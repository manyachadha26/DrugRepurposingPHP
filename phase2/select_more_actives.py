#!/usr/bin/env python3
"""
Pick 110 more PDE4 actives from the pool already on disk, to reach 150.

WHY EXPAND
----------
The interaction fingerprint passed at AUC 0.727, but with two weaknesses that
both come from having only 40 actives:

  - the 95% confidence interval runs 0.637 to 0.817, so the lower bound sits
    below the 0.70 pass bar
  - the classifier uses 44 features to separate 40 actives, so it is
    over-parameterised and the feature weights are unstable - they reshuffled
    completely between the dry and wet runs

150 actives fixes both. Nothing about the first 40 is wasted: they stay in
the benchmark and this only adds to them.

NO NETWORK NEEDED
-----------------
The 5,349-compound PDE4 activity pool is already cached in chembl_raw.pkl, so
this runs offline. No ChEMBL calls, no HTTP 500s, no waiting on their server.

HOW THE 110 ARE CHOSEN
----------------------
Same rules as the first 40, so the combined set stays homogeneous:

  - pIC50 between 7.0 and 10.0
  - molecular weight 250 to 550
  - one compound per Bemis-Murcko scaffold, and no scaffold already used by
    the existing 40
  - a RANDOM member of each scaffold group, never the most potent one.
    Taking the most potent gave a first draft spanning pIC50 9.7 to 11.1,
    all picomolar, which is not what a screen encounters
  - spread evenly across five potency bands, so the set is not all
    sub-nanomolar at one end

USAGE
-----
    python select_more_actives.py chembl_raw.pkl actives_chembl.smi \\
           [n_wanted] [seed]

OUTPUT
------
    actives_extra.smi       the 110 new ones, SMILES + ChEMBL ID
    actives_extra.csv       with pIC50, MW and scaffold, for the paper
    actives_150.smi         all 150 together, ready to prepare and dock
"""

import os
import sys
import pickle
import numpy as np
import pandas as pd

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors
    from rdkit.Chem.Scaffolds import MurckoScaffold
    RDLogger.DisableLog("rdApp.*")
except ImportError as exc:
    sys.exit(f"ERROR: rdkit import failed ({exc}). Activate your env.")

PIC50_MIN, PIC50_MAX = 7.0, 10.0
MW_RANGE = (250.0, 550.0)
N_BANDS = 5


def find_col(df, *keys):
    """Locate a column by substring, since ChEMBL field names vary."""
    for k in keys:
        for c in df.columns:
            if k in c.lower():
                return c
    return None


ASSAY_TYPES = {"ic50", "ki", "kd", "ec50"}


def to_pic50(row, pchembl_col, value_col, units_col, type_col,
             rel_col=None, assume_nm=False):
    """
    Convert an activity record to pIC50.

    Takes pchembl_value when the row has one, and otherwise converts
    standard_value using its units. Preferring one column outright and
    giving up when it is blank throws away every record that carries only a
    raw value - in testing that lost 30% of the molecules for no reason.

    ChEMBL mixes nM, uM and pre-converted values in the same table, so units
    are read rather than assumed. Assuming wrong moves a compound by three
    orders of magnitude.
    """
    if pchembl_col:
        v = row.get(pchembl_col)
        if v is not None and not pd.isna(v):
            try:
                return float(v)
            except (TypeError, ValueError):
                pass

    if not value_col:
        return None

    # A censored value is not a measurement. ">8000 nM" means "weaker than
    # 8000", so reading it as exactly 8000 invents a potency the assay never
    # established - and these rows are common in ChEMBL.
    if rel_col:
        if str(row.get(rel_col, "=")).strip() not in ("=", ""):
            return None

    # Only direct binding or inhibition endpoints. Percent-inhibition and
    # similar rows cannot be converted and would otherwise slip through.
    t_raw = str(row.get(type_col, "")).strip().lower()
    if t_raw and t_raw not in ASSAY_TYPES and not t_raw.startswith("p"):
        return None

    v = row.get(value_col)
    if v is None or pd.isna(v):
        return None
    try:
        v = float(v)
    except (TypeError, ValueError):
        return None
    if v <= 0:
        return None

    t = str(row.get(type_col, "")).lower()
    if t.startswith("pic") or t.startswith("pki"):
        return v                      # already a log value
    u = str(row.get(units_col, "")).strip().lower() if units_col else ""
    factor = {"nm": 1e-9, "um": 1e-6, "µm": 1e-6, "mm": 1e-3,
              "m": 1.0, "pm": 1e-12}.get(u)
    if factor is None:
        if not assume_nm:
            return None
        # No units column was cached. ChEMBL stores standard_value for
        # IC50/Ki/Kd in nM, so nM is the right assumption - but it IS an
        # assumption and belongs in the methods section.
        factor = 1e-9
    return -np.log10(v * factor)


def scaffold_of(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None, None, None
    try:
        s = Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(m))
    except Exception:
        return None, None, None
    if not s:
        s = Chem.MolToSmiles(m)          # acyclic: its own group
    return s, Descriptors.MolWt(m), Chem.MolToSmiles(m)


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    pkl_path, existing_path = sys.argv[1], sys.argv[2]
    n_want = int(sys.argv[3]) if len(sys.argv) > 3 else 110
    seed = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    for p in (pkl_path, existing_path):
        if not os.path.exists(p):
            sys.exit(f"ERROR: not found: {p}")

    with open(pkl_path, "rb") as fh:
        raw = pickle.load(fh)
    df = raw if isinstance(raw, pd.DataFrame) else pd.DataFrame(raw)
    print(f"Cached pool: {len(df)} activity records, "
          f"{len(df.columns)} columns")

    smi_col = find_col(df, "canonical_smiles", "smiles")
    id_col = find_col(df, "molecule_chembl_id", "chembl_id")
    pch_col = find_col(df, "pchembl_value")
    val_col = find_col(df, "standard_value")
    uni_col = find_col(df, "standard_units", "units")
    typ_col = find_col(df, "standard_type", "type")
    rel_col = find_col(df, "standard_relation", "relation")
    if not smi_col or not id_col or not (pch_col or val_col):
        sys.exit(f"ERROR: could not find the needed columns.\n"
                 f"  smiles={smi_col} id={id_col} pchembl={pch_col} "
                 f"value={val_col}\n  available: {list(df.columns)}")
    print(f"Using columns: smiles={smi_col}, id={id_col}, "
          f"pchembl={pch_col}, value={val_col}, units={uni_col}, "
          f"type={typ_col}, relation={rel_col}")
    assume_nm = (pch_col is None and uni_col is None)
    if assume_nm:
        print("  NO units column in the cache. Assuming nM, which is what")
        print("  ChEMBL stores standard_value in for IC50/Ki/Kd. State this")
        print("  assumption in the methods section.")
    if rel_col:
        print("  Censored rows (>, <, >=) are dropped: they are bounds, not")
        print("  measurements.")
    else:
        print("  WARNING: no relation column, so censored values cannot be")
        print("  dropped and some potencies may be bounds, not measurements.")

    # ---- the 40 already chosen, which must be excluded ----
    have_ids, have_scaffolds = set(), set()
    with open(existing_path) as fh:
        for line in fh:
            p = line.split()
            if len(p) < 2:
                continue
            a, b = p[0], p[1]
            cid = b if b.startswith("CHEMBL") else a
            smi = a if cid == b else b
            have_ids.add(cid)
            # Take the scaffold from the existing file's OWN smiles. Looking
            # it up in the filtered candidate table instead misses any
            # existing active that the pIC50 or MW filter removed, and its
            # scaffold would then be free to reappear in the new picks -
            # which silently breaks the leave-one-scaffold-out split.
            sc, _, _ = scaffold_of(smi)
            if sc:
                have_scaffolds.add(sc)
    print(f"Already selected: {len(have_ids)} actives to exclude, "
          f"{len(have_scaffolds)} of their scaffolds resolved")
    if len(have_scaffolds) < len(have_ids):
        print(f"  WARNING: {len(have_ids) - len(have_scaffolds)} of the "
              f"existing actives' SMILES would not parse. Their scaffolds")
        print(f"  cannot be excluded, so a collision is possible. Check "
              f"{existing_path}.")

    # ---- pIC50 per record, then best record per molecule ----
    df = df.dropna(subset=[smi_col, id_col]).copy()
    df["_pic50"] = df.apply(
        lambda r: to_pic50(r, pch_col, val_col, uni_col, typ_col,
                           rel_col, assume_nm), axis=1)
    df = df.dropna(subset=["_pic50"])
    best = df.sort_values("_pic50", ascending=False).groupby(
        id_col, as_index=False).first()
    print(f"Unique molecules with a usable pIC50: {len(best)}")

    in_band = best[(best["_pic50"] >= PIC50_MIN)
                   & (best["_pic50"] <= PIC50_MAX)]
    print(f"  within pIC50 {PIC50_MIN}-{PIC50_MAX}: {len(in_band)}")

    # ---- scaffolds and weights ----
    rows = []
    for _, r in in_band.iterrows():
        sc, mw, can = scaffold_of(r[smi_col])
        if sc is None or mw is None:
            continue
        if not (MW_RANGE[0] <= mw <= MW_RANGE[1]):
            continue
        rows.append({"chembl_id": r[id_col], "smiles": can,
                     "pIC50": round(float(r["_pic50"]), 2),
                     "MW": round(mw, 1), "scaffold": sc})
    cand = pd.DataFrame(rows)
    print(f"  parseable and MW {MW_RANGE[0]:.0f}-{MW_RANGE[1]:.0f}: "
          f"{len(cand)}")
    if cand.empty:
        sys.exit("ERROR: nothing left after filtering. Check the columns "
                 "printed above are the right ones.")

    # ---- exclude the scaffolds the existing 40 already occupy ----
    print(f"  scaffolds used by the existing actives: "
          f"{len(have_scaffolds)}")
    cand = cand[~cand.chembl_id.isin(have_ids)
                & ~cand.scaffold.isin(have_scaffolds)]
    print(f"  candidates on new scaffolds: {len(cand)} "
          f"across {cand.scaffold.nunique()} scaffolds")

    # ---- one random member per scaffold ----
    rng = np.random.default_rng(seed)
    picks = (cand.sample(frac=1.0, random_state=seed)
                 .groupby("scaffold", as_index=False).first())
    print(f"  one per scaffold: {len(picks)}")

    if len(picks) < n_want:
        print(f"\nWARNING: only {len(picks)} distinct new scaffolds are "
              f"available, fewer than the {n_want} requested.")
        print("Taking all of them. The benchmark will be smaller than "
              "planned; report the real number, do not reuse scaffolds.")
        n_want = len(picks)

    # ---- spread across potency bands ----
    edges = np.linspace(PIC50_MIN, PIC50_MAX, N_BANDS + 1)
    picks["_band"] = np.clip(
        np.digitize(picks["pIC50"], edges[1:-1]), 0, N_BANDS - 1)
    per = n_want // N_BANDS
    chosen = []
    for b in range(N_BANDS):
        grp = picks[picks._band == b]
        take = min(per, len(grp))
        if take:
            chosen.append(grp.sample(take, random_state=seed + b))
        print(f"  band {edges[b]:.1f}-{edges[b+1]:.1f}: "
              f"{len(grp)} available, taking {take}")
    out = pd.concat(chosen) if chosen else picks.head(0)

    # top up from whatever is left if some bands were thin
    if len(out) < n_want:
        left = picks[~picks.chembl_id.isin(out.chembl_id)]
        need = n_want - len(out)
        if len(left):
            out = pd.concat([out, left.sample(min(need, len(left)),
                                              random_state=seed + 99)])
            print(f"  topped up with {min(need, len(left))} from other bands")

    out = out.drop(columns=["_band"]).reset_index(drop=True)

    out[["smiles", "chembl_id"]].to_csv(
        "actives_extra.smi", sep="\t", header=False, index=False)
    out.to_csv("actives_extra.csv", index=False)

    # ---- the combined 150 ----
    with open("actives_150.smi", "w") as fh:
        with open(existing_path) as src:
            for line in src:
                if line.strip():
                    fh.write(line if line.endswith("\n") else line + "\n")
        for _, r in out.iterrows():
            fh.write(f"{r['smiles']}\t{r['chembl_id']}\n")

    n_total = len(have_ids) + len(out)
    print(f"\nWrote actives_extra.smi   {len(out)} new actives")
    print(f"Wrote actives_extra.csv   with pIC50, MW and scaffold")
    print(f"Wrote actives_150.smi     {n_total} actives in total")
    print(f"\npIC50 of the new ones: min {out.pIC50.min():.2f}  "
          f"median {out.pIC50.median():.2f}  max {out.pIC50.max():.2f}")
    print(f"MW of the new ones   : min {out.MW.min():.0f}  "
          f"median {out.MW.median():.0f}  max {out.MW.max():.0f}")
    print(f"Distinct scaffolds   : {out.scaffold.nunique()} "
          f"(should equal {len(out)})")
    if out.scaffold.nunique() != len(out):
        print("  WARNING: a scaffold repeats. That breaks the "
              "leave-one-scaffold-out split. Tell me before docking.")

    print(f"\nNext: generate 15 decoys for each of these {len(out)}, then "
          f"dock.")


if __name__ == "__main__":
    main()
