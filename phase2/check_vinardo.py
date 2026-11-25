#!/usr/bin/env python3
"""
Is Vinardo's early enrichment real, or is it the size artifact again?

BACKGROUND
----------
Gate 2 round 1 looked like it had weak signal, but a diagnostic showed the
scores were really tracking molecular weight (r = -0.739): bigger molecule,
more contacts, better score. That is not drug discovery, that is a ruler.

Vinardo now gives 7.7x enrichment in the top 1% (3 of 40 actives where 0.39
is expected, p = 0.0056) while Vina gives zero. Before we build a screen on
that, we have to rule out the same confound.

FOUR CHECKS
-----------
1. MW vs score correlation. If |r| is large and negative, the score is a
   size proxy and the enrichment is an artifact. We want |r| small.
2. Are the top-12 actives unusually heavy compared to the other actives?
   If the 3 hits are the 3 biggest actives, that is the artifact.
3. Does Vinardo rank differently from Vina at all, or is it the same order
   with different numbers? Spearman close to 1.0 means nothing changed.
4. Are the top hits the most potent actives, or a random potency draw?
   Real signal should favour potent compounds.

USAGE
-----
    python check_vinardo.py rescore_results.csv actives_chembl.smi \\
           decoys_chembl.smi [actives_chembl.csv]

The last argument is optional - it supplies pIC50 for check 4.
"""

import os
import sys
import numpy as np
import pandas as pd

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors
    RDLogger.DisableLog("rdApp.*")
except ImportError as exc:
    sys.exit(f"ERROR: rdkit import failed ({exc}). "
             f"conda activate your project env.")

from scipy.stats import spearmanr, pearsonr, mannwhitneyu


def read_smi(path):
    """
    Read a .smi file into {name: smiles}.

    Handles both column orders, because .smi files are written either way
    and guessing wrong would silently produce no matches.
    """
    out = {}
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 2:
                continue
            a, b = parts[0], parts[1]
            # whichever token parses as a molecule is the SMILES
            if Chem.MolFromSmiles(a) is not None:
                out[b] = a
            elif Chem.MolFromSmiles(b) is not None:
                out[a] = b
    return out


def main():
    if len(sys.argv) < 4:
        print(__doc__)
        sys.exit(1)

    df = pd.read_csv(sys.argv[1])
    smiles = {}
    for p in sys.argv[2:4]:
        smiles.update(read_smi(p))
    print(f"SMILES loaded for {len(smiles)} compounds\n")

    # --- attach MW and heavy atom count ---
    mw, hac = [], []
    for n in df["name"]:
        s = smiles.get(n)
        m = Chem.MolFromSmiles(s) if s else None
        mw.append(Descriptors.MolWt(m) if m else np.nan)
        hac.append(m.GetNumHeavyAtoms() if m else np.nan)
    df["MW"] = mw
    df["HAC"] = hac

    miss = int(df["MW"].isna().sum())
    if miss:
        print(f"WARNING: no SMILES matched for {miss} of {len(df)} names. "
              f"Checks below use the {len(df)-miss} that matched.\n")
    d = df.dropna(subset=["MW"]).copy()

    # ---------- CHECK 1: size confound ----------
    print("=" * 62)
    print("CHECK 1  Is the score just measuring molecular size?")
    print("=" * 62)
    for col in ("vina_score", "vinardo_score"):
        if col not in d:
            continue
        r_mw, p_mw = pearsonr(d["MW"], d[col])
        r_ha, _ = pearsonr(d["HAC"], d[col])
        flag = ("ARTIFACT - score is a size proxy" if abs(r_mw) > 0.60
                else "borderline" if abs(r_mw) > 0.40
                else "OK - not size driven")
        print(f"  {col:<15} r(MW) = {r_mw:+.3f}   "
              f"r(heavy atoms) = {r_ha:+.3f}   -> {flag}")
    print("  (round 1 Vina was r = -0.739, which was the artifact)\n")

    # ---------- CHECK 2: are the top hits just the heavy actives? ----------
    print("=" * 62)
    print("CHECK 2  Are the top-ranked actives the heaviest actives?")
    print("=" * 62)
    col = "vinardo_score" if "vinardo_score" in d else "vina_score"
    ranked = d.sort_values(col).reset_index(drop=True)
    ranked["rank"] = np.arange(1, len(ranked) + 1)
    k = max(1, int(round(len(ranked) * 0.01)))
    acts = ranked[ranked.label == 1]
    top_acts = acts[acts["rank"] <= k]
    rest_acts = acts[acts["rank"] > k]

    print(f"  all {len(acts)} actives      median MW {acts.MW.median():.0f}")
    if len(top_acts):
        print(f"  the {len(top_acts)} in the top {k}   median MW "
              f"{top_acts.MW.median():.0f}   "
              f"(individually: {', '.join(f'{v:.0f}' for v in top_acts.MW)})")
        if len(rest_acts) > 1:
            u, pu = mannwhitneyu(top_acts.MW, rest_acts.MW,
                                 alternative="greater")
            print(f"  heavier than the other actives? p = {pu:.3f}"
                  f"  -> {'YES, suspicious' if pu < 0.05 else 'no, fine'}")
    else:
        print(f"  no actives in the top {k}")
    print()

    # ---------- CHECK 3: did the ranking actually change? ----------
    if "vina_score" in d and "vinardo_score" in d:
        print("=" * 62)
        print("CHECK 3  Did Vinardo actually reorder anything?")
        print("=" * 62)
        rho, _ = spearmanr(d["vina_score"], d["vinardo_score"])
        print(f"  Spearman between the two scores: {rho:+.3f}")
        print("  1.00 would mean identical order and no new information.")
        if rho > 0.95:
            print("  -> near-identical ranking. The enrichment difference is")
            print("     then driven by a handful of swaps at the very top.")
        else:
            print("  -> genuinely different ranking, not a rescaling.")
        v_top = set(d.nsmallest(12, "vina_score")["name"])
        d_top = set(d.nsmallest(12, "vinardo_score")["name"])
        print(f"  top-12 overlap: {len(v_top & d_top)} of 12 compounds\n")

    # ---------- CHECK 4: potency ----------
    if len(sys.argv) > 4 and os.path.exists(sys.argv[4]):
        print("=" * 62)
        print("CHECK 4  Are the top hits the potent actives?")
        print("=" * 62)
        a = pd.read_csv(sys.argv[4])
        namecol = next((c for c in a.columns
                        if "chembl" in c.lower() and "id" in c.lower()), None)
        pcol = next((c for c in a.columns if "pic50" in c.lower()), None)
        if namecol and pcol:
            pot = dict(zip(a[namecol], a[pcol]))
            acts = acts.copy()
            acts["pIC50"] = acts["name"].map(pot)
            got = acts.dropna(subset=["pIC50"])
            if len(got):
                rho_p, p_p = spearmanr(got["rank"], got["pIC50"])
                print(f"  rank vs pIC50 across {len(got)} actives: "
                      f"rho = {rho_p:+.3f}, p = {p_p:.3f}")
                print("  negative rho = potent compounds rank better "
                      "(what we want)")
                t = got[got["rank"] <= k]
                if len(t):
                    print(f"  top-{k} actives pIC50: "
                          f"{', '.join(f'{v:.1f}' for v in t.pIC50)}")
                print(f"  all actives median pIC50: {got.pIC50.median():.1f}")
        else:
            print("  could not find ChEMBL ID / pIC50 columns in that csv; "
                  "skipping")
        print()

    print("=" * 62)
    print("HOW TO READ THIS")
    print("=" * 62)
    print("  Check 1 OK and Check 2 'no'  -> the Vinardo early enrichment is")
    print("      real. Use Vinardo to rank the screen, still confirm affinity")
    print("      with MM-GBSA, and report both scoring functions.")
    print("  Check 1 ARTIFACT or Check 2 'YES' -> it is the size effect again.")
    print("      Do not rank on score at all; triage on pocket occupancy and")
    print("      interactions, and let MD do the ranking.")


if __name__ == "__main__":
    main()
