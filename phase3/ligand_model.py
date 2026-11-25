#!/usr/bin/env python3
"""
A ligand-based model: predict PDE4 activity from 2D structure alone.

WHY A SECOND MODEL
------------------
The interaction fingerprint classifier reached AUC 0.822, but it needs a
docked pose for every compound it scores. This one uses only the molecule's
2D structure - no protein, no docking - so it is an independent line of
evidence built from a different feature space with a different algorithm.

If both models rank sotagliflozin and sitagliptin highly, that is real
corroboration. If they disagree, we need to know before writing anything.

It is trained on the SAME 146 actives and 2,246 matched decoys, so its AUC
sits directly beside the structure-based one in a single table.

WHAT IT LEARNS, AND WHAT THAT IS WORTH
--------------------------------------
The decoys are PRESUMED inactive - property-matched compounds with no
recorded PDE4 data - not measured inactives. So strictly the model learns
"known PDE4 actives vs property-matched compounds without PDE4 data", which
is the same thing the structure-based model learns. That is the right
comparison, and it is not the same as a true activity model. Say so in the
methods.

VALIDATION
----------
Scaffold-grouped cross-validation: every compound sharing a Bemis-Murcko
scaffold goes into the same fold, so a prediction never comes from a model
that saw that chemotype. A random split leaks badly here - analogues land
on both sides and the model learns the scaffold, not the chemistry.

A label-permutation control runs alongside. Shuffle the labels, repeat
everything, and a sound pipeline collapses to AUC near 0.5.

USAGE
-----
    python ligand_model.py actives_150.smi decoys_150.smi [screen.smi]
                           [--loso] [--folds N]

    screen.smi   optional: a library to rank once the model is trained
    --loso       leave-one-scaffold-out instead of grouped k-fold. Matches
                 the structure-based model exactly, but 146 folds of a
                 random forest takes a while.
    --folds N    number of grouped folds (default 10)

OUTPUT
------
    ligand_model_oof.csv       out-of-fold prediction per training compound
    ligand_model_screen.csv    the screen library, ranked (if one was given)
"""

import os
import sys
import numpy as np
import pandas as pd

try:
    from rdkit import Chem, RDLogger, DataStructs
    from rdkit.Chem import Descriptors, Crippen, Lipinski
    from rdkit.Chem.Scaffolds import MurckoScaffold
    from rdkit.Chem import rdFingerprintGenerator as fpg
    RDLogger.DisableLog("rdApp.*")
except ImportError as exc:
    sys.exit(f"ERROR: rdkit import failed ({exc}). Activate your env.")

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GroupKFold

N_BITS = 2048
RADIUS = 2          # ECFP4
ALPHA = 20.0


# ----------------------------------------------------------------- metrics
def metrics(scores, labels, alpha=ALPHA):
    """Higher score = more likely active, so negate for the shared ranking."""
    order = np.argsort(-np.asarray(scores))
    y = np.asarray(labels)[order]
    n = len(y); n_act = int(y.sum()); n_dec = n - n_act
    if n_act == 0 or n_dec == 0:
        raise ValueError("need both classes")
    ranks = np.arange(1, n + 1)
    auc = 1 - ((ranks[y == 1].sum() - n_act * (n_act + 1) / 2) / (n_act * n_dec))
    efs = {}
    for pct in (1, 5, 10):
        k = max(1, int(round(n * pct / 100)))
        found = int(y[:k].sum()); exp = n_act * k / n
        efs[pct] = (found, k, found / exp if exp else 0.0)
    ri = np.where(y == 1)[0] + 1
    ra = n_act / n
    s = np.sum(np.exp(-alpha * ri / n)) / n_act
    s_rand = (1.0 / n) * (1 - np.exp(-alpha)) / (np.exp(alpha / n) - 1)
    rie = s / s_rand
    hi = (1 - np.exp(-alpha * ra)) / (ra * (1 - np.exp(-alpha)))
    lo = (1 - np.exp(alpha * ra)) / (ra * (1 - np.exp(alpha)))
    bedroc = float(np.clip((rie - lo) / (hi - lo), 0, 1))
    return auc, efs, bedroc


def se_auc(A, n1, n0):
    Q1 = A / (2 - A); Q2 = 2 * A * A / (1 + A)
    return float(np.sqrt((A * (1 - A) + (n1 - 1) * (Q1 - A * A)
                          + (n0 - 1) * (Q2 - A * A)) / (n1 * n0)))


# ------------------------------------------------------------------ inputs
def read_smi(path):
    """[(name, smiles)], tolerant of either column order."""
    out = []
    with open(path) as fh:
        for line in fh:
            p = line.split()
            if len(p) < 2:
                continue
            a, b = p[0], p[1]
            if Chem.MolFromSmiles(a) is not None:
                out.append((b, a))
            elif Chem.MolFromSmiles(b) is not None:
                out.append((a, b))
    return out


DESCRIPTORS = [
    ("MW", Descriptors.MolWt),
    ("logP", Crippen.MolLogP),
    ("TPSA", Descriptors.TPSA),
    ("HBD", Lipinski.NumHDonors),
    ("HBA", Lipinski.NumHAcceptors),
    ("RotB", Lipinski.NumRotatableBonds),
    ("AromRings", Lipinski.NumAromaticRings),
    ("Rings", Lipinski.RingCount),
    ("FracCsp3", Lipinski.FractionCSP3),
    ("HeavyAtoms", Descriptors.HeavyAtomCount),
    ("MolMR", Crippen.MolMR),
    ("BalabanJ", Descriptors.BalabanJ),
    ("BertzCT", Descriptors.BertzCT),
    ("NHOHCount", Lipinski.NHOHCount),
    ("NOCount", Lipinski.NOCount),
    ("HeteroAtoms", Lipinski.NumHeteroatoms),
]


def featurise(smiles_list):
    """ECFP4 bits plus physicochemical descriptors."""
    gen = fpg.GetMorganGenerator(radius=RADIUS, fpSize=N_BITS)
    X, keep = [], []
    for i, smi in enumerate(smiles_list):
        m = Chem.MolFromSmiles(smi)
        if m is None:
            continue
        fp = gen.GetFingerprint(m)
        arr = np.zeros((N_BITS,), dtype=np.int8)
        DataStructs.ConvertToNumpyArray(fp, arr)
        desc = []
        for _, fn in DESCRIPTORS:
            try:
                v = fn(m)
            except Exception:
                v = 0.0
            desc.append(0.0 if v is None or np.isnan(v) else float(v))
        X.append(np.concatenate([arr, np.array(desc)]))
        keep.append(i)
    names = [f"bit_{i}" for i in range(N_BITS)] + [d[0] for d in DESCRIPTORS]
    return np.array(X, dtype=float), keep, names


def scaffold_of(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    try:
        s = Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(m))
    except Exception:
        return None
    return s if s else Chem.MolToSmiles(m)


def new_forest(seed=0):
    return RandomForestClassifier(
        n_estimators=500, max_depth=None, min_samples_leaf=2,
        max_features="sqrt", class_weight="balanced_subsample",
        n_jobs=-1, random_state=seed,
    )


def cross_validate(X, y, groups, use_loso, n_folds, seed=0, quiet=False):
    """Out-of-fold probabilities, never predicting a scaffold it trained on."""
    proba = np.full(len(y), np.nan)
    if use_loso:
        folds = [(np.array([g != gg for gg in groups]),
                  np.array([g == gg for gg in groups]))
                 for g in sorted({g for g, l in zip(groups, y) if l == 1})]
    else:
        gk = GroupKFold(n_splits=n_folds)
        folds = [(np.isin(np.arange(len(y)), tr),
                  np.isin(np.arange(len(y)), te))
                 for tr, te in gk.split(X, y, groups)]
    for i, (tr, te) in enumerate(folds, 1):
        if y[tr].sum() < 2 or (1 - y[tr]).sum() < 2:
            continue
        clf = new_forest(seed + i)
        clf.fit(X[tr], y[tr])
        proba[te] = clf.predict_proba(X[te])[:, 1]
        if not quiet and (i % 10 == 0 or i == len(folds)):
            print(f"    fold {i} / {len(folds)}")
    return proba


def report(tag, proba, y):
    ok = ~np.isnan(proba)
    auc, efs, bedroc = metrics(proba[ok], y[ok])
    n1 = int(y[ok].sum()); n0 = int(ok.sum()) - n1
    s = se_auc(auc, n1, n0)
    print(f"\n=== {tag} ===")
    print(f"  {n1} actives, {n0} decoys, {int(ok.sum())} predicted")
    print(f"  AUC-ROC   {auc:.3f}   95% CI "
          f"[{max(0,auc-1.96*s):.3f}, {min(1,auc+1.96*s):.3f}]")
    for pct, (found, k, ef) in efs.items():
        print(f"  EF {pct:>2}%     {ef:>5.1f}x   ({found} of {n1} in top {k})")
    print(f"  BEDROC    {bedroc:.3f}")
    return auc


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    use_loso = "--loso" in sys.argv
    n_folds = 10
    for a in sys.argv:
        if a.startswith("--folds"):
            try:
                n_folds = int(a.split("=")[1]) if "=" in a else n_folds
            except (IndexError, ValueError):
                pass
    if len(args) < 2:
        print(__doc__)
        sys.exit(1)
    for p in args[:2]:
        if not os.path.exists(p):
            sys.exit(f"ERROR: not found: {p}")

    acts = read_smi(args[0])
    decs = read_smi(args[1])
    names = [n for n, _ in acts] + [n for n, _ in decs]
    smis = [s for _, s in acts] + [s for _, s in decs]
    labels = np.array([1] * len(acts) + [0] * len(decs))
    print(f"Read {len(acts)} actives and {len(decs)} decoys")

    X, keep, feat_names = featurise(smis)
    y = labels[keep]
    names = [names[i] for i in keep]
    smis = [smis[i] for i in keep]
    groups = [scaffold_of(s) for s in smis]
    bad = [i for i, g in enumerate(groups) if g is None]
    if bad:
        mask = np.array([i not in set(bad) for i in range(len(groups))])
        X, y = X[mask], y[mask]
        names = [n for i, n in enumerate(names) if mask[i]]
        groups = [g for i, g in enumerate(groups) if mask[i]]
        print(f"  dropped {len(bad)} compounds with no scaffold")

    n_act = int(y.sum())
    n_sc = len(set(groups))
    print(f"Featurised {len(y)}: {n_act} actives, {len(y)-n_act} decoys")
    print(f"Features   : {X.shape[1]} "
          f"({N_BITS} ECFP4 bits + {len(DESCRIPTORS)} descriptors)")
    print(f"Scaffolds  : {n_sc}")
    print(f"Validation : "
          + ("leave-one-scaffold-out" if use_loso
             else f"{n_folds}-fold, grouped by scaffold"))
    print("\nCross-validating (this is the number that counts):")

    proba = cross_validate(X, y, groups, use_loso, n_folds)
    auc = report("ligand-based random forest, out-of-fold", proba, y)

    pd.DataFrame({"name": names, "label": y, "scaffold": groups,
                  "probability": proba}).to_csv(
        "ligand_model_oof.csv", index=False)

    # ---------------- permutation control ----------------
    print("\n=== label permutation control ===")
    print("  Shuffled labels should collapse to AUC near 0.5.")
    perm = []
    for seed in range(3):
        rng = np.random.default_rng(seed)
        ysh = y.copy(); rng.shuffle(ysh)
        p = cross_validate(X, ysh, groups, use_loso, n_folds,
                           seed=100 + seed, quiet=True)
        ok = ~np.isnan(p)
        try:
            perm.append(metrics(p[ok], ysh[ok])[0])
        except ValueError:
            pass
    if perm:
        m, s = float(np.mean(perm)), float(np.std(perm))
        print(f"  shuffled AUC: {m:.3f} +/- {s:.3f} over {len(perm)} runs")
        if m > 0.60:
            print("  *** LEAKAGE. The real number above cannot be trusted.")

    # ---------------- what the model uses ----------------
    clf = new_forest().fit(X, y)
    imp = clf.feature_importances_
    print("\n=== most informative features ===")
    print("  Fit on all data, for interpretation only - NOT performance.")
    desc_idx = list(range(N_BITS, X.shape[1]))
    top_desc = sorted(desc_idx, key=lambda i: -imp[i])[:8]
    print("  descriptors:")
    for i in top_desc:
        print(f"    {feat_names[i]:<14} {imp[i]:.4f}")
    bits = sorted(range(N_BITS), key=lambda i: -imp[i])[:5]
    print(f"  top ECFP4 bits: "
          + ", ".join(f"{feat_names[i]}({imp[i]:.4f})" for i in bits))
    print("  A descriptor dominating here means the model may be keying on")
    print("  bulk properties rather than chemistry - compare against the")
    print("  structure-based model before drawing conclusions.")

    # ---------------- rank the screen ----------------
    if len(args) > 2 and os.path.exists(args[2]):
        lib = read_smi(args[2])
        Xl, keepl, _ = featurise([s for _, s in lib])
        lnames = [lib[i][0] for i in keepl]
        p = clf.predict_proba(Xl)[:, 1]
        pa = clf.predict_proba(X[y == 1])[:, 1]
        out = pd.DataFrame({
            "name": lnames, "probability": p,
            "pct_of_actives_below": [float((pa < v).mean() * 100) for v in p],
        }).sort_values("probability", ascending=False).reset_index(drop=True)
        out.insert(0, "rank", np.arange(1, len(out) + 1))
        out.to_csv("ligand_model_screen.csv", index=False)

        print(f"\n=== the screen, ranked by the LIGAND-BASED model ===")
        print("  (the model has seen none of these compounds)")
        print(f"  {'#':>3}  {'compound':<26} {'prob':>6}  {'%act':>5}")
        print("  " + "-" * 46)
        for _, r in out.iterrows():
            print(f"  {int(r['rank']):>3}  {str(r['name'])[:26]:<26} "
                  f"{r['probability']:>6.3f}  "
                  f"{r['pct_of_actives_below']:>4.0f}%")
        print("\n  Compare these ranks against the structure-based ranking.")
        print("  Agreement on a compound is corroboration from an")
        print("  independent method. Disagreement must be reported, not")
        print("  resolved by picking the ranking you prefer.")


if __name__ == "__main__":
    main()
