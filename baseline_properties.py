#!/usr/bin/env python3
"""
Experiment B - what does trivial physicochemistry alone score on this benchmark?

The interaction-fingerprint model reaches AUC 0.822. That number only means
something if simple molecular properties cannot reach it. This fits a logistic
regression on SIX descriptors only - no interaction features, no pose, no
docking at all - using exactly the same leave-one-scaffold-out folds, so the
two numbers are directly comparable.

Report whatever it gives. If this baseline is high, the fingerprint result is
weaker than it looks, and that has to be said in the manuscript.

    python baseline_properties.py fingerprints_150.csv actives_150.smi decoys_150.smi

Writes baseline_properties_oof.csv (per-compound out-of-fold probabilities,
needed by delong.py) and prints the metrics.
"""
import sys

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger, DataStructs
from rdkit.Chem import Descriptors, Crippen
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline

RDLogger.DisableLog("rdApp.*")
ALPHA = 20.0


def metrics(score_lower_better, labels, alpha=ALPHA):
    order = np.argsort(score_lower_better)
    y = np.asarray(labels)[order]
    n, na = len(y), int(np.sum(labels))
    nd = n - na
    ranks = np.arange(1, n + 1)
    auc = 1 - ((ranks[y == 1].sum() - na * (na + 1) / 2) / (na * nd))
    efs = {}
    for pct in (1, 5, 10):
        k = max(1, int(round(n * pct / 100)))
        found = int(y[:k].sum())
        efs[pct] = (found, k, found / (na * k / n))
    ri = np.where(y == 1)[0] + 1
    ra = na / n
    s = np.sum(np.exp(-alpha * ri / n)) / na
    s_rand = (1.0 / n) * (1 - np.exp(-alpha)) / (np.exp(alpha / n) - 1)
    rie = s / s_rand
    rie_max = (1 - np.exp(-alpha * ra)) / (ra * (1 - np.exp(-alpha)))
    rie_min = (1 - np.exp(alpha * ra)) / (ra * (1 - np.exp(alpha)))
    bedroc = float(np.clip((rie - rie_min) / (rie_max - rie_min), 0, 1))
    return auc, efs, bedroc


def read_smi(path):
    out = {}
    for line in open(path):
        p = line.split()
        if len(p) < 2:
            continue
        a, b = p[0], p[1]
        if Chem.MolFromSmiles(a) is not None:
            out[b] = a
        elif Chem.MolFromSmiles(b) is not None:
            out[a] = b
    return out


def scaffold_of(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    s = Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(m))
    return s if s else Chem.MolToSmiles(m)


def props(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    return [Descriptors.MolWt(m), Crippen.MolLogP(m), Descriptors.TPSA(m),
            Descriptors.NumRotatableBonds(m), Descriptors.NumHDonors(m),
            Descriptors.NumHAcceptors(m)]


def make_folds(y, groups):
    by = {}
    for i, g in enumerate(groups):
        by.setdefault(g, []).append(i)
    act = sorted({g for g, l in zip(groups, y) if l == 1})
    fold = np.full(len(y), -1, int)
    for f, g in enumerate(act):
        for i in by[g]:
            fold[i] = f
    others = sorted((g for g in by if g not in set(act)),
                    key=lambda g: -len(by[g]))
    sizes = np.array([(fold == f).sum() for f in range(len(act))], float)
    for g in others:
        f = int(np.argmin(sizes))
        for i in by[g]:
            fold[i] = f
        sizes[f] += len(by[g])
    return fold


def loso(X, y, fold, seed=0):
    proba = np.full(len(y), np.nan)
    for f in np.unique(fold[fold >= 0]):
        te = fold == f
        tr = ~te
        if y[tr].sum() < 2 or (1 - y[tr]).sum() < 2:
            continue
        clf = make_pipeline(StandardScaler(),
                            LogisticRegression(max_iter=2000,
                                               class_weight="balanced",
                                               random_state=seed))
        clf.fit(X[tr], y[tr])
        proba[te] = clf.predict_proba(X[te])[:, 1]
    return proba


def main():
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    fp = pd.read_csv(sys.argv[1])
    smiles = {}
    for p in sys.argv[2:4]:
        smiles.update(read_smi(p))

    rows, keepidx = [], []
    for i, n in enumerate(fp["name"]):
        if n not in smiles:
            continue
        pr = props(smiles[n])
        sc = scaffold_of(smiles[n])
        if pr is None or sc is None:
            continue
        rows.append(pr + [sc, n, int(fp["label"].iloc[i])])
        keepidx.append(i)

    cols = ["MW", "logP", "TPSA", "RotB", "HBD", "HBA", "scaffold",
            "name", "label"]
    d = pd.DataFrame(rows, columns=cols)
    dropped = len(fp) - len(d)
    if dropped:
        print(f"WARNING: {dropped} of {len(fp)} compounds had no SMILES match "
              f"and were dropped. Check the .smi files are the right ones.")

    X = d[["MW", "logP", "TPSA", "RotB", "HBD", "HBA"]].values.astype(float)
    y = d["label"].values.astype(int)
    fold = make_folds(y, d["scaffold"].tolist())

    print(f"Compounds   : {len(d)}   actives {int(y.sum())}   "
          f"decoys {int((1 - y).sum())}")
    print(f"Features    : 6 physicochemical descriptors only "
          f"(MW, logP, TPSA, RotB, HBD, HBA)")
    print(f"Folds       : {len(np.unique(fold[fold >= 0]))}, "
          f"leave-one-scaffold-out, identical scheme to fp_classify.py")

    p = loso(X, y, fold)
    ok = ~np.isnan(p)
    auc, efs, bedroc = metrics(-p[ok], y[ok])
    print(f"\n=== property-only baseline, out-of-fold ===")
    for pct, (found, k, ef) in efs.items():
        print(f"  EF at {pct:>2}% : {found:>3} in top {k:<5} = {ef:>5.1f}x")
    print(f"  AUC-ROC       : {auc:.3f}")
    print(f"  BEDROC (a=20) : {bedroc:.3f}")

    perm = []
    for seed in range(5):
        rng = np.random.default_rng(seed)
        ysh = y.copy()
        rng.shuffle(ysh)
        pp = loso(X, ysh, fold, seed=seed)
        o = ~np.isnan(pp)
        perm.append(metrics(-pp[o], ysh[o])[0])
    print(f"  permutation   : {np.mean(perm):.3f} +/- {np.std(perm):.3f}")

    out = d[["name", "label"]].copy()
    out["probability"] = p
    out.to_csv("baseline_properties_oof.csv", index=False)
    print("\nWrote baseline_properties_oof.csv")

    print("\n" + "=" * 62)
    print("HOW TO READ THIS - decide before you look, then follow it:")
    print("  baseline AUC < 0.65  -> the interaction features carry the signal.")
    print("                          Report the baseline as a control row.")
    print("  baseline 0.65 - 0.75 -> a real share of the 0.822 is molecular")
    print("                          size and polarity. Say so explicitly and")
    print("                          run delong.py to test the difference.")
    print("  baseline > 0.75      -> the fingerprint advantage is small. The")
    print("                          manuscript must be rewritten around that,")
    print("                          not around 0.822. Do not bury this.")
    print("Whatever it says, this number goes in the paper.")
    print("=" * 62)


if __name__ == "__main__":
    main()
