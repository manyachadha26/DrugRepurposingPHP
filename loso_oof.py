#!/usr/bin/env python3
"""
Run the interaction-fingerprint classifier with leave-one-scaffold-out and
SAVE the per-compound out-of-fold probabilities.

Same model and same fold construction as fp_classify.py - this adds nothing
to the method, it just writes the predictions to disk so they can be compared
statistically (delong.py) instead of by eyeballing confidence intervals.

    python loso_oof.py fingerprints_150.csv actives_150.smi decoys_150.smi out.csv

Writes out.csv with columns: name, label, probability.
Also prints AUC / EF / BEDROC so you can confirm it matches fp_classify.py.
If the AUC it prints does NOT match what fp_classify.py printed on the same
inputs, stop and work out why before using either number.
"""
import sys

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
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


def main():
    if len(sys.argv) < 5:
        sys.exit(__doc__)
    fp_csv, act_smi, dec_smi, out_csv = sys.argv[1:5]
    with_hac = "--with-hac" in sys.argv

    fp = pd.read_csv(fp_csv)
    smiles = {}
    smiles.update(read_smi(act_smi))
    smiles.update(read_smi(dec_smi))
    fp["scaffold"] = [scaffold_of(smiles[n]) if n in smiles else None
                      for n in fp["name"]]
    missing = int(fp["scaffold"].isna().sum())
    if missing:
        print(f"WARNING: {missing} of {len(fp)} compounds had no SMILES match "
              f"and were dropped.")
        print("         If this number is large you have passed the wrong "
              ".smi files for this fingerprint set.")
    fp = fp.dropna(subset=["scaffold"]).reset_index(drop=True)

    drop = {"name", "label", "scaffold"}
    if not with_hac:
        drop.add("hac")
    Xdf = fp[[c for c in fp.columns if c not in drop]].apply(pd.to_numeric,
                                                             errors="coerce")
    Xdf = Xdf.fillna(Xdf.median(numeric_only=True)).fillna(0.0)
    Xdf = Xdf[[c for c in Xdf.columns if Xdf[c].std() > 0]]
    X = Xdf.values.astype(float)
    y = fp["label"].values.astype(int)

    fold = make_folds(y, fp["scaffold"].tolist())
    proba = np.full(len(y), np.nan)
    for f in np.unique(fold[fold >= 0]):
        te = fold == f
        tr = ~te
        if y[tr].sum() < 2 or (1 - y[tr]).sum() < 2:
            continue
        clf = make_pipeline(StandardScaler(),
                            LogisticRegression(max_iter=2000,
                                               class_weight="balanced",
                                               random_state=0))
        clf.fit(X[tr], y[tr])
        proba[te] = clf.predict_proba(X[te])[:, 1]

    ok = ~np.isnan(proba)
    auc, efs, bedroc = metrics(-proba[ok], y[ok])
    print(f"compounds {len(fp)}  actives {int(y.sum())}  "
          f"decoys {int((1 - y).sum())}  features {X.shape[1]}")
    print(f"AUC {auc:.3f}   EF1% {efs[1][2]:.1f}x   BEDROC {bedroc:.3f}")

    pd.DataFrame({"name": fp["name"], "label": y,
                  "probability": proba}).to_csv(out_csv, index=False)
    print(f"Wrote {out_csv}")


if __name__ == "__main__":
    main()
