#!/usr/bin/env python3
"""
Can the interaction fingerprint separate actives from decoys, where the
docking score could not?

THE CIRCULARITY TRAP, AND HOW THIS AVOIDS IT
--------------------------------------------
Fitting a model on 40 actives and then scoring those same 40 would produce
a beautiful number that means nothing. Worse, a random split still leaks:
two analogues of the same chemotype land on both sides, so the model
recognises the scaffold rather than the interaction pattern.

So we use LEAVE-ONE-SCAFFOLD-OUT. Every compound gets its Bemis-Murcko
scaffold. For each scaffold that contains an active, we hold out EVERY
compound with that scaffold (actives and decoys alike), train on the rest,
and predict the held-out set. Each prediction is therefore made by a model
that never saw that chemotype. We pool those out-of-fold predictions and
score them with exactly the metrics used on the docking scores, so the
comparison is like for like.

We also run a LABEL PERMUTATION check. Shuffle the labels and repeat the
whole procedure: a sound pipeline collapses to AUC near 0.5. If it does
not, there is leakage and the real number cannot be trusted either.

THE PASS BAR, fixed before looking
----------------------------------
  out-of-fold AUC >= 0.70 AND EF at 1% >= 5x   -> Path 1 succeeds, expand
                                                   to 150 actives, confirm
  AUC 0.60 - 0.70                              -> inconclusive, 40 is too
                                                   few; expand, then re-run
  AUC < 0.60                                   -> Path 1 fails, go to Path 2

USAGE
-----
    python fp_classify.py fingerprints.csv actives_chembl.smi \\
           decoys_chembl.smi [rescore_results.csv] [--with-hac]

--with-hac adds heavy-atom count as a feature. Off by default: raw size is
what confounded Gate 2 round 1, and we do not want it back by the side door.
Run it both ways and report both.
"""

import os
import sys
import numpy as np
import pandas as pd

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem.Scaffolds import MurckoScaffold
    RDLogger.DisableLog("rdApp.*")
except ImportError as exc:
    sys.exit(f"ERROR: rdkit import failed ({exc}). Activate your env.")

from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline

ALPHA = 20.0


def metrics(scores, labels, alpha=ALPHA):
    """Lower score = better, same convention as the docking energies."""
    order = np.argsort(scores)
    y = np.asarray(labels)[order]
    n = len(y)
    n_act = int(y.sum())
    n_dec = n - n_act
    if n_act == 0 or n_dec == 0:
        raise ValueError("need both classes")

    ranks = np.arange(1, n + 1)
    auc = 1 - ((ranks[y == 1].sum() - n_act * (n_act + 1) / 2)
               / (n_act * n_dec))

    efs = {}
    for pct in (1, 5, 10):
        k = max(1, int(round(n * pct / 100)))
        found = int(y[:k].sum())
        expected = n_act * k / n
        efs[pct] = (found, k, found / expected if expected else 0.0)

    ri = np.where(y == 1)[0] + 1
    ra = n_act / n
    s = np.sum(np.exp(-alpha * ri / n)) / n_act
    s_rand = (1.0 / n) * (1 - np.exp(-alpha)) / (np.exp(alpha / n) - 1)
    rie = s / s_rand
    rie_max = (1 - np.exp(-alpha * ra)) / (ra * (1 - np.exp(-alpha)))
    rie_min = (1 - np.exp(alpha * ra)) / (ra * (1 - np.exp(alpha)))
    bedroc = float(np.clip((rie - rie_min) / (rie_max - rie_min), 0, 1))
    return auc, efs, bedroc


def read_smi(path):
    """{name: smiles}, tolerant of either column order."""
    out = {}
    with open(path) as fh:
        for line in fh:
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
    try:
        s = MurckoScaffold.GetScaffoldForMol(m)
        sm = Chem.MolToSmiles(s)
        # a molecule with no ring system yields an empty scaffold; fall back
        # to the whole molecule so it forms its own group rather than
        # collapsing every acyclic compound into one giant fold
        return sm if sm else Chem.MolToSmiles(m)
    except Exception:
        return None


def make_folds(y, groups):
    """
    Assign every compound to exactly one fold, without splitting a scaffold.

    One fold per active-bearing scaffold. A scaffold that contains an active
    defines its own fold. Every other scaffold - the decoys that share no
    chemotype with any active, which is most of them - is handed round-robin
    to those folds, kept whole.

    This matters: an earlier version only ever predicted compounds sitting on
    an active's scaffold, so nearly all decoys were never scored and there was
    no negative class left to compute an AUC against. Here each compound is
    predicted exactly once, by a model that never saw its scaffold.
    """
    by_group = {}
    for i, g in enumerate(groups):
        by_group.setdefault(g, []).append(i)

    act_groups = sorted({g for g, lab in zip(groups, y) if lab == 1})
    if not act_groups:
        raise ValueError("no actives have a scaffold; cannot build folds")

    fold = np.full(len(y), -1, dtype=int)
    for f, g in enumerate(act_groups):
        for i in by_group[g]:
            fold[i] = f

    # the rest, biggest groups first so the folds stay balanced in size
    others = sorted((g for g in by_group if g not in set(act_groups)),
                    key=lambda g: -len(by_group[g]))
    sizes = np.array([(fold == f).sum() for f in range(len(act_groups))],
                     dtype=float)
    for g in others:
        f = int(np.argmin(sizes))
        for i in by_group[g]:
            fold[i] = f
        sizes[f] += len(by_group[g])
    return fold


def loso(X, y, fold, seed=0):
    """
    Out-of-fold probabilities using the fold assignment from make_folds.
    nan for any compound whose fold could not be fitted.
    """
    proba = np.full(len(y), np.nan)
    for f in np.unique(fold):
        if f < 0:
            continue
        test = fold == f
        train = ~test
        if y[train].sum() < 2 or (1 - y[train]).sum() < 2:
            continue
        clf = make_pipeline(
            StandardScaler(),
            LogisticRegression(max_iter=2000, class_weight="balanced",
                               random_state=seed),
        )
        clf.fit(X[train], y[train])
        proba[test] = clf.predict_proba(X[test])[:, 1]
    return proba


def report(tag, proba, y):
    ok = ~np.isnan(proba)
    # the classifier outputs P(active); metrics expect lower = better
    auc, efs, bedroc = metrics(-proba[ok], y[ok])
    print(f"\n=== {tag} ===")
    print(f"  scored {int(ok.sum())}   actives {int(y[ok].sum())}   "
          f"decoys {int((1 - y[ok]).sum())}")
    for pct, (found, k, ef) in efs.items():
        print(f"  EF at {pct:>2}% : {found:>3} of {int(y[ok].sum())} "
              f"in top {k:<5} = {ef:>5.1f}x")
    print(f"  AUC-ROC       : {auc:.3f}")
    print(f"  BEDROC (a=20) : {bedroc:.3f}")
    return auc, efs[1][2], bedroc


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    with_hac = "--with-hac" in sys.argv
    if len(args) < 3:
        print(__doc__)
        sys.exit(1)

    fp = pd.read_csv(args[0])
    smiles = {}
    for p in args[1:3]:
        smiles.update(read_smi(p))
    baseline_csv = args[3] if len(args) > 3 else None

    # ---------- scaffolds ----------
    fp["scaffold"] = [scaffold_of(smiles[n]) if n in smiles else None
                      for n in fp["name"]]
    missing = int(fp["scaffold"].isna().sum())
    if missing:
        print(f"WARNING: no scaffold for {missing} of {len(fp)} compounds "
              f"(no SMILES match). They are dropped.")
    fp = fp.dropna(subset=["scaffold"]).reset_index(drop=True)

    drop = {"name", "label", "scaffold"}
    if not with_hac:
        drop.add("hac")
    feat_cols = [c for c in fp.columns if c not in drop]
    # drop constant and all-blank columns: they carry nothing and make the
    # scaler divide by zero
    Xdf = fp[feat_cols].apply(pd.to_numeric, errors="coerce")
    Xdf = Xdf.fillna(Xdf.median(numeric_only=True)).fillna(0.0)
    keep = [c for c in Xdf.columns if Xdf[c].std() > 0]
    dropped = sorted(set(Xdf.columns) - set(keep))
    Xdf = Xdf[keep]

    X = Xdf.values.astype(float)
    y = fp["label"].values.astype(int)
    groups = fp["scaffold"].tolist()

    n_act = int(y.sum())
    n_sc = len({g for g, lab in zip(groups, y) if lab == 1})
    print(f"Compounds   : {len(fp)}   actives {n_act}   "
          f"decoys {len(fp) - n_act}")
    print(f"Features    : {len(keep)}"
          + (f"  ({len(dropped)} constant columns dropped)" if dropped else ""))
    print(f"Folds       : {n_sc} active-bearing scaffolds, "
          f"leave-one-scaffold-out")
    print(f"Heavy-atom count as a feature: {'YES' if with_hac else 'no'}")
    if n_act < 25:
        print(f"\nNOTE: {n_act} actives is few. Expect wide error bars; a "
              f"borderline result here means expand, not conclude.")

    # ---------- the real run ----------
    fold = make_folds(y, groups)
    nfold = len(np.unique(fold[fold >= 0]))
    sizes = [int((fold == f).sum()) for f in range(nfold)]
    print(f"Fold sizes  : {min(sizes)} to {max(sizes)} compounds, "
          f"{nfold} folds, every compound predicted once")
    proba = loso(X, y, fold)
    auc, ef1, bedroc = report("interaction fingerprint, out-of-fold",
                              proba, y)

    # ---------- permutation sanity check ----------
    print("\n=== label permutation check ===")
    print("  A sound pipeline collapses to AUC near 0.5 on shuffled labels.")
    perm_aucs = []
    for seed in range(5):
        rng = np.random.default_rng(seed)
        ysh = y.copy()
        rng.shuffle(ysh)
        try:
            p = loso(X, ysh, fold, seed=seed)
            ok = ~np.isnan(p)
            perm_aucs.append(metrics(-p[ok], ysh[ok])[0])
        except ValueError:
            continue
    if perm_aucs:
        m, s = float(np.mean(perm_aucs)), float(np.std(perm_aucs))
        print(f"  shuffled AUC: {m:.3f} +/- {s:.3f} over {len(perm_aucs)} runs")
        if m > 0.60:
            print("  *** LEAKAGE. The pipeline scores well on random labels,")
            print("  *** so the real number above is not trustworthy. Stop")
            print("  *** and work out where information is crossing folds.")

    # ---------- baseline comparison ----------
    if baseline_csv and os.path.exists(baseline_csv):
        base = pd.read_csv(baseline_csv)
        sub = fp[["name", "label"]].merge(base, on="name", how="inner",
                                          suffixes=("", "_b"))
        print("\n=== the docking scores, on the same compounds ===")
        for col in ("vina_score", "vinardo_score"):
            if col in sub.columns:
                s = sub.dropna(subset=[col])
                a, e, b = metrics(s[col].values, s["label"].values)
                print(f"  {col:<15} AUC {a:.3f}   EF1% {e[1][2]:>5.1f}x   "
                      f"BEDROC {b:.3f}")

    # ---------- interpretation, clearly labelled ----------
    clf = make_pipeline(
        StandardScaler(),
        LogisticRegression(max_iter=2000, class_weight="balanced"),
    ).fit(X, y)
    coef = clf[-1].coef_[0]
    top = sorted(zip(keep, coef), key=lambda t: -abs(t[1]))[:10]
    print("\n=== which contacts the model leans on ===")
    print("  Fit on ALL data, for interpretation only. These numbers are")
    print("  NOT evidence of performance - the out-of-fold AUC above is.")
    for nm, c in top:
        print(f"  {nm:<18} {c:+.3f}   "
              f"{'favours actives' if c > 0 else 'favours decoys'}")

    # ---------- verdict ----------
    print("\n" + "=" * 62)
    # The highest EF reachable at any cut-off is n_total / n_actives. On a
    # small set that ceiling can sit below the 5x bar, so applying the bar
    # there would fail a perfect result. Check before judging on it.
    ok = ~np.isnan(proba)
    ef_ceiling = int(ok.sum()) / max(1, int(y[ok].sum()))
    if ef_ceiling < 5.0:
        print(f"EF at 1% cannot exceed {ef_ceiling:.1f}x on a set this size,")
        print("so the 5x bar does not apply. Judging on AUC alone.")
        if auc >= 0.70:
            print("\nPASS on AUC. Confirm on the full-size benchmark before")
            print("relying on it.")
        else:
            print(f"\nAUC {auc:.3f} is below 0.70. Not usable as it stands.")
    elif auc >= 0.70 and ef1 >= 5.0:
        print("PASS. The interaction pattern discriminates where the score")
        print("did not. Expand to 150 actives and confirm, then use this to")
        print("rank the antidiabetic screen.")
    elif auc >= 0.60:
        print("INCONCLUSIVE. Better than the score, short of the bar, and")
        print(f"{n_act} actives cannot settle it. Expand to 150 and re-run")
        print("before drawing any conclusion.")
    else:
        print("FAIL. The interaction pattern does not separate these")
        print("compounds either. Do not expand. Move to Path 2 (MM-GBSA),")
        print("one week, then write Path 3 if that fails too.")
    print("=" * 62)


if __name__ == "__main__":
    main()
