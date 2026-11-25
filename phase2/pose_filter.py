#!/usr/bin/env python3
"""
Test the pose filter on the full benchmark, not on 7 drugs.

WHY
---
A separate analysis on 7 marketed PDE4 drugs reported AUC 0.864 for a
"pose-filtered score": keep compounds whose top pose hydrogen-bonds the
invariant Gln and sits in the hydrophobic clamp. All 7 actives passed and
23% of decoys passed.

Two reasons that number cannot be taken at face value. A binary filter with
those pass rates gives AUC 0.885 by itself, so the Vina score adds nothing -
the filter is the whole signal. And with 6 to 7 actives, several of them
protonation variants of the same molecule, the confidence interval runs from
about 0.68 to 1.0, which settles nothing.

This applies the SAME filter to the 40-active, 1,197-decoy benchmark. Same
criteria, six times the actives, properly property-matched decoys, 40
distinct scaffolds. If the filter holds up here it is a real finding. If it
does not, the 0.864 was the small sample.

It reports three things so the filter and the score can be told apart:
  1. the filter alone, as a binary classifier
  2. the score alone
  3. the score within the passing set only, which is what "pose-filtered
     score" actually means

USAGE
-----
    python pose_filter.py fingerprints.csv rescore_results.csv
"""

import os
import sys
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, hypergeom


def metrics(scores, labels, alpha=20.0):
    order = np.argsort(scores)
    y = np.asarray(labels)[order]
    n = len(y); n_act = int(y.sum()); n_dec = n - n_act
    if n_act == 0 or n_dec == 0:
        raise ValueError("need both classes")
    ranks = np.arange(1, n + 1)
    auc = 1 - ((ranks[y == 1].sum() - n_act * (n_act + 1) / 2) / (n_act * n_dec))
    efs = {}
    for pct in (1, 5, 10):
        k = max(1, int(round(n * pct / 100)))
        found = int(y[:k].sum())
        exp = n_act * k / n
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
    """Hanley & McNeil 1982. Wide intervals at small n are the point."""
    Q1 = A / (2 - A); Q2 = 2 * A * A / (1 + A)
    return float(np.sqrt((A * (1 - A) + (n1 - 1) * (Q1 - A * A)
                          + (n0 - 1) * (Q2 - A * A)) / (n1 * n0)))


def show(tag, scores, labels):
    auc, efs, bedroc = metrics(scores, labels)
    n1 = int(np.sum(labels)); n0 = len(labels) - n1
    s = se_auc(auc, n1, n0)
    print(f"\n=== {tag} ===")
    print(f"  {n1} actives, {n0} decoys")
    print(f"  AUC      {auc:.3f}   95% CI "
          f"[{max(0,auc-1.96*s):.3f}, {min(1,auc+1.96*s):.3f}]")
    for pct, (found, k, ef) in efs.items():
        print(f"  EF {pct:>2}%   {ef:>5.1f}x   ({found} of {n1} in top {k})")
    print(f"  BEDROC   {bedroc:.3f}")
    return auc


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    for p in sys.argv[1:3]:
        if not os.path.exists(p):
            sys.exit(f"ERROR: not found: {p}")

    fp = pd.read_csv(sys.argv[1])
    sc = pd.read_csv(sys.argv[2])
    df = fp.merge(sc, on="name", suffixes=("", "_s"))
    lab = "label" if "label" in df.columns else "label_s"
    for c in ("gln_hbond", "n_clamp"):
        if c not in df.columns:
            sys.exit(f"ERROR: no {c} column. Re-run fingerprint.py.")

    df["gln_hbond"] = pd.to_numeric(df["gln_hbond"], errors="coerce").fillna(0)
    df["passes"] = ((df["gln_hbond"] == 1) & (df["n_clamp"] >= 1)).astype(int)

    y = df[lab].values.astype(int)
    n_act = int(y.sum())
    print(f"Compounds {len(df)}   actives {n_act}   decoys {len(df)-n_act}")

    # ---------- pass rates, the thing actually being claimed ----------
    pa = int(df[(df[lab] == 1)]["passes"].sum())
    pd_ = int(df[(df[lab] == 0)]["passes"].sum())
    na, nd = n_act, len(df) - n_act
    print(f"\nPass rates for 'H-bonds Gln AND contacts the clamp':")
    print(f"  actives {pa}/{na} = {100*pa/na:.0f}%        "
          f"(the 7-drug set reported 100%)")
    print(f"  decoys  {pd_}/{nd} = {100*pd_/nd:.0f}%        "
          f"(the 7-drug set reported 23%)")
    odds, pf = fisher_exact([[pa, na - pa], [pd_, nd - pd_]],
                            alternative="greater")
    print(f"  Fisher exact p = {pf:.2e}   odds ratio {odds:.2f}")
    if pa + pd_:
        enr = (pa / (pa + pd_)) / (na / (na + nd))
        print(f"  an active is {enr:.2f}x enriched in the passing set")

    # ---------- 1. the filter on its own ----------
    # ties broken at random so the AUC is not inflated by file order
    rng = np.random.default_rng(0)
    show("the filter alone, as a binary classifier",
         -df["passes"].values + rng.uniform(0, 1e-6, len(df)), y)

    # ---------- 2. the score on its own ----------
    for col in ("vina_score", "vinardo_score"):
        if col in df.columns:
            s = df.dropna(subset=[col])
            show(f"{col} alone", s[col].values, s[lab].values.astype(int))

    # ---------- 3. the score inside the passing set ----------
    # This is what "pose-filtered score" means: non-passing compounds are
    # pushed to the bottom, passing ones ordered by score.
    for col in ("vina_score", "vinardo_score"):
        if col not in df.columns:
            continue
        s = df.dropna(subset=[col]).copy()
        worst = s[col].max() + 1.0
        # Every non-passing compound gets the SAME placeholder, so they form
        # one large tie block. numpy's argsort breaks ties by row order, and
        # if actives happen to be listed before decoys the AUC inflates for
        # free - no signal required. Jitter breaks the ties at random, which
        # is the honest treatment of a tie.
        base = np.where(s["passes"] == 1, s[col], worst)
        jit = np.random.default_rng(0).uniform(0, 1e-6, len(s))
        show(f"filter + {col}", base + jit, s[lab].values.astype(int))

    # ---------- does the score add anything beyond the filter? ----------
    print("\n" + "=" * 62)
    print("DOES THE SCORE ADD ANYTHING?")
    print("=" * 62)
    for col in ("vina_score", "vinardo_score"):
        if col not in df.columns:
            continue
        s = df[(df["passes"] == 1) & df[col].notna()]
        a = int(s[lab].sum()); d = len(s) - a
        if a >= 3 and d >= 3:
            auc, efs, _ = metrics(s[col].values, s[lab].values.astype(int))
            se = se_auc(auc, a, d)
            print(f"  {col} WITHIN the passing set only: AUC {auc:.3f} "
                  f"[{max(0,auc-1.96*se):.3f}, {min(1,auc+1.96*se):.3f}]")
            print(f"    {a} actives, {d} decoys. If this interval contains")
            print(f"    0.5, the score is not ranking inside the passing set")
            print(f"    and the filter is doing all the work.")
        else:
            print(f"  {col}: too few compounds pass to test "
                  f"({a} actives, {d} decoys)")


if __name__ == "__main__":
    main()
