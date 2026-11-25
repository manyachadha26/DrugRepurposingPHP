#!/usr/bin/env python3
"""
DeLong's test for two correlated ROC curves.

Use this instead of saying "the confidence intervals overlap". When two models
or two receptor preparations are scored on THE SAME compounds, their AUCs are
correlated, and comparing them through independent intervals is the wrong test -
it is far too conservative and will hide a real difference, or let you imply a
difference that is not there.

    python delong.py a.csv b.csv ["label for A"] ["label for B"]

Each CSV needs columns: name, label (1 active / 0 decoy), probability.
Higher probability must mean "more likely active" in both files.
Only compounds present in BOTH files are used, matched on name.

Reference: DeLong, DeLong & Clarke-Pearson, Biometrics 44:837 (1988);
fast algorithm from Sun & Xu, IEEE Signal Process Lett 21:1389 (2014).
"""
import sys

import numpy as np
import pandas as pd
from scipy import stats


def _midrank(x):
    J = np.argsort(x)
    Z = x[J]
    N = len(x)
    T = np.zeros(N, dtype=float)
    i = 0
    while i < N:
        j = i
        while j < N and Z[j] == Z[i]:
            j += 1
        T[i:j] = 0.5 * (i + j - 1) + 1
        i = j
    T2 = np.empty(N, dtype=float)
    T2[J] = T
    return T2


def fast_delong(predictions_sorted_transposed, label_1_count):
    """predictions_sorted_transposed: (k, n) with all positives first."""
    m = label_1_count
    n = predictions_sorted_transposed.shape[1] - m
    positive = predictions_sorted_transposed[:, :m]
    negative = predictions_sorted_transposed[:, m:]
    k = predictions_sorted_transposed.shape[0]

    tx = np.empty([k, m], dtype=float)
    ty = np.empty([k, n], dtype=float)
    tz = np.empty([k, m + n], dtype=float)
    for r in range(k):
        tx[r, :] = _midrank(positive[r, :])
        ty[r, :] = _midrank(negative[r, :])
        tz[r, :] = _midrank(predictions_sorted_transposed[r, :])
    aucs = tz[:, :m].sum(axis=1) / m / n - (m + 1.0) / 2.0 / n
    v01 = (tz[:, :m] - tx[:, :]) / n
    v10 = 1.0 - (tz[:, m:] - ty[:, :]) / m
    sx = np.cov(v01)
    sy = np.cov(v10)
    delongcov = sx / m + sy / n
    return aucs, np.atleast_2d(delongcov)


def delong_test(y, score_a, score_b):
    order = np.argsort(-y, kind="mergesort")
    y = y[order]
    preds = np.vstack((score_a[order], score_b[order]))
    m = int(y.sum())
    aucs, cov = fast_delong(preds, m)
    var = cov[0, 0] + cov[1, 1] - 2 * cov[0, 1]
    if var <= 0:
        return aucs, 1.0, 0.0, cov
    z = (aucs[0] - aucs[1]) / np.sqrt(var)
    p = 2 * stats.norm.sf(abs(z))
    return aucs, p, z, cov


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    fa, fb = sys.argv[1], sys.argv[2]
    la = sys.argv[3] if len(sys.argv) > 3 else fa
    lb = sys.argv[4] if len(sys.argv) > 4 else fb

    A = pd.read_csv(fa)[["name", "label", "probability"]]
    B = pd.read_csv(fb)[["name", "label", "probability"]]
    M = A.merge(B, on="name", suffixes=("_a", "_b")).dropna()
    if len(M) == 0:
        sys.exit("No compounds in common. Are these the same benchmark?")
    if not (M.label_a == M.label_b).all():
        sys.exit("Labels disagree between the two files for some compounds. "
                 "Fix that before testing.")

    y = M.label_a.values.astype(int)
    sa = M.probability_a.values.astype(float)
    sb = M.probability_b.values.astype(float)

    print(f"matched compounds {len(M)}   actives {int(y.sum())}   "
          f"decoys {int((1 - y).sum())}")
    dropped_a, dropped_b = len(A) - len(M), len(B) - len(M)
    if dropped_a or dropped_b:
        print(f"  (dropped {dropped_a} from A and {dropped_b} from B "
              f"as unmatched or missing)")

    aucs, p, z, cov = delong_test(y, sa, sb)
    se = np.sqrt(np.diag(cov))
    for lab, a, s in ((la, aucs[0], se[0]), (lb, aucs[1], se[1])):
        print(f"  {lab:<34} AUC {a:.4f}  95% CI "
              f"[{a - 1.96 * s:.4f}, {a + 1.96 * s:.4f}]")
    diff = aucs[0] - aucs[1]
    sd = np.sqrt(cov[0, 0] + cov[1, 1] - 2 * cov[0, 1])
    print(f"\n  difference (A - B) {diff:+.4f}  95% CI "
          f"[{diff - 1.96 * sd:+.4f}, {diff + 1.96 * sd:+.4f}]")
    print(f"  DeLong z = {z:.3f},  p = {p:.4g}")

    print("\n" + "=" * 62)
    if p < 0.05:
        better = la if diff > 0 else lb
        print(f"Significant at p < 0.05: {better} is the better ranker on")
        print("these compounds. Report the difference and its interval.")
    else:
        print("NOT significant. Do not describe either as better. The correct")
        print("sentence is that the two are statistically indistinguishable on")
        print(f"this benchmark (difference {diff:+.3f}, p = {p:.3g}) - which")
        print("is a result in itself, not a failed experiment.")
    print("Report the p-value whichever way it goes. A non-significant")
    print("difference that gets described as an improvement is the single")
    print("easiest thing for a referee to catch.")
    print("=" * 62)


if __name__ == "__main__":
    main()
