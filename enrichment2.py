#!/usr/bin/env python3
"""
Compare two scoring functions on the SAME Gate 2 benchmark.

WHAT IT DOES
------------
Reads rescore_results.csv (name, label, vina_score, vinardo_score) and for
each score column reports:

  AUC-ROC            0.5 = coin flip, 1.0 = perfect. We want >= 0.70.
  EF at 1/5/10%      "how many times better than random in the top slice".
  BEDROC alpha=20    like AUC but only rewards hits near the very top,
                     which is what matters in a real screen.

It prints the two side by side so the comparison is honest: identical
molecules, identical decoys, only the scoring function changed.

PASS BAR (fixed in advance, do not move it after seeing the numbers)
--------------------------------------------------------------------
  AUC >= 0.70  AND  EF1% >= 5x    -> the scoring function is usable for
                                      ranking; proceed with it.
  AUC 0.60-0.70                   -> weak. Usable only to triage which
                                      compounds enter the pocket, never
                                      to claim predicted potency.
  AUC < 0.60                      -> not usable for ranking at all.

USAGE
-----
    python enrichment2.py rescore_results.csv
"""

import sys
import numpy as np
import pandas as pd


def metrics(scores, labels, alpha=20.0):
    """
    scores: lower = better (docking energies are negative)
    labels: 1 active, 0 decoy
    """
    order = np.argsort(scores)            # most negative first
    y = np.asarray(labels)[order]
    n = len(y)
    n_act = int(y.sum())
    n_dec = n - n_act
    if n_act == 0 or n_dec == 0:
        raise ValueError("need both actives and decoys")

    # --- AUC-ROC via the rank-sum identity (exact, no binning) ---
    ranks = np.arange(1, n + 1)
    rank_sum = ranks[y == 1].sum()
    auc = (rank_sum - n_act * (n_act + 1) / 2) / (n_act * n_dec)
    auc = 1 - auc                         # because low score = good

    # --- enrichment factor at cut-offs ---
    efs = {}
    for pct in (1, 5, 10):
        k = max(1, int(round(n * pct / 100)))
        found = int(y[:k].sum())
        expected = n_act * k / n
        efs[pct] = (found, k, found / expected if expected else 0.0)

    # --- BEDROC, Truchon & Bayly 2007, closed form ---
    ri = np.where(y == 1)[0] + 1          # 1-based ranks of the actives
    ra = n_act / n
    # RIE = observed exponential sum / the sum expected at random
    s = np.sum(np.exp(-alpha * ri / n)) / n_act
    s_rand = (1.0 / n) * (1 - np.exp(-alpha)) / (np.exp(alpha / n) - 1)
    rie = s / s_rand
    # map RIE onto 0-1
    rie_max = (1 - np.exp(-alpha * ra)) / (ra * (1 - np.exp(-alpha)))
    rie_min = (1 - np.exp(alpha * ra)) / (ra * (1 - np.exp(alpha)))
    bedroc = float(np.clip((rie - rie_min) / (rie_max - rie_min), 0.0, 1.0))

    return auc, efs, bedroc, n_act, n_dec


def verdict(auc, ef1):
    if auc >= 0.70 and ef1 >= 5.0:
        return ("USABLE - this scoring function ranks PDE4D binders. "
                "Proceed with it for the screen.")
    if auc >= 0.60:
        return ("WEAK - use only to triage pocket occupancy, never to "
                "claim predicted potency. Must be stated in the paper.")
    return "NOT USABLE for ranking."


def report(df, col):
    sub = df.dropna(subset=[col])
    auc, efs, bedroc, n_act, n_dec = metrics(sub[col].values,
                                             sub["label"].values)
    print(f"\n=== {col} ===")
    print(f"  scored {len(sub)}   actives {n_act}   decoys {n_dec}")
    for pct, (found, k, ef) in efs.items():
        print(f"  EF at {pct:>2}% : {found:>3} of {n_act} actives "
              f"in top {k:<5} = {ef:>5.1f}x")
    print(f"  AUC-ROC          : {auc:.3f}")
    print(f"  BEDROC (a=20)    : {bedroc:.3f}")
    print(f"  -> {verdict(auc, efs[1][2])}")
    return auc, efs[1][2], bedroc


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    df = pd.read_csv(sys.argv[1])
    if "label" not in df.columns:
        sys.exit("ERROR: csv needs a 'label' column (1 active, 0 decoy)")

    cols = [c for c in ("vina_score", "vinardo_score") if c in df.columns]
    if not cols:
        sys.exit("ERROR: no score columns found")

    print("Gate 2 re-test: same molecules, same decoys, "
          "different scoring function.")
    out = {c: report(df, c) for c in cols}

    if len(out) == 2:
        a_vina = out["vina_score"][0]
        a_vinardo = out["vinardo_score"][0]
        print("\n" + "=" * 60)
        print(f"Vina AUC {a_vina:.3f}  ->  Vinardo AUC {a_vinardo:.3f}  "
              f"({a_vinardo - a_vina:+.3f})")
        if a_vinardo >= 0.70:
            print("Vinardo clears the bar. Use Vinardo scores for the screen,")
            print("and report both numbers in the paper.")
        elif a_vinardo > a_vina + 0.05:
            print("Vinardo is better but still under 0.70. Report both, and")
            print("treat the screen as pocket-occupancy triage, not ranking.")
        else:
            print("No real improvement. Both scoring functions fail to rank")
            print("this target. That is the finding - report it plainly and")
            print("let MD / MM-GBSA do the affinity work in Phase 4.")
        print("=" * 60)


if __name__ == "__main__":
    main()
