#!/usr/bin/env python3
"""
Do the hits bind the same way in a second crystal structure?

WHY
---
Every conclusion so far rests on one receptor conformation, 1XOQ. A pocket
is not rigid: a compound that fits one crystal structure can fail in another
of the same protein. That is the most obvious thing a reviewer will ask, and
the cheapest way to be wrong.

So the 41 screen compounds were docked again against 1XMU, an independent
PDE4D structure that passed the same redocking gate (0.927 A). 7CBQ was
excluded for failing that gate at 3.511 A, on both the dry and the wet
preparation.

WHY RAW FEATURES AND NOT THE TRAINED MODEL
------------------------------------------
The classifier's per-residue features are named by residue number, and the
two structures number their residues differently - 1XOQ's Gln369, Phe340 and
Phe372 are Gln443, Phe414 and Phe446 in 1XMU, a constant offset of 74. The
columns cannot be matched, and training a second classifier would mean
docking all 2,392 benchmark compounds again.

The pharmacophore features carry across unchanged, because the script locates
the invariant glutamine and the clamp phenylalanines GEOMETRICALLY in each
structure rather than by number. So the question becomes the right one
anyway: does the compound make the same contacts in both conformations?

THE CONSENSUS TEST, fixed before looking
----------------------------------------
A compound is CONSISTENT when, in BOTH structures, it:
  - hydrogen-bonds the invariant glutamine (gln_hbond = 1), and
  - contacts at least one clamp phenylalanine (n_clamp >= 1)

Anything that satisfies this in only one structure is conformation-dependent
and must be reported as such, not quietly promoted on its better result.

USAGE
-----
    python compare_structures.py screen_results.csv library_fp.csv \\
                                 library_fp_1XMU.csv
"""

import os
import sys
import numpy as np
import pandas as pd

KEY = ["d_gln", "gln_hbond", "n_clamp", "d_metal", "f_buried"]


def main():
    if len(sys.argv) < 4:
        print(__doc__)
        sys.exit(1)
    for p in sys.argv[1:4]:
        if not os.path.exists(p):
            sys.exit(f"ERROR: not found: {p}")

    rank = pd.read_csv(sys.argv[1])
    a = pd.read_csv(sys.argv[2])     # 1XOQ
    b = pd.read_csv(sys.argv[3])     # 1XMU

    miss = [c for c in KEY if c not in a.columns or c not in b.columns]
    if miss:
        sys.exit(f"ERROR: missing columns in one of the fingerprint "
                 f"files: {miss}")

    df = (rank[["rank", "name", "probability", "pct_of_actives_below",
                "is_control"]]
          if "is_control" in rank.columns
          else rank[["rank", "name", "probability",
                     "pct_of_actives_below"]]).copy()
    df = df.merge(a[["name"] + KEY], on="name", how="left")
    df = df.merge(b[["name"] + KEY], on="name", how="left",
                  suffixes=("_oq", "_mu"))

    for c in ("gln_hbond_oq", "gln_hbond_mu", "n_clamp_oq", "n_clamp_mu"):
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0)

    df["ok_oq"] = (df.gln_hbond_oq == 1) & (df.n_clamp_oq >= 1)
    df["ok_mu"] = (df.gln_hbond_mu == 1) & (df.n_clamp_mu >= 1)
    df["consistent"] = df.ok_oq & df.ok_mu

    n = len(df)
    print(f"Compounds compared: {n}")
    print(f"  satisfy the pharmacophore in 1XOQ : {int(df.ok_oq.sum())}")
    print(f"  satisfy it in 1XMU               : {int(df.ok_mu.sum())}")
    print(f"  satisfy it in BOTH               : "
          f"{int(df.consistent.sum())}")

    agree = int((df.ok_oq == df.ok_mu).sum())
    print(f"  the two structures agree on      : {agree} of {n} "
          f"({100*agree/n:.0f}%)")

    print(f"\n{'#':>3} {'compound':<24} {'1XOQ':>18}   {'1XMU':>18}   both")
    print(f"{'':>3} {'':<24} {'dGln  HB clamp':>18}   "
          f"{'dGln  HB clamp':>18}")
    print("-" * 78)
    for _, r in df.iterrows():
        ctrl = "*" if r.get("is_control", False) else " "

        def fmt(suf):
            d = r[f"d_gln_{suf}"]
            d = "  --- " if pd.isna(d) else f"{float(d):5.2f}"
            return (f"{d} {int(r[f'gln_hbond_{suf}']):>3} "
                    f"{int(r[f'n_clamp_{suf}']):>5}")

        mark = "YES" if r.consistent else ("one" if (r.ok_oq or r.ok_mu)
                                           else " - ")
        print(f"{int(r['rank']):>3}{ctrl}{str(r['name'])[:23]:<23} "
              f"{fmt('oq'):>18}   {fmt('mu'):>18}   {mark}")

    print("\n  * = positive control.  'YES' = pharmacophore satisfied in "
          "both structures.")

    # ---- how the controls did, which is what licenses the rest ----
    if "is_control" in df.columns:
        c = df[df.is_control.astype(bool)]
        if len(c):
            print(f"\nControls consistent across both structures: "
                  f"{int(c.consistent.sum())} of {len(c)}")
            for _, r in c.iterrows():
                print(f"  {r['name']:<20} "
                      f"{'both' if r.consistent else ('one' if (r.ok_oq or r.ok_mu) else 'neither')}")

    # ---- the verdict on the hits ----
    print("\n" + "=" * 70)
    hits = df[(~df.get("is_control", pd.Series(False, index=df.index))
               .astype(bool)) & (df["pct_of_actives_below"] >= 50)]
    if hits.empty:
        print("No library compound cleared the 50% mark in the 1XOQ "
              "ranking, so there is nothing to confirm here.")
    else:
        for _, r in hits.iterrows():
            if r.consistent:
                print(f"{r['name']}: CONSISTENT. Same binding mode in two "
                      f"independent\n  crystal conformations. The hit "
                      f"survives the robustness check.")
            elif r.ok_oq and not r.ok_mu:
                print(f"{r['name']}: CONFORMATION-DEPENDENT. It satisfies "
                      f"the pharmacophore in\n  1XOQ but not in 1XMU. "
                      f"Report that. Do not present it as a hit on the "
                      f"strength of one structure.")
            else:
                print(f"{r['name']}: does not satisfy the pharmacophore in "
                      f"1XOQ.")
    print("=" * 70)

    df.to_csv("cross_structure_comparison.csv", index=False)
    print("\nWrote cross_structure_comparison.csv")


if __name__ == "__main__":
    main()
