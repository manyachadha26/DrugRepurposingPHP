#!/usr/bin/env python3
"""
Rank the antidiabetic library by interaction pattern, not by docking score.

WHY IT WORKS THIS WAY
---------------------
The docking score cannot rank PDE4D binders: AUC 0.604 for Vina and 0.586
for Vinardo on 40 actives against 1,197 property-matched decoys, both barely
above a coin flip. Classifying on WHERE the pose sits instead reached 0.727
under leave-one-scaffold-out cross-validation, with a label-permutation
control at 0.465.

So this trains that classifier on the whole benchmark and applies it to the
library. Training on all 1,237 is correct here: the library compounds are
different molecules, so there is no leakage to guard against. The 0.727 from
cross-validation is what tells you how much to trust the output - this script
does not re-estimate it, and nothing it prints is evidence of performance.

WHAT THE OUTPUT IS FOR
----------------------
A probability alone is an opaque number no reviewer will accept. So each
compound is reported with the interactions behind it: the hydrogen bond to
the invariant glutamine, clamp contacts, distance to the bimetal centre,
burial. The ranking and its reasons arrive together.

WHAT IT IS NOT
--------------
Not a potency prediction. Not an affinity estimate. There is no claim here
that a high-probability compound inhibits PDE4D, only that its pose
reproduces the interaction pattern shared by known PDE4 inhibitors. Write it
that way in the paper.

USAGE
-----
    python screen_library.py benchmark_fp.csv library_fp.csv [controls.txt]

    benchmark_fp.csv   fingerprints.csv from the 1,237-compound benchmark
    library_fp.csv     fingerprints of the docked library
    controls.txt       optional, one control name per line (theophylline,
                       rolipram), which are flagged in the output

OUTPUT
------
    screen_results.csv   ranked, with every interaction column kept
"""

import os
import sys
import numpy as np
import pandas as pd

from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)

    bench_path, lib_path = sys.argv[1], sys.argv[2]
    ctrl_path = sys.argv[3] if len(sys.argv) > 3 else None
    for p in (bench_path, lib_path):
        if not os.path.exists(p):
            sys.exit(f"ERROR: not found: {p}")

    bench = pd.read_csv(bench_path)
    lib = pd.read_csv(lib_path)

    controls = set()
    if ctrl_path and os.path.exists(ctrl_path):
        with open(ctrl_path) as fh:
            controls = {l.strip().lower() for l in fh if l.strip()}

    if "label" not in bench.columns:
        sys.exit("ERROR: the benchmark file needs a label column.")

    # ---- features must be the SAME columns in the same order ----
    # The library was fingerprinted against the same receptor, so the
    # per-residue columns should match. Any that do not are dropped rather
    # than filled, because a silently zero-filled contact column would look
    # like "no contact" and quietly bias every library compound.
    drop = {"name", "label", "scaffold", "hac"}
    bcols = [c for c in bench.columns if c not in drop]
    lcols = [c for c in lib.columns if c not in drop]
    feat = [c for c in bcols if c in lcols]
    missing = [c for c in bcols if c not in lcols]
    if missing:
        print(f"WARNING: {len(missing)} feature columns are in the benchmark "
              f"but not the library and are dropped:")
        print("  " + ", ".join(missing[:12])
              + (" ..." if len(missing) > 12 else ""))
        print("  If that list is long, the two fingerprint runs used "
              "different receptors. Stop and check.")

    def numeric(df, cols):
        X = df[cols].apply(pd.to_numeric, errors="coerce")
        return X.fillna(X.median(numeric_only=True)).fillna(0.0)

    Xb = numeric(bench, feat)
    keep = [c for c in Xb.columns if Xb[c].std() > 0]
    Xb = Xb[keep].values.astype(float)
    Xl = numeric(lib, keep).values.astype(float)
    y = bench["label"].values.astype(int)

    print(f"Training on   : {len(bench)} compounds "
          f"({int(y.sum())} actives, {len(bench)-int(y.sum())} decoys)")
    print(f"Scoring       : {len(lib)} library compounds")
    print(f"Features used : {len(keep)}")
    print(f"Trust level   : AUC 0.727 from cross-validation on the "
          f"benchmark.\n                Treat this as a triage ordering, "
          f"not a potency ranking.\n")

    clf = make_pipeline(
        StandardScaler(),
        LogisticRegression(max_iter=4000, class_weight="balanced"),
    ).fit(Xb, y)

    lib = lib.copy()
    lib["probability"] = clf.predict_proba(Xl)[:, 1]

    # where each library compound sits relative to the known actives
    act = bench[bench.label == 1]
    pa = clf.predict_proba(numeric(act, keep).values.astype(float))[:, 1]
    lib["pct_of_actives_below"] = [
        float((pa < p).mean() * 100) for p in lib["probability"]]

    lib["is_control"] = [n.lower() in controls for n in lib["name"]]
    lib = lib.sort_values("probability", ascending=False).reset_index(drop=True)
    lib.insert(0, "rank", np.arange(1, len(lib) + 1))
    lib.to_csv("screen_results.csv", index=False)

    # ---- the readable table ----
    cols = [c for c in ("d_gln", "gln_hbond", "n_clamp", "d_metal",
                        "f_buried", "n_res") if c in lib.columns]
    print(f"{'#':>3}  {'compound':<24} {'prob':>5}  {'%act':>5}  "
          + "  ".join(f"{c:>9}" for c in cols))
    print("-" * (42 + 11 * len(cols)))
    for _, r in lib.iterrows():
        tag = " *" if r["is_control"] else "  "
        vals = []
        for c in cols:
            v = r[c]
            vals.append(f"{'':>9}" if pd.isna(v) else
                        (f"{int(v):>9}" if c in ("gln_hbond", "n_clamp",
                                                 "n_res")
                         else f"{float(v):>9.2f}"))
        print(f"{int(r['rank']):>3}{tag}{str(r['name'])[:24]:<24} "
              f"{r['probability']:>5.3f}  {r['pct_of_actives_below']:>4.0f}%  "
              + "  ".join(vals))
    if controls:
        print("\n  * = positive control")

    # ---- how to read it, including the outcome nobody wants ----
    print("\n" + "=" * 62)
    print("HOW TO READ THIS")
    print("=" * 62)
    print("  'prob' is how much the pose resembles a known PDE4 inhibitor's")
    print("  interaction pattern. 'pct' is the share of the 40 known actives")
    print("  this compound scores above.")
    print()
    ctrl = lib[lib.is_control]
    if len(ctrl):
        print("  CHECK THE CONTROLS FIRST. Rolipram and theophylline are")
        print("  known PDE4 inhibitors and should rank high. Where they land:")
        for _, r in ctrl.iterrows():
            print(f"    {r['name']}: rank {int(r['rank'])} of {len(lib)}, "
                  f"prob {r['probability']:.3f}, "
                  f"above {r['pct_of_actives_below']:.0f}% of known actives")
        if (ctrl["pct_of_actives_below"] < 50).any():
            print("    At least one control ranks below half the known")
            print("    actives. That is a problem with the protocol, not")
            print("    with the control. Do not report the library ranking")
            print("    until it is explained.")
    else:
        print("  NO CONTROLS FOUND. Add theophylline and rolipram and re-run.")
        print("  Without them there is no way to tell a working screen from")
        print("  a broken one.")

    top = lib[lib.pct_of_actives_below >= 50]
    nc = top[~top.is_control]
    print()
    if len(nc) == 0:
        print("  NO library compound scores above half the known actives.")
        print("  That is a legitimate result: on this evidence these")
        print("  antidiabetics do not reproduce the PDE4 interaction")
        print("  pattern. Report it as a negative screen. Do not lower the")
        print("  threshold to manufacture hits.")
    else:
        print(f"  {len(nc)} library compounds score above at least half the")
        print(f"  known actives. These are candidates for experimental")
        print(f"  testing - not predicted inhibitors. Inspect each pose by")
        print(f"  eye before putting any of them in the paper.")


if __name__ == "__main__":
    main()
