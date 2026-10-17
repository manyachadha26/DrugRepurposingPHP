#!/usr/bin/env python3
"""
Describe the bimetal centre of each PDE4D structure, and test whether the
actives the docking score misranks are the metal-dependent ones.

WHY THIS IS THE IMPORTANT SCRIPT
--------------------------------
We believe the scoring functions fail because the PDE4D catalytic site
coordinates Zn2+ and Mg2+ through a bridging water network, and neither Vina
nor Vinardo models metal coordination explicitly. Right now that is a story.
This script turns it into a measurement, which is the difference between a
discussion paragraph and a result.

Two parts:

PART 1, per structure: the coordination shell of each metal - which residue
atoms and waters complete it, at what distance - and how far the native
ligand sits from each metal. Run this on the RAW downloaded PDB, not the
cleaned one, because cleaning removed the waters and the bridging water is
half the argument.

PART 2, the test: split the docked actives by whether their pose reaches the
metal shell, then compare docking rank between the two groups. If the
metal-reaching actives are the ones ranked worst, the hypothesis has
evidence behind it. If rank is the same either way, it does not, and we say
so - a negative result here is still a result, and it would send us to look
for a different explanation rather than publishing a guess.

USAGE
-----
  Part 1:
    python metal_site.py shell 7CBQ.pdb 1XOQ.pdb 1XMU.pdb
    (raw PDBs, waters intact)

  Part 2:
    python metal_site.py test fingerprints.csv rescore_results.csv
"""

import os
import sys
import numpy as np

METALS = {"ZN", "MG", "MN", "FE"}
WATER = {"HOH", "WAT", "DOD"}
SHELL = 2.8          # first coordination shell, A
SHELL_WIDE = 3.2     # a looser pass, since crystal geometry varies


def parse_pdb(path):
    out = []
    with open(path) as fh:
        for line in fh:
            if line[:6].strip() not in ("ATOM", "HETATM"):
                continue
            try:
                xyz = (float(line[30:38]), float(line[38:46]),
                       float(line[46:54]))
            except ValueError:
                continue
            name = line[12:16].strip()
            el = line[76:78].strip().upper()
            if not el:
                letters = "".join(c for c in name if c.isalpha()).upper()
                el = letters[:2] if letters[:2] in METALS else letters[:1]
            out.append({
                "name": name, "resn": line[17:20].strip().upper(),
                "chain": line[21:22].strip(), "resi": line[22:27].strip(),
                "alt": line[16:17].strip(), "el": el, "xyz": xyz,
                "hetatm": line[:6].strip() == "HETATM",
            })
    return out


def shells(path):
    atoms = parse_pdb(path)
    if not atoms:
        print(f"  ERROR: no atoms parsed from {path}")
        return
    xyz = np.array([a["xyz"] for a in atoms])

    metals = [i for i, a in enumerate(atoms)
              if a["el"] in METALS or a["resn"] in METALS]
    if not metals:
        print("  no Zn/Mg/Mn/Fe in this file")
        return

    # the native ligand: the largest non-water, non-metal HETATM residue
    groups = {}
    for i, a in enumerate(atoms):
        if a["hetatm"] and a["resn"] not in WATER and a["resn"] not in METALS:
            groups.setdefault((a["chain"], a["resn"], a["resi"]), []).append(i)
    native = max(groups.values(), key=len) if groups else []
    if native:
        a0 = atoms[native[0]]
        print(f"  native ligand    : {a0['resn']} {a0['chain']}{a0['resi']}"
              f"  ({len(native)} atoms)")
    else:
        print("  native ligand    : none found (apo structure?)")

    n_wat = sum(1 for a in atoms if a["resn"] in WATER)
    print(f"  waters in file   : {n_wat}"
          + ("   <-- none; use the RAW pdb, the bridging water matters"
             if n_wat == 0 else ""))

    for mi in metals:
        m = atoms[mi]
        d = np.linalg.norm(xyz - xyz[mi], axis=1)
        print(f"\n  {m['resn']}{m['resi']} ({m['el']}) coordination:")
        partners = []
        for j in np.argsort(d):
            if j == mi or d[j] > SHELL_WIDE:
                continue
            a = atoms[j]
            if a["el"] in ("C",) or a["el"].startswith("H"):
                continue          # carbon does not coordinate a metal here
            kind = "water" if a["resn"] in WATER else (
                "metal" if (a["el"] in METALS or a["resn"] in METALS)
                else "protein")
            partners.append((d[j], a, kind))
        if not partners:
            print("    nothing within 3.2 A - check the file")
        for dd, a, kind in partners:
            flag = "" if dd <= SHELL else "  (outer, 2.8-3.2 A)"
            print(f"    {a['resn']}{a['resi']:<5} {a['name']:<4} "
                  f"{dd:5.2f} A   {kind}{flag}")

        if native:
            dn = np.linalg.norm(xyz[native] - xyz[mi], axis=1)
            print(f"    closest native-ligand atom: {dn.min():.2f} A "
                  f"({atoms[native[int(np.argmin(dn))]]['name']})")

    # is there a water bridging the two metals?
    if len(metals) >= 2:
        wat = [i for i, a in enumerate(atoms) if a["resn"] in WATER]
        if wat:
            best = None
            for w in wat:
                ds = [np.linalg.norm(xyz[w] - xyz[mi]) for mi in metals]
                if max(ds) <= SHELL_WIDE + 0.4:
                    tot = sum(ds)
                    if best is None or tot < best[0]:
                        best = (tot, w, ds)
            if best:
                a = atoms[best[1]]
                print(f"\n  BRIDGING WATER   : {a['resn']}{a['resi']} at "
                      + ", ".join(f"{d:.2f} A" for d in best[2])
                      + " from the metals")
                print("  This is the species neither scoring function "
                      "represents.")
            else:
                print("\n  no single water bridges both metals within 3.6 A")


def test(fp_csv, score_csv):
    import pandas as pd
    fp = pd.read_csv(fp_csv)
    sc = pd.read_csv(score_csv)
    df = fp.merge(sc, on="name", suffixes=("", "_s"))
    lab = "label" if "label" in df.columns else "label_s"
    acts = df[df[lab] == 1].copy()
    if acts.empty:
        sys.exit("ERROR: no actives after the merge. Check the name columns "
                 "match between the two files.")
    if "metal_shell" not in acts.columns:
        sys.exit("ERROR: no metal_shell column. Re-run fingerprint.py on a "
                 "receptor that still has its metals.")

    from scipy.stats import spearmanr, mannwhitneyu

    print(f"Actives merged: {len(acts)}")
    if "d_metal" not in acts.columns:
        sys.exit("ERROR: no d_metal column in the fingerprint file.")

    dm = acts["d_metal"].dropna()
    if dm.empty:
        sys.exit("ERROR: d_metal is empty. The receptor had no metals.")
    print(f"\nHow far the actives sit from a metal:")
    print(f"  min {dm.min():.2f}   25% {dm.quantile(.25):.2f}   "
          f"median {dm.median():.2f}   75% {dm.quantile(.75):.2f}   "
          f"max {dm.max():.2f} A")
    if dm.max() - dm.min() < 1.5:
        print("  They are all at essentially the same distance, so there is")
        print("  no contrast to test. Report the distribution and stop here.")
        return
    print("  The binary metal_shell flag is NOT used: a 6 A cut-off put all")
    print("  40 actives on one side of it, which is no test at all. The")
    print("  continuous distance below carries the information instead.")

    for col in ("vina_score", "vinardo_score"):
        if col not in df.columns:
            continue
        d = df.dropna(subset=[col]).copy()
        d["rank"] = d[col].rank(method="min")
        a = d[(d[lab] == 1) & d["d_metal"].notna()].copy()
        if len(a) < 8:
            print(f"\n=== {col} ===\n  only {len(a)} actives; too few.")
            continue

        print(f"\n=== {col} ===")
        # The hypothesis: actives that bind THROUGH the metal centre are the
        # ones the score gets wrong. Closer to the metal (small d_metal)
        # should therefore mean a WORSE rank (large rank number), i.e. a
        # NEGATIVE correlation between d_metal and rank.
        rho, p = spearmanr(a["d_metal"], a["rank"])
        print(f"  distance to metal vs rank: rho = {rho:+.3f}, p = {p:.4f}")
        print(f"  (negative rho = the metal-proximal actives rank worse, "
              f"which is what the hypothesis predicts)")

        med = a["d_metal"].median()
        near = a[a["d_metal"] <= med]
        far = a[a["d_metal"] > med]
        print(f"  closer half (<= {med:.2f} A): {len(near):>2} actives, "
              f"median rank {near['rank'].median():.0f}")
        print(f"  farther half (> {med:.2f} A): {len(far):>2} actives, "
              f"median rank {far['rank'].median():.0f}")
        if len(near) >= 3 and len(far) >= 3:
            u, pm = mannwhitneyu(near["rank"], far["rank"],
                                 alternative="greater")
            print(f"  closer half ranked worse? p = {pm:.4f}")
            if p < 0.05 and rho < 0:
                print("  -> SUPPORTED. The actives sitting nearest the metal")
                print("     centre are the ones the score misranks.")
            elif p < 0.05 and rho > 0:
                print("  -> OPPOSITE of the hypothesis. Metal-proximal")
                print("     actives rank BETTER. The explanation is wrong as")
                print("     stated; do not publish it. Report this instead.")
            else:
                print("  -> NOT SUPPORTED. Metal distance does not explain")
                print("     which actives are misranked. Report that plainly")
                print("     and look for another explanation. Do not assert")
                print("     the metal story without this evidence.")


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    mode = sys.argv[1]
    if mode == "shell":
        for p in sys.argv[2:]:
            print("=" * 62)
            print(os.path.basename(p))
            print("=" * 62)
            if not os.path.exists(p):
                print("  file not found")
                continue
            shells(p)
            print()
    elif mode == "test":
        if len(sys.argv) < 4:
            sys.exit("usage: python metal_site.py test fingerprints.csv "
                     "rescore_results.csv")
        test(sys.argv[2], sys.argv[3])
    else:
        print(__doc__)
        sys.exit(1)


if __name__ == "__main__":
    main()
