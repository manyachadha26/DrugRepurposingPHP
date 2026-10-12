"""
Gate 1: how far is the docked pose from the crystal pose?

Usage:
    python calc_rmsd.py redock_1XOQ.pdbqt 1XOQ_native_ligand.pdbqt

The first file is Vina's output, the second is the ligand taken out of the
crystal. Under 2.0 A means your setup reproduces reality and you can trust it.

Symmetry is handled. Roflumilast's dichloropyridine ring can be flipped 180
degrees and look identical, so a naive atom-by-atom comparison reports a large
error for two poses that are chemically the same. RDKit's CalcRMS tries the
equivalent atom matchings and takes the best one.

The molecules are NOT superimposed first. That matters: superimposing would
report a perfect score for a pose docked into entirely the wrong part of the
protein, because the shapes would still match once moved on top of each other.
We want the distance where the pose actually sits.
"""

import os
import subprocess
import sys
import tempfile


def to_sdf(path, tmpdir, tag):
    """Convert pdbqt/pdb to sdf with Open Babel so RDKit can read bond orders."""
    out = os.path.join(tmpdir, f"{tag}.sdf")
    r = subprocess.run(["obabel", path, "-O", out], capture_output=True, text=True)
    if not os.path.exists(out) or os.path.getsize(out) == 0:
        sys.exit(f"Open Babel could not convert {path}\n{r.stderr[:400]}")
    return out


def naive_rmsd(a, b):
    """Fallback: paired heavy-atom RMSD, atom order as written."""
    import math

    def heavy(p):
        out = []
        for line in open(p):
            if line.startswith(("ATOM", "HETATM")):
                el = line[76:78].strip() or line[12:16].strip()[0]
                if el.upper().startswith("H"):
                    continue
                out.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
        return out

    x, y = heavy(a), heavy(b)
    if len(x) != len(y):
        return None
    s = sum((p[0]-q[0])**2 + (p[1]-q[1])**2 + (p[2]-q[2])**2 for p, q in zip(x, y))
    return math.sqrt(s / len(x))


def main():
    if len(sys.argv) != 3:
        sys.exit("Usage: python calc_rmsd.py <docked.pdbqt> <native_ligand.pdbqt>")
    docked, native = sys.argv[1], sys.argv[2]
    for f in (docked, native):
        if not os.path.exists(f):
            sys.exit(f"Not found: {f}")

    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem import AllChem, rdMolAlign
        RDLogger.DisableLog("rdApp.*")
    except ImportError:
        v = naive_rmsd(docked, native)
        print(f"RDKit not installed - naive RMSD (symmetry NOT handled): {v}")
        print("Install it for the correct number:  conda install -c conda-forge rdkit")
        return

    with tempfile.TemporaryDirectory() as td:
        ref = Chem.SDMolSupplier(to_sdf(native, td, "ref"), removeHs=True)[0]
        poses = [m for m in Chem.SDMolSupplier(to_sdf(docked, td, "poses"), removeHs=True) if m]
        if ref is None or not poses:
            sys.exit("Could not read the molecules. Check both files are the right ligand.")

        print(f"Reference (crystal): {ref.GetNumAtoms()} heavy atoms")
        print(f"Docked poses found:  {len(poses)}\n")

        best = None          # lowest RMSD anywhere in the list
        top = None           # RMSD of the pose Vina ranked FIRST
        for i, p in enumerate(poses, start=1):
            try:
                # CalcRMS, not GetBestRMS: GetBestRMS superimposes the two
                # molecules first, which would report 0 for a pose sitting in
                # completely the wrong place. We need the distance as docked.
                v = rdMolAlign.CalcRMS(Chem.Mol(p), Chem.Mol(ref))
            except Exception:
                v = None
            if v is not None:
                if i == 1:
                    top = v
                if best is None or v < best[1]:
                    best = (i, v)
            print(f"   pose {i}: RMSD = {v:.3f} A" if v is not None
                  else f"   pose {i}: could not compare")

        if best is None:
            sys.exit("\nNo pose could be compared. Are both files the same molecule?")

        bi, bv = best
        print("\n" + "=" * 62)
        print("TWO DIFFERENT QUESTIONS - both matter, and they can disagree:\n")
        print(f"  1. Did the TOP-RANKED pose get it right?   {top:.3f} A"
              if top is not None else "  1. Top-ranked pose: could not compare")
        print(f"     (this is the one that counts for screening)")
        print(f"  2. Was the right pose found ANYWHERE?      {bv:.3f} A  (pose {bi})")
        print(f"     (this only shows the search explored the pocket properly)")
        print()

        if top is not None and top < 2.0:
            print("PASS. The best-scoring pose is the correct one, so the scoring")
            print("function can be trusted to rank unknown compounds.")
        elif bv < 2.0:
            print("PARTIAL. The search FOUND the correct pose but did not rank it")
            print(f"first - it came {bi}th. Sampling works; scoring does not.")
            print("For a screen you only ever look at the top pose, so a ranking")
            print("that puts the right answer third cannot be relied on here.")
            print("Try: raise exhaustiveness to 64, confirm the receptor is the")
            print("PDBFixer-repaired one, and check the pocket is complete.")
            print("If it still fails, use this structure for comparing poses")
            print("only - not for scoring - and say so in the methods.")
        else:
            print("FAIL. The correct pose was never found. Do not screen with this.")
            print("Check the metals are still present, the box is centred on the")
            print("ligand, and protonation was handled.")
        print("=" * 62)

        print("\nReport Vina's own score for pose 1 as well - it is in the docked file:")
        print(f"   grep 'VINA RESULT' {docked} | head -1")


if __name__ == "__main__":
    main()
