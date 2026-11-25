#!/usr/bin/env python3
"""
Rebuild incomplete side chains, WITHOUT losing the metals or the waters.

WHY
---
Meeko refused 7CBQ_wet.pdb: eight residues have side chains the
crystallographer could not resolve. Seven sit 21-31 A from the ligand and
could be deleted harmlessly. Arg342 sits 9.0 A away, inside the grid box -
deleting it would carve exactly the kind of void we just spent a day
removing from this project.

So we rebuild the missing atoms instead.

WHAT IT DOES NOT DO
-------------------
  - does not build missing RESIDUES (whole unresolved loops). Those are
    absent for a reason and inventing a loop conformation is worse than
    leaving the gap
  - does not call removeHeterogens, which would strip the Zn, the Mg and
    the coordination waters we deliberately kept
  - does not replace non-standard residues. CME is a real modified cysteine
    in this structure, not an error to be corrected away

It counts metals and waters before and after and refuses to report success
if either changed, because losing them silently is the exact failure this
whole re-preparation was fixing.

USAGE
-----
    python fix_sidechains.py 7CBQ_wet.pdb

OUTPUT
------
    7CBQ_wetfix.pdb
"""

import os
import sys

METALS = {"ZN", "MG", "MN", "FE", "CA", "CO", "NI"}
WATER = {"HOH", "WAT", "DOD"}


def census(path):
    """Count metals, waters and total heavy atoms in a PDB."""
    metals, waters, atoms = [], 0, 0
    with open(path) as fh:
        for line in fh:
            if line[:6].strip() not in ("ATOM", "HETATM"):
                continue
            atoms += 1
            resn = line[17:20].strip().upper()
            if resn in WATER:
                waters += 1
            elif resn in METALS:
                metals.append(f"{resn}{line[22:26].strip()}")
    return metals, waters, atoms


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    src = sys.argv[1]
    if not os.path.exists(src):
        sys.exit(f"ERROR: not found: {src}")

    try:
        from pdbfixer import PDBFixer
        from openmm.app import PDBFile
    except ImportError as exc:
        sys.exit(f"ERROR: pdbfixer/openmm not available ({exc}).\n"
                 f"  conda install -c conda-forge pdbfixer openmm")

    stem = os.path.basename(src).rsplit(".", 1)[0]
    out = f"{stem}fix.pdb" if stem.endswith("_wet") else f"{stem}_fixed.pdb"

    before_metals, before_waters, before_atoms = census(src)
    print(f"{os.path.basename(src)}")
    print(f"  atoms  : {before_atoms}")
    print(f"  metals : {', '.join(before_metals) if before_metals else 'NONE'}")
    print(f"  waters : {before_waters}")

    fixer = PDBFixer(filename=src)

    # Missing RESIDUES: find them so the library knows about the gaps, then
    # clear the list so none are built. An invented loop is a worse artifact
    # than an honest gap.
    fixer.findMissingResidues()
    if fixer.missingResidues:
        print(f"\n  {len(fixer.missingResidues)} gap(s) in the chain - "
              f"left as gaps, not rebuilt")
        fixer.missingResidues = {}

    fixer.findMissingAtoms()
    n_res = len(fixer.missingAtoms)
    if n_res:
        print(f"\n  rebuilding side-chain atoms in {n_res} residue(s):")
        for res, atoms in list(fixer.missingAtoms.items())[:12]:
            names = ",".join(a.name for a in atoms)
            print(f"    {res.name}{res.id}: {names}")
        if n_res > 12:
            print(f"    ... and {n_res - 12} more")
    else:
        print("\n  no missing atoms found - nothing to rebuild")

    if fixer.missingTerminals:
        print(f"  {len(fixer.missingTerminals)} terminal(s) completed")

    fixer.addMissingAtoms()

    with open(out, "w") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh,
                          keepIds=True)

    after_metals, after_waters, after_atoms = census(out)
    print(f"\nWrote {out}")
    print(f"  atoms  : {before_atoms} -> {after_atoms} "
          f"(+{after_atoms - before_atoms})")
    print(f"  metals : {', '.join(after_metals) if after_metals else 'NONE'}")
    print(f"  waters : {after_waters}")

    ok = True
    if sorted(after_metals) != sorted(before_metals):
        print("\n  *** METALS CHANGED. Do not dock against this file.")
        ok = False
    if after_waters != before_waters:
        print(f"\n  *** WATERS CHANGED: {before_waters} -> {after_waters}.")
        print("  *** The coordination shell is what we were protecting.")
        print("  *** Do not dock against this file.")
        ok = False
    if after_atoms < before_atoms:
        print("\n  *** ATOMS WERE LOST, not just added. Check before using.")
        ok = False

    if ok:
        print("\n  Metals and waters intact. Next:")
        print(f"    mk_prepare_receptor.py -i {out} -o {stem} -p")
        print(f"    python patch_waters_pdbqt.py {out} {stem}.pdbqt")
    else:
        sys.exit(1)


if __name__ == "__main__":
    main()
