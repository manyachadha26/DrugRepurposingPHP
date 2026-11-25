#!/usr/bin/env python3
"""
Make sure the coordination waters actually reached the receptor pdbqt.

WHY THIS EXISTS
---------------
Receptor preparation tools routinely strip waters, because for most targets
that is the right default. Here it is not. If mk_prepare_receptor drops them,
the pdbqt we dock against has the same empty cavity we are trying to fix and
the whole re-run is wasted - silently, because nothing errors.

So: count the waters in the pdbqt. If they are missing, append them as rigid
receptor atoms typed OA.

WHAT AN OA WATER CAN AND CANNOT DO
----------------------------------
Typed OA, the water oxygen is a hydrogen-bond ACCEPTOR and it occupies space.
Occupying space is the fix: a ligand atom can no longer sit where the water
sits. What it cannot do is DONATE a hydrogen bond, because crystal waters
carry no hydrogens in the PDB and their orientation is unknown, so we are not
inventing one. State that as a limitation in the methods; do not describe
these as fully modelled waters.

Vina scores by atom type and ignores receptor partial charges, so the charge
column is written as 0.000 rather than a made-up value.

USAGE
-----
    python patch_waters_pdbqt.py <wet_receptor.pdb> <receptor.pdbqt>

Writes <receptor>_wet.pdbqt if a patch was needed, and tells you which file
to dock against either way.
"""

import os
import sys

WATER = {"HOH", "WAT", "DOD"}


def pdb_waters(path):
    """Water oxygens from the wet receptor PDB."""
    out = []
    with open(path) as fh:
        for line in fh:
            if line[:6].strip() not in ("ATOM", "HETATM"):
                continue
            if line[17:20].strip().upper() not in WATER:
                continue
            name = line[12:16].strip()
            el = line[76:78].strip().upper()
            # keep only the oxygen; a stray D or H on a water is dropped
            if el and el != "O":
                continue
            if not el and not name.startswith("O"):
                continue
            out.append({
                "name": name or "O",
                "resn": line[17:20].strip().upper(),
                "chain": line[21:22].strip() or "A",
                "resi": line[22:26].strip(),
                "x": float(line[30:38]), "y": float(line[38:46]),
                "z": float(line[46:54]),
            })
    return out


def pdbqt_water_count(path):
    n = 0
    with open(path) as fh:
        for line in fh:
            if (line[:6].strip() in ("ATOM", "HETATM")
                    and line[17:20].strip().upper() in WATER):
                n += 1
    return n


def last_serial(path):
    n = 0
    with open(path) as fh:
        for line in fh:
            if line[:6].strip() in ("ATOM", "HETATM"):
                try:
                    n = max(n, int(line[6:11]))
                except ValueError:
                    pass
    return n


def fmt(serial, w):
    """
    One PDBQT receptor atom line.

    Column layout, which has to be exact or Vina misreads the file:
      1-6 record · 7-11 serial · 13-16 name · 17 altloc · 18-20 resname
      22 chain · 23-26 resseq · 31-38 x · 39-46 y · 47-54 z
      55-60 occupancy · 61-66 b-factor · 71-76 partial charge · 78-79 type
    """
    # field widths spelled out so the off-by-one that put the atom type in
    # columns 77-78 instead of 78-79 cannot come back: Vina reads the type
    # by column, so a one-character shift silently mistypes every water
    return (
        "ATOM  "                      # 1-6   record
        f"{serial:>5}"                # 7-11  serial
        " "                           # 12
        f"{w['name']:<4}"             # 13-16 atom name
        " "                           # 17    altloc
        f"{w['resn']:>3}"             # 18-20 residue name
        " "                           # 21
        f"{w['chain']:1}"             # 22    chain
        f"{w['resi']:>4}"             # 23-26 residue number
        "    "                        # 27    icode + 28-30 blank
        f"{w['x']:8.3f}{w['y']:8.3f}{w['z']:8.3f}"   # 31-54
        f"{1.00:6.2f}{0.00:6.2f}"     # 55-66 occupancy, b-factor
        "    "                        # 67-70
        f"{0.0:6.3f}"                 # 71-76 partial charge, unused by Vina
        " "                           # 77
        "OA"                          # 78-79 AutoDock type
    )


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    pdb_path, pdbqt_path = sys.argv[1], sys.argv[2]
    for p in (pdb_path, pdbqt_path):
        if not os.path.exists(p):
            sys.exit(f"ERROR: not found: {p}")

    waters = pdb_waters(pdb_path)
    have = pdbqt_water_count(pdbqt_path)

    print(f"coordination waters in {os.path.basename(pdb_path)} : "
          f"{len(waters)}")
    print(f"waters already in {os.path.basename(pdbqt_path)}    : {have}")

    if not waters:
        print("\nThere are no waters in the wet PDB. Either the cleaning step "
              "kept none,\nor you pointed at the wrong file. Nothing to do.")
        sys.exit(1)

    if have >= len(waters):
        print("\nThe waters survived the conversion. Dock against "
              f"{pdbqt_path} as it is.")
        return

    if have:
        print(f"\nOnly {have} of {len(waters)} survived. Appending the rest "
              f"would risk duplicates,\nso this writes all {len(waters)} "
              f"fresh after removing the water lines present.")

    out = pdbqt_path.rsplit(".", 1)[0]
    out = f"{out}_wet.pdbqt" if not out.endswith("_wet") else f"{out}.pdbqt"
    if os.path.abspath(out) == os.path.abspath(pdbqt_path):
        out = pdbqt_path.rsplit(".", 1)[0] + "_withwat.pdbqt"

    serial = last_serial(pdbqt_path)
    written = 0
    with open(pdbqt_path) as src, open(out, "w") as dst:
        for line in src:
            if (line[:6].strip() in ("ATOM", "HETATM")
                    and line[17:20].strip().upper() in WATER):
                continue                      # drop any partial water lines
            if line.startswith("END"):
                continue                      # re-added at the end
            dst.write(line)
        for w in waters:
            serial += 1
            dst.write(fmt(serial, w) + "\n")
            written += 1
        dst.write("END\n")

    print(f"\nWrote {out} with {written} waters typed OA.")
    print(f"Dock against {out}, NOT {pdbqt_path}.")
    print("\nSanity check before you launch 1,237 jobs:")
    print(f"  grep -c HOH {out}        # should print {written}")
    print(f"  grep -c ' OA$' {out}     # at least {written}")


if __name__ == "__main__":
    main()
