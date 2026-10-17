#!/usr/bin/env python3
"""
Re-prepare the receptor KEEPING the waters that coordinate the metals.

WHY
---
The first preparation kept the Zn and Mg but deleted every water. In the real
structure those metals are water-coordinated: in 1XOQ, five of Mg1002's six
coordinating partners are waters, and HOH1010 bridges the two metals at
2.07 and 2.18 A. Deleting them left an empty cavity at the bimetal centre.

The consequence, measured: all 40 docked actives ended up 2.28 to 5.58 A from
a metal (median 2.86) where crystallographic roflumilast sits 4.00 to 4.99 A
away. Docking pushed ligand atoms into a hole that does not exist, and those
poses scored BETTER for it (rank vs metal distance rho = +0.326, p = 0.040).

So we keep the first coordination shell. Not every water in the file - a
receptor packed with 360 waters would block the pocket - just the ones
actually bonded to a metal.

WHAT IT KEEPS
-------------
  - chain A protein, original numbering preserved
  - the metals in chain A
  - every water within 3.0 A of a retained metal (the first shell)
  - modified amino acids that are part of the chain
  - altloc A only, where alternates exist

WHAT IT DROPS
-------------
  - all other chains
  - all other waters
  - cryoprotectants and buffer (EDO, GOL, PEG, SO4, and similar)
  - the native ligand, which goes to its own file for the grid box and for
    the redocking check

USAGE
-----
    python clean_with_waters.py raw/1XOQ.pdb [chain] [shell_A]

    chain    defaults to A
    shell_A  defaults to 3.0

OUTPUT
------
    1XOQ_wet.pdb                 receptor with its coordination waters
    1XOQ_wet_native_ligand.pdb   the crystal ligand, on its own
"""

import os
import sys
import numpy as np

METALS = {"ZN", "MG", "MN", "FE", "CA", "CO", "NI"}
WATER = {"HOH", "WAT", "DOD"}
# buffer, cryoprotectant and crystallisation junk
JUNK = {"EDO", "GOL", "PEG", "PG4", "SO4", "PO4", "CL", "NA", "K", "IOD",
        "MES", "TRS", "ACT", "ACY", "DMS", "MPD", "FMT", "BME", "NAG"}

STANDARD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}
# modified residues that sit inside the chain and must be kept, not dropped
MODIFIED_AA = {"CME", "MSE", "CSO", "CSD", "SEP", "TPO", "PTR", "KCX",
               "MLY", "HIC", "OCS", "CAS"}


def parse(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            rec = line[:6].strip()
            if rec not in ("ATOM", "HETATM"):
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
            rows.append({
                "line": line.rstrip("\n"), "rec": rec, "name": name,
                "alt": line[16:17].strip(), "resn": line[17:20].strip().upper(),
                "chain": line[21:22].strip(), "resi": line[22:27].strip(),
                "el": el, "xyz": xyz,
            })
    return rows


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    path = sys.argv[1]
    chain = sys.argv[2] if len(sys.argv) > 2 else "A"
    shell = float(sys.argv[3]) if len(sys.argv) > 3 else 3.0

    if not os.path.exists(path):
        sys.exit(f"ERROR: not found: {path}")

    atoms = parse(path)
    if not atoms:
        sys.exit(f"ERROR: no atoms parsed from {path}. Is it really a PDB?")

    stem = os.path.basename(path).rsplit(".", 1)[0]
    out_rec = f"{stem}_wet.pdb"
    out_lig = f"{stem}_wet_native_ligand.pdb"

    n_wat_total = sum(1 for a in atoms if a["resn"] in WATER)
    print(f"{os.path.basename(path)}: {len(atoms)} atoms, "
          f"{n_wat_total} waters in file")

    # altloc: blank or A only
    atoms = [a for a in atoms if a["alt"] in ("", "A")]

    this = [a for a in atoms if a["chain"] == chain]
    if not this:
        chains = sorted({a["chain"] for a in atoms})
        sys.exit(f"ERROR: no chain {chain}. Chains present: {chains}")

    # ---- the metals of this chain ----
    metals = [a for a in this
              if a["el"] in METALS or a["resn"] in METALS]
    if not metals:
        print("  WARNING: no metals in this chain. Nothing to keep waters "
              "around, and the whole point of this script is the metals.")
    else:
        print("  metals kept     : " + ", ".join(
            f"{m['resn']}{m['resi']}" for m in metals))

    m_xyz = np.array([m["xyz"] for m in metals]) if metals else None

    # ---- the native ligand: biggest HETATM residue that is not water,
    #      metal or junk ----
    groups = {}
    for a in this:
        if (a["rec"] == "HETATM" and a["resn"] not in WATER
                and a["resn"] not in METALS and a["resn"] not in JUNK
                and a["resn"] not in STANDARD_AA
                and a["resn"] not in MODIFIED_AA):
            groups.setdefault((a["resn"], a["resi"]), []).append(a)
    native = max(groups.values(), key=len) if groups else []
    native_key = (native[0]["resn"], native[0]["resi"]) if native else None
    if native:
        print(f"  native ligand   : {native_key[0]}{native_key[1]} "
              f"({len(native)} atoms) -> {out_lig}")
    else:
        print("  native ligand   : none found")

    # ---- decide what to keep ----
    keep, kept_waters = [], set()
    for a in this:
        resn = a["resn"]

        if resn in WATER:
            if m_xyz is None:
                continue
            d = np.linalg.norm(m_xyz - np.array(a["xyz"]), axis=1).min()
            if d <= shell:
                keep.append(a)
                kept_waters.add(a["resi"])
            continue

        if resn in METALS or a["el"] in METALS:
            keep.append(a)
            continue

        if (resn, a["resi"]) == native_key:
            continue                      # goes to its own file

        if resn in JUNK:
            continue

        if resn in STANDARD_AA or resn in MODIFIED_AA or a["rec"] == "ATOM":
            keep.append(a)
            continue
        # anything else unrecognised is dropped, and named so it is not silent
        print(f"  dropped residue : {resn}{a['resi']} (unrecognised)")

    print(f"  waters kept     : {len(kept_waters)} within {shell} A of a "
          f"metal, out of {n_wat_total} in the file")
    if metals and not kept_waters:
        print("  WARNING: no waters were within range of a metal. Check the")
        print("  input really is the raw PDB and not an already-cleaned file.")

    def blank_altloc(line):
        """
        Clear column 17.

        We have already chosen altloc A above, so leaving the letter in place
        makes the file look like it still offers a choice. Meeko then refuses
        the structure with "Residues with alternate location" and nothing is
        written. Blanking it states the choice we made.
        """
        return line[:16] + " " + line[17:] if len(line) > 17 else line

    n_alt = sum(1 for a in keep if a["alt"])
    with open(out_rec, "w") as fh:
        for a in keep:
            fh.write(blank_altloc(a["line"]) + "\n")
        fh.write("END\n")
    if native:
        with open(out_lig, "w") as fh:
            for a in native:
                fh.write(blank_altloc(a["line"]) + "\n")
            fh.write("END\n")
    if n_alt:
        print(f"  altloc          : kept variant A for {n_alt} atoms, "
              f"column 17 blanked")

    print(f"  wrote {out_rec}  ({len(keep)} atoms)")

    # ---- the number that matters: how close can a ligand now get? ----
    if metals and kept_waters:
        wat = [a for a in keep if a["resn"] in WATER]
        w_xyz = np.array([a["xyz"] for a in wat])
        dmin = np.linalg.norm(
            w_xyz[:, None, :] - m_xyz[None, :, :], axis=2).min(axis=1)
        print(f"\n  kept waters sit {dmin.min():.2f} to {dmin.max():.2f} A "
              f"from a metal")
        print("  Those positions are now occupied, so docking can no longer")
        print("  put a ligand atom there. That is the whole fix.")

    print(f"\nNext: convert to pdbqt, then CHECK the waters survived:")
    print(f"  mk_prepare_receptor.py -i {out_rec} -o {stem}_wet -p")
    print(f"  python patch_waters_pdbqt.py {out_rec} {stem}_wet.pdbqt")


if __name__ == "__main__":
    main()
