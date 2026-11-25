"""
Step 4: prepare the PDE structures for docking.

For each PDB file this writes two things:
    <ID>_clean.pdb          the receptor: chain A protein + Zn + Mg
    <ID>_native_ligand.pdb  the drug that came in the crystal

What it removes, and why:
    - chains other than A      several copies of the protein; docking needs one
    - water (HOH)              standard for docking
    - ethanediol (EDO)         antifreeze from crystallisation, not biology
    - alternate conformations  side chains modelled in two places at once;
                               conformation A is kept

What it deliberately KEEPS:
    - ZN and MG                both sit in the active site and are part of how
                               the enzyme works. Delete them and the pocket
                               collapses. This is why you must never strip all
                               HETATM records in one go.
    - CME and other modified   these look like HETATM but are amino acids inside
      amino acids              the chain. Deleting one leaves a hole in the protein.

Residue numbering is left exactly as deposited, so the glutamine numbers from
Step 3 (Gln369 / Gln443 / Gln466) stay valid.

Run:  python clean_structures.py 1XOQ.pdb 7CBQ.pdb 1XMU.pdb
"""

import os
import sys
from collections import Counter

KEEP_CHAIN = "A"
METALS = {"ZN", "MG", "MN", "CA", "FE", "CO", "NI", "CU"}
DROP = {"HOH", "WAT", "DOD", "EDO", "GOL", "PEG", "SO4", "PO4", "ACT", "MPD",
        "TRS", "DMS", "IOD", "FMT", "ACY", "NO3", "EPE", "MES", "CL", "BR"}
# HETATM records that are really amino acids sitting in the polypeptide chain
MODIFIED_AA = {"MSE", "CME", "CSO", "CSD", "OCS", "SEP", "TPO", "PTR", "KCX",
               "LLP", "MLY", "M3L", "HYP", "PCA", "CAS", "CSS", "SMC"}


def clean(path):
    pdbid = os.path.splitext(os.path.basename(path))[0].upper()
    protein, ligand = [], []
    dropped, kept_het = Counter(), Counter()
    alt_dropped = 0
    lig_resn = None

    for line in open(path, errors="replace"):
        rec = line[:6]
        if rec not in ("ATOM  ", "HETATM"):
            continue
        chain = line[21]
        resn = line[17:20].strip()
        alt = line[16]

        if chain != KEEP_CHAIN:
            continue
        if alt not in (" ", "A"):
            alt_dropped += 1
            continue
        line = line[:16] + " " + line[17:]          # blank the altloc flag

        if rec == "ATOM  " or resn in MODIFIED_AA:
            protein.append(line)
            if resn in MODIFIED_AA:
                kept_het[resn] += 1
            continue
        if resn in METALS:
            protein.append(line)
            kept_het[resn] += 1
            continue
        if resn in DROP:
            dropped[resn] += 1
            continue
        # anything left that is big enough to be a drug is the native ligand
        ligand.append(line)
        lig_resn = resn

    if not protein:
        print(f"!! {pdbid}: nothing kept - is chain {KEEP_CHAIN} present?")
        return

    rec_path = f"{pdbid}_clean.pdb"
    with open(rec_path, "w") as fh:
        fh.write(f"REMARK   Receptor prepared from {os.path.basename(path)} "
                 f"by clean_structures.py\n")
        fh.write(f"REMARK   Chain {KEEP_CHAIN} only; waters and crystallisation "
                 f"additives removed; metals retained;\n")
        fh.write("REMARK   alternate conformations reduced to A; original "
                 "residue numbering preserved.\n")
        fh.writelines(protein)
        fh.write("TER\nEND\n")

    msg = [f"### {pdbid}",
           f"   receptor      -> {rec_path}  ({len(protein)} atoms)"]
    if ligand:
        lig_path = f"{pdbid}_native_ligand.pdb"
        with open(lig_path, "w") as fh:
            fh.write(f"REMARK   Native ligand {lig_resn} extracted from "
                     f"{os.path.basename(path)}\n")
            fh.writelines(ligand)
            fh.write("END\n")
        xs = [float(l[30:38]) for l in ligand]
        ys = [float(l[38:46]) for l in ligand]
        zs = [float(l[46:54]) for l in ligand]
        n = len(ligand)
        msg.append(f"   native ligand -> {lig_path}  ({lig_resn}, {n} atoms)")
        msg.append(f"   grid centre      center_x = {sum(xs)/n:.3f}   "
                   f"center_y = {sum(ys)/n:.3f}   center_z = {sum(zs)/n:.3f}")
    else:
        msg.append("   !! no native ligand found")

    if kept_het:
        msg.append("   KEPT  " + ", ".join(f"{k} x{v}" for k, v in sorted(kept_het.items()))
                   + "   <- metals and modified residues")
    if dropped:
        msg.append("   removed " + ", ".join(f"{k} x{v}" for k, v in sorted(dropped.items())))
    if alt_dropped:
        msg.append(f"   removed {alt_dropped} atoms in alternate conformations")

    has_zn = any(k == "ZN" for k in kept_het)
    has_mg = any(k == "MG" for k in kept_het)
    if not (has_zn and has_mg):
        msg.append("   !! WARNING: a catalytic metal is missing. Do not dock this.")
    print("\n".join(msg) + "\n")


def main():
    files = sys.argv[1:]
    if not files:
        sys.exit("Usage: python clean_structures.py 1XOQ.pdb 7CBQ.pdb 1XMU.pdb")
    for f in files:
        clean(f)
    print("Next: convert each receptor to pdbqt with Meeko, e.g.")
    print("    mk_prepare_receptor.py -i 1XOQ_clean.pdb -o 1XOQ -p")


if __name__ == "__main__":
    main()
