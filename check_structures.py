"""
Step 3 structure checks, computed from the downloaded PDB files themselves.

Reports per structure:
  - chains present, and how many copies of the protein
  - the co-crystallised ligand and every metal ion
  - residues within 5 A of that ligand, computed from coordinates
  - the invariant glutamine, with ITS NUMBER IN THIS STRUCTURE
  - residues missing from the model (REMARK 465), flagged if near the pocket
  - alternate conformations

Why computed and not read off RCSB: the SITE records in a PDB file are written
by whoever deposited it. Different depositors list different residues, so two
structures of the same protein can look like they have different pockets when
they do not. Distances do not have that problem.

Put the four .pdb files in the same folder and run:
    conda activate phppde4 && python check_structures.py
"""

import glob
import math
import os
import sys
from collections import defaultdict, OrderedDict

CUTOFF = 5.0          # angstroms, standard for "lining the pocket"
METALS = {"ZN", "MG", "MN", "CA", "FE", "NA", "K", "CO", "NI", "CU"}
SOLVENT = {"HOH", "WAT", "DOD"}
# common crystallisation additives - present because of how the crystal was
# grown, not because they mean anything biologically
CRYO = {"EDO", "GOL", "PEG", "SO4", "PO4", "ACT", "CL", "MPD", "TRS", "DMS",
        "IOD", "BR", "FMT", "ACY", "NO3", "EPE", "MES"}

AA3 = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU",
       "LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL","MSE","SEC","PYL"}


def parse(path):
    atoms, het, missing, altlocs = [], [], [], set()
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            rec = line[:6]
            if line.startswith("REMARK 465"):
                body = line[10:].strip()
                parts = body.split()
                # data lines look like:  ALA A   123   (optionally with a model no.)
                if len(parts) >= 3 and parts[0] in AA3 and len(parts[1]) == 1:
                    try:
                        missing.append((parts[0], parts[1], int(parts[2])))
                    except ValueError:
                        pass
                continue
            if rec not in ("ATOM  ", "HETATM"):
                continue
            try:
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            except ValueError:
                continue
            alt = line[16].strip()
            resn = line[17:20].strip()
            chain = line[21].strip() or "_"
            try:
                resi = int(line[22:26])
            except ValueError:
                continue
            if alt:
                altlocs.add((chain, resi, resn))
            entry = {"resn": resn, "chain": chain, "resi": resi, "xyz": (x, y, z),
                     "name": line[12:16].strip()}
            (atoms if rec == "ATOM  " else het).append(entry)
    return atoms, het, missing, altlocs


def centre_and_contacts(atoms, lig_atoms):
    """Residues with any atom within CUTOFF of any ligand atom."""
    near = {}
    for la in lig_atoms:
        lx, ly, lz = la["xyz"]
        for a in atoms:
            ax, ay, az = a["xyz"]
            if abs(ax - lx) > CUTOFF or abs(ay - ly) > CUTOFF or abs(az - lz) > CUTOFF:
                continue
            d = math.sqrt((ax - lx) ** 2 + (ay - ly) ** 2 + (az - lz) ** 2)
            if d <= CUTOFF:
                key = (a["chain"], a["resi"], a["resn"])
                if key not in near or d < near[key]:
                    near[key] = d
    return near


def report(path):
    name = os.path.basename(path)
    pdbid = os.path.splitext(name)[0].upper()
    atoms, het, missing, altlocs = parse(path)
    if not atoms:
        print(f"\n### {pdbid}: no ATOM records found - is this really a PDB file?")
        return None

    chains = sorted({a["chain"] for a in atoms})
    per_chain = {c: len({a["resi"] for a in atoms if a["chain"] == c}) for c in chains}

    groups = defaultdict(list)
    for h in het:
        if h["resn"] in SOLVENT:
            continue
        groups[(h["chain"], h["resi"], h["resn"])].append(h)

    metals = {k: v for k, v in groups.items() if k[2] in METALS}
    cryo = {k: v for k, v in groups.items() if k[2] in CRYO}
    ligands = {k: v for k, v in groups.items()
               if k[2] not in METALS and k[2] not in CRYO}

    print("\n" + "=" * 68)
    print(f"### {pdbid}   ({name})")
    print("=" * 68)
    print(f"Chains: {', '.join(chains)}  ({len(chains)} copies of the protein)")
    if len(chains) > 1:
        print("   -> keep ONE chain before docking:  remove not chain A")
    print("   residues per chain: " + ", ".join(f"{c}={n}" for c, n in per_chain.items()))

    mcount = defaultdict(int)
    for (_, _, rn) in metals:
        mcount[rn] += 1
    if mcount:
        print("Metal ions: " + ", ".join(f"{k} x{v}" for k, v in sorted(mcount.items())))
        if not ({"ZN", "MG"} & set(mcount)):
            print("   !! PDE4 needs a zinc AND a magnesium. Check before you strip anything.")
    else:
        print("Metal ions: NONE FOUND  !! PDE4 cannot be docked without its metals.")

    if cryo:
        cc = defaultdict(int)
        for (_, _, rn) in cryo:
            cc[rn] += 1
        print("Crystallisation additives to delete: "
              + ", ".join(f"{k} x{v}" for k, v in sorted(cc.items())))

    if not ligands:
        print("Co-crystallised ligand: NONE FOUND")
        return None

    # the real ligand = the largest non-metal, non-additive group
    key = max(ligands, key=lambda k: len(ligands[k]))
    lig = ligands[key]
    print(f"Co-crystallised ligand: {key[2]}  (chain {key[0]}, residue {key[1]}, "
          f"{len(lig)} atoms)")
    if len(ligands) > 1:
        others = ", ".join(sorted({k[2] for k in ligands if k != key}))
        print(f"   other non-solvent groups present: {others}")

    cx = sum(a["xyz"][0] for a in lig) / len(lig)
    cy = sum(a["xyz"][1] for a in lig) / len(lig)
    cz = sum(a["xyz"][2] for a in lig) / len(lig)
    print(f"GRID BOX CENTRE (ligand centroid):  "
          f"center_x = {cx:.3f}   center_y = {cy:.3f}   center_z = {cz:.3f}")

    same_chain = [a for a in atoms if a["chain"] == key[0]]
    near = centre_and_contacts(same_chain, lig)
    ordered = OrderedDict(sorted(near.items(), key=lambda kv: kv[0][1]))
    print(f"\nPocket residues within {CUTOFF:.0f} A of the ligand "
          f"({len(ordered)} residues, chain {key[0]}):")
    print("   " + ", ".join(f"{rn}{ri}" for (_, ri, rn) in ordered))

    glns = [(ri, d) for (_, ri, rn), d in near.items() if rn == "GLN"]
    if glns:
        ri, d = min(glns, key=lambda t: t[1])
        print(f"\n>>> INVARIANT GLUTAMINE IN THIS STRUCTURE: Gln{ri}  "
              f"(closest glutamine, {d:.2f} A from the ligand)")
        print("    Use THIS number in your pose filter for THIS structure.")
    else:
        print("\n!! No glutamine within the cutoff - unexpected for PDE4. Check by eye.")

    if missing:
        in_chain = [m for m in missing if m[1] == key[0]]
        print(f"\nResidues missing from the model (REMARK 465): "
              f"{len(missing)} total, {len(in_chain)} in chain {key[0]}")
        pocket_nums = {ri for (_, ri, _) in ordered}
        if pocket_nums:
            lo, hi = min(pocket_nums) - 8, max(pocket_nums) + 8
            nearby = [f"{rn}{ri}" for (rn, ch, ri) in in_chain if lo <= ri <= hi]
            if nearby:
                print("   !! MISSING NEAR THE POCKET: " + ", ".join(nearby))
                print("   Decide whether to model these back in, and say so in methods.")
            else:
                print("   None of them are near the pocket - fine to proceed.")
    else:
        print("\nResidues missing from the model: none listed.")

    if altlocs:
        print(f"\nAlternate conformations on {len(altlocs)} residues - "
              "keep conformation A only.")

    return {"pdb": pdbid, "gln": (min(glns, key=lambda t: t[1])[0] if glns else None),
            "centre": (cx, cy, cz), "ligand": key[2], "metals": dict(mcount),
            "nres": len(ordered)}


def main():
    files = sorted(glob.glob("*.pdb"))
    files = [f for f in files if "native_ligand" not in f.lower()]
    if not files:
        sys.exit("No .pdb files here. Put the four downloaded structures in this folder.")
    print(f"Checking {len(files)} structure(s): {', '.join(files)}")
    results = [r for r in (report(f) for f in files) if r]

    if len(results) > 1:
        print("\n" + "=" * 68)
        print("### THE NUMBERING WARNING - read this before writing any filter")
        print("=" * 68)
        for r in results:
            print(f"   {r['pdb']}: invariant glutamine = Gln{r['gln']}")
        nums = {r["gln"] for r in results}
        if len(nums) > 1:
            print("\n   These numbers are DIFFERENT. The same glutamine, doing the same")
            print("   job, is numbered differently in each structure, because the")
            print("   crystallographers used different constructs.")
            print("\n   A pose filter that hardcodes one number will silently pass or")
            print("   fail everything in the other structures. Key your filter on the")
            print("   per-structure number above, never on a single global one.")
    print()


if __name__ == "__main__":
    main()
