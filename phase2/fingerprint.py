#!/usr/bin/env python3
"""
Build an interaction fingerprint for every docked pose.

WHY
---
The docking SCORE cannot rank PDE4D binders (AUC 0.612 Vina, 0.624 Vinardo).
But the POSES are accurate (Gate 1 redocked at 0.77 A). So we stop asking
"what energy did Vina assign" and start asking "does this pose make the
contacts that every real PDE4 inhibitor makes".

Known PDE4 pharmacophore, from the crystallography literature:
  - a hydrogen bond to the invariant glutamine, which anchors the catechol
    ether oxygens of rolipram-like inhibitors
  - the ligand sandwiched in a hydrophobic clamp between two phenylalanines
  - proximity to the bimetal centre (Zn2+ / Mg2+) and its bridging waters

NOTHING IS HARDCODED BY RESIDUE NUMBER. We already hit that trap once: the
invariant Gln is 369 / 443 / 466 depending on the construct, so a hardcoded
number silently zeroes the filter. Here the invariant Gln and the clamp
phenylalanines are FOUND GEOMETRICALLY, from where the native crystal ligand
actually sits. The script prints which residues it picked so you can check
them against the literature before trusting anything downstream.

USAGE
-----
    python fingerprint.py <receptor.pdb> <native_ligand.pdb> \\
                          <docked_folder> <actives_folder> [out.csv]

Example:
    python fingerprint.py 1XOQ_clean.pdb 1XOQ_native_ligand.pdb \\
                          docked2 actives2 fingerprints.csv

OUTPUT
------
fingerprints.csv - one row per pose, columns:
    name, label, hac
    d_gln          A    closest ligand heavy atom to the invariant Gln amide
    gln_hbond      0/1  that distance <= 3.5 A
    d_metal        A    closest ligand heavy atom to Zn or Mg
    metal_shell    0/1  within 6.0 A of a metal (the water-bridged shell)
    n_clamp        int  how many clamp Phe residues are contacted
    n_res          int  distinct residues within 4.0 A
    f_apolar       0-1  fraction of ligand carbons making an apolar contact
    f_polar        0-1  fraction of ligand N/O making a polar contact
    f_buried       0-1  fraction of ligand atoms with >=8 receptor atoms near
    r_<RES><num>   0/1  one column per pocket residue, contacted or not

The f_ columns are FRACTIONS on purpose. Raw counts scale with molecular
size, and a size-driven score is exactly the artifact that invalidated
Gate 2 round 1. Fractions do not.
"""

import os
import sys
import glob
import math
import numpy as np

# --- distance cut-offs, in Angstrom -----------------------------------
CUT_HBOND = 3.5     # donor-acceptor heavy atom distance
CUT_CONTACT = 4.0   # generic residue contact
CUT_APOLAR = 4.0    # carbon-carbon
CUT_CLAMP = 4.5     # slightly looser, stacking is not a point contact
CUT_METAL = 6.0     # metal + its bridging water shell
CUT_POCKET = 8.0    # which residues count as "the pocket" at all
BURIED_N = 8        # receptor atoms within CUT_CLAMP to call an atom buried

METALS = {"ZN", "MG", "MN", "FE", "CA"}
POLAR = {"N", "O"}


def parse_pdb(path, want_hetatm=True):
    """
    Minimal fixed-column PDB parser. Returns a list of dicts.

    We parse by column positions, not by splitting on whitespace: PDB files
    routinely run fields together (a -100.123 coordinate touches its
    neighbour) and split() silently mangles those lines.
    """
    out = []
    with open(path) as fh:
        for line in fh:
            rec = line[:6].strip()
            if rec == "ATOM" or (want_hetatm and rec == "HETATM"):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except ValueError:
                    continue
                name = line[12:16].strip()
                el = line[76:78].strip().upper()
                if not el:
                    # fall back to the atom name: first alphabetic character,
                    # except for two-letter metals which we check whole
                    el = "".join(c for c in name if c.isalpha())[:2].upper()
                    if el not in METALS:
                        el = el[0] if el else "X"
                out.append({
                    "name": name,
                    "resn": line[17:20].strip().upper(),
                    "chain": line[21:22].strip(),
                    "resi": line[22:27].strip(),
                    "el": el,
                    "xyz": (x, y, z),
                })
    return out


def parse_pose(path):
    """
    Read MODEL 1 of a docked pdbqt and return its heavy atoms.

    Only the top-ranked pose matters: that is the one a real screen acts on.
    AutoDock element types sit in columns 77-79; anything starting with H
    is a hydrogen (H, HD) and is dropped.
    """
    atoms, inside, seen_model = [], False, False
    with open(path) as fh:
        for line in fh:
            if line.startswith("MODEL"):
                if seen_model:
                    break
                seen_model, inside = True, True
                continue
            if line.startswith("ENDMDL"):
                break
            if not line.startswith(("ATOM", "HETATM")):
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            t = line[77:79].strip().upper() or line[12:16].strip()[:1].upper()
            if t.startswith("H"):
                continue
            el = "C" if t in ("C", "A") else t[0]
            atoms.append({"el": el, "xyz": (x, y, z)})
    return atoms


def coords(atoms):
    return np.array([a["xyz"] for a in atoms], dtype=float)


def dmat(a, b):
    """Pairwise distances between two coordinate arrays."""
    if len(a) == 0 or len(b) == 0:
        return np.empty((len(a), len(b)))
    return np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2)


def main():
    if len(sys.argv) < 5:
        print(__doc__)
        sys.exit(1)

    rec_path, nat_path = sys.argv[1], sys.argv[2]
    docked_dir, actives_dir = sys.argv[3], sys.argv[4]
    out_csv = sys.argv[5] if len(sys.argv) > 5 else "fingerprints.csv"

    for p in (rec_path, nat_path, docked_dir, actives_dir):
        if not os.path.exists(p):
            sys.exit(f"ERROR: not found: {p}")

    rec = parse_pdb(rec_path)
    nat = [a for a in parse_pdb(nat_path) if not a["el"].startswith("H")]
    if not rec:
        sys.exit(f"ERROR: no atoms parsed from {rec_path}")
    if not nat:
        sys.exit(f"ERROR: no atoms parsed from {nat_path}")
    print(f"Receptor atoms      : {len(rec)}")
    print(f"Native ligand atoms : {len(nat)}")

    rec_xyz = coords(rec)
    nat_xyz = coords(nat)

    # ---------- find the metals ----------
    metal_idx = [i for i, a in enumerate(rec)
                 if a["el"] in METALS or a["resn"] in METALS]
    if metal_idx:
        print("Metals found        : " + ", ".join(
            f"{rec[i]['resn']}{rec[i]['resi']}" for i in metal_idx))
    else:
        print("Metals found        : NONE")
        print("  WARNING: no Zn/Mg in this receptor file. Either the cleaning")
        print("  step dropped them, or you are pointing at the wrong file.")
        print("  The metal columns will be meaningless. Check before using.")

    # ---------- group receptor atoms by residue ----------
    residues = {}
    for i, a in enumerate(rec):
        residues.setdefault((a["chain"], a["resn"], a["resi"]), []).append(i)

    # ---------- which residues line the pocket ----------
    d_rec_nat = dmat(rec_xyz, nat_xyz)
    min_to_nat = d_rec_nat.min(axis=1)
    pocket = sorted(
        (k for k, idx in residues.items()
         if min_to_nat[idx].min() <= CUT_POCKET and k[1] not in METALS
         and k[1] != "HOH"),
        key=lambda k: (k[0], int("".join(c for c in k[2] if c.isdigit() or c == "-") or 0)),
    )
    print(f"Pocket residues     : {len(pocket)} within {CUT_POCKET} A "
          f"of the native ligand")

    # ---------- the invariant glutamine, found geometrically ----------
    gln_atoms = []
    best = None
    for k in pocket:
        if k[1] != "GLN":
            continue
        amide = [i for i in residues[k] if rec[i]["name"] in ("NE2", "OE1")]
        if not amide:
            continue
        d = min_to_nat[amide].min()
        if best is None or d < best[0]:
            best = (d, k, amide)
    if best:
        print(f"Invariant Gln       : {best[1][1]}{best[1][2]} "
              f"chain {best[1][0]}, amide {best[0]:.2f} A from native ligand")
        gln_atoms = best[2]
        if best[0] > 5.0:
            print("  WARNING: that is far from the native ligand. Check it "
                  "is really the invariant Gln before trusting d_gln.")
    else:
        print("Invariant Gln       : NOT FOUND")
        print("  WARNING: no glutamine near the ligand. d_gln and gln_hbond")
        print("  will be blank. Do not use them.")

    # ---------- the hydrophobic clamp phenylalanines ----------
    clamp = [k for k in pocket
             if k[1] == "PHE" and min_to_nat[residues[k]].min() <= 5.0]
    print("Clamp Phe           : " + (", ".join(
        f"{k[1]}{k[2]}" for k in clamp) if clamp else "none within 5 A"))
    print()

    # ---------- labels ----------
    active_names = {
        os.path.basename(f).rsplit(".", 1)[0]
        for f in glob.glob(os.path.join(actives_dir, "*"))
    }

    poses = sorted(glob.glob(os.path.join(docked_dir, "*.pdbqt")))
    if not poses:
        sys.exit(f"ERROR: no .pdbqt files in {docked_dir}")
    print(f"Poses to fingerprint: {len(poses)}")

    res_cols = [f"r_{k[1]}{k[2]}" for k in pocket]
    header = (["name", "label", "hac", "d_gln", "gln_hbond", "d_metal",
               "metal_shell", "n_clamp", "n_res", "f_apolar", "f_polar",
               "f_buried"] + res_cols)

    rows, skipped = [], 0
    for n, path in enumerate(poses, 1):
        name = os.path.basename(path).rsplit(".", 1)[0]
        lig = parse_pose(path)
        if not lig:
            skipped += 1
            continue
        lig_xyz = coords(lig)
        D = dmat(lig_xyz, rec_xyz)                 # ligand x receptor
        lig_el = [a["el"] for a in lig]

        # Gln hydrogen bond: ligand N/O to the Gln amide N/O
        if gln_atoms:
            pol = [i for i, e in enumerate(lig_el) if e in POLAR]
            sub = D[np.ix_(pol, gln_atoms)] if pol else np.empty((0, 0))
            d_gln = float(sub.min()) if sub.size else float("nan")
        else:
            d_gln = float("nan")

        # metal proximity
        if metal_idx:
            d_metal = float(D[:, metal_idx].min())
        else:
            d_metal = float("nan")

        # clamp
        n_clamp = sum(1 for k in clamp
                      if D[:, residues[k]].min() <= CUT_CLAMP)

        # distinct residues contacted
        n_res = sum(1 for k in pocket
                    if D[:, residues[k]].min() <= CUT_CONTACT)

        # apolar fraction: ligand carbons near a receptor carbon
        rec_c = [i for i, a in enumerate(rec) if a["el"] == "C"]
        lig_c = [i for i, e in enumerate(lig_el) if e == "C"]
        if lig_c and rec_c:
            f_apolar = float((D[np.ix_(lig_c, rec_c)].min(axis=1)
                              <= CUT_APOLAR).mean())
        else:
            f_apolar = 0.0

        # polar fraction: ligand N/O near a receptor N/O
        rec_p = [i for i, a in enumerate(rec) if a["el"] in POLAR]
        lig_p = [i for i, e in enumerate(lig_el) if e in POLAR]
        if lig_p and rec_p:
            f_polar = float((D[np.ix_(lig_p, rec_p)].min(axis=1)
                             <= CUT_HBOND).mean())
        else:
            f_polar = 0.0

        # burial
        f_buried = float(((D <= CUT_CLAMP).sum(axis=1) >= BURIED_N).mean())

        row = [name,
               1 if name in active_names else 0,
               len(lig),
               "" if math.isnan(d_gln) else round(d_gln, 3),
               "" if math.isnan(d_gln) else int(d_gln <= CUT_HBOND),
               "" if math.isnan(d_metal) else round(d_metal, 3),
               "" if math.isnan(d_metal) else int(d_metal <= CUT_METAL),
               n_clamp, n_res,
               round(f_apolar, 4), round(f_polar, 4), round(f_buried, 4)]
        row += [int(D[:, residues[k]].min() <= CUT_CONTACT) for k in pocket]
        rows.append(row)

        if n % 200 == 0 or n == len(poses):
            print(f"  {n} / {len(poses)}")

    import pandas as pd
    df = pd.DataFrame(rows, columns=header)
    df.to_csv(out_csv, index=False)

    n_act = int(df.label.sum())
    print(f"\nWrote {out_csv}")
    print(f"  rows {len(df)}   actives {n_act}   decoys {len(df) - n_act}")
    if skipped:
        print(f"  skipped {skipped} files with no parsable MODEL 1")
    print(f"  features {len(header) - 3} "
          f"({len(res_cols)} of them per-residue contacts)")
    print("\nCheck the residues named above against the PDE4 literature "
          "before you trust the numbers.")
    print(f"\nNext:  python fp_classify.py {out_csv} actives_chembl.smi "
          f"decoys_chembl.smi rescore_results.csv")


if __name__ == "__main__":
    main()
