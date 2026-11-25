"""
check_3d.py - automatic sanity check (and optional repair) of ligands.sdf made by Open Babel.

Usage:
  python check_3d.py library.smi ligands.sdf          # just check
  python check_3d.py library.smi ligands.sdf --fix    # check, then rebuild the bad ones with RDKit

Checks: (1) every input molecule made it into the SDF, (2) no molecule is flat,
(3) no atoms sitting on top of each other, (4) every stereocentre you specified in the
SMILES has the same R/S label in the 3D structure (catches flipped chirality).
"""
import sys, random, shutil
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")
smi_file, sdf_file = sys.argv[1], sys.argv[2]
FIX = "--fix" in sys.argv

ref = {}
for line in open(smi_file):
    p = line.split()
    if len(p) >= 2: ref[p[1]] = p[0]

def problems_of(m, name):
    out = []
    conf = m.GetConformer()
    pos = conf.GetPositions()
    if np.ptp(pos[:, 2]) < 0.05 and m.GetNumAtoms() > 3: out.append("FLAT (all atoms in one plane)")
    d = np.linalg.norm(pos[:, None] - pos[None], axis=-1) + np.eye(len(pos)) * 9
    if d.min() < 0.5: out.append("atoms overlap")
    if name in ref:
        a = Chem.MolFromSmiles(ref[name])
        want = {i: c for i, c in Chem.FindMolChiralCenters(a, includeUnassigned=False, useLegacyImplementation=True)}
        if want:
            m2 = Chem.RemoveHs(m); Chem.AssignStereochemistryFrom3D(m2)
            hit = m2.GetSubstructMatch(a)
            if not hit: out.append("3D structure does not match the SMILES")
            else:
                got = {i: c for i, c in Chem.FindMolChiralCenters(m2, includeUnassigned=False, useLegacyImplementation=True)}
                if any(got.get(hit[i]) != c for i, c in want.items()): out.append("CHIRALITY FLIPPED vs your SMILES")
    return out

mols = {}
for m in Chem.SDMolSupplier(sdf_file, removeHs=False):
    if m is not None: mols[m.GetProp("_Name")] = m

bad = {n: problems_of(m, n) for n, m in mols.items()}
bad = {n: p for n, p in bad.items() if p}
missing = [n for n in ref if n not in mols]
print(f"Input molecules: {len(ref)}   in SDF: {len(mols)}")
if missing: print("MISSING from SDF (Open Babel dropped them):", ", ".join(missing))
for n, p in bad.items(): print(f"  PROBLEM  {n}: {'; '.join(p)}")
if not bad and not missing: print("No automatic problems found.")

if FIX and (bad or missing):
    shutil.copy(sdf_file, sdf_file.replace(".sdf", "_raw.sdf"))
    for n in list(bad) + missing:
        a = Chem.AddHs(Chem.MolFromSmiles(ref[n]))
        ps = AllChem.ETKDGv3(); ps.randomSeed = 7
        if AllChem.EmbedMolecule(a, ps) != 0:
            ps.useRandomCoords = True
            if AllChem.EmbedMolecule(a, ps) != 0: print("  could not rebuild", n); continue
        AllChem.MMFFOptimizeMolecule(a, maxIters=2000)
        a.SetProp("_Name", n); mols[n] = a
        print("  rebuilt", n, "-> remaining problems:", problems_of(a, n) or "none")
    w = Chem.SDWriter(sdf_file)
    for n in ref:
        if n in mols: w.write(mols[n])
    w.close()
    print(f"\nRewrote {sdf_file} (old copy kept as {sdf_file.replace('.sdf','_raw.sdf')}).")

names = [n for n in mols]
random.seed(1)
flagged = [n for n in bad]
rest = [n for n in names if n not in flagged]; random.shuffle(rest)
print("\nOpen these 20 in PyMOL and look at them by eye (flagged ones first):")
print(", ".join((flagged + rest)[:20]))
