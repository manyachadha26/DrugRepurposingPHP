"""
Gate 2 round 2: 3D structures and pdbqt files for the ChEMBL actives and
their decoys.

    python prep_gate2_v2.py actives_chembl.smi decoys_chembl.smi

Takes SMILES files (tab or space separated: SMILES<TAB>name), generates 3D
coordinates, and writes pdbqt into actives2/ and decoys2/.

Molecules Open Babel cannot embed are retried with RDKit, and anything that
still fails is reported rather than silently dropped - a missing active would
quietly distort the enrichment number.
"""
import os
import subprocess
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
RDLogger.DisableLog("rdApp.*")

if len(sys.argv) != 3:
    sys.exit("Usage: python prep_gate2_v2.py <actives.smi> <decoys.smi>")

JOBS = [(sys.argv[1], "actives150"), (sys.argv[2], "decoys150")]


def embed(smi, name, path):
    """3D coordinates. Open Babel first, RDKit as the fallback."""
    tmp = path + ".in.smi"
    open(tmp, "w").write(f"{smi}\t{name}\n")
    sdf = path + ".sdf"
    subprocess.run(["obabel", tmp, "-O", sdf, "--gen3d", "--minimize",
                    "--ff", "MMFF94"], capture_output=True, timeout=300)
    os.remove(tmp)
    ok = os.path.exists(sdf) and os.path.getsize(sdf) > 0
    if ok:
        m = next(iter(Chem.SDMolSupplier(sdf, removeHs=False)), None)
        if m is None or m.GetNumConformers() == 0:
            ok = False
    if not ok:
        m = Chem.MolFromSmiles(smi)
        if m is None:
            return None
        m = Chem.AddHs(m)
        ps = AllChem.ETKDGv3(); ps.randomSeed = 42
        if AllChem.EmbedMolecule(m, ps) != 0:
            return None
        AllChem.MMFFOptimizeMolecule(m, maxIters=2000)
        m.SetProp("_Name", name)
        w = Chem.SDWriter(sdf); w.write(m); w.close()
    return sdf


for smi_file, outdir in JOBS:
    if not os.path.exists(smi_file):
        sys.exit(f"Missing {smi_file}")
    os.makedirs(outdir, exist_ok=True)
    entries = []
    for line in open(smi_file):
        p = line.split()
        if len(p) >= 2:
            entries.append((p[0], p[1]))
        elif len(p) == 1:
            entries.append((p[0], f"{outdir}_{len(entries)+1}"))
    print(f"\n{smi_file}: {len(entries)} molecules -> {outdir}/")

    ok = failed = 0
    bad = []
    for i, (smi, name) in enumerate(entries, start=1):
        safe = "".join(c if c.isalnum() or c in "-_" else "_" for c in name)
        base = os.path.join(outdir, safe)
        pdbqt = base + ".pdbqt"
        if os.path.exists(pdbqt) and os.path.getsize(pdbqt) > 0:
            ok += 1
            continue
        sdf = embed(smi, safe, base)
        if sdf is None:
            failed += 1; bad.append(name); continue
        subprocess.run(["mk_prepare_ligand.py", "-i", sdf, "-o", pdbqt],
                       capture_output=True)
        if os.path.exists(sdf):
            os.remove(sdf)
        if os.path.exists(pdbqt) and os.path.getsize(pdbqt) > 0:
            ok += 1
        else:
            failed += 1; bad.append(name)
        if i % 100 == 0:
            print(f"   {i}/{len(entries)}")
    print(f"   done: {ok} written, {failed} failed")
    if bad:
        print(f"   failed: {', '.join(bad[:15])}{' ...' if len(bad) > 15 else ''}")
        if outdir == "actives150":
            print("   NOTE: a missing ACTIVE distorts enrichment. Fix these or")
            print("         remove them from the actives list before docking.")

print("\nNext: edit run_gate2.sh to use actives2/ and decoys2/, then run it.")
