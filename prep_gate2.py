"""
Gate 2, part 1: turn the actives and decoys into pdbqt files Vina can read.

Actives = the SEVEN DUD-E protomers, not the six neutral controls.
Why: DUD-E matched each decoy's charge to a specific protonation state. If you
dock neutral cilomilast against decoys matched to the cilomilast anion, the
charge property is no longer matched and the enrichment number is not clean.
Using the protomers makes actives and decoys consistent by construction, and
gives exactly 50 decoys per active.

Run:  python prep_gate2.py
Needs: controls_dude_protomers.sdf, decoys.sdf, and Meeko installed.
"""

import os
import subprocess
import sys
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")

JOBS = [("controls_dude_protomers.sdf", "actives"), ("decoys.sdf", "decoys")]

for sdf, outdir in JOBS:
    if not os.path.exists(sdf):
        sys.exit(f"Missing {sdf}")
    os.makedirs(outdir, exist_ok=True)
    mols = [m for m in Chem.SDMolSupplier(sdf, removeHs=False) if m]
    print(f"{sdf}: {len(mols)} molecules -> {outdir}/")
    ok = fail = 0
    for i, m in enumerate(mols, start=1):
        name = m.GetProp("_Name") if m.HasProp("_Name") else f"{outdir}_{i}"
        name = "".join(c if c.isalnum() or c in "-_" else "_" for c in name)
        one = os.path.join(outdir, name + ".sdf")
        w = Chem.SDWriter(one); w.write(m); w.close()
        out = os.path.join(outdir, name + ".pdbqt")
        r = subprocess.run(["mk_prepare_ligand.py", "-i", one, "-o", out],
                           capture_output=True, text=True)
        os.remove(one)
        if os.path.exists(out) and os.path.getsize(out) > 0:
            ok += 1
        else:
            fail += 1
            print(f"   FAILED: {name}  {r.stderr.strip()[:120]}")
        if i % 50 == 0:
            print(f"   {i}/{len(mols)}")
    print(f"   done: {ok} written, {fail} failed\n")

print("Next:  bash run_gate2.sh")
