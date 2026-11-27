import glob, math
import numpy as np

M = [(float(l[30:38]), float(l[38:46]), float(l[46:54]))
     for l in open("3HR1_clean.pdb")
     if l.startswith("HETATM") and l[17:20].strip() in ("ZN", "MG")]

def top_pose(path):
    xyz, started = [], False
    for ln in open(path):
        if ln.startswith("MODEL"):
            if started:
                break
            started = True
        elif ln.startswith(("ATOM", "HETATM")) and started:
            if ln[76:78].strip() not in ("H", "HD"):
                xyz.append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
    return xyz

cry = [(float(l[30:38]), float(l[38:46]), float(l[46:54]))
       for l in open("3HR1_native_ligand.pdb") if l.startswith(("ATOM", "HETATM"))]
print("crystal ligand to metal: %.2f A" % min(math.dist(a, m) for a in cry for m in M))

for mode in ("dry", "wet"):
    ds = []
    for f in sorted(glob.glob("taskC/docked_3HR1_%s/*.pdbqt" % mode)):
        p = top_pose(f)
        if p:
            ds.append(min(math.dist(a, m) for a in p for m in M))
    ds = np.array(ds)
    print("%s: n=%d  min %.2f  median %.2f  max %.2f"
          % (mode, len(ds), ds.min(), np.median(ds), ds.max()))
