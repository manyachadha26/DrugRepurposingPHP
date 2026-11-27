import sys, math
from collections import defaultdict

METALS = {"ZN","MG","MN","FE","CA","NI","CO","CU","CD"}
JUNK = {"HOH","GOL","EDO","SO4","PO4","CL","NA","K","ACT","MES","TRS",
        "DMS","PEG","IMD","BME","FMT","ACY","NO3","CAC","TRIS","MPD","IPA"}

for path in sys.argv[1:]:
    metals, groups = [], defaultdict(list)
    for ln in open(path):
        if not ln.startswith("HETATM"):
            continue
        res = ln[17:20].strip()
        ch, seq = ln[21], ln[22:26].strip()
        xyz = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
        if res in METALS:
            metals.append((res + seq, xyz))
        elif res not in JUNK:
            groups[(res, ch, seq)].append(xyz)
    name = path.split("/")[-1]
    if not groups:
        print(f"{name}: no ligand found"); continue
    key = max(groups, key=lambda k: len(groups[k]))
    lig = groups[key]
    print(f"\n{name}")
    print(f"  ligand      : {key[0]} chain {key[1]} {key[2]}, {len(lig)} atoms")
    print(f"  metals      : {', '.join(m for m,_ in metals) if metals else 'NONE'}")
    for mname, mxyz in metals:
        d = min(math.dist(a, mxyz) for a in lig)
        verdict = "CONTACTS METAL" if d <= 3.0 else "no metal contact"
        print(f"  d(lig,{mname:<7}) = {d:5.2f} A   {verdict}")
