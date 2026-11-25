"""
Re-pick the actives from the already-downloaded ChEMBL data.
No re-download - reads chembl_raw.pkl.

Why: taking the single most potent compound from each scaffold gives you 40
picomolar clinical candidates. That is not a fair test. Your real question is
whether docking can pick out a MODEST binder, because a repurposed diabetes
drug will never be picomolar. This samples across the potency range instead.

    python select_actives.py
"""
import math, os, pickle, random, sys
from collections import defaultdict
import statistics as st

N_WANTED = 40
PIC50_MIN, PIC50_MAX = 7.0, 10.0   # 100 nM to 0.1 nM - real but not freakish
MW_RANGE = (250, 550)
SEED = 42

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, inchi
from rdkit.Chem.Scaffolds import MurckoScaffold
RDLogger.DisableLog("rdApp.*")
random.seed(SEED)

if not os.path.exists("chembl_raw.pkl"):
    sys.exit("chembl_raw.pkl missing - run get_chembl_actives.py first.")
raw = pickle.load(open("chembl_raw.pkl", "rb"))
print(f"Records in cache: {len(raw)}")

best = {}
for r in raw:
    try:
        v = float(r.get("standard_value"))
    except (TypeError, ValueError):
        continue
    if v <= 0 or r.get("standard_relation") not in (None, "=", "'='"):
        continue
    cs = r.get("confidence_score")
    if cs is not None and int(cs) < 8:
        continue
    smi, cid = r.get("canonical_smiles"), r.get("molecule_chembl_id")
    if not smi or not cid:
        continue
    p = 9.0 - math.log10(v)
    if cid not in best or p > best[cid][1]:
        best[cid] = (smi, p)

ctrl = set()
if os.path.exists("controls.smi"):
    for line in open("controls.smi"):
        m = Chem.MolFromSmiles(line.split()[0])
        if m:
            ctrl.add(inchi.MolToInchiKey(m))

cands, seen = [], set()
for cid, (smi, p) in best.items():
    if not (PIC50_MIN <= p <= PIC50_MAX):
        continue
    m = Chem.MolFromSmiles(smi)
    if m is None or "." in smi:
        continue
    mw = Descriptors.MolWt(m)
    if not (MW_RANGE[0] <= mw <= MW_RANGE[1]):
        continue
    k = inchi.MolToInchiKey(m)
    if k in seen or k in ctrl:
        continue
    seen.add(k)
    try:
        s = MurckoScaffold.MurckoScaffoldSmiles(mol=m)
    except Exception:
        s = "none"
    cands.append(dict(id=cid, smiles=Chem.MolToSmiles(m), pic50=p, mw=mw, scaf=s))

print(f"In range pIC50 {PIC50_MIN}-{PIC50_MAX}, MW {MW_RANGE[0]}-{MW_RANGE[1]}: {len(cands)}")

# one compound per scaffold, then spread those scaffolds across potency bands
by_scaf = defaultdict(list)
for c in cands:
    by_scaf[c["scaf"]].append(c)
reps = [random.choice(v) for v in by_scaf.values()]   # not the most potent - a random member
print(f"Distinct scaffolds available: {len(reps)}")

BANDS = [(7.0, 7.5), (7.5, 8.0), (8.0, 8.5), (8.5, 9.0), (9.0, 10.0)]
bins = defaultdict(list)
for c in reps:
    for lo, hi in BANDS:
        if lo <= c["pic50"] < hi or (hi == 10.0 and c["pic50"] == 10.0):
            bins[(lo, hi)].append(c); break
for b in bins:
    random.shuffle(bins[b])

picked, i = [], 0
per_band = N_WANTED // len(BANDS)
while len(picked) < N_WANTED:
    progressed = False
    for b in BANDS:
        if len(bins[b]) > i and len(picked) < N_WANTED:
            picked.append(bins[b][i]); progressed = True
    if not progressed:
        break
    i += 1

print(f"\nSelected {len(picked)} actives from {len({c['scaf'] for c in picked})} scaffolds")
print("\npotency spread:")
for lo, hi in BANDS:
    n = sum(1 for c in picked if lo <= c["pic50"] < hi or (hi == 10.0 and c["pic50"] == 10.0))
    ic = f"{10**(9-hi):.2g}-{10**(9-lo):.2g} nM"
    print(f"   pIC50 {lo}-{hi}  ({ic:>16s}): {n:3d}  {'#'*n}")
mws = [c["mw"] for c in picked]
print(f"\nMW: {min(mws):.0f} to {max(mws):.0f}, median {st.median(mws):.0f}")
print(f"pIC50: {min(c['pic50'] for c in picked):.2f} to {max(c['pic50'] for c in picked):.2f}, "
      f"median {st.median(c['pic50'] for c in picked):.2f}")

with open("actives_chembl.csv", "w") as f:
    f.write("name,chembl_id,smiles,pic50,mw\n")
    for c in sorted(picked, key=lambda c: -c["pic50"]):
        f.write(f'{c["id"]},{c["id"]},"{c["smiles"]}",{c["pic50"]:.2f},{c["mw"]:.1f}\n')
with open("actives_chembl.smi", "w") as f:
    for c in sorted(picked, key=lambda c: -c["pic50"]):
        f.write(f'{c["smiles"]}\t{c["id"]}\n')
print("\nRewritten: actives_chembl.csv, actives_chembl.smi")
