"""
Pull real PDE4 inhibitors from ChEMBL to use as actives in Gate 2.

Seven actives is too few to measure enrichment. This fetches measured PDE4
inhibitors, keeps the potent and well-measured ones, and picks a structurally
DIVERSE subset - ChEMBL's PDE4 data is dominated by a few chemical series, and
40 near-identical rolipram analogues would tell you nothing.

    python get_chembl_actives.py

ChEMBL's server is sometimes down (HTTP 500) or rate limited (429). This
retries with backoff, saves partial progress, and can be rerun safely - it
caches what it already downloaded.

Writes actives_chembl.csv and actives_chembl.smi
"""

import math
import os
import pickle
import sys
import time
from collections import defaultdict

N_WANTED = 40          # how many diverse actives to keep (25 is fine too)
PIC50_MIN = 7.0        # 100 nM or better
MW_RANGE = (200, 600)  # dockable size range
CACHE = "chembl_raw.pkl"
MAX_TRIES = 5

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors, inchi
    from rdkit.Chem.Scaffolds import MurckoScaffold
    RDLogger.DisableLog("rdApp.*")
except ImportError:
    sys.exit("RDKit missing. conda activate phppde4")


def find_targets(new_client):
    """Human PDE4 only. Verified by name - never trust a hardcoded ID."""
    found = {}
    for term in ["phosphodiesterase 4", "PDE4"]:
        try:
            hits = new_client.target.search(term)
        except Exception as e:
            print(f"   search '{term}' failed: {str(e)[:80]}")
            continue
        for r in hits:
            name = (r.get("pref_name") or "")
            low = name.lower()
            if (r.get("organism") == "Homo sapiens"
                    and r.get("target_type") == "SINGLE PROTEIN"
                    and "phosphodiesterase 4" in low
                    and "phosphodiesterase 4" == low.split("'")[-1].strip()[:23].strip()
                    or (r.get("organism") == "Homo sapiens"
                        and r.get("target_type") == "SINGLE PROTEIN"
                        and "phosphodiesterase 4" in low)):
                found[r["target_chembl_id"]] = name
    # final guard: the name must mention PDE4, never anything else
    return {k: v for k, v in found.items() if "phosphodiesterase 4" in v.lower()}


def fetch(res, label):
    """Page through a ChEMBL query, retrying on server errors."""
    out = []
    for attempt in range(1, MAX_TRIES + 1):
        try:
            out = list(res)
            return out
        except Exception as e:
            msg = str(e)[:100].replace("\n", " ")
            wait = min(60, 5 * 2 ** (attempt - 1))
            if attempt == MAX_TRIES:
                print(f"   {label}: gave up after {MAX_TRIES} tries - {msg}")
                return []
            print(f"   {label}: attempt {attempt} failed, retrying in {wait}s")
            time.sleep(wait)
    return out


# ---------- 1. fetch ----------
if os.path.exists(CACHE):
    print(f"Using cached download ({CACHE}). Delete it to re-fetch.")
    raw = pickle.load(open(CACHE, "rb"))
else:
    try:
        from chembl_webresource_client.new_client import new_client
    except ImportError:
        sys.exit("pip install chembl_webresource_client")

    print("Finding human PDE4 targets in ChEMBL...")
    targets = find_targets(new_client)
    if not targets:
        sys.exit("Could not reach ChEMBL or found no PDE4 targets.\n"
                 "Their server returns HTTP 500 when it is down - try again later,\n"
                 "or use the manual route in the message from Claude.")
    for tid, name in targets.items():
        print(f"   {tid}  {name}")

    raw = []
    for tid, name in targets.items():
        print(f"Downloading activities for {tid} ({name})...")
        q = new_client.activity.filter(
            target_chembl_id=tid,
            standard_type__in=["IC50", "Ki"],
            standard_units="nM",
            assay_type="B").only(
            "molecule_chembl_id", "canonical_smiles", "standard_type",
            "standard_value", "standard_relation", "confidence_score",
            "target_chembl_id")
        got = fetch(q, tid)
        print(f"   {len(got)} records")
        for g in got:
            g["_target_name"] = name
        raw.extend(got)
        if got:
            pickle.dump(raw, open(CACHE, "wb"))   # save after each target

    if not raw:
        sys.exit("\nChEMBL returned nothing. Their server is likely down.\n"
                 "Rerun this later - it will pick up where it left off.")
    print(f"Cached to {CACHE}\n")

print(f"Total records: {len(raw)}")

# ---------- 2. filter ----------
best = {}
for r in raw:
    try:
        v = float(r.get("standard_value"))
    except (TypeError, ValueError):
        continue
    if v <= 0:
        continue
    if r.get("standard_relation") not in (None, "=", "'='"):
        continue            # '>' and '<' are limits, not measurements
    cs = r.get("confidence_score")
    if cs is not None and int(cs) < 8:
        continue
    smi, cid = r.get("canonical_smiles"), r.get("molecule_chembl_id")
    if not smi or not cid:
        continue
    pic50 = 9.0 - math.log10(v)
    if cid not in best or pic50 > best[cid][1]:
        best[cid] = (smi, pic50)

print(f"Compounds with a usable measurement: {len(best)}")

cands, seen = [], set()
for cid, (smi, pic50) in best.items():
    if pic50 < PIC50_MIN:
        continue
    m = Chem.MolFromSmiles(smi)
    if m is None or "." in smi:
        continue
    mw = Descriptors.MolWt(m)
    if not (MW_RANGE[0] <= mw <= MW_RANGE[1]):
        continue
    k = inchi.MolToInchiKey(m)
    if k in seen:
        continue
    seen.add(k)
    cands.append(dict(chembl_id=cid, smiles=Chem.MolToSmiles(m),
                      pic50=pic50, mw=mw, mol=m))

print(f"Potent (pIC50 >= {PIC50_MIN}), dockable, unique: {len(cands)}")
if len(cands) < 10:
    sys.exit("Too few actives. Lower PIC50_MIN to 6.5 and rerun.")

# ---------- 3. drop anything already a control ----------
ctrl = set()
if os.path.exists("controls.smi"):
    for line in open("controls.smi"):
        m = Chem.MolFromSmiles(line.split()[0])
        if m:
            ctrl.add(inchi.MolToInchiKey(m))
before = len(cands)
cands = [c for c in cands if inchi.MolToInchiKey(c["mol"]) not in ctrl]
if before != len(cands):
    print(f"Removed {before - len(cands)} already present in controls.smi")

# ---------- 4. diverse subset ----------
scaf = defaultdict(list)
for c in cands:
    try:
        s = MurckoScaffold.MurckoScaffoldSmiles(mol=c["mol"])
    except Exception:
        s = "none"
    scaf[s].append(c)
print(f"Distinct Murcko scaffolds: {len(scaf)}")

for s in scaf:
    scaf[s].sort(key=lambda c: -c["pic50"])
picked, rnd = [], 0
while len(picked) < N_WANTED:
    added = False
    for s in sorted(scaf, key=lambda s: -scaf[s][0]["pic50"]):
        if len(scaf[s]) > rnd:
            picked.append(scaf[s][rnd]); added = True
            if len(picked) >= N_WANTED:
                break
    if not added:
        break
    rnd += 1

nscaf = len({MurckoScaffold.MurckoScaffoldSmiles(mol=c["mol"]) for c in picked})
print(f"Selected: {len(picked)} actives from {nscaf} scaffolds")
print(f"pIC50 range: {min(c['pic50'] for c in picked):.1f} to "
      f"{max(c['pic50'] for c in picked):.1f}")
print(f"MW range: {min(c['mw'] for c in picked):.0f} to "
      f"{max(c['mw'] for c in picked):.0f} Da")

with open("actives_chembl.csv", "w") as f:
    f.write("name,chembl_id,smiles,pic50,mw\n")
    for c in sorted(picked, key=lambda c: -c["pic50"]):
        f.write(f'{c["chembl_id"]},{c["chembl_id"]},"{c["smiles"]}",'
                f'{c["pic50"]:.2f},{c["mw"]:.1f}\n')
with open("actives_chembl.smi", "w") as f:
    for c in sorted(picked, key=lambda c: -c["pic50"]):
        f.write(f'{c["smiles"]}\t{c["chembl_id"]}\n')

print("\nWritten: actives_chembl.csv, actives_chembl.smi")
print("Next: decoys for these, then rerun Gate 2.")
