"""
Build property-matched decoys without DUD-E.

For each active, ask ChEMBL for compounds with similar physical properties,
throw out anything with any reported PDE4 activity, throw out anything
structurally similar to any active, and keep a random sample of what's left.

That is the same logic DUD-E uses: look alike by properties, differ by
structure, presumed inactive. The difference is the source pool (ChEMBL rather
than purchasable ZINC) and that you can see exactly what it did.

    python make_decoys.py

Needs chembl_raw.pkl and actives_150.smi (both already in your folder).
Writes decoys_150.smi and decoy_report_150.csv

Safe to rerun - it caches each active's candidate pool.
"""
import csv, os, pickle, random, sys, time
from collections import defaultdict
import statistics as st

N_DECOYS = 15          # per active
TANIMOTO_MAX = 0.35    # above this a "decoy" is too close to a real active
SEED = 42
POOL_CACHE = "decoy_pool_cache_150.pkl"

from rdkit import Chem, RDLogger, DataStructs
from rdkit.Chem import Descriptors, inchi
from rdkit.Chem import rdFingerprintGenerator as fpg
RDLogger.DisableLog("rdApp.*")
random.seed(SEED)

try:
    from chembl_webresource_client.new_client import new_client
except ImportError:
    sys.exit("pip install chembl_webresource_client")

# ---- the actives ----
actives = []
for line in open("actives_150.smi"):
    p = line.split()
    m = Chem.MolFromSmiles(p[0])
    if m:
        actives.append(dict(id=p[1], smiles=p[0], mol=m,
                            mw=Descriptors.MolWt(m),
                            logp=Descriptors.MolLogP(m),
                            rtb=Descriptors.NumRotatableBonds(m),
                            hbd=Descriptors.NumHDonors(m),
                            hba=Descriptors.NumHAcceptors(m)))
print(f"Actives: {len(actives)}")

# ---- everything with ANY measured PDE4 activity is banned as a decoy ----
banned = set()
if os.path.exists("chembl_raw.pkl"):
    for r in pickle.load(open("chembl_raw.pkl", "rb")):
        if r.get("molecule_chembl_id"):
            banned.add(r["molecule_chembl_id"])
print(f"Compounds excluded for having PDE4 data: {len(banned)}")

gen = fpg.GetMorganGenerator(radius=2, fpSize=2048)
afps = [gen.GetFingerprint(a["mol"]) for a in actives]
akeys = {inchi.MolToInchiKey(a["mol"]) for a in actives}

pool = pickle.load(open(POOL_CACHE, "rb")) if os.path.exists(POOL_CACHE) else {}

def fetch_pool(a):
    """Compounds from ChEMBL whose properties bracket this active."""
    for attempt in range(4):
        try:
            q = new_client.molecule.filter(
                molecule_properties__full_mwt__gte=a["mw"] - 25,
                molecule_properties__full_mwt__lte=a["mw"] + 25,
                molecule_properties__alogp__gte=a["logp"] - 1.0,
                molecule_properties__alogp__lte=a["logp"] + 1.0,
                molecule_properties__rtb__gte=max(0, a["rtb"] - 2),
                molecule_properties__rtb__lte=a["rtb"] + 2,
                molecule_properties__hbd__gte=max(0, a["hbd"] - 1),
                molecule_properties__hbd__lte=a["hbd"] + 1,
            ).only("molecule_chembl_id", "molecule_structures")
            return list(q[:400])
        except Exception as e:
            wait = 5 * 2 ** attempt
            print(f"      ChEMBL error, retrying in {wait}s ({str(e)[:60]})")
            time.sleep(wait)
    return []

report, chosen, used_keys = [], [], set()
for i, a in enumerate(actives, start=1):
    if a["id"] not in pool:
        print(f"[{i}/{len(actives)}] {a['id']}  MW {a['mw']:.0f} logP {a['logp']:.1f} ...")
        pool[a["id"]] = fetch_pool(a)
        pickle.dump(pool, open(POOL_CACHE, "wb"))
    cands = pool[a["id"]]

    keep = []
    for c in cands:
        cid = c.get("molecule_chembl_id")
        ms = c.get("molecule_structures") or {}
        smi = ms.get("canonical_smiles")
        if not cid or not smi or cid in banned or "." in smi:
            continue
        m = Chem.MolFromSmiles(smi)
        if m is None:
            continue
        k = inchi.MolToInchiKey(m)
        if k in akeys or k in used_keys:
            continue
        sim = max(DataStructs.BulkTanimotoSimilarity(gen.GetFingerprint(m), afps))
        if sim >= TANIMOTO_MAX:
            continue
        keep.append((cid, Chem.MolToSmiles(m), m, sim))

    random.shuffle(keep)
    take = keep[:N_DECOYS]
    for cid, smi, m, sim in take:
        used_keys.add(inchi.MolToInchiKey(m))
        chosen.append((f"decoy_{cid}", smi, a["id"], sim))
    report.append((a["id"], len(cands), len(keep), len(take),
                   st.median([Descriptors.MolWt(t[2]) for t in take]) if take else 0,
                   a["mw"]))
    print(f"[{i}/{len(actives)}] {a['id']}: pool {len(cands)} -> usable {len(keep)} "
          f"-> kept {len(take)}")

print(f"\nTotal decoys: {len(chosen)}")
short = [r for r in report if r[3] < N_DECOYS]
if short:
    print(f"{len(short)} actives got fewer than {N_DECOYS} decoys:")
    for r in short[:10]:
        print(f"   {r[0]}: {r[3]}")
    print("   (uneven decoy counts are acceptable - the report records them)")

with open("decoys_150.smi", "w") as f:
    for name, smi, parent, sim in chosen:
        f.write(f"{smi}\t{name}\n")
with open("decoy_report_150.csv", "w") as f:
    w = csv.writer(f)
    w.writerow(["decoy", "smiles", "matched_active", "max_tanimoto_to_any_active"])
    for name, smi, parent, sim in chosen:
        w.writerow([name, smi, parent, f"{sim:.3f}"])

sims = [c[3] for c in chosen]
print(f"\nTanimoto to nearest active: max {max(sims):.3f}, median {st.median(sims):.3f}")
print("(all below 0.35 by construction - decoys look similar by properties,")
print(" not by structure, which is the point)")
print("\nWritten: decoys_150.smi, decoy_report_150.csv")
