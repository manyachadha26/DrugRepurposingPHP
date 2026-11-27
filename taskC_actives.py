import random
import pandas as pd
from chembl_webresource_client.new_client import new_client
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")
random.seed(42)
N_WANTED = 100

tgts = new_client.target.filter(pref_name__icontains="phosphodiesterase 10",
                                organism="Homo sapiens")
ids = [t["target_chembl_id"] for t in tgts if t["target_type"] == "SINGLE PROTEIN"]
print("targets:", ids)

rows = []
for tid in ids:
    q = new_client.activity.filter(
        target_chembl_id=tid,
        standard_type__in=["IC50", "Ki", "Kd", "EC50"]).only(
        ["molecule_chembl_id", "canonical_smiles", "standard_value",
         "standard_units", "standard_relation", "standard_type"])
    for a in q:
        rows.append(a)
print("raw records:", len(rows))

df = pd.DataFrame(rows).dropna(subset=["canonical_smiles", "standard_value"])
df = df[df["standard_relation"].fillna("=") == "="]
df = df[df["standard_units"].fillna("nM") == "nM"]
df["standard_value"] = pd.to_numeric(df["standard_value"], errors="coerce")
df = df[df["standard_value"] > 0]
import numpy as np
df["pact"] = 9 - np.log10(df["standard_value"])
df = df[(df.pact >= 7.0) & (df.pact <= 10.0)]
print("after potency filter:", len(df), "rows,",
      df.molecule_chembl_id.nunique(), "compounds")

best = df.groupby("molecule_chembl_id").agg(
    smiles=("canonical_smiles", "first"), pact=("pact", "max")).reset_index()

keep = []
for _, r in best.iterrows():
    m = Chem.MolFromSmiles(r.smiles)
    if m is None:
        continue
    mw = Descriptors.MolWt(m)
    if not (250 <= mw <= 550):
        continue
    sc = Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(m)) or r.smiles
    keep.append((r.molecule_chembl_id, r.smiles, r.pact, mw, sc))
print("after size filter:", len(keep))

by_sc = {}
for rec in keep:
    by_sc.setdefault(rec[4], []).append(rec)
picked = [random.choice(v) for v in by_sc.values()]
print("distinct scaffolds:", len(picked))

picked.sort(key=lambda r: r[2])
if len(picked) > N_WANTED:
    step = len(picked) / N_WANTED
    picked = [picked[int(i * step)] for i in range(N_WANTED)]

with open("taskC/pde10a_actives.smi", "w") as fh:
    for cid, smi, pact, mw, sc in picked:
        fh.write(f"{smi}\t{cid}\n")
print(f"\nwrote taskC/pde10a_actives.smi with {len(picked)} compounds")
print(f"  pAct range {min(p[2] for p in picked):.2f} to {max(p[2] for p in picked):.2f}")
print(f"  MW range  {min(p[3] for p in picked):.0f} to {max(p[3] for p in picked):.0f}")
