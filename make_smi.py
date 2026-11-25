"""
make_smi.py - turn your spreadsheet / CSV / DUD-E decoy file into a clean .smi file for Open Babel.

Usage (run in the folder that holds your files, inside the phppde4 environment):
  python make_smi.py drug_library.xlsx library.smi --sheet library
  python make_smi.py drug_library.xlsx controls.smi --sheet controls
  python make_smi.py decoys.csv decoys.smi --prefix decoy_

Input can be .xlsx, .csv (columns 'name' and 'smiles'), or a plain text file with
'SMILES ID' on each line (this is what DUD-E gives you).
What it does: skips blank SMILES, makes names safe (no spaces/brackets), checks every
SMILES with RDKit, drops exact duplicates, and prints a report of anything wrong.
"""
import argparse, re, sys
import pandas as pd
from rdkit import Chem
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

ap = argparse.ArgumentParser()
ap.add_argument("infile"); ap.add_argument("outfile")
ap.add_argument("--sheet", default=None, help="sheet name if input is .xlsx")
ap.add_argument("--prefix", default="", help="text added to the start of every name, e.g. decoy_")
a = ap.parse_args()

f = a.infile.lower()
if f.endswith((".xlsx", ".xlsm")):
    df = pd.read_excel(a.infile, sheet_name=a.sheet or 0)
else:
    df = pd.read_csv(a.infile, sep=None, engine="python")
    if "smiles" not in [str(c).lower() for c in df.columns]:   # headerless 'SMILES ID' file
        df = pd.read_csv(a.infile, sep=r"\s+", header=None, names=["smiles", "name"], usecols=[0, 1])
df.columns = [str(c).strip().lower() for c in df.columns]
if "name" not in df.columns:
    df["name"] = [f"mol{i+1}" for i in range(len(df))]

ok, seen, problems = [], set(), []
for i, r in df.iterrows():
    name = str(r["name"]).strip() if pd.notna(r["name"]) else ""
    smi = str(r["smiles"]).strip() if pd.notna(r["smiles"]) else ""
    if not name or name.lower() == "nan":
        continue
    if not smi or smi.lower() == "nan":
        problems.append((name, "NO SMILES - go copy it from PubChem/ChEMBL")); continue
    m = Chem.MolFromSmiles(smi)
    if m is None:
        problems.append((name, "RDKit cannot read this SMILES: " + smi)); continue
    safe = a.prefix + re.sub(r"[^A-Za-z0-9_.+-]+", "_", name).strip("_")
    key = Chem.MolToSmiles(m)
    if key in seen:
        problems.append((name, "duplicate structure, skipped")); continue
    seen.add(key)
    ok.append((smi, safe))

with open(a.outfile, "w") as fh:
    for smi, nm in ok:
        fh.write(f"{smi} {nm}\n")
print(f"Wrote {len(ok)} molecules to {a.outfile}")
if problems:
    print(f"\n{len(problems)} problem(s) - fix these in your sheet and re-run:")
    for n, p in problems: print(f"  - {n}: {p}")
