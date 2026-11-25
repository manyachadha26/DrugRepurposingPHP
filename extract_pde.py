"""
Pull proximal-tubule expression for the PDE4 family and its competitors
out of the Human Protein Atlas single-cell-type file.

Needs rna_single_cell_type.tsv in the same folder (unzipped from
rna_single_cell_type.tsv.zip on the HPA single cell type data page).

Run:  conda activate phppde4 && python extract_pde.py
"""

import sys
import pandas as pd

TSV = "rna_single_cell_type.tsv"

GENES = {
    "PDE4A": "PDE4", "PDE4B": "PDE4", "PDE4C": "PDE4", "PDE4D": "PDE4",
    "PDE1A": "Competitor", "PDE1B": "Competitor", "PDE1C": "Competitor",
    "PDE2A": "Competitor", "PDE3A": "Competitor", "PDE3B": "Competitor",
    "PDE7A": "Competitor", "PDE8A": "Competitor", "PDE8B": "Competitor",
}

# The nephron, in anatomical order. Proximal tubule is where PTH acts and
# where Gsa is imprinted, so it is the one that decides the project.
KIDNEY_CELLS = [
    "proximal tubule cells",
    "loop of henle epithelial cells",
    "distal convoluted tubule cells",
    "renal connecting tubule cells",
    "renal collecting duct principal cells",
    "renal collecting duct intercalated cells",
    "podocytes",
]
KEY = "proximal tubule cells"

try:
    df = pd.read_csv(TSV, sep="\t")
except FileNotFoundError:
    sys.exit(f"Could not find {TSV}. Unzip it next to this script first.")

print("Columns in the file:", list(df.columns))

name_cols = [c for c in df.columns if "name" in c.lower() and "gene" in c.lower()]
gene_cols = [c for c in df.columns if c.lower() in ("gene", "gene_name", "gene name")]
gene_col = name_cols[0] if name_cols else gene_cols[0]
cell_col = next(c for c in df.columns if "cell type" in c.lower())
val_col = next(c for c in df.columns if "ncpm" in c.lower() or "ntpm" in c.lower())

sub = df[df[gene_col].isin(GENES)].copy()
if sub.empty:
    for alt in df.columns:
        if alt != gene_col and df[alt].isin(GENES).any():
            gene_col = alt
            sub = df[df[alt].isin(GENES)].copy()
            break
    else:
        sys.exit("No PDE genes matched in any column.")
print(f"Using: gene={gene_col!r}  cell type={cell_col!r}  value={val_col!r}")

# HPA stores cell types lowercase; match case-insensitively so capitalisation
# differences between releases cannot silently empty the table.
sub["_cell"] = sub[cell_col].str.strip().str.lower()
present = set(sub["_cell"])
missing = [c for c in KIDNEY_CELLS if c not in present]
for c in missing:
    print(f"WARNING: cell type {c!r} not found - skipped")
wanted = [c for c in KIDNEY_CELLS if c in present]
if not wanted:
    sys.exit("None of the kidney cell types were found. Check the file.")

sub = sub[sub["_cell"].isin(wanted)]
wide = sub.pivot_table(index=gene_col, columns="_cell", values=val_col, aggfunc="max")
wide = wide.reindex(columns=wanted)
wide.insert(0, "Group", [GENES[g] for g in wide.index])

if KEY in wide.columns:
    wide = wide.sort_values(KEY, ascending=False)
    top = wide.index[0]
    best4 = wide[wide["Group"] == "PDE4"][KEY].max()
    bestc = wide[wide["Group"] == "Competitor"][KEY].max()
    best4_gene = wide[(wide["Group"] == "PDE4") & (wide[KEY] == best4)].index[0]
    bestc_gene = wide[(wide["Group"] == "Competitor") & (wide[KEY] == bestc)].index[0]
    print("\n" + "=" * 62)
    print(f"Highest in proximal tubule:  {top}  ({wide.loc[top, KEY]:.1f} nCPM)")
    print(f"Best PDE4 gene:              {best4_gene}  ({best4:.1f})")
    print(f"Best competitor:             {bestc_gene}  ({bestc:.1f})")
    if best4 >= bestc:
        print("VERDICT: a PDE4 gene leads. Rationale holds - carry on.")
    elif best4 >= 0.5 * bestc:
        print("VERDICT: PDE4 is present but not dominant. Project works;")
        print("         discuss the competitor and consider adding it as a target.")
    else:
        print(f"VERDICT: {bestc_gene} is far above every PDE4 gene.")
        print("         Raise this with Dr. Upadhayay this week.")
    print("=" * 62 + "\n")

wide.round(1).to_csv("pde_proximal_tubule.csv")
print("Written: pde_proximal_tubule.csv\n")
print(wide.round(1).to_string())
