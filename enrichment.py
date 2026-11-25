"""
Gate 2, part 3: does the docking actually separate real PDE4 inhibitors
from look-alike molecules that do not work?

This is the test that decides whether any later ranking means anything.
A protocol that cannot pull rolipram above its property-matched decoys cannot
be trusted to judge a diabetes drug it has never seen.

Run after run_gate2.sh:   python enrichment.py
Writes gate2_results.csv and prints the verdict.
"""

import glob
import math
import os
import sys

def best_score(pdbqt):
    """First REMARK VINA RESULT line = the top-ranked pose's affinity."""
    try:
        for line in open(pdbqt):
            if line.startswith("REMARK VINA RESULT"):
                return float(line.split()[3])
    except OSError:
        pass
    return None

ADIR = sys.argv[1] if len(sys.argv) > 1 else "actives"
DDIR = sys.argv[2] if len(sys.argv) > 2 else "docked"
actives = {os.path.basename(f)[:-6] for f in glob.glob(f"{ADIR}/*.pdbqt")}
if not actives:
    sys.exit(f"No {ADIR}/ folder found.")

rows = []
for f in glob.glob(f"{DDIR}/*.pdbqt"):
    name = os.path.basename(f)[:-6]
    s = best_score(f)
    if s is not None:
        rows.append((name, s, 1 if name in actives else 0))

if not rows:
    sys.exit(f"No results in {DDIR}/.")

missing = len(actives) - sum(r[2] for r in rows)
if missing:
    print(f"WARNING: {missing} active(s) have no docked result.\n")

# more negative = better, so sort ascending
rows.sort(key=lambda r: r[1])
N = len(rows)
n_act = sum(r[2] for r in rows)
if n_act == 0:
    sys.exit("No actives among the docked results - check the folder names.")
n_dec = N - n_act

print(f"Docked: {N}   actives: {n_act}   decoys: {n_dec}")
print(f"Ratio: 1 active per {n_dec/n_act:.0f} decoys\n")

ranks = [i for i, r in enumerate(rows, start=1) if r[2] == 1]

def EF(frac):
    k = max(1, int(round(N * frac)))
    hits = sum(1 for i in range(k) if rows[i][2] == 1)
    expected = n_act * k / N
    return hits, k, (hits / expected if expected else float("nan"))

# AUC-ROC by rank sum (Mann-Whitney), no sklearn needed
rank_sum = sum(ranks)
auc = 1.0 - (rank_sum - n_act * (n_act + 1) / 2) / (n_act * n_dec)

# BEDROC, Truchon & Bayly 2007, alpha = 20
alpha = 20.0
ra = n_act / N
s = sum(math.exp(-alpha * i / N) for i in ranks)
rie_den = (ra * (1 - math.exp(-alpha)) / (math.exp(alpha / N) - 1))
rie = s / rie_den if rie_den else float("nan")
bedroc = (rie * ra * math.sinh(alpha / 2) /
          (math.cosh(alpha / 2) - math.cosh(alpha / 2 - alpha * ra)) +
          1 / (1 - math.exp(alpha * (1 - ra))))

print("Where each active landed (1 = best score of all):")
for name, sc, lab in rows:
    if lab == 1:
        print(f"   rank {rows.index((name, sc, lab))+1:4d} of {N}   {sc:7.3f}   {name}")
print()

for frac, label in [(0.01, "1%"), (0.05, "5%"), (0.10, "10%")]:
    hits, k, ef = EF(frac)
    print(f"Enrichment at {label:>3s}: {hits} of {n_act} actives in the top {k}"
          f"   EF = {ef:.1f}x")
print(f"\nAUC-ROC: {auc:.3f}")
print(f"BEDROC (alpha=20): {bedroc:.3f}")

with open("gate2_results.csv", "w") as fh:
    fh.write("rank,name,vina_score,is_active\n")
    for i, (name, sc, lab) in enumerate(rows, start=1):
        fh.write(f"{i},{name},{sc},{lab}\n")

print("\n" + "=" * 62)
if auc >= 0.8:
    print(f"PASS. AUC {auc:.3f} - the protocol separates real inhibitors from")
    print("decoys well. Proceed to the screen.")
elif auc >= 0.7:
    print(f"ACCEPTABLE. AUC {auc:.3f} - usable, but say so plainly in the paper")
    print("and treat borderline hits with caution.")
elif auc >= 0.6:
    print(f"WEAK. AUC {auc:.3f} - barely better than chance. Rankings from this")
    print("protocol will not support strong claims. Consider rescoring the top")
    print("hits with a second method before drawing conclusions.")
else:
    print(f"FAIL. AUC {auc:.3f} - no better than random. Do not screen with this.")
    print("Check the receptor kept its metals, the box is on the pocket, and")
    print("the ligand protonation states are right.")
print("=" * 62)
print("\nReport AUC, EF at 1% and BEDROC together. AUC describes the whole")
print("ranking; EF and BEDROC describe the top, which is all a screen uses.")
