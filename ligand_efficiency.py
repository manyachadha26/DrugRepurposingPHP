"""
Does normalising for molecular size fix the ranking?

Vina scores correlate strongly with how many atoms a molecule has, so a big
weak molecule can outscore a small tight-binding one. Ligand efficiency
divides the score by the number of heavy atoms, which removes most of that.
It is a standard medicinal-chemistry measure, not a trick.

    LE = score / heavy atom count       (more negative is still better)

If enrichment improves under LE, rank your screen by LE and say so in the
methods. If it does not, the problem is the scoring function itself and
size was not the main issue.

    python ligand_efficiency.py
"""
import csv, math, os, sys
from collections import defaultdict
import statistics as st
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")

if not os.path.exists("gate2_results.csv"):
    sys.exit("Run enrichment.py first.")

smi = {}
matched = {}
for r in csv.DictReader(open("decoys.csv")):
    smi[r["name"]] = r["smiles"]; matched[r["name"]] = r["matched_active"]
for r in csv.DictReader(open("controls_dude_protomers.csv")):
    smi[r["name"]] = r["smiles"]

rows = []
for r in csv.DictReader(open("gate2_results.csv")):
    n = r["name"]
    if n not in smi: continue
    m = Chem.MolFromSmiles(smi[n])
    if not m: continue
    ha = m.GetNumHeavyAtoms()
    s = float(r["vina_score"])
    rows.append((n, s, s/ha, ha, int(r["is_active"])))

def stats(data, key, label):
    data = sorted(data, key=key)
    N = len(data); na = sum(1 for d in data if d[4] == 1); nd = N - na
    if na == 0 or nd == 0: return None
    ranks = [i for i, d in enumerate(data, 1) if d[4] == 1]
    auc = 1.0 - (sum(ranks) - na*(na+1)/2)/(na*nd)
    k = max(1, int(round(N*0.01)))
    hits = sum(1 for i in range(k) if data[i][4] == 1)
    ef = hits/(na*k/N) if na*k/N else 0
    alpha = 20.0; ra = na/N
    s_ = sum(math.exp(-alpha*i/N) for i in ranks)
    den = ra*(1-math.exp(-alpha))/(math.exp(alpha/N)-1)
    rie = s_/den if den else 0
    bed = (rie*ra*math.sinh(alpha/2)/(math.cosh(alpha/2)-math.cosh(alpha/2-alpha*ra))
           + 1/(1-math.exp(alpha*(1-ra))))
    print(f"{label:28s} AUC {auc:.3f}   EF1% {ef:5.1f}x   BEDROC {bed:.3f}")
    return auc

POTENT = {"roflumilast", "rolipram", "apremilast", "cilomilast_anion"}

print("=" * 72)
print("ALL SEVEN ACTIVES, pooled")
print("=" * 72)
stats(rows, lambda d: d[1], "raw Vina score")
stats(rows, lambda d: d[2], "ligand efficiency")

print("\n" + "=" * 72)
print("FOUR POTENT ACTIVES ONLY, pooled")
print("=" * 72)
sub = [d for d in rows if d[4] == 0 or d[0] in POTENT]
sub = [(d[0], d[1], d[2], d[3], 1 if d[0] in POTENT else 0) for d in sub]
stats(sub, lambda d: d[1], "raw Vina score")
stats(sub, lambda d: d[2], "ligand efficiency")

print("\n" + "=" * 72)
print("EACH ACTIVE IN ITS OWN DECOY SET, by ligand efficiency")
print("=" * 72)
groups = defaultdict(list)
for n, a in matched.items():
    groups[a].append(n)
byname = {d[0]: d for d in rows}
fr = []
print(f"{'active':20s} {'rank (raw)':>12s} {'rank (LE)':>12s}")
for a in sorted(d[0] for d in rows if d[4] == 1):
    g = groups.get(a) or groups.get(a.replace("_anion", "")) or []
    g = [x for x in g if x in byname]
    if not g: continue
    ar = byname[a]
    r_raw = sum(1 for x in g if byname[x][1] < ar[1]) + 1
    r_le = sum(1 for x in g if byname[x][2] < ar[2]) + 1
    n = len(g) + 1
    fr.append(r_le/n)
    print(f"{a:20s} {r_raw:5d}/{n:<6d} {r_le:5d}/{n:<6d}")
if fr:
    print(f"\nmean fractional rank by LE: {st.mean(fr):.3f}   (0.5 = random)")

xs = [d[3] for d in rows if d[4] == 0]
ys = [d[2] for d in rows if d[4] == 0]
mx, my = st.mean(xs), st.mean(ys)
num = sum((a-mx)*(b-my) for a, b in zip(xs, ys))
den = math.sqrt(sum((a-mx)**2 for a in xs)*sum((b-my)**2 for b in ys))
print(f"\nsize correlation after LE: r = {num/den if den else 0:+.3f}  "
      f"(was -0.739 on raw score)")
