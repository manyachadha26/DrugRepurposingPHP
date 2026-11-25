"""
Why did Gate 2 come out weak? Three checks on the results you already have.

    python diagnose_gate2.py

Needs gate2_results.csv (from enrichment.py) and decoys.csv.
"""
import csv, math, os, sys
from collections import defaultdict
import statistics as st

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors
    RDLogger.DisableLog("rdApp.*")
    HAVE_RDKIT = True
except ImportError:
    HAVE_RDKIT = False

if not os.path.exists("gate2_results.csv"):
    sys.exit("Run enrichment.py first.")

scores = {}
labels = {}
for r in csv.DictReader(open("gate2_results.csv")):
    scores[r["name"]] = float(r["vina_score"])
    labels[r["name"]] = int(r["is_active"])

# map each decoy to the active it was matched to
matched = {}
smiles = {}
for r in csv.DictReader(open("decoys.csv")):
    key = r["matched_protomer"]
    matched[r["name"]] = r["matched_active"] if key == "" else r["matched_active"]
    smiles[r["name"]] = r["smiles"]
prot_smiles = {}
if os.path.exists("controls_dude_protomers.csv"):
    for r in csv.DictReader(open("controls_dude_protomers.csv")):
        prot_smiles[r["name"]] = r["smiles"]

# DUD-E's matched_active column uses the base name; protomer rows may differ
def group_of(decoy):
    return matched.get(decoy)

groups = defaultdict(list)
for d, a in matched.items():
    if d in scores:
        groups[a].append(d)

actives = [n for n, l in labels.items() if l == 1]

print("=" * 70)
print("CHECK 1 - each active against ONLY its own 50 matched decoys")
print("=" * 70)
print("Pooling all 7 actives into one list compares a 180 Da molecule against")
print("decoys built for a 460 Da one. Judging each active inside its own")
print("matched set removes that confound.\n")
print(f"{'active':20s} {'rank':>10s} {'score':>9s} {'decoy median':>13s}  verdict")
per = []
for a in sorted(actives):
    base = a
    g = groups.get(base) or groups.get(a.replace("_anion", "")) or []
    if not g:
        print(f"{a:20s}  no matched decoy group found")
        continue
    ds = sorted(scores[d] for d in g)
    rank = sum(1 for v in ds if v < scores[a]) + 1
    n = len(ds) + 1
    frac = rank / n
    per.append((a, rank, n, frac))
    verdict = ("top 2%" if frac <= 0.02 else "top 10%" if frac <= 0.10
               else "top 25%" if frac <= 0.25 else "poor")
    print(f"{a:20s} {rank:4d}/{n:<5d} {scores[a]:9.3f} {st.median(ds):13.3f}  {verdict}")
if per:
    print(f"\nmean fractional rank: {st.mean(p[3] for p in per):.3f}  "
          f"(0.5 = random, lower is better)")
    good = sum(1 for p in per if p[3] <= 0.10)
    print(f"actives in the top 10% of their own decoy set: {good}/{len(per)}")

print("\n" + "=" * 70)
print("CHECK 2 - is Vina just ranking by molecular size?")
print("=" * 70)
if not HAVE_RDKIT:
    print("RDKit not available, skipping.")
else:
    xs, ys = [], []
    for d, smi in smiles.items():
        if d in scores:
            m = Chem.MolFromSmiles(smi)
            if m:
                xs.append(Descriptors.MolWt(m)); ys.append(scores[d])
    if len(xs) > 10:
        mx, my = st.mean(xs), st.mean(ys)
        num = sum((a-mx)*(b-my) for a, b in zip(xs, ys))
        den = math.sqrt(sum((a-mx)**2 for a in xs) * sum((b-my)**2 for b in ys))
        r = num/den if den else 0
        print(f"correlation between molecular weight and Vina score: r = {r:.3f}")
        print(f"(negative = bigger molecules score 'better'; |r| above about 0.6")
        print(f" means the pooled ranking is largely sorting by size, not by fit)")
        print(f"\ndecoy MW range: {min(xs):.0f} to {max(xs):.0f} Da")

print("\n" + "=" * 70)
print("CHECK 3 - the potent actives on their own")
print("=" * 70)
print("Theophylline is roughly 100 uM against PDE4 and zardaverine is modest.")
print("Only four of these seven are nanomolar PDE4 inhibitors. An enrichment")
print("test is meant to use genuine potent binders.\n")
POTENT = {"roflumilast", "rolipram", "apremilast", "cilomilast", "cilomilast_anion"}
sub = [(n, scores[n], labels[n]) for n in scores
       if labels[n] == 0 or n in POTENT]
sub.sort(key=lambda t: t[1])
N = len(sub); na = sum(1 for t in sub if t[2] == 1)
if na:
    ranks = [i for i, t in enumerate(sub, 1) if t[2] == 1]
    nd = N - na
    auc = 1.0 - (sum(ranks) - na*(na+1)/2) / (na*nd)
    k = max(1, int(round(N*0.01)))
    hits = sum(1 for i in range(k) if sub[i][2] == 1)
    ef = hits / (na*k/N) if na*k/N else 0
    print(f"actives used: {na}  decoys: {nd}")
    print(f"AUC-ROC (potent actives only): {auc:.3f}")
    print(f"EF at 1%: {ef:.1f}x   ({hits} of {na} in the top {k})")
    print("\nranks:")
    for i, t in enumerate(sub, 1):
        if t[2] == 1:
            print(f"   {i:4d}/{N}   {t[1]:7.3f}   {t[0]}")
