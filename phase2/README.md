# Phase 2 — docking protocol validation

Repurposing antidiabetic drugs as PDE4 inhibitors for pseudohypoparathyroidism.
This folder covers everything after Phase 1: validating whether the docking
protocol can actually rank PDE4D binders, and building a method that can.

Target: PDE4D catalytic domain, PDB **1XOQ**, chain A, with its Zn and Mg.

---

## What Phase 2 established

| # | Finding | Evidence |
|---|---|---|
| 1 | Vina reproduces PDE4D crystal geometry | redocking RMSD 0.756 Å, top-ranked pose |
| 2 | Docking **scores** cannot rank PDE4D binders | AUC 0.604 (Vina) and 0.586 (Vinardo) on 40 actives vs 1,197 matched decoys; confirmed on two independently prepared receptors |
| 3 | Deleting metal-coordinating waters distorts poses | docked actives sat median 2.86 Å from a metal vs 4.00 Å for the crystal ligand; after retaining the waters, 4.27 Å |
| 4 | Redocking validation **cannot detect** finding 3 | 0.770 Å on the dry receptor vs 0.756 Å on the wet one — indistinguishable, because the native ligand never reaches the affected region |
| 5 | Interaction-pattern classification works | **AUC 0.822, 95% CI [0.780, 0.864]** on 146 actives vs 2,246 matched decoys, leave-one-scaffold-out |

Finding 4 is the cautionary one: a protocol validated by redocking alone,
on a metalloenzyme prepared by stripping waters, carries this error with no
way to see it.

---

## Benchmark results in full

| Benchmark | Actives | Decoys | Method | AUC | 95% CI | EF 1% | BEDROC |
|---|---|---|---|---|---|---|---|
| round 1 | 7 | 210 | Vina score | 0.628 | — | — | — |
| round 2, dry receptor | 40 | 1,197 | Vina score | 0.612 | [0.518, 0.706] | 0.0× | 0.122 |
| round 2, dry receptor | 40 | 1,197 | Vinardo score | 0.624 | [0.530, 0.718] | 7.7× | 0.168 |
| round 2, dry receptor | 40 | 1,197 | interaction fingerprint | 0.806 | [0.724, 0.888] | 5.2× | 0.273 |
| round 3, wet receptor | 40 | 1,197 | Vina score | 0.604 | [0.510, 0.698] | 5.2× | 0.157 |
| round 3, wet receptor | 40 | 1,197 | Vinardo score | 0.586 | [0.492, 0.680] | 5.2× | 0.197 |
| round 3, wet receptor | 40 | 1,197 | interaction fingerprint | 0.727 | [0.637, 0.817] | 12.9× | 0.324 |
| **round 4, wet receptor** | **146** | **2,246** | **interaction fingerprint** | **0.822** | **[0.780, 0.864]** | **10.2×** | **0.458** |

Round 1 was invalid: 7 actives, and the score tracked molecular weight at
r = −0.739, so it measured size rather than binding.

Round 4 enrichment, hypergeometric: 15 of 146 in the top 24 (p = 2×10⁻¹³),
43 in the top 120 (p = 2×10⁻²⁴), 73 in the top 239 (p = 3×10⁻³⁸).

Label-permutation control for round 4: AUC 0.494 ± 0.011 over 5 shuffles,
confirming no information crosses the cross-validation folds.

---

## Benchmark construction

- **Actives**: ChEMBL PDE4 activity data, 8,806 records. Censored values
  (`>`, `<`) dropped — they are bounds, not measurements. Percent-inhibition
  rows dropped. pIC50 restricted to 7.0–10.0, MW to 250–550 Da. One compound
  per Bemis-Murcko scaffold, chosen at random within the scaffold group
  rather than by potency, and sampled evenly across five potency bands.
  146 actives on 146 distinct scaffolds.
- **Decoys**: 15 per active, property-matched on MW (±25), logP (±1.0),
  rotatable bonds (±2) and H-bond donors (±1). Every compound with any
  recorded PDE4 activity excluded. Tanimoto to the nearest active capped at
  0.35 — observed max 0.348, median 0.212. 2,246 decoys.
- **Units caveat**: the cached ChEMBL pull has no `standard_units` column.
  Values are taken as nM, which is what ChEMBL stores `standard_value` in
  for IC50/Ki/Kd. This assumption belongs in the methods section.

---

## Why interaction pattern instead of score

The poses are accurate (finding 1) but the energies do not rank (finding 2).
So the classifier uses *where the pose sits* rather than what Vina scored it.

Features per pose: distance to the invariant glutamine and whether it
hydrogen-bonds, contacts with the hydrophobic clamp phenylalanines, distance
to the bimetal centre, fractions of ligand atoms making apolar and polar
contacts, burial, and one binary contact column per pocket residue.

The invariant Gln and clamp Phe are located **geometrically**, from where the
native crystal ligand sits — never by residue number. The invariant Gln is
369, 443 or 466 depending on the construct, and a hardcoded number silently
zeroes the filter. On 1XOQ the script resolves Gln369, Phe340 and Phe372.

Fractions rather than raw counts, because raw counts scale with molecular
size and a size-driven score is the artifact that invalidated round 1.

Validation is **leave-one-scaffold-out**: every compound sharing a scaffold
with the held-out active is held out with it, so each prediction comes from a
model that never saw that chemotype. A random split leaks, because analogues
land on both sides and the model learns the scaffold.

---

## Scripts

| Script | Purpose |
|---|---|
| `clean_with_waters.py` | Prepare the receptor keeping waters within 3.0 Å of a metal; drop bulk waters, buffer and other chains; blank the altloc column |
| `patch_waters_pdbqt.py` | Verify the waters survived conversion to pdbqt; append them typed OA if not |
| `rescore.py` | Re-score existing poses with Vinardo, no re-docking |
| `enrichment2.py` | AUC, EF at 1/5/10%, BEDROC for two score columns side by side |
| `check_vinardo.py` | Four checks for the molecular-size confound |
| `fingerprint.py` | Build the interaction fingerprint for every pose |
| `fp_classify.py` | Leave-one-scaffold-out classification, with a label-permutation control |
| `metal_site.py` | `shell`: metal coordination geometry from the raw PDB. `test`: whether metal-proximal actives are the misranked ones |
| `pose_filter.py` | Test a binary H-bond-plus-clamp filter, separating what the filter contributes from what the score contributes |
| `select_more_actives.py` | Select additional actives from the cached ChEMBL pool, offline |
| `screen_library.py` | Rank a compound library by interaction pattern, with the interactions shown alongside |

---

## Reproducing Phase 2

```bash
conda activate phppde4

# receptor, keeping the metal-coordinating waters
python clean_with_waters.py raw/1XOQ.pdb
mk_prepare_receptor.py -i 1XOQ_wet.pdb -o 1XOQ_wet -p
python patch_waters_pdbqt.py 1XOQ_wet.pdb 1XOQ_wet.pdbqt

# gate 1 — pose reproduction
vina --receptor 1XOQ_wet.pdbqt --ligand 1XOQ_native_ligand.pdbqt \
     --config config_1XOQ.txt --out redock_1XOQ_wet.pdbqt
python calc_rmsd.py 1XOQ_native_ligand.pdb redock_1XOQ_wet.pdbqt

# benchmark — 146 actives, 2,246 matched decoys
python select_more_actives.py chembl_raw.pkl actives_chembl.smi 110 0
python make_decoys_150.py
python prep_150.py actives_150.smi decoys_150.smi
bash run_150.sh                      # ~31 h, 8 cores, restartable

# gate 2 — can anything rank these compounds
python fingerprint.py 1XOQ_wet.pdb 1XOQ_wet_native_ligand.pdb \
       docked150 actives150 fingerprints_150.csv
python fp_classify.py fingerprints_150.csv actives_150.smi decoys_150.smi
```

`make_decoys_150.py`, `prep_150.py` and `run_150.sh` are the Phase 1 scripts
with three settings changed: 15 decoys per active instead of 30, the 150-active
input, and the `actives150` / `decoys150` / `docked150` folders.

---

## Data files to add alongside these scripts

Not included here — add from the working folder:

```
actives_150.smi            actives_extra.csv          decoys_150.smi
decoy_report_150.csv       fingerprints_150.csv       fingerprints.csv
rescore_results.csv        1XOQ_wet.pdb               1XOQ_wet.pdbqt
1XOQ_wet_native_ligand.pdb redock_1XOQ_wet.pdbqt      config_1XOQ.txt
```

The `docked150/` folder is 2,392 poses plus logs and is too large for git —
it regenerates from `run_150.sh`.

---

## Open items

- The `name` column in `actives_chembl.csv` needs checking: the control
  independence check assumed it holds compound names rather than ChEMBL IDs
- `calc_rmsd.py` reads only the first MODEL of a multi-pose file. Harmless
  for the top-ranked pose, which is what Gate 1 uses, but worth fixing
- `metal_site.py shell` lists chain A and chain B metals together, and its
  bridging-water check requires one water near *every* metal found, so it
  misses the real bridge. In 1XOQ, HOH1010 bridges Zn1001 and Mg1002 at
  2.18 Å and 2.07 Å
- The metal-coordination hypothesis for the scoring failure was **tested and
  rejected**: rank vs distance to metal gave rho = +0.326 (p = 0.040), the
  opposite direction. Do not state it as the explanation
- 3 of 149 prepared actives and 4 of 2,250 decoys failed 3D generation, and
  one active failed pdbqt conversion (CHEMBL593656). Final counts are 146
  and 2,246

---

## Not established

- No antidiabetic has been screened yet
- No potency or affinity prediction is possible from this work. The docking
  score does not rank for this target, so scores must never be reported as
  predicted potency
- AUC 0.822 measures separation of known actives from property-matched
  decoys. It does not establish that a high-ranking compound inhibits PDE4D
