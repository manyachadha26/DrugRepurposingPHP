# Phase 3 — the screen, and whether its hits survive

Screening 34 antidiabetic drugs plus 6 known PDE4 inhibitors against PDE4D,
using the interaction-pattern method validated in Phase 2 (AUC 0.822).

**Result: one candidate, sitagliptin.**

---

## The screen

41 compounds: 34 antidiabetics, 6 known PDE4 inhibitors as positive
controls, plus cilomilast re-docked as its anion. Docked against the wet
1XOQ receptor, ranked by the Phase 2 classifier trained on 146 actives and
2,246 matched decoys.

### The controls validate it

| Control | Rank of 41 | Above this % of known actives |
|---|---|---|
| zardaverine | 1 | 92% |
| roflumilast | 2 | 90% |
| apremilast | 3 | 88% |
| rolipram | 4 | 87% |
| theophylline | 9 | 40% |
| cilomilast | 10 | 40% |

All 6 in the top 10 (p = 5.5×10⁻⁵); 4 of 6 took ranks 1–4 (p = 1.6×10⁻⁴).
Controls versus library, Mann-Whitney p = 1.7×10⁻⁵.

**Theophylline ranking 9th is correct, not a failure.** Its PDE4 IC50 is
around 100 µM, roughly a thousandfold weaker than the pIC50 7–10 compounds
the model was trained on. A method that placed it beside roflumilast would
be the broken one.

**Cilomilast at 40% is a genuine false negative.** It is a potent inhibitor
(~100 nM) that the method ranked mid-pack. Re-docking it as the anion made
it *worse* (rank 13, 27%), so protonation is not the explanation — the
hypothesis was tested and rejected. Its pose is fine; the scoring is what
missed it. The model leans heavily on burial, and cilomilast is smaller and
less enclosed than the compounds it was trained on. That is the method's
known weakness and belongs in the limitations.

### The two initial hits

| | Rank | % of actives below | d(Gln369) | Clamp | d(metal) |
|---|---|---|---|---|---|
| sotagliflozin | 5 | 71% | 3.24 Å | 2 | 3.58 Å |
| sitagliptin | 6 | 65% | 3.13 Å | 2 | 4.00 Å |

---

## Visual inspection against the crystal reference

Poses compared against crystallographic roflumilast in PyMOL
(`inspect_hits.pml`).

| | Atoms within 2 Å of roflumilast | Contacts < 2.5 Å from metals/waters |
|---|---|---|
| sotagliflozin | 25 of 31 — 81% | 1 |
| sitagliptin | 24 of 30 — 80% | 1 |
| cilomilast *(known inhibitor)* | 22 of 26 — 85% | — |

Both hits occupy essentially the same pocket volume as a marketed PDE4
inhibitor. Cilomilast has the best overlap of the three and was still
mis-ranked, which confirms the weakness is in scoring, not in pose
generation.

---

## Cross-structure robustness — the test that changed the answer

Every result above rests on one receptor conformation. The 41 compounds were
re-docked against a second PDE4D structure.

### Receptor preparation and gating

Both candidate structures were re-prepared keeping the metal-coordinating
waters. 7CBQ additionally needed 8 incomplete side chains rebuilt
(`fix_sidechains.py`) — seven sat 21–31 Å from the ligand, but Arg342 at
9.0 Å was inside the grid box, so deletion was not acceptable.

| Structure | Redocking RMSD (wet) | Redocking RMSD (dry, Phase 1) | Used |
|---|---|---|---|
| 1XOQ | 0.756 Å | 0.770 Å | yes |
| 1XMU | 0.927 Å | — | yes |
| 7CBQ | 3.511 Å | 3.415 Å | **excluded** |

7CBQ failed a gate set in advance, on both preparations. The water fix did
not rescue it, which is itself informative: that fix corrected the
metal-proximity artifact specifically and did not paper over an unrelated
problem.

### Why raw features, not the trained model

The classifier's per-residue features are named by residue number, and the
two structures number differently: 1XOQ's Gln369, Phe340, Phe372 are
Gln443, Phe414, Phe446 in 1XMU — a constant offset of 74. The script found
the equivalent residues **geometrically in each structure**, with nothing
hardcoded, which is also an independent check that the detection works.

### The consensus criterion, fixed before looking

A compound is consistent if, in **both** structures, it hydrogen-bonds the
invariant glutamine (≤ 3.5 Å) **and** contacts at least one clamp
phenylalanine.

| | 1XOQ | 1XMU | Both |
|---|---|---|---|
| Compounds satisfying the pharmacophore | 18 | 15 | 14 |

The two structures agree on 36 of 41 compounds (88%). Five of six controls
are consistent; only zardaverine flips, at 3.49 Å versus 3.68 Å — a
sharp-threshold artifact, not a difference in binding.

### Outcome

| | 1XOQ | 1XMU | Verdict |
|---|---|---|---|
| **sotagliflozin** | 3.24 Å, H-bond | **4.47 Å, no H-bond** | conformation-dependent — **dropped** |
| **sitagliptin** | 3.13 Å, H-bond | 3.03 Å, H-bond | **consistent — retained** |

Sotagliflozin loses the glutamine hydrogen bond entirely in the second
structure, a full Ångström past the cutoff. Without this test it would have
been reported as a hit.

---

## A ligand-based model, and why it could not corroborate

A random forest on ECFP4 fingerprints plus 16 descriptors, trained on the
same 146 actives and 2,246 decoys, scaffold-grouped 10-fold
cross-validation.

Out-of-fold AUC 0.988, permutation control 0.507.

**That 0.988 is inflated by the benchmark design and must not be reported
as better than the structure-based 0.822.** The decoys were selected to have
Tanimoto below 0.35 to every active, so a fingerprint model separates them
almost perfectly by rediscovering the selection rule. The structure-based
model never saw 2D structure and was unaffected.

It ranked all six controls 1–6, then placed every antidiabetic at the
bottom: sitagliptin 14th, sotagliflozin 20th, all at 0%.

This is **not** evidence against the hits. A ligand-based model cannot score
a novel chemotype highly — it recognises only what it has seen, and with
decoys chosen for dissimilarity it scores every unfamiliar scaffold near
zero. It is structurally incapable of providing evidence either way on a
scaffold hop.

**Consequence, stated plainly: the sitagliptin result has no independent
corroboration.** It rests on the structure-based evidence alone.

The contrast is worth reporting: ligand-based screening put all six known
PDE4 drugs in the top six and every novel-chemotype candidate at the bottom,
while structure-based screening found candidates in unrelated drug classes.
That is a concrete demonstration of why repurposing across chemotypes needs
structure.

---

## What must be said carefully

**14 of 41 compounds satisfy the pharmacophore in both structures** —
including metformin, rosiglitazone, acetohexamide and lobeglitazone.
Pharmacophore consensus is a necessary condition, not a sufficient one. What
separates sitagliptin (rank 6, 65%) from metformin (rank 21, 12%) is the
ranking, not the pharmacophore. The 14 are not 14 hits.

**No potency is predicted.** The docking score does not rank for this target.
Nothing here estimates affinity.

**One candidate, for experimental testing.** Sitagliptin reproduces the
canonical PDE4 binding mode in two independent crystal conformations and
ranks above two marketed PDE4 inhibitors. That is a hypothesis for an enzyme
assay, not a finding about inhibition.

---

## Scripts

| Script | Purpose |
|---|---|
| `fix_sidechains.py` | Rebuild incomplete side chains with PDBFixer, verifying metals and waters survive |
| `compare_structures.py` | Cross-structure pharmacophore consensus |
| `ligand_model.py` | Ligand-based random forest, 2D structure only |
| `inspect_hits.pml` | PyMOL session for visual inspection against the crystal reference |

`screen_library.py` and `pose_filter.py` are in `phase2/`.

---

## Reproducing Phase 3

```bash
conda activate phppde4

# the screen
cat library.smi controls.smi > screen_set.smi
python prep_screen.py screen_set.smi empty.smi
bash run_screen.sh
python fingerprint.py 1XOQ_wet.pdb 1XOQ_wet_native_ligand.pdb \
       docked_screen screen_lig library_fp.csv
python screen_library.py fingerprints_150.csv library_fp.csv controls.txt

# second structure
python clean_with_waters.py raw/1XMU.pdb
mk_prepare_receptor.py -i 1XMU_wet.pdb -o 1XMU_wet -p
python patch_waters_pdbqt.py 1XMU_wet.pdb 1XMU_wet.pdbqt
vina --receptor 1XMU_wet.pdbqt --ligand 1XMU_native_ligand.pdbqt \
     --config config_1XMU.txt --out redock_1XMU_wet.pdbqt
python calc_rmsd.py 1XMU_wet_native_ligand.pdb redock_1XMU_wet.pdbqt
bash run_screen_1XMU.sh
python fingerprint.py 1XMU_wet.pdb 1XMU_wet_native_ligand.pdb \
       docked_screen_1XMU screen_lig library_fp_1XMU.csv
python compare_structures.py screen_results.csv library_fp.csv \
       library_fp_1XMU.csv

# ligand-based model
python ligand_model.py actives_150.smi decoys_150.smi screen_all.smi
```

---

## Still to do

- Cross-validate the PDE4D expression finding against a second human kidney
  dataset (KPMP or a published single-nucleus atlas). It currently rests on
  HPA alone
- Check against pharmacokinetic sources whether sitagliptin reaches the
  proximal tubule at relevant concentrations. Not yet verified
- MD and MM-GBSA on the sitagliptin complex — needs a GPU, declared as
  future work
- Author order, agreed in writing
- The eight Vancouver references in the introduction need checking against
  the actual papers; two citations are still missing (PHP prevalence, and a
  review of PDE4 inhibitor clinical development covering the emesis dose
  limit)
