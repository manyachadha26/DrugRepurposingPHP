# Gate 1: redocking validation

Date: 10 October 2026
AutoDock Vina 1.2.7 (macOS arm64). Receptors prepared with Meeko 0.8.0 from
PDB files cleaned by `clean_structures.py`. 7CBQ side chains repaired with
PDBFixer. Grid boxes centred on each co-crystallised ligand's centroid.
RMSD computed with RDKit `CalcRMS` (symmetry-aware, no superposition).

| Structure | Target | Ligand | Top-ranked pose RMSD | Best pose RMSD (rank) | Top score | Outcome |
|---|---|---|---|---|---|---|
| 1XOQ | PDE4D | roflumilast | **0.770 Å** | 0.770 Å (1) | −9.488 | Pass |
| 1XMU | PDE4B | roflumilast | **0.942 Å** | 0.942 Å (1) | −9.597 | Pass |
| 7CBQ | PDE4D | apremilast | 3.415 Å | 0.694 Å (3) | −9.173 | Sampling pass, scoring fail |

## Interpretation

1XOQ and 1XMU reproduce the crystallographic binding mode, and in both the
top-scored pose is the correct one. Those two are used for scoring.

For 7CBQ the correct pose is sampled (0.694 Å) but ranked third. Repeating at
exhaustiveness 64 gave 3.420 Å for the top pose and 0.694 Å at rank 3 again,
with scores unchanged to two decimal places, so the search is converged and the
discrepancy is a scoring-function limitation rather than insufficient sampling.
Apremilast is larger and has more rotatable bonds than roflumilast, which makes
the ranking harder. 7CBQ is therefore used for pose comparison only.

## Benchmark for the screen

Roflumilast in PDE4D (1XOQ): **−9.488 kcal/mol**. Candidate antidiabetic
compounds are compared against this value.

## Methods sentence

Redocking of roflumilast into PDE4D (PDB 1XOQ) reproduced the crystallographic
binding mode with a heavy-atom RMSD of 0.77 Å for the top-ranked pose, and
0.94 Å for PDE4B (1XMU). For PDE4D structure 7CBQ the crystallographic pose was
recovered among the sampled conformations (0.69 Å) but ranked third; this
structure was consequently used for pose comparison rather than scoring.

## Settings

exhaustiveness 32 (64 tested for 7CBQ), num_modes 9, seed 42, grid spacing 0.375 Å.
Box dimensions per structure are in the `config_*.txt` files.
