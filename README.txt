STEP 4 OUTPUT - receptor preparation
Produced 10 October 2026 by clean_structures.py + Meeko 0.8.0

READY TO DOCK
  1XOQ.pdbqt                 primary target (PDE4D)        grid: config_1XOQ.txt
  7CBQ.pdbqt                 second PDE4D                  grid: config_7CBQ.txt

STILL NEEDS ONE COMMAND ON YOUR MACHINE
  1XMU.pdbqt                 not produced here - see below

CLEANED RECEPTORS (input to the pdbqt step)
  1XOQ_clean.pdb  7CBQ_clean.pdb  1XMU_clean.pdb
  Chain A only. Waters and ethanediol removed. Zinc and magnesium KEPT.
  Alternate conformations reduced to A. Original residue numbering preserved,
  so Gln369 / Gln443 are still Gln369 / Gln443.

NATIVE LIGANDS (for the Week 4 redocking test)
  1XOQ_native_ligand.pdb/.pdbqt   roflumilast  (26 atoms)
  7CBQ_native_ligand.pdb/.pdbqt   apremilast   (32 atoms)
  1XMU_native_ligand.pdb/.pdbqt   roflumilast  (26 atoms)

GRID BOXES - sized to each ligand plus 8 A clearance, not guessed
  config_1XOQ.txt   centre 14.117, 28.703, 53.487    box 22 x 26 x 28
  config_7CBQ.txt   centre 84.199, 157.634, 27.011   box 26 x 24 x 26
  config_1XMU.txt   centre 2.335, 0.725, 41.711      box 28 x 22 x 28

TWO THINGS TO FINISH
1. 1XMU: its residue 194 is CME, a chemically modified cysteine. Meeko has to
   fetch a definition file for it from rcsb.org, and this sandbox has no
   internet. On your machine, with internet, run:
       mk_prepare_receptor.py -i 1XMU_clean.pdb -o 1XMU -p

2. 7CBQ: eight residues have side chains the crystallographers could not
   resolve, so they were dropped (-a flag). Seven sit 25-36 A from the pocket
   and do not matter. ARG342 is 12.9 A from the centre, which is inside the
   search box, so for the final run rebuild it instead of dropping it:
       conda install -c conda-forge pdbfixer openmm
       pdbfixer 7CBQ_clean.pdb --add-atoms=heavy --output=7CBQ_fixed.pdb
       mk_prepare_receptor.py -i 7CBQ_fixed.pdb -o 7CBQ -p
   The current 7CBQ.pdbqt is fine for testing the pipeline meanwhile.

NEXT: Week 4 Gate 1, the redocking test.
       vina --receptor 1XOQ.pdbqt --ligand 1XOQ_native_ligand.pdbqt \
            --config config_1XOQ.txt --out redock_1XOQ.pdbqt
   Then measure RMSD against the native ligand. Under 2 A and you proceed.
