GATE 2 - the enrichment test
============================

Your teammates' Phase 1 files, with one fix, plus the scripts to run Gate 2.

WHAT I CHANGED
  decoys.sdf   one decoy (decoy_C94816686) had a torn 3D structure - a 4.71 A
               "single bond". Rebuilt with RDKit ETKDGv3 + MMFF94. The other
               349 were fine. Everything else is exactly as they sent it.

THE ONE DECISION THEY FLAGGED - and the answer
  They noticed the controls are neutral but DUD-E matched some decoys to
  charged forms. Use the SEVEN PROTOMERS as the actives for this test, not
  the six neutral controls. DUD-E matched each decoy's charge to a specific
  protonation state, so pairing neutral cilomilast with decoys matched to the
  cilomilast anion breaks the very property matching the decoys exist to
  provide. Using the protomers makes it consistent and gives exactly 50
  decoys per active. prep_gate2.py already does this.

  Chemically this is also the right call: cilomilast's carboxylic acid is
  ionised at pH 7.4, so the anion is the real species.

HOW TO RUN
  Copy these files into ~/phppde4 (where 1XOQ.pdbqt and config_1XOQ.txt are).

    conda activate phppde4
    python prep_gate2.py        # makes actives/ and decoys/ as pdbqt
    bash run_gate2.sh           # docks all 357 - a few hours, restartable
    python enrichment.py        # the verdict

  run_gate2.sh skips anything already docked, so you can stop it with Ctrl+C
  and start it again later.

WHAT THE NUMBERS MEAN
  AUC-ROC    whole-ranking quality. 0.5 is random, 0.7 acceptable, 0.8 good.
  EF at 1%   how many times better than chance the very top of the list is.
             This is what a screen actually uses.
  BEDROC     like EF but smoother, weights the top heavily. Standard to report.

  Report all three. AUC alone can look respectable while the top of the list
  is useless, which is the part that matters.

IF IT FAILS
  Do not screen anyway. Check the receptor still has its zinc and magnesium,
  that the box is centred on the pocket, and the protonation states. A failed
  Gate 2 reported honestly is a finding; a screen run on a failed protocol is
  three months wasted.
