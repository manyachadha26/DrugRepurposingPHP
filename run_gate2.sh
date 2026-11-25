#!/bin/bash
# Gate 2, part 2: dock every active and decoy against the primary target.
#
# 357 molecules at exhaustiveness 32. Expect roughly 2-4 hours on an M4.
# Safe to stop and restart - it skips anything already done.

set -u
RECEPTOR=1XOQ.pdbqt
CONFIG=config_1XOQ.txt
VINA=./vina            # change to "vina" if you installed it through conda

mkdir -p docked
total=$(ls actives/*.pdbqt decoys/*.pdbqt 2>/dev/null | wc -l | tr -d ' ')
n=0
for f in actives/*.pdbqt decoys/*.pdbqt; do
    n=$((n+1))
    base=$(basename "$f" .pdbqt)
    out="docked/${base}.pdbqt"
    if [ -s "$out" ]; then continue; fi
    $VINA --receptor "$RECEPTOR" --ligand "$f" --config "$CONFIG" \
          --out "$out" --num_modes 9 > "docked/${base}.log" 2>&1
    echo "[$n/$total] $base"
done
echo "Docking finished. Now: python enrichment.py"
