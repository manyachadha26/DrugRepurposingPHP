#!/bin/bash
# Gate 2 round 2 - dock the ChEMBL actives and their decoys.
# Restartable: skips anything already in docked150/.
set -u
RECEPTOR=1XOQ_wet.pdbqt
CONFIG=config_1XOQ.txt
VINA=./vina

mkdir -p docked150
total=$(ls actives150/*.pdbqt decoys150/*.pdbqt 2>/dev/null | wc -l | tr -d ' ')
n=0
start=$(date +%s)
for f in actives150/*.pdbqt decoys150/*.pdbqt; do
    n=$((n+1))
    base=$(basename "$f" .pdbqt)
    out="docked150/${base}.pdbqt"
    [ -s "$out" ] && continue
    $VINA --receptor "$RECEPTOR" --ligand "$f" --config "$CONFIG" \
          --out "$out" --num_modes 9 > "docked150/${base}.log" 2>&1
    if [ $((n % 25)) -eq 0 ]; then
        el=$(( $(date +%s) - start ))
        echo "[$n/$total]  ${el}s elapsed"
    fi
done
echo "Done. Now: python enrichment.py actives150 docked150"
