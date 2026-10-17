#!/bin/bash
# Gate 2 round 2 - dock the ChEMBL actives and their decoys.
# Restartable: skips anything already in docked_screen/.
set -u
RECEPTOR=1XOQ_wet.pdbqt
CONFIG=config_1XOQ.txt
VINA=./vina

mkdir -p docked_screen
total=$(ls screen_lig/*.pdbqt  2>/dev/null | wc -l | tr -d ' ')
n=0
start=$(date +%s)
for f in screen_lig/*.pdbqt ; do
    n=$((n+1))
    base=$(basename "$f" .pdbqt)
    out="docked_screen/${base}.pdbqt"
    [ -s "$out" ] && continue
    $VINA --receptor "$RECEPTOR" --ligand "$f" --config "$CONFIG" \
          --out "$out" --num_modes 9 > "docked_screen/${base}.log" 2>&1
    if [ $((n % 25)) -eq 0 ]; then
        el=$(( $(date +%s) - start ))
        echo "[$n/$total]  ${el}s elapsed"
    fi
done
echo "Done. Now: python enrichment.py screen_lig docked_screen"
