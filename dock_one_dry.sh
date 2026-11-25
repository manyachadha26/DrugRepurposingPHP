#!/bin/bash
f="$1"
base=$(basename "$f" .pdbqt)
out="docked150_dry/${base}.pdbqt"
[ -s "$out" ] && exit 0
./vina --receptor 1XOQ.pdbqt --ligand "$f" --config config_box.txt \
       --cpu 1 --out "$out" > "docked150_dry/${base}.log" 2>&1
