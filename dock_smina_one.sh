#!/bin/bash
f="$1"
base=$(basename "$f" .pdbqt)
o="docked_screen_smina/${base}.pdbqt"
[ -s "$o" ] && exit 0
smina --receptor 1XOQ_wet.pdbqt --ligand "$f" --config config_box.txt \
      --cpu 1 --out "$o" > /dev/null 2>&1
