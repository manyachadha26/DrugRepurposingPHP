#!/bin/bash
MODE="$1"; f="$2"
base=$(basename "$f" .pdbqt)
o="taskC/docked_3HR1_${MODE}/${base}.pdbqt"
[ -s "$o" ] && exit 0
./vina --receptor "3HR1_${MODE}.pdbqt" --ligand "$f" --config config_3HR1.txt \
       --cpu 1 --out "$o" > /dev/null 2>&1
