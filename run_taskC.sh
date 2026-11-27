#!/bin/bash
set -u
for MODE in wet dry; do
  mkdir -p taskC/docked_3HR1_$MODE
  ls taskC/pde10a_lig/*.pdbqt | xargs -P 8 -n 1 ./dock_taskC_one.sh $MODE
  echo "$MODE: $(ls taskC/docked_3HR1_$MODE | grep -c pdbqt)"
done
