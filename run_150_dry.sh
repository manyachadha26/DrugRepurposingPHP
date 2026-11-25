#!/bin/bash
set -u
mkdir -p docked150_dry
ls actives150/*.pdbqt decoys150/*.pdbqt | xargs -P 8 -n 1 ./dock_one_dry.sh
