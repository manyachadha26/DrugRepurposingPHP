#!/bin/bash
# Gate 2 round 2, run in parallel.
#
# Vina scales poorly across cores, so one job using 8 cores is much slower
# than 8 jobs using 1 core each. This usually cuts total time 3-4x.
#
# Restartable: skips anything already finished.
# Stop with Ctrl+C; rerun to continue.
set -u

RECEPTOR=1XOQ.pdbqt
CONFIG=config_1XOQ.txt
VINA=./vina
JOBS=${1:-8}          # bash run_gate2_parallel.sh 6   to use 6 at a time

mkdir -p docked2
total=$(ls actives2/*.pdbqt decoys2/*.pdbqt 2>/dev/null | wc -l | tr -d ' ')
done_already=$(ls docked2/*.pdbqt 2>/dev/null | wc -l | tr -d ' ')
echo "$total to dock, $done_already already finished, running $JOBS at a time"
echo "started $(date '+%H:%M:%S')"

dock_one() {
    f="$1"
    base=$(basename "$f" .pdbqt)
    out="docked2/${base}.pdbqt"
    [ -s "$out" ] && return 0
    ./vina --receptor "$RECEPTOR_X" --ligand "$f" --config "$CONFIG_X" \
           --out "$out" --num_modes 9 --cpu 1 > "docked2/${base}.log" 2>&1
}
export -f dock_one
export RECEPTOR_X="$RECEPTOR" CONFIG_X="$CONFIG"

ls actives2/*.pdbqt decoys2/*.pdbqt 2>/dev/null \
  | xargs -P "$JOBS" -I {} bash -c 'dock_one "$@"' _ {} &

PID=$!
while kill -0 $PID 2>/dev/null; do
    sleep 60
    n=$(ls docked2/*.pdbqt 2>/dev/null | wc -l | tr -d ' ')
    echo "  $(date '+%H:%M:%S')  $n / $total"
done
wait $PID
echo "Finished $(date '+%H:%M:%S')"
echo "Now: python enrichment.py actives2 docked2"
