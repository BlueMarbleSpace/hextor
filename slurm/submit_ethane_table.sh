#!/bin/bash
# Submit (or resubmit) the v5 ethane table build to Slurm at the lowest
# runnable priority.  Resume-safe: the manifest lists only columns not yet in
# the cache, so running this again after a partial build submits just the rest.
# Do NOT re-run while a previous array from this script is still queued or
# running -- the two would duplicate work (harmless, but wasteful).
#
#   NICE     default 9999.  Priority here is 10000 x partition factor - nice,
#            and priority 0 means HELD, so 9999 -> priority 1 is the floor.
#   MAXCONC  default 4.  ExoColumn is memory-bandwidth bound (~240 MB/process);
#            throughput on this 6-core/12-thread box peaks at 3-4 concurrent and
#            falls beyond, so running more tasks at once finishes no sooner.
set -eo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
source "$HERE/ethane_table_common.sh"
# Deliberately NOT sourcing setvars.sh here: the manifest step only parses text
# files, and whatever this shell exports is inherited by every job it submits.
NICE=${NICE:-9999}
MAXCONC=${MAXCONC:-4}
MAXTASKS=1000          # MaxArraySize is 1001 here, so indices must stay <= 1000
mkdir -p "$LOGDIR"

# The array tasks read the manifest at run time, so each submission gets its
# own file rather than overwriting one a queued array still depends on.
MANIFEST=$LOGDIR/manifest_$(date +%Y%m%d_%H%M%S).txt
python3 "$HEXTOR/tools/make_radiation_table.py" "${TABLE_ARGS[@]}" \
        --list-uncached "$MANIFEST"
N=$(wc -l < "$MANIFEST")

asm_dep=()
if [ "$N" -gt 0 ]; then
  PER_TASK=$(( (N + MAXTASKS - 1) / MAXTASKS ))
  NTASKS=$(( (N + PER_TASK - 1) / PER_TASK ))
  ARR=$(sbatch --parsable --nice="$NICE" \
          --array="0-$((NTASKS - 1))%${MAXCONC}" \
          --export=ALL,MANIFEST="$MANIFEST",PER_TASK="$PER_TASK" \
          "$HERE/ethane_table_array.sbatch")
  echo "array job $ARR: $N columns as $NTASKS tasks x $PER_TASK, <= $MAXCONC at once, nice $NICE"
  asm_dep=(--dependency="afterany:$ARR")
else
  echo "every column is already cached; submitting assembly only"
fi
ASM=$(sbatch --parsable --nice="$NICE" "${asm_dep[@]}" "$HERE/ethane_table_assemble.sbatch")
echo "assembly job $ASM${asm_dep:+ (after $ARR)}"
echo "manifest: $MANIFEST"
