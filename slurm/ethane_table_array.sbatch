#!/bin/bash
# One array task = PER_TASK consecutive lines of MANIFEST (flat column indices),
# computed serially into the shared cache.  Writes no table; re-running is safe.
#SBATCH --job-name=eth_table
#SBATCH --partition=main
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1000M
#SBATCH --time=03:00:00
#SBATCH --chdir=/models/hextor
#SBATCH --output=/models/hextor/model/radiation/slurm_logs/eth_%A_%a.out
set -eo pipefail
source /models/hextor/slurm/ethane_table_common.sh
source /opt/intel/oneapi/setvars.sh >/dev/null 2>&1   # ExoColumn is ifx-built (libimf)
: "${MANIFEST:?MANIFEST not set}" "${PER_TASK:?PER_TASK not set}"
first=$(( SLURM_ARRAY_TASK_ID * PER_TASK + 1 ))
last=$(( first + PER_TASK - 1 ))
cols=$(sed -n "${first},${last}p" "$MANIFEST" | paste -sd, -)
if [ -z "$cols" ]; then echo "task $SLURM_ARRAY_TASK_ID: no columns"; exit 0; fi
echo "task $SLURM_ARRAY_TASK_ID on $(hostname): columns $cols"
exec python3 -u tools/make_radiation_table.py "${TABLE_ARGS[@]}" \
     --only-columns "$cols" --timeout "$COLUMN_TIMEOUT"
