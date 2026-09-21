#!/bin/bash
# One array task = PER_TASK consecutive lines of MANIFEST (flat column indices),
# computed serially into the shared cache.  Writes no table; re-running is safe.
#SBATCH --job-name=eth_table
#SBATCH --partition=main
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1000M
#SBATCH --time=01:10:00
# 2 columns x COLUMN_TIMEOUT (1800 s) + overhead.  Typical is 10-20 min; keeping
# the limit realistic lets a task use merlin's 4 night-only CPUs, which Slurm
# grants only to jobs that will end before the 07:00 daytime reservation.
#SBATCH --chdir=/models/hextor
#SBATCH --output=/models/hextor/model/radiation/slurm_logs/eth_%A_%a.out
set -eo pipefail
source /models/hextor/slurm/ethane_table_common.sh
# --force: sbatch hands the job the SUBMITTING shell's environment, and if that
# shell had already sourced setvars.sh it carries SETVARS_COMPLETED=1, which makes
# setvars.sh return 3 without doing anything.  Under `set -e` with its output
# discarded, that killed every task of array 55 silently (exit 3, empty log).
# Check it explicitly so any future failure here says why.
if ! source /opt/intel/oneapi/setvars.sh --force >/dev/null 2>&1; then
    echo "ERROR: /opt/intel/oneapi/setvars.sh failed; ExoColumn needs its libimf" >&2
    exit 1
fi
: "${MANIFEST:?MANIFEST not set}" "${PER_TASK:?PER_TASK not set}"
first=$(( SLURM_ARRAY_TASK_ID * PER_TASK + 1 ))
last=$(( first + PER_TASK - 1 ))
cols=$(sed -n "${first},${last}p" "$MANIFEST" | paste -sd, -)
if [ -z "$cols" ]; then echo "task $SLURM_ARRAY_TASK_ID: no columns"; exit 0; fi
echo "task $SLURM_ARRAY_TASK_ID on $(hostname): columns $cols"
exec python3 -u tools/make_radiation_table.py "${TABLE_ARGS[@]}" \
     --only-columns "$cols" --timeout "$COLUMN_TIMEOUT"
