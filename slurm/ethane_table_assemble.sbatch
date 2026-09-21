#!/bin/bash
# Runs after the array (afterany).  Normally every column is cached and this just
# writes the HDF5 table; any column an array task failed on is retried here,
# serially.  Its output is the build log.
#SBATCH --job-name=eth_table_asm
#SBATCH --partition=main
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2000M
#SBATCH --time=06:00:00
# If every column is cached this takes minutes.  If it runs out of time while
# retrying failed columns, just resubmit: each finished column is already cached.
#SBATCH --chdir=/models/hextor
#SBATCH --output=/models/hextor/model/radiation/slurm_logs/eth_assemble_%j.out
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
exec python3 -u tools/make_radiation_table.py "${TABLE_ARGS[@]}" \
     --workers 1 --timeout "$COLUMN_TIMEOUT"
