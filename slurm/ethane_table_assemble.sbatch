#!/bin/bash
# Runs after the array (afterany).  Normally every column is cached and this just
# writes the HDF5 table; any column an array task failed on is retried here,
# serially.  Its output is the build log.
#SBATCH --job-name=eth_table_asm
#SBATCH --partition=main
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2000M
#SBATCH --time=1-00:00:00
#SBATCH --chdir=/models/hextor
#SBATCH --output=/models/hextor/model/radiation/slurm_logs/eth_assemble_%j.out
set -eo pipefail
source /models/hextor/slurm/ethane_table_common.sh
source /opt/intel/oneapi/setvars.sh >/dev/null 2>&1
exec python3 -u tools/make_radiation_table.py "${TABLE_ARGS[@]}" \
     --workers 1 --timeout "$COLUMN_TIMEOUT"
