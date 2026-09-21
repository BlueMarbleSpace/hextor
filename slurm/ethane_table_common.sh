# Shared settings for the v5 ethane radiation-table build under Slurm.
# Sourced by submit_ethane_table.sh and both sbatch scripts, so the manifest,
# the array tasks and the assembly job can never disagree about the grid.
# These match run_ethane_table.sh (the pre-Slurm nohup launcher) exactly, so
# the cache it built is reused as-is.
HEXTOR=/models/hextor
OUT=$HEXTOR/model/radiation/radiation_N2_CO2_CH4_C2H6_Sun_p_rh0.8.h5
CACHE=$OUT.cache
LOGDIR=$HEXTOR/model/radiation/slurm_logs
TABLE_ARGS=(--out "$OUT" --cache "$CACHE" --star G2V_SUN_n68.nc
            --ch4-axis --rh 0.8 --c2h6-photochem)
# Per-column ExoColumn timeout.  The script default (600 s) is sized for a
# quiet node; a full column is 2700 records, and while this shares the node
# with other ExoRT work it measured ~4.5 records/s, i.e. right at 600 s.  The
# Slurm --time limit is the real guard against a hung process.  1800 s is 3x the
# contended estimate; a column that still times out is retried by the assembly job.
COLUMN_TIMEOUT=1800
