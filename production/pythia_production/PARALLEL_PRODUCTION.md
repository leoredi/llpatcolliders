# Parallel Pythia Production Guide

## Overview

The Pythia production stage generates ~150 HNL samples (different masses and flavours) which takes **~20 hours on a single core**. The parallel production script reduces this to **~2-3 hours** on an 8-core machine.

## Quick Start

### 1. Use Bash (fish/zsh users)

The script uses Bash-only features (`read -p`, arrays). If your login shell is fish or zsh, launch it via Bash:

```bash
cd production/pythia_production
bash -lc 'source activate llpatcolliders && ./run_parallel_production.sh muon'
```

(Or open a Bash shell, `conda activate llpatcolliders`, then run the script.)

### 2. Configure Parallelism

Edit `run_parallel_production.sh` line 15:

```bash
MAX_PARALLEL=8  # Adjust based on your CPU cores
```

**Recommended settings:**
- **4 cores:** `MAX_PARALLEL=3` (leave 1 for system)
- **8 cores:** `MAX_PARALLEL=7`
- **16 cores:** `MAX_PARALLEL=14`
- **Mac M1/M2:** Check `sysctl hw.ncpu` and use N-2

### 3. Run Full Production

```bash
cd production/pythia_production
./run_parallel_production.sh
```

**What it does:**
- Runs the meson production grid (electron, muon, tau)
- Runs jobs in parallel (respects `MAX_PARALLEL` limit)
- Each job: 100k events → ~10 minutes
- Total time: ~2-3 hours (vs ~20 hours sequential)

**Output:**
- CSV files: `../../output/csv/simulation/HNL_*.csv` (from `production/pythia_production`)
- Logs: `../../output/logs/simulation/HNL_*.log`
- Summary log: `../../output/logs/simulation/production_run_TIMESTAMP.log`

### 4. Monitor Progress

In a separate terminal:

```bash
# Watch job count
watch -n 5 "jobs -r | wc -l"

# Watch CSV files being created
watch -n 10 "ls ../../output/csv/simulation/*.csv | wc -l"

# Check recent completions
tail -f ../../output/logs/simulation/production_run_*.log
```

### 5. Verify Results

After completion:

```bash
# Count output files (should be ~150)
ls ../../output/csv/simulation/HNL_*.csv | wc -l

# Check for failures
grep -l "FAILED\|ERROR" ../../output/logs/simulation/*.log

# Check file sizes (should all be > 1 MB)
find ../../output/csv/simulation -name "*.csv" -size -1M

# Check for overlaps with MadGraph EW files
cd ../../analysis_pbc
python limits/combine_production_channels.py --dry-run
```

## Script Details

### run_parallel_production.sh

**Key features:**
- Job queue with `MAX_PARALLEL` limit
- Background job execution with `&`
- Automatic log file management
- Progress tracking
- Summary statistics

**Safety features:**
- Waits for available slot before starting new job
- Each job has isolated log file (no conflicts)
- Failed jobs logged separately
- Estimated time calculation

### Job Control Functions

```bash
# Count running jobs
count_jobs() {
    jobs -r | wc -l | tr -d ' '
}

# Wait for free slot
wait_for_slot() {
    while [ $(count_jobs) -ge $MAX_PARALLEL ]; do
        sleep 1
    done
}
```

### Log File Naming

Each job creates a unique log:
```
HNL_2p6GeV_muon_direct_20251202_143022.log
    ^mass   ^flav  ^mode   ^timestamp
```

## Performance Estimates

| CPU Cores | MAX_PARALLEL | Time (est.) | Speedup |
|-----------|--------------|-------------|---------|
| 1 (single)| 1           | ~20 hours   | 1×      |
| 4 cores   | 3           | ~7 hours    | 3×      |
| 8 cores   | 7           | ~3 hours    | 7×      |
| 16 cores  | 14          | ~1.5 hours  | 13×     |

**Notes:**
- Time per job varies: kaon (5 min) < charm (8 min) < beauty (12 min)
- I/O can bottleneck at very high parallelism
- Leave 1-2 cores for system responsiveness

## Common Issues

### Issue: Jobs fail with "library not found"

**Solution:** Check Pythia library path
```bash
export DYLD_LIBRARY_PATH="/path/to/pythia/lib:$DYLD_LIBRARY_PATH"
./run_parallel_production.sh
```

### Issue: Too many jobs slow down system

**Solution:** Reduce `MAX_PARALLEL`
```bash
# Edit line 15 in run_parallel_production.sh
MAX_PARALLEL=4  # Lower value
```

### Issue: Disk full during production

**Check available space:**
```bash
df -h ../../output/csv/simulation
```

**Each CSV is ~5-50 MB**; the full meson set is ~5 GB. If you also add EW files later, plan for more headroom.

### Issue: Jobs hang or freeze

**Kill all background jobs:**
```bash
killall main_hnl_production
jobs -p | xargs kill
```

## After Production

### 1. Combine Overlapping Channels

If both meson and EW files exist at same mass (4-8 GeV):

```bash
cd ../../analysis_pbc
python limits/combine_production_channels.py
```

This creates unified files combining B-meson and W-boson production.

### 2. Run Analysis

```bash
cd ../../analysis_pbc
python limits/run_serial.py
```

### 3. Generate Plots

```bash
cd ../money_plot
python plot_money_island.py
```

## Comparison: Sequential vs Parallel

### Sequential (run_full_production.sh)

```bash
for mass in "${MUON_MASSES[@]}"; do
    ./main_hnl_production $mass muon 100000  # Wait for completion
done
```

✅ **Pros:** Simple, safe, guaranteed order
❌ **Cons:** Very slow (~20 hours)

### Parallel (run_parallel_production.sh)

```bash
for mass in "${MUON_MASSES[@]}"; do
    wait_for_slot
    run_production_job $mass muon &  # Background execution
done
wait  # Wait for all jobs
```

✅ **Pros:** Fast (~3 hours), efficient CPU usage
❌ **Cons:** More complex, requires monitoring

## Best Practices

1. **Test first** (run a small subset or `test_parallel.sh` if present)
2. **Monitor disk space** during production
3. **Don't exceed** `(N_cores - 2)` for `MAX_PARALLEL`
4. **Check logs** after completion for failures
5. **Verify CSV files** before deleting logs
6. **Keep summary log** for documentation

## Troubleshooting Checklist

- [ ] Pythia executable compiled? (`ls -lh main_hnl_production`)
- [ ] Library path set? (`echo $DYLD_LIBRARY_PATH`)
- [ ] Output directory writable? (`touch ../../output/csv/simulation/test.txt`)
- [ ] Enough disk space? (`df -h output/`)
- [ ] Test script passes? (`./test_parallel.sh`, if available)

## References

- Sequential script: `run_full_production.sh`
- Test script: `test_parallel.sh`
- Analysis docs: `../../analysis_pbc/README.md`
- Configuration: `../../config_mass_grid.py`

---

**Created:** 2025-12-02
**Status:** Ready for production use
