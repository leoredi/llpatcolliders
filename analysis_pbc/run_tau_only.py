#!/usr/bin/env python
"""Run tau analysis only and update results."""
import sys
import argparse
import multiprocessing
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description="Run tau-only HNL analysis")
    parser.add_argument("--parallel", action="store_true", help="Use parallel processing")
    parser.add_argument("--workers", type=int, default=None, help="Number of workers (default: all CPU cores)")
    args = parser.parse_args()

    sys.path.insert(0, str(Path(__file__).parent.resolve()))
    from limits.run import run_flavour, GEOM_CACHE_DIR, ANALYSIS_OUT_DIR
    import pandas as pd

    GEOM_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    ANALYSIS_OUT_DIR.mkdir(parents=True, exist_ok=True)

    n_workers = args.workers if args.workers else multiprocessing.cpu_count()
    df = run_flavour(
        'tau', '001', 3000.0,
        use_parallel=args.parallel,
        n_workers=n_workers,
        separation_m=0.001,
        decay_seed=12345,
    )
    df['separation_mm'] = 1.0

    # Load existing results and update tau only
    results_file = ANALYSIS_OUT_DIR / 'HNL_U2_limits_summary.csv'
    if results_file.exists():
        existing = pd.read_csv(results_file)
        existing = existing[existing['flavour'] != 'tau']
        final = pd.concat([existing, df], ignore_index=True)
    else:
        final = df

    final.to_csv(results_file, index=False)
    print(f'Saved {len(df)} tau mass points to {results_file}')

if __name__ == '__main__':
    main()
