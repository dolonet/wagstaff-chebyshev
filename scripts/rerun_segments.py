#!/usr/bin/env python3
"""
Re-run a specific list of inert-factor survey segments.

This is a re-entrant, fsync-strict companion to ``survey.py``. It
re-runs an explicit list of segment ids — useful if an interrupted
survey left any segment with mismatched row counts in the audit report,
or for spot-checking a subrange deterministically.

Design points (differ deliberately from ``survey.py``):

  * Opens the output CSV in ``'w'`` mode only — never appends.
  * Calls ``os.fsync()`` on the CSV and stats file descriptors after
    every segment, so an abrupt SIGKILL cannot leave buffered rows
    unflushed.
  * Takes an explicit ``--segments`` list (e.g. ``"4001,4002,4003"``).
  * Writes to a fresh ``--output-dir`` that must not already exist,
    which guarantees it cannot pollute an existing survey output.

Internally it imports ``_check_p`` from ``survey.py``, so the
arithmetic is bit-for-bit identical to the original run.

Usage::

    python3 rerun_segments.py \
        --segments 4001,4002,4003 \
        --segment-size 100000000 \
        --max-r 1e12 \
        --cores 128 \
        --output-dir rerun/
"""
import argparse
import json
import math
import multiprocessing as mp
import os
import sys
import time
from datetime import datetime, timezone

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from survey import (  # noqa: E402
    _check_p, _worker_init, primes_in_range, sieve_small,
    WAGSTAFF_PRP, CSV_HEADER, VERSION,
)


def _fsync(f):
    f.flush()
    os.fsync(f.fileno())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--segments', required=True,
                    help='Comma-separated segment ids to re-run (e.g. "4001,4002,4003")')
    ap.add_argument('--segment-size', type=int, default=100_000_000)
    ap.add_argument('--max-r', type=float, default=1e12)
    ap.add_argument('--min-r', type=float, default=0)
    ap.add_argument('--max-k-iter', type=int, default=1_000_000)
    ap.add_argument('--cores', type=int, default=mp.cpu_count())
    ap.add_argument('--output-dir', required=True,
                    help='Fresh output directory (must not exist)')
    args = ap.parse_args()

    segments = sorted({int(s) for s in args.segments.split(',') if s.strip()})
    seg_size = args.segment_size
    max_r = int(args.max_r)
    min_r = int(args.min_r)
    cores = args.cores
    outdir = args.output_dir

    if os.path.exists(outdir):
        print(f"ERROR: output dir {outdir!r} already exists — refuse to overwrite.",
              file=sys.stderr)
        sys.exit(2)
    os.makedirs(outdir)

    csv_path = os.path.join(outdir, 'inert_factors.csv')
    stats_path = os.path.join(outdir, 'segment_stats.jsonl')
    params_path = os.path.join(outdir, 'params.json')
    log_path = os.path.join(outdir, 'rerun.log')

    max_p = (max_r - 1) // 2

    params = {
        'version': VERSION,
        'driver': 'rerun_segments.py',
        'segments': segments,
        'segment_size': seg_size,
        'max_r': max_r,
        'min_r': min_r,
        'max_p_effective': max_p,
        'max_k_iter': args.max_k_iter,
        'cores': cores,
        'started': datetime.now(timezone.utc).isoformat(),
    }
    with open(params_path, 'w') as f:
        json.dump(params, f, indent=2)

    max_seg_hi = max(s + 1 for s in segments) * seg_size
    sieve_limit = math.isqrt(max_seg_hi) + 1

    logf = open(log_path, 'w')
    def log(msg):
        ts = datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M:%S')
        line = f"[{ts}] {msg}"
        print(line, flush=True)
        logf.write(line + '\n')
        logf.flush()

    log(f"rerun_segments.py  segments={segments}  cores={cores}  max_r={max_r:,}")
    log(f"Sieving small primes up to {sieve_limit:,}...")
    small_primes = sieve_small(sieve_limit)
    log(f"  {len(small_primes):,} small primes ready")

    cfg = {'max_r': max_r, 'min_r': min_r, 'max_k_iter': args.max_k_iter}

    csvf = open(csv_path, 'w')
    csvf.write(CSV_HEADER)
    _fsync(csvf)
    statsf = open(stats_path, 'w')

    pool = mp.Pool(cores, initializer=_worker_init, initargs=(cfg,))
    t0 = time.time()

    try:
        for seg in segments:
            seg_lo = seg * seg_size
            seg_hi = min(seg_lo + seg_size, max_p + 1)
            if seg_lo > max_p:
                log(f"  seg {seg}: above max_p, skipping")
                continue

            t_seg = time.time()
            primes = primes_in_range(seg_lo, seg_hi, small_primes)
            primes = [p for p in primes if p >= 5 and p not in WAGSTAFF_PRP]

            chunksize = max(1, min(32, len(primes) // (cores * 16)))

            seg_inert = 0
            seg_pure = 0
            seg_passes = 0
            seg_g4_le_43 = 0
            seg_g4_max = 0
            obs_counts = {}

            for rows, g4_le43_count in pool.imap_unordered(_check_p, primes, chunksize=chunksize):
                seg_g4_le_43 += g4_le43_count
                seg_inert += g4_le43_count
                if g4_le43_count:
                    obs_counts['proved'] = obs_counts.get('proved', 0) + g4_le43_count
                for row in rows:
                    p_val, r_val, G4, d, d_str, pure_iii, obs, bq, bc = row
                    seg_inert += 1
                    if pure_iii and d > 1:
                        seg_pure += 1
                    if obs == 'PASS':
                        seg_passes += 1
                    if G4 > seg_g4_max:
                        seg_g4_max = G4
                    obs_counts[obs] = obs_counts.get(obs, 0) + 1
                    csvf.write(f"{p_val},{r_val},{G4},{d},{d_str},"
                               f"{pure_iii},{obs},{bq},{bc}\n")
                    if obs == 'PASS':
                        log(f"*** LOCAL PASS: p={p_val} r={r_val} G4={G4} d={d} d_facs={d_str} ***")

            _fsync(csvf)

            seg_stat = {
                'segment': seg,
                'p_range': [int(seg_lo), int(seg_hi)],
                'primes': len(primes),
                'inert': seg_inert,
                'pure_iii': seg_pure,
                'passes': seg_passes,
                'g4_le_43': seg_g4_le_43,
                'g4_max': seg_g4_max,
                'obstruction_counts': obs_counts,
                'elapsed_s': round(time.time() - t_seg, 1),
            }
            statsf.write(json.dumps(seg_stat) + '\n')
            _fsync(statsf)

            log(f"Seg {seg}: primes={len(primes):,} inert={seg_inert:,} "
                f"g4_le_43={seg_g4_le_43:,} g4_max={seg_g4_max} "
                f"rows_written={seg_inert - seg_g4_le_43:,} "
                f"[{seg_stat['elapsed_s']}s]")
    finally:
        pool.terminate()
        pool.join()
        csvf.close()
        statsf.close()

    summary = {
        **params,
        'completed': datetime.now(timezone.utc).isoformat(),
        'elapsed_s': round(time.time() - t0, 1),
    }
    with open(os.path.join(outdir, 'summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)

    log(f"DONE in {time.time() - t0:.1f}s")
    logf.close()


if __name__ == '__main__':
    main()
