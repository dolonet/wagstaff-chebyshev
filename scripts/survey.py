#!/usr/bin/env python3
"""
Wagstaff inert-factor survey.

For every odd prime p for which W_p = (2^p+1)/3 is composite, this script
finds all inert prime factors r ≡ 3 (mod 8) with r | W_p in the specified
range and records the following per factor:

  1. G_r/4 = gcd((r+1)/4, W_{p-2})  — computed by modular arithmetic only.
  2. If G_r/4 ≤ 43 the factor is recorded as "proved" (Theorem 7.5 /
     N_d exclusion rules out any local pass of Condition (II) for d ≤ 43).
  3. Otherwise a full structural analysis is performed: order of ω₃ = 3+2√2
     modulo r, factorization of d = ord(ω₃)/4, Class I/II/III classification
     of each prime in d, and identification of the obstruction that blocks
     Condition (II) locally.

The G_r/4 ≤ 43 shortcut is valid because any local pass of Condition (II)
requires d | G_r/4, so d ≤ G_r/4. The N_d exclusion theorem rules out any
pass for d ≤ 43; therefore every potential pass has G_r/4 > 43 and only
~2.3% of factors need the expensive order computation.

Output (inside --output-dir):
  inert_factors.csv    Per-factor: p,r,G4,d,d_factors,pure_iii,obstruction,blocking_q,blocking_class
  segment_stats.jsonl  Per-segment aggregates
  checkpoint.json      Resume state
  survey.log           Timestamped progress + any local-pass alerts
  params.json          Run parameters (for reproducibility)
  summary.json         Final totals

Usage:
  python3 survey.py --cores 128 --max-r 1e12 --output-dir survey_output/
  python3 survey.py --cores 128 --output-dir survey_output/ --resume

The script is pure-Python and has no external dependencies beyond the
standard library.
"""

import argparse
import json
import math
import multiprocessing as mp
import os
import sys
import time
from datetime import datetime, timezone

VERSION = "2.1.0"

# ═══════════════════════════════════════════════════════════════
# Known Wagstaff PRP exponents — skip these (W_p is prime/PRP)
# ═══════════════════════════════════════════════════════════════

WAGSTAFF_PRP = frozenset({
    3, 5, 7, 11, 13, 17, 19, 23, 31, 43, 61, 79, 101, 127,
    167, 191, 199, 313, 347, 701, 1709, 2617, 3539, 5807,
    10501, 10691, 11279, 12391, 14479, 42737, 83339, 95369,
    117239, 127031, 138937, 141079, 267017, 269987, 374321,
    986191, 4031399, 13347311
})

# ═══════════════════════════════════════════════════════════════
# Primality testing and factoring (pure Python, no dependencies)
# ═══════════════════════════════════════════════════════════════

_SMALL_PRIMES = (2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97)
_SMALL_SET = frozenset(_SMALL_PRIMES)

def _miller_rabin(n, witnesses=(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)):
    """Deterministic for n < 3.3×10^24."""
    d, s = n - 1, 0
    while d % 2 == 0: d >>= 1; s += 1
    for a in witnesses:
        if a >= n: continue
        x = pow(a, d, n)
        if x == 1 or x == n - 1: continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1: break
        else:
            return False
    return True

def is_prime(n):
    if n < 2: return False
    if n in _SMALL_SET: return True
    if any(n % p == 0 for p in _SMALL_PRIMES): return False
    return _miller_rabin(n)

def _pollard_rho(n):
    if n % 2 == 0: return 2
    for c in range(1, 256):
        x = y = c + 1; d = 1
        while d == 1:
            x = (x * x + c) % n
            y = (y * y + c) % n; y = (y * y + c) % n
            d = math.gcd(abs(x - y), n)
        if d != n: return d
    return n

def factor(n):
    """Complete factorization → {prime: exp}."""
    if n <= 1: return {}
    result = {}
    for p in _SMALL_PRIMES:
        while n % p == 0: result[p] = result.get(p, 0) + 1; n //= p
    if n == 1: return result
    d = 101
    while d * d <= n and d < 100_000:
        while n % d == 0: result[d] = result.get(d, 0) + 1; n //= d
        d += 2
    if n == 1: return result
    stack = [n]
    while stack:
        m = stack.pop()
        if m <= 1: continue
        if is_prime(m): result[m] = result.get(m, 0) + 1; continue
        f = _pollard_rho(m)
        if f == m: result[m] = result.get(m, 0) + 1
        else: stack.append(f); stack.append(m // f)
    return result

# ═══════════════════════════════════════════════════════════════
# Segmented prime sieve
# ═══════════════════════════════════════════════════════════════

def sieve_small(limit):
    """Primes up to limit."""
    sieve = bytearray(b'\x01') * (limit + 1)
    sieve[0] = sieve[1] = 0
    for i in range(2, math.isqrt(limit) + 1):
        if sieve[i]:
            sieve[i*i::i] = bytearray(len(sieve[i*i::i]))
    return [i for i in range(2, limit + 1) if sieve[i]]

def primes_in_range(lo, hi, small_primes):
    """Primes in [lo, hi) via segmented sieve."""
    if lo >= hi: return []
    sz = hi - lo
    sieve = bytearray(b'\x01') * sz
    for p in small_primes:
        if p * p >= hi: break
        start = (-(lo % p)) % p
        if lo <= p: start = p * p - lo
        if start < sz:
            sieve[start::p] = bytearray(len(sieve[start::p]))
    for i in range(max(0, 2 - lo)):
        if i < sz: sieve[i] = 0
    return [lo + i for i in range(sz) if sieve[i]]

# ═══════════════════════════════════════════════════════════════
# T₂ arithmetic and order computation
# ═══════════════════════════════════════════════════════════════

def pow_T2(n, r):
    """ω₃^n = (3+2√2)^n mod r → (a, b) with a+b√2."""
    a, b = 1, 0
    ca, cb = 3 % r, 2 % r
    while n > 0:
        if n & 1:
            a, b = (a*ca + 2*b*cb) % r, (a*cb + b*ca) % r
        ca, cb = (ca*ca + 2*cb*cb) % r, 2*ca*cb % r
        n >>= 1
    return a, b

def ord_omega3(r):
    """Order of ω₃ in norm-1 subgroup of (Z[√2]/rZ)*, divides r+1."""
    rp1 = r + 1
    facs = factor(rp1)
    order = rp1
    for p, e in facs.items():
        for _ in range(e):
            t = order // p
            if pow_T2(int(t), int(r)) == (1, 0):
                order = t
            else:
                break
    return int(order)

# ═══════════════════════════════════════════════════════════════
# Classification cache (per worker process)
# ═══════════════════════════════════════════════════════════════

_cls_cache = {}

def classify(q):
    """Classify prime q by v₂(ord_q(2)).
    Returns (class, ord_q(2), m).
      I:   v₂=0 (odd order)
      II:  v₂≥2
      III: v₂=1, m=ord/2
    """
    if q in _cls_cache: return _cls_cache[q]
    if q == 2:
        r = ("special", 1, 0)
    else:
        g = q - 1; facs = factor(g); o = g
        for p, e in facs.items():
            for _ in range(e):
                if pow(2, o // p, q) == 1: o //= p
                else: break
        v2 = 0; tmp = o
        while tmp % 2 == 0: v2 += 1; tmp >>= 1
        if v2 == 0:   r = ("I", o, 0)
        elif v2 == 1:  r = ("III", o, o >> 1)
        else:          r = ("II", o, 0)
    _cls_cache[q] = r
    return r

# ═══════════════════════════════════════════════════════════════
# G_r/4 via modular arithmetic (no big-integer W_{p-2})
# ═══════════════════════════════════════════════════════════════

def compute_G4(p, r):
    """G_r/4 = gcd((r+1)/4, W_{p-2}).
    Uses: 3·W_{p-2} = 2^{p-2}+1 and 3·W mod 3m = 3·(W mod m)."""
    m = (r + 1) // 4
    if m <= 1: return m
    three_m = 3 * m
    t = (pow(2, p - 2, three_m) + 1) % three_m
    # t = 3·W_{p-2} mod 3m, which equals 3·(W_{p-2} mod m)
    # Safe because 3 | (2^{p-2}+1) for odd p-2, i.e. all p≥5
    return math.gcd(m, t // 3)

def _vq(n, q):
    if n == 0: return 99
    v = 0
    while n % q == 0: v += 1; n //= q
    return v

# ═══════════════════════════════════════════════════════════════
# Worker
# ═══════════════════════════════════════════════════════════════

_CFG = {}

def _worker_init(cfg):
    global _CFG, _cls_cache
    _CFG = cfg
    _cls_cache = {}

def _check_p(p):
    """For prime p with composite W_p, find inert factors in [min_r, max_r]
    and return list of result tuples."""
    max_r = _CFG['max_r']
    min_r = _CFG['min_r']
    max_k_iter = _CFG['max_k_iter']

    k_lo = max(1, (min_r // (2 * p)) + (1 if min_r > 0 else 0))
    k_hi = min(max_r // (2 * p), max_k_iter)
    if k_hi < k_lo: return ([], 0)

    # Optimization: only iterate k values giving r ≡ 3 mod 8.
    # r = 2pk+1, need 2pk ≡ 2 mod 8, i.e. pk ≡ 1 mod 4.
    # p odd prime: p%4 ∈ {1,3}. Need k ≡ 1/p mod 4.
    p4 = p % 4
    k_res = 1 if p4 == 1 else 3  # k ≡ k_res mod 4 gives r ≡ 3 mod 8
    k_start = k_lo + (k_res - k_lo % 4) % 4
    if k_start < k_lo: k_start += 4

    rows = []
    g4_le43 = 0
    for ki in range(k_start, k_hi + 1, 4):
        r = 2 * p * ki + 1
        if r > max_r: break
        if not is_prime(r): continue
        if pow(2, p, r) != r - 1: continue

        # ── Inert factor found ──
        G4 = compute_G4(p, r)

        # G4 ≤ 43 shortcut: N_d exclusion proves no local pass possible
        # when d ≤ 43, and any local pass requires d | G4 ≤ 43.
        if G4 <= 43:
            g4_le43 += 1
            continue

        # G4 > 43: quick test before expensive ord computation.
        # For a local pass we need d | G4, hence ord(ω₃) = 4d | 4·G4.
        # One cheap exponentiation rules out ~99.99% of cases.
        if pow_T2(4 * G4, r) != (1, 0):
            rows.append((p, r, G4, 0, '', False, 'order_excess', 0, ''))
            continue

        # ω₃^{4·G4} = 1, so d | G4. Full analysis to find exact d and obstruction.
        order = ord_omega3(r)
        d = order // 4

        if d <= 1:
            rows.append((p, r, G4, d, '', False, 'trivial', 0, ''))
            continue

        d_facs = factor(d)
        d_str = '*'.join(f'{q}^{e}' if e > 1 else str(q)
                         for q, e in sorted(d_facs.items()) if q > 1)

        # Classify all odd prime factors
        classes = {}
        for q in d_facs:
            if q > 2: classes[q] = classify(q)

        pure_iii = all(c[0] == "III" for c in classes.values())

        # Determine obstruction (first blocking prime)
        obs = None; bq = 0; bc = ''
        for q, a in sorted(d_facs.items()):
            if q <= 2: continue
            cls, oq, mq = classes[q]
            if cls in ("I", "II"):
                obs = f'class_{cls}'; bq = q; bc = cls; break
            # cls == "III": check if q | W_{p-2}
            if pow(2, p, q) != (-4) % q:
                obs = 'congruence'; bq = q; bc = 'III'; break
            if a >= 2:
                if q == 3:
                    if _vq(p - 2, 3) < a:
                        obs = 'valuation'; bq = q; bc = 'III'; break
                else:
                    if (p - 2) % mq != 0:
                        obs = 'valuation'; bq = q; bc = 'III'; break
                    if 1 + _vq((p - 2) // mq, q) < a:
                        obs = 'valuation'; bq = q; bc = 'III'; break

        if obs is None:
            obs = 'PASS'  # Local pass of Condition (II). Conjecturally empty.

        rows.append((p, r, G4, d, d_str, pure_iii, obs, bq, bc))
    return (rows, g4_le43)

# ═══════════════════════════════════════════════════════════════
# Checkpoint and I/O
# ═══════════════════════════════════════════════════════════════

CSV_HEADER = "p,r,G4,d,d_factors,pure_iii,obstruction,blocking_q,blocking_class\n"

def load_checkpoint(path):
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    return None

def save_checkpoint(path, data):
    tmp = path + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(data, f, indent=2)
    os.replace(tmp, path)

def log(logfile, msg):
    ts = datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M:%S')
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    if logfile:
        logfile.write(line + '\n')
        logfile.flush()

# ═══════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description='Wagstaff inert factor survey')
    parser.add_argument('--cores', type=int, default=mp.cpu_count())
    parser.add_argument('--max-r', type=float, default=1e10)
    parser.add_argument('--min-r', type=float, default=0)
    parser.add_argument('--max-k-iter', type=int, default=1_000_000,
                        help='Cap on k iterations per prime (limits work for small p)')
    parser.add_argument('--segment-size', type=int, default=10_000_000)
    parser.add_argument('--output-dir', type=str, default='survey_output')
    parser.add_argument('--resume', action='store_true')
    args = parser.parse_args()

    max_r = int(args.max_r)
    min_r = int(args.min_r)
    cores = args.cores
    seg_size = args.segment_size
    outdir = args.output_dir

    # Effective max_p: r = 2p·k+1 ≥ 2p+1, so p < max_r/2
    max_p = (max_r - 1) // 2

    os.makedirs(outdir, exist_ok=True)
    csv_path = os.path.join(outdir, 'inert_factors.csv')
    stats_path = os.path.join(outdir, 'segment_stats.jsonl')
    ckpt_path = os.path.join(outdir, 'checkpoint.json')
    log_path = os.path.join(outdir, 'survey.log')
    params_path = os.path.join(outdir, 'params.json')

    logf = open(log_path, 'a')

    # Save parameters
    params = {
        'version': VERSION,
        'max_r': max_r, 'min_r': min_r, 'max_p_effective': max_p,
        'max_k_iter': args.max_k_iter, 'segment_size': seg_size,
        'cores': cores, 'started': datetime.now(timezone.utc).isoformat()
    }
    with open(params_path, 'w') as f:
        json.dump(params, f, indent=2)

    # Resume state
    ckpt = load_checkpoint(ckpt_path) if args.resume else None
    resume_seg = ckpt['last_segment'] + 1 if ckpt else 0
    totals = ckpt.get('totals', {}) if ckpt else {}
    total_primes = totals.get('primes_checked', 0)
    total_inert = totals.get('inert_found', 0)
    total_pure_iii = totals.get('pure_iii', 0)
    total_passes = totals.get('passes', 0)
    total_g4_le_43 = totals.get('g4_le_43', 0)

    num_segments = math.ceil(max_p / seg_size)

    log(logf, f"=== Wagstaff Inert Factor Survey v{VERSION} ===")
    log(logf, f"Range: min_r={min_r:,.0f}  max_r={max_r:,.0f}  max_p={max_p:,.0f}")
    log(logf, f"Cores: {cores}  Segments: {num_segments}  Segment size: {seg_size:,}")
    if ckpt:
        log(logf, f"Resuming from segment {resume_seg} "
            f"(primes={total_primes:,} inert={total_inert:,} passes={total_passes})")

    # ── Precompute small primes for sieve ──
    sieve_limit = math.isqrt(max_p) + 1
    log(logf, f"Sieving small primes up to {sieve_limit:,}...")
    small_primes = sieve_small(sieve_limit)
    log(logf, f"  {len(small_primes)} small primes ready")

    # ── Open output files ──
    csv_mode = 'a' if args.resume else 'w'
    csvf = open(csv_path, csv_mode)
    if not args.resume:
        csvf.write(CSV_HEADER)
        csvf.flush()

    statsf = open(stats_path, 'a' if args.resume else 'w')

    # ── Worker config ──
    cfg = {
        'max_r': max_r,
        'min_r': min_r,
        'max_k_iter': args.max_k_iter,
    }

    # ── Main loop ──
    pool = mp.Pool(cores, initializer=_worker_init, initargs=(cfg,))
    t_start = time.time()
    t_survey_start = t_start

    try:
        for seg_idx in range(num_segments):
            seg_lo = seg_idx * seg_size
            seg_hi = min(seg_lo + seg_size, max_p + 1)

            if seg_idx < resume_seg:
                continue
            if seg_lo > max_p:
                break

            t_seg = time.time()

            # Generate primes in this segment
            primes = primes_in_range(seg_lo, seg_hi, small_primes)
            # Filter: skip <5, skip Wagstaff PRPs
            primes = [p for p in primes if p >= 5 and p not in WAGSTAFF_PRP]

            if not primes:
                # Save checkpoint even for empty segments
                save_checkpoint(ckpt_path, {
                    'last_segment': seg_idx, 'totals': {
                        'primes_checked': total_primes, 'inert_found': total_inert,
                        'pure_iii': total_pure_iii, 'passes': total_passes,
                        'g4_le_43': total_g4_le_43,
                    }
                })
                continue

            # Small chunksize for good load balancing (small p are 1000× heavier)
            chunksize = max(1, min(32, len(primes) // (cores * 16)))

            # Process
            seg_inert = 0
            seg_pure = 0
            seg_passes = 0
            seg_g4_le_43 = 0
            seg_g4_max = 0
            seg_rows = 0
            obs_counts = {}

            for result in pool.imap_unordered(_check_p, primes, chunksize=chunksize):
                rows, g4_le43_count = result
                seg_g4_le_43 += g4_le43_count
                seg_inert += g4_le43_count
                obs_counts['proved'] = obs_counts.get('proved', 0) + g4_le43_count
                for row in rows:
                    # row = (p, r, G4, d, d_str, pure_iii, obs, bq, bc)
                    p_val, r_val, G4, d, d_str, pure_iii, obs, bq, bc = row
                    seg_inert += 1
                    if pure_iii and d > 1: seg_pure += 1
                    if obs == 'PASS': seg_passes += 1
                    if G4 > seg_g4_max: seg_g4_max = G4
                    obs_counts[obs] = obs_counts.get(obs, 0) + 1
                    seg_rows += 1

                    # Write CSV row
                    csvf.write(f"{p_val},{r_val},{G4},{d},{d_str},"
                               f"{pure_iii},{obs},{bq},{bc}\n")

                    # Alert on local pass
                    if obs == 'PASS':
                        log(logf, f"*** CONDITION II LOCAL PASS: p={p_val} r={r_val} "
                            f"G4={G4} d={d} d_facs={d_str} ***")

            csvf.flush()

            # Update totals
            total_primes += len(primes)
            total_inert += seg_inert
            total_pure_iii += seg_pure
            total_passes += seg_passes
            total_g4_le_43 += seg_g4_le_43

            # Segment stats
            seg_elapsed = time.time() - t_seg
            seg_stat = {
                'segment': seg_idx,
                'p_range': [int(seg_lo), int(seg_hi)],
                'primes': len(primes),
                'inert': seg_inert,
                'pure_iii': seg_pure,
                'passes': seg_passes,
                'g4_le_43': seg_g4_le_43,
                'g4_max': seg_g4_max,
                'obstruction_counts': obs_counts,
                'elapsed_s': round(seg_elapsed, 1),
            }
            statsf.write(json.dumps(seg_stat) + '\n')
            statsf.flush()

            # Checkpoint
            save_checkpoint(ckpt_path, {
                'last_segment': seg_idx,
                'totals': {
                    'primes_checked': total_primes,
                    'inert_found': total_inert,
                    'pure_iii': total_pure_iii,
                    'passes': total_passes,
                    'g4_le_43': total_g4_le_43,
                },
            })

            # Progress
            elapsed_total = time.time() - t_survey_start
            segs_done = seg_idx - resume_seg + 1
            segs_left = num_segments - seg_idx - 1
            if segs_done > 0 and segs_left > 0:
                eta_s = elapsed_total / segs_done * segs_left
                eta_str = f"  ETA: {eta_s/3600:.1f}h"
            else:
                eta_str = ""

            if (seg_idx + 1) % 5 == 0 or seg_inert > 0:
                log(logf, f"Seg {seg_idx}/{num_segments}: "
                    f"p=[{seg_lo:,},{seg_hi:,}) "
                    f"primes={len(primes):,} inert={seg_inert:,} "
                    f"pure_iii={seg_pure} passes={seg_passes} "
                    f"g4_max={seg_g4_max} "
                    f"[{seg_elapsed:.1f}s]"
                    f"  TOTAL: {total_inert:,} inert, {total_passes} passes"
                    f"{eta_str}")

    except KeyboardInterrupt:
        log(logf, f"Interrupted at segment {seg_idx}. Resume with --resume.")
    finally:
        pool.terminate()
        pool.join()
        csvf.close()
        statsf.close()

    # ── Final summary ──
    elapsed = time.time() - t_survey_start
    summary = {
        'version': VERSION,
        'max_r': max_r,
        'min_r': min_r,
        'max_p_effective': max_p,
        'total_primes_checked': total_primes,
        'total_inert_factors': total_inert,
        'total_pure_class_iii': total_pure_iii,
        'total_condition_ii_passes': total_passes,
        'total_g4_le_43': total_g4_le_43,
        'fraction_g4_le_43': round(total_g4_le_43 / max(total_inert, 1), 6),
        'fraction_pure_iii': round(total_pure_iii / max(total_inert, 1), 6),
        'elapsed_hours': round(elapsed / 3600, 2),
        'completed': datetime.now(timezone.utc).isoformat(),
    }
    summary_path = os.path.join(outdir, 'summary.json')
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)

    log(logf, "")
    log(logf, "═" * 72)
    log(logf, "SURVEY COMPLETE")
    log(logf, f"  Primes checked:        {total_primes:>15,}")
    log(logf, f"  Inert factors found:   {total_inert:>15,}")
    log(logf, f"  Pure Class III:        {total_pure_iii:>15,}")
    log(logf, f"  Condition II passes:   {total_passes:>15,}")
    log(logf, f"  G_r/4 ≤ 43:           {total_g4_le_43:>15,} "
        f"({100*total_g4_le_43/max(total_inert,1):.2f}%)")
    log(logf, f"  Elapsed:               {elapsed/3600:>14.2f}h")
    log(logf, "═" * 72)

    if total_passes > 0:
        log(logf, f"!!! WARNING: {total_passes} CONDITION II LOCAL PASSES DETECTED !!!")
        log(logf, "Check survey.log and inert_factors.csv for details.")
    else:
        log(logf, "CONFIRMED: 0 Condition II local passes.")

    logf.close()
    return 0 if total_passes == 0 else 1

if __name__ == '__main__':
    sys.exit(main())
