#!/usr/bin/env python3
"""cp358: PPPC falsification scan — p-indexed pre-window extension.

Extends the cp354/sec:empirical pre-window from k <= 500 to k <= KMAX
(default 1e6) over N_PRIMES exponents p* starting at START (default
5e11, the Platinum boundary, contiguous with the paper's window).

For each p* and k <= KMAX: candidate r = 2 p* k + 1 with
  (o) 2^{p*} == -1 (mod r)        [cheap powmod FIRST: implies r | 2^{p*}+1,
                                   and for prime r, ord_r(2) = 2p*]
  (i) r prime (BPSW, on (o)-survivors only)
  (ii) r == 3 (mod 8)
  (iii) d_r computed exactly (factor r+1, order descent of omega_3 in
        F_{r^2}), d_r > 43
  (iv) d_r odd, 4 d_r | r+1
  (v)  d_r | W_{p*-2}  <=>  2^{p*-2} == -1 (mod 3 d_r)   [THE event]

Any candidate passing (i)-(v) is a REAL danger triple (falsification-
relevant!); two at one p* would be a PPPC counterexample. Either
outcome is progress: a pass falsifies, a clean sweep extends the
empirical record by three orders of magnitude in k.

Pairs of (i)-(iv) candidates at one p* additionally get the
pair-separation asserts (thm:pair-separation 1-3) as a free real-data
test.

usage: cp358_falsification_prewindow.py [N_PRIMES] [KMAX] [START] [NPROC]
"""
import sys
from math import gcd
from multiprocessing import Pool

import sympy
from sympy import factorint, isprime, nextprime

N_PRIMES = int(sys.argv[1]) if len(sys.argv) > 1 else 50000
KMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 1_000_000
START = int(sys.argv[3]) if len(sys.argv) > 3 else 5 * 10**11
NPROC = int(sys.argv[4]) if len(sys.argv) > 4 else 12

WHEEL = [3, 5, 7, 11, 13, 17, 19, 23]


def d_r_of(r):
    """Exact ord(omega_3)/4 in F_{r^2} for inert r = 3 mod 8 dividing
    2^p+1: descent on (r+1)/4 using arithmetic in Z[sqrt2]/r."""
    def mul(a, b):
        return ((a[0]*b[0] + 2*a[1]*b[1]) % r, (a[0]*b[1] + a[1]*b[0]) % r)
    def powm(a, e):
        res = (1, 0)
        while e:
            if e & 1:
                res = mul(res, a)
            a = mul(a, a)
            e >>= 1
        return res
    n = (r + 1) // 4
    if powm((3, 2), 4 * n) != (1, 0):
        return None  # ord does not divide r+1 — not the danger shape
    g = n
    for q in factorint(n):
        while g % q == 0 and powm((3, 2), 4 * (g // q)) == (1, 0):
            g //= q
    return g


def scan_prime(p):
    """Return list of (i)-(iv) candidates (k, r, d_r, condv) at p."""
    out = []
    pm4 = p % 4
    for k in range(1, KMAX + 1):
        if k % 4 != pm4:        # C1: r = 3 mod 8 <=> k = p (mod 4)
            continue
        r = 2 * p * k + 1
        bad = False
        for q in WHEEL:
            if r % q == 0:
                bad = True
                break
        if bad:
            continue
        if pow(2, p, r) != r - 1:   # (o)
            continue
        if not isprime(r):          # (i)
            continue
        if r % 8 != 3:              # (ii) paranoia (C1 should ensure)
            continue
        d = d_r_of(r)               # (iii)
        if d is None or d <= 43 or d % 2 == 0:
            continue
        if (r + 1) % (4 * d) != 0:  # (iv)
            continue
        condv = pow(2, p - 2, 3 * d) == 3 * d - 1   # (v)
        out.append((k, r, d, condv))
    return (p, out)


def main():
    primes = []
    p = START - 1
    for _ in range(N_PRIMES):
        p = nextprime(p)
        primes.append(p)
    print(f"falsification pre-window: {N_PRIMES} primes from {primes[0]} "
          f"to {primes[-1]}, k <= {KMAX}", flush=True)
    n_cand = n_condv = n_pairs = 0
    with Pool(NPROC) as pool:
        for i, (p, cands) in enumerate(pool.imap_unordered(
                scan_prime, primes, chunksize=20)):
            if cands:
                n_cand += len(cands)
                for (k, r, d, condv) in cands:
                    line = (f"CAND p={p} k={k} r={r} d_r={d} "
                            f"condv={condv}")
                    print(line, flush=True)
                    if condv:
                        n_condv += 1
                        print(f"!!! (v)-PASS — REAL DANGER TRIPLE: {line}",
                              flush=True)
                if len(cands) >= 2:
                    n_pairs += 1
                    # free real-data separation asserts (thm:pair-separation)
                    for a in range(len(cands)):
                        for b in range(a + 1, len(cands)):
                            k1, r1, d1, _ = cands[a]
                            k2, r2, d2, _ = cands[b]
                            if k1 > k2:
                                k1, r1, d1, k2, r2, d2 = k2, r2, d2, k1, r1, d1
                            assert (k2 - k1) % 4 == 0 or (k1 % 4 == k2 % 4), \
                                f"SEP1 FAIL p={p}"
                            assert (k2 - k1) % gcd(p*k1+1, p*k2+1) == 0, \
                                f"SEP2 FAIL p={p} {k1} {k2}"
                            assert (k2 - k1) % (4 * gcd(d1, d2)) == 0, \
                                f"SEP3 FAIL p={p} {k1} {k2}"
                    print(f"PAIR at p={p}: {[(c[0], c[2]) for c in cands]} "
                          f"(separation asserts PASS)", flush=True)
            if (i + 1) % 2000 == 0:
                print(f"  progress {i+1}/{N_PRIMES} primes; cands={n_cand} "
                      f"condv={n_condv} pairs={n_pairs}", flush=True)
    print(f"FALSIFICATION SCAN DONE: {N_PRIMES} primes x k<={KMAX}; "
          f"(i)-(iv) candidates={n_cand}; (v)-passes={n_condv}; "
          f"multi-candidate exponents={n_pairs}", flush=True)


if __name__ == "__main__":
    main()
