#!/usr/bin/env python3
"""cp358: corner-census extension of cp354_diag_corner_gcd.py part C:
admissible odd d in (124, 1000], condition-(v) primes p <= 1e8.
Factorization-free diagonal census via the four corner gcds
G(eps,eps') = gcd(V_d - eps*2^{(p+1)/2}, 2U_d - eps') (cp354 L2,
PROVED UNCONDITIONALLY; gate-verified). Reports every corner hit,
QUALifying danger events (r = 3 mod 8, primitive at d, d > 43), and
any diagonal PAIR (>= 2 qualifying r at one (d, p)).

Speedups vs cp354: (1) condition (v) as a single residue class
p == 2 + ord/2 (mod ord), ord = ord_{3d}(2), instead of a powmod per
prime; (2) multiprocessing over d.

usage: cp358_corner_sweep_d1000.py [D_LO D_HI P_MAX NPROC]
"""
import json
import sys
from math import gcd
from multiprocessing import Pool

import sympy
from sympy import factorint, isprime, primerange

D_LO = int(sys.argv[1]) if len(sys.argv) > 1 else 125
D_HI = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
P_MAX = int(sys.argv[3]) if len(sys.argv) > 3 else 100_000_000
NPROC = int(sys.argv[4]) if len(sys.argv) > 4 else 12


def pell_VU(n):
    V0, V1, U0, U1 = 2, 2, 0, 1
    if n == 0:
        return V0, U0
    for _ in range(n - 1):
        V0, V1 = V1, 2 * V1 + V0
        U0, U1 = U1, 2 * U1 + U0
    return V1, U1


def in_class_A(q):
    o, t = 1, 2 % q
    while t != 1:
        t = (t * 2) % q
        o += 1
    return o % 2 == 0 and (o // 2) % 2 == 1 and pow(2, o // 2, q) == q - 1


def admissible(d):
    return all(in_class_A(q) for q in factorint(d))


def corners(V_d, U_d, m):
    out = {}
    for epsp in (+1, -1):
        Hr = 2 * U_d - epsp
        t = pow(2, (m + 1) // 2, Hr)
        for eps in (+1, -1):
            g = gcd((V_d - eps * t) % Hr, Hr)
            while g % 2 == 0:
                g //= 2
            while g % 3 == 0:
                g //= 3
            out[(eps, epsp)] = g
    return out


def classify_corner(g, m):
    if g == 1:
        return []
    if g < (2 * m + 1) ** 2 or isprime(g):
        return [(g, 1)]
    return sorted(factorint(g).items())


def primitive_at_d(r, d):
    for e in range(1, d):
        if d % e == 0:
            Ve, _ = pell_VU(e)
            if (Ve * Ve + 2) % r == 0:
                return False
    return True


PRIMES = None  # built once in main() before the Pool fork (COW-shared)


def sweep_d(d):
    """Full corner census at one d; returns (d, npairs_tested, hits)."""
    V_d, U_d = pell_VU(d)
    n3d = 3 * d
    # condition (v): 2^{p-2} == -1 mod 3d  <=>  p == 2 + ord/2 (mod ord)
    ord3d = sympy.n_order(2, n3d)
    if ord3d % 2 or pow(2, ord3d // 2, n3d) != n3d - 1:
        return (d, 0, [])
    c = (2 + ord3d // 2) % ord3d
    hits = []
    tested = 0
    for p in PRIMES:
        if p % ord3d != c:
            continue
        # paranoia cross-check on a small sample
        if tested < 3:
            assert pow(2, p - 2, n3d) == n3d - 1
        tested += 1
        for key, g in corners(V_d, U_d, p).items():
            if g > 1:
                for r, e in classify_corner(g, p):
                    qual = (r % 8 == 3) and primitive_at_d(r, d)
                    hits.append((d, p, key, str(r), e, r % 8, qual,
                                 r > 10**12))
    return (d, tested, hits)


def main():
    global PRIMES
    adm = [d for d in range(D_LO + (D_LO % 2 == 0), D_HI + 1, 2)
           if admissible(d)]
    print(f"admissible odd d in ({D_LO},{D_HI}]: {len(adm)} values",
          flush=True)
    print("sieving primes...", flush=True)
    PRIMES = list(primerange(5, P_MAX))
    print(f"  {len(PRIMES)} primes < {P_MAX}", flush=True)
    total_pairs = 0
    all_hits = []
    with Pool(NPROC) as pool:
        for d, tested, hits in pool.imap_unordered(sweep_d, adm):
            total_pairs += tested
            all_hits.extend(hits)
            print(f"  d={d}: {tested} cond-v primes <= {P_MAX}, "
                  f"{len(hits)} corner hits", flush=True)
    print(f"\nTOTAL: {total_pairs} (d,p) rectangles, {len(all_hits)} hits",
          flush=True)
    qual_by_dp = {}
    for (d, p, key, r, e, m8, qual, big) in all_hits:
        print(f"  HIT d={d} p={p} corner({key[0]:+d},{key[1]:+d}) "
              f"r={r}^{e} mod8={m8} qual={qual} big={big}")
        if qual:
            qual_by_dp.setdefault((d, p), []).append(r)
    pairs = {k: v for k, v in qual_by_dp.items() if len(v) >= 2}
    print(f"\nQUAL events: {sum(len(v) for v in qual_by_dp.values())}; "
          f"DIAGONAL PAIRS: {len(pairs)}")
    for k, v in pairs.items():
        print("  PAIR!!", k, v)
    with open("scripts/output/cp358_corner_sweep_summary.json", "w") as f:
        json.dump({"D_LO": D_LO, "D_HI": D_HI, "P_MAX": P_MAX,
                   "admissible": adm, "rectangles": total_pairs,
                   "hits": [(h[0], h[1], list(h[2]), h[3], h[4], h[5],
                             h[6], h[7]) for h in all_hits],
                   "pairs": [[list(k), v] for k, v in pairs.items()]},
                  f, indent=1)
    print("summary written; DONE")


if __name__ == "__main__":
    main()
