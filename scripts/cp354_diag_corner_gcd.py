#!/usr/bin/env python3
"""
GATE NOTE (cp354): the corner-prime classification
("all prime factors of a corner are = 1 mod 2m", hence "corner <
(2m+1)^2 is 1 or prime") requires m PRIME (after stripping 3s);
for composite m there are counterexamples (e.g. d=15, m=95, r=11).
All uses in this project have m = p prime.

cp354_diag_corner_gcd.py — Corner-gcd localization of the PPPC diagonal,
and a (d, p)-sweep diagonal census that needs NO factorization of V_{2d}.

MATH (cp354 session; claimed PROVED UNCONDITIONALLY, elementary):

  d odd >= 3, m odd, V_{2d} = V_d^2 + 2 = 2(2U_d-1)(2U_d+1).

  (L1) HALF-PINNING: r > 3 prime, r | V_{2d}, ord_r(2) = 2m, m odd
       ==> V_d = eps*2^{(m+1)/2} mod r for a unique eps in {+-1}.

  (L2) CORNER LOCALIZATION: with G(eps,eps') = gcd(V_d - eps*2^{(m+1)/2},
       2U_d - eps'), a prime r > 3 divides V_{2d} with 2^m = -1 mod r
       IFF r | G(eps,eps') for exactly one corner. All prime factors of
       every corner are = 1 mod 2m; so a corner < (2m+1)^2 is 1 or prime.

  (L3) PRODUCT CAPS for same-(d,m) pairs: same eps' => r1 r2 <= 2U_d+1;
       same eps => r1 r2 <= |V_d - eps 2^{(m+1)/2}|; anti-aligned: none.

  ALGORITHM: corner = gcd(2U_d - eps', (V_d - eps*pow(2,(m+1)//2, 2U_d-eps'))).
  Diagonal census at (d, p) = four modpows; no factorization of V_{2d}.

OUTPUT: A. cross-check vs direct factorization (small d).
        B. recover the d=57, p=11, r=683 danger triple via corners.
        C. sweep admissible odd d in [45,121] x cond-v primes p <= P_MAX;
           report all corner hits and any diagonal pair.
        D. resonance-band minimum of |V_d - 2^{(p+1)/2}|.
"""

import sys
from math import gcd, log2
from sympy import isprime, primerange, factorint

P_MAX = 1_000_000
D_LO, D_HI = 45, 121

def pell_VU(n):
    V0, V1 = 2, 2
    U0, U1 = 0, 1
    if n == 0:
        return V0, U0
    for _ in range(n - 1):
        V0, V1 = V1, 2 * V1 + V0
        U0, U1 = U1, 2 * U1 + U0
    return V1, U1

def in_class_A(q):
    # v2(ord_q(2)) = 1: ord = 2m, m odd, 2^m = -1
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
    """factor a corner using: every prime factor = 1 mod 2m."""
    if g == 1:
        return []
    if g < (2 * m + 1) ** 2 or isprime(g):
        return [(g, 1)]
    return sorted(factorint(g).items())

print("=" * 100, flush=True)
print("A. cross-check corner localization vs direct factorization (odd d <= 35)")
mismatch = total = 0
for d in range(3, 37, 2):
    V_d, U_d = pell_VU(d)
    direct = {}
    for half in (2 * U_d - 1, 2 * U_d + 1):
        for r in factorint(half):
            if r <= 3:
                continue
            o = r - 1
            for q in factorint(r - 1):
                while o % q == 0 and pow(2, o // q, r) == 1:
                    o //= q
            if o % 2 == 0 and pow(2, o // 2, r) == r - 1 and (o // 2) % 2 == 1:
                direct[r] = o // 2
    for r, m in direct.items():
        total += 1
        cs = corners(V_d, U_d, m)
        hits = [key for key, g in cs.items() if g % r == 0]
        if len(hits) != 1:
            mismatch += 1
            print("   MISMATCH d=%d r=%d m=%d hits=%s" % (d, r, m, hits))
print("   checked %d (r, m) cases; mismatches: %d (expected 0)" % (total, mismatch), flush=True)

print("=" * 100, flush=True)
print("B. danger triple d=57, p=11 via corners (no factorization)")
V57, U57 = pell_VU(57)
for key, g in corners(V57, U57, 11).items():
    if g > 1:
        print("   corner (eps=%+d, eps'=%+d): odd part = %d, factors %s" %
              (key[0], key[1], g, classify_corner(g, 11)))
print("   contains 683:", any(g % 683 == 0 for g in corners(V57, U57, 11).values()), flush=True)

print("=" * 100, flush=True)
adm_ds = [d for d in range(D_LO, D_HI + 1, 2) if admissible(d)]
print("C. sweep: admissible odd d in [%d,%d] = %s" % (D_LO, D_HI, adm_ds))
print("   cond-v primes p <= %d; corner census" % P_MAX, flush=True)

primes = list(primerange(5, P_MAX))
VU = {d: pell_VU(d) for d in adm_ds}
prim_cache = {}

def primitive_at_d(r, d):
    key = (r, d)
    if key in prim_cache:
        return prim_cache[key]
    ok = True
    for e in range(1, d):
        if d % e == 0:
            Ve, _ = pell_VU(e)
            if (Ve * Ve + 2) % r == 0:
                ok = False
                break
    prim_cache[key] = ok
    return ok

total_dp = 0
hit_lines = []
qual_events = []
res_min = (10.0, None)

for d in adm_ds:
    V_d, U_d = VU[d]
    lv = log2(V_d)
    mod3d = 3 * d
    condv_primes = [p for p in primes if pow(2, p - 2, mod3d) == mod3d - 1]
    print("   d=%-3d cond-v primes <= %d: %d" % (d, P_MAX, len(condv_primes)), flush=True)
    for p in condv_primes:
        total_dp += 1
        cs = corners(V_d, U_d, p)
        for key, g in cs.items():
            if g > 1:
                for r, e in classify_corner(g, p):
                    m8 = r % 8
                    prim = primitive_at_d(r, d)
                    big = r > 10**12
                    qual = (m8 == 3) and prim
                    hit_lines.append("   d=%-3d p=%-6d corner(%+d,%+d): r=%d^%d mod8=%d prim=%s >1e12=%s QUAL=%s"
                                     % (d, p, key[0], key[1], r, e, m8, prim, big, qual))
                    if qual:
                        qual_events.append((d, p, r, big))
        if abs((p + 1) / 2 - lv) < 3:
            H = abs(V_d - (1 << ((p + 1) // 2)))
            if H > 0 and log2(H) / lv < res_min[0]:
                res_min = (log2(H) / lv, (d, p))

print("   (d,p) pairs tested: %d ; corner hits: %d" % (total_dp, len(hit_lines)), flush=True)
for line in hit_lines:
    print(line)
print("\n   QUALIFYING danger events (r=3 mod 8, primitive, d>43): %d" % len(qual_events))
bydp = {}
for (d, p, r, big) in qual_events:
    bydp.setdefault((d, p), []).append(r)
npair = sum(1 for v in bydp.values() if len(v) >= 2)
print("   DIAGONAL PAIRS (>=2 qualifying r at one (d,p)): %d" % npair)
for k, v in bydp.items():
    if len(v) >= 2:
        print("   PAIR!!", k, v)

print("\nD. resonance minimum log2|H|/log2(V_d) = %.4f at (d,p)=%s" % res_min)
print("done.")
