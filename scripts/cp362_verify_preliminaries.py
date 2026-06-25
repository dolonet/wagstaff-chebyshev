#!/usr/bin/env python3
"""
cp360: Adversarial numerical verification of FOUNDATIONAL PRELIMINARIES of
wagstaff_chebyshev_reduction_v3.6.tex.

Per-factor structural checks (lem:ord2, lem:mod8, lem:v2-class, thm:condI,
thm:ord-sqrt2) need only KNOWN prime factors of W_p, not full factorization.
We collect prime factors via factorint with an ECM/limit budget and run the
per-factor checks on every prime factor found (whether or not the cofactor
splits fully).  Pure-arithmetic claims (W_p mod 8/16, fundamental, pndivW,
necessity) need no factoring.
"""
import sys
from sympy import factorint, isprime, n_order, primerange, sqrt_mod, gcd

def mulZ2(x, y, m):
    (a, b) = x; (c, d) = y
    return ((a*c + 2*b*d) % m, (a*d + b*c) % m)

def powZ2(x, e, m):
    res = (1 % m, 0); base = (x[0] % m, x[1] % m)
    while e > 0:
        if e & 1: res = mulZ2(res, base, m)
        base = mulZ2(base, base, m); e >>= 1
    return res

def ordZ2(x, m, group_order):
    o = group_order
    for q in factorint(group_order):
        while o % q == 0 and powZ2(x, o // q, m) == (1 % m, 0):
            o //= q
    return o

def W(p): return (2**p + 1)//3
OMEGA3 = (3, 2)
def v2(n): return (n & -n).bit_length() - 1

# collect prime factors of n with a bounded budget (don't insist on full factoring)
def known_prime_factors(n, budget=20000):
    out = set()
    fac = factorint(n, limit=budget)   # trial up to budget + pollard rho/p-1
    for r, e in fac.items():
        if isprime(r):
            out.add((r, e))
    return out

test_ps = list(primerange(5, 130))

print("=== Per-factor checks on KNOWN prime factors of composite W_p ===")
fails = []; nfac = 0; pchecked = []
for p in test_ps:
    n = W(p)
    if isprime(n):
        continue
    pchecked.append(p)
    for r, e in known_prime_factors(n):
        nfac += 1
        o2 = n_order(2, r)
        if o2 != 2*p: fails.append(f"p={p} r={r}: ord_r(2)={o2} != 2p")
        if r % (2*p) != 1: fails.append(f"p={p} r={r}: r%2p={r%(2*p)} !=1")
        if r % 8 not in (1,3): fails.append(f"p={p} r={r}: r%8={r%8} not in {{1,3}} !!")
        if r % 8 == 3 and r % 16 not in (3,11): fails.append(f"p={p} r={r}: r%16={r%16} (mod16b)")
        if pow(2,(n-1)//2,r) != r-1: fails.append(f"p={p} r={r}: condI fail")
        grp = r-1 if r%8 in (1,7) else r*r-1
        ow = ordZ2(OMEGA3, r, grp); vw = v2(ow)
        exp = {3:2,5:1,7:0}.get(r%8, None)
        if exp is not None and vw != exp: fails.append(f"p={p} r={r}(={r%8}mod8): v2(ord w3)={vw}!={exp}")
        if r%8==1 and vw > v2(r-1)-1: fails.append(f"p={p} r={r}(1mod8): v2={vw}>v2(r-1)-1")
        if r%8==1:
            s = sqrt_mod(2,r); os = n_order(s,r)
            if os != 4*p: fails.append(f"p={p} r={r}(1mod8): ord(sqrt2)={os}!=4p")
print(f"  composite W_p covered: {pchecked}")
print(f"  prime factors checked: {nfac}")
print(f"  FAILURES: {len(fails)}")
for f in fails[:60]: print("   !!", f)

print("\n=== W_p mod 8 / mod 16 + parity (FULL factorization for parity where feasible) ===")
af = []
for p in test_ps:
    n = W(p)
    if n % 8 != 3: af.append(f"W_{p}%8={n%8}!=3")
    if n % 16 != 11: af.append(f"W_{p}%16={n%16}!=11")
print(f"  W_p mod 8/16 over p in [5,130): FAILURES={len(af)}")
for f in af[:10]: print("   !!", f)

# parity of 3mod8 multiplicity -- need full factorization; do it for fully-factorable p only
print("\n  parity of 3mod8-factor multiplicity (fully factored W_p only):")
parity_fail = []; pf_done = []
for p in test_ps:
    n = W(p)
    if isprime(n):
        pf_done.append((p,'prime, 1 factor'))
        # prime W_p = 3 mod 8 itself: single factor =3mod8, mult 1 odd OK
        if n % 8 == 3 and 1 % 2 != 1: parity_fail.append(p)
        continue
    fac = factorint(n)  # full
    if all(isprime(r) for r in fac):
        mult3 = sum(e for r,e in fac.items() if r%8==3)
        pf_done.append((p, f'{mult3} (3mod8 mult)'))
        if mult3 % 2 != 1: parity_fail.append(f"p={p}: 3mod8 mult={mult3} even !!")
print(f"    fully resolved p: {[x[0] for x in pf_done]}")
print(f"    parity FAILURES: {len(parity_fail)}", parity_fail if parity_fail else '')

print("\n=== lem:fundamental + lem:pndivW (p in [5,200)) ===")
ff = []
for p in primerange(5, 200):
    n = W(p); nm2 = W(p-2)
    if n+1 != 4*nm2: ff.append(f"p={p}: W_p+1!=4W_{{p-2}}")
    if gcd(n, nm2) != 1: ff.append(f"p={p}: gcd!=1")
    if nm2 % p == 0: ff.append(f"p={p}: p|W_{{p-2}} !!")
print(f"  FAILURES: {len(ff)}", ff if ff else '')

print("\n=== prop:necessity (prime W_p) + discriminating-on-composite sanity ===")
nf = []; primeWp = []
for p in primerange(5, 600):
    n = W(p)
    if isprime(n):
        primeWp.append(p)
        if powZ2(OMEGA3, (n+1)//2, n) != ((-1)%n, 0):
            nf.append(f"p={p}: necessity fails on prime W_p")
print(f"  prime W_p (p): {primeWp}")
print(f"  necessity FAILURES on prime W_p: {len(nf)}", nf if nf else '')
comp_pass = []
for p in primerange(5, 200):
    n = W(p)
    if not isprime(n) and n>1 and powZ2(OMEGA3, (n+1)//2, n) == ((-1)%n, 0):
        comp_pass.append(p)
print(f"  composite W_p passing Cond II globally (expect none): {comp_pass if comp_pass else 'NONE'}")

print("\n=== mod-8 dichotomy brute: any r=5/7 mod8 | W_p, p in [5,200) (KNOWN factors) ===")
viol = []
for p in primerange(5, 200):
    for r,e in known_prime_factors(W(p)):
        if r%8 in (5,7): viol.append((p,r,r%8))
print(f"  r=5or7 mod8 factors found: {viol if viol else 'NONE'}")

tot = len(fails)+len(af)+len(parity_fail)+len(ff)+len(nf)
print(f"\n{'ALL PRELIMINARY CHECKS PASSED' if tot==0 else 'TOTAL FAILURES: '+str(tot)}")
