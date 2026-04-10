#!/usr/bin/env python3
"""Verify that no Class III prime q < 5000 is a Wieferich-at-minus-1 prime.

Paper reference: Section 6.2, lines 1371-1375.

A prime q is Class III if v_2(ord_q(2)) = 1 (i.e., ord_q(2) = 2m with m odd).
The paper claims: v_q(2^m + 1) = 1 for all 194 Class III primes q < 5000.

Equivalently: 2^m = -1 (mod q) but 2^m != -1 (mod q^2). This matters because
if q^2 | (2^m + 1), the LTE (Lifting the Exponent) lemma would give extra
q-divisibility of W_{p-2}, potentially allowing d to absorb a higher power of q.

Requires: sympy (for isprime)
Usage:    python3 class_iii_wieferich.py [LIMIT]
"""
import sys
from sympy import isprime


def multiplicative_order(a, m):
    """Return ord_m(a) by trial."""
    x, o = a % m, 1
    while x != 1:
        x = (x * a) % m
        o += 1
    return o


def v2(n):
    """2-adic valuation of n."""
    if n == 0:
        return float('inf')
    v = 0
    while n % 2 == 0:
        v += 1
        n //= 2
    return v


def main():
    limit = int(sys.argv[1]) if len(sys.argv) > 1 else 5000

    print(f"Class III Wieferich-at-minus-1 check (q < {limit})")
    print("=" * 60)

    class_iii = []
    wieferich = []

    for q in range(3, limit):
        if not isprime(q):
            continue
        o = multiplicative_order(2, q)
        if v2(o) != 1:
            continue
        # q is Class III: ord_q(2) = 2m, m odd
        class_iii.append(q)
        m = o // 2
        # Check: 2^m = -1 mod q (should always hold)
        assert pow(2, m, q) == q - 1, f"2^m != -1 mod q for q={q}"
        # Check: 2^m != -1 mod q^2 (Wieferich-at-minus-1 test)
        q2 = q * q
        val = pow(2, m, q2)
        if val == q2 - 1:
            wieferich.append(q)
            print(f"  WIEFERICH at q={q}: 2^{m} = -1 (mod {q}^2)")

    print(f"\nClass III primes < {limit}: {len(class_iii)}")
    print(f"Wieferich-at-minus-1 primes: {len(wieferich)}")

    if not wieferich:
        print(f"\n{'=' * 60}")
        print(f"No Wieferich-at-minus-1 prime among {len(class_iii)} "
              f"Class III primes < {limit}.")
        print(f"v_q(2^m + 1) = 1 for all {len(class_iii)} primes.")
        print(f"{'=' * 60}")
    else:
        print(f"\n*** {len(wieferich)} Wieferich-at-minus-1 primes found ***")


if __name__ == "__main__":
    main()
