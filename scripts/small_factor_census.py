#!/usr/bin/env python3
"""Census of small inert prime factors of composite Wagstaff numbers.

Paper references: Section 8.5 (Table 1), Section 6.4 (Corollary small-Gr).

Enumerates all primes r = 3 (mod 8) with r < LIMIT (default 200000),
computes ord_r(2), identifies those where ord_r(2) = 2p for prime p >= 5
with W_p composite, and computes G_r/4 = gcd((r+1)/4, W_{p-2}).

Expected output for LIMIT = 200000 (the paper's census):
  - 869 inert factors of composite W_p
  - 11 distinct G_r/4 values: {1, 3, 9, 11, 19, 57, 129, 177, 321, 2217, 2451}
  - Prime factors of all G_r/4 values: {3, 11, 19, 43, 59, 107, 739}
  - 845/869 (97.2%) have G_r/4 <= 43 (Corollary small-Gr)
  - 24/869 (2.8%) have G_r/4 > 43 (verified by order excess)

Requires: sympy (for isprime, factorint)
Usage:    python3 small_factor_census.py [LIMIT]
"""
import sys
from math import gcd
from sympy import isprime, factorint


# Known Wagstaff prime/PRP exponents — W_p is prime or probable prime.
# From survey.py; matches the paper's list of 36 known Wagstaff primes/PRPs.
WAGSTAFF_PRIME_EXPONENTS = {
    3, 5, 7, 11, 13, 17, 19, 23, 31, 43, 61, 79, 101, 127,
    167, 191, 199, 313, 347, 701, 1709, 2617, 3539, 5807,
    10501, 10691, 11279, 12391, 14479, 42737, 83339, 95369,
    117239, 127031, 138937, 141079, 267017, 269987, 374321,
    986191, 4031399, 13347311,
}


def multiplicative_order(a, m):
    """Return ord_m(a) by dividing down from the group order.

    For prime m, ord_m(a) | (m-1). Factor m-1 and remove prime
    factors while the power remains 1. Much faster than trial.
    """
    order = m - 1
    for p, e in factorint(order).items():
        for _ in range(e):
            if pow(a, order // p, m) == 1:
                order //= p
            else:
                break
    return order


def main():
    limit = int(sys.argv[1]) if len(sys.argv) > 1 else 200000

    print(f"Small inert factor census (r < {limit})")
    print("=" * 60)

    G4_values = []
    rows = []

    for r in range(11, limit + 1):
        if not isprime(r) or r % 8 != 3:
            continue
        o = multiplicative_order(2, r)
        if o % 2 != 0:
            continue
        p = o // 2
        if not isprime(p) or p < 5:
            continue
        # Exclude known Wagstaff primes (W_p must be composite).
        # For r < 200000, the Wagstaff prime exponents in range are
        # all listed below. We skip isprime(W_p) for large p since
        # W_p has thousands of digits.
        if p in WAGSTAFF_PRIME_EXPONENTS:
            continue

        # G_r/4 = gcd((r+1)/4, W_{p-2})
        # Use modular arithmetic: compute W_{p-2} mod (r+1)/4.
        q = (r + 1) // 4
        W_p2_mod_q = (pow(2, p - 2, 3 * q) + 1) // 3
        G4 = gcd(q, W_p2_mod_q)
        G4_values.append(G4)
        rows.append((p, r, G4))

    total = len(G4_values)
    unique = sorted(set(G4_values))

    print(f"\nTotal inert factors of composite W_p: {total}")
    print(f"Distinct G_r/4 values ({len(unique)}): {unique}")

    # Prime factors of all G_r/4 values
    all_prime_factors = set()
    for v in unique:
        if v > 1:
            all_prime_factors.update(factorint(v).keys())
    print(f"Prime factors of all G_r/4: {sorted(all_prime_factors)}")

    # Distribution
    print(f"\nDistribution:")
    for threshold in [1, 3, 43]:
        n = sum(1 for g in G4_values if g <= threshold)
        print(f"  G_r/4 <= {threshold:3d}: {n:4d}/{total} ({100*n/total:.1f}%)")
    n_large = sum(1 for g in G4_values if g > 43)
    print(f"  G_r/4 >  43: {n_large:4d}/{total} ({100*n_large/total:.1f}%)")

    # Verification
    ok = True
    if total != 869 and limit == 200000:
        print(f"\n*** UNEXPECTED: expected 869 factors, got {total}")
        ok = False
    if len(unique) != 11 and limit == 200000:
        print(f"\n*** UNEXPECTED: expected 11 distinct G_r/4 values, got {len(unique)}")
        ok = False

    if ok and limit == 200000:
        print(f"\n{'=' * 60}")
        print(f"Census matches paper: 869 factors, 11 distinct G_r/4 values.")
        print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
