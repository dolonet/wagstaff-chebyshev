#!/usr/bin/env python3
"""NCT verification for d in [83, 200].

Paper reference: Section 8.4.

For each parity-unblocked d in {83, 99, 107, 121, 129, 131, 139, 163,
171, 177, 179}, searches primes r = 8dk + 1 with r < R_LIMIT for
primitive divisors of U_{4d}, then checks whether any compatible triple
(d, r, p) exists: ord_r(2) = 2p with p prime and d | W_{p-2}.

A compatible triple would refute Conjecture R3 at that factor; finding
zero compatible triples extends NCT beyond the unconditional d <= 81.

The paper's full computation used R_LIMIT = 10^11 on a 128-core machine.
Default here is 10^9 (scales to full range with more time).

Requires: sympy (for factorint, isprime)
Usage:    python3 nct_extended_verify.py [R_LIMIT]
"""
import sys
import time
from sympy import factorint, isprime


def U_mod(n, m):
    """Pell U_n mod m via matrix exponentiation."""
    if m == 1 or n <= 0:
        return 0
    def mm(A, B):
        return [[(A[0][0]*B[0][0]+A[0][1]*B[1][0]) % m,
                 (A[0][0]*B[0][1]+A[0][1]*B[1][1]) % m],
                [(A[1][0]*B[0][0]+A[1][1]*B[1][0]) % m,
                 (A[1][0]*B[0][1]+A[1][1]*B[1][1]) % m]]
    R = [[1, 0], [0, 1]]
    B = [[2, 1], [1, 0]]
    e = n - 1
    while e > 0:
        if e & 1:
            R = mm(R, B)
        B = mm(B, B)
        e >>= 1
    return R[0][0] % m


def multiplicative_order(a, r):
    """Compute ord_r(a) by dividing down from r-1."""
    order = r - 1
    for p, e in factorint(order).items():
        for _ in range(e):
            if order % p == 0 and pow(a, order // p, r) == 1:
                order //= p
    return order


# The 11 parity-unblocked odd d values in [83, 200].
PARITY_UNBLOCKED = [83, 99, 107, 121, 129, 131, 139, 163, 171, 177, 179]


def main():
    r_limit = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10**9

    print(f"NCT extended verification (d in [83, 200], r < {r_limit})")
    print("=" * 60)

    total_primes = 0
    total_primitives = 0
    total_compatible = 0

    for d in PARITY_UNBLOCKED:
        t0 = time.time()
        step = 8 * d
        primitives = 0
        compatible = 0
        primes_checked = 0

        r = step + 1
        while r < r_limit:
            if isprime(r):
                primes_checked += 1
                # Check primitivity: rho(r) = 4d
                # r | U_{4d} is necessary
                if U_mod(4 * d, r) == 0:
                    # Verify it's primitive: r does not divide U_n
                    # for any proper divisor n of 4d
                    is_primitive = True
                    for q, e in factorint(4 * d).items():
                        if U_mod(4 * d // q, r) == 0:
                            is_primitive = False
                            break
                    if is_primitive:
                        primitives += 1
                        # Check compatibility: ord_r(2) = 2p, p prime,
                        # d | W_{p-2}
                        o = multiplicative_order(2, r)
                        if o % 2 == 0:
                            p = o // 2
                            if p >= 5 and isprime(p):
                                        # d | W_{p-2} = (2^{p-2}+1)/3
                                val = (pow(2, p - 2, 3 * d) + 1) % (3 * d)
                                if val % 3 == 0 and (val // 3) % d == 0:
                                    compatible += 1
                                    print(f"  *** COMPATIBLE TRIPLE: "
                                          f"d={d}, r={r}, p={p}")
            r += step

        elapsed = time.time() - t0
        total_primes += primes_checked
        total_primitives += primitives
        total_compatible += compatible
        print(f"d={d:3d}: {primes_checked:7d} primes, "
              f"{primitives:2d} primitives, "
              f"{compatible} compatible  ({elapsed:.1f}s)")

    print(f"\n{'=' * 60}")
    print(f"Total: {total_primes} primes checked, "
          f"{total_primitives} primitive divisors, "
          f"{total_compatible} compatible triples")
    if total_compatible == 0:
        print(f"No compatible triple found for any d in [83, 200].")
    else:
        print(f"*** {total_compatible} COMPATIBLE TRIPLES FOUND ***")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
