#!/usr/bin/env python3
"""Wagstaff Chebyshev primality test (paper, section 8.6).

Usage: python3 chebyshev_test.py <p>

Tests whether W_p = (2^p + 1) / 3 is a probable prime using the
Chebyshev congruence omega_3^{(W_p+1)/2} = -1 (mod W_p).
"""
import sys


def chebyshev_test(p):
    """Return True if W_p passes the Chebyshev test (probable prime)."""
    N = (pow(2, p) + 1) // 3

    # Work in Z[sqrt(2)]/(N), elements are pairs (a, b) for a + b*sqrt(2).
    def mul(x, y):
        a, b = x
        c, d = y
        return ((a * c + 2 * b * d) % N, (a * d + b * c) % N)

    # Compute omega_3^{(N+1)/2} mod N by binary exponentiation.
    exp = (N + 1) // 2
    result = (1, 0)   # identity
    base = (3, 2)     # omega_3 = 3 + 2*sqrt(2)
    while exp > 0:
        if exp & 1:
            result = mul(result, base)
        base = mul(base, base)
        exp >>= 1

    return result == (N - 1, 0)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 chebyshev_test.py <p>")
        print("Tests whether W_p = (2^p+1)/3 is a probable prime.")
        sys.exit(1)

    p = int(sys.argv[1])
    if p < 5 or p % 2 == 0:
        print(f"Error: p must be an odd prime >= 5, got {p}")
        sys.exit(1)

    if chebyshev_test(p):
        print(f"W_{p} is a PROBABLE PRIME (Chebyshev test passed)")
    else:
        print(f"W_{p} is COMPOSITE (Chebyshev test failed)")
