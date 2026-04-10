#!/usr/bin/env python3
"""Verify NCT (No Compatible Triples) for d <= 81 and Pell exclusion for d <= 43.

Paper references: Theorem 5.12 (NCT for d <= 81), Theorem 6.5 (Pell exclusion).

For each admissible odd d:
  1. Compute U_{4d} and factor it.
  2. Find primitive divisors r = 1 (mod 8) of U_{4d}.
  3. Check that no compatible triple (d, r, p) exists.
  4. (For d <= 43) Compute N_d = gcd((V_{4d}+2)/2, U_{4d}) and verify
     that no inert prime r = 3 (mod 8) dividing N_d can yield a valid
     configuration.

Requires: sympy (for factorint)
Usage:    python3 nct_verify.py
"""
import sys
from sympy import factorint, isprime


def pell_UV(n):
    """Compute (U_n, V_n) for the Pell sequence: alpha = 1 + sqrt(2)."""
    if n == 0:
        return (0, 2)
    U_prev, V_prev = 0, 2
    U_curr, V_curr = 1, 2
    for _ in range(n - 1):
        U_prev, U_curr = U_curr, 2 * U_curr + U_prev
        V_prev, V_curr = V_curr, 2 * V_curr + V_prev
    return (U_curr, V_curr)


def ord2(n):
    """Return v_2(n), the 2-adic valuation of n."""
    if n == 0:
        return float("inf")
    v = 0
    while n % 2 == 0:
        v += 1
        n //= 2
    return v


def multiplicative_order(a, m):
    """Return ord_m(a), assuming gcd(a, m) = 1."""
    if m == 1:
        return 1
    r = 1
    x = a % m
    while True:
        x = (x * a) % m
        r += 1
        if x == 1:
            return r


def is_class_iii(q):
    """Check if prime q has v_2(ord_q(2)) = 1 (Class III)."""
    o = multiplicative_order(2, q)
    return ord2(o) == 1


def primitive_divisors_mod8_1(d):
    """Find primitive divisors r = 1 (mod 8) of U_{4d}."""
    U4d, _ = pell_UV(4 * d)
    factors = factorint(abs(U4d))
    result = []
    for r, _ in factors.items():
        if r < 3 or r % 8 != 1:
            continue
        # Check that r is a *primitive* divisor: rho(r) = 4d.
        # rho(r) = min n >= 1 with r | U_n.
        # Since r | U_{4d}, we need r does not divide U_n for any proper
        # divisor n of 4d.
        is_primitive = True
        n = 4 * d
        for q, e in factorint(n).items():
            if pell_UV(n // q)[0] % r == 0:
                is_primitive = False
                break
        if is_primitive:
            result.append(r)
    return result


def check_nct(d, r):
    """Check whether a compatible triple (d, r, p) could exist.

    Returns (blocked, reason).
    """
    o = multiplicative_order(2, r)
    v2_o = ord2(o)
    if v2_o != 1:
        return True, f"v_2(ord_r(2)) = {v2_o} != 1"
    p_candidate = o // 2
    if not isprime(p_candidate):
        return True, f"ord_r(2)/2 = {p_candidate} is composite"
    # Check d | W_{p-2}
    Wp2 = (pow(2, p_candidate - 2, d * d) + 1) // 3
    if Wp2 % d != 0:
        return True, f"d does not divide W_{{p-2}} (p = {p_candidate})"
    return False, f"COMPATIBLE TRIPLE FOUND: d={d}, r={r}, p={p_candidate}"


def admissible_d_values(limit):
    """Return odd d in [1, limit] with all prime factors in Class III."""
    result = []
    for d in range(1, limit + 1, 2):
        if d == 1:
            result.append(d)
            continue
        factors = factorint(d)
        if all(is_class_iii(q) for q in factors):
            result.append(d)
    return result


def main():
    print("=" * 60)
    print("NCT verification for d <= 81 (Theorem 5.12)")
    print("Pell exclusion verification for d <= 43 (Theorem 6.5)")
    print("=" * 60)

    admissible = admissible_d_values(81)
    print(f"\nAdmissible d values in [1, 81]: {admissible}")
    print(f"Count: {len(admissible)}")

    # --- NCT for d <= 81 ---
    print("\n--- NCT (split branch) ---\n")
    total_prims = 0
    all_blocked = True

    for d in admissible:
        if d <= 3:
            U4d, _ = pell_UV(4 * d)
            factors = factorint(abs(U4d))
            has_r1mod8 = any(r % 8 == 1 and r > 2 for r in factors)
            if not has_r1mod8:
                print(f"d = {d:3d}: U_{4*d} = {U4d}, "
                      f"no prime r = 1 (mod 8) exists. BLOCKED.")
                continue

        prims = primitive_divisors_mod8_1(d)
        total_prims += len(prims)

        if not prims:
            print(f"d = {d:3d}: no primitive r = 1 (mod 8). BLOCKED.")
            continue

        for r in prims:
            blocked, reason = check_nct(d, r)
            status = "BLOCKED" if blocked else "*** FAIL ***"
            print(f"d = {d:3d}: r = {r}, {reason}. {status}")
            if not blocked:
                all_blocked = False

    print(f"\nTotal primitive divisors r = 1 (mod 8) checked: {total_prims}")

    # --- Pell exclusion for d <= 43 ---
    print("\n--- Pell exclusion (inert branch) ---\n")
    admissible_43 = [d for d in admissible if d <= 43]

    for d in admissible_43:
        U4d, V4d = pell_UV(4 * d)
        Nd = abs((V4d + 2) // 2)
        from math import gcd
        Nd = gcd(Nd, abs(U4d))
        factors = factorint(Nd)
        inert_primes = [r for r in factors if r % 8 == 3 and r > 3]
        print(f"d = {d:3d}: N_d = {Nd}, "
              f"inert primes r = 3 (mod 8), r > 3: {inert_primes or 'none'}",
              end="")
        blocked = True
        for r in inert_primes:
            o = multiplicative_order(2, r)
            if ord2(o) != 1:
                print(f" [r={r}: v_2(ord)={ord2(o)}, not Class III]", end="")
                continue
            p_cand = o // 2
            if not isprime(p_cand):
                print(f" [r={r}: ord/2={p_cand} composite]", end="")
                continue
            # r is inert and divides W_p for prime p; but is W_p composite?
            # W_p = (2^p+1)/3; if W_p = r then p is a Wagstaff prime exponent.
            Wp = (pow(2, p_cand) + 1) // 3
            if Wp == r:
                print(f" [r={r}: W_{p_cand}={r} is prime]", end="")
                continue
            blocked = False
            print(f" *** FAIL at r={r}, p={p_cand} ***", end="")
        if blocked:
            print(" BLOCKED.")
        else:
            all_blocked = False

    if all_blocked:
        print("\n" + "=" * 60)
        print("ALL CHECKS PASSED. No compatible triple exists for d <= 81.")
        print("Pell exclusion verified for all admissible d <= 43.")
        print("=" * 60)
    else:
        print("\n*** VERIFICATION FAILED ***")
        sys.exit(1)


if __name__ == "__main__":
    main()
