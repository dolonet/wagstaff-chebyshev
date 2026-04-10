#!/usr/bin/env python3
"""Census of known split Wagstaff factors and case classification.

Paper references: Section 5.3 (Cases A/B/4U), Section 8.2 (case breakdown).

Enumerates all primes r = 1 (mod 8) with r < R_LIMIT and ord_r(2) = 2p
for prime p <= P_LIMIT, with W_p composite. For each split factor:
  1. Classify: Case A (r | U_p), Case B (r | V_p), 4U (4U_p^2 = 1 mod r),
     or general Case C.
  2. Verify Conjecture R3: p | ord_r(alpha), alpha = 1 + sqrt(2).
  3. Compute d = ord_r(omega_3)/4 and check NCT coverage.

For Case C factors, d = ord_r(omega_3)/4 determines whether Condition II
could potentially hold. Condition II requires d odd and 4 | ord(omega_3).
Factors where this fails are automatically safe (no compatible triple).

Expected output for R_LIMIT = 10^7, P_LIMIT = 5000:
  185 split factors, R3 verified at all 185
  Case A: 1, Case B: 6, Case 4U: 46, Case C: 132

Requires: sympy (for factorint)
Usage:    python3 split_factor_census.py [R_LIMIT] [P_LIMIT]
"""
import sys
from sympy import factorint, isprime


# Known Wagstaff prime/PRP exponents (W_p prime, skip these).
WAGSTAFF_PRP = {
    3, 5, 7, 11, 13, 17, 19, 23, 31, 43, 61, 79, 101, 127,
    167, 191, 199, 313, 347, 701, 1709, 2617, 3539, 5807,
    10501, 10691, 11279, 12391, 14479, 42737, 83339, 95369,
    117239, 127031, 138937, 141079, 267017, 269987, 374321,
    986191, 4031399, 13347311,
}


def multiplicative_order(a, r):
    """Compute ord_r(a) by dividing down from r-1."""
    order = r - 1
    for p, e in factorint(order).items():
        for _ in range(e):
            if order % p == 0 and pow(a, order // p, r) == 1:
                order //= p
    return order


def U_mod(n, m):
    """Compute Pell U_n mod m via matrix exponentiation."""
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


def V_mod(n, m):
    """Compute Pell V_n mod m via trace of matrix power."""
    if m == 1:
        return 0
    if n == 0:
        return 2 % m
    def mm(A, B):
        return [[(A[0][0]*B[0][0]+A[0][1]*B[1][0]) % m,
                 (A[0][0]*B[0][1]+A[0][1]*B[1][1]) % m],
                [(A[1][0]*B[0][0]+A[1][1]*B[1][0]) % m,
                 (A[1][0]*B[0][1]+A[1][1]*B[1][1]) % m]]
    R = [[1, 0], [0, 1]]
    B = [[2, 1], [1, 0]]
    e = n
    while e > 0:
        if e & 1:
            R = mm(R, B)
        B = mm(B, B)
        e >>= 1
    return (R[0][0] + R[1][1]) % m


def mul_z2(a, b, r):
    return ((a[0]*b[0] + 2*a[1]*b[1]) % r, (a[0]*b[1] + a[1]*b[0]) % r)


def pow_z2(base, exp, r):
    result = (1, 0)
    b = (base[0] % r, base[1] % r)
    while exp > 0:
        if exp & 1:
            result = mul_z2(result, b, r)
        b = mul_z2(b, b, r)
        exp >>= 1
    return result


def order_z2(base, r, group_order):
    order = group_order
    for p, e in factorint(group_order).items():
        for _ in range(e):
            if order % p == 0 and pow_z2(base, order // p, r) == (1, 0):
                order //= p
    return order


def main():
    r_limit = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10**7
    p_limit = int(float(sys.argv[2])) if len(sys.argv) > 2 else 5000

    print(f"Split Wagstaff factor census (r < {r_limit}, p <= {p_limit})")
    print("=" * 60)

    factors = []  # (p, r)

    for r in range(17, r_limit, 8):  # r = 1 mod 8 starts at 17
        if not isprime(r):
            continue
        o = multiplicative_order(2, r)
        if o % 2 != 0:
            continue
        p = o // 2
        if p < 5 or p > p_limit or not isprime(p):
            continue
        if p in WAGSTAFF_PRP:
            continue
        # Verify r | W_p: 2^p = -1 mod r
        if pow(2, p, r) != r - 1:
            continue
        factors.append((p, r))

    total = len(factors)
    print(f"\nTotal split factors: {total}")

    case_A = 0
    case_B = 0
    case_4U = 0
    case_C = 0
    case_C_auto = 0       # d inadmissible (4 nmid k or d even)
    case_C_nct81 = 0      # d odd, d <= 81
    case_C_nct200 = 0     # d odd, 83 <= d <= 200
    case_C_open = 0       # d odd, d > 200
    r3_fail = 0

    for p, r in factors:
        up = U_mod(p, r)
        vp = V_mod(p, r)
        four_up2 = (4 * up * up) % r

        if up == 0:
            case_A += 1
        elif vp == 0:
            case_B += 1
        elif four_up2 == 1:
            case_4U += 1
        else:
            case_C += 1
            # d = ord_r(omega_3)/4; Condition II requires d odd, 4|k
            k = order_z2((3, 2), r, r - 1)
            if k % 4 != 0 or (k // 4) % 2 == 0:
                case_C_auto += 1
            else:
                d = k // 4
                if d <= 81:
                    case_C_nct81 += 1
                elif d <= 200:
                    case_C_nct200 += 1
                else:
                    case_C_open += 1

        # R3: p | ord_r(alpha), alpha = (1,1) in Z[sqrt(2)]
        ord_alpha = order_z2((1, 1), r, r - 1)
        if ord_alpha % p != 0:
            r3_fail += 1
            print(f"  R3 FAIL: p={p}, r={r}")

    print(f"\nCase breakdown:")
    print(f"  Case A (r | U_p):      {case_A:4d} ({100*case_A/total:.0f}%)")
    print(f"  Case B (r | V_p):      {case_B:4d} ({100*case_B/total:.0f}%)")
    print(f"  Case 4U (4U_p^2 = 1):  {case_4U:4d} ({100*case_4U/total:.0f}%)")
    print(f"  Case C (general):      {case_C:4d} ({100*case_C/total:.0f}%)")

    uncond = case_A + case_B + case_4U + case_C_auto + case_C_nct81
    print(f"\nCase C detail (d = ord(omega_3)/4):")
    print(f"  d inadmissible (auto-safe): {case_C_auto}")
    print(f"  d odd, d <= 81 (NCT-81):    {case_C_nct81}")
    print(f"  d odd, 83 <= d <= 200:      {case_C_nct200}")
    print(f"  d odd, d > 200:             {case_C_open}")

    print(f"\nR3 verification: {total - r3_fail}/{total} passed")
    print(f"Unconditionally covered: {uncond}/{total}")

    if r3_fail == 0:
        print(f"\n{'=' * 60}")
        print(f"All {total} split factors: R3 verified.")
        if case_C_open + case_C_nct200 > 0:
            print(f"Note: {case_C_open + case_C_nct200} Case C factors have d > 81;")
            print(f"these are covered by R3 computational verification, not NCT-81.")
        print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
