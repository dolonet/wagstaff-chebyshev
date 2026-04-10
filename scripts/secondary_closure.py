#!/usr/bin/env python3
"""Verify the secondary-factor closures for d=57 and d=67.

Paper reference: Section 8.8, lines 1980-1999.

The three local passes of Condition II in the extended survey occur at
d=57 and at two structured cases. For each, the paper finds a SECOND
inert prime factor r_2 of the same composite W_p and shows that
Condition II fails at r_2 (so the global pass is blocked).

d=57 closure:
  p = 85684865933, r_1 = 171369731867 (Sophie Germain: r_1 = 2p+1)
  r_2 = 7,675,646,348,774,307,083 = 2pk+1 with k = 44,789,977
  r_2 = 3 (mod 8) (inert), ord_{r_2}(omega_3) = r_2 + 1
  (r_2+1)/4 has prime factor 487; ord_487(2) = 243 is odd (Class I)
  Therefore d_{r_2} is not Class-III-smooth, Cond II fails at r_2.

d=67 closure:
  p = 24822855876938710967
  r_1 = 148937135261632265803
  r_2 = 5,474,333,343,676,555,561,818,313 = 2pk+1 with k = 110268
  r_2 = 1 (mod 8) (split), and R3 holds at r_2.

Requires: sympy
Usage:    python3 secondary_closure.py
"""
from sympy import isprime, factorint
from math import gcd


def v2(n):
    if n == 0:
        return float('inf')
    v = 0
    while n % 2 == 0:
        v += 1
        n //= 2
    return v


def multiplicative_order_fast(a, m):
    """Compute ord_m(a) by dividing down from m-1."""
    order = m - 1
    for p, e in factorint(order).items():
        for _ in range(e):
            if order % p == 0 and pow(a, order // p, m) == 1:
                order //= p
    return order


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


def classify_prime(q):
    """Classify q by v_2(ord_q(2))."""
    o = multiplicative_order_fast(2, q)
    v = v2(o)
    if v == 0:
        return "I", o
    elif v == 1:
        return "III", o
    else:
        return "II", o


def verify_factor(label, p, r):
    """Verify that r divides W_p and characterize."""
    print(f"\n  {label}: r = {r}")
    print(f"    isprime(r) = {isprime(r)}")
    print(f"    r mod 8 = {r % 8}", end="")
    if r % 8 == 1:
        print(" (split)")
    elif r % 8 == 3:
        print(" (inert)")
    else:
        print(" (excluded)")
    # Check r | W_p via 2^p = -1 mod r
    check = pow(2, p, r)
    print(f"    2^p mod r = {check} (expect r-1 = {r-1}): "
          f"{'r | W_p' if check == r-1 else 'r DOES NOT divide W_p'}")
    # Check Sophie Germain
    if r == 2*p + 1:
        print(f"    Sophie Germain: r = 2p+1")
    else:
        k = (r - 1) // (2 * p)
        print(f"    r = 2pk+1 with k = {k}")
    return check == r - 1


def analyze_inert_closure(label, p, r2):
    """For inert r2, check that Cond II fails via Class I/II obstruction."""
    print(f"\n  Analyzing Cond II at {label}...")
    # For inert r, group order is r+1
    rp1 = r2 + 1
    q4 = rp1 // 4
    print(f"    (r+1)/4 = {q4}")
    facs = factorint(q4, limit=10**7)
    print(f"    Partial factorization of (r+1)/4: {dict(facs)}")
    # Check each prime factor for Class III status
    for q in sorted(facs.keys()):
        if q <= 3:
            continue
        cls, o = classify_prime(q)
        print(f"    q={q}: ord_q(2)={o}, v_2={v2(o)}, Class {cls}", end="")
        if cls != "III":
            print(f" -> OBSTRUCTION (not Class III)")
            return True
        else:
            print()
    return False


def main():
    print("Secondary-factor closures for d=57 and d=67")
    print("=" * 60)

    # === d=57 closure ===
    print("\n--- d=57 ---")
    p57 = 85684865933
    r57_1 = 171369731867
    r57_2 = 7675646348774307083

    print(f"  p = {p57}")
    ok1 = verify_factor("r_1 (danger factor)", p57, r57_1)
    ok2 = verify_factor("r_2 (closure factor)", p57, r57_2)

    if ok1 and ok2:
        blocked = analyze_inert_closure("r_2", p57, r57_2)
        if blocked:
            print(f"\n  d=57: CLOSED. Cond II fails at r_2.")
        else:
            print(f"\n  d=57: *** NOT CLOSED ***")

    # === d=67 closure ===
    print("\n--- d=67 ---")
    p67 = 24822855876938710967
    r67_1 = 148937135261632265803
    r67_2 = 5474333343676555561818313

    print(f"  p = {p67}")
    ok1 = verify_factor("r_1 (danger factor)", p67, r67_1)
    ok2 = verify_factor("r_2 (closure factor)", p67, r67_2)

    if ok1 and ok2:
        # r_2 is split (r_2 = 1 mod 8), so R3 applies
        print(f"\n  r_2 = 1 mod 8 (split). Checking R3...")
        # ord_{r_2}(alpha) where alpha = 1+sqrt(2)
        # For split r, compute in F_r
        rm1 = r67_2 - 1
        rm1_facs = factorint(rm1, limit=10**9)
        print(f"    Partial factorization of r_2-1: {dict(rm1_facs)}")
        # Check if p | ord_{r_2}(alpha) via alpha^{(r-1)/p}
        # This is equivalent to checking that alpha is NOT a p-th power
        if rm1 % p67 == 0:
            exp = rm1 // p67
            val = pow_z2((1, 1), exp, r67_2)
            r3_holds = val != (1, 0)
            print(f"    alpha^((r-1)/p) = {val} "
                  f"({'!= 1, R3 holds' if r3_holds else '= 1, R3 FAILS'})")
            if r3_holds:
                print(f"\n  d=67: CLOSED. R3 blocks at r_2 (split factor).")
            else:
                print(f"\n  d=67: *** NOT CLOSED ***")
        else:
            print(f"    p does not divide r_2-1; checking full order...")
            print(f"    (r_2-1) mod p = {rm1 % p67}")


if __name__ == "__main__":
    main()
