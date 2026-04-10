#!/usr/bin/env python3
"""Primitive divisor survey and QPRO analysis.

Paper reference: Section 8.3.

For each prime q <= Q_LIMIT with q = 3 (mod 8), finds primitive divisors
r = 1 (mod 8) of U_{4q} by factoring V_{2q} = V_q^2 + 2. A primitive
divisor r satisfies rho(r) = 4q (i.e., r | U_{4q} but r does not divide
U_n for any proper divisor n of 4q).

For each primitive divisor, checks whether 2 is a qth power residue
modulo r: 2^{(r-1)/q} = 1 (mod r).

Expected result for q <= 1000 (paper's survey):
  62 primitive divisors among small q (where V_{2q} can be factored).
  Exactly one qth-power case: q=43, r=14449.
  QPRO counterexample: q=163, r=12036824328615957881 (found by ECM;
  not reproduced here due to size).

Requires: sympy (for factorint, isprime)
Usage:    python3 primitive_divisor_survey.py [Q_LIMIT]
"""
import sys
import time
from sympy import factorint, isprime


def V_pell(n):
    """Compute V_n (companion Pell) exactly. V_0=2, V_1=2, V_n = 2V_{n-1} + V_{n-2}."""
    if n == 0:
        return 2
    if n == 1:
        return 2
    a, b = 2, 2
    for _ in range(n - 1):
        a, b = b, 2 * b + a
    return b


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


def v2(n):
    """2-adic valuation."""
    if n == 0:
        return float('inf')
    v = 0
    while n % 2 == 0:
        v += 1
        n //= 2
    return v


def multiplicative_order(a, r):
    """Compute ord_r(a) by dividing down from r-1."""
    order = r - 1
    for p, e in factorint(order).items():
        for _ in range(e):
            if order % p == 0 and pow(a, order // p, r) == 1:
                order //= p
    return order


def main():
    q_limit = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    factor_limit = int(sys.argv[2]) if len(sys.argv) > 2 else 0

    print(f"Primitive divisor survey (q <= {q_limit}, q = 3 mod 8)")
    print("=" * 60)

    q_values = [q for q in range(3, q_limit + 1)
                if isprime(q) and q % 8 == 3]

    total_primitives = 0
    qth_power_cases = []

    for q in q_values:
        t0 = time.time()
        Vq = V_pell(q)
        V2q = Vq * Vq + 2

        digits = len(str(V2q))
        if digits > 60 and factor_limit == 0:
            # Skip: V_{2q} too large for sympy to factor
            print(f"q={q:4d}: V_{{2q}} has {digits} digits, skipping "
                  f"(use factor_limit>0 for partial)")
            continue

        try:
            if factor_limit > 0:
                factors = factorint(V2q, limit=factor_limit)
            else:
                factors = factorint(V2q)
        except Exception:
            print(f"q={q:4d}: factoring failed ({digits} digits)")
            continue

        primitives = []
        for r in factors:
            if r < 3 or not isprime(r) or r % 8 != 1:
                continue
            # Check primitivity: rho(r) = 4q
            # Need: r | U_{4q}, r nmid U_{2q}, r nmid U_4
            if U_mod(4 * q, r) != 0:
                continue
            if U_mod(2 * q, r) == 0 or U_mod(4, r) == 0:
                continue

            # Primitive divisor found
            o = multiplicative_order(2, r)
            is_qth = pow(2, (r - 1) // q, r) == 1
            primitives.append((r, o, is_qth))
            if is_qth:
                qth_power_cases.append((q, r, o))

        elapsed = time.time() - t0
        n_prim = len(primitives)
        total_primitives += n_prim
        if n_prim > 0:
            for r, o, is_q in primitives:
                tag = " [QTH POWER]" if is_q else ""
                print(f"q={q:4d}: r={r}, ord_r(2)={o}, "
                      f"v2(ord)={v2(o)}{tag}  ({elapsed:.1f}s)")
        else:
            if digits <= 60 or factor_limit > 0:
                print(f"q={q:4d}: 0 primitives (V_{{2q}}: {digits}d, "
                      f"{elapsed:.1f}s)")

    print(f"\n{'=' * 60}")
    print(f"Total primitive divisors (r = 1 mod 8): {total_primitives}")
    print(f"Q-th power cases: {len(qth_power_cases)}")
    for q, r, o in qth_power_cases:
        print(f"  q={q}, r={r}: 2^{{(r-1)/{q}}} = 1 mod r, "
              f"ord_r(2)={o}, v2(ord)={v2(o)}")
        half = o // 2
        if half > 1 and not isprime(half):
            print(f"    ord/2 = {half} is composite -> NCT unaffected")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
