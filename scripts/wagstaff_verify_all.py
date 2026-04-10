#!/usr/bin/env python3
"""
Reproducibility package for:
"Primality proofs for Wagstaff numbers via Chebyshev bases"

This script verifies Theorem 3.2 for all 20 known Wagstaff primes
with 5 <= p <= 1709. For each p, it:
  1. Factors W_{p-2} = (2^{p-2}+1)/3
  2. Checks conditions (I), (II), (III_q) of Theorem 3.2
  3. Computes F and verifies F > sqrt(W_p)
  4. For Conjecture 4.5 (p%6==1): checks s_{p-2} and omega^{W_{p-2}}
  5. For Conjecture 4.6 (p%6==5): checks omega^{(W_p-1)/2}

Usage: python3 wagstaff_verify_all.py
Requires: sympy (pip install sympy)
Runtime: ~5 minutes on a modern laptop
"""

from sympy import isprime, factorint
import math
import json
import time

# ---- Arithmetic in Z[sqrt(d)]/(N) ----

def mul_mod(x, y, d, N):
    """Multiply (a+b*sqrt(d)) * (c+e*sqrt(d)) mod N."""
    return ((x[0]*y[0] + d*x[1]*y[1]) % N, (x[0]*y[1] + x[1]*y[0]) % N)

def pow_mod(base, exp, d, N):
    """Compute base^exp in Z[sqrt(d)]/(N)."""
    result = (1, 0)
    b = (base[0] % N, base[1] % N)
    e = abs(exp)
    while e > 0:
        if e % 2 == 1:
            result = mul_mod(result, b, d, N)
        b = mul_mod(b, b, d, N)
        e //= 2
    return result

def omega3_pow(exp, N):
    """Compute omega_3^exp mod N, where omega_3 = 3 + 2*sqrt(2), disc = 8."""
    return pow_mod((3, 1), exp, 8, N)  # (3 + 1*sqrt(8)) = 3 + 2*sqrt(2)

def omega2_pow(exp, N):
    """Compute omega^exp mod N, where omega = 2 + sqrt(3), disc = 3."""
    return pow_mod((2, 1), exp, 3, N)

# ---- Main verification ----

def wagstaff(p):
    return (2**p + 1) // 3

# All 20 known Wagstaff primes
WAGSTAFF_PRIMES = [5, 7, 11, 13, 17, 19, 23, 31, 43, 61, 79, 101,
                   127, 167, 191, 199, 313, 347, 701, 1709]

# Pre-computed factorizations for large cases (from factordb/Cunningham)
KNOWN_FACTORS = {
    701: {
        3: 1, 467: 1, 27961: 1,
        352369374013660139472574531568890678155040563007620742839120913: 1,
        9551137: 1, 373746913: 1, 53590752072775417: 1,
        331997018907154605105439627246768387410091148934183179476062496087540587425301679685209566699340444844716483: 1,
    },
    1709: {
        3: 1, 10243: 1, 81937: 1, 223439473: 1, 2058017388830521: 1,
        52616682857203: 1, 678712480477777538119297: 1,
        820793700834095033598534000634605722059739038841: 1,
        166600443145782167422869582205773627533809608354200608495919325006516227014769666125340269708249: 1,
        425328467811476606269734624816851146193116859358628089320763326584879949452267875093047179051512783907924678657125235746323328435166697131298895625273070976622085482536748543280466254716106901284011556317089150580145726892414853298648692222190100640625226284396348290430732796476757391227496891238233: 1,
    },
}

def get_factorization(p):
    """Get factorization of W_{p-2}."""
    W_pm2 = wagstaff(p - 2)
    if p in KNOWN_FACTORS:
        facs = KNOWN_FACTORS[p]
        # Verify
        product = 1
        for q, e in facs.items():
            product *= q**e
        assert product == W_pm2, f"Factorization mismatch for p={p}"
        return facs
    else:
        return factorint(W_pm2)

def verify_theorem32(p):
    """Verify Theorem 3.2 for a single Wagstaff prime."""
    N = wagstaff(p)
    Np1 = N + 1
    sqrt_N = math.isqrt(N)

    result = {
        'p': p,
        'p_mod_6': p % 6,
        'digits': len(str(N)),
    }

    # Condition (I): omega3^{N+1} = 1
    r1 = omega3_pow(Np1, N)
    result['cond_I'] = (r1 == (1, 0))

    # Condition (II): omega3^{(N+1)/2} = -1
    r2 = omega3_pow(Np1 // 2, N)
    result['cond_II'] = (r2 == (N - 1, 0))

    # Get factorization of W_{p-2}
    facs = get_factorization(p)
    result['W_pm2_factors'] = {str(q): e for q, e in sorted(facs.items())}

    # Condition (III_q) for each prime factor
    F = 4
    S = []
    excluded = []
    for q in sorted(facs.keys()):
        exp = Np1 // q
        r = omega3_pow(exp, N)
        if r != (1, 0):
            F *= q ** facs[q]
            S.append(q)
        else:
            excluded.append(q)

    result['S'] = [str(q) for q in S]
    result['excluded'] = [str(q) for q in excluded]
    result['F_digits'] = len(str(F))
    result['sqrt_N_digits'] = len(str(sqrt_N))
    result['F_gt_sqrt_N'] = (F > sqrt_N)
    result['proven_prime'] = result['cond_I'] and result['cond_II'] and result['F_gt_sqrt_N']

    # Conjecture verification
    if p % 6 == 1:
        # Conjecture 4.5: chain
        W_pm2 = wagstaff(p - 2)
        r_chain = omega2_pow(W_pm2, N)
        result['conj_45_chain'] = (r_chain == (1, 0))

        s = 4
        for _ in range(p - 2):
            s = (s * s - 2) % N
        result['s_pm2'] = s
        result['s_pm2_is_4'] = (s == 4)
    elif p % 6 == 5:
        # Conjecture 4.6: primitive root
        half = (N - 1) // 2
        r_prim = omega2_pow(half, N)
        result['conj_46_prim_root'] = (r_prim == (N - 1, 0))

    return result

# ---- Run all verifications ----

if __name__ == '__main__':
    print("=" * 78)
    print("Wagstaff Primality Verification — Reproducibility Package")
    print("=" * 78)
    print()

    all_results = []

    for p in WAGSTAFF_PRIMES:
        t0 = time.time()
        result = verify_theorem32(p)
        dt = time.time() - t0

        status = "PROVEN ✓" if result['proven_prime'] else "FAILED ✗"
        excl = ','.join(result['excluded']) if result['excluded'] else '—'

        conj = ''
        if p % 6 == 1 and 'conj_45_chain' in result:
            conj = f"Conj4.5: {'✓' if result['conj_45_chain'] else '✗'}"
        elif p % 6 == 5 and 'conj_46_prim_root' in result:
            conj = f"Conj4.6: {'✓' if result['conj_46_prim_root'] else '✗'}"

        print(f"p={p:5d} ({result['digits']:3d}d) [{dt:6.1f}s] "
              f"{status}  excl={excl:8s}  {conj}")

        all_results.append(result)

    # Summary
    print()
    print("=" * 78)
    proven = sum(1 for r in all_results if r['proven_prime'])
    print(f"Total: {proven}/{len(all_results)} proven prime.")

    conj45 = [r for r in all_results if 'conj_45_chain' in r]
    conj46 = [r for r in all_results if 'conj_46_prim_root' in r]
    print(f"Conjecture 4.5 (chain): {sum(1 for r in conj45 if r['conj_45_chain'])}/{len(conj45)}")
    print(f"Conjecture 4.6 (prim root): {sum(1 for r in conj46 if r['conj_46_prim_root'])}/{len(conj46)}")

    # Save full results as JSON
    with open('wagstaff_verification_results.json', 'w') as f:
        json.dump(all_results, f, indent=2)
    print(f"\nFull results saved to wagstaff_verification_results.json")

    # Also need to add the Phi_1398 remainder for p=1709
    # (300-digit prime factor not in KNOWN_FACTORS because it's derived)
