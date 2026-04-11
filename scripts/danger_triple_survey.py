#!/usr/bin/env python3
"""Danger-triple survey for the multi-factor reduction (Theorem 6.3).

Paper reference: Section 6.7.

For each admissible d <= D_MAX (odd, all prime factors q satisfy
v_2(ord_q(2)) = 1, i.e. Class III primes):
  1. Compute V_{2d} and extract its primitive part.
  2. Factor the primitive part (via sympy).
  3. Identify primitive inert divisors: primes r with r = 3 (mod 8)
     and r = -1 (mod 4d).
  4. For each such r, check whether ord_r(2)/2 is prime p with
     d | W_{p-2}.  If so, (d, r, p) is a danger triple.

A danger triple is 'phantom' if W_p is prime (so it can never appear
as a factor of a composite W_p).

Expected results for D_MAX = 1000:
  79 admissible d-values.
  5 danger triples: d in {57, 57, 67, 331, 683}.
  3 phantom (W_p prime), 2 real.

NOTE: For d > ~200, V_{2d} has hundreds of digits and sympy's
factorint may not complete in reasonable time.  The pre-computed
results in data/danger_triple_data.json include ECM-assisted
factorisations for these larger d-values.

Requires: sympy
Usage:    python3 danger_triple_survey.py [D_MAX]
"""
import sys
import json
import time
from sympy import factorint, isprime, gcd, divisors


def V(n):
    """Companion Pell number V_n(2,-1): V_0=2, V_1=2, V_n=2V_{n-1}+V_{n-2}."""
    if n == 0:
        return 2
    if n == 1:
        return 2
    a, b = 2, 2
    for _ in range(n - 1):
        a, b = b, 2 * b + a
    return b


def primitive_part(d):
    """Primitive part of V_{2d}: remove all factors shared with V_m for proper divisors m of 2d."""
    val = V(2 * d)
    n = 2 * d
    for m in sorted(divisors(n)):
        if m == n:
            continue
        g = gcd(val, V(m))
        while g > 1:
            val //= g
            g = gcd(val, g)
    return val


def is_admissible(d):
    """d is admissible if odd, > 1, and every prime factor q has v_2(ord_q(2)) = 1."""
    if d % 2 == 0 or d <= 1:
        return False
    for q in factorint(d):
        q = int(q)
        # ord_q(2): find smallest k with 2^k = 1 mod q
        x, k = 2 % q, 1
        while x != 1:
            x = x * 2 % q
            k += 1
            if k > q:
                return False
        # v_2(k) must be exactly 1
        if k % 2 != 0 or (k // 2) % 2 == 0:
            return False
    return True


def check_danger(d, r):
    """Check whether (d, r) yields a danger triple."""
    r = int(r)
    # Compute ord_r(2) via factorisation of r-1
    fac_r1 = factorint(r - 1)
    order = r - 1
    for p, e in fac_r1.items():
        p = int(p)
        for _ in range(e):
            if pow(2, order // p, r) == 1:
                order //= p
            else:
                break

    if order % 2 != 0:
        return None
    p_val = order // 2
    if not isprime(p_val):
        return None

    # Check d | W_{p-2}, i.e. for each prime power q^e || d,
    # 2^{p-2} = -1 (mod q^e)  [with adjustment at q = 3]
    for q, e in factorint(d).items():
        q = int(q)
        mod = (3 ** (e + 1)) if q == 3 else (q ** e)
        if pow(2, p_val - 2, mod) != mod - 1:
            return None

    # Check if W_p is prime (phantom).  For small p, W_p = r itself;
    # for large p, materialising 2^p is infeasible — mark as unknown.
    KNOWN_WAGSTAFF_PRIME_EXPONENTS = {
        3, 5, 7, 11, 13, 17, 19, 23, 31, 43, 61, 79, 101, 127,
        167, 191, 199, 313, 347, 701, 1709, 2617, 10501, 12391,
    }
    if p_val in KNOWN_WAGSTAFF_PRIME_EXPONENTS:
        phantom = True
    elif p_val < 10000 and isprime((pow(2, p_val) + 1) // 3):
        phantom = True
    else:
        phantom = False  # conservatively treat as real
    return {"d": d, "r": r, "p": p_val, "phantom": phantom}


def analyze_d(d, verbose=False, digit_limit=45):
    """Full analysis of one d-value.  Returns (d, status, dangers).
    Skips factorisation if primitive part exceeds digit_limit digits."""
    t0 = time.time()
    prim = primitive_part(d)
    digits = len(str(prim))
    if verbose:
        print(f"  d={d}: primitive part has {digits} digits", flush=True)

    if digits > digit_limit:
        if verbose:
            print(f"    skipped (>{digit_limit} digits; see data/danger_triple_data.json)", flush=True)
        return d, f"SKIPPED({digits}d)", [], digits, time.time() - t0

    fac = factorint(prim)

    product = 1
    for p, e in fac.items():
        product *= int(p) ** e
    complete = (product == prim)
    if not complete and isprime(prim // product):
        fac[prim // product] = 1
        complete = True

    # Identify primitive inert divisors
    dangers = []
    for r in fac:
        r = int(r)
        if r < 3:
            continue
        if r % 8 != 3:
            continue
        if r % (4 * d) != 4 * d - 1:
            continue
        hit = check_danger(d, r)
        if hit is not None:
            dangers.append(hit)

    cofactor_digits = len(str(prim // product)) if not complete else 0
    if complete:
        status = f"EXCLUDED" if not dangers else f"DANGER({len(dangers)})"
    else:
        status = f"PARTIAL({cofactor_digits}d)" if not dangers else f"DANGER+PARTIAL({len(dangers)})"

    return d, status, dangers, digits, time.time() - t0


def main():
    d_max = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    print(f"Danger-triple survey for admissible d <= {d_max}")
    print(f"(For the full d <= 1000 survey, see data/danger_triple_data.json)\n")

    admissible = [d for d in range(3, d_max + 1, 2) if is_admissible(d)]
    print(f"Admissible d-values: {len(admissible)}")

    all_dangers = []
    results = []
    for d in admissible:
        d_val, status, dangers, digits, elapsed = analyze_d(d, verbose=True, digit_limit=45)
        results.append({"d": d_val, "status": status, "digits": digits, "time": round(elapsed, 1)})
        for dt in dangers:
            ptype = "PHANTOM" if dt["phantom"] else "REAL"
            print(f"    *** DANGER TRIPLE: (d={dt['d']}, r={dt['r']}, p={dt['p']}) [{ptype}]")
            all_dangers.append(dt)

    print(f"\n{'='*60}")
    print(f"Survey complete: {len(admissible)} admissible d-values, {len(all_dangers)} danger triples")
    for dt in all_dangers:
        ptype = "PHANTOM" if dt["phantom"] else "REAL"
        print(f"  (d={dt['d']}, r={dt['r']}, p={dt['p']}) [{ptype}]")

    # Write JSON output
    out = {
        "d_max": d_max,
        "n_admissible": len(admissible),
        "admissible_d": admissible,
        "n_danger_triples": len(all_dangers),
        "danger_triples": all_dangers,
        "results": results
    }
    outfile = f"danger_triple_survey_d{d_max}.json"
    with open(outfile, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nResults written to {outfile}")


if __name__ == "__main__":
    main()
