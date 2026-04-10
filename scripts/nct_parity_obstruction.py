#!/usr/bin/env python3
"""
NCT Parity Obstruction Analysis
================================

The "No Compatible Triple" (NCT) conjecture states: there is no triple (d, r, p) with
  - d > 1 odd
  - r prime, r ≡ 1 mod 8, r | U_{4d} (Pell sequence)
  - ord_r(2) = 2p, p prime
  - d | (2^{p-2} + 1) / 3

PARITY OBSTRUCTION:
  For d | W_{p-2} = (2^{p-2}+1)/3 with p >= 5 prime (so p-2 is odd),
  we need 2^{p-2} ≡ -1 mod each prime power q^a || d (q != 3).

  This requires ord_{q^a}(2) = 2e with e ODD (since e | p-2 and p-2 is odd).

  So d is BLOCKED if any prime q | d (q != 3) has:
    (a) ord_q(2) odd  =>  -1 not in <2> mod q, OR
    (b) v_2(ord_q(2)) >= 2  =>  e = ord_q(2)/2 is even, can't divide odd p-2

  d is UNBLOCKED iff every prime q | d with q != 3 has v_2(ord_q(2)) = 1.
"""

from math import gcd
from sympy import factorint, isprime, nextprime


def multiplicative_order(a, n):
    """Compute ord_n(a), the multiplicative order of a modulo n."""
    if gcd(a, n) != 1:
        return None
    order = 1
    current = a % n
    while current != 1:
        current = (current * a) % n
        order += 1
    return order


def v2(n):
    """2-adic valuation of n: largest k such that 2^k | n."""
    if n == 0:
        return float('inf')
    k = 0
    while n % 2 == 0:
        n //= 2
        k += 1
    return k


def classify_prime_by_mod8(q):
    """
    Classify a prime q != 2, 3 by its residue mod 8 and the 2-adic valuation
    of ord_q(2).

    Returns (q_mod_8, v2_ord, ord_q_2, status)
    where status is 'blocked' or 'unblocked'.
    """
    if q == 2 or q == 3:
        return None
    ord_q = multiplicative_order(2, q)
    v = v2(ord_q)
    qmod8 = q % 8
    if v == 1:
        status = 'unblocked'
    else:
        status = 'blocked'
    return (qmod8, v, ord_q, status)


def classify_d(d):
    """
    Classify odd d > 1: is it blocked or unblocked by the parity obstruction?

    Returns (is_blocked, blocking_info) where blocking_info lists the
    blocking primes and reasons.
    """
    factors = factorint(d)
    blocking = []
    unblocked_info = []

    for q, a in factors.items():
        if q == 3:
            continue  # skip q = 3
        if q == 2:
            continue  # d is odd, but just in case

        ord_q = multiplicative_order(2, q)
        v = v2(ord_q)

        if v == 0:
            blocking.append((q, a, ord_q, v, 'ord_q(2) is odd => -1 not in <2> mod q'))
        elif v >= 2:
            blocking.append((q, a, ord_q, v, f'v_2(ord_q(2)) = {v} >= 2 => e even'))
        else:
            # v == 1: unblocked contribution
            e = ord_q // 2  # this is odd
            # For prime powers q^a with a >= 2, check lifting
            if a >= 2:
                ord_qa = multiplicative_order(2, q**a)
                v_qa = v2(ord_qa)
                if v_qa != 1:
                    blocking.append((q, a, ord_qa, v_qa,
                                     f'v_2(ord_{{q^{a}}}(2)) = {v_qa} != 1 at prime power'))
                    continue
                e = ord_qa // 2
            unblocked_info.append((q, a, ord_q, e))

    is_blocked = len(blocking) > 0
    return is_blocked, blocking, unblocked_info


def compute_congruence_condition(unblocked_info):
    """
    For an unblocked d, compute the required congruence class for p.

    Each prime power q^a in d requires e_i | p-2 where e_i = ord_{q^a}(2)/2.
    So p-2 ≡ 0 mod lcm(e_1, ..., e_k).

    Returns (L, description) where L = lcm of all e_i.
    """
    from math import lcm as math_lcm
    if not unblocked_info:
        return 1, "no non-3 prime factors, p arbitrary (mod trivial)"

    L = 1
    details = []
    for q, a, ord_q, e in unblocked_info:
        L = math_lcm(L, e)
        if a == 1:
            details.append(f"e({q}) = {e}")
        else:
            details.append(f"e({q}^{a}) = {e}")

    return L, f"p ≡ 2 mod {2*L} (i.e., p-2 ≡ 0 mod {L}); components: {', '.join(details)}"


def main():
    print("=" * 80)
    print("NCT PARITY OBSTRUCTION ANALYSIS")
    print("=" * 80)

    # =========================================================================
    # PART 1: Verify mod-8 classification of primes
    # =========================================================================
    print("\n" + "=" * 80)
    print("PART 1: CLASSIFICATION BY q mod 8")
    print("=" * 80)

    mod8_classes = {1: [], 3: [], 5: [], 7: []}

    q = 5  # start from first prime > 3
    while q < 500:
        if q != 2 and q != 3:
            result = classify_prime_by_mod8(q)
            if result:
                qmod8, v, ord_q, status = result
                mod8_classes[qmod8].append((q, v, ord_q, status))
        q = nextprime(q)

    for r in [3, 5, 7, 1]:
        primes_in_class = mod8_classes[r]
        print(f"\n--- q ≡ {r} mod 8 ---")

        # Check if all have the same v2 pattern
        v2_vals = [v for _, v, _, _ in primes_in_class]
        statuses = [s for _, _, _, s in primes_in_class]

        if r == 3:
            expected = "v_2(ord_q(2)) = 1 (UNBLOCKED)"
            all_match = all(v == 1 for v in v2_vals)
        elif r == 5:
            expected = "v_2(ord_q(2)) >= 2 (BLOCKED)"
            all_match = all(v >= 2 for v in v2_vals)
        elif r == 7:
            expected = "v_2(ord_q(2)) = 0, odd order (BLOCKED)"
            all_match = all(v == 0 for v in v2_vals)
        elif r == 1:
            expected = "varies (need individual check)"
            all_match = True  # no uniform expectation

        print(f"  Expected: {expected}")
        print(f"  Number of primes checked: {len(primes_in_class)}")

        if r != 1:
            print(f"  All match expected pattern: {all_match}")
            if not all_match:
                exceptions = [(q, v, o) for q, v, o, _ in primes_in_class
                              if (r == 3 and v != 1) or (r == 5 and v < 2) or (r == 7 and v != 0)]
                print(f"  EXCEPTIONS: {exceptions}")
        else:
            # For q ≡ 1 mod 8, show detailed breakdown
            v2_dist = {}
            for q_val, v, ord_q, status in primes_in_class:
                v2_dist.setdefault(v, []).append(q_val)

            print(f"  Distribution of v_2(ord_q(2)):")
            for v_val in sorted(v2_dist.keys()):
                qs = v2_dist[v_val]
                status_str = "UNBLOCKED" if v_val == 1 else "BLOCKED"
                print(f"    v_2 = {v_val} ({status_str}): {len(qs)} primes: {qs}")

        # Show first few examples
        print(f"  First examples:")
        for q_val, v, ord_q, status in primes_in_class[:8]:
            print(f"    q={q_val}: ord_q(2)={ord_q}, v_2={v}, {status}")

    # =========================================================================
    # PART 2: Classify all odd d from 5 to 500
    # =========================================================================
    print("\n" + "=" * 80)
    print("PART 2: CLASSIFY ALL ODD d FROM 5 TO 500")
    print("=" * 80)

    blocked_list = []
    unblocked_list = []

    for d in range(5, 501, 2):  # odd d from 5 to 500
        if d % 2 == 0:
            continue

        is_blocked, blocking, unblocked_info = classify_d(d)

        if is_blocked:
            blocked_list.append((d, blocking))
        else:
            L, desc = compute_congruence_condition(unblocked_info)
            unblocked_list.append((d, unblocked_info, L, desc))

    total = len(blocked_list) + len(unblocked_list)
    pct_blocked = 100.0 * len(blocked_list) / total

    print(f"\n  Total odd d in [5, 500]: {total}")
    print(f"  Blocked:   {len(blocked_list)} ({pct_blocked:.1f}%)")
    print(f"  Unblocked: {len(unblocked_list)} ({100-pct_blocked:.1f}%)")

    # =========================================================================
    # PART 3: List all unblocked d
    # =========================================================================
    print("\n" + "=" * 80)
    print("PART 3: ALL UNBLOCKED d VALUES")
    print("=" * 80)

    print(f"\n  {'d':>5}  {'factorization':>20}  {'lcm(e_i)':>10}  congruence condition on p")
    print("  " + "-" * 75)

    for d, unblocked_info, L, desc in unblocked_list:
        factors = factorint(d)
        fact_str = " * ".join(f"{q}^{a}" if a > 1 else str(q) for q, a in sorted(factors.items()))
        print(f"  {d:>5}  {fact_str:>20}  {L:>10}  {desc}")

    # =========================================================================
    # PART 4: Detailed analysis of q ≡ 1 mod 8 primes
    # =========================================================================
    print("\n" + "=" * 80)
    print("PART 4: DETAILED ANALYSIS OF q ≡ 1 mod 8 PRIMES < 500")
    print("=" * 80)

    print(f"\n  {'q':>5}  {'q mod 8':>7}  {'ord_q(2)':>10}  {'v_2':>4}  {'e=ord/2':>8}  status")
    print("  " + "-" * 55)

    q1mod8_unblocked = []
    q1mod8_all = []

    q = 5
    while q < 500:
        if q % 8 == 1 and isprime(q):
            ord_q = multiplicative_order(2, q)
            v = v2(ord_q)
            status = "UNBLOCKED" if v == 1 else "BLOCKED"
            e = ord_q // 2 if v >= 1 else "N/A"
            print(f"  {q:>5}  {q%8:>7}  {ord_q:>10}  {v:>4}  {str(e):>8}  {status}")
            q1mod8_all.append((q, ord_q, v, status))
            if v == 1:
                q1mod8_unblocked.append(q)
        q = nextprime(q)

    print(f"\n  Total q ≡ 1 mod 8 primes < 500: {len(q1mod8_all)}")
    print(f"  Unblocked (v_2 = 1): {len(q1mod8_unblocked)}")
    if q1mod8_unblocked:
        print(f"  Unblocked q values: {q1mod8_unblocked}")
    else:
        print(f"  ALL q ≡ 1 mod 8 primes < 500 are BLOCKED!")

    # =========================================================================
    # PART 5: Structural summary
    # =========================================================================
    print("\n" + "=" * 80)
    print("PART 5: STRUCTURAL SUMMARY")
    print("=" * 80)

    # Categorize unblocked d by their prime factorization structure
    pure_3_power = []
    only_3mod8_primes = []
    has_1mod8_prime = []

    for d, unblocked_info, L, desc in unblocked_list:
        factors = factorint(d)
        non3_primes = [q for q in factors if q != 3]

        if not non3_primes:
            pure_3_power.append(d)
        elif all(q % 8 == 3 for q in non3_primes):
            only_3mod8_primes.append(d)
        else:
            has_1mod8_prime.append(d)

    print(f"\n  Unblocked d categories:")
    print(f"    Pure powers of 3:           {len(pure_3_power):>3}  values: {pure_3_power[:20]}{'...' if len(pure_3_power) > 20 else ''}")
    print(f"    Only q ≡ 3 mod 8 factors:   {len(only_3mod8_primes):>3}  values: {only_3mod8_primes[:20]}{'...' if len(only_3mod8_primes) > 20 else ''}")
    print(f"    Has q ≡ 1 mod 8 factor:     {len(has_1mod8_prime):>3}  values: {has_1mod8_prime[:20]}{'...' if len(has_1mod8_prime) > 20 else ''}")

    # Check: which q ≡ 1 mod 8 primes appear in unblocked d?
    q1mod8_in_unblocked = set()
    for d, unblocked_info, L, desc in unblocked_list:
        for q, a, ord_q, e in unblocked_info:
            if q % 8 == 1:
                q1mod8_in_unblocked.add(q)

    if q1mod8_in_unblocked:
        print(f"\n  q ≡ 1 mod 8 primes appearing in unblocked d: {sorted(q1mod8_in_unblocked)}")
    else:
        print(f"\n  NO q ≡ 1 mod 8 primes appear in any unblocked d!")

    # =========================================================================
    # PART 6: Key theorem summary
    # =========================================================================
    print("\n" + "=" * 80)
    print("PART 6: KEY THEOREMS")
    print("=" * 80)

    print("""
  THEOREM (Parity Obstruction by Residue Class):
    Let q be an odd prime, q != 3. Then:

    (1) q ≡ 7 mod 8  =>  ord_q(2) is ODD  =>  d with q | d is BLOCKED
        Proof: (2/q) = 1 when q ≡ ±1 mod 8, but q ≡ 7 mod 8 means
               q-1 ≡ 6 mod 8, so v_2(q-1) = 1. Since ord_q(2) | q-1
               and (2/q) = +1, we have 2^{(q-1)/2} ≡ 1 mod q.
               But (q-1)/2 is odd (since v_2(q-1) = 1), so
               ord_q(2) | (q-1)/2 which is odd  =>  ord_q(2) is odd.

    (2) q ≡ 5 mod 8  =>  v_2(ord_q(2)) >= 2  =>  d with q | d is BLOCKED
        Proof: (2/q) = -1 when q ≡ 3,5 mod 8, so 2^{(q-1)/2} ≡ -1 mod q.
               Thus ord_q(2) does not divide (q-1)/2, meaning
               v_2(ord_q(2)) > v_2((q-1)/2) = v_2(q-1) - 1.
               Since q ≡ 5 mod 8, v_2(q-1) >= 2, so v_2(ord_q(2)) >= 2.

    (3) q ≡ 3 mod 8  =>  v_2(ord_q(2)) = 1  =>  UNBLOCKED
        Proof: (2/q) = -1 so 2^{(q-1)/2} ≡ -1 mod q, hence ord_q(2)
               does not divide (q-1)/2. Since v_2(q-1) = 1 (as q ≡ 3 mod 8),
               (q-1)/2 is odd, and ord_q(2) = 2 * (odd), so v_2 = 1.

    (4) q ≡ 1 mod 8  =>  VARIES (depends on individual q)
        Since v_2(q-1) >= 3, the 2-adic valuation of ord_q(2) can be
        0, 1, 2, 3, ..., up to v_2(q-1).

  COROLLARY: An odd d > 1 is AUTOMATICALLY BLOCKED if it has any prime
  factor q with q ≡ 5 or 7 mod 8 (and q != 3).

  The only unblocked d are those whose prime factors are:
    - 3, or
    - primes q ≡ 3 mod 8, or
    - primes q ≡ 1 mod 8 with v_2(ord_q(2)) = 1 (rare/nonexistent < 500)
""")

    # =========================================================================
    # PART 7: Density estimate
    # =========================================================================
    print("=" * 80)
    print("PART 7: ASYMPTOTIC DENSITY")
    print("=" * 80)

    # Count primes by residue class
    counts = {1: 0, 3: 0, 5: 0, 7: 0}
    q = 5
    while q < 10000:
        counts[q % 8] += 1
        q = nextprime(q)

    print(f"\n  Primes < 10000 by residue mod 8:")
    for r in [1, 3, 5, 7]:
        print(f"    q ≡ {r} mod 8: {counts[r]}")

    print(f"""
  Since primes are equidistributed mod 8, about 50% of odd primes > 3
  fall in the classes q ≡ 5 or 7 mod 8 (which are always blocked).

  For a random odd d > 1 with k distinct prime factors (all > 3):
    - Prob(all factors ≡ 3 mod 8) ~ (1/4)^k  (if factors random mod 8)
    - So the fraction of unblocked d decreases exponentially with k.

  This explains the high blocking rate observed above.
""")

    print("=" * 80)
    print("COMPUTATION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    main()
