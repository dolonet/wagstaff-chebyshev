#!/usr/bin/env python3
"""
CP164L: Extended R3 sweep at r ≤ 10¹¹, tracking sum(1/p).

Goal: decide whether the cluster's 0/1.835M R3-failure pattern at r ≤ 10¹⁰
is structural or random under the naive Chebotarev null.

Naive Chebotarev (CP145 + linear disjointness) predicts χ_r(α) uniform on
μ_p, so each split factor independently has P(R3 fails) ≈ 1/p. Expected
rogue count is therefore E = Σ 1/p over split factors. CP164f extrapolated
E(r ≤ 10¹⁰) ≈ 0.4 from a small-R sample with log-log slope ~0.17.

This sweep computes E exactly at r ≤ 10¹¹ (10x cp161). If we still see 0
fails AND E now exceeds ~1.5–2.0, the case for structural R3 starts to
firm up. If a fail appears, the structural premise is dead and Angle B
(p-th power reciprocity for α=1+√2) has no target.

For each prime r in this shard's range:
  - Compute ord_r(2). If ord = 2p with p prime, then r | W_p.
  - If r ≡ 1 mod 8 (split):
        compute ord_r(α) where α = 1+√2 (in F_r);
        test R3 (p | ord);
        accumulate inv_p_sum += 1/p.
  - If r ≡ 3 mod 8 (inert): test v₂(ord_r(ω₃)) = 2 (cheap free check).

Sharding: contiguous chunks of [5, R_MAX) across N shards.

Key difference from cp161_r3_sweep_cluster.py:
  - R_MAX = 10¹¹ (was 10¹⁰).
  - Tracks sum_inv_p per shard.
  - Uses α = 1+√2 directly (cp161 used ω₃ = 3+2√2 = α²; ω₃ has same
    p-th power character but α is the natural element for R3).

Output line: DONE shard=<i> ... split(yes/no)=<y>/<n> inv_p_sum=<x.xxxxxx>
"""
import sys, time

R_MAX = 10**11
SEGMENT = 50_000_000

# ===== Self-contained primality + factoring =====

WITNESSES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)

def is_prime(n):
    if n < 2: return False
    for p in (2,3,5,7,11,13,17,19,23,29,31,37):
        if n == p: return True
        if n % p == 0: return False
    d = n - 1; s = 0
    while d % 2 == 0:
        d //= 2; s += 1
    for a in WITNESSES:
        if a >= n: continue
        x = pow(a, d, n)
        if x == 1 or x == n - 1: continue
        for _ in range(s - 1):
            x = (x * x) % n
            if x == n - 1: break
        else:
            return False
    return True

def factor_small(n):
    """Trial division for n < ~10¹¹; returns dict {p: e}."""
    fs = {}
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47):
        while n % p == 0:
            fs[p] = fs.get(p, 0) + 1
            n //= p
        if n == 1: return fs
    p = 53
    import math
    lim = int(math.isqrt(n)) + 1
    while p <= lim and n > 1:
        while n % p == 0:
            fs[p] = fs.get(p, 0) + 1
            n //= p
            lim = int(math.isqrt(n)) + 1
        p += 2
    if n > 1:
        fs[n] = fs.get(n, 0) + 1
    return fs

# ===== Quadratic ring arithmetic over Z[√D]/N =====

def mul_mod(a1, b1, a2, b2, D, N):
    return ((a1*a2 + b1*b2*D) % N, (a1*b2 + a2*b1) % N)

def pow_quad_mod(a, b, exp, D, N):
    ra, rb = 1, 0
    a, b = a % N, b % N
    while exp > 0:
        if exp & 1:
            ra, rb = mul_mod(ra, rb, a, b, D, N)
        a, b = mul_mod(a, b, a, b, D, N)
        exp >>= 1
    return ra, rb

def order_int(g, N, bound, bound_factors):
    order = bound
    for q, e in bound_factors.items():
        for _ in range(e):
            if order % q == 0:
                test = order // q
                if pow(g, test, N) == 1:
                    order = test
                else:
                    break
    return order

def order_in_quad(a, b, D, N, bound, bound_factors):
    order = bound
    for q, e in bound_factors.items():
        for _ in range(e):
            if order % q == 0:
                test = order // q
                ra, rb = pow_quad_mod(a, b, test, D, N)
                if ra == 1 and rb == 0:
                    order = test
                else:
                    break
    return order

def v2(n):
    if n == 0: return 99
    v = 0
    while n % 2 == 0:
        v += 1; n //= 2
    return v

# ===== Segmented sieve =====

def primes_in_range(lo, hi):
    if lo < 2: lo = 2
    if hi <= lo: return
    import math
    sqrt_hi = int(math.isqrt(hi)) + 1
    base_sieve = bytearray([1]) * (sqrt_hi + 1)
    base_sieve[0] = base_sieve[1] = 0
    for i in range(2, int(math.isqrt(sqrt_hi)) + 1):
        if base_sieve[i]:
            for j in range(i*i, sqrt_hi + 1, i):
                base_sieve[j] = 0
    base_primes = [i for i in range(2, sqrt_hi + 1) if base_sieve[i]]
    seg_size = SEGMENT
    seg_lo = lo
    while seg_lo < hi:
        seg_hi = min(seg_lo + seg_size, hi)
        size = seg_hi - seg_lo
        seg = bytearray([1]) * size
        for p in base_primes:
            start = max(p * p, ((seg_lo + p - 1) // p) * p)
            for j in range(start - seg_lo, size, p):
                seg[j] = 0
        for i in range(size):
            if seg[i]:
                v = seg_lo + i
                if v >= 2:
                    yield v
        seg_lo = seg_hi

# ===== Main shard worker =====

def main():
    if len(sys.argv) != 3:
        print("usage: cp164L_extended_split_sweep.py <shard_id> <num_shards>", file=sys.stderr)
        sys.exit(1)
    shard_id = int(sys.argv[1])
    num_shards = int(sys.argv[2])

    chunk_size = R_MAX // num_shards + 1
    lo = max(5, shard_id * chunk_size)
    hi = min(R_MAX, (shard_id + 1) * chunk_size)
    if lo >= hi:
        print(f"DONE shard={shard_id} (empty range)", flush=True)
        return

    t0 = time.time()
    split_yes = 0
    split_no = 0
    inert_v2_2 = 0
    inert_other = 0
    inv_p_sum = 0.0
    last_report = t0
    count_r = 0
    fail_examples = []

    for r in primes_in_range(lo, hi):
        count_r += 1
        rm8 = r % 8
        if rm8 not in (1, 3): continue
        rm1 = r - 1
        rm1_factors = factor_small(rm1)
        o = order_int(2, r, rm1, rm1_factors)
        if o % 2 != 0: continue
        p_cand = o // 2
        if p_cand < 5 or not is_prime(p_cand): continue
        p = p_cand
        if pow(2, p, r) != r - 1: continue
        if r == 3: continue

        if rm8 == 1:
            # SPLIT: test R3 via ord_r(α) where α = 1+√2.
            # In F_r, √2 exists; the order of α in (Z[√2]/r)* equals the
            # order of α as element of F_r* under the splitting map.
            # We compute it in the quadratic ring (avoids needing √2 mod r),
            # which gives the same order divisor.
            order_alpha = order_in_quad(1, 1, 2, r, rm1, rm1_factors)
            inv_p_sum += 1.0 / p
            if order_alpha % p == 0:
                split_yes += 1
            else:
                split_no += 1
                fail_examples.append(('SPLIT_FAIL', p, r, order_alpha))
                if len(fail_examples) <= 20:
                    print(f"  R3 FAIL shard={shard_id}: p={p}, r={r}, ord_α={order_alpha}", flush=True)
        elif rm8 == 3:
            # INERT: free v₂ check on ω₃ = 3+2√2 in (Z[√2]/r)*
            rp1 = r + 1
            rp1_factors = factor_small(rp1)
            order_omega = order_in_quad(3, 2, 2, r, rp1, rp1_factors)
            v2_ord = v2(order_omega)
            if v2_ord == 2:
                inert_v2_2 += 1
            else:
                inert_other += 1
                fail_examples.append(('INERT_v2', p, r, order_omega, v2_ord))
                if len(fail_examples) <= 10:
                    print(f"  INERT v₂≠2 shard={shard_id}: p={p}, r={r}, ord={order_omega}, v₂={v2_ord}", flush=True)

        now = time.time()
        if now - last_report > 120:
            elapsed = now - t0
            total = split_yes + split_no + inert_v2_2 + inert_other
            print(f"  shard={shard_id} r~{r} count={count_r} hits={total} split(yes/no)={split_yes}/{split_no} inv_p_sum={inv_p_sum:.6f} inert(v₂2/other)={inert_v2_2}/{inert_other} {elapsed:.0f}s", flush=True)
            last_report = now

    elapsed = time.time() - t0
    total = split_yes + split_no + inert_v2_2 + inert_other
    print(f"DONE shard={shard_id} range=[{lo},{hi}) primes={count_r} hits={total} split(yes/no)={split_yes}/{split_no} inv_p_sum={inv_p_sum:.6f} inert(v₂2/other)={inert_v2_2}/{inert_other} elapsed={elapsed:.1f}s", flush=True)

if __name__ == "__main__":
    main()
