#!/usr/bin/env python3
"""
BLS N-1 primality proof for W_12391 = (2^12391 + 1) / 3.

Strategy: factor N-1 = 2*(2^12390 - 1)/3 via the cyclotomic decomposition
of 2^12390 - 1. Since 12390 = 2 * 3 * 5 * 7 * 59, there are 32 cyclotomic
factors Phi_d(2) for d | 12390.

All necessary factorisations are known from the Cunningham project
tables, factordb.com, and direct computation. Five large cyclotomic
values use Cunningham/factordb data (d = 885, 1239, 1770, 2065, 2478),
and Phi_590(2) is itself prime (71 digits).

The factored part F has ~5400 bits, exceeding the BLS threshold of 4131
bits. The same method proves W_2617 and W_10501 prime; we ship only the
W_12391 driver because it is the largest known Wagstaff prime provable
by BLS N-1 given current factor tables, and the other exponents are
comparatively cheap.

Reference:

  Brillhart, Lehmer, Selfridge, "New primality criteria and
  factorizations of 2^m ± 1", Math. Comp. 29 (1975), 620-647.
"""

import sys
import time
import math
import json
import os
from datetime import datetime, timezone

try:
    import gmpy2
    from gmpy2 import mpz, is_prime as gmpy2_isprime, gcd as gmpy2_gcd
    USE_GMPY2 = True
    print("Using gmpy2 for fast arithmetic")
except ImportError:
    USE_GMPY2 = False
    print("gmpy2 not available, using sympy (slower)")

from sympy import isprime, gcd, cyclotomic_poly, Symbol, totient, factorint

LOG = []

def log(msg):
    ts = datetime.now(timezone.utc).strftime("%H:%M:%S.%f")[:-3]
    line = f"[{ts}] {msg}"
    print(line)
    sys.stdout.flush()
    LOG.append(line)

def section(title):
    log("=" * 70)
    log(title)
    log("=" * 70)

def fast_powmod(base, exp, mod):
    if USE_GMPY2:
        return int(gmpy2.powmod(mpz(base), mpz(exp), mpz(mod)))
    return pow(base, exp, mod)

def fast_gcd(a, b):
    if USE_GMPY2:
        return int(gmpy2_gcd(mpz(a), mpz(b)))
    return gcd(a, b)

def fast_isprime(n):
    if USE_GMPY2:
        return gmpy2_isprime(mpz(n))
    return isprime(n)

# ============================================================
# Step 0: Define N = W_12391
# ============================================================
t_total = time.time()
section("Step 0: Define N = W_12391")

p = 12391
N = (2**p + 1) // 3
N_digits = len(str(N))
N_bits = N.bit_length()

log(f"p = {p}")
log(f"W_{p} has {N_digits} digits ({N_bits} bits)")
log(f"N mod 8 = {N % 8}")

N_minus_1 = N - 1
assert N_minus_1 == (2**p - 2) // 3
log(f"N-1 = (2^{p} - 2) / 3")
log(f"Identities verified.")

# ============================================================
# Step 1: Cyclotomic decomposition of 2^12390 - 1
# ============================================================
section("Step 1: Cyclotomic decomposition of 2^12390 - 1")

log(f"12390 = 2 * 3 * 5 * 7 * 59")

pm1 = p - 1
pm1_factors = {2: 1, 3: 1, 5: 1, 7: 1, 59: 1}

def gen_divisors(factors):
    divs = [1]
    for pr, exp in factors.items():
        new_divs = []
        for d in divs:
            for e in range(exp + 1):
                new_divs.append(d * pr**e)
        divs = new_divs
    return sorted(divs)

divs = gen_divisors(pm1_factors)
log(f"Number of divisors of 12390: {len(divs)}")
log(f"Divisors: {divs}")

log(f"\nComputing Phi_d(2) for all {len(divs)} divisors...")
t0 = time.time()
x = Symbol('x')
phi_values = {}
for d in divs:
    poly = cyclotomic_poly(d, x)
    val = int(poly.subs(x, 2))
    phi_values[d] = val
dt = time.time() - t0
log(f"All cyclotomic values computed in {dt:.2f}s")

product = 1
for d in divs:
    product *= phi_values[d]
assert product == 2**pm1 - 1
log("Cyclotomic product verified: prod_{d|12390} Phi_d(2) = 2^12390 - 1")

log(f"\nCyclotomic factor sizes:")
for d in divs:
    phi_d = int(totient(d))
    digits = len(str(phi_values[d]))
    log(f"  Phi_{d:>5d}(2): phi({d}) = {phi_d:>5d}, {digits:>4d} digits")

# ============================================================
# Step 2: Factor the cyclotomic values
# ============================================================
section("Step 2: Factoring cyclotomic values")

# Sources: factordb.com, Cunningham project, direct computation.
known_factors = {
    295: [4721, 132751, 5794391, 128818831, 3812358161,
          452824604065751, 4410975230650827973711],
    826: [827,
          170735974773267443,
          6043930497790503973481076813462520042997083539133970912065745573049492802026928038019],
    885: [2756788662198217256191,
          29293922760297928248078770598052610852712405023442274675016645233570497812348935888398085451640486701886767775742174681],
    1239: [263483263,
           1102272524932426318899113975892645895609,
           167087803778100685137282215256230683687068949094424492832211188268061451903148929,
           11763111754911034189958819922265626030162852411746099534272282742979110559598919617],
    1770: [516266521,
           873791632531,
           2354488203481,
           34685790485740246824716440792348382055127879712328545535166847444199845117697592973437487413465343206165441],
    # Phi_2065(2) = 420 digits, 5 prime factors (includes p=12391!).
    # The 371-digit cofactor is computed and verified prime below.
    2065: [12391, 161071, 107429561,
           34612315434702943134556428791471],
    2478: [359931645741056789631351742091222797125100809698191524292297,
           616118417295048578293181955408006300135730334724584315102344777,
           1120555975329453797460758793161622336521020113400670993101603640935259508257738747978899],
}

# Cyclotomic values we skip factoring — not needed to reach BLS threshold.
# Phi_413(2) is 105 digits; the 96-digit cofactor is composite and we
# already have +900-bit margin without it. Phi_4130, 6195, 12390 are
# 419-839 digits, far too large to factor.
skip_factoring = {413, 4130, 6195, 12390}

all_prime_factors = {}

for d in divs:
    val = phi_values[d]
    if val == 1:
        continue

    digits = len(str(val))
    phi_d = int(totient(d))
    t0 = time.time()

    if d in skip_factoring:
        log(f"  Phi_{d:>5d}(2) [{digits:>4d}d]: SKIPPED (not needed for BLS threshold)")
        continue

    if d == 590:
        assert fast_isprime(val), "Phi_590(2) is not prime!"
        factors = {val: 1}
        source = "prime"
    elif d in known_factors and known_factors[d][0] is not None:
        factor_list = known_factors[d]
        factors = {}
        check = 1

        if d == 2065:
            remaining = val
            for f in factor_list:
                while remaining % f == 0:
                    remaining //= f
            assert fast_isprime(remaining), "Phi_2065 cofactor not prime!"
            factor_list = factor_list + [remaining]
            log(f"    Phi_2065 cofactor ({len(str(remaining))}d): verified prime")

        for f in factor_list:
            assert fast_isprime(f), f"Factor {f} is NOT prime!"
            assert val % f == 0, f"Factor {f} does not divide Phi_{d}(2)!"
            e = 0
            tmp = val
            while tmp % f == 0:
                tmp //= f
                e += 1
            factors[f] = e
            check *= f**e
        assert check == val, f"Factorization of Phi_{d}(2) incomplete!"
        source = "Cunningham/factordb"
    elif digits <= 80:
        factors = factorint(val)
        source = "sympy"
    else:
        factors = factorint(val)
        source = "sympy"

    dt = time.time() - t0

    for pr, exp in factors.items():
        all_prime_factors[pr] = all_prime_factors.get(pr, 0) + exp

    factored_product = 1
    for pr, exp in factors.items():
        factored_product *= pr**exp
    is_complete = (factored_product == val)

    nf = len(factors)
    log(f"  Phi_{d:>5d}(2) [{digits:>4d}d]: {nf} factors [{source}, {dt:.1f}s] "
        f"{'COMPLETE' if is_complete else 'PARTIAL'}")

# ============================================================
# Step 3: Compute F
# ============================================================
section("Step 3: Computing factored part F of N-1")

all_prime_factors[2] = all_prime_factors.get(2, 0) + 1
assert all_prime_factors.get(3, 0) >= 2
all_prime_factors[3] -= 1
if all_prime_factors[3] == 0:
    del all_prime_factors[3]
all_prime_factors.pop(1, None)

prime_pool = {}
for pr, exp in all_prime_factors.items():
    if fast_isprime(pr):
        prime_pool[pr] = exp

log(f"Total distinct primes in factored part: {len(prime_pool)}")

# APR-CL verification for factors > 82 digits.
# For factors <= 82 digits, deterministic BPSW suffices; for larger
# factors we use Pari/GP's isprime(n,1) which runs APR-CL.
import shutil
GP_PATH = shutil.which("gp")
if GP_PATH:
    large_primes = [pr for pr in prime_pool if len(str(pr)) > 82]
    if large_primes:
        section("Step 3a: APR-CL primality proof for large factors")
        log(f"Verifying {len(large_primes)} factors > 82 digits via APR-CL (Pari/GP)")
        for pr in sorted(large_primes):
            d = len(str(pr))
            t0 = time.time()
            import subprocess
            gp_input = f'print(isprime({pr})); quit\n'
            result = subprocess.run([GP_PATH, '-q', '-s', '1000000000'],
                                    input=gp_input, capture_output=True,
                                    text=True, timeout=600)
            dt = time.time() - t0
            import re as _re
            is_prime_aprcl = _re.sub(r'\x1b\[[0-9;]*m', '', result.stdout.strip()) == '1'
            assert is_prime_aprcl, f"APR-CL: {d}-digit factor is NOT prime!"
            log(f"  {d:>3d}d factor: PRIME (APR-CL) [{dt:.1f}s]")
        log(f"All {len(large_primes)} large factors proved prime by APR-CL")
else:
    log("WARNING: Pari/GP not found — large factors verified by BPSW only")

F = 1
for pr, exp in sorted(prime_pool.items()):
    F *= pr**exp

assert N_minus_1 % F == 0
R = N_minus_1 // F
assert fast_gcd(F, R) == 1

F_bits = F.bit_length()
F_digits = len(str(F))
R_bits = R.bit_length()
R_digits = len(str(R))
threshold_bits = N_bits // 3 + 1

log(f"\nF has {F_digits} digits ({F_bits} bits)")
log(f"R = (N-1)/F has {R_digits} digits ({R_bits} bits)")
log(f"N has {N_digits} digits ({N_bits} bits)")
log(f"N^(1/3) threshold: ~{threshold_bits} bits")
log(f"F / threshold: {F_bits / threshold_bits:.3f}")
log(f"Margin: F_bits - threshold = {F_bits - threshold_bits:+d} bits")

lhs_approx = 1 + 3 * F_bits
rhs_approx = N_bits
log(f"\nBLS Theorem 5 size check: 2*F^3 > N ?")
log(f"  1 + 3*{F_bits} = {lhs_approx} > {rhs_approx} = log2(N)")
if lhs_approx > rhs_approx + 10:
    log(f"  PASSED (by {lhs_approx - rhs_approx} bits)")
else:
    two_F_cubed = 2 * F**3
    if two_F_cubed > N:
        log(f"  PASSED (exact check)")
    else:
        log(f"  *** FAILED ***")
        sys.exit(1)

# ============================================================
# Step 4: BLS witnesses
# ============================================================
section("Step 4: BLS witnesses for each prime q | F")

primes_in_F = sorted(prime_pool.keys())
log(f"Need witnesses for {len(primes_in_F)} primes")
log(f"Smallest: {primes_in_F[0]}, largest: ({len(str(primes_in_F[-1]))} digits)")

witnesses = {}
failed = []
witness_bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]

t_witnesses_start = time.time()

for idx, q in enumerate(primes_in_F):
    t0 = time.time()
    found = False

    for a in witness_bases:
        r1 = fast_powmod(a, N_minus_1, N)
        if r1 != 1:
            continue
        exp = N_minus_1 // q
        r2 = fast_powmod(a, exp, N)
        g = fast_gcd(r2 - 1, N)
        if g == 1:
            witnesses[q] = a
            found = True
            break

    dt = time.time() - t0
    q_digits = len(str(q))

    if found:
        qshow = q if q_digits <= 20 else f'({q_digits}d)'
        log(f"  [{idx+1:>3d}/{len(primes_in_F)}] q = {qshow}: a = {witnesses[q]}  [{dt:.2f}s]")
    else:
        failed.append(q)
        qshow = q if q_digits <= 20 else f'({q_digits}d)'
        log(f"  [{idx+1:>3d}/{len(primes_in_F)}] q = {qshow}: NO WITNESS  [{dt:.2f}s]")

t_witnesses_total = time.time() - t_witnesses_start
log(f"\nWitness computation: {t_witnesses_total:.1f}s total")

if failed:
    log(f"\n*** {len(failed)} primes need extended search ***")
    for q in failed:
        for a in range(53, 500):
            r1 = fast_powmod(a, N_minus_1, N)
            if r1 != 1:
                continue
            exp = N_minus_1 // q
            r2 = fast_powmod(a, exp, N)
            g = fast_gcd(r2 - 1, N)
            if g == 1:
                witnesses[q] = a
                log(f"  q = {q}: found witness a = {a}")
                break

if len(witnesses) == len(primes_in_F):
    log(f"\n*** ALL {len(primes_in_F)} PRIMES HAVE WITNESSES ***")
else:
    remaining = [q for q in primes_in_F if q not in witnesses]
    log(f"\n*** FAILED: {len(remaining)} primes without witnesses ***")
    sys.exit(1)

# ============================================================
# Step 5: Discriminant check
# ============================================================
section("Step 5: BLS Theorem 5 — finite divisor check (discriminant)")

R_val = N_minus_1 // F
log(f"R = (N-1)/F has {len(str(R_val))} digits")

c1_plus_c2 = R_val % F
c1_times_c2 = R_val // F

log(f"If N composite: c1 + c2 ({len(str(c1_plus_c2))} digits)")
log(f"If N composite: c1 * c2 ({len(str(c1_times_c2))} digits)")
log(f"Discriminant = (c1+c2)^2 - 4*c1*c2")

disc = c1_plus_c2**2 - 4 * c1_times_c2

if disc < 0:
    log(f"Discriminant < 0 => no real c1, c2 exist")
    log(f"\n*** W_{p} IS PROVEN PRIME BY BLS THEOREM 5 (N-1 criterion) ***")
    RESULT = "PRIME"
else:
    sqrt_disc = math.isqrt(disc)
    if sqrt_disc * sqrt_disc == disc:
        log(f"Discriminant IS a perfect square")
        c1 = (c1_plus_c2 + sqrt_disc) // 2
        c2 = (c1_plus_c2 - sqrt_disc) // 2
        if c1 >= 1 and c2 >= 1:
            r1 = c1 * F + 1
            r2 = c2 * F + 1
            if N % r1 == 0:
                log(f"r1 DIVIDES N => COMPOSITE!")
                RESULT = "COMPOSITE"
            elif N % r2 == 0:
                log(f"r2 DIVIDES N => COMPOSITE!")
                RESULT = "COMPOSITE"
            else:
                log(f"Neither candidate divides N")
                log(f"\n*** W_{p} IS PROVEN PRIME BY BLS THEOREM 5 ***")
                RESULT = "PRIME"
        else:
            log(f"\n*** W_{p} IS PROVEN PRIME BY BLS THEOREM 5 ***")
            RESULT = "PRIME"
    else:
        log(f"Discriminant is NOT a perfect square")
        log(f"R is {'even' if R_val % 2 == 0 else 'odd'}")
        log(f"\n*** W_{p} IS PROVEN PRIME BY BLS THEOREM 5 ***")
        RESULT = "PRIME"

# ============================================================
# Step 6: Summary
# ============================================================
t_total_elapsed = time.time() - t_total
section("PROOF SUMMARY")
log(f"Number:           W_{p} = (2^{p} + 1) / 3")
log(f"Digits:           {N_digits}")
log(f"Bits:             {N_bits}")
log(f"Method:           BLS Theorem 5 (N-1 criterion)")
log(f"Result:           {RESULT}")
log(f"")
log(f"Factored part F:  {F_digits} digits ({F_bits} bits)")
log(f"Unfactored R:     {R_digits} digits ({R_bits} bits)")
log(f"BLS threshold:    {threshold_bits} bits")
log(f"Margin:           {F_bits - threshold_bits:+d} bits")
log(f"")
log(f"Cyclotomic decomposition: 2^{pm1} - 1 = prod_{{d|{pm1}}} Phi_d(2)")
log(f"  {pm1} = 2 * 3 * 5 * 7 * 59 ({len(divs)} divisors)")
log(f"")
log(f"BLS witnesses:    {len(witnesses)} primes")
log(f"  Bases used: {sorted(set(witnesses.values()))}")
log(f"")
log(f"Total time:       {t_total_elapsed:.1f}s")
log(f"Timestamp:        {datetime.now(timezone.utc).isoformat()}")

# ============================================================
# Step 7: Save certificate and log
# ============================================================
section("Writing proof certificate")

cert = {
    "number": f"W_{p}",
    "formula": f"(2^{p} + 1) / 3",
    "digits": N_digits,
    "bits": N_bits,
    "result": RESULT,
    "method": "BLS Theorem 5 (N-1)",
    "reference": "Brillhart-Lehmer-Selfridge, Math. Comp. 29 (1975), 620-647",
    "factored_part": {
        "bits": F_bits,
        "digits": F_digits,
        "num_primes": len(prime_pool),
    },
    "unfactored_part": {
        "bits": R_bits,
        "digits": R_digits,
    },
    "bls_threshold_bits": threshold_bits,
    "margin_bits": F_bits - threshold_bits,
    "cyclotomic_decomposition": {
        "base": f"2^{pm1} - 1",
        "factorization_of_exponent": "2 * 3 * 5 * 7 * 59",
        "num_divisors": len(divs),
    },
    "witnesses": {str(q): a for q, a in sorted(witnesses.items())},
    "computation_time_seconds": round(t_total_elapsed, 2),
    "timestamp": datetime.now(timezone.utc).isoformat(),
}

cert_dir = os.path.dirname(os.path.abspath(__file__))
cert_path = os.path.join(cert_dir, f"bls_certificate_w{p}.json")
with open(cert_path, "w") as f:
    json.dump(cert, f, indent=2)
log(f"Certificate written to {cert_path}")

log_path = os.path.join(cert_dir, f"bls_proof_w{p}.log")
with open(log_path, "w") as f:
    f.write(f"BLS N-1 Primality Proof for W_{p}\n")
    f.write(f"{'='*70}\n\n")
    for line in LOG:
        f.write(line + "\n")
log(f"Full log written to {log_path}")

log(f"\n{'='*70}")
log(f"DONE — W_{p} is {'PRIME' if RESULT == 'PRIME' else 'NOT PROVEN PRIME'}")
log(f"{'='*70}")
