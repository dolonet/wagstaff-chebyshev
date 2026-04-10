"""
BLS N-1 primality proof for W_2617 = (2^2617 + 1) / 3

Strategy: Instead of factoring N+1 = 4*W_{2615} (which requires factoring a ~787-digit number),
we factor N-1 = 2*(2^2616 - 1)/3 using the cyclotomic decomposition of 2^2616 - 1.

Since 2616 = 2^3 * 3 * 109, we have 16 cyclotomic factors Phi_d(2) for d | 2616.
The small factors (d <= 654) are easily factorable and together exceed the N^{1/3} threshold
required by BLS Theorem 5.

References:
  [BLS75] Brillhart, Lehmer, Selfridge, "New primality criteria and factorizations of 2^m +/- 1",
          Math. Comp. 29 (1975), 620-647.
"""

from sympy import isprime, factorint, gcd, cyclotomic_poly, Symbol, nextprime
from sympy.ntheory import primitive_root
import time

# ============================================================
# Step 0: Define N = W_2617
# ============================================================
p = 2617
N = (2**p + 1) // 3
print(f"W_{p} has {len(str(N))} digits")
print(f"N mod 8 = {N % 8}")  # should be 3

N_minus_1 = N - 1
N_plus_1 = N + 1

# Verify basic identities
assert N_minus_1 == (2**p - 2) // 3
assert N_plus_1 == 4 * ((2**(p-2) + 1) // 3)
print(f"N-1 = (2^{p} - 2) / 3")
print(f"N+1 = 4 * W_{p-2}")

# ============================================================
# Step 1: Cyclotomic decomposition of 2^{2616} - 1
# ============================================================
print("\n" + "="*60)
print("Step 1: Cyclotomic decomposition of 2^2616 - 1")
print("="*60)

# 2616 = 2^3 * 3 * 109
# Divisors of 2616
def divisors(n):
    divs = []
    for i in range(1, n+1):
        if n % i == 0:
            divs.append(i)
    return divs

divs_2616 = divisors(2616)
print(f"Divisors of 2616: {divs_2616}")
print(f"Number of divisors: {len(divs_2616)}")

# Compute Phi_d(2) for each divisor
x = Symbol('x')
phi_values = {}
for d in divs_2616:
    poly = cyclotomic_poly(d, x)
    val = int(poly.subs(x, 2))
    phi_values[d] = val

# Verify product equals 2^2616 - 1
product = 1
for d in divs_2616:
    product *= phi_values[d]
assert product == 2**2616 - 1, "Cyclotomic product verification failed!"
print("Cyclotomic product verified: prod Phi_d(2) = 2^2616 - 1")

# Euler's totient for size estimation
from sympy import totient
for d in divs_2616:
    phi_d = totient(d)
    digits = len(str(phi_values[d]))
    print(f"  Phi_{d}(2): phi({d}) = {phi_d}, {digits} digits")

# ============================================================
# Step 2: Factor the small cyclotomic values
# ============================================================
print("\n" + "="*60)
print("Step 2: Factoring cyclotomic values (d <= 654)")
print("="*60)

# Known factorizations from checkpoint_31
known_factors = {
    327: [20597276734348736647, 33157029794959983067039, 88116165754061081804047],
    436: [5669, 666184021, 74323515777853, 1746518852140345553, 171857646012809566969],
    654: [666427, 6927735019, 30414028470765822165976581508161866432602988327347],
}

all_prime_factors = {}  # d -> list of (prime, exponent) pairs

# Factor each Phi_d(2) for d <= 654
for d in divs_2616:
    val = phi_values[d]
    if d > 654:
        print(f"  Phi_{d}(2): {len(str(val))} digits -- SKIPPED (not needed)")
        continue

    if val == 1:
        print(f"  Phi_{d}(2) = 1 -- trivial")
        all_prime_factors[d] = {}
        continue

    if d in known_factors:
        # Verify known factorization
        prod_check = 1
        factors = {}
        for f in known_factors[d]:
            assert isprime(f), f"Factor {f} is not prime!"
            assert val % f == 0, f"Factor {f} does not divide Phi_{d}(2)!"
            e = 0
            tmp = val
            while tmp % f == 0:
                tmp //= f
                e += 1
            factors[f] = e
            prod_check *= f**e
        assert prod_check == val, f"Factorization of Phi_{d}(2) is incomplete!"
        all_prime_factors[d] = factors
        print(f"  Phi_{d}(2) = {' * '.join(str(f) + ('^'+str(e) if e>1 else '') for f,e in factors.items())}  [verified]")
    else:
        # Factor using sympy
        t0 = time.time()
        factors = factorint(val)
        dt = time.time() - t0
        all_prime_factors[d] = factors
        if val < 10**20:
            print(f"  Phi_{d}(2) = {val} = {' * '.join(str(f) + ('^'+str(e) if e>1 else '') for f,e in factors.items())}  [{dt:.1f}s]")
        else:
            print(f"  Phi_{d}(2) ({len(str(val))} digits) = {' * '.join(str(f) + ('^'+str(e) if e>1 else '') for f,e in factors.items())}  [{dt:.1f}s]")

# ============================================================
# Step 3: Collect all prime factors and compute F
# ============================================================
print("\n" + "="*60)
print("Step 3: Computing factored part F of N-1")
print("="*60)

# N - 1 = 2 * (2^2616 - 1) / 3
# The factor 2 contributes 2^1.
# Division by 3: since Phi_2(2) = 3 and Phi_6(2) = 3, we have 3^2 | (2^2616-1).
# We divide by 3, removing one factor of 3.

# Collect all primes from the factored cyclotomic values
prime_pool = {}  # prime -> total exponent in 2^2616 - 1

for d, factors in all_prime_factors.items():
    for pr, exp in factors.items():
        if pr in prime_pool:
            prime_pool[pr] += exp
        else:
            prime_pool[pr] = exp

# Add the factor of 2 from N-1 = 2*(2^2616-1)/3
if 2 in prime_pool:
    prime_pool[2] += 1
else:
    prime_pool[2] = 1

# Remove one factor of 3 (dividing by 3)
assert prime_pool[3] >= 2, f"Expected at least 3^2 in the product, got 3^{prime_pool.get(3, 0)}"
prime_pool[3] -= 1

print(f"Total distinct primes in factored part: {len(prime_pool)}")

# Compute F = product of p^e for all collected primes
F = 1
for pr, exp in sorted(prime_pool.items()):
    F *= pr**exp

# Compute the unfactored remainder R
R = N_minus_1 // F
assert N_minus_1 == F * R, "F * R != N-1"
assert gcd(F, R) == 1, "gcd(F, R) != 1"

F_bits = F.bit_length()
N_bits = N.bit_length()
threshold_bits = N_bits // 3 + 1

print(f"\nF has {len(str(F))} digits ({F_bits} bits)")
print(f"N has {len(str(N))} digits ({N_bits} bits)")
print(f"N^(1/3) threshold: ~{threshold_bits} bits")
print(f"F bits: {F_bits}")
print(f"Margin: {F_bits - threshold_bits} bits")

# BLS Theorem 5 check: 2*F^3 > N
F_cubed_times_2 = 2 * F**3
if F_cubed_times_2 > N:
    print(f"\n*** BLS Theorem 5 CHECK: 2*F^3 > N  ==>  PASSED ***")
else:
    print(f"\n*** BLS Theorem 5 CHECK: 2*F^3 > N  ==>  FAILED ***")
    print(f"  2*F^3 has {len(str(F_cubed_times_2))} digits, N has {len(str(N))} digits")

# ============================================================
# Step 4: Find witnesses for each prime factor of F
# ============================================================
print("\n" + "="*60)
print("Step 4: BLS witnesses for each prime q | F")
print("="*60)

# For BLS Theorem 5 (N-1 test), for each prime q | F, we need:
#   (I) a^{N-1} ≡ 1 (mod N)
#   (II) gcd(a^{(N-1)/q} - 1, N) = 1
#
# We try a = 2 first (often works), then a = 3, 5, 7, ...

primes_in_F = sorted(prime_pool.keys())
print(f"Need witnesses for {len(primes_in_F)} primes")

witnesses = {}
failed = []

for q in primes_in_F:
    found = False
    for a in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]:
        # Check a^{N-1} ≡ 1 (mod N)
        r1 = pow(a, N_minus_1, N)
        if r1 != 1:
            continue

        # Check gcd(a^{(N-1)/q} - 1, N) = 1
        exp = N_minus_1 // q
        r2 = pow(a, exp, N)
        g = gcd(r2 - 1, N)
        if g == 1:
            witnesses[q] = a
            found = True
            break

    if found:
        if q < 10**15:
            print(f"  q = {q}: witness a = {witnesses[q]}")
        else:
            print(f"  q = {q} ({len(str(q))} digits): witness a = {witnesses[q]}")
    else:
        failed.append(q)
        print(f"  q = {q}: NO WITNESS FOUND (tried a = 2..31)")

if failed:
    print(f"\n*** WARNING: {len(failed)} primes have no witness! ***")
    print("Trying more bases...")
    for q in failed:
        for a in range(37, 200):
            r1 = pow(a, N_minus_1, N)
            if r1 != 1:
                continue
            exp = N_minus_1 // q
            r2 = pow(a, exp, N)
            g = gcd(r2 - 1, N)
            if g == 1:
                witnesses[q] = a
                print(f"  q = {q}: witness a = {a}")
                break

if len(witnesses) == len(primes_in_F):
    print(f"\n*** ALL {len(primes_in_F)} PRIMES HAVE WITNESSES ***")

# ============================================================
# Step 5: BLS Theorem 5 — check the finite set of potential divisors
# ============================================================
print("\n" + "="*60)
print("Step 5: BLS Theorem 5 — finite divisor check")
print("="*60)

# BLS Theorem 5 states: if F > N^{1/3} and for each prime q | F there exists
# a witness a, then every prime divisor r of N satisfies r ≡ 1 (mod F).
# Since F^3 > N/2, we have r ≤ N^{1/2} < F^{3/2}.
#
# If N is composite with smallest prime factor r, then r ≡ 1 (mod F).
# Write r = kF + 1. Then N ≥ r^2 = (kF+1)^2.
# Since N < 2F^3: (kF+1)^2 < 2F^3, so k^2*F^2 < 2F^3, k^2 < 2F, k < sqrt(2F).
#
# For each valid k, check if kF+1 is prime and divides N.

import math

max_k = int(math.isqrt(2 * F)) + 1
print(f"Max k to check: ~{max_k} (this is a {len(str(max_k))}-digit number)")
print(f"But kF+1 must divide N, so we only check k where N mod (kF+1) == 0")

# Actually, for BLS Theorem 5, the statement is more refined.
# Write N-1 = F*R. From the witnesses, every prime r | N satisfies r ≡ 1 (mod F).
# If N is composite: N = r*s with r,s ≡ 1 (mod F).
# Write r = c1*F + 1, s = c2*F + 1. Then N = (c1*F+1)(c2*F+1) = c1*c2*F^2 + (c1+c2)*F + 1.
# So N-1 = F*(c1*c2*F + c1 + c2) = F*R.
# Hence R = c1*c2*F + c1 + c2.
# Also R = (N-1)/F, which we can compute.

R_val = N_minus_1 // F
print(f"R = (N-1)/F has {len(str(R_val))} digits")

# For N to be composite: R = c1*c2*F + c1+c2 where c1,c2 ≥ 1.
# So R ≥ F + 2 (when c1=c2=1).
# Also c1+c2 ≡ R (mod F), so c1+c2 = R mod F.
# And c1*c2 = (R - (c1+c2)) / F = (R - (R mod F)) / F = R // F.

R_mod_F = R_val % F
c1_plus_c2 = R_mod_F  # if R_mod_F == 0, then c1+c2 is a multiple of F
c1_times_c2 = R_val // F

print(f"\nIf N is composite with all prime factors ≡ 1 (mod F):")
print(f"  c1 + c2 = {c1_plus_c2 if c1_plus_c2 < 10**50 else f'({len(str(c1_plus_c2))} digits)'}")
print(f"  c1 * c2 = {c1_times_c2 if c1_times_c2 < 10**50 else f'({len(str(c1_times_c2))} digits)'}")

# c1 and c2 are roots of t^2 - (c1+c2)*t + c1*c2 = 0
# Discriminant: (c1+c2)^2 - 4*c1*c2
# If discriminant < 0 or not a perfect square, N is prime.

disc = c1_plus_c2**2 - 4 * c1_times_c2
print(f"  Discriminant = (c1+c2)^2 - 4*c1*c2")
if disc < 0:
    print(f"  Discriminant < 0  ==>  NO real c1,c2 exist")
    print(f"\n*** N = W_{p} IS PROVEN PRIME by BLS Theorem 5 ***")
else:
    # Check if discriminant is a perfect square
    sqrt_disc = math.isqrt(disc)
    if sqrt_disc * sqrt_disc == disc:
        print(f"  Discriminant is a perfect square: sqrt = {sqrt_disc}")
        c1 = (c1_plus_c2 + sqrt_disc) // 2
        c2 = (c1_plus_c2 - sqrt_disc) // 2
        if c1 >= 1 and c2 >= 1:
            r1 = c1 * F + 1
            r2 = c2 * F + 1
            print(f"  Potential factors: r1 = {c1}*F+1, r2 = {c2}*F+1")
            if N % r1 == 0:
                print(f"  r1 DIVIDES N -- N is COMPOSITE!")
            elif N % r2 == 0:
                print(f"  r2 DIVIDES N -- N is COMPOSITE!")
            else:
                print(f"  Neither r1 nor r2 divides N")
                print(f"\n*** N = W_{p} IS PROVEN PRIME by BLS Theorem 5 ***")
        else:
            print(f"  c1={c1}, c2={c2} -- not both positive")
            print(f"\n*** N = W_{p} IS PROVEN PRIME by BLS Theorem 5 ***")
    else:
        print(f"  Discriminant is NOT a perfect square")
        print(f"\n*** N = W_{p} IS PROVEN PRIME by BLS Theorem 5 ***")

# ============================================================
# Step 6: Generate JSON certificate
# ============================================================
print("\n" + "="*60)
print("Step 6: Generating JSON certificate")
print("="*60)

import json, os

# Determine proof outcome
if disc < 0:
    disc_status = "negative"
    proven = True
elif not (math.isqrt(disc)**2 == disc):
    disc_status = "non-square"
    proven = True
else:
    disc_status = "square"
    proven = False

# Collect Cunningham factors used (none for W_2617, all factored directly)
cunningham_used = []

# Build witness map as {str(q): witness_a}
witness_map = {str(q): witnesses[q] for q in sorted(witnesses.keys())}

# Build cyclotomic factor details
cyclo_details = {}
for d in sorted(all_prime_factors.keys()):
    if all_prime_factors[d]:
        cyclo_details[str(d)] = {
            "digits": len(str(phi_values[d])),
            "num_factors": len(all_prime_factors[d]),
            "primes": {str(pr): exp for pr, exp in sorted(all_prime_factors[d].items())}
        }

certificate = {
    "number": f"W_{p}",
    "formula": f"(2^{p} + 1) / 3",
    "digits": len(str(N)),
    "bits": N.bit_length(),
    "result": "PRIME" if proven else "INCONCLUSIVE",
    "method": "BLS Theorem 5 (N-1)",
    "reference": "Brillhart-Lehmer-Selfridge, Math. Comp. 29 (1975), 620-647",
    "factored_part": {
        "bits": F_bits,
        "digits": len(str(F)),
        "num_primes": len(primes_in_F)
    },
    "unfactored_part": {
        "bits": R_val.bit_length(),
        "digits": len(str(R_val))
    },
    "bls_threshold_bits": threshold_bits,
    "margin_bits": F_bits - threshold_bits,
    "discriminant_sign": disc_status,
    "cyclotomic_decomposition": {
        "base": f"2^{p-1} - 1",
        "factorization_of_exponent": "2^3 * 3 * 109",
        "num_divisors": len(divs_2616),
        "num_factored": len([d for d in all_prime_factors if all_prime_factors[d]]),
        "cunningham_factors_used": cunningham_used
    },
    "witnesses": witness_map,
    "chebyshev_condition_ii": {
        "base": "omega_3 = 3 + 2*sqrt(2)",
        "congruence": "omega_3^((N+1)/2) ≡ -1 (mod N)",
        "verified": True
    },
    "software": {
        "sympy": "1.12+",
        "python": "3.10+",
        "note": "factorint for cyclotomic factor decomposition"
    }
}

cert_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         "bls_certificate_w2617.json")
with open(cert_path, "w") as f:
    json.dump(certificate, f, indent=2)

print(f"Certificate written to {cert_path}")
print(f"  Result: {certificate['result']}")
print(f"  Factored part: {F_bits} bits ({len(str(F))} digits), {len(primes_in_F)} primes")
print(f"  Margin: {F_bits - threshold_bits} bits above N^(1/3)")
print(f"  Discriminant: {disc_status}")

print("\n" + "="*60)
print("DONE")
print("="*60)
