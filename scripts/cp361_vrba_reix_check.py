"""cp361: numerically confirm the Vrba-Reix recurrence test on W_p matches
Condition (II) omega_3^{(W_p+1)/2} == -1 (mod W_p) in Z[sqrt2], for both
prime and composite Wagstaff numbers. De-risks the motivational equivalence
claim added to Paper 02 v3.6 intro (Vrba-Reix == Condition (II)).

Vrba-Reix test: S_0 = seed (mod W_p), S_{i+1} = S_i^2 - 2 (mod W_p),
W_p prime  <=>  S_{p-1} == S_0 (mod W_p).  Seeds tried: 3/2 and 1/4.
"""
import sympy


def wagstaff(p):
    return (pow(2, p) + 1) // 3


# Z[sqrt2] arithmetic mod N: element (a,b) = a + b*sqrt2
def mul(x, y, N):
    a, b = x
    c, d = y
    return ((a * c + 2 * b * d) % N, (a * d + b * c) % N)


def powz(x, e, N):
    r = (1, 0)
    while e:
        if e & 1:
            r = mul(r, x, N)
        x = mul(x, x, N)
        e >>= 1
    return r


def condition_II(p):
    N = wagstaff(p)
    omega3 = (3 % N, 2 % N)
    res = powz(omega3, (N + 1) // 2, N)
    return res == ((N - 1) % N, 0)   # == -1


def vrba_reix(p, seed_num, seed_den):
    N = wagstaff(p)
    S0 = (seed_num * pow(seed_den, -1, N)) % N
    S = S0
    for _ in range(p - 1):
        S = (S * S - 2) % N
    return S == S0


# prime-exponent Wagstaff: prime for p in {3,5,7,11,13,17,19,23,31,43,...},
# composite for p in {29,37,41,47,53,59,...}
PRIME_WP = [5, 7, 11, 13, 17, 19, 23, 31, 43, 61]
COMPOSITE_WP = [29, 37, 41, 47, 53, 59, 67, 71, 73]

print(f"{'p':>4} {'W_p prime?':>11} {'CondII':>7} {'VR(3/2)':>8} {'VR(1/4)':>8}")
ok = True
for p in sorted(PRIME_WP + COMPOSITE_WP):
    N = wagstaff(p)
    is_prime = sympy.isprime(N)
    cII = condition_II(p)
    vr32 = vrba_reix(p, 3, 2)
    vr14 = vrba_reix(p, 1, 4)
    # the equivalence claim: CondII pass <=> W_p prime, and VR matches CondII
    agree = (cII == is_prime)
    print(f"{p:>4} {str(is_prime):>11} {str(cII):>7} {str(vr32):>8} "
          f"{str(vr14):>8}  {'OK' if agree else 'MISMATCH'}")
    if not agree:
        ok = False

print()
print("CONCLUSION: Condition (II) tracks primality of W_p exactly on this "
      "sample" if ok else "CONCLUSION: Condition (II) does NOT track primality!")
# report which VR seed matches Condition II across the whole sample
for sn, sd in [(3, 2), (1, 4)]:
    match = all(vrba_reix(p, sn, sd) == condition_II(p)
                for p in PRIME_WP + COMPOSITE_WP)
    print(f"  VR seed {sn}/{sd} matches Condition (II) on all samples: {match}")
