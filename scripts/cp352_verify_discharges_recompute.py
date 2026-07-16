"""cp352_verify_discharges_recompute.py -- ADVERSARIAL independent recompute of
the witness-exponent discharges X3 -> X2 -> X1 (cp351 Result 4 + cp352 hit).

Fresh code, no reuse of cp351/cp352 helpers. Verifies:

(1) p* = 85684865933, r2 = 7675646348774307083 (cp351 discharge):
    r2 prime (SymPy AND PARI/GP), r2 = 3 mod 8, 2^p* = -1 mod r2 (=> r2 | W_p*),
    ord_{r2}(2) = 2p* exactly (certified against the full factorization of
    r2 - 1), the EXACT order of omega_3 = 3 + 2 sqrt(2) in F_{r2^2} (own
    F_r[sqrt 2] arithmetic, order certified against the full factorization of
    r2 + 1), v_2(ord) = 2, d_{r2} = odd part of ord/4, full factorization of
    d_{r2}, Class-III status (v_2(ord_q(2)) = 1?) of every prime factor q,
    the parity obstruction (some q not in A), the direct condition-(v) test
    2^{p-2} mod 3d, AND a fully direct evaluation of Condition (II) at r2:
    omega_3^{(W_p+1)/2} mod r2 with the exponent reduced modulo r2 + 1.
(2) The same full recompute for p* = 1251488009, r = 33829735778964491
    (the cp352 C-scanner factor, claimed k = 13515805).
(3) The three Platinum witness data lines of Paper B v3.3 sec:inert-survey:
    r prime, r = 3 mod 8, r = 2 p k + 1 with the stated k, r <= 10^12,
    the claimed d equals the exact computed ord(omega_3)/4, and the
    local-pass condition (v) HOLDS at the witness (both routes).

Every primality assertion is re-certified in PARI/GP (APR-CL), and every
exact-order claim is independently re-certified in PARI/GP polmod arithmetic.
"""

import subprocess
import sympy

FAILURES = []


def check(cond, msg):
    tag = "OK  " if cond else "FAIL"
    print(f"    [{tag}] {msg}")
    if not cond:
        FAILURES.append(msg)
    return cond


def v2(n):
    return (n & -n).bit_length() - 1


def exact_order_int(a, r, M, M_fact):
    """Exact multiplicative order of a mod r, given a^M = 1 mod r and the
    complete factorization of M."""
    assert pow(a, M, r) == 1, "a^M != 1: M is not an exponent multiple"
    o = M
    for q in M_fact:
        while o % q == 0 and pow(a, o // q, r) == 1:
            o //= q
    return o


# ---- own F_r[sqrt 2] arithmetic (r = 3 mod 8 => x^2 - 2 irreducible) ----

def x2_mul(u, w, r):
    a, b = u
    c, d = w
    return ((a * c + 2 * b * d) % r, (a * d + b * c) % r)


def x2_pow(u, e, r):
    res = (1, 0)
    while e > 0:
        if e & 1:
            res = x2_mul(res, u, r)
        u = x2_mul(u, u, r)
        e >>= 1
    return res


def exact_order_x2(u, r, M, M_fact):
    assert x2_pow(u, M, r) == (1, 0), "u^M != 1 in F_{r^2}"
    o = M
    for q in M_fact:
        while o % q == 0 and x2_pow(u, o // q, r) == (1, 0):
            o //= q
    return o


def fact_str(f):
    return " * ".join(f"{q}^{e}" if e > 1 else f"{q}"
                      for q, e in sorted(f.items()))


PARI_PRIMES = set()   # integers to certify prime via PARI APR-CL
PARI_ORDS = []        # (r, fourd, [prime factors of fourd]) order certificates
PARI_CONDV = []       # (p, d, expected_bool)


def analyze_inert_factor(label, p, r, expect_pass,
                         claimed_d=None, claimed_k=None, claimed_d_primes=None):
    print("=" * 78)
    print(label)
    print("-" * 78)

    check(sympy.isprime(p), f"p = {p} is prime (SymPy)")
    check(sympy.isprime(r), f"r = {r} is prime (SymPy)")
    PARI_PRIMES.update([p, r])

    check(r % 8 == 3, f"r mod 8 = {r % 8} (inert class, need 3)")

    k = (r - 1) // (2 * p)
    check(r == 2 * p * k + 1, f"r = 2*p*k + 1 exactly, with k = {k}")
    if claimed_k is not None:
        check(k == claimed_k, f"k = {k} matches the stated k = {claimed_k}")

    check(pow(2, p, r) == r - 1,
          "2^p = -1 (mod r)  => r | 2^p + 1, r != 3 => r | W_p "
          "(and r << W_p => W_p composite)")

    # exact ord_r(2), certified against the complete factorization of r-1
    rm1_fact = sympy.factorint(r - 1)
    print(f"    r - 1 = {fact_str(rm1_fact)}")
    PARI_PRIMES.update(rm1_fact)
    o2 = exact_order_int(2, r, r - 1, rm1_fact)
    check(o2 == 2 * p, f"ord_r(2) = {o2} = 2p exactly")

    # 2 must be a non-residue mod r so that F_r[sqrt 2] is a field
    check(pow(2, (r - 1) // 2, r) == r - 1,
          "2^((r-1)/2) = -1 (mod r): 2 is a non-residue, F_r[sqrt2] = F_{r^2}")

    # sanity of the quadratic arithmetic
    check(x2_mul((1, 1), (1, 1), r) == (3 % r, 2 % r),
          "(1 + sqrt2)^2 = 3 + 2 sqrt2 in own arithmetic")
    rp1 = r + 1
    check(x2_pow((1, 1), rp1, r) == (r - 1, 0),
          "alpha^{r+1} = -1 in F_{r^2} (norm/Frobenius identity)")

    # exact order of omega_3, certified against the full factorization of r+1
    rp1_fact = sympy.factorint(rp1)
    print(f"    r + 1 = {fact_str(rp1_fact)}")
    PARI_PRIMES.update(rp1_fact)
    om = (3, 2)
    oo = exact_order_x2(om, r, rp1, rp1_fact)
    vv = v2(oo)
    check(vv == 2,
          f"ord(omega3) = {oo}; v_2(ord) = {vv} = 2 (Lemma v2-class), "
          f"ord = r+1 (full torus): {oo == rp1}")
    if vv != 2:
        print("    !! v_2(ord) != 2: Condition (II) fails trivially; "
              "d_r undefined as ord/4 odd")
        return None

    d = oo // 4
    odd_part = d >> v2(d) if d % 2 == 0 else d
    check(d % 2 == 1 and odd_part == d,
          f"d_r = ord/4 = {d} is odd (= its own odd part)")
    if claimed_d is not None:
        check(d == claimed_d, f"computed d_r = {d} matches claimed d = {claimed_d}")

    d_fact = sympy.factorint(d)
    print(f"    d_r = {fact_str(d_fact)}")
    PARI_PRIMES.update(d_fact)
    if claimed_d_primes is not None:
        check(sorted(d_fact) == sorted(claimed_d_primes) and
              all(e == 1 for e in d_fact.values()),
              f"claimed prime factorization {sorted(claimed_d_primes)} matches")

    # Class III membership of every prime factor of d
    bad = []
    for q in sorted(d_fact):
        oq = exact_order_int(2, q, q - 1, sympy.factorint(q - 1))
        vq = v2(oq)
        in_a = (vq == 1)
        print(f"      q = {q}: ord_q(2) = {oq}, v_2 = {vq}, "
              f"Class III (in A): {in_a}")
        if not in_a:
            bad.append(q)

    # condition (v): d | W_{p-2}  <=>  3d | 2^{p-2}+1  <=>  2^{p-2} = -1 mod 3d
    cv = pow(2, p - 2, 3 * d) == 3 * d - 1
    print(f"    condition (v) [d | W_(p-2)] via 2^(p-2) mod 3d: {cv}")
    PARI_CONDV.append((p, d, cv))

    # DIRECT Condition (II) at r: omega3^{(W_p+1)/2} mod r.
    # (W_p+1)/2 = (2^{p-1}+2)/3; reduce mod r+1 (ord(omega3) | r+1):
    M = 3 * rp1
    c = pow(2, p - 1, M)
    num = (c + 2) % M
    check(num % 3 == 0, "exponent reduction: 3 | (2^{p-1}+2 mod 3(r+1))")
    e0 = num // 3                      # = (W_p+1)/2 mod (r+1)
    val = x2_pow(om, e0, r)
    cond2 = (val == (r - 1, 0))
    print(f"    DIRECT Condition (II) at r: omega3^((W_p+1)/2) mod r "
          f"{'= -1' if cond2 else '!= -1 (value ' + str(val) + ')'}")
    check(cond2 == cv,
          "descent consistency: direct Condition-(II) test agrees with (v)")

    PARI_ORDS.append((r, 4 * d, sorted(set([2] + list(d_fact)))))

    if expect_pass:
        check(cond2, "EXPECTED LOCAL PASS: Condition (II) holds at r")
        check(not bad, "all prime factors of d are Class III "
                       "(necessary for a pass)")
        print("    => LOCAL PASS confirmed: genuine Platinum witness behaviour")
    else:
        check(not cond2, "EXPECTED FAIL: Condition (II) fails at r")
        check(bool(bad),
              f"parity obstruction present: q in {bad} divide d but are not "
              f"Class III => q never divides W_(odd) => d does not divide W_(p-2)")
        print(f"    => Condition (II) FAILS at the factor r of W_{p}")
        print(f"    => W_{p} CANNOT pass Condition (II) globally "
              f"(the congruence mod W_p restricts mod r)")
    print()
    return d


# ============================ the five cases ============================

print("#" * 78)
print("# (1) cp351 discharge of p* = 85684865933 via secondary factor r2")
print("#" * 78)
d1 = analyze_inert_factor(
    "secondary factor r2 = 7675646348774307083 of W_85684865933 "
    "(claimed: FAILS Condition II; d = 3*487*1313423399858711)",
    p=85684865933, r=7675646348774307083, expect_pass=False,
    claimed_d=3 * 487 * 1313423399858711,
    claimed_d_primes=[3, 487, 1313423399858711])

print("#" * 78)
print("# (2) cp352 discharge of p* = 1251488009 via new factor (k = 13515805)")
print("#" * 78)
d2 = analyze_inert_factor(
    "new factor r = 33829735778964491 of W_1251488009 "
    "(claimed: FAILS Condition II; d = 3*41*68759625567001; k = 13515805)",
    p=1251488009, r=33829735778964491, expect_pass=False,
    claimed_d=3 * 41 * 68759625567001,
    claimed_k=13515805,
    claimed_d_primes=[3, 41, 68759625567001])

print("#" * 78)
print("# (3) the three Platinum witness data lines (paper sec:inert-survey)")
print("#" * 78)
w1 = analyze_inert_factor(
    "witness 1: p = 1251488009, r = 2p+1 = 2502976019, claimed d = 41716267",
    p=1251488009, r=2502976019, expect_pass=True,
    claimed_d=41716267, claimed_k=1)
w2 = analyze_inert_factor(
    "witness 2: p = 10916765939, r = 152834723147 (= 2*p*7+1), "
    "claimed d = 38208680787",
    p=10916765939, r=152834723147, expect_pass=True,
    claimed_d=38208680787, claimed_k=7)
w3 = analyze_inert_factor(
    "witness 3: p = 85684865933, r = 2p+1 = 171369731867, claimed d = 57",
    p=85684865933, r=171369731867, expect_pass=True,
    claimed_d=57, claimed_k=1)

for (rw, name) in [(2502976019, "witness 1"), (152834723147, "witness 2"),
                   (171369731867, "witness 3")]:
    check(rw <= 10**12, f"{name} r = {rw} <= 10^12 (inside Platinum range)")
for (rd, name) in [(7675646348774307083, "discharge factor r2"),
                   (33829735778964491, "discharge factor r_new")]:
    check(rd > 10**12, f"{name} = {rd} > 10^12 (outside Platinum range, "
                       f"consistent with Platinum completeness)")

# ===================== PARI/GP independent re-certification ==============

print("#" * 78)
print("# PARI/GP independent re-certification (APR-CL primality + orders)")
print("#" * 78)

lines = []
lines.append('allok = 1;')
for n in sorted(PARI_PRIMES):
    lines.append(f'if(!isprime({n},2), allok = 0; '
                 f'print("PARI: NOT PRIME (APR-CL): ", {n}));')
lines.append('if(allok, print("PARI: all listed integers are prime (APR-CL)"));')
for (r, fourd, qs) in PARI_ORDS:
    qlist = ",".join(str(q) for q in qs)
    lines.append(
        f'r = {r}; w = Mod(Mod(1,r)*(3 + 2*x), x^2 - 2); fd = {fourd}; '
        f'ok = (w^fd == 1); foreach([{qlist}], q, '
        f'if(fd % q == 0 && w^(fd/q) == 1, ok = 0)); '
        f'if(ok, print("PARI: ord(omega3) in F_(", r, "^2) = ", fd, '
        f'" certified"), allok = 0; '
        f'print("PARI: ORDER CERT FAILED at r = ", r));')
for (p, d, expected) in PARI_CONDV:
    lines.append(
        f'got = (Mod(2, 3*{d})^({p}-2) == -1); '
        f'if(got == {1 if expected else 0}, '
        f'print("PARI: cond (v) at p={p}, d={d} -> ", got, " as expected"), '
        f'allok = 0; print("PARI: COND-(V) MISMATCH at p={p}"));')
# order pinning: 2^p = -1 mod r at the five factors
for (p, r) in [(85684865933, 7675646348774307083),
               (1251488009, 33829735778964491),
               (1251488009, 2502976019),
               (10916765939, 152834723147),
               (85684865933, 171369731867)]:
    lines.append(
        f'if(Mod(2,{r})^{p} == -1, '
        f'print("PARI: 2^{p} = -1 mod {r}"), '
        f'allok = 0; print("PARI: DIVISIBILITY FAILED {r}"));')
lines.append('if(allok, print("PARI: ALL CHECKS PASSED"), '
             'print("PARI: SOME CHECKS FAILED"));')
lines.append('quit;')

res = subprocess.run(["gp", "-q"], input="\n".join(lines),
                     capture_output=True, text=True, timeout=300)
print(res.stdout)
if res.stderr.strip():
    print("PARI stderr:", res.stderr)
pari_ok = ("PARI: ALL CHECKS PASSED" in res.stdout)
check(pari_ok, "PARI/GP re-certification: all primality, order, and "
               "condition-(v) checks passed")

# ================================ summary ================================

print("#" * 78)
print("# SUMMARY")
print("#" * 78)
if FAILURES:
    print(f"{len(FAILURES)} CHECK(S) FAILED:")
    for m in FAILURES:
        print("  -", m)
else:
    print("ALL CHECKS PASSED.")
    print()
    print("(1) r2 = 7675646348774307083 is a prime inert factor of "
          "W_85684865933, ord_{r2}(2) = 2p*, d_{r2} = 3*487*1313423399858711;")
    print("    487 and 1313423399858711 are NOT Class III (odd ord_q(2)) "
          "=> parity obstruction => Condition (II) FAILS at r2")
    print("    => W_85684865933 is NOT a global Condition-(II) passer. "
          "DISCHARGED.")
    print("(2) r = 33829735778964491 is a prime inert factor of "
          "W_1251488009 (k = 13515805), ord_r(2) = 2p*, "
          "d_r = 3*41*68759625567001;")
    print("    41 and 68759625567001 are NOT Class III "
          "=> parity obstruction => Condition (II) FAILS at r")
    print("    => W_1251488009 is NOT a global Condition-(II) passer. "
          "DISCHARGED.")
    print("(3) All three Platinum witness data lines check out exactly "
          "(r, k, d, r <= 10^12), and condition (v) HOLDS at each witness")
    print("    (so the witnesses are genuine local passes; "
          "p* = 10916765939 remains undischarged: X1).")
