"""cp352: INDEPENDENT re-verification of the cp351 fixed-d NCT closures at
d in {129, 139, 171, 177} (thm-cp351-nct-closures-129-139-171-177).

Fresh code throughout -- no reuse of cp351 implementations. Checks:

A. Pell machinery from scratch: U, V by their own recurrences; internal
   consistency via U_{2n} = U_n V_n, V_n = 2(U_n + U_{n-1}),
   V_n^2 - 8 U_n^2 = 4(-1)^n.
B. Primitive parts two independent ways:
     (i)  Psi_n = prod_{e|n} U_e^{mu(n/e)}      (Mobius over U; exact int)
     (ii) Psi_{4m} = V_{2m} / prod_{m'|m, m'<m} Psi_{4m'}   (V-route)
   plus full re-derivation checks  V_{2d} = prod_{m|d} Psi_{4m}  and
   U_{4d} = prod_{e|4d} Psi_e  as exact integer identities, and the
   gcd-structure check gcd(Psi_{4d}, U_e) = 1 for every proper e | 4d
   (no intrinsic prime present; certificate-completeness sanity).
C. Claimed factor lists: distinctness + product == Psi_{4d} EXACTLY.
D. Primality of all 15 factors: sympy.isprime (BPSW) AND PARI
   isprime(.,2) (APR-CL) AND isprime(.,1) (Pocklington-Lehmer).
E. Rank of apparition of every factor recomputed from scratch
   (first n with r | U_n); must be exactly 4d (primitive divisor).
F. For every split factor r = 1 mod 8: factor r-1 with PARI, re-verify
   that factorization in Python (product + APR-CL primality of each
   prime), run the standard order descent from r-1 (exact ord_r(2) with
   minimality witnesses), p_r = ord/2, primality of p_r, and -- were any
   p_r prime -- the exact condition d | W_{p_r-2} <=> 2^(p_r-2) = -1
   (mod 3d).
G. Cross-check outcomes against the cp351-claimed exclusion reasons.
H. Context: recompute the parity-unblocked (admissible) odd d in
   (99,200] and the residual open set after these four closures.

Exit: prints PASS/FAIL verdict per d and overall; any assertion failure
means a defect in the claim under test (or in this script -- investigate).
"""

import subprocess
import sys
from math import gcd

import sympy

D_SET = [129, 139, 171, 177]

# The claim under test (factor lists from cp351 / FactorDB):
CLAIMED = {
    129: [75390435818587404514891047039001, 276431598001487177038457370080041],
    139: [6568027, 574593323881,
          197671787019026728856728443573405844788899,
          5764333243492029733503944458218529327657872457],
    171: [8209, 10259, 1186057, 82658931465203, 300531150542299499,
          19269878832774783314803453706733003299],
    177: [6362089, 21015437931105774861319763619250473457,
          490240983802792962966192466615950084458511337],
}
CLAIMED_DIGITS = {129: 65, 139: 106, 171: 83, 177: 89}

FAILURES = []


def check(label, ok):
    print(f"  [{'OK' if ok else 'FAIL'}] {label}")
    if not ok:
        FAILURES.append(label)


def gp_run(script):
    """Run a PARI/GP script, return stdout lines."""
    res = subprocess.run(["gp", "-q"], input=script, capture_output=True,
                         text=True, timeout=600)
    if res.returncode != 0:
        raise RuntimeError(f"gp failed: {res.stderr}")
    return [ln.strip() for ln in res.stdout.splitlines() if ln.strip()]


def gp_factor(n):
    """Factor n with PARI; returns list of (prime, exponent). Verified by caller."""
    lines = gp_run(
        f"fa = factor({n}); "
        f"for(i=1, matsize(fa)[1], print(fa[i,1], \" \", fa[i,2]));"
    )
    return [(int(a), int(b)) for a, b in (ln.split() for ln in lines)]


# ----------------------------------------------------------------------
# A. Pell sequences from scratch
# ----------------------------------------------------------------------
print("=" * 72)
print("A. Pell U/V recurrences + internal identities")
N = max(4 * d for d in D_SET)            # 708
U = [0, 1]
V = [2, 2]
for n in range(2, N + 1):
    U.append(2 * U[-1] + U[-2])
    V.append(2 * V[-1] + V[-2])

check("U[1..8] = 1,2,5,12,29,70,169,408",
      U[1:9] == [1, 2, 5, 12, 29, 70, 169, 408])
check("V[1..8] = 2,6,14,34,82,198,478,1154",
      V[1:9] == [2, 6, 14, 34, 82, 198, 478, 1154])
check("U_{2n} = U_n * V_n  for all n <= 354",
      all(U[2 * n] == U[n] * V[n] for n in range(1, N // 2 + 1)))
check("V_n = 2(U_n + U_{n-1})  for all n <= 708",
      all(V[n] == 2 * (U[n] + U[n - 1]) for n in range(1, N + 1)))
check("V_n^2 - 8 U_n^2 = 4(-1)^n  for all n <= 708",
      all(V[n] ** 2 - 8 * U[n] ** 2 == 4 * (-1) ** n for n in range(N + 1)))

# ----------------------------------------------------------------------
# B. Primitive parts two ways + decomposition identities
# ----------------------------------------------------------------------
print("=" * 72)
print("B. Primitive parts Psi: Mobius-over-U vs V-recursive; decomposition")


def mobius(n):
    f = sympy.factorint(n)
    if any(e > 1 for e in f.values()):
        return 0
    return -1 if len(f) % 2 else 1


def psi_mobius(n):
    """Psi_n = prod_{e|n} U_e^{mu(n/e)} as an exact integer."""
    num = den = 1
    for e in sympy.divisors(n):
        m = mobius(n // e)
        if m == 1:
            num *= U[e]
        elif m == -1:
            den *= U[e]
    assert num % den == 0, f"Psi_{n} not an integer?!"
    return num // den


PSI_U = {}
for d in D_SET:
    for e in sympy.divisors(4 * d):
        if e not in PSI_U:
            PSI_U[e] = psi_mobius(e)

_psi_v_memo = {}


def psi_via_V(m):
    """Psi_{4m} for odd m via the V_{2m} cyclotomic decomposition."""
    if m in _psi_v_memo:
        return _psi_v_memo[m]
    val = V[2 * m]
    for mp in sympy.divisors(m):
        if mp < m:
            q = psi_via_V(mp)
            assert val % q == 0, f"V_{2*m} not divisible by Psi_{4*mp}"
            val //= q
    _psi_v_memo[m] = val
    return val


for d in D_SET:
    odd_divs = [m for m in sympy.divisors(d)]      # d odd => all divisors odd
    check(f"d={d}: all divisors odd", all(m % 2 == 1 for m in odd_divs))
    check(f"d={d}: Psi_{4*d} via U-Mobius == via V-recursion",
          PSI_U[4 * d] == psi_via_V(d))
    check(f"d={d}: V_{2*d} == prod_{{m|{d}}} Psi_{{4m}}  (decomposition)",
          V[2 * d] == sympy.prod(PSI_U[4 * m] for m in odd_divs))
    check(f"d={d}: U_{4*d} == prod_{{e|{4*d}}} Psi_e  (Mobius inversion)",
          U[4 * d] == sympy.prod(PSI_U[e] for e in sympy.divisors(4 * d)))
    # Psi_e | U_e for EVERY e | 4d -- the link used in the soundness proof
    # that a primitive divisor of U_{4d} must divide Psi_{4d}:
    check(f"d={d}: U_e == prod_{{f|e}} Psi_f for every e | {4*d}",
          all(U[e] == sympy.prod(PSI_U[f] for f in sympy.divisors(e))
              for e in sympy.divisors(4 * d)))
    check(f"d={d}: Psi_{4*d} has {CLAIMED_DIGITS[d]} digits",
          len(str(PSI_U[4 * d])) == CLAIMED_DIGITS[d])
    check(f"d={d}: gcd(Psi_{4*d}, U_e) = 1 for all proper e | {4*d} "
          f"(no intrinsic prime)",
          all(gcd(PSI_U[4 * d], U[e]) == 1
              for e in sympy.divisors(4 * d) if e < 4 * d))

# ----------------------------------------------------------------------
# C. Product identity of the claimed factor lists
# ----------------------------------------------------------------------
print("=" * 72)
print("C. Claimed factor lists multiply EXACTLY to Psi_{4d}")
for d in D_SET:
    facs = CLAIMED[d]
    check(f"d={d}: {len(facs)} claimed factors pairwise distinct",
          len(set(facs)) == len(facs))
    check(f"d={d}: product of claimed factors == Psi_{4*d} EXACTLY",
          sympy.prod(facs) == PSI_U[4 * d])

# ----------------------------------------------------------------------
# D. Primality of all 15 factors (sympy BPSW + PARI APR-CL + Pocklington)
# ----------------------------------------------------------------------
print("=" * 72)
print("D. Primality certification of all claimed factors")
ALL_FACTORS = [r for d in D_SET for r in CLAIMED[d]]
check("total claimed factors = 15", len(ALL_FACTORS) == 15)
check("sympy.isprime (BPSW) on all 15",
      all(sympy.isprime(r) for r in ALL_FACTORS))

gp_script = "".join(f"print(isprime({r},2), \" \", isprime({r},1));\n"
                    for r in ALL_FACTORS)
lines = gp_run(gp_script)
aprcl_ok = all(ln.split()[0] == "1" for ln in lines)
pock_ok = all(ln.split()[1] == "1" for ln in lines)
check("PARI isprime(r,2) [APR-CL] = 1 on all 15", aprcl_ok and len(lines) == 15)
check("PARI isprime(r,1) [Pocklington-Lehmer] = 1 on all 15", pock_ok)

# ----------------------------------------------------------------------
# E. Rank of apparition of every factor (fresh modular recurrence)
# ----------------------------------------------------------------------
print("=" * 72)
print("E. Rank of apparition rho(r) recomputed from scratch")
for d in D_SET:
    for r in CLAIMED[d]:
        a, b = 0, 1                      # U_0, U_1 mod r
        rho = None
        for n in range(1, 4 * d + 1):
            if b % r == 0:
                rho = n
                break
            a, b = b, (2 * b + a) % r
        check(f"d={d}: rho({r}) == {4*d} (primitive divisor of U_{4*d})",
              rho == 4 * d)

# ----------------------------------------------------------------------
# F. Split factors: exact order descent, half-order, W-condition
# ----------------------------------------------------------------------
print("=" * 72)
print("F. Split factors r = 1 (mod 8): ord_r(2) from scratch, certificate")
outcomes = {}            # d -> list of (r, reason, ...)
for d in D_SET:
    outcomes[d] = []
    print(f"--- d = {d} ---")
    for r in CLAIMED[d]:
        cls = r % 8
        if cls != 1:
            check(f"d={d}: r={r} = {cls} (mod 8): inert, not a triple candidate",
                  cls == 3)
            continue
        # (F1) factor r-1 with PARI, verify completely in Python
        fac = gp_factor(r - 1)
        prod_check = 1
        for q, a in fac:
            prod_check *= q ** a
        check(f"d={d}: r={r}: PARI factorization of r-1 multiplies back",
              prod_check == r - 1)
        check(f"d={d}: r={r}: each prime of r-1 passes sympy BPSW",
              all(sympy.isprime(q) for q, _ in fac))
        qlines = gp_run("".join(f"print(isprime({q},2));\n" for q, _ in fac))
        check(f"d={d}: r={r}: each prime of r-1 passes PARI APR-CL",
              all(ln == "1" for ln in qlines))

        # (F2) exact multiplicative order descent
        assert pow(2, r - 1, r) == 1     # Fermat (r certified prime)
        o = r - 1
        rem = {q: a for q, a in fac}
        for q in list(rem):
            while rem[q] > 0 and pow(2, o // q, r) == 1:
                o //= q
                rem[q] -= 1
        rem = {q: a for q, a in rem.items() if a > 0}
        check(f"d={d}: r={r}: 2^ord == 1 (mod r)", pow(2, o, r) == 1)
        check(f"d={d}: r={r}: minimality witnesses 2^(ord/q) != 1 for all q|ord",
              all(pow(2, o // q, r) != 1 for q in rem))

        # (F3) the certificate steps
        if o % 2 == 1:
            print(f"    r={r}: ord_r(2)={o} ODD -> ord != 2p, excluded")
            outcomes[d].append((r, "ord_odd"))
            continue
        p_r = o // 2
        pr_fac = dict(rem)
        pr_fac[2] -= 1
        pr_fac = {q: a for q, a in pr_fac.items() if a > 0}
        pr_prime_sympy = sympy.isprime(p_r)
        pr_prime_gp = gp_run(f"print(isprime({p_r}));")[0] == "1"
        check(f"d={d}: r={r}: sympy/PARI agree on primality of p_r",
              pr_prime_sympy == pr_prime_gp)
        # consistency: descent-derived factorization vs primality verdict
        recon = 1
        for q, a in pr_fac.items():
            recon *= q ** a
        check(f"d={d}: r={r}: descent factorization of p_r multiplies back",
              recon == p_r)
        is_one_prime = (len(pr_fac) == 1 and list(pr_fac.values())[0] == 1)
        check(f"d={d}: r={r}: descent factorization consistent with primality",
              is_one_prime == pr_prime_sympy)
        if not pr_prime_sympy:
            fstr = " * ".join(f"{q}^{a}" if a > 1 else f"{q}"
                              for q, a in sorted(pr_fac.items()))
            print(f"    r={r}: ord={o}, p_r=ord/2={p_r} COMPOSITE = {fstr}"
                  f" -> excluded")
            outcomes[d].append((r, "half_composite", p_r, pr_fac))
            continue
        # p_r prime: the W-condition (would be a SURVIVOR if it held)
        cond = (p_r >= 5) and (pow(2, p_r - 2, 3 * d) == 3 * d - 1)
        print(f"    r={r}: ord={o}, p_r={p_r} PRIME; d|W_(p_r-2): {cond}")
        if cond:
            outcomes[d].append((r, "SURVIVOR", p_r))
        else:
            outcomes[d].append((r, "prime_but_W_fails", p_r))

# ----------------------------------------------------------------------
# G. Conclusion per d + cross-check against cp351's claimed reasons
# ----------------------------------------------------------------------
print("=" * 72)
print("G. Certificate conclusion and cp351 cross-check")
EXPECTED = {        # from cp351 checkpoint/node
    129: {"ord_odd": 1, "half_composite": 1},
    139: {"ord_odd": 1, "half_composite": 1},
    171: {"ord_odd": 1, "half_composite": 1},
    177: {"half_composite": 3},
}
for d in D_SET:
    surv = [t for t in outcomes[d] if t[1] == "SURVIVOR"]
    check(f"d={d}: NO surviving factor -> no compatible triple at d={d}",
          len(surv) == 0)
    counts = {}
    for t in outcomes[d]:
        counts[t[1]] = counts.get(t[1], 0) + 1
    check(f"d={d}: exclusion-reason profile matches cp351 claim {EXPECTED[d]}",
          counts == EXPECTED[d])

# specific factor-level claims from the cp351 node text
hc = {d: [t for t in outcomes[d] if t[1] == "half_composite"] for d in D_SET}
check("d=129: composite half-order divisible by 3 * 13^2 * 43",
      any(t[2] % (3 * 13 ** 2 * 43) == 0 for t in hc[129]))
check("d=139: composite half-order divisible by 3 * 139 * 157",
      any(t[2] % (3 * 139 * 157) == 0 for t in hc[139]))
check("d=171: composite half-order == 1026 = 2 * 3^3 * 19",
      any(t[2] == 1026 for t in hc[171]))
check("d=177: exactly three split factors, all half-orders composite",
      len(hc[177]) == 3)

# ----------------------------------------------------------------------
# H. Context: parity-unblocked (admissible) d in (99, 200]
# ----------------------------------------------------------------------
print("=" * 72)
print("H. Admissible (parity-unblocked) odd d in (99,200] and residual set")


def in_A(q):
    """q in A  <=>  v_2(ord_q(2)) = 1."""
    o = sympy.n_order(2, q)
    return o % 2 == 0 and (o // 2) % 2 == 1


admissible = [d for d in range(101, 200, 2)
              if all(in_A(q) for q in sympy.factorint(d))]
print(f"  admissible d in (99,200]: {admissible}")
check("admissible set == {107,121,129,131,139,163,171,177,179} (paper list)",
      admissible == [107, 121, 129, 131, 139, 163, 171, 177, 179])
check("each target d in {129,139,171,177} is admissible (non-vacuous closure)",
      all(d in admissible for d in D_SET))
residual = [d for d in admissible if d not in (107, 121) and d not in D_SET]
check("residual open parity-unblocked set == {131, 163, 179}",
      residual == [131, 163, 179])

# ----------------------------------------------------------------------
print("=" * 72)
if FAILURES:
    print(f"OVERALL: {len(FAILURES)} FAILURE(S):")
    for f in FAILURES:
        print(f"  - {f}")
    sys.exit(1)
print("OVERALL: ALL CHECKS PASS -- NCT closures at d = 129, 139, 171, 177 "
      "independently reconfirmed.")
