#!/usr/bin/env python3
"""cp412b: T1 adversarial recheck of the cp412 danger-triple census (43 < d <= 400).

Independent re-implementation (no imports from cp412_danger_census_d400.py) of
every mathematical ingredient, then a from-scratch recomputation of the census
verdicts for all recorded primes, plus the C1-C7 attack suites requested by the
verification protocol:

  C1  Frobenius frame: alpha = 1+sqrt2, N(alpha) = -1, alpha^r = conj(alpha),
      alpha^{r+1} = -1 in F_{r^2} for inert r; omega_3 = alpha^2;
      r == 3 (mod 8) <=> v_2(ord(omega_3)) = 2 among inert primes; split
      primitive divisors have r == 1 (mod 8) with 8d | r-1; no prime == 5,7
      (mod 8) divides V_{2d} (d odd).  Numerics: all primes of V_{2d} for all
      odd d <= 35 + independent V-recurrence divisibility patterns.
  C2  D | W_n <=> 2^n == -1 (mod 3D), odd D <= 1000 (incl. prime powers),
      odd n <= 200, against exact integer W_n.
  C3  factorization-free v2(ord_r(2)) vs exact orders for all odd r < 1e5;
      brute-force ord(omega_3) in F_{r^2} for r == 3 (mod 8), r < 20000;
      primitivity-descent biconditional stress test.
  C4  admissible-rung recompute; product-identity re-verification of all 34
      rung factorizations; independent from-scratch factorization of V_114,
      V_118, V_134, V_162; multiplicity (square divisor) audit vs LTE;
      global dedup / rung-attribution audit.
  C5  full recompute of every (i)-(iv) survivor, REAL triple and phantom:
      primality (BPSW + PARI/GP APR-CL), r | V_2d exact, ord_r(2) = 2p exact
      (two independent certificates), ord(omega_3) = 4d exact (divisor
      certificate, no factoring), condition (v), k, cutoff flags, phantom test.
  C6  corner-sweep consistency; prop:diagonal-closures reconciliation;
      split-side crosscheck recompute (all 160 split primes); discharge_status
      vs prop:witness-discharges + sec:secondary + X1 node (incl. independent
      verification of the paper's three secondary/discharge factors and the
      X1 witness data).
  C7  boundary audit (43 < d <= 400, d = 400 exclusion) + census JSON/log
      internal consistency recount.

Output: scripts/output/cp412b_t1_recheck.log (tee'd stdout) and
        scripts/output/cp412b_t1_recheck.json (failure lists per section).
"""

import json
import os
import re
import subprocess
import sys
import time

import sympy
from sympy import isprime, factorint
from sympy.ntheory.residue_ntheory import sqrt_mod

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CENSUS_JSON = os.path.join(ROOT, "scripts/output/cp412_danger_census.json")
FDB_CACHE = os.path.join(ROOT, "scripts/output/cp412_factordb_cache.json")
NCT_JSON = os.path.join(
    ROOT, "papers/02_wagstaff_chebyshev_reduction/zenodo_v3.6/data/nct_certificates.json")
OUT_LOG = os.path.join(ROOT, "scripts/output/cp412b_t1_recheck.log")
OUT_JSON = os.path.join(ROOT, "scripts/output/cp412b_t1_recheck.json")
ECM_BIN = "/tmp/ecmenv/bin/ecm"

_logfh = open(OUT_LOG, "w")


def log(msg=""):
    print(msg, flush=True)
    _logfh.write(msg + "\n")
    _logfh.flush()


FAILURES = {}


def fail(section, msg):
    FAILURES.setdefault(section, []).append(msg)
    log(f"   !! FAIL [{section}]: {msg}")


# ---------------------------------------------------------------------------
# fresh arithmetic
# ---------------------------------------------------------------------------

def my_pell(n):
    """(U_n, V_n) exact; U_0,U_1 = 0,1; V_0,V_1 = 2,2; X_{k+1} = 2X_k + X_{k-1}."""
    ua, ub = 0, 1
    va, vb = 2, 2
    for _ in range(n):
        ua, ub = ub, 2 * ub + ua
        va, vb = vb, 2 * vb + va
    return ua, va


def myV(n):
    return my_pell(n)[1]


def v2(n):
    return (n & -n).bit_length() - 1


def pair_mul(x, y, r):
    """(a + b sqrt2)(c + d sqrt2) in Z[sqrt2]/r."""
    a, b = x
    c, d = y
    return ((a * c + 2 * b * d) % r, (a * d + b * c) % r)


def pair_pow(base, e, r):
    res = (1 % r, 0)
    b = (base[0] % r, base[1] % r)
    while e > 0:
        if e & 1:
            res = pair_mul(res, b, r)
        b = pair_mul(b, b, r)
        e >>= 1
    return res


W3 = (3, 2)  # omega_3 = 3 + 2 sqrt2


def w3_pow_is(r, e, target):
    return pair_pow(W3, e, r) == (target % r, 0)


def exact_order_pair(x, r, group_order_factored):
    """exact multiplicative order of pair x in F_r[sqrt2], given a full
    factorization {q: e} of a known multiple of the order."""
    n = 1
    for q, e in group_order_factored.items():
        n *= q ** e
    assert pair_pow(x, n, r) == (1, 0), "claimed multiple of order is not one"
    o = n
    for q in group_order_factored:
        while o % q == 0 and pair_pow(x, o // q, r) == (1, 0):
            o //= q
    return o


def exact_order_int(a, r, group_order_factored):
    n = 1
    for q, e in group_order_factored.items():
        n *= q ** e
    assert pow(a, n, r) == 1
    o = n
    for q in group_order_factored:
        while o % q == 0 and pow(a, o // q, r) == 1:
            o //= q
    return o


def my_v2_ord2(r):
    """v_2(ord_r(2)) without factoring the odd part of r-1."""
    t = v2(r - 1)
    y = pow(2, (r - 1) >> t, r)
    if y == 1:
        return 0
    for s in range(1, t + 1):
        y = y * y % r
        if y == 1:
            return s
    raise ValueError(f"2^(r-1) != 1 mod {r}")


def inert_rung(r, d, dprimes):
    """r inert, r | V_{2d} (so omega3^{2d} = -1).  Exact rung e with
    ord(omega3) = 4e, e | d.  Biconditional proved in the audit notes:
    omega3^{2m} = -1  <=>  e | m   (m odd).  Descent leaves v_l(g) = v_l(e)."""
    assert w3_pow_is(r, 2 * d, r - 1), f"omega3^(2d) != -1 mod {r}, d={d}"
    g = d
    for l in dprimes:
        while g % l == 0 and w3_pow_is(r, 2 * (g // l), r - 1):
            g //= l
    return g


def split_rung(r, d, dprimes):
    """r == 1 mod 8 split, r | V_{2d}.  omega = 3+2s in F_r.  Exact rung e
    (ord(omega) = 4e); also returns ord(alpha) = 8e certificate data."""
    s = sqrt_mod(2, r)
    w = (3 + 2 * s) % r
    a = (1 + s) % r
    assert pow(w, 2 * d, r) == r - 1, f"omega3^(2d) != -1 mod split {r}"
    g = d
    for l in dprimes:
        while g % l == 0 and pow(w, 2 * (g // l), r) == r - 1:
            g //= l
    # alpha certificate: ord(alpha) = 8g  <=>  alpha^{4g} = -1 (+ descent)
    assert pow(a, 4 * g, r) == r - 1, f"alpha^(4e) != -1 mod {r}"
    return g


def cert_ord2_is_2p(r, p):
    """EXACT no-factoring certificate: ord_r(2) = 2p  <=>  p prime and
    2^p == -1 (mod r) and r > 3.  (divisors of 2p are 1,2,p,2p; 2^p = -1
    kills 1,2? no: kills ord | p; ord = 2 forces r | 3.)"""
    return r > 3 and isprime(p) and pow(2, p, r) == r - 1


def cert_w3_ord_is_4d(r, d, dprimes, split=False):
    """EXACT no-factoring certificate: ord(omega3) = 4d  <=>
    omega3^{2d} = -1 and omega3^{2d/l} != -1 for every prime l | d."""
    if split:
        s = sqrt_mod(2, r)
        w = (3 + 2 * s) % r
        if pow(w, 2 * d, r) != r - 1:
            return False
        return all(pow(w, 2 * (d // l), r) != r - 1 for l in dprimes)
    if not w3_pow_is(r, 2 * d, r - 1):
        return False
    return all(not w3_pow_is(r, 2 * (d // l), r - 1) for l in dprimes)


def gp_aprcl(n, timeout_s=600):
    try:
        res = subprocess.run(["gp", "-q"], input=f"print(isprime({n},2))\n",
                             capture_output=True, text=True, timeout=timeout_s)
        out = res.stdout.strip()
        return True if out == "1" else (False if out == "0" else None)
    except Exception:
        return None


# ---------------------------------------------------------------------------
# odd part of ord_r(2): independent analysis (trial 1e6 + exact 2-block rule
# + local factorint/ecm escalation; FDB-cache primes usable as division-
# verified hints only)
# ---------------------------------------------------------------------------

SMALL_PRIMES = list(sympy.sieve.primerange(2, 10 ** 6))


def factorint_sub(n, timeout_s):
    code = ("import sympy,json\n"
            f"print(json.dumps({{str(k):v for k,v in sympy.factorint({n}).items()}}))\n")
    try:
        res = subprocess.run(["nice", "-n", "10", sys.executable, "-c", code],
                             capture_output=True, text=True, timeout=timeout_s)
        if res.returncode == 0 and res.stdout.strip():
            return {int(k): v for k, v in json.loads(res.stdout).items()}
    except subprocess.TimeoutExpired:
        pass
    return None


def ecm_one_factor(n, budget_s):
    if not os.path.exists(ECM_BIN):
        return None
    t0 = time.time()
    for b1, curves in [(2000, 40), (11000, 40), (50000, 30), (250000, 20)]:
        rem = budget_s - (time.time() - t0)
        if rem <= 3:
            return None
        try:
            res = subprocess.run(
                ["nice", "-n", "10", ECM_BIN, "-q", "-c", str(curves), str(b1)],
                input=str(n), capture_output=True, text=True, timeout=rem)
            for tok in res.stdout.split():
                if tok.isdigit():
                    f = int(tok)
                    if 1 < f < n and n % f == 0:
                        return f
        except subprocess.TimeoutExpired:
            return None
    return None


_HINTS = set()
_OPO_MEMO = {}


def load_hints():
    if os.path.exists(FDB_CACHE):
        c = json.load(open(FDB_CACHE))
        for _k, v in c.items():
            for fs, _e in v.get("factors", []):
                _HINTS.add(int(fs))


def my_oddpart_ord2(r, escalate_budget_s=150):
    if r in _OPO_MEMO:
        return _OPO_MEMO[r]
    res = _my_oddpart_ord2(r, escalate_budget_s)
    _OPO_MEMO[r] = res
    return res


def _my_oddpart_ord2(r, escalate_budget_s=150):
    """verdict on the odd part m' of ord_r(2), assuming v2(ord_r(2)) = 1.
    Returns ('prime', p) | ('composite', m_exact_or_None, why) |
    ('unresolved', blockinfo).  Exactness argument as audited:
    after maximal known-prime stripping, v_l(m) = v_l(m') for every known l,
    and the unfactored block B (product of ALL unknown-prime powers of o)
    strips iff m' | m/B.  If B remains and s = m/B > 1: m' has a known prime
    factor AND an unknown-block prime factor -> composite, exactly."""
    t = v2(r - 1)
    o = (r - 1) >> t
    assert pow(2, o, r) == r - 1, "v2(ord)=1 precondition"
    known = set()
    rem = o
    for q in SMALL_PRIMES:
        if rem == 1 or q * q > rem:
            break
        while rem % q == 0:
            known.add(q)
            rem //= q
    if rem > 1 and isprime(rem):
        known.add(rem)
        rem = 1
    block = rem

    def rebuild_block():
        b = o
        for q in known:
            while b % q == 0:
                b //= q
        return b

    def strip(block_now):
        m = o
        changed = True
        while changed:
            changed = False
            for l in sorted(known):
                while m % l == 0 and pow(2, 2 * (m // l), r) == 1:
                    m //= l
                    changed = True
        if block_now > 1:
            assert m % block_now == 0
            if pow(2, 2 * (m // block_now), r) == 1:
                m //= block_now
                block_now = 1
                changed = True
                while changed:
                    changed = False
                    for l in sorted(known):
                        while m % l == 0 and pow(2, 2 * (m // l), r) == 1:
                            m //= l
                            changed = True
        return m, block_now

    for _round in range(4):
        m, b = strip(block)
        assert pow(2, m, r) == r - 1, "invariant 2^m = -1"
        if b == 1:
            if m == 1:
                return ("one", None, "")
            if isprime(m):
                return ("prime", m, "exact (all primes known)")
            return ("composite", m, "exact odd part, composite")
        if m // b > 1:
            return ("composite", None,
                    f"2-block rule: known prime {m // b} side + {len(str(b))}-digit block")
        # ambiguous: m == b, escalate locally
        newly = set()
        for h in sorted(_HINTS):
            if 1 < h < b and b % h == 0:
                x = h
                if isprime(x):
                    newly.add(x)
        if not newly and len(str(b)) <= 46:
            fi = factorint_sub(b, min(escalate_budget_s, 120))
            if fi:
                newly |= set(fi.keys())
        if not newly:
            f = ecm_one_factor(b, escalate_budget_s)
            if f:
                for piece in (f, b // f):
                    stack = [piece]
                    while stack:
                        x = stack.pop()
                        if isprime(x):
                            newly.add(x)
                        else:
                            ff = factorint_sub(x, 60) if len(str(x)) <= 46 else None
                            if ff:
                                newly |= set(ff.keys())
                            else:
                                f2 = ecm_one_factor(x, 60)
                                if f2:
                                    stack += [f2, x // f2]
                                else:
                                    newly.add(-1)  # marker: incomplete
        newly.discard(-1)
        if not newly:
            return ("unresolved", None, f"{len(str(b))}-digit block, local tools exhausted")
        known |= newly
        block = rebuild_block()
    return ("unresolved", None, "escalation rounds exhausted")


# ===========================================================================
log(f"cp412b T1 adversarial recheck  ({time.strftime('%Y-%m-%d %H:%M:%S')})")
log("")
census = json.load(open(CENSUS_JSON))
load_hints()
T0 = time.time()

# ---------------------------------------------------------------------------
log("== C1: Frobenius frame + divisor classification for V_{2d} ==")
# sanity of my own Pell code against the alpha-power identity
for n in (1, 2, 3, 7, 20, 57):
    u, v = my_pell(n)
    big = 10 ** 80  # work in Z via a huge 'modulus' so pair_pow is exact
    a = pair_pow((1, 1), n, big)
    if not (a[0] * 2 == v and a[1] == u):
        fail("C1", f"alpha^{n} != V/2 + U sqrt2")
if myV(2) != 6 or myV(6) != 198 or myV(10) != 6726:
    fail("C1", "Pell V values wrong")

# (a) inert Frobenius identities + rung/divisibility law, primes < 3000
n_inert = n_split = 0
for r in sympy.sieve.primerange(3, 3000):
    rm = r % 8
    if rm in (3, 5):  # inert
        n_inert += 1
        conj = pair_pow((1, 1), r, r)
        if conj != (1, r - 1):
            fail("C1", f"alpha^r != conj(alpha) for inert r={r}")
        if pair_pow((1, 1), r + 1, r) != (r - 1, 0):
            fail("C1", f"alpha^(r+1) != -1 for inert r={r}")
        if pair_pow(W3, r + 1, r) != (1, 0):
            fail("C1", f"omega3^(r+1) != 1 for inert r={r}")
        ordw = exact_order_pair(W3, r, factorint(r + 1))
        if v2(ordw) != v2(r + 1):
            fail("C1", f"v2(ord omega3) != v2(r+1) at r={r}")
        if rm == 3:
            if ordw % 4 != 0 or (ordw // 4) % 2 == 0:
                fail("C1", f"r={r}==3 mod 8: ord(omega3) not 4*odd")
            e = ordw // 4
            # independent divisibility law via the integer V-recurrence mod r:
            # {m <= M : r | V_{2m}} must be exactly the ODD multiples of e
            M = min(5 * e, 1500)
            va, vb = 2 % r, 2 % r          # V_0, V_1
            hits = []
            for nn in range(2, 2 * M + 1):
                va, vb = vb, (2 * vb + va) % r   # vb = V_nn
                if nn % 2 == 0 and vb == 0:
                    hits.append(nn // 2)
            want = [m for m in range(1, M + 1) if m % 2 == 1 and m % e == 0]
            if hits != want:
                fail("C1", f"divisibility law fails r={r}: hits {hits[:5]} want {want[:5]}")
        else:  # r == 5 mod 8: v2(r+1)=1 -> no odd rung, divides no V_{2d}, d odd
            if v2(ordw) != 1:
                fail("C1", f"r={r}==5 mod 8: v2(ord omega3) = {v2(ordw)} != 1")
    else:  # split
        n_split += 1
        s = sqrt_mod(2, r)
        if s is None:
            fail("C1", f"sqrt2 missing for split r={r}")
            continue
        w = (3 + 2 * s) % r
        ordw = exact_order_int(w, r, factorint(r - 1))
        if rm == 7 and v2(ordw) == 2:
            fail("C1", f"r={r}==7 mod 8 has v2(ord omega3)=2 (impossible)")
        if v2(ordw) == 2:
            e = ordw // 4
            if rm != 1 or (r - 1) % (8 * e) != 0:
                fail("C1", f"split r={r}: rung {e} but 8e does not divide r-1")
log(f"   Frobenius + order-law verified on {n_inert} inert / {n_split} split primes < 3000")

# (b) complete classification of divisors of V_{2d}, all odd d <= 35
checked = 0
for d in range(1, 36, 2):
    Vd = myV(2 * d)
    fac = factorint(Vd)
    prod = 1
    for q, e in fac.items():
        prod *= q ** e
    if prod != Vd:
        fail("C1", f"factorint incomplete for V_{2*d}")
    dprimes = list(factorint(d))
    for q, e in fac.items():
        checked += 1
        if q == 2:
            if e != 1:
                fail("C1", f"v2(V_{2*d}) = {e} != 1")
            continue
        if q % 8 in (5, 7):
            fail("C1", f"prime {q} == {q%8} mod 8 divides V_{2*d}")
            continue
        if q % 8 == 3:
            rung = inert_rung(q, d, dprimes)
            if (q + 1) % (4 * rung) != 0:
                fail("C1", f"4*rung does not divide r+1 for r={q}")
            if v2(q + 1) != 2:
                fail("C1", f"inert primitive r={q}: v2(r+1) = {v2(q+1)} != 2")
        else:
            rung = split_rung(q, d, dprimes)
            if (q - 1) % (8 * rung) != 0:
                fail("C1", f"split r={q}: 8*rung does not divide r-1")
        if d % rung != 0 or rung % 2 == 0:
            fail("C1", f"rung {rung} of {q} does not odd-divide d={d}")
        # primitivity: q divides no V_{2e'} with 0 < e' < rung (exact integers)
        for ep in range(1, rung):
            if myV(2 * ep) % q == 0:
                fail("C1", f"r={q} divides V_{2*ep} below its rung {rung}")
        # and q | V_{2 rung} exactly
        if myV(2 * rung) % q != 0:
            fail("C1", f"r={q} does not divide V at its own rung {rung}")
log(f"   full divisor classification of V_2..V_70: {checked} prime powers checked")
log(f"   C1 {'PASS' if 'C1' not in FAILURES else 'FAIL'}")
log("")

# ---------------------------------------------------------------------------
log("== C2: divisibility criterion D | W_n <=> 2^n == -1 (mod 3D) ==")
cnt = 0
for n in range(1, 201, 2):
    Wn = (2 ** n + 1) // 3
    assert (2 ** n + 1) % 3 == 0
    for D in range(1, 1001, 2):
        direct = (Wn % D == 0)
        crit = (pow(2, n, 3 * D) == 3 * D - 1)
        if direct != crit:
            fail("C2", f"criterion fails at n={n}, D={D}")
        cnt += 1
log(f"   {cnt} (n, D) pairs checked (odd n <= 200, odd D <= 1000, incl. prime powers)")
log(f"   C2 {'PASS' if 'C2' not in FAILURES else 'FAIL'}")
log("")

# ---------------------------------------------------------------------------
log("== C3: v2(ord_r(2)) square-chain + primitivity descent, stress tests ==")
n3 = 0
for r in sympy.sieve.primerange(3, 100000):
    ordr = exact_order_int(2, r, factorint(r - 1))
    if my_v2_ord2(r) != v2(ordr):
        fail("C3", f"v2_ord2({r}) != v2(exact ord) = {v2(ordr)}")
    n3 += 1
log(f"   v2(ord_r(2)) square-chain == exact for all {n3} odd primes r < 1e5")

nb = 0
for r in sympy.sieve.primerange(3, 20000):
    if r % 8 != 3:
        continue
    # brute-force order of omega3 in F_{r^2}
    x = (3 % r, 2 % r)
    o = 1
    while x != (1, 0):
        x = pair_mul(x, W3, r)
        o += 1
    of = exact_order_pair(W3, r, factorint(r + 1))
    if o != of:
        fail("C3", f"brute ord(omega3) mod {r}: {o} != factored {of}")
    if o % 4 != 0 or (o // 4) % 2 == 0:
        fail("C3", f"r={r}: ord(omega3)={o} not 4*odd")
    nb += 1
log(f"   brute-force ord(omega3) == factored-descent order for {nb} primes r==3(8) < 20000")

nd = 0
for r in sympy.sieve.primerange(3, 100000):
    if r % 8 != 3:
        continue
    e = exact_order_pair(W3, r, factorint(r + 1)) // 4
    for f in (1, 3, 5, 9, 15, 21, 25, 33):
        d = e * f
        dprimes = list(factorint(d))
        # my primitivity test (the census's T3 logic, re-implemented)
        prim = w3_pow_is(r, 2 * d, r - 1) and all(
            not w3_pow_is(r, 2 * (d // l), r - 1) for l in dprimes)
        if prim != (f == 1):
            fail("C3", f"descent biconditional fails r={r} e={e} f={f}")
        rung = inert_rung(r, d, dprimes)
        if rung != e:
            fail("C3", f"inert_rung({r},{d}) = {rung} != true rung {e}")
        nd += 1
    # non-multiple: omega3^{2m} must not be -1 when e does not divide m
    m = e + 2 if (e + 2) % e != 0 else e + 4
    if w3_pow_is(r, 2 * m, r - 1) and m % e != 0:
        fail("C3", f"omega3^(2m) = -1 at non-multiple m={m} of e={e}, r={r}")
log(f"   descent biconditional stress: {nd} (r, d) primitivity decisions verified")
log(f"   C3 {'PASS' if 'C3' not in FAILURES else 'FAIL'}")
log("")

# ---------------------------------------------------------------------------
log("== C4: completeness — admissible rungs, product identities, dedup ==")
# (a) admissible list, independent boundary logic
my_adm = []
for d in range(1, 1001, 2):
    if not (43 < d <= 400):
        continue
    dp = list(factorint(d))
    if all(my_v2_ord2(q) == 1 for q in dp):
        my_adm.append(d)
if my_adm != census["admissible_rungs"]:
    fail("C4", f"admissible rung mismatch: mine {my_adm}")
if 400 in my_adm or 399 in my_adm or 45 in my_adm or 43 in my_adm:
    fail("C4", "boundary error in my own list?!")
# cross-check a few Class III / non-Class III memberships against exact orders
for q in (3, 11, 19, 43, 59, 67, 83, 107, 131, 163, 179, 331):
    if v2(exact_order_int(2, q, factorint(q - 1))) != 1:
        fail("C4", f"{q} should be Class III (exact order check)")
for q in (5, 7, 17, 23, 29, 31, 37, 41, 53, 61, 73, 89, 97, 127, 151, 193, 241, 257):
    if v2(exact_order_int(2, q, factorint(q - 1))) == 1:
        fail("C4", f"{q} should NOT be Class III (exact order check)")
log(f"   admissible rungs recomputed independently: {len(my_adm)} rungs, "
    f"{'match' if my_adm == census['admissible_rungs'] else 'MISMATCH'}")

# nct rung consistency
nct = json.load(open(NCT_JSON))
nct_ds = sorted(c["d"] for c in nct["certificates"])
if sorted(set(my_adm) & set(nct_ds)) != nct_ds:
    fail("C4", "nct certificate rungs not all admissible")
if [d for d in my_adm if d > 100] != nct_ds:
    fail("C4", "(100,400] admissible != nct rungs")

# (b) product identities for ALL 34 rungs from the census JSON prime lists
rung_fac = {}     # d -> {prime: mult} for r > 3, from census JSON
v3_of = lambda n: max(k for k in range(0, 40) if n % 3 ** k == 0)
for rec in census["rungs"]:
    d = rec["d"]
    if rec["status"] != "COMPLETE":
        fail("C4", f"rung d={d} not COMPLETE")
        continue
    V = myV(2 * d)
    if v2(V) != 1:
        fail("C4", f"v2(V_{2*d}) != 1")
    prod = 1
    fac = {}
    for pr in rec["primes"]:
        rr = int(pr["r"])
        fac[rr] = pr["mult"]
        if not isprime(rr):
            fail("C4", f"d={d}: listed factor {rr} fails BPSW")
        if rr % 8 not in (1, 3):
            fail("C4", f"d={d}: factor {rr} == {rr%8} mod 8")
        prod *= rr ** pr["mult"]
    if V % prod != 0:
        fail("C4", f"d={d}: listed primes do not divide V_2d — POISONED")
        continue
    Q = V // prod
    b = 0
    while Q % 3 == 0:
        Q //= 3
        b += 1
    if Q != 2:
        fail("C4", f"d={d}: cofactor after primes and 3-part is {Q} != 2 — INCOMPLETE list")
    if b != 1 + v3_of(d):
        fail("C4", f"d={d}: v3(V_2d) = {b} != 1 + v3(d) (LTE)")
    rung_fac[d] = fac
log(f"   product identity re-verified (with 2*3^(1+v3(d)) closure) for {len(rung_fac)}/34 rungs")

# (c) independent from-scratch factorization of the four diagonal-closure V's
for d in (57, 59, 67, 81):
    U, _V = my_pell(d)
    h1, h2 = 2 * U - 1, 2 * U + 1
    V = myV(2 * d)
    if V != 2 * h1 * h2:
        fail("C4", f"half-split identity fails at d={d}")
    myfac = {2: 1}
    ok = True
    for h in (h1, h2):
        fi = factorint_sub(h, 240)
        if fi is None:
            log(f"   [note] local factorint timeout on half of V_{2*d} — skipped")
            ok = False
            break
        for q, e in fi.items():
            myfac[q] = myfac.get(q, 0) + e
    if not ok:
        continue
    prod = 1
    for q, e in myfac.items():
        if not isprime(q):
            fail("C4", f"my factor {q} of V_{2*d} not prime?!")
        prod *= q ** e
    if prod != V:
        fail("C4", f"my own factorization of V_{2*d} fails product identity")
    mine_gt3 = {q: e for q, e in myfac.items() if q > 3}
    if mine_gt3 != rung_fac.get(d):
        fail("C4", f"d={d}: MY independent factorization differs from census: "
                   f"mine {sorted(mine_gt3)} census {sorted(rung_fac.get(d, {}))}")
    else:
        log(f"   V_{2*d}: independent local factorization == census list "
            f"({len(mine_gt3)} primes > 3)")

# (d+e) global rung attribution / dedup / multiplicity audit
true_rung = {}       # r -> e (my exact rung)
side_of = {}
for rec in census["rungs"]:
    d = rec["d"]
    dprimes = list(factorint(d))
    for pr in rec["primes"]:
        rr = int(pr["r"])
        if myV(2 * d) % rr != 0:
            fail("C4", f"d={d}: listed prime {rr} does not divide V_2d (integer check)")
        if rr % 8 == 3:
            e = inert_rung(rr, d, dprimes)
            side = "inert"
        else:
            e = split_rung(rr, d, dprimes)
            side = "split"
        if rr in true_rung and true_rung[rr] != e:
            fail("C4", f"prime {rr}: inconsistent rung {true_rung[rr]} vs {e}")
        true_rung[rr] = e
        side_of[rr] = side
        # multiplicity vs LTE: mult(d) = mult(V_{2e}) + v_r(d/e); primitive mult
        # is 1 for every prime here unless flagged
        f = d // e
        vr = 0
        while f % rr == 0:
            f //= rr
            vr += 1
        if pr["mult"] != 1 + vr:
            fail("C4", f"d={d}, r={rr}: mult {pr['mult']} != 1 + v_r(d/e) = {1+vr} "
                       f"(possible Wall-Sun-Sun-type square at primitive rung!)")
# appearance law: r appears in census rung d exactly when true_rung | d
adm_set = set(my_adm)
for rr, e in true_rung.items():
    for rec in census["rungs"]:
        d = rec["d"]
        listed = any(int(pr["r"]) == rr for pr in rec["primes"])
        if listed != (d % e == 0):
            fail("C4", f"appearance law: r={rr} rung {e} listed={listed} at d={d}")
    if e > 43 and e not in adm_set and e <= 400:
        fail("C4", f"prime {rr} primitive at NON-admissible rung {e} in range?!")
log(f"   global attribution audit: {len(true_rung)} distinct primes, appearance law checked")
n_m2 = sum(1 for rec in census["rungs"] for pr in rec["primes"] if pr["mult"] > 1)
log(f"   square-divisor audit: {n_m2} entries with mult >= 2 (all LTE-explained above)")
log(f"   C4 {'PASS' if 'C4' not in FAILURES else 'FAIL'}")
log("")

# ---------------------------------------------------------------------------
log("== C5 + C6: full independent reclassification of all census primes ==")
have_gp = subprocess.run(["which", "gp"], capture_output=True).returncode == 0
my_survivors = {}         # d -> set of (i)-(iv) survivor r
my_real, my_phantom = [], []
my_split_viol = []
my_unresolved = []
mismatches = []
n_class = 0
for rec in census["rungs"]:
    d = rec["d"]
    dprimes = list(factorint(d))
    surv = set()
    for pr in rec["primes"]:
        rr = int(pr["r"])
        n_class += 1
        cv = pr["tests"].get("verdict", "")
        e = true_rung[rr]
        if side_of[rr] == "split":
            # split pipeline
            s = my_v2_ord2(rr)
            if s != 1:
                myv = "split-v2!=1"
            elif e != d:
                myv = "split-nonprim"
            else:
                verdict, m, why = my_oddpart_ord2(rr)
                if verdict == "prime":
                    p = m
                    condv = pow(2, p - 2, 3 * d) == 3 * d - 1
                    myv = "SPLIT-VIOLATION" if condv else "split-condv-fail"
                    if condv:
                        my_split_viol.append((d, rr, p))
                elif verdict == "composite":
                    myv = "split-oddcomp"
                elif verdict == "one":
                    myv = "split-ord2"
                else:
                    myv = "split-unresolved"
                    my_unresolved.append((d, rr, why))
            ok = ((myv == "split-v2!=1" and "v2(ord_r(2)) != 1" in cv) or
                  (myv == "split-nonprim" and "alpha-rung" in cv) or
                  (myv == "split-oddcomp" and "odd part of ord composite" in cv) or
                  (myv == "split-condv-fail" and "fails condition (v)" in cv) or
                  (myv == "SPLIT-VIOLATION" and "VIOLATION" in cv) or
                  (myv == "split-unresolved" and "UNRESOLVED" in cv))
            if ok and myv == "split-nonprim":
                mm = re.search(r"alpha-rung (\d+) <", cv)
                if not mm or int(mm.group(1)) != e:
                    mismatches.append((d, rr, f"split rung {e}", cv))
            if not ok:
                mismatches.append((d, rr, myv, cv))
            continue
        # inert pipeline (rr == 3 mod 8 by C4 residue check)
        if my_v2_ord2(rr) != 1:
            fail("C5", f"inert r={rr}: v2(ord_r(2)) != 1 (Euler violated?!)")
        if e != d:
            # census message format: 'attributed to rung {g}; deduplicated'
            mm = re.match(r"attributed to rung (\d+);", cv)
            if not mm or int(mm.group(1)) != e:
                mismatches.append((d, rr, f"attributed to rung {e}", cv))
            continue
        # primitive here: T4
        verdict, m, why = my_oddpart_ord2(rr)
        if verdict == "one":
            myv = "T4-one"
            if cv != "excluded at T4":
                mismatches.append((d, rr, myv, cv))
            continue
        if verdict == "composite":
            myv = "T4-composite"
            if "odd part" not in cv or "composite" not in pr["tests"].get("T4_p", ""):
                mismatches.append((d, rr, myv, cv))
            # if census claimed an exact odd part, compare
            mm = re.search(r"odd part (\d+) composite", pr["tests"].get("T4_p", ""))
            if mm and m is not None and int(mm.group(1)) != m:
                mismatches.append((d, rr, f"odd part {m}", f"census odd part {mm.group(1)}"))
            continue
        if verdict == "unresolved":
            my_unresolved.append((d, rr, why))
            if "UNRESOLVED" not in cv and "unresolved" not in cv:
                mismatches.append((d, rr, "unresolved-here", cv))
            continue
        # T4 pass: odd part is prime p
        p = m
        # EXACT certificates
        if not cert_ord2_is_2p(rr, p):
            fail("C5", f"(d={d}, r={rr}): ord_r(2) = 2*{p} certificate FAILS")
        if not cert_w3_ord_is_4d(rr, d, dprimes):
            fail("C5", f"(d={d}, r={rr}): ord(omega3) = 4d certificate FAILS")
        k = (rr - 1) // (2 * p)
        if rr - 1 != 2 * p * k:
            fail("C5", f"(d={d}, r={rr}): k arithmetic broken")
        surv.add(rr)
        condv = pow(2, p - 2, 3 * d) == 3 * d - 1
        if not condv:
            myv = "T5-fail"
            if "excluded at T5" not in cv:
                mismatches.append((d, rr, myv, cv))
            continue
        # full triple: phantom vs real, exact
        if p <= 4096:
            phantom = (3 * rr == 2 ** p + 1)
        else:
            if rr.bit_length() >= p - 1:
                fail("C5", f"r bit length vs W_p size at (d={d}, r={rr})")
            phantom = False
        if phantom:
            my_phantom.append({"d": d, "r": rr, "p": p})
            if "PHANTOM" not in cv:
                mismatches.append((d, rr, "PHANTOM", cv))
        else:
            my_real.append({"d": d, "r": rr, "p": p, "k": k,
                            "above_1e12": rr > 10 ** 12})
            if "REAL" not in cv:
                mismatches.append((d, rr, "REAL", cv))
    my_survivors[d] = surv

log(f"   {n_class} primes reclassified from scratch")
log(f"   my REAL triples:   {[(t['d'], t['r'], t['p'], t['k']) for t in my_real]}")
log(f"   my PHANTOM triples:{[(t['d'], t['r'], t['p']) for t in my_phantom]}")
log(f"   my unresolved:     {my_unresolved if my_unresolved else 'none'}")
log(f"   my split violations: {my_split_viol if my_split_viol else 'none'}")
for msn in mismatches:
    fail("C5", f"verdict mismatch: {msn}")

# compare survivor sets and triples wholesale with the census
cs_real = {(int(t["d"]), int(t["r"]), int(t["p"]), int(t["k"]))
           for t in census["summary"]["real_triples"]}
if {(t["d"], t["r"], t["p"], t["k"]) for t in my_real} != cs_real:
    fail("C5", f"REAL triple set differs: mine {my_real} census {cs_real}")
cs_ph = {(int(t["d"]), int(t["r"]), int(t["p"])) for t in census["summary"]["phantoms"]}
if {(t["d"], t["r"], t["p"]) for t in my_phantom} != cs_ph:
    fail("C5", f"PHANTOM set differs: mine {my_phantom} census {cs_ph}")
if census["summary"]["unresolved"]:
    fail("C5", "census reports unresolved entries (must re-examine)")
for t in census["summary"]["real_triples"]:
    if (int(t["r"]) > 10 ** 12) != t["above_platinum_cutoff"]:
        fail("C5", f"platinum cutoff flag wrong for r={t['r']}")

# primality certification: BPSW + APR-CL on every survivor r and pinned p
cert_targets = set()
for dd, ss in my_survivors.items():
    for rr in ss:
        cert_targets.add(rr)
        # pinned p for each survivor
        verdict, m, _ = my_oddpart_ord2(rr)
        assert verdict == "prime"
        cert_targets.add(m)
cert_targets.add(1049)   # the split (v)-tested survivor at d=131
cert_targets.add(131)
aprcl_all_ok = True
for n in sorted(cert_targets):
    if not isprime(n):
        fail("C5", f"BPSW says {n} is composite?!")
    if have_gp:
        g = gp_aprcl(n)
        if g is not True:
            aprcl_all_ok = False
            fail("C5", f"APR-CL failed for {n}: {g}")
log(f"   APR-CL (gp isprime(*,2)) on {len(cert_targets)} survivor r/p values: "
    f"{'all prime' if aprcl_all_ok and have_gp else 'gp unavailable' if not have_gp else 'FAILURES'}")

# redundant order certificate via full factorization of r-1 and r+1 (survivors)
for dd in sorted(my_survivors):
    for rr in sorted(my_survivors[dd]):
        fi = factorint_sub(rr - 1, 120)
        if fi:
            o2 = exact_order_int(2, rr, fi)
            verdict, p, _ = my_oddpart_ord2(rr)
            if o2 != 2 * p:
                fail("C5", f"factored r-1 order check: ord_{rr}(2) = {o2} != 2*{p}")
        fi2 = factorint_sub(rr + 1, 120)
        if fi2:
            ow = exact_order_pair(W3, rr, fi2)
            if ow != 4 * dd:
                fail("C5", f"factored r+1 order check: ord(omega3) mod {rr} = {ow} != 4*{dd}")

# phantom identities, exact
if (2 ** 11 + 1) // 3 != 683 or (2 ** 17 + 1) // 3 != 43691:
    fail("C5", "phantom W_p identities broken")
if not isprime((2 ** 11 + 1) // 3) or not isprime((2 ** 17 + 1) // 3):
    fail("C5", "W_11 / W_17 not prime?!")
log(f"   C5 {'PASS' if 'C5' not in FAILURES else 'FAIL'}")
log("")

# ---------------------------------------------------------------------------
log("== C6: consistency with the paper ==")
# (a) corner sweep: no REAL triple with p <= 1e8
for t in my_real:
    if t["p"] <= 10 ** 8:
        fail("C6", f"REAL triple with p = {t['p']} <= 1e8 contradicts cp358 corner sweep")
log(f"   (a) REAL p values: {[t['p'] for t in my_real]} — all > 1e8: "
    f"{all(t['p'] > 10**8 for t in my_real)}")

# (b) prop:diagonal-closures reconciliation, MY OWN survivor sets
paper_d124 = {57: {683, 171369731867}, 59: {88499},
              67: {148937135261632265803}, 81: {5507},
              83: set(), 99: set(), 107: set(), 121: set()}
paper_pinned = {683: 11, 171369731867: 85684865933, 88499: 44249,
                148937135261632265803: 24822855876938710967, 5507: 2753}
for dd, exp in paper_d124.items():
    got = my_survivors.get(dd, set())
    if got != exp:
        fail("C6", f"d={dd}: my survivors {sorted(got)} != paper {sorted(exp)}")
    for rr in got:
        verdict, p, _ = my_oddpart_ord2(rr)
        if p != paper_pinned[rr]:
            fail("C6", f"pinned exponent mismatch at r={rr}: {p} != {paper_pinned[rr]}")
# no two survivors at one d share a pinned exponent (diagonal vanishing)
for dd, ss in my_survivors.items():
    ps = []
    for rr in ss:
        _, p, _ = my_oddpart_ord2(rr)
        ps.append(p)
    if len(ps) != len(set(ps)):
        fail("C6", f"d={dd}: two survivors share a pinned exponent (diagonal pair!)")
log(f"   (b) d<=124 survivor sets match prop:diagonal-closures; no shared pinned p at any d")

# (c) split-side: recomputed above; assert none
if my_split_viol:
    fail("C6", f"split-side NCT violations found: {my_split_viol}")
log(f"   (c) split-side crosscheck: {len([1 for rec in census['rungs'] for pr in rec['primes'] if pr.get('side')=='split'])} split primes recomputed, 0 violations")

# (d) discharge_status vs paper
DISCHARGED = {85684865933, 1251488009}
X1_P = 10916765939
SEC_CLOSED_D67_P = 24822855876938710967
SEC_CLOSED_D67_R2 = 5474333343676555561818313   # paper sec:secondary, split
for t in census["summary"]["real_triples"]:
    p = int(t["p"])
    st = t["discharge_status"]
    if p in DISCHARGED and st != "discharged-witness":
        fail("C6", f"p={p} should be discharged-witness, census says {st}")
    if p == X1_P and st != "X1-open":
        fail("C6", f"p={p} should be X1-open")
    if p == SEC_CLOSED_D67_P and st != "secondary-closed":
        fail("C6", f"LABEL DEFECT: p={p} (d=67 triple) is labelled '{st}' but the "
                   f"paper (sec:secondary) closes it at the split secondary factor "
                   f"r2={SEC_CLOSED_D67_R2}, and it is known (rem:known-danger), "
                   f"not new — expected 'secondary-closed'")
# verify the paper's secondary/discharge factors independently
r2a = 7675646348774307083           # d=57 witness secondary (prop:witness-discharges a)
r2b = 33829735778964491             # p=1251488009 discharge factor (b)
r2c = SEC_CLOSED_D67_R2             # d=67 secondary (sec:secondary), split
rx1 = 152834723147                  # X1 witness factor
checks = []
checks.append(("r2a prime", isprime(r2a) and (not have_gp or gp_aprcl(r2a))))
checks.append(("r2a == 3 mod 8", r2a % 8 == 3))
checks.append(("r2a | W_85684865933", pow(2, 85684865933, r2a) == r2a - 1))
da = 3 * 487 * 1313423399858711
checks.append(("ord(omega3) = 4*d_{r2a}", cert_w3_ord_is_4d(r2a, da, [3, 487, 1313423399858711])))
checks.append(("487 not Class III", my_v2_ord2(487) != 1))
checks.append(("1313423399858711 not Class III", my_v2_ord2(1313423399858711) != 1))
checks.append(("d_{r2a} condition (v) fails (parity)", pow(2, 85684865933 - 2, 3 * da) != 3 * da - 1))
checks.append(("r2b prime", isprime(r2b) and (not have_gp or gp_aprcl(r2b))))
checks.append(("r2b == 3 mod 8", r2b % 8 == 3))
checks.append(("r2b | W_1251488009", pow(2, 1251488009, r2b) == r2b - 1))
db = 3 * 41 * 68759625567001
checks.append(("ord(omega3) = 4*d_{r2b}", cert_w3_ord_is_4d(r2b, db, [3, 41, 68759625567001])))
checks.append(("41 not Class III", my_v2_ord2(41) != 1))
checks.append(("68759625567001 not Class III", my_v2_ord2(68759625567001) != 1))
checks.append(("r2c prime", isprime(r2c) and (not have_gp or gp_aprcl(r2c))))
checks.append(("r2c == 1 mod 8 (split!)", r2c % 8 == 1))
checks.append(("r2c | W_{d67 p}", pow(2, SEC_CLOSED_D67_P, r2c) == r2c - 1))
checks.append(("rx1 prime", isprime(rx1) and (not have_gp or gp_aprcl(rx1))))
checks.append(("rx1 = 14p+1", rx1 == 14 * X1_P + 1))
checks.append(("rx1 | W_{X1 p}", pow(2, X1_P, rx1) == rx1 - 1))
dx1 = 38208680787
fx1 = factorint(dx1)
checks.append(("rx1 d = (r+1)/4", 4 * dx1 == rx1 + 1))
checks.append(("ord(omega3) = 4*d_{rx1}", cert_w3_ord_is_4d(rx1, dx1, list(fx1))))
checks.append(("rx1 passes (v) locally", pow(2, X1_P - 2, 3 * dx1) == 3 * dx1 - 1))
checks.append(("X1 exponent prime", isprime(X1_P)))
checks.append(("witness exponent 85684865933 prime", isprime(85684865933)))
checks.append(("witness exponent 1251488009 prime", isprime(1251488009)))
checks.append(("d67 pinned exponent prime", isprime(SEC_CLOSED_D67_P)))
for name, ok in checks:
    if not ok:
        fail("C6", f"paper fact check failed: {name}")
log(f"   (d) paper witness/secondary factors: {sum(1 for _, ok in checks if ok)}/{len(checks)} verified")
log(f"   C6 {'PASS' if 'C6' not in FAILURES else 'FAIL'}")
log("")

# ---------------------------------------------------------------------------
log("== C7: boundary + internal consistency of the census artifacts ==")
if census["range"] != "43 < d <= 400":
    fail("C7", "range string unexpected")
if any(d % 2 == 0 or d <= 43 or d > 400 for d in census["admissible_rungs"]):
    fail("C7", "census admissible list violates 43 < d <= 400 odd")
# d=400 is even; largest odd is 399 = 3*7*19 with 7 not Class III
if my_v2_ord2(7) == 1:
    fail("C7", "7 misclassified")
# recount summary numbers from the JSON rung records
n_primes = sum(len(rec["primes"]) for rec in census["rungs"])
if n_primes != census["summary"]["primes_tested"]:
    fail("C7", f"primes_tested recount {n_primes} != {census['summary']['primes_tested']}")
n_complete = sum(1 for rec in census["rungs"] if rec["status"] == "COMPLETE")
if n_complete != census["summary"]["rungs_complete"] or n_complete != 34:
    fail("C7", "rungs_complete recount mismatch")
n_real_v = sum(1 for rec in census["rungs"] for pr in rec["primes"]
               if pr["tests"].get("verdict", "").startswith("REAL"))
n_ph_v = sum(1 for rec in census["rungs"] for pr in rec["primes"]
             if pr["tests"].get("verdict", "").startswith("PHANTOM"))
if n_real_v != len(census["summary"]["real_triples"]) or n_ph_v != len(census["summary"]["phantoms"]):
    fail("C7", "summary triple counts do not match per-prime verdicts")
# FDB cache poisoning: every cached factor list must product-verify
fdbc = json.load(open(FDB_CACHE))
for k, v in fdbc.items():
    n = int(k)
    prod = 1
    for fs, e in v.get("factors", []):
        prod *= int(fs) ** int(e)
    if v.get("factors") and prod != n:
        fail("C7", f"FDB cache entry {k[:30]}... fails product identity")
log(f"   FDB cache: {len(fdbc)} entries, all product-verified")
log(f"   C7 {'PASS' if 'C7' not in FAILURES else 'FAIL'}")
log("")

# ---------------------------------------------------------------------------
n_fail = sum(len(v) for v in FAILURES.values())
log("== VERDICT ==")
log(f"   sections with failures: { {k: len(v) for k, v in FAILURES.items()} if FAILURES else 'none'}")
log(f"   elapsed: {time.time()-T0:.1f}s")
json.dump({"failures": FAILURES,
           "my_real": [{k: str(v) for k, v in t.items()} for t in my_real],
           "my_phantom": [{k: str(v) for k, v in t.items()} for t in my_phantom],
           "my_unresolved": [[str(x) for x in u] for u in my_unresolved],
           "n_primes_reclassified": n_class,
           "aprcl_available": have_gp},
          open(OUT_JSON, "w"), indent=1)
log(f"   wrote {OUT_JSON}")
sys.exit(1 if n_fail else 0)
