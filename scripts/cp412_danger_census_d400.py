#!/usr/bin/env python3
"""cp412: complete danger-triple census for 43 < d <= 400 (Paper 02, def:danger-triple).

For every admissible odd d in (43, 400] (all prime factors Class III, i.e. in
A = {q odd prime : v_2(ord_q(2)) = 1}), assemble the COMPLETE verified
factorization of V_{2d} (Pell-Lucas), then run the danger tests on every prime
factor r > 3:

  (T1) r == 3 (mod 8)
  (T2) v_2(ord_r(2)) = 1               [factorization-free square chain]
  (T3) ord_r(omega_3) = 4d in F_{r^2}  [primitivity at rung d; else attribute
                                        r to its true rung e = d/f]
  (T4) odd part of ord_r(2) is a PRIME p; record p, k = (r-1)/(2p)
  (T5) d | W_{p-2}  <=>  2^(p-2) == -1 (mod 3d)

Full T1-T5 passers: PHANTOM if r == W_p = (2^p+1)/3 (W_p prime), else REAL.

Cross-checks:
  (a) reproduce prop:diagonal-closures survivor sets at d <= 124 exactly,
      decide condition (v) for each survivor;
  (b) corner-sweep consistency: no REAL triple with p <= 10^8 (cp358 sweep);
  (c) split-side NCT: at each complete rung no r == 1 (mod 8) with
      {v_2(ord_r(2)) = 1, odd part prime p, ord_r(alpha) = 8d, d | W_{p-2}}.

Non-admissible odd d in (43,400] cannot carry a danger triple BY DEFINITION
(def:danger-triple requires every prime factor of d in A) -- equivalently by
the parity obstruction (thm:parity-obstruction): d | W_{p-2} forces every
prime factor of d into A.  Their V_{2d} are not factored here.

Factorization sources (priority order):
  1. papers/02_wagstaff_chebyshev_reduction/zenodo_v3.6/data/nct_certificates.json
     (28 pinned Psi_{4d} distinct-prime lists, (100,400], APR-CL certified)
  2. scripts/output/cp398_factordb_raw.json (complete V_{2d} FF, verified here)
  3. scripts/output/cp412_factordb_cache.json (FactorDB pulls, >=2s apart)
  4. local factoring: coprime halves V_{2d} = 2*(2U_d-1)*(2U_d+1), sympy
     factorint, gmp-ecm fallback (/tmp/ecmenv/bin/ecm), ~10 min cap/cofactor.

Every assembled factorization is verified by the exact product identity
(with multiplicities) against the locally recomputed V_{2d}, and every listed
prime is BPSW-checked (sympy.isprime).  REAL-triple r and p and phantom r are
APR-CL certified via PARI/GP isprime(*,2).

Outputs: scripts/output/cp412_danger_census.json, .log
"""

import json
import os
import subprocess
import sys
import time

import sympy
from sympy import isprime, factorint, primerange

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_JSON = os.path.join(ROOT, "scripts/output/cp412_danger_census.json")
OUT_LOG = os.path.join(ROOT, "scripts/output/cp412_danger_census.log")
FDB_CACHE_PATH = os.path.join(ROOT, "scripts/output/cp412_factordb_cache.json")
NCT_JSON = os.path.join(
    ROOT, "papers/02_wagstaff_chebyshev_reduction/zenodo_v3.6/data/nct_certificates.json")
CP398_RAW = os.path.join(ROOT, "scripts/output/cp398_factordb_raw.json")
ECM_BIN = "/tmp/ecmenv/bin/ecm"

D_MIN, D_MAX = 43, 400          # census range: 43 < d <= 400
PLATINUM_CUTOFF = 10**12
CORNER_P_LIMIT = 10**8
DISCHARGED_WITNESS_P = {85684865933, 1251488009}   # prop:witness-discharges
X1_OPEN_P = 10916765939                            # conj:x1
# sec:secondary: the known d=67 real triple's exponent is closed by a SPLIT
# secondary factor r2 of W_p (the all-inert hypothesis fails at p; "R3
# verified directly" in the paper).  r2 is re-verified at runtime (prime,
# == 1 mod 8, 2^p == -1 mod r2) before the label is applied.
SECONDARY_CLOSED_P = {24822855876938710967: 5474333343676555561818313}

# prop:diagonal-closures expected survivor sets of tests (i)-(iv), d <= 124
EXPECTED_D124_SURVIVORS = {
    57: {683, 171369731867},
    59: {88499},
    67: {148937135261632265803},
    81: {5507},
    83: set(), 99: set(), 107: set(), 121: set(),
}
EXPECTED_D124_PINNED_P = {
    683: 11, 171369731867: 85684865933, 88499: 44249,
    148937135261632265803: 24822855876938710967, 5507: 2753,
}


class Log:
    def __init__(self, path):
        self.fh = open(path, "w")

    def __call__(self, msg=""):
        print(msg, flush=True)
        self.fh.write(msg + "\n")
        self.fh.flush()


log = Log(OUT_LOG)


# ----------------------------------------------------------------------------
# Pell arithmetic (exact bigint)
# ----------------------------------------------------------------------------

def pell_UV(n):
    """(U_n, V_n): U_0=0,U_1=1; V_0=2,V_1=2; X_{n+1} = 2X_n + X_{n-1}."""
    u0, u1 = 0, 1
    v0, v1 = 2, 2
    for _ in range(n):
        u0, u1 = u1, 2 * u1 + u0
        v0, v1 = v1, 2 * v1 + v0
    return u0, v0


def V(n):
    return pell_UV(n)[1]


def U(n):
    return pell_UV(n)[0]


def zs2_mul(a, b):
    """multiply (x + y*sqrt2) pairs exactly in Z[sqrt2]."""
    return (a[0] * b[0] + 2 * a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def zs2_pow(a, e):
    r = (1, 0)
    while e:
        if e & 1:
            r = zs2_mul(r, a)
        a = zs2_mul(a, a)
        e >>= 1
    return r


def v2(n):
    return (n & -n).bit_length() - 1


# ----------------------------------------------------------------------------
# F_{r^2} = F_r[sqrt2] arithmetic (2 a non-residue mod r), and F_r analog
# ----------------------------------------------------------------------------

def fr2_pow(base, e, r):
    """(x + y*sqrt2)^e mod r as pairs; requires (2|r) = -1 for field."""
    res = (1, 0)
    b = (base[0] % r, base[1] % r)
    while e:
        if e & 1:
            res = ((res[0] * b[0] + 2 * res[1] * b[1]) % r,
                   (res[0] * b[1] + res[1] * b[0]) % r)
        b = ((b[0] * b[0] + 2 * b[1] * b[1]) % r, (2 * b[0] * b[1]) % r)
        e >>= 1
    return res


OMEGA3 = (3, 2)   # 3 + 2*sqrt2 = alpha^2, alpha = 1 + sqrt2


def omega3_2m_is_minus1(r, m):
    """omega_3^{2m} == -1 in F_{r^2}?"""
    return fr2_pow(OMEGA3, 2 * m, r) == (r - 1, 0)


def inert_true_rung(r, d, d_prime_factors):
    """Given inert r | V_{2d} (so omega_3^{2d} = -1), return g0 with
    ord_r(omega_3) = 4*g0, g0 | d.  Asserts the -1 sanity condition."""
    assert pow(2, (r - 1) // 2, r) == r - 1, f"2 is a QR mod {r}: not inert"
    assert omega3_2m_is_minus1(r, d), f"omega3^(2d) != -1 mod {r} at d={d}"
    g = d
    for l in d_prime_factors:
        while g % l == 0 and omega3_2m_is_minus1(r, g // l):
            g //= l
    # invariant: omega3^{2m} = -1 <=> g0 | m (m odd | d); descent is exact
    assert (r + 1) % (4 * g) == 0, f"4*{g} does not divide r+1 for r={r}"
    return g


def split_alpha_true_rung(r, d, d_prime_factors):
    """r == 1 mod 8, r | V_{2d}: alpha in F_r, alpha^{4d} = -1.  Return g0
    with ord_r(alpha) = 8*g0, g0 | d."""
    s = sympy.sqrt_mod(2, r)
    a = (1 + s) % r
    assert pow(a, 4 * d, r) == r - 1, f"alpha^(4d) != -1 mod split r={r}"
    g = d
    for l in d_prime_factors:
        while g % l == 0 and pow(a, 4 * (g // l), r) == r - 1:
            g //= l
    assert (r - 1) % (8 * g) == 0
    return g


def v2_ord2(r):
    """v_2(ord_r(2)) factorization-free: y = 2^o mod r (o = odd part of r-1);
    v2 = number of squarings of y to reach 1 (0 if y == 1)."""
    t = v2(r - 1)
    o = (r - 1) >> t
    y = pow(2, o, r)
    if y == 1:
        return 0
    for s in range(1, t + 1):
        y = y * y % r
        if y == 1:
            return s
    raise AssertionError(f"2^(r-1) != 1 mod {r}: r not prime?")


# ----------------------------------------------------------------------------
# FactorDB (cached, throttled)
# ----------------------------------------------------------------------------

_fdb_cache = {}
_fdb_last_t = [0.0]


def fdb_load():
    global _fdb_cache
    if os.path.exists(FDB_CACHE_PATH):
        _fdb_cache = json.load(open(FDB_CACHE_PATH))


def fdb_save():
    json.dump(_fdb_cache, open(FDB_CACHE_PATH, "w"), indent=1)


def fdb_query(n):
    """Raw FactorDB API answer for n, cached.  >= 2s between live requests."""
    key = str(n)
    if key in _fdb_cache:
        return _fdb_cache[key]
    import urllib.request
    import urllib.parse
    wait = 2.2 - (time.time() - _fdb_last_t[0])
    if wait > 0:
        time.sleep(wait)
    url = "http://factordb.com/api?query=" + urllib.parse.quote(key)
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=30) as resp:
                data = json.load(resp)
            _fdb_last_t[0] = time.time()
            _fdb_cache[key] = {"fetched": time.strftime("%Y-%m-%d %H:%M:%S"),
                               "status": data.get("status"),
                               "factors": data.get("factors", [])}
            fdb_save()
            return _fdb_cache[key]
        except Exception as e:
            log(f"    [fdb] attempt {attempt+1} failed for {str(n)[:40]}...: {e}")
            _fdb_last_t[0] = time.time()
            time.sleep(5)
    return None


# ----------------------------------------------------------------------------
# local factoring helpers
# ----------------------------------------------------------------------------

def factorint_subprocess(n, timeout_s):
    """sympy.factorint in a subprocess with a hard timeout; None on timeout."""
    code = ("import sympy,json,sys\n"
            f"f=sympy.factorint({n})\n"
            "print(json.dumps({str(k):v for k,v in f.items()}))\n")
    try:
        res = subprocess.run(["nice", "-n", "10", sys.executable, "-c", code],
                             capture_output=True, text=True, timeout=timeout_s)
        if res.returncode == 0 and res.stdout.strip():
            return {int(k): v for k, v in json.loads(res.stdout).items()}
    except subprocess.TimeoutExpired:
        pass
    return None


def ecm_find_factor(n, budget_s=600):
    """gmp-ecm schedule under a total time budget; returns one nontrivial
    factor or None."""
    if not os.path.exists(ECM_BIN):
        return None
    t0 = time.time()
    schedule = [(2000, 50), (11000, 60), (50000, 60), (250000, 40), (10**6, 30)]
    for b1, curves in schedule:
        rem = budget_s - (time.time() - t0)
        if rem <= 5:
            return None
        try:
            res = subprocess.run(
                ["nice", "-n", "10", ECM_BIN, "-q", "-c", str(curves), str(b1)],
                input=str(n), capture_output=True, text=True,
                timeout=min(rem, budget_s))
            for tok in res.stdout.split():
                if tok.isdigit():
                    f = int(tok)
                    if 1 < f < n and n % f == 0:
                        return f
        except subprocess.TimeoutExpired:
            return None
    return None


def full_local_factor(n, budget_s=600, label=""):
    """Complete factorization dict or None.  Trial -> factorint(sub) -> ECM."""
    if n == 1:
        return {}
    out = {}
    rem = n
    for q in primerange(2, 100000):
        if rem == 1 or q * q > rem:
            break
        while rem % q == 0:
            out[q] = out.get(q, 0) + 1
            rem //= q
    if rem > 1 and (rem < 10**10 or isprime(rem)):
        # q*q > rem exit means rem is prime; also catches small residues
        out[rem] = out.get(rem, 0) + 1
        rem = 1
    if rem == 1:
        return out
    if isprime(rem):
        out[rem] = out.get(rem, 0) + 1
        return out
    t0 = time.time()
    stack = [rem]
    while stack:
        m = stack.pop()
        if isprime(m):
            out[m] = out.get(m, 0) + 1
            continue
        left = budget_s - (time.time() - t0)
        if left <= 5:
            log(f"    [local] budget exhausted on {label} cofactor ({len(str(m))}d)")
            return None
        f = None
        if len(str(m)) <= 46:
            fi = factorint_subprocess(m, min(left, 120))
            if fi:
                for p_, e_ in fi.items():
                    for _ in range(e_):
                        stack.append(p_)
                continue
        f = ecm_find_factor(m, budget_s=left)
        if f is None:
            log(f"    [local] cannot factor {label} cofactor ({len(str(m))}d)")
            return None
        stack.append(f)
        stack.append(m // f)
    return out


# ----------------------------------------------------------------------------
# ord_r(2) odd-part analysis with partial factorization (exact where decided)
# ----------------------------------------------------------------------------

def _harder_factor_pass(n, budget_s, use_fdb):
    """Best-effort factor pass on composite n: FactorDB, factorint (small),
    gmp-ecm.  Returns the set of prime factors found (exponents irrelevant
    for the order descent)."""
    newp = set()
    t0 = time.time()
    stack = [n]
    seen = set()
    while stack:
        x = stack.pop()
        if x == 1 or x in seen:
            continue
        seen.add(x)
        if isprime(x):
            newp.add(x)
            continue
        if use_fdb and len(str(x)) >= 18:
            ans = fdb_query(x)
            if ans and ans.get("factors"):
                y = x
                progressed = False
                for fs, _e in ans["factors"]:
                    f = int(fs)
                    if f == x:
                        continue
                    while y % f == 0:
                        stack.append(f)
                        y //= f
                        progressed = True
                if progressed:
                    if y > 1:
                        stack.append(y)
                    continue
        left = budget_s - (time.time() - t0)
        if left <= 5:
            continue
        if len(str(x)) <= 46:
            fi = factorint_subprocess(x, min(left, 120))
            if fi:
                newp.update(fi.keys())
                continue
        f = ecm_find_factor(x, budget_s=left)
        if f is not None:
            stack.append(f)
            stack.append(x // f)
    return newp


def ord2_oddpart(r, ecm_budget_s=600, use_fdb=True):
    """For prime r with v2(ord_r(2)) = 1 (T2 already passed): decide whether
    the odd part p of ord_r(2) = 2p is prime.

    Returns (verdict, payload):
      ('prime', p)            -- odd part is the prime p (exact)
      ('composite', detail)   -- odd part provably composite (exact)
      ('one', None)           -- odd part 1 (r | 3; excluded upstream)
      ('unresolved', detail)  -- cannot decide (factoring failed)

    Method: o = odd part of r-1; the odd part m' of ord divides o.  Descend
    over the known prime factors of o; keep the unfactored part of o as one
    block m1 (coprime to the knowns by construction): if 2^{2*(m/m1)} == 1
    the whole block strips (m' | m/m1, exact); if the block does not strip,
    some prime of m1 divides m', and if additionally a known prime remains
    in m, m' has >= 2 distinct prime factors -> composite (exact, no
    factorization of m1 needed).  Only the ambiguous corner (m == m1
    composite) escalates to real factoring (FactorDB / factorint / ECM)."""
    t = v2(r - 1)
    o = (r - 1) >> t
    assert pow(2, o, r) == r - 1, "T2 precondition violated"
    known = set()
    rem = o
    for q in primerange(2, 10**6):
        if rem == 1 or q * q > rem:
            break
        while rem % q == 0:
            known.add(q)
            rem //= q
    if rem > 1 and isprime(rem):
        known.add(rem)
        rem = 1
    m1 = rem  # composite block, coprime to known primes, m1 | o

    def rebuild_block():
        r2 = o
        for q in known:
            while r2 % q == 0:
                r2 //= q
        return r2

    def descend(block):
        m = o
        for l in sorted(known):
            while m % l == 0 and pow(2, 2 * (m // l), r) == 1:
                m //= l
        if block > 1:
            # block | m always (block coprime to knowns, block | o)
            assert m % block == 0
        if block > 1 and pow(2, 2 * (m // block), r) == 1:
            m //= block
            block = 1
            for l in sorted(known):
                while m % l == 0 and pow(2, 2 * (m // l), r) == 1:
                    m //= l
        return m, block

    escalated = False
    while True:
        m, m1r = descend(m1)
        if m > 1:
            assert pow(2, m, r) == r - 1, "descent invariant (2^m = -1)"
        if m1r == 1:
            if m == 1:
                return ("one", None)
            if isprime(m):
                return ("prime", m)
            return ("composite",
                    f"odd part {m} composite (fully factored descent)")
        s_min = m // m1r
        if s_min > 1:
            return ("composite",
                    "odd part has a known prime factor and a factor of the "
                    f"{len(str(m1r))}-digit unfactored block -> composite")
        # ambiguous corner: m == m1r, composite, unfactored
        if escalated:
            return ("unresolved",
                    f"odd part divides unfactored {len(str(m1r))}-digit "
                    "composite block")
        escalated = True
        newp = _harder_factor_pass(m1r, ecm_budget_s, use_fdb)
        if not newp:
            return ("unresolved",
                    f"odd part divides unfactored {len(str(m1r))}-digit "
                    "composite block (FDB/factorint/ECM found nothing)")
        known |= newp
        m1 = rebuild_block()


# ----------------------------------------------------------------------------
# APR-CL via PARI/GP
# ----------------------------------------------------------------------------

def aprcl_certify(n, timeout_s=1800):
    try:
        res = subprocess.run(["nice", "-n", "10", "gp", "-q"],
                             input=f"print(isprime({n},2))\n",
                             capture_output=True, text=True, timeout=timeout_s)
        return res.stdout.strip() == "1"
    except Exception:
        return None


# ----------------------------------------------------------------------------
# admissibility
# ----------------------------------------------------------------------------

def is_class_iii(q):
    """q odd prime in A  <=>  v_2(ord_q(2)) = 1  <=>  2^(odd part of q-1) == -1."""
    return v2_ord2(q) == 1


def admissible_rungs():
    out = []
    for d in range(D_MIN + 2 if D_MIN % 2 == 1 else D_MIN + 1, D_MAX + 1, 2):
        if d <= D_MIN:
            continue
        fac = factorint(d)
        if all(is_class_iii(q) for q in fac):
            out.append(d)
    return out


# ----------------------------------------------------------------------------
# self-tests of every mathematical fact used
# ----------------------------------------------------------------------------

def selftests():
    log("== self-tests ==")
    # Pell recurrences vs alpha-power identities
    for n in [1, 2, 3, 5, 10, 17, 40]:
        u, v = pell_UV(n)
        a = zs2_pow((1, 1), n)              # alpha^n = V_n/2 + U_n sqrt2
        assert (2 * a[0], a[1]) == (v, u), f"alpha^{n} mismatch"
        # V_n = alpha^n + (1-sqrt2)^n : conjugate sum = 2*a[0]
        assert v == 2 * a[0]
    # omega_3 = alpha^2
    assert zs2_pow((1, 1), 2) == (3, 2)
    # identities for odd d
    for d in [1, 3, 5, 9, 11, 57]:
        u, v = pell_UV(d)
        v2d = V(2 * d)
        assert v2d == v * v + 2, f"V_2d != V_d^2+2 at d={d}"
        assert v2d == 2 * (2 * u - 1) * (2 * u + 1), f"half split fails d={d}"
        assert v * v - 8 * u * u == -4, f"V^2-8U^2 != -4 at odd d={d}"
    # W divisibility test: D | W_n <=> 2^n == -1 (mod 3D), incl. prime powers
    for n in [3, 5, 7, 9, 11, 15, 21, 27]:
        Wn = (2**n + 1) // 3
        for D in [3, 5, 9, 11, 19, 27, 33, 43, 57, 81, 121, 171, 243, 331]:
            direct = (Wn % D == 0)
            viatest = (pow(2, n, 3 * D) == 3 * D - 1)
            assert direct == viatest, f"W-div test fails n={n} D={D}"
    # class III examples: 3,11,19,43,59,67,83,107 in A; 5,7,17,23,29,31,41 not
    for q in [3, 11, 19, 43, 59, 67, 83, 107, 131, 163, 179, 211, 227, 251,
              281, 283, 307, 331, 347, 379]:
        assert is_class_iii(q), f"{q} should be Class III"
    for q in [5, 7, 17, 23, 29, 31, 37, 41, 53, 61, 73, 89, 97, 101, 103, 109,
              113, 127, 137, 151, 157, 173, 181, 193, 197, 199]:
        assert not is_class_iii(q), f"{q} should NOT be Class III"
    # inert primitivity spot checks: 11 at rung 3 (11 | V_6 = 198)
    assert V(6) % 11 == 0
    assert inert_true_rung(11, 3, [3]) == 3
    # 683 | V_114 (from the paper d=57 factorization) -- rung determined below
    assert V(114) % 683 == 0
    # 43691 = W_17 | V_662 (from nct d=331 list)
    assert V(662) % 43691 == 0
    # v2_ord2 examples
    assert v2_ord2(11) == 1 and v2_ord2(19) == 1
    assert v2_ord2(17) == 3 and v2_ord2(41) == 2 and v2_ord2(7) == 0
    # T2 for r == 3 mod 8 is automatic (Euler: (2|r) = -1): sample check
    for r in [11, 19, 43, 59, 67, 83, 107, 131, 683, 5507, 88499]:
        if r % 8 == 3:
            assert v2_ord2(r) == 1
    # ord2_oddpart on knowns: 43691 -> p = 17
    assert ord2_oddpart(43691, use_fdb=False) == ("prime", 17)
    v_, p_ = ord2_oddpart(683, use_fdb=False)
    assert (v_, p_) == ("prime", 11)
    assert ord2_oddpart(88499, use_fdb=False) == ("prime", 44249)
    assert ord2_oddpart(5507, use_fdb=False) == ("prime", 2753)
    # split-side arithmetic: r = 89 == 1 mod 8 divides V_22; true alpha-rung?
    assert V(22) % 89 == 0
    g89 = split_alpha_true_rung(89, 11, [11])
    assert g89 in (1, 11) and (89 - 1) % (8 * g89) == 0
    log("   all self-tests PASS")


# ----------------------------------------------------------------------------
# factorization assembly
# ----------------------------------------------------------------------------

_leaf_cache = {}     # m -> (factor dict of V_{2m}, source string)


def leaf_V_factorization(m, cp398, allow_fdb=True):
    """Complete verified factorization of V_{2m} for small/leaf m (m <= 100).
    Returns (dict prime->exp, source) or (None, reason)."""
    if m in _leaf_cache:
        return _leaf_cache[m]
    v = V(2 * m)
    src = None
    fac = None
    key = str(2 * m)
    if key in cp398:
        entry = cp398[key]
        fac = {}
        for fs, e in entry["factors"]:
            f = int(fs)
            fac[f] = fac.get(f, 0) + int(e)
        src = "cp398_factordb_raw.json"
    if fac is None and v < 10**36:
        fac = full_local_factor(v, budget_s=300, label=f"V_{2*m}")
        src = "local (sympy/trial)"
    if fac is None and allow_fdb:
        ans = fdb_query(v)
        if ans and ans.get("status") in ("FF", "P", "PRP", "PP"):
            fac = {}
            for fs, e in ans["factors"]:
                f = int(fs)
                fac[f] = fac.get(f, 0) + int(e)
            src = "FactorDB (cp412 cache)"
        elif ans and ans.get("factors"):
            # partial: take prime factors, factor the rest locally
            fac = {}
            remv = v
            for fs, e in ans["factors"]:
                f = int(fs)
                if isprime(f):
                    while remv % f == 0:
                        fac[f] = fac.get(f, 0) + 1
                        remv //= f
            sub = full_local_factor(remv, budget_s=600, label=f"V_{2*m} cofactor")
            if sub is None:
                _leaf_cache[m] = (None, f"V_{2*m} partial on FactorDB, local factoring failed")
                return _leaf_cache[m]
            for p_, e_ in sub.items():
                fac[p_] = fac.get(p_, 0) + e_
            src = "FactorDB partial + local"
    if fac is None:
        fac = full_local_factor(v, budget_s=600, label=f"V_{2*m}")
        src = "local (sympy/ecm)"
    if fac is None:
        _leaf_cache[m] = (None, f"V_{2*m} could not be fully factored")
        return _leaf_cache[m]
    # verify exact product and BPSW primality
    prod = 1
    for p_, e_ in fac.items():
        assert isprime(p_), f"leaf V_{2*m}: listed factor {p_} not prime (BPSW)"
        prod *= p_ ** e_
    assert prod == v, f"leaf V_{2*m}: product identity FAILS"
    _leaf_cache[m] = (fac, src)
    return _leaf_cache[m]


def assemble_rung(d, nct_by_d, cp398):
    """Complete verified factorization of V_{2d}.
    Returns (dict prime->exp or None, source string, note)."""
    v = V(2 * d)
    assert v % 2 == 0 and v % 4 != 0, "v_2(V_2d) must be exactly 1"
    odd_divs = [m for m in sorted(sympy.divisors(d)) if m % 2 == 1]
    pool = {2, 3}
    src_parts = []
    for m in odd_divs:
        if m == 1:
            continue
        if m in nct_by_d:
            pool.update(int(f["r"]) for f in nct_by_d[m]["factors"])
            src_parts.append(f"nct_certificates.json Psi_{4*m}")
            # cp398 cross-check pool enrichment (harmless, division-verified)
            if str(2 * m) in cp398:
                pool.update(int(fs) for fs, e in cp398[str(2 * m)]["factors"])
                src_parts.append(f"cp398 V_{2*m} (cross)")
        else:
            fac, src = leaf_V_factorization(m, cp398)
            if fac is None:
                return None, "", f"leaf V_{2*m} incomplete: {src}"
            pool.update(fac.keys())
            src_parts.append(f"{src} V_{2*m}" if m < d else f"{src} V_{2*d}")
    fac = {}
    rem = v
    for r in sorted(pool):
        while rem % r == 0:
            fac[r] = fac.get(r, 0) + 1
            rem //= r
    note = ""
    if rem > 1:
        if isprime(rem):
            fac[rem] = fac.get(rem, 0) + 1
            note = f"residual {len(str(rem))}-digit prime added (BPSW)"
        else:
            sub = full_local_factor(rem, budget_s=600, label=f"V_{2*d} residual")
            if sub is None:
                return None, "; ".join(src_parts), \
                    f"INCOMPLETE: composite residual of {len(str(rem))} digits"
            for p_, e_ in sub.items():
                fac[p_] = fac.get(p_, 0) + e_
            note = "residual factored locally"
    # verification: exact product with multiplicity + BPSW on every prime
    prod = 1
    for p_, e_ in fac.items():
        assert isprime(p_), f"d={d}: listed factor {p_} fails BPSW"
        prod *= p_ ** e_
    assert prod == v, f"d={d}: product identity FAILS"
    # structural sanity: every odd prime factor == 1 or 3 (mod 8)
    for p_ in fac:
        if p_ > 2:
            assert p_ % 8 in (1, 3), f"d={d}: prime {p_} == {p_%8} mod 8 in V_2d!"
    return fac, "; ".join(sorted(set(src_parts))), note


# ----------------------------------------------------------------------------
# main census
# ----------------------------------------------------------------------------

def main():
    t_start = time.time()
    log(f"cp412 danger-triple census, 43 < d <= {D_MAX}   "
        f"({time.strftime('%Y-%m-%d %H:%M:%S')})")
    log("")
    selftests()

    fdb_load()
    nct = json.load(open(NCT_JSON))
    nct_by_d = {c["d"]: c for c in nct["certificates"]}
    cp398 = json.load(open(CP398_RAW))
    # verify cp398 entries once (product identity)
    for k, entry in cp398.items():
        n = int(k)
        prod = 1
        for fs, e in entry["factors"]:
            prod *= int(fs) ** int(e)
        assert prod == V(n), f"cp398 V_{n} product identity fails"
    log(f"   sources: nct_certificates.json ({len(nct_by_d)} rungs: "
        f"{sorted(nct_by_d)}),")
    log(f"            cp398_factordb_raw.json (V_2m for 2m in {sorted(map(int, cp398))}, all verified)")
    log("")

    adm = admissible_rungs()
    log(f"== admissible rungs 43 < d <= {D_MAX} (every prime factor Class III) ==")
    log(f"   {len(adm)} rungs: {adm}")
    log("   Non-admissible odd d cannot carry a danger triple: def:danger-triple")
    log("   requires every prime factor of d in A (equivalently, d | W_{p-2}")
    log("   forces this by the parity obstruction thm:parity-obstruction).")
    log("   Their V_{2d} are not factored.")
    d124 = [d for d in adm if d <= 124]
    assert d124 == [57, 59, 67, 81, 83, 99, 107, 121], \
        f"admissible d<=124 mismatch with prop:diagonal-closures: {d124}"
    assert sorted(set(adm) & set(nct_by_d)) == sorted(nct_by_d), \
        "nct certificate rungs must all be admissible"
    assert [d for d in adm if d > 100] == sorted(nct_by_d), \
        "(100,400] admissible list must equal the 28 nct rungs"
    log("")

    rung_records = []
    real_triples = []
    phantoms = []
    unresolved = []
    incomplete = []
    survivors_by_d = {}      # d -> set of T1-T4 passers (tests i-iv)
    split_violations = []
    primes_tested = 0
    dedup_notes = []

    for d in adm:
        v = V(2 * d)
        dfac = sorted(factorint(d))
        log(f"== d = {d}  (V_{2*d}: {len(str(v))} digits) ==")
        fac, source, note = assemble_rung(d, nct_by_d, cp398)
        rec = {"d": d, "V2d_digits": len(str(v)), "source": source,
               "note": note, "status": None, "primes": []}
        if fac is None:
            rec["status"] = "INCOMPLETE"
            incomplete.append(d)
            log(f"   INCOMPLETE: {note}")
            rung_records.append(rec)
            continue
        rec["status"] = "COMPLETE"
        rec["distinct_primes"] = len(fac)
        log(f"   complete factorization verified: {len(fac)} distinct primes, "
            f"product identity OK, all BPSW prime")
        log(f"   source: {source}" + (f"; {note}" if note else ""))

        surv = set()
        for r in sorted(fac):
            if r <= 3:
                continue
            primes_tested += 1
            pr = {"r": str(r), "digits": len(str(r)), "mult": fac[r],
                  "mod8": r % 8, "tests": {}}
            rec["primes"].append(pr)
            T = pr["tests"]
            if r % 8 == 1:
                # split side: cross-check (c)
                pr["side"] = "split"
                s = v2_ord2(r)
                T["split_v2ord2"] = f"v2(ord_r(2))={s}"
                if s != 1:
                    T["verdict"] = "excluded: split with v2(ord_r(2)) != 1 (no prime half-order)"
                    continue
                g = split_alpha_true_rung(r, d, dfac)
                T["split_alpha_rung"] = f"ord_r(alpha)=8*{g}"
                if g != d:
                    T["verdict"] = f"excluded: split, alpha-rung {g} < {d} (non-primitive here)"
                    continue
                verdict, payload = ord2_oddpart(r)
                if verdict == "prime":
                    p = payload
                    T["split_p"] = str(p)
                    condv = pow(2, p - 2, 3 * d) == 3 * d - 1
                    T["split_condv"] = "holds" if condv else "fails"
                    if condv:
                        T["verdict"] = "SPLIT NCT VIOLATION"
                        split_violations.append((d, r, p))
                        log(f"   !!! SPLIT NCT VIOLATION at d={d}, r={r}, p={p}")
                    else:
                        T["verdict"] = "excluded: split survivor fails condition (v)"
                elif verdict == "composite":
                    T["verdict"] = f"excluded: split, odd part of ord composite ({payload})"
                elif verdict == "one":
                    T["verdict"] = "excluded: split, ord_r(2) = 2"
                else:
                    T["verdict"] = f"UNRESOLVED split candidate: {payload}"
                    unresolved.append({"d": d, "r": str(r), "side": "split",
                                       "detail": payload})
                    log(f"   ?? unresolved split candidate r={r}")
                continue

            # inert side: T1-T5
            pr["side"] = "inert"
            T["T1_mod8"] = "pass (r == 3 mod 8)" if r % 8 == 3 else f"fail (r == {r%8} mod 8)"
            if r % 8 != 3:
                T["verdict"] = "excluded at T1"
                continue
            s = v2_ord2(r)
            T["T2_v2ord2"] = f"pass (v2 = 1)" if s == 1 else f"fail (v2 = {s})"
            assert s == 1, f"r={r} == 3 mod 8 must have v2(ord_r(2)) = 1 (Euler)"
            g = inert_true_rung(r, d, dfac)
            if g != d:
                T["T3_primitive"] = f"fail: true rung {g} (ord_r(omega3) = 4*{g})"
                T["verdict"] = f"attributed to rung {g}; deduplicated"
                if g > 43:
                    dedup_notes.append(
                        f"d={d}: r={r} is non-primitive here; true rung {g} "
                        f"({'in census' if g in adm else 'outside census'})")
                continue
            T["T3_primitive"] = f"pass (ord_r(omega3) = 4*{d})"
            verdict, payload = ord2_oddpart(r)
            if verdict == "one":
                T["T4_p"] = "fail (ord_r(2) = 2, r | 3 impossible)"
                T["verdict"] = "excluded at T4"
                continue
            if verdict == "composite":
                T["T4_p"] = f"fail ({payload})"
                T["verdict"] = "excluded at T4 (odd part of ord_r(2) composite)"
                continue
            if verdict == "unresolved":
                T["T4_p"] = f"UNRESOLVED ({payload})"
                T["verdict"] = "unresolved-survivor of T1-T3 (T4 undecided)"
                unresolved.append({"d": d, "r": str(r), "side": "inert",
                                   "detail": payload})
                log(f"   ?? unresolved inert T1-T3 survivor r={r}")
                continue
            p = payload
            k = (r - 1) // (2 * p)
            assert (r - 1) == 2 * p * k, "k arithmetic"
            T["T4_p"] = f"pass (p = {p}, k = {k})"
            surv.add(r)
            condv = pow(2, p - 2, 3 * d) == 3 * d - 1
            T["T5_condv"] = "pass (d | W_(p-2))" if condv else "fail (d does not divide W_(p-2))"
            log(f"   survivor of (i)-(iv): r = {r}, p = {p}, k = {k}, "
                f"condition (v) {'HOLDS' if condv else 'fails'}")
            if not condv:
                T["verdict"] = "excluded at T5"
                continue
            # full danger triple: phantom vs real (r = W_p iff 3r = 2^p + 1).
            # W_p has bit length p-1; every census r divides V_{2d} <= V_800
            # < 2^1020, so p > 4096 forces r < W_p: exact-compare only small p.
            if p <= 4096:
                is_phantom = (3 * r == 2**p + 1)
            else:
                assert r.bit_length() < p - 1, "r >= W_p size: impossible"
                is_phantom = False
            if is_phantom:
                T["verdict"] = "PHANTOM danger triple (r = W_p, W_p prime)"
                phantoms.append({"d": d, "r": str(r), "p": str(p)})
                log(f"   -> PHANTOM triple ({d}, {r}, {p}): r = W_p")
            else:
                above = r > PLATINUM_CUTOFF
                if p in DISCHARGED_WITNESS_P:
                    st = "discharged-witness"
                elif p == X1_OPEN_P:
                    st = "X1-open"
                elif p in SECONDARY_CLOSED_P:
                    r2 = SECONDARY_CLOSED_P[p]
                    assert isprime(r2) and r2 % 8 == 1 \
                        and pow(2, p, r2) == r2 - 1, \
                        f"secondary-closure factor for p={p} fails verification"
                    st = "secondary-closed"
                else:
                    st = "new-undischarged"
                T["verdict"] = f"REAL danger triple ({st})"
                real_triples.append({"d": d, "r": str(r), "p": str(p),
                                     "k": str(k), "above_platinum_cutoff": above,
                                     "discharge_status": st})
                log(f"   -> REAL danger triple ({d}, {r}, {p}), k={k}, "
                    f"r > 10^12: {above}, status: {st}")
        survivors_by_d[d] = surv
        rung_records.append(rec)
        log("")

    # ------------------------------------------------------------------
    # cross-check (a): reconcile with prop:diagonal-closures at d <= 124
    # ------------------------------------------------------------------
    log("== cross-check (a): prop:diagonal-closures reconciliation (d <= 124) ==")
    diag_checks = []
    recon_ok = True
    for d in sorted(EXPECTED_D124_SURVIVORS):
        exp = EXPECTED_D124_SURVIVORS[d]
        got = survivors_by_d.get(d, set())
        ok = (exp == got)
        recon_ok &= ok
        log(f"   d={d}: survivors of (i)-(iv) expected {sorted(exp)}, got "
            f"{sorted(got)} -> {'MATCH' if ok else 'MISMATCH'}")
        assert ok, f"d={d} survivor set mismatch vs prop:diagonal-closures"
        for r in sorted(got):
            p = EXPECTED_D124_PINNED_P[r]
            condv = pow(2, p - 2, 3 * d) == 3 * d - 1
            diag_checks.append({"d": d, "r": str(r), "p": str(p),
                                "v_holds": bool(condv)})
            log(f"      (d={d}, r={r}, p={p}): condition (v) "
                f"{'HOLDS -> danger triple' if condv else 'FAILS -> not a triple'}")
    log("")

    # reconciliation notes: phantom rung attribution vs corner sweep
    log("== reconciliation: phantom rows vs corner sweep ==")
    r683_rung = inert_true_rung(683, 57, [3, 19])
    log(f"   683 = W_11: primitive rung by T3 = {r683_rung} "
        f"(ord_683(omega3) = {4*r683_rung}).")
    log(f"   683 | V_2e iff e == {r683_rung} (mod {2*r683_rung}): appearances at "
        f"e = 57, 171, 285, 399, ... The cp358 corner sweep counts rectangle")
    log("   hits (non-primitive appearances included): its phantom rows 683@d=57")
    log("   (original d<=124 sweep) and 683@d=171 (extended sweep; 171 = 57+114,")
    log("   admissible) are the SAME prime; this census records it once, at its")
    log(f"   primitive rung {r683_rung}. (285 = 3*5*19 and 399 = 3*7*19 are")
    log("   non-admissible, so no rectangle exists there.)")
    assert r683_rung == 57
    r43691_rung = inert_true_rung(43691, 331, [331])
    log(f"   43691 = W_17: primitive rung by T3 = {r43691_rung}; corner-sweep rows")
    log("   at d=331 (primitive) and d=993 = 331 + 2*331 (non-primitive appearance,")
    log("   outside this census range).")
    assert r43691_rung == 331
    log("   2731 = W_13 has primitive rung 683 > 400: outside this census; its")
    log("   corner-sweep row d=683 is the primitive one.")
    o2731 = inert_true_rung(2731, 683, [683])
    assert o2731 == 683 and V(1366) % 2731 == 0
    log("")

    # ------------------------------------------------------------------
    # cross-check (b): corner-sweep consistency (no REAL triple, p <= 1e8)
    # ------------------------------------------------------------------
    log("== cross-check (b): corner-sweep consistency ==")
    corner_ok = True
    for t in real_triples:
        p = int(t["p"])
        if p <= CORNER_P_LIMIT:
            corner_ok = False
            log(f"   !!! REAL triple with p = {p} <= 10^8 contradicts cp358 corner sweep")
    log(f"   REAL triples with p <= 10^8: "
        f"{[t for t in real_triples if int(t['p']) <= CORNER_P_LIMIT]} "
        f"-> {'PASS (none)' if corner_ok else 'FAIL'}")
    assert corner_ok
    log("")

    # ------------------------------------------------------------------
    # cross-check (c): split-side NCT
    # ------------------------------------------------------------------
    log("== cross-check (c): split-side NCT consistency ==")
    split_ok = (len(split_violations) == 0)
    log(f"   split primes r == 1 (mod 8) satisfying all of "
        f"{{v2(ord)=1, odd part prime, ord_r(alpha)=8d, d | W_(p-2)}}: "
        f"{split_violations if split_violations else 'NONE'} -> "
        f"{'PASS' if split_ok else 'FAIL (contradicts proved NCT d<=400!)'}")
    assert split_ok, "split-side NCT violation -- must not happen"
    log("")

    # ------------------------------------------------------------------
    # APR-CL certification of census survivors
    # ------------------------------------------------------------------
    log("== APR-CL certification (PARI/GP isprime(*,2)) ==")
    have_gp = subprocess.run(["which", "gp"], capture_output=True).returncode == 0
    aprcl_results = {}
    if have_gp:
        to_certify = set()
        for t in real_triples:
            to_certify.add(int(t["r"]))
            to_certify.add(int(t["p"]))
        for ph in phantoms:
            to_certify.add(int(ph["r"]))
        for n in sorted(to_certify):
            ok = aprcl_certify(n)
            aprcl_results[str(n)] = ok
            log(f"   isprime({n}, 2) = {ok}")
            assert ok is True, f"APR-CL failed/timed out for {n}"
    else:
        log("   gp NOT FOUND -- BPSW-only (caveat)")
    log("")

    # ------------------------------------------------------------------
    # summary + outputs
    # ------------------------------------------------------------------
    n_complete = len([r for r in rung_records if r["status"] == "COMPLETE"])
    log("== summary ==")
    log(f"   admissible rungs: {len(adm)}; complete: {n_complete}; "
        f"incomplete: {incomplete if incomplete else 'none'}")
    log(f"   per-rung prime tests (r > 3): {primes_tested}")
    log(f"   REAL danger triples: {len(real_triples)}")
    for t in real_triples:
        log(f"      (d={t['d']}, r={t['r']}, p={t['p']}) k={t['k']} "
            f"above_cutoff={t['above_platinum_cutoff']} {t['discharge_status']}")
    log(f"   PHANTOM triples: {len(phantoms)}")
    for ph in phantoms:
        log(f"      (d={ph['d']}, r={ph['r']}, p={ph['p']})  [r = W_p prime]")
    log(f"   unresolved: {unresolved if unresolved else 'none'}")
    if dedup_notes:
        log("   dedup/attribution notes:")
        for n_ in dedup_notes:
            log(f"      {n_}")
    log(f"   elapsed: {time.time()-t_start:.1f}s")

    out = {
        "generated_by": "scripts/cp412_danger_census_d400.py",
        "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "range": f"43 < d <= {D_MAX}",
        "admissible_rungs": adm,
        "rungs": rung_records,
        "summary": {
            "rungs_complete": n_complete,
            "rungs_incomplete": incomplete,
            "primes_tested": primes_tested,
            "real_triples": real_triples,
            "phantoms": phantoms,
            "unresolved": unresolved,
            "diagonal_survivor_v_checks": diag_checks,
            "split_crosscheck_pass": split_ok,
            "corner_consistency_pass": corner_ok,
            "aprcl": aprcl_results,
            "dedup_notes": dedup_notes,
        },
    }
    with open(OUT_JSON, "w") as fh:
        json.dump(out, fh, indent=1)
    log(f"   wrote {OUT_JSON}")
    log(f"   wrote {OUT_LOG}")


if __name__ == "__main__":
    main()
