"""cp356: fixed-d NCT certificates for the six (200,300] targets whose
Psi_{4d} is already fully factored on FactorDB (sweep
cp356_nct_200_300_factordb_sweep.py): d in {201, 209, 249, 251, 281, 297}.

Certificate = cp351/cp353 pattern (Cor. nct-fixed-d-certificate):
  1. recompute Psi_{4d} locally (independent of FactorDB),
  2. verify the product identity against the FactorDB factor list,
  3. BPSW-certify every factor prime (APR-CL flag-2 re-cert separately),
  4. mod-8 filter (only r = 1 mod 8 can be split / triple-relevant),
  5. per split factor: ord_r(2) via explicit descent (r-1 factored via
     FactorDB, then gp/ECM fallback), half-order primality, condition
     (v) d | W_{p-2} via pow(2, p-2, 3d) == 3d-1.
NCT closed at d iff no survivor.

Order descent needs the FULL factorization of r-1; if both FactorDB
and gp (with timeout) fail on a cofactor, the factor is reported
STALLED and the d stays open (no false closure possible).
"""
import json
import subprocess
import sys
import time
import urllib.parse
import urllib.request

import sympy

import sys
TARGETS = [int(x) for x in sys.argv[1:]] or [201, 209, 249, 251, 281, 297]
GP_TIMEOUT = 1800  # seconds per stubborn r-1


def V(n):
    a, b = 2, 2
    for _ in range(n - 1):
        a, b = b, 2 * b + a
    return b


def psi_4d(d):
    odd = sorted(m for m in sympy.divisors(d) if m % 2 == 1)
    psi = {}
    for m in odd:
        val = V(2 * m)
        for mm in odd:
            if mm < m and m % mm == 0:
                val //= psi[mm]
        val //= 6  # Psi_4 = V_2 = 6 divides every V_{2m}, m odd... no:
        psi[m] = val
    return psi[d]


def psi_4d_correct(d):
    """Psi_{4d} via V_{2d} = prod over odd m | d of Psi_{4m}, Psi_4 = 6."""
    odd = sorted(m for m in sympy.divisors(d) if m % 2 == 1)
    psi = {}
    for m in odd:
        val = V(2 * m) // 6  # divide out Psi_4 (m=1 part is V_2 = 6)
        for mm in odd:
            if mm < m and m % mm == 0 and mm > 1:
                val //= psi[mm]
        psi[m] = val
    return psi[d]


def fdb_factors(n, retries=3):
    """FactorDB factor list for n, or None."""
    url = "http://factordb.com/api?query=" + urllib.parse.quote(str(n))
    for _ in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=30) as resp:
                data = json.load(resp)
            if data.get("status") in ("FF", "P", "PRP"):
                out = []
                for f, e in data["factors"]:
                    out.extend([int(f)] * int(e))
                return out
            return None
        except Exception:
            time.sleep(5)
    return None


def gp_factor(n, timeout=GP_TIMEOUT):
    """Full factorization via PARI gp, or None on timeout/failure."""
    try:
        res = subprocess.run(
            ["gp", "-q"], input=f"f=factor({n}); for(i=1,matsize(f)[1], "
            f"print(f[i,1],\" \",f[i,2]))\n",
            capture_output=True, text=True, timeout=timeout)
        out = []
        for line in res.stdout.strip().splitlines():
            p, e = line.split()
            out.extend([int(p)] * int(e))
        if out and sympy.prod(out) == n:
            return out
        return None
    except Exception:
        return None


def factor_full(n, label):
    """Cascade: sympy small -> FactorDB -> gp. Returns list or None."""
    fac = []
    rem = n
    for q in sympy.primerange(2, 100000):
        while rem % q == 0:
            fac.append(q)
            rem //= q
    if rem == 1:
        return fac
    if sympy.isprime(rem):
        return fac + [rem]
    f = fdb_factors(rem)
    if f is None:
        print(f"    [{label}] FactorDB lacks cofactor C{len(str(rem))}; gp/ECM (<= {GP_TIMEOUT}s)...")
        sys.stdout.flush()
        f = gp_factor(rem)
    if f is None:
        return None
    assert sympy.prod(f) == rem
    for q in f:
        assert sympy.isprime(q), f"gp/FactorDB returned composite piece {q}"
    return fac + f


def order2(r, rm1_factors):
    """ord_r(2) by explicit descent from the factorization of r-1."""
    o = r - 1
    for q in set(rm1_factors):
        while o % q == 0 and pow(2, o // q, r) == 1:
            o //= q
    assert pow(2, o, r) == 1
    return o


overall = {}
for d in TARGETS:
    print(f"\n=== d = {d} ===")
    psi = psi_4d_correct(d)
    factors = fdb_factors(psi)
    if factors is None:
        print("  FactorDB no longer FF?! skipping")
        overall[d] = "FDB-MISSING"
        continue
    prod = 1
    for r in factors:
        prod *= r
    assert prod == psi, f"product mismatch at d={d}!"
    print(f"  product of {len(factors)} factors == Psi_{4*d} "
          f"({len(str(psi))} digits): OK")
    for r in factors:
        assert sympy.isprime(r), f"composite factor {r}"
    print("  all factors BPSW-prime: OK")
    survivors, stalled = [], []
    for r in factors:
        cls = r % 8
        if cls != 1:
            print(f"  r ({len(str(r))}d) = {cls} mod 8 -> not split, irrelevant")
            continue
        rf = factor_full(r - 1, f"d={d}, r~10^{len(str(r))-1}")
        if rf is None:
            print(f"  r ({len(str(r))}d) = 1 mod 8: r-1 factorization STALLED")
            stalled.append(r)
            continue
        o = order2(r, rf)
        if o % 2:
            print(f"  r ({len(str(r))}d) split; ord odd -> excluded")
            continue
        p_r = o // 2
        if not sympy.isprime(p_r):
            print(f"  r ({len(str(r))}d) split; half-order composite -> excluded")
            continue
        cond_v = pow(2, p_r - 2, 3 * d) == 3 * d - 1
        print(f"  r ({len(str(r))}d) split; p_r = {p_r} PRIME; "
              f"d | W_(p-2): {cond_v}")
        if cond_v:
            survivors.append((r, p_r))
    if survivors:
        overall[d] = f"SURVIVOR {survivors}"
    elif stalled:
        overall[d] = f"OPEN (stalled r-1: {[len(str(r)) for r in stalled]}-digit)"
    else:
        overall[d] = "CLOSED"
    print(f"  d = {d}: {overall[d]}")

print("\n===== SUMMARY =====")
for d in TARGETS:
    print(f"  d = {d}: {overall.get(d)}")
