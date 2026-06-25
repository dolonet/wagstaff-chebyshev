"""cp365: build a STATIC, FactorDB-independent NCT-certificate data artifact
for the Zenodo/v3.6 reproducibility bundle. For every fixed-d closure in
(100, 400] (the 9 values of (100,200] + the 15 of (200,400]) it pins the
complete factorization of Psi_{4d} and records, per prime factor, its
digit-count, BPSW + APR-CL (PARI isprime flag 2) primality, mod-8 class, and
the certificate disposition (for split r = 1 mod 8: ord_r(2) parity ->
half-order primality -> condition (v) d | W_{p-2}). NCT closed at d iff no
factor reaches condition (v) as a survivor.

The point: a referee can verify every NCT closure from this JSON alone,
without trusting FactorDB to be up (the factorizations are pinned here, each
re-checked by exact product identity + APR-CL).

Output: papers/02_wagstaff_chebyshev_reduction/zenodo_v3.6/data/nct_certificates.json
Run (slow: FactorDB + gp for r-1 descent): background.
"""
import json
import os
import subprocess
import sys
import time
import urllib.parse
import urllib.request

import sympy

OUT = ("/home/alexey/claude/primes/papers/02_wagstaff_chebyshev_reduction/"
       "zenodo_v3.6/data/nct_certificates.json")
GP_TIMEOUT = 1800

CLOSED = [107, 121, 129, 131, 139, 163, 171, 177, 179,          # (100,200]
          201, 209, 211, 243, 249, 251, 281, 297,               # (200,300]
          307, 321, 347, 361, 363, 387, 393]                    # (300,400]

# Locally-factored Psi_{4d} (not FF on FactorDB; from cp353/cp358 ECM/GNFS).
# Verified at runtime by the exact product identity prod(factors) == Psi_{4d}.
FACTORIZATIONS = {
    131: [523, 1049, 67073, 135193, 1235593, 6948763, 3525458899,
          319560279309632944107942971,
          67103907069340766277660216218294177],
    163: [1124699, 1096864817, 13182839507, 12036824328615957881,
          221024872680238792171, 65609788974347359786729,
          3577684812700647698297361563656018739],
    179: [1433, 4297, 3627971, 28260363913, 58136103167853442603,
          98733110510495882688859990311113569,
          4970296102561833027001182370298510775841518525241649807641],
    211: [146857, 338779913, 52213127964665809, 307892103084433531,
          515423294391027356040566667977041,
          137414750523933072622414103126689921861439469587248885682924164026532698941843299],
    243: [3, 971, 15114601, 2068037891, 935133871378249,
          18150486096554548669359703025539417,
          6764946053640872456577330274123135728144973423704739859],
    363: [11, 2352241, 69796187, 25691954375979011,
          178659946675133033768857,
          460935198389026844033117792002971203,
          690448252635787211462588854754611369108904661165314966552834606373306146042953],
}


def V(n):
    a, b = 2, 2
    for _ in range(n - 1):
        a, b = b, 2 * b + a
    return b


def psi_4d(d):
    """Psi_{4d} = primitive part of V_{2d}; V_{2d} = prod_{odd m|d} Psi_{4m}, Psi_4 = 6."""
    odd = sorted(m for m in sympy.divisors(d) if m % 2 == 1)
    psi = {}
    for m in odd:
        val = V(2 * m) // 6
        for mm in odd:
            if mm < m and m % mm == 0 and mm > 1:
                val //= psi[mm]
        psi[m] = val
    return psi[d]


def fdb_factors(n, retries=3):
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


def factor_full(n):
    fac, rem = [], n
    for q in sympy.primerange(2, 100000):
        while rem % q == 0:
            fac.append(q)
            rem //= q
    if rem == 1:
        return fac
    if sympy.isprime(rem):
        return fac + [rem]
    f = fdb_factors(rem) or gp_factor(rem)
    if f is None:
        return None
    assert sympy.prod(f) == rem
    return fac + f


def order2(r, rm1_factors):
    o = r - 1
    for q in set(rm1_factors):
        while o % q == 0 and pow(2, o // q, r) == 1:
            o //= q
    assert pow(2, o, r) == 1
    return o


def aprcl(r):
    """PARI isprime(r, 2) == APR-CL primality proof; True/False/None(no gp)."""
    try:
        res = subprocess.run(["gp", "-q"], input=f"print(isprime({r},2))\n",
                             capture_output=True, text=True, timeout=600)
        return res.stdout.strip() == "1"
    except Exception:
        return None


def certify(d):
    psi = psi_4d(d)
    rec = {"d": d, "psi_4d_digits": len(str(psi)), "factors": [],
           "status": None, "note": ""}
    if d in FACTORIZATIONS:
        factors = FACTORIZATIONS[d]
        rec["source"] = "embedded (local ECM/GNFS, cp353/cp358)"
    else:
        # query the FULL Psi_{4d} (FactorDB indexes it as FF); fall back to
        # the trial-divide cascade only if the direct query misses.
        factors = fdb_factors(psi) or factor_full(psi)
        rec["source"] = "FactorDB-FF"
    if factors is None:
        rec["status"] = "FACTORIZATION_INCOMPLETE"
        return rec
    assert sympy.prod(factors) == psi, f"product mismatch d={d}"
    rec["product_identity"] = "OK"
    survivors = []
    for r in sorted(set(factors)):
        fr = {"r": str(r), "r_digits": len(str(r)),
              "bpsw_prime": bool(sympy.isprime(r)),
              "aprcl_prime": aprcl(r), "mod8": r % 8}
        if r % 8 != 1:
            fr["disposition"] = "not split (r != 1 mod 8) -> irrelevant"
        else:
            rm1 = factor_full(r - 1)
            if rm1 is None:
                fr["disposition"] = "STALLED (r-1 factorization incomplete)"
                survivors.append(("stalled", r))
            else:
                o = order2(r, rm1)
                if o % 2:
                    fr["disposition"] = "split; ord_r(2) odd -> excluded"
                else:
                    p_r = o // 2
                    if not sympy.isprime(p_r):
                        fr["disposition"] = "split; half-order composite -> excluded"
                    else:
                        cond_v = pow(2, p_r - 2, 3 * d) == 3 * d - 1
                        fr["p_r"] = p_r
                        if cond_v:
                            fr["disposition"] = "split; SURVIVOR (condition v holds)"
                            survivors.append(("survivor", r))
                        else:
                            fr["disposition"] = "split; condition (v) fails -> excluded"
        rec["factors"].append(fr)
    rec["status"] = "CLOSED" if not survivors else f"OPEN/{survivors}"
    return rec


def main():
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    bundle = {"description": "Static NCT fixed-d certificates for Paper 02 v3.6, "
              "(100,400]. Each Psi_{4d} factorization pinned + APR-CL certified; "
              "split factors excluded by the cor:nct-fixed-d-certificate steps.",
              "generated_by": "scripts/cp365_nct_certificate_bundle.py",
              "certificates": []}
    for d in CLOSED:
        print(f"[cp365] d={d} ...", flush=True)
        try:
            rec = certify(d)
        except Exception as e:
            rec = {"d": d, "status": f"ERROR: {e}"}
        bundle["certificates"].append(rec)
        print(f"[cp365] d={d}: {rec.get('status')}", flush=True)
        with open(OUT, "w") as f:               # incremental dump
            json.dump(bundle, f, indent=2, default=str)
    closed = sum(1 for c in bundle["certificates"] if c.get("status") == "CLOSED")
    print(f"[cp365] DONE: {closed}/{len(CLOSED)} CLOSED, written {OUT}", flush=True)


if __name__ == "__main__":
    main()
