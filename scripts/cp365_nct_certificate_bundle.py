"""cp365: build a STATIC, FactorDB-independent NCT-certificate data artifact
for the Zenodo/v3.7 reproducibility bundle. For every fixed-d closure in
(100, 400] (the 9 values of (100,200] + all 19 of (200,400]; 28 total since
cp377) it pins the complete factorization of Psi_{4d} and records, per prime
factor, its
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

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "..", "data", "nct_certificates.json")
GP_TIMEOUT = 1800

CLOSED = [107, 121, 129, 131, 139, 163, 171, 177, 179,          # (100,200]
          201, 209, 211, 227, 243, 249, 251, 281, 283, 297,     # (200,300]
          307, 321, 331, 347, 361, 363, 379, 387, 393]          # (300,400]

# Locally-factored Psi_{4d} (not FF on FactorDB; from cp353/cp358 ECM/GNFS,
# and cp366-cp377 half-split + ECM + GNFS for d = 227/283/331/379).
# Verified at runtime by the exact product identity prod(factors) == Psi_{4d}.
FACTORIZATIONS = {
    227: [907,
          57203,
          5580569,
          16552486912925363,
          103897484432089424574671548054573458540642269367990262779403,
          201760911884372460299372338734188130826364691502781503178121274772114731708580885801],
    283: [6793,
          12451,
          11489975459,
          21850712277340062427,
          806304582166957550068693465373407792418260097,
          49679742606807434161308646677418138917837852718337698366097,
          877319209569244846810168253907948677873416655650083562956858621782177415571],
    331: [6619,
          43691,
          4510867,
          224936385721,
          1732445993139049817801,
          5264282380581441760369,
          104620710651125217917754500294169403,
          19496626790588131600768193365585307203602374799879171175562014606473433,
          762557006971255766502991992708120614676358835210615414018541037281740581273457],
    379: [4547,
          351030947611,
          55408612133213734588648530107,
          783398649198311086301698856177,
          547727124112532774042307292211984248239243330440030687686315134398603465411,
          611819812080262238556029279191272259222584574755509428287193285422780631900357459059587512290065179360322811461486058125545862408366246029937],
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


def v2_ord2(r):
    """v2(ord_r(2)), factorization-free (repeated squaring; cp372 method).
    Lets the certificate dispose of v2 != 1 split factors without factoring
    r - 1 — essential for the P71-P141 split factors of d = 283/331/379."""
    e2 = ((r - 1) & -(r - 1)).bit_length() - 1
    t = pow(2, (r - 1) >> e2, r)
    if t == 1:
        return 0
    for j in range(1, e2 + 1):
        t = t * t % r
        if t == 1:
            return j
    raise AssertionError("2^(r-1) != 1 mod r -- r not prime?")


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
            # factorization-free order-parity first (cp372 method): only the
            # v2 == 1 case needs the exact order, hence an r-1 factorization.
            v2 = v2_ord2(r)
            fr["v2_ord2"] = v2
            if v2 == 0:
                fr["disposition"] = "split; ord_r(2) odd -> excluded"
            elif v2 >= 2:
                fr["disposition"] = ("split; v2(ord_r(2)) >= 2 -> half-order "
                                     "even -> composite -> excluded")
            else:
                rm1 = factor_full(r - 1)
                if rm1 is None:
                    fr["disposition"] = "STALLED (r-1 factorization incomplete)"
                    survivors.append(("stalled", r))
                else:
                    o = order2(r, rm1)
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
    bundle = {"description": "Static NCT fixed-d certificates for Paper 02 (v3.7), "
              "(100,400] complete. Each Psi_{4d} factorization pinned + APR-CL "
              "certified; split factors excluded by the "
              "cor:nct-fixed-d-certificate steps.",
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
