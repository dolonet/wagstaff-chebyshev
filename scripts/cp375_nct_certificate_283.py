"""cp375: complete NCT fixed-d certificate for d = 283.

Closure chain: cp360 trial+ECM on the unsplit C208 (factors 6793, 12451,
11489975459, 21850712277340062427), cp366 Pell half-split
Q_566 = 2 (2 P_283 - 1)(2 P_283 + 1), cp368 ECM t35 on the L-half residual
C104 (no factor), cp374 GNFS on the C104 (cado-nfs, this script reads the
factors from scripts/output/cp374_gnfs_results.txt line "d283_L_C104 ...").
The M-half residual is BPSW/APR-CL prime (P75), so with the GNFS split
Psi_1132 is completely factored:

  Psi_1132 = 12451 * f1 * f2  *  6793 * 11489975459 * 21850712277340062427 * P75
             (L-half, f1*f2 = C104)   (M-half)

(The rational factor 3 of L = 2 P_283 - 1 is not in Psi_1132 = V_566/6.)
Certificate steps per cor:nct-fixed-d-certificate, same trust base as the
cp353/cp358/cp365/cp367/cp372 closures: pinned factorization verified by
exact product identity, APR-CL (PARI isprime flag 2) for every factor,
factorization-free v2(ord_r(2)) for split factors, exact-order descent only
in the v2 = 1 case.
"""
import json
import math
import os
import subprocess
import sympy

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "output", "cp375_certificate_283.json")
GNFS_RESULTS = os.path.join(HERE, "output", "cp374_gnfs_results.txt")

D = 283
C104 = 40057004104743866669602345140100256037457556974869233878102970879973584214152054621594927825184472731409
KNOWN_L = [12451]
KNOWN_M = [6793, 11489975459, 21850712277340062427]


def gnfs_factors(label, target):
    """Parse 'label N f1 f2 ...' from the cp374 results file."""
    with open(GNFS_RESULTS) as f:
        for line in f:
            tok = line.split()
            if len(tok) >= 4 and tok[0] == label:
                N = int(tok[1])
                fs = [int(x) for x in tok[2:]]
                assert N == target, "cp374 line does not match the pinned target"
                assert math.prod(fs) == N, "GNFS factor product mismatch"
                return fs
    raise SystemExit(f"no '{label}' line in {GNFS_RESULTS} — run cp374 first")


def pell(n):
    if n == 0:
        return 0, 2
    P, Q = pell(n // 2)
    s = -2 if n // 2 % 2 else 2
    P2, Q2 = P * Q, Q * Q - s
    if n % 2 == 0:
        return P2, Q2
    P1 = P2 + Q2 // 2
    Q1 = 2 * P1 + 2 * P2
    return P1, Q1


def V(n):
    a, b = 2, 2
    for _ in range(n - 1):
        a, b = b, 2 * b + a
    return b


def aprcl(r):
    res = subprocess.run(["gp", "-q"], input=f"print(isprime({r},2))\n",
                         capture_output=True, text=True, timeout=7200)
    return res.stdout.strip() == "1"


def v2_ord2(r):
    """v2(ord_r(2)), factorization-free (repeated squaring)."""
    e2 = ((r - 1) & -(r - 1)).bit_length() - 1
    t = pow(2, (r - 1) >> e2, r)
    if t == 1:
        return 0
    for j in range(1, e2 + 1):
        t = t * t % r
        if t == 1:
            return j
    raise AssertionError("2^(r-1) != 1 mod r -- r not prime?")


def main():
    fs = gnfs_factors("d283_L_C104", C104)
    psi = V(2 * D) // 6                     # 283 prime: no nested division
    P, _ = pell(D)
    L, M = 2 * P - 1, 2 * P + 1
    assert V(2 * D) == 2 * L * M
    gL, gM = math.gcd(psi, L), math.gcd(psi, M)
    assert gL * gM == psi
    assert gL == math.prod(KNOWN_L) * C104, "L-half known-factor ledger wrong"
    pM = gM // math.prod(KNOWN_M)
    factors = sorted(KNOWN_L + KNOWN_M + [pM] + fs)
    assert math.prod(factors) == psi, "product identity fails"
    print(f"[cp375] Psi_1132 C{len(str(psi))} = product of {len(factors)} "
          f"factors, identity OK")

    rec = {"d": D, "psi_4d_digits": len(str(psi)),
           "source": "cp360 trial+ECM + cp366 Pell half-split "
                     "(L=2P_283-1, M=2P_283+1) + cp368 ECM t35 (no factor) "
                     "+ cp374 GNFS on the L-half C104 (cado-nfs)",
           "product_identity": "OK", "factors": [], "status": None}
    survivors = []
    for r in factors:
        fr = {"r": str(r), "r_digits": len(str(r)),
              "bpsw_prime": bool(sympy.isprime(r)),
              "aprcl_prime": aprcl(r), "mod8": r % 8,
              "half": "L" if L % r == 0 else "M"}
        assert fr["bpsw_prime"] and fr["aprcl_prime"], f"factor {r} not prime!"
        if r % 8 != 1:
            fr["disposition"] = "not split (r != 1 mod 8) -> irrelevant"
        else:
            v2 = v2_ord2(r)
            fr["v2_ord2"] = v2
            if v2 == 0:
                fr["disposition"] = "split; ord_r(2) odd -> excluded"
            elif v2 >= 2:
                fr["disposition"] = ("split; v2(ord_r(2)) >= 2 -> half-order "
                                     "even -> composite -> excluded")
            else:
                # v2 == 1: need exact order (descent)
                o = sympy.n_order(2, r)
                p_r = o // 2
                if not sympy.isprime(p_r):
                    fr["disposition"] = "split; half-order composite -> excluded"
                else:
                    cond_v = pow(2, p_r - 2, 3 * D) == 3 * D - 1
                    fr["p_r"] = int(p_r)
                    fr["disposition"] = ("split; SURVIVOR (condition v holds)"
                                         if cond_v else
                                         "split; condition (v) fails -> excluded")
                    if cond_v:
                        survivors.append(r)
        print(f"[cp375] r ({fr['r_digits']}d, mod8={fr['mod8']}, "
              f"half={fr['half']}, aprcl={fr['aprcl_prime']}): "
              f"{fr['disposition']}", flush=True)
        rec["factors"].append(fr)
    rec["status"] = "CLOSED" if not survivors else f"OPEN/{survivors}"
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(rec, f, indent=2)
    print(f"[cp375] d=283: {rec['status']} -> {OUT}")


if __name__ == "__main__":
    main()
