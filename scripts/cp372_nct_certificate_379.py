"""cp372: complete NCT fixed-d certificate for d = 379.

Closure enabled by the Pell half-split identity (cp366) + cp360/cp368 ECM:
Q_758 = 2 (2 P_379 - 1)(2 P_379 + 1). Known factors 4547 (L-half; cp360
trial), 351030947611 (M-half; cp360 ECM p12), 55408612133213734588648530107
and 783398649198311086301698856177 (M-half; cp368 ECM on the C134), leave
prime remainders P141 (L) and P75 (M), completing the factorization of
Psi_1516 (C290):

  Psi_1516 = 4547 * P141  *  351030947611 * p29 * p30 * P75.

Assembled independently in-session (gcd split + trial over 1516k +- 1 up to
2e8 found nothing new; both remainders BPSW-prime). Certificate steps per
cor:nct-fixed-d-certificate, same trust base as cp353/cp358/cp365/cp367
closures: pinned factorization verified by exact product identity, APR-CL
(PARI isprime flag 2) for every factor, factorization-free v2(ord_r(2)) for
split factors, exact-order descent only in the v2 = 1 case.
"""
import json
import math
import os
import subprocess
import sympy

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "output", "cp372_certificate_379.json")

D = 379
KNOWN = [4547, 351030947611,
         55408612133213734588648530107,
         783398649198311086301698856177]


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
    psi = V(2 * D) // 6                     # 379 prime: no nested division
    P, _ = pell(D)
    L, M = 2 * P - 1, 2 * P + 1
    assert V(2 * D) == 2 * L * M
    gL, gM = math.gcd(psi, L), math.gcd(psi, M)
    assert gL * gM == psi
    pL = gL // 4547
    pM = gM // (351030947611 * KNOWN[2] * KNOWN[3])
    factors = sorted(KNOWN + [pL, pM])
    assert sympy.prod(factors) == psi, "product identity fails"
    print(f"[cp372] Psi_1516 C{len(str(psi))} = product of {len(factors)} "
          f"factors, identity OK")

    rec = {"d": D, "psi_4d_digits": len(str(psi)),
           "source": "cp360 trial+ECM + cp366 Pell half-split "
                     "(L=2P_379-1, M=2P_379+1) + cp368 ECM on the C134",
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
        print(f"[cp372] r ({fr['r_digits']}d, mod8={fr['mod8']}, "
              f"half={fr['half']}, aprcl={fr['aprcl_prime']}): "
              f"{fr['disposition']}", flush=True)
        rec["factors"].append(fr)
    rec["status"] = "CLOSED" if not survivors else f"OPEN/{survivors}"
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(rec, f, indent=2)
    print(f"[cp372] d=379: {rec['status']} -> {OUT}")


if __name__ == "__main__":
    main()
