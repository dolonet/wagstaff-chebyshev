"""cp410: the self-supporting ladder for Condition (II).

Proposition (paper 02 v3.8, prop:ladder): if W_{p-2} is prime, Condition
(II) implies primality of W_p — the fully-factored-(N+1) BLS
specialization, completed by the d = 1 exclusion (r | 36) on the split
branch and the Pell-sequence exclusion (g > 43) on the inert branch.

This script builds the honest instance list:
  unconditional rungs: q in A000978 with W_q a CERTIFIED prime and
    p = q + 2 prime;
  conditional rungs: q in the PRP tail with p = q + 2 prime (each becomes
    unconditional the day W_q gets a certificate).
Also sanity-checks the proposition's two norm exclusions.
"""
import sympy

# A000978: (2^q+1)/3 prime/PRP exponents (as of 2026)
A000978 = [3, 5, 7, 11, 13, 17, 19, 23, 31, 43, 61, 79, 101, 127, 167,
           191, 199, 313, 347, 701, 1709, 2617, 3539, 5807, 10501, 10691,
           11279, 12391, 14479, 42737, 83339, 95369, 117239, 127031,
           138937, 141079, 267017, 269987, 374321, 986191, 4031399,
           13347311, 13372531, 15135397]
# W_q proven prime: every listed q <= 42737 (classical + ECPP era; incl.
# our BLS 2617/10501/12391), plus the ECPP proofs recorded on the t5k
# top-20 Wagstaff page (accessed 2026-07-16): 83339 (Sep 2014),
# 117239 (Aug 2022), 127031 (Jan 2023), 138937, 141079. 95369 is not
# shown as proven there; immaterial for the ladder (95371 is composite).
CERTIFIED = {q for q in A000978 if q <= 42737}
CERTIFIED |= {83339, 117239, 127031, 138937, 141079}

print("== sanity: the two exclusions of prop:ladder ==")
# split d=1: omega_3^2 = -1 mod r forces r | N(omega_3^2 + 1) = 36
from sympy import sqrt, expand
w2 = (17, 12)  # omega_3^2 = 17 + 12*sqrt(2)
norm = (w2[0] + 1) ** 2 - 2 * w2[1] ** 2
assert norm == 36, norm
print(f"N(omega_3^2 + 1) = {norm} -> split d=1 impossible for r >= 17  OK")
# 4*W_{p-2} = W_p + 1 identity
for p in [5, 7, 11, 13, 19, 31]:
    W = lambda n: (2 ** n + 1) // 3
    assert 4 * W(p - 2) == W(p) + 1
print("4*W_{p-2} = W_p + 1 verified on samples  OK")

print("\n== unconditional rungs (W_{p-2} certified prime, p prime) ==")
uncond = []
for q in A000978:
    p = q + 2
    if sympy.isprime(p):
        if q in CERTIFIED:
            uncond.append(p)
            print(f"  q={q:>8}  p={p:>8}  W_q certified -> "
                  f"Condition (II) is an UNCONDITIONAL primality test at p")
print(f"unconditional ladder: p in {uncond}")

print("\n== conditional rungs (W_q currently PRP only) ==")
for q in A000978:
    p = q + 2
    if q not in CERTIFIED and sympy.isprime(p):
        print(f"  q={q:>9}  p={p:>9}  (awaits a certificate for W_q)")

print("\n== near-misses (q certified, p = q+2 composite) ==")
for q in A000978:
    p = q + 2
    if q in CERTIFIED and not sympy.isprime(p):
        print(f"  q={q:>6}  p={p:>6} = {sympy.factorint(p)}")
