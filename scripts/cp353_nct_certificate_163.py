"""cp353: fixed-d NCT certificate for d = 163 from the ECM full
factorization of Psi_652 (YC campaign, cp352/R3.4). Same certificate as
cp351 (Cor. nct-fixed-d-certificate): verify product, certify primality
(BPSW here; APR-CL flag 2 in PARI below), filter split factors, check
half-order primality and condition (v)."""
import sympy

def V(n):
    a, b = 2, 2
    for _ in range(n - 1):
        a, b = b, 2 * b + a
    return b

d = 163
# Psi_652 = V_326 / Psi_4  (odd divisors of 163: 1, 163; Psi_4 = V_2 = 6)
psi = V(2 * d) // 6
FACTORS = [1124699, 1096864817, 13182839507, 12036824328615957881,
           221024872680238792171, 65609788974347359786729,
           3577684812700647698297361563656018739]
prod = 1
for r in FACTORS:
    prod *= r
assert prod == psi, "product mismatch!"
print(f"product of {len(FACTORS)} factors == Psi_652 ({len(str(psi))} digits): OK")
for r in FACTORS:
    assert sympy.isprime(r), f"composite {r}"
print("all factors BPSW-prime: OK")
survivors = []
for r in FACTORS:
    cls = r % 8
    if cls != 1:
        print(f"  r={r}  = {cls} mod 8 -> not split, irrelevant")
        continue
    o = sympy.n_order(2, r)
    if o % 2:
        print(f"  r={r}  split; ord_r(2)={o} odd -> excluded")
        continue
    p_r = o // 2
    if not sympy.isprime(p_r):
        print(f"  r={r}  split; half-order {p_r} = {sympy.factorint(p_r)} composite -> excluded")
        continue
    cond_v = pow(2, p_r - 2, 3 * d) == 3 * d - 1
    print(f"  r={r}  split; p_r = {p_r} PRIME; d | W_(p-2): {cond_v}")
    if cond_v:
        survivors.append((r, p_r))
print()
print("NCT CLOSED at d = 163" if not survivors else f"SURVIVORS: {survivors}")
