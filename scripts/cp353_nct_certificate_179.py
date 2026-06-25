"""cp353: fixed-d NCT certificate for d = 179 from the complete ECM
factorization of Psi_716 (YC campaign; the B1=11e6 stage finished the
number: final prime cofactor P58)."""
import sympy

def V(n):
    a, b = 2, 2
    for _ in range(n - 1):
        a, b = b, 2 * b + a
    return b

d = 179
psi = V(2 * d) // 6
FACTORS = [1433, 4297, 3627971, 28260363913, 58136103167853442603,
           98733110510495882688859990311113569,
           4970296102561833027001182370298510775841518525241649807641]
prod = 1
for r in FACTORS:
    prod *= r
assert prod == psi, "product mismatch"
print(f"product of {len(FACTORS)} factors == Psi_716 ({len(str(psi))} digits): OK")
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
print("NCT CLOSED at d = 179" if not survivors else f"!!! SURVIVORS: {survivors}")
