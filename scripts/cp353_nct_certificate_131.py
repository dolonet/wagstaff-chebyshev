"""cp353: fixed-d NCT certificate for d = 131 from the cado-nfs (GNFS)
complete factorization of Psi_524 (YC campaign). Certificate logic as in
cp351/cp353_163."""
import sympy

def V(n):
    a, b = 2, 2
    for _ in range(n - 1):
        a, b = b, 2 * b + a
    return b

d = 131
psi = V(2 * d) // 6  # odd divisors of 131: {1, 131}; Psi_4 = 6
FACTORS = [523, 1049, 67073, 135193, 1235593, 6948763, 3525458899,
           319560279309632944107942971,
           67103907069340766277660216218294177]
prod = 1
for r in FACTORS:
    prod *= r
assert prod == psi, f"product mismatch! prod has {len(str(prod))} digits, psi {len(str(psi))}"
print(f"product of {len(FACTORS)} factors == Psi_524 ({len(str(psi))} digits): OK")
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
print("NCT CLOSED at d = 131" if not survivors else f"!!! SURVIVORS: {survivors}")
