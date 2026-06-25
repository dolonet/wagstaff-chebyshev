"""cp358: fixed-d NCT certificate for d = 363 from the local ECM full
factorization of Psi_1452 (gmp-ecm 7.0.6 driver, B1 <= 11e6, ONE round — the C169 fell instantly; ECM-friendliness lesson #5).
Certificate = cp351/cp353/cp356 pattern."""
import sympy


def V(n):
    a, b = 2, 2
    for _ in range(n - 1):
        a, b = b, 2 * b + a
    return b


d = 363
odd = [m for m in sympy.divisors(d) if m % 2 == 1]
psi = {}
for m in odd:
    val = V(2 * m) // 6
    for mm in odd:
        if mm < m and m % mm == 0 and mm > 1:
            val //= psi[mm]
    psi[m] = val
PSI = psi[d]

FACTORS = [11, 2352241, 69796187, 25691954375979011,
           178659946675133033768857,
           460935198389026844033117792002971203,
           690448252635787211462588854754611369108904661165314966552834606373306146042953]
prod = 1
for r in FACTORS:
    prod *= r
assert prod == PSI, "product mismatch!"
print(f"product of {len(FACTORS)} factors == Psi_1452 "
      f"({len(str(PSI))} digits): OK")
for r in FACTORS:
    assert sympy.isprime(r), f"composite {r}"
print("all factors BPSW-prime: OK")
survivors = []
for r in FACTORS:
    cls = r % 8
    if cls != 1:
        print(f"  r ({len(str(r))}d) = {cls} mod 8 -> not split, irrelevant")
        continue
    o = sympy.n_order(2, r)
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
print()
print("NCT CLOSED at d = 363" if not survivors
      else f"SURVIVORS: {survivors}")
