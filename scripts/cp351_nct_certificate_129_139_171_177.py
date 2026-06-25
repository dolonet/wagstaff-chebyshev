"""cp351: fixed-d NCT certificates for d in {129, 139, 171, 177} via the
complete FactorDB factorizations of the primitive parts Psi_{4d}.

Certificate (Paper B v3.3, Cor. nct-fixed-d-certificate):
  1. verify  prod(factors) == Psi_{4d}  and each factor is prime
     (sympy BPSW here; APR-CL via PARI in cp351_aprcl_verify.gp);
  2. keep split factors r = 1 mod 8;
  3. p_r = ord_r(2)/2 must be prime, else no compatible triple at r;
  4. if p_r prime, test d | W_{p_r-2}  <=>  2^(p_r-2) = -1 mod 3d.
No surviving factor => NCT holds at this d.
"""
import sympy

def V(n):
    a, b = 2, 2
    for _ in range(n - 1):
        a, b = b, 2 * b + a
    return b if n >= 1 else 2

def psi_parts(dmax_set):
    needed = sorted(set(m for d in dmax_set for m in sympy.divisors(d) if m % 2 == 1))
    psi = {}
    for d in needed:
        val = V(2 * d)
        for m in sympy.divisors(d):
            if m % 2 == 1 and m < d:
                val //= psi[m]
        psi[d] = val
    return psi

FACTORS = {
    129: [75390435818587404514891047039001, 276431598001487177038457370080041],
    139: [6568027, 574593323881,
          197671787019026728856728443573405844788899,
          5764333243492029733503944458218529327657872457],
    171: [8209, 10259, 1186057, 82658931465203, 300531150542299499,
          19269878832774783314803453706733003299],
    177: [6362089, 21015437931105774861319763619250473457,
          490240983802792962966192466615950084458511337],
}

psi = psi_parts(FACTORS.keys())
all_ok = True
for d, facs in FACTORS.items():
    print(f"=== d = {d}  (Psi_{4*d}, {len(str(psi[d]))} digits) ===")
    prod = 1
    for r in facs:
        prod *= r
    assert prod == psi[d], f"product mismatch at d={d}"
    print(f"  product of {len(facs)} factors == Psi_{4*d}: OK")
    for r in facs:
        assert sympy.isprime(r), f"composite factor {r}"
    print(f"  all factors BPSW-prime: OK")
    survivors = []
    for r in facs:
        cls = r % 8
        if cls != 1:
            print(f"  r={r}  = {cls} mod 8  -> not split, irrelevant")
            continue
        o = sympy.n_order(2, r)
        if o % 2 != 0:
            print(f"  r={r}  split; ord_r(2)={o} odd -> excluded")
            continue
        p_r = o // 2
        if not sympy.isprime(p_r):
            print(f"  r={r}  split; ord_r(2)/2 = {p_r} = "
                  f"{sympy.factorint(p_r)} composite -> excluded")
            continue
        # exact condition d | W_{p_r - 2}  <=> 2^(p_r-2) == -1 mod 3d
        cond_v = pow(2, p_r - 2, 3 * d) == 3 * d - 1
        print(f"  r={r}  split; p_r = {p_r} PRIME; d | W_(p-2): {cond_v}")
        if cond_v:
            survivors.append((r, p_r))
    if survivors:
        all_ok = False
        print(f"  !!! COMPATIBLE TRIPLE CANDIDATE(S): {survivors}")
    else:
        print(f"  NCT CLOSED at d = {d}: no compatible triple exists.")
print()
print("ALL FOUR d CLOSED" if all_ok else "SURVIVORS FOUND — NOT CLOSED")
