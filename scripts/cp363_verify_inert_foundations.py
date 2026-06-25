"""
cp361: Independent adversarial verification of the inert-branch foundations
of Paper B (wagstaff_chebyshev_reduction_v3.6.tex). Lean/incremental version.

Re-derives and numerically checks on REAL Wagstaff inert factors:
  thm:lte, thm:pell-divisibility, prop:pell-closed-form, lem:inert-primitivity,
  thm:pell-exclusion, thm:parity-obstruction, lem:classIII-division.
Census 845/869 is confirmed separately by scripts/cp152_inert_G4_scan.py.

Enumeration capped at r < ENUM_LIMIT for speed (spot-check, not full census).
"""
import sys
import sympy
from sympy import factorint, isprime, primerange, n_order

ENUM_LIMIT = 60000
PR = lambda *a: (print(*a), sys.stdout.flush())

def pell_UV(n):
    U0, U1, V0, V1 = 0, 1, 2, 2
    if n == 0: return 0, 2
    if n == 1: return 1, 2
    for _ in range(2, n + 1):
        U0, U1 = U1, 2 * U1 + U0
        V0, V1 = V1, 2 * V1 + V0
    return U1, V1

def ord2(r):
    return int(n_order(2, r))

def mul2(x, y, r):
    (a, b), (c, d) = x, y
    return ((a * c + 2 * b * d) % r, (a * d + b * c) % r)

def pw2(x, e, r):
    res = (1, 0)
    while e > 0:
        if e & 1: res = mul2(res, x, r)
        x = mul2(x, x, r); e >>= 1
    return res

def ord_omega3(r):
    om = (3, 2); o = r + 1
    for q in factorint(r + 1):
        while o % q == 0 and pw2(om, o // q, r) == (1, 0):
            o //= q
    return o

def W(n): return (pow(2, n) + 1) // 3
def vq(x, q):
    v = 0
    while x % q == 0: x //= q; v += 1
    return v

# Class III membership
def in_A(q):
    if q == 3: return True
    o = ord2(q)
    return o % 4 == 2

# ---- build A and check |A| ----
A = set([3])
for q in primerange(5, 5000):
    if in_A(q): A.add(q)
PR(f"|A cap [3,5000)| Class III primes = {len(A)} (paper says 194 with q<5000; includes 3? -> {len(A)})")

# ============ CHECK 1: LTE exactness ============
PR("\n=== CHECK 1: thm:lte v_q(W_n) exact (off-by-one audit) ===")
fail = 0; aq_fail = 0
for n in range(3, 600, 2):
    if vq(W(n), 3) != vq(n, 3): fail += 1; PR(f"  q=3 FAIL n={n}")
for q in sorted(A - {3})[:40]:
    m = ord2(q) // 2
    aq = vq(pow(2, m) + 1, q)
    if aq != 1: aq_fail += 1; PR(f"  alpha_q!=1: q={q} a_q={aq}")
    for j in range(1, 12, 2):  # n=m*j, j odd
        n = m * j
        if n % 2 == 0: continue
        if vq(W(n), q) != aq + vq(j, q): fail += 1; PR(f"  q={q} FAIL n={n}")
PR(f"  LTE failures={fail}; alpha_q!=1 count={aq_fail}")

# ============ ENUMERATE inert factors r<ENUM_LIMIT ============
PR(f"\n=== Enumerating inert factors r<{ENUM_LIMIT} (spot-check subset) ===")
fac_list = []
for r in range(11, ENUM_LIMIT, 8):
    if not isprime(r): continue
    o2 = ord2(r)
    if o2 % 2: continue
    p = o2 // 2
    if not (isprime(p) and p >= 5): continue
    Wp = W(p)
    if Wp == r or Wp <= r or Wp % r: continue  # composite W_p, r proper
    oo = ord_omega3(r)
    assert oo % 4 == 0
    g = oo // 4
    Wpm2 = W(p - 2)
    Gr = int(sympy.gcd((r + 1) // 4, Wpm2))
    fac_list.append((r, p, g, Gr, Wpm2))
PR(f"  inert factors found (r<{ENUM_LIMIT}): {len(fac_list)}")
PR(f"  distinct g: {sorted(set(g for _,_,g,_,_ in fac_list))}")
PR(f"  distinct G_r: {sorted(set(G for _,_,_,G,_ in fac_list))}")
PR(f"  G_r<=43: {sum(1 for *_,G,_ in [(r,p,g,G,W) for (r,p,g,G,W) in fac_list])}" )
le43 = sum(1 for (r,p,g,G,W) in fac_list if G <= 43)
PR(f"  G_r<=43 count = {le43} / {len(fac_list)} (these are the cor:small-Gr blocked ones)")

# ============ CHECK 2: pell-divisibility + closed form on real r ============
PR("\n=== CHECK 2: thm:pell-divisibility + prop:pell-closed-form ===")
pd = cf = 0; tested = 0
for (r, p, g, Gr, Wpm2) in fac_list:
    d = g
    if d > 350: continue  # keep Pell numbers small
    U4, V4 = pell_UV(4 * d); U2, V2 = pell_UV(2 * d); U1, V1 = pell_UV(d)
    gg = int(sympy.gcd((V4 + 2) // 2, U4))
    if gg % r != 0: pd += 1; PR(f"  PD FAIL r={r} d={d}")
    if not (gg == V2 == V1**2 + 2): cf += 1; PR(f"  CF FAIL d={d}")
    tested += 1
PR(f"  tested {tested}; pell-div fails={pd}, closed-form fails={cf}")

# ============ CHECK 3: inert-primitivity ============
PR("\n=== CHECK 3: lem:inert-primitivity (r|V_2d, r∤V_2e for 0<e<d) ===")
pf = 0; tested = 0
for (r, p, g, Gr, Wpm2) in fac_list:
    d = g
    if d > 300: continue
    _, V2d = pell_UV(2 * d)
    if V2d % r: pf += 1; PR(f"  PRIM r∤V2d r={r} d={d}"); continue
    for e in range(1, d):
        _, V2e = pell_UV(2 * e)
        if V2e % r == 0: pf += 1; PR(f"  PRIM early-div r={r} d={d} e={e}"); break
    tested += 1
PR(f"  tested {tested}; primitivity fails={pf}")

# ============ CHECK 4: pell-exclusion admissible set + per-d ============
PR("\n=== CHECK 4: thm:pell-exclusion ===")
adm_claim = {1, 3, 9, 11, 19, 27, 33, 43}
recomp = set(d for d in range(1, 44, 2) if all(in_A(q) for q in factorint(d)))
PR(f"  admissible claimed   : {sorted(adm_claim)}")
PR(f"  admissible recomputed: {sorted(recomp)}  MATCH={adm_claim==recomp}")
for d in sorted(recomp):
    _, Vd = pell_UV(d)
    Nd = Vd**2 + 2
    fac = factorint(Nd)
    danger = []
    surv = []
    for r in fac:
        if r == 3 or r % 8 != 3: continue
        o2 = ord2(r)
        if o2 % 4 != 2: continue          # cond 1: v2(ord)!=1 excluded
        m = o2 // 2
        if not isprime(m): continue       # cond 2: m not prime excluded
        pp = m
        Wpp = W(pp)
        if Wpp == r: surv.append((r, pp, "phantom-Wp-prime")); continue
        if W(pp - 2) % d == 0: danger.append((r, pp)); surv.append((r,pp,"DANGER"))
        else: surv.append((r, pp, "excluded-cond3 d∤W_{p-2}"))
    PR(f"  d={d}: N_d={Nd} ({len(fac)} primefac); r≡3mod8 survivors past c1,c2: {surv}; REAL DANGER={danger}")

# ============ CHECK 5: parity-obstruction on real factors ============
PR("\n=== CHECK 5: thm:parity-obstruction (g has q∉A => g∤W_{p-2}) ===")
contra = 0; passers = []
for (r, p, g, Gr, Wpm2) in fac_list:
    fac = factorint(g)
    nonA = any(not in_A(q) for q in fac)
    div = (Wpm2 % g == 0)
    if nonA and div:
        contra += 1; PR(f"  CONTRADICTION r={r} p={p} g={g}: nonClassIII factor but g|W_{{p-2}}")
    if div and g > 1: passers.append((r, p, g))
PR(f"  parity contradictions = {contra}")
PR(f"  inert factors with g>1 and g|W_{{p-2}} (Cond II passers) = {len(passers)} -> {passers[:8]}")

# ============ CHECK 6: classIII-division lem ============
PR("\n=== CHECK 6: lem:classIII-division (q>3 in A: q|W_{p-2} iff m_q|p-2) ===")
cd = 0; n = 0
for q in sorted(A - {3})[:30]:
    mq = ord2(q) // 2
    for p in primerange(5, 3000):
        lhs = (W(p - 2) % q == 0)
        rhs = ((p - 2) % mq == 0)
        if lhs != rhs: cd += 1; PR(f"  CD FAIL q={q} p={p}")
        n += 1
PR(f"  tested {n} (q,p) pairs; classIII-division fails={cd}")
PR("\n=== DONE ===")
