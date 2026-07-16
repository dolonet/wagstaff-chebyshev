"""cp354: same-d diagonal exclusion certificates for d in
{59,81,83,99,107,121} from complete FactorDB factorizations of V_{2d}.
A diagonal PPPC pair at d needs TWO distinct primes r1 != r2, both | V_{2d},
both = 3 mod 8, both PRIMITIVE (ord_r(omega3) = 4d), with the SAME pinned
prime p = ord_{r1}(2)/2 = ord_{r2}(2)/2 and d | W_{p-2}. With V_{2d}
completely factored we enumerate ALL such r and group by pinned p:
no group of size >= 2 => S_diag(d) = 0 for ALL p, unconditionally."""
import json, subprocess, sympy

def V(n):
    a, b = 2, 2
    for _ in range(n - 1):
        a, b = b, 2 * b + a
    return b

def fdb(n):
    out = subprocess.run(['curl', '-s', '--max-time', '20',
        f'http://factordb.com/api?query={n}'], capture_output=True, text=True)
    d = json.loads(out.stdout)
    return d['status'], [(int(f), int(m)) for f, m in d['factors']]

def ord_omega3(r, d):
    # exact order of omega3 in F_{r^2}; returns the odd part g with ord = 4g
    R = r
    def mul(a, b):
        return ((a[0]*b[0] + 2*a[1]*b[1]) % R, (a[0]*b[1] + a[1]*b[0]) % R)
    def pw(a, e):
        res = (1, 0)
        while e:
            if e & 1: res = mul(res, a)
            a = mul(a, a); e >>= 1
        return res
    om = (3, 2)
    n = (R + 1) // 4
    if pw(om, 4 * n) != (1, 0):
        return None  # ord doesn't divide r+1 pattern (split r); handle below
    g = n
    for q, e in sympy.factorint(n).items():
        for _ in range(e):
            if pw(om, 4 * (g // q)) == (1, 0): g //= q
            else: break
    return g

all_closed = True
for d in [59, 81, 83, 99, 107, 121]:
    v = V(2 * d)
    status, facs = fdb(v)
    assert status == 'FF', f"d={d} not fully factored: {status}"
    prod = 1
    for f, m in facs:
        prod *= f ** m
    assert prod == v, f"d={d} product mismatch"
    danger = {}
    for f, m in facs:
        if f <= 3: continue
        assert sympy.isprime(f), f"composite listed factor {f}"
        if f % 8 != 3: continue            # need inert
        g = ord_omega3(f, d)
        if g is None or g != d: continue   # need primitivity: ord(omega3) = 4d
        o2 = sympy.n_order(2, f)
        if o2 % 2: continue
        v2 = (o2 & -o2).bit_length() - 1
        if v2 != 1: continue               # need ord = 2 * odd
        p_r = o2 // 2
        if not sympy.isprime(p_r): continue
        cond_v = pow(2, p_r - 2, 3 * d) == 3 * d - 1
        danger.setdefault(p_r, []).append((f, cond_v))
    pairs = {p: rs for p, rs in danger.items() if len(rs) >= 2}
    singles = {p: rs for p, rs in danger.items() if len(rs) == 1}
    print(f"d={d}: V_{2*d} {len(str(v))}d FF; danger-grade primitive inert divisors "
          f"with prime half-order: {sum(len(x) for x in danger.values())}; "
          f"singles={ {p: rs[0][0] for p, rs in singles.items()} }; "
          f"SAME-p PAIRS: {pairs if pairs else 'NONE'}")
    if pairs:
        all_closed = False
print()
print("DIAGONAL CLOSED (S_diag(d) = 0 for ALL p) at every tested d"
      if all_closed else "!!! same-p pair found")
