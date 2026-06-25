"""cp356: APR-CL flag-2 re-certification of all prime factors of the
six Psi_{4d}, d in {201, 209, 249, 251, 281, 297} (cp353 evidentiary
standard). Factor lists are pulled from FactorDB and verified against
the locally recomputed Psi by exact product identity BEFORE the APR-CL
run (so FactorDB is only a hint source, never trusted bare). Each
factor goes to PARI isprime(r, 2) = APR-CL.
"""
import json
import subprocess
import time
import urllib.parse
import urllib.request

import sympy


def V(n):
    a, b = 2, 2
    for _ in range(n - 1):
        a, b = b, 2 * b + a
    return b


def psi_4d(d):
    odd = sorted(m for m in sympy.divisors(d) if m % 2 == 1)
    psi = {}
    for m in odd:
        val = V(2 * m) // 6
        for mm in odd:
            if mm < m and m % mm == 0 and mm > 1:
                val //= psi[mm]
        psi[m] = val
    return psi[d]


def fdb_factors(n):
    url = "http://factordb.com/api?query=" + urllib.parse.quote(str(n))
    with urllib.request.urlopen(url, timeout=30) as resp:
        data = json.load(resp)
    out = []
    for f, e in data["factors"]:
        out.extend([int(f)] * int(e))
    return out


total = 0
import sys
for d in ([int(x) for x in sys.argv[1:]] or [201, 209, 249, 251, 281, 297]):
    psi = psi_4d(d)
    factors = fdb_factors(psi)
    prod = 1
    for r in factors:
        prod *= r
    assert prod == psi, f"product mismatch d={d}"
    gp_in = "".join(
        f"if(!isprime({r},2), print(\"FAIL {r}\"); quit(1));\n"
        for r in factors) + f'print("d={d}: {len(factors)} factors APR-CL OK");\n'
    res = subprocess.run(["gp", "-q"], input=gp_in, capture_output=True,
                         text=True, timeout=3600)
    print(res.stdout.strip())
    assert "FAIL" not in res.stdout, f"APR-CL failure at d={d}"
    total += len(factors)
    time.sleep(2)
print(f"ALL APR-CL PASS ({total} factors)")
