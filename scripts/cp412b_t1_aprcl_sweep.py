#!/usr/bin/env python3
"""cp412b: APR-CL sweep over EVERY distinct prime factor recorded in the
cp412 danger census (all 34 rung factorizations, r > 3), via PARI/GP
isprime(*, 2).  Purpose: eliminate the BPSW-pseudoprime residual risk in the
rung-completeness certificates (the census APR-CL-certified only the
survivor-related values; the paper's toolchain standard, sec:toolchain, is
APR-CL for every asserted prime factor used in finite exhaustion).

Output: scripts/output/cp412b_t1_aprcl_sweep.log
"""

import json
import os
import subprocess
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CENSUS_JSON = os.path.join(ROOT, "scripts/output/cp412_danger_census.json")
OUT_LOG = os.path.join(ROOT, "scripts/output/cp412b_t1_aprcl_sweep.log")

census = json.load(open(CENSUS_JSON))
primes = sorted({int(pr["r"]) for rec in census["rungs"] for pr in rec["primes"]})
t0 = time.time()

gp_in = ["default(parisizemax, 1G);", "v = [" + ",".join(map(str, primes)) + "];",
         'for(i=1, #v, if(isprime(v[i], 2) != 1, print("COMPOSITE ", v[i])));',
         'print("SWEPT ", #v);']
res = subprocess.run(["nice", "-n", "10", "gp", "-q"], input="\n".join(gp_in),
                     capture_output=True, text=True, timeout=3600)
out = res.stdout.strip()

with open(OUT_LOG, "w") as fh:
    fh.write(f"cp412b APR-CL sweep  ({time.strftime('%Y-%m-%d %H:%M:%S')})\n")
    fh.write(f"distinct primes r > 3 across 34 rungs: {len(primes)}\n")
    fh.write(f"largest: {len(str(primes[-1]))} digits\n")
    fh.write(out + "\n")
    fh.write(f"elapsed: {time.time()-t0:.1f}s\n")
print(out)
print(f"n={len(primes)}, max_digits={len(str(primes[-1]))}, "
      f"elapsed={time.time()-t0:.1f}s")
assert f"SWEPT {len(primes)}" in out and "COMPOSITE" not in out, "APR-CL sweep FAILED"
print("APR-CL sweep PASS: every listed factor is certified prime")
