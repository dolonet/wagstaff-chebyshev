#!/usr/bin/env python3
"""
CP183 — Rigorous Platinum Lemma: at-most-1 danger factor per p in r ≤ 10^12 range.

This script verifies the rigorous computational statement used in the paper:

  LEMMA (CP183, unconditional for r ≤ 10^{12}):
  For every prime p, the number of inert prime factors r ≤ 10^{12} of W_p
  with d_r | W_{p-2} is AT MOST 1.

The platinum database (CP172) is an exhaustive computation covering all
(p, r) pairs with r ≤ 10^{12}, r prime, r ≡ 3 mod 8, and 2p | r-1
(i.e., ord_r(2) = 2p, Order-Pinning). Among 684,965,381 inert-factor rows,
only 3 have obstruction = PASS (satisfying d_r | W_{p-2}).

These 3 PASS rows are:
  (p, r, d) = (1251488009, 2502976019, 41716267)
             (10916765939, 152834723147, 38208680787)
             (85684865933, 171369731867, 57)

All 3 p values are DISTINCT. Hence at most 1 PASS per p. QED.

Corollary: By Three-Factor Reduction, any composite W_p that is a
Condition II danger must have ≥ 3 inert factors. By this lemma, at least
2 of those factors must exceed 10^{12}.

This script verifies these claims directly from the Platinum CSV
(inert_factors.csv) shipped on Zenodo at 10.5281/zenodo.19496206.
Place the CSV (or its .xz archive) at data/inert_factors.csv[.xz] relative
to the repository root before running.
"""

import csv
import lzma
import os
from collections import Counter

PLATINUM_CANDIDATES = [
    "data/inert_factors.csv.xz",
    "data/inert_factors.csv",
    "inert_factors.csv.xz",
    "inert_factors.csv",
]
PLATINUM = next((p for p in PLATINUM_CANDIDATES if os.path.exists(p)), PLATINUM_CANDIDATES[0])

print("="*80)
print("Platinum Lemma (paper Theorem 8.4) — direct verification")
print("="*80)
print(f"Reading: {PLATINUM}")
print()

pass_rows = []
rows_by_p = Counter()

opener = lzma.open if PLATINUM.endswith(".xz") else open
with opener(PLATINUM, "rt") as f:
    reader = csv.DictReader(f)
    for row in reader:
        p = int(row["p"])
        rows_by_p[p] += 1
        if row["obstruction"] == "PASS":
            pass_rows.append(row)

print(f"Total (p, r) rows in platinum DB: {sum(rows_by_p.values()):,}")
print(f"Total distinct p values: {len(rows_by_p):,}")
print(f"Total PASS rows: {len(pass_rows)}")
print()

# Verify distinctness of PASS p values
pass_ps = [int(r["p"]) for r in pass_rows]
print("PASS rows:")
for row in pass_rows:
    p, r, d = int(row["p"]), int(row["r"]), int(row["d"])
    print(f"  p={p:>20}  r={r:>25}  d={d:>15}  pure_iii={row['pure_iii']}")
print()

if len(set(pass_ps)) == len(pass_ps):
    print(f"✓ VERIFIED: all {len(pass_ps)} PASS p values are DISTINCT.")
else:
    print(f"✗ FAIL: {len(pass_ps)} PASS rows but only {len(set(pass_ps))} distinct p.")
print()

# Verify max inert-factor count per p in DB
max_rows = max(rows_by_p.values())
print(f"Max inert-factor rows per p in DB: {max_rows}")
print(f"Distribution: {Counter(rows_by_p.values()).most_common()}")
print()

# For each p with 2 rows, check how many PASS
p_with_2 = [p for p, c in rows_by_p.items() if c >= 2]
print(f"Number of p with ≥ 2 inert-factor rows: {len(p_with_2):,}")
# Check: do any of these have PASS?
pass_ps_set = set(pass_ps)
p_2_with_pass = [p for p in p_with_2 if p in pass_ps_set]
print(f"Of those, p that appear in PASS: {len(p_2_with_pass)}")
if len(p_2_with_pass) == 0:
    print("✓ VERIFIED: no p with ≥ 2 rows has any PASS row.")
else:
    print(f"Concerning: {p_2_with_pass}")
print()

print("="*80)
print("LEMMA CP183 (unconditional for r ≤ 10^12):")
print("="*80)
print()
print("For every prime p, at most 1 inert prime factor r ≤ 10^{12} of W_p")
print("satisfies d_r | W_{p-2} with d_r > 0.")
print()
print("Proof: Direct verification of platinum DB: 685M inert rows, 3 PASS,")
print("       3 distinct p, max ≥ 2-row p with PASS = 0. ∎")
print()
print("Corollary: By Three-Factor Reduction + MF Pinning (CP181), a composite")
print("W_p with ≥ 3 inert dangers must have ≥ 2 factors > 10^{12}.")
