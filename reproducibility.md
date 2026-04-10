# Reproducibility walkthrough

This document lists every step required to reproduce the paper's
computational claims, from scratch.

## Prerequisites

- **Python** ≥ 3.10 (standard library only for the survey pipeline).
- **sympy** and optionally **gmpy2** for the BLS primality driver.
- **Pari/GP** (`gp`) for APR-CL verification of large cyclotomic
  factors in the BLS driver. Optional: if absent, the BLS script falls
  back to BPSW for every factor.
- **TeX Live** (`pdflatex`) to rebuild the paper PDF.
- **~16 MB RAM** for the audit script; **~128 cores** recommended for a
  full-scale survey re-run.

Python dependencies:

```sh
pip install sympy gmpy2
```

## Step 1 — download the full survey CSV

The canonical clean CSV (802 MB, 15 587 021 rows) is hosted on Zenodo.
Download `inert_factors.csv` and `segment_stats.jsonl` into a working
directory of your choice.

Verify the download:

```sh
cd /path/to/downloads
shasum -a 256 -c /path/to/repo/data/SHA256SUMS
```

You should see:

```
inert_factors.csv: OK
```

(The file `sample_1000.csv` is shipped inside this repo at
`data/sample_1000.csv` and its hash is included in `SHA256SUMS`.)

## Step 2 — audit the CSV

```sh
python3 scripts/audit.py inert_factors.csv segment_stats.jsonl \
    --output audit_report.json
```

Expected output (abridged):

```
Segments in stats: 5000
Expected CSV rows: 15,587,021
Total CSV data rows: 15,587,021
Bad-parse rows:      0
Net diff:            +0
Bad segments: 0
Total duplicates within segments: 0
```

A non-zero `net_diff`, any `Bad segments`, or any duplicates means the
download is corrupt. The full audit runs in ≈30 seconds and peaks at
~16 MB of RAM.

## Step 3 — independent recomputation of a sample slice

```sh
python3 scripts/verify_sample.py data/sample_1000.csv
```

`verify_sample.py` recomputes `p`, `r`, `G_r/4`, and (where applicable)
the full order `ord_r(ω₃)/4` and the obstruction column for every row.
Expected:

```
Total rows: 1000
Passed:     1000
Failed:     0

OK — sample CSV verified.
```

You can run the same script on the full 802 MB CSV as well — expect it
to take several hours because it does not exploit any shortcut.

## Step 4 — re-run survey segments (optional)

To re-run a specific survey segment and compare row-by-row against the
published CSV:

```sh
python3 scripts/rerun_segments.py \
    --segments 1234 \
    --segment-size 100000000 \
    --max-r 1e12 \
    --cores 32 \
    --output-dir rerun_1234/
```

Each segment covers one billion prime-exponent candidates; a single
segment takes tens of seconds to a few minutes on a 32-core machine.
Then diff the re-run CSV against the corresponding slice of the
published CSV:

```sh
awk -F, -v lo=123400000000 -v hi=123500000000 \
    'NR>1 && $1>=lo && $1<hi' inert_factors.csv \
    | sort > slice_pub.csv
tail -n +2 rerun_1234/inert_factors.csv | sort > slice_rerun.csv
diff slice_pub.csv slice_rerun.csv  # should be empty
```

The arithmetic is deterministic and the row order is sorted after
`build_clean.py`, so a re-run should produce a bit-identical slice.

## Step 5 — full survey re-run (expensive)

A full `p < 5×10¹¹`, `r < 10¹²` survey takes ≈4 hours on a 128-core
workstation with 64 GB RAM. Segments are checkpointed, so the run can
be resumed across machine restarts.

```sh
python3 scripts/survey.py \
    --cores 128 \
    --max-r 1e12 \
    --segment-size 100000000 \
    --output-dir full_rerun/
```

Output lands in `full_rerun/`:

```
inert_factors.csv    — raw survey output (~840 MB)
segment_stats.jsonl  — per-segment totals
checkpoint.json      — resume state
params.json          — run parameters
summary.json         — final totals
survey.log           — timestamped log
```

Then canonicalise (dedupes, sorts, validates against `segment_stats.jsonl`):

```sh
python3 scripts/build_clean.py \
    full_rerun/inert_factors.csv \
    full_rerun/segment_stats.jsonl \
    --output-dir full_rerun/clean/
```

The resulting `full_rerun/clean/inert_factors.csv` should be
byte-identical to the Zenodo download (the SHA-256 hash is pinned in
`data/SHA256SUMS`).

## Step 6 — BLS primality proofs

The paper proves three Wagstaff primes via BLS N-1 (§3.6--3.8).
Dedicated scripts reproduce each proof and write JSON certificates
to `data/`:

```sh
python3 scripts/bls_n_minus_1_w2617.py     # W_2617,  §3.6, ~5 s
python3 scripts/bls_n_minus_1_w10501.py    # W_10501, §3.7, ~5 min
python3 scripts/bls_n_minus_1_w12391.py    # W_12391, §3.8, ~70 s (with gmpy2)
```

Pre-computed certificates are shipped in `data/bls_certificate_w*.json`.

The original generic driver `scripts/primality_bls.py` also reproduces
the W_12391 proof:

```sh
python3 scripts/primality_bls.py
```

## Step 7 — NCT and Pell exclusion verification

Verify the No Compatible Triples theorem (d ≤ 81, §5.6) and the Pell
exclusion result (d ≤ 43, §6.3):

```sh
python3 scripts/nct_verify.py
```

Runtime: ≈1 minute. Requires `sympy`. Expected output ends with:

```
ALL CHECKS PASSED. No compatible triple exists for d <= 81.
Pell exclusion verified for all admissible d <= 43.
```

## Step 8 — split factor census and R3 verification

Reproduce the 185-factor split census and verify Conjecture R3 at
every known split factor (§8.2):

```sh
python3 scripts/split_factor_census.py
```

Runtime: ~5 minutes. Requires `sympy`. Expected output ends with:

```
R3 verified: 185/185
```

The case breakdown should be: Case A: 1, Case B: 6, Case 4U: 46,
Case C: 132 (43 d-inadmissible + 89 general d>200).

## Step 9 — primitive divisor survey

Reproduce the primitive divisor survey for small Class III primes (§8.3):

```sh
python3 scripts/primitive_divisor_survey.py 200
```

Runtime: ~3 seconds. Expected: 10 primitive divisors among q <= 200,
exactly one qth-power case (q=43, r=14449).

## Step 10 — NCT extended verification

Search for compatible triples at parity-unblocked d in [83, 200] (§8.4):

```sh
python3 scripts/nct_extended_verify.py 1e8
```

Runtime: ~6 seconds at r < 10^8. Expected: 0 compatible triples.
The paper's full computation used r < 10^11; this script scales to that
range with more time.

## Step 11 — auxiliary verifications

These scripts verify specific claims from individual paper sections:

```sh
# BLS ceiling analysis (§3.5): which W_p admit BLS proofs
python3 scripts/bls_ceiling.py

# Verify Theorem 3.2 for all 20 known Wagstaff primes (§3.4)
python3 scripts/wagstaff_verify_all.py

# Parity obstruction analysis (§5.5): which d are blocked
python3 scripts/nct_parity_obstruction.py

# Class III Wieferich check (§6.2): v_q(2^m+1)=1 for all q<5000
python3 scripts/class_iii_wieferich.py

# Secondary-factor closures for d=57, d=67 (§8.8)
python3 scripts/secondary_closure.py

# 869-factor inert census (§8.5)
python3 scripts/small_factor_census.py
```

## Step 12 — Chebyshev test smoke-check

Run the standalone Chebyshev probable-prime test (§8.6) on the three
BLS-proved exponents as a sanity check:

```sh
python3 scripts/chebyshev_test.py 2617
python3 scripts/chebyshev_test.py 10501
python3 scripts/chebyshev_test.py 12391
```

Each should print `PROBABLE PRIME (Chebyshev test passed)`.

## Step 13 — rebuild the paper PDF

```sh
cd paper/
make
```

This invokes `pdflatex` three times (for internal references) and
produces `wagstaff_chebyshev.pdf`.

## Expected hardware

| Step | CPU | RAM | Disk |
|---|---|---|---|
| Audit (step 2) | 1 core | 16 MB | 802 MB (the CSV) |
| Verify sample (step 3) | 1 core | 50 MB | — |
| Re-run one segment (step 4) | 32 cores | 1 GB | tens of MB |
| Full survey (step 5) | 128 cores | 64 GB | 840 MB |
| BLS proofs (step 6) | 1 core | 500 MB | — |
| NCT verification (step 7) | 1 core | 100 MB | — |
| Split census (step 8) | 1 core | 200 MB | — |
| Primitive divisors (step 9) | 1 core | 100 MB | — |
| NCT extended (step 10) | 1 core | 100 MB | — |
| Auxiliary (step 11) | 1 core | 200 MB | — |
| Chebyshev test (step 12) | 1 core | 50 MB | — |
| Paper build (step 13) | 1 core | — | — |

## Provenance

All scripts in `scripts/` are pure Python and self-contained. The
survey pipeline has no network, subprocess, or C extension
dependencies. The only script that calls an external binary is
`primality_bls.py`, which optionally invokes Pari/GP's `gp` for APR-CL
verification of large cyclotomic factors.
