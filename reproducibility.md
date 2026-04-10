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

## Step 6 — BLS primality proof

To re-verify that `W_12391` is prime by BLS N-1 (the current ceiling of
this approach):

```sh
python3 scripts/primality_bls.py
```

Runtime: ≈70 seconds on a modern CPU with `gmpy2`, several minutes
without. The script writes a certificate file
`bls_certificate_w12391.json` and a detailed log alongside itself.

## Step 7 — rebuild the paper PDF

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
| BLS proof (step 6) | 1 core | 500 MB | — |
| Paper build (step 7) | 1 core | — | — |

## Provenance

All scripts in `scripts/` are pure Python and self-contained. The
survey pipeline has no network, subprocess, or C extension
dependencies. The only script that calls an external binary is
`primality_bls.py`, which optionally invokes Pari/GP's `gp` for APR-CL
verification of large cyclotomic factors.
