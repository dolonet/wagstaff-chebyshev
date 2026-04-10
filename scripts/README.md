# Scripts

All scripts are pure Python 3 and depend only on the standard library,
except for `primality_bls.py` which uses `sympy` (and optionally
`gmpy2` / Pari/GP).

| Script | What it does | Paper section | Runtime | Peak RAM |
|---|---|---|---|---|
| `survey.py` | Runs the inert-factor survey end to end | §9 | ~4 h @128 cores | ~64 GB |
| `rerun_segments.py` | Re-runs an explicit list of survey segments (fsync-strict) | §9 | tens of seconds / segment | ~1 GB |
| `build_clean.py` | Dedupes + sorts + validates the raw survey CSV | §9 | ~2 min | <1 GB |
| `audit.py` | Streams the clean CSV and checks it against `segment_stats.jsonl` | §9 | ~30 s | ~16 MB |
| `verify_sample.py` | Independently recomputes every column on a CSV slice | §9 | ~2 min / 1000 rows | ~50 MB |
| `primality_bls.py` | BLS N-1 primality proof that W_12391 is prime | §3, §9 | ~70 s | ~500 MB |

## Entry points

```sh
# Full survey (128-core machine, ~4 hours wall-clock)
python3 survey.py --cores 128 --max-r 1e12 --output-dir out/

# Re-run specific segments (deterministic, fsync-strict)
python3 rerun_segments.py --segments 4001,4002,4003 \
    --max-r 1e12 --cores 32 --output-dir rerun/

# Canonicalise raw survey output
python3 build_clean.py out/inert_factors.csv out/segment_stats.jsonl \
    --output-dir out/clean/

# Audit a clean CSV
python3 audit.py out/clean/inert_factors.csv out/segment_stats.jsonl \
    --output audit_report.json

# Independent row-by-row verification of a sample
python3 verify_sample.py ../data/sample_1000.csv

# BLS primality proof
python3 primality_bls.py
```

## Dependencies

```sh
pip install sympy        # required for primality_bls.py
pip install gmpy2        # optional; 10× faster BLS proof
```

Pari/GP (`gp`) is optional; `primality_bls.py` uses it for APR-CL
verification of factors larger than 82 digits. Without Pari/GP the
script falls back to BPSW primality testing.

See `../reproducibility.md` for the end-to-end walkthrough.
