# Scripts

All scripts are pure Python 3 and depend only on the standard library,
except where noted below.

## Survey pipeline

| Script | What it does | Paper section | Runtime | Peak RAM |
|---|---|---|---|---|
| `survey.py` | Runs the inert-factor survey end to end | §9.5 | ~4 h @128 cores | ~64 GB |
| `rerun_segments.py` | Re-runs an explicit list of survey segments (fsync-strict) | §9.5 | tens of seconds / segment | ~1 GB |
| `build_clean.py` | Dedupes + sorts + validates the raw survey CSV | §9.5 | ~2 min | <1 GB |
| `audit.py` | Streams the clean CSV and checks it against `segment_stats.jsonl` | §9.5 | ~30 s | ~16 MB |
| `verify_sample.py` | Independently recomputes every column on a CSV slice | §9.5 | ~2 min / 1000 rows | ~50 MB |

## Primality proofs and tests

| Script | What it does | Paper section | Runtime | Dependencies |
|---|---|---|---|---|
| `chebyshev_test.py` | Chebyshev probable-prime test for W_p | §9.6 | seconds | stdlib only |
| `primality_bls.py` | BLS N-1 primality proof driver (W_12391) | §3.8 | ~70 s | sympy, gmpy2 (opt) |
| `bls_n_minus_1_w2617.py` | BLS N-1 proof for W_2617 | §3.6 | ~5 s | sympy |
| `bls_n_minus_1_w10501.py` | BLS N-1 proof for W_10501 | §3.7 | ~5 min | sympy |
| `bls_n_minus_1_w12391.py` | BLS N-1 proof for W_12391 | §3.8 | ~70 s | sympy, gmpy2 (opt) |
| `bls_ceiling.py` | BLS N+1 ceiling analysis for all known Wagstaff primes | §3.5 | ~5 min | sympy |
| `wagstaff_verify_all.py` | Verifies Theorem 3.2 for all 20 known W_p with 5 <= p <= 1709 | §3.4 | ~5 min | sympy |

## Split factor analysis (Section 5 / §9.2)

| Script | What it does | Paper section | Runtime | Dependencies |
|---|---|---|---|---|
| `split_factor_census.py` | Census of 185 split factors, case breakdown, R3 verification | §9.2 | ~5 min | sympy |
| `r3_split_sweep.py` | Extended R3 sweep at r <= 10^11 with sum(1/p) tracking | §9.8 | hours @128 cores | stdlib only |

## NCT and Pell exclusion (Sections 5-6)

| Script | What it does | Paper section | Runtime | Dependencies |
|---|---|---|---|---|
| `nct_verify.py` | NCT verification (d <= 81) and Pell exclusion (d <= 43) | §5.7, §6.3 | ~1 min | sympy |
| `nct_parity_obstruction.py` | Parity obstruction analysis: which d are blocked | §5.5 | seconds | stdlib only |
| `nct_extended_verify.py` | NCT search for parity-unblocked d in [83, 200] | §9.4 | ~6 s @r<10^8 | sympy |
| `primitive_divisor_survey.py` | Primitive divisor survey and qth-power residue analysis | §9.3 | ~3 s @q<=200 | sympy |

## Inert factor analysis (Sections 6-7 / §9.5)

| Script | What it does | Paper section | Runtime | Dependencies |
|---|---|---|---|---|
| `small_factor_census.py` | Census of 869 inert factors with r < 2x10^5 | §9.5 | ~30 s | sympy |
| `class_iii_wieferich.py` | Verify no Class III prime q < 5000 is Wieferich-at-minus-1 | §6.2 | seconds | sympy |
| `secondary_closure.py` | Verify secondary-factor closures for d=57, d=67 | §9.8 | ~1 min | sympy |
| `danger_triple_survey.py` | Danger-triple survey for admissible d (multi-factor reduction) | §6.7 | ~10 s @d≤57 | sympy |

## Unconditional reductions from pinning (Section 8)

| Script | What it does | Paper section | Runtime | Dependencies |
|---|---|---|---|---|
| `multi_factor_pinning.py` | Order-Pinning + Multi-Factor Pinning + Phantom Exclusion checks on the full danger catalogue | §8.1-§8.3 | ~30 s | sympy |
| `platinum_lemma.py` | Direct verification of Platinum Lemma against the 684,965,381-row Zenodo CSV | §8.4 | ~25 min streaming | stdlib only |
| `exact_ap_density.py` | Exact AP-density theorem verification on class-A primes; rigorous 1/phi(ord) pair-count sharpening | §8.5 | ~1 min | sympy |
| `second_moment_reduction.py` | Deterministic Second-Moment Reduction bookkeeping; heuristic pair-count estimate | §8.6 | seconds | sympy |

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

# Chebyshev primality test for a single W_p
python3 chebyshev_test.py 2617

# NCT and Pell exclusion verification
python3 nct_verify.py

# BLS primality proofs (certificates written to ../data/)
python3 bls_n_minus_1_w2617.py
python3 bls_n_minus_1_w10501.py
python3 bls_n_minus_1_w12391.py

# BLS ceiling analysis
python3 bls_ceiling.py

# Verify Theorem 3.2 for all known Wagstaff primes
python3 wagstaff_verify_all.py

# Split factor census and case breakdown
python3 split_factor_census.py

# NCT parity obstruction analysis
python3 nct_parity_obstruction.py

# NCT extended verification (d in [83,200])
python3 nct_extended_verify.py

# Primitive divisor survey (q <= 200)
python3 primitive_divisor_survey.py

# Small inert factor census (Table 1 in §9.5)
python3 small_factor_census.py

# Class III Wieferich check
python3 class_iii_wieferich.py

# Secondary closure for d=57, d=67
python3 secondary_closure.py

# Danger-triple survey (d <= 57 feasible in ~10 s)
python3 danger_triple_survey.py 57

# Pinning reductions (Section 8 of the paper)
python3 multi_factor_pinning.py
python3 exact_ap_density.py
python3 second_moment_reduction.py
# Platinum Lemma: requires inert_factors.csv from Zenodo placed at data/
python3 platinum_lemma.py
```

## Dependencies

```sh
pip install sympy        # required for most scripts
pip install gmpy2        # optional; 10x faster BLS proof
```

Pari/GP (`gp`) is optional; `primality_bls.py` uses it for APR-CL
verification of factors larger than 82 digits. Without Pari/GP the
script falls back to BPSW primality testing.

See `../reproducibility.md` for the end-to-end walkthrough.
