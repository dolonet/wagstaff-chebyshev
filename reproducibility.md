# Reproducibility walkthrough (v3.6)

Every step required to reproduce the paper's computational claims, from
scratch. The BLS primality proofs are in the companion paper (Zenodo
`10.5281/zenodo.19645478`) and are not reproduced here.

## Prerequisites

- **Python** ≥ 3.10 with **sympy** (the survey pipeline uses the standard
  library only).
- **PARI/GP** (`gp`) for APR-CL primality certification (`isprime(n, 2)`) of
  the `Psi_{4d}` factors and the order-descent factorizations. If absent, the
  certificate generator falls back to BPSW.
- **TeX Live** (`pdflatex`) to rebuild the paper PDF.

```sh
pip install sympy
```

## Step 1 — NCT fixed-d certificates (no FactorDB, no network)

The central unconditional claim — NCT closed for every odd `d <= 200` and at
fifteen further parity-unblocked values in `(200, 400]` — is verifiable from
the shipped static artifact `data/nct_certificates.json`. For each closed `d`
it pins the complete factorization of the Pell primitive part `Psi_{4d}` and,
per prime factor, its APR-CL primality and the certificate disposition (mod-8
class; for split `r = 1 mod 8`: `ord_r(2)` parity → half-order primality →
condition (v) `d | W_{p-2}`). NCT is closed at `d` iff no factor reaches
condition (v) as a survivor.

Re-generate and re-verify it:

```sh
python3 scripts/cp365_nct_certificate_bundle.py
```

Expected: `DONE: 24/24 CLOSED`. The locally-factored values
`d in {131,163,179,211,243,363}` carry embedded factorizations checked by exact
product identity `prod(factors) == Psi_{4d}`; the rest are pulled once from the
public factor tables and pinned into the JSON. Every factor is APR-CL-certified
(154 factors total, all prime). Individual per-d certificates are reproduced by
`scripts/cp351/cp353/cp356/cp358_nct_certificate*.py` and the APR-CL re-cert by
`scripts/cp356_aprcl_recert.py`, `scripts/cp358_aprcl_range.py`.

## Step 2 — Vrba–Reix equivalence

```sh
python3 scripts/cp361_vrba_reix_check.py
```

Confirms, on the prime and composite Wagstaff numbers up to `p = 73`, that the
recurrence `S_{i+1} = S_i^2 - 2 (mod W_p)` return test (both seeds `3/2` and
`1/4`) matches Condition (II) and tracks primality of `W_p` exactly.

## Step 3 — independent verification of the foundations

The preliminaries and the inert-branch foundations are independently
re-derived and numerically checked:

```sh
python3 scripts/cp362_verify_preliminaries.py        # mod-8/16, ord2, condI, ...
python3 scripts/cp363_verify_inert_foundations.py    # LTE, Pell-exclusion, parity
```

Both print PASS on all checks (mod-8 dichotomy airtight; the 845-of-869 inert
census reproduced; LTE valuations exact).

## Step 4 — survey CSV (the Platinum enumeration)

The 684,965,381-row inert-factor enumeration (`r <= 10^12`, `p < 5×10^11`) is
the one computational input of the inert-branch reduction. The canonical clean
CSV (802 MB, 15,587,021 rows) is on Zenodo `10.5281/zenodo.19496206`.

```sh
# verify the download against the pinned hash
shasum -a 256 -c data/SHA256SUMS            # inert_factors.csv: OK
# integrity audit (~30 s, ~16 MB RAM)
python3 scripts/audit.py inert_factors.csv segment_stats.jsonl
# independent recomputation of a slice
python3 scripts/verify_sample.py data/sample_1000.csv   # Passed: 1000 / Failed: 0
```

Direct Platinum-Lemma verification over the full CSV:

```sh
# place the Zenodo CSV at data/inert_factors.csv, streaming, ~25 min
python3 scripts/platinum_lemma.py
```

Expected: 684,965,381 rows, exactly 3 local-pass rows, all at distinct
exponents, zero `p` with ≥ 2 passes.

To re-run the survey from scratch (≈4 h on 128 cores) use `scripts/survey.py`
then `scripts/build_clean.py`; a single segment is reproduced by
`scripts/rerun_segments.py` (deterministic — a re-run yields a bit-identical
slice).

## Step 5 — cross-case reductions and censuses

```sh
python3 scripts/multi_factor_pinning.py        # Order/Multi-Factor Pinning, phantom exclusion
python3 scripts/exact_ap_density.py            # exact-AP characterization + pair-count
python3 scripts/second_moment_reduction.py     # E_3 <= Pi + W bookkeeping
python3 scripts/cp354_diag_corner_gcd.py       # diagonal corner census (d <= 124)
python3 scripts/cp358_corner_sweep_d1000.py    # corner sweep d <= 1000, p <= 1e8
python3 scripts/cp358_falsification_prewindow.py  # falsification pre-window scan
```

## Step 6 — rebuild the paper PDF

```sh
cd paper/ && make        # pdflatex x3 -> wagstaff_chebyshev_reduction_v3.6.pdf
```

## Expected hardware

| Step | CPU | Notes |
|---|---|---|
| NCT certificates (1) | 1 core + gp | minutes (APR-CL on 154 factors) |
| Vrba–Reix / verification (2,3) | 1 core | seconds–minutes |
| CSV audit / sample (4) | 1 core, ~16 MB | ~30 s / minutes |
| Platinum direct (4) | 1 core | ~25 min, needs 802 MB CSV |
| Full survey re-run (4) | 128 cores, 64 GB | ~4 h, 840 MB out |
| Reductions / censuses (5) | 1 core | seconds–minutes |
| Paper build (6) | 1 core | seconds |

## Provenance

All scripts are pure Python + sympy and self-contained, except the APR-CL
calls (`cp365_nct_certificate_bundle.py`, `cp356/cp358_aprcl_*`) which invoke
PARI/GP's `gp`, and the FactorDB lookups in the per-`d` certificate scripts
(the consolidated `nct_certificates.json` pins those results so the closures
are verifiable without network access).
