# Scripts

Reproducibility scripts for the paper (v3.6). Pure Python 3 + `sympy`, except
the APR-CL calls (PARI/GP `gp`) and the FactorDB lookups noted below. See the
repo-root `reproducibility.md` for the end-to-end procedure and `../README.md`
for the high-level map.

## NCT fixed-d certificates (cor:nct-fixed-d-certificate)

| Script | What it does |
|---|---|
| `cp365_nct_certificate_bundle.py` | Builds the static `data/nct_certificates.json`: pinned `Psi_{4d}` factorizations + APR-CL + dispositions for all 24 closures in (100,400]. FactorDB-independent once generated. Needs `gp`. |
| `cp351_nct_certificate_129_139_171_177.py` | per-d closures in (100,200] |
| `cp353_nct_certificate_{131,163,179}.py` | per-d closures (local ECM/GNFS) |
| `cp356_nct_certificates_200_300.py` | (200,300] closures (FactorDB-FF) |
| `cp358_nct_certificate_{211,243,363}.py` | (200,400] closures (local ECM) |
| `cp358_nct_certificates_range.py` | parameterized (300,400] closures |
| `cp356_aprcl_recert.py`, `cp358_aprcl_range.py` | APR-CL re-certification (PARI flag 2) |

## Inert branch, pinning, exact-AP (Section: cross-case reductions)

| Script | What it does |
|---|---|
| `platinum_lemma.py` | Platinum Lemma: 684,965,381-row enumeration (needs the Zenodo CSV) |
| `second_moment_reduction.py` | `E_3 <= Pi + W` bookkeeping |
| `multi_factor_pinning.py` | Order-Pinning + Multi-Factor Pinning + phantom exclusion |
| `exact_ap_density.py` | exact-AP characterization + pair-count |
| `nct_*.py`, `primitive_divisor_survey.py`, `small_factor_census.py` | NCT / inert censuses |

## Corner / falsification censuses

| Script | What it does |
|---|---|
| `cp354_diag_corner_gcd.py` | diagonal corner census (d <= 124) |
| `cp358_corner_sweep_d1000.py` | corner sweep d <= 1000, p <= 1e8 |
| `cp358_falsification_prewindow.py` | falsification pre-window scan |

## Vrba–Reix equivalence & independent verification

| Script | What it does |
|---|---|
| `cp361_vrba_reix_check.py` | confirms the recurrence test == Condition (II), tracks primality |
| `cp362_verify_preliminaries.py` | re-derives mod-8/16, ord2, Condition (I), order √2 |
| `cp363_verify_inert_foundations.py` | re-derives LTE, Pell-exclusion, parity; reproduces 845/869 census |

## Survey pipeline (the Platinum enumeration data)

| Script | What it does | Runtime |
|---|---|---|
| `survey.py` | full inert-factor survey | ~4 h @128 cores |
| `rerun_segments.py` | re-run explicit segments (deterministic) | tens of s / segment |
| `build_clean.py` | dedupe + sort + validate the raw CSV | ~2 min |
| `audit.py` | stream-audit the clean CSV vs segment stats | ~30 s |
| `verify_sample.py` | independently recompute every column on a slice | ~2 min / 1000 rows |

Dependencies: `pip install sympy`; PARI/GP `gp` for APR-CL. The per-d
certificate scripts query FactorDB for order-descent factorizations; the
consolidated `nct_certificates.json` pins those results so the closures are
verifiable offline.
