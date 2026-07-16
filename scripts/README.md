# Scripts

Reproducibility scripts for the paper (v4.0). Pure Python 3 + `sympy`, except
the APR-CL calls (PARI/GP `gp`) and the FactorDB lookups noted below. See the
repo-root `reproducibility.md` for the end-to-end procedure and `../README.md`
for the high-level map. Scripts named `cpNNN_*` are pinned as they were run
in the research tree; a few of the census/verification drivers read
research-tree paths and are shipped as provenance for the corresponding
static artifacts in `data/` (noted per table below) — everything else runs
in-place from this directory.

## NCT fixed-d certificates (cor:nct-fixed-d-certificate)

| Script | What it does |
|---|---|
| `cp365_nct_certificate_bundle.py` | Builds the static `data/nct_certificates.json`: pinned `Psi_{4d}` factorizations + APR-CL + dispositions for all 28 closures in (100,400] (with the `d <= 99` theorem: NCT for every odd `d <= 400`). FactorDB-independent once generated. Needs `gp`. |
| `cp350_nct_107_121_verify.py` | the first two fixed-d closures, d = 107, 121 |
| `cp351_nct_certificate_129_139_171_177.py` | per-d closures in (100,200] |
| `cp353_nct_certificate_{131,163,179}.py` | per-d closures (local ECM/GNFS) |
| `cp356_nct_certificates_200_300.py` | (200,300] closures (FactorDB-FF) |
| `cp358_nct_certificate_{211,243,363}.py` | (200,400] closures (local ECM) |
| `cp358_nct_certificates_range.py` | parameterized (300,400] closures |
| `cp367_nct_certificate_227.py`, `cp375_nct_certificate_283.py`, `cp376_nct_certificate_331.py`, `cp372_nct_certificate_379.py` | the last four (200,400] closures (Pell half-split `V_{2d} = 2(2P_d-1)(2P_d+1)` + local ECM/GNFS) |
| `cp352_verify_nct_recompute.py` | independent end-to-end re-check of the NCT closures |
| `cp356_aprcl_recert.py`, `cp358_aprcl_range.py` | APR-CL re-certification (PARI flag 2) |

## Danger census through d = 400 (prop:danger-census-400, cor:inert-floor-400)

| Script | What it does |
|---|---|
| `cp412_danger_census_d400.py` | the census driver: 34 admissible rungs `43 < d <= 400`, complete `V_{2d}` factorizations, tests (i)–(v) per primitive prime factor; writes `data/cp412_danger_census.json`. Provenance script (reads research-tree sources); the shipped JSON is the self-contained artifact. |
| `cp412b_t1_recheck.py` | adversarial re-implementation reproducing all 340 verdicts (provenance) |
| `cp412b_t1_aprcl_sweep.py` | APR-CL sweep of all 226 distinct primes of the census (provenance) |

## Witness discharges and X1 (sec:x1)

| Script | What it does |
|---|---|
| `cp352_verify_discharges_recompute.py` | self-contained end-to-end recomputation: the two discharged witness exponents carry an explicit failing factor (parity obstruction, direct congruence evaluation bypassing the descent); the three Platinum witnesses are genuine local passes. Runs in-place, no network. |

## Diagonal, corner and falsification censuses

| Script | What it does |
|---|---|
| `cp354_samed_diagonal_certificates.py` | same-d diagonal certificates (danger-grade divisors grouped by pinned exponent) |
| `cp354_factordb_v2d_diagonal.py` | `V_{2d}` factorization harvest for the diagonal closures |
| `cp354_diag_corner_gcd.py` | diagonal corner census (d <= 124) |
| `cp358_corner_sweep_d1000.py` | corner sweep d <= 1000, p <= 1e8 |
| `cp358_falsification_prewindow.py` | falsification pre-window scan (k <= 1e6 deep-k complement) |

## Heuristics and the ladder

| Script | What it does |
|---|---|
| `cp395_pppc_heuristic_constants.py` | the M1–M3 model constants and the bracketed `E[Pi]` figures (heuristic only; proves nothing). Runs in-place. |
| `cp410_twin_ladder.py` | the self-supporting exponent ladder (`W_{p-2}` prime rungs). Runs in-place. |

## Cross-case reductions (pinning, Platinum, exact-AP)

| Script | What it does |
|---|---|
| `platinum_lemma.py` | Platinum Lemma: 684,965,381-row enumeration (needs the Zenodo CSV) |
| `multi_factor_pinning.py` | Order-Pinning + Multi-Factor Pinning + phantom exclusion |
| `exact_ap_density.py` | exact-AP characterization + pair-count |
| `second_moment_reduction.py` | **historical**: the `E_3 <= Pi + W` counting bookkeeping of pre-v4.0 versions. v4.0 routes the chain through the branch decomposition instead; every check this script performs remains true, and it is kept as a cross-check on the counting remark (`N_in <= Pi + 1`). |
| `nct_*.py`, `primitive_divisor_survey.py`, `small_factor_census.py` | NCT / inert censuses |

## Vrba–Reix agreement & independent verification

| Script | What it does |
|---|---|
| `cp361_vrba_reix_check.py` | global agreement of the `S^2-2` return tests with Condition (II) on tested exponents (empirical agreement of global verdicts, **not** per-factor equivalence — see the Erratum in `../README.md`) |
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
consolidated `nct_certificates.json` and `cp412_danger_census.json` pin those
results so the closures and the census are verifiable offline.
