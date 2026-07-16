# Reproducibility walkthrough (v4.0)

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

## Step 1 — NCT fixed-d certificates (pinned; offline-verifiable)

The central unconditional split-branch claim — NCT closed for **every odd
`d <= 400`** (the `d <= 99` theorem plus 28 fixed-`d` closures in
`(100, 400]`) — is verifiable from the shipped static artifact
`data/nct_certificates.json`. For each closed `d`
it pins the complete factorization of the Pell primitive part `Psi_{4d}` and,
per prime factor, its APR-CL primality and the certificate disposition (mod-8
class; for split `r = 1 mod 8`: `ord_r(2)` parity → half-order primality →
condition (v) `d | W_{p-2}`). NCT is closed at `d` iff no factor reaches
condition (v) as a survivor.

The shipped JSON is the verification object: every factorization carries an
exact product identity (`prod(factors) == Psi_{4d}`, with `Psi_{4d}`
recomputable from scratch by the recurrence in the generator script), and
every listed factor can be re-certified offline with `gp`'s `isprime(., 2)`.
To re-**generate** the artifact from scratch:

```sh
python3 scripts/cp365_nct_certificate_bundle.py
```

Expected: `DONE: 28/28 CLOSED`. Note this **overwrites**
`data/nct_certificates.json` in place, requires `gp`, and — for the values
without embedded factorizations — FactorDB network access. The
locally-factored values (e.g.
`d in {131,163,179,211,227,243,283,331,363,379}`) carry embedded
factorizations; the rest are pulled once from the public factor
tables and pinned into the JSON. Every factor is APR-CL-certified
(182 factors total, all prime); if `gp` is absent the generator falls back
to BPSW, and the output is then only BPSW-checked, **not** certified. The 28
values are exactly the admissible odd `d` in `(100, 400]` — for every other
odd `d` there, no compatible triple can exist for congruence/parity reasons
and no certificate is needed. Individual per-d certificates are reproduced
by `scripts/cp350_nct_107_121_verify.py` and the
`cp351`/`cp353`/`cp356`/`cp358`/`cp367`/`cp372`/`cp375`/`cp376` per-`d`
certificate scripts,
an independent end-to-end re-check by `scripts/cp352_verify_nct_recompute.py`,
and the APR-CL re-cert by `scripts/cp356_aprcl_recert.py`,
`scripts/cp358_aprcl_range.py`.

## Step 1b — danger census through d = 400 (static artifact)

The complete danger-triple census for admissible `43 < d <= 400`
(Proposition `danger-census-400` / Corollary `inert-floor-400`) is shipped as
`data/cp412_danger_census.json`: 34 rungs, complete `V_{2d}` factorizations,
340 (rung, prime) tests, 226 distinct primes, exactly two real triples (both
at exponents closed by secondary factors) plus two phantoms at prime `W_p`.
The generating and adversarial-recheck drivers
(`scripts/cp412_danger_census_d400.py`, `scripts/cp412b_t1_recheck.py`,
`scripts/cp412b_t1_aprcl_sweep.py`) are pinned provenance — they were run in
the research tree and read its layout; the JSON itself is self-contained
(every factorization carries an exact product identity, and each listed prime
re-checks with `gp`'s `isprime(., 2)`).

## Step 1c — witness discharges (self-contained)

```sh
python3 scripts/cp352_verify_discharges_recompute.py
```

Re-verifies end-to-end, with no network and by direct congruence evaluation
(bypassing the order descent), that the two discharged witness exponents
carry an explicit prime factor failing Condition (II), and that the three
Platinum witnesses are genuine local passes.

## Step 2 — Vrba–Reix global agreement (empirical)

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
CSV (802 MB, 15,587,021 rows) is on Zenodo `10.5281/zenodo.19496206`. The
clean CSV retains exactly the 15,587,021 factors with `G_r/4 > 43` — by the
Pell-sequence exclusion theorem no `d <= 43` factor can pass, so these are
the only factors at which a local pass is arithmetically possible; the
remaining 669,378,360 enumerated factors are excluded wholesale by that
theorem and tallied per-segment in `segment_stats.jsonl` (shipped in the same
Zenodo record), against which `scripts/audit.py` reconciles the row count.

```sh
# verify the download against the pinned hash (hashes are relative to data/)
(cd data && shasum -a 256 -c SHA256SUMS)    # inert_factors.csv: OK
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

Expected: 15,587,021 CSV rows (the `G_r/4 > 43` survivors of the
684,965,381-factor enumeration, per the note above), exactly 3 PASS rows,
all at distinct exponents, zero `p` with ≥ 2 passes.

To re-run the survey from scratch (≈4 h on 128 cores) use `scripts/survey.py`
then `scripts/build_clean.py`; a single segment is reproduced by
`scripts/rerun_segments.py` (deterministic — a re-run yields a bit-identical
slice).

## Step 5 — cross-case reductions and censuses

```sh
python3 scripts/multi_factor_pinning.py        # Order/Multi-Factor Pinning, phantom exclusion
python3 scripts/exact_ap_density.py            # exact-AP characterization + pair-count
python3 scripts/cp354_diag_corner_gcd.py       # diagonal corner census (d <= 124)
python3 scripts/cp358_corner_sweep_d1000.py    # corner sweep d <= 1000, p <= 1e8 (long: hours)
python3 scripts/cp358_falsification_prewindow.py  # falsification pre-window scan (long: hours)
python3 scripts/cp395_pppc_heuristic_constants.py # M1-M3 heuristic constants (heuristic only)
python3 scripts/cp410_twin_ladder.py           # self-supporting exponent ladder
python3 scripts/second_moment_reduction.py     # historical: pre-v4.0 counting bookkeeping
```

(`second_moment_reduction.py` verifies the `E_3 <= Pi + W` bookkeeping of
pre-v4.0 versions; the v4.0 chain routes through the branch decomposition,
and the script is kept as a cross-check on the counting remark.)

## Step 6 — rebuild the paper PDF

```sh
cd paper/ && make        # pdflatex x3 -> wagstaff_chebyshev_reduction_v4.0.pdf
```

(The last pre-restructure version, v3.8, is archived in `paper/archive/`;
that directory has no Makefile — build it with a direct
`pdflatex -interaction=nonstopmode wagstaff_chebyshev_reduction_v3.8.tex`
run, three times.)

## Expected hardware

| Step | CPU | Notes |
|---|---|---|
| NCT certificates (1) | 1 core + gp | minutes (APR-CL on 182 factors) |
| Census artifact check / discharges (1b, 1c) | 1 core (+ gp for APR-CL) | seconds–minutes |
| Vrba–Reix / verification (2,3) | 1 core | seconds–minutes |
| CSV audit / sample (4) | 1 core, ~16 MB | ~30 s / minutes |
| Platinum direct (4) | 1 core | ~25 min, needs 802 MB CSV |
| Full survey re-run (4) | 128 cores, 64 GB | ~4 h, 840 MB out |
| Reductions / censuses (5) | 1 core | seconds–minutes (corner sweep and pre-window scan: hours) |
| Paper build (6) | 1 core | seconds |

## Provenance

All scripts are pure Python + sympy and self-contained, except: the APR-CL
calls (`cp365_nct_certificate_bundle.py`, `cp356/cp358_aprcl_*`) which invoke
PARI/GP's `gp`; the FactorDB lookups in the per-`d` certificate scripts
(the consolidated `nct_certificates.json` pins those results so the closures
are verifiable without network access); and the census drivers (`cp412*`,
Step 1b), which were run in the research tree and are shipped as provenance
for the self-contained census artifact.
