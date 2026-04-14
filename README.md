# Chebyshev primality criteria for Wagstaff numbers

This repository contains the LaTeX source, reproducibility scripts, and
data pointers for the paper **"Chebyshev primality criteria for
Wagstaff numbers"** by Alexey Dolotov.

## Abstract

The Chebyshev base `omega_3 = 3 + 2*sqrt(2)` and the identity
`W_p + 1 = 4*W_{p-2}` yield an unconditional `N+1` primality criterion
for the Wagstaff numbers `W_p = (2^p + 1)/3`; combined with a companion
`N-1` criterion and the Brillhart-Lehmer-Selfridge framework, this gives
fully verified primality proofs of `W_2617`, `W_10501`, and `W_12391`,
independent of ECPP. For the sufficiency conjecture -- whether the
single congruence `omega_3^{(W_p+1)/2} = -1 (mod W_p)` implies
primality -- we show that Condition~(I) holds automatically, exclude
factors `r = 5, 7 (mod 8)`, reduce the split branch `r = 1 (mod 8)`
to the nonexistence of compatible Pell triples (proved for rank
parameters `d <= 81`), and prove that the inert branch `r = 3 (mod 8)`
is impossible for order parameter `d <= 43`. An Order-Pinning Lemma
shows that every inert prime factor `r` of a composite `W_p` satisfies
`ord_r(2) = 2p`, and a Multi-Factor Pinning Theorem forces all inert
factors of a common `W_p` to share a single order parameter; combined
with a multi-factor bound, these unconditionally exclude the entire
finite catalog of small-`d` danger configurations produced by the
parity obstruction. An exhaustive Platinum Lemma rules out every inert
factor with `r <= 10^12` via a 684,965,381-row enumeration, and an
elementary deterministic Second-Moment Reduction transfers the
residual problem to the single claim that no composite `W_p` admits a
pair of distinct residual inert danger factors `r_1, r_2 > 10^12` at
the same `p`.

## Layout

```
.
├── paper/                         LaTeX source of the paper
│   ├── wagstaff_chebyshev.tex
│   └── Makefile                   pdflatex targets
├── scripts/                       Reproducibility scripts (27 Python scripts)
│   ├── README.md                  script -> paper section map
│   │
│   │  Survey pipeline (§9.5)
│   ├── survey.py                  inert-factor survey (128 cores, ~4 h)
│   ├── rerun_segments.py          re-run explicit survey segments
│   ├── build_clean.py             canonicalise the raw survey CSV
│   ├── audit.py                   CSV integrity audit
│   ├── verify_sample.py           independent verifier for a CSV slice
│   │
│   │  BLS primality proofs (§3.6-3.8)
│   ├── bls_n_minus_1_w2617.py     BLS proof for W_2617
│   ├── bls_n_minus_1_w10501.py    BLS proof for W_10501
│   ├── bls_n_minus_1_w12391.py    BLS proof for W_12391
│   ├── bls_ceiling.py             BLS N+1 ceiling analysis
│   ├── primality_bls.py           BLS N-1 proof driver (generic)
│   ├── wagstaff_verify_all.py     verify Theorem 3.2 for all known W_p
│   ├── chebyshev_test.py          standalone Chebyshev primality test
│   │
│   │  Split factor analysis (§5, §9.2)
│   ├── split_factor_census.py     185-factor census and case breakdown
│   ├── r3_split_sweep.py          extended R3 sweep at r <= 10^11
│   │
│   │  NCT and Pell exclusion (§5-6, §9.3-9.4)
│   ├── nct_verify.py              NCT (d<=81) and Pell exclusion (d<=43)
│   ├── nct_parity_obstruction.py  parity obstruction analysis
│   ├── nct_extended_verify.py     NCT search for d in [83, 200]
│   ├── primitive_divisor_survey.py primitive divisor survey (§9.3)
│   │
│   │  Inert factor analysis (§6-7, §9.5, §9.8)
│   ├── small_factor_census.py     869-factor census (r < 2x10^5)
│   ├── class_iii_wieferich.py     Class III Wieferich check (§6.2)
│   ├── secondary_closure.py       d=57, d=67 closure (§9.8)
│   ├── danger_triple_survey.py    danger-triple survey for admissible d (§6.7)
│   │
│   │  Unconditional reductions from pinning (§8)
│   ├── multi_factor_pinning.py    Order-Pinning + Multi-Factor Pinning + Phantom Exclusion (§8.1-8.3)
│   ├── platinum_lemma.py          Platinum Lemma direct verification (§8.4)
│   ├── exact_ap_density.py        Exact AP density (§8.5) and rigorous pair-count sharpening
│   └── second_moment_reduction.py Deterministic Second-Moment Reduction bookkeeping (§8.6)
│
├── data/                          Sample data, hashes, and BLS certificates
│   ├── README.md                  data dictionary + schema
│   ├── sample_1000.csv            first 1000 rows of the clean CSV
│   ├── SHA256SUMS                 hashes of full and sample CSVs
│   ├── bls_certificate_w2617.json
│   ├── bls_certificate_w10501.json
│   ├── bls_certificate_w12391.json
│   └── danger_triple_data.json    V_{114}, V_{134}, V_{662} factorisations + 5 danger triples
├── reproducibility.md             end-to-end reproduction walkthrough
├── LICENSE                        MIT
├── CITATION.cff                   machine-readable citation
└── README.md                      this file
```

## Cite this work

If you use the survey data or the computational scripts, please cite:

- **Paper (arXiv):** *pending submission -- arXiv ID will be inserted here*
- **Survey CSV (Zenodo):** [10.5281/zenodo.19496206](https://doi.org/10.5281/zenodo.19496206)

A machine-readable citation is provided in `CITATION.cff`.

## Reproduce the central computational claims in 4 steps

1. **Download the full CSV** from Zenodo (802 MB, 15 587 021 rows).
2. **Verify integrity:** `shasum -a 256 -c data/SHA256SUMS` (expected
   hash is pinned in `data/README.md`).
3. **Run the audit:** `python3 scripts/audit.py inert_factors.csv
   segment_stats.jsonl`. Expect `Bad segments: 0` and
   `Total duplicates: 0`.
4. **Spot-check a slice:** `python3 scripts/verify_sample.py
   data/sample_1000.csv`. Expect `Passed: 1000 / Failed: 0`.

For the full reproduction procedure, including re-running the survey
and the four pinning-reduction verifications (§8.1-8.6),
see `reproducibility.md`.

## Contact

Issues and questions: please open a GitHub issue on this repository.

## License

MIT. See `LICENSE`.
