# Chebyshev primality criteria for Wagstaff numbers

This repository contains the LaTeX source, reproducibility scripts, and
data pointers for the paper **"Chebyshev primality criteria for
Wagstaff numbers"** by Alexey Dolotov.

## Abstract

The Chebyshev base `omega_3 = 3 + 2*sqrt(2)` and the identity
`W_p + 1 = 4*W_{p-2}` yield an unconditional `N+1` primality criterion
for the Wagstaff numbers `W_p = (2^p + 1)/3`; combined with a companion
`N-1` criterion and the Brillhart-Lehmer-Selfridge framework, this gives
self-contained primality proofs of `W_2617`, `W_10501`, and `W_12391`,
independent of the ECPP proofs in the literature. For the sufficiency of
the single Chebyshev congruence, we reduce the split branch
`r = 1 (mod 8)` to the nonexistence of compatible Pell triples (NCT),
proved for all rank parameters `d <= 81`; Conjecture R3 is verified
computationally at all 185 known split factors. For the inert branch we
prove that `r = 3 (mod 8)` cannot pass Condition (II) whenever
`d <= 43`.

## Layout

```
.
├── paper/                         LaTeX source of the paper
│   ├── wagstaff_chebyshev.tex
│   └── Makefile                   pdflatex targets
├── scripts/                       Reproducibility scripts (22 Python scripts)
│   ├── README.md                  script -> paper section map
│   │
│   │  Survey pipeline (§8.5)
│   ├── survey.py                  inert-factor survey (128 cores, ~4 h)
│   ├── rerun_segments.py          re-run explicit survey segments
│   ├── build_clean.py             canonicalise the raw survey CSV
│   ├── audit.py                   CSV integrity audit
│   ├── verify_sample.py           independent verifier for a CSV slice
│   │
│   │  BLS primality proofs (§3.5-3.8)
│   ├── bls_n_minus_1_w2617.py     BLS proof for W_2617
│   ├── bls_n_minus_1_w10501.py    BLS proof for W_10501
│   ├── bls_n_minus_1_w12391.py    BLS proof for W_12391
│   ├── bls_ceiling.py             BLS N+1 ceiling analysis
│   ├── primality_bls.py           BLS N-1 proof driver (generic)
│   ├── wagstaff_verify_all.py     verify Theorem 3.2 for all known W_p
│   ├── chebyshev_test.py          standalone Chebyshev primality test
│   │
│   │  Split factor analysis (§5, §8.2)
│   ├── split_factor_census.py     185-factor census and case breakdown
│   ├── r3_split_sweep.py          extended R3 sweep at r <= 10^11
│   │
│   │  NCT and Pell exclusion (§5-6, §8.3-8.4)
│   ├── nct_verify.py              NCT (d<=81) and Pell exclusion (d<=43)
│   ├── nct_parity_obstruction.py  parity obstruction analysis
│   ├── nct_extended_verify.py     NCT search for d in [83, 200]
│   ├── primitive_divisor_survey.py primitive divisor survey (§8.3)
│   │
│   │  Inert factor analysis (§6-7, §8.5, §8.8)
│   ├── small_factor_census.py     869-factor census (r < 2x10^5)
│   ├── class_iii_wieferich.py     Class III Wieferich check (§6.2)
│   └── secondary_closure.py       d=57, d=67 closure (§8.8)
│
├── data/                          Sample data, hashes, and BLS certificates
│   ├── README.md                  data dictionary + schema
│   ├── sample_1000.csv            first 1000 rows of the clean CSV
│   ├── SHA256SUMS                 hashes of full and sample CSVs
│   ├── bls_certificate_w2617.json
│   ├── bls_certificate_w10501.json
│   └── bls_certificate_w12391.json
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

For the full reproduction procedure, including re-running the survey,
see `reproducibility.md`.

## Contact

Issues and questions: please open a GitHub issue on this repository.

## License

MIT. See `LICENSE`.
