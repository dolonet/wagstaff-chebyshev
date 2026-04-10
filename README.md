# Chebyshev primality criteria for Wagstaff numbers

This repository contains the LaTeX source, reproducibility scripts, and
data pointers for the paper **"Chebyshev primality criteria for
Wagstaff numbers"** by Alexey Dolotov.

## Abstract

A Wagstaff number is `W_p = (2^p + 1) / 3` with `p` an odd prime. The
paper defines a Chebyshev-style criterion (Condition (II)) that every
known Wagstaff prime satisfies, proves it sufficient for primality, and
conjectures it is necessary as well (Conjecture 5.2: no composite
Wagstaff number satisfies Condition (II)).

Conjecture 5.2 is proven unconditionally in ≈98% of cases for all
`r < 10¹²` and all composite `W_p` with `p < 5×10¹¹`. The remaining
open cases sit inside explicit parameter bands that the paper describes
precisely, and an extended computational survey of 684 965 381 inert
prime factors found only three local passes of Condition (II), each
attributable to a structural coincidence analysed in the paper.

## Layout

```
├── paper/                 LaTeX source of the paper
│   ├── wagstaff_chebyshev.tex
│   └── Makefile            pdflatex targets
├── scripts/                Reproducibility scripts (pure Python)
│   ├── README.md           script → paper section map
│   ├── survey.py           inert-factor survey
│   ├── audit.py            CSV integrity audit
│   ├── build_clean.py      canonicalise the raw survey CSV
│   ├── rerun_segments.py   re-run explicit survey segments
│   ├── primality_bls.py    BLS N-1 proof that W_12391 is prime
│   └── verify_sample.py    independent verifier for a CSV slice
├── data/                   Small sample data and hash file
│   ├── README.md           data dictionary + schema
│   ├── sample_1000.csv     first 1000 rows of the clean CSV
│   └── SHA256SUMS          hashes of full and sample CSVs
├── reproducibility.md      end-to-end reproduction walkthrough
├── LICENSE                 MIT
├── CITATION.cff            machine-readable citation
└── README.md               this file
```

## Cite this work

If you use the survey data or the computational scripts, please cite:

- **Paper (arXiv):** *pending submission — arXiv ID will be inserted here*
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
