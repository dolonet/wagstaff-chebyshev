# Data

This directory contains a small sample slice of the inert-factor survey
that underpins the computational results in the paper. The full dataset
is hosted on Zenodo (see below).

## Files in this directory

| File | Description |
|---|---|
| `sample_1000.csv` | First 1000 rows of the canonical clean CSV. Smoke test for `scripts/verify_sample.py`. |
| `SHA256SUMS` | SHA-256 hashes of `inert_factors.csv` (full) and `sample_1000.csv`. |

## The full dataset

The full clean CSV is **`inert_factors.csv`**, 802 MB, 15 587 021 rows.
It is distributed separately via Zenodo:

> **Zenodo DOI:** *(pending upload — see the paper for the final citation)*

Verify the download with:

```sh
shasum -a 256 -c SHA256SUMS
```

The expected hash of the full CSV is

```
db2b08d8542b0e5c46cf5c6b54077d3122f297eadb2f3735cf900da89eb63e37
```

## Schema

Every row in `inert_factors.csv` describes one inert prime factor
`r` ≡ 3 (mod 8) of some composite Wagstaff number `W_p = (2^p+1)/3`
with `p < 5×10¹¹` and `r < 10¹²`. Columns:

| # | Column | Type | Meaning |
|---|---|---|---|
| 1 | `p` | int | Exponent. Prime, and `W_p` is composite. |
| 2 | `r` | int | Inert prime factor of `W_p`, `r ≡ 3 (mod 8)`. |
| 3 | `G4` | int | `G_r/4 = gcd((r+1)/4, W_{p-2})`. Computed by modular arithmetic. |
| 4 | `d` | int | `ord_r(ω₃)/4`, with `ω₃ = 3+2√2`. `0` means the row short-circuited via an order-excess check (no full order computed). |
| 5 | `d_factors` | str | Factorisation of `d`, e.g. `3*5419*2350291`. Empty if `d = 0`. |
| 6 | `pure_iii` | bool | All prime factors of `d` are in Class III (`v₂(ord_q(2)) = 1`). |
| 7 | `obstruction` | str | Why Condition (II) fails locally. One of: `order_excess`, `class_I`, `class_II`, `congruence`, `valuation`, `trivial`, `PASS`. |
| 8 | `blocking_q` | int | The prime factor of `d` that caused the obstruction (if applicable, else `0`). |
| 9 | `blocking_class` | str | Class of `blocking_q`: `I`, `II`, or `III`. |

### The `order_excess` shortcut

A row with `d = 0` and `obstruction = order_excess` used the following
shortcut: the survey first checks whether `ω₃^(4·G_r/4) = 1 (mod r)`. If
not, the order of ω₃ cannot divide `4·G_r/4`, hence no local pass of
Condition (II) is possible for this factor and the full order
computation is skipped. The vast majority of inert factors take this
branch.

### The three local passes

Exactly three rows in the full CSV have `obstruction = PASS`:

| `p` | `r` | `G_r/4 = d` | mechanism |
|---|---|---|---|
| 1 251 488 009 | 2 502 976 019 | 41 716 267 | primitive-root coincidence |
| 10 916 765 939 | 152 834 723 147 | 38 208 680 787 | structured cofactor |
| 85 684 865 933 | 171 369 731 867 | 57 | small-`d` gap at `d = 57` |

Each row is analysed in the paper (§9.8).

## How the data was generated

Run `scripts/survey.py` (pure-Python, no dependencies). The survey is
fully reproducible at segment granularity; see `reproducibility.md`
in the repo root.
