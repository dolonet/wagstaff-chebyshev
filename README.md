# Reduction of the Wagstaff Chebyshev sufficiency conjecture to the Pell Primitive Pair Conjecture

LaTeX source, reproducibility scripts, and certified data for the paper
**"Reduction of the Wagstaff Chebyshev sufficiency conjecture to the Pell
Primitive Pair Conjecture"** by Alexey Dolotov (v3.6).

## Erratum (2026-07-15)

The v3.6 abstract and introduction call the conjecture "equivalently the open
'only if' direction of the Vrba–Reix Wagstaff primality test". **This is
incorrect as stated.** The literal Vrba–Reix prize iteration (`S_0 = 1/4
(mod 2^p+1)`, `S_{i+1} = S_i^2 - 2`) is a Lucas-type iteration in
`Q(sqrt(-7))` — its base is `theta = (1+3*sqrt(-7))/8 = (pibar/pi)^2` with
`2 = pi*pibar` — while Condition (II) lives in `Z[sqrt(2)]`. The two return
conditions decouple at the level of individual prime factors **in both
directions**, so the two sufficiency statements are logically incomparable:
proving either would not formally settle the other. They agree empirically at
every tested exponent. The corrected framing ("the Chua-form counterpart of
the Vrba–Reix test") appears in the v3.7 revision, and a companion note works
out the exact relation. The mathematical content of the paper — the reduction
of Condition-(II) sufficiency to PPPC + NCT-general + `W = 0` — is unaffected.

## Abstract

For `W_p = (2^p + 1)/3` and `omega_3 = 3 + 2*sqrt(2)`, the congruence
`omega_3^{(W_p+1)/2} = -1 (mod W_p)` (Condition II) holds for every prime
`W_p`; whether it implies primality is the open **Wagstaff Chebyshev
sufficiency conjecture** — the Chua-form counterpart of the **Vrba–Reix**
Wagstaff primality test (`S_{i+1} = S_i^2 - 2`; open since 2008, standing
500-euro reward; see the Erratum above on the exact relation). This paper
reduces the conjecture to a quantitative
residual: the **Pell Primitive Pair Conjecture (PPPC)**, a general-`d`
nonexistence-of-compatible-triples (NCT) statement on the split branch, and a
mild Wieferich–Wagstaff exclusion `W = 0` (unconditional for all factors below
`2^64`). On the way it proves, unconditionally: NCT for every odd `d <= 200`
and at fifteen further parity-unblocked values in `(200, 400]`; the
Order-Pinning, Multi-Factor Pinning, and Second-Moment reductions; a Pair
Separation Theorem; quartic residuosity of 2 at primitive Pell divisors; and an
exact-AP characterization. The inert branch collapses, via these reductions and
one computational input (the **Platinum Lemma**, a 684,965,381-row enumeration),
to PPPC. The reduction is conditional, not an unconditional proof.

This is the companion to the Brillhart–Lehmer–Selfridge primality-proofs paper
(see *Cite this work*).

## Layout

```
.
├── paper/
│   ├── wagstaff_chebyshev_reduction_v3.6.tex   LaTeX source (62 pp)
│   ├── wagstaff_chebyshev_reduction_v3.6.pdf
│   └── Makefile                                pdflatex targets
├── scripts/                                    reproducibility scripts
│   ├── README.md                               script -> paper section map
│   │   NCT fixed-d certificates (cor:nct-fixed-d-certificate):
│   │     cp351/cp353/cp356/cp358_nct_certificate*.py  per-d closures
│   │     cp356_aprcl_recert.py, cp358_aprcl_range.py  APR-CL re-cert
│   │     cp365_nct_certificate_bundle.py        consolidated static artifact
│   │   Inert/cross-case (Platinum, Second-Moment, pinning, exact-AP):
│   │     platinum_lemma.py, second_moment_reduction.py,
│   │     multi_factor_pinning.py, exact_ap_density.py, nct_*.py
│   │   Corner / falsification censuses:
│   │     cp354_diag_corner_gcd.py, cp358_corner_sweep_d1000.py,
│   │     cp358_falsification_prewindow.py
│   │   Vrba–Reix global agreement:  cp361_vrba_reix_check.py
│   │   Independent verification (preliminaries + inert foundations):
│   │     cp362_verify_preliminaries.py, cp363_verify_inert_foundations.py
│   │   Survey pipeline:  survey.py, build_clean.py, audit.py, verify_sample.py
├── data/
│   ├── nct_certificates.json   static fixed-d certificates: pinned Psi_{4d}
│   │                           factorizations + APR-CL + dispositions for all
│   │                           24 closures in (100,400] (FactorDB-independent)
│   ├── danger_triple_data.json V_{114}, V_{134}, V_{662} factorisations
│   ├── sample_1000.csv         first 1000 rows of the survey CSV
│   ├── SHA256SUMS              hashes
│   └── README.md               data dictionary
├── reproducibility.md          end-to-end reproduction walkthrough
├── CITATION.cff   .zenodo.json   LICENSE   README.md
```

## Cite this work

- **Paper (arXiv):** *pending submission — arXiv ID will be inserted here*
- **This bundle (Zenodo):** *DOI minted on release*
- **Inert-factor survey CSV (Zenodo):**
  [10.5281/zenodo.19496206](https://doi.org/10.5281/zenodo.19496206)
  (802 MB, 15,587,021 rows — the Platinum enumeration data)
- **Companion BLS primality-proofs paper (Zenodo):**
  [10.5281/zenodo.19645478](https://doi.org/10.5281/zenodo.19645478)

A machine-readable citation is in `CITATION.cff`.

## Reproduce the central claims

**NCT closures (no FactorDB needed).** `data/nct_certificates.json` pins, for
every fixed-`d` closure in `(100, 400]`, the complete factorization of the Pell
primitive part `Psi_{4d}` with each factor APR-CL-certified and its certificate
disposition. Re-generate / re-verify it with
`python3 scripts/cp365_nct_certificate_bundle.py` (PARI/GP `gp` required for
APR-CL; the embedded local factorizations are checked by exact product
identity).

**Vrba–Reix global agreement.** `python3 scripts/cp361_vrba_reix_check.py`
checks that, at every tested exponent, the `S^2-2` return tests and
Condition (II) give the same global verdict and track primality of `W_p`.
This is empirical agreement of global outcomes, **not** equivalence of the
tests — the per-factor conditions differ (see the Erratum above).

**Survey CSV (Platinum enumeration).** Download the 802 MB CSV from Zenodo
(`10.5281/zenodo.19496206`), verify `shasum -a 256 -c data/SHA256SUMS`, audit
with `python3 scripts/audit.py`, spot-check with
`python3 scripts/verify_sample.py data/sample_1000.csv`.

See `reproducibility.md` for the full procedure.

## License

MIT (code) / CC-BY-4.0 (paper, data). See `LICENSE`.
