# Reduction of the Wagstaff Chebyshev sufficiency conjecture to the Single-Pass Conjecture

LaTeX source, reproducibility scripts, and certified data for the paper
**"Reduction of the Wagstaff Chebyshev sufficiency conjecture to the
Single-Pass Conjecture"** by Alexey Dolotov (v4.0; formerly *"... to the Pell
Primitive Pair Conjecture"*, v3.6 — see *What changed in v4.0* below; the
last pre-restructure version, v3.8, is archived in `paper/archive/`).

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
the Vrba–Reix test") appears from v3.7 onwards, and a companion note works
out the exact relation. The mathematical content of the reduction is
unaffected.

## Abstract (v4.0)

For `W_p = (2^p + 1)/3` and `omega_3 = 3 + 2*sqrt(2)`, the congruence
`omega_3^{(W_p+1)/2} = -1 (mod W_p)` (Condition II) holds for every prime
`W_p`; whether it implies primality is the open **Wagstaff Chebyshev
sufficiency conjecture** — the Chua-form counterpart of the **Vrba–Reix**
Wagstaff primality test (open since 2008, standing 500-euro reward; see the
Erratum above on the exact relation).

Call a prime `r` with `ord_r(2) = 2p` a **local passer** at `p` if
`omega_3^{(W_p+1)/2} = -1` in `Z[sqrt(2)]/(r)`, and write `P_p` for the set
of local passers; every member of `P_p` divides `W_p`, and `P_p = {W_p}`
whenever `W_p` is prime. The **Single-Pass Conjecture (SPC)** asserts
`|P_p| <= 1` for every prime `p >= 5`. The headline theorem is that **SPC
implies the sufficiency conjecture, unconditionally and with zero
computational inputs** (the only outside ingredient is a Nagell–Ljunggren
result: `W_p` is never a perfect power, so a composite `W_p` has two distinct
prime factors, and under a global pass both would lie in `P_p`).

Double passes decompose by residues mod 8: the inert pair (above the Platinum
cutoff this is the **Pell Primitive Pair Conjecture, PPPC**), the mixed pair
(**MPC**), and the split pair (subsumed by general-`d` NCT). The working
chain proves the sufficiency conjecture from `PPPC + MPC^> + X1`, taking two
stated computational inputs (the **Platinum Lemma** — a 684,965,381-row
enumeration of all inert factors `r <= 10^12` — and the witness discharges);
that conjunction is in turn implied by SPC. Proved unconditionally along the
way: NCT for **every odd `d <= 400`** (28 fixed-`d` certificates in
`(100,400]`), a **complete danger-triple census through `d = 400`** (exactly
two real triples exist, both at exponents closed by secondary factors), the
Order-Pinning and Multi-Factor Pinning theorems, a Pair Separation theorem,
diagonal vanishing through `d = 400`, quartic residuosity of 2 at primitive
Pell divisors, and an exact-AP characterization. The reduction is
conditional; the sufficiency conjecture remains open. Across every record,
exactly four composite-`W_p` exponents with `P_p` nonempty are known — all
singletons, all inert — and no configuration with `|P_p| >= 2` has ever been
observed.

This is the companion to the Brillhart–Lehmer–Selfridge primality-proofs
paper (see *Cite this work*).

## What changed in v4.0

- **v4.0 (2026-07-16).** The paper is restructured around the Single-Pass
  Conjecture: new section "Local passes and the Single-Pass Conjecture"
  (the passer set `P_p`, SPC, and the input-free headline reduction), a
  branch-decomposition proposition, and a rewritten reduction section. The
  working chain `PPPC + MPC^> + X1` keeps hypotheses and computational
  inputs **identical to v3.8**, with a much shorter proof. The
  second-moment/bucket bookkeeping of earlier versions is removed; its
  surviving content is restated in sharper form (quantitative bound
  `N_in <= Pi + 1`; input-free pair-level form; triple-count erosion with a
  Wieferich-paired corner; Wieferich/multiplicity constraints). NCT is now
  closed for every odd `d <= 400` (v3.6 had fifteen of the nineteen values
  in `(200,400]`), and the complete danger census through `d = 400` is new.
  The last pre-restructure version is archived at
  `paper/archive/wagstaff_chebyshev_reduction_v3.8.{tex,pdf}`.
- **v3.6 (2026-06-25).** Previous bundle cut ("... to the Pell Primitive
  Pair Conjecture"), plus the 2026-07-15 Erratum above.

## Layout

```
.
├── paper/
│   ├── wagstaff_chebyshev_reduction_v4.0.tex   LaTeX source (67 pp)
│   ├── wagstaff_chebyshev_reduction_v4.0.pdf
│   ├── archive/                                v3.8 (last pre-SPC version)
│   └── Makefile                                pdflatex targets
├── scripts/                                    reproducibility scripts
│   └── README.md                               script -> paper statement map
│       NCT fixed-d certificates: cp350/cp351/cp353/cp356/cp358/cp367/
│         cp372/cp375/cp376 per-d closures; cp365 consolidated artifact;
│         cp352_verify_nct_recompute.py independent re-check;
│         cp356_aprcl_recert.py, cp358_aprcl_range.py APR-CL re-cert
│       Danger census d <= 400: cp412_danger_census_d400.py +
│         cp412b_t1_recheck.py, cp412b_t1_aprcl_sweep.py (adversarial pass)
│       Witness discharges / X1: cp352_verify_discharges_recompute.py
│       Diagonal + corner + falsification: cp354_samed_diagonal_certificates.py,
│         cp354_factordb_v2d_diagonal.py, cp354_diag_corner_gcd.py,
│         cp358_corner_sweep_d1000.py, cp358_falsification_prewindow.py
│       Heuristic constants (M1–M3): cp395_pppc_heuristic_constants.py
│       Self-supporting ladder: cp410_twin_ladder.py
│       Cross-case: platinum_lemma.py, multi_factor_pinning.py,
│         exact_ap_density.py; second_moment_reduction.py (historical)
│       Vrba–Reix agreement: cp361_vrba_reix_check.py
│       Independent verification: cp362_verify_preliminaries.py,
│         cp363_verify_inert_foundations.py
│       Survey pipeline: survey.py, build_clean.py, audit.py, verify_sample.py
├── data/
│   ├── nct_certificates.json      static fixed-d certificates: pinned Psi_{4d}
│   │                              factorizations + APR-CL + dispositions for
│   │                              all 28 closures in (100,400]
│   ├── cp412_danger_census.json   complete danger-triple census 43 < d <= 400
│   │                              (34 rungs, 340 tests, 226 APR-CL primes)
│   ├── danger_triple_data.json    V_{114}, V_{134}, V_{662} factorizations
│   ├── sample_1000.csv            first 1000 rows of the survey CSV
│   ├── SHA256SUMS                 hashes
│   └── README.md                  data dictionary
├── reproducibility.md             end-to-end reproduction walkthrough
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
every fixed-`d` closure in `(100, 400]` (28 values — with the `d <= 99`
theorem this closes NCT for every odd `d <= 400`), the complete factorization
of the Pell primitive part `Psi_{4d}` with each factor APR-CL-certified and
its certificate disposition. Re-generate / re-verify with
`python3 scripts/cp365_nct_certificate_bundle.py` (PARI/GP `gp` required;
expected: `DONE: 28/28 CLOSED`; 182 factors).

**Danger census `d <= 400`.** `data/cp412_danger_census.json` pins the
complete census of danger triples for admissible `43 < d <= 400`: 34 rungs,
complete `V_{2d}` factorizations, 340 (rung, prime) tests, 226 distinct
primes all APR-CL-certified; exactly two real triples (both at exponents
closed by secondary factors) plus two phantoms at prime `W_p`. The `cp412*`
scripts are the pinned provenance (they were run in the research tree; the
artifact is self-contained, and every listed prime can be re-checked with
`gp`'s `isprime(., 2)`).

**Witness discharges.** `python3 scripts/cp352_verify_discharges_recompute.py`
re-verifies end-to-end (self-contained, no network) that the two discharged
witness exponents carry an explicit prime factor at which Condition (II)
fails, and that the three Platinum witnesses are genuine local passes.

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
