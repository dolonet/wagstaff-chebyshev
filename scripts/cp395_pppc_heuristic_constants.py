#!/usr/bin/env python3
"""cp395: formalized PPPC pair-count heuristic -- explicit constants.

*** HEURISTIC ONLY. THIS SCRIPT PROVES NOTHING. ***
Every number below is a probabilistic model of deterministic objects.
PPPC remains NOT PROVED regardless of these numbers.

PURPOSE (track-2 draft cp395, paper 02 v3.7 sec:heuristic upgrade).
The paper states: "For pair counts the same heuristic, combined with the
Multi-Factor Pinning constraint that both factors share the pinned p
exactly, gives a total pair heuristic of O(10^-6)."  This script upgrades
that one-line figure to an explicit model-based expected-count computation
with numbered assumptions M1/M2/M3, per-exponent violation probabilities,
and a summed expectation with a ceiling/central/best bracket -- the same
presentation style as the v3.4 MPC joint heuristic (cp352), so the two
brackets are directly comparable.

MODEL (mirrors paper v3.7 sec:heuristic + cp352 machinery).

M1 (danger-triple production, d-indexed).  For admissible odd d > 43 the
expected number of danger triples (d, r, p) is the paper's own E(D)
summand with the inert fraction included:

  E1(d) = s_in * (lnln V_{2d})^2 / (2 ln V_{2d} * phi(L_d)),
  s_in = 1/2,  ln V_{2d} = 2 d ln(1+sqrt2),
  L_d  = valuation-adjusted modulus (paper Rem. valuation-adjusted-Ld).

  Decomposition: [# inert primitive divisors ~ s_in lnln V]
               * [Pr(pinned p_r = ord_r(2)/2 prime) ~ lnln V/(2 ln V)]
               * [Pr(condition (v): d | W_{p-2}) = 1/phi(L_d)].
  Dominant-factor convention ln p_bar(d) = ln V_{2d} - ln 4  (r = 2pk+1,
  k = O(1); enters only through slowly-varying factors).

M2 (condition-(v) density anchors, from the paper's measured records).
  delta_B (central) = 3 / 15,587,021 = 1.925e-7
      the Platinum enumeration's measured pass rate among the rows that
      survive the G_r <= 43 exclusion (15,587,018 order-excess rows + 3
      local passes; r <= 10^12).
  delta_A (ceiling) = 3 / 376,435 = 7.970e-6
      95%-CL rule-of-three upper bound from the above-cutoff zero-pass
      record: 369,876 (k<=500 window at 37.1M exponents) + 6,559 (deep-k
      window, 30,000 exponents) tested (i)-(iv) candidates, 0 passing (v).
  delta_C(r) (best) = delta_B * min(1, r0 / r), r0 = 10^12
      structural decay: the alignment d_r | W_{p-2} concentrates on
      candidates near the cutoff; extrapolating the measured constant to
      factors of size e^{2d} is the deliberate ceiling reading (same
      device as cp352 cols A/B vs C/D).

M3 (independence + partner count).  Distinct factors' (v)-events are
independent; the number of *other* locally passing inert factors at the
pinned exponent p is  q_part = min(1, omega_in(p) * delta),
omega_in(x) = (1/2)(ln x - lnln x)  (normal order, half inert; identical
to cp352's partner device).  Unordered-pair double count corrected by 1/2.
Unconditional repulsions the paper PROVES (pair separation Thm 6.10,
diagonal closures d <= 124, sub-Poisson variance ratio 0.977) are NOT
credited -- every neglected constraint pushes the true count DOWN.

SUMMED EXPECTATION:
  E[Pi]_x = (1/2) * sum_{d admissible > 43} E1(d) * q_part_x(d),
  x = A (ceiling, flat delta_A), B (central, flat delta_B),
  C (best: size-resolved anchor spectrum u in [ln 1e12, ln V_{2d}] with
     spectrum du/u, prime factor lnln/(2u) at ln r = u, partner density
     delta_C at partner scale max(r0, 2p)).

P-INDEXED PER-EXPONENT TABLE (independent cross-check, calibrated on the
paper's measured window means lambda(500) = 0.1148 and lambda(1e6) =
6559/30000 = 0.2186 at p* ~ 5e11):
  nu_p = per-exponent expected (v)-pass count among above-cutoff
         candidates;  Pr[pair at p] = nu_p^2/2  (Poisson).
"""

import math
from array import array

LN_ALPHA = math.log(1.0 + math.sqrt(2.0))       # 0.881373587...
LN2 = math.log(2.0)
N_MAX = 500_000                                  # d-sum range
D_MIN = 43
S_IN = 0.5                                       # inert fraction
R0 = 1.0e12                                      # Platinum cutoff
LN_R0 = math.log(R0)

# ---- M2 anchors (paper-stated records) --------------------------------
N_PLAT_TESTED = 15_587_018 + 3                   # order-excess rows + passes
DELTA_B = 3.0 / N_PLAT_TESTED                    # central (measured)
N_ZERO_PASS = 369_876 + 6_559                    # above-cutoff record
DELTA_A = 3.0 / N_ZERO_PASS                      # ceiling (95% rule of 3)

print("=" * 78)
print("cp395  PPPC pair heuristic, explicit constants -- HEURISTIC ONLY,")
print("       PROVES NOTHING.  Model M1/M2/M3; see docstring.")
print("=" * 78)
print(f"\n[M2] condition-(v) density anchors")
print(f"     delta_B (central) = 3/{N_PLAT_TESTED:,} = {DELTA_B:.3e}"
      f"   (Platinum, r <= 1e12)")
print(f"     delta_A (ceiling) = 3/{N_ZERO_PASS:,} = {DELTA_A:.3e}"
      f"     (95% UB, 0-pass record above cutoff)")
print(f"     delta_C (best)    = delta_B * min(1, 1e12/r)   (1/r decay)")

# ---------------------------------------------------------------- sieve
spf = array('i', range(N_MAX + 1))
for i in range(2, int(N_MAX ** 0.5) + 1):
    if spf[i] == i:
        for j in range(i * i, N_MAX + 1, i):
            if spf[j] == j:
                spf[j] = i

def fac(n):
    f = {}
    while n > 1:
        p = spf[n]
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        f[p] = e
    return f

_ord2 = {}
def ord2(q):
    if q in _ord2:
        return _ord2[q]
    o = q - 1
    for ell in fac(q - 1):
        while o % ell == 0 and pow(2, o // ell, q) == 1:
            o //= ell
    _ord2[q] = o
    return o

def v2(n):
    v = 0
    while n % 2 == 0:
        n //= 2
        v += 1
    return v

_inA = {}
def in_A(q):
    if q not in _inA:
        _inA[q] = (q == 3) or (v2(ord2(q)) == 1)
    return _inA[q]

_aq = {}
def lift_b(q, e):
    key = (q, e)
    if key in _aq:
        return _aq[key]
    m = ord2(q) // 2
    qe = q ** e
    x = (pow(2, m, qe) + 1) % qe
    if x == 0:
        b = 0
    else:
        a = 0
        while x % q == 0:
            x //= q
            a += 1
        b = max(0, e - a)
    _aq[key] = b
    return b

def phi_L(dfac):
    """phi of the valuation-adjusted modulus L_d; None if inadmissible."""
    L = {}
    for q, e in dfac.items():
        if not in_A(q):
            return None
        if q == 3:
            L[3] = max(L.get(3, 0), e)
            continue
        m = ord2(q) // 2
        b = lift_b(q, e)
        L[2] = max(L.get(2, 0), 1)
        for ell, c in fac(m).items():
            L[ell] = max(L.get(ell, 0), c)
        if b > 0:
            L[q] = max(L.get(q, 0), b)
    out = 1
    for ell, c in L.items():
        out *= (ell - 1) * ell ** (c - 1)
    return out

def omega_in(lnp):
    """expected # inert prime factors of W_p (cp352 device)."""
    return 0.5 * max(lnp - math.log(max(lnp, 3.0)), 1.0)

# ---------------------------------------------------------- [0] gate
print("\n[0] CALIBRATION GATE -- reproduce paper sec:heuristic E(1000)")
cal_raw = cal_half = 0.0
q = 1009
while q <= N_MAX:
    if spf[q] == q and in_A(q):
        lnV = 2.0 * q * LN_ALPHA
        lnln = math.log(lnV)
        m = ord2(q) // 2
        pL = 1
        for ell, c in fac(2 * m).items():
            pL *= (ell - 1) * ell ** (c - 1)
        t = lnln * lnln / (2.0 * lnV * pL)
        cal_raw += t
        cal_half += 0.5 * t
    q += 2
print(f"     prime-level sum, admissible q in (1000, {N_MAX}]:")
print(f"     raw (no inert fraction)  = {cal_raw:.4f}   paper: ~0.0076")
print(f"     with s_in = 1/2          = {cal_half:.4f}   paper: ~0.003")
gate0 = abs(cal_raw - 0.0076) < 0.002 and abs(cal_half - 0.003) < 0.001
print(f"     gate: {'PASS' if gate0 else 'FAIL'}")

# ------------------------------------------------- [1] p-indexed nu_p
print("\n[1] P-INDEXED CALIBRATION (window means, p* ~ 5e11)")
P_STAR = 5.0e11
LAM_500 = 0.1148          # paper sec:empirical, 37,123,274 exponents
LAM_1E6 = 6559.0 / 30000  # deep-k record, 30,000 exponents

def lam_raw(K, p):
    """sum over k <= K, k in one class mod 4 (avg over the two odd
    classes), of (1/2k) * 1/ln(2pk+1)."""
    tot = 0.0
    for a in (1, 3):
        s = 0.0
        k = a
        while k <= K:
            s += 1.0 / (2.0 * k * math.log(2.0 * p * k + 1.0))
            k += 4
        tot += 0.5 * s
    return tot

raw500 = lam_raw(500, P_STAR)
C_CAL = LAM_500 / raw500
pred1e6 = C_CAL * lam_raw(1_000_000, P_STAR)
print(f"     model lambda(K) = C_cal * avg_a sum_(k<=K, k=a mod 4)"
      f" 1/(2k ln(2pk))")
print(f"     C_cal = {LAM_500}/{raw500:.5f} = {C_CAL:.3f}"
      f"   (absorbs AP/Mertens + mod-8 concentration)")
print(f"     validation: predicted lambda(1e6) = {pred1e6:.4f}"
      f"   measured = {LAM_1E6:.4f}"
      f"   ratio {pred1e6/LAM_1E6:.2f}")
gate1 = 0.7 < pred1e6 / LAM_1E6 < 1.4
print(f"     gate ([0.7,1.4]): {'PASS' if gate1 else 'FAIL'}")

def nu_p(p, mode):
    """per-exponent expected (v)-pass count among above-cutoff candidates.
    mode: 'A'/'B' flat delta over the full factor spectrum (Erdos-Kac
    count omega_in with the k >= k_min cutoff), 'C' decay (calibrated
    k-sum, converges)."""
    lnp = math.log(p)
    if mode in ('A', 'B'):
        delta = DELTA_A if mode == 'A' else DELTA_B
        # factors r = 2pk+1 > 1e12: spectrum lnln W_p - lnln(max(2p,1e12))
        lnlnW = math.log(p * LN2)
        floor = math.log(max(math.log(2.0 * p), LN_R0))
        mu = 0.5 * max(lnlnW - floor, 0.0)
        return mu * delta, mu
    # mode C: sum_k q1(k) * delta_C(2pk+1); explicit head + integral tail
    kmin = max(1, int((R0 - 1) / (2 * p)) + 1)
    tot = 0.0
    NTERMS = 20_000
    for a in (1, 3):
        k = kmin + ((a - kmin) % 4)      # least k >= kmin, k = a mod 4
        s = 0.0
        for _ in range(NTERMS):
            r = 2.0 * p * k + 1.0
            s += (1.0 / (2.0 * k * math.log(r))) * min(1.0, R0 / r)
            k += 4
        # tail: sum_{t>=k, t=a(4)} ~ (1/4) int R0/(4 p t^2 ln 2pt) dt
        s += R0 / (16.0 * p * k * math.log(2.0 * p * k))
        tot += 0.5 * s
    return C_CAL * DELTA_B * tot, None

print("\n     per-exponent nu_p and pair probability nu^2/2:")
print(f"     {'p':>10s} {'nu_A':>10s} {'nu_B':>10s} {'nu_C':>10s}"
      f" {'pair_A':>10s} {'pair_B':>10s} {'pair_C':>10s}")
for p in (5.0e11, 1.0e13, 1.0e15, 1.0e18):
    nA, muA = nu_p(p, 'A')
    nB, _ = nu_p(p, 'B')
    nC, _ = nu_p(p, 'C')
    print(f"     {p:10.0e} {nA:10.2e} {nB:10.2e} {nC:10.2e}"
          f" {nA*nA/2:10.2e} {nB*nB/2:10.2e} {nC*nC/2:10.2e}")
_, mu_star = nu_p(P_STAR, 'B')
print(f"     (mu_p at p* = 5e11: {mu_star:.1f} above-cutoff factors"
      f" expected; EK/normal-order model)")

# consistency with the zero-pass record
exp_500win = 369_876 * DELTA_A, 369_876 * DELTA_B
print(f"\n     consistency: expected (v)-passes among the 369,876 tested")
print(f"     candidates: ceiling {exp_500win[0]:.2f} (the 95% bound's own")
print(f"     definition), central {exp_500win[1]:.3f}; observed 0.")

# sum over p of nu_C^2/2 (best-estimate p-indexed total)
print("\n     summed p-indexed best estimate  sum_p nu_C(p)^2 / 2 :")
# integrate over primes via dx/ln x, log grid; two regimes around 5e11
def sum_pairs_C():
    tot = 0.0
    lo, hi, ng = math.log(5.0), math.log(1.0e18), 300
    h = (hi - lo) / ng
    prev = None
    x_prev = None
    for i in range(ng + 1):
        x = math.exp(lo + i * h)
        n, _ = nu_p(x, 'C')
        f = (n * n / 2.0) / math.log(x)     # prime density
        if prev is not None:
            tot += 0.5 * (prev + f) * (x - x_prev)
        prev, x_prev = f, x
    return tot

SUMC = sum_pairs_C()
print(f"     integral p in [5, 1e18] (negligible tail beyond: nu_C ~ 1/p):")
print(f"     E[Pi]_p-indexed,best = {SUMC:.2e}")
print(f"     NOTE: under flat delta (A/B) the same sum DIVERGES like")
print(f"     sum_p (lnln W_p)^2 -- the model face of the paper's remark")
print(f"     that even the single-pass first moment is divergent")
print(f"     (rem:triple-status).  Flat readings are therefore reported")
print(f"     per-exponent and through the d-indexed capped sum below,")
print(f"     never as an uncapped p-sum.")

# --------------------------------------------- [2] d-indexed E[Pi]
print("\n[2] D-INDEXED SUMMED EXPECTATION  E[Pi] = (1/2) sum_d E1(d)*q_part")
totA = totB = totC = 0.0
winA = winB = winC = 0.0            # tail-window sums
sum_E1 = 0.0
WIN_LO = 400_000
n_adm = 0
rows = []
d = D_MIN + 2
while d <= N_MAX:
    pL = phi_L(fac(d))
    if pL is not None:
        n_adm += 1
        lnV = 2.0 * d * LN_ALPHA
        lnln = math.log(lnV)
        E1 = S_IN * lnln * lnln / (2.0 * lnV * pL)
        sum_E1 += E1
        lnp = lnV - math.log(4.0)
        w = omega_in(lnp)
        qA = min(1.0, w * DELTA_A)
        qB = min(1.0, w * DELTA_B)
        # best: size-resolved anchor spectrum u in [ln 1e12, lnV]
        # integrand (per u): [du/u spectrum] * [lnln r/(2u) prime]
        #                    * [min(1, omega_in * delta_C(partner))]
        # partner scale max(r0, 2p), ln 2p = u - ln 2  (k = O(1))
        EC = 0.0
        u_lo, u_hi = LN_R0, lnV
        if u_hi > u_lo * 1.0001:
            ng = 120
            h = (u_hi - u_lo) / ng
            prev = None
            for i in range(ng + 1):
                u = u_lo + i * h
                ln2p = u - LN2
                dpart = DELTA_B * min(1.0, math.exp(LN_R0 - max(ln2p, LN_R0)))
                f = (math.log(u) / (2.0 * u * u)) * \
                    min(1.0, omega_in(max(ln2p - LN2, 2.0)) * dpart)
                if prev is not None:
                    EC += 0.5 * (prev + f) * h
                prev = f
        EC *= S_IN / pL
        cA = 0.5 * E1 * qA
        cB = 0.5 * E1 * qB
        cC = 0.5 * EC
        totA += cA; totB += cB; totC += cC
        if d >= WIN_LO:
            winA += cA; winB += cB; winC += cC
        rows.append((cA, d, pL, E1, qA, qB))
        if len(rows) > 4000:
            rows.sort(reverse=True)
            del rows[40:]
    d += 2

rows.sort(reverse=True)
n_odd_win = (N_MAX - WIN_LO) // 2
NBAR = 0.5 * (WIN_LO + N_MAX)

# tails: A/B summands ~ ln^2 t/t until the cap saturates (A: saturated
# well below N_MAX so tail ~ E1-tail ~ 1/t^2; B: cap at d* ~ 1.2e7).
# crude device (cp352 style): window mean per odd d * shape * x4 slack.
tailA = 4.0 * (winA / n_odd_win) * NBAR ** 2 / (2.0 * N_MAX)
d_star_B = 1.0 / (0.44 * DELTA_B)     # omega_in(1.76 d) delta_B ~ 1
# B pre-cap growth ~ ln^2: window-mean linear-log extrapolation to d*_B,
# then post-cap E1-tail beyond d*_B:
meanB = winB / n_odd_win
tailB_precap = meanB * (d_star_B - N_MAX) / 2.0 * \
    (math.log(d_star_B) / math.log(NBAR)) ** 2 / (d_star_B / NBAR)
# summand_B ~ c*ln^2(t)/t: integral N..d* = c/3 (ln^3 d* - ln^3 N);
# with c ~ meanB * NBAR / ln^2 NBAR:
cB_fit = meanB * NBAR / math.log(NBAR) ** 2
tailB_precap = cB_fit / 3.0 * (math.log(d_star_B) ** 3 -
                               math.log(N_MAX) ** 3)
# post-cap: (1/2) sum_{d>d*} E1 ~ (1/2) * 0.28-admissible * S_IN *
# ln^2(1.76d*)/(2*1.76*d*) * (2/avg phiL/d ~ 2/d) integrated:
adm_frac = n_adm / ((N_MAX - D_MIN) / 2)
lnVs = 2.0 * d_star_B * LN_ALPHA
postB = 0.5 * adm_frac * S_IN * math.log(lnVs) ** 2 / (2.0 * LN_ALPHA) \
    * 2.0 / d_star_B
tailB = tailB_precap + postB
tailC = 4.0 * (winC / n_odd_win) * NBAR ** 2 / (2.0 * N_MAX)

EA = totA + tailA
EB = totB + tailB
EC_ = totC + tailC
print(f"     admissible d in ({D_MIN}, {N_MAX}]: {n_adm}"
      f" ({adm_frac:.1%} of odd)")
print(f"     sum E1(d) (single danger triples, d > 43) = {sum_E1:.4f}")
print(f"     top contributors:")
for cA, dd, pL, E1, qA, qB in rows[:8]:
    print(f"       d={dd:<7d} phi(L)={pL:<8d} E1={E1:.2e}"
          f"  qA={qA:.2e} qB={qB:.2e}  contribA={cA:.2e}")
print(f"\n     E[Pi] ceiling (A, flat {DELTA_A:.1e}):"
      f" {totA:.3e} + tail {tailA:.1e} = {EA:.3e}")
print(f"     E[Pi] central (B, flat {DELTA_B:.1e}):"
      f" {totB:.3e} + tail {tailB:.1e} = {EB:.3e}")
print(f"       (B tail: pre-cap ln^2 growth to d* ~ {d_star_B:.1e},"
      f" then E1-tail {postB:.1e})")
print(f"     [diagnostic] d-indexed size-resolved C:"
      f" {totC:.3e} + tail {tailC:.1e} = {EC_:.3e}")
print(f"       REJECTED as a headline: anchoring every d at the cutoff")
print(f"       scale double-counts primes across d (each inert prime has")
print(f"       ONE d_r; the per-d du/u spectrum may not be summed over d")
print(f"       at fixed size).  The p-indexed sum of block [1], which")
print(f"       counts each candidate once via its k-position, is the")
print(f"       correct decay-model total: {SUMC:.2e}")

# --------------------------------------------- [3] headline + context
print("\n[3] HEADLINES (heuristic only; PROVES NOTHING)")
print(f"     E[Pi] ceiling  = {EA:.1e}   (d-indexed, flat delta_A, capped)")
print(f"     E[Pi] central  = {EB:.1e}   (d-indexed, flat delta_B, capped;")
print(f"                      tail-dominated, read as ~1e-4)")
print(f"     E[Pi] best     = {SUMC:.1e}   (p-indexed, structural decay;")
print(f"                      dominated at the anchor scale p ~ 5e11)")
print(f"     paper's existing figure: 'total pair heuristic of O(1e-6)'")
print(f"     inside bracket [best, ceiling]:"
      f" {'yes' if SUMC <= 1e-6 <= EA else 'NO -- check'}")
print(f"     MPC joint bracket (cp352, for comparison): 8.7e-4 / 4.4e-4 /"
      f" 9.8e-10")
print(f"     per-exponent at the Platinum frontier p* = 5e11:")
nA, _ = nu_p(P_STAR, 'A'); nB, _ = nu_p(P_STAR, 'B'); nC, _ = nu_p(P_STAR, 'C')
print(f"       Pr[pair] ceiling {nA*nA/2:.1e} / central {nB*nB/2:.1e}"
      f" / best {nC*nC/2:.1e}")

print("""
[4] CAVEATS (load-bearing; heuristic, NOT a proof)
 1. Randomness model for deterministic objects (divisor counts lnln V,
    1/ln prime density, AP density 1/phi(L_d)); nothing is rigorous.
 2. Flat delta (A/B) extrapolates a density measured at r <= 1e12 to
    factors of size e^(2d); with it even the single-pass first moment
    diverges over p (consistent with paper rem:triple-status), so flat
    readings are meaningful only per-exponent or inside the d-indexed
    capped sum.  The six-orders ceiling/best spread is genuine model
    uncertainty, reported bracketed, not resolved (same as the MPC
    heuristic's spread).
 3. Uncredited PROVED repulsions (pair separation, diagonal closures at
    admissible d <= 124, sub-Poisson ratio 0.977) all push the true count
    below the model.  Diagonal share within the model is not separated.
 4. min(1, .) caps make the flat d-indexed sums converge; the B tail
    crosses its cap at d* ~ 1.2e7 > N_MAX, so the B headline carries the
    extrapolated pre-cap + post-cap tail estimate (crude, x2-ish).
 5. Independence across the two factors at one p (M3): both divisibilities
    d_i | W_{p-2} live on the same integer; common prime factors of the
    d_i correlate them (worst case ~x2), and gcd constraints (Thm
    pair-separation (3)) anti-correlate.  Orders of magnitude unaffected.
 6. C_cal calibration at p* ~ 5e11 validated at 1.9x window growth
    (gate [1]); using it at p far from 5e11 assumes the same AP/Mertens
    constant.  The p-indexed best total is dominated within one decade of
    the anchor scale, so this extrapolation is mild exactly where it
    matters.
 7. The d-indexed size-resolved decay variant (diagnostic line, block
    [2]) is inflated by cross-d double counting near the cutoff and is
    not a headline; the flat d-indexed readings (A/B) use the paper's
    dominant-divisor convention (anchor at size ~ V_{2d}, distinct per
    d), which does not suffer from this.
END (heuristic only; nothing above is PROVED)""")
