# CONTRACT: fs-eproc

## Purpose and layer
Anytime-valid inference (plan §9.6, Bet 5): betting e-processes, mixture
confidence sequences, e-value arithmetic, e-BH — the statistical layer that
makes optional stopping VALID. Layer: L4.

## Public types and semantics
- `BettingEProcess::new(null_mean)` — one-sided test of H₀: mean ≤ m for
  outcomes in [0,1] (WSR betting lineage): log-wealth via strict fs-math ln;
  PREDICTABLE plug-in bet (regularized empirical mean/var, clipped one-sided
  to `aggressiveness / max(m, 1−m)`); `observe`, `e_value`, `log_e_value`,
  `rejects_at(alpha)` (Ville threshold 1/α).
- `LossSpan::new(span)` — checked finite-positive, fixed-before-observation
  support for `abs(loss_b - loss_a)`. The span is part of the statistical
  claim and replay identity; it must not be estimated from the race data.
- `PairwiseRace::new(loss_span)` — the e-racing primitive for
  lower-is-better losses. Construction requires the scale explicitly; there
  is no default unit span. It maps
  `d = ((loss_b - loss_a) / span + 1) / 2` without clipping and feeds a
  betting process at null 1/2; `a_beats_b(alpha)` reports evidence that A has
  lower conditional mean loss than B. The required null is
  `E[loss_b,t - loss_a,t | F_(t-1)] <= 0`, not marginal-mean equality alone.
- `GaussianMixtureCs::new(sigma, rho, alpha)` — Robbins normal-mixture CS
  for sub-Gaussian means: time-uniform interval with closed-form radius
  √((tσ²+ρ)(ln((tσ²+ρ)/ρ) + 2ln(1/α)))/t; also `e_value_for(m)` (the
  mixture martingale as a two-sided e-value).
- Arithmetic: `combine_product` (independent), `combine_average` (valid
  under ARBITRARY dependence; max-shifted log-sum-exp), `e_to_p` (Markov).
- `e_benjamini_hochberg(log_e, alpha)` — Wang–Ramdas e-BH: reject k̂ largest
  where e_(k) ≥ m/(αk); FDR ≤ α under arbitrary dependence; deterministic
  tie-breaking (descending e, ascending index).

- `hardening` module (patch Rev M, bead 7tv.9; [F], behind
  `conformal-hardening`): the anytime-valid layer's assumptions as
  operational contracts. `MondrianConformal` (per-bucket split
  conformal at the ⌈(n+1)(1−α)⌉ rule; buckets below the calibration
  floor REFUSE with a teaching count instead of extrapolating;
  `marginal_band` kept for the comparison), `DriftMonitor` (sequential
  two-sample PIT betting pair against the COMPOSITE null
  mean ∈ 1/2 ± δ, δ = max(1/√n_train, 0.02) — the finite-calibration
  tolerance, added after an unslacked monitor false-fired on the
  training sample's own noise; detection shrinks `validity_scale`),
  `CoverageClaim`/`fcr_flag` (per-claim miscoverage e-processes,
  e-BH over them spends the explicit false-coverage budget),
  `admission_alpha` (the Bonferroni reservation for study admission),
  `ExchangeabilityCard` (assumptions declared, ledger-ready).

## Invariants
- Validity is STRUCTURAL: any predictable bet in the admissible range yields
  a supermartingale under H₀ — strategy affects power only.
- Boundedness contract enforced: observations outside [0,1] are REFUSED
  (feeding unbounded data would silently void the guarantee).
- Pairwise support is enforced before wealth changes: non-finite losses,
  subtraction overflow, and differences outside the declared `LossSpan`
  return `PairwiseInputError`. Clipping is forbidden because it changes the
  estimand and can create betting drift under an equal-raw-mean null.
- Trajectories are bit-replayable: strict fs-math functions + caller-supplied
  fs-rand streams (tested: identical tournament decision, time, and bits).
- log-space wealth: no overflow at extreme evidence.

## Error model
Malformed scalar process parameters (null mean outside (0,1), out-of-range
bounded-process outcomes, non-positive sigma/rho, alpha outside (0,1)) panic
with teaching messages. Pairwise loss data and span violations are fallible
`PairwiseInputError` values and leave pairwise wealth unchanged, because a
production race must be able to surface a no-claim result rather than panic.

## Determinism class
Deterministic CROSS-ISA (strict fs-math + pure arithmetic).

## Cancellation behavior
O(1) per observation; the CALLER owns kill-handles (fs-exec) — this crate
supplies the decision signal e-racing acts on.

## Unsafe boundary
None.

## Feature flags
`conformal-hardening` enables the [F] Mondrian conformal, drift-monitor, and
false-coverage-budget machinery. It is default-off.

## Conformance tests (empirical, seeded, release-mode)
Ville validity under ADVERSARIAL stopping (4000 sims × 2000 horizon:
type-I 0.0305 ≤ 0.05); power scaling (median stop 241 @ δ=.05, 32 @ δ=.15);
time-uniform CS coverage (miss 0.0487 ≤ 0.05; radius 0.573→0.086);
e-BH FDR 0.001 ≤ 0.1 at power 1.00 (10 signals / 30 nulls × 300 sims);
race decided + bit-replayable; arithmetic laws; bounded-input refusal;
declared-span boundaries and one-ULP refusal; equal-mean skew counterexample.

## No-claim boundaries
- Empirical-Bernstein/hedged closed-form CS for bounded means (mixture CS
  covers sub-Gaussian; EB variant is follow-up scope with its consumers).
- Conformal e-prediction (fs-surrogate + conformal-hardening beads).
- Nightly adversarial-rotation trials (certify-certifiers bead).
- Two-sided betting tests; asymptotic variants; sub-exponential extensions.
- A caller-supplied finite span is a checked assumption, not a certificate that
  the underlying stochastic process has that support. Runtime breach refuses
  the observation; justifying the bound remains the caller's obligation.

## No-claim boundaries (hardening)

- The drift monitor tests the CANDIDATE SCORE distribution (1-D PIT);
  multivariate covariate drift needs a score projection chosen by the
  caller — the projection's blind spots are the caller's declaration.
- The composite-null slack δ trades sensitivity below δ for
  no-false-alarm robustness: shifts smaller than the calibration
  noise floor are undetectable BY DESIGN.
- FCR flagging inherits e-BH's guarantee under arbitrary dependence;
  the admission bound is Bonferroni (conservative on purpose).
- Weighted-conformal handling of optimization-induced shift (beyond
  detect-and-refuse/recalibrate) is the growth path the bead names.
