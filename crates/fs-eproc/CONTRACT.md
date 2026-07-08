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
- `PairwiseRace` — the e-racing primitive: paired losses → d = midpoint of
  clipped difference and 1 → betting process at null ½; `a_beats_b(alpha)`.
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
- Trajectories are bit-replayable: strict fs-math functions + caller-supplied
  fs-rand streams (tested: identical tournament decision, time, and bits).
- log-space wealth: no overflow at extreme evidence.

## Error model
Contract violations (null mean outside (0,1), out-of-range outcomes,
non-positive σ/ρ, α∉(0,1)) panic with teaching messages — these are
guarantee-voiding programmer errors, not data errors.

## Determinism class
Deterministic CROSS-ISA (strict fs-math + pure arithmetic).

## Cancellation behavior
O(1) per observation; the CALLER owns kill-handles (fs-exec) — this crate
supplies the decision signal e-racing acts on.

## Unsafe boundary
None.

## Feature flags
None.

## Conformance tests (empirical, seeded, release-mode)
Ville validity under ADVERSARIAL stopping (4000 sims × 2000 horizon:
type-I 0.0305 ≤ 0.05); power scaling (median stop 241 @ δ=.05, 32 @ δ=.15);
time-uniform CS coverage (miss 0.0487 ≤ 0.05; radius 0.573→0.086);
e-BH FDR 0.001 ≤ 0.1 at power 1.00 (10 signals / 30 nulls × 300 sims);
race decided + bit-replayable; arithmetic laws; bounded-input refusal.

## No-claim boundaries
- Empirical-Bernstein/hedged closed-form CS for bounded means (mixture CS
  covers sub-Gaussian; EB variant is follow-up scope with its consumers).
- Conformal e-prediction (fs-surrogate + conformal-hardening beads).
- Nightly adversarial-rotation trials (certify-certifiers bead).
- Two-sided betting tests; asymptotic variants; sub-exponential extensions.

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
