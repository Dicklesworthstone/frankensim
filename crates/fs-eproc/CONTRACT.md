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
  under ARBITRARY dependence; extended-real-safe max-shifted log-sum-exp),
  `e_to_p` (Markov). NaN and indeterminate zero-times-infinity products are
  refused before they can become evidence.
- `e_benjamini_hochberg(log_e, alpha)` — Wang–Ramdas e-BH: reject k̂ largest
  where e_(k) ≥ m/(αk); FDR ≤ α under arbitrary dependence; deterministic
  tie-breaking (descending e, ascending index).

- `hardening` module (patch Rev M, bead 7tv.9; [F], behind
  `conformal-hardening`): the anytime-valid layer's assumptions as
  operational contracts. `MondrianConformal` (per-bucket split
  conformal at the ⌈(n+1)(1−α)⌉ rule; buckets below the calibration
  floor REFUSE with a teaching count instead of extrapolating;
  `marginal_band` kept for the comparison), `DriftMonitor` (sequential
  two-sample PIT betting pair, combined as one equal-weight mixture e-process
  so the two-sided decision spends alpha once, against the COMPOSITE null
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
- Confidence-sequence observations, nulls, calibration residuals, and drift
  samples are finite; malformed alpha/budget values are refused. A rejected
  observation cannot mutate running state.
- Pairwise support is enforced before wealth changes: non-finite losses,
  subtraction overflow, and differences outside the declared `LossSpan`
  return `PairwiseInputError`. Clipping is forbidden because it changes the
  estimand and can create betting drift under an equal-raw-mean null.
- Trajectories are bit-replayable: strict fs-math functions + caller-supplied
  fs-rand streams (tested: identical tournament decision, time, and bits).
- log-space wealth: no overflow at extreme evidence.

## Error model
Malformed scalar process parameters (null mean outside (0,1), out-of-range
bounded-process outcomes, non-finite confidence-sequence observations,
non-finite/non-positive sigma/rho, alpha outside (0,1), malformed calibration
data) panic with teaching messages. Pairwise loss data and span violations are fallible
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

The five empirical Gauntlet aggregates use canonical fs-obs
`ConformanceCase` rows under `fs-eproc`, with distinct `ville-validity`,
`power`, `cs-coverage`, `e-bh`, and `race-replay` scopes. Every row records the
actual shared fs-rand input root `0xB7E1` (`0xE90C ^ 0x5EED`) in its `seed`
field. Kernel 1 and the logical input tiles are: 0..3999 for Ville validity;
1000..1199 and 2000..2199 for power; 50000..51499 plus shrinkage tile 99999
for CS coverage; `200000 + sim*64 + hypothesis` for e-BH; and tile 777 twice
for race replay. These are logical RNG coordinates, never worker, thread, Cx,
or scheduler seeds. Each fresh emitter produces one sequence-zero event and
runs the failure-record lint plus fs-obs wire validation. The rows remain after
their existing statistical and replay assertions, so a failed gate aborts
before emission as an ordinary Rust diagnostic rather than a structured
failure record.

`tests/hardening.rs`, behind `conformal-hardening`, emits canonical aggregate
`ConformanceCase` verdicts for ch-001..ch-005 and validated object-shaped
`Custom` measurement companions for ch-001..ch-004. The literal input roots
are `0xa11ce`, `0xd21f7`, `0xfc12`, `0x0b71`, and `0x5eed`, respectively;
the suite has no execution/Cx seed. Assertions and expectations reached before
an aggregate verdict remain ordinary Rust test diagnostics. The fixed
regressions ch-006..ch-008 remain assertion-only and do not claim aggregate
event coverage. Central proof must explicitly enable
`fs-eproc/conformal-hardening`; a default-feature pass skips the target.

`tests/adversarial_rotation.rs` is the fixed seeded e-process
certifier-of-certifier tranche in `frankensim-epic-ascent-7tv.21.9`, advancing
the e-process scope of `frankensim-epic-gauntlet-6nb.7`. Four
strategies choose predictably among four dyadic bounded null laws whose
conditional means are exactly one half, then stop at the Ville threshold,
their declared history-dependent rule, or a 1024-observation horizon. Each
strategy uses 1024 disjoint fs-rand logical streams and must record at most 80
level-0.05 crossings (a greater-than-four-binomial-standard-deviation finite
gate), exact law/stop/observation accounting, and a non-colliding tile range.
One canonical configuration identity binds the laws, strategies, stopping
rules and numeric selector thresholds, alpha, seed, Philox kernel/tile/index
layout, draw method/count, stream-semantics version, and crate versions; a
result identity and object-shaped `fs-obs` receipt retain all four empirical
outcomes before five linted, wire-validated
`ConformanceCase` verdicts are emitted.

`tests/ebh_dependence.rs` is the G0 e-BH certifier-of-certifier tranche in
`frankensim-epic-ascent-7tv.21.10`. It exhaustively enumerates three exact
finite laws over 4,096 simultaneous hypotheses. In the global-null law, 64
perfectly dependent blocks of 64 hypotheses are mutually exclusive across
1,024 equiprobable phases; an active block has e-value 1,024, so every null's
mean e-value is exactly one and FWER = FDR = 64/1024 = 1/16. In the mixed
perfect-dependence law, 4,076 true nulls share one e=17 shock among 17 phases
beside 20 fixed alternatives, giving exact FDR
`4076 / (4096 * 17) <= 1/16`. A third exact law gives each of 63 mutually
exclusive 64-null blocks e=63 once beside 64 alternatives: production e-BH
rejects no nulls, while a test-local seeded-broken comparator that omits the
family-size factor from `m/(alpha*k)` rejects the active block and incurs FDP
1/2 in every phase. The latter implementation is mutation-test code only and
is never used by production paths.

One fs-rand seed chooses a one-draw cyclic rotation of original hypothesis
IDs; all probability phases are then enumerated rather than sampled. A
canonical configuration identity binds the laws, phase order, evidence and
threshold values, family/block geometry, alpha, rotation coordinates and
algorithm, stream-semantics version, comparator mutation, and crate versions.
The result identity and object-shaped `fs-obs` receipt retain exact rejection,
false-rejection, marginal-validity, and first-shape-failure fields before four
linted, wire-validated `ConformanceCase` verdicts are emitted.

`tests/conformal_shift_receipt.rs`, behind `conformal-hardening`, is the exact
finite conformal-shift and surrogate-policy tranche in
`frankensim-epic-ascent-7tv.21.19`. It exhausts all 20 leave-one-out ranks in
each of two separated residual regimes at alpha 0.1. Per-regime Mondrian bands
must cover exactly 18/20 ranks in both buckets. The pooled band covers 36/40
overall while covering only 16/20 shock ranks, so nominal marginal coverage
cannot conceal that conditional failure. Every pooled width is required to be
bit-identical to `fs_surrogate::conformal_band`; the production
`certify_or_escalate` seam must authorize finite calibrated bands and escalate
both unseen/unbounded and drift-invalid queries. A fixed in-support control
sequence must not fire the scalar PIT drift monitor, a fixed above-support
shift must fire within 64 observations, and recalibration with 19 shock
residuals must restore an admissible bucket. The caller mapping from a drift or
refusal verdict to `in_validity_domain = false` is exercised and recorded here;
it is not automatic library wiring.

One fs-rand draw chooses a cyclic regime/rank presentation rotation; the full
rank laws are still enumerated rather than sampled. Canonical configuration,
campaign, policy, replay, and result identities bind the residual lattices,
quantile rule, alpha, bucket geometry, RNG coordinates and semantics, drift
grid/order/budget, refusal and refresh policy, implementation versions, every
retained band and policy decision, exact coverage counts, and first mismatch.
An object-shaped `fs-obs` receipt is emitted before five statistical rows and
five linted, wire-validated `ConformanceCase` verdicts.

## No-claim boundaries
- Empirical-Bernstein/hedged closed-form CS for bounded means (mixture CS
  covers sub-Gaussian; EB variant is follow-up scope with its consumers).
- Conformal e-prediction (fs-surrogate + conformal-hardening beads).
- Both-ISA nightly rotation over a refreshed null/adversary corpus remains
  certify-certifiers follow-up scope. The landed rotation battery is one fixed
  finite seeded pseudorandom campaign: it can catch implementation drift but
  is not an exhaustive proof over null laws, stopping times, dependence
  structures, seeds, horizons, or ISAs. No cross-ISA execution evidence is
  produced by this code-first slice.
- The e-BH dependence battery exhausts its three declared finite laws, not all
  possible dependence structures or alternative configurations. It is an
  executable cross-examination of the production step-up rule and one seeded
  mutation, not a replacement for the arbitrary-dependence theorem. It makes
  no nightly-refresh, cross-ISA, scheduler/Cx, wall-clock, or throughput claim;
  central Cargo/RCH/DSR proof is pending for the code-first landing.
- The conformal-shift receipt is one exact two-bucket, no-tie, alpha-0.1 finite
  orbit and one scalar drift projection. It does not prove arbitrary
  distributions, ties, alpha/sample-size combinations, bucket schemes, shift
  magnitudes, weighted conformal validity, or automatic caller wiring. It
  makes no cross-ISA, Cx/cancellation, or performance claim; central batch
  proof is pending for the code-first landing.
- Two-sided betting tests; asymptotic variants; sub-exponential extensions.
- A caller-supplied finite span is a checked assumption, not a certificate that
  the underlying stochastic process has that support. Runtime breach refuses
  the observation; justifying the bound remains the caller's obligation.
- The canonical Gauntlet rows retain evidence for the fixed seeded empirical
  schedules above; they do not turn those observed rates or stopping times into
  an exhaustive finite-sample proof for untested data-generating processes.
- Under `docs/CERTIFICATE_REGIMES.md`, time-uniform and sequential evidence can
  support a declared long-horizon observable or reliability route only with
  the named population, sampling/dependence law, stopping policy, failure
  event, and separate model/field evidence. Anytime validity does not validate
  the simulator, mint population transfer, or certify one exact long chaotic
  trajectory.

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
