# fs-race CONTRACT

## Purpose and layer

Layer: L4 (ASCENT). e-RACING (plan §9.6, Bet 8 [M]): anytime-valid
sequential tests DRIVE structured candidate cancellation — pairwise
fs-eproc races with e-BH false-discovery-rate control eliminate dominated
candidates mid-evaluation, firing their fs-exec kill-handles. The [M]
payoff claim is measured, never assumed.

## Public types and semantics

- `LossSpan` is a checked finite-positive maximum absolute paired-loss
  difference. It is fixed before the race, carried in `RaceSettings` and
  `RaceOutcome`, and therefore belongs in replay provenance.
- `RaceSettings::new(loss_span)` supplies the standard alpha and round budgets
  only after receiving an explicit scale. There is no `Default` implementation
  that can silently assume unit-normalized inputs.
- `race_field(loss, n, settings, kills)` → `Result<RaceOutcome, RaceError>`:
  rounds are
  the ONLY clock — every survivor consumes exactly one observation per
  round in canonical index order, e-value crossings are evaluated only
  at round boundaries, so the elimination sequence is a pure function
  of (seed, logical stream identities), never wall-clock arrival.
  Full pairwise `PairwiseRace` matrix fed in BOTH directions; per-
  candidate elimination evidence = the arithmetic mixture, computed in log
  space, over the fixed family of all original opponents;
  `e_benjamini_hochberg` at alpha across the surviving population per round;
  kills dispatched ascending (deterministic).
  `min_rounds` delays the first check (skipped, never peeked).
- `RaceOutcome`: survivors, elimination events `(round, candidate)`,
  winner (lowest running mean, index tie-break), evaluations used vs
  `fixed_n_equivalent`, and `savings()` — the falsifiable ledger.
- `successive_halving(...)` → `BracketLedger`: rank-based kills at
  budget milestones (standard SH semantics — does NOT carry the
  e-guarantee; documented), bracket schedule ledgered.
- Kill wiring: callers register candidate ids `0..n` in a
  `KillRegistry` to hold gates; eliminated candidates' whole
  evaluation trees drain at their next poll point.
- `RaceError` is a no-verdict outcome for malformed settings, non-finite
  losses, arithmetic overflow, or a paired difference outside `LossSpan`.
  Earlier-round kill requests cannot be rolled back if a later observation
  breaches support; callers must treat the whole returned race as no-claim.

## Invariants

1. Bitwise replay: identical inputs give identical elimination
   sequences, winners, and counters (race-001).
2. Ground truth: on a separated field the true best wins and every
   dominated candidate is eliminated within budget (race-002).
3. ANYTIME VALIDITY, empirically: across 200 seeded replays the true
   best was eliminated 0 times against an α = 0.05 budget of 10 ± 9.2
   (3σ binomial) — zero excess false elimination (race-003).
4. The MEASURED payoff: 11.7× evaluations saved vs fixed-N on the
   separated fixture (the stated 2–5× claim, exceeded and gated at
   ≥ 2×); the INSEPARABLE field reports 1.03× — no fake payoff — with
   elimination α-controlled (race-004).
5. Kill gates fire exactly for eliminated candidates; survivors'
   gates stay clean (race-005).
6. Successive halving follows its declared bracket schedule and beats
   fixed-N while the true best survives (race-006).
7. The statistical scale is structural: support boundaries map to exactly
   0/1, one ULP outside refuses without changing pairwise wealth, changing
   loss units and `LossSpan` together leaves decisions unchanged, and the
   equal-mean skew family that clipping misclassified stays calibrated
   (race-010/011 plus fs-eproc G0 tests).
8. Non-finite e-race observations abort without a verdict or a newly inferred
   candidate-specific kill; rank-only successive halving keeps its separately
   documented structural-invalid behavior (race-009).
9. Winner means use overflow-safe online updates, so finite values near
   `f64::MAX` retain their mathematical ordering (race-012).

## Elimination-evidence validity (bead 7tv.7.1, the derivation)

Candidate i's elimination evidence is the MIXTURE (arithmetic mean,
computed in log space via `fs_eproc::combine_average`) of the pairwise
e-processes e_ji("j beats i") over the FIXED, predeclared family of all
n−1 ORIGINAL opponents. Validity: (1) each e_ji is a test
supermartingale for the pairwise conditional null
`E[L_i,t - L_j,t | F_(t-1)] <= 0` ("i is not worse than j"; betting
process, predictable lambda, bounded outcomes); marginal equality alone is
insufficient; (2) under candidate i's
composite null "i is not worse than ANY opponent" every family member's
null holds, so each e_ji has expectation ≤ 1 at every stopping time;
(3) a dead opponent's process is frozen at its elimination round — a
stopped supermartingale is a supermartingale (optional stopping), so
freezing preserves validity; (4) a convex combination of
supermartingales is a supermartingale, hence the mixture is itself an
anytime-valid e-process for i's null; (5) e-BH (Wang–Ramdas) controls
the elimination FDR at alpha under ARBITRARY dependence among the input
e-values. At each batch, the live-family threshold is at least as strict as
the full-family threshold for the cumulative rejection count; rejected
processes stay frozen. Thus the cumulative set remains e-BH self-consistent
for the final vector of stopped e-values. Under the global null only, FDR is
the probability of any rejection; no general family-wise-error claim is made.
The REJECTED former
construction — the maximum over currently-surviving opponents — fails
(4) (a max of e-values is not an e-value) and additionally selected its
family from the same data (survivor-dependence); the battery's
certifier test (race-008) demonstrates the max exceeding the E[e] ≤ 1
Markov budget (measured 1.51 ± 0.27 vs the mixture's 0.79 ± 0.15 under
the exchangeable null) and race-007 pins the shipped construction's
any-elimination rate under the global null within α plus binomial
slack, with optional stopping and adaptive elimination active. The
certifier's first catch was the battery's own former noise fixture,
whose per-candidate persistent offsets made it a non-null (recorded in
the battery header).

Every pair uses `d = ((L_b - L_a) / LossSpan + 1) / 2` with no clipping.
Clipping is invalid because it changes the estimand: an equal-mean skew
distribution can have a nonzero clipped-difference mean. Non-finite losses,
subtraction overflow, and out-of-support differences abort with `RaceError`
and no `RaceOutcome`. The current round cannot emit elimination evidence.

## Error model

`race_field` returns structured `RaceError` values and no verdict for fewer
than two candidates, invalid alpha/round settings, non-finite losses,
subtraction overflow, and support breaches. A late failure can follow valid
earlier eliminations whose kill requests have already fired; returning an
error revokes the aggregate race claim but cannot undo external cancellation.
`successive_halving` retains explicit panics for invalid field size/eta and
for an all-invalid field because it is a separate rank-based primitive.

## Determinism class

Bit-deterministic by construction (see Public types): rounds as the
only clock. Parallel evaluation of a round cannot change the result
because crossings are checked only at round boundaries — the same
read-parallel/apply-canonical discipline as fs-mesh's coloring.

## Cancellation behavior

The crate IS the cancellation driver: eliminations request the
candidate's gate; everything running under that gate drains at its
next poll. The tournament loop itself is bounded and synchronous.

## Unsafe boundary

`#![deny(unsafe_code)]` via workspace lints; no capsules.

## Feature flags

None.

## Conformance tests

`tests/battery.rs`: race-001 replay; race-002 domination; race-003
false-elimination calibration (200 replays); race-004 measured
savings, separated and inseparable; race-005 kill wiring; race-006
successive halving; race-007 adaptive global-null calibration; race-008
mixture-vs-invalid-max e-value certifier; race-009 no-verdict non-finite
refusal; race-010 clipping counterexample and span guard; race-011 scale
covariance and malformed-setting refusal; race-012 overflow-safe means.

## No-claim boundaries

- Reclaim-LATENCY histograms (the ≤ 200 µs systems gate) need the
  real async tile-pool lanes under load — perf-CI scope; the smoke
  tier proves the wiring, not the latency.
- CMA-ES/NES/bayesopt integration APIs: the driver ships; optimizer
  glue lands with the ornithoid flagship's step 2.
- fs-uq CS-stopping cross-wiring (per-candidate MLMC streams stopping
  on their own confidence sequences): demonstrated independently in fs-frame's
  fragility stage; the joint API here is a successor.
- Elimination-order OPTIMALITY (racing theory regret bounds): the
  battery gates validity and measured savings, not minimax rates.
- `LossSpan` proves only that the supplied number is finite and positive.
  Establishing that it bounds the process almost surely is the caller's
  scientific obligation; runtime checks catch observed breaches, not future
  or unobserved tail mass.
- No transactional rollback of already-dispatched kill requests after a later
  support failure. Production orchestration must ledger the returned error and
  discard the aggregate tournament claim.
