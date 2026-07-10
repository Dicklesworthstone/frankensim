# fs-race CONTRACT

## Purpose and layer

Layer: L4 (ASCENT). e-RACING (plan §9.6, Bet 8 [M]): anytime-valid
sequential tests DRIVE structured candidate cancellation — pairwise
fs-eproc races with e-BH family-wise control eliminate dominated
candidates mid-evaluation, firing their fs-exec kill-handles. The [M]
payoff claim is measured, never assumed.

## Public types and semantics

- `race_field(loss, n, settings, kills)` → `RaceOutcome`: rounds are
  the ONLY clock — every survivor consumes exactly one observation per
  round in canonical index order, e-value crossings are evaluated only
  at round boundaries, so the elimination sequence is a pure function
  of (seed, logical stream identities), never wall-clock arrival.
  Full pairwise `PairwiseRace` matrix fed in BOTH directions; per-
  candidate elimination evidence = the strongest surviving opponent's
  log e-value; `e_benjamini_hochberg` at α across the surviving
  population per round; kills dispatched ascending (deterministic).
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

## Elimination-evidence validity (bead 7tv.7.1, the derivation)

Candidate i's elimination evidence is the MIXTURE (arithmetic mean,
computed in log space via `fs_eproc::combine_average`) of the pairwise
e-processes e_ji("j beats i") over the FIXED, predeclared family of all
n−1 ORIGINAL opponents. Validity: (1) each e_ji is a test
supermartingale for the pairwise null "i is not worse than j" (betting
process, predictable λ, bounded outcomes); (2) under candidate i's
composite null "i is not worse than ANY opponent" every family member's
null holds, so each e_ji has expectation ≤ 1 at every stopping time;
(3) a dead opponent's process is frozen at its elimination round — a
stopped supermartingale is a supermartingale (optional stopping), so
freezing preserves validity; (4) a convex combination of
supermartingales is a supermartingale, hence the mixture is itself an
anytime-valid e-process for i's null; (5) e-BH (Wang–Ramdas) controls
the elimination FDR at α under ARBITRARY dependence among the input
e-values, which round-by-round re-testing and adaptive elimination do
not violate because each input is anytime-valid. The REJECTED former
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

Non-finite losses are rejected STRUCTURALLY (fail closed): the
offending candidate is condemned in the round the value appears
(`invalid` + `eliminated`, kill fired), the poisoned value never
reaches e-processes or running means, and winner selection uses
`total_cmp` over means that are finite by construction — the former
panic at winner selection is gone.

## Error model

Structured panics on programmer contracts (field size, eta); loss
streams are caller-owned pure functions — non-finite losses are the
caller's contract violation and panic at the ordering comparison with
a teaching message.

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
successive halving.

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
