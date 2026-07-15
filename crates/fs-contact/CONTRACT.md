# CONTRACT: fs-contact

> Status: ACTIVE (Stage 1, increment 1 — bead tqag). Capability-routed
> body-to-body contact over certified motion.

## Purpose and layer

Blocker B3 (expansion plan, phase E2): body-to-body contact detection
with certificates instead of sampled heuristics. Layer: **L3** (deps:
fs-motion L2, fs-query L2, fs-geom L2, fs-ivl L1, fs-exec L0).
Explicitly NOT a dependency of-or-on fs-solid/fs-mbd solver internals —
those consume adapters; reusable contact protocols live here.

## Public types and semantics

- `SpacetimeBody` — a finite body-frame support box bound to a
  `fs_motion::CertifiedMotorTube` (body-to-world). Validation refuses
  non-finite/inverted supports.
- `spacetime_candidates(bodies, window, max_pairs, cx)` →
  `BroadPhaseReport`: the conservative spacetime broad phase. Each
  body's windowed world box is `CertifiedMotorTube::box_action_over`
  over the WHOLE window — an enclosure for every `t`, so a
  non-overlapping pair provably cannot touch inside the window (no
  sampled instants, no tunneling between samples). Deterministic
  sweep-and-prune on world `x` (`total_cmp`, index tie-breaks); output
  pairs sorted; report carries checked/pruned counts and the worst
  motion versor-defect bound, which consumers must carry forward.
- `NarrowRoute` / `narrow_phase(pair, route_a, route_b, iters, cx)` →
  `NarrowVerdict`: capability routing. Stage 1 routes Convex×Convex
  through fs-query's certified `convex_separation` (its semantics pass
  through unchanged: `separation_proven ⇔ lo > 0`, overlap never
  claimed). Any pairing without a compatible declared route refuses
  with `MissingCapability` naming the pair and capability — never a
  guess.
- `ContactError` — typed refusals throughout;
  `CandidateBudgetExhausted` (program risk #2) lists every unresolved
  overlapping pair beyond the budget so the resolved prefix is never
  mistaken for the complete candidate set.

## Invariants

- Broad-phase candidacy is conservative over the query window: a pair
  absent from `pairs` has certifiably disjoint windowed enclosures.
- Output ordering is a pure function of the inputs (deterministic
  sort keys everywhere; no HashMap iteration).
- Refusals leave no partial claim: budget exhaustion returns the
  unresolved remainder, capability gaps name the pair.

## Error model

`ContactError` wraps `fs_motion::MotionError` and
`fs_query::QueryError` unchanged (their teaching text passes through)
and adds contact-specific refusals: body-count/support/window
validation, candidate budget exhaustion with the unresolved list,
missing narrow-phase capability, cancellation.

## Determinism class

Bit-deterministic given deterministic inputs: sorted sweeps, fixed
tie-breaks, fs-motion/fs-query deterministic enclosures underneath.

## Cancellation behavior

`Cx` checkpoints per body enclosure and per sweep row; narrow phase
inherits fs-query's cancellation strides.

## Unsafe boundary

None. Workspace lints; no `unsafe` blocks.

## Feature flags

None yet. CCD lanes will gate under features when they land.

## Conformance tests

`tests/contact.rs`, cases ct-001..ct-004: analytic screw-motion broad
phase (approach window overlaps, retreat window disjoint, both against
hand-computed enclosure geometry); determinism replay; budget
exhaustion listing exact unresolved pairs; capability refusal; convex
narrow-phase distance containment at a frozen time against the
analytic value.

## No-claim boundaries

- CERTIFIED CCD is NOT yet claimed: no time-of-impact enclosures, no
  conservative advancement, no feature-pair spacetime inclusion.
  Increment 2 lands them over prescribed `CertifiedMotorTube`s;
  Stage 2 consumes simulated-flow tubes through a tube-source-agnostic
  interface. Per the bead's load-bearing subtlety, CCD will be
  feature-pair-wise — a global `separation(t)` root guard is refused
  as unsound during persistent contact, and a fixture will prove that
  refusal.
- Narrow-phase routes: Stage 1 is Convex×Convex only. SDF-pair local
  gaps (fs-query `ImplicitGapOracle`), nonconvex decomposition,
  interval global optimization, and mixed-route pairings all refuse
  as `MissingCapability` today.
- Penetration depth is never claimed (fs-query's convex overlap
  no-claim passes through); EPA-class certificates arrive with
  fs-query bead hk8f5.
- Rep Router conversion/motion errors do not yet inflate contact
  bounds (fs-query bead fugfk); claims apply to the presented charts,
  not to abstract regions behind conversions.
- The broad phase prunes on certified geometry enclosures, but the
  motion versor defect is REPORTED, not folded into the boxes; the
  fold lands with the CCD increment.
