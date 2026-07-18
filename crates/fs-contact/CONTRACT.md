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

## Certified CCD (bead tqag, increment 2)

`certified_ccd(a, b, window, time_tolerance, max_windows, cx)` proves
clearance or localizes possible contact by conservative window
bisection over `CertifiedMotorTube::box_action_over` enclosures:

- SOUNDNESS (the Sev-0 no-tunneling claim): a subwindow is cleared only
  when the two whole-subwindow image enclosures are disjoint along a
  coordinate axis — no instant inside it can produce contact, with no
  sampling anywhere. Everything not proven clear subdivides to the time
  tolerance and is reported as a possible-contact window; the union of
  reported windows contains every true contact instant. ct-005 drives a
  bullet fully through a thin plate INSIDE the window (both endpoint
  enclosures disjoint — the exact trap endpoint sampling falls into)
  and requires a possible window containing the true crossing,
  localized to under 1% of the window.
- HONESTY: contact is never CLAIMED (box overlap is necessary, not
  sufficient); `ClearWindow` carries a certified lower bound on the
  axis gap; budget exhaustion is a refusal carrying the exact pending +
  unresolved windows in time order, never a truncated verdict (ct-008).
- THE ROOT-GUARD REFUSAL, EXECUTABLE (ct-007): bodies overlapping the
  whole window have no separation sign change for a global-root guard
  to find; certified CCD reports one possible window covering the whole
  domain instead of a false clear. This is why the design bisects
  enclosures rather than guarding roots of `separation(t)`.
- Determinism: LIFO bisection with the earlier half examined first;
  reports replay bit-identically (ct-006).

## Swept-vertex-hull refinement (bead tqag, increment 3)

`refine_possible_windows(a_vertices, a_tube, b_vertices, b_tube,
windows, max_iterations, cx)` re-tests each `PossibleContact` window
for POLYTOPE bodies: every body-frame vertex trajectory is enclosed by
`point_action_over`, and the convex hull of all trajectory-box corners
contains the body's image at every instant (a rigid image of a hull is
the hull of the vertex images). Certified separation of the two swept
hulls (`fs_query::convex_separation`; corner selection is the support
trait's documented exact case) PRUNES the window with a certified gap —
tight exactly where per-instant axis-aligned boxes are structurally
loose. ct-009 pins the trap: two 45°-rotated cubes passing on the
diagonal whose AABBs overlap at EVERY instant (the box route can never
clear them at any tolerance) prune with the analytic edge-to-edge gap.
Soundness: pruning requires `separation_proven` over a SUPERSET of each
swept body, so a window containing true contact can never be pruned
(ct-010 keeps the bullet's crossing window Retained); retention claims
nothing, exactly like the box verdict.

## SDF-obstacle route (bead tqag, SDF increment)

`refine_windows_against_sdf(vertices, tube, obstacle, windows,
time_tolerance, max_windows, cx)` prunes possible-contact windows for a
polytope body against a STATIC exact-distance chart. Soundness: exact
Euclidean signed distance is 1-Lipschitz (the theorem carried by
`TraceStepClaim::ExactDistance` — weaker claims refuse at entry by
capability name), so with the swept-vertex-hull corners enclosed in a
ball of radius `r` around any center `c`, every swept point satisfies
`φ(q) ≥ φ_lo(c) − r`; a positive bound proves the whole subwindow clear
with a certified gap. The radius computation self-covers its f64
rounding (multiplicative inflation, documented in-line); the center
choice affects only tightness. Because certified_ccd MERGES adjacent
possible windows, the route bisects internally with the same
LIFO/tolerance/budget discipline (budget exhaustion refuses with exact
partial state). ct-011 pins the value (a corner-region pass the AABB
route retains forever prunes against the curved surface; a lying chart
with a weaker claim refuses at entry) and ct-012 the soundness arm (a
through-shot's sphere-entry window survives as Retained).

## No-claim boundaries

- Certified CCD verdicts remain ENCLOSURE verdicts: `PossibleContact` /
  `Retained` windows localize in time but never adjudicate contact;
  time-of-impact enclosures tighter than window bounds, refinement for
  non-polytope MOVING bodies (needs swept support maps from fs-motion),
  and MOVING SDF obstacles (the chart is static in this increment) are
  later work in this bead's staging plan. Stage 2 consumes
  simulated-flow tubes through a tube-source-agnostic interface.
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
- The broad phase and CCD prune on certified geometry enclosures, but
  the motion versor defect is REPORTED (`BroadPhaseReport::max_defect`,
  `CcdReport::max_defect`), not folded into the boxes; the fold is
  still open in this bead's staging plan.
