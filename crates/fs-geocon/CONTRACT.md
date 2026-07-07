# fs-geocon — CONTRACT

## Purpose and layer

L4 (ASCENT-facing geometry). Differentiable/certifiable geometric
constraint primitives (plan §7.6): the manufacturability and
design-intent constraints ASCENT enforces, each FIRST-CLASS — a value,
a derivative story, a certificate story, and its fs-constraint kind
mapping. Geometry supplies values and derivatives; ASCENT owns the
constraint SEMANTICS.

## Public types and semantics

- `min_thickness_soft` — the anti-paperclip constraint: mean p-norm
  aggregation of fs-query thickness samples. `soft_min` is the C¹
  optimizer value (an over-approximation converging DOWN to the
  minimum as p grows; exact on uniform samples, which keeps lever
  derivatives clean); `hard_min` and the LOCALIZED violation list are
  the ledger/verdict values. The two are reported side by side, never
  conflated. Oracle-skipped samples are counted.
- `draft_violations` — normals versus the pull direction for the mold
  half being assessed: smooth hinge² penalty (C¹) plus EXACT violating
  regions; normals opposing the pull within the half's own reach are
  UNDERCUTS, flagged separately; faces beyond `n·d < −0.5` belong to
  the other mold half (the v1 parting model is the plane perpendicular
  to the pull).
- `QuotientChart` + `SymmetryGroup` (`ReflectX`, `Cyclic{n}`,
  `Periodic{period}`) — symmetry ENFORCED BY CONSTRUCTION: evaluation
  points fold into the fundamental domain, so the shape is invariant
  for ARBITRARY inner designs and lever values; gradients chain back
  through the fold. Reflection folds bitwise; cyclic/periodic at fp
  scale. The seam makes the chart C0 (declared).
- `envelope_violation` — containment (`φ_allowed ≤ 0` on design
  boundary samples) and keep-outs (`flip = true`): the sampled worst
  plus a SUM-FORM log-sum-exp aggregate that is `≥ worst` within
  `ln(n)/β` — CONSERVATIVE, so driving the smooth value to zero drives
  the true worst to zero (never stops short).
- `volume_certified` / `volume_smooth` — over an EXPLICIT integration
  domain (fixed independently of levers, so derivatives see the shape
  change, not grid realignment): a RIGOROUS enclosure (sure-inside
  cells vs the Lipschitz uncertainty band; truth in `[lo, hi]` for
  1-Lipschitz-certified charts) beside a smoothed-Heaviside estimate
  whose lever derivative matches Hadamard on fixtures.
- `GeoPrimitive::descriptor()` — the declared table: differentiability
  class (fs-opt `Class`), certificate story (`Enclosure` /
  `SmoothEstimate` / `ExactByConstruction`), and the fs-constraint
  `ConstraintKind` mapping; `proof_escalation()` names the interval
  obligations where they exist.

## Invariants

1. Thickness: soft over-approximates converging down with p;
   violations localize to exactly the thin samples; the lever FD
   derivative matches the analytic 2 on the neck; a toy descent DRIVES
   the neck to feasibility (gcp-001).
2. Draft: a 10° wall passes 5° and fails 15° with all samples
   localized; vertical walls violate any positive draft; a mushroom
   shoulder is flagged as an UNDERCUT, not low draft; the smooth
   penalty's FD derivative matches the analytic hinge slope (gcp-002).
3. Symmetry: quotient shapes are invariant under their groups for
   ARBITRARY asymmetric inner designs (bitwise for reflection,
   fp-scale for cyclic/periodic) across random levers × groups ×
   orbits; folded gradients match finite differences off-seam —
   violation is structurally impossible (gcp-003).
4. Envelopes: containment/keep-out match analytic penetrations; the
   LSE aggregate is conservative within `ln(n)/β`; the FD derivative
   is right; the descent returns an escaping design fully inside
   (gcp-004).
5. Volume: certified enclosures bracket the analytic sphere volume at
   two resolutions and tighten with h; the smoothed volume's lever
   derivative matches Hadamard `4πr²` within 2%; a descent meets a
   volume cap (gcp-005).
6. The descriptor table is total, symmetry is `ExactByConstruction`,
   volume is `Enclosure`, thickness maps to Fabrication, and proof
   escalations are declared exactly where they exist (gcp-006).

## Error model

fs-query's teaching errors carry through (`NoGradient`,
`NotOnBoundary`, `Cancelled`, …); envelope assessment is total.
Honest gaps (skipped thickness samples) are COUNTED, not hidden.

## Determinism class

Fully deterministic: fixed sample sets, canonical iteration, no
randomness. Identical inputs give identical reports bitwise.

## Cancellation behavior

Volume integrations poll `cx.checkpoint()` per grid slab; thickness
aggregation polls through the fs-query oracle; both return the carried
`Cancelled` teaching error.

## Unsafe boundary

None. `#![forbid(unsafe_code)]` via workspace lints; no capsules.

## Feature flags

None.

## Conformance tests

`tests/conformance.rs`, cases gcp-001..gcp-006 — JSON-line verdicts,
seeded LCG randomness, fs-obs events for the volume/Hadamard table and
the descriptor table. Any reimplementation must pass the suite
unchanged.

## No-claim boundaries

- Derivatives here are FD demonstrations against analytic truths;
  exact adjoints join the gradient-stack bead.
- The draft parting model is the plane perpendicular to the pull;
  parting-line OPTIMIZATION and multi-pull molds are follow-ups.
- Symmetry groups are the plan's named trio about fixed axes/planes;
  arbitrary group generators and lattice groups follow with fs-ga.
- Envelope/volume rigor rests on the 1-Lipschitz chart contract;
  interval-arithmetic per-cell bounds (fs-ivl) would replace the
  Lipschitz band with tighter enclosures — the declared escalation.
- Mass is volume × constant density; density fields join fs-material.
