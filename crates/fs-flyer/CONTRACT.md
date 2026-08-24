# fs-flyer — CONTRACT

Bead: frankensim-wf-root-guzez.4.1 (E3.1, Wright Flyer program).
Spec: COMPREHENSIVE_PLAN_FOR_REAL_TIME_WRIGHT_FLYER_SIM_WITH_FRANKENSIM.md §5.1
(ROUND 6 steady state). Evidence: data/wright-flyer/{flyer-reference,
canard-mechanics,geometry-conventions,frame-conventions}.

## Purpose and layer

L4. The aircraft-level Wright Flyer model. E3.1 ships the `FlyerDesign`
schema (Round-2 canard/mechanics/structure fields), control-topology enums,
admission with typed refusals, the component mass/inertia build-up
(cross-checked against the published single-lineage inertias), and the
derived-quantity panel. The later E5 lifecycle engine now also lives here;
the browser boundary remains in `fs-flyer-wasm`.

## Public types and semantics

- `FlyerDesign { wing, canard, rudder, lateral, warp, masses… }` — frd-body-v1
  coordinates, +x forward from the wing leading edge, SI units, radians.
- `LateralControlTopology::{WarpWithSlavedRudder{ratio}, WarpIndependentRudder}`;
  `WarpStructureMode::{FlexibleTrussed, Rigidized}`.
- `MassComponent { mass, position, gyration_sq }` — six-bucket engine
  convention included; gyration radii carry each component's own spread.
- `mass_build_up()` → gross mass, gross CG, diagonal body inertia about the
  gross CG.
- `derived_panel(hinge_ratio)` → canard volume, naive two-surface neutral
  point, fixed AND free-control static margins, thin-plate hinge-moment
  gradient 2π·(x_h − ¼).
- `digest()` — canonical content identity under
  `org.frankensim.fs-flyer.design.v1` (a PhysicalScenarioId ingredient).
- `FlyerDesign::reference_1903()` — the dossier-valued reference config.
- `simloop::SimLoop` / `SimStateOut` — 120 Hz lifecycle state including
  longitudinal state plus Estimated reduced roll/yaw state. Snapshot v2 is a
  14-float layout whose first 12 words are the unchanged v1 prefix and whose
  appended words are roll attitude and heading.

## Invariants

- Admission precedes every computation; an unadmitted design cannot produce
  a build-up or panel.
- The declared empty mass must equal the component sum within 1 kg
  (`mass-spec-mismatch` otherwise — no hidden mass, ever).
- The hinge-axis admission domain IS the E1.5 prior [0.25, 0.50] x/c.
- Reference-config positions/gyrations are a CALIBRATED reconstruction
  (documented in-source): positions place the gross CG at the dossier's
  29.7% chord; gyrations land the inertias inside the published ±15% band.

## Error model

Typed `Refusal { code, message, ranked_repairs }`. Codes:
`non-finite-input`, `span-outside-domain`, `chord-outside-domain`,
`camber-outside-domain`, `area-inconsistent`, `hinge-axis-outside-prior`,
`pilot-mass-outside-domain`, `component-count-exceeded`,
`component-mass-invalid`, `mass-spec-mismatch`, `hinge-ratio-invalid`.
Caps tested at cap AND cap+1 (workspace law).

## Determinism class

Deterministic: pure arithmetic + fs-blake3 digests; the reference design
digest is pinned (golden-bump protocol; one bump recorded when travel_rad
moved from the 0.5236 literal to exact π/6).

## Cancellation behavior

Synchronous pure functions; nothing to cancel.

## Unsafe boundary

Workspace `deny(unsafe_code)`; no unsafe.

## Feature flags

None.

## Conformance tests

`tests/design_battery.rs`: mass/CG pinned to dossier values (gross 340.2 kg,
CG 29.7% chord); inertia reproduction inside the ±15% band vs the Jex-Culick
lineage with slender-biplane ordering; admission caps at cap AND cap+1
(pilot both edges, fabricated area, hinge prior, hidden-component
mass-spec falsifier, NaN/zero-span); derived panel pinned against
independent hand calculations (canard volume, naive NP, margin arithmetic,
zero gradient at x_h = ¼, self-driving sign, free-worse-than-fixed) plus
the documented naive-vs-Culick NP comparison; digest determinism +
sensitivity + pinned golden. JSONL receipts per case.

## No-claim boundaries

- The published inertias are a single-lineage reconstruction; agreement
  inside the band is a CROSS-CHECK of our decomposition, not a validation
  of either lineage.
- The derived panel's aerodynamic formulas are documented simple models
  (two-surface, no interference, thin plate). Their divergence from the
  Culick neutral point (3.9% c) is recorded data. Real stability claims are
  V-02a territory (E4.6a) and need the full aero model.
- The hinge-moment gradient is a SIGN/SHAPE device over the E1.5 prior; its
  quantitative level inherits the Estimated ceiling (A7a promotion path).
- The lateral state uses `ReducedAeroelasticWarp` and
  `ReducedLateralBuildUp`, with declared inertia and rudder-moment constants.
  It supports bank/heading presentation and reduced control-response claims;
  it is not a six-DOF structural-margin, post-stall, or calibrated-history
  claim.
- No aerodynamic force model lives here; sections are fs-airfoil, planform
  effects are fs-wing.
