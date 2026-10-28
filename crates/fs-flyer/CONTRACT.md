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
- `spine::step` / `spine::advance` — deterministic six-DoF fixed-step
  integration. The array boundary is admitted into canonical `fs-ga::So3`
  before either load callback runs; body angular velocity advances through the
  `fs-time` body/right Lie step, and returned quaternions use the canonical
  representative of the rotation double cover. Zero-horizon advance and replay
  still perform body/state/timestep/group admission and return that canonical
  state without invoking the load callback.

## Canard hinge friction and stops

`canardmech::step` solves one monotone scalar midpoint-velocity equation
with a bounded safeguarded Newton/bisection iteration. A successful step
uses that same velocity for friction, stop damping, and applied work;
`StepReceipt` exposes all three. The quadratic stop spring uses its discrete
energy gradient across engagement/release, and stop damping acts over the
fraction of the configuration path beyond either stop. Mechanical energy
change plus both losses equals held-drive work to solve/roundoff tolerance.
The smooth Coulomb law still permits sub-Coulomb creep; exact static sticking
is not claimed. Non-finite derived steps and unresolved bounded solves refuse
with `mech-step-non-finite` and `mech-step-unresolved`, respectively.

## Invariants

- Admission precedes every computation; an unadmitted design cannot produce
  a build-up or panel.
- The declared empty mass must equal the component sum within 1 kg
  (`mass-spec-mismatch` otherwise — no hidden mass, ever).
- The hinge-axis admission domain IS the E1.5 prior [0.25, 0.50] x/c.
- Reference-config positions/gyrations are a CALIBRATED reconstruction
  (documented in-source): positions place the gross CG at the dossier's
  29.7% chord; gyrations land the inertias inside the published ±15% band.

## Spatial wind in the flight lifecycle

`SimLoop::init_with_spatial_wind(scenario, SpatialWindSpec)` opts into the
wall-compatible evolving `fs-atmo` field. The specification declares the mean
log law, spatial strength and length, advection speed, and every modal OU
correlation time; the scenario supplies the seed and load density. Its mean
reference speed must match the scenario headwind. This finite-mode spectrum
is synthetic and is not a validated reconstruction of historical weather.

The longitudinal force and propeller solves consume CG-sampled body air-relative
velocity. Ground velocity is separate inertial state: temporal OU changes and
spatial wind gradients change relative airspeed without fabricating an aircraft
impulse. Forces use the accepted tick-n field; the returned tick-n+1 state and
read-only `sample_spatial_wind` queries share the advanced field. The old scalar
OU angle increment is disabled in this mode. Gravity acts vertically in world
coordinates at the actual pitch attitude. The sampler exposes the full velocity
and analytic gradient in flight-line-forward/right/up coordinates, but the
current longitudinal force interface consumes only forward and vertical inflow;
distributed span/disk gust loads and crosswind/sideslip dynamics are not claimed.
Mechanical initialization retains the baseline steady trim, not gust-equilibrated
prelaunch dynamics.

Spatial runs save V3 checkpoints containing the entire wind specification,
both actual OU quadratures, ground velocity, and run identity. The existing
`restore_checkpoint` resumes this state exactly and checks the relative-air
velocity against the restored ground/wind state. Legacy initialization and V2
checkpoint bytes retain their existing behavior. Focused tests live in
`simloop/spatial_wind/tests.rs`: actual load/propeller response, zero-amplitude
gust exclusion, observer independence, rail/airborne replay, and gravity-frame
consistency. Run with `cargo test -p fs-flyer --lib simloop::spatial_wind`.

### Distributed strip and propeller inflow

`SimLoop::init_with_distributed_spatial_wind` samples every wing/canard strip
reference point and both propeller disk centers from the same frozen tick.
`OpenLoopDesign::force_buildup_with_incident` admits the corresponding explicit
body-relative incident field. Rigid pitch velocity `q cross (position - CG)`
enters each local sample once; the legacy pitch-as-twist approximation is
disabled on this path. Local chord-plane speed and incidence drive the nonlinear
strip solve, and its accepted local velocity drives Kutta-Joukowsky forces,
profile drag and section pitching moments. Converged propeller wash is included
in both circulation and force velocity. The separate disk-center axial inflows
drive the two BEMT thrust/torque solves and the existing engine dynamics.

`coupled_prop_airframe_step_with_incident` carries disk-induced velocity `vi`,
with `T = 2 rho A vi (V_disk + vi)`, and applies `wash_factor * vi` at each
strip. Returned strip velocities belong to the actual evaluated iterate that
produced circulation. Distributed acceptance checks each disk's relative
fixed-point residual separately; a strong rotor cannot mask the weak rotor's
error. The rationalized momentum root retains weak-loading accuracy.
Existing legacy and CG modes retain their previous wash
scaling and arithmetic. Invalid or reversed local flow and below-ground surface
samples refuse through the existing envelope receipt without a CG-flow fallback.

Distributed runs select V4 checkpoint integrity and a distinct run identity;
the complete OU state, wind specification and inertial velocity retain their
existing encoding. V3 explicitly selects CG loading and V2 selects legacy
loading. Changing a distributed payload to V3 cannot silently downgrade its
scenario identity. Focused coupling, force and lifecycle regressions cover
independent rotor inflow, the momentum/wash relation, canard-only gust pitch
loads, local-force reconstruction, single pitch-motion ownership and replay.

This remains quasi-steady strip aerodynamics with the existing spanwise prop-wash
overlap approximation, a common straight wake and
reference chord-row circulation ratios, axial disk-center BEMT, reference-flow
parasite/classical induced-drag estimates and reduced longitudinal motion.
It does not supply chordwise gust memory, crossflow dynamics, disk harmonics,
or a validation of historical turbulent flight.

## Error model

Typed `Refusal { code, message, ranked_repairs }`. Codes:
`non-finite-input`, `span-outside-domain`, `chord-outside-domain`,
`camber-outside-domain`, `area-inconsistent`, `hinge-axis-outside-prior`,
`pilot-mass-outside-domain`, `component-count-exceeded`,
`component-mass-invalid`, `mass-spec-mismatch`, `hinge-ratio-invalid`.
Caps tested at cap AND cap+1 (workspace law).

The six-DoF spine additionally reports `orientation-outside-group` before
invoking a load callback, `non-finite-loads` for a non-finite force or moment,
`non-finite-state` when finite admitted inputs overflow during a derived step,
and `lie-step-refused` if the canonical integration owner refuses a derived
step. Non-finite start/end time is `non-finite-input` and refuses before load
callback effects.

## Determinism class

Deterministic: pure arithmetic + fs-blake3 digests; the reference design
digest is pinned (golden-bump protocol; one bump recorded when travel_rad
moved from the 0.5236 literal to exact π/6). Spine digests are computed from
canonical output states, so quaternion representatives `q` and `-q` describe
the same callback inputs, trajectory, and digest.

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
