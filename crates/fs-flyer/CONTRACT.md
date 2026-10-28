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
  The second load evaluation receives a first-order endpoint velocity predictor
  including both halves of the initial moment impulse; smooth angular drag and
  forces coupled to angular velocity therefore retain second-order convergence.
  The accepted rotational update delivers the trapezoidal world angular impulse
  `dt/2 (R_start tau_start + R_end tau_end)` to roundoff. This is a finite-step
  balance, not an exact solution for general varying loads. The spine battery
  checks exponential angular decay, its integrated angle and coupled translation,
  and a constant world-frame torque on an asymmetric rotating body (G0/G1).

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

Ordinary, CG spatial-wind and distributed-wind airborne paths share one coupled
longitudinal/rotor midpoint step. The first native aerodynamic solve predicts
midpoint velocity, pitch attitude/rate and rotor speed. A second solve evaluates
the actual loads at that predicted state, including spatial wind at the stage's
CG and material points. Both stages use the same held controls and temporal wind
realization; the legacy scalar gust's incidence increment uses stage airspeed.
The first solution warm-starts the second, and only the accepted midpoint wash
is cached for the next tick. There is no duplicate aerodynamic owner or extra
atmosphere/control clock advance.

The force point uses `X_mid = X_n + dt V_n / 2`; accepted centre-of-mass drift
uses the mean of initial and final ground velocities. This closes both linear
impulse and orbital angular impulse about any translated origin to roundoff:
`Δ(X × mV) = dt X_mid × F_mid`. Pitch and rotor inertia receive their respective
midpoint torques. Smooth coupled longitudinal/rotor dynamics are second order
under the held inputs; this does not extend the claim to OU innovations,
discrete controls, existing load-envelope clamps or the reduced lateral model.
Constant force/gravity and constant angular accelerations integrate exactly.
A refused midpoint publishes no predicted motion or warm cache. A positive
midpoint CG reaches the native force solve even when the held-load endpoint
predictor would cross ground. Actual accepted midpoint kinematics determine
GroundContact; the terminal state retains that midpoint velocity, pitch/rotor
impulse and warm cache, with height placed on the wall. Distributed material
points remain subject to native atmosphere admission: a positive CG does not
certify that every surface is above ground. Only a nonpositive midpoint CG uses
the bounded held-load terminal fallback without requesting loads underground.
These event branches claim no smooth order or continuous collision/localization.

The ordinary path retains body airspeed as stored state, converting through the
steady headwind offset before and after the inertial step. Force-free attitude
changes therefore preserve world momentum and kinetic energy to roundoff instead
of creating energy through explicit rotating-frame transport. Airborne ground
motion may be backward. Both paths remove the force assembly's level-flight
gravity convention and apply vertical world gravity once, independent of angle
of attack, pitch attitude or a gust. Focused analytic regressions live in
`simloop/motion.rs`. Ordinary and spatial launches also share one constrained
world-momentum step: stationary pitched aircraft have both body airspeed
components, the loaded rail holds vertical ground velocity at zero, and actual
rotated non-gravity loads determine support, friction and release. If tick-start
upward load meets or exceeds weight, the same full airborne midpoint interval
begins immediately, including vertical/pitch/rotor and lateral evolution, with
no preceding rail friction or tow impulse. Release introduces no velocity jump
and advances controls/atmosphere only once. Positively loaded rail ticks retain
the existing single-force-solve constraint/tow law; within-tick loss of support
is not localized.

Pilot member 3 retains its registered gains. Its earlier crash-class calibration
predates these momentum/gravity corrections and is no longer a validation claim
for the current plant. With consistent stationary/rail velocity, the unchanged
seed-1903 Dec-17 reference reaches RailEndWithoutLift at tick 975, approximately
18.316 m: the 18.3 m rail is too short for this currently calibrated plant.
The battery pins that outcome rather than fabricating lift or retuning gains.
Airborne lifecycle tests explicitly add the existing registered catapult as a
test fixture, distinct from historical Dec-17. With same-tick support release,
member 3 flies 11.41 s with five undulations and 44.90 m of travel before the
native outer coupling solver refuses midpoint residual growth at tick 1443.
Fixed controls fly 5.02 s, reach 11.60 m forward and turn back to 4.09 m before
the propeller station solver exhausts its 120-evaluation budget at tick 676.
Stable station momentum arithmetic and local residual-certified finishing solve
the earlier tick-673 slow-convergence failure; the later operating point remains
unresolved within that bounded local solve. It is not a station wake-domain
refusal, and no broader station-robustness claim is made. These are
numerical nonconvergence receipts carried by EnvelopeExceeded, not established
physical envelope boundaries. The substantive flight duration/undulation and
travel gates remain, along with exact terminal code/tick and replay evidence.
Neither assisted fixture nor original rail shortfall establishes a revised
historical launch or pilot calibration.

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
velocity against the restored ground/wind state. Legacy initialization and the V2
checkpoint layout retain their existing encoding. Focused tests live in
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

### Recording and scrubbing the selected wind model

`replayenv::{record_replay_with_spatial_wind,
record_replay_with_distributed_spatial_wind}` retain the exact V3/V4 tick-zero
checkpoint in the existing envelope layout. Legacy `record_replay` keeps its
checkpoint cadence and serialized layout. `scrub_to_tick(0)` reads the actual
initial state without advancing either physics or the wind; positive ticks
resume the nearest strictly earlier checkpoint and reproduce that step's output.
Both tick-zero digest and run intent must match. Spatial replay requires its
original initial checkpoint even when a later checkpoint is available. Missing,
downgraded, duplicate or wrongly indexed initialization refuses instead of
selecting legacy gusts. The existing deterministic-pilot recording restriction
remains; this does not add an applied human-input trace. The already emitted
damage-model-unavailable terminal code also round-trips through the parser.

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
