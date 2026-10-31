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
- `mass_build_up()` → gross mass, gross CG and full body inertia about the
  gross CG. `cross_inertia_kgm2` retains the negative parallel-axis products;
  `inertia_tensor()` reconstructs the symmetric matrix and `rigid_body()`
  prepares the original six-DOF body directly from it.
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

`spine::RigidBody` owns the shared, immutable `InertiaTensor3`. `new(mass,
principal_moments)` explicitly declares principal body axes;
`from_inertia_tensor(mass, tensor)` retains arbitrary products of inertia.
Principal-frame preparation occurs once, outside all tick loops. Ordinary
flight torque kicks use the full inverse tensor, and native free rotation
transports angular momentum through the same fixed body-frame transformation.
Replay and partitioned flight carry that same body; no geometry or load axes
are relabelled. The diagonal path retains its original arithmetic and digest.
G1/G3 fixtures independently check component raw moments, the actual
design-to-flight path, body-coordinate covariance and cross-axis motor recoil.
Component gyration inputs still describe body-aligned diagonal COM tensors;
this does not infer unprovided intrinsic cross terms or changing distributions.

## Finite-strip rotational fluid inertia

`addedmass::assemble_analytic_strip` integrates each finite uniform strip's
normal translation, pitch about its explicit span axis, and normal-velocity
variation along its span. With chord `c`, span `b` and normal added mass
`m_a = rho*pi*c²*b/4`, its centroidal rotary terms are `m_a*c²/32` about the
span axis and `m_a*b²/12` about `span_axis cross normal`. The pitch coefficient
is the flat-plate limit of the ellipse coefficient in
[MIT 13.021, supplemental problem A1, page 3](https://web.mit.edu/13.021/13021_2004/suppnotes/supp04.pdf),
with ellipse semiaxis `a=c/2`; the span coefficient is the exact uniform second
moment. These are coefficients of the existing independent-strip potential-flow
approximation, not a full three-dimensional rectangular-plate solution.

`Strip` supplies orthonormal `normal`/`span_axis` directions plus the existing
centroid `control_gain` and explicit `control_angular_gain`. All three positive
kinetic-energy modes contribute to rigid, hinge and cross blocks together;
their shared added momentum supplies the original workless gyroscopic bias.
The joint solve also retains the rigid body's complete inertia tensor.
The actual extraction path uses these terms and records algorithm version
`config-resolved-strip-v2`; its existing wire schema is unchanged.

G0/G2/G3 tests in `addedmass_battery` cover the thin-ellipse pitch coefficient,
hinge/rigid recoil, finite-span subdivision, reference shifts, oblique frame
rotation, length-to-the-fifth rotary scaling, PSD, work closure and invalid
orientation. The intentionally changed joint and extraction goldens distinguish
the added physical storage. Geometry/control Jacobians remain frozen during
each assembly; configuration derivatives for finite hinge motion, cross-strip
fluid interaction, free-surface effects and finite-aspect-ratio corrections are
outside this approximation.

## Hydrostatic pressure and material atmosphere loads

`atmosphere::MaterialAirProfile` composes the L2 atmospheric wind with the
existing L3 material gas owner at L4. It either retains one coherent arbitrary
gas snapshot or uses `StandardAtmosphere1976` with an explicit geometric ground
elevation. Wind heights remain above-ground distances; their sum with ground
elevation is the profile's geometric sea-level altitude. Queries preserve the
source's 0..80 km domain and retain its complete pressure, temperature, density,
sound speed and transport state. No duplicate EOS or hydrostatic law is added.

Frozen and evolving samples use the original mean/turbulence flow. They
transform both vector and gradient into flight-line forward/right/up, without
advancing an OU clock or consuming randomness. Legacy `AirScenario` constants
cannot override the material state. Homogeneous gas and the existing SimLoop
initialization/checkpoint behavior remain unchanged.

`OpenLoopDesign::force_buildup_in_atmosphere` and its evolving counterpart
subtract the sampled local wind from explicit inertial ground velocity, use
the shared longitudinal frame transform, and send this sample's density into
the original wing, drag and propeller load solve. G0/G3 tests compare the full
native result, verify hydrostatic `dp/dz = -rho*g`, pressure/density scaling of
dynamic pressure and Reynolds number, real lift/thrust scaling with altitude,
and exact read-only wind sampling. No extra constant reference density appears.

This is prescribed CG-sampled wind and a gas state frozen across the aircraft
for one load evaluation. The standard atmosphere is a dry-air reference, not
local-weather evidence. The original reduced longitudinal loads still omit
sideslip and Reynolds/Mach-dependent polar corrections; their gravity/frame
convention is unchanged. The bridge does not model pressure work, buoyancy,
weather evolution or validate high-altitude aircraft performance, and does not
make the historical SimLoop stratified. The prescribed solenoidal velocity
field does not imply divergence-free mass flux in a stratified gas; this is
not an anelastic or compressible atmospheric-flow solve.

## Six-DOF rotor momentum and gyroscopic motion

`spine::rotors::{RotorAssembly, RotorState, RotorTorques, step_held}` couples
up to 64 independent axisymmetric rotors to the native six-DOF spine.
`reduced_body` includes the whole assembly's mass and transverse/offset rotor
inertia but excludes rotor axial spin inertia. Its full tensor and all rotor
axes use the same declared body frame, which need not be principal. Signed absolute axial momenta
are `p = J (relative_spin + axis dot omega)`; `state_from_relative_spin`
constructs them from physical relative speeds. This is the reaction-wheel
reduction documented by Basilisk's `mrpSteering`, equations (2)--(6).

Motor torque changes rotor momentum and produces equal/opposite carrier recoil
without becoming external torque. Rotor aerodynamic torque changes both rotor
momentum and total external torque. `loads.moment_nm` excludes those rotor
aerodynamic torques, which the owner adds once. The native `fs-time` gyrostat
map includes precession from all stored rotor momenta. Translation uses a
world-force drift/kick/drift update, preserving its discrete orbital impulse.
The endpoint state includes unwrapped relative rotor phases and absolute
momenta; the receipt returns relative speeds, fixed-origin linear/angular
impulse residuals, separate motor/external work, actual kinetic energies and
the numerical work defect. No numerical defect is relabeled as dissipation.

The successful path is allocation-free, O(rotors) plus at most five elementary
gyrostat rotations, without a nonlinear solve. State, inertia, axes, loads,
count and finite derived quantities are admitted before publication. G0/G1/G3
tests in `tests/rotor_spine.rs` cover analytic motor recoil/phase/work,
external aerodynamic braking, multiple-axis impulse balance, balanced
counter-rotation and zero-rotor reduction, independent locked-inertia energy,
gyroscopic response, second-order phase/work refinement, replay and refusal.
The existing spine battery stays unchanged.

This is second order for the declared held loads, fixed axes and constant
inertia; coupled state-dependent loads need appropriate stage evaluation.
Physical rotor inertia must not be inferred from reflected engine/gearbox
shaft inertia. Gears, gimbals, imbalance, flexible rotors, changing mass and
contact are outside this owner. No historical calibration, complete rotor
checkpoint wire format or integration into the reduced longitudinal SimLoop
is claimed by this additional six-DOF API.

### Ideal gears and real propeller/airframe motion

`spine::rotors::geared` admits one ideal relative gear coordinate and 1–64
axisymmetric rotors with explicit signed ratios `g_a`, physical axial inertias
`J_a`, and body-fixed unit axes `e_a`. The caller supplies the whole assembly's
mass and full locked inertia tensor, including all rotor inertia. Assembly
computes `J_eq = sum(J_a g_a²)`, `c = sum(J_a g_a e_a)`, and the SPD Schur
reduction `I_red = I_locked - c cᵀ/J_eq`. Products of inertia are retained by
`fs_time::lie::gyrostat::InertiaTensor3`; its principal-axis assembly occurs once.
State is common relative phase `s`, canonical momentum
`p = J_eq s_dot + c.dot(omega)`, and the full carrier six-DOF state.
Total body momentum is `I_red omega + c p/J_eq`. Physical kinetic energy is
`m |V|²/2 + omega.dot(I_red omega)/2 + p²/(2 J_eq)`; the battery independently
reconstructs it from the original locked inertia and relative rotor rates.

Internal generalized motor torque changes `p`; aerodynamic torque changes it
by `sum(g_a tau_a)` and contributes `sum(e_a tau_a)` to total external torque
exactly once. Carrier moments must exclude those axial shaft torques. Motor
work is `Q_motor delta_s`; aerodynamic shaft work uses absolute axial rate
`g_a s_dot + e_a.dot(omega)`. A perfectly balanced counterrotating pair has zero
`c` and no motor recoil; unequal external shaft torque still changes total
angular momentum. Signed ratios are never replaced by reflected inertia in
that cancellation law.

`geared::step_held` is second order for fixed world force/body torques;
`geared::step` evaluates a smooth load law at the initial state and a half-step
predictor before accepting a full step with the midpoint loads. The two calls
occur in deterministic start/midpoint order. The work residual is the actual
endpoint kinetic change minus declared work, with no artificial dissipation.
The global angular impulse combines the translation kick's orbital impulse
with trapezoidal endpoint-frame torque impulses; its residual measures the
accepted discrete momentum balance. Successful mechanics use fixed-size
storage. Invalid models/states refuse; input state is immutable, while side
effects inside caller-owned callbacks remain the caller's responsibility.

`spine::rotors::propulsion::GearedPropulsion` is a concrete consumer: both
stages run the existing distributed two-way wing/BEMT solve with actual stage
attitude, velocity, all three angular rates and generalized shaft momentum.
World-position material-wind samples are converted to body relative incident
velocities at every strip and both disk centers. The converged circulations
and accepted local wash drive full vector Kutta–Joukowsky forces, profile drag,
section moments, existing classical induced drag and parasite ledger, and
thrust at both disk arms. Gravity is added once in world axes. Force moments
retain all three components about the same assembly COM. Propeller mounts can
refer to two rotors inside a larger gear train, so additional engine/flywheel
inertia is supplied physically by the caller rather than invented.

Gear ratio, rotor axis, blade handedness and positive thrust direction are
distinct inputs. Handedness maps signed ABSOLUTE spin into the current BEMT
positive-speed branch; the physical aerodynamic torque reverses with that
handedness. The current actuator-disk coupling only supports body-x shafts and
explicit +x thrust. Stopped/opposite blade rotation returns
`propeller-spin-outside-domain`; reversed net disk flow retains the existing
`disk-inflow-invalid` refusal. Neither spin nor inflow is floored. Assembly
and stage defects use `propulsion-assembly-invalid`, `propulsion-state-invalid`
or existing native refusals. No unconverged/one-way aerodynamic fallback runs.

G0/G1/G3 evidence in `tests/geared_propulsion.rs` covers nonunit gear-ratio
recoil/phase/work, balanced and unequal torque cases, arbitrary-axis dense
inertia and independent physical momentum reconstruction, admitted-count and
state boundaries, current symmetric aerodynamic parity, real asymmetric
six-DOF motion, absolute-speed/handedness refusal, and second-order refinement
of both an analytic motor/drag law and real coupled aerodynamic dynamics.
The `GearedPropulsion` aerodynamic tier remains quasi-steady, with existing 2-D section closures,
axial disk inflow, straight wakes, spanwise wash overlap and Estimated drag.
It does not add rudder dynamics, propeller crossflow/harmonics, apparent mass,
shaft compliance/failure, gimbals, contact, checkpoint serialization, or
historical calibration. Historical preset mass/inertia and SimLoop behavior
are unchanged. Full coupled accuracy remains bounded by the declared wing,
BEMT and coupling tolerances; refinement evidence does not certify physical
validity or arbitrary operating branches.

### Spatially mounted propellers

`propcoupling::spatial::{DiskAxes, coupled_step}` and
`GearedPropulsion::new_spatial` extend that consumer to canted, body-fixed
shafts. Positive-thrust axes are admitted unit vectors, collinear with their
physical rotor axes; neither gear sign nor blade handedness chooses the
thrust direction. Native wing-induced velocity and the actual local rigid-body
incident velocity are projected onto each disk's axis before BEMT. The returned
thrust acts along that axis at its actual arm, while aerodynamic shaft torque
acts along the separately signed spin axis, exactly once in total momentum.

Spatial wash measures radius in the full disk plane, including vertical
separation. Each disk's scalar induced velocity becomes a vector increment
along its thrust axis at every strip; circulation and full-vector force
assembly use the same accepted local velocities. Geometry factors are prepared
once per coupling solve. The existing Aitken solver, per-disk convergence gate,
momentum root and refusal behavior are shared. The spatial recipe and actual
axes have a distinct content digest; `new` and the historical axial APIs keep
their original spanwise model and numerical behavior.

G0/G1/G3 checks in `tests/spatial_propulsion.rs` independently reconstruct
disk-plane overlap, projected native BEMT thrust/torque, vector wash and global
angular momentum; verify rotation covariance, the coplanar axial limit, and
invariance to reversing a mechanical axis and its coordinate signs; and reject
invalid mounts or a late failed stage without returning partial motion. This
remains a disk-average, quasi-steady axial model. Transverse disk loads,
azimuthal harmonics, swirl, contracted/skewed wakes, direct disk-to-disk wake
interference and gimbal motion are absent. The finite 2R wash window remains
an explicit model truncation, not a resolved far wake. Canted-mount support
does not establish a crossflow validity envelope or historical calibration.

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
18.314 m: the 18.3 m rail is too short for this currently calibrated plant.
The battery pins that outcome rather than fabricating lift or retuning gains.
Airborne lifecycle tests explicitly add the existing registered catapult as a
test fixture, distinct from historical Dec-17. With consistent disk-induced
wash and same-tick support release, member 3 flies 11.42 s with five undulations
and 44.90 m of travel before the native wing solver exhausts its 600-evaluation
budget at tick 1444 (residual 5.92825e-10). Fixed controls fly 4.84 s, reach
11.61 m forward and turn back to 5.41 m before the outer coupling exhausts its
12-correction budget at tick 655 (residual 0.00584813, tolerance 0.005).
The native strip-cycle safeguard preserves the unchanged Huffman 4 m/s catapult
liftoff gate: lift occurs at tick 197 and 10.229 m, followed by 4.66 s airborne
before an outer-coupling second-growth refusal at tick 756. The unchanged
3 m/s case passes its initial wing solve but reaches the 120-evaluation
propeller station cap at tick 8. These are numerical nonconvergence receipts
carried by EnvelopeExceeded, not established physical envelope boundaries.
The substantive flight duration/undulation and travel gates remain, along with
exact terminal code/tick and replay evidence. Neither assisted fixture nor
original rail shortfall establishes a revised historical launch or pilot calibration.

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

Both native propeller/airframe paths carry the same disk-induced axial velocity
`vi`. Ideal actuator-disk momentum gives mass flux `rho A (V + vi)`, far-wake
increment `2 vi`, and thrust `2 rho A vi (V + vi)`. The rationalized disk root
therefore returns `vi`, not `2 vi`. Each strip receives `wash_factor * vi` once;
the factor already supplies the disk-plane value 1 and the upstream/downstream
spatial variation. The former additional factor 1/2 in ordinary loading was
inconsistent with its own momentum root. Both loading modes now use the same
wash sum and root. The coupling spec digest uses v2 to distinguish the corrected
physical convention; candidate tuples, acceptance tolerances and budgets are
unchanged. The independent native regression reconstructs mass-flux thrust and
the washed wing circulation, and rejects a second halving of the wash.

This corrects the induction amplitude within the existing model: the ordinary
path retains freestream Kutta–Joukowsky force evaluation, while the distributed
path evaluates forces with local strip velocities. The finite 2R influence
window, spanwise overlap approximation and omission of wake contraction, swirl
and vertical disk overlap remain explicit model limitations. References: NASA
Glenn, [Propeller Thrust](https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/propeller-thrust/),
and Deters, Ananda and Selig, AIAA 2015-2265, equations (3)–(7).

Both native propeller/airframe paths judge residual growth against the evaluated
predecessor of the current correction. A first growth above 25% rejects the
entire trial and retries `x_previous + omega_min * residual_previous`; the
rejected state never becomes the retry's origin or Aitken history. Subsequent
corrections retain the relaxation floor, and a second genuine growth refuses.
Both the rejected trial and retry count against the unchanged candidate cap;
there are at most `cap + 2` native map evaluations including initialization.
The same acceptance tolerance remains in force. An admitted off-trim native
fixture checks rollback at cap 5 and exhaustion at the registered cap 12; the
old retry falsely classified its recovery as a second growth. Correct disk
wash also resolves the former member-3 tick-1443 force input within the unchanged
candidate-B tolerance; that local recovery does not establish a complete
flight envelope. No one-way fallback is introduced.

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

Material-bound spatial constructors (`init_with_material_spatial_wind` and
`init_with_material_distributed_spatial_wind`) admit the existing complete
`fs-material::GasState` and replace scalar scenario density **before** trim.
The original trim, wing/canard loads, propellers, rotor dynamics and lateral
loads consequently share its density. `scenario()` exposes this canonical
scenario; `sample_material_air` returns the same accepted wind and gas, with
dynamic-pressure and Reynolds diagnostics relative to an explicit observer
velocity. The complete thermodynamic and transport tuple remains available
without re-evaluation or a second reference-air assumption.

Material runs use V5 checkpoints, retaining the whole numerical gas and an
explicit CG/distributed sampling mode. Every coefficient binds tick-zero and
run identity, including properties not consumed by the current force laws.
The material replay recorders preserve the canonical density and original
initial checkpoint. A contradictory density at restore refuses. Historical
V2/V3/V4 bytes are unchanged. Focused G0/G3/G5 tests in
`simloop/spatial_wind/material/tests.rs` cover actual pressure-dependent flight,
native force/trim equivalence, q/Re, observational noninterference and replay.

The material gas is homogeneous and prescribed: no heat/composition evolution,
altitude stratification or fluid/structure thermodynamics is introduced. The
reference polars remain Reynolds-independent and incompressible; retaining
viscosity, temperature and sound speed does not claim those corrections.
Numerical checkpoints do not serialize material-source evidence or qualify a
gas composition for historical air; that authority stays with the producing
model. The existing flight density envelope still applies.

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


## Transient wake, shaft and aircraft propulsion

`spine::rotors::propulsion::dynamic::DynamicPropulsion` couples two original
uniform axial inflow states to the existing physical geared rotor/aircraft
assembly. `DynamicState` retains both disk-induced velocities alongside native
six-DOF motion and canonical shaft momentum. The apparent-air-mass model is
owned by `fs-airscrew::dynamic_inflow`; its density remains fixed throughout
a trajectory and its air mass is excluded from aircraft inertia.

Each actual aerodynamic evaluation solves the original wing under the
current vector wash, projects incident plus wing-induced velocities onto
the declared disk axes, and evaluates original blade forces at the stored
induction. The shared carrier-wrench owner applies thrust arms, counterrotating
shaft loads, internal motor torque, and fixed/canted geometry exactly once.
The quasi-steady and dynamic models reuse the same section forces and
body/shaft wrench assembly while retaining distinct public result types.

The initial real loads predict aircraft, shaft and both wake states to
half time. A second real evaluation supplies the accepted native gear map
and the two implicit-midpoint wake-flux steps. The complete result is published
only after every motion and wake stage succeeds; borrowed input state is
unchanged on late refusal. Callback side effects remain caller-owned.
World aircraft/rotor impulse balances and axial reduced-wake impulse/work
balances remain separately inspectable, with no resolved-wake world balance
or whole-fluid/body energy theorem inferred from them.

Four new native coupled cases establish canted startup momentum accounting,
state-dependent motor/gust/body/shaft/wake refinement, exact shaft-sign
convention invariance, and late density refusal. Together with six dynamic
inflow and 37 original BEMT/station/geared/spatial/coupling cases, 47 focused
native cases pass and legacy goldens remain unchanged. Smooth admitted
coupled refinement is second order. The tier retains normal axial,
constant-density, disk-average flow, quasi-steady blade/wing sections and
the existing wash window; it supplies no skew/crossflow loading, swirl,
harmonics, transported wake geometry, blade memory, gimbal motion, historical
calibration, browser wiring or cross-platform qualification.


## Oblique propeller loads in coupled aircraft motion

`spine::rotors::propulsion::dynamic::DynamicPropulsion::new_oblique` explicitly
selects native revolution-averaged oblique blade loads over the existing geared
aircraft and uniform axial wake. The normal constructor keeps the original
axial model. Each actual aerodynamic stage supplies incident plus wing-induced
velocity at each disk center and the carrier angular rate to the native blade
owner. Fixed body-frame disk references and physical spin directions are
prepared independently of the signs chosen for mechanical shaft coordinates.

The shared aerodynamic assembly applies the full force at the disk center,
its moment arm about the aircraft COM, and the transverse hub moment. Axial
aerodynamic torque goes only to the original physical shaft port. Consequently
shaft reaction, carrier tilt, transverse loading and aircraft world momentum
use the same real force evaluation without duplicated shaft torque. Both start
and predicted-midpoint stages include oblique loads before advancing the
original aircraft/gear/wake state. No independent mechanics or wing solver is
introduced. The section power account includes hub translation, tilt and
absolute shaft spin; relative shaft power and carrier axial rotation sum to
that same physical absolute-spin work.

The wake remains the named constant-density uniform AXIAL apparent-air-mass
model and its straight axial wash. The scalar momentum/energy receipts are not
a full oblique-fluid impulse or passivity theorem. Incident flow is uniform
across each disk apart from actual carrier tilt kinematics; intra-disk wind
gradients, skew-wake geometry/induction, blade-phase harmonics, flapping and
unsteady airfoil memory remain outside this model. The existing normal wake
and new section-domain refusals both apply; no negative-speed clipping or
reversed-flow fallback occurs.

Three focused native cases in `tests/oblique_propulsion.rs` reconstruct the
actual carrier wrench difference from blade forces/hub moments, advance a
canted two-propeller aircraft with independent physical world momentum
reconstruction, and reverse all mechanical shaft-axis sign conventions while
retaining an identical physical trajectory. The existing dynamic, geared,
spatial and coupled aerodynamic suites also remain applicable. No historical,
full-workspace, browser, real-time or performance qualification is implied.

## First harmonic disk wakes in native aircraft and shaft motion

`spine::rotors::propulsion::dynamic::first_harmonic::HarmonicPropulsion`
explicitly selects a normal-transport mean-plus-two-mode wake for each
physical propeller. `HarmonicState` retains all six induction velocities
alongside the original carrier pose, twist and generalized shaft momentum.
`disk_frames` exposes the prepared body-fixed mode references so callers
can specify the actual physical field. Mechanical shaft-axis sign choices
do not change these mode directions or their physical forcing.

One original wing solve per actual aerodynamic stage now consumes the
local modal induction at each real wing-strip location inside the existing
straight axial wash window. The same field reaches every native blade
quadrature point. Blade thrust and actual transverse hub moments drive the
three modal air-inertia ports; all body forces, lever arms, transverse hub
moments and aerodynamic shaft torques pass through the original shared
wrench and geared mechanics owners exactly once. No duplicate section,
wing, mechanical or nonlinear wake solver is introduced.

The initial force stage predicts aircraft, shaft and all wake modes to
half time; the actual predicted-midpoint stage advances the original
mechanics and each implicit transport step. The wake input work and the
section induced work use their actual respective midpoint states.
`coupling_work_residual_j` explicitly exposes their difference rather than
calling a partitioned quadrature defect dissipation or silently claiming
coupled passivity. The modal frame-rate term uses carrier rotation about
the disk normal, never physical shaft spin. Original immutable-state and
late-failure atomicity remain; callback effects are caller-owned.

The selected wake closure is still a normal-flow, zero-skew transport
approximation even when the blade law resolves incident crossflow. The
existing straight wash geometry is a finite-window transfer approximation;
it is not a free wake, skewed induction solution, full world-fluid momentum
account, wake-distortion model or total fluid/body energy theorem. Incident
air/wing-induced flow is still sampled at each disk center apart from actual
carrier tilt kinematics. Finite blade phase, higher modes, flapping, cyclic
pitch, swirl, direct inter-disk wake coupling and intradisk ambient wind
gradients remain unimplemented. Apparent air inertia remains excluded from
the aircraft's rigid inertia.

Four native cases in `tests/harmonic_propulsion.rs` establish zero-mode
parity with original oblique aircraft loads; nonuniform induction reaching
real blades and wing strips; independent physical aircraft/shaft world
momentum reconstruction over 24 canted-rotor ticks; full-state and coupling
work refinement; mechanical-axis sign parity; and late-density refusal with
deterministic retry. Measured 24-tick residuals are `8.73145e-13 N s` linear,
`1.58448e-11 N m s` angular, `2.94253e-15 N s` modal and `2.10943e-15 J`
within-wake storage/work. At 24/48/96 steps, the complete-state error against
384 steps is `3.56045e-4`, `8.81038e-5`, `2.09938e-5`; accumulated absolute
induced coupling-work defects are `1.45322e-3`, `3.65155e-4`, `9.15029e-5 J`.
The state error is a numerical fixture norm, not a common physical unit or
trajectory certificate. All ten new and 55 original propulsion/airscrew
cases pass in the focused actual-source native build. This is not a full
Cargo/DSR, cross-platform, browser, real-time or physical validation gate.

## Skew-flow harmonic wakes in actual turning aircraft motion

`HarmonicPropulsion::new_skewed` selects the full Pitt--Peters first harmonic
transport from the `fs-airscrew` source owner. `new` retains its original
normal-transport choice. The dimensional source gain, moment signs and
normal-throughflow limits are specified in the airscrew contract and Chen,
NASA TM-102219, equations (13b), (14), (17), (18):
<https://ntrs.nasa.gov/api/citations/19900006622/downloads/19900006622.pdf>.
This extends the earlier explicitly zero-skew transport assumption; it does
not replace the shared section, wing or geared mechanical owners.

Each actual aerodynamic stage projects the same wing/ambient incident
velocity into its fixed physical disk frame. The current transverse
direction changes only the source operator; stored wake modes retain their
body-fixed meaning. `disk_incident_mps` and optional `skew_transport` expose
these actual stage inputs and coefficients. Mean and both harmonics reach
the original radial/azimuth blade points and existing wing-strip wash
window. Their actual thrust and transverse hub moments drive the modal
ports. Rotor body/shaft forces, torque, carrier reaction, lever arms and
motor work still use the original mechanical update exactly once.

The optional `skew_inflow` entries contain the same accepted modal steps
also exposed in `inflow`; they are additional accounting detail, not extra
updates. Source cross work is signed and does not become invented heat.
The original `coupling_work_residual_j` still measures implicit wake input
work minus actual predicted-midpoint blade induced work. A late wake,
blade, wing or input refusal returns no changed borrowed state; callbacks
remain caller-owned. Source iteration caps and tolerances are explicit.

Three additional actual-source aircraft cases exercise a time-varying
crosswind and turning disk velocity; exact shaft-axis/ratio sign parity;
unchanged aerodynamic loads at identical induction but changed derivatives;
independent physical world momentum reconstruction; signed source work;
full-state and partitioned-work refinement; exhausted iteration budgets;
and late refusal/retry. A 24-tick, 96 ms canted-rotor trajectory measures
`2.67193e-13 N s` linear and `1.84520e-11 N m s` angular impulse defects,
with largest modal defect `4.76659e-11 N s`, within-wake energy defect
`1.89884e-11 J`, and at most four coefficient iterations. Its accumulated
signed cross transport work is `-2.21195e-3 J`.

The complete coupled refinement check uses the initial smooth 6 ms interval
of that fixture: 8/16/32 steps versus a 256-step reference give state-norm
errors `1.79455e-4`, `4.49075e-5`, `1.10028e-5` and accumulated absolute
coupling-work defects `5.12082e-5`, `1.22387e-5`, `2.97900e-6 J`. This norm
mixes native coordinates and is a numerical convergence measure, not a
common physical unit. Longer refined versions encounter nonuniform
convergence or the existing quasi-steady wing solver's 600-iteration
refusal; no wing tolerance, cap or closure was relaxed. The 96 ms momentum
ledger is therefore not a convergence certificate for that entire flight.
Independent source-level nonlinear refinement remains separately tested.

All nine new and 65 prior airscrew/propulsion cases pass in the focused
pinned native build (74 cases across twelve suites). This does not claim
full Cargo/DSR, cross-platform, real-time or physical validation. Existing
positive section throughflow/angle domains, disk-center incident sampling,
straight finite-window wash geometry and constant air density still apply.
The extension supplies the skew transport matrix; it does not supply a
transported/skewed far-wake geometry, direct inter-disk wake coupling,
intradisk ambient wind gradients, reverse-flow regimes, higher modes,
swirl, flapping/cyclic pitch, or a complete fluid/body energy/momentum model.

## Actual stage-by-stage wing continuation in skewed flight

`first_harmonic::ContinuedHarmonicState` pairs the existing complete physical
`HarmonicState` with an explicit `fs-wing::nonlinear::continuation::WingContinuation`.
`initialize_continuation`, `evaluate_continued` and `step_continued` select the new
wing continuation lane. The original entry points retain their original cold
wing solve. The actual initial force stage continues from the caller's anchor;
its accepted circulation seeds the actual midpoint force stage. The returned
physical endpoint carries the accepted midpoint circulation as the next initial
guess. Both stages resolve current section equations and apply current loads;
neither force is lagged. The original aircraft/shaft mechanics, blade section
quadrature, skew-inflow owner, force/moment arms and impulse/work accounts remain
the only physical owners. The anchor adds no aerodynamic momentum state.

The exposed long-flight failure had two distinct causes. On the original
`vel_z=0.8 m/s` positive-gust trajectory, the old Picard iteration refuses at
8 ms after its 600-attempt cap even while the section solution is attached and
the row split is well conditioned. Caller-owned warm starts alone did not cure
it. The new native Newton--Krylov continuation crosses that failure at the same
`1e-10` actual componentwise section gate. Continuing the same positive-gust
trajectory later approaches cancellation of a signed linear chord-row total;
the declared row-condition bound correctly refuses. The original
`vel_z=1.5 m/s` fixture also has near-filament and very large effective-incidence
roots; the admitted continuation API explicitly refuses it. A low algebraic
residual on that legacy fixture is not evidence of realistic long-trajectory
accuracy.

The full coupled refinement regression uses the original attached
`vel_z=0.8 m/s` initial state and the opposite smooth crosswind turn,
`wind(t,p)=[0.02*p_y, -0.5*sin(20*t), 0.1*cos(20*t)] m/s`. It completes the actual
96 ms interval at 24, 48, 96 and 384 physical steps with canted counterrotating
disks, moving attitude, local wing velocities, first-harmonic skew transport and
state-dependent aircraft/shaft loads. The section equations at every real
initial and midpoint stage are independently reconstructed from the public full
horseshoe field minus an independent finite-bound field. Their maximum measured
normalized residual is `3.069e-11`; all effective incidences remain attached.
The run does not weaken the force residual or shorten the physical interval.

Against the 384-step reference, the existing full-state distance errors at
24/48/96 steps are `3.4904523143e-4`, `8.6389460900e-5`, `2.0587886173e-5`.
The accumulated absolute induced-work coupling defects are `1.4235880765e-3`,
`3.5729731447e-4`, `8.9494495170e-5 J`. Both refine at the expected second order.
Across those runs and the reference, independent solid-body/rotor world linear
momentum defects stay below `7.48e-12 N s`, angular momentum defects below
`9.03e-11 N m s`, reduced-wake impulse defects below `6.92e-11 N s`, and wake
energy-balance defects below `3.59e-11 J`. No stage uses more than six Newton
attempts or 44 Jacobian products. Reversing equivalent shaft-axis/ratio/blade
handedness conventions gives exactly the same complete continued run and work
account. Late source and exhausted-Newton refusals leave both physical state and
wing anchor unchanged; exact replay is tested.

These measured bounds concern this admitted 96 ms convergence fixture. They do
not certify arbitrary maneuver duration, separated-wing roots, a moving far
wake, unsteady section hysteresis or full fluid/body momentum closure. The
previous explicit induced-work partition defect remains a quadrature defect,
not dissipation or a passivity guarantee. Row-cancellation, incidence and
filament domains remain explicit physical/numerical admission limits rather
than hidden regularization.

The focused final sweep passes 107 cases: 30 wing/continuation cases and 77
airscrew/propulsion cases, including the ten-case harmonic flight battery.
The native harness compiles current wing and propulsion sources and the actual
Newton/FGMRES/Armijo numerics; it omits the solver's two `fs_exec::Cx` convenience
wrappers and unrelated aircraft surfaces. This does not substitute for a full
workspace Cargo/clippy/CI gate, which was not run in this environment.

## Gravity in longitudinal perturbations

`longitudinal::linearize` holds pitch fixed while perturbing body air
velocities. The original `aircraft::force_buildup` retains its level-flight
trim convention, where pitch equals incidence. Its inserted weight is removed
before taking velocity derivatives; actual attitude-dependent gravity enters
only the pitch column. Angle-of-attack changes therefore cannot create
spurious gravity terms in the speed or incidence derivatives. Aircraft force
assembly, linearization, airborne motion and rail force projection share the
same terrestrial gravity constant and weight-removal operation.

The focused G1/G3 checks distinguish incidence from attitude using an exact
gravity-only manufactured field and compare full native aerodynamic force
derivatives after changing mass at a fixed operating point. The latter must
remain mass-independent before division by inertial mass. This corrects the
linearized poles and all augmented/pilot models consuming that matrix; it does
not change the underlying aerodynamic solve or the original level-flight trim.
The historical model still uses uniform standard Earth gravity. Spatial
planetary gravity is not inferred from this correction.

The focused x86 native sweep passes 27 cases: 11 motion tests, seven V-02a
tests, seven augmented-spectrum tests, and the actual assistance/pilot flight
behavior regressions. The original longitudinal source reproduces the previous
V-02a digest and fails the new mass-scaling regression. Correcting its gravity
derivatives changes nominal doubling time from 0.6204376902 to 0.6670863177 s;
all physical acceptance thresholds remain unchanged. Only measured x86 golden
pins were updated; other ISA pins still require native remeasurement.

This source-selected build uses the original aircraft/wing/propeller and
Newton/FGMRES numerics, excluding only unrelated module surfaces and the two
unused solver `Cx` convenience wrappers. It does not establish a full Cargo or
DSR gate; full-runtime builds were interrupted while compiling asupersync.


## Prescribed material atmosphere in the evolving flight lifecycle

`SimLoop::init_with_profile_spatial_wind` and its distributed counterpart bind
an immutable `MaterialAirProfile` to the original native trajectory. Initial
trim and prelaunch state use the source gas at the 2.4 m rail height. Each
beginning/midpoint load stage samples its own CG height, and the same density
drives native wing/canard forces, propeller thrust and torque, aeroelastic
dynamic pressure, and reduced/distributed lateral response. The value in
`scenario().rho_kg_m3` remains the initial rail density. `sample_material_air`
returns the complete original gas at the requested height without advancing
wind or trajectory. Predicted endpoints are profile-admitted before publishing
motion and its corresponding warm start.

V6 checkpoints retain the profile definition, geometric ground datum, sampling
mode and complete existing wind/flight state. Uniform profiles retain all twelve
gas fields; the standard atmosphere retains its model discriminator and ground
altitude. Profile identity binds into tick-zero and run identities. Restore
admits the initial density and saved height before publishing state. Matching
replay recorders retain the initial V6 checkpoint; existing V2–V5 modes keep
their layouts and homogeneous paths retain their existing arithmetic.

The original material owner supplies hydrostatic pressure and coherent gas
thermodynamics/transport. Density is frozen across the aircraft within each
load stage; distributed mode resolves incident velocities across strips and
disks. Existing reduced flight dynamics, constant flight gravity, engine torque
law and Reynolds-independent polar remain in force. This adds no weather
evolution, compressible/anelastic atmospheric dynamics, buoyancy, pressure work
or high-altitude aircraft validation claim. A prescribed solenoidal velocity
field with stratified density does not imply divergence-free mass flux.

Six focused regressions cover native trajectory/rotor changes with altitude,
stage-density coherence through midpoint/warp/lateral loads, homogeneous
equivalence, profile-domain refusal before motion publication, complete source
encoding, and checkpoint/replay continuation. Source/parser review is distinct
from native execution, which remains unverified for this addition.
