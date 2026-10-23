# fs-mbd-material

L3 composition boundary from immutable material data and volume geometry to
`fs-mbd::articulated::SpatialInertia`. The numerical `fs-mbd` owner remains
independent of material databases; no replacement dynamics solver is introduced.

`resolve_body` requires a unique key, volume moments, material card, complete
state point and source selection for every region. Only density is required.
It retains the original resolved bundles and composes their mass, centre of
mass and full centroidal tensors, including products of inertia. Region keys
fix the reduction order. The parallel-axis merge works about the centroids,
not by subtracting origin tensors. Final spatial-matrix admission is performed
by the same downstream inertia owner used for momentum and kinetic energy.

`VolumeMoments` accepts moments from a geometry integration owner and supplies
analytic cuboid and straight-tetrahedron constructors. SE(3) placement uses
`fs-ga`; tensors remain centroidal after placement. Lengths are metres, volume
is m3, volume inertia is m5, density is kg/m3 and mass inertia is kg m2.

The caller supplies a non-overlapping rigid partition. Numeric admission does
not establish geometry coverage, mesh accuracy, rigidity or measured material
accuracy. Source uncertainty is retained, not propagated to a certified inertia
bound. Density at a changed state must be paired with the intended geometry:
compiling a new specimen is not a mass-conserving thermal-expansion update, a
live state transfer, or an implicit source of energy/momentum.

Compilation has a 4096-region hard cap and a caller-selected lower cap. Missing,
ambiguous, wrong-unit, nonpositive or unsupported density refuses. Cancellation
is checked before work, between fixed-size region operations, and before
publication. Sorting and each underlying material query are bounded cold-path
operations, not a real-time cancellation-latency guarantee.

`articulated::resolve_articulated` compiles complete material-defined link
geometry into the existing `ArticulatedModel` without placeholder inertias.
The region budget is global across all links. The existing owner still validates
root-first topology, unique names and physical joints. Bodies and source bundles
remain aligned with the unchanged link order. The resulting model feeds existing
forward/inverse dynamics, forward kinematics and unconstrained free-flight APIs.

`articulated::rebind_articulated` requires exactly one named material binding
for every link of an existing model. It preserves link order, parents, zero
transforms, joint axes and limits. Missing, duplicate or foreign bindings refuse
instead of retaining stale mass values. Publication is atomic and never mutates
the template or a live generalized state. This is model substitution, not an
energy-conserving material-exchange event in a moving system.

The `material_lever` executable evaluates the same geometry, joint and loads
using two explicit synthetic density cards through actual articulated solvers:
`cargo run -p fs-mbd-material --example material_lever -- 1000 3000`.
It reports mass, pivot inertia, torque-driven acceleration and gravitational
holding torque. The densities are demonstration inputs, not material measurements.
`cargo test -p fs-mbd-material --all-targets` includes the source/geometry,
full-mechanism and executable regressions. No independent time integrator,
contact law, stiffness law or empirical rigid-body qualification is claimed.

`ResolvedBody::principal_body` uses the existing L1 Jacobi eigensolver on a
scaled 3x3 inertia tensor and checks reconstruction through the actual returned
quaternion frame before supplying the existing free-body `MassProperties`.
Principal moments are not silently clipped. The centre-relative principal
frame is retained together with the original source-bound body. Its measured
relative tensor residual is distinct from source and geometric uncertainty;
repeated moments have non-unique eigenvectors, not a smooth frame guarantee.

`state_from_kinematics` maps the authored pose and origin twist to principal
centre-of-mass pose and momenta, including omega cross centroid. The reverse
`world_from_body` recovers the original geometry pose after free-body stepping.
`arm_for_body_point` maps original geometric/contact points to the existing
point-kinematics, effective-mass and impulse APIs. This retains one mechanics
owner while material-induced centroid and principal-axis changes propagate
without moving the geometry or misplacing contact forces. The source model
must match the state; the adapter does not infer a live material-exchange law.


## Heat transport to mechanical configuration

`thermal::heating::HeatedBody` binds each original region to its source-resolved
solid heat-capacity chart and an explicit thermal control volume. The compiler
copies the original mechanical mass into the thermal cell, checks the same
material card, density claim, state context and reference volume, and creates
the existing `fs-thermal` enthalpy timeline. Contact/boundary events are retained;
new mass is refused. `advance_to` returns the thermal history, configured body
and original energy audit only after all actual thermal substeps and material
checks succeed. A hot unsupported intermediate state cannot hide behind a
subsequent cool endpoint.

`advance_with` accepts the existing `DepositionObserver` interface for generated
heat and step-local transport. It uses `CoupledObservers`, not a second thermal
solve, so source heat is deposited and accounted once. Source observer history,
enthalpy history and mechanical configuration fail or succeed together. Original
caloric data and thermal-expansion receipts remain accessible.

This is one-way uniform-region solid free expansion, not constrained
thermoelasticity. Thermal areas, reference volumes and contact conductances are
not deformed automatically; this needs a justified small-geometry-change/Biot
approximation. The output inertia feeds existing dynamics, but updating live
momentum, mechanical work and stress is outside this adapter. Source claims,
uniform-region assumptions and refinement still bound accuracy.

Focused native regressions: `cargo test -p fs-mbd-material --test heating`.

## Loaded material heat and shape (MR21/MR70)

`TetrahedralMass::solve_thermal_with_nodal_forces` composes the existing source
resolution and compatible thermal solve with the original solid owner's dead
nodal forces. Current E(T), integrated expansion, reference density, constraints
and exact mechanical loads remain attached to the accepted shape and its
`initial_stress_pencil`. `support_reactions_n` subtracts applied support loads
from the original internal forces; free components remain zero.

`gravity_nodal_forces` integrates original per-region density through
`TetLoadGeometry3` on the same reference mesh and region/element map. Heating
changes actual shape and density without changing reference matter or weight.
`ConstrainedHeatedBody::new_with_nodal_forces` keeps the immutable force vector
through every actual thermal substep, generated-heat observer composition,
cooling and restart. Empty forces retain the thermal-only operator identity.

This is quasistatic, fixed-topology small-strain mechanics with reference-frame
dead loads, not follower-load or inertial dynamics. The heat audit stays thermal;
mechanical work and elastic storage are not silently added to or claimed by it.
Focused synthetic cases exercise density/weight, loaded softening and modes,
correct support reactions, heterogeneous element order, heat/cool/restart and
refusal. Native Rust execution remains a separate verification requirement.

`ThermalCell` and `ConstrainedRegion` select the same explicit
`ThermalStiffnessModel`: isotropic E/nu or nine oriented orthotropic engineering
constants. The existing material-law owner admits the orthotropic compliance;
the existing solid owner rotates it once into the mesh's Mandel coordinates.
Both consume the original reference density, preserving closed specimen mass
without demanding a new density claim at the heated temperature. Source data,
rotation and orientation identity remain attached to the accepted state and
its prestressed vibration pencil. Missing orthotropic properties or invalid
axes refuse; there is no isotropic fallback. Every accepted heat substep uses
that same model, including after restart and cooling.

Both stiffness rungs currently require a separately sourced **isotropic**
thermal-expansion path. Orthotropic elasticity alone does not establish that
this expansion approximation applies to a specimen. This adapter does not
claim anisotropic expansion, full tensor material support, large-strain
mechanics or experimental qualification. Focused G1/G3 cases in
`tests/thermoelastic_body/orthotropic.rs` compare the isotropic limit, analytic
directional compliance, heat/shape/inertia/vibration propagation, restart and
atomic refusal of an unsupported intermediate temperature.
