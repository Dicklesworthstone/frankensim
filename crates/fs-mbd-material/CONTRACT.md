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
