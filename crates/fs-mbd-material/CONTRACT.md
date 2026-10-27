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

Thermal symmetry is now selected independently through `ThermalExpansionModel`
on both `ThermalCell` and `ConstrainedRegion`: the original isotropic path or
three source-integrated principal coefficients with explicit fixed axes.
`ThermalExpansionSource` retains the original integrated source and, for the
directional model, the existing solid owner's complete frame-transform receipt.
Its accepted mesh-frame `TetThermalStrainState` is used unchanged by equilibrium,
stress recovery and the retained initial-stress vibration operator. The scalar
route preserves the original arithmetic and identities. The two material
models may have different explicitly authored principal frames.

`ThermalCellSource::expansion()` exposes the selected source variant;
`thermal_strain()` exposes the exact accepted state. The latter avoids callers
reconstructing a scalar approximation. Geometry and full inertia follow the
actual compatible displacement while reference mass remains conserved. Every
accepted constrained heat substep, restart and cooling endpoint retains the
same expansion model, axes and source selection. Unsupported expansion paths
reject the entire proposed advance even when the caloric source remains valid.

Focused cases in `tests/thermoelastic_body/orthotropic.rs` cover the isotropic
limit, analytic rotated strain and fully restrained stress, compatible free
deformation and inertia, changed prestressed operators, heat/cool/restart and
atomic refusal. The constitutive fixtures are synthetic. This is small-strain,
fixed-principal-frame expansion and one-way quasistatic heat-to-mechanics
coupling; it supplies no general off-diagonal coefficient law, finite-strain
evolution, deformation-dependent heat geometry or experimental qualification.
Native Rust execution remains a separate verification requirement.

The public stiffness/expansion model names remain available through the same
`thermoelastic` API, but now re-export the common
`fs-solid::linear3::source` implementations used by spatial enthalpy coupling.
There is no second source resolver, rotation formula or scalar expansion path.
The reference density passed to that shared lowering is the same retained
density as before; existing prescribed-temperature arithmetic and identities
are unchanged.

## Continuing free-body motion under finite material heating

`PrincipalFrame` is the shared full-inertia lowering used by the original
`PrincipalBody` and by changed configurations. `transfer_state_from` retains
world centre-of-mass position, linear and spin angular momentum, and authored
body orientation when the centroid and principal axes change. The geometry
origin translates to retain the physical mass centre. Equal, bit-identical
mass is required; the composing owner must additionally establish unchanged
matter and a common body frame. The returned signed kinetic parameter work is
explicit and bounded. It is neither heat nor damping, and principal-axis
relabeling is not a new physical orientation.

`thermal::heating::motion::HeatedRigidBody` establishes that lineage from one
`HeatedBody`, retaining its original source/material inventory and thermal
timeline. Every actual thermal substep advances the existing free-body
integrator using the current inertia, resolves the next source expansion and
transfers the continuing momentum to that configuration. A mechanical step
limit tightens the same thermal scheduler; events and cooling retain one clock.
`advance_with` composes the existing generated-heat/transport observers without
another heat solve. Source failures, cancellations and work-budget failures
discard the entire proposed thermal and mechanical history.

The audit separates caloric heat, signed configuration work and fixed-inertia
numerical energy change. Limits apply to sums of absolute work/error during an
advance, so heating followed by cooling cannot hide an exceeded allowance.
This is one-way thermal forcing of free flight under uniform gravity. Internal
deformation velocities, thermal stresses, support/contact work and reciprocal
thermoelastic heat are omitted; configuration work is not silently debited from
the thermal cells. The retained thermal fixed points are in the moving body
frame, not world-space clamps. Uniform-region and small-geometry-change
assumptions remain. The split is first order in changing inertia and requires
time refinement independently of energy closure.

G1/G3 tests in `principal` and `heating_motion` cover full-tensor/centroid changes,
world angular momentum including its orbital contribution, source heat-capacity
changes to actual spin, nonuniform expansion, heating/cooling, restart equality,
and atomic source/work/cancellation refusal. Synthetic material fixtures test
the implemented mechanics; they are not empirical material qualification.

## Continuing articulated motion under material heating

`thermal::heating::articulated_motion::HeatedArticulatedMotion` binds the
existing `HeatedArticulated` inventory, source heat timeline and accepted joint
positions/velocities. Each actual heat substep advances `fs-mbd` ABA dynamics
through the existing checked `fs-time` RK45 owner, with stationary base,
constant declared joint efforts and uniform gravity. Position and rate error
scales are independent per coordinate, including mixed prismatic/revolute trees.

At a thermal configuration change, the native RNEA/ABA mass actions preserve
`p = M(q) qdot` at unchanged joint positions. This transfers the coupled
generalized momentum, not each link's angular velocity in isolation. The native
joint-position and velocity limits apply to the resulting state. Momentum is
not checked as actuator force; actual applied efforts retain their own limits.
No dense mass matrix or independent dynamics/time solver is introduced.

Thermal mass and original source ownership remain unchanged. Mechanical kinetic
and gravitational storage, exact held-effort work, signed thermal-configuration
work and fixed-configuration integration error are separate accounts. Total
absolute parameter work and numerical defects have caller-specified bounds;
opposite heating/cooling signs cannot conceal excessive unresolved work.
Parameter work is external to this one-way caloric model, never a second heat
credit. Joint momentum residuals use componentwise dimensional tolerances.

The original heat scheduler owns both clocks, generated sources and contact
events. `advance_with` composes ordinary thermal observers; source, material,
joint-limit, work-budget and late cancellation refusals discard the complete
proposed trajectory. Native recursions are bounded by the source mechanism's
global region/link cap; cancellation is checked at every RK stage and thermal
substep. Clone retains the full restart state.

Scope is the existing free-expansion rigid-link approximation with prescribed
joint placements and transport geometry. The thermal/mechanical split is first
order and requires step refinement. No moving-base or closed-loop mechanism,
contact/impact, shape-change kinetic energy, reciprocal thermoelastic heat or
experimental qualification is claimed.

Focused G1/G3 tests: `cargo test -p fs-mbd --test articulated_mass` and
`cargo test -p fs-mbd-material --test heated_articulated_motion`. Synthetic
coaxial two-joint cases independently specify the full coupled inertia law,
heat-capacity-only swaps, applied work, cooling/restart and transactional
source/cancellation/work-budget refusal. Full workspace qualification remains
a separate gate.

## Material-driven journal loss in continuing articulated motion (MR20)

`HeatedArticulatedMotion::with_journal_bearings` binds explicit revolute links,
receiving thermal cells, concentric geometry and source-pinned
`LiquidTemperatureLaw` instances. Native ABA/RK stages integrate the relative
resisting torque and its conjugate nonnegative power using the same coefficient.
Each accepted interval's integrated loss enters its receiving enthalpy cell
once. Endpoint viscosity/density refresh from the actual wall temperature and
original claims; the next interval uses that material state. All retained
liquid properties require continuous temperature-path source support.

`Audit.viscous_heat_j` is internal mechanical-to-thermal transfer already present
in the thermal owner's supplied heat. The combined audit cancels it, retaining
applied joint work, inertia-remap parameter work and numerical energy defect
as separate accounts. Native stages enforce a creeping-flow bound; endpoints
enforce the declared viscosity-change limit. Every source, budget, mechanical,
thermal or cancellation refusal leaves the accepted coupled history intact.

The film is massless, isothermal and quasisteady with prescribed concentric
geometry. No pressure-supported load, eccentricity, cavitation, turbulence,
end effects, film thermal storage or spin-up is claimed. The thermal/mechanical
split remains first order. Four G1/G3/G4 tests in
`tests/heated_articulated_viscous_motion.rs` cover analytic spin-down/local heat,
two-link angular momentum, source-only viscosity feedback and refinement,
continuation, and late rollback. Native execution remains required for
numerical qualification.

### Material dry joint brakes and two-body thermal feedback

`articulated_motion::dry::DryJointBrake` joins the existing ordered dry-interface,
tribology, articulated dynamics and finite-enthalpy owners. `with_dry_joint_brakes`
attaches prescribed parent-child sliding contacts to `HeatedArticulatedMotion`.
The existing ABA/RK45 stages apply native `FrictionLaw::regularized_traction_1d`
and integrate its mechanical power. Those same accepted loss increments pass
through native `DissipationStep`/`HeatPartition` into the two actual finite thermal
regions, using the existing shared substeps. No separate motion or heat solver is
introduced, and journal films and dry brakes can act together.

For a revolute joint, slip is `v = r qdot` and resisting generalized effort is
`r F`; for a prismatic joint they are `qdot` and `F`. This conjugate mapping
preserves the full articulated mass coupling and the opposite parent reaction.
Normal load, nominal contact area, revolute effective radius and heat fractions
are explicitly prescribed. The preload has no normal displacement/work in the
admitted coordinate. The source pressure is exactly load divided by nominal area.
The ordered interface is bound by `fs-contact`; surface A must be a thermal region
of the named child, B a region of its direct parent, retaining their original
bulk cards and caloric objects. There is no unnamed third heat destination.

Static/kinetic coefficients come from the existing canonical interface property
resolver. A scalar must explicitly declare `ConstantWithinValidity`; an exact
point sample is not promoted to a law over evolving temperature.
`ThermalInterfacePoint` keeps the two body temperatures distinct and
only supplies a mean temperature when explicitly declared. Each accepted endpoint
refreshes the original selected claims. The next mechanical interval therefore
uses the actual heated interface response, including thermally weakened friction.
Rate, frequency, field, strain and history-dependent responses/applicability
refuse because this slice has no owner for those evolving coordinates. Coefficient
curves and sampled surfaces may use one shared evolving temperature coordinate
(A, B or declared mean), with additional held pressure/composition coordinates.
Two evolving source-response temperatures refuse instead of treating point samples
as proof of support on a diagonal. Native source-line queries check continuous
support across every crossed knot/cell and the union of both coefficient breakpoints
checks `mu_k <= mu_s`. Per-step coefficient-change limits require split refinement.

Native regularization is an explicit positive speed and gives zero force at rest;
this does not claim static sticking, a contact onset/normal-compliance solve, wear,
texture-direction dependence, flash temperature, or a derived heat partition.
The two finite regions remain spatially lumped and solid. Thermal contact exchange
between them can be supplied separately by the existing thermal contact owner.
Normal force, area and radius do not change with expansion. Existing fixed-topology,
stationary-base and thermal-parameter-work limits still apply; reciprocal
thermoelastic work is not invented. The split remains first order and requires
time refinement independently of RK tolerance.

`Audit::friction_heat_j` is the mechanical loss already included in the thermal
owner's external heat counter; it is an internal transfer in the combined balance.
Numerical integration defects and inertia-parameter work remain separate accounts.
Each brake also retains the native cumulative two-destination `WorkLedger`.
There are at most 64 journal/brake loss channels in total. Every attachment and
advance is a candidate; source support, phase, mapping, coefficient, energy,
integration-budget or cancellation failure preserves motion, both heat stores,
source receipts and work history together. Source-line work polls the existing
checkpoint callback at bounded property/knot boundaries.

`tests/heated_articulated_dry_motion.rs` contains independent cube-inertia constant
deceleration with interface-data swaps and reversed motion, a coupled two-rotor
momentum/heat balance, prismatic translation, low-speed exponential decay including
rest, thermal-fade convergence against the analytically reduced coupled ODE,
aligned restart, wrong-body refusal, source-domain refusal and late cancellation.
These are synthetic constitutive fixtures, not measured brake validation. This
delivery has parser/source-review checks only; execution of these new regressions
is not claimed without a current-source compile/runtime result.
