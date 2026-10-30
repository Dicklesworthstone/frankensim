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

`with_expanding_journal_bearings` additionally binds each cylinder to its exact
child/parent material region. Its reference cylinders are coaxial with, and
centered on, the revolute joint; their free-expansion fixed points must match
that joint in their own frames. Each endpoint takes the original regions'
integrated expansion receipts, scales the original radii and takes the shorter
expanded axial overlap. The same geometry drives resistance, torque, Reynolds
admission and the next interval's integrated heat. A closed gap or excessive
geometric resistance change at an endpoint refuses the entire coupled step.
Geometry never compounds across updates, and each cylinder retains its own
temperature. These sampled checks do not localize intermediate contact events;
the caller's timestep and model regime must exclude unresolved gap closure.

The film is massless, isothermal and quasisteady. Geometry is either explicitly
prescribed or follows this declared free-expansion model; it is not inferred
from inertia moments. No pressure-supported load, eccentricity, cavitation,
turbulence, end effects, film thermal storage or spin-up is claimed. Radial
motion/squeeze-film work and constrained bearing stress are outside this
circumferential-shear model; callers must justify neglecting them. The
thermal/mechanical split remains first order. G1/G3/G4 tests in
`tests/heated_articulated_viscous_motion.rs` cover analytic spin-down/local heat,
two-link angular momentum, source-only viscosity feedback and refinement,
source expansion through clearance/motion/heat, independent cylinder heating,
continuation, gap closure, coefficient budgets and late rollback. Link inertia
parameter work remains separate from viscous heat; expanding a film does not
create a second heat or mechanical-work account.

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

## Material interface damping and finite-body thermal feedback

With `dissipative-coupling`, `thermal::heating::motion::damping::DampedRigidPair`
retains two original `HeatedRigidBody` trajectories and their finite enthalpy
networks. An ordered `ResolvedLinearInterface` supplies a pure axial viscous
connector: explicitly sourced zero stiffness, positive damping over its fixed
declared area, and no constant-loss anchor. Attachments name original material
regions and reference body points. Existing source expansion moves those points;
their actual world separation determines the central impulse direction.
The caller establishes that points belong to the named regions and that the
lumped interface law represents the physical connection.

`MaterialDamper::temperature_side` explicitly selects the receiving cell whose
temperature supplies the interface law. Other original coordinates, surface
roles, property roles and selected claims are retained. `first_heat_fraction`
routes physical heat to side 0; side 1 receives the exact remainder. Neither an
average temperature nor an effusivity partition is invented. The connector has
no separate heat store. Original links and reservoirs remain with each body's
thermal network. `new` defaults to insulation between the bodies;
`with_thermal_contact` adds the ordered source-bound exchange described below.

Each step solves both public `fs-thermal::BackwardEulerStep` residuals using the
existing Newton-Krylov and fallible finite-difference owners. Unknowns are cell
energy increments in each body's original sorted order. Every trial starts from
the same accepted bodies and refreshes the source at the proposed temperatures.
A fresh final evaluation must pass local, global, and mechanical/thermal energy
gates before either history is returned. Unsupported sources, pending events,
nonlinear budgets, physical step limits and cancellation publish no endpoint.
The returned pair retains the accepted source for subsequent continuation.

For starting axial relative velocity v and summed native inverse effective mass g,
the midpoint impulse on the second body is J=-c*dt*v/(1+c*dt*g/2). The native
equal-and-opposite point-impulse owner updates both momenta. Independently
calculated heat D=c*dt*(v/(1+c*dt*g/2))^2 must agree with both its kinetic-energy
change and midpoint impulse work before it enters thermal storage once.
Central impulses conserve linear and total angular momentum in exact arithmetic.
The admitted c*dt*g ceiling is at most 2, refusing scalar slip reversal from an
excessively coarse damping step. The original frequency domain must include
zero and the resolved relaxation rate; this does not establish broadband
physical validity or remove time-refinement obligations.

After the impulse, the original integrator advances actual positions and
orientations. Existing expansion/full-inertia transfer then updates configuration
while preserving world momenta. This is a first-order split. Its audit separates
physical dissipation, independently evaluated reservoir/authored heat, expansion
parameter work and native numerical drift. Mechanical heat generation is exactly
the physical damping loss. Expansion parameter work remains the existing explicit
unresolved thermoelastic-work boundary: this closes dissipative feedback, not a
general thermoelastic free-energy, entropy or coupled-passivity model.
Spring storage, dry friction, normal contact, connector mass and deformation
kinetics require their own owners.

The focused regression target is `heating_damping` with
`--features dissipative-coupling`. Manufactured sources exercise midpoint
impulse/heat, material-temperature feedback, continuing trajectories and joint
refusal. These tests specify implementation checks, not measured validation;
native execution must be reported separately from source review.

### Ordered thermal exchange through the same damping interface

`DampedRigidPair::with_thermal_contact` binds the shared `ThermalInterfacePoint`
to the existing complete interface bundle, which must already contain the
canonical area-specific thermal-contact resistance alongside stiffness and
damping. Admission retains the same interface card, ordered actual receiving
bulk cards, original claim pins, coordinate schema and fixed external
conditions. The declared side selected by `temperature_side` supplies the
mechanical layer-temperature coordinate `T`; both side temperatures, optional
phase fractions and any explicitly declared mean follow their actual caloric
states. No mean temperature is inferred.

The original mechanical interface area A and sourced resistance R'' [m² K/W]
supply G=A/R'' [W/K]. At every joint backward-Euler trial the complete original
bundle refreshes together, and Q=dt G (T_A-T_B) enters the first network as -Q
and the second as +Q. Each network retains its original cells, conduction,
reservoirs and authored sources. Generated dashpot heat keeps its explicit
partition. `contact_heat_to_second_j` and `contact_conductance_w_k` expose
internal transfer separately from dissipation, external heat, expansion
parameter work and numerical drift.

This is a fixed-area, negligible-storage connection between lumped thermal
nodes. It infers no bulk half-cell resistance, evolving gap, contact active set
or pressure law. Source refusal, cancellation or failed energy closure returns
neither proposed body. Cloning and the third `into_parts` component preserve
the contact binding. The `heating_damping` regressions cover unequal starting
temperatures in both side orders, data-only nonlinear resistance swaps that
alter both thermal histories and motion, lossless restart and atomic refusal.

## External impulses on continuing heated material bodies

`thermal::heating::motion::loads::MaterialAttachment` identifies an original
material region and a point in reference body coordinates. The current source
expansion maps that point about the region's original fixed point; the existing
principal-frame owner supplies its actual centre-of-mass-relative arm.
`HeatedRigidBody::material_point_kinematics` exposes that same world point,
velocity and angular velocity for an external load or contact owner.

`apply_material_impulses` applies an explicitly bounded, ordered batch of native
point-impulse events to a private candidate. Each impulse uses the existing
`fs-mbd::RigidBodyState::apply_impulse_at_body_point` owner and the current
source-derived mass and full inertia. Signed external work is independently
evaluated as impulse dot midpoint point velocity and compared with native work
and the actual kinetic-storage change, at both event and batch level. Extraction
of mechanical work remains signed; it is not relabeled as generated heat.

Accepted events retain the original thermal cells, sources, event cursor,
configuration, pose and physical clock. The returned body continues through the
same thermal/mechanical advance, so subsequent heating updates its inertia
without discarding the applied momentum. Failed attachment admission, native
dynamics, work bounds or cancellation publish no partial batch. The audit keeps
the ordered native receipts and separate external-work/storage accounts.

This is an instantaneous mechanical-input boundary. The caller supplies the
physical impulse law, the other body's reaction and any separately resolved
impact/friction heat. It does not infer a contact active set, a force-holding
duration or a continuous-force time integrator. Material coordinates do not prove
geometric contact or membership inside a region; the geometry owner supplies
that physical correspondence. Focused tests in `heating_loads` exercise analytic
off-centre work/momentum, continued heating, source-only density and expansion
changes, and atomic refusal. Native execution is reported separately.

## Source-bound creeping-flow motion of a sphere (MR57)

`VolumeMoments::sphere` supplies analytic volume and centroidal moments to the
existing material-to-inertia compiler. `creeping_flow::CreepingSphere` retains
the original `FluidTransport` source, resolving its coherent density and
dynamic viscosity once for a frozen isothermal reservoir. Its declared radius
must agree with a single uniform spherical material region's volume and full
inertia. This comparison rejects incompatible bodies but does not certify a
shape from moments alone; the spherical boundary is explicit model data.

The Stokes force is -6*pi*mu*a*(v-U), rotary torque is -8*pi*mu*a^3*omega,
and fully immersed buoyancy is -rho_fluid*V*g. Body density therefore changes
inertia and settling, fluid density changes buoyancy and regime admission,
and viscosity changes both translational and rotational response. No material
name, empirical drag multiplier, reference mass or independent damping setting
selects those quantities.

`step` closes the existing `RigidBodyIntegrator::step` midpoint force relation
algebraically and passes the resulting wrench to that original owner. It never
advances a second momentum, position or attitude state. Endpoint and midpoint
loads are evaluated from the corresponding actual velocities. Viscous heat is
dt*(6*pi*mu*a*|v-U|^2 + 8*pi*mu*a^3*|omega|^2) at the midpoint. Maintaining a
nonzero uniform background velocity supplies dt*F_drag dot U. An independent
balance checks kinetic, gravitational and buoyancy-potential change plus
viscous heat minus background-flow work; numerical residual is neither
physical heating nor source work.
The fluid is an external reservoir, so these reported transfers must not be
deposited a second time by a caller claiming the same accounting boundary.

This bounded reduction admits translational and rotational Reynolds limits
no greater than 0.1. The fluid diffusion time a^2/nu relative to the fastest
body relaxation time must fit a separate caller-declared limit no greater
than 0.1. The physical step is bounded by a declared fraction of that body
relaxation time. These are explicit approximation/admission limits, not
certified continuum errors. Time refinement is required even when energy
closure is small. The conservative timescale gate excludes, in particular,
many density-matched solid/liquid transients that need fluid inertia/history.

The model is fully immersed, unbounded, no-slip, incompressible, Newtonian and
quasi-steady. It supplies no Basset initial-layer/history force, added mass,
walls, particle interactions, free-surface forces, rarefaction, turbulence,
thermal evolution or other exterior/gap resistance. Continuum validity and
actual spherical geometry remain caller obligations. Source/geometry/regime
failure, work-budget failure or cancellation returns no proposed state.

G1/G3/G4 synthetic tests in `tests/creeping_flow.rs` compare actual repeated
native trajectories with independent settling and rotary-decay formulae,
demonstrate source-only body/fluid substitution, check second-order time
refinement, separate maintained-flow work from heat, and exercise geometry,
regime and transactional refusals. This is a bounded MR57 capability; it does
not close the broader exterior/gap-flow bead or claim empirical validation.

## Continuing hydrostatics during material-resolved thermal expansion

`HeatedRigidBody::advance_hydrostatic` composes the original heat scheduler,
clipped-volume hydrostatics, native kick-drift mechanics and momentum-preserving
thermal configuration transfer. Named `ThermalDisplacementRegion` envelopes
follow their own original integrated free strain and fixed point. Each rebuild
starts from reference vertices; closed-body mass remains unchanged. Fluid
density is the exact `FluidTransport` snapshot, and reservoir gravity must
match the original mechanical integrator.

`HeatedHydrostaticAudit` includes pressure-potential storage and separately
reports its change caused by expansion as external configuration work, never
caloric heat. Absolute parameter-work sums count inertia and displacement
contributions before cancellation. The numerical budget bounds the native
kinetic-plus-gravity-plus-pressure-potential defect. Failed bindings, sources,
work limits and late cancellation publish neither proposed trajectory.

The fluid reservoir has frozen homogeneous density and a prescribed surface;
the coupling uses a first-order split. It adds no waves, drag, added mass,
flooding, expansion-induced fluid kinetic energy, constrained deformation or
reciprocal thermoelastic heat. The caller owns envelope coverage, nonoverlap
and compatibility of free regional expansions. Replacing the reservoir or
envelope between calls changes the model; no exchange work is inferred.

G1/G3/G4 definitions in `tests/heated_hydrostatic_motion.rs` check material-alpha
changes reaching flotation, independent regional center of buoyancy and torque,
source-density acceleration, reference rebuilding and restart, and late
cancellation/work-budget rollback. Synthetic coefficients test the implementation;
they do not establish measured material or hydrodynamic validity.

## Continuing creeping flow during material-resolved thermal expansion

`HeatedRigidBody::advance_creeping_sphere` composes the original thermal
timeline, source free expansion, `CreepingSphere` load/midpoint owner and
momentum-preserving configuration transfer. A single uniform spherical region
is required. Its current radius is rebuilt from the declared ORIGINAL sphere
and the same integrated strain that supplies its volume and inertia. Geometry
is never compounded and the original closed-body mass stays bit-identical.
The original fluid snapshot supplies all density and viscosity values without
another source query when radius changes; its gravity must match native gravity.

Consequently Stokes translation, rotary resistance and displaced-fluid force
change as r, r^3 and r^3 while the body's moment of inertia changes as r^2.
The original owner retains its force closure, native position/momentum/attitude
step and all Reynolds, diffusion-time and relaxation-step gates. Admission
also runs at each UPDATED thermal endpoint, even when there is no subsequent
mechanical step in that advance. `advance_creeping_sphere_with` composes an
existing heat/transport observer on these exact substeps; its generated heat,
the finite material history and continuing motion publish atomically.

The audit separates external-fluid viscous heat, maintained-background-flow
work, inertia-remap work and the change of buoyancy potential caused only by
configuration. The latter uses U_b = rho_fluid V g dot x_COM, with the native
world-origin gravitational datum; this unbounded reservoir has no declared
absolute pressure, so that parameter account is NOT a thermodynamic p dV law.
Absolute parameter-work budgets sum inertia and displacement contributions
before cancellation. The native fixed-configuration mechanical/fluid defect
remains numerical error and is never turned into heat. Neither Stokes loss nor
configuration work is silently added to the finite solid's caloric storage.

This is a first-order thermal/mechanical split with the original uniform-cell,
reference-thermal-geometry and slow free-expansion assumptions. It adds no
resolved expanding-surface fluid flow, radial viscous heating, fluid kinetic
energy, added mass/Basset history, reciprocal thermoelastic work, free surface,
phase change or constrained deformation. The original low-Reynolds and fluid
diffusion gates remain model-admission limits, not certified error bounds.
Changing the frozen reservoir between calls changes the model; its exchange
work is not inferred. Source, geometry, flow-regime, work or cancellation
refusal discards the full proposed heat, motion and optional source history.

G1/G3/G4 definitions in `tests/heated_creeping_motion.rs` compare actual
translation, spin and external heat/work against an independent midpoint
recurrence with source-only alpha, viscosity and fluid-density substitutions.
They also cover exact frozen-geometry native-step parity, reference rebuilding,
restart, shared generated-heat substeps, updated-radius regime refusal, late
rollback and retained fluid snapshots. Synthetic coefficients establish only
implementation checks; execution and physical validation are reported separately.

## Material-fluid heat exchange on the continuing sphere trajectory

`advance_creeping_sphere_in_thermal_fluid` and its `_with` composition use the
original enthalpy observer and `CreepingExpansion` to close ambient heat,
temperature, source expansion, radius, inertia and fluid-loaded motion. The
steady exterior spherical-conduction limit supplies G = 4*pi*k_fluid*r and
signed heat into the body G*(T_fluid-T_body). Current expanded radius enters
every exchange; the ambient rho/Cp/k/mu tuple is identical to the motion view.
The original body card and prepared conductivity path supply the Biot check.

Admission requires a wholly exposed spherical region, no overlapping authored
Robin/contact surface, low translational/rotational/expansion Peclet, low
Rayleigh and Biot numbers, slow thermal/radius evolution and small temperature
contrast for frozen fluid properties. Gas buoyancy expansion comes from its
original ideal EOS. Unknown liquid expansion refuses nonzero gravity with
nonzero thermal contrast, including contrast generated inside a candidate.

Signed reservoir heat enters the existing external-heat account once. Fluid
viscous loss and inertia/buoyancy configuration work retain separate accounts.
The heat-exchange layer evaluates prepared solid curves and cached caloric
slope, adding no runtime material queries.
The additional explicit relaxation limit composes with the original timeline
and any generated-heat observer, whose history publishes atomically with heat
and motion. Source/regime/boundary-event/cancellation failure discards all of
those proposed histories.

Seven focused integration definitions cover independent exponential heating
and cooling with refinement, k/Cp/expansion substitutions, shared heater work,
pressure/Rayleigh effects, high-Peclet refusal at low Reynolds, source/boundary
ownership and late rollback. Their native test target passes `cargo check`;
runtime execution remains pending. This is an established exterior conduction
field and an infinite prescribed reservoir. Initial fluid diffusion, finite
fluid storage, advective or nonlinear convection, phase change and reciprocal
thermoelastic heat are not supplied.
