# fs-solid CONTRACT

## Current-area pressure on finite-deformation membranes

`MembranePatch::solve_with_pressure` uses the original membrane material,
kinematics and Newton--Krylov owner with a uniform held pressure jump. Pressure
acts on current oriented triangle areas, including the exact matrix-free
follower-load tangent and loads at prescribed nodes. Positive pressure acts
along the reference-oriented normal (initially +z); it is a signed pressure
difference, not absolute thermodynamic pressure. Nonzero-pressure equilibrium
requires a consistently oriented surface and a fully prescribed perimeter.

`pressure_work_j` integrates current-area virtual work exactly along a straight
nodal-displacement segment. At fixed perimeter this is pressure times swept
volume. The equilibrium potential uses zero free displacement at the same
prescribed perimeter as its pressure-work datum. Moving-edge segment work is
path-dependent and is not advertised as a global potential.

`solve_vibration_with_pressure` solves that same equilibrium and uses
`K_material + K_geometric - dF_pressure/du` with the conserved reference mass.
Its additional pressure/volume port reuses the exact current-area load owner.
This holds the background pressure fixed; it does not infer gas compressibility,
reservoir volume, flow, radiation or a finite-amplitude dynamic trajectory.
The existing material-domain, cancellation, bounded solve and full-spectrum
admission rules still apply. G1/G3 tests in `tests/membrane/pressure.rs` check
analytic pyramid deflection/reactions/work, deformed-area frame covariance,
every pressure-force derivative and the loaded modal pencil with two interior
nodes (nonzero follower stiffness). No experimental membrane validation or
global injectivity/continuum stability certificate is claimed.

## Variable-amplitude critical-plane Miner damage on native J2 histories

The quasistatic and backward-Euler J2 fatigue owners now admit
`J2VariableAmplitudeGauge` alongside fixed-plane CMB, representative-cycle SWT,
and Fatemi-Socie gauges. Every accepted native endpoint contributes its actual
stress and total mechanical-strain tensor to the material-level
`CriticalPlaneMiner`; no external tensor history or reconstructed equivalent
stress is accepted.

Each gauge retains its original strain-life calibration/source binding and
declared finite plane set. Uncorrected, Morrow and Smith-Watson-Topper policies
reuse the existing bounded rainflow/Palmgren-Miner owner. Paired policies keep
their existing refusal for phase-lag, relaxation, or nested reductions that
would hide unmatched stress extrema. Sample and samples-times-planes work
budgets remain explicit.

Native mechanical advancement and every fatigue bank are one transaction.
Quasistatic block failure, inertial-step failure, cancellation, source mismatch,
or a late variable-amplitude gauge refusal restores the previously accepted J2
state and every fatigue history. Finalization is likewise atomic across all
banks. The variable-amplitude bank does not require or manufacture representative
cycle closure and is not reset when a separate SWT/Fatemi-Socie cycle window is
restarted.

This is a finite-plane multiaxial variable-amplitude life estimate under linear
Miner summation. It does not establish cyclic stabilization, optimize the
continuum plane sphere, infer unsampled extrema, model load-sequence interaction,
or evolve fatigue damage into stiffness, heat, fracture topology, or crack
growth. The existing fixed-temperature small-strain J2 and source-validity
limits remain authoritative.

## Shear-strain Fatemi–Socie life from native J2 cycles

The shared quasistatic and inertial J2 fatigue observer now admits
`J2FatemiSocieGauge` alongside fixed-plane CMB and critical-plane SWT gauges,
or as the only fatigue bank. Both owners consume the original native stress
and total mechanical strain, including plastic history and unforced inertial
motion. Existing cycle closure, conditioning restart and finish operations
include both critical-plane banks; no alternative constitutive or motion
solver is introduced.

Admission binds the actual immutable material identity, temperature, native
shear modulus and initial yield stress. Source-backed calibrations additionally
match the original card, named state, typed operating coordinates and exact
selected E, nu and initial-yield answers. The evolving hardened yield surface
cannot replace the fitted initial-yield reference. Original material and
calibration evidence remain inspectable.

A native block and every fatigue bank publish together only after requested
closure succeeds. Failure restores mechanics, fixed/SWT snapshots, each
Fatemi–Socie accepted prefix and any provisional result. Ordinary advances do
not copy the retained shear tensor histories. Restart samples the actual
accepted mechanical state and preserves fixed-plane counters and the global
sample budget. A decaying or ratcheting path that fails closure cannot silently
become a stabilized periodic-life result.

This extension predicts bounded empirical shear-controlled life for a declared
resolved representative cycle at fixed source conditions. It does not establish
stabilization or experimental component validity, count multiaxial
variable-amplitude cycles, or evolve fatigue stiffness loss or crack geometry.

## Findley shear-path fatigue from native elastic and Newmark responses

`fatigue::TetFatiguePath::new_with_findley` and
`TetTransientRun3::with_fatigue_and_findley` attach Findley gauges alongside
fixed-plane S-N gauges, or use Findley gauges alone. Both feed the same shared
observer from the original native tensor stress recovery, including actual
inertial stress under zero external load and prescribed-support motion.
Physical tensor shear is retained in the declared mesh frame. The calculation
does not replace the elastic solve, Newmark step or material constitutive law.

A sourced Findley gauge requires the original immutable resolved material
bundle, with identity equal to the actual native elastic material identity.
Its original card, named material state and every calibration query coordinate
and quantity convention must match. The retained thermal state and fatigue
calibration must have the same temperature. A label or caller-supplied hash
alone cannot authorize a source-backed transfer; derived/oriented material
identities need their own explicit original receipt mapping before admission.

Explicit window closure evaluates the full stress path. A combined endpoint
advance and closure accepts native mechanics, S-N history and all Findley
gauges together. Any late failure or cancellation restores every observation
and leaves native motion unchanged. Ordinary steps append bounded tensor
samples without cloning the entire preceding Findley history. A restart begins
at the last actual accepted native stress, preserves S-N and mechanical history,
and does not reset the global sample budget. Completed Findley usage factors
are exposed independently from S-N/Miner summaries.

This is fixed-temperature, small-strain elastic fatigue-limit usage on resolved
closed representative blocks and declared candidate planes. It does not infer
periodicity for decaying vibration, count multiaxial variable-amplitude cycles,
predict finite life or couple fatigue to stiffness loss or crack propagation.
Synthetic native regressions cover rotating shear, unforced vibration,
source binding, restart and all-gauge rollback against the original solvers.

## Critical-plane SWT from accepted native J2 cycles

The quasistatic and inertial J2 fatigue owners can carry bounded critical-plane
cycle gauges alongside fixed-plane rainflow gauges, or use critical-plane
gauges alone. Both use the same private observer, original material-state
binding, CMB source checks and actual native stress/total-mechanical-strain
tensors. Plane projection and SWT inversion remain in `fs-material`; no
alternative constitutive or motion solver is introduced.

An explicit cycle-window restart samples the current accepted native solution
and discards only the preceding critical-plane window and its optional result. It permits
mechanical conditioning before the representative observation, preserving
displacement, velocity, plastic history and fixed-plane fatigue counters.
A closed stress/strain window is an admitted empirical input, not proof that
the evolving constitutive state has stabilized.

Ordinary advances retain every observation across processing chunks. Explicit
`advance_cycle` operations close the requested observation on a private native
candidate before publishing any motion or fatigue state. Nonclosure, nested
strain cycles, unsupported life, a late gauge failure and cancellation preserve
the previously accepted state. Finalization closes fixed and critical gauges
together. Critical-plane results are exposed separately from the existing
fixed-plane summaries.

This fixed-temperature, small-strain tensile-mode SWT path admits only resolved
representative cycles and declared candidate planes. Decaying free vibration,
virgin plastic transients or ratcheting paths that do not close cannot silently
become periodic-life estimates. Plane/time refinement and calibration
applicability remain the caller's obligations; there is no multiaxial
variable-amplitude counting, thermal-fatigue transfer, crack growth or
fatigue-to-stiffness/heat coupling in this extension.

## Conservative gravity on the current finite-strain configuration

With cellwise reference density, `hyper3::gravity` integrates a stationary
world-field potential, acceleration and analytic spatial gradient at four
current P1 tetrahedral quadrature points.
`PreparedHyperTetProblem::gravity_load` retains the original reference geometry;
the element-ordered densities are the same fixed reference densities consumed
by native `reference_mass_matrix`. Current volume never creates or removes mass.
Force, external-force tangent and potential use the same degree-two quadrature.
Uniform and affine acceleration reproduce the full consistent reference-mass
projection, including prescribed coordinates and material interfaces.

The problem and prepared `gravity_load_with_nodal_density` entries admit one
positive P1 reference density per original node. Their eight positive cubic
points (vertices of weight V/40 and face centroids of weight 9V/40) integrate
rho0*N_a*N_b exactly. Force, potential and analytic tangent share those fixed
reference mass weights, preserving graded center of mass and inertia moments.
Under admitted positive-J deformation this represents rho_current = rho0/det F: changed current
volume does not change inventory or replace the retained field by cell means.
Globally uniform nodal input retains the original four-point arithmetic.
`solve_with_nodal_density_gravity` and its combined pressure entry use the same
native Newton, line search and reaction owner; the bound mass weights are
prepared once. Current vertices and face centroids must lie in the field domain.

`HyperTetProblem::solve_with_gravity` uses that owner at every Newton and trial
configuration, subtracts its analytic tangent from internal stiffness, includes
its potential in conservative globalization, and recovers reactions with all
supported weight. The existing load factor ramps gravity. The combined pressure
entry preserves current-area pressure and its conservative/nonsymmetric split.
If a large absolute gravity datum obscures a small potential increment, an
actual residual decrease is required; equilibrium tolerance is not relaxed.

Callbacks must supply derivatives of one stationary scalar potential. Finite
values and gradient symmetry are checked; derivative consistency remains the
field owner's contract. Domain exits refuse or backtrack, never clamp. This is
an external field, not self-gravity, moving-source dynamics or a new material
law. Smooth nonlinear fields require spatial refinement; no quadrature error
certificate, singular-source resolution or mesh-independent result is inferred.
Focused G0/G1/G3 regressions in `tests/hyper3_gravity.rs` exercise heterogeneous
native mass, finite-difference derivatives, rotated geometry, analytic finite
stretch with pressure and support reactions, potential-datum invariance, and a
three-rung nonlinear quadrature ladder. These are synthetic numerical checks.
The nodal-density unit regressions additionally compare all cubic simplex
moments, actual shared reference-mass action, deformed gravity force/energy/
tangent and center-of-mass torque, nonlinear discrete-potential derivatives,
and scalar-path bits. `tests/hyper3_gravity/nodal_density.rs` adds the native
finite-stretch/pressure/reaction and public-admission regressions. The bounded
native source validation executes the three new assembly regressions and the
three shared-mass regressions; the full Newton integration regressions remain
unexecuted in the current dependency-constrained environment.

## Fatigue on source-bound inertial J2 motion

`plastic3::dynamics::fatigue::J2DynamicFatigueRun` owns accepted physical time,
displacement, velocity, native plastic history and bounded fatigue observation
while borrowing the original source-resolved J2 problem. It initializes through
the real zero-load solve, admits explicit velocity compatible with the original
supports, and freezes reference density from each original material bundle.
Every strictly increasing endpoint calls the existing `advance_backward_euler`
with actual endpoint loads. Native stress, reactions and inertial energy remain
inspectable; there is no alternate constitutive or dynamic solve.

Single steps and requested endpoint blocks publish only after every mechanical
increment and every fatigue gauge accepts. A late gauge failure, invalid time,
nonconvergence, cancellation or original strain-domain refusal leaves all prior
motion, plastic state, cycles and damage unchanged. Checkpoint clones retain
native motion and unresolved fatigue cycles together. Only the latest native
solution/step and bounded gauge residues are retained.

Quasistatic `J2FatiguePath` and inertial motion share the same private observer.
Uncorrected strain-life uses the scalar counter; named Morrow/SWT policies feed
actual simultaneous signed stress and total mechanical strain to the paired
counter, including its coincident-extrema and source-domain checks. A sourced
CMB curve must match the original material card, named state, external query
conventions and modulus claim. Caller-declared curves keep their explicit
identity/chemistry/process/fixed-temperature/modulus admission.

Backward-Euler numerical dissipation remains separately reported in the native
inertia receipt; it is not material heat or fatigue damage. Physical timestep
refinement and resolved turning points remain necessary. Miner accumulation
does not change stiffness, crack topology or thermodynamic state. The six native
regressions in `tests/plastic_dynamic_fatigue.rs` compare unforced vibration and
time refinement, cyclic J2/Prager response and inertial support reactions, named
mean-stress policies, continuation, source admission and atomic rollback.

## Strain-life from the original cyclic J2 constitutive path

`plastic3::fatigue::J2FatiguePath` owns the latest native solution and one
scalar or paired strain-life counter per fixed-plane gauge while borrowing the original
`SmallStrainTetProblem<ResolvedJ2Plasticity>`. The real virgin zero-load solve
is the first observation. Every original absolute load target goes through
`solve_from`; its plastic strain, accumulated isotropic history and Prager
backstress continue from the preceding accepted state. The observer projects
the returned tensor-Voigt total mechanical strain and signed normal stress.
It does not substitute nondecreasing equivalent plastic strain, unsigned von
Mises stress or a second constitutive approximation.

Each calibration must match its original resolved material-state identity,
chemistry, process, fixed temperature and Young's modulus. Original material
receipts remain inspectable through the path. The adapter accepts only the
source-resolved J2 owner, which has no thermal eigenstrain; thermal-wrapped
laws and changing-temperature histories require a separately admitted transfer.
Every accepted element, including ungauged ones, must remain within the native
law's finite-state and componentwise small-strain envelope.

`advance(count)` observes every intermediate target and commits the whole
requested block only after mechanics and all gauges accept. Cancellation,
nonconvergence, a later gauge refusal or an unsupported strain leaves the prior
mechanical and fatigue histories intact. Clones retain the complete current
state and unresolved cycles, without storing the full response trace. `finish`
requires completion of the declared observation and includes final half cycles
once. Processing windows do not create artificial cycles.

The CMB model and selected Uncorrected, Morrow or Smith-Watson-Topper policy
retain their empirical, fixed-plane scope. No stiffness, heat, fracture, failure topology, ratcheting
model or validated service life is inferred from Miner damage. The five focused
G0/G1/G3/G4 tests in `tests/plastic_fatigue.rs` compare the actual cyclic J2 coupon
with independent one-dimensional Prager equations and an analytical strain-life
damage sum, and check continuation, source admission, budgets and rollback.
Their data are synthetic numerical references, not held-out material evidence.


## Fatigue from accepted inertial stress histories

`loads3::transient::fatigue::TetTransientFatigueRun3` owns the original sparse
Newmark run and the same fixed-plane rainflow/Miner observer used by
`fatigue::TetFatiguePath`. `TetTransientRun3::with_fatigue` admits only step zero,
requires an explicit retained thermal field, checks every gauge's material-state
identity and calibration temperature, and samples the actual initial stress.
Every successful stationary- or moving-support advance recovers native endpoint
stress and samples every gauge once. Inertia, support motion, damping, thermal
eigenstrain, reactions and work remain in their original mechanical owners;
zero applied load does not imply zero vibration stress or zero fatigue.

The wrapper exposes no mutable run, external stress injection or material/thermal
reconfiguration. A late gauge refusal, stress-domain failure or cancellation
commits neither native mechanics nor any fatigue counter. Cloning retains both
accepted mechanics and unresolved cycles, with shared immutable matrices and
bounded state/residue copies. Read-only stress recovery adds no sample.
`finish` closes endpoint half cycles exactly once and prevents further steps;
processing chunks must not finalize their intermediate residues.

This is an elastic, fixed-temperature, fixed-plane empirical stress-life
estimate. Caller-supplied S-N/process/environment applicability is not inferred
from matching identity/temperature. Time and spatial refinement remain necessary
to resolve physical stress extrema. There is no multiaxial criterion, plastic
strain-life, changing-temperature law, remesh transfer, crack propagation,
stiffness degradation, damage-derived heat or service-life validation claim.

G0/G1/G3/G4 regressions in `tests/transient_fatigue.rs` cover unforced oscillator
stress and independently accumulated damage, a three-rung physical-period
reference, exact native mechanics/reactions, moving supports, checkpoint/chunk
continuation, complete rollback and source/history refusals. All fatigue data
in those tests is explicitly synthetic; the regressions are implementation
checks, not experimental material calibration.

## Finite-strain traction-free plane stress

`HyperProblem::solve_plane_stress(PlaneStressSettings, &Cx)` runs the original
load-stepped Newton solver on P1 triangles and Q1 quadrilaterals while the native
`fs-material` energy relaxes all three transverse nominal tractions at each
quadrature point. Thickness contraction and anisotropic transverse shear are
solved, not prescribed. The in-plane force and Schur-complement tangent derive
from the same reduced energy; supports, dead tractions and reactions retain the
original reference-thickness convention. `plane_stress_states` exposes the
actual relaxed 3-D deformation and local equilibrium residual.

Reference gradients and quadrature are prepared once per solve. Energy-only
line searches use the material's three-direction local Newton path without
forming full constitutive Hessians. The global solver checks the true linear
residual, polls cancellation between quadrature points and bounded Krylov
chunks, and refuses nonconvergence, folded/inverted elements, invalid controls
and material-domain violations without returning a partial solution.

G0/G1/G3 tests in `tests/hyper_plane_stress.rs` cover an independently derived
finite Neo-Hookean uniaxial traction solve with lateral/thickness contraction,
nodal energy/force/tangent consistency on P1 and skew Q1 elements with tilted
fibres, rigid-motion invariance, and domain/cancellation refusals. This is a
thin in-plane continuum: no shell bending, wrinkling, global injectivity,
fully incompressible formulation, general stability or convergence certificate
is implied. Existing `solve()` and continuation retain plane-strain semantics.
`cargo run -p fs-solid --example hyperelastic_plane_stress` solves a thin coupon
under nominal traction and reports displacement, relaxed thickness and support
reaction. Its constitutive coefficients are explicitly authored model inputs.

## Exactly incompressible finite-strain plane stress

`HyperProblem::solve_incompressible_plane_stress(PlaneStressSettings, &Cx)`
uses the same P1/Q1 geometry, integration, Newton solver and boundary loads as
the finite-bulk plane-stress path. At each point the native material reduction
enforces `det F = 1`, relaxes transverse shear and recovers the traction-free
pressure multiplier. Its reduced energy, forces and tangent do not use a bulk
penalty. The stored constitutive bulk coefficient can span finite scales
without changing this explicitly selected incompressible problem.

`incompressible_plane_stress_states` recovers physical first Piola stress,
volume-preserving 3-D deformation and the original material remainder's
pressure multiplier. This multiplier's gauge is not mean Cauchy pressure.
The corresponding potential and residual/tangent APIs retain reference-thickness
work conjugacy and reaction recovery. Original material stretch bounds apply to
the actual relaxed thickness; an allowed in-plane stretch does not override them.

G0/G1 tests in `tests/hyper_plane_stress.rs` independently recover incompressible
Neo-Hookean uniaxial contraction, stress, energy and support reaction on triangles
and quads; show bulk-coefficient independence through `1e300`; differentiate
assembled reduced energy/forces for tilted fibres on skew elements; and check
thickness-domain/cancellation refusals. The existing example accepts
`--incompressible`. This is the exact thin-continuum volume constraint, not a
general 3-D mixed formulation, shell bending model or arbitrary-mesh stability
certificate.

## Finite-strain Q1 mean dilatation

`HyperProblem::solve_mean_dilatation(&Cx)` adds an opt-in plane-strain Q1
formulation for quasi-incompressible compression. Its condensed element energy
integrates the original material remainder at all four Gauss points and the
original bulk energy at reference-volume-mean det F. Eight-DOF nested automatic
differentiation supplies the complete residual/Hessian, including geometric and
cross-quadrature terms. NH, MR, Ogden and full-invariant HGO retain their original
material split; no alternate constitutive coefficients are substituted.

Reference geometry is prepared once. The original load stepping/Newton owner
is reused, with bounded cancellation and an independently recomputed linear
residual. `mean_dilatation_states` distinguishes local det F, mean det F and
tension-positive condensed bulk Cauchy stress from the full constitutive stress.
Reference/current corners and every local material state must remain admissible;
an admissible mean cannot hide local inversion or a stretch-domain violation.

Focused G0/G1 regressions in `tests/mean_dilatation.rs` cover energy derivatives,
analytic finite compression/reactions on two meshes, an independent bulk-stiffness
identity, a nonuniform compressed-block reference, objectivity and refusal.
This mitigates the finite bulk penalty's excess local constraints. It supplies
no 3D/P1/TDNNS, fully incompressible, hourglass-stability, continuum-convergence
or contact guarantee. Existing `solve()` and continuation retain standard behavior.

## Static compression by spatial follower pressure

`HyperTetProblem::solve_with_face_pressures` consumes the existing canonical
face-pressure field in native finite-deformation equilibrium. Explicit zero
values select unloaded faces; positive values act inward on CURRENT area.
Loads, prescribed displacements and pressures share the original load ramp,
constitutive law, Newton assembly and support-reaction owner. The original
scalar uniform-pressure entry point retains its conservative solver behavior.

A partially loaded boundary generally has a nonsymmetric follower tangent and
path-dependent work. The new path retains the full analytic pressure Jacobian,
uses cancellable FGMRES with ILU(0) or absolute-diagonal preconditioning, and
globalizes against the force residual. It does not manufacture a pressure-volume
potential. The explicit linear iteration budget and independently recomputed
linear residual apply. All-zero face fields use the dead-load path exactly.

G1 tests in `tests/hyper3.rs` reproduce independent Neo-Hookean uniaxial
compression on two meshes, including lateral expansion, current-area reaction,
volume and stored energy. `tests/hyper3_pressure.rs` exercises a compressed,
sheared equilibrium with an actually nonsymmetric free tangent, reaction loads,
uniform-field equivalence, exact zero-field identity, invalid face coverage,
nonfinite values, cancellation and an exhausted linear budget. The existing
P1 volumetric-locking and stability limitations remain. Face loads are supplied
data, not a fluid equation of state or contact/platen law; no mesh-convergence
claim follows from an exactly representable affine test.

## Spatial follower pressure on finite-deformation exterior faces

`hyper3::pressure::TetPressure::PerFace` supplies finite signed pressures in the
canonical order of `ClosedTetPressureSurface::faces`. Explicit zeros select
unloaded faces. The original prepared geometry and force/tangent kernel serve
both a scalar uniform pressure and a face field. Positive pressure acts inward
on current face area; a patch keeps its own resultant, torque and generally
nonsymmetric external-force Jacobian. It is not assigned a volume potential.

Held face pressures have exact work along a straight nodal-displacement segment
to rounding, using Simpson integration of the quadratic P1 area vector. Common
translation work on an open patch remains present. Curved physical trajectories
require time refinement; this is not path independence or an injectivity proof.
Uniform pressure retains its original exact endpoint-volume work and static
solver path. Missing/extra/nonfinite face values and cancellation refuse before
returning a load. G0/G1/G3/G4 tests check localized outward orientation, independent
finite-difference tangents, path-dependent analytic work, uniform compatibility
and refusal. Dynamic consumption belongs to `fs-couple::finite_solid`.

## Finite-deformation static follower pressure

`hyper3::HyperTetProblem::solve_with_pressure` couples dead referential loads
and prescribed displacement with uniform follower pressure on the exact matching
complete closed boundary. The existing Newton owner evaluates current-area
forces, the analytic pressure tangent and conservative pressure-volume work at
every accepted or trial state. Force and volume share one geometry evaluation;
Armijo compares energy/work increments. Pressure ramps with load factor, and
support reactions include current pressure and applied dead loads. Positive
pressure compresses; negative pressure denotes outward traction. Invalid trial
geometry backtracks, cancellation aborts, and nonfinite reactions refuse.

This does not add an equation of state, open-patch loading, nonuniform pressure,
mixed incompressible elements or limit-point continuation. The focused static
tests in `tests/hyper3_pressure.rs` compare finite hydrostatic compression and
expansion, frame transformations, fixed-area traction differences, zero-pressure
equivalence and finite-deformation support reactions with independent formulas.
These synthetic regressions require native execution; they do not establish
experimental material calibration or mesh convergence.

## Prepared finite-deformation reference geometry

`hyper3::prepared::PreparedHyperTetProblem` retains admitted reference shape
gradients and volumes while immutably borrowing the original mesh and material
field. Repeated force, energy and tangent evaluations use the same native
constitutive and assembly kernels; current displacement, positive Jacobian,
material-domain and finite-response checks remain active. New interval dead
loads are validated on every call and use the original reference volumes.
Preparation takes O(elements) storage and removes repeated reference geometry
construction; it does not change the P1 discretization or mitigate volumetric
locking. The finite-solid dynamic adapter shares this preparation across
candidate steps. G0/G1 checks cover changing finite stretches against closed-form
Neo-Hookean energy/force/tangent, exact native assembly, new load admission,
inversion, invalid dimensions, budgets and cancellation.

### Objective energy/momentum force on the prepared mesh

`PreparedHyperTetProblem::energy_momentum_internal_force` and
`energy_momentum_tangent` evaluate an objective strain-space discrete gradient
between two configurations on the original P1 mesh. The mean right
Cauchy--Green metric is `Cbar=(C0+C1)/2`, represented by
`F*=chol(Cbar)^T`; it is not the metric of the mean deformation gradient.
The original material supplies energy, Piola stress, tangent and endpoint-domain
admission. A symmetric second-Piola work correction, multiplied by midpoint F,
retains partition-of-unity force balance, midpoint internal torque balance and
the potential-work identity. No physical damping or constitutive replacement
is introduced.

The analytic generally nonsymmetric 9-by-9 endpoint tangent differentiates the
same Cholesky representative, material pullback and work quotient. A consumer
using midpoint displacement must multiply its action by two. Near
`||Delta C||^2 <= 64 epsilon ||Cbar||^2`, the stress takes the continuous limit;
measured nonlinear and energy tolerances remain authoritative. Original
reference geometry, material/reference adapters and consistent mass are retained.

Four exact native stress/tangent tests pass, covering work, internal force and
torque, all nine tangent directions, objectivity, endpoint reversal, the
continuous limit and source-domain refusal. Two 250-step native free-tumbling
P1 tetrahedron trajectories with the original generalized-alpha/Newton-FGMRES
owner retain energy and both momenta near roundoff for the admitted fixtures.
Full prepared/Cx/material-card integration regressions are authored but remain
unexecuted because the required sibling runtime dependencies are unavailable.

### Uniform closed pressure in the objective discrete gradient

`energy_momentum_force_with_uniform_pressure` and its endpoint tangent combine
the original frozen material with the pressure potential `p det(F)` before
the same objective strain-space algorithm. Positive pressure is inward. The
native prepared closed-pressure surface owns reference matching and finite-load
admission; its original endpoint-volume work remains independently measurable.
The analytic cofactor stress and Hessian preserve the original material/source
domain. Zero pressure delegates exactly to the unpressurized algorithm.

The symmetric algorithmic second-Piola stress closes pressure work and midpoint
torque together. Exact straight-segment averaging of face forces closes work
alone and can produce a midpoint torque defect; a nontrivial regression exposes
that difference. This is a finite-step algorithmic volume load. It approaches
native boundary traction as endpoints coalesce, but its finite-step nodal field
need not equal midpoint traction or cancel on interior nodes. Nonuniform face
pressure remains outside this energy/momentum entry.

All six production stress/tangent groups pass, including signed pressure,
endpoint tangent directions, objectivity/reversal, native face-traction limits
and the work-versus-torque counterexample. Native time integration also closes
energy and both momenta and shows resolved second-order refinement. Public
prepared/Cx integration remains syntax/static reviewed, not an executed Cargo
qualification in the dependency-incomplete environment.

## Source-pack mechanical coupons and nonthermal dead loads

`linear3::eigenstrain` admits explicit sorted unique dead nodal forces in
assembly, equilibrium solve and input-bound stress recovery. It shares the
existing force validator and physical solver, while retaining generic initial
strain and requiring no fictitious temperature. Fixed-DOF loads remain in the
full force vector for support reactions. An empty force list preserves the
original eigenstrain-only identity and response exactly; changed loads refuse
recovery against an older solution.

`examples/material_coupon.rs` reads a normalized material pack and an explicit
bounded TSV scenario, then runs native tetrahedral face-load integration,
equilibrium and stress recovery. It resolves independent isotropic parameter
pairs or complete source-bound FPL orthotropy with explicit L/R/T orientation.
Source property roles, quantity-typed coordinates, source/derivation/orientation
identities and strain ceilings remain explicit. Templates exercise aluminum
5005, C52100 bronze and loblolly pine without inline material presets.

The aligned uniform traction fixture allows Poisson contraction and has an
independent axial compliance reference. Output reports source-bound mass,
displacement, applied/support forces, strain energy and residuals. Its
load-displacement product is F dot u, not a loading-history work integral.
Rotated general anisotropy is a coarse FE result. Neither a source data lookup
nor this linear model supplies strength, creep, fracture or mesh-error bounds.

G0/G1 regressions cover combined initial strain and dead force, exact affine
displacement/stress/energy, fixed-support reactions, empty-load equivalence,
malformed forces, stale-load refusal, signed coupon loads and strain admission.

## Material-bound spatial gravity loads

`loads3::history::VolumeHistory3::FixedAcceleration` accepts uniform, nodal P1,
or element-local P1 accelerations through `fields::AccelerationField3`. Values
are world-frame m/s² on the retained reference mesh. The borrowed solve problem
supplies the authoritative element density. Full consistent P1 mass acts before
support elimination, preserving interface density jumps, restrained-node loads,
resultants and first moments. The same scalar mass kernel supplies the original
elastic inertia. Typed uniform acceleration histories use this same integration
owner without allocating a second per-element force-density field.

The accepted sample enters the existing static, stress/reaction, support-history,
and sparse transient owners. G0/G1/G3 tests in `load_history3` compare against
full M*g, heterogeneous material mass, an analytic gravitational-potential
gradient, rotated reactions, initial dynamic acceleration and refusal behavior.
This is prescribed reference-sampled loading, not self-gravity, follower gravity,
or an accelerating-coordinate model. An inverse-square field sampled at vertices
retains P1 spatial approximation error and requires explicit resampling if its
physical owner changes position or time.

`TetLoadGeometry3::assemble_acceleration_with_nodal_reference_density` and
`TetLoadHistory3::at_with_nodal_reference_density` explicitly bind a positive
P1 reference-density field instead. They use the same exact cubic mass kernel
as graded elastic inertia, preserving density-weighted gravity moments and
the product of spatially varying density and acceleration. Full supported-node
forces remain available to the original static stress/reaction owner; surface
and point forces are superposed once. Uniform fields retain the original bits.
This entry projects a held force for static use; it does not replace the mass
in the separate `loads3` transient constructors. Three focused G1/G3 regressions
in `load_history3` cover independent barycentric integration, graded `M*a`,
weight/torque/reactions, scalar parity and refusals.

## Moving force/couple ports on the native shell trajectory

`moving_load::MovingShellModel` binds the original `shell::ShellPlate` geometry,
six-DOF spatial operators and stationary support map once. A point or explicit
patch quadrature transfers Cartesian forces and free couples with the transpose
of its P1 material-point motion, including `offset cross force`. Full restrained
rows survive for reactions. The same rows sample displacement, director rotation
and material velocity; moving a load's attachment across triangles does not add
the convective traversal velocity to that material velocity or its power.

`MovingShellMotion` initializes acceleration from the original mass equilibrium
and reuses `shell_time::PreparedNewmark` with beta=1/4 and gamma=1/2. Retained old
and supplied new endpoint loads provide trapezoidal force and reference-moment
impulses. Support and damping exchanges, supplied work, and measured momentum
and energy residuals stay separate. Mass-proportional Rayleigh damping is a
momentum sink; it is not silently attributed to internal material heat. Force
projection is linear in point count, factors are reused, reaction evaluation
visits only restrained rows, and impulse reductions need no extra dense vectors.
The caller's point/step/horizon bounds and cancellation checks precede atomic
publication; duplicate, stale or future endpoint timestamps refuse.

This is the existing flat, small-rotation, estimate-only shell surrogate. Its
angular inventory and wrench moments use original reference positions and the
linearized director field. It is not a finite-rotation continuum shell,
follower-force/contact solver, automatic patch tracker, moving-support law, or
modal reduction. Time-varying load paths require temporal refinement. The
physical subset of bead `frankensim-b8bxd.9.3` is implemented; that broader
adapter's remaining features are not claimed complete. G0/G1/G3/G4 tests in
`tests/moving_load.rs` exercise independent wrench resultants, rigid motion,
reciprocal work, analytic accelerated translation, native damped/undamped moving
loads, support reactions, common-edge traversal, and atomic refusal/retry.

## Momentum-stable shell integration

`shell_time::step_newmark` solves the effective mass equation for acceleration,
then adds complete displacement and velocity increments. Absolute positions are
never differenced to reconstruct acceleration; large rigid drift and short time
steps therefore do not inject artificial momentum. The original `fs-la`
partial-pivot LU owns the solve, without a dimensional absolute pivot cutoff.
For average-acceleration parameters and initially equilibrated acceleration,
the full consistent mass carries the trapezoidal force impulse and the discrete
force-work balance includes internal viscous dissipation. Other Newmark
parameters do not inherit those conservation or unconditional stability claims.
The G0/G3 `shell_time` tests cover free drift, applied impulses, off-diagonal
mass, internal spring/damper reactions, work and physical coefficient scaling.

`shell_time::PreparedNewmark` borrows fixed operators and factors the effective
mass once through the same `fs-la::LuWorkspace` owner. Repeated steps reuse that
factorization with quadratic-cost matrix-vector products and triangular solves;
the one-shot API shares its complete-increment advancement and admission rules.
The prepared system must be rebuilt when operators, timestep or Newmark
parameters change. It can continue a restored checkpoint or another trajectory
without refactoring, with bit-identical results to fresh one-shot solves.
Additional `shell_time` tests exercise full nonproportional damping with
force-work balance, prepared/fresh equality over a forced trajectory, exact
checkpoint continuation after rebuilding factors, and second-order convergence
to a manufactured forced response.

## Purpose and layer

Layer: L3 (FLUX). The elasticity core (plan §8.2, bead tfz.13):
small-strain elasticity through finite-strain hyperelasticity, with a
MEASURED locking-free formulation, on both the body-fitted and the
CutFEM-on-SDF frontends. Constitutive laws live in fs-material; this
crate owns body-fitted kinematics, weak forms, boundary operators, and
the globalized Newton loop. Vector CutFEM assembly and solution are
owned canonically by fs-cutfem; fs-solid retains only a legacy
constructor facade.

The `linear3` module is the generic body-fitted three-dimensional
small-strain rung. It assembles a sparse stiffness/consistent-mass pencil from
P1 tetrahedra plus uniform or per-element `fs-material` state points. Its
constitutive seam is a finite, exactly symmetric, positive-definite 6x6 tensor
in an explicitly documented Mandel basis; isotropic and oriented-orthotropic
constructors lower into that shared representation. Object names and material
presets are deliberately absent: downstream vibration and vibroacoustic
consumers receive the same `(K,M)` artifact for any admitted
geometry/material composition. Anisotropy is therefore supported at the
operator level, but material-card resolution of engineering/crystal constants
and texture/orientation fields remains a separate evidence-bearing admission
step; absence of that data is not permission to assume isotropy.

`TetLinearElasticProblem::assemble_with_nodal_reference_density` admits an
explicit positive P1 density field on the original shared nodes. It produces
the same native sparse pencil for `fs-modal`, with exact cubic integration of
`rho_0 N_a N_b`; a cell mean preserves total mass but not these inertia moments.
The scalar element kernel is shared with finite-strain `P1ReferenceMass`.
Density controls reference mass while the original material field controls
stiffness. Uniform nodal input preserves the original mass bits; connectivity,
support reduction, quality, budgets and cancellation share the original
assembly path. Invalid or unrepresentable density/mass refuses before any
pencil is returned. Total inventory includes supported nodes. The density
field is explicit authored input, not material-source resolution or a coupled
thermal/phase refresh; a closed body's reference field must be retained.
Native modal regression moves density toward the free vertex at identical
total mass/stiffness and verifies the analytic eigenvalue ratio `2/3`.

`HyperTetProblem::reference_mass_matrix_with_nodal_density` carries that same
exact P1 density integration into finite-strain mechanics. It assembles in the
original node/element order, retains supported coordinates, and admits finite
positive nodal density and representable mass under the native mesh, budget
and cancellation checks. The matrix is reference inventory, held through
deformation. It does not refresh density from current volume or infer a
constitutive field from inertia. G0/G1/G3 regressions in
`tests/hyper3_reference_density.rs` cover cubic moments, centroid, kinetic
energy against independent quadrature, positive definiteness, scalar bit
parity, shared-vertex/permuted assembly and refusal boundaries.

`tet10::QuadraticTetProblem::assemble_with_reference_density` carries graded
reference inertia into the same native quadratic stiffness/mass pencil. Each
element supplies four positive corner densities, interpolated with P1 reference
barycentric weights, including on curved physical maps. Material interfaces may
be discontinuous. The ten quadratic displacement shapes retain the full
consistent mass; their negative corner lobes are not used to interpolate density
or to lump mass. Total inventory is integrated before constrained rows are
removed. `assemble_thermoelastic_with_reference_density` changes stiffness with
the admitted thermal material field while retaining that original density.
No temperature-dependent current density refreshes a closed body's mass.

Four mass Duffy points per direction integrate affine graded quadratic mass
exactly (degree five). A non-default profile on a curved map requires six
(degree eight including the cubic Jacobian determinant); this guarantees the
mass quadrature degree, not stiffness accuracy or global map injectivity.
Existing geometry, work budgets and cancellation still apply. A profile equal
to the original scalar material density preserves the existing mass bits.
Native regression checks every mass entry against independent barycentric
moments, translation/rotation inertia, positive definiteness, curved quadrature,
and frozen mass under changed thermal tensors. The actual generalized modal
solve verifies that relocating density at identical total mass and stiffness
increases a supported edge's mass by 25% and reduces its eigenvalues by 20%.

`hex::HexElasticProblem::assemble_with_reference_density` extends original
reference inertia to HEX20 and HEX27. Each element supplies eight positive
corner densities; nonnegative Q1 weights interpolate that reference field
while the original quadratic displacement shapes retain the full consistent
mass. Interfaces may be discontinuous. The resulting native sparse pencil
retains total inventory before support reduction, density-weighted first
moments and rotational inertia. Thermal assembly can change local stiffness
through `assemble_thermoelastic_with_reference_density` while the same graded
reference inventory remains fixed.

The explicit mass quadrature is independent of stiffness quadrature. Equal
rules reuse one prepared point array. Three Gauss points per direction
integrate affine graded quadratic mass exactly; other new-profile geometries
conservatively require six for coordinatewise degree ten, including the
quadratic map's Jacobian determinant. Exact nodal affine relations determine
that admission. These bounds certify mass polynomial integration, not global
map injectivity or stiffness accuracy. An identical scalar profile with the
original rule preserves prior mass bits; a different requested mass rule is
honored. Native regressions compare every mass entry with independent monomial
integrals on affine and curved maps, retain first/rotational moments, expose a
five-point underintegration case, and check exact supported eigenvalue shifts
at equal total mass and stiffness. Existing geometry, budget and cancellation
admission remains in force.

`linear3::TetElasticMaterial::acoustic_response` exposes the shared
`fs-material::elastic::acoustic` law using the exact native Mandel stiffness,
density and material identity. Source-resolved isotropy, oriented orthotropy and
full tensors therefore derive directional wave speeds/impedances and hydrostatic
compressibility/strain from the same coefficients as assembled mechanics. The
caller declares an explicit strain ceiling, as for stress recovery, without
extending source validity. This is a frozen small-strain model query; static
moduli do not become measured ultrasonic data. Native regressions connect
source-selected bulk/shear constants and oriented full-tensor states to these
derived responses while preserving original identities.

The shared `linear3::source` adapter resolves explicitly selected isotropic
Young/Poisson, Young/shear or bulk/shear source pairs, or oriented
orthotropic elastic constants, and independently selected isotropic or three-axis
integrated expansion. Both `fs-mbd-material` and spatial
`fs-couple::solid_enthalpy` use this one lowering into the original Mandel
tensor and eigenstrain owners. Current stiffness does not require current
density: the caller supplies conserved reference mass per reference volume.
Reference resolution can select density and the full elastic schema atomically.
Every public lowering checks a complete matching property schema before
accessing engineering constants. Source integration, compliance admission and
proper-frame transforms remain in their original owners.

All three isotropic source choices lower through
`fs-material::state_point::ResolvedIsotropicElasticModuli`; derived E/nu values are not
inserted into the card or relabeled as measurements. The selected pair, source
receipts and reference inventory remain bound to the native material identity.
The source owner's Mandel tensor carries the retained K/G or E/G moduli directly
into mechanical, modal and thermal assembly. It never reconstructs them through
a rounded Poisson ratio near incompressibility or the auxetic limit. Canonical
E/nu input preserves its original tensor arithmetic. The same direct K/G tensor
construction serves the acoustic owner; this numerical correction does not
remove the existing displacement elements' volumetric locking.
Changing temperature re-resolves the same selected source schema before
rebuilding stiffness, thermal stress and vibration inputs, while reference mass
and inertia remain unchanged. G1/G3 tests exercise equivalent source pairs,
temperature-dependent compliance, thermal stress, exact mass retention and
source/identity refusals through the existing tetrahedral solve. This does not
extend the separate reciprocal affine thermoelastic potential's E/nu derivative
law: that consumer explicitly refuses alternative pairs until its thermal
derivatives are supported.

TET10 and HEX20/HEX27 expose `assemble_thermal_load` using the same
isoparametric geometry, Mandel strain operator, elastic tensors, constraints
and fully integrated stiffness rule as their mechanical assembly. Existing
`TetThermalStrainField` values supply total eigenstrain, uniform or constant per
element, bound to the exact elastic material state. The sign is
`K u = f_external + integral B^T C epsilon_free dV`. The existing `TetThermalLoad`
carrier retains both full support rows and the canonical reduced free-DOF order,
so the existing solver can compute compatible thermal displacement. Load identity
binds geometry, connectivity, source states, actual elastic coefficients,
quadrature and supports. Both stiffness and mass quadrature locations retain the
original geometry-admission checks. Cancellation and numerical/material refusal
publish no partial load. This does not add a thermal transport solver, within-cell
temperature interpolation, geometric prestress, or
an accuracy certificate for curved-map quadrature. `tests/quadratic_thermal.rs`
checks uniform anisotropic eigenstrain including Mandel shear, resultant/moment
balance, thermal work, actual free-expansion solves on all three element families,
source/state mismatch and constrained loads against analytic references.

The same quadratic APIs also accept `QuadraticThermalStrainField::AtQuadrature`.
Each element supplies one source-integrated Mandel eigenstrain at each of its
original stiffness quadrature points, in explicit element-major/rule order.
The field checks both the points-per-element value and complete sample count;
every sample must bind the exact element elastic material identity. The
elementwise input path retains its original arithmetic, load identity and
one-constitutive-contraction-per-element behavior. No temperature interpolation,
constant-alpha approximation or averaging of nonlinear expansion is inserted.
Elastic tensors remain frozen per element. Changing quadrature requires new
samples, and source-data and quadrature error remain the caller's obligations.

Pointwise thermal loading is shared by equilibrium, stress/energy recovery and
the initial-stress vibration pencil. A zero mean strain does not erase stored
energy or the correlation between stress and shape gradients. G1 tests in
`tests/quadratic_thermal_gradient.rs` verify exact compatible thermal-gradient
bending on TET10/HEX20/HEX27, the analytic restrained energy of a zero-mean
through-thickness gradient, and its signed geometric-stiffness quadratic form.
Count/rule and late source mismatches refuse without publishing a partial load.

`resolve_nodal_thermal_strains` connects actual nodal reference/current absolute
temperatures to those same TET10/HEX20/HEX27 mechanics. It uses each family's
native geometry/shapes and stiffness quadrature, interpolates temperature first,
then calls the existing source-integrated isotropic or oriented orthotropic
expansion owner. Nonconstant alpha is integrated over each local temperature
path; nodal eigenstrains or endpoint alpha times temperature change are not
substitutes. The original resolved expansion reference pins the card, selected
claims, nonthermal coordinates and quantity conventions. A source selection is
prepared once per element, and uniform nodal fields resolve their integral once
per element. Returned original source receipts and accepted strains remain in
native element/quadrature order and borrow directly into load/response/tangent
methods. There is no second thermal law, shape kernel, reference-mass update or
implicit elastic-tensor refresh. Invalid counts, source mismatch, unsupported
temperature paths, nonpositive/nonfinite temperatures and higher-order
interpolation overshoot outside source support refuse the whole result; no
clamping hides those exits. Cancellation is checked per element and point.
Temperature transfer from a different thermal mesh/space remains explicit,
and this prescribed-temperature bridge does not claim caloric or reciprocal
energy closure. `tests/quadratic_temperature_source.rs` checks nonconstant-alpha
integration with a nonuniform stress-free reference, real gradient bending and
joint geometry/material-frame rotation, exact constant-field load equivalence,
source/count/domain refusal and positive-nodal temperature overshoot.

`resolve_nodal_thermoelasticity` explicitly couples source-resolved local
elasticity to the same nodal-temperature bridge. `QuadraticThermoelasticSource`
binds the original elastic and expansion references to one card and physical
context, with independent declared tensor symmetries and orientations. Pinned
current elastic properties are resolved before assembly; current density is
never queried. Uniform element temperatures resolve each source law once.
`ThermomechanicalAtQuadrature` carries matching local materials and integrated
eigenstrains at the native stiffness points. The exact local C(T) enters the
actual displacement stiffness, thermal load, recovered stress/energy and signed
initial-stress stiffness. `assemble_thermoelastic` exposes that same stiffness
and original consistent mass. Every local density must equal the original
element density; heating cannot redefine the reference inventory. Local sample
counts, material/strain identity and source support are checked, and the existing
load identity includes every local tensor's identity and coefficient bits.
Legacy elementwise and strain-only field paths retain their previous behavior.

This is the mechanical tangent at prescribed, frozen temperature, not a
coupled temperature derivative, transient thermoelastic feedback or a caloric
energy closure. Source and quadrature validity remain explicit: variable C(T)
and prestress can require a higher integration rule. The independent polynomial
oracles in `tests/quadratic_thermoelastic_material.rs` check stiffness, thermal
work, stored energy, signed prestress, compatible thermal bending and unchanged
mass on all three quadratic families. They also check a stiffness-domain exit
despite admitted expansion data and late local material/count/density refusals.

The same TET10/HEX20/HEX27 problems expose `solve_thermal_response` for source-bound
thermal equilibrium on a frozen reference mesh. The existing `fs-solver::CgState`
and Jacobi preconditioner solve the assembled stiffness with bounded
iteration/cancellation windows and an independently recomputed true Euclidean
residual. Every connected component must have all six rigid motions restrained.
Fully fixed bodies have an exact empty solve, while their nonzero thermal stress,
elastic energy and support force are retained.

Recovery integrates actual `sigma = C (B u - epsilon_free)` and
`epsilon_elastic:sigma/2` at the same stiffness quadrature points, checking both
total and elastic Mandel strain against the caller's finite positive ceiling.
`QuadraticThermalResponse` contains compatible nodal displacement, explicitly
named volume-mean element strains/stresses, integrated energy, and full nodal
internal forces. Those forces are reactions **on the solid** at fixed DOFs and
unmodified equilibrium residuals at free DOFs, with no external mechanical load.
Element means do not estimate pointwise or peak stress; energy integrates the
varying strain field before averaging. The immutable result binds the thermal
load, actual material coefficients, mesh frame, controls and observations;
failure publishes no partial solution or field. Geometry, mass and density are
unchanged. This does not add live thermoelastic transport, geometric stiffness,
finite-strain behavior, stress extrapolation or quadrature-error certification.
`tests/quadratic_thermal_response.rs` exercises free expansion, a piecewise
nonuniform expansion field, clamped stress/energy/support reactions, repeatable
results and refusals for rigid modes, controls and applicability limits.

Quadratic problems also expose `assemble_thermal_load_with_nodal_forces` and
`solve_thermal_response_with_nodal_forces`. Their dead nodal forces use the same
strictly ordered, unique, finite full-DOF validation and combined-load identity
owner as P1. Thermal and mechanical forces enter one stiffness solve. Actual
stress still integrates `C(Bu - epsilon_free)`; `nodal_internal_force_n` retains
that full internal force, while `nodal_force_balance_n` subtracts applied dead
nodal forces. The latter gives reactions on the solid at fixed DOFs and unchanged
residuals at free DOFs, including applied loads directly on supports. The result
retains both the original thermal-load identity and the complete load identity.
`tests/quadratic_thermal_loads.rs` checks analytic thermal expansion plus end
traction on all three families, support-load accounting, and malformed loads.

`linear3::visco` (IRP M1.3-W2a, in progress) prepares P1, TET10, HEX20 or
HEX27 elements with caller-supplied, mesh-frame `fs-material::TensorProny`
laws and density.
It reuses the existing elastic mesh admission, equilibrium stiffness and
consistent mass, and the shared Mandel B operator. Complex stiffness retains
the full storage and loss tensors, including off-diagonal couplings. Every
frequency query must lie in every element's declared law band. Empty branches
reproduce the canonical elastic stiffness bitwise. Positive quadrature weights
propagate pointwise relative quadratic-form bounds to the assembled forms;
this alone is not a frequency-dependent modal-eigenvalue certificate.


`ViscoElements::prepare_harmonic_response` connects the complete material
spectrum and the original consistent mass to the shared `fs-modal` harmonic
owner. At each positive requested frequency it solves
`(K'(omega) - omega^2 M - i K''(omega)) u = f` under `exp(-i omega t)`;
the exact single-frequency identity `C_eff = K''/omega` does not replace
hereditary histories with memoryless damping. The full coupled tensor entries
remain present, and every element law must admit the frequency. Preparation
reuses one dense factor for multiple force phasors: O(n^3) factorization,
O(n^2) retained storage and solves, plus original quadrature strain checks.
Before publication, the scaled 2x2 Gram problem of real/imaginary Mandel
strain admits the maximum norm over the complete cycle at each source point.
Outputs retain displacement/velocity phasors, numerical backward error,
independent mean input work and constitutive heat, and their signed defect.
No dispersive stored-energy formula is inferred from storage stiffness.
Budgets include the three dense input matrices and shared harmonic payload;
source sparse assembly retains its existing limits. Unrepresentable effective
loss, source-domain violations, cancellation and singular pencils refuse.
This bounded dense route does not claim spatial convergence, unsampled strain
extrema, transient response, moving supports or experimental validation.
`tests/visco_harmonic.rs` specifies analytic Prony/Kelvin response and heat,
nonproportional coupled tensor/consistent-mass phasors, complete-cycle strain
refusal, deterministic factor reuse, lossless response, budget/band and
cancellation checks (G0/G1/G3/G4).

`SolidMeshRef::Hex` retains the original hexahedral connectivity and positive
tensor-product quadrature (at least three points per axis). The original
`HexElasticProblem` admits topology and sampled Jacobians and assembles the
reference elastic pencil; the same hex shape gradients feed full-tensor
spectra, projected operators, constitutive histories and integrated heat.
This includes curved HEX20/HEX27 geometry without a tetrahedral remesh.
Jacobian and small-strain checks apply at declared integration points;
they do not certify unsampled geometry or strain extrema.

The required `ViscoStateBudget` counts six internal strain scalars per branch
and stiffness quadrature point, before history allocation. A history borrows
its immutable prepared mesh/laws and starts relaxed at zero displacement.
Prescribed full nodal displacement ramps advance the constitutive owner's
exact exponential update; B-transpose returns complete force rows, including
support reactions, and reduced forces use the canonical free-DOF ordering.
Stored energy, exact ramp work and dashpot heat are integrated in joules.
Work is the material owner's time-integrated work, not endpoint force times
displacement. Candidate histories commit only after all points and assembled
forces/energies succeed; cancellation and invalid/nonfinite inputs preserve
the previous state. Cx is checked per element/point, including state copying.
No damping clipping, diagonal loss approximation or secondary constitutive
law is introduced. Determinism is same-host at fixed element/point order.

`tests/visco_elements.rs` names G0 uniform-strain tensor/element energy,
elastic parity, 4096 counter-stream multiaxial steps, one-branch Zener response,
band/state-cap/invalid-history refusals; G1 harmonic endpoint force on three
time rungs (observed order 2 +/- 0.2); G3 rotation invariance; and a G4
cancellation-preserves-state check. These fixtures are synthetic software
checks and require native execution.
The HEX20/HEX27 fixtures additionally check rotated multiaxial spectra,
work-conjugate full forces, exact material work/heat integration, canonical
elastic mass/stiffness parity, quadrature/state-budget refusals, and curved
HEX27 storage/loss convergence against an analytic spatial integral.

`ViscoElements::quasistatic_operator(dt, cx)` prepares a reusable sparse
consistent tangent for one positive duration; `ViscoElementState::step_quasistatic`
solves endpoint nodal loads and prescribed support motion through the shared
Jacobi-preconditioned CG owner. The original full Mandel Kelvin/Prony tensors,
quadrature and accepted branch histories remain intact. Each step solves free
displacement increments, independently checks the original internal-force
residual, and returns complete displacements, support reactions and constitutive
storage/work/heat. Support motion contributes its actual integrated material
work. All connected components must pass the existing rigid-motion constraint
admission, and every original material band must include DC for sustained-load
creep. Failed solves, source-domain violations and cancellation preserve the
accepted displacement and memory. This is fixed-material quasistatic evolution
without inertia: each affine displacement ramp has exact constitutive
integration, while a continuously prescribed load trajectory requires time-step
convergence. DC admission does not certify that trajectory's transient spectrum.
`tests/visco_quasistatic.rs` checks analytic Maxwell/Kelvin creep convergence,
independent heat, nonproportional coupling, moving supports, tangent derivatives
and transactional refusals.

`ViscoElements::dynamic_operator(dt, config, cx)` prepares consistent-mass
midpoint dynamics on the same P1/TET10/HEX20/HEX27 tensor histories. Its
`start(history, full_velocity, cx)` retains an accepted constitutive preload;
`ViscoDynamicState::advance(mean_force, cx)` consumes the actual interval-mean
nodal load. The solve uses `M + dt² K_mean/2`, where `K_mean` differentiates
the original exact strain-ramp mean force. Geometry, tensor assembly and CG
are shared with the existing owners. Endpoint stresses remain observations;
their trapezoidal average does not replace the work-conjugate material force.
The accepted interval independently checks momentum and
`external_work = change(kinetic + material_storage) + physical_heat`.
No energy defect is relabeled as heat. Exact per-step work and loss survive
large cumulative counters. Explicit step/free-DOF/iteration budgets, source
strain limits, finite arithmetic, cancellation and energy allowances are
checked before publishing motion, memory, work or clock changes.

Supports are stationary at their accepted displacements and have zero velocity
and zero work. Full dynamic support reactions and moving-support inertial
couplings are not reconstructed from the reduced mass. Free bodies are admitted
through their original consistent mass, without artificial constraints. Every
material band must include DC; transient spectral admission and temporal
convergence remain caller obligations. The material ramp is exact, while the
coupled inertial trajectory is a time discretization. `tests/visco_dynamics.rs`
checks elastic energy conservation, analytical Kelvin-oscillator convergence,
independent heat, nonproportional memory, applied work, free-body motion and
transactional refusals. The existing element suite also carries a shared
TET10/HEX20/HEX27 kinetic/material-energy regression.

`ViscoElements::modes` solves the storage-dispersive ordered eigenproblem
`lambda(K'(omega), M) = omega^2` using the existing dense fs-modal eigensolver
and an explicit frequency bracket, iteration/tolerance budget, and a dense
reference dimension cap (at most 256 free DOFs). For positive equilibrium and
PSD Prony branches, `K'(omega)/omega^2` decreases in Loewner order, so ordered
positive-mode roots can be bisected without assuming fixed-point contraction.
Endpoints must bracket the selected mode inside the common material band;
missing roots, exhausted iterations, and an eigenvalue whose certified
quadratic (Kato–Temple) radius exceeds the budget refuse. The radius comes
from `fs_modal::recentred_dense_radii`, because a dense eigenvector's linear
residual stays near `ε·λ_max` on a stiff mesh and is no measure of the
eigenvalue; the linear residual is still reported.
The output retains the mass-normalized shape, full loss-energy quotient,
numerical frequency interval and nonlinear residual. Numerical bisection is
not an outward-rounded modal certificate, and the loss quotient does not
certify a diagonal damping approximation. A one-coordinate P1/TET10 Zener
fixture compares the self-consistent frequency with an independent analytic
quadratic (2e-10 relative), and deliberately shifts it by more than a factor
of two from the equilibrium result.

This module does not yet supply the fractional-target modal-loss ladder,
certified modal block reduction, valve-card integration, moving-support inertial
evolution, or material/specimen validation. Arbitrary
time-domain input is not a certification that its spectrum lies in the law's
band. Intrinsic material heat excludes support, contact and radiation losses.

`linear3::tet10` adds a separate isoparametric quadratic displacement/geometry
operator (IRP M1.3-W1a). `QuadraticTetProblem` consumes ten-node connectivity,
the existing uniform/per-element material tensor, fixed DOFs and assembly
budget, and a required `TetQuadrature`. Canonical order is four vertices, then
edges `(0,1),(1,2),(0,2),(0,3),(1,3),(2,3)`. Shared edges must share exactly one
midside node; a midside node cannot also be a corner or represent another edge.
The signed Jacobian must be finite, positive and above the scaled quality floor
at every stiffness and mass quadrature point. The P1 arithmetic, thermal loads,
thermal prestress and stress-recovery paths remain unchanged.

`linear3::hex` (IRP M1.3-W1b) adds the fully integrated isoparametric HEX20
serendipity and HEX27 tensor-Q2 families. They use the **same** const-generic
Mandel/consistent-mass/COO assembler extracted from TET10, the same material
field, fixed-DOF map, `TetElasticAssembly` carrier and fs-modal eigen path.
Canonical node coordinates and Gmsh-12/17 permutations live in L2 fs-mesh's
`structured_hex`; fs-solid now has a runtime dependency on that owner. The
historical `maximum_tetrahedra` budget counts hex elements too. Family and
connectivity cardinality must agree. Shared corners/edges/HEX27 face centres
have one global node per topological entity; duplicate cells and conflicting
node roles refuse. This checks conformity, not arbitrary cell intersections.

`HexQuadrature` requires 3..=12 Gauss-Legendre points per direction for both
operators. A 2x2x2 rule refuses `UnderIntegrated`: a test-only contraction
counts its HEX20 hourglass modes and cross-checks full-rule stiffness. Mass
is always consistent, SPD on admitted test meshes; HEX20 corner row sums are
negative and row-sum lumping is not available. Every stiffness/mass Gauss
point checks signed Jacobian and the declared scaled floor. General curved
maps still require a quadrature ladder; 3x3x3 is exact for affine mass/stiffness,
not a curved-map accuracy certificate. The common kernel retains TET10's
ordered local triangle arithmetic and polls Cx per element and quadrature point.

`tests/hex.rs` specifies G0 partition/Kronecker/gradient identities (1e-14),
six rigid modes, exact mass/first moments (1e-14), quadratic bending completeness
(1e-10), topology/Jacobian/underintegration/cap/cap+1/cancellation refusals, and
Gmsh permutation checks. HEX27 reproduces quadratic physical fields under the
tested trilinear distortion; HEX20's loss there is logged, not qualified away
(Arnold, Boffi & Falk 2002, https://arxiv.org/abs/math/0005036).
G1 roller-face P waves use an fs-mms OrderGate on three x-refinement rungs.
G2 square free-free bars (L/h=50, nu=.3) gate the first three bending pairs and
first axial mode to Timoshenko/Cowper and c/(2L), respectively (0.1 percent),
after a three-rung fs-mms-fitted Richardson frequency estimate <=0.02 percent.
The Timoshenko characteristic is solved by deterministic Newton and checked
against independent free-end shear/moment residuals; Euler-Bernoulli is a
diagnostic only. The scalar oracle follows the cosh-scaled free-free equation
in Khasawneh & Segalman 2019, Table 7, https://arxiv.org/abs/1811.03222.
G3 compares the first six HEX27/TET10 modes on the same bar within the sum of
their three-rung discretization estimates plus propagated eigen residuals.
These are estimated discretization bounds, not certified continuum enclosures.
All new numerical gates require central native execution before qualification.
No mixed/incompressible stabilization, reduced-integration stabilization,
finite-strain hexes, inside-element material variation or global injectivity
claim is introduced. Distorted bending-dominated parts require HEX27/TET10;
HEX20's affine benchmark does not qualify it on arbitrary distorted geometry.

Stiffness and mass use the existing `fs-feec` positive-weight collapsed Duffy
rule. Required counts are respectively 3..=12 and 4..=12 per direction: affine
stiffness has degree two and consistent mass degree four, while the Duffy map
adds degree two. Mass is never row-sum lumped. Curved-element stiffness is
rational and needs an independent quadrature ladder; specifying a count does
not certify its error. `tet10_shapes` and `tet10_point` expose the same mapping
for consistent loads and independent norm integration. Assembly is ordered by
element and quadrature point, stores symmetric local triangles identically,
and polls the caller's `Cx` before allocation and at each quadrature point.
Mesh/free-DOF caps refuse before allocating maps or matrix storage; no partial
matrix is published after an error or cancellation. Replay is same-host
deterministic; no new cross-ISA claim is made.

`promote_p1_to_p2` gives each sorted global endpoint pair one midpoint ID,
independent of element order. It canonicalizes element orientation and applies
an optional caller-supplied source-surface projection only to boundary edges.
Straight-sided approximation of curved geometry is a separate error axis.
Positive sampled Jacobians do not prove injectivity between samples; projection
does not supply source-surface or specimen authority. Finite-strain and mixed incompressible formulations remain outside this
quadratic operator's scope. Its explicit initial-stress extension is described below. Thermal stress observations are limited to
the explicit volume averages returned by `solve_thermal_response`.

Valve reduction requires `ValveCardRequest::element_order`: `Tet4` or explicit
`Tet10 { elements, quadrature }`. Quadratic corner sets must match the original
per-element topology, and every observation/pressure-load node must refer to
the supplied complete quadratic mesh. No boundary condition or load is
implicitly interpolated from P1. `ValveCard` v2 retains the discretization;
its identity binds order, quadrature and full quadratic connectivity. The v2
decoder refuses v1 bytes rather than reinterpreting them. This changes card
identity, not the legacy P1 numerical operator or the rl-001 authored band.
Reduction remains an explicitly truncated exploratory consumer, not an IRP
acceptance rendering engine.

### Thermal equilibrium in pressure-loaded valve reduction

`reduce::thermal::reduce_pressure_loaded_valve` carries the original P1 or
TET10 thermal/dead-load equilibrium into the existing valve/reed consumer.
The material field, integrated eigenstrain and explicit dead nodal forces
use the original solid owners. The common reducer receives that equilibrium's
`K_elastic + K_geometric` pencil; both its modal spectrum and exact opening/
pressure compliance therefore change with the same actual stress field.
TET10 retains pointwise quadrature stress, not element-mean stress.

The existing LDLT inertia owner must admit the WHOLE tangent as positive
definite. Selecting a positive modal window cannot hide an unstable direction
outside it. A singular/unstable pencil refuses instead of converting negative
eigenvalues into frequencies. This is numerical admission within the original
factorization tolerances, not a continuum buckling certificate.

Reference geometry gradients, mass and homogeneous incremental supports remain
unchanged. The aperture gap, width and thickness are measured from the accepted
displaced nodes against the original fixed observation axes and opposing plane.
The original P1 geometry-update owner checks its displacement bound and moved
cell quality. TET10 uses that declared displacement bound and the existing
isoparametric point evaluator at both original quadrature rules. Positive sampled
Jacobians do not establish injectivity between points. Closed apertures, excessive
strain/displacement and incompatible source states refuse; no contact model is
selected implicitly.

The pressure-force-per-pascal field remains an explicit FIXED reference load.
Its work-conjugate modal areas use the same new modes; surface-normal updates,
follower stiffness, pressure-dependent equilibrium and changed load geometry
are not inferred. The loss factor remains authored. This is the conservative
small-displacement initial-stress approximation, not finite-strain dynamics,
a source loss-spectrum model or live modal-state transfer.

`ThermalPressureLoadedValve` retains the original equilibrium/tangent,
displaced observation nodes, ordinary pressure-loaded valve and a combined
identity. The existing `BeatingReedSource` ingress compiles its actual pressure
port after the original first-mode compliance admission. The ValveCard codec
is unchanged; its ordinary material/mesh fields keep their original meaning.
The complete thermal receipt must accompany that card to retain thermal input
provenance. Every refusal or cancellation returns no partially minted reduction.

G1 definitions in `reduce/tests/thermal.rs` use a six-tet unit cube with one
free corner displacement. Independent barycentric integrals give closed-form
P1 and TET10 mass, stiffness, thermal/dead-load displacement, signed geometric
stiffness and pressure compliance; the same result reaches the playing-reed
ingress. G3 checks zero-preload parity and an unstable disconnected direction
outside the selected positive window. G4 covers displacement, closed-aperture,
source mismatch and cancellation refusals. These are synthetic implementation
checks, not measured valve/material validation.

The focused `tests/tet10.rs` target specifies G0 shape, quadrature, consistent
mass, rigid-mode, distorted bending patch, rotated-tensor, refusal and budget
checks; G1 MMS (L2 order 3 and energy order 2, ±0.2), projected-cylinder geometry
and free-bar eigenvalue ladders; and the Cook/MacNeal–Harder G2 bands. The
`reduce::reduce_tests::rl_005_tet10_cantilever_card` regression specifies a
fivefold static-compliance improvement over P1 at matched DOFs on the rl-001
slab geometry. `scripts/e2e/irp_s1_structures.sh` retains their structured
receipts. These new probes require central native execution before any DEF-08
closure claim. The Euler-disc structural consumer and Bar-S1 descriptions
require their own explicit migration; adding this operator alone does not
establish those consumers' convergence or remove their P1 limitations.

The free-bar per-commit test computes r=1,2,4 and explicitly compares against
the historical ts1 r=8 diagnostic identified in the .58 bead; it does not
qualify that historical rung. `g1_free_free_bar_eigen_ladder_r8_qualification`
computes all four rungs (the first diagnostic cost 2101 s on ts1). Its gates
require rising observed order toward the theoretical four, finest order at
least 3.5, and eigen residual uncertainty at most 2% of each finest difference.
Successive Richardson limits are compared with an explicitly estimated
pre-asymptotic truncation allowance plus propagated eigen residual uncertainty.
Mass scaled by `(1+h²)` must fail the same trend/limit gate. MacNeal–Harder
alone uses `8 eps kappa_1_est ||b||_1` for the ordinary-SpMV residual, with a
five-step inverse-norm estimate and a local 1e-6 displacement falsifier. This
is a conditioning estimate, not a verified inverse bound; all other static
fixtures retain their original relative Euclidean residual limit of 1e-9.

`TetElasticMaterial::try_new_oriented_tensor` admits complete anisotropic
stiffness data into that same operator. `ElasticTensorBasis` requires explicit
tensor-shear, engineering-shear or Mandel notation, one of the documented
`xx,yy,zz,xy,yz,zx` / `xx,yy,zz,yz,zx,xy` orders, and a nonzero source-frame
identity. The caller also supplies a target-frame identity, a proper rotation
whose columns are source axes in the target frame, density and an upstream
material-state identity. A frame mapped to itself requires identity rotation.
The constructor checks full major symmetry in work-conjugate coordinates and
positive definiteness of the entire source and target law. Tensor-shear
matrices with normal/shear coupling are generally not ordinarily symmetric;
engineering matrices transform as `C_M = S C_E S`, with
`S = diag(1,1,1,sqrt(2),sqrt(2),sqrt(2))`. Negative off-diagonal coefficients
are admissible when the whole law is stable. This does not apply a stiffness
positivity rule to arbitrary off-diagonal blocks of coupled constitutive laws.

`ElasticTensorTransformReceipt` retains the exact source matrix and descriptor,
both frame identities, rotation and upstream material identity alongside the
executable target-frame material. The material's identity binds those fields,
density and every output coefficient under
`org.frankensim.fs-solid.elastic-tensor-frame-transform.v1`. Rotation reuses the
existing Mandel basis transform (equivalent to the four-index tensor law) and
its deterministic triangle storage; source asymmetry is refused before that
step. No missing coefficient, noisy-source symmetry projection or frame
calibration is inferred. G0/G3 tests compare all three notations and both
orders against a separate four-index rotation with a fully coupled synthetic
tensor, check whole-law instability and frame/identity refusals, and retain
the existing orthotropic tests. A G1 affine single-tet test compares actual
nodal forces and strain energy against direct tensor contraction. This is
software evidence, not measured crystal/wood validation. The receipt has no
portable codec yet. `from_resolved_elastic_tensor` now connects the complete
`fs-material::state_point::ElasticTensorStatePoint` directly to this constructor.
The shared basis descriptors are owned by L1 `fs-matdb` and re-exported here.
Source compiler v6 and normalized-pack v5 preserve each coefficient context;
the state resolver requires all 36 entries plus density. The transform binds
the bundle identity, hence every selected receipt and the actual query point.
The G1 tet regression also decodes a portable material card, resolves the
complete tensor, and reproduces the independently checked assembled operator.
Coefficient uncertainties remain in the upstream receipts; no independent
coefficient assumption, propagated tensor bound or empirical authority is added.

The committed `silicon-cubic-25c-nasa-rp1057` source bundle now supplies a
separate sourced witness in `xtask/tests/matdb_pack_cli.rs`. It compiles all
36 explicit engineering coordinates, persists and reloads the material card,
and resolves its density/tensor into the same oriented tetrahedral operator.
The G1/G3 test checks mass, nodal forces and strain energy at zero and 45-degree
orientations against an independent three-constant cubic stress calculation.
The matrix entries derive from three handbook constants at 25 C; density is
from a separate pure-crystal reference. This cross-source engineering example
retains unknown pressure and specimen/dopant conditions. It neither qualifies
a particular wafer nor supplies a thermal/electrical silicon model.

`TetThermalStrainState` and `assemble_thermal_load` add the matching generic
small-strain thermomechanical load boundary. An upstream material law supplies
the already-integrated, possibly anisotropic Mandel thermal strain between an
explicit positive reference and current temperature. Each state is bound to
the exact elastic material-state identity and a nonzero thermal-law identity;
the isotropic binder consumes fs-material's exact signed
`integral(alpha(T) dT)` only when its current endpoint and tangent elasticity
share the exact card and complete state point. Otherwise it refuses before
assembly.
the operator assembles `integral B^T C epsilon_free dV` on the same tetrahedra
and publishes both full nodal forces and the constrained free-DOF projection.
`solve_thermal_displacement` advances a Jacobi-PCG solve in cancellation-bounded
chunks and accepts only an independently recomputed true Euclidean residual;
it reconstructs the full zero-Dirichlet nodal displacement field for an
explicit later geometry-update transaction. That transaction is
`update_geometry_from_displacement`: it identity-checks the reference mesh,
enforces a caller displacement ceiling, preserves connectivity, and refuses
every degenerate or orientation-flipped tetrahedron before publishing the
updated fixed-topology chart.
Fully fixed homogeneous-Dirichlet meshes now admit their exact zero-displacement
thermal response with zero iterations and zero reduced residual. The thermal
load still validates all geometry/material/state inputs and retains its full
force vector; recovered stress, support reactions and stored elastic energy
can be nonzero. Matrix/modal `assemble` continues to refuse `NoFreeDofs`.
A G1 fully clamped isotropic heating test checks hydrostatic compression,
reactions, energy and unchanged geometry against closed-form values; invalid
thermal counts and degenerate geometry still refuse before the empty solve.
No material name, constant-alpha assumption, or Euler-disc geometry appears in
this layer.

`TetThermalStrainState::try_from_oriented_strain` connects a source-frame
symmetric strain to that same thermal operator. `StrainTensorBasis` declares
physical tensor shear, engineering strain shear (`2 epsilon_ij`), or Mandel
shear and either existing six-component order. These are strain conventions,
not stress-vector or elastic-matrix scaling. The numerical transform is the
second-order law `epsilon' = Q epsilon Q^T`. The admitted
`ElasticTensorTransformReceipt` supplies the target frame and elastic-state
identity; a source frame mapped to itself requires identity rotation. Identity
rotation preserves normalized Mandel bits. Nonfinite input, improper rotation,
shear overflow, complete nonzero-to-zero rotation underflow, and target strains
outside the declared Mandel component ceiling refuse.

`ThermalStrainTransformReceipt` retains the source six-vector, basis, both
frames, rotation and executable thermal state. State identity domain
`org.frankensim.fs-solid.thermal-strain-frame-transform.v1` binds the ordinary
thermal-state identity (temperatures, law, elastic state, target strain and
validity ceiling) plus the source vector and complete transform. The existing
constructor's identity preimage is unchanged. G1/G3 compares all six source
conventions against independent second-order rotation, trace/Frobenius
invariance, and four-index elastic stress contraction; actual thermal loads,
constrained displacement and fixed-topology geometry reproduce affine
anisotropic expansion. This is nominal, synthetic numerical evidence. The
upstream law must supply the integrated strain and admit its temperature path;
no coefficient integration, frame calibration, portable transform codec or
uncertainty propagation is added here.

`try_from_resolved_strain` connects the complete six-component
`fs-material::StrainTensorStatePoint` (in its `state_point` module) to this same
operator ingress. The selected card/claim/query-point bundle identity becomes
the upstream thermal-strain identity; callers still declare the strain to be
stress-free and integrated over the supplied temperature path. The G1/G3 test
now loads every convention through a portable material card before checking
the same independent stress, nodal force, displacement and geometry oracles.

`try_from_joint_strain` additionally retains admitted source covariance or its
exact unknown-correlation reason. With a fixed proper rotation Q and fixed
admitted stiffness C, it computes `Sigma_e = A Sigma_source A^T` and
`Sigma_stress = C Sigma_e C^T`, where A includes source shear scaling, order
and the second-order frame transform. Both outputs use target Mandel order;
strain covariance is dimensionless squared and free-stress covariance is Pa
squared. No off-diagonal covariance is discarded. The fixed 6x6 calculations
publish symmetric results or refuse nonfinite values, negative variances and
complete loss of a nonzero covariance. No eigenvalue clipping manufactures
positive semidefiniteness.

`JointThermalStrainTransformReceipt` holds the executable nominal receipt,
conditional covariance outcome and source joint-state identity. Domain
`joint-thermal-strain-transform.v1` binds all of them and the output covariance
bits. G3 compares all six conventions and both source policies with a correlated
finite ensemble rotated by an independent full-tensor calculation, including
covariance of actual thermal nodal loads. G0 retains NoBlock, PartialMembership
and UnstatedMarginal outcomes and checks arithmetic loss/overflow. These are
conditional linear second moments, not confidence intervals or guarantees that
all realizations remain in the admitted small-strain/material domain. Uncertain
stiffness, orientation, temperatures, their cross-correlations and nonlinear
uncertainty propagation remain unsupported; unknown covariance is never zero.

`recover_thermal_stress` consumes the accepted thermal displacement and its
exact mesh, thermal field, material coefficients and constraints. It computes
the actual P1 elastic strain `B u - epsilon_free`, Cauchy stress
`C (B u - epsilon_free)`, stored elastic energy `V/2 epsilon_elastic:sigma`,
and full nodal internal force `integral B^T sigma dV`. These forces are support
reactions on fixed DOFs and equilibrium residuals on free DOFs; residual entries
are retained. Free expansion produces zero stress up to numerical error, while
constrained expansion produces the corresponding physical stress and energy.
The caller names the reference mesh frame and supplies a finite positive
ceiling applied to both total and elastic Mandel strains. This is a declared
small-strain applicability bound, not a yield or material-validation certificate.

Recovery rechecks the original thermal-load identity and the actual density
and stiffness bits retained under `thermal-response-input.v1`; changing a law
while reusing its caller-supplied material ID therefore refuses. Existing v1
load and displacement identity formats are unchanged. The immutable
`TetThermalStressReport` binds those inputs, solution, frame, strain ceiling,
element outputs, nodal forces and energy under `thermal-stress-recovery.v1`.
Cancellation is checked before allocation and at every element boundary; no
partial field is published. Existing mesh/element budgets bound report storage.

`stress_in_frame` expresses an element's recovered stress by `Q sigma Q^T`
with an explicit `StressTensorBasis`. Physical and engineering stress shear
are both `sigma_ij`; Mandel shear is `sqrt(2) sigma_ij`. Engineering-strain
factor two must never be applied to stress. Both existing component orders
are supported. Proper rotation and nonzero frames are required; a frame mapped
to itself requires identity rotation. Nonfinite and complete nonzero-to-zero
transforms refuse. `StressTensorTransformReceipt` retains source stress,
report/element identity, both frames, target convention/order, rotation and
output, binding them under `stress-frame-transform.v1`.

G1/G3 checks free expansion and two supported bodies with distinct thermal
strains/stiffness against closed-form displacement, stress, reaction and
energy formulas. All six output conventions are compared with independent
full-tensor rotation and stress/strain work conjugacy. G0/G4 checks reused
material IDs with changed stiffness, changed thermal fields/constraints,
strain bounds, invalid frames and pre-cancellation. This is nominal linear
elastic recovery in the reference geometry, not stress-source pack import,
finite-deformation/Piola stress, yield/buckling admission, evolved-state
recoupling, propagated solved-stress uncertainty or measured validation.

`JointStressTensorTransformReceipt::try_from_resolved` consumes the immutable
joint stress state from fs-material. It transforms source stress and full
covariance through stress-specific order/shear normalization and a fixed proper
rotation into canonical target Mandel coordinates. Tensor and engineering stress
both store physical shear; the transform never applies engineering strain's
factor two. The existing bounded 6x6 covariance kernel computes `A Sigma A^T`
and preserves output symmetry without clipping or assuming independence.
Nonfinite normalization/rotation/covariance, negative variance and complete
nonzero-to-zero numeric loss refuse. Nonzero target frames and identity rotation
for a frame mapped to itself are required; identity rotation retains normalized
nominal bits. The work is one bounded tile and allocates no growing buffers.

The receipt retains source joint identity, original coordinates and nominal
stress, target frame, rotation and both numerical outputs under
`joint-stress-frame-transform.v1`. `StressTensorUncertainty` holds the complete
target covariance [Pa^2] or the unchanged source unknown reason. G3 compares
all six source conventions and both joint policies against a full-rank
twelve-realization ensemble transformed as full tensors and passed through
actual clamped thermal displacement/stress recovery and checks receipt tampering.
G0 retains NoBlock, PartialMembership and UnstatedMarginal and checks frame
errors, shear overflow and complete subnormal loss. These are synthetic numerical tests.
This is source-data covariance conditional on fixed orientation, not solved-field
uncertainty, a confidence/domain bound, material/frame calibration, propagation
of uncertain temperature/stiffness/orientation, or experimental validation.

## Material replacement on an accepted P1 trajectory

`loads3::transient::events::rebind_material_candidate` permits source-data
changes to density, stiffness and stress-free strain on the same reference
mesh and fixed-DOF set. It retains the accepted physical displacement, clock
and step count, solves `M_after_ff v_after = M_before_ff v_before`, and uses
the original initializer to rebuild operators and re-equilibrate acceleration
under the explicitly supplied new load/thermal endpoint. Subsequent ordinary
Newmark steps and stress recovery consume the actual retained new problem.
The event does not consume a time-advancement slot.

The conserved quantity is FREE canonical momentum. Stationary constraints
may exchange full nodal or total physical momentum; this operation reports
no support impulse. Moving-support trajectories and changed meshes/restraints
refuse. This is externally supplied material reconfiguration, not mass
transport or a continuous material transition.

Required signed work is `Delta(1/2 v'Mv) + Delta(sum V/2 e'C e)`, where
`e = strain - strain_free`. Kinetic and physical shifted-elastic changes
remain separate. Positive work must be supplied by the reconfiguration owner,
and negative work returned to it; neither is automatically caloric heat or
dissipation. An explicit absolute work budget and independent SI momentum
residual gate admit the candidate. Per-element old/new material and thermal
identities and reference masses remain in its report. Failure of source,
solve, momentum, strain, work or cancellation checks preserves the old run.

`tests/transient_events3.rs` exercises immutable source-card substitution
through the native solver: doubled density, tripled stiffness, changed
eigenstrain, exact scalar momentum/acceleration/work, subsequent motion,
zero-work equivalent replacement, retained nonzero clock/step history,
fully clamped storage, and late work/strain or cancellation refusal.
This does not differentiate the replacement map or validate a real material
transition.

## Instantaneous impulses on an accepted P1 trajectory

`loads3::transient::TetTransientRun3::apply_nodal_impulse` applies
explicit nodal impulses through the same sparse consistent reference mass as
Newmark stepping. It solves `M_ff delta_v_f = J_f`, retains supported velocity,
and recovers full support impulses as `M delta_v - J` on restrained rows.
Momentum gates use the actual representable velocity change. Moving supports
retain their impulse work; `average(v) dot M delta_v` is checked against applied
plus support work without subtracting large endpoint kinetic energies.

The event preserves displacement, time, step count, material, thermal state
and finite loads. Native acceleration and finite-force reactions are refreshed
for the new velocity and damping force before transactional publication.
It does not advance a strain/history law or generate damping heat at zero time.
Collision timing, restitution and the opposite partner remain with the caller.
`tests/transient_impulse3.rs` supplies direct consistent-mass, support-work,
continued-motion and refused-candidate regressions. These are local numerical
checks, not calibration of an impact law.

## Thermal material evolution in sparse P1 dynamics

`loads3::transient::thermal::advance_thermal_candidate` refreshes the actual
source-bound thermal strain, material identity and global Mandel stiffness.
Changed tensors rebuild the native endpoint stiffness using the existing P1
assembly; unchanged tensors reuse the original sparse operators. Both paths
preserve the accepted Newmark displacement, velocity, acceleration and clock.
Reference density, sparse mass, reference inventory, mesh connectivity and
stationary zero restraints stay fixed. Refreshed source counts, identities,
geometry, quality and strain limits pass the original assembly and response
checks. A fully clamped trajectory uses the existing zero-free-DOF dynamic
assembly and still reports stress and stored energy.

The ordinary Newmark audit retains its unshifted quadratic-energy convention,
averages old/new thermal nodal forcing, and records signed stiffness parameter
work `W_K = q_old' * (K_next-K_old) * q_next / 2`.
`ThermalTransientStep3` separately audits physical kinetic plus shifted elastic
energy, with `e = epsilon-epsilon_free`. Its independent element work is
`W_T = -sum V * (sigma_old+sigma_next)/2 : delta(epsilon_free)` and
`W_C = sum V/2 * e_old : (C_next-C_old) : e_next`. Both terms are needed when a
fully restrained solid expands and softens with zero displacement work.
Material parameter work is signed prescribed input; it is not dissipated heat.

Changing stiffness requires zero stiffness-proportional Rayleigh damping;
fixed reference mass damping remains supported. Changing endpoint Rayleigh
stiffness damping is not unconditionally dissipative and is refused.
Unchanged-stiffness trajectories retain their existing damping behavior.
The native and physical audits use the declared run tolerances and refreshed
source identities. These are exact discrete endpoint work identities;
continuum accuracy still requires resolving the prescribed material path and
checking timestep refinement. Failed assembly, motion, stress, energy or
cancellation checks leave the original trajectory unchanged.

G1/G3/G4 tests in `tests/transient3.rs` cover an analytical thermally driven
oscillator with timestep refinement, source-identity refresh, scalar softening
and stiffening with retained motion and mass damping, fully clamped thermal
storage with distinct expansion/stiffness work, changed-density and variable
stiffness-damping refusal, late strain refusal and cancellation. This is
one-way prescribed thermal/material loading. Reciprocal caloric feedback,
changing reference mass, moving thermal supports, an adjoint of the new
material inputs and measured validation remain outside this path.

## Prescribed support motion in sparse P1 dynamics

`loads3::transient::motion` extends the existing sparse average-acceleration
Newmark trajectory with absolute world-component support displacement. The
fixed DOF set, reference mesh, materials, thermal eigenstrain and Rayleigh
damping remain frozen. `begin_transient_with_supports` accepts full initial
displacement/velocity fields and explicit constrained acceleration; free
acceleration entries in that input must be zero because they are solved.
Initialization and every step include the complete partitioned load
`f_f + fT_f - M_fc a_c - C_fc v_c - K_fc u_c`, with
`C = eta M + zeta K`. Consistent-mass coupling is retained.

`advance_with_supports` accepts sorted, unique fixed-DOF displacement values.
Omitted fixed components mean zero at that endpoint. Support acceleration
and velocity follow the same beta=1/4, gamma=1/2 recurrence as the free state;
they are not finite differences with an unrelated integration convention.
Recovering acceleration from rounded displacement samples amplifies their
roundoff by the inverse squared step size; a smaller step is not an unlimited
accuracy guarantee for prescribed derivatives.
Full element actions reuse the canonical P1 stiffness/mass kernel. No dense
global operator is assembled. Fully prescribed meshes use empty reduced
matrices and zero solve iterations after full geometry/material validation;
public modal assembly still refuses `NoFreeDofs`.

`nodal_motion` includes the accepted support u/v/a. State energy includes all
nodes, and `TransientStep3::support_work_j` is the trapezoidal work of support
reactions **on the solid**. Mechanical work includes external loads on both
free and prescribed components; the frozen full thermal equivalent load has
its own work term. The acceptance identity is
`delta(kinetic + quadratic_strain_energy) + damping_loss = mechanical_work +
thermal_equivalent_work + support_work`. Quadratic strain energy remains
`u'Ku/2`; `recover_response` separately reports the actual thermally shifted
elastic energy and stress, full inertia/damping and support reactions.

Input admission, full nodal displacement limits, element limits, sparse solve
limits, per-node/element cancellation checkpoints and the energy/residual
gates precede atomic publication. A clone retains the prescribed kinematics,
reaction endpoint and immutable full thermal load as well as the free state.
Force-only advancement and reduced-only state adjoints
explicitly refuse this lane; they cannot omit
support terms while claiming a valid trajectory or gradient. Displacement
impulses, changing restraint axes/topology, adaptive event discovery,
temporal/continuum error certification and measured validation are not supplied.
Prescribed motion requires a physically resolved time-step ladder.

`loads3::history::dynamics` connects typed `SupportMotion3` displacement
histories to this same inertial trajectory. `begin_transient_with_support_history`
requires the initial restrained displacement to match the retained history
sample, with explicit initial velocity and constrained acceleration.
`advance_history_with_supports` samples all loads and support values together
at the run's next physical time before delegating one coupled Newmark step.
Value-table clamping or wrapping does not replace physical time for moving
frames, and rotating a payload does not rotate the fixed world restraint axes.
Malformed bindings, sampling failures and cancellation leave the accepted
trajectory intact. These adapters do not infer analytic derivatives, impulses
or an event schedule from displacement tables.

`EventLoadHistory3::new_with_supports` explicitly indexes those original
support, load and frame knots together under one schedule budget. Its support
initializer requires the same full initial kinematics. Event stepping uses
the exact left force before a knot and the right force after it, retaining
the native support Newmark recurrence. Finite force jumps keep support u/v/a,
include the original inertial/damping/stiffness lift in free acceleration,
and refresh reactions for subsequent support-work accounting; they add no
impulse or energy. A support-displacement jump, including one caused by a
frame jump, refuses without advancing accepted motion. The existing adaptive
two-half-step comparison now includes constrained-node kinematics as well as
free motion, including fully prescribed bodies. Support-driven runs compare
raw full/half disagreement for all kinematics: prescribed derivative errors
also enter free acceleration through the consistent mass. No Richardson order
is assumed for that recurrence, which can retain alternating acceleration
error. This local consistency indicator is not a derivative-error bound;
neither unknown events, exact support derivatives nor a global temporal bound
are inferred. Accepted work, clock and fork/restart semantics
remain those of the original trajectory.

`tests/transient_supports3.rs` supplies G1 exact mass/stiffness/both-Rayleigh
cross-term initialization, a closed-form accelerating-base oscillator with
second-order refinement, fully prescribed thermal/external/support work and
reaction checks, and G3/G4 fork, malformed-input and cancellation checks.
`tests/transient_support_history3.rs` covers rigid translation with damping,
agreement with explicit core endpoints, physical frame time beyond a clamped
table, a base-excited oscillator with both held-force pulse edges, full-support
adaptive rejection, and transactional sampling, binding and cancellation refusals.

### Prescribed-support motion adjoints

`loads3::transient::adjoint::support` supplies a separate full-state transpose
companion through `prescribed_adjoint`. `NodalMotionCotangent3` includes q/v/a
on every node, so direct objective dependence on constrained motion is retained.
`pullback_step` returns previous full-state cotangents, endpoint free-force
cotangents, and absolute support-displacement control cotangents in sorted
`fixed_dofs()` order. An omitted primal fixed row is still a zero displacement
control; varying it requires explicitly adding that row at the endpoint.
`pullback_initialization` reverses the actual initial coupled mass solve and
returns initial full q/v and constrained-a derivatives. Free initial a was
solved and its unused input cotangent is exactly zero.

The reduced solve reuses the original sparse transpose owner. A bounded
canonical element sweep supplies actual `M_fc^T lambda` and `K_fc^T lambda`
without assuming bitwise matrix symmetry. For `b=h²/4`, `d=h/2`, the support
endpoint seeds are `q_bar-K_fc^T lambda`,
`v_bar-(eta M_fc^T+zeta K_fc^T)lambda`, and `a_bar-M_fc^T lambda`.
The pullback then reverses the original prescribed Newmark recurrence, including
the previous support acceleration. Fully prescribed meshes require no sparse
solve but retain the kinematic transpose and direct objective seeds.

Full-DOF and element caps precede scratch allocation; independent sparse
transpose/solve limits remain explicit. Failed or cancelled pullbacks publish
no partial gradient and leave both the companion and primal run reusable. These
are derivatives of fixed-step discrete motion equations, not of PCG iterations,
acceptance, floating-point rounding, moving constraint axes, impulses, or step
times. Direct stress/reaction/work/energy-observation and material derivatives
are not supplied by this motion pullback. Prescribed derivatives retain the
forward recurrence's inverse-squared-step conditioning.

`tests/transient_support_adjoint3.rs` checks the analytical one-free-DOF
mass/stiffness/both-Rayleigh transpose, every initial/force/support control
group through three native steps against centered differences with full nodal
objectives, the fully prescribed kinematic limit, fixed-index ordering and
malformed-input, budget and cancellation refusals.

`material_scale_pullback` now accepts prescribed-support runs at their exact
accepted endpoints, using actual full nodal q/v/a. For the corresponding
step/initialization free-force cotangent lambda, it contracts
`-lambda' d(M*a + (eta*M + zeta*K)*v + K*q - fT)` with zero constrained
equation multipliers. Canonical P1 element matrices and Mandel B are reused;
the complete `M_fc*a_c`, `C_fc*v_c` and `K_fc*q_c` terms are retained, including
the explicitly supplied initial support acceleration. Returned independent
per-element stiffness/density scale derivatives and global Rayleigh eta/zeta
derivatives compose with the existing full-motion pullback. A common material
scale sums its element entries; fully prescribed motion has zero implicit
material sensitivity.

The free thermal strain is frozen, but the stiffness derivative includes its
equivalent-force dependence. Do not count that stiffness-dependent thermal
load a second time. Geometry, support axes/authored kinematics and external
nodal forces are fixed. Parameter-dependent external forces, initial/support
controls and direct objective terms must be chained separately. This is not a
tensor-entry, thermal-law, shape, reaction or stress-observation derivative.
Independent full-DOF and element bounds precede scratch allocation; cancellation
or refusal returns no partial result and never changes accepted motion. The
reduced-state-only internal helper continues to refuse driven supports because
reduced q/v/a cannot recover their physical cross terms.

`tests/transient_support_material3.rs` checks the analytical initial coupled
mass/stiffness/damping terms, a three-step anisotropic heterogeneous thermal
trajectory against six independent material/damping finite-difference
directions at two parameter points, and the fully prescribed zero-gradient,
budget, cancellation and retry limits.

`adjoint::support::reaction` supplies physical support-force observation
pullbacks through `reaction_pullback`. The observed force is the force ON the
solid, `R_c=(M*a+(eta*M+zeta*K)*v+K*q-fT-f_external)_c`, in fixed world
Cartesian components. A caller supplies finite objective derivatives with
respect to sorted unique fixed DOFs, typically from a recovered reaction
residual. The result includes full motion seeds `(K' w,C' w,M' w)`, direct
fixed-force rows `-w`, and direct material/Rayleigh terms `+w' dR` at fixed
motion. The full canonical transpose and material contraction owners are
reused; no global full matrix, new equilibrium solve or finite difference is
introduced. `reaction_motion_pullback` on the prepared prescribed companion
supplies the same motion seed without a primal-state tape.

For a complete reaction objective, chain the motion seed through the existing
initialization/step/trajectory reverse, then add the direct material and
fixed-force terms. Fully prescribed motion has no implicit material derivative,
but its physical reactions retain nonzero direct thermal, inertia and damping
derivatives. The same `-w` applies independently to a full frozen thermal
equivalent force. Stiffness-dependent thermal effects already present in the
material term must not be counted again. The endpoint API supports stationary
and driven restraints, keeps full consistent-mass and damping coupling, and
checks independent full-DOF/element allowances before scratch allocation.
Malformed rows, nonfinite arithmetic, refusal and cancellation publish no
partial result and never change accepted motion. No moving-axis, shape,
thermal-law, stress-observation, physical-work or acceptance-gate derivative
is inferred.

`tests/transient_reaction_adjoint3.rs` checks analytic direct and implicit
initial reaction terms, the same observation seed on stationary and driven
supports, three-step reaction-target finite differences in seven independent
material/damping/support/free-force/fixed-force directions, fully prescribed
nonzero direct material sensitivities, and malformed/budget/cancellation retry
boundaries.

`adjoint::support::stress` adds pointwise elastic stress objectives through
`stress_pullback`. Rows contain sorted unique tetrahedron IDs and finite
`dJ/dsigma` in mesh-frame Mandel `[xx,yy,zz,xy,yz,zx]`; the final three stress
entries are `sqrt(2)` times physical tensor shears. The canonical P1 strain
operator and admitted constitutive tensor supply `B' C' w`, including fixed
DOFs, and the direct per-element stiffness-scale term `w' sigma`, where
`sigma=C*(B*q-epsilon_free)`. Velocity, acceleration, direct density and direct
Rayleigh terms are zero. No element-volume or time-quadrature factor is added;
the objective owner includes any such weights. For complete design gradients,
chain the full displacement seed through the matching accepted-state reverse
and add its implicit material terms to the direct stiffness contribution.

Geometry, material axes and free thermal strain remain fixed. Stationary,
driven and fully prescribed endpoints use their actual full displacements.
Full-DOF and element caps apply before scratch allocation; both total and
elastic Mandel strains must obey the same finite positive ceiling as stress
recovery, even for unobserved elements. Invalid rows, nonfinite arithmetic and
cancellation yield no partial gradient and preserve accepted motion.
`von_mises_limit` supplies the smooth dimensionless inequality
`(sigma_vm/limit)^2-1 <= 0` and its Mandel stress seed, with a defined zero
derivative at hydrostatic stress and a caller-supplied positive finite limit.
It does not infer a yield law, failure probability or physical admission.
`tests/transient_stress_adjoint3.rs` checks analytic thermal/shear terms without
volume weighting, anisotropic nodal and material finite differences at two
mesh scales, von-Mises derivatives, and refusal/cancellation/retry boundaries.

`adjoint::support::trajectory` records complete prescribed-support simulations
through `record_prescribed_controls`. For controls alone, its owning recorder
retains only accepted interval times: the frozen linear motion Jacobian needs
no full primal tape or finite-difference replay. Each successful original coupled step is recorded
atomically; failed sampling, support admission, solves, cancellation or storage
admission preserves both the accepted state and schedule. Cloning forks both.
The typed-history adapter samples one retained force/support endpoint before
advancing, using the same physical time and frame binding as the original run.

`reverse_with_observations` adds sparse full-nodal objective derivatives at any
accepted endpoint, including initialization, and reverses the complete schedule
once. It returns every initial/endpoint force row, every absolute prescribed
displacement row, independent initial free q/all v/fixed a, and the sum of force
cotangents for the one frozen thermal equivalent load. Initial fixed q appears
only in the initial displacement-control row, avoiding duplicate control
contributions. Sparse observations preserve caller summation order with
O(steps + observations) indexing. Explicit step, returned-gradient and observation
caps precede allocation. Refused/cancelled reversals return no partial gradient.

`enable_material_recording` additionally retains the original full nodal q/v/a
at initialization and each accepted endpoint. It must be enabled once, before
the first finite step. Independent element/full-DOF limits and a checked
`samples * nodes * 9` scalar cap precede recording; failed advancement preserves
both the accepted motion and material tape. Cloning forks the retained tape.
The typed-history adapter uses the same bounded atomic recording path.

`reverse_design_with_observations` returns all existing control derivatives and
per-element stiffness/density plus global Rayleigh damping derivatives in one
reverse sweep. It contracts each residual multiplier against the exact retained
full motion, including the solved initial acceleration. No additional transpose
solve, replay or finite difference is needed. Material derivatives hold geometry,
restraint axes, authored support kinematics, free thermal strain and independent
external nodal loads fixed. The stiffness-induced thermal equivalent force is
already included; its dependence must not be chained again through the separate
frozen-material thermal-force cotangent. Parameter-dependent external loads and
direct objective terms remain caller-owned chains. Refusal/cancellation returns
no partial combined gradient and leaves the tape reusable.

`reverse_design_with_stress` adds sparse element-stress observations at any
recorded endpoint to those motion objectives. `ElementStressObservation3`
contains objective/Pa seeds in mesh-frame Mandel `[xx,yy,zz,xy,yz,zx]`; duplicate
element/state rows add in caller order, and callers include quadrature weights.
The same reverse sweep accumulates full `B' C' w` displacement seeds, direct
`w' sigma` stiffness terms, and implicit trajectory material/load/support terms,
including initialization. No additional equilibrium solves or stress tape are
needed. The combined motion/stress row cap and existing element/full-motion
budgets apply. Every state with stress rows obeys the finite positive total and
elastic strain ceiling, even when its seeds cancel to zero. Thermal strain,
geometry, frames, restraint axes and accepted times remain fixed. The smooth
`von_mises_limit` seed composes with this path for strength objectives; it does
not add plasticity or infer material failure limits. The three regressions in
`tests/transient_stress_trajectory3.rs` cover six material/damping/support/load
directions for mixed stress/strength/motion objectives, fully prescribed
intermediate-state controls, and observation-budget/strain/cancellation retries.

Callers supply actual residual/quadrature weights for their motion objectives
and chain returned world-component controls through authored waveforms. This
does not infer derivatives of frame transforms, knot/step times, interpolation,
geometry, direct control costs, physical work reports or acceptance.
`tests/transient_support_trajectory3.rs` supplies a multi-step least-squares
objective checked against all twelve independent initial/force/support control
directions, an analytical typed fully-prescribed ramp-amplitude chain rule,
and bounded recording, cloning, cancellation and refusal/retry checks.
`tests/transient_support_design3.rs` adds a thermal multi-step least-squares
objective observed at every endpoint, compared with six material/damping/support/
force directions at two parameter points, plus opt-in recording, storage refusal,
atomic failure, cancellation/retry, clone identity and fully prescribed limits.

## Public types and semantics

- `Mesh2` / `Patch`: structured body-fitted meshes — P1 triangles, Q1
  quads, mapped Q1 panels (`cooks_membrane`), with the four sides as
  named patches. `shapes_at`/`quad_points` expose the isoparametric
  machinery (chain rule with the transposed inverse Jacobian — the
  conformance battery pinned this exact slip on triangles and mapped
  panels).
- `LinearProblem` (plane strain/stress via `PlaneKind`, `lame`):
  small-strain elasticity with strong Dirichlet and dead-load traction
  patches; `Formulation::Standard` vs `Formulation::BBar` (Hughes'
  dilatation projection — element-average volumetric gradient,
  equivalent to condensed constant-pressure mixed elements).
- `HyperProblem` + `NewtonSettings`/`NewtonReport`: plane-strain
  finite strain; the 2D displacement gradient embeds into the 3D
  deformation gradient (F₃₃ = 1); residual `∫ P : ∇δu` and consistent
  tangent `∫ ∇δu : A : ∇Δu` come from fs-material's exact AD Piola
  stress and 9×9 tangent. Newton globalizes by Armijo backtracking on
  the potential energy (material refusals — det F ≤ 0 — halve the
  step) under linear load stepping; systems solve with MINRES
  (symmetric, indefinite-safe near buckling). `residual_and_tangent` /
  `potential_energy` are the public consistency probes.
- `CutElasticity`/`CutSolution`: source-compatible legacy facade over
  `fs_cutfem::CutElasticity` and its canonical solution type. The facade
  constructs `IsotropicElastic::new(E, ν, 1.0)`; the fixed strain limit
  supplies compatibility metadata only and is not evidence that a solved
  field remains inside that domain. Historical stabilization knobs retain
  their physical coefficients by translating
  `β_canonical = β_legacy(λ+2μ)/μ` and
  `γ_canonical = γ_legacy(λ+2μ)/μ` before delegation. Canonical assembly,
  graded-tree constraints, Nitsche and ghost terms, deterministic CG,
  retained topology, and error integration are not reimplemented here.
  The facade fixes the legacy solver controls at tolerance `1e-12` and
  `60_000` iterations. The removed implementation's post-hoc
  `1e-8` acceptance after targeting `1e-12` is not reproduced: a miss
  of the canonical `1e-12` gate fails closed. `BoundaryTraction`,
  `DesignBoxEdge`, and opaque checked `EdgeBand` are direct fs-cutfem
  re-exports. `solve_with_boundary_traction` delegates the typed support
  descriptor without inspecting or reimplementing it; the legacy callback
  field must be `None`, and an SDF crossing through supported load still
  refuses canonically.
  Physical-coefficient equivalence is claimed only when each active term's
  dimensionless translated coefficient is finite. Canonical validation refuses
  an overflowing Nitsche translation when `traction_free_interface` is false
  and any overflowing ghost translation. When `traction_free_interface` is
  true, the unused translated Nitsche value is deliberately ignored and
  carries no physical-coefficient claim. A future physical-coefficient input
  surface is required to represent active terms without a non-representable
  intermediate.
- `select_formulation(RegimeIndicators)`: the documented
  element-selection policy — B-bar when ν ≥ 0.45 or slenderness ≥ 5
  (fs-regime's indicators feed this; the locking battery measures the
  boundary).
- `stability` (bead tfz.15): `reduced_pencil` (free-DOF K and
  geometric stiffness K_G(σ₀) from the linear prebuckling state),
  `buckling_loads` — the pencil `(K + λK_G)φ = 0` solved MATRIX-FREE
  by fs-la LOBPCG on the Cholesky-reduced operator L⁻¹(−K_G)L⁻ᵀ
  (dense reduction fixture-gated at n ≤ 4096), modes −K_G-normalized;
  `lambda_indicator` (Richardson) + `evidence_row`;
  `eigenvalue_derivative` (direct pencil derivative at frozen
  prebuckling stress), `group_stiffness` (per-group Young's-scale
  lever), and the clustered-eigenvalue policy: `ks_aggregate` /
  `ks_aggregate_derivative` — the conservative smooth KS lower
  envelope with softmax-weighted derivatives, because individual
  eigenvalues are not differentiable where branches cross (measured
  in the battery).
- `continuation` (bead tfz.15): `PathResidual` (generic equilibrium
  residual R(u,λ) = R_int − λF_ext; `HyperProblem` implements it for
  homogeneous clamps + ramped dead tractions), `advance` — Keller
  bordered pseudo-arclength with pivoted dense LU for total bordered dimension
  at most 1024 and matrix-free FGMRES above that cap. Both predictor and
  corrector solve the complete, row-equilibrated nonsymmetric border, so a
  simple fold does not require an inverse of the singular structural tangent.
  Corrector acceptance requires both physical force balance and the tangent-plane
  arclength constraint; a lost displacement increment cannot pass by returning
  to the previous equilibrium. Step reduction/growth respects the declared
  minimum/maximum, and failure preserves the complete pre-step checkpoint,
  limit-point events (λ̇ sign flips) and branch-point candidates
  (Cholesky definiteness flips without a λ̇ flip, gated n ≤ 1024);
  `PathState` — the checkpointable trajectory (clone to checkpoint,
  hand back to resume; split runs are BITWISE identical to straight
  runs); `switch_branch` — null-direction predictor pinning
  (`pending_switch`): the arc constraint forces the first step along
  the buckling mode so the corrector lands on the bifurcated branch
  (a state-only perturbation relaxes back to the fundamental branch;
  a basin-scale jump inverts elements — both measured failure modes).
- `rod` (bead tfz.14): geometrically exact Cosserat rods — nodal
  positions plus canonical `fs-ga::So3` frames updated
  MULTIPLICATIVELY through body/right group operations; strains Γ = Rᵀr′ − e₁ and
  κ = Log(Rᵢ⁻¹Rᵢ₊₁)/L₀ from RELATIVE quantities (rigid motions give
  exactly zero strain — energy-invariance battery). Statics: analytical native
  residual/tangent Newton with residual backtracking and load
  stepping, fixture-gated; fixture slenderness matters — a
  penalty-stiff EA puts a quadratic spurious-stretch wall around
  every Newton iterate (measured, documented in the battery).
  Arbitrary energy callbacks retain an element-local finite-difference fallback;
  SE(3) dynamics under fs-time integrators remain successor scope.
- `fiber` (bead tfz.14): fiber sections over fs-material `Uniaxial`
  cards (Menegotto–Pinto steel; Mander concrete with the
  compression-positive convention flipped at the section boundary);
  `respond`/`commit` split (pure response vs state commitment);
  `update_sections_batched` — the §15.2 hot loop through fs-la's
  batched Cholesky with a per-outlier direct fallback, throughput
  measured and ledgered; `rc_section` fixture (confined core +
  unconfined cover + steel layers).
- `fiber::material`: source-selected or explicitly authored Mander/MP laws
  share the existing state-point resolver and fiber engine. Every density and
  constitutive parameter is required; selected source bundles and receipts
  survive section construction and area scaling. Static source admission
  refuses unowned strain, stress, rate, frequency or history dependence.
  Selecting coefficients does not calibrate the cyclic model. `RcSection`
  partitions four-sided cover, core and two symmetric reinforcement bands;
  steel displaces concrete in both response and mass. Steel response is lumped
  at the declared band centroids, so this is smeared reinforcement, not resolved
  round-bar geometry. Geometry and source bindings stay immutable; consumer
  clones own independent constitutive histories. Area scaling changes mass and
  every fiber area together. The legacy `rc_section` remains an authored fixture
  with its original gross-concrete-plus-rebar convention. Focused G0/G1/G3 tests
  are in `tests/fiber_material.rs`; no measured-building validation is implied.
- `beamcol` (bead tfz.14): force-based (Spacone) cantilever element —
  EXACT equilibrium interpolation, five-point Gauss–Lobatto section
  integration (a point AT the hinge), per-section Newton with
  residual backtracking, virtual-work tip deflection; pushover and
  cyclic histories with committed, cloneable state (G4).
- `frame2d` (plan §15.2 step 3): planar frames of force-based
  distributed-plasticity beam-columns (`Frame2d`, `FrameElement`,
  `SectionModel::{Elastic, Fiber}`). Basic system `q = [N, Mᵢ, Mⱼ]` with
  EXACT equilibrium interpolation `D(ξ) = [N, (ξ − 1)Mᵢ + ξMⱼ]`, 3–7
  Gauss–Lobatto sections, Neuenhofer–Filippou state determination
  (compatibility correction + section Newton to compatibility AND section
  equilibrium; sub-incremented `v_c → v` in up to 64 substeps when the
  direct iteration fails), tangent `K_b = F⁻¹`. `Geometry::Linear` or
  `Geometry::Corotational` (exact chord rotation, consistent geometric
  stiffness, chord angle unwrapped against the committed angle). Analyses:
  `static_load` (load control), `displacement_control` (bordered pushover
  through softening), `natural_frequencies` (massless DOFs statically
  condensed), `newmark` (average acceleration, uniform base motion,
  Rayleigh damping on the initial stiffness via `Rayleigh::from_modes`) with
  a `DynamicHistory` energy ledger (input, kinetic, damping, internal work)
  and base shear. `newmark_with_load` retains an explicit constant nodal load
  throughout the base-motion record: `M a + C v + f_int = load - M r a_ground`.
  It includes that load in the initial acceleration and external-work ledger,
  so gravity-preloaded columns keep their compression and P–Δ response.
  `newmark` continues to mean zero sustained load; a preload is never inferred
  from internal forces. The energy audit uses the change in kinetic energy
  from the supplied initial state, alongside incremental internal work.
  Newmark solves directly for endpoint acceleration with tangent
  `M + gamma*h*C + beta*h^2*K_t`; position and velocity use the corresponding
  increments. It never recovers acceleration from two absolute displacements.
  Backtracking and the residual gate act on that same physical acceleration
  equation. Force work uses the resolved kinematic displacement increment,
  retaining momentum/work balance when the stored position cannot resolve it.
  `velocities()` and `accelerations()` expose the accepted base-relative
  generalized state, including angular components. Absolute displacement
  remains `f64`; compensated position or sub-ULP constitutive strain is not
  claimed. `tests/newmark_acceleration.rs` covers lumped and consistent mass,
  sub-ULP vibration, changed base forcing, damping and kinetic/internal work.
  Activating a fixed coordinate freezes its accepted displacement and clears
  its velocity/acceleration. This is a kinematic boundary reset; no resolved
  support-engagement impact or associated impulse/energy receipt is claimed.
  The initial Rayleigh reference tangent is frozen at the first admitted
  transient initialization (after any gravity/static preload) and retained
  across record boundaries, clones and retries after a partially accepted
  history. Static loading does not silently replace it. Adding a node or
  member invalidates the reference; `reset_rayleigh_reference` explicitly
  requests a new one without resetting mechanical or material state. Current
  point masses still enter the mass-proportional term, and supports restrict
  the same full reference operator. Explicit coefficients may change between
  records; that is a caller-requested change of the damping model.
  Every frame, including one without material viscosity, refuses nonfinite
  forcing/directions and negative/nonfinite Rayleigh coefficients before state
  mutation. For nonzero stiffness-proportional damping, the complete free-DOF
  C (including mass and material loss) must have nonnegative work within a
  floating-point sign-resolution check. Positive diagonal congruence scaling
  followed by the existing Cholesky owner with a scratch-only
  `64*eps*n_active` diagonal allowance admits semidefinite loss directions;
  negative diagonals and nonzero rows with zero diagonal are refused.
  The original C is never shifted, clipped or replaced in forces or work.
  This numerical admission is not an interval passivity certificate and does
  not suppress physical buckling in the stiffness/internal-force model.
  The focused frame tests cover nonlinear whole/split-history agreement,
  explicit reference changes, structural edits, negative-damping refusal,
  mass-stabilized damping and retained zero-loss directions.
  The `static_load` and `newmark` Newton steps backtrack
  (step 1, ½, … down to 2⁻⁹) on the free residual norm, because fiber laws
  are non-smooth: Mander concrete has zero tangent at exactly zero strain, so
  the first tangent of an unstressed RC section is steel-only and full steps
  from rest overshoot and cycle.
- `koiter` [F], feature `koiter-asymptotics` (off by default): FD
  energy expansion along the buckling mode at the critical state →
  a/b coefficients and `Bifurcation` classification, with the
  SAMPLED-CONTINUATION fallback oracle (the imperfect-geometry path)
  cross-checking imperfection tolerance.
- `contact` [F], feature `contact` (off by default): SDF-native
  log-barrier contact around `HyperProblem`, with `Barrier`,
  optional lagged smoothed-Coulomb `Friction`, `ContactProblem`,
  `ContactSolution`, and the frictionless `translation_gradient`
  adjoint probe. The solver requires a finite, intersection-free
  initial state; bad parameters, non-finite SDF values/gradients,
  invalid pins, and unsupported contact-adjoint inputs return
  `SolidError::InvalidInput`. `solve_from` holds the previous physical
  displacement fixed as the slip reference; `solve` uses the zero state.
  `solve_from_with_surface_displacement` additionally admits one finite
  world-coordinate material-displacement vector per node. Friction uses the
  relative increment `u - reference_u - surface_displacement`; the prescribed
  surface increment stays fixed and is not scaled by the load factor. The
  supplied SDF independently defines current obstacle geometry, so a conveyor
  can transport the solid without changing its planar chart. Reference-state
  feasibility is still checked against that chart; no surface trajectory or
  material correspondence is inferred.
  Numerical lag rounds refresh normal loads and tangent frames only,
  and success requires the refreshed friction residual to meet tolerance.
  The frozen-frame Huber energy, gradient, and Hessian agree in both
  sticking and sliding branches. Exact Newton directions are retained when
  finite, descending, and admitted by a true max-norm relative linear
  residual at most `sqrt(f64::EPSILON)`. Otherwise sliding search directions
  first add the Huber quadratic majorizer curvature
  `mu * normal_load / abs(slip)`. This positive search metric changes neither
  the physical Hessian nor the force law. A stiffness-scaled diagonal shift
  remains the fallback if the majorized system also fails direction admission;
  it preserves the physical force law, tangent, pin equations, and final
  residual gate. Both line-search acceptance branches
  require admitted material stresses as well as positive sampled gaps;
  an energy decrease alone cannot bypass the material's stretch domain.
  At satisfied pin equations, the residual-decrease shortcut also requires
  energy nonincrease within a scaled floating-point allowance, preventing
  alternating energy/residual merit cycles. A tiny admitted Newton correction
  may terminate despite unresolved potential differences only if its actual
  force residual already meets the unchanged tolerance; its coordinate size
  is bounded by `sqrt(EPSILON)` times the shortest reference edge in max norm.
  This numerical-resolution rule is not an energy-error certificate. Pending
  prescribed displacement equations retain their boundary-work path.
  `eps_v` is a displacement smoothing
  width; this regularized law permits small creep and does not claim exact
  set-valued static friction. Solutions expose separate normal and
  tangential reactions plus nonnegative endpoint relative friction work.
  Signed `friction_support_work = sum(f_t dot surface_displacement)` records
  supplied support work, with `f_t` the force on the solid. Endpoint products
  close `body_friction_work + friction_work = friction_support_work`.
  These products are not trajectory-integrated heat or actuator work;
  thermal deposition and rigorous continuous-trajectory CCD remain separate
  capabilities.
  `solve_timed_from` additionally admits a positive physical interval and a
  declared slip-speed regularization. Its displacement smoothing width is
  `dt * regularization_speed`, so subthreshold creep follows elapsed time
  rather than the number of equilibrium solves. Existing entry points retain
  their displacement-width semantics. `ContactTrajectory` binds one unchanged
  problem, holds accepted displacement and physical time, and accumulates body,
  relative and support work with backward-Euler endpoint quadrature. A failed
  solve or nonrepresentable time/work update changes none of that state.
  The cumulative ports retain `body + relative = support` up to roundoff;
  body work is independently evaluated from actual displacement increments.
  This is quasistatic velocity-regularized friction: it adds no inertia,
  recoverable tangential state, ideal static-friction constraint, thermal
  partition or continuous-trajectory accuracy claim. Changing physical loads
  and slip directions still require temporal refinement.
- `accept_scenario_bc`: fs-scenario integration — `Dirichlet` and
  `Traction` under `Physics::Elasticity` accepted; dimensioned-value
  resolution stays with the scenario consumer.
- `shell` (bead frankensim-b8bxd.9.1): a safe-Rust, deterministic,
  flat, single-patch triangular **estimate-only surrogate** for the Euler-disc
  base. It combines lumped inertia, edge-spring in-plane resistance, and
  normal-slope/rotation penalties; it is not a CST plane-stress membrane or a
  Kirchhoff--Love/isotropic continuum discretisation.
  `ShellPlate::assemble` returns type-separated dense
  `ShellMatrix<Mass>`, `ShellMatrix<Stiffness>`, and optional Rayleigh
  `ShellMatrix<Damping>` operators in node-major `(x,y,z,rx,ry,rz)` order.
  Each triangle contributes edge-spring, normal-slope, and rotation-gradient
  surrogate energy; translational and rotary inertia are lumped from the
  material card. `ShellSupport` applies three distinct, normal-aligned point
  pinned supports by reducing all translations at those nodes while retaining
  rotations; its normal identifies the plate side and it does not claim a
  frictional or rotational support law. Inputs carry model,
  source, and state IDs and fixed SI/right-handed-Cartesian interpretation.
  Assembly bounds node, element, dense-memory, deterministic
  validation-and-assembly work, and diagnostic sizes before allocation; it
  returns no partial operator. Validation charges finite-node scans,
  scale-relative coincident-node pairs, canonical face insertion, oriented
  edge incidence/adjacency, connectivity traversal, used-node marking, and
  triangle assembly against the declared work budget. Small assemblies report
  raw algebraic mass-array definiteness, stiffness-array nullity, and symmetry
  residue via a deterministic Jacobi eigensweep. Its raw mixed-unit spectral
  spread is diagnostic only, not a physical condition number or continuum
  certificate; larger valid assemblies explicitly report that those bounded
  algebraic diagnostics were not computed.
  Scale-relative coincident nodes, duplicate index faces, geometrically
  duplicated faces through coincident nodes, non-manifold edges (more than two
  incident faces), same-orientation shared edges, and components disconnected
  by an edge are refused before assembly; callers requiring several components
  must model separate plates.

- `reduce` module (music bead 3ez8g.3.3): the offline lip/reed reduce
  lab. `reduce_valve(request, cx) -> ValveCard` runs
  `TetLinearElasticProblem` -> `fs_modal::slice_window` -> per-mode
  re-normalization from the solver's `phi^T M phi = 1` convention to
  UNIT MEAN TRANSLATION of the orifice face (mass normalization would
  destroy the units). The mean includes every declared face node, with
  constrained components contributing zero; duplicate face nodes refuse.
  This mints physical `m_eff = phi^T M phi` [kg],
  `k_eff = lambda m_eff` [N/m], and `c = eta sqrt(k m)` from the card's
  AUTHORED loss factor (the lab does not invent Q values). The
  retained-compliance law: the modal compliance
  `sum (phi^T r)^2 / lambda` of the retained modes — `phi^T r` is the
  FORCE-pattern projection; the first build used the base-excitation
  `phi^T M r` and retained 1e-11 of the compliance, an executed
  falsification — is compared against the EXACT static compliance
  `r^T K^{-1} r` (sparse LDLT), and a reduction below the authored
  floor REFUSES (`RetainedComplianceTooLow`). Orifice geometry is
  MEASURED from the mesh: the face is the surface strip facing the
  opposing plane, gap = closest approach, width/channel-thickness =
  face extents. Opening and width axes must be perpendicular (normalized
  dot-product tolerance 1e-12); nonfinite opposing-plane offsets refuse
  before assembly. An oblique frame is not silently projected or rescaled.
  The card carries the provenance chain (mesh digest,
  material state-point receipt identity, source id) under a
  domain-separated schema-versioned content identity with
  canonical-bytes round-trip and tamper refusal. Loading requires at least
  one mode, bounds mode allocation by available input rows, and refuses
  trailing records after the identity. `beating_reed()`
  requires an explicit effective pressure area `A` [m²], conjugate to the
  mean-opening coordinate, and mints `fs_scenario::BeatingReed` with
  `P_c = k H / A`. No fixed pressure-face length is assumed. Distributed
  pressure-load projection is available through `reduce_pressure_loaded_valve`;
  authored areas supplied directly to the plain card remain research inputs.
  Empty cards and invalid converted parameters refuse.
  It also passes the card loss as explicit reed damping
  `zeta = eta/2`, preserving `c = eta sqrt(k m)` in the playing dynamics.
  This is the existing reduced viscous model, not a broadband loss-law fit.
  `reduce_pressure_loaded_valve` accepts one finite mesh-axis force-per-pascal
  vector [m²] per node. It projects the load onto each retained mode using
  `A = phi^T load / mean_face(phi)`. This is the same work-conjugate coefficient
  for generalized force `p A` and swept-volume flow `A v`; constrained DOFs
  do no virtual work. The loaded surface need not equal the opening-observation
  face. The returned `PressureLoadedValve` retains the card, signed per-mode
  areas and an identity covering card/load/projections, and converts its first
  mode without a separately authored area. Nonpositive first-mode area refuses
  the positive-closing reed conversion instead of taking an absolute value.
  The caller owns consistent surface-load assembly, pressure orientation and
  units. This does not discover wetted surfaces, solve a pressure field, or
  certify pressure-load modal convergence from the opening-compliance fraction.
  The existing card serialization remains structural; it does not serialize
  this pressure port. The cane integration assembles uniform-pressure loads
  from its planar underside triangles and exercises the resulting bore response.

### Three-dimensional finite-strain hyperelasticity (`hyper3`)

`HyperTetProblem { nodes_m, tetrahedra, material, prescribed_m,
nodal_forces_n, body_force_n_m3, budget, settings }.solve(cx)` is the 3-D
Total-Lagrangian counterpart of `hyper2d`: P1 tetrahedra (constant `F` per
element, so one card evaluation is exact quadrature), the `fs-material`
energy cards (Neo-Hookean, Mooney–Rivlin, decoupled Ogden, HGO/GOH;
AD Piola stress and consistent 9x9 tangent, with coalescence-safe spectral
derivatives for Ogden), DEAD loads (nodal forces, body force per reference
volume lumped to vertices), and linear ramping over `load_steps`. Each
Newton step solves `K_ff d = -r_f - K_fp dp` so pending prescribed motion
enters through the consistent linearization; trial states that invert an
element are halved; once the boundary data are applied an Armijo condition
on the potential `sum V W - f_ext . u` is enforced (residual decrease on a
non-descent direction). The gate is `||r_f|| <= relative_tolerance *
max(||f_int||, ||f_ext||, max|K_dd| max|u|)`; the stiffness term keeps it
meaningful at stress-free states (a rigid rotation, where every force
vanishes). Linear solves: ILU(0)-PCG, falling back to |diag|-Jacobi PMINRES
for an indefinite tangent; acceptance on the recomputed residual.
`HyperTetSolution` carries displacements, reactions in `prescribed_m`
order, stored energy, the smallest `det F`, and per-step residual histories,
backtracks, linear iterations and fallbacks. `internal_force` and
`tangent_matrix` are exposed for verification. Element quality and sizes are
admitted against `TetAssemblyBudget`; refusals are `HyperTetError`
(`InvalidInput`, `DegenerateElement`, `MaterialRefused`, `NewtonStalled`,
`LinearSolveFailed`, `Cancelled`).

`HyperelasticLawField::Uniform` and `PerElement` use this same owner. A
per-element field follows tetrahedron connectivity order and must cover the
mesh exactly (`LawCountMismatch` otherwise). Each local law supplies every
trial energy, Piola stress and consistent tangent, including line-search and
final evaluations. Regional materials are not averaged. Fields may borrow
native `Hyperelastic` laws or `ResolvedHyperelastic` values retaining the exact
admitted constitutive model cards, including parameters, validity, source hashes
and provenance. Its explicit fixed external query must satisfy every source
bound and exact quantity descriptor. Static temperature, pressure, composition,
moisture and concentration conditions are supported; evolving mechanical axes
require their own owner and refuse. The context-free constructor only accepts
unconstrained source domains. The source wrapper delegates constitutive parameter
and native stretch checks to the existing hyperelastic owner; it introduces no
additional energy model.
The focused heterogeneous regression solves two source-defined Ogden layers in
series against an independent axial stress/energy formula and traction-balance
root, changes only the second layer's material data, and checks actual interface
motion and reaction. Uniform native/source/element fields agree bit for bit;
wrong cardinality and a single element's unsupported stretch refuse. Material
conditions remain fixed over the load path; thermal evolution, dynamic inertia
and source uncertainty propagation remain outside this static reduction.

Conformance (`tests/hyper3.rs`, all executed): the affine finite-strain
patch test is reproduced to 1.2e-14 m with the exact stored energy; the
compressible Neo-Hookean uniaxial stretch to 1.6 (rollers) matches the
closed-form lateral stretch and nominal stress to 1e-9 in 4-5 Newton
iterations per step; the tangent matches a central difference of the
internal force to 4e-11 (Neo-Hookean and Mooney–Rivlin); at the reference
state it equals the independently written `linear3` stiffness to 4e-16; a
rigid rotation of the boundary stores 1e-17 J with reactions below 1e-15 N;
inverted elements, repeated prescribed DOFs and cancellation refuse.

No-claims (`HyperTetProblem::solve`): static, isothermal, rate-independent; constant-strain
tetrahedra (volumetric locking of nearly incompressible cards is not
mitigated; no mixed or F-bar formulation); dead loads only (no follower
pressure); no contact; no limit-point continuation (drive by displacement
or refuse); no mesh-convergence claim from one run.

#### Complete-boundary current-area pressure (`hyper3::pressure`)

`ClosedTetPressureSurface::new` derives the entire exterior of the original
tetrahedral mesh through `TetLoadGeometry3`'s incidence, connectivity and quality
admission. It orients each exterior triangle outward and requires every
boundary edge to have exactly two oppositely directed incidences. Interior
faces receive no load. Duplicate/same-side tetrahedra, nonmanifold faces or
edges, open/inconsistent boundaries, unused vertices and work-budget overflow
refuse. The prepared owner retains exact reference-coordinate bits and ordered
connectivity for binding to a consumer's original geometry. Either reference
tet winding is admitted; every current tet must preserve its original sign.

For each outward triangle `(a,b,c)`, `force_n` places
`-p/6 * ((x_b-x_a) cross (x_c-x_a))` on each of its three nodes. Positive signed
pressure acts inward on the current area. The analytic `force_tangent` is the
external force Jacobian: perturbing a face vertex by `h` varies the cross
product by `h cross (x_b-x_c)`, `h cross (x_c-x_a)` or
`h cross (x_a-x_b)`. A consumer subtracts this Jacobian from its internal
residual tangent. On the complete closed boundary these assembled loads equal
`-p grad(V)` and the tangent equals `-p Hessian(V)`.

`volume_m3` sums current oriented tetrahedron volumes with finite arithmetic
and inversion/collapse refusal. `work_j` returns `-p (V_after - V_before)` for
pressure held between the two endpoints. This exact discrete volume difference
is not replaced by midpoint traction work; P1 volume is cubic in position.
Finite signed pressure, geometry and all returned arithmetic are admitted,
work is bounded by the reference node/tet/face limits, and cancellation is
checked during geometry and load operations. The owner changes no material,
reference mass or accepted trajectory and introduces no time solver.

`tests/hyper3_pressure.rs` specifies independent scalar-triple-product force
and volume oracles, finite-difference tangent and symmetry, rotated and
translated force/resultant/moment checks, a two-tet union with its interior
face removed, signed work, invalid/budget and cancellation refusal. Native
parsing passed. Full compilation and execution remain pending because the
shared executor killed the real `asupersync` prerequisite even during a
metadata-only check. No selected open patch, spatially varying pressure,
fluid equation of state, global injectivity certificate or continuous-time
inversion guarantee is claimed. Static `HyperTetProblem::solve` remains the
dead-load owner; finite-solid dynamics consumes this separate pressure kernel.

### Incremental small-strain inelasticity (`plastic3`)

`SmallStrainTetProblem { nodes_m, tetrahedra, laws, prescribed_m,
nodal_forces_n, load_path, budget, settings }.solve(cx)` follows a load path
on P1 tetrahedra for ANY `fs-material::SmallStrainLaw` (J2 radial return
with the Simo–Hughes algorithmic tangent, isotropic/orthotropic
elasticity). Tensor-component Voigt strain `eps = B u`, internal force
`sum V B^T W sigma`, tangent `sum V B^T W C_alg B` with `W = diag(1,1,1,2,2,2)`.
Each load factor (which may decrease or reverse) is one Newton solve from
the last converged state, prescribed increments through the consistent
linearization, residual-increase halving, and the same linear-solve policy
as `hyper3`. The static force gate uses element stress-force bounds on free
rows and actual free external loads, retaining the largest accepted scale
through each unloading increment. Fully fixed regions and rigid travel do
not weaken another region's equilibrium tolerance. Constitutive nonfinite
stress, force or tangent refuses before publishing state. Element states are updated by the law and
committed only after convergence. The solution carries displacements,
element strain/stress/states, and per-increment reactions, residual
histories, and counts of evolving elements. `internal_force` and
`tangent_matrix` take explicit committed states for verification.
The convergence check precedes tangent assembly. The element tangent kernel
precomputes its twelve constitutive strain-column contractions once on the
stack and reuses them for all test columns; no changed convergence order,
speedup factor or mesh-accuracy claim is inferred from this operation reduction.

`SmallStrainLawField::Uniform(&law)` and
`SmallStrainLawField::PerElement(&laws)` use this same numerical owner.
Per-element entries follow tetrahedron connectivity order and may come from
different material cards or spatially varying resolved conditions. Each
tetrahedron initializes its own law's state and uses that law for every trial
stress, consistent tangent, converged history update, and final stress.
A wrong field cardinality returns `LawCountMismatch` before assembly.
No effective-property averaging or separate heterogeneous solver is introduced.
Constitutive coefficients remain fixed throughout the load path.
Reported endpoint stress is evaluated from the same previous committed history
as the converged force and tangent, before publishing the next internal state.
An explicit-duration overstress law therefore advances exactly once per load
increment; reading the result does not perform an extra relaxation interval.
The same publication rule serves the backward-Euler inertial adapter.
`ResolvedJ2Plasticity` implements the native `SmallStrainLaw` interface directly,
so fields can borrow complete resolved material inputs and retain all source
receipts without extracting bare coefficients. This field alone does not add
thermal feedback, state remapping, or finite-strain mechanics.

The homogeneous numerical measurements below predate this field extension.

The source-selected `MixedPlasticStrainModuli` J2 extension uses this same
P1 owner and retained per-element plastic strain. Its Prager backstress changes
reverse yielding while the supplied isotropic and kinematic moduli can retain
the same monotonic tangent. A focused authored regression in
`tests/plastic3_material.rs` drives tension, unloading and reversal on the
original six-tetrahedron coupon and compares reactions and accumulated/signed
plastic strains with independent uniaxial equations. This introduces no
alternate solver, fatigue-life estimate or finite-strain plasticity claim;
the test's presence is not a claim of native execution.

Conformance (`tests/plastic3.rs`, executed): homogeneous uniaxial tension
with linear hardening (E = 1000, nu = 0.3, sigma_y = 1, H = 100) follows the
closed form through yield to eps = 0.006 (force 1.4545454547 vs
1.4545454545) and elastic unloading to the plastic strain leaves zero force
(2.6e-11), alpha = eps_p and the isochoric lateral contraction eps_p / 2 to
1e-8 relative, in 2-5 Newton iterations per increment; the algorithmic
tangent matches central differences of the internal force to 4.6e-10 at a
plastic state; the elastic-law tangent equals `linear3` to 4e-13 relative; a
cantilever driven by tip displacement yields progressively (1, 20, 33 newly
yielding elements) with a falling secant stiffness.

No-claims (`plastic3`): small strain and rotation; no implicit clock
(rate-dependent laws require explicit physical intervals);
constant-strain tetrahedra lock under confined plastic incompressibility
(confined limit loads are over-predicted; no B-bar/mixed formulation); dead
loads; no limit-load continuation; no mesh-convergence claim.

### Sub-yield creep on the native solid (`plastic3::creep`)

`SmallStrainTetProblem::solve_creep_from` consumes the accepted nodal
displacement and per-element `CreepState` through the original sparse Newton
owner. Uniform `NortonCreep` and source-bound `ResolvedNortonCreep` fields
retain the same constitutive flow, tangent, source receipts and reference
geometry. Exactly one absolute load endpoint represents one positive common
physical duration; the adapter refuses multi-entry paths or unequal element
durations. Newton and backtracking always borrow the old history. The published
stress is the converged trial stress evaluated before committing that history,
so observing a rate-dependent endpoint cannot advance another interval.

The adapter admits finite old/new material histories and bounds the full
displacement gradient by each law's small-strain limit, including its skew
part. Each checked constitutive increment is compared with the actual native
stress and committed state. Per-element accounts in J/m³ and reference
volumes are returned for conservative physical-heat transfer. The whole-solid
audit independently compares applied-force and support-reaction endpoint work
with integrated stress work, then with elastic storage change, physical creep
dissipation and backward-Euler numerical loss. Numerical loss is not heat.
All local and global gates and cancellation precede publication; caller-owned
displacements and histories are never mutated.

Elastic coefficients and geometry must remain fixed across an interval.
Prescribed temperature may change Norton kinetics between intervals while
retaining creep history; elastic parameter work, thermal expansion, reciprocal
thermal evolution and inertia require their separate owners. This is a
small-strain, secondary-creep model, with no primary/tertiary creep, rupture,
confined-flow accuracy or material-population qualification claim.

`tests/plastic3_creep.rs` checks a loaded six-tetrahedron bar against analytical
constant-stress Norton creep, permanent set after unloading, held-length
relaxation against the backward-Euler recurrence and three timestep rungs,
per-element temperature-dependent histories and heat, and transactional
refusals including cancellation. These are manufactured implementation checks,
not measured validation.
The source-bound regression warms an accepted loaded trajectory, preserves
every pinned claim, and substitutes a material with twice the reference rate;
the actual bar extension follows the independent accumulated creep law.

### Simultaneous creep and yielding on native tetrahedra

`SmallStrainTetProblem::solve_creep_from` consumes the material owner's
`CreepLaw` and preserves its complete associated history. The same native
Newton solve, small-strain/rotation checks, endpoint stress/state comparison,
support/applied work and physical/numerical energy audit serve Norton creep
and simultaneous isotropic Norton/J2. `CreepTetStep` defaults to the original
creep state; combined steps retain separate creep and plastic memories.
Per-element accounts expose both physical heat channels and their total;
energetic hardening storage and backward-Euler loss remain excluded from heat.
The single physical interval, bounded work and immutable proposal semantics
are unchanged. `tests/plastic3_creep_plastic.rs` exercises an actual loaded
tetrahedron yielding and continuing to creep, unloading, held-shape relaxation,
analytical deformation/energy and preservation of both histories on refusal.
Finite deformation, inertia and noncoaxial kinematic hardening remain outside
this checked creep path.

### Permanent creep in native fiber members and frames

`FiberLaw::Norton` and `FiberState::Norton` consume the original material
owner's free-sided axial reduction, with signed permanent strain, accumulated
creep and optional original source receipts. The existing section integration,
force-based member compatibility and global frame Newton solver remain shared.
No effective modulus, alternate beam solver or inferred creep coefficient is
introduced. Explicit authored reference density in `FiberMaterial::norton`
does not acquire source authority from the law's separate coefficient receipts.

`Frame2d::solve_creep_step` admits one absolute load endpoint and one common
positive physical duration across the original Norton fiber points. Section
Newton, compatibility retries and global backtracking see accepted history;
the returned frame commits once. Stored forces and `creep_fiber` observations
do not evaluate another physical interval. `with_creep_laws` refreshes time
and prescribed temperature while retaining native geometry, history, elastic
and kinetic parameters, original card/claim identities and nonthermal context.

The existing untimed static, pushover and modal routes, and the inertial
Newmark route, refuse Norton fibers. They do not own the explicit creep clock
or its checked work account; an algorithmic creep tangent cannot silently
become an instantaneous modal stiffness or advance on a mismatched time step.

The checked step integrates each constitutive energy account using its original
fiber area, Lobatto weight and member length. Applied/support work is checked
independently against this material work and against storage change plus
physical creep loss and numerical loss. Force and moment gates use separate SI
tolerances. Physical heat and point volumes remain available to thermal
consumers; numerical loss is not heat. Immutable proposals, explicit fiber-work
budgets and cancellation checks preserve accepted history on refusal.

This path admits stationary inertia, linear geometry and Norton or joint
Norton/J2 sections.
It bounds nodal rotations, chord displacement and absolute quadrature curvature;
unsampled spatial extrema require refinement. Corotational geometric work,
inertia, thermal expansion, reciprocal temperature evolution, unsupported
inelastic fiber laws and measured-life qualification are outside this checked account.
`tests/frame_creep.rs` exercises constant-force extension/permanent set,
held-length tension relaxation, actual cantilever curvature/tip shape under
sustained moment, source-temperature continuation and atomic refusal/cancellation.
It also verifies that the untimed and inertial entry points cannot age fibers.

### Rod couples held in laboratory coordinates

`rod::WorldTipLoad` declares a world-held tip force in newtons and a free couple
about that moving tip in newton metres. `Rod`, `MaterialRod` and
`SpatialMaterialRod::solve_static_world` use the original energy, Lie-group
Newton, load stepping and residual backtracking. Each trial evaluates
`M_body = R_tip^T M_world`; the exact external contribution to the body/right
residual tangent is `-hat(M_body)`. It is included alongside both native
analytical internal tangents and the arbitrary-energy finite-difference fallback.
The constant couple is not replaced by the gradient of `-M dot Log(R)`.

`TipLoad` continues to express a world-held force plus a body-held follower
couple. These are distinct physical load laws when a tip turns. World loading
composes with the existing source-resolved shear/torsional stiffness, spatial
materials, thermal sections, intrinsic geometry, gravity, tensile links and
fixed supports. `support_reactions_world` uses the identical load transform at
the accepted pose. The material owners retain atomic publication, source
admission and existing cancellation boundaries; no separate rod solver,
material averaging or load-induced change of reference mass is introduced.

G1/G3 tests check the load tangent against independent rotation differences,
laboratory-frame covariance, actual whole-rod force/torque balance after finite
bending, distinction from follower loading, fixed-support frame conversion,
source-bound torsional compliance and temperature-dependent shear modulus.
These remain bounded static equilibrium results. A residual is not a stability
or continuum-error certificate; no dynamic trajectory, accumulated torque work,
general conservative torque potential or material calibration is claimed.

### Intrinsically curved material rods

`SpatialMaterialRod::try_from_reference` binds `RodReferenceGeometry` to the
existing per-segment material sections. Both present and unloaded strain measures
use the same midpoint directors/SO(3) logarithm; elastic energy and resultants use
their difference. A cambered or twisted unloaded part therefore carries no
invented bending/torsional preload. The original Newton, support reactions,
reference-mass gravity and local thermal-section owners remain authoritative.
`RodMaterialSource::with_reference_geometry` preserves this datum, source pins,
and mass through temperature refresh; gravity uses the actual reference nodes.

Inputs declare an unloaded shape with equal material-coordinate intervals and
resolved neighboring rotations below pi. A loaded scan cannot supply unloaded
geometry without an independently justified reconstruction. Refine nodes,
frames and section assignments together. This remains bounded static slender-rod
mechanics with circular concentric sections and finite-difference tangents;
it does not supply dynamic stick/ribbon coupling, scan accuracy, contact,
self-intersection certification, or experimentally validated instrument audio.
G0/G1/G3 tests in `rod_material_field/reference.rs` cover intrinsic zero stress,
a three-rung circular-arc end-moment limit, material changes, gravity datum and
objectivity; `rod_material_source.rs` covers source/temperature continuity.

### Source-bound tensile members on material rods

`SpatialMaterialRod::with_tension_links` attaches up to 64 immutable
`MaterialTensionLink` bundles to existing nodes, including body-frame offsets.
Each keeps its original `ConcentricElasticSection`, natural length, strand count
and explicit tensile-strain limit. At the actual deformed endpoint distance L,
the complete bundle stores EA/(2 L0) max(L-L0,0)^2; slack links exert no compressive
force. The original rod energy/Newton/support-reaction owner sees this energy
once. Rotating offsets therefore transmit reciprocal forces and moments through
the same equilibrium, including on intrinsically cambered rods. Public link
observations report span, tension, storage and endpoint force/moment pullbacks.

The source-bound reference mass uses natural length, not the stretched span.
Explicit gravity includes half that mass at each endpoint and the corresponding
rotating offset potential and moment. Replacing a link rebuilds these weights
from retained inventories. `reference_mass()` remains the rod-only receipt;
`total_reference_mass_kg()` includes links and each link retains its own source.
`gravity_nodal_forces_n()` and `gravity_nodal_moments_body_nm()` expose the
matching loads. The separate link is held at its explicitly resolved material
state; rod heating does not silently refresh its properties or natural length.

This is bounded quasistatic straight-link mechanics with tensile linear strain,
not a sagging cable, dynamic bow ribbon, thermal link, rosin/contact law, damage,
or validated instrument model. Nonstatic material claims, invalid attachments,
unrepresentable inventory and trial strains outside the declared domain refuse.
G1/G3 tests in `rod_material_field/tension.rs` cover analytic parallel compliance,
energy-derived offset reactions, camber/material preload sensitivity, conserved
mass, gravity work, slack and refusal preservation. The existing
`fs-couple` `heated_material_rod` example accepts an explicit `--tension-link`
natural-length/offset/strand-count record and reports its physical observations;
its separately cold 300 K member is not a horsehair material preset.

### Shared material fluid loading on submerged rods

`SpatialMaterialRod::with_submerged_fluid` consumes the existing
`fs-material::FluidTransport` interface and the specimen's explicit gravity.
The uniform hydrostatic pressure resultant is `-rho_fluid V g`; its reference
volume uses the original layered section geometry rather than an effective
material density. `ReferenceDisplacementVolume` explicitly distinguishes flooded
central voids (annular material volume) from sealed cavities (outer envelope).
Solid mass, density receipts and native motion inertia remain unchanged.

The original shape-function integration retains total displaced volume and its
first moment, including clamped nodes. Attached strands use their natural-length
volume and strand count. Their rotating offsets contribute a separate buoyant
potential, force, moment and analytical tangent through the same rod Newton,
support-reaction and motion owners. Replacing fluid, gravity or links rebuilds
the load without accumulation; zero gravity retains an explicitly bound fluid.
The original material sources retain measurement authority; the shared numerical
transport snapshot does not invent a source receipt.

All represented bodies must be fully immersed with nonoverlapping displaced
volumes. This is an explicit fixed-reference-volume approximation: current
solid density does not redefine mass or displacement volume during thermal
resolution. No partial immersion, free surface, radial expansion, fluid
inventory in sealed cavities, capillarity, drag, added mass or compressible
hydrostatic profile is inferred. Material-fluid replacement is static specimen
construction, not a conservative time-dependent environment update.
Focused `rod_gravity` regressions cover density-driven equilibrium and neutral
buoyancy, layered/hollow volume and first moments, gas-pressure changes, thermal
reference inventory, eccentric strand work/moments, motion inertia, and atomic
load replacement/cancellation. These synthetic checks make no measured-fluid
or experimentally validated material-response claim.

## Invariants

1. Patch tests exact: linear displacement fields are reproduced to
   solver tolerance on P1, Q1, B-bar, and the CutFEM frontend
   (all-embedded fixture) — sol-001.
2. MMS optimal orders per implemented family × frontend: L2 ≈ 2,
   H1 ≈ 1 for P1/Q1/CutFEM-Q1 (sol-002).
3. Objectivity through the fs-material interface: W(RF) = W(F) and
   P(RF) = R·P(F) for both cards (sol-003).
4. Consistency merge gate: the assembled tangent equals the FD
   directional derivative of the residual, and the residual is the
   exact gradient of the potential energy (rel ≤ 1e-6) — sol-004.
5. Locking is MEASURED, not assumed: on thin near-incompressible
   bending, the standard element's error degrades by ≥10× from
   ν = 0.3 to ν = 0.4999 while B-bar stays within 2× (sol-005).
6. Cook's membrane envelope: converged tip at the literature reference
   point (48, 52) inside [23.5, 24.5] (plane stress, ν = 1/3), ≤2%
   self-convergence deviation; near-incompressible B-bar
   self-converges while the standard element's shortfall is logged
   (sol-006).
7. Newton robustness: large-strain bending and a buckling-adjacent
   compressed strip converge under load stepping + line search with
   fast terminal contraction; histories and backtracks are evidence
   (sol-007).
8. Determinism: BTree/insertion-ordered assembly, deterministic
   Krylov solvers; bit-identical across runs on a platform —
   including continuation trajectories (checkpoint/resume proven
   bitwise-equal in stab-002).
9. Buckling: the Euler strut's EXTRAPOLATED critical load lands
   within 3% of the analytic value with the Richardson indicator
   covering the raw fine-mesh gap (Q1 parasitic-shear inflation,
   reported, stab-001); pencil symmetry to 1e-10 and mode
   K-orthogonality to 1e-6.
10. Continuation: validated against the CLOSED-FORM von Mises truss
    (limit points within 5–10% of analytic, manifold deviation
    ≤ 1e-7, load control jumps where arclength traces, stab-002);
    the continuum tent traces the full snap-through Z-curve on one
    path — snap and snap-back limit points, load recovery beyond the
    load-control failure point (stab-003).
11. Branch handling: bifurcation detected within 15% of the pencil
    prediction on the compressed column; null-direction switching
    lands on the bent branch (transverse growth ≥ 20× measured,
    stab-004).
12. Flat plate-surrogate energy consistency: each assembled stiffness contribution is a
    positive weighted outer product, so the free operator is symmetric
    positive-semidefinite while positive lumped translational/rotary inertia
    makes the raw mass array positive-definite. A single free triangular panel has the six
    infinitesimal rigid motions in its stiffness null space; the focused
    `shell_operator` battery verifies this and independently reconstructs the
    quadratic energy from `Kx`.

## Error model

`SolidError` teaching errors: `SolveFailed` (linear gate missed),
`NewtonStalled` (carries the residual history; repair = more load
steps), `MaterialRefused` (det F ≤ 0 escaped globalization; repair =
smaller steps), `UnknownPatch`, `UnsupportedBc` (fs-scenario mapping
outside Dirichlet/Traction × Elasticity), and `InvalidInput`
(non-finite/inconsistent caller data such as invalid contact
parameters, infeasible starts, or malformed adjoint seeds).
`InternalInvariant` distinguishes a corrupted lower-level solver graph from
repairable caller input.

The CutFEM facade maps only `CutFemError::SolveNotConverged` to
`SolidError::SolveFailed`, preserving its iteration count and relative
residual. A canonical `ConstraintCycle` maps to `InternalInvariant`; every
material, parameter, geometry, callback, topology-admission, and
unsupported-regime refusal, including an invalid typed edge band or a cut
supported segment, maps to `SolidError::InvalidInput` with the canonical
diagnostic. No canonical refusal is converted into a solution.

## Determinism class

Bit-deterministic across runs on a fixed platform (ordered assembly,
fs-solver CG/MINRES, no threading, no ambient state). CutFEM determinism
and topology ordering are inherited directly from fs-cutfem rather than
duplicated in this crate. Cross-ISA golden hashes are not yet recorded
(follow-up).

## Cancellation behavior

Bounded synchronous loops (assembly, capped Krylov iterations, capped
Newton with capped backtracking). Chunked Cx polling belongs to the
fs-exec driver (the L3 discipline; fs-feec/fs-cutfem precedent).

## Unsafe boundary

`#![deny(unsafe_code)]` via workspace lints; no capsules.

## Feature flags

`koiter-asymptotics` (off by default): the plan-flagged [F] Koiter
post-buckling module and its battery, per the Ambition-Tag gating
rule.

`contact` (off by default): the plan-flagged [F] SDF-native barrier
contact prototype and its battery, per the Ambition-Tag gating rule.

### Structural invariants (tfz.14)

12. Rod objectivity: superposed rigid motion leaves strain energy
    invariant to 1e-12 (str-001).
13. Elastica: large-deflection tip matches an in-test SHOOTING oracle
    (RK4 + bisection on the elastica BVP) within 0.015·L at
    PL²/EI = 2 (str-002).
14. Helical family: pure bending gives the EI/M circle (2%), pure
    torsion stays straight (1e-6) at rate M/GJ (2%), combined moments
    give the constant helix strain state (3%) — str-003.
15. RC hysteresis: per-cycle dissipation positive and growing with
    amplitude, peak moment inside the hand capacity band, section
    states bitwise-deterministic across runs (str-004).
16. Batched section updates are scalar-consistent to 1e-8 with
    ledgered throughput (str-005).
17. Force-based pushover softens past yield (secant < 0.5× initial),
    dissipates under reversal, and resumes from a checkpoint
    BITWISE (G4) — str-006.
18. Contact barrier calculus and solve hygiene (feature `contact`):
    barrier derivatives match finite differences, invalid/penetrating
    contact inputs return structured `InvalidInput` instead of
    panicking, accepted iterates keep positive gap, reactions balance
    the applied normal load, tight-clearance and curved-obstacle
    fixtures retain positive gap, lagged friction separates stick/slip
    cases, and the frictionless translation adjoint matches FD on a
    plane — cnt-001..007.

11. Reduce-lab gates (rl-001..rl-004): the thick cantilever-slab
    oracle (f above the slender Euler-Bernoulli line by the plate-strip
    factor, measured +11.4% inside an authored 0-35% band;
    strip-normalized effective-mass ratio 0.355 in (0.30, 0.42);
    `k = lambda m` at machine precision; the
    `k_ratio = f_ratio^2 * mass_ratio/0.25` self-consistency identity;
    retained compliance 0.98 > 0.95), the minted card driving the
    existing fs-phs valve island to its quasi-static equilibrium at
    1e-3 with supply defect < 1e-8, refusals (empty/out-of-range face,
    the retained-compliance floor via a window that excludes the
    compliant fundamental, a face behind the opposing plane), and the
    canonical-bytes round-trip + tamper refusal.

## Conformance tests

`tests/battery.rs`: sol-001 patch tests (body-fitted max nodal error
≤ 1e-9; CutFEM at the CG-tolerance floor ≤ 1e-7); sol-002 MMS orders
(P1 2.00/1.00, Q1 2.00/1.00, CutFEM 2.19/1.04); sol-003 objectivity
(energy rel ≤ 1e-10, stress rotation ≤ 1e-8, 8 random states × 2
cards); sol-004 tangent/energy consistency (≤ 1e-6); sol-005 locking
battery (standard degrades ≥10×, B-bar ≤2×); sol-006 Cook's envelope;
sol-007 Newton robustness (5-step large bending, 8-step compression,
terminal contraction ≥ 100×); the selection-policy probe. Unit tests:
mapped-mesh areas, selection thresholds. The CutFEM ownership probes
compare facade and direct canonical coefficients, nodal values,
residuals, iteration counts, active cells, and error norms bitwise after
the documented stabilization translation. A typed right-edge band likewise
compares coefficients, nodal values, and assembled-load compliance bitwise and
proves that an occupied legacy callback slot refuses. The probes also require
invalid material data, the uncertified near-incompressible regime, and an
empty domain to fail closed as `InvalidInput`.

`tests/stability.rs` (bead tfz.15): stab-001 Euler pencil (analytic
G2 + Richardson indicator); stab-002 von Mises truss closed-form
continuation oracle (+ bitwise resume); stab-003 continuum tent
snap-through with the self-calibrating load-control-failure probe;
stab-004 branch detection + null-direction switching; stab-005
eigenvalue-derivative gradient gate (rel ≤ 1e-3 vs frozen-K_G FD) and
the clustered-eigenvalue trap (min() kink ≈ 2 measured, KS aggregate
smooth with derivative matching FD to 1e-6, conservative).
`tests/koiter.rs` (feature-gated): stab-006 symmetric-stable
classification + imperfect-geometry sampled oracle.

`tests/contact.rs` (feature-gated): cnt-001 barrier calculus and
force divergence; cnt-002/003 intersection-free hard press and
reaction balance; cnt-004 tight-clearance and curved-obstacle gap
audits; cnt-005 lagged friction cone separation; cnt-006
frictionless contact translation adjoint versus FD; cnt-007
structured invalid-input returns; cnt-008 fixed-reference friction balance;
cnt-009 unpinned conveyor transport, common-translation invariance, motion
reversal and signed endpoint work closure, zero-motion equivalence, and invalid
surface-displacement refusal (G1/G3, force/translation tolerance 1e-9 and
endpoint work balance 1e-14); cnt-010 fixed-duration 1/2/4/8-step friction
against an independent force-balance creep oracle (1e-10), timed conveyor
reversal and cumulative work closure (1e-14), held-state zero work and
prescribed boundary work;
cnt-011 timed/default equivalence, invalid time/speed/product admission,
representable-clock progression and atomic nonlinear-failure behavior.

`tests/structural.rs` (bead tfz.14): str-001 objectivity; str-002
elastica vs shooting oracle; str-003 circle/twist/helix; str-004 RC
moment-curvature hysteresis + determinism; str-005 batched
consistency + throughput ledger; str-006 force-based pushover,
reversal dissipation, G4 bitwise resume.

`tests/frame2d.rs`: a single elastic element reproduces the Euler–Bernoulli
cantilever tip deflection, rotation and axial shortening to 1e-12; a
fixed-base portal matches the closed-form lateral stiffness
`(24EI_c/h³)(6ρ+1)/(6ρ+4)` to 1e-6; the corotational cantilever under tip
moment `πEI/L` rolls into the exact half circle (tip within 0.5% L, rotation
π to 1e-9); tip-mass frequency `√(3EI/mL³)` to 1e-6 and its softening
`ω²/ω₀² ≈ 1 − P/P_cr` at half the Euler load; RC pushover with exact
`Mᵢ = H·L` equilibrium bounded by an independent section moment envelope and
base-concentrated curvature; cyclic pushover dissipates energy over closed
cycles; sustained-gravity equilibrium for linear and corotational members,
loaded/unloaded linear superposition with the exact sustained-load work,
continued-state kinetic-energy accounting, and compressed-column vibration
against the gravity-loaded tangent oscillator; elastic SDOF Newmark against
the exact Duhamel response (< 0.1% of
peak) with energy balance < 1e-6; a 3-story 2-bay RC fiber frame with gravity
preload, modal Rayleigh damping and an 8 s synthetic record (energy balance <
2%); bit-identical replay; malformed-input refusals.

## No-claim boundaries

- `frame2d` members are planar Euler–Bernoulli (no shear deformation, torsion,
  warping or out-of-plane buckling); reference mass is lumped by default or
  consistently integrated when explicitly selected; no
  shear failure, bond-slip, joint-panel flexibility, P-Δ of leaning columns,
  or soil–structure interaction; the corotational formulation is large
  rotation / small strain. The dense LU global solve is sized for building
  frames of a few hundred DOFs, not for sparse large models. Fiber-section
  response inherits fs-material's Mander / Menegotto–Pinto cards and their
  no-claims.

- `linear3` claims infinitesimal-strain elasticity and work-conjugate thermal
  eigenstrain loading on a fixed, conforming, non-degenerate tetrahedral mesh.
  The thermal operator does not infer expansion coefficients, solve a transient
  response, remesh, establish mesh/stress convergence, or claim finite-strain
  validity. Its coordinate update is fixed-topology and inversion-gated, not a
  finite-strain constitutive update. Its static solve establishes only its reported linear
  residual on the supplied discretization. It does not claim plasticity,
  physical damping, phase change, fracture, moving interfaces, or topology
  evolution. Callers must escalate or refuse when the resolved material state
  leaves this solid rung.
- Its free-body stiffness correctly contains rigid modes. Removing, fixing,
  or slicing around those modes is a boundary-condition/eigenproblem decision;
  the assembler never adds hidden grounding stiffness.

- The facade inherits fs-cutfem's no-claim for clipping a supported traction
  segment cut by the SDF and for quadrature error bounds on arbitrary traction
  callbacks. It only delegates checked support; it does not certify callback
  smoothness or grow a second boundary assembler.
- Extreme `(E, beta, gamma)` combinations whose active physical stabilization
  products are finite but whose fs-cutfem dimensionless coefficients are not
  representable are not claimed. Nitsche-active and ghost-active terms fail
  closed pending a canonical physical-coefficient API; an unused Nitsche value
  in traction-free mode is ignored rather than validated.
- TDNNS-proper / weakly-imposed-symmetry FEEC stress elements: awaits
  the simplicial H(curl)/H(div) family bead (dcng); B-bar is the
  shipped, measured locking-free path with the same acceptance metric.
- 3D and shells; higher-order families (Qk/Pk, k ≥ 2); IGA frontend
  (bead tfz.9).
- Hyperelasticity on the CutFEM frontend (linear ships there; the
  finite-strain cut path composes this crate's Newton loop with
  fs-cutfem quadrature in a successor bead).
- The legacy CutFEM facade makes no independent coercivity,
  near-incompressible, quadrature, graded-tree, convergence, or
  cancellation claim. In particular, material ratios beyond
  fs-cutfem's certified compressible regime are refused rather than
  routed to fs-solid's body-fitted B-bar formulation.
- fs-opdsl-generated residual/adjoint paths (the DSL coverage bead);
  the hand path here passes the same consistency gates the DSL output
  must pass.
- Contact no-claims: mesh-mesh point-triangle fallback, rigorous
  interval CCD over curved SDFs, adaptive barrier-stiffness policy,
  frictional adjoints/VJPs, dynamics/time stepping, and fs-scenario
  contact-law coefficient wiring. Plasticity flow remains successor
  scope (fs-material's plastic module wires in later).
- Ogden cards use the existing HyperTetProblem material interface unchanged;
  `tests/hyper3.rs::ogden_card_solves_uniaxial_bar` compares reactions and
  lateral motion against an independent compressible material-point root on
  three affine P1 mesh rungs (1e-8 relative reaction band). This establishes
  card consumption and the affine patch response, with no P1 locking claim
  or intrinsic soft-tissue calibration.
- Rod SE(3) DYNAMICS (fs-time symplectic wiring) and analytic rod
  tangents; 3D frames/lattices assembled from rods (topo-truss's
  consumer pass); J2 continuum element wiring (the fs-material card
  ships with algorithmic moduli; the element battery is successor
  scope); published El Centro fiber-model plot digitization (the
  shipped envelope is the capacity band + self-consistency).
- fs-scenario dimensioned-value resolution (units plumbing stays with
  the scenario consumer; kind/physics validation ships here).
- Production-scale pencil solves (shift-invert/sparse factorizations
  beyond the fixture-gated dense reduction) and the
  prebuckling-adjoint chain dλ/ds through σ₀(s) (fs-ad/ASCENT
  integration successor); LOBPCG preconditioning hooks.
- Plate/shell buckling references (needs bending elements — the
  shells bead); the 2D canonical here is the strut.
- CST plane-stress patch accuracy, Kirchhoff--Love/isotropic-continuum claims,
  physical condition-number authority, curved shells, multi-patch/mortar joins,
  NURBS geometry, transverse-shear
  (thick/Reissner--Mindlin) behavior, nonlinear/finite-rotation response,
  loads, time integration, constraint multipliers for oblique support normals,
  shell error enclosures, and NAFEMS qualification. The `shell` operator is a
  flat, linear spatial base only; fs-iga's dedicated surface-shell frontend
  remains the required successor for those claims.
- Koiter coefficients versus LITERATURE tables (v1 validates
  classification + the sampled oracle; quantitative a/b benchmarking
  awaits the shell fixtures where the classic tables live).
- The reduce lab is Estimate-authority end to end until measured
  tissue/cane packs exist: the material card in the battery is
  synthetic-authored (no soft-tissue elastic pack exists in
  data/matdb — the recorded upgrade path), the loss factor is an
  authored card field, and the "lip-shaped" fixture is an authored
  tapered slab. The lab claims reduction correctness, never tissue
  truth; claim-side promotion happens in the per-track gates beads.

### Continuation of accepted 3D inelastic history

`SmallStrainTetProblem::solve_from` starts the original native increment/Newton
loop from supplied accepted nodal displacements and element histories.
`solve` supplies virgin history to that same implementation. Every nonlinear
trial uses the previous converged element history, and the method returns a
complete immutable candidate only after the requested load path converges.
Wrong array sizes and non-finite displacements refuse before history copying;
cancellation or a later failed increment leaves the supplied state reusable.

The new problem's load factors are absolute multipliers of its own prescribed
displacements and dead forces. A continuation may release constraints, change
loads or rebind per-element coefficients. The caller owns mesh/material lineage
and law-specific history and parameter-work admission: this generic numerical
entry point neither transfers constitutive families nor certifies a thermal
coefficient change. It retains the existing small-strain and P1 limitations.

Focused G1/G3 regression sources cover split-path equivalence through reversal
and reloading, closed-form permanent strain after releasing the driven end,
elementwise yield reduction at fixed displacement, and retry after a later
increment refusal without consuming the accepted history.

## Quadratic thermal prestress and vibration

`linear3::QuadraticPrestressedAssembly` extends the reference initial-stress
approximation to TET10, HEX20 and HEX27. Each problem's
`solve_thermal_tangent` or `solve_thermal_tangent_with_nodal_forces` first solves
its source-bound thermal/dead-load equilibrium and returns that immutable response
with the existing `TetElasticAssembly` vibration pencil. An independently supplied
stress array or stale solution cannot enter this API.

The correction is
`K_g[(a,i),(b,j)] = delta_ij integral grad(N_a)^T sigma_0 grad(N_b) dV`.
The same isoparametric map and stiffness quadrature evaluate the actual
`sigma_0 = C(Bu - epsilon_free)` at every point, through the shared response
kernel. Element-average stress is not substituted into this product. The
original reference consistent mass, geometry, material state and free-DOF order
remain unchanged. The pencil `(K_e + K_g, M)` feeds the existing modal/spectral
owner directly; its original elastic and signed geometric matrices remain
separately inspectable.

Tension stiffens and compression softens the relevant incremental motions.
Neither negative stiffness nor nonpositive eigenvalues may be clamped,
absolutized or presented as real vibration frequencies. This symmetric pencil
does not itself certify stability. It omits finite-strain kinematics, follower
loads, geometry-dependent boundary loads, live modal state transfer and
post-buckling. Fully fixed bodies have a valid static thermal response but no
vibration pencil and return `NoFreeDofs`.

Existing assembly, source-state, rigid-support, residual and pointwise
small-strain gates precede publication. Element/quadrature and sparse-row
checkpoints preserve cancellation. Identity binds the admitted equilibrium and
the complete signed stiffness/mass pencil.

G1/G3/G4 definitions in `tests/quadratic_thermal_tangent.rs` use all three
families to check zero-stress recovery, signed uniform thermal prestress,
unchanged mass/support maps, spectral shifts, and a manufactured varying-stress
quadratic field whose geometric energy differs from an element-average shortcut.
They also check the fully constrained and cancellation boundaries. These are
numerical approximation checks, not experimental material validation.

## Generic initial strain on the native tetrahedral owner

`linear3::eigenstrain` admits source-bound mesh-frame Mandel initial strain
without manufacturing a temperature. Principal strains rotate through the
existing second-order Mandel transform with explicit nonzero source/mesh frames.
Generic and thermal inputs share the original equivalent-load, PCG displacement,
stress recovery, geometry admission and prestress kernels; thermal source types
and thermal identity formats remain intact. Generic output carriers now name
the shared eigenstrain response; existing thermal names denote those same results.
`SolidStiffnessModel` is the existing source schema and lowering shared across
thermal and other evolving material conditions. No duplicate constitutive or
linear solver is introduced. Deformed charts expose their admitted tetrahedral
volumes for mass-conserving current density.

Generic initial strain is checked against exact current elastic identities;
stress/tangent recovery refuses a stale material/source/load/solution. Fully
constrained equilibrium returns zero displacement and actual reactions; a modal
pencil still requires free DOFs. The tangent uses the existing signed
reference-configuration initial-stress approximation, not a finite-deformation
or post-buckling theory. Input source work is external to this mechanical energy
account; no combined chemo/thermo-mechanical passivity claim is added. G1 tests
cover free directional deformation, constrained stress and energy, thermal/generic
numerical equivalence, frame rotation, stale input and cancellation.


## Source viscoelastic pressure-valve dynamics

The material-valve path in fs-solid::reduce::visco connects the existing
TetViscoField / TensorProny owner to the ordinary valve reduction, pressure
ports and BeatingReedSource interface. Uniform or per-element anisotropic
Kelvin-Voigt/Prony data supplies density, storage and loss. The same TET4 or
explicit-quadrature TET10 mesh, constraints, modal shapes, face normalization
and pressure work projection are used throughout. No material-name branch or
independent damping constant is inserted.

This is an explicitly admitted frozen-frequency reduction:
K = K'(omega_ref), C = K''(omega_ref)/omega_ref. All modes share ONE storage
pencil and its mass-orthogonal basis. The existing fs-modal damping certifier
examines the full projected symmetric C; the legacy `MaterialValve` refuses
a required coupled block instead of exporting independent decay rates. Full retained storage
and loss matrices are also compared with the frozen model at EVERY retained
modal frequency. Controls explicitly bound the resulting storage-energy and
absolute loss-factor operator errors, the modal loss factor and allocation.
These are finite-subspace, discrete-frequency checks. They do not establish
broadband/hereditary accuracy, coupling to omitted modes, nonlinear mechanics,
or a thermally coupled constitutive history.

The static opening and pressure comparisons still use the source EQUILIBRIUM
stiffness, not K'(omega_ref). Thus even a complete retained frozen model must
disclose any lost equilibrium compliance. Conversion to the legacy one-mode
playing image must additionally pass its existing playing-mode admission. Every law must admit DC, the reference frequency and each modal
sample. The wrapper retains the original field and approximation report;
identities bind constitutive tensors, branch times, density, source identities,
field ordering and the selected reference frequency.

Authored loss must be zero on this source-material path. The existing card
codec preserves each mode's actual damping_n_s_m, and reed lowering reads that
coefficient as c/(2*sqrt(k*m)). The legacy loss_factor field denotes only the
authored contribution; it is not total material loss. Exact legacy eta/2
rounding is retained only when the authored formula reproduces that mode's c.
The same source loss consequently changes the existing mechanical runtime's
energy decay after serialization, rather than merely changing metadata.

The focused visco_valve integration target covers analytic single-coordinate
mass/stiffness/pressure normalization, source-viscosity changes to actual
port-Hamiltonian decay, Prony changes to frequency and damping, equilibrium
compliance disclosure, codec-to-reed continuity, uniform/per-element identity,
and refusal of incompatible bands, double-counted loss, excessive dispersion
or modal loss requiring a coupled runtime.

### Complete selected material valve modes

`reduce::visco::reduce_modal_pressure_loaded_valve` explicitly exports a
`MaterialModalValve` for a coupled runtime. It reuses the same source assembly,
modal window, frozen-frequency dispersion checks and equilibrium comparators,
retains the original `ViscoElements`, every selected mass-normalized mode and
the complete projected viscosity, and permits the certifier's full blocks.
A mode with zero mean opening remains available to carry pressure work,
external force or coupled material loss. The legacy card and one-mode path
retain their refusal; no independent-mode conversion is exposed by this owner.

Signed mean-opening and pressure projections stay separate. Nodal force
projection and nodal displacement/velocity expansion use the same original
shapes and support map. Static diagnostics include all selected modes while
preserving the original equilibrium reference. Identity binds the complete
basis, full viscosity and signed ports. This is still an admitted frozen K/C
image, with no claim to broadband hereditary, omitted-mode, contact or
large-strain accuracy.

The added `visco_valve` regressions exercise an exactly zero-opening,
pressure-carrying mode in an anisotropic two-coordinate tetrahedron, full-C
work against the original element assembly, virtual-work reciprocity, legacy
refusal, source-dispersion admission and cancellation.

### Hereditary complete material valve modes

`reduce::visco::reduce_hereditary_pressure_loaded_valve` exports a separate
`HereditaryModalValve` from the original equilibrium pencil. It retains all
selected mass-normalized modes, including zero-opening modes, their signed
pressure/force observations, original tetrahedral laws and complete Kelvin C.
The same `ViscoElements::project_basis` supplies every quadrature strain map,
physical volume and tensor-Prony factor. This path does not freeze storage or
loss at a reference frequency, use modal loss ratios, or loosen the legacy
frozen/independent-mode admission.

Physical coordinates satisfy q=sqrt(total_mass)*x. Explicit projection and
factor-work caps apply before workspaces are allocated. The projected mass
and equilibrium stiffness must match the original modal pencil within the
caller-declared basis error. DC and the conservative complete instantaneous
frequency bound must lie in the original material bands.

Bit-identical relaxation times share one existing Householder-QR row reduction.
Column equilibration and original Gram reconstruction preserve all coupled and
weak physical directions without a relative rank cutoff. At most n channels
remain per distinct source time, independent of the number of mesh points.
This is energy/work/heat equivalent for wholly relaxed or unrelaxed initial
histories driven by the same retained coordinates; arbitrary independently
authored quadrature memories and per-element heat reconstruction are not
provided. Original material strain limits and maps remain attached and bound
into image identity. A conservative Mandel-norm bound admits ordinary runtime endpoints in O(n);
near a limit the original quadrature maps decide. Endpoint strain admission
is not an interval enclosure.
Source/modal/spatial convergence and sampled coupling accuracy remain separate.

Regressions compare full frequency-dependent matrices with the original
sparse material assembly, and compact-channel stress, energy, work and heat
with the original tensor history under signed shear ramps and holds. They
also cover strong dispersion refused by a frozen fit, dark modes, very weak
loss directions, source bands and explicit capacity refusals.

## Source tetrahedra projected into full-memory structural dynamics

`linear3::visco::ViscoElements::project_basis` projects an explicit displacement
basis using the original Tet4/Tet10 strain matrices, physical quadrature volumes,
material laws and constrained free-DOF map. The resulting
`ProjectedViscoSolid` retains the complete consistent mass matrix
`Phi^T M Phi`, each original quadrature-point tensor constitutive contribution,
and its owning element. It converts Mandel shear to tensor shear at the existing
`TensorProny::project_strain` boundary. All Kelvin cross terms and all original
Maxwell factor rows survive; no mass lumping, material averaging, fitted damping
or silent basis truncation is introduced.

The same retained basis maps nodal loads by its transpose and expands
displacement/velocity with the original supports. These maps preserve discrete
virtual work; fixed rows contribute no work. `SolidProjectionBudget` bounds
coordinate count, stored basis entries, quadrature points and scalar projection
work before constructing the result. Overflow, nonfinite maps, source refusal
and cancellation return errors without a partial projection. The immutable
projection borrows its original source assembly.

A dynamic consumer must independently admit positive-definite, representable
projected inertia. Geometry, supports, basis and material state remain fixed.
This API does not claim spatial/modal convergence, arbitrary frequency validity,
nonlinear strain or thermally evolving geometry. Focused G0/G3 integration
regressions in `fs-couple/tests/tensor_solid.rs` cover Tet4/Tet10 shear and
consistent mass, nodal virtual work, nonorthogonal constrained bases, original
Maxwell prehistory and regional material-loss swaps. They use synthetic material
fixtures and do not establish empirical calibration.

## Continuing thermal-valve motion across source updates

`reduce::thermal::motion::ThermalValveModalMotion` owns the actual thermal
reduction and the existing `fs-modal` time-domain state. The reducer retains
its original mass-normalized mode shapes and exact pressure-load projections.
Pressure steps use the original exact ZOH owner; applied work is conjugate to
swept volume, and nodal observations reconstruct equilibrium plus vibration.

Candidate source updates preserve the exact sample clock and use the shared
mass-inner-product transfer on the same reference mesh, frame, supports and
consistent mass. New equilibrium, stiffness, mode count/order and supplied
pressure observations are admitted with explicit state-loss and signed
incremental-storage budgets. Updates and pressure steps publish only complete
candidates. Discarded state is not damping or heat; bore feedback, radiation
reaction and reciprocal thermoelastic heat are not inferred. Three focused
one-tetrahedron tests use independent consistent mass, stiffness, thermal
equilibrium and oscillator references to check motion/work, source-update
continuity, changed pressure response and refused/cancelled candidates.


### Native finite-strain storage and reference-mass access

The `hyper3` owner exposes full referential dead loads, actual stored energy and
a consistent reference mass for dynamic consumers. Energy uses the same
deformation gradient, determinant/stretch admission and constitutive owner as
the native Piola force and tangent. Reference mass delegates to
`linear3::p1_tet_consistent_mass`, with one positive finite density per element;
unrepresentable element/assembled mass refuses. Consumers must retain that
mass during deformation instead of multiplying updated volume by a fixed
reference density. These accessors introduce no constitutive or time solver.

`FinitePlasticTetMotion::try_new_with_nodal_density` also admits an authored
positive P1 reference density on the original shared vertices. Its consistent
matrix integrates rho*N_a*N_b analytically through degree three and retains
the complete nodal inventory across the existing finite-plastic momentum
steps. The cell-density accessor reports cell means; the separate nodal
accessor preserves the actual distribution. `reference_acceleration_load_n`
uses this SAME mass action to project held nodal accelerations, including
uniform gravity, into the existing nodal-load input. No current-volume mass
refresh or one-quarter-per-node lumping is introduced. Uniform nodal density
delegates to the original cellwise mass path. Independent simplex moments,
rigid kinetic energy and native backward-Euler work checks cover the mass
operator; the finite-plastic regression exercises the native free-body Newton
step. This is not material-source inference or a new nonlinear-gravity
quadrature, deformation-field solver, or thermal mass-transfer model.


### Impulses with complete finite constitutive memory

`FinitePlasticTetMotion::impulse_candidate` and
`FiniteViscoTetMotion::impulse_candidate` share the original `ReferenceMass`
owner for zero-duration momentum transfer. Native checked CG solves
`M_ff delta_v_f = J_f` through borrowed sparse rows with diagonal
preconditioning; it builds neither a dense inverse nor a second mass matrix.
Both cellwise and graded nodal reference density retain their exact existing
mass. Caller-supplied held DOFs keep their current velocity, including nonzero
support motion. Repeated applied DOFs accumulate in input order.

`SolidImpulseReport` resolves linear and angular momentum at the actual
deformed positions, applied and support impulses, kinetic change, and signed
applied/support work. Physical gates use the actual representable velocity
change and include each free-DOF residual, so cancellation of net impulses
cannot conceal lost local momentum. Kinetic change is evaluated through
`average(v) dot M delta_v`; the work residual uses the independently
recomputed impulse defect instead of subtracting large work totals.

Only velocity changes. Original deformation, plastic metrics and accumulated
flow, every Maxwell branch, stress, storage, material time, law bindings and
reference inventory survive the event. Candidates publish only after solve,
workspace, momentum, work and cancellation gates. Subsequent native motion
continues the accepted material history. Impulses do not themselves generate
relaxation/plastic heat or select a contact/restitution law. The caller owns
the opposite impulse and any contact energy loss; mixed incompressibility
constraints are outside these displacement-only P1 entry points.
The `hyper_plastic_dynamics/impulse.rs` and `hyper_visco_dynamics/impulse.rs`
regressions use independent cubic-simplex mass, deformed torque and support-work
oracles, retained inelastic prehistory and subsequent native motion.

`FiniteThermoplasticTetSolid::impulse_candidate` extends this same account to
an already inertia-bound thermal tetrahedral solid; a quasistatic owner without
original reference inertia refuses. `FiniteThermoplasticHexMotion` supplies
the same event for both standard and mean-dilatation Hex8 formulations, using
the original tensor-product quadrature mass and physical current positions.
The thermal family bindings and reference temperatures, current temperature,
all plastic/creep and quadrature history, physical/assumed deformation, stress,
free energy and internal energy remain unchanged. The event advances no
material interval and deposits no heat. `tests/thermoplastic_impulse.rs`
checks hot plastically set bodies, independent reference-mass/torque/work
oracles, moving supports, state retention and subsequent native thermal motion.

### Material-point mobility with inelastic memory

The finite plastic and viscoelastic motion owners expose `TetImpulsePoint`,
`point_kinematics`, `point_directional_inverse_mass_per_kg` and
`point_impulse_candidate`. One admitted P1 barycentric trace B determines
current position and velocity; its exact transpose distributes a world
impulse. Nonnegative unit-sum weights admit interior, face, edge and vertex
points, with only roundoff-sized sum normalization. Extrapolation refuses.
The directional mobility is `n^T B M_ff^-1 B^T n` for a unit world direction,
using the original cellwise or graded reference density and held velocity DOFs.
It borrows the existing checked-CG mass response without cloning constitutive
history or forming a dense inverse. A fully held trace has zero mobility.
Only the accepted impulse candidate copies the complete unchanged material
state. The native visco impulse regression checks an interior point against
independent cubic-simplex mass inversion, including negative consistent-mass
nodal responses, supported mobility, deformed torque and conjugate point work.

The same owners expose `point_mobility_matrix_per_kg` and
`point_impulses_candidate` for simultaneous contact. Every row/column pair of
`B M_ff^-1 B^T` retains its full coupling through the original supported mass.
One checked response is retained at a time; no dense mass inverse or material
history copy is needed to build mobility. Caller row and workspace limits
bound the matrix, trace maps and solve reports before allocation. The returned
matrix is not silently symmetrized or clamped; the contact-law owner admits it
with its own numerical tolerances. All point loads reach the compensated nodal
accumulator together and produce one native velocity, support and work account.

Hot Tet4 and Hex8 motion expose the same five point operations. Tet4 uses the
original barycentric trace; `HexImpulsePoint` supplies reference coordinates
in `[-1,1]^3` for the original trilinear trace. Physical
deformed coordinates determine position, velocity and torque even when the
Hex8 constitutive owner uses mean dilatation. Original density, thermal source
families, temperatures, quadrature histories and free/internal storage survive
unchanged. Tet4 requires its reference inertia to be bound first. Focused native
regressions compare complete two-point matrices with independent simplex and
tensor-product inverses, including oblique directions, nonvertex points,
supported responses, summed work and unchanged thermal/material state.

### Exact-incompressible pressure impulses

`mixed::conservative::ConservativeTetSolid::impulse_candidate` solves the
original P2/P1 velocity and pressure-impulse KKT equations at the current
configuration. The endpoint cofactor supplies the same weak volume derivative
B as the native mixed geometry. Original consistent mass supplies M; native
sparse LDL handles `[M_ff B_f^T; B_f 0]`. Explicit velocity, pressure-impulse
and work scales preserve a symmetric dimensionless system, with assembly and
factor budgets and a recomputed true residual. No second mass solve occurs.

The multiplier is an impulse in Pa s, separate from the stored algorithmic
pressure. Held velocities remain unchanged; full support impulses and work
remain explicit. The report separates external, internal pressure and support
contributions. Actual represented velocities must satisfy every original weak
volume-rate row, including a removed gauge row. Pressure contributes neither
net linear/angular momentum nor admitted work beyond explicit tolerances.
A gauge is allowed only for a verified connected constant-pressure null mode;
other rank deficiencies refuse. Incoming incompatible velocity refuses rather
than silently receiving a dissipative projection. In particular, an earlier
volume-constrained time step does not automatically guarantee endpoint tangent
velocity admissibility. The impulse preserves geometry, density, material
points, stored pressure and time, and adds external/support work to the
existing conservative trajectory ledger. `tests/mixed_conservative.rs` checks
independent P2 mass moments, pressure loading, rotated geometry, supported
gauge reactions and continuation. The exact-incompressible branch uses this
pressure solve; finite-compliance conservative formulations retain their
original mass-only response and supply no pressure impulse.

### Reference-aware finite-deformation response

`hyper3` evaluates energy, Piola force and tangent through the original
`HyperelasticMaterial` interface, preserving source reference-deformation
adapters. Bare laws retain the default native response. The actual total
deformation must have a finite positive Jacobian independently of the
material's elastic-gradient admission. Reference geometry and the consistent
reference mass remain fixed; a thermal adapter does not alter the mass mesh.

### Source-defined cohesive load paths

`cohesive::CohesivePatch::from_interface` consumes
`fs-material::cohesive::source::ResolvedBilinearCohesive`, checks the actual
ordered surface material assignments, and retains the original source bundle
through every accepted displacement/damage state. `p1::MaterialFacePair` and
`P1Interfaces::from_interfaces` lower the same source through the existing
matching-face quadrature and work-conjugate maps; both full and reduced patches
retain their source. No constitutive, quadrature or nonlinear solver is copied.

Swapping an interface dataset changes the traction/separation envelope, failure
point and area-integrated fracture work in the native equilibrium solver.
Conditions remain fixed during each mechanical increment. The explicit
`cohesive::temperature` transfer below advances the same sourced interface to
new temperatures without resetting history. Chemical evolution, adaptive
fracture topology and complete structural/thermal energy accounting are not
introduced. The `cohesive_source` integration tests
check an independently solved compliant pull, material substitution, exact P1
area/fracture accounting, surface-order refusal and atomic cancellation.

### Nonmatching bonded P1 surfaces

`cohesive::p1::nonmatching` binds independently triangulated planar surfaces
through their triangle intersections. Each intersection is triangulated and
integrated with positive degree-two quadrature. Distinct barycentric weights
on the two sides build one shared jump map; the existing cohesive operator uses
its transpose for forces and its consistent tangent for Newton iterations.
Undamaged linear traction is integrated exactly to floating-point precision,
and nonlinear damage retains discrete virtual work and force balance. Moment
and rigid-motion consistency hold on the common plane; any admitted off-plane
projection contributes error bounded by the supplied geometric tolerance.
The source-bound constructor retains the original ordered interface source at
every full and reduced quadrature point.

Construction checks exterior-face incidence, disjoint side DOFs, opposite-side
bodies, face quality, planarity, within-side overlap and full area coverage on
both sides. Work and storage limits bound face count, pair checks and quadrature
points, with cancellation throughout. A deterministic bounding-box sweep skips
disjoint pairs before they consume the pair budget. Worst-case preprocessing
remains quadratic for heavily overlapping boxes. Runtime
solves reuse the fixed sparse maps and original constitutive/solver path.

This is an explicit initially bonded, fixed small-strain reference interface.
The caller admits the relative common-plane projection/coverage tolerance.
Curved interfaces require separate locally planar regions; large sliding,
dynamic rebonding and adaptive fracture topology are not inferred. Nonlinear
softening still requires independent surface and load-path refinement. Tests
cover crossed-diagonal area and analytic linear energy, virtual work and rigid
motion, matching-mesh stiffness recovery, force-controlled softening/fracture/compressive
reclosure, source retention, malformed coverage and construction budgets.


### Exact native transient continuation across material lifetimes

`loads3::transient::checkpoint::TetTransientCheckpoint3` detaches accepted native
motion from its borrowed material problem. It retains original geometry, shared
assembled operators, motion/acceleration, thermal field, applied loads, support
state, controls and accepted step count. Binding requires the same mesh, fixed
components and complete material values/source identities. It neither
reassembles nor restarts dynamics; changed coefficients still use the existing
thermal/rebind owners. This immutable in-process continuation is not portable
serialization. G3 continuation, step-budget preservation and changed-material
refusal are exercised by the reciprocal dynamic thermoelastic tests in
`fs-couple`.

## Consistent reference inertia for planar frame vibration

`frame2d::FrameMassFormulation::ConsistentReference`, selected by
`Frame2d::with_mass_formulation`, integrates linear axial and cubic Hermite
transverse centerline motion on each original Euler–Bernoulli member. The
result retains translation/rotation and endpoint cross terms; separately
added point masses remain diagonal. `Frame2d::new` retains the original
lumped default, and a model with only nodal mass retains the diagonal path
even when consistent element inertia is selected.

Cross-section rotary inertia uses the same owner: material members derive
`sum(rho I)` from each original layer's density and centroidal area moment,
including the layer's own material state. `add_element_with_rotary_inertia`
accepts an explicit section mass moment per length [kg m] for raw-coefficient
members; `add_element` retains its declared zero-section-spin approximation.
Consistent mass integrates `rho I (d(v_dot)/dx)^2` exactly with the Hermite
slope (Rayleigh beam kinematics), including its endpoint and translation-spin
cross terms. Endpoint lumping assigns half the section spin inventory to each
node. Neither path counts centerline orbital inertia a second time.

One reference mass owner supplies the generalized vibration pencil, the
initial acceleration solve, every Newmark inertial residual/tangent, Rayleigh
mass damping, uniform-base excitation, kinetic energy and support reactions.
Massless free coordinates still use the original static stiffness condensation;
the full massive block is reduced with the existing fs-la Cholesky owner before
the existing symmetric eigenvalue solve. Mass is prepared once per analysis,
not rebuilt within Newton iterations. `inertial_forces` applies the same M to
an explicit full acceleration vector, including consistent gravity/body loads
and their generalized nodal moments. Reaction records include supported inertia,
damping and applied loads, including the initial endpoint.

This is fixed reference-configuration inertia for linear or moderate-rotation
vibration. It does not include configuration-dependent mass, centrifugal or
gyroscopic terms, or shear deformation. Distributed body loads are
projected through the reference interpolation; one force-based element does
not certify spatial convergence of dynamic stress or distributed-load section
recovery. The existing dense building-scale solver and constitutive/kinematic
limits remain. Selecting a mass model does not alter the material law.

`tests/frame2d_mass.rs` checks exact rigid translation/rotation inertia in
rotated frames, two cantilever frequencies under mesh refinement, free-body
base acceleration and work, full supported inertia reactions, a loaded
two-coordinate history against independently solved generalized eigenvalues,
and unchanged nodal-only behavior. The one-element first bending-frequency
error falls from about 30% for endpoint lumping to below 0.5% with consistent
mass; higher-mode accuracy still requires member refinement.
Additional G0/G1 tests compare section spin against independent slope
quadrature, exact rigid-body rotational inertia and the analytic simply
supported Rayleigh spectrum under mesh refinement. Source-layer tests check
the full density-weighted spin inventory and its effect on the original
Newmark trajectory.

## Signed, mass-normalized frame modes and participation

`Frame2d::modal_analysis` exposes the complete signed current-tangent spectrum,
full node-major mode shapes and horizontal/vertical uniform-base participation.
Shapes satisfy `phi^T M phi = 1`; fixed coordinates are zero and massless
coordinates use the original static condensation map
`phi_r = -K_rr^-1 K_rt phi_t`. The massless stiffness must be positive definite:
an unstable or unconstrained massless subsystem cannot silently disappear
through condensation. Each massless and mass factor is prepared once.
The original fs-la Cholesky and symmetric Jacobi owners supply the algebra;
the former frame-only eigenvalue kernel is retired.

With `r_x/r_y` the full unit translational influence, including supported
coordinates, participation is `Gamma = phi^T M r` and effective mass is
`Gamma^2`. The same reference M is used by the frame transient and body-force
owners. Supported inertia may also contribute directly to reactions, so the
modal sum need not equal the complete `r^T M r` inventory.

Every mode retains its unmodified signed eigenvalue, a freshly recomputed
free-pencil backward residual and a mass-whitened residual. The reported
sign-resolution diagnostic is the latter plus
`64*eps*n*||L^-1 K_cc L^-T||_inf`. It is a floating-point estimate, not a
certified interval or stability theorem. Modes are explicitly oscillatory,
unstable, or unresolved near zero; unresolved does not certify a neutral mode.
Resolved negative eigenvalues expose their linearized exponential growth rate.
`natural_frequencies` now requires every mode to be resolved and oscillatory;
rigid or unstable frames must use the signed analysis instead of receiving
clipped zero frequencies.

The analysis is frozen at the current admitted material/displacement state;
it does not locate a nonlinear buckling load or refresh during a history.
Repeated eigenvalues define a mass-orthonormal subspace, not unique vectors.
`tests/frame2d_modes.rs` covers exact massless beam-shape restoration and
participation, full-mass orthogonality and reconstruction, free rigid motion,
and an axially compressed column crossing its Euler instability.

### Material-defined frame members

`Frame2d::add_material_element` consumes the existing
`fs-material::section::ConcentricElasticSection`. It derives axial rigidity,
bending rigidity and reference line mass together from the original layer
integrals and retains the full section through `material_section`. There are
no independent stiffness/mass overrides on this entry point. Heterogeneous
member assignments and perfectly bonded concentric layers retain their own
material cards, query states and receipts; no average material replaces them.

The original force-based static, signed-modal, reference-mass, gravity and
Newmark owners consume those coefficients. A sampled frequency/rate law is
refused, as is a finite frequency-applicability band that this general time
history cannot enforce. Constant or temperature-only coefficients are frozen
at the declared state. This constructs a new stress-free reference specimen;
it does not implement thermal expansion, live mass exchange, plasticity,
interlayer slip or material-derived damping. Existing planar Euler–Bernoulli
and reference-inertia assumptions still apply. Missing data and source/geometry
admission failures publish no partial member or added mass.

`tests/frame2d_material.rs` exercises independent cantilever compliance,
`sqrt(E/rho)` frequency scaling, integrated weight, layered and spatial
assignments, and a guided-tip oscillator's actual Newmark history and energy
balance. `cargo run -p fs-solid --example material_frame -- 1000 2e9 4000 2e9`
demonstrates data-only substitution through the same geometry and solvers.
The example inputs and regression data are synthetic, not measured alloy
qualification.

### Physical bond rates in native cohesive equilibrium

`CohesiveState::solve_over_time[_controlled]` uses the same accepted history,
Newton/FGMRES owner, work-conjugate jump/scatter map, force gate and atomic
publication as the original equilibrium path. Explicit dt gives the physical
constant jump rate; source-selected bond viscosity contributes to total endpoint
traction and its complete consistent tangent. Untimed solves refuse viscous
bonds. Intrinsic cohesive work/storage/fracture and physical viscous loss are
reported separately; no structural stability or unique-root claim is added.

The existing restrained thermal joint passes its actual interval to both bulk
Maxwell and cohesive owners. Its total bond path work includes analytically
integrated viscous loss; endpoint-force quadrature corrections and the combined
energy audit use that same work. Bulk loss, intrinsic fracture and bond viscosity
stay distinct. Loss is not silently deposited as heat, and source coefficients
are not refreshed across a damaged bond history.

### Source-bound frame strain-rate dissipation

`FrameViscousSection::try_new` binds explicit axial and bending viscosity
properties to the SAME original layered elastic section.
`Frame2d::add_material_viscous_element` retains those sources, stiffness and
reference inertia together. Every native Newmark call includes the complete
physical damping matrix. Heterogeneous member losses retain directional,
rotational and inter-node terms; no common Rayleigh ratio is inferred.

`FrameMaterialDamping` retains the original sections and three positive
work-conjugate resistance ports per member. The native full matrix is assembled
from those exact ports, which are also available to reduced receivers. Each
port maps physical nodal displacement to an extension in metres; rotational
weights therefore carry metres. `material_damping_energy` and
`rayleigh_damping_energy` are separate parts of the existing total damping work,
not additional heat sources.

This path requires fixed linear geometry and elasticity throughout the frame,
positive inertia on every free coordinate, explicit viscosity support from DC
through the actual sampled Nyquist band, and the complete original spectrum
below Nyquist even when a receiver retains fewer modes. Temperature and geometry
remain frozen. It adds no hereditary loss, thermal deposition or acoustic
radiation. `tests/frame2d_viscosity.rs` supplies continuous oscillator and
heterogeneous two-coordinate references, physical-port work equivalence,
continuation, and source/admission regressions.

## Complete frame material inventory for live thermoelastic consumers

`FrameMaterialDamping::covers_all_members` records whether every original frame
member has its material-bound elastic and viscous section, including members
fixed at both ends. A complete material thermal consumer must account for those
members' stored energy and physical layer inventories even when they do not
appear in a free modal coordinate. This flag adds no mass, stiffness or damping
law. The original operator still rejects massless free-coordinate condensation
and retains its complete reference geometry, source sections and native ports.

## Forced resonance on the native linear frame

`Frame2d::prepare_harmonic_response` connects the original reference mass,
elastic stiffness, full section viscosity and explicitly authored Rayleigh
matrix to `fs_modal::harmonic::PreparedHarmonicResponse`. The prepared factor
is reusable for force cases at one frequency. Its `free_dofs()` map preserves
node-major translations/rotations and their conjugate forces/moments; no mode
or off-diagonal inertia/loss term is discarded. The source frame is borrowed
immutably, including its accepted displacement and retained Rayleigh reference.

Peak phasors use `exp(-i omega t)`. Responses contain physical amplitude,
phase, velocity, numerical backward residual and independent cycle-averaged
input/viscous powers. Source viscosity and Rayleigh damping appear once in the
total loss; this is not automatic thermal deposition. Original source bands
admit the actual queried frequency, including explicit DC when supported.
Singular lossless poles refuse without a damping floor or artificial mass.

The scope is fixed linear geometry, elastic sections, stationary homogeneous
supports and positive inertia on every free coordinate. Consistent member mass
or source-derived section spin supplies physical rotational inertia. Original
rigid modes remain intact at positive forcing frequency. Source viscosity
assembly is shared with the native transient owner, without repeating its
natural-mode eigensolve or inventing a Nyquist clock for a harmonic query. The
public transient operator retains its complete-spectrum/timestep admission.
Massless viscous coordinates, moving supports, support reactions,
nonlinear periodic orbits and startup transients are not supplied by this API.
Roundoff-sized native stiffness asymmetry is admitted in diagonal-scaled
coordinates and averaged in scratch; larger skew refuses. Source coefficients
and accepted frame state are not rewritten.

Controls bound the full coordinate count and conservative numerical payload
before dense assembly. Cancellation polls surround native bounded preparations,
per-member stiffness assembly and the shared factorization/solve; an interrupted
native mass/material call finishes before its next checkpoint. G1/G3 tests in
`tests/frame2d_harmonic.rs` and `tests/frame2d_viscosity.rs` cover analytic axial
resonance, full consistent bending inertia, heterogeneous nonproportional loss,
rotational covariance, reciprocity, source-band boundaries, physical power,
singular limits, budgets, cancellation and unchanged accepted state.

## Original material strain domains on frozen modal solid trajectories

`ViscoElements::validate_displacement` evaluates the existing element Mandel
strain maps at every original material quadrature point and invokes the source
law's strain admission. `MaterialModalValve::validate_modal_displacement` uses
the existing modal-to-node expansion before that same check. These observations
advance no material history and construct no alternate elastic operator. They
poll cancellation at source points and retain all accepted state on refusal.
TET10 observes its admitted quadrature; unsampled spatial extrema and maxima
between time samples remain separate accuracy obligations.

## Original exterior geometry on material valve images

Frozen and hereditary `MaterialModalValve` images retain one immutable
`MaterialValveSurface` prepared from the original Tet4 nodes and exterior
`TetLoadGeometry3` faces. Canonical node triples, physical face areas and
outward normals use the original solid owner, including reversed tetrahedron
winding. Preparation respects its existing cancellation and load budgets.
No geometry is rebuilt during stepping and no material response is fitted.
Tet10 images deliberately carry no P1 exterior surface: a consumer requiring
that trace must refuse rather than invent midside geometry or interpolate an
unqualified contact surface. Mechanical pressure and swept-volume maps retain
their original ownership independently of a consumer's contact-face selection.

## Continuing material-fluid rod loads

`PreparedRodMotion::body_loads` exposes the external conservative wrench at an admitted dynamic pose using the same gravity, hydrostatic displacement, attached-strand forces, and reference potential as the native force assembly. Nodal forces are expressed in world axes; pure attachment couples are expressed in the current body axes. Rotating attachment arms retain their actual mass and displaced-volume distributions, including cases where their centroids differ. `rebind_submerged_fluid` replaces an existing uniform fluid snapshot through `FluidTransport` while retaining the original solid reference inventory, inertial data, displaced-volume policy, reference geometry, and stationary supports. This model-level operation neither advances time nor accounts for environmental work; the continuing time owner must account for the potential change at its actual pose. These operations retain the existing full-immersion, nonoverlapping-volume, fixed-reference-volume approximation and infer no drag, added mass, partial immersion, free surface, solid compression, or reservoir dynamics.

## Continuing source-temperature rod sections

`PreparedRodMotion::rebind_material_temperatures` refreshes section response through the original `RodMaterialSource`, retaining reference inventory, intrinsic geometry, supports, attached links, gravity and hydrostatic displacement. It admits the exact reference material description and verifies replay of current section claims, expansion paths and strain bounds before resolving a new temperature field. Original unheated reference sections may acquire their explicit thermal binding. The source owner re-evaluates E/rho/G and integrates expansion from the original stress-free reference; current density does not redefine closed-body mass or rotary inertia. `strain_energy_j` evaluates the original conservative strain-storage owner at an admitted dynamic pose. Section geometry, reference displaced volume and separately bound tension-link material states remain fixed.

## Native rod strains and source thermoelastic lowering

`RodMaterialSource::thermoelastic_laws` lowers each original concentric segment/layer through the same pinned mechanical and thermal-expansion data into the complete material thermoelastic law. `PreparedRodMotion::segment_strain_states` returns the actual-pose axial, bending, shear and torsional strains and section energy through the existing native kinematic and energy owners. Neither path introduces a replacement rod model, effective material, or extra trajectory.

The reciprocal finite-enthalpy consumer in `fs-couple::rod_motion::enthalpy` uses these views to evaluate thermal exchange and refresh the original source-derived stiffness and integrated expansion while retaining reference mass/inertia, pose, momentum, supports, attached links, loads and clock. Existing source and small-strain refusals remain authoritative. Independent native-energy, nonuniform conduction and torsional-reaction regressions were added under `fs-couple/tests/rod_motion/enthalpy.rs`; native execution remains pending the unchanged dependency/environment build blocker.

## Original bulk and interface sources on thermal joints

`cohesive::thermal_joint::source::SourceJoint` binds the existing restrained-member
and cohesive solver to original `ResolvedThermoMaxwell` bulk sources and one ordered
`ResolvedBilinearCohesive` interface. Receiving member identities must match both
interface sides. Source temperature, strain and rate domains remain enforced by
the native material law; member memory, cohesive damage, force/extension
compatibility and work/storage checks remain owned by the existing joint.

The ordinary advance explicitly admits the unchanged interface query, including quantity
conventions, so evolving bulk temperature cannot silently freeze an unsupported
cohesive-state change. Three source integration tests cover independent series and
thermal-work references, data-only substitutions with relaxation and fracture
continuation, and source/domain/cancellation refusal. This remains a quasistatic
fixed-geometry series model with prescribed bulk temperatures and fixed external
interface state; no thermal transport, fracture-heat partition or
temperature-dependent cohesive-history transfer is inferred from that path.

`SourceJoint::advance_with_interface_temperature[_controlled]` explicitly
accepts a temperature successor from the same card, envelope and pinned claims.
Its normal series geometry consumes common-mode or mixed-mode pure-opening
data; the material extension does not supply shear equilibrium in this owner.
It stages the material owner's fixed-jump history transfer before the original
native equilibrium solve, retaining every bulk Maxwell strain, reference
geometry and source assignment. `TemperatureJointReport` separates signed
interface parameter work, fixed-jump damage release and the native mechanical
accounts. An independent old-to-new energy balance includes both stages before
joint state, source receipts, damage and physical time publish together.

Three focused source/native regressions add an analytical series-force and
parameter-work reference, hot strength/toughness loss followed by cooling
without healing, and source-domain/late-cancellation refusal with identical
retry. This is first-order prescribed-temperature composition; subdivide
combined thermal/load paths for accuracy. It does not solve thermal transport,
assign fracture heat, validate an unobserved continuous temperature path, or
refresh frozen-onset mixed-mode history.

The 34 cohesive constitutive tests, including the four temperature-transfer
tests, passed in an exact-physics slice. Source/native integration targets remain
unexecuted in this worker: compiling the unchanged pinned `asupersync` dependency
was SIGKILLed even with single-threaded metadata-only compilation.

## Original thermal eigenstrain in tensor-Prony projection

The existing tetrahedral tensor-Prony projector carries each explicitly supplied
thermal free strain into its original material quadrature instead of discarding
it or inventing an external force. Original geometry, physical strain columns,
consistent inertia, material frame and supports remain unchanged. Eigenstrain
admission uses the existing small-strain material domain. Its constitutive
energy and continuing Maxwell memory belong to the original projected material
and native coupled PHS owners; projection adds no heat or independent solver.
Source-integrated isotropic/oriented-orthotropic expansion and reciprocal heat
are exercised by the tensor coupling regressions, currently authored but unrun
in this environment.

## Shared J2 plasticity in fibers, beam-columns and frames

`FiberLaw::J2` consumes the original material owner's exact uniaxial stress
reduction, including plastic lateral contraction and separate isotropic/Prager
history. The existing section tangent, force-based beam-column and planar frame
solvers therefore support metal yielding, elastic unloading and permanent
curvature/deflection without a second structural plasticity law. Raw section
kernels retain their caller-admission convention; `Section::try_respond` and the
beam/frame nonlinear paths check finite response and the full compatible tensor
strain against the original J2 small-strain domain.

`FiberMaterial::resolve_j2` uses the same source requirements, explicit hardening
basis and coefficient conversion as tetrahedral solids, retaining density and
all original receipts at one fixed external state. RC reinforcement may select
the existing steel law or native J2; its replaced concrete area and reference
mass remain coherent. A source post-yield tangent does not invent a kinematic
hardening split. Temperature feedback, rate effects and calibrated cyclic steel
behavior are not inferred for this fixed-state axial reduction.

A force-based cantilever commits every section only after all sections converge;
a refused step preserves the whole accepted element for retry or reversal.
Focused regressions compare independent bilinear section and beam formulas,
source/mass consistency, full-tensor domain refusal and transactional history.
These additions do not provide finite-strain forming, plastic collapse
continuation, fatigue/damage or a complete nonlinear plastic shell solver.

Complete unloading in the native frame iteration retains its independent
compatibility and relative force checks while admitting a local constitutive
roundoff floor derived from section tangent and strain. This avoids rejecting
a physically equilibrated permanent bend solely because endpoint load is zero.

## Saturating metal fibers

`FiberLaw::Voce` and its retained axial history use the exact native
`UniaxialVoce` restriction in existing sections, force-based beam-columns, and
frame elements. `FiberMaterial::resolve_voce` admits independent saturation
magnitude/rate, linear isotropic/Prager coefficients and reference density in
one source bundle. Reinforced sections preserve the same concrete replacement
and constituent mass calculation. Checked section response admits the full
tensor strain including lateral contraction before a structural commit.

`fiber_plastic` checks the independent nonlinear moment/curvature relation,
consistent section tangent, source retention, unloading curvature, and matching
beam/frame permanent deflection. These fixed-condition, small-strain material
checks do not establish cyclic calibration, thermal recoupling, damage, or
finite-strain forming. Existing geometric choices keep their own model limits.

## Incremental plane continua

`plastic2::SmallStrainPlaneProblem` uses the existing `Mesh2` P1 triangles and
Q1 quadrilaterals with a native uniform or per-element `SmallStrainLaw`.
Plane stress solves all transverse stress constraints with the original
material reduction and its condensed algorithmic tangent. Plane strain
retains the constrained transverse stress. Every quadrature point carries its
complete three-dimensional strain and committed constitutive history.

Absolute chronological load factors support loading, unloading, reversal and
repeated physical intervals for timed laws. All Newton and line-search trials
reuse the incoming history. Only a converged increment returns updated history
and its matching endpoint stresses. `solve_from` leaves its borrowed accepted
state untouched on refusal. Reactions and free-row force convergence retain
their original mechanical meaning. Geometry, quadrature, matrix-entry and
load-path budgets are explicit, as are full tensor strain and small-rotation
bounds. Prepared geometry and B columns, reused C B products, and one accumulated
element matrix avoid repeated geometry work and duplicate quadrature insertion.
The existing tetrahedral sparse linear correction is shared without changing
its numerical acceptance criteria.

G0/G1/G3 `plastic2` tests cover P1/Q1 cyclic analytical references, confined
plane strain, skew-element global tangents, heterogeneous yielding, timed
relaxation, atomic continuation/refusal, geometry and rotation admission, and
Voce saturation. Standard displacement elements can lock in confined nearly
incompressible plastic flow; the opt-in Q1 plane-strain B-bar path below reduces
volumetric locking. Shell bending, finite-strain forming and limit-load
continuation remain outside this owner.

## Nonlinear backstress recovery in metal fibers and members

`FiberLaw::Chaboche` and its matching state use the exact free-sided axial
restriction of the native tensor law. Signed plastic strain, accumulated
equivalent plastic strain and independent recovering backstress remain separate
through section curvature reversals and member unloading. Section forces and
their tangent come from the same checked candidate; rejected domain exits leave
accepted fiber histories unchanged. Existing force-based beams and frame fiber
sections consume this law through their original section interface.

`FiberMaterial::resolve_chaboche` resolves density, Voce/Prager coefficients and
one independently sourced C/gamma branch from the same fixed-condition bundle.
All original property receipts remain available and density retains its original
mass meaning. The adapter does not infer nonlinear recovery from a monotonic
total-strain tangent or silently freeze evolving material coordinates.

The native fiber battery compares beam deflection and permanent set against an
independently stress-parametrized implicit-recovery curve, checks a reversed
section tangent and retained history, and exercises complete sourced properties
and mass. These are numerical constitutive/member checks, not calibration of
real cyclic behavior, finite-strain forming or a fatigue-life prediction.

## B-bar plastic flow and general constitutive Newton tangents

`SmallStrainPlaneProblem::formulation` admits `Standard` for existing P1/Q1
behavior and `BBar` for Q1 plane strain. The latter replaces only the full
three-dimensional volumetric strain by its element-quadrature mean, adding
`(mean(div u) - div u)/3` to each normal component. Local deviatoric strain,
plastic flow and separate quadrature histories survive unchanged. The same
projected B drives constitutive strain, internal force and algorithmic tangent;
the zz stress contributes to projected virtual work. Returned B-bar strains are
assumed material strains, whose zz component can be nonzero even though the
compatible displacement is plane strain. Both raw and projected strains, and
the raw rotation, must satisfy their independently declared admission limits.
B-bar plane stress and P1 are explicit refusals. History cannot be transferred
between formulations without a separately admitted projection.

Both native plane and tetrahedral increments preserve nonsymmetric constitutive
tangents. A law declaring `tangent_symmetric=false` selects existing FGMRES with
ILU(0) or absolute-diagonal preconditioning, bounded restart/iteration work,
in-Arnoldi cancellation checks and a recomputed true residual. The original
symmetric PCG/PMINRES path remains available. Nonlinear backstress recovery is
therefore usable under nonproportional load without false symmetrization; timed
tetrahedral dynamics inherits the same correction owner.

The focused plane battery includes constant plastic patches, skew-element
deviator/tangent checks, raw-domain refusal and a nearly incompressible pressure
cylinder against the independent plane-strain Lame displacement solution.
`tests/plastic_chaboche.rs` exercises actually nonsymmetric free tangents,
finite-difference force derivatives, triaxial/plane-stress/B-bar equilibrium,
retained nonproportional backstress and exhausted linear-budget refusal.
These tests do not establish a finite-strain, shear-locking or plastic-collapse
mesh-convergence claim.

## True-residual correction for confined continuum solves

The shared small-strain symmetric correction checks the actual Euclidean
residual independently of MINRES's preconditioner-norm recurrence estimate.
If the estimate stops early but the true residual misses the existing gate,
the solver may perform one fresh MINRES correction on `b - A x`. Its iteration
allowance is only the unused part of the original MINRES-path cap; no extra
allowance is created. Acceptance recomputes `b - A (x + correction)` relative
to the original right-hand side, retains only a better finite candidate, and
preserves the original linear gate and physical Newton force tolerance.
Cancellation or an exhausted allowance returns no accepted material history.

The nearly incompressible B-bar cylinder regression now also resolves a
12-by-48 mesh: a single bounded correction removes the observed mismatch
between the MINRES estimate and the independently checked residual. Its
displacement is compared with the same analytical Lame reference and the two
coarser meshes. A separate finite-load refusal test exercises an insufficient
linear allowance. This is a bounded residual correction, not a replacement
preconditioner or a claim that every ill-conditioned mesh can converge.



## Eight-node small-strain plasticity with full 3D B-bar

`plastic3::hex::SmallStrainHexProblem<L>` uses the original hexahedral geometry/Jacobian owner, FEEC two-point Gauss rule in each direction, and the shared native continuum increment/Newton solver. Uniform and per-element laws retain eight independent histories per cell. Absolute load paths support force/displacement control, unloading, reversal and continuation. The tetrahedral and inertial paths retain the same force scaling, corrections and old-history endpoint stress semantics after extraction into `plastic3::increment`.

Standard integration uses compatible strain. B-bar replaces only volumetric strain by the reference-volume-weighted element mean, with the full three-dimensional one-third correction. Identical projected B columns enter strain, virtual work and tangent; no material entries are symmetrized. Existing FGMRES handles genuinely nonsymmetric recovery tangents. Deviatoric quadrature remains fully integrated; there is no one-point stiffness or artificial hourglass penalty. Local 24-by-24 matrices accumulate on the stack before sparse insertion; no per-point matrix allocation is introduced.

Explicit node/element/point/free-DOF/matrix-entry/load-step budgets and sampled Jacobian checks precede assembly. Raw compatible strain and rotation are checked independently of projected strain. Histories commit only at converged increments; a refused candidate leaves caller-owned state unchanged. A rate law contributes one physical interval per increment, including endpoint stresses.

Focused verification covers cyclic axial response and permanent set, complete distorted-element nonproportional Chaboche tangent and force equilibrium, distinct point histories, heterogeneous yielding, actual-time overstress, refused-step retry, and geometry/domain/budget admission. The nearly incompressible axially constrained Lamé-cylinder displacement error falls from 98.3257% (standard 2x8x1) to 1.86313% (B-bar 2x8x1) and 0.500276% (B-bar 4x16x1) under the original physical force gate. This is small-strain volumetric-locking mitigation, not a general stable mixed pressure element, global injectivity certificate, finite-strain model, general shear-locking fix or plastic-limit convergence claim.


## Exact-branch Chaboche fibers and structural source receipts

`FiberLaw::chaboche_branches` and `FiberMaterial::chaboche_branches` admit one through eight native axial Chaboche branches. One branch retains the existing variant; a bounded typed dispatcher selects the actual two-through-eight branch law and matching history. Response/commit allocate no heap storage and do not execute padded inactive branches. Enum storage is bounded by its largest variant; this is not a branch-proportional memory claim. Native return, tangent, tensor-compatible lateral flow and source coefficients are unchanged.

`FiberMaterial::resolve_chaboche_branches` resolves the ordered independent C_i/gamma_i roles, Voce/Prager data and density together. It retains the complete original bundle and an explicit immutable branch-index-to-source-role mapping, including N=1. Density names must differ from every constitutive role; a conflicting binding returns a typed refusal before lookup. The material-card owner independently rejects same-name claims with incompatible dimensions. Existing Norton, single-branch, Voce, J2 and RC variants remain available.

Checked section responses refuse material-domain and branch/history-shape mismatches without changing accepted history. Raw multi-branch responses return NaNs on these failures, including newly introduced single/multi Chaboche mismatches. Existing force-based beam and static-frame owners use the same retained fiber histories. RC reinforcement retains the original geometry-consistent density and line mass.

Seven focused regressions cover every supported branch count through cycling, independent three-scale bending/tangent/backstress formulas and reversal, failed-response isolation, complete ordered source metadata and RC mass, custom-density alias refusal, sourced beam/frame loaded and permanent deflection with reversal, and an independent one-interval Norton relaxation oracle. Existing fiber plasticity regressions remain required. These are fixed-condition small-strain cyclic laws, not thermal Chaboche history transfer, new frame kinematics, damage or fatigue calibration.


## Enclosed gas-cavity geometry on deforming solids

`hyper3::pressure::cavity::ClosedTetCavity` selects one connected closed inner
boundary of an already admitted tetrahedral solid. Enclosed gas volume is
separate from total material volume. Gas-outward orientation supplies the volume
and its derivatives while canonical solid-face pressure signs are preserved.
Evaluated geometry stores compact boundary-node `dV/du` and analytic
O(boundary faces) Hessian-vector actions. Positive differential gas pressure
produces `+p grad(V)`; held-pressure endpoint work is exactly `p ΔV`.
The generated canonical face field feeds the existing static and dynamic
pressure owners without a second solid solver or equation of state.
The cavity also exposes `adiabatic_gas_load`, binding an original
`ClosedGasParcel` at the undeformed cavity volume to the exact current cavity
geometry. Every evaluation follows the source-owned isentrope to the trial
volume and pulls the same pressure and acoustic bulk modulus through
`GasMechanicalLoad`. The returned full-solid force Jacobian combines the
matrix-free geometry Hessian with dp/dV; its conservative contribution is
`Delta U_gas + p_ext Delta V`. This bridge uses no fixed gamma, does not freeze
pressure or mutate the caller parcel, and adds neither heat nor time integration.
The existing 256 compact-coordinate material pullback bound still applies.


Preparation rejects duplicate, out-of-range, open, disconnected and outer-shell
selections. Each evaluation admits every material tetrahedron and requires
positive current cavity volume. Relative reference and displacement coordinates
preserve a large common translation. These local checks do not certify
self-intersection, global injectivity or contact. Geometry supports the admitted
boundary size; `GasMechanicalLoad` retains its separate 256-coordinate chart
bound. Gas evolution, heat exchange, mixed incompressibility and time integration
remain with their respective owners.

Five focused regression groups cover independent volume/cofactor derivatives,
work and pressure signs, matrix-free tangent action and symmetry, objectivity,
actual native static-solver compression of a hollow twelve-tetrahedron specimen,
and invalid geometry/selection/cancellation. Independent analytical and NumPy
calculations verify the fixture and its tangent; the authored native Rust tests
have not been executed in this disk-constrained session.


### Simultaneous creep and yielding in native frames

`FiberLaw::NortonJ2` and `FiberState::NortonJ2` retain the original free-sided
joint axial law, separate signed creep/plastic strains and both accumulated
paths. `Frame2d::solve_creep_step` shares its existing section/member/global
Newton solve, reference quadrature, physical-duration checks, small-geometry
admission and work gates across Norton and joint Norton/J2 fibers. Repeated
section and equilibrium trials borrow the same accepted history; acceptance
advances both memories once.

`FrameCreepEnergy` reports separate physical creep and plastic dissipation.
Energetic hardening remains storage and numerical loss never becomes heat.
`creep_plastic_fiber` reads accepted histories without a constitutive update.
`with_creep_laws` accepts either axial family and preserves the original family,
elastic/hardening storage, kinetic coefficients, complete card/claim/query
identities and fixed external coordinates. Yield may refresh with temperature
only through its original pinned source. Authored and sourced histories cannot
be silently exchanged.

Untimed static, displacement, modal and inertial routes refuse both timed
families. Corotational geometry, finite deformation, thermal strain and reciprocal
heating remain outside this frame account. The `frame_creep_plastic` target
contains independent axial load/hold/unload references, nonuniform-temperature
cantilever curvature and tip-displacement references, and immutable-history
checks for unsupported routes, source replacement, inconsistent clocks and
cancellation. These manufactured fixtures are implementation checks, not
measured service-life or structural-safety qualification.

Verification for this increment: the complete native material/cut-cell/solid
source chain and both frame integration targets passed pinned-compiler metadata
checks. The joint material behavior passed its 22 focused executable regressions;
the subsequent expanded material suite passed 56/56. Native frame tests were
typechecked, not numerically executed. This is not the full Cargo/DSR gate.


### Thermal expansion and reciprocal heat in native creep intervals

The same generic `SmallStrainTetProblem::solve_creep_from` accepts
`ThermalCreep` with its complete temperature-bearing history. It retains the
original sparse equilibrium solve, common physical duration, total displacement
gradient bounds, frozen trial histories and endpoint constitutive consistency
gates. `CreepTetEnergy` additionally integrates the material's change of T times
mechanical entropy and signed reversible caloric heat over the original reference
tetrahedral volumes. The complete applied/support work balance includes these
terms exactly once. Plain Norton and Norton/J2 laws retain exact zero values for
both new accounts. Thermal transfer uses the per-element total caloric heat,
including reciprocal exchange; numerical loss and hardening storage remain separate.

The complete native source chain and the original/joint frame and tetrahedral
integration targets typechecked successfully. Those four native integration targets
were not numerically executed in this validation; the material behavior suite was.


## Native tetrahedra with thermally transported Chaboche history

The existing generic SmallStrainTetProblem directly consumes ThermalChaboche with its complete ChabocheState and nonsymmetric tangent. Each resolved old-to-new temperature law represents one physical solve_from increment; all Newton trials retain old state. For a rate-independent source, a new-to-new law recovers accepted endpoint tractions without repeating thermal transport. Timed sources retain the accepted stress; reevaluation represents another physical interval. tests/plastic_thermal_chaboche.rs verifies a six-tetrahedron held member through heating, cooling and compressive reversal, retained separate branch histories, free-row equilibrium and native per-element thermal energy. A nonproportional traction case independently checks the full nonsymmetric assembled tangent. Prescribed temperature changes reaction and loaded lateral displacement through actual constitutive flow. These tests establish a bounded small-strain prescribed-temperature capability, not an automatically heat-coupled Chaboche solid, finite-deformation or experimental qualification.


## Timed cyclic plasticity in the native tetrahedral consumer

The existing generic tetrahedral owner also advances source-pinned timed thermal Chaboche histories through heating, tensile/compressive holds and reversal. The timed regression in tests/plastic_thermal_chaboche.rs uses affine Voce Q(T), independent affine branch moduli and n=3 Perzyna kinetics. Independent endpoint flow, branch-recovery and viscous heat, coefficient storage and the first law are checked alongside force equilibrium and support work. Repeated frozen-history solves are identical; a shorter physical interval relaxes less. The test uses prescribed temperature and does not claim automatic Chaboche-state coupling in ThermoplasticSolid's heat owner, finite deformation or experimental validation.

## Retained finite-strain viscoelastic tetrahedra

`hyper3::visco::FiniteViscoTetSolid` consumes uniform or per-element original
`ViscoHyperelastic` laws. `try_new` admits the native reference geometry and
relaxed history. `advance(dt_s, supports, loads, body_load, settings, controls,
cx)` solves one actual quasistatic endpoint through the existing generic
`HyperTetProblem`, with its original force-equilibrium gate, consistent tangent,
linear solver and globalization. Every trial borrows the same accepted
`ViscoState` and physical interval; no trial or load-search iteration advances
viscous memory. Original reference geometry is prepared once and shared between
accepted candidates. Regional spectra are never averaged.

Displacement, complete element histories, Piola stress, physical storage and
elapsed time change together only after the native solve and physical work
gates pass. Dropping or refusing a candidate leaves the previous solid intact.
Absolute endpoint supports and loads permit continued loading, unloading,
stress relaxation and creep; `dt_s = 0` is the native elastic jump. One native
load endpoint is required. The caller controls physical interval subdivision
and explicit finite work tolerances.

The accepted energy record distinguishes physical storage before/after,
nonnegative native relaxation loss, native frozen-state loading work, reduced
incremental potential, actual internal/dead-load/support endpoint work, native
material balance and signed endpoint integration defect. The potential is used
only by Newton. Native work minus storage change and relaxation loss is checked
independently; assembled potential and support/load equilibrium work are also
checked. The signed endpoint defect is actual internal work minus physical
storage change and native loss. It has its own explicit allowance and is never
deposited as heat or hidden by clipping dissipation.

`tests/hyper_visco.rs` covers an affine relaxation patch with a free interior
node against native material-point histories, maintained-load creep with
physical work/loss accounts, recovery of the original elastic solver when there
are no Maxwell branches, and heterogeneous compatible histories with failed-step
rollback and deterministic retry. `cargo run -p fs-solid --example
viscoelastic_creep` reports extension, maintained force, storage, interval loss
and endpoint defect for an explicitly synthetic constrained-lateral coupon.

The base `FiniteViscoTetSolid::advance` path is isothermal quasistatic finite-bulk
P1 tetrahedra. Volumetric locking is not mitigated; this is not the Taylor–Hood
mixed/incompressible formulation from IRP M1.3-W2b. The motion owner below adds
inertia. Contact, evolving temperature, source-parameter resolution, general
time/mesh convergence and experimental qualification require their respective
owners and evidence. No closure of that full mixed-element bead is claimed.

### Native motion with retained viscoelastic memory

`hyper3::visco::dynamics::FiniteViscoTetMotion::try_new` binds the complete
accepted solid to explicit original reference density and initial nodal velocity.
Consistent `P1ReferenceMass` is prepared once and shared between candidates;
current volume and viscous deformation do not change mass. Previously relaxed
or deformed material histories and elapsed time are retained.

`try_new_with_nodal_density` uses the same native reference-mass owner for a
positive P1 nodal density. It retains the complete original field and exposes
cell means only as inventory diagnostics. `reference_acceleration_load_n`
projects acceleration through this exact matrix, including supported rows.
Candidates share the mass and preserve all accepted Maxwell histories; neither
viscous deformation nor a new temperature query changes the closed inventory.
G1/G3 regressions in `tests/hyper_visco_dynamics/nodal_density.rs` cover graded
gravity/free flight, a relaxed-memory oscillator with analytic recurrence,
kinetic-loss rollback, uniform bit parity and invalid-input refusal.

`advance` uses the same positive physical interval for the native branch update
and `BackwardEulerInertia`. The original hyper3 Newton owner assembles material
and inertial forces/tangents together, uses its physical free-row momentum
stopping rule, and includes consistent-mass coupling in support reactions.
Free bodies require no artificial supports. There is no alternate motion solver
or lumped-mass approximation. The quasistatic endpoint remains unchanged when
inertia is absent; its zero-time elastic jump remains available.

Tiny-strain energy subtraction can obscure a final Newton decrement. An
otherwise rejected inertial trial may be accepted only when its supports are
fully applied and its actual free-force residual already satisfies the original
convergence gate using the trial's native material/inertia/load normalization.
The next Newton iteration still certifies that same gate. No residual tolerance,
material-domain condition or physical work allowance is relaxed.

The returned step distinguishes physical stored-energy change, native viscous
loss, signed material endpoint defect, kinetic-energy change, nonnegative
backward-Euler numerical kinetic loss, and momentum/kinetic balance residuals.
Their sum closes actual applied/support endpoint work. Numerical kinetic loss
has its own finite caller allowance and is never material heat. A refused solve,
work gate, numerical-loss budget or cancellation returns no next motion and
leaves displacement, velocity, time, original mass and all viscous histories
unchanged. Dead loads are held over each interval; no follower/contact load,
mass transfer or automatic thermal coupling is introduced.

`tests/hyper_visco_dynamics.rs` checks unconstrained rigid acceleration against
the independent discrete momentum solution; a deforming oscillator against
the small-strain scalar Maxwell/backward-Euler recurrence with the original
consistent mass; three timestep rungs against the continuous elastic oscillator;
and post-solve kinetic-budget rollback with deterministic retry. The comparison
explicitly bounds the small finite-strain correction. These are local model/time
checks, not a general finite-strain trajectory or spatial-convergence certificate.
`cargo run -p fs-solid --example viscoelastic_creep -- --dynamic` releases a
synthetic initially stretched coupon and reports displacement, velocity,
physical storage/kinetic energy, viscous loss and both numerical work terms.

## Finite-strain Taylor–Hood displacement/pressure solids

`hyper3::mixed::MixedUpProblem` solves P2 Tet10 displacements with one continuous
P1 pressure per original corner vertex. Uniform, per-element and per-quadrature
mixed material fields retain their original laws; the material owner supplies
all four consistent saddle-point tangent blocks. Exact incompressibility is
weakly imposed in that pressure space. Original finite-bulk mode preserves the
material's admitted volume/pressure dual, including its domain restrictions.

The original Tet10 topology/map and explicit Duffy rules prepare shared immutable
reference geometry. Four through twelve points per collapsed direction are
admitted: the nonlinear pressure constraint needs the quartic rule, not the
linear-stiffness minimum. Curved geometry and nonlinear energies still require
quadrature refinement. Tangents are accumulated as 34×34 element blocks and
scattered once per element; repeated Newton actions reuse the assembled matrix,
its original sparse LDL factor, and one symbolic analysis for the fixed pattern.
The previous numeric factor is released before its replacement is allocated.

The native Newton–Krylov owner supplies residual Armijo globalization. Required
length and stress scales define positive work-conjugate congruence scaling;
force [N] and constraint-volume [m³] residuals are not directly mixed. The same
physical force and complete pressure residuals are checked at acceptance.
Native LDL uses 2×2 pivots and requires inertia (n_free_u positive, n_free_p negative),
including for an initially equilibrated or zero-iteration endpoint. An unstable
or singular candidate refuses. Numeric factor payload, sparse assembly payload,
mesh/free-DOF counts and quadrature-point count have explicit caller caps.
Cancellation is checked during assembly and around the bounded native factor.

Pressure normally follows traction equilibrium. An optional exact-mode pressure
pin is accepted only when the selected component's constant-pressure field is a
true null mode on free displacement rows; the omitted constraint is also checked
at the accepted endpoint. Finite-bulk pins and physically overconstraining pins
refuse. The public raw evaluation retains constrained rows, reactions, physical
residuals, the full tangent and mixed saddle potential. That potential is not
silently interpreted as viscous stored energy.

`tests/mixed_up.rs` verifies independent uniaxial traction/stretch/pressure,
original bulk/shear ratios 1 through 1e8, full tangent derivatives and objectivity,
quartic constraint integration, pressure gauges, budget/cancellation boundaries,
and rejection of an already equilibrated unstable compressed bar. A non-affine
Cook panel has a three-rung refinement diagnostic and compares bulk/shear 1e6
against exact incompressibility. This bounded evidence is not the complete
Treloar/MMS/inf-sup campaign or a universal locking-free certificate; bead
`frankensim-irp-v5-instrument-realism-pexfp.61` remains open. The displacement-only
P1 hyper3 path retains its original volumetric-locking limitation.

## Retained finite-strain mixed viscoelastic solids

`hyper3::mixed::visco::MixedViscoTetSolid` owns the accepted displacement,
continuous P1 pressure, time, Piola stresses and one distinct
`ViscoRemainderState` per TET10 quadrature point. `MixedViscoLawField`
accepts uniform, per-element or per-quadrature original spectra without
averaging. Prepared P2/P1 geometry is shared across accepted endpoints.

`advance` solves one physical endpoint with accepted metrics and interval frozen
throughout the native mixed Newton solve, warm-starting displacement and pressure.
A zero interval means an instantaneous elastic jump. Numerical load continuation
must have exactly one endpoint; physical substeps remain explicit caller actions.
A returned candidate owns the complete new state. Cancellation, Newton refusal
or a failed work gate cannot mutate the accepted solid.

Physical storage is the native remainder plus `U(J(p))`, reconstructed from
pressure. It is not the original bulk energy at the quadrature deformation:
the mixed pressure constraint holds weakly. Exact incompressibility never
evaluates the unused bulk penalty, even when publishing accepted branch updates.
The report separates remainder/bulk storage, native viscous loss, frozen-history
loading work, incremental remainder and mixed saddle potentials, pressure
constraint work, actual internal/dead-load/support endpoint work, equilibrium
work residual and signed endpoint-work defect. Explicit
`FiniteViscoStepControls` gate the work residual and defect; neither numerical
defect nor saddle potential is counted as physical dissipation or heat.

`tests/mixed_visco.rs` verifies held-stretch relaxation against the native
material point, maintained-traction creep with an independent scalar equilibrium,
original-bulk storage and the non-aging elastic limit, heterogeneous spectra and
transactional work-gate refusal. Exact-mode relaxation is bitwise independent
of unused bulk moduli `1e3` and `1e300`. The runnable
`mixed_viscoelastic_creep` example applies exact quadratic-face loads and reports
stretch, pressure, storage, interval viscous loss and signed endpoint defect.
The standalone solid advance is isothermal and quasistatic; the retained motion
owner below adds consistent inertia. Neither path supplies thermal closure, a
general mesh inf-sup certificate or the full M1.3-W2b acceptance ladder.

## Consistent TET10/P1 mixed viscoelastic motion

`hyper3::mixed::visco::dynamics::MixedViscoTetMotion` binds an accepted mixed
solid to explicit original per-element density and initial nodal velocity.
Its immutable consistent reference mass is integrated from the same cached
quadratic shapes and positive reference weights. Signed corner/edge entries
are retained; row-sum lumping and deformation-dependent mass changes are not
used. The matrix is exact for affine TET10 geometry and uses the caller's
explicit quadrature approximation for curved geometry. Mass assembly has a
separate byte cap and cancellation checkpoints.

The existing native backward-Euler inertia enters the same mixed Newton solve:
`M*a` augments displacement forces and `M/dt^2` augments only the displacement
tangent. The pressure equations are unchanged and pressure carries no mass.
Residuals, final admission, support reactions and the KKT inertia check use the
complete incremental system. The reported mixed potential remains the material
saddle potential. Static entry points continue to omit inertia.

The original reference mass owner also supports signed consistent entries in
its tangent and uses `abs(M)` for its positive momentum reference bound; existing
positive P1 and Hex8 actions retain their accumulation order. Density-weighted
mass, not current volume, defines kinetic energy and momentum.

A returned motion candidate advances displacement, pressure, velocity, material
history and time together. Physical branch loss, signed material endpoint-work
defect and nonnegative backward-Euler kinetic loss remain separate. Existing
`FiniteViscoMotionControls` independently bound material work and numerical
kinetic loss. Refusal leaves the complete accepted endpoint unchanged.

Exact pressure enforces endpoint volume weakly. Velocity is the backward-Euler
displacement chord; there is no separate tangent-velocity projection. This is
a first-order convergent motion path with explicitly budgeted numerical damping,
not the instrument-render doctrine's zero-algorithmic-damping or discrete-gradient
integrator. Closed work accounts do not establish timestep accuracy or general
mesh inf-sup stability. Thermal evolution, follower loads and contact retain
their separate owners.

`tests/mixed_visco_dynamics.rs` checks exact-incompressible free-body
acceleration, a genuine quadratic shear mode against its independent scalar
backward-Euler recurrence and continuous solution on 20/40/80 intervals,
released finite-bulk versus exact-pressure viscoelastic motion, signed moving-
support reactions and complete rollback after numerical-loss refusal. Native
inertia unit tests check rational P2 coefficients, linear-velocity kinetic
moments, tangent finite differences, the absolute-mass bound, original positive
P1/Hex8 arithmetic, budget and cancellation. The existing creep example also
offers `--dynamic` to release the deformed body and report physical and
numerical losses independently.

### Instantaneous impulses with retained mixed viscoelastic memory

`MixedViscoTetMotion::impulse_candidate` transfers nodal impulses through the
original signed consistent P2 mass for `VolumetricTreatment::OriginalMaterial`.
The same shared mass owner supplies bounded CG, repeated-load accumulation,
actual free-row/linear/angular momentum gates, kinetic change and signed
support work. Held DOFs retain their current velocity, including nonzero
support motion. Angular momentum uses the current deformed geometry.

The complete mixed solid survives unchanged: continuous pressure, displacement,
all original quadrature Maxwell metrics and accumulated relaxation loss,
stresses, remainder and pressure-reconstructed bulk energy, physical clock,
source laws, reference mesh and density. Only velocity changes. A later native
`advance` continues these accepted histories. An instantaneous event creates no
material relaxation loss or heat, and the caller retains ownership of its
contact/restitution law and opposite body's impulse.

This mass-only entry returns the typed
`MixedViscoTetError::ExactIncompressibleImpulse` for an exact pressure
constraint. It does not assign pressure a mass or silently project an
incompatible velocity. All solve, input, budget, balance and cancellation
refusals leave the accepted motion unchanged.

`MixedViscoTetMotion::constrained_impulse_candidate` supplies the separate
exact-incompressible path through the existing shared P2/P1 pressure-impulse
owner. Each original accepted quadrature deformation supplies its endpoint
cofactor; the shared native mass/pressure KKT solve retains the original
supports, weak volume rows, physical scales and sparse budgets. Internal
pressure impulse has independently gated zero resultant, torque and work.
Its Pa s multipliers are reported separately and never overwrite physical
pressure or any material history. The finite-compliance formulation returns
the typed `FiniteCompliancePressureImpulse` refusal from this entry.

The incoming velocity must already meet the caller's weak-volume-rate allowance.
A velocity accepted through backward-Euler position chords is not automatically
tangent-compatible; the event refuses a large pre-existing rate defect instead
of concealing it through a dissipative projection. A pressure-impulse gauge is
admitted only for a true connected constant null mode. Shared pressure and
mass accounting gate the actual rounded successor before the same atomic
complete-history transfer used by the finite-compliance entry.

`Tet10ImpulsePoint` names a point by its original cell and four barycentric
coordinates. The native quadratic shape owner supplies all ten signed weights
for current isoparametric position, velocity and transpose impulse distribution;
negative corner weights are retained. `point_kinematics` samples either pressure
formulation. `point_directional_inverse_mass_per_kg` uses the original supported
consistent mass response without a material trial or history clone, while
`point_impulse_candidate` delegates the complete-history nodal event. These two
mass-only operations explicitly refuse exact incompressibility. The shared P1
and P2 trace helpers retain conjugate point work and torque on curved geometry;
they do not replace the original mass by a point mass. The finite-compliance
entry is connected to the existing rigid impact owner in `fs-couple`.

The focused `mixed_visco_dynamics/impulse.rs` regressions use independent
rational P2 integrals and signed corner loads, deformed angular impulse,
moving-support work, nonzero retained pressure and relaxed Maxwell histories,
subsequent native motion, typed formulation refusal and transactional retry.
The exact-mode extension checks an independently integrated pressure impulse
on an already stretched, relaxed body, retained translation and Maxwell
continuation, and incompatible velocity/gauge/factor-budget refusals. A genuine
pressure-gauge regression also checks cancelling applied and moving-support
work with unchanged kinetic energy and retained relaxation history.
The material-point extension checks signed quadratic interpolation and loads,
supported mobility against independent P2 integrals, torque, work and retained
history. The point owner's polynomial fixture also checks curved reference
geometry and virtual work independently of the finite-element solve.
Native execution remains pending the existing dependency-build blocker; this
addition is not a completed Cargo or DSR gate.

## Deformable cohesive dynamics with retained inertia

`cohesive::dynamics::CohesiveMotion` adopts an existing `CohesiveState`, its
original reduced consistent mass operator, explicit velocity and clock. It
retains the same elastic operator, sparse jump maps, patch areas, ordered source
bundles and accepted fracture histories. There is no second interface assembly
or independent damage state. Matching and nonmatching P1 traces use their
existing work-conjugate maps. The mass must use exactly the same reduced DOFs,
be symmetric positive definite and remain unchanged; the elastic operator must
be symmetric positive semidefinite and unchanged. These are requirements on
the original operator owner, not conclusions from dimension or energy samples.
Rigid-body modes do not require artificial elastic restraints in this path.
The bulk potential is `u^T K u/2` in the admitted reduced coordinates; an
additional prestress potential or moving essential boundary data requires its
own mechanics, inertia and support-work treatment.

For a physical interval `h`, the existing Newton/FGMRES owner solves
`M (v1-v0)/h + K (u0+u1)/2 + sum B^T A t_mean = f_mean`, with
`v1 = 2 (u1-u0)/h - v0`. Each trial starts from the same accepted material
history. `t_mean` is the original analytic vector integral along the straight
three-component jump path and includes surviving-bond viscosity. The endpoint
traction tangent is not substituted for the derivative of this path mean.
The existing fallible numerical-difference adapter supplies Jacobian actions
with an explicit displacement stencil; all actual residual evaluations count
against the declared solve budget. A source or stencil refusal refuses the
whole increment, and the complete physical endpoint is independently replayed
before acceptance.

The step reports original kinetic, bulk and cohesive endpoint energies,
intrinsic fracture loss, physical bond-viscosity loss, optional damaged-area
friction loss, external force work and
prescribed-jump work separately. Prescribed work is
`sum A t_mean dot (jump1-jump0-B(u1-u0))`, so adoption of a previously loaded
solid neither guesses its old actuator offset nor drops that actuator's work.
The full gate is external plus prescribed work minus kinetic/bulk/cohesive
storage changes, fracture, viscosity and friction. No numerical residual is relabeled as
physical dissipation or routed to a thermal field. Signed roundoff-scale bulk
energies are retained rather than clamped, with materially negative values
refused under the explicit energy allowance.

The accepted endpoint replay also reports `DynamicsPatchDissipation` in original
patch order, with actual area-integrated fracture, viscosity and friction losses.
Its sums use the same arithmetic as the global loss accounts. A downstream
thermal owner can identify the original patch/source through the read-only
`CohesiveState::system` accessor; these reports alone assign no heat fraction.

Complete failure removes tensile and shear resistance through the original
material law while retaining mass, velocity, history and clock. Original normal
compression remains available on reclosure. The default is the existing
frictionless penalty law; the explicit sourced damaged-area extension below
adds residual shear resistance. Neither is exact nonpenetration or newly
discovered contact. Immutable
candidate publication makes late cancellation, material-domain refusal, failed
work gates and exhausted solves leave the entire accepted trajectory intact.

This is fixed-reference small-strain dynamics with fixed material/source
conditions during each mechanical increment, fixed frames and fixed topology.
The explicit temperature transfer below may advance the cohesive source and
history between increments. It does not update bulk coefficients, impose
moving reduced constraints, discover contact or provide a unique/stable
nonlinear root.
Prescribed interface offsets and reduced external forces are explicit inputs;
their physical construction and resolution of reversals, crack events and
softening branches remain caller responsibilities. Smooth elastic motion has
the midpoint second-order limit, but a closed energy account alone does not
establish temporal accuracy through fracture.

`tests/cohesive_dynamics.rs` specifies original P1 consistent-mass frequency
and refinement, independent scalar actuator/viscosity work, release of an
already softened solid, and source-bound three-component debonding of two
unrestrained deformable tetrahedra. Its postfailure checks retain each body's
momentum, drive physical compressive reclosure and retain the same sources and
damage. Late cancellation, a material-domain refusal and work-budget exhaustion
also retain the complete accepted trajectory. These focused native tests do
not constitute a dynamic crack-propagation, finite-strain or service-validation
campaign.

### Instantaneous cohesive momentum transfer

`CohesiveMotion::impulse_candidate` applies an explicitly supplied generalized
impulse in the original free-coordinate ordering by solving `M delta_v = J`
with the native checked CG owner and supplied SPD preconditioner. It uses
operator actions and bounded Krylov storage, without a dense inverse or a
replacement mass. The original immutable motion supplies displacement, damage,
source-bound viscosity, friction/contact bindings, clock, step count and controls.
Only an accepted velocity and its incremental event account are published.

The event independently rechecks `M (v_after - v_before) - J` using the actual
representable returned velocity. Its signed work is
`((v_before + v_after)/2) dot J`; kinetic change uses the corresponding mass
increment. Compensated accumulation and a work tolerance scaled to the event
avoid hiding a lost impulse in a large unchanged background energy. Numerical
workspace, solver iterations, cancellation and the original velocity limit are
explicit gates. Refusals leave the accepted state intact.

No displacement path or physical time elapses, so the event does not invent
fracture, viscous or frictional loss. The next physical interval retains those
original owners. An affected Coulomb tangential trace reconsiders inherited
sliding occupancy using that interval's actual incoming velocity, including its
prescribed rate; prior reactions and material history survive the impulse.
An unaffected trace retains its existing occupancy. Accepted finite-time replay
clears this pending classification, and temperature transfer preserves it.

The caller owns impulse generation, other-body reaction, restitution and any
impact energy loss, and adds the signed event work to its trajectory account.
These reduced-coordinate balances do not establish Cartesian angular momentum
without the original geometric trace. Moving support impulses, new collision
topology and automatic heat routing are outside this operation.

### Damaged-area friction in the original deformable dynamics

`CohesiveMotion::with_damaged_area_friction` binds one explicit regularization
speed per selected original patch. The kinetic coefficient is read from the
canonical `cohesive_residual_kinetic_friction` property already retained in
that patch's complete `ResolvedBilinearCohesive` source bundle. The receiving
owner requires a nonnegative dimensionless coefficient at fixed external
conditions; it never infers friction from cohesive strength or fracture energy.
Unselected patches keep their original response, and zero mu retains the
original arithmetic. Immutable shared configuration avoids per-step copying of
the patch coefficient array.

The existing material owner supplies the analytic path measure
`P = integral_0^1 d(s) Kc max(-gn(s),0) ds`, including changing damage and
compressive reclosure. With the actual straight-path tangential velocity `v`,
the mean resisting shear is `mu P v/max(|v|,v_reg)`. It is assembled through the
same `B^T A` map as the original bond and solved with the original mass,
elasticity, cohesive history, viscosity and Newton/FGMRES owner. A signed net
normal adhesive force or endpoint `d*p` is not substituted for that integral.
The original material's damage variable remains the sole damage history.

Friction work is the mean added shear dotted with the actual tangential jump,
integrated over original area. It is nonnegative, reported separately from
fracture and surviving-bond viscosity, and debited exactly once in the full
mechanical energy audit. Prescribed-interface work uses the total bond plus
friction traction, so a moving support supplies its actual contribution. No
heat destination or entropy claim is assigned implicitly. Source-temperature
transfer refreshes the same selected coefficient in the original successor
bundle while preserving regularization speed, damage, mass, velocity, time and
step count. Changing mu at fixed jump creates no friction work or heat.

This is velocity-regularized kinetic Coulomb on a declared damaged-area mixture,
not ideal static friction, a tangential elastic spring, roughness calibration or
new contact discovery. The speed is a numerical closure in m/s, independent of
the timestep; timestep refinement still resolves reversals and stopping. The
fixed-reference and original material-domain limitations remain in force.

`tests/cohesive_dynamics_friction.rs` compares actual native inertial response
against independent scalar equations for intact, partially damaged and fully
failed compressed bonds, including surviving-bond viscosity. It checks an
independent analytic integral while shear initiates and grows damage under
compression, rotated prescribed-slip work, zero-mu equivalence, low-speed
physical-time regularization, source-temperature substitution, missing/invalid
data and cancellation before final publication. These manufactured-data tests
establish implementation behavior, not measured adhesive-friction calibration.

### Static sticking and kinetic breakaway of reclosed deformable bonds

`CohesiveMotion::with_damaged_area_coulomb` selects unregularized static/kinetic
friction on the same original damaged area. It requires both canonical residual
friction properties in the original source bundle, dimensionless at fixed
external conditions, with `0 <= mu_kinetic <= mu_static`. Initial closed-contact
occupancy is explicit; accepted endpoint motion and occupancy determine whether
the next interval may test static capture. This owner is mutually exclusive
with the velocity-regularized extension. Neither introduces a second cohesive
history or changes the original damage, pressure or surviving-bond viscosity.

Two scaled tangential traction coordinates per nonzero-capacity selected patch
join the existing midpoint Newton/FGMRES system. The natural map is
`t - project_disk(t + k*delta_g, mu*P) = 0`, using the common `fs-material`
projection already used by rigid and finite-patch damaged-area owners. `P` is
the original analytic integral of damage times compressive pressure over the
actual straight cohesive path, including changing damage and closure. The
positive numerical `k` conditions the map and scaled unknowns; it is not an
interface spring and stores no energy. The native body mass, elasticity,
work-conjugate trace map and full coupling between patches remain in the solve.
An initially static candidate with resolved slip or endpoint velocity is re-solved at its sourced
kinetic capacity. All passes share the original nonlinear iteration and
physical-residual budgets; Krylov admission counts the augmented dimension.

The accepted physical audit checks the static reaction disk or the direction
and magnitude of kinetic traction. Sticking has exactly zero physical friction
loss only when slip and both endpoint velocities are unresolved; finite
tolerance-scale traction times slip remains in the work/energy
defect. Sliding loss is `A * mu_kinetic * P * |delta_g|`, and the same added
traction contributes to prescribed-interface work. Original fracture and bond
viscosity losses remain separate. Temperature transfer refreshes both original
coefficients while retaining damage, source identity, geometry, inertia,
velocity, clock and occupancy. Zero coefficients preserve original arithmetic.

This is a bounded straight-path rigid-friction specialization on fixed damaged
area, without tangential elastic slip storage or newly discovered contact.
Resolved reversal of tangential endpoint velocities while the pressure path is
active refuses the whole increment, including nonzero-net-slip reversals;
time-step/event refinement must resolve stopping before physical heat is
accepted. Nonlinear convergence does not establish time accuracy or uniqueness
of redundant static reactions, and singular or exhausted native solves refuse.

The Coulomb cases in `tests/cohesive_dynamics_friction/coulomb.rs` check static
hold and kinetic breakaway against independent inertial equations for partial
and complete damage, a rotated two-component reaction disk, changing damage
and opening-path pressure, actual prescribed work, source-temperature threshold
changes, physical kinetic heat below the displacement tolerance, zero-coefficient
parity, source/budget/cancellation refusal and both
zero- and nonzero-net-slip stopping refusals. Source data are manufactured;
these tests do not establish specimen calibration.

### Native cohesive intervals across Coulomb stopping and resticking

`CohesiveMotion::advance_interval` completes a requested duration under constant
external force and affine prescribed jumps using the original mechanical owner.
Each substep retains the original maximum step, constitutive laws, histories,
momentum equation and physical acceptance tolerances. Only the typed
`UnresolvedFrictionReversal { patch }` refusal triggers bracketing: shorter
candidates are solved from the same accepted predecessor until an admitted
candidate stops inside the original velocity-resolution band or resolves the
declared time bracket, then the remaining interval is advanced. Other refusals
propagate unchanged. A rejected trial's damage, motion
and losses are never adopted. No velocity projection, stopping impulse, extra
friction spring or numerical heating is introduced.

The original endpoint velocity tolerance determines static capture. Residual
kinetic energy stays in the actual native state and first-law account. This is
bounded continuation of the discrete straight-path model, not an exact stopping
time, earliest-event certificate, monotone-guard proof or unique softening-path
claim. Source conditions stay fixed over the interval; affine prescribed motion
and constant force are explicit restrictions. Native time resolution, conditioning,
source domains and unresolved events can still refuse.

Independent caps bound attempted solves, accepted substeps, cumulative physical
residual evaluations and nonlinear iterations. Each attempt reserves its native
allowance; successful work refunds unused counts while a refused solve charges
the whole reserved allowance. Per-residual material work retains the original
patch bound. Cancellation or exhaustion returns no partially accepted interval.
Memory retains only the current motion and best bracket candidate, not an entire
state trajectory. Original fracture, bond viscosity, friction and prescribed
work are accumulated separately, including per-patch totals. The interval's
endpoint first law is checked against the sum of original accepted tolerances.

The interval also retains each original patch's maximum normal jump over its
incoming endpoint and all accepted straight local paths. Rejected trials never
enter these bounds; a final closed endpoint cannot erase earlier opening. The
final remainder uses the exact requested clock coordinate and the unchanged
positive native timestep, avoiding an extra rounded time addition after a stop.
Ordinary single-step clock arithmetic is unchanged.

For resolved incoming slip, the augmented Coulomb solve starts its reaction
guess on the sourced kinetic disk at the incoming damaged compression and in
the actual displacement predictor's slip direction. This avoids a zero-reaction
line-search barrier for tiny sliding increments without another material trial
or solve. The original integrated path pressure, nonlinear equations, work
limits and physical audits still determine the accepted reaction and losses;
stationary capture retains its previous reaction predictor.

`tests/cohesive_dynamics_friction/interval.rs` compares an unaligned scalar stop
with independent constant-friction stopping distance and kinetic-energy loss,
then holds above kinetic capacity and breaks away above static capacity. It
also checks accumulated original fracture/viscosity/friction on prescribed
paths, native replay, work bounds, late cancellation and immutable rollback.
Path maxima are checked against native replay and incoming opening followed by
closure; nonzero-origin stopping retains the exact requested final clock.
A small native kinetic step is checked against constant-force displacement,
endpoint velocity and frictional-work formulas before completing the stop.

### Temperature-dependent adhesion in deformable equilibrium and dynamics

`CohesiveSystem::with_temperature_sources` prepares temperature successors on
the original patches. `CohesiveState::transfer_temperature` accepts those
coefficients only after matching the original elastic operator, DOFs, patch
identities, areas and complete sparse trace maps. Every ordered source retains
the original card, selected claims, material surfaces and nonthermal query
coordinates through `ResolvedBilinearCohesive::require_temperature_successor`.
Unsourced interfaces and substituted geometry/operators refuse.

Each quadrature point uses the existing material owner's fixed-jump coefficient
change followed by irreversible relaxation. Common-mode and frozen-physical-mode
BK histories retain their original interpretation; cooling cannot heal bonds or
rescale past fracture work. Area-integrated signed parameter work, additional
damage loss and storage change are reported separately and gated by
`parameter work - storage change - damage loss`. No displacement or dashpot
work occurs during transfer, and neither parameter work nor damage loss is
automatically assigned to heat. Bulk elasticity and thermal expansion remain
outside this transfer.

`CohesiveMotion::transfer_temperature` retains the original mass, velocity,
physical time, accepted step count and numerical controls. The next native
mechanical solve responds to the resulting interface tractions. The transferred
state is not asserted to be equilibrated at its new temperature. All sources,
histories and integrated accounts publish atomically; per-site and bounded
vector polls preserve cancellation without partial accepted state.

`tests/cohesive_temperature.rs` checks independent heated-bond storage/damage
identities, area integration across two original sites, cooling irreversibility,
native continuation to failure, changed-frequency inertial motion with retained
mass, trajectory-budget retention, source/geometry/operator refusal and late
cancellation. These manufactured source curves verify the integration; they
are not physical validation of a temperature-dependent adhesive specimen or a
coupled thermal-transport/heat-partition model.

## Conservative finite-strain elastic motion

`hyper3::mixed::conservative::{ConservativeTetDefinition, ConservativeTetSolid}`
retain quadratic TET10 displacement and velocity, physical stored energy, elapsed
time, original element material identities, and a shared consistent reference
mass. `FullMaterial` evaluates the original finite-strain energy through the
objective C-based discrete gradient in `fs-material`. Initial displacement and
velocity are explicit; preparation does not silently relax a supplied state.

`ConservativeTetSolid::impulse_candidate` applies zero-duration nodal impulses
through the same full consistent reference mass for `FullMaterial` and
`OriginalMaterial`. Current deformed positions supply the angular-momentum
arms. Held velocity coordinates retain their actual velocity, and the shared
impulse owner reports their reaction impulses and signed work. Displacement,
physical time, material points, original density, retained bulk volumes and
algebraic pressure remain unchanged. Applied and support impulse work enter
the original compensated cumulative account; the initial energy is not reset.
Subsequent native intervals therefore retain the entire physical work history.

The same public method retains the native pressure-impulse constraint solve
for `ExactIncompressible`. Compressible formulations use only `controls.balance`
and require no pressure gauge; pressure/rate arrays are empty and pressure
momentum/work is zero. `MixedImpulseReport::kkt_inertia` is `Some` only when the
native constrained factor exists, and `None` for the checked mass solve.
The caller owns restitution/contact physics and any opposite body's impulse.
G1/G3 fixtures in `tests/mixed_conservative_bulk/impulse.rs` compare independent
quadratic-simplex momentum and kinetic integrals, signed support impulses,
moving-support work, retained state, native continuation and refusal behavior.

Each interval uses centered kinematics, `v_after = 2 Delta u / dt - v_before`,
and force `2 M (Delta u - dt v_before) / dt^2`. The native consistent P2 mass
retains its signed entries and quadratic body-force weights. The inertia work
equals the kinetic-energy increment with zero algorithmic loss. This is a
separate `MidpointInertia` path; existing backward-Euler motion keeps its
original dissipative semantics.

`ConservativeTetDefinition::prepare_with_reference_density` binds four positive
P1 corner densities per original tetrahedron to that same signed consistent P2
mass. The cached barycentric pressure basis supplies the density samples; no
geometry reconstruction, row lumping or new time integrator is introduced.
The immutable field and assembled mass remain shared through accepted motion
for full-material, exact-incompressible and original finite-bulk formulations.
Four Duffy points per direction admit graded affine inertia; a varying profile
on a curved quadratic map requires six. Original constant cell fields retain
their scalar quadrature semantics and accumulation bits.

`reference_corner_density_kg_m3` exposes the retained field. For a graded body,
`reference_density_kg_m3` reports reference-volume-weighted cell means; physical
inertia and loads use the full field and mass operator. The body's
`reference_acceleration_load_n` projects nodal acceleration through its actual
full `M`, including prescribed coordinates and signed quadratic corner loads.
Uniform gravity and affine acceleration therefore retain the body's original
mass, force resultant and torque. The returned loads are held over a native
interval; evolving position-dependent force laws retain their own owner.

Independent barycentric polynomial integrals check every affine and curved
mass entry. Native deforming, translating and spinning trajectories check
energy and both momenta; acceleration fixtures exercise all three formulations
and constrained support reactions. Admission, scalar parity and cancellation
checks retain the original transactional boundary.

Assembly contracts each point's true generally nonsymmetric endpoint tangent
into one local displacement block before one element scatter. The native
Newton-FGMRES owner uses that actual Jacobian; only its symmetric part supplies
the sparse LDL preconditioner. The preconditioner is neither a substituted
elastic tangent nor an inertia-based stability certificate. Geometry, mass and
accepted material points are retained between candidate evaluations.

Dead nodal/body loads, moving Dirichlet supports and returned signed reactions
have explicit work accounts. Each accepted interval checks constitutive work,
equilibrium work, local total-energy balance and cumulative total-energy balance
against the caller's allowance in joules. Kinetic work and numerical loss remain
separate fields. Invalid domains, exhausted budgets, failed solves, failed work
gates or cancellation return no changed accepted state. One physical interval
requires one load step; Newton backtracking does not advance material history.

G0/G1 evidence in `tests/hyper_conservative.rs` includes independent consistent
mass and original-energy checks for a freely translating, rotating and deforming
body; exact constant-acceleration motion; moving-support work; and rejected or
cancelled candidate rollback. A genuinely nonlinear quadratic shear benchmark
uses an independent quartic Ogden oscillator, independent scalar discrete steps,
and a time-refinement ladder against a converged continuous reference. The
focused native-source lane passed four tests: the three-rung time order was
1.99 and 2.00, and the 48-step free-body energy/momentum discrepancies were at
floating-point roundoff. The native inertia unit lane passed nine tests.

`FullMaterial` is displacement-only and a very large bulk modulus may still
lock. The exact and original finite-bulk mixed formulations below supply
separate pressure treatments. Nonlinear modal reduction, contact, viscoelastic
dissipation and the frequency-exact finite-rotation nonlinear split remain
separate. Centered nonlinear time integration has phase error even when its
discrete work balance closes; the complete-coordinate exact linear owner below
has its own scope. The bounded nonlinear fixtures do not satisfy the separate
long-horizon, nonlinear frequency-matching and complete port/checkpoint
obligations of bead `frankensim-irp-v5-instrument-realism-pexfp.62`.

## Conservative exact-incompressible P2/P1 motion

`ConservativeFormulation::ExactIncompressible` extends the retained elastic
motion owner with continuous corner-node P1 interval pressure and endpoint
constraints `integral N_p (J_after - 1) dV = 0`. The original bulk-free
remainder supplies physical storage and elastic work. No bulk penalty is
evaluated, subtracted or approximated; the pressure carries neither mass nor
stored energy. Original full-invariant fibre dependence is retained.

Momentum uses the objective discrete volume derivative, while the endpoint
constraint uses the endpoint cofactor. The two Jacobian cross-blocks are therefore
distinct at finite motion. The pressure contribution to the displacement block
includes the full derivative of the discrete cofactor. Native Newton-FGMRES
retains all these nonsymmetric derivatives. Its symmetric saddle preconditioner
uses the endpoint constraint derivative in both cross-blocks, keeping redundant
endpoint constraints singular instead of hiding them by averaging the crosses.

Preparation admits both the initial weak volume and the initial weak volume
rate `integral N_p cof(F_initial):grad_X(v_initial) dV`. Their explicit
controls are `maximum_relative_volume_residual` and
`maximum_initial_relative_volume_rate_per_s`, normalized by total reference
volume. Inconsistent initial velocities are refused without a pressure impulse
or silent projection. This is an initial compatibility check; subsequent
centered velocities are not claimed to satisfy a separate exact endpoint
tangency constraint.

`maximum_point_volume_work_defect` bounds each dimensionless local
discrete-volume identity. The accepted step separately reports actual pressure
work, its discrepancy from pressure dotted with the weak-volume change, and
the endpoint pressure residuals in cubic metres. Pressure work, local and
cumulative physical energy, and weak volume are all admitted before publishing
the new state. `algorithmic_pressure_pa()` returns the last interval's
algebraic work multipliers, initialized at zero; these are not endpoint Cauchy
pressures. `pressure_nodes()` identifies their sorted corner-node layout.

Independent endpoint pressure rows on the free displacement space are required.
The final factorization checks this even for a zero-iteration Newton solution.
There is no pressure-pin/gauge API; closed fully prescribed displacement cases
with redundant volume equations refuse. This method enforces weak P1 volume
constraints, not pointwise J=1; a coarse element can retain local volume
variation even when every weak moment closes. It is not a conservative
finite-bulk mixed formulation.

G0/G1 regressions in `tests/mixed_conservative.rs` use independently integrated
Mooney-Rivlin remainder energy, rational consistent P2 mass, momenta and weak
volume. The free nonaffinely strained, spinning body changes shape substantially
over 32 intervals while preserving energy and both momenta to roundoff and
normalized weak volume within `2e-15` in this fixture. Separate checks cover
bitwise trajectory independence from unused bulk moduli `60` and `1e300`,
incompatible initial volume/velocity, invalid rate controls, and redundant
pressure-row refusal without state mutation. Full-material regression tests
remain green. `examples/conservative_elastic_motion.rs` runs both formulations
and prints deformation, energy, work, sampled J, weak volume and interval
pressure with explicit SI controls. These bounded tests do not establish a
mesh-independent inf-sup constant or the remaining long-horizon/render contract.

## Exact complete-coordinate linear solid motion

`dynamics::ExactLinearSolid` advances every free coordinate of an admitted
`TetElasticAssembly` using the existing
`fs-phs::linear_exact::prepared::PreparedQuadraticFlow`. Full consistent mass
and stiffness are retained; no modal truncation, mass lumping or frequency
retuning is introduced. State is `[u,v]`, with original physical storage
`(u^T K u + v^T M v)/2`. Native Cholesky prepares the inverse mass and its
residual is independently admitted against the original matrix.

Preparation accepts explicit initial free-coordinate displacement/velocity,
one fixed positive clock, and `ExactLinearSolidControls`: complete-coordinate
capacity, dense preparation bytes, exponential multiplication work,
operator-error allowance, work allowance in joules and accepted-step cap.
The existing complete augmented dimension is bounded by `3*N+1 <= 256`.
The prepared exponential is reused for warm intervals. Zero-stiffness directions
remain in the state and retain rigid motion and constant-acceleration response.

Force inputs, initial fields and returned mean velocities follow the source
assembly's full `free_dofs()` order. The mean velocity is `Delta u/dt`.
Held-force work is independently `f dot Delta u`; it is not inferred from
the energy difference. Interval and accumulated original physical work balances
gate publication. Invalid operators, exhausted budgets, invalid force, work
refusal and cancellation preserve every accepted state field.

The supported reduced assembly represents homogeneous stationary constraints.
This API neither accepts a full nodal-force vector that silently drops
constrained forces nor reconstructs missing support reactions or moving-support
coupling. It supplies no material damping, nonlinear finite-rotation claim or
sparse exponential scalability claim. Exactness refers to the admitted spatial
pencil and held input, not continuum accuracy or arbitrary time-varying forcing.
The nonlinear conservative owner retains its separate phase limitation.

G0/G1/G5 evidence in `tests/solid_linear_dynamics.rs` passed all five native
tests. A 5 kHz oscillator at a 48 kHz clock runs 1,000,000 accepted intervals
with pitch error below `1e-6` cent and maximum relative energy drift
`5.03e-10`; the independent midpoint formula gives `-59.1076` cents.
A coupled non-diagonal pencil retains both 5 kHz and 11 kHz analytic frequencies
and held-force work. Rational P1 and signed P2 consistent-mass pencils, zero-K
motion, capacity/domain refusal, cancellation and repeatability are covered.
These tests construct independent pencils in the native assembly type; they
do not execute the original geometry assembly frontends. The native lane
compiled the complete original fs-phs crate and was rerun against the current
shared propagator. No full-workspace Cargo result is claimed.

## Original finite-bulk conservative P2/P1 motion

`ConservativeFormulation::OriginalMaterial` extends the same quadratic solid
owner with the ORIGINAL finite-compliance bulk energy. It accepts
`DiscreteBulkSettings`, an explicit relative weak-volume allowance and a local
discrete-volume-work allowance. Quadratic and log-squared laws use their own
native secant inverses. Actual geometric and retained bulk volumes must both
remain inside the caller's positive-volume branch; no clipping, averaged modulus
or surrogate quadratic law is used.

At each reference quadrature point, accepted `Jbar_before` is frozen.
The interval P1 pressure determines `Jbar_after` through
`p_alg = Delta U / Delta Jbar`. The spatial equation is
`integral N_p (det(F_after)-Jbar_after) dV = 0`. Actual geometric pressure
work uses the objective discrete cofactor. Its two distinct pressure/displacement
cross-blocks remain unchanged, and the pressure block is the negative original
secant compliance. The symmetric preconditioner retains that pressure block,
so a finite-compliance body with fully prescribed displacement can still solve
its pressure equations. This does not add a gauge to the exact-constraint mode.

Physical storage is the original remainder at actual F plus original U at
retained Jbar. `bulk_volume_ratios()` exposes those values in reference
element/quadrature order; `bulk_stored_energy_j()` exposes their storage.
The energy record separately reports bulk storage before/after, local bulk work
and its storage defect. Its constraint-work defect subtracts both bulk work and
pressure dotted with the weak-volume change. Bulk constitutive work, pressure
storage work, weak volume, complete physical energy and cumulative work must
all pass before the new displacement, velocity, pressure and volume history
are published together.

`prepare_with_pressure` takes a compatible P1 initial ENDPOINT bulk pressure.
Original `pressure_potential` reconstructs Jbar, and the initial weak-volume
equations must close against the supplied displacement. `prepare` uses zero
initial pressure and refuses incompatible initial geometry. No initial
displacement, velocity or volume is silently projected. In particular, arbitrary
geometric quadrature volumes are not copied into a mixed state: unresolved
pressure-space components would alternate under the quadratic secant.
Compressible initial velocity is admitted. Subsequent accepted volumes follow
the solved original history and are never reset to geometric determinants or
reprojected to endpoint pressure. Interval pressure thereafter remains an
algorithmic work multiplier, not endpoint Cauchy pressure.

G0/G1 tests in `tests/mixed_conservative_bulk.rs` passed all six native
cases, alongside all four full-material and three exact-constraint regressions.
Independent original-energy, rational consistent-mass, momentum and weak-volume
calculations cover freely deforming/spinning bodies, heterogeneous quadrature
moduli, fully prescribed dilation/compression, compatible prestressed
initialization and immutable domain/budget/cancellation refusal. Over 32 free
steps, maximum energy drift was `1.14e-13 J` for quadratic bulk and
`2.09e-13 J` for logarithmic bulk, with nonzero bulk-energy exchange.
The logarithmic trajectory ladder compares displacement, velocity AND retained
volume at 8/16/32 steps against a 128-step reference; observed orders were
1.961 and 2.057. This is bounded time-convergence evidence for the retained
nonlinear history, not a pointwise incompressibility or mesh-independent
stability claim.

`examples/conservative_elastic_motion.rs` now executes FullMaterial,
OriginalMaterial and ExactIncompressible with explicit controls and reports
physical bulk storage separately. All three 32-step native example runs passed.
Contact, material relaxation, moving-support linear exponentials and a
frequency-exact nonlinear finite-rotation split remain separate capabilities.

## Exact affine-force histories on linear solids

`ExactLinearSolid::prepare_affine_forces` prepares the same complete physical
pencil for a fixed clock and an explicit positive characteristic force in N.
`advance_affine(force_start_n, force_end_n, cx)` evolves the force linearly
over each interval. Native `PreparedQuadraticFlow::new_extended_input_controlled`
retains source coordinates `[u,v,f]`; a held force increment enters through
`I/dt`, avoiding a per-step force-slope calculation. These added force
coordinates carry no structural storage or inertia, and no new exponential
implementation or modal approximation is introduced.

Supplied work is independently integrated as `integral f(t) dot v(t) dt`
through the native quadratic observation. It is not inferred from energy change
or replaced by endpoint-average force dotted with displacement. Returned
`mean_velocity_mps` remains the kinematic `Delta u/dt`; only held-force
exchange can use it alone as a work-conjugate flow. The affine source endpoint,
physical interval/cumulative work and cancellation gates all pass before the
physical `[u,v]` state is accepted. Force jumps between intervals are explicit
and change acceleration without introducing a velocity impulse.

The affine complete dimension is `4*N+1 <= 256` (at most 63 free coordinates),
with a distinct explicit cold-memory estimate and the same multiplication cap.
The held route retains its original capacity. Wrong prepared mode, malformed
forces, unrepresentable finite force increments, source error, exhausted budgets
and cancellation refuse without changing accepted state. Exactness applies to
the declared affine force history; approximating another continuous history
still requires temporal refinement.

G0/G1 tests in `tests/solid_linear_force_history.rs` passed all four cases.
Coupled non-diagonal pencils are checked against analytic oscillator responses,
with work evaluated independently by integration by parts. Zero-stiffness
constant-jerk motion uses exact cubic/quadratic trajectories and a polynomial
work integral. One interval and four subintervals agree through continuous
ramps and finite force jumps. Affine-specific dimension/overflow boundaries,
wrong mode, cancellation and deterministic replay are covered. All five
existing held-force tests, including the million-step phase case, remain green.

## Exact coupled viscosity on linear solids

`dynamics::ExactLinearSolid::prepare_viscous` and
`prepare_viscous_affine_forces` retain the original complete free-coordinate
stiffness K and consistent mass M while admitting an explicit symmetric
positive-semidefinite viscosity C [N s/m] in the SAME free-DOF order. They solve
M vdot + C v + K u = f through the existing fixed-clock quadratic exponential.
C can include non-proportional coupling; no modal diagonalization or Rayleigh
assumption replaces its entries. This is a memoryless constitutive resistance
or boundary impedance, not the hereditary tensor-Prony material owner.

The original-dimensional C is checked for symmetry and negative eigenvalues
against the caller's explicit relative operator allowance, normalized by C's
own magnitude with no unit-sized floor. No eigenvalues are clipped. Dynamics
use R_vv = M^-1 C M^-1, retaining Q = diag(K,M) and the original physical
energy. An independent quadratic observation integrates D = integral v^T C v dt.
Held-force work remains f dot Delta u; affine work remains integral f(t) dot v(t) dt.
Interval and compensated cumulative balances enforce Delta H + D = W before
state publication. Negative integrated loss beyond the stated energy allowance
refuses. Energy reports expose interval and cumulative physical dissipation.

Zero viscosity retains the conservative trajectory and existing fixed-clock
semantics. Observation storage/work are included in the same explicit preparation
budgets. Clock, forces, overflow, cancellation and work/loss refusals leave all
accepted physical state, time and cumulative history unchanged. Stationary
homogeneous supports and the existing dense augmented-dimension limits remain;
this adds neither support reactions nor moving-support forcing.

Native `tests/solid_linear_viscous.rs` passes five focused tests: underdamped,
critical and overdamped analytic motion with independently integrated loss;
rigid drag and terminal velocity; coupled consistent M with noncommuting K/C;
affine force work and heat against continuous integrals and partitioned motion;
and invalid/indefinite viscosity, resource limits, zero-C equivalence and
refused/cancelled-step rollback. The five held-force and four affine-force
regressions also pass, including the existing million-step phase check.

## Complete-spectrum linear solid motion

`dynamics::spectral::CompleteSpectralSolid::prepare_with_basis` admits all n
generalized eigenpairs of the original n-coordinate TetElasticAssembly. The
n*n basis is row-major with modes in columns. Every supplied eigenvalue must
be finite and nonnegative and is retained verbatim. Missing modes, duplicate
columns, invalid source matrices or incompatible dimensions refuse. This is a
complete coordinate transformation; no spectral window, modal truncation,
tail substitute, diagonal mass or fitted stiffness enters the motion.

Preparation recomputes K phi - lambda M phi against the original matrices,
using the caller's explicit absolute eigen-residual allowance plus a relative
scale from that mode's actual K phi and lambda M phi. There is no unit-sized
floor. The full row-sum norm of Phi^T M Phi - I must satisfy a declared
tolerance strictly below one; initial u/v reconstruction has a separate
relative allowance. Symmetry is admitted and rounded as in the dense exact
owner. These are numerical residual checks, not interval-certified spectral
enclosures or experimental qualification.

Each fixed-clock mode uses exact oscillator/affine-source coefficients with
stable entire-function series near zero. Zero-frequency modes retain the
exact cubic displacement and quadratic velocity under force ramps without
inverse stiffness. Physical force projection uses Phi^T; all displacement
and velocity coordinates are reconstructed. The independent original-coordinate
work is f_start dot Delta u + Delta f dot (u_after - mean(u)), obtained from
integration by parts along the continuous trajectory, not from the energy
difference or endpoint-average force. Original K/M energy, interval and
compensated cumulative work are checked before accepting modal state, physical
state, time or history. Force jumps between intervals apply no impulse.

The 3n+1/4n+1 augmented-exponential limit does not apply to this owner. Explicit
admission caps cover free coordinates, simultaneous preparation bytes,
preparation and interval scalar products, and accepted steps. Preparation is
O(n^3); retained storage and warm intervals are O(n^2). Cancellation is checked
at rows and vector-sized arithmetic tiles; capacity, overflow, energy and
cancellation refusals preserve the accepted state. This conservative path has
zero physical dissipation and does not silently diagonalize a coupled viscous
operator. General coupled viscosity retains the separate dense exact owner.

The caller supplies the complete basis. The existing
`fs_modal::eigh_gen_dense` can produce all n modes and is used explicitly by
`examples/spectral_elastic_motion.rs`; its small, bounded cold eigensolve has
no internal cancellation polls. No automatic unpolled eigensolve is hidden in
the spectral owner, and no sparse million-coordinate scalability claim is made.

Native verification passes four `tests/solid_spectral_dynamics.rs` tests:
128 coupled coordinates against continuous analytical affine motion and work;
rigid and tiny-positive modes against an independent ODE power series and
partitioned ramps; invalid/incomplete bases including a 1e-30 stiffness
falsifier; and exact capacity boundaries, cancellation, energy refusal and
retry. `tests/solid_spectral_tet10.rs` runs actual 12-cell TET10 assembly,
retains signed consistent mass, obtains all 108 native modes, exercises the
lowest and highest modes, and independently checks original sparse K/M energy
and force work. Measured free-motion energy error is 1.345e-20 J in that
fixture. The example compiles and runs a synthetic nodal pulse/release through
the same production geometry, assembly, spectrum and dynamics. This is
represented-pencil temporal evidence, not mesh convergence or measured-material
validation. The previous nine held/affine and five viscosity regressions pass
alongside these new tests.


## Exact cubic force histories on the complete elastic spectrum

`CompleteSpectralSolid::prepare_cubic_forces_with_basis` explicitly prepares
piecewise cubic physical forces while retaining every generalized mode and the
original consistent mass and stiffness. `advance_cubic([F0, F1, F2, F3], cx)`
uses `F(s) = F0 + F1*s + F2*s^2 + F3*s^3` for normalized interval time
`s = t/dt`; every coefficient has units of newtons, not newtons per second.
Initial state, mode admission, the fixed time grid and force coordinates retain
the complete-spectrum contract above. No mode truncation, mass lumping, physical
damping, or nonlinear stiffness is introduced by this owner.

Cold preparation caches the entire oscillator functions through E8. A compensated
series handles zero and small phase, including exact rigid polynomial motion;
trigonometric recurrence handles larger phase. The fixed-clock step projects the
four force coefficients, advances all modal positions and velocities, and
reconstructs physical endpoint acceleration and three normalized displacement
moments `Ik = integral_0^1 s^k*u(dt*s) ds`, k = 0, 1, 2, in metres. These moments
allow support work and other polynomial coupling observations to be integrated
without subdividing the physical interval or evaluating extra force samples.

The applied work is independently evaluated as
`F0 dot (u1-u0) + sum_(r=1..3) Fr dot (u1-r*I_(r-1))`. It is compared with
original-coordinate strain plus kinetic energy both per interval and
cumulatively. The accepted state, moments, work accumulation and clock publish
only after all finite-value, work-budget and cancellation checks pass. Internal
read-only trial evaluation permits an enclosing physical coupling to impose its
own full-system gates before the same infallible state publication.

The explicit cubic capacities are `64*n^2 + 2048*n` bytes,
`4*n^3 + 32*n^2 + 4096*n` cold multiplication terms and
`24*n^2 + 1024*n` warm multiplication terms. Existing held and affine preparation
capacities and their arithmetic are unchanged. Source coefficients, candidate
moments and acceleration are finite-checked; refusal preserves accepted history.

`solid_spectral_polynomial` checks the zero/small/series-switch/multiple-period
phase regimes against independent continuous solutions, all 96 coupled
coordinates with a non-diagonal mass, independent physical work, affine
agreement, interval partitioning of all three moments, exact capacity boundaries,
and invalid, overflowing, cancelled and energy-refused candidate rollback.
Focused native source-closure validation also reruns the held, affine, viscous and
complete-spectrum tests. This is floating-point exact linear propagation subject
to admitted source/basis and balance tolerances, not an interval certificate or a
finite-rotation nonlinear dynamics claim.


## Moving supports with original consistent inertia and reaction work

`dynamics::supported::SupportedSpectralSolid` advances the conservative linear
solid with explicitly prescribed support motion. Supply an original assembly
that still contains the prescribed rows and columns. Support labels must be a
nonempty, strictly increasing proper subset of the assembly's represented global
DOF labels; at least one free coordinate is required. Initial displacement,
velocity and applied-force endpoints use the original represented order; support
accelerations and reactions use the supplied support order. Coordinates already
omitted by assembly remain homogeneous stationary supports, and their omitted
reaction rows cannot be recovered by this owner.

`partition_free_pencil` extracts the entire free/free K and M blocks for a native
complete generalized eigensolve. `prepare_with_basis` requires all free modes,
checks the whole source mass for positive definiteness and the whole source
stiffness for numerical positive semidefiniteness, and admits the same source
symmetry allowance as the complete-spectrum owner. A valid free block alone
cannot hide indefinite support inertia or a zero-diagonal indefinite stiffness
coupling. The full K check uses complete diagonal pivoting and the infinity norm
of the entire unresolved Schur residual; no eigenvalues are clipped and retained
operators are not replaced by the admission factors. These floating-point
operator checks are not interval-certified positivity claims.

Each `advance_affine_supports` linearly interpolates retained support acceleration
to the supplied endpoint and integrates it to continuous support velocity and
cubic displacement. It linearly interpolates the full physical applied force.
The free equation uses
`M_ff*a_f + K_ff*u_f = f_f - K_fs*u_s - M_fs*a_s`, retaining the ORIGINAL
consistent cross mass. This effective cubic force is advanced through every free
mode by the exact polynomial spectral kernel. There is no mass lumping, omitted
mode, velocity impulse or repeated dense eigensolve in an interval.

Reported endpoint reactions are the forces exerted BY prescribed supports ON the
solid, `r_s = (M*a + K*u - f)_s`. The original free-row equation residual must fit
its explicit absolute newton allowance. Strain and kinetic energies use full
original-coordinate K/M, including free/support cross terms. Applied work is
integrated from the continuous trajectory. Reaction work is separately integrated
using the exact displacement moments and twice integrating its inertial term by
parts; it is never inferred from the energy difference. Applied loads on prescribed
rows correctly change the reaction and the two work accounts without changing the
prescribed trajectory. The physical closure is
`Delta H = W_applied + W_support` both per interval and cumulatively.

The inner free candidate, original-row residual, independent full work balance,
finite values and physical endpoint energy all pass before either owner's accepted
state changes. Materially negative original strain energy is refused even if a
numerically tolerated stiffness residual produces matching negative support work.
The final cancellation checkpoint precedes infallible joint publication of modal
state, full motion, support acceleration, compensated work histories and clock.
No prepared matrix is cloned during an interval.

Separate whole-system caps admit n represented coordinates, preparation storage
`192*n^2 + 8192*n + 8*element_count` bytes, cold multiplication work
`12*n^3 + 128*n^2 + 16384*n` and warm multiplication work
`64*n^2 + 8192*n`. The inner cubic spectral controls also apply. Preparation is
O(n^3), while retained matrices and intervals are O(n^2). Cancellation is polled
through preparation, positivity factors and interval rows. The basis remains
caller-supplied; the example's small native dense eigensolve is an explicitly
bounded cold operation without internal cancellation polls.

Five `solid_supported_dynamics` regressions check an independently solved
moving-base oscillator with non-diagonal consistent mass, separately integrated
applied and actuator powers, uniform accelerating translation, partition invariance,
constrained-load invariance, original full-pencil rejection, capacities, invalid
sources, cancellation, energy refusal and bitwise rollback. The independent
oscillator's cumulative balance defect is below 3.886e-16 J.
`examples/supported_elastic_motion.rs` assembles an actual 30-coordinate TET10
specimen, prescribes 18 base-face coordinates and retains all 12 free modes. Two
opposite triangular acceleration pulses move the base 62.5 micrometres and stop it;
480 intervals at 48 kHz retain independently observed actuator work and original
physical energy within 9.672e-17 J. Geometry and material are synthetic. This is a
native assembly/solve/motion consumer, not spatial convergence or experimental
calibration evidence. Viscosity, hereditary loss, finite rotations, changing
supports/topology and impulsive support motion are outside this conservative owner.

## Physical boundary loading of hereditary harmonic solids

`ViscoElements::prepare_harmonic_response_with_port_impedance` assembles the
original material spectrum and consistent mass together with a caller-supplied
physical load `Z(omega)` and real power-conjugate force map `B`. The complete
pencil is `K' - omega^2 M - i K'' - i omega B Z B^T`. Reactive and mutual
entries remain present; a loaded resonance is solved before any bare-system
inversion. One factor is reused at that frequency.

The shared harmonic owner supplies pointwise Hermitian-PSD admission and
separate port/projection/payload limits. The load owner retains responsibility
for its frequency domain, units, causality and source accuracy. Returned
`port_load` contains the original flow, opposing effort and net absorbed
boundary work. This work is separate from the original Prony/Kelvin material
heat and is not relabeled as radiation or heat. The signed balance subtracts
both contributions. Original source-band and full-cycle strain admission,
structural payload budgets and cancellation still apply.

`tests/visco_harmonic.rs` checks a load-regularized bare resonance, an independent
coupled spring/mass/resistance boundary on a hereditary solid, separate work and
heat, exact zero-load response parity, and material/load/capacity/cancellation
refusals (G0/G1/G3/G4). These bounded numerical checks do not establish broadband
load passivity, a spatial convergence certificate or measured material accuracy.

### Original hereditary free-decay poles

`ViscoElements::prepare_resonance_poles` prepares the complete homogeneous
fixed-support spatial pencil in `[free displacement; original Maxwell strain]`
coordinates. Every six-component Mandel strain at every original quadrature
point and branch remains present, including equal-time and zero-energy
directions. The existing consistent mass, equilibrium stiffness and complete
Kelvin tensor are retained. Under `exp(s*t)`, the original equations are
`(K_eq+s*C_kv+s^2*M)u + sum(V*B^T*C_j*q_j)=0` and
`(1+s*tau_j)q_j-s*tau_j*B*u=0`. Keeping the branch coordinates makes this
pencil polynomial through material relaxation rates; eliminating them would
introduce poles into the operator and invalidate holomorphic contour counts.

Three real coefficient matrices and fixed two-sided scales are cached once.
Subsequent operator/analytic-derivative evaluations cost O(N^2), with all dense
factorization and contour iterations delegated to the existing `fs-la`
Newton/Beyn owners. Source coordinate, numerical payload, preparation-work and
per-query evaluation limits are admitted explicitly. The source budgets are
separate from the shared solver's iteration/storage/work limits.

Returned modes retain original coordinates and separately check mechanical
force and material-strain residual families against their original coefficient
envelopes. Negative real Laplace rate means decay; a purely real pole is
relaxation, and an unobservable memory state is not an audible resonance.
Eigenvectors have arbitrary common amplitude and are not admitted finite-strain
histories. The original real-frequency material-band intersection is retained
as a conservative `|s|` gate, including the complete contour region. Analytic
continuation of the declared Prony model is not experimental validation at
complex frequency or a material-data extrapolation certificate.

Local refinement needs a supplied structural seed and refuses initialization
exactly at an uncoupled material relaxation rate. Contour discovery keeps the
shared numerical conditioning, nested-quadrature, count/rank and simple-root
admissions. Two-moment Beyn cannot resolve more independent roots than its
probe rank; a quadratic augmented pencil may therefore require contour
subdivision even when every retained coordinate fits the resource budget.
Unresolved multiplicities/defective roots refuse. Numerical residual and
quadrature agreement do not certify spatial convergence, spectral completeness
or rigorous damping enclosures.

`tests/visco_poles.rs` checks an independently factorized scalar Prony cubic,
its oscillatory pair and real relaxation pole, a live contour quadrature ladder,
coupled tensor/Kelvin/consistent-mass roots against an independent determinant,
source/resource/cancellation refusals, lossless roots and full-spectrum rank
refusal (G0/G1/G3).

### Physical boundary stiffness, damping and inertia in hereditary solid modes

`ViscoElements::prepare_resonance_poles_with_port_load` adds a complete real
power-conjugate port map B and symmetric numerical-PSD boundary coefficients
K_port, C_port and M_port to the original hereditary spatial pencil. The
opposing effort is `K_port*d+C_port*d'+M_port*d''`, where `d=B^T*u`; the
body receives `-B*effort`. The mechanical pencil therefore receives
`B*(K_port+s*C_port+s^2*M_port)*B^T` before balancing. Its causal impedance
is `K_port/s+C_port+s*M_port`. No single-frequency impedance sample is
continued, and no boundary mass, mutual entry or source Maxwell coordinate is
discarded. This represents kinematically attached mass, springs and dampers;
an attachment with independent moving coordinates requires its own state.

The caller supplies the complete original free-coordinate force map, all
three port matrices, a declared real-frequency band, and additional port,
payload and preparation-work limits. Each coefficient is checked for exact
symmetry, finite values and numerical positive semidefiniteness under a
diagonal congruence, reusing the bounded shared fs-la Jacobi eigensolver.
Negative diagonal entries and nonzero rows with zero diagonals refuse before
numerical tolerance could conceal an active direction. The original matrices
are retained unchanged; singular passive directions are allowed. These
floating-point checks are not interval-certified passivity statements.

The load's additional payload envelope is `24*n*p+80*p^2+256*p` bytes. Its
preparation-work envelope is `3*(720*p^3+256*p^2)+6*n*p*(n+p)` scalar-work
units; every individual PSD admission also obeys fs-la's unpolled-work cap.
These limits are separate from the original source-pencil and Newton/Beyn
budgets. Port projection is prepared once; repeated pole and contour queries
reuse the original quadratic coefficients and fixed scales. Cancellation is
polled before/after bounded PSD calls and at projection/observation rows.

Each returned pole includes optional original port displacement, flow and
opposing effort, load storage `(d^H*K_port*d+flow^H*M_port*flow)/2`, boundary
heat `Re(flow^H*C_port*flow)`, net absorbed power `Re(flow^H*effort)`, and the
signed defect `power-2*Re(s)*storage-heat`. These are complex-envelope forms
at the mode's arbitrary normalization, not finite-cycle averages or admitted
finite-amplitude trajectories. Boundary storage and loss remain separate
from intrinsic material loss. The original mechanical and memory residuals
include the added load. The intersection of original material and load bands
gates `|s|` and the whole contour; analytic continuation is not measured-data
validation at complex frequency. Existing numerical-root, multiplicity,
completeness and spatial-accuracy limits remain unchanged.

`tests/visco_loaded_poles.rs` checks an independently factorized loaded Prony
cubic, original material memory ratios, independent boundary power/storage/heat,
two mutually coupled ports against an independently assembled rational
determinant with consistent cross mass, exact zero-load pole/vector parity,
singular passive directions, source-band admission, active and malformed
loads, exact resource boundaries, and cancellation (G0/G1/G3/G4). The full
native fs-solid library builds, all four loaded-pole tests pass, and the five
existing hereditary-pole regressions pass against the updated library.

## Full hereditary solid poles with independent dynamic physical loads

`ViscoElements::prepare_resonance_poles_with_dynamic_load` assembles the complete
original body/material/load pencil. The supplied `ViscoPoleDynamicLoad` retains
every coordinate of an admitted exact quadratic `PortHamiltonian`, including
mechanically dark states, full mutual storage/resistance and all original port
columns. For body displacement `u`, the same signed physical map `B` gives load
flow `v=B^T*u'` and opposing nodal force `-B*effort`. With the original load law
`z'=(J-R)Qz+Gv`, `effort=G^T*Qz`, the added equations are
`(sI-(J-R)Q)z-s*G*B^T*u=0` and the body receives `B*G^T*Qz`. No standalone load
resolvent, sampled impedance, material denominator, projection or modal damping
ratio replaces those rows. The joint polynomial is regular to evaluate at bare
load and Maxwell relaxation poles. Its bounded dense coefficients and fixed
scales are prepared once, while all root iterations remain with `fs-la`.

The existing complete-state dimension/payload/evaluation/solver budgets include
the appended load states. `ViscoPoleDynamicBudget` separately admits original
load state/port counts, payload, projection work and the shared numerical PSD
check. Exact finite quadratic storage, symmetric numerically PSD Q/R and the
original admitted pHS structure are required; nonlinear storage is refused.
Source and load bands are intersected and applied conservatively to `|s|` and
the entire contour. This is a declared causal continuation, not experimental
evidence outside measured frequencies. Preparation and observation poll `Cx`,
and no accepted material/attachment history is read, advanced or reset.

`physical_vector` is original free displacement, all original six-component
Mandel branch states, then the independent load's original state order.
`dynamic_load_range()` exposes the suffix. Each pole must pass separate
original coefficient-scaled mechanical, Maxwell and load-state residuals.
`resonance_pole_from_state` accepts a full authored seed, including dark load
or memory modes, without a displacement-only resolvent initialization.
The convenience displacement seed initializes load coordinates to zero and
retains its existing exact Maxwell initialization/refusal behavior.

The dynamic-load observation independently reports `H=z^H Qz/2`,
`D=(Qz)^H R(Qz)`, incoming port power, and `P-2 Re(s)H-D`. Body equilibrium,
consistent-inertia and every Maxwell storage term are evaluated separately,
as are Kelvin and Maxwell dissipation. The joint signed defect is
`2 Re(s)*(H_body+H_load)+D_body+D_load`. These are complex-envelope quantities
at the returned arbitrary mode amplitude, not cycle averages or fitted heat.
Signed numerical defects remain visible. Finite pencil roots need not be
observable mechanical resonances. Local/contour results remain numerical,
without interval enclosures, global completeness, nonlinear-load support,
measured material validation or a new spatial-convergence certificate.

G1/G3/G4 specimens in `visco_dynamic_poles.rs` verify analytical attached-mass
mode splitting, a seed exactly at the bare attachment pole, 128/256-node contour
agreement and honest 64-node count refusal, hereditary loaded roots against an
independent scalar eliminant, separate original body/load energy and heat,
signed mutual ports, mechanically dark original load states, nonlinear-law
refusal, resource/band refusal and cancellation. Native fallback verification
compiles the complete `fs-solid` crate against aligned cached dependencies;
it does not claim a DSR or full-workspace run.

### Forced response of the complete hereditary solid and dynamic load

`PreparedViscoPoles::harmonic_factor` factors its original conditioned
polynomial at `s = -i omega`, retaining every original free displacement,
six-component Maxwell strain state and independent quadratic pHS load state.
It accepts positive finite frequencies in the intersection of the original
material and attached load bands, including regular coupled problems at a bare
host resonance or bare attachment resonance. It never evaluates an isolated
load inverse, condenses away a load/material state, fits a replacement damping
law, advances accepted histories or introduces a homogeneous state.
`ViscoJointHarmonicFactor::solve` reuses that immutable shared `fs-la` complex LU
for multiple forcing columns and phases. Original mechanical, memory and load
rows are independently checked with the prepared physical residual tolerance;
up to three corrective solves reuse the same factor without changing that
tolerance. Genuine joint singularities and unrepresentable states refuse.

Force and state phasors are peak amplitudes for `exp(-i omega t)`. Every original
material quadrature point undergoes the same full-cycle strain-envelope
admission as `ViscoHarmonicOperator`; separate real/imaginary amplitude checks
cannot admit an excessive intermediate phase. The returned body storage uses
the original equilibrium stiffness, consistent mass and Maxwell strain energy;
material heat uses the original Kelvin tensors and every Maxwell relaxation
branch once. Static K/C/M and independent dynamic loads retain separately
observed original storage and dissipation. Body energy is never obtained by
subtracting a much larger attachment energy from a combined account.

All reported energies, heat rates and powers are cycle means. Applied work is
`Re(force^H velocity)/2`; absorbed interface work is
`Re(flow^H effort)/2`, positive into the load. Signed body, total and individual
load power defects remain visible. Power admission uses the unchanged physical
tolerance and the sum of apparent applied/interface work, actual dissipation
magnitudes and `2 omega` times original mean storage magnitudes. This physical
oscillatory work scale admits exact antiresonances with zero interface flow and
finite internal load motion; it adds no constant floor, clips no signed power,
and does not classify reactive storage as heat. Passivity restrictions remain
those of the original source and load preparation.

`ViscoJointHarmonicBudget` bounds complete coordinate count, numerical payload
and conservative factor/solve/refinement/observation work, separately from the
already admitted immutable source cache. Polling occurs during assembly, shared
LU/solves, original-row checks, quadrature admission and observation. Refusals
leave both the factor and accepted source histories unchanged. Determinism is
fixed-order same-build/ISA behavior; numerical budgets exclude allocator
metadata and process RSS. This is a dense bounded response path, not a sparse
large-model performance claim, an interval enclosure, a nonlinear harmonic
balance method, a zero-frequency solver or support for descriptor/nonquadratic
independent loads.

Focused G1/G3/G4 evidence in `tests/visco_joint_harmonic.rs` covers analytic
coupled-oscillator states at bare attachment antiresonance and bare host
resonance; full Maxwell/Kelvin and mutual load dissipation across five forced
frequencies; resolving body energy beside an attachment more than ten orders
larger; phase/repeated-RHS reuse; original full-cycle strain refusal; singular,
invalid-band, resource and cancellation refusal; and subsequent factor reuse.
Native compilation of the complete `fs-solid` crate and the four new tests,
nine existing `visco_harmonic` tests and four existing `visco_dynamic_poles`
tests passed. These are focused native checks, not a workspace or DSR gate.


### Original constrained dynamic attachments

`prepare_resonance_poles_with_descriptor_load` accepts original quadratic `DescriptorPortHamiltonian` loads through `ViscoPoleDescriptorLoad`. It retains `E=diag(I_diff,0)`, effort `e=[Q_diff*z_diff,z_alg]`, every differential state, every zero-storage multiplier, and every original Maxwell strain. Velocity input `B^T*u'` and opposing force `-B*G^T*e` preserve power conjugacy. No descriptor matrix, isolated load impedance or constraint is inverted or regularized.

Finite poles and forced harmonic response reuse the existing prepared pencil and shared numerical owners. Original state ranges distinguish differential states and reactions. Multipliers contribute to port power and full mutual resistance heat but never acquire stored energy. Every original algebraic row independently passes the physical residual gate, including same-factor harmonic refinement. An all-differential descriptor preserves ODE arithmetic and state order. Full source bands, state/port/payload/work limits and cancellation apply. Exact quadratic storage is required; descriptor index and infinite-spectrum completeness are not certified; redundant or singular constraints remain explicit refusals.

Native full `fs-solid` compilation and all 12 focused regressions passed: four descriptor tests, four original joint harmonic tests and four original dynamic-pole tests, using the actual source and cached production dependency artifacts. This focused validation is not a full workspace Cargo gate or a descriptor-index/infinite-spectrum certificate.


### Original continuous state-space boundaries in hereditary solid resonance

`ViscoPoleStateSpaceLoad` and `prepare_resonance_poles_with_state_space_load` retain the complete original real law `z'=A*z+B*v`, `effort=C*z+D*v+E*v'` with `v=force_map^T*u'` and opposing body force `-force_map*effort`. Full mutual D and derivative E terms enter the existing first- and second-degree mechanical coefficients. Every original state, including conjugate-pair coordinates and mechanically dark states, remains in the full polynomial; zero-state static D/E boundaries are supported. No isolated load impedance, source resolvent or derivative/direct matrix is inverted, fitted or condensed. Shared Newton/Beyn and reusable complex LU remain the only numerical solver owners.

The caller owns the original realization, state/port units, real-frequency support and stability/passivity evidence. Source and boundary bands are intersected and applied conservatively to |s|. Complete state/port/payload/preparation budgets and cancellation apply before copying the load. Each original boundary-state equation is checked independently within the physical residual/refinement gate. Original body Maxwell/Kelvin storage and heat remain separate from the boundary.

A generic state-space realization supplies no physical storage/dissipation witness. Pole and forced observations therefore expose original states, flows, efforts and signed per-port/total absorbed work, including negative mutual exchanges. They do not label that work as heat or invent boundary energy. Forced `average_total_storage_j` and `total_power_balance_residual_w` are now `Option<f64>`: the existing complete physical pHS/KCM paths return `Some`; an unqualified state-space boundary returns `None`. The original body energy/heat and body-plus-boundary-work balance remain independently computed and enforced. Numerical roots do not certify source passivity, experimental complex-frequency support or spectral/spatial completeness.

G1/G3/G4 regressions in `visco_state_space_loads` cover bare oscillator antiresonance with original derivative inertance, exact loaded free modes, signed mutual D/E work beside all Maxwell/Kelvin strains, original dark-state finite poles, finite-amplitude admission, resources, support, cancellation and factor reuse. Existing fully physical response tests explicitly unwrap complete energy accounts.

### Instantaneous momentum in complete linear-solid trajectories

`ExactLinearSolid::apply_momentum_impulse` and
`CompleteSpectralSolid::apply_momentum_impulse` accept resolved impulses in
the original free-DOF order. The dense owner retains the inverse already
computed and admitted for its source map; the complete-spectrum owner uses
`Phi^T J` and all basis columns, including zero-stiffness rigid modes. Each
event uses O(n²) operations and O(n) scratch, without a new factorization or
propagator rebuild. The dense owner retains one additional n² matrix within
its conservative preparation envelope. These are operation/storage bounds,
not measured performance results.

Both owners audit `M * (v_after - v_before)` using the actual represented
endpoint and an explicit absolute per-component impulse allowance. An admitted
approximate basis cannot bypass this physical audit. Signed work is
`J dot (v_before + v_after)/2`; its independent residual uses the momentum
residual directly, with compensated product accounting. Work enters the
original cumulative account once. Displacement, clock, accepted-step count,
viscosity/source preparation and earlier physical losses remain intact.
Zero impulses preserve physical/modal bits; cancellation, unrepresentable
kicks and failed momentum/work gates publish nothing.

These reduced pencils cannot recover discarded support reactions or claim a
full-body angular inventory from an incomplete DOF map. The API applies a
resolved impulse; it does not select contact laws, restitution or event times.
G1/G3/G4/G5 regression definitions cover off-diagonal mass, viscous/affine
continuation, a free tetrahedral mass's independent linear/angular balance,
mixed rigid/elastic spectral response and analytic continuation, zero identity,
rounded-basis work refusal and cancellation. Native execution remains outstanding.

### Impulses and instantaneous support velocity changes

`SupportedSpectralSolid::apply_momentum_impulse` accepts applied impulses
in the original represented-coordinate order and the prescribed support
velocities immediately after the event. Positions, clock, accepted-step count,
prepared operators and retained support acceleration remain unchanged.
Subsequent continuous intervals start from the new support and free velocities.

The original consistent mass supplies
`J_eff = J_f - M_fc * Delta v_c` to the complete free spectral owner.
Support reactions are `R_c = M_cf * Delta v_f + M_cc * Delta v_c - J_c`.
Both cross blocks remain active; applied impulses on support rows change the
actuator reaction without changing the prescribed motion. The actual represented
full velocity is independently checked against the original mass in every row,
including whether large opposing applied/support impulses can represent their
small required net transfer. A failed outer audit also discards the unpublished
inner spectral candidate.

Applied and support work are separately computed from their impulses dotted
with average pre/post velocity. Their signed cumulative accounts continue through
later native intervals. Compensated product accounting and direct momentum
increments avoid subtracting large momentum totals. The existing rounded
total-energy reference convention remains unchanged, so an exact zero event
does not invent a debt when one store is below the ULP of another.
The event borrows all prepared matrices and uses O(n²) operations with O(n)
scratch. These are algorithmic bounds, not measured performance.

G1/G3/G4 regression definitions cover independent consistent-mass reactions,
actuator work, analytical post-event continuation, uniform boosts, support-load
partition invariance, zero identity across separated energy scales, and late
full-row refusal/cancellation followed by an unchanged native continuation.
This API prescribes resolved impulses/support velocities; it does not choose
contact, restitution, event timing or a compliant actuator law. Native test
execution remains outstanding.
