# fs-solid CONTRACT

## Strain-life from the original cyclic J2 constitutive path

`plastic3::fatigue::J2FatiguePath` owns the latest native solution and one
shared `StrainLifeMiner` per fixed-plane gauge while borrowing the original
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

The CMB model and named uncorrected mean-stress policy retain their empirical,
fixed-plane scope. No stiffness, heat, fracture, failure topology, ratcheting
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

`linear3::visco` (IRP M1.3-W2a, in progress) prepares P1 or TET10 elements
with caller-supplied, mesh-frame `fs-material::TensorProny` laws and density.
It reuses the existing elastic mesh admission, equilibrium stiffness and
consistent mass, and the shared Mandel B operator. Complex stiffness retains
the full storage and loss tensors, including off-diagonal couplings. Every
frequency query must lie in every element's declared law band. Empty branches
reproduce the canonical elastic stiffness bitwise. Positive quadrature weights
propagate pointwise relative quadratic-form bounds to the assembled forms;
this alone is not a frequency-dependent modal-eigenvalue certificate.

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
certified modal block reduction, valve-card integration, HEX histories,
displacement time integration, or material/specimen validation. Arbitrary
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
Force-only advancement, finite force jumps and reduced-only state adjoints
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

`tests/transient_supports3.rs` supplies G1 exact mass/stiffness/both-Rayleigh
cross-term initialization, a closed-form accelerating-base oscillator with
second-order refinement, fully prescribed thermal/external/support work and
reaction checks, and G3/G4 fork, malformed-input and cancellation checks.
`tests/transient_support_history3.rs` covers rigid translation with damping,
agreement with explicit core endpoints, physical frame time beyond a clamped
table, and transactional sampling, binding and cancellation refusals.

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
  exactly zero strain — energy-invariance battery). Statics: FD
  residual/tangent Newton with residual backtracking and load
  stepping, fixture-gated; fixture slenderness matters — a
  penalty-stiff EA puts a quadratic spurious-stretch wall around
  every Newton iterate (measured, documented in the battery).
  Analytic tangents and SE(3) DYNAMICS under fs-time symplectic
  integrators are recorded successor scope.
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
  Numerical lag rounds refresh normal loads and tangent frames only,
  and success requires the refreshed friction residual to meet tolerance.
  The frozen-frame Huber energy, gradient, and Hessian agree in both
  sticking and sliding branches. When sliding exposes an indefinite
  tangent, a stiffness-scaled diagonal shift produces a descent search
  direction while preserving the physical force law, tangent, pin
  equations, and final residual gate. `eps_v` is a displacement smoothing
  width; this regularized law permits small creep and does not claim exact
  set-valued static friction. Solutions expose separate normal and
  tangential reactions plus nonnegative endpoint incremental friction work.
  The latter is not trajectory-integrated heat; thermal deposition and
  rigorous continuous-trajectory CCD remain separate capabilities.
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

No-claims (`plastic3`): small strain and rotation; rate independence;
constant-strain tetrahedra lock under confined plastic incompressibility
(confined limit loads are over-predicted; no B-bar/mixed formulation); dead
loads; no limit-load continuation; no mesh-convergence claim.

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
structured invalid-input returns.

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
Conditions remain fixed for a damage history. Source refresh across thermal or
chemical evolution, adaptive fracture topology and complete structural/thermal
energy accounting are not introduced. The `cohesive_source` integration tests
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

Every advance explicitly admits the unchanged interface query, including quantity
conventions, so evolving bulk temperature cannot silently freeze an unsupported
cohesive-state change. Three source integration tests cover independent series and
thermal-work references, data-only substitutions with relaxation and fracture
continuation, and source/domain/cancellation refusal. This remains a quasistatic
fixed-geometry series model with prescribed bulk temperatures and fixed external
interface state; no thermal transport, fracture-heat partition or
temperature-dependent cohesive-history transfer is inferred. New modules pass
formatting and whitespace checks. Native tests remain unexecuted: a single-threaded
metadata-only check was SIGKILLed in unchanged asupersync in the authoring worker.

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
