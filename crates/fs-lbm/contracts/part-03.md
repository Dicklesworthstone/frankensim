<!-- fs-lbm contract, part 3 of 3; a higher number is newer. Index: ../CONTRACT.md. -->

## Physical backward Euler with original tensor conductivity laws

`advance_nonlinear_tensor_conduction` advances original-law full-tensor
conduction with complete fixed sensible capacities C[J/K]. Frozen and nonlinear
drivers share the same time/storage admission and original-state energy-record
owner. Each nonlinear trial retains original cell-to-face conductivity integrals,
endpoint derivatives, source support, finite-contact heat and the immutable
physical old-temperature field. Positive storage admits insulated bodies.
Prescribed-cell reservoirs supply their conductive reaction and storage change
exactly once. Newton, backtracking, shared Krylov and cancellation limits remain
active, and refusal returns no partial interval. Conductivity data does not
infer latent heat, variable capacity or material-history evolution. Physical
time accuracy is first order.

Seven native pipeline regressions cover analytic nonlinear cooling, unequal
capacities with finite contact, temporal refinement, equivalence to frozen
tensor storage, insulated rotated heating, thermostat storage, and invalid
time/storage/support/cancellation. An additional steady regression forces
backtracking inside a narrow original source domain and checks zero-budget
transactional refusal. These pipeline regressions are authored and syntax
checked; the complete dependency graph was not compiled or executed here.
Analytical fixture calculations validate the reference values only.

## Original-law nonlinear full-tensor conduction

`solve_nonlinear_tensor_conduction` accepts one immutable
`TensorConductivityLaw` per voxel. `RotatedSolidConductivity` binds the
original native directional tables to explicitly declared principal axes.
Every cell-to-face heat uses the original conductivity integral across the
whole temperature interval; Newton derivatives use original endpoint tensors,
without a duplicate conductivity derivative. Fixed common-factor tensors
recover the scalar Kirchhoff transformation of native HMM; grid-diagonal
laws recover each integrated directional half-cell. Arbitrary full K(T) is a
nonsymmetric columnwise primitive extension with smooth local flux consistency,
not a general discrete entropy, maximum-principle or unconditional-convergence
guarantee.

The native face/contact topology and physical heat ledger are shared with
frozen tensor conduction. Exact local Newton elimination retains only a small
cell-recovery record; the sparse system contains face increments. Original
cell/interface/contact/boundary heat must pass at the unchanged endpoint.
Bounded line search checks positive Kelvin and complete original source support;
one cumulative Krylov budget and explicit Newton/backtracking budgets bound
the endpoint. No trial advances material history or replaces original source
support. The last linear correction residual is reported separately from
original physical heat acceptance.

Evidence: 25 actual-native selected-source Rust tests passed (13 new
constitutive/integral/tangent/condensation checks and 12 existing conductance/HMM
checks). Five public pipeline regressions cover transformed-affine full-tensor
heat, independent rotated laws with finite contacts, constant-law equivalence,
budget/support/cancellation refusal and original-source bounded backtracking. Full pipeline compilation and
execution remain unverified in the current dependency/capacity-limited environment.

## Physical SDIRK2 on the full-tensor conduction operator

`advance_tensor_conduction_sdirk2` provides two-stage, second-order, L-stable
conduction for frozen SPD tensors and fixed per-cell sensible capacities.
Independent stage setups are sampled at c=[gamma,1], gamma=1-1/sqrt(2).
Both native spatial stages retain their actual face temperatures and heat and
share one Krylov iteration budget. Stage 2 retains the original interval
temperature and carries stage 1's native heat through an internal storage RHS,
excluded from external source accounting. Tableau quadrature includes each
prescribed cell's reservoir storage and corrects already-carried heat when the
stage-2 prescribed mask changes. One whole-interval `TransientRecord` reports
physical storage, source and reservoir heat. Smooth forcing and fixed
constraints qualify the temporal order; discontinuities require interval
boundaries. No latent heat, variable capacity, nonlinear tensor law or
maximum-principle guarantee is inferred.

Evidence: five exact-production temporal numerical tests passed with pinned
native rustc. Six native pipeline regressions cover contact decay, temporal
order, rotated-tensor forcing, changing constraints/contacts, thermostat
storage and shared-budget/transactional refusal. Those pipeline regressions
were authored but not executed in this capacity/dependency-limited environment.

## Transient full-tensor conduction with physical cell storage

`advance_tensor_conduction` advances the same conservative hybrid face
operator by one backward-Euler interval. Its explicit per-cell sensible
capacities are in J/K and are counted once, without another volume factor.
The original previous temperature contributes C/dt and (C/dt)*T_old only to
cell equations before static condensation. Faces and zero-thickness contacts
carry no invented storage. Positive cell storage admits insulated bodies.

`TensorConductionStep` retains the actual endpoint face temperatures and
outward heat together with the native `TransientRecord`. Native local
heat-residual checks include storage. The shared physical record accounts
for source heat, external/sink exchange and prescribed-cell reservoir
reactions; a prescribed-temperature jump also receives its C*delta(T) from
that reservoir exactly once. Failed/cancelled trials cannot mutate the
caller's original temperature, capacities or authored setup.

This entry has fixed SPD tensors and fixed sensible capacities over each
interval, and uses the supplied endpoint boundary/source sample. It does not
infer latent heat, varying calorics, advected enthalpy or nonlinear tensor
updates. Backward Euler remains first order in time; a general discrete
maximum principle for rotated anisotropy is not claimed.

Two exact production-kernel native Rust checks passed for statically condensed
unequal-capacity contact equilibration, interface flux/jump and stored-energy
conservation, and convergence to the independent diffusion decay mode.
Four full native pipeline regressions cover those contacts, non-unit-volume
uniform heating with rotated tensors, fixed-cell reservoir storage, and time
refinement; they are defined but unexecuted in this workspace's incomplete
native build.

## Conservative full-tensor conduction on the native box grid

The opt-in `solve_tensor_conduction` owner admits one frozen finite symmetric
positive-definite Cartesian conductivity tensor per voxel on the native graded
grid. A hybrid face-temperature energy retains every tensor component,
reproduces affine fluxes, and reduces to native harmonic half-cell/contact
resistances for diagonal tensors. Local cell condensation leaves shared
perfect-interface traces or two finite-contact traces joined by A/R''.
Temperature, heat-flux, adiabatic and film boundaries, fixed cells, sources,
and linear sinks use their actual physical faces and loads.

`TensorConductionSolution` publishes cell-outward heat and both contact surface
temperatures. Its solid-to-fluid heat observer sums those actual face fluxes;
no scalar conductance is invented for a rotated tensor. Native k(T), source
contacts, eddy additions, compact components and advection remain with their
explicit owners and are refused by this frozen-field entry.

At most eight same-matrix residual refinements share the declared native
Krylov budget. Acceptance checks the recomputed returned-trace algebraic
residual and local cell/interface/contact/boundary heat closure, with native
arithmetic roundoff allowances, positive finite Kelvin endpoints and bounded
cancellation. General rotated anisotropy has no discrete maximum-principle
guarantee. The private cell-only storage seam retains original old
temperatures for a physical transient owner; face traces store no energy.

The Cartesian energy follows Droniou, Eymard, Gallouet and Herbin,
M3AS 20 (2010), doi:10.1142/S0218202510004222, equations 2.17--2.21.
Four exact production-kernel native Rust tests passed on the final local
operator: rotated affine cross-flux, diagonal/contact reduction, condensed
energy/source consistency, and tensor/geometry admission. Full native
pipeline regressions cover graded heterogeneous/contact patches, actual
solid-to-fluid heat, native diagonal parity and unsupported-owner refusals.
Those integration definitions have not run in this incomplete native build.

## Kinetic storage across tetrahedral liquid projection

Tetrahedral liquid ingress evaluates its kinetic change from the original P1
consistent-mass integral and the receiving native cell masses/velocities. Each
face-connected component retains one fixed reference velocity `U`, then reduces
`delta K_rel + U dot delta P_rel + |U|² delta M/2` with double-double products.
An unchanged common translation therefore cannot erase an unresolved strain
or spin mode. The mass and momentum terms remain present: an actual native
readback defect is not removed by changing the arithmetic reference frame.
The original overlap, linear/angular momentum correction and physical model
are unchanged; neither sign of remap energy becomes heat.

The predicted velocity projection, initialized populations and subsequent
thermal/rheology/buoyancy attachment are checked independently against the same
source inventory. Component kinetic changes are kept until the shared sum of
absolute mechanical errors is tested, so disconnected pieces cannot cancel
their remap errors. The existing rounded source/target kinetic totals remain
diagnostics; subtracting them need not reproduce the retained change.
Native readback uses the original tracked interface mass or full-cell
population mass, never density times a rounded/clamped fill fraction. The
same mass supplies the post-attachment component enthalpy audit.

The completed transfer receipt exposes a retained native kinetic increment and
the independently integrated original P1 inventory in a caller-specified fixed
velocity reference. Reframing carries mass and relative momentum as well as
component kinetic storage. Rigid release can therefore compare its original
momentum/inertia storage and compose this transfer stage without reconstructing
a difference from rounded endpoint diagnostics. Thermal attachment refreshes
the native increment from its actual final readback.

The storage kernel uses fixed local arrays and one retained source/candidate
inventory per component; prediction scratch is reused for native readback.
Seven executed production-leaf tests cover P1 variance beneath common motion,
consistent-mass cross terms, actual mass/momentum defects, separate component
errors, reference changes, component merging and invalid inventories. Full
native ingress/mass-owner regression
sources are also present; their execution remains pending an admitted full
dependency build. This is a projection-accounting extension, not a
conservative kinetic remap, subcell energy model or physical heating law.

## Original spin impulses and retained immersed-body kinetic storage

The D3Q19 moving single-body, assembly and finite/shared film-carrier paths
transport the original authored, fluid and film COM spin impulses directly
into the body's momentum. They do not reconstruct those impulses from a
subtraction of rounded angular rates followed by a second inertia product.
The resulting spin drives the original final native free rotation and pose;
fluid collision and the independent film resistance/work laws remain their
existing owners. Source impulses are reduced before their one frame transform.

All three trajectory paths evaluate actual endpoint kinetic increments from
the fixed material/body principal-inertia form and represented velocity and
angular-momentum changes. Compensated products retain small changes beside
unchanged translation/spin. Both free-half-step increments are measured
separately; their numerical rotational defect is removed from the independent
work residual before rounding. No energy residual defines drag, film loss or
heat, and native settling observations consume these same motion receipts.

Four executed actual-production-leaf regressions cover direct small-couple
transport, hidden transverse storage, a rotated full inertia and separate
free-rotation work. Native single/assembly flow regressions additionally cover
accepted spin and fluid acceleration with transverse kinetic storage hidden
by a large common translating inventory; these full integration tests remain
unexecuted because the real dependency build is memory limited. This extends
the existing retained-fluid/dry-correction models, not their physical scope or
finite-step free-rotation energy accuracy claim.

## Variable-Cp phase-volume mixing and bath stability

Source sensible intervals with affine rho(T) and nonconstant Cp(T) retain their
original nonlinear h coordinate in thermal-volume coupling. The selected
single-interval entry uses the same original-chart mixing-volume remainder as
the full-chart entry, with one bounded interval per query and the selected
enthalpy support enforced. The affine-in-h density shortcut applies only to
that actual density law. Full-chart calls retain their aggregate traversal
budget and cancellation checks.

Explicit conductive baths bound dT/dh by 1/min(Cp) over every intersected
source interval, including clipped conductivity support. Authored h/T
intervals retain their existing slope bound. Current source conductivity,
reference inventory, transfer accounts and the explicit mixing-weight
acceptance gate remain the existing owners.

Native regressions exercise a source-Cp volume remainder whose sign differs
from the old density-only shortcut, selected/full-chart equivalence,
calorimetry, traversal limits, cancellation and a variable-Cp bath overshoot.
They are included but were not executed here; full native compilation remains
blocked by the real asupersync dependency's memory-limit termination.

## Native D3Q19/D3Q27 open velocity and pressure faces

`lattice3::open_boundaries::OpenBoundaryGrid3` binds existing face specifications
to the original dense native collision, pull stencil and three population
buffers. Physical donor moments include the uniform half force. The prescribed
velocity or density combines with the first interior cell's other moment and
physical nonequilibrium stress, then the native equilibrium, second-moment
projection and complete half source reconstruct the accepted populations.
Quadratic/product equilibrium and collision rates remain explicit choices.

Open faces occupy one axis, with a fluid donor one cell inward. Fixed exterior
and voxel wall links retain their original populations at open/wall seams;
only open-owned rows change there, so exact full prescribed moments apply to
face interiors. Nonzero wall motion and incompatible topology refuse. Both
grid and selected cell-collision Mach envelopes apply before publication.

Serial and pooled steps share the original 64-cell kernels. No provisional
open reflection is admitted as a final reservoir state. Populations, step count
and the measured mass/linear/angular reservoir exchange publish together after
the final cancellation checkpoint. The receipt includes temporary open-link
reflection and every actual regularization change, with canonical compensated
reduction; it is not a pressure-work or thermodynamic-energy balance. Target
rebinds alter held values atomically without advancing or resetting fluid state.

`OpenBoundaryGrid3::{step_serial_with_open_boundaries,
step_pooled_with_open_boundaries}` accepts new uniform reservoir targets in the
same native transaction as fluid advancement. Existing topology, boundary-kind,
profile-ownership and collision-envelope admission still apply. The original
collision/stream/measurement path uses the candidate face plan; the plan publishes
at the population/step-count commit point after the final cancellation checkpoint.
A late failure therefore retains both the old controls and the old fluid state.
The held-target step and explicit rebind APIs share the same preparation and
stepping implementation without a second solver or population copy.

`OpenBoundaryGrid3::set_pressure_profile` supplies a held scalar density field
on either pressure face for transverse hydrostatic heads and localized pressure
loads, including full D3Q27 collisions. Two reusable face buffers preserve the
accepted trace on source, density or cancellation refusal. Serial/pooled steps,
wall ownership and actual linear/angular exchange retain the original kernels;
sampling never occurs in a worker. An installed profile cannot be silently
overridden by a uniform-target rebind. `tests/native_pressure_profiles.rs`
exercises two-direction hydrostatics in every reservoir axis, off-axis pressure
impulse/torque, wall seams, worker parity and transactional face-buffer reuse.

`tests/lattice3_open_boundaries.rs` covers uniform flow for native collision
families, physical forced velocities/stress, independent total momentum/mass
balance, wall/open ownership, moment-origin covariance, worker invariance and
late cancellation/admission rollback. This owner currently retains uniform
body forces and fixed geometry; spatial/endpoint force fields, moving walls,
spatial velocity profiles and nonreflecting characteristics remain separate.

`OpenBoundaryGrid3::{step_serial_with_material_rheology_workspace,
step_pooled_with_material_rheology_workspace}` resolves each original material
at the incoming physical strain rate and applies its molecular `mu/rho0` to
the actual native shear relaxation. Velocity reservoirs use quadratic on-site
one-sided derivatives, with a linear fallback if the second interior donor is
solid. Pressure reservoirs retain their zero-normal-velocity-gradient donor
condition. Interior/periodic and tangential halfway-wall derivatives preserve
the shared grid's stencil, including open/wall seams. Uniform forced plug flow
therefore does not acquire false wall shear at a reservoir. Bulk and higher
rates, TRT magic and the equilibrium choice remain explicit native settings.

The reusable `RheologyWorkspace3` reserves velocity scratch and two original
source epochs; successful fluid/reservoir publication swaps the accepted epoch
without adapter-owned step allocations. Source, scaling, native admission or
cancellation failure preserves accepted populations and the prior report.
Sources own temperature, composition, history and any external side effects.
The adapter is an explicit generalized-Newtonian update; it does not add
implicit stress iteration, SGS viscosity or thermal evolution. Focused module
regressions cover physical plug and quadratic shear in every reservoir axis,
wall seams, solid-donor fallback, actual collision coefficients, pooled parity,
buffer reuse and failed-trial retention.

## Physical velocity initialization on shared D3Q19/D3Q27 grids

`Grid3::try_new_equilibrium` accepts physical lattice density and velocity.
It stores equilibrium minus the full half source, retaining the requested
velocity and `rho (cs² I + u u)` physical stress. D3Q19 and quadratic D3Q27
use Guo's source; factorized D3Q27 BGK/TRT reuse their complete product
derivative. D3Q27 moment families start from the product equilibrium with
the same hydrodynamic initial moments; their finite-force accuracy scope
is unchanged. Both grid and collision Mach limits apply to requested and
represented physical velocities. Initial populations are never clipped.

`ForcedGrid3::try_new_preserving_velocity` attaches a spatial initial force
field to an unstepped zero-force grid while preserving density, physical
velocity and existing nonequilibrium stress. This changes source coordinates
without an impulse, elapsed time or heat. Force shape, finiteness, solid
ownership and collision capability are checked before publication. Candidates
reuse existing population scratch and poll cancellation at fixed tile bounds.
Zero-force populations retain their bits. The existing raw constructors keep
their explicit population-frame behavior.

Native regressions cover every collision family, full physical stress and
nonequilibrium preservation, all 27 product-source moments, first-step momentum,
zero-force identity, the physical Mach convention and invalid/cancelled inputs.

## Spatial reservoir pressures on native SIMPLE faces

`SimpleConfig::pressure_profiles_pa` optionally replaces each Pressure/Outlet
plane's uniform value with the static pressure at every original Cartesian
boundary face. The profiles retain the original staggered geometry, ghost
pressure momentum load, SIMPLEC correction, body-force coupling and transient
state. Plane values use x-fast indexing with the normal coordinate omitted;
solid faces receive no load. No interpolation or plane averaging is implicit.

The numerical datum comes from a fluid-exposed profile face; unused solid-face
values must remain finite but cannot set a gauge or enter pressure differences.
The datum is removed before momentum assembly and restored only for pressure
readback. Restart admission checks active profile values against the retained
datum. The existing Stokes reservoir-traction heating mode samples
pressure times volume flow at each original face, preserving spatial covariance
and the same relative-gauge work account. Uniform profiles reproduce the existing
uniform solve, and omitting profiles preserves its API behavior.

Focused G1/G3 regressions cover every axis, graded-grid hydrostatic loading with
throughflow, a discrete hydrostatic rest state, independently driven counterflow
channels separated by solid voxels, pressure-work/physical-heat closure and gauge
shifts. Malformed, nonfinite or overconstraining profiles refuse. These are fixed
spatial pressure data for the existing incompressible model; material density,
thermodynamic pressure, time-dependent reservoirs and exterior interpolation
remain with their respective owners. Exact rest-state verification is a spatial
operator check, not a new near-zero-flow convergence guarantee.

## Original-source conduction with phase/resistance feedback

The conductive phase-resistance thermal-volume momentum entries prepare
paired original-state isotropic conduction once through the existing
harmonic half-cell operator, then call the coupled phase/resistance owner.
Authored net heat and conduction enter each original M/h/P trial once;
physical drag heat, endpoint source viscosity/permeability, phase density
and compatible volume flux converge together before momentum transport.

Conductivity and permeability retain the same original caloric chart and
material card. Initial, incoming, accepted and drag-heated caloric
intermediate states must remain inside the conductivity source. A final
in-range mixed endpoint cannot conceal an out-of-source intermediate.
The retained conductivity identity, paired heat, explicit enthalpy-secant
mixing-weight bound and endpoint resistance source receipts remain visible.

Conduction assembly/balance, caloric representation/transport and coupled
heat/constitutive residuals share one original enthalpy budget. Their
aggregate plus native resistance work defects also shares the mechanical
energy budget. Controls derive conduction's volume clock/gates from the same
resistance/thermal controls, so a second divergent configuration is absent.
Ordinary and spatial-pressure variants retain one accepted pressure field
and one accepted set of material/momentum routes; no numerical loss is heat.

Focused tests combine sourced conduction with nonzero mushy drag, actual
density/volume change and conservative enthalpy transport, plus identity,
stability, source-bound and cancellation refusals. Conduction remains an
explicit original-state first-order split with its original monotonicity
limit; it is not reevaluated inside the nonlinear drag loop. Exterior and
fixed-solid thermal faces remain insulated. Exact-solid momentum, moving
interfaces and a pressure-dependent thermodynamic EOS remain outside scope.

## Endpoint-source phase, resistance and thermal-volume feedback

The resolved phase-resistance thermal-volume momentum entries close the
original enthalpy chart, actual-density continuity, physical Darcy heating
and material-source resistance in one bounded accepted step. Original and
trial endpoint phases resolve through the existing source owner. True liquid
volume fraction sets permeability; dynamic viscosity retains its current
source density, including kinematic-viscosity conversion.

Every returned endpoint source must match the exact phase, law and density.
The applied native beta field and endpoint source beta are retained as
distinct data. Their difference is bounded on accepted native velocities by
non-cancelling half-cell impulse and heat sums; these share original force,
caloric and local heat budgets. No arbitrary coefficient-relative tolerance
replaces the physical residual. One accepted pressure field supplies volume,
momentum and material routes. Numerical trial heat/impulses are never replayed
into the original inventories.

Original, incoming, returned and accepted caloric-intermediate states must
remain mobile. Exact solids refuse before their liquid-property query.
The intermediate check prevents hot incoming material from concealing an
unsupported solid pretransport split. Endpoint receipts are retained;
current M/V, the net-enthalpy pressure-work convention, bounded iterations,
explicit source shear coordinates and atomic cancellation remain unchanged.

Focused regressions cover an independent unequal-density coupled endpoint
(h=60 J/kg, rho=2.5 kg/m³, beta=5/27 Pa s/m², outlet flux=2/15 m³/s,
physical heat=1/1215 J),
source identity and phase refusals, concealed caloric immobility and retry.
The source model is still a bounded equilibrium/mobile fixed-grid law.
This adds no exact-solid momentum constraint, moving skeleton, kinetic
phase law, source shear inference or pressure-dependent thermodynamic EOS.

## Coupled physical Darcy heating and thermodynamic volume

The resistance-aware piecewise thermal-volume momentum entries solve a
bounded nonlinear heat/volume problem on the existing native owners. Every
trial starts from the original M, P and chart state. The once-authored net
external caloric increment plus a physical drag-heat guess advances h; its
chart density change and original enthalpy-mixing volume drive continuity.
The same accepted implicit resistance projection supplies actual-density
pressure/body/drag impulses and physical cell heat.

Acceptance requires agreement of caloric and native drag heat, including
addition roundoff, local heat limits and shared non-cancelling mechanical
and enthalpy budgets. Current M/P is reconstructed once. Only the accepted
impulses are scattered; the accepted mass routes then transport kicked P.
Iterations do not accumulate heat, impulses or elapsed physical time.
Numerical backward-Euler, reconstruction, scatter and transport losses
retain their separate accounts. Actual M/V remains the inertial density.

Spatial reservoir pressure preserves native relative gauges and explicit
incoming material/velocity. Required source-pressure and opening work retain
the original NET-enthalpy convention. Prescribed resistance is frozen in this
entry; the source-response seam also measures endpoint constitutive changes
as physical impulse and heat residuals, without a second accepted pressure
solve. Coupling/volume iteration and cancellation bounds remain explicit.
Exhaustion or any early/late refusal returns no candidate.

Focused tests cover an independent quadratic heat/expansion solution on all
axes, actual M/h/P and phase volume, zero-resistance reduction under large
pressure offsets, and failed-step/retry behavior. No exact-solid locking,
moving skeleton, free surface or pressure-dependent thermodynamic law is
introduced by this fixed-grid coupling.

## Phase/source-bound stationary resistance mechanics

The phase-resistance momentum entries resolve each admitted original fluid
state through its existing material source and `MushyPhaseLaw`, then set
`beta = mu K^-1` for the native implicit resistance owner. Permeability uses
the law's true liquid volume fraction, including unequal constituent
densities. Kinematic viscosity is converted by the original source resolver
at its chart-matching current density. Complete mobility/law identity, exact
h-resolved phase and source/transport density are checked; one original
source receipt is retained per fluid cell.

Original `M/rho_chart(h)` must fill its cell within an explicit volume
budget. That admission never replaces authoritative M or native M/V inertia.
Exactly solid constitutive states refuse before liquid-property queries;
fixed grid solids retain zero inventory and need no source. Fully liquid
states yield exactly zero resistance. Nonzero beta overflow or underflow
refuses. Ordinary and reservoir-pressure entries perform one fresh native
momentum step and retain physical cell drag heat for the caloric owner.

Focused tests cover unequal-density volume fractions, source nu-to-mu
conversion, analytic graded-column response, actual-mass retention,
fully liquid reduction, source provenance and immobility/refusal boundaries.
This is an explicit original-state mechanical coupling: h, phase and
thermodynamic volume do not advance, source shear is caller-declared, and no
coupled heat/volume fixed point or exact-solid momentum owner is claimed.

## Material transport with physical resistance heating

`advance_material_momentum_with_resistance` and its reservoir-pressure form
consume current M, h and P, perform one fresh resistance-aware mechanical
step, and reuse that accepted native field for all original material and
momentum routes. Accepted drag impulses enter original cell P once. Only the
native cell partition of physical `dt B u_final²` heat is inserted into h
before those same mass routes transport M, M h and kicked momentum. Incoming
source and reservoir material/velocity remain explicit, including backflow.

Actual `sum M (h_heated - h_original)` and absolute local insertion defects
are checked, so unrepresentable heat cannot silently disappear under a large
enthalpy datum. Heat insertion shares the existing mechanical energy and
material enthalpy budgets with their native stages. Failed reconstruction,
mechanics, caloric insertion, material or momentum transport, or cancellation
returns no candidate and does not mutate caller state.

Focused G0/G2/G3 tests cover graded density/resistance, source and boundary
flow, accepted heat and momentum routes, zero-resistance reduction, held
controller work and atomic refusal/retry. This is fixed-geometry transport
with frozen isotropic resistance. It does not update an EOS or phase law,
solve heat-dependent resistance or thermal-expansion feedback, or convert
numerical damping or pressure work into additional enthalpy.

## Density-aware implicit stationary resistance

The staggered `FlowField` resistance entries solve
`(M + dt B) u1 = M u0 + I_body + dt A delta_p`, with both M and B integrated
over the original graded half-cells. Caller-supplied nonnegative isotropic
`beta` has units Pa s/m²; B is its volume integral in kg/s. The pressure
operator uses `dt A²/(M + dt B)`. Actual density supplies inertia, kinetic
energy and material impulse shares; resistance never becomes fictitious mass.

The accepted endpoint provides physical `dt B u1²` heat, partitioned by
each contributing cell's own beta and half-volume. Separate backward-Euler,
reconstruction and scatter losses are numerical and never heat. Fixed
normal flux retains its controller reaction and endpoint work against drag.
Relative pressure readback preserves small heads under large common gauges;
volume sources and reservoirs retain the original pressure-work convention.
Cell feedback scatters accepted body, pressure and drag impulses once into
the original extensive momentum using authoritative half-cell mass.

Shape, finite coefficients, positive representable mobility, native volume
and residual budgets, non-cancelling impulse/work defects, heat-partition
roundoff and cancellation all gate publication. Zero beta retains the old
face arithmetic. G0/G2/G3 regressions cover analytic transient/Darcy limits,
graded density/resistance, source flow, held-flux work and gauge changes.
This is frozen linear resistance against a stationary matrix, with no
constitutive permeability inference, exact solid locking, moving skeleton,
phase evolution or nonlinear resistance convergence claim.

## Physical impulses on source-bound liquid trajectories

The original `MaterialFreeSurface2`, `MaterialFreeSurface3`,
`MaterialFreeSurfaceEnthalpy3` and `MaterialFreeSurfaceRheology3` expose
`apply_impulse_density` in N s/m³ per whole geometric cell. Conversion uses
the retained material reference density and velocity scale. The native
located liquid mass supplies the impulse inertia, including partial fill;
the input is never multiplied by a timestep or silently weighted by fill.

The 3-D receipt contains total momentum [kg m/s], mass [kg] and signed work
[J]. The 2-D receipt uses the original mass per depth, giving momentum [kg/s],
mass [kg/m] and signed work [J/m]. Explicit SI defect budgets and the original
scenario speed ceiling remain binding. Nonfinite or nonzero-underflow
conversions refuse. Actual SI momentum/work readback, including conversion
roundoff, is checked through a native receipt-admission hook before the
final cancellation poll and population publication. A failed conversion
cannot leave an accepted native state behind.

The thermal and source-rheology trajectories delegate to that same native
liquid transaction. Exact canonical density and retained force-corrected
stress preserve pressure/source coordinates to the native model's roundoff.
Enthalpy, phase and material inventories, accepted constitutive samples,
wetting history, held forces and physical clocks do not change. Mechanical
impulse work remains a separate kinetic account; no heat or constitutive
relaxation is fabricated at zero elapsed time. Subsequent native evolution
retains the existing transport, molecular-heating and material owners.

Physical conversion allocates one temporary input field. Sparse publication
stages only affected wet cells; dense publication reuses existing scratch.
Neither path clones a complete liquid, thermal or material owner.

## Instantaneous loading of retained liquid mass

The sparse D3Q19 owner scans the complete x-fastest input so loads on inactive,
gas, solid or piston-occupied cells cannot disappear from the account. It stages
only nonzero wet-cell proposals in canonical Morton/lane order and publishes
only after the full requested count and balances agree. An all-zero input
allocates no replacement storage; no full sparse owner clone is required.

`FreeSurface::apply_impulse_density` and `FreeSurface3::apply_impulse_density`
apply resolved zero-duration inputs to the original dense D2Q9 and sparse D3Q19
populations. Inputs are lattice impulse per whole geometric cell. Integrated
momentum/work is per out-of-plane depth in 2-D and a total in 3-D. Fluid
inertia is the canonical population sum; interface inertia is its original
tracked liquid mass. Equal impulses therefore produce larger velocity changes
at partial fill. A source specified per occupied liquid volume must supply its
fill factor; the API never silently multiplies or clamps the requested impulse.

The shared Guo-frame kernel preserves intrinsic density and corrected viscous
stress. Liquid density closes to the exact original ordered sum before momentum
and work readback, preserving fill, conversion thresholds and the native mass
ledger. The applied impulse, actual liquid momentum change, signed midpoint
work and stable kinetic change are independently reduced and budgeted.
Nonpositive, nonfinite or overfilled addressed interface mass refuses; no
unlocated conversion carry supplies missing local inertia.

All cells are admitted into reusable native scratch before one publication.
Unsupported gas/wall inputs, malformed fields, allocation, low-Mach, balance
or cancellation refusal retain accepted populations, forces, mass lanes,
topology, carry, capillarity and wetting history. No time, streaming, pressure
projection, conversion, heat or damping is fabricated. Subsequent native
evolution resolves the resulting pressure and interface motion.

## Dimensional free-surface laminar film

`film::NusseltFilm` evaluates the constant-property Newtonian no-slip,
shear-free planar base flow from actual SI flux, viscosity, density, tension
and effective acceleration. It exposes film thickness, mean/surface speed,
Reynolds number, wall shear and the wall-normal velocity/stress profiles.
Flux is conserved analytically; there is no Reynolds floor.

`long_wave_mode` uses the isothermal Benney surface-wave approximation with
caller-declared `kh` and `Re*kh` budgets. Its onset uses mean-speed
`Re=q/nu`, with critical `Re=(5/6) g_n/g_t`; it does not substitute the
two-wall Poiseuille threshold. Invalid inputs, dry flow, unrepresentable
derived quantities and modes outside the budgets refuse explicitly.
Base construction and mode evaluation allocate no buffers. Profile/mode
arithmetic has fixed order; thickness uses the existing deterministic math.

This is a locally planar, fully developed, laterally wide, near-onset
long-wave model, not a turbulence, contact-line, arbitrary-wavelength,
non-Newtonian or transient-pour certificate. The module cites its primary
derivation and includes independent flux, stress/power, threshold, capillary,
material-scaling and refusal regressions. Native Rust execution was
unavailable in the implementing session; numerical comparisons are separate.


## Instantaneous momentum input on the original fluid populations

`Grid::apply_impulse_density`, `BoundaryGrid3::apply_impulse_density`, and
the shared `Grid3`/`ForcedGrid3` D3Q19 and D3Q27 entries
accept one prescribed impulse density per original x-fastest cell. The input
is already time-integrated: density is retained and `rho * delta_u = J`.
The operation advances no collision, stream, geometry or time. Existing held
forces, relaxation parameters and boundary data remain attached to the flow.

One shared pointwise kernel evaluates the second-order equilibrium difference
and the matching Guo half-source change algebraically. It retains the full
force-corrected nonequilibrium second moment instead of equilibrating away
viscous stress. Actual population readback supplies momentum, mass roundoff,
signed midpoint impulse work and stable kinetic-energy change. Explicit cell,
low-Mach speed, momentum-defect and work-defect limits gate publication.
The 2-D path reuses caller scratch; the tiled 3-D path reuses collision scratch;
the shared 3-D grids reuse their next-population buffer. No new population or
force allocation is needed. All paths cancel at bounded traversal boundaries
and publish once, after the complete balance check. The shared grid retains
its original, possibly tighter Mach limit; `ForcedGrid3` uses each accepted
cell's source frame and retains the source history. The impulse is delivered
once, without being replayed by subsequent collision/forcing steps. These
entries retain hydrodynamic stress for every native collision family; they
do not assert translation or preservation of all higher kinetic moments.

`MaterialScaling::{impulse_density_lattice, impulse_density_si}` converts
N s/m³ using the same source density and dx/dt as the flow. Total 3-D impulses
include dx³; 2-D totals include dx² and are per unit depth. There is no second
timestep factor. Only fully occupied `Cell::Fluid` receives a dense impulse;
free-surface interface inertia belongs to its original mass owner. Static
solid entries must be zero. Partial-saturation fictitious-fluid subtraction
and reciprocal rigid-body reaction remain responsibilities of the enclosing
coupling owner. This is a low-Mach hydrodynamic velocity jump, not an
incompressible projection, compressible shock model or thermal source.
G0/G3 regressions in `tests/momentum_impulse.rs` cover independent moments and
stress, signed work, subsequent native collision, SI scaling and atomic refusal,
including all D3Q27 collision families with uniform and spatial force histories.

## Implicit anisotropic resistance on native D3Q19 momentum

`BoundaryGrid3::step_with_linear_resistance` and its pooled/selected-wall
counterparts admit the local law `F = b + D(v_s-u)` with a prescribed skeleton.
`LinearResistance3` owns finite nonnegative principal force-density/velocity
coefficients and orthonormal material axes; zero coefficients admit free
directions. Coefficients, skeleton velocity and independent body force are
held at their supplied midpoint values for one lattice tick. They replace
both old force components after preserving the initial physical velocity and
full corrected stress through the original Guo reframe.

Native collision uses the initial law. Native pull streaming solves
`(rho I + D/2)u = m_raw + (b + D v_s)/2` with three principal-axis divisions,
then publishes populations and endpoint forces together. The same grouped
kernels serve serial and pooled paths; no second collision or post-stream
population scan is introduced. Invalid data, cancellation and late refusal
publish neither state nor force. Velocity/pressure reservoirs use the joint
open-face/source solve described below. Nonlinear drag shares the extension
below; evolved skeleton inertia and heat deposition belong to coupled owners.

The receipt separates independent forcing, constitutive skeleton reaction,
their angular impulses/work, nonnegative midpoint relative dissipation and
selected native wall exchange. Complete momentum balances require every wall
link and close to floating-point precision; measured reservoir exchange and
periodic angular transport remain separate. Tiny constitutive drag is retained in the receipt even if it rounds
out of the total force, without a promise to accumulate sub-ULP fluid impulses.
The local drag/work identity is not a transported kinetic-energy equality.
Uniform resolved relaxation is second order and has bounded relative source
energy; stiff trapezoidal updates can reverse velocity and are not L-stable.

### Vector-speed Darcy–Forchheimer extension

`DarcyForchheimer3` and the serial/pooled `step_with_darcy_forchheimer`
methods add `-beta |u-v_s| (u-v_s)` to the same anisotropic linear resistance
owner. The finite nonnegative coefficient is force density per speed squared,
including any physical density factor. The quadratic term uses vector speed,
not separate principal-component magnitudes. Isotropic linear resistance has
a stable analytic quadratic endpoint; anisotropy solves one monotone speed
equation with safeguarded Newton and a 64-iteration cap. Individual positive
quadratic roots bracket its speed within a factor of sqrt(3), independently
of principal resistance contrast. Zero beta uses the exact linear arithmetic.

The original collision, pull stream, wall exchange and atomic force/population
publication remain shared. Skeleton reaction is the endpoint-averaged drag.
The nonlinear relative dissipation is
`beta/8 (r0+r1) (|v0+v1|²+(r0-r1)²)`, a nonnegative expression equal to its
trapezoidal force work, rather than the different midpoint-speed cubic.
Analytic quadratic endpoints, resolved nonlinear decay order, rotated
anisotropy, large principal contrast, wall/source/skeleton momentum and work,
serial/pooled parity, exact beta-zero reduction and nonlinear atomic refusal
are exercised in `tests/boundary_forcing.rs`. The preceding prescribed-skeleton,
floating-point, no-heat and stiff-time-resolution limits apply.

### Joint porous resistance and open reservoirs

Both resistance laws admit the native velocity/pressure reservoir geometry and
spatial traces. Face interiors evaluate the constitutive law at the original
regularized physical target. Mixed wall/open cells retain every wall-owned and
tangential population and solve `F = law((m(F)+F/2)/rho(F))` on the actual
reconstructed state. The native Guo transform supplies its affine row
derivative; the existing prepared 3x3 LU solves a bounded damped Newton step.
Force-scaled final admission checks the represented populations and force,
including weak drag; failure publishes neither populations nor forces. Mixed
seam macros retain the native partial-reconstruction convention and need not
equal full-face targets.

The shared receipt independently measures open population mass, linear and
angular transfer and actual periodic coordinate angular transport. With all
wall links selected, `delta_P = external - skeleton + wall_population + open`;
the angular balance additionally includes periodic transport, to floating-point
precision. Reservoir pressure-work and heat are not supplied by this account.
Zero resistance retains original open-step bits, and the existing closed
endpoint solve and arithmetic remain unchanged. Native grid clones preserve
complete kinetic state, topology and profiles; `clone_from` reuses candidate
buffers without reconstructing populations or rebuilding links.

Five additional physical/refusal cases in `tests/boundary_forcing.rs` exercise
all six face normals with spatial traces and corrected stress, mixed wall/open
edges on all three axes, mass/linear/angular closure, serial/pooled parity,
weak drag, exact zero-drag reduction, and late refusal with atomic retry.

## Evolving resistance with retained physical force history

`lattice3::sources::resistance::ResistanceGrid3::advance_serial` and
`advance_pooled` accept the complete resistance field at the next physical
endpoint. The existing accepted force and separately retained drag drive native
collision; the next Darcy tensor, nonlinear coefficient, skeleton velocity and
independent body force enter the original implicit endpoint solve after native
pull streaming. Laws, populations, force, drag and elapsed steps publish
together. Updating a law never silently changes the already accepted physical
momentum frame. Fixed-law `step_*` retains its original arithmetic and receipt.

`ResistanceEndpointExchange3` retains actual endpoint-averaged source impulses,
native wall impulses and measured periodic-cut angular transport. Each
`ResistanceWork3` separately reports impulse work, trapezoidal endpoint power
and their directly evaluated signed temporal difference. For drag d ON fluid
and relative velocity v=u-v_s, endpoint fluid-drag plus skeleton work equals
minus the nonnegative endpoint resistance-power quadrature. Their impulse-work
sum additionally contains R=-(d1-d0).(v1-v0)/4. R can have either sign when the
law changes and is never physical heating; independent body-force timing has
its own work remainder.

An explicit finite nonnegative caller limit bounds the sum of absolute local
relative-work remainders in lattice energy units before publication. Opposite
cell defects cannot cancel this admission check. This matters for rapid changes
in resistance: positive endpoint dissipation alone does not prevent an
underresolved trapezoidal update from amplifying kinetic energy. Callers must
resolve source histories in time and choose an application-appropriate budget;
this local source criterion is not a complete fluid-energy certificate.

Complete endpoint-law admission precedes collision, and all original collision,
Mach, endpoint representability, local force/work and cancellation refusals
remain atomic, including the law/drag history. Serial and pooled paths share
fixed64-cell arithmetic and ordered receipts on the original D3Q19/D3Q27
collision families. This does not evolve porous geometry, material state,
skeleton inertia or thermal energy.

## Implicit porous resistance on the shared D3Q19/D3Q27 grid

`lattice3::sources::resistance::ResistanceGrid3<L,Q>` adds frozen spatial
anisotropic Darcy/Brinkman and rotationally invariant vector-speed Forchheimer
resistance to the original shared lattice grid. Both the D3Q19 adapter and
D3Q27 quadratic/factorized BGK and TRT, central-moment and cumulant collisions
reuse the existing `DarcyForchheimer3` admission, endpoint solve and local
force/work arithmetic. Each collision keeps its own native source convention.

Construction requires an unstepped, zero-uniform-force grid and explicitly
declares its initial population frame. The initial force is solved from supplied
raw moments; attachment does not preserve a previously read unforced velocity,
reset populations or rebind an evolved state. Resistance, prescribed skeleton
velocity and independent body force remain immutable. Read physical macros from
the retained force-aware owner.

Accepted force and separately retained constitutive drag remain at the same
epoch as the populations. Native collision uses that force; the endpoint solve
is fused into the original pull pass. Population, force/drag and completed-step
buffers publish atomically. Fixed 64-cell arithmetic and ordered receipt
reductions retain identical serial and pooled results across worker counts.
The original wall probe measures every stationary exterior/voxel wall link.
Periodic angular jumps are measured from actual fluid pulls, including D3Q27
corner wraps and excluding reflected links. Neither boundary term is inferred
from a global momentum residual.

Initial and endpoint admission check that represented physical velocity recovers
the solved drag within 4096 machine epsilons relative, without a unit-sized or
body-force tolerance floor. Unresolved half-force cancellation refuses rather
than reporting a recomputed source different from the retained collision force.
Work uses actual retained endpoint drags and is independently checked against
uncancelled endpoint work. Constitutive loss remains distinct from total fluid
kinetic transport; no heat or finite skeleton inertia evolves here.

Seven new native numerical cases cover all collision families, rotated anisotropy,
nonlinear second-order decay, resolved stiff reversal, extreme resolved
resistance contrast, unresolved cancellation refusal, an independent push
transport oracle, complete wall/periodic exchange, exact zero-resistance native
parity, worker bits and atomic refusal/retry. All 38 existing boundary-forcing
cases also pass, preserving the specialized D3Q19 behavior.

The original integration target also typechecks against the full actual
production dependency graph.

The original low-Mach/collision bounds and spatial/time refinement requirements
remain. This frozen-law owner adds no moving/open boundaries, evolving skeleton
motion, force-frame conversion or complete fluid-energy guarantee.

## D3Q27 TRT for wall-bounded laminar flow

`CollisionModel27::Trt { tau, magic, equilibrium }` extends the original
`Grid27` collision/stream owner with even/odd population relaxation. Molecular
viscosity stays `(tau-1/2)/3`; `tau_odd=1/2+magic/(tau-1/2)`. Each force
parity uses its own trapezoidal prefactor. Quadratic equilibrium retains the
Guo source; factorized equilibrium retains its full product derivative.
Thirteen opposite pairs and the rest population require no moment transform
or heap allocation. Existing collision models and arithmetic stay unchanged.

Material collision fields, material rheology and SGS read/retune the same
model while retaining its selected equilibrium and magic product. Invalid
magic, nonfinite tau and an unrepresentable odd relaxation time refuse before
stepping. The original grid admission, cancellation and publication rules apply.

The `magic=3/16` wall-placement claim is restricted to uniformly forced,
fully developed planar Poiseuille flow on stationary halfway walls. It does
not establish curved-wall, variable-viscosity-interface, open-boundary or
turbulent-flow accuracy. Added native regressions cover both equilibria, all
three axes at tau 0.65/0.9/1.7, nonunit density, reversed flow, a BGK slip
control, force/stress moments, BGK reduction and material/SGS retuning.
Independent JavaScript lattice-equation checks exercised those planar cases;
Rust regressions remain unexecuted in this session because no shell/runtime
was available. They are tests to run, not retained native acceptance evidence.


## Moving shared film carriers through one fluid/body step

`MovingCarrierSurface3::Attached` binds body-coordinate `LubricationPlanePatch3`
geometry to an explicit carrier registry index. Several faces share one physical
mass, complete inertia, accepted pose and external center load. Equal body values
at different registry indices remain independent. Each registry entry must own
at least one patch; at most 64 carriers and 64 patches are admitted.

`step_with_lubricated_spheres_and_carriers` and
`step_with_dry_lubricated_spheres_and_carriers` reuse the original moving-sphere,
fluid and dry-interior correction owners. The existing shared-carrier resistance
solver sees all faces at their actual common midpoint pose and includes their
cross-face mobility. Each carrier receives one free half-drift, the sum of its
patch spin impulses plus one authored center couple, and one final half-drift.
The solved linear response already includes its one external impulse. Prescribed
surfaces retain their original translated geometry, axial spin and work port.

Per-patch reports retain separate reaction, couple, work and active-film counts.
One carrier report owns their aggregate transfer, actual trajectory momentum and
kinetic change, external work, free-rotation energy defect and work residual.
Carrier storage must be counted once regardless of the number of attached faces.
The independent moving-surface path uses the same rigid momentum/trajectory helper.

Every dry-correction probe starts from the same accepted carrier registry. The
chosen carrier candidate remains paired with the accepted sphere/fluid candidate.
Endpoint construction, report allocation and both conservative swept-gap checks
finish before fluid publication; subsequent carrier/patch publication is infallible.
Invalid associations, native budgets, inadmissible sweeps and cancellation preserve
all accepted owners. The empty registry retains the original singleton fast path.

This remains the original local planar missing-film model: geometry does not
certify a finite face extent or introduce a resolved wall, contact, sampled carrier
fluid interior, displaced volume, inferred gravity/buoyancy or thermal destination.
Full-grid regressions check an independent three-mass midpoint oracle, coupled
full-inertia wrench/work, unique moving storage and actual population momentum,
legacy single-face/empty equivalence, rejected trials and exact retry.

## Pressure-aware dense liquid initialization

`MaterialFreeSurfaceScaling::free_surface2_with_pressure` initializes the
original closed D2Q9 classifier with positive absolute liquid pressure at
occupied cell centres. The field includes any hydrostatic or capillary
contribution and is sampled only during construction; subsequent gas traction
remains a separate prescribed field using the same fixed datum. Numerical
density and the caller's excursion bound use the existing rho0*(dx/dt)^2
conversion. The original native force transform establishes physical rest with
each cell's actual rho*g before interface insertion and mass initialization.

Mass per depth includes that numerical density and the existing half-filled
boundary convention; this is neither a physical liquid EOS nor a mass-preserving
projection. No relaxation, damping, second solver or stored pressure map is
introduced. Row cancellation or refusal returns no incomplete flow. Regression
definitions cover loaded equilibrium, gauge and material scaling, four discrete
hydrostatic orientations, zero-load replay and late refusal/cancellation.
Native execution of these new regressions remains unverified.

## Conservative material pressure defaults in D2Q9

The established `MaterialFreeSurface2::step_with_gas_pressure` and
`step_with_gas_link_pressure` now admit the same native conservative proposal
and physical-unit receipt as their `_and_momentum` counterparts. They discard
only the returned receipt, so both linear and selected angular redistribution
remain active when existing callers add pressure. One private candidate serves
pressure, conversion and SI admission before publication. Original authored
pressure errors and cancellation retain `MaterialLbmError::FreeSurfacePressure2`;
local conversion refusals retain their native reason. Source and scaling errors
are unchanged. Native pressure convenience APIs use the same conservative
candidate and retain their existing error type, including the original
conversion refusal reason without a pressure-cell coordinate.

## Prescribed gas traction with conservative D2Q9 conversion

`step_with_gas_pressure_and_momentum` and
`step_with_gas_link_pressure_and_momentum` compose the original dense pressure
reconstruction with signed local mass-and-momentum conversion. Ordinary and
pressure-loaded conservative steps share one admission/streaming/conversion
sequence. Each accepted step uses one staged surface; a refused pressure load,
unsupported numerical state, unavailable conversion capacity or cancellation
publishes no liquid state. Pressure callbacks are sampled only on the actual
pre-conversion exposed cells or links; newly exposed links are sampled next step.

The same entries on `MaterialFreeSurface2` and `MaterialFreeSurfaceWetting2`
accept absolute pressure in Pa with an explicit fixed reference and geometry in
metres. SI reporting and sourced contact-line history admit the same native
candidate before publication. A source-domain refusal preserves both liquid
state and all wetting history. Uniform cell/link loads use identical native
traction arithmetic, including the existing Guo force-frame reconstruction.

The caller's combined gas/capillary density-excursion bound remains a numerical
admission limit. This composition does not add a gas equation of state, gas
storage, pressure-work heating or an angular-momentum conservation claim. The
returned physical receipt covers conversion only, with mass, linear momentum
and signed numerical kinetic change per depth in kg/m, kg/s and J/m.

Dense gas reconstruction evaluates the even Guo correction with the accepted
endpoint force `rho_next*g + fext`. The held external force stays in the
existing gas-equilibrium cache. At each actual missing gas link, gravity adds
the coefficient `A_q = w_q*(-u_n.g/CS2 + (e_q.u_n)*(e_q.g)/CS2^2)`;
one scalar closure gives `rho_next = sum(f_base)/(1 + sum(A_q))`, followed by
`f_next[q] = f_base[q] - rho_next*A_q`. The existing explicit boundary velocity
is retained. This adds no nonlinear iteration or heap buffer and preserves
callback order and cache reuse. Invalid positive-density closure refuses the
transaction. Tracked mass publishes after all gas closures pass, reusing the
original post-collision populations and cached fill without another buffer.
This is force evolution across a tick, without a held-load reframe.
Gravity-driven D2Q9 trajectories intentionally change; zero gravity retains the
historical population arithmetic, including nonzero external force. Regression
sources restore the actual endpoint pressure pairs, distinguish source-dependent
traction from conversion momentum, and pin exact zero-gravity gas populations.

## Source-backed conservative D2Q9 momentum conversion

`MaterialFreeSurface2::step_with_momentum` and
`MaterialFreeSurfaceWetting2::step_with_momentum` expose the original dense
solver's signed, capacity-bounded local mass-and-momentum conversion without
changing its material density, viscosity, surface tension or force bindings.
A single native proposal is shared with physical-unit reporting and the sourced
contact-line update; native, SI-range, source-domain and cancellation refusals
publish neither flow nor wetting history.

The returned `MaterialFreeSurfaceMomentumReport2` retains the native receipt
and converts mass, linear momentum and signed kinetic remap change per unit
out-of-plane depth using `rho_ref*dx^2`, `dx/dt` and their products. The units are
kg/m, kg/s and J/m. Numerical remap energy is not deposited as physical heat.
This receipt covers conversion only; it does not certify full-step force/traction
balance or a finite-depth three-dimensional flow. The linear-only entries do not
preserve angular momentum; the opt-in angular mode below adds that conversion
account. The original native numerical-domain, capacity and physical scaling
limits still apply.

## D2Q9 cell-centre angular conversion

`FreeSurface::step_with_angular_momentum` and the source-bound
`MaterialFreeSurface2` / `MaterialFreeSurfaceWetting2` methods retain their
explicit local angular step. The native and material
`set_conversion_redistribution` selectors additionally admit
`ConversionRedistribution2::ConservativeLocalAngular` through the common
momentum candidate, including prescribed cell/link gas pressure and subsequent
source-owned wetting composition. `ConservativeLocal` remains the default.
Native ordinary and pressure convenience steps also honor this selection.

The original signed, capacity-bounded adjacent transport still owns mass and
linear momentum. Sparse components contain only actual nonzero transfer routes;
disconnected fronts cannot exchange corrective impulse or torque. Original
signed sources are captured during the class cascade, without an additional
full-grid snapshot, and survive later reuse as recipients or debt payers.
Positive final masses determine a centroid and scalar planar inertia. A
translation plus rotation supplies the minimum mass-weighted velocity
correction for simultaneous linear/angular closure.

Existing population retargeting preserves density and physical viscous stress
in the held Guo force frame. Actual population readback must close mass, linear
momentum and angular momentum to their scaled roundoff allowances. Original
transport defects remain accounted for alongside final component residuals.
Unresolved spin, inadmissible Mach number, nonpositive corrected populations,
unit overflow or cancellation refuses the complete candidate. The original
single native candidate also owns pressure loading, SI accounts and wetting
source/history acceptance; pressure traction and its physical torque precede
the conversion inventory and are not removed by the numerical correction.

Angular momentum is conserved at cell centres in the selected inertial lattice
frame. The remap does not conserve first mass moments and therefore does not
promise general boost-equivariance. Periodic coordinates refuse because
transport across a periodic cut has no angular owner here. The optional angular
defect sums absolute residuals about each component's fixed local origin;
translating their bounds also requires the separately reported linear
residuals. SI angular quantities are per depth [kg m/s]. Signed numerical
kinetic changes include the actual corrected velocities and remain separate
from physical heat. These are conversion accounts, not full-step force/torque
certificates.

Regression definitions cover unequal off-axis fronts under held force, a
linear-only torque-producing control, separate selected frames, translated
coordinates/origins, detached signed-debt components and retained sources,
unresolved-spin rollback, periodic/cancellation refusal, unchanged
no-conversion dynamics, pressure traction with positive and negative
conversion mass, SI angular units and gauge shifts, and dynamic-wetting source
refusal/retry. Route processing polls cancellation without introducing a dense
angular factorization or another full-flow clone.

## Rigid-frame momentum forcing on the native D2Q9 core

`core2::frame::RigidFrame2::prepare_step` supplies prescribed translation,
centrifugal, Euler and Coriolis acceleration to the original D2Q9 collision,
streaming and free-surface paths. All coordinates and rates use lattice units
(`dx = dt = 1`), positive spin is counterclockwise, and physical gravity is
supplied already resolved in the moving frame. The owner explicitly replaces
`g` and `fext`; pressure, viscosity and boundary traction keep their native
owners.

The local Coriolis source uses the exact two-by-two implicit-midpoint solve.
Its isolated force step rotates momentum without increasing kinetic energy;
there is no velocity rescaling, nonlinear iteration or additional flow solve.
The complete spatial/free-surface scheme does not acquire an energy-exact
claim. Motion coefficients are held over one step; scheduled callers provide
midpoint coefficients and refresh after every native step.

Changing forces first transforms the existing Guo half-source populations at
fixed physical velocity and stress. Preparation therefore applies no impulse
and advances no time. Two linear passes admit every live cell before publishing
anything, without a dense scratch allocation. Nonfinite motion, inconsistent
public storage and unrepresentable cell states refuse atomically. The existing
low-Mach, geometry and timestep-resolution requirements remain; frame motion is
prescribed, with no inferred finite-inertia support/vessel dynamics.

`tests/rotating_frame.rs` checks native collision/streaming response, 2,000-step
pure-Coriolis kinetic preservation, centrifugal and Euler impulses, free-fall
equivalence, unchanged physical stress during force replacement, and atomic
refusal. The acceleration decomposition follows
[MIT rotating-frame mechanics](https://ocw.mit.edu/courses/res-tll-004-stem-concept-videos-fall-2013/resources/rotating-frames-of-reference/);
the force representation follows [Guo et al. (2002)](https://doi.org/10.1103/PhysRevE.65.046308).

## Native terrain velocity evolution

`terrain_wind::run_v04c_terrain_wind` now evolves the existing D3Q19 BGK
boundary grid above a synthetic Gaussian hill. The authored-heightfield entry
accepts physical elevations, the caller's native runner and cancellation gate.
Cubic physical spacing and explicit `timestep_s` determine the velocity unit;
`tau` determines the reported laminar viscosity. Grid, memory, finite input,
open-face-neighbor and low-Mach admission precede a successful receipt.

The x inlet imposes a roughness-shifted logarithmic profile, x outlet fixes
reference density, y is periodic, the lower wall and voxel terrain are still,
and the upper wall moves at the reference top speed. Initial fluid is at rest;
`num_steps` advances actual populations. The original partial reconstruction at
wall/open seams remains. Final vector discrepancy, mean velocity, maximum speed
and downstream reverse flow come from interior fluid observations, excluding
the two x boundary columns on each end. Reverse flow requires a negative
streamwise velocity, rather than merely slow positive flow. Both RMS and
maximum discrepancy must meet the declared band. Schema v2 binds exact input
bits, including timestep, viscosity control, terrain and step count.

This is a finite-time laminar numerical experiment. The logarithmic input is
not a turbulence/wall-function model; no atmospheric LES, open-sky condition,
steady convergence, matched fs-atmo statistics or experimental validation is
claimed. The default resolution is 48×48×12 on 200×200×50 m so all three
spacings are equal. The existing battery covers rest, live time/viscosity/terrain
response, physical-velocity rescaling, admission and cancellation.

## Prepared analytic rheology in the native sparse free surface

`d3q19::freesurface3::rheology::FreeSurfaceRheology3` supplies prepared
Newtonian, power-law, Carreau and regularized Herschel--Bulkley laws to the
existing sparse VOF wet-cell BGK/TRT relaxation hook. This complements the
pre-existing source-resolved material rheology owner. A uniform law requires
no dense field; a spatial `RheologyField<3>` covers the full domain in fixed
Eulerian global coordinates, including newly wet cells after sparse-slot
reshuffles. SI laws prepared through `MaterialScaling::bind_rheology` enter
this same constitutive path.

The shared force-corrected stress inversion supplies persistent Morton-key/
lane relaxation warm starts and endpoint floor/cap statistics. Each step
uses one native candidate clone and admits every endpoint wet cell before
publishing populations, phase, mass, relaxation and statistics together.
Law rebinds and typed constitutive, native and cancellation refusals are
atomic. The returned collision statistics and retained next-endpoint
statistics describe distinct physical states. Gas and unlocated carry have
no constitutive stress; zero-mass wet collision ghosts do.

Fallible preparation and inversion reuse the original local kernel, while
legacy APIs keep their panic contract and successful arithmetic. A positive-
infinite trial strain may recover a finite limiting root; an unrepresentable
accepted strain refuses. Five focused regressions cover manufactured stress
and actual forced collision, Newtonian/uniform-field reduction, phase/slot
association, rollback and fallible-kernel limiting behavior. No advected
material label, inferred mixture law, thermal owner or yield history is added.

The prepared owner accepts native TRT without changing its authored even/odd
relaxation product. The full force-corrected even stress enters the same local
constitutive inversion as BGK. Every resolved even relaxation also admits the
actual derived odd relaxation before the source epoch can publish; finite even
viscosity alone is insufficient. Unsupported collision families still refuse,
and the original constitutive floor/cap statistics and numerical limits remain.

`FreeSurfaceRheology3::step_with_planar_wall_velocity_field` connects prepared
Newtonian, power-law, Carreau and regularized Herschel--Bulkley viscosity to the
original tangential face-field wall transaction. Accepted local rates govern
its actual collision and reflected populations. The actual endpoint stress,
including wall-generated shear and the Guo half-force correction, then refreshes
the same law before populations, mass, phase, relaxation and statistics publish.
No stationary exterior or gas velocity is substituted into a gradient. The
returned native face/total impulses, moments, couples and actuator work belong
to that same accepted flow step; applied and next-endpoint constitutive
statistics remain distinct.

Two Morton-keyed relaxation maps retain existing entries across steps and law
rebinds, pruning gas/retired keys and admitting newly wet cells at their global
field coordinates. Accepted rates supply warm starts while the other map is
private scratch; a failed native or endpoint constitutive solve cannot replace
the accepted epoch. One native candidate is cloned, with no unused thermal
mass-route trace. Wall/domain budgets and the existing fixed, completely flooded,
tangent-only plane, Mach, positivity and mass restrictions remain. This is
isothermal explicit constitutive splitting, without inferred heat, normal
geometry motion or an implicit wall/material solve.

Focused regressions cover independently reconstructed forced TRT collision with
an odd nonhydrodynamic mode, all-face Newtonian/native population and wall-account
identity, wall-generated Carreau shear and work, retained map entries, derived
odd-rate refusal and complete endpoint rollback with exact retry. They make no
new measured accuracy, convergence or performance claim.

## Exploratory Smagorinsky transport on the native sparse free surface

`FreeSurface3::prepare_sgs` and `step_with_sgs` accept explicit shared
`SgsConfig` controls and a reusable `FreeSurfaceSgsWorkspace3`. This path
requires Smagorinsky with native BGK, TRT or independent `CentralMomentBulk`;
other central/cumulant collisions refuse. The full-gradient WALE/Vreman path
admits the same three families. Independent central-moment collision preserves
each actual uniform or source-local molecular bulk rate and higher-order rate
while modeled viscosity changes only deviatoric shear relaxation.
WALE and Vreman need information absent from a symmetric local moment and
refuse on this stress path; the separate full-gradient path is described below.
Acceptance use, invalid controls, and excess eddy/molecular viscosity
refuse without changing accepted flow. The zero coefficient follows the
original native collision path bit for bit when the endpoint bulk field remains complete.

`FreeSurface3::set_trt_magic` selects TRT with an explicit positive even/odd
relaxation product while preserving the exact molecular even relaxation and
all populations. The checked selection is atomic; `collision_model` reports
the authoritative molecular law. Native local material/SGS relaxation uses
the shared `D3q19::with_shear_relaxation` owner, preserving the selected family
and authored TRT product as local even viscosity changes. SGS preparation
also refuses unrepresentable derived odd relaxation, before cache acceptance.
Fixed product is not a viscosity-independent wall-error guarantee under
arbitrary forcing, variable viscosity or free-surface geometry; the native
Guo source retains its second-order force moment (Postma and Silva, Physical
Review E 102, 063307, doi:10.1103/PhysRevE.102.063307).

The shared Guo half-source correction uses each wet cell's actual retained
force and physical velocity. Removing the stress trace yields `tau*S_dev`;
the shared Smagorinsky closure then supplies `K = tau*nu_t`. The rationalized
positive root of `tau*(tau-tau_molecular) = 3*K` determines current relaxation
without a lagged viscosity or a stencil through gas. This follows Chiappini,
Water 12 (2020), 1212, equations (6)-(8), doi:10.3390/w12041212, restricted to
deviatoric shear. TRT's stress is even, so the same root supplies its even
relaxation; the native collision applies each Guo parity at its own rate.
Both BGK and TRT still tie shear/bulk relaxation. Independent central moment
uses its original bulk rate to reconstruct the spherical stress separately.

Morton-key/lane cache entries survive sparse slot changes; newly wet cells
are resolved on the next preparation. Only one native candidate is cloned.
Its original collision, gas reconstruction, capillarity, mass exchange and
conversion finish before publication, so native or cancellation refusal
preserves every accepted component. Scratch is consumable only on success.
`FreeSurfaceSgsReport3` describes the pre-step input: instantaneous resolved
and modeled shear dissipation are weighted by native bounded liquid fill.
`resolved_bulk_dissipation_lu` separately reports molecular `zeta*(div u)^2`.
BGK/TRT's extra tied eddy volumetric damping is excluded from molecular power
and modeled deviatoric power. A bound local bulk field must exactly cover all
endpoint wet cells; native SGS refuses a topology change that needs a source
response it cannot resolve. The accepted flow remains unchanged on refusal.
Empty interface ghosts and unlocated carry add no liquid dissipation. No
diagnostic power is added as heat to this isothermal owner.

`MaterialFreeSurface3::{step_with_sgs_si,step_with_gradient_sgs_si}` expose
these ordinary steps with an authored physical filter width, original frozen
molecular coefficients and one native candidate. Both endpoint scenario-speed
bounds, SI power conversion and cancellation finish before publication. The
three reported SI powers distinguish molecular shear, molecular bulk and
modeled shear; they add no heat deposition or pressure-work account.

Six focused G0/G1/G4/G5 regressions in `d3q19::freesurface3::sgs::tests`
cover independently reconstructed forced populations and actual collision,
trace exclusion, all three zero-control conversion trajectories for BGK/TRT,
active mass conservation/replay, viscosity/usage refusal, cancellation and new
wet cells after sparse activation. Variable-viscosity TRT populations include
an odd nonhydrodynamic mode so the oracle detects either BGK substitution or
a frozen odd rate that changes the authored product. Selection and derived
odd-rate refusal are atomic. This is an exploratory constitutive integration,
not a turbulent-flow calibration, subgrid interface or entrainment model,
wall-damping law, resolved heat equation, or free-surface energy certificate.

## Source-caloric molecular viscous self-heating

With `enthalpy-coupling`, `MaterialChtFluid` exposes
`march_natural_convection_with_enthalpy_and_viscous_heat` and
`march_natural_convection_with_enthalpy_and_material_solids_and_viscous_heat`.
Both use the existing backward-Euler unsteady-Stokes work owner on the
original projected MAC faces and source-start local molecular viscosity.
Its cell heat enters the existing thermal power field once, after authored
power scheduling, and stays fixed during every caloric/radiative trial.
The source law h(T), joint solid storage/conductivity, contacts and compact
storage keep their existing owners. An empty authored source field is reset
between intervals, so previous molecular deposits do not accumulate as power.

Generated heat at fixed-temperature cells is reported as equal physical
source and thermostat outflow. Duplicate fixed declarations do not double
count it. The original nonlinear energy-closure residual is retained across
these cancelling additions. Accepted mechanical receipts separately retain
wall/body supply, native face kinetic storage, backward-Euler numerical loss
and pressure-constraint work; only molecular work enters thermal storage.

`MaterialViscousHeatMarch<M>` returns the original thermal result, the native
`UnsteadyHeatingAudit` and one mechanical receipt per accepted interval.
Audit acceptance follows thermal and source endpoint acceptance. The native
closed-wall/symmetry Stokes admission, Reynolds ceiling at most 0.1, physical
work tolerance, and no fan/resistance/SGS restrictions remain. Heated intervals
restart BE history from exact accepted faces; unheated entries retain their
existing history path. Source feedback remains explicitly first order.
No second momentum solve, gravitational store, numerical-loss heat or inferred
turbulent heating is introduced. Four physical regressions cover evolving
source Cp/mu, thermostat rejection, kinetic/time-loss history and joint solids.

Lattice Boltzmann core with a D2Q9 BGK path, a tile-major D3Q19 BGK/Guo path,
link-wise D3Q19 wall/open-boundary topology, the lattice-scaling assistant,
and dense-grid extension scaffolding for vector forcing, local rheology,
thermal double-population fixtures, and a free-surface mass-ledger prototype.

## Source-backed surface contact in native conjugate heat transfer

`ThermalSetup::source_face_contacts` binds original ordered interface cards to
explicit cell-face pairs. Each law updates only its two declared absolute-
temperature coordinates; pressure, medium, finish, history and the selected
claim retain their original context. The resistance is queried at the solved
surface temperatures, not the bulk cell centres. The existing half-cell owner
integrates native authored k(T) curves and recovers both surfaces from one
common heat flux. Explicit resolved conductivity fields keep their declared
frozen-bulk semantics, including the fluid's eddy contribution exactly once.

Native and resolved-field energy entry points share one bounded contact
iteration. Each trial reuses the original conservative fixed-resistance solve;
transient trials retain the original stored-energy history. There is no second
time integrator, interfacial storage or additional source of heat. Explicit
compact case/board joints add to the package branch resistance; their source
queries use the actual package and outside surfaces after the junction drop.
Duplicate ownership, compact sides/internal aliases and through-flow remain
refused. A failed query, cancellation or exhausted iteration publishes no
candidate and changes no caller-owned setup or thermal history.

Acceptance checks both the relative difference between source and assembled
resistance and the local Kelvin defect `(T_surface_A-T_surface_B)-R_source*q`.
`EnergySolution::source_contact_states` retains the actual surfaces, ordered
heat, applied/source resistance, residuals and original source receipts.
Explicit iteration, relaxation and tolerance controls govern this fixed point;
arbitrary positive source laws do not imply uniqueness or global convergence.
The no-source path retains its existing behavior. Focused G1/G3/G4 tests cover
closed transient energy against an independent scalar flux equation, nonlinear
solid/fluid half-cell integration, ordered source reversal, surface-only source
support, the heat observer, and transactional refusals.

## Purpose and layer

Layer L3 (FLUX). Depends on `fs-evidence` (the `Color` for the Evidence-typed
scaling plan), deterministic `fs-math` primitives, `fs-matdb` for immutable
interface-system queries and usage receipts, `fs-qty` for runtime
dimension checks, `fs-exec` for cancel gates and pooled sweeps, and
`fs-sparse` (COO/CSR assembly and ILU(0)) for the conjugate energy and
projection solves. Pure, deterministic (fixed tile/cell/link order).

### Material-owned thermal buoyancy in two and three dimensions

`MaterialScaling::thermal_slab_from_material` and `thermal_slab3_from_material`
consume `ThermalFluidTransport` from the original material owner. One coherent
state supplies momentum/heat diffusivity, isobaric volumetric expansion and the
actual reference temperature used by the body force. A separate expansion
constant or a wall-midpoint temperature cannot silently replace that source.
The caller's `MaterialBuoyancyBoundary` bounds the linearized fractional density
departure at both walls. Unknown expansion refuses; explicit zero and negative
expansion retain their physical meaning.

The corresponding `rebind_thermal_slab_from_material` methods also check every
current cell before the existing conservative rebind. Both relaxation fields,
buoyancy and reference temperature update together, preserving momentum and
temperature populations exactly. Density/Cp or lattice-unit changes still
require conservative remapping and refuse. Failed admission changes no field.
These remain frozen-property Boussinesq solvers: the density-departure check
does not establish phase/EOS support or bound future external heating. Focused
`material_buoyancy` tests exercise real 2-D/3-D evolution, pressure/composition
substitution, signed forces, off-center source temperature and atomic refresh.

The native `ThermalLbm` and `ThermalLbm3` also admit a finite nonzero
`set_buoyancy_direction` vector, normalized once before stepping. The vector
points opposite physical gravity; the existing signed `gbeta` and temperature
coordinate still determine the force magnitude. Defaults remain +y in 2-D
and +z in 3-D. Inclined layers and heated vertical channels therefore use the
same thermal populations, material scaling, local/tensor transport and Guo
force integration as the horizontal slab. Changing direction preserves the
current physical state; the next step rebinds the half-source frame without
an additional impulse. Checked source snapshots retain the direction, and
external-force admission checks every component of the total buoyancy/load.
The Rayleigh diagnostic retains its full coefficient; multiplying by the
wall-normal direction component gives its normal projection. The original
Rayleigh–Bénard onset bracket is not an inclined-channel stability criterion.
`directional_buoyancy` tests compare vertical-channel counterflow with the
analytic cubic viscous profile, grid refinement, rotated/reversed force
directions and exact per-step impulse through checked source advances.

`ThermalLbm` and `ThermalLbm3` expose opt-in
`BuoyancyQuadrature::EndpointAverage`; `OpeningTemperature` retains the original
default. The shared scalar operator advances once at the input physical
velocity. Flow then receives `gbeta * ((T_before + T_after)/2 - t_ref)` in the
authored direction through the original Guo rebind and collision. Symmetric
prescribed-source half-kicks make this mean equal the full interval's endpoint
mean. Old temperatures reuse nonoverlapping spare flow scratch, without an
additional per-step allocation or a second heat/flow solve.

Checked source, local-transport and tensor candidates retain the selection and
refuse a nonfinite transported total force before publication. Low-level direct
steps remain nontransactional after scalar transport. This native quadrature setting alone does not center caller-supplied external
loads. The ordinary density-resolved material adapter retains its held-force
sampling; the separate explicit material endpoint method below resolves that
force from its own original source epochs. This does not establish the time order of the complete partitioned
solver or a total mechanical/thermal energy closure. `buoyancy_quadrature`
tests use actual wall/source thermal ramps, independent fluid plus reflected-wall
momentum, time refinement, rotating loads, retained thermal history and
late-refusal rollback.

The opt-in 2-D/3-D
`step_with_boussinesq_material_transport_endpoint_averaged` methods return
`BoussinesqEndpointStep`. Initial material mu/k continue to set local relaxation
with fixed rho0/cp0. After one actual scalar transport, the original bounded
source is queried at the admitted scalar prediction
`T_transport + heat_per_step/2`; initial and predicted density departures are
averaged, prescribed external force is added once, and the original Guo
frame/collision advances flow once. The trailing heat half-kick retains the new
physical velocity and history. `force_endpoint` retains the distinct force-query
epoch; `transport.endpoint` is independently re-queried from actual final
populations, preserving exact source identity despite possible ulp differences.
All three query epochs, source/force conversion, heat/wall balance and
cancellation gates precede publication. No replay, inferred expansion law,
density/storage rebind or global coupled time-order claim is introduced.

The private source-owner regressions execute the actual one-pass 2-D/3-D
kernels: an independent quadratic-density thermal-ramp impulse has the analytic
trapezoidal defect and fourfold error reduction under timestep halving. Typed
force-query and final-endpoint refusals, nonfinite replacement loads and
continued accepted thermal history are checked. These two new checks and
eighteen prior selected native checks passed in isolated native-source execution.
`material_buoyancy_endpoint` additionally specifies nonlinear original-gas
density, one-read source tickets, all three exact query epochs, initial
molecular-rate ownership and material-source/cancellation rollback. Those two
full material integration tests remain unexecuted; selected kernel execution
does not establish full-crate integration or experimental accuracy.

`BoundaryGrid3::rebind_shear_relaxation` extends uniform transport refresh to
the original tiled boundary/moving-body solver. It delegates to the existing
`ViscousLattice3` rate owner, preserving the collision family, higher-order
rates, populations, forces, topology and work buffers. D3Q19's original
second-order bulk/shear tie remains. Invalid rates refuse before mutation;
the operation does not change mass units or add thermal energy.
`SphereLubrication3::with_dynamic_viscosity` preserves original film geometry,
tangential modes and work limits for the same material refresh. A coupled
consumer must bind both rates from one physical snapshot, as the material
immersed-sphere adapter does; neither numeric setter supplies material
identity or a conservative density change.

### Shared solid/fluid material feedback in natural convection

`MaterialChtFluid::solve_natural_convection_with_material_solids` resolves both
the fluid response and every solid's conductivity on the same steady temperature
field. `SolidConductivityLaw` retains the complete original material bundle and
pins its sources while changing only its declared absolute temperature. Scalar
or explicitly grid-aligned orthotropic coefficients feed the original conservative
half-cell/contact heat assembly; its temperature feeds the same buoyancy/viscosity
momentum owner. Solid conductivity defects participate in the existing nonlinear
acceptance gate, alongside fluid transport, force and temperature defects.

Applied and endpoint solid receipts are retained. Unsupported source states,
cancellation and exhausted coupling budgets publish no partial result. The
fixed-solid entry point retains its explicit constant-coefficient semantics.
This remains laminar, steady, fixed-geometry Boussinesq conjugate heat transfer,
without solid heat storage, melting, deformation or off-grid conductivity tensors.
Focused G1/G3 tests check independent nonlinear two-layer heat balances, data-only
material and principal-axis changes; G0 tests check source/cancellation refusal
and a sharp conductivity change hidden by a small temperature defect. A coarse
heated-cavity G3 regression checks that a solid-only source substitution changes
the actual fluid buoyancy and circulation while conserving heat.

### Source-caloric steady conjugate convection

`MaterialChtFluid::solve_natural_convection_with_enthalpy` and its
`with_material_solids_and_enthalpy` counterpart reuse the existing steady
material/momentum/heat fixed point with the original fluid's full h(T).
`ConvectionCaloricLaw` adapters borrow `ResolvedHeatCapacityCurve` or
`GasEnthalpyLaw` or the retained `IsobaricGasLaw`; the caller binds that same fluid and fixed source context to
the transport resolver. This correspondence is a caller contract; the chart's
chord Cp can differ from raw source Cp within its admitted interpolation
approximation, so no arbitrary pointwise equality is imposed. Every incident
fluid cell is also queried through that transport resolver at each inlet or
potential-backflow temperature, admitting its source domain and both density
and Cp departures before a solve. Direct original-law enthalpy increments establish one
common temperature datum, including every inlet/backflow state, so arbitrary
formation-energy offsets do not become numerical heat. The same donor h and
positive chart Cp enter both sides of each native upwind face. Reference rho
still owns continuity/inertia, and the existing rho AND Cp departure bounds
remain enforced; no physical approximation gate is bypassed.

`EnergyConfig::solve_with_enthalpy` supplies the steady affine operator.
`EnergyConfig::evaluate_with_enthalpy` applies that same native assembly to an
unchanged candidate supplied with its source enthalpy/Cp and molecular
conductivity fields. It checks the original algebraic and local physical heat
residuals and returns boundary, sink and fixed-cell heat accounting without a
Krylov solve. The caller owns original-law evaluation and nonlinear convergence;
enthalpy anchors differing from their endpoint temperatures, unresolved source
contacts, differing prescribed cell temperatures and compact aliases that would
be rewritten refuse. Field, flow and cancellation admission
remain shared. The public `source_caloric_endpoint` regressions cover nonlinear
enthalpy heat, stale-law refusal and unchanged prescribed/compact temperatures.

`EnergyConfig::evaluate_step_with_conductivity` and
`evaluate_step_with_enthalpy` extend that unchanged-state audit to physical
backward-Euler storage. Ordinary cells supply their original-law energy secants
from the fixed old time level; compact junctions refresh exact secants through
their existing storage owner and never borrow capacity from aliases. The same
native local heat residual includes storage, and the original step ledger
retains exact compact energy, prescribed-temperature jumps and physical sinks.
Source enthalpy anchors and limited traces are evaluated at the candidate;
no Krylov solve or constitutive update occurs. The public
`source_caloric_step_endpoint` regressions specify independent nonlinear
storage/advection balances for upwind and limited transport, compact C(T),
reservoir jumps and unchanged-state/refusal behavior.

The material owner re-evaluates the true caloric chart at the returned state
and requires both the maximum individual face-heat defect and the maximum
cell net-heat defect to meet the explicit tolerance, alongside its original
temperature, conductivity, viscosity and force checks. Interior cancellation
in the global balance cannot manufacture nonlinear acceptance. The result
retains the borrowed caloric owner, applied/endpoint fields and exact boundary
enthalpies. Source-domain, local mass-conservation, cancellation or iteration
failure publishes no solution. Existing source-solid conduction/contact paths
remain shared; latent plateaux, variable-density mass transport and pressure
work remain outside this steady single-phase entry. `material_conjugate`'s
`source_enthalpy` regressions specify actual channel heat, source Cp(T) changes,
reversed flow, datum invariance, local closure and atomic refusal.

### Radiation in steady material convection

`MaterialBuoyancyConfig::radiation` optionally borrows the original
`RadiationConfig` for all four steady material-convection entries, including
source-defined solids and source h(T). One cached `Radiators` geometry supplies
the first and every subsequent thermal tangent. Original fluid rho/mu/k, solid
conductivity, caloric transport and radiation are re-evaluated at the same
temperature endpoint. The smaller material/radiation pass budget applies.

Acceptance retains the original material and caloric gates and additionally
checks every cell's physical radiative power against the assembled tangent,
using the same normalization and roundoff allowance as native transient heat.
Relaxation or a loose temperature tolerance cannot waive this heat gate.
`MaterialRadiationReport` retains the original exposed-face count, heat escaping
to the surroundings, local power defects and radiation-corrected global heat residual.
The native `EnergySolution` report remains the actually solved tangent system;
the separate net-outflow correction excludes fixed-cell radiation because its
matching thermostat correction cancels. Surface exchange remains reciprocal.

No second heat solver, extra material query context or radiation ray trace per
pass is introduced. The absent-radiation path retains its original solves.
Existing laminar/Boussinesq/source-domain bounds remain; compact junction aliases
refuse this option because they do not own resolved radiating cell rows.
Focused G0/G1 regressions specify multiple faces per heat row, nonlinear-power
refusal, thermostat accounting, an independent Stefan–Boltzmann equilibrium,
and a source-solid/source-enthalpy channel carrying 6 W alongside radiation.

### Source-caloric native transient convection

`MaterialChtFluid::march_natural_convection_with_enthalpy` extends the SAME
native unsteady material march with source-law storage and upwind transport.
Reference `rho0 V` remains each fluid cell's mass; direct original-law
`h(T_new)-h(T_old)` replaces frozen-Cp sensible storage. Positive source
secants enter the original backward-Euler storage diagonal, while complete
donor h and positive tangents enter the existing conservative face operator.
Every nonlinear trial retains the original accepted old temperature and the
same projected native flow: a caloric iteration never advances physical time.
The native MAC momentum/pressure state and BE/BDF2 history remain intact.

Acceptance checks unrelaxed temperature change, EACH cell's true storage-rate
defect, EACH face and cell's nonlinear advective heat defect, and the actual
source-law heat closure. Global cancellation cannot hide local errors.
Fixed-cell reservoirs include their exact source-law storage jump, contacts
and fixed solid heat capacity retain their ordinary owners, and endpoint
radiation uses the original converged radiation law. The shared native heat
record reports true stored enthalpy and physical boundary/reservoir heat.
Explicit per-step caloric-pass and whole-march linear-energy-solve ceilings
bound work; the latter counts radiation trials as well. Domain/boundary,
conservation, budget, solve or cancellation refusal publishes no partial march.

Material density, viscosity, conductivity and buoyancy still refresh at each
accepted time level, with their original density AND raw/chart-Cp departure
limits. This is FIRST ORDER coupled time accuracy, even with BDF2 momentum;
the source coefficients retain the existing explicit lag. Caloric/transport
source identity and fixed pressure/composition correspondence remain the
caller's contract. No latent plateau, variable mass, pressure work, viscous
heating or new physical energy store is introduced. Focused `source_enthalpy`
transient cases specify independent heating/channel enthalpy recurrences,
reversed flow and datum shifts, native buoyant Darcy time history with exact
thermostat storage, and transactional refusals. These are regression definitions,
not an assertion that a native validation run or measured campaign completed.

### Compact junctions in source-caloric transient conduction

The source-caloric fluid, source-integral solid, and joint fluid/solid native
marches accept compact components with an authored positive
`heat_capacity_j_k`. Each nonlinear trial keeps that capacity in one appended
junction row of the original backward-Euler solve. Collapsed voxel aliases
have exactly zero cell capacity and are excluded from source-material queries,
source-integral storage, and source-solid energy subtotals. Ordinary cells still
require their original material assignment. The resolved-field operator admits
omitted material slots only for validated compact aliases; it does not create
placeholder material laws or relabel the domain.

The total native heat record includes actual fluid and resolved-solid enthalpy
increments plus exactly one `Cj * delta(Tj)` per junction. Prescribed-cell
reservoirs still supply their own actual source-law storage jumps. The shared
source schedule also scales junction power from the authored value each step.
The source-solid-only record and `solid_energy_j` continue to describe resolved
source solids, excluding compact storage. Junction temperatures participate in
the common endpoint temperature-convergence check; positive source secants,
local enthalpy defects, source support, cancellation and work budgets retain
their existing admission.

This remains a single isothermal junction with constant authored capacity, not
a package-internal temperature field or a nonlinear junction enthalpy model.
The separate radiation geometry still requires declared surface-material
slots; omitted bulk slots do not supply emissivity or imply radiating compact
surfaces. No material source, latent heat, or variable mass is inferred.

`tests/material_unsteady_buoyancy/compact_caloric.rs` defines independent G1/G3
quadratic heating/cooling and reservoir balances, unequal-area alias ownership,
source-solid subtotals and G0/G4 source-domain/cancellation refusal checks.
These are regression definitions; source inspection and parsing do not establish
an executed native validation result.

The cooling-CHT CLI admits these components in temperature-dependent
unsteady coupled-energy scenes. The source-completeness check covers authored
resolved solids, excluding the trailing compact material slot. Each component
still requires positive `heat_capacity_j_k`; the existing source domains,
fluid enthalpy controls and solid-storage controls remain enforced.
`cooling_cht`'s compact caloric regression checks the integrated CLI heat
ledger and verifies that its source receipts contain only resolved materials.

### Source calorimetry in native unsteady material convection

`SolidConductivityLaw::with_heat_storage` prepares the original selected Cp(T)
claim through `PreparedHeatStorage` on an explicitly declared temperature span.
`march_natural_convection_with_material_solids` then uses the SAME native
momentum history, projected fluxes, backward-Euler heat operator and time steps.
Every solid assignment supplies this original source law. Fixed reference
rho times voxel volume defines invariant mass; a changing source density cannot
create mass. Original source knots, conditions, selected claims and storage
identity remain retained, with no resampled enthalpy or invented Cp average.

Within each native heat interval, exact source-integral secants and grid-aligned
k(T) converge to a common supported endpoint. The acceptance gate checks the
unrelaxed temperature defect, conductivity defect, each cell's exact storage
defect, and the corrected whole-step heat closure, with explicit roundoff
allowance. The existing native contact, radiation, source and boundary owners
remain authoritative. Prescribed-temperature reservoirs supply actual integrated
storage once even with repeated declarations. Result heat records contain true
integrated storage; solid records retain nonlinear defects and solid increments.
Per-cell solid sensible energy uses the declared lower-span temperature datum.

Unsupported linear predictors may backtrack only their constitutive trial;
accepted temperatures are never clamped or extrapolated. Initial and endpoint
source failures, nonconvergence and cancellation return no candidate march.
This is single-phase fixed-geometry solid calorimetry at the source's declared
mechanical conditions, without latent plateaus, density-driven deformation or
inferred Cp/Cv conversion. Fluid properties and buoyancy retain their existing
explicit first-order coupling, even with BDF2 momentum. The old frozen-solid
entry retains its original behavior. Focused tests independently reduce a
piecewise-Cp two-cell energy balance, material substitutions, fixed-reservoir
storage/radiation, bounded large-step iteration, and source/iteration refusals.

### Exact source Cp and local solid storage admission

The original `IsobaricHeatCapacityLaw` implements `ConvectionCaloricLaw`
directly through its exact source-Cp segment integral and tangent. Both existing
steady and transient enthalpy entries can therefore retain every original Cp
breakpoint without constructing a resampled enthalpy chart.

The joint fluid/solid transient owner admits EACH solid cell's constitutive
storage defect against that cell's own floating-point allowance. Another
voxel's large reference mass or temperature cannot waive a local energy error.
Summed roundoff is used only by the whole-domain heat ledger; nonfinite
allowances refuse. Focused source-based regressions cover an independently
reduced coupled energy balance, repeated thermostats with radiation, source and
work-budget refusals, and a large unheated specimen that cannot hide a small
voxel's 30 J nonlinear caloric defect.

### Solved material flow on the shared enthalpy timeline

With `enthalpy-coupling`, `MaterialChtFluid::solve_flow` retains the actual
SIMPLEC face fluxes, immutable voxel domain and original material snapshot.
`MaterialChtFlow::advect_enthalpy` converts those same volume fluxes to mass
fluxes with the density that drove momentum. It reuses `fs-thermal::Advected`
for upwind enthalpy, boundary accounting and the combined transport CFL.
No cell-velocity reconstruction, replacement Cp law or second heat integrator
is introduced. The returned observer composes with existing current, material
interface and source histories on their common accepted thermal substeps.

Every fluid voxel maps one-to-one to an active thermal cell with the same
volume/centre, `mass = rho * volume` and actual borrowed caloric law. Inflow,
including outlet backflow, requires explicitly supplied supported enthalpy;
outflow carries the current local h. Internal advective watts never enter the
external-heat audit. The shared roundoff mass-balance admission is checked
before an observer is returned; a larger momentum divergence refuses without
repairing the flux or creating an artificial mass source. Geometry/storage
bindings are checked again at source epochs, steps and events. Refusal or
cancellation keeps the accepted thermal and source histories unchanged.

This is stationary, constant-density, frozen-momentum enthalpy transport.
A latent plateau remains caloric storage, not permission to infer the flow
mobility of a partly solid material. No phase-volume change, skeleton support,
free-surface motion, pressure work, viscous heat or eddy heat conductivity is
inferred. Those require their existing or explicitly admitted separate owners.
`tests/flow_enthalpy.rs` exercises actual solved flow, viscosity-data changes,
enthalpy/source balances and rollback; the module tests cover graded geometry,
all three flow directions, closed circulation and admission failures.

### Material-bound enthalpy transport

With `enthalpy-coupling`, `conjugate::enthalpy::transport::material::MaterialCell`
binds an active thermal cell to its actual borrowed caloric law and original
material card. One atomic source query supplies thermal conductivity and the
explicitly chosen electrical resistivity/conductivity basis at the same evolving
temperature and liquid mass fraction. The common material conductor owns any
reciprocal; original source claims remain on `MaterialTransportState`. Initial
claims are pinned for the binding's lifetime; evolving state cannot silently
select another measurement after leaving the original source's support. Explicit
insulation differs from missing data. The binding plugs directly into the
existing explicit and implicit current/enthalpy callbacks; link geometry and
thermal interfaces retain their existing owners. Endpoint source refusals prevent
publication. Fixed mass/geometry, scalar DC conduction and the declared caloric
chart remain the model limits; this does not implement moving melt or induction.
`tests/joule_enthalpy_material.rs` exercises source-driven phase/current/heat
changes on a latent plateau, equivalent electrical source bases and rollback.

### Moving-liquid enthalpy on native free-surface mass routes

`material::free_surface::enthalpy::MaterialFreeSurfaceEnthalpy3` binds the
original `FreeSurface3` state, material card, pinned `LiquidTemperatureLaw`
and `EquilibriumEnthalpyPhaseCurve`. Native collision emits its actual pair
mass exchanges; the original interface conversion cascade emits its exact
carry/redistribution transfers. The companion transports extensive `mass*h`
over those same routes with implicit donor upwinding and the existing sparse
BiCGStab/ILU owner. The solve uses specific-enthalpy differences from one
deterministic old-storage reference per connected route component, so an
arbitrary caloric datum or disconnected reservoir cannot suppress real mixing
through the relative residual tolerance. Born cells and empty conversion carry
inherit their connected storage reference. Any admitted roundoff mass-closure
defect retains its reference-energy contribution in the same donor equation;
neither masses nor fluxes are repaired. Datum invariance is therefore to the
native continuity and stored-enthalpy rounding precision. After restoring the
datum, the reported residual is recomputed from the enthalpies actually stored
against the same centered system. Endpoint acceptance requires the actual L2
residual [J/kg] to satisfy `absolute_transport_tolerance_j_kg +
relative_tolerance * norm(centered_rhs)`. The absolute allowance is an explicit
finite nonnegative control, defaulting to `1e-8 J/kg`; setting it to zero opts
into strict relative acceptance. Both actual absolute and relative residuals
are reported. This permits explicitly budgeted representation roundoff without
allowing a large datum to hide a materially inaccurate mixing result. The
Krylov iteration itself retains its original relative convergence criterion.
Extensive enthalpy diagnostics retain the original caloric datum. It does not
reconstruct a separate flux from velocities.
Positive old storage must anchor every participating route component; both
local mass closure and the final extensive-enthalpy defect are checked. The
conversion carry retains its own enthalpy. No regularization, balancing heat
or density repair is added.

Declared cell heat acts at the beginning of the step. The original caloric
chart supplies T, and every wet cell's pinned material query supplies viscosity
to its original native collision. All endpoint states are queried again before
publication. A private flow candidate makes source-support, negative-storage,
solver, budget and cancellation refusals atomic for populations, classifications,
carry, heat and material history. Zero-mass interface ghosts store no energy;
new unconnected ghosts receive a mass-weighted one-ring extension only for their
collision coefficients. Their enthalpy is never a heat reservoir.

`step_with_heat_and_gas_pressure` combines declared cell heat and an imposed
absolute gas-pressure field on this same timeline. An explicit pressure
reference fixes the Pascal datum; physical cell-center samples lower through
the same checked pressure-unit conversion as the isothermal material owner.
Temperature-dependent viscosity and the admitted local force field enter the
original native pressure reconstruction, while enthalpy follows that step's
actual mass routes. The pressure callback does not introduce a separate flow
solve, caloric source or heat term. Invalid pressure, density-excursion bounds,
source paths and observed cancellation leave flow, enthalpy, material states
and clock unpublished together. Reference pressure retains the ordinary heat
step's arithmetic; a common pressure-datum shift retains its physical load.
Focused regressions cover heated pressure-driven motion, the temperature's
effect on viscosity, mass and extensive-enthalpy budgets, reference/gauge
identity and late callback/cancellation refusal. This remains incompressible
thermal enthalpy transport: compressive heating, gas-reservoir evolution and a
total mechanical-plus-thermal energy closure are not supplied by this entry.

With `enthalpy-coupling`, `step_with_thermal_boundaries` accepts one borrowed
`FreeSurfaceThermalBoundaries` input combining contact resistance, prescribed
gas-film heat, direct radiation, solved interface temperature and optional
`FreeSurfaceGasPressure`. The original thermal coordinator still advances
finite solid bodies, attached solid radiation and bound electrical/Joule/Lorentz
physics; optional molecular heating retains its original account. The same
native flow candidate receives the imposed pressure alongside temperature-
dependent surface stress, with no extra thermal or momentum solve. Existing
entry points pass no prescribed pressure and retain their arithmetic.
Pressure is sampled in absolute Pascals at physical liquid-cell centers, with
the same fixed reference and combined density-excursion bound. It does not
rebind gas-film temperature/pressure source contexts, infer an evolving gas
volume or enter the thermal heat account. The caller owns compatibility of
that mechanical load with the declared source conditions. Finite-atmosphere
storage remains with its existing separate coupled entry points. A late
pressure or cancellation refusal discards the already-computed thermal stage
and finite-body candidate along with native flow. Focused regressions cover
conducting pressure-driven liquid with finite solid heat, reference-pressure
identity, late rollback and pressure acting beside solved surface-temperature
capillarity without weakening the complete-boundary continuation guard.

This bounded rung requires one pure liquid, unchanged physical density,
stationary walls and explicitly neglected capillarity. Negative accepted native
interface/carry mass refuses. Phase conversion, volume change, conduction,
viscous heating, pressure work and surface-energy evolution remain separate
missing closures. `h` is enthalpy and is never described as internal or total
energy. Implicit donor transport and local conversion redistribution introduce
numerical mixing within neighboring liquid cells; no physical thermal-diffusion
or mechanical-energy claim is made. The focused `material_free_surface_enthalpy` tests exercise actual
interface conversion, material-only viscosity changes, heat conservation and
atomic physical/source-domain refusals.

### Finite solid heat exchange with moving liquid

With `enthalpy-coupling`, `step_with_body_conduction` adds source-bound
`FreeSurfaceThermalBody` cells to the existing moving-liquid thermal timeline.
Each body keeps finite mass and specific enthalpy from its original immutable
caloric chart. Its scalar conductivity is resolved at the same current
temperature, with the initial source claim pinned throughout its lifetime.
Changing material data changes the actual contact heat, subsequent liquid
viscosity/buoyancy and native VOF populations.

`with_radiation_faces` also binds canonical total hemispherical emissivity from
that original solid card. Its selected claim joins the pinned conductivity
bundle and refreshes at each actual caloric state; known zero is reflecting,
while missing or unsupported data refuses. Explicit dry wall patches use the
same diffuse-gray fourth-power kernel as liquid radiation. Their areas fit the
remaining solid surface after contact coverage; gas films may share a patch.
All paths use the same source timeline and summed capacity-based timestep
bound. `body_radiation` and `heat_from_radiation_j` report the external radiant
supply separately from internal contact and finite-gas heat. The gas never
receives an opposite debit for the prescribed radiant bath. Source-only wall
swaps therefore reach liquid enthalpy, viscosity and native populations.
G1/G3/G4 regressions in `conduction/atmosphere/solid/radiation.rs` check an
independent four-store balance, finite-capacity subcycling and atomic refusal.
The fixed-wall, uniform-cell, opaque-gray approximation does not add melting,
resolved solid boundary layers, spectral/participating-gas radiation or visibility.

Finite wall cells may also retain explicit `FreeSurfaceBodyLink` geometry
through `with_conduction_links`. Each unordered solid pair is declared once,
with its actual shared area and two thermal-centroid distances. All existing
finite-wall, environmental and solved-surface entry points use those links.
The common observer reuses each body's pinned k(T) evaluation to form the
series bulk conductance; the original thermal owner applies equal and opposite
pair heat and includes every solid link in its capacity-based time-step bound.
Solid/solid and liquid/solid patches share each body's surface-area budget.
Missing/duplicate/self links, invalid geometry and work-budget excess refuse.

This permits nonuniform temperature and heat spreading through a discretized
wall without a second thermal integrator. Each cell retains its original mass,
caloric chart and solid-phase gate. Solid-internal exchange cancels from the
reported wall-to-liquid heat; applied and endpoint samples retain both original
material sources. All endpoint and native-flow refusals still publish none of
the liquid, wall or gas candidates. Body geometry is explicitly supplied and
fixed; perfect solid/solid thermal contact is the declared approximation.
Spatial/contact resolution and finite-volume accuracy remain caller obligations.
G1/G4 tests in `tests/material_free_surface_enthalpy/conduction/body/network.rs`
check independent three-store heat spreading, conductance-limited substeps,
combined energy, surface coverage, source-domain refusal and atomic rollback.

Explicit `FreeSurfaceBodyContact` patches must join a positive liquid cell to
a face-adjacent permanent solid tile or outer domain wall. Inactive gas tiles
cannot become thermal solids. Contact area, liquid/solid centroid distances
and stationary solid geometry are caller-supplied physical data; each native
wall face occurs once. The original thermal owner checks total surface
coverage and applies the two material half-cell resistances in series to both
the heat exchange and CFL bound. Undeclared surfaces remain insulated. Contact
resistance is explicitly neglected by this entry point; the source-bound
interface entry below supplies it without replacing either bulk layer.

Solid and liquid storage share the SAME substeps. Their exchange is internal
in the combined thermal audit; `flow.body_heat_j` records the signed transfer
into the liquid, separately from external and molecular heat. The corresponding
body report exposes both finite stores and the combined balance defect. The
viscous-heating entry point composes the original molecular source with this
same timeline. After thermal advancement, the original native mass routes
transport liquid enthalpy. All material, thermal, native flow, conservation,
budget and cancellation checks precede publishing either participant.

This remains a first-order thermal/flow split with fixed solid geometry and
constant caloric density, pure liquid, stationary native walls and neglected
capillarity. Solid cells refuse melting because this adapter supplies no
mechanical phase handoff. It is a finite-volume/uniform-cell model, not a new
mechanical solver, inferred wall geometry, interface-resistance model or
experimental validation. Focused tests in
`tests/material_free_surface_enthalpy/conduction/body.rs` compare independent
two-material heat/storage calculations, data-only substitutions and exact
insulation reduction; they also exercise source/cancellation and late native
flow refusals of the complete solid/liquid transaction.

### Ordered interface materials between finite walls and moving liquid

`step_with_body_interfaces` and its viscous-heating counterpart add the original
ordered `InterfaceSystemCard` through `FreeSurfaceBodyInterface`. The common
`ThermalInterfacePoint::resolve_resistance` selects one positive area-specific
resistance [m2 K/W] at the actual initial liquid/solid caloric states. Each
thermal substep and endpoint refreshes that same selected source with each
side's own temperature (and declared phase coordinates), retaining the complete
surface/medium/history context. Source support cannot silently switch claims.
Pressure and other external conditions stay explicitly prescribed.

The shared `fs-thermal::enthalpy::series_contact_conductance_w_k` composes
`A / (d_liquid/k_liquid + R'' + d_solid/k_solid)`. This SAME conductance supplies
both thermal exchange and the existing explicit stability limit. The original
finite stores, native mass-route enthalpy transport, viscosity and optional
material buoyancy therefore respond to interface-data substitutions. Pair heat
remains internal; molecular heat and external supply retain separate accounts.
The interface entry points also accept the existing environmental gas patches,
so wall contact, gas exchange and molecular heating share one source timeline
and combined thermal step bound. Gas heat remains external; wall exchange is
internal. Applied/endpoint samples keep the original interface receipt alongside
both bulk material states. A zero-area/insulated patch still queries its interface;
omitting an interface binding explicitly retains the historical zero-resistance
approximation. Duplicate or absent contact indices refuse before advancement.

All thermal, source, endpoint, native-flow and cancellation gates precede joint
publication. Fixed geometry, negligible contact storage, pure liquid and
neglected capillarity remain the physical scope. This does not infer contact
pressure, wetting area, gas heat exchange or a melting-solid handoff. Focused
G1/G3/G4 tests in the existing finite-body suite compare independent three-layer
heat calculations, actual interface-data changes in VOF populations, pinned
source-domain exits and rollback; synthetic fixtures are not experimental data.

### Environmental gas heat reaches the moving liquid

With `enthalpy-coupling`, `step_with_gas_environment` composes explicitly
exposed `FreeSurfaceGasContact` patches with the same liquid conduction,
finite-wall contacts, optional molecular heating and native VOF advance.
Each patch names a positive-mass liquid cell, an adjacent gas cell, physical
area, liquid-centroid distance and gas-film thickness. Duplicate faces,
permanent walls, wet gas neighbors and excessive coverage refuse. Geometry
must be supplied from the currently accepted interface at every call.

The existing series-conductance owner combines the original liquid k(T) with
the supplied `ThermochemicalGasModel` conductivity at the mean liquid/bath
temperature. Its species fractions, NASA9 laws, pressure and transport cards
remain source-owned data. Every thermal substep and endpoint reevaluates the
same sources, with the initial film transport claims retained. The thermal
step bound includes the SUM of gas and internal contact conductances at each
affected cell. Gas heat changes the actual liquid enthalpy, viscosity and
native populations; changing the gas data requires no solver-code change.

The gas bath is an explicit external reservoir. Its signed heat appears in
the environmental report and is already counted once in `external_heat_j`.
Finite-wall heat stays internal and molecular heating keeps its separate
account. A late source, thermal, flow, budget or cancellation error publishes
neither liquid nor wall state. This is a continuum, quasi-steady Fourier-film
closure and the existing first-order thermal/flow split; it does not derive
the film thickness, resolve gas convection, evolve finite gas storage, or
supply evaporation, radiation, rarefied-gas transfer or gas pressure forces.
Focused tests in `tests/material_free_surface_enthalpy/conduction/environment.rs`
compare independent layer/storage calculations, data-driven native population
changes, combined heat accounting, stiff exchange and transactional refusal.

### Surface radiation on the native free-surface thermal timeline

`step_with_radiative_environment` adds gray diffuse radiation on
explicitly exposed liquid/gas patches. The nonlinear Stefan–Boltzmann law is
reevaluated every thermal substep, using a factored fourth-power difference
to preserve fluxes near equilibrium. Its positive secant conductance joins
gas-film and internal-contact conductances in the same monotonicity bound;
no temperature or energy clipping supplies stability. Radiant and gas-film
heat may share a physical surface and are each counted exactly once. Finite
walls, source-bound interfaces, optional molecular heat and attached electrical
sources retain their original owners. The resulting enthalpy feeds the actual
material viscosity and native VOF momentum/transport transaction.

`FreeSurfaceEmissivity::Prescribed` retains a declared constant in [0,1].
`Material` consumes `hemispherical-total-emissivity` already selected in the
original complete `LiquidTemperatureLaw` bundle; it never selects another claim
or falls back to a supplied constant. Source values must be static, dimensionless
and in [0,1]. A source-known zero retains a reflecting boundary, independently
of absent data. The original pressure/composition/finish context and source receipt
travel with each sample. Changing that material DATA changes radiant heat, finite
enthalpy, viscosity and the native momentum trajectory without solver changes.

This is a gray opaque surface with view factor one to black isothermal
surroundings, and cell temperature used as
surface temperature. It needs a resolved thermal boundary layer or a justified
lumped-cell approximation, plus time refinement for the explicit first-order
split. It does not infer emissivity, visibility, multiple reflections, spectral
transparency, gas radiation, evaporation or finite surrounding storage. Source,
geometry, budget and cancellation refusals publish no partial state. Focused
G1/G3/G4 tests compare the analytic nonlinear cooling integral and its time
refinement, combined gas/radiant energy, native viscosity/momentum response,
stiff subcycling and atomic refusal on synthetic constitutive fixtures.
The shared solved-interface path samples the same pinned emissivity at its
actual surface temperature, independently of bulk liquid properties. Its
nonlinear residual, finite-gas opposite heat and optional implicit enthalpy
evaluations use that value, retaining their original conservation/rollback
boundaries. Source bounds narrow the surface bracket; unsupported trials refuse.
Temperature-dependent emissivity does not imply a spectral radiation model.

### Condensation/evaporation inventory on existing native liquid interfaces

`MaterialFreeSurfaceEnthalpy3::exchange_interface_mass` accepts integrated
signed mass, carried liquid specific enthalpy and separately partitioned heat.
It changes the actual original VOF interface fill and extensive liquid enthalpy
in one transaction. Positive mass grows liquid; negative mass removes it. The
external phase-change owner supplies the flux and latent-heat partition and
explicitly maps any gas formation reference to the native liquid caloric datum.
The actual represented transfer is reported separately from the request, and
the external owner must debit that represented amount within its own balance
budget. Latent heat is never inferred a second time from the carried enthalpy.

Material crosses at the accepted local liquid velocity. Its carried momentum
is independently checked against actual endpoint momentum after the original
buoyancy/thermocapillary force frame is refreshed. No population-equilibrium
reset, recoil impulse, numerical dissipation heat or hidden mass carry is added.
Energy is audited using centered local increments, so a large common caloric
datum cannot provide artificial tolerance for lost thermal exchange. Every
source, geometry, mass, energy, momentum and cancellation check precedes state
publication. Source viscosity, rheology and supported capillary geometry remain
owned by the original liquid; its accepted clock is unchanged for the outer
coupling owner to advance the normal thermal/flow interval exactly once.

The bounded rung requires occupied native interface cells with positive final
fill at most one. Nucleation in gas, complete interface exhaustion, overfill,
relative gas/liquid velocity and phase-front remapping are outside this API.
An active solved-surface stress or dynamic contact history requires its complete
coupling owner and refuses this instantaneous update. This is a conservative
physical ingress seam, not an independently inferred condensation rate or gas
energy model. G0/G1/G3/G4 regressions in `material_free_surface_enthalpy/mass_exchange.rs`
cover real fill growth/recession, caloric-source refresh, carried momentum,
subsequent native mass/enthalpy transport, reference shifts and atomic refusal.
The shared condensation regression exports the original wet-gas proposal's film
into the native liquid exactly once, retaining only its gas endpoint. It checks
combined gas/liquid mass, maps IAPWS and native caloric datums explicitly, and
reconciles the phase owner's internal energy with native enthalpy using the
reported saturation-pressure volume work. The isothermal phase owner's latent
heat duty is retained rather than deposited a second time into the liquid.

`project_interface_mass_exchange` selects a reproducible increment from the
actual tracked native mass before an external gas caloric solve; it does not
recover mass by multiplying a rounded fill fraction. The original pure-liquid
caloric inverse and bounded face-adjacent gas-patch admission are available to
that coupling owner. `fs-couple::liquid_gas::evaporation` uses these read-only
operations, the shared finite-rate phase owner and this native exchange to
publish finite gas plus real VOF liquid together on one native flow interval.

### Kinetic account for co-moving native mass exchange

`MaterialFreeSurfaceEnthalpy3::exchange_interface_mass_with_kinetic` applies the original local-velocity signed exchange and returns its actual resolved kinetic-energy change, carried kinetic energy and local enthalpy-plus-kinetic defect. Positive, negative and zero-mass heat-only updates use the same source, geometry, material refresh and publication rules. It adds no mixing heat or powered recoil. The existing `exchange_interface_mass` API and arithmetic are unchanged. Finite donor/receiver owners can close a shared mechanical account without reconstructing kinetic energy from pressure density or losing force-refresh changes.

### Positive liquid inflow with explicit local complete mixing

`MaterialFreeSurfaceEnthalpy3::mix_interface_inflow` adds positive mass at a
declared incoming velocity to an existing occupied interface. The caller selects
a perfectly inelastic local mixing law: with actual represented incoming mass
`dm`, old mass `m` and final mass `M`, velocity becomes
`u + (dm/M)*(v_in-u)` and `0.5*dm*(m/M)*|v_in-u|^2` becomes liquid enthalpy.
This is the authored loss of relative kinetic energy under complete mixing;
numerical remap defects are never converted to heat. Incoming enthalpy must
invert to pure liquid on the SAME original caloric chart; its density need
not equal receiver storage before mixing. The accepted mixture still passes
the original phase, source and reference-density gates. An unsupported final
thermal density/volume change refuses the whole transaction. Positive
throughflow streams use this same incoming admission; outgoing restrictions
remain unchanged. Every accepted mass transaction invalidates a retained
planar thermal-source continuation binding. The incoming
and accepted velocities must both fit the original scenario's speed domain.

The native update retains canonical population density exactly, preserving
pressure and the mass-to-fill relation. It translates the physical equilibrium
in the actual Guo force frame while retaining physical nonequilibrium stress.
Original material, force, rheology and surface owners re-admit the complete
candidate before publication. The report separates external heat, mixing heat,
carried kinetic energy and actual kinetic change; independent local mass,
momentum, enthalpy and enthalpy-plus-kinetic defects must all fit their budgets.
The co-moving signed exchange API and its arithmetic remain unchanged.

Several streams can enter the same cell in one atomic native update. Each
stream independently fits the original liquid caloric chart and speed bound;
averaging cannot admit an invalid source or imply a new composition law.
Mass-weighted enthalpy and velocity retain a central kinetic-energy variance
`q = sum(m_i*|v_i-v_mean|^2)/(2*sum(m_i))`. Represented incoming kinetic energy
includes `dm*q`, and the same amount becomes additional local mixing heat
exactly once. Opposing jets therefore retain their energy even at zero mean
velocity. Stable relative-velocity reduction follows the authored stream order;
permutations agree to floating-point roundoff. The existing stream-count bound
continues to bound work.

Each shared cell yields one sample with its summed requested and actual
represented mass. All its source owners use that sample's common
represented/requested ratio to debit their own mass, momentum and carried
energy. Separately supplied heat is summed without mass scaling. Single-stream
and co-moving arithmetic are retained, with no intermediate material refresh
between co-located streams. Focused numerical tests cover independent stream
budgets, boosts, density/datum shifts, permutations and overflow refusal; the
public native fixture covers opposing streams and whole-transaction refusal.

This bounded filling boundary retains positive interface capacity and the
original liquid density/phase assumptions. It does not resolve jet impact,
nucleate cells, supply arbitrary outflow recoil, advance the clock or establish
global momentum/kinetic conservation of VOF conversion. Native regressions cover
unequal density, force-frame stress and carried momentum; public material tests
exercise filling, thermal storage, normal subsequent transport, co-moving
identity and transactional refusal.

### Powered liquid outflow with reciprocal recoil

`MaterialFreeSurfaceEnthalpy3::eject_interface_outflow` removes actual native
interface mass at a prescribed outlet velocity through an explicitly selected
ideal actuator. For removed mass `r`, old mass `m`, retained mass `M=m-r`, old
velocity `u` and outlet velocity `v`, the remaining velocity is
`u' = u - (r/M)*(v-u)`. Liquid plus stream retain their linear momentum. The
external actuator supplies `W = 0.5*r*(m/M)*|v-u|^2` with zero net impulse.
The caller must authorize a finite nonnegative batch work allowance and debit
the returned actual work from its source. Work exceeding that allowance refuses
the complete transaction. Required work uses stable relative velocities; actual
native readback independently closes enthalpy plus kinetic energy against signed
stream enthalpy, separately supplied heat, signed stream kinetic energy and W.

Mechanical work never becomes heat or a negative mixing loss. The outgoing
enthalpy itself must fit the same original pure-liquid fixed-density caloric
chart. Outlet and accepted liquid speeds must both satisfy the original physical
scenario. Canonical pressure density, native mass-to-volume conversion, physical
nonequilibrium stress, source refresh and transactional publication use the same
owners as inflow. Co-moving outflow requires exactly zero work and retains the
original signed-exchange arithmetic. Exhaustion, excess recoil, unrepresented
mass, coupled contact/surface restrictions and source-domain exits still refuse.
This boundary does not infer evaporation recoil, nozzle pressure work, gas-volume
displacement or thermodynamic phase expansion, and it does not advance time.

Several outlets may remove mass from the same interface in one native update.
The shared stream reducer uses positive removed-mass magnitudes to retain mean
velocity, mean liquid enthalpy and central specific kinetic variance `q`.
The combined signed transfer is still negative. Its carried kinetic energy
includes `dm*q`, while actuator work includes `r*q` for `r=-dm`; this energy
leaves with the streams and never becomes heat in the remaining liquid.
Opposing outlet velocities therefore require work even when their mean causes
no recoil. Work admission covers the sum across every cell and every outlet.

Each outlet independently satisfies the original liquid caloric and speed
domains before reduction. The native owner projects the total removal once
and returns one sample per cell. Its common represented/requested mass ratio
applies to each outlet's credited mass, momentum and carried energy. Authored
heat sums without mass scaling. Positive remaining inventory and accepted
recoil limits remain binding, and every refusal leaves the accepted owner
unchanged. The original single-outlet arithmetic is retained; co-moving
exchange and inflow use their existing laws. Focused tests check opposing
outlets, all three-source permutations, frame/datum/density covariance,
represented fractions, dispersion-work overflow and single-stream identity;
public fixtures cover immediate/prepared publication and shared-cell refusals.

### Prepared interface mass transfers

`prepare_interface_inflow` and `prepare_interface_outflow` return typed handles
that exclusively borrow the original thermal-flow owner and retain one accepted
native lattice/enthalpy candidate. `report()` supports joint multi-owner balance
checks; dropping a handle discards its candidate. Consuming `publish()` moves
only the already accepted native fields into their original owner, with no new
allocation, callbacks, cancellation poll or fallible work. Coupled callers must
perform their final cancellation poll before publishing any participant. The
exclusive borrow prevents stale-state publication; immediate inflow/outflow
entries use the same preparation and publication arithmetic.

### Finite atmosphere on the native free-surface thermal timeline

`step_with_finite_gas_environment` gives the same moving-liquid thermal owner
an actual `IsochoricGasParcel`. Each accepted thermal substep transfers equal
and opposite film heat between the liquid and the original NASA-9 gas internal
energy. Fixed gas mass, composition, volume and formation reference are retained;
its new temperature and EOS pressure re-evaluate the same source-backed films.
Bulk and film transport claims are checked again at the actual endpoint.

The shared explicit step bound includes the sum of all exposed conductances
against the finite gas capacity, together with the original liquid, internal
conduction and finite-wall contact rows. A small atmosphere therefore cannot
silently act as an infinite heat reservoir. Source-dependent liquid viscosity
and rheology, ordered wall interfaces, optional molecular heat and native VOF
mass/enthalpy transport continue through their original owners. All proposed
gas, finite-wall and liquid states publish only after the native flow and
advection gates succeed.

`FreeSurfaceConductionReport::finite_gas` records actual gas internal-energy
change, reciprocal heat defect and the combined liquid/wall/gas thermal
account. Environmental heat appears once as a supply to the liquid subsystem
and cancels against gas storage in their union. Authored heat and the existing
physical molecular-heating account remain distinct; no full VOF kinetic-energy
or pressure-work closure is inferred.
The gas storage audit integrates the original frozen-composition Cv law between
the actual initial and final temperatures, independently of the accumulated heat
receipts. It does not subtract large absolute NASA9 energies. Changing a common
formation-energy reference therefore cannot reject an otherwise admitted liquid,
wall and atmosphere step or alter its heat, pressure or native liquid motion.

The caller supplies currently exposed geometry at each interval; the film is
held during the thermal stage. Gas volume remains explicitly fixed. These
thermal entries do not themselves perform reaction, phase transfer, gas-volume
work or a phase-density handoff. Original attached surface/contact-history
admission and first-order thermal/flow limits remain binding. Five focused
synthetic regressions cover independent nonlinear NASA energy/pressure limits,
data-driven native population changes, finite wall/interface/molecular heat,
summed small-gas-capacity timestep admission, and late source/cancellation
rollback. Their addition does not establish executed native validation.

`step_with_radiative_finite_gas_environment` connects the existing liquid
radiation input to this same finite-atmosphere kernel. All old finite-gas
entries delegate with empty radiation. Radiation uses the original factored
fourth-power law, source or prescribed emissivity, exposed geometry and combined
thermal timestep bound. It is external supply to liquid; a gas film's opposite
heat remains internal. Radiation-only exposure is admitted without a gas heat
debit. The original finite-gas caloric account is still present.

An optional density-excursion allowance selects the existing post-thermal gas
pressure load, retaining its fixed reference and stationary geometry checks.
The accepted EOS, liquid, finite solids and reports publish atomically after
one native flow step. This is the existing cell-temperature gray-radiation
approximation, requiring resolved boundary layers or justified lumped cells;
it does not add solved-interface latent feedback or relax wetting restrictions.

## Public types and semantics

- D2Q9 constants: `Q` (9), `CS2` (`1/3`), `MACH_LIMIT` (`0.3`).
- `equilibrium(rho, ux, uy) -> [f64; 9]` — the D2Q9 equilibrium distribution
  (recovers `Σf = ρ`, `Σeᵢfᵢ = ρu`).
- `Lbm::channel(nx, ny, tau, gx)` — a body-force-driven channel (periodic x,
  bounce-back y-walls). `step`/`run` (collide + Guo forcing + stream +
  bounce-back); `density`, `velocity` (Guo-corrected), `total_mass`,
  `viscosity` (`ν = (τ−½)/3`), `x_velocity_profile`.
- `core2::CollisionModel2` (bead 3zkcr) — the D2Q9 collision-operator
  selector on `Grid::collision`: `Bgk` (default; byte-for-byte the
  historical arithmetic, so every existing pin/golden is unchanged),
  `Regularized` (Latt–Chopard projected), `RecursiveRegularized`
  (Coreixas-style RR3), `Trt { magic }`, and `CentralMoment`
  (cascaded, the D2Q9 analog of `CollisionModel3::CentralMoment`).
  All modes keep the exact shear-viscosity law `nu = (tau - 1/2)/3`
  (executed shear-wave identity, all under 0.1% deviation), conserve
  mass/momentum exactly, and are bitwise deterministic. EXECUTED
  stability landscape (2026-08-10, `examples/stab_scan` and the
  fs-aeroac jet-labium probe): `CentralMoment` is the robust
  high-cell-Reynolds operator (rig runs through Re_cell = 192 where
  BGK dies at ~48); plain `Regularized` destabilizes on unresolved
  (~1-cell) shear layers where BGK survives, `RecursiveRegularized`
  repairs that but hits a temporal-Nyquist parity artifact with
  bounce-back at high shear, and every `Trt` magic above the BGK
  point REDUCED the stable window at tau near 1/2 (the "magic = 1/4
  optimal" folklore is false on this family). Body forces in the regularized
  and central-moment modes act through the physical `f + S/2` frame and
  return `C(f + S/2) + S/2`; TRT splits the Guo source into even/odd parts
  with their corresponding relaxation rates. All modes deliver the imposed
  momentum increment independently of viscosity. Forced uniform-flow stress
  and laminar Poiseuille profiles are checked in the collision battery.
  Combining non-BGK modes with the partial-saturation or free-surface
  operators remains mechanically possible but unvalidated.
- `core2::Grid` / re-exported `Grid` — a general dense D2Q9 grid with cell
  flags, vector gravity, per-cell relaxation time, per-cell external force,
  periodicity flags, deterministic collide/stream steps, and gas/wall
  boundary bounce handling for the plain non-free-surface step.
  `Grid::rebind_gravity` and `Grid::rebind_forcing` replace held acceleration
  and local loads without an instantaneous physical half-force impulse, using
  the caller's reusable collision scratch and atomic accepted-state publication.
  `core2::MomentumExchange2`, `stream_from_with_wall_momentum`, and
  `step_with_wall_momentum` add an opt-in raw lattice-impulse receipt for an
  explicitly selected subset of stationary wall cells without changing the
  legacy step path. `core2::MovingWallMomentumExchange2` plus the paired
  moving-wall stream/step methods apply per-wall-cell halfway-bounce velocity
  corrections and return boundary-relative impulse, torque, work, resolved
  population momentum, and moving-mass transfer terms for an explicitly
  selected wall subset. `core2::CurvedMovingWallLink2` plus the paired curved
  moving-wall stream/step methods consume an exhaustive per-fluid-cell,
  per-direction wall-link map and apply the linear Bouzidi-Firdaouss-Lallemand
  off-lattice interpolation. Link-local velocity and the declared wall
  fraction determine both the incoming population and the exact intersection
  used for torque. `core2::WallTopologyTransition2` and
  `Grid::transition_wall_topology` atomically replace a fluid/wall cell mask,
  initialize newly uncovered cells from unique surviving one-ring donors, and
  receipt the exact active-population mass/momentum removal and insertion.
  `core2::PartialSaturationCell2`, `core2::PartialSaturationExchange2`, and
  the paired Grid collision/step methods add a BGK/Guo dense-D2Q9
  Noble-Torczynski partially saturated cell rung. Caller-supplied solid
  fractions and velocities blend fluid and non-equilibrium-bounce-back solid
  collisions while the receipt retains aggregate equal-and-opposite
  fluid/solid impulse, angular impulse, solid work, and coupling-weight sum.
  `core2::PartialSaturationForcing2` separately accounts the complete imposed
  fluid impulse, torque and midpoint work, including occupied cells.
  `core2::PartialSaturationBody2`, `core2::PartialSaturationBodyExchange2`,
  and the paired finite-body collision/step methods couple that same operator
  to finite planar translational and rotational inertia at the midpoint.
  `core2::UnderResolvedCylinderPair2` and
  `core2::UnderResolvedLubricationExchange2` add a separate pair-local,
  per-unit-depth normal-resistance correction for smooth circular cylinders
  below a caller-declared cutoff. The immutable descriptor owns checked gap,
  reduced radius, viscosity, common normal, line-of-action point, and surface
  velocities; its receipt returns exact equal-and-opposite solid impulse and
  angular impulse, dissipative pair work, the effective gap, and the applied
  resistance increment without mutating the lattice or a rigid-body state.
  `core2::VelocityPressureX2` plus the paired Grid step methods provide a
  low-Mach regularized velocity inlet at x-min and density outlet at x-max with
  periodic y closure. Density/velocity and physical Guo nonequilibrium stress
  are extrapolated from the first-interior columns. Held gravity and local
  force densities use `F = rho*g + fext` at each source and destination;
  reconstruction restores the destination half-source after copying
  `Pi_neq + (u F + F u)/2`. Original collision applies force once and original
  streaming owns the selected-wall impulse. All face candidates are admitted
  before publication using O(ny) reusable scratch. Scratch retains the original
  collision-field extent on success or refusal.
  `core2::PressureX2` and `Grid::step_pressure_x` prescribe pressure densities
  at both x faces, so the pressure head determines the flow in either direction.
  Both faces use the same first-interior physical Guo velocity/stress and
  destination half-source reconstruction; neither imposes an inlet velocity.
  The measured variant retains the original stationary-wall impulse, separate
  from open-face transport. Face densities represent isothermal lattice
  pressure `p = CS2*rho`; SI inputs use the material pressure-perturbation
  conversion, not absolute thermodynamic pressure. The native battery checks
  uniform pressure stationarity, barotropic hydrostatic gravity balance,
  reversed pressure-head flow around an obstacle, and exact discrete momentum
  balance between open-face transport and wall drag. Prescribed-density force
  overflow is refused before collision; late source refusal remains atomic.
- `rheology::Rheology`, `rheology::update_tau`, and
  `rheology::channel_flow` — local apparent-viscosity laws with a shared,
  self-consistent constitutive solve in D2Q9 and D3Q19 `update_tau3`.
  Both dimensions infer strain from `Pi_neq + (u F + F u)/2` in the
  physical Guo frame. D2Q9 includes both `rho*g` and the local external
  force density; uniform acceleration is not shear. The shared estimators
  correct three/six tensor entries directly, avoiding source-array
  construction and an extra population pass; material-driven rheology paths
  reuse these estimators.
  The population stress supplies `A = gamma*tau`; a safeguarded bracketed
  Newton solve closes `tau = 1/2 + nu(A/tau)/CS2` for the same local state.
  The previous tau is only a warm start. Unclipped solutions require a
  relative tau defect at most `1e-12`, with a 96-iteration bound and an
  explicit failure if convergence is exhausted. Constant Newtonian and
  constant-limit power-law/Carreau laws are evaluated once per grid update,
  avoiding per-cell iteration. The existing power-law shear regularization
  at `gamma = 1e-12` remains explicit in the evaluator.
  `TAU_FLOOR = 0.505` and `TAU_CAP = 3.0` bound representable relaxation;
  activating either changes the requested constitutive law and is counted
  in `RheologyStats`. The reported shear rate uses the applied tau. A local
  constitutive solve does not make the collide/stream time integration
  implicit. Force-reframing regressions cover zero/full strain, density and
  relaxation changes, and independent force directions; analytic channel
  gates remain the flow-level checks.
- `thermal::ThermalLbm` and `thermal::gbeta_for_rayleigh` — D2Q9 flow plus
  D2Q5 temperature populations for Rayleigh-Bénard-style onset fixtures with
  fixed-temperature wall rows.
  D2Q9/D3Q19 slab forces are held over each flow interval. Before replacing
  buoyancy plus an external (including Lorentz) load, the existing populations
  are rebound by half the difference of the raw Guo sources, at the accepted
  physical velocity. Mass, physical momentum and the force-corrected stress
  stay unchanged by this coordinate change; only the ensuing native collision
  supplies the physical impulse. Unchanged force leaves population bits alone.
  This prevents a spurious extra half impulse when a magnetic field is enabled,
  reversed or switched off, including material/midpoint/local-transport paths.
  The operation adds no heat, work, solver or physical state. G0 tests preserve
  all stress components through arbitrary source changes; G1 tests check
  `rho*delta_u=F*dt` through switch-on/off/reversal and independent SI Lorentz
  acceleration in both dimensions. These checks do not supply a total-energy,
  unconditional-stability or temporal-order certificate for the split solver.
  Onset perturbations shift the equilibrium velocity in the physical source
  frame and retain the existing nonequilibrium stress, density and transverse
  velocity; zero amplitude leaves population bits unchanged. Thermal collision
  uses the input physical velocity and the original momentum scratch buffer,
  before the ensuing flow advance overwrites that scratch. Each slab retains
  one advective-flux vector per cell for the scalar correction described below.
- `MaterialSlab::step_with_material_rheology` (2-D and 3-D) resolves shared
  material transport at actual temperature and local physical shear rate.
  The existing non-equilibrium stress estimators include the Guo half-force
  second-moment correction; a bounded local root solve closes
  `tau=1/2+3*nu(T,rate)*dt/dx^2` with `rate=numerator/(tau*dt)`.
  `RheologySolve::relaxation` controls the original damped search until two
  consecutive material evaluations give opposite tau-residual signs. A
  residual-scaled secant then stays strictly inside that sampled tau bracket;
  bisection takes over if the interval has not halved since the preceding
  trial or rounding prevents an interior secant. No global tau endpoints are
  probed, and each query retains the original source-domain admission.
  `maximum_iterations` counts every material evaluation, with cancellation
  checked before each. Acceptance requires the actual returned material
  query's tau defect to meet the original absolute/relative tolerances;
  bracket width alone cannot accept a discontinuous response. Exhausting the
  query budget or the representable bracket interior refuses convergence.
  The retained background relaxation field is only an initial iterate, never
  an assumed viscosity for the evolving material. Both initial and candidate endpoint
  source evaluations must converge before the existing thermal/flow step is
  published. Original source payloads, physical queries, iteration counts and
  residuals remain on `transport::rheology::RheologyMaterial`. Source failures,
  exhausted iteration and observed cancellation preserve all accepted state.
  No viscosity floor/cap or material-name law is substituted. This is a
  fixed-density/Cp, low-Mach generalized-Newtonian BGK path with explicit
  initial-state transport, not a viscoelastic model, implicit thermal coupling,
  viscous-heating account or arbitrary-law convergence guarantee.
  `tests/material_rheology.rs` checks independent prescribed stress moments
  with forcing, both-dimensional constant-material equivalence, a source-card
  viscosity curve against its analytic constitutive root, actual material-swap
  momentum response, and transactional refusals.
  Source-authored `Surface2D` laws additionally carry simultaneous temperature
  and shear-rate dependence through the same complete material bundle. The
  focused manufactured affine law checks the independent constitutive quadratic,
  source pins and actual momentum changes under data-only substitutions. With
  viscous heating enabled, the identical evaluated viscosity supplies both
  momentum and physical heat, and the heated endpoint changes the next flow
  step. These fixtures establish plumbing and numerical identities, not
  measured non-Newtonian material accuracy.
- `MaterialSlab::step_with_material_rheology_and_viscous_heating` (2-D and
  3-D) explicitly adds molecular `2 mu D:D` from those same cached
  initial material/shear samples to the existing source/transport step.
  Authored heat must exclude this molecular loss. Per-cell physical power,
  integrated viscous heat, authored heat and their checked source partition
  accompany the original material receipts and independently measured
  source/storage/wall balance. Endpoint material evaluation sees the heated
  temperature, and repeated calls feed it back into viscosity and flow.
  No second constitutive query or parallel thermal solver is introduced.
  All conversion, temperature, accounting, endpoint and cancellation refusals
  precede publication. This remains fixed-density/Cp, incompressible,
  generalized-Newtonian, low-Mach BGK with explicitly frozen initial loss;
  numerical momentum dissipation, pressure work, independently sourced bulk
  heating and elastic storage are excluded. Tied BGK/TRT bulk work is retained.
  A sensible-heat balance is not
  a transient total-energy conservation proof or an implicit-time claim.
  The material rheology tests compare both dimensions with independently
  prescribed shear power, check the exact zero-stress limit, exercise actual
  source-card temperature/viscosity feedback, and cover atomic refusals.
- `MaterialSlab::step_with_boussinesq_material_transport` (2-D and 3-D) —
  evaluates actual local material density for `(rho(T)-rho0)*gravity`, and
  the same evaluation's dynamic viscosity and conductivity for
  `nu=mu/rho0`, `alpha=k/(rho0*cp0)`. The shared material sampler, thermal
  source step, lattice kernels, independent wall heat accounting and endpoint
  checks own the advance. Original source payloads and unmodified transport
  remain on `transport::boussinesq::BoussinesqMaterial`; its explicit projected
  transport retains the original reference density and heat capacity.
  `BoussinesqConfig` requires finite signed gravity and caller-selected relative
  density/Cp departure bounds in `[0,1)`, enforced at initial and proposed
  endpoint states. No mass, heat storage or physical units are rebound.
  Density-resolved force replaces the independent background beta load for
  this step; prescribed external forces remain additive. The background beta
  and relaxation fields are retained for other explicit solver entry points.
  Any refusal or observed cancellation preserves the entire accepted state.
  These are explicit initial-state coefficients, with support checked again
  at the endpoint; the numerical bounds do not establish physical validity.
  The caller must justify small neglected density/Cp departures and timestep
  refinement. No compressible work, phase/species evolution or gravitational
  energy account is introduced. G1/G3 ideal-gas pressure/temperature checks,
  three-dimensional source/force ownership and G0/G4 refusal/retry coverage
  live in `tests/material_transport_feedback.rs`.
- `MaterialSlab::step_with_boussinesq_material_rheology_and_viscous_heating`
  composes those existing source owners in both dimensions. A single converged
  material response at actual temperature/shear rate supplies density-derived
  buoyancy, molecular viscosity/conductivity and `mu * rate^2` heat. The shared
  sampler, heat partition, tentative transport step and endpoint admission
  publish together. Original source receipts and unprojected density survive;
  reference density/Cp still own mass and sensible storage, and the independent
  background beta load is replaced for the step. The same explicit Boussinesq,
  generalized-Newtonian and sensible-energy limits above apply; this does not
  claim total mechanical/gravitational energy conservation. G1/G3 tests in
  `tests/material_rheology.rs` compare both dimensions against independently
  prescribed force and stress power, change density/shear source data in the
  actual repeated flow, and verify initial/endpoint refusal and cancellation
  rollback.
- `freesurface::FreeSurface`, `ContactModel`, `dam_break`, and `surge_front`
  — dense-grid VOF-style mass tracking with interface/gas/fluid conversion
  bookkeeping, conservative neighbor-local redistribution, contact-model bracketing,
  and qualitative dam-break / jet-fragment fixtures.
- `freesurface::query_dynamic_wetting2`, `DynamicWettingAnswer2`,
  `ContactLineRegime2`, and `DynamicWettingError2` — a read-only adapter from
  an ordered, history-bearing `InterfaceSystemCard` to a receipt-bearing
  receding/advancing contact-angle bracket. Both speed inputs are `QtyAny`
  values carrying coherent-SI velocity dimensions. It queries both endpoints
  at the absolute contact-line speed, selects by signed speed outside a
  symmetric pinning deadband, and returns no selected angle while pinned. The
  answer retains both `MaterialAnswer` receipts and the interface-system
  content identity.
- `d3q19::{Duct, equilibrium3, duct_analytic}` — the frozen D3Q19 BGK + Guo
  body-force core on aligned 4x4x4 SoA tiles, with stationary halfway
  bounce-back x/y walls and periodic z for the rectangular-duct fixture. Its
  axial-z BGK collision dispatches once to a bit-neutral NEON/AVX2 tile capsule
  or the scalar twin. Its separate pull-stream dispatcher moves four-cell x
  rows through a bit-neutral NEON/AVX2 capsule or the scalar stencil; SIMD
  lanes are independent cells, and streaming performs no floating arithmetic.
- `d3q19_bgk_simd_tier` / `D3q19BgkSimdTier` — truthful operation-specific
  receipt for the Duct collision kernel. AVX-512 hosts report AVX2 because this
  capsule deliberately uses four-lane AVX2; x86 without AVX2+FMA reports scalar.
- `d3q19_stream_simd_tier` / `D3q19StreamSimdTier` — separate receipt for the
  Duct pull-stream kernel. AArch64 reports NEON; x86 requires only AVX2 for this
  pure-move capsule and otherwise reports scalar.
- `d3q19::{CollisionModel3, collide_cell3, CollisionError3}` — the checked
  per-cell collision authority used by `BoundaryGrid3` and public cell calls.
  `Bgk { tau }` retains the frozen general-force projection; the Duct tile
  scalar/SIMD twin is separately bit-locked to its test-only axial oracle.
  `Trt { tau, magic }` relaxes even populations with `1/tau` and odd
  populations with `1/(1/2 + magic/(tau-1/2))`. Each Guo source parity uses
  its own `1-omega/2` factor. Opposite pairs share their increments, with no
  heap allocation or moment transform. The viscosity stays `(tau-1/2)/3`;
  `magic=(tau-1/2)^2` recovers BGK to roundoff. Positive finite magic and
  representable even/odd relaxation times are admitted before collision.
  `magic=3/16` removes the viscosity-dependent wall offset for uniformly
  forced, fully developed planar Couette-Poiseuille flow; arbitrary forcing,
  curved/voxel geometry and open faces retain their own error models.
  The same model is available through `BoundaryGrid3`, sparse D3Q19 and
  `Grid3<D3q19>`, and via `LbmCollisionChoice::Trt { magic }` on the original
  physical duct-to-conjugate-heat path. Material/rheology/SGS retuning keeps
  magic fixed while updating even viscosity and its corresponding odd rate.
  Focused regressions cover the three-axis planar profile across viscosities,
  a BGK wall-slip control, force/stress moments, invalid parameters, physical
  duct selection, source retuning and portable/FMA equivalence. This adds no
  turbulent-transition or general high-Reynolds stability claim.
  `CentralMoment` relaxes a full-rank centered D3Q19 monomial basis with
  independent second/higher-order rates. `ReducedCumulant` applies the
  nonlinear `C220`/`C202`/`C022` corrections to the three independent
  fourth-order rows represented by that basis. Both moment families admit
  uniform and per-cell body forces through `C(f + S/2) + S/2`, using the
  original Guo source `S` and force-corrected velocity. Central moments thus
  receive `(I - Lambda/2) S` with each mode's own rate; reduced cumulants
  are evaluated in the physical frame with zero first central moments.
  **`ReducedCumulant` is linearly unstable at rest for any field that varies
  along all three axes** (pexfp.138). With third- and fourth-order rates 1
  and `second_order_rate` 1.98807 (nu = 1e-3), with no wall and no force, a
  seed grows at +4.09e-2 per step on a periodic 4 x 4 x 4 box, at
  k = (pi, pi, pi/2), and at +5.11e-2 on 20 x 4 x 4, at k = (3pi/5, pi, pi).
  From 1e-8 it reaches negative density within about 420 steps. At the
  rates checked, fields uniform along any one axis don't grow.
  - **Cause.** At rest the relaxed rows are Hermite polynomials. On D3Q19,
    H220 and H002, and their permutations, have weighted inner product
    -1/27: these are the eight corner velocities D3Q19 lacks. It is the only
    overlap between rows relaxed at different rates. Without it the linear
    step can't grow in the 1/w norm, at any wavenumber or any rates in (0, 2).
  - **`CentralMoment` at `higher_order_rate` 1 also grows at rest** once
    `second_order_rate` is above about 1.988, more slowly and at longer
    wavelengths. It grows at +8.31e-4 per step at 1.98807 on 12 x 4 x 16.
  - **`CollisionModel27::Cumulant` has no such overlap and shows no growth**
    at rest at any wavenumber. This was checked with omega
    [r, r, 1 x 8] and with ten unlike rates, and on two periodic boxes in
    Rust.
  - **Not covered:** shear, walls and mean flow.
  - `tests/rest_wave_three_axes.rs` checks this: a three-axis density wave
    grows on `ReducedCumulant` and decays on `Cumulant`, and its one-axis
    twin decays on both.
  All paths reject invalid parameters or state before publication. On an
  x86-64 CPU with FMA, every cell runs the same body compiled under
  `target_feature(enable = "fma")` (the `d3q19/fma` capsule) instead of one
  `fma()` call per `mul_add`. `mul_add` rounds once on both paths, so the
  results and refusals are bitwise those of the portable path.
- `D3Q19_MOMENT_COLLISION_SEMANTICS_VERSION` — independent version for the
  optional moment-space basis, rate grouping, `mul_add` chains, cumulant
  projection, and deterministic back-transform. Version 3 uses binomial
  central-to-raw shifts and the fixed D3Q19 inverse instead of per-cell dense
  factorization. This changes moment-collision roundoff and the unregistered
  reduced-cumulant replay candidate, including unforced cells. The registered
  BGK surfaces retain version 1 and their existing golden values.
- `Face3`, `FaceBoundary3`, and `BoundarySpec3` — six-face axis-aligned
  boundary declarations: paired periodic faces, tangential moving/stationary
  halfway walls, regularized velocity faces, and isothermal pressure/density
  faces. Open faces currently require zero body force.
- `LinkMaskTile3` and `BoundaryLink3` — aligned per-tile D3Q19 effective-wall
  masks and canonical `(tile, lane, direction)` enumeration. Separate open
  masks distinguish populations reconstructed after streaming.
- `BoundaryGrid3` — tile-major D3Q19 grid whose compiled masks drive the
  runtime pull stencil. The legacy `new(..., tau, ...)` constructor remains
  BGK, while `with_collision_model` selects an explicit checked model for the
  whole grid. `voxelize_sdf` samples a scalar field at cell centers, commits
  solid occupancy/masks atomically, and then locks topology. Solid and planar
  walls use link-wise halfway bounce-back. Moving walls add the standard
  `2 w_i rho (c_i dot u_wall) / c_s^2` incoming-link correction.
  Velocity/pressure face-interior cells use second-order Hermite regularized
  non-equilibrium-stress reconstruction from the first interior cell.
- `WallMomentumExchange3` and `BoundaryGrid3::{step_with_wall_momentum,
  step_pooled_with_wall_momentum}` supply hydrodynamic impulse, three-axis
  angular impulse, and work on any canonical subset of the existing wall
  links. The caller can reuse the selection across steps. Loads read the same
  post-collision populations and halfway-bounce expression as streaming;
  the default path's population arithmetic is unchanged. Measurement costs
  O(selected links), allocates no per-step buffer, and sums density only once
  per selected boundary cell. Quantities use lattice units, with cell centers
  at `(x+0.5,y+0.5,z+0.5)` and torque evaluated at the halfway intersection.
- `PeriodicMomentumTransport3` closes the periodic-cell angular account for
  single/multiple immersed bodies and their moving/material/lubricated wrappers.
  It measures the actual post-collision populations' coordinate jumps in
  O(periodic boundary area), without allocation or changed population arithmetic.
  The fluid/body balance includes its signed `angular_momentum_change` alongside
  external and wall moments. The value is independent of moment origin; diagonal
  double wraps retain both jumps, while wall-owned links are excluded. One
  assembly-level receipt owns the transport, never each body's nested response.
  Periodicity preserves linear momentum, but does not imply conservation of
  angular momentum inside the represented cell. This receipt is not a body
  torque or an open-boundary flux model.
- `plan_scaling(reynolds, char_length_lu, u_lattice) -> ScalingPlan { tau,
  viscosity, u_lattice, mach, tau_margin, stable }` — the lattice-scaling
  assistant. `ScalingPlan::color()` (verified when comfortably stable, else
  estimated). Panics on non-positive Reynolds / length.
- `poiseuille_analytic(gx, viscosity, ny, y)` — the analytic reference profile.

### Absorbing layer (`sponge`, bead 9ok02)

- `sponge::Sponge2::new(side, width, sigma_max, target_rho, target_u)`
  — a checked absorbing layer over a contiguous x-range of the D2Q9
  `core2::Grid`: after each step, fluid cells in the layer blend
  toward the equilibrium of the authored far-field target in the physical
  Guo frame, `h <- h + sigma (f_eq - h)` with `h = f + S(u,F)/2`,
  then restore the endpoint source at the actual blended density/velocity.
  The held acceleration and spatial external forces stay attached; target
  velocity, momentum and second moments receive no half-force bias. The
  unforced path retains its original population arithmetic, with a quadratic ramp
  `sigma = sigma_max ((d+1)/width)^2` growing toward the boundary (an
  abrupt sigma is its own reflector). `SpongeSide::{LeftX, RightX}`.
  Construction-time asserts (crate boundary convention): positive
  width, `sigma_max` in (0, 1], positive finite target density,
  finite low-Mach target velocity. `apply` skips non-Fluid cells and
  panics if the layer exceeds the grid or a profile mismatches the
  height. `Sponge2::with_profile` takes a per-y-row `(rho, u)`
  profile, turning the layer into a spectral-DNS style FRINGE region
  (periodic-x domains: outflow re-conditioned to an authored inflow
  profile — the jet-labium fixture in fs-aeroac consumes it).
  This imposed reservoir exchanges mass and momentum with the fluid; it is
  not a closed-domain conservative or thermal-energy model. Focused G0/G3
  tests retain forced matched states, independently check mass, physical
  momentum and the full second-moment blend, and pin the unforced arithmetic.
- Purpose: outlet-reflection suppression for aeroacoustic source
  extraction (bead 9ok02 — the jet pilot's 6% Re-200 flux imbalance
  was outlet-reflection sensitivity; source spectra are untrusted
  without a MEASURED reflection coefficient).

### Sparse active-tile infrastructure (`d3q19::sparse`, bead sjro)

- `morton3(tx, ty, tz)` / `demorton3(key)` — 3×21-bit Morton (Z-order) tile
  keys; `morton3` panics past the 21-bit range (a compile-shape bound, not a
  runtime input).
- `SparseGrid3::new(nx, ny, nz, tau, force)` — an EMPTY sparse D3Q19 grid
  over a domain of positive 4-multiples; state exists only for active tiles.
  `activate_tiles(&[(tx,ty,tz)])` is idempotent, atomic (refuses whole-batch
  on any out-of-domain tile), preserves existing tile state bit-exactly, and
  re-sorts the active list into ascending Morton order — iteration, sweep,
  and reduction order depend only on the active SET, never on activation
  history. New tiles start at physical rest in the retained Guo force frame,
  `f = f_eq(rho,0) - S(0,F)/2`; this also gives the native free-surface
  constructor zero initial physical momentum under nonzero forcing.
- `step_serial()` / `step_pooled(&runner, &gate)` — one sweep step as two
  passes (collide from the shared per-cell `collide_cell3` kernel, then
  pull-stream) over three retained buffers. The pull-stream pass writes a
  private transactional destination; only a completely drained successful
  pass swaps that destination into the published state. A link whose source
  lies in an inactive tile or outside the domain bounces back (inactive space
  is solid wall in this increment). `step_pooled` runs both passes through any
  `fs_exec::KernelRunner` and is bitwise-identical to `step_serial` for
  every worker count: kernel tiles own disjoint Morton-contiguous slot
  groups sized by a worker-count-independent constant, and all inputs are
  read-only within a pass.
- `perturb(seed, amplitude)` — splitmix64 hash-seeded density perturbation
  at physical rest in the same force frame (the hash IS the seed schedule);
  `cell_macros` retains its raw population-velocity convention so local-force
  free-surface consumers add their own half source exactly once.
  G0/G1 regressions verify activation, density changes, unchanged existing
  tiles, full resting stress, and three complete force impulses from rest.
  `state_bits()` — canonical-order exact
  bits for golden hashing; `total_mass()`, `cell_macros(slot, lane)`,
  `active_tiles()`, `allocated_state_bytes()` (exactly three population
  buffers per active tile), `state_bytes_per_tile()`.
- `SparseError3` — typed refusals: `Dims`, `TileOutOfDomain`, `Collision`
  (Morton key + lane + the underlying `CollisionError3`), `Cancelled`
  (pre-step state intact; re-issue is deterministic), `Pool`.

#### Full-gradient SGS on active sparse tiles

`SparseGrid3::{sgs_workspace, prepare_sgs, step_serial_with_sgs,
step_pooled_with_sgs}` apply the shared explicit Smagorinsky, WALE or Vreman
closure through the original sparse collision and pull-stream kernels. The
workspace retains only active cells and six neighbor slots per active tile;
SGS preparation performs no step allocation or per-cell sparse-map lookup.
The existing pooled owner's chunk descriptors remain its scheduling cost.
Domain dimensions and canonical Morton keys must match the workspace exactly;
activation/deactivation requires a newly prepared workspace, even when the
active-cell count stays the same. Empty active sets refuse preparation.

Gradients use the actual force-centered velocity, centered differences across
tile boundaries and the dense owner's quadratic stationary-halfway-wall
stencil at inactive tiles and external links. Inactive space retains its solid
wall meaning. Native molecular collision settings are the input on every step;
eddy viscosity is never accumulated from the prior step. The shared viscosity
ratio bound refuses rather than clips. Zero coefficient uses the original
scalar-collision path with unchanged model bits. Both serial and pooled
publication retain populations and step count on closure, native collision or
cancellation refusal, including a request observed after streaming.

`d3q19::sparse::sgs::tests` checks dense/sparse local-collision and population
parity for all three closures with forcing and internal/external walls,
analytic quadratic wall and halo shear, zero controls, sparse memory/topology
binding, mass conservation, real worker-count identity and cancellation/retry.
These are local/integration G0/G1/G4/G5 checks. Diagnostics remain instantaneous
shear dissipation in lattice units; no heat, free-surface closure, wall model,
DNS calibration, developed-turbulence validation or acceptance render is added.

#### Source-bound molecular viscosity under sparse SGS

`MaterialScaling::bind_sparse_sgs_materials` samples one original
`FluidTransport` response per active sparse cell, retaining that response, its
admitted transport snapshot and its applied molecular relaxation in
`sgs::material::SparseMaterialSgs3`. Canonical Morton topology and the native
collision law bind the field; inactive cells are neither sampled nor allocated.
The existing momentum-scaling routine supplies the physical conversion.
Reference density must remain unchanged; thermal coefficients are retained
without imposing a thermal relaxation bound on this isothermal path.

`SparseMaterialSgs3::{prepare_sgs, step_serial, step_pooled}` add the shared
closure's eddy viscosity to each cell's own molecular viscosity and apply the
ratio refusal against that same local value. Zero SGS retains the heterogeneous
molecular field. Repeated steps allocate no new material field, resample no
source and never compound a previously augmented viscosity. Original collision
families and non-shear controls remain intact. Topology, source, numerical or
cancellation refusals publish no flow update or partial binding.

Focused tests cover original-source identity, physical scaling, heterogeneous
zero-SGS native equivalence, analytic local eddy addition, local ratio refusal,
pooled identity and stale-epoch/cancellation refusal. This is a frozen source
epoch: temperature/composition evolution, source refresh, nonlinear shear-law
closure, variable-density remapping and developed-turbulence validation remain
outside this adapter.

### 3-D free surface on sparse tiles (`d3q19::freesurface3`, bead sxnm)

- `FreeSurface3::new(nx, ny, nz, tau, gravity, sigma, contact, fluid)` /
  `with_solid_tiles(..., solid_tiles, fluid)` — the 2-D Körner VOF
  machinery ported onto the sparse grid: per-lane `Cell3`
  (fluid/interface/gas) and tracked mass keyed by MORTON KEY
  (activation-reshuffle-proof), interface inserted between fluid and gas,
  closure asserted. Solid tiles are permanent walls that advance-driven
  activation refuses to open (fluid declared inside one panics).
- `step()` — refresh smoothed fill, collide wet lanes via the shared
  `collide_cell3`, pairwise-antisymmetric mass exchange, pull-stream with
  gas reconstruction at the `σκ/cs²` reference density, conversion
  cascade with the conservation carry. ACTIVATION INVARIANT: every newly
  created interface cell immediately activates every inactive in-domain
  tile in its 18-neighborhood (activate-on-fluid-arrival deadlocks: an
  interface backed by the activation frontier pins as fluid-against-wall).
- Probes: `fill`, `ledger_mass` (fluid Σf + interface m + carry),
  `region_mass(pred)` + `carry()` (ledger partitioning), `cell`,
  `wet_extent`, `surge_front_x(z_below)`, `conversions()`,
  `worst_step_violation()` (per-step ledger violation, logged not
  assumed), `retire_gas_tiles()` (ledger-neutral memory maintenance).
- `FreeSurfaceError3` typed refusals: `Grid`, `Collision` (Morton key +
  lane + source), `NoFluid`, `InvalidContactAngle`, `Cancelled`.
- Goldens: `fs-lbm:d3q19-pour` frozen four-quadrant (see
  golden-couplings.json). Batteries: G0 ledger algebra, G2 Martin-Moyce
  envelope (reference data in `data/reference/martin-moyce-1952.jsonl`,
  report-only curve comparison), G5 bitwise replay.
- No-claims: no throughput (712t), turbulence, calibrated sessile-drop or
  dynamic contact-line accuracy, thermocapillary traction, phase evolution,
  moving substrates or closed surface-energy balance. The original curvature
  clipping and smoothed-fill discretization remain. Quantitative capillary
  oscillation/spurious-current validation and Martin-Moyce central-band
  agreement require further resolution studies. Internal solids are tile-granular.

## Invariants

- Sponge: the acoustic reflection coefficient is MEASURED, not
  assumed — density pulses at two spectral contents read
  R = 1.1e-4 / 2.3e-4 (authored ceiling 5e-3, ~20x headroom;
  re-measure in the same commit if the fixture is re-dimensioned)
  while the bounce-back-wall control measured identically reads 0.75
  and must exceed 20x the sponge's R (the disabled-sponge mutation
  cannot pass); a target-matched sponge is a NO-OP on the matched
  state to 1e-14; left/right sides symmetric; bitwise deterministic.
- The equilibrium recovers its density + momentum moments exactly.
- Unforced collision preserves density and all three momentum components to
  roundoff. The boundary grid delegates to the public cell authority; the Duct
  tile scalar/SIMD twin is differential-tested against its frozen axial oracle.
- Equal-rate central-moment collision reduces to BGK within deterministic
  transform/solve roundoff, including nonzero body force. Split higher-order
  relaxation changes nonequilibrium modes without relaxing degree-zero/one
  invariants. Both moment families deliver the imposed momentum increment
  and relax physical viscous stress at the declared second-order rate.
- Reduced-cumulant collision preserves degree-zero/one invariants, fixes the
  discrete equilibrium to solve roundoff, and is covariant under Cartesian
  axis permutation when rates are equal within each represented order. Its
  nonlinear fourth-order relaxation does not reduce to BGK at equal rates.
- MASS is conserved by a closed-domain step (collision, forcing, streaming,
  and bounce-back all conserve mass). Prescribed velocity/pressure faces are
  open-system flux boundaries and do not claim global mass conservation.
- Duct SIMD streaming is a pure population-bit permutation. Its x/y wall
  crossings read the opposite population at the destination cell, including
  corners and simultaneous z wraps; every non-wall z source is periodic.
- Steady Poiseuille channel flow matches the analytic parabola to a few percent
  (halfway bounce-back resolves the quadratic profile).
- `plan_scaling` derives `τ = 3ν + ½`, flags `stable` iff `τ > ½` AND
  `Mach < MACH_LIMIT`.
- General dense-grid constructors reject zero dimensions and nonphysical
  relaxation times before arithmetic can produce NaNs.
- Gas cells do not act as fluid population sources in the plain dense-grid
  stream step; absent gas-side populations bounce at the fluid boundary until
  explicit free-surface bookkeeping lands.
- D2Q9 wall momentum exchange sums `-2 f_post c_q` in fixed row-major/link
  order for selected stationary halfway-bounce-back links. Gas and exterior
  bounces are excluded; an isolated wall in equilibrium fluid has zero net
  impulse, and selecting one obstacle cannot silently include another.
- D2Q9 moving-wall streaming adds
  `2 w_q rho_post (c_q dot u_wall) / c_s^2` to each incoming wall link. Its
  selected-wall impulse uses the boundary-relative Wen-style exchange
  `(c_out - u_wall) f_out - (c_in - u_wall) f_in`. The receipt independently
  retains the resolved fluid-population impulse and `u_wall (f_in - f_out)`;
  their linkwise balance closes in fixed row-major/direction order. Torque uses
  destination-local halfway-link midpoints about the caller's finite origin.
  The force convention is from Wen et al., *Galilean Invariant Fluid-Solid
  Interfacial Dynamics in Lattice Boltzmann Simulations* (2014,
  <https://arxiv.org/abs/1303.0625>).
- D2Q9 curved moving-wall streaming uses the linear BFL rule with fraction
  `r = distance(fluid center, wall intersection) / link length`. Writing
  `delta = 2 w_q rho_post (c_q dot u_wall) / c_s^2`, the short-link branch is
  `f_in = 2r f_out + (1-2r) f_far + delta` for `r < 1/2`; the long-link branch
  is `f_in = (f_out + delta)/(2r) + (2r-1) f_opposite/(2r)` for `r > 1/2`.
  At `r = 1/2` the established moving-halfway operation tree is used exactly.
  The complete geometry map is deterministic and exhaustive: every actual
  fluid-wall pull link has one entry and every other slot is empty. Short links
  require the next lattice node away from the wall to be fluid or interface.
  Receipt torque uses the declared intersection `x_fluid - r c_q`, while its
  impulse/mass/work balance retains the existing boundary-relative convention.
  The interpolation rule follows Bouzidi, Firdaouss, and Lallemand (2001,
  <https://doi.org/10.1063/1.1399290>).
- D2Q9 partial saturation uses the Noble-Torczynski non-equilibrium-bounce-back
  variant with lattice time step one. For solid fraction `epsilon_s`,
  `B = epsilon_s (tau - 1/2) / ((1 - epsilon_s) + (tau - 1/2))`, and
  the unforced law is `f_post = f + (1-B) Omega_f + B Omega_s`, where
  `Omega_s[q] = f[opp(q)] - f[q] + feq[q](rho,u_s) - feq[opp(q)](rho,u)`.
  With imposed Guo force `F = rho*g + fext`, the same operators act on physical
  populations `h = f + S(u,F)/2`; `B*(Omega_s(h)-Omega_f(h))` is added to the
  native forced BGK result. The complete fluid force is applied once, without
  an occupancy multiplier, and its account stays separate from solid reaction.
  The solid-coupling correction conserves population mass to roundoff; its
  fluid momentum and angular-momentum changes are equal and opposite to the
  solid receipt. Cells and directions reduce in fixed row-major order. The
  construction follows Noble and Torczynski (1998,
  <https://doi.org/10.1142/S0129183198001084>).
- D2Q9 under-resolved cylinder lubrication is an **Estimated** local subgrid
  closure. With unit normal `n` directed from solid 1 to solid 2, signed gap
  rate `s = (u_2 - u_1) dot n`, reduced radius
  `R* = (1/R_1 + 1/R_2)^-1`, supplied gap `h`, positive floor `h_min`, and
  cutoff `h_c`, it uses `h_eff = max(h, h_min)` and
  `lambda = 3 pi sqrt(2) mu (R*)^(3/2)` times
  `(h_eff^(-3/2) - h_c^(-3/2))` when `h < h_c`; otherwise `lambda = 0`.
  The one-step correction is `J_1 = lambda s n`, `J_2 = -J_1`, with pair
  work `-lambda s^2 <= 0`. Kromkamp et al. derive the leading equal-cylinder
  force per unit length (2005, <https://doi.org/10.1017/S0022112005003551>);
  the unequal-radius expression is the local parabolic-gap Reynolds inference,
  not a separately validated unequal-cylinder result. Subtracting the
  resistance evaluated at `h_c` follows Nguyen and Ladd's general cutoff rule
  (2002, equations 34-35,
  <https://doi.org/10.1103/PhysRevE.66.046708>); their three-dimensional sphere
  exponent and calibrated cutoff values are not transferred. Powers of
  three-halves use `x * sqrt(x)` and pair evaluation has no reduction-order
  dependence.
- D2Q9 wall-topology replacement is a two-phase publication. Fresh cells use
  an equal-weight average of unique, surviving pre-transition one-ring fluid
  populations, relaxation times, and external forces; newly covered
  populations are cleared. The receipt's fresh-minus-removed mass and momentum
  equal the active-grid change to roundoff. Donor discovery and all receipt
  reductions use fixed row-major/direction order.
- D2Q9 regularized x faces impose the declared inlet velocity and outlet
  density to roundoff while copying the complementary moment and independently
  measured physical Guo non-equilibrium stress under each cell's `rho*g + fext`
  from the first interior column. Measured and unmeasured paths are bit-identical.
- Rheology laws reject inadmissible physical parameters. An unclipped local
  update closes its stress/viscosity relation to a relative tau defect of
  `1e-12`; the previous tau only initializes the solve. Every update reports
  floor/cap counts when the requested law leaves the representable tau window,
  and `max_shear` is consistent with the applied relaxation time.
- Thermal wall populations encode the declared wall temperatures, so the
  public `temperature` query is consistent on wall and fluid rows.
- Free-surface steps conserve the tracked ledger mass (fluid `Σf` plus
  interface mass plus carry) to the test tolerance, and gas/interface/fluid
  conversions are counted rather than hidden.
- Dynamic wetting admission requires both endpoint claims to answer under the
  same explicit `SelectionPolicy`. Each endpoint must be dimensionless,
  finite, and strictly inside `(0, pi)` radians, with receding no greater than
  advancing. Both caller inputs and any selected curve whose abscissa is
  `contact_line_speed` must carry SI velocity dimensions. The cloned query
  point binds that axis to the absolute coherent-SI speed; the sign only
  chooses the endpoint. Exact replay reproduces both usage receipts, and the
  retained system hash moves when ordered surfaces, history, claims, or models
  move.
- D3Q19 wall and open masks are disjoint and compiled in tile/lane/direction
  order. Where an open face meets a wall, directions crossing both faces are
  stationary wall-owned while pure open-face directions are reconstructed;
  where moving exterior walls disagree at an edge/corner, the whole cell is
  stationary so diagonal corrections cannot create an unpaired mass source.
- Voxelized topology is initialization-only and two-phase: failed/non-finite,
  all-solid, or open-neighbor-invalid proposals leave occupancy and masks
  unchanged. A committed topology cannot be silently changed after
  perturbation or stepping.
- Stationary planar and voxel bounce-back conserve mass to roundoff. Tangential
  moving-wall corrections sum to zero over the represented link set.
- D3Q19 measured wall loads use the same Wen-style boundary-relative exchange
  as D2Q9, retaining the moving-mass term so `J_wall + delta_p_fluid =
  sum(u_wall * delta_mass)` and its three-axis angular counterpart close.
  Body-force impulse and open-boundary momentum flux remain separate.
  Global angular-momentum closure additionally requires accounting for
  periodic coordinate-wrap flux; the closed-box regression uses no periodic
  faces. This surface does not move voxel topology or integrate rigid bodies.
- On face-interior cells, regularized velocity/pressure reconstruction
  preserves the declared zeroth and first moments and copies the independently
  measured force-corrected non-equilibrium second moment from the first interior cell. Mixed
  wall/open rim cells preserve their per-link wall ownership instead of
  claiming exact target moments.

## Error model

Most operations are total over physically admissible inputs. Constructors and
parameter helpers panic on nonsensical requests: zero dimensions, non-finite
forces/relaxation times, non-positive viscosities/rheology indices, non-positive
Rayleigh height, or non-positive Reynolds/length in the scaling assistant.
D2Q9 momentum measurement additionally requires a full-grid boolean mask whose
selected entries are all `Cell::Wall`; this is checked before a measured step
can mutate populations.
D2Q9 moving-wall calls additionally require one finite velocity per grid cell,
zero velocity on every non-wall cell, speed squared below 0.03, and a finite
moment origin. Every moving wall adjacent to fluid requires positive finite
post-collision density, and every outgoing wall-link population must be finite.
Request fields and post-collision state are admitted before streaming mutates
the grid.
D2Q9 curved moving-wall calls require a full-grid link map that is exactly
populated on fluid/interface-to-wall pulls, fractions in `(0, 1]`, finite
link-intersection velocities below the same low-Mach envelope, and a finite
moment origin. Short-link interpolation requires an in-domain fluid/interface
far donor with a finite outgoing population; long-link interpolation requires
the finite local opposite population. Moving links require positive finite
local post-collision density. Geometry is admitted before collision, and all
required post-collision values are admitted before grid publication.
D2Q9 partial-saturation descriptors require a finite solid fraction in
`[0,1]`, exactly zero velocity at zero fraction, and finite solid velocity
inside the same low-Mach envelope. The full-grid field may assign positive
fraction only to `Cell::Fluid`; all populations must be finite, every active
cell needs positive finite density, finite momentum, and `tau > 0.5`. BGK is
required; imposed acceleration, active-cell external force, the combined force
and physical velocity must be finite, as must the moment origin. Collision is
built in a private proposal and the caller's output/scratch is published only after
every proposed population and the aggregate receipt are finite.
D2Q9 under-resolved cylinder lubrication descriptors require a finite
nonnegative gap; positive finite cutoff, floor, reduced radius, and dynamic
viscosity; `h_min < h_c`; a finite unit normal within `128 * f64::EPSILON`;
finite surface velocities individually inside the squared-speed `< 0.03`
low-Mach envelope; and a finite line-of-action point. Evaluation requires a
finite moment origin and refuses non-finite derived resistance, gap rate,
impulse, work, offset, or angular impulse. Negative gaps belong to contact and
are refused; exactly zero gap is admitted only through the explicit positive
floor. The descriptor and evaluation mutate no state, so refusal is atomic.
D2Q9 wall-topology replacement currently requires a fluid/wall-only domain, a
full-grid target mask that leaves at least one fluid cell, positive finite
covered/donor population mass, finite donor populations/momentum/external
force, and donor relaxation times above 0.5. Every fresh cell needs at least one
unique surviving one-ring fluid donor. The complete proposal and receipt are
validated before flags, populations, relaxation times, or forces are
published.
D2Q9 regularized x flow requires at least three columns, non-periodic x,
periodic y, fluid face/first-interior columns, complete finite gravity/external
forcing, representable current-cell and prescribed-outlet total forces,
a positive finite outlet density, and finite inlet speed squared below 0.03.
These setup conditions are checked before collision. Predicted source velocities
must also stay inside that low-Mach envelope and all reconstructed face
populations must be finite with positive density; every row is admitted before
streaming, so late refusal preserves the grid and restores collision scratch extent.
Dynamic wetting queries return typed errors when either `QtyAny` input is not a
velocity; for a non-finite signed speed or negative/non-finite pinning speed;
when a selected `contact_line_speed` curve declares non-velocity abscissa
dimensions; for either material-query refusal, non-dimensionless or
out-of-range angles, or a reversed hysteresis bracket. They query both
endpoints before selecting one, so an incomplete card refuses even when the
current direction would consume the endpoint that is present. The adapter
mutates neither the card nor free-surface state.
D3Q19 boundary construction additionally rejects non-4-multiple dimensions,
tile-count overflow, invalid collision parameters, a nonzero force with an
unsupported collision model, unpaired periodic faces, non-tangential/non-finite or
outside-low-Mach wall velocities, non-finite or outside-low-Mach inlet
velocities, non-positive pressure density, and open faces on more than one axis.
Voxelization rejects non-finite samples, an all-solid domain, an obstructed
first-interior open-face layer, and topology mutation after initialization.
Boundary-grid perturbation rejects non-finite amplitudes or magnitudes at least
one before changing populations or locking topology.

D3Q19 load measurement rejects non-finite origins, links outside the actual
wall topology, inconsistent cell/tile/lane coordinates, duplicates, and
out-of-order selections before collision. Non-finite loads are refused before
population publication. The pooled load path returns `BoundaryStepError3`,
polls cancellation at intervals of at most 256 measured links and before
publication, and returns no loads or new fluid populations on refusal.
`collide_cell3` returns typed errors rather than publishing non-finite or
non-positive cell states. Both moment-space rungs additionally refuse rates
outside `(0, 2)`. Their fixed lattice inverse requires no numerical pivot or
rank-admission threshold; non-finite outputs still refuse. Body
forcing is supported by all three collision families; regularized open
velocity/pressure faces reconstruct physical moments in each cell's Guo frame.

## Determinism class

Fully deterministic: fixed cell iteration order for D2Q9 and fixed logical
tile/cell/direction identities plus face/z/y/x reconstruction order for D3Q19;
no RNG. Duct collision retains the frozen direction-ascending reductions and
scalar operation tree in fixed lane blocks. Duct pull-streaming uses a pinned
direction/z/y/x-tile row schedule, but only moves bits into slot-exclusive
destinations, so schedule changes cannot create arithmetic reassociation.

## Cancellation behavior

Dense-grid steps are synchronous (no polling). `SparseGrid3::step_pooled`
polls `Cx::checkpoint` once per kernel tile: a tripped gate drains the pass,
surfaces `SparseError3::Cancelled`, and leaves the pre-step state intact. The
collide pass writes its retained scratch buffer and the stream pass writes a
third retained transactional destination; neither becomes published state
until both passes complete. A re-issued step is bitwise-identical to an
uncancelled one.

## Unsafe boundary

Five registered leaf capsules opt out of the workspace default at module scope:
the collision pair under `src/d3q19/simd/{neon,x86}/mod.rs`, the pull-stream
pair under `src/d3q19/simd/stream/{neon,x86}/mod.rs`, and the cell-collision
FMA capsule `src/d3q19/fma/mod.rs`. Each stays below the
capsule line cap with an adjacent `SAFETY.md`. Safe façades validate tile shapes,
Rust borrows prove aliasing and lifetimes, and x86 one-shot selectors admit
private thunks only after their operation-specific feature checks. The FMA
capsule's one-shot selector admits its thunk only after
`is_x86_feature_detected!("fma")`, hot loops hold the selected kernel across the
loop, and the thunk runs the safe portable body.
Miri selects safe scalar twins.

## Feature flags

None.

## Conformance tests

`tests/sponge_reflection.rs` (8): measured reflection coefficients at
two pulse widths + wall control; matched-state no-op; left-side
symmetry; four constructor refusal panics; bitwise determinism.

`tests/lbm.rs` covers the v0 core: equilibrium moments; mass conservation;
Poiseuille flow matches the analytic parabola (symmetric, centered); the
scaling assistant derives τ + flags stability + colors the plan; it rejects a
high-Mach plan and nonsense inputs; determinism; exact one-link wall-impulse
sign/magnitude, obstacle selection, equilibrium cancellation, replay
determinism, and pre-step mask refusal for D2Q9 momentum exchange; an
independently enumerable moving-wall link that pins bounce correction,
boundary-relative force, torque, work, and moving-mass balance; exact
zero-velocity compatibility with the stationary API; moving-step bit replay;
and fail-closed moving-field admission. A periodic-x, two-moving-wall Couette
deck checks the manufactured linear relative-velocity profile, exact replay,
separate top/bottom impulse and receipt balance, and a measured common-boost
profile/force residual inside one low-Mach lattice envelope. A separately
enumerated partial-saturation cell pins every population of the published
blend, mass and equal-and-opposite impulse/torque/work receipts, origin
translation, exact zero-fraction compatibility, the full-fraction stationary
fixed point, replay, descriptor admission, and atomic refusal. Linear
curved-wall fixtures independently pin both BFL branches, off-lattice torque
arms, full receipt balance, exact halfway compatibility, replay, checked link
construction, exhaustive geometry, far-donor requirements, finite-population
admission, and atomic refusal. A separate topology-transition fixture
pins unique-donor fresh-cell initialization, exact covered/fresh counts,
mass/momentum delta closure, replay, idempotence, and atomic no-donor/mixed-domain
refusal. The file also covers force-free and forced regularized x-face
moment/stress reconstruction, nonunit outlet density with unequal source and
destination forces, measured-path bit equivalence and unchanged wall impulses,
nonfinite/overflow forcing refusal, and late-source refusal without partial
publication.

`tests/extensions.rs` covers the current extension scaffolding: power-law and
Newtonian-limit channel profiles, Carreau plateaus, Rayleigh-Bénard onset
bracketing and Nusselt heat transport, gas-neighbor streaming behavior, thermal
wall-temperature queries, invalid-parameter rejection, free-surface mass-ledger
conservation, qualitative dam-break front advance, rotation equivariance,
contact-model bracketing, qualitative jet fragmentation, and dynamic wetting
selection, receipt replay, interface-system identity binding, and typed
refusal paths.

The rheology unit battery independently checks strongly thinning and
thickening constitutive roots, previous-tau independence in both dimensions,
repeated-update identity, and clipping at the applied strain rate.
For D2Q9 central-moment collision, the trace uses its independent unit bulk
relaxation time; only deviatoric strain varies with the solved shear tau.
The manufactured tensor battery covers mixed dilation/shear and Guo forcing,
and the local-root battery includes fixed trace contributions.
`strongly_thickening_channel_matches_analytic_profile` adds an n = 5 channel
with 33 fluid rows: it must reach steady state before 80,000 steps and keep
worst pointwise velocity error below 3% of the analytic peak.

`Rheology::HerschelBulkley` extends the same D2Q9/D3Q19 local solver with
`nu = k max(gamma, gamma_min)^(n-1) + Y (1-exp(-m gamma))/gamma`, where
`Y` is yield stress divided by density. `n=1` recovers the
Bingham–Papanastasiou law. The regularization time `m` and power-law minimum
rate `gamma_min` are explicit inputs; the yield contribution at zero strain
is exactly `Y*m`. Deterministic `expm1` and a small-rate series preserve that
limit and the derivative without a second constitutive iteration. This is a
finite-viscosity regularization with sub-yield creep, not an exact rigid plug,
cessation criterion, thixotropic memory or experimentally qualified material.
The existing tau floor/cap remains a separately counted resolution limit;
capped cells cannot establish physical yield behavior. The formula follows
Papanastasiou (1987), DOI:10.1122/1.549926, as written for Herschel–Bulkley by
Jeon & Hodges (2018), equation (6), DOI:10.5194/nhess-18-303-2018; the explicit
power-law plateau additionally regularizes its otherwise singular `n<1`
zero-strain term. The focused unit battery checks plateaus, monotone stress,
zero-yield reduction, invalid inputs and constitutive closure.
`regularized_bingham_channel_matches_stress_integral_in_both_dimensions`
specifies native 2-D/3-D channel profiles against an independent integrated
momentum balance, a nearly flat creeping core, mass balance and absence of
numerical viscosity clipping. Full-crate execution of that added channel test
is not established by a standalone constitutive-kernel check.

`rheology::RheologyField<2>` and `<3>` retain one prepared constitutive law
per native D2Q9 or D3Q19 plate-grid cell, with x-fastest full-grid extents
including walls. Newtonian, power-law, Carreau and regularized Herschel–Bulkley
laws can coexist. Construction checks dimensions, checked cell multiplication,
exact law count and all parameters before allocation. Complete `refresh`
admits every proposal before changing any accepted law and reuses storage;
`set_law` admits one local replacement. Callers may refresh parameters from
temperature/composition between steps; the field does not transport material
labels, infer mixture/interface laws, evolve temperature or add stress memory.

`update_tau_field` and `d3q19::update_tau3_field` share the uniform updater's
implicit local constitutive solve, Guo correction and clipping statistics.
D2Q9 central moments retain separate bulk/shear strain reconstruction. Static
dispatch preserves the scalar path, and prepared-field updates allocate nothing.
Shape/storage mismatch refuses before relaxation writes; physical populations
retain the original updater's validity preconditions. Exact production-source
native tests cover mixed roots, clipping, storage reuse and atomic admission;
full-grid regressions specify force/trace invariance, inactive cells and
scalar/field bit parity in both dimensions, with execution still pending a
resource-capable native build.

`rheology::physical::PhysicalRheology` binds caller-declared SI Newtonian,
power-law, Carreau and regularized Herschel–Bulkley parameters through the
existing `MaterialScaling::bind_rheology`. `ScaledRheology` retains the original
material scale, source, frozen reference density, physical law and declared
`LiquidViscosityBasis`; dynamic coefficients divide by density exactly once.
For spacing dx and step dt, kinematic viscosity scales by dt/dx², power-law
consistency by dt^(2-n)/dx², yield stress/density by dt²/dx², constitutive times
by 1/dt, and minimum shear rate by dt. Physical power laws require an explicit
SI low-rate cutoff and use the existing zero-yield Herschel–Bulkley image so
that grid refinement does not change that material plateau.

Normal conversions retain the original material-scaling operation order;
deterministic logarithmic fallback admits representable results despite
intermediate overflow/underflow. Invalid parameters or a lost nonzero term
refuse. Uniform and spatial 2-D/3-D updates consume the same converted law;
existing tau floor/cap accounting still bounds numerical resolution. No source
fit, new material evolution or stability guarantee is inferred. Four focused
regressions specify physical stress/local-root equivalence under refinement,
dynamic/kinematic agreement, extreme-scale refusal and source retention through
both field dimensions. Syntax validation is established; native execution of
these new whole-crate tests remains pending sufficient host resources.

`tests/extensions.rs` (tfz.19): lbm-101 power-law/Carreau profiles vs
analytic (0.14% for the n = 0.8 profile with τ floor+cap ledger);
lbm-102 Rayleigh–Bénard onset
bracket (decay Ra=1200 / growth Ra=2500, Nu>1); lbm-104 STRICT free-
surface mass ledger (5e-14 over 600 dam-break steps, conversions
counted); lbm-105 dam-break front envelope (coarse honesty band);
lbm-106 G3 rotation equivariance (3e-14); lbm-107 contact-model
bracket band + Plateau–Rayleigh jet fragmentation with strict ledger;
lbm-108 level-jump refinement (Poiseuille through the interface +
shear-wave decay-rate transmission, first-order labels); lbm-110 through
lbm-112 dynamic wetting interpolation/pinning, G5 receipt replay and ordered
system identity, and typed incomplete/ambiguous/out-of-domain/unit/bracket
refusals.

`tests/d3q19_battery.rs` covers the frozen D3Q19 core: exact integer lattice
moments/opposites, equilibrium moments, mass conservation, analytic
rectangular-duct flow, replay determinism, the registered core golden, shared
kernel bit-equivalence, BGK/central-moment equal-rate equivalence, split-rate
conservation, reduced-cumulant equilibrium/nonlinearity/axis-covariance laws,
and fail-closed inputs. The D3Q19 `forcing_tests` unit module adds independent
uniform-acceleration populations, forced equal-rate central/BGK equivalence,
and prescribed physical-stress relaxation for both moment families.

The D3Q19 SIMD module adds G0 seeded differential batteries for both operations.
Collision compares every lane to its scalar tile twin and directly locks that
twin to the frozen axial cell authority. Pull-streaming compares every
population bit on a single-tile fixture and an asymmetric 2x3x4-tile fixture;
independent route anchors cover inter-tile x/y/z moves, periodic z, corner
bounce, and wall precedence over simultaneous z wrap. Both batteries retain
first-divergence JSONL receipts, hostile dispatch-selector cases, and one-shot
table identity. The existing Duct golden remains the end-to-end bit-surface
gate.

`tests/d3q19_sparse.rs` drives both entry and mid-stream cancellation. Its
deterministic one-worker runner lets collision finish, requests cancellation
after the first stream kernel tile, and pins unchanged published bytes, zero
completed steps, drained structured refusal, and bitwise serial-equivalent
fresh-gate reissue.

`tests/d3q19_boundaries.rs` (bead 40p2) covers all six hand-enumerated planar
link masks, aligned deterministic mask ordering, the exact 18 links around one
voxel, atomic immutable SDF topology and rejection paths, axis-generic
regularized face density/velocity plus independent inlet/outlet stress
reconstruction, grid-level collision selection and equal-rate BGK/central
agreement, split-rate cumulant conservation, uniform plus per-cell force
impulse, and force-driven planar Poiseuille flow for both moment families.
The forced channel gate requires less than 0.4% peak-normalized profile error
at 16 fluid rows and better than 2.5-fold error reduction from 8 to 16 rows.
The battery also covers open/wall and moving/open rim ownership, stationary
planar and voxel mass conservation, exact one-step
moving-lid momentum oracles, qualitative primary cavity circulation,
pressure-driven Poiseuille shape, and a boundary replay-hash candidate. Ignored
release fixtures carry the full 10,000-step leak and 32x32 full-rim 3%
pressure-Poiseuille gates.

`tests/cylinder_re100.rs` (frankensim-wghy) keeps the inexpensive evidence
machinery in the default suite: a checked linear-detrended, Hann-window
`fs-fft` lift-frequency estimator with peak-prominence, exact-bin,
affine-trend, replay, flat/non-finite/ambiguous/Nyquist, and
resolution-refusal oracles; and an exactly projected stair-step cylinder mask
pinned by occupied cells, centroid, reflection symmetry, fluid-wall link count,
and FNV-1a fingerprint. The ignored release-scale `lbm-109` fixture runs Re=100
at diameter 10, nominal inlet speed 0.1, 32D streamwise extent, 12D and 16D
periodic lateral spans, 8,192 warm-up steps, and 32,768 retained force samples.
It normalizes drag by measured mean inlet density, requires raw-domain (`[1.20, 1.55]`)
and split-window (`<= 2%` raw, `<= 3.5%` under `(4, -3)` two-width extrapolation)
Cd guards before admitting a bounded (`<= 0.10`) empirical two-width
zero-blockage Cd intercept, and requires the full-window and both half-window
intercepts to stay in `[1.25, 1.45]`. It does
not extrapolate spectral bins: both widths and both half-windows must put St in
`[0.155, 0.175]`, split windows must agree within `0.0062`, and the reported 16D
St must agree with the 12D sensitivity run within `0.01`. The deterministic
symmetry-breaking seed is disclosed as offset `(10, 6)` at transverse speed
`1e-4`. Roshko NACA TR-1191, Posdziech-Grundmann 2007, Behr et al. 1995, and
Maskell ARC R&M 3400 are cited at the executable gate; no publisher artifact is
redistributed.

### Collision-operator battery (`tests/regularized_collision.rs`)

Hermite-state BGK equivalence of the projected regularization (and
RR3 at zero velocity); TRT-at-BGK-magic equivalence; exact
mass/momentum conservation of every operator; force-independent impulse
and stress under homogeneous forcing; the shear-wave and body-force-driven
Poiseuille viscosity identities for every operator; the executed
unresolved-shear-layer stability asymmetry (BGK/RR3/CM survive,
projected Regularized destabilizes — also the mode-discrimination
control); bitwise determinism of every non-BGK path.

## D3Q27 and shared reference transport (IRP M3.2-W1)

`d3q27` represents all independent central moments/cumulants through (2,2,2).
Geier et al. (2015), doi:10.1016/j.camwa.2015.05.001 section4 is the collision
reference. A separable central transform and normalized multivariate
generating-function recurrence retain all fourth/fifth/sixth couplings. This
is not D3Q19's ReducedCumulant projection. All ten publication-order relaxation
rates are explicit, including separate bulk/shear groups and the diagonal
cubic-alias corrections (58--63). The central-moment operator uses Gaussian
central targets; cumulant equilibria above degree2 are zero. Cumulant midpoint
forcing changes the first cumulants as in85--87. Central-moment forcing uses
the Gaussian source moments and grouped trapezoidal factors of Fei, Luo & Li
(2018), doi:10.1103/PhysRevE.97.053309, equations5,20--21: the six mixed third
moments receive `(1-omega3/2) cs^2 F/rho`, and the three fifth moments receive
`(1-omega9/2) cs^4 F/rho`. The third-order pair differences have zero source;
their separate omega4 does not set the source factor. Bit semantics v2 corrects
these previously omitted force terms. An explicit zero-force branch preserves
all unforced arithmetic; BGK and cumulant collision arithmetic are unchanged.
Finite-force flow accuracy remains a downstream qualification, separate from
conservation and these local constitutive moment balances.

Bit semantics v3 additionally matches factorized BGK forcing to
`S = (F/rho) . d_u feq`. Cached marginal derivatives retain the full represented
Maxwellian source tensor, including mixed third/fifth moments omitted by the
quadratic source in a moving fluid. Quadratic BGK, every unforced BGK path,
central-moment and cumulant arithmetic remain unchanged. The G0/G3
`factorized_bgk_force_matches_maxwellian_moments_in_moving_frames` regression
checks all 27 moments and independent product derivatives across densities,
moving frames, force signs/axes and relaxation times 0.5005, 0.8 and 1.7.
This establishes local kinetic consistency, not full Galilean invariance or
finite-force continuum convergence. The unforced cumulant golden is unchanged.

CollisionConfig27 requires model/rates or BGK tau and a Mach cap in (0,0.3].
Tau must exceed0.5 and each rate must be strictly in (0,2). Nonfinite state,
force or parameters, nonpositive density and excessive Mach refuse. No
population clipping, adaptive rates or unspecified defaults are introduced.
BGK explicitly chooses a quadratic or factorized equilibrium; like-polynomial
D3Q19 BGK is the lattice control.

`lattice3::Grid3<L,Q>` shares checked dense reference transport between D3Q27
and an adapter calling the existing D3Q19 cell owner. Collision, force,
periodic axes, immutable solid mask, initial populations and cell/Mach caps
are required inputs. Nonperiodic exterior and solid links use stationary
halfway bounce-back. Fixed64-cell tiles share arithmetic between serial Cx
and pooled KernelRunner execution. Three buffers preserve published state and
step count on collision, pool or cancellation refusal. Cancellation is polled
at each tile and before publication. The library spawns no threads itself;
specialized D3Q19 sparse/boundary/SIMD owners are not redirected or duplicated.

`Grid3::wall_momentum_probe` compiles selected solid or nonperiodic exterior
halfway links once. `step_serial_with_wall_momentum` and its pooled counterpart
return the wall impulse and torque from the actual outgoing post-collision and
incoming streamed populations, including pressure and viscous traction. Each
corner/edge link is counted once. The optional collision field retains the same
material or SGS viscosity binding as the ordinary per-cell step. Only selected
links are revisited for measurement; compensated sums use canonical link order
for serial/worker-count agreement and resolve drag against cancelling pressure.
The receipt and populations publish together: topology, arithmetic, collision
or cancellation refusal leaves accepted state untouched. Positions are halfway
intersections with cell centres at integer coordinates plus one half; torque
uses that domain image about the explicit origin. Force scales by
`rho_scale*dx^4/dt^2`, torque by `rho_scale*dx^5/dt^2`. This stationary-wall
measurement does not add moving-wall work or a continuous-image torque model
for bodies crossing periodic seams. The `wall_momentum_*` native battery checks
channel pressure, force/torque balance with applied forcing and an obstacle,
material-field stepping, unchanged populations and worker counts 1/2/7.

AcousticScalingPlan::from_gas uses the supplied GasState nu=mu/rho and c plus
explicit dx, maximum physical speed and Mach cap. It returns dt=sqrt(1/3)dx/c,
nu_lu=nu*dt/dx^2, tau=0.5+3nu_lu and the existing ScalingPlan. Invalid,
excessive-Mach and unrepresentable scalings refuse.

`tests/d3q27_battery.rs` specifies exact rational G0 lattice laws, equilibrium
and10000-state conservation at1e-14, independent set-partition cumulant
oracles through degree6, refusal boundaries, G4 cancellation and G5 actual
workers1/2/7. G0/G3 force regressions independently differentiate a product
equilibrium and check all27 central moments across force signs/axes, moving
frames, densities and unequal grouped rates. A separately prescribed Hermite
stress checks that forced shear and bulk relaxation retain their independent
rates, alongside the third/fifth source balances. The product-state collision
hash is unchanged by v2 because its fixture is unforced cumulant collision.
An independent dyadic oracle and focused exact-source x86_64 debug/release
replay reproduce `b256c77f0c5c1bf7`; no cross-ISA qualification is added.
G1 requires >=1.8 shear-decay error
order on200/100/50micrometre acoustic rungs and live, shrinking higher-rate
spread. G2 uses four Gaussian axial-jet cross-section rungs (n = 16 to 128;
the central-moment operator converges slowest at tau - 1/2 = 0.003) and Fourier
angular interpolation: all D3Q27 collisions decrease m4/m0 below0.5%, while
same-polynomial D3Q19 BGK must be at least2x worse, or the fixture is refused.
The rotational fixture uses diffusive scaling, separate from the fixed-c
acoustic ladder. `scripts/e2e/music/lbm_qualification.sh` retains both release
lanes. G1 and G2 pass natively in release; G2's 128 rung takes about ten
minutes there and does not finish within 90 minutes in a debug build.

This is not a White--Chong nozzle reproduction, full jet/edge-tone campaign,
total-energy/NSF solver, 3-D multiblock refinement, transparent boundary,
high-Re stability/performance, cross-ISA identity or instrument acceptance.
M3.2-W2a/W3a and M1.1-W8 must consume the explicit API and qualify their own
physics. Historical D3Q19 arithmetic and semantic/golden versions are unchanged.

## Exploratory three-dimensional SGS (`sgs3`, IRP M3.2-W5 partial)

`sgs3::eddy_viscosity` is an allocation-free, unit-consistent Smagorinsky/WALE
closure for a supplied full velocity gradient. Smagorinsky uses
`nu_t = (Cs Delta)^2 sqrt(2 S:S)`. WALE uses the traceless symmetric part of
the **matrix square** of the velocity gradient, including rotation, following
Nicoud & Ducros (1999), DOI `10.1023/A:1009995426001`. Tensor normalization
avoids high-power overflow/underflow; no arbitrary denominator epsilon or
viscosity floor changes its near-wall or refinement asymptotes. Zero coefficient
is an explicit disabled control; negative/nonfinite coefficients, nonpositive
filter width and nonfinite gradients refuse.

`Grid3::prepare_sgs`, `step_serial_with_sgs` and `step_pooled_with_sgs` are
explicit opt-in exploration APIs. `SgsConfig` requires the model, filter width
in lattice cells, maximum eddy/molecular viscosity ratio and intended use.
`SgsUse::AcceptanceRender` refuses with a typed error. Exceeding the supplied
ratio cap refuses rather than silently clipping the closure. No constructor,
ordinary step, instrument render or acceptance path enables SGS by default.
This does not complete the instrument program's render-admission integration
or grant a turbulent-DNS/acceptance claim.

The integration caches one macroscopic velocity/density field, using the
grid's existing uniform force and half-force velocity shift. Interior and
periodic derivatives are centered; stationary halfway wall derivatives use
the quadratic through the actual zero wall value and two fluid centres.
Solid/exterior walls use the same rule; a one-cell gap uses the symmetric
zero derivative. Grid topology supplies these boundary semantics explicitly.
WALE cannot be reconstructed from symmetric nonequilibrium stress alone.
Moving/open walls and independently supplied force fields are outside these
`Grid3` entry points. The force-history owner below supplies its separate SGS
methods without rerouting the original owner.

An optional molecular collision field composes with existing per-cell material
viscosity. The original equilibrium and independent bulk/higher rates survive:
D3Q27 modifies only omega1; D3Q19 retains its pre-existing tied second-order
bulk/shear rate. Each step starts from molecular viscosity, never the previous
SGS-augmented field. Zero SGS retains the original collision bits. A reusable
workspace removes SGS sweep allocations; preparation is an O(N) serial pass
with cancellation every 64 cells, followed by the existing serial or pooled
collision/streaming owner. Refusal/cancellation never publishes populations or
increments the step counter. Per-step memory and computation scale linearly;
no measured throughput or speedup is claimed.

`SgsReport` returns fixed-order, instantaneous pre-step lattice-unit sums
`sum(2 rho nu S_dev:S_dev)` and `sum(2 rho nu_t S_dev:S_dev)`, their ratio
and maxima. These diagnose resolved and modeled **shear** dissipation; bulk
dissipation is excluded (including D3Q19's tied-rate contribution). They are
not accumulated thermal energy and do not close a total-energy or scalar-flux
equation.

`examples/turbulent_taylor_green.rs` runs the same D3Q27 cumulant step owner
with `none`, `smagorinsky`, `wale`, `vreman`, `dynamic-local`, or
`dynamic-domain`. The dynamic variants fit the squared Smagorinsky coefficient
from the current density and velocity using the 27-cell Favre test box; they
reject an authored `--coefficient`. An explicit coefficient admission bound,
eddy/molecular ratio bound and workspace byte budget constrain that fit.
CSV output includes its maximum applied and minimum raw coefficient, negative
fits excluded by the no-backscatter policy, and undetermined cells. Dynamic
fit columns are NaN for fixed models. Explicit model, Reynolds-number,
coefficient and resolution controls accompany CSV energy and instantaneous
dissipation histories. Its `--help` gives a fixed-Re, matched-time refinement
sequence; the example does not substitute for a developed-turbulence benchmark.

`sgs3::model` unit tests check analytic simple-shear behavior, WALE's cubic
near-wall asymptote, rotational/unit scaling and invalid-input/acceptance
refusals. `tests/sgs3_rungs.rs` checks quadratic halfway-wall profiles,
zero-coefficient bitwise identity, preserved independent moment rates,
ratio-cap/acceptance/cancellation refusals and actual 1/2/4-worker identity.
The 8/16/32-cell Taylor-Green ladder at fixed Re1600 verifies the resolved
analytic field's SGS/molecular dissipation ratio decreases with order 2 +/-0.3
and exercises conservative short-time SGS evolution. It does **not** qualify
developed turbulent TGV dissipation histories, a channel DNS profile, jet
physics or finest-rung acoustic equivalence. Those full M3.2-W5/G2 campaigns
remain open. Neither SGS nor LVEL is a resolved-turbulence acceptance render.

### Cartesian geometry for the dynamic LES field

`PeriodicDynamicSgs::new_cartesian` admits uniform Cartesian cells with
independent positive finite spacings `[hx, hy, hz]`. The momentum and scalar
fits use centered derivatives in each physical direction and the equivalent
filter length `Delta = (hx*hy*hz)^(1/3)`, evaluated without forming the volume
product. The 27-cell test box has physical widths `[3hx, 3hy, 3hz]`; its
equivalent squared filter-length ratio remains nine. Unrepresentable geometry
metrics refuse. The isotropic constructor and equal-spacing Cartesian calls
retain the original arithmetic, numerical payload and allocation-free sweeps.

This remains a scalar eddy-viscosity/gradient-diffusion model on a uniform,
fully periodic grid with at least three cells per axis. Correct directional
metrics do not qualify arbitrary cell aspect ratios, tensor-valued SGS
transport, nonuniform grids, walls or turbulent DNS agreement. The geometric
mean scale follows Scotti, Meneveau & Lilly (1993), DOI `10.1063/1.858537`,
whose anisotropic-grid analysis supports it for moderate resolution anisotropy.

### Dynamic scalar transport on the periodic LES field

`sgs3::dynamic::scalar` supplies an independently fitted scalar diffusivity.
`PeriodicDynamicSgs::prepare_with_scalar` prepares momentum and scalar fits
from the same density, velocity and scalar fields, reusing the Favre-filtered
velocity and strain. The scalar Leonard vector and modeled gradient vector
are contracted by least squares, following Squires (1991), equations
(2.13)--(2.19), with the implementation's explicit Cartesian test box. The scalar
model has no stress-tensor factor of two. It returns diffusivity directly;
it never divides by a possibly zero momentum coefficient to form Prandtl.

Domain or local-box averaging operates on numerator and denominator
separately. Scalar offsets and units are removed before the fit; negative
fits follow the explicit no-backscatter policy, excessive positive fits
refuse, and undetermined fits return zero. The consuming solver enforces the
authored diffusivity/molecular bound against its original material source.
The reusable scalar workspace adds 64 nominal numerical bytes per cell,
with actual capacity admission and bounded cancellation. Constant scalars
and uniform velocity use zero diffusivity without consuming stale scratch.
This exploratory gradient-diffusion model supplies no wall or shock model,
SGS energy equation, or DNS qualification.

`scalar::components::PeriodicDynamicComponents` fits one common diffusivity
to a bank of same-unit scalars, including every species mass fraction. It
minimizes the sum of the scalar Germano residuals: numerator and denominator
are summed across components before spatial averaging and division. One
shared normalization preserves their physical relative weight; independently
normalizing each species would give trace species artificial influence.
Constant and absent components contribute zero. Species permutations and
orthogonal same-unit component transformations preserve the exact-arithmetic
fit, with floating reductions following the declared component order.

`prepare_with_scalar_and_components` fits momentum, heat and common species
diffusivity from one preparation of the current flow. Its immutable snapshot
accessor reads the existing composition layout without a transposed copy.
All three transports are independent; positive species diffusivity remains
possible when momentum viscosity and heat diffusivity are zero. The component
workspace adds 80 nominal numerical bytes per cell, independent of species
count, with no allocation during preparation and bounded cancellation.
The native consumer must enforce physical diffusivity limits and preserve
discrete zero-sum mass flux and partial-enthalpy transport. This constrained
vector extension of the scalar fit does not replace molecular differential
diffusion, model unresolved reaction rates, or establish reacting-flow/DNS
accuracy. The existing periodic, grid-aspect and exploratory-use limits apply.

### Dynamic momentum and scalar transport in planar channels

`sgs3::dynamic::channel::{ChannelDynamicSgs, ChannelDynamicScalar}` extends
the shared normalized Favre--Germano kernels to a uniform Cartesian channel
with two periodic directions and one axis bounded by planar halfway walls.
There must be at least three fluid cells in each periodic direction and two
in the wall-normal direction. The equal-volume test box contains nine cells
within the same wall-normal plane; it never samples through a wall or wraps
to the opposite plate. Numerator and denominator are averaged over each
homogeneous plane before division, so the fitted coefficients remain functions
of wall distance and time. Channel configurations make this plane averaging
intrinsic rather than accepting a misleading periodic averaging selection.

The resolved filter length remains `Delta = (hx*hy*hz)^(1/3)`. Expanding only
the two homogeneous widths by three gives the squared test/grid length ratio
`3^(4/3)`, rather than the fully periodic box's nine. This volume-width
convention and homogeneous-direction filtering follow the alternatives
described in Cabot (1991), *Large eddy simulations of passive and buoyant
scalars with dynamic subgrid-scale models*, CTR Annual Research Briefs,
section 2.1. It does not copy that paper's spectral discretization or claim
its channel-flow accuracy.

Declared velocity and scalar wall derivatives use the quadratic through the
two nearest fluid centres and the halfway boundary condition. Each preparation supplies both
uniform wall velocities and independent scalar `Value`, `ZeroGradient`, or
`Interior` policies. A scalar zero-gradient wall requires three normal fluid cells:
the filtered scalar uses a second-order one-sided derivative of its actual
Favre-filtered interior values there. Density weighting need not preserve the
physical scalar's zero normal gradient, so that boundary condition is not
imposed on the filtered quotient. Uniform Dirichlet values do survive the
tangential Favre filter. Prescribed wall values participate in normalization:
a uniform fluid field is not mistaken for zero strain or zero scalar gradient when its
walls differ. Filtering, strain and scalar gradients share the same geometry;
scalar diffusivity remains independently fitted even when momentum transport
is zero. The existing explicit no-backscatter policy, coefficient admission
bounds, finite-value checks and exploratory-use restriction remain in force.

`Interior` reconstructs the wall-adjacent cell derivative from three current
fluid centres at both physical and Favre-filtered levels. It is exact for a
quadratic and second-order for smooth resolved fields; it supplies neither a
wall temperature nor a transport boundary law. Flux/film and coupled-flux
owners can use it without deriving a gradient from a diffusivity still being
fitted, substituting a reservoir for the surface, or consuming an old wall
report. The opposite wall retains its own independent policy. Three normal
fluid cells are admitted before uniform-field shortcuts whenever interior
reconstruction is needed; two cells remain sufficient for two fixed values.
Kernel regressions cover mixed policies along every normal axis, analytic
Favre density covariance, affine and reversed scalar coordinates, and early
refusal of insufficient normal resolution.

The reused numerical storage remains 136 bytes per fluid cell for momentum
and 64 for one scalar, with actual-capacity admission. Plane reduction is
linear in cell count, and preparation adds no allocation. The caller retains
ownership of physical wall fluxes, molecular transport and the accepted step;
scratch is usable only after a successful preparation. This is a wall-resolved
planar closure, not a wall function, rough-wall model, arbitrary solid/open
boundary filter, nonuniform mesh method or DNS qualification.

`channel::components::ChannelDynamicComponents` fits one common component
diffusivity from that same current channel flow, alongside independent momentum
and heat. All species, including a closure species, enter one same-unit
normalization and joint least-squares objective. Raw component contractions
are summed before a single plane reduction, retaining physical relative
weights without a component-by-cell copy. Storage adds 80 numerical bytes per
fluid cell, independent of component count; preparation remains allocation-free.

Both physical and Favre-filtered component gradients use second-order
three-fluid-centre interior reconstruction at the walls. No composition wall
value or zero gradient is inferred: an impermeable Maxwell--Stefan flux can
balance composition and thermal-diffusion driving terms. The native owner must
retain its actual species-wall flux, absent-species rule, discrete zero-sum
mass flux, partial-enthalpy transport and physical diffusivity limit. This
component path therefore requires three normal fluid cells. Momentum and heat
retain their separate declared wall conditions; existing periodic arithmetic
and scalar wall APIs remain unchanged. Focused tests cover quadratic profiles,
same-unit basis/species-order invariance, trace weighting, a zero-flux thermal
diffusion balance with nonzero composition gradient, and bounded refusal/retry.

`Grid3::{channel_dynamic_sgs_workspace, prepare_channel_dynamic_sgs,
step_serial_with_channel_dynamic_sgs, step_pooled_with_channel_dynamic_sgs}`
connects this closure to the original isothermal lattice collision/streaming
owner, including D3Q27 cumulant collision. It infers the unique nonperiodic
wall axis, requires an all-fluid rectangle and passes the original stationary
halfway wall values. The two remaining axes must be periodic. The existing
uniform body force enters through the same Guo-centered macroscopic velocity.
Spatially varying or time-dependent forces use the `ForcedGrid3` methods below.

The explicit moving-wall entries
`prepare_channel_dynamic_sgs_with_moving_walls`,
`step_serial_with_channel_dynamic_sgs_and_moving_walls` and their pooled
counterpart derive both wall velocities from the complete native
`WallMomentumProbe3` and its prescribed velocity field. Every reflected link
on both walls is required, including a stationary plate when only one moves.
Each plane must have one uniform tangential velocity. Native topology,
finite-value, Mach and discrete impermeability checks remain authoritative;
partial probes and nonuniform drives refuse even when their net mass source
vanishes.

That exact wall pair enters the dynamic filter and the quadratic wall-relative
dissipation stencil. The prepared molecular-plus-SGS field then feeds the
original moving-wall stream, which returns its actual impulse and wall-work
receipt alongside the SGS report. Actuator work supplied to the fluid is the
negative of reported wall work. No new arrays, stochastic force, thermal
source or persistent boundary drive are introduced. Zero drive retains the
stationary arithmetic. The existing channel resolution and exploratory-use
limits apply; arbitrary moving geometry and turbulent-flow qualification
remain outside this API.

The local Smagorinsky, WALE and Vreman closures use the same complete-wall
admission through `prepare_channel_sgs_with_moving_walls`,
`step_serial_with_channel_sgs_and_moving_walls` and the pooled counterpart.
Their ordinary `SgsWorkspace3` is reused without a test filter or additional
allocation. The actual wall-relative full gradient drives both the closure
and its dissipation diagnostic, so translating the fluid and both plates
together does not introduce artificial wall shear. Source molecular fields,
independent collision rates, native impulse/work and atomic publication use
the same owners as the dynamic path. A single normal fluid layer is admitted
by the local quadratic stencil; there is no nonlocal-filter resolution claim.
Zero wall drives retain the original stationary local-SGS arithmetic. Focused
local regressions cover analytic Couette transport, co-motion, active-model
frame changes, native/pooled equivalence and failure/retry boundaries.

Five native regression definitions cover Couette dissipation and co-motion
for both lattices and all wall normals, active-fit frame changes, material
rates and original impulse receipts, serial/pooled and zero-drive identity,
and incomplete/nonuniform drive, cancellation and late-stream refusal/retry.
These new integration tests await native execution; the implementation
session's shared executor exhausted storage and memory before that build.

The channel and periodic native adapters share molecular-ratio admission,
collision conversion and dissipation reporting. Every preparation starts from
the original collision or supplied molecular field. A zero fitted viscosity
preserves that collision exactly; D3Q27 central/cumulant changes only the shear
rate while retaining independent trace and higher rates. Prepared collisions
feed the original serial or pooled transaction, so refusal or cancellation
publishes neither populations nor the accepted step count.

`ChannelDynamicSgsWorkspace3` includes the shared filter, current density and
velocity, and the collision field, with nominal payload
`N * (136 + 32 + size_of::<Collision>())` bytes and actual-capacity admission.
Preparation allocates no additional arrays. Native regressions exercise
wall/body-force transport, nonshear rates, local molecular fields, prepared
collision equivalence, serial/pooled identity and failed-trial retry. This
adds a native resolved-channel route for the cumulant solver, without claiming
a developed-turbulence accuracy/stability envelope or calibrated wall law.

### SGS with spatial and time-dependent force history

`ForcedGrid3` supplies local Smagorinsky/WALE/Vreman, periodic dynamic and
planar-channel dynamic preparation and stepping. Each family offers serial
and pooled steps that hold the accepted force, plus advances that accept the
next force endpoint. They reuse `SgsWorkspace3`, `DynamicSgsWorkspace3` and
`ChannelDynamicSgsWorkspace3`, their existing reports and the original or
explicitly supplied molecular collision field; previously augmented rates
must not be resubmitted as molecular transport.

Every closure reconstructs incoming physical velocity as
`u_n = (j_n + F_n/2)/rho_n` using the accepted per-cell force. The supplied
`F_(n+1)` never reinterprets that incoming field or its strain. Collision uses
`F_n`; endpoint admission and velocity readback use `F_(n+1)` through the
original force transaction. On a periodic domain the physical momentum
increment remains `sum_cells((F_n + F_(n+1))/2)`. Populations, accepted force
history and step count publish together only after successful completion;
closure/native refusal or cancellation retains all three accepted values.

Preparation adds no arrays or force-history snapshot and retains the existing
bounded workspaces and cancellation checkpoints. Ordinary `Grid3` paths keep
their original uniform-force arithmetic. Local solid/stationary-halfway-wall
stencils, all-fluid periodic/channel dynamic topology, independent nonshear
rates, molecular-ratio limits, byte limits and exploratory-use admission stay
in force.

The moving-wall entries on `ForcedGrid3` now use the same native wall probe,
correction and impulse/work receipt as `Grid3`, while retaining both body-force
epochs. They offer held-force steps and endpoint-force advances in serial and
pooled execution. Physical momentum changes by the trapezoidal body impulse
plus the receipt's `fluid_population_impulse`; the original wall-relative
mass terms and sign of wall work remain unchanged. The receipt is admitted
from the actual candidate buffers before populations, clock or next forces
publish. No second collision or population buffer is introduced.

Local and dynamic channel SGS compose with these entries through
`prepare_channel_sgs_with_moving_walls` and
`prepare_channel_dynamic_sgs_with_moving_walls`, plus their combined
`step_*`/`advance_*` methods. The current force frame, complete uniform
tangential wall pair and original molecular field determine each incoming
fit. A next force endpoint cannot alter that fit. The same prepared
collisions and exact wall drive then enter the source-aware native stream.
Existing channel topology, Mach, viscosity-ratio and exploration-only
admission remain authoritative. This adds no open-boundary model, thermal
energy owner or turbulent-DNS qualification.

Focused regression definitions cover uniform-force native equivalence,
changing-force momentum balance, force-frame and next-endpoint independence,
local/dynamic prepared-step equivalence, worker identity and atomic
refusal/cancellation/retry. Native compilation and execution were unavailable
in this implementation session after the shared executor disconnected.

### Resolved vorticity on native three-dimensional flow owners

`Grid3::kinematics_with_moving_walls` and the corresponding `ForcedGrid3`
entry bind the resolved gradient to the actual prescribed native wall drive.
They share moving-stream topology, Mach, tangential-motion and local mass-flux
admission. The wall velocity for each missing Cartesian neighbor comes from
that axial reflected link; diagonal contacts are not averaged into it.
Unselected links remain stationary, matching native streaming. This supports
spatially varying drives and fixed staircase solids, without a channel-only
restriction. A single ordered sweep of selected links uses fixed six-face
scratch; repeated sampling adds no heap allocation. The ordinary stationary
path retains its arithmetic, and spatial forcing retains the accepted epoch.
`MovingWallKinematicsError` separates native wall refusals from gradient
refusals. No moving geometry or curved-wall-distance model is introduced.
The moving-wall regressions in `tests/kinematics.rs` check both lattices and
all axes, co-translation, linear and quadratic shear, a one-cell gap,
spatial forcing, a nonuniform wall drive, and refusal/retry without mutation.
`Grid3::kinematics_with_finite_wall` binds the same query to the accepted
finite-inertia wall endpoint, reconstructing its translation/spin field rather
than reusing the previous stream's midpoint scratch. A regression advances the
coupled wall/fluid state and distinguishes these two velocity epochs explicitly.

`Grid3::kinematics`, `ForcedGrid3::kinematics`,
`BoundaryGrid3::kinematics` and `NsfGrid3::kinematics` expose the full resolved
velocity gradient independently of an SGS closure. Every adapter samples the
accepted half-force-corrected velocity once and reuses its original native
SGS/rheology stencil: centered periodic/interior differences, quadratic
halfway-wall reconstruction, and the owner's existing on-site or characteristic
open-face/reservoir convention. `ForcedGrid3` retains its accepted spatial
force epoch. NSF uses accepted characteristic traces and distinguishes
full-spacing reservoir ghosts from half-cell walls.

`kinematics::FlowKinematics3` defines `gradient[i][j] = du_i/dx_j`, signed
`curl(u)`, `|curl(u)|²/2` enstrophy density, full divergence,
`Q = (Omega:Omega - S:S)/2` rotation excess and `(omega dot grad)u`
stretching/tilting. Curl is twice rigid-body angular velocity; pure strain has
zero curl, and simple shear has nonzero curl but zero Q. Enstrophy is neither
density-weighted nor volume-integrated. These invariants are local resolved
quantities, not a coherent-vortex certificate or a full material vorticity
derivative: baroclinic generation, compressible dilation and viscous terms are
not folded into the stretching value. Native output uses lattice units;
`to_physical(dt_s)` divides the tensor by the actual positive timestep, giving
curl/divergence in s^-1 and quadratic invariants in s^-2. Cell size cancels.

`KinematicsWorkspace3::new(dimensions, maximum_bytes)` admits both numerical
buffers before allocation. Repeated sampling allocates no buffers, polls
cancellation at most every 64 cells and returns a canonical x-fastest field
with `None` for solid cells, independent of the tile-major D3Q19 storage.
Wrong workspace shape, invalid states, nonfinite invariants and cancellation
refuse without changing the flow. Scratch is consumable only from a successful
sampling call and is rebuilt on retry. No new collision, transport, closure,
resolution or physical-validation envelope is claimed.

`tests/kinematics.rs` checks G0 rotation/strain/shear/stretching identities and
SI conversion, G1 all-axis Couette and quadratic Poiseuille wall curl with
spatial forcing, a three-rung three-dimensional Taylor--Green curl convergence
ladder on both D3Q19 and D3Q27, analytic viscous enstrophy decay of the exact
two-dimensional Taylor--Green solution extruded in z, NSF moving-wall and
characteristic-throughflow sampling, solid exclusion and bounded refusal/retry.

`MaterialScaling::viscous_stress_from_lattice_kinematics` evaluates the shared
`fs-material::fluid_transport::viscous` law on that original resolved gradient,
divided once by the retained physical timestep. The retained source transport
supplies molecular shear viscosity; the caller supplies independent bulk
viscosity. Tensor stress, oriented shear/normal traction, dissipation and the
constitutive tangent use SI units and xx/yy/zz/xy/yz/zx tensor-shear ordering.
The companion `thermal_viscous_stress_from_lattice_kinematics` evaluates the
authored thermal Carreau-Yasuda law at the actual temperature and deviatoric
rate, retaining its source support, shift and derivatives. Dilation cannot
masquerade as shear thinning. These read-only observations allocate no buffers,
do not resample the retained transport, and do not insert forces or heat.
They are continuum constitutive stress on the native spatial stencil, distinct
from collision stress and discrete wall-link momentum exchange; no SGS stress,
thermodynamic pressure or additional validation claim is included.
`tests/material_shear_stress.rs` checks both lattices and all six shear-plane
orientations against prescribed Couette/Poiseuille profiles with moving walls
and current spatial half-force correction, signed tractions, physical power,
thermal law evaluation, rigid rotation, independent bulk loss and source refusal.

### Resolved vorticity on the native planar flow owner

`core2::Grid::kinematics` and `kinematics_with_moving_walls` expose the full
D2Q9 spatial gradient using `KinematicsWorkspace3` with `[nx, ny, 1]`. The
planar tensor has homogeneous z, signed `omega_z = dv/dx - du/dy`, enstrophy,
divergence, Q, and exactly zero stretching; the symmetric collision stress
cannot supply this rotation information. Each native physical velocity,
including gravity and spatial Guo half-force correction, is cached once.
Interior/periodic derivatives are centered second order; halfway boundaries
use the same quadratic wall/center/opposite-neighbor reconstruction as 3D.
Moving-wall arrays use the native streaming operator's wall-only full-grid
layout. Missing nonperiodic neighbors remain stationary halfway faces.
`to_physical(dt_s)` supplies SI scaling with the accepted state's timestep.

The adapter samples Fluid/Wall topology and explicitly refuses Gas/Interface
cells; it does not invent wall conditions at a free surface. Curved-link
bounce-back is outside this halfway-wall geometry. One-cell periodic axes
are homogeneous; two-cell periodic axes refuse an unresolved aliased
stencil. Walls return `None`. Repeated calls allocate no buffers, poll at
most every 64 cells, and return data only on success. Invalid mutable grid
buffers, invalid wall velocities and cancellation leave flow state unchanged.
NSF inviscid vorticity budgets also refuse two-cell aliased source stencils,
matching their existing body-acceleration-curl query.

`tests/kinematics2.rs` checks analytic rotation/strain/shear/dilation with
spatial forcing and Galilean shifts; both orientations of linear/quadratic
moving-wall profiles including one-fluid-cell gaps; stationary exterior
Poiseuille walls; a three-rung periodic Taylor--Green curl convergence ladder;
and viscous Taylor--Green enstrophy decay through native D2Q9 stepping.
Free-surface, cancellation, invalid state and underresolved stencil refusals
are covered separately. No collision, forcing or transport law is changed.

### Source-resolved rheology with local and dynamic SGS

`lattice3::rheology::RheologySgs3` selects a local closure, periodic dynamic
fit or channel dynamic fit together with its existing reusable workspace.
The native `step_*_with_material_rheology_and_sgs` methods on `Grid3` and
`ForcedGrid3`, and the latter's endpoint-force advances, first resolve the
original material producer at each incoming physical shear rate. That exact
source-derived molecular field supplies SGS viscosity-ratio admission and
molecular dissipation. Each call starts from fresh molecular transport;
previously augmented collisions are never recycled as material coefficients.

Source records keep molecular relaxation times and original transport
responses. The accepted `RheologySgsStep3` also retains the actual applied
molecular-plus-SGS collision field and typed local/dynamic diagnostics.
Independent nonshear rates, reference density and source-domain checks keep
their original owners. Current F_n sets both source and SGS preparation;
F_(n+1) enters only the native force-history advance.

The source and SGS workspaces are reused without new adapter-owned arrays.
Copying the applied collision snapshot occurs in bounded chunks before native
stepping. Successful population publication is followed only by the original
infallible source-report swap. Refusal or cancellation retains the previous
accepted source epoch, populations, force history and clock. Source callbacks
retain responsibility for their own allocation and external side effects.

This is an explicit resolved-shear generalized-Newtonian LES composition.
It retains each closure's stationary-wall/periodic/channel topology and
exploration-only policy. It does not close correlations between unresolved
strain and nonlinear molecular stress, or add constitutive memory, heat
evolution, moving-wall rheology or turbulent-flow qualification. Focused
regression definitions cover original-source and zero-SGS equivalence,
nonzero local/dynamic transport, incoming force timing, source epochs,
worker identity and rejected-step retry; native execution was unavailable
after the shared executor disconnected.

### Moving/open D3Q19 boundaries and turbulent wall loads

`BoundaryGrid3::{sgs_workspace, prepare_sgs, step_with_sgs,
step_pooled_with_sgs}` apply the same explicit Smagorinsky/WALE kernel to the
existing tiled boundary owner. Molecular transport comes from its constructor
collision law on every step; no previous eddy viscosity is compounded. The
cached velocity includes both uniform and per-cell Guo force corrections.
Workspace collision entries follow native `tile * TILE_CELLS + lane` order.

Interior and periodic derivatives are centred. Halfway-wall derivatives use
the actual moving-wall value and the streaming owner's stationary wall/open
seam and corner conventions. On-site velocity faces use a quadratic one-sided
derivative through two interior fluid cells, with a linear fallback when the
second is solid; pressure faces use their zero-normal-velocity-gradient rule.
The original regularized open-face population reconstruction remains the
boundary model. No synthetic turbulent inflow, transparent outlet, developed
jet qualification or improvement of that model's existing accuracy envelope
is implied by adding SGS.

Optional canonical wall-link selections return `BoundarySgsReport3` with the
actual same-step wall impulse, angular impulse and work, alongside pre-step
shear-dissipation diagnostics. Both serial and pooled paths reuse the original
collision, link streaming, open reconstruction and wall-momentum kernels.
Closure sweeps allocate no buffers after workspace creation; publication uses
the existing reusable stream destination. Input, ratio, collision, load and
cancellation refusals publish neither new populations nor a load result.
Ordinary stepping does not enable SGS. The zero coefficient retains original
collision bits and ordinary trajectories; independent higher rates survive.

`tests/boundary_sgs.rs` covers moving/open/periodic zero controls, analytic
Couette and quadratic inlet gradients, force-frame consistency, actual
SGS wall loads and worker-count identity, plus refusal/retry behavior. These
are focused integration regressions, not a turbulent-flow validation campaign.

`BoundaryGrid3::{dynamic_sgs_workspace, prepare_dynamic_sgs,
step_with_dynamic_sgs, step_pooled_with_dynamic_sgs}` also connect the shared
Favre Germano--Lilly fit to this native collision/streaming and wall-load owner.
`BoundaryDynamicSgsConfig3` explicitly selects a fully periodic box or a planar
channel with one impermeable moving-wall pair and two periodic directions.
Channel fitting uses the existing tangential test filter and separate
wall-parallel plane averages, including both current physical wall velocities.
Open faces, interior solids, extra wall pairs and mismatched workspaces refuse;
there is no truncated-filter fallback at unsupported geometry.

The reusable workspace includes native macros/collisions, x-fastest filter
inputs and the original dynamic filter. Its complete numerical byte budget is
checked before allocation and against actual capacities. Each preparation
reconstructs the current uniform-plus-local half-force velocity, fits once,
and maps eddy viscosity back to native tile order. Coefficient and
eddy/molecular limits refuse instead of clipping; the constructor's molecular
collision and independent higher rates are retained. Dynamic diagnostics,
instantaneous shear dissipation and optional same-step wall loads publish with
the accepted step. `tests/boundary_dynamic_sgs.rs` covers shared-filter oracles,
all-axis moving walls, force centering, laminar Couette, physical wall-momentum
balance, serial/worker parity, bounded workspace and refusal/retry. These APIs
remain exploratory, with the existing D3Q19 tied shear/bulk rate and no wall
function, source-material evolution or turbulent-DNS qualification.

### Spatial and evolving D3Q19 inlet velocity

`BoundaryGrid3::set_velocity_profile` samples an explicit vector trace on an
existing velocity face in canonical z/y/x order. Coordinates are lattice-node
indices; the caller owns physical units, time and the coherent incident field.
The same accepted trace feeds the original open-population reconstruction in
ordinary, pooled and SGS steps. Sampling does not advance populations; the
trace supplies the next boundary endpoint. Uniform face specifications remain
the fallback, and `clear_velocity_profile` restores them.

`InletFlux3::MeanNormal` applies one uniform normal correction to attain a
declared signed outward-normal mean velocity (inflow is negative on either
end). Tangential components are unchanged. This constrains the prescribed
volumetric velocity trace, not density-weighted or reconstructed population
mass flux. Partially reconstructed wall seams refuse this control. Every raw
and corrected sample must satisfy the existing velocity bound. Invalid source
values, geometry, cancellation or callback errors preserve the previous
complete trace and populations. Two face buffers are reused after admission;
the first successful profile locks geometry as ordinary stepping does.

`tests/boundary_profiles.rs` covers all six orientations, normal-mean control,
exact uniform/clear trajectories, coherent evolving serial/pooled SGS input,
buffer reuse and refusal atomicity. `examples/sgs_spatial_inlet.rs` drives the
actual WALE boundary step with two advected, divergence-free periodic Fourier
modes and reports transport, shear dissipation and mass. It demonstrates the
spatial inlet path; it supplies neither a measured turbulent spectrum nor a
developed turbulent inflow or a newly qualified transparent outlet.

### Full-tensor Mann inlet observations on D3Q19

With `turbulent-inflow`, `BoundaryMannInflow3` prepares an explicit SI/lattice
frame and point or exact square-face observation for the original
`BoundaryGrid3::set_velocity_profile` owner. The same accepted velocity trace
feeds ordinary, pooled and SGS collision/streaming. Node positions use the
existing on-site boundary coordinates, not exterior NSF ghost centers. The
full Mann cross-component spectrum comes from fs-atmo; there is no variance
restoration, inferred mean or independent random generator. The prepared
`MannBoxFilter` binds the complete immutable model and uses no gain storage for
point observations. Face gains are prepared once under an explicit byte bound.

One borrowed source selector consumes a Taylor-frozen snapshot with convection
in grid Cartesian axes at the authored physical time, a current
`FixedStepMannField` endpoint, or its accepted stationary
`MannStepAverage`. Live and averaged sources retain their own advection exactly
once, must use the declared CFD dt, and are never advanced by the adapter. An
interval mean must end at that source's accepted endpoint. It averages velocity,
not nonlinear momentum/energy fluxes; it adds no higher-order boundary claim.
Spatial and advective phase bounds must not exceed pi, and active-node/mode work
is bounded before sampling. The temporal bound does not bandlimit OU innovations.

Arbitrary Mann wavevectors do not wrap periodically: transverse periodic faces
refuse. Native transverse-wall/open-seam rules remain, including `InletFlux3`'s
mixed-seam refusal for mean-normal control. Unconstrained sampling retains the
finite aperture's fluctuating mean. Source, profile and population state remain
unchanged on refusal; successful application only publishes the original profile.
A later fluid-step refusal retains that accepted control. No viscous inlet wall
layer, developed turbulent-flow calibration or transparent outlet is claimed.
The focused `boundary_mann_inflow` tests exercise all six face mappings, exact
spatial averaging and native serial/pooled SGS, original endpoint/interval-mean
consumption, and source/clock/resolution/topology/work/cancellation refusals.

### Transient finite-volume momentum and heat

`FlowStepRecord::kinetic_energy_j` uses the original staggered face velocities
and their fluid half-cell masses. Interior faces own the adjacent half cells;
exposed inlet and outlet faces own only the inner half cell, excluding mirrored
ghost volume. Stationary solid and impermeable faces contribute zero. Squaring
cell-centred averaged velocities would discard resolved kinetic storage, even
for a divergence-free circulation. Both the ordinary and source-material
marches now use the same allocation-free face norm. In closed domains its
interior inertia and summation order agree with the existing physical
viscous-heating account. This diagnostic does not promote the ordinary open
flow discretization to a total-energy-conservative scheme or create heat from
numerical kinetic losses. Analytical vortex, graded-boundary, uniform-motion
and solid-occupancy regressions exercise the norm and its march integration.

`conjugate::simple_unsteady_with_sgs` and `march_conjugate_with_sgs` select the
same full-gradient Smagorinsky/WALE kernel through explicit `FvSgsConfig`.
The physical filter width is the declared factor times `cell_volume^(1/3)`;
neither `SimpleConfig` nor existing march entry points enable SGS. Acceptance
use and simultaneous equilibrium LVEL refuse. A zero coefficient keeps the
ordinary laminar arithmetic and empty closure state. These raw entries use
the declared constant molecular viscosity. The source-material entry below
retains a separate evolving molecular-viscosity field.

One reusable velocity workspace freezes the gradient field before each
nonlinear momentum sweep. Normal derivatives use the owning staggered faces;
off-diagonal derivatives use quadratic nonuniform centred stencils. Domain
walls/inlets supply their actual prescribed velocity at the face, voxel solids
supply stationary no-slip values, and symmetry/outlet faces impose zero
tangential derivative at the face rather than at the adjacent cell centre.
The closure supplies the existing complete symmetric stress and is refreshed
again after the step's pressure projection. Invalid gradients or an exceeded
eddy/molecular ratio cap refuse; no clipping or silent laminar fallback occurs.
The filter widths are cached once; gradient/closure sweeps allocate no buffers
and poll cancellation every 64 cells.

Every coupled energy step consumes `rho cp nu_t / Pr_t` from that step's
projected flow. `SimpleConfig::turbulent_prandtl` independently declares the
positive finite heat-diffusion ratio (default 0.9), without changing the
momentum model or molecular material. `FvFlow` retains it for subsequent heat
solves and exports. The cooling scene exposes `solver.turbulent_prandtl` for
LVEL and transient SGS and reports the actual retained value. The constant
ratio remains an authored model assumption, not a universal material law;
see [NASA Wind-US Prandtl controls](https://www.grc.nasa.gov/www/winddocs/user/keywords/prandtl.html).
A caller-supplied nonzero eddy-conductivity field refuses instead of being
combined or overwritten. Resolved-SGS marches select
`ThermalSetup::eddy_wall_treatment = EddyWallTreatment::Molecular`: fixed and
film domain walls, ordinary fluid/solid interfaces, compact package surfaces
and explicit stationary thermal joints retain their molecular fluid half-cell
resistance. Fluid/fluid mixing and open-face diffusion still use the effective
conductivity. This is a resolved-wall boundary condition, requiring adequate
near-wall resolution, not a new coarse-grid heat-transfer correlation. The
default `EffectiveConductivity` policy preserves existing equilibrium LVEL
behavior; an explicitly requested equilibrium thermal wall function and the
resolved molecular policy cannot be combined.

The accepted molecular wall coefficients remain separate from effective cell
conductivity. Source-contact surface queries, nonlinear solid half-cell laws,
compact links, heat observers and material endpoint flux corrections all use
the same accepted wall coefficients. In particular, a molecular conductivity
change cannot disappear from the endpoint audit merely because an opposing
eddy-conductivity change keeps their sum constant. No SGS dissipation heating, turbulent scalar-variance
equation, equilibrium thermal wall law, calibrated turbulent channel profile,
or resolved turbulent spectrum is claimed. The existing first-order thermal
coupling and power-law/upwind convection limits still apply; the caller must
establish temporal and spatial convergence for the intended exploration.

Focused `conjugate::simple::sgs` tests cover graded Couette and voxel walls,
quadratic Neumann gradients, native MAC normal derivatives, WALE's pure-shear
zero, physical Smagorinsky viscosity, exact disabled flow identity, admission
and cap refusals, and a moving-wall flow/heat step compared against an
independently supplied turbulent-conductivity energy solve.
The thermal face-contact battery additionally checks graded molecular wall and
film resistances, retained interior/open-face mixing, ordinary and compact
solid interfaces, and source queries at the molecular interface surfaces.
An endpoint regression retains the changed wall heat at equal total cell
conductivity; the heated moving-wall transient checks the automatic SGS policy.

`MaterialChtFluid::march_natural_convection_with_sgs` uses the same native
transient owner with source-dependent molecular viscosity and conductivity.
Each accepted thermal endpoint refreshes the original source field; each
nonlinear momentum sweep checks `nu_t/(mu_local/rho)` against the explicit
cap and applies `mu_local + rho*nu_t` in the reciprocal symmetric stress.
Material edges retain harmonic transmission, and no-slip half-cell wall
traction remains molecular. Source snapshots never absorb the eddy field,
so a constitutive refresh cannot compound turbulent viscosity. The projected
eddy field supplies exactly one separate gradient heat flux with reference
`rho*Cp` and the original `flow.turbulent_prandtl`. Neither molecular nor modeled mechanical loss is
deposited as heat. Source density/Cp bounds, first-order feedback, cancellation
and exploratory SGS restrictions remain in force; equilibrium LVEL and a
second authored eddy-conductivity field still refuse. The focused source/SGS
tests cover affine shear, reciprocal edges, molecular walls, local viscosity
caps, and refresh without compounding. The `material_conjugate` heated-cavity
regression checks live viscosity/momentum response, zero-closure identity,
heat conservation and an independent once/omitted/twice eddy-heat replay.

`march_natural_convection_with_enthalpy_and_sgs`,
`march_natural_convection_with_material_solids_and_sgs`, and
`march_natural_convection_with_enthalpy_and_material_solids_and_sgs` extend
that same explicit native SGS owner through fluid source h(T), solid source
Cp(T)/k(T), and joint fluid/solid endpoint evolution. Existing signatures
retain their behavior. One accepted eddy field and projected flow are reused
through every nonlinear thermal trial; momentum is not stepped again. Eddy
conductivity is added only to fluid cells, once, and never enters molecular
snapshots or source storage. Fluid h(T) owns actual storage and conservative
donor heat; solid source integrals, contacts, reservoirs, radiation, local
nonlinear gates, work bounds and atomic refusal remain in force. In the source-enthalpy
entries, gradient diffusion uses the current source caloric derivative and
the original `flow.turbulent_prandtl`. The linear-Cp fluid entry retains its declared reference conversion. These entries select no new dissipation-heating or
wall model. First-order feedback and exploratory SGS restrictions apply;
zero SGS recovers the corresponding original march. The caloric-solids
regressions specify joint source-mu(T) momentum/heat response, independent
enthalpy conservation, constant-Cp reductions and late source refusal; these
new whole-crate tests remain unexecuted under the current host memory limit.

For the explicit source-enthalpy SGS marches, each native caloric trial
rebuilds `k_t = rho0 * Cp_source(T_trial) * nu_t / Pr_t` from the immutable
original Prandtl control and
accepted eddy field. It never compounds a prior thermal trial's conductivity.
Only fluid cells receive this conversion; original h(T) continues to own
storage and donor-cell advective heat. After each candidate, native individual
face and net cell/junction conductive defects compare applied and endpoint
source Cp against the existing `maximum_flux_defect_w` tolerance. Both sides
reuse the actually accepted fixed/source contact resistances and original
compact-link laws. The signed boundary/thermostat correction enters the same
heat ledger exactly once. A globally small energy residual cannot hide a large
local heat-flux mismatch.

`maximum_sgs_heat_flux_defect_w` records the admitted local defect;
`SgsHeatNotConverged` refuses an exhausted endpoint without publishing state.
No SGS and matching constant source/reference Cp retain their original
reductions. An independent moving-wall two-cell regression checks analytic
SGS shear, true variable-Cp enthalpy storage, harmonic half-cell heat flux,
thermostat accounting and the explicit one-pass flux refusal. Both source
fluid-only and joint fluid/solid caloric owners use this closure. First-order
coupling, exploratory SGS scope and separate dissipation/wall-model selection
remain unchanged. All changed sources pass syntax parsing; full native
integration execution remains unverified under the host resource limits.

## Compressible Navier-Stokes-Fourier hybrid (IRP M3.2-W3a, `nsf3`)

`nsf3::NsfGrid3<L: Lattice3<27>>` advances density, momentum and total energy of an ideal gas on the D3Q27 lattice.
The formulation and its conservation class are stated here because the bead requires the choice to be disclosed.

**Formulation.** Mass and momentum are carried by the populations of an ISOTHERMAL lattice: equilibrium at the lattice
reference `theta0 = 1/3` (`THETA0`, the lattice `c_s^2`), a regularised collision whose shear rate is
`1 / (mu(T) / (rho theta0) + 1/2)`, a trace (bulk) rate derived from a SOURCED bulk-to-shear viscosity ratio
(`NsfSpec::bulk_to_shear_ratio` with a mandatory citation string; `1/w_b = 1/2 + 1.5 ratio mu / (rho theta0)`, from the
bulk viscosity `(2/3) rho theta0 (1/w_b - 1/2)` of the lattice, checked against the measured sound attenuation of the
scheme) and the third-order non-equilibrium part rebuilt recursively from the deviatoric stress. The thermodynamic pressure difference `p' = rho (theta - theta0)` of `p = rho theta` enters as the Guo body force
`F = -grad p'` in conservative face form with the half-force velocity `u = (sum c f + F/2) / rho`; the first-moment
non-equilibrium `-F/2` that the regularised projection drops is added back so that the momentum receives exactly `F`.
Total energy `E = rho (c_v theta + |u|^2 / 2)` is a conservative finite-volume field on the same cells with the face flux
`(E + p) u - sigma . u - kappa grad T` (`sigma` from the second moment of the pre-collision populations about `feq`, shifted by
Guo's `(F u + u F) / 2` as described below; fourth-order face values away from walls, two-point beside walls). **Conservation class: total-energy conservative.** Interior face fluxes are single-valued,
so the sum of `E` changes only by boundary heat/work and explicit source supplies, to round-off; a source-free closed box
with stationary adiabatic walls conserves total energy. Mass and
momentum are the lattice's. The temperature is `(E - KE) / (rho c_v)` and feeds back only through the force.
Coupling in time: trapezoid on the flux divergences of the two levels with `COUPLING_PASSES` fixed
predictor-corrector passes. `NSF3_BIT_SEMANTICS_VERSION` moves with any change of order, coefficient or coupling.

**Optional endpoint accuracy control (`NSF3_BIT_SEMANTICS_VERSION` 9).**
`with_coupling_control` / `set_coupling_control(Some(NsfCouplingControl))`
adds four explicit positive tolerances and `maximum_endpoint_evaluations`.
`None` remains the exact original three-pass default, with no probes or changed
corrector arithmetic. The minimum budget is four: three original correctors,
then a probe of that same retained candidate. A failed probe permits another
original corrector only while a further verification evaluation remains. The
report exposes `corrector_passes`, `endpoint_evaluations` and `residual_probes`;
a loose control reports 3/4/1, while total native energy-flux evaluations equal
correctors plus probes. Mass collision/streaming still occurs exactly once.

The probe evaluates the candidate's actual stored E, temperature, physical
velocity and Guo force. It checks the thermodynamic/source force and its
half-force velocity separately, re-evaluates source caloric closure, then
recomputes the native local conservative energy equation with actual endpoint
stress, material/SGS transport, boundary data and already sampled supplies.
Gravity and changing body acceleration use their original work quadratures;
manufactured callbacks are not resampled. The SGS/enthalpy-contact probe
resolves its boundary at the retained physical velocity. No candidate force,
energy, face or divergence history is overwritten by the probe.

All four residuals are infinity norms. Force uses `rho*THETA0`, velocity uses
`sqrt(THETA0)`, and both caloric and energy defects use
`rho*THETA0/(gamma_ref-1)`. In physical units these correspond per unit mass
to `R*T_ref/dx`, `sqrt(R*T_ref)` and `cv_ref*T_ref`. These fixed native reference
scales remain defined at zero force/velocity and never grow with an arbitrary
energy datum; unrepresentable scales refuse rather than use an epsilon floor.
The energy norm is a local equation defect, distinct from the existing global
telescoping conservation ledger. The accepted conservative candidate retains
that ledger's exact flux ownership.

`last_coupling_report` gives the four accepted residuals, actual work counts,
worst component and largest residual/tolerance ratio. Budget exhaustion returns
`NsfError::CouplingNotConverged` with the requested controls and final unaccepted
report. Refusal/cancellation publishes no populations, drive, thermal/flux
history, report, clock or ledger. Characteristic interval copies retain both
controls and the latest accepted report. Material caloric inversion can impose
a residual floor; the control does not relax it or widen any physical domain.
This is an accuracy-control capability, not a finding that the fixed-pass
method fails, a truncation-error estimator or new physical validation.
Focused tests cover the published-state residual, untouched three-pass
acceptance, work-budget/cancellation rollback, unchanged source sampling,
reference-normalized zero flow, source caloric/formation invariance and
SGS thermal-contact interval ownership. Their native execution remains pending
the host dependency-build resource limit.

**Stress of the energy flux (found by the Mach map; `NSF3_BIT_SEMANTICS_VERSION` 3).** Guo's scheme stores the populations half a
force source behind the physical ones (`f = f_phys - S / 2`, and the source has second moment `F u + u F`), so the viscous
stress that enters `sigma . u` is `-(1 - w_s / 2)` times the deviatoric part (and `-(1 - w_b / 2)` times the third of the trace)
of `Pi_neq + (F u + u F) / 2`. The first version took `Pi_neq` from `f - feq` alone, so the thermal-pressure force `F` acting on a
flow `u` appeared as a spurious viscous stress: an entropy wave advected at Mach 0.2 decayed slower than the quiescent wave by
-2.2e-2 (-0.57 Mach^2, equal to `-gamma (gamma - 1) Mach^2`, and unchanged from 32 to 64 points per wavelength, so not a
truncation error); with the shift the defect is +2.3e-3 on 32 points and +8.8e-4 on 64. Removing `sigma . u` altogether gave
+9.8e-4, which located the term. The momentum equation was never affected (the lattice collision carries its own stress); the
Taylor-Green numbers of the battery did not move (kinetic-energy decay rate over `4 nu k^2` 0.99945 on 16 and 0.99964 on 32
points; accumulated heating +1.388e-2 and +1.138e-2).

**Force projection by viscous mode (`NSF3_BIT_SEMANTICS_VERSION` 4).** The Guo
source's isotropic second moment now uses `1 - w_b/2`; its deviatoric second
moment retains `1 - w_s/2`. Previously both used the shear prefactor, creating
an artificial isotropic stress during acceleration whenever the bulk and shear
rates differed. The correction has zero mass and momentum moments and adds no
allocation. The local collision regressions `guo_source_respects_bulk_relaxation_for_uniform_acceleration`
and `bulk_and_shear_stress_relax_independently_in_forced_flow` verify independent
trace/deviatoric relaxation against prescribed physical moments, including
unequal rates, nonparallel velocity/forcing, and the equal-rate limit. These
local identities do not extend the existing Mach or transport envelopes.

**Reservoir density and forced stress (`NSF3_BIT_SEMANTICS_VERSION` 7).**
Incoming reservoir links now extrapolate `post - feq(rho,u) - S(u,F)/2`
and add `feq(rho_boundary,u_boundary) + S(u_boundary,F)/2`. Post-collision
velocity is recovered with the negative half-force correction. This retains
the physical shear/bulk nonequilibrium without the artificial
`F tensor F * (1/rho_boundary - 1/rho)/4` stress caused by shifting equilibrium
velocities. Each encountered cell/face builds its ghost populations once;
existing face selection, streamed mass and finite-volume energy accounting
remain the owners. Zero-force reconstruction retains its previous bits.
Focused tests cover unequal densities and viscous rates, all six prescribed
faces, worker counts, and native mass/energy closure. This changes forced
reservoir trajectories; it does not qualify a wider flow or acoustic regime.

**Physical initialization and boundary rebinding (`NSF3_BIT_SEMANTICS_VERSION`
5).** Initial populations now use `feq(rho,u) - S(u,F)/2`. A shifted-velocity
equilibrium alone introduced an artificial `F tensor F/(4 rho)` second moment.
Prescribed reservoir changes and characteristic boundary changes share the
same complete Guo half-source rebase, preserving physical normal/shear stress
as well as mass, velocity and total energy. The shared source helper allocates
no memory. G0 regressions exercise nonparallel pressure/body forcing, both
viscous modes, every boundary orientation and reservoir attach/change/remove.

Calorically perfect gas admission also checks the actual local kinematic
viscosity `mu(T)/rho` at construction, boundary admission, predictors and accepted
endpoints. Checking the reference-density transport table alone admitted dilute
states outside the existing verified viscosity limit. Such states now refuse
transactionally, with no clipping, extra material-cache allocation or expansion
of the declared flow regime. A focused rarefaction regression checks rollback.

**Why not a temperature-corrected equilibrium (measured).** The first formulation put `theta - theta0` into the equilibrium
(Hermite through third order). Its Bloch-linearised one-step matrix (28 x 28: 27 populations and the energy) about a
quiescent uniform state has a REAL eigenvalue of 1.4 to 2.0 at oblique wavevectors near the Nyquist limit at air
viscosity (dx 20 um): Taylor-Green at Mach 0.03 left the temperature window in 55 to 110 steps and a quiescent gas with
1e-5 noise in under 20, while one-axis acoustic tests were stable. Its eigenvector is the rest population and the
energy: at fixed `E`, `p = (gamma - 1)(E - KE)` does not depend on `rho`, so a density perturbation enters the
equilibrium with +2.5 on the rest population and -2 on the corners and the trace relaxation amplifies it. The eigenvalue
grew with the bulk rate (the physical bulk viscosity of air needs the largest rate) and was not removed by the flux order,
upwind dissipation, Adams-Moulton order, 1 to 8 coupling passes, the time level of the collision temperature, constant
viscosity, an internal-energy equation or dropping the viscous work. The energy-free equilibrium used here has
`|lambda| < 1` on the same scan (worst 0.9977 at the physical bulk rate).

**Units and receipts.** `theta = T / (3 T0)`; `dt = dx / sqrt(3 R T0)`; lattice density 1 is the reference density of
`GasState::try_new(gas, T0, p0)`; `mu(T)` and `kappa(T)` are tabulated from `GasState` receipts across the DECLARED
temperature window (`transport_points`, linear interpolation) and converted with `mu_lat = mu dt / (rho0 dx^2)`,
`kappa_lat = kappa dt / (R rho0 dx^2)`. No transport value is a literal. `NsfLaw::energy_unit_j()` converts lattice
energy to joules; the API reports joules.

**Bulk viscosity (sourced, dx-dependent like the shear rate).** The ratio `mu_b / mu` is a declared, cited input, not a lattice
constant: dry air near 293 K, 0.73 (M. S. Cramer, Phys. Fluids 24, 066102 (2012), doi 10.1063/1.4729611: "mu_b about
0.73 mu for N2 at 293 K", from the acoustic absorption of Prangsma, Alberga and Beenakker; shock-thickness value 0.67;
curve fit `mu_b = 0.778e-5 (T/200 K)^1.376 kg/(m s)` for N2, 77 to 1073 K). A Rayleigh-Brillouin measurement at
hypersonic frequency gives 1.527 / 1.826 = 0.84 for air at 295.8 K (Z. Gu and W. Ubachs, arXiv:1409.3014, Table II), so the
range carried is 0.67 to 0.84; the audio-frequency N2 value is used because instrument acoustics is far below the
relaxation frequency. The ratio is taken temperature-independent over the declared window (Cramer: it rises gradually with
T). The lattice rate therefore follows dx exactly like the shear rate: `w_b` is 1.9922 at dx 50 um (lattice viscosity 6.0e-4) and 1.9804 at dx 20 um (1.5e-3).
Stability scan (Bloch, 28 x 28, eight oblique modes, ratio 0.73): max |lambda| = 1.0000, 1.0000, 0.9999, 0.9998, 0.9984,
0.9953 for lattice viscosity 1e-5, 1e-4, 5e-4, 1.5e-3, 1e-2, 3e-2 and 1.146 at 0.1; ratio 0.84 at 5e-4 gives 0.9999.
The verified envelope is therefore `MAX_VERIFIED_LATTICE_VISCOSITY` = 0.03 and a finer dx is REFUSED (`NsfError::Viscosity`),
not compensated with inflated bulk viscosity.
Round-off scan (the four-digit values above hide the last digits). Complex-step Bloch matrix (`Im f(x + i 1e-30 e) / 1e-30`,
no subtraction error, so the 28 x 28 matrix is exact to machine precision) of the numpy prototype of this scheme (same
algorithm, not the Rust code), at EVERY wavevector of a 16 x 16 x 4 periodic box (242 modes: Nyquist and the short axis
included), ratio 0.73, eight coupling passes, self-conjugate Nyquist modes projected with 1/N. `max |lambda| - 1`:
-2.2e-16 (lattice viscosity 1e-5), +8.9e-16 (1e-4), +6.7e-16 (5e-4), -3.3e-16 (1.5e-3), -1.1e-16 (1e-2), -2.2e-16 (3e-2).
The extremum is a self-conjugate Nyquist mode whose eigenvalue is exactly neutral (`lambda = 1`), so +8.9e-16 and +6.7e-16 are
3 to 4 ULP of eigen-solver round-off, not growth; were it real it would be 9e-9 relative amplitude after 1e7 steps, nine
orders below the 1e-6 energy-ledger band. The smallest-k modes damp at -1.562e-6, -1.562e-5, -7.809e-5, -2.342e-4, -1.558e-3
and -4.656e-3 per step for the same six viscosities. No range is refused beyond the 0.03 envelope.

**Refusals.** `new` and `step` refuse: a gas state outside the receipt window (T in [50, 2000] K: 50 K admitted, 49.999 K
refused); a temperature outside the declared run window; more cells than `maximum_cells` (cap admitted, cap + 1
refused); a local Mach number `|u| / sqrt(gamma theta)` above `maximum_mach`; `|theta / theta0 - 1|` above
`maximum_temperature_ratio`; non-finite or non-positive density, energy or temperature; malformed bounds, bulk rate
empty citation or non-positive ratio, a lattice viscosity above the verified envelope (above), zero workers; a declared
`maximum_mach` above `MAX_QUALIFIED_MACH` = 0.2 (`NsfError::MachBound`: 0.2 admitted, 0.2 x (1 + 1e-9) refused; the bound is the
largest Mach number of the Galilean-invariance map in `tests/nsf3_qualification.rs` and nothing is claimed above it); a wall
velocity on a periodic axis, non-finite, wall-normal or above the Mach bound (`set_edge_velocity`). A failed `step` leaves the published state unchanged.

**Determinism.** Every cell pass is a pure function of read-only arrays written by `par_fill` over fixed contiguous
chunks; every reduction is serial in cell order. The result is bitwise independent of `NsfSpec::workers`
(`nsf3_bitwise_worker_count_invariance`: 1, 2, 3 and 7 workers). No transcendental function is evaluated in the module
(the receipts use the repository's `det` functions). No wasm32 build or golden hash is claimed.

**Moving walls.** `set_edge_velocity(axis, high, u)` gives the bounce-back wall at the low or high edge of a non-periodic axis a
tangential velocity: each link that crosses it gains the momentum `2 w_q rho (c_q . u_w) / theta0` (Ladd), and the energy flux
through the wall face gains the viscous work `-sigma_{a b} u_{w b}` with `sigma` of the cell beside the wall (exact for a uniform
stress, first order in `dx` for a developing flow). Solid cells inside the domain stay at rest, a velocity with a wall-normal
component is refused (no mass crosses a wall), and the speed is held to the declared Mach bound.
`wall_heat_j` retains Fourier heat alone; `wall_work_j` separately retains the
stress work of each wall at the same native face and time quadrature. Both are
positive into the fluid, with work negative when the wall extracts energy.
Heat and work are accumulated independently before their possible cancellation
in the total-energy flux. Their history refreshes with the original boundary
binding and publishes only with accepted steps. No extra field sweep, source
query or numerical allocation is introduced.

The internal-energy dissipation account subtracts thermal wall heat, without
subtracting external mechanical work as though it were heat. A steady
isothermal Couette flow can therefore retain positive wall work, equal
negative wall heat and positive viscous heating while net boundary energy is
zero. Version 8 records this corrected physical accounting; original combined
energy-flux arithmetic, fluid trajectories and transport envelopes are unchanged.
The local stress closure still has the original developing-flow boundary
accuracy and does not assert exact equality with bounce-back link impulse work.
Focused G0/G1 regressions exercise all six wall orientations and both work
signs, cancelling heat/work, an evolving adiabatic Couette field and refused
wall-law rollback. The independent steady temperature-profile qualification
remains `couette_viscous_heating_profile`.

**Spatial thermal walls (`NSF3_BIT_SEMANTICS_VERSION` 10).**
`with_thermal_walls` / `set_thermal_walls` attach a sparse list of
`NsfThermalWallFace { cell, axis, high, thermal }` overrides on the original
impermeable FV faces. A face belongs to an adjacent fluid cell and either a
solid neighbor or a non-periodic domain edge. Different hot/cold walls,
localized heaters and ambient cooling can therefore drive the native coupled
flow/total-energy equations. Unlisted faces retain `NsfSpec::wall`; an empty
list preserves the original global wall arithmetic and trajectory. Periodic
fluid neighbors, solid donor cells, duplicate/out-of-range faces and open
reservoir/characteristic faces are refused. Opening a face already carrying a
thermal override is also refused; ownership cannot be silently replaced.

`NsfWallThermal` admits adiabatic, isothermal [K], prescribed inward heat flux
[W/m²], and `Film { coefficient_w_m2_k, ambient_temperature_k }` conditions.
Positive heat enters the gas independently of face orientation. Isothermal
faces retain the original half-cell discretization and actual cell pressure:
`q_in = k((T_cell+T_wall)/2, p_cell) (T_wall-T_cell)/(dx/2)`.
The finite film solves this same molecular law together with
`q_in = h (T_ambient-T_wall)`, so the half cell and ambient film are in series.
It does not use conductivity at an assumed bulk/ambient mean. The native gas
source supplies conductivity; no inferred wall transport or additional heat
owner is introduced. Film coefficients must be finite and nonnegative; zero
is insulating. Both prescribed and ambient temperatures must be in the
declared run window, and native material-domain failures propagate.

The film surface solve brackets the interface between the cell and ambient
temperatures. Safeguarded secants, with every third iteration a bisection,
have a hard 192-iteration budget after at most two endpoint queries. The
wall-temperature fraction residual must be at most `32*f64::EPSILON`; failure
returns `ThermalWallConvergence` without publication. A stable series
conductance avoids large coefficient products and subtraction of nearly equal
interface temperatures in the insulating/isothermal limits. This numerical
surface tolerance is distinct from the optional coupled-flow endpoint controls.

One sparse surface pass evaluates each local face flux once per original FV
evaluation; divergence and the wall-heat account reuse that same value, with
the original deterministic reduction order. There is no additional full-cell
sweep or modeled SGS wall heat. Wall stress work remains a separate account,
including at a face that is simultaneously heated/cooled and moving. Local
molecular, Smagorinsky, WALE and Vreman flow paths retain their existing owners.
Channel dynamic SGS now accepts these thermal tables in either attachment
order. A sparse scan determines each wall plane's actual scalar derivative
policy. Uniform isothermal planes retain their own prescribed value, including
different values on the two sides; uniform insulated planes retain the original
zero physical gradient and separate Favre derivative. Full-plane overrides can
replace the global default. Partial overrides retain a uniform policy only
when they agree with the remaining default faces. Zero flux and zero film
coefficient are exact insulation.

Heterogeneous, nonzero-flux and finite-film planes instead use the shared
channel filter's `ScalarWall::Interior`: second-order three-fluid-center normal
derivatives at both physical and Favre-filtered levels. A spatially varying
prescribed value or native material-dependent wall flux does not authorize an
invented wall datum on the Favre quotient. No averaged wall temperature,
additional material query or SGS wall heat is supplied. The original molecular
faces still impose every actual local boundary law. This extends the existing
exploratory plane-averaged closure, retaining its derivative truncation and
plane-averaging assumptions; it does not calibrate a localized turbulent
thermal boundary layer. Empty tables keep the original global path exactly.

Changing the table at an accepted time preserves populations, temperature,
force, energy, clocks and cumulative ledgers, and refreshes the current native
FV history under the new wall data. Channel scalar transport is refitted from
the current accepted fields in the same transaction, including when the table
is cleared; failures preserve the previously accepted coefficients. The next original trapezoid therefore uses
the new condition at both time levels, without an instantaneous energy impulse.
`set_thermal_walls_cancellable` checks cancellation before/after admission and
after the complete flux evaluation; any refusal restores controls and history.
Characteristic interval copies retain the same sparse table. Focused tests in
`nsf3::thermal_wall::tests` exercise six orientations and an interior solid,
independent film series-resistance/interface algebra, heterogeneous heating
and cooling with moving-wall work, original/default equivalence, refused and
late-cancelled updates, characteristic trials and actual-pressure source refusal.
These conditions prescribe external thermal data; they do not introduce solid
heat storage, a turbulent wall function or a new natural-convection qualification.

`preview_thermal_wall_step` exposes one reversible original fluid step for
finite-wall composition; `preview_thermal_wall_interval` accepts a positive,
explicitly bounded count of those same native steps. Its `NsfThermalWallTrial` retains the exact grid,
complete native energy ledger, original prior energy and canonical per-face
`NsfThermalWallExchange` joules. Positive local heat enters the fluid on either
orientation. The starting boundary-history refresh returns its already
evaluated sparse heat fluxes; the accepted conservative corrector retains the
same endpoint face fluxes used by its energy equation. Their original trapezoid
supplies each exchange. Residual probes do not replace that consumed history,
and no accepted-temperature property requery reconstructs the heat afterward.
Only explicitly listed local faces appear in this account; unlisted global
walls and mechanical wall work retain their distinct native accounts.

The candidate refreshes starting history under its proposed held table once,
including an unchanged table. Later steps carry the exact prior corrector's
sparse flux vector together with its existing native divergence history;
there is no intermediate accepted-temperature property requery or history
restart. The grid is copied once for the complete interval. Per-face heat and
all native energy-transfer terms accumulate in step order, final storage comes
from the last step, and a multi-step residual closes the entire interval.
One-step calls preserve their original native ledger and residual arithmetic.
Ordinary public setters retain their unchanged-table no-op, and ordinary steps
retain no extra face history.

The trial retains the exact source-free physical state, supports existing held
openings, gravity, body forces and molecular/local/dynamic SGS, and still refuses
manufactured volume callbacks whose external state cannot be rolled back.
Periodic and channel dynamic closures share their already-admitted mutex
scratch through `Arc`, without allocating a second 200-byte/cell filter payload.
Each fit overwrites its required scratch from the current candidate density,
Guo velocity and temperature before consumption; these closures have no
Lagrangian or other temporal model history. The exclusive trial borrow prevents
concurrent access to the accepted physical owner. Accepted transport
coefficients, controls, histories and reports remain separately retained.
A refused/dropped trial may leave only private filter scratch overwritten, as
an ordinary refused native step can; retry performs a complete fresh fit.
The trial helper retains original filter geometry, controls and poisoned-lock
refusal without a second scratch allocation. The native wall update refits
transport for each proposed boundary; molecular heat alone crosses the wall.

`cell_count()` reports the original allocated grid including masked solids.
The preview's full-grid cell cap bounds candidate storage, and a separate
`maximum_cell_steps` gate admits the checked product of that count with the
requested positive step count before copying. Zero steps, work/clock overflow
and insufficient budgets refuse. Original corrector limits remain in force;
cancellation is polled between steps, at native work boundaries and while
accumulating faces. A late failure discards the entire candidate, including
any successfully computed prefix. `accept` publishes the already admitted
complete interval without additional fallible work.

This supplies the fluid side of `fs-couple::nsf_thermal`, whose original finite
thermal cells receive the opposite actual joules and own their storage and
spatial conduction. That adapter separately gates local and joint energy and
retains first-order held-endpoint thermal coupling. Native regression sources
cover all heat-flux orientations, original direct-step equivalence for newly
installed controls, unchanged-table continuation, optional endpoint residual
control, work/source refusal and late cancellation. Interval regressions compare
the unchanged direct native trajectory and all energy accounts, distinguish
continuous corrector history from per-step restarts, preserve physical/dynamic
coefficients after a dropped or late-cancelled active channel interval, exercise
shared workspace reuse and native heat-cap refusal, and replay periodic fits.
The same native topology, molecular wall heat and exploratory turbulence
qualification remain authoritative. This is a coupling seam, not an
independently calibrated turbulent or conjugate benchmark.

**Prescribed open faces (pexfp.138).**
`set_open_face(axis, high, Some(NsfOpenState))` declares the density [kg/m³],
Cartesian velocity [m/s] and temperature [K] of a uniform exterior ghost plane.
The face lies halfway between its ghost centre and the adjacent fluid centre.
Only one non-periodic axis with periodic transverse domain edges is admitted;
internal solid masks may describe a nozzle/channel but must not touch the open
plane. Solid/open junctions, intersecting open/wall domain edges, moving open
faces and manufactured volume sources are refused. Prescribed reservoir and
characteristic boundary modes are mutually exclusive in either attachment order.
`None` restores a stationary wall. Incoming populations use non-equilibrium
extrapolation (Guo, Zheng & Shi 2002, Chinese Physics 11, 366–374,
doi:10.1088/1009-1963/11/4/310), including the post-collision half-force shift.
The same exterior state enters the centred thermal-pressure face force and the
total-energy face flux (enthalpy, extrapolated stress work, and conduction).
Caloric energy and face transport use the grid's gas owner, including NASA9
when supplied. Changing a boundary restarts flux history and adjusts the
pre-collision half-force moment without changing the fluid's physical velocity
or energy. Refused updates leave the state and prior boundary intact.

`EnergyLedger::open_energy_j` reports the signed net energy entering prescribed
faces, separately from `wall_heat_j` and characteristic `open_boundary_energy_j`;
the conservation residual subtracts all three. With prescribed faces the legacy
`viscous_heating_j` remainder also contains kinetic energy transport and is
**not** a physical viscous-heating measurement.
The `nsf3::open::tests` checks cover uniform throughflow in both directions on
all three axes (G0, error < 1e-11, sealed-wall negative control), heat input
and the total-energy balance (G0, relative residual < 1e-12), boundary-update
refusals and velocity preservation, and worker-count equality (G5).
The existing thermochemical suite also exercises open throughflow away from the
caloric reference temperature, arbitrary formation-energy offsets and refusal
outside source pressure support. A fixed-domain pulse at 64/128/256 cells reports
the prescribed boundary's reflection; it does not assert a transparency band.
The native 20 mm-domain, sigma=1 mm Gaussian diagnostic measured amplitude
ratios 0.68695 / 0.81251 / 0.90254 at dx=0.3125 / 0.15625 / 0.078125 mm
(-3.26 / -1.80 / -0.89 dB): reflection increases with refinement. This broadband
packet result is not frequency-resolved R(f) over an instrument's declared band.
These are boundary-implementation checks, not a spatial convergence ladder,
characteristic zone-to-DG coupling, < -40 dB reflection qualification, jet
validation, or instrument/audio evidence. A prescribed exterior state generally
reflects acoustic waves. The separate characteristic mode has its own contract
and qualification limits below; it is not replaced by this primitive.

**Prescribed open faces (pexfp.138).**
`set_open_face(axis, high, Some(NsfOpenState))` declares the density [kg/m³],
Cartesian velocity [m/s] and temperature [K] of a uniform exterior ghost plane.
The face lies halfway between its ghost centre and the adjacent fluid centre.
Only one non-periodic axis with periodic transverse domain edges is admitted;
internal solid masks may describe a nozzle/channel but must not touch the open
plane. Solid/open junctions, intersecting open/wall domain edges, moving open
faces and manufactured volume sources are refused. Prescribed reservoir and
characteristic boundary modes are mutually exclusive in either attachment order.
`None` restores a stationary wall. Incoming populations use non-equilibrium
extrapolation (Guo, Zheng & Shi 2002, Chinese Physics 11, 366–374,
doi:10.1088/1009-1963/11/4/310), including the post-collision half-force shift.
The same exterior state enters the centred thermal-pressure face force and the
total-energy face flux (enthalpy, extrapolated stress work, and conduction).
Caloric energy and face transport use the grid's gas owner, including NASA9
when supplied. Changing a boundary restarts flux history and adjusts the
pre-collision half-force moment without changing the fluid's physical velocity
or energy. Refused updates leave the state and prior boundary intact.

`EnergyLedger::open_energy_j` reports the signed net energy entering prescribed
faces, separately from `wall_heat_j` and characteristic `open_boundary_energy_j`;
the conservation residual subtracts all three. With prescribed faces the legacy
`viscous_heating_j` remainder also contains kinetic energy transport and is
**not** a physical viscous-heating measurement.
The `nsf3::open::tests` checks cover uniform throughflow in both directions on
all three axes (G0, error < 1e-11, sealed-wall negative control), heat input
and the total-energy balance (G0, relative residual < 1e-12), boundary-update
refusals and velocity preservation, and worker-count equality (G5).
The existing thermochemical suite also exercises open throughflow away from the
caloric reference temperature, arbitrary formation-energy offsets and refusal
outside source pressure support. A fixed-domain pulse at 64/128/256 cells reports
the prescribed boundary's reflection; it does not assert a transparency band.
The native 20 mm-domain, sigma=1 mm Gaussian diagnostic measured amplitude
ratios 0.68695 / 0.81251 / 0.90254 at dx=0.3125 / 0.15625 / 0.078125 mm
(-3.26 / -1.80 / -0.89 dB): reflection increases with refinement. This broadband
packet result is not frequency-resolved R(f) over an instrument's declared band.
These are boundary-implementation checks, not a spatial convergence ladder,
characteristic zone-to-DG coupling, < -40 dB reflection qualification, jet
validation, or instrument/audio evidence. A prescribed exterior state generally
reflects acoustic waves. The separate characteristic mode has its own contract
and qualification limits below; it is not replaced by this primitive.

**Manufactured-solution hooks.** `NsfGrid3::new_with_sources(spec, solid, initial, Some(source))` takes a `SourceFn`: for a cell
index and a lattice time (steps) it returns `(S_rho, S_m, S_E)` in lattice units for `d rho/dt + div(rho u) = S_rho`,
`d(rho u)/dt + div(rho u u + p I - sigma) = S_m` and `dE/dt + div((E + p) u - sigma . u - kappa grad theta) = S_E`. The momentum
source joins the force at both time levels (so the Guo half-force velocity and the stress shift see it), the energy source enters
the trapezoid with the flux divergences, and the mass source is added to the populations with the lattice weights (so it carries
no momentum) HALF BEFORE AND HALF AFTER the stream: all of it before the stream leaves the stream's pressure gradient one source
increment ahead of the momentum source and a first-order momentum error of `theta0 grad(S) / 2` per step (seen as a velocity
error of that size after the first step and a density error that did not shrink with refinement). With the half/half split the
residual per step is `theta0 laplacian(S) / 4`, which is second order globally under acoustic scaling. A grid without sources
is bitwise unchanged. The hooks exist for the qualification (`mms_full_nsf_order`) and are not an instrument source API.

**Ledger.** `EnergyLedger` reports internal, kinetic, pressure-work, viscous-heating, wall-heat, wall-work, boundary-energy and
volume-source terms in joules. `source_energy_j` is the exact trapezoid of the same `S_E` samples applied to fluid cells;
solid cells and callback reevaluation cannot add an audit-only supply. Negative supplies remove energy. The total-energy
residual subtracts wall, boundary and volume supplies independently; no internal-exchange remainder is a balancing credit.
Cumulative source energy is published only with an accepted step and uses the medium's fixed sensible-energy chart.

`pressure_work_j = -integral(p div u) dV dt` is the reversible exchange (trapezoid in time, fourth-order divergence).
For source-free and heat-only steps, `viscous_heating_j` is DEFINED as
`internal change - wall heat - pressure work - source energy`, with the existing boundary terms where applicable.
With no source and no moving wall, it equals the kinetic-energy loss net of pressure work term by term to round-off.
THAT IDENTITY IS BOOKKEEPING. The physical content (the heating is continuum dissipation) is the Taylor-Green decay rate.

Manufactured mass/momentum sources also alter kinetic storage, and this owner has no separately qualified local
source-work split for its streamed mass source. Such steps report the signed residual in
`unresolved_internal_exchange_j` and zero in `viscous_heating_j`; the unresolved term includes both dissipation and
unseparated source mechanics and must not be deposited as heat. Uniform acceleration can give a negative unresolved
term while the exact supplied energy becomes kinetic energy. This does not weaken the independent total-energy audit.
This volume-source distinction does not change the fluid dynamics.

Focused regressions in `tests/nsf3_thermochemical.rs` cover affine signed heating/cooling in both gas branches,
solid-cell exclusion, one/three-worker replay, uniform source-driven acceleration, mass/energy supply, and existing
source-domain rollback. They are authored behavioral checks, not a native execution or empirical qualification claim.

**Conformance (`tests/nsf3_battery.rs`).** `equilibrium_recovers_conserved_moments` (G0: mass, momentum, the second moment
and every third moment the lattice can represent, 1e-13), `refusals_cap_and_cap_plus_one` (G0, including the viscosity envelope, the anonymous ratio and the Mach-bound refusal; the
spec declares `MAX_QUALIFIED_MACH`),
`closed_box_energy_budget_closes` (G0: random velocity + temperature state, six bounce-back walls, 1e4 steps at the
viscosity of dx 50 um, total energy <= 1e-11 relative), `viscous_heating_equals_kinetic_loss` (G0 identity <= 1e-13 of
the total energy per step; G1: kinetic-energy decay rate within 5e-3 of `4 nu k^2` and accumulated heating within 2e-2,
shrinking with refinement), `isothermal_wall_heat_ledger` (G0: wall heat equals the energy change to 1e-11),
`adiabatic_sound_speed` (G1: 48 / 96 / 192 / 288 points per wavelength, relative error <= 1e-4 at the finest rung,
second order observed >= 1.8, and the isothermal speed fails; the 32-point value is printed, not asserted),
`nsf3_bitwise_worker_count_invariance` (G5). Each physics test prints one execution receipt (JSON line: lane, stage, verdict,
formulation, conservation class, gas receipt with blake3, rungs `{dx_m, dt_s}`, ledger, sound-speed error, observed values,
tolerance, elapsed); `scripts/e2e/music/lbm_qualification.sh --lane nsf-energy --out DIR` runs the battery natively and requires
one passing receipt per stage (closed-box, viscous-heating, wall-heat, sound-speed). Native evidence of the first landing:
`wa.sh test 295c9c4d5 --release -p fs-lbm --test nsf3_battery` on trj, 7 passed, exit 0.

**No-claim boundaries (`nsf3`).** The sound-speed band of 1e-4 is met at 288 points per wavelength, NOT at 32:
dispersion is second order (observed order 2.00 between 192 and 288), and the athermal D3Q27 scheme alone is already
about -1.8e-4 at 32 points. Not turbulence-resolved, not shock-capturing, no open boundary and no moving boundary other than a
tangential velocity of a non-periodic domain edge, no entropic or
energy-population thermal LBM (the thermal pressure is a force), no temperature-dependent equilibrium. The force-coupled
lattice resolves the sound speed `sqrt(gamma theta)` only while it stays below the lattice speed: the Mach and
temperature-ratio bounds are refusals, not tuning. The bulk-to-shear ratio is a cited input, not a module result; its
uncertainty (0.67 to 0.84) is not propagated. No recording, graded or instrument claim.

### Exploratory NSF3 subgrid momentum and heat transport

`NsfGrid3::with_sgs(NsfSgsConfig)` explicitly enables the shared full-gradient
Smagorinsky, WALE or Vreman closure on the existing low-Mach NSF owner. Only the
deviatoric relaxation receives `nu_t`; the source-owned molecular bulk rate is
unchanged. The matching stress enters the conservative total-energy flux.
Interior heat conduction additionally receives `rho (cp/R) nu_t / Pr_t` in
lattice units, with explicit positive turbulent Prandtl number and current-cell
NASA9 heat capacity when a thermochemical source is attached. The same face
coefficient is shared by both cells. No second viscous-heating source is added.

Transport is refreshed from the force-corrected velocity in the existing
predictor and accepted endpoint. Periodic derivatives are centred; solid and
domain-wall derivatives use actual halfway wall values, including prescribed
wall motion. Impermeable wall heat transfer retains its molecular half-cell
resistance, with zero modeled wall heat flux. Characteristic faces use their
original reconstructed half-cell velocity, while a localized transverse
reservoir aperture or complete prescribed face uses its original ghost velocity
one full spacing away. A one-cell layer between a prescribed ghost and an
internal solid uses the quadratic derivative with distances one and one-half,
respectively; it is not interpreted as two stationary halfway walls.
The same predictor faces feed SGS, thermal-pressure forcing and energy flux.
Adjacent algebraic conductivity is extrapolated to either opening (zero
normal gradient of the modeled transport coefficient); temperature differences
retain their original half-cell or full-cell distances. Modeled conduction
and stress work enter the existing total-energy and open-boundary accounts,
including selected characteristic patches. Rebinding controls refreshes SGS
and conservative history atomically; rejected/dropped trials retain both.
Complete prescribed-reservoir faces admit active SGS in either configuration
order; `set_open_face` refreshes the closure, matching shear stress and heat
flux together while retaining the original physical state. Removing a face
restores its actual stationary wall stencil. The existing
Mach, temperature and total-viscosity envelope remains in force; ratio-cap
violations refuse rather than clip. Zero coefficient preserves molecular
arithmetic, and acceptance-render use remains forbidden.

`sgs_transport` exposes accepted-cell eddy viscosity and conductivity. Focused
tests cover moving-wall gradients, independent bulk relaxation, conservative
heat diffusion, closed energy, worker identity, exact zero control, fixed-Re
closure refinement, NASA9 heat capacity, uniform throughflow without spurious
open-face wall shear, reservoir gradients, open heat flux, complete
mass/energy exchange and transactional admission failures. The prescribed-face
tests cover all three axes and closures, exact quadratic ghost/wall gradients
on both orientations, heated throughflow with conservative total energy,
1/3-worker identity, zero controls and viscosity-cap rebind rollback. These are
algebraic low-Mach transport integrations, not a new compressible-turbulence
qualification: no SGS kinetic-energy equation, isotropic SGS pressure,
pressure-dilatation, shocklet closure or thermal wall function is supplied.
The modeled shear and Prandtl heat-flux form follows NASA report 20200002378,
equations (8)-(10), without its additional compressible terms.

`NsfGrid3::with_dynamic_thermal_sgs(NsfDynamicThermalSgsConfig, maximum_bytes)`
instead fits momentum and temperature diffusivity independently through the
shared Favre--Germano producer. This path requires a fully periodic, all-fluid
domain with at least three cells per axis. Each original predictor and accepted
endpoint uses its current density, Guo-corrected physical velocity, temperature
and material response. Fitted viscosity modifies only deviatoric relaxation;
`rho (cp/R) alpha_t` enters the same conservative heat flux, including cases
with zero modeled stress and nonzero modeled heat. No second heating source or
energy integrator is introduced.

Momentum-ratio and total-viscosity limits retain the original molecular source.
Heat-ratio admission checks both the cached cell source and the original shared
face conductivity, without a second face material query. Positive modeled heat
with zero molecular conductivity refuses a finite ratio bound. Local and
dynamic controls are mutually exclusive; wall, solid and open-face topologies
are outside this dynamic path. Refusal preserves accepted transport, fields
and energy accounts. The two reused filter buffers require 200 nominal bytes
per cell, with actual-capacity admission; the existing NSF state, predictor
and transport vectors are outside that filter budget. A single serial fit lock
preserves the original parallel cell/face loops. This does not make the NSF
step allocation-free or add whole-step cancellation. The existing low-Mach,
temperature and stability qualifications remain in force; no DNS, wall-model,
shock or complete compressible SGS-energy claim is added.

`NsfGrid3::with_channel_dynamic_thermal_sgs(NsfChannelDynamicThermalSgsConfig,
maximum_bytes)` extends that same native transport/energy owner to an entirely
fluid plane channel: exactly one halfway-wall axis, two periodic axes and at
least three cells on every axis. Its paired channel momentum/scalar workspace
also requires 200 nominal bytes per cell. Current stationary or tangentially
moving wall velocities and the actual thermal table enter each fit. Uniform
isothermal or adiabatic planes retain their declared values/derivatives;
heterogeneous or native flux/film planes use the original three-fluid-center
reconstruction at both scalar levels, as specified above. Adiabatic physical
scalar derivatives and the Favre-filtered scalar's one-sided normal derivatives
remain distinct, as specified by the shared channel model.

The original material heat capacity, Guo velocity, independent molecular bulk
rate, interior shared-face conductivity and conservative total-energy owner
remain authoritative. Modeled wall heat flux remains zero under the original
molecular half-cell conduction law. Moving-wall traction/work instead retains
the original adjacent-cell total-stress approximation, which is first order
for developing flow; this does not explicitly impose zero SGS viscosity at
the geometric wall or establish its near-wall asymptotic accuracy. Wall-speed
changes refresh transport and energy history transactionally. Attaching open
faces or characteristic boundaries after either dynamic closure refuses before
mutation, retaining the admitted filter topology.

The focused native channel regressions cover current source heat capacity and
force/wall data, heat-only planes, molecular wall conduction, shared-face
admission, moving-wall energy, worker identity, workspace reuse, refused-step
rollback, zero controls and topology/budget refusals. The periodic regressions
retain their original physical assertions. Additional local-wall regressions
compare both attachment orders and exact uniform trajectories, preserve
distinct uniform side values, exercise heterogeneous temperature/flux/film
planes with one molecular heat account, verify current-endpoint refits and
coefficient/cancellation rollback, and cover all three wall axes and both sides.
The shared channel derivative MMS remains the independent physical/Favre
gradient oracle. Interior solids, wall functions, spatially varying wall
velocity and channel DNS calibration remain outside this path.

### Evolving Mann turbulence on an original reservoir aperture

With `turbulent-inflow`, `NsfEvolvingMannInflow` consumes the original full-tensor
`FixedStepMannField` at exactly the NSF grid's physical time step. Construction
binds the current CFD step to the current atmosphere step, so a retained
snapshot can initialize a running grid without inventing past evolution.
Later updates use checked relative step increments. The atmosphere owns phase
advection once; the adapter supplies the explicit physical frame, mean velocity
and original gas-source pressure/temperature at actual exterior ghost centres.
Point, face-average and cell-average observations reuse the existing Mann
spatial filters without restoring unresolved variance.

Each apply advances a trial temporal state and fills a reused canonical tile
buffer. The original atomic reservoir-profile rebind admits gas states and
refreshes its Guo force and conservative flux history; only success publishes
the atmosphere. Refusal or cancellation preserves accepted boundary and source
state. Reapplying at the same CFD time consumes no noise. Native fluid stepping
remains a separate operation, as with `NsfFrozenMannInflow`; a later refused
fluid step retains the already accepted boundary control. No full-grid clone
or additional fluid/energy solver is introduced.

The caller budgets both temporal owners, retained profile capacity, optional
filter gains, mode samples and OU transitions. Spectral geometry is immutable;
the original boundary rebind and NSF step retain their existing allocations.
Grid spacing, dimensions, dt, aperture orientation and tile membership stay
bound. Spatial and advective phase bounds retain their original refusals;
the temporal bound covers advection only because OU innovations are not
bandlimited. Reports expose both clocks, authored correlation times and actual
work. The default observation is a synthetic homogeneous turbulent inlet
with held boundary controls during each native step. Authored OU evolution
supplies no resolved eddy-turnover, nonreflecting-boundary or temporal-order
qualification. The explicit planar observation below adds its own boundary
condition.

### Plane-blocked Mann observations on frozen and evolving inlets

Both `NsfFrozenMannInflow` and `NsfEvolvingMannInflow` expose
`with_plane_wall(MannPlaneWall)`. The original atmosphere supplies a
harmonic no-penetration correction before the existing frame rotation,
physical-unit conversion, gas admission and atomic profile publication.
The plane belongs to the Mann frame; every actual exterior ghost observation
footprint must lie in its upper half-space. Mean velocity must be tangential, and the frozen
adapter also requires tangential convection because it translates query
coordinates. Tangency uses the admitted rotation's roundoff scale without
altering the supplied velocities.

The wall observation reads the same frozen snapshot or evolving trial
quadratures. Point sampling adds no numerical arrays. Face and cell averages
prepare fs-atmo's complete `MannWallBoxFilter`, retaining both the homogeneous
wave and the evanescent correction over the original rotated ghost footprint.
The original native gas/force/energy owner consumes this resolved mean velocity;
nonlinear momentum/energy fluxes are not averaged and missing variance is not
restored. `wall_spatial_filter()` exposes the prepared observation for the
source's two-position filtered-covariance query. `spatial_filter()` continues
to expose only homogeneous observations.

Both attachment orders use one selected filter. Adding a plane after an
average rebuilds it within the retained byte allowance; old gains are released
before replacement. `with_plane_wall_cancellable` exposes cancellation during
that preparation. The evolving owner also admits the actual total numerical
capacity against its existing workspace limit. Per-query samples allocate no
filter arrays or random draws. The complete footprint, not just its centre,
must remain above the plane: partially submerged support refuses instead of
being clipped or renormalized. Rejected samples and cancellation preserve the
accepted profile and atmosphere; retries use the original temporal owner.

This remains outer-flow slip-plane conditioning, without a viscous no-slip
layer, wall-generated turbulence or calibrated near-ground spectra. Regression
coverage includes independent physical-space quadrature of actual rotated
ghost-cell/face averages, both attachment orders, finite capacity, footprint
refusals, native energy exchange, idempotent clocks, and evolving retry after
cancellation or boundary refusal. Existing point-sampling regressions retain
their wall-normal and source-history checks.

### Source-owned NASA9 energy and local transport in the NSF grid

`NsfSpec<G = GasSpec>` admits `ThermochemicalNsfGas` through the same
`NsfGrid3<L>` constructors. Existing GasSpec call sites retain their constant-cv
branch and arithmetic order. The source owner retains the NASA9 mixture, original
transport cards/selections, typed reference query and bounded `GasEnergySolve`.
Only T/p evolve; composition and every other material coordinate stay fixed.

Each fluid-cell predictor resolves pressure as rho*R*T using that mixture's
constant R and actual conserved density. The resulting source viscosity drives
the existing collision and stress rates. Shared energy faces resolve source
conductivity at mean T and mean physical pressure; isothermal wall faces use
mean fluid/wall T at the adjacent fluid pressure. No source transport is clamped
or replaced with the reference-pressure diagnostic table. Every source query,
inversion, local Mach/gamma and kinematic-viscosity check can refuse before
publication. Accepted cell responses, populations, energy, flux history and
ledgers publish together. The original fixed predictor-corrector passes remain;
this is not a new compressible solver.

The conserved internal-energy chart is u_s(T) = u_NASA9(T) - u_offset, with the
single retained offset u_offset = u_NASA9(T_ref) - cv(T_ref)*T_ref. It matches
cv*T at the reference state and removes arbitrary formation-energy offsets
from the separate lattice mass / finite-volume energy fluxes. NASA9 absolute
energy may be negative. Lattice specific energy converts through (dx/dt)^2.
Cell and face energies evaluate this chart as cv(T_ref)*T_ref plus the original
NASA9 integral of Cv from T_ref to T. Density/energy reconstruction inverts that
same temperature-anchored increment. Common formation references therefore cancel
before floating-point evaluation, including when they dwarf the thermal change;
no subtraction of large absolute NASA9 energies participates in evolution.
The existing material caloric branch is prepared once per native timestep and
reused across every cell and corrector, while endpoint transport is still resolved
at the actual state. No new polynomial evaluator or caloric integrator participates.
Construction admits this same source branch over the entire inversion bracket,
including positive Cv. A raw NASA9 bracket must stay within one region per species;
crossing a region boundary requires the source owner's explicit
`with_caloric_continuity` policy. The NSF owner neither reconciles source offsets
implicitly nor waits until the first timestep to refuse an unsupported bracket.
`NSF3_BIT_SEMANTICS_VERSION=6` records this source-path arithmetic change;
the legacy GasSpec branch is unchanged. `energy_reference_j_kg()` exposes the offset.
Reported energies and manufactured `SourceFn` energy supplies use this chart;
an absolute-energy source carrying mass must subtract u_offset times its mass
supply. Composition-changing/reaction sources are unsupported.

The source reference T/p must equal the grid reference, and its inversion
bracket must equal the declared run window. Source properties are re-evaluated
at actual local states, including final energy inversion, so a valid reference
point does not authorize an unsupported density/temperature trajectory.
`thermochemical_state(cell)` exposes the original source receipt at the
accepted physical state; `medium()` retains the source owner.
`dispersion::GasPoint::from_law` uses the exact reference coefficients for
a LOCAL uniform-state small-signal linearization. Its gamma is not a global
nonlinear caloric approximation.

The existing Mach/viscosity bounds remain necessary admission guards; the
legacy perfect-gas convergence and stability evidence does not establish
nonlinear NASA9 qualification. No reactions, species transport, molecular
relaxation, phase change or variable bulk/shear ratio are supplied. Synthetic
regressions in `tests/nsf3_thermochemical.rs` exercise caloric heating and local
EOS/transport feedback, exact formation-energy-offset invariance during density
changes and signed small-heat evolution, data-only conductivity substitution, worker invariance, local-source
rollback, temperature-dependent Mach admission and legacy ingress equivalence.
Their presence does not claim native execution or empirical gas validation.


## No-claim boundaries

- Tiled D3Q19 grids default to BGK + Guo on a dense set of aligned SoA
  tiles. The selectable central-moment and `ReducedCumulant` cell operators
  support body force with deterministic `O(Q^2)` forward transforms and an
  `O(Q)` algebraic inverse. They remain scalar references with no new high-Re
  stability claim. A bounded complete-cell microbenchmark on Xeon Platinum
  8573C, rustc 1.100-nightly (908501772), native opt-level 2, measured
  2.50–2.93 times lower collision time than the version-2 dense solve across
  forced/unforced central and reduced-cumulant modes (medians of three
  30,000-call samples). These shared-host measurements cover cell collision;
  they are not stream/collide GLUP/s, worker-scaling or cross-ISA evidence.
  D3Q27 has the separate
  checked full collision/reference-grid surface above, pending native qualification.
  Normalized drag/lift histories and bandwidth roofline /
  fs-tilelang kernels remain staged. The sparse active-tile layer
  (`d3q19::sparse`) claims determinism and active-set memory proportionality
  only — NOT throughput (the GLUP/s perf lane is bead 712t): its collision is
  the scalar per-cell kernel, occupancy is whole-tile (per-cell masks and the
  free-surface activation/deactivation rules are WS1-E's increment), tile
  size is pinned at 4³ (the 4³→8³ autotune sweep is follow-on work through
  the fs-exec tune path). The frozen serial sparse golden covers only its
  pinned L-shaped active set and step count; a separate bitwise test covers
  the worker-count sweep. The retained transactional stream destination raises
  population storage from two to three buffers per active tile; this is the
  explicit memory cost of failure-atomic mid-stream cancellation, not a
  throughput optimization. The raw D2Q9 stationary-wall momentum receipt is
  not a normalized aerodynamic
  coefficient. Geier
  et al.'s primary derivation (doi:10.1016/j.camwa.2015.05.001) explicitly
  restricts itself to D3Q27 after identifying non-refining D3Q19 anisotropy.
  `ReducedCumulant` therefore implements only the general cumulant definitions
  for the independent D3Q19 `C220`/`C202`/`C022` projection; it is not a
  "Geier D3Q19" operator and cannot borrow the paper's high-Re evidence.
- The current SIMD slice accelerates the all-fluid Duct axial-BGK collision and
  its pure-move pull stream. `BoundaryGrid3` mixed masks/general forcing and the
  central-moment/cumulant paths remain scalar; collision and streaming are not
  fused. Therefore this slice alone makes no stream+collide GLUP/s claim. Both
  native ISA batteries, the unchanged Duct golden, and the four-quadrant
  performance lane remain required before promotion evidence is green.
- SDF voxelization is midpoint classification followed by stair-step halfway
  bounce-back. It is second-order for the represented flat, lattice-aligned
  halfway wall, not a second-order certificate against the original continuous
  curved SDF. D3Q19 interpolation remains staged; the separate D2Q9 linear BFL
  path consumes caller-supplied link fractions and does not derive intersections
  from an SDF or certify geometric/convergence order.
- D2Q9 `MomentumExchange2` is a raw stationary-wall lattice impulse over the
  caller's exact cell mask. `MovingWallMomentumExchange2` adds a static-topology
  boundary-relative force evaluation for both halfway and caller-described
  linear-BFL links, but it does not move cell topology or initialize newly
  uncovered fluid. Neither
  receipt applies physical-unit conversion, reference-area normalization,
  blockage correction, averaging, or shedding-frequency estimation; the BFL
  path additionally has no quadratic/multireflection option or measured curved-
  geometry convergence deck. Therefore this is not yet the Re=100 cylinder
  Cd/St validation. Boundary-relative exchange improves the force receipt's
  frame behavior. The moving-Couette conformance deck now measures relative-
  profile and shear-impulse residuals for one flat, periodic-x, steady, low-Mach
  common-boost pair; it does not prove Galilean invariance for curved walls,
  topology changes, transients, finite-Re bodies, or the full coupled solver.
- `PartialSaturationCell2` consumes a caller-supplied cell volume fraction and
  one aggregate solid velocity; it does not rasterize geometry, resolve
  overlapping bodies, or attribute the aggregate receipt among bodies. This
  first rung is dense BGK D2Q9 with an additive Guo-force extension and selects
  the published non-equilibrium-bounce-back variant. It has no stresslet-accuracy,
  convergence-order, free-surface coexistence or high-Re stability claim.
  Its prescribed-velocity entry does not integrate rigid-body inertia; the
  finite-inertia coupling below adds a bounded velocity update, while pose
  evolution, geometry rasterization and complete moving-body FSI remain open.
- `UnderResolvedCylinderPair2` is a leading normal squeeze-film correction for
  smooth convex circular or osculating-cylinder geometry, Newtonian fully
  wetted creeping-gap flow, and raw lattice quantities per unit out-of-plane
  depth. Its caller supplies geometry, the gap and cutoff calibration, and
  enumerates each unordered pair exactly once. It neither estimates the
  already-resolved lattice resistance nor mutates fluid populations, so using
  an uncalibrated cutoff or applying the pair twice will double count force.
  The positive gap floor is a disclosed numerical regularization, not a
  roughness, contact, or continuum-validity model. There is no tangential or
  rotational lubrication, many-body hydrodynamics, cavitation, asperity/contact,
  elastic-hydrodynamic, thermal/non-Newtonian, gap-generation, rigid-body
  integration, cutoff-calibration, convergence, or validated-bearing claim.
  Full Reynolds/Elrod-Adams/EHL films remain owned by `fs-tribo`.
- D2Q9 `transition_wall_topology` consumes an already-discretized next wall
  mask; it does not integrate geometry motion, infer covered cells, or couple a
  rigid-body state. Equal-weight one-ring population averaging is a
  deterministic first fresh-cell rung, not a conservative remap or an
  accuracy-order certificate. The receipt deliberately exposes any net active
  mass/momentum change; callers must ledger or correct it rather than assuming
  conservation. Gas/interface transitions and newly covered-body impulse/work
  coupling remain staged.
- The separate `lbm-109` release fixture encodes the intended normalization,
  warm-up, detrended FFT, raw/split-window guards, primary-source envelopes,
  and empirical two-width Cd sensitivity treatment. Maskell's closed-tunnel
  law is not claimed to transfer to a periodic-y fixture; the linear-beta Cd
  intercept merely eliminates one assumed leading coefficient from two
  disclosed widths, and its magnitude is bounded. Strouhal uses the 16D result
  plus 12D sensitivity rather than a frequency-bin extrapolation. Until that
  explicitly ignored release lane is executed green on the combined batch
  snapshot, the crate still makes no completed Re=100 cylinder validation
  claim.
- `VelocityPressureX2` and `PressureX2` are low-Mach regularized fixture boundaries with periodic
  lateral closure. It is not a characteristic or non-reflecting far-field
  condition and makes no unbounded-domain or blockage
  claim. The cylinder battery must disclose its domain, resolution, warm-up,
  averaging window, and lateral-domain sensitivity separately.
- The cylinder fixture's diameter-10 stair step and `tau = 0.53` have no grid-
  convergence certificate, curved-wall accuracy claim, or seed-independent
  shedding-phase claim. Split-window and two-width agreement detect selected
  transients and domain sensitivity; they do not replace a resolution family.
- The face-generic regularized closure was selected instead of six
  face-specialized Zou-He tables because its Hermite stress projection has an
  independent second-moment oracle and preserves arbitrary tangential target
  components under one implementation. This is not evidence that a correct
  Zou-He/Hecht-Harting implementation is unstable. The current closure is a
  low-Mach fixture, not a high-Re characteristic/non-reflecting boundary.
  Prescribed velocity traces and simultaneous Guo force use that same closure.
- Open/wall intersections use an explicit per-link mixed policy: wall-crossing
  links bounce stationary, pure open-normal links are regularized, and links
  crossing both faces are wall-owned. This is not a geometric corner
  reconstruction, and exact target moments are claimed only on face-interior
  cells. Moving-wall cavity evidence is qualitative until the separate G2
  benchmark bead lands.
- Solid topology is static after initialization. Moving/deforming SDF
  boundaries require a future explicit population initialization, mass delta,
  and topology-transition receipt rather than reusing `voxelize_sdf`.
- Interface and gas flags exist so the plain core can share the future data
  model, but free-surface mass/VOF bookkeeping is not implemented in
  `Grid::step`; gas-side pulls currently bounce rather than reconstructing
  missing free-surface populations.
- Thermal and rheology fixtures are dense-grid correctness scaffolding, not
  validated LES, cumulant, sparse-tile, or production multiphase solvers.
- The free-surface implementation is a dense prototype with ledger and
  metamorphic gates. It does not yet claim quantitative dam-break agreement,
  contact-angle calibration, production wetting physics, or validated
  surface-tension breakup rates.
- The dynamic wetting adapter proves system- and endpoint-provenance-bound
  hysteresis admission and directional selection only. Its caller-supplied
  pinning deadband is recorded but has no independent calibration evidence or
  MatDB usage receipt. The adapter does not impose a returned angle in
  `phi_at`, evolve roughness or contamination, implement Cox-Voinov,
  Hoffman-Tanner, or another contact-line law, or claim capillary-rise or
  contact-line validation. The existing binary neutral/wetting solver path is
  behaviorally unchanged.
- The scaling assistant covers the `τ`/`ν`/`Mach` core; consuming fs-regime's
  dimensionless groups and emitting a full `dx`/`dt` unit conversion with
  Evidence provenance is the fuller deliverable.
- Grid refinement is the TWO-LEVEL 1:2 channel coupling with Dupuis–Chopard
  non-equilibrium rescaling and a FIRST-ORDER temporal interface handoff.
  Actual postcollision populations crossing the ownership plane now enter
  a reused directional flux register. Both fine substeps carry their original
  one-quarter cell-volume weight; the resulting mismatch changes only the
  owned coarse interface row before restriction. Interpolation to the two
  adjacent coarse centres preserves the line of action of every directional
  transfer, closing mass, linear momentum and interface torque without a
  global velocity correction. A scalar residence-time register also corrects
  coarse acceleration for mass crossing at the first and second fine substeps'
  quarter- and three-quarter-interval times. Its local Guo-frame update preserves
  physical nonequilibrium stress; fine collision impulses and prescribed
  external force keep their native owners. Wall loads and angular transport
  through the periodic cut remain external exchanges. Compact-pulse/origin/
  locality, variable-density acceleration and independent force-plus-wall
  regressions cover these accounts. A focused actual-kernel
  execution of the unchanged lbm-108 test gives about 0.33% peak Poiseuille
  error and 0.26% shear-decay-rate discrepancy at its original resolution;
  its acceptance bands are unchanged. No second-order, energy-conservation
  or arbitrary-flow stability claim follows. Post-collision second-order
  transfer, general octree topologies and DWR-driven adaptation remain
  separate successors.
- Legacy dense contact-line fixtures remain MODEL-BRACKETED (neutral vs wetting
  fill ghosts; lbm-107 reports the sensitivity band). The opt-in dense 2-D and
  sparse 3-D planar equilibrium-angle closures retain their separate, bounded
  conventions and do not claim calibrated dynamic wetting.
- ADJOINT HONESTY (plan §8.7 [M]): free-surface LBM adjoints are NOT
  promised — cell-conversion events make the map non-differentiable;
  gradients for free-surface objectives go through surrogate or
  gradient-free lanes. The model card is this paragraph.
- Pouring scenarios: tilt schedules enter as the rotating gravity vector
  (lbm-106 pins 90-degree equivariance at 3e-14); full fs-scenario moving-
  frame integration and Plateau–Rayleigh breakup SCORING (beyond the
  qualitative fragment gate) are staged with the vessel flagship.

## WS1 3-D feature preservation: thermal + non-Newtonian (bead eg62)

`d3q19::coupled` ports the two VALIDATED 2-D capabilities to the 3-D
path without forking constitutive code. `PlatesGrid3` is a plain-storage
D3Q19 plate channel (rigid z-plates via halfway bounce-back, periodic
x/y) with per-cell relaxation time and per-cell external force, colliding
through the SHARED `collide_cell3` BGK + Guo kernel; the tiled SoA
`Duct`/`BoundaryGrid3` frozen BGK arithmetic remains unchanged.
`ThermalLbm3` adds a D3Q7 advection-diffusion temperature
population (c_s² = ¼) with anti-bounce-back fixed-temperature plates and
Boussinesq coupling into the Guo term. `update_tau3` closes tau and shear
rate from the force-corrected D3Q19 non-equilibrium second moment through the
SAME safeguarded local constitutive solver, `rheology::Rheology` law,
`TAU_FLOOR`/`TAU_CAP`, and `RheologyStats` ledger as 2-D. Its previous tau is
only a warm start, and the flow time integration remains explicit. The
plate-channel profile is certified against the SHARED
`powerlaw_poiseuille_analytic` (the geometry is 1-D in the wall-normal
axis). Battery (tests/d3q19_coupled.rs, verdict-JSON): Rayleigh–Bénard
onset bracket (decay at Ra = 1200, growth at Ra = 2500, Nu > 1), n ∈
{0.5, 1.0, 1.5} profiles to the 2-D 3% bar, 1e-11 mass ledger, exact
plate temperatures, linear conduction, ledger honesty, and a
replay-stable seeded determinism CANDIDATE.

### Corrected convection in both thermal slabs

`ThermalLbm` and `ThermalLbm3` share the D2Q5/D3Q7 scalar collision and
streaming operator. Linear equilibrium retains the zero-mass correction
`w_i*(1-1/(2*tau_g))*e_i dot delta(T*u)/c_T^2` from
[Chopard, Falcone and Latt (2009), equation (14)](https://doi.org/10.1140/epjst/e2009-01035-5).
The actual local thermal relaxation sets both BGK relaxation and its correction;
the diffusivities remain `(tau_g-1/2)/3` and `(tau_g-1/2)/4`, respectively.
Temperature collision samples the input physical velocity; force-frame
rebinding preserves that velocity, and flow advances afterward under the
selected held buoyancy/external load.

Each retained state stores one previous advective-flux vector per cell.
Transactional source snapshots copy that history. Fresh history, an explicit
`reset_advection_history`, or a nonzero onset perturbation uses one startup
step with zero correction. Direct population or discontinuous velocity
replacement requires that explicit restart.

Prescribed heat uses equilibrium half-kicks at the current physical velocity,
with the same `delta_T*u` shift in history. This preserves kinetic
nonequilibrium and the transport derivative through each source substep;
heating is not mistaken for advective-flux evolution. Opposite correction
sources cancel locally, and wall exchange remains independently measured from
anti-bounce-back links.

Six focused G0/G1 tests execute the production operator: axial/oblique Fourier
decay and phase at several relaxation times, time-varying velocity,
source-maintained diffusion, source/history invariants, local-rate wall heat
balance and startup. The full slab regression additionally targets late
rejection followed by accepted continuation against an independently warmed
ordinary solver. Existing low-Mach and lattice-stability limits remain;
no unconditional stability or temporal-order guarantee for the complete
partitioned buoyancy/source solver is added.

### Thermal wall laws in the original thermal slabs

`ThermalLbm::set_thermal_walls` and `ThermalLbm3::set_thermal_walls`
select independent, spatially uniform laws at the existing stationary,
impermeable halfway plates: prescribed temperature, insulation, prescribed
inward heat flux, or a film to a prescribed ambient. With no override, the
original public `t_bottom`/`t_top` temperatures retain their original
anti-bounce-back arithmetic. Resetting the override restores that behavior.

`set_spatial_thermal_walls` supplies one law per surface link on both plates,
in x order for 2-D and x-fastest `x + nx*y` order for 3-D. Heaters, cooling
patches, films and insulation may share a plate. Both complete arrays are
admitted before publication; ordinary streaming borrows the retained fields
without allocation or copies. `thermal_wall_fields` exposes the active laws,
while the existing `thermal_walls` getter retains the uniform defaults.
Setting uniform laws or resetting walls clears spatial overrides. For a
nonuniform field, Rayleigh/Nusselt diagnostics use the last measured mean wall
contrast and are unavailable before its first step; the existing steady,
source-free interpretation of the volume Nusselt diagnostic still applies.

For each reflected normal pair, `q = g_in - g_out` is the actual scalar
exchange and `T_wall = (g_in + g_out)/(2w)` is the surface trace. The film
closure solves `q = H*(T_ambient - T_wall)` locally, with bounded arithmetic
in both insulating and prescribed-temperature limits. Flux units are scalar
increment per link per step; `H = h*dt/(rho*Cp*dx)`. Collision already owns
diffusivity, so no second relaxation-time factor multiplies this exchange.
The reports sum actual exchange with compensated reductions and retain the
mean and local extrema of the surface trace. Ghost planes remain metadata,
not extra fluid heat storage. Reports describe heat already transported;
they are not sources to apply again.

The checked source owner snapshots both the active laws and last accepted
report. It admits active fluid temperatures, prescribed wall temperatures,
and every resolved wall trace against the caller's scalar interval. A late
refusal preserves the original populations, histories, wall laws and report.
The channel dynamic SGS filter retains declared temperatures and zero-flux
conditions independently at each wall. For nonzero flux/film or nonuniform
spatial walls it uses
three-fluid-centre reconstruction of both physical and Favre scalar gradients.
The original reflected-pair operator still applies actual surface heat, and
the same source transaction admits temperature, Mach, cancellation and heat
balance before publishing state. The filter neither invents a wall trace nor
adds modeled mechanical dissipation as a heat source. Its zero-diffusion path
retains ordinary transport and local molecular rates; physical-unit and
source-resolved material adapters use the same coupling. Other ordinary/source
thermal steps share the original wall streaming operator in both dimensions.

`MaterialThermalWall` provides the same laws in SI through the retained
material slab, including local fields through `set_spatial_thermal_walls`.
It converts signed `q''` by `dt/(rho*Cp*dx*DeltaT)` and `h`
by `dt/(rho*Cp*dx)`; reversing the original temperature coordinate does not
reverse physical heating. The wall report converts actual link exchange to
joules using the same captured material heat capacity and physical cell
volume, with explicit extrusion depth in 2-D and `dx^3` in 3-D. Surface
mean/minimum/maximum temperatures are reported in kelvin, including proper
extrema ordering for negative `DeltaT`. The existing frozen-storage and
Boussinesq assumptions remain. Films use authored coefficients and prescribed
ambients; correlation selection and finite external storage remain with their
existing owners.

Focused kernel regressions cover the original temperature arithmetic, exact
insulation, signed imposed flux, local film balance, and the analytic series
resistance of a conductive slab between two films across several relaxation
times. Source and channel-SGS regressions cover boundary publication and
refusal. Spatial regressions compare constant fields bit for bit with uniform
laws, preserve wall patterns and derivative history after rejected steps, and
check opposed buoyant plumes driven by balanced local surface heating/cooling
in 2-D and 3-D. This closure is limited to the existing planar, uniform-cell,
constant-capacity temperature slabs; it adds neither curved boundaries nor
a finite ambient reservoir or a general unsteady accuracy guarantee.

### Full tensor diffusion in the original thermal slabs

`thermal::tensor::DiffusivityTensor<D>` prepares the full symmetric positive
definite Cartesian diffusivity for the existing D2Q5/D3Q7 temperature owners.
Its units are lattice-cell squared per lattice-step. The physical conversion
is `alpha = K/(rho Cp) * dt/dx^2` for a square/cubic grid and common volumetric
heat capacity. This conversion supplies neither a material source nor
variable-capacity storage. Exact input symmetry, positive diffusivity,
representable inversion, and both `S > 0` and `2I-S > 0` are checked without
clipping. Existing fs-la Cholesky solves prepare
`S=(I/2+alpha/c_T^2)^(-1)`, with `c_T^2=1/3` in 2-D and `1/4` in 3-D.

`step_with_sources_and_tensor_diffusion` uses the original force frame,
scalar equilibrium, source kicks, retained advective-flux history, streaming,
wall measurements and transactional acceptance. Only the vector flux
relaxation becomes a matrix; even kinetic moments retain the original
`tau_g`. The correction applies `I-S/2` to the existing time difference of
`T u`, preserving the intended diffusion tensor during advection.
The matrix form follows Zhao et al. (2020), equations (5) and (13),
<https://arxiv.org/abs/1902.09813>; the D3Q7 tensor construction follows
Yoshida and Nagaoka (2010), <https://doi.org/10.1016/j.jcp.2010.06.037>.

A uniform field reuses one prepared tensor; a per-cell field retains the
supplied axes and spatial variation. No factorization, lattice allocation or
tensor re-preparation is added to the collision loop. The existing checked
source owner still allocates its transactional candidate. Tensor fields are
borrowed for one interval and do not replace background scalar coefficients.
Malformed fields and rejected endpoints retain the original populations,
forces and derivative history.

The focused tests exercise rotated Fourier decay against `k^T alpha k`,
spatial refinement with nonzero advection, reflected axes, nonpassive and
unrepresentable tensor refusal, per-cell/uniform parity, source and wall heat
balance, failed-candidate retry and subsequent scalar evolution. SPD
admission is not a population-positivity or coupled-buoyancy stability theorem.
The original low-Mach, startup-history, uniform-storage and wall-model
limitations remain; sharp material interfaces, variable heat capacity and
unrestricted high-anisotropy accuracy require their own physical treatment.

### Original material sources for thermal tensor diffusion

`material::transport::tensor::ThermalMaterialField` supplies uniform or
per-cell original `PreparedOrthotropicThermalTransport<D>` paths to
`MaterialSlab::step_with_material_tensor_diffusion` in two and three
dimensions. Principal directions are expressed in the square/cubic lattice
frame. The retained SI mapping supplies `alpha_lattice=alpha_SI*dt/dx^2`;
the existing tensor preparation and collision remain the numerical owners.

Both original density and specific-heat paths must be constant at every
canonical source knot and match the slab's retained physical storage exactly.
Endpoint equality alone is insufficient. Conductivity may vary with actual
temperature and declared cell assignment. Initial coefficients are frozen
for one explicit step; endpoint samples qualify source support and tensor
representation before adoption. They do not constitute an implicit solve.

Power, force, wall heat and joule balances use the original physical adapter.
Its existing acceptance callback performs endpoint qualification and
cancellation checks inside the single native candidate transaction; there is
no additional full-grid snapshot. Source records borrow their original owner,
retaining every receipt without per-cell source clones. Adjacent identical
evaluated tensors reuse their factorization. No factorization is added to
the population collision loop.

Five integration tests cover original-source authority, scalar limits in both
dimensions, temperature-dependent/per-cell transport and SI heat balance,
rotated 3-D flux, and endpoint/storage/cancellation refusal with replay.
The source-only native subset passes three tests; the full LBM integration
battery remains unexecuted in the constrained development environment.

The caller's declared orthotropic symmetry and fixed orientation remain
explicit. Variable storage, latent heat, contact resistance and evolving
material axes require their own physical models. Momentum and buoyancy retain
the existing background setup; scalar background diffusivity/Nusselt
diagnostics do not become tensor diagnostics.

### Independently fitted dynamic transport on the thermal slab

`ThermalLbm3::{thermal_channel_sgs_workspace, prepare_thermal_channel_sgs,
step_with_sources_and_channel_sgs}` explicitly attaches the shared planar
dynamic momentum/scalar model to the original D3Q19/D3Q7 slab. Only active
fluid cells enter the filters. Its current Guo-centered velocity, density,
temperature, stationary halfway walls and actual fixed plate temperatures
determine both fits. Momentum and thermal diffusion are independent, including
planes with zero modeled viscosity and positive modeled heat diffusion.
The original molecular rates or supplied local molecular fields remain the
baseline; zero fitted contributions preserve each corresponding rate's bits.

The original transactional source/flow/temperature owner consumes the total
rates for one step and measures the anti-bounce-back wall heat exchange.
`MaterialSlab<ThermalLbm3>::step_with_power_density_and_channel_sgs` uses the
same retained material, SI power/force conversion and normalized temperature
map, including either sign of the physical plate temperature difference.
Source and wall heat, temperature admission and endpoint Mach checks all
precede publication. Rejected or cancelled trials retain the accepted physical
state, and a retry overwrites failed scratch. The original serial collision
kernel is checked before/after, without claiming finer cancellation inside it.

The reusable numerical workspace requires `240*N_fluid + 16*N_slab` nominal
bytes, with actual-capacity admission. This includes current rho/u/T, both
filters and full step-local rates, but excludes the existing source owner's
complete trial snapshot. Preparation allocates no additional arrays. Explicit
coefficient, viscosity/diffusivity-to-molecular and Mach bounds are refusals.
Shear dissipation diagnostics reuse the original thermal velocity stencil and
do not add mechanical dissipation to the sensible-temperature equation.

The original wall discretization applies adjacent-fluid transport through
halfway bounce-back/anti-bounce-back; it is not a new molecular-only wall
resistance or wall function. The six focused native channel tests cover
independent heat, original source/wall transport, force-centering, zero-fit
equivalence, budgets, cancellation/retry, endpoint refusal and both SI thermal
orientations. This remains exploratory, resolved-wall Boussinesq transport;
it supplies no total-energy or developed turbulent-channel/DNS qualification.

### Source-resolved molecular transport with thermal SGS

`MaterialSlab<ThermalLbm3>` now composes source-resolved transport with both
local SGS and independently fitted channel momentum/heat. The
`step_with_material_transport_and_sgs` and
`step_with_material_transport_and_channel_sgs` methods use one original
material query per fluid cell at the accepted temperature and one at the
candidate endpoint. Their concrete source records, SI heat balance and SGS
report return together as `MaterialSgsTransportStep`.

Fresh `mu/rho0` and `k/(rho0*Cp0)` rates enter the existing SGS preparation
and provide its molecular-ratio bounds. Augmented rates are step-local and
never replace retained material properties or accumulate across steps. The
accepted Guo force frame supplies the incoming velocity. Original stationary
wall transport, source splitting, independently measured wall heat, endpoint
Mach admission and final material publication retain their existing owners.

The corresponding `step_with_boussinesq_material_transport_and_sgs` and
`step_with_boussinesq_material_transport_and_channel_sgs` methods also use
that same sample's density in `(rho_source(T)-rho0)*gravity`. This replaces
the legacy beta contribution on the tentative step. Both material epochs
must satisfy the authored density/Cp departure bounds; reference inertia and
sensible storage stay fixed. Refusal or cancellation publishes no partial
state. Endpoint sampling checks source support, not implicit convergence.
The coefficients remain explicit, and no modeled dissipation heat or
total-energy qualification is introduced.

The explicit `step_with_boussinesq_material_transport_and_sgs_endpoint_averaged`
and `step_with_boussinesq_material_transport_and_channel_sgs_endpoint_averaged`
variants return `MaterialSgsEndpointStep`. The original SGS owner prepares total
rates once from accepted Guo-centered velocity, accepted scalar state where
needed, and fresh initial molecular mu/k. The existing density endpoint
projector consumes the one scalar transport performed under those total rates;
its original bounded source queries and density departures drive the one
subsequent flow step. The report retains initial, force-predictor and
independently sampled actual-final material epochs. Initial molecular and SGS
coefficients remain frozen for that step; no augmented rate becomes a
material property.

Both families retain original ratio/usage limits, reusable scratch,
independent dynamic momentum/heat fits, wall/source heat checks and post-flow
Mach admission. Typed source/SGS failures and cancellation publish no candidate.
Fixed rho0/Cp0, sensible storage and the existing no-global-time-order/no-total-
energy boundaries remain unchanged.

Three additional G0/G1/G4 regressions cover zero-SGS parity, nonlinear source
density driven by nonzero modeled thermal transport, exact query/report
identity, original force-frame preparation, independent positive channel heat,
scratch reuse and late-source/Mach/cancellation rollback with nonuniform bulk
relaxation. These three, the five existing material SGS checks and two native
force-resolution checks passed in a bounded original-source build with actual
material, SGS and numerical owners. Unrelated registries/executor scheduling
and the full Cargo dependency graph were not built.

Five focused regressions are included for native-step equivalence, evolving
molecular transport, independent channel heat, shared source/force ownership,
ratio bounds and endpoint/cancellation rollback. All five pass in the native
debug build after execution was restored.

### No-claim boundaries (3-D coupled ports)

- Single-threaded correctness-first reference rung: no SIMD/tile/parallel
  claim, no performance claim against the tiled fast path.
- The baseline battery establishes the fixture-scale onset band and profile
  bars it actually gates; the explicit SGS paths require their separate
  exploratory controls and do not qualify developed turbulence.
- The seeded determinism hash is a replay-checked CANDIDATE, not a frozen
  golden: freezing awaits the GOLDEN_POLICY four-quadrant ceremony
  (committed tree, M4 + x86, debug + release), per the 40p2 precedent.
- Wall layers hold equilibrium state and are excluded from physics
  claims; only the fluid interior is certified.

## Steady conjugate heat transfer (`conjugate`, bead rc-root-q61wp.34)

`conjugate` couples the existing `BoundaryGrid3` airflow to ONE conservative
finite-volume energy equation over fluid AND solid voxels, for forced-air
electronics cooling (the wedge's second physics after steady conduction).

- `VoxelDomain` (cubic cells of edge `dx` metres, index `(z ny + y) nx + x`,
  `Voxel::Fluid | Voxel::Solid(material)`), `FluidProperties` (constant
  `rho`, `c_p`, `k`, `nu`; `dry_air_300k()` cites Incropera Table A.4) and
  `SolidMaterial` (labelled isotropic `k`).
- `FlowField`: staggered volumetric face fluxes (m^3/s). Producers:
  `quiescent`, `from_face_velocity` (analytic/measured face-normal velocity,
  no projection), and `from_cell_velocities`, which interpolates cell
  velocities and PROJECTS them onto discretely divergence-free fluxes with
  one SPD correction-potential solve on the fluid cells (`FlowFace::Fixed`
  fluxes held, potential zero on `FlowFace::Free` faces, components without
  a free face pinned). `ProjectionReport` retains divergence before/after
  and the largest correction. Faces touching a solid carry exactly zero.
  Before pinning a component, the projection checks that its prescribed
  boundary inflow and outflow balance. It refuses incompatible fixed fluxes
  per connected component, even if separate pockets' imbalances cancel
  globally: a potential gauge is not a physical mass source. The check uses
  a compensated boundary-only sum with a `128 epsilon * sum |Q|` floating
  point allowance, independently of the iterative solver tolerance. G0/G3
  tests exercise both flow directions, uniform and graded grids, separated
  pockets, and preservation of uniform temperature under balanced forced
  convection with both upwind and power-law energy schemes.
- `FlowField::project_density_volume_source` projects an existing native face
  field with a positive frozen density per fluid cell and an authored net
  outward volume rate per cell (m^3/s). It uses the same component graph,
  native ILU(0)/BiCGStab solve, and face-correction path. Face inertia is the
  actual staggered dual mass `M = A (rho_L h_L/2 + rho_R h_R/2)`; the pressure
  coefficient is `dt A^2/M`, including boundary half-cells and graded grids.
  Uniform density retains the original geometric-potential arithmetic exactly
  and converts the solved potential to physical correction pressure afterward.
  `Wall`/solid flux is zero, `Fixed` flux stays prescribed and `Free` correction
  pressure is zero. Each component without a free face must independently
  balance its volume source with its prescribed boundary outflow; pinning only
  removes the pressure gauge and never subtracts an incompatible source mean.
  The result contains actual corrected fluxes, cell pressure correction, and
  independent per-cell volume-defect and mechanical-work evidence. Before
  publication it admits `delta_K + sum(M delta_u^2/2) = dt sum(p source) -
  dt sum_fixed_boundary(p Q_out)` within explicit caller budgets; kinetic
  change is evaluated from the actual flux increment without subtracting
  large energy totals. The nonnegative projection loss remains numerical and
  never becomes heat. The public method leaves its input unchanged on invalid
  source/density/flux, incompatible components, failed solve, exceeded budgets
  or cancellation. Dimensions/minimum spacing and held-face semantics are
  checked; original domain geometry/occupancy identity remains the caller's
  obligation because `FlowField` does not retain a complete domain snapshot.
  This is a working frozen-density pressure/continuity component, not a phase
  or density timestep: the caller still owns density evolution, mass/momentum
  sources and source/boundary energy. At constant density a volume source
  implies mass source `rho source`; closed-material phase expansion needs a
  compatible evolving-density and enthalpy owner. Existing fixed-density
  phase/thermal gates remain unchanged. Focused native regressions cover
  analytic graded two-density pressure/flux/work on all axes, exact homogeneous
  legacy parity, fixed-boundary work, sealed source/sink compatibility,
  disconnected components, independent readback budgets and refusal rollback.
- `project_density_volume_source_with_boundary_pressure` extends that same
  native solve with finite authored pressure on each `Free` domain plane.
  Pressure acts through the boundary half-cell inertia and can accelerate
  inflow or outflow; specifying pressure on a velocity-prescribed `Wall` or
  `Fixed` plane is refused. Each connected open region solves relative to its
  first reservoir pressure, restores the authored gauge in cell readback, and
  measures all pressure work in the same relative gauge. This retains small
  pressure differences at a large common pressure without tolerance-driven
  false convergence or catastrophic pressure-work cancellation. The additional
  signed reservoir work is `dt sum(p_boundary Q_out)`; the independent identity
  is `delta_K + projection_loss - source_work + fixed_work + reservoir_work = 0`.
  The default zero-pressure entry retains its original flux arithmetic. Analytic
  graded two-density tests cover both directions on all axes, homogeneous and
  heterogeneous large-datum invariance, independent disconnected reservoirs,
  and refusal rollback. This is a finite-duration mechanical pressure impulse;
  it neither supplies a viscous steady-flow law nor updates a thermodynamic EOS.
- `StaggeredBoundaryPressure` and
  `project_density_volume_source_with_pressure_profile` resolve a separate
  reservoir pressure at each native x/y/z boundary face. The caller-owned
  optional-pressure buffers must match the existing staggered layout exactly;
  `Some` is finite and allowed only on a fluid-exposed `Free` exterior slot.
  `None` supplies zero there; interior, solid, `Wall` and `Fixed` slots must
  remain `None`. Admission precedes field cloning, with no profile allocation.
  The original half-cell inertia, solve and work audit use each actual slot;
  component references remain independent even for disconnected regions on
  the same plane. Constant/zero profiles retain the original plane API bits.
  G1/G3 regressions reproduce affine-gradient acceleration and reservoir work,
  sharp-stratification hydrostatic rest through the actual gravity predictor
  on graded grids in all axes, and simultaneous inflow/outflow on one plane.
  G0/G4 checks cover native shape/role/finite admission, independent large
  references and transactional refusal. Continuous-profile face averaging,
  time-dependent reservoirs and thermodynamic state remain caller-owned.
- `FlowField::predict_density_body_force` supplies native staggered momentum
  prediction from frozen physical density, uniform acceleration and an explicit
  cell body-force density. Each face integrates the same two half-cells as the
  pressure dual mass: `F = M*g + A*(h_L*b_L/2 + h_R*b_R/2)` and
  `delta_Q = A*dt*F/M`. Boundary and solid interfaces use their fluid half-cell.
  Wall/solid and prescribed `Fixed` fluxes remain held, with the balancing
  constraint impulse and its work reported explicitly. Actual returned fluxes
  pass independent impulse and kinetic-work budgets; no density mean or
  Boussinesq substitution replaces the admitted inertias.
  `force_and_project_density` composes that predictor with the original density
  pressure/volume constraint. It reports predictor-stage and endpoint-midpoint
  force/constraint work separately, plus the signed force-pressure split term.
  At resolved hydrostatic rest the endpoint force work vanishes and this cross
  term cancels the intermediate projection loss; neither is material heat.
  Original source-pressure work, fixed-boundary pressure work and projection
  loss retain their separate native accounts. The predictor can directly feed
  material or thermal-volume transport. This is a frozen-density force/pressure
  split: advective face inertia, cell/face reconstruction, force constitutive
  closure and gravitational-potential transport remain with their physical
  owners. Native G1/G3/G4 tests cover graded stratified hydrostatics, sharp
  density contrasts, independently derived density-driven circulation, held
  boundary reactions/work, zero-force parity and atomic budget/cancel refusal.
- `FlowField::project_and_advect_material` consumes that pressure projection
  through conservative backward-Euler donor transport on the SAME native face
  volume fluxes. Cell mass is the authoritative inventory; `rho=M/V` supplies
  the pressure coefficient without reconstructing an accepted mass from rho.
  Explicit volume injection supplies its own incoming density and caloric h;
  extraction and boundary outflow use the accepted local donor state. Incoming
  and backflow boundaries require an authored rho/h pair in the same datum.
  The existing sparse row scaling and Krylov owner solve the density equation;
  its donor densities define transferred masses once. Paired transfers create
  actual endpoint masses, independently read back in that original equation.
  A second native donor solve carries `mass*h` on exactly those same transfers.
  Component-wise caloric centering excludes arbitrary enthalpy datums from the
  solve tolerance while retaining mass-closure terms. Actual stored h and all
  local/global mass and enthalpy balances pass explicit physical-unit budgets.
  No negative inventory clipping, mass renormalization, or numerical-to-heat
  conversion occurs. Allocation/connection/iteration limits and cancellation
  precede candidate publication. Focused native G1/G3/G4 cases check independent
  graded two-cell mass/energy solutions on every axis, volume-driven material
  injection/removal, large datum invariance, continuing closed circulation,
  loose-solve readback refusal and atomic input preservation. This is an
  operator-split pressure/material-advection component; the enclosing physical
  owner supplies momentum prediction, source momentum/work and compatible
  material/EOS closure. It does not infer phase expansion, mixing heat,
  pressure heat or a complete Navier–Stokes/enthalpy timestep. First-order
  donor diffusion and splitting error require spatial/temporal refinement.
- `MaterialAdvectionStep::advect_momentum` consumes the exact accepted mass
  transfers of that step, including phase-volume results. Private retained
  routes and the already assembled mass*h operator serve all three momentum
  components without another density solve or a reconstructed mass flux.
  Original extensive momenta and prescribed pre-advection impulses remain
  explicit, as do incoming boundary/source velocities. Component-centered
  native solves provide donor velocities; paired transfers form authoritative
  returned momenta, whose actual P/M independently passes the original donor
  equations. Zero transfers and zero impulses retain original momentum exactly.
  Local/global momentum and kinetic-energy gates include impulse representation,
  external carry, nonnegative backward-Euler/upwind loss and the separate
  energy contribution of retained mass roundoff. Numerical losses never change
  enthalpy. The staggered projection retains pressure correction and its work;
  cell momentum does not receive that correction again. Force splitting, phase
  mobility and feedback from staggered pressure to cell momentum remain with
  the enclosing owner.
  Altered public masses, missing incoming velocities, numerical budget failure
  or cancellation returns no candidate. Native G1/G3/G4 regressions check
  analytic graded transport on all axes, injection/extraction, impulse work,
  exact zero-flow continuation, closed circulation, velocity-datum invariance,
  actual phase-volume consumption and independent roundoff-energy refusal.
- `FlowField::reconstruct_material_momentum` transfers authoritative cell mass
  and extensive momentum onto the original staggered pressure inertias. Each
  cell contributes `M/2` and `P_axis/2` to each corresponding face; internal
  faces merge those quantities and divide by the native half-cell dual mass
  computed from `rho=M/V`. Every fluid half-cell participates, including both
  solid-interface orientations. Returned densities accompany the unprojected
  fluxes for the enclosing owner's compatible force/pressure/volume solve.
  Actual mass, momentum and kinetic-energy readbacks use independent budgets.
  Stable relative velocities expose nonnegative averaging loss. Held fluxes
  report support/controller impulse, reservoir work `I*u_held`, and the
  nonnegative constraint-projection loss; stationary walls perform zero work.
  Native dual-mass roundoff has a separate signed energy account whose absolute
  magnitude is still bounded, so naming it cannot bypass admission. None of
  these numerical losses becomes heat. The transfer is linear in cell count,
  with no additional sparse solve. It consumes actual transported cell momenta
  but supplies no reverse pressure feedback or complete Navier-Stokes timestep;
  repeated reconstruction's numerical diffusion remains explicit. Native
  regressions cover graded weighted velocities and mixing loss, large velocity
  datums, held/solid reactions, mass-representation refusal, actual accepted
  material-momentum transport followed by native pressure, and an independent
  density-weighted circulation solution.
- `FlowField::force_project_material_momentum` composes reconstruction and ONE
  fresh native force/pressure correction, then returns the actual impulses to
  the original cell momenta. Each contributing half-cell receives its share
  in proportion to authoritative `M/2`; original sub-face velocity differences
  remain intact when there is no impulse. Held reconstruction reactions reach
  the cells, including both orientations of a solid interface. Actual returned
  momentum independently passes linear and orbital angular-impulse budgets,
  using original face and cell centers on graded rectilinear geometry.
  Authored forces and returned pressure supply impulses; `Q/A` representation
  errors are audited numerical defects, never invented forces. The original
  reconstruction, force and pressure accounts remain separate. Incremental
  cell energy adds a nonnegative half-cell scatter loss and a signed native
  versus authoritative mass-representation correction; its absolute size is
  budgeted. The original reconstruction averaging loss cancels from this
  incremental balance and is not charged again. Centered local work accounts
  retain small changes under large common velocities. Numerical losses never
  alter material enthalpy. No arbitrary public projection receipt is accepted:
  all stages execute here with the same admitted mass, density and geometry.
  The caller supplies a FRESH volume target, compatible material/EOS evolution
  and phase mobility; already applied phase-volume pressure must not be applied
  again. Reconstructing returned cell momentum need not reproduce the final
  divergence-constrained face field. This force/pressure split does not provide
  staggered advective inertia or a complete Navier-Stokes timestep. Native
  regressions cover analytic density-driven circulation, graded angular
  impulse, held reactions, large velocity datums, zero-impulse preservation,
  independent representation-work refusal and repeated accepted material
  transport followed by fresh forcing, pressure and cell feedback.
- `force_and_project_density_with_pressure_profile` and
  `force_project_material_momentum_with_pressure_profile` carry the original
  borrowed `StaggeredBoundaryPressure` through native forcing, projection and
  cell feedback. Every exposed reservoir uses its actual half-cell inertia
  and its own pressure difference. Reservoir work joins the existing combined
  force/pressure identity, and physical pressure impulses reach original cell
  momentum once. The native solve's relative pressure, component references
  and existing cell-to-row map survive privately through readback; restoring
  a large public pressure gauge cannot erase a small driving head. This adds
  no pressure solve, profile copy, index-map allocation or component walk.
  The original zero-pressure entrypoints retain their readback arithmetic.
  Native regression sources specify unequal-density pressure acceleration,
  spatial hydrostatic heads, large common pressure references, cell impulse
  and energy feedback, zero-profile parity and invalid-profile refusal.
- `FlowField::advance_coexistence_volume` supplies a bounded physical consumer
  of that conservative transport: one mobile, common-velocity material mixture
  on ONE adjacent isothermal solid/liquid coexistence interval of its original
  `EquilibriumEnthalpyPhaseCurve`. The original chart gives affine specific
  volume `v(h)=a+b*h` on that interval, independent of its sensible-density
  interpolation rule. Authoritative original masses and phase states must
  fill their cells within an explicit occupied-volume budget. Caller-prescribed
  NET material enthalpy increments give `h_star=h+delta_H/M` and the physical
  expansion/contraction rate `S=M*(v(h_star)-v(h))/dt`. This source enters native
  pressure continuity with ZERO local material injection. The same projected
  donor fluxes then carry actual mass and mass*h. Incoming/backflow densities
  come from original-chart states in the same interval, preserving the affine
  closure at the endpoint. Actual returned mass and chart density independently
  pass the occupied-volume gate; tolerated initial EOS roundoff is propagated,
  never converted into a physical volume source or repaired by renormalization.
  Requested-versus-represented caloric increments and native transport defects
  share an explicit enthalpy budget. The returned `dt*sum(p*S)` source pressure
  work is an ADDITIONAL mechanical account in the native correction-pressure
  gauge; fixed-boundary work and numerical projection loss remain separate in
  the native projection report. The prescribed increments are not raw supplied
  heat, and no numerical loss becomes latent or sensible heat. The caller owns
  the thermal/work closure, momentum predictor and admissible phase mobility;
  thermodynamic pressure, composition and source applicability remain fixed.
  This substep does not infer pressure-dependent EOS feedback, solid mechanics,
  mushy drag or free-surface geometry. Other charts, sensible-branch crossing,
  incompatible sealed expansion, failed solves, budgets or cancellation refuse
  without publishing a candidate. Native G1/G3/G4 cases check analytic melting
  outflow and freezing inflow, continued inventories, balanced sealed transfer,
  uniform-density parity, caloric-datum invariance, zero-input continuation with
  admitted roundoff, independent endpoint EOS readback and atomic refusals.
- `FlowField::advance_thermal_volume` extends the same native owner to sensible
  heating/cooling on one adjacent original chart interval. The chart's explicit
  density-law query distinguishes affine specific volume from sourced affine
  sensible density. Affine-volume intervals retain the original one-pass
  closure, including exact coexistence parity. For affine density, conservative
  enthalpy mixing changes occupied volume: each original/incoming donor mass
  contributes `-m/rho_d*((rho_d-rho_r)/rho_r)^2` at its receiver. The exact
  accepted native mass routes supply these weights. Bounded relaxed iteration
  updates only the continuity source; original mass, flow, boundary states and
  prescribed net enthalpy increments remain fixed, with the caloric increment
  applied once. Independent endpoint mass, enthalpy and EOS-volume budgets
  admit the result. Trial EOS queries allocate/hash no new phase identities;
  only the accepted endpoint receives original-chart phase states and pressure
  work. The mixing source is not heat, and initial EOS roundoff is never an
  expansion source. Fixed thermodynamic pressure and one original interval
  remain explicit limits: incompatible sealed volume, interval crossing,
  exhausted iterations, numerical budgets or cancellation refuse atomically.
  Thermal/work closure, momentum prediction, phase mobility and free-surface
  geometry remain with their existing owners. Native G1/G3/G4 tests include
  an independent nonlinear inflow quadratic solution, graded internal mixing,
  sensible heating/cooling on every axis, continuation, affine parity,
  caloric-datum/relaxation invariance and bounded refusals. First-order donor
  mixing still requires spatial and temporal refinement.
- `FlowField::advance_piecewise_thermal_volume` admits the full original chart
  through the same pressure, mass and extensive-enthalpy transport owner.
  Heating, cooling and donor mixing may cross sensible/latent boundaries and
  distinct original sensible-density slopes. A receiver's mixing source sums
  `m*(v(hr) + v'(hr)*(hd-hr) - v(hd))` over original inventory and the actual
  accepted incoming mass routes. Conservative enthalpy balance cancels the
  common receiver tangent; the chart evaluates signed slope changes and local
  reciprocal curvature without fitting another EOS. Original inputs remain
  frozen across trials and the prescribed net caloric increment is applied
  exactly once. `PiecewiseThermalVolumeControls` bounds both nonlinear trials
  and TOTAL crossed intervals across every donor query and every source trial;
  cancellation and independent occupied-volume/conservation gates publish no
  partial state. The existing single-interval fast paths remain intact.
  Chart conditions and source identities remain fixed. The caller still owns
  phase mobility and thermal/work closure; this does not supply a thermal PDE,
  solid mechanics, pressure-dependent EOS or free-surface geometry. Native
  regressions cover heating/cooling through both boundaries on every axis,
  independent cross-branch mass/enthalpy/EOS solutions, graded internal mixing,
  both signs of knot-volume corrections, distinct affine-density branches,
  large caloric datums, exact total-work limits and atomic refusals.
- `FlowField::advance_isotropic_conductive_thermal_volume` connects an original
  `PreparedPhaseConductivity` to that full-chart volume owner. Calling it
  explicitly declares the scalar source isotropic. One source sample per fluid
  cell supplies original `T(h)`, phase and conductivity; the same borrowed
  caloric chart supplies density throughout the accepted transport. Each
  internal fluid/fluid face exchanges paired heat through its physical series
  resistance `G=A/(d_left/k_left+d_right/k_right)`, with graded half-cell
  distances. Zero conductivity, solid interfaces and exterior thermal faces
  insulate. The explicit diffusion step requires each cell's sum of actual
  nonnegative enthalpy-secant weights `Q_ij/(M_i*(h_j-h_i))` to be at most one;
  true latent plateaus need no apparent heat capacity, clipping or hidden
  thermal substeps. Excess duration or bounded face-work exhaustion refuses.
  Conductive heat is added once to the caller's net caloric increments, then
  native pressure, mass, enthalpy and nonlinear occupied-volume closure run
  unchanged. Paired-heat imbalance and all heat-addition residuals share the
  same enthalpy budget as caloric representation and native transport. Original
  source/chart domains admit initial, incoming, heated and accepted states.
  Required source-pressure work remains explicit and separate; the caller
  retains momentum, phase mobility and thermomechanical closure. This is an
  explicit-conduction/implicit-donor split requiring refinement, with no new
  pressure-dependent EOS, solid mechanics or free-surface evolution. Native
  regressions check analytic sealed heat-driven density flow on every axis,
  graded resistances, liquid/latent crossing and continued inventories,
  constant-density reduction, shared heat budgets and atomic refusals.
- `lbm_duct_flow` -> `LbmFlow { field, velocity_m_s, pressure_pa, report }`:
  velocity inlet on `x-min`, pressure outlet on `x-max`,
  no-slip elsewhere and on solid voxels; derives `nu_lat`/`tau` from the
  declared physical state and refuses (`LatticeResolution`, with the voxel
  size that would admit it) below `min_tau`; `Auto` collision is BGK (on
  the plate-fin heatsink at `tau = 0.519` the central-moment operator
  diverged at step 85 while BGK stayed bounded, so `Auto` never switches to
  it). A lattice that leaves its physical regime refuses (`FlowDiverged`)
  instead of panicking: `BoundaryGrid3::step_pooled` checks the density of
  every open-face source cell on the unpublished state
  (`BoundaryStepError3::Unphysical`).
  The steady criterion is the relative L2 change of the cell MASS flux
  `rho u` between checks. Mass flux, not velocity, feeds the projection:
  the weakly compressible lattice conserves `rho u` per layer (measured
  uniform to 1e-6 over interior layers) while `u` rises along the
  pressure drop. Both on-site faces are `Free`, because boundary-cell
  momentum is a reconstruction; `LbmFlowReport` retains the NOMINAL and
  REALIZED inflow (inlet rim cells share wall links, so a coarse grid
  admits less than the plug). `pressure_pa` is the gauge pressure
  `c_s^2 (rho_lat - 1) rho S^2`; `report.pressure_drop_pa` spans the first
  to the last interior layer and `mean_pressure_x` averages one layer.
- `solve_energy`: every face flux is `J = d T_P - a T_N - r`; interior
  faces use the harmonic face conductivity (exact series resistance of two
  half cells, so temperature and normal flux are continuous across any
  fluid/solid or solid/solid contrast) and power-law, exponential, or upwind
  convection in conservative form. `ThermalFace` rules: `Adiabatic`,
  `Temperature`, `HeatFlux` (into the domain), `Convective { h, ambient }`,
  `Inflow { temperature }`, `Outflow { backflow_temperature }`; solid cells
  on open faces are adiabatic. Sources via `ThermalSetup::power_w` (W per
  cell, `add_uniform_power`) and fixed-temperature cells. The M-matrix is
  solved by ILU(0)-right-preconditioned BiCGStab on the Jacobi-scaled
  system; acceptance is the RECOMPUTED true residual.
- `ConvectionScheme::Exponential` uses the same conservative face assembly in
  both energy and staggered momentum, with `A(P) = |P| / expm1(|P|)`. It
  reproduces source-free constant-coefficient 1-D profiles across graded
  centre distances and boundary half cells. Stable small/large-Peclet
  evaluation retains nonnegative coefficients without exponential overflow or
  premature underflow. Power-law remains the cheaper default; exponential
  fitting makes no multidimensional higher-order or turbulence claim. The
  numerical-leaf tests check exact face fluxes and limiting coefficients;
  `exponential_convection` checks the original CHT boundary-value solve.
  All schemes form owner and neighbour coefficients from one nonnegative
  fitted conductance, so reversing a heat face swaps them bit-for-bit rather
  than cancelling inflow diffusion. Exponential fitting evaluates the
  dimensioned conductance before underflow; a vanishing dimensionless weight
  does not erase representable physical transport. G0/G1/G3 tests cover
  paired heat cancellation, counterflow at large Peclet number, dimensional
  rescaling, and the existing graded analytic profile. CHT semantics is 4.
- `EnergyBalance` is computed from boundary and fixed-cell fluxes alone:
  `source + fixed_injection - boundary_outflow` closes to roundoff because
  the face fluxes telescope; `EnergySolution` adds `bulk_temperature_x`,
  `solid_to_fluid_heat_w`, `max_where`.

`BoundaryGrid3` streaming reads a pull-source table compiled with the link
masks (one `u32` per cell and direction) instead of re-classifying every
link on every step. The table selects exactly the source the classifier
would, so results are bit-identical: `BOUNDARY_GOLDEN_HASH`
(`d3q19_boundaries`) and `GOLDEN_HASH` (`d3q19_battery`) are unchanged and
`D3Q19_BOUNDARY_BIT_SEMANTICS_VERSION` stays 1. Memory cost: 76 bytes per
cell beside the 304 bytes of populations.

`BoundaryGrid3::step_pooled(pool, gate)` runs collide and pull-stream as
`fs_exec::TileKernel`s over groups of 8 lattice tiles (512 cells); the
open-face reconstruction (face cells only) stays serial. Serial and pooled
streams call the same `StreamView::pull_cell`, and collision the same
`collide_cell3`, so the published state is bit-identical to `step()` for
any worker count (`lbm3bc-pool`: workers 1/2/3 on an obstructed
velocity/pressure duct and on a moving-lid cavity). The stream writes a
scratch buffer that is swapped in only after both passes complete:
cancellation (`BoundaryStepError3::Cancelled`), a collision refusal (first
cell in canonical order), or a pool fault leaves the pre-step state
published (`lbm3bc-pool-cancel`). The scratch buffer adds 304 bytes per cell
on first pooled use. `lbm_duct_flow` steps through one parked crew
(`LbmFlowConfig::workers`, `0` = every host core) and takes moments from
the populations, so a diverging run is a `FlowDiverged` refusal, not a
panic.

### Natural convection (`natural_convection`)

Two-way Boussinesq coupling in a CLOSED enclosure (walls, or periodic
lattice axes whose thermal faces are adiabatic — exact for solutions
invariant along that axis). `BoundaryGrid3::set_force_field` selects the
initial per-cell Guo force frame, added to the uniform force, without changing
stored populations. Its initial physical momentum therefore includes the
selected half force. Evolved states instead use `rebind_force_field`: the
existing Guo source transform preserves mass, physical momentum
`sum(f_q c_q) + F/2`, and full force-corrected nonequilibrium stress at fixed
physical velocity. Rebinding advances no time and adds no impulse or work;
the next collision supplies one step of the new force. Zero per-cell force
disables that contribution while retaining any declared uniform force.
`force_density` exposes their total for physical-moment sampling.

Both field APIs support BGK, central moments, and reduced cumulants. Open
velocity/pressure faces admit uniform and spatial forcing, including different
forces at a face and its first-interior neighbour. Reconstruction uses physical
velocity `u=(sum_q f_q c_q + F/2)/rho` and physical nonequilibrium stress
`Pi_neq + (u F + F u)/2` at the interior cell. It projects that stress with the
prescribed face velocity or density, then uses the original `reframe3` to enter
the destination's own Guo source frame. Zero force retains the old arithmetic.
Mixed wall/open cells still replace only open-owned links; exact target moments
are claimed only on face interiors. There is no new boundary lattice allocation.
`set_force_field_scale(s)` holds exactly `s * selected` in every cell for a
fixed spatial shape with a caller-owned time history (an imposed acoustic
pressure gradient). Like rebinding, it preserves accepted mass, physical
momentum and full force-corrected stress before publishing the scaled force:
changing scale supplies no elapsed-time impulse. The following collision
applies the newly held force once. Smooth histories sampled at each interval
midpoint supply midpoint impulses. `step_with_endpoint_force_scale` and its
pooled counterpart own a separate trapezoidal endpoint convention: collision
uses `F_start`, then populations and `F_end` publish together through the private
endpoint publisher without a held-source reframe. Scaling itself does not
advance or stage time. The shape is retained once;
subsequent scales reuse collision scratch without population/force allocation.
Open traces are reconstructed in unpublished storage in the final endpoint
force frame, preserving the prescribed physical velocity when the force changes.
A nonfinite or nonpositive reconstructed face refuses before either populations
or the endpoint force is published.
All fluid samples and total forces must be
finite. Rebinding preflights every changed cell's physical velocity and new
population frame in the existing collision scratch before publishing either
populations or field; callback failure or late invalid-cell refusal leaves
both accepted values intact. Solid cells are skipped and unchanged-force
population bits are retained. No second population lattice is allocated.
Without a field, the existing step arithmetic is unchanged. Each natural
convection coupling rebinds the buoyancy
`-beta (T - T_ref) g` in lattice units, advances up to `coupling_interval`
pooled steps, projects the physical cell mass fluxes including `F/2`, and takes one implicit
energy step over that actual batch duration, including a final batch shortened
by the original step budget. A shorter nonzero change cannot supply a
geometric-tail estimate by comparison with a longer batch. Its storage uses
the fluid's volumetric heat capacity in
every cell (solids get it as a PSEUDO capacity: it shapes the path, never
the fixed point). Convergence is judged on the geometric extrapolation
`delta r / (1 - r)` of the remaining per-coupling change (no estimate
until the sequence contracts), for the relative mass flux and for the
temperature over its range. A candidate STEADY energy solve on those same
projected fluxes must additionally satisfy the original `steady_tolerance`:
`coupling_residual_k` is its maximum distance from both the temperature that
drove the lattice batch and the pseudo-time endpoint, divided for admission
by the candidate's actual temperature range (no one-kelvin or absolute-temperature
floor). Failure continues the original pseudo-time trajectory with restarted
contraction history, within the unchanged `max_steps` budget. Success returns
that exact field, energy solution and projection without another solve.
The endpoint acceptance and actual batch duration advance
`CHT_SEMANTICS_VERSION` to 7.
The focused endpoint tests check the lagged-force defect, thermal-scale
invariance, bounded refusal and returned thermal-field consistency; the
existing cavity checks also require the normalized endpoint tolerance.
A lattice speed above `max_lattice_speed` refuses as
`LatticeResolution`; this guard and the reported speed use the actual
`u = (sum(f_q c_q) + F/2)/rho`, while flux projection retains `rho*u` over the
reference density. These corrections advance `CHT_SEMANTICS_VERSION` to 3.

`tests/boundary_forcing.rs` independently checks mass/momentum/full-stress
invariance, force switch/reversal/off impulses for all three collision laws,
finite-inertia wall balances on all three face normals, unchanged-force bits,
and atomic field/frame refusal. The convection benchmark measurements below
predate version 3 and have not been requalified after this correction.

| Fixture | Reference | Measured |
|---|---|---|
| de Vahl Davis cavity, Ra = 1e3, Pr = 0.71, 16 x 16 (periodic z) | Nu = 1.118 | 1.1127 (0.5 %); remaining-change estimates 8e-7 / 2.8e-7; coupling residual 3.5e-7 K |
| Conducting block heated in a cold enclosure | balance and buoyancy response | balance closes to 1e-9 W; peak rise 0.3621 K -> 0.3546 K as beta x5 |

Release lane (`--ignored`, executed) on 32 x 32: Ra = 1e4 (tau 0.8) Nu
2.2308 vs 2.243 (-0.55 %, 3120 steps); Ra = 1e5 (tau 0.56) Nu 4.4087 vs
4.519 (-2.4 %, 12880 steps); gated at 5 %.

### Resolved-material transient storage

`EnergyConfig::step_with_conductivity` and `step_with_enthalpy` use the
same physical backward-Euler storage rows and energy ledger as the native
march. Explicit cell capacities remain in J/K and retain one entry per cell.
Ordinary cells require positive capacity; every collapsed compact alias must
supply exactly zero so that the component's authored `heat_capacity_j_k`
is counted once. Junction aliases must share the accepted old temperature.
Fixed-temperature cells are supported: their reservoir supplies any imposed
stored-energy jump as well as the conductive and sink reaction.

The returned `EnergySolution::report.balance` contains physical sources,
reservoir reactions and sinks; storage is no longer injected into
`sink_outflow_w`. `TransientRecord` separately accounts for ordinary and
junction storage and closes the physical step. Callers of the resolved-field
stepping APIs must use that transient record for the step residual. As before,
an external nonlinear enthalpy owner must check its actual `m * delta(h)`;
a frozen capacity secant is not by itself convergence of a nonlinear chart.

`tests/compact_material_storage.rs` independently checks resolved conductivity
with closed junction equilibration, duplicate-capacity refusals, a prescribed
reservoir temperature jump with physical sinks, and source-enthalpy advection
past a storing compact component. The operator retains one common interface
heat flux and counts internal component/fluid exchange only once.

### Transient conjugate conduction (`march_energy`)

Backward-Euler steps of `C dT/dt + A T = b(t)` over a FROZEN flow field
(forced convection), with `C = (rho c) V` from the fluid and each solid's
declared `SolidMaterial::with_heat_capacity` (an ordinary solid without one
refuses). Compact components instead require an authored total
`CompactComponent::heat_capacity_j_k`; their collapsed voxel aliases carry
zero additional capacity. Both cell and component sources use
`power_w * power_schedule(t)`, with per-step records of the peak
(solid) temperature and an energy closure computed two independent ways
(stored-energy change vs `dt (source - boundary outflow)`). Positive `C / dt`
storage anchors the absolute temperature even with completely insulated or
prescribed-flux boundaries; it adds no artificial heat sink. A steady solve
still requires a physical temperature anchor. Unconditionally stable for
the M-matrix, so electronics-scale conductivity and capacity
contrasts impose no step limit; first order in time. Declared radiation
(`Some(&RadiationConfig)`, escape or surface exchange) converges emission
and irradiation at each backward-Euler endpoint, retaining the original
time level in every nonlinear trial. Both temperature change and local
physical radiation-power defects must satisfy `RadiationConfig::tolerance`
within `max_iterations`; invalid absolute temperatures or an unconverged
endpoint refuse. The per-step balance uses the physical face law.
Prescribed-cell reservoir reactions include any imposed stored-energy
jump, and their physical radiation reactions cancel their matching sink
correction in the net boundary balance. `TransientRecord::radiated_w` is
the actual external radiation loss; internal exchange in a sealed
enclosure conserves stored energy. `march_conjugate` shares this endpoint
solve. These radiation and reservoir corrections move
`CHT_SEMANTICS_VERSION` to 2.

| Fixture | Reference | Measured |
|---|---|---|
| Near-isothermal cube (Biot 4e-7) cooled by convective faces, 100 steps of tau/50 | backward-Euler recursion of the lumped ODE | 0.1380329637 vs 0.1380329672 |
| Same, dt halved | continuous exponential | error 2.698e-3 -> 1.351e-3 (ratio 2.00) |
| Heated plate under channel flow, 400 steps from cold | steady solve; energy budget | stored 2.616326564 J = net input; worst step closure 4e-12 J; within 2.5e-9 K of steady |
| Radiating floor plate (the Stefan–Boltzmann fixture with `rho c` 2.4e6), 30 steps to t = 1800 s | independent bisection of `C(T+ - T) = dt [P - eps sigma A(T+^4 - T_amb^4)]` | 379.448207585 vs 379.448207595 K; every step closes within 1e-8 of input |
| Same, 15 vs 30 steps | RK4 solution of `C dT/dt = P - eps sigma A (T^4 - T_amb^4)` | temporal errors 1.145121 -> 0.577911 K (ratio 1.9815) |
| Same, 200 steps of 600 s | closed-form steady state | 433.501643679 vs 433.501645484 K; 4.999999892 of 5 W radiated |

Additional regressions cover coarse-step heating and cooling in both
radiation modes, the zero-velocity coupled march, exhausted nonlinear
budgets, closed-enclosure stored-energy conservation, and prescribed
reservoir heating/cooling jumps with and without radiation.

### Steady finite-volume flow (`simple_flow`)

Steady laminar incompressible Navier–Stokes by finite volumes on the
staggered (MAC) faces of the same `VoxelDomain`, coupled by SIMPLEC (Van
Doormaal & Raithby 1984). This is the forced-convection path with no
lattice constraint: no `tau`, Mach or cell-Reynolds admission, so operating
points the LBM would need millions of cells for run on the declared voxels.
The returned `FlowField` holds the face fluxes the energy solver consumes,
so `simple_flow` -> `solve_energy` / `march_energy` is the conjugate chain.

- Faces touching a solid voxel are blocked; a blocked transverse neighbour
  is a no-slip wall half a cell away (staircase). Domain faces are
  `FvBoundary::{Wall { velocity }, Symmetry, Inlet { velocity }, Outlet,
  Pressure { pressure_pa }}`;
  an outlet face is a momentum unknown whose outside twin mirrors the inside
  cell (zero normal gradient) and carries the ghost pressure
  `2 p_boundary - p_inside`, so the face holds its declared static pressure
  and the pressure correction uses the matching
  linearization `u' = 2 d p'_inside`. The steady state is one fixed point of
  the discrete equations, independent of the relaxation, also for
  undeveloped (even partly reversed) outflow. (The first version
  re-extrapolated outlet velocities each iteration and stalled at a 3e-4
  mass residual on the heatsink, whose outflow is not developed.)
- `Pressure` holds a finite signed static gauge pressure [Pa]; `Outlet` retains
  exactly zero pressure. All openings and reported cell pressures share the
  caller's reference (the removed reference hydrostatic state in Boussinesq
  runs). This incompressible load is not absolute gas pressure and does not
  update density or a material EOS. Opening velocities are solved in either
  direction. `SimpleConfig::backflow_tangential_velocity_m_s` optionally supplies
  the external tangential velocity carried by incoming fluid on each opening;
  absence retains zero incoming tangential momentum. Vectors are finite, have
  zero normal component, and are allowed only on Pressure/Outlet faces. Pressure
  and continuity still determine the normal velocity. Each of the two half-face
  portions of a transverse staggered control volume selects its own donor:
  simultaneous inflow and outflow exchange momentum even when their net volume
  flux is zero. Outgoing fluid carries the solved interior velocity. This is
  advected momentum, not an imposed viscous trace, shear traction, or an audited
  reservoir kinetic-energy/heat supply. The shared steady/transient momentum
  operator uses the same condition; the normal pressure equation is unchanged.
  The existing energy owner supplies backflow temperature/enthalpy; no reservoir
  thermal state is inferred from pressure. Momentum is iterated relative to one
  fixed numerical datum and that datum is restored on public pressure output, so
  common pressure offsets do not contaminate small pressure-gradient residuals.
  Native material momentum
  checkpoints retain the datum. Face fans retain their ambient-zero reference.
  Nonfinite pressure and unrepresentable pressure differences refuse at admission.
  Focused G1/G3 `pressure_opening` tests cover planar Poiseuille refinement,
  pressure datum/offset invariance, reversed flow on all axes, source-viscosity
  scaling, the existing backflow heat transport, and transient start-up.
  Backflow-momentum regressions independently close the actual three-row MAC
  predictor's momentum budget, including mixed-sign half-face fluxes, axis/flow
  reversal, zero-default equivalence, native steady/transient response and
  unsupported boundary/work-model refusal. Test presence alone is not an
  executed validation or mesh-convergence claim.
- Momentum uses shared power-law, exponential, or upwind coefficients with `a_P / alpha`
  under-relaxation (alpha in (0, 1); SIMPLEC's `d` degenerates at 1). Closed
  components of the pressure correction are pinned.
- `SimpleConfig::momentum_reconstruction = MomentumReconstruction::Minmod`
  adds a conservative limited linear correction to the existing upwind momentum
  operator. It requires `scheme = Upwind`; the default `DonorCell` retains the
  previous operator. Both steady and transient nonlinear sweeps refresh the
  correction from actual native face velocities. Each shared side uses one
  canonical donor and equal-and-opposite correction, with physical staggered
  distances on graded meshes. Smooth interior affine velocity is reproduced;
  extrema, walls, pressure openings, solids and incomplete stencils retain
  donor transport. Molecular stress, pressure and creeping-flow owners are
  unchanged. No per-face allocation, thermal reconstruction, numerical-heat
  deposit, global multidimensional maximum principle or unconditional kinetic
  stability is introduced. Focused regressions cover graded affine recovery,
  paired conservation, velocity offsets, unsupported stencils, scheme admission
  and the original momentum operator's manufactured transported-shear residual.
  These are regression definitions; native execution and mesh/time refinement
  remain distinct evidence.
- Convergence is two steady residuals of the same iterate: the largest cell
  mass imbalance before correction over the largest face mass flux, and
  per component `||b - A u||` of the Jacobi-scaled momentum system at the
  velocities entering the iteration, over the flow's momentum scale
  `sqrt(rows) max_c rms(b_c)` (a component that vanishes by symmetry has a
  round-off right-hand side; normalizing it by its own norm left plug flows
  through a porous block stalled at a "residual" of 0.3). Inner solves only reduce their
  entry residuals (by `momentum_tolerance`, default 0.1, and
  `pressure_tolerance`, default 0.01), so a converged report never rests on a
  skipped inner solve; tight final corrections (1e-8 of each entry imbalance)
  then hand the energy equation a projected flux field. Budget
  exhaustion refuses as `FlowNotSteady`; non-finite velocities or pressures
  after any momentum solve or correction refuse as `FlowDiverged` before they
  reach an incomplete factorization.
  Both steady handoff and transient checkpoints require every fluid cell's
  native volume flux to satisfy `abs(sum Q) <= 128 * f64::EPSILON * sum(abs Q)`.
  At most four corrections reuse the same pressure operator; unresolved local
  continuity refuses without clipping fixed faces or quiet circulation. The
  added `pressure_projection` regressions cover the multiscale failure and
  caloric handoff, but have not been executed in this change's environment.
- `SimpleConfig::fan` puts a piecewise-linear `FanCurve` (static pressure rise
  against flow) on one `Inlet` face: the face's uniform normal inflow is the
  operating point where the curve meets the mean pressure of the fluid layer
  behind the face (fan pulling from ambient at pressure zero). The delivery
  moves by secant steps on settled flows only (the inlet-layer pressure just
  after an inflow change is a pressure-correction transient; a per-sweep
  update was measured to oscillate with growing amplitude), and the fan
  mismatch joins the momentum residual. The report carries the operating
  point and mismatch.
- The pressure correction is symmetric positive definite (equal face
  coefficients both ways; outlets add only to the diagonal), so
  `PressureSolver::AmgCg` (the default) solves it by conjugate gradients
  preconditioned with fs-sparse's smoothed-aggregation AMG V-cycle (strength
  0.08, three smoothing sweeps), the hierarchy rebuilt every
  `AMG_REBUILD_SWEEPS` = 10 sweeps because the coefficients drift slowly; a
  failed CG solve falls back to ILU(0)-BiCGStab, which
  `PressureSolver::IluBicgstab` selects outright. The report counts momentum
  and pressure Krylov iterations. On the 1 mm heatsink AMG-CG took a third
  of the ILU pressure iterations at no wall-time gain (16 800 cells); at
  0.5 mm the AMG flow solve took 3 285 s (debug) where an ILU run spent
  4 871 s (below; not a like-for-like timing, the projection tolerance
  differed).

| Fixture | Reference | Measured |
|---|---|---|
| Plane Poiseuille, n = 4, 8, 16 cells across | exact discrete law `12 mu U / H^2 * n^2 / (n^2 + 2)` (half-cell wall stencil; second order to the continuum) | 10.666666608, 11.636363639, 11.906976745: within 1e-8 relative |
| Lid-driven cavity, Re 100, 16^2 and 32^2 | Ghia, Ghia & Shin (1982) centreline u | worst deviation 0.0268 -> 0.0059 |
| Square duct (quarter, symmetry), half-side 4 and 8 cells | Darcy f Re = 56.91 (Shah & London) | 53.749, 56.069 (error ratio 3.8); Richardson 56.843 (0.12 %) |
| Aspect-0.5 duct, 8 x 16, Re 10, same thermal problem as the LBM rung | analytic developed profile and its Nu | profile within 0.4 %; local Nu 3.5078 vs 3.5005 (0.2 %); dp/dx -2.98 % vs f Re = 62.19 |
| Block two cells upstream of the outlet (backflow -0.099 m/s through the outlet plane) | convergence; relaxation independence | converged to 1e-10; alpha 0.5 vs 0.8 velocities within 2.6e-10 |
| Channel driven by a linear fan curve (3 Pa shut-off, 2.5e-3 m^3/s free delivery), operating mid-curve | the operating point is on the curve AND on the system curve (a fixed-velocity solve at the solved flow) | Q = 1.6917e-3 m^3/s at 0.96997 Pa; curve mismatch 6e-10; fixed-velocity inlet pressure equal to 1e-6 |
| Refusals | outward inlet, wall with normal velocity, alpha = 1, budget, all-solid domain, cancellation | structured errors |

### Internal fans and flow resistances (`InternalFan`, `FlowResistance`)

The enclosure ingredients of electronics cooling, all inside the SIMPLEC
momentum equations (momentum only; no heat):

- `InternalFan { patch, blows_positive, curve }`: across each open face of
  a `FacePatch` (an interior face plane and a transverse cell range) the
  static pressure rises by the `FanCurve` value at the flow through the
  whole patch. The source is linearized about each sweep's flow,
  `S(u) = s A dp(Q0) + A dp'(Q0) A_fan (u - u0)`: exact at the fixed point,
  and the curve's own (non-positive) slope damps the fan-system loop. The
  report lists each fan's operating point.
- `FlowResistance::Planar { patch, loss_coefficient }`: a grille or
  perforated plate dropping `1/2 rho K |u| u` across each face;
  `FlowResistance::perforated_plate_loss(f)` gives Idelchik's thin
  sharp-edged plate `K = (1 + 0.707 sqrt(1 - f) - f)^2 / f^2`.
- `FlowResistance::Volume { lo, hi, permeability_m2, inertial_per_m }`: an
  orthotropic Darcy–Forchheimer block, `-dp/dx_a = mu u_a / kappa_a +
  1/2 rho C_a |u_a| u_a` on the superficial velocity.

Resistances enter `a_P` (Picard in `|u|`), hence SIMPLEC's `d`.

| Fixture (plug flow, symmetry sides) | Reference | Measured |
|---|---|---|
| Grille K = 3.9994 (50 % plate) at U = 0.5 | jump `1/2 rho K U^2`, flat elsewhere | within 1e-8 relative |
| Internal fan (10 Pa, 8 m^3/s free delivery) against a K = 6 grille, both ends open | root of `p0 (1 - Q/Qmax) = 1/2 rho K (Q/A)^2` | Q 4.694396385 vs 4.694396386 in 91 sweeps; reversed fan same magnitude |
| 10-cell porous block, kappa 2, C 3 | `(mu U / kappa + rho C U^2 / 2) L` | within 1e-8 relative |
| Refusals | boundary plane, negative K, non-positive permeability, fan on solid faces | structured errors |

No-claims: a fan is a pressure jump on a plane (no swirl, hub blockage or
blade aerodynamics); planar and volume losses are the declared coefficients
(the perforated-plate formula is Idelchik's turbulent thin-plate
correlation, not a measurement of a particular grille).

### LVEL turbulence (`turbulence`, `SimpleConfig::turbulence`)

`Turbulence::Lvel` adds the algebraic LVEL closure (Agonafer, Gan-Li &
Spalding 1996), the electronics-cooling model for fan-driven flows on
coarse voxel meshes: per fluid cell, the wall distance `L` (Euclidean distance
to the nearest solid voxel box or wall-type domain plane; openings,
inlets and symmetry planes do not count) and the local speed relative to that
wall give `Re_L = V L / nu`, which Spalding's unified law of the wall (`k` 0.41, `E`
8.6) inverts for `u+`. Interior diffusion uses the profile's tangent
`nu_eff = nu dy+/du+`; the half-cell flux to a wall uses its secant
`nu_w = nu y+/u+`, so the discrete wall shear `mu_w |u_P - u_wall| / L` equals
`rho u_tau^2` of the law of the wall at any first-cell `y+` (using the
tangent at the wall over-predicted the friction by 61 % and the Nusselt
number by 43 % at Re_Dh 2e4). Both are under-relaxed by one half per sweep.
`FvFlow::eddy_conductivity` gives the energy equation `rho c_p nu_t / Pr_t`
(the secant value in wall-adjacent cells); `fv_natural_convection` applies
the same closure and the positive `SimpleConfig::turbulent_prandtl` (default
0.9). The cooling-cht scene selects it with `solver.turbulence: "lvel"` and
may declare `solver.turbulent_prandtl` independently.

The separable linear-time distance transform minimizes squared distance to
the actual cell intervals, including on graded grids and at obstacle corners.
Selecting the nearest solid centre and then measuring only that box is not
equivalent: a regression with unequal widths previously returned 2.5 m for a
surface 1.5 m away. A direct box/plane oracle covers graded and uniform cells,
all domain-wall selections and unit scaling. The same minimizing feature
retains the selected domain-wall identity and its prescribed velocity;
internal solid voxels remain stationary. A native closure regression checks
that a nearer wide solid retains its stationary frame even when a moving
domain wall's ghost centre is nearer. This advances `CHT_SEMANTICS_VERSION`
to 5. LVEL evaluates
`|u - u_wall| L / nu` in that wall's translating frame, consistent with
[Spalding wall functions](https://cpp.openfoam.org/v13/nutUSpaldingWallFunctionFvPatchScalarField_8C_source.html).
Both tangent momentum viscosity and secant wall/heat transport use this same
relative speed. A uniform fluid/wall translation therefore adds no modeled
turbulent viscosity or conductivity. The original deterministic nearest-feature
tie breaking and distance geometry are retained, and the wall identity is
cached at construction without a per-sweep geometry search. Focused G0/G3
tests in `simple/turbulence_tests.rs` cover all three translation axes, graded
Couette Galilean boosts, mixed moving domain walls and stationary voxel solids,
and a public moving-wall plug-flow solve. This extends the equilibrium closure
to the existing prescribed tangential wall velocities; it supplies no moving
geometry, rotating-frame or turbulence-history model.

The varying effective viscosity enters the complete symmetric momentum stress,
`div(mu_eff (grad(u) + grad(u)^T))`, including normal and reciprocal shear terms.
Both reciprocal components share one positive coefficient at each staggered
shear edge: harmonic transmission for molecular material jumps, and the native
arithmetic interpolation for smooth eddy fields. Shear edges and wall secants
use actual half-cell area weights on graded grids. The isotropic Reynolds
stress is part of the pressure; this algebraic closure does not resolve
turbulent kinetic energy. This is the incompressible Boussinesq stress
convention described by the
[Turbulence Modeling Resource](https://tmbwg.github.io/turbmodels/noteonrunning.html).
G1/G3 tests in `simple/turbulence_tests.rs` check a three-grid manufactured
variable-eddy shear balance and equality of reciprocal stress coefficients.
LVEL admission refuses domains with no no-slip wall or solid: an infinite wall
distance is outside this wall model, not a licence to invent a finite length.

The historical channel verification below used the earlier scalar
eddy-viscosity diffusion operator. It compares against the continuum LVEL model itself (the developed
half-channel ODE `rho nu (1 + nu_t/nu) du/dy = tau_w (1 - y/h)`, solved in
the test), validation against Dean's channel friction and Gnielinski's Nu
(a pipe correlation on `D_h = 4h`, Pr 0.71). Half-channel, wall at y-,
symmetry at y+, `h` = 10 mm, 120 h long, developed values over 80..110 h:

| Cells across, Re_Dh | Darcy f | continuum LVEL | Dean | Nu (x = 100 h) | Gnielinski |
|---|---|---|---|---|---|
| 6, 2e4 (first cell y+ ~ 25; 208 sweeps) | 0.03342 | 0.03395 (-1.6 %) | 0.02920 (+14.4 %) | 46.19 | 51.77 (-10.8 %) |
| 12, 2e4 (`--ignored`; first cell y+ ~ 13; 523 sweeps)* | 0.03493 | 0.03395 (+2.9 %) | 0.02920 (+19.6 %) | 52.49 | 51.77 (+1.4 %) |

\* Measured before the momentum residual took the flow-wide scale (same
1e-6 tolerance; the 6-cell case moved by 0.03 % in f across that change).

The friction excess over Dean is the model's own (+16.3 % at Re_Dh 2e4,
+12.8 % at 1e5 for the continuum LVEL): extending the log-law mixing to the
centreline flattens the profile (centreline/bulk 1.098 against Dean's
1.150). No-claims: equilibrium algebraic model only (no turbulence
transport, history, separation physics or transition prediction; it
applies its eddy viscosity wherever `Re_L` is large, laminar regions
included); the bands above are the evidence, not a general accuracy.

### Unsteady flow and transient conjugate marches (`simple_unsteady`, `march_conjugate`)

The momentum equations gain `rho V du/dt` over each staggered volume,
implicit by backward Euler or BDF2 (`TimeScheme`; BDF2 starts with one
backward-Euler step), inside `a_P` so SIMPLEC's `d` stays consistent. Each
step iterates sweeps to `inner_tolerance` on the step's own mass and
momentum residuals (a step that misses it in `inner_iterations` refuses as
`SolverNotConverged { system: "unsteady flow step" }`), then one tight
pressure correction leaves its fluxes divergence-free. Inlets follow an
`inlet_schedule(t)` (ramps, fan start-up); a face fan's operating-point
iteration is steady and refuses (internal fans march). Records carry the
sweeps, residuals, a probe velocity and the kinetic energy per step, plus a
time-averaged velocity and time-averaged face fluxes (`mean_field`,
divergence-free as an average of projected fields) over a declared window;
the steady energy equation may run on `mean_field` when the flow's time
scales are far below the solids' thermal ones (no-claim: that neglects the
unsteady correlation `<u' T'>`). `march_conjugate` advances
flow and energy together: each step's flow (with the Boussinesq force of
the previous step's temperature when declared, and LVEL's eddy
conductivity), then the backward-Euler energy step of `march_energy` on
those fluxes with its two-way closure.

| Fixture | Reference | Measured |
|---|---|---|
| Smoothly started channel (32 x 8, inlet `1 - exp(-t/0.05)`), probe at t = 0.05 viscous times | 256-step BDF2 reference | errors 1.894e-5 -> 4.100e-6 (32 -> 64 steps): order 2.21 (BDF2); 1.330e-4 -> 6.799e-5: order 0.97 (backward Euler) |
| Same channel, 40 steps of one viscous time | steady SIMPLEC solution | largest velocity difference 2.5e-10 |
| Plug flow heated mid-channel, 10 steps | `march_energy` on the steady plug flow | 8.2e-11 K; every step closes |
| Refusals | face fan, empty march, compact junction without declared capacity | structured errors |

An impulsively started flow (inlet on at t = 0 over fluid at rest) is
non-smooth initial data: BDF2 measured order 1.35 at t = 0.04 there;
smooth starts (`inlet_schedule`) recover second order.

### Graded grids (`VoxelDomain::graded`)

A domain may declare its cell widths per axis (a non-uniform Cartesian
grid: fine around small features and boundary layers, coarse elsewhere, the
electronics-cooling grid). Every finite-volume path generalizes:

- energy: face conductance `A / (w_P/(2 k_P) + R'' + w_N/(2 k_N))` with the
  cell's own face area and half widths; boundary half cells `A k / (w/2)`;
  heat capacity `rho c V_c`; sources by volume (`add_uniform_power`);
- SIMPLEC (MAC): each face velocity's control volume spans the two adjacent
  centres (`(w_- + w_+)/2` long, the face's area across); its transverse
  sides take area-weighted fluxes from the two half cells, diffusion uses
  node-to-node distances (half the row width to a wall), and body forces,
  porous terms and the time derivative weight the two half cells by width;
  fluxes, mass residuals and the pressure correction use each face's area;
- radiation: rays cross non-uniform cell boundaries (each exit time from its
  own face coordinate); escape and exchange use each face's area;
- LVEL wall distance: an exact separable feature transform over the cell
  centres gives each cell's nearest seed, and the distance is to that
  seed's box (exact for axis-aligned walls; at a solid corner the corner
  distance, not the earlier centre distance less half a voxel).

- cell-velocity projection: distance-weighted interpolation to the actual
  face plane, followed by the same conservative pressure correction with
  `A / distance` on each interior/open face. Fixed inlet fluxes and zero
  solid-face fluxes retain their semantics. Focused `graded_cell_projection`
  tests check linear interpolation, removal of a zero-pressure-drop gradient
  along each axis, disconnected closed pockets, and downstream heat balance.

The lattice-Boltzmann paths (`lbm_duct_flow`, `natural_convection`)
need uniform voxels and refuse a graded domain;
`VoxelDomain::dx` is NaN there so any uniform-spacing assumption fails
loudly.

| Fixture | Reference | Measured |
|---|---|---|
| Two-material slab on 8 arbitrary widths (5 to 30 mm) | exact series profile | every cell within 1e-9 K |
| Plane Poiseuille on a smooth wall-clustered map (walls 3x finer), n = 4, 8, 16 | 12 mu U / H^2 | errors 1.682, 0.488, 0.127: order 1.94 |
| 10-cell porous block on graded x | `(mu U / kappa + rho C U^2 / 2) L` | within 1e-8 |
| Wall distance on a clustered column | each centre's distance to the wall | within 1e-14 |
| Exposed plate on graded x/y | `F = 1`, total-area Stefan–Boltzmann law | 5 W radiated to 1e-6; temperature within 0.05 K |

### Orthotropic solids and contact resistance

`SolidMaterial::with_orthotropic([k_x, k_y, k_z])` declares grid-aligned
principal conductivities (a PCB laminate: about 30 W/(m K) in-plane, 0.3
through the board); each face uses the conductivity along its normal
axis. `ThermalSetup::contacts` lists `ContactResistance { materials,
resistance_m2_k_w }` between two distinct solid materials (a thermal
interface material, a bonded or pressed joint), added in series on every
face the two materials share: the face conductance is
`A / (dx/2 / k_P + R'' + dx/2 / k_N)`, exact for piecewise-constant
properties. `EnergySolution::conductivity` reports the per-axis values.

| Fixture | Reference | Measured |
|---|---|---|
| Laminate (30, 30, 0.3) in series with k = 3, along x and along z | exact series profile | every cell within 1e-9 K |
| k 10 / k 2 bar with a 1e-4 m^2K/W joint | exact series profile; temperature step q'' R'' | every cell within 1e-9 K; step exact to 1e-9 |

### Temperature-dependent conductivity (`SolidMaterial::conductivity_table`)

`SolidMaterial::axis_conductivity_tables` and `with_axis_conductivity_tables`
declare independent absolute `kx(T), ky(T), kz(T)` in W/(m K), on the fixed
grid axes. All three tables are nonempty, positive and strictly ordered in
temperature; a one-knot table is constant. Each axis owns its own knots and
constant endpoint extensions. The legacy scalar table retains its common-factor
semantics; declaring both forms refuses, as does combining either with an
independently resolved molecular/source conductivity field. No coordinate
rotation or off-diagonal tensor component is inferred.

The selected axis law supplies the same exact Kirchhoff mean, nonlinear Robin
half cell, two-material/contact common flux, source-contact surface recovery
and compact-component link. Cell reports and native buoyancy endpoint
coefficients evaluate all three laws at the actual local temperature. Initial
frozen coefficients retain their original temperature epoch; later nonlinear
proposals and endpoint probes borrow one candidate temperature field for
conductivity and reconstruction, without a separate per-cell factor buffer.
Existing source ownership, storage history and residual acceptance remain
authoritative.

The `cooling-cht` material JSON accepts
`axis_conductivity_tables: [kx_points, ky_points, kz_points]`, where each points
array contains `[temperature_k, conductivity_w_m_k]` pairs. Authored axis
tables cannot override a source-derived material. Six direct actual-production
numerical tests passed, including independent directional integrals and Robin
heating/cooling; native regression definitions cover graded profiles, contact
jumps, Robin/transient heat closure and invalid/competing laws. Full native/CLI
integration execution remains pending the complete dependency build.

A solid may declare `k(T)` as a piecewise-linear table; every axis
conductivity scales by `k(T) / conductivity_w_m_k` at the cell's
temperature. Within one continuous material with zero contact resistance,
each face instead uses the positive Kirchhoff secant
`mean(k; Ta,Tb) = integral(Ta..Tb, k(T) dT)/(Tb-Ta)` over the centre distance.
Solid Dirichlet half cells use the same integral to the boundary temperature.
For Robin films, the physical surface temperature lies between cell and ambient
and solves `integral(Tsurface..Tcell, k(T) dT) = d*h*(Tsurface-Tambient)`.
The same table traversal locates its piece; a scaled, rationalized quadratic
solves that piece analytically, with no separate iterative tolerance or budget.
The resulting half-cell secant enters the original series-film conductance.
Integration is exact across every linear table piece, including its declared
constant endpoint extensions; near equilibrium it integrates temperature
intervals directly rather than subtracting two large potentials. Table lookup
is logarithmic, followed only by the crossed knots, with no face allocations.
Reciprocal rows use identical coefficients. Cell-reported conductivity remains
the local value. Different-material and finite-contact faces integrate each
bulk half cell separately and solve for their single common heat flux. The
surface temperature jump equals the authored area-specific resistance times
that flux. Positive tables give a monotone scalar equation. Two knot cursors
merge the relevant events in flux order, invert each linear piece analytically,
and use bounded safeguarded Newton iterations only in the final piece pair.
Normalized temperature drops avoid subtracting large conductivity potentials;
constant laws and equal-temperature limits retain a direct series calculation.
No face allocations or new numerical dependencies are introduced. The same law
applies to native solid/fluid walls, including the existing fluid eddy term once.
Resolved molecular fields retain their explicitly frozen coefficients.

`EnergySolution::solid_fluid_links` retains the actual ordinary wall coefficients
as a contiguous boundary-only list. The solid-to-fluid observer sums those links
and the separate compact links instead of reconstructing a different law from
cell-centre values or rescanning every component. Fixed-temperature solid cells
are included. `tests/conjugate/nonlinear_contact.rs` checks different nonlinear
materials with known interface temperatures, graded axes, reversal, contact
jumps, solid/fluid observations and an independent transient storage balance.

Compact-component links integrate the native conductivity law of the outside
solid's half cell in series with the declared junction-to-board or
junction-to-case resistor and any explicit case/board contact. The combined
area-specific resistance is `R_junction * A_total + R_contact''`; the existing
analytic film law resolves the outside surface without a nested bulk iteration.
Each face receives its correct parallel area share, and the collapsed package
adds no bulk resistance. Frozen material fields and fluid links retain
their declared coefficients. Junction powers, fixed-cell reactions and observers
all use the same stored links. Ordinary boundary conditions on collapsed component
sides contribute no heat, matching the assembly's adiabatic side ownership.
`tests/conjugate/nonlinear_compact.rs` checks independently derived heating and
cooling profiles with unequal face areas, branch powers, and the compact-side
heat-balance regression. Physical marches use explicitly declared compact
junction heat capacity; no capacity is inferred from a resistor. The analytic film is the declared stationary
scalar Robin law, not a resolved fluid boundary layer.

The steady and transient Picard/Krylov owner proposes a full temperature update
from the previous pass. With authored k(T) and native minmod convection, one
outer iteration owns both laws: an upwind seed supplies the initial field, then
each proposal freezes conductivity and reconstruction together for one solve.
No nested undamped reconstruction loop can bypass the conductivity safeguard;
constant-Cp-only solves retain their existing reconstruction owner.
A bounded residual line search evaluates the original
assembly and physical accounting on that update segment, rebuilding conductivity,
nonlinear boundary half cells, material/contact/compact links and native limited
convection from each trial's own full temperature field. Prescribed reservoirs,
held flow, source/contact/radiation outer owners and the original previous-time
storage state stay fixed. Each probe performs no Krylov solve or time advance.
The merit is the maximum of the actual row-scaled residual and every free
cell/junction's local heat residual; each heat residual uses its own throughput
and explicit product-roundoff allowance. Opposing defects cannot cancel globally.
A full step is retained when it decreases that merit; otherwise at most forty
halved steps seek sufficient decrease. A small step alone never permits return:
both residuals must meet the original `EnergyConfig::tolerance`. Stagnating inner
solves may receive a tighter algebraic tolerance, with the configured Krylov
budget unchanged. The existing 200-pass bound refuses unresolved endpoints.
Returned conductivity, wall links, junction branches and reservoir heat are
evaluated at the accepted field; the iteration report sums actual Krylov work.

The production step-selector's direct Rust tests reproduce and resolve the
undamped two-cycle of a sharp positive conductivity transition using an
independent Fourier integral, and check bounded no-descent/error propagation.
Native regression definitions in `tests/conjugate/nonlinear_contact.rs` add the
350.5 K steady/implicit storage endpoint, graded composite contact in both
directions on all axes, and heterogeneous-fluid face Peclet numbers. These
integration definitions still require a complete native dependency build.
The Peclet diagnostic uses the same two-half-cell face conductance as the
energy equation, including graded widths and molecular/eddy conductivity.
`tests/conjugate/nonlinear_contact.rs` also refuses steady and transient fields
whose opposing local heat errors cancel globally and pass a loose linear gate.
G1/G3 regressions in `tests/conjugate.rs` check the source-free analytic
Kirchhoff profile at 1/2/10/20 cells, crossed knots on arbitrary widths along
all orthotropic axes under both heat-flow signs, a finite same-material contact
jump, the equal-temperature limit, and a two-cell backward-Euler contrast-decay
solution with independently closed stored energy. The profile acceptance is
1e-7 K; this is numerical model verification, not measured-material validation.
Additional Robin regressions cover single-cell heating and cooling with an
independent quadratic reference, two opposing films, crossed table knots on a
graded slab in both heat directions, and an analytic backward-Euler cooling
step with actual boundary heat and stored energy checked separately.

### DC Joule heating (`solve_electric`, `solve_electrothermal`)

Solids may conduct electricity (`ElectricSetup::conductors`, per
material: `rho(T) = rho_0 (1 + alpha (T - T_ref))`). Electrodes are boxes
of cells treated as ideal equipotentials, held at a voltage or carrying a
declared total current (at least one voltage electrode is the reference).
The potential solves the current balance on every conductor cell with
`G_f = A_f / (rho_P w_P/2 + rho_N w_N/2)` between conductor cells and
`A_f / (rho_P w_P/2)` to an electrode surface (AMG-preconditioned CG, ILU
BiCGStab fallback); fluid, insulators and the domain boundary carry no
current. Each face dissipates `G_f dphi^2`, split by half resistance, so
the Joule heat equals the electrodes' `sum V I` (discrete Tellegen) to
solver precision. `solve_electrothermal` adds the Joule heat to the
sources and iterates potential and temperature. Its accepted electrical
field is re-evaluated at the returned thermal temperature, with that solve
cached between passes. Acceptance requires both a temperature change below
`1e-9 * max(span, 1 K)` and every cell's Joule-source mismatch below
`1e-12 W + 1e-9 * max(old, new local power)`. Source differences retain
watts per cell, independently of heat capacity or other cells' power;
`ElectrothermalReport::joule_source_change_w` reports the largest mismatch.
Temperature-independent conductors retain one electrical and one thermal
solve. After 200 thermal passes the loop refuses as
`SolverNotConverged { system: "electrothermal (thermal runaway?)" }`,
the answer for a current-driven conductor past runaway.
Refusals: no voltage electrode, overlapping or touching electrodes, a
current electrode with no conducting path to a voltage electrode, a
non-positive resistivity at the solved temperature.

`solve_electrothermal_with_contacts` includes explicit `ElectricalContact`
area resistances in series with bulk half-cell resistances. Contact Joule
heat follows each declared surface partition, including heat assigned to
electrode cells, and joins bulk heat and authored power exactly once.
Electrical contacts do not infer or replace thermal contact resistances.
The contact values stay fixed throughout this native-law solve; the existing
source-resolved material coupling owns temperature-dependent contact laws.
Focused G1 regressions in `electric/electrothermal_tests.rs` compare the
three-cell thermal resistance network, an independent voltage-driven cubic,
and a strongly cooled resistor whose small temperature update would otherwise
hide a resolved Joule-source mismatch.

| Fixture | Reference | Measured |
|---|---|---|
| 10 mm x 1 mm^2 bar, rho 1e-6 ohm m, 10 mV, ends held at 300 K by copper electrodes | `I = V A / (rho L)` = 1 A; uniform 1e6 W/m^3; peak `300 + q L^2 / (8 k)` | 1.000000000000001 A; Joule = delivered = 0.0100000000 W; peak 301.250000005 K |
| Same bar current-driven, `alpha` 3.9e-3, `m L / 2 = 0.5`, 10/20/40 cells | `theta = (cos(m (x - L/2)) / cos(m L/2) - 1) / alpha` | peak errors 3.33e-2, 8.32e-3, 2.08e-3 K: order 2.00; 11 passes; delivered = Joule to 1e-9 |
| Same, 3.5 x the current (`m L / 2 = 1.75 > pi / 2`) | no steady state | refused (electrothermal) |

No-claims: DC only (no skin effect, inductance or capacitance); ohmic,
isotropic conductors; ideal electrodes; no inferred contact resistance,
surface heat partition, or contact pressure/dynamics model.

### Two-resistor compact components (`CompactComponent`)

`ThermalSetup::compact_components` holds JEDEC two-resistor models
(JESD15-3), the standard board-level representation of a package from its
datasheet `R_jc` and `R_jb`: the component's box of solid cells blocks flow
like any solid, but its interior is collapsed. One extra unknown per
component, the junction, receives the power and reaches every cell beyond
the case-top face through `R_jc` and every cell beyond the board face
through `R_jb`. Each resistor is split over its face in proportion to area
and placed in series with the outside cell's half width. An explicit
`face_contacts` or `source_face_contacts` entry may add a package/TIM joint:
`G = A_cell / (R A_face + R_contact'' + dx / (2 k))`, with zero added
contact resistance when no joint is declared. The sides are adiabatic; the box
cells report the junction temperature through identity rows. `EnergySolution::junctions` gives the junction temperature and the
case and board heats; `solid_to_fluid_heat_w` counts the case links into
fluid instead of the collapsed cells. Refusals: a box that holds fluid,
overlaps or sits on another component, or puts its case or board face on
the domain boundary; power, sinks or fixed temperatures inside a box; a
physical transient without an explicitly declared positive junction capacity.

For an outward package heat flux `q`, the case/board surface is
`T_package = T_junction - q R A_face`; the outside surface is
`T_outside = T_package - q R_contact''`. Source laws query these surfaces in
their declared order, retain the original source context, and never query a
collapsed alias as bulk material. Fixed and source contacts share the stored
branch links, nonlinear outside half-cell law, conductivity-defect audit and
native storage equations. Contact layers add no heat or capacity. Sides,
internal aliases, two compact components and conflicting face owners refuse.
`tests/compact_contacts.rs` checks independent graded anisotropic series
networks, surface-only nonlinear source support and reversed order, and closed
fixed/source transients through both resolved-field and native marches.

`heat_capacity_j_k: Some(C_j)` adds a constant lumped junction storage law to
the same resistive model. `march_energy` and `march_conjugate` solve
`(C_j/dt + sum G) T_j_new - sum G T_neighbor_new = P_j_new + (C_j/dt) T_j_old`
in the existing implicit energy system, including nonlinear outside bulk
conductivities. Capacities follow the actual unknown layout: ordinary cells,
zero-capacity collapsed aliases, then one capacity per junction. The independent
step audit adds `C_j * (T_j_new - T_j_old)` exactly once. Package voxel size or
subdivision cannot multiply its capacity. A component's initial aliases must
share one temperature; accepted junction temperatures are copied to aliases
exactly, so a frozen-flow march can restart from its returned temperature field.
The caller offsets the schedule clock on restart. Both native marches apply
`power_schedule` to cell and component powers from the authored setup each step.
Artificial steady pseudo-time explicitly omits physical junction storage.

A component may instead bind one `CompactHeatStorageLaw` through
`ThermalSetup::compact_heat_storage`. `TabulatedHeatStorage` admits an authored
2..4096-knot total C(T) table; `PreparedHeatStorage` retains source-qualified
constituents and their explicit fixed amounts. Constant capacity and a bound
law are mutually exclusive. The shared caloric kernel integrates every crossed
piecewise-affine breakpoint directly, with no extrapolation or latent-heat
inference. Steady energy solves validate ownership but do not query transient
storage support.

The common transient endpoint owner refreshes positive enthalpy secants in the
junction rows while retaining the original physical time level. It serves
frozen-flow and coupled native marches, public resolved-field steps, and
source-caloric fluid/solid marches. Radiation shares the same endpoint passes;
source-caloric whole-march budgets count every additional energy solve.
Capacity-derivative damping and per-junction support limits change only the
constitutive trials. A junction at its span boundary may stay there while
coupled neighbors adjust. Only an actual supported solved endpoint can pass.

Each bound junction must satisfy both its exact storage/secant consistency and
`delta H_j = dt * (P_j - Q_case - Q_board)`, independently of global closure.
`CompactStorageControls` declares temperature, absolute/relative energy,
relaxation and iteration limits. The energy tolerance uses locally transferred
heat plus an explicit floating-point resolution allowance; absolute kelvin is
not the relative heat scale. An insufficient scaled linear solve is tightened
within the same finite budget. Constant tables skip unnecessary temperature
fixed-point passes, while retaining their local heat and support checks.
Accepted ledgers always evaluate the original energy integral once per
junction; source-solid subtotals and compact voxel aliases exclude it.

`tests/conjugate/compact_transient.rs` checks an independent pulsed RC recurrence,
unequal-area aliases, closed-system equilibration, backward-Euler refinement,
nonlinear outside half cells, admission, cancellation and split restart.
Its C(T) regressions cover heating/cooling across knots, the public resolved-field
step, sharp supported laws with out-of-span predictors, local heat errors that
cancel globally, and coupled junctions starting at a support boundary.
`tests/material_unsteady_buoyancy/compact_caloric.rs` includes a separate
prepared-material junction alongside original source-caloric fluid/solid storage.
`tests/conjugate/compact_coupled.rs` checks the same storage and power schedule
through the actual coupled flow/energy march. Energy remains backward Euler
even when momentum selects BDF2. Radiating ordinary surroundings retain their
existing owner; an exposed compact alias cannot receive a second ordinary
radiation or wall sink because it is not a resolved package surface.

| Fixture | Reference | Measured |
|---|---|---|
| 4 x 4 x 2-cell component (R_jc 50, R_jb 20 K/W, 1 W) between k = 5 slabs held at 300 K / 310 K | two-path network `R_b = R_jb + 2 dx/(kA)`, `R_t = R_jc + 2 dx/(kA)` | junction temperature and both resistor heats within 1e-9; balance 8e-12 |

No-claims: the two-resistor model's own limits (boundary-condition dependence
of a 2R reduction, no spreading inside the package, adiabatic sides). The
optional single-junction storage is an explicit lumped sensible-energy
reduction, with either constant or bounded temperature-dependent capacity.
It is not a measured dynamic package impedance, a distributed die/case model,
or a latent-heat model.

### Surface radiation to the surroundings (`radiation`)

`solve_energy_radiating` (and `FvBuoyancyConfig::radiation` inside the
natural-convection coupling) adds gray diffuse emission from every solid
voxel face that borders fluid, `q = eps sigma F A (T^4 - T_amb^4)`, to the
surroundings seen through declared domain faces (`RadiationConfig::
surroundings_k`; openings and inlets, typically). The escape factor `F` is
a Monte Carlo estimate: `rays_per_face` cosine-weighted rays per face,
marched voxel by voxel (Amanatides–Woo), with counter-based streams keyed by
`(seed, cell, face, ray)`, so estimates are bit-reproducible and carry the
standard error `sqrt(F (1 - F) / rays)`. The energy equation takes `q` as a
cell sink. Mixed-temperature surroundings are reduced using
`T_amb = <T_surroundings^4>^(1/4)` over escaping rays, so incident radiant
powers add before the effective temperature is formed. Exit counts by
domain face make this reduction at most six terms per emitting face;
normalizing by the hottest reached surroundings avoids intermediate
fourth-power overflow while preserving uniform surroundings exactly.
The sink is
Newton-linearized about the previous iterate
(`ThermalSetup::cell_sinks`, reported in `EnergyBalance::sink_outflow_w`),
starting from the linearization about the surroundings temperature.

For prescribed-flow steady convection, the same accepted endpoint must satisfy
both the relative temperature-change criterion and each cell's nonlinear
physical radiation power versus its solved tangent at `RadiationConfig::tolerance`.
The original radiation owner supplies the heat scale and floating-point
allowance; a large global temperature span cannot waive local physical heat
closure. The existing `max_iterations` bounds the shared loop, and accepted
flow/temperature fields are returned without another solve. `RadiationReport`
retains the maximum cell defect, relative power defect, net physical outflow
correction and corrected energy residual in watts. Fixed-temperature cells'
corrections also change their thermostat heat, so they cancel from that net
correction. The native energy report continues to describe its solved tangent.
The focused `tests/radiative_convection_closure.rs` target exercises escape-only
and surface exchange on prescribed advective flow, budget refusal despite a
passing temperature gate, independent Stefan–Boltzmann heat, and physical
fixed-temperature reservoir accounting in a sealed exchange enclosure.

| Fixture | Reference | Measured |
|---|---|---|
| Slab filling a box floor, five faces surroundings at 300 K, non-conducting air, 5 W | `F = 1` exactly; `T = (P / (eps sigma A) + T_amb^4)^(1/4)` | every escape factor 1.0; 433.501645484 K vs 433.501645484 K; 5.000000000 W radiated; 6 Newton iterations |
| Two opposed 8 x 8 plates four cells apart, sides open | `1 - F_12`, `F_12` the closed-form aligned-rectangle view factor (X = Y = 2: 0.41525) | mean escape 0.58231 vs 0.58475 with 64 x 4096 rays; reruns bit-identical |

Surface-to-surface exchange (`RadiationConfig::surface_exchange`, the
default; `SurfaceExchange`): the emitting faces are grouped into planar
tiles of up to `patch_size x patch_size` faces of one material and
orientation, and every ray ends where it is absorbed, on a patch or in the
surroundings; non-emitting solids and opaque domain faces reflect it
diffusely (re-radiating surfaces, at most `max_bounces` times). The Monte
Carlo exchange `A_i F_ij` is averaged with `A_j F_ji` and then balanced
symmetrically (`G = D S D`, symmetric Sinkhorn) so the factors are
reciprocal AND every row sums to the emitter's own absorbed share:
pairwise averaging alone broke the row sums, and clipping the negative
self-views it produced biased the black-shell fixture below by 1 % in
`T^4`. Patch irradiation solves the gray radiosity system by Gauss–Seidel;
each face loses `eps A (sigma T^4 - H)`, Newton-linearized in its own
temperature with `H` lagged. The net heat the surfaces lose equals what the
surroundings receive (zero in a sealed box). The escape-only model remains
available (`surface_exchange: false`); `fv_natural_convection` uses either.

Transient, steady, and finite-volume natural-convection endpoints
independently check the original irradiation
equations and the net surface-emission versus external-escape balance.
They continue iterating when a small radiosity change alone is insufficient,
poll cancellation each sweep, and refuse unresolved reflection at the
100,000-sweep budget. Extremely small emissivity therefore cannot silently
publish fictitious heat loss from a sealed enclosure. The public
`SurfaceExchange::irradiation` sampling API retains its legacy return type
and change criterion; all physical conjugate radiation solvers use the
checked path and reuse irradiation for each endpoint's sinks and heat loss.

| Fixture (S2S) | Reference | Measured |
|---|---|---|
| Convex 2 x 2 x 2-cell cube (eps 0.8, 1 W) in a sealed 8^3 shell at 300 K, non-conducting air, 1024 rays per face, 48 patches | two-surface enclosure `P = sigma A1 (T1^4 - T2^4) / (1/eps1 + (A1/A2)(1/eps2 - 1))` | black shell: 1.5e-11 relative in T1; eps2 = 0.5: 2.3e-6 (the formula's uniform-radiosity assumption); radiated to surroundings 0 |

No-claims: one irradiation per patch tile; opaque domain faces and
non-emitting solids reflect perfectly (declare radiating enclosure walls as
emissive solid cells); transparent air; gray diffuse opaque surfaces;
Monte Carlo factors.

### Finite-volume natural convection (`fv_natural_convection`)

Steady Boussinesq natural (or mixed) convection on the SIMPLEC flow and the
conjugate energy equation, including OPEN boundaries: `FvBoundary::Outlet`
faces are pressure openings that pass flow both ways (paired with
`ThermalFace::Outflow { backflow_temperature }` for the ambient), with
pressure measured from the reference hydrostatic state. The momentum source
is `-rho beta (T - T_ref) g` on fluid cells (cell-sized staggered volumes
average the two adjacent cells). The expansion coefficient is signed, as in
the transient and material-owned paths: negative expansion reverses thermal
buoyancy, and zero expansion reduces to forced flow or conduction. The focused
cavity regression compares a negative-expansion fluid with reversed gravity
at positive expansion and independently checks the zero-expansion limit.
The coupling starts from conduction
through still fluid with the openings held at ambient (a uniform reference
field carries no buoyancy, and openings without flow anchor no temperature).
Each coupling runs `sweeps_per_coupling`
SIMPLEC iterations under the current force. A batch that meets
`flow.tolerance` is projected before the energy re-solve, so
`temperature_tolerance` times the span bounds the change of the actual
returned projected pair. Acceptance returns that pair without another
pressure correction or thermal solve. This closes the LBM path's
"closed enclosures only" boundary.

When radiation is declared, its positive finite `tolerance` also bounds
temperature change over the span and each cell's physical radiation power
minus the solved tangent, normalized by the original radiation owner's heat
scale and floating-point allowance. A loose buoyancy temperature criterion
cannot waive this physical heat gate. The shared iteration limit is the
smaller of `max_couplings` and `radiation.max_iterations`; an unresolved
radiative endpoint refuses instead of returning a fresh radiation total
beside an unclosed tangent balance. The report retains the relative power
defect, net physical heat correction excluding fixed-temperature
thermostats, and corrected physical energy residual in watts. Focused
regressions exercise a single projected coupling epoch and a heated
radiative enclosure under deliberately loose thermal convergence controls.

| Fixture | Reference | Measured |
|---|---|---|
| De Vahl Davis cavity, 32^2, Ra 1e3 | Nu = 1.118 | 1.1201 (0.2 %), 78 couplings |
| Same, Ra 1e4 | Nu = 2.243 | 2.2671 (1.1 %), 36 couplings |
| Open vertical channel (chimney), isothermal plates 8 cells apart, Ra_b = 20, openings top and bottom | fully developed Elenbaas limit g beta dT b^3 / (12 nu) of the discrete stencil ((n^2 + 2) / n^2 x the continuum) | 0.9795 (L/b = 10), 0.9931 (L/b = 30): approaches from below; energy closure 1e-9; induced in/outflow equal to 1e-12 |

No-claims: steady Boussinesq only (an unsteady configuration refuses
as `FlowNotSteady` rather than returning a time average); radiation only
through `FvBuoyancyConfig::radiation`; the SIMPLEC no-claims apply.

### Worked example (`examples/heatsink_cht.rs`)

`cargo run --release -p fs-lbm --example heatsink_cht -- fv [U] [voxel_mm]`
(finite-volume SIMPLEC flow, the default) or `-- lbm [U] [u_lat]
[auto|bgk|central] [max_steps]`: a ducted aluminium plate-fin heatsink
(30 mm long, five 1 mm fins 10 mm tall at 4 mm pitch on a 2 mm base, 2 W
chip under the fins) in a 60 x 20 x 14 mm duct at 0.5 mm voxels, air at 300 K
and 0.25 m/s. Executed 2026-10-07 on a 4-core host (release, lto off):

| Quantity | Value |
|---|---|
| Lattice | 134 400 cells, tau 0.5191 (BGK), 45 000 steps, 755 s |
| Realized / nominal inflow | 7.306e-5 / 7.000e-5 m^3/s |
| Pressure drop | 0.621 Pa |
| Energy solve | 158 BiCGStab iterations, 2.8 s, balance closure 6.8e-12 |
| Junction temperature | 335.31 K (thermal resistance 17.66 K/W) |
| Outlet bulk temperature | 323.41 K (= inlet + 2 W / (rho c_p Q)) |
| Effective film coefficient | 23.2 W/m^2K over 3.78e-3 m^2 wetted area |

The same case with the SIMPLEC flow, executed 2026-10-07 through
`frankensim cooling-cht examples/cooling-cht/heatsink-duct.json` (debug
build):

| Quantity | FV, 1 mm voxels (16 800 cells) |
|---|---|
| Flow | 76 SIMPLEC iterations, mass residual 2.6e-7, momentum residual 9.4e-7, max cell Re 39, 89 s |
| Inflow / outflow | 7.000e-5 / 7.000e-5 m^3/s (divergence 3e-22) |
| Energy | 61 BiCGStab iterations; 2 W in, 2 W advected out; balance 8e-13 |
| Junction temperature | 334.98 K (17.49 K/W) |

The same scene at 0.5 mm voxels (134 400 cells), `heatsink_cht fv 0.25 0.5
amg`, debug build, executed 2026-10-07: 167 SIMPLEC iterations (mass
residual 3.5e-8, momentum 9.5e-7), 1 765 AMG-CG pressure and 499 momentum
Krylov iterations, flow 3 285 s; energy 145 BiCGStab iterations, balance
2e-11; junction 334.36 K (17.18 K/W), pressure drop 0.468 Pa. (An earlier
ILU-BiCGStab run at this resolution took 4 871 s and then failed the
since-relaxed 1e-12 final projection.)

| Junction temperature | K |
|---|---|
| FV 1 mm | 334.98 |
| FV 0.5 mm | 334.36 |
| LBM 0.5 mm (realized inflow +4.4 %, Ma ~0.28 in the gaps) | 335.31 |

The FV junction moves 0.62 K on halving the voxel and lands 0.95 K (0.3 %
of absolute, 3 % of the rise) below the LBM run at the same resolution: a
cross-method band, not a mesh-converged value. At 1 mm
the fin faces fall exactly on voxel centres; the example and the
`cooling-cht` scene parser now resolve such ties deterministically (the
voxel at a box's min edge is in, the one at its max edge is out), and an
earlier 1 mm run that let rounding decide carried one fin two voxels thick.

Estimated numerical evidence at one resolution. The peak lattice speed in
the fin gaps reached 0.163 (Mach ~0.28): compressibility error is O(Ma^2)
there, so treat the pressure drop as indicative; a lower lattice velocity
(second argument) trades runtime for that error. No experimental
validation is claimed.

### Invariants

- Discrete conservation: the two sides of every interior face carry equal
  and opposite `J` for any temperatures; the balance residual is the sum of
  solver residuals.
- Maximum principle: off-diagonals are non-positive and the operator is
  weakly diagonally dominant for divergence-free fluxes, so with no sources
  no temperature leaves the hull of the boundary data (tested).
- Flow never crosses a non-open thermal face (`FlowThroughClosedFace`).
- A problem with no temperature anchor refuses instead of stalling.

### Conformance (`tests/conjugate.rs`)

| Fixture | Reference | Measured |
|---|---|---|
| Two-material slab, `k` 1 and 10 | exact series profile | every cell within 1e-9 K |
| Parallel plates, isothermal, Pe_Dh 100 | Nu_T = 7.541 | n=8: 7.5015, n=16: 7.5313; error ratio 4.07 (second order) |
| Plates with conducting walls (k_s/k_f = 10), uniform outer flux | Nu_H = 8.235 | 8.2448 (0.12 %) |
| Rectangular duct, aspect 0.5, isothermal, x* = 0.125 | Nu_T = 3.391 (Shah & London) | 3.4326 on 12 x 24 (1.2 %) |
| LBM duct (8 x 16, Re 10, Pr 0.7) vs the analytic-profile solve | same thermal problem | local Nu 3.446 vs 3.503 (1.6 %); profile within 4.3 % of the developed shape; interior layer mass flux uniform to 2e-3 |
| Same LBM duct, developed pressure gradient | Darcy f Re = 62.19 (Shah & London, aspect 0.5) | 2.574e8 vs 2.496e8 Pa/m (3.1 %) |
| Metamorphic: wall `k` x10 under a hot spot | peak must fall | 309.70 K -> 304.40 K |
| Refusals | closed-face flow, no anchor, unknown material, lattice resolution, non-tile dims, cancellation | structured errors |

Release lane (`--ignored`, executed): simultaneously developing duct at
Pr = 0.72, Re ~ 49.5, BGK, against the `fs-convection` Shah–London Table 52
card (`convection.rectangular-duct-laminar-cwt-developing-pr072`):

| Rung | Gz | Nu_m (LBM + FV) | Card | Discrepancy |
|---|---|---|---|---|
| 12 x 24 x 40 | 14.165 | 4.6521 | 4.4457 | +4.64 % |
| 16 x 32 x 56 | 13.590 | 4.6060 | 4.4118 | +4.40 % |

The measured band (positive, shrinking with refinement) is the Level-B
cross-code evidence for that card at this Péclet number (36), where axial
conduction toward the Dirichlet inlet — absent from the boundary-layer
table — raises the entrance heat transfer. Gated at 10 %. At Re = 100 with
the central-moment operator the same duct diverged; it now refuses as
`FlowDiverged` instead of panicking.

### No-claim boundaries (conjugate)

- Constant-property convection: steady or unsteady (`march_conjugate`)
  forced convection
  through `lbm_duct_flow` or `simple_flow`, Boussinesq natural convection in
  closed enclosures through `natural_convection` and in closed or open
  domains through `fv_natural_convection`; turbulence only through the
  algebraic LVEL closure of the FV path; no open-boundary LBM natural
  convection or temperature-dependent properties;
  radiation only between gray diffuse surface patches and to the
  surroundings (no participating media). A natural-convection run that does not settle (for example
  above the transition Rayleigh number) refuses as `FlowNotSteady`.
- Staircase voxel geometry at the declared `dx`; one run makes no
  mesh-convergence claim. Power-law convection is first order where the
  cell Péclet number exceeds about 10. All results are Estimated numerical
  evidence; no enclosure or Verified colour is produced.
- `lbm_duct_flow` runs `BoundaryGrid3` pooled but scalar (no SIMD path);
  steps to steady scale with the slowest viscous mode. No throughput claim.
- `simple_flow` is sequential (assembly, AMG-CG or ILU(0)-BiCGStab), its
  only turbulence closure is the optional algebraic LVEL model (a laminar
  solution above transition is a laminar idealization, not a prediction),
  and power-law momentum convection is first order at high cell Péclet
  numbers; its numerical diffusion is not bounded here. No throughput
  claim.
- No product (`.fsim`) stage consumes this module yet; Journey A budgets are
  unchanged by it.

