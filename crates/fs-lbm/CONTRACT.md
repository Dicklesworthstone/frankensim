## Acceleration endpoints during native sparse rigid-wall motion

`FreeSurface3::planar_rigid_wall_motion_with_body_acceleration_candidate` and
`advance_body_acceleration_with_planar_rigid_wall_motion` compose uniform/affine
body-acceleration endpoints with prescribed tangential translation and axial
spin of the existing completely flooded fixed plane. One native collision,
stream, gas reconstruction and conversion uses the accepted old force and
`F_next = rho_next * a_next(x)`. No zero-time reframe removes the changing load's
half impulse. Newly wet cells retain their own density and position. Actual
reflected populations still own wall impulse, angular impulse and mechanical
work; body acceleration is not added to those boundary receipts.

Every acceleration-owned native step now admits positive finite final wet
populations and actual force-corrected Mach after conversion, including cells
whose raw velocity alone would pass. All native/source/boundary/cancellation
refusals preserve the accepted state. Equal endpoints retain held-field bits.
Focused regressions cover heterogeneous transported density, endpoint impulse,
rigid-wall work, conversion, and refused high-Mach candidate replay. Motion is
prescribed over the interval; this fixed-plane isothermal path adds no normal
geometry motion, finite-inertia wall dynamics or fluid heat/energy equation.

## Actual-density acceleration through native open reservoirs

`ForcedOpenBoundaryGrid3::try_new_with_body_acceleration_preserving_velocity`
attaches spatial lattice acceleration without changing physical velocity or
stress. Serial and pooled `advance_*_with_body_acceleration` accept a borrowed
next-epoch field. Collision retains `F_n`; each streamed donor and reconstructed
endpoint uses its own actual `F_(n+1)=rho_(n+1)*a_(n+1)`. At wall/open seams,
retained wall rows and the native half-source give a scalar affine density
closure, solved directly and checked against the final population mass.
Pressure/velocity profiles and D3Q19/D3Q27 source families retain their native meaning.

A reusable force vector keeps the authored acceleration immutable while other
cells still query their donors. The original two executor passes remain; each
stream tile resolves its endpoint force and physical Mach before publication.
The forced-open owner also exposes pooled held and evolving force-density steps.
Accepted populations, source history, step count and native wall receipt publish
together. Invalid fields, source closure, endpoint Mach or cancellation preserve
all accepted state. Body impulse is trapezoidal; reservoir mass/linear/angular
transport and prescribed-wall impulse/work remain distinct native contributions.

Focused regressions use heterogeneous density, all three reservoir axes,
profiled pressure/velocity, moving-wall seams/corners and native collision families;
independently converged prescribed-force steps check the direct closure.
Mass, linear/angular momentum, wall work, serial/pool identity, zero-source
bit preservation and refusal/retry are exercised. This is the existing fixed
geometry, low-Mach isothermal model; it adds no moving topology, implicit
velocity-dependent force law, impact acceleration or fluid energy equation.

## Scheduled differential pressure in native transient flow

`simple_unsteady_with_pressure_schedule` and its explicit SGS companion sample
one finite signed factor at each native time-step endpoint. Every authored
Pressure/Outlet value, including spatial profiles, becomes
`p_ref + factor * (p_authored - p_ref)`, using the existing fixed reference of
its own connected fluid component. Factor one preserves the original march;
zero removes reservoir head and negative values reverse it. This is scalar
modulation of the differential profile, not independent per-port pressure
control or thermodynamic evolution of a common pressure datum.

The original solver retains its BE/BDF2 velocity history, face inertia,
pressure projection, physical half-cell opening masses and backflow rules.
Internal fan curves, body loads and tangential reservoir velocities keep their
original owners. No per-step allocation or topology traversal is added:
changed samples only admit the original exposed pressure planes. Invalid or
unrepresentable samples and cancellation return no partial march and leave
caller-owned data unchanged. Discontinuous schedules require time steps aligned
with their authored transitions. This flow-only entry does not advance heat.

Four native regression definitions exercise all-axis pressure-driven inertia
and reversal with BE/BDF2, oscillatory temporal refinement and conservation,
spatial profiles with independent large component offsets and constant-one
parity, and rejected samples/cancellation. The changed Rust sources passed the
pinned compiler toolchain's syntax/format parser. Native execution remains
pending: the shared filesystem is full and no cached LBM test binary or complete
dependency outputs are available. No executed convergence claim is made here.

## Independent pressure datums in disconnected finite-volume flow

The native SIMPLE owner centers pressure independently in each face-connected
fluid component. Its first exposed pressure reservoir, in native face/plane
order, selects the numerical datum; a component reaching a face fan retains
that fan's ambient zero. Solids and unrelated components cannot choose its
datum. Pressure differences drive the original momentum and pressure-correction
equations. Reported cell pressure restores each component's authored reference,
while pressure-work heat uses the same local reference as momentum.

A large constant offset in a separate duct therefore cannot launch a spurious
startup flow, erase a small driving head, or overflow pressure admission across
regions that exchange no fluid. Connected/common-datum cases retain scalar
reference storage. Connectivity is resolved at initialization; no new
per-step graph traversal, pressure solve or matrix is introduced. Native
momentum checkpoints retain the reference map and rebase their pressure guess
when a changed phase mask joins or splits components.

Five focused regression definitions cover all-axis Poiseuille flow with
independent offsets, exact rest at unequal component pressures, pressure-work
and physical-heat invariance, checkpoint restoration/topology changes, and
fan-ambient ownership. Execution of these new regressions is pending because
the development executor became unavailable. These are incompressible
mechanical-pressure coordinates, not a gas EOS or a new pressure boundary law.
Result-changing conjugate semantics advance to version 10.

## Source-owned nonisothermal turbulent species mixing

`FvFlow::diffuse_turbulent_gas_species` connects the original momentum eddy
viscosity to finite `ClosedGasParcel` inventories with an explicit turbulent
Schmidt number. Canonical source species and typed mole-to-mass conversion
supply the fractions; the shared native graded-cell half-resistance supplies
`G`. All faces sample one original snapshot with
`dm_k = dt G (Y_left,k - Y_right,k)`. Exact absent species remain absent until
an actual incoming stream introduces their original source.

The existing material species-network owner transports each donor's partial
enthalpy, combines all streams, and reconstructs each affected fixed-volume
gas once. Its NASA caloric increments cancel arbitrary formation references
before thermal arithmetic. The returned source state supplies the updated
temperature, thermodynamic pressure, heat capacity, acoustic state, viscosity
and conductivity. No mixture-average enthalpy replacement, second EOS or
additional formation-energy heat is introduced.

Native cell IDs, dimensions, fluid mask, exact parcel volumes and accepted
carrier masses are admitted before stream construction. Common species on a
conducting face must retain compatible caloric and transport sources. The
original global withdrawal, local/global species and energy budgets remain
authoritative; source partitions, face zero-net-mass defects and aggregate
carrier-mass drift are also checked. Solid, exterior and zero-diffusivity
faces are impermeable. Refusal/cancellation publishes no partial endpoint.

This is a first-order explicit frozen-state component, retaining the original
256-parcel and 4096-stream limits. Excessive withdrawals refuse; incoming
material cannot replenish a donor within the same interval. Face conductances
and typed compositions are prepared once, and no per-face gas reconstruction
is performed. Time/mesh refinement remains necessary. Molecular diffusion,
turbulent heat flux, reactions and subsequent pressure/momentum evolution
retain their own operators. SIMPLEC gauge pressure is separate from the
returned thermodynamic pressure, and the carrier's old passive enthalpy is
not an alternate caloric state.

Focused regressions exercise independent unequal-temperature NASA integrals
on graded cells, formation-reference invariance, native Smagorinsky coupling,
Schmidt/heat-policy separation, introduction from exact species absence,
impermeable solids, exact null controls and immutable refusal.

## Spatial acceleration on shared D3Q19/D3Q27 transport

`ForcedGrid3::try_new_with_body_acceleration_preserving_velocity` initializes
`F = rho*a` in the native half-source frame without changing physical velocity
or stress. `advance_serial_with_body_acceleration` and its pooled/moving-wall
counterparts accept the next spatial acceleration in lattice units (`a_SI*dt²/dx`).
Collision uses the accepted force; each completed native stream tile resolves
`F_(n+1) = rho_(n+1)*a_(n+1)` using its actual transported, wall-corrected density.
Both the grid and native collision Mach limits apply to that physical endpoint.

The caller's vector becomes force storage in place. The original collision and
stream remain the only executor passes; workers perform endpoint admission and
retain canonical error ordering. Native wall receipts measure the same corrected
populations and retain their impulse/work meaning. Failure or cancellation
publishes neither populations, force history nor step count. Existing force-density
steps retain their meaning; a held acceleration requires explicit next endpoints.
The selected molecular/SGS collision field, fixed wall topology and low-Mach
isothermal scope remain native. These APIs do not add open or moving-geometry
boundary physics, implicit velocity-dependent loads, or a fluid energy equation.

## Plane-wall Mann turbulence at native D3Q19 inlets

With `turbulent-inflow`, `BoundaryMannInflow3::new_with_plane_wall` binds the
existing harmonic Mann slip-plane observation to the native D3Q19 velocity
profile. Point samples and exact face averages consume the complete corrected
field, including its evanescent wall term. Frozen snapshots, accepted OU
endpoints and exact accepted interval means retain the atmosphere owner's
modes, stochastic path and clock; applying a profile never advances the source.

The corrected field is divergence-free and has zero wall-normal velocity at
the authored plane. Mean flow and frozen translation must be tangent to that
plane. A mean-normal inlet-flux correction is admitted only when its direction
is also tangent. Every complete sampling footprint must remain above the
plane. Failed geometry, source matching, work/storage admission or cancellation
preserves the previously accepted profile.

Only one complete spatial filter is retained: point observations allocate no
gain array; finite faces store the original filter's required per-mode gains.
The same accepted profile feeds the original ordinary, pooled and SGS boundary
stepping. `spatial_filter()` now returns `Option<&MannBoxFilter>`; wall
observations instead expose `wall_spatial_filter()` for their complete
two-position covariance operator. A homogeneous sinc filter alone does not
represent wall-blocked covariance.

Focused regressions cover wall no-penetration and divergence through native
SGS stepping, independent spatial quadrature of finite faces, original OU
endpoint/interval clocks, and atomic refusal. This is outer-flow slip-plane
blocking, not a viscous no-slip layer or a wall-generated spectrum. Spatial
and temporal averages apply to velocity, not nonlinear momentum/energy fluxes.
The native low-Mach boundary and the authored turbulence model retain their
original resolution and applicability limits.

## Spatial body acceleration on native dense D3Q19

`BoundaryGrid3::rebind_body_acceleration_field` selects a prescribed local
acceleration in lattice units while retaining physical density, velocity and
force-corrected stress. `step_with_endpoint_body_acceleration_field` and its
pooled/selected-wall counterparts advance the original collision and stream
with `F_n`, then use each cell's transported density in `F_(n+1) = rho_(n+1)*a`.
The sampled local value follows that cell through donor observation, regularized
velocity/pressure reconstruction and the existing scalar wall-seam closure.
Body impulse is the trapezoidal sum of endpoint forces; actual wall/reservoir
transport remains separate.

Callbacks sample active cells once in x-fast order before mutation and skip
static solids. One acceleration vector becomes accepted force storage in place;
native population scratch is reused, with no extra collision or cloned grid.
Failed sampling, numerical admission or cancellation publishes neither source
nor populations. The prescribed field replaces both previous force components;
an ordinary step after rebinding holds force density and does not automatically
refresh acceleration. Velocity-dependent loads require their own implicit law.

G0/G1/G3/G4/G5 tests cover uniform-field parity, analytic transported-density
response, stress-preserving rebinding, local sampling, linear/angular wall
impulse, all six profiled reservoir orientations against independently
converged prescribed-force steps, serial/pooled identity, refusal and exact
retry. This extends the existing low-Mach lattice model; it does not establish
general compressible-flow validity or a whole-solver convergence certificate.

## Physical momentum at forced D3Q19 initialization

`BoundaryGrid3` constructors initialize unit-density fluid at physical rest in
the held Guo source frame, `f = f_eq(rho, 0) - S(0, F)/2`. Density perturbation
uses each cell's actual uniform-plus-spatial force. Thus initialization adds
neither a half-step body impulse nor the artificial `F F/(4 rho)` pressure
produced by shifting the equilibrium velocity. The native collision and stream
then supply their ordinary full-step momentum exchange. Boundary semantics
version 2 records this correction; zero-force population arithmetic and the
unforced boundary golden remain unchanged. `set_force_field` still explicitly
selects source coordinates; `rebind_force_field` preserves an evolved physical
state.

G0/G1/G3 regressions cover both force signs, all five D3Q19 collision families,
the initial isotropic pressure tensor, and integrated body impulse on periodic
grids with nonuniform density and spatial forcing. These establish the stated
initialization and balance laws, not a new whole-solver convergence claim.

The public tiled `Duct` constructor and its seeded density perturbation now
use the same physical-rest half-source convention. Both force signs and
nonuniform density therefore start with zero physical velocity; an interior
uniform-fluid probe gains one full force-density-over-density increment per
native step before wall signals arrive. Invalid perturbation amplitudes refuse
before mutation. Zero-force arithmetic, collision, streaming and wall rules
are unchanged. The historical core golden remains a kernel sentinel with its
original explicitly authored raw initial populations in a native unit test;
its hash is unchanged and does not certify the new initialization behavior.

## Positive Kelvin in physical conjugate heat transfer

The scalar and axis-aligned CHT entry points require positive finite Kelvin
for prescribed boundaries, ambient and backflow reservoirs, fixed cells, sinks
and initial storage states. Steady solves, transient steps and converged native
buoyancy results also admit their complete physical temperature fields and
compact junction temperatures before returning an accepted result. Cooling
that reaches zero or negative Kelvin refuses instead of publishing a clamped
or unphysical endpoint. CHT semantics version 9 records this admission.

Signed temperature increments and intermediate native k(T) linear proposals
remain numerical work. Their final physical endpoint must satisfy the Kelvin
domain; an intermediate proposal alone does not establish failure. Native
buoyancy pseudo-time and coupling candidates follow the same distinction.
The native lattice passive-scalar API retains its own scalar contract.
Existing normalized CHT benchmarks now declare an explicit positive Kelvin
datum and subtract it only from their analytic reference observations.

Focused regression definitions cover boundary and storage admission, steady
and transient overcooling, valid cooling with the original heat balance, and a
nonlinear conductive root whose first frozen-coefficient proposal is negative.

## Centered temperatures in conjugate conduction

The native scalar/axis-aligned finite-volume energy solve uses temperature
increments about a reference selected from held bath, prescribed-cell, sink or
old-time data. Boundary, sink and storage loads use temperature differences;
material conductivity, contacts, radiation and caloric queries keep physical
Kelvin. The original sparse operator and iteration budget remain in use.
Native constant-Cp advection retains its zero-K enthalpy convention through the
actual local volume-divergence term; source-owned affine enthalpy evaluates its
original tangent at the reference without replacing its caloric datum.

Physical face accounting centers the same face equation at its actual cell
temperature. Pure conduction therefore vanishes exactly at equal temperatures.
Nonlinear endpoint algebra checks subtract only the actual row-propagated
uncertainty from Kelvin rounding and coordinate subtraction, then apply the
requested increment residual tolerance. Independent physical heat checks remain
active; an absolute Kelvin-sized residual denominator cannot hide a small load.
These result-changing semantics are version 8 of `CHT_SEMANTICS_VERSION`.

Six focused regression definitions cover small heat loads across conductivity
contrasts and Kelvin datums, exact uniform equilibrium, finite-contact heat and
fixed-cell reactions, closed-body transient storage, and nonzero-divergence
open-flow accounting, including backflow-only reservoirs. Pinned Rust syntax
checks and five native tests of the original face kernels and reference selection
passed, including reciprocal affine-enthalpy transport. Full integration execution
remains pending after the dependency build exhausted shared resources; the native
kernel checks do not establish full solver, convergence or performance claims.

## Conducting walls in native open-flow heat transport

`ThermalOpenField3::set_wall_conduction(axis, links, controls)` attaches the existing implicit finite-wall graph to an original lower/upper external wall pair. Multiple axis pairs can coexist, including both transverse pairs of a rectangular duct. Link indices omit the selected normal coordinate and retain the original x-fast surface order. Every endpoint must face an active fluid cell; blocked wall-array slots cannot become hidden solid storage or consume sources.

Each graph prepares its reflections from actual post-collision D3Q7 populations before streaming. The original trial wall cells receive the accepted endpoints once; final joint publication retains the native flow/heat transaction. Rejected attempts alter only reusable graph scratch, which is rebuilt on retry. Internal conductive transfers enter the solid first-law account with opposite endpoint signs and do not become external heat. Prepared factors and workspace are reused; each pair has separate explicit graph resource bounds in addition to the carrier storage limit.

This enables forced throughflow with laterally conducting finite boundaries using one wall-temperature owner. Cross-axis graph connections, resolved solid depth, variable caloric properties, automatic buoyancy feedback and the full developing-channel correlation benchmark remain outside this implementation.

## Native open-flow temperature transport

`OpenBoundaryGrid3::prepare_thermal_field` binds a bounded D3Q7 scalar carrier to the original native D3Q19/D3Q27 geometry, solid mask, physical velocity and accepted time. `step_serial_with_thermal` stages the original momentum step and thermal transport with one final joint publication point. Sources use the existing equilibrium half-kicks and retained advective-flux history; optional per-cell thermal relaxation is admitted for the step. Persistent candidate arrays and histories are reused.

Every velocity/pressure face requires explicit reservoir temperature, including nominal outlets. Accepted local force-corrected normal velocity selects Danckwerts total-flux inflow or zero-conductive-flux outflow independently at each boundary link. Receipts measure actual incoming-minus-outgoing scalar populations. Stationary external walls retain independent finite cells per active fluid link; unassigned internal obstacle laws remain stateless and do not create hidden solid storage. Temperature, velocity, source, resource, cancellation and independently measured fluid/solid balance gates precede publication of populations, wall state, histories, receipts and counters.

The conserved scalar uses fixed reference volumetric heat capacity. This is passive low-Mach temperature transport, not compressible mass-specific enthalpy or automatic buoyancy feedback. Its heat receipt must not be replaced by the momentum reservoir mass receipt times temperature. The scalar component velocity envelope is additional to native Mach admission. No nonreflecting-outlet or general coupled time-order claim is made.

### Prescribed spatial and evolving force in native convection

`ForcedOpenBoundaryGrid3::prepare_thermal_field`, `step_serial_with_thermal`
and `advance_serial_with_thermal` carry the same original D3Q7 thermal owner
alongside native spatial force epochs. Initial thermal equilibrium and heat
advection use each cell's accepted physical velocity, including its local
half-force correction. Collision retains the opening force; reservoir
reconstruction and the final source half-kick use the actual endpoint force.
The tensor-diffusion variants retain the original full SPD diffusion kernel.

The final publication advances momentum, force, temperature, finite walls,
flux history, exterior-wall momentum/work receipts and counters together.
An accepted changing-force step exchanges the caller's endpoint buffer with
the old force buffer for reuse; refusal or cancellation preserves accepted
states, the last wall receipt and the caller's buffer. There is no per-step
adapter allocation, extra force reframe, inferred body-work heat or duplicate
thermal source. Reference heat capacity and the original opening-velocity
thermal split remain unchanged; no new coupled time-order claim is made.
Focused G1/G3 regressions cover analytic affine acceleration and source moments
on both native lattices, spatial heat-packet displacement, independent combined
fluid/wall heat closure, late refusal with exact retry, equivalent tensor and
held-force representations, and moving-wall receipt publication/cancellation.

## Finite immersed bodies in native convection

`ThermalOpenField3::set_immersed_thermal_bodies` groups explicitly selected
original solid voxels into finite lumped thermal bodies. Each body has ONE
temperature, heat capacity and held external source across every exposed
D3Q7 link. The original native flow mask and no-slip obstacle stay unchanged.
The same capability is available through `BoussinesqOpenGrid3`, retaining
native buoyancy/heat publication and the separate fluid density envelope.

Each body's backward-Euler balance sums the original reflected-pair resistance
law over actual post-collision outgoing populations. All faces use the common
endpoint; neither capacity nor source is multiplied by voxel/face count.
Per-body receipts measure actual reflected heat and independently computed
storage. Internal-wall and combined fluid/solid receipts include each body
once. Refusal or cancellation preserves body temperatures, accepted receipts,
thermal populations/history and native flow. Interface lists are prepared in
streaming order, with explicit additional storage bounds and no step allocation
or geometry lookup.

The independent-body installer uses an isothermal-interior approximation with fixed heat capacity and a
common authored bulk-to-surface resistance per exposed link. Its applicability
requires negligible internal temperature gradients; it does not resolve solid
diffusion, body contact, moving interfaces, variable caloric properties or
surface-area error of stair-stepped native voxels. Focused regressions check
the analytic ten-face/two-voxel step, capacity-weighted closed-box equilibrium,
separate hot/cold body inventories, periodic seams, exact refused-step retry,
and ownership, storage-budget and temperature-domain admission.

## Spatial immersed-solid conduction in native convection

`set_immersed_thermal_network` installs original immersed thermal nodes and
authored reciprocal conductances together. Different nodes retain distinct
temperatures, capacities and sources; a fully buried node may have zero fluid
interfaces. Each positive-conductance component must reach an exposed node.
The same installation is available through `BoussinesqOpenGrid3` without
detaching its force history or thermal publication.

One simultaneous backward-Euler solve couples all node storage, external
sources, solid conduction and actual post-collision fluid traces. The existing
`fs-sparse` LDLT implementation factors the diagonally scaled operator once;
later steps reuse the factor and retained linear workspace. No temperature
copy becomes a separate physical owner. Each graph edge supplies the same
endpoint heat with opposite signs to its two original nodes. Actual reflected
fluid heat, original capacity times representable temperature change, source
and conductive heat must close per node and over the complete graph. The
existing final thermal/flow publication is the only acceptance point.

Node/interface and graph/factor storage have separate explicit budgets.
Repeated undirected paths must be combined before installation; invalid
graphs, exhausted factor budgets, unsupported temperatures, numerical balance
failure and cancellation preserve the previous accepted state and receipts.
Geometry and fixed material/contact conductances remain caller-owned. This is
spatial lumped-cell conduction with first-order solid time accuracy; it does
not introduce nonlinear caloric laws, moving geometry, an automatic accuracy
certificate or a new stability guarantee for the fluid scheme. Focused tests
cover independent two-node balances including stiff coupling, heat release
from a fully buried core through an exposed shell, whole-capacity equilibrium,
late cancellation/exact retry and graph admission.

### Cartesian material voxels on the immersed network

`set_immersed_thermal_voxels` prepares the same native thermal network directly
from selected original solid cells. Each voxel retains its own fixed capacity,
temperature, isotropic conductivity and external source. Neighboring selected
voxels share the harmonic conductance of their two half cells; the fluid-facing
resistance adds the solid half-cell path to an explicitly authored contact
resistance. Interior solid contacts are ideal. No temperature or heat inventory
is maintained by the returned cell-to-node mapping.

The installer uses canonical native cell order, the original periodic axes and
actual Cartesian faces. Distinct periodic faces between the same two nodes add
in parallel; self-neighbors carry no heat. Unselected solid neighbors and
nonperiodic exterior solid faces insulate. The geometry preparation has its own
payload bound in addition to retained node/interface and sparse-factor bounds.
All installation checks precede replacement of the accepted network. The same
entry point is exposed through `BoussinesqOpenGrid3`.

Four executed native tests cover the independent heterogeneous two-cell balance,
half-solid surface resistance, periodic parallel faces and input-order invariance,
a sourced 3x3x3 solid with a buried center, and transactional topology/budget
refusal. Together with the graph, body and original owner tests, 38 native
tests passed against the actual sparse solver. Momentum-wrapper and complete
Cargo integration remain separate verification obligations. Spatial/time
refinement, nonlinear caloric laws, rotated anisotropy and subcell surface
accuracy are not established by these local balances.

The `immersed_body_cooling` example now accepts `resolved` (default) or `lumped`
solid treatment on the same native flow geometry. The resolved case applies
heat to the interior voxel and reports capacity-weighted mean, minimum and
maximum solid temperatures alongside actual fluid heat and combined balance.

## Spatial finite-wall conduction on native thermal slabs

`ThermalLbm::set_wall_conduction` and `ThermalLbm3::set_wall_conduction` prepare authored reciprocal conductances between existing finite wall cells. The original wall values remain the only physical solid temperature/storage owner. The shared sparse LDLT factor is prepared with diagonal scaling and explicit node/link/factor bounds; each subsequent boundary solve reuses retained workspace. Unconnected wall links retain their original scalar law.

A single backward-Euler boundary solve couples solid storage, sources, lateral or explicitly authored cross-plane conduction and the actual post-collision scalar reflection. Each edge is evaluated once with opposite endpoint contributions. Independently measured reflected heat, actual representable endpoint storage and internal conductive transfer must satisfy local and whole-plane first-law checks. Wall replacement invalidates the operator; failed preparation preserves the accepted plan. Checked source snapshots retain the plan and roll back wall state and transport history with the native populations.

This is fixed-capacity, fixed-resistance conduction on lumped wall cells. It does not resolve solid depth, infer geometry/conductivity, establish nonlinear material behavior, or remove the fluid solver's stability/refinement requirements. Time integration is first order for this solid solve.

## Finite solid cells on native thermal walls

`ThermalWall::FiniteSolid` couples an explicitly sized solid cell to each
normal D2Q5/D3Q7 slab link. Capacity, bulk temperature, through-wall/contact
resistance and external heat per original lattice step are authored in the
same scalar/energy coordinates as the fluid. A uniform declaration expands
once into the existing per-link spatial state; live bulk temperatures are
read through `thermal_wall_fields`. Ordinary steps allocate no wall scratch.

The local backward-Euler solve closes both `C dT = Q - q` and
`T_bulk - T_surface = R q` against the original reflected-pair heat. The
fluid receives that heat once. `ThermalWallExchange::solid` retains independently
measured solid storage, its external source, the actual fluid exchange and
the largest per-link first-law defect. Checked source/material steps admit
surface and solid-bulk temperatures and both summed and local balance before
publishing; existing spatial snapshots roll wall state back with flow, scalar
populations and advective history on any endpoint refusal. At extreme capacity
ratios, an unrepresentable solid temperature increment can leave a first-law
defect. Checked source/material steps refuse defects above the caller's balance
tolerance; ordinary unchecked steps expose the defect in their wall report.

Finite-wall kernel tests exercise the analytic interface/storage identities,
heating/cooling, conductivity changes and tiny/large capacity limits.
The shared scalar integration regression evolves both lattices under tangential
advection and spatial heating/cooling, checking fluid/solid heat accounts each
step and the independent capacity-weighted equilibrium after the source stops. Native
`finite_thermal_walls` regressions cover both lattice stencils, capacity-weighted
equilibrium, combined energy, one-time source deposition, independent wall
cells and endpoint rollback; the source-history test also rejects after full
candidate advancement. This is lumped first-order solid storage behind a
stationary planar resistance: no lateral wall conduction, resolved solid
diffusion, temperature-dependent solid law or coupled second-order time claim.



## Midpoint material tensor diffusion

The D2Q5/D3Q7 material slabs expose `step_with_material_tensor_diffusion_midpoint`
for source-bound, rotated K(T). The original native step is solved at the actual
old/new temperature midpoint with bounded under-relaxed iteration. Temperature
consistency and a fresh relative check of every principal conductivity must pass;
a strong axis cannot hide a mismatch on a weak one. The original caloric storage,
physical source/force conversion, wall exchange and energy acceptance remain
authoritative. No directional average replaces the tensor.

Each trial uses the original source-step transaction and accepted populations,
including force frames and advective history. Nonconvergence, source-domain or
endpoint refusal, and cancellation publish no partial time step. Density/Cp remain
constant over their source domains. This closes conductivity lag in the stated
fixed-point equation; it does not claim general second-order time accuracy or
unconditional stability of the split LBM. G1/G3/G4 integration cases compare an
independent nonlinear two-cell tensor reference, exact constant-law replay in
both dimensions, weak-axis admission, and transactional failure/replay.

## Planar fluid force-and-couple impulses

`Grid::apply_wrench_impulse` and `FreeSurface::apply_wrench_impulse` accept a total
linear impulse and out-of-plane angular impulse on strictly increasing original
wet-cell indices. Actual population mass, or retained partial liquid mass at an
interface, determines the centroid and scalar moment of inertia. The resulting
translation plus rotation is the least mass-weighted velocity increment realizing
the wrench. Cell centers are `(x+1/2,y+1/2)` in lattice lengths; native planar
impulse and work are per unit depth. A single cell cannot carry a free couple.
Periodic topology refuses because no unwrapped angular coordinate chart is owned.

The shared population impulse retains the original Guo force frame, intrinsic
density and physical viscous stress. Independent linear/angular momentum and
signed kinetic-work readback must meet separate budgets before publication;
invalid inertia, excessive speed, unsupported cells and cancellation preserve
accepted state. Partial liquid mass, fill and classification are retained exactly.
The implementation uses existing population scratch, one compact stencil and a
scalar inertia solve, without a dense impulse field or additional collision.
No time, heat or reciprocal body reaction is assigned. Six focused integration
tests cover the physical quadratures, partial inertia, stress, origin covariance,
Galilean work, braking, unresolved couples and transactional refusal/zero input.

## Compact borrowed free-surface impulse staging

`FreeSurface3::prepare_impulse_stencil` retains a mutable borrow of the original native owner while sampling a bounded, strictly increasing wet-cell stencil. Native populations, located retained liquid mass, physical Guo-frame velocity and original cell centers feed the shared population impulse kernel. Gas, inactive/static solid, piston, unresolved interface inventory, excessive speed and unrepresentable states refuse. Staging preserves original populations until the final commit; dropping a stage/candidate or final cancellation leaves the original unchanged. Candidate readback exposes each original mass/position and actual endpoint velocity so a finite receiving owner can independently gate reciprocal angular momentum and work before joint publication. No fluid clock, held force, fill/classification, geometry or unlocated carry is changed.

`MaterialFreeSurface3::prepare_native_impulse_stencil` reuses original physical impulse admission and scenario speed, while explicitly returning native-unit staging for a composing owner’s SI transaction. It does not silently reinterpret lattice readback as physical units.

The native stencil's three focused tests execute actual retained populations, force-frame stress, readback and transactional refusal. A further native refinement regression checks retained solid keys outside the physical domain and beyond representable child Morton coordinates: refinement refuses before allocation, preserving source mass, geometry and populations. All six native refinement cases pass after this guard. The prior three SI receipt checks remain unchanged; full streaming/material integration was not executed.

## Temperature-activated molecular wetting

`MolecularKineticLaw::try_new_arrhenius` optionally evolves the equilibrium jump
frequency as `K(T)=K_ref exp[(E_a/k_B)(1/T_ref-1/T)]` inside the original calibrated
temperature domain. Activation energy is explicit joules per molecular event;
the frequency in `MolecularKineticControls` is its value at the declared reference
temperature. This constant-prefactor Arrhenius closure adds no inferred adhesion
measurement, Eyring prefactor, extrapolation or activation-heat source. The
original constructor still selects constant frequency. Reference-temperature and
zero-activation responses preserve that original model exactly.

The same law feeds the existing dense, sparse and enthalpy wall-angle owners,
retaining original interface claims and actual local temperature. Accepted samples
report the local jump frequency separately from static material measurements.
The rate and velocity scales must be finite and positive across both temperature
bounds; thermal domain refusal retains native flow, enthalpy and wetting history.
Per-length molecular loss remains diagnostic, without claiming a resolved line
measure, a global surface-energy balance or heat deposition. G1/G3 unit cases
check independent forward kinetics, low-speed friction and constant-model limits;
the existing thermal integration case also exercises activated rates and rollback.

## Isothermal kinetics from reversible equilibrium adhesion

The dense D2Q9 and sparse D3Q19 material owners can attach an original equilibrium
angle or Young--Dupre work with `with_equilibrium_planar_wetting`, then explicitly
select Cox or molecular kinetics with `with_dynamic_equilibrium_wetting`. The
same frozen liquid tension feeds both adhesion and native pressure. Static
coefficients sampled at fixed external conditions retain their source axes and
receipts; time, loading and rate restrictions cannot be silently frozen. No
advancing/receding measurements or contact-speed query are fabricated from work.

`with_equilibrium_wetting_hysteresis` optionally binds independent static pinning
thresholds from that interface, preserving the equilibrium potential and separate
resistance receipt. The existing full-angle rate domains, geometric speed history
and transactional native publication still apply. Pure equilibrium alone bypasses
the velocity quotient when fill changes on a locally flat wall stencil; adding
kinetics or pinning retains the unresolved-speed refusal. The backward geometric
speed drives the next step's capillary ghosts. Line losses remain diagnostics;
these owners do not infer a line measure, heat deposition or closed energy balance.

Four public integration groups exercise both source origins and rate laws in
both dimensions, independent forward equations and native response, separately
sourced pinning, domain rollback, and the pressure-loaded flat-fill distinction.

## Heated wetting at neutral transverse corners

The height-graph equilibrium-wetting owner now closes the two ends of its
sourced side wall against the explicitly selected neutral Cartesian transverse
walls. The simultaneous conditions `n.n_wall=cos(theta)` and
`n.n_transverse=0` give zero transverse height slope and the existing sourced
normal slope. Graph normals, area, Marangoni projection and public heat patches
share this metric. The original raw/smoothed VOF pressure ghost uses even
reflection at that same neutral wall, including diagonal corner samples.
The internal boundary selection survives thermal angle-field replacement;
ordinary planar wetting retains its original stencil.

This closure requires one sourced stationary wall, no internal solids or piston,
complete monotone liquid/gas columns, one angle through each occupied interface
band, and the original graph-slope envelope. It does not supply an angle for a
second substrate, arbitrary wedges, corner films, corner pinning, capillary heat
or a closed surface-energy balance. Invalid source work, geometry and cancellation
still refuse the complete thermal/momentum candidate. G1 tests cover both
graph orientations on all Cartesian wall pairs, raw affine-plane preservation,
zero artificial Laplace pressure, independent corner traction and an actual
heated Young--Dupre corner response through the public material owner.

## Material-flow enthalpy reconstruction

`MaterialFlowAdvection::with_reconstruction` attaches the original thermal owner's shared `LimitedEnthalpyReconstruction` to solved material face fluxes. Selected contiguous axial cells reconstruct total specific enthalpy, including a declared passive equal-density latent inventory; paired face transfers, source-law inversion, accepted history and the combined conduction/advection timestep restriction remain with `fs-thermal`. Unselected cells and external streams retain donor enthalpy. Empty plans retain the donor-cell arithmetic. The observed timeline retains first-order time integration on stationary control volumes; it does not infer multidimensional junction slopes, phase mobility, density-change flow, pressure work or moving interfaces. The `advection()` accessor exposes the underlying graph and actual transport audit; implicit clients must attach the same `reconstruction()` plan to their own native implicit transport owner.

Pure frozen `MaterialFlowAdvection<_, ()>` owners additionally expose native
`advance_ssprk2_to` and `advance_ssprk2_adaptive_to`, including corresponding
`with_coefficients` entries for repeatable stage sources and link/bath
conductances. These methods retain the solved `rho * face_volume_flux` graph,
original material/density snapshot and selected reconstruction in the returned
accepted wrapper. Original fluid-cell mass, geometry and caloric-law identity
are admitted at stage queries and before publication. The native thermal owner
alone integrates storage and heat/mass accounts. Adaptive comparisons are one
full step versus two fine steps, and publish only the conservative fine state
and its ledgers. Mandatory datum-invariant h error resolves latent storage at
unchanged temperature; optional T error queries the original caloric law.
Native events, CFL/caloric admission, shared work budgets and cancellation
remain authoritative. The smooth-solution local estimate gives no spatial or
certified global accuracy bound. Stateful inner observers and live momentum
feedback are not repeatable coefficient queries and do not expose these methods.
Integration regressions use actual material-driven momentum for analytic
latent-throughflow temporal refinement, reconstruction in both flow directions,
separate supplied/stream heat, map/law admission and transactional refusal.

`MaterialFlowFeedback::with_reconstruction` retains that same shared plan across
live source viscosity, bounded Boussinesq buoyancy, mushy resistance and retained
momentum updates, including both native implicit entries. Each actual source
epoch rebuilds the original mass graph before binding the same reconstruction;
no additional momentum solve or alternate caloric law is introduced. Required
axial interior connections remain present at zero flow, allowing rest and flow
reversal. They count against the original connection budget. Dormant transverse
faces add no artificial junction, while nonzero crossflow remains visible and
refuses an unsupported axial stencil. A missing neighbor at an exact-solid lock
also refuses the candidate; external or cross-solid connectivity is not inferred.

The stationary-skeleton liquid-carrier attachment reconstructs phase-carried
liquid enthalpy, leaving solid latent inventory in original mixture storage.
Its equal-density, pure-substance, source-support and physical Stokes-work gates
still apply. The optional bulk entropy face limiter remains incompatible with
liquid-only carriage. Reconstruction can be selected before the first advance;
accepted feedback preserves it with the original material pins, physical history
and transactional candidate. Empty selection preserves donor arithmetic and CFL.

## Adaptive implicit thermoviscous flow

`MaterialFlowFeedback::advance_implicit_adaptive_to` controls time error in the
existing quasisteady Stokes/material/enthalpy endpoint solve. Each comparison
uses one original backward-Euler interval and two fine intervals; only the fine
thermal states and their actual endpoint flow, source and physical work are
accepted. The original `ViscousHeating` owner integrates both accepted intervals,
retaining wall/pressure/body supply and molecular/Darcy/Forchheimer heat without
turning momentum residuals or temporal error into heat. All coarse, fine,
rejected and finite-difference evaluations share the blueprint's one
`max_solves` budget. Individual one-step reports retain per-interval work counts.

Mandatory enthalpy error scales changes from the comparison's starting state,
so latent inventory remains resolved at fixed temperature and absolute
formation-enthalpy datums do not inflate the allowance. Optional temperature
error queries the original caloric laws. First-order step-doubling uses the
actual representable fine spans and accepts no extrapolated state. The caller's
whole-advance joule allowance bounds cumulative absolute accepted algebraic
residuals and the final compensated, datum-relative storage/heat balance;
native actuator/physical-loss power tolerance remains separate. Reports retain
the final original coupled endpoint and accumulated accepted physical accounts.

This entry retains fixed topology, stateless `()` material feedback, and the
original stationary-heating restrictions. Interior or unprocessed native events
refuse before momentum work; target-time events remain pending. The native
accepted-substeps budget counts both fine intervals. Minimum timestep, trial
budget, source/work/nonlinear refusal and cancellation publish no partial
thermal or flow prefix. The estimate controls local time error on admitted
smooth solutions, not spatial/model error or a certified global bound. Focused
integration cases compare continuous source-viscosity heating, latent mushy
mobility under shifted h datums, and shared speculative-work/transaction gates.

## Hydrostatic rest under material density-dependent gravity

`MaterialFreeSurfaceScaling::free_surface3_hydrostatic` initializes a planar
column using the original material-bound Guo force, collision, streaming and
gas-pressure reconstruction. Gravity points toward one declared Cartesian
supporting wall; zero gravity is also admitted. Complete liquid layers retain
the original stationary solid tiles and at least one layer of gas. The caller
supplies a positive absolute pressure reference and a density-excursion bound.
Capillarity and wetting must be explicitly absent.

`free_surface2_hydrostatic` applies the same original density recurrence to
the dense D2Q9 owner. Its X/Y supporting wall is halfway between the wall row
and first interior row; the declared depth counts only complete interior
liquid layers, with at least one interior gas layer. Boundary interface cells
retain full mass so their volume agrees with the native gas-pressure face.
This explicit constructor does not change the generic classifier's half-filled
interface convention. Observations and mass remain per unit out-of-plane depth.

`FreeSurface::new_hydrostatic_column` applies that same profile to a fresh
native D2Q9 grid with an existing wall/weir mask. Complete occupied columns
run from the floor to one declared upper gas face; vertical/submerged gas
exposure refuses. The input is closed and nonperiodic, with downward held
gravity, zero external force and raw unit-density BGK rest populations.
The constructor retains the mask, initializes the original physical rest
frame and full interface mass, and introduces no per-step work. Its masked
weir regressions share the dense integration seam below. Historical generic
initialization and its half-filled inventory remain unchanged.

For lattice acceleration magnitude `a`, set `h = a/(2*cs^2)`. The nearest wet
layer has `rho = 1/(1-h)`; successive deeper layers multiply by
`(1+h)/(1-h)`. This surface factor follows the native missing-gas-link rule.
A continuum pressure profile evaluated at cell centers alone does not supply
the same discrete rest state. Layer densities are computed once and reused
while the original constructor owns liquid mass, sparse activation and the
physical-zero-velocity Guo frame. Ordinary `step()` retains the same gas datum.

This balances the weakly compressible numerical operator. It is not a liquid
equation of state or an incompressible mass projection: native mass contains
the retained density perturbation. Pressure and submerged traction approach
their incompressible hydrostatic limits as `dt/dx` decreases with the physical
configuration fixed; refining `dx` at fixed `dt/dx` does not remove that
compressibility error. The constructor adds no moving-solid or partially
wetted-body force closure.

G1/G3/G4 regressions in `tests/material_free_surface3/hydrostatic.rs` exercise
all six gravity directions, zero-gravity equality, pressure-datum invariance,
actual reflected-link traction and torque with decreasing compressibility,
mass retention and transactional refusal. The declared density, source-domain,
geometry and cancellation bounds remain with the original owners.

The dense `tests/material_free_surface/gas_pressure2.rs` cases exercise all
four directions including a single wet layer, full-depth mass retention,
ordinary native rest evolution, the `dt/dx` pressure convergence ladder,
zero gravity, density/viscosity similarity, pressure-datum invariance and
invalid/cancelled construction. Both dimensions use one temporary profile
value per wet layer and add no work or retained storage to time stepping.

## Original-law shear-free Carreau film base

`film::carreau::CarreauFilm` resolves a positive prescribed flux per width at
fixed temperature under tangential gravity. It owns the original thermal
Carreau-Yasuda law, source-rate ceiling and accepted `CarreauSlitProfile`.
The lower half of a symmetric stationary-wall slit of gap `2h` supplies the
shear-free film exactly: `tau(y)=rho*g_t*(h-y)`. Native constitutive inversion,
positive quadrature, temperature shifting, source support and viscous heat
remain authoritative. No fitted effective viscosity enters the profile.

A bounded thickness solve brackets against the original zero-shear Nusselt
thickness and supported wall stress. Safeguarded Newton uses the exact
`dq/dh=h*gamma_wall` sensitivity and admits both flux residual and relative
thickness correction. Requested and accepted flux remain distinct. Point
samples retain apparent and differential viscosities at the same local shear;
`mu_t*u''=-rho*g_t` supplies curvature. Surface shear is zero, wall speed is
zero, and original heat per wall area equals gravitational work within the
native quadrature tolerance. Density and surface tension are explicit actual
SI inputs, not inferred from a constitutive curve or fabricated material card.

This fully developed, isothermal, planar base does not itself claim stability,
entry-region accuracy, thermal feedback or a generalized turbulent threshold.
All support/work failures and cancellation refuse without a partial film.
Focused regressions compare the exact constant-source Nusselt limit at shifted
temperature, an independently integrated nonlinear profile, gravity/heat balance,
source limits and cancellation.

## Lagrangian turbulence with source-resolved rheology

`RheologySgs3::Lagrangian` composes the original source-resolved molecular
viscosity with pathline-averaged dynamic SGS on periodic or stationary channel
grids. `ChannelRheologySgs3::Lagrangian` uses the same history with prescribed
moving channel walls. Source shear, turbulence fitting and pathline departures
use the same incoming Guo force frame and physical wall velocities. Material
queries retain SI shear and time; the authored SGS memory time remains in
lattice timesteps. Each step samples molecular transport afresh, without
feeding an earlier augmented viscosity back into the source law.

The original source and filter workspaces own all numerical storage. The
history trial remains unpublished through collision, streaming, endpoint-force
admission and moving-wall work measurement. Only a successful native step
commits its populations, force/clock, source report and turbulence memory;
rejected candidates preserve the elapsed history interval for a retry. The
`Lagrangian` report variant distinguishes this model from instantaneous fits.

Focused G0/G3/G4 regressions compare several evolving periodic/channel steps
with the independently composed native operators, serial/pooled results,
source/force refusal and retry, and thermal Carreau--Yasuda moving channels
in all three orientations. This retains the exploratory closure's topology,
coefficient, molecular-ratio and departure bounds. It does not supply
unresolved non-Newtonian stress correlations, a thermal energy equation or
turbulent-DNS validation. Locally adaptive memory is an explicit option below.

## Physical strain in planar material rheology

Native and source-backed D2Q9 thermal and free-surface rheology use the same physical
`sqrt(2 dev_3(D):dev_3(D))` invariant as `fs-material`, embedding the original
planar tensor with `D_zz = D_xz = D_yz = 0`. Equal in-plane expansion therefore
retains its actual three-dimensional deviator. Guo half-source correction and
independent central-moment trace relaxation precede this embedding; the latter's
planar trace contributes `div(u)^2/3` to the squared material shear rate.

Thermal power separately retains `2 mu dev_3(D):dev_3(D) + zeta div(u)^2`.
Ordinary tied collision keeps its existing implicit `zeta = 2 mu / 3`, so its
total planar loss remains `2 mu D:D`; the material query no longer counts that
bulk term as shear. Independent planar bulk relaxation still requires a physical
bulk-heating owner and is refused by the shear-only heating entry. Native
lattice-unit laws and SI-converted `PhysicalRheology` use this same material
invariant, so a planar flow does not acquire another apparent viscosity when
embedded in D3Q19. This corrects the former native full-strain convention;
incompressible simple shear and Poiseuille retain their original law.
Manufactured original-population tests cover native and source
owners, force correction, trace relaxation, actual viscosity, and heat partition
against the material tensor and independently prescribed thermal source.

## Density-aware coarse/fine interpolation

`RefinedChannel` interpolates spatial and temporal ghost data in the physical
Guo population frame `h=f+S(u,F)/2`, whose first moment is `rho*u`.
Each donor's actual density and force determine its source transform; the
receiving force frame is applied only after the physical interpolation.
Restriction likewise blends rescaled physical populations before rebinding
the coarse force. This preserves the existing mass-weighted velocity and
physical second moments without a force/velocity covariance artifact at
nonuniform density. Unforced transfer arithmetic is unchanged. Existing
rescaling, reflux, timestep order and low-Mach model limits still apply;
this does not introduce a variable-density incompressible equation of state.

## Spatial pressure reservoirs and sealed D2Q9 channels

`core2::PressureProfileX2` borrows one positive finite lattice-density value
per y row at each x face. `Grid::step_pressure_profile_x` and its selected-wall
receipt variant apply these held pressure profiles through the same native
regularized reconstruction as `PressureX2`. No callback runs inside a step and
there is no profile allocation. Uniform profiles reproduce the uniform-pressure
arithmetic. All candidates reuse collision scratch; its O(ny) extra capacity
is retained while the caller-visible extent returns to the collision field.

Pressure faces now admit paired stationary-wall face/interior rows. Non-periodic
y requires complete top and bottom wall rows, and each x face needs at least
one open fluid row. Velocity/pressure inlet mode retains periodic y and fluid
face requirements. A profile-length, topology, nonfinite force, or late donor
refusal leaves grid populations unchanged. Wall populations are never rebuilt
as pressure cells. Wall receipts retain the original bounce-back link impulse;
the subsequent reservoir reconstruction is a separate momentum exchange.

Each fluid face prescribes `p = cs² rho`, copies its streamed first-interior
physical velocity and non-equilibrium stress, and restores the destination's
own `rho*g + fext` half-source. The pressure profile can therefore sustain
transverse hydrostatics with the discrete relation
`rho[y+1]/rho[y] = (1 + g_y/(2 cs²))/(1 - g_y/(2 cs²))`.
This follows directly by equating opposite Guo half-source populations across
a streaming link; the continuum exponential alone is not the discrete balance.
Pressure remains a weakly compressible isothermal lattice quantity. SI loads
require a shared material pressure reference. These regularized reservoirs
reflect acoustics and do not supply a characteristic or thermodynamic boundary.

The formulation retains the equilibrium/non-equilibrium split of Guo, Zheng
and Shi, [Chinese Physics 11 (2002), 366–374](https://doi.org/10.1088/1009-1963/11/4/310),
with the existing second-order Hermite stress projection of Latt and Chopard,
[Mathematics and Computers in Simulation 72 (2006), 165–168](https://doi.org/10.1016/j.matcom.2006.05.017).
Focused `pressure_driven2` regressions cover both signs of hydrostatic gravity,
the plane Poiseuille pressure/viscosity law, local Guo/stress recovery,
uniform-profile replay, native wall receipts, scratch reuse and atomic refusal.

## Original-law nonlinear tensor SDIRK2

`advance_nonlinear_tensor_conduction_sdirk2` advances the existing full-tensor
voxel conduction equations with two physical SDIRK2 stages. Every nonlinear
trial reevaluates the original conductivity interval integrals and endpoint
Newton tensors. Both stages retain the physical starting temperature in
storage; stage two uses the accepted first-stage temperature only as its seed.
Full anisotropy, finite-contact traces, sinks and independently authored
stage source/boundary samples retain their original spatial treatment.

The shared interval ledger integrates actual external heat with the tableau
weights. Internal stage heat is excluded from generated power, and prescribed
reservoir storage is counted once, including changing prescribed masks.
Krylov iterations, Newton corrections and rejected line-search candidates
each have one budget shared across the complete interval. Refusal or
cancellation returns no partial interval and leaves inputs unchanged.

Each cell has a fixed, complete sensible capacity C in J/K, without a second
volume factor. This entry does not supply temperature-dependent caloric
storage or latent heat. Second-order time accuracy requires smooth laws and
correctly sampled smooth forcing with fixed constraints; forcing jumps must
bound intervals. The L-stable tableau does not establish unconditional
nonlinear convergence, arbitrary-step positivity or a discrete maximum
principle for general rotated anisotropy.

## Locally adaptive Lagrangian momentum memory

`sgs3::dynamic::lagrangian::AdaptiveLagrangianMemory` is an explicit opt-in
configuration for the existing signed Germano history. Before its first
accepted fit, `LagrangianDynamicSgs::with_adaptive_memory` and the native
`LagrangianSgsWorkspace3` can author a positive dimensionless coefficient
and a minimum time. The step configuration's original `memory_time` is the
upper bound. Existing owners without this selection retain their fixed-time
arithmetic. Scalar and component owners expose the same momentum selection;
their independently authored heat and species memory times remain unchanged.

The local law is `T = theta Delta (I_LM I_MM / rho^4)^(-1/8)`, based on
Meneveau, Lund and Cabot (1996), DOI `10.1017/S0022112096007379`.
It uses the advected accepted history and the original grid-filter width.
The density normalization removes the squared-density units of these Favre
contractions; this variable-density extension is not a calibrated claim.
Both time bounds and theta remain explicit. Zero, undetermined or nonpositive
signed numerator history selects the upper bound without discarding its sign.
This bounded policy differs from the paper's positivity mechanism.

The first fit remains the original local spatial average. Each later local
time is frozen during the existing first-order pathline/source relaxation.
Logarithmic invariant arithmetic avoids products, fourth powers and premature
underflow; old-history decay enters the shared scale before exponentiation.
No per-cell arrays are added. Accepted minimum/maximum times publish only
when the original turbulence trial commits, including native collision and
material-rheology consumers. Cancelled, dropped or refused trials preserve
accepted histories and times. Native-grid times remain in lattice timesteps.

Six new regressions verify unit similarity, advected-history relaxation,
time bounds, signed and extreme-scale histories, all channel orientations,
transactional retries and independent thermal/species memories. The exact
production SGS module graph passes all 80 tests with the pinned compiler and
actual native math dependencies. The full fs-lbm grid/material adapters were
source-reviewed but not compiled in the constrained environment. This adds no
adaptive scalar memory, backscatter-energy equation, calibrated theta, DNS
qualification or higher-order pathline integration.

## Source-local planar bulk viscosity

The D2Q9 thermal material owner now admits separately authored physical
`mu(T,gamma)` and `zeta_3(T,gamma)` through the existing `BulkViscousMaterial`.
Planar reduction uses `zeta_2=zeta_3+mu/3`, so even constant physical bulk
viscosity requires a local trace relaxation when shear viscosity varies.
The original central collision reads an optional prepared per-cell trace-tau
field; default/native scalar selections retain their original arithmetic.

The bounded material solve uses the original Guo-corrected stress numerators
`A=gamma_dev2*tau_s`, `B=sqrt(2/3)|div u|*tau_trace`, and
`N=hypot(A,B/sqrt(2))`. Its auxiliary relaxation satisfies
`gamma=N/(lambda*dt)` and
`1/lambda=hypot((A/N)/tau_s,(B/(sqrt(2)*N))/tau_trace)`.
Normalized weights avoid squared-component overflow/underflow; zero stress
queries the zero-rate source. The report's relaxation defect is in lambda.
Collision receives the actual shear and trace taus from that same source,
never lambda. Shear-only and trace-only inputs recover their respective taus.

The unchanged thermal transaction converts `mu*gamma²+zeta_3*(div u)²` once,
retains original material receipts, checks the endpoint, and publishes the
applied trace field with the accepted state. Explicit `zeta_3=0` is allowed
in planar flow through the positive `mu/3` trace coefficient. The 3-D trace
collision retains its strict positive-relaxation admission; negative bulk
viscosity and unrepresentable positive losses refuse. This remains explicit,
fixed-density/Cp low-Mach sensible-heat transport, without a total-energy claim.
Manufactured tests cover the physical tensor/heat identity, source-dependent
coefficients, zero bulk, original native collision, and full rollback.

## Lagrangian turbulent heat with original material sources

`thermal::sources::sgs3::lagrangian` composes the shared joint Favre/Germano
momentum and scalar histories with the original D3Q19/D3Q7 source transaction.
The resolved slab retains periodic tangent directions and its actual rigid
halfway walls. Current prescribed-temperature, insulating, flux, film and
spatial wall derivative policies are shared with the instantaneous channel
closure; the native thermal owner alone imposes and measures their wall heat.

Stress and heat keep independent signed contractions, memory times, coefficient
limits and molecular-ratio limits. A zero momentum fit cannot disable a
positive heat fit. The existing optional invariant momentum memory is exposed
without altering scalar memory. Sampling reads one accepted Guo-centered
density/velocity/temperature state before replacement forces. The caller
supplies actual sampling time in lattice timesteps; ThermalLbm3 has no separate
authoritative clock. Elapsed sample time controls bounded pathline departures.

One prepared token borrows only the reusable workspace and retains both
candidate histories plus admitted total relaxation fields. Native source, heat
balance, temperature range, endpoint Mach and cancellation checks precede
publication. Material wrappers additionally retain that token through every
original source epoch and outer publication check. Only successful completion
is followed by infallible joint history commit. Dropped, cancelled and refused
steps preserve both histories, including their last sample time for retry.

Material APIs cover ordinary source-resolved transport, source-density
Boussinesq forcing and endpoint-averaged source-density forcing. The last keeps
all three original material epochs: initial transport, thermal-predictor force
quadrature and the actual accepted endpoint. Original SI source samples, fixed
reference inertia/sensible storage, density/Cp bounds and force conversions
remain intact. Fresh molecular mu/rho0 and k/(rho0 Cp0) own their ratio limits;
augmented rates never feed back into a later molecular source query.

Field arrays are prepared once under an explicit numerical byte budget and
reused by sampling and history publication. Original native snapshots retain
their independent cell budget. The instantaneous and Lagrangian adapters share
the same rate augmentation and shear-power arithmetic. No alternate collision,
wall solver or material law is introduced. This remains explicit, low-Mach
sensible-temperature transport; modeled kinetic transfer is diagnostic rather
than heat deposition, and no total-energy or turbulent-DNS claim is added.

Five new native regressions compare multiple evolving steps against the
independent shared history owner and original transport, including fixed and
adaptive momentum memory, independently positive heat, actual wall policies,
Guo frames, exact zero controls and late refusal/replay. Three material groups
cover all wrappers, fresh source epochs and endpoint buoyancy, then warmed
history rollback after late source refusal and cancellation. All changed Rust
files passed pinned parser checks and bounded source review. These eight new
integration regressions have not been compiled or executed because required
dependency artifacts and build storage are unavailable; earlier shared-kernel
and atmospheric test results do not establish this native integration gate.

## Physical planar liquid wrench impulses

`MaterialFreeSurface2::apply_wrench_impulse` applies an explicitly authored
linear impulse and independent out-of-plane couple to retained liquid cells.
Inputs and receipts use the existing per-unit-depth convention: linear
impulse in kg/s, angular impulse in kg m/s, liquid mass in kg/m, and signed
work in J/m. Positive angular impulse is counterclockwise. No thickness or
extra timestep is inferred. The original source density, cell length and
actual interface fill determine the selected mass, centroid and inertia;
the native minimum-mass-weighted translation plus rotation realizes the wrench.

Conversion transports the requested moment to the representable native
origin and transports actual deposited and observed moments back to the
original authored SI origin. Separate physical linear, angular and work
budgets admit the complete SI receipt before native population publication.
Scenario speed, topology, singular-couple and cancellation gates retain the
original transaction; a late unit or balance refusal also leaves the accepted
liquid unchanged. Held Guo forcing, viscous stress, density, fill and material
identity remain with their original owners. No collision, thermal advance or
reciprocal body reaction is implicit.

Focused regressions cover independent physical impulse/inertia/work sums for
bulk and half-filled cells across distinct length/time scales, covariance
under a shifted SI origin, signed braking work and state-preserving native
and physical refusals. Dependent-receipt rejection is checked at the native
publication boundary.

## Original caloric storage in nonlinear tensor conduction

With `enthalpy-coupling`, `advance_caloric_tensor_conduction` and `advance_caloric_tensor_conduction_sdirk2` add backward Euler and SDIRK2 storage in the original caloric law's specific-energy coordinate. Each cell retains its supplied reference mass and original energy basis/datum. Newton solves centered energy increments with actual storage `m * (e - e0)`; face traces remain Kelvin. The original analytic `dT/de` scales only the cell column of the full-tensor Kirchhoff tangent. Exact zero slopes on latent plateaus remain regular through the positive mass/storage diagonal, without apparent-capacity smoothing, numerical differentiation, or a prescribed spatial phase front.

`TensorCaloricLaw` extends the canonical `ThermalLaw` with original bounded energy support, an analytic temperature derivative, a fused constitutive query, and reference-cell admission. Direct implementations retain the equilibrium phase curve, restricted pure-solid view, and original heat-storage assembly; amount-bound storage enforces its original mass normalization. Every prescribed-temperature cell requires an explicit original specific energy at each stage, validated against its authored temperature, so an isothermal phase fraction is never selected by a temperature inverse. Prescribed boundary-face temperatures do not require a cell-energy prescription.

Both stages reevaluate original caloric and tensor laws, share nonlinear/Krylov/backtracking budgets, and reuse the independently audited interval heat ledger with actual representable stored-energy changes and prescribed-cell reservoir corrections. Conductivity remains the independently supplied tensor law `K(T)`; independent phase-fraction dependence, moving geometry, expansion work, and advection are outside these entries. Backward Euler is first order. SDIRK2 is second order on smooth caloric/forcing segments with fixed constraints; unresolved sharp phase transitions or forcing changes require time refinement or interval splitting. Either stage's refusal or cancellation returns no partial interval.

## Lagrangian momentum and heat in conservative NSF3

`NsfGrid3::with_lagrangian_thermal_sgs` binds the original joint Favre--Germano
pathline histories to native D3Q27 mass/momentum and conservative total energy.
Periodic fluid boxes and resolved channels with two periodic axes are admitted;
actual translating-wall velocities and original thermal-wall derivative policies
feed the same filters. Independent momentum and temperature memories, optional
bounded invariant momentum times, coefficient bounds and molecular-ratio
admissions are explicit. Times are lattice timesteps from the native grid clock.
Source-owned viscosity, heat capacity, conductivity and bulk relaxation retain
their original state and domain. Only deviatoric collision receives modeled
viscosity; shared interior heat faces receive rho*(cp/R)*alpha_t. Mechanical
transfer stays in the original total-energy stress work, without an additional
heating source. Molecular half-cell conduction remains the wall heat law.

Every endpoint corrector/probe starts from the same accepted signed histories.
The final source/force/temperature state retains a borrowing joint history token
through all later native work, energy and cancellation admissions; one infallible
commit accompanies physical publication. Failed/retried steps cannot age memory.
Same-time wall/force changes re-emit the accepted contractions with zero elapsed
transport, without publishing a new model sample. Regular timesteps reuse the
264 nominal bytes/cell filter/history payload and do not copy history or repeat
the final fit. Thermal-wall preview intervals allocate an independent workspace
within their cell and original per-workspace bounds, copying only accepted
histories/scales/policy. A rejected multi-step prefix cannot modify its parent.

G1/G3/G4 regressions cover the one-update endpoint against the original shared
model, conservative closed-box energy, final-checkpoint rollback/replay,
moving-channel same-time refresh, independent multi-step thermal trials,
workspace/geometry/usage and representable-clock refusals. Core regressions
cover zero-time refresh with both fixed and adaptive momentum memory, all
three channel normals, current-strain re-emission and history-copy independence
from rejected scratch. An exactly zero interval is an identity; a positive
interval whose adaptive relaxation underflows still refuses.
The pathline/source update remains first order even though native energy uses
its original endpoint trapezoid. No open boundary, interior-solid filter,
backscatter-energy model, complete compressible LES closure, turbulent-wall
function, shock qualification or DNS claim is added. Native integration tests
are authored but require execution in a complete build environment.

## Original source contacts in nonlinear tensor conduction

The native nonlinear tensor owner consumes the existing ordered
SourceFaceContact assignments directly. This applies to steady conduction,
sensible backward Euler and SDIRK2, and both original-caloric time drivers.
Every finite source contact retains two physical face temperatures. The
original declared interface query seeds these traces; unequal bulk cell
temperatures may lie outside the interface source's temperature box, while
their conductivity intervals must remain supported by their original laws.

Every physical Newton trial refreshes the same selected original resistance
claim at those actual A/B surface temperatures. Card identity, surface order,
pressure, separation, finish, medium and other supplied conditions remain
with the existing source owner. Supported tangent payloads are a declared
constant scalar, an original linear curve, or original affine surface cells.
The analytical chain rule for q = A*(Ta-Tb)/R includes both resistance
temperature derivatives. Source corners use original one-sided primal
tangents and do not establish a unique adjoint derivative.

The two trace equations receive exact opposite contact heat and derivative
rows. Cell increments still condense locally into the existing sparse face
system; the contact adds no storage, external source, boundary heat or outer
fixed-resistance thermal solve. Numerical source responses and derivatives
are cached per trial, and source derivative preparation occurs once per
spatial solve. Original source receipts and actual bulk face heat are returned
in TensorConductionSolution.source_contact_states, with independent local
interface, contact and global energy checks.

TensorNonlinearConfig owns Newton, Krylov and rejected-trial work budgets,
shared across SDIRK2 stages. Scalar SourceContactControls govern only the
separate scalar outer iteration. Duplicate or overlapping contact owners,
unsupported original source states, failed convergence and cancellation
return no partial endpoint or interval. Frozen tensor entries continue to
require fixed contacts; existing fixed-resistance and perfect-interface
arithmetic retain their original paths.

## Original directional sources in tensor conduction

`RotatedSourceConductivity` borrows three original `ResolvedConductionLayer` owners in an explicitly declared proper material frame. The directions must share the original caloric identity, pure phase, temperature-axis convention and all fixed source context. The adapter uses the intersection of their original temperature and energy domains, preserves original source lines and receipts, and refuses extrapolation. Preparation compares fixed context once with bounded cancellation; evaluating a point or interval performs no material selection, table copy or source re-resolution.

Every crossed original conductivity knot contributes to each directional interval mean before rotation. The canonical `ConductivityCurve::mean_conductivity_w_m_k` accumulates positive segment means with convex width weighting, so representable means survive overflow or underflow of the dimensional Kirchhoff integral. Equal endpoints use the original point law; reversed intervals use identical arithmetic.

With `enthalpy-coupling`, the same adapter implements the original caloric storage trait. Its temperature, density and entropy queries retain the narrower source energy support. The binding original source owns each limiting energy endpoint, preserving its admitted one-sided temperature tangent at melting onset. No apparent heat capacity, new caloric identity, phase-mixture conductivity or inferred orientation is introduced. The adapter composes with the existing steady, sensible and caloric tensor drivers and original interface contacts.

Two actual production conductivity-mean tests passed under the pinned Rust compiler, including extreme representability and crossed-knot/support cases. Added adapter regressions cover source intersections, explicit rotation, original-context and chart admission, cancellation, and the pure-solid caloric endpoint; these integration definitions were parsed but not executed in the dependency-constrained workspace.

## Adaptive original-caloric tensor time integration

The original caloric tensor owner admits bounded adaptive SDIRK2(1) windows through `advance_caloric_tensor_conduction_adaptive`. Every completed attempt uses the same two native spatial stages; the embedded first-order energy discrepancy controls dt with the existing square-root thermal controller. Its scale uses changes in original h/u, including latent plateaus, never absolute energy datums. The embedded coordinate is never passed to the material inverse. Only the native second-order endpoint and native accepted stage heat are physical; the discrepancy is a local estimator, not a certified endpoint, global or channel-quadrature bound.

Forcing is sampled at actual native stage times. The callback is deterministic and read-only with respect to physical history and may borrow constant setup/prescription arrays. Strictly increasing declared transitions bound every attempted interval and mark endpoint samples as explicit left limits. Prescribed masks remain constant across all attempts in each smooth segment. Minimum/maximum dt, a maximum attempt count, and the native shared Newton/Krylov/backtracking caps bound work; only temporal rejection and native SolverNotConverged permit retries. Invalid source data, cancellation and other native refusals are fatal.

Rejected state and heat never enter retained physical records. Final stored energy is independently recomputed from initial/final original energies and fixed reference masses. Signed net accepted heat is compensated independently of the displayed source and outflow totals; closure subtracts the compensated storage and net-heat accounts before rounding away either correction. Both net whole-window closure and the sum of absolute accepted native closures must meet the single configured joule allowance. Any late error or cancellation discards the private accepted prefix and leaves caller inputs unchanged.

Five tests of the actual production embedded-error and heat-account arithmetic passed under the pinned Rust compiler. Five added native integration definitions cover diffusion refinement, latent energy and datum shifts, transition left limits and thermostat changes, bounded attempts and minimum dt, and fatal source/mask/cancellation outcomes. These integration definitions and public exports were parsed; full native integration execution was unavailable in the dependency-constrained workspace.

### Full-tensor diffusion on native open 3D flow

`OpenBoundaryGrid3::step_serial_with_thermal_tensor_diffusion` transports heat
with the original prepared `DiffusivityTensor<3>` and borrowed uniform or
per-cell `DiffusivityField`. It uses the existing D3Q7 flux block
`S = (I/2 + alpha_lattice/(1/4))^-1`, retaining off-diagonal transport in all
three Cartesian planes. Physical input scales as `K/(rho*Cp) * dt/dx^2`
with the carrier's fixed, common reference volumetric heat capacity. The
carrier's scalar tau continues to relax even kinetic moments. No tensor is
installed persistently or substituted into later scalar calls.

The original source half-kicks, advective-flux history, native force-corrected
velocities, actual wall/reservoir heat, finite-wall storage and final joint
publication remain shared with scalar transport. Solid grid cells carry no
thermal population inventory; explicitly configured immersed bodies retain
separate shared storage, while unassigned obstacle links retain the original
boundary law. Per-cell field dimensions are admitted before any candidate staging.
Refused or cancelled attempts publish no heat, momentum, wall state or time.
The original scalar collision path retains its arithmetic. G0 tests compare
the isotropic limit for D3Q19/D3Q27 flow with reversed reservoirs and sources,
and exercise late rejection/cancellation plus unchanged scalar retry. G1 tests
resolve signed oblique Fourier decay and spatial refinement in each Cartesian
plane against `exp(-t k^T alpha k)`.

This is anisotropic temperature diffusion, not variable-capacity conjugate
enthalpy or an inferred contact-resistance law. SPD alone establishes neither
population positivity nor full coupled-flow stability. Native wall and outlet
approximations retain their original scope; spatial/time refinement remains
the caller's responsibility.

### Tensor diffusion with native Boussinesq feedback

`BoussinesqOpenGrid3::step_serial_with_tensor_diffusion` carries the same
borrowed full thermal tensor into the original native temperature/force/flow
transaction. The transported temperature, including the final source half-kick,
drives the actual endpoint Boussinesq force. Off-diagonal conduction therefore
changes buoyant motion through the existing force owner and consistency audit.

Accepted and candidate force epochs, trapezoidal body impulse, reservoir heat,
finite-wall storage and wall conduction retain their existing accounting.
Temperature, momentum, walls, force history and time publish jointly; a rejected
or cancelled attempt retains the complete accepted state. The scalar entry
point retains its original arithmetic and no tensor persists into later calls.

Three regression definitions compare the isotropic limit on both native
momentum lattices, verify anisotropy changes endpoint force while closing the
fluid-plus-wall heat balance, and cancel at the final publication checkpoint
before an exact retry. This keeps the existing opening-velocity split and
small-density-departure Boussinesq assumptions. Frozen tensor coefficients and
a common reference heat capacity do not supply variable-storage enthalpy,
total-energy evolution, or a new temporal-order guarantee.

### Binary molecular diffusion on accepted carrier mass

`MaterialAdvectionStep::diffuse_binary_species` advances two original
constituent inventories on the SAME accepted carrier mass. It can consume
`advect_species::<2>`'s endpoint or a previous diffusion endpoint. Prescribed
per-cell molecular diffusivity is held fixed; it is never inferred from
viscosity or from a species index. The mass-frame law is
`j_A = -rho D grad(Y_A)`, `j_B = -j_A`, with density from accepted mass
and native voxel volume.

Every conducting fluid-fluid face uses the two half-cell resistances in
series, including graded Cartesian widths and spatial jumps in `rho D`.
Stagnant faces still diffuse. Exterior, solid and zero-diffusivity interfaces
are impermeable. Uniform mass fractions are stationary even when density
varies. The domain must be the original geometry: dimensions and fluid mask
are checked, while identical physical widths remain the caller's obligation.

Each pair's exponential relaxation approaches its mass-weighted equilibrium.
Forward/reverse half sweeps form a symmetric second-order diffusion splitting
for frozen coefficients. Both original constituent inventories receive
opposite transfers; neither clipping nor complementary-species reconstruction
is used. The less abundant constituent supplies the gradient to preserve
trace transport when majority fractions round to one. Prepared face factors
keep exponentials outside repeated sweeps. Work is
`O(cells + faces*substeps)`; storage is `O(cells + faces)`.

Independent constituent-mass and carrier-partition budgets admit the complete
result. Work/storage limits and bounded cancellation polling precede publication;
carrier, inputs and accepted state remain unchanged on refusal. Eight regression
definitions cover exact unequal-density pair relaxation, uniform fractions,
barriers/solids, tiny traces and species permutation, temporal/spatial convergence,
strong-step positivity/conservation, and invalid-data/budget/cancellation refusal.

This is an isothermal, fixed-carrier binary Fick component. It supplies no EOS,
partial-enthalpy transport, nonideal thermodynamic factor, Soret/barodiffusion or
arbitrary multicomponent closure. Large stable sweeps still require accuracy
refinement. The original donor advection followed by diffusion remains an
operator split with first-order overall accuracy; second order refers only
to the frozen diffusion sweep composition.

`MaterialAdvectionStep::diffuse_species_adaptive::<N>` adds automatic uniform
temporal refinement of that original operator. It compares complete n- and
2n-substep solutions from the same original inventories, accepts the finer
positive conservative state, and reuses a rejected fine state as the next
coarse comparison. Each constituent has its own absolute mass-fraction
tolerance; relative scaling uses the local initial/fine concentration, so
trace species do not inherit a bulk absolute allowance. The per-cell estimate
is `|Y_fine - Y_coarse| / 3`, from the second-order asymptotic regime. This is
not a certified global error bound and does not measure spatial, material or
advection-splitting errors. No Richardson extrapolation alters physical masses.
Aggregate work includes every discarded resolution; storage admission includes
the live coarse comparison alongside the original solver. Refusal or
cancellation returns no partial interval. Independent original constituent and
carrier-partition budgets remain mandatory at every resolution.

### Conservative multicomponent turbulent mixing

`MaterialAdvectionStep::diffuse_species::<N>` extends the same accepted-carrier
diffusion owner to 1..=128 constituents sharing a common diffusivity. Each
original constituent inventory evolves independently under
`j_k = -rho D grad(Y_k)`. One prepared face conductance and exponential factor
serve every constituent. The original binary entry and trace-safe complementary
arithmetic remain available. All components retain separate mass-defect budgets;
the sum is audited against the immutable accepted carrier, with no clipping,
normalization or reconstructed remainder species.

`FvFlow::diffuse_turbulent_species` consumes the actual momentum eddy viscosity
from the original LVEL or transient SGS flow, using an explicit positive
turbulent Schmidt number: `D_t = nu_t / Sc_t`. This is the gradient-diffusion
mass-flux closure in NASA/TM-2022-0008776, equation (46). The thermal Prandtl
number and LVEL heat-wall secant do not select species transport. Source-resolved
molecular viscosity remains with the momentum owner. The original accepted
mass divided by native cell volume supplies density, including graded grids and
spatial density jumps. Coefficients are frozen for the requested mixing interval.

`FvFlow::diffuse_turbulent_species_adaptive` connects the native momentum
eddy field to the original conservative temporal-refinement owner.
`FvFlow::diffuse_molecular_and_turbulent_species_adaptive` additionally consumes
one explicit nonnegative common molecular diffusivity per cell. The law is
`j_k = -rho (D_molecular + nu_t/Sc_t) grad(Y_k)`: molecular transport survives
zero SGS viscosity. The two cellwise contributions are added BEFORE the
original graded half-cell resistance is formed; adding independently
harmonized molecular and turbulent face coefficients is a different interface
model. Molecular D is never inferred from viscosity or the thermal Prandtl
number. Both coefficients must vanish in solids; a zero total coefficient
blocks the face. Negative/nonfinite inputs and unrepresentable positive
transport refuse without clipping or an invented diffusion floor.

Both adaptive entries retain per-species concentration tolerances, all
coarse/fine work in the aggregate budget, and the original mass/partition
accounts. Storage includes the combined D field and the simultaneously live
coarse comparison. The finer conservative state is authoritative; no
Richardson extrapolation alters its inventories. The temporal estimate is
asymptotic and does not certify spatial or model error. The molecular entry
admits binary Fick diffusion or an independently justified COMMON-D molecular
model; it does not stand in for unequal-species or Maxwell–Stefan transport.

The common face map preserves each species' global mass and the cell carrier
partition to the caller's declared budgets. Uniform mass fractions are
stationary at variable density. Zero eddy viscosity is an exact identity for
the turbulent-only entries; molecular diffusion can remain active.
The default entries keep exterior and solid faces impermeable; explicitly
authored surface reservoirs use the boundary entries below. Advective
inlet/backflow composition remains with `advect_species`. Work is
`O(N * (cells + faces * substeps))`, storage `O(N * cells + faces)`;
the checked component-work count and temporary diffusivity buffer are included
in the work/storage account. Cancellation or failed admission publishes no
partial inventory and does not modify the carrier or flow.

G0/G1/G3 regressions cover a four-species analytic pair on graded unequal-density
cells, trace transport and species permutation, the common-diffusivity Neumann
mode's second-order temporal refinement, null controls and resource refusals.
The native transient Smagorinsky-to-three-species test exercises the actual
momentum closure and accepted material flow. Separate coupling checks establish
Schmidt scaling, independence from the heat-wall policy and refusal atomicity.

This is an isothermal common-diffusivity mixing component. It does not supply
unequal molecular diffusion, partial-enthalpy transport, reaction-rate closure,
scalar variance or a species wall function. Matching original physical grid
widths, flow mask and the flow epoch remains the caller's responsibility. The
original turbulence model's applicability and resolution limits remain in force.


### Reservoir-driven species diffusion

`MaterialAdvectionStep::diffuse_species_with_boundaries` adds explicitly
authored fixed-composition reservoirs at exterior or fluid-solid faces.
Each face supplies unit-sum mass fractions and a common finite nonnegative
film resistance in m² s/kg. With the original accepted density and native
half-cell width, `G = A / (w/(2 rho D) + R_film)` in kg/s. Zero film resistance
prescribes composition at the actual face; zero cell D still blocks exchange.
Boundary identities must be unique and ordered by `(cell, face)`; interior
fluid-fluid faces cannot acquire a second transport owner. Unlisted faces
remain impermeable. A solid-facing reservoir does not evolve solid storage.

Exact positive cell/reservoir relaxation maps surround the original interior
half sweeps in canonical forward/reverse order. Signed, gross inward and gross
outward constituent transfers are retained at each face. Original constituent
mass admission subtracts the actual signed external transfers; every boundary
map's nonzero summed carrier transfer is charged to the partition budget.
Both boundary and interior updates count toward work, and prepared coefficients,
transfer accumulators and face readbacks count toward storage.

`diffuse_species_with_boundaries_adaptive` reuses the original complete-interval
coarse/fine refinement. Besides per-cell concentration tolerances, each species
has an explicit integrated-transfer allowance in kg. The estimate sums, over
faces, the largest signed/inward/outward transfer difference divided by three;
opposing reservoirs cannot cancel unresolved transport. The retained fine
inventories and transfer accounts are authoritative, without extrapolation.
These are asymptotic temporal estimates, not certified spatial/model bounds.
This is an isothermal common-D reservoir law, with no inferred selective
membrane, sorption equilibrium, finite reservoir/solid storage or enthalpy flux.

### Native open 3D local SGS momentum and heat

`lattice3::open_boundaries::sgs` prepares Smagorinsky, WALE or Vreman
on the original D3Q19/D3Q27 reservoir geometry. The accepted, force-centered
velocity field uses the native full gradient: quadratic one-sided velocity
reservoir slopes (linear when the second donor is solid), the pressure
reservoir's zero-normal-gradient donor convention, and unchanged interior,
periodic and stationary-wall stencils. Open faces are never interpreted as
impermeable halfway walls. Shared local SGS admission and collision-family,
bulk/higher-mode preservation are retained. The original molecular collision
or a fresh explicit per-cell molecular field is the baseline on each call;
previously augmented collisions must not be reused as molecular transport.

`OpenSgsWorkspace3` retains the actual prepared eddy viscosity rather than
recovering it by subtracting rounded collision viscosities. Its explicit
numerical payload budget covers five f64 values and one native collision per
cell, with actual Vec capacities checked. Serial and pooled native momentum
steps consume these prepared collisions without adding persistent closure
state. Preparation performs no adapter-owned allocation.

`step_serial_with_thermal_sgs` couples that same local eddy viscosity to the
original native open-flow/D3Q7 heat/finite-wall transaction. Its explicit
positive turbulent Prandtl number gives alpha_t = nu_t/Pr_t. Since the actual
D3Q7 scalar sound-speed squared is 1/4, it adds 4*alpha_t to each supplied
molecular thermal tau. The caller supplies a fresh complete molecular thermal
tau field. A separate nonnegative finite alpha_t/alpha_molecular maximum is
an admission limit, never clipping. Nonfinite or underflowed positive modeled
transport refuses; ordinary finite floating-point addition may round a tiny
positive contribution back to the molecular tau. Exact zero SGS preserves
both molecular collision and thermal-tau bits. The coupled workspace adds
one f64 per cell to the momentum workspace.

The existing native owner retains reservoir reversal, scalar-population heat
exchange, symmetric source kicks, advective history, finite-wall storage and
wall conduction. All coefficient admissions precede staging, and its final
heat balance and cancellation gates precede joint publication. A refusal may
change scratch only: accepted momentum/thermal populations, finite walls,
receipts and step counters remain unchanged. Modeled mechanical transfer
remains diagnostic, with no invented viscous-heating source. This is passive
fixed-reference-capacity low-Mach temperature transport with frozen incoming
coefficients, not compressible enthalpy, a dynamic exterior test filter,
nonreflecting-outlet physics or a qualified turbulent-jet result.

Four native regression groups cover D3Q19 TRT and D3Q27 cumulant/full-gradient
reservoir preparation for all three models; exact zero-model joint trajectories
with sources and conducting finite walls; active native composition and an
independently measured fluid/solid heat inventory; and warmed-state budget,
ratio, late thermal, cancellation and retry behavior including positive
sub-ulp modeled transport. These regressions are authored and parsed/formatted;
the native integration suite was not compiled or executed in the constrained
implementation workspace.

## Reciprocal momentum on evolving material liquids

The sealed `MaterialFreeSurfaceMomentum3` interface binds the existing
`MaterialFreeSurface3`, `MaterialFreeSurfaceRheology3` and
`MaterialFreeSurfaceEnthalpy3` owners to the same compact native impulse
transaction. It exposes the original reference-density, length, velocity and
scenario scaling together with a borrowed native wet-cell stage. Implementations
delegate to the original physical speed/unit admission and retained liquid
inventory; composing solvers cannot substitute another material or flow owner.

This instantaneous mechanical source preserves canonical density, actual partial
mass, held force frame and physical nonequilibrium stress. The thermal and
rheological owners retain their accepted source samples, caloric/phase/rheology
history and chronology. No collision, constitutive query, heat deposition,
geometry update or clock advance is inferred. Dropping a refused stage publishes
nothing, and no complete material owner or population field is cloned.

The interface enables reciprocal finite-body force/couple and dissipation
checks in `fs-couple` before native publication, with the same source binding
for frozen, non-Newtonian and thermally evolving liquids.

## Tangential moving walls in native open flow

`OpenBoundaryGrid3` now admits prescribed tangential exterior-wall velocities
for both D3Q19 and D3Q27. The existing native moving reflection and compensated
relative-velocity impulse/work reduction serve closed and open grids. Reservoir
donors and reservoir exchange use the same moving-wall candidate. At open/wall
seams, reflected links stay wall owned; the reservoir cannot conceal wall mass
injection. Every reflected set passes a local zero-mass-source check before
admission. Streamwise moving-wall through-flow is supported, including its seams.
Normal motion and incompatible velocities at intersecting wall faces refuse;
the fixed geometry does not represent swept volume or fresh cells.

Boundary rebinding and serial/pooled steps accept wall-speed updates without
resetting populations. `last_wall_exchange` measures all exterior wall links
after an accepted moving interval, separately from reservoir transport. It
retains wall, population and moving-mass impulses and their angular counterparts,
plus prescribed wall work. Refusal preserves the accepted populations, targets,
time and last receipt. Spatial-force and thermal publication retain this same
transaction. Rheology and SGS gradients sample the actual wall velocities.
The thermal carrier supports stateless normal wall laws on tangentially moving
faces; finite moving-wall storage refuses because solid heat advection is not
represented. Mechanical wall work is not silently converted to heat.

Focused tests cover uniform co-motion with through-flow in all orientations,
separate wall/reservoir momentum closure and seam populations, physical wall
shear, force-epoch publication, and native-kernel parity through analytic
Couette convergence. The linear moving reflection retains the native low-Mach
model; no new high-Mach or moving-geometry accuracy claim is added.

## Spatial velocity traces on native open flow

`OpenBoundaryGrid3` and `ForcedOpenBoundaryGrid3` accept spatial velocity
profiles on all six native face orientations for D3Q19 and D3Q27. They reuse
the specialized owner's `InletFlux3`, profile report and refusal types. Callers
may resample at successive physical times without resetting populations or
force history; two face-sized vector buffers are reused, with no profile
callback or allocation during stepping. Sampling/refusal leaves accepted
populations and time unchanged. A successfully installed control remains held
if a later fluid step refuses.

Optional outward-normal mean control uses compensated deterministic reduction.
Both sampled and corrected vectors respect the native grid/collision Mach
envelope; an optional inclusive speed limit applies to the final corrected
trace. Mixed wall/open seams retain their original per-link ownership and
refuse mean-normal control because their complete velocity is not imposed.
Unrelated wall/reservoir updates preserve profiles. An attached profile must
be cleared before changing its uniform fallback; clearing restores that held
fallback without changing physical state.

The original force-corrected donor stress and regularized reconstruction consume
the held vectors directly. Momentum receipts measure the actual operator, while
thermal advection, rheology and SGS consume the accepted physical velocities.
No separate flow or turbulence source model is introduced. Focused cases cover
time updates across orientations/collision families, momentum and mass closure,
uniform-profile bit parity, native Mach/final-speed refusal, late cancellation,
moving-wall seam ownership, and spatial force epochs. The requested mean is a
volume-velocity control, not an assertion of exact mass flux.

## Native open-grid spatial force epochs

`lattice3::open_boundaries::forcing::ForcedOpenBoundaryGrid3` binds an explicit per-cell force density to the original native D3Q19/D3Q27 open-grid populations. Initial attachment changes only the native half-source coordinates: density, physical velocity, nonequilibrium stress, accepted time and reservoir/pressure controls survive. There is no public inner-grid escape with stale uniform-force macroscopic readback.

A step collides with the accepted opening force and interprets streamed donors and targets with their own endpoint force. Donor physical stress uses the donor source; reconstruction attaches the target source. The fluid body impulse is the trapezoidal sum of the two force epochs, separate from measured reservoir and wall transport. Ordinary steps do not perform an extra source-coordinate reframe. Both source fields must satisfy original collision capability, finite-value, shape and solid-zero admission. Cancellation, Mach or reconstruction failure leaves accepted populations, force history and time unchanged. The wrapper retains original native collision kernels and preallocated grid scratch; holding a spatial force allocates no per-step adapter storage.

Focused checks execute the actual generic-poll native owner cores across five D3Q19 and six D3Q27 collision configurations: physical velocity/stress attachment, periodic trapezoidal impulse, bitwise constant-field equivalence, all three open axes with distinct donor/target sources, measured reservoir closure, and late refusal/exact retry. These extracted-production checks exclude the thin `Cx`/pool wrappers and do not substitute numerical kernels. They do not constitute a full Cargo integration run or a new finite-force accuracy claim for the native collision families.

## Native open-flow Boussinesq temperature feedback

`lattice3::open_boundaries::thermal::buoyancy::BoussinesqOpenGrid3` owns native D3Q19/D3Q27 momentum, D3Q7 temperature, finite thermal walls/conduction and both force epochs as one accepted state. Its explicit affine law is `F_b = -rho_ref * beta_theta * (theta - theta_ref) * gravity_lu`, plus an authored per-cell mechanical background that is held or supplied at the next accepted time. Here `beta_theta = beta_SI * delta_T` and `gravity_lu = gravity_SI * dt^2 / dx`; the law uses fixed reference density while native macroscopic velocity retains the actual population density. Attachment preserves physical momentum/stress and existing thermal history. This is fixed-reference-capacity heat transport with Boussinesq momentum feedback, not variable-density mass transport or a total-energy equation.

Each step transports heat at the accepted force-corrected velocity, previews the endpoint force from transported temperature plus the remaining source half-kick, and performs one original native momentum step. Collision uses the retained opening force; donor/target reconstruction and endpoint readback use their local candidate force. The final scalar source kick uses the admitted native endpoint velocity. Its actual temperature must reproduce the applied force within the caller's absolute component tolerance; the accepted representation retains that applied force. Scalar splitting retains its original temporal-order limits. The separately reported trapezoidal body/buoyancy impulse is not deposited as sensible heat.

The declared relative-density-departure envelope admits cell endpoints/source previews, authored reservoir/fixed-temperature wall data and actual candidate fluid traces at reservoirs and external/internal walls. Reservoir receipts retain the extrema of their reconstructed incoming/outgoing population pairs, including outflow/backflow. Finite-solid bulk temperatures remain governed by the thermal model, not the fluid density law. A refusal or cancellation publishes none of the flow, force, thermal history, wall state, time or receipts. Three explicitly budgeted force fields and the existing native/scalar/solid scratch are reused without per-step adapter allocation; no iterative heat re-solve is needed.

Focused checks execute five actual coupled-owner cases: analytic uniformly heated boxes on both native lattices with reference density distinct from population density; thermal reversal with zero net buoyancy impulse; a conducting-wall open channel with measured two-way feedback and independent fluid-plus-solid thermal closure; final-poll/domain refusal with bitwise exact retry; and fluid-boundary domain enforcement despite admissible mixed cells. Two direct law tests cover signs, temperature-datum shifts, anomalous expansion, zero-expansion background preservation and explicit domain refusal. The numerical-owner tests use extracted production cores and real native/math/sparse kernels, excluding only execution-context wrappers. Full Cargo/public-context integration and the planned developing-channel Nusselt/refinement reference remain unverified by these tests.

### Mixed buoyant and prescribed-force convection

`BoussinesqOpenGrid3::advance_serial` and `advance_serial_with_tensor_diffusion`
accept a caller-owned `background_at_end: &mut Vec<[f64; 3]>`. The field is the
mechanical force density at the next accepted lattice time, before adding
thermal buoyancy. It covers every original cell, is finite and vanishes on
solid cells. Spatially nonuniform and time-varying backgrounds use the same
native D3Q19/D3Q27 momentum and D3Q7 heat owners as held-background steps.

Collision uses the retained total force `F_n = B_n + F_b(theta_n)`.
Reservoir reconstruction, endpoint physical velocity and the final scalar
source half-kick use `F_(n+1) = B_(n+1) + F_b(theta_transport + Q/2)`.
The original final force-consistency gate checks the actual accepted
temperature against this endpoint background. Thermal advection retains its
opening physical velocity and existing source split; no additional temporal
order is claimed.

`BoussinesqOpenStep3::background_impulse` reports the fluid-cell sum of
`(B_n + B_(n+1))/2`. `body_impulse` uses the corresponding total-force
endpoints, and `buoyancy_impulse` subtracts the mechanical background using
compensated accumulation. Wall and reservoir momentum remain separate.
Mechanical and buoyant work are not deposited as heat; the original reservoir,
finite-wall, immersed-solid and thermal-source accounts remain authoritative.

After the final cancellation gate, accepted flow, total force, thermal state,
wall histories and mechanical background publish together. Success exchanges
the caller's endpoint buffer with the old accepted background for reuse;
refusal or cancellation preserves both the input buffer and accepted epoch.
The owner retains its existing three force fields without a per-step force
allocation. Caller endpoint storage is outside `maximum_force_storage_bytes`,
which bounds retained field payload. An unchanged endpoint background retains
held-background arithmetic. Focused regression definitions cover affine
forcing/heating on both lattices, identical-endpoint equivalence, spatial
forcing with tensor/finite-wall heat closure, and refusal with exact retry.

## Constituent mass on accepted material transfers

`MaterialAdvectionStep::advect_species` transports 1..=128 fixed-order
constituent inventories on the retained native three-dimensional material
advection routes. Per-cell constituent **mass** is authoritative; its initial
partition is checked against the accepted step's original mass and its endpoint
partition against that same step's represented endpoint mass. Incoming boundary
planes and local sources have separately authored mass fractions selected by
each retained route's actual origin. Extraction and backflow use their original
donor/origin rules. No second density solve, inferred partial volume, or
composition-dependent EOS is introduced.

Each species reuses the original row-scaled mass-times-specific-quantity donor
matrix and the existing BiCGStab recurrence. One native ILU(0) factor is reused
across species, with only one scalar solve workspace live at a time. The caller
bounds component work and payload storage. Per connected component, fraction
amplitude scaling keeps an entirely trace species representable in the native
Krylov norm; absent species retain the original exact zero-RHS path. This does
not promise arbitrary relative accuracy for a tiny local tail of a bulk species.

Paired accepted transfers determine endpoint constituent masses. Independent
readback substitutes those actual masses into the unscaled extensive donor
equations, checks each species' global boundary/source balance, and audits the
initial, endpoint and route partitions against the accepted carrier. Each
species has its own absolute mass-defect budget, including trace species;
partition errors have a separate declared kilogram allowance. Pure-material
derived fractions can differ from one by admitted carrier/constituent roundoff.
There is no clipping, seeded trace or composition renormalization. Negative or
nonrepresentable inventories/transfers and exceeded budgets return no candidate.
Cancellation, including after completed components, leaves all caller data and
the accepted carrier step unchanged.

This is the native first-order implicit donor advection scope. Species diffusion,
reactions, mixing heat, thermodynamic composition feedback and a higher-order
advection scheme remain with their respective future/material owners.

## Thermal evolution coupled to resolved Carreau film momentum

With `implicit-enthalpy-coupling`,
`film::transient::thermal::ThermalCarreauFilm` binds the retained original
wall-normal velocity to the original identity-bearing caloric timeline.
Ordered fluid cells must form the same straight, uniform wall-normal mesh:
each cell retains volume `A*dy`, mass `rho*A*dy`, the same original caloric
material and the same accepted time. Existing solid or wall cells may share
the native thermal timeline. No capacity, conductivity or liquid phase
occupancy is inferred from a viscosity law.

Each native backward-Euler thermal residual, derivative probe and independent
final evaluation restarts the SAME accepted momentum profile. Interior shear
face temperatures interpolate the original native caloric cell temperatures;
the callback supplies the actual wall and free-surface fluid-side traces
from its original interface/conduction model. A solid bulk temperature is not
automatically an interface temperature. Those endpoint source temperatures
are held during the original midpoint momentum solve. This composition is
first order in coupled time, even though momentum at a prescribed fixed
temperature retains its midpoint scheme.

The original local viscous-loss partition supplies each fluid cell's heat
exactly once. Signed body, moving-wall and free-surface work remain mechanical
exchange; neither kinetic-energy change nor algebraic defects are deposited
as heat. Original native thermal acceptance precedes a separate combined
`delta(thermal + kinetic) - other heat - mechanical work` gate. Other thermal
power is measured directly by the native owner before viscous heat is added;
the combined account uses the independently measured stored-energy change.
Reports retain actual deposited heat, its partition residual, kinetic change,
mechanical work and thermal/nonlinear residuals separately.

The owner reuses its momentum trial arrays, temperature buffer and local heat
buffer. Every thermal trial resets only the original accepted velocity/time;
rejected Newton candidates do not become physical histories. Endpoint and
aggregate momentum constitutive budgets bound nonlinear work. Native thermal
and Krylov allocations remain. Cancellation is latched for the whole advance.
No accepted velocity, thermal state or time changes unless all final gates
pass; an explicit `into_parts` handoff permits native events and later rebinding.

This is a locally parallel, fixed-thickness, fixed-reference-density liquid
film. The thermal source must remain in its admitted liquid domain.
Streamwise transport, evolving interfaces, thermal expansion and phase
occupancy are not introduced by this coupling. Nonlinear convergence alone
does not establish temporal or spatial accuracy.

The coupling's default Newton controls use an inner forcing floor of `1e-4`,
appropriate for its finite-difference Jacobian and caloric/constitutive roundoff.
The original outer residual and physical energy gates remain independent;
explicit caller-supplied solver controls are unchanged.

All 37 focused source-linked film tests passed, including four coupled
regressions for viscous self-heating and original-source feedback, an independent
two-cell conduction/bath balance, first-order refinement against independent
coupled momentum/thermal integration, and joint refusal/cancellation/retry.
Executed code includes the original material laws, math, tridiagonal, thermal
backward-Euler, finite-difference and Newton/Krylov kernels. The harness omits
unrelated crate modules and uses an uninhabitable context type for unexercised
runtime entrypoints. These results do not establish full-crate Cargo or
asupersync runtime integration.

## Resolved transient Carreau films

`set_face_temperatures` admits an explicitly supplied wall-normal temperature
field at all N+1 original shear faces. Every face retains the SAME thermal
Carreau source and physical rate ceiling, with its own prepared temperature
shift. Interior stresses/tangents and the free-surface stress inverse use those
actual local preparations; no mean temperature or fitted apparent viscosity is
substituted. Uniform fields retain the original compact path. Spatial updates
reuse two prepared buffers after first use and publish only after complete
source/scaling admission and the final cancellation check. Updating temperature
never replaces velocity history, changes density/thickness, or advances time.

The Newtonian, original-source-bound Newtonian and Carreau `step_with_heat`
entries share one shear-segment partition. Each exports a nonnegative
viscous-energy increment per original
wall-normal cell, in J/m² of wetted wall area. Interior shear-segment work is
split equally between its two adjacent cells; the wall and surface half-segment
losses enter only their adjacent cells. The sum is the original momentum work
account's viscous heat, within summation roundoff. The surface contribution uses
the imposed traction and retains the separate constitutive-inverse defect.
The caller supplies its actual represented wall area and original caloric
owner; neither a heat capacity nor a temperature rise is inferred. These
increments already exclude kinetic storage, so body, wall and surface work
must not be deposited again. Heat output, velocity and time publish together;
shape, source, solve, speed and cancellation refusals retain every accepted
value. Reusable O(N) scratch avoids momentum-step workspace allocation.

G1/G3/G4 spatial-field tests compare against independent original-law stress
bisection and continuum integration under a linear temperature gradient and
opposing surface traction; they check spatial refinement, uniform-path parity,
local heat partition and global work, original face/source retention, thermal
rebind refusal, and cancellation/retry including the caller's heat buffer.
Native Newtonian checks also recover exact signed Couette cell heating and
bit-identical momentum/report behavior with and without heat export.
The temperature field is prescribed over each momentum interval. These APIs
provide spatial source feedback and conservative heat transfers, not a second
thermal inventory, conductivity closure or automatic thermal time integrator.

`film::transient::TransientCarreauFilm` retains the full wall-normal velocity
history under signed body acceleration, translating-wall velocity and imposed
surface traction. Each implicit-midpoint face evaluates the original thermal
Carreau stress and positive stress tangent. The shared `fs-material`
`PreparedCarreau` owner retains the original temperature shift and constitutive
inverse; no fitted apparent viscosity or summed steady profiles are substituted.
Bounded Newton/backtracking uses the native tridiagonal solve and reusable O(n)
arrays. Accepted trials retain their constitutive/Jacobian evaluation.

Newtonian and nonlinear films share one finite-volume momentum and work account,
including both boundary half-cells. The nonlinear surface inverse exposes its
stress defect separately. Initial, midpoint and final rates remain within the
authored source-supported ceiling; velocity, source, iteration, work and final
cancellation refusals leave accepted history and time unchanged. A prescribed
temperature update retains momentum, density, geometry and the physical rate
ceiling while re-preparing the original source and diffusion scale.

G1/G3/G4 tests cover the original Newtonian limit, independent constitutive
integration and spatial refinement, nonlinear time refinement, reversing and
moving-wall work, exact rest/Couette limits, temperature continuation and atomic
retry. This is fixed-thickness, locally parallel generalized-Newtonian flow.
It adds no elastic memory, thermal evolution, free-surface geometry, entrance
flow, capillary waves or turbulent-transition prediction. Midpoint time and
wall-normal resolution remain accuracy requirements.

## Resolved transient laminar films

`film::transient::TransientFilm` evolves the wall-normal velocity of the same
constant-property Newtonian film used by `NusseltFilm` and `ShearedFilm`.
The local equation is rho u_t = rho g + mu u_yy at fixed thickness, with
prescribed tangential no-slip wall speed and surface traction. Uniform finite
volumes and implicit midpoint retain startup, decay, reversal and counterflow;
drive changes do not replace the accepted profile by a steady parabola.

The original native tridiagonal solver uses O(n) work/storage, reusable scratch
and a cached factorization for unchanged durations. Every numerical traversal
polls cancellation at least every 128 cells. The caller supplies physical-speed
and original-equation residual limits. Invalid source/geometry/drive, failed
solve, exceeded limits or cancellation publishes neither velocity nor time.

Receipts separate body, wall and surface work, represented kinetic change,
nonnegative viscous heat, and momentum/work defects. The prescribed surface
stress reconstructs the actual boundary trace across its half-cell; that
half-cell's traction work and viscous loss occur exactly once. Heat remains a
reported transfer, not a second thermal inventory or inferred temperature rise.
`MaterialFilmLiquid::transient_film` retains the original liquid source,
viscosity basis, query and receipts beside this native history. No additional
material selection or apparent-viscosity fit is introduced.

G1 tests compare startup with a continuous mixed-boundary sine expansion,
measure second-order spatial and temporal refinement, and recover the driven
Nusselt/sheared steady limit. G0/G3/G4 checks cover signed work and momentum,
common-frame boosts, reversal, rest, domain refusal and exact cancellation
retry. These tests exercise the numerical owner; the material adapter also
has source-to-native parity regressions for both viscosity bases. Verification
of the numerical owner uses exact production film/math/tridiagonal sources;
full Cargo and source-adapter integration remain separate verification gates.

This is fixed-thickness, locally parallel, laterally wide Newtonian laminar
flow. It does not evolve a free surface, infer turbulent transition, resolve
entrance effects or change material temperature/composition. Midpoint may
oscillate on unresolved stiff modes; conservation and linear stability do not
replace temporal or spatial resolution checks for the forcing bandwidth.

## Source-resolved Carreau viscosity in nonlinear slender jets

`PeriodicSlenderJet::new_carreau` retains an original
`ThermalCarreauYasuda` law and fixed physical temperature, density and surface
tension. Each actual explicit or implicit RK stage queries the material at
`gamma=sqrt(3)*|du/dz|`, the invariant of incompressible uniaxial extension.
The conservative axial stress is `3*mu*du/dz`; the same sample supplies
`3*mu*(du/dz)^2` to the existing integrated viscous-loss account. No lattice
tau cap, fitted Newtonian substitute or second constitutive implementation enters.

The explicit viscous step bound retains both apparent and differential
viscosity. The apparent conductance is necessary near strongly thinning states:
a small tangent alone can allow an Euler stage to overshoot through zero strain.
Every explicit trial stage rechecks the bound and the original source support.
Source refusal or cancellation preserves accepted area, momentum, geometry and
time. The explicit path retains its original numerical workspace.

Regression definitions cover exact Newtonian-degeneracy trajectories, shifted
temperature, independent tensor stress and heat, conservative momentum/volume,
strong-thinning dissipation, the viscous timestep bound, and refusal/rollback.
These definitions require execution in the repository's Rust validation lane;
they are not a claim that a new convergence or performance campaign has run.

`ImplicitViscosity` retains the original IMEX-SSP3(4,3,3) tableau for Carreau
jets. Each fixed-area nonlinear stage solves the original conservative stress
with its exact positive differential tangent. The native cyclic LDL/rank-one
kernel resolves the symmetric positive definite Newton Jacobian in O(N) work.
Five additional arrays are prepared once for a non-Newtonian selection; successful
stepping allocates nothing. Exact Newtonian degeneracies retain the native direct
solve and arithmetic. `SlenderJetImplicitControls` bounds corrections, backtracks
and relative stage residual; the outer step budget bounds the number of stages.

The source-weighted mean is separated before nonlinear residual scaling, so a
uniform translation does not become the convergence tolerance for weak shear.
Backtracking re-queries the same original material; unsupported algebraic trial
rates may be rejected without clipping or extending source support. A search
with no supported candidate preserves its original source diagnostic. Actual RK
stage refusal, exhausted nonlinear work and cancellation roll back the full
advance. Returned momentum and source rates are rechecked after mean restoration.

Positive RK weights integrate physical viscous loss at the solved stage rates.
The fixed-stage kinetic balance distinguishes molecular loss from implicit
numerical damping; full IMEX endpoints retain a signed truncation remainder.
Third-order temporal behavior requires resolved dynamics and nonlinear residual
error below time-discretization error. The embedded second-order estimate remains
a local indicator, not an unconditional accuracy or monotonic-energy guarantee.
Focused regressions use an independent nonlinear alternating-face decay law,
original-stress/tangent differences, stage energy balance, large translations,
stiff shear-thinning evolution, source-domain backtracking and transactional retry.

Temperature does not evolve, loss is not deposited as heat, and there is no
elastic memory or independently calibrated extensional response. Original
slenderness/resolution limits and pre-pinch stopping remain.

## Adaptive source-owned turbulent gas heat and species transport

`FvFlow::advance_turbulent_gas_species` advances a complete physical window
through the original material network's coarse-versus-two-half-step controller.
Native geometry, accepted carrier masses and momentum `nu_t/Sc_t` face
conductances remain fixed for that flow epoch. Every probe refreshes actual
composition, donor partial enthalpy and the original caloric/EOS/transport
endpoint. Only accepted fine steps contribute to the compensated stream ledger.
The result `TurbulentGasSpeciesAdvance` retains local species/temperature error
ratios, bounded work and accepted-step statistics, original-to-final physical
audits and carrier drift against the same original accepted masses.

`diffuse_turbulent_gas_transport` and `advance_turbulent_gas_transport` add
turbulent heat to the snapshot and adaptive entries respectively. Source-derived
`k_t = rho cp(T,Y) nu_t,heat / Pr_t` uses the retained thermal eddy field and
positive Prandtl number; momentum `nu_t/Sc_t` remains the independent species
law. Graded half-cell series resistances give antisymmetric internal heat
`Q = dt K (T_left - T_right)`. Each probe samples heat and species from the
same original state and sends both to one material reconstruction. Donor partial
enthalpy and formation references are carried once; additional heat has its own
accepted-only signed ledger. Solid, exterior and zero-conductivity faces supply
no heat in this component. Molecular conduction and external wall heat retain
their separate owners.

The original 256-parcel and 4096-stream bounds apply; heat and species share the
stream/whole-window report bound. Typed withdrawal and caloric-bracket overshoot
may shorten adaptive trials. Other source-domain errors are fatal. Local-error,
step and evaluation limits bound refinement; rejection/cancellation never
publishes a partial candidate. Fixed whole-window mass, each-species, actual
source-relative energy and endpoint entropy budgets are never multiplied by
the number of accepted steps. The source-relative energy audit removes common
formation references before arithmetic. Carrier drift is also measured against
the initial flow state, not reset after every substep.

Snapshots and the retained adaptive fine path remain first order in time.
Local control is not a global-error or turbulence-model validation certificate.
The held eddy/geometry epoch does not refresh pressure or momentum; subsequent
flow consumers must use the returned material endpoints. The original planar
or finite-volume turbulence closure limitations remain.

Focused regressions are authored for analytic species relaxation, refreshed
partial enthalpy, typed withdrawal recovery, variable-heat-capacity heat
exchange, formation-reference invariance, an independent combined RK4
reference, accepted-only heat accounting, independent Prandtl/Schmidt response,
thermal-bracket recovery and atomic late work/source/cancellation failure.
Scoped pinned formatting and whitespace checks passed. Full native integration
execution remains blocked by the pinned asupersync dependency exceeding this
environment's build resources before the changed crates compile.


### Partial enthalpy on the same constituent diffusion transfers

`MaterialAdvectionStep::diffuse_species_with_enthalpy` adds thermal storage to
the original common-diffusivity owner. One callback supplies the original
partial specific enthalpies in J/kg at each conducting face, frozen for the
interval. Every actual signed constituent transfer contributes `sum(h_k dm_k)`;
one cell loses exactly the energy its neighbor receives. The existing positive
mass relaxations, graded-grid conductances, ordered sweeps and mass/partition
budgets remain the only species solve. Source face queries occur once before
the sweeps, and both inventories publish together only after final admission.

Because the barycentric diffusive carrier flux is zero, a shared enthalpy datum
is subtracted before products. Species-specific formation/latent differences
remain. Error-free subtraction and product transforms retain their low parts
until opposing species contributions cancel; compensating already-rounded
products alone is insufficient. Compensated cell increments update `m*h` storage
in centered coordinates; an absolute joule budget includes local unrepresented
work, closed-domain storage/transfer defects, and rounded carrier-flux defects
weighted by the partial-enthalpy range. The endpoint audit also retains low
parts of mass times the represented enthalpy change; opposing local storage
errors cannot disappear through rounded products or a zero global balance. Large formation-energy totals cannot
enlarge that budget. Thermal face/scratch/output storage is included in the
original allocation ceiling. Source error, cancellation or failed mass/heat
admission returns no partial candidate.

The returned endpoint enthalpy already contains the reported diffusion energy.
Consumers must invert it with the returned composition using their original
caloric law, re-evaluate face properties as needed, and accept both inventories
together. This component supplies no EOS, source identity, temperature solve,
reaction or pressure-work source, arbitrary multicomponent molecular law, Soret
term or coupled temporal-order guarantee. Freezing partial enthalpies is a
declared approximation for a changing state. Focused G0/G1/G3 regression
definitions cover analytic binary relaxation/energy flux, multicomponent
isothermal caloric inventory, common-datum invariance, retained small heat amid
opposing formation terms, unchanged mass-only results, and unrepresented-work,
source, storage and cancellation refusals. These
new regressions require native execution; the authoring environment went
offline before that check completed.

## Original gas calorics on constituent diffusion

`GasSpeciesEnthalpySource` binds the full canonical species bank, original
NASA-9 formation references, pressure context and transport data to the native
species/enthalpy diffusion owner. Every cell preparation converts the actual
constituent masses to molar amounts, checks their partition of the immutable
carrier to machine roundoff, and inverts the original composition-dependent
enthalpy law. The same resolved state supplies temperature, Cp, conductivity,
viscosity and EOS density. Zero local abundance retains the full source bank.
The full bank must admit a continuous increasing caloric branch, including
species absent from a particular cell; discontinuous seams do not get hidden
by a diffusion step.

Face partial specific enthalpies come from each original species law at the
physical face temperature, interpolated with the actual graded cell widths.
`GasSpeciesDiffusivity::Binary` reuses an explicit `BinaryDiffusionMaterial`
for exactly its two species, refreshing its m²/s law at the reconstructed T/p
while preserving auxiliary source coordinates and property selection. The
`Frozen` variant retains an explicitly authored common mixing field, including
zero-diffusion barriers. Viscosity never substitutes for molecular diffusivity,
and one binary coefficient cannot silently represent a larger molecular mixture.

`diffuse_thermochemical_species` advances the existing positive species kernel
and its exact constituent-carried energy account, then reconstructs the endpoint
from both returned inventories. Initial/endpoint caloric residuals and the
original transport defect share an absolute joule budget. Owned prepared arrays
are reserved from the native storage ceiling before allocation; immutable source
banks and temporary original material-evaluation scratch remain additional.
Source, identity, domain, storage or cancellation refusal returns no partial
candidate. The returned h already contains diffusion energy.

This entry freezes D and face partial h over its declared interval; repeated
steps must refresh them. It is an isobaric fixed-carrier constitutive closure:
EOS density is returned for subsequent flow/pressure coupling and never replaces
the accepted mass. It introduces no pressure work, reactions, arbitrary
multicomponent molecular friction, Soret term or coupled temporal-order claim.
Regression definitions cover variable-Cp analytic binary exchange, formation
reference invariance, graded physical face interpolation, absent species,
source-dependent diffusivity, validity axes and atomic refusals.

## Enthalpy through original species reservoirs

`diffuse_species_with_enthalpy_and_boundaries` observes each actual signed
constituent transfer from the existing positive reservoir/interior maps.
One original partial-specific-enthalpy vector in J/kg is supplied per active
boundary face. The same local face law applies to counterdiffusing species;
the common carrier flux is zero. The original film/half-cell conductance,
ordered sweeps, positive species inventory and mass/partition gates remain
authoritative. Empty boundaries reproduce the closed enthalpy entry exactly.

Interior and boundary energy share the original centered partial differences
and error-free products. A boundary credits the same retained energy components
to its cell and the independent external account once. The shared endpoint
audit compares represented storage and transported energy with external heat,
rather than incorrectly requiring a closed-domain zero. Per-cell storage
mismatches and every boundary mass-frame defect still cost absolute joules.
Internal cumulative energy and species receipts retain low parts; public
scalar readbacks round only at the reporting boundary.

Per-face and aggregate receipts retain signed energy and positive/negative
NET map energy after cancellation among species. Large opposing reservoirs
or a reversal during an interval cannot disappear behind net domain heat.
This energy is already contained in the returned enthalpy. The generic
reservoir observer forwards both kinds of transfers; the mass-only consumer
is a zero-sized no-op. Active thermal face caches, retained accounts and
readbacks are admitted against storage before allocation, and source refusal
or cancellation publishes no candidate.

`GasSpeciesEnthalpyBoundarySource` binds the exact canonical reservoir list to
explicit PHYSICAL face temperatures and the original gas source. Its partial
enthalpies use the same NASA evaluations as interior faces, including species
absent from a cell. `diffuse_thermochemical_species_with_boundaries` resolves
the endpoint composition and total enthalpy together, refreshes gas transport,
and charges initial/endpoint caloric residuals alongside native open-energy
defects. Species-specific formation shifts transform external chemical energy
and stored enthalpy together; a common mass-specific datum cancels from the
zero-carrier diffusion flux.

Reservoir composition and a finite mass-film resistance do not prescribe a
thermal interface temperature. The supplied face trace must come from its
own thermal boundary/interface owner; no donor-temperature, Fourier-film,
reaction, extra mixing-heat or pressure-work source is inferred here. The
reservoir is fixed and owns no evolved storage. Source coefficients are frozen
over this entry's interval; coupled temporal refresh has its own driver.

Focused regressions include exact finite-film exchange, simultaneous interior
and boundary transport, empty-boundary equivalence, signed/gross heat reversals,
common-reference invariance, retained formation and species-product remainders,
and source/storage/cancellation refusals. Original gas integration fixtures
exercise an explicit face temperature distinct from cell temperature and
formation-reference covariance.

## Adaptive species and enthalpy with original material feedback

`SpeciesEnthalpyDiffusionModel` provides one original state-dependent caloric
and transport closure. `diffuse_species_with_enthalpy_adaptive` prepares the
original state, advances a positive half-step predictor, and prepares its
midpoint properties. A full-step corrector then starts from the original
substep inventories using those properties. Predictor heat is never deposited
again. The endpoint source image supplies the next substep and is returned
as `state` on acceptance; a gas caller obtains refreshed Cp, conductivity,
viscosity, EOS density and diffusivity without repeating the final inversion.
`GasSpeciesEnthalpySource` implements this same generic closure directly.

Complete intervals at n and 2n physical refresh steps estimate smooth
second-order time error through |fine-coarse|/3. Temperature in K, specific
enthalpy in absolute J/kg, and each constituent mass fraction have independent
allowances. Formation-energy totals cannot enlarge a thermal tolerance.
The accepted candidate is the original positive fine trajectory, with no
extrapolation, clipping or composition normalization. The prior fine solution
becomes the next coarse comparison. This is estimated temporal accuracy,
not a certified global error bound; sharp constitutive changes and phase corners
need their own source validity and resolution controls.


Complete relaxation can make coarse and fine endpoints agree after both have
skipped a thermal transient. The driver therefore compares sampled D and
error-free centered species partial-enthalpy differences at the start and
midpoint. For changing coefficients, the accepted outer face relaxation must
not exceed one half. That rate is recovered over the full physical refresh
step, so increasing the inner frozen sweeps cannot bypass the check. Exact
unchanged coefficients retain the large-step path. This anti-saturation check
uses counted face queries without another source inversion or face buffer;
it supplements the temporal estimate rather than certifying a constitutive law.

Every predictor, corrector and discarded trajectory contributes to total
pair/component work and source-query counts. The storage gate reserves live
comparison states and three source preparations before allocation; obsolete
coarse source images are dropped before finer work. Original source-table
storage and temporary source-evaluation scratch remain additional. A refused
source state, depleted budget or cancellation publishes no candidate.

Accepted physical correctors accumulate their original constituent-mass and
enthalpy defects. A separate initial-to-final audit uses error-free subtraction
and products for represented storage. Internally, the frozen kernel returns its
existing compensated energy increments by move; both components survive
substep composition, even when a public scalar heat readback rounds away a
small term beside formation energy. Initial/final caloric-inversion residuals
remain explicit in the joule budget, and the largest accepted-stage residual
is reported separately. Rejected trajectories consume resources but are not
added to the accepted physical history.

Focused original-kernel regressions cover the isothermal frozen limit,
independent variable-Cp analytic convergence, a smooth nonlinear heat-of-mixing
case that defeats unguarded step doubling, simultaneous thermal/diffusivity
feedback against an independently refined ODE, source/work/storage/cancellation
refusals, and exact conservation of small heat between large opposing formation
references. A separate integration regression couples the original NASA source,
binary material diffusivity and the native carrier to an independent two-cell
mass/energy ODE. These remain fixed-carrier constitutive diffusion steps;
pressure/volume evolution, reactions and nonideal molecular cross-effects
require their original physical owners.

## Local SGS with general prescribed moving walls

`Grid3::prepare_sgs_with_moving_walls` and the corresponding serial/pooled
step entries use local Smagorinsky, WALE or Vreman closures on selected fixed
halfway-wall links, including interior solids and spatially varying tangential
drives. The complete wall-link field remains with the native moving bounce-back
owner. SGS and kinematics share each missing Cartesian neighbor's actual axial
wall-link velocity; diagonal-only contacts do not invent an axial wall.

The ordered link mapping uses fixed storage and no per-cell probe search or
additional field allocation. Existing collision-family rates, supplied molecular
viscosities, Mach/admission checks and eddy-to-molecular limits remain intact.
The collision field is prepared from fresh molecular values on each call.

`ForcedGrid3` exposes the same preparation, held-force steps and endpoint-force
advances. The accepted physical velocity uses the current `F_n/2` correction
before gradient fitting; `F_(n+1)` admission and publication remain in the original
force owner. Native population, force, clock and wall-load transactions are reused.
SGS diagnostics report instantaneous modeled shear transfer and add no heat.

Regression groups are provided for analytic wall gradients on both lattices and
all three models, zero-drive/disabled controls, serial/pooled/native equivalence,
force histories, cancellation and refusal. This remains the exploratory local
closure on fixed staircase geometry; it does not add wall functions, geometry
remapping, a turbulent-flow qualification or an integrated energy balance.


### Molecular stress with phase resistance and changing density

`advance_viscous_phase_resistance_thermal_volume_momentum` composes the
existing full-stress molecular owner with source-resolved mushy resistance,
physical heating, thermal volume flow and conservative M/Mh/P transport.
Original material receipts supply molecular viscosity; actual M/V supplies
inertia. Darcy resistance resolves the accepted endpoint phase on the same
material chart. Molecular and Darcy heat each enter caloric storage once;
numerical losses, reservoir work and source-pressure work retain their separate
accounts. Noncancelling caloric and mechanical defects share the existing
energy budgets. Source, solve, cancellation or budget refusal publishes no
partial inventory.

This is a first-order split with frozen original molecular viscosity and
endpoint Darcy coupling. It provides neither simultaneous nonlinear shear
coupling nor exact solid locking, moving skeletons or free surfaces. Initial,
incoming, heated and transported phases must remain mobile.

Four new focused regressions cover closed kinetic/caloric balance with both
heat sources, unequal phase densities and volume response on every axis, exact
fully-liquid parity with the original viscous route, and atomic source/solid/
cancellation refusals with retry. These and three existing phase-resistance
tests passed in a native production-source harness; this is not a full workspace
build or a temporal-convergence certificate.

## Original material rheology with curved-wall SGS

`MaterialScaling::step_curved_walls_with_material_rheology_and_sgs` couples
the original source response to the native off-lattice wall owner. Each fluid
cell is sampled once at incoming physical deviatoric shear. One accepted
Guo-centered full gradient supplies both molecular rheology and the explicit
Smagorinsky, WALE or Vreman closure. Selected axial wall links use their actual
intersection distances and supplied wall velocities; diagonal intersections
retain the original BFL population/load geometry.

The exact same frozen wall velocities drive gradient, population exchange,
impulse, torque and wall work. Original local molecular viscosity supplies the
eddy/molecular ratio gate and molecular diagnostic. Collision family, TRT
product, independent bulk and higher-order controls remain unchanged. There
is no second gradient, full-grid clone, implicit material solve or heat deposit.
The original resolved-shear nonlinear-material approximation and fixed
reference-density/scaling limits remain explicit.

This reuses `BoundaryRheologySgsWorkspace3` and its two source buffers.
Applied material samples and the SGS/selected-wall report publish only after
the original curved collision/stream/load transaction succeeds. The complete
returned native receipt also retains other walls, body force, BFL mass change,
open and periodic transport. Refusal/cancellation leaves accepted populations
and material samples unchanged.

Three focused regressions are authored for affine/quadratic wall profiles on
every axis and all five collision families, independent material/eddy
viscosities and momentum/angular-momentum/work accounting, exact constant-source
native recovery, single source reads, and late-refusal replay. Scoped pinned
formatting and whitespace checks passed; native execution remains pending the
pinned dependency build's resource requirement. This composition adds no
turbulent-flow, unresolved roughness or wall-resolution validation claim.

## Source molecular conductivity in turbulent gas transport

`FvFlow::diffuse_turbulent_gas_transport_with_molecular_heat` and
`advance_turbulent_gas_transport_with_molecular_heat` explicitly compose
original gas molecular conduction with the existing turbulent heat/species
update. Each current material state supplies
`k_total = k_source + rho cp nu_t,heat / Pr_t`. These total cell conductivities
enter the graded half-cell series resistance; independently harmonic molecular
and eddy face conductances are not added. Molecular heat therefore persists
when either cell's thermal eddy viscosity vanishes.

Every adaptive probe refreshes original conductivity, density, heat capacity
and donor partial enthalpy at its own state. One original caloric update receives
all heat and species streams. Its accepted-only heat ledger reports total
molecular-plus-turbulent paired heat; callers must not add another molecular
conduction step for the same interval. Existing turbulent-only entries retain
their original component semantics. This supplies neither molecular species
diffusion nor Soret/Dufour, external-wall or reservoir transport.

The original bounded first-order snapshot/adaptive owner, formation-reference
invariance, fixed mass/energy/entropy allowances and cancellation behavior
remain. Focused regressions cover the zero-SGS molecular limit, total half-cell
resistance, source retention and adaptive accuracy. New numerical code/tests
are scoped-formatted and statically reviewed; native execution remains blocked
before the changed crate by the pinned dependency build's resource requirement.

## Transient native thermal boundary inputs

The existing `ThermalOpenField3::set_reservoir_temperature` supports changing inlet and backflow reservoir temperatures without reconstructing the scalar field. `set_wall_boundary(face, law)` now updates an existing external stateless wall uniformly over its original face; prescribed temperature, signed inward heat flux, ambient-film transfer/data and insulation can change between accepted intervals. `set_internal_wall_boundary(law)` changes only the fallback law on unassigned fluid-to-solid links. Explicit immersed bodies retain their link ownership, stored heat and conduction graph. These setters preserve native face geometry, scalar populations, flux history and accepted time; they cannot create or replace finite thermal storage.

`set_finite_wall_heat(face, heat_per_step)` changes the external source of each retained finite wall cell in original x-fast surface order. The field covers the complete wall array, is finite, and is zero in slots adjacent to solid grid cells. Negative source values represent cooling. Stored temperatures, capacities, contact resistances and prepared wall-conduction factors remain unchanged; no energy is deposited by the setter. The original finite-wall solve credits each active source once and balances it against fluid heat, solid storage and reciprocal wall conduction. `wall_states` exposes current sources, while `boundaries` retains the original finite-wall declaration.

All inputs are admitted before mutation. A refused update leaves both inputs and physical state unchanged. A subsequent failed or cancelled step retains the newly authored inputs and the prior physical epoch, matching reservoir-update semantics; the next attempt rebuilds scratch from the retained state. The Boussinesq owner forwards all three APIs, pre-admits prescribed wall temperature through its density-departure envelope, and retains the original actual-trace checks for flux and ambient-film laws. Native held-force, varying-force, scalar-diffusion and tensor-diffusion steps use these same boundary owners without a second heat solve or extra source allocation.

## Spatial-force and buoyant open-boundary SGS

`ForcedOpenBoundaryGrid3` now prepares and advances the existing local
Smagorinsky, WALE and Vreman closures with each cell's accepted spatial Guo
half-force. Held-force and next-force entries share the native reservoir/wall
gradient and original collision/stream transaction. SGS samples only the
opening physical velocity; the supplied endpoint force retains its separate
native reconstruction and trapezoidal impulse semantics.

`BoussinesqOpenGrid3::step_serial_with_sgs` and `advance_serial_with_sgs`
compose that same momentum closure with `alpha_t = nu_t / Pr_t` and the
original D3Q7 molecular thermal relaxation. Fresh molecular collision/thermal
fields are explicit, augmentation is frozen for the interval, zero SGS retains
their exact bits, and eddy/molecular limits refuse rather than clip. One shared
coefficient preparation serves passive and buoyant thermal owners without
additional per-step arrays or population clones.

The original joint owner retains heat-source half-kicks, reservoir/backflow
exchange, finite conducting walls, temperature-dependent endpoint buoyancy,
and separate body/background/buoyancy impulse accounts. Successful publication
alone advances flow, heat, walls, force epochs and clocks or exchanges a supplied
endpoint-background buffer. Late source, ratio, force or cancellation refusal
does not publish a partial physical state. Modeled shear transfer remains
diagnostic; neither it nor mechanical-force work is added to heat.

Four focused regression groups cover analytic spatial-force gradients and
native composition, exact zero-SGS trajectories, active buoyant feedback with
independent fluid/solid heat and impulse accounts, and warmed late-refusal
replay. The first three exercise D3Q19 TRT and D3Q27 cumulant lattices with all
three closures. Scoped pinned formatting and whitespace checks passed before the execution
environment went offline; the final source was recovered from those reviewed
edits through GitHub. Native execution remains pending the pinned dependency
build's resource requirement.
This is the existing explicit fixed-reference-capacity Boussinesq split with
local exploratory closures, not a variable-density total-energy model,
buoyancy-production turbulence closure or new temporal-order qualification.

## Forced open-flow SGS momentum and heat

`ForcedOpenBoundaryGrid3::step_serial_with_thermal_sgs` and
`advance_serial_with_thermal_sgs` combine explicit local Smagorinsky,
WALE or Vreman momentum/heat transport with held or advancing spatial force
fields on D3Q19 and D3Q27. Preparation uses the accepted physical velocity,
including each cell's Guo half-force correction F_n/2; a proposed endpoint
force does not enter that incoming fit. The existing shared coefficient
preparation supplies one eddy-viscosity field to native shear relaxation
and to D3Q7 heat diffusivity through `alpha_t=nu_t/Pr_t`.
Callers supply fresh molecular transport each time, preserving the original
closure, coefficient-ratio and disabled-model admissions.

The original joint forced-thermal transaction retains force-buffer exchange,
reservoir transport, prescribed thermal-source half-kicks, finite stationary
wall storage, wall-work receipts, histories and counters. One existing SGS
workspace supplies all coefficient storage; the wrapper introduces no
per-step allocation or alternate physical state. Publication occurs only
after the original native transaction succeeds. Refusal or cancellation
leaves accepted physical fields and caller-owned endpoint forces unchanged;
scratch may change.

This serial API remains fixed-reference-capacity temperature transport:
modeled mechanical transfer and body work are not deposited as heat. It
does not add a total-energy model, wall function, pooled force/thermal owner
or empirical turbulence validation.

Three regression groups are authored for both lattices and all three
closures: constant physical velocity under spatial forcing versus a
wrong-force-frame negative control; coefficient/Prandtl agreement, native
joint-transaction equivalence, disabled controls and storage reuse; and
invalid endpoint force, late thermal failure, diffusivity-ratio refusal,
cancellation and retry. Pinned-nightly syntax/format checks and source
review passed. Native compilation and execution of this extension remain
unverified: real-graph builds exhausted memory upstream in `asupersync`
before reaching LBM. No test-pass or performance claim is made for those
unexecuted cases.

## Analytic full-gradient local SGS response

`sgs3::eddy_viscosity_tangent` returns the original Smagorinsky, WALE or
Vreman viscosity plus all nine analytic derivatives with respect to
`G[i][j]=d(u_i)/d(x_j)`. `vreman_eddy_viscosity_tangent` additionally
retains three directional filter widths, assigned to derivative coordinates.
Both return `SgsViscosityTangent { viscosity, derivative }`.
The original primal closure supplies the viscosity bit-for-bit. Fixed-size
invariant differentiation introduces no allocation or production finite
differences. WALE retains the full gradient's strain and rotation dependence;
Vreman retains compensated weighted minors and exponent-scaled arithmetic.

Derivatives have dimensions of length squared, with coefficients and filter
widths held fixed. A stress Jacobian must include BOTH terms of
`d(2 nu_t S)=2 nu_t dS+2 S(dnu_t/dG:dG)`; the returned viscosity derivative
alone is not the stress Jacobian. Overflow or underflow of an otherwise
nonzero derivative explicitly refuses instead of clipping.

At Smagorinsky zero strain, all-zero WALE gradient and Vreman zero/rank-one
states, the response uses a documented zero generalized-derivative
selection. These are not claims of classical differentiability: a retained
rank-one Vreman example has nonzero one-sided derivatives. At nonzero WALE
gradients with zero traceless symmetric square, the ordinary derivative is
zero. No new nonlinear convergence, turbulent-flow validation, wall-function
or transition claim is introduced.

Twenty native numerical tests passed on pinned nightly: ten original
closure tests plus ten tangent groups covering all nine entries against
independent centered differences, exact primal equality, isotropic rotation
covariance, parity and homogeneity, unit/filter scaling, WALE near-wall cubic
behavior, anisotropic Vreman axes/aspect ratios/rank limits, complete modeled
stress response, disabled models and unrepresentable derivatives. The build
compiled the complete original `fs_math` crate and the actual
`model.rs`, `vreman.rs` and `tangent.rs` sources without runtime
adapters, numerical stubs or substituted kernels. This establishes the
numerical slice, not full `fs-lbm` or downstream solver integration.

## Implicit conjugate heat with molecular and phase resistance

`FlowField::advance_isotropic_implicit_conjugate_viscous_phase_resistance_thermal_volume_momentum_with_baths`
combines the native full molecular stress stage, implicit spatial conduction,
finite solid storage, bath/contact exchange and phase-dependent Darcy drag on
the original material chart. Molecular inertia uses actual M/V; original checked
phase sources supply its viscosity once. Its physical heat joins authored heat
once, then the existing resistance owner iterates Darcy heat, implicit thermal
state and density-compatible volume flow from the same original fluid and solid
inventories. Only the accepted trial transports M, M*h and P and publishes a
finite-wall endpoint. Fixed-solid mass and momentum stay outside fluid routes.

Source conductivity uses the native analytic enthalpy Newton/Krylov graph.
Whole-call residual-evaluation and Krylov ceilings include rejected outer
trials. Molecular assembly, implicit heat, conservative caloric transport and
Darcy heat/constitutive defects enter the shared caloric and both mechanical
budgets exactly once. Physical molecular/Darcy heat, bath/contact transfers,
pressure/source work and numerical mechanical losses remain separate accounts.
Cancellation or any source, heat, mechanical or work-budget refusal publishes
no partial momentum, material inventory or heated wall.

G1/G3 checks independently reconstruct finite-wall/contact Fourier heat,
unequal-density latent volume expansion and combined fluid/wall/advection
mass/enthalpy balance while both physical dissipations are active. Zero Fourier
exchange preserves the prior composed molecular/resistance endpoint exactly.
G4 checks cancel after a completed heat trial, exhaust aggregate residual work,
and reproduce the unchanged accepted result on scratch reuse. Focused execution
uses actual material, native sparse, deterministic reduction, LU, FGMRES,
enthalpy/conduction and flow sources; full Cargo/execution-context integration
remains a separate verification gate.

The split retains original-state molecular viscosity, implicit pre-advection
conductivity and transported-endpoint phase resistance. It is first order and
requires refinement; it does not infer simultaneous endpoint shear feedback.
Fluid phases must remain mobile. Exact solid locking, moving skeletons and
free surfaces remain outside this owner.

## Actual-material-density SGS momentum and separated energy transfer

`FlowField::advance_material_momentum_with_viscosity_and_sgs` composes the
original source-viscosity predictor, force/pressure update and accepted
mass-route transport on closed fixed geometry. Original cell M/V supplies
every SGS dynamic-viscosity conversion and local molecular-ratio limit;
native half-cell masses still supply inertia. Source reference density does
not replace either authority. Smagorinsky, WALE and Vreman retain their shared
filter, native reciprocal-edge and explicit exploration policies.

One retained backward-Euler predictor refreshes SGS on its current native
velocity before the first solve and after every coupled sweep, before the
same-velocity residual. Original molecular viscosity is sampled once per
fluid cell and held for the interval. The added stress is deviatoric
`2 rho nu_t dev(D)`; the original molecular `2 mu D` response retains its
existing isotropic part. SGS supplies no bulk viscosity. Zero model coefficient
validates the configuration and then uses the original molecular path exactly.

Normal-cell work and each reciprocal shear-edge work use the SAME total
dynamic viscosity as momentum. At a material edge, total cell viscosities
enter harmonic transmission first. Resistance shares distribute the common
stress work, then each cell's original molecular and SGS fractions divide its
local loss. Wall traction and wall heat retain the existing molecular closure.
No independently averaged SGS stress, source-reference-density conversion or
second momentum solve is introduced.

`MaterialViscosityReport::cell_sgs_transfer_j` and `sgs_transfer_j` retain
modeled removal from resolved kinetic energy separately from molecular heat.
Only molecular heat enters the original specific-enthalpy inventory. The
existing kinetic/work audit includes both physical channels while retaining
backward-Euler loss, scatter loss and representation work separately. Source,
cap, solve, heat-insertion or cancellation refusal returns no partial candidate.

This remains a first-order split with frozen material density/viscosity during
the predictor. The modeled transfer is not an SGS energy-storage, cascade or
heating equation. Active pressure-mouth SGS remains outside this entry, as do
a variable-density total-energy qualification and resolved-turbulence accuracy
claims. Ten focused regression groups cover analytic graded-density decay,
deviatoric stress, harmonic work partition, molecular-only walls/heat, every
adjacent edge cap, zero-model identity and refusal/replay. Source parsing,
scoped formatting and independent arithmetic were checked; native execution
remains blocked before fs-lbm by the pinned asupersync dependency build.
