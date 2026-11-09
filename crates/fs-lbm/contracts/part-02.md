<!-- fs-lbm contract, part 2 of 3; a higher number is newer. Index: ../CONTRACT.md. -->
## D3Q19 sparse-sweep performance evidence model (bead 712t)

`perf` is the normative semantic substrate for the ignored GLUP/s runner. It
freezes these workload rows before any timing result is inspected:

- one 128³ sparse-grid domain with every 4³ tile active (`dense-active`);
- one 128³ sparse-grid domain whose whole-tile active set rounds upward to at
  least ten percent (`sparse-ten-percent`: 3,277 tiles / 209,728 cells);
- a serial row with exactly one worker and an all-core row retaining the exact
  normalized worker/placement identity.

The version-2 logical traffic model is source-counted, not a compulsory-miss
fantasy. One active BGK update performs two population reads and two population
writes across collide and pull-stream: `19 * 8 * 4 = 608 B/cell`. Splitting
same-tile pulls from cross-tile/domain halo pulls removes sparse-map lookup from
784 of the 1,216 links in one 4³ tile. The remaining 432 halo links are exactly
6.75 lookups/cell, adding `6.75 * 16 = 108 B/cell` of logical source key/slot
traffic and `16 / 64 = 0.25 B/cell` of amortized active-tile identity, for
exactly `716.25 B/cell`. BTree node/cache-line traffic can only increase real
traffic; it is not invented as a portable constant. The generic three-axis
BGK+Guo path counts 142 FLOPs for density/momentum/forced velocity, 271 for
equilibrium, and 611 for force projection plus relaxation, totaling 1,024
FLOPs/update and approximately 1.4293 FLOP/B. Any kernel/model change must bump
`D3Q19_PERF_MODEL_VERSION` and restate this arithmetic.

The plan's 1.0 GLUP/s Apple-M-class and 0.6 GLUP/s Threadripper-class values are
retained as informational `plan_target_glups` fields. They are never silently
used as initial gates. A positive/negative gate requires `EvidenceClass::Citable`,
an exact admission-receipt hash, at least three timed repetitions, stable
critical-path signatures, and a separately authorized positive `floor_glups`.
Before the first honest quiet-window calibration, rows are explicitly
`report_only` with a named refusal and a null floor. Contaminated pre/post axes
are `environment_invalid`, neither floor pass nor floor miss.

Per-repetition tile/halo timing forms a bounded task DAG.
`attribute_critical_path` computes the longest path in the max-plus semiring,
independent of task and predecessor enumeration. Exact path ties choose the
smaller stable task id; exact dominant-class ties choose the stable class order.
Receipts retain the source-to-sink path, makespan, per-class critical wall,
dominant class from `activation|collide|halo|stream`, and its parts-per-million
share. Citable rows require the same path and dominant class across
deterministic repetitions; varying wall values remain telemetry, while
topology/class drift fails closed.

`SparseGrid3::step_pooled_observed` is the opt-in measurement seam for that
driver. It retains the exact executor reports and whole-pass wall envelopes for
the existing collide and combined stream/halo pool invocations. The stream
kernel uses the same worker-count-independent groups as the ordinary path and
partitions its writes into two disjoint loops: same-tile pulls first, then
cross-tile/domain halo lookup and pulls. The ordinary path runs those same
loops without clocks or timing buffers; the observed path adds one clock pair
per completed group and returns canonical `(group, first_tile_slot, tiles,
local_stream_wall_ns, halo_wall_ns)` rows. Group walls overlap across workers
and therefore MUST enter max-plus task-DAG composition; summing them as serial
time is invalid. A missing group observation refuses before the private
destination is published, and cancellation retains the pre-step state exactly.

`sparse_sweep_task_samples` validates successful executor completion, exact
worker/group geometry, non-cancelled reports, and nested wall envelopes before
allocating a bounded DAG. Its stable identities encode
`optional activation -> collide -> {local-stream[g] -> halo[g]} ->
stream-residual + publication`. Activation is present only when the measured
repetition actually rebuilt the active set; fixed-active-set throughput lanes
omit it instead of charging one-time setup to every sweep. The final join
carries only the measured stream-pass wall not already on the selected group
path plus publication, so the max-plus makespan is exactly the sum of optional
activation, collide-pass, stream-pass, and publication rather than a sum of
overlapping worker walls. The join may name every bounded halo terminal; the
global task and edge caps remain the allocation authority.

### No-claim boundaries (D3Q19 performance model)

- `perf` performs no timing, machine probe, authority validation, baseline
  promotion, durable-ledger write, SIMD dispatch, scheduler change, or tuning
  publication. It validates the record the ignored runner must emit.
- Logical byte and FLOP accounting defines a portable roofline denominator; it
  does not claim measured DRAM traffic, cache behavior, instruction count, or
  target attainment on any host.
- No floor is authorized by source code alone. Both reference-ISA rows,
  pre/post admitted axes, quiet-host dispersion, retained external receipts,
  and the sealed tune-row path remain required measurement evidence.
- The observed sweep does not itself time active-set construction,
  authenticate machine axes, emit a GLUP/s row, or write the ledger. Its task
  lowering consumes an optional caller-owned activation wall and exposes
  truthful primitive max-plus inputs for the ignored driver; neither surface
  is itself citable performance evidence.

### Ignored D3Q19 throughput driver

`tests/perf_lane.rs` is the opt-in hardware driver for bead `712t`. It measures
the existing sparse-grid sweep at the fixed 128³ memory-resident shape for the
dense-active and smallest whole-tile ten-percent-active sets, each with one
worker and all reported logical workers. The active set is a canonical Morton
prefix, the executor seed and complete placement identity are retained, one
warmup precedes five timed repetitions, and each timed observation is lowered
with no fabricated activation task. The driver emits the raw fs-roofline
attainment, every per-repetition critical path, the typed aggregate GLUP/s row,
and pre/post machine-axis admission receipts through fs-obs.

The driver applies the same bounded baseline-store, firmware, promotion-policy,
and retained-source-receipt checks as the existing fs-feec/fs-fft lanes. Axis
reprobe contamination is environment-invalid even when baseline configuration
is absent. An authority-admitted baseline is still downgraded to report-only:
no LBM anti-collapse floor has been calibrated and the centralized external
gate recorder does not yet own an LBM lane. Therefore this driver can collect
the first honest rows, but cannot mint a citable pass or failure. Both reference
ISA rows, floor authorization, durable LBM ledger recording, and sealed tile
tuning remain explicit follow-on acceptance work.

## Density-proportional acceleration on the sparse free surface

`FreeSurface3::rebind_body_acceleration` selects a held lattice acceleration
with actual cell force `F=rho*a`. The complete accepted force history travels
with the original populations. Replacing the acceleration preserves canonical
density, physical velocity and force-corrected stress without advancing time.
During a native step, collision uses the initial force and the streamed endpoint
uses its new density times the same acceleration, without a fictitious reframe.
This retains the trapezoidal source impulse of the original Guo convention.

Missing gas-link populations retain the existing explicit physical boundary
velocity and solve their density-dependent even source through one scalar
linear closure. New interfaces enter their own density-dependent source frame.
Existing force entries are updated in place; there is no second fluid solver or
iterative pressure solve. Normal and pressure-loaded steps are transactional in
this mode. Independently bound local forces and LegacyGlobal replay refuse;
the historical reference-density constructors keep their original behavior.
The original low-Mach, capillary, topology and work-budget limitations remain.
This adds no liquid EOS, thermal expansion or phase-density volume source.

The `body_acceleration` integration regressions specify equal acceleration at
unequal densities, load switching, an independent affine-density endpoint
impulse and discrete atmosphere/wall hydrostatics in all six directions.
Focused native tests also cover physical stress through rebinding, gas-pressure
moments, interface birth, continuing mass/strain/motion operations and refusal.

## Material-bound isothermal free surface

`material::free_surface::plan_free_surface_from_material` lowers a frozen
`fs_material::liquid::surface::ResolvedLiquidSurfaceState` into the existing
`d3q19::freesurface3::FreeSurface3`. The source bundle requests density and the
explicit dynamic/kinematic viscosity basis, plus surface tension only for
`LiquidCapillarity::Sourced`. Thermal conductivity and heat capacity are not
invented for an isothermal momentum solve. The complete source point and usage
receipts retain temperature, pressure, composition and surface-environment
conditions, including stated unknowns. An already resolved multiphysics bundle
can be adopted intact; the surface consumer does not query independently or
drop other properties. Constants, supported temperature coefficients and
explicitly typed held composition/moisture/concentration curves or surfaces
are admitted by the original material owner. That owner retains scalar/T-only
interface selectors and requires typed semantics for new environmental laws.
Pressure response and dynamic frequency/rate/strain dependence still need
their evolving owners. The frozen scenario supplies no composition transport,
mixing closure or live change of reference density.

The existing material momentum-scaling helper owns nu_lattice and tau for both
thermal and free-surface consumers. Explicit dx/dt convert gravity to native
acceleration and sigma to sigma/(rho0*(dx/dt)^2*dx). Construction
bounds the full possible voxel domain, validates tile coordinates and refuses
initial liquid inside a permanently solid tile. The source-bound wrapper returns
native ledger, regional mass and unregioned carry in kg, and wet-cell density,
Guo-corrected velocity and incompressible pressure perturbation in SI. Gas
populations are never presented as liquid samples.

The native `rho_lattice*g_lattice` law is the default for both binary and
supplied-state sparse material constructors, including supplied pressure.
`MaterialFreeSurfaceScaling::with_density_dependent_gravity` explicitly selects
that same default. Binary initialization starts at physical rest in the selected
source frame; supplied velocities are retained.
`MaterialFreeSurface3::rebind_gravity` replaces this held SI acceleration without
instantaneous impulse, retaining the original material, pressure reference and
dx/dt units. All components are admitted before native publication; nonzero
loads that disappear in conversion refuse. The current numerical density sets
gravity, while sourced reference density still owns the SI mass/pressure units.
The native owner retains a private material-gravity descriptor. Existing thermal,
Marangoni and Lorentz owners explicitly bind their complete force field with
current-density gravity. Nonzero extra forces retain the original first-order
held-force split and endpoint Guo reframe; zero extra force keeps or restores
native density-aware gravity, including its endpoint pressure closure. Source
Boussinesq acceleration is `rho(T)/rho0*g` at any numerical pressure density,
while additional Lorentz and surface loads remain force densities. Complete
rigid-frame forcing takes over material gravity once, using current density
for its complete acceleration. These transfers preserve physical velocity,
stress and extensive inventories without an extra impulse or heat source.
Independently authored native/affine acceleration and LegacyGlobal replay still
refuse incompatible composition. Native handoff regressions cover unequal
densities, zero-extra restoration, force moments, cancellation and ownership;
the thermal regression varies pressure density independently of source rho(T).
The SI regressions
specify equal acceleration across source and pressure-density changes, refreshed
gas-pressure trajectories, no-impulse/atomic load replacement, and zero-load
compatibility. This does not supply a compressible or thermally expanding liquid.

The capillary reconstruction must remain positive at both signs of the native
smoothed-fill curvature clamp; its identical bound is exposed as
`FreeSurface3::CURVATURE_LIMIT`. This is a necessary numerical gate, not a
capillary accuracy or stability certificate. The existing neutral/wetting ghosts
remain explicit numerical scenario choices; no measured contact angle, surface
energy balance, Marangoni force, gas evolution, thermal update, phase transfer or
moving-wall law is inferred. Lattice sound speed is not physical liquid/gas sound
speed. `tests/material_free_surface.rs` covers dimensional Young-Laplace
bookkeeping, material similitude, separate viscosity/capillary substitutions in
real native evolution, mass accounting, Guo samples and domain refusals.

An independent conservative explicit capillary-wave admission now requires
`dt <= sqrt(rho0*dx^3/(4*pi*sigma))`, equivalently native sigma at most
`EXPLICIT_CAPILLARY_LATTICE_LIMIT = 1/(4*pi)`. The pressure-only outer phase
does not contribute gas inertia. This bounds the combined travel of two
counterpropagating shortest resolved waves (`lambda = 2*dx`) to one cell per
timestep; it remains a conservative wave-resolution policy, not a proof of
complete LBM/VOF stability or convergence. Positive reconstructed density and
zero advective Mach alone do not satisfy it. The dense D2Q9 constructor and
every public-sigma step, sparse D3Q19 construction, physical source scaling,
and local/extrapolated thermocapillary source updates share the same bound.
The native owners additionally require `sigma_lattice/rho_lattice <= 1/(4*pi)`
at each actual wet-cell density; the reference-density coefficient ceiling
remains in force. Liquid density, not tracked mass or fill fraction, sets this
inertia. Transferred initial states, mutable dense states, local source fields,
height-column averages and local interface extrapolations all enforce it.
Streaming checks both the current and candidate endpoint density before
publishing any populations or tracked mass, so prescribed pressure rarefaction
cannot silently leave the accepted liquid outside its capillary-wave domain.
No collision/curvature/pressure evolution is changed for admitted states.
The helper `maximum_capillary_timestep_s` exposes the physical maximum before
constructing a plan; `capillary_timestep_ratio` reports its consumed fraction.
Explicitly neglected tension has no capillary limit. Source-temperature or
interface extrapolation increases that exceed the bound refuse before accepted
fluid/thermal state publication. The focused G0/G1 tests check independent
wave dispersion, dx^(3/2) and material scaling, positive-pressure refusal,
zero-tension behavior and unchanged state after rejected coefficient updates.
Nonunit-density and partial-fill regressions exercise the native constructors;
independent planar pressure reconstruction (`rho_next = 5/6`) exercises atomic
2-D/3-D rarefaction refusal. Thermocapillary regressions retain the old source
map when local density or a column average makes the new coefficient unsafe.
These checks do not establish measured capillary-wave convergence.

### Prescribed gas pressure on the native free surface

`FreeSurface3::step_with_gas_pressure` adds an explicitly supplied gas-pressure
perturbation to the original capillary anti-population reconstruction:
`rho_boundary = 1 + (p_gas_gauge + sigma*kappa)/cs²`. One callback sample per
exposed liquid cell supplies every missing gas population at that cell. It is
sampled anew each step, including newly exposed interfaces; no stale map or
second momentum integrator is retained. Zero load preserves ordinary evolution.
The combined density excursion must fit a finite caller bound in (0,1).
Pressure refusal and observed cancellation precede population/mass publication.

The default version-4 boundary uses the physical Guo velocity
`u = (sum(e_i f_i) + F/2)/rho` and the FULL local force. For a reconstructed
opposite-link pair, its stored even moment is
`feq_i(rho_boundary,u) + feq_opp(rho_boundary,u) - (S_i(u,F)+S_opp(u,F))/2`.
This follows the native population source transformation; it is not a new
published boundary algorithm. It removes the artificial force-squared stress
introduced by putting raw momentum/rho into equilibrium, and makes pressure
independent of the chosen uniform/local force split. Explicit `LegacyGlobal`
replay retains its historical pressure and conversion arithmetic. Force-free
reconstruction and the missing-gas-links-only choice are unchanged. Focused
regressions require planar halfway hydrostatic columns to stay at rest for 64
steps in all six gravity orientations, with the even Guo moment and equivalent
force frames checked separately.
This does not repair the general free-interface viscous-stress, geometry or
transient force-balance limitations of the underlying reconstruction.
For the forcing convention see Guo, Zheng and Shi (2002),
<https://doi.org/10.1103/PhysRevE.65.046308>; for free-surface force-balance
limitations see Schwarzmeier and Rüde (2023),
<https://doi.org/10.1002/fld.5173>.

`MaterialFreeSurfaceScaling::with_pressure_reference` fixes the absolute Pascal
datum before constructing a trajectory. Its material-bound `step_with_gas_pressure`
accepts absolute pressures at physical cell centers and lowers only their
difference from that datum through the original rho0*(dx/dt)² scale.
`absolute_pressure_si` and the existing perturbation samples share this datum.
The load reuses original source viscosity, capillarity and mass-conversion
routes. It does not evolve an external gas, refresh the frozen liquid source,
account pressure work, or turn artificial lattice compressibility into a
physical liquid/gas equation of state. Tests cover actual direction/reversal on
all axes, pressure-reference invariance, material-density similitude, capillary
zero-load equivalence, mass preservation and atomic refusal/cancellation.

The dense D2Q9 `FreeSurface` now uses the same pressure law through
`step_with_gas_pressure` and `step_with_gas_link_pressure`. Cell fields are
sampled once per exposed cell; link fields independently load each missing
gas population. Curvature and equal-density equilibria are cached per cell.
The Guo half-source correction applies to ordinary and prescribed-pressure
reconstruction, removing spurious moving-force stress from the dense boundary.
One complete native candidate supplies atomic pressure/geometry/cancellation
refusal while preserving the original mass and local conversion owner.
The pressure and conservative-momentum wrappers share one collision/streaming
kernel, including the same Guo correction. Pressure loading retains the legacy
local conversion policy; `step_with_momentum` retains its separately admitted
capacity-bounded momentum conversion. Neither wrapper silently selects the
other's conversion policy.
`MaterialFreeSurface2` exposes the same absolute-Pascal adapters and retains
its mass-per-depth convention. Physical links carry both centres and their
midpoint, without inferring a surface normal. G1/G3/G4 regression definitions
cover loaded static equilibrium, independent link traction, forced even stress, zero-load
replay, gauge/density similarity, mass closure and late failure rollback.
These are implementation checks, not measured free-interface validation;
the isothermal single-liquid and pressure-work exclusions above still apply.

`MaterialFreeSurfaceScaling::free_surface3_from_cells_with_pressure` also
admits a positive absolute liquid-pressure field at occupied initial cell
centers, using that same fixed reference and pressure scale. The caller can
initialize hydrostatic or capillary pressure together with transferred fill
and Guo-corrected velocity; the pressure is sampled only at construction.
Every initial density must satisfy an explicit excursion bound in (0,1).
Native mass remains `fill*rho_lattice*rho_source*dx³`, so pressure-induced
numerical density contributes to initial mass; no physical equation of state
or mass-preserving pressure projection is asserted. Frozen source properties
are unchanged. The ordinary zero-perturbation initializer is bit-preserved.
Tests check hydrostatic cancellation of initial interior acceleration in all
six gravity directions, a pressure-matched flat surface at rest, the supplied
spherical Laplace jump, gauge invariance and refused/cancelled ingress. These
are initialization checks, not an exactly well-balanced hydrostatic or
static-drop scheme, and do not establish long-time or refinement accuracy.

### Sourced equilibrium angle at one planar domain wall

`MaterialFreeSurfaceScaling::with_planar_wetting` binds a
`ResolvedLiquidWettingState` to `PlanarWall3` before constructing flow.
It requires the same liquid phase/process/revision, complete query and quantity
descriptors as the existing surface bundle; only the angle comes from this
ordered substrate/liquid interface. The original source/scaling owner continues
to supply rho, viscosity and sigma. The caller assigns substrate and named
environment to the receiving geometry; this does not resolve a gas EOS.

The same binding now reaches the dense D2Q9 owner through `free_surface2`.
`MaterialFreeSurfaceDomain2` supplies an explicitly bounded, stationary closed
rectangle including its wall rows. X/Y faces map to `freesurface::PlanarWall2`;
Z faces and a conflicting legacy Wetting ghost refuse before classification or
allocation. The source still owns rho, viscosity and sigma; the geometric
boundary can consume either a direct angle or the existing Young–Dupré
work-of-adhesion resolution. The original liquid/query consistency checks are
unchanged. Dense gravity retains its native density-dependent rho*g convention.

`FreeSurface::set_planar_contact_angle` extends frozen fill before smoothing
and extends smoothed fill through the selected complete wall row and the second
ghost layer used by curvature. This preserves a planar fill field at the wall
instead of introducing curvature by truncating its smoothing stencil.
With inward wall normal it imposes `dphi/dn_wall=-|grad_t(phi)|cot(theta)` without
clipping that derivative or changing tracked mass. Other walls remain neutral;
corner intersections do not acquire inferred contact angles. Exactly pi/2
preserves the original dense neutral trajectory bit for bit. Invalid angles,
periodic normal directions, incomplete wall rows or missing inward stencils
refuse assignment; changing the admitted topology or fallback model later
panics before native step mutation, as documented by that dense API.

`MaterialFreeSurface2` retains the original native evolution and read-only SI
sampling. Its mass ledger and carry are **kg per metre of out-of-plane depth**,
using rho*dx², rather than treating a 2-D cell as a 3-D volume. This is a
translationally invariant section, with no out-of-plane curvature or friction.
The dense owner still has its original panic/conversion semantics and no
transactional cancellation API. G1 planar-gradient/zero-curvature checks cover
all four orientations; G3 source-to-pressure trajectory, mass conservation,
replay, exact neutrality and pre-allocation refusal checks reuse
`tests/material_free_surface3/adhesion.rs`. They establish this bounded
equilibrium closure, not contact-line kinetics, calibrated sessile-drop
accuracy, a surface-energy conservation law or a completed MR26 delivery.

`FreeSurface3::set_planar_contact_angle` imposes the equilibrium angle through
liquid at that one domain face. With wall normal pointing into fluid, its
smoothed-fill ghost satisfies `dphi/dn_wall = -|grad_t(phi)| cot(theta)`.
The ghost changes the existing curvature and actual gas-side pressure
reconstruction; it never edits tracked liquid mass. Exactly pi/2 takes the
existing Neutral path bit for bit. All other walls remain neutral, including
internal tile solids; arbitrary staircase normals are not inferred.
Combining the sourced planar angle with the legacy Wetting ghost refuses before
geometry or allocation; the remaining walls must be explicitly Neutral.
Unrepresentable near-endpoint stencil magnitudes refuse before assignment.

`step_with_gate` reuses the native step with tile-level cancellation checks
and a final check before population publication. Observed cancellation preserves
populations, mass, carry and conversion counters. After publication begins the
original conversion cascade completes. The derived fill cache refreshes on retry.
The original no-argument `step`, sampling and mass APIs remain available.

G1 manufactured checks cover all six wall-normal conventions, zero planar
curvature and spherical Laplace-curvature refinement. The focused integration
tests in `tests/material_free_surface3.rs` exercise source-angle trajectory
changes, exact 90-degree neutrality, omitted-capillarity equivalence, mass
conservation, source/state/unit/dynamic-support refusal, cancellation and retry.
These tests establish coupling and operator geometry; they do not establish
calibrated sessile-drop dynamics, hysteresis, thermocapillary stress, capillary
oscillation/spurious-current convergence or a closed surface-energy account.

### Original liquid sources on solved-flow enthalpy transport

`conjugate::enthalpy::advection::liquid` adds
`MaterialChtFlow<ResolvedLiquidState>::advect_liquid_enthalpy` to the existing
native flow-to-enthalpy owner. The same solved face fluxes, fixed control-volume
map, donor-cell transport, combined CFL, energy audit and thermal timeline own
the evolution. The adapter retains the exact resolved reference state used by
momentum and the original borrowed caloric chart used by thermal storage.

Admission requires common material-card identity, identical selected Cp/density
claims, and all non-temperature coordinates and quantity descriptors. Density
and viscosity must be validity-wide scalar sources under the existing
instantaneous Newtonian admission. Cp and conductivity use the same shared
support guard as the existing surface/momentum coefficients: non-temperature
response curves and unevolved dynamic source bounds refuse. Every knot in the
complete caloric chart
must remain fully liquid at the fixed source density. Momentum and thermal eddy
viscosities must be zero. A matching card identity or equal endpoint densities
alone do not establish a compatible material binding.

The original pinned liquid law is re-evaluated at actual thermal source epochs,
endpoints and events. Each accepted endpoint retains its source receipts and
caloric state. Leaving source support or refusing an inner physical-history step
discards the complete thermal/advection/history candidate. Conduction and other
physical histories remain with the composed inner observer. This is passive
sensible-enthalpy transport on frozen laminar momentum; it supplies no phase
mobility, variable-density volume flux, viscous or pressure work, buoyancy,
turbulent heat closure or free-surface motion. The native strict local mass
balance gate remains unchanged and can require further momentum refinement.

The seven regression sources exercise the actual native duct and internal-fan
solves, independent upwind energy accounting, a viscosity-only material swap,
original-claim support loss, candidate rollback, and context/selection/object
identity, complete-chart and dynamic thermal-source refusals.

### Local liquid viscosity and molecular conduction on native substeps

`LiquidFlowEnthalpyInput::feedback` binds the original pinned liquid/caloric
sources to `MaterialFlowFeedback`. Each actual thermal source epoch supplies
the native full-stress laminar MAC viscosity field and conservative face mass
fluxes. Momentum is quasisteady; density and reference cell mass stay fixed.
Thermal storage, advection, source histories and energy accounting keep their
existing owners.

`MaterialFlowFeedback::with_fluid_conduction` optionally assigns existing full
fluid-face thermal link IDs. The same applied material samples supply k(T);
the shared thermal series-conductance helper reduces both native half-cell
resistances. The original pairwise link exchange and combined CFL use those
values, without adding sources or replacing authored contacts. Both patch
coverages and total cell surfaces must match native voxel geometry. Parallel
links on an owned face, competing coefficient owners, stale h/time epochs and
changed geometry refuse the whole candidate. Missing k never implies insulation.

G1/G3 source-only conductivity and viscosity swaps exercise actual native
momentum, heat propagation, graded-cell resistance, cell mass and boundary
energy. G0/G4 ownership/source/geometry refusals preserve accepted state.
This fixed-density, scalar Newtonian bulk-fluid model supplies no additional
contact resistance, phase mobility, pressure/viscous heat or moving surface.

### Material-density natural convection on the enthalpy timeline

`material::conjugate::MaterialChtFluid::solve_natural_convection` also exposes
the shared source-driven Boussinesq projection directly to steady CHT consumers,
without requiring a transient enthalpy timeline. One local source query supplies
`(rho(T)-rho0)*g`, dynamic viscosity and molecular conductivity. The existing
full-symmetric-stress SIMPLEC owner and conductivity-field energy owner solve
on the same original domain and conservative volume fluxes. Reference rho/Cp
stay in inertia, continuity and sensible heat-capacity advection; the supplied
density/Cp departure bounds apply at every trial and proposed endpoint.

Acceptance requires converged native flow/heat solves plus independent
unrelaxed temperature, relative mu/k and absolute force-defect tolerances.
Underrelaxation cannot shrink the defect used for acceptance. Returned fields
retain both the actual applied samples and fresh endpoint samples. Negative or
nonmonotone thermal expansion follows the source density law without an
independent beta. Fixed solid conductivities, ordinary contacts and authored
thermal boundaries are retained; stateful solid tables, LVEL and independent
eddy-conductivity fields refuse. This is a numerical steady Boussinesq model,
not variable-density mass transport, transient dynamics or experimental
validation. Source queries must be repeatable; budget, domain or cancellation
failure publishes no partial solution. Focused G1/G3/G4 cases in
`tests/material_conjugate.rs` cover nonlinear heat balance, Darcy material swaps,
negative expansion, endpoint source refusal and independent force convergence.
Low-conductivity warm and fan-driven starts keep the caller's initial field;
conductive initialization is reserved for unanchored, exactly zero-flow startup.

The same `MaterialChtFluid::march_natural_convection` now supplies local source
density, viscosity and conductivity to the original native unsteady CHT march.
Each accepted temperature is sampled once per fluid cell; the three coefficient
buffers are reused, and the last applied and endpoint source epochs and times
are retained. Reference rho/Cp still own inertia, continuity, sensible storage
and heat-capacity advection. Source support and density/Cp bounds are checked at
the initial state and every endpoint, including the final returned temperature.
A refused source endpoint or cancellation returns no partial solution.

This path retains MAC velocity, pressure and BE/BDF2 history. Local molecular
viscosity also refreshes geometric Darcy resistance. The exact original
backward-Euler thermal storage, fixed-cell reservoir account and converged
endpoint radiation remain in use with the admitted conductivity field. Solids
retain explicit fixed heat capacity and grid-aligned conductivity; scalar beta,
SGS/LVEL, independently authored eddy conductivity and dynamic solid tables
cannot compete with the material owner. Material/thermal coupling is FIRST
ORDER even when momentum uses BDF2; time-step refinement remains necessary.
No variable-density enthalpy, viscous-heating or gravitational-storage account
is implied. Added `tests/material_unsteady_buoyancy.rs` cases cover independent
nonlinear-conductivity storage and signed Darcy BE/BDF2 recurrences, unchanged
radiation/reservoir accounting, and final source/departure/cancellation refusal.
These are regression definitions; adding them does not claim numerical execution.

`LiquidFlowEnthalpyInput::feedback_boussinesq` admits the same pinned liquid
and caloric sources for bounded density variation. It delegates to
`MaterialFlowFeedback::try_new_boussinesq` and the existing `BoussinesqMaterial`
projection. One source evaluation supplies actual rho for `(rho-rho0)*g` and
the same mu/k for native stress/conduction. Reference rho0 stays in inertia,
fixed cell mass and the original conservative face-mass graph; caloric density
metadata cannot reassign mass. The original h/T law and its source pins/context
remain authoritative. Pressure openings use the removed reference hydrostatics.
No independent expansion coefficient or duplicate buoyancy load participates.

`simple_flow_with_viscosity_and_force` applies finite cell force densities to
the same staggered momentum operator, requiring zero force on solid voxels.
Source/endpoint density and Cp departures must satisfy explicit bounds; native
convergence, energy, CFL and all-or-nothing history owners remain unchanged.
Applied epochs retain the exact buoyancy field. This is quasisteady Boussinesq
natural/mixed convection, not variable-density continuity or a total-energy
model: gravitational/kinetic storage, pressure/viscous heat, phase mobility and
free surfaces remain outside this path. G1/G3 source-density swap regressions
exercise a heated porous column without a fan and independent Darcy/heat
balances. G4 endpoint-bound and binding refusals preserve accepted history.


### Conservative initial fill and physical velocity

`FreeSurface3::from_cells` and
`MaterialFreeSurfaceScaling::free_surface3_from_cells` initialize the original
sparse VOF from supplied cell fills and physical velocities. They preserve
initial mass and Guo-corrected momentum to floating-point tolerance instead of
creating a half-filled resting shell. Full wet boundary cells retain fill one;
fractional interior cells remain interfaces. Deterministic Morton/lane
initialization refuses duplicate, invalid, out-of-domain or wall cells. The SI
adapter retains the original sourced rho/mu/sigma and planar-wetting binding,
uses reference density one, and checks the admitted speed and domain budgets.

Construction uses the same `Cx` and returns a new candidate; it never resets an
accepted trajectory. Positive populations and finite conserved moments are
necessary admission checks, not a stability certificate. Equilibrium population
projection does not transfer non-equilibrium stress, subcell moments or caloric
state, and does not establish subsequent momentum or surface-energy conservation.
The focused material-free-surface regressions exercise material substitution,
Guo velocity, exact initial fill, order independence, a native step, and refusal
and cancellation. Their presence does not imply an executed native test pass.


### Enthalpy-dependent mushy mobility in native flow

`MaterialFlowFeedback::with_mushy_mobility` connects the original equilibrium
h chart and declared phase closure to the existing momentum/enthalpy timeline.
It requires the exact registered caloric object and equal constant phase
densities, multiplies inverse permeability by the same resolved local viscosity,
and applies that resistance implicitly through the existing MAC/SIMPLE momentum
and pressure correction. The original face fluxes transport total h once,
including latent energy. Applied drag is retained in `MaterialFlowEpoch`, and
source/endpoint support is checked before accepting the coupled candidate.
The existing liquid-only ingress remains unchanged; the explicit mixed-phase
attachment uses the existing generic `FlowEnthalpyInput` path.

This is bounded quasisteady, homogeneous single-velocity enthalpy/porosity flow
in a stationary drag frame. It supplies no phase-volume source, exact solid
locking, moving skeleton/free surface or drag-dissipation heating. Permeability
parameters are declared model data. Tests cover isothermal latent-plateau flow,
data-only permeability changes, independent series-Darcy/enthalpy balances,
source retention and transactional support refusal; native execution remains a
separate verification requirement.


### Source-bound orthotropic porous momentum and heat

`MaterialFlowFeedback::with_porous_materials` binds one original
`PermeabilityTemperatureLaw` per declared porous fluid voxel. Clear-fluid voxels
are explicit `None` entries. Each source epoch samples the skeleton's intrinsic
principal permeability in m2 at the same local caloric temperature used by the
fluid transport source; all other original material conditions, keys and claim
receipts stay pinned. The existing `IntrinsicPermeability` owns the material
frame. The complete `B = mu R diag(1/K_i) R^T` enters the native world/grid
coordinates, including oblique frames and signed cross-coordinate resistance.
An isotropic material is independent of its arbitrary frame. Positive finite
principal permeability is required. Zero permeability needs an exact velocity
constraint and is not replaced by a large finite resistance.

The actual source coefficients enter the original MAC momentum rows and
pressure correction. Grid-diagonal inputs retain the exact original half-cell
arithmetic. If P averages each cell's two faces per component, cross terms use
the symmetric native action `P^T V (B-diag(B)) P`. Its local power plus the
original diagonal half-face power is
`V [ubar^T B ubar + sum_a B_aa (u_a+ - u_a-)^2/4]`, nonnegative for SPD B.
This retains face coercivity instead of introducing a checkerboard null mode
through cell-centred drag alone. The symmetric coefficients are compiled with
the source epoch; no material query or tensor rotation occurs in an iteration.
Scalar mushy resistance adds once to the diagonal; its existing API is retained.

The existing bounded component iteration resolves the cross coupling. Before
accepting convergence, all three coupled momentum residuals are reassembled on
one frozen velocity/pressure iterate. The final tight pressure projection must
also preserve that residual bound. Failure refuses instead of reporting a
converged solve with stale cross-coordinate forces. This is not an unconditional
iteration-convergence claim for arbitrarily ill-conditioned anisotropy.
Cells with nonzero cross terms may touch closed normal walls/symmetry, but
cannot touch pressure/inlet openings: the existing opening rows mirror an
exterior staggered volume, which requires a separate open-port work reduction.
A clear-fluid layer between a porous insert and those openings is admitted.

With quasisteady or unsteady heating, the SAME native action supplies the
original cell's Darcy heat once; a cross term alone may be negative and is not
a separate dissipation channel. Independent face-work and cell-heat traversals
must agree. Native kinetic energy and backward-Euler numerical loss retain
their separate accounts. Without a heating attachment this changes momentum
only. `MaterialFlowEpoch::darcy_resistance_pa_s_m2` exposes the full total matrix
and source-start skeleton receipts; `endpoint_porous_samples` retains accepted
endpoint sources.

This closure has a fixed, zero-velocity rigid skeleton and fixed geometry,
linear Darcy drag and the existing density/caloric limits. Fluid and skeleton
share the registered local temperature; that caloric chart must represent the
intended participating heat capacity. No second skeleton heat/mass store,
deformation, saturation, gas-slip closure or moving-frame work is inferred.
Bindings attach once before the first advance. Source, frame, coefficient,
endpoint and cancellation checks remain inside the existing coupled candidate
transaction, including native momentum, heat, clocks and inner history.

Focused tests in `liquid_flow_enthalpy/viscous/porous.rs` independently reduce
the four-face MAC equations to check directional drag and heat, signed axis
permutations, exact compatibility with the scalar native APIs, real
two-interval `mu(T)/K(T)` feedback, and source-domain refusal with retry.
The `porous/rotated.rs` cases independently reduce an alternating-frame 2x2
MAC circulation in all three coordinate planes. They check signed off-diagonal
motion, local and total heat, source feedback, native backward-Euler kinetic
and numerical losses, and transactional source/cancellation failure.
The supported-freezing regression also requires the original porous receipts
to match the final temperature after kinetic capture heat. These numerical
fixtures are authored regressions, not experimental material validation; native
execution remains a separate verification requirement.

### Source-bound calibrated porous inertia and heat

The original `with_porous_materials` binding also consumes an explicitly
attached `PermeabilityTemperatureLaw::with_forchheimer` source. Its three beta
claims share the original complete material bundle, fixed conditions and
temperature continuation with permeability. The material convention is
`rho beta_i |u_i| u_i`, without a factor one half. The native adapter passes
`rho beta_i` directly to its existing quadratic MAC momentum coefficient;
the older authored `FlowResistance::Volume` convention uses `C_i = 2 beta_i`.
No coefficient is inferred from permeability, porosity or particle diameter.

This is a calibrated separable material-principal-direction closure, not the
rotation-invariant vector law `rho beta |u| u`. The original orthonormal frame
is retained for arbitrary rotations, including equal principal beta values.
`MaterialFlowEpoch::forchheimer_resistance` exposes optional native
`ForchheimerResistance` values with the principal rho*beta coefficients and
their material axes; clear or immobile cells are explicit None entries. This
nonlinear law is not represented as a constant world-frame resistance matrix.

For each voxel, form the eight velocity tuples `u=(u_x±,u_y±,u_z±)` from its
SAME six MAC normal face values. The convex potential is
`Phi = V/8 sum_tuples sum_i (rho beta_i/3) |r_i.u|^3`.
Differentiating this potential supplies the original face momentum forces;
their physical work is `3 Phi >= 0`. This positive cubature retains face
coercivity instead of discarding opposing face velocities in a cell average.
Each component assembly freezes the current principal speeds, retains the
nonnegative Picard diagonal implicitly and puts cross-coordinate forces in
the existing RHS. Per-cell face rows reuse allocated storage and are computed
once per component assembly; no constitutive source query occurs in a solve.

Signed grid-axis permutations retain the EXACT original half-cell force and
heat arithmetic `sum_a V/2 rho beta_a (|u_a-|^3 + |u_a+|^3)`, with no oblique
row storage or cubature work. Exactly zero beta preserves the original linear
Darcy path without nonlinear storage or extra residual probes. Each source's
maximum superficial principal-component speed is checked on all eight tuples
at source/accepted endpoints, including retained transient motion. An admitted
grid face speed or cancelled cell average does not waive this principal limit.

An independent face-force traversal checks local cubic heat against actual
nonlinear momentum work. `ViscousPower`
exposes `forchheimer_heat_w` and `cell_forchheimer_heat_w`; both steady and
unsteady heat audits expose `forchheimer_heat_j`. These are subaccounts of the
existing total physical heat, not additional thermal inputs. Linear Darcy,
molecular viscosity and backward-Euler numerical loss keep their distinct
meanings. Nonlinear solves check the complete residual at the final frozen
velocity, including after pressure projection, before accepting motion/heat.

All updates use the original caloric storage and candidate transaction.
Source, speed, frame, arithmetic, convergence or cancellation failure publishes
neither momentum nor heat. The fixed, zero-velocity, locally thermally
equilibrated skeleton, admitted constant density, creeping-flow heat ceiling
and existing opening restrictions remain. An oblique nonlinear cell must have
a clear-fluid layer before pressure/inlet/outlet openings because their native
rows mirror an exterior staggered volume. This also applies when the linear
permeability is isotropic. Closed normal walls and symmetry remain supported.
This does not establish turbulent,
compressible, moving-skeleton or uncalibrated high-Reynolds-number flow.
The focused material-flow fixtures check source-only beta substitution,
direction permutation, signed quadratic circulation, exact zero-beta recovery,
temperature feedback, physical/numerical energy partition and rollback.
Rotated nonlinear fixtures independently reduce eight-tuple circulation in
all three coordinate planes, and exercise principal-speed and open-port
refusal without changing accepted heat or motion.

### Source-bound isotropic vector-speed porous inertia

`PermeabilityTemperatureLaw::with_isotropic_forchheimer` selects the distinct
force law `-rho beta |u| u`, with Euclidean superficial speed. The same original
material/temperature binding and momentum/heat transaction consume its scalar
beta; an anisotropic Darcy tensor remains independent. Neither equal
directional beta entries nor isotropic permeability select this convention.
`ForchheimerResistance::try_isotropic` retains the scalar coefficient without
a material frame. Its scalar accessor and the optional principal/frame
accessors identify the selected law without inventing directional coefficients.

The same eight native face tuples use the convex potential
`Phi = V/8 sum_tuples rho beta |u_tuple|^3/3`. Each momentum row freezes the
current full speed in a nonnegative implicit diagonal; all components still
participate in that speed. Independent face-force and cell-power traversals
close `sum face_force*face_velocity = 3 Phi`. This physical heat enters the
original caloric store once, with Darcy loss and backward-Euler numerical
loss retaining their separate meanings. The actual maximum tuple norm is
checked against the source's Euclidean-speed limit. Opposing face values
cannot cancel that limit or the dissipated heat through a cell average.

Exactly zero beta takes the unchanged Darcy path. Positive vector-speed cells
retain the existing coupled residual checks and require clear fluid before
inlet/pressure openings; the stationary skeleton, density and heat-model
limits above still apply. Focused regressions cover arbitrary uniform-vector
rotation, an independent potential derivative, native circulation in all
coordinate planes, source-only beta changes, exact zero-beta recovery,
temperature feedback, physical/numerical work balance and atomic full-speed
refusal. These are discrete synthetic oracles, not experimental validation.

### Tetrahedral liquid geometry and momentum enter native VOF

`material::free_surface::transfer` connects the same resolved liquid scaling to
`fs_mesh::overlap` and the native supplied-state initializer. It intersects each
source tetrahedron with actual Cartesian cells, integrates its P1 velocity field,
and uses the existing `fs_la` LU/SVD owners to reconstruct momentum.
Mass-weighted translation and rotation corrections preserve total linear and
angular momentum at receiving cell centres, subject to explicit error bounds.
The final checks read actual native fill, density and Guo-corrected velocity.

The angular solve uses the represented range of the receiving inertia. A
collinear cell row can carry both transverse angular-momentum components;
only unresolved spin about its line requires refinement. The fixed 3-by-3
native SVD is scaled before rank detection and its result is independently
reapplied to the original inertia with the caller's physical angular budget.
This is the same mass-weighted least-change correction, with no null-space
spin invented. The existing per-fragment and global projection-energy gates
still bound either sign of numerical kinetic change. Focused tests cover
straight and oblique rows, unsupported spin, full-rank recovery and unit
rescaling; the thin-prism integration case independently bounds its required
kinetic projection change and refuses a tighter budget.

Source mass, first mass moment, linear/angular momentum and kinetic energy use
analytic simplex integrals. Signed kinetic and uniform-gravity potential changes
from projection are separately reported numerical errors, never physical heat.
Source coverage, target overfill, mass/momentum, centroid displacement, projection
energy, native speed, work and occupied-cell caps all gate publication. Unsupported
spin on a rank-deficient receiving grid requires refinement. Original Cx polls
cover geometry, integration, correction, initialization and final publication.

This is an already-liquid, uniform-material isothermal transfer without stored
elastic/plastic history. Its P1 velocity field supplies the initial viscous strain.
Source tetrahedra must be nonoverlapping; target overfill is not a full mesh
intersection certificate. Capillarity must be explicitly neglected because a
surface-energy projection is not represented. There is no enthalpy advection,
solid-to-liquid history release, initial pressure projection, moving wall, or
per-fragment lineage/constraint claim. The returned accounting describes ingress,
not subsequent native momentum/energy conservation. Tests include translating
and spinning cubes with analytic inertia and grid-refined projection error,
cut tetrahedra, material density similitude, an actual native step, budget and
cancellation refusals. Native Cargo execution remains separately required.

### Prescribed body work, material heat and latent feedback

`conjugate::creeping_flow_with_body_force_heat` uses the original closed
Stokes solver, including fixed no-slip solid voxels, with one explicit
force-density vector per voxel (zero in solids).
The original MAC half-cell interpolation supplies both the momentum force and
its signed work against the returned face velocities. An independent cell
half-volume partition must agree. Negative body work represents an absorbing
actuator. Only the original molecular and stationary Darcy dissipation enters
heat; wall/body supply, pressure-constraint work and numerical residuals are
separately accounted and never deposited again.

Fixed solid voxels are exact impermeable momentum constraints, not large finite
Darcy coefficients. Normal viscous heat is sampled only in fluid cells. The
reciprocal shear quadrature requires four fluid cells at its MAC edge; a blocked
transverse face instead uses the native stationary half-cell wall coefficient,
including staircase corners. Momentum and heat share that coefficient. The
original parallel half-cell conductances partition wall dissipation into the
adjacent fluid cells; an internal stationary wall supplies zero mechanical work.
Ignored solid viscosity slots neither receive fluid heat nor set the fluid
Reynolds bound. A moving exterior-wall patch without an active tangential fluid
face is refused. Solid heat conduction remains with its existing thermal owner.

The obstacle regressions include a wall-driven cavity above a fixed solid bed
with an independent scalar Stokes reduction, graded heterogeneous half-cell
resistances, and the one-loop circulation around a central solid voxel. They
check zero obstacle flux/heat, exact actuator-work balance, and material-driven
deposition through the existing enthalpy observer. This is fixed geometry; it
does not admit moving solids, all-solid flow, phase activation, or stress release
during melting.

`MaterialFlowFeedback::with_body_force_heating` retains this force field over
an accepted advance. Every actual thermal substep reuses the same viscosity,
phase/permeability, conservative advection and caloric owners. Source-data
changes therefore alter actual velocities, dissipation and subsequent
temperatures or latent fraction. Accepted-duration audits integrate body and
wall work; source, endpoint or cancellation refusal preserves accepted state.

The original closed creeping-domain, fixed-density and fixed-reference-mass
restrictions remain. Force is prescribed external actuator data; this does not
advance a finite electrical/mechanical source, gravitational potential, kinetic
store, moving skeleton or free-surface energy. Any enclosing source must debit
the reported work. Reynolds/power gates and thermal-step refinement still apply.
Focused regressions cover wall/body/Darcy balance, absorbing and pressure-balanced
forces, material substitution, thermal thinning, latent feedback and rollback.
Their presence does not establish native execution or experimental validation.

### Reservoir pressure work feeds native material flow and enthalpy

`conjugate::creeping_flow_with_pressure_work_heat` and
`MaterialFlowFeedback::with_pressure_work_heating` add one explicit opposed
Pressure/Outlet plane pair to the existing quasisteady source-viscosity owner.
These mouths prescribe normal traction `-p_res*n` and ZERO TANGENTIAL VELOCITY
trace. Tangential reaction shear performs zero work. This is a distinct opt-in
from legacy static-pressure/zero-gradient openings; it does not prescribe zero
tangential traction. Other planes remain walls/symmetry. The original closed
and inertial heating modes retain their existing boundary admission.
Nonzero `backflow_tangential_velocity_m_s` refuses in this zero-trace work model;
an explicitly zero vector is equivalent to the existing absent value. The
ordinary pressure-opening momentum option does not silently add unaccounted
reservoir shear work to this separate heat-producing mode.

Native mouth momentum rows mirror one interior cell, but their physical work
volume is HALF that cell. The normal stress coefficient is doubled in that
mirrored row; physical half-row work therefore matches the original cell
`2*mu*(du_n/dn)^2*V`. A mouth shear half-edge uses
`gamma=(u_n_plus-u_n_minus)/l_t + inward*u_t/(h_n/2)` and heat
`mu_edge*(h_n/2)*l_t*w_other*gamma^2`. Both complementary momentum rows use
the derivatives of half this nonnegative quadratic: the normal row has an odd
tangential ghost, and the tangential row has a zero-trace half-cell coefficient
AND reciprocal normal shear. The same two-cell harmonic source viscosity
supplies both derivatives and partitions the heat by resistance fractions.
Original wall, body and grid-aligned porous work likewise use physical mouth
half volumes. Fully rotated porous cells still require a clear-fluid layer
before an opening; their existing admission is not relaxed.

`ViscousPower::pressure_power_w` independently sums
`-(p_res-p_ref)*Q_out` at each mouth using the actual projected face flux that
enters the existing enthalpy advection owner. `supplied_pressure_power_w` is
signed mechanical supply, including absorbing reservoirs. It joins wall and
body supply in the power check, never in the heat source map. The original
local molecular/Darcy/Forchheimer dissipation is deposited exactly once.
`ViscousHeatingAudit::pressure_work_j` and `supplied_pressure_work_j` integrate
that supply over actual accepted substeps, while incoming/outgoing stream
enthalpy remains the original advection account. Source viscosity,
permeability and calibrated nonlinear loss refresh on that same timeline.
`with_body_force_heating` preserves an already selected pressure-work mode.

This remains fixed geometry, fixed density/reference mass and incompressible
caloric storage, with momentum convection omitted and the existing Reynolds
ceiling at most 0.1. It has no transient/open kinetic-energy, compressibility,
buoyancy-work, turbulent or freely turning reservoir-mouth claim. The final
frozen coupled momentum residual, continuity, independent supply/dissipation
balance and source domains all gate the original candidate transaction.
Cancellation or rejection publishes neither motion nor heat. Focused
regressions cover discrete developed channels, source-only viscosity and
nonlinear porous swaps, nonuniform/graded mouth shear, gauge/reversal,
temperature feedback, stream enthalpy and rollback; native execution remains
separate from source/parser validation.

### Shared material buoyancy and electromagnetic transport

`MaterialSlab::step_magnetothermal_with_boussinesq_material_transport` and its
electrode-free counterpart compose the existing current/Joule/Lorentz owner
with the existing Boussinesq material projection in both 2D and 3D. One local
source evaluation supplies resistivity, density-driven buoyancy, viscosity and
thermal conductivity. The original source object remains attached to each
sample. Reference density and heat capacity retain inertia and sensible storage;
explicit departure bounds apply at initial and endpoint states. Material buoyancy
replaces the independent beta load for the step, and Joule heat is deposited once.
The original source accounting, heat balance and atomic refusal gates apply.

Coefficients and sources are explicit initial-state values, with endpoint support
checks rather than within-step convergence. Buoyancy work is not Joule heat;
gravitational/kinetic total-energy conservation, induction, compressible work,
species evolution and phase change remain outside this model. Focused regressions
cover direct force superposition, existing-owner equivalence, data-only material
changes, depth scaling, retained source evaluations and transactional refusals.

### Shared phase-material current, flow and enthalpy evolution

`conjugate::enthalpy::advection::feedback::electrical::advance_joule_enthalpy_with_material_flow`
feeds the existing Ohmic current/account owner from the EXACT applied samples of
`MaterialFlowFeedback`. Each actual thermal substep resolves h-derived phase/T
once per flow sample; those objects drive viscosity, declared buoyancy/mushy
mobility, owned fluid-face conductivity, current and irreversible electrical heat.
The original MAC face mass fluxes advect total enthalpy, including latent heat.
Electrical contact heat keeps its original partition and is deposited once.
Independent inner sources and admitted viscous/Darcy heating remain on the same
conservative timeline. The electrical accessor reads retained material response;
it must not resolve a separate state or mutate external physical history.

The electrical and flow maps must use the same borrowed all-fluid domain and
identical active receiving cells, including ideal electrode voxels. Source and
endpoint h/T/cell identities are checked; final current support uses already
resolved endpoint objects and does not perform another momentum solve. The
original whole-advance current budget includes this final solve. Thermal state,
flow/inner history and electrical work/deposition reports publish together; any
map, material, solver, endpoint, budget or cancellation refusal preserves both
accepted inputs. Joule energy is already inside the thermal external-heat audit.

This retains the selected material-flow momentum mode and fixed-mass/fixed-volume
domain. It adds no magnetic/motional EMF, induction, charge convection, phase-volume
source, free surface or kinetic-plus-thermal conservation claim. Thermal timestep
refinement remains required; the default quasisteady mode additionally requires
rapid momentum equilibration. Retained native momentum has the restrictions below.
Focused synthetic-card regressions compare independent series-current,
contact-heat, Darcy, conduction and donor-cell enthalpy balances; exercise
constant-temperature phase changes, data swaps and zero-current reduction; and
check exact sample identity/query counts plus refused endpoint/retry behavior.

### Shared phase conductivity on the live momentum/enthalpy timeline

`MaterialFlowFeedback::with_phase_conductivity` binds the original registered
equilibrium h chart and shared `LiquidMassFractionConductivity` to owned fluid
faces. The actual local liquid mass fraction scales each source-backed k(T)
before the existing half-cell series resistance. The same conductances feed
pairwise heat exchange and the existing conductive/advective CFL. At a latent
plateau h can therefore change both this conduction and separately attached
mushy Darcy mobility while temperature stays fixed. Base material samples and
receipts remain intact; `phase_conductivity_law()` exposes the applied model.

Admission requires the exact caloric object and constant reference density at
every chart knot. Source and endpoint evaluations refuse unsupported phase
states or unrepresentable scaled conductivity before publishing a coupled
candidate. The ordinary no-law path is unchanged; an identity-law regression
compares accepted heat, flow and history bitwise. A four-cell test independently
closes series-Darcy, half-cell conduction, donor-cell total-h and external-energy
balances; endpoint overflow/cancellation retain the accepted state. Native Rust
execution remains a separate verification requirement.

This is a declared mass-fraction scaling law, not an inferred universal mixture
rule. It does not supply variable-density continuity, distinct phase velocities,
solid locking, moving geometry or a free surface, and it does not add latent heat
or replace the original caloric/source/history owners.

### Source-material enthalpy with retained native momentum

`MaterialFlowFeedback` defaults to the quasisteady mode described above.
`with_unsteady_momentum` instead initializes the existing native MAC momentum
owner at rest at the supplied thermal time. Each accepted substep first advects
total h with the accepted old face field, then performs exactly one native
backward-Euler momentum update over the actual thermal/CFL-limited duration.
Both stages use the starting source epoch: local viscosity, material Boussinesq
forces, declared h-dependent mushy drag and owned thermal/phase coefficients.
Caloric storage, latent heat, conservative advection and current/source histories
retain their original owners. This is a first-order split requiring refinement.

Native face velocities, pressure and the exact accepted clock survive advances.
Clock, starting samples and accepted enthalpy must agree before reuse. The
original thermal clone-and-publish transaction includes this momentum state;
thermal, source, native-solver or cancellation refusal publishes none of the
candidate. `applied().flow` remains the field actually used for heat advection;
`unsteady_flow`, `unsteady_time_s` and `unsteady_step` separately expose endpoint
momentum and its existing native convergence/kinetic-energy diagnostics.

Admission requires closed laminar wall/symmetry boundaries, no fans, fixed
reference mass/density and fixed geometry. Initial nonzero flow is not inferred.
Without the explicit unsteady Stokes heating attachment below, the diagnostic
kinetic energy is not a closed mechanical-plus-thermal work account; pressure,
molecular and Darcy power are not deposited as heat. Steady viscous/body-force
heating refuses in either attachment order.
Existing phase-support, equal-density mushy and no-free-surface limits remain.
`try_new_phase` also supports retained momentum on its original fixed grid.
An exact solid start needs no liquid viscosity or Cp; all-solid intervals
advance the accepted clock with zero native solves. Melting opens the original
faces at their retained zero speed without restarting other moving faces.
Freezing can close a native face only if its actual staggered velocity is
already exactly zero. Otherwise the entire candidate refuses a missing
momentum/kinetic-energy handoff; no speed clipping, pressure projection or
unreported capture energy is introduced. This preserves native kinetic energy
under the admitted zero-speed occupancy changes. It supplies no elastic stress
release, shrinkage, moving skeleton or resolidification reference-state law.
Regressions exercise heated and isothermal latent-plateau transients, retained
history, the native steady limit, actual CFL durations and atomic refusals.
Their addition does not establish an executed native test pass.


### Source-material unsteady Stokes heat and mechanical storage

`MaterialFlowFeedback::with_unsteady_viscous_heating` selects the existing native backward-Euler momentum owner with convection explicitly omitted under `ViscousHeatControls` (Reynolds ceiling at most 0.1). `with_unsteady_body_force_heating` additionally admits prescribed signed actuator work and permits force updates only between accepted advances without restarting momentum. Source-start viscosity and mobile mushy Darcy coefficients remain those already used by the material/enthalpy owner. Current sourced viscosity must admit both retained starting and computed endpoint speeds.

The shared thermal owner chooses the actual combined conductive/advective CFL and event-clipped interval first. Its `step_sources` hook then prepares exactly one native momentum solve and the existing stress/Darcy quadratures supply endpoint physical heat over that duration. Accepted old flow still advects total h; material coefficients are frozen at the source start, and all endpoint material checks see the generated heat. This first-order split requires time and space refinement. `applied().flow` remains the old advection field; `unsteady_mechanical_step` separately exposes the endpoint power and exact native face-mass kinetic account.

The discrete account is `W_wall + W_body + W_pressure = Q_molecular + Q_Darcy + (K_end-K_start) + N_BE`, where `K=sum(M_face*u_face^2)/2` uses exactly the inertia row's staggered half-cell volumes, and `N_BE=sum(M_face*(u_end-u_start)^2)/2`. The nonnegative backward-Euler numerical loss is reported, never turned into heat. Pressure-constraint work must independently remain negligible; it is not heat or a residual correction. Only physical molecular/Darcy watts enter the original caloric kernel, once. `UnsteadyHeatingAudit` accumulates that physical heat, signed supplies, native kinetic change, numerical loss and residual over the accepted advance; generated heat is already inside the original thermal audit. Flow diagnostics use the same native face mass norm, while work/heat balances retain its actual increments.

`UnsteadyMechanicalStep::kinetic_energy_change_j` retains the actual native
face increment `sum(M*(u1-u0)*(u1+u0)/2)` with compensated endpoint differences,
products and face reduction. The independently evaluated backward-Euler loss
uses that same represented difference, retaining a small starting speed even
when native `u1-u0` would round it away. Both quantities retain their reduction
tails through accepted substeps. `UnsteadyHeatingAudit::kinetic_energy_change_j`
also includes actual removed-face storage during stationary-support capture.
Absolute initial/final K remains a descriptive, bit-matched history; subtracting
those rounded totals is not the work/heat gate. The existing implicit caloric
and MHD combined balances consume the retained increment. Physical stress,
porous and capture heat keep their independent original laws and partitions;
no numerical loss or residual is thermalized. Scalar native regressions cover
the endpoint kinetic/BE identity and opposite large transfers. The supported
freezing regression advances the original thermal owner beside an unchanged
large circulation, retaining the small cavity's decay and capture heat even
when its absolute total kinetic diagnostics are identical.

Source-derived Boussinesq forcing also composes with this retained-momentum heat mode. The same applied material epoch supplies density departures, viscosity and conductivity; its buoyancy field adds to any prescribed actuator field in the native momentum equation. The original staggered body-work quadrature evaluates both total force and the buoyancy component on the same endpoint faces. `UnsteadyMechanicalStep::supplied_buoyancy_power_w` and `cell_buoyancy_power_w` retain that signed component; `UnsteadyHeatingAudit::supplied_buoyancy_work_j` accumulates it over accepted intervals. Buoyancy is a SUBSET of total body work, never an additional supply or heat source. Actuator work is total body work minus that subset. Physical molecular/Darcy losses alone enter enthalpy, while numerical loss remains separate.

Admission is closed fixed geometry, laminar source transport, constant reference mass/density and the original fixed-pressure caloric chart. The default source path is Newtonian; the explicit retained-flow rheology binding below owns shear dependence. Material buoyancy retains its explicit small density/Cp-departure bounds and source/endpoint checks. Its signed work represents exchange with the declared Boussinesq forcing model; no resolved gravitational-potential reservoir is introduced. A stationary Darcy frame must share the bulk enthalpy at local thermal equilibrium. Exact phase occupancy uses the same zero-speed face admission above; an all-solid interval has zero mechanical work, physical heat, kinetic change and numerical loss. Locking a moving face, openings, fans, configured extra resistance, moving skeletons and compressible pressure work refuse or remain outside this mode. Prescribed actuator work does not debit an unspecified finite electrical/mechanical reservoir. Duration-dependent histories refuse instantaneous `power_balance_observed` queries.

All state, source snapshots, prepared momentum, heat and accounts remain inside the original thermal deep-clone transaction. A pending interval can be prepared and accepted only once; solver, caloric, endpoint, history, cancellation or budget refusal publishes none of it. Focused G1/G3/G4 regressions cover the independently reduced 2x2 MAC startup/work identity, negative actuator work, drive-off kinetic decay without numerical heating, material feedback, CFL clipping, latent Darcy heating and late rollback. Exact-phase cases additionally cover cold-solid startup, quiescent freeze/remelt, the independent eight-face moving-ring kinetic/work account while its central solid melts, data-only viscosity changes, and atomic refusal to lock a moving face.

Additional G1/G3/G4 cases in `tests/liquid_flow_enthalpy/viscous/unsteady.rs` independently reduce source-density-driven circulation and its work/heat/storage account, compare combined buoyancy and absorbing actuators/walls with the identical prescribed total force, and retain accepted thermal and momentum histories through late refusal/cancellation. These regressions exercise the original native solver and source cards; adding or compiling them is not a claim of executed numerical or experimental validation.

`examples/material_buoyant_flow.rs` exposes this same composition as a runnable
material-swap entry. A bounded synthetic input table becomes one immutable
`MaterialCard`; the existing liquid and caloric owners supply all interpolation
and source evaluation. Local source density drives buoyancy, local viscosity
drives stress, and local conductivity drives the four declared thermal faces.
The accepted native momentum and enthalpy history continue across reports.
CSV output separates physical heat, native kinetic storage, numerical loss,
signed buoyancy work and the combined energy defect. The two supplied material
files change only the density curve. They are numerical examples, not measured
substance data or a spatially converged benchmark. The example retains the
closed single-phase Boussinesq/Stokes restrictions and introduces no additional
solver, state clock or gravitational-potential reservoir. Its actual entry tests
cover material source identity, bounded support and data-only coupled response.

### Local rotating flow in the conjugate finite-volume momentum owner

`SimpleConfig::coriolis_angular_velocity_rad_s` configures constant local
rotation [rad/s] for the original steady and native transient SIMPLE solver.
Velocities and wall motions are relative to the grid axes; centrifugal
potential belongs to caller-supplied effective gravity or pressure. No Euler
acceleration, spherical metric or moving-grid transform is inferred.

The momentum equations include the skew action
`C = P^T V [2*rho*Omega]x P`, where P averages the two native normal faces
of each cell and V is its physical volume. The transpose distributes each
cell's force back through its actual half volumes, including graded meshes.
The density is the SAME constant reference density used by native inertia;
material Boussinesq departures continue to weight gravity only. The three
mass-rate coefficients are cached once, with no rotation-sized field buffer.

The action is separate from refreshed material/external body loads, Darcy
resistance and physical dissipation. Each component solve uses the current
coupled velocity. Prospective convergence and final pressure projection both
check the frozen coupled momentum residual, so stale cross-component forces
cannot be mistaken for convergence. Zero rotation skips the additional source
arithmetic and coupled-residual work.

Supported normal boundaries are closed walls and symmetry, including internal
solid voxels and moving tangential walls. Normal pressure/inlet/outlet rows
need a separate adjoint boundary-work treatment and are refused with nonzero
rotation. The distinct molecular-viscosity split with heterogeneous face
masses retains zero spin; combining nonzero rotation with that inertia metric
is explicitly refused rather than assigning the wrong inertial frequency.

Coriolis contributes no physical work or heat. Native backward Euler/BDF2 keep
their endpoint time convention; backward-Euler kinetic damping remains in the
existing numerical-loss receipt and is never deposited into enthalpy. Strong
rotation may require shorter timesteps, tighter relaxation or more coupled
iterations; exhausted convergence refuses the candidate. The shared solver
carries the same source through its original material convection and thermal
accounts, without a second flow step or altered material-query ownership.

Rotating transient convergence is admitted on the pressure-projected endpoint
using the stricter of the requested inner momentum tolerance and the original
flow tolerance. The refreshed frozen operator supplies the reported momentum
residual; an excessive defect returns to the existing sweep budget. Accepted
faces, pressure and receipts share the same time-level equations. Final march
reporting reuses that accepted endpoint after history advancement, without
auditing an unadvanced next step or applying another pressure correction.

Focused regression definitions cover a manufactured rotating free-slip vortex
and its pressure through the real momentum solve, arbitrary-axis graded-grid
transpose/no-work identities, invalid/open/heterogeneous-inertia admission,
zero configuration, native transient mechanical accounting and refused material
advance followed by clean retry.

### Implicit material flow closes on the thermal endpoint

With the existing `implicit-enthalpy-coupling` feature,
`MaterialFlowFeedback<..., (), ...>::step_implicit` closes one stationary
Stokes/material/enthalpy interval through the ORIGINAL
`BackwardEulerStep::evaluate_with_transport` residual and existing
`NewtonKrylovState`. Unknowns are cell-energy increments in stable CellId order.
Every trial starts from a clone of the same accepted material blueprint and
resolves actual h/T, phase/mushy/porous response, optional quasisteady source
rheology, native flow and physical power. The same endpoint's face mass fluxes
reach `Advected::endpoint_transport`; directed donor enthalpy is not replaced
by a symmetric conductance. Assigned fluid-face conductance uses the same
accepted source samples. The original caloric law retains latent storage.

Native Stokes work is bounded across all Newton, numerical-Jacobian,
line-search and final acceptance evaluations by the original `max_solves`.
Caller-declared joule difference scale and nonlinear controls are explicit;
these are numerical Jacobian actions, not analytic tangents or adjoints.
Source-domain refusal latches and terminates the nonlinear attempt. No failed
trial, source, coefficient or diagnostic is published as physical history.
Local and global native thermal residual gates are independently rechecked.

Only the converged endpoint's molecular/Darcy/Forchheimer heat enters the
native thermal update, once. The existing physical-power owner integrates the
matching wall/reservoir/body supply over that same duration. Net open-stream
enthalpy and independently evaluated authored heat/boundary supply remain
separate report accounts. Internal transfer cancels. Combined actual storage
minus these external supplies must also close. Target-time events stay pending;
events inside the frozen interval refuse through the original thermal owner.

This entry requires the existing Stokes heating mode and an empty wrapped
observer. Retained momentum, buoyancy and mutable observer histories remain
outside its domain. It returns an accepted thermal state plus endpoint report,
not an explicit-start observer history; subsequent calls reuse that thermal
state and the same material blueprint. No explicit advective/conductive CFL
is imposed, but time/space/Jacobian/constitutive refinement remains necessary.
No uniqueness, unconditional nonlinear convergence, inertial/free-surface
flow, moving skeleton or full melt-shape prediction is claimed. Focused native
regressions compare analytic thermal/rate-flow roots, pressure-driven donor
and conductive exchange, latent mobility, constant-source behavior, shared
work budgets, source refusal and deterministic cancellation/retry.

### Implicit material feedback with retained native inertia

Under `implicit-enthalpy-coupling`,
`MaterialFlowFeedback<..., (), ...>::step_unsteady_implicit` closes total
enthalpy, source response and the ORIGINAL unsteady Stokes momentum at one
endpoint. Each thermal Newton trial and constitutive rate iterate restores the
same accepted native face velocity/pressure history. It advances that history
over the exact proposed interval with endpoint viscosity, optional sourced
Darcy/Forchheimer resistance, material buoyancy and explicit actuator force.
The endpoint flow supplies the original directed enthalpy graph; the same
material epoch supplies fluid-face conduction and molecular/porous heat.
Latent fractions can change coefficients while temperature stays constant.
Neither momentum nor caloric storage is reimplemented.

`ImplicitUnsteadyFlowControls` uses the original joule residual/Newton controls
and the shared `RheologyControls` when a native-rate source is attached.
Stationary and transient closure share the actual-rate source admission and
viscosity/conductivity/Cp mismatch calculation. An active source buoyancy
density departure must close at that same rate. `QuasisteadyRheologyControls`
and its report name denote these same numerical types. Only source-rate
coordinates relax; velocities, heat and history never relax or accumulate
across candidates. Every real native solve, including finite-difference and
final acceptance trials, consumes the original `max_solves` ceiling.

The method returns both accepted thermal state and its retained feedback owner.
Resume with those two together, including after changing the prescribed drive.
The original `BackwardEulerStep` independently checks local and total thermal
energy defects. Native work closes physical heat, actual face-mass kinetic
change and nonnegative backward-Euler numerical loss. The latter NEVER heats
material. A combined thermal/kinetic/work gate independently bounds the total
defect, including the separately limited pressure-constraint work. Buoyancy
work remains a subset of body work, not another heat source or a newly modeled
gravitational reservoir. `applied()` describes the implicit endpoint epoch;
`ImplicitUnsteadyFlowReport` owns this interval's accounting, while `audit()`
remains the explicit-advection observer account. Authored sources/contacts stay
unchanged. No trial state or prepared mechanical receipt is published after
source refusal, nonconvergence, exhausted work, cancellation or failed balance.

Admission retains closed walls/symmetry, creeping-flow/source bounds, fixed
reference density, grid and caloric object. Mobile mushy/porous response may
evolve; a change in mobile/solid occupancy refuses until an admitted phase
event/kinetic handoff is supplied. Mutable wrapped observers, moving skeletons,
free surfaces, compressible storage and constitutive memory remain outside this
entry. Target-time events stay pending and interior events refuse. No explicit
thermal CFL is imposed, but time/space/constitutive/Jacobian refinement and
source support remain necessary; nonlinear convergence is not a physical
accuracy certificate. Five focused G1/G3/G4 groups independently reduce
temperature/shear-dependent startup and drive-off, directed donor transport
with Fourier exchange, latent Darcy response, constant-source native parity
and observer continuation, and shared-budget/source/late-cancellation rollback.
Their addition is not a claim that native tests or physical validation ran.

### Quasisteady source rheology closes native pressure flow and heat

`MaterialFlowFeedback::with_quasisteady_rheology` and the original pinned
`LiquidFlowRheology::quasisteady_feedback` bind a declared per-voxel rate seed
to the existing source/flow/enthalpy owner. Zero shear is not presumed to lie
in source support. Before advancing, select the original closed viscous/body
heating or opposed reservoir pressure-work heating mode. This is an additional
explicit constitutive closure; the retained, time-lagged mode below is unchanged.

At fixed source-start enthalpy and temperature, each trial resolves the SAME
source, rebuilds viscosity and viscosity-dependent mushy/porous resistance,
and runs the original native Stokes and physical-power gates. Unit molecular
quadrature observes cell-volume RMS `sqrt(2 dev(D):dev(D))` from the returned
MAC fluxes, including graded reciprocal shear, no-slip half cells and the
pressure-work mouth half-edges. Darcy loss, forcing and numerical residuals
are not kinematic strain. Re-querying the immutable material at this ACTUAL
rate must succeed before any relaxation; unsupported physical states cannot
be hidden by damping the next source coordinate.

`QuasisteadyRheologyControls` bounds constitutive iterations, absolute/relative
rate error, relative viscosity/conductivity/Cp mismatch, and convex relaxation
of the next rate query. Both rate and transport closure are required. Every
actual native trial also consumes the original `max_solves` budget. No trial
calls wrapped source, conductance or thermal-step hooks. Only the final source
epoch, its original physical heat and actual advective flux reach those owners.
Rate-dependent conductivity therefore enters the original half-cell series
conductance and conductive/advective CFL using the SAME accepted source query.

`applied().shear_rate_per_s` retains the coordinates actually queried for those
applied coefficients; `quasisteady_rheology_report()` exposes their bounded
discrepancy from solved rates. `endpoint_shear_rate_per_s()` retains the actual
solved rate. The thermal endpoint is re-admitted at that rate and its actual
new temperature, without claiming a second equilibrium solve at the endpoint.
The next thermal source epoch closes that new equilibrium. Source, nonlinear
or native budgets, arithmetic/power refusal and cancellation preserve the
entire original accepted thermal, flow and wrapped-observer state.

This remains the admitted fixed-geometry/reference-mass incompressible caloric
Stokes approximation. It adds no material memory, yield evolution, open
kinetic-energy store, universal convergence or uniqueness certificate. The
cell-RMS constitutive reduction and explicit thermal split require spatial and
temporal refinement. Focused fixtures use independent affine-source channel
and circulation roots, data-only material swaps, distinct seeds/directions,
rate-dependent heat conduction, thermal feedback, stream enthalpy, constant
source parity, observer counts and rollback. Source/parser validation does not
substitute for native test execution or measured rheology validation.

### Source rheology on the retained momentum and enthalpy timeline

`LiquidFlowRheology` binds the original pinned liquid, density and caloric
sources to actual temperature and physical shear rate. It creates the existing
`MaterialFlowFeedback` with retained native MAC momentum; it adds no fluid or
thermal stepper. The general `with_native_rheology` attachment also composes
with the existing phase-aware source callback and its original admission rules.

The cell-volume RMS `sqrt(2 dev(D):dev(D))` is observed from the retained face
velocities using the original symmetric-stress quadrature, including graded
shear edges and no-slip half-cell walls. A unit molecular coefficient removes
viscosity weighting from this kinematic observation. Darcy drag, LES losses,
cell speed, collision parameters and numerical time loss are not strain rate.
The original scalar Newtonian sampling path remains unchanged.

Each accepted-start enthalpy/temperature/rate resolves one immutable source
epoch. The identical evaluated viscosity drives native full-stress momentum
and optional physical molecular heat. Its conductivity drives the existing
fluid conduction and finite-atmosphere exchange; caloric storage retains its
original owner. `MaterialFlowEpoch::shear_rate_per_s` records the rates that
actually supplied those applied coefficients. Actual endpoint temperature and
native rate are re-admitted before publication, and
`endpoint_shear_rate_per_s` exposes that endpoint separately.

The pure-liquid convenience requires fixed density, original material and
Cp/rho pins, and rate-independent h(T). The general phase attachment queries
an exact solid at zero constitutive rate according to its current caloric
state, without modifying live face velocity or bypassing the original
zero-speed freezing handoff. All source samples, thermal and native momentum
histories, physical heat and wrapped finite-gas state remain inside the same
clone-and-publish transaction. Late source or phase refusal and cancellation
publish none of the candidate.

Coefficients are explicitly lagged from the accepted start: this first-order
generalized-Newtonian split requires time/space refinement and does not add an
implicit nonlinear rheology solve, yield or viscoelastic memory, moving phase
geometry, or a mass/density handoff. Six focused regressions cover an
independent two-interval circulation/work/heat reference, graded-grid strain,
exact constant-Newtonian parity, endpoint-source/cancellation rollback,
rate-dependent conductivity with finite gas exchange, and phase-first solid
admission before the native kinetic handoff. Their addition is not a claim of
executed numerical or experimental validation.




### Open NSF axis with source-bound local linear characteristics

`NsfGrid3::with_characteristic_axis` opens two opposed faces of the original
NSF owner. One non-periodic axis (at least four cells) and a two-cell all-fluid
collar are required. Transverse axes may be periodic or stationary halfway
bounce-back walls; interior stationary obstacles remain supported. A diagonal
population link touching both an open plane and a transverse wall belongs to
the wall, independently of axis ordering. The original pull streaming and
outgoing/incoming mass/momentum quadrature share this ownership, including
one-cell-wide wall channels. The normal FV face keeps its original pressure,
enthalpy/conduction/viscous-work flux; stationary transverse wall heat remains
in the separate wall-heat account. Topology binds before stepping; manufactured
body sources and moving edge walls are excluded. This supersedes the original
closed-grid "no open boundary" limitation only within this narrow domain.

Each face prescribes one incoming outward-normal acoustic characteristic
`J- = delta u_n - delta p/(rho_ref*c_ref)`; the outgoing `J+` comes from the
adjacent live cell. For incoming-wave and pressure-only controls, entropy and
tangential velocity follow the reference normal speed: inward uses the reference,
outward or stationary uses the adjacent cell. Acoustic flow reversal therefore
does not reset those convected modes. The reference density, impedance and caloric tangent
come from the SAME retained `NsfGasMedium`. Reconstructed pressure and temperature
query its actual NASA9 caloric energy and local transport, or the original
GasSpec path; no independent reservoir gamma, density, sound speed, viscosity
or conductivity is accepted. The declared bound on pressure, temperature and
velocity perturbations must lie in (0, 0.1], with the original material support,
Mach, temperature-window and viscosity refusals still active. Composition stays
frozen. The entropy/acoustic boundary reconstruction is explicitly linearized
about the authored reference; NASA9 energy inversion inside the grid remains
nonlinear and source-owned.

This fixed-reference donor is the local linear characteristic approximation:
entropy and shear travel at the reference flow speed while the acoustic modes
travel at that speed plus/minus sound speed. The explicit `PressureTemperature`
source instead retains its declared temperature and reference tangential velocity
on actual inflow, and the adjacent entropy/tangential state on outflow or at rest.
That material-source switch can remain discontinuous at reversal; it is not an
acoustic-only control or a finite exterior material history. No donor smoothing
or nonlinear entropy/enthalpy transport through an attached acoustic duct is
implied. Both paths retain the same original face material and conservative flux.

An explicit cold-only `with_characteristic_thermal_policy([low, high])`
can select `NsfCharacteristicThermalPolicy::EnthalpyFlux` independently on
each original face, including its selected patch. The default `UpwindEntropy`
retains the legacy trajectory and hash. For incoming outward-normal flow, the
new face temperature balances source enthalpy influx with the same half-cell
Fourier flux used by `open_energy_flux`:
`abs(rho*u_n)*(h_face-h_upstream) + 2*kappa*(theta_face-theta_cell) = 0`
in the original lattice chart. Density, conductivity and thermochemical
sensible enthalpy come from the retained material; the calorically perfect
branch uses the original frozen `law.gamma`. Conductivity is evaluated at
the original mean face/cell temperature and pressure. A bounded solve stays
between the cell and upstream temperatures and propagates material refusals.

Outgoing and stationary flow use the adjacent lattice temperature exactly,
with zero boundary conduction. An explicit pressure-temperature drive provides
the upstream temperature; other drives retain the existing isentropic reference
temperature at the reconstructed pressure. Normal acoustic and tangential
characteristics are unchanged. Policy, native face states, trial clones and
rollback remain with the original owner. An explicit zero SGS coefficient
retains molecular transport. This is a thermal contact, not full nonlinear NSCBC, general
flow-reversal admission, finite player heat capacity or species transport.
Focused independent-flux and native-input regressions accompany the selection;
these checks do not establish measured instrument fidelity.

Active Smagorinsky, WALE and Vreman closures also admit `EnthalpyFlux` contacts.
Each acoustic face velocity supplies the same half-cell gradient and adjacent
modeled conductivity as the conservative transport field. The inlet solves
`abs(rho_f*u_n)*(h_f-h_up) + 2*(k_m+k_t)*(theta_f-theta_i) = 0` using the
original source enthalpy and molecular conductivity. Matching deviatoric stress
work remains in the original total-energy flux. A bounded 32-pass face/Guo-force
iteration closes the velocity dependence before each thermal predictor;
nonconvergence refuses instead of using inconsistent heat. Only boundary
closure values are recomputed during the iteration, followed by one complete
SGS field. Outward/stationary conduction remains exactly zero. Attaching SGS
in either order, changing face/patch drives and changing transverse wall
velocity rebind force and conservative history together. Molecular and
exact-zero paths retain their original arithmetic. Six focused regressions
cover the total-conductivity root, axes/orientations/models, attachment order,
patches, worker identity, mass/energy balance, drop/cancel/retry and failed
rebinds. This adds no turbulent inflow spectrum, SGS kinetic-energy equation
or new boundary-transparency qualification; the original finite outer thermal
coupling accuracy remains in force.

The same reconstructed face supplies incoming populations (adjacent
post-collision non-equilibrium extrapolation), the conservative thermal-pressure
force and the original FV enthalpy/conduction/viscous-work flux. Face energy is
not injected as a manufactured source and boundary cells are not overwritten.
The two-ended exchange reports actual streamed mass, streamed momentum plus the
trapezoidal thermal-pressure impulse, FV total energy, the integral of normal
volume flow and the actual pressure-work integral. Positive exchange leaves the
grid. Area-averaged pressure times total flow is not substituted for integrated
pressure work on a nonuniform face. Per-cell pressure/velocity/area traces remain
available for a later explicit spatial projection. Lattice mass flux and
continuum volume flow remain distinct, exposing their discretization difference.

`set_incoming_characteristics` changes the two incoming acoustic data only at an
accepted time boundary. It preserves physical density, velocity and total energy
and rebases the full Guo half-source, including its second moment, together with
the current face-flux history. A refusal publishes nothing. Every accepted step
publishes face exchange with the existing fluid/material state; failed predictor
passes publish neither. The optional boundary has no effect on the closed-grid
population/energy arithmetic. Fixed cell-order reductions retain worker-count
determinism. Face storage is proportional to face cells, not volume cells.

`EnergyLedger::open_boundary_energy_j` is signed total energy entering these
faces; its residual now subtracts this term separately from wall heat. The
viscous-heating bookkeeping subtracts open internal-energy advection and heat,
not the enthalpy/kinetic flux. Its remainder still is not independent evidence
of physical dissipation. Momentum exchange closes the periodic-transverse,
no-obstacle domain;
stationary internal obstacle and transverse-wall reactions are additional wall
impulses. A finite walled passage does not claim periodic-domain momentum closure.

Focused tests in `nsf3_open_boundaries.rs` exercise analytic uniform throughflow
on all axes, separately summed mass/momentum/energy budgets, incoming-wave work,
outgoing-pulse refinement, NASA9 formation-offset invariance and source support,
atomic refusals and workers. Finite-passage cases exercise independent mass and
FV total-energy balance, cyclic coordinate permutations, opposed wall corners,
and separate isothermal-wall heat. `nsf_dg::tests` additionally composes a finite
walled constriction with the original live DG bore. These are numerical-owner
tests, not a measured nozzle or instrument campaign.

The ratio citation in `NsfSpec<G, S>` accepts both literal and owned source text
(`S: Into<String>`). The grid retains its own string through characteristic
previews and acceptance, available from `bulk_ratio_source()`. Data-file input
therefore does not leak buffers or replace a supplied source with a static label.

`nsf_open_channel` is an executable throughflow and
incoming-pulse example with pressure/temperature/speed inputs; dry-air GasSpec
owns its thermodynamics and transport. NASA9 source replacement is exercised
through the same API in the focused tests. These tests/example have not been
executed natively in the authoring environment.

**No-claim boundary.** This is a first-order normal-face, local linear acoustic
characteristic rung, not full nonlinear NSCBC, a qualified obliquely transparent
boundary, moving interfaces, a DG/flow multirate mortar, a resolved flute jet or
instrument-accuracy evidence. An independent NumPy plane-wave calculation of
these formulas (not Rust execution) gave returned/residual pressure L2 ratios
0.08120, 0.04191 and 0.01920 for 64/128/256 cells. That is a refinement check;
it does not meet or establish the plan's -40 dB interface target. Native and
multidimensional qualification remain necessary before stronger claims.

### Moving a completely flooded planar wall beneath a free surface

The native `step_with_planar_wall_motion` and frozen-material
`step_with_planar_wall_motion_si` advance the SAME VOF collision, streaming,
capillary and conversion owner with a prescribed tangential halfway-wall
velocity. Reflected populations and their outgoing partners supply reciprocal
liquid/wall impulse and velocity-dot-impulse mechanical work. Every selected
plane cell must remain fully liquid before and after the step; gas, partial
wetting, a moving contact line, or normal wall motion refuses. The plane has
zero swept normal volume, so this entry does not claim the cover/vacate or
conservative moving-geometry remap required by the broader MR69 bead.

Admission checks absolute and wall-relative physical Mach, positive wet
populations, the explicit relative mass-defect tolerance, and cancellation on
one private native candidate. The material wrapper additionally retains its
original source, dx/dt units and scenario speed bound, completing SI output
admission before installing that candidate. Velocity rebinds are checked on
every call. Zero prescribed motion leaves stationary population bits unchanged.

Raw momentum exchange retains the numerical density-one pressure for its
momentum account. Reports separately identify this reference impulse and the
pressure-perturbation force. An absolute SI wall load requires the existing
explicit pressure reference; the lattice pressure datum is never presented as
an absolute physical pressure. Tangential work is independent of that datum.
The actuator is prescribed, with no finite-wall inertia, thermal deposition,
kinetic-energy closure, moving-contact calibration or curved moving-body claim.

`step_with_planar_wall_motion_and_gas_pressure` and its material
`step_with_planar_wall_motion_and_gas_pressure_si` adapter combine that wall
exchange with the original prescribed gas reconstruction in ONE private native
candidate. The SI field uses the shared physical cell-center/Pascal conversion,
requires an explicit fixed pressure reference, and admits the combined gas plus
capillary density excursion. Applied gas load propagates through the actual
liquid populations into the submerged wall reaction and velocity-dot-impulse
work; no instantaneous pressure-force correction or second fluid step is added.
Zero/reference-pressure loading preserves the ordinary moving-wall bits; zero
wall velocity preserves the pressure-only step. Pressure refusal, late speed or
mass failure, unrepresentable SI output and cancellation leave the whole
accepted liquid state unchanged. The same tangent-only fixed geometry limits
apply: this does not supply piston swept volume, an evolved gas chamber, gas
pressure-volume energy, or an enthalpy moving-wall solve. Uniform loaded-state
traction, raw SI/native parity, pressure propagation to actual wall reaction,
gauge-invariant tangential work and late pressure rollback have four additional
authored regressions; their native execution remains unverified.

Five native regressions specify analytical six-plane shear impulse/work, true
VOF motion, zero-motion bit parity, mechanical power reversal, inertial/relative
speed admission and transactional early/late refusal. Three sourced-material
regressions specify SI/native parity, pressure-reference shifts, viscosity-only
changes in actual flow/reaction, and late scenario-bound rollback. Native test
execution remains a separate gate; authored tests are not a passing claim.

## Source-resolved rheology in moving liquid (MR23/MR27)

The explicit `free_surface3_from_cells_with_rheology` ingress carries shared
temperature/shear-rate liquid sources into the actual native VOF collision.
It reuses the fixed-domain constitutive solver with bracketed safeguards and
force-corrected molecular stress, including the same material query budget,
cancellation and actual-query residual acceptance. The resulting viscosity
supplies collision tau, molecular viscous heating and temperature-dependent liquid, finite-wall and gas-film
thermal transport. Accepted endpoint material/rate observations use the new
native populations; ordinary frozen constructors refuse evolving-rate bundles.
`rheology_at` exposes the physical rate, collision coefficient and local
constitutive residual at an accepted cell.

Thermal substeps hold the initial native stress while resolving temperature and
rate together. This is the existing first-order split, not exact discrete
kinetic/thermal conservation. Numerical conversion carry stores enthalpy without
an invented shear rate, viscosity or heat source. Numerical damping is not
thermalized. Source-domain, constitutive-convergence, cancellation, thermal and
native-flow checks precede publication of both liquid and finite-wall states.

The current domain is constant-density pure liquid, stationary walls and
explicitly omitted capillarity. Shear-dependent caloric storage refuses. No
viscoelastic/yield history, phase mobility, pressure-work or surface-energy
claim is added. Five focused runtime tests plus source-admission coverage check
independent tensor/constitutive roots, body-force correction, data-only
thinning/thickening, native populations/heat, Newtonian bit parity, thermal
conductance and rollback. The existing fixed-domain rheology tests also cover
the shared solver; full native execution remains required.



### Dynamic wetting data reaches the native planar capillary boundary

`material::free_surface::wetting::MaterialFreeSurfaceWetting3` binds both
original advancing/receding angle claims of an ordered substrate/liquid
`InterfaceSystemCard` to one planar domain wall. The liquid phase/process/revision
must match the actual source. Every fixed liquid query coordinate and quantity
descriptor is retained, with only the explicit SI `contact_line_speed` axis
added. Both original claims must support stationary initialization. Subsequent
queries pin those SAME selections, including when another claim could answer
outside their support. Temperature and other dynamic source dependencies are
not silently evolved; the liquid and environment remain the admitted frozen
isothermal context. Source medium, preparation and history remain on the card.
The shared `ResolvedScalarProperty::has_only_owned_and_fixed_external_conditions`
admission exempts exactly the typed speed owned by this runtime. Every other
response and applicability axis uses the existing canonical fixed-external
policy, including typed pressure/composition; their complete values stay fixed
through the trajectory. Thus a pressure-dependent angle law changes native
capillary motion when the prescribed pressure data changes between scenarios.
The existing static-only admission and scalar constant-support rules are intact.

Physical VOF fractions supply the wall history independently of imposed-angle
curvature ghosts. A cached tangential average uses only complete opposite
neighbor pairs, retaining centered affine geometry even at transverse walls.
The first inward physical layer supplies the normal derivative. Inactive
interior gas contributes zero through the existing sparse geometry owner;
explicit solid, outside-domain and piston cells remain excluded. Merely
changing the prescribed angle or activating an empty gas halo cannot move
this measured field. Only boundary tile lanes are scanned after walking the
active tile keys, and the wall-cell budget precedes stencil memoization.
The capillary ghost/smoothing/pressure path remains its original owner.

For outward liquid normal `-grad(phi)/|grad(phi)|`, the
signed geometric speed is `(phi_new - phi_old) * dx / (dt * |grad_t(phi_new)|)`;
the gradient is per LATTICE CELL, not per metre. Positive speed advances liquid.
It is a first-cell spatial estimate and a backward time difference, not a
subcell contact-line reconstruction or the fluid velocity at a no-slip wall.
Newly activated wall history was gas. Nonzero fill motion with a zero tangential
gradient, nonfinite/unrepresentable speed, missing source support or an exceeded
wall-cell budget refuses. Initial history is explicitly stationary.

Above/below the caller's nonnegative speed deadband, the original advancing or
receding endpoint supplies the angle. Inside the closed deadband, the measured
geometric angle is projected into the sourced hysteresis interval. The resulting
PER-CELL field enters the same geometric wall ghost, curvature, Laplace-pressure
reconstruction and actual native momentum/mass step. Endpoint angles drive the
next step, an explicit one-step-lagged closure. A collapsed bracket reproduces
the original homogeneous equilibrium boundary; omitted capillarity leaves the
angle data inert in momentum. Other walls and internal solid tiles are neutral.
Existing static wetting and the legacy Wetting ghost cannot be combined with
this binding. The scalar query-only `query_dynamic_wetting2` remains unchanged.

Native populations, fill, carry, conversion counters, contact history, original
source receipts and the next-step angle field publish together only after the
complete native candidate and all endpoint/source/cancellation gates succeed.
No thermal or surface-energy closure, exact kinematic pinning, calibrated
dynamic-contact-line accuracy, arbitrary substrate geometry, moving wall or
thermocapillary force is claimed. The thermal free-surface capillarity refusal
remains in force. Synthetic G1 tests cover speed/dx/dt scaling, direction and
closed-deadband selection, and all six wall-normal/ghost conventions; G1/G3/G4
integration tests cover collapsed-bracket native equivalence, rate-law data
effects on populations, the sigma-zero limit and transactional source exits.

## Finite atmosphere coupled to material phase and native momentum

`conjugate::enthalpy::advection::feedback::environment::advance_finite_gas_with_material_flow`
composes the original fixed-density material/phase/MAC timeline with a closed,
fixed-volume `IsochoricGasParcel`. Each accepted thermal substep resolves the
condensed phase conductivity and the original gas-film conductivity at the
mean interface temperature and current EOS pressure. A prescribed two-layer
film transfers equal/opposite joules into condensed total enthalpy and gas
NASA-9 internal energy. Gas temperature, pressure and transport refresh through
the same source model. Latent heat, exact solid/mushy occupancy, Darcy mobility,
liquid viscosity, retained face momentum and optional molecular heating remain
with the original owners; no second phase or flow solver is added.

The source-derived gas `cv` bound and original condensed caloric slope combine
with conductive/advective timestep limits. Exposed-area coverage, original gas
transport claims, explicit film geometry and duplicate ambient-convection
refusal prevent hidden parallel heat paths. Gas and combined thermal ledgers
check actual storage against the transferred energy; that internal exchange is
removed from the union's external heat. Mechanical heating retains its original
work/kinetic account and actual accepted `step_sources` interval. The returned
gas, thermal and flow candidates publish together; late source, phase, momentum
or cancellation failures preserve all inputs.

This is a frozen-composition, stationary-interface, well-mixed gas/film model.
Gas pressure is observed but does not impose traction or volume work. Reaction,
evaporation, species transfer, moving solid geometry and nonzero-momentum freeze
handoff remain outside this path. Four runtime regressions independently check
NASA internal-energy cooling and latent melting, gas-data-only changes that open
native flow faces, late pressure-domain/cancellation rollback and small-reservoir
capacity limits. Coefficients are synthetic. MR25/MR29 remain open for their
broader moving solid/liquid integration.

### Source convection on stationary material/phase contacts

`environment::advance_finite_gas_with_material_convection` reuses the same
finite-gas/material-flow observer and original thermal timeline.
`FiniteGasConvectionContact` adds an explicitly authored physical unit normal
to the existing exposed cell, area and half-cell geometry. Active area and
material distance are positive; gas thickness is zero because convection
replaces the prescribed film. Unit normals admit only normalization roundoff;
no voxel face or flow projection is invented.

Each source epoch retains the original current scalar conductivity and existing
phase correction across the condensed half-cell. A constant-k representation
of that explicit finite-volume approximation feeds the shared material surface
root. Every actual surface trial queries the same gas at mean surface/bulk
temperature and current EOS pressure, retains its selected source claims and
calls `FiniteGasConvection` for original h evidence. Driven h is positive;
`None` denotes exact equilibrium. Declared zero conductivity remains insulating.
This held-cell-phase path does not resolve an internal phase front or a
source-integrated k(T) profile.

Actual signed gas heat enters condensed enthalpy and the opposite finite-gas
energy once. Native latent fraction, exact solid/mushy face occupancy, retained
momentum, advection, molecular heating, combined capacity bounds and rollback
remain with their original owners. Reported secant conductance is a frozen
heat-rate representation, not a nonlinear derivative. Surface temperature,
material power, continuity residual, iteration count and coefficient evidence
remain in `FiniteGasContactSample::convection`. Provider failures retain their
original cell/material error; exhausted root work yields the typed
`FiniteGasFlowError::SurfaceBalance`.

Four focused native specifications cover independent half-cell/film and NASA
internal-energy equations in both heat directions, latent heat, actual phase
face/momentum response, equilibrium and transactional refusals. Syntax and
whitespace checks passed; native compilation/execution remain pending an
admitted constellation build.

## Material MHD drives native momentum and enthalpy reciprocally

`conjugate::enthalpy::advection::feedback::mhd::advance_mhd_enthalpy_with_material_flow`
joins the existing moving-conductor current operator to retained native MAC
unsteady Stokes momentum and the original conservative enthalpy timeline.
Resistivity, viscosity, conductivity, buoyancy and admitted mobile mushy drag
read the same material objects at each actual starting h/T. A bounded implicit
iteration solves current with `u × B`, applies the native reciprocal Lorentz
force, and advances momentum from the same accepted history over the actual
thermal interval. Its endpoint cell velocity uses the same half-cell force
interpolation as the electrical work identity. Independent velocity and work
gates refuse incomplete coupling; every real current/momentum trial is budgeted.
An inner duration-dependent observer publishes one interval, not one per trial.

Only conductor-frame `I²R` and native molecular/Darcy dissipation enter thermal
storage. Signed Lorentz work transfers energy between current and momentum;
magnetic braking can heat material while electrode power is zero. Native
kinetic storage, backward-Euler numerical loss and bounded pressure-constraint
work retain separate accounts. Electrical sampled work, force work on actual
endpoint velocity, the native face-force lift and mapped heat are independently
checked. Absolute accumulated transfer defects cannot cancel between intervals.
The combined thermal/kinetic/work gate includes original thermal input, external
body and buoyancy work, wall work and electrode supply without double-counting
physical loss. The final current support solve uses actual endpoint h/T and
velocity and supplies no extra heat.

This entry requires a complete shared fixed-volume mobile map and the original
unsteady physical-heating admission. The imposed magnetic source is stationary;
the caller owns solenoidal-field validity and low-induction applicability. There
is no induced-field energy store, Hall response, moving boundary or force
handoff into a conducting frozen solid. Exact solid locking is explicitly
refused; existing mobile mushy drag remains available. Electrode-connected and
fully isolated current topologies are explicit choices. Fixed original contact
resistance and heat partition are preserved. Changing drive or B between
accepted advances retains native momentum. Material, budget, convergence,
endpoint and cancellation failures preserve both original states.

Synthetic regressions reduce a driven 3x2 MAC box independently to its current,
Lorentz-force, viscous-loss and face-mass equations. They cover resistivity-data
substitution, B reversal, external-force composition, zero-voltage magnetic
braking, actual temperature feedback, one material epoch across repeated
current trials, and late support-solve/nonconvergence rollback. This closes the
bounded MR43 moving-current/fluid/thermal connection; general phase transition,
induction and moving solid/fluid mechanical coupling remain broader work.

### Source-owned electrical interfaces on live MHD

`conjugate::enthalpy::advection::feedback::mhd::interfaces` adds
`advance_mhd_enthalpy_with_material_interfaces` and
`advance_mhd_enthalpy_with_anchored_interfaces`. Both reuse the original MHD
current/momentum/enthalpy driver. Ordered interface callbacks receive the exact
retained adjacent bulk/phase samples used by the same current trial. The shared
interface response owns finite resistance and explicit heat partition, or an
explicit supported opening; no interface temperature, pressure, surface law or
missing-data fallback is inferred. Original source payloads need not be Clone.
The existing source binder checks actual material/process/texture surface roles.

Resolved contacts and openings enter the native moving-current topology solve
before electrical reachability and assembly. That ONE solution supplies bulk
and contact heat, reciprocal Lorentz force and its work. An opening removes the
face from current, motional voltage, contact heat and Lorentz transfer, while
independent material, thermal and fluid owners retain their state. Zero
resistance remains perfect contact. Authored coefficients are not added to a
resolved source resistance; unlisted faces keep perfect-contact semantics.
Wired components still require their native voltage references; the isolated
owner retains its per-cell current-continuity admission.

Accepted interface epochs accompany the actual accepted MHD intervals; rejected
Picard trials cannot publish contact or observer history. Fresh endpoint source
support uses returned h/T and velocity without depositing another heat increment.
Source, topology, current/momentum budget, convergence and cancellation failures
preserve the original state. Existing combined energy gates, fixed-field
low-induction and fixed-geometry assumptions remain unchanged. Interfaces own no
separate momentum or energy storage. Source changes are sampled at actual thermal
interval starts; discontinuities need explicit advance endpoints and refinement.

Synthetic tests in `liquid_flow_enthalpy/mhd/interfaces.rs` compare native motion,
current, ordered contact heat and total energy against an independent MAC/circuit
reduction. They exercise original interface-card substitution and temperature
feedback, actual surface order, non-Clone receipts, opening with retained viscous
coasting, reclosing, and late source-support refusal. These are numerical
regressions, not measured material or contact validation.

## Isothermal source rheology with native capillarity

`material::free_surface::rheology::IsothermalRheologyLaw::try_new` binds one immutable
material card and complete resolved bundle to fixed temperature, pressure,
composition and source environment. Only the explicitly named physical shear
rate `sqrt(2 dev(D):dev(D))` changes. The dynamic/kinematic basis and sourced
versus neglected capillarity are explicit. No caloric data is requested. The
shared material admission permits temperature/rate source response and refuses
unowned dynamic axes; density and surface tension must be shear-independent.
All runtime queries retain the original selected claims and full context.

`MaterialFreeSurfaceRheology3` uses the existing physical-unit and cell-state
ingress, molecular force-corrected stress and shared bounded `solve_local`
owner. Its sampled-bracket safeguards preserve source admission, evaluation
counts and acceptance from the actual residual. Each wet collision lane,
including zero-mass interface ghosts, closes rate against `tau = 1/2 + 3 nu dt/dx^2`. Those exact material coefficients drive the
existing native collision, streaming, VOF mass exchange and conversion. The
native capillary reconstruction retains the sourced uniform sigma; frozen and
evolving-viscosity paths share the same dimensional and positive-boundary-density
admission. No clipping, material-name dispatch or duplicate viscosity fit occurs.

The accepted per-cell source states drive one explicit step on a private native
candidate. Its new stress and speed fields must pass source, convergence and
scenario bounds before populations, fill, carry, conversion counters and source
receipts are published together. Cancellation is polled between constitutive
queries and before publication. The unregioned mass carry has no invented shear
rate or material collision cell. Original native mass accounting is unchanged.

`from_cells_with_pressure_reference` fixes the positive absolute-pressure datum
through the existing material scaling owner before the trajectory is returned.
`step_with_gas_pressure` then accepts absolute Pascals at exposed physical cell
centers, using the shared conversion and original gas-plus-capillary
reconstruction. Pressure and current source-resolved local relaxation act in
one native step without allocating unused thermal transport routes. The actual
new stress closes the next constitutive samples before joint publication;
pressure, endpoint source, convergence, speed or cancellation refusal rolls
back flow and all material/rate receipts together. Zero/reference loading keeps
the ordinary path's population bits; with the default fixed-context law, a
common pressure/reference shift keeps the same dynamics. Numerical density
retains its original incompressible meaning: this default law does not rebind
the material's fixed source-pressure context,
evolve a gas reservoir, or add pressure-volume heat/work. Three focused
regressions specify Newtonian pressure-loaded parity, gauge invariance, actual
pressure-driven thinning/thickening with independent stress/constitutive roots,
mass retention, reference-load identity and late pressure/source rollback.
Their native execution remains a separate, unverified gate.

`IsothermalRheologyLaw::try_new_with_pressure` explicitly owns viscosity
dependence on absolute mechanical pressure as well as shear. It uses the shared
material admission's pressure marker and the SAME selected source claims. The
trajectory's fixed pressure datum must equal the initial source pressure; every
wet-cell query then uses `p_ref + cs²*(rho_lattice-1)*rho_source*(dx/dt)²` from
the actual liquid populations. External gas loading is not substituted for
interior pressure. That pressure is held within each local shear root; the
resulting viscosity drives the next native collision. Pressure is reconstructed
once per wet endpoint, with no second momentum solve or thermal route storage.
This is explicit coefficient feedback, not an implicit pressure/shear time solve.

Density and active surface tension must have no pressure/shear response axes
and must remain exactly invariant at every source sample. This excludes physical
liquid compression and variable-sigma traction. Positive absolute pressure,
original source validity, local rheology convergence and speed must all pass
before flow and material receipts publish together. Temperature, composition
and other context remain fixed. Shifting actual pressure and its datum can now
change viscosity and motion, as required by a pressure-dependent source.
Frozen consumers refuse this explicit pressure owner. The thermal composition
below owns pressure-dependent transport with pressure-independent caloric storage.
G0/G1/G3/G4 regressions specify the actual interior pressure coordinate, an
independent pressure/shear constitutive root, constant-law native bit parity,
absolute-pressure-dependent motion, source-domain rollback, cancellation and
unsupported density/sigma or missing thermal pressure-reference refusal. No new runtime execution,
continuum convergence, physical calibration or total-energy claim is made.
For the incompressible pressure/shear-dependent constitutive model family see
Bulíček, Majdoub and Málek (2010), <https://arxiv.org/abs/1004.0835>; its analytic
existence results are not asserted for this discrete LBM/VOF implementation.

This is isothermal generalized-Newtonian response, with the explicit numerical
contact ghost inherited from the domain. Stationary walls remain the default;
the distributed tangential-wall entry below owns prescribed wall motion.
Thermal storage, viscous heat, viscoelastic or yield history, variable-sigma
traction, phase change and conservative surface energy remain outside this owner.
Source and numerical admission are not stability or experimental-accuracy
certificates. Focused tests in `tests/material_free_surface_rheology.rs` exercise
an independent force-corrected tensor and analytical constitutive root, native
Newtonian bit parity, source-only viscosity/capillary motion changes with mass
conservation and deterministic replay, and endpoint source/convergence rollback.

### Source viscosity with distributed wall motion

`MaterialFreeSurfaceRheology3::step_with_planar_wall_velocity_field_si` admits
the original nonuniform tangential wall field in SI units. Each wet cell's
accepted source viscosity supplies relaxation in the same native collision and
reflection that measure wall impulse, torque, face couples and actuator work.
This lets the existing isothermal material assignment drive a shear wall without
inventing heat capacity or conductivity. The original fixed source context and
selected claims remain attached to the trajectory.

Actual new liquid pressure and force-corrected strain close the next constitutive
queries before publication. Pressure-dependent laws retain their fixed initial
datum and interior pressure reconstruction. Wall pressure is not used as a
substitute for the local source coordinate. The constitutive update is explicit
in time; endpoint admission does not imply an implicit wall/material solve.
The unchanged `MaterialPlanarWallFieldReport3` carries the measured SI exchange.
Mechanical work is not converted into an undeclared heat source or a claim of
total-energy closure; the thermal rheology owner remains the finite-heat path.

Material and native field budgets precede one private candidate clone. The
shared native candidate accepts the existing local relaxation map and allocates
no thermal mass-route trace. Native mass/positivity/Mach, source validity,
convergence, speed, SI conversion and cancellation all precede joint publication
of populations, fill, carry, conversion counters and material/rate samples.
The geometry remains a fixed, completely flooded plane with held tangential
face velocities. Normal wall motion, contact-line evolution and simultaneous
gas-pressure actuation are not added by this entry.

`tests/material_free_surface_rheology/motion.rs` specifies exact constant-source
population/work parity with the existing wall owner, independent pressure/shear
roots and face-impulse work, shear-thinning/thickening changes to actual flow and
actuator work, and late source refusal with cancellation and safe retry. These
are synthetic constitutive regressions, not measured material validation.

## Pressure-dependent viscosity on the thermal VOF timeline

`free_surface3_from_cells_with_rheology` also accepts the existing explicitly
pressure-marked source bundle. Its `with_pressure_reference` value must equal
the original positive absolute source pressure. The same scaling owner supplies
`p_ref + cs²*(rho_lattice-1)*rho0*(dx/dt)²` from each actual wet cell; prescribed
external gas pressure never substitutes for this interior state. That pressure
is held with the original force-corrected stress while the existing bounded
constitutive root and thermal substeps query the SAME selected material claims
at their current temperature and physical shear rate. The actual new native
pressure, temperature and shear are re-admitted before joint publication.

Viscosity and conductivity may respond to these owned transport coordinates.
Pressure-response axes on density, specific heat or surface tension refuse,
even when sampled values happen to coincide. The original fixed-reference
caloric h(T), mass transport, optional buoyancy/capillary owners and thermal
balance remain unchanged. Pressure-dependent stored energy, physical liquid
compression and pressure-volume heating are not inferred. This is the
constant-density, pressure-independent caloric model with pressure/temperature/
shear-dependent transport described by Bulíček, Málek and Rajagopal (2009),
§2.1–2.2, DOI:10.1137/07069540X; its existence theorem is not a discrete claim.

Molecular heating retains its original first-order source split: accepted-start
mu and strain supply held shear-plus-tied-bulk heat; heated material then refreshes
the collision viscosity before the single native step. Both use this same law
and actual held pressure, at their distinct thermal epochs. No new viscosity
fit, solver, numerical-loss heating or source-format extension is introduced.
Original two-axis surfaces can express T/p or p/shear response while the query
retains all three coordinates. Unmarked T/shear paths keep their original bits.

The focused `material_free_surface_enthalpy/rheology/pressure.rs` cases compare
an independent pressure/shear quadratic root and physical heat, exercise actual
T/p feedback and constant-law bit parity, and refuse unsupported caloric axes,
late source-domain failures and cancellation without changing accepted state.

### Anchored conducting phase changes retain reciprocal MHD

`conjugate::enthalpy::advection::feedback::mhd::advance_mhd_enthalpy_with_anchored_phase`
explicitly selects an external stationary anchor for exact frozen cells. The
existing mobile MHD entry still refuses implicit exact-solid support. Both
entries share the same native current, MAC momentum and conservative enthalpy
algorithm. Original cell IDs, source snapshots, conductor contacts and thermal
mass remain unchanged through solid, mushy and liquid epochs. Frozen conductors
carry source-resolved current and receive I²R without querying liquid-only
viscosity or Cp.

The native current owner's Lorentz density remains available for every voxel.
In a frozen voxel its negative times original volume is reported as
`frozen_lorentz_reaction_n`, the anchor-on-material force in N; the anchor receives
the opposite force. Only this supported magnetic load is excluded from the
mobile momentum RHS. Accepted intervals integrate that force over actual
durations into `frozen_lorentz_reaction_impulse_n_s`; endpoint melting does not
erase the preceding support impulse. This covers direct Lorentz support, not
native fluid/solid or outer-wall traction. Other frozen body/buoyancy loads
refuse. Supported material must have exactly zero velocity, so
`frozen_support_work_j` is exactly zero. Original reciprocal force/work and
combined energy gates remain active on actual accepted mobile velocity. No
support work or numerical loss is deposited as heat.

All-solid intervals solve current and advance heat on the original timeline
without consuming a momentum solve or fabricating liquid transport. The existing
phase owner releases melted native faces at zero speed and retains all other
momentum. Freezing moving faces still refuses: this supplies no capture impulse,
kinetic-to-thermal handoff, deforming/moving support, shrinkage, phase-density jump
or free-surface transition. Static imposed-field, low-induction, closed unsteady
Stokes, fixed-contact, equal phase-density and time-refinement assumptions apply.

Four focused regressions in `liquid_flow_enthalpy/mushy/continuity/anchored_mhd.rs`
cover independently predicted frozen Ohmic heat and signed support under material
resistivity/field changes; melting into viscosity-dependent native MHD flow;
moving-face freezing and endpoint current-budget rollback; and a conducting
frozen center with a driven liquid ring. These are synthetic numerical fixtures.
Reconstructed sources received source/API/physics review; native parsing,
typechecking and execution were unavailable in this session.

### Spatial enthalpy through tetrahedral liquid handoff (MR28 increment)

`material::free_surface::transfer::enthalpy` carries nodal P1 specific enthalpy
through the original exact tetrahedron/cell overlap moments used for momentum.
The existing mass, linear/angular momentum and mechanical projection gates remain
active. Kinetic and gravitational projection discrepancies remain numerical
accounts; neither becomes heat. An independent source integral
`sum(rho*V*sum(h_i)/4)` is checked against the actual native `sum(m*h)` with a
caller-selected absolute energy bound.

All referenced nodes and intervening caloric knots must remain fully liquid at
the original density. The original pinned source-line resolver also covers the
complete temperature interval for density, viscosity, heat capacity and
conductivity, including each property's own breakpoints. Source holes, unsupported
line topology, claim switching and interior density changes refuse. Source
queries and geometric work retain explicit bounds. The same immutable material
card and original thermal attachment own caloric inversion, source receipts and
local viscosity; no default temperature is invented for a missing thermal cell.

The result is the existing `MaterialFreeSurfaceEnthalpy3`, ready for its native
mass transport, heating, conduction and viscosity-dependent motion. Source,
projection, attachment, budget and cancellation failures return no accepted
candidate. Four focused regressions cover independently integrated spinning/cut
geometry, continued native flow, material-only viscosity substitutions and
hidden caloric/source-density changes with late refusal and deterministic retry.
These synthetic fixtures received source/API and independent arithmetic review;
Rust parsing, compilation and test execution were unavailable in this session.

MR28 remains open. This handoff requires already-liquid material without stored elastic/plastic history,
one constant-density caloric chart and explicitly omitted capillarity. It does
not discard solid elastic/plastic history, implement melting release, change
phase volume, supply pressure/capillary work or certify total fluid mechanical
energy conservation.


### Source rheology through conservative tetrahedral liquid transfer

`free_surface3_from_tetrahedra_with_enthalpy_and_rheology` carries the original temperature/shear material law into the existing thermal VOF owner. Each source tetrahedron supplies its actual affine velocity gradient and `sqrt(2 dev(D):dev(D))`. Original pinned property-line admission covers that tetrahedron's complete P1 temperature span at its actual rate; a reference rate cannot conceal unsupported source shear. Shared vertices are admitted separately when adjacent tetrahedra have different strain. Source-query and original geometric budgets remain explicit.

The original exact overlap, per-fragment mass, linear/angular momentum and enthalpy transfer remain authoritative. The receiving native owner resolves its own force-corrected shear and source viscosity, retaining the same law through continuing motion and molecular shear-plus-tied-bulk heat. Source rates are reported separately from receiving constitutive samples. The receiving viscous stress is initialized from the overlap-projected symmetric P1 strain and the original source at the receiving enthalpy and rate. Nonlinear stress is evaluated at that projected state; an average of source stress tensors is not implied, and elastic memory remains outside this handoff. Kinetic/potential projection discrepancies remain numerical accounts, never heat. The existing `with_time_origin` preserves an explicitly supplied source clock before continuation.

Two synthetic regressions specify data-only rheology changes in actual native motion/dissipation, independent mass/enthalpy and both-momentum checks, source-versus-receiver rates, rigid-rotation objectivity, out-of-domain source shear, bounded admission, late receiver refusal and cancellation. Constant reference density, fully liquid source, absent stored elastic/plastic history and neglected capillary surface energy remain required. MR28, MR29 and MR69 remain open for their broader phase, stress-history and moving-geometry obligations.

### Source-density buoyancy across conservative tetrahedral liquid transfer

`free_surface3_from_tetrahedra_with_boussinesq_enthalpy` carries original
source-density feedback into the existing thermal VOF trajectory. Tetrahedral
mass, geometric volume, inertia and the caloric chart retain the original
reference density rho0. Every nodal/caloric-knot/property-line sample uses the
shared `BoussinesqMaterial` admission for source density and heat-capacity
departures, including interior source breakpoints. Actual source rho(T) and
its receipts remain available; no material-name expansion coefficient is
introduced. Optional original rheology controls compose through the same
thermal attachment and constitutive samples.

The native receiving owner applies its material-owned gravity to current carrier
density, with the bounded source-density departure in the acceleration.
Its Guo-frame rebind preserves physical
velocity and stress without adding an impulse, time or heat. Transfer rereads
actual receiving mass and both momenta after attachment, using the original
per-fragment projection and shared absolute error budgets. Mechanical
projection differences remain separate from thermodynamic storage. Later
heating, advection and molecular dissipation refresh that same material/force
trajectory; local dynamic viscosity is divided by rho0 for collision.

Three synthetic regressions specify independently moving fragments, analytic
gravity impulses before/after heating, zero-departure parity, simultaneous
source shear and buoyancy, reference-inertia invariance, interior density
bounds and late gravity/cancellation refusal. This bounded Boussinesq
feedback does not add volume expansion, phase-density change, pressure work or a new
gravitational-energy store. Fully liquid material, absent stored elastic/plastic
history and neglected capillary surface energy remain required. The broader
MR28/MR29/MR69 deliveries remain open.

### Constant sourced capillarity in the original heated free surface

`MaterialFreeSurfaceEnthalpy3` now retains sourced capillarity when sigma
is declared `ConstantWithinValidity`, has no response axes and is bound
into the same complete liquid bundle as density, heat capacity, conductivity
and the selected viscosity basis. Original card, query meanings,
requirements, selected claims and coefficient identity are preserved through
every thermal/rheology stage and native endpoint. A varying curve with equal
sampled values does not establish a constant law; a missing active sigma in
the old four-property thermal bundle refuses.

The original native capillary scaling and pressure reconstruction continue
unchanged. Temperature/shear-dependent transport, finite gas storage and
finite wall heat use the same retained source trajectory. Leaving the
selected sigma claim's validity domain refuses the complete candidate even
if a fresh query could select a different fallback claim. Liquid populations,
enthalpy, finite wall stores and finite gas state publish together.

Four focused regressions in
`tests/material_free_surface_enthalpy/conduction/capillarity.rs` specify
bitwise native capillary parity at uniform enthalpy, sigma-data substitution
in actual populations, independent four-substep NASA gas/nonlinear-wall heat
references, real temperature/shear feedback into continuing momentum,
complete-bundle admission and late source-domain rollback/retry/cancellation.
The independent reference places the late validity crossing after one
admitted internal thermal substep. These are synthetic numerical fixtures,
not material calibration or experimental validation. Source/API review and
independent scalar arithmetic completed; the reconstructed Rust sources and
tests could not be parsed, compiled or executed in the unavailable executor.

The enthalpy ledger does not invent capillary heat or claim conservative
surface energy or total mechanical energy. This ordinary thermal binding
continues to require constant sigma; the explicit planar thermocapillary
entry point below admits a bounded temperature-varying law. Evaporation and
phase-density change remain unsupported. Existing single-liquid,
constant-density, stationary-wall, native pressure-boundary and low-Mach
assumptions remain in force.


### Thermal free-surface material wetting

`MaterialFreeSurfaceEnthalpy3::with_dynamic_planar_wetting` attaches the
original ordered interface-card closure to an already admitted pure-liquid
trajectory. Frozen and thermal flow share one `WettingBinding`, geometric
contact-speed estimator, hysteresis selection and native contact-angle-field
installation. Original frozen isothermal arithmetic is retained. Sigma remains
the separately admitted source-declared constant, not a hidden function of
contact angle.

Each wet/interface sample uses the original liquid material tuple's actual
caloric temperature. Zero-mass wet/interface ghosts retain the enthalpy owner's
existing temperature extension. Genuinely dry wall stencil cells retain the
explicit initial source context and create no thermal mass. Other substrate
and environmental coordinates remain fixed; no resolved interface-temperature
model is inferred.

Before collision, thermal feedback re-evaluates both pinned advancing/receding
angle claims while holding the last accepted signed contact-line speed. Those
angles drive the SAME interval's native capillary curvature/pressure update.
Endpoint admission uses advected temperature and the original backward
geometric speed, retaining one accepted history. Conductive substeps also check
interface source support, so leaving and reentering a source interval cannot
escape admission. Populations, mass, h, angle field, source receipts, contact
history, finite wall heat and finite gas state publish only after the complete
native candidate succeeds.

This is a first-order local-cell wetting closure with the original spatial and
contact-line limits. It supplies no moving substrate, calibrated pinning,
Marangoni force, capillary heat, or conservative surface/total mechanical
energy claim. Changing contact angle is not credited as thermal energy.

Four synthetic regression definitions cover original frozen/native parity;
same-interval temperature -> sourced angle -> population sensitivity and
neglected-capillary inertness; finite wall/gas/rheology feedback with independent
wall heat balance; and late thermal/endpoint source failure, rollback, retry
and cancellation. The independently calculated first two wall-heating endpoints
are 400.961538461538 K and 401.864956392677 K, placing the 401.2 K source failure
after one valid internal substep. Final changed files passed Rust syntax parsing,
source/API review and bounded independent ownership review. Native compilation
and execution were unavailable (no Rust/Cargo, DSR or RCH); these are authored
regressions, not passing native or experimental evidence.
Focused command:
`cargo test -p fs-lbm --test material_free_surface_enthalpy conduction::capillarity::wetting`.

`with_equilibrium_planar_wetting` adds an explicitly equilibrium alternative
on that same thermal/native transaction. It selects either an original
equilibrium-angle claim or the existing Young--Dupre law from reversible
`equilibrium-work-of-adhesion` W(T) and the SAME liquid/gas tension retained
by the thermal material bundle. Ordinary attachment requires constant sigma;
the explicit uniform thermal-capillarity entry below or the bounded sidewall
height graph can evolve it. The actual local
liquid caloric temperature refreshes the selected interface claim before
collision, at conductive substeps, and after advection. Source substitution,
loss of partial wetting (`0 < W < 2 sigma`), and the original tension source's
domain exit refuse the complete step. Dry numerical ghosts retain reference
conditions. Unchanged reference-temperature sources are reused verbatim,
without per-cell claim resolution.

Equilibrium samples retain their explicit angle origin and
`ContactLineRegime2::Equilibrium`; they contain the original work or angle
receipt, without fabricated hysteresis claims or a contact-speed query.
This boundary does not need a resolved contact-line velocity, and therefore
does not reject locally flat filling stencils merely for a zero tangential
gradient. An implicit Cox extension refuses rather than turning equilibrium work
into a measured advancing/receding bracket. The explicit rate-law selector below
retains the equilibrium source instead. The separate planar and local
thermocapillary geometries remain incompatible with imposed wall angles;
the height-graph sidewall domain is specified below. The local-cell
interface-temperature approximation, stationary homogeneous-wall assumption
and absence of closed surface-energy/heat accounting remain explicit.

Four focused regressions in
`tests/material_free_surface_enthalpy/conduction/wetting/equilibrium.rs`
cover constant-source native parity, local W(T) feedback into the same
momentum step with conserved mass and explicit heat, original-source and
partial-wetting rollback/retry, and refusal of an implicit Cox extension.
They are synthetic numerical tests, not experimental wetting validation.

`MaterialFreeSurfaceEnthalpy3::with_dynamic_equilibrium_wetting` explicitly
extends that original equilibrium binding with either
`EquilibriumWettingRateLaw::CoxVoinov` or `MolecularKinetic`. The existing
full-angle Cox or Blake--Haynes law uses the sourced equilibrium angle, or
Young--Dupre angle from the original W(T) and SAME native sigma(T), as its
microscopic endpoint. Cutoffs and molecular calibration remain separate model
inputs. Accepted VOF history supplies signed contact-line speed; this smooth
homogeneous-wall mode introduces no hysteresis interval or pinning deadband.
Samples retain one original equilibrium claim with no fabricated speed axis,
plus explicit dynamic-model diagnostics. Both thermal preparation and the
advected endpoint refresh the selected source and native boundary together.

Only constant or spatially uniform thermal tension is admitted for this moving
boundary; the equilibrium height-graph extension remains static. Cox requires
a Newtonian liquid and its declared low-Ca domain. Unresolved moving geometry,
source/calibration exits, duplicate rate laws and cancellation refuse the whole
candidate. Molecular loss remains a per-length diagnostic, not caloric heat.
First-cell spatial and backward-time speed approximations, stationary-wall
assumptions and the absence of closed surface/mechanical energy accounting
remain. Focused synthetic regression definitions check independent forward Cox
and molecular equations, original thermal adhesion receipts, native population
response, and rollback; they do not establish experimental calibration or
unconditional stability.

`MaterialFreeSurfaceEnthalpy3::with_dynamic_equilibrium_gas_wetting` attaches
the original source-resolved two-phase Cox owner to that equilibrium binding.
The microscopic angle follows the original angle or Young--Dupre work claim,
or the independently sourced pinning thresholds when explicitly selected below.
The wedge ratio uses actual outer-gas viscosity divided by local liquid viscosity. Gas and liquid keep distinct temperatures, source identities and
state coordinates. No contact-speed axis is added to the equilibrium receipt.

Ordinary steps retain the explicit gas snapshot. Existing finite-atmosphere
heat-transfer entry points instead refresh the same original gas parcel at
thermal preparation and accepted endpoints, including its admitted temperature,
pressure and transport coefficients. Original source matching and both phases'
low-capillary-number limits remain mandatory. This selector shares rate-law
exclusivity and the constant/uniform-tension, Newtonian and planar-wall domain;
it adds no separate gas evolution or contact-line heat model. Three focused
synthetic integration regressions specify fixed-snapshot parity, the independent
two-wedge equation, finite-gas feedback into native momentum, original adhesion
receipts, combined thermal balance and transactional refusal. Their native
execution remains unverified in the current environment.

`MaterialFreeSurfaceEnthalpy3::with_equilibrium_wetting_patches` resolves
explicit, fixed substrate chemistry/finish regions on that same planar wall.
Each `EquilibriumWettingPatch` names native wall-cell footprints and its own
original equilibrium-angle or reversible adhesion-work source. Surface A and
its fixed finish/history may differ; the original liquid surface, medium and
named gas environment must agree. Unlisted cells retain the explicit base
substrate. Overlaps, repeated cells, off-wall coordinates and the original
wall-cell budget refuse deterministically; ordering does not select a winner.

One source/physics binding is shared per patch, and an immutable shared index
routes cell history and temperature without cloning cards per cell. The
pattern is assigned once before physical time advances and before rate-law
or independent-threshold selection. Static Young--Dupre, subsequent Cox/MKT and sourced-gas refresh
all use the selected local original claims, actual liquid temperature and the
same native tension. Existing finite-atmosphere paths update the same outer
gas on every patch. Samples retain distinct original receipts; `wetting_source`
continues to identify the base substrate.

The spatial angle field drives the original capillary ghosts, pressure and
native populations. Static height-graph coupling rebuilds the same boundary
normal, area and Marangoni stress and retains column/corner admission. Patch
edges are resolved only at wall-cell scale: no subcell Cassie/Wenzel mixture,
chemical-edge pinning law, moving substrate, evolving coating chemistry or
surface-energy/heat transfer is inferred. Focused synthetic regressions
specify original-source/temperature routing, native momentum response,
homogeneous parity, order independence, molecular/two-phase Cox forward
relations, evolving-gas refresh and complete rollback on patch-source exit.

`MaterialFreeSurfaceEnthalpy3::with_equilibrium_wetting_hysteresis` explicitly
adds independent static receding/advancing claims from each original base or
patch card. The original equilibrium angle or Young--Dupre W(T)/sigma(T)
remains a separate reversible source; it must lie inside that card's current
threshold bracket at attachment and every thermal/native endpoint. Original
threshold identities stay pinned, with their own source receipts and no
fabricated speed axis. A source exit or loss of this ordering refuses the
complete coupled candidate.

Resolved accepted contact motion outside the declared numerical deadband
selects its direction's threshold. Inside the deadband, the actual geometric
angle supplies the bounded static reaction or triggers depinning at the
appropriate endpoint. Modeled speed is zero there, even when a tiny observed
speed has the opposite sign. Unresolved contacts use the equilibrium angle;
a zero-width bracket exactly recovers smooth equilibrium kinetics. The
original physical-VOF history is independent of imposed curvature ghosts.

Optional Cox, molecular and sourced-gas models consume that selected
microscopic angle and modeled speed. The actual observed speed remains in
their admission checks and diagnostics, including the observed and modeled
Cox capillary numbers. Original local T, sigma and gas transport refresh
through their existing owners on every selected patch. The same angle field
drives native curvature, capillary pressure and momentum. Shared patch owners
remain bounded and no source card is cloned per cell.

Signed static resistance is sigma*(cos(theta_Y)-cos(theta_branch)); its product
with modeled speed is a nonnegative per-length pinning-loss diagnostic.
Molecular kinetic loss remains separate. Neither becomes heat without an
owned contact-line measure and energy exchange. This is a first-cell geometric
boundary closure, not an exact contact-position constraint, experimental
calibration, moving-wall model or closed surface-energy budget. Assignment is
before time advances and after patch selection; constant or spatially uniform
tension is required. Native synthetic regressions specify sticking and both
depinning directions in wall-adjacent liquid mass, thermal/patch source routing,
independent Cox/MKT relations, finite-gas feedback, smooth zero-width parity and
transactional source/order refusal. These definitions passed syntax and
source/API review; native execution remains unverified in this environment.

`free_surface3_from_cells_with_uniform_thermal_capillarity` extends that same
owner to time-varying sigma(T) in the exact zero-Marangoni domain. Every actual
wet-cell source value must agree exactly in SI units before lattice conversion,
at initial admission, every conductive thermal substep, and both native step
endpoints. No tolerance, averaging or suppressed tangential stress admits a
spatial mismatch. Complete source bundles and original selected claims remain
pinned. Native capillary-wave/density/pressure admission still applies to the
resulting field before publication.

The same local sigma(T) reaches native Laplace pressure and Young--Dupre
W(T)/sigma(T). Storage-free wet-interface ghosts refresh from adjacent real
liquid using the existing enthalpy extension; difference-first interpolation
preserves a uniform donor temperature exactly without adding mass or heat.
Genuinely dry wall ghosts retain their declared reference context. Uniform
capillarity needs no graph or local PLIC reconstruction and retains ordinary
native conversion/mass routes. It adds no tangential Guo load and no separate
surface-energy or heating model. Separately solved surface temperatures still
require their original geometric stress owner. Nonuniform tension remains
unsupported with an imposed wall angle and refuses atomically.

Three additional native regressions in the same equilibrium test module
specify homogeneous W(T)/sigma(T) heating with the actual pressure coefficient,
nonconstant-ratio hot-state angle, constant-sigma trajectory/receipt parity,
and complete rollback for local/interior-substep gradients, original tension
source exit, loss of partial wetting and cancellation. These are authored
numerical regressions; native execution remains unverified in this environment.


### Source-bound DC current on the free-surface thermal trajectory

`MaterialFreeSurfaceEnthalpy3::with_electrical_core` retains a bounded rectangular
native current domain, ideal electrode drive and explicit contact laws. Every
selected cell must remain a fully occupied native `Fluid` cell; partial VOF
current geometry refuses before publication. Core exterior faces impose declared
zero normal current. They neither infer surrounding material insulation nor
create mechanical walls. Electrical regions and contact indices are local to
the declared core.

The same complete liquid T/shear bundle supplies `ResolvedConductorState`;
exactly one scalar conductivity or resistivity source is required, and changing
that source basis is a data change. All original requirements, query meanings
and pinned source receipts remain. The original voxel current solver supplies
half-cell/contact I^2 R to the existing thermal observer at each actual accepted
interval. Diagnostic endpoint solves consume budget but deposit no heat. The
report independently integrates electrode power and deposited Joule power,
admits both instantaneous W and integrated J defects, and retains applied and
thermal-endpoint current/source epochs.

Every existing conduction, finite-wall, bath, finite-atmosphere, rheology and
wetting path automatically uses the bound current. Instantaneous-heat/pure-flow
paths refuse an attached drive rather than omitting it. The original thermal
timeline updates conductivity, viscosity, film and contact sources; its endpoint
h feeds the original collision and VOF enthalpy transport. Current/source/budget/
cancellation failure or candidate core dewetting leaves liquid, walls and
atmosphere unchanged.

Scope: zero imposed magnetic field, quasi-static isotropic Ohmic transport,
fixed complete current-cell geometry, explicit thermal-current feedback and the
existing first-order thermal/flow split. No cut-cell conductivity interpolation,
motional EMF/Lorentz work, electrochemistry, phase-density closure, surface-work
heating or total kinetic-energy claim. Refine time for nonlinear electrical
feedback. Native Joule heat is external heat to this thermal system; it is
counted once and compared separately with measured electrode supply.
`combined_thermal_balance_defect_j` in the finite-gas report includes electrical
Joule heat among external inputs.

Five focused regression definitions cover independent R(T)/current/heat,
conductivity/resistivity equivalence and actual population response, coupled
liquid/wall/gas storage, late source/budget refusal and partial geometry/drive
admission. All five changed Rust files passed syntax parsing, and independent
scalar recurrences were evaluated. Native compilation/tests were not run:
Rust/Cargo and DSR/RCH were unavailable in the editing environment.


### Imposed-field motional current and Lorentz forcing of moving liquid

`MaterialFreeSurfaceEnthalpy3::with_imposed_magnetic_field` extends the retained
fully occupied electrical core with a prescribed stationary Cartesian B field.
It calls the original moving-conductor face operator: the same half-cell
interpolation produces u cross B voltage and its transpose Lorentz force.
Accepted native `sample_at` velocity already includes the Guo half-force
correction; no second correction is added. Temperature/shear-dependent
conductivity or resistivity remains in the original complete pinned liquid bundle.

Every actual thermal interval deposits native conductor-frame I^2 R once.
Electrical power and accumulated energy independently satisfy electrode supply
= Joule heat + signed sampled Lorentz work to the declared W/J tolerances.
Endpoint current diagnostics deposit neither heat nor impulse. Integrated native
cell forces become one interval-average load through the original rho0/dx/dt
scaling and Guo collision; original gravity and source buoyancy remain.
Endpoint reframing restores the original nonmagnetic force field while preserving
physical velocity and force-corrected stress, so continuation does not inherit a
stale magnetic half-force.

The explicit current/thermal/flow split holds accepted velocity during thermal
current solves. `sampled_lorentz_work_j` is that sampled electrical/mechanical
exchange, not an exact LBM kinetic balance. `flow.magnetic_motion` separately
observes actual core control-volume kinetic change and the largest endpoint
velocity change, including advection, pressure, viscosity and other forces.
A caller-declared velocity-change threshold refuses excessive lag before the
liquid/wall/gas transaction commits; it is an error-control observation, not a
rigorous error bound or balancing heat. Neither kinetic change nor a work
discrepancy is thermalized.

Scope remains low induction, imposed static B, scalar Ohmic transport, fully
occupied fixed current geometry and declared current boundaries. No partial-cell
current closure, magnetic field storage, induction, Hall force, magnetization,
Maxwell traction or global kinetic-energy closure is claimed. Four focused
regressions cover motional braking, force reversal, conductivity substitution,
zero-B parity with original gravity, accepted-state continuation and late
split-gate rollback. The seven changed Rust files passed Tree-sitter syntax
parsing, and independent scalar current/force/energy references were executed.
Native compilation, Cargo tests and rustfmt remain unexecuted in this environment.

### Isolated current loops in the moving material liquid core

`MaterialFreeSurfaceEnthalpy3::with_isolated_electrical_core` selects the existing
isolated moving-current owner explicitly. The original electrode-bound entry
retains its previous behavior. The isolated input must have no electrodes and a
finite nonnegative per-cell current-continuity tolerance. Each current epoch
records the selected `FreeSurfaceCurrentBoundary`; no empty drive is silently
interpreted as isolation.

The native owner gauges each conducting component without replacing its gauge
voxel with an ideal electrode. All original resistivity, contact resistance,
Joule heat and reciprocal Lorentz force remain present, including at the gauge.
Every cell's Kirchhoff balance is checked. The same accepted physical velocity
and original T/shear material bundle drive each thermal substep's current, heat
and force. With no electrical supply, conductor-frame Joule heat is balanced by
negative sampled Lorentz work. Integrated forces enter native Guo forcing once;
accepted-state continuation does not retain a stale magnetic half-force.

Nonuniform imposed B can therefore brake and heat an electrically isolated
moving liquid. Material conductivity/resistivity substitutions change both
responses. Uniform translation through uniform B instead develops a potential
gradient with no circulating current; omitting B also supplies no current or
invented leakage. No additional empirical magnetic loss is added.

The core remains a fixed, fully occupied Cartesian current region with explicitly
declared insulating/symmetry exterior boundaries. This is not a partial-VOF-cell
electrical closure or time-dependent field-induction model. Original source,
current/thermal budget, power/energy, velocity-change and atomic liquid/wall/gas
gates remain active. Sampled mechanical work remains distinct from the observed
LBM control-volume kinetic change; no discrepancy becomes compensating heat.

Synthetic regressions in `material_free_surface_enthalpy/conduction/electrical/`
`magnetic/isolated.rs` independently reconstruct a four-edge current loop, every
cell's charge balance, gauge-cell heat, magnetic impulse and Joule/work exchange.
They cover equivalent source bases, conductivity substitution, B reversal,
actual native motion and continuation, uniform-field polarization, zero-B
isolation, forbidden physical electrodes and late current-budget rollback.

### Original material interfaces on the moving-liquid current timeline

`MaterialFreeSurfaceEnthalpy3::with_material_electrical_interfaces` binds one
original ordered interface card and complete selected source bundle to each
authored contact. The interface condition callback receives the same adjacent
specific enthalpy, temperature, bulk source identities and declared texture
frames used by current. It supplies an explicit physical query, connectivity
and closed-contact heat partition. Resistance is refreshed by the original
`ResolvedElectricalInterface::at_query_point` owner; copied coefficients or a
freshly available fallback claim cannot replace the retained source selection.
Actual surface material/process/phase identities and texture frames must match
in order. Unsupported original source data refuses even for an open interface.

The resolved contact or open face enters the existing stationary or moving
current topology owner, including isolated component gauges. One current
solution supplies bulk/contact Joule heat and reciprocal Lorentz force.
Opening a face interrupts electrical loops without deleting liquid, mass,
enthalpy, thermal links or retained physical velocity. Ordered interface heat
is already included in native cell power and is deposited once. Applied and
thermal-endpoint current reports retain their exact interface/source epochs;
endpoint diagnostics add no heat or impulse.

All existing liquid conduction, finite wall/atmosphere, radiation, rheology and
wetting steps retain the binding. Callbacks must be pure across retries and
endpoint queries. Original source, current, thermal, work and split-motion
gates admit the whole liquid/wall/gas transaction or reject it unchanged.
Geometry remains a fixed fully occupied current core; interface temperature,
pressure, history, partition and switching conditions are explicit caller
models. There is no new interface energy store, contact mechanics, induction,
partial-cell electrical model or switching-event locator. Discontinuous laws
must resolve their timing by declared substep refinement.

Three synthetic regressions in
`material_free_surface_enthalpy/conduction/electrical/interfaces.rs` cover
an independent temperature-dependent series circuit and ordered heat recurrence,
actual native motion under source substitution, static/moving current agreement,
an isolated induced loop opening during heating and reclosing on retained flow,
and atomic refusal of reversed surface assignments or late original-source
failure despite an available fallback.

### Source-driven planar thermocapillary flow (MR26/MR27 increment)

`MaterialFreeSurfaceScaling::free_surface3_from_cells_with_thermocapillarity`
binds the original liquid, caloric chart, optional buoyancy/rheology and one
explicit `PlanarThermocapillarity` envelope. The full original material bundle
retains a positive surface-tension claim in N/m; only a source-declared constant
or the existing law's absolute-temperature response axis is admitted. Every
thermal evaluation pins the same claim, and the actual accepted local caloric
temperature drives both the native Laplace-pressure coefficient and tangential
surface stress. No duplicate source epoch or material-name rule is introduced.

The geometry reduction is a Cartesian tank/slab with a single nearly planar
surface in every column, neutral stationary walls and no internal walls,
overhangs, disconnected liquid, multiple partial layers or unlocated carry.
The declared surface-cell layer cannot change during the admitted trajectory.
Actual raw VOF fill reconstructs each column height. A caller-declared departure
smaller than half a cell and a positive surface-fill floor gate both endpoints.
Centered tangential sigma differences (one-sided at actual tank walls) implement
the declared flat-interface stress. Each projected patch has area dx² and its
surface cell has occupied volume fill*dx³; their ratio supplies the physical
force density without a diffuse-interface normalization. Normal gradients do
not become tangential forces. Normal capillary pressure remains solely in the
existing smoothed-fill reconstruction with that cell's original sigma.

The extra traction composes with the original Guo gravity, source buoyancy and
any integrated Lorentz load. Source-frame changes retain physical velocity and
force-corrected stress without a second impulse. Original thermal substeps,
mass routes, enthalpy transport, finite wall/gas histories and endpoint source
admission remain transactional. A vanishing gradient restores the original
force frame. Constant-source ordinary curved trajectories keep their original
arithmetic; this opt-in does not relax their source restrictions.

Scope is a first-order held-force split with cell-centered source temperature
and a caller-admitted flat-interface approximation. It is not arbitrary curved
thermocapillarity, a resolved interface-temperature law, surfactant transport,
contact-line calibration, capillary heat, surface entropy or a conservative
surface/total-mechanical-energy model. Refine time and interface resolution;
no experiment-backed accuracy or unconditional stability is claimed.

Focused G1/G3/G4 regressions cover all normal-axis/sign choices, independent
traction magnitude and occupied-volume weighting, normal-gradient exclusion,
original Laplace-pressure arithmetic, constant-source parity, source-only
native flow direction, heat/mass/replay and transactional refusal/retry.
Test targets are the native `d3q19::freesurface3::thermocapillary` unit module
and `material_free_surface_enthalpy`'s
`conduction::capillarity::thermocapillary` module with `enthalpy-coupling` enabled.
The broad MR26/MR27 deliveries remain open.

### Hydrodynamic adhesion from static hysteresis (MR26 increment)

`MaterialFreeSurfaceWetting3::with_cox_voinov` and
`MaterialFreeSurfaceEnthalpy3::with_cox_voinov_wetting` add the single-liquid
Cox relation to an existing sourced static advancing/receding binding. The
through-liquid apparent angle satisfies
`g(theta_app) - g(theta_micro) = (mu U / sigma) ln(L/l)`, using the full angular
integral `g(theta) = integral_0^theta (u - sin(u)cos(u))/(2sin(u)) du`, rather
than extending the small-angle cubic approximation to hydrophobic surfaces.
Reference: Chan et al., *Journal of Fluid Mechanics* 900 A8 (2020), equations
(1.1)-(1.3), <https://doi.org/10.1017/jfm.2020.499>.

The caller supplies the effective microscopic cutoff, apparent-angle observation
scale, maximum capillary number (at most 0.01) and angular domain within
[1e-4, 3] radians. The cutoff is not silently interpreted as a measured Navier
slip length. A shared 2048-interval fixed Gauss/Hermite table is constructed
once with cancellation checkpoints; each subsequent prediction uses bounded
arithmetic and table lookup, without allocation or transcendental evaluation.
Static source receipts stay at zero contact speed. Measured speed-dependent
angle curves refuse, preventing duplicate rate effects. Predicted angle and
signed capillary number are explicit model outputs, separate from the claims.

Both native flow paths consume the angle in their original curvature and
pressure reconstruction. The thermal path uses each liquid cell's actual
source viscosity, including conduction substeps and advected endpoint states;
dry numerical ghosts retain the reference coefficients. Dynamic planar wetting
also composes with the existing uniform thermal-capillarity owner: source
`sigma(T)` then drives BOTH native Laplace pressure and the Cox capillary
number at each admitted thermal and advected endpoint. Every actual wet-cell
tension must still agree exactly; unresolved spatial gradients refuse rather
than silently discarding Marangoni stress. Without that owner, tension remains
the native constant. Non-Newtonian rheology, missing tension, high capillary
number and absent solutions within the declared angular domain refuse before
publishing candidate flow or heat. No angle is clipped to conceal a film
transition. The existing pinning deadband and first-cell speed estimate retain
their prior limitations. Resolved gas motion, slip, entrainment, contact-line
calibration and conservative surface-energy accounting are not added.

`CoxVoinovLaw::try_new_with_viscosity_ratio` additionally admits a caller-declared
fixed ratio `M = mu_outer/mu_liquid`, using the full two-viscous-wedge angular
function in Chan et al. (2020), equation (2.13), with
`g'(theta,M) = -sin(theta)^2/(3 F(theta,M))`. The original constructor retains
the inviscid outer phase. Both liquid and outer-phase capillary numbers obey
the declared low-Ca bound. The table clusters nodes at both angular endpoints
to resolve narrow wedges and viscosity contrast, and large-ratio evaluation
uses reciprocal scaling. Runtime predictions still require no quadrature,
allocation or transcendental evaluations.

The ratio is explicit model data in every `CoxVoinovSample`, separate from the
original interface claims and liquid viscosity receipts. The same boundary
feeds native isothermal and thermal flow. In thermal use a fixed ratio means
the modeled outer viscosity varies proportionally with the local liquid
viscosity; independent gas temperature/viscosity evolution is not inferred.
This is an outer-wedge contact-angle closure, not a second-fluid momentum
solver. G1/G3 checks cover independent two-phase quadrature, the analytic
right-angle relation, phase exchange, small-angle resolution, native flow
changes and thermal coefficient use. G0/G4 tests retain transactional refusal
when either phase leaves the capillary-number envelope.

`MaterialFreeSurfaceEnthalpy3::with_cox_voinov_gas_wetting` instead binds the
original resolved thermochemical gas transport state. Each contact sample uses
`M = mu_gas(T_gas,p,composition) / mu_liquid(T_liquid)` with separate phase
temperatures and source receipts. Changing gas transport data therefore changes
the apparent angle and its subsequent native capillary response without a
material-specific fluid branch. The original fixed-ratio API is unchanged.
An explicitly bound gas snapshot remains fixed during ordinary liquid-only
steps. Finite-atmosphere conduction and solved-interface entry points consume
their actual evolving gas parcel, including callers that already change its
composition through the existing reaction owner. They admit the same original
transport claims at thermal source/endpoint checks, then pass that accepted gas
into the original native wetting preparation and advected endpoint acceptance.
Gas heat and composition remain with their existing owners; no second gas
momentum or thermal model is introduced.

`wetting_outer_gas` exposes the accepted outer source. Each wall sample retains
the same immutable gas receipt through an `Arc`, alongside its distinct local
liquid viscosity, applied ratio and interface claim bundle. The declared static
microscopic-angle claims and their fixed context remain authoritative; this
hydrodynamic connection does not invent a new equilibrium contact angle for a
different chemical environment. Surface tension follows the same native
constant or explicitly uniform thermal-capillarity owner.
Both phase capillary bounds, source support, native transport and cancellation
must pass before any liquid, finite gas or wall candidate publishes.

The native thermal-wetting regressions verify fixed-sigma trajectory equivalence,
heated source tension shared by pressure and moving Cox response, and complete
rollback of flow, enthalpy, tension, source history and time on a spatial-gradient
or source-domain refusal. This composition does not admit the incompatible
planar/local/height-graph thermocapillary geometry owners for dynamic wetting,
or add contact-line heat, a surface-energy balance or new calibration data.

The original reference-ratio table is reused exactly when ratio bits match.
Changed moving ratios first invert the SAME angular relation over only the
actual angle increment, without allocating another table or retaining a ratio cache.
A safeguarded scalar Newton solve uses a fixed stack and adaptive 4/8-point
Gauss quadrature, bounded by 64 inverse iterations, 24 subdivision levels and
32,768 angular-slope evaluations. Both sides are scaled by `max(1,M)` before
integration, retaining representable response at extreme viscosity contrasts.
Quadrature discrepancy and roundoff must fit the local angular tolerance.
If the local inverse is unresolved, the original bounded table is constructed
once for that sample, preserving its existing support for strongly receding
motion near small angles. This fallback retains the original table's accuracy
and physical-domain checks; the stricter local tolerance is not a global
certificate. There is no persistent ratio cache. Cancellation applies to changed,
unchanged and zero-speed sourced queries. Pinned/zero-speed samples retain their
original stationary behavior. The focused varying-ratio matrix verifies fewer
than one quarter of a rebuilt table's angular-slope evaluations per sample;
this is an operation-count bound for those cases, not a general timing claim
or a bound for fallback cases.
These evaluations remain bounded by the existing wall-cell and thermal-substep
limits. This supplies a two-viscous-wedge
closure at the declared local phase conditions, with the existing geometric
contact-speed and scale-separation assumptions; it adds no resolved gas motion,
contact-line heat, surface-energy theorem or calibrated entrainment model.

G1/G3 tests compare full-angle predictions against independent Simpson
quadrature, recover the cubic small-angle limit, check length-scale similarity,
and exercise actual capillary flow from static source data. Thermal integration
checks local viscosity and unchanged enthalpy accounting; G0/G4 tests cover
double-counting refusal, model exits and transactional state preservation.
Run `cargo test -p fs-lbm --lib cox_voinov`,
`cargo test -p fs-lbm --test material_free_surface3 dynamic_wetting`, and
`cargo test -p fs-lbm --features enthalpy-coupling --test material_free_surface_enthalpy conduction::capillarity::wetting`.

### Solved liquid/gas/radiative surface temperature (MR19/MR20 increment)

`MaterialFreeSurfaceEnthalpy3::step_with_surface_exchange` solves one shared
zero-storage interface temperature for liquid conduction, a declared gas film
and gray opaque radiation to black isothermal surroundings. The companion
`step_with_finite_gas_surface_exchange` attaches the existing fixed-mass,
fixed-volume NASA-9 atmosphere. These entrypoints run on the original thermal
substeps before the same native rheology, VOF enthalpy and momentum update.

The caller supplies currently exposed face area, liquid-centroid distance,
gas-film thickness, emissivity and radiative-bath temperature. Liquid half-cell
conductivity is sampled at its current bulk temperature; gas conductivity is
resolved at the mean of solved surface and actual gas-bulk temperatures, with
the gas's current EOS pressure and original composition/source conditions.
The surface balance is `Gl (Tl - Ts) + Gg (Tg - Ts) + Gr (Tr - Ts) = 0`, where
`Gr = emissivity sigma A (Ts + Tr) (Ts^2 + Tr^2)` is the exact fourth-power
secant. The finite-volume liquid half-cell approximation remains explicit;
there is no inferred convective coefficient or parallel radiation path that
bypasses its resistance.

A bounded safeguarded scalar solve requires both temperature and power
residual tolerances. Its trials intersect the current pure-liquid caloric
interval, connected NASA temperature support and original gas transport
claim domains. A hot radiative reservoir is never queried as gas, so a
supported interior root can exist when that reservoir is outside the gas
source range. Every film trial still passes through the original complete
source evaluator; invalid interpolation trials backtrack within the same work
budget, and switched source claims cannot supply an accepted root. Neither an
unsupported phase/film solution nor exhausted work publishes partial state.

With `S = Gl + Gg + Gr`, eliminating the surface produces nonnegative pair
conductances `Gl Gg / S`, `Gl Gr / S`, and `Gg Gr / S`. The resulting pair
powers conserve exchange independently of the bounded scalar residual. Liquid
and finite-gas CFL rows each include all incident conductances, together with
the liquid's original conduction/contact losses. In particular, the finite
gas limit includes its radiative-reservoir path through the shared surface.
Every accepted interval stores the exact opposite gas-supplied joules in the
original internal-energy owner. Radiation can therefore warm initially
equal-temperature liquid and gas, with only the radiant supply external to
their union. Separate surface reports already enter thermal/flow heat once;
endpoint diagnostics supply no additional energy.

This is a first-order explicit thermal/flow split with a quasi-steady film,
bulk-sampled liquid half-cell conductivity and prescribed view factor one.
It requires temporal and spatial refinement for accuracy. No evaporation,
surface heat capacity, gas momentum/pressure work, participating-gas radiation,
moving-solid handoff or liquid phase conversion is supplied. The existing
resolved-cell radiation entrypoint retains its cell-temperature semantics.
Temperature-dependent planar and height-graph thermocapillarity can consume
this surface field under the complete-column admission below. Source-declared constant capillarity
remains supported, including its planar limit. Original resolved-cell
thermal/thermocapillary entrypoints retain their behavior until a trajectory
explicitly enters the solved-surface stress mode.

Six focused regressions in `material_free_surface_enthalpy`'s
`conduction::environment::surface` module cover independent nonlinear film
roots and heat, NASA gas storage/pressure, small-gas radiative stiffness,
supported roots under hotter out-of-domain radiant baths and partially
supported gas-source interpolation, and actual native
viscosity/population changes, and atomic liquid/wall/gas refusal/retry.
Their source coefficients are authored synthetic data, not experimental
validation. The broader MR19/MR20 multiphysics deliveries remain open.
The existing planar capillarity regressions additionally check constant-source
native parity and refusal of incomplete solved-temperature stress fields.

### Solved interface temperature drives capillary flow

For temperature-dependent planar or height-graph surface tension, the shared
heat-boundary entrypoints require exactly one full outward patch per admitted tank
column. Missing columns, partial patches, duplicate columns and side-facing
patches refuse. The original sigma claim alone is queried at each solved
interface temperature, with its original selection, complete source context
and source-domain checks. Bulk material bundles keep their actual local
caloric temperature, enthalpy, viscosity, rheological state and density.
Neither thermal expansion nor a temperature substitution changes reference mass.
The original complete bulk material bundle still admits its sigma claim at
bulk temperature before the surface override; this conservative source-domain
requirement remains. The interface query has its own retained Ts receipt and
is the coefficient actually used for native pressure and traction.

At the end of the original thermal stage, this same sigma(T_surface) field
drives both the native tangential Marangoni traction and Laplace-pressure
coefficient. Each physical column's surface value extends normally onto its
native wet and zero-mass ghost lanes within the admitted planar or height-graph
columns. The local detached-interface reconstruction has no column extension.
The original geometry, occupied-volume traction scaling,
Guo force-frame conversion, buoyancy and optional Lorentz loading remain
with their existing owners. No force contribution is added twice.

After the native mass/enthalpy advection, the shared scalar surface closure
resolves again against actual accepted-candidate liquid temperature/rheology
and the finite gas's accepted thermal T/p. Original gas-source and sigma
support still gate this diagnostic epoch. The resulting native endpoint
updates pressure and traction without another heat interval, source debit or
gas-energy change. Surface reports retain this native endpoint separately
from the thermal endpoint that drove the step. Caller-declared normal conduction
distances are held over the interval. Graph patch owners, normals and areas
are reconstructed after motion; planar height/layer admission and the graph's
strict fill/slope envelope gate their respective native configurations.
Callers resupply current patches on continuation. This remains a first-order split requiring
time and interface-resolution refinement.

Once this mode accepts a step, subsequent steps must provide its complete
solved heat boundary; they cannot silently revert to bulk-temperature stress.
The persistent mode flag, native populations/forces/sigma, enthalpy/source
states and finite gas/walls publish together after all gates. Failed mode
entry does not latch the requirement, and source, geometry, budget or
cancellation refusal preserves the prior accepted state. Surface free energy,
capillary-work heating, surface entropy and conservative total mechanical
energy remain outside this model.

Focused regressions in the thermocapillary test module independently compare
nonuniform gas baths, series surface temperatures, original heat storage,
sigma gradients, actual native traction and capillary pressure; preserve a
uniform solved-surface/no-Marangoni control despite a bulk temperature
gradient; and check surface-only source exits, incomplete fields and mode
continuation/rollback with finite gas. Their synthetic source laws test
implementation and physical coupling rather than experimental accuracy.

## Finite-mass planar D3Q19 wall coupling

`BoundaryGrid3::step_with_finite_mass_walls` and its pooled counterpart advance
fluid populations and the tangential velocities of two opposite complete planar
walls in one transaction. `FiniteMassWall3::mass` is the total mass of one face
in lattice units (`dx = dt = 1`), not an areal density. The other two axes must
be periodic; voxel solids, open faces, and additional walls refuse. The grid's
boundary velocities hold the accepted endpoint wall state.

The existing collision runs once. Two face-only reductions use its actual
post-collision populations: the first obtains the affine tangential reaction
`J(u) = J(0) - D*u`, with `D = sum(rho_post)/3`; a direct midpoint solve couples
that reaction to finite wall inertia and the declared external impulse. The
second evaluates the shared link-wise momentum-exchange owner at the solved
midpoint velocity. Streaming uses that same velocity. This avoids repeated
fluid solves, full-state cloning on the serial path, and explicit added-mass
instability for light walls. The ordinary prescribed-wall step is unchanged.

`FiniteMassWallReceipt3` reports before/midpoint/after velocities, measured
hydrodynamic impulse and moment, wall kinetic-energy change, and hydrodynamic
and external work. Tangential fluid-plus-wall momentum includes the declared
external impulses. Wall work uses the midpoint velocity and closes its
kinetic-energy increment to floating-point roundoff. External impulse acts at
the geometric face center. Normal translation and all rotation remain
constrained: the receipt explicitly reports the support impulse and angular
impulse about the caller's moment origin. Periodic transport contributes angular
momentum flux, so this is not a closed-domain global angular-momentum claim.

Geometry, finite positive mass, finite inputs, collision state, load arithmetic,
and midpoint/endpoint low-Mach wall speeds are checked before publication.
Pooled cancellation preserves the old populations and wall velocities; the
unpublished destination follows the existing pooled stream path. This capability
has fixed planar geometry: it does not cover rotating immersed bodies, changing
voxel occupancy, fresh-cell reconstruction, or arbitrary moving interfaces.

Focused serial regressions cover every normal axis, an independent implicit
impulse solution, very light walls, common translation, passive relaxation,
fluid-plus-wall momentum, external work/support reactions, and refusal without
publication. A separate pooled regression covers multi-group replay and
cancellation/refusal. The broad moving-boundary bead remains open.

## Finite-inertia immersed D2Q9 translation and rotation

`Grid::collide_into_with_partial_saturation_body` and
`Grid::step_with_partial_saturation_body` couple one planar rigid velocity
state to the original Noble-Torczynski partial-saturation operator. The caller
supplies total mass, polar inertia about the center of mass, the current rigid
velocities, frozen center/fractions, and external impulse/couple in lattice
units. Linear external impulse acts at the center of mass. Both rigid
velocities and actual fluid populations advance; pose/fraction updates stay
with the geometry owner.

At each coupled cell, the original collision supplies solid impulse
`B rho (u - G v_mid)`, where `G v = V + omega cross (x - center)`. A single
3-by-3 positive-definite `fs-la` solve determines the rigid velocity increment
with inertia `diag(m, m, I) + sum(B rho G^T G)/2`. This retains off-center
translation/rotation coupling and avoids dividing an explicitly sampled load
by a small body mass. One ordinary collision evaluation at that midpoint
then supplies the measured reciprocal linear/angular impulse and solid work.
The receipt reports body kinetic-energy change, midpoint work, and the
momentum residual using the actual representable endpoint
velocity differences. There is no repeated fluid solve.

The coupled entry admits BGK/Guo fluid/wall domains. It checks
positive finite inertia, the matrix response, and initial/midpoint/endpoint
body point speeds before publishing body state or caller scratch. The
collision exchange conserves fluid-plus-body linear momentum and total
orbital-plus-spin angular momentum after accounting for the imposed fluid
force and external body loading, to floating-point roundoff. Streaming
can additionally transfer momentum to stationary walls or angular momentum
through periodic seams; those are distinct from immersed-body exchange.

This is a frozen-geometry velocity-coupling substep, not a complete moving-body
solver or a fluid total-energy model. The inertia belongs to the solid degree
of freedom of this discrete fluid/solid model. Geometry-derived displaced- or
interior-fluid mass/inertia corrections, overlapping bodies, pose integration,
free surfaces and physical FSI validation remain outside this entry. Tests in
`tests/partial_saturation_body.rs` independently reduce an off-center point
mobility, sum actual population/body momentum, verify passive light-body
collision energy and replay, and exercise common rigid motion, the zero-solid
limit and atomic early/late refusals.

The physical Guo velocity supplies both the finite-inertia solve and the
prescribed solid-collision law. `PartialSaturationForcing2` sums the incoming
held force `rho*g + fext` on all active cells, its moment about the requested
origin, and its work at mean physical pre/post-collision velocity. Interior
fluid is retained, so this force is independent of solid fraction. Its work
does not enter the body's hydrodynamic work or the external solid-work account.
The force account uses the same held source for both collision endpoints;
changing a force field or a density-dependent source between steps requires
the usual Guo half-source accounting. No fluid total-energy or global time-order
claim follows from the local mechanical work receipt. Zero force retains the
original arithmetic and zero occupancy retains native forced BGK bits.

The forced cases in `tests/partial_saturation_body.rs` independently reconstruct
the Guo polynomial and every blended population with nonzero stress, then close
repeated fluid/body/wall linear and angular momentum, force work and body work.
They also reject nonfinite sources, combined-force overflow and a late native
collision overflow without publishing body state, populations or caller scratch.

### Moving disk geometry on the same finite-inertia operator

`core2::moving_disk::MovingDisk2` advances a disk's center, material angle,
translation, spin and actual native fluid populations on each lattice tick.
A half-drift places the disk at the collision configuration; one original
finite-inertia solve supplies the reciprocal impulse; the accepted velocity
supplies the final half-drift. Each free drift is parallel to its own linear
momentum and therefore preserves orbital angular momentum. The step receipt
also reports the collision center and the external angular impulse about the
fixed lattice origin. The existing midpoint work and kinetic account remain
those of the original velocity-coupling operator.

Circle/square intersection areas provide analytic partial-cell fractions,
including subcell disks and exact tangencies. One reusable fraction buffer
clears only the previous bounding box and fills the new one. Interior and
exterior cells avoid sector integration. Native populations retain their mass
and momentum as occupancy changes; no fresh-cell population reconstruction or
additional full fluid candidate is introduced. Body mass and polar inertia
are explicitly the additional solid degree of freedom: retained interior
fluid remains part of the discrete model. This is not a physical-density or
displaced-fluid inertia closure.

The whole swept bounding box of each half-drift must contain fluid cells and
remain within the interior domain halo. The endpoint is admitted before native
streaming or accepted body publication; a late refusal changes only reusable
scratch, which the next call rebuilds. Fixed radius/inertia, one disk,
BGK/Guo and lattice `dx = dt = 1` are the supported scope. Multiple bodies, wall
contact, periodic body seam crossing and free surfaces require further owners.
Native outer-wall reactions and periodic fluid angular flux remain distinct
from immersed exchange. No global temporal order, physical drag calibration or
complete moving-boundary qualification is inferred.

`tests/moving_disk.rs` checks analytic translated area, common translation over
many actual steps without stale-cell coupling, independent fluid-plus-body
linear/orbital/spin balances before outer-boundary transport, imposed force
through changing occupied cells, separate fluid and body work, bitwise replay,
and late-endpoint rollback with scratch recovery. The broad moving-boundary
bead remains open.

### Physical dry-disk momentum on the native D2Q9 operator

`core2::moving_disk::DryDisk2` retains original positive physical dry mass and
polar inertia at light, neutral and heavy density. It reuses analytic circle
fractions, the native finite-inertia BGK/Guo collision and its two free
half-drifts around collision and streaming. A bounded three-impulse solve
closes `J_c = Delta Pi_inside,physical - J_inside,F` using actual candidate
populations and endpoint geometry. No nominal rigid-interior mass or inertia
is subtracted from the body.

D2 acceleration is density dependent: `F = rho*g + fext`. Full-fluid and
collision-midpoint-occupancy source accounts therefore use
`(F_old + F_streamed)/2`, recomputed inside every correction candidate.
Authored body impulse/couple work, fictitious-fluid correction work, native
held-source work, wall reaction and periodic angular transport remain distinct.
Both center and requested fixed-origin correction residuals must converge;
independent dry-body and exterior-fluid momentum balances also gate acceptance.

A shared native stationary pull observes streaming without allocating a second
fluid grid and preserves existing pull/reflection arithmetic. Fluid, dry body
and angle publish together after final checks. Refused calls may change private
scratch, which a retry rebuilds. New geometry, numerical, evaluation-budget and
cancellation refusals are typed errors. Inherited low-Mach, numerical-inertia
and population assertions remain prepublication panics. Default work is capped
at 96 actual native collision evaluations; each numerical Jacobian costs three
additional probes, followed by bounded backtracking.

Seven new physical groups and the six existing moving-disk/eight existing
finite-inertia PSM groups pass against the actual native source closure.
They cover symmetric translation at density ratios 0.15/1/3, rising and
sedimenting disks, spin/orbital transfer about a translated origin, an
independently integrated density-dependent endpoint source, interior-only
forcing, wall/exterior bounce and rollback/exact retry after actual probes.
Full Cargo remains unexecuted where the constellation runtime closure is absent.

`disk_cell_moments` extends the original analytic circle/unit-square integral
with first moments relative to the cell center. The area-only specialization
retains the same arithmetic and skips moment work. Dry interior angular
momentum and midpoint interior-source torque integrate these exact cut-cell
moments against the piecewise-constant cell momentum/force. Full-cell native
collision and source torque remain unchanged; exterior quantities are their
cut-interior complements. Only midpoint cut cells enter an additional reusable
sorted cache, without another full-grid buffer.

The previously demonstrated oblique uniform-translation spin artifact is
removed. At dry density ratios 0.15/1/3, sixteen native ticks retain both velocity
components, spin and material angle within 3e-14, with zero Newton corrections
and one candidate per tick. Eight geometry/moving-disk and eight physical
dry-disk groups pass against rebuilt actual sources. Independent Cartesian
cap/corner integrals and full/half/quarter-disk moments, translated area/centroid
sums and area-only bit agreement check the geometric change. The retained
nonuniform-density/source, wall, angular-origin and rollback cases also pass.
Exact cut-cell first moments do not establish a general moving-boundary GCL,
higher-order fluid reconstruction or exact no-slip.
One fixed-radius disk, lattice dx=dt=1 and a fluid-only swept interior domain
are required. Multiple disks, body seam crossing, contact and free surfaces
remain outside this owner. The existing `MovingDisk2` retains its explicit
additional-inertia interpretation.

### Simultaneous dry-disk and prescribed-plane normal lubrication

`DryDisk2::step_with_lubrication_planes` adds up to 64 distinct prescribed
translating `LubricationPlane2` reservoirs to the original physical dry-disk
operator. The original `UnderResolvedCylinderPair2` normal resistance supplies
the positive gap floor and cutoff subtraction. Its dynamic viscosity is in
lattice units per unit depth; the caller matches the fluid viscosity and
calibrates the unresolved-film cutoff.

Every active film uses the same solved body midpoint velocity and frozen
half-drift geometry. The existing three-unknown correction closes
`J_supplemental = Delta Pi_inside - J_inside,F + sum(J_film)` with the
original positive dry mass and polar inertia. It does not add sequential kicks
or replace the native fluid collision. Empty or inactive plane sets retain the
previous arithmetic and exact accepted state.

Normal film impulses act through the disk center, with zero center couple and
their actual moment about the requested fixed origin. Each plane receives the
opposite impulse and moment. The returned `DryLubricatedDiskExchange2` separates
film body work, each prescribed plane's reservoir work and nonnegative film
dissipation from authored loads, the fictitious-interior correction, represented
exterior forcing and hydrodynamic exchange. No thermal destination is implied.

Start, midpoint and accepted endpoint gaps must remain positive. Intermediate
nonlinear probes may cross a plane so that a sufficiently strong implicit film
can prevent crossing in the accepted solution. Geometry and final receipt
checks, numerical budgets and cancellation gates run before fluid/body
publication. Input plane descriptors are immutable; the caller advances them
only after an accepted tick. Normal-only descriptors add no tangential film;
the optional cylinder extension below supplies sliding and spin resistance.
Contact, finite-plane recoil, plane rotation and a resolved fluid-wall boundary
remain separate owners.

The final actual-source closure passes six new film groups plus eight existing
physical dry-disk and eight moving-disk/analytic geometry groups. These exercise
exact inactive parity, simultaneous oblique films with independent resistance
and translated-origin balances, driven-plane work, density-dependent Guo source
and analytic interior moments, refusal/exact retry, and an implicit film that
prevents the uncorrected candidate's crossing. Full Cargo remains blocked where
the required sibling runtime dependency closure is absent.

### Cylinder/plane tangential translation and spin

`LubricationPlane2::with_tangential(cutoff_gap, minimum_gap)` explicitly adds
the leading small-gap Stokes resistance of an infinite circular cylinder
sliding and spinning beside a no-slip plane. The constructor alone retains
the original normal-only model. Tangent/spin scales are independent of its
normal cutoff and floor; all lengths and velocities are lattice quantities,
and impulses/work are per unit out-of-plane depth for one lattice tick.

With cylinder radius R, dynamic viscosity mu, relative center tangent speed
u_t and body spin omega, the additional law is

```text
h_eff  = max(h, minimum_gap)
zeta_t = 2*pi*mu*sqrt(2*R)*(1/sqrt(h_eff) - 1/sqrt(cutoff_gap))
J_t    = -zeta_t*u_t*t
tau_c  = -zeta_t*R^2*omega
```

The coefficient is exactly zero when disabled or h_eff is at/above cutoff.
A rationalized root difference retains positive drag immediately inside the
cutoff. Translation and spin have zero cross coefficients for this cylinder
asymptotic law. Pure sliding has zero center couple, pure spin has zero net
force, and rolling with zero contact-point slip still dissipates
zeta_t*(u_t^2 + (R*omega)^2). These are the Jeffrey--Onishi cylinder
coefficients, reproduced in [Reboux, Richardson and Jensen, equation (2.1)](https://personalpages.manchester.ac.uk/staff/oliver.jensen/OpenArchive/RebouxRichardsonJensen.pdf);
the sphere's logarithmic mixed resistance is a different law.

For example, a radius-2 cylinder may use independent normal and tangent scales:

```rust
let wall = LubricationPlane2::new(
    [0.0, 0.0], [0.0, 1.0], [0.0, 0.0], 0.125, 1.0 / 1024.0,
)?.with_tangential(1.0 / 32.0, 1.0 / 512.0)?;
```

At solve time, tangent cutoff must not exceed 0.02*R, a declared small-gap
admission policy rather than a finite-gap accuracy bound. Positive floors and
strictly positive actual gaps remain mandatory. This Estimated subgrid law
does not resolve finite-span end effects, cavitation, contact, a finite-inertia
wall or the destination of dissipated heat.

Every selected film evaluates BOTH its translation and spin at the same
solved native midpoint. Its center couple enters the existing three-unknown
dry-inertia solve, scaled by radius with the original angular unknown. The
fictitious-interior correction and its work exclude the complete film wrench.
`PlaneLubricationExchange2.native` remains the original normal-only receipt;
the added total body/plane impulses use the caller's same fixed origin.
The plane's free couple at its actual midpoint footpoint p is
-tau_c + (p-c) cross J_t. It accompanies the force reaction, so body plus plane
angular impulse closes even for pure translation and pure spin. Plane work
uses its prescribed translation; body work includes tau_c*omega.

Three focused native law groups and three new coupled groups pass alongside
all six existing normal-film groups. They cover the independent coefficients,
rolling passivity, origin/common-velocity transforms, separate floor/cutoff
behavior, a tangent-only active film outside the normal cutoff, simultaneous
oblique films, solved-midpoint spin, body/exterior-fluid/wall momentum and work,
exact inactive parity, typed tangent overflow/admission refusal, cancellation,
and exact retry after an exhausted evaluation budget. The original normal
pair-local assertion contract is unchanged. Validation compiles the complete
current D2Q9/moving/dry numerical owners with the original Guo and cancellation
leaves and genuine numerical dependencies; it is not a full workspace/runtime
or three-dimensional validation.

### Electrically insulating interfaces in moving conductors

`ElectricSetup::solve_with_cell_motion_topology` and
`solve_isolated_cell_motion_topology` admit explicit electrically open faces
alongside finite resistive contacts. One validated face map removes the same
interface from current continuity, motional EMF, Joule heat and reciprocal
Lorentz-force transfer, while preserving both cells and their bulk material,
mechanical and thermal laws. Existing motion entry points delegate with no open
faces and preserve their previous results.

Isolated-body potential gauges follow the actual conducting graph after cuts;
each separate lamination receives its own reference and can retain an independent
circulating current. Breaking a closed loop extinguishes its eddy current and
braking loss rather than leaving a fictitious path across insulation. Wired
components still require voltage references. Invalid or duplicate faces,
overlapping open/resistive definitions, unreferenced wired islands and cancellation
refuse without mutating the inputs.

This remains a frozen-state, imposed-field low-induction solve. It does not add
dielectric charging, displacement current, arcing, changing material occupancy or
self-induced magnetic storage. The original `electric_isolated_motion` and
`electrical_open_faces` targets cover slit opening/reclosing, independent
laminations, reversed motion, coupled power balance and input refusals.

### Insulating faces in coupled magnetic midpoint flow

`MaterialSlab<ThermalLbm>` and `MaterialSlab<ThermalLbm3>` expose
`step_magnetothermal_midpoint_with_open_faces` and
`step_magnetothermal_midpoint_isolated_with_open_faces`. Explicit adjacent
cell pairs interrupt electrical current while retaining the original fluid
momentum, material volume and thermal transport. The same face topology enters
every nonlinear trial and independent final source recheck. Isolated components
receive separate numerical gauges; gauge cells retain their resistance, Lorentz
force and heat. Wired components retain the existing voltage-reference
requirement. Empty openings preserve existing entry-point behavior.

Topology and thermal transport coefficients remain fixed during a call; bulk
resistivity is re-evaluated at the proposed temperature. Existing imposed-field,
low-Mach, source-work and transactional-refusal boundaries apply. The actual slab
regression opens a four-edge current loop into a zero-current tree, compares its
flow and heat with the unforced slab, and checks independently gauged 3D
conducting sheets against four-resistor motional-EMF power balances at the
accepted midpoint velocities. Explicit and local material-transport entry points
retain their prior topology behavior.

## Incoming-time material rheology on the reference 3-D grid

`lattice3::rheology` supplies material-rheology steps for `Grid3` and
`ForcedGrid3`. Each incoming velocity field is interpreted with its same-time
uniform or spatial half-force shift. The shared centered/quadratic-halfway-wall
stencil and shared material deformation-rate kernel evaluate
`gamma = sqrt(2 dev(sym(grad u)):dev(sym(grad u))) / dt` in 1/s. Rigid spin and
isotropic expansion/compression do not change the shear-viscosity query. The
original producer resolves one material response per fluid cell, its transport
view is frozen once, and `tau = 1/2 + (mu/rho0) dt / (cs^2 dx^2)` drives that
collision/stream step. Source objects and query coordinates, time and rate are
retained. Reference dx/dt/rho0 stay fixed; nonmatching density refuses. D3Q27
independent bulk/higher rates and all equilibrium choices survive; D3Q19 models
retain their configured independent bulk rate or existing bulk/shear tie.

This is explicit incoming-time generalized-Newtonian feedback, not a
self-consistent nonequilibrium-stress closure or implicit time integration.
No source-domain extension, viscosity floor/cap, material fit, heat deposition,
thermal evolution or hereditary memory is added. The preparation sweep caches
velocities in O(N) storage/work and queries sources in fixed serial cell order.
Failed source, scaling, collision or cancellation checks publish no step.
Time-varying force endpoints remain owned by the existing trapezoidal force
step and are not used to reinterpret incoming strain.

`Grid3::rheology_workspace` and `ForcedGrid3::rheology_workspace` reserve two
source/collision epochs and one velocity buffer. The corresponding
`*_with_material_rheology_workspace` steps and force-endpoint advances reuse
the same preparation kernel as the one-shot APIs, resampling every current
fluid cell rather than caching its previous viscosity. Successful steps reuse
all adapter vector allocations. Original material objects and callbacks,
native pooled scheduling and caller-owned force vectors retain their own
allocation behavior. Exact grid extents are admitted before source queries;
matching cell counts alone are insufficient. Source, scaling, native endpoint
and cancellation refusals preserve both the accepted flow and the previous
source report. Only native success permits the infallible report-buffer swap.

`Grid3::step_{serial,pooled}_with_channel_material_rheology_and_moving_walls`
bind that same reusable material preparation to prescribed native moving walls.
The admitted channel has one fixed halfway-wall pair, two periodic tangents,
no interior solids and a complete native probe containing both wall planes.
Each plane has one uniform tangential lattice velocity. The material strain
stencil uses those exact wall velocities, including the two-wall difference
for a one-cell gap; a common translation does not become artificial shear.
The same local molecular collision field and wall velocities drive one native
stream, whose actual impulse/work accompanies the accepted source epoch.
Wall admission precedes source queries; source, streaming, measurement and
cancellation refusals preserve both accepted flow and the previous source
report. The ordinary workspace supplies all adapter storage, with no wall map
or per-step allocation. This remains incoming-time generalized-Newtonian
rheology on fixed geometry; wall work is not deposited as thermal energy.

The moving-material regressions check analytic Couette strain and the original
Carreau-Yasuda material producer on D3Q19 central moments and D3Q27 cumulants,
including one-cell gaps, two-component shear, common translation and the
same-time half-force shift. Repeated changed drives compare serial/pooled
population bits and native wall-load receipts, retain independent rates and
reuse allocation addresses. Missing-wall and post-stream cancellation cases
retain the previous source epoch and permit an exact successful retry.

Five native exact-source serial tests exercise the actual grid owners and real
execution context: all seven Newtonian collision variants reproduce existing
population bits; force-centered Carreau queries ignore rigid rotation;
quadratic wall shear is recovered; original source/query ownership and rollback
are retained; and a CM power-law channel with index 0.75 reduces maximum
analytic-profile error from 0.723% at 12 fluid rows to 0.190% at 24. Every source
query stays outside the fixture's explicitly declared low-rate regularization.
The workspace regressions cover repeated serial/pooled source refresh across
D3Q19/D3Q27 collision families, buffer reuse, source and force-endpoint refusals,
post-stream cancellation and exact retry. The Carreau force-frame oracle and
power-law channel refinement also use the reusable workspace. No cross-ISA or
measured-throughput claim is added.

### Local capillary-interface conversion preserves detached drops

Both dense D2Q9 and sparse D3Q19 free surfaces redistribute each converted
cell's signed excess only to adjacent surviving interface cells after the full
flag cascade. A conversion no longer changes every distant drop's fill and
Laplace pressure through a global pool. Recipients follow fixed lattice-link
order, use fixed-size neighbor storage and take the final rounding remainder
locally. The redistribution cost depends on converted cells and their stencil,
without another whole-domain interface scan.

New D3Q19 constructors select `ConversionRedistribution3::ConservativeLocal`.
The named `Local` policy explicitly retains its former mass-only arithmetic.
`with_conversion_redistribution(LegacyGlobal)` explicitly reproduces historical
version-1 trajectories, including their nonlocal pooling. The frozen pour
sentinel selects this replay mode and retains its original four-quadrant hash;
new locality and thermal regressions exercise the corrected default. The
independent `FREE_SURFACE3_BIT_SEMANTICS_VERSION` is 8 for the current default;
version 2 introduced local redistribution, and version 3 corrects the density
and force handling of fresh-cell initialization described below. Version 4
uses the physical Guo frame in gas-pressure reconstruction as described above.
Version 5 retains actual positive density in fill fractions and actual total
mass in relative conservation diagnostics, without an absolute numerical floor.
Version 6 selects bounded physical momentum conversion for both binary-occupancy
and transferred-cell constructors. Ordinary and cell/link-pressure steps thereby
reuse the existing atomic conservative candidate without requiring an opt-in.
Explicit `LegacyGlobal` replay retains its historical conversion and
gas-pressure arithmetic, leaving the unchanged D3Q19 core, boundary and sparse
goldens untouched.

Both native free-surface owners compute interface volume as `mass/density` at
every admitted positive density. Rescaling both quantities by the same factor
therefore preserves fill, resolved/smoothed geometry, curvature and the fill
weights of mass exchange. The sparse relative mass diagnostic likewise uses
the actual initial inventory, so a small total mass cannot conceal a large
fractional loss. An unchanged empty inventory reports zero; creation from zero
or nonfinite inventory reports infinite defect. Initial zero-density wet cells
remain refused by their existing constructors. Both streaming paths also refuse
nonpositive or nonfinite endpoint carrier density before publishing populations
or mass, including at zero surface tension. A dry interface has zero tracked
liquid mass and positive carrier density; its fill and velocity remain defined.
These corrections do not change
the imposed gas-pressure datum, which remains reference density one, and do not
claim a physical low-pressure or rarefied-flow regime.

Focused G0/G3 regressions cover disk/sphere initialization and geometry under
power-of-two density rescalings, mass/volume conversion, zero-density refusal,
dry-interface readback, pressure-driven density collapse before publication,
and fractional mass defects across small, ordinary and large inventories.

`FreeSurface3::advance_body_acceleration` evolves prescribed acceleration from
its retained endpoint to the next one. Collision consumes the old `rho*a`;
the shared native stream, gas-pressure closure and interface initialization
use the actual next density and acceleration. No zero-time source reframe
removes the changed load's half impulse. Homogeneous constant-density bulk flow
receives `delta_u = (a_old + a_next)/2` per lattice step at any admitted density.
Equal acceleration endpoints reproduce the held-force trajectory bit for bit.
The gas-link-pressure variant composes this evolution with the original frozen
pressure field, and material `advance_gravity` entries retain the same SI units
and pressure reference. Native flow, conversions, force history and the next
acceleration publish together after all refusals/cancellation checks. This
adds no vessel dynamics, thermal work or total-energy closure. Focused G0/G1/G3
regressions cover acceleration ramps across density/unit choices, pressure
source moments at changing endpoints, held-load parity and atomic retry.

### Spatial actual-density acceleration on the original free surface

`BodyAccelerationField3` shares uniform and affine `a(x) = a0 + G(x - x0)`
loading on the native D3Q19 free surface. The affine owner stores constant-size
coefficients, with no dense spatial buffer or per-cell material callback.
Ordinary steps hold the field; endpoint advances retain old collision force
and use each cell's next actual density times its next local acceleration in
accepted velocity and gas-pressure reconstruction. Newly wet cells evaluate
the field at their own centres. Native centres are `i + 1/2`; SI centres are
`(i + 1/2) dx` in the original zero-origin Cartesian grid. SI conversion maps
`a0` by `dt²/dx`, `G` by `dt²`, and `x0` by `1/dx` once per field binding.
Piston relocation holds that world-coordinate field fixed and zero-time
reframes transported populations at their destinations before the physical
momentum audit. Uniform fields retain their existing population arithmetic.
The shared owner excludes independent force and rigid-frame ownership;
candidate refusal or cancellation retains the prior field and trajectory.
This prescribed position-dependent load does not supply Coriolis coupling or
thermodynamic density expansion. Eight regression definitions cover independent
density/momentum endpoints, local pressure moments, new interfaces, relocation,
uniform parity, ownership/rollback and SI scaling. Scoped formatting and
whitespace checks passed; native compilation and execution remain outstanding.

D3Q19 records direct donor/recipient mass routes, including reversed routes for
conversion debt. The original extensive enthalpy transport consumes these same
routes, avoiding temperature mixing through a shared zero-storage pool. The
explicit historical `Local` policy retains excess or debt without an adjacent
recipient as unregioned carry; without spatial provenance it is not injected
into another component on a later step. That mass-only policy does not preserve
conversion momentum or resolve orphan geometry. The default conservative policy
refuses unlocated carry or unavailable adjacent capacity before publication.
Neither policy supplies a conservative capillary surface-energy discretization.

Focused G0 tests cover both conversion signs, disconnected-drop locality,
per-node mass-route balances, orphan carry retention and scalar D3Q19 tension
admission. The G1 enthalpy regression compares two independent conversion
components and reverse-debt cycles against analytic local energy balances.

### Bounded local conversion carries physical liquid momentum

`ConversionRedistribution3::ConservativeLocal` is the native constructor and
timestep default, also selectable on `MaterialFreeSurface3` or the already bound
`MaterialFreeSurfaceEnthalpy3`. Explicit Local and LegacyGlobal conversion
rules retain their arithmetic. The conservative mode completes signed
over/undershoots at zero conversion tolerance and
routes each excess/debt only to adjacent interfaces with actual available
capacity/mass. Dry ghosts cannot pay debt, full cells cannot overfill, and no
unlocated carry or distant component supplies a missing recipient.

The converted cell's physical pre-conversion velocity accompanies every signed
represented mass delta. This includes negative old `m*u`: reversing a debt's
scalar route does not substitute the receiving cell's velocity. Contributions
are accumulated relative to each recipient's old velocity, then its population
velocity is updated once in the actual Guo frame. Canonical density and pressure
remain exact, and physical nonequilibrium stress is retained to roundoff.
The SAME actual routes feed the original extensive enthalpy transport.

`conversion_momentum_report()` records mass and component momentum defects,
redistributed mass, route count and signed numerical kinetic change. Remap
kinetic change is never physical heating. Local representation defects must fit
population-roundoff bounds; insufficient adjacent storage, inadmissible Mach,
nonpositive populations or unrepresentable residual momentum refuse. Public
normal and prescribed-pressure steps publish the complete candidate atomically.

This closes conversion's represented linear momentum, not the full streaming
integrator's momentum/energy balance, angular momentum, capillary surface energy
or phase-density volume sources. Native regressions cover unequal densities, Guo
forces, advancing/receding and disconnected fronts, actual before/routes/after
mass continuity, a real collision/stream/conversion step, zero-pressure replay,
co-moving identity and late public-step rollback. These regressions use the
actual constructor defaults; binary ingress separately checks advancing and
receding momentum under force, and both pressure entries retain conservative
conversion and late-refusal atomicity without a policy override.

`ConversionRedistribution3::ConservativeLocalAngular` extends that native
conversion with angular closure. It preserves signed tracked mass, Guo physical
linear momentum and cell-centre angular momentum separately for each component
of actual adjacent conversion routes. Original source inventories are captured
before class changes, including negative mass debt and retained source fluid
that later receives a transfer. Detached fronts cannot exchange corrective
impulse or torque through a global pool.

The final positive masses define a component centroid and inertia. A
mass-weighted least-change translation and rotation close the two momenta;
the same scaled SVD owner already used for tetrahedral liquid ingress admits
only represented inertia directions. The actual native population readback
must then satisfy the existing local roundoff policy. Unresolved spin,
inadmissible Mach, nonpositive populations or failed closure refuse before a
public normal/pressure step publishes its candidate. Original transport and
readback defects remain accounted for alongside final component defects.

The conversion receipt supplies `absolute_angular_momentum_defect = Some(...)`
only for this mode, about fixed local component origins; other modes return
`None`. Signed numerical kinetic change includes the corrective motion and
never becomes physical heat. The remap is in the declared wall-bounded lattice
frame: relocating cell mass changes its first moment, so a uniform velocity
boost can require a different angular correction. This is neither a geometric
projection that preserves mass centroids nor a claim of whole-timestep angular
or energy conservation. Existing modes retain their arithmetic.

Focused regression sources independently check asymmetric forced fronts,
detached components, a retained source paying signed debt, unrepresentable spin
with one surviving cell, density preservation, numerical kinetic change and
public normal/pressure rollback. Their presence is not native execution or
experimental validation evidence.

### Conservative conversion on the resolved thermocapillary thermal trajectory

`MaterialFreeSurfaceEnthalpy3::with_conversion_redistribution` selects the native
mode after attaching planar, height-graph or local sourced thermocapillarity.
The original thermal step still advances one flow/enthalpy clock, reconstructs
surface traction and pressure, and re-admits the same material sources.

Resolved-fill completion checks actual interface mass against its canonical
population density. An already bounded ConservativeLocal endpoint retains its
populations, exact native routes and signed conversion-energy receipt. An
unresolved private candidate delegates to the same momentum-carrying native
conversion, appends its actual routes, then checks the resulting bounds. It
does not apply a second mass-only projection or overwrite the step's receipt
with an empty conversion. The original enthalpy solve consumes the complete
before/routes/after account; numerical kinetic change never becomes heat.

Occupied surface volume remains `mass/density`, so integrated Marangoni force
is the patch traction even with unequal numerical densities. Native regressions
exercise that balance, three real thermocapillary timesteps, preserved geometry
receipts, and unresolved-candidate momentum/route closure. Full material tests
cover the original source/caloric owner and transactional thermal refusal.
Existing slope, cut-volume and source/phase-density gates remain mandatory.
Vanishing cut films still refuse; this does not supply thermal expansion,
phase-density volume flux, a liquid EOS or pressure-work closure.

### Conservative D2Q9 interface momentum

`FreeSurface::step_with_momentum` advances the original D2Q9 collision and
streaming operators, then transports each signed conversion mass with its
donor's physical velocity. Excess enters only adjacent interface capacity;
negative debt draws only from adjacent positive inventories. Periodic aliases
are deduplicated and disconnected droplets cannot exchange through a global
pool. An unlocated carry or insufficient local mass/capacity refuses.

All source transfers accumulate before one velocity update per recipient.
The original Guo source-frame transformation and an equilibrium difference
retain physical nonequilibrium stress; the native density reduction remains
bit-identical. No equilibrium reset erases viscous stress. A zero-mass recipient
with unresolved momentum refuses. New interface populations are admitted before
they become initialization donors, including finite forces, relaxation and
the existing low-Mach domain.

One staged surface makes cancellation and admission failure atomic. The returned
`ConversionMomentumReport2` records signed conversion kinetic change and local
mass/linear-momentum defects. Numerical conversion loss is not physical heat.
The receipt excludes collision/streaming, pressure traction and imposed forces.
The linear-only policy does not preserve first mass moments or angular momentum;
the opt-in angular policy adds local angular closure without preserving first
mass moments. Native `step` and the cell/link gas-pressure convenience entries
use this same conservative conversion owner.

Seven focused G0/G1/G3 regressions cover unequal-density advancing and receding
fronts, signed debt, physical boosts and forcing, independent droplets, periodic
capacity aliases, retained physical stress, native zero-conversion replay,
cancellation and newly wet-cell refusal.

### Native and material D2Q9 defaults carry conversion momentum

`FreeSurface::step`, `MaterialFreeSurface2::step` and
`MaterialFreeSurfaceWetting2::step_with_gate` delegate to their existing
`step_with_momentum` owners. Ordinary native and material
stepping therefore carries signed donor momentum through actual adjacent
conversion capacity and publishes one complete candidate. Material stepping
also admits the existing SI conversion account. Dynamic wetting publishes the
same candidate together with its accepted source queries, contact geometry and
angle history.

Native `step` replaces its former mass-only conversion and reuses one staged
surface. Its convenience signature and the material `step` signatures remain;
conservative refusal panics without publishing the proposal. Explicit momentum
entries provide typed cancellation and refinement refusals. Native cell/link
pressure convenience entries similarly discard only the conversion receipt,
retaining physical gas traction, the selected angular policy, and their existing
error type. No extra clone is added to the already staged pressure path.
The wetting convenience entry retains its original error type and canonical
cancellation classification.

`resolved_linear_momentum_kg_per_s` sums actual liquid mass times Guo-corrected
velocity and converts to SI momentum per unit out-of-plane depth (kg/s). Gas,
walls and unlocated carry have no represented liquid momentum. External forces
and wall/pressure traction can change the observation; it does not certify
whole-step or angular-momentum conservation.

Native regressions independently measure signed liquid momentum before and
after advancing/receding conversions, including body forcing and real pressure
force/torque. They exercise ordinary and both pressure defaults, selected angular
closure and complete rollback when local capacity is unavailable. Material
regressions compare default and explicit conservative evolution through actual
nonzero conversion routes, check SI momentum against physical
cell observations, and preserve wetting-source refusal and cancellation
atomicity. Native capillary fixtures cover all four wall orientations and the
existing material parameters; the gravity fixture starts at physical rest in
the held-force frame before evolving real drop conversions.

### Density-consistent physical state at free-surface ingress and conversion

`FreeSurface3::from_cells` closes the canonical population sum to each supplied
density before applying the Guo force frame. Tracked mass remains the authored
`fill*density`; a full interface therefore has no rounding-only excess to
convert. The existing bounded density-closure owner is shared with strain
initialization. Forced/unforced ingress regressions check exact density and
mass, physical momentum, and absence of spurious conversion across density scales.

Newly wet D2Q9 cells and default D3Q19 cells average their neighbors' physical
velocities, including each neighbor's own `F/(2 rho)` correction. They then
initialize `f = feq(rho,u) - S(u,F)/2` in the receiving cell's force frame.
D2Q9 uses that cell's `rho*g + fext`; D3Q19 retains its uniform reference-density
force or the existing averaged local-force field. Averaging raw velocities and
then dividing a force by averaged density would create a false velocity jump
at unequal density. An unshifted equilibrium instead adds an extra half-force
impulse; shifting only the equilibrium velocity leaves a false `F F/(4 rho)`
physical stress. The same full source transform now initializes explicit D3Q19
`from_cells` ingress, preserving its supplied fill, density and velocity.

The existing mass redistribution and source-force owners remain unchanged.
Focused G1/G3 tests independently check physical velocity and the full second
moment `rho*(cs^2 I + u u)` at nonunit and unequal densities, with uniform and
spatial forces, and close the original mass ledger after conversion. This
removes force-frame artifacts; it does not claim a globally momentum-conserving
VOF remap, transferred nonequilibrium stress or a physical liquid equation of
state. Frozen version-1 replay continues to expose its historical behavior.

### Physical momentum at dense force changes and rest initialization

The D2Q9 `Grid::rebind_gravity` and `Grid::rebind_forcing` operations change
piecewise-held acceleration and local force density in the original Guo chart.
For every active fluid/interface cell they use its own density and physical
velocity and replace `f` by `f + (S(u,F_old)-S(u,F_new))/2`, where
`F = rho*g + fext`. Physical mass, momentum and the force-corrected second
moment are preserved to population roundoff, including evolved nonequilibrium
stress. No time elapses and no physical impulse, work or heat is added. The
next original collision/stream step supplies the held force's actual impulse.
Direct writes to the public force fields retain the previous stored populations
and therefore have a different half-source/time convention; they are not the
physical-state-preserving operation.

Admission checks the full field lengths and finite forces, plus each active
cell's density, physical velocity and proposed populations before publication.
The caller's reusable collision buffer holds the proposal. Gas/wall population
bits, topology, relaxation and boundary settings are retained. Identical inputs
leave both accepted populations and scratch unchanged. No second force-field
copy or additional persistent buffer is allocated.

Ordinary dense channel, rheology, two-level refinement, dam-break and material
free-surface constructors now initialize physical rest in the declared source
frame. Previously they encoded zero raw population momentum and hence started
with an unintended half-force velocity; the refined levels even started with
different physical velocities. Focused regressions check full physical moments
at nonunit densities, changing spatial loads, atomic late-cell refusal, all
five D2Q9 collision families through switch/reversal/off histories, periodic
whole-grid impulse closure, and sourced SI gravity independent of liquid density.
This does not change the existing force, pressure, boundary or free-surface
transport laws, or establish temporal order for arbitrary force histories.

### Material-owned dense force histories

`FreeSurface::{rebind_gravity,rebind_forcing}` reuse the native collision scratch
to change held loads without stepping or altering tracked mass, topology or
conversion history. `MaterialFreeSurface2` exposes the same operations with SI
inputs: gravity in m/s² and a complete row-major external field in N/m³. The
original source density and dx/dt convert `g_lat = g_SI*dt²/dx` and
`F_lat = (F_SI/rho0)*dt²/dx`; collision uses `rho_lat*g_lat + F_lat`. Gravity-only
replacement retains the external load. Physical density, momentum, stress and
fill survive to native population roundoff, with the material source unchanged.

All conversion precedes atomic native publication. Nonfinite inputs,
disappearing nonzero loads and malformed fields refuse without changing accepted
state. One temporary converted field is used, with no persistent duplicate.
Three authored integration regressions cover SI force histories under density
and scaling substitution, pressure-evolved nonunit density and stress, retained
external loads and atomic refusal. These integration tests remain unexecuted
under the dependency-build resource limit; focused native force-frame tests are
separate evidence and do not constitute a full material/free-surface gate.

### Dense D3Q19 plate force histories

`PlatesGrid3::rebind_forcing` changes a held volume-force field through the
shared Guo half-source transform. The original density, physical momentum and
complete force-corrected stress survive to population roundoff. Wall
populations and relaxation times are retained; the next native collision
supplies the actual new impulse. Existing private collision storage stages the
proposal with no allocation or persistent duplicate. Invalid layouts, nonfinite
forces and unrepresentable live cells refuse before accepted publication;
identical input preserves population and scratch bits.

`plate_channel_flow3` and the plate-grid density perturbation initializer now
start at physical rest under held forcing. Uniform channel initialization uses
a private selector without allocating another force field. The raw public
`fext` assignment convention and the separate legacy Duct/Boundary constructors
are not redefined. Four focused native tests cover nonuniform density and all
stress components, switched/reversed/removed force through original collision
and streaming, late-cell refusal/retry, and actual Newtonian/Carreau startup.
These selected-production-source checks do not replace the full Cargo gate.

### Finite-slope thermocapillary height graphs

`free_surface3_from_cells_with_height_graph_thermocapillarity` binds the original
material sigma(T), caloric state and native VOF/Guo trajectory to a complete,
single-valued tank interface. `HeightGraphThermocapillarity` declares its column
axis, liquid side, maximum slope (at most one) and minimum occupied interface
volume. A graph may cross lattice layers and have multiple partial cells per
column. Both orientations and all three Cartesian column axes are supported.

Each monotone liquid-to-gas column supplies a liquid-volume height and
interface-band-weighted sigma. Centered interior and second-order edge
differences reconstruct both height and sigma derivatives. The full graph
metric projects the gradient tangentially, including its axial component;
actual patch area sqrt(1+|grad h|²) and occupied interface volume set force
density. This avoids a small-slope approximation and does not double-count
normal capillary pressure. The same column sigma extends onto interface
pressure cells, including empty ghosts; bulk caloric temperatures remain
unchanged. A dense column buffer and two bounded cell passes avoid a global
surface solve or per-cell allocation.

`with_equilibrium_planar_wetting` also admits one planar SIDE wall on this
graph. The original interface W(T) or equilibrium-angle claim and original
liquid sigma(T) use the same actual local caloric temperature. Every positive
occupied interface cell in a wall column must give exactly one angle; separate
columns may differ. Missing angle samples, through-band angle differences and
a graph-normal wall refuse. The two end columns touching another side wall
must remain exactly neutral (pi/2); nonneutral two-wall corners are unresolved.
At an XMin wall with liquid below a z graph, the imposed slope is
`h_x = -cot(theta)*sqrt(1+h_y^2)`, with the corresponding signs for the opposite
wall, liquid side and axis rotations. This slope supplies the SAME graph
normal, actual area and tangential gradient projection. Its total slope still
obeys the declared finite-slope envelope. Native normal pressure keeps its
existing smoothed-VOF angle-ghost curvature and original column sigma; this is
the same continuum angle, not a shared discrete variational surface-energy law.

Initial attachment, each heated flow preparation and the advected endpoint
install the current angle BEFORE reconstructing graph traction. The original
single candidate lattice and native force rebind retain physical momentum;
the next collision supplies the actual impulse. Storage-free wet ghosts use
the same enthalpy extension at attachment and heated endpoints. A separately
solved interface-temperature stress boundary is refused with graph wetting,
because it has no matching solid/liquid temperature closure for W(T). Source,
geometry or endpoint refusal rolls back populations, heat, force, pressure
coefficients and the original wetting receipts together.

Source/geometry/slope/carry/cancellation refusal discards the entire candidate.
Raw interface fractions must remain monotone and in [0,1]. Before endpoint
geometry, the material candidate resolves native conversion hysteresis through
actual adjacent-interface mass routes. The original enthalpy transport solve
uses these same routes; geometric clipping supplies neither mass nor heat.
This remains a resolved height-graph reduction with held-force time splitting;
it does not admit overhangs, detached drops, internal solid walls, dynamic
contact lines, surfactants or a conservative surface-energy closure. The curved
solved-temperature heat boundary uses the same admitted graph geometry.

Native manufactured tests independently check graph-metric projection, area
and volume scaling, quadratic mixed slopes, constant sigma, all orientations,
pressure/traction consistency and geometry refusal. Material regressions check
SI traction, conservative mass/enthalpy evolution, exact uniform-sigma baseline
population bits, source-validity refusal and cancellation rollback.
The sidewall extension adds analytical normal/area/traction checks across all
axis and wall orientations, an independent affine VOF-pressure check, and
thermal source, same-interval impulse and atomic corner/band/slope regressions.
The six direct native curvature-source tests pass; the broader native graph
and thermal integration regressions are authored but unexecuted under the
dependency-build resource limit. These checks are numerical, not experimental
wetting validation.

### Local thermocapillary traction on detached interfaces

`free_surface3_from_cells_with_local_thermocapillarity` retains the original
material, enthalpy, native pressure and Guo-force owners while removing the
global tank-column restriction. `LocalThermocapillarity` bounds the smallest
positive cut-cell fill and the conditioning of the local source reconstruction.
Complete monotone raw-VOF height stencils share their graph selection and
derivatives with native curvature. Both try bounded axial depths of 7, 11 and
15 cells in descending raw-fill direction order, including alternate axes
when a preferred graph is incomplete. Unallocated in-domain nonsolid gas is
zero fill, so gas-tile activation alone cannot switch resolved height geometry.
The selected graph must have a normal-axis magnitude of at least one half.
Native pressure retains its smoothed-fill contact-angle fallback at walls and
where no complete bounded height graph is available.
This changes the native free-surface bit semantics to version 8. Direct
regressions require spherical Young-Laplace pressure convergence and bitwise
pressure invariance when empty gas tiles are retired without moving liquid.
Cartesian direction magnitudes that differ only within unit-fill arithmetic
uncertainty use a deterministic coordinate tie rule. This prevents a force-frame
rebind from switching equivalent height columns and rotating the local normal;
it changes neither mass nor the surface/source admission envelope.
Force rebinding preserves the canonical population density exactly, preventing
full cells from becoming phantom partial owners through a reduction ulp. A
bounded rest correction and verified final-add fallback retain physical
velocity and force-corrected stress to population roundoff without mass routes.

A volume-matching plane/cube intersection supplies each partial cell's surface
area and centroid. A weighted affine fit uses only face-connected wet samples
in a 3-by-3-by-3 neighborhood. Its tangential projection times actual area,
divided by occupied liquid volume, supplies traction. The same fit extrapolates
sigma to the patch centroid for native Laplace pressure. Empty/full closure
ghosts receive coefficients from connected partial patches and own no duplicate
surface area or force. The reconstruction uses fixed-size per-cell arrays and
bounded plane inversion; it performs no global surface solve.

Source fitting, cut geometry, tiny volumes, pressure-domain violations, missing
support or cancellation refuse before publishing a candidate field. Walls,
contact lines, ambiguous disconnected ghost neighborhoods and all-full/all-empty
unresolved interfaces are not admitted. Ordinary evolution reconstructs sigma
from sourced bulk cell values and requires spatial refinement. A complete
solved-temperature boundary can instead sample source sigma at each actual
PLIC centroid, as described below. This extension does not add
surfactant transport, capillary heating or a conservative surface-energy step.

### Curved solved-interface heat exchange

With `enthalpy-coupling`, `height_graph_interface_patches` returns the current
physical centroid, outward normal, actual area and liquid/gas owner of every
admitted graph column. These patches use the same volume-height reconstruction
as native traction. `FreeSurfaceHeightGraphPatch::with_exchange` attaches
explicit positive liquid-normal and gas-film distances, emissivity and the
radiative surroundings. The surface quadrature point is not presented as an
occupied-liquid centroid; a column spanning multiple cut cells declares its
outermost positive cell as the thermal storage owner.

Surface exchange admits exactly one current patch per column, with unique
owners and the reconstructed area, including for constant-sigma source claims.
The thermal area budget replaces the projected outer face with its actual
curved area while preserving the independent Cartesian contact-face limits.
The existing bounded zero-storage surface solve supplies liquid, gas and
radiative heat transfer. Original sigma(T_surface) drives both native pressure
and tangential stress; bulk material temperatures remain their caloric values.

After native advection, the adapter reconstructs live graph geometry and
re-solves interface temperatures using the advected liquid states and retained
gas boundary. This endpoint refresh supplies no second thermal interval and no
extra heat. It retains the graph's strict raw-fill and slope envelope, explicit
normal distances and held-force time split. Detached local interfaces use their
own PLIC patch geometry and connected surface gradient, described below.

Focused tests check all six orientations, a curved cap crossing lattice layers,
the corrected exposed-area budget with five internal contacts, shared source
sigma/pressure/traction, mass and heat conservation, and transactional malformed
patch refusal for both variable and constant sigma. Existing planar solved-T
and finite-gas atomicity regressions remain unchanged and pass.

### Conservative completion of resolved interface fills

Material thermocapillary candidates call `resolve_interface_fills` after native
streaming/conversion and before enthalpy advection. Excess or missing interface
mass is shared across adjacent interface capacities using bounded D3Q19 local
storage. If that first ring is insufficient, the original conversion cascade
exposes another layer at an exact threshold. Every transfer enters the original
mass trace and consumes its route budget; unresolved carry or lack of progress
refuses the private candidate. The initial mass snapshot is retained, and the
updated endpoint snapshot and complete trace drive the original implicit
enthalpy transport solve.

This completes physical fills rather than clipping their geometric appearance.
The ordinary native trajectory retains its original conversion hysteresis and
LegacyGlobal behavior. Resolved thermocapillary trajectories still enforce
their slope, stencil, minimum-volume and material-source envelopes. Conversion
does not gain a new momentum-conservation or capillary-energy claim.

Regressions exercise overfill and underfill, disconnected-drop locality,
interface-front expansion, bounded work and cancellation, and advancing curved
graphs of both curvature signs crossing lattice layers. A separate regression
passes the amended trace through the actual enthalpy transport solver and
compares extensive heat and local temperatures with an analytic two-cell case.

### Solved interface temperature on detached local surfaces

`local_interface_patches` returns each actual partial-cell PLIC patch's centroid,
normal, area and occupied-liquid boundary area. `FreeSurfaceLocalPatch::with_exchange`
binds explicit liquid-normal and gas-film distances, emissivity and surroundings.
The interface divides liquid and gas subvolumes of the same partial voxel;
equal liquid/gas coordinates are admitted only by this resolved Local entry.
Empty/full closure ghosts own no patch and receive no duplicate thermal supply.
The thermal area budget uses the actual cut surface plus clipped wet Cartesian
faces, while independent Cartesian contact limits remain unchanged.

The original zero-storage surface solver supplies one interface temperature for
liquid conduction, gas-film exchange and radiation. Original source sigma at
that temperature supplies both native capillary pressure and tangential stress.
A weighted tangent-plane fit uses connected PLIC centroids, excludes opposing
surface normals and admits its conditioning explicitly. It does not reuse bulk
temperature gradients or extrapolate solved surface sigma a second time.

`FreeSurfaceInterfaceSolve::local_boundary_evolution` defaults to `FixedOwners`:
each reconstructed owner must retain a declared boundary. `UniformReservoir`
explicitly permits newly exposed patches when every contact has identical
source binding, complete context, normal distances and radiative data. The
number of live patches remains bounded by the declared thermal face budget.
This prescribed-bath mode refuses finite gas; it does not invent a remapping
of finite gas storage or spatially varying boundary data.

After native motion, live PLIC geometry and advected material states re-solve
the endpoint interface temperatures for pressure/traction only. Heat remains
accounted once by the original thermal interval, and bulk material temperatures
retain their caloric values. Endpoint ownership and geometry validation also
apply to constant-sigma sources; those sources do not acquire a persistent
temperature-dependent force-boundary requirement. These closures retain the
existing held-force time split and do not
add surfactant transport, surface entropy, contact-line laws or capillary heat.

### Implicit source-resolved surface exchange for stiff thermal capacities

With `implicit-enthalpy-coupling`, `FreeSurfaceInterfaceSolve::implicit` selects
the original thermal owner's backward-Euler residual. The default remains the
existing explicit trajectory. The optional path supports prescribed solved
surface reservoirs or one finite isochoric atmosphere, together with declared
liquid-liquid, liquid-solid and solid-solid conduction links. Original finite
solid calorics, ordered contact sources, direct solid/gas films and solid
radiation participate in that same residual. Electrical coupling and separate
liquid gas/radiation APIs still require additional coupled residuals and refuse.

Liquid and solid coordinates are specific-enthalpy increments in J/kg, mapped
to actual cell energy increments with the original masses. A finite atmosphere appends
one original specific-internal-energy increment to this same nonlinear system. A small cut cell therefore does not receive
the same finite-difference joule probe as a large cell. Frozen positive row
scales combine cell mass and the interval's conductance/caloric slope bound.
The existing bounded Newton-Krylov and fallible directional-difference owners
solve these scaled coordinates. The derivative is numerical, not an analytic
material tangent; unsupported probes retain their original source refusals.

Every trial privately re-evaluates surface temperature, source conductivity,
gas-film properties, radiation and liquid-link conductances at its endpoint.
Trial evaluations never integrate heat. Only the accepted endpoint payload
supplies heat for that interval through the original conduction observer, once.
The resulting interface sigma continues into the same native pressure/traction
refresh and conservative mass/enthalpy transport.

Acceptance independently checks every unscaled cell energy residual and the
total energy residual in joules. The whole-call storage change is also checked
against accumulated external supply. Unknown count, Krylov dimension, Newton
attempts, residual evaluations and interval count are explicitly bounded;
cancellation or any refusal publishes no accepted fluid/material state.
`maximum_dt_s` still limits each interval for accuracy, and `maximum_substeps`
still bounds total intervals. This removes the explicit thermal stability limit
without enlarging capacity, clipping geometry or inferring heat. Backward Euler
and the thermal/flow split remain first-order; temporal refinement is required.

### Inelastic freezing into an explicit stationary support

With `enthalpy-coupling`, the exact-phase feedback owner's
`advance_with_supported_freezing` permits a bounded moving-to-solid endpoint.
It requires native unsteady physical heating, equal constant phase density,
fixed geometry, a no-op inner observer, and an explicit infinitely massive
zero-velocity `StationaryPhaseSupport`. Only the requested final endpoint
defers moving-face admission; earlier substep phase events still refuse.

Capture uses the SAME staggered face mass as native kinetic accounting.
Removed linear/angular momentum reaches the support, whose zero velocity does
no mechanical work. Removed kinetic energy heats the original mapped caloric
cells through `apply_heat_j`, with native half-cell weights. The separate
`support_capture_heat_j` channel is part of generated physical heat; molecular,
Darcy and backward-Euler numerical losses remain distinct. Final kinetic
history includes capture and continues on the same accepted clock.

Surviving native velocities remain unchanged. A nonzero face between newly
solid and surviving liquid refuses because its momentum redistribution needs
a pressure/impulse solve. After depositing capture heat, all material sources
are queried again and the final exact phase mask must agree. Remelting,
unsupported capture, source exits and cancellation publish none of the thermal,
flow or support candidate. This first-order endpoint closure admits complete
moving components and zero-flux phase boundaries; it does not supply moving
solid mechanics, shrinkage or arbitrary partial freezing.

The `liquid_flow_enthalpy::mushy::occupancy::locking` regressions and native
`phase_capture` unit test specify independent energy/torque references,
graded half-cell transfers, data-only viscosity/latent changes, freeze/remelt
continuation, preserved separate flow, and deterministic refused retries.

## Complete melting into the original moving liquid owner

`MaterialFlowFeedback<ResolvedPhaseTransportState, _, (), PhaseFlowTransport>
::release_fully_liquid` transfers an accepted, completely molten fixed-grid body
into the ORIGINAL `MaterialFreeSurfaceEnthalpy3`. It takes no independent
temperature or velocity seed. The original registered caloric object, selected
transport claims, physical query context, whole-cell mass, specific enthalpy
and accepted native face motion determine the receiver. Release returns that
native owner plus a transfer report, so every existing finite-wall, finite-gas,
electrical, magnetic, surface and thermal method remains directly accessible.
No parallel liquid integrator or restricted forwarding wrapper is introduced.

This boundary requires the same uniform cubic grid and original fixed density,
every mapped material cell exactly liquid, original retained unsteady momentum
and kinetic audit, and no extra finite thermal bodies or observer history.
`StationaryConstraintRelease::Withdraw` explicitly retires stationary support,
its old boundary schedule and thermal contacts at the accepted time. Receiving
permanent walls and gravity are declared for the new domain. Old source
rheology, phase-conductivity and buoyancy extensions require an explicit separate
handoff and refuse here. Surface tension must be explicitly neglected because
this rung has no owned surface-creation-energy account. Partial melting,
moving-mould work and residual solid stress are not silently discarded.

MAC face-to-cell projection on this shared grid preserves linear and angular
momentum. Source momentum is independently integrated over the original
half-cell face masses; source kinetic energy is checked against the original
native audit. Target quantities use ACTUAL force-corrected initialized native
samples and mass ledgers, not an authored velocity claim. The gate checks mass,
total caloric enthalpy, both momenta and a separately declared maximum kinetic
projection deficit. This representation loss is reported and bounded; it never
becomes physical heat or external work. Initial numerical pressure is zero,
not an equilibrated free surface or a pressure-projection solve. Grid/time
refinement and an appropriate representation allowance remain necessary.

The physical clock lives in `MaterialFreeSurfaceEnthalpy3` itself.
`with_time_origin` may assign a finite nonnegative origin only before the
thermal/native owner advances. `accepted_steps` and `physical_time_s` follow
every original stepping API. The existing sole native acceptance point checks
integer overflow and representable forward time before any publication, then
increments exactly once with accepted populations, enthalpy, material state and
wetting after the final cancellation check. Thermal subcycles, failed source
queries, force refusals and rejected candidates never advance this clock.
Gas, finite-body, current/Lorentz and heat paths use that same acceptance point.

New regressions exercise actual latent heating followed by gravity-driven native
VOF, a nonzero circulating MAC field with independent linear/angular-momentum
and kinetic-projection formulas, viscosity-only material substitution in
continued motion, molecular heat without double-counting projection loss,
original source/caloric/geometry/phase/clock refusal and cancellation. The
existing combined gas/wall/current/Lorentz late-refusal and retry regression now
also checks the native clock. A nonrepresentable time origin refuses with
unchanged populations, enthalpy and count. Native Rust execution is separate
from source review and formatting validation.

## Source rheology through phase transport and complete-melt release

The phase-first material resolver can own an explicitly typed inverse-time shear
coordinate. It samples the original pinned viscosity and conductivity sources at
the native MAC rate; solid cells do not request liquid properties, and
rate-dependent density or caloric storage still requires a separate owner.

`MaterialFlowFeedback::release_fully_liquid_with_rheology` retains that source
through the existing complete-melt conservation path. Each accepted MAC cell's
enthalpy-derived temperature and actual shear rate must replay the target's
original material identity, selected claims, context and values. The receiver is
the original `MaterialFreeSurfaceEnthalpy3`, initialized by its native rheology
constructor and carrying the accepted physical clock. Its existing heat,
conduction, gas, finite-body, electrical and magnetic entry points remain usable.
A different reference shear coordinate cannot freeze the receiving response.

The shared transfer still checks mass, enthalpy, linear and angular momentum,
and separately bounds the kinetic projection deficit. Equilibrium population
initialization rebuilds receiving viscous stress; it does not transfer the MAC
stress tensor. This memoryless generalized-Newtonian rung needs startup and time
refinement and does not claim elastic memory, capillary energy, shrinkage or
moving-solid geometry. The original phase-conductivity multiplier and bounded
buoyancy projection can accompany the same native rheology owner.

The integrated regression melts the specimen, drives native circulation,
releases it and advances native free-surface motion with data-only rheology
variants. It independently checks molecular heat against the volume sum of
mu*gamma_dot^2*dt and face conduction against the same source k(T,gamma_dot),
alongside the original transfer budgets, clock and refusal rollback.


## Spatial finite-inertia immersed-body coupling on D3Q19

`BoundaryGrid3::step_with_partial_saturation_body` couples the three translational
and three angular velocities of `PartialSaturationBody3` to the native D3Q19
populations and then runs the existing pull stream. The supplied symmetric,
positive-definite inertia is expressed in world axes about the body's center;
off-diagonal products and asymmetric sampled fractions are retained. Fractions
are row-major `(z*ny+y)*nx+x`; cell centers are `[x+0.5,y+0.5,z+0.5]`.

At offset `r`, `Gv = V + omega cross r`. The same Noble-Torczynski
non-equilibrium reflection law as D2Q9 supplies solid impulse
`B*rho*(u-Gv)`, with `B=epsilon*(tau-1/2)/(1-epsilon+tau-1/2)`.
For TRT, `tau` is the even/shear relaxation time; the fluid collision retains
its independently configured odd rate. The solid reflection law is unchanged.
For central-moment and reduced-cumulant models, `tau=1/second_order_rate`;
the original native nonlinear collision and every higher-order rate remain
active in the fluid fraction. A finite reciprocal is required.
The operator reduces `K=sum(B*rho*G^T*G)` through a scalar weight, first
moments and the six independent rotary moments; one deterministic 6-by-6
Cholesky solve advances `(M+K/2)*delta_v=sum(B*rho*G^T*(u-Gv_n))+J_external`.
It never feeds an explicit fluid load back through a tiny inverse body mass.
The native collision runs once at the solved midpoint velocity. No grid-sized
matrix, per-step field allocation, nonlinear coupling iteration, or duplicated
fluid owner is introduced. Existing post/next buffers are reused; first-use
next storage and fixed 3-by-3/6-by-6 Cholesky factors allocate.

The receipt measures the actual population correction to the selected fluid collision,
equal-and-opposite body impulse and angular impulse about a caller-selected
origin, center torque reduced without distant-origin subtraction, and midpoint
body work. Actual rounded velocity changes independently determine body
kinetic-energy change and a six-component momentum residual. External linear
impulse acts at the center; the final three external components are a couple.
All stationary exterior/voxel wall exchanges are separately measured from the
blended populations using the existing native link receipt. Thus the immersed
body's action/reaction is not confused with wall reaction. Global angular
balance on a periodic domain additionally needs periodic transport flux.

This is a D3Q19/Guo, fixed-geometry/world-inertia velocity substep in
lattice units `dx=dt=1`. The caller owns pose drift and gyroscopic evolution and
must rotate world inertia when orientation changes. Fluid populations and their
mass/momentum remain in the immersed volume: body inertia is **additional** to
that retained fluid, not the physical body density with displaced interior
fluid removed. This rung does not claim moving geometry/remap, physical-density
closure, multiple overlapping bodies, body contact, lubrication, free surfaces,
high-Re validity, calibrated 3-D drag, or complete rigid-body trajectories.
It also makes no suspension-stress/rheology claim: the inherited
non-equilibrium reflection law can mispredict stresslets even when torque is
accurate (Najuch and Sun, *Computers & Fluids* 189, 2019,
<https://doi.org/10.1016/j.compfluid.2019.05.004>). Momentum balance does not
validate that constitutive response.

Admission requires finite positive mass, exactly symmetric numerically SPD
inertia, fractions in `[0,1]` positive only on fluid cells, finite held fluid forcing,
periodic or stationary-wall faces, and body center/covered-point speeds squared
below `0.03` at the initial, midpoint and final velocities. Unsupported models,
invalid requests, numerical overflow and cancellation return
`PartialSaturationError3` without publishing populations, body state or topology
admission. Scratch may change. Polling occurs at every 64-cell tile in fraction/
response assembly, collision, first-use scratch initialization and streaming,
plus immediately before publication; allocation itself is not interruptible.

G0/G3 coverage in `tests/partial_saturation_body3.rs` independently reconstructs
streamed linear and angular impulse for all six velocities with off-diagonal
inertia, off-center loading and arbitrary moment origin; checks midpoint work;
exercises mass/inertia down to `1e-18`; closes repeated full-domain momentum with
nonzero stationary-wall reaction; verifies common-translation invariance, proper
axis rotation and origin covariance; reproduces ordinary BGK bits at zero
fractions; and checks invalid-model/input/cancellation refusal and bit-identical
retry. These are discrete balance/conformance checks, not a continuum error or
performance certification, and do not close the full moving-boundary bead.

### Moving spatial bodies on the native fluid and rotation owners

`BoundaryGrid3::step_with_moving_partial_saturation_body` advances a body's
center, proper `So3` orientation, world linear velocity and body-frame spin
angular momentum, together with the actual native D3Q19 populations. Fixed
body inertia is an admitted `fs_time::lie::gyrostat::InertiaTensor3`; its native
principal-frame representation is cached at construction, not diagonalized on
every tick. The integrator composes a free half-tick, the existing six-velocity
implicit midpoint exchange, and a second free half-tick. Each free translation
is parallel to its own momentum; native rotation transports body momentum with
the inverse attitude increment and preserves spatial spin to roundoff. At the
midpoint, the body tensor is rotated into world axes and its rounded triangle
is mirrored for the existing exactly symmetric Cholesky admission. The native
coupling's actual representable angular-velocity increment supplies the spin
impulse; that impulse is rotated back and added to the retained body momentum.
Zero exchange does not gratuitously rebuild momentum from angular velocity.

`BodyFrameShape3` supplies membership in a fixed body frame and a positive
bounding radius about the center. Custom shapes must honor that bound and
return deterministic, unchanged geometry throughout a step.
`BodyFrameEllipsoid3` supplies checked positive semiaxes, including a sphere
constructor, with the largest semiaxis as its enclosing radius.
`BodyFrameSampler3` reuses one row-major fraction field and samples equally
spaced subcell midpoints, with 1 through 8 samples per axis. Only the midpoint
pose's fractions drive collision. The receipt also measures raw sampled solid
volumes at the initial and final poses: quadrature variation remains visible
and is never hidden by rescaling fractions, inserting mass or changing density.
Standalone sampling clips to the domain; an accepted motion step cannot clip.
Each half-step's swept enclosing ball must remain within the interior domain
halo and its enclosing box must avoid every fixed solid voxel. This deliberately
conservative admission encloses intermediate rotations too, can refuse a
non-contacting elongated shape near a wall, and supplies no collision response.

Hydrodynamic impulses and work remain measurements of the original population
correction. All force arms, including the supplied external world linear
impulse's arm, use the midpoint center. The final receipt independently measures
actual endpoint body linear and total orbital-plus-spin angular momentum about
the requested fixed origin. Stationary wall reactions remain in the existing
separate receipt; periodic fluid angular transport must be accounted for by a
caller making a global angular-balance claim. Body kinetic-energy change is
measured from its endpoints. Energy changes in the two native free rotations
are reported separately from midpoint coupling work; the symmetric free map is
not an exact energy integrator. A work residual exposes the remaining rounded
frame/velocity-conversion defect rather than asserting an exact energy identity.

The original fluid coupling now has a private final-validation callback before
its population swap. The moving owner completes the second free drift, swept
support admission, endpoint low-Mach and geometry checks, all receipts and a
final cancellation check inside that transaction. Any refusal publishes neither
fluid nor body nor topology admission; reusable scratch may change. Sampling,
buffer clearing and swept-support checks poll at most every 64 cells, with at
most 512 membership queries per cell. Native collision and streaming retain
their existing tile checkpoints. Scratch allocation itself is not interruptible.

This is one D3Q19/Guo immersed body in lattice `dx=dt=1`, with solid inertia
additional to all retained interior fluid. It supplies no displaced-fluid
physical-density closure, moving-mask/fresh-cell remap, discrete GCL, exact
cell-intersection volume, wall contact, periodic body wrapping, overlapping
bodies, lubrication, free surfaces, calibrated three-dimensional drag or
full-coupled temporal-order qualification. The broad moving-boundary bead stays
open. `tests/moving_partial_saturation_body3.rs` covers body-frame quadrature and
proper-axis rotation, actual moving ellipsoid geometry over repeated native
steps, sampled-volume variation, exact replay, independent full-domain
fluid/body linear and angular balances with nonzero wall reaction, body kinetic
accounts, zero-occupancy equivalence to the native free map and ordinary BGK,
quadrature/primitive admission, endpoint swept-support refusal, late cancellation,
fixed-obstacle refusal and deterministic retry after refusal (G0/G3/G4/G5).

### Simultaneous immersed rigid assemblies

`BoundaryGrid3::step_with_partial_saturation_bodies` advances `N > 0`
frozen rigid velocity states with one collision and one streaming pass on the
original populations. For total cell occupancy `epsilon = sum_a epsilon_a`,
the declared multi-solid extension uses
`B_a = epsilon_a (tau-1/2)/(1-epsilon+tau-1/2)` and population correction
`sum_a B_a (Omega_s,a - Omega_f)`. Each fraction must be in `[0,1]`, positive
only in fluid cells, and the sum must not exceed one. Fractions represent
disjoint occupied subvolumes; arbitrary callers own that geometry obligation.
There is no clipping or overlap renormalization. The combined weight is the
original Noble–Torczynski weight of the total occupied volume, including its
fully occupied limit. The multi-solid collision sum follows the usual PSM
decomposition; the total-occupancy weight and implicit body update are the
specific law implemented here, not a claim to reproduce every PSM variant
(see Kemmler et al., 2025, equations 18–22,
<https://doi.org/10.1177/10943420241313385>).

Holding incoming fluid populations fixed gives one independent symmetric
six-by-six midpoint block per body, with the same full world inertia,
translation/rotation cross terms and additional-inertia interpretation as the
single-body owner. Hydrodynamic interaction passes through the shared fluid;
there is no artificial dense cross-body inertia or sequential per-body stream.
Each body reports its actual population-correction impulse, center and world
torque, midpoint work and accepted momentum/kinetic-energy defects. Stationary
wall reaction is measured once from the composite outgoing populations. The
assembly additionally measures composite-population mass, linear and angular
roundoff relative to the sum of attributed corrections. A singleton delegates
to the existing owner, preserving its numerical operation order and bits.

`step_with_moving_partial_saturation_bodies` composes these blocks with each
body's existing native SO(3) free half-steps and midpoint geometry sampler.
The fixed inertia may be a full body-frame tensor. Each half-step has affine
center motion; pair admission minimizes the relative center distance on the
same interval and requires strictly separated enclosing spheres, including a
floating-point grazing margin. This conservative numerical check can reject
separated elongated shapes and is not interval-certified CCD. Existing
wall/halo support restrictions remain. All final shapes, velocities, work and
momentum accounts and cancellation checks precede fluid or body publication;
a late failure in the last body rolls back the entire assembly. Nested
per-body trajectory wall and fluid-forcing fields are zero: the assembly owns
both accounts.

There is one reusable fraction field per body and one shared fluid scratch
pair, with fixed-size `O(N)` body work and two small native Cholesky allocations
per body. Sampling/clearing and fraction traversal are `O(N * grid cells)`;
pair separation is `O(N^2)`. This is a small-assembly path, not a sparse
million-particle implementation. Input order pins reductions; permutation
invariance is to roundoff, not a bitwise claim. Cancellation is checked at
native lattice tiles, body iterations and bounded pair/fraction batches.
`tests/partial_saturation_bodies3.rs` checks an independent shared-cell formula,
full-tensor world momentum/work with nonzero wall reaction, permutation and
zero occupancy, moving-body replay and spatial balances, exact singleton
compatibility, overfill refusal, late pair collision and finalization
cancellation. No displaced-fluid closure, contact/lubrication model, remap/GCL,
full-system order, suspension rheology or continuum-force accuracy is claimed.

## Prescribed pressure on original NSF characteristic faces

`NsfBoundaryDrive` distinguishes an incoming acoustic characteristic [m/s]
from a pressure offset [Pa] relative to the original face reference gas.
The pressure variant sets the incoming characteristic from the live outgoing
wave at each original face reconstruction: `J- = J+ - 2 dp/(rho_ref c_ref)`.
Normal velocity, inflow entropy, tangential velocity, thermochemical state,
population reconstruction and conservative mass/energy flux retain their
existing owners. It is an ideal prescribed pressure within the existing local
linear characteristic domain, not a finite reservoir or stagnation pressure.

`set_boundary_drives` and `preview_boundary_interval` retain both controls and
all accepted fields on refusal or cancellation. The existing characteristic
methods delegate to the typed path without changing their reconstruction
arithmetic. Source pressure, temperature, Mach and relative-perturbation limits
still apply; prescribed pressure does not bypass material admission. The
focused pressure-boundary tests exercise all three axes, returning-wave and
material impedance dependence, original integrated work/mass/energy accounts,
typed characteristic identity and atomic refusals. They do not establish a
reflection bound or measured-instrument accuracy.

`NsfBoundaryDrive::ResistiveImpedance` adds a local passive acoustic load in
specific SI impedance [Pa s/m]: `p-p_ref = R*(u_n-u_ref,n) + p_drive`, with
outward normal velocity on either face and finite `R >= 0`. The live outgoing
characteristic and the original reference gas's `Z0=rho*c` determine both face
pressure and velocity at every predictor. The pressure reflection coefficient
is `(R-Z0)/(R+Z0)`: zero resistance is pressure release, `R=Z0` is matched,
and large resistance approaches the positive rigid-wall reflection. Scaled
weights retain the small velocity near that rigid limit without subtracting
opposite characteristics. They are prepared once per face or patch; evaluating
the original quadrature adds no allocation or constitutive query.

The same load reaches complementary faces, localized patches and reversible
intervals, preserving source/Mach/temperature bounds and atomic refusal.
Its source-free perturbation work is `R*(u_n-u_ref,n)^2 >= 0`; the original
boundary mass, momentum, pressure-work and total-energy accounts already carry
the exchange. There is no extra heat sink or separately evolved exterior.
`nsf3_open_boundaries::impedance` checks source-dependent reflection on all
axes, signed returning packets and refinement, native conservation, patches,
trial replay and refusals. These are plane-normal local-linear numerical
checks, not qualified oblique transparency or measured load data. Reactive
storage, frequency-dependent impedances and finite exterior inventories remain
outside this instantaneous resistive model. Resolved NSF-DG input files can
select `impedance PA_S_M PA` for a source/coupled patch's complementary face;
the final pressure value is the optional ideal source offset, not another
incoming wave or an absolute pressure.

`reservoir_aperture_state_at(P, T)` resolves held absolute pressure and
temperature through the grid's original material source while retaining the
attached transverse aperture's declared Cartesian ghost velocity. The warm
`set_reservoir_aperture_state` path preserves stored fluid density, velocity,
stress and total energy and rebases the original force and flux history.
Source-domain or boundary refusal restores the previous held conditions.

`preview_boundary_interval_with_reservoir` stages those conditions together
with the two original characteristic controls and optional selected patch.
No reservoir, field, material or ledger state is published on a dropped,
cancelled or refused trial. The mutable accepted ghost participates in the
original state hash. G0/G3 cases in `nsf3::open::tests` compare stored
physical state, warm-source refusal, exact replay and independent mass/energy
balances. The source remains a prescribed reservoir of the grid's original
composition, with no finite inventory, moving topology or new radiation owner.


## Original phase conductivity through moving-liquid release

The original `MaterialFreeSurfaceEnthalpy3` can retain `LiquidMassFractionConductivity`. Both ordinary and rheological complete-melt release transfer that exact law before publishing the receiving liquid and its clock. One shared helper keeps the pinned raw source, declared phase multiplier and effective conductivity separate. Liquid/liquid, finite-solid contact, gas-film and solved gas/radiation-interface half-cell resistances all consume that same effective coefficient, including thermal CFL admission. Source overflow and cancellation refuse before native publication. Focused regressions independently close each heat balance, preserve identity-law behavior and source claims, and carry the law through real melt/release and viscosity-dependent motion. This owner still admits only pure liquid at fixed caloric density; moving mushy phase, solidification, shrinkage and surface-creation energy remain outside this feature.

### Atomic uniform flow and held fluid loading

`BoundaryGrid3::initialize_uniform_equilibrium` prepares the requested physical
density and velocity using the original D3Q19 equilibrium and the current
uniform/per-cell Guo source frame. It deliberately replaces initial fluid
storage and stress; it is not an impulse or an evolved-flow remap. Solids,
body states, viscosity, forcing and boundary topology remain unchanged.
A mismatch with stationary walls or bodies describes a startup transient.

`BoundaryGrid3::rebind_uniform_force` changes only the uniform part of held
volume forcing. The original `thermal::forcing::reframe3` transformation retains
mass, physical momentum and force-corrected stress at every fluid cell. It
therefore applies no artificial half-force impulse when an external load
switches, and supplies no elapsed-time work. Existing spatial fields and their
selected scaling shape remain intact; the next step applies the new held force.
`set_force_field_scale` uses the same piecewise-held convention for a retained
spatial shape. Smooth endpoint histories use `step_with_endpoint_force_scale`
and its pooled counterpart, which own the old/new force timing explicitly.

`BoundaryGrid3::rebind_force_field_with_gate` replaces the spatial held load
through that same Guo transformation while retaining the uniform contribution.
It samples in x-fastest order, skips static solids and includes partially
occupied moving-body cells. The proposal allocates one native force field and
reuses collision scratch. A late invalid sample, unphysical population or
cancellation publishes neither field nor populations. Unchanged cell forces
retain population bits. The new field also becomes the shape selected by the
existing force-scale API; the legacy field setter and rebind entry delegate to
the same checked sampler. Focused native-source tests cover nonunit-density
stress preservation, actual switched/reversed collision impulses, ordering,
solid exclusion and exact retry after late refusal.

The uniform operations use existing collision scratch, allocate no new population or
force field, poll cancellation at each native tile and publish only after the
last check. A late unphysical cell preserves all accepted populations/forcing.
An unchanged uniform force preserves population bits. Regularized open boundaries
admit the same uniform and spatial forcing, with force-corrected velocity/stress
reconstruction at both the source and destination cells.
The initializer's caller supplies scenario Mach/physical-regime admission.
These entries add no pressure-reservoir, variable-density or continuum-order claim.

Focused native-source regression groups cover physical density/velocity
initialization under a held force, stress-preserving uniform switches with an
existing spatial field, elapsed impulse, identical-force bits, and atomic late
invalid-density/cancellation refusal. The forced-open cases in
`tests/boundary_forcing.rs` additionally check every face orientation, independent
physical momentum and stress, mixed wall-link ownership, serial/pooled
endpoint-force parity, and atomic late reconstruction overflow with exact retry.
These added cases were not executed in the constrained
editing environment. The original moment observer
is reached only after typed density admission. The full Cargo workspace remains
unavailable with missing constellation runtime siblings.

### Imposed fluid forces through immersed rigid assemblies

The single, simultaneous and moving partial-saturation owners admit the original
uniform plus per-cell Guo fluid force, alongside each solid's external impulse.
This enables body-force-driven flow with moving, rotating full-inertia bodies.
All non-solid cells retain fluid, including occupied fractions, so the complete
held fluid force acts once per cell. No saturation factor suppresses it and no
fluid force is silently attributed to a solid. Additional solid loading remains
in the original explicit impulse/couple input. Changing an evolved force field
uses `BoundaryGrid3::rebind_force_field` to retain physical momentum and stress.

Stored Guo populations are `f = g - S(u,F)/2`, with physical velocity
`rho*u = sum_q(f_q*c_q) + F/2`. The existing native source-frame transform
recovers `g`. The immersed correction is
`sum_a B_a (Omega_s,a(g) - Omega_fluid(g))`, added to the unchanged native
`collide_cell3(f, model, F)` for the selected collision model. Its first moment remains
`sum_a B_a*rho*(v_a_mid-u)`; the original six-velocity solve uses this same
physical velocity. The remaining imposed collision impulse is `F`. This is
the declared additive-force extension of the retained-fluid PSM law, not a
new displaced-fluid density closure or a full coupled temporal-order claim.
Existing BGK/TRT paths retain their numerical arithmetic, and zero occupancy
reproduces the native selected forced populations bit for bit.

Central-moment and reduced-cumulant models use the same physical-frame
coupling. Their native forced output is `C(g)+S/2`; the immersed owner recovers
the actual fluid increment as `(ordinary-g)-(g-f)`, equal to `C(g)-g` up to
the roundoff of the two half-source additions. Thus no second moment transform,
collision, source evaluation, heap allocation or population field is needed.
Because this increment has zero first moment, its replacement retains the
same six-velocity implicit body response and counts the imposed force once at
every solid fraction. It also retains the selected higher-order relaxation
and nonlinear cumulant products instead of substituting a BGK fluid increment.
Single, simultaneous, moving, dry-inertia and lubricated owners share this
collision path. This extends model compatibility; it does not establish a
high-Re continuum, stresslet, suspension-viscosity or moving-body accuracy
certificate for either moment model.

The focused unit tests in `boundary/partial_saturation/collision.rs` compare
recovered increments against a separate original unforced collision at three
viscous rates, with unequal higher-order rates, nonzero stress and third/fourth
content. Independent source and moment sums check mass, full external impulse
and body recoil at empty, partial, nearly full and full saturation. Existing
BGK/TRT increments are checked against their native owners. The original
`forced_partial_saturation3` integration target additionally covers moment-model
single/shared cells, all fill limits, repeated rotating/moving bodies, wall and
periodic angular transport, singleton replay and atomic refusal/retry.

`PartialSaturationForcing3` reports the sum of the actual configured fluid forces,
their moment about the requested world origin, and imposed-force work at the
mean of physical pre/post-collision velocities. This local mechanical work
account excludes streaming and pressure/internal-energy accounting; it is not
a complete fluid energy balance. Wall and hydrodynamic reactions remain
separate. For a moving assembly, only the assembly coupling receipt owns fluid
forcing and walls; nested per-body versions are zero. The existing atomic
publication and cancellation boundaries apply, with no new population storage
or second collision/stream pass.

The six focused cases in `tests/forced_partial_saturation3.rs` check an
independent Guo source and shared-cell population formula with nonzero stress,
native zero-occupancy bits, repeated full-tensor force/torque/work balances with
walls and arbitrary origin, moving geometry with shared cells, singleton
equality, late force overflow, cancellation and exact retry. Existing unforced
fluid/body suites remain applicable. No sedimentation calibration, rheology,
lubrication, contact, remap/GCL or continuum-force accuracy is established.


### Simultaneous normal lubrication for moving spherical assemblies

`SphereLubrication3` supplies the unresolved normal Stokes resistance between
strictly separated spheres. For actual positive midpoint gap `h`,
`R_eff = Ra Rb/(Ra+Rb)` and `h_eff = max(h, roughness_gap)`, the added
resistance is `zeta = 6 pi mu R_eff^2 (1/h_eff - 1/cutoff_gap)` below the
cutoff and zero above. `mu` is constant dynamic viscosity in lattice
mass/(length*time), not kinematic viscosity. The caller must match the fluid
viscosity and calibrate the resolved-grid cutoff. The explicit positive
roughness floor must not exceed the cutoff; it does not admit solid overlap.
The normal correction follows the leading Stokes gap asymptotic (Janoschek,
Harting and Toschi, 2013, equation 11 for equal spheres,
<https://arxiv.org/abs/1308.6482>; Nguyen and Ladd, 2002,
<https://doi.org/10.1103/PhysRevE.66.046708>). This implementation does not
claim those papers validate its partial-saturation cutoff or trajectories.

`step_with_lubricated_spheres` samples the supplied spherical geometry at the
original native midpoint pose and adds the correction inside the simultaneous
fluid/body solve. The geometry-free
`step_with_lubricated_partial_saturation_bodies` exposes the same operator;
its caller owns consistency between sphere radii and occupied fractions.
Each original six-by-six body block is factored once. Only three translational
impulse columns are needed to form the active-pair Schur matrix
`diag(1/zeta) + J A^-1 J^T/2`; native Cholesky solves all normal impulses
together, including pairs that share a body. There is no sequential pair kick,
new numerical solver, dense six-times-body-count matrix, grid-sized buffer or
additional collision/streaming pass. At most 64 bodies and 64 active pairs are
admitted, with an explicit caller pair cap. Pair search is quadratic in body
count and Schur storage is quadratic in active-pair count.

`coupling.lubrication` owns each body's central impulse, its orbital torque
about the requested origin, work from the actual returned midpoint velocity,
nonnegative sum `zeta * slip_mid^2`, and measured constitutive/work defects.
Its impulses are not attributed again to fluid populations. Native Guo forces,
solid external loads, full-tensor spin, wall exchange and moving-body balances
retain their owners; all per-body residuals include the additional lubrication
impulse/work. No thermal destination is implied for the model viscous loss.
No active pairs and a singleton retain the original uncorrected arithmetic.
Every endpoint and cancellation check completes before joint body/fluid
publication; a refusal may change scratch only.

This is normal sphere-sphere creeping-flow correction with retained interior
fluid and additional solid inertia. It adds no wall/tangential lubrication,
cavitation, contact, displaced-fluid closure or calibrated suspension
rheology. The midpoint law is passive in real arithmetic but not L-stable;
stiff under-resolved relaxation can reverse endpoint slip. Positive-gap
admission and balance diagnostics do not certify temporal resolution or
prevent every discretization artifact. The native tests retain this stiff
case explicitly, verify unequal-radius asymptotics and both three-body chain
modes, roughness/cutoff behavior, rotation/permutation/boost covariance,
forced moving full-tensor world momentum/work and atomic failure/retry.

### Simultaneous sphere/plane squeeze films and prescribed surface work

`BoundaryGrid3::step_with_lubricated_partial_saturation_planes` and
`step_with_lubricated_spheres_and_planes` extend the original body/fluid
response with prescribed infinite planar surfaces. `LubricationPlane3` carries
a world point, an admitted unit normal into the permitted half-space, a constant
low-Mach translation, and independently specified wall cutoff and roughness
scales. The moving entry point interprets the point at the tick start, samples
the plane at the same midpoint as the original spherical bodies, and admits
both complete affine relative half-drifts by checking their endpoint gaps.
`translated` provides the next window's prescribed geometry. No body pose or
velocity is projected to repair contact.

For a sphere of radius R and actual positive plane gap h, the additional
normal resistance is `zeta = 6*pi*mu*R^2*(1/max(h,h_rough)-1/h_cut)` below the
cutoff and zero above it. This is the leading Stokes squeeze term, not the
complete Brenner resistance or its logarithmic corrections. Dynamic lattice
viscosity mu comes from `SphereLubrication3`; each wall has independent cutoff
and roughness inputs because wall and sphere-pair resolution errors need not
agree. This normal leading term and cutoff/roughness construction follow the
sphere-wall asymptotics in Costa, Westerweel, Breugem and Boersma (2015),
section III, equation 33, <https://arxiv.org/pdf/1506.01880>. Their higher-order
terms and collision model are not implemented by this owner. Calibration to
the native partial-saturation grid remains the caller's responsibility.

Sphere-pair and sphere-plane rows enter the SAME original small Schur solve.
For a plane row the relative midpoint speed is `(v_body_mid-v_plane).n`; only
the body's response contributes to that row's mobility. A single factorization
of each native six-velocity response supplies the three translational columns.
The shared-body off-diagonal terms are retained across wall and particle films;
there are no sequential plane kicks, giant generalized mass matrix, extra fluid
collision pass, or population correction that counts the film twice. At most
64 bodies, 64 declared planes, and 64 TOTAL active films are admitted, with the
existing caller-selected active-film cap. Discovery is O(N^2+N*planes), small
Schur storage O(active_films^2), and per-plane reaction storage O(planes).
Empty plane lists retain the original APIs' arithmetic and results.

`SphereLubricationExchange3` now includes aggregate plane impulse, angular
impulse, prescribed-plane work and the active sphere-plane count. New wrapper
results also return `PlaneLubricationExchange3` in exactly the input plane
order, including zero reactions for inactive planes. Normal impulses pass
through each sphere center; the opposite plane reaction is applied at its
actual projected footpoint. Both angular impulses use the caller's same fixed
world origin. Plane work is work RECEIVED by the prescribed surface, so its
driving source supplies the negative. The measured identity is total body film
work plus total plane work plus nonnegative modeled viscous loss equals the
reported work residual. A driven plane can supply positive body work without
being misclassified as negative dissipation. Individual body and moving-body
momentum/work accounts include the complete film impulse once.

Planes define ONLY the unresolved-film component. They neither create solid
voxel cells nor impose or move a resolved no-slip fluid boundary. The complete
stationary-wall integration uses an existing matching halfway voxel wall;
translating-plane results describe the prescribed film port, not a complete
moving-wall fluid-displacement model. Existing stationary wall/halo admission,
retained interior fluid plus additional solid inertia, and actual sphere
sampling remain authoritative. All physical gaps must be positive even outside
the correction cutoff. Exact duplicate planes, overlap, exhausted limits,
nonfinite solves and cancellation refuse atomically, including final geometry
failure after fluid streaming has been prepared. Only private fluid/sampler
scratch may change before successful finalization.

Six new focused native groups pass alongside all 32 existing related groups.
Independent checks include the analytic forced sphere/translated-plane impulse
and source work under scene rotations and common translation; an independent
two-body velocity solve with one sphere pair and two different driven plane
films; independent wall cutoffs/roughness; exact empty-plane continuation;
24 ticks of two sampled full-tensor bodies beside an actual voxel wall with
Guo forcing and external loads; and strict overlap/combined-film/plane-count
admission, cancellation, late translated-plane crossing and deterministic
retry. The sampled trajectory's independently reconstructed maximum world
linear/angular momentum defect is `3.748246288001744e-13` lattice units, with
modeled viscous loss `1.3701980928724615e-5`. These source-selected native checks
do not establish the full Cargo/DSR gate, ISA parity or experimental calibration.

Normal-only creeping-flow resistance, fixed plane normals and constant
translation are the admitted model. Tangential/rolling film forces, plane
rotation, finite wall recoil, solid contact, cavitation, deformable surfaces,
thermal allocation and complete moving resolved-fluid boundaries remain
outside this increment. Midpoint damping retains its existing passivity
interpretation and unresolved endpoint reversal; it is not L-stable or a claim
that arbitrarily small-gap trajectories are temporally resolved. The broad
moving-boundary bead remains open.

### Reciprocal tangential and parallel-spin wall lubrication

`LubricationPlane3::with_tangential(cutoff_gap, roughness_gap)` adds an explicit
leading-logarithmic sphere-wall resistance to the same frozen and moving
immersed-body entry points. The original normal sphere-pair and plane films
remain present and simultaneous. The tangential scales are independent of both
normal cutoffs and retain `0 < roughness <= cutoff`. Admission also requires
the tangential cutoff to be at most `0.02*radius` for every supplied sphere.
This is a declared small-gap model domain, not a certified accuracy bound.
Actual gaps remain positive even below the roughness scale.

For either unit tangent t, let `s = n cross t`,
`u = (v_body - v_plane).t`, and `w = R*omega.s`. The missing resistance in
these power-conjugate coordinates is
`k * [[4,-1],[-1,4]]`, with
`k = (4*pi*mu*R/5) * ln(h_cut/max(h,h_rough))` below the cutoff, zero above.
The coefficients follow the Goldman-Cox-Brenner leading logarithms, reproduced
in Ryu, Hommel, Roberts and Frechette (2019), equations 7-10,
<https://link.aps.org/accepted/10.1103/PhysRevE.99.022802>. The constant terms
cancel in the cutoff difference; higher-order gap terms are not modeled.
The code evaluates the logarithm through the existing deterministic
`fs_math::det::ln` owner. Normal-axis spin has no leading singular resistance
in this model and receives no invented correction.

The two eigenvalues are `3*k` and `5*k`, acting on `(u+w)/sqrt(2)` and
`(u-w)/sqrt(2)`. Both tangential directions therefore add four scalar modes to
the original impulse-space Schur solve. This retains translation-rotation
reciprocity, pressure-induced moments and nonzero rolling-mode resistance;
zero material-point slip alone does not eliminate the hydrodynamic loss.
All original normal modes and all new modes participate in the same solve,
including off-diagonal terms for shared bodies and arbitrary admitted full
inertia tensors. The tangent basis is an internal coordinate choice, not an
anisotropic friction direction.

Each original native six-velocity response is still factored once. Its three
rotational impulse columns are evaluated only for bodies touched by an active
tangential block; normal-only and inactive-tangent paths retain the previous
three-column arithmetic. There is no second velocity solver, additional fluid
collision pass, or extra film population source. The existing maximum of
64 bodies, 64 planes and 64 physical films remains. Scalar solve work has its
own `SphereLubrication3::with_resistance_mode_limit` control: default 64,
explicit opt-in through 320. A normal film consumes one mode; an active
tangential block consumes four additional modes. Storage is O(modes squared)
and factorization O(modes cubed), with bounded polling and fallible reservation.

`SphereLubricationExchange3` now reports center spin impulses separately,
includes spin in each body's world angular impulse and midpoint wrench work,
and reports actual mode dimension and tangential-plane count. Every frozen
body's angular balance subtracts the film spin impulse exactly once; the
original moving rigid-body owner propagates that changed momentum through its
unchanged half-drifts. Each plane report includes the free reaction couple at
its actual footpoints. For body force F and center torque T at center c,
the wall receives force -F at footpoint p and couple
`-T - (c-p) cross F`. The resulting world angular impulses cancel, including
the moment of tangential traction. Dropping that couple would violate the
global angular balance.

The loss is the nonnegative sum of each scalar resistance times its ACTUAL
returned midpoint mode speed squared. Body film work includes force times
translation plus center torque times spin. Since these planes translate but
do not rotate, prescribed-plane work remains reaction force dotted with plane
velocity. The work residual includes both body and plane work plus modeled
loss. Spin impulses are never mislabeled as external actuator torque or
counted again in the fluid populations. All fields, mode limits, endpoint
geometry, native speed admission and cancellation complete before publication.

Six new focused native groups pass alongside all 38 previous related groups.
They check an independent physical two-coordinate translation/spin solve and
reaction couple; rolling and opposed-spin eigenmodes, including second-order
time refinement; a separately assembled 12-by-12 full-tensor physical system
against the native 16-mode mixed sphere/plane solve under scene rotation,
common translation and plane permutation; actual sampled close-wall motion;
independent tangential cutoff/roughness and exact inactive-option results;
and mode/asymptotic-domain admission, cancellation, late driven-plane crossing
and atomic deterministic retry. Rolling-mode errors at 4/8/16 steps are
`1.339430650594553e-7`, `3.348109667462851e-8`, and
`8.369982382530616e-9`. Over 24 sampled ticks, two full-tensor bodies beside an
actual voxel wall with Guo forcing, external loads and normal particle-pair
films have a maximum independently reconstructed world momentum defect of
`3.037859740881953e-13` lattice units and modeled loss
`7.877813066050252e-10`.

These are focused selected-source native checks, not full Cargo/DSR, ISA parity
or experimental calibration. The surfaces still supply missing-film ports;
they do not create or move resolved fluid boundaries. Fully flooded Newtonian
creeping flow, fixed wall normals, constant prescribed translation and the
small-gap leading asymptotics are the admitted model. General tangential
sphere-pair resistance, finite-gap complete drag, cavitation, solid contact,
rotating or recoiling walls, thermal allocation and complete moving-wall fluid
displacement remain outside it. The original midpoint damping is not L-stable
and does not establish temporal resolution of arbitrarily stiff gap dynamics.
The broad moving-boundary bead remains open.


## Surface-reservoir heat through flow and feedback adapters

The native material-flow advection, sourced-liquid, momentum-feedback,
electrical-feedback, magnetohydrodynamic-feedback and environment-feedback
observers forward the original inner observer's `reservoir_conductances`
hook. State-dependent half-cell and film resistance therefore reaches the
same thermal substep as internal conductance and advective enthalpy. Authored
reservoir temperature/area and original error mappings remain intact; a
missing prepared flow epoch refuses through the existing error path.

The original thermal owner installs these coefficients only for actual heat
and entropy evaluation and restores authored boundaries before publication.
No added heat source, duplicate bath credit, independent flux approximation
or new thermal integrator is introduced. Source-only electrical adapters
continue to compose physical observers through `CoupledObservers`.

A focused native-flow regression uses actual face mass routes, `Advected`,
the geometric `ConductionObserver` and the original thermal timeline. An
independent accepted-interval donor recursion verifies state-dependent bath
conductance, final enthalpy/temperature, bath heat minus outlet enthalpy export,
actual-bath entropy, fixed mass, authored-boundary restoration and rollback
after late cancellation. Four original-source tests pass in the bounded
native composition. The other five forwarding methods retain their existing
adjacent delegation/error patterns; this is not a full fs-lbm/workspace build.

### Reciprocal tangent and spin lubrication between unequal spheres

`SphereLubrication3::with_tangential_pairs(cutoff_gap, roughness_gap)` opts the
existing frozen and sampled-moving sphere owners into a leading-logarithmic
unequal-sphere resistance block. The previous normal-only path is retained when
this option is absent or inactive. Pair tangent/roughness lengths are independent
of the normal-film and plane scales. All actual gaps must remain positive;
roughness regularizes resistance and never grants permission to intersect.

For radii `a,b`, actual center distance `d`, unit normal `n` from A to B, tangent
`t`, and `s=n cross t`, define `u=(v_a-v_b).t`,
`z_a=a*(omega_a.s+u/d)`, `z_b=b*(omega_b.s+u/d)`, and
`k=4*pi*mu*(a*b/(a+b))/5*ln(cutoff/max(gap,roughness))` below cutoff. The
reciprocal resistance in `(z_a,z_b)` is `k*[[4,1],[1,4]]`. Its modes
`(z_a+z_b)/sqrt(2)` and `(z_a-z_b)/sqrt(2)` have positive resistances `5k` and
`3k`. Both orthogonal tangent directions enter the shared solve. Counterrotation
therefore retains a singular dissipative mode even when the two surface
velocities coincide. Rotation about the center normal has no divergent leading
term and receives no fabricated correction.

At contact this block reproduces all leading `YA/YB/YC` force/torque
coefficients of Jeffrey--Onishi with Townsend's corrected dimensionalization:
Townsend, *Physics of Fluids* 35, 127126 (2023),
<https://arxiv.org/abs/1802.08226v4>, equations 19-25 and 33-34a; coefficient
functions `G2ya`, `G2yb`, `G2yc`, `G4yc` in the author's
`Pecnut/stokesian-dynamics/find_resistance_scalars/functions_nearfield.py`,
observed source blob `bab8a0e4d281db4d6fc2e0c7965ec0220e1e7bbe`.
The normal coefficient remains the separately calibrated existing `1/gap` law.

At positive gap the rotational velocities are measured relative to the actual
line-of-centers angular velocity. This is an explicitly objective finite-gap
completion of the leading logarithmic block: replacing the contact distance
`a+b` by actual `d` in the frame rate changes omitted-order `h*ln(h)` terms.
Forces and both center torques are the power-conjugate derivatives of those
same modes. Thus every pair conserves linear and orbital-plus-spin angular
momentum and has a null common rigid motion, including common rotation at
finite gap. This is not the complete higher-order two-sphere resistance, a
far-field drag model, or an accuracy certificate. The tangent cutoff must be
at most `0.02` times every supplied sphere radius; this is a conservative
declared small-gap domain, not a proven error bound. Arbitrary positive radius
ratios are admitted only when their actual geometry and resistance remain
representable. No contact, friction, cavitation, background-strain stresslet,
normal-axis finite spin drag, or dense-suspension closure is inferred.

The original 6-velocity response factors are reused once per body. Only bodies
touched by active tangent modes need their three angular-response columns.
Unequal pair modes carry distinct angular covectors for both original bodies;
they join the same simultaneous solve as normal pairs and prescribed-plane
films. There is no sequential pair kick, extra fluid sweep, substitute inertia,
or dense global physical mass matrix. The limits remain 64 bodies, 64 planes,
64 active physical films and an explicitly admitted scalar-mode limit of
1 through 320 (default 64); an active tangent block adds four modes to any
normal mode. Dense response storage and factorization remain O(modes^2) and
O(modes^3). `active_tangential_pairs` reports the active sphere-pair tangent
blocks independently of `active_pairs` and plane counts.

Both center spin impulses flow through the existing `body_spin_impulses`,
`body_angular_impulses`, full-tensor moving-body half-drifts, and per-body balance
receipts. Work uses actual returned midpoint translation and angular velocity;
mode loss is nonnegative and the actual constitutive/work defects remain
visible. Midpoint damping is passive but not L-stable; time refinement is still
required for stiff relaxation. Existing private finalization preserves atomic
body/fluid publication on invalid geometry, exhausted modes, cancellation or a
late swept-sphere intersection.

`pair_tangential_lubrication3.rs` supplies six focused groups: independent
dimensional `YA/YB/YC` probes at equal, 1:3 and 1:10 radii in both endpoint
orders and convergence to contact coefficients; zero-slip counterrotation with
second-order time refinement; common-rigid-motion nullspace, scene covariance
and endpoint exchange; an independent full 18-velocity solve for shared unequal
pairs and translating-wall films; a 24-tick sampled unequal/full-inertia body
trajectory beside an actual voxel wall with Guo forcing, independent global
momentum and midpoint wrench/work reconstruction; and independent cutoffs,
roughness, normal-only equality, mode/domain refusal, cancellation, exact late
swept-sphere refusal and deterministic retry. The six new and 44 prior focused
native fluid cases pass. These are source-selected native tests, not the
unavailable repository DSR quality verdict.

## Shared physical acoustic scaling

`lattice3::AcousticScalingPlan::from_material` consumes both checked
`FluidAcoustics` and `FluidTransport` from the same original material. Density
must match bitwise across the two views; the existing scaling kernel then uses
physical `c` and `nu = mu/rho` to determine `dt`, relaxation and Mach limits.
`from_gas` delegates to this admission, refusing stale public gas fields.
Source-resolved liquids use their own acoustic response without ideal-gas
inference or a numerical sound-speed default.

Three focused regressions cover source-only liquid substitutions changing
physical timestep and relaxation, accepted-gas bit parity and stale-state
refusal, and inconsistent acoustic/transport density. This mapping does not
extend the original low-Mach lattice solver to shock, cavitation or multiphase
physics, nor certify its spatial or temporal discretization error.

### Shared finite atmosphere for moving liquid and finite solid walls

`FreeSurfaceThermalBody::with_gas_faces` attaches explicit stationary solid/gas
patches to the SAME finite atmosphere used by the original moving-liquid
conduction and surface-exchange entry points. `FreeSurfaceBodyGasFace` supplies
native wall/gas coordinates, patch area, solid centroid distance and continuum
gas-film thickness. No independent bath temperature, conductivity or thermal
inventory is introduced.

The original solid source and caloric chart supply k(T), finite mass and
temperature; `ThermalMaterialPoint::at` preserves their material-card identity.
The original gas model evaluates film transport at the mean solid/gas temperature
and the finite parcel's current EOS pressure and composition. Selected film
claims remain pinned across each interval. All solid and liquid rates use one
thermal substep state, and their aggregate opposite heat enters the original gas
internal-energy owner exactly once.

Native wall/gas adjacency, dry gas occupancy, unique physical patches and the
shared face budget are checked before advancement. Solid film areas must fit
the original thermal owner's remaining exposed area after BOTH solid/solid and
solid/liquid contact coverage. The timestep bound combines film and native
contact conductances with each body's own caloric temperature slope; the gas
capacity bound includes all liquid and solid films together.

`body_environment` reports direct gas/solid heat and actual source samples.
`FreeSurfaceBodyExchange::heat_from_gas_j` is distinct from its
`heat_to_liquid_j = heat_from_gas_j + heat_from_radiation_j - solid_storage_change`. Direct gas/solid
exchange never enters the liquid's external-heat account. The finite-gas audit
closes actual liquid, solid and gas storage, retaining molecular and electrical
heat in their existing accounts. All native VOF, material, thermal, geometry,
budget and cancellation gates precede publication of liquid, solid and gas
candidates.

This remains a fixed-geometry, finite-cell solid and single-liquid,
constant-density VOF approximation with an explicitly declared well-mixed,
fixed-volume gas and quasi-steady continuum films. It supplies no gas momentum,
pressure traction, evaporation, reaction, moving-wall or phase-density handoff,
nor a total kinetic/capillary-energy closure. Film geometry is held during the
thermal stage and must describe the accepted configuration on each call.
Three focused G1/G3/G4 regressions cover independent finite-store balances,
solid-data changes reaching actual liquid populations through the common gas,
aggregate small-gas timestep control, geometry/coverage/domain refusal,
cancellation, retained histories and successful retry.


### Conducting solid surfaces with source-backed gas convection

`FreeSurfaceThermalBody::with_gas_convection` replaces the explicitly prescribed
gas film on that body's exposed faces with a borrowed
`FreeSurfaceBodyGasConvection` provider. The original conductivity key and
selected claim prepare one continuously supported pure-solid k(T) line inside
an explicit temperature span. The same original caloric chart, material
context and fixed mass remain authoritative. Body clones and thermal trials
share this preparation; they do not rebuild conductivity charts.

Each face requires positive area and solid centroid distance, with zero gas
thickness. The shared material kernel solves the supported conductivity
integral against convection at the actual surface temperature. Each provider
query receives the original signed wall/dry-gas geometry, outward normal,
finite bulk gas and original film sources at mean surface/bulk temperature and
current EOS pressure. Driven coefficients must be finite and positive. `None`
is admitted only at exact equilibrium. A source, correlation, arithmetic,
surface-budget or cancellation refusal aborts the original transaction.

The actual gas-side power enters finite solid storage and the same aggregate
opposite gas debit once. Existing explicit and implicit timelines, gas energy
residuals, capacity bounds and rollback retain their owners. The reported
secant conductance is a frozen heat-rate representation, not a nonlinear
derivative. `FreeSurfaceBodyGasSample::convection` retains surface temperature,
original coefficient evidence, shared source line, integral-coordinate power,
heat-continuity residual and iteration count.

Radiation sharing the same native face uses the coupled surface owner below.
Distinct radiating faces retain the original cell-temperature law. Prescribed gas-film behavior remains
available. This is a fixed-geometry, quasi-steady, pure-solid surface closure
coupled to the original well-mixed finite gas. Spatial, temporal and correlation
validity require their existing assessment; no new gas momentum or phase
handoff is supplied.

Four focused native regressions cover an independent affine-k surface and
finite-store energy oracle in both heat directions, source sharing, convection
feedback into native populations, equilibrium, distinct-face radiation and
atomic refusal. Source parsing and whitespace checks passed; native compilation
and execution remain pending the required constellation dependencies.

### Shared convection and radiation on conducting solid faces

A radiating patch with the same original body, native face and exact area as a
conducting convection patch now shares one solved surface temperature. Unequal
areas and conflicting body ownership refuse. The two-bath material kernel
balances the original k(T) integral against separate gas and radiant powers,
including opposing drives, equal cell/gas temperatures and zero net solid heat.
Original emissivity is queried at actual surface temperature with its existing
claim and material context; no new conductivity chart is built during trials.

The gas sample retains the accepted radiant sample. The radiant timeline reuses
that exact sample at applied and endpoint states, without another surface root.
Actual gas heat and external radiation enter their original energy accounts
once. Material/gas and material/radiant pair conductances supply the solid row;
the gas row additionally includes the gas/radiant-bath relay. Explicit capacity
limits and implicit gas conditioning both retain that path. These are frozen
secant scales, not derivatives of changing coefficients.

`FreeSurfaceBodyConvectionSample::solid_power_to_surface_w` retains the signed
material integral power; it is distinct from either bath's heat. Reports also
retain actual emitting temperature, original emissivity and coefficient
evidence, heat-continuity residual and the gas/radiant pair conductance.
The original fixed geometry, pure-solid profile, well-mixed gas and diffuse-gray
unit-view-factor radiation assumptions remain explicit.

Three native specifications cover independent affine-k/fourth-power surfaces
and finite liquid/solid/gas balances in both heat directions, the radiant relay,
exact equilibrium, actual surface-source refusal and late-provider rollback.
They passed source/API and syntax review; native compilation/execution remain
pending an admitted constellation build. The shared numerical kernel's 23-test
bounded execution is documented separately in fs-material.

### Independent fragment conservation in tetrahedral liquid ingress

The original `free_surface3_from_tetrahedra` and
`free_surface3_from_tetrahedra_with_enthalpy` paths identify face-connected
pieces through the original shared vertex indices. Each piece uses the SAME
tet/Cartesian overlap geometry, P1 integrals and native dense momentum
correction independently. A detached component cannot supply another's
linear or angular impulse during remapping.

`TetLiquidTransferReport::components` retains source tetrahedron indices,
receiving native cell coordinates, actual source/target moments, centroid
displacement and signed kinetic/gravitational projection changes for each
piece. Components are ordered by their lowest source tetrahedron index.
A receiving cell occupied by two disconnected pieces refuses: unresolved
coalescence cannot replace separate states. Repeated tet vertices and more
than two cells incident to one indexed face also refuse. This topological
admission is not a global geometric nonintersection certificate.

Mass and both momentum budgets bound sums of absolute component errors
(momentum uses Euclidean norms). The mechanical error budget bounds the sum
over pieces of absolute kinetic change plus absolute gravitational change.
Each component also meets the original centroid limit. Opposite errors cannot
cancel between fragments or independently spend the full shared allowance.

The enthalpy extension retains actual source and receiving mass*h accounts
per component, in the same order as the mechanical reports. Their summed
absolute discrepancies meet the original thermal allowance. Mechanical
projection changes remain numerical representation error and never become
heat. Source-card, pinned-claim, caloric, density, native scaling, work-budget
and cancellation gates remain in force.

These constructors still admit already-liquid, uniform-material snapshots
without stored elastic/plastic history, with explicitly neglected capillarity. They do not release solid
elastic energy, transport plastic history, create a pressure projection,
support different fluids in one cell or claim a complete melting transition.
The component mapping is ingress lineage, not persistent fragment tracking
after native VOF evolution or coalescence. Failed admission returns no
candidate and leaves supplied fields unchanged.

Three G1/G3/G4 regressions cover independently translating/counter-rotating
heated cubes with actual continued native flow, isolated-component
equivalence, equal/opposite numerical-energy errors that cannot cancel,
unresolved component aliasing, cancellation and retry.

## Source buoyancy through melting and complete native liquid release

`MaterialFlowFeedback::try_new_phase_boussinesq` attaches original source density
to the SAME native solid/mushy/liquid momentum and enthalpy timeline. The original
phase chart, reference mass/volume, exact-solid occupancy and sourced Darcy
resistance remain authoritative. Mobile material receives (rho_source-rho0)*g;
exact solids retain the opposite stationary support-on-material reaction, zero
velocity and zero support work. Liquid-only properties are never fabricated for
exact solids. Source density/Cp departures retain the existing Boussinesq bounds;
raw dynamic viscosity and conductivity remain physical source values.

The original unsteady momentum owner retains signed buoyancy work, kinetic
storage, molecular/Darcy heat and backward-Euler numerical loss separately. Only
physical dissipation heats the material. Every endpoint source or bound refusal
uses the existing atomic thermal/momentum publication and deterministic retry.

Complete melting can release the original phase source into the native VOF
Boussinesq owner while retaining the same caloric object, card, query conditions,
selected claims, reference transport, gravity and approximation bounds. Existing
mass, enthalpy, momentum, angular momentum, clock and projection-loss admission
remain in force. Removing stationary support exposes native reference gravity
plus its density perturbation; force-frame rebind preserves initial velocity.
The joint phase resolver admits original h/T/fraction, native physical shear
and bounded source density in one query. Shear-dependent density or Cp still
refuses. `BoussinesqPhaseMaterial::at_source` refreshes that source while retaining
the original reference, bounds, material, caloric chart and mobility law; exact
solids continue to omit liquid-only properties. Native viscosity drives the same
momentum, Darcy resistance and molecular heat, while the original k drives heat.

`release_fully_liquid_with_rheology` now carries this combined source into
`free_surface3_from_cells_with_boussinesq_rheology`, which uses the existing VOF
force and constitutive owners together. Source MAC rates admit the original
source; receiving rates are solved from its own force-corrected stress. Mass and
inertia retain rho0, departing source rho drives buoyancy, and force-frame
rebinding precedes all receiving momentum/kinetic checks. Projection loss is
still reported separately from heat. The joint regression heats exact solid
through the latent plateau, drives native shear, releases complete liquid and
checks data-only motion changes, original transfer budgets and independently
summed mu*gamma_dot^2*volume*dt heating on the evolving VOF owner.

This is a bounded Boussinesq stationary-skeleton model, not moving partial-solid
geometry, shrinkage, elastic-energy release, separate phase velocities or a
closed gravitational-potential inventory. Three native regressions cover cold
startup to mushy circulation/remelting, endpoint refusal, and source-preserving
complete release followed by further heated motion. They are authored but unrun
here. Independent exact-rational ring force/work arithmetic passed; it does not
replace native execution, refinement or experimental qualification.

## Conservative affine enthalpy ingress on the original thermal operator

`EnergyConfig::step_with_enthalpy` uses the SAME finite-volume conductivity,
contact, boundary, sink, storage and Krylov owners with an explicit frozen
`AdvectiveEnthalpyField`. Each fluid cell supplies h0, T0 and positive Cp from
one admitted material chart; solids have no fluid entry. Both sides of each
interior face use the same donor's complete h0+Cp*(T-T0), including its offset,
so transported energy is equal and opposite even when local Cp differs. Inlet
and backflow values are explicit enthalpies at their declared temperature.
Outflow has zero diffusion; Inflow retains the original half-cell diffusion.

This entry supports first-order upwind and opt-in limited source-enthalpy
reconstruction, with one constant positive density, locally volume-conservative
flow, finite coefficients and impermeable solid faces. Reference FluidProperties Cp does not replace supplied caloric tangents.
The original conductivity-only and constant-Cp entries retain their existing
numerical scheme. Public VoxelDomain adjacency permits coupled material owners
to audit the same physical faces without constructing another grid topology.

The returned native balance describes the frozen linearization. Its open-face
audit evaluates centered enthalpy and temperature differences to reduce affine
cancellation. The operator subtracts the first fluid cell's h0 from every donor
and inlet/backflow enthalpy before forming affine intercepts, limiter corrections
and boundary fluxes. Local continuity makes this one common reference cancel
from each physical balance; its admission is unchanged. This avoids cancellation
of large formation-energy fluxes hiding small heat sources and uses no additional
per-cell allocation. Net boundary transport is reference-independent when net
boundary mass flow vanishes. With fixed-temperature fluid cells on open faces,
the separately reported fixed-cell injection and nonfixed boundary outflow use
this same centered convention; their combined physical balance is invariant.
The nonlinear material audit likewise subtracts the common h0 before its Cp*dT
prediction, so a small actual face-law defect cannot round to a false zero.
Focused regressions retain a small source and detect a sub-ULP affine mismatch
under positive/negative common datum shifts, both flow directions and both
supported enthalpy reconstructions.
Nonlinear callers must independently check true chart storage and
face fluxes before accepting state; neither raw coefficients nor a small global
linear residual establish material-source validity or nonlinear conservation.
Three native regressions cover differing donor Cp, reverse flow, common-datum
invariance, real outlet backflow, nonconservative-flow/invalid-field refusal and
cancellation. They are authored but unrun here; pinned syntax checks, static API
review and independent three-cell arithmetic passed.

## Physical storage admission in combined buoyancy and rheology

The shared moving-liquid rheology binding rejects a density or Cp source with
the evolving physical-rate axis, even if its value at the reference rate equals
the original storage coefficient. A coincident value cannot grant a missing
shear-induced volume or caloric-state model. This one admission applies to
direct VOF initialization, complete-melt release and conservative tetrahedral
ingress through their existing shared binding. Temperature-dependent density
remains supported through the original bounded Boussinesq projection.

Three additional native VOF regressions cover source-density and shear-viscosity
substitutions in continued motion and molecular heat, an independent
force-corrected constitutive root, density-bound and cancellation rollback, and
rejection of rate-dependent storage at equal reference values. The separate
constitutive/force/heat arithmetic closes to roundoff. Native Rust compilation
and execution of these added tests remain unavailable in this environment.

### Ordinary immersed dry body with measured interior-fluid feedback

The optional one-body `step_with_dry_partial_saturation_body` uses the original
positive dry mass and complete tensor for light, neutral and heavy solids. It
solves the six-component retained-interior momentum correction implicitly through
the unchanged NT collision, streaming and native rigid split. Endpoint interior
inventories come from actual streamed populations and sampled occupancy; no
analytic rigid-fluid tensor subtraction, epsilon mass, or density-matched bypass
is used. Inventories use physical Guo momentum `sum(f*c) + F/2`, including its
endpoint-occupancy weighting. Uniform and per-cell fluid forces are admitted and
held during each candidate solve. The correction is
`J_c = Delta(Pi_inside) - J_inside,F`, where `J_inside,F = sum(epsilon_mid * F)`;
angular components use the same fixed world origin. This excludes forcing of
retained fictitious fluid from the physical body's load. The midpoint source is
reduced once from the first native candidate's existing sampled geometry and
cached across Newton probes, without extra quadrature or population storage.

Cut-cell angular inventories retain the occupied subvolume's first spatial
moment `s = integral_inside (x - x_cell) dV` from exactly those midpoint samples.
For cell momentum density `p`, the interior contribution is
`epsilon * ((x_cell - origin) cross p) + s cross p`; midpoint imposed-force
torque uses the identical geometry with `F` replacing `p`. Thus an asymmetric
partial cell does not relocate its fluid momentum or force to the cell center.
The first moment is a volume integral, not a centroid to multiply by occupancy
again. Mass fractions and linear-momentum quadrature are unchanged. Integer
sample-offset sums give exactly zero first moments for full and symmetric
cells. Only nonzero cut-cell moments are stored, sorted into the native
tile/lane order and consumed by cursors; no second grid-sized field or per-cell
lookup is needed. The shared dry-body owner supplies this convention to
single-body, simultaneous-body and lubricated-sphere paths. Changing the fixed
origin shifts each moment by minus the origin displacement crossed with its
linear inventory. These are moments of the declared finite-resolution
quadrature, not an exact curved-surface or moving-boundary GCL claim.

Full-pressure gravity uses the dry-body load `M_s*g` alongside the configured
fluid source. Reduced-pressure gravity uses `(M_s-rho_f V)g` with the hydrostatic
fluid source removed; these conventions must not be combined. An evolved force
field uses the existing `rebind_force_field` to preserve its Guo source frame.
The represented fluid source is full-fluid forcing minus interior forcing;
authored body work and the solved interior-correction work remain separate.

Finite-difference probes, Newton corrections and backtracking have explicit hard
work limits and condition admission. Probes read population inventories only in
occupied cells; independent full-grid linear and fixed-origin angular inventories
run before and after convergence. Body and all-fluid minus interior plus body
momentum close against authored impulse, represented fluid forcing,
stationary-wall exchange and periodic angular transport before the single
publication. Failure/cancellation preserves populations, dry pose/momenta and
topology admission; only scratch may change.

Controls bound initial/midpoint/endpoint sampled-volume variation (default 1%,
at most 10% per step), and refuse zero resolved interior volume. Actual interior
mass, slip, quadrature volumes, correction work and macroscopic kinetic-energy
changes remain explicit. These diagnostics do not certify exact no-slip, a GCL,
complete thermodynamic energy, calibrated drag, coupled temporal order or
unconditional density-ratio stability. The isothermal, fully immersed lane has
no moving free surface. Focused regressions cover repeated nonuniform forcing
and torque for light, neutral and heavy bodies with periodic and stationary-wall
boundaries, source subtraction for force confined to the fictitious interior,
physical source-frame inventories, and forced rollback with exact retry.

### Endpoint fluid-force histories through one physical dry body

`step_with_dry_partial_saturation_body_and_endpoint_force_scale` carries the
original selected spatial force shape through the SAME positive-dry-inertia
Newton solve, native rigid split and streamed fluid candidate. Collision and
reciprocal reflection retain the accepted starting force. Endpoint physical
Guo inventories and original open-face reconstruction use the actual next force;
the unchanged uniform contribution remains held. Full-fluid source impulse is
`(F_old + F_next)/2` per cell, and its local work uses physical velocities in
their corresponding source frames.

The fictitious-interior subtraction integrates that SAME trapezoidal source
over the existing collision-midpoint occupancy and cut-cell first moments.
It is cached once across probes. Candidate endpoint inventories use the next
force before it is published, so its half-source momentum cannot masquerade as
an additional physical body load. Independent body and exterior-fluid-plus-body
linear/angular accounts retain wall, periodic and reservoir contributions.
All candidates use the original accepted field and reusable post/next storage;
no full-grid backup, second collision, or force-field clone occurs per probe.
Force, populations, pose and momenta publish together after final admission.
Legacy held callers pass no endpoint scale and retain their original arithmetic.

The frozen-geometry `step_with_partial_saturation_body_and_endpoint_force_scale`
uses the same source/reconstruction path. These are prescribed smooth fixed-shape
histories, initialized with momentum-preserving field rebinding. They do not
supply an implicit velocity-dependent load, a discontinuous force switch,
complete coupled temporal order, or a thermodynamic fluid-energy/heat balance.
Held-force entry points retain their original semantics; the simultaneous
endpoint-force extension is specified below.

Focused regressions specify zero-occupancy native parity across four collision
models and periodic/wall/reservoir boundaries, independently reconstructed
physical source work, light/neutral/heavy moving-body ramps and reversal through
zero, cut-cell interior subtraction, exact held-path parity and atomic
overflow/budget/late-cancellation refusal with retry. The native owners now
compile and execute in the assembly selection below; the standalone
`dry_partial_saturation_body3` target and full workspace were not run.

### Endpoint fluid-force histories through simultaneous immersed assemblies

`step_with_partial_saturation_bodies_and_endpoint_force_scale`, the moving
assembly/sphere-film entries, and the physical dry assembly/sphere/prescribed-
plane-film entries carry one selected endpoint source through the original
shared collision, rigid split and resistance solve. Collision and reciprocal
body response retain the accepted old force; final physical Guo inventories
and velocity/pressure reservoir reconstruction use the next force. The shared
full-fluid source integrates their trapezoidal impulse and angular impulse
once. Dry assemblies subtract that same source over each distinct collision-
midpoint interior, including cut-cell first moments, cached across the full
coupled 6N Newton probes. The original positive masses, full tensors,
conditioning, candidate limits, geometry and film law remain authoritative.

Candidate probes never publish force or populations. All bodies, fluid and
endpoint source publish only after final trajectory, film, population and
independent momentum accounts. None callers and singleton delegation preserve
their original arithmetic; no field clone or additional fluid stream is
introduced per probe. New material assembly entries cover retained and physical
dry sphere/ellipsoid consumers. Moving plane-film and finite/shared-carrier
consumers are covered below; this adds no complete coupled-order,
implicit constitutive force, GCL or thermodynamic work claim.

Focused selected native execution passed 39 tests across
`dry_partial_saturation_bodies3`, `forced_partial_saturation3` and
`dry_lubricated_spheres3`. Eight new checks reconstruct actual fluid-plus-body
linear/angular inventories, distinct cut-cell source subtraction, open/wall/
periodic transport, analytic ramp-source history, active pair/plane films,
constant/zero/singleton parity and late refusal with exact retry. Original serial
numerical owners and actual cached dependencies compile; unrelated pooled/SGS/
frontend graphs are excluded. A subsequent current-source build retained the
concurrent local-viscosity curved-wall owner and passed seven endpoint/history
cases plus the original curved-wall halfway parity case. The full Cargo
workspace and material source-card constructors were not executed.

### Endpoint force histories with moving finite and shared film carriers

The moving `*_and_surfaces_and_endpoint_force_scale` and
`*_and_carriers_and_endpoint_force_scale` entries carry the original prescribed
source history through retained-fluid and physical dry-body transactions.
Collision uses the accepted old force; endpoint Guo inventories and open
velocity/pressure reconstruction use the next force. Every dry probe starts
from the same accepted surface/registry state and reuses the cached midpoint
fictitious-interior source. Finite recoil and cross-patch resistance therefore
participate inside the accepted coupled candidate.

A registry carrier owns one original mass/full tensor, authored load, aggregate
film wrench and kinetic inventory regardless of its attached face count.
Film-only carriers acquire no sampled-fluid interior or extra source
subtraction. Prescribed reactions remain external ports; finite reaction work
belongs to actual carrier motion. All final trajectory, sweep, allocation and
momentum checks precede the shared source/fluid/body publication; subsequent
surface and registry copies are infallible. Held callers preserve their None
arithmetic.

Sixteen focused selected-native checks passed, including seven new groups.
Independent population/Guo inventories close linear and angular momentum for
a dry sphere, two unique carriers, three attached faces and a prescribed
spinning surface under source ramps/reversal across periodic, wall and
pressure-reservoir boundaries. Retained-body controls prove the first
collision retains old-force recoil and later source history changes finite
motion. Constant/zero parity, candidate refusal, late finite-surface crossing
and exact retry preserve source and every accepted owner. These tests use
the current original serial numerical owners and genuine cached dependencies;
full Cargo/DSR and the material constructor graph remain unexecuted.
Frozen-geometry finite entries and the original film-only geometric/thermal
limitations are unchanged.

### Simultaneous physical dry-body assemblies

`step_with_dry_partial_saturation_bodies` advances one to eight bodies with
their original positive dry masses and complete inertia tensors. The one-body
case delegates exactly to the existing dry-body path. For multiple bodies, one
native simultaneous collision/streaming candidate supplies the full coupled
`6N` interior-momentum correction; off-diagonal derivatives retain interaction
through shared streamed fluid. Bodies are not advanced by sequential fluid
steps. Physical Guo inventories and each body's collision-midpoint interior
source retain the one-body convention. Each fictitious interior source is
subtracted once from the represented exterior-fluid source.

Correction convergence is checked in both body-COM and fixed-world-origin
coordinates, so a long moment arm cannot hide an angular defect. Independent
body and combined exterior-fluid-plus-dry-body momentum inventories gate a
single publication of all bodies and populations. Failed probes, exhausted
budgets, ill conditioning and late cancellation leave accepted state and
topology admission unchanged. Shared wall, forcing and periodic receipts belong
to the assembly, rather than being duplicated in each nested body receipt.

The bounded lane requires disjoint swept enclosing spheres, nonoverlapping
sampled ownership, resolved interiors and the existing volume-variation limits.
Each Newton iteration uses `6N` full assembly probes plus trial/backtracking
evaluations, bounded by the candidate budget; at most a 48-by-48 correction
Jacobian is retained. Existing population scratch is reused without cloning the
grid. Contact/lubrication, a moving free surface, exact GCL, calibrated drag and
unconditional added-mass stability remain outside this increment.

Nine focused native regressions cover same-candidate cross-body transmission,
forced full-tensor dynamics at a translated origin, distinct source subtraction,
permutation covariance, three-body translation, exact singleton delegation and
atomic refusal/retry. All six existing single-dry regressions also pass.

### Physical dry spheres with simultaneous missing-film resistance

`step_with_dry_lubricated_spheres` and
`step_with_dry_lubricated_spheres_and_planes` compose the existing normal and
optional reciprocal tangent/spin resistance with the bounded physical dry-body
correction. Every full `6N` candidate retains original positive dry masses and
complete tensors, solves the native films simultaneously, and uses the resulting
shared collision/streaming state for actual interior Guo momentum. Each midpoint
fictitious-fluid source is still subtracted exactly once.

Native film impulses are separate from resolved-fluid hydrodynamic loads and
are not reapplied to populations. Individual physical body balances include the
film wrench once. The combined exterior-fluid-plus-dry account includes the
opposite prescribed-plane reaction about the same fixed origin, including
reaction couples; internal sphere-pair wrenches cancel. Individual plane
reactions are exposed by `DryLubricatedPlanesExchange3`. Plane work is work
received by the prescribed surface, with its reservoir supplying the negative.
The native nonnegative film-loss receipt remains separate from any thermal
destination.

Plane geometry remains prescribed and supplies only unresolved-film resistance;
it neither moves resolved fluid boundaries nor provides finite-inertia recoil.
Existing film/mode limits, positive-gap admission, tangent asymptotic limits,
dry candidate budgets, conditioning, volume checks and atomic publication remain
authoritative. Empty planes preserve the sphere-only path; a singleton without
planes delegates to the original single-dry owner. No contact law, GCL, exact
impermeability or unconditional added-mass stability is added.

Five new native groups pass, including unequal-density reciprocal resistance,
translated-origin tangent/spin plane reactions, nonuniform forcing with distinct
interior sources and 15 simultaneous pair/plane modes, inactive-film exact
compatibility, and post-stream endpoint refusal with exact retry. All nine
existing dry-assembly and six single-dry groups also pass. Full Cargo
qualification remains blocked by the unavailable sibling runtime dependencies.


## Integer-cell planar piston geometry epochs (MR69 foundation)

\`FreeSurface3::remap_planar_piston_one_cell\` adds a bounded moving-geometry
primitive for the exact plug-translation limit. One complete planar piston can
advance or retreat by exactly one Cartesian cell. Every wet VOF classification,
tracked interface mass, native D3Q19 population vector, local force frame and
local surface coefficient is translated one-to-one. Newly uncovered cells
therefore receive the conserved prior liquid state; they are never reset to an
arbitrary equilibrium. Covered piston slabs are wall geometry for neighbor,
sampling and fill-ghost queries even while their sparse tiles remain allocated.

The remap is transactional. It activates every destination first, refuses domain
escape, permanent-solid overlap, destination aliasing, incompatible wall/contact
geometry, mass drift, momentum drift, cancellation and epoch overflow, then
publishes one candidate. \`geometry_epoch\` increments only on an accepted geometry
publication; ordinary VOF timesteps leave it unchanged. Advance followed by
retreat restores the physical VOF mass/classification snapshot, while gas
population scratch remains intentionally unspecified.

This operation is a conservative geometry remap, not a normal-wall pressure
solve. Population momentum is retained up to floating reduction error, so its
reported remap impulse and work are zero. A finite-speed piston still requires a
pressure/traction dynamics owner to supply physical impulse/work. Subcell cut
volumes, deforming nonplanar solids, moving contact lines and simultaneous
phase-mass transfer remain open MR69/MR29 work. Focused G1/G3/G4 regressions
cover all six wall orientations, nonzero carried momentum, exact advance/retreat
state, stale-wall refusal, cancellation and full-domain overflow.

## Endpoint-coupled material buoyancy on the native transient clock

`MaterialChtFluid::march_natural_convection_implicit` and
`march_natural_convection_with_enthalpy_implicit` close original source viscosity,
conductivity and density-derived buoyancy against each accepted native thermal
endpoint. The first retains reference-Cp thermal storage; the second composes
the existing source h(T) storage and conservative donor-enthalpy transport.
Original density/Cp departure limits and source-domain admission still apply.

Each bounded nonlinear pass solves the same original BE/BDF2 momentum interval
and backward-Euler heat interval. Accepted previous/older velocity levels and
old temperature remain fixed throughout those trials. Current velocity and
pressure may warm-start the next native solve, but do not become new physical
history. Flow projection, thermal assembly, caloric iteration, compact storage,
contacts, fixed solid properties, thermostats, scheduled sources and endpoint
radiation remain on their existing owners.

Acceptance independently bounds the actual unrelaxed temperature change, every
applied-versus-endpoint viscosity/conductivity value and every buoyancy-force
component. A small temperature increment cannot hide a sharp material response.
Temperature relaxation affects only the next source-temperature iterate. One
accepted physical step adds one flow/heat record and one average sample, then
advances native momentum time levels once. The returned original result retains
applied and endpoint source samples and one closure report per physical step.

`MaterialEndpointConfig` bounds passes per step and whole-march momentum solves,
actual linear energy solves and transport-source refreshes. Discarded material,
caloric, radiation and compact-storage trials consume these limits. The
source-caloric entry debits both its original heat budget and the outer budget
by actual solves; the smaller remaining allowance governs the next native call.
Repeated outer trials accumulate their caloric work into one physical-step
report. Source refusal is terminal; native failure, cancellation, exhausted
work or failed closure returns no partial march or accepted history.

Regression sources use an independent nonlinear-conductivity backward-Euler
solution, signed Darcy BE/BDF2 recurrences, true h(T) thermostat storage,
open-channel enthalpy recurrence and datum shifts. They check single history/
average advancement, sharp material-response refusal, late source failure and
joint caloric/radiation work limits. The new entries are opt-in laminar source
coupling, with first-order coupled time accuracy even when momentum uses BDF2.
Implicit SGS, source-solid calorics, wall closures and viscous-heating
composition are not selected by these entries. No unconditional nonlinear
convergence, measured performance or experimental accuracy is claimed.

## Prescribed gas pressure resolved on actual free-surface links

`FreeSurface3::step_with_gas_link_pressure` supplies a distinct pressure to each
actual missing incoming gas population. `GasPressureLink3` identifies the liquid
cell, signed gas source coordinates and incoming native D3Q19 direction. The
source is exactly `liquid - E3[q]`; callbacks visit native cells deterministically
and gas directions in increasing q. Inactive, covered-piston and out-of-domain
wall sources never become gas pressure links. Separate neighbors of one liquid
cell can therefore carry separate prescribed loads without assigning one
pressure to that entire cell.

The original reconstruction remains the only momentum boundary owner. Every
link adds its gas-pressure perturbation to the SAME cell capillary density and
retains the existing physical Guo source correction. The combined gas/capillary
density must remain positive and inside the explicit excursion bound. Curvature
is evaluated once per exposed cell; repeated identical boundary densities reuse
the same equilibrium without a map or per-link allocation. Existing cell fields
are still sampled exactly once per exposed cell, and uniform link loading keeps
their population arithmetic and bits. Interior cells do not evaluate curvature
or initialize a gas equilibrium.

The frozen SI, source pressure/shear rheology, heat/enthalpy, complete conducting
thermal and native/SI tangential moving-wall entries share this carrier and
original collision/transport step. `MaterialGasPressureLink3` adds both cell
centers and their midpoint in metres. That midpoint is not an asserted VOF
surface intersection; the streaming direction is not a geometric surface
normal. Absolute positive Pascals use the original fixed datum and pressure
unit. The composed thermal link entry refuses a simultaneous cell-pressure
owner. All existing source, surface/electrical, body heat, wall work and atomic
publication gates remain in force. The original native Local/LegacyGlobal
conversion-error semantics are unchanged; ConservativeLocal and the existing
material transactions publish only complete candidates.

Four native regression sources independently check every reconstructed
population and pressure-traction sign on all axes, exact uniform-load parity
including force/capillarity and callback cadence, wall exclusion, and last-link
pressure/cancellation rollback. Four physical regression sources check SI
coordinates and gauge, pressure/shear source updates, heated mass-route parity,
and complete finite-body thermal rollback; the existing wall-pressure regression
also compares both link entrypoints and reciprocal reports. Native execution
remains unverified in the resource-constrained delivery environment. This is
prescribed gas traction, not evolved gas topology, gas-volume reconstruction,
gas inertia, compressive heating or a pressure-volume work closure.

## Implicit material convection with native SGS transport

The opt-in `march_natural_convection_with_sgs_implicit` and
`march_natural_convection_with_enthalpy_and_sgs_implicit` entries compose
`MaterialEndpointConfig` with the original projected momentum, SGS and thermal
owners. Every outer trial starts from the same accepted thermal and momentum
histories. Its current projected eddy viscosity replaces the reusable eddy heat
field; it is never accumulated across trials or written into molecular material
coefficients. Reference-caloric heat uses rho0*Cp0*nu_t/Pr_t. Original source
calorics use the existing actual-Cp(T) SGS heat owner and its independent face,
cell, storage and enthalpy-flux closure tests.

The physical local eddy-viscosity limit is checked against the actual refreshed
endpoint molecular viscosity before acceptance, even when a loose nonlinear
material tolerance would otherwise accept a large source change. This check
does not advance SGS or alter its accepted field. Every actual momentum and
heat solve consumes its existing whole-march allowance, including discarded
outer and caloric trials. Coupled time accuracy remains first order with
BE or BDF2 momentum. The original exploration-only SGS admission remains;
these entries add neither a turbulent wall model nor dissipation heat.

Regression sources exercise nonzero SGS, exact zero-SGS specialization,
independent native heat replay with one eddy field, source-Cp(T) heat transfer,
and physical-cap/heat-flux/work refusal. They were syntax-checked; native
execution was unavailable in the delivery environment.

## Implicit source-solid conjugate convection

The `march_natural_convection_with_material_solids_implicit` and
`march_natural_convection_with_enthalpy_and_material_solids_implicit` entries,
and their `and_sgs_implicit` counterparts, close original solid k(T), integrated
solid Cp(T), fluid transport and buoyancy on the same accepted native endpoint.
The enthalpy entry also retains the original fluid h(T) and its source-owned
storage/advection defects. These wrappers reuse the existing solid and joint
caloric heat solvers; no independent time integrator or synthetic solid source
is introduced. Solid conductivity is retained for endpoint-defect evaluation
while the fluid source refresh replaces only fluid coefficients.

Every nonlinear trial measures storage from the same previous accepted
temperature. Repeated trials at one physical time aggregate constitutive and
linear work while replacing the provisional storage record with the accepted
record. Source energy datums, reference masses, contacts, radiation, compact
storage, fixed-temperature reservoirs and reciprocal interfaces retain their
original owners. The plain and joint caloric paths consume the shared outer
energy allowance on every actual solve. A source/native refusal, cancellation,
failed endpoint closure or exhausted allowance returns no partial march.

SGS variants use the current outer trial's projected eddy field only in fluid
heat transport; joint source calorics convert it with actual fluid Cp(T).
Original SGS limits and flux-defect gates remain active. Coupled time accuracy
remains first order; these APIs do not select viscous or SGS dissipation heating.
Four regression sources cover an independent nonlinear two-cell heat balance,
same-step buoyancy from solid heating, exact work-budget/source refusal, and
joint SGS/solid energy accounting with constant-Cp specialization. They were
syntax-checked; native execution was unavailable in the delivery environment.

## Dense source-pinned moving wetting

`MaterialFreeSurface2::with_dynamic_planar_wetting` binds the original ordered
solid/liquid interface to one complete stationary X/Y wall. The existing dense
curvature and capillary-pressure reconstruction consume a complete spatial
angle field with direct O(1) tangential lookup. Advancing/receding endpoint
selection, immutable source claims, speed dimensions, fixed source conditions,
deadband and optional full-angle two-phase Cox law share their owner with sparse
3-D and thermal wetting. Measured speed-dependent angle curves cannot also use
Cox rate dependence; stationary receipts stay at zero speed for Cox predictions.

Contact-line speed comes from consecutive accepted first-row VOF fields in
physical m/s. A local tangential average preserves affine planes and takes
O(wall cells) work. It is independent of imposed angle ghosts, so changing an
angle cannot itself fabricate line motion. The geometric normal uses the next
inward row. This first-cell/backward-time estimate drives the next native step;
the deadband projects the geometric angle into the sourced interval without
claiming exact pinning or calibrated slip. Bulk flow is isothermal, single-liquid
and invariant along the missing dimension; mass remains kg per metre of depth.

`MaterialFreeSurfaceWetting2::step_with_gate` works on a private native candidate
and publishes populations, mass, carry, conversions, angle field and source
receipts together. Source-domain exit, degenerate moving geometry, Cox limits
and cancellation retain the last accepted state. Cancellation is polled around
the dense clone/kernel and at each wall query; the original dense sweep is not
internally interruptible, and retains its numerical panic semantics. No surface
energy closure, heat partition, moving substrate or second-fluid momentum is
invented. Constant-angle trajectories reproduce the original equilibrium ghost;
zero surface tension leaves the angle field inactive in momentum.

## Finite-atmosphere EOS pressure driving free-surface momentum

The opt-in `step_with_pressurized_finite_gas_environment` and
`step_with_pressurized_finite_gas_surface_exchange` close the original
finite-atmosphere heat-to-pressure-to-liquid-momentum path. The SAME
`IsochoricGasParcel` receives the aggregate opposite liquid/solid heat on
each existing thermal substep. After all substeps, its actual endpoint EOS
pressure supplies the original native missing-gas-link reconstruction over
one flow/advection stage. The initial gas pressure is not reused as a stale
mechanical load. Prescribed and finite-atmosphere pressure owners cannot
both supply that boundary.

No EOS, caloric source, momentum solver or second physical clock is added.
Finite solid conduction/radiation, solved interface temperature, electrical
and thermocapillary bindings retain the existing shared transaction and
source-domain checks. Original gas/capillary density-excursion and native
flow/transport gates run before publication; late pressure refusal therefore
preserves gas, liquid, finite solids, sources and clock together.
`FreeSurfaceFiniteGasPressureLoad` reports the original initial and actual
applied EOS pressures, fixed parcel volume and declared native pressure bound.

All exposed gas links must communicate with one well-mixed atmosphere.
The caller owns that physical connectivity and the parcel's fixed physical
volume, which may include an explicitly connected rigid plenum. A fixed
rigid vessel and fixed amount of incompressible liquid are the intended
geometry; moving-solid epochs refuse. Numerical lattice density and
`mass/rho` fill never evolve physical gas volume. The original first-order
thermal/flow split and fixed source-pressure liquid caloric scope remain.
There is no gas momentum, inferred pocket topology, compression work,
evaporation or complete mechanical-energy/passivity claim. The combined
thermal storage audit does not become a whole-fluid kinetic-energy audit.

Four authored regressions independently calculate caloric pressure change
and planar-link traction for heating/cooling, exercise radiant surface heat
feeding gas pressure, preserve formation-energy gauges and finite-volume
causality, and check late pressure refusal, requested cancellation and exact
continuation of gas/liquid/solid storage. The test vessel's declared gas
volume comes from its initial rigid geometry and physical liquid amount.
Native compilation, rustfmt and test execution were unavailable after the
worker filesystem/exec service failed during this implementation.

## Implicit molecular viscous heating with material calorics

The opt-in
`march_natural_convection_with_enthalpy_and_viscous_heat_implicit` and
`march_natural_convection_with_enthalpy_and_material_solids_and_viscous_heat_implicit`
entries compose the original molecular viscous-heating owner with implicit
material/buoyancy endpoint feedback. They retain the original backward-Euler
unsteady-Stokes trajectory, fluid h(T) owner and source-solid caloric/exchange
owners. Accepted time histories stay frozen across nonlinear trials. Each
trial replaces its provisional physical heat deposit; one accepted mechanical
receipt enters each physical step's audit/history.

Acceptance independently samples the actual endpoint viscosity and force
through the original Reynolds, physical heat-partition and kinetic/work laws.
The non-cancelling L1 local thermal-credit/mechanical-debit defect must fit the
original power tolerance, the mechanical balance must fit that tolerance times
the duration, and the original frozen momentum-row residual must meet the
native inner tolerance. Deposited and accepted endpoint mechanical heat agree
within the original power tolerance times duration, not necessarily bitwise.
Loose outer material tolerances cannot waive these physical checks.

Discarded thermal trials remain charged to existing actual-work budgets.
Fixed-temperature molecular heat retains reciprocal reservoir accounting.
Source solids receive heat through their original exchanges and receive no
direct fluid viscous deposit. The original explicit path shares the same
strict mechanical sampler and physical-deposit helper.

Existing closed-wall/symmetry, Reynolds, mechanical-power and backward-Euler
restrictions remain. These APIs add no SGS viscous heating, numerical-loss
heating, new momentum integrator, higher coupled time order or buoyancy
potential-energy storage claim. Four regression sources cover independent
endpoint heating, physical admission under loose outer tolerances, native
history/work budgets and joint fluid/source-solid energy accounting.
They were source/API reviewed; parsing, compilation and native execution
were unavailable in the delivery environment.

## Limited source-enthalpy convection on the native finite-volume operator

`ConvectionScheme::EnthalpyMinmod` opts an original source-caloric energy solve
into conservative limited linear reconstruction. It is available on
`solve_with_enthalpy`, `step_with_enthalpy`, and the existing steady and
transient material-convection entrypoints. Native constant-Cp temperature solves
also accept it, including steady, transient, buoyant and resolved-conductivity
entrypoints. In that lane Cp factors out of minmod, so the existing temperature
operator transports the limited donor temperature with rho Cp Q. It does not
infer a nonlinear caloric chart. Momentum owners still refuse this selection.

For each fluid cell and coordinate axis, minmod limits the two specific-enthalpy
gradients using actual adjacent cell-centre distances. A donor half-width gives
the face correction. Two fluid neighbours are required; boundaries, solid
interfaces and unresolved one-cell directions retain donor transport. Both
rows of an interior face use the same signed mass flux, donor tangent and
correction, so the transfer has one credit and one debit. The existing positive
upwind Cp matrix, molecular/contact diffusion, source terms, reservoir laws,
reference mass and native caloric storage remain authoritative.

Three offsets per cell are cached once for each frozen energy operator, outside
its matrix and Krylov work. One graded-grid sampler serves native temperature
and source-enthalpy reconstruction. Native temperature transport starts from
the admitted upwind solution, refreshes its face correction for at most 200
bounded solves, and accepts only when the unchanged endpoint passes the original
algebraic and independent local heat residuals. Every trial retains the same
old-time storage, sources, reservoirs, contacts and compact junctions; internal
transfer remains paired. A small iterate change alone cannot publish a result.
The original material iteration refreshes source h,
Cp and the limiter. Independent acceptance re-evaluates the full reconstructed
endpoint heat on every active face and net cell row, including the mismatch
from a stale limiter when the caloric tangent itself is exact. Internal errors
cannot cancel into a falsely accepted global energy balance. Original per-pass
and whole-march work limits, cancellation, source support and physical gates
remain in effect. Nonfinite reconstruction or face coefficients refuse.

The same check is used by steady material convection and by native transient
fluid/source-solid caloric owners, including implicit material feedback, SGS
transport and molecular viscous heating where those owners are already
admitted. Boundary enthalpy keeps the original inlet/backflow source and common
datum. Existing thermostat and thermal-ledger corrections retain the original
physical accounting. No new thermal clock, storage law or momentum solver is
introduced.

This is a limited spatial reconstruction of specific enthalpy; backward-Euler
heat and existing coupled temporal orders do not change. Adjacent face traces
are bounded by their source enthalpies at the reconstruction state. This is
not a general maximum-principle, unconditional nonlinear-convergence or
certified-error claim for multidimensional coupled sources. First-order
boundary treatment and caller-selected mesh/time refinement remain explicit.
Four regression sources cover an independent graded-mesh linear-enthalpy
balance with signed flow, rest and common-datum shifts; a native nonlinear
caloric march and independent stored-energy accounting; native constant-Cp
admission with momentum refusal; and stale-limiter rejection even with an
exact affine caloric law. The native temperature regression source additionally
specifies oblique manufactured-solution refinement, reversed flow, insulated
circulation energy, uniform equilibrium and cancellation. The cooling-cht CLI
accepts energy_convection=enthalpy-minmod without fluid_enthalpy and reports
the constant-Cp scope. Rust parsing, compilation and native numerical execution
were unavailable in the delivery environment; source/API review and separate
numerical stencil calculations do not establish executed Rust results.

## Thermal density change with authoritative material momentum

`FlowField::advance_piecewise_thermal_volume_momentum` couples the original
chart-derived thermal/phase volume source to current material mass and
extensive momentum. The pressure-profile variant retains native reservoir
admission and component-relative pressure. The original density chart owns
sensible expansion, latent density change and nonlinear mixing volume; there
is no additional material mass source or fitted density law.

Current M/P reconstruction and body-force prediction run once. Every bounded
thermal-volume trial uses that same predictor and once-applied net caloric
increment. Only the accepted trial supplies pressure impulse to the original
cell momenta. The same accepted mass routes then transport the kicked momentum
with zero additional impulse. Rejected trials do not advance physical time,
and the composition adds no second pressure solve or projection clone.

The returned mass, specific enthalpy, phase state and extensive momentum are
the next inventory. Original controls, shared duration, source bounds,
incoming material and velocity, local/global balances and cancellation remain
required. A late refusal returns no partial endpoint. Pressure-source and
reservoir work, projection/scatter losses and donor transport losses retain
their separate native accounts; none silently becomes material heat.

Phase mobility, solid support and the caloric/work convention remain explicit
caller responsibilities. This first-order split supplies no solid stress,
free-surface geometry, viscosity or pressure-dependent thermodynamic EOS.
Analytic regression sources cover sealed expansion/contraction on all axes,
continued evolution from actual endpoint momentum, reservoir gauge invariance
and atomic late refusal. Compiler checks do not constitute numerical execution.

## Native shear and source rheology on evolving liquid density

`FlowField::advance_piecewise_thermal_volume_momentum_with_source_rheology`
and its pressure-profile variant observe cell RMS `sqrt(2 dev(D):dev(D))`
on the same original M/P reconstruction used by molecular viscosity. The
existing unit-molecular-coefficient stress quadrature retains native graded
edges, solid obstacles, walls and reservoir mouths. The observation solves
no momentum equation and adds no independent strain discretization.

The resolver receives each original phase state and physical rate once.
The standard material binding is
`resolve_phase_transport_state_with_shear_for_volume_change`; the caller's
resolver must consume those coordinates and retain its fixed source context.
Returned phase state and source density must match the original chart exactly.
Kinematic viscosity therefore converts with the actual admitted density.
The result retains original rates and resolved source snapshots alongside the
existing viscous thermal-volume momentum endpoint.

This entry admits fully liquid mobile material with zero declared Darcy
resistance. Original cells, declared incoming material, the actual molecular-
heated intermediate states and transported endpoints must remain liquid.
A phase needing solid or mushy support refuses atomically. Sensible expansion,
the original volume constraint, mass/momentum routes, shared heat-error gate
and separate numerical/mechanical accounts retain their existing owners.
Each mobile cell is queried once and the native viscosity step runs once;
thermal-volume trials repeat neither. Source viscosity is frozen at the
original rate, so time refinement is required for nonlinear rheology.

The focused native rheology regressions use a manufactured source-backed
melting chart with changing sensible density and kinematic shear viscosity,
constant-source equivalence, phase/source refusals, cancellation and retry.

## Molecular viscosity driving thermal density and phase volume

`FlowField::advance_piecewise_thermal_volume_momentum_with_viscosity` and its
pressure-profile variant route original molecular dissipation through the
same chart that owns material density and occupied volume. The existing
symmetric-stress viscosity owner samples each transport source once, uses
actual M/V face inertia and returns its native momentum impulse and physical
cell heat. The sampled source density must exactly match the original chart
state; other source coordinates and admitted phase mobility remain explicit
caller responsibilities.

Physical molecular heat is added once to authored NET caloric increments.
Original phase states remain the initial state until the thermal owner advances
h and phase and derives compatible expansion/contraction flux. Viscously
updated momentum then receives the accepted pressure impulse and follows the
same accepted mass routes as M and M*h. No source is sampled again during
thermal trials, and no second viscosity solve or momentum pressure solve is
introduced.

Heat-addition residuals share the original thermal enthalpy budget and the
viscosity energy gate. Wall work, required source-pressure work, temporal loss
and scatter/transport losses retain separate accounts; numerical dissipation
does not become heat. A sealed geometry incompatible with its chart's net
expansion refuses instead of inventing compression. The existing closed
wall/symmetry and opposed pressure-mouth restrictions remain. This is a
first-order frozen-source split with caller-owned thermodynamic work and phase
support, not a pressure-dependent EOS or deformable-solid closure. All stages
must accept before any endpoint is returned.

## Implicit original-source heat driving thermal volume and momentum

With `implicit-enthalpy-coupling`,
`advance_isotropic_implicit_conductive_thermal_volume_with_baths` and its
momentum counterpart solve one backward-Euler heat interval before the existing
phase-density/pressure/material transport. The caller supplies bounded native
Newton/Krylov work and a reusable workspace sized to the actual fluid cells.
The original flow duration is the heat interval; there are no hidden thermal
substeps or retries, and no explicit conductive mixing-weight restriction.

The adapter uses the existing `ImplicitConductionGraph` and the SAME
`PreparedPhaseConductivity`/borrowed caloric objects. Actual mass and graded
geometry are held during this heat interval. T(h), k(h), dT/dh and dk/dh are
evaluated at each nonlinear trial; latent storage retains zero dT/dh and its
actual phase-dependent conductivity. An admitted zero conductivity must retain
the source's zero selected tangent. Native sparse analytic linearization,
line search, source-domain checks and work limits decide acceptance.

Every selected exterior or embedded-wall face retains its actual fluid half
distance, area, bath temperature and resistance. Multiple baths on one cell
contribute separate Fourier powers and their complete analytic derivatives.
They are attached as local native feedback powers without blending bath
temperatures or adding thermal inventory. Their native entropy contribution is
explicitly unclassified; the native report is not a complete fluid/bath entropy
audit. Original authored heat is additive over the same interval.

Accepted endpoint internal/bath heat is independently reassembled and compared
with each cell's native implicit storage change. Absolute local mismatches,
heat-addition defects, boundary balance and the original transport defect share
the existing total enthalpy budget. Only the accepted physical heat enters
native volume/momentum once; pressure work and numerical losses keep their
existing accounts. Fixed solids retain zero material inventory. The optional
mechanical pressure profile keeps its original owner and actual mass routes.
This remains a first-order conduction/flow split requiring time refinement.

Focused definitions cover stiff pair/cycle/cube backward-Euler limits, decimal
graded face coverage, nonlinear k(T) and k(f) latent storage, opposed baths on
an embedded wall, independent caloric/momentum transfer and transactional
source/work/cancellation refusal. Explicit APIs retain their original behavior.

### Molecular heat in the same implicit thermal-volume endpoint

`advance_isotropic_implicit_conductive_thermal_volume_momentum_with_viscosity_and_baths`
prepares the existing molecular impulse and physical cell heat once. Its
original addition/work gates combine authored heat with that molecular heat
before the shared implicit conduction/bath solve. The accepted caloric endpoint
therefore includes molecular heating when evaluating k(T,f) and every bath
flux. Native volume/momentum then uses the already updated cell momentum;
neither viscosity nor its heat is applied again. Original viscosity and phase
sources remain frozen at the original state, without a simultaneous
thermal-viscosity constitutive iteration. Molecular numerical losses and wall
work retain their existing separate accounts. The original viscosity report,
heat-addition residual and total thermal budget remain available. The focused
one-sided-bath regression checks an independent two-cell backward-Euler balance
and the original momentum evolution, including one sample per molecular source.

## Prescribed thermal baths driving phase-volume heat and momentum

`advance_isotropic_conductive_thermal_volume_with_baths` and
`advance_isotropic_conductive_thermal_volume_momentum_with_baths` attach explicit
`ThermalVolumeBathFace` laws to the original conductive volume timeline. Each
unique fluid boundary face supplies its bath temperature and nonnegative
area-specific resistance. Faces at the grid exterior and adjoining fixed-solid
voxels are admitted; fluid/fluid interfaces are not bath boundaries. The
original graded grid supplies actual face area and the fluid half-cell distance;
the same original k(T,f) sample used by internal conduction gives `G=A/(d/k+R)`.
An embedded solid defines the boundary geometry without contributing its width,
conductivity, heat capacity or material inventory. The caller explicitly owns
the prescribed ideal wall/bath temperature and additional resistance. Unlisted
faces retain insulation. The same bath input is
available alongside the original molecular-viscosity and phase-resistance
momentum owners, including their existing mechanical pressure profiles. Each
owner keeps one bath entry point and delegates legacy insulation to empty baths.

Bath heat enters the existing net-enthalpy/phase-density/pressure/momentum
advance exactly once. Signed heat by face normal and gross absolute transfers remain
separate from internal exchange, advective enthalpy and required mechanical
source work. The original local heat-addition defects and total enthalpy budget
still gate publication. Sorted physical cell/face order removes input-list
ordering from the floating-point result. Combined internal/bath face work is
bounded by the existing conductive face limit.

Bath weights use the original piecewise caloric chart's maximum dT/dh over
the admitted source range, jointly with the original internal enthalpy secants.
Sharp latent storage is not replaced by an apparent heat capacity; an entirely
latent source interval has zero thermal slope and still enforces actual
enthalpy-domain admission. Coefficients are held at the original cell state
for this explicit interval. This is first-order split heat transport, not a
nonlinear through-wall solve, evolving finite bath, or closed pressure-work
model. Focused G1/G3/G4 regression definitions cover six graded physical faces,
source coefficients and caloric datum shifts, embedded-wall calorimetry with
zero solid inventory, opposed bath heat, latent energy,
native material momentum, additive molecular heat, work/stability refusal and
cancellation/replay.

## Source conductivity driving thermal-density momentum

`FlowField::advance_isotropic_conductive_thermal_volume_momentum` and its
pressure-profile variant carry original-source heat through the same accepted
thermal-volume and material-momentum step. They reuse the existing conductivity
owner and its borrowed original caloric chart; source identity, density,
temperature and phase share that chart throughout the endpoint.

Original-state source samples and harmonic physical half-cell resistances
produce paired face heat once. The original enthalpy-secant monotonicity gate
also covers latent plateaus without an invented heat capacity. Those increments
are added once to authored net heat, and the same shared enthalpy-defect budget
covers assembly, caloric representation and accepted native transport. The
accepted pressure impulse reaches cells once before their original mass routes
carry momentum; thermal trials do not repeat conduction or body forcing.

The source is explicitly isotropic; exterior and fixed-solid thermal faces
remain insulated. Source-pressure work and numerical mechanical losses are
reported separately and never deposited as heat. Existing first-order split
and refinement requirements, phase mobility and caller-owned thermodynamic
work convention remain. Source, heat, chart, volume or momentum refusal exposes
no partial candidate. The original conduction-only API retains its behavior
through the same private heat-preparation and acceptance runner.

Independent analytic regression sources cover sensible-expansion heat transfer
and exact cell impulse on all axes and enthalpy datums, latent heat without
spurious plateau conduction, and bidirectional reservoir transport under large
ambient references with missing-inflow-velocity refusal and exact retry.
These are compiler-checked regression sources; they have not been executed here.

## Conductive and molecular heat on one evolving density chart

`FlowField::advance_isotropic_conductive_thermal_volume_momentum_with_viscosity`
and its pressure-profile variant combine original source conductivity with
molecular stress and its reciprocal physical heat before the accepted
thermal-density/momentum evolution. Both source evaluations use the same
original phase states. Conductive face heat, the viscosity solve and body-force
prediction each run once; bounded phase-volume trials consume their fixed
caloric increments and viscously updated momentum.

The original conduction preparation and final source-domain checks are shared
with the conduction-only and conductive-momentum owners. The actual caloric
state after adding molecular heat must also stay inside the original source's
support, even if later advection would bring the endpoint back inside. Its
enthalpy gate
includes conductive assembly, molecular-heat assembly and native caloric/
transport defects together, preventing two individually accepted heat stages
from exceeding the common budget. Source-pressure work, reservoir/wall work
and numerical losses retain their native separate accounts.

The original isotropic-conductivity and insulated exterior/fixed-solid thermal
boundary assumptions remain, as do the viscosity boundary restrictions,
first-order frozen-source splitting and caller-owned phase/work conventions.
A source, heat, chart, volume, momentum or cancellation refusal returns no
partial heated or accelerated state.

## Evolving material momentum on one native pressure and transport step

`FlowField::advance_material_momentum` composes the existing current-mass
reconstruction, native body-force/pressure impulse feedback and conservative
material transport. Each original injection/extraction rate supplies its own
signed volume target. The same pressure-corrected field produces the accepted
mass routes, and those immutable routes transport the mechanically updated
extensive momentum with zero additional impulse. There is one pressure solve;
momentum creates neither another mass solve nor an independent material flux.

The returned cell mass, specific enthalpy and extensive momentum are the
endpoint inventory; extensive enthalpy is mass times specific enthalpy. The
retained `transport.projection.field` is explicitly the force/project transport
stage, not a reconstruction of endpoint P/M. Every subsequent call reconstructs
interior/free faces from the supplied current mass and momentum. The receiver
supplies held boundary fluxes only, so stale transport velocities cannot replace
the authoritative endpoint momentum.

Opposed volume sources retain their gross mass, enthalpy and momentum exchanges
even when their net volume is zero. Every incoming material route requires its
authored velocity; extraction and boundary outflow use the actual donor.
Mechanics and material must share identical original projection controls,
including duration and acceptance gates. Original work limits, local/global
balance checks and cancellation remain in force. No stage returns a candidate
until mass, enthalpy, momentum and every original stage gate have accepted; all
input inventories remain unchanged on refusal.

This is a first-order force/pressure/advection split. It adds no viscosity,
thermal expansion, EOS or phase mobility and does not claim a complete
Navier-Stokes integrator. Native numerical losses remain separate from physical
enthalpy. Three authored regression sources cover independently calculated
opposed-source inventories and one pressure impulse, repeated evolution from
endpoint momentum with deliberately unrelated receiver velocities, and
mismatched clocks/gates, missing inflow velocity, native budgets, cancellation
and exact retry. Independent arithmetic checks were evaluated separately;
Rust compilation, rustfmt and native test execution were unavailable.

`advance_material_momentum_with_pressure_profile` extends that same step to
authored reservoir heads on native open faces. The fresh current-density
pressure impulse drives one accepted material transport, retaining its
original mass, enthalpy and momentum routes. Incoming/backflow material and
velocity remain explicit. Mechanical correction pressure does not become a
thermodynamic EOS, pressure heat or an extra enthalpy source; its reservoir
work retains the native gauge and the enclosing owner's flow-work convention.
No stage returns partial inventory on profile, source, solve or budget refusal.
Two regression sources independently derive unequal-density pressure-driven
inventories and hydrostatic gravity cancellation, including large pressure
references, exact empty-profile parity and late-failure retry.

### Uniform reservoirs through the existing material pressure owner

`ReservoirPressure` adds six optional uniform plane pressures in Pa alongside
the existing borrowed `StaggeredBoundaryPressure` profile. The generalized
`force_and_project_density_with_reservoir_pressure`,
`force_project_material_momentum_with_reservoir_pressure` and
`advance_material_momentum_with_reservoir_pressure` entries accept either input.
Uniform reservoirs need no staggered pressure buffers or profile allocation.
Each supplied plane value must be finite and its face role must be `Free`;
`None` retains zero pressure. Existing profile and zero-pressure APIs remain.

Both representations feed the existing native projection, `PressureReadback`,
force/pressure feedback and accepted material routes. Component-relative
impulse/work arithmetic, the authored public pressure gauge and original
transactional admission are retained. No second pressure map, solve or impulse
owner is introduced. Reservoir work remains in the existing mechanical accounts;
material transport still carries exactly `M*h`, with explicit incoming material
and velocity. Pressure work is not added again to enthalpy, and no EOS or new
total-energy closure is inferred from the reservoir input.

Additional native material regressions check analytic pressure-driven flow and
reversal, actual mass/enthalpy/momentum transfer, a graded density jump with
bit-identical physical results under a `1e12 Pa` common datum, and hydrostatic
rest at a stratified side opening. Uniform plane and published spatial-profile
entries give identical complete material results. Invalid plane data, missing
inflow state, invalid profiles, cancellation and retry retain refusal semantics.
All 64 focused native tests passed, including the previously published profile
tests; this is not a full workspace build.

### Endpoint-sampled D3Q19 force histories

`BoundaryGrid3::step_with_endpoint_force_scale` and its pooled counterpart
advance the original collision/stream operator with the force at `t_n`, then
publish the force at `t_(n+1)` with the accepted populations. The next scalar
multiplies the field's original selected spatial shape. Physical momentum
`sum(f_q c_q) + F/2` consequently receives `(F_n + F_(n+1))/2` per lattice tick;
the unchanged uniform force is included once. Periodic transport cancels in
the global balance; wall reaction and open-boundary transport remain separate.
This is second-order trapezoidal integration of a prescribed smooth force
history, not an implicit endpoint constitutive solve. Initialize its first
force with `rebind_force_field` to preserve the initial physical velocity.
Invalid endpoint data, a collision refusal or cancellation preserve both
accepted populations and force. The serial and pooled paths reuse the same
cell kernels and publication scratch; worker count preserves population bits.

The forced-jet example now uses physical force amplitude `omega * v_ac`,
momentum-preserving initial rebinding, and the pooled endpoint API. Its smooth
sinusoidal response has the declared second-order time discretization error.
`boundary_forcing` tests cover spatial and uniform forcing, all four forcing-
compatible collision models, physical-momentum closure, sinusoidal convergence,
the incorrect-order negative control, serial/pooled bits and atomic refusals
(G0/G1/G3/G4/G5). No new turbulence or open-boundary accuracy claim is made.

### Actual-density body acceleration on the dense lattice

`BoundaryGrid3::rebind_body_acceleration` uses the original zero-motion Guo
preparation to select `F = rho*a` without elapsed impulse or a stress reset.
Its private open-boundary admission does not broaden `RigidFrame3`'s public
rotating-frame restrictions.
`step_with_endpoint_body_acceleration` and its pooled/selected-wall variants
then use the original collision, streaming and wall-load kernels. Collision
retains `F_n`; the streamed density supplies `F_(n+1) = rho_(n+1)*a_(n+1)`.
There is no endpoint population reframe, so physical momentum receives the
native trapezoidal source impulse `(F_n+F_(n+1))/2`. Passing the same acceleration
at successive endpoints gives a held physical acceleration at every density.
An ordinary `step` still holds force density, as its existing contract states.

The endpoint force/velocity admission runs on the populations already held by
each streaming kernel. Publication reuses the accepted force allocation after
first use; it needs no new grid clone, source solver or second collision.
Cancellation, nonrepresentable endpoint forces and late collision/stream refusal
leave accepted populations and sources unchanged. Closed walls, periodic
domains and native velocity/pressure faces are supported. Open reconstruction
uses the interior's actual endpoint density and force-corrected stress. A full
face retains the regularized template's density. At a wall/open seam, only
open-owned rows change: the native Guo source gives an affine mass law
`rho = B + C*rho`, solved directly as `B/(1-C)` before reconstructing in the
`F = rho*a` source frame. This avoids iterative density guesses and retains
wall-owned populations. Singular, nonpositive or nonrepresentable closures
refuse atomically; the actual final population sum is checked before publication.

Spatial velocity/pressure profiles remain authoritative. As with the original
masked reconstruction, nominal full-face density/velocity targets need not hold
at a mixed wall seam. Reservoir mass and momentum transport stay separate from
body-force impulse; neither sealed-domain conservation nor reservoir work is
claimed. Immersed finite-body feedback and full-flow energy conservation still
need their own coupled owners.

`tests/dense_body_acceleration.rs` checks unequal-density response and changing
acceleration through four collision laws, retained non-equilibrium stress,
independent affine-density endpoint moments, sinusoidal second-order time
convergence, periodic total/local momentum, wall linear/angular impulse,
serial/pooled population bits and transactional refusal (G0/G1/G3/G4/G5).
`tests/open_body_acceleration.rs` adds all-axis/full-face analytical response,
wall-seam closure against independently iterated prescribed-force reference
steps, reservoir mass/momentum accounting, pooled parity and failed open-source
publication. The scalar closure's inline tests cover retained rows, independent
fixed-point agreement, noncontractive cases, density scaling and refused roots.

### Shared momentum reduction for immersed finite bodies

Single- and multiple-body partial-saturation coupling share one native spatial
response assembly. Fixed-size double-double accumulators retain weak net forces
and couples under opposing cell tractions before the original f64 Cholesky solve.
A shared reaction reduction preserves actual fluid/body impulses, moments about
the requested origin and body center, and work; moving, dry-body and lubrication
owners consume the same corrected path. Collision and streaming arithmetic is
unchanged, with no per-cell heap allocation or additional fluid solve.

This does not manufacture sub-ULP population feedback: the measured momentum
residual remains exposed. `partial_saturation_cancellation3` checks the actual
single-/multiple-body response to analytically balanced incoming fluid loads;
the internal `momentum_reduction_tests` checks a weak reaction couple and work.

## Prescribed three-dimensional rigid-frame flow

`d3q19::RigidFrame3::prepare_step` supplies translating and arbitrarily rotating
co-moving `BoundaryGrid3` domains with the full lattice acceleration
`g - a_origin - alpha cross r - omega cross (omega cross r) - 2 omega cross u`.
Native cell centres are `[x+0.5, y+0.5, z+0.5]`, matching geometry and wall torque
positions. The Coriolis source uses an analytic 3-by-3 midpoint inverse; its
isolated source step preserves speed and axial velocity without clipping or
iteration. Physical body acceleration and relative wall/fluid velocities are
already expressed in the moving frame.

The owner replaces uniform and per-cell forces together. The original Guo
half-source transform preserves current mass, physical momentum and corrected
stress, so preparing or switching frame loads introduces no impulse. Existing
collision scratch and per-cell force storage are reused after first allocation;
late admission failure or cancellation before publication changes neither
accepted populations nor forces. Call the original step after preparation and
refresh the frame load each step. BGK, TRT, central-moment and reduced-cumulant
collision laws retain their original numerical kernels.

Closed finite domains are supported. A periodic direction requires both spin
and angular acceleration parallel to that direction; incompatible periodic
rotations and open pressure/velocity faces refuse. This prescribed-motion
owner does not solve vessel inertia, transform an inertial initial velocity,
or assert full-flow energy conservation or engineering validation.

Six focused native regressions exercise independent arbitrary-axis rotation
and Coriolis work, physical Euler/centrifugal signs, source-frame momentum and
stress preservation, compatible periodicity, transactional refusals, and the
actual closed-box collision/stream/wall linear and angular impulse account.
The runnable
`rotating_tank3` example supplies a smooth spin-up schedule to this same owner.

### Local Coriolis flow on periodic domains

`core2::frame::Coriolis2` and `d3q19::Coriolis3` select the local source
`a_effective - 2*omega cross u` explicitly. Centrifugal acceleration is already
absorbed into the authored effective potential/gravity; Euler and origin
acceleration are not inferred. The two-dimensional f parameter is twice the
out-of-plane angular velocity. Three-dimensional angular velocity may have
any direction, including the horizontal component of planetary rotation.

The original two-/three-dimensional midpoint inverse supplies the source,
with the same physical Guo reframe and atomic field publication. The local
source is translation invariant and therefore admits arbitrary periodic axes;
the full `RigidFrame3` owner retains its incompatible-periodicity refusal.
Open density/velocity faces remain unsupported. Constant specialization removes
the unused position-dependent rigid-frame terms without adding an iterative
solve, per-cell trigonometry or a second dense scratch allocation.

Refresh before every native step. A homogeneous force-free inertial oscillation
preserves kinetic energy and axial velocity to roundoff; its per-step angle is
`-2*atan(|omega|)`, so phase is second order rather than exact. The homogeneous
geostrophic state `a_effective = 2*omega cross u` is stationary. This does not
claim exact energy conservation for transport/walls, spherical geometry,
beta-plane variation, or resolved atmospheric thermodynamics. Parameters use
the existing lattice units, with `omega_lattice = omega_SI * dt`.

`tests/coriolis.rs` covers both spin signs, zero spin, arbitrary-axis periodic
collision/streaming in all four admitted D3Q19 collision models, independent
Rodrigues phase, geostrophic balance, current-velocity preservation, full-frame
periodicity refusal, closed/open faces, invalid cells and cancellation.

### Local Coriolis on sparse and material liquid trajectories

`Coriolis3::prepare_free_surface_step` uses that same source on active sparse
wet cells, including newly created interfaces on the next refresh. It retains
canonical density, tracked mass/fill and the original Guo velocity/stress
representation. Historical replay and independent body-acceleration ownership
keep the full-frame owner's explicit refusal. Source specialization shares
one sparse traversal/rebind with the existing rigid-frame path.

`MaterialCoriolis3` and `MaterialFreeSurface3::step_with_coriolis_si` convert
spin by the original dt and effective acceleration by dt²/dx. Unrepresentable
nonzero conversions refuse. Source preparation and the native timestep share
one private candidate and atomic publication; source identity, pressure,
capillarity, topology and mass stay with the original material owner.

`FreeSurfaceRheology3::step_with_coriolis` and
`MaterialFreeSurfaceRheology3::step_with_coriolis_si` apply the same force while
retaining the original local relaxation and endpoint constitutive admission.
Coriolis does not become shear heating or a new constitutive source. These
isothermal entries add no thermal coupling or finite-inertia vessel dynamics.
The native force/mass/moment, new-interface, SI conversion, constitutive and
cancellation regression definitions accompany the implementations.

### Local rotating thermal and Boussinesq liquid

`MaterialFreeSurfaceEnthalpy3::step_with_coriolis_si` and
`step_with_coriolis_and_viscous_heating_si` compose the local rotating source
with the existing single thermal/VOF transaction. Material density departure
weights authored effective gravity only; Coriolis uses reference inertia,
retained nongravitational loads remain force densities, and Coriolis supplies
no thermal source. SI conversion, midpoint momentum, heated source viscosity,
native mass routes, extensive enthalpy transport and endpoint/clock acceptance
share one sparse candidate. The accepted endpoint restores the original
thermal force datum without another impulse; rotation and gravity are supplied
on each interval.

Full rigid-frame rotating-Boussinesq refusal remains, as do first-order thermal
splitting and complete-boundary ownership limits. This adds neither
density-dependent accelerating-vessel inertia nor a total mechanical-energy
theorem. Independent dimensional midpoint, material feedback, no-heat,
interface-transport and atomic-refusal regression definitions are in
`tests/material_free_surface_enthalpy/frame.rs`.

### Persistent local Coriolis in dense thermal and material slabs

`ThermalLbm::set_coriolis` and `ThermalLbm3::set_coriolis` retain admitted
local rotation for ordinary, external-force, prescribed-heat, transport and
rheology steps. Their `MaterialSlab::set_coriolis_si` methods convert physical
angular velocity with the slab's retained timestep,
`Omega_lattice = Omega_SI*dt`. Invalid or unrepresentable configuration leaves
the setting unchanged; setting or clearing rotation changes no populations,
physical velocity, stress or thermal history. Zero disables the source on the
next step.

The common force-binding seam composes rotation with the complete selected
non-Coriolis load immediately before the single native momentum collision.
Opening, endpoint-average and original-material buoyancy keep their existing
sampling epochs. The analytic midpoint inverse is admitted once per setting;
no grid-sized rotation field is allocated. The two-dimensional grid's separate
uniform acceleration participates in the same solve and retains its original
owner. Guo reframing preserves physical velocity and stress at force refresh,
without an extra half impulse. Material density departures weight effective
gravity only, never the Coriolis term.

Checked source/material snapshots retain the setting and existing rollback.
Wall exchange, heat accounting, material query ownership and constitutive
stress reconstruction stay on the original trajectory. Coriolis supplies no
molecular or sensible heat. Isolated homogeneous source motion preserves speed
and axial velocity with the midpoint phase; spatial transport, thermal
splitting, boundaries and their stability limits remain unchanged. Fluid and
wall velocities are relative to the local axes, with background centrifugal
effects represented in effective gravity.

Regressions in `thermal_sources.rs`, `material_buoyancy_endpoint.rs` and
`material_thermal_rheology.rs` exercise signed/arbitrary-axis phase, both
buoyancy quadratures and native collision families, independent combined-force
oracles, material SI scaling and source epochs, refusal/cancellation replay,
and unchanged initial molecular strain and heat. This extends the dense thermal
owners; the sparse liquid rotating-source contracts above remain separate.

### Acceleration of sparse free surfaces and material-owned vessels

`RigidFrame3::prepare_free_surface_step` applies the SAME admitted rigid-frame
midpoint source to original sparse `FreeSurface3` wet cells. The field includes
physical body acceleration minus origin acceleration, Euler and centrifugal
terms, and the energy-neutral local Coriolis solve. All vectors use moving
axes, velocities are relative, and centres remain `[x+0.5,y+0.5,z+0.5]`.
Current local density converts acceleration to lattice force density. The
complete force field replaces initial gravity; pressure, capillarity and wall
traction are not counted twice. Historical `LegacyGlobal` replay refuses.
An independently selected `rebind_body_acceleration` owner also refuses,
including an owned zero acceleration: its density-dependent endpoint source
law cannot silently become a held midpoint force field. The frozen material
constructor retains its original reference-force convention and does not
implicitly select that independent owner.

The existing local Guo rebind preserves canonical density exactly and physical
velocity/force-corrected stress to population roundoff, without elapsed time,
impulse, mass, fill or topology changes. Storage follows active wet cells,
including zero-mass interface ghosts, not the full domain. Native conversions
inherit the neighbor force frame; callers refresh the prescribed source each
step so newly wet cells use their actual coordinates and accepted velocities.
Cancellation and invalid cells refuse before publication.

`MaterialFreeSurface3::step_in_rigid_frame_si` consumes `MaterialRigidFrame3`
and body acceleration in SI using the original liquid's dx/dt. Origin scales
by 1/dx, angular rate by dt, angular acceleration by dt² and linear acceleration
by dt²/dx. Nonfinite, overflowing or disappearing nonzero conversions refuse.
One private candidate joins preparation and the original native step, so
failure/cancellation preserves accepted forces, populations, mass and topology.
The caller supplies midpoint motion and owns its external clock. This frozen
liquid entry does not silently combine thermal buoyancy or solve a finite-inertia
vessel/gas, and no exact energy theorem is asserted for the full VOF step.

Focused G0/G1/G3/G4 regressions in `rotating_frame3/free_surface.rs` and
`material_free_surface3/forcing.rs` exercise independent Coriolis rotation and
work, density-independent translational acceleration/free fall, physical
Euler/centrifugal signs, exact mass/fill and physical moment retention, real
wet-cell conversion, SI unit/source identity, changing load histories and
transactional/independent-owner refusal. Original dense-grid rigid-frame
arithmetic is retained.

`FreeSurfaceRheology3::step_in_rigid_frame` and
`MaterialFreeSurfaceRheology3::step_in_rigid_frame_si` advance that SAME source
law with the original analytic or source-resolved non-Newtonian viscosity.
The native uniform/spatial assignment and material card/claim identities remain
attached. The SI entry shares the frozen-flow unit conversion. Its actual
pressure-dependent source still queries interior liquid pressure at the fixed
construction datum, not an externally prescribed gas value.

The Guo preparation preserves physical strain, so accepted local relaxation
enters the same forced collision without a redundant constitutive solve.
Actual endpoint force-corrected stress and pressure then refresh viscosity,
including new wet cells, before any state publishes. One private sparse
candidate contains source preparation, native collision/transport/conversion
and endpoint admission; no unused enthalpy routes or extra flow solve are
allocated. Refusal preserves force history, populations, mass, local viscosity
and source samples. Independent body-acceleration owners and historical replay
retain their explicit refusal. Motion coefficients are refreshed every step;
heat, finite-inertia vessel dynamics and simultaneous prescribed gas loading
remain outside these isothermal entries.

Native rheology tests and `material_free_surface_rheology/frame.rs` specify
independent incoming BGK/Guo populations with local relaxation, the midpoint
SI source equation, a changing source-viscosity trajectory against its closed
quadratic constitutive root, free fall and uniform acceleration without false
shear, and atomic source/convergence/ownership/cancellation refusal.

`MaterialFreeSurfaceEnthalpy3::step_in_rigid_frame_si` now joins prescribed SI
frame motion and authored heat to that same source-owned VOF trajectory.
The original thermal split resolves heated viscosity and optional rheology
before one native collision; the actual mass/conversion routes carry extensive
enthalpy once. The existing frame kernel admits locally resolved acceleration,
so no second frame integrator, constitutive query or sparse fluid candidate is
introduced. Source pressure, caloric mass and collision scaling remain owned by
the original material binding.

Constant-density thermal/rheological liquid admits arbitrary-axis rotation.
Bounded Boussinesq liquid admits translation: the original full gravity plus
density-departure load is removed exactly, and its departure multiplies
`g - a_origin`. The frame kernel subtracts origin acceleration once, yielding
`(1 + density_departure) * (g - a_origin)` before numerical-density conversion.
Already-resolved non-gravitational forces retain their force-density convention.
Rotating Boussinesq liquid refuses either nonzero angular velocity or nonzero
angular acceleration, including startup from rest: this reference-inertia
model does not own the additional density-dependent inertial couplings.

At acceptance, the original thermal Guo force representation is restored
without a second impulse or strain change; callers supply frame motion again
on every driven interval. `step_in_rigid_frame_with_viscous_heating_si` deposits
only the original accepted molecular strain power. Frame work, Coriolis work,
pressure work and numerical kinetic loss do not become heat. Electrical and
surface-heat bindings still require their complete original boundary timelines.
The existing first-order thermal/upwind split, low-Mach limits and prescribed
vessel/gas assumptions remain. No total mechanical-plus-thermal energy theorem
is asserted. All material, transport, mass, enthalpy, speed and cancellation
gates precede joint force/population/source/enthalpy/clock publication.

Four focused definitions in `material_free_surface_enthalpy/frame.rs` check an
independent SI midpoint solve and all incoming heated BGK/Guo populations,
repeated free fall without spurious molecular heat, changed sourced density and
effective gravity, exact endpoint force restoration, actual interface heat
transport, and complete rollback including a late native mass-route refusal.

### Finite-inertia planar missing-film surfaces

`FiniteLubricationPlane3` and
`BoundaryGrid3::step_with_lubricated_partial_saturation_finite_planes`
bind the existing normal and optional leading-logarithmic tangential film
laws to finite mass and a full world-frame inertia tensor. Plane mobility
joins the original simultaneous scalar resistance solve. Every sphere film
sharing a plane therefore sees that plane's solved recoil and spin during
the same native fluid tick.

For a sphere mode `(n, a)`, the opposing plane wrench uses
`(n, (sphere_center - plane_center) cross n + a)` with negative sign.
The local arm retains total world angular momentum and makes common rigid
translation/rotation a null mode. Plane film work uses the solved center
and spin midpoint velocities. It enters the plane's kinetic state;
authored external impulses and couples remain separate. Original film
dissipation has no inferred thermal destination.

Returned finite-plane responses contain accepted velocities, directly
accumulated center-angular impulse, kinetic change, external work and
momentum/work residuals. These residuals measure numerical balance; they
are not additional numerical-tolerance acceptance gates. Finiteness,
positive-definite inertia, native low-Mach admission, coupled footpoint
speeds, existing distinct-plane/positive-gap rules and film/mode budgets
must all pass before any sphere, plane or population state publishes.
Even with no active film, finite planes receive their authored external
impulses. Cancellation and refused candidates preserve accepted state.

This is a frozen-geometry exchange, like the existing spatial immersed
body response. Centers, plane geometry and world inertia remain fixed.
Plane extent and mass are caller-owned; no resolved moving wall, solid
contact, pose drift, gyroscopic step, displaced-fluid closure or GCL is
provided. Existing stationary resolved-wall reactions remain separate
external ports. Regression definitions cover the independent two-mass
normal response, shared-plane angular/work balances and origin shifts,
common rigid rotation, forced sampled-body BGK/TRT coupling and atomic
refusals. Native compilation and test execution remain outstanding.
No measured physical-validation or performance claim is added.

### Prescribed axial spin in the original sphere/plane film solve

`LubricationPlane3::with_axial_spin` adds constant rotation about the plane
normal through its translating `point`. Each active film uses
`V + Omega cross (footpoint-point)` at the original exchange pose; the moving
body owner supplies its original midpoint plane and midpoint sphere geometry.
The same cached local velocity drives both the native simultaneous Schur solve
and the accepted-slip dissipation. There is no additional body kick, fluid
collision, contact projection, or clock. Existing normal/tangent cutoff,
roughness, resistance-mode and positive-gap restrictions remain authoritative.

Plane reaction work includes translation and the complete moment about the
spin axis, including the existing free couple. Rotational work is accumulated
from local lever arms with compensated summation, without subtracting large
world-origin torques. Reports retain equal/opposite linear and total angular
impulse, body work, nonnegative modeled loss, and the work residual. Zero spin
retains the original arithmetic. Nonfinite or high-Mach active-footpoint
velocity refuses before any body/fluid publication; no global Mach bound is
claimed for an infinite rotating plane. Resolved-fluid boundaries remain a
separate required owner and must share the prescribed surface motion.

This extends only the original divergent squeeze and tangent/parallel-spin
resistance. It does not add finite axial torsional resistance: a stationary
sphere on the axis receives no added axial torque from this correction.
The finite small-gap torsional limit is discussed by Liu and Prosperetti,
JFM (2010), section 4, equation 4.1,
<https://pages.jh.edu/aprospe1/publications/PapersPublished/Rotating/LiuJFM_2010.pdf>.
General tilting planes, solid contact and complete resolved rotating-wall flow
remain outside this correction. The separate finite-inertia plane owner keeps
its original solved recoil, full-wrench work and Schur mobility; a finite
carrier refuses an additional nonzero prescribed spin. Zero remains an identity.

Six selected-source native Rust tests execute the original resistance solve,
six-velocity body response, Cholesky, cancellation gate and deterministic math.
They cover local rigid motion, rotation covariance, rolling-mode resistance,
opposing footpoints and actuator torque, origin-independent work, translation/
rotation cancellation, explicit absent axial drag, zero-spin continuation and
local-speed refusal. Two public sampled-body integration tests additionally
cover midpoint kinematics and atomic failure/retry; these source tests parsed
but were not executed in the delivery environment. This is not a full crate
build, resolved-wall validation or experimental calibration.

### Evolving material momentum with sourced molecular viscosity

`FlowField::advance_material_momentum_with_viscosity` borrows one original
`FluidTransport` evaluation per mobile cell, with no source in fixed solids.
Only its frozen dynamic viscosity is consumed; actual evolving M/V supplies
inertia through the original density projection's half-cell face masses.
The original full symmetric MAC rows and Krylov solver advance a
backward-Euler molecular predictor without convection or pressure. Reciprocal
shear, harmonic material-jump transmission, tangentially moving walls and
stationary no-slip solids retain their original operators. Existing time-level
paths keep their original constant-density arithmetic.

Native impulses J_f = M_f*(u_new-u_old) reach original extensive cell momentum
through half-cell mass shares. Independent native and cell accounts check
delta_K_face + BE_loss + Q_molecular - W_wall and
delta_K_cell + scatter_loss - delta_K_face - representation_work.
BE loss is sum M_f*delta_u_f^2/2; scatter loss is
sum M_c*(delta_u_left-delta_u_right)^2/8. Physical unscaled frozen-row residual,
actual linear/angular impulse transfer, absolute per-face representation work,
molecular heat partition and actual M*(h_heated-h_original) insertion have
physical-unit acceptance gates. The native active-face impulse is not a complete
fixed-wall reaction measurement.

Only duration times the original molecular cell power enters specific enthalpy.
Backward-Euler, scatter, projection and material-transport losses remain
numerical accounts. The existing force/pressure/accepted-mass-route advance then
runs once with the updated momentum and heated enthalpy. Returned flow remains
the transport-stage field; endpoint state is the actual returned M, h and P.
A later call reconstructs current momentum and consumes newly evaluated sources.

Scope is closed, fixed geometry and a first-order frozen-viscosity split.
The original entry retains `2 mu D`, with no independently supplied bulk law.
There is no added EOS, phase mobility, turbulence,
pressure-opening work model or implicit viscosity-temperature feedback.
Cell and sweep ceilings bound work: each sweep has three existing momentum
solves, each capped at 20,000 Krylov iterations, and a frozen coupled residual
readback. Numerical storage is linear in original cells/faces. Any source,
solve, balance, heat-insertion or cancellation refusal returns no candidate.

`advance_material_momentum_with_bulk_viscosity` and
`advance_piecewise_thermal_volume_momentum_with_bulk_viscosity`, including their
pressure-profile variants, opt into `2 mu dev(D) + zeta div(u) I`. The caller
supplies one finite nonnegative bulk coefficient per mobile cell at the same
state as the original shear source; fixed-solid slots are ignored. Admission
uses the existing `NewtonianViscosity` law, including explicit zero bulk.
Normal rows use `4 mu/3 + zeta` and reciprocal `zeta - 2 mu/3` cross derivatives.
The matching physical heat is `2 mu dev(D):dev(D) + zeta div(u)^2`, entering the
same enthalpy or phase-density chart once. Existing shear-edge transmission,
actual M/V inertia, mirrored reservoir half cells, separate numerical losses
and atomic refusal gates remain unchanged. Both coefficients are frozen during
the original split; no EOS or bulk law is inferred from shear viscosity.

The rheology-resolved thermal-volume and implicit conductive/conjugate entries
also provide explicit bulk-viscosity variants. The rheology path observes the
original physical deviatoric rate before resolving shear from its source; bulk
is independently supplied at that same original state. The implicit thermal
paths add the predictor's matching shear-plus-bulk heat once to the existing
fluid, bath and finite-wall caloric graph. The molecular impulse is not applied
again after conduction. Both viscosities remain frozen over this first-order
split; implicit conduction does not imply implicit viscosity-temperature feedback.

The endpoint-coupled entry
`advance_isotropic_implicit_conjugate_thermoviscous_thermal_volume_momentum_with_baths`
adds bounded thermal feedback for both shear and independent bulk viscosity.
Every trial solves the original fluid and finite-solid caloric histories with
authored heat plus a guessed molecular heat vector. The original material source
is resolved at that exact native thermal endpoint, before advection; its shear
and bulk responses drive the existing native momentum kick using original M/P
and actual M/V inertia. Fluid states must remain fully liquid with no
unrepresented Darcy resistance; finite solids keep their own caloric owners.
The sum of absolute cellwise differences between heat used by the graph and
physical heat produced by momentum must fit the feedback tolerance and the
shared thermal and molecular energy gates. Accepted thermal transport retains
the graph's used heat, while the molecular report retains actual physical heat;
their signed and non-cancelling mismatch is explicit, never silently redeposited.
Only the accepted trial advances pressure and material transport. Rejected
trials consume bounded computation without accumulating physical history.
Thermal residual evaluations, thermal linear iterations and molecular sweeps
have whole-call limits, and unchanged coefficient fields reuse their original
momentum solve. This couples temperature-dependent viscosity within the existing
first-order split; it does not claim simultaneous endpoint shear rheology,
unconditional fixed-point convergence or an inferred bulk constitutive law.

The endpoint shear-rheology entry
`advance_piecewise_thermal_volume_momentum_with_endpoint_rheology` resolves the
original phase/material source on each native coupled momentum sweep. Its
physical RMS deviatoric rate is measured on the retained MAC velocity before
impulse scatter; source shear and independent bulk viscosities are refreshed
before all three nonlinear momentum residuals and the molecular heat audit.
Original time levels and actual M/V face inertia remain fixed throughout.
Acceptance requires both the relative row residual and the absolute physical
impulse residual, plus the existing angular-impulse and energy budgets. A small
coefficient change alone is not an acceptance criterion. Physical heat and wall
work use the same final coefficients; backward-Euler and scatter losses stay
separate. The original thermal phase is held for this constitutive solve, then
physical heat, pressure and conservative material transport advance once.
Retained endpoint rates, original source receipts and total source-update counts
make the operating point explicit. The original viscosity sweep limit bounds
the complete nonlinear solve, and source, support, work or cancellation refusal
publishes no endpoint. No rate cutoff, extrapolation or phase resistance is
introduced; this entry requires fully liquid mobile material.

The joint implicit entry
`advance_isotropic_implicit_conjugate_thermorheological_thermal_volume_momentum_with_baths`
closes both temperature and native shear-rate dependence. Each outer heat trial
samples its exact pretransport fluid phase once, then closes the native momentum
rows using shear and independent bulk responses at that phase and the current
native rate. Original M/P, actual inertia and fluid/finite-solid histories stay
fixed across both iterations. The same shared heat driver compares the inner
solve's physical molecular heat with the heat supplied to the implicit graph;
only acceptance permits one pressure/material advance. The temperature-only
entry retains its existing exact-coefficient cache, while the joint entry
performs a fresh nonlinear momentum solve at each thermal trial. Total source
queries, coefficient-field updates, native sweeps/Krylov work and thermal work
include the original predictor and rejected trials. The existing sweep and
thermal-work ceilings apply to the whole call. Returned rates, original material
receipts, used/actual heat and separate residuals describe the same accepted
constitutive endpoint. Fully liquid fluid support and the original first-order
pressure/advection split remain explicit; no hidden time steps or convergence
guarantee is introduced.

Four source-bound regressions cover graded actual inertia distinct from source
reference density, original liquid-viscosity changes, moving-wall work and
physical heat, current-state reconstruction, quiescent compatibility and atomic
refusals. All eight changed Rust sources passed pinned-nightly rustfmt parsing
through stdin; the three new files were formatted. Independent scalar response
and 101 unequal-mass scatter identities were evaluated, with the latter closing
within 2.23e-15 J. Native Rust typechecking, tests and runtime benchmarks remain
unexecuted because the shared filesystem is full and memory is constrained.

`advance_material_momentum_with_viscosity_and_pressure_profile` admits one
opposed pair of pressure reservoirs with the original normal-traction and
zero tangential-velocity-trace convention. Remaining planes retain closed
no-slip or symmetry conditions. `FvBoundary::Outlet` identifies the two
viscous mouths; the staggered pressure profile is their sole pressure owner.
Actual evolving M/V supplies both interior and boundary half-cell inertia.
Native mirrored mouth rows mirror that inertia by the same factor, and their
existing half-volume weight restores physical impulse, work and molecular heat.
One-sided mouth impulses reach their adjacent cell once. The viscous predictor
uses zero reservoir pressure; the following original pressure/material step
applies the authored profile once. Only molecular dissipation enters h.
Temporal, scatter, projection and advective losses remain numerical accounts.
The split remains first order with frozen source viscosity and fixed geometry.
Regression sources specify an independent unequal-density three-face diffusion
solution, pressure-driven viscous throughflow with a complete split energy
identity, and incompatible-pressure/trace/boundary refusal. General opening
configurations and nonzero tangential reservoir traces are not admitted.

## Rigid motion and actual link torque on a flooded free-surface wall

`FreeSurface3::planar_rigid_wall_motion_candidate` and the material owner's
`planar_rigid_wall_motion_candidate_si` admit two tangential translations and
normal-axis spin on one fixed, completely flooded rectangular domain plane.
The original halfway reflection evaluates the rigid velocity at each actual
link midpoint. The same population exchange supplies linear impulse, moment
about the declared origin, and local velocity dotted with impulse. The old
uniform-wall entry retains its original population arithmetic.

One link traversal also accumulates the symmetric 3-by-3 derivative of the two
tangential impulses and axial angular impulse with respect to the prescribed
twist. This affine response belongs to that unchanged input fluid and origin;
it is available to finite-inertia coupling without finite-difference trial
solves. Actual native replay still owns population and endpoint acceptance.
Accumulation needs constant auxiliary storage beyond the original candidate.

The SI adapter retains the original material source and scaling. Raw numerical
pressure, reference-pressure impulse and its moment, perturbation loads, and
optional absolute physical loads remain distinct. Work is mechanical exchange
at an isothermal boundary, not deposited heat or a liquid total-energy claim.
Corner speeds bound the finite plane; actual link-relative and both wet-cell
endpoint speeds retain their original low-Mach and material-scenario limits.
Geometry, mass, arithmetic and cancellation failures publish no candidate.
Normal motion, tilting, finite plate edges, cover/vacate and moving contact
lines are outside this fixed-plane operation.

Focused native regressions specify all-six-plane uniform-path parity, an
independent rotational impedance and liquid angular-momentum increment,
affine response replay, origin-shift covariance, and atomic refusal.

### Mixed prescribed and finite-inertia film surfaces

`LubricationSurface3::{Prescribed, Finite}` and
`step_with_lubricated_partial_saturation_surfaces` compose driven and finite
planar film ports in one original simultaneous sphere/fluid resistance solve.
Finite recoil changes the same shared film impulses received by prescribed
surfaces; no sequential wall kick or second fluid stream is introduced. One
bounded descriptor array preserves stable surface order. Original film/mode
caps, positive gaps, local speed admission and missing-resistance laws apply.

The external-load slice matches the surface count exactly. Each finite slot
owns one independent spatial carrier and accepts its center impulse/couple.
Prescribed slots require zero authored impulse: their known motion and received
film work already define the actuator port. Reports retain one reaction per
surface; `finite_body: Some(_)` identifies actual finite kinetic change and
external work, while `None` identifies prescribed reservoir work. Aggregate
`coupling.lubrication.plane_work` includes both and is not wholly external work.

All sphere, finite-carrier and population updates publish through the original
checked owner. The final finite assignments are infallible; refusal/cancellation
preserve accepted state. Prescribed descriptors are unchanged. Homogeneous and
empty inputs retain native arithmetic. Geometry and world inertia stay frozen;
there is no moving-wall geometry, dry-body trajectory, gyroscopic pose update,
solid contact, shared multi-face carrier identity or thermal destination added.

Four new native production-source G1/G3/G4 tests pass with six existing spin
regressions. Coverage includes an independent two-inertia velocity oracle with
fluid response, full-wrench finite recoil beside axial spin, frame covariance,
momentum/work, homogeneous/empty bit continuation and refusal/retry. Three
public BGK/TRT grid integration tests additionally cover actual population
momentum and two-tick state, native parity and joint rollback; these parsed but
were not executed in the delivery environment. No full workspace/grid pass or
experimental calibration is claimed.

### Immersed bodies with native velocity/pressure reservoirs

The single-body, simultaneous, moving and physical dry-body D3Q19 owners now
admit the grid's existing prescribed-velocity and prescribed-pressure faces.
The original boundary law still limits open faces to one grid axis, retains
paired periodic directions and stationary resolved walls, and reconstructs
physical Guo moments in the destination cell's own held force frame. All four
admitted collision models keep their original native fluid increment.

Open-face reconstruction runs on the unpublished streamed candidate, in the
same canonical face/cell order as an ordinary grid step. Positive immersed
fractions on reconstructed open-face cells are refused before collision.
Moving-body support still stays inside the existing domain halo; bodies
entering/leaving the reservoir, moving topology and contact are not introduced.
Zero occupancy recovers the ordinary native open step bit for bit.

`OpenBoundaryMomentumTransport3`, exposed as the shared `open` receipt, measures
mass, linear momentum and fixed-origin angular momentum delivered to the fluid
by the actual boundary operator. Its opposite belongs to the external reservoir.
During the existing collision traversal each open incoming link contributes
`2*c_q*f_post,opp`, the exact exchange from the original temporary reflection.
The subsequent native reconstruction contributes its actual population
differences, including replacement of tangential rows on open-face interiors.
Wall-owned corner links retain their original rows and wall account. Periodic
coordinate transport remains separately counted. The account is measured from
those local operators, not inferred from a whole-domain balance residual.

Reservoir reductions retain double-double low bits until the shared receipt is
rounded. After first-use native face-cache allocation, the added work scales
with open links and reconstructed boundary cells, with no additional grid-sized
field or whole-grid inventory scan. The single/assembly/dry owners reuse their
existing collision and next-population scratch. Failure during reconstruction,
receipt overflow, endpoint checks or cancellation restores the scratch/cache
owners without publishing populations, bodies or topology admission.

Physical dry-body combined linear and angular momentum gates include reservoir
exchange once. Nested body receipts do not duplicate shared wall, periodic,
forcing or reservoir terms. Full and retained-interior fluid masses remain
explicit inventories; this extension does not renormalize moving occupancy or
claim exact displaced-fluid mass/GCL closure. Reservoir mass is the measured
full-grid population transfer. This is an isothermal momentum extension, with
no pressure-work, thermal/internal-energy or complete fluid-energy claim, and
does not by itself establish calibrated continuum drag or high-Re validity.

`tests/open_partial_saturation3.rs` checks repeated exact native-step parity for
BGK/TRT/central-moment/reduced-cumulant empty overlays; independent forced
shared-body mass, linear and angular balances with wall/open and periodic
corners; translated-origin covariance; occupied-face refusal and exact retry;
late reservoir overflow rollback; cancellation after reconstruction; and actual
exterior-fluid-plus-dry-body inventories for a single body, delegated singleton
and simultaneous pair. The earlier unsupported-model regression now refuses
an unrepresentable reciprocal viscous time and moving walls, reflecting the
newly admitted moment collision and open boundary capabilities.

### Live reservoir targets without topology changes

`BoundaryGrid3::rebind_open_boundaries` replaces held velocity/density values
on existing open faces. Every face retains its boundary kind, and wall and
periodic data must stay unchanged. The original fallible face admission runs
before publication. A changed velocity target with an installed spatial
profile refuses until that profile is explicitly cleared or updated; an
unchanged inlet profile remains authoritative during outlet-pressure changes.

The update inspects six records, allocates no grid storage and preserves
compiled links, cached open cells, populations, held forces and their physical
Guo frame. It applies no instantaneous impulse and advances no time. The next
accepted native step applies the targets through the original regularized
reconstruction, whose shared `open` receipt records actual mass, linear and
angular momentum transfer. Invalid data, kind/wall changes and cancellation
retain all previous controls. Targets are piecewise held per tick; no waveform
interpolation, finite reservoir, moving topology or pressure-work law is added.

The live-target cases in `tests/open_partial_saturation3.rs` target warmed
solver parity against independently constructed native grids for all four
admitted collision models, physical force-frame preservation, independent
mass/linear/angular inventory, spatial-profile precedence and atomic refusal.

`BoundaryGrid3::set_velocity_profile_with_speed_limit` adds an inclusive
caller limit on the final corrected velocity trace. The ordinary profile entry
retains its previous arithmetic and admission. Both entries share the same
deterministic sampling, original low-Mach checks and two reusable plane
buffers. The added bound is checked after mean-normal correction and before
publication; invalid limits or corrected overspeed preserve the full previous
profile and populations. No new per-step field, collision law or momentum
receipt is introduced. The spatial-shear regression in
`tests/open_partial_saturation3.rs` exercises actual body spin and independent
mass, linear and angular inventories with all four admitted collision models.

`BoundaryGrid3::set_pressure_profile` holds spatial positive finite density
targets on declared pressure faces. It reuses the profile owner's deterministic
sampling and two scalar plane buffers. Source, density, geometry or cancellation
refusal retains the complete previous trace and accepted populations. Successful
profiles lock solid topology. `clear_pressure_profile` restores the held uniform
density; changing that uniform target while a profile is installed refuses.

The original regularized pressure operator changes only its per-node density
selection. First-interior velocity/stress, Guo framing and mixed wall/open link
ownership remain authoritative. At mixed seams actual density need not equal
the complete prescribed trace. Updates advance no populations or physical time;
the following common step's shared `open` receipt measures actual mass, linear
and angular transport. The report contains sampled count and minimum/maximum
lattice density, without a pressure-work or thermodynamic interpretation.

`tests/pressure_profiles.rs` covers all six face orientations, uniform
trajectory identity including mixed seams, spatial density reconstruction,
buffer reuse, target precedence and atomic refusal. The pressure-history case
in `tests/open_partial_saturation3.rs` drives an elongated body's translation
and spin, compares a uniform-pressure control, and independently closes mass
and all six momentum components through the original four collision models.

### Moving finite film carriers on the dry sphere/fluid clock

`MovingLubricationSurface3::{Prescribed, Finite}` and the moving/dry
`step_with_*lubricated_spheres_and_surfaces` entries attach each finite plane
to the original `MovingPartialSaturationBody3` using a body-frame point and
normal. Original cached full inertia and momentum-preserving free half-drifts
advance both pose and world mobility. Actual midpoint geometry enters the
existing simultaneous film/fluid solve. The direct received film spin impulse
and authored center couple kick the carrier's momentum before its second free
half-drift. Prescribed translation also advances once to the accepted endpoint.

The dry 6N retained-fluid correction evaluates every candidate from the same
accepted finite-carrier state. Its mobility participates in every shared
resistance response; probes never publish or advance a surface. Accepted
spheres, surfaces and the single streamed fluid candidate publish together,
after all allocations, endpoint geometry checks and diagnostics. Exact surface
load counts, the 64-surface ceiling and original physical-film/scalar-mode
budgets remain authoritative; prescribed load slots must be zero.

Conservative angular-travel bounds exclude sphere/plane contact and duplicate
planes throughout both numerical half-drifts, including the exchanged second
half rate. They can refuse valid near-gap trajectories and are not exact CCD.
Finite reports separate actual fixed-origin momentum change, authored loads,
film work, kinetic change and native free-rotation energy defect. The existing
dry sphere/exterior-fluid balance remains an open subsystem containing all
received surface reactions; finite carrier work is not prescribed reservoir
work or deposited heat.

Each finite slot is an independent film-only carrier. No resolved wall,
sampled carrier interior, displaced volume, buoyancy, inferred weight,
carrier-carrier shape contact, shared multi-face body identity, moving-wall
cover/vacate conservation or thermal destination is supplied. Refusal preserves
all accepted physical states; reusable numerical scratch may change.

Focused native tests cover moving full-inertia recoil beside a spinning
prescribed surface, independent population/body momentum closure across ticks,
origin covariance, exact legacy prescribed/empty continuation, dry probe and
late-carrier-refusal rollback/retry, and swept collisions hidden by clear
endpoints. These tests were authored but could not be compiled or executed in
the delivery environment: the pinned runtime compiler was killed before this
crate, and the execution environment subsequently became unavailable.

### P1 viscous strain through native liquid ingress

Existing tetrahedral-liquid handoffs carry the volume-weighted symmetric P1
velocity gradient through their original tet/cell overlaps. The same native
Hermite reconstruction and Guo force-frame owner initialize the receiving
viscous second moment without changing represented density or momentum.
The required pre-collision nonequilibrium tensor uses molecular tau, while
physical viscous stress includes the original tau-minus-half correction.
Zero strain preserves the original populations exactly; actual velocity and
stress moments are read back before each cell is accepted.

Frozen viscosity uses its original admitted relaxation. Thermal and
generalized-Newtonian ingress query the original source at projected enthalpy
and strain before attaching the existing constitutive owner. Receiving-cell
queries consume the original explicit source-query budget. An admitted local
relaxation can seed that same constitutive solve after coverage, coordinate,
finiteness and relaxation-domain checks. The final source solve replaces every
provisional sample before publication. Ordinary cell constructors and default
thermal attachments retain their existing initialization semantics.

These local seeds avoid false source-domain queries caused by an unrelated
global reference guess. No second viscosity law or retained redundant state
map is introduced. Actual per-fragment mass, linear/angular momentum and
mechanical projection accounts are rechecked after thermal/force attachment.
Numerical remap energy remains separate from enthalpy; no elastic release or
stress-history transfer is implied. Nonlinear stress is evaluated at the
projected strain rather than claimed to equal an average source stress.

Focused regression definitions independently recover forced population
momentum and the viscous tensor, translating/rotating/shearing P1 fields and
continued native motion, population-positivity refusal, and a tight material
rate domain whose valid initialization would fail from the reference guess.
Nonpositive populations, unsupported sources, exhausted query bounds and
cancellation return no candidate. The broader deforming/melting boundary,
pressure-projection and elastic-history obligations remain open.

### Simultaneous liquid throughflow on the original native inventory

`MaterialFreeSurfaceEnthalpy3::{exchange_interface_throughflow,
prepare_interface_throughflow}` accepts incoming and powered outgoing streams
in one original mass/enthalpy/population transaction. Equal gross inlet and
outlet masses can transfer nonzero momentum and caloric energy while keeping
the native fill unchanged. The shared owner creates one lattice/enthalpy
candidate, applies each cell's velocity once, refreshes the original material
and force sources once, and gates actual native endpoint balances before
publication. A prepared proposal borrows the accepted owner exclusively and
can be discarded or published with the finite source owners.

The selected local law holds the streams' velocity, enthalpy and relative
mass-flow proportions over integrated progress `s` in `[0,1]`. For gross
incoming mass `a` and outgoing `b`, `m(s)=m0+s*(a-b)` and
`p(s)=m0*u0+s*(a*v_in-b*v_out)`. Both directions use this SAME `u(s)=p(s)/m(s)`.
Inlet mixing heat integrates `a*|v_in-u(s)|²/2`; outlet actuator work integrates
`b*|v_out-u(s)|²/2`. Each co-located stream group's central velocity variance
adds its own nonnegative gross-mass energy. Thus
`delta_K + mixing_heat = carried_K + actuator_work`; mixing heat enters the
original caloric owner once, and outlet work remains an external debit.
No inlet-first/outlet-first split chooses a different finite-transfer loss.

Analytic Gram weights evaluate the relative-speed integral as a sum of
nonnegative squares. A bounded 64-term positive series handles near-equal
initial/final mass; deterministic `fs-math` logarithms handle separated masses.
The equal-mass limit is regular. Double-double signed stream sums preserve
small net momentum/enthalpy beside common velocity or caloric references.
Single-direction cells delegate the original mixing/ejection state arithmetic.

The native net mass projection remains authoritative. Its rounding adjusts
only the larger directional mass batch; the smaller batch is unchanged.
Each per-cell receipt reports requested and represented incoming/outgoing
masses separately, so external owners debit/credit their own gross streams.
Gross representation defects and their net closure consume the original mass
budget. No effective velocity or enthalpy is formed by division by net mass.
Original stream-count, source, speed, capacity, work, cancellation and complete
state-publication bounds apply. The accepted clock advances only through the
original subsequent flow/thermal step.

This is a local held-stream complete-mixing/ideal-actuator model. It does not
resolve jets, infer nozzle pressure work, evolve finite donor mechanics or
remove existing phase-front/contact/surface-history restrictions. Positive
initial and final native mass bound the linear mass path; source evolution
and detailed spatial transport remain with their original owners.

Six new numerical regressions compare the simultaneous path with an
independent RK4 control-volume integration and an analytic equal-mass pump,
and check one-way parity, individual stream dispersion, frame/datum/order
covariance, near-balanced projection, extreme mass-ratio limits and refusal.
All six and eleven existing stream numerical tests pass in a bounded native
build of the actual numerical source with genuine `fs-math`. Three additional
material/population regression definitions exercise zero-net continuation,
complete original one-way results and atomic refusal/publication; those full
fixtures could not be compiled or run with the available dependency products.
Full Cargo/RCH/DSR validation remains unrun; numerical-kernel execution is not
claimed as end-to-end material or experimental validation.

### Multiple missing-film faces on one finite carrier

`CarrierLubricationFace3`, `SharedLubricationSurface3` and
`BoundaryGrid3::step_with_lubricated_partial_saturation_shared_carriers`
attach several original planar film geometries to one separately supplied
`PartialSaturationBody3`. A carrier owns one mass, center, full world inertia,
translation/spin state and authored external impulse/couple. Face identity is
explicit; equal floating-point body values do not infer shared inertia.
Each carrier must own at least one face. Prescribed translation/axial-spin
surfaces can participate in the same solve.

Modes on different faces of one carrier now share the off-diagonal spatial
mobility terms in the original simultaneous resistance matrix. All local face
wrenches are accumulated before one carrier endpoint is solved. Its external
load is applied once, and each face receives work at that common midpoint
translation/spin. Canonical private first-face references must agree on the
complete carrier and external load before matrix assembly. The original
ungrouped response retains its arithmetic. Carrier factorization is reused
across its faces; no additional physical film modes or second fluid stream
are introduced.

The public result separates `CarrierFaceExchange3` surface reactions from
unique `CarrierLubricationExchange3` mechanical inventories. Per-face force,
world-origin moment, free couple and work remain attributable, while carrier
kinetic change, external work and complete momentum/work residuals appear
once. Center-spin impulse is accumulated directly from local arms. The
original composite plane-work field includes both finite-carrier work and
prescribed actuator work; it is not all external reservoir work.

All fallible face/carrier receipt assembly runs inside the original
prepublication callback. Fluid populations, sphere states and carrier states
publish together. Surface geometry remains immutable; refused geometry,
inertia, bounds, allocation or cancellation yields no partial publication.
The existing 64-surface/carrier ceilings, actual positive-gap checks, footpoint
low-Mach admission and total film/scalar-mode budgets remain authoritative.

This is the original frozen-geometry missing-film correction. It does not
create resolved fluid boundaries, finite face extents, moving carrier
trajectories, gyroscopic pose evolution, displaced-fluid closure, contact or
a thermal destination for viscous loss. The original midpoint damping remains
passive but is not L-stable; sharing inertia does not establish time accuracy
for an unresolved stiff film.

The final selected-source native run passed 18 tests: five new shared-inertia
mechanics cases, three new public descriptor/unique-receipt cases, six existing
spin cases and four existing mixed-surface cases. The independent cases cover
cross-face resistance, full-inertia momentum/work, rigid-motion covariance,
origin changes, prescribed surfaces, exact legacy continuation and malformed
owner/cancellation refusal. Original response, Cholesky, DD arithmetic and
cancellation code executed. The public entry point compiled against the exact
original transition signature; populated-grid streaming/publication was not
executed in this seam. No full Cargo/workspace verdict is claimed.


### Nonuniform in-plane wall fields on the original free surface

`FreeSurface3::planar_wall_velocity_field_candidate` and the material
`planar_wall_velocity_field_candidate_si` accept a tangent-velocity pair at
every Cartesian face of one fixed, completely flooded domain plane. The
cyclic tangent axes and row-major field order are explicit. Every crossing
link arriving in one liquid cell uses that face's same quadrature velocity;
opposed diagonals therefore cancel added mass locally even for nonuniform
shear or in-plane stretching. Actual rounded incoming-minus-outgoing
population mass is measured and checked per face, alongside the original
whole-liquid mass gate.

Original reflection, collision, streaming, capillarity and conversion remain
the only liquid evolution operators. Each face returns actual reciprocal
impulse, original link angular impulse about a declared origin, a separately
accumulated residual couple about its face center, and velocity paired with
actual impulse. Force at the face center alone generally loses a transverse
moment; downstream transfer must retain the couple. Its constrained tilt has
zero angular speed in this fixed-normal face kinematics. Mechanical work is
not deposited heat or a whole-fluid energy certificate.

The same traversal measures the two diagonal tangent-impulse derivatives.
Cross-face and cross-tangent entries are zero for this single original
D3Q19 reflection from the unchanged input state; no steady-state or future
fluid response is claimed. A changed field still requires native replay and
endpoint admission. Added storage is linear in wall faces, with explicit
face and complete-domain ceilings checked before the native fluid clone.
Local Mach, relative wall/liquid speed, population positivity, mass,
representability and cancellation checks precede candidate publication.

The material adapter retains the original source, units and scenario speed
limit at every wall face and actual liquid endpoint. Raw population impulses,
reference-pressure impulses/moments, perturbation loads and optional absolute
loads remain distinct. Geometry has fixed normal position; this does not add
cover/vacate, normal motion, moving contact lines, finite plate edges, a solid
constitutive law or complete MR69 coupling. The existing rigid-link and
uniform-wall operators retain their population arithmetic.

Six regression definitions cover independent nonuniform link and liquid
moment balances in all six orientations, per-face couples/work/mass, affine
replay, origin covariance, continued native motion, exact uniform-path parity,
physical source/pressure units and atomic refusal/retry. Targeted Rust syntax,
format and whitespace checks passed. Native typechecking and test execution
remain outstanding in the delivery environment.

## Finite-inertia planar wall rotors

`WallRotor3` and `BoundaryGrid3::step_with_finite_inertia_rotor` advance a flush
circular disk's actual spin together with native D3Q19 collision and wall
streaming. The grid stays inertial. Disk center and Cartesian spin axis are
fixed by a bearing; disk inertia uses lattice mass-length-squared units.
Positive spin uses the positive Cartesian axis on either face. Cell centers
are at `[x+0.5,y+0.5,z+0.5]`, and wall planes at zero or the grid dimension.

Each selected face cell owns its complete five-link stencil. Open/periodic
domains, moving housing walls, side-wall seams and obstructed patch cells
refuse. Complete tangent pairs cancel moving-wall mass and moving-mass
momentum corrections. Native half-link impulse yields
`J_h=A-D*omega_mid`, where `D=sum(rho*r_perp^2/3)`. The analytic midpoint
inertia solve has no explicit added-inertia step restriction. The held
midpoint spin drives the original bounce; bearing force and transverse
torque are returned separately, with all moments about the requested origin.

The receipt closes fluid, rotor and housing momentum against the applied
rotor couple and held fluid body force. Rotor kinetic change equals midpoint
external work plus measured fluid-on-rotor work to roundoff; fluid
viscous/numerical dissipation is not assigned to that mechanical identity.
Reusable patch metadata and next-state storage avoid dense per-step
allocations after warmup. Cancellation and late refusal publish neither
fluid nor rotor state. Call the coupled method each step; ordinary grid
stepping treats the underlying face as stationary.

`tests/finite_inertia_rotor3.rs` covers native impulse slopes on all six faces
and two radii, small inertia, all four collision models, arbitrary-origin
source/bearing balance, fluid-driven spin-up, topology admission and atomic
late refusal. `examples/finite_rotor3.rs` runs unforced spin-down and reports
actual fluid/rotor/housing angular exchange.

### Simultaneous independent wall rotors

`BoundaryGrid3::step_with_finite_inertia_rotors<const N>` advances a nonempty
fixed array of rotors against one shared native collision and streaming step.
All midpoint responses use the same post-collision fluid state. Complete
patch cells must be disjoint; overlap refuses before collision. Rotors may
share one face, occupy opposite faces, or rotate about different Cartesian
axes. A fixed-size cursor merge caches the next sparse patch address, so
ordinary fluid cells do not scan the rotor list and no dense owner map is
allocated.

`WallRotorsExchange3` reports per-rotor `WallRotorResponse3` entries in input
order, with stationary-housing and fluid body-force exchanges counted once
globally. Each rotor retains its own applied couple, bearing reactions and
midpoint work identity. The single-rotor API delegates to the same kernel and
preserves its receipt. All proposed fluid and rotor states publish together;
a later rotor refusal or cancellation accepts none. Tests cover bitwise
one-element equivalence, source-free counterrotation, three spatial axes with
body force/applied couples, same-face disjoint patches, overlap and atomic
late refusal.


### Thermal distributed wall momentum on the original VOF owner

`MaterialFreeSurfaceEnthalpy3::step_with_planar_wall_velocity_field_si` and its explicit molecular-heating variant advance the original native free-surface candidate with the selected face velocity field, current cellwise material relaxation and retained Guo force frame. Native collision/streaming/conversion records the same bounded mass routes used by the existing extensive-enthalpy transport; no reconstructed companion flux or second fluid/thermal engine is introduced.

The original thermal source, optional rheology/buoyancy, transport, wetting and enthalpy gates precede publication. The actual final native candidate is then readmitted for flooded-wall topology, population/Mach and ledger mass, followed by material scenario speed and SI conversion. Populations, enthalpy, material/rheology state, wetting and clock accept together only after cancellation polling.

Global and local reciprocal impulse, angular impulse, face couple, response, mechanical work and pressure-reference reports are retained. Wall work is not deposited as caloric heat. The optional heating method reuses the existing held initial-state molecular stress-power source; its existing first-order splitting and thermal-domain limits remain.

Scope is a fixed, completely flooded planar wall with tangential face velocities. Normal cover/vacate, finite moving edges, moving contact lines, prescribed gas-pressure actuation and a solved conduction boundary are not supplied by these new calls. Existing electrical and thermal surface-boundary ownership requirements remain operative.

Five synthetic regression definitions cover zero-field parity with thermal and buoyant flow, actual mass-route enthalpy transport and work separation, one-time molecular source reuse, native-stress rheological closure, and transactional source refusal/cancellation/retry. Formatting/parsing and scoped whitespace checks passed. Native Rust regressions, convergence and performance were not executed in this resource-constrained session; these definitions are not a runtime validation or performance claim.


