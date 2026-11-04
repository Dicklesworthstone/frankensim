# CONTRACT: fs-lbm

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
are unchanged. Native pressure-only APIs retain their documented legacy scope.

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
Legacy native mass-only stepping retains its historical behavior.

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
existing sparse VOF wet-cell BGK relaxation hook. This complements the
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

## Exploratory Smagorinsky transport on the native sparse free surface

`FreeSurface3::prepare_sgs` and `step_with_sgs` accept explicit shared
`SgsConfig` controls and a reusable `FreeSurfaceSgsWorkspace3`. This path
requires native BGK or TRT and Smagorinsky; central/cumulant collisions refuse.
WALE and Vreman need information absent from a symmetric local moment and
refuse. Acceptance use, invalid controls, and excess eddy/molecular viscosity
refuse without changing accepted flow. The zero coefficient follows the
original native BGK or TRT collision path bit for bit.

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
Both BGK and TRT still tie shear/bulk relaxation.

Morton-key/lane cache entries survive sparse slot changes; newly wet cells
are resolved on the next preparation. Only one native candidate is cloned.
Its original collision, gas reconstruction, capillarity, mass exchange and
conversion finish before publication, so native or cancellation refusal
preserves every accepted component. Scratch is consumable only on success.
`FreeSurfaceSgsReport3` describes the pre-step input: instantaneous resolved
and modeled shear dissipation are weighted by native bounded liquid fill.
Empty interface ghosts and unlocated carry add no liquid dissipation. No
diagnostic power is added as heat to this isothermal owner.

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
time integrator, interfacial storage or additional source of heat. Duplicate
face ownership, compact aliases and through-flow use the existing face-contact
admission. A failed query, cancellation or exhausted iteration publishes no
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
itself invert to pure liquid at the original fixed caloric density. The incoming
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
  3-D) explicitly adds `mu * rate^2 = 2 mu D:D` from those same cached
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
  numerical momentum dissipation, pressure work, compressible/bulk-viscosity
  heating and elastic storage are excluded. A sensible-heat balance is not
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
  toward the equilibrium of the authored far-field target,
  `f <- f + sigma (f_eq - f)`, with a quadratic ramp
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
  history.
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
  (the hash IS the seed schedule); `state_bits()` — canonical-order exact
  bits for golden hashing; `total_mass()`, `cell_macros(slot, lane)`,
  `active_tiles()`, `allocated_state_bytes()` (exactly three population
  buffers per active tile), `state_bytes_per_tile()`.
- `SparseError3` — typed refusals: `Dims`, `TileOutOfDomain`, `Collision`
  (Morton key + lane + the underlying `CollisionError3`), `Cancelled`
  (pre-step state intact; re-issue is deterministic), `Pool`.

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

Normal derivatives use the quadratic through the two nearest fluid centres
and the actual halfway boundary condition. Each preparation supplies both
uniform wall velocities and independent scalar `Value` or `ZeroGradient`
conditions. A scalar zero-gradient wall requires three normal fluid cells:
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
in force. This adds no moving/open wall model, material-wrapper API, thermal
energy owner or turbulent-DNS qualification.

### Resolved vorticity on native three-dimensional flow owners

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

### Transient finite-volume momentum and heat

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
combined or overwritten. No SGS dissipation heating, turbulent scalar-variance
equation, equilibrium thermal wall law, calibrated turbulent channel profile,
or resolved turbulent spectrum is claimed. The existing first-order thermal
coupling and power-law/upwind convection limits still apply; the caller must
establish temporal and spatial convergence for the intended exploration.

Focused `conjugate::simple::sgs` tests cover graded Couette and voxel walls,
quadratic Neumann gradients, native MAC normal derivatives, WALE's pure-shear
zero, physical Smagorinsky viscosity, exact disabled flow identity, admission
and cap refusals, and a moving-wall flow/heat step compared against an
independently supplied turbulent-conductivity energy solve.

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
moving wall velocities and the original uniform isothermal or adiabatic wall
condition enter each fit. Adiabatic physical scalar derivatives and the
Favre-filtered scalar's one-sided normal derivatives remain distinct, as
specified by the shared channel model.

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

Four focused native channel regressions cover current source heat capacity and
force/wall data, heat-only planes, molecular wall conduction, shared-face
admission, moving-wall energy, worker identity, workspace reuse, refused-step
rollback, zero controls and topology/budget refusals. The periodic regressions
retain their original physical assertions. Interior solids, spatially varying
wall data, wall functions and channel DNS calibration remain outside this path.

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
The plane belongs to the Mann frame; every actual exterior ghost centre must
lie in its upper half-space. Mean velocity must be tangential, and the frozen
adapter also requires tangential convection because it translates query
coordinates. Tangency uses the admitted rotation's roundoff scale without
altering the supplied velocities.

The wall observation reads the same frozen snapshot or evolving trial
quadratures and adds no numerical arrays. Rejected samples and cancellation
preserve the accepted profile and atmosphere; retries use the original
temporal owner. Box/face averaging combined with wall blocking refuses in
both attachment orders because the homogeneous sinc gain does not average
the evanescent wall correction. This is point-sampled slip-plane conditioning,
without a viscous no-slip layer, wall-generated turbulence or calibrated
near-ground spectrum. Four native regressions cover rotated profiles, exact
wall-normal values, incompatible settings and evolving rollback/replay.
All four pass in the native debug build, together with the ten existing Mann
inlet tests covering original, filtered and evolving boundary behavior.

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
- `VelocityPressureX2` is a low-Mach regularized fixture boundary with periodic
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
Temperature collision samples the input physical velocity after force rebinding
and before flow advances under the frozen buoyancy/external load.

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
`-beta (T - T_ref) g` in lattice units, advances `coupling_interval`
pooled steps, projects the physical cell mass fluxes including `F/2`, and takes one implicit
energy step whose storage uses the fluid's volumetric heat capacity in
every cell (solids get it as a PSEUDO capacity: it shapes the path, never
the fixed point). Convergence is judged on the geometric extrapolation
`delta r / (1 - r)` of the remaining per-coupling change (no estimate
until the sequence contracts), for the relative mass flux and for the
temperature over its range. The answer is a fresh STEADY energy solve on
the final fluxes; `coupling_residual_k` retains its distance from the last
pseudo-time iterate. A lattice speed above `max_lattice_speed` refuses as
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
junction-to-case resistor. Its area-specific resistance is `R_junction * A_total`;
the existing analytic film law resolves the outside surface without a nested
iteration. Each face receives its correct parallel area share, and the collapsed
package adds no bulk resistance. Frozen material fields and fluid links retain
their declared coefficients. Junction powers, fixed-cell reactions and observers
all use the same stored links. Ordinary boundary conditions on collapsed component
sides contribute no heat, matching the assembly's adiabatic side ownership.
`tests/conjugate/nonlinear_compact.rs` checks independently derived heating and
cooling profiles with unequal face areas, branch powers, and the compact-side
heat-balance regression. Physical marches use explicitly declared compact
junction heat capacity; no capacity is inferred from a resistor. The analytic film is the declared stationary
scalar Robin law, not a resolved fluid boundary layer.

The existing steady and transient Picard/Krylov owner evaluates these
coefficients from the previous pass and warm-starts the next solve until the
largest temperature change is below 1e-10 of the span. This update only proposes
acceptance. The original assembly and physical accounting are then evaluated
at the unchanged candidate with its own conductivity, nonlinear half cells,
contact and compact links. No extra Krylov solve or time advance occurs.
The actual row-scaled residual and every free cell/junction's heat closure must
meet `EnergyConfig::tolerance`; each heat residual uses its own throughput and
an explicit local product-roundoff allowance. Opposing defects cannot cancel
through the global ledger. Storage retains the original previous-time state.
If the local heat gate needs more accuracy, the next counted linear solve uses
a tighter tolerance; the requested physical endpoint allowance stays fixed.
The same 200-pass bound refuses an unresolved endpoint. Returned conductivity,
wall links, junction branches and reservoir heat describe the accepted field;
the iteration report sums the actual native conductivity Krylov work.
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
sources and iterates potential and temperature (resistivity at the
previous pass) until the largest change is 1e-9 of the span; 200 passes
refuse as `SolverNotConverged { system: "electrothermal (thermal
runaway?)" }`, the answer for a current-driven conductor past runaway.
Refusals: no voltage electrode, overlapping or touching electrodes, a
current electrode with no conducting path to a voltage electrode, a
non-positive resistivity at the solved temperature.

| Fixture | Reference | Measured |
|---|---|---|
| 10 mm x 1 mm^2 bar, rho 1e-6 ohm m, 10 mV, ends held at 300 K by copper electrodes | `I = V A / (rho L)` = 1 A; uniform 1e6 W/m^3; peak `300 + q L^2 / (8 k)` | 1.000000000000001 A; Joule = delivered = 0.0100000000 W; peak 301.250000005 K |
| Same bar current-driven, `alpha` 3.9e-3, `m L / 2 = 0.5`, 10/20/40 cells | `theta = (cos(m (x - L/2)) / cos(m L/2) - 1) / alpha` | peak errors 3.33e-2, 8.32e-3, 2.08e-3 K: order 2.00; 11 passes; delivered = Joule to 1e-9 |
| Same, 3.5 x the current (`m L / 2 = 1.75 > pi / 2`) | no steady state | refused (electrothermal) |

No-claims: DC only (no skin effect, inductance or capacitance); ohmic,
isotropic conductors; no electrical contact resistance; ideal electrodes.

### Two-resistor compact components (`CompactComponent`)

`ThermalSetup::compact_components` holds JEDEC two-resistor models
(JESD15-3), the standard board-level representation of a package from its
datasheet `R_jc` and `R_jb`: the component's box of solid cells blocks flow
like any solid, but its interior is collapsed. One extra unknown per
component, the junction, receives the power and reaches every cell beyond
the case-top face through `R_jc` and every cell beyond the board face
through `R_jb`. Each resistor is split over its face in proportion to area
and placed in series with the outside cell's half width:
`G = 1 / (R A_face / A_cell + dx / (2 k A_cell))`. The sides are
adiabatic; the box cells report the junction temperature through identity
rows. `EnergySolution::junctions` gives the junction temperature and the
case and board heats; `solid_to_fluid_heat_w` counts the case links into
fluid instead of the collapsed cells. Refusals: a box that holds fluid,
overlaps or sits on another component, or puts its case or board face on
the domain boundary; power, sinks or fixed temperatures inside a box; a
physical transient without an explicitly declared positive junction capacity.

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
SIMPLEC iterations under the current force, then re-solves energy on the
current fluxes; convergence needs the SIMPLEC residuals below
`flow.tolerance` AND the energy re-solve's largest temperature change below
`temperature_tolerance` times the span, at the same coupling. The returned
energy solution is computed on exactly the returned (tightly projected)
fluxes. This closes the LBM path's "closed enclosures only" boundary.

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
drop other properties. Constants and temperature-only curves are supported,
including an explicitly typed absolute-temperature axis under any source name.
Other response curves and dynamic frequency/rate/strain bounds refuse.

The existing material momentum-scaling helper owns nu_lattice and tau for both
thermal and free-surface consumers. Explicit dx/dt convert rho0*g to the native
constant Guo force density and sigma to sigma/(rho0*(dx/dt)^2*dx). Construction
bounds the full possible voxel domain, validates tile coordinates and refuses
initial liquid inside a permanently solid tile. The source-bound wrapper returns
native ledger, regional mass and unregioned carry in kg, and wet-cell density,
Guo-corrected velocity and incompressible pressure perturbation in SI. Gas
populations are never presented as liquid samples.

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
No collision/curvature/pressure evolution is changed for admitted states.
The helper `maximum_capillary_timestep_s` exposes the physical maximum before
constructing a plan; `capillary_timestep_ratio` reports its consumed fraction.
Explicitly neglected tension has no capillary limit. Source-temperature or
interface extrapolation increases that exceed the bound refuse before accepted
fluid/thermal state publication. The focused G0/G1 tests check independent
wave dispersion, dx^(3/2) and material scaling, positive-pressure refusal,
zero-tension behavior and unchanged state after rejected coefficient updates.
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

The discrete account is `W_wall + W_body + W_pressure = Q_molecular + Q_Darcy + (K_end-K_start) + N_BE`, where `K=sum(M_face*u_face^2)/2` uses exactly the inertia row's staggered half-cell volumes, and `N_BE=sum(M_face*(u_end-u_start)^2)/2`. The nonnegative backward-Euler numerical loss is reported, never turned into heat. Pressure-constraint work must independently remain negligible; it is not heat or a residual correction. Only physical molecular/Darcy watts enter the original caloric kernel, once. `UnsteadyHeatingAudit` accumulates that physical heat, signed supplies, native kinetic change, numerical loss and residual over the accepted advance; generated heat is already inside the original thermal audit. The older cell-centred kinetic diagnostic is not substituted for native mechanical storage.

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
quadrature observes actual cell-volume RMS `sqrt(2 D:D)` from the returned
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

The cell-volume RMS `sqrt(2 D:D)` is observed from the exact retained face
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

The existing native smoothed VOF field supplies the actual wall-adjacent fill
and in-wall gradient. For outward liquid normal `-grad(phi)/|grad(phi)|`, the
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
rate `sqrt(2 D:D)` changes. The dynamic/kinematic viscosity basis and sourced
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
mu and strain supply the held `mu*gamma_dot²` heat; heated material then refreshes
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

`free_surface3_from_tetrahedra_with_enthalpy_and_rheology` carries the original temperature/shear material law into the existing thermal VOF owner. Each source tetrahedron supplies its actual affine velocity gradient and `sqrt(2 D:D)`. Original pinned property-line admission covers that tetrahedron's complete P1 temperature span at its actual rate; a reference rate cannot conceal unsupported source shear. Shared vertices are admitted separately when adjacent tetrahedra have different strain. Source-query and original geometric budgets remain explicit.

The original exact overlap, per-fragment mass, linear/angular momentum and enthalpy transfer remain authoritative. The receiving native owner resolves its own force-corrected shear and source viscosity, retaining the same law through continuing motion and molecular `mu*gamma_dot^2` heat. Source rates are reported separately from receiving constitutive samples. The receiving viscous stress is initialized from the overlap-projected symmetric P1 strain and the original source at the receiving enthalpy and rate. Nonlinear stress is evaluated at that projected state; an average of source stress tensors is not implied, and elastic memory remains outside this handoff. Kinetic/potential projection discrepancies remain numerical accounts, never heat. The existing `with_time_origin` preserves an explicitly supplied source clock before continuation.

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

The native receiving owner retains its reference gravity load and adds the
bounded source-density departure. Its Guo-frame rebind preserves physical
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
dry numerical ghosts retain the reference coefficients. Source tension remains
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
different chemical environment. Surface tension remains the native constant.
Both phase capillary bounds, source support, native transport and cancellation
must pass before any liquid, finite gas or wall candidate publishes.

The original reference-ratio table is reused exactly when ratio bits match.
Changed moving ratios use one temporary table from the SAME bounded Gauss/
Hermite owner; pinned or zero-speed samples do not construct it. There is no
growing ratio cache. These constitutive evaluations remain bounded by the
existing wall-cell and thermal-substep limits. This supplies a two-viscous-wedge
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
stencil evaluates `gamma = sqrt(2 sym(grad u):sym(grad u)) / dt` in 1/s. The
original producer resolves one material response per fluid cell, its transport
view is frozen once, and `tau = 1/2 + (mu/rho0) dt / (cs^2 dx^2)` drives that
collision/stream step. Source objects and query coordinates, time and rate are
retained. Reference dx/dt/rho0 stay fixed; nonmatching density refuses. D3Q27
independent bulk/higher rates and all equilibrium choices survive; D3Q19 moment
models retain their existing second-order bulk/shear tie.

This is explicit incoming-time generalized-Newtonian feedback, not a
self-consistent nonequilibrium-stress closure or implicit time integration.
No source-domain extension, viscosity floor/cap, material fit, heat deposition,
thermal evolution or hereditary memory is added. The preparation sweep caches
velocities in O(N) storage/work and queries sources in fixed serial cell order.
Failed source, scaling, collision or cancellation checks publish no step.
Time-varying force endpoints remain owned by the existing trapezoidal force
step and are not used to reinterpret incoming strain.

Five native exact-source serial tests exercise the actual grid owners and real
execution context: all seven Newtonian collision variants reproduce existing
population bits; force-centered Carreau queries ignore rigid rotation;
quadratic wall shear is recovered; original source/query ownership and rollback
are retained; and a CM power-law channel with index 0.75 reduces maximum
analytic-profile error from 0.723% at 12 fluid rows to 0.190% at 24. Every source
query stays outside the fixture's explicitly declared low-rate regularization.
Pooled wrappers reuse the same preparation and existing pool owner, but pooled
execution was not included in this qualification; no cross-ISA or throughput
claim is added.

### Local capillary-interface conversion preserves detached drops

Both dense D2Q9 and sparse D3Q19 free surfaces redistribute each converted
cell's signed excess only to adjacent surviving interface cells after the full
flag cascade. A conversion no longer changes every distant drop's fill and
Laplace pressure through a global pool. Recipients follow fixed lattice-link
order, use fixed-size neighbor storage and take the final rounding remainder
locally. The redistribution cost depends on converted cells and their stencil,
without another whole-domain interface scan.

New D3Q19 constructors select `ConversionRedistribution3::Local`.
`with_conversion_redistribution(LegacyGlobal)` explicitly reproduces historical
version-1 trajectories, including their nonlocal pooling. The frozen pour
sentinel selects this replay mode and retains its original four-quadrant hash;
new locality and thermal regressions exercise the corrected default. The
independent `FREE_SURFACE3_BIT_SEMANTICS_VERSION` is 4 for the current default;
version 2 introduced local redistribution, and version 3 corrects the density
and force handling of fresh-cell initialization described below. Version 4
uses the physical Guo frame in gas-pressure reconstruction as described above.
Explicit `LegacyGlobal` replay retains its historical conversion and
gas-pressure arithmetic, leaving the unchanged D3Q19 core, boundary and sparse
goldens untouched.

D3Q19 records direct donor/recipient mass routes, including reversed routes for
conversion debt. The original extensive enthalpy transport consumes these same
routes, avoiding temperature mixing through a shared zero-storage pool. Excess
or debt with no adjacent recipient remains explicit unregioned carry; without
spatial provenance it is not injected into another component on a later step.
This preserves the total numerical mass ledger but does not claim resolved
geometry for orphan carry, momentum conservation during cell conversion or a
conservative capillary surface-energy discretization.

Focused G0 tests cover both conversion signs, disconnected-drop locality,
per-node mass-route balances, orphan carry retention and scalar D3Q19 tension
admission. The G1 enthalpy regression compares two independent conversion
components and reverse-debt cycles against analytic local energy balances.

### Bounded local conversion carries physical liquid momentum

`ConversionRedistribution3::ConservativeLocal` is an explicit native timestep
option, also selectable on `MaterialFreeSurface3` or the already bound
`MaterialFreeSurfaceEnthalpy3`. Existing Local and LegacyGlobal trajectories retain their arithmetic.
The new mode completes signed over/undershoots at zero conversion tolerance and
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
co-moving identity and late public-step rollback.

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
mass moments. The historical `step` retains its existing conversion behavior.

Seven focused G0/G1/G3 regressions cover unequal-density advancing and receding
fronts, signed debt, physical boosts and forcing, independent droplets, periodic
capacity aliases, retained physical stress, native zero-conversion replay,
cancellation and newly wet-cell refusal.

### Material D2Q9 defaults carry conversion momentum

`MaterialFreeSurface2::step` and `MaterialFreeSurfaceWetting2::step_with_gate`
delegate to their existing `step_with_momentum` owners. Ordinary material
stepping therefore carries signed donor momentum through actual adjacent
conversion capacity, admits the existing SI conversion account, and publishes
one complete candidate. Dynamic wetting publishes the same candidate together
with its accepted source queries, contact geometry and angle history.

The native `FreeSurface::step` remains the historical replay entry. Material
`step` preserves its convenience signature and panics on a conservative refusal
without publishing the proposal; the existing explicit momentum entry provides
typed cancellation and refinement refusals. The wetting convenience entry
retains its original error type and canonical cancellation classification.
No additional lattice clone or duplicate conversion-report owner is introduced.

`resolved_linear_momentum_kg_per_s` sums actual liquid mass times Guo-corrected
velocity and converts to SI momentum per unit out-of-plane depth (kg/s). Gas,
walls and unlocated carry have no represented liquid momentum. External forces
and wall/pressure traction can change the observation; it does not certify
whole-step or angular-momentum conservation.

The material regressions compare default and explicit conservative evolution
through actual nonzero conversion routes, check SI momentum against physical
cell observations, and preserve wetting-source refusal and cancellation
atomicity. Native capillary fixtures cover all four wall orientations and the
existing material parameters; the gravity fixture starts at physical rest in
the held-force frame before evolving real drop conversions.

### Density-consistent physical state at free-surface ingress and conversion

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

### Local thermocapillary traction on detached interfaces

`free_surface3_from_cells_with_local_thermocapillarity` retains the original
material, enthalpy, native pressure and Guo-force owners while removing the
global tank-column restriction. `LocalThermocapillarity` bounds the smallest
positive cut-cell fill and the conditioning of the local source reconstruction.
Complete monotone raw-VOF height stencils share their derivatives with native
curvature; local traction may try bounded axial depths of 7, 11 and 15 cells.
The selected graph must have a normal-axis magnitude of at least one half.
Native pressure keeps its existing seven-cell curvature stencil and fallback.
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
liquid-liquid conduction links and the same material, caloric, gas-film,
convection and radiation laws. Finite bodies, electrical coupling and separate
gas/radiation APIs still require additional coupled residuals and refuse.

Liquid coordinates are specific-enthalpy increments in J/kg, mapped to actual
cell energy increments with the original masses. A finite atmosphere appends
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
Multi-body and film entry points retain their existing held-force semantics.

Focused regressions specify zero-occupancy native parity across four collision
models and periodic/wall/reservoir boundaries, independently reconstructed
physical source work, light/neutral/heavy moving-body ramps and reversal through
zero, cut-cell interior subtraction, exact held-path parity and atomic
overflow/budget/late-cancellation refusal with retry. Rustfmt parsing and source
review are available; native compilation and execution remain unverified.

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
transient material-convection entrypoints. Temperature-only energy and momentum
owners refuse this selection; no caloric chart is inferred from reference Cp.

For each fluid cell and coordinate axis, minmod limits the two specific-enthalpy
gradients using actual adjacent cell-centre distances. A donor half-width gives
the face correction. Two fluid neighbours are required; boundaries, solid
interfaces and unresolved one-cell directions retain donor transport. Both
rows of an interior face use the same signed mass flux, donor tangent and
correction, so the transfer has one credit and one debit. The existing positive
upwind Cp matrix, molecular/contact diffusion, source terms, reservoir laws,
reference mass and native caloric storage remain authoritative.

Three offsets per cell are cached once for each frozen energy operator, outside
its matrix and Krylov work. The original material iteration refreshes source h,
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
caloric march and independent stored-energy accounting; explicit refusal by
noncaloric energy and momentum paths; and stale-limiter rejection even with an
exact affine caloric law. Rust parsing, compilation and native numerical
execution were unavailable in the delivery environment; source/API review
does not establish those executed results.

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
There is no added EOS, phase mobility, bulk-viscosity law, turbulence,
pressure-opening work model or implicit viscosity-temperature feedback.
Cell and sweep ceilings bound work: each sweep has three existing momentum
solves, each capped at 20,000 Krylov iterations, and a frozen coupled residual
readback. Numerical storage is linear in original cells/faces. Any source,
solve, balance, heat-insertion or cancellation refusal returns no candidate.

Four source-bound regressions cover graded actual inertia distinct from source
reference density, original liquid-viscosity changes, moving-wall work and
physical heat, current-state reconstruction, quiescent compatibility and atomic
refusals. All eight changed Rust sources passed pinned-nightly rustfmt parsing
through stdin; the three new files were formatted. Independent scalar response
and 101 unequal-mass scatter identities were evaluated, with the latter closing
within 2.23e-15 J. Native Rust typechecking, tests and runtime benchmarks remain
unexecuted because the shared filesystem is full and memory is constrained.

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


## Rigid wall momentum on the original thermal VOF trajectory

The immutable thermal rigid-wall candidate methods use the existing link-resolved rigid reflection with current cellwise material relaxation and the original Guo force frame. That same native tick records the actual mass routes consumed by extensive-enthalpy transport. Final source/rheology, wetting, wall mass, population/Mach, material speed, SI conversion and cancellation checks precede return.

A private evolved-state tick separates preparation from publication without cloning prepared source or electrical operators. Each thermal rigid-wall trial clones the original native fluid once and retains the same material/caloric bindings and physical clock. The existing electrical/conduction ownership gate remains operative. Optional molecular heating reuses the original held initial-state stress-power source; actuator work and numerical loss are not additional caloric inputs. Fixed flooded-plane and single-phase limits remain. Focused regression definitions exercise isothermal parity, heat and buoyancy, actual rigid replay, original molecular source reuse, and refusal/retry. Formatting, parsing and scoped whitespace checks passed; native compilation and execution remain outstanding.


### Finite guided plane sliders on the physical D2Q9 disk

The native dry-disk owner now accepts FiniteLubricationSlider2 through
DryDisk2::step_with_lubrication_sliders. Each source supplies a positive
physical mass per out-of-plane depth, a distinct center of mass, a signed
tangential velocity and the original stationary plane/film geometry.
The guide fixes normal position and orientation. Infinite-plane geometry
therefore stays exactly unchanged under travel; no finite face extent,
rotary inertia, resolved-fluid wall or free normal recoil is inferred.

Every original three-unknown dry-fluid candidate condenses the finite
slider midpoint with k_eff = k/(1+k/(2m)) and
u_free = u_before + J_authored/(2m). The stiff branch avoids overflow of
k/(2m). The disk receives the film load inside its existing native solve;
the slider receives the opposite tangent impulse and actual kinetic work.
Normal traction and the complete distributed-traction couple go to a
separate ideal guide port, whose work is zero. The original cylinder spin
resistance is retained; translating recoil does not silently remove it.

Slider centers undergo both free half-drifts around the same impulse.
Receipts separately retain authored work, true endpoint momentum/kinetic
change, guide impulse/couple and their measured balance. The combined
exterior-fluid/dry-disk/slider momentum account includes guide reactions
as external inputs. Every slider endpoint, work and momentum check is
prepared before the original sole body/fluid publication, followed only
by infallible assignments of the already-admitted slider states.

Focused native-source checks cover simultaneous oblique sliders with
Guo acceleration and authored loads, independent fixed-origin momentum
and kinetic-work reconstruction, stiff recoil, original inactive/empty
bit parity, late refusal and cancellation with exact retry. The bounded
stiff limit is also checked without forming the overflowing raw ratio.
Implicit midpoint remains passive but is not L-stable: unresolved stiff
films may reverse relative velocity during a tick, so physical relaxation
still needs time resolution. Existing model cutoffs, positive-gap
admission, 64-plane limit, no contact and no film-heat destination remain.

## Native curved D3Q19 wall geometry and momentum

`CurvedWall3` binds canonically ordered original wall links to their actual
fluid-side fractions and caches the opposite fluid donor for fractions below
one half. `from_level_set` checks all native occupancy signs and brackets every
voxel-wall intersection with explicit root tolerance/iteration limits. The
callback's continuity and single crossing on each link remain caller geometry
obligations; no SDF/topology certificate is inferred. An unresolved short-link
gap, stale donor, mismatched mask or exhausted root budget refuses.

`BoundaryGrid3::step_with_curved_walls` uses the original four collision laws,
Guo source frame and pull table, then applies linear Bouzidi-Firdaouss-Lallemand
interpolation only at selected links. Surface velocity is sampled once per
link at its actual intersection. Exactly-halfway geometry retains the original
population and wall-exchange arithmetic. Sparse geometry/cursors and reusable
post/next fields avoid a dense distance/velocity field or fresh per-step
population allocation after warmup.

The SAME incoming population is used for both streaming and the shared wall
momentum convention. Actual population impulse, boundary-relative impulse and
velocity-times-mass correction remain separate; torque uses the off-lattice
intersection, and work uses that same wall velocity. Unselected walls, held
fluid forcing and native periodic angular transport are separately counted.
Interpolation mass drift is measured, never renormalized into a false swept
volume. Every candidate population and load passes finite-state/cancellation
checks before one final swap. Late refusal publishes neither fluid nor loads.

This is fixed-occupancy curved geometry with prescribed surface kinematics,
including tangential rotation of invariant spheres/cylinders. It neither
advances a solid pose nor remaps newly covered/uncovered fluid. Existing
moving-body owners remain authoritative for evolving immersed geometry.
Regularized reservoir faces have their own measured exchange below. There is no general
GCL, full-solver Galilean-invariance, calibrated body-drag, complete fluid-energy
or global second-order convergence claim.

Six native regressions pass, covering original halfway bit compatibility and
global momentum accounts for BGK/TRT/central-moment/reduced-cumulant collision,
both interpolation branches, actual generated sphere torque/work, periodic
transport, malformed/missing donors and late cancellation. Fractional planar
Couette recovers its physical-wall linear profile within 2e-10. Resolved
cylindrical Couette relative velocity errors are 0.01472427871037479 on 24x24x4
and 0.007795637138649035 on 48x48x4, under respective 2%/1% budgets; this measured
pair does not establish asymptotic order. A 2.3-cell inner radius is outside
that accuracy envelope (about 12.4% error). The ten existing native finite
wall-rotor regressions also pass after sharing the incoming-population momentum
helper. These are direct original numerical-owner checks with only unused
runtime/pool and unrelated module wrappers omitted; full workspace Cargo,
DSR, cross-ISA and performance gates were not executed.

### Local material viscosity on the original curved-wall step

`BoundaryGrid3::step_with_curved_walls_and_relaxation` accepts a step-local
x-fastest shear-relaxation field. It uses the SAME native collision, held Guo
force, interpolated pull populations and wall/periodic momentum accounts as
the uniform curved-wall entry. BGK/TRT/central-moment/reduced-cumulant families
retain their original non-shear controls, including TRT's magic product and
D3Q19's existing bulk/shear tie. Unchanged rates retain original bits; a later
ordinary step still uses its constructor model. Solid padding is ignored.
Invalid field dimensions/rates or native endpoint/load/cancellation refusal
leave accepted populations and topology unchanged.

`material::MaterialScaling::step_curved_walls_with_material` samples
each original `FluidTransport` response once, retaining that concrete value
and its exact applied tau. Borrowed source values permit deliberately frozen
fields without repeating property interpolation. The existing momentum
scaling owner converts nu=mu/rho0 using the reference dx/dt. Every local source
must retain the same physical reference density; mass, inertia and fluid
inventories are not remapped. Cp/k remain source-admitted but impose no
thermal-lattice restriction on this momentum-only operation. The explicit
cell budget precedes sampling/allocation. All source and conversion checks
precede native stepping, with no fallible publication after native acceptance.

Temperature/composition evaluation, source epochs and any thermal evolution
remain with the caller. The method does not reconstruct source pressure,
solve implicit shear rheology, deposit viscous heat or substitute actuator
work for heat. Geometry, wall velocity, moment origin, held force and returned
momentum/work accounts retain the native lattice-unit convention. The caller
must use the pinned density/dx/dt mapping and maintain its scenario speed
bound. Curved-wall interpolation mass drift and fixed-occupancy restrictions
are unchanged.

Focused regressions exercise constant-source population/report identity for
all four collision families, independent local collision/pull reconstruction,
two-temperature piecewise Couette fluidity and actual wall work, original
source retention and single evaluation, and budget/source/density/rate/late
wall/cancellation rollback. These new regressions await native execution;
no new measured convergence or performance claim is made here.

## Source-backed gas convection on planar moving-liquid interfaces

`MaterialFreeSurfaceEnthalpy3::step_with_surface_convection` and
`step_with_finite_gas_surface_convection` evaluate a call-local
`FreeSurfaceGasConvection` inside the original zero-storage surface solve.
Each trial receives its actual interface temperature, current original gas
bulk state, and the same selected gas source claims at the mean interface/bulk
temperature and current pressure/composition. A positive finite coefficient
replaces only `Gg = kg A/dg` with `Gg = h A`; it does not add another gas film or
remove the liquid half-cell resistance or gray radiation. The exact coefficient
evidence and original film source remain in `FreeSurfaceInterfaceSample`.
`None` is accepted only at exact interface/gas thermal equilibrium. A failed
correlation, unsupported source, invalid coefficient or unrepresentable
conductance refuses the complete transaction.

The same borrowed provider remains attached to every cloned thermal observer,
nonlinear trial and endpoint evaluation. Existing positive pairwise surface
elimination supplies the signed liquid/gas/radiant heat and the explicit
capacity-based interval bound. The finite atmosphere stores exactly the
opposite gas heat on each original thermal substep; its new temperature, EOS
pressure and source properties feed subsequent trials. Optional pressure loading
uses this same parcel through the existing native pressure owner. Accepted
liquid, wall and gas states still publish together. The opt-in implicit residual
calls the provider at each trial endpoint for prescribed or single finite-gas
surfaces. Finite bodies, electrical and separate environmental paths retain
their original implicit exclusions.

This bounded entry admits explicit Cartesian planar faces. Height-graph and
local reconstructed surfaces refuse. A higher owner must declare the physical
plate geometry, whole-plate characteristic length and held relative gas flow;
a grid face is not a new boundary-layer entrance. The working plate adapter is
`fs-airflow::graph::thermal::medium::convection::PlanarFreeSurfaceConvection`.
The held-geometry thermal/native-flow split, isothermal-patch and held-property
assumptions still require time, mesh and model assessment. No gas momentum,
evaporation, arbitrary curved-surface law or total mechanical-energy closure is
introduced. Focused regressions in
`tests/material_free_surface_enthalpy/conduction/surface/convection.rs` exercise
the actual surface root and liquid resistance, finite-gas reciprocal energy and
step admission, exact equilibrium, late refusal with identical accepted retry,
and the original implicit endpoint residual.

## Curved D3Q19 walls with native open reservoirs

`step_with_curved_walls` also admits the original velocity/pressure open faces.
It streams the actual curved populations first, then invokes the SAME cached,
regularized open-face owner used by immersed-body coupling. The first-interior
stress and prescribed velocity/density use their native Guo force frame,
including spatially varying held forces and existing inlet profiles. At mixed
wall/open nodes, the native mask retains every curved or halfway wall-owned
population; only the open directions reconstruct.

`CurvedWallStep3::open` reports reservoir mass, linear impulse and angular
impulse TO the fluid. The shared compensated accumulator counts the temporary
open-link reflection and every actual reconstruction delta, including the
replacement of tangential populations on face interiors. These are direct
operator contributions, not a residual inferred from a whole-domain balance.
They remain distinct from curved-body reaction, other walls, held volume-force
input and periodic coordinate transport. BFL population mass drift remains
measured without renormalization. Reservoir exchange is not attributed to the
body, and is not labeled pressure work or a complete fluid-energy balance.

The candidate reservoir populations and accounts finish before the existing
single publication swap. A late reconstruction, arithmetic or cancellation
refusal preserves accepted populations and geometry. The open-face cache,
post/next population buffers and compensated boundary owner are reused; no
second collision law, reservoir law or whole-volume balance scan is added.

Four focused G0/G1/G3 native regressions pass: exact halfway recovery for both
velocity/pressure and pressure/pressure ends on all three axes and all four
collision families; rotating curved-body mass/linear/angular balance with
separate reservoir, periodic, stationary-wall, uniform and spatially varying
force contributions; wall-corner population preservation and steady
pressure-driven flow between fractional walls; and late reservoir-account
overflow without publication. The resolved 96x16x4 duct with a 14.6-cell gap
has relative axial-profile error 0.009887809049680325 under a 3% budget. A
6.6-cell gap gave 0.046973215227596746 and is outside that envelope; no global
order claim follows. The six original curved-boundary regressions also pass.
These are direct original numerical-owner checks with unused runtime/pool,
body-motion and unrelated API wrappers omitted. Full Cargo, DSR, cross-ISA
and performance gates were not run. The fixed-occupancy and prescribed-wall-
kinematics limits of the curved owner remain unchanged.


### Finite translating plane carriers on the physical D2Q9 disk

FiniteLubricationPlane2 and
DryDisk2::step_with_translating_lubrication_planes advance finite normal
AND tangential plane momentum within the original dry-disk/fluid transaction.
The supplied mass is per unit out-of-plane depth; center and surface point
are distinct. Orientation remains fixed. The guide receives only a pure
reaction couple, with zero work and no supported linear impulse.

The original UnderResolvedCylinderPair2 owner supplies the normal
resistance, cutoff and gap floor at the actual first-half-drift geometry.
The shared finite-mass elimination condenses normal and tangent resistance
independently with k/(1+k/(2m)), using one carrier mass. Both enter every
existing three-unknown dry-fluid candidate, with no later impulse kick or
extra fluid collision. Original cylinder spin resistance is unchanged.

Each accepted carrier moves its center AND plane point through both free
half-drifts. The second uses its solved endpoint velocity. Admission checks
the actual carrier/disk endpoint gap, never a prescribed-wall extrapolation.
With fixed normal, the relative gap is affine within each drift; positive
start/midpoint/endpoint gaps exclude crossing along these two trajectories.
Complete vector low-Mach admission applies before, at midpoint and after.

Actual carrier momentum and kinetic increments, authored loads/work and
the orientation-guide couple are independently checked before native
publication. The aggregate exterior-fluid/dry-disk/carrier momentum account
includes guide couples as external sources. Refusal/cancellation preserves
all accepted states. Only infallible assignments publish the already
admitted finite endpoints.

Focused native tests reconstruct three successive two-plane steps from
actual changing geometry, coefficient laws, momentum and work. They check
recoil that rescues a falsely predicted held-wall crossing, true first/last
drift crossings, complete vector speed refusal, exact retry and inactive
parity. This supplies missing-film physics, not resolved moving walls,
finite face extents, contact or a thermal destination. The inherited
implicit midpoint method is not L-stable; physical stiff-film relaxation
still needs time resolution.

## Coupled implicit internal energy of one finite atmosphere

Selecting `FreeSurfaceInterfaceSolve::implicit` on the existing finite-gas
surface or surface-convection entry couples one gas energy coordinate to the
original liquid backward-Euler residual and bounded Newton-Krylov driver.
Every trial reconstructs a private parcel through its original isochoric caloric
inversion. Its actual temperature, EOS pressure, composition and selected source
claims feed the same interface root as the trial liquid endpoint. No physical
history advances in a residual or directional-difference evaluation.

The added gas row is original-law `m_g Delta u_g + dt sum(Q_from_gas) = 0`.
It uses actual represented storage, with formation references canceled by the
original temperature-anchored caloric law. The liquid rows retain the original
`BackwardEulerStep` heat-rate and local/total energy gates. A separate union
gate compares actual liquid-plus-gas storage with independently accumulated
authored/molecular and radiant heat. Gas exchange is internal to that union.
The accepted observer validates the already prepared gas endpoint and accounts
for its surface heat once; it never applies that heat to the parcel a second
time. Native flow, source, pressure and final endpoint failures still discard
all candidates before publication.

`maximum_unknowns` includes the gas coordinate. Existing residual, Newton,
Krylov, interval, source-inversion and interface-root budgets remain operative.
`cell_energy_tolerance_j` bounds each liquid row and the gas row per accepted
interval. `GasEnergySolve` retains its separate original caloric-inversion
controls. `FreeSurfaceImplicitReport::maximum_gas_residual_j` exposes the actual
accepted gas defect; `FreeSurfaceFiniteGasExchange::implicit_energy_tolerance_j`
identifies the accumulated declared gas-row allowance within the reported
combined gas budget. Conditioning scales do not replace any raw joule gate.

This removes the explicit gas-capacity CFL restriction; it does not establish
timestep accuracy or unconditional nonlinear convergence. The original
first-order thermal/native-flow split and fixed gas mass, composition and volume
remain. Finite bodies, current, and separate environment/radiation paths remain
excluded from this implicit rung. Focused public-API regressions compare the
stiff coupled endpoint with an independent caloric/series-resistance equation,
check large formation-reference invariance and evolving gas/film source states,
close the radiating liquid-plus-gas union, and verify late-domain, cancellation
and work-budget refusal with an identical accepted retry.
