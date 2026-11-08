## Conducting walls in native open-flow heat transport

`ThermalOpenField3::set_wall_conduction(axis, links, controls)` attaches the existing implicit finite-wall graph to an original lower/upper external wall pair. Multiple axis pairs can coexist, including both transverse pairs of a rectangular duct. Link indices omit the selected normal coordinate and retain the original x-fast surface order. Every endpoint must face an active fluid cell; blocked wall-array slots cannot become hidden solid storage or consume sources.

Each graph prepares its reflections from actual post-collision D3Q7 populations before streaming. The original trial wall cells receive the accepted endpoints once; final joint publication retains the native flow/heat transaction. Rejected attempts alter only reusable graph scratch, which is rebuilt on retry. Internal conductive transfers enter the solid first-law account with opposite endpoint signs and do not become external heat. Prepared factors and workspace are reused; each pair has separate explicit graph resource bounds in addition to the carrier storage limit.

This enables forced throughflow with laterally conducting finite boundaries using one wall-temperature owner. Cross-axis graph connections, resolved solid depth, variable caloric properties, automatic buoyancy feedback and the full developing-channel correlation benchmark remain outside this implementation.

## Native open-flow temperature transport

`OpenBoundaryGrid3::prepare_thermal_field` binds a bounded D3Q7 scalar carrier to the original native D3Q19/D3Q27 geometry, solid mask, physical velocity and accepted time. `step_serial_with_thermal` stages the original momentum step and thermal transport with one final joint publication point. Sources use the existing equilibrium half-kicks and retained advective-flux history; optional per-cell thermal relaxation is admitted for the step. Persistent candidate arrays and histories are reused.

Every velocity/pressure face requires explicit reservoir temperature, including nominal outlets. Accepted local force-corrected normal velocity selects Danckwerts total-flux inflow or zero-conductive-flux outflow independently at each boundary link. Receipts measure actual incoming-minus-outgoing scalar populations. Stationary external walls retain independent finite cells per active fluid link; internal obstacle laws remain stateless and do not create hidden solid storage. Temperature, velocity, source, resource, cancellation and independently measured fluid/solid balance gates precede publication of populations, wall state, histories, receipts and counters.

The conserved scalar uses fixed reference volumetric heat capacity. This is passive low-Mach temperature transport, not compressible mass-specific enthalpy or automatic buoyancy feedback. Its heat receipt must not be replaced by the momentum reservoir mass receipt times temperature. The scalar component velocity envelope is additional to native Mach admission. No nonreflecting-outlet or general coupled time-order claim is made.

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

Source-backed D2Q9 thermal and free-surface rheology use the same physical
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
lattice-unit rheology retains its explicitly documented full two-dimensional
strain convention. Manufactured original-population tests cover both source
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
cover zero-time refresh and history-copy independence from rejected scratch.
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
