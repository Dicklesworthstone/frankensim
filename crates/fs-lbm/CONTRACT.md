# CONTRACT: fs-lbm

Lattice Boltzmann core with a D2Q9 BGK path, a tile-major D3Q19 BGK/Guo path,
link-wise D3Q19 wall/open-boundary topology, the lattice-scaling assistant,
and dense-grid extension scaffolding for vector forcing, local rheology,
thermal double-population fixtures, and a free-surface mass-ledger prototype.

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

`step_with_radiative_environment` adds prescribed gray diffuse radiation on
explicitly exposed liquid/gas patches. The nonlinear Stefan–Boltzmann law is
reevaluated every thermal substep, using a factored fourth-power difference
to preserve fluxes near equilibrium. Its positive secant conductance joins
gas-film and internal-contact conductances in the same monotonicity bound;
no temperature or energy clipping supplies stability. Radiant and gas-film
heat may share a physical surface and are each counted exactly once. Finite
walls, source-bound interfaces, optional molecular heat and attached electrical
sources retain their original owners. The resulting enthalpy feeds the actual
material viscosity and native VOF momentum/transport transaction.

This is a gray opaque surface with prescribed constant hemispherical emissivity,
view factor one to black isothermal surroundings, and cell temperature used as
surface temperature. It needs a resolved thermal boundary layer or a justified
lumped-cell approximation, plus time refinement for the explicit first-order
split. It does not infer emissivity, visibility, multiple reflections, spectral
transparency, gas radiation, evaporation or finite surrounding storage. Source,
geometry, budget and cancellation refusals publish no partial state. Focused
G1/G3/G4 tests compare the analytic nonlinear cooling integral and its time
refinement, combined gas/radiant energy, native viscosity/momentum response,
stiff subcycling and atomic refusal on synthetic constitutive fixtures.

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

This retains the existing pure-liquid, neglected-capillarity and first-order
thermal/flow approximation. The caller supplies currently exposed geometry at
each interval; the film is held during the thermal stage. Gas volume is an
explicit fixed-volume approximation, with no gas momentum, pressure traction,
evaporation, reaction, radiation or phase-density handoff. Five focused
synthetic regressions cover independent nonlinear NASA energy/pressure limits,
data-driven native population changes, finite wall/interface/molecular heat,
summed small-gas-capacity timestep admission, and late source/cancellation
rollback. Their addition does not establish executed native validation.

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
  the paired Grid collision/step methods add a force-free dense-D2Q9
  Noble-Torczynski partially saturated cell rung. Caller-supplied solid
  fractions and velocities blend fluid and non-equilibrium-bounce-back solid
  collisions while the receipt retains aggregate equal-and-opposite
  fluid/solid impulse, angular impulse, solid work, and coupling-weight sum.
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
  periodic y closure; density/velocity and non-equilibrium stress are
  extrapolated from the respective first-interior columns.
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
  `CentralMoment` relaxes a full-rank centered D3Q19 monomial basis with
  independent second/higher-order rates. `ReducedCumulant` applies the
  nonlinear `C220`/`C202`/`C022` corrections to the three independent
  fourth-order rows represented by that basis. Both moment families admit
  uniform and per-cell body forces through `C(f + S/2) + S/2`, using the
  original Guo source `S` and force-corrected velocity. Central moments thus
  receive `(I - Lambda/2) S` with each mode's own rate; reduced cumulants
  are evaluated in the physical frame with zero first central moments.
  All paths reject invalid parameters or state before publication.
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
  `f_post = f + (1-B) Omega_f + B Omega_s`, where
  `Omega_s[q] = f[opp(q)] - f[q] + feq[q](rho,u_s) - feq[opp(q)](rho,u)`.
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
  measured non-equilibrium stress from the first interior column. The measured
  and unmeasured open-step paths are bit-identical.
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
  measured non-equilibrium second moment from the first interior cell. Mixed
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
cell needs positive finite density, finite momentum, and `tau > 0.5`, forcing
must be exactly zero, and the moment origin must be finite. Collision is built
in a private proposal and the caller's output/scratch is published only after
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
periodic y, fluid face/first-interior columns, zero gravity/external forcing,
a positive finite outlet density, and a finite inlet speed squared below 0.03.
Every condition is checked before collision.
Dynamic wetting queries return typed errors when either `QtyAny` input is not a
velocity; for a non-finite signed speed or negative/non-finite pinning speed;
when a selected `contact_line_speed` curve declares non-velocity abscissa
dimensions; for either material-query refusal, non-dimensionless or
out-of-range angles, or a reversed hysteresis bracket. They query both
endpoints before selecting one, so an incomplete card refuses even when the
current direction would consume the endpoint that is present. The adapter
mutates neither the card nor free-surface state.
D3Q19 boundary construction additionally rejects non-4-multiple dimensions,
tile-count overflow, invalid collision parameters, moment-space collision with
nonzero force, unpaired periodic faces, non-tangential/non-finite or
outside-low-Mach wall velocities, non-finite or outside-low-Mach inlet
velocities, non-positive pressure density, open faces on more than one axis,
and body force combined with the current regularized open-face closure.
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
forcing is supported by all three collision families; the regularized open
velocity/pressure boundary closure retains its separate force-free admission.

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

Four registered leaf capsules opt out of the workspace default at module scope:
the collision pair under `src/d3q19/simd/{neon,x86}/mod.rs` and the pull-stream
pair under `src/d3q19/simd/stream/{neon,x86}/mod.rs`. Each stays below the
capsule line cap with an adjacent `SAFETY.md`. Safe façades validate tile shapes,
Rust borrows prove aliasing and lifetimes, and x86 one-shot selectors admit
private thunks only after their operation-specific feature checks. Miri selects
safe scalar twins.

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
refusal. The file also covers regularized x-face moment/stress reconstruction,
measured-path bit equivalence, and pre-step topology/forcing refusal.

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
`strongly_thickening_channel_matches_analytic_profile` adds an n = 5 channel
with 33 fluid rows: it must reach steady state before 80,000 steps and keep
worst pointwise velocity error below 3% of the analytic peak.

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
Moving/open walls and independently supplied force fields are outside this
SGS step surface; their owners are not silently rerouted.

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
with `none`, `smagorinsky`, or `wale`. Explicit model, Reynolds-number,
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
ordinary laminar arithmetic and empty closure state. These entries currently
use the declared constant molecular viscosity, not a material-viscosity field.

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

Every coupled energy step consumes `rho cp nu_t / Pr_t`, with fixed `Pr_t=0.9`,
from that step's projected flow. Returned eddy fields expose the same closure.
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
so the sum of `E` over the fluid cells changes only by the wall heat flux, to round-off, in a closed adiabatic box; mass and
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
component is refused (no mass crosses a wall), and the speed is held to the declared Mach bound. With a moving wall
`wall_heat_j` is the net energy through the wall faces (heat plus work), and `viscous_heating_j`, defined as the remainder of
the internal-energy budget, excludes that work: in a steady Couette flow the net wall energy vanishes (work in equals heat out),
so the ledger carries no heating information there and the physical check is the temperature profile
(`couette_viscous_heating_profile`).

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

**Ledger.** `EnergyLedger` reports internal, kinetic, pressure-work, viscous-heating, wall-heat, boundary-energy and
volume-source terms in joules. `source_energy_j` is the exact trapezoid of the same `S_E` samples applied to fluid cells;
solid cells and callback reevaluation cannot add an audit-only supply. Negative supplies remove energy. The total-energy
residual subtracts wall, boundary and volume supplies independently; no internal-exchange remainder is a balancing credit.
Cumulative source energy is published only with an accepted step and uses the medium's fixed sensible-energy chart.

`pressure_work_j = -integral(p div u) dV dt` is the reversible exchange (trapezoid in time, fourth-order divergence).
For source-free and heat-only steps, `viscous_heating_j` is DEFINED as
`internal change - wall heat - pressure work - source energy`, with the existing boundary terms where applicable.
With no source, it equals the kinetic-energy loss net of pressure work term by term to round-off.
THAT IDENTITY IS BOOKKEEPING. The physical content (the heating is continuum dissipation) is the Taylor-Green decay rate.

Manufactured mass/momentum sources also alter kinetic storage, and this owner has no separately qualified local
source-work split for its streamed mass source. Such steps report the signed residual in
`unresolved_internal_exchange_j` and zero in `viscous_heating_j`; the unresolved term includes both dissipation and
unseparated source mechanics and must not be deposited as heat. Uniform acceleration can give a negative unresolved
term while the exact supplied energy becomes kinetic energy. This does not weaken the independent total-energy audit.
Source-free arithmetic, dynamics and bit-semantics version are unchanged.

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
Smagorinsky or WALE closure on the existing low-Mach NSF owner. Only the
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
resistance, with zero modeled wall heat flux. Characteristic and prescribed
reservoir faces refuse active SGS in either configuration order. The existing
Mach, temperature and total-viscosity envelope remains in force; ratio-cap
violations refuse rather than clip. Zero coefficient preserves molecular
arithmetic, and acceptance-render use remains forbidden.

`sgs_transport` exposes accepted-cell eddy viscosity and conductivity. Focused
tests cover moving-wall gradients, independent bulk relaxation, conservative
heat diffusion, closed energy, worker identity, exact zero control, fixed-Re
closure refinement, NASA9 heat capacity and admission failures. These are
algebraic low-Mach transport integrations, not a new compressible-turbulence
qualification: no SGS kinetic-energy equation, isotropic SGS pressure,
pressure-dilatation, shocklet closure or thermal wall function is supplied.
The modeled shear and Prandtl heat-flux form follows NASA report 20200002378,
equations (8)-(10), without its additional compressible terms.

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
The original material owner reconstructs density/energy with its stable
(u_NASA9 - u_offset) - u_s residual; no new polynomial evaluator or caloric
integrator participates. `energy_reference_j_kg()` exposes the offset.
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
EOS/transport feedback, formation-energy-offset invariance during density
changes, data-only conductivity substitution, worker invariance, local-source
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
  first rung is force-free dense D2Q9 and selects the published
  non-equilibrium-bounce-back variant. It has no stresslet-accuracy,
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
  condition, accepts no body force, and makes no unbounded-domain or blockage
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
  low-Mach fixture, not a high-Re characteristic/non-reflecting boundary; it
  accepts only constant per-face targets and refuses simultaneous Guo force.
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
  non-equilibrium rescaling and a FIRST-ORDER interface handoff: measured
  2.5% steady Poiseuille deviation at the level jump and 5.6% extra shear-
  wave decay rate — honesty-labeled in lbm-108; the post-collision
  (Filippova–Hänel-style) second-order transfer, general octree topologies,
  and dwr-adaptivity-driven refinement signals are recorded successors.
- Dense contact-line fixtures remain MODEL-BRACKETED (neutral vs wetting
  fill ghosts; lbm-107 reports the sensitivity band). The opt-in sparse 3-D
  planar equilibrium-angle closure retains its separate, bounded convention
  and does not claim calibrated dynamic wetting.
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

### No-claim boundaries (3-D coupled ports)

- Single-threaded correctness-first reference rung: no SIMD/tile/parallel
  claim, no performance claim against the tiled fast path.
- No turbulence model and no claim beyond the fixture-scale onset band
  and profile bars actually gated in the battery.
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
  fluid/solid or solid/solid contrast) and Patankar power-law (or upwind)
  convection in conservative form. `ThermalFace` rules: `Adiabatic`,
  `Temperature`, `HeatFlux` (into the domain), `Convective { h, ambient }`,
  `Inflow { temperature }`, `Outflow { backflow_temperature }`; solid cells
  on open faces are adiabatic. Sources via `ThermalSetup::power_w` (W per
  cell, `add_uniform_power`) and fixed-temperature cells. The M-matrix is
  solved by ILU(0)-right-preconditioned BiCGStab on the Jacobi-scaled
  system; acceptance is the RECOMPUTED true residual.
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

Both field APIs support BGK, central moments, and reduced cumulants and refuse
open velocity/pressure faces. All fluid samples and total forces must be
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

### Transient conjugate conduction (`march_energy`)

Backward-Euler steps of `C dT/dt + A T = b(t)` over a FROZEN flow field
(forced convection), with `C = (rho c) V` from the fluid and each solid's
declared `SolidMaterial::with_heat_capacity` (a solid without one refuses),
sources `power_w * power_schedule(t)`, and per-step records of the peak
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
  `FvBoundary::{Wall { velocity }, Symmetry, Inlet { velocity }, Outlet}`;
  an outlet face is a momentum unknown whose outside twin mirrors the inside
  cell (zero normal gradient) and carries the ghost pressure `-p_inside`, so
  the face holds pressure zero and the pressure correction uses the matching
  linearization `u' = 2 d p'_inside`. The steady state is one fixed point of
  the discrete equations, independent of the relaxation, also for
  undeveloped (even partly reversed) outflow. (The first version
  re-extrapolated outlet velocities each iteration and stalled at a 3e-4
  mass residual on the heatsink, whose outflow is not developed.)
- Momentum uses power-law (or upwind) coefficients with `a_P / alpha`
  under-relaxation (alpha in (0, 1); SIMPLEC's `d` degenerates at 1). Closed
  components of the pressure correction are pinned.
- Convergence is two steady residuals of the same iterate: the largest cell
  mass imbalance before correction over the largest face mass flux, and
  per component `||b - A u||` of the Jacobi-scaled momentum system at the
  velocities entering the iteration, over the flow's momentum scale
  `sqrt(rows) max_c rms(b_c)` (a component that vanishes by symmetry has a
  round-off right-hand side; normalizing it by its own norm left plug flows
  through a porous block stalled at a "residual" of 0.3). Inner solves only reduce their
  entry residuals (by `momentum_tolerance`, default 0.1, and
  `pressure_tolerance`, default 0.01), so a converged report never rests on a
  skipped inner solve; one tight final correction (1e-8 of the converged
  imbalance) then hands the energy equation a projected flux field. Budget
  exhaustion refuses as `FlowNotSteady`; non-finite velocities or pressures
  after any momentum solve or correction refuse as `FlowDiverged` before they
  reach an incomplete factorization.
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
coarse voxel meshes: per fluid cell, the wall distance `L` (exact Euclidean
transform to the nearest solid voxel or wall-type domain face; openings,
inlets and symmetry planes do not count) and the local speed give
`Re_L = V L / nu`, which Spalding's unified law of the wall (`k` 0.41, `E`
8.6) inverts for `u+`. Interior diffusion uses the profile's tangent
`nu_eff = nu dy+/du+`; the half-cell flux to a wall uses its secant
`nu_w = nu y+/u+`, so the discrete wall shear `mu_w u_P / L` equals
`rho u_tau^2` of the law of the wall at any first-cell `y+` (using the
tangent at the wall over-predicted the friction by 61 % and the Nusselt
number by 43 % at Re_Dh 2e4). Both are under-relaxed by one half per sweep.
`FvFlow::eddy_conductivity` gives the energy equation `rho c_p nu_t / 0.9`
(the secant value in wall-adjacent cells); `fv_natural_convection` applies
the same closure. The cooling-cht scene selects it with
`solver.turbulence: "lvel"`.

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
| Refusals | face fan, empty march, compact components | structured errors |

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

The lattice-Boltzmann paths (`lbm_duct_flow`, `natural_convection`, the
cell-velocity projection) need uniform voxels and refuse a graded domain;
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
temperature. The energy solve (steady, and each transient step) iterates
Picard passes, each evaluating `k` at the previous pass and warm-starting
the Krylov solve, until the largest temperature change is below 1e-10 of
the span (200 passes, else `SolverNotConverged { system: "energy
conductivity" }`). Face conductances stay harmonic in the cell values.

| Fixture | Reference | Measured |
|---|---|---|
| Slab, `k = 1 + (T - 300)/100`, faces 400 K / 300 K, 10 and 20 cells | Kirchhoff transform: `int k dT` linear in x | max error 0.229 K -> 0.063 K of a 100 K span: order 1.86 |

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
the domain boundary; power, sinks or fixed temperatures inside a box; any
transient march (the model has no heat capacity).

| Fixture | Reference | Measured |
|---|---|---|
| 4 x 4 x 2-cell component (R_jc 50, R_jb 20 K/W, 1 W) between k = 5 slabs held at 300 K / 310 K | two-path network `R_b = R_jb + 2 dx/(kA)`, `R_t = R_jc + 2 dx/(kA)` | junction temperature and both resistor heats within 1e-9; balance 8e-12 |

No-claims: the two-resistor model's own limits (boundary-condition
dependence of a 2R reduction, no spreading inside the package, adiabatic
sides); steady only.

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
average the two adjacent cells). The coupling starts from conduction
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

### Sourced equilibrium angle at one planar domain wall

`MaterialFreeSurfaceScaling::with_planar_wetting` binds a
`ResolvedLiquidWettingState` to `PlanarWall3` before constructing flow.
It requires the same liquid phase/process/revision, complete query and quantity
descriptors as the existing surface bundle; only the angle comes from this
ordered substrate/liquid interface. The original source/scaling owner continues
to supply rho, viscosity and sigma. The caller assigns substrate and named
environment to the receiving geometry; this does not resolve a gas EOS.

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


### Tetrahedral liquid geometry and momentum enter native VOF

`material::free_surface::transfer` connects the same resolved liquid scaling to
`fs_mesh::overlap` and the native supplied-state initializer. It intersects each
source tetrahedron with actual Cartesian cells, integrates its P1 velocity field,
and uses the existing `fs_la` LU/Cholesky owners to reconstruct momentum.
Mass-weighted translation and rotation corrections preserve total linear and
angular momentum at receiving cell centres, subject to explicit error bounds.
The final checks read actual native fill, density and Guo-corrected velocity.

Source mass, first mass moment, linear/angular momentum and kinetic energy use
analytic simplex integrals. Signed kinetic and uniform-gravity potential changes
from projection are separately reported numerical errors, never physical heat.
Source coverage, target overfill, mass/momentum, centroid displacement, projection
energy, native speed, work and occupied-cell caps all gate publication. Unsupported
spin on a rank-deficient receiving grid requires refinement. Original Cx polls
cover geometry, integration, correction, initialization and final publication.

This is an already-liquid, stress-free, uniform-material isothermal transfer.
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
adjacent live cell. Entropy and tangential velocity use the inflow reference or
the actual outflow donor. The reference density, impedance and caloric tangent
come from the SAME retained `NsfGasMedium`. Reconstructed pressure and temperature
query its actual NASA9 caloric energy and local transport, or the original
GasSpec path; no independent reservoir gamma, density, sound speed, viscosity
or conductivity is accepted. The declared bound on pressure, temperature and
velocity perturbations must lie in (0, 0.1], with the original material support,
Mach, temperature-window and viscosity refusals still active. Composition stays
frozen. The entropy/acoustic boundary reconstruction is explicitly linearized
about the authored reference; NASA9 energy inversion inside the grid remains
nonlinear and source-owned.

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

## Isothermal source rheology with native capillarity

`material::free_surface::rheology::IsothermalRheologyLaw` binds one immutable
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

This is isothermal generalized-Newtonian response on stationary walls, with the
explicit numerical contact ghost inherited from the domain. It does not add
thermal storage, viscous heat, viscoelastic or yield history, moving walls,
variable-sigma traction, phase change or conservative surface-energy evolution.
Source and numerical admission are not stability or experimental-accuracy
certificates. Focused tests in `tests/material_free_surface_rheology.rs` exercise
an independent force-corrected tensor and analytical constitutive root, native
Newtonian bit parity, source-only viscosity/capillary motion changes with mass
conservation and deterministic replay, and endpoint source/convergence rollback.

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

MR28 remains open. This handoff requires already-liquid, stress-free material,
one constant-density caloric chart and explicitly omitted capillarity. It does
not discard solid elastic/plastic history, implement melting release, change
phase volume, supply pressure/capillary work or certify total fluid mechanical
energy conservation.


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

The coupled entry admits only force-free BGK fluid/wall domains. It checks
positive finite inertia, the matrix response, and initial/midpoint/endpoint
body point speeds before publishing body state or caller scratch. The
collision exchange conserves fluid-plus-body linear momentum and total
orbital-plus-spin angular momentum to floating-point roundoff. Streaming
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
scratch, which the next call rebuilds. Fixed radius/inertia, one disk, force-free
BGK and lattice `dx = dt = 1` are the supported scope. Multiple bodies, wall
contact, periodic body seam crossing and free surfaces require further owners.
Native outer-wall reactions and periodic fluid angular flux remain distinct
from immersed exchange. No global temporal order, physical drag calibration or
complete moving-boundary qualification is inferred.

`tests/moving_disk.rs` checks analytic translated area, common translation over
many actual steps without stale-cell coupling, independent fluid-plus-body
linear/orbital/spin balances before outer-boundary transport, bitwise replay,
and late-endpoint rollback with scratch recovery. The broad moving-boundary
bead remains open.

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
independent `FREE_SURFACE3_BIT_SEMANTICS_VERSION` is 2 for the new default,
leaving the unchanged D3Q19 core, boundary and sparse goldens untouched.

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
surface reservoirs and declared liquid-liquid conduction links, with the same
material, caloric, gas-film and radiation laws. Finite gas, finite bodies,
electrical coupling and the separate gas/radiation APIs require additional
coupled residuals and refuse this initial implicit entry.

Unknowns are specific-enthalpy increments in J/kg, mapped to actual cell energy
increments with the original masses. A small cut cell therefore does not receive
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
refinement and does not claim elastic memory, capillary energy, phase-conductivity
multiplier or buoyancy transfer. Unsupported owners refuse explicitly.

The integrated regression melts the specimen, drives native circulation,
releases it and advances native free-surface motion with data-only rheology
variants. It independently checks molecular heat against the volume sum of
mu*gamma_dot^2*dt and face conduction against the same source k(T,gamma_dot),
alongside the original transfer budgets, clock and refusal rollback.


## Spatial finite-inertia immersed-body coupling on D3Q19

`BoundaryGrid3::step_with_partial_saturation_body` couples the three translational
and three angular velocities of `PartialSaturationBody3` to the native BGK
populations and then runs the existing pull stream. The supplied symmetric,
positive-definite inertia is expressed in world axes about the body's center;
off-diagonal products and asymmetric sampled fractions are retained. Fractions
are row-major `(z*ny+y)*nx+x`; cell centers are `[x+0.5,y+0.5,z+0.5]`.

At offset `r`, `Gv = V + omega cross r`. The same Noble-Torczynski
non-equilibrium reflection law as D2Q9 supplies solid impulse
`B*rho*(u-Gv)`, with `B=epsilon*(tau-1/2)/(1-epsilon+tau-1/2)`.
The operator reduces `K=sum(B*rho*G^T*G)` through a scalar weight, first
moments and the six independent rotary moments; one deterministic 6-by-6
Cholesky solve advances `(M+K/2)*delta_v=sum(B*rho*G^T*(u-Gv_n))+J_external`.
It never feeds an explicit fluid load back through a tiny inverse body mass.
The native collision runs once at the solved midpoint velocity. No grid-sized
matrix, per-step field allocation, nonlinear coupling iteration, or duplicated
fluid owner is introduced. Existing post/next buffers are reused; first-use
next storage and fixed 3-by-3/6-by-6 Cholesky factors allocate.

The receipt measures the actual population correction to ordinary BGK,
equal-and-opposite body impulse and angular impulse about a caller-selected
origin, center torque reduced without distant-origin subtraction, and midpoint
body work. Actual rounded velocity changes independently determine body
kinetic-energy change and a six-component momentum residual. External linear
impulse acts at the center; the final three external components are a couple.
All stationary exterior/voxel wall exchanges are separately measured from the
blended populations using the existing native link receipt. Thus the immersed
body's action/reaction is not confused with wall reaction. Global angular
balance on a periodic domain additionally needs periodic transport flux.

This is a BGK/Guo, fixed-geometry/world-inertia velocity substep in
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

This is one BGK/Guo immersed body in lattice `dx=dt=1`, with solid inertia
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


## Original phase conductivity through moving-liquid release

The original `MaterialFreeSurfaceEnthalpy3` can retain `LiquidMassFractionConductivity`. Both ordinary and rheological complete-melt release transfer that exact law before publishing the receiving liquid and its clock. One shared helper keeps the pinned raw source, declared phase multiplier and effective conductivity separate. Liquid/liquid, finite-solid contact, gas-film and solved gas/radiation-interface half-cell resistances all consume that same effective coefficient, including thermal CFL admission. Source overflow and cancellation refuse before native publication. Focused regressions independently close each heat balance, preserve identity-law behavior and source claims, and carry the law through real melt/release and viscosity-dependent motion. This owner still admits only pure liquid at fixed caloric density; moving mushy phase, solidification, shrinkage and surface-creation energy remain outside this feature.

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
`sum_a B_a (Omega_s,a(g) - Omega_BGK(g))`, added to the unchanged native
`collide_cell3(f, BGK, F)`. Its first moment remains
`sum_a B_a*rho*(v_a_mid-u)`; the original six-velocity solve uses this same
physical velocity. The remaining imposed collision impulse is `F`. This is
the declared additive-force extension of the retained-fluid PSM law, not a
new displaced-fluid density closure or a full coupled temporal-order claim.
Unforced paths retain their numerical arithmetic, and zero occupancy reproduces
the native forced BGK populations bit for bit.

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
