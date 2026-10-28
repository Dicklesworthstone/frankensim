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
BiCGStab/ILU owner. It does not reconstruct a separate flux from velocities.
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
energy. Implicit transport and native global conversion redistribution introduce
numerical mixing; no physical thermal-diffusion or mechanical-energy claim is
made. The focused `material_free_surface_enthalpy` tests exercise actual
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
  optimal" folklore is false on this family). Body forces under the
  non-BGK modes reuse the BGK-form Guo term (exact for unforced
  cells; disclosed approximation otherwise); combining non-BGK modes
  with the partial-saturation or free-surface operators is
  mechanically possible but unvalidated.
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
  `rheology::channel_flow` — local apparent-viscosity laws and explicit
  τ updates with floor/cap counts for cells outside the representable
  relaxation window. Both D2Q9 `update_tau` and D3Q19 `update_tau3` infer
  strain from `Pi_neq + (u F + F u)/2` in the physical Guo frame. D2Q9
  includes both `rho*g` and the local external force density. Uniform
  acceleration is not shear. The shared estimators correct three/six tensor
  entries directly, avoiding source-array construction and an extra population
  pass; the material-driven rheology paths reuse them. The force-reframing
  regressions cover zero/full strain, density and relaxation changes, and
  independent force directions; the existing analytic channel gates remain
  the flow-level checks.
- `thermal::ThermalLbm` and `thermal::gbeta_for_rayleigh` — D2Q9 flow plus
  D2Q5 temperature populations for Rayleigh-Bénard-style onset fixtures with
  fixed-temperature wall rows.
- `MaterialSlab::step_with_material_rheology` (2-D and 3-D) resolves shared
  material transport at actual temperature and local physical shear rate.
  The existing non-equilibrium stress estimators include the Guo half-force
  second-moment correction; a bounded damped local solve closes
  `tau=1/2+3*nu(T,rate)*dt/dx^2` with `rate=numerator/(tau*dt)`. The retained
  background relaxation field is only an initial iterate, never an assumed
  viscosity for the evolving material. Both initial and candidate endpoint
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
  bookkeeping, conservative carry redistribution, contact-model bracketing,
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
  scalar/SIMD twin is separately bit-locked to its test-only axial oracle. The unforced
  `CentralMoment` reference rung relaxes a full-rank centered D3Q19 monomial
  basis with independent second/higher-order rates. `ReducedCumulant` applies
  the nonlinear `C220`/`C202`/`C022` corrections to the three independent
  fourth-order rows represented by that basis. All paths reject invalid
  parameters or state before publication.
- `D3Q19_MOMENT_COLLISION_SEMANTICS_VERSION` — independent version for the
  optional moment-space basis, rate grouping, `mul_add` chains, cumulant
  projection, and deterministic back-transform. It does not refreeze the
  established BGK grid golden.
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
  transform/solve roundoff. Split higher-order relaxation changes
  nonequilibrium modes without relaxing degree-zero/one invariants.
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
- Rheology laws reject non-finite or non-positive physical parameters, and
  every update reports floor/cap counts when viscosity leaves the representable
  τ window.
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
`collide_cell3` returns typed errors rather than publishing non-finite or
non-positive cell states. Both moment-space rungs additionally refuse rates
outside `(0, 2)`, nonzero body forcing, or a numerically rank-deficient moment
transform.

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

`tests/extensions.rs` (tfz.19): lbm-101 power-law/Carreau profiles vs
analytic (0.12% with τ floor+cap ledger); lbm-102 Rayleigh–Bénard onset
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
and fail-closed inputs.

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
agreement, split-rate cumulant conservation, open/wall and moving/open rim
ownership, stationary planar and voxel mass conservation, exact one-step
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
mass/momentum conservation of every operator; the shear-wave
viscosity identity for every operator; the executed
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
central targets; cumulant equilibria above degree2 are zero. Midpoint forcing
changes the first cumulants/central moments as in85--87; finite-force flow
accuracy remains a downstream qualification, separate from conservation.

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
workers1/2/7. Its new product-state collision hash has an independent dyadic
oracle; native reproduction is pending. G1 requires >=1.8 shear-decay error
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

- D3Q19 grids remain BGK + Guo on a dense set of aligned SoA tiles. The
  selectable central-moment and `ReducedCumulant` cell operators are
  deterministic `O(Q^3)` correctness references: they are unforced and have no
  performance or high-Re stability claim. D3Q27 has the separate checked full
  collision/reference-grid surface above, pending native qualification.
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
  convergence-order, free-surface coexistence, rigid-body integration,
  high-Re stability, or FSI claim.
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
`Duct`/`BoundaryGrid3` fast paths and their frozen bit-semantics goldens
are untouched. `ThermalLbm3` adds a D3Q7 advection-diffusion temperature
population (c_s² = ¼) with anti-bounce-back fixed-temperature plates and
Boussinesq coupling into the Guo term. `update_tau3` adapts τ from the
local D3Q19 non-equilibrium second moment through the SAME
`rheology::Rheology` law, `TAU_FLOOR`/`TAU_CAP`, and `RheologyStats`
ledger as 2-D; the plate-channel profile is certified against the SHARED
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
invariant along that axis). `BoundaryGrid3::set_force_field` installs a
per-cell Guo force added to the uniform force (refused with open faces or a
collision model without body forcing; without a field every result is
bit-identical to before). Each coupling sets the buoyancy
`-beta (T - T_ref) g` in lattice units, advances `coupling_interval`
pooled steps, projects the cell mass fluxes, and takes one implicit
energy step whose storage uses the fluid's volumetric heat capacity in
every cell (solids get it as a PSEUDO capacity: it shapes the path, never
the fixed point). Convergence is judged on the geometric extrapolation
`delta r / (1 - r)` of the remaining per-coupling change (no estimate
until the sequence contracts), for the relative mass flux and for the
temperature over its range. The answer is a fresh STEADY energy solve on
the final fluxes; `coupling_residual_k` retains its distance from the last
pseudo-time iterate. A lattice speed above `max_lattice_speed` refuses as
`LatticeResolution`.

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
(stored-energy change vs `dt (source - boundary outflow)`). Unconditionally
stable for the M-matrix, so electronics-scale conductivity and capacity
contrasts impose no step limit; first order in time. Declared radiation
(`Some(&RadiationConfig)`, escape or surface exchange) enters each step as
sinks Newton-linearized about the previous step's temperature, implicit in
each face's own emission with the irradiation lagged one step; the
records carry the power radiated to the surroundings
(`TransientRecord::radiated_w`) and the closure counts it in the outflow.
`march_conjugate` takes the same argument.

| Fixture | Reference | Measured |
|---|---|---|
| Near-isothermal cube (Biot 4e-7) cooled by convective faces, 100 steps of tau/50 | backward-Euler recursion of the lumped ODE | 0.1380329637 vs 0.1380329672 |
| Same, dt halved | continuous exponential | error 2.698e-3 -> 1.351e-3 (ratio 2.00) |
| Heated plate under channel flow, 400 steps from cold | steady solve; energy budget | stored 2.616326564 J = net input; worst step closure 4e-12 J; within 2.5e-9 K of steady |
| Radiating floor plate (the Stefan–Boltzmann fixture with `rho c` 2.4e6), 30 steps to t = 1800 s | linearized backward-Euler recursion `T+ = T + dt (P - q(T)) / (C + dt q'(T))` | 379.460985216 vs 379.460985226 K; every step closes to 1e-8 of its input |
| Same, 15 vs 30 steps | RK4 solution of `C dT/dt = P - eps sigma A (T^4 - T_amb^4)` | error 1.095 -> 0.565 K (ratio 1.94) |
| Same, 200 steps of 600 s | closed-form steady state | 433.501643767 vs 433.501645484 K; 4.9999999 of 5 W radiated |

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

Verification is against the continuum LVEL model itself (the developed
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
cell sink Newton-linearized about the previous iterate
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
NSF owner. One non-periodic axis (at least four cells), periodic transverse
axes and a two-cell all-fluid collar are required. Interior stationary
obstacles remain supported. Topology binds before stepping; manufactured body
sources and moving edge walls are excluded. This supersedes the original
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
of physical dissipation. Momentum exchange closes the no-obstacle domain;
stationary internal obstacle reactions are additional wall impulses.

Focused tests in `nsf3_open_boundaries.rs` exercise analytic uniform throughflow
on all axes, separately summed mass/momentum/energy budgets, incoming-wave work,
outgoing-pulse refinement, NASA9 formation-offset invariance and source support,
atomic refusals and workers. `nsf_open_channel` is an executable throughflow and
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
It reuses the fixed-domain local constitutive solver and force-corrected
molecular stress. The resulting viscosity supplies collision tau, molecular
viscous heating and temperature-dependent liquid, finite-wall and gas-film
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
ingress, molecular force-corrected stress and shared bounded `solve_local` owner.
Each wet collision lane, including zero-mass interface ghosts, closes rate
against `tau = 1/2 + 3 nu dt/dx^2`. Those exact material coefficients drive the
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
