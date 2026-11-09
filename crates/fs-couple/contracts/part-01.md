<!-- fs-couple contract, part 1 of 8; a higher number is newer. Index: ../CONTRACT.md. -->
## Generalized-Newtonian film resistance in impact mechanics

`render::plate::impact::squeeze::power_law` binds the native
`fs_tribo::resistive_film::power_law::PowerLawFilm` to the existing
`ImpactSystem` force, directional tangent, and endpoint-acceptance paths.
`with_power_law_film` accepts explicit SI channel laws and bounded nonlinear
flow/work controls. Read-only count and indexed access preserve the complete
constitutive owner; its insertion indices are independent of normal and
sliding film indices.

Both reference and prepared execution evaluate the full shear-dependent
pressure law at the original configuration and actual discrete-gradient
velocity. They add its work-conjugate force once through the original
nonlinear dissipation callbacks; the analytic Newton image uses the same
native law's directional tangent. The original time owner advances motion,
retains material history, and accounts physical loss. Constitutive laws and
their declared solver budgets remain immutable during mechanical trials.

Normal-pressure, sliding, and generalized-Newtonian films share the existing
four-island ceiling through one internal count. All attachment orders enforce
this work bound and compressible-gas exclusion. Mechanical-mode dimensions,
mapped support-memory restrictions, initial zero-rate domain admission, and
final zero-rate endpoint admission are retained. The new field is initialized
in the single mechanical constructor, and later relaxation/radiation memory
attachment preserves its original mode coordinates and film ownership.

The native owner refuses sealed pockets, invalid geometry, excessive
pressure, unrepresentable terms, or failed flow/work convergence before
the mechanical transaction publishes. Original cancellation and mechanical
force/energy gates continue to govern accepted state and time. This is a
stationary-wall, full-film reduction with prescribed laws. Physical film work
can be routed into finite material enthalpy as described below; temperature
feedback, cavitation, inertia, and non-Newtonian sliding remain separate.

### Generalized-Newtonian film work in finite material enthalpy

`HeatedImpactSystem::new_with_power_law_losses` adds explicit
`PowerLawFilmHeatRoute` destinations for a complete original power-law pressure
island or selected original channels. Squeeze, sliding, and power-law insertion
indices remain distinct. Duplicate channels, whole-island/channel overlap,
absent constitutive owners, and invalid thermal destinations or partitions
refuse before stepping.

Each selected island uses its original nonlinear evaluator once during heat
staging, at the same midpoint geometry and actual discrete-gradient velocity
as mechanics. Physical heat is `dt * Q * delta_p`; the separate numerical work
defect is not deposited. Original iteration/work controls and bounded
cancellation checks remain active. The common finite-enthalpy, coupled-transport,
and thermal-network paths retain one accepted mechanical/thermal clock and
publish all owners together.

Unselected work retains its mechanical-loss ownership. Selected heat transfers
internally once and is deducted once from remaining mechanical loss. The heat
route alone holds rheology coefficients fixed; temperature feedback requires
the explicit material binding below. It does not infer wall/fluid partition,
fluid inventory, sliding entrainment or cavitation.

### Finite-temperature feedback into the original power-law pressure law

`LiquidFilmImpactSystem::new_with_power_law` adds complete canonical per-channel
`LiquidPowerLawFilmChannelBinding` coverage in the original power-law island
namespace. It explicitly selects `K(T)=K_ref*mu_source(T)/mu_source(T_ref)`, with
positive caller-supplied reference consistency in Pa s^n, pinned
`LiquidTemperatureLaw` claims/context, and the unchanged native flow index and
physical core shear rate. The original viscosity supplies a dimensionless
ratio; it is not reinterpreted as power-law consistency. This separable
thermo-rheological approximation requires caller justification. Sourced
viscosity alone does not establish K_ref or temperature-independent n/core.

Each original finite liquid store supplies the selected channel's temperature.
The existing material-card, fully-liquid, constant-density/mass/volume and
per-sample relative-source-viscosity gates remain authoritative. Every channel
of a selected original island must be covered once; distinct passages retain
their own source law and temperature without averaging constitutive values.
Power-law, sliding and normal film insertion indices stay independent.

The native `PowerLawFilm` compiles the proposed consistency and core pressure
operator, retaining geometry, domains, solver budgets, nonlinear force and
directional tangent. A completed mechanical sample charges its original
physical Q*delta_p once at the previously accepted coefficients. All proposed
sources and operators stage inside the existing contact/enthalpy/mechanical
transaction before infallible accepted-boundary swaps. Failure or cancellation
preserves motion, heat, clock, laws and source receipts. Thermal-network
conduction uses the same transaction and can change local resistance.
`into_power_law_parts` retains fixed reference coefficients so rebinding an
already warm system does not compound the temperature shift.

Six focused regression groups cover independent shear-profile work for
thinning/Newtonian/thickening laws with both Newton images, local channel
flow/heat, canonical binding order and warm rebind, source/thermal/step-change
refusals with exact retry, and network-conduction feedback. Source parsing and
scoped whitespace checks passed; native compilation and execution remain
unverified. Temperature-dependent n/core, non-Newtonian Barus pressure coupling,
sliding power-law films and a spatial gap-temperature PDE are not supplied.

Four focused regression groups cover independently integrated thinning,
Newtonian and thickening shear work, finite caloric and momentum balances,
selected-channel ownership and canonical ordering, thermal-carrier continuation,
and late-cancellation/domain rollback. New files passed source formatting and
scoped whitespace checks before an executor outage. Recovered manual patches
were reviewed and read back through GitHub; native compilation and execution
remain unverified.

Focused regressions exercise an above-core thinning law against independently
integrated stress profiles and a scalar discrete-gradient step, reference
versus prepared finite-difference/analytic execution, directional derivatives,
work and endpoint/attachment/refusal behavior. Exact coupling sources parse;
a bounded independent owner review covers wiring and state preservation.
The underlying exact Reynolds sources and 61 native film plus seven real
math test definitions pass native type and borrow checking. An independent
98-case numerical profile check bounds flux relative error below 6.3e-15 and
pressure/aperture derivative error below 3.4e-10. Coupling runtime tests and a
full crate build remain unexecuted in the full-storage environment.

<!-- porous-covariance-feedback:start -->
### Partially coherent pressure spectra on the native thermal trajectory

`porous_heating::feedback::covariance::PorousCrossSpectrumHeating` binds one
admitted `PressureCovariance` to each original physical frequency. Each
`PorousCovarianceBin` has a fixed positive real exterior reference and the same
full specimen area. Distinct frequencies remain canonical and statistically
independent. Duplicate frequencies refuse; independent same-frequency emitters
must first add covariance matrices, while coherent emitters require their
phasor sum or the appropriate cross correlation. Components never masquerade
as extra frequencies.

The covariance is bin-integrated RMS pressure squared. Native coherent peak
amplitudes have `C = p p*/2`; the angular-density constructor applies a positive
`d omega` quadrature weight exactly once. One changed cell evaluates its native
gas state once across the spectrum. Each affected frequency then rebuilds its
changed porous coefficients and composes one stack before at most two original
coherent excitations. Both statistical components share the same physical
layers, interface map and admitted source state.

Positive native viscous/thermal cell work adds before one accepted-start
thermal step. Boundary deficits never manufacture heat. The same `PreparedHeat`
and `HeatAccount` owners retain total and per-frequency exterior-energy
accounts without duplicate credit. Native caloric identities, thermal masses,
geometry, conduction, authored heat, events, latent storage and continuation
remain bound; late domain/cancellation refusal publishes neither candidate.
Existing frequency/cell-product limits apply, with at most two component fields
per cached cell/frequency pair and no second stack composition.

This retains stationary ensemble/long-time averaging, prescribed isobaric gas,
rigid-frame normal incidence and justified common gas/skeleton temperature.
Finite-time beats, stochastic thermal fluctuations, frequency conversion and
transient acoustic storage remain outside the model. Four integration groups
cover exact coherent-spectrum parity, uncorrelated illumination, zero drive,
shared gas/composition work, latent heat, events/continuation and atomic
refusals, plus independent direct `H C H*` spatial work and native temperature/
energy refinement against a two-frequency RK4 reference.
<!-- porous-covariance-feedback:end -->

<!-- porous-oblique-feedback:start -->
### Oblique porous absorption on the shared native thermal trajectory

`TemperaturePorousHeating`, `PorousSpectrumHeating` and
`PorousCrossSpectrumHeating` expose opt-in `try_new_oblique` constructors.
A single drive supplies one conserved tangential wavenumber [rad/m]; a spectrum
supplies one per INPUT frequency bin. Frequency canonicalization carries each
wavenumber with its original drive/covariance, and duplicate frequencies still
refuse. Public drive/bin literals retain their existing shape. Tangential
wavenumbers are finite and nonnegative, including for zero incident power, and
all original bin, cell-product, composition and cancellation budgets apply.

Every frequency uses the original `PorousLayer::try_new_oblique` and existing
stack/positive-work owners. Both normal and tangential viscous work reach the
original native heat source; pressure-relaxation heat and the harmonic factor
are unchanged. Covariance components at a frequency share the same tangential
mode and one composed stack. Changed-cell gas evaluation is still shared across
frequencies, and accepted-start deposition and endpoint refresh retain the
original thermal timeline, per-bin diagnostics and atomic candidate semantics.

The exterior tangential wavenumber and positive real NORMAL wave reference stay
fixed during internal temperature feedback. The caller supplies a physical
propagating exterior mode; no angle is recomputed from a changing internal gas
speed. Full specimen area and normal thickness are unchanged. Source support,
rigid isotropic frame, prescribed isobaric gas, common gas/skeleton temperature
and original finite thermal mass/caloric identities remain binding. This does
not infer angular diffusion, finite-face diffraction, anisotropy, poroelastic
frame motion, exterior evanescence or transient acoustic storage.

The native spectrum tests retain exact zero-tangent arithmetic, per-bin angle
association, coherent/covariance equivalence, independent two-port RMS work,
shared refresh counts and atomic refusal. The independent
`porous_oblique_heating/refinement.rs` witness reconstructs normal and transverse
velocities from original constitutive coefficients and interface continuity,
checks separate volume losses against normal exterior flux, and refines the
native temperatures, enthalpy and per-bin/local energy against independent RK4
for two distinct frequency/angle bins. Its transverse viscous term is material,
so replacing vector work with the normal component alone fails the witness.
<!-- porous-oblique-feedback:end -->

<!-- porous-angular-feedback:start -->
### Independent angular ensembles and diffuse porous heating

`PorousCrossSpectrumHeating::try_new_independent_modes` accepts explicit
`angular::PorousAngularMode` entries with zero cross-covariance between entries.
Distinct tangential wavenumbers may share a physical frequency; each mode keeps
its original two-face Hermitian pressure covariance. Canonical order is
`(omega, nonnegative k_t)`. Duplicate keys refuse: independent emitters with the
same key, including isotropic equal-angle azimuths, must first combine their
covariances. The existing frequency-only constructors continue to reject any
repeated frequency. No coherent cross-angle interference is silently discarded.

The original oblique layers, one stack solve per mode, positive normal plus
transverse viscous work, pressure-relaxation work, and compensated native heat
sum are unchanged. The shared gas still evaluates once per changed thermal
cell, independent of the mode count. Original `bin_power`, `bin_cells` and
`bin_energy` now refer to canonical mode entries for this opt-in constructor;
per-mode accounts remain diagnostics, not extra deposited energy. Existing
mode/cell-product budgets, accepted-start rates, endpoint refresh, conduction,
latent heat, events, continuation and transactional refusal remain native.

`angular::DiffuseHemisphere` binds a checked original exterior acoustic state,
one frequency/band, full face area and independent left/right normal incoming
irradiances E [W/m2]. The exterior is modeled as lossless. Positive quadrature
nodes use `mu=cos(theta)` and weights w in `dmu`. Full-hemisphere weights must
integrate both 1 and `2*mu` to unity within 64 binary64 epsilons, checked with
compensated sums. The builder applies no empirical grazing cutoff, hidden
azimuth factor or renormalization. The caller refines angular quadrature.

An isotropic hemisphere has normal flux measure `2*mu*dmu`: a mode receives
`area*E*2*mu*w` watts and normal reference `R=Z0/(area*mu)`. RMS pressure
variance is formed directly as `2*Z0*E*w`, cancelling the grazing cosine before
arithmetic. Thus pressure variance stays finite as incidence approaches grazing.
Exact grazing, invalid moments/weights, nonfinite inputs, unrepresentable normal
references and numerically lost positive variances refuse, including zero-drive
geometry checks. Irradiance is already integrated over its frequency band; no
extra bandwidth, peak/RMS factor or `2*pi` is inferred. The external mode and
reference stay fixed while the internal pore gas changes temperature.

This is a stationary independent plane-wave ensemble on laterally homogeneous
isotropic rigid-frame layers. It does not certify a finite reverberant room,
finite-face diffraction, nearby-source absorption, coherent angular speckle,
anisotropy, frame motion or fluctuating temperature. The random-incidence
plane-wave definition is discussed by Sgard, Atalla and Robin (2024), section
2.2.3, [doi:10.3389/facou.2024.1414356](https://doi.org/10.3389/facou.2024.1414356).

Focused tests retain exact normal-mode/native parity, frequency and mode
association, duplicate refusal, hemispherical irradiance, finite near-grazing
pressure, shared refresh work, continuation and atomic cancellation. An
independent original-rho/K transfer-matrix and vector-velocity volume integral
checks per-angle R/T and separate viscous/thermal losses, angular quadrature
refinement, and two-frequency many-angle thermal/energy convergence against RK4.
<!-- porous-angular-feedback:end -->

### Heterogeneous layered moisture and compatible solid response

`moisture::solid::MoistureLayeredSolid` binds one original ordered
`MoistureLayeredSlab` to a conforming prism mesh with six native tetrahedra
per actual moisture cell. Shared interface nodes impose perfect mechanical
bonding independently of the stack's water-transfer contact or sealing laws.
Each layer retains its own pinned elastic and hygroscopic sources, original
moisture datum, principal directions, dry inventory and local accepted
moisture. No material, moisture or dry capacity is averaged across a join.

The homogeneous and layered consumers share the same private source/cell
evaluation and native FE equilibrium owner. One globally compatible solve
returns displacement, stress, reactions, elastic energy, admitted deformed
geometry, conserved wet inertia and the signed initial-stress vibration
pencil. Cell/source arrays remain in physical layer order, with explicit
layer-cell ranges.

Fixed-duration and finite-vapor candidates publish water and mechanics
together. Adaptive wrappers reuse the original layered water controller and
project its exact final endpoint into native mechanics; late source,
equilibrium, geometry, mass or cancellation refusal discards the entire
proposed window. Temporal estimates concern moisture/activity and reservoir
water, not mechanical interpolation error.

This remains fixed-reference, isothermal, quasistatic small-strain mechanics.
Mechanical bonding does not imply ideal water contact. Deformation does not
silently change diffusion geometry; stress-assisted transport, sorption heat,
hysteresis, finite deformation and dynamic mode-state remapping remain
separate models. Authored G1/G3 fixtures cover homogeneous splitting,
heterogeneous rotated free swelling, constrained stress/energy, native
generalized eigenvalue changes, finite-vapor coupling and whole-window
atomic refusal; they are not measured laminate or wood validation.

### Quasi-steady interference feedback in one original film

With the existing layered-wave-heating feature,
layered_wave_heating::interference::QuasiSteadyFilmHeating rebuilds a selected
Debye, Ohmic or DebyeWithOhmic film at the accepted native thermal temperature.
It retains the original multilayer order, thicknesses, incident and exit source
objects, periodic transverse grid and authored physical depth quadrature.
Only the selected homogeneous film changes state. Its original complete material
bundle is refreshed through the existing same-card query owner, including
permeability; the combined Debye/Ohmic law uses one identical refreshed bundle.
Fixed source conditions and typed absolute-temperature conventions must agree
with retained thermal source receipts. No constitutive data are fitted or inferred.

Both this owner and `QuasiSteadySpatialInterferenceHeating` expose
`try_new_with_carrier` for one authored transverse wavevector in rad/m.
Both face arrays, every frequency and every statistical factor describe periodic
tangential-electric envelopes with the same `exp(-i k_carrier.x)` removed.
The shared angular-spectrum solver preserves that wavevector through each
frequency solve and material-temperature rebuild. Incidence angle therefore
follows the original exterior wavelength; this is not fixed-angle broadband
illumination or an independent carrier for each line. Existing constructors
retain zero carrier. Continuous cell/rectangle heat keeps the actual spatial
sum/difference interference; sampled heat quadrature still requires its own
spatial resolution. The original two-port light cones, causal reservoir heat,
thermal energy and absolute DD clock remain authoritative.

The thermal receiver is the original LumpedEnthalpyBody. Its complete chart
must be fully solid with fixed density; its material card must match the film.
The declared quadrature covers the complete film thickness and original full
periodic cell. Partial masks are refused. Thermal mass divided by density must
match both geometric and actually weighted quadrature volume within 64*f64::EPSILON
relative allowance. Weights remain unchanged. Native Biot admission remains
mandatory; it does not alone establish internal equilibration under nonuniform
standing-wave/beam heating. The uniform-body approximation additionally requires
internal thermal equilibration faster than the relevant thermal-envelope evolution.

ThermalInterferenceDrive retains both physical-face peak E arrays on one source
clock for every line. The optional nonempty deterministic mean and each explicitly
real/proper factor preserve the original joint second-moment assumptions. All
frequencies in a factor remain coherent. Frequency caches are rebuilt once per
distinct frequency from the refreshed original film law; original angular-spectrum,
two-port, exposure, reservoir and integrated-heat owners supply all field physics.
The output spectrum is already W and reaches the native enthalpy owner at scale 1.
Reversible field storage, polarization energy and parameter work are never
reclassified as thermal dissipation.

advance is an EXPLICIT quasi-steady split. It holds the initial-temperature
coefficients over the requested feedback interval, even if native thermal/ripple
controls subdivide it. The accepted clock is DD and is never reset by a refresh.
The native solver keeps its analytic interval source joules, nonlinear phase/chart
and boundary laws, separate ripple/thermal/root controls and energy ledger.
Every accepted sample must satisfy the caller's maximum temperature excursion;
every sample's material refresh is admitted before a joint result is returned.
The endpoint response is admission for the next rebuild, not a second heat credit.
A failed source query, excursion, native solve, budget or cancellation returns no
candidate and leaves all original inputs unchanged. Native samples are moved out
without duplicating their trajectory vector.

This extension does not preserve transient polarization/field history through
coefficient changes, close reversible parameter work, solve switch-on Maxwell
dynamics, perform an implicit/converged thermo-optical iteration, or support
changing geometry/density/phase. It therefore does not complete MR07. Feedback
interval refinement and the validity of the quasi-steady approximation remain
caller-owned. For statistical inputs, the nonlinear thermal result is the
trajectory driven by expected heat, not an ensemble-mean temperature.

ThermalInterferenceLimits bounds total authored lines and conservatively sums
coexisting native numerical workspace before each large field preparation.
Existing limits separately bound field, frequency-pair, depth, reservoir, FFT and
spectral work. Material-source payloads, borrowed data and allocator overhead
remain outside these numerical workspace figures. Native thermal output/work caps
remain independent.

Three focused G1/G3/G4 integration regressions in
tests/interference_film_heating.rs compiled, linked and executed natively:
3 passed, 0 failed. The exact source query/refresh, conductor, Debye and phase
owners, original Maxwell film/two-port/interior/angular/exposure/reservoir/heat
modules, shared harmonic math and native lumped enthalpy engine all execute.
The scoped executable supplies a one-point identity FFT boundary (all fixtures
use a 1x1 periodic grid), unreachable optical-family refusal types, and a nominal
ConductionError envelope. It does not execute the original FFT implementation or
build the complete crates. Both new Rust files pass native rustfmt.

The independent oracle uses direct native two-port interior fields, original
Debye ODE currents and the exact eight mixed real/proper coefficient states.
Regressions check temperature-dependent permeability/conductivity/polarization
refresh, native thermal energy, subsequent-interval feedback, a 2^53+0.25 DD
source clock, and source-domain, temperature-excursion, mass/volume, resource and
cancellation refusal. These are authored constitutive model tests, not
experimental material qualification, a spatial thermal-uniformity proof, or
verification of an implicit/transient electromagnetic feedback scheme.

<!-- porous-terminated-feedback:start -->
### Passive-backed porous heating on the shared thermal trajectory

The existing harmonic and spectral porous heating owners provide opt-in
`try_new_terminated` and `try_new_oblique_terminated` constructors. Covariance
heating also provides `try_new_independent_modes_terminated`, including the
original flux-normalized diffuse hemisphere. The authored right-face pressure
or covariance component must be zero: a passive rear load determines the
counter-wave. A diffuse source therefore has zero right irradiance.

Each finite passive pressure/volume-flow load or exact rigid cap remains paired
with its INPUT frequency/angular mode before canonical sorting. Original
normal-mode references, tangential wavenumbers, and loads stay fixed while
accepted cell temperatures update the original gas and JCA(L) coefficients.
Each changed mode solves its terminal feedback again through the existing
conditioned stack. The original positive local work and physical cell binding
remain authoritative. There is still one gas evaluation per changed physical
cell, independent of frequency/angular count, and one stack solve per changed
mode. Unchanged temperatures reuse the complete cache.

The frozen `PorousStackHeating` adapter and all native feedback paths consume
the same prescribed-source account. Raw two-face boundary diagnostics remain
available, but reflected rear pressure is not counted as fresh source power.
`PorousSpectrumPower::load_absorbed_power_w` and
`PorousHeatingEnergy::load_absorbed_j` report work exported into a passive load.
The native acoustic balance is incident minus outgoing minus load work minus
deposited pore work. The load is outside the bound porous control volumes;
its loss never enters their source map or enthalpy. Existing open-face APIs
retain a zero load channel and their original fields, source powers and heat.

Native accepted-start timing remains unchanged: source, reflected, load and
positive cell watts integrate over the same accepted substep, then endpoint
temperature refresh prepares the next rates. Original thermal storage,
conduction, authored heat, events, continuation and transactional refusal stay
native. The backing is a fixed quasisteady acoustic load, not a new thermally
evolving reservoir, structural wall model or transient acoustic energy store.

Focused tests compare native source/load/cell powers with an independent
rho/K transfer solve and direct vector-velocity volume integration, verify
accepted-start load joules and refreshed endpoint rates, and refine the coupled
two-cell thermal trajectory against independent RK4. Frozen/harmonic/spectral/
covariance/angular one-mode paths agree, load-to-mode association survives
sorting, and rigid-backed diffuse heat converges with angular refinement while
preserving exact hemispherical input flux and shared gas work.
<!-- porous-terminated-feedback:end -->


## Prepared full tensor mechanics with exact shared relaxation history

`tensor_material::prepared::PreparedTensorBody` provides a fixed-clock cached
route for the original full three-dimensional tensor material, consistent inertia
and force/velocity coordinates. It consumes the existing `ProjectedTensorProny`
points, retaining every equilibrium tensor, Kelvin factor, Maxwell factor,
positive quadrature volume and exactly represented relaxation time. It does not
fit modal damping, lump mass, merge nearby poles, remove weak material directions,
or change the caller's displacement basis.

For each BIT-IDENTICAL positive source time tau, one complete displacement-memory
vector satisfies `y_tau' = M^-1 p - y_tau/tau`. Original point/branch springs are
`z_pb = A_pb y_tau`, and the mechanical force is
`p' = f - K_eq u - C M^-1 p - sum_tau K_tau y_tau`, where
`K_tau = sum_(pb at tau) A_pb^T A_pb`. Every source factor remains available for
original storage. Thus the dynamic state has `(2 + distinct_times)*coordinates`
entries instead of six extra states for every branch at every quadrature point.
This is an exact shared-coordinate realization for the admitted source history,
not a reduction of mechanical modes or a constitutive approximation.

The constructor binds the source and clock; borrowing histories explicitly choose
`InitialMemory::Relaxed` (all y=0) or `Instantaneous` (each y=u0). No prior work or
heat is invented for that initial state. Nonzero free eigenstrain refuses because
its initial offsets need not belong to the displacement-driven source image.
Arbitrary imported point memories, changing laws/geometry/supports, plane-stress
descriptors, finite rotation and thermal/source rebinding retain their separate
owners. Endpoint source-strain admission is enforced; a within-interval strain
or frequency-band enclosure is not claimed.

`PreparedQuadraticFlow` caches the SAME native coupled generator and propagation
kernel once. Each declared heat group uses `PreparedSourceLoss` on the ORIGINAL
positive Kelvin and Maxwell source squares. Diagonal similarity conditioning
changes only numerical coordinates. When source rows exceed the state dimension,
the existing native Householder QR supplies a full factor without rank truncation;
its normalized Gram matrix must match the original source Gram, including cross
terms. Group labels are explicit, contiguous and nonempty; grouping changes
observation cost and attribution only. One original element can collect all its
quadrature heat in one cache.

A warm interval performs cached contractions, original source-domain and energy
checks, and independently integrates held-force work `f dot (u1-u0)`. Equilibrium
and kinetic storage use original K_eq and M; Maxwell storage is the sum of the
squares of EVERY original A_pb*y_tau. Stable quadratic increments and compensated
work/heat histories check `W - D - Delta H` per interval and cumulatively before
publishing state, clock or ledgers. Heat is a positive source observation, never
an energy-deficit fill. Finite force jumps apply no impulse. Histories clone only
O(state) accepted data and scratch; immutable dense caches are borrowed.

Explicit capacities cover coordinates, points, complete states, heat groups,
ALL preparation/cache bytes and multiply terms, warm work and accepted steps.
Preliminary grouping is admitted before allocation. The existing native dense
limits remain: the base split, extended flow and doubled positive-loss state
must each fit their own 256-coordinate bound. No capacity escalation or hidden
truncation occurs. Preparation includes original Cholesky/eigensolve/QR and
bounded flow-integral calls; cancellation is checked around those existing
atomic native kernels and throughout source/cache work. No hard cold-path
cancellation-latency or real-time timing claim is made.

Six `tests/tensor_material_prepared.rs` tests pass natively: an independently
solved standard-linear-solid with exact roots -1,-2,-3 and analytical dashpot heat
for both prehistories; separate analytical heat for unequal elements with reversed
channel labels; coupled anisotropic, nonproportional Kelvin and consistent-mass
agreement with the original uncompressed owner (including every reconstructed
spring coordinate); exact versus adjacent-f64 source times; capacity/eigenstrain
admission; and bitwise rollback for malformed input, material-domain refusal,
cancellation, exhausted steps and work-budget refusal. The native source-closure
build uses the original material law/projection/history, full fs-phs and native
factor/eigen/GEMM/SIMD kernels. It is not a full workspace Cargo or DSR pass.

## Prepared tensor-memory motion on complete original solid coordinates

`tensor_material::solid::prepared::PreparedTensorSolidBody` binds the prepared
shared-history runtime to `ProjectedViscoSolid`. Its force projection, nodal
motion expansion, complete consistent mass, tensor strain maps and positive
quadrature all come from that original source. An identity basis retains every
represented nodal coordinate. The adapter does not create a second eigensolve,
constitutive law, damping fit or geometry approximation. Other supplied bases
retain the existing projection owner's explicit approximation obligations.

Heat groups are the ORIGINAL element indices, so all quadrature contributions
from one element share one positive source-loss cache. The public step exposes
that element heat and the original mean nodal velocity. Original held nodal work
is checked against retained-coordinate work using the same absolute energy
allowance. An O(state) candidate history lets force projection, source checks,
physical work/heat checks, nodal expansion and the final cancellation poll finish
before publication. No dense prepared matrix or immutable source is cloned per
interval. The source projection retains its own admission; prepared controls
bound the additional runtime/caches. Material changes, eigenstrain, moving support
maps and imported arbitrary point memories remain outside this fixed-source path.

Two native `tests/tensor_solid_prepared.rs` tests exercise actual source assembly.
The complete TET10 case uses all 30 original nodal coordinates and all 27 positive
material points with stiffness quadrature 3 and consistent-mass quadrature 4.
It retains signed consistent-mass entries. One Maxwell branch per point requires
222 states and a 506-coordinate block in the old point-memory path, which still
correctly refuses its existing bound. The exact shared history retains the SAME
source mechanics in 90 states. The test independently recomputes original sparse
K/M energy, every original Maxwell spring energy, nodal applied work and material
heat over 480 intervals; maximum observed physical balance error is 2.436e-17 J.
A separate two-element test swaps Kelvin and Maxwell materials and verifies heat
follows the original element, fixed-row loads do no work, and malformed outer
input preserves accepted state.

`examples/prepared_tensor_solid.rs` runs a complete free TET10 with anisotropic,
nonproportional Kelvin and Maxwell source tensors, original consistent inertia,
explicit 1e-5 affine initial strain and instantaneous memory. Its balanced nodal
pair now uses exact affine intervals: force magnitude rises from 0 to 0.2 N over
12 intervals, falls to zero over the next 12, and then releases freely through
interval 480 at 48 kHz. Each side's pulse impulse is 5e-5 N*s. All 90 physical
states remain; affine force records require 120 propagation coordinates.
CSV includes extension, equilibrium/kinetic/Maxwell storage, external work,
physical heat and balance errors.

The native affine example reaches 10 ms with extension 2.175098106e-6 m,
positive material heat 5.432572103e-5 J, maximum cumulative balance error
3.557e-17 J and maximum original nodal-work mapping defect 1.805e-20 J.
Source data are synthetic; this is complete assembly/projection/motion evidence,
not experimental calibration, spatial convergence or a real-time guarantee.
Focused native validation passes six existing prepared-material tests, three
new affine-material tests, five solid tests and the updated example. It does
not assert a full workspace Cargo, DSR or RCH result.

<!-- porous-load-thermal-receiver:start -->
### Opt-in native thermal destination for passive backing absorption

`PorousStackHeating`, `TemperaturePorousHeating`, `PorousSpectrumHeating` and
`PorousCrossSpectrumHeating` expose `with_load_receiver`. The latter includes
independent angular and diffuse modes. An explicit existing native cell receives
the original positive terminal work. It must be distinct from every porous
cell, active, and retain its original identity, geometry, mass, material,
storage basis and caloric identity. Its temperature must remain within its
native caloric domain. Binding is rechecked during source preparation and after
the final accepted events before publishing the returned candidate.

This optional binding declares that the complete real-impedance work thermalizes
in the identified lumped body. No material, heat fraction or destination is
inferred. Escaping radiation remains exported unless a separate receiving
absorber justifies that thermal destination. The receiver's native storage,
latent heat, boundary conditions, authored sources, events and conduction links
remain unchanged. Backing-to-pore conduction can change pore temperatures and
the subsequent acoustic response. The backing's temperature does not change
the prescribed impedance or independently trigger a new acoustic solve.

The default remains the export-only load behavior above. `load_absorbed_j` and
`load_absorbed_power_w` always retain the original total load work. New
`load_deposited_j`/`load_deposited_power_w` report the work credited to a bound
receiver; `load_exported_j`/`load_exported_power_w` report the remainder outside
the native thermal system. Exactly one of these channels receives the load
work. `deposited_j` and `cell_energy_j` retain their pore-only meanings.
`load_energy_j` separately identifies the receiving cell and its credited heat.
The native thermal audit includes pore plus bound-load heat exactly once;
per-bin accounts are diagnostics, not additional sources. Acoustic closure
still subtracts total load work once, regardless of the thermal destination.

One shared source-map path serves all owners; independent spectral/angular
load powers are summed before a single native advance. Accepted-start rates
integrate over actual native step durations, including event clipping. Fixed
or latent pore temperatures reuse their acoustic caches. Tests use an
independent density/bulk-modulus transfer solve, volume quadrature and three
caloric ODEs with conduction to check refinement and conservation. They also
check frozen/harmonic/spectral/covariance/angular parity, one-bin fast paths,
latent backing storage, continuation, export-only defaults, exact zero loss,
distinct binding, changed mass, caloric-domain refusal and cancellation.
<!-- porous-load-thermal-receiver:end -->

## Regional continuous modal acoustics and original mouth memory

`mm_line::continuous::ContinuousModalLine::new_regional` consumes one original
`GasState` per uniform axial element through fs-duct's full regional Galerkin
owner. Local density and compressibility determine inertia, interface impedance
and modal cutoff without averaging gas states or changing the retained radial
basis. Every region has the same static pressure; this is frozen stratification,
with no mean flow, mixing or thermal evolution.

`new_regional_loaded` attaches a complete `ModalAdmittanceLoad` to that same
regional operator before preparing the shared exact flow. The original mouth
radius, mode coordinates and equilibrium pressure must agree; the load's own gas
may have a different density or composition. Every mutual term, dynamic load
state and direct conductance remains in the interconnected system. Only inlet
wave ports remain external. `new_loaded` shares the same original interconnection
and retains its existing admission order and uniform arithmetic.

External wave references use each face's own original gas and area. Stored
energy, load heat and incident/outgoing wave energies retain their separate
physical accounts. The duct interior remains lossless; this constructor does
not add a radiation law, fit damping or establish the attached source's
frequency accuracy. Regional loaded regressions cover exact repeated-gas parity,
mutual-load memory and heat under three clock subdivisions, source/capacity
refusals and preservation of accepted history after malformed input.

## Exact affine forcing of prepared tensor mechanics

`PreparedTensorBody::from_mass_matrix_affine` prepares the original tensor
mechanics for F(s)=(1-s)F0+sF1, s=t/h. A positive caller-supplied force scale
conditions numerical coordinates only. It neither clips forces nor changes
material parameters. Each interval seeds its explicit force record; finite
force jumps between intervals apply no mechanical impulse. Accepted state
remains the original displacement, momentum and shared Maxwell histories.

The existing native extended-input flow propagates the complete coupled
trajectory. Its signed quadratic observation integrates F(t)^T M^-1 p(t)
independently of original physical energy and positive material heat. The
existing positive source-loss owner integrates every original Kelvin and
Maxwell channel on that same affine trajectory. No endpoint-force average or
energy residual substitutes for applied work or physical heat.

`advance_affine` returns the ordinary motion/heat report and the first normalized
velocity moment integral_0^1 s*v(h*s) ds. One separately prepared native
one-way displacement-integral exponential supplies all coordinates' moments.
Its observers do not feed back into motion or material history. Endpoint-weighted
moment work must agree with independently integrated force/velocity work before
publication, in addition to the original interval and cumulative energy gates.
Preparation is cached; warm steps use bounded contractions and source checks.

`maximum_state_coordinates` bounds physical plus force-record propagation
coordinates; `state_dimension()` reports physical coordinates and
`propagation_state_dimension()` additionally includes the force records.
All flow, source-loss and observer native dimensions, storage, preparation work
and warm work are admitted explicitly. Existing held preparation retains its
previous cache and evaluation route. An affine-prepared history also accepts a
held force through equal endpoints. Held-only preparation refuses affine input.

The fixed-source, fixed-clock, compatible-prehistory and zero-eigenstrain scope
of prepared mechanics is unchanged. Source validity remains an endpoint check,
not a certificate over all intermediate strains. This prescribed-force path
does not supply an implicit contact solve, evolving material law or nonlinear
finite-deformation integrator.

Native `tests/tensor_material_prepared_affine.rs` passes three focused tests:
an independently derived standard-linear-solid solution under a global ramp,
with work from integral F*v and heat from the branch-force square; constant
endpoint agreement with the original held runtime; and capacity/failed-input/
cancellation/exhausted-step/wrong-preparation atomicity. The analytical test uses
both compatible prehistories, three conditioning scales and 40 consecutive
intervals. It explicitly falsifies average-force times displacement as work
(maximum observed difference 2.485e-8 J). All six existing prepared-material
tests remain green. This is original-source native numerical verification,
not experimental material qualification or a full workspace Cargo/DSR/RCH pass.

## Affine nodal loading of complete original solids

`PreparedTensorSolidBody::new_affine` prepares the same original
`ProjectedViscoSolid` for endpoint-interpolated nodal force fields.
`PreparedTensorSolidHistory::advance_affine` projects both endpoints through
the original basis, advances the prepared tensor trajectory, and expands both
mean velocity and the first normalized temporal velocity moment through the
same nodal map. Fixed rows have zero velocity and work. The adapter checks
original nodal moment work against the independently integrated material-runtime
force/velocity work before publishing its O(state) candidate history.

Three native tests in `tests/tensor_solid_prepared/affine.rs` exercise the real
source assembly and full prepared consumer. A complete TET10 retains all 30
mechanical coordinates, 27 material points and 90 physical states; its affine
propagation has 120 coordinates and its separate native observer has 180.
The original consistent-mass body force produces the analytical rigid
translation under linearly changing acceleration, with independently checked
displacement, velocity, impulse, work and both velocity moments. Its measured
nodal-work discrepancy is 4.236e-19 J in magnitude and physical material loss
stays within the declared zero-strain tolerance.

A deforming TET10 triangular pulse has the same physical endpoint, external
work and original Kelvin/Maxwell heat when its intervals are split into 1, 2
and 4 parts. Original sparse K/M and original spring-factor storage are
recomputed independently; maximum measured physical balance error across these
runs is 4.089e-17 J. All three heat totals agree within 3.721e-18 J. This is an
exact-affine subdivision identity up to floating-point error, not evidence of
spatial convergence or arbitrary-load resolution. A supported TET4 test checks
fixed-row loads do no work and malformed endpoint fields preserve accepted
state. The original two held-solid tests remain green.

## Loaded hereditary shell resonance on original physical ports

`shear_shell::coating::harmonic::loaded` connects the original equilibrium-basis
film and substrate spectrum to `fs_modal::harmonic::PreparedHarmonicResponse`.
Both `MaterialShellImage` and `MaterialShellAssemblyImage` can prepare a reusable
`prepare_coating_harmonic_with_impedance` factor at one admitted frequency.
The original film/substrate laws are evaluated once at that drive frequency;
their complete coupled storage/loss and the unchanged image inertia enter
`K - omega^2 M - i omega C - i omega B Z B^T` before factorization. No inversion
of an unloaded resonance, per-mode loss approximation or extra material fit is
introduced. The existing unloaded solve retains its original operation order.

Port selection uses original force/velocity footprints or the original signed
exterior/cavity pressure-volume rows. The entire complex mutual impedance is
retained. Its Hermitian part receives the shared pointwise PSD admission; the
caller remains responsible for its units, source band, physical boundary law
and broadband causality. Internal background, substrate, interface and film
heat remain separate from the load's net absorbed boundary work, with an
independent combined mean-power residual. Absorbed boundary work alone is not
a radiation or material-heat certificate.

Each returned motion passes the existing original whole-film spatial/cycle
strain-envelope observer. An insufficient envelope refuses the response and
consumes that completed strain query; it does not establish material failure.
Preparation and solve poll both the original source gate and shared execution
context. Source, factor, port and strain work/storage limits remain explicit.
This is a bounded fixed-basis periodic response, without a continuum/modal
truncation certificate or startup, moving-material or nonlinear response claim.

`tests/coating_harmonic.rs` covers zero-load parity and factor reuse, signed
fluid ports, an independent full two-port feedback oracle, material/load power
separation, single-component assembly parity, whole-film strain refusal,
pointwise passivity, source/resource refusals and cancellation (G0/G1/G3/G4).

## Prepared tensor material continuation at accepted instants

`tensor_material::prepared::PreparedTensorHistory::rebind_materials`
continues an accepted complete-coordinate state on a replacement prepared
operator. `PreparedTensorMaterialChangeBudget` gives an explicit maximum
absolute externally supplied coefficient work. The candidate reports signed
equilibrium and Maxwell storage changes and a replacement history; the original
accepted history remains unchanged on success and refusal.

Reference mass, original point order and volumes, strain columns, heat ownership,
ordered branch counts and the positive step duration remain fixed. Both endpoints
must have zero eigenstrain. Original branch identity and material applicability
belong to the source-binding caller; the numerical API cannot authorize replacing
an unrelated internal variable merely because its dimensions match. The new
material must admit the accepted displacement and remaining step capacity.

The state retains complete displacement-history vectors, so changed tensor
factor rank requires no inverse, pseudoinverse, spectral cutoff or strain reset.
One old relaxation-time group may split by copying its exact history. Several
old groups may acquire the same new time only when all their complete accepted
histories compare exactly equal. An incompatible collision returns the typed
`MemoryCollision` refusal rather than averaging history or fitting damping.
Held/affine preparation and its numerical force scale may change because those
caches carry no accepted physical history.

Displacement, momentum, velocity, time, accepted steps, initial energy and
compensated force-work/heat accounts survive. A distinct compensated signed
parameter-work account records the fixed-state storage change. Stiffness and
original factor-row differences evaluate that work without subtracting unchanged
kinetic storage. Pure relaxation-time or Kelvin-viscosity changes have zero
instantaneous coefficient work and create no instantaneous heat. Both held and
affine cumulative balances include parameter work. A finite source allowance,
replacement-domain admission, original physical energy balance and cancellation
are checked before returning a complete candidate.

Focused native evidence: three new continuation tests plus the six held and
three affine prepared regressions pass. The continuation cases compare real
post-change motion, spring history and independently integrated physical heat
with the analytic three-exponential Zener solution; exercise positive/negative
coefficient work and changed factor rank; split and merge equal histories; refuse
unequal-history collisions; and preserve accepted state under mass/map/clock,
work-budget and cancellation refusals. This is prescribed endpoint coefficient
evolution, not a solved thermal expansion, entropy or heat-feedback model.


## Complete causal material resonance of shell images

`shear_shell::coating::poles` prepares one immutable continuous model for
`MaterialShellImage` or `MaterialShellAssemblyImage` through
`prepare_coating_resonance`. It retains the original mass-normalized equilibrium
basis, complete background/glue viscosity, original interface dashpots, and
every already-selected substrate and coating memory row. The existing
`fs_phs::PortHamiltonian::with_relaxation_branches` owns the internal-strain
storage and loss; `fs_phs::harmonic::HarmonicModel` owns the prepared generator,
local pole refinement and contour discovery. No new eigensolver, material fit,
branch aggregation, projection or modal-damping approximation is introduced.

State order is interleaved equilibrium modal displacement/momentum, followed
by all substrate memory coordinates and then all coating coordinates in their
original selected order. A mechanical seed can be lifted with `p=s*q` and
`z=c.q/(1+s*tau)`; complete caller-supplied seeds also admit material relaxation
directions, including at a source branch rate. Singular seed lifting refuses
instead of removing those states. The full polynomial/linear state model is
used during iterations and contours, without eliminating Maxwell denominators.

Frozen films and unselected anchored constant-loss substrate images refuse:
the existing explicit causal selections are required. Coating Prony lowering,
source-fit tolerances, constant-loss substrate continuation and the original
fixed spatial/bubble/mode basis remain independent approximation boundaries.
The prepared image retains the intersection of real-frequency source bands;
each mode reports whether its oscillation frequency lies there. The explicitly
selected global continuation also has real relaxation modes and can have poles
outside that real-axis band. Neither flag establishes experimental accuracy at
complex frequency or widens source authority.

Each pole carries original force-port motion and signed exterior/cavity volume
flow, the complete state, original-equation residual, and independently evaluated
mechanical storage, background/interface heat and per-memory-branch storage and
heat. Equal-time coating factors may mix source films, so no per-film allocation
is invented. `H=x*Qx/2`, `D`, and `2 Re(s) H+D` are complex-envelope accounts,
not cycle averages of a decaying real trajectory. Independent energy/loss
partition and decay residuals must pass the declared physical tolerance. Shape
normalization is not finite-amplitude strain admission or acoustic observability.

Complete coordinates, conservative preparation payload/work and separate shared
query budgets remain explicit. Physical observations reserve their resources
before spectral work. Preparation polls source rows and before/after synchronous
shared physical admissions; the shared spectral iterations and observations poll
the source cancellation gate. Source/history is immutable on success or refusal.
Count, rank, contour conditioning and nested quadrature diagnostics are numerical
evidence, not rigorous spectral completeness or damping enclosures.

`tests/coating_poles.rs` exercises independent scalar rational dispersion and
physical decay, a real material relaxation root, contour quadrature refinement,
full cross-mode memory and single-component assembly parity, source/resource
refusals, cancellation and exact replay (G0/G1/G3/G4).


### Continuous harmonic scattering of the original multimodal line

`ContinuousModalLine::prepare_harmonic_scattering` exposes the same complete
physical duct/load system through the shared `fs_phs::harmonic::scattering`
owner. Every external inlet mode and every unbound mouth mode uses exactly
the original runtime's pressure-wave reference and port order. The reference
feedback is closed in the original continuous operator before factorization;
the accepted runtime state and time step do not enter that calculation.

The returned model retains all original duct and attached load states and
physical resistances. Responses report original complex state, pressure,
physical volume-flow inputs and independently calculated material/load heat,
together with incident/reflected wave amplitudes and powers. Wave-reference
exchange remains a boundary account, not added material dissipation. One
frequency factorization can serve several incident-wave vectors. Preparation,
factorization and query use explicit shared harmonic/scattering budgets and
cooperative cancellation without changing accepted runtime history.

A static mouth conductance retains its original ordinary-input source slots;
set those slots to zero for the unforced termination. A dynamic mouth load
retains all of its original state and consumes the interconnected mouth slots.
Regional gas density/compressibility and mutual load coupling are unchanged.
The scalar positive reference is a wave-coordinate convention, not a matched
evanescent characteristic. This adapter supplies no new radiation law, physical
material validation, spatial accuracy band, steady-state convergence guarantee
or continuum mode-completeness claim.

`tests/continuous_modal_scattering.rs` exercises a two-mode cylinder's analytic
reflection/transmission on three axial rungs with port, rigid and pressure-
release ends; reverse-port reciprocity and shared-factor reuse; a regional
two-mode dynamic mouth with independent original oscillator/direct-conductance
heat; bitwise harmonic-query independence from accepted history and three
runtime clocks; and a static load's original ordinary source slots with no
artificial memory states (G1/G3). The synthetic axial ladder is numerical
evidence for that declared source, not measured-material validation.


## Native rigid cohesive loss and finite-body heat

`friction_heat::rigid_cohesive::apply_rigid_cohesive_heat` advances the original
`RigidCohesiveSystem` and original finite enthalpy network through one common
native timestep. Every original site has exactly one matching-area thermal
route; explicit ordered body/material mappings retain A/B identity. Several
sites may share finite receiving cells. This is a fixed material-reference
patch/conduction model: native finite rotations keep their objective mechanical
work, but the thermal geometry is not remapped into current world coordinates.
Moving footprints, changing area/coverage and finite adhesive-layer thermal
resistance/storage are outside this owner.

The actual audited native per-site damaged-area friction and surviving-bond
viscosity supply heat exactly once. Fracture remains a distinct nonthermal
surface-creation account. The original common-temperature, massless-face network
solves signed heat into both bodies through their own k*A/d conductances;
ordinary hot-to-cold exchange can make one body's net heat negative. No prescribed
partition, thermal reservoir or numerical-loss heating is added. Both ends of
the original straight local material path must remain closed. Mechanical source
coefficients stay frozen during this increment; source temperature continuation
and temporal refinement remain explicit.

The native stage/publish boundary admits the complete thermal candidate before
any body, material history, observation or clock changes. A reusable frame buffer
keeps ordinary native stepping allocation-free after construction. Source,
nonlinear, domain, energy and cancellation refusals preserve both accepted
owners; final accepted native and caloric publication is infallible. Independent
whole-system accounting checks the change of kinetic/gravity/bond plus thermal
storage, adding fracture and subtracting applied mechanical work and other heat.
The shared finite thermal solver is also used by the deformable cohesive owner.

`tests/rigid_cohesive_friction/heat.rs` exercises independent vector-friction and
bond-viscosity work, analytic unequal finite-cell heating with signed ordinary
conduction, continued exact native-mechanics parity, progressive fracture kept
out of heat, and late thermal-source/cancellation rollback (G1/G3/G4).

`friction_heat::third_body::cohesive::propose_deformable_cohesive_layer_heat`
adds real finite adhesive storage to one original sourced deformable cohesive
patch. The source card's named third body maps explicitly to an existing layer
material; ordered surface/body/material identities, both contact areas and the
full layer cross-section must match. The native mechanical candidate supplies
only surviving-bond viscosity to the existing `ThirdBodyHeatStep`; its original
mass/caloric law retains heat, while two source contact resistances and geometric
bulk resistances control delayed asymmetric body heating. Fracture remains
nonthermal. Nonzero damaged-surface friction refuses because its localized
surface destination is absent. Mechanical coefficients remain frozen and no
source parameter work is thermalized. The small-deformation thermal geometry is
a fixed material-reference layer with a closed normal path; the original source
must be calibrated for that thickness. Open cracks and evolving thickness are
outside this split. Optional entropy covers caloric
storage, conduction and viscosity, without claiming a fracture/surface-entropy
closure. Both immutable candidates remain unpublished on source, nonlinear,
energy or cancellation failure. The `friction_heat/cohesive/layer.rs` regression
uses independent three-capacity conduction and bilinear damage/viscosity
integrals, plus native friction, identity and late-refusal cases.

`friction_heat::third_body::cohesive::propose_deformable_cohesive_layer_heat_with_temperature`
adds prescribed thermal weakening to that same finite adhesive transaction.
Caller-owned successors retain the original selected claims, nonthermal query
coordinates, patch geometry, elastic operator and accepted irreversible history.
The existing native transfer runs once at held geometry before motion; thermal
iterations see actual subsequent damage without repeating the source change.
Signed coefficient work remains explicit external work, and both fixed-jump
thermal fracture and later mechanical fracture stay outside caloric heat.
Only the actual subsequent surviving-bond viscosity heats the finite adhesive
and conducts through its original two resistances. An independent audit reads
original/final bond storage, accrued fracture and finite thermal storage, with
the unchanged initial native kinetic/bulk energy and signed source work; the
three original owners' residuals provide a separate check. Cooling cannot heal
damage. Either both final candidates are returned or both accepted owners stay
unchanged, including on foreign-source and late cancellation refusals. This is
explicit prescribed temperature with frozen coefficients during motion, without
cohesive entropy or closed thermoelastic feedback. The layer temperature tests
compare analytic inertial unloading and three-capacity diffusion after thermal
weakening, continued cooling without healing, and late source/cancellation
rollback (G1/G3/G4).

`propose_deformable_cohesive_layer_heat_with_deformation` and its
`with_temperature_and_deformation` counterpart admit tensile and compressive
extension of that same surviving adhesive layer. `CohesiveLayerDeformation`
declares the absolute normal and Euclidean shear strain bounds relative to
original layer thickness. Both endpoints are checked; the native straight
constitutive path and irreversible damage preserve those bounds within the
step. Complete rupture refuses before thermal trials, requiring a separate
crack/contact transition. Partial damage retains its fracture account without
being interpreted as missing thermal area.

The fixed-area, laterally uniform layer has endpoint thickness `L0 + g_n`.
`ThirdBodyHeatStep::with_layer_normal_stretch` uses the physical half-thickness
in each original series resistance, retaining physical spatial conductivity,
both contact resistances, body half-distances and original mass/caloric storage.
The report exposes `layer_normal_stretch`. Ambient side-area convection refuses
on this path; explicitly owned total sources and reservoirs retain their units.
Edge overlap, lateral-area changes, open cracks and closed thermoelastic
feedback remain outside the model. Existing closed-reference entries preserve
their original behavior. Focused `friction_heat/cohesive/layer/stretch.rs`
regressions compare native inertial motion and viscous heat with independent
physical-length three-cell balances, zero-extension parity, progressive damage,
source-temperature weakening and transactional strain/rupture refusals.

### Prescribed source temperature with native cohesive heating

`apply_rigid_cohesive_heat_with_temperature` stages one successor of each
original cohesive source, plus its paired damaged-area friction source when
enabled, before the same native mechanical increment and finite-body heat solve.
Selected claims, ordered surfaces, held query coordinates and source history
remain pinned by the original transfer owners. The transfer runs once at
accepted geometry; nonlinear heat trials never repeat thermal damage or advance
the mechanical clock. The complete temperature/mechanical/caloric transaction
publishes only after the original audits and final cancellation gate succeed.

The prescribed coefficient change contributes signed `parameter_work_j` as an
external input to the combined energy balance. It is neither deposited as heat
nor silently taken from caloric storage. Fixed-jump thermal damage remains in
the fracture account; only actual bond viscosity and residual sliding heat the
finite receivers. Mechanical coefficients remain frozen during each increment.
This explicit temperature control does not provide cohesive free-energy
temperature derivatives, surface entropy, or closed thermomechanical feedback;
the optional entropy audit covers the original caloric transport and physical
heat destinations. The fixed-reference, closed-path geometric restrictions of
the ordinary heat operation remain in force.

The native heat regression target verifies exact parity with the original
paired-source mechanical step, independent hot friction/viscosity and finite
cell heating, heating/cooling with fixed-jump damage excluded from heat, and
late conductivity/cancellation rollback of both source owners and both physical
states (G1/G3/G4).

### Native rigid adhesive stopping and interval continuation

`RigidCohesiveSystem::advance_interval` continues the original finite-rotation
pair under constant loads, using its existing nominal timestep as the maximum
substep. Only a typed unresolved local Coulomb reversal triggers refinement
from the same accepted predecessor. Actual objective endpoint velocity can
resolve stopping and restore the original static disk; no projected velocity,
new impulse, compliance, or numerical heat changes the physics. This brackets
an admissible boundary of the discrete path model, without certifying the first
continuous event or a unique softening branch. Source coefficients remain at
their currently accepted temperature; this operation performs no heat solve or
additional temperature transfer.

Reusable native solver matrices are retained. Three bounded accepted-state
buffers preserve both original body motions, every site history, source
reference, friction branch/effort, observations and clock. Any outer failure
restores the incoming state, including after successful internal substeps.
Separate cumulative trial, substep, evaluation, material and iteration ceilings
bound rejected work as well as accepted work. Physical acceptance tolerances
are unchanged. Whole-interval endpoint energy and momentum are checked against
the sums of actual accepted work, impulses and per-site fracture, viscosity and
friction losses. Final native observations still describe the last substep.

Ordinary stepping retains its original multiplication clock until an explicit
interval is used; subsequent nominal steps start from that interval's endpoint.
`next_step_time_s` exposes this shared endpoint to coupled owners. The focused
native regressions cover an unaligned analytic translation-limit stop at cold
and previously transferred hot sources, resticking and renewed sliding, exact
rotating multisite native-path parity, and late outer rollback with retry.

`friction_heat::rigid_cohesive::interval::apply_rigid_cohesive_interval_heat`
composes that native interval with the existing finite-body heat network.
Accepted per-site friction and surviving-bond viscosity sums supply heat once,
including sliding before a final static substep; rejected stopping trials and
fracture never become caloric heat. Original body/material/area assignments are
shared with the fixed-step adapter. The largest normal jump over the initial
and every accepted local endpoint must stay closed, so intermediate opening
cannot be hidden by final reclosure. The native interval's rollback checkpoint
encloses the thermal/source solve and final cancellation gate; both complete
owners publish together. Whole-interval mechanical and caloric energy close
independently. Thermal transport uses interval-average physical power and one
backward-Euler endpoint: shorter complete intervals refine heating timing, and
this frozen-coefficient split does not imply closed temperature feedback or
surface-entropy closure. The focused heat/interval regression covers analytical
stopping heat with signed body conduction, static hold and breakaway, late
thermal/cancellation rollback, and native opening followed by reclosure.

## Independently dynamic loads in causal shell resonance

`MaterialShellImage::prepare_coating_resonance_with_load` and the corresponding
monolithic `MaterialShellAssemblyImage` method compose the unchanged causal shell
with a complete `ShellResonanceLoad`. The caller supplies an existing exactly
quadratic `PortHamiltonian` and one original `ShellHarmonicPort` per load input.
Force ports retain their physical velocity/force units; exterior and cavity
ports retain their independently signed volume-flow/opposing-pressure units.
Repeated, signed and nonorthogonal original projections are permitted. Mutual
load coupling is retained through the original full J/R/Q/G matrices.

The existing `with_parallel_load_columns` owner connects load input to shell
flow and applies minus the load output through the same shell columns. All load
states follow the original shell, substrate and coating states. The shared
fs-phs pole and contour solvers see this complete joint generator. No isolated
impedance inverse, load condensation, new mode projection, damping ratio or
material fit is introduced. In particular, an independently evolving load mode
survives when signed shell projections cancel and its port is unobservable.

`ShellResonanceLoadResponse` reports every original load coordinate, physical
flow and effort, signed absorbed power per port, load storage and actual R loss.
Those quantities use complex-envelope conventions, not finite-cycle averages.
Independent load rows check `s*z = (J-R)*Q*z + G*flow`. Independent accounts check
`2 Re(s) H_load + D_load - P_load`,
`2 Re(s) H_shell + D_shell + P_load`, and the complete joint decay balance.
Source mechanical/material channels remain distinct from load dissipation;
absorbed port work is never relabeled as radiation or as irreversible heat.
All these accounts must pass the declared physical tolerance before publication.

The complete dimension, including every load coordinate, is admitted before
allocation. Preparation payload/work include original load matrices, full joint
interconnection and shared physical admissions. Observation resources reserve
both load-state and load-port work before spectral queries. Cancellation is
polled around synchronous shared admissions and through source/load rows.
The mechanical seed helper initializes appended load coordinates to zero;
callers can supply a complete load-dominated seed or use contour discovery.
Source images and accepted histories are immutable on success or refusal.

Nonquadratic or undeclared load storage refuses. This API accepts an ordinary
continuous pHS owner, not a descriptor load or a separately authored algebraic
feedthrough law. The load's finite authored real-frequency band is intersected
with the shell source bands and retained with results. Selecting the owner
selects its global causal continuation; out-of-band/real-relaxation results
retain the existing explicit support flag and no experimental accuracy claim.
The original source Prony/constant-loss approximations, fixed shell basis and
shared numerical count/rank limitations remain unchanged.

`tests/coating_dynamic_poles.rs` exercises independent hereditary dispersion of
a moving-support mass/spring/damper attachment, the two split resonances and
their joint contour discovery, original mutual force-port coupling, a pure load
pole at an isolated component resonance with zero net signed volume coupling,
and dimension/band/port/cancellation refusals (G0/G1/G3).

## Continuous original-state rational modal resonance

`mm_line::continuous::rational::analysis::PreparedRationalModalAnalysis`
provides continuous forced response, isolated finite pole refinement and
contour discovery for the original complete modal Galerkin duct and an
already admitted `ModalRationalLoad`. Its direct constructor supports zero,
invertible and mixed-rank direct impedance matrices; the existing time-domain
runtime's mixed-rank limitation does not apply to this analysis. The existing
`RationalModalLine::prepare_analysis` adapter reads its original physical
sources without reading, changing or resetting accepted trajectory state.

The retained descriptor coordinates are every duct displacement and canonical
pressure momentum, every original real rational pole coordinate, and every
outward mouth-flow multiplier. The last group has no storage. The original
equations are `x' = A_d x + G_in (2a-P_in)/Z0 - G_out U`,
`z' = A_l z + B_l U`, and `P_out = C_l z + D_l U`. No direct matrix, bare
impedance or descriptor mass is inverted, and no epsilon resistance is added.
Mutual mode coupling and zero-transfer-residue internal states remain present.

All harmonic queries use peak `exp(-i omega t)` phasors, with `s=-i omega`,
and obey the original real-frequency load band. One shared `fs-la` complex LU
serves all incident modal drives at a frequency. Poles and spectra use the
same entire descriptor pencil and the shared `fs-la` Newton/Beyn owners.
Their complex continuation is the selected rational model, not an independent
exterior measurement or an exact fractional/material law. A pole separately
reports whether `|Im(s)|` lies in the original source band. Infinite roots and
unresolved repeated/defective modes are not reported as resolved finite poles.
The contour's argument count, rank, nested quadrature and conditioning reports
are numerical evidence, not interval bounds or a completeness certificate.

Constant reversible row/column conditioning changes no physical state.
Acceptance independently evaluates the unclosed original duct, memory and
pressure-continuity equations. Its tolerance is a rowwise backward error in
the fixed numerical coordinate chart: the denominator includes the original
right-hand side and the **split** absolute coefficient envelope times the
scaled state infinity norm. Opposing body/interface coefficients and `sE`
are not cancelled in that envelope. This definition remains meaningful when
contour discovery gives roundoff-sized components in an exactly dark physical
block; no state is clipped or projected to zero. Nonfinite accumulations,
resource violations, failed original residuals and cancellation refuse.

Observations independently retain all inlet/mouth pressure and volume-flow
phasors, duct storage, incident/reflected wave power, duct-side mouth work and
original rational-load work. Harmonic quadratic quantities are cycle means;
free-pole quantities are complex-envelope forms. The independent duct balance
is `2 Re(s) H_duct - inlet_work + mouth_work`; the interface defect compares
duct pressure work against `C_l z + D_l U` work. Generic rational coordinates
have no supplied physical storage matrix: their signed work is never called
heat and duct energy is never described as total coupled energy.

The prepared image admits the complete dimension (at most 1024), numerical
payload and preparation/factor/observation work before allocation. Pole and
contour numerical work has its own shared-owner limits; their explicit
iteration/node controls also bound callback assembly counts. Source identity,
fit/passivity reports, spatial/radial resolution and runtime clock remain
unchanged. Neither the runtime clock nor its held-wave states enter the pencil.

G0/G1/G3/G4 evidence in `tests/rational_modal_analysis.rs` covers zero and
mixed-rank two-mode mouths against independent continuous cylindrical transfer
equations, including evanescent propagation and mutual impedance; second-order
refinement at 4/8/16 axial elements; the loaded decay pole against the analytic
`Z(s)+Z0=0` pole over 1/2/4 elements; a finite pole exactly at an unobservable
original material pole; regular joint response at an exact bare duct mode;
bitwise clock/history independence; the existing runtime adapter; factor reuse;
and source-band, resource and cancellation refusals. This evidence does not
qualify nonlinear flow, evolving gas, continuum BEM accuracy, arbitrary fit
passivity between samples, or the full repository build.

### Forced response of complete causal shell and dynamic-load resonance images

`PreparedShellResonance::harmonic_factor` evaluates the same complete selected shell/substrate/coating/load state image as its existing pole APIs. `ShellForcedFactor::response` accepts peak complex forces at the original force ports, followed by signed exterior and cavity pressures, using `exp(-i omega t)`. Internal attachment ports remain internally connected. One shared physical LU serves all external drive vectors at each admitted positive angular frequency; the common authored shell/load frequency band is enforced before factorization. Bare component resonances need no separate impedance inverse.

`ShellForcedResponse` retains the complete joint harmonic state, original port outputs and per-port supplied work, plus physical modal and force-port displacement, velocity and signed volume flow. It independently accounts for cycle-mean mechanical storage, background/interface heat, every actual selected substrate/coating memory branch, and original load storage and resistance heat. Mutual load-port work remains signed per port. Shell and load balances, original forced mechanical/kinematic/material/load equations, and full-owner energy/heat are cross-checked. These forced cycle means are distinct from the existing free-pole complex-envelope quantities.

`ShellForcedBudget` admits the additional observation payload and work before the solve. Each branch/interface projection is reused across original mechanical equation checks. Cancellation, malformed input, insufficient resources or failed original physical accounting returns no partial response and preserves the factor/source. The fixed source basis and already selected causal material approximation remain the source model's limits. This is a linear image response; finite-amplitude coating strain admission remains the source image's separate geometric strain check.

`coating_dynamic_response` checks independently solved causal body/load amplitudes, actual source memory and load heat, a finite joint response at an exactly singular standalone attachment resonance, two-coordinate cross-memory and dynamic-load coupling, mutual signed port work, both oriented pressure drives, shared-factor reuse, and source-band/resource/cancellation admission. Existing shell pole and free-mode semantics remain unchanged.

### Mixed-rank continuous rational terminals

The continuous DG and modal rational consumers admit reciprocal positive-semidefinite mixed-rank direct impedance through one shared constrained-terminal owner. Existing invertible-direct and exactly-zero-direct paths retain their established solve order. The mixed path keeps every original duct coordinate, original rational pole coordinate, signed physical port, accepted source history, source clock, source band and source admission report.

For the original equations `x' = A x + F q + b u` and `H x = D q`, a fixed complete orthogonal range/null chart changes only the coordinates of zero-storage flow multipliers. The actual original direct range is solved through the existing terminal solver. The remaining original null-pressure constraints are differentiated once through the existing zero-direct owner. Returned flow rows and inlet factors are reconstructed in the original port basis, then checked against original pressure closure and the original state equation. No material or duct coordinate is removed, no history is projected or reset, and no resistance is added to regularize a constraint. Modal held-input coordinates let all incident modes participate in the same elimination.

A spectral cutoff alone is not authority to discard a fast physical direction. The shared Jacobi chart must agree with independent original-coordinate principal-minor checks through the shared LU owner: a resolved original principal range block is extended by every one- and two-coordinate combination, and the first exactly-zero LU pivot must agree with the proposed rank. Original null rows, componentwise reconstruction, orthogonality and closure are checked separately. Small nonzero positive directions, hidden negative pairs, unresolved/rounded-rank ambiguity, nonreciprocal direct matrices and singular higher-index constraints refuse. This is conservative numerical admission, not an exact-arithmetic rank or passivity certificate; a rounded source that acquires additional nonzero pivots may refuse even when its author intended an exact null.

The existing prepared exact generator and integral-observation owner advances the complete admitted finite state and computes interval means and work. Original load pressure remains `C z + D q`, including any computed inlet-flow contribution. Independent duct-side and load-side signed work must agree. Generic rational coordinates have no declared physical storage, so generic terminal work remains signed exchange rather than being labeled total stored energy or heat. Attachment and every update check the retained original pressure manifold; inconsistent history, nonfinite inputs, cancellation and failed updates do not mutate accepted history or accounts.

The shared mixed cold path requires a `256*n²` numerical byte allowance, with complete dimension `n <= 256`, bounded shared spectral/principal-minor work and cancellation polls. Each consumer separately admits its retained and exact-propagation storage; the modal consumer adds this workspace to its existing preparation envelope. The range/null chart and flow maps are prepared once; playback reuses the existing exact flow/observation operators. Spatial/radial truncation, source-fit error, source-band coverage and held-input clock error remain distinct from exact propagation of the admitted finite operator.

Focused evidence covers analytic common/differential compliance flow, signed port covariance, several positive and null channels, isolated and rotated tiny nonzero directions, hidden indefinite pairs, higher-index refusal, unchanged zero-direct regressions, mutual two-mode real/conjugate-pole history and work under exact substepping, clock convergence to independent weak harmonic and original joint-descriptor responses, preserved compatible history, inconsistent-history refusal and atomic rejected updates. These checks do not imply a full-workspace test run or certify arbitrary source fits.

### Mixed-rank rational terminations in the original DG trajectory

The continuous DG terminal path shares its reciprocal mixed-rank direct-term
closure with the modal duct path. It retains every original DG, compliant-wall
and rational pole coordinate. Only terminal flow multipliers are eliminated:
resolved direct-range flows satisfy the original impedance law, while the
remaining null-pressure constraints are differentiated once. Existing full-rank
and exactly-zero direct paths preserve their original arithmetic. Unresolved
rank, nonreciprocal sources, higher-index constraints and capacity or
cancellation failures refuse before changing accepted history.

The exact mechanical/DG renderer carries the resulting inlet-dependent flows
through its original mechanical velocity columns. Independent load-pressure
observations use the actual original `C z + D q`, including the computed
`D * inlet_flow_factor`; that term is mathematically zero on the admitted null
space, but is not erased from floating-point observations. Its mechanical dual
norm contribution participates in the existing interval pressure check. Each
side retains its signed interface work; rational coordinates are not assigned
invented storage or called heat.

Cold assembly admits at least `256 * n * n` numerical bytes for a mixed-rank
terminal, with complete original `n <= 256` and the shared spectral work cap.
The existing exact propagator separately admits its complete generator and
work forms. Caller-authored terminal ordering is checked against every original
pressure and flow column. Inconsistent initial or candidate histories remain
subject to the existing constraint checks; they are never projected or reset.

The four focused `dg_port::continuous::terminal_tests` exercise the actual cold
constructor with original DG and wall operators, a signed mutual two-port
rational relaxation law, and the actual prepared quadratic flow. They compare
separate integrated duct/wall and load work/heat accounts across zero, mixed
and full-rank direct terms; original-state semigroup agreement across clocks;
physical terminal permutation; and resource, cancellation and history refusal.
The fixture's load storage and heat are known independently from its authored
relaxation law, not inferred for arbitrary rational fits. Native checks compile
the exact production cold assembly with those real owners. They do not claim a
complete renderer, workspace or DSR validation run.


## Spatial moisture diffusion through compatible axial swelling

`MoistureLayeredSolid::advance_spatial_diffusion` explicitly interprets the selected scalar source as current-frame isotropic dry-basis Fick diffusivity: `j = -rho_d,current D grad_x M`. This is not the concentration-gradient law `-D grad_x(rho_d M)`. Conserving the dry skeleton gives reference flux `-rho_d,reference F^-1 D F^-T grad_X M`; the admitted separable prism reduces this to a cellwise `1/lambda_x^2` metric. Existing area-specific activity-transfer laws consume the actual exposed/contact area. Dry inventory, original geometry and source selections remain unchanged.

Every nonlinear moisture endpoint invokes the original compatible native solid solve. The original layered transport owner then evaluates its true residual and measured water balance with that endpoint's geometry, inside its existing bounded implicit iteration. Source, mechanics, work-limit and cancellation refusals publish no changed specimen. The result retains the original water accounts, measured axial stretches/current area, geometric residual and identity, and mechanical evaluation count.

Focused G1 tests independently solve the one-cell swelling polynomial, finite current-area surface resistance, and unlike-layer inventory exchange. G3 tests cover the undeformed limit, identical-material layer splitting, and late work/source/non-axial refusals. The five new groups and six retained layered-solid tests passed against the current source and native dependencies. This is focused execution evidence, not a repository-level DSR verdict or experimental validation.

Scope: explicitly opted-in spatial diffusivity, isothermal quasistatic small-strain mechanics, and axis-aligned rectangular common sections within the caller's nodal tolerance. Bending, shear and taper outside that tolerance refuse. Stress chemical potentials, strain-dependent diffusivity, sorption heat, general moving meshes, finite-strain constitutive mechanics and liquid transport are separate capabilities.


### Finite atmospheric water coupled to swelling geometry

`MoistureLayeredSolid::advance_spatial_diffusion_finite_atmosphere` joins the original spatial transport and native solid solve to the original finite-vapor atmosphere balance. Each endpoint-humidity trial starts from the same accepted moisture field. Its actual axial stretch supplies the bulk `1/lambda_x^2` metric, and its actual current area supplies surface and contact transfer. The atmosphere owner alone accepts water uptake, reconstructs the gas, and checks water, carrier, and total-mass conservation.

`MoistureSpatialAtmosphereControls` bounds complete transport evaluations and actual native solid solves cumulatively across every humidity and nonlinear geometry trial, including discarded trials. The solid allowance is reserved before each native solve. Reports retain original transport/source/geometry identities and cumulative work counts; failed admission, budgets, balances, or cancellation publish no changed input.

This is an isothermal, quasistatic small-strain extension for the existing admitted separable prism. The atmosphere retains its fixed gas free volume: specimen expansion is not deducted from that volume. Displacement-volume coupling, gas-pressure traction, sorption heat, and pressure work are outside this API. Five public-API regression groups cover independent coupled inventory/geometry equations, finite current-area surface resistance, zero-swelling reduction, sealed and wet/dry finite reservoirs, cumulative trial budgets, source/geometry refusal, replay, and cancellation. Rustfmt parsing and source/API review passed; native execution of these new integration tests remains unverified in this environment.


### Source-bound coating admission for complete forced shell response

`PreparedShellResonance::harmonic_factor_with_coating_strain` pairs the complete forced factor with its own original shell or assembled-shell strain owner. Each accepted response retains all selected material memory and independent load states, the source-band checks, signed physical work and original equation residuals. The actual coherent complex displacement must pass every original film's sufficient strain bound over the whole coated triangle, both film faces and the complete oscillation cycle. No phase sampling, modal normalization, foreign observer or copied geometry enters the check.

The caller admits separate physical-observation and complete strain budgets. A successfully completed strain query consumes one allowance even when its sufficient envelope refuses the response; failed solves and failed strain queries consume none. Exhaustion refuses before solving. Reuse leaves the original factor, immutable sources and accepted histories unchanged. This admits only the modeled coating envelopes: the original shell regions author no substrate strain limit, and unmodeled preload or removed static strain, material failure and continuum convergence remain outside the claim.

Two focused regressions in `coating_dynamic_response` cover a phase-interior exceedance under actual dynamic-load forcing and an original assembly with distinct film limits, including source identity, response parity, refusal, reuse, cancellation and resource exhaustion. The current shared `fs-phs` owner and the actual shell source-owner harness compiled; all 14 tests in the focused `coating_dynamic_response` target passed, including both new admission regressions. The harness contains the actual complete shell tree and original supporting physical owners; unrelated application modules are omitted. This is not a full workspace Cargo gate.


### Regional gas sources for continuous rational modal radiation

`RationalModalLine::new_regional` and
`PreparedRationalModalAnalysis::new_regional` consume the original
`assemble_lossless_regional` source with one authored `GasState` per uniform
axial element. Every gas state, local density and sound speed, consistent mass,
stiffness, transverse cutoff, and inlet/outlet reference survives. Repeated
uniform states reproduce the uniform path; no averaged gas, effective impedance,
resonance correction, or eliminated physical coordinate is introduced.

The original rational mouth keeps its exterior gas, signed radial basis, source
band, source clock, and accepted pole history. Runtime analysis uses the same
retained regional source without consulting or changing accepted history.
Continuous response, pole refinement and spectrum discovery retain the complete
joint descriptor and original mouth multipliers. Existing preallocation,
cancellation, source-pressure, pressure-manifold, and transactional work checks
remain in force.

These regions are frozen at one common static pressure, with interfaces on the
authored element grid. Refinement must preserve each physical interval. Mean
flow, mixing, thermal/species evolution, and causal higher-mode viscothermal
losses are outside this lossless interior model. Rational mouth work remains
signed exchange; it is not reclassified as stored energy or heat.

Native compilation of the edited runtime and analysis with the original mixed/zero-direct owners and production physical/numerical dependencies succeeded. All 10 focused `rational_modal_analysis` integration tests passed, including five regional tests for uniform equivalence, local cutoff/transfer convergence, loaded finite-pole convergence, exact-clock work/state subdivision, and source/refusal/history preservation. The bounded harness supplies BEM producer metadata for type linkage only; BEM producer execution and a full `fs-couple` Cargo gate are not claimed.


### Original constrained attachments in complete shell resonance

`ShellResonanceDescriptorLoad` and `prepare_coating_resonance_with_descriptor_load` on both material shell images and material shell assembly images retain every original shell/material coordinate, load differential state, and algebraic reaction. The shared physical interconnection and descriptor harmonic/pole owner preserve the original signed force/volume-flow maps and internal work. No constraint condensation, artificial inertia, or duplicate eigensolver is introduced; the existing ordinary-load path is preserved.

Load observations use the original effort coordinates `[Q_diff*z_diff, z_alg]`. Reactions participate in the complete mutual resistance and port effort, while carrying zero stored energy and no derivative term. Every differential and algebraic load equation is independently checked. `load_differential_range()` and `load_algebraic_range()` identify the original coordinate blocks. The same complete load observer feeds forced response and the existing original coating-strain admission, including assemblies and source identities.

Source-band intersections, complete coordinate/payload/work budgets and cancellation remain required. Singular or redundant constraints refuse without regularization. Finite-pole calculations do not newly certify descriptor index, infinite-spectrum completeness, source-basis convergence, or experimental continuation accuracy. Four regressions in `coating_dynamic_response/descriptor.rs` exercise constrained inertance, mutual algebraic resistance with all material memory, all-differential ODE parity with coating admission, and original assemblies with resource/support/cancellation refusals.


### Original rational radiation coupled to hereditary solid resonance

`solid_rational::PreparedSolidMimoResonance` borrows the admitted original `PreparedMimoLoad` and materializes its exact continuous real A/B/C/D/E under a separate allowance. All source coordinates, mutual terms, source band, clock, accepted history and continuous/sampled report classifications survive. The complete signed force map couples original solid velocity to original source effort through `ViscoPoleStateSpaceLoad`. Dimension/map checks precede realization allocation; the solid's complete-state, work, memory, source-band and physical residual owners remain in force.

Harmonic forcing and free modes do not advance or reset the borrowed accepted history, use the digital mean-direct matrix, or infer source heat/storage from passivity reports. Signed per-port boundary work is explicit. Without an original physical storage witness, complete total-storage and total-dissipation accounts remain `None`; body storage, material heat and body-plus-boundary-work balance remain enforced. The generic solver owner is reused unchanged, and original numerical-continuation, physical-source and spatial-accuracy limits remain applicable.

### Scheduled molecular leakage during finite-gas performance playback

`AperturePerformance::with_finite_gas_effusion` attaches an owned
`FiniteGasEffusionProgram` to the existing finite transport owner before
playback. Fixed, source-qualified interfaces connect its original finite gas
inventories; `GasEffusionControl` schedules their exposed areas on the accepted
acoustic source clock. The existing nozzle program remains usable and may be
empty. Molecular apertures and nozzles must describe distinct physical paths.

Both directional species beams and their molecular crossing energy enter the
same original-inventory update as other paths. The same source/thermal/acoustic
transaction rebuilds each gas once, updates the source from its actual gas,
pays parameter work and checks final qualification. Candidate areas and control
cursors commit only with an accepted source sample. Late physical refusal or
cancellation preserves accepted gas, source state, areas, controls and output.
Original interface definitions are borrowed during preview; enabling this
feature adds one bounded vector of borrowed path descriptors per sample.

The source must satisfy the original fixed-thin-aperture, free-molecular and
negligible-reservoir-drift assumptions. Area scheduling exposes qualified pores;
it does not change pore diameter, add moving-wall work, reinterpret the sound
source's Bernoulli aperture or establish a hard real-time guarantee. Focused
tests compare every accepted state/report with direct native transport, cover
composition and cooling, closing/reopening, shared control/work budgets, cold
preflight and atomic late refusal/cancellation. Rustfmt parsing and source/API
review passed; native test execution remains blocked by missing constellation
sibling checkouts in the current environment.

### Reacting native liquid on the original gas transport network

`GasReaction::transporting_evaporating_liquid_candidate` connects the original
complete gas transport snapshot to the existing radiative reactive phase and
conduction owner. The exact transported atmosphere supplies chemical inventory
and caloric state. Chemistry, phase, finite-wall heat and native liquid flow
share one retained timestep; the original complete gas map, liquid and walls
publish together after the actual final map passes transport driving-force and
source admission.

The nonreactive and reactive entries share one network/candidate/accounting
implementation. Their generic report retains the original coupled report, including reaction
stoichiometry. Chemical constituent changes
are not treated as phase or transport mass defects. Original total mass and
source-relative energy residuals join with absolute local defects before the
whole-interval gate. Formation energy, latent heat and transferred partial
enthalpy receive no additional heat credit.

Network endpoint and chemical caloric controls remain separate explicit stage
inputs. Reactive mass/energy tolerances bound the enclosing interval; its
temperature-departure limit remains relative to the transported atmosphere.
This is an explicit first-order transport/chemistry/phase/conduction split,
with fixed gas volumes and the existing geometry, source and attached-history
restrictions. It adds no simultaneous reaction-transport solution, compression
work, jet momentum or combined entropy/mechanical-energy theorem.

Authored regressions in `tests/evaporating_liquid_gas/transport/reactive.rs`
check transported composition reaching chemistry, complete species changes
against original stoichiometry, whole-inventory thermal energy, finite walls,
radiation and pressure, chemical driving-force reversal rollback, and
empty-path/caloric-control agreement with the original reactive owner.

### Finite-vessel effusive recoil with inertial exhaust accounts

`rocket_blowdown::effusion` connects original angular-qualified molecular
outlets to the native variable-mass exhaust owner. `GasRocketState` retains one
gas inventory and one dry mass; simultaneous outlets spend their original
species allowances and share one caloric reconstruction. Their accepted
species losses determine molecular momentum and energy without resampling an
endpoint flow or adding a continuum pressure-thrust term. Opposed outlets can
cancel recoil while the reservoir continues to cool and fractionate.

Explicit unit outward normals share the state's inertial frame. The exported
momentum and energy include the changing vessel velocity. Accepted-stream rate
rounding and the remap to the actual gas inventory are included in absolute
momentum and energy budgets before either endpoint is published. The raw native
report remains available; corrected physical exports belong to the coupled
step. This closes the gas, vessel and exhaust account without subtracting large
formation-energy baselines. Cancellation or a late physical refusal leaves all
accepted input state unchanged.

The scope is a free translating, locally mixed vessel with fixed inertial
aperture normals and an absorbing vacuum. Rotation, gravity, external work,
slosh, wall accommodation and finite receivers require other owners. The gas
flux is held on each explicit step; changing temperature and composition need
caller-controlled time refinement. Focused G0/G1/G3/G4 fixtures cover analytic
half-pressure recoil, mixtures, opposed apertures, inertial-frame covariance,
full mass/momentum/energy accounts and atomic refusal, including loss of small
recoil under an enormous common boost. Native execution remains unverified
because the required constellation siblings are absent from this checkout.

The `effusion_recoil` example accepts the existing bounded synthetic gas file,
dry mass and authored thin-aperture geometry. It emits accepted gas state,
position/velocity, species fractions and cumulative inertial exhaust accounts
through the same public physical owner; it does not introduce another solver.

### Finite-inertia guided walls on original free-surface liquid momentum

`free_surface_wall::GuidedFreeSurfaceWall` retains the complete original
material-bound D3Q19 free-surface liquid and a native finite-mass rigid state.
One fully flooded planar wall translates along an explicit unit tangent guide,
with no normal displacement or spin. Each accepted step solves the shared
midpoint velocity against the ACTUAL native wall impulse and an optional held
guide force. One native rigid-wall pilot at zero spin supplies the original
reflection's affine scalar resistance; one native replay at the resulting
midpoint supplies the accepted populations and complete wrench. An already
balanced pilot needs no replay. The authored midpoint-speed interval remains
an admission bound. No fitted drag, added-mass coefficient or replacement fluid
algorithm is used, and the final native impulse must independently close.

All trials start from the same accepted fluid. A new immutable native wall
candidate API shares exactly the existing moving-wall admission and receipts,
avoiding an additional full-population clone per trial. Complete domain cells,
aggregate trial-cell work and trial count are bounded before the respective
scan or allocation. Only the accepted candidate advances the liquid, rigid
momentum, midpoint displacement, elapsed time and accepted-step count together.

The native rigid impulse receipt, actual midpoint mismatch, full angular
momentum and liquid boundary work must meet separate physical tolerances.
Small representable impulse/work
defects are retained without subtracting a large kinetic-energy baseline.
The raw guide reaction closes native numerical momentum; pressure-perturbation
and explicitly pressure-referenced physical support forces AND couples remain
separate. The full native link torque is taken about the free half-drift COM;
the guide cancels that couple to hold zero spin. Liquid, guide and applied force
retain their orbital moment about the caller's fixed balance origin. Actual
stored endpoint increments close that angular account without subtracting large
absolute angular momenta or allowing numerical pressure to inflate its budget.
An absolute load is not invented when the source has no pressure reference.

This guided translation owner does not evolve moving occupancy or wall
rotation. Its native link-moment and affine response owner is the SAME one used
by the rotor and joint planar wall paths; no new reflection kernel is added.
Original isothermal VOF viscosity, forcing, capillarity and phase
conversion keep their numerical/physical limits; boundary work balance does
not turn that discretization into an energy-conserving thermal fluid. The wall
angular account does not certify the liquid's other boundaries or VOF conversion.
Six
integration regressions cover independent D3Q19/finite-inertia momentum and
work oracles, per-link moments at an off-center COM, shifted angular origins,
bitwise native population replay and continuation, pressure-datum force/couple
shifts, bounded geometry/work/trials, cancellation and late native refusal.

## Finite native liquid-to-liquid momentum transfer

`liquid_transfer::transfer_interface_mass` atomically moves a finite packet
between two original `MaterialFreeSurfaceEnthalpy3` interface inventories at the
same accepted physical time. Explicit Cartesian grid origins and validated
grid-to-world orientations must map both cell centers to one point, and both
trajectories must retain the same borrowed
surface source, liquid temperature law and caloric curve. Native withdrawal
carries the actual donor velocity with zero actuator work; native receiver
mixing changes its populations and deposits the physical relative kinetic loss
as enthalpy once. Receiver credit uses the donor's represented withdrawal.
Actual two-owner mass, linear momentum and enthalpy-plus-kinetic residuals must
fit the supplied native budgets before either owner publishes. Donor velocity
is rotated into receiver axes before native mixing, and both actual momentum
increments are rotated into world axes before the joint balance. Native
participant reports keep their original local axes. Identical grid axes skip
the velocity round trip. These fixed coordinate rotations do not introduce
moving-reference-frame work, resample either lattice or deposit extra heat.
Cancellation,
source/clock/geometry mismatch, capacity or speed refusal preserves both states.
The original timelines continue from the accepted fills, motion and enthalpy.

The finite liquid transfer prepares one native lattice/enthalpy candidate per
participant and retains unchanged material and coupled-history bindings in the
original owners. `prepare_interface_inflow` and `prepare_interface_outflow`
return typed handles that exclusively borrow their originating trajectory;
read-only reports support joint acceptance without stale-state publication.
Dropping a handle discards its proposal. After joint checks and a final
cancellation poll, both `publish` operations are infallible and allocate no new
state. The existing immediate inflow/outflow APIs use the same preparation and
publication. This removes the coupling wrapper's extra full-owner clones while
retaining the original native conservation and admission arithmetic. Focused
tests cover discarding both prepared proposals after cancellation and exact
parity of published reports/states with the immediate native entry points.

This is coincident-interface transfer between overlapping representations, not
a resolved connector, pipe, jet, phase change or pressure-work model. It adds no
compressible internal-energy, capillary-energy or whole-flow energy theorem.
G1/G3/G4 integration tests independently read complete native mass, linear and
angular momentum, and enthalpy-plus-kinetic inventories, and exercise unequal
velocities, mixing heat, independent grid rotations and world angular momentum
about an offset origin, datum/boost covariance, chained continuation and atomic
refusals. These are synthetic conservation tests, not experimental validation.

### Effusive vessels under gravity and external force

`GasRocketState::advance_effusion_with_environment` retains the same original
molecular outlets, shared species withdrawals and one gas reconstruction while
the native physical-exhaust owner advances explicit uniform gravity and another
held inertial force. Gravity follows the changing vehicle mass; the other force
is held in newtons. Their native integrated impulse and mechanical work enter
the coupled momentum/energy balances once. Closed or absent outlets preserve
the gas and follow the same loaded constant-mass trajectory. The original
`advance_effusion` entry delegates with exactly zero loads.

Gas is the authoritative starting and ending propellant inventory. A local
native starting copy uses its actual gas mass; an existing native/gas mismatch
must pass the mass budget and its remap's momentum/kinetic magnitudes consume
the corresponding physical budgets. Endpoint remapping and accepted-stream
rate rounding retain their explicit signed corrections. Unresolved inventory
rounding during gravity is bounded separately: an absolute mass uncertainty
epsilon_m reserves `epsilon_m*abs(g_j)*dt` per momentum component and
`epsilon_m*|g|*dt*sqrt(mean(|v|²))` in work. The native powered trajectory supplies
that velocity moment. These reserves accompany the signed residuals and do not
claim an exact interpolation of rounded gas mass or become fictitious heat.
No-outflow coast has no mass-history uncertainty. A late budget, source or
cancellation refusal publishes neither owner.

The compact gas remains locally Maxwellian and mixed; internal gravity
stratification, rotating apertures and resolved feed flows remain outside this
model. Loads and normals are fixed during a call. Six added regression groups
cover explicit/free parity, independent changing-mass trajectories, complete
external-work and Galilean accounts, loaded closed-outlet motion, canonical
inventory budgets, extreme-scale gravity reserves and atomic refusal. Rustfmt parsing and source/API review
passed; numerical execution remains unverified while required constellation
sibling dependencies are unavailable.

## Directional momentum exchange in solids with material history

`inelastic_impact::multiple::elliptic` extends the native plastic, Maxwell,
thermal Tet4 and thermal Hex8 manifold entries with authored directional
friction. Its four `*_mixed_elliptic_contact_candidate` functions retain the
same material-point traces, full consistent cross-point mobility, rigid
full-tensor response, support handling and single summed impulse update per
body. A shared law callback selects the original
`fs_contact::frictional_multicontact::elliptic` owner without copying mechanical
assembly, contact search, material state or endpoint audits.

Each contact supplies two nonnegative principal coefficients and its first
principal direction in the admitted tangent frame. The complete work-conjugate
ellipse transformation retains normal/tangent and intercontact blocks. Only
original physical impulses and velocities enter the bodies; transformed
determinacy witnesses remain labeled in the report. Equal coefficients recover
the original circular response exactly, both zero recover frictionless contact,
and a single zero coefficient retains exact free slip without small substitute
traction.

Unique transformed impulses determine their physical image and admit the
original signed work at driven solid supports. A velocity-only witness still
refuses unresolved moving-support work; a representative redundant impulse is
not promoted into authoritative local work or heat. Native caloric storage,
plastic history, every Maxwell branch, configuration and time remain unchanged.
Impact dissipation remains separate from material dissipation and is not
automatically deposited as heat.

Regressions extend the existing nonzero plastic/Maxwell-history fixtures with
circular/frictionless parity, collapsed ellipses, the physical principal-slip
law, heterogeneous-contact permutation, joint texture/tangent-frame rotation,
driven-support uniqueness and redundant-row refusal, retry and cancellation.
The existing bounded solver, source-domain, native mass-solve and physical
acceptance gates remain authoritative. This is a frozen-geometry mechanical
event; evolving texture/source laws and sustained-contact stepping retain their
separate original owners.

## Retained gas drift in native phase transfer

`liquid_gas::evaporation::momentum::advance_liquid_with_moving_gas` advances the original Stefan/enthalpy phase solve and one native VOF interval with explicit retained gas bulk momentum. `MovingGasParcel` requires an authored Cartesian velocity and preserves the original finite gas composition, internal energy, volume and source axes; bulk kinetic energy is derived separately. Evaporation carries the actual liquid mean velocity into gas; condensation carries the retained gas mean velocity into the original native complete-mixing owner. The receiver's physical relative kinetic loss heats that same receiver once through its original caloric law. Gas momentum updates by the reciprocal accepted mass impulse, and actual represented liquid plus gas source momentum and kinetic-plus-caloric energy must close within the declared budgets.

Optional sensible heat retains the original source solve; optional pressure uses the actual mixing-heated gas EOS. Neither pressure traction nor a molecular vapor-jet recoil is inferred as an additional phase impulse. Source domains, latent-datum compatibility, temperature-increment limits, native capacity and speed bounds remain authoritative after mixing. The complete private gas/liquid candidate publishes atomically after one retained flow timestep; zero-flux sensible heating preserves the same mechanical account. This is a homogeneous mean-drift, first-order phase/mixing/flow split, with fixed connected gas volume. Its closed source ledger excludes subsequent VOF wall/pressure forces and numerical kinetic truncation; those are not mixing heat.


## Two moving superfluid pistons with authored fluid and body forces

With `superfluid`, `superfluid::piston_pair::advance_balanced_piston_pair`
composes the original Landau channel with two original native rigid bodies.
Both ends move; there is no stationary support absorbing unreported momentum.
Left/right body states, mass properties and balanced-pressure references are
explicit. Bodies share one world axis and remain centered, spin-free and
within the declared small-displacement reference configuration.

The original reciprocal two-end midpoint response supplies a mass-scaled
symmetric 2-by-2 Cholesky solve for both velocities. A force-only anchor
preserves uniform translation and common acceleration without artificial
blocked-pressure loads. Native `RigidBodyIntegrator` advances positions and
momenta without changing native masses or inertias. Each endpoint fluid
half-cell carries separate momentum and kinetic storage, distinct from the
resolved channel interior.

Authored component body forces and independent outward piston loads retain
their physical impulse and work. Reports distinguish pressure impulse, actual
native momentum change, fluid contact impulse, external body impulse,
endpoint-fluid source impulse and actual midpoint work. Whole-system momentum
and energy are checked before any owner publishes. Friction heat is exported
once. Nonlinear advection, moving meshes, rotation and temperature-driven
material updates remain outside this reference-grid model. The existing
single-piston API shares native admission and motion preparation helpers.

`superfluid_piston_pair` independently checks isolated translation/recoil,
common acceleration of all fluid mass and both bodies, spatial fluid/body
forcing (including disparate native masses), and complete rollback on refusal
or cancellation. Its cumulative balances include both fluid half-cells.
All five pair, five existing piston and seven heat tests passed native
source-module validation; a full workspace Cargo build is not claimed.

### Adaptive source-gas recoil trajectories

`GasRocketState::advance_effusion_adaptive` advances a complete requested
interval through the original loaded effusion owner. Each attempted interval
compares one full step with two half steps, accepting only the two half steps.
The gas discretization remains first order. Explicit dimensional tolerances and
a relative tolerance bound the estimated local difference in temperature,
pressure, every original species inventory, position and velocity. Mechanical
relative scales use velocity increments and displacement beyond initial drift;
a large coordinate origin or common boost cannot enlarge those scales.

The controller deterministically halves rejected intervals and may double
small-error accepted intervals. Only the original typed gross-withdrawal refusal
is retried at a shorter step. Source-domain, frame, physical balance and other
owner failures propagate. Minimum/maximum steps, attempts and accepted interval
counts bound all work and retained history; each accepted interval stores two
native reports. The final interval is shortened or split into a representable
pair without silently discarding time. Nonadvancing half-step clocks, unmet
accuracy at the minimum step, exhaustion or cancellation publish no partial
trajectory endpoint.

Compensated totals contain only accepted molecular mass, inertial exhaust
momentum/energy, external impulse/work, signed conservation residuals and
nonnegative gravity reserves. Original physical budgets apply per accepted
native half step. Adaptive tolerances estimate local truncation error, with no
certified or fixed global error budget. The complete accepted reports remain
available for consumers needing their original source and species accounts.
Outlets, normals and environment are held over the requested interval; the
controller does not infer control events, resolve rotating gas, or replace the
original thermochemical and mechanical laws.

The `effusion_recoil` example retains fixed stepping by default and accepts
optional adaptive mode plus world gravity and force vectors. Its CSV includes
all three position, velocity and momentum components, cumulative external
impulse/work, gravity reserves and accepted/rejected step counts. It uses the
same explicitly synthetic gas fixture and authored aperture qualification.


## Atomic two-piston superfluid heating

`superfluid::piston_pair::advance_balanced_piston_pair_with_heat` composes
the same reciprocal fluid/two-body proposal with the original finite thermal
timeline. Active unique destinations receive the current interval's per-cell
friction joules exactly once. Original thermal sources, events and contacts
advance through `DepositionSimulation::advance_to`; the original caloric
owner accepts heat through `apply_heat_j`. No alternate solver, repeated
mechanical step, copied damping update or cumulative-heat replay is used.

The complete storage gate adds the pair's mechanical-fluid residual to the
original thermal march, caloric intake and friction-transfer residuals. It
counts neither friction nor the original wave residual twice. Channel fields,
heat and clock, both native body positions/momenta, and thermal cells remain
unchanged on caloric refusal or final cancellation. Publication follows the
original private-candidate boundary, with only infallible assignments after
acceptance. Frozen material coefficients and the first-order endpoint heat
split are unchanged; native thermal storage conventions remain authoritative.

The seven-test `superfluid_piston_pair` suite includes repeated driven
heating, independent per-cell intake and whole-system source/storage balances,
caloric refusal, late cancellation and exact retry. Thirty-six distinct
channel/coupling tests passed native source-module validation in this increment.

Run the complete example with:

```sh
cargo run -p fs-couple --features superfluid --example superfluid_piston_pair
```

The example uses the original sourced 1.8 K helium reference, two unequal
native pistons, and physical reference-cell masses. Its explicit local
isochoric caloric tangent uses the reference cv and permits at most 1 mK of
dissipative heating; it is not a nonlinear helium EOS or coefficient-feedback
model. Geometry, excitation and mutual-friction coefficient are authored demo
inputs. At default 64 cells, 926 intervals cover 5 ms. The native executable
completed independent cumulative momentum, represented-energy, heat-transfer
and clock checks. Observed relative energy defect was 3.66e-13, momentum defect
9.84e-22 N s, and heat-transfer defect 1.85e-26 J. These are numerical balance
results for this fixture, not accuracy bounds against an experiment. Full
workspace Cargo verification remains unavailable under the current resource
limits.


### Prescribed axial rotation in material film surfaces

`ResolvedImmersedBodies::with_spinning_prescribed_planes` accepts one signed
angular speed in rad/s per unchanged `PlaneFilmSi` source descriptor. Rotation
is about each plane normal through its translating source point. Rates convert
with `Omega_lattice = Omega_SI * timestep_s`, with no length factor; local rigid
wall velocity uses the original cell-size/time conversion. Source SI rates are
retained separately from accepted native plane state. Read-only observations
expose current SI angular speed and local velocity without mutating the owner.

All film resistance, equal/opposite linear and angular impulse, and translational
plus rotational wall work stay inside the original simultaneous native
fluid/body solve. Translation, material rebind and cancellation/retry retain the
same plane spin and accepted clock. Zero rates follow the original prescribed
plane path. Existing `PlaneFilmSi` construction remains source-compatible, and
the spatial fluid-force forwarding API is preserved. No separate force update,
post-step impulse, wall trajectory or finite carrier is introduced.

The resolved fluid boundary must independently match the prescribed motion.
The native positive-gap, local-speed, missing-resistance and asymptotic rules
remain authoritative. This correction does not supply finite axial torsional
resistance, tilting planes or solid contact. The native six-test production-source
suite exercises the coupled resistance/work model. Three material integration
tests cover SI conversion, actuator work and momentum, shared-clock rebind,
zero-rate continuation and refusal/retry; they parsed but were not executed in
the delivery environment. No full material/grid integration pass is claimed.

## Finite axial rotor on the original flooded free-surface wall

`free_surface_wall::rotor::FreeSurfacePlanarRotor` retains one original material
liquid and native rigid-body state for continued fixed-center rotation. The
plane normal must be a principal axis of the original inertia. The guide fixes
the center of mass and carries transverse torques; the original MBD integrator
updates axial spin and attitude. Unsupported initial motion is refused rather
than projected onto the supported subspace.

A pilot native reflection measures the axial angular impedance. Its candidate
storage is discarded before the final replay. The affine midpoint equation
couples actual liquid torque, finite rotor inertia and held external torque
analytically; exactly two native liquid candidates and one native body step are
used per accepted liquid tick. The final native link wrench, actual represented
body momentum, midpoint velocity, attitude, kinetic increment and work must all
pass the declared bounds. No apparent inertia replaces the physical inertia.

Reports distinguish raw numerical-pressure reactions, physical perturbation
loads, optional absolute pressure loads, transverse guide moments and actuator
work. All moments use the fixed center. Guide work is measured, not silently
discarded. The liquid remains isothermal: mechanical boundary work is not an
invented heat deposit. Work, angular advance and full-grid trial cells are
explicitly bounded; cancellation and late endpoint-speed refusal publish neither
owner. Fixed flooded geometry excludes normal translation, tilting, finite
plate edges, cover/vacate and moving contact lines. Time and spatial refinement
remain necessary for trajectory accuracy, including very light rotors.

Focused integration regressions specify independent light/heavy rotor response,
native liquid angular momentum and body attitude, driven continuation,
off-center pressure-reference guide loads, and atomic early/late refusal.

### Bidirectional molecular sources on a finite translating vessel

`GasRocketState::advance_effusion_exchange` composes prescribed molecular-source
ports and true-vacuum outlets on the SAME original gas inventory, gross species
withdrawal allowance and single caloric reconstruction. The original reciprocal
effusion law and explicit uncollimated angular qualification supply each beam;
the native particle-exchange owner advances the actual vessel with simultaneous
positive intake and discharge. Transferable mass can grow, shrink or remain
constant. Incoming mass cannot fund another port's withdrawal during the same
finite proposal, and complete original species/source admission remains binding.

An `EffusionReservoirPort` supplies a prescribed thermochemical source that
remains instantaneously comoving with the vessel-mounted aperture throughout
the held interval. Its declared source volume does not imply a depleted finite
reservoir. A stationary inertial ambient, two independently moving vessels,
shifted Maxwellian arrival distributions and particle flight geometry are not
inferred. Exterior wall-pressure and support loads remain explicit external
forces. The two opposite molecular beams BOTH contribute to transmitted normal
stress; equal counterflow does not erase force or imply a pressure-difference
wall-force model. No additional continuum pressure thrust is added.

Each original accepted species mass carries its original crossing energy and
qualified directional momentum. Incoming relative kinetic/internal energy
enters the gas caloric solve exactly once, without a second thermalization heat.
Canonical per-port receipts and their totals retain separate incoming/outgoing
mass, relative momentum/energy and inertial momentum/energy. Stable signed net
mass and energy receipts use original boundary coordinates before gross reports
are subtracted, preserving heat exchange beside large formation references.
Both stream directions use the SAME native time-mean velocity and squared speed
along the changing-mass trajectory, including held gravity and external force.

The gas remains the authoritative inventory. Initial native mass normalization
reserves its momentum/energy cost; gross source-to-rate differences are bounded
separately, and final native mass remapping contributes its actual signed
momentum and kinetic correction. Source-relative gas closure and native kinetic
closure combine without subtracting formation energies. Conservative gravity
mass-history reserves consume the remaining physical budgets without claiming
an unqualified interpolation of gas reconstruction rounding. Refusal or final
cancellation publishes neither a gas endpoint nor a mechanical endpoint.

Six focused genuine-source fixture groups cover hot-source filling, simultaneous
vacuum discharge and zero-net-mass thermal counterflow, selective addition of a
previously absent species, common boosts and loaded full accounts, shared gross
withdrawal despite replenishment, and source/domain/frame/cancellation rollback.
Production and fixture files passed Rustfmt parsing and direct source/API review;
native Cargo execution remains unavailable with the current missing constellation
sibling manifests. No experimental accuracy or full workspace pass is claimed.


### Adaptive molecular feeding with evolving species and shared time control

`GasRocketState::advance_effusion_exchange_adaptive` advances prescribed
comoving reservoir ports and vacuum outlets over a complete requested interval.
It shares one step-doubling controller with the existing vacuum-only adaptive
entry point. Each trial compares a complete native transaction with two native
half steps, and only the two accepted half-step candidates become history.
The original coupled gas, momentum, energy, source qualification and shared
gross-withdrawal checks remain active in every proposal. Typed withdrawal
limits can request a shorter trial; other physical/source refusals remain
atomic failures.

The local error comparison follows temperature, pressure, individual species
mass, velocity increment and displacement beyond the original inertial drift.
Its species axis is the canonical union of the vessel and geometrically active
source catalogs, retaining original caloric identity and assigning zero mass
to absent species. A newly arriving species therefore participates immediately
in error control. The shared axis uses the original thermochemical bound of
128 species and bounded stack scratch, preserving closed-vessel coasting;
active apertures retain their separate original 32-species qualification.
Geometrically closed source ports do not inspect inactive source catalogs.

`EffusionExchangeTrajectoryStep` retains original per-port receipts in its
accepted native history. Compensated whole-interval accounts separately sum
gross incoming/outgoing streams, authoritative net mass and molecular energy,
inertial momentum/energy, external impulse/work and physical residual/reserve
reports. Rejected trials contribute no transfers or totals. Exact final-time
handling, minimum/maximum intervals, attempt/accepted-history limits and late
cancellation are owned by the same controller for both entry points.

The finite-gas scheme remains first order. Its step-doubling estimate controls
local error and does not certify a global error bound; physical conservation
budgets apply to each accepted half step. Source thermochemical states, aperture
normals and external loads remain held across the requested window. A stationary
ambient, rotating vessel or independently moving finite receiver is outside
this common-frame source closure.

Focused genuine-source regressions cover introduced species, temporal
refinement, accepted-only gross/net accounting, vacuum-only parity, original
closed-mixture limits, bounded clocks/work and atomic source/cancellation
refusals. Rustfmt parsing and source/API review are available; native compilation
and execution remain unverified with the missing constellation dependencies.

### Endpoint-sampled fluid forcing on one material dry body

`ResolvedImmersedBody::step_forces_si_with_endpoint_fluid_force_scale` advances
the original sourced positive dry mass/tensor, fluid populations and accepted
clock with a smooth prescribed spatial load history. The dimensionless endpoint
scale multiplies the original SI field selected by
`rebind_fluid_force_field_si`; the uniform fluid load stays held. Existing
source density, cell size, timestep, authored body forces/couples and reduced
gravity retain their original conversion and ownership.

The native dry-body solver uses the next Guo source frame for actual streamed
inventories and open faces, and subtracts the same trapezoidal source from the
retained fictitious interior. Complete native momentum, work and boundary
receipts return unchanged. No fallible conversion follows native publication:
the body, fluid, force and material clock are adopted together. Refused source
data, numerical budgets and cancellation preserve every accepted owner.

A genuine source fixture compares ramp, zero crossing and reversal against
independently constructed original native motion, reconstructs the SI volume
impulse, and checks physical source identity and rollback/retry. Rustfmt parsing
and source/API review passed; compilation and execution remain unverified.
This entry supplies no extra pressure-gradient body traction, endpoint-dependent
constitutive load, finite heat storage or complete coupled temporal-order claim.

### Endpoint-sampled fluid forcing on material assemblies

`RetainedFluidBodies::step_forces_si_with_endpoint_fluid_force_scale` and
`ResolvedImmersedBodies::step_forces_si_with_endpoint_fluid_force_scale` carry
the selected SI spatial force field through the original simultaneous
sphere/ellipsoid owner, including sphere-pair films. The uniform fluid force
stays held; original reference density, dx/dt, authored center forces/couples
and reduced weight retain the existing conversions. Physical dry assemblies
measure endpoint Guo inventories and subtract each midpoint fictitious-interior
trapezoidal source once within the full coupled solve. Returned native receipts
retain explicit lattice units and shared boundary/source ownership.

All SI conversion and next-clock admission precede native execution.
Populations, all body poses/momenta, endpoint force and accepted clock advance
together; no fallible post-publication conversion was added. Held entry points
retain their original path. The moving plane-film and finite/shared-carrier
continuation below extends this source history to those owners. The native implementation passed the
focused 39-check selection documented in fs-lbm; these two typed source-card
wrapper implementations were reviewed and rustfmt-parser checked, but their
complete constructor/dependency graph was not executed.


### Endpoint force histories on material film surfaces and carriers

`ResolvedImmersedPlanes`, `ResolvedImmersedSurfaces` and
`ResolvedImmersedSharedSurfaces` expose
`step_forces_si_with_endpoint_fluid_force_scale` through their original native
transaction. The selected SI spatial force receives the next dimensionless
scale; uniform fluid force, reference density, dx/dt, original material inertia,
authored center loads and reduced-weight ownership remain unchanged. Shared
registries convert and apply one complete load per carrier, independently of
patch count.

The prescribed-plane owner prepares every next position before native
acceptance. Finite/shared owners prepare and admit their actual coupled
trajectories inside native candidates. Populations, dry-body and carrier/patch
states, endpoint force and accepted material clock advance together; no fallible
conversion follows publication. Refusal preserves the accepted source and all
clocks/states, while reusable scratch may change. Native receipts retain lattice
units and distinct prescribed work/finite storage. The 16 native checks are
documented in fs-lbm; these three typed SI implementations passed source/API
review and parser checks, but their complete source-card constructor graph
was not executed.

### SI reservoir boundaries on material immersed bodies

`ImmersedDomain::with_boundaries` creates an explicit
`ImmersedBoundaryDomain`. Its six `FaceBoundarySi` entries follow `Face3::ALL`
and admit stationary walls, paired periodic faces, prescribed material-fluid
velocity in m/s, and pressure perturbation in Pa. All eight retained-fluid,
settling, and physical dry sphere/ellipsoid constructors accept this setup
through their shared original construction path. Existing `ImmersedDomain`
literals retain their periodic/stationary defaults. Resolution, allocation,
material inertia, body geometry and time integration keep their original owners.

The common converter uses the original source-owned `MaterialScaling`.
Velocity components use `velocity_lattice`, and the complete prescribed vector
must respect the scenario speed bound. The inverse pressure mapping reuses
`pressure_perturbation_si` and the same reference-one density: in the original
numerical isothermal convention,
`rho_lu = 1 + p_prime / (CS2 * rho_reference * (dx/dt)^2)`.
Zero gauge pressure maps exactly to density one; negative perturbations must
retain positive density. Nonfinite values, unrepresentable scales or densities,
and nonzero perturbations too small to change reference density refuse. This is
not an absolute gas-pressure or physical acoustic equation of state.

`BoundarySpec3::try_validate` exposes the existing native face admission as a
fallible operation, and the original asserting constructor delegates to that
same implementation. SI conversion and boundary admission happen before grid
allocation. Periodic faces stay paired, open faces occupy at most one axis,
and the native finite/tangential/low-Mach checks remain. The immersed SI surface
does not admit prescribed moving walls. There is no live topology mutation.

Accepted construction feeds the same original open-boundary immersed solver.
Reservoir mass and momentum remain in its shared native receipt, convertible
through the existing mass/impulse/angular-impulse methods; no duplicate body
reaction or heat destination is introduced. Source rebinding still requires
unchanged reference density and fixed dx/dt, preserving these unit mappings.
The original body halo, pressure-work and thermodynamic limitations remain.

Accepted owners expose `rebind_open_boundaries_si` through retained-fluid,
settling, ordinary dry-body/assembly, prescribed-plane and moving-surface
interfaces. Construction and live updates share the same converter and source
density, dx/dt and scenario speed bound, including after viscosity-only source
rebinding. All face kinds and stationary/periodic data remain fixed.

A successful update changes only held reservoir targets, preserving populations,
physical force-corrected moments, body/surface histories and the accepted clock.
The next common native step owns the actual reservoir mass and linear/angular
momentum transfer exactly once. Invalid controls, topology changes or
cancellation preserve the previous targets. Callers sample their histories at
native tick boundaries; no interpolation, pressure-work or finite-reservoir
storage is inferred. The actual source-backed regression in
`tests/material_immersed/flow.rs` targets positive/reversed/zero inlet histories,
pressure changes after viscosity rebinding, unchanged pre-step physical state,
independent physical momentum/mass balances and atomic refusal.

The focused cases in `material_immersed/boundary/tests.rs` exercise pressure
round trips across scales, legacy defaults, typed pressure/speed/topology
refusal, and six actual native body/fluid steps. Independently authored native
faces reproduce every population, body state and reservoir receipt exactly,
and the prescribed inlet accelerates the immersed body.

Accepted material immersed owners also expose `set_velocity_profile_si` and
`clear_velocity_profile` through the same retained-fluid, settling, dry-body,
assembly, prescribed-plane and moving-surface interfaces, including shared
finite-carrier surfaces. Sampling receives fixed world/grid cell centers in
metres and returns velocity in m/s. Original source-owned scaling, native
regularized reconstruction and the two reusable plane buffers remain
authoritative. No converted volume field is allocated.
`InletFluxSi::MeanNormalVelocity` specifies a signed outward plane mean in m/s;
inflow is negative on either face. This is a prescribed velocity constraint,
not density-weighted mass-flow or exact population-flux control. Mixed
wall/open seams retain the original mean-control refusal.

The source-scaled speed bound applies to the final corrected trace before
the native buffer swap. Source, conversion, low-Mach, corrected-speed and
cancellation refusals preserve the previous complete profile. Updates and
clearing advance no accepted state; the following interval's shared `open`
receipt owns actual reservoir mass and linear/angular momentum once.
Profile diagnostics return the native explicitly named lattice units so no
fallible conversion follows publication. Callers retain coherent-field time,
spectral resolution and physical inflow validation. Focused native tests cover
all face orientations, corrected-speed refusal, source/cancellation retry and
nonuniform shear driving body translation/spin with independent six-component
momentum balance; the full source-backed continuation/clear regression is in
`tests/material_immersed/flow.rs`.

The same material immersed interfaces also expose `set_pressure_profile_si`
and `clear_pressure_profile`, including shared finite-carrier surfaces.
Sampling receives fixed world/grid cell centers in metres and returns gauge
pressure perturbations in Pa. The original source-owned
`density_lattice_from_pressure_perturbation` converts each sample to the native
reference-one isothermal density target. Its positive-density, finite-scale and
representability checks remain binding; no new pressure or Mach law is added.

The common native pressure-profile owner reuses two scalar plane buffers and
publishes only a complete admitted trace. Source/conversion/cancellation refusal
retains the preceding profile. Updating or clearing profiles preserves accepted
fluid/body/surface history, force frames and time. Subsequent original steps
measure reservoir mass and linear/angular transport once in their shared
`open` receipt. Diagnostics retain explicit lattice units, so no fallible
conversion follows publication. These controls do not supply absolute gas
pressure, finite reservoir storage, pressure work or a calibrated inflow model.

Focused original-source SI cases compare physical cell-center sampling and
pressure conversion against independently authored native profiles on all
six faces at two source densities and two cell sizes. Invalid or unrepresentable
pressure, retained source errors, cancellation and successful retry are covered.

### Inertial ambient capture and complete molecular momentum

`GasRocketState::advance_effusion_exchange_with_ambient` adds prescribed
`EffusionAmbientPort` sources with fixed world velocities in the vessel's
declared inertial frame. It extends the SAME vacuum/comoving material and native
mechanical owners. The existing comoving entry point delegates with no ambient
ports; all three port families share canonical IDs, one original gross species
withdrawal budget and one finite-gas caloric solve.

Each ambient samples its original Maxwellian relative to the actual starting
vessel velocity through an `EffusiveDriftInterface`. Its accepted incoming
species masses carry the complete conditional velocity vector, including
tangential slip, and their original internal plus translational crossing energy.
Outgoing gas keeps its original isotropic effusion law. Both directions enter
the same changing-mass native momentum and kinetic-energy accounts; incoming
crossing energy heats the gas once, without an added thermalization term or
continuum pressure thrust. External wall pressure remains a separately authored
load. Source volume does not impose a finite ambient inventory.

The finite native step holds the sampled relative moments, so this remains a
first-order approximation to an inertially fixed source. Subsequent steps must
resample actual relative velocity. Original source identity, active species,
temperature/pressure support, frame and authored relative-speed qualification
remain binding. The final actual relative speed is admitted before publication;
sampled and endpoint checks do not certify a continuous trajectory enclosure.
Geometrically closed ports leave their source distribution and velocity inactive.
Refusal or terminal cancellation exposes neither endpoint.

Focused original-source regressions cover mixed-port zero-slip parity,
independent half-range headwind moments, tangential momentum and energy,
Galilean boosts of BOTH source and vessel, complete physical accounts and
final-drift/frame/source/cancellation refusal. Rustfmt parsing and source/API
review passed; native compilation and execution remain unverified with the
missing constellation sibling manifests.

### Adaptive inertial ambient exchange and executable trajectory input

`GasRocketState::advance_effusion_exchange_with_ambient_adaptive` refreshes each
fixed-world source's relative velocity at every coarse and fine physical trial.
It reuses the common step-doubling controller and accepted-only gross/net
accounts, without a second mechanical integrator or species solver. Active
ambient source species join the original estimator axis even when absent from
the initial vessel. Original source, speed, gross withdrawal and conservation
admission remain active on every native candidate. The source velocity stays
fixed while the vessel accelerates; rejected trials spend no material and emit
no receipts. This is a first-order local error estimate, not a certified global
accuracy bound or a continuous source-domain enclosure.

The existing `effusion_feed` executable accepts `--inertial-source VX VY VZ
MAX_RELATIVE_SPEED`, `--vessel-velocity VX VY VZ`, and `--adaptive` while retaining
the original positional gas/geometry inputs and default comoving fixed steps.
Bounded parsing refuses duplicate, incomplete or nonfinite velocity options.
CSV rows retain separate gross streams, authoritative net accounts, original
species fractions and residuals, and append accepted native substeps, attempted
intervals and rejected intervals. Its explicitly authored synthetic aperture
and drift qualification make no claim of empirical validation.

Focused regressions specify changing-velocity resampling, refinement against
512 genuine native fixed steps, full interval accounts, new ambient species,
inactive foreign sources and late cancellation rollback. Rustfmt parsing and
direct production/fixture/example API review passed; native execution remains
unverified with the missing constellation dependencies.

## Free two-component translation on the original nonlinear ALE owner

With `experimental-stokes-motion`,
`stokes_motion::translating::planar::PlanarTranslationMotion` frees both world
translation components on the SAME original nonlinear translating flow owner.
Each actual native rigid impulse and its pressure/viscous reaction close one
coupled two-component endpoint. The entire original domain follows the body's
finite translation; reference operators and material transport are retained.
No fitted drag, apparent rigid mass or independent x/y fluid update is used.

A bounded local damped Newton search obtains the full cross-coordinate response
from two original flow probes. Authored velocity boxes and dimensional probe
sizes are explicit. The reusable native two-by-two LU owner solves the scaled
response; backtracking uses fixed-base physical residual scaling. Actual probe
velocity changes, including represented cross-coordinate rounding, determine
the correction. Every base, probe and rejected trial counts against the same
native-flow-solve budget. The local method does not guarantee root existence,
uniqueness or global convergence; original native nonlinear gates remain active.

Acceptance closes both ACTUALLY stored body momentum increments, and separately
checks represented body/mesh drift and mechanical work. Springs and dashpots
act in fixed world coordinates. Their stored energy, physical damping, signed
interface/external work and backward-Euler numerical loss remain separate;
numerical loss is not heat. Only pressure/viscous stress acts on the impermeable
body wall; advected momentum and kinetic energy retain their native accounts.
Original component-local pressure datums and material density remain binding.
Orientation is fixed, with the complete zero-work restraint couple about the
moving COM reported explicitly. Initial spin/out-of-plane momentum refuses.

The accepted native ALE candidate publishes once, together with body momentum
and pose, after final cancellation admission. Budgets, bounds, source/datum
refusal, work failure and cancellation preserve both continuing physical owners
and their shared clock. Exterior callbacks describe the translating outer
boundary at endpoint world coordinates; independently stationary obstacles,
rotation, remeshing and a whole-fluid energy/passivity theorem are not implied.

Four focused regression groups specify changing-direction momentum with finite
fluid inertia and common boosts; an independent off-diagonal world-force/spring
oracle and complete uniform-container energy balance; simultaneous open pressure
and tangential viscous recoil; and budget/datum/bounds/final-cancellation rollback.
New source and tests passed pinned Rustfmt parsing. Native compilation and test
execution remain unverified in the delivery environment.

### Conservative molecular exchange between two finite moving gases

`molecular_exchange::step_candidate` joins two ORIGINAL `MovingGasParcel`
owners through one explicitly shared thin aperture. Both gas inventories,
compositions, temperatures and bulk momenta evolve. The common material sampler
supplies two distinct drifting-Maxwellian streams; each sampled species mass,
full conditional vector impulse and complete crossing energy is used with
opposite signs on the two owners. Neither side is an infinite prescribed bath.

Both donors' gross original-species withdrawals are admitted before either
caloric solve. The retained small-drift net mass has its own gross-rounding
check and supplies the common frame transform once. Each original moving-gas
boundary owner uses its actual final mass and represented momentum to remove
the corresponding kinetic change before ONE original caloric reconstruction.
There is no second EOS, fitted heat capacity or added receiver-mixing heat.

Joint admission checks actual total and individual-species mass changes,
componentwise momentum, reciprocal original-reference energy and the combined
internal-plus-bulk energy account. The aperture-frame and laboratory-frame
energy accounts remain distinct, including the kinetic transform of actual
mass/momentum representation residuals. Tolerances scale physical exchange and
relative caloric work, not absolute formation energies or unrelated large
momentum components. Both final source/transmission and relative-speed domains
are admitted without resampling. Terminal cancellation publishes neither gas.

This is a first-order homogeneous fixed-volume gas/control-volume exchange;
time refinement resamples changing gas states. A shared interface does not
infer separated-vessel flight, moving container geometry, remaining-wall
pressure/support forces or a resolved flow field. Equal counterflow retains
normal momentum transfer; the corresponding finite bulk kinetic change is
funded by the original internal energies. Zero geometric area leaves inactive
source distributions and aperture velocity unqueried.

Six genuine-source test groups specify independent finite mass/pressure moments,
zero-net-mass stress, tangential kinetic loss, common gas/aperture boosts and
formation-reference invariance, shared gross withdrawal and newly arriving
species, and atomic source/drift/cancellation refusal. Scoped Rustfmt parsing
and direct API/source review passed. Native compilation and execution remain
unverified because the required sibling manifests are unavailable.

### Adaptive finite molecular exchange with original inventory reuse

`molecular_exchange::trajectory::advance_candidate` advances a complete pair
window through the original fixed-step owner. Each trial compares one full
step with two native half steps and retains only the accepted half steps.
Both finite temperatures, compositions and bulk velocities are resampled in
each proposal. Local error follows both gas temperatures/pressures, every
original interface-catalog species (including zero initial mass) and each
velocity increment. Common inertial drift cannot enlarge the velocity error
allowance. The canonical material inventory is queried once per gas/state;
fixed-step withdrawal and joint mass checks reuse the same public inventory
owner instead of rebuilding it inside per-species loops.

Separate finite limits control initial/minimum/maximum intervals, attempted
intervals and accepted history. A shorter final remainder must advance the
physical clock exactly. Only local error or typed duration/gross-withdrawal
limits request a shorter trial; source, kinetic/caloric and conservation
refusals retain their physical meaning. Cancellation or budget exhaustion
publishes neither partial gas nor partial accepted history.

Whole-window accounts use compensated accepted-only sums for both gross
directions, original relative caloric coordinates, bulk kinetic changes,
signed transferred mass/momentum and mass/momentum/energy residuals. These
residuals never subtract gross formation-energy totals. The estimator remains
first order and local; no global error or continuous source-domain enclosure
is claimed. It converges the shared-aperture operator. Entropy/passivity and
isolated closed-vessel equilibration require the remaining wall/face and gas
transport operators; timestep refinement alone cannot supply their physics.

Focused native-owner fixtures specify temporal refinement, accepted-only
accounts, changing source samples, newly arriving species and bounded work,
minimum-step and late-cancellation refusal. Scoped Rustfmt parsing and source/API
review are available; native compilation and execution remain unverified with
the missing constellation sibling manifests.

### Runnable finite molecular exchange example

`cargo run -p fs-couple --example molecular_exchange --
crates/fs-material/examples/material_gas_spring/mixture-a.tsv --adaptive`
uses the existing seven-column material source loader and one original species
catalog for two finite moving-gas owners. Explicit SI inputs set each initial
temperature, pressure, volume and bulk velocity, the common aperture velocity
and oriented normal, its area and qualified molecular geometry, the maximum
admitted relative speed, interval duration/count and optional native adaptive
advance. Defaults are authored synthetic values. Arguments, source bytes,
species count, total duration and adaptive work are bounded; malformed,
duplicate, unknown and nonfinite options refuse.

CSV rows carry both accepted gas states and species fractions, separate gross
directional and signed net masses, transferred vector momentum, interval and
compensated cumulative mass/momentum/energy residuals, and adaptive work counts.
A refused interval emits no new endpoint. The original source, caloric and
conservation owners remain authoritative in both fixed and adaptive modes.
The example explicitly describes an interface stage: gas bulk velocities are
not moving-container coordinates, no separated-vessel flight model is supplied,
and the remaining wall/face momentum and gas-transport operators are required
for complete closed-vessel equilibration or entropy/passivity conclusions.
Reducing the timestep alone cannot supply those omitted physical operators.

Scoped Rustfmt parsing, concrete API review and CSV column accounting passed.
Native compilation/execution remains unverified because the constellation
sibling manifests are unavailable in this checkout; no numerical run or
performance result is asserted.

### Joint tangential translation and spin of finite planar liquid walls

`free_surface_wall::planar::FreeSurfacePlanarWall` extends the existing guided
wall and rotor family with simultaneous finite-mass translation along both
plane tangents and principal normal-axis spin. The original material wall
reflection supplies its full three-mode affine impedance, including both
translation/spin cross terms. For the midpoint change q_mid-q_0, one
mass-whitened three-by-three Cholesky solve supplies the shared twist; final
native replay supplies the actual wrench to the original MBD integrator.
Neither the solved twist nor a fitted fluid-added mass replaces the actual
native position, attitude and momentum endpoint.

Both fluid candidates use the same free half-drift COM as their moment origin.
This retains full orbital as well as spatial spin momentum under the native
midpoint position update. Actual linear momentum, spin, orbital-plus-spin,
midpoint velocity/angular velocity, affine-response replay and mechanical work
must pass separate explicit bounds before the liquid, body or clock escapes.
Kinetic and source work use impulse/midpoint pairings; no residual is assigned
to heat or repaired by an energy rescaling.

The guide carries normal force and transverse moment. Raw numerical-pressure,
pressure-perturbation and explicitly absolute physical reactions remain
separate. Their measured work is retained. The complete-domain work allowance
covers exactly two original fluid candidates: the discarded pilot and final
replay. There is no coordinatewise root iteration or finite-difference sweep.

Focused regression definitions cover independent off-centre reflection
impedance on all six planes and three mass/inertia ratios, native replay and
driven continuation with noncollinear COM motion, full angular-impulse arms,
pressure-reference covariance, and source/geometry/budget/cancellation/late
endpoint refusal with unchanged-state retry. Fixed flooded geometry excludes
normal motion, tilt, finite plate edges, cover/vacate, moving contact lines and
caloric heat. The implementation does not add a whole-liquid energy or complete
coupled temporal-order claim.


## Original rigid matter released into thermal free-surface motion

`liquid_transfer::rigid::release_rigid_liquid` removes a no-storage rigid
kinematic constraint from an already-liquid tetrahedral material body. The
existing `TetrahedralDensityBody`, its actual `PrincipalFrame`, accepted
`RigidBodyState`, nodal specific enthalpy and physical clock supply the
source. The full original inertia and COM must match exactly. Native point
kinematics derives every world vertex and `P/m + omega cross r` velocity;
nominal inertia, caller-substituted velocity fields and zero-spin fallbacks
are not used.

Every original nodal density retains the receiving liquid's immutable card,
selected density claim, nonthermal context, quantity conventions and constant
reference density. Equal total mass alone cannot admit a different density
distribution. The original `fs-lbm` tetrahedral overlap and thermal transfer
owners then preserve per-piece mass, both momenta, caloric state and ingress
lineage, including their P1 viscous-strain initialization and optional native
rheology/Boussinesq attachments. Their complete source-domain, positivity,
Mach, geometric-resolution, allocation and projection gates remain binding.

Independent checks compare integrated P1 and actual native receiver M/P/L/h
against the original rigid matter. Original cell masses supply a compensated
source enthalpy integral. The report distinguishes rigid/P1 identity error
from signed receiver kinetic and gravitational projection changes. Neither
mechanical discrepancy is deposited as heat; the original enthalpy datum is
retained without an invented internal-energy or pressure-work conversion.
The receiver starts at the original accepted time and can continue through
its ordinary thermal VOF step.

Preparation is linear in original vertices and tetrahedra before the existing
overlap work. It adds no eigensolve, constitutive law, fluid solver, trajectory
clone or physical timestep. Refusal and cancellation return no receiver and
leave the original source available for retry. Successful callers replace the
rigid representation rather than evolving both representations of the matter.

This entry requires material already fully liquid and a rigid model with no
stored elastic, plastic, damage or other mechanical memory. It does not melt
a stressed solid, detach stored cohesive bonds, supply density-change volume
work, model melting kinetics, or introduce capillary energy or a pressure
projection. Four focused G1/G3/G4 regression definitions cover independent
rotated-cuboid inertia, world and per-piece momenta, enthalpy and clock
continuation, frame/material/phase refusal, density-distribution refusal,
budgets, cancellation and retry. Native compilation and test execution remain
outstanding in the delivery environment; no measured performance or physical
validation claim is added.


## Simultaneous hard impacts on continuing modal and cavity surfaces

`render::schedule::force::coupled::impulse::impact::multiple` prepares complete
relative attachment rows and their full mass-normalized Delassus matrix.
`PreparedModalImpacts` retains shared-mode and cross-contact terms; one original
`fs-contact::multicontact` solve supplies the simultaneous normal impulses, and
one compensated transpose application updates the original modal momenta. The
solver's redundant-row velocity determinacy remains distinct from individual
contact-load authority. Restitution choices that create global energy refuse.

`CoupledModalSystem::apply_normal_impacts` independently checks actual native
momentum increments, all active contact velocities and complementarity, net
contact work and kinetic-energy change before one atomic publication. Failed
geometry, cancellation, limits or balance checks retain every accepted state.
Positions, clock, held inputs, preceding frames, material relaxation and proper
radiation histories remain continuous. Prepared rows and projection scratch are
reused; the bounded native contact solver and returned receipts still allocate.
No allocation-free multicontact or benchmark speedup is claimed.

`SharedCavity::apply_surface_impacts` admits only original physical surfaces and
checks actual volume flows before publication. Finite strikers may have zero
cavity-volume shape. Neck and auxiliary acoustic velocities, compression,
viscous/thermal/molecular memory and radiation history remain continuous and
respond during subsequent original evolution. These APIs consume an explicitly
supplied instantaneous contact set: they do not discover collisions, solve
friction or sustained contact, infer a thermal loss partition, or establish
spatial momentum conservation absent the represented mechanical symmetry.

`tests/modal_impacts.rs` covers analytic finite-mass chains, contact permutation,
redundant rows, global energy refusal, inactive cross-response, atomic retry,
small momentum beside large kinetic backgrounds and continued live-memory
motion. Its `modal_impacts/cavity.rs` child retains the full shared-cavity
integration fixture separately from the core modal tests.

## Finite planar body momentum on the original rigid ALE fluid

With `experimental-stokes-motion`, `stokes_motion::rigid::PlanarRigidMotion`
closes both world translations AND principal-z angular momentum through the
original nonlinear rigid-domain ALE owner. Geometry and moment arms follow the
body's finite motion. The body z axis must be principal and world aligned;
transverse momentum is refused rather than silently projected. The whole
reference domain translates/rotates with the body; independently moving
obstacles, deformation and remeshing remain outside this owner.

The original native momentum impulse constructs each actual trial rate. A
read-only native ALE kinematics preview supplies the SAME represented endpoint
pivot, orientation and mesh rates to wall callbacks and to the subsequent flow
candidate. Initial strong normal means and linear moments include rigid spin.
The new pose, unwrapped angle and each actual native rate must agree with the
moving domain before acceptance. Rotation increments stay below pi; accumulated
angle can cross quaternion/principal-angle branches without resetting springs.

One bounded coupled endpoint search measures three native flow response columns
and uses the original reusable small LU solver. Its algorithm is shared with
the existing two-coordinate translation owner; all base/probe/backtracking
solves consume the declared budget. Dimensional force/torque gates and physical
probe scales remain distinct. Actual represented cross-coordinate rate changes
enter the correction; fixed-base residual scaling controls descent. This is a
local method with native nonlinear refusal, not a global convergence theorem.

Original physical stress alone supplies force and torque about the moving COM.
No-slip and finite Navier-slip body walls are supported. Advected momentum and
physical slip loss retain separate native accounts. Body kinetic/spring change,
physical mechanical damping, signed external/interface work and backward-Euler
numerical loss are admitted separately. Slip loss is reported in joules and
is not added again to the body-only work identity or assigned to heat.

The unique flow candidate commits once, together with body pose/momentum and
unwrapped angle after final cancellation admission. Refusal preserves body,
fluid history, geometry and clock. The body work gate is not an exact finite-step
whole-fluid angular-conservation or passivity theorem: native covariant BDF1
geometry/time defects, fluid dissipation, streams and sources retain their own
accounting and require refinement.

Four focused G1/G3/G4 regressions cover an independently manufactured rotating
Couette field with nonzero force and torque, springs/damping and no-slip/finite/
perfect-slip walls; common boosts; the independent closed-fluid added-mass
limit and existing translation parity; continued finite rotation across angle
branches; actual normal-moment admission; and budget/datum/final-cancellation
rollback with exact retry. Pinned Rust parsing is recorded separately from
native compilation and numerical execution, which remain unrun in this session.

## Simultaneous Coulomb impacts on continuing modal and cavity surfaces

`impulse::impact::frictional::PreparedModalFrictionalImpacts` retains all three
reciprocal contact-frame rows and their complete mass-normalized Gram mobility.
One original `fs-contact::frictional_multicontact` solve accounts for normal/shear
and cross-contact response; one compensated transpose application transfers its
impulses to the original modal momenta. Normal, sticking and sliding reactions
are not split into separate corrections. Locked tangent rows remain exactly
locked, and the original solver's velocity-only determinacy does not authorize
a unique local load or heat allocation.

`CoupledModalSystem::apply_frictional_impacts` checks actual native momentum
increments, actual normal and tangential velocities, complementarity, the
circular Coulomb cone, final sliding direction, contact work, net passivity and
kinetic-energy closure before one atomic publication. Failure or cancellation
retains the accepted network. Normal and frictional impact owners share the
same original relative-row construction and native momentum admission.
Positions, clock, held forcing, previous frames, material relaxation, acoustic
memory and radiation state continue across the instantaneous momentum jump.

`SharedCavity::apply_surface_frictional_impacts` restricts endpoints to original
physical surfaces and admits the candidate's actual cavity volume flows before
publication. It does not impart a direct impulse to neck or auxiliary acoustic
coordinates; those evolve through their original continuing dynamics.

Preparation reuses contact rows, full mobility and projection scratch. The
bounded native solve and result still allocate. Frames are explicit frozen
orthonormal physical inputs; their projected modal rows need not be orthogonal.
The API admits at most 21 contacts and retains the existing setup-work budget.
Collision discovery, moving-frame transport, persistent friction histories,
material inference and a per-contact thermal partition remain outside this
instantaneous owner. Spatial momentum conservation requires the corresponding
symmetry in the retained modal basis. No global nonsmooth existence theorem or
measured performance improvement is asserted.

`tests/modal_frictional_impacts.rs` checks analytic finite-mass stick/slip under
tangent rotation, coupled off-centre spin and orbital angular momentum, shared
body cross-contact mobility and stable-ID reduction, locked tangents and native
determinacy, invalid frames, global energy refusal, cancellation, penetration,
late native-state refusal, and continued material/acoustic/forced histories.


## Partial release of original rigid liquid with a moving retained body

`liquid_transfer::rigid::partial::release_rigid_liquid_cells` partitions selected original tetrahedra through the material geometry owner and returns both an ordinary thermal VOF liquid and a source-bound rigid remainder. Original cell masses, density receipts, full tensors, topology, nodal enthalpy, and physical clock remain attached through explicit selected/retained index maps. Child velocities come from the original native rigid velocity field sampled at each actual center of mass; the existing principal-frame and full-inertia momentum owners construct the child states without release impulse or supplied work.

The complete original source must already be liquid on the same constant reference-density caloric chart and have no stored elastic, plastic, cohesive, or history energy. Every retained tetrahedron must have a world bounding box strictly outside the receiver grid plus explicit clearance. Touching or ambiguous overlap refuses, including overlap with dry or stationary-solid voxels. This is a handoff-time exclusion and does not authorize later grid entry without an admitted moving occupancy/coupling owner.

Native child M/P/L/h and zero-work kinetic identity are checked against the original before transfer. The actual receiving liquid plus retained native body are checked again against the original with one global budget. Numerical kinetic and gravitational projection changes never become heat. Whole-source caloric admission and selected-piece transfer share one source-query ceiling. Any source, topology, phase, occupancy, projection, budget, or cancellation refusal returns neither child and leaves the borrowed source unchanged.

Three focused regressions cover analytic rotated unequal-body accounting, actual liquid readback, native continuation, successive releases and composed lineage, cold retained matter, missing occupancy, strict projection budgets, cancellation, and deterministic retry. Source parsing, formatting, whitespace, and API checks passed; native compilation and tests were not executed in this resource-constrained session.

## Spatial rigid gravity through mixed finite-solid contact trajectories

`FiniteSolidBody::rigid_free_flight_candidate_in_field` and
`mixed_sphere_trajectory_candidate_in_field` connect the original
`RigidGravityField` to the existing free-flight and transient mixed-sphere
owners. The caller supplies an explicit field and frame and may retain an
explicit uniform pressure on the solid's original closed boundary. A solid's
separately bound stationary field is neither inferred from nor replaced by
the rigid field.

Each rigid step reuses native endpoint-load Verlet and retains both actual
field samples and the state after free drift. Central forces and explicitly
selected quadrupole force/torque enter the original full-inertia rigid owner;
there is no extra field query, alternate integrator or inferred finite-body
gravity model. The original uniform-gravity entry and the uniform helper used
by inelastic trajectories retain their previous map and accounting semantics.

Rigid potential decrease is physical source work once. Held world-force
displacement work and body-torque kick work stay separate from free-rotation
drift error, gravity half-kick work minus potential decrease, and the difference
between held-force kick work and displacement work under changing field loads.
Signed defects close the measured work identity; the sum of their individual
absolute magnitudes must independently pass the integration allowance. Opposite
errors cannot cancel through admission, and none becomes contact heat.

The same field reaches every original geometry retry and accepted free phase.
Only accepted phases accumulate time, source work, impulses or numerical
defects. Actual endpoint forces and rotated body torques supply the external
linear/angular impulse account; original simultaneous rigid/solid and
rigid/rigid contact retain their own unallocated mechanical loss. Final
whole-trajectory work, momentum and absolute-defect gates remain binding.
Native field/domain/frame refusal, exhausted work or cancellation returns
neither accepted mechanical candidate.

This remains the original transient contact owner with its declared affine
numerical collision path. Endpoint gravity samples do not enclose the field
domain between them; close approaches require caller refinement and admission.
The sustained first-order contact lane and inelastic material trajectories
retain their existing uniform rigid-gravity interfaces. No source recoil,
gravity inference from an undeclared shape, continuum path certificate or
global impact-accuracy claim is added.


## Finite planar wall recoil with thermal liquid

`free_surface_wall::planar::thermal::FreeSurfaceThermalPlanarWall` retains the original thermal VOF owner and native full-inertia rigid body. The isothermal and thermal variants share one planar mechanical implementation. One original rigid-link pilot and one final native replay feed the existing mass-whitened three-mode midpoint solve, preserving translation/spin coupling. Actual replay loads drive the native MBD endpoint, followed by the original linear, spin, orbital-plus-spin, midpoint, affine-response and work audits before either owner escapes.

Discarded trials advance no accepted source history, enthalpy or clock. The final thermal receipt distinguishes authored and optional molecular heat; mechanical wall, applied and guide work remain separate. Existing body, full-grid trial, motion and work limits remain binding. Geometry is a fixed flooded plane with tangential translation and principal spin; cover/vacate, tilt, finite edges and a solved conduction boundary are not supplied. Four regression groups cover all-plane isothermal parity, heated/buoyant finite recoil, molecular heat applied once, and transactional late refusal/cancellation/retry. Source parsing, formatting and scoped whitespace checks passed; native build and test execution remain outstanding in this resource-constrained session.

## Native rigid bodies striking continuing modal and cavity surfaces

`impulse::impact::rigid::PreparedRigidModalImpacts` prepares complete world-vector
modal material-point traces and their double-double Gram response. Each event
recomputes the original rigid bodies' current-pose translation and full-inertia
response. Shared bodies, shared modes, off-centre spin and normal/tangent terms
enter one original simultaneous Coulomb solve; the same rows transfer reciprocal
impulses through native rigid momentum application and continuing modal staging.

Actual per-coordinate momentum, endpoint contact laws, finite-witness impulse
couples, signed contact work, net passivity, kinetic closure and original modal
state budgets all precede publication of either receiver. Tiny impulses cannot
disappear behind a large background momentum. Stable contact IDs order the
projection and work reductions. Redundant rows retain the original solver's
velocity determinacy and unresolved local load/heat allocation.

`SharedCavity::apply_surface_rigid_frictional_impacts` admits only physical
surface components and validates actual cavity volume flows before either
rigid or modal publication. Failed late admission and cancellation retain both
accepted receivers and allow prepared scratch to be reused. Rigid poses, modal
displacements, clocks, previous interval records, held inputs, material
relaxation, radiation and auxiliary acoustic histories remain continuous.

This is an instantaneous event on explicitly supplied material witnesses and
frozen physical frames, with linear modal point geometry. Conservation of
spatial momentum requires the corresponding represented modal symmetry. The
API adds no collision discovery, persistent-contact integrator, shock model or
thermal partition. Prepared modal rows and scratch are reused, while the native
contact solver and returned receipts still allocate; no measured speedup is
claimed.

G1/G3/G4 regressions in `tests/modal_rigid_impacts.rs` cover analytic finite-mass
stick/slip, tangent rotation, full-tensor off-centre spin, a distant world origin,
shared-body/mode simultaneous response, native redundancy, current rigid
pose/inertia, small losses beside common kinetic motion, refused source and
geometry cases, lost rigid momentum, and continued material/radiation evolution.
The owner-local test exercises late receiving-state refusal and cancellation
immediately before joint publication.

## Material rheology in finite planar fluid-body momentum

With `experimental-stokes-motion`,
`stokes_motion::rigid::PlanarRigidMotion::step_with_rheology` closes the same
three actual body momentum equations with shear-dependent constitutive stress
inside each original rigid ALE nonlinear solve. Newtonian and rheological
entries share one mechanical loop, original native impulse construction,
represented wall/mesh kinematics, response search, stress force/torque, work
gates and joint candidate commit. A statically dispatched trial policy selects
only the original native constitutive solve; no duplicate body solver, boxed
per-cell law or frozen initial-viscosity surrogate is introduced.

The resolver receives the endpoint WORLD centroid, invariant physical strain
rate and immutable reference transport snapshot. Each accepted response must
retain that source density. All native constitutive, nonlinear, boundary and
source refusals propagate. Temperature/composition can be captured in a spatial
resolver, but thermal state and composition are held during this mechanical
interval; mechanical pressure is never substituted for an EOS query. The
callback must be repeatable for identical input. Its external side effects
remain outside transactional state ownership.

`PlanarRigidStep<G = ()>` retains the original material owners in its generic
`RigidAleStep<G>`; `flow.rheology_samples()` exposes only the accepted endpoint
responses in cell order without querying the source again. The ordinary
Newtonian result keeps the default type and an empty response collection. The
same accepted body and fluid history may continue with either constitutive
entry. No-slip and finite Navier slip retain physical stress and dissipation;
a Newtonian slip-length formula is not used as a nonlinear constitutive law.

Three G1/G3/G4 regression groups independently manufacture changing rigid
translation/spin and changing Carreau shear, including thinning/thickening,
finite slip, shared boosts, force/torque/work and world-position source owners;
check the closed-fluid added-mass/zero-strain Newtonian limit and continuation;
and refuse source failure, changed density and final cancellation with exact
retry. Pinned nightly Rust parsing is separate from native compilation and
numerical execution, which remain unrun in this session.

Scope remains constant-density generalized-Newtonian flow on one globally rigid
two-dimensional domain with a principal-z body. Viscoelastic memory, density
evolution, thermal deposition, independent obstacles, exact finite-step whole-
fluid angular conservation and a full coupled passivity theorem are not added.

## Finite-inertia planar induction rotor

`planar_induction::PlanarInductionRotor` couples the original linear
`fs-em::PlanarEddy` moving-conductor midpoint to the original
`fs-mbd::articulated` prescribed-base inverse-inertia/impulse owner. The
mechanical model is one original revolute-z root link on a stationary support,
with its original centered axisymmetric spatial inertia. Field history,
unwrapped angle and angular velocity form one immutable accepted state.
The API returns both field and mechanical histories together only after
coupling, field, work, impulse and cancellation admission.

The admitted physical conductor is a homogeneous isotropic disk or annulus
whose analytic occupancy is invariant under rotation. Binding checks all
original conductor cells: one source-resolved material and isolated current
region, positive triangle orientation, opposite interior edges, one connected
component and complete outer/inner boundary cycles winding once about the
native axis. Every boundary chord must fit the declared radial approximation
budget. A finite square, offset rotor, sector, extra hole or varying conductor
source cannot silently become an invariant rotor. The fixed polygonal mesh
still approximates circles: this admission does not certify exact polygon
no-penetration, geometry convergence or electromagnetic field accuracy.
No interface, material source or mesh actually rotates.

The complete field domain must have uniform scalar permeability. This avoids
claiming full bearing-force closure from Lorentz volume forces while omitting
permeability-interface traction. Every conducting cell belongs to this rotor;
separate conducting stators require another composition. Impressed current
inside the rotor refuses at binding and every endpoint, because the extracted
loads represent the native conduction current. Fixed external impressed-current
coils in insulating cells retain their original source-work owner.

Each bounded scalar trial uses the original accepted magnetic history and
rigid velocity `omega * [-dy, dx]`. The native nonsymmetric induction solve
supplies the actual midpoint current, field, Joule heat and Lorentz work.
The original exact P1 load integrator supplies the full force/torque wrench,
which is transformed to the actual midpoint link frame. The native ABA
impulse response then applies its integrated body load and any separately
declared held shaft torque. Its original constrained mass, support reaction,
representability and joint-limit checks remain active. The impulse API is
used as an inverse-inertia map inside this finite-duration composition;
the composition advances the rotor angle itself.

The converged field speed must agree with the mean of the actual old/new
mechanical speeds. Actual kinetic increment minus native field mechanical
work and shaft work is independently bounded, as is the combined magnetic
plus kinetic storage change plus Joule heat minus impressed-source and shaft
work. Native support impulses, external impulses and their angular/linear
residuals remain separately available. A stationary ideal support has zero
work. Numerical field, mechanical and coupling defects are exposed; none is
reclassified as heat, removed from torque or repaired by changing state.
Nonconvergence has no prescribed-speed fallback and returns no partial state.

The implementation reuses the prepared original magnetic operator, native
Krylov/storage/Peclet controls, original P1 force transfer and original linear
storage ABA owner. The immutable native model, conductor selection, unit rotational velocities
and geometry receipt share one private reference-counted binding, so accepting
a new history does not copy mesh-sized prepared arrays. Bounded field trials
reuse their velocity buffer. No performance, nonlinear-material,
moving-mesh, full-machine, thermal-feedback or cross-platform claim is made.
Focused native tests exercise genuine free magnetic DOFs, signed braking,
finite-inertia feedback versus frozen-speed forcing, continuation, driven
annular rotation, source/geometry refusals and final-publication cancellation.



## Native articulated mechanisms against continuing modal surfaces

`PreparedArticulatedModalImpacts` and
`CoupledModalSystem::apply_articulated_frictional_impacts` connect original
free/prescribed articulated trees to complete modal world-point traces through
one simultaneous 3N Coulomb solve. Original ABA mobility retains root joints,
cross-link coupling, full inertias and COM offsets. Mobility probes do not
apply fictitious unit events or enforce their speeds; actual combined native
responses retain physical endpoint joint limits and prescribed reactions.

Rigid and articulated consumers share the same modal frame/trace/Gram
preparation. Native geometry and inverse-mass responses are rebuilt at every
event. Actual DD link-COM kinetic/spatial increments, per-coordinate velocity
changes, modal momentum, endpoint contact conditions, finite-witness couples
and contact/support work are checked before either receiver publishes.
Prescribed support work is an explicit source, distinct from passive contact
work. Redundant contact/support allocations retain native determinacy and
representative-work qualifications.

`SharedCavity::apply_surface_articulated_frictional_impacts` permits original
physical surface incidence and checks actual volume flow before publishing
modal or mechanism state. Configuration, clocks, prior interval receipts and
material/radiation histories remain continuous. The caller selects states at
one physical event; independent clocks are not resynchronized. This adds no
collision discovery, force duration, contact compliance or thermal allocation.
Modal spatial conservation extends only to represented basis symmetries.

Six new native integration groups cover a rotated free-chain mass-matrix
oracle, driven offset-COM hinge/source work, shared-modal free/prescribed
recoil and permutation, mobility probes versus actual joint limits, redundant
determinacy and atomic modal-budget refusal, and live Maxwell/radiation
continuation. All six passed, as did twelve existing rigid/modal and continuing
impulse groups. Two owner-local tests passed for late-admission/final-cancellation
atomicity in both bridges. Validation used bounded original-source Rust graphs;
the full cavity wrapper received source review but its complete fixture and
the full workspace/DSR lane were not executed.

## Ideal supports in finite planar fluid-body momentum

With `experimental-stokes-motion`,
`stokes_motion::rigid::PlanarRigidMotion::set_constraints` selects stationary
world x/y COM coordinates and/or principal-z rotation on the original native
fluid-body owner. `PlanarRigidConstraints::SHAFT` fixes the COM while allowing
spin; `FIXED` fixes all three coordinates; public flags also admit guided
translation. Lock activation requires an exactly zero actual native rate.
No momentum is projected away and no stopping impulse or lost energy is
inferred. Unlocking changes no accepted body, fluid, geometry or clock.

Newtonian and rheological entries retain one original mechanical loop and
constitutive owner. Unlocked coordinates solve the same finite momentum
equations. Locked coordinates retain exactly zero native rate and represented
drift, and return the missing physical force/torque as `PlanarRigidSupport`.
Reactions act ON the body, angular quantities are about its COM, and finite
reaction impulses are checked independently. The stationary support performs
zero work; its opposite load belongs to the external support, whose dynamics
are not evolved. Existing fluid stress, slip loss, mechanical work and joint
cancellation/publication gates remain binding.

The shared response search probes only free coordinates and embeds them in
the original small reusable LU with inactive identity rows/columns. No penalty
stiffness, extra body solver or additional per-cell constitutive dispatch is
introduced. A fully fixed body performs one native fluid solve and needs no
response probe or LU solve. Locked coordinates may use zero-width bounds at
zero and zero probe steps; the ordinary free path retains its prior controls.

Four focused regressions cover x/y-guided closed-fluid added mass and reduced
flow-solve count, a shear-thinning/thickening shaft with finite slip and actual
support force, a fixed body carrying nonzero fluid stress in one solve, and
exact activation/refusal, overflowing support impulse, final cancellation and
deterministic retry. Five changed Rust files passed pinned nightly parsing.
Native compilation and numerical execution remain unrun because storage is
full and remote build tools are unavailable.

These are ideal stationary world-aligned supports and a shaft through the COM.
Off-COM hinges, moving supports, general Pfaffian constraints, independent
obstacle motion, thermal deposition and whole-system passivity remain outside
this owner. The underlying constant-density rigid-domain and finite-step
angular-conservation limits are unchanged.

## Physical modal damping rates in coupled consumers

The modal bank's `damping_rates_per_s()` is the physical coefficient in
`q'' + c q' + omega^2 q = F`, including free coordinates with `omega = 0`.
Coupled exact and prepared propagators, harmonic response, frequency tangents and pole pencils,
continuous broadband observation, inverse contact responses and acceleration
readouts consume that same coefficient. A damping ratio is diagnostic metadata
and is not sufficient to reconstruct resistance at zero or very small natural
frequency. Independent ratio-authored source models retain their original law.
The linear impact and cavity exchange compilers preserve their supplied physical
rates directly. Intrinsic body, acoustic and neck drag consumes no synthetic
bilateral connection slots; actual volume springs and spatial dampers retain
their existing connection budgets and complete signed columns.

Component heat remains `integral(c v^2 dt)`; harmonic peak phasors dissipate
`c |v|^2/2`. Original signed connection columns and separate material, fluid
and connection losses retain their ownership. Material heat caches bind the
actual rates. A source viscosity update replaces only its retained resistance
rows, preserves other material/ambient drag, and changes no accepted state,
stored energy, mode ordering or clock. Changed rates invalidate old heat and
forcing preparations even when a free mode's diagnostic ratio is unchanged.

`CoupledModalSystem::set_component_damping_rates` also exposes those continuing
intrinsic rates to callers of the generic coupled owner. A sorted, unique list
of `damping::ComponentDampingRate` selects original `(component, mode)` pairs
and replaces their complete finite nonnegative rates, including drag on free
coordinates. It uses the existing material-mobility transaction and original
ordinary/prepared/harmonic-only construction budgets. Mechanical coordinates,
Maxwell memory, radiation state, connection coefficients, sample clock and
previous interval reports survive unchanged; the instantaneous change supplies
no impulse, work or heat. Empty and unchanged batches rebuild nothing and
allocate nothing. Invalid input, cancellation and failed operator preparation
leave the accepted owner intact. The caller owns the complete intrinsic rate
and its physical/temporal applicability; selective source-temperature updates
continue to use the source-bound thermal owner.

`tests/coupled_damping_update.rs` exercises underdamped, critical and overdamped
continuation against closed-form motion, independently integrated component
and connection heat, warm Maxwell/radiation histories against the continuous
physical ODE (ordinary) and the existing average-velocity exchange (prepared),
selective free-coordinate updates across components, harmonic-only response,
and atomic refusal. These are numerical tests of the supplied law, not measured
material-source validation.

Focused regressions cover a damped free coordinate's separate component/link
heat and harmonic force response, source-selective viscosity refresh, and
held/ramp/harmonic broadband convolutions. These extend the existing finite
modal approximation; they introduce no damping floor, inferred drag or new
constitutive source.

## Reciprocal finite-solid / point-source gravitational kicks

`FiniteSolidBody::reciprocal_gravity_kick_candidate` applies actual Newtonian
tidal impulses to the original deforming solid and equal/opposite recoil to a
finite native rigid source. `CentralGravity` owns the field; hyper3 owns the
current-position quadrature with the same original cellwise or nodal reference
density as inertia. Its force/potential-only entry shares the full gravity
evaluation and skips sparse Newton-tangent assembly. The original consistent
mass and rigid momentum owners perform both velocity updates. The source is
a monopole at its COM: its spin is unchanged, without an invented compensating
couple for quadrature error.

Native per-owner residuals and a complete pair check bound actual linear and
orbital-plus-spin impulse and impulse work. Solid support reactions remain
external exchanges, including when the source recoils against an anchored
solid. The balance origin is explicit; no large background momentum enlarges
the absolute tolerances. Failure or cancellation publishes neither candidate.
The frozen pair potential is counted once. The kick preserves configuration,
material history, acceleration and physical time: it is a force subflow for
caller-owned splitting, not a complete moving-source trajectory or a claim of
energy conservation during this isolated kinetic change. Source multipoles,
solid self-gravity, near-source quadrature accuracy and timestep convergence
require their respective physical models or refinement. The same interaction
must not also be supplied through the stationary external-gravity binding.

### Moving reciprocal solid/source trajectories

`reciprocal_gravity_step_candidate` composes the original reciprocal half-kick,
native objective solid energy/momentum step plus symmetric rigid free drift,
and a second half-kick evaluated on both moved configurations. The native solid
clock/history advances once. Its retained acceleration belongs to the material
middle flow; gravity remains an explicit split impulse. Held solid loads,
stationary external fields, optional uniform pressure, source free spin and all
stationary support reactions retain their native owners and work conventions.

The gravitational splitting defect is kick work plus the change of the single
pair potential. It has its own absolute gate. A separate complete energy gate
uses the sum of absolute middle, gravity-splitting and impulse/work defects, so
opposing numerical errors cannot conceal each other. Full pair linear and
angular increments use native measured residual identities, without another
mass action or subtraction of large endpoint momenta. Individual native gates
and noncancelling stage residuals remain active; support exchanges are external.

Explicit duration/node limits bound the wrapper's interval and extra storage.
Per-element swept source-relative boxes enforce a declared minimum separation
before and after the middle flow. Outward-rounded coordinate gaps give a
conservative lower bound along linear interpolation of the accepted endpoint
geometry and source COM; a source cannot tunnel through that numerical path.
This is not continuous physical-path CCD, near-field quadrature certification,
finite-source contact or a global error bound. Smooth trajectory accuracy still
requires spatial/time refinement. Every late separation, energy, momentum or
cancellation refusal publishes neither state.

## Spatial molecular absorption and modal wall diffusion on the exact clock

`ContinuousModalLine::new_field_molecular` and `new_field_molecular_loaded`
retain the original `AxialGasField`, its local caloric reservoirs and both
actual endpoint gas references. `new_viscothermal` and
`new_viscothermal_loaded` retain the complete uniform-cylinder modal field
and the original positive wall-diffusion histories. Signed modal wall traces
preserve coherent pressure and tangential-velocity interference. All acoustic,
thermal, molecular, wall and attached mouth-load coordinates participate in
the existing prepared exact scattering flow and its complete resource limits.
No retained mode or positive reservoir is removed to satisfy capacity.

`duct_state()` includes every physical duct-memory coordinate, and
`load_state()` begins at `source.state_dimension()`. Original mouth-load
binding, face work, independent physical heat, cancellation and exact clock
subdivision remain binding. Empty molecular laws retain the existing numerical
operator and time samples bitwise. The gas background is prescribed and
frozen. Wall physics requires the original cylindrical shape, uniform gas and
declared thin-boundary-layer frequency band; generic time evolution does not
establish accuracy outside that band. Molecular and wall attachments refuse
each other until a shared thermal realization is supplied.

Focused native checks exercise the original spatial molecular law, wall
trace/lift/loss kernels and source-row resistance reconstruction. The authored
loaded/unloaded exact-runtime regressions remain unexecuted in this workspace
because the full constellation dependencies are unavailable.

## Finite-band thermo-optical feedback

`SpectralThermoOpticalSource` extends the same sealed source interface with
explicit finite-band `SpectralOpticalStatePoint` families. Construction admits
the original n/k projection or n/alpha table, a caller-declared maximum knot gap
in vacuum metres, and the same-state magnetic permeability through
`fs-em::plane_wave::SpectralOpticalResponse`. The existing thermo-optical driver
refreshes each source at its actual accepted cell temperature, retains its
original family, claims and gap limit, then readmits the current frequency,
source band and permeability before applying coherent absorbed work and before
publishing the thermal endpoint. Original source wavelength properties therefore
reach the same coherent film and implicit heat owners without a visible-band
clamp, substituted spectrum, nonmagnetic default or second energy debit.

Focused infrared regression sources exercise mixed n/k and original n/alpha
layers with nonuniform native cell temperatures, local caloric/conductive energy
balance, successive accepted feedback steps, missing permeability, retained gap
and frequency admission, and atomic endpoint-source refusal after the native
thermal candidate is otherwise admissible. The existing fixed-geometry,
quasi-steady split and timestep-convergence responsibility remain unchanged.

### Instantaneous momentum on source-backed gas pistons

`cavity_phs::piston::momentum` applies resolved wall impulses to the original
`GasPiston` and `DoubleActingGasPiston`, including their material-bound wrappers.
Only the one wall velocity changes. Both gas volumes, conserved masses,
compositions, calorics, source records and preceding finite-time port accounts
remain intact. No gas pressure work or heat is assigned to zero elapsed time.
Free walls use their actual conserved mass; fixed supports return the opposite
reaction impulse with zero support work. Endpoint ceilings, represented
momentum and signed work are admitted before publication, using impulse-local
scales and double-double polarized work rather than total-energy subtraction.
Callers account the returned event receipt separately from finite-time windows.

The axial rigid-impact path uses the original rigid point kinematics, full
inertia mobility and unpublished native point-impulse candidate, together with
`fs_contact::impact::solve_normal_impact`. A declared world axis and origin
locate the actual piston face from its original volume/area coordinate; the
rigid witness must lie on that axis. Equal and opposite normal impulses act on
the two native mechanical owners, with actual momentum, restitution and
combined kinetic-loss checks before the final shared cancellation/publication
boundary. Off-centre striker rotation is retained. Impact loss is reported by
the contact law without an invented gas or solid heat partition.

This is one selected frictionless axial event, with frozen configuration and
no gas/material query. It supplies no collision search, oblique guide dynamics,
simultaneous-contact selection, shock closure or sustained-contact integrator.
Focused kernel and source-gas regressions cover kicks, braking, fixed-support
reaction, reciprocal recoil, unchanged gas/port histories, cancellation and
continuation through the original finite-time equations.



## Reflected moving-source histories on the original audio clock

`MovingCompactRenderer::with_ideal_reflections` consumes the original borrowed `PreparedReflectedSourceTrajectory` owners and accepted Q/F histories. Every source shares one original `IdealReflectionObserver`; a separately bound microphone or boundary owner refuses. The existing cubic strength interpolant is compiled once per physical source, then queried at each direct/reflected path's own retarded time. Pressures add coherently in the existing compensated sum before decimation or PCM conversion. No source history or trajectory is copied to create images.

`source_count` continues to count physical sources; `source_path_count` reports all propagation paths. The existing `maximum_source_receiver_pairs` budget counts every path/output pair, and retained root-iteration bounds include that multiplicity. Full history coverage and Doppler bandwidth admission apply separately to every path. Missing reflected prehistory or tails refuse the complete observation instead of becoming silence. Existing global sample indices, no-allocation callbacks, cancellation rollback, pressure units and `PressureRenderer` sinks remain the owners.

Focused G1/G3 tests reconstruct the reflected waveform from independent source-image polynomials, verify original physical-source accounting and exact callback partitioning, retain coherent source cancellation, and exercise path-work, history, boundary identity and transactional cancellation refusals. The underlying ideal stationary-wall, quiescent lossless compact-source model limits remain; this adds no source feedback or experimental acoustic qualification.

## Extended gravity source recoil and tidal spin

`reciprocal_quadrupole_gravity_kick_candidate` and
`reciprocal_quadrupole_gravity_step_candidate` extend the original finite-solid
reciprocal interaction with the native `QuadrupoleSourceGravity` field. The
source's original mass, full inertia and actual COM pose determine one
monopole-plus-second-moment potential, acceleration and analytic spatial
derivative. The original cellwise/nodal reference-mass quadrature and consistent
mass remain the solid's owners. The original monopole APIs and isotropic-inertia
reduction retain the original arithmetic.

The rigid source receives the opposite net nodal impulse and the negative
source-relative moment of those same impulses as a COM spin impulse. Affine
shape completeness makes this the orientation derivative of the same discrete
quadrature potential, up to rounding; it is not a momentum projection. Native
rigid impulse work includes spin, and the joint angular account includes actual
source spin, orbital momentum, and the solid's separate stationary support
reaction. Source-relative compensated moment accumulation keeps the physical
update independent of the reporting origin. Each native participant and the
combined exchange retain their original residual gates.

The trajectory uses the same reciprocal half-kicks and native material/free
rigid middle flow. Both endpoint source positions and orientations are sampled,
so tidal spin changes the subsequent field. The original independent splitting
defect, total work/energy and noncancelling momentum gates remain; no complete
trajectory energy-conservation claim is added. Native acceleration/history
continue to describe the material middle flow.

An explicit source COM bounding radius must be compatible with its native
inertia, and radius/distance must remain below the declared ratio in (0,1).
The frozen kick admits quadrature sample positions; the trajectory additionally
requires every per-element swept relative bounding box to clear radius/ratio
and the caller's minimum separation. These conservative linear-path checks do
not certify curved physical CCD, source contact, neglected higher multipoles,
interior/near-field gravity, or spatial quadrature accuracy. A source/domain,
work, native-step or cancellation refusal publishes neither participant.

Focused tests compare field potential, force, spatial derivative and source
orientation torque; exercise exact isotropic reduction; measure actual pair
linear/spin/orbital momentum including supports; replay the original native
subflows; refine both source models; shift the reporting origin; and refuse
inadmissible extents atomically.
### Simultaneous coherent pressure channels and calibrated WAV output

`pcm_wav::observation::harmonic::HarmonicPressureFrameRenderer` binds up to 64
named receiver/readout channels to one original `PreparedHarmonicArrayRealization`
and absolute integer/DD sample clock. Real and proper-complex factor draws stay
shared across every channel and frequency. Directional readouts reuse the
original pressure/velocity projection, signed weights and causal delays once.
Preparation budgets cover the complete selected bank; original retained lines,
including dark channels, must satisfy strict output Nyquist. Callback storage is
bounded, sampling allocates nothing, and cancellation or physical refusal
publishes no partial frame or counter advance.

The original `Pcm16WavStream` and shared pressure frame loop retain complete
interleaved frames, one authored Pa-to-PCM full scale, whole-request admission,
callback-boundary cancellation and poisoned partial-I/O refusal. Frame counts
are distinct from scalar channel samples. Mono/stereo bytes remain compatible;
3–64 channels use WAVEFORMATEXTENSIBLE direct outputs with a zero speaker mask.
The original PCM24/float32 batch writer supports the same ordered channels.
Applications retain channel identity/geometry/calibration metadata beside WAV
bytes. This is settled harmonic observation and offline file export; it adds
no source process, propagation, transient onset, per-channel normalization,
automatic speaker mapping or inference of ensemble statistics from one draw.

### Recorded generator coherence through the original loaded exterior

`driving_point::network::radiation::PreparedRadiatingNetwork::prepare_welch_bin`
binds one measured interior frequency bin to its exact original loaded network
and retained exterior operator. Recorded channel labels **and** SI source
quantities must equal the canonical `NetworkSourceCoordinate` order. Pressure
channels are open-circuit generator pressures before the complete original
source impedance; prescribed flow channels are graph-inward m³/s. Frequency,
order and units are checked even for silent records. DC and Nyquist bins refuse.

`PreparedRadiatingWelchBin` borrows the original `WelchCrossSpectrum`, retaining
its acquisition and window provenance. It maps each segment factor through the
original loaded source-to-port flow bank, changes graph inflow to outward
opening flow exactly once, and sums complex contributions before any quadratic
observation. The original mutual radiation matrix, scalar/coupled generator
impedances and conventional/regional/molecular duct propagation remain bound.
Internal pressure is never used as an exterior pressure proxy.

Cold `RadiatingWelchLimits` admits the complete segment count, logical factor
storage and mapping work before allocation. The cache contains only
`segments * openings` complex peak flow values; its fixed ceiling is 4096
segments and 32 openings, or 2 MiB logical storage. The admitted product count is
`segments * sources * (openings + 2)`, including the shared peak conversion's
two scaling products per source. Checked counts and allocation refusal publish
no partial preparation. Cancellation polls each segment, with at most the
original 32-source/32-opening bounded composition between polls. These are
logical ceilings, not RSS or latency guarantees. Authored coherence retains its
separate rank-32 limit.

Repeated `evaluate_into` calls reuse the original duct measured contraction for
all p/p, U/U and p/U moments and signed generator/interface/dissipation/network
work. Cached outward-flow factors reach the original BEM p/u field owner.
`PreparedRadiatingMicrophones::welch_bin` retains the complete complex Hermitian
microphone bin CSD and actual near-field pressure/velocity directional response.
`PreparedRadiatingFarField::{observe_welch_bin,evaluate_welch_bin}` reuse the
original positive angular quadrature, physical exterior medium and requested
directions. All measured factors enter one convergence decision for the actual
covariance, including weak angular residuals exposed by source cancellation.
Individual source self-powers cannot certify that convergence.

The shared Welch owner applies `sqrt(2*df)` once; each original harmonic owner
applies the peak time-average half once. Outputs are bin-integrated moments and
watts, not per-Hz densities. The original signed exterior surface work is already
included in network external-load power; observing it never adds a second
radiation debit. A quiet microphone does not imply zero source dissipation or
zero radiated power elsewhere. Overlapping recorded segment factors do not
establish independent trials, stationarity, confidence, a time realization or
cross-frequency coherence.

Each runtime owner retains its own explicit work cap; the cold opening mapping
is not charged or repeated during observations. Exact owner-pointer binding
rejects otherwise same-shaped substitute preparations. `evaluate_into`
invalidates all statistics getters before work and after any refusal or
cancellation. Microphone failures return no slice; caller-provided angular
output must be discarded on any error, as in the original angular API.

Validation boundary: a native retained-operator harness executes the production
adapter and original Welch, duct work, BEM field/microphone/angular and
`CancelGate` methods with independently specified analytical duct maps and
manufactured exterior columns. Six passing cases cover 40 mixed-SI factors,
rectangular/Hann normalization, cross-receiver phase, signed heat/load work,
pressure-null reactive velocity, covariance-specific angular convergence,
identity/frequency/units, limits, cancellation and replay. This does **not**
execute graph/material/BEM assembly or establish mesh convergence. Two additional
repository regression tests exercise the original loaded graph and exterior
preparation when the complete constellation is available. Full Cargo/DSR was not
run here because seven required sibling repositories are absent.

### Coherent finite-footprint heat on original spatial thermal cells

`layered_wave_heating::spatial_interference::SpatialInterferenceHeating` binds
one original angular-spectrum exposure ensemble to `LayeredWaveHeating`'s
existing cell-split optical stack. Every native depth cell is prepared from
that same ensemble; the anchor stack must be the exact original object. The
original optical/thermal owner retains material cards, fixed source conditions,
native cell IDs and equal-medium cell faces. No independent geometry or
material-law mapper is introduced.

The declared continuous rectangle must match the original planar thermal area.
Every original cell integrates its complete thickness with its own authored
positive depth quadrature. Actual continuous cross-mode, sum/difference and
reservoir-filter interference is retained by the original electromagnetic heat
owner. Actual weights are not renormalized. This is a through-thickness thermal
column of the selected footprint; it does not spread patch watts over a larger
body. The optical field still observes the original periodic continuation.
Lateral temperature gradients, lateral heat transfer and patch-edge scattering
are not resolved. The native one-temperature cross section requires its own
physical justification; area matching alone does not establish equilibration.

Prepared cell spectra are already watts. The actual native thermal interval
and retained DD source clock determine analytic per-cell joules. Accepted
continuation advances that clock without resetting source phase. A separately
bounded positive sum of local primitive-ripple bounds prevents opposite cell
beats from hiding temporal forcing. Maximum thermal step, source ripple and
native nonlinear/energy tolerances remain separate controls. No carrier or beat
is dropped because it is faster than the thermal time step.

Each source credit enters once through the original
`LayeredThermalStep::with_cell_energy_j`. The existing thermal residual, sparse
solve, material-domain admission, nonlinear conductivity/capacity, contacts,
baths, authored sources and energy acceptance remain authoritative. The native
audit includes other authored heat and rounding in its mean-power conversion.
Negative or unrepresentable irreversible heat refuses rather than being clipped.
All candidate source/thermal data are withheld on cancellation or any refusal;
the accepted input state is immutable. Original source and solve work bounds
remain enforced, and cumulative preparation storage is conservatively bounded.

Optical coefficients and geometry remain frozen. This does not add a lateral
temperature law, thermo-optical refresh, phase-change chart, transient
polarization history or reversible parameter-work transfer. Statistical inputs
drive the native nonlinear thermal problem with expected heat, not an
ensemble-mean temperature. Callers must not separately credit these same optical
loss joules. The existing uniform-body complete-volume guard is unchanged.

Five exact-source native source-integration/ripple/clock tests passed. They run
the production cell-integration helpers with the original native `fs_fft`
spectral engine. The full adapter also passed a Rust type/lifetime check with
signature-only external-owner stubs; that check supplies no physical-owner
validation. The authored public Maxwell/material/thermal integration test was
not run: the full Cargo/DSR workspace remains unavailable with missing sibling
checkouts. These limits do not establish full optical-depth, thermal-time or
material-data convergence.

## Three-axis porous rotation with native gyroscopic recoil

`porous_momentum::spherical::SphericalPorousBody` couples the original D3Q19
fluid to a finite native rigid inertia whose center is fixed and whose three
rotational degrees of freedom remain free. Radial isotropic Darcy/Forchheimer
coefficients inside a fully contained sphere make the porous support invariant
under every attitude. Authored mass, principal inertia, initial orientation and
spin remain native data; anisotropic inertia can precess. The six box walls are
stationary, with no voxel solids. The guide cancels force at the center only,
so all components of the measured distributed fluid torque reach the body.

Each three-coordinate Newton--Krylov trial invokes the original rigid
Euler-top/Cayley integrator and original nonlinear fluid collision, streaming
and endpoint forcing. The residual uses actual endpoint world spin change
against measured fluid torque and the supplied world external couple. The
accepted candidate is independently replayed through both original owners.
Radial resistance and independent fluid-force callbacks are sampled once per
eligible cell per tick, and native trial/cell budgets bound the solve.

For the native Cayley update, Rbar=(R0+R1)/2 maps body torque into spatial
impulse, while Omega=R0*omega_body_mid=R1*omega_body_mid is its work-conjugate
world angular velocity. The skeleton field Omega cross (x-center) uses that
same velocity. A normalized midpoint quaternion or averaged endpoint world
angular velocity is not substituted. Native attitude phase still requires time
refinement; the declared rotation bound limits h*|Omega|.

The receipt retains native mechanics and fluid exchange, independently reduced
fluid momentum, actual body spin change, wall and guide impulses, external work,
body kinetic change, drag work and relative dissipation. Participant momentum
residual norms and constitutive/interface/rigid work defects are summed in
absolute value before admission, so opposing defects cannot mask one another.
The remaining fluid kinetic change is explicitly transport/pressure/collision
work, not physical heat or a complete fluid-energy certificate. Failed trials
leave the original body and populations untouched. This extension does not
remap a general rotating mask, model moving exterior walls, or infer apparent
fluid inertia or evolving porosity.

## Finite axisymmetric porous-rotor momentum

`porous_momentum::rotor::AxisymmetricPorousRotor` couples a native finite
principal-axis rotor to the original D3Q19 Darcy–Forchheimer flow. The COM and
axis remain fixed while native angular momentum and attitude evolve. A finite
cylinder must fit completely inside a nonsolid box with six stationary walls,
including for an oblique axis. Scalar isotropic Darcy and vector-speed
Forchheimer coefficients depend only on cylindrical radius and axial position;
the material support and law therefore remain invariant under rotation.
Dry inertia excludes the fluid already represented by the original populations.

The coefficient and independent fluid-force callbacks are sampled once per
eligible cell for each tick and retained through the solve. One scalar midpoint
angular speed sets the skeleton field `omega * axis x (x-COM)`. Every residual,
numerical derivative probe and final replay executes the actual native fluid
operator from the same accepted populations. One reused native scratch grid
avoids reconstructing populations or recompiling topology. The original
Newton–Krylov owner solves the measured nonlinear torque response; no affine
pilot response or replacement effective inertia is assumed.

The original rigid midpoint integrator receives the axial drag and actuator
impulses and advances the actual attitude. The stationary bearing cancels the
measured skeleton force and transverse couple, independently of the resulting
momentum residual. Fluid/body linear and angular changes are checked separately,
including measured wall exchange and explicit shifts to the reporting origin.
Opposing participant defects cannot cancel their acceptance gates.

Native constitutive loss, fluid drag work, rotor interface work, actuator work
and the bearing's zero-work condition have independent checks. Remaining fluid
kinetic change includes pressure, transport and collision; it is not deposited
as heat or classified as a numerical defect. Native trial/cell budgets charge
all probes and replay. Fluid populations, rigid state and physical clock publish
only as one accepted candidate after cancellation and physical gates.

The retained fluid source quadrature can reverse slip under stiff unresolved
steps. Spatial, relaxation-time and rotor-phase refinement remain necessary.
This is finite spin in invariant porous support, not general six-degree motion,
moving occupancy, nonaxisymmetric blades, evolving porosity or a caloric model.

Four native numerical integration cases cover finite-inertia slowdown and real
attitude, independent cell impulse/torque/dissipation reconstruction, wall and
bearing balance, oblique-axis motion, reporting-origin invariance, exact driven
zero-drag phase, native population replay and atomic refusal. The original
integration target typechecks against these owners' full actual dependencies.

## Finite periodic porous-skeleton recoil

`porous_momentum::PeriodicPorousSkeleton` couples the original fully periodic
D3Q19 Darcy–Forchheimer flow to a native finite dry-skeleton mass per
representative domain. Homogeneous resistance and fixed material axes make
support invariant under translation; voxel solids, wall/open faces and initial
spin refuse. `PorousLatticeUnits3` explicitly maps the existing cell size, tick
and fluid mass at density one into SI. Dry mass excludes that represented fluid;
no porosity factor, displaced-fluid correction or additional inertia is inferred.

The three shared midpoint velocities are solved by native Newton–Krylov with
fallible numerical Jacobian actions. Every residual resets one reusable native
grid candidate to the accepted state and executes the original nonlinear
resistance step. Native trial/cell ceilings include stencils and the final
accepted replay. Original fluid populations, force fields, rigid momentum, COM
position and one physical tick publish together after cancellation and physical
balance gates.

The nonrotating skeleton receives actual translational recoil. Its ideal
rotational guide cancels the distributed native drag moment about the actual
represented midpoint COM, separately from the held external COM drive. Actual
fluid and body momentum changes are checked independently; angular balance
includes that guide and the native periodic-cut transport receipt. Constitutive
dissipation and rigid interface work close without cancelling participant
defects. Remaining fluid kinetic change includes pressure, transport and
collision; it is not automatically heat or numerical error.

This is a homogeneous periodic representative-domain model, not motion of a
finite bounded porous object or evolving occupancy/material axes. The held
midpoint skeleton velocity retains the original endpoint fluid-force quadrature
and its stiff-step reversal; spatial and time refinement remain required. No
heat partition, complete energy-conservation claim or trajectory-error
certificate is added.

Four focused tests use independent two-mass linear relaxation and kinetic loss,
quadratic continuous-relaxation refinement, anisotropic nonuniform native-flow
replay, cross-axis external loading with independent physical guide moments and
periodic angular balance, shifted reporting origins, and atomic trial-budget,
cancellation and late-endpoint refusal. The original integration target also
typechecks with the complete real dependency graph of these numerical owners.

## Reflected pressure at an original moving microphone

`MovingCompactRenderer::with_moving_receiver_reflections` borrows `PreparedReflectedMovingReceiverSource` histories sharing the same original `MovingReflectionObserver`. Every physical strength interpolant is compiled once; each direct/reflected path queries it at its own retarded source time and contributes signed physical pressure to the existing compensated sum. The microphone's trajectory changes arrival times and source-clock bandwidth, never the Eulerian pressure coefficient law.

The existing source/path counts, complete finite-history checks, source/receiver-pair budgets, retained root work, global sample clock and `PressureRenderer` interface remain the owners. Bandwidth admission includes every image's reflected receiver velocity; a direct-path relative-motion bound is insufficient. A foreign microphone/wall owner, incomplete reflected prehistory, a tighter unsupported residual policy or an out-of-window clock refuses before rendering. Callback partitioning remains bit-identical, coherent opposing sources cancel, and cancelled blocks preserve both output and accepted sample count. The underlying stationary ideal-wall, homogeneous lossless compact-source model boundaries remain explicit; no trajectory, radiation feedback or PCM encoder semantics are introduced.

The same pressure sink accepts reflected moving-source/microphone owners
prepared in one uniform flow frame when the mean flow is tangent to every
stationary laboratory wall. No source history or force vector is copied into
an image history. A full waveform regression compares the rendered samples to
independent laboratory source images with convected propagation, mixed
rigid/release parity and a nonzero spatial/time frame origin. Normal flow,
foreign trajectory frames and moving walls remain outside this specialization.

### Per-cell temperature feedback for finite-footprint interference heat

`layered_wave_heating::spatial_interference::feedback::QuasiSteadySpatialInterferenceHeating`
connects accepted original thermal-cell temperatures to the existing continuous
finite-rectangle heat path. It retains the original `LayeredWaveHeating` body,
material-layer cards, cell order, geometry, footprint and authored quadratures.
Each cell independently refreshes its complete selected Debye/Ohmic material
bundle at that cell's native accepted temperature, including permeability. A
combined Debye/Ohmic law refreshes one shared complete snapshot. No layer-average
temperature, replacement source selection or inferred n+k reservoir is used.

One rebuilt original cellwise stack per distinct authored frequency retains the
actual source-backed incident and exit responses. The existing two-face drive
preparer now serves both the uniform and spatial feedback owners; mean, empty,
real and proper factors retain their original order, joint moments, frequency
coherence and source clock. All cell heat comes from one common ensemble across
the complete rebuilt stack. Distinct frequencies are never tolerance-merged.
The original cell-area/volume and work admission, including the uniform-body
complete-volume guard, remain in force.

One explicit interval holds the initially refreshed optical laws and calls the
existing spatial heat step once. Its analytic source joules, DD phase epoch,
local source-ripple bounds, nonlinear thermal solve and native energy acceptance
are unchanged. After that candidate solve, EACH cell must satisfy the declared
temperature-excursion limit. Every original endpoint material source and actual-
frequency stack is admitted before the joint result is returned. Applied and
endpoint source responses remain inspectable in native cell order; endpoint
admission supplies no second heat credit. Cancellation, source/stack refusal,
excursion failure or native solve failure publishes no partial accepted state.

The endpoint excursion is a control on an explicit quasi-steady split, not a
coupled residual or a temporal-convergence certificate. The existing thermal
column still does not resolve lateral gradients or patch-edge scattering.
Exterior source states, geometry and native thermal laws remain held. This adds
no phase-change chart, finite-band-table-to-causal-law inference, transient
polarization history, reversible parameter work or ensemble-mean temperature.
Preparation bounds cover coexisting original field caches and retained cell
sources; dropped applied caches are not charged as coexisting with endpoint
stack admission. Source payloads and allocator overhead remain excluded.

Seven focused native tests passed: the five existing exact source/clock/ripple
tests plus two new per-cell excursion/domain/cancellation tests. The complete
new adapter and actual shared drive helper typechecked together with original
native math/FFT libraries and signature-only optical/material/thermal owner
declarations. That supplies no full physical-owner execution claim. Two public
actual-owner tests are authored for mixed-bank temperature feedback, independent
two-cell thermal balance, source-joule change, clock continuation, and atomic
endpoint optical-domain refusal. They remain unrun while the full Cargo/DSR
workspace is blocked by the missing sibling checkouts.

### Stationary selected-bin observations of the loaded exterior

`driving_point::network::radiation::PreparedRadiatingWelchBand::stationary`
binds an explicitly selected, strictly increasing set of original
`PreparedRadiatingWelchBin` values. It borrows the exact same
`WelchCrossSpectrum` acquisition and original `RadiatingNetwork` across every
bin. Canonical source labels/SI quantities, graph port order, opening source
kinds/order and receiver identities must agree, including for silence.
Conventional propagation retains one identical `LossModel`; molecular
propagation must borrow the **exact same** `MolecularNetwork` reservoir owner.
Mixed propagation kinds or replacement molecular owners refuse. Original
frequency-dependent source/load impedances remain explicitly authored per bin.

This API declares a stationary-bin sum in which cross-frequency products
average to zero. Neither Welch segmentation nor matching array shapes proves
that assumption. The consumer does not synthesize a waveform, construct
cross-frequency coherence, interpolate missing source/transfer data, infer
uniform spacing of selected frequencies, or infer unselected bins to be silent.
`bins()` retains exact discrete indices and their original preparations;
`coverage()` reports selected versus available interior-bin counts. DC and
Nyquist always remain absent, including when every interior bin is selected.
This metadata does not establish measured bandwidth, leakage control or
stationarity.

The original bin adapters have already applied shared `sqrt(2*df)` and the
original harmonic time-average half. Band observations sum those integrated
moments once, with no additional df, peak/RMS, or one-sided multiplier.

`statistics` prepares bounded reusable original per-bin buffers and aggregate
port p/p, U/U, p/U and exterior p/u statistics. Evaluation preserves the original
signed generator, source dissipation, interface, network and external-load work.
`microphones` prepares the **same physical patterns** at each original receiver
for every bin, preserving frequency-dependent near-field pressure/velocity and
complex inter-microphone phase. Its selected-bin Hermitian CSD remains in Pa².
`far_field` requires original angular preparations at all selected frequencies,
with identical requested directions and exterior fluid. Every actual measured
covariance must pass its own original positive quadrature convergence. Original
per-bin convergence reports remain accessible; summing powers does not create
a broadband convergence certificate. Exterior surface work already belongs to
network external-load power; outgoing observations never add another debit.

The band admits at most 4096 selected bins plus explicit totals for original
opening-factor bytes and already performed mapping products. It never copies
those factors. Observation constructors admit complete retained logical bytes
before allocation, including original per-bin buffers and port-label payloads
where applicable; borrowed original operators and allocator overhead are
excluded. Repeated evaluations allocate nothing and repeat no propagation,
boundary solve or source-to-opening composition. Original per-bin budgets stay
independently enforced. The additional whole-band work gate conservatively sums
those per-bin allowances plus the documented reduction count; unused allowances
cannot weaken another bin's gate. All arithmetic counts are checked.

Stationary reductions use componentwise compensated sums. Microphone and
directional observations accumulate while processing each bin, avoiding extra
bins-by-microphones² or bins-by-directions result banks. Cancellation polls each
bin and bounded reduction chunks. Any identity, finite, work, convergence or
cancellation refusal hides **all** band observation getters, including after a
previous success. There is no partially published band result.

Validation boundary: five passing native tests execute the production band
modules and original retained Welch/duct/field/microphone/angular owners with
independently specified analytical passive source/outlet maps and separated
monopole columns. Sparse frequencies 256, 640 and 1152 Hz test full reciprocal
mutual loading, delays, signed source heat, exterior work, physical p/u and
Hermitian microphone moments; analytical sphere integrals independently match
outgoing power and directional moments, including active boundary work.
Rectangular and Hann normalization, exact acquisition/SI/model refusals,
coverage, storage/work limits, cancellation/replay, missing angular bins,
direction mismatch and one-bin convergence failure are covered. A targeted
propagation-identity test covers conventional mismatch, mixed kinds and exact
molecular-owner identity using opaque identity markers; it does not execute a
molecular constitutive solve. The native proof does not execute original
material, graph or BEM assembly. One additional public original-assembly band
regression is authored but unrun. Full Cargo/DSR remains unavailable because
seven constellation sibling repositories are absent.


### Loaded generator coordinates for physical exterior energy control

`PreparedRadiatingNetwork::source_control` binds the original canonical
generator bank to the original BEM source-control owner. Open-circuit pressure
controls remain Pa before the full source impedance; prescribed volume flows
remain graph-inflow m³/s. The retained complete loaded port-flow map supplies
each actual opening as `Q_out = -U_graph`. Full mutual exterior loading,
source impedance, mutual generator impedance, original propagation, receiver
geometry, frequency and gas remain those of the borrowed preparation.
Neither an independent outlet actuator nor an additional radiation-power
debit is introduced.

`PreparedRadiatingSourceControl` retains the exact original owner, canonical
source identities and defaults. Its `observation()` is the composed original
p/u map; `energy_control()` returns the unchanged
`PreparedHarmonicEnergyControl`. Thus original BEM value/Cartesian complex
gradient, Hermitian pullback and matrix-free HVP operate in actual generator
coordinates, including unselected disturbances and dark inputs. Quadrature
retains original receiver indices and m³ weights. The objective is exterior
residual cycle-mean energy in J with the original peak-phasor factor 1/4;
it is not generator heat, an actuator penalty, or a geometry/material
derivative. Source limits/scales use each original coordinate's SI unit.

The narrow `PreparedHarmonicObservation::compose_sources` factory accepts an
explicit new-source-major, old-source-minor complex matrix at one frequency.
It performs DD composition before one final round-down, preserving weak
product cancellation and all four p/u components. Original frequency,
medium, physical receivers and nondispersive-storage capability survive.
The new producer hash binds the original producer, both ordered source
identities, receiver geometry, medium, frequency, capability and every complex
matrix entry. It never silently inherits an identity omitting the source law.
Existing source-kind hash tags are unchanged; new
`SourceId::OpenCircuitPressure` and `PrescribedVolumeFlow` represent upstream
coordinates truthfully and are explicitly refused as raw panel traces.

Cold admission bounds complete composition products, source/transfer counts,
logical input/retained identity bytes and field storage before composition;
fixed caps remain 4096 sources, 64 million products and 256 MiB logical cold
storage. The loaded-network adapter retains its original 32-source/32-port
caps and bounded stack scratch. Silent columns still undergo complete
identity, shape, finite-value and budget checks. No partial composed owner is
published on refusal/cancellation. Repeated physical field, energy-gradient,
pullback and HVP applications use the original allocation-free numerical
owners and repeat no boundary or network solve.

Validation boundary: five passing native tests execute the production BEM
source-composition factory, including complex Hermitian chain rules, a
`(1+d)(1-d)-1 = -d²` product-cancellation residual, transformation identity,
silent-input admission and preservation of the relaxing-exterior storage
refusal. Three more passing native tests execute the production loaded-source
adapter with retained original duct/BEM operators and independently specified
passive mutual outlet loading plus separated-monopole p/u. Mixed Pa/m³/s
controls agree with analytic energy, all four Cartesian finite differences,
analytic HVP, the physical opening sign and the loaded-source Hermitian
pullback; signed network work is unchanged. Budget, opening identity,
cancellation and replay are covered. This native proof does not execute
original graph, material or BEM assembly. Two original-assembly regressions
are also authored (mixed-SI fresh-field derivatives and full mutual generator
impedance/raw-trace refusal) but unrun. Full Cargo/DSR remains unavailable
because seven constellation sibling repositories are absent.

## Changing-cutoff aspiration colour

`turbulence_noise::ColouredGaussian` carries its normalised first-order
Gaussian state through changing flow cutoffs. For a cutoff held over a step,
`b = exp(-2*pi*fc*dt)` and `y_next = b*y + sqrt(1-b*b)*xi`; the innovation
uses `expm1` so a small positive interval does not clear a warm history.
This is the exact Ornstein--Uhlenbeck transition (Gillespie, Phys. Rev. E
54, 2084, 1996, DOI 10.1103/PhysRevE.54.2084). Replacing a cutoff changes the
correlation time without rescaling the retained sample by a new stationary
filter gain. Each call consumes exactly one original keyed Philox normal.

At a fixed cutoff the law and autocorrelation equal the previous normalised
matched-pole filter up to rounding; previous sample bits are not promised.
A prescribed changing cutoff retains unit equilibrium variance, while the
zero-start and closed-gate reset retain finite warmup. A cutoff coupled to
the random flow state need not preserve an unconditional stationary variance.
The aspiration amplitude law, source placement, authored gain and Reynolds
gate are unchanged. This remains optional stochastic colour, without resolved
turbulence, measured aspiration calibration or acceptance-render authority.

Focused regressions cover the fixed-cutoff limit, independent exact
transitions after a warm history, a tiny positive relaxation interval,
changing-cutoff ensemble variance, keyed replay and gate clearing. Four focused
native Rust tests pass using the complete original math and RNG libraries and
the pinned nightly compiler. Full fs-couple and workspace execution remain
unrun.


### Absorbing exit media on the original layered thermal timeline

`layered_wave_heating::substrate::PreparedSubstrateWaveHeating` binds an original
coherent angular-spectrum beam or monochromatic coherence bank to a finite
planar column of its homogeneous exit medium. Native thermal faces define each
depth interval; the original Maxwell owner integrates whole-periodic-cell
absorption analytically over that interval. Positive watts are cached per
original CellId. No depth quadrature, artificial optical interface, field FFT
or per-step material solve is added, and the finite-film 64-layer bound does
not constrain this homogeneous-substrate thermal subdivision.

The original exit card and shared fixed query context must agree with EVERY
conductivity and caloric source. The thermal cross section must equal the
complete periodic optical cell within roundoff. Explicit cell, cell-mode-factor
work and logical storage bounds precede allocation. Preparation is
O(cells * grid points * factors); repeated held-source intervals add O(cells)
work to the original thermal solver. Coherence is contracted by the existing
source owner; mutually uncorrelated factors are never summed as fields.

`advance_to` converts cached watts to joules over the actual native interval
and calls `LayeredThermalStep::with_cell_energy_j` once. Native nonlinear
capacity/conduction, authored heat, contacts, reservoirs, events and acceptance
remain authoritative. Entering/leaving flux and a closure residual are reported
separately; heat is not recovered from their subtraction or rescaled to close
the residual. Loss outside the represented column and in finite films is not
credited to this thermal body. Failure/cancellation leaves the initial timeline
unchanged.

This is a held cycle/ensemble-mean electromagnetic source. Optical coefficients,
geometry and source statistics do not refresh with the changing cell temperature.
Lateral thermal gradients, carrier transients, multifrequency beats, phase change,
electromagnetic storage changes and constitutive parameter work are not supplied.
Expected heat drives one nonlinear temperature trajectory, not an ensemble mean.

Five authored G1/G3/G0/G4 regressions cover independent oblique Fresnel and
exponential cell integrals against native finite-volume temperatures, 96-cell
refinement and partial coherence, lossless/TIR refusal to invent heat, weak
positive heat below endpoint-flux roundoff, and source/geometry/work/cancellation
admission. Rust execution remains unverified: the authoring runtime is unavailable.

## Two finite vapor chambers across one native moisture field (MR46)

`moisture::atmosphere::pair` couples one original slab or Cartesian field to
two independent original `FiniteMoistureAtmosphere` inventories. Every exposed
face explicitly names its first or second chamber and retains its source-pinned
dry-basis sorption and optional finite surface-transfer law; absent faces are
sealed. The gas owners retain separate volumes, carriers, composition, pressure,
transport and caloric source catalogs. Water species/molar mass and prescribed
temperature must agree. No second gas or dry-skeleton inventory is fabricated.

Each chamber balance uses only the signed native face transfers assigned to it.
Consequently, water may traverse a specimen whose net uptake is zero. The final
actual vapor changes plus represented solid-water change, and the actual total
gas-mass changes plus solid-water change, must each fit one explicit global
budget. Original per-chamber carrier/mass gates and native diffusion/source
gates remain in force. Failed trials and cancellation publish no endpoint.

The original scalar finite-reservoir solver is nested: first-chamber RH is
solved at every second-chamber RH probe, and every native field proposal starts
from the accepted field. A cumulative evaluation cap spans all nested trials;
the original native linear/nonlinear budgets still bound each proposal. All
inner roots and source states must remain admissible throughout the outer
bracket. This deterministic construction can refuse when a coupled solution
exists outside those admitted probe paths; it is not a general globally
convergent two-dimensional root solver. Worst work is the product of the scalar
trial counts, explicitly bounded by the cumulative cap.

This is isothermal, fixed-reference-geometry transport at prescribed gas
volumes. Carrier permeation, condensation, moving boundaries and temperature
feedback are not added; reported gas internal-energy changes are not total
sorption heat or thermostat duty. Focused authored G1/G3/G4 regressions cover
the independently eliminated three-inventory backward-Euler equations,
through-transfer with zero uptake, port/owner reversal, continuation, all slab
schemes, six-face Cartesian accounting, sealed/single-chamber limits, and
budget/source/cancellation/global-conservation refusal. These added native
Rust tests have not been executed in this session: DSR and the executor are
unavailable. Source review and independent floating-point balance algebra do
not establish a native compilation or test pass.

## Dispersive impedance resonance realization

`driving_point::realize_impedance_line` retains explicit resonance states only
for the analytic lossless rigid cylinder. A real-frequency viscothermal wave
number does not identify a complex pole or its residue: Bessel dispersion and
decay remain in the original TMM samples, rather than being replaced by an
oscillator at the lossless frequency. Nested FFT rungs reuse even samples and
measure the residual FIR's last-quarter energy and the first three peaks of
the actual FIR-plus-compliance response. At least two measured rungs are needed
for a finite peak-shift result. The original geometry, wall, radiation and gas
refusals remain in force. The exact lossless modal path is unchanged.

Both adaptive and explicit-size Bessel impedance callers retain the isothermal
static gain `rho*c^2/(gamma*V)` for supported sealed rigid volumes. The public
geometric helper still reports its documented adiabatic gain. The finite real
part uses the original TMM's 1e-3 rad/s quasi-static probe; mouth-pressure
transfers never receive a sealed-volume compliance state. Fixed-size callers
retain an explicitly unqualified FFT length. The adaptive caller refuses
unresolved high-Q decay or Nyquist truncation at its tap cap; no damping, gain
repair or spectral taper makes that refusal disappear. `WideTube` is a
high-shear approximation without admitted DC support and refuses this full-band
realization; the Bessel owner is required instead. These sampled numerical
checks do not certify passivity between samples, continuum error, all resonances
or experimental acoustic validity.

Focused G0/G1/G3 regressions cover a causal resonator's independent impulse
response and peak, even-bin reuse, actual tail measurements, dispersive
high-Q refusal, original Bessel isothermal storage and lossless ring-down.


## Changing-band turbulent jet-noise colour

The shared `ButterworthBandGaussian` uses a covariance-normalized two-state
realization for the flue, voice, lip and constriction noise sources. Each
zero-input step contracts the state norm even when the band changes. The
identity `F Fᵀ + G Gᵀ = I` preserves unit equilibrium state covariance for
prescribed band schedules; a zero-start or cleared source retains its warmup.
The high-pass truncation shares the band-pass's Nyquist-limit state. This fixes
the direct form's parametric instability: a native regression against the
previous core reached amplitude 41.52 after seven unit-scale noise samples.

At fixed edges the original Butterworth transfer and in-band calibration are
preserved up to rounding. Wide-band normalization uses whole/split 16-point
quadrature with at most 255 panels and refuses unresolved normalization. This
removes the former roughly 13.4% gain jump at the 1–24 kHz Nyquist limit on a
48 kHz clock. Repeated identical edges/step reuse a prepared section, removing
repeated spectral quadrature and coefficient square roots from steady-flow
sampling. Sampling allocates nothing and consumes exactly one original keyed
Philox normal per call, including silent calls.

Source placement, Reynolds gates, authored gains, Strouhal bands and retained
below-Nyquist power fractions remain unchanged. Changed realizations change
sample bits. Rapidly varying bands have no stationary spectrum or instantaneous
in-band RMS guarantee, and feedback-driven random edges need not preserve
unconditional covariance. This is a stable stochastic colour model, not a
resolved turbulence model or additional experimental noise calibration.


### Full-state driven resonance frequency tangents

`BemLoadedModes::harmonic_response_batch_with_frequency_derivative` evaluates
T*y=f and T*y_omega=f_omega-T_omega*y with one original boundary assembly and
one balanced LU. The analytic pencil derivative includes mechanical inertia,
complete Kelvin coupling, original acoustic A/B dispersion and every retained
Maxwell memory state. Geometry, medium parameters, retained basis and attachment
maps remain fixed; source phasor derivatives are optional explicit inputs.
No acoustic-block inverse, frequency perturbation or pole fit is introduced.

All returned displacement, velocity, boundary-pressure, boundary-velocity and
memory tangents retain their original physical coordinates. Differentiated
mechanical, acoustic and material rows are audited against the original
unscaled equations, with the response owner's arithmetic rounding policy.
Independent source-work, Kelvin-heat, Maxwell-heat and signed surface-work
derivatives must close. These derivatives may be negative; they are not loss
coefficients. Source and tangent amplitude normalization preserves the physical
units of every returned quantity. Original and additional work/payload limits
are separate; cancellation or refusal publishes no partial batch.

Fixed real displacement/velocity and surface-pressure/velocity receivers use
fs-modal's shared compensated `HarmonicTransferSample::project_real`, preserving
coherent cancellation for local peak/notch/linewidth refinement. This supplies
actual driven receiver slopes, not a spectral-completeness or continuum-error
claim. G1/G3/G4/G5 regressions exercise an analytic dark Maxwell mode, full
coupled frequency differences, differentiated physical work, unchanged primal
responses, phase-dependent forcing, zero-force crossings, budgets and cancellation.

## Nonlinear helium–solid thermal interfaces

The `superfluid` single-contact, paired-contact and distributed-dissipation
compositions accept `CounterflowContactLaw::Kapitza` alongside explicitly linear
interfaces. `CounterflowThermalContact::kapitza` consumes the original immutable
`fs_material::superfluid::kapitza::KapitzaInterface`: actual area, coefficient,
exponent and absolute-temperature domain remain caller-supplied. Each endpoint
may select its own law. `CounterflowThermalContact::linear` retains constant
conductance, including an exactly closed zero-conductance port.

The nonlinear heat powers satisfy the original law at BOTH finite midpoint
temperatures, together with the full reciprocal channel response and the native
caloric charts, including latent plateaux. Root probes use a monotone flat
continuation at interface temperature bounds; this is numerical bracketing only.
Within the admitted domain, the material law consumes reference temperature
and the original perturbations separately, preserving heat transport even when
their difference is smaller than one absolute-temperature floating-point ULP.
The original unmodified material law must admit both sides at the initial,
midpoint and final states before either owner publishes. No clipped temperature
or extrapolated coefficient becomes accepted physics. Scalar root work remains
bounded and uses O(1) caloric/interface probes; it introduces no repeated wave
or thermal-network solves. The paired convergence bound uses the largest
channel-side interface tangent over its declared temperature domain.

Boundary heat and entropy remain the original reference-state `Q` and `Q/T0`.
Interface wave-availability loss is `H(T_channel,T_cell)*(T_channel-T_cell)*dt/T0`,
where `H` is the selected constitutive law; it is nonnegative and is not deposited
as additional heat. Distributed viscous/vortex dissipation still enters only its
explicit native caloric destinations. Frozen small-signal helium properties,
first-order network splitting and the separate wave/reference-heat accounts
remain unchanged; this does not supply a nonlinear bulk helium EOS or phase
transition dynamics.

Focused G1/G3/G4 tests independently eliminate the two-cell wave equations,
recover manufactured nonlinear and mixed-interface heat powers, close finite
heat/availability balances, exercise bounded-domain probes and final-state
refusals, and verify cancellation/retry with profile memory and both heat
intakes. `superfluid_heat_link` now runs a sourced 1.8 K helium channel between
finite caloric cells through authored quartic interfaces, with normal-profile
memory and explicit distributed heat receivers. Its surface coefficient and
geometry are illustrative inputs, not an experimental apparatus calibration.

### Lip tension actuator work

The linear and nonlinear brass control paths retain each lip coordinate and
momentum when changing the admitted tension multiplier. The instantaneous
actuator work is the signed spring increment `0.5 (K_after-K_before) q^2`;
damping changes do not produce instantaneous heat. This increment is evaluated
independently of unchanged kinetic storage, so a large kinetic energy cannot
erase a smaller spring-work transfer. Nonfinite retained/candidate storage or
an overflowing work account refuses before publishing the rebuilt island,
control value or pending transfer. Subnormal decay whose energy rounds to zero
retains that represented zero without introducing a control refusal.

`lip_tension_work_j` retains one signed joule total per lip mode. Pending
control changes enter the next successful nonempty `BrassBlockDiag` once;
voice/performance getters also expose accepted work before rendering. The
optional second mode has no admitted nonunit tension law and therefore keeps
its original state and zero parameter-work increment. Linear voice reset starts
a new work-accounting horizon while retaining historical block diagnostics.
Performance windows and the `music_render lip` sidecar expose the same separate
actuator account; control-clock order makes it independent of audio partitioning.

Focused regressions cover spring-energy changes of both signs, rest and repeated
controls, retained states, independent modes, rejected controls/blocks, staged
changes, and scheduled output partitioning. Pneumatic work and physical losses
retain their existing meanings. This control account does not close the still
approximate lip/cup/gas swept-volume coupling or establish whole-system passivity.

### Resolved electrostatic geometry with a finite material circuit

The default-off `resolved-electrostatic-circuit-motion` feature exposes
`electrostatic::motion::circuit::PlanarCapacitorCircuitMotion`. It joins the
existing source-bound planar dielectric path, native material RLC circuit,
native articulated body impulse and bounded displacement feedback owners.
The actual field is rebuilt and solved at every changed trial geometry, so
heterogeneous permittivity and fringing inside the authored planar domain
change capacitance, charge/current and mechanical motion together.

For `H(q,s)=q²/(2C(s))`, the circuit uses the harmonic-mean endpoint capacitance.
The symmetric charge/geometry discrete gradient supplies mechanical work
`-(q0²+q1²) delta(1/C)/4`. Its zero-displacement limit uses the original field's
analytic shape derivative. This closes the variable-charge identity as well as
the original fixed-charge limit; replacing the charge-square mean with q0*q1
would violate that identity. Native source work, bulk resistor and contact heat,
unrelated capacitor/mutual-inductor storage and actual body kinetic work remain
separate observations. The physical electrical endpoint account includes the
resolved field once. The native circuit subproblem's `delta_h` uses its auxiliary
constant capacitance and is not the physical changing-geometry storage change.

Construction supplies an initial joint assembly, replacing one unbound scalar
capacitor placeholder while retaining its original charge coordinate. All other
material snapshots, graph topology, contacts, storage laws and inertia remain
fixed during evolution. A competing source-owned capacitor, nonlinear law or
coupled-capacitance block on the active branch refuses. Native constraints must
admit unchanged charge/flux at every coefficient binding; incompatible ideal
clamps refuse without projection. The existing held-voltage motion path serves
direct ideal voltage control. Cancellation, domain exits and failed independent
electrical/mechanical/total work checks return no joint candidate.

The original fixed-thickness, explicitly isotropic, strain-independent dielectric
mesh path and single guided translating body restrictions remain. Heat is reported
for explicit receiving owners, not silently deposited into a thermal state.
This does not add electrostriction, dielectric charge transport, relaxation,
breakdown, general remeshing, 3-D motion, contact events or thermal feedback.
Work admission is not a time/spatial error bound or experimental device validation.

Focused G1/G3/G4 regressions are in `tests/electrostatic_motion/circuit.rs`:
isolated-control parity, independently reduced driven RC/RLC plate equations,
separate passive storage/heat, original branch/material identity, continuing
tapered heterogeneous-field feedback, competing source and clamp refusal, and
late cancellation/retry. Run with
`cargo test -p fs-couple --features resolved-electrostatic-circuit-motion --test electrostatic_motion`.

## Complete-field finite-opening interior/exterior acoustics

`cavity_modal::interior::exterior::RadiatingInterior` couples the complete
original pressure FE field to the native finite-opening exterior BEM matrix.
It constructs `InteriorMultiport::new_full` (or its original boundary-layer
variant) from the original FE mesh, signed P1 boundary modes and coherent source
gas. Every original pressure DOF remains, subject to the caller's dense cap and
the FE hard ceiling. Original-coordinate field recovery uses the existing
Robin-loaded factorization without an eigenmode cutoff or unloaded-pole
impedance inversion.

Opening triangles must match the retained outward closed exterior exactly and
have planar support. All modes on radiating faces must be included; driven
interior feeds cannot share those faces. Different signed radiation modes may
overlap, including zero-net-flow modes. An exterior P0 mesh that erases an entire
supplied mode refuses rather than silently dropping it. Native BEM velocity
shapes remain dimensionless P1 panel means, so the unchanged source artifact
retains its generalized force/velocity units [N s/m]. The FE reaction matrix is
`conj(Z_native[r,s])/(A_r A_s)` [Pa s/m^3], using TWO geometric support areas,
never signed net-flow gains. All self and mutual entries are retained, and a
resolved native passivity deficit refuses without clipping or matrix repair.

With FE flow positive inward and exterior flow outward, the interface is
`P + Z Q_inward = 0`. Positive real power references only set numerical
coordinates. Source changes reuse the full pressure factor, radiation matrix
and small port closure. Returned amplitudes retain original FE coordinate
order. Source work pairs solved pressure with PRESCRIBED feed flow; radiated
work and native viscous/thermal wall loss are independently evaluated. Original
feed-flow and opening-pressure constraints are checked in power coordinates,
with default relative ceiling 1e-8 and an explicit tolerance entry point.
Port and source-minus-radiation-minus-wall-power residuals remain observable;
no residual is assigned to a physical heat channel.

This is a bounded harmonic model with one coherent homogeneous gas, planar
openings, a rigid remaining exterior and the selected BEM route's existing
mesh/band admission. Pressure, trace and exterior convergence are independent
obligations. Geometry/support, source identity, finite data, explicit work caps,
physical residuals and cooperative cancellation produce refusals. There is no
new dependency direction or external numerical runtime. A causal time-domain
radiation realization and the broader pexfp.99/.103/.104 in-loop acceptance
requirements remain open.

Four native regressions compare against an independently assembled original
CSR FE/BEM mixed block and separately computed BEM matrix; they exercise full
P2 coordinates, mutual loading, live zero-net modes, reference invariance,
finite-radiation versus ideal-release behavior, native wall-loss balance and
binding/refusal paths. Deliberately corrupting the cached closure must fail
original physical constraint admission. Pinned Rust parser/format checks and
independent three-port closure/unit/work algebra passed. Native execution is
still pending: the current environment kills the unchanged upstream dependency
build before these crates. Source checks are not native numerical validation.

### Typed velocity histories on native open flow faces

`scenario_velocity::sample_velocity3` executes the existing typed scenario
velocity in its declared Cartesian axes at physical time, then converts all
three signed components using the supplied frozen `MaterialScaling`. The full
Euclidean speed must fit that scaling's declared scenario bound. No material
is re-evaluated, history copied, direction inferred from a scalar, or carrier
velocity added by a coordinate-frame rotation.

`rebind_velocity_face3` sends the checked sample to the original D3Q19
`BoundaryGrid3::rebind_open_boundaries`. The caller explicitly maps the named
region to a face and supplies the same scaling used for the grid; lattice axes
are WORLD Cartesian axes. Other face data, populations, force coordinates and
time remain unchanged. The next native step applies the held velocity through
the existing force-corrected boundary reconstruction. Sampling, bounds, native
kind/profile refusal or cancellation publish no target. This is a uniform
prescribed velocity on a fixed open face, not a moving boundary, a spatial
profile, a finite reservoir or an additional time integrator.

`tests/scenario_velocity.rs` checks signed rotated/interpolated histories on a
real forced inlet, all three physical velocity components after native steps,
full-vector speed admission, unchanged state on refusal and cancellation.

## Causal exterior reaction on complete pressure fields

`cavity_modal::interior::causal::CausalInteriorTime` interconnects the original
`PreparedExactTime` FE field with an already admitted `PreparedMimoLoad` or
`PreparedDelayedMimoLoad`. Every free pressure/potential coordinate and every
native exterior pole/delay history is retained. Dimensionless P1 opening
profiles use outward velocity [m/s] and conjugate integrated pressure force [N];
the impedance is N s/m, matching native BEM aperture work coordinates. Signed
and zero-net modes never divide by net volume gain. Source rows/columns map
explicitly to original profiles; remaining feed profiles and local resistive
faces cannot overlap an exterior opening. Original volume injection and exact
pressure-release admission remain active.

Preparation reuses one complete-field exponential and caches only the small
`(H + D_mean)` port factor. Each interval matches original interval-mean forces
with held opposite interior/exterior velocities, retaining mutual reaction and
native fractional arrivals. Both owners are exact for their held inputs; this
finite-clock exchange is a discretization of continuous coupling and requires
an independent time ladder. Warm steps perform no factorization, exponential,
BEM/material query or exterior-history clone. Native finite-state admission runs
last, followed by an infallible fluid/history commit; refusals preserve both
accepted histories and step count.

Prescribed feed work, signed volume-source work, independently integrated FE
Robin loss and signed exterior exchange work are accounted once. Exterior
exchange is not inferred heat or radiated energy; returning waves can make it
negative. Original continuous/sampled passivity reports remain unchanged,
including GridOnly. Supplied realization geometry, gas, mode normalization,
causality/fit accuracy and band applicability remain caller obligations. This
interconnection does not itself bake or qualify an exterior source or establish
a coupled storage certificate.

Five native regressions specify a continuum-moment tetrahedron with independent
exact pressure, Robin heat, source work and pole history; a three-rung exchange
ladder against a separately integrated continuous coupled ODE; full P2 pressure
with signed zero-net modes and independently integrated delayed mutual
response; atomic refusal; and source/clock/face/budget admission. Rust parsing
and formatting were checked; native tests remain unrun under the dependency
memory and workspace-space blockers. Independent analytical-oracle exchange
errors were 6.500772e-6, 1.618632e-6 and 4.042733e-7 on 20/40/80 steps over 0.4 s;
these are oracle results, not native implementation evidence.

### Atomic scenario velocity face updates

`scenario_velocity::rebind_velocity_faces3` samples all explicitly selected
typed velocity histories at one physical time and publishes them through one
native open-boundary transaction. Its fixed six-slot input and output follow
`Face3::ALL`; absent slots retain their current targets and spatial profiles.
Every selected history, frame, full-vector speed bound and profile precedence
check completes before publication. A later refusal cannot leave an earlier
face at a different control time. The single-face API uses this same path.

No histories, fields or populations are copied, and no physical step occurs.
The native restriction to open faces on one axis remains authoritative. This
is atomic prescribed-velocity control, not an inferred mass-flow balance.
The native integration regressions exercise paired signed histories in distinct
frames, later-history and native-kind refusals, selected versus unselected
profile authority, empty selections and cancellation.

### Quadratic opening motion in causal full-field reaction

`CausalInteriorTime::new_with_boundary_sources` appends the original prepared
`BoundaryVelocitySource` columns after linear velocity ports, preserving the
exact FE owner's ordering and work coordinates. Both families can supply
prescribed feeds or exterior reaction modes. Midside-only motion and signed
zero-net quadratic profiles retain their complete original pressure projection,
with no re-integration, corner substitution, area normalization or additional
exterior states. Exterior indices address the combined order. Original
source/model identity, pressure-release restrictions and complete augmented-state
limits remain with the exact FE owner. Canonical physical support is checked
even for zero feed shapes; an exterior source with zero original pressure
projection refuses. The old constructor delegates with an empty prepared-source
list. Warm stepping, atomic history acceptance and original load reports remain
unchanged.

Three focused regressions specify independently integrated midside-only
inlet/opening work with simultaneous linear and volume sources, complete P2
linear-profile recovery and live zero-net midside modes with fractional delayed
mutual reaction, and identity/support/release/combined-count refusals plus
bit-identical empty-source compatibility. Parsing and formatting passed; native
tests remain unrun in this environment. This extends interconnection only:
supplied exterior realizations must describe the same original physical
quadratic modes and retain their original qualification limits.

### Typed velocity controls on shared D3Q19/D3Q27 reservoirs

`scenario_velocity::lattice3` applies the same typed SI velocity histories to
the original `OpenBoundaryGrid3<L, Q>`. The specialized D3Q19 and shared native
adapters use one fixed-size preparation function for history sampling, frame
rotation, dimensional conversion, full-vector speed bounds and profile policy.
The shared grid has no spatial profiles; the specialized owner still refuses
selected profiles and retains unselected ones.

One native target transaction follows preparation. Native grid/collision Mach
limits and boundary-kind admission remain authoritative, and `Error::Boundary`
retains the complete `GridError3<L::Error>` instead of reducing it to text.
No collision model, equilibrium, time step, force or population is replaced by
the adapter. The subsequent native step owns physical force correction and
measured reservoir exchange. Regression cases drive signed rotated histories
through actual factorized D3Q27 steps and check later-face refusal, native Mach
and kind errors, cancellation and unchanged accepted state.


## Analytic scenario frames drive fixed-topology moving walls

`scenario_wall_velocity::FrameWallVelocity3` binds a prepared analytic frame path, the grid's admitted `MaterialScaling`, and an explicit Cartesian WORLD grid origin in metres. At a caller-selected physical sample time it evaluates the frame once, maps each original native wall-link position into WORLD metres, evaluates the carrier velocity at that fixed Eulerian point, and converts all three signed velocity components through the material scaling. Query work allocates nothing: the caller supplies one scratch vector per selected wall link, and long traversals poll cancellation every 64 links. The complete speed norm must fit the supplied scaling's admitted maximum speed. Invalid samples, history corners, unit conversion failures, insufficient scratch, overspeed, or cancellation enter no fluid step; scratch is not a published control and may be partially filled on refusal.

Serial and pooled adapters then call the original shared D3Q19/D3Q27 moving-wall operator. That owner remains authoritative for probe topology, wall tangency, discrete impermeability, collision Mach limits, actual momentum/work receipts and atomic population/step-count publication. Returned `MovingWallMomentum3` retains native lattice units. The time sample is explicit, allowing midpoint sampling without assuming a solver clock. This is prescribed surface motion on fixed geometry; it does not implement normal wall displacement, swept cells or remeshing.

The integration target `scenario_wall_velocity` checks translating drives against the original moving-wall populations, impulse and mechanical work on both native lattices, independently checks rotation at fixed WORLD wall-link locations, and checks overspeed, forbidden normal motion, cancellation and scratch-length refusal. Source review and parsing passed. These integration tests are not claimed as executed while the unchanged `asupersync` dependency blocks the full native crate build.


## Distributed dielectric loss on the original finite thermal timeline

The default-off `resolved-dielectric-thermal` feature couples native planar,
axisymmetric and tetrahedral dielectric candidates to the original finite
enthalpy and complete heat-network owners. Polarization and Ohmic losses enter
thermal cells once; electric/polarization storage remains separate. Source
claims are refreshed at actual trial temperatures under fixed, validity-wide
constant electrical coefficients. Physical cell volumes, source ownership,
clocks and joint electrical/thermal energy must agree before either endpoint
is returned. Fourier transport, thermal contacts, layers and closed-enclosure
radiation retain their existing owners. See
[the dielectric thermal contract](docs/resolved-dielectric-thermal.md) for the
API, first-order coupled thermal accuracy and bounded physical assumptions.
Seven focused regressions are in `tests/dielectric_field_thermal.rs`.

## Source-bound causal radiation for complete interior fields

`cavity_modal::interior::exterior::causal::BoundCausalRadiation` binds original
planar P1 opening shapes, the complete pressure FE space, one coherent
homogeneous source gas, native full BEM aperture matrices, and the native proper
MIMO exact-ZOH load. The same opening binding serves harmonic construction.
Malformed, duplicate, empty or nonfinite original traces refuse before geometry
indexing. The bound field is assembled from the source gas itself; an arbitrary
heterogeneous coefficient field is never assumed to equal its reference medium.

Original BEM matrices remain in N s/m and are conjugated once for the FE/fitting
convention. Every self and mutual reaction and signed zero-net mode retains its
original force/velocity work coordinates. Native causality is checked before
fitting, under explicit finite-band tail assumptions. Exactly three increasing
proper shared-pole orders retain their original reports and must decrease
measured matrix-relative error up to roundoff. Independently evaluated,
interleaved BEM frequencies remain outside the fit and must also meet the
declared final tolerance. Actual sampled endpoints must equal the load's source
band and cover Nyquist. The positive first frequency leaves DC fit accuracy
unqualified; the native sampled passivity check still includes DC.

Sample storage, finite-band causality work and native relocation QR scratch/work
are bounded before baking. Native BEM mesh/route limits, pole fitting,
continuous/sampled passivity admission and complete-field time preparation
remain authoritative. The immutable owner constructs the existing time closure
from its own field, mode order and zero-history load. Sources created through
`field()` retain that exact FE identity. Cancellation precedes history cloning
during preparation; repeated simulations have independent native load histories.

This construction extracts zero prescribed delay. It retains phase propagation
in finite-band BEM data without inferring sharp arrivals from aperture centres
or bounding boxes. It applies no gain scaling or passivity repair and promotes
no native certificate class. GridOnly remains GridOnly. It does not establish
spatial convergence, a uniform unsampled error bound, sharp first-arrival
accuracy, exterior stored energy or closed-loop decay qualification. Planar P1
opening traces remain its explicit source family; the lower-level time owner
continues to accept separately qualified quadratic and delayed sources.

Six focused native regressions specify original signed BEM coordinates,
source/trace and resource refusals, a three-order reciprocal relaxation model
with independent validation, analytic phasor/causality signs, and the complete
public box bake -> original volume source -> full-field time-step route with
separate BEM and FE-work checks. Pinned parsing/format and source/API reviews
passed. Native compilation and numerical execution remain unverified under the
existing upstream dependency-memory and workspace-space blocker. In particular,
successful numerical admission of the full BEM fitting fixture is pending; no
fit, causality or passivity gate was relaxed to manufacture it.

## Isothermal moisture slab SDIRK2 diffusion

`MoistureSlab::advance_sdirk2` integrates the original conservative slab flux
with Alexander's two-stage, second-order L-stable SDIRK method,
gamma = 1 - 1/sqrt(2). Both stages use the original pinned transport source,
endpoint-dependent nonlinear iteration and native tridiagonal solve. The
selected Kirchhoff path retains its exact source-integrated face law, including
zero dry diffusivity. Ordinary `advance` retains the caller's original method.

The second equation retains the physical original field and adds the first
stage's explicitly evaluated conservative face transfers with weight
(1-gamma)*dt. It never reconstructs a flux by subtracting rounded stored
moisture. The endpoint receipt spans the full physical interval; its signed
boundary water is the weighted sum of both stage fluxes. The separate
`first_stage` receipt is diagnostic and must not be added to that account.
Endpoint nonlinear residual and water conservation checks include the known
first-stage load.

L-stability does not imply unconditional positivity. A second-stage candidate
outside the original field/active-reservoir range or admitted moisture domain
returns `MoistureError::StageRange` before querying unsupported source states.
A zero-transfer reservoir cannot enlarge that range. Failed stages publish no
state, clip no moisture and do not fall back to first order.

`advance_sdirk2_adaptive_to` reuses the original whole-duration controller.
Its full-versus-two-half comparison estimates second-order fine-step error;
unequal representable fractions u and v use (u^3+v^3)/(3*u*v), equal halves
giving one third. Accepted state and boundary water come only from the two
native fine steps. Stage-range refusals permit shorter retries; cancellation,
source-domain and inventory failures remain transactional. Each attempted
SDIRK call reserves two stage evaluations, including early refusals, and a
complete comparison costs six. Ordinary first-order callers retain their
existing error scaling and work counts.

Temperature and boundary environments remain constant during a call; callers
split at environmental events. This path covers the native slab only. It does
not establish second-order finite-reservoir, layered, Cartesian, or
nonisothermal coupling, nor a global temporal error certificate.


## Periodic driving of nonlinear physical-port feedback

`phs_orbit::solve_forced_feedback_harmonic_balance_*scaled_seeded` combines an
authored `PeriodicPortInput` with the original instantaneous
`AutonomousPortFeedback` law in the same full pHS ODE. Each residual uses
`(J-R) grad(H) + G [u_periodic(t) + u_feedback(x,G^T grad(H))]`. The externally
prescribed frequency and phase remain fixed, with no autonomous phase row or
unknown frequency. Original storage, native resistance, complete port order
and every explicit material-memory coordinate remain present.

The existing feedback tangent owner supplies `(J-R) H'' + G du_total/dx` to
the shared scaled/adaptive AFT solver. Its declared total feedback derivative
includes the physical flow/Hessian chain once; the periodic drive is
state-independent. Native complete storage tangents are discovered, with
explicit complete Hessian actions available for other laws. A refusing or
incomplete declared action never falls back to finite differences. Forcing
bandwidth is admitted before source preparation, including arbitrarily small
nonzero lines; zero padding does not consume resolution. All callbacks are
deterministic and read-only, with physical memory in explicit states.

`audit_forced_feedback_orbit` evaluates the original storage gradient once per
quadrature sample and retains separate signed external-drive and feedback
work through the original compensated physical input map. Negative feedback
work is net absorption by that authored law, separate from native R loss; it
is not automatically material heat. The combined original input supplies the
cycle power residual. Separately rounded large opposing component averages
need not sum bit-for-bit to that combined account.

`solve_forced_feedback_shooting_variational` and its `_analytic` variant expose
that same complete ODE and total derivative to the original discrete RK4
variational owner. Period is exactly the external source period, without an
autonomous phase row or period column. The complete state, monodromy and Floquet
spectrum retain material memory. Adapter and native solver resource costs share
the caller's byte/work caps; original forcing bandwidth is checked before
source preparation, and cancellation precedes preparation and is polled at
native stage/product boundaries. A declared tangent refusal never changes to
finite differences. No second integrator or constitutive law is introduced.

The bounded finite-order HB lane has no cooperative Cx cancellation. Neither
its sampled balance nor one discrete shooting map establishes an infinite
harmonic tail, physical-law validity, global stability or continuous-orbit
accuracy. Audit quadrature must resolve drive and retained state orders;
nonlinear work additionally needs sampling refinement. Focused G1/G3 tests
cover exact phased multiport quadratic/cubic harmonic conversion under unit
rescaling, retained relaxation memory, passive asymmetric drag versus native
time-domain shooting at different phases, the analytic discrete scalar Floquet
multiplier, nonlinear storage and signed power splits, and strict bandwidth,
source, tangent, resource and cancellation refusals.

## Finite wall heat on ideal and nonideal elastic gas cavities

`cavity_phs::modal::GasValveMotion::{thermal_candidate,thermal_layers_candidate}`
joins the original nonlinear gas/elastic pressure-work step to one finite
material enthalpy body. The existing bounded scalar feedback owner solves one
signed transfer Q into the gas and out of the body. Every trial restarts the
accepted gas and original modal motion, with gas heat `Q_gas+Q` and body heat
`Q_body-Q`. Its interface sees the actual gas energy/volume midpoint and body
enthalpy midpoint, including any latent fraction. No independent piston, EOS,
heat-capacity approximation, modal reduction or thermal root finder is added.

Both native ideal and Peng--Robinson gas parcels retain their original
composition, mass, source conditions, branch admission, pressure, density,
calorics, compressibility and transport. The series-layer entry point uses the
same gas transport state and the body's original material card, condition axes,
conductivity query and receipts. An ideal reference state cannot substitute
for the actual dense-gas midpoint.

A `ModalGasThermalCandidate` returns the complete original motion/gas owner,
`UniformEnthalpyStep`, and `ModalGasThermalStep`. Acceptance independently
checks the native gas/solid pressure work, the body storage residual, the heat
law, and the enclosing balance
`dU_gas+dE_modal+D_modal+dH_body+W_exterior-Q_gas-Q_body`.
Internal Q cancels once. Native modal dissipation remains external to both
thermal stores. All heat trials use the same physical sample period; source,
clock, bracket, work-budget or cancellation failures return no partial state.
The caller publishes both candidate endpoints together.

The thermal participant follows the original uniform-temperature,
fixed-geometry, negligible-pV-work enthalpy convention. It does not silently
rebuild the source-fixed elastic basis, change stiffness, add thermal expansion
or transfer mechanical state across a phase transition. The existing gas-loaded
sample-rate and midpoint-resolution guards remain active. Temporal refinement
is still required; no unconditional convergence, monotone relaxation, entropy
theorem, measured-fluid calibration or automatic exterior-radiation refresh is
claimed.

Six G1/G3/G4 regression definitions reuse the original six-tetrahedron modal
and Peng--Robinson fixtures: independent constant-Cv two-storage heat and
compression elimination; continuing dense-gas/material-layer exchange and
whole-system energy; latent-enthalpy midpoints; exact insulation reduction;
late source/domain/clock/budget refusal with deterministic retry; and
cancellation after actual interface evaluation. Focused native command:
`cargo test -p fs-couple --test solid_enthalpy valve::gas_motion::nonideal::thermal`.

## Causal reflected volume stresses on the native microphone/PCM path

`pcm_wav::observation::quadrupole::TransientQuadrupoleRenderer` borrows an
original `PreparedTransientQuadrupoleObserver` and read-only physical source
histories. The observer retains each direct/image path's actual retarded time,
full tensor transformation and coherent pressure/particle-velocity sum. The
existing ideal `FirstOrder` microphone reduces that physical field, and the
existing `PressureRenderer`, decimator and WAV stream retain their own clocks,
filter delay and explicit Pa-to-PCM scaling. No second propagation, interpolation,
source evolution, gain normalization or audio encoder is introduced.

`QuadrupoleSourceWindow` declares each original history's closed DD-clock
support and bandwidth. The history reader supplies consistent stress moment,
rate, acceleration and integral, including its original integration constant;
the same source can be read at different times for its images. Preparation
checks the complete receiver window against every source's shortest/longest
propagation delay, source and receiver bandwidths against strict Nyquist, and
physical-source, all-image work and scratch-storage budgets before reading any
history. A bandwidth declaration is not an interpolation-convergence or source
compactness certificate. Missing reflected prehistory/tails and unresolved
timestamps refuse; no zero history or artificial tail is synthesized.

One bounded field buffer serves all callbacks. Integer-index timestamps use
the native `HarmonicSampleClock` without imposing a harmonic source assumption,
and emission times retain DD precision through the original transient observer.
`receiver_time_at` exposes that full clock; the field sample's f64 time is an
inspection value. Field and microphone callbacks share one accepted cursor and
are partition invariant. History refusal or cancellation leaves both caller
output and accepted sample count unchanged. Readers must be read-only: hidden
history mutation is outside this adapter's atomicity contract. Runtime observation
adds no allocation; borrowed history-query costs remain with their original owner.

Five focused regression groups cover original harmonic full-tensor and vector
field parity, ideal-microphone directionality and instantaneous `p*u`, split
field/pressure callbacks, existing PCM streaming and causal decimation, large
epoch translation, callback rollback and complete reflected-history admission.
These controls are authored; a complete native dependency build has not been
executed in the current constrained environment. The original stationary ideal
rigid/release boundary and homogeneous lossless propagation scope is retained;
this observer does not add radiation feedback, moving geometry or finite wall
impedance.


## Native joint friction with finite-cell thermal feedback

`friction_heat::joints::propose_joint_friction_heat` returns one private candidate containing the original native rigid bodies, constrained joint-friction flow, original finite enthalpy endpoint and source receipt. Every nonlinear trial starts from the same accepted states. Source capacities/drives and the two original conductivities are resolved at the same candidate interface temperatures; Newton-Krylov and fallible directional differences retain their original convergence and refusal semantics.

Each selected hinge or slider has one distinct existing thermal link. Canonical link order and original joint identities are checked independently. The source explicitly owns the reduced bearing/rail-to-thermal correspondence, load distribution and applicability: an aggregate bilateral reaction is not treated as a bearing preload. Multiple interfaces can share finite cells, whose original storage is updated once. Their common-temperature coordinates and material/geometry conductances determine signed heat into each body, including ordinary hot-to-cold transfer.

Heat comes from each original joint's continuous friction-work integral through arrest and reversal. An outer endpoint impulse pairing is not a per-source heat measure. The source declares a heat fraction in [0,1]; the remaining physical loss is an explicit nonthermal energy export, with separate per-port and total allocation residuals. Neither constraint work nor numerical residual is thermalized. The combined gate closes actual native kinetic change plus caloric change and exported nonthermal work against axial drive work, external heat and separately reported bilateral residual work. An optional original entropy audit covers the caloric/converted-heat subsystem; it does not assign entropy or state to the nonthermal exports.

This is a fixed-pose, fixed-thermal-geometry phase with held efforts, common static/kinetic Coulomb bounds and backward-Euler temperature continuation. Time/spatial refinement and the surrounding native pose owner remain required. No wear evolution, surface recession, thermal remap, distinct static coefficient or new heat store is inferred. Source refusal, invalid identity/domain, exhausted bounds and cancellation return no published body or caloric state.

Direct regressions cover the reversing slider's 7/18 J continuous loss (whose separate outer friction impulse pairing has the wrong positive sign), an independently solved temperature/friction quadratic, heat fractions 0/0.25/1 with signed conduction, two coupled hinges sharing one finite middle cell, zero heat during holding, and atomic source/identity/budget/cancellation refusals.
