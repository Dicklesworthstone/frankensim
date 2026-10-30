# fs-gas contract

`fs-gas` owns L3 evolving gas systems. Constitutive and species laws remain in
`fs-thermochem`; this crate does not replace those sources with an independent
caloric model. Production dependencies remain pure Rust and Franken-owned.

## Experimental frozen-composition Euler flow

Feature `experimental-euler` exposes `dynamics::{FrozenGas, Euler1d}`. It is an
Estimated engineering baseline for one-dimensional, inviscid, stationary variable-area flow
of a single frozen ideal-gas mixture. It advances actual compressible waves,
contacts, shocks and rarefactions without the low-Mach lattice assumption.
It does not grant the separate entropy/invariant-domain lane's theorem claim.
This implements the first transport slice of
`frankensim-ext-gas-1d-network-vvnx`; that broad network bead remains open.

### Source and thermodynamic state

The immutable `FrozenGas` owns the complete original NASA-9 mixture and a
caller-declared temperature interval. All active species must support that
interval without gaps. Crossing source regions requires the original source's
explicit caloric-continuity construction. Bounded directed-rounding interval
subdivision of each component's original Cv/R polynomial admits positive Cv
over the whole path. Failure to establish that bound refuses construction.
The subdivision uses a fixed-size stack and at most 4095 work items per region.

Pressure is `rho R T`; frozen sound speed is `sqrt((Cp/Cv) R T)` at the same
temperature and composition. Temperature inversion uses the original source's
reference-invariant internal-energy difference and analytic Cv, with safeguarded
Newton/bisection and explicit absolute/relative energy tolerances and iteration
limits. There is no constant-gamma substitution for a variable-Cp source.
Full source-bearing evaluations remain available through `thermodynamics`.

Conserved cells contain mass density [kg/m³], axial momentum density
[kg/(m² s)] and referenced total-energy density [J/m³]. The energy variable is

`rho * (e(T) - e(Tmin) + u²/2)`.

The original signed `e(Tmin)` remains in the owner. Physical total energy is
recovered by adding `e(Tmin)*rho`, and physical energy flux by adding the same
reference times mass flux. This exact linear change of conservative variables
avoids rounding away sensible energy under large formation references. It is
valid for the declared frozen composition and cannot silently transfer to a
reacting or composition-changing cell. No positivity requirement is imposed on
the original absolute formation/internal energy.

### Conservative transport and boundaries

The numerical flux is HLLC using both states' actual frozen acoustic speeds.
Contact intermediate states must have positive density/pressure and internal
energy inside the admitted source interval. Degenerate or inadmissible contact
states use an HLL fallback; its intermediate state is checked, with bounded
wave-speed expansion if required. These local checks are not a proof of global
entropy stability or positivity for arbitrary EOS data.

`NormalConserved` and `FrozenGas::hllc_normal_flux` expose the same closure to
two-dimensional spatial owners in an orthonormal normal/tangent frame. They
carry both momentum components and full referenced total energy. Transverse
momentum and kinetic energy follow the same HLLC/HLL waves as the normal
equations; HLL intermediate-state admission includes heat produced by numerical
shear mixing. No extra constant-gamma source or independent passive-scalar flux
is introduced. `normal_conserved` refuses loss of sensible-energy precision
when adding transverse kinetic energy, and `normal_primitive` removes that
energy exactly once before the original caloric inversion. The stationary slip
wall adapter reuses the thermally perfect shock/rarefaction wall solve and
returns zero mass, transverse momentum and energy flux. Spatial reconstruction,
mesh evolution and timestep acceptance remain with the consuming spatial owner.
Focused G0/G3/G4 tests cover axial correspondence, moving shear contacts,
transverse-boost covariance on both flux branches, shear-heating source-domain
refusal, slip-wall pressure and cancellation. This adds no entropy theorem,
viscous wall model or physical shock-thickness claim.

`NormalPrimitive` and the `*_normal_flux_from_primitive` entries admit
reconstructed density, normal/tangential velocity and pressure through the
same source and conservative-energy precision checks. They reuse the original
source-computed thermodynamics in the unchanged HLLC/HLL and wall kernels,
avoiding an immediate inverse temperature solve on energy just constructed
from those same primitives. Callers cannot supply a cached sound speed or
substitute a caloric model. Two added regressions compare conservative and
primitive entries for contacts, supersonic flow, HLL shear mixing and both wall
branches, plus source/precision/cancellation refusals. All 154 focused native
gas library tests pass on the pinned nightly.

`SpatialConserved`, `SpatialPrimitive` and `SpatialFaceFlux` extend that same
closure to a normal and two orthogonal tangent directions. The state carries
`rho * (e(T)-e(Tmin) + (u_n²+u_t1²+u_t2²)/2)`. All three momenta and the full
energy participate in one HLLC/HLL solve. A private const-generic kernel is
shared by the planar and volumetric entry points, including star-state
admission, bounded HLL wave-speed expansion, and physical flux evaluation;
it does not fork the Riemann model or append independent scalar advection.
The three-dimensional constructor, decoder and primitive fast path preserve
the same original source and caloric-precision admission. Wall pressure uses
the existing thermally perfect shock/rarefaction closure, with zero mass,
both tangential fluxes and work at a stationary adiabatic slip wall.

The spatial face tests check bitwise planar reduction when the second tangent
vanishes; double-shear contacts; rotations and Galilean boosts of both tangent
components on HLLC and HLL branches; combined shear-heating domain refusal;
tangentially invariant wall pressure; source-owned primitive/conservative
agreement; and precision/domain/cancellation refusal. This is the closure for
the multi-D owner in `fs-flux`, not a three-dimensional mesh integrator here.
The Estimated claim and the entropy, vacuum, wall-model and physical-validation
limitations remain unchanged. All 162 focused native gas library tests passed
on the pinned nightly through the original dependency sources; the complete
workspace DSR gate was unavailable and is not implied by this narrower result.

MUSCL uses monotonized-central primitive slopes, reduced together when needed
to retain positive density/pressure and the declared source interval. SSP RK2
uses one shared oriented flux per face. Periodic end faces share the same flux
bits; reflecting stationary walls impose exact zero mass and energy flux;
transmissive boundaries use zero-gradient ghost states and are not claimed to
be nonreflecting subsonic outlets.

Every stage checks its actual face-speed CFL against the caller's value, at most
0.45. Failed trials halve dt within an explicit retry cap. CFL-only retries
retain the requested spatial reconstruction; actual admissibility failures use
a first-order retry. Every accepted stage and final cell must remain in the
source domain. No conservative variable, temperature or pressure is clipped.
`StepReport` retains the actual time increment, maximum stage CFL, fallback
counts and extensive mass/momentum/energy balances against the two RK stages'
integrated boundary fluxes.

Workspace is allocated once, with a hard cap of 1,048,576 cells and fallible
allocation. Successful steps allocate no numerical workspace. Cancellation is
polled through cells, faces, source admission, inversion and retries. A refusal
or cancellation preserves the complete original cell field and time. Operations
are ordered deterministically on one arithmetic target; cross-ISA bit identity
and parallel scaling have not been measured.

### Focused native evidence

Eight tests of the complete source compiled and executed through the actual
pure thermochemistry dependency closure. They cover source-gauge invariance
under ±1e18 formation offsets, variable-Cv/high-speed EOS round trips, HLLC
stationary contacts, uniform state and reflecting-wall conservation, strong
rarefaction positivity, source-domain refusal and early/mid/late cancellation.

Measured regression results on the authored fixtures:

| Case | Resolution | Result |
| --- | --- | --- |
| Smooth contact advection, density L1 | 40 / 80 / 160 cells | 7.98344e-4 / 1.96127e-4 / 4.72099e-5 |
| Sod tube, independent exact gamma-law oracle, density L1 | 120 / 240 cells | 4.01893e-3 / 2.20519e-3 |
| Lax tube, same independent oracle, density L1 | 120 / 240 cells | 1.48205e-2 / 9.28277e-3 |
| Variable-Cp acoustic wave | 128 cells | Measured speed 1.15993219 vs source speed 1.15966702 m/s |

The independent oracle uses shock Hugoniot and rarefaction-invariant relations,
not production HLLC/source-inversion code. A nonlinear variable-Cp shock also
executes with conserved boundary accounting. These are manufactured/conformance
tests, not experimental validation or a performance benchmark.

`cargo run -p fs-gas --features experimental-euler --example shock_tube -- 240`
runs a one-metre N2 tube with the 200–1000 K NASA-9 coefficient region already
retained in the thermochemistry fixtures, and emits the resulting field as CSV.
The native run advanced 292 steps to 4.739365351855e-4 s, with zero reported mass
balance residual and -7.68e-15 J referenced-energy balance residual. The material
association is explicitly declared by the example; this is not a validated
nitrogen shock-tube experiment.

### Scope boundary

No evolving species transport, chemistry/flow splitting,
junctions, valves, viscous stresses, heat conduction, moving mesh, vacuum state,
real-gas EOS, multidimensional flow, adjoint, entropy certificate or experimental
validation is supplied by this Euler slice. A source-domain refusal is an actual
refusal and does not extend a NASA fit to hotter or colder shocked gas.

References: Toro, Spruce and Speares, *Shock Waves* 4 (1994), 25–34,
DOI 10.1007/BF01414629; Ketcheson, LeVeque and del Razo, *Riemann Problems and
Jupyter Solutions* (2020), Euler and approximate-Euler chapters.

## Stationary area and nonuniform-mesh extension

`Mesh1d` retains ordered physical face positions, positive face areas, and
integrated cell volumes. `Euler1d::new_on_mesh` consumes that geometry without
resampling or substituting a midpoint volume. The uniform constructor remains
available. `cell_width()` and `cross_section()` return `Option` so a nonuniform
or tapered field cannot advertise a nonexistent uniform dimension.

MUSCL slopes use the actual distances between cell centres. Each physical face
flux is multiplied by its area exactly once and shared between neighbouring
volumes. The pressure-area momentum source is paired algebraically with the
pressure traction, preserving a resting uniform-pressure state even across
density contacts in a taper. Stationary area changes do no energy work.

Both stages enforce a local area/volume-scaled face CFL, including any HLL
wave-speed expansion. `StepReport::area_momentum_impulse` records the geometric
force integral. The residual subtracts both boundary transfer and this impulse.
Periodic boundaries require equal end areas; they share the same flux bits.

Conservative state construction also checks that stored momentum and energy
retain the requested sensible energy within the existing caloric tolerance.
Finite values alone do not admit a state whose kinetic energy has rounded away
its thermal energy.

Fifteen focused native gas tests pass, including the original eight tests and:
nonuniform-taper rest with density contacts; finite-amplitude tapered pulse
mass/energy/momentum closure; mesh and periodic-boundary refusals; cancellation;
caloric-round-trip refusal; and a moving variable-Cp Mach-3 shock whose field and
reports are bitwise invariant under a +1e18 formation-reference shift.
The latter checks front position against independent Rankine-Hugoniot motion
and checks shock heating and entropy increase from the original source.

For smooth periodic contact advection on a stretched mesh, density L1 error
falls from 1.66346e-4 at 40 cells to 3.67860e-5 at 80 cells (4.52-fold).
This is a manufactured convergence result, not a throughput benchmark.
Source geometry/quadrature accuracy and the quasi-one-dimensional smooth-area
model remain explicit assumptions; abrupt junction losses are not supplied.

## Physical reflecting-wall Riemann boundary

Reflecting ends now solve a one-sided thermally perfect Riemann problem
directly. Incoming gas connects to a stopped wall state by the source-owned
shock Hugoniot; the velocity jump is factored from the root equation so weak
shocks do not converge to an unchanged spurious root. Retreating gas connects
to the wall by an isentropic rarefaction, integrating Cp(T)/c(T) in sqrt(T)
coordinates and recovering pressure from the source entropy increment.
Quadrature is bounded, cancellable, and split at original caloric seams.

The wall imposes exactly zero mass and energy transfer and positive normal
pressure. It does not reuse an HLL ghost-state momentum flux that can become
artificially tensile during expansion. Unsupported hot or cold wall states
refuse; the source is not extended to an unadmitted vacuum. The actual wall
wave and acoustic speeds enter the same stage CFL calculation.

Twenty-one focused native gas tests pass. Six added checks compare weak and
strong constant-Cp shock/rarefaction pressures with analytic relations, check
variable-Cp Hugoniot and entropy balance, compare expansion with an independent
2,000-panel quadrature, exercise orientation/formation-reference invariance,
and verify source-window refusal and cancellation. These qualify the
one-dimensional inviscid boundary implementation, not viscous wall layers,
heat transfer, physical shock thickness or a fully discrete entropy theorem.

## Adaptive sealed finite-rate reactors

`chemistry::{advance_adiabatic, advance_adiabatic_energy}` advances an actual
caller-declared duration at fixed volume and total internal energy. The former
obtains the energy from the original feed temperature; the latter accepts the
actual stored energy from another conservative owner. The complete original
`IdealGasKineticsV1` mechanism, source receipts, canonical species amounts,
elapsed time, final pressure/temperature and balances remain in the result.
This implements a bounded homogeneous reactor slice of
`frankensim-ext-gas-chemistry-ladder-paqh`, not the complete combustion ladder.

The solver advances reaction extents with damped backward Euler (the default)
or optional second-order L-stable SDIRK2, both with adaptive step doubling. Every candidate reconstructs amounts from the original feed
plus the exact admitted stoichiometric matrix; there is no independent species
clipping, positive trace seeding, or equilibrium projection. Temperature solves
the original total internal energy at each candidate composition, with the
formation contribution present exactly once. No separate reaction heat source
is added. The source-path/positive-Cv interval admission is the same L1 helper
used by the Euler EOS, including possible products absent from the feed.

The Newton Jacobian combines the analytic concentration and temperature rate
derivatives with `dT/dxi_j = -sum_i(nu_ij u_i)/sum_i(c_i Cv_i)`. A bounded scaled
partial-pivot solve and positivity/energy backtracking support stiff reactions.
Stage, Newton and root scratch is allocated once and reused. Local error control
compares the complete step with two half steps, covering species and temperature;
the two-half-step result is accepted without extrapolating species. Backward
Euler is first order; SDIRK2 is second order for smooth homogeneous reactor
ODEs. Both use estimated local errors rather than global accuracy certificates.

Caller controls separately bound absolute energy, absolute amount, relative
species/temperature error, relative element/mass balance, numerical entropy
decrease, minimum/initial/maximum time step and all nonlinear/root/attempt work.
Every accepted step and the complete integration check actual ideal-mixture
entropy. Final element, molecular-mass and energy balances are checked again
before publication. A budget exhaustion, unsupported source state, arithmetic
failure or cancellation returns no partial integration. Zero duration preserves
the admitted initial state. No ignition delay or empirical rate accuracy is
inferred from successful integration of a caller-authored mechanism.

The focused reactor fixtures cover analytic reversible-isomer time evolution,
stiff relaxation to the independent Gibbs/UV endpoint, exothermic energy and
entropy balance, molecular-count-changing dissociation, extensive scaling,
one-correction Newton convergence, overflow/refusal and cancellation. Source
interval admission also rejects negative interior Cv despite positive endpoint
values. These are numerical/constitutive conformance tests, not experimental
validation of a real-fuel mechanism.

## Common material gas integration

Optional feature `material` exposes
`material::advance_closed_parcel_adiabatic`. It evolves an existing
`fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel`, using its
actual stored internal energy as the one target. The complete product bank
must match the mechanism's original NASA-9 models and retain every active
feed transport card, source-pack identity and declared water identity.

Generated amounts pass through `ThermochemicalGasModel::with_mole_amounts`,
which uses the existing exact-simplex closure and refuses unrepresentable
positive traces. The existing material composition/energy replacement resolves
all final transport and acoustic properties at the preserved mass and volume.
The final material state is checked against both the original stored energy
and the kinetic endpoint, including mass, density, temperature, pressure and
composition-roundoff discrepancies. Mass and volume retain their original
bits; zero duration retains the entire original parcel. Publication follows
the last source, balance and cancellation check.

Five material-consumer regressions exercise finite elapsed reaction against an
analytic solution, reaction-driven heating and pressure rise, updated common
transport/acoustics, exact storage preservation, source/water identity refusal,
zero-time identity, early/late cancellation, exhausted work and final transport
domain refusal. The numerical reactor and source model are the production
implementations; the bounded native material test closure selects original
material storage/resolver/caloric bodies and excludes unrelated runtime and
mechanical modules. It does not establish a full `fs-material` or workspace
build. The existing workspace/runtime resource limits remain applicable.

Across the final focused native runs, all eight reactor tests, fifteen Euler
tests and five material-reactor tests passed. Eighteen selected original
material tests also passed, including the two generated-amount regressions;
the complete original thermochemistry library passed all 111 tests. The
material-consumer run compiled the actual `fs-gas` library with both features
enabled against the bounded original material source closure described above.

The Euler field remains frozen-composition. The zero-dimensional reactor does
not add species advection or chemistry/flow splitting to that field. Neither
path supplies diffusion flames, generic third-body/falloff mechanisms, spray,
knock, real-gas phases, moving volumes or a validated ignition model.

### Material cards in the shared molecular transport law

With both `material` and `experimental-euler` enabled,
`material::transport::MaterialTransport<L>` implements the existing
`multicomponent::dissipation::TransportLaw`. This connects the original
`ThermochemicalGasModel` transport cards to every consumer of that shared law,
including molecular species flow, reacting flow, finite-plenum networks and
the frozen Cartesian `fs-flux` consumer. No separate mixture formula or Euler
material registry is introduced.

Coefficient ownership is explicit: the original material evaluator supplies
viscosity and thermal conductivity using its source selection and Wilke/WMS
owners; the separately supplied `L` supplies only the complete canonical
binary-diffusivity matrix; bulk viscosity is an explicit nonnegative input.
The scalar coefficients returned by `L` are deliberately unused and never
serve as missing-card fallbacks. Pair transport retains its own source domains
and refusal behavior. Material cards do not imply binary diffusivity or bulk
viscosity. Admission compares the entire original NASA-9 species models, not
their labels, and requires the material bank to cover the complete pair axis.

The adapter retains the original material bank, cards, selection policies,
water/pack identities, complete external `QueryPoint`, and binary law. Each
actual face temperature and pressure replaces only those two coordinates in
the query, with canonical absolute-temperature/pressure meanings. Other axes
and their quantity conventions remain unchanged. Exactly normalized canonical
mole fractions pass through the original `with_composition` owner unchanged.
Native conservative faces may close a compensated sum while their canonical
ordered sum differs from one by roundoff. Only the native `64*epsilon*n`
unit-sum envelope is admitted; the original `with_mole_amounts` generated-simplex
owner reconciles those values. The actual resulting canonical fractions are
then supplied to BOTH material and pair-transport owners. Genuinely nonunit
inputs refuse. Exact-zero species are omitted from material queries; no trace
floor makes an absent species appear. The full pair matrix remains present because the
original diffusion owner supports evolution into those species. A newly active
species must have its real transport card support at the actual state.

`resolve` exposes the original resolved material state and query receipts for
the same face inputs. Source dimension/state admission precedes both material
and diffusion queries. The caller's binary output is copied only after complete
material resolution, successful original pair evaluation, shared matrix
validation and a final cancellation poll. Failed or cancelled evaluations leave
that output unchanged. Material refusals remain structured as
`DissipationError::Material(MaterialError)`, including through the existing
reacting-flow transaction. The adapter contains no accepted fluid state.

The original material resolver allocates and evaluates source cards at each
call; this path does not inherit the allocation-free claim of `KineticTransport`.
Five focused G0/G1/G3/G4 tests exercise complete query/source retention, zero
species, coefficient ownership, independent Fourier/viscous mode responses to
data-only card changes, no-query input admission, late output cancellation and
binary refusal, canonical three-active-species/one-zero roundoff reconciliation
through the same material and binary composition, and whole reacting-flow rollback when a new product's card
cannot support the final temperature. Fixtures are manufactured engineering
data, not measured gas validation. Existing source domains, neutral ideal-gas
assumptions, diffusion omissions and each consuming solver's numerical limits
remain in force.

## Reusable finite-rate reactor for cell-wise chemistry

`chemistry::IsochoricReactor` owns one original immutable kinetic mechanism and
one declared set of numerical controls. Construction admits the complete
caloric source interval and allocates the original source, stoichiometric,
Newton, adaptive-stage, concentration and endpoint workspace once. Subsequent
`advance_adiabatic` and `advance_adiabatic_energy` calls take dense `Amount`
values in canonical source order and reuse that storage without successful-path
allocation. This removes repeated identity lookup, source-path admission and
workspace construction when processing a finite-volume field.

The returned `IsochoricKineticsView` borrows the accepted endpoint and its
source owner. Explicit `to_owned(cancelled)` alone constructs the independent
result and original NASA-9 evaluation receipts. The existing allocating free
functions delegate to the same evolution algorithm. Each new solve resets
extents, extensive inputs and statistics; cancellation, an unsupported energy
target or a work-limit refusal cannot contaminate the next cell.

The previously declared implicit-Euler/step-doubling accuracy, positivity,
energy, element, mass and entropy boundaries remain in force. Focused native
verification passes all eight existing reactor regressions and four prepared
workspace tests covering analytic independent parcels, retained buffer
addresses, explicit-energy/source identity and recovery after refusals. This
is not a full-workspace or DSR verdict.

## Conservative transport of evolving ideal-gas species

Behind `experimental-euler`, `multicomponent::MixtureGas` owns one complete
original `IdealGasEquilibriumV1` species bank and an explicitly admitted
positive-Cv caloric interval. Every cell evolves all species partial mass
densities, axial momentum and total energy. Density is the compensated species
sum; temperature, pressure, Cv and frozen sound speed use the current
composition and original NASA-9 laws. Zero species stay admissible without
clipping, seeding or normalization.

The stored linear energy gauge is `rho E - sum(rho_i e_i(Tmin))`. Shared species
fluxes carry their corresponding original energy references, preserving total
physical energy while avoiding loss of sensible heat to large formation
constants. A two-species regression confirms bit-identical numerical transport
under independent formation-reference shifts of magnitude 1e18.

`MixtureEuler1d` uses the stationary `Mesh1d` owner, including nonuniform spacing,
face areas and integrated volumes. Its MUSCL/HLLC, bounded HLL fallback, SSP RK2,
actual stage CFL, pressure-area force, periodic/transmissive/reflecting boundaries
and transactional retries follow the existing engineering flow lane. Species,
momentum and energy use common face wave speeds. `MixtureGas::hllc_flux` exposes
the same source-bound Riemann closure with caller-supplied output storage.
Contiguous state, stage, slope and face workspaces are allocated once; successful
numerical steps allocate no workspace. Source or numerical refusal and
cancellation publish neither a partial field nor a new clock.

Six focused native regressions cover variable-Cp mixture thermodynamics and
sound speed; stationary differing-caloric contacts on a tapered nonuniform mesh;
smooth species-contact refinement; agreement with the frozen-source shock
solver at uniform composition; formation-reference invariance; and rarefaction,
structural zeros and late cancellation. Smooth contact L1 errors at 40/80/160
cells are 5.266259985e-4 / 1.198517512e-4 / 2.972969663e-5. The complete focused
pure-gas library passed 33 tests when this transport increment was admitted.

This remains an Estimated inviscid engineering lane. Fully conservative moving
material interfaces with differing caloric laws can exhibit pressure
oscillations; exact moving-interface pressure-equilibrium preservation is not
claimed. There is no fully discrete entropy-stability certificate, physical
species diffusion, viscosity, moving mesh, arbitrary junction or multi-D flow
claim.

### Shared planar faces for evolving species

`MixtureGas` also admits borrowed `NormalCell` and `NormalPrimitive` traces.
The canonical output axis is every original species, normal momentum,
transverse momentum, then full species-referenced energy. Its axial and planar
consumers share one HLLC/HLL implementation; transverse momentum is not appended
to an independently solved one-dimensional flux. Both kinetic components enter
star-state and HLL intermediate-state admission through the original caloric
owner. In particular, shear mixing cannot hide heat outside the source interval.

The reconstructed-primitive entry reuses thermodynamics computed by the original
source constructor, avoiding an immediate second temperature inversion. All
successful face evaluations use bounded stack scratch and caller output storage.
Stationary slip walls retain the original composition-local shock/rarefaction
pressure; species, tangential momentum and energy wall fluxes are zero. Output
is scratch until success; neither input trace nor the species source is mutated.

Five focused native regressions cover exact axial parity, analytical material
contact transport, transverse-boost covariance on HLLC and HLL, formation-gauge
invariance, source refusal under shear heating, primitive/conservative parity,
wall pressure and cancellation. These are source-kernel claims, not a spatial
convergence, entropy-stability, chemistry or physical diffusion theorem.

## Coupled finite-rate multicomponent flow

`multicomponent::reactive::ReactingEuler1d` binds one conservative mixture field
to an identical original species bank in `IdealGasKineticsV1`. It composes an
isochoric chemistry half step, a complete conservative transport step, and a
second chemistry half step. If heat release changes acoustic speed enough for
transport to shorten the trial duration, the complete composition restarts
from the original field using the actual shorter duration. Failed chemistry,
source reconstruction, work budgets and cancellation publish no cell or clock.

Each chemistry solve uses the actual physical cell volume and total internal
energy after removing kinetic energy. Momentum and physical total energy are
held fixed. The stored energy changes only by the negative species-reference
change: `Delta Eg = -sum(e_i(Tmin) Delta rho_i)`. Original formation energy is
therefore accounted for once; no additional heat-release term is added.
The accepted reactor temperature is independently checked against the original
species sensible-energy functions at the reconstructed cell composition. Lost
absolute-energy precision, molecular-mass inconsistency that changes kinetic
energy, or an endpoint outside the source/energy budget is refused.

Two complete preallocated flow instances provide accepted and candidate state;
one prepared reactor and dense amount buffer are reused across cells. Source
identity/path admission and numerical workspace allocation occur at
construction. Chemistry tolerances are per physical cell; callers must choose
absolute amount and energy tolerances appropriate to the actual mesh volumes.
`ReactiveStepReport` retains transport boundary changes, area-wall impulse,
actual chemistry species/gauge-energy changes, chemistry entropy production,
solver work and mass/chemical-total-energy residuals. Chemistry entropy alone
is not a total transport/source/boundary entropy certificate.

Four coupled regressions cover a moving exothermic periodic field against the
homogeneous reactor at the actual accepted duration (including forced acoustic
restarts), independent sinusoidal advection/reaction decay and refinement,
late cancellation followed by reuse, and source-identity/chemical-energy
precision refusal. Reactive-contact L1 errors at 40 and 80 cells are
2.3273859838e-4 and 5.1358044353e-5. The complete current pure-gas native library
passes 38 tests. The five existing material-adapter regressions also pass using
actual fs-gas source with the selected original fs-material source bodies;
this is not a full-workspace or DSR verdict.

The runnable `reacting_oxygen` example uses the O/O2 high-temperature NASA7
coefficients and isolated `2 O + M <=> O2 + M` reaction from GRI-Mech 3.0 as
distributed in Cantera v3.2.0 `data/gri30.yaml`. NASA7 embeds exactly in NASA9
with zero inverse-temperature Cp terms. It retains the original 1000..3500 K
region, 1 atm standard-state pressure, explicit SI rate conversion and unit
O/O2 collider efficiencies. This isolated oxygen subsystem is not a complete
air or fuel-combustion mechanism and carries no experimental-validation claim.
A 40-cell periodic run reaches 1 ms in 113 accepted steps, with reported mass
residual 0 kg and chemical-total-energy residual -1.718625e-13 J; atomic oxygen
falls from 4.237578e-5 to 3.380236e-5 kg. Local temperatures reach
1831.691..1912.266 K from an initially uniform 1800 K.

Run with `cargo run -p fs-gas --features experimental-euler --example
reacting_oxygen -- 80`. The example emits coherent-SI CSV and a balance summary.
The symmetric composition does not upgrade the underlying adaptive
implicit-Euler chemistry to a certified second-order time method. Local
chemistry error estimates are not global coupled error bounds. Physical
species diffusion, heat conduction, viscosity, flames, sprays, real-gas phases,
moving meshes and arbitrary reacting junction networks remain separate work.

## Characteristic compact-reservoir boundaries

Frozen Euler `Boundary::Reservoir` retains a stationary conservative state in
the original gas owner's energy gauge. Its pressure and temperature supply
stagnation conditions for inflow; its pressure supplies the outflow backpressure.
A one-sided classical-wave compatibility relation determines the actual port.
A stationary two-sided ghost would add an artificial acoustic impedance and is
not used for this compact-volume model.

Incoming gas accelerates isentropically using the source's temperature-dependent
enthalpy and entropy. A safeguarded velocity solve matches the outgoing duct
wave; an admitted sonic state limits choked inflow. The energy flux carries the
reservoir's original total enthalpy exactly once. Outgoing gas follows the
thermally perfect shock Hugoniot or rarefaction invariant and carries its
resolved total enthalpy. The fixed interface samples departing supersonic waves
and sonic fans rather than overprescribing an outlet.

Scaled pressure-wave roots retain weak-shock/rarefaction increments. All caloric
differences come from the original source, without constant-gamma replacement.
Actual states and classical wave ordering must be admitted; cancellation and
bounded root work remain explicit. Unused stopped shocks or downstream fan
tails may lie outside the caloric interval without disqualifying an admitted
accelerating inlet or sonic outlet. A genuinely unsupported port state refuses.

`set_boundaries` updates prescribed end conditions atomically between steps.
Reservoir conditions are held fixed during each SSPRK2 step. The accepted
`left_boundary_change` and `right_boundary_change` are each signed into the
interval and retain the two-stage mass, momentum and energy exchange, alongside
the existing aggregate report. This supports driven inlets/outlets without
inferring exchange from a sampled pressure or flow value.

Seven focused native boundary regressions compare gamma-law shock/rarefaction
jumps and pressure derivatives (including 1e-10 pressure increments), analytic
sonic inlet/outlet limits, independent variable-Cp Hugoniot and quadrature,
reference-gauge invariance, source-window edge cases and driven-field balances.
All 28 frozen-dynamics tests passed with the integrated source. A separate
regression verifies that multicomponent flow refuses this frozen-state reservoir:
it requires a species-bearing state and energy-reference closure and must not
silently reinterpret this boundary as transmissive.

This estimated inviscid model assumes an ideal compact reservoir/entrance.
It supplies neither an empirical discharge coefficient nor finite neck inertia,
wall loss, multidimensional mixing, moving-volume work or an entropy theorem.

## Synchronized finite-plenum gas networks

`dynamics::network::GasNetwork1d` owns existing Euler branches, finite rigid
`Plenum` volumes, and explicit `PlenumConnection` endpoints. Every branch
retains its original mesh and field. All branches must share the exact frozen
source, composition, caloric interval, energy reference and accepted clock.
Duplicate, out-of-range and periodic endpoint connections refuse. Branch and
plenum counts are each capped at 4096; combined cell/plenum storage shares the
existing `MAX_EULER_CELLS` cap.

Each connected port evaluates the same physical characteristic reservoir
closure once at each stage. Its area-weighted mass and total-enthalpy exchange
feeds the branch and plenum with opposite signs. A compact plenum stores mass
and internal energy, without resolved bulk momentum; incoming kinetic energy
remains in that energy account. Its fixed support receives the reported port
impulse in the branch's axial coordinate and does no mechanical work.
Different branch axes are not added into a fictitious scalar network momentum.

One global SSPRK2 duration advances every owner. Actual branch face speeds and
the summed port area/speed divided by plenum volume limit both stages. A small
plenum can therefore limit the whole network timestep. Every retry restarts all
owners; state, clock and accepted exchange/report slices commit together only
after the final admission and cancellation check. Numerical scratch and report
storage are allocated at construction and reused.

Global reports integrate all branch and plenum mass/energy and count only
external-end exchange. Accepted per-port exchanges are the exact negation of
the corresponding branch-end report. Each branch also retains its full
momentum/area-force balance. Immutable access prevents independent branch
stepping while the network owns the common state.

Ten focused native network regressions cover original isolated-Euler identity,
tapered rest, compact-volume CFL, topology/source/clock admission, final-stage
cancellation and failure rollback, finite-pulse splitting, direct extensive
storage and port exchange reconstruction, orientation symmetry and bitwise
variable-Cp formation-gauge invariance. A 64/128-cell compact-plenum acoustic
check compares against `p'(t)/p'(0)=exp(-c*A*t/V)` before returning reflections;
normalized errors are 2.97668e-6 and 3.61846e-6, both within the fixture's
weak-wave tolerance. This is an independent impedance check, not an assertion
of an asymptotic refinement order.

The `fs-duct` `nonlinear_manifold` consumer uses original cylinder/cone
geometry and a 30 mL plenum. Runs with 16 and 64 cells per branch reach 4 ms
in 131 and 527 accepted steps. Both report zero mass residual; physical energy
residuals are -1.421085e-14 J and zero respectively. All temperatures remain in
the declared N2 200--1000 K source window, and both daughter passages respond.
The full staged pure-gas native run passed 58 tests and the existing eight
focused duct tests passed. This is not a workspace Cargo/DSR result.

A finite plenum is an explicitly compact, perfectly mixed storage model, not
a zero-volume lossless junction or a resolved three-dimensional shock/mixing
region. Species-changing/reacting junctions, valves/discharge coefficients,
wall heat, viscous loss, compliant volumes and entropy certificates remain
separate models. This contributes to the open broad gas-network bead without
claiming its complete reactive/MoC scope.

## Physical reflecting walls for evolving mixtures

Reflecting multicomponent boundaries solve the one-sided stopped shock or
isentropic rarefaction with the incident trace's current species composition.
They share the thermally perfect wall kernel with frozen Euler instead of
retaining the momentum component of a mirrored HLLC/HLL ghost flux. Every
species mass flux and the stationary adiabatic wall energy flux are exactly
zero; the physical wave speed enters the existing stage CFL calculation.

The borrowed local caloric law uses the admitted original species models,
temperature-dependent Cp/Cv, molar masses, energy/entropy differences and
source seams. Composition is fixed across each local inviscid wall wave but
may change between faces and accepted steps, including through the existing
reactive transport owner. Successful wall evaluation allocates no heap storage
and constructs no replacement gas/source owner. Source failures remain typed;
unsupported hot/cold wave paths and numerical budgets refuse without
extrapolation or clipping.

Five focused regressions compare arbitrary constant-Cp compositions against
analytic shock/rarefaction pressures, variable-Cp compression against independent
mass/momentum/energy jumps and entropy production, and expansion against a
separate Simpson integral. They also check orientation and independent species
formation-reference invariance, actual finite-volume wall wiring, and evolving
species/energy conservation with wall impulse accounting. In the retreating
manufactured fixture, the previous mirrored flux gives -0.722497216 Pa traction;
the physical wall pressure is +0.0110814482 Pa.

The complete integrated pure-gas native run passes 67 tests, including existing
reacting flow, characteristic reservoirs and finite-plenum networks. Eight
existing nonlinear duct tests pass against the refactored kernel, and the
64-cell-per-branch manifold still reaches 4 ms in 527 steps with zero reported
mass and physical-energy residuals. These are focused production-code checks,
not a workspace Cargo/DSR result or experimental validation. Reactive/catalytic
surface kinetics, wall heat transfer and moving-wall work remain separate models.


## Conservative molecular, thermal and viscous field evolution

`multicomponent::dissipation::DissipationOperator<L>` binds the original gas,
stationary constant-area mesh, boundary conditions and explicit transport law.
`advance_for` completes an exact requested duration through bounded SSP RK2
parabolic substeps. It publishes the complete field atomically and leaves its
clock unchanged for the composing integrator. Nonuniform spacing is admitted
with the actual distances and physical cell volumes; variable-area viscous
geometry requires additional stress terms and is refused.

Every shared face carries barycentric Maxwell-Stefan mass fluxes Ji, momentum
flux `-tau`, and gauged-energy flux
`sum_i (h_i/M_i - e_ref_i) Ji - lambda grad(T) - u*tau`.
The longitudinal stress is `tau=(4mu/3+mu_bulk) grad(u)`; bulk viscosity is an
explicit authored value, with zero selecting Stokes' hypothesis. The source
enthalpy difference is evaluated as `[Delta u_i(Tmin,T)+RT]/M_i`, avoiding
large-formation-reference subtraction. No independent viscous-heating source
is added to total energy. Viscous loss of kinetic energy therefore appears as
internal heating through the conservative update.

Periodic faces are shared bit-for-bit. Transmissive ends use zero diffusive
gradients. Stationary reflecting ends are impermeable and adiabatic, retain
viscous traction, and perform exactly zero wall energy work. Frozen reservoir
states have no species axis and are refused by the mixture boundary owner.
Parabolic estimates are checked at both stages; positivity, source support,
work limits and cancellation use bounded refusal/retry, with no clipping of
conserved states. Complete species, momentum, physical-energy and boundary
balances are returned.

`TransportLaw` accepts actual face T,p and composition. `ConstantTransport`
retains explicit coefficients and validity for authored approximations and
analytic cases. `KineticTransport` adapts the source-bound Lennard-Jones,
Chapman-Enskog/Wilke/Eucken model with reusable property storage. Matrix,
face, stage and candidate buffers are prepared once. Derived fraction/gradient
simplex identities use a checked roundoff-sized remainder on the largest
component; no conserved species value or physical flux is projected.

Seven focused native tests pass: Fourier heat-mode decay/refinement; analytic
binary diffusion; viscous velocity decay and kinetic-to-internal conversion;
unequal-mass partial-enthalpy transport and bit-identical evolution after a
1e18 formation-reference shift; nonuniform geometry; wall traction with zero
energy work; and cancellation, work limits, zero duration and reuse. Fourier
L1 errors at 40/80 cells are 1.91446548e-5 / 4.77807730e-6. This remains an
Estimated continuum operator: the parabolic estimate is not a spectral or
fully discrete entropy/invariant-domain proof. Soret/Dufour effects, nonideal
phases, Knudsen transport and variable-area viscous stresses remain separate
closures.


## Reacting Navier-Stokes composition with original-source molecular transport

`ReactingEuler1d::step_with_dissipation` composes exact-duration molecular and
chemistry half steps around conservative advection:
`D(dt/2), C(dt/2), A(dt), C(dt/2), D(dt/2)`. Each D stage may subcycle its
parabolic update. Only A advances the clock. If heat release changes acoustic
stability and A accepts a shorter duration, all five stages restart from the
accepted field at that duration. Cancellation, source refusal and exhausted
retries cannot publish a partly reacted or diffused field.

The final D endpoint is checked through the original caloric owner against
the reactor's declared temperature interval before publication. Viscous
heating or species enthalpy transport may otherwise leave a narrower chemistry
interval even while the molecular/flow source still admits the state; this
condition triggers bounded complete-step halving without temperature clipping.

The returned optional `SplitDissipationReport` retains both accepted molecular
halves. Whole-step species, momentum, mass and physical total-energy residuals
include both molecular boundary transfers, advection/area terms, and chemical
species/reference-energy changes. No extra chemical or viscous heat source
is added to conservative total energy. The existing `step` method retains
the chemistry/advection/chemistry path with no molecular stages.

Five coupled native tests pass: analytic advection/diffusion/isomerization with
mesh refinement; heat-release-driven duration restarts; closed-wall molecular
traction and zero wall energy work; late cancellation and reused workspaces;
and viscous heating against a narrower chemistry interval, including atomic
refusal and successful shorter replay. The analytic concentration L1 errors
are 6.60993650e-5 / 2.03342713e-5 at 40/80 cells. The complete isolated native
source-linked suites pass 132 thermochemistry and 79 gas unit tests; these
results are scoped to those owner builds, not the full runtime dependency graph.

`cargo run -p fs-gas --features experimental-euler --example reacting_oxygen -- 40 molecular`
enables the complete original-source path in the O/O2 example. Calorics,
reaction and Lennard-Jones parameters retain Cantera v3.2.0 GRI3 attribution
and original blob 9a62004ec17bbec8ac11e97785932d426dd58610. At 40 cells over
1 ms the native run accepted 113 coupled steps / 226 molecular substeps,
with mass residual -4.336809e-19 kg, physical total-energy residual
-6.274981e-13 J, and zero reported momentum residual. Atomic oxygen mass
falls from 4.237578e-5 to 3.380238e-5 kg. These are numerical consistency
checks of an isolated reaction, not an experimentally validated air mechanism.

The transport law's dilute-neutral-gas, Eucken-conductivity, explicit bulk
viscosity and constant-area restrictions remain in force. This symmetric
composition and either selected chemistry method do not by themselves
establish global second-order accuracy or a fully discrete entropy inequality.

## Characteristic ports with distinct local compositions

`MixtureGas::reservoir_flux` connects a conservative duct trace to a stationary
species-bearing reservoir in the same original source bank and per-species
energy gauges. Local compositions may differ. One shared borrowed-caloric
solver now supplies both frozen and mixture ports: the duct composition owns
the outgoing pressure wave, while the reservoir composition owns isentropic
acceleration and the sonic inlet limit. The advected donor supplies every
species flux and its total enthalpy; structural species zeros remain zero.

The public output buffer is updated only after the complete port has passed
source, numerical and cancellation checks. Successful evaluation uses bounded
stack scratch and no replacement gas owner or heap workspace. A strong inflow
whose auxiliary velocity bound reaches the full cooling interval locates its
supported sonic throat before constructing an unused cold endpoint; actual
port admission remains strict.

Five independent regressions cover two different gamma laws across a port,
analytic choked inflow and sonic outflow, branch acoustic impedance, donor
enthalpy, original formation-reference/orientation invariance and unchanged
output on late cancellation. All 72 focused staged gas tests passed, including
the existing frozen networks, wall reflections and reacting-flow checks.
This run precedes the concurrently added molecular-dissipation module and is
not a full workspace Cargo/DSR verdict. The closure remains an inviscid compact
reservoir model, without resolved junction mixing or empirical discharge loss.

## Conservative networks with evolving species

`multicomponent::network::MixtureNetwork1d` connects the existing conservative
mixture branches to finite rigid `Plenum` volumes. Every owner shares the exact
original species bank, canonical order, caloric path and individual energy
references; its local composition and temperature may differ. The shared
characteristic port supplies donor species, pressure impulse and total
enthalpy once. The corresponding branch-end and plenum accepted exchanges
are exact opposites for each species and energy.

One global SSPRK2 transaction advances all branches and plenums. Actual port
speeds constrain branch and compact-volume CFL at both stages; every retry
restarts all owners. State, clocks and reports publish together after final
admission/cancellation. Flat numerical/report storage is prepared once and
reused, and compact-volume temperature inversions use cached guesses. Compact
storage has no bulk momentum; its fixed support receives the reported local
port impulse while incoming kinetic energy remains in the energy balance.

Four focused regressions cover different-composition rest in tapered geometry,
finite-pulse splitting and actual species transport, independent extensive
species/physical-energy accounting, accepted port exchange, compact-volume CFL,
source/topology refusal, cancellation/reuse, orientation and formation gauges.
The integrated native gas build, including the concurrently landed molecular
and reacting transport modules, passes 88 tests.

`cargo run -p fs-gas --features experimental-euler --example mixture_manifold -- 64`
runs pressurized N2 into an air surrogate containing N2/O2 through three passages
and a 30 mL plenum. It retains the original 200--1000 K NASA-9 coefficients from
Cantera v3.2.0 `data/airNASA9.yaml`, blob
`5975d1303819332c1bd77a783d9873b66f65c952`, with the declared 1 atm reference.
The 16/64-cell runs reach 4 ms in 149/593 steps, with no HLL fallbacks or
halvings. Physical-energy residuals are -3.170523e-14 / 1.919073e-15 J;
the largest species mass residual is 1.084202e-19 / 1.355253e-20 kg.
Both daughters receive nitrogen and the final plenum N2 mass fractions are
0.973409707 / 0.979168255.

All sampled states remain inside physical species bounds and the source path.
A small composition-extremum undershoot remains in this conservative MUSCL
baseline (64-cell daughter minimum N2 fraction about 0.766145 versus initial
air 0.767082); no composition maximum principle or mesh-converged prediction
is claimed. This network transport owner has no chemical source, resolved
junction mixing, molecular port diffusion, wall loss or entropy certificate.
The required `dsr quality --tool frankensim` could not run because DSR is
unavailable; focused native checks are not the workspace gate.

## Finite-rate reacting shockwave networks

`multicomponent::network::reactive::ReactingNetwork1d` owns a complete mixture
network and one prepared original-source isochoric reactor. Every branch cell
and finite rigid plenum participates in `C(dt/2), A(dt), C(dt/2)`, using its
actual volume. The existing reacting single-branch owner and this network
share the same per-cell chemistry update. Momentum and physical chemical total
energy are fixed during C; the stored referenced-energy change is exactly the
negative change in species reference energy, without a second heat-release
source. Endpoint sensible energy is independently reconstructed from the
original species models.

The initial duration respects branch and compact-volume acoustic scales. If
the first chemical half raises wave speeds and transport accepts less time,
the whole composition restarts from accepted branch and plenum states at the
actual duration. Exhausted retries, source refusal and cancellation cannot
publish a partly reacted network. Fields, clocks, combined reports and port
exchanges publish together after final admission. Trial fields, amount storage,
reactor and numerical workspaces are reused; reset work polls in bounded tiles.

Combined reports distinguish chemical species/reference-energy changes from
external transport and reconstruct physical total-energy balance. Accepted
port exchanges remain opposite to the branch transport reports. Those reports
explicitly describe the intermediate transport stage; releasing the raw network
clears its stale transport-only whole-step report. Compact plenums retain the
same perfectly mixed, zero-bulk-momentum, fixed-support interpretation.

Four new regressions compare homogeneous reversible isomerization with analytic
decay, and strong exothermic steps with independent constant-volume reactors at
the accepted duration after a heat-driven restart. They separately reconstruct
branch/plenum chemistry and actual port-energy budgets, reject mismatched
original sources and plenum temperatures outside the chemistry interval, and
verify late cancellation plus workspace reuse. The complete integrated native
gas suite passes 92 tests, including existing reacting molecular-transport
checks. Eight existing nonlinear duct tests pass against the shared port and
mixture-network kernel. These are production-source native checks, not a full
workspace Cargo/DSR result; DSR is unavailable in this environment.

This owner composes inviscid network transport with finite-rate chemistry.
Molecular junction diffusion, wall heat/loss, valves and resolved mixing remain
separate closures. The chemistry entropy diagnostic is not a network entropy
certificate, and the selected chemistry method combined with symmetric
splitting does not establish global second-order convergence.


## Two-layer compartments with conservative moving-interface work

`multicomponent::compartment::TwoZone` adds a reusable two-zone control volume
for issue #2. It retains the original `MixtureGas` species bank and source-path
admission. Each layer stores canonical species masses, sensible internal energy
in that bank's gauge, physical volume and temperature. A layer may initially
be absent exactly; real entrainment and heat can create it without artificial
mass, species or temperature seeds.

A source increment carries heat and advected partial enthalpies. The fixed room
volume is partitioned by common pressure and the actual variable-Cp ideal-gas
EOS. The implicit midpoint interface work is
`W = (p_old+p_new)/2 * (V_lower_new-V_lower_old)`; the other layer receives
exactly `-W`. A bounded scalar energy-partition solve recovers both temperatures
through the original gas inversion. Its derivative is positive on the admitted
positive-Cv source paths. Neither fixed gamma nor a new caloric surrogate is
introduced. The total sensible-energy source is conserved by construction;
physical-energy reports retain the original species formation references.

`inflow` and `carried_rate` provide canonical enthalpy-bearing source rates.
The same transfer is added to the receiver and negated for the donor.
`advance` completes exactly one supplied-source duration atomically; a network
or time-quadrature owner must recompute its rates when retrying. Prepared rooms
can copy accepted physical state only when their source and geometry agree.
Cancellation, negative inventory, unsupported calorics, an out-of-room interface,
or a bounded solve failure publishes nothing and leaves the clock unchanged.

Seven native tests pass: the constant-Cp CFAST pressure reduction and adiabatic
interface refinement; independently integrated variable-Cp heating; conservative
unequal-species enthalpy exchange; formation of an initially absent upper layer;
original EOS, composition and interface checks; cancellation/refusal/reuse;
and representability admission including nonbinary room volumes. The governing
mass, energy, ideal-gas and fixed-volume equations are documented in
[NIST TN1889v1, chapter2](https://doi.org/10.6028/NIST.TN.1889v1).
This is an Estimated well-mixed low-Mach compartment owner, not a flame,
ceiling-jet, acoustic, soot-radiation or experimental-validation claim. Prescribed
source rates are frozen over a step; the midpoint work alone does not make the
complete source/network integration second order.


## Inventory-limited compartment burning

`multicomponent::compartment::burning` binds an original admitted ideal-gas
kinetics mechanism and one conserved global reaction to the exact room species,
NASA calorics, and per-species energy references. The requested nonnegative
forward extent is capped by every net reactant inventory, including oxygen,
with optional nonnegative incoming species rates over the explicit interval.
Incoming rates are availability only: the caller must add their physical mass
and enthalpy exactly once. A driver with simultaneous reactant outflow must
reserve that withdrawal first or transport before burning; timestep halving
alone cannot make an additional withdrawal compatible with complete depletion.

The adapter returns chemistry-only species rates and the energy-reference
correction `-sum(dm_s/dt * e_ref_s)`. Original formation energies therefore
produce sensible heating while preserving physical internal energy. Reported
reaction enthalpy and power are diagnostics, never a second heat source.
Power-to-extent conversion uses the original reaction enthalpy at an explicit
reference temperature. One common extent may be rounded down by a bounded
number of adjacent floating-point values to admit actual update arithmetic;
individual species are never clipped, seeded, or independently repaired.

This is an estimated global, inventory-limited burning closure, not an
Arrhenius ignition/extinction model, equilibrium solver, soot model, or proof of
second-order temporal accuracy. Reactions with species on both sides are
explicitly refused by this single-pass adapter. Six focused tests cover oxygen
exhaustion, atom/mass conservation, analytic variable-Cv formation-energy
heating, inflow counted exactly once, a one-ULP depletion case, temperature-
dependent reaction enthalpy, cancellation, and source identity/refusals.


### Clock-preserving transport/burning composition

`InventoryLimitedReaction::apply_split` evaluates the actual transported
inventory with no repeated incoming source, delegates chemical species/energy
and pressure-interface relaxation to the original atomic TwoZone advance, and
then restores the physical clock infallibly. The source duration controls the
reaction increment but does not advance elapsed time twice. This first-order
split allows prior ventilation withdrawals to reduce oxygen availability before
burning. The returned chemical source is recorded once by a composing network;
reported reaction heat is still never an extra energy source. All eight
original burning tests pass, including two added regressions for post-outflow
oxygen depletion, fixed clock, paired interface work, physical-energy balance,
late cancellation, reuse, and caloric-domain refusal of excessive source heating.

## Molecular transport through finite-plenum networks

`multicomponent::network::dissipation::NetworkDissipationOperator` extends the
existing original-source Maxwell--Stefan, Fourier and longitudinal-viscous
operator across branch/plenum connections. Every branch retains the existing
stationary constant-area restriction. A caller supplies one strictly positive
plenum-side `Length` per connection. Together with the branch half-cell width
it defines the constitutive gradient and interpolated face state; no mixing
length is inferred from plenum volume.

Each connected face uses the same molecular kernel as an interior branch face.
Species diffusion is mass-averaged, carries each original partial enthalpy,
and shares total-energy flux with conduction and viscous work. The branch and
plenum receive equal-and-opposite species and energy exchange. A stationary
plenum retains zero bulk momentum; its fixed support receives the signed local
axial impulse. No scalar momentum is summed over differently oriented branches,
and no separate viscous-heating source duplicates conservative energy.

One global explicit SSPRK2 subcycle advances all participating physical
volumes. Its estimated parabolic bound includes both branch cells and the sum
of all port rates divided by each plenum volume. Every requested molecular
interval completes its exact duration before fields publish; the flow clock
is unchanged. Law, constitutive solver, stages and diagnostic arrays are prepared
once and reused. Candidate port/branch reports remain separate from accepted
reports, allowing a composing owner to reject an entire split step atomically.

`multicomponent::network::dissipating::DissipatingNetwork1d` owns a complete
`D(dt/2), A(dt), D(dt/2)` step. If advection shortens the duration, both molecular
halves restart from the original accepted state. Its report closes species and
physical energy against the combined external boundary transfers, and its
accepted molecular-port slice is the sum of both D halves. The A-only exchange
and branch reports explicitly describe the intermediate advection stage.
Releasing the raw network clears its stale intermediate whole-step report.

Five molecular regressions check two-volume binary/Fourier exponential
relaxation and timestep refinement, source-consistent partial-enthalpy transfer,
viscous support impulses, formation-gauge/orientation invariance, compact-volume
subcycling and atomic cancellation/reuse. Three composing-owner regressions
check zero-law identity, independent A+D port/storage budgets, late cancellation
and exhausted molecular work without partial publication. All eight new checks
passed in the integrated production-source native build.

The compact gradient-length closure is an Estimated continuum model. It does
not resolve multidimensional mixing, radial wall layers, variable-area viscous
stresses or rarefied gas. Soret/Dufour effects and a fully discrete entropy or
invariant-domain theorem remain outside this operator. Native source-linked
checks are supplementary: `dsr quality --tool frankensim` is unavailable here
(`dsr: command not found`), so no workspace DSR verdict is claimed.

## Reacting networks with molecular transport

`ReactingNetwork1d::step_with_dissipation` composes
`D(dt/2), C(dt/2), A(dt), C(dt/2), D(dt/2)` for every branch cell and finite
plenum. It reuses the same prepared molecular network operator and original
constant-volume chemistry update; each molecular half may subcycle, while
only A advances time. If heat release or diffusion changes the admissible
acoustic duration, all five stages restart from accepted physical storage at
the actual duration.

Molecular operator reports remain unpublished through coupled retries and
cancellation. Complete fields, clocks, reaction/transport/molecular reports and
the sum of both molecular port transfers publish only after the final check.
The original C/A/C method remains available without molecular stages. Combined
balances include both external molecular transfers and original chemical
species/reference-energy changes; no duplicate heat-release or viscous-heating
source is introduced. Final molecular temperatures in every branch cell and
plenum must remain inside the mechanism's possibly narrower caloric interval.

Three new coupled regressions pass: a two-pool analytic reversible
reaction-plus-diffusion exponential, heat-release-driven complete-duration
restart with homogeneous chemistry identity, and late cancellation preserving
both network and separately owned molecular reports followed by successful
reuse. All 110 staged native gas tests passed, including the earlier
single-branch molecular/reacting paths and the original compartment tests,
before the subsequent Maxwell--Stefan trace-gradient scaling correction.
After that correction, all ten exact-production-module Maxwell--Stefan tests
and rebuilt 16/64-cell molecular manifold consumers passed. This scoped
validation does not claim tests added concurrently elsewhere on main or the
unavailable workspace DSR gate.

The selected chemistry method and symmetric composition do not by themselves
establish global second-order accuracy or a discrete entropy certificate. The
molecular operator's constant-area, continuum, explicit mixing-length and
transport-law applicability limits remain unchanged.

## Molecular N2/O2 shockwave manifold consumer

`cargo run -p fs-gas --features experimental-euler --example mixture_manifold -- 64 molecular`
evolves the pressure/composition pulse with `DissipatingNetwork1d<KineticTransport>`.
The existing default remains inviscid and retains its conical daughter. Molecular
mode explicitly selects three stationary cylinders of lengths 0.6/0.5/0.7 m and
radii 0.010/0.008/0.008 m, a 30 mL rigid plenum, and a declared 0.010 m plenum-side
mixing length at each port; each gradient includes the endpoint half-cell width.

Calorics, masses and reference energies remain the original airNASA9 source bank.
The dilute nonpolar kinetic law uses original Cantera v3.2.0 `data/gri30.yaml`
collision data (blob `9a62004ec17bbec8ac11e97785932d426dd58610`): N2 sigma
3.621 angstrom and epsilon/k 97.53 K; O2 sigma 3.458 angstrom and epsilon/k 107.4 K.
It declares Wilke viscosity, estimated Eucken conductivity, zero bulk viscosity,
200--1000 K and 1000--1e6 Pa. This does not resolve radial wall friction/heat
transfer, dense gas, Soret/Dufour or rotational relaxation, and the compact
mixing distance is a constitutive input rather than a resolved plenum flow.

Fresh native runs with Maxwell--Stefan evaluator version 2 completed 4 ms:

| Cells per branch | Wave steps | Molecular substeps | Maximum species residual (kg) | Physical-energy residual (J) |
| --- | --- | --- | --- | --- |
| 16 | 168 | 336 | 1.36e-20 | +5.19e-15 |
| 64 | 669 | 1338 | 1.09e-19 | +1.43e-14 |

Both runs had zero rejected molecular substeps, coupling retries, HLL fallbacks
and timestep halvings. Every sampled field was finite, with positive pressure,
source-admitted temperature, and species fractions in [0,1] summing to one at
printed precision. The sampled minimum N2 fraction was 0.765013831926 at 16 cells
and 0.766338161781 at 64 cells, below initial air's 0.76708: these are not
composition maximum-principle, entropy, convergence-rate or experimental
accuracy certificates. Both daughters receive pressure/composition waves and
nonzero conservative molecular port transfers. The rebuilt default 16-cell
inviscid run retained its prior 149 steps and -3.17e-14 J physical-energy residual.

## Stationary wall drag and finite wall heat storage

`multicomponent::wall_exchange` applies a source interval to one original
mixture cell without changing species or the flow clock. Linear drag uses its
exact frozen-rate exponential; Darcy drag uses its exact frozen-factor
quadratic decay with an explicit hydraulic diameter. The gas retains total
energy while lost kinetic energy becomes sensible energy, and the stationary
support receives the opposite signed axial impulse. No second viscous-heating
term is added to total energy.

Thermal contacts may be adiabatic, a prescribed-temperature reservoir, or a
finite lumped wall with explicit heat capacity in J/K. A contact's conductance
is the complete exposed h*A in W/K. A safeguarded backward-Euler temperature
solve evaluates the original mixture's NASA-9 energy and Cv rather than a
frozen-gamma substitute. Finite-wall storage is eliminated in the scalar solve
and its actual represented temperature change is included in the energy audit.
Reservoir heat is external; gas/finite-wall exchange is internal. The step
refuses source-domain or residual-budget violations rather than clipping gas
temperature or creating heat.

`network::wall_exchange::NetworkWallOperator` binds one distinct contact to
every branch cell and finite plenum. It stages all momentum/energy and wall
temperature updates in prepared buffers, preserves each species exactly, and
publishes the entire source interval only after its final cancellation check.
Failure leaves accepted gas, walls, clocks and reports intact. Its accepted
support impulses remain per branch axis. A successful standalone source update
clears an old advection report that no longer describes the changed fields.
The source operator accepts the existing stationary mesh geometry; the separate
molecular operator retains its constant-area restriction.

The network has a shared finite-wall trial workspace for complete coupled
retries. Candidate walls start from accepted storage once per trial, retain
both wall half steps, and publish only when the composing owner accepts every
stage. The report separates heat into gas, actual finite-wall energy change,
external reservoir heat, kinetic energy thermalized, and numerical residuals.

These are Estimated lumped/averaged closures with explicitly supplied
coefficients. They do not resolve developing, turbulent or unsteady shock
boundary layers; solid conduction between contacts; shared multi-cell walls;
radiation; moving-wall work; or a thermodynamic/entropy theorem. Backward-Euler
heat integration is first order even when surrounded by symmetric source
composition. No throughput or experimental-accuracy claim is made.

Focused Rust regressions accompany the source/network implementation for an
independent two-body heat solution and timestep refinement, exact linear/Darcy
drag and support impulses, reservoir exchange/source refusal, and late
cancellation followed by reuse. Native compilation and execution are blocked in
this session: the execution environment has no usable connection. Consequently
these new regressions have not been run, and no DSR, Cargo, native-suite or
physical-consumer pass is claimed for this extension.

## Wall sources composed with molecular and reacting networks

`DissipatingNetwork1d::step_with_wall` performs W/2, D/2, A, D/2, W/2.
`ReactingNetwork1d::step_with_wall` performs W/2, C/2, A, C/2, W/2;
`step_with_dissipation_and_wall` surrounds the existing D/2, C/2, A, C/2, D/2
composition with both wall halves. The existing methods continue without W.
Every half covers the actual accepted wave duration. A shorter advection step
or a refused source stage resets gas fields and finite wall storage before
retrying the complete composition.

Both wall heat transfers enter the gas boundary-energy balance. The wall report
separately retains actual finite-wall storage changes and external reservoir
heat, so a gas-plus-wall audit can distinguish internal exchange from external
energy supply. Stationary-wall kinetic thermalization is diagnostic only and
does not enter the total-energy source a second time. Final branch and plenum
temperatures must still fit the chemistry mechanism's admitted interval after
the final wall half. Gas, molecular reports, wall temperatures and both wall
halves' aggregate support impulses publish only after the final cancellation
check.

The manifold example accepts `molecular-wall`:

```bash
cargo run -p fs-gas --features experimental-euler --example mixture_manifold -- 64 molecular-wall
```

It retains the original N2/O2 caloric and kinetic sources, and adds explicitly
illustrative fixed inputs: Darcy factor 0.02, h=50 W/(m² K), wall areal capacity
2000 J/(m² K), initial wall temperature 300 K, cylindrical sidewall areas and an
authored 0.005 m² plenum wall area. It independently integrates final retained
wall temperatures to report a complete gas-plus-wall energy balance alongside
gas-only and molecular-port budgets. These coefficients are not a named
material calibration or an unsteady shock-layer model.

Focused coupled regressions exercise closed gas-plus-wall storage, zero-wall
identity, exact half-duration accounting and late cancellation/reuse. Additional
cases check an independently integrated nonconstant mixture Cv polynomial and
simultaneous reaction, molecular transport and finite-wall energy storage.
These new Rust tests and the new example mode are committed but have not been compiled
or executed in this session: `dsr quality --tool frankensim` could not launch
because the execution environment is unavailable. No native/DSR pass, measured
wall attenuation, global second-order accuracy or full entropy certificate is
claimed.

## Finite-rate chemistry in common-pressure gas layers

`multicomponent::compartment::kinetics::ZoneKinetics` binds the complete original
gas/kinetic source and a prepared room geometry. Each selected present layer
uses the existing `IsochoricReactor` at its actual volume and physical internal
energy, including formation terms. The resulting time-integrated species
changes become canonical `LayerRates`; their only energy term converts the
original numerical species references. No second combustion-heat source is
added. Both layers then use the original simultaneous pressure/interface-work
solve. The enclosing transport owner advances the physical clock once.

The prepared reactor and room candidate are reused. Actual rate reconstruction,
mass and physical-energy residuals are checked before publication. Gas/source
and geometry mismatches, exhausted work, domain failures and late cancellation
leave the accepted room unchanged. The kinetic temperature interval is checked
again after mechanical relaxation: compression can invalidate a layer even
when its homogeneous reactor endpoint was admissible. A shorter retry starts
from the supplied accepted state, independent of previous workspace contents.

The chemistry report retains original per-layer entropy, conservation, local
error and adaptive-work diagnostics. The selected fixed-volume chemistry
method is followed by first-order mechanical/source composition; the homogeneous entropy check is not a
fully discrete theorem for the moving-interface or ventilated system. Supplied
mechanisms retain their own source and validity limitations. No ignition,
flame-structure, CO/soot-yield or experimental accuracy claim is inferred.

Eight new regressions cover analytic finite-time reversible relaxation,
unequal reacting layers with variable heat capacities and actual interface
work, disabled/absent layers, trial replay, formation-reference invariance,
source/domain/work refusals, post-relaxation domain rejection and late
cancellation. All 23 focused native compartment tests passed; that is separate
from the unavailable full-workspace DSR gate.

## Reduced finite-rate Newton system and independent hydrogen reactors

The original extent integrator selects the smaller of the full reaction
matrix and an equivalent species-plus-temperature Schur system. With ten
species and 29 reactions it factors 11 rather than 29 unknowns. The additional
thermal coordinate retains the original compensated temperature derivative
with respect to each reaction extent, preserving elemental energy-reference
invariance. Every reaction extent is reconstructed, and bounded residual
correction checks the solution against the original full linear equations.
The original pivoting, positive-amount, caloric-domain, nonlinear, energy,
element/charge, entropy, local-error, retry and cancellation checks still apply.
No measured wall-clock speedup or higher integration order is implied.

Five focused reduced-system cases compare independent full/reduced equations,
stiff cancellation recovery, exothermic conservation, elemental datum shifts,
and deterministic workspace reuse. The original gas library passes 151 tests.
A further independent Cantera 3.2.0 test uses the complete sourced hydrogen
mechanism, continuous calorics, no radical seed, 1 litre at 101325 Pa and
H2:O2:N2 = 1:1:3.76. Initial 1000 K advanced 1 ms and 1200 K advanced 0.1 ms
agree with Cantera within 3 K, 0.2% pressure and 0.5% + 1e-12 mol per species;
physical internal energy agrees within 2e-8 J and atom/mass inventories close.
For the 1200 K case, tightening local relative tolerance from 2e-4 to 2e-6
reduces the combined temperature/species endpoint error by over 70%.
`data/reference/thermochem/hydrogen_cantera32.py` reproduces the independent
reference using a separately installed, pinned development-only Cantera.
These are numerical source/integration comparisons, not experimental ignition
validation, a fully discrete coupled entropy theorem, or a second-order claim.

## Optional second-order L-stable homogeneous chemistry

`IsochoricKineticsSolve::method` explicitly selects
`IsochoricKineticsMethod::BackwardEuler` (unchanged default) or `Sdirk2`.
The latter implements Alexander's two-stage stiffly accurate method with
gamma = 1 - 1/sqrt(2), A = [[gamma, 0], [1-gamma, gamma]], b = A[1].
Reference: R. Alexander, SIAM J. Numer. Anal. 14(6), 1006–1021 (1977),
https://doi.org/10.1137/0714068. This is a second-order L-stable method for the
smooth homogeneous reactor ODE, not an order claim for any surrounding split
transport, radiation, pyrolysis or pressure-interface composition.

Both stages reuse the original constrained extent Newton owner, source fixed
U,V, smaller factorization, caloric inversion, positive-amount and work gates.
The second stage keeps its known reaction forcing separate from the physical
baseline and starts Newton from the admitted first-stage state; an explicit
shifted baseline may be nonphysical and is never used as a state. Each stage
must pass an entropy comparison with its physical step baseline. Stage two is
not incorrectly required to exceed stage one's entropy. These are actual
acceptance gates, not unconditional Runge–Kutta positivity/entropy theorems.

One reaction-sized forcing vector is reserved only for SDIRK2 at construction;
advances allocate no new stage storage. The two-half-step result is retained
without Richardson extrapolation, with the error difference divided by three
and a deterministic cube-root step controller. Backward-Euler behavior and its
square-root controller remain the default.

Four new native tests demonstrate fixed-step second-order convergence against
analytic reversible relaxation and its independent stability function, stiff
damping, physical energy/elemental datum invariance, cancellation with replay,
and two independent Cantera 3.2.0 hydrogen endpoints. All 156 native gas tests
pass. At identical local relative tolerance 1e-5, the 1200 K/0.1 ms hydrogen
fixture uses 322 accepted SDIRK2 steps and 3714 Newton iterations versus
4485 and 19983 for backward Euler; temperature errors against Cantera are
0.01556 K and 0.04168 K respectively. These are fixture-specific work counts,
not a universal speed or accuracy guarantee.

## Shared chemical sources for axial and planar conserved cells

`multicomponent::reactive::CellChemistry` binds the identical original
`MixtureGas` species bank and `IdealGasKineticsV1` mechanism, retaining the
reactor's selected time method and admitted source interval. Its axial
`[rho_i..., px, E_ref]` and planar `[rho_i..., px, py, E_ref]` updates share
the same conversion, original fixed-U,V reactor and endpoint reconstruction
as existing one-dimensional and network chemistry. Both planar kinetic terms
enter internal-energy admission; both momenta are preserved. Reaction changes
the stored reference gauge by the negative original species-reference change.
There is no additional heat-release source on physical chemical total energy.

The conservative endpoint is independently reconstructed at the original
reactor temperature. Caloric precision and the kinetic-energy consequence of
any molecular-mass residual must fit the declared energy tolerance. Each cell
and its accumulated report publish together only after the final cancellation
poll. Source, energy, work-counter, dimension and cancellation refusals leave
both unchanged. Numerical scratch is reusable after a refusal. Successful
advances allocate no heap storage.

The prepared reactor and cell source expose actual numerical heap capacities,
including original kinetics, constrained Newton, stage and amount buffers.
Immutable mechanisms, gas caloric preparation, fixed-size owners and bounded
stack scratch are explicitly outside that numerical budget. Higher-level
spatial owners can add these capacities to their own field workspaces before
reserving a second complete field.

Four focused regressions compare exothermic planar cells with the original
homogeneous reactor for backward Euler and SDIRK2, check the exact axial limit,
rotation and energy-datum invariance, and exercise late cancellation/replay,
source/path mismatch and complete cell/report rollback. This chemical-source
adapter alone supplies no spatial transport, molecular diffusion, ignition or
flame validation, and no global splitting-order or entropy certificate.
