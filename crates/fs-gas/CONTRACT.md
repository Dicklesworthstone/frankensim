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

The solver advances reaction extents with damped implicit Euler and adaptive
step doubling. Every candidate reconstructs amounts from the original feed
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
the two-half-step result is accepted. This is a first-order, L-stable implicit
base method with an estimated local error, not a global accuracy certificate.

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
