# fs-gas contract

## Constant-pressure finite-rate chemistry

`chemistry::isobaric::IsobaricReactor` advances a sealed ideal-gas parcel at
prescribed pressure with total enthalpy conserved. It shares the original
extent integrator, backward Euler/SDIRK2 stages, adaptive comparisons, reduced
Newton factorization, source admission and conservation gates. The rigid U,V
entry points retain their original behavior. The prepared adapter reuses its
admitted storage across independent parcels; owning convenience results retain
the original mechanism and endpoint thermodynamic receipts.

Every trial computes its physical volume from nRT/p and evaluates kinetics at
the resulting actual concentrations. Reference-volume extent coordinates keep
stoichiometry exact; their rate Jacobian includes both volume and concentration
dilution derivatives. Temperature recovery uses original source enthalpy
increments, coefficient-combined reaction Delta h and Cp. Chemical formation
energy enters that single enthalpy constraint, never an additional heat source.
Entropy and elemental, charge and mass balances remain extensive at the
physical variable volume. Results expose initial/final V, H, U and outward
p Delta V, so the adiabatic first law can be checked independently.

The explicit-H entry point first resolves the initial thermodynamic state at
the authored composition. That resolved state defines initial volume and zero
boundary work; the temperature hint is not an additional energy source.
Unresolved energy targets, source domains, invalid pressure, cancellation or
bounded work exhaustion publish no accepted result. This is a homogeneous,
single-temperature ideal-gas reactor with prescribed constant pressure, not a
piston inertia, heat-transfer, nonideal or multiphase model.

Focused G1/G3 fixtures cover analytic isothermal dissociation, concentration-
sensitive bimolecular recombination, both stiff methods, full constrained
Jacobian differences, original reduced Newton work, actual binding-energy
changes and exact elemental energy-reference shifts. G4 exercises prepared
cancellation/reuse and zero-duration explicit-H semantics.


`fs-gas` owns L3 evolving gas systems. Constitutive and species laws remain in
`fs-thermochem`; this crate does not replace those sources with an independent
caloric model. Production dependencies remain pure Rust and Franken-owned.

## Material wall conduction in gas networks

With `experimental-euler,material`, `NetworkWallOperator::with_material_conduction`
binds original temperature or phase-dependent conductivity sources and explicit
planar length/area to each selected branch or plenum thermal contact. The gas
film and solid resistances add in series. Existing source-bound circular gas
transport still supplies the film from the actual gas mixture, density, T and p;
other contacts retain their originally authored film. The operator preserves
that separate film so repeated steps do not compound the same solid resistance.

In the default `WallConductionEvaluation::InitialState` mode, material k is
evaluated at the actual starting wall state of each source half, then frozen
for the original conservative implicit heat solve. The final
wall state must remain inside the original source domain. Each next half reads
the evolved wall temperature, including the existing reaction/wave composition.
Finite wall storage remains exclusively with the original contact and optional
material caloric binding; conductivity adds no second heat capacity or heat
supply. Any source refusal or late cancellation preserves all accepted gas,
wall temperatures and reports. Source lines, card identity and geometry remain
inspectable. The prepared temperature-only source refreshes without database
queries or allocation; phase-source refresh uses the original scalar resolver
with its selected claim pinned, retaining each query's evidence and conditions.

`MaterialWallConduction::with_evaluation(ImplicitEndpoint)` instead resolves
the original conductivity at every proposed caloric wall endpoint inside the
same bounded scalar heat equation. It requires the existing material heat-storage
binding and supports both temperature-dependent and phase-dependent sources.
The retained gas film remains separate and frozen over the source half. Candidate
storage, endpoint conductivity and actual represented heat must satisfy the same
conservative backward-Euler transaction. Increasing conductivity does not inherit
the fixed-conductance heat-flux bracket. Without an exact conductivity tangent the
solver uses bounded bisection; nonlinear conductivity carries no uniqueness claim.
Every source evaluation in the required heat bracket must be supported. Source,
work-budget and cancellation refusal preserve all accepted state.

`MaterialWallConduction::from_phase_source` uses `PhaseWallConductionSource`
and the same `MaterialWallCaloric::from_phase_curve` chart to derive BOTH T(h)
and liquid mass fraction. Thus equal temperatures on a melting plateau can
produce different source k(T,phi); no phase interpolation is invented by gas
code. The source card must match the chart, the storage chart must match the
conduction chart, and the original density must remain compatible with the
fixed reference mass/volume. Latent heat stays entirely in h. Missing phase
coordinates, a different chart, source extrapolation or an unavailable pinned
claim refuse before state publication. Temperature-only sources continue to
refuse energy-coordinate walls. Inspectors distinguish the two source forms.

The initial-state mode retains constitutive lag. Implicit endpoint evaluation
closes that local lag but remains a first-order, fixed-geometry, lumped-wall
backward-Euler approximation with frozen gas-film/drag coefficients; neither
mode resolves a through-wall temperature field.
For centre-temperature storage, the conduction length is the surface-to-centre
distance. Adiabatic contacts refuse. Original phase-aware spatial and mechanical
owners remain responsible for moving fronts, wall support, shrinkage and liquid
motion; a thermal wall cannot remain a mechanical support after melting merely
because this heat-transfer operator still exists. Source endpoint admission is
not a continuous source-domain certificate between sampled states. Focused
G1/G3 tests check independent two-body heat balances, source-only conductivity
and capacity swaps, combined gas-film/solid response, actual chemical-rate and
pressure feedback, endpoint-domain refusal, cancellation and deterministic replay.
Phase regressions additionally check equal-T distinct latent inventories,
sensible-to-latent crossing, three-rung time refinement against an independently
integrated latent-plateau energy equation, k refresh without compounding resistance, exact
original-chart matching, pinned-source refusal and native reacting-network
energy/rate feedback. These manufactured fixtures are not measured-material
or wall-geometry qualification.

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
junctions, valves, viscous stresses, heat conduction, general moving geometry, vacuum state,
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

## Gravity-balanced isothermal atmospheres

`Euler1d::new_isothermal_atmosphere` adds signed constant axial gravity to the
existing frozen-mixture Euler owner. Its `IsothermalAtmosphere` contains typed
acceleration, temperature, pressure and position. The source's original gas
constant determines the exponential hydrostatic profile; original NASA calorics
determine energy and sound speed. Nonuniform constant-area cells receive
analytically integrated density/energy averages, with a regular small-gravity
series. Unsupported temperatures, unrepresentable profiles, nonconstant area,
inconsistent volumes and periodic boundaries refuse before publication.

Both first-order and MUSCL reconstruction operate on perturbations about the
retained reference. Source-evaluated reference faces are prepared once. Pressure
traction and the reference gravitational force cancel before the update; the
remaining density perturbation receives its actual gravitational force. The
gas evolves adiabatically after initialization: the reference is not a
thermostat and does not freeze disturbances or substitute constant gamma.
Successful steps reuse the original RK workspace and allocate no new numerical
arrays. Original source-domain admission, CFL checks, retries and cancellation
remain in force.

The shared one-, two- and three-dimensional HLLC kernel returns the physical
flux directly for identical admitted traces and the exact pressure-only flux
for a stationary pressure contact. This avoids roundoff transport from a
redundant star-state multiplication/division, including across a stationary
tangential shear jump. Ordinary shocks, moving contacts and HLL fallback still
use the original source-admitted wave solver.

Gravity's energy work uses the same face mass flux as continuity, making the
discrete sum of gas energy and `rho*phi(cell_center)` conservative up to boundary
transport and rounding. `StepReport` exposes the gravitational impulse and
work separately from area force and includes both in the balance residual.
`gravitational_potential_energy` exposes the retained midpoint potential
quadrature; this is not an exact within-cell gravitational-energy integral.
The focused G1/G3/G4 cases cover exact atmospheric averages, stationary balance,
disturbance evolution with closed energy accounting, analytic free fall under
refinement, formation-reference/coordinate covariance and cancellation rollback.

This reference-state method is bounded to a frozen ideal mixture, constant
acceleration and constant area. It does not claim balance for arbitrary
nonisothermal atmospheres, self-gravity, multicomponent/reaction transport or
an entropy/positivity theorem. Large disturbances may still refuse source
admission. Frozen networks use the explicit shared-potential constructor described
below. Arbitrary caller-owned end ports still refuse without their own declared
gravitational-potential and energy convention.

## Isentropic hydrostatic atmospheres

`Euler1d::new_isentropic_atmosphere` extends the same gravity-balanced owner to
uniform-entropy references with nonuniform temperature. It solves
`h(T)-h(T_ref) = g*(x-x_ref)` using the original positive-Cp NASA source and
recovers pressure from its source entropy difference. The lapse is `g/Cp(T)`;
neither a constant gamma nor a prescribed linear temperature fit replaces the
source. Internal-energy and entropy differences cancel their original reference
constants before arithmetic. `hydrostatic_reference()` distinguishes both
families; the existing `atmosphere()` accessor remains specific to isothermal
references.

Source-evaluated face states are retained once. Conserved density and referenced
energy are integrated with a cancellable composite four-point Gauss rule, using
`AtmosphereQuadrature`'s relative convergence tolerance and explicit panel cap.
Successive-rule differences are estimates, not certified integration bounds.
Caloric inversion retains the source owner's work cap and tolerance, tightening
the reference solve toward floating-point accuracy where the requested tolerance
allows. Source-domain escape or exhausted quadrature/inversion budgets refuses
before returning an owner. This avoids confusing spatially averaged primitive
fields with the decoded state of conserved cell averages.

The same prepared reference, perturbation reconstruction and shared-flux gravity
work preserve the discrete atmosphere and evolve real disturbances. There is
no added thermodynamic solve, integration or allocation in time stepping. Physical
gravity-force consistency converges with the reference quadrature and caloric
solve. G1/G3/G4 cases check analytic constant-Cp polytropes, independent 70-digit
variable-Cp conserved integrals, force balance, gas-plus-potential energy, formation
reference invariance and cancellation/restart. The existing constant-area,
constant-acceleration, frozen ideal-mixture and physical-boundary restrictions
apply; network coupling requires the shared-potential declaration below. This is a convectively neutral source isentrope, not arbitrary
entropy stratification, an evolving chemistry/phase equilibrium, or a global
entropy/positivity theorem.

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
or optional second-order L-stable SDIRK2, both with adaptive step doubling.
Each candidate uses the last accepted composition plus local reaction extents
through the exact admitted stoichiometric matrix. The origin is refreshed only
after acceptance; the original feed and fixed total energy remain untouched for
whole-call mass, elemental, charge, energy and entropy audits. No species clipping,
positive trace seeding or equilibrium projection is introduced. Temperature solves
the original total internal energy at each candidate composition, with the
formation contribution present exactly once. No separate reaction heat source
is added. The source-path/positive-Cv interval admission is the same L1 helper
used by the Euler EOS, including possible products absent from the feed.

Caloric roots use `sum(c_i * Delta u_i(T_initial,T))` plus the compensated
original reaction energy accumulated through accepted extents, instead of
subtracting two gross formation-energy totals. Source reaction energies are
combined before polynomial evaluation by the same L1 kinetic owner. This
preserves actual reaction binding energies and variable Cv while canceling
conserved reference offsets. Rejected stages do not change the accumulated
chemical energy. The prepared workspace retains one species and one reaction
reference vector; nonlinear trials allocate no additional storage.

The physical `energy_residual_j` uses these source-integrated differences.
`absolute_energy_residual_j` independently subtracts the original endpoint
source energy from its target. Each accepted species update also records the
original-energy equivalent of `c_new-c_old-nu*xi` as
`composition_reference_roundoff_j`. It is a numerical inventory diagnostic,
not heat. The independent absolute-energy check allows this measured term and
a magnitude-scaled floating-point allowance; this is not an interval certificate.
Element, charge and mass gates retain their original tolerances. Explicit
absolute energy targets whose floating-point spacing exceeds the requested
Joule tolerance refuse; temperature-initialized evolution retains its precise
source-relative equation. G1/G3 tests independently recover temperature from
the published product amounts, including changed binding energy and molecular
count under positive and negative common elemental references.

The Newton Jacobian combines the analytic concentration and temperature rate
derivatives with `dT/dxi_j = -sum_i(nu_ij u_i)/sum_i(c_i Cv_i)`. A bounded scaled
partial-pivot solve and positivity/energy backtracking support stiff reactions.
Stage, Newton and root scratch is allocated once and reused. Local error control
compares the complete step with two half steps, covering species and temperature;
the two-half-step result is accepted without extrapolating species. Backward
Euler is first order; SDIRK2 is second order for smooth homogeneous reactor
ODEs. Both use estimated local errors rather than global accuracy certificates.

Caller controls separately bound absolute energy, absolute amount, relative
species/temperature error, relative element/charge/mass balance, numerical entropy
decrease, minimum/initial/maximum time step and all nonlinear/root/attempt work.
Every accepted step and the complete integration check actual ideal-mixture
entropy. Final element, charge, molecular-mass and energy balances are checked again
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

### Ionized species and numerical charge conservation

The same homogeneous reactor accepts charged mechanisms admitted by the
original `fs_thermochem::kinetics::IdealGasKineticsV1` owner. Integer stoichiometry must
conserve both original elemental and charge vectors, and the declared source
molar masses must agree with every reaction. Ions and electrons contribute to
the original species pressure and internal energy. Ionization energy belongs
to those source formation references and is never added as a second heat term.
For reversible channels, the original source Gibbs energies determine reverse
rates and the concentration-normalized equilibrium constants.

Exact stoichiometric conservation does not by itself bound floating-point
charge drift. Before each accepted step is rebased, the reactor compares its
actual concentrations against the untouched initial concentrations. Final
publication also compares the actual returned mole amounts against the input
mole amounts, including the effects of volume conversion. Both checks use
`IsochoricKineticsSolve::relative_conservation_tolerance` and the residual

`abs(sum_i z_i (n_i - n_i_initial)) / max(sum_i abs(z_i) n_i, sum_i abs(z_i) n_i_initial)`.

An identically zero carrier inventory has zero residual. Charges and carrier
amounts are scaled before multiplication, and compensated sums preserve the
canonical species order. Neutral species never enter the denominator. The
audit takes linear work in species count and allocates no storage.
`IsochoricKineticsBalance::maximum_relative_charge_residual` retains the largest
accepted-step or final-publication residual. An excessive or nonfinite residual
refuses the integration with `IsochoricKineticsError::Balance`; the solver does
not clip, renormalize, neutralize or publish an inconsistent endpoint.

The conserved target is the authored feed charge, not necessarily zero. This
is a homogeneous constituent/energy evolution law: admitting a nonzero charge
does not establish an electrostatic equilibrium or supply its field energy.
All species share one source temperature and ideal-gas convention. Separate
electron and heavy-particle temperatures, electric-field acceleration, ambipolar
diffusion, non-Maxwellian electron distributions, radiation-driven ionization,
and nonideal plasma corrections require separate owners and are not supplied
by this reactor.

Both owned `IsochoricKineticsResult` and borrowed `IsochoricKineticsView` expose
`frozen_state`. This explicit owning conversion retains the accepted actual
amounts, source bank, volume and common temperature in L1's
`FrozenIonizedGasStateV1`, for nonequilibrium charged-transport consumers. It
performs no Saha/equilibrium projection. The original reactor still preserves
an authored nonzero feed charge, but conversion to this quasineutral material
lane refuses it. State cloning shares the immutable sources and populations;
the prepared reactor's ordinary advance remains allocation-free.

The ionization regressions use `H + e- <-> H+ + 2 e-`, ideal monatomic
Sackur--Tetrode thermodynamics and NIST's H I ionization energy. They compare
detailed balance and the final adiabatic state with independent Saha formulas,
and finite-time evolution with an independent scalar integration. They check
electron pressure, ionization cooling, mass/energy/charge balance, preserved
nonzero feed charge, and refusal of trace-carrier drift over extensive scales.
Exactly zero initial electrons remain zero for this electron-impact channel;
no artificial seed is inserted. Its forward rate is explicitly synthetic, so
these numerical tests do not validate an empirical ionization time or discharge
mechanism. Source reference: NASA RP-1311, Part I, section 3.7,
https://ntrs.nasa.gov/citations/19950013764, and NIST H I ionization energy,
https://physics.nist.gov/cgi-bin/ASD/ie.pl?spectra=H+I.

An additional source-backed regression runs the original stiff SDIRK2 reactor
with L1's NIST H BEB electron-impact law, including its Maxwellian temperature
derivative and the original Gibbs-derived reverse law. Ionization increases
the electron count and consumes thermal energy while mass, nuclei, charge and
total internal energy remain conserved. This proves that the retained physical
rate reaches the existing evolution path; it is not empirical validation of
ionization time or an electrical discharge.

## Frozen electron/heavy energy evolution

`chemistry::energy::advance_frozen_energy` evolves the actual source-bound
`FrozenIonizedGasStateV1` at fixed species amounts and volume, with independent
constant electron and heavy-particle powers over the requested interval.
Positive power deposits energy and negative power removes it. An absent
electron population cannot accept electron power or acquire a seed. Electrical
consumers supply their actual, separately resolved carrier work; this owner
does not invent a field or add total Joule work a second time.

Elastic energy exchange uses the same source-retained electron-neutral
momentum cross sections as material drag:
`G = 3 R n_e sum_target[N_target (m_e/m_target) k_m(Te)]`, with electron amount
`n_e` in mol, neutral number density `N_target` in m^-3 and
`Q_e_to_h = G (Te-Th)`. Every present neutral requires exactly one original
momentum source. The declared target/electron speed limit also bounds
`sqrt(m_e/m_target)`, so cold but unrealistically light targets cannot pass a
massive-target approximation. Charged fraction, source temperature intervals,
cross-section coverage and numerical work remain explicit. These limits do
not establish that neglected Coulomb or inelastic processes are insignificant
for an arbitrary experiment. The model is Banks (1966), Eqs. 16–20,
https://ntrs.nasa.gov/api/citations/19660020172/downloads/19660020172.pdf.

Each first-order substep freezes G at its starting state and solves a monotone
implicit scalar exchange. The same signed exchange enters the two pool energy
targets with opposite signs. Full/two-half comparison controls an estimated
local error, including the changing collision rate; it is not a global error
certificate. Every caloric inversion uses original-source internal-energy
increments relative to the original temperatures, canceling formation offsets
before rounding. Energy and temperature residuals both gate convergence, so
a trace electron pool cannot hide a large Te error under a parcel-wide energy
tolerance. Nonlinear trials allocate no heap storage. Source gaps, invalid
domains, exhausted work and cancellation publish no partial endpoint.

The result retains the initial and final source-bound states, original
collision sources, original controls, separately deposited energies, signed
exchange and actual per-pool energy residuals. This is a frozen-composition
thermal step. It does not supply two-temperature molecular detailed balance,
ionization population evolution, photon absorption, field energy or charged
spatial transport; those owners must exchange their own conservative sources.

## Common material gas integration

Optional feature `material` exposes
`material::advance_closed_parcel_adiabatic`. It evolves an existing
`fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel`, using its
actual stored internal energy as the one target. The complete product bank
must match the mechanism's original NASA-9 models and retain every active
feed transport card, source-pack identity and declared water identity.
This adapter refuses mechanisms containing charged species before integrating:
its Wilke viscosity and Wassiljewa--Mason--Saxena thermal-conductivity closure
belongs to neutral dilute gas. Source identity alone does not make that closure
an electron/ion transport model. A charged-mechanism refusal leaves the
original parcel unchanged, including for a zero-duration request.

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

### Caller-owned end ports and local cell sources

`Euler1d::step_with_ports(max, controls, ports, cancelled)` advances a
non-periodic interval whose end fluxes are owned by the caller. The SSP-RK2
predictor offers stage 0 at `t^n` and the corrector stage 1 at `t^n + dt`, each
as a `PortStage { stage, time, dt, left, right }` holding the end cells'
conserved and primitive states; the closure returns `[Option<FaceFlux>; 2]`
(`None` keeps the declared boundary). A returned flux is the whole end
exchange: the accepted report's boundary changes are `dt * area * flux`. A
refusing port refuses the step (subject to the declared halvings) and leaves
field and time unchanged; periodic intervals refuse ports.

The relations a port needs are public and source-owned:

- `FrozenGas::end_wave(interior, normal, pressure, cancelled)` returns the end
  state reached from `interior` along its outgoing wave curve (isentrope in
  expansion, Hugoniot in compression), its velocity derivative with respect to
  the end pressure and its signal speed. Its weak limit is `du = -n dp/(rho c)`
  and it stops gas at the wall solver's pressure.
- `FrozenGas::nozzle_exit(donor, back, cancelled)` is the steady isentropic
  exit of a stagnation donor into a back pressure, with choking; its
  incompressible limit is Bernoulli and its stagnation enthalpy is the donor's.
- `stagnation_enthalpy(T)`, `face_flux(face, signal)` and
  `inflow_flux(face, T0, signal)` (energy flux = mass flux times the donor's
  stagnation enthalpy) close a face from one state.

The shared nozzle closure resolves the imposed pressure difference even near
equal reservoir pressures. It factors the temperature decrement by the original
source's isothermal entropy increment `R ln(p0/pb)` and evaluates the source's
Cp/T and enthalpy secants in that coordinate. The residual is relative to the
driving increment; an absolute Cp-scale tolerance cannot introduce a flow floor.
The derivative limit retains finite nozzle kinetic energy when the temperature
decrement is smaller than one representable absolute-temperature step. Returned
temperature and EOS state still have ordinary floating-point rounding.

Positive source Cv gives a conservative subsonic shortcut: when
`2 ln(p0/pb) < pb/p0` and the bounded cooling interval is source-supported, the
closure avoids the preliminary sonic solve. Stronger pressure drops retain the
existing source-owned sonic limit and choking decision. Neither route freezes
gamma or extends the temperature domain. The same closure serves public frozen
gas end ports and finite-volume network valves. Focused regressions compare
constant/variable-Cp analytic isentropes down to adjacent pressures, the exact
zero-drop limit, reversed donor flow, energy-reference invariance, cold-source
refusal and cancellation. A counter over the actual FrozenGas source compares
the complete weak solve with the previously mandatory preliminary sonic solve;
this is a source-query work comparison, not a throughput or scaling benchmark.

`Euler1d::apply_cell_sources(update, cancelled)` replaces every cell average
through a caller-owned local update between accepted steps (wall friction, wall
heat exchange). Every candidate is decoded through the original caloric law
before any is published; one refusal leaves the whole field unchanged. Time
does not advance.

Tests (`dynamics::evolve::port_tests`, 6): wall ports reproduce reflecting ends
bit for bit with ordered stages; an injecting port is the whole end exchange and
conserves mass to 1e-15 and energy to 1e-14 relative; refusal and periodic
refusal; the end wave's acoustic limit and wall pressure (1e-9); the nozzle's
Bernoulli limit (1e-4 at a 1e-4 pressure drop) and face-flux consistency with
the Riemann interface; transactional and exact cell sources. Measured on trj
lane t, log 20261030T190307Z (61 `dynamics` tests pass).

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

`MaterialBinaryTransport` supplies the complete molecular law directly from the
original material gas and explicit `BinaryDiffusionMaterial` pair assignments.
The original caloric source is compared in full; each unordered species pair
must occur exactly once, including currently absent species. Pair cards retain
their original external contexts and pin the selected binary-diffusivity claims
at their declared reference T,p. Actual face T,p refresh those same claims;
missing support refuses rather than switching to another source or estimating
D from viscosity. Binary data have units m²/s and must be positive. Pair contexts
cannot freeze mixture-composition coordinates or aliases of evolving T,p.
Source axes and applicability bounds require supported static external meanings;
untyped restrictions and frequency/field/rate support refuse. Unused context
labels remain in the original query receipts, and declared external coordinates
such as humidity remain fixed while T,p evolve.

This law enters the existing axial, planar, spatial, finite-parcel and reacting
owners through `TransportLaw`; it does not introduce another transport solver.
The shared continuum owner retains its mass frame, pressure-driven separation,
original partial-enthalpy flux and full energy accounting. The original
`MaterialSoretDufourTransport` may wrap this law with independent cross data.
Scalar viscosity/conductivity still come from the original material gas, and
bulk viscosity remains explicit. Failed evaluations preserve the output and
last complete pair receipts. Source queries allocate and carry their existing
validity limits; this path is not an empirical validation or performance claim.
G0/G1/G3/G4 tests in `material::tests::transport::binary` check unequal-mass
composition/pressure diffusion and partial enthalpy against independent binary
answers, conservative spatial mode refinement, data-only diffusivity changes,
pair symmetry/coverage, source pinning and cancellation.

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
species diffusion, viscosity, general moving geometry, arbitrary junction or multi-D flow
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
general moving meshes and arbitrary reacting junction networks remain separate work.

### Finite-rate chemistry driving a resolved inertial piston

`ReactingEuler1d::step_with_piston` composes the original cell chemistry half
steps around `MixtureEuler1d::step_with_piston`. The first source uses the old
cell volumes; the second uses their actual moved volumes. Changes in species
and temperature therefore drive the original pressure waves and inertial wall.
The reactor's existing physical-energy target and species-reference correction
supply chemical conversion once; no additional heat-release source is added.
The constant-area interval has a fixed left wall and one inertial right wall.

Every complete attempt starts from the accepted field, geometry and piston
velocity. A shorter mechanical/acoustic duration restarts both chemistry half
steps at that duration. Second-half source refusal, failed whole-step balance
or cancellation cannot publish the successful middle pressure/motion stage.
The existing two field workspaces, reactor and amount buffer are reused. A
bounded, cancellable same-owner geometry reset copies only private trial
coordinates/volumes after motion; the ordinary stationary scratch-copy and
molecular-owner geometry guards remain intact.

`ReactivePistonStep` retains the original reactive and piston reports and adds
actual physical chemical-gas energy change plus whole-step momentum/energy
residuals. Endpoint species-reference changes participate in the physical
energy balance, together with gas energy, piston kinetic energy and exterior
pressure/additional-force work. The fixed left support carries its explicit
opposite impulse with zero work. The same caller-authored piston tolerances
gate whole-step balances as well as the pressure solve; callers must account
for their per-cell chemistry tolerances when choosing that global budget.

Four focused G1/G4 regressions cover unforced exothermic pressure-driven
expansion and species/work/support balances; acoustic restart with independent
stage replay at the accepted duration and moved volumes; a second chemistry
refusal after successful compression; and cancellation/recovery after prior
motion. They are authored but unexecuted under the session's worker resource
limit. This entry does not advance the separately bound fixed-geometry
molecular, finite-wall heat or network owners. Original chemistry method,
splitting-error, source precision and unequal-caloric-contact limits remain;
there is no general second-order, entropy-stability or combustion-validity claim.

### Molecular stress and heat transport through a reacting inertial piston

`ReactingEuler1d::step_with_piston_and_dissipation` evaluates the original
`D(dt/2), C(dt/2), A(dt), C(dt/2), D(dt/2)` composition inside the shared
midpoint-speed piston root. Every proposal restarts all five stages from the
same accepted field at one duration and wall speed. The molecular stages use
the old and moved physical widths/volumes, respectively; only the original
extensive Euler A stage changes geometry and time. Chemical conversion,
Maxwell--Stefan species enthalpies, Fourier/Soret/Dufour transport and the
longitudinal Newtonian stress remain with their existing source owners.

Both pressure and molecular wall stresses drive piston inertia. Molecular
wall gradients use the candidate wall velocity, and the very same signed
momentum flux supplies laboratory wall work `w*(-tau)` in total energy.
`DissipationReport` now retains separate signed left/right conservative
transfers, including network branch reports. The returned
`PistonStep<ReactiveStepReport>` combines those impulses/work with the moving
Euler pressure exchanges. Actual endpoint chemical gas energy, piston kinetic
energy, fixed-left support impulse and exterior-pressure/additional-force work
must satisfy the caller's original physical budgets before publication. No
separate heat-release or viscous-heating source and no energy correction is added.

The constant-area interval has an impermeable adiabatic fixed left wall and
one impermeable adiabatic inertial right wall. Conduction and species diffusion
act inside the gas; there is no conductive wall-reservoir heat to count as
mechanical work. Finite wall storage, prescribed thermal end reservoirs,
variable-area or multidimensional motion remain outside this entry.

The molecular operator retains its accepted geometry separately from one
fallibly allocated candidate mesh. Trial widths are refreshed for every
stage, so a refused moved candidate cannot contaminate a later stationary
step. A shortened acoustic duration restarts the complete composition and
root. Gas, mesh, clock, piston and the molecular geometry binding publish
together after final acceptance and cancellation polling. Source and molecular
scratch are reused; no complete field is cloned per velocity proposal. The
original stationary molecular entry still requires its exact accepted mesh.

Four authored G1/G3/G4 regressions cover an independent one-cell RK2 viscous
wall force/work oracle; reacting expansion with nonzero viscous load, physical
energy/support balance and source-gauge changes; complete replay at the solved
speed after acoustic shortening; and final cancellation after prior motion,
scratch reuse and ordinary stationary continuation on the accepted geometry.
They are unexecuted in this resource-limited session. Original parabolic,
source, chemistry/splitting and unequal-caloric-contact limits remain; no
experimental combustion, general second-order, entropy or FSI-stability claim
is made.

### Finite end-wall heat driving reacting piston pressure

`ReactingEuler1d::step_with_piston_dissipation_and_heat` adds the original
caloric heat exchange around the same molecular/reacting midpoint-speed root:
`H(dt/2), D(dt/2), C(dt/2), A(dt), C(dt/2), D(dt/2), H(dt/2)`.
`PistonThermalEnds` retains two distinct finite wall stores or prescribed
thermal reservoirs, their original caloric bindings and accepted gas geometry.
Every velocity proposal resets both stores from accepted values; its first
heat half uses old end-cell volumes and its last uses the moved volumes.
Heat therefore changes the actual caloric state used by pressure fluxes and
continuing motion. End order is left then right, including when both ends
contact the same gas cell. Only A advances geometry and time.

Pressure and longitudinal viscous stress remain the only gas-side mechanical
impulse/work terms. Conductive heat never enters the stress-work ledger. The
independent complete energy admission is
`Delta E_gas + Delta U_finite_walls + Delta K_piston + P_ext*A*Delta x
- F_ext*Delta x - Q_prescribed_reservoir = 0`.
Gas chemical energy includes actual endpoint species-reference changes once.
Finite-wall storage is reconstructed from accepted endpoints using the original
energy difference or nonlinear capacity integral; a cached capacity tangent is
not multiplied by the whole temperature change. Arbitrary wall energy datums
cancel in differences. Prescribed reservoir heat is external, whereas heat
exchanged with a finite store remains internal to this energy balance.

`ThermalReactiveStepReport` retains complete gas endpoints and heat transfers
alongside the original middle D/C/A/C/D report. That middle report's balances
exclude the outer heat halves. `PistonStep::total_energy_residual_j` includes
the finite stores and subtracts external reservoir heat on this entry;
`mechanical_energy_residual_j` continues to contain mechanical work alone.
The end owner's last report, contacts and accepted geometry publish atomically
with gas, clock, piston and molecular geometry after the final cancellation
poll. Acoustic shortening restarts all seven stages and the whole root.
Two small contact arrays and a prepared candidate mesh are reused; no thermal
source tables or complete fields are cloned per velocity proposal.

Source-backed storage must have an explicit internal-energy or admitted
negligible-expansion solid convention. Bare source enthalpy refuses because
its pressure/deformation work has no owner in this rigid-body balance. Radial
drag also refuses, avoiding a duplicate or stationary-wall momentum sink on
the moving piston. End conductance is caller-prescribed hA; finite thermal
masses represent distinct fixed amounts and do not add a second inertial mass.
Thermal expansion, axial wall conduction, latent material work not represented
by the admitted source, changing area and a fully implicit conjugate thermal
solution are outside this entry. Backward-Euler heat halves do not establish
global second-order accuracy or an entropy/combustion-validity theorem.

Four authored regressions cover an independent two-body heat oracle and hot/
cold finite-wall pressure response with full storage/work balance; nonlinear
original wall storage with simultaneous chemical and wall energy-gauge shifts
and a separately counted external reservoir; final cancellation after prior
motion, exact replay and original source-interval refusal; and admission of
the rigid thermal-work/drag convention. These tests are unexecuted under this
session's worker resource limit; no native build or solver result is claimed.

### Two resolved gas chambers sharing an inertial piston

`MixtureEuler1d::step_with_shared_piston` couples the original evolving
multicomponent fields on BOTH sides of one impermeable inertial wall. Each
chamber retains its own source, canonical species axis, caloric domain and
resolution. Their constant cross sections, shared interface coordinates and
physical clocks must match. The wall has negligible axial thickness; both
outer ends are fixed and reflecting. Positive axial motion expands the left
volume and compresses the right volume.

Each proposal in the original midpoint-speed root evaluates both original
moving-volume RK2 fields at one common duration from their accepted states.
The pressure on each side is its resolved moving-wall wave pressure, never
an imposed opposing reservoir or a guessed whole-chamber pressure law.
Piston impulse is the opposite of the summed interface gas impulses. The
very same face exchanges supply gas work; total gas energy plus piston
kinetic change equals the applied additional-force work. Fixed support
reaction is the opposite of BOTH outer gas impulses, with zero support work.
`SharedPistonGasReport` retains the two separate original conservation and
boundary reports, including their independent species balances and gauges.
The outer `PistonStep` sums fixed support reactions; exterior-pressure work
is zero because both gas sides are already represented.

The initial timestep obeys both relative acoustic limits and both possible
mesh-contraction limits over the authored speed bracket. A refused stage or
root halves and restarts the complete two-field proposal. Each original
field reuses its candidate/RK/geometry workspace; no source tables or full
fields are cloned per velocity proposal. After the last cancellation poll,
both extensive fields, meshes, equal clocks and piston velocity publish
together. Shared candidate endpoints must represent exactly one position.
`gas_trials` counts paired proposals, each requiring two original field solves.

Four authored G1/G3/G4 regressions cover equal-pressure rest with different
gases and independent opposing p*A*dt impulses; unequal-pressure waves in
both directions, source-energy gauge invariance and full work/support balance;
fine-side CFL shortening with replay of both prescribed moving fields; and
late cancellation, continuation plus incompatible geometry/clock refusals.
Native compilation, formatting and runtime execution remain unverified after
the delivery executor became unavailable. This is sealed adiabatic Euler
transport, without advancement of separately bound chemistry, molecular or
thermal owners. Original unequal-caloric contact and boundary-order limits
remain; there is no general FSI stability or multidimensional-rung claim.

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
surface kinetics and wall heat transfer remain separate models.

### Resolved moving-wall pressure and mechanical work

Frozen and evolving-mixture planar/spatial face owners additionally expose
`moving_wall_normal_flux` and `moving_wall_spatial_flux`, each with a
source-admitted primitive fast path. These are moving-face ALE fluxes. The
original one-sided shock/rarefaction owner uses normal velocity relative to the
wall and retains the same density, temperature, pressure, actual species
composition, temperature-dependent calorics and source-domain refusals.

In the supplied normal axis the flux contains exactly zero relative mass and
species transfer, normal momentum `p_star`, zero tangential momentum, and
laboratory-frame energy `p_star * wall_velocity`. Multiplication by the existing
-1/+1 outward-normal sign gives outward transfer. Thus a piston doing work on
the gas adds internal and kinetic energy through its pressure traction, with
no duplicated heat or species-enthalpy transfer. Tangential translation has no
effect on the inviscid slip pressure or work. Independent species formation
references cancel because no species crosses the wall.

The reported signal-speed bound is relative to the moving face, not the
laboratory frame. Zero wall velocity retains stationary arithmetic and exactly
zero work. Both adapters use the existing shared caloric wave kernel; successful
face work creates no new gas/source owner and needs only bounded stack scratch.

Seven added regressions cover independent constant-Cp piston shock/rarefaction
pressures, variable-Cp swept-slab mass/momentum/energy balances, evolving-mixture
species and formation-reference invariance, normal/tangential Galilean changes,
relative wave-speed semantics, planar/spatial and conservative/primitive
agreement, stationary identity, source limits and cancellation without partial
work publication. They have not been executed in the current resource-limited
session. This remains an instantaneous inviscid face closure: moving geometry,
cell-volume evolution/geometric conservation, temporal quadrature, piston
inertia and viscous wall work are the corresponding spatial/mechanical owners'
responsibility. The separate moving-interval entry below supplies that spatial
evolution for a bounded axial geometry; ordinary stationary updates are unchanged.

### Spatial piston dynamics and conservative moving volumes

`Euler1d::step_with_moving_walls` advances the original resolved Euler field in
a constant-area interval with two impermeable prescribed end walls. Their two
axial speeds are fixed over each accepted step; all interior faces move by
affine interpolation, retaining the original nonuniform cell fractions.
The existing MUSCL reconstruction, source admission and HLLC/HLL wave owner are
shared with stationary dynamics. An interior Riemann problem is solved in its
face's translating frame, then its flux is returned to laboratory conserved
variables as `[f_mass, f_momentum+w*f_mass,
f_energy+w*f_momentum+(w*w/2)*f_mass]`. Thermodynamic traces are reconstructed
through the original source in that frame, avoiding subtraction of large
laboratory kinetic energies. Wall pressure and `p_star*w` work use the existing
one-sided caloric shock/rarefaction closure.

Both RK2 stages advance extensive `Q=V*U`, with the predictor resolved at the
new volume and the corrected endpoint on that same new geometry. Cell volume
changes are exactly the discrete swept volume `dt*A*(w_right-w_left)` used by
the moving-face flux. Actual face-relative CFL is checked at both old and new
stage volumes, and contracting geometry also bounds the initial trial. The
existing bounded retry/source-refusal policy never clips conserved fields.
`StepReport` retains the same stage-integrated left/right pressure impulses
and wall work that changed the gas; its conservation residual subtracts those
boundary exchanges once. Geometry scratch is fallibly allocated at the first
nonzero-motion entry and reused on subsequent steps. Cancellation/refusal
preserves the accepted mesh, cells and clock together; zero speeds retain the
ordinary stationary path and its exact arithmetic.

Six focused regressions cover independent Reynolds-transport and Galilean flux
relations, zero-speed identity, a co-moving nonuniform pressure-contact field,
resolved variable-Cp piston compression and formation-gauge invariance, a
16/32/64-cell homologous-compression comparison against the independent smooth
isentropic Euler solution, and early/final cancellation or geometry/source
refusal. They are authored but unexecuted in the resource-limited session.
The original unreconstructed wall-cell traces can limit global boundary order;
second-order global convergence is not asserted. This lane does not evolve
piston inertia, heat transfer, composition, radial/variable-area deformation,
open ends, bound gravitational atmospheres or multidimensional moving meshes.
The face-frame conservation convention follows Springel (2010), §3.3,
<https://doi.org/10.1111/j.1365-2966.2009.15715.x>; it does not claim that paper's
Voronoi mesh algorithm or its full accuracy properties.

`MixtureEuler1d::step_with_moving_walls` supplies the same bounded axial motion
for the original evolving-composition field. It shares the geometry workspace,
affine face velocities and swept-volume preparation with the frozen owner,
while retaining its own compatible mass-fraction reconstruction, canonical
species storage, source calorics and shared HLLC/HLL waves. Every partial
density advances extensively on the same RK2 volumes. Face-frame species
fluxes are unchanged by the laboratory transformation; their compensated sum
is the mass flux used in both momentum and kinetic-energy transformations.
There is no independently rounded mass-flux field, composition normalization
or absent-species seeding. Each species keeps its own formation-energy datum.
All wall species fluxes are zero, and the returned original mixture report
accounts for the pressure impulse and mechanical work at both ends.

The shared geometry's cancellation is translated into the mixture's native
`Cancelled` refusal before retry handling. Accepted mesh, all species, energy
and clock remain atomic. Five focused regressions cover independent species
Reynolds-transport fluxes, stationary identity, co-moving unequal-caloric
contacts, nonuniform-composition piston work and independent species gauges,
single-species agreement with the original frozen owner, structural zeros,
and early/late geometry or source refusal. These are authored but unexecuted
under the same worker resource limit. This standalone entry does not advance
the separately bound stationary reaction, diffusion, heat or network owners.
Generally advecting unequal-caloric material contacts retain the original
conservative scheme's pressure-oscillation limitation; no exact general
pressure-equilibrium or entropy-stability claim follows from this extension.

### Simultaneous inertia of both resolved gas end walls

`Euler1d::step_with_inertial_walls` advances a closed constant-area field and
two finite `Piston` masses together. Both wall velocities use the increasing
axial coordinate. Each `InertialWallLoad` holds exterior pressure and an
additional signed force: exterior pressure pushes the left wall rightward and
the right wall leftward. Gas momentum remains in the original resolved field.

A bounded two-variable Newton solve selects both midpoint wall speeds.
Every trial runs the original moving-volume RK2 operator from the same
accepted state and duration. Each wall receives the negative of its actual
gas-side pressure impulse; the identical boundary work enters its kinetic
balance. Both wall momentum/work residuals and combined gas-plus-walls
momentum/energy are admitted independently, using represented mesh
displacements for exterior/load work. No fixed support, reduced wall mass,
uniform cavity pressure or posterior energy correction replaces an owner.

Both velocity intervals bound relative-wave and contraction CFL. Row-scaled
finite differences and bounded backtracking reuse the native candidate and
geometry buffers. An exhausted solve retries the complete interval under
the original gas halving budget. Only the last evaluated native candidate
can publish; gas, mesh, clock and both wall velocities commit after the
final cancellation check. No field or source table is cloned per trial.

Six regressions cover balanced outside pressure, symmetric vacuum recoil,
unequal walls with nonuniform pressure and signed loads, a common velocity
boost through the original first-order ALE route, the massive-left-wall
limit of the original one-piston owner, and cancellation/refusal rollback.
They are authored but unexecuted under the session's worker resource limit.
This frozen-composition adiabatic lane does not yet compose the two moving
walls with mixture chemistry, heat transfer, gravity, contact stops or
transverse dynamics; it establishes no general FSI stability theorem.

### Evolving mixtures with two finite end walls

`MixtureEuler1d::step_with_inertial_walls` now uses the same simultaneous
wall solve and physical acceptance with its original species-resolved ALE
operator. Local partial densities and original caloric sources determine
pressure at every stage. The complete native gas report retains canonical
species masses, both boundary exchanges and independently referenced total
energy. Each impermeable wall has zero flux of every species; absent species
remain absent. No nominal frozen composition replaces the evolving field.

`InertialWallsStep<R>` retains the frozen report as its public default. A
small internal exchange record passes only actual gas momentum/energy changes
and each boundary impulse/work to the shared mechanical equations; no species
or source state enters a duplicate wall model. Trials reuse the native mixture
candidate and geometry workspaces. Gas, all species, mesh, time and both walls
publish only after the common physical checks and final cancellation gate.
Original mixture failures retain their typed causes; cancellation and
allocation failures propagate without numerical retries.

Three focused regression groups cover unequal walls with heterogeneous
composition and independent formation-gauge shifts; agreement with the
original frozen owner for a fixed two-species composition; and early/late
cancellation, failed roots and unsupported-boundary rollback. These tests
are authored and formatting checked, but have not been compiled or executed
in this constrained environment. Chemistry, molecular diffusion, viscosity
and wall heat remain separate owners. The underlying conservative mixture
scheme retains its pressure-oscillation limitation at generally advecting
unequal-caloric contacts.

### Finite-rate chemistry with two inertial gas end walls

`ReactingEuler1d::step_with_inertial_walls` composes the original source
half steps around the simultaneous two-wall mixture pressure solve. Both
wall momenta remain finite and gas momentum stays spatially resolved. The
first chemistry half uses accepted cell volumes and the second uses actual
moved volumes. A shortened native interval restarts both source halves and
both walls from accepted state at that same duration.

`ReactiveInertialWallsStep` retains complete reactive accounting and the
original middle-stage wall report. Whole-system momentum subtracts both
signed exterior-pressure/additional-force impulses, with no fixed-support
term. Whole-system energy uses actual chemical gas-energy change, both wall
kinetic changes and represented exterior/additional-force work. Each changed
species reference participates once through the original coupled report;
no separate chemical heat or posterior energy correction is added.

All source, field and geometry workspaces are reused. Species, mesh, clock
and both wall velocities publish together after complete physical acceptance
and the final cancellation poll. Four authored G1/G3/G4 groups specify
unequal-wall reactive recoil with an independent energy-gauge comparison,
accepted-duration replay using moved source volumes, postcompression source
refusal and cancellation/retry after prior motion. They parsed with the
pinned formatter but have not been compiled or run in the constrained
editing environment. Molecular transport and thermal walls remain separate
owners. Existing chemistry/splitting and unequal-caloric-contact limits
remain; no general accuracy, entropy or combustion-validity claim is added.

### Pressure-driven inertia of a resolved gas end wall

`Euler1d::step_with_piston` couples the existing spatial Euler field to
`dynamics::piston::Piston`, a finite inertial mass at the right endpoint. The
left wall remains fixed, both ends are impermeable, and area is constant.
The mesh alone owns wall position and area. Authored exterior pressure and
additional force act on the piston; neither a uniform chamber law nor a
second pressure model replaces the source-owned wave field. The existing
zero-dimensional `fs-couple::cavity_phs::piston` owner remains independent.

Each bounded scalar trial chooses midpoint velocity `w`, then runs the
original moving-volume RK2 operator from the accepted field at the SAME
duration. Its accepted right-wall pressure impulse `I` enters
`M*(v1-v0)=I+dt*(F-A*p_ext)`, with `v1=2*w-v0`. The identical two-stage wall
pressure supplies gas work. Consequently `Delta K=w*M*Delta v` pairs the
mechanical update with the discrete gas exchange; no subsequent energy
correction is applied. The report independently measures actual gas energy,
wall kinetic change, exterior-pressure and force work, and both mechanical
and combined momentum/energy defects. The fixed-left-wall reaction is
explicit in the combined momentum balance and does zero work. External work
uses represented mesh displacement, so a rounding discrepancy is admitted
only within the same caller-specified physical energy tolerance.

`PistonControl` requires a source-supported midpoint-speed bracket containing
the accepted velocity, positive momentum/energy tolerances, and a bounded
secant/bisection budget. Bracket speeds bound the initial relative-wave CFL;
both gas stages retain original source, moving-geometry and CFL admission.
Any retry halves and resets the complete coupled attempt. Only the last
successfully evaluated candidate can publish gas, geometry, clock and piston
velocity, after a final cancellation poll. Original gas workspace is reused;
geometry scratch is allocated fallibly on first use. There are no per-trial
field clones, and cancellation is never converted into a numerical retry.

Five focused G0/G1/G3/G4 regressions cover balanced pressure, variable-Cp
compression with source-energy-gauge invariance and global work/momentum
closure, exact-duration replay of the final solved motion through the
prescribed-wall owner, three mesh refinements against the independent
small-amplitude acoustic piston-impedance solution before reflections, and
early/late cancellation or source/bracket/geometry refusal. These tests are
authored but unexecuted under the session's worker resource limit. The
acoustic oracle is a linear small-amplitude limit; neither it nor residual
closure establishes experimental validity or a general FSI stability theorem.
Existing spatial-order and unequal-caloric-interface limits remain unchanged.
This entry does not supply piston collisions, friction, springs, crank
linkages, multiple moving inertial walls, heat or reactive/multicomponent
coupling. Context: Lefrançois and Boufflet, SIAM Review 52(4), 747–767 (2010),
<https://doi.org/10.1137/090758313>; the discrete work identity above is derived
directly for this implementation, not inherited from that paper.

`MixtureEuler1d::step_with_piston` now supplies the same inertial end wall for
the original evolving-composition field. It uses the identical `Piston`,
`PistonLoad`, `PistonControl`, scalar root and mechanical balance code. A small
source-neutral exchange carries actual gas momentum/energy changes, both
wall impulses, right-wall work and elapsed duration into that shared code.
`PistonStep<R>` retains the native gas report: the public default remains the
original frozen `StepReport`, while the mixture method returns its original
species-bearing report. No species, total-mass or energy-reference axis is
collapsed into a frozen surrogate.

Every mixture trial calls the existing extensive species ALE operator at the
common nonlinear duration. Actual partial densities, independent formation
references and source calorics determine pressure throughout the field.
Closed walls have zero flux of every species, and structural zeros remain
zero. The complete accepted field, geometry, clock and piston publish only
after shared momentum/energy acceptance and final cancellation polling.
Original mixture errors remain distinguishable as `PistonError::MixtureGas`;
native source/geometry cancellation is converted to `PistonError::Cancelled`
before retry handling. The prescribed-motion mixture path shares only its
original report extraction and trial visibility; its arithmetic is unchanged.

Four additional G0/G1/G4 regressions cover heterogeneous-composition inertial
compression with every species balance, independent source-gauge shifts and
global work/support exchange; single-species agreement with frozen inertia;
replay of the final solved midpoint speed through the original prescribed
mixture motion; and cancellation plus distinctly typed source, bracket and
geometry rollback. They are authored but unexecuted under the same worker
resource limit. This entry does not advance separately bound reaction,
diffusion, heat or network owners. The existing pressure-oscillation limit at
general moving unequal-caloric contacts remains explicit; shared mechanical
balance acceptance is not an entropy or pressure-equilibrium theorem.


## Conservative molecular, thermal and viscous field evolution

`multicomponent::dissipation::DissipationOperator<L>` binds the original gas,
stationary constant-area mesh, boundary conditions and explicit transport law.
`advance_for` completes an exact requested duration through bounded SSP RK2
parabolic substeps. It publishes the complete field atomically and leaves its
clock unchanged for the composing integrator. Nonuniform spacing is admitted
with the actual distances and physical cell volumes; variable-area viscous
geometry requires additional stress terms and is refused.

Every shared face carries barycentric Maxwell-Stefan mass fluxes Ji, momentum
flux `-tau`, and, without thermal cross coupling, gauged-energy flux
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
fully discrete entropy/invariant-domain proof. Explicit source-owned Soret/Dufour
coupling is described below; nonideal phases, Knudsen transport and variable-area
viscous stresses remain separate closures.

### Explicit reciprocal thermal diffusion

`TransportLaw::heats_of_transport` optionally supplies canonical finite kinetic
heats of transport `H_i [J/kg]` at the actual common face T,p,composition.
`SoretDufourTransport<L>` supplies a constant authored approximation on an
explicit T/p domain and the full composition simplex, with source identity,
version and typed species coefficients matched to the underlying caloric bank.
The wrapped molecular law retains its own source/domain checks. No coefficient
is inferred from species names, molecular masses, heat capacities or collision
diameters. State-dependent or composition-restricted laws implement that same
source query directly. Existing molecular laws omit the extra closure by default.

With `material`, `MaterialSoretDufourTransport<L>` supplies this closure from
original material-card scalar, curve or two-coordinate surface claims. Each
canonical species is explicitly assigned a signed J/kg kinetic coefficient
key, card and selection policy, and an explicit typed mole-fraction query axis.
The caller declares the wrapped law's conductivity at zero diffusive species
flux; partial enthalpies cannot substitute for these kinetic coefficients.
Construction resolves a reference state, then pins the exact selected claims.
Every face refreshes the original sources at actual T,p and composition,
preserving fixed external conditions and full usage receipts. Original static
source admission rejects unowned frequency, field or rate response and
unsupported scalar constancy. Missing/ambiguous claims, type/domain/source
refusal and cancellation publish no partial coefficients or new retained cross
receipts. A successful face's coefficient and heat requests reuse one resolved
cross query; source resolution allocates and is not an allocation-free kernel
claim. Existing material transport omits cross coupling by default.
Manufactured surface fixtures independently verify the resulting Soret/Dufour
flux and entropy across changing T,p,x, plus source pinning, absent species,
typed admission and atomic refusal; they are not measured-material qualification.

The shared axial/network/planar/spatial face adds the balanced Soret force
`d_i^T = x_i M_i (H_i-Hbar)/(R T) grad(T)/T`, with `Hbar=sum_i Y_i H_i`,
to the original ideal chemical/pressure forces before calling the same
Maxwell--Stefan owner. Its reciprocal reduced heat flux is
`q_reduced = -lambda grad(T) + sum_i (H_i-Hbar) Ji`.
Here lambda means conductivity at zero diffusive species flux. These kinetic
cross coefficients are distinct from original partial enthalpies: total energy
contains partial-enthalpy diffusion, this reduced heat flux and stress work
exactly once each. A common shift of every H is removed before arithmetic;
zero coupling recovers the original molecular law. Exact absent species are
neither seeded nor assigned spurious flux. At impermeable walls, zero species
flux makes the reciprocal heat contribution zero.

For positive fractions, the local conjugate entropy balance is
`-R sum_i (Ji/M_i) d_chem_i/x_i - q_reduced grad(T)/T²
 = sigma_pair + lambda [grad(T)/T]² >= 0`.
This is the Onsager Schur-complement construction (Van-Brunt et al.,
<https://arxiv.org/pdf/2109.05082>, equations 3.5–3.9). The existing parabolic
engineering scale includes both thermal/composition coupling scales and their
product; no variable-coefficient spectral or fully discrete entropy theorem is
claimed. Bounded scratch remains allocation-free on successful numerical faces.
Focused fixtures independently reconstruct unequal-mass Soret-only, Dufour-only
and combined fluxes and entropy; common-heat gauge/zero-coupling/spatial reduction;
source-axis/domain/cancellation refusal; and periodic composition evolution with
closed original species and total-energy balances. Authored fixture coefficients
are not empirical transport validation.

### Finite material gas inventories with thermal cross diffusion

With `material` and `experimental-euler`,
`material::transport::parcel::ParcelMolecularTransport<L>` consumes the same
prepared `PlanarTransport<L>`, including `MaterialSoretDufourTransport`, on
explicit continuum links between rigid `ClosedGasParcel` inventories. Every
link uses zero barycentric bulk velocity and the original canonical mass-frame
species flux. These are separately exposed continuum paths; they do not add a
second molecular law to the same dusty-gas pore or supply thermal transpiration.
The original caloric species must match exactly, including initially absent
species when they enter a receiving parcel.

Each distinct participating original and proposed endpoint state is qualified
through the shared zero-gradient face kernel, in addition to every actual face
query. The kinetic, material and thermal-cross source domains therefore apply
jointly to the admitted states. The law retains its last successful query,
which may be an endpoint qualification. Per-link reports separately retain the
actual transported species, face enthalpy temperature and reduced heat flux.

All continuum fluxes and native `GasTransportPaths` sample one unchanged
inventory and enter the existing material snapshot reconstruction together.
Each species carries its original donor partial enthalpy once. Added energy is
`q_reduced*A*dt + sum_i dm_i [h_i(T_face)-h_i(T_donor)]`, with the differences
evaluated by the original analytic caloric owner before rounding. Neither large
formation references nor kinetic heat-of-transport datums become extra heat.
The native combined work bound, gross original-species withdrawal, source catalog
join, fixed absolute mass/energy budgets and closed endpoint entropy gate remain
authoritative. Source refusal, cancellation or an overlarge explicit step returns
no candidate state. This is a first-order physical candidate operator; it adds no
separate physical solver. Its whole-window accuracy owner is described below.

`tests/material_parcel_transport.rs` covers independent binary Soret/Dufour
temperature and inventory predictions, unequal-mass frame/orientation,
formation/kinetic-datum invariance, shared nozzle withdrawal, absent species,
initial and final source-domain refusal, and cancellation. Its coefficients are
manufactured integration fixtures, not measured-material qualification.

`ParcelMolecularTransport::{advance_candidate,advance_to}` completes a positive
`[start,target]` interval using bounded full-versus-two-substep comparisons of
those same native physical candidates. Only the fine path enters accepted
states and accounts. Per-parcel species inventory (kg) and temperature (K)
errors use absolute tolerances plus a relative allowance on actual changes,
so large carriers or absolute temperatures cannot conceal a small response.
The fine first-order estimate includes the correction for unequal represented
half intervals. Endpoint clocks define every physical duration; the target is
exact, a shorter final remainder still needs accuracy, and an unresolved
midpoint refuses. An adjacent representable final tail is coalesced into the
proposed span before its physical probes, never by snapping an already solved
clock. The estimate is local and asymptotic, not a global error bound.

The evaluation limit counts every full/fine probe, including discarded and
failed probes. Only typed original withdrawal, pressure-order, entropy or
caloric bracket/root-convergence refusals, and temporal error excess, can reduce
the interval. Original temperatures must lie inside the supplied caloric
brackets. Actual material, kinetic and thermal-cross source/domain refusals are
hard failures. Every physical trial retains initial/final source qualification,
the original enthalpy convention and simultaneous native `GasTransportPaths`.
No extrapolated state, source substitution or clipped inventory is introduced.

Accepted gross species streams and donor partial enthalpies reuse the native
transfer account; additional internal energy retains both directions on each
canonical parcel pair. Rejected probes contribute only to work counts. Detailed
source reports are retained for the final accepted substep. The original material
`audit_closed_window` independently evaluates actual original-to-final total and
source-qualified species mass, NASA internal-energy increments with common
species references cancelled, and caloric/mixing entropy. It and the compensated
sum of accepted native energy residuals must fit the unchanged absolute physical
budgets, with no multiplication by step count. Late refusal or cancellation
publishes neither a gas-map prefix nor an advanced caller clock.

The public adaptive fixtures exercise an independent exact coupled two-mode
Soret/Dufour relaxation, separate species/temperature refinement, stiff pure-gas
Fourier relaxation, accepted-only replay with native nozzle paths, reference
invariance, typed physical retries, late source/work/cancellation rollback and
ordinary plus large represented clocks.

`ParcelMolecularTransport::reservoir_candidate` connects independently exposed
held material baths to finite gas traces through the original zero-barycentric
Maxwell-Stefan/pressure/Soret-Dufour/Fourier face. Every bath reads one original
finite inventory. Signed species mass, combined gross withdrawal and separate
reduced heat enter the existing open species/enthalpy inversion once. Both
stream directions carry original bath partial enthalpy via analytic increments;
there is no fictitious finite bath volume or additional caloric solver.

Original bath, gas, sampled face and final gas source domains qualify. Nonzero
exchange at zero bath species activity refuses without introducing a trace.
Acceptance uses actual endpoint mass, source-qualified species, reference-free
energy residual and gas-plus-bath entropy against fixed absolute budgets. Gas
entropy may decrease when the enclosing bath account compensates it. The
candidate remains first order and requires an admissible caller-supplied step;
it does not add bulk advection, finite bath depletion or an adaptive timeline.
Public fixtures check independent flux/temperature/entropy equations, unequal
molar masses, datum invariance, pure heat, absent species, simultaneous baths,
supported entropy-decreasing overshoot, and late source/cancellation refusal.

`ParcelMolecularTransport::network_candidate` combines internal continuum links
and held-bath faces on that same original snapshot. Its internal face sampler is
shared with the closed finite-gas consumer. All signed species streams, original
donor/midpoint enthalpy corrections and separate bath reduced heat enter one
native gross-withdrawal budget, source join and caloric reconstruction. Bath
arrivals cannot replenish an internal outflow allowance. Internal transfers
cancel in the final open mass/energy/entropy audit; only actual bath streams
cross its accounting boundary. All participating final gas states, including
cells with no direct bath contact, must satisfy the original constitutive law.

Both face families share globally unique physical IDs and one native work bound
before rate allocation. Reports retain internal and bath faces separately;
`reservoir_candidate` delegates with no internal links. Empty bath/internal
families preserve the corresponding physical specialization. This remains a
first-order continuum candidate with an explicit duration. Public network
fixtures cover independent two-cell/bath Fourier balances, simultaneous binary
species/Dufour increments, closed/open parity, shared gross withdrawal, joint
work/identity limits and late source/cancellation atomicity.

`advance_network_candidate` and `advance_network_to` evolve that simultaneous
network across a complete represented time interval. They share the existing
full-step/two-half-step controller with closed molecular transport. The fine
trajectory alone enters accepted states and accounts; no extrapolated state or
rejected probe enters a physical ledger. Local error scales separately bound
each species inventory, temperature, directional species stream and additional
heat stream, so throughflow and reversals cannot hide behind a small net change.
The stream-relative allowance scales transported amounts, not carrier inventory.

Accepted bath streams retain compensated signed sums and separate gross inward
and outward totals. The final original-to-final audit uses unchanged original
bath models, partial enthalpies/entropies and integrated reduced heat; changing
per-step relative caloric charts are not integrated as external energy. Internal
streams cancel at this boundary. Accumulated native caloric defects and actual
endpoint mass/species/energy/entropy each keep one fixed whole-window budget.
Retained accounts share the native work bound, and graph size is admitted before
allocating the held-bath map. Typed physical/time-error refusals refine the step;
source-domain failures and cancellation remain terminal and atomic with the
caller clock. These are local temporal estimates, not certified global errors.
Public fixtures compare two independent thermal modes and a direct bath-heat
integral, isolate stream accuracy from state accuracy, retain balanced-bath gross
species/heat and mixing entropy, and check face-order, work, source, cancellation
and representable-clock behavior. Finite reservoir depletion and native bulk or
dusty-gas path composition remain outside this continuum-plus-held-bath API.


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

`TwoZone::apply_split` composes a supplied heat/stream source substep with an
already timed outer update. It uses the same atomic inventory, caloric,
pressure and interface-work solve as `advance`, then restores the exact
original clock after success. Only the enclosing time owner advances the
physical interval and records those source rates once. A focused regression
compares its complete physical state and formation-energy-aware stream balance
with an ordinary source advance while checking the retained clock and atomic
cancellation. This operation does not add chemistry or a higher-order split.

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

The adapter returns chemistry-only species rates and the source-coefficient
energy-reference correction `-extent_rate * Delta u(Tmin)`, evaluated at the
original `MixtureGas` reference temperature. Conserved elemental formation
references cancel before absolute-energy rounding; distinct formation energies,
original caloric corrections and molecular-count pressure work remain active.
Reported reaction enthalpy and power are diagnostics, never a second heat
source. Power-to-extent conversion uses source-coefficient Delta h at an explicit
reference temperature. `chemical_reference_roundoff_w` separately exposes the
discrepancy from independently rounded species-reference powers and is never
added as heat. The original moving-interface solve retains all paired p dV work. One common extent may be rounded down by a bounded
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
stresses or rarefied gas. Thermal cross coupling requires the explicit source
law above; a fully discrete entropy or invariant-domain theorem remains outside
this operator. Native source-linked
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

Opt-in `WallDrag::MovingLinear` and `WallDrag::MovingDarcy` prescribe constant
axial wall velocity U over a source interval at fixed gas volume. They exactly
relax relative velocity u-U using the same frozen linear/Darcy laws. Gas receives
mechanical work U times its actual accepted momentum change; the wall/drive
receives the opposite impulse and work. Loss of relative-slip kinetic energy
becomes sensible gas energy once, even during entrainment or braking.
`mechanical_work_to_gas_j` remains separate from `heat_to_gas_j` and
`drag_heating_j`, and frozen branch reports retain `wall_work` separately from
`wall_heat`. Both frozen and reacting/molecular mixture network balances include
drive work. Original caloric wall storage and material-derived coefficients use
the same owner; coefficient refresh preserves wall motion and admits Reynolds
and Mach numbers in the wall frame. Stationary cases retain their original
arithmetic. This is a prescribed sliding-wall source, not resolved wall inertia,
moving geometry or an unsteady boundary-layer model; zero-momentum plenums
continue to refuse drag. Four numerical-kernel and eight native integration
regressions cover analytic slip relaxation, frame covariance, weak/strong drag,
finite calorics, network composition and transactional refusal.

`WallDrag::LinearQuadratic` applies authored linear and inertial resistance
simultaneously to slip s = u-U: ds/dt = -a*s - b*abs(s)*s, with nonnegative
a [1/s], b [1/m] and a prescribed constant wall velocity U. Its single exact
frozen-coefficient update is s(t) = s(0)*exp(-a*t) /
[1 + b*abs(s(0))*(1-exp(-a*t))/a], using the continuous a=0 limit. Stable
exponential differences retain weak impulses; scaled rational decay admits
strong resistance without splitting the two forces or clipping the residual
slip. The existing accepted-impulse work and slip-heating account is shared.
The zero-coefficient limits are the original linear and quadratic laws.
G1 regressions compare an independent reciprocal-slip solution, subdivision,
frame transformations and weak/strong limits; the mixture-network source test
also checks actual species, temperature, support impulse and drive-work balances.
No coefficient, transition law or applicability range is inferred.

`exchange_cell_with_recoil` couples one original mixture cell to a finite
`SlidingSupport` on the same axis at fixed gas volume and species mass. The
existing exact linear/quadratic slip map uses the reduced mass m*M/(m+M), so
both bodies respond to the same interval. Accepted gas impulse I gives the
support -I, wall work I*(U0-I/(2*M)), and nonnegative slip loss
-I*(s0+I/(2*reduced_mass)). An authored `gas_heat_fraction` partitions this
loss; a nonzero wall share requires finite `WallThermal::Lumped` storage.
The original caloric heat owner then advances the contact without resetting a
flow clock. Returned gas/support endpoints and actual finite-wall storage must
pass independent momentum and complete energy budgets before return.
`RecoilUpdate` exposes both slip-heat shares; the original `WallCellReport`
retains total slip loss and accounts for the wall share as heat leaving gas,
so its sensible-energy identity remains additive. Inputs are immutable on
refusal or cancellation. This owner binds one support to one gas amount;
shared multi-cell inertia and moving geometry require their own coupled owner.
G1 regressions cover reduced-mass analytic decay, subdivision, large mass ratios,
frame covariance, heat fractions 0/0.5/1 with subsequent conduction, actual
combined storage, and late cancellation. Friction coefficients and the heat
fraction are authored; the thermal step retains its original first-order law.

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

### Current material state in developed pipe-wall losses

`NetworkWallOperator::with_circular_laminar_transport` binds the existing
`TransportLaw` (including source-card `MaterialTransport` and original kinetic
transport) to selected circular branches. Full caloric/species-source identity
must match the network. At every wall source half it obtains temperature,
pressure, density and composition from the actual conservative cell, reuses the
original molecular mole-fraction closure, and queries the same source-owned
viscosity and conductivity. Thus new reaction products, changed inlet mixtures,
pressure or thermal evolution change wall drag and cooling without a caller
re-entering property scalars or a separate reacting-wall solver.

The declared `CircularLaminarWall` supplies diameter, Nusselt number, supported
Reynolds interval (within 0..2300), and maximum Mach (at most 0.3). Diameter and
cell volumes must agree with the original constant-area circular mesh. Drag is
`32 mu/(rho D²)`; complete cell conductance is `Nu k (4 V/D)/D`, where `4 V/D`
is that cell's wetted area. The Nusselt number and fully developed applicability
are explicit engineering inputs, not values inferred from a gas name. This
closure does not establish entrance length, negligible axial conduction,
high-frequency oscillatory boundary layers, turbulent flow or a measured pipe
correlation. Source coefficients are evaluated at bulk gas temperature and
frozen during each source half; that constitutive lag remains first order.

Both source-step endpoints must satisfy the original property law and declared
Reynolds/Mach limits. A changed product outside its material-card domain rejects
the complete trial, including earlier cells, wall temperatures, reaction updates
and reports. The existing exponential drag, implicit finite-wall thermal solve,
species conservation, support impulse and energy accounting remain unchanged.
No molecular or chemical heat source is added twice. Finite wall heat capacity
is still an explicitly authored constant; this binding updates gas-side
transport, not solid-wall calorics. Unbound branches and stationary plenums
retain the existing explicitly authored contacts. Original law workspaces may
change during refused trials, but carry no accepted physical state.

Focused G1/G3 regressions check independent drag/finite-wall heat equations,
density and transport-data substitution, exact caloric-source and geometry
admission, heated-endpoint rollback, and positive reaction-product refresh
through `ReactingNetwork1d::step_with_wall`. An original material-card refusal
after reaction also checks whole-network/finite-wall rollback. These are
manufactured numerical integration tests, not experimental combustion or
convection validation.

### Existing source composition

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

## Three-dimensional original-mixture Riemann and wall closures

`multicomponent::SpatialCell` and `SpatialPrimitive` extend the original
mixture owner to an orthonormal normal/two-tangent face frame. Canonical
species, all three momenta and full species-referenced total energy use the
same const-generic HLLC/HLL implementation as axial and planar consumers.
Star and HLL source admission includes both transverse kinetic contributions;
mixing their shear energy cannot escape the original caloric interval. The
original planar arithmetic and source-owned shock/rarefaction wall pressure
are retained. Reconstructed primitives retain source-computed temperature,
heat capacity and sound speed without another energy inversion.

Successful face calls use bounded stack scratch and caller output. Source and
precision refusals never clip temperature or renormalize species. Stationary
adiabatic slip walls have exactly zero species, tangential momentum and energy
flux and use the incident cell's original composition and calorics. Eight
focused tests cover bitwise planar reduction, material contacts with both
shears, tangential rotation/normal reversal/Galilean covariance, combined-shear
source refusal, variable-Cp and formation-gauge retention, original wall
pressure, precision/cancellation and the full 128-species bound.

This remains the Estimated conservative engineering lane. No moving
unequal-caloric interface pressure-equilibrium property, fully discrete
entropy/invariant-domain theorem, chemistry, spatial convergence or molecular
shock-thickness validation follows from a successful face evaluation.

## Complete methane chemistry with trace-preserving local extents

The existing prepared reactor now admits the complete 53-species, 325-channel
`kinetics::gri30::cantera_gri30` source, including its 16 originally irreversible
channels. `ContinuousAt300K` uses the source's explicitly declared upper-region
enthalpy/entropy continuity adjustment, retaining every original Cp law and
each species' own NASA seam. Its common admitted interval is 300--3000 K.
No reduced global methane reaction or added combustion heat source is used.

Reaction coordinates are rebased only after an accepted complete step.
Candidates reconstruct from that accepted composition plus local stoichiometric
extents; the original feed and fixed total internal energy remain untouched for
whole-call mass, elemental, charge, energy and entropy audits. This prevents a trace
intermediate from becoming the subtraction of large accumulated reaction-cycle
extents. It introduces no species clipping, positive trace floor, radical seed
or equilibrium projection, and rejected attempts retain their accepted origin.

The smaller Newton system remains algebraically equivalent to the full extent
equation. Exact nonzero stoichiometric rows and columns are prepared once,
without thresholds or rank truncation; compensated sums retain source order.
For GRI, the common factorization is 54 by 54 instead of 325 by 325. Its
stoichiometric matrix has 1227 nonzero entries out of 17225, reducing reduced
assembly product terms from 947700 to 83808, including the unchanged thermal
row. A focused release measurement of 40 identical source-derived equations
took 13.151 ms with the sparse traversal versus 90.728 ms with the dense
traversal; their 325 corrections were bit-identical. These are one-run
operation timings, not a whole-reactor or portable speed guarantee.

Every reduced correction must satisfy the original unscaled componentwise
linear backward-error gate. Up to four residual corrections are bounded.
If trace-channel cancellation still fails that gate, prepared storage permits
the existing scaled-partial-pivot full-extent solver, followed by the same gate
and at most two direct residual corrections. Physical and numerical acceptance
tolerances are unchanged. `IsochoricKineticsBalance::direct_fallbacks` reports
successful recovered Newton corrections, including those in later-rejected
time-step attempts. For GRI the fallback matrix/RHS/scales reserve 850200
bytes; the sparse column cache adds 22240 bytes beyond the sparse-row path.
All numerical buffers are prepared before advancing; no hot-path allocation
is introduced.

The independent development-only reproducer
`data/reference/thermochem/gri30_cantera32.py` verifies Cantera 3.2.0 and the
pinned source blob, then integrates the complete mechanism separately with
CVODES. Two native SDIRK2 regressions use CH4:O2:N2 = 1:4:15.04, one litre,
101325 Pa and no radical seeds. At 2 ms, the 1300 K induction case agrees with
the 1300.5483131534243 K reference within 0.00000181 K; the 1500 K ignition
case agrees with the 2600.310562321808 K reference within 0.0000561 K.
Both compare all 53 species within 0.5% plus 1e-12 mol, original physical
energy within 2e-8 J, pressure, independent atom/mass inventories and entropy.

The hot case uses 911 accepted steps, 154 rejected attempts, 11997 Newton
iterations and 1296 successful direct fallbacks at relative tolerance 2e-5.
Its observed mass residual is 1.11e-19 kg, maximum relative elemental residual
1.19e-15 and energy residual magnitude 2.85e-14 J. The induction case uses
110 accepted steps and no rejections. These are fixture-specific work and
accuracy observations, not certified error bounds. All 25 focused chemistry
tests pass, including original hydrogen, integration-order, conservation,
datum-invariance, cancellation and sparse/dense-equivalence regressions. The
complete focused native gas library also passes all 158 tests.

This extends homogeneous methane/natural-gas/nitrogen chemistry. It supplies
neither a fitted wood/PMMA pyrolysate nor soot, char oxidation, spray, real-gas
or experimentally validated flame/emission prediction. A finite-step entropy
check does not turn the source's irreversible channels into detailed-balanced
reactions, establish a fully discrete coupled entropy theorem, or give higher
order to surrounding pressure, transport, radiation or solid-source splitting.

## Shared original-mixture planar molecular faces

`multicomponent::dissipation::PlanarTransport` prepares the original
Maxwell--Stefan workspace and transport law on the exact evolving species
bank. Its shared private face kernel also serves the original axial and gas
network operators; there is no second species conversion, diffusion solver or
enthalpy convention. A planar face transports every diffusive species flux,
both momentum components and full referenced energy together. Its Newtonian
stress includes normal and transverse derivatives under the three-dimensional
trace convention with zero out-of-plane velocity. Energy includes partial
referenced species enthalpy, Fourier heat and both components of stress work.
No separate chemical or viscous heat source duplicates those transfers.

Conservative inputs use the original caloric inversion. Already reconstructed
primitive inputs use the original source admission without another inversion.
Stationary adiabatic slip walls retain normal viscous traction and exactly zero
species, shear and energy flux. Output is scratch until success; callers own
spatial gradients, conservative divergence, requested-duration transaction and
the complete boundary impulse. The reported diffusion scale is a conservative
engineering step-size estimate, not a spectral or invariant-domain certificate.

Prepared numerical capacity includes the binary matrix and all original
Maxwell--Stefan solve buffers. Successful kernels allocate no new scratch.
Immutable gas/species storage, bounded stack scratch and the caller's transport
law storage/evaluation allocations are outside that budget. Material-card
providers retain their original source-query allocation and validity contract.

Five focused regressions cover independent binary diffusion with unequal Cp,
Fourier/full-tensor/enthalpy flux, primitive/conservative agreement, tangential
frame and formation-gauge covariance, exact axial reduction, slip-wall loads,
exact-zero species, source/axis refusal and final cancellation/replay. Original
axial/network regressions remain authoritative for their geometry and balances.
Soret/Dufour coupling uses the explicit source-law extension above. These face
fixtures establish no nonideal transport, spatial accuracy or combustion
validation. Explicit held
reservoir faces use the source-owned boundary energy convention below.

### Shared original-mixture molecular faces in three dimensions

`multicomponent::dissipation::SpatialTransport` extends the same prepared
source and const-generic molecular kernel to all three momentum components.
The original axial and planar owners retain their public interfaces and use
that same kernel with zero or one transverse component; no diffusion solver,
species conversion, caloric reference or constitutive law is duplicated.
Three-dimensional conservative inputs remove all kinetic components through
the original source inversion. Primitive inputs use the original complete
three-momentum admission without repeating the temperature solve.
Every face receipt retains `dynamic_viscosity_pa_s` from that original law
evaluation. A composing SGS owner can check its eddy/molecular ratio against
the actual face coefficient without querying a stateful source a second time.

Each face supplies its positive normal separation and both tangential rows of
the velocity gradient in the original orthonormal face frame. The Newtonian
traction includes both transverse dilatations in the normal component and
both normal/tangential derivatives in each shear component. All three stress
work terms accompany original referenced species enthalpy diffusion and
Fourier heat in one total-energy flux. Stationary impermeable adiabatic slip
walls retain both transverse dilatations in the normal viscous load while
returning exactly zero species, shear and energy transfer. A field owner must
still supply physical gradients, divergence, time integration and complete
boundary exchanges.

The prepared binary matrix and Maxwell--Stefan numerical heap capacities are
unchanged by the extra momentum direction. Successful kernels allocate no
additional scratch; caller-supplied law allocations keep their own contract.
The reported engineering diffusivity estimate includes the complete
three-dimensional Newtonian row couplings and must be combined with all
spatial inverse lengths by the field owner. It is not a spectral,
invariant-domain or fully discrete entropy certificate.

Six focused regressions independently reconstruct binary species/Fourier/full
stress-work fluxes, compare conservative and primitive inputs, rotate both
tangent axes and the complete derivative tensor, check Galilean and formation
reference covariance, retain the planar reduction bit for bit, and verify
slip-wall loads, a resolved `1e-180` trace species, an exact absent species,
source/axis/gradient/caloric refusal and final cancellation/replay. They do
not establish spatial convergence, physically resolved shock/flame thickness,
Soret/Dufour effects, no-slip walls or experimental combustion validity.

### Original-source energy effort at physical reservoir faces

Both planar and spatial face owners expose `reservoir_face_from_primitives`
for an explicitly selected exterior trace at a physical boundary. The original
source admits both traces. Maxwell--Stefan, traction and Fourier transport
retain the shared two-trace midpoint approximation, while species energy uses
the bath's referenced enthalpies and stress work uses the bath's velocity.
A stationary bath therefore receives zero mechanical work even with nonzero
viscous traction. The supplied distance is cell-center to physical boundary,
not the spacing between two cells. The existing three-dimensional gradient
layout and all interior-face source conventions remain unchanged.

The kernel does not evolve the exterior reservoir. The field owner must retain
that state, its geometry and signed boundary transfers consistently with its
advective boundary condition. A focused regression checks both exterior-side
orientations, original unequal-Cp bath enthalpy with large formation datums,
all three work components, unchanged species/traction fluxes, and zero work at
a stationary bath in both dimensions. Spatial accuracy and finite-reservoir
storage require their respective field and coupling owners.

### Three-momentum chemistry on the same prepared cell source

`CellChemistry::advance_spatial` extends the same prepared axial/planar owner
to canonical species, px/py/pz and referenced total energy. All three momentum
components remain fixed, and all three kinetic terms enter the original-source
internal-energy target and independent endpoint reconstruction. The original
formation-reference change is applied once, with no extra heat-release source
on chemical total energy. Dimension-specific public entry points share one
const-generic conversion and acceptance kernel.

The existing atomic cell/report publication, exact source identity, caloric
interval and bounded-work requirements remain unchanged. Actual prepared
reactor capacities are available to the enclosing spatial memory budget.
Focused dimension-reduction/rotation and late-cancellation/replay checks cover
complete chemical energy, all momenta and reusable prepared storage. This is
the local fixed-volume chemical source; it does not itself establish a spatial
entropy theorem, resolved flame structure or experimental combustion validity.

### Compatible composition through density and compression waves

`multicomponent::reconstruction::mass_fraction_slopes` supplies the same
bounded, allocation-free composition limiter to axial and Cartesian mixture
consumers. Traces reconstruct density and every canonical mass fraction
separately, then recover partial densities as their product. Conservative
cell averages and their original species equations are unchanged. A single
MC factor scales the entire centered fraction-slope vector; constant inert
fractions and other linear composition constraints survive density variation
to floating-point roundoff. No dependent remainder species, normalization,
trace seed or independent species correction is introduced.

The relative-only 32-epsilon roundoff guard prevents division noise in an
otherwise constant fraction from flattening all composition slopes, while
retaining resolved gradients even for a 1e-200 trace constituent and preserving
exact absent species. The axial owner uses actual nonuniform cell-centre
distances. Its source-domain contraction scales density, every fraction,
velocity and pressure together before the original thermodynamic admission.
The additional density reconstruction scalar is prepared once per cell and
does not add a second conserved density equation or per-step allocation.

Focused compression-wave tests on uniform and nonuniform 80-cell meshes
retain a uniform inert fraction to 2.23e-16 and an independent weighted
composition constraint to 4.45e-16 while closing every species and energy
balance. These verify compatible reconstruction; they do not establish a
fully discrete composition-maximum principle or entropy theorem.


## Quasistatic moving compact gas volumes

\`multicomponent::network::moving\` adds the first explicit moving-boundary
thermodynamic owner for species-bearing compact gas storage. A
\`Plenum::quasistatic_adiabatic_to\` candidate preserves every canonical species
mass while changing physical volume. It solves the original fixed-composition
NASA-9 ideal-gas isentrope, \`Delta s(T,p)=0\` with \`p=rho R_mix T\`, through
the same borrowed caloric source used by the wall/characteristic owners.
Safeguarded Newton/bisection is bounded by the existing caloric-domain work and
temperature limits; unsupported compression/expansion refuses rather than
clipping temperature or pressure.

The receipt reports source temperature/pressure endpoints, mechanical boundary
work on the gas, physical internal-energy change, extensive entropy defect,
per-species mass residual and the referenced-versus-absolute energy-gauge
residual. Because the closed species inventories are unchanged, formation
references cancel from work exactly in the continuum. Boundary work is not
reclassified as heat and no \`p*Delta V\` endpoint approximation is used.

`MixtureNetwork1d::quasistatic_adiabatic_plenum_to` accepts that same endpoint
directly into a running connected chamber between transport steps. Every branch
field, other plenum, connection and clock stays unchanged. The next transport
step reads the accepted chamber volume and source state for both port fluxes
and its acoustic CFL restriction. A single branch-field pass admits old/new
combined inventories before publication, including refusal if total network
storage would overflow. The returned `NetworkVolumeReport` separates mechanical
work from transport/heat and reports the combined species and energy balance.
An actual accepted volume change invalidates the obsolete flow/port reports;
an exact no-op, refusal or cancellation retains them. No network reconstruction
or source-gauge conversion is required.

These operations are deliberately clock-free and quasistatic. Network use
holds ports closed during the volume change and leaves the schedule and
accuracy of composing that change with transport to its caller. They do not
implement a finite-speed moving-wall Riemann solver. Piston inertia, prescribed
motion, shock generation during boundary motion, heat leakage and moving
geometry/GCL remain explicit separate owners. Three focused endpoint tests cover an independent
constant-Cp isentrope and reversible round trip, variable calorics with +/-1e18
formation-energy gauges, and source-domain/cancellation atomicity. Three
additional network regressions cover compression-driven branch flow with an
independently integrated physical work balance, reversible volume cycles, and
state/clock/report preservation under early or late cancellation, source
refusal and exact no-op.

## Finite-time compression coupled to network transport

`MixtureNetwork1d::step_with_volume_rates` evolves prescribed per-plenum dV/dt
simultaneously with branch and characteristic-port flow through the original
SSP RK2 stages. Species inventory and referenced internal energy advance
extensively, with source-owned pressure work `-p dV/dt` accounted separately
from transported energy. Accepted endpoint volumes and pressure work are
available per plenum, and `NetworkStepReport::boundary_work_on_gas_j` closes
the combined network energy balance with that same pressure quadrature.

Positive representable volumes and acoustic/volume-rate CFL bounds gate atomic
publication. All-zero rates preserve the rigid-step arithmetic; failure retains
fields, volumes, clocks and accepted receipts. Work/report buffers are allocated
at construction and reused. G1 tests compare 20/40/80-step compression and
expansion ladders with independent constant/linear-Cp isentropes. Other focused
tests check simultaneous port enthalpy/species/work closure, exact rigid
identity, geometry/source-domain refusal and early/intermediate/final cancellation.

The compact perfectly mixed chamber approximation remains. Rates are constant
over each accepted step. Spatial piston shocks, branch-wall motion and piston
inertia require other physical owners; this is not a spatial ALE solver.

### Moving chambers with reacting source operators

`ReactingNetwork1d` composes prescribed plenum volume rates with finite-rate
chemistry, molecular species/heat transport, viscosity and finite-wall storage.
The existing split is W/2 -> V/2 -> D/2 -> C/2 -> moving transport -> C/2
-> D/2 -> V/2 -> W/2, omitting absent operators. Source halves before transport use the
original volume; halves afterward use the actual accepted endpoint volume.
Boundary work has its own report field and is subtracted exactly once from
referenced and physical-energy residuals.

Prepared wall, molecular and valve operators keep accepted geometry separate
from unpublished interval volumes. Complete-step acceptance publishes the new
bindings, while shortened-duration retries restore original fields AND volumes.
Final cancellation preserves accepted geometry, reports, valve exchanges and
wall storage. All-zero rates retain the original rigid path, including valves.

Valve halves use the actual original/final chamber volumes in the unchanged
source-owned nozzle and finite donor-blowdown solve. Each half holds its stage
volumes fixed and transfers identical species and finite enthalpy to the
receiver; only moving transport supplies mechanical pressure work. Source
receipts sum the two valve intervals' storage defects, excluding intervening
compression, chemistry and heat. Public standalone valve calls still require
the accepted geometry, and failed trials cannot rebind it.

Authored wall conductances, contact areas and compact mixing lengths remain
fixed inputs; a volume rate does not determine their changing geometry.
The compact mixing and existing reaction/source-splitting approximations remain.
Focused regressions exercise live reaction/transport/compression/wall energy,
duration restart, final cancellation and reuse, accepted geometry rebinding and
exact zero-motion identity. Moving-valve cases additionally reconstruct each
tank's species and energy change from finite enthalpy exchanges and pressure
work, close physical gas-plus-wall storage with all sources active, and check
complete-duration retries and cancellation/restart. These new valve-composition
regressions are added but have not been executed in this resource-constrained
session. No new global temporal-order or entropy-stability claim is made.

## Original material calorics in finite gas walls

`MaterialWallCaloric` binds the existing `PreparedHeatStorage` to finite
branch-cell and plenum walls. Its explicit caloric constraint declares
fixed-volume internal-energy capacity or the stationary-solid
negligible-expansion approximation; Cp is never silently converted to Cv.
Original source lines, fixed constituent amounts and temperature support remain
attached. Source integrals determine actual storage; a temperature contact's
capacity is only its refreshed tangent cache.

The same binding's `from_thermal_law` retains an original
`fs-thermal::ThermalLaw`, its physical reference cell, fixed mass and explicit
supported energy interval. Both owned phase charts and adapters borrowing
prepared material sources can be bound without copying a fitted heat capacity.
`NetworkWallOperator::with_material_heat_storage` supplies either source to the existing
ordinary, molecular, reacting and valve-composed network transactions.
Developed pipe-wall transport continues to query current gas composition,
viscosity and conductivity; it replaces only drag/conductance, preserving the
original wall source and evolving stored energy.

The source step solves one bracketed backward-Euler equation in transferred heat:
`Q = dt G [T_wall(e0-Q/m)-T_gas(E0+Q/V)]`. The original gas calorics recover gas
temperature and the original wall law supplies T(h/u), including exact latent
plateaus. Monotonicity belongs to the fixed-conductance case; a source-dependent
endpoint conductance makes no uniqueness claim. The same safeguarded scalar owner uses exact source tangents when
available and bisection otherwise. No apparent capacity or separate latent-heat
source is added. Actual
represented gas energy, wall storage and heat-flux residuals must all meet the
existing joule allowance. Source/geometry/range/work/cancellation failure returns
no candidate; the original complete-network transaction publishes every gas and
wall state together. Source storage is immutable and numerical work is reused.

The law's storage basis remains explicit: internal energy describes an admitted
fixed-volume source; enthalpy retains the material's own pressure-work convention.
The combined report sums the selected thermal coordinates and does not relabel
wall enthalpy as rigid-vessel internal energy. Reference heat-transfer geometry,
uniform-wall applicability, imposed conductance and any external pressure work
remain caller/source responsibilities. Melting-wall motion, loss of mechanical
support, wall sharing between cells, nucleation, and spatial wall conduction are
not inferred by this finite-store coupling.

Focused tests check independent nonlinear calorimetric references, original
source-only substitutions, reaction/domain/cancellation rollback, latent storage
and full plateau traversal in both heat directions, accepted-history rollback
and replay, continuing molecular/wave energy conservation, and the original
borrowed solid material chart in a reacting network with current gas-card losses.
Native execution remains unverified in the authoring worker: the unchanged
asupersync dependency was SIGKILLed during code generation under its memory cap.

Four further G1/G3/G4 cases exercise implicit material conductivity against
independent latent-plateau quadratic roots in both heat directions, a valid root
beyond the frozen-flux bound, variable-capacity endpoint energy, source-only
changes, constant-conductance and disabled-film limits, continuation and
source/work/cancellation refusal. The reacting-network regression runs both time
evaluations. These new Rust cases remain unexecuted; six changed Rust files
passed pinned syntax parsing and scoped whitespace checks.


## Pressure-actuated network valves

`Orifice::with_pressure_opening(PressureOpening)` adds a reversible quasi-static
differential-pressure aperture to the existing finite-plenum valve owner.
Finite typed thresholds satisfy `0 <= cracking < full_open`: aperture is zero
through cracking, linear between thresholds, and capped at one. It multiplies
the declared maximum `Cd * opening * area` exactly once. Existing
`ValveDirection` enforces exact reverse closure.

Initial and implicit endpoint evaluations use current source-decoded pressures.
The actual nozzle pressure ratio, source-owned variable calorics and choking,
adiabatic donor blowdown, and equal-and-opposite species and energy transfers
remain authoritative. Cracking pressure is never subtracted from nozzle
backpressure. Closed endpoint trials bracket the transfer instead of accepting
threshold overshoot under a loose absolute mass tolerance. Invalid thresholds,
nonfinite differences, and unresolved nonzero aperture or area refuse.

This is an Estimated algebraic aperture law without hysteresis, poppet inertia,
spring-energy storage, leakage, calibration, or a certified safety-relief claim.
Fixed-opening construction retains its previous behavior. Focused
`pressure_actuated` regressions exercise analytical choking/un-choking,
thresholds and directions, finite-tank blowdown and energy/species balance,
source rollback, and pressure-wave propagation in the original gas network.
Their presence is not an executed validation claim.


## Prescribed species-bearing pressure reservoirs

`MixtureEuler1d::step_with_reservoirs` accepts stationary species-bearing end
cells in the same source species order and per-species energy-reference gauges.
It recomputes the existing characteristic fluxes from the actual interior cell
at both SSP RK stages and every halved retry. Reservoir inflow carries the
source composition and stagnation enthalpy; outflow carries the evolving
interior composition and total enthalpy. Reservoir pressure supplies inflow
stagnation pressure or outflow backpressure through the original variable-
caloric shock/fan and choking owner.

An absent endpoint preserves its declared boundary; two absent reservoirs
execute the original step exactly. Periodic intervals with reservoirs, moving
reservoirs, invalid dimensions and source-domain violations refuse before cell
or clock publication. Borrowed `Cell` values carry no owner identity token:
callers retain responsibility for matching canonical species order and energy
gauges, as with `MixtureGas::reservoir_flux`. Stack endpoint buffers add no
per-step heap workspace. Prescribed states are fixed during one step and are
not depleted; finite-reservoir feedback remains with the existing network owner.

Six focused `reservoir_boundary_tests` regressions exercise independently
reconstructed predictor/corrector fluxes, variable-Cp donor enthalpy, pressure/
temperature/composition drives, end reversal, formation-reference invariance,
CFL retry replay, exact absent-reservoir parity and cancellation/refusal rollback.
Their presence does not establish executed native validation.


## Prescribed mixture pressure reservoirs with finite-volume motion

`MixtureNetwork1d::step_with_reservoirs` accepts one optional stationary left/right source pair per branch, in canonical species order and energy reference gauges. Raw `Cell` source compatibility remains caller-owned; species dimensions, stationarity and source-state validity are admitted. Only unconnected physical ends may have prescribed reservoirs. Characteristic fluxes refresh at both SSP RK2 stages and every retry, preserving reservoir inflow species and stagnation enthalpy, and interior outflow species and total enthalpy. External fluxes enter the global boundary balance; finite-plenum exchanges remain exact opposites. An all-None source slice keeps the original step. Fields, clocks and accepted reports publish atomically. Construction-time primitive caches and existing endpoint workspaces avoid per-step numerical heap allocation.

`step_with_reservoirs_and_volume_rates` applies prescribed exterior sources and finite-plenum motion in the same native SSP RK2 transaction. Each stage uses the evolving chamber volume, composition and energy alongside fresh branch-end characteristic states; mechanical `-p*dV/dt` updates extensive energy beside conservative port transport. Reports separate exterior enthalpy transfer from accepted plenum boundary work. All-None sources preserve the volume-rate lane; all-zero rates preserve rigid-reservoir arithmetic, and every retry replays both drives before one shared publication point. Reservoir controls are not yet propagated through reactive or wall-source composed network wrappers.


## Pressure reservoirs through coupled reacting and molecular network owners

The prescribed species-bearing pressure reservoirs of `MixtureNetwork1d` now
propagate through the complete reacting and molecular network transactions
instead of requiring a caller-side split. `ReactingNetwork1d::step_with_reservoirs`
and `step_with_reservoirs_and_sources` compose the original C/A/C chemistry
with reservoir-driven characteristic transport and optional wall/valve halves.
`step_with_reservoirs_and_volume_rates` evaluates the same exterior pressure
drive simultaneously with finite-plenum pressure work. The molecular owner
likewise exposes `DissipatingNetwork1d::step_with_reservoirs` and
`step_with_reservoirs_and_sources`, retaining D/2-A-D/2 and optional source
halves in one complete-duration retry.

Reservoir states remain prescribed and undepleted during one call. Every
transport retry recomputes their characteristic flux from the actual branch
stage state. Reservoir mass/enthalpy appears exactly once in the transport
boundary report; chemistry reference changes, molecular boundary transfer,
wall heat, valve-internal exchange and moving-volume work keep their existing
separate accounts. A shortened acoustic step restarts all chemistry, molecular,
wall and valve owners from accepted storage before replaying the pressure drive.
Cancellation or source refusal publishes none of those states or receipts.

Focused regressions drive a new species from a high-pressure reservoir through
the reacting network, including simultaneous plenum compression, and through
the molecular D/A/D owner. They reconstruct species and physical-energy changes
from the complete reports. These tests are authored but native execution remains
subject to the repository's existing build-resource limits.

## Gravity-aware finite-plenum networks

`GasNetwork1d::new_with_gravity` couples the original isothermal and isentropic
branches through one explicit `NetworkGravity` convention. Each branch retains
its signed local acceleration and reference position, with
`phi(x) = branch_offset - g*(x-x_ref)`. A branch without gravity has constant
potential. Every rigid compact plenum has one declared specific potential;
all incident physical faces must meet it. Admission subtracts branch/plenum
datums before adding the local head, with tolerance based on that physical
difference rather than a large arbitrary absolute zero. Missing declarations,
mismatched heights, non-finite or unrepresentable storage refuse.

Connected characteristic ports reconstruct the physical hydrostatic face from
the original retained reference and current cell perturbation. Volume-averaged
end-cell pressure is not substituted. Identical admitted resting traces reuse
the original source evaluation. `Euler1d::hydrostatic_reference_face` exposes
the retained face state for consistent compact-volume initialization. Each
port still evaluates one source-owned characteristic flux transferred with
opposite signs. Gravity work uses the original branch face-mass-flux identity;
no extra junction work term is introduced.

Gas reports retain branch kinetic/internal energy and compact internal energy
and now explicitly account for branch gravity work. `NetworkGravityReport`
separately balances gas plus midpoint branch/compact potential storage against
external gas energy and `phi*mass` transport. All potential terms are accumulated
after subtracting the first branch offset. The exposed reference and centered
storage/transport preserve resolved work under a large common datum; the
absolute-potential residual is `R_centered + phi_ref*R_mass`. Absolute source
formation energy additionally uses the existing `e_ref*mass` convention. Direct
subtraction of huge absolute storage levels is not used for the conditioned
residual. No per-cell arrays are allocated during stepping.

Six native regressions cover exact joined isothermal/isentropic rest with both
reconstruction orders, transient exchange among three elevations with independent
storage reconstruction, external potential-energy transport, bitwise flow and
conditioned-balance invariance under a representable `2^40 J/kg` potential shift
and `1e18` formation datum, physical height/refusal checks, and rollback at the
final publication poll. All six new tests and all 80 native dynamics tests pass.

A compact node cannot represent a hydrostatic column spanning different
elevations. Original constant-area/constant-acceleration frozen-ideal-mixture
assumptions remain. This adds neither arbitrary entropy stratification, gravity
for reacting or diffusive networks, moving nodes, self-gravity, nor an
entropy/positivity certificate.

## Coupled atomic ionization and collisional thermal evolution

`chemistry::ionization::advance_atomic_ionization_with_collisions` evolves the
original ground-state `e + A <=> A+ + 2e` channel and electron/heavy thermal
exchange simultaneously in the existing rigid-parcel owner. It consumes
`energy::ElectronHeavyCollisionSources`: original neutral momentum cross
sections, optional Coulomb model, and explicit charged-fraction, massive-target
and Maxwellian quadrature controls. The original Coulomb receipt must match
the initial population. Later evaluations explicitly refresh the same source
law at the current amounts, density, and separate temperatures. Recombination
requires the neutral source even when that neutral is initially absent.

For reaction extent xi and cumulative electron-to-heavy transfer Q, the
electron budget is `Cv_e*(ne0+xi)*Te=Cv_e*ne0*Te0+Pe*t-chi*xi-Q`, with
`Cv_e=3R/2`. The heavy source-difference caloric budget is `H(Th)=Ph*t+Q`.
Thus chemical storage, changing electron heat capacity, external work and
internal heat transfer enter once. Inert species retain their original heavy
calorics. Exactly absent electrons remain absent and cannot receive power.

Each backward-Euler substep refreshes actual neutral and Coulomb conductances
at its start, then solves a scalar reaction residual with an inner monotone
thermal solve. The analytic constrained derivative includes both changing
electron population and heavy-pool thermal response. Safeguarded brackets
admit supported interior states without clipping temperatures or abundances.
Full/two-half-step estimates control this first-order coefficient lag and the
combined reaction/thermal evolution. Both temperatures, all reaction amounts,
and separate neutral/Coulomb transfers enter the error norm. Source admission,
physical mixture entropy and independent pool-energy balances gate publication.
Nonlinear trials allocate no state or species workspace; source snapshots are
refreshed between substeps. Cancellation leaves caller-owned input unchanged.

G1/G3/G4 regressions compare three refinement levels with an independent
NIST-BEB/Landau/analytic-neutral ODE calculation, exercise signed exchange and
held powers under formation-energy gauge changes, create an initially absent
neutral through recombination, and refuse stale source bindings or cancellation.
The native EM integration carries cell-volume heating through this same joint
update and refreshes transport from actual changed populations and temperatures.
These are numerical checks of the declared classical Maxwellian ground-state
model. They do not establish experimental discharge validation, radiation or
photoionization, excited-state kinetics, a continuously fixed-voltage transient,
spatial charged transport, or a general multireaction plasma network.

## Finite-plenum thermal equilibrium in gravity networks

`GasNetwork1d::step_with_thermal` couples the original frozen-gas transport to
a prepared `PlenumThermalOperator`. Each compact volume has one explicit
`PlenumThermalContact`: adiabatic, a prescribed-temperature reservoir,
authored finite capacity, or retained original `MaterialWallCaloric` storage.
The latter keeps its nonlinear h/u integral, original identity/geometry and
latent inventory; a plateau is not replaced by apparent heat capacity.
Distinct contacts represent distinct physical wall masses.

The existing implicit gas/wall source now borrows either the original frozen
or species-carrying gas calorics. Both adapters retain their actual source,
energy gauge and inversion; no equivalent species bank or constant-gamma
substitute is constructed. Existing mixture entry points and source-dependent
wall-conductance behavior remain on the same equation owner. Exact
source-representable thermal fixed points avoid introducing inversion
tolerance into zero heat flux. Passive source brackets, original source
domains, actual represented heat, energy residuals and bounded work remain
authoritative.

The wave owner first prepares its ordinary SSPRK2 candidate. Backward-Euler
plenum heat then uses that trial's actual accepted duration, including CFL and
source retries. This is first-order thermal/transport splitting. All gas
states, finite-wall states, clock, port exchanges and reports publish together
only after the final poll. Failed later contacts and final cancellation retain
both accepted owners; prepared source buffers are reused. Amount and volume
are fixed during heat exchange, so it adds neither momentum nor gravitational
work.

Gas balances subtract heat exactly once. `NetworkThermalReport` separately
records the original wall source and closes gas plus finite-wall storage
against external gas/reservoir transfer; with gravity it also includes
centered potential storage and external potential transport. The complete
absolute-potential residual is reconstructed from the complete centered
residual plus `phi_ref*mass_residual`. Formation energy retains the original
separate convention.

Seven focused native regressions cover exact stiff hydrostatic/thermal
equilibrium, independent variable-Cv finite-wall roots and long-time
equilibration, heating/cooling without temperature crossing, actual accepted
duration, nonlinear original wall storage, latent inventory, pressure-driven
gravity flow with independently reconstructed energy, large potential/formation
datum invariance, final cancellation/reuse and later-source rollback.
The seven new cases, all 87 dynamics tests and 14 existing shared-wall
regressions passed against the original production crate. The material-enabled
production path compiled. Its supplemental gate initially rejected three shared
fixtures with undeclared T/phi quantity types; their validity declarations now
match the original source support and strict query types. The corrected
material-only cases remain unverified: the rerun was killed by the shared-memory
limit during test compilation.

This remains a compact-volume, frozen ideal-mixture engineering model with
explicitly authored conductance. It supplies neither a validated heat-transfer
correlation, resolved wall conduction or shared wall mass, moving-volume work,
gas phase change, nor a second-order or entropy-stability certificate.

## Distributed wall equilibrium in gravity-aware gas networks

`NetworkWallOperator` extends the same prepared wall owner to explicitly
selected finite-volume branch cells and compact plenums. The original
`PlenumThermalOperator` name is an alias of this owner; existing plenum calls
retain their behavior. `BranchWallContact` binds original linear/Darcy drag
and reservoir, finite sensible, nonlinear or latent thermal storage. Original
source, complete mesh geometry, cell volume and distinct wall amount remain
explicit. A shared physical wall requires its own thermal owner rather than
duplicated stored energy.

`GasNetwork1d::step_with_walls` uses the existing transactional transport
step. The shared caloric wall kernel acts on unpublished branch and plenum
candidates over the actual accepted retry duration. Transport/thermal coupling
is first-order splitting. Drag converts kinetic to internal energy while
retaining total gas energy; only accepted wall heat changes that total.
Fixed gas amount and cell volume mean no additional gravitational work during
the wall source. Original variable heat capacity and latent h/u inventory are
retained without constructing an equivalent gas or apparent wall heat capacity.

Each accepted branch `StepReport` exposes signed axial wall impulse and
wall heat and includes both in its residual. Impulses along differently
oriented branches are not summed into a fictitious scalar network momentum.
The shared source ledger closes gas plus finite-wall storage and, where
declared, centered gravitational potential. Absolute-potential residuals
retain the existing reference-times-mass-residual conversion. Numerical
candidate arrays are allocated at preparation and reused. Gas, walls, clocks,
port exchanges and reports publish together after the final cancellation
poll; a later source refusal retains all accepted state.

Seven native regressions cover exact hydrostatic/thermal rest, analytic linear
and Darcy drag, independent variable-Cv finite-wall exchange, heat-driven
gravity flow with potential/formation-reference invariance, retained latent
inventory, final cancellation and later-source failure. The complete dynamics
suite passed 94 tests, including those seven; all 14 original shared-wall
regressions also passed (108 distinct tests). These used the original
production crate with `experimental-euler`, without the optional
`material` feature.

This remains a frozen ideal-mixture passage/plenum model with explicit wall
laws. It adds no gas phase change, resolved wall conduction, shared thermal
wall network, moving-volume work, calibrated heat-transfer correlation,
second-order splitting or entropy-stability certificate.

## Shared finite walls across frozen gas networks

With `experimental-euler`, `dynamics::network::NetworkWallOperator::with_shared_walls`
binds one finite original solid inventory to several branch or plenum gas
control volumes. A `SharedWall` owns that storage once; its `SharedWallPort`
values identify the original gas cells and their complete nonnegative
conductances in W/K. Sharing heat never mixes their mass or species. Local
contacts must be adiabatic without competing finite storage, and a gas target
may belong to only one shared group. Existing local linear/Darcy drag remains
active and converts kinetic energy to internal energy before thermal exchange.
The shared storage's own conductance is zero: every actual film belongs to its
declared interface.

Each group solves the simultaneous backward-Euler equations
`Q_i = dt G_i (T_wall,new - T_gas,i,new)` and
`Delta E_wall + sum(Q_i) = 0`. The unknown is an energy increment about the
initial wall state, not a difference of large absolute energy totals. Every
gas response reuses the original frozen-source caloric adapter and existing
safeguarded reservoir solver. Original source capacity integrals or T(h/u)
supply the wall storage and temperature, including latent plateaus; no fitted
gas Cv, apparent latent capacity, or repeated copy of the physical wall mass
is introduced. When original tangents exist, the scalar derivative uses the
individual gas heat responses and wall tangent; otherwise the same safeguarded
owner bisects. All interfaces retain the same solved endpoint wall state.
Reprojecting that state after convergence would change stiff interface fluxes
and is not performed.

Acceptance checks the actual represented wall and gas storage changes and
each original endpoint heat equation against the joule allowances. Large gas
formation-energy gauges cancel through the existing sensible-energy owner.
An original wall coordinate that cannot represent the required exchange must
meet the same actual-storage residual; a small solver update alone is never
acceptance. Original wall enthalpy retains its source's pressure-work convention
and is not relabelled as rigid-vessel internal energy.

The source uses the actual duration accepted by the transport retry.
Per-branch signed heat is added to the existing branch residual, while the
single wall inventory enters the global gas/wall and gravity-aware accounts
once. Fixed gas mass and volume produce no extra gravitational work during
heating. All branch/plenum fields, wall inventories, clocks and reports publish
together after the final cancellation poll. Source-domain refusal, insufficient
nested work, or cancellation preserves the preceding accepted history.
Prepared gas-response arrays and compensated reductions are reused.

`SharedWallControl` bounds the total original gas primitive/caloric and wall
state queries across all groups, brackets and final checks in one source
attempt; each original inversion retains its own bounds, and complete network
retries retain their existing cap. `SharedWallReport` retains the accepted
coupled storage residual, maximum individual wall and interface residuals,
actual source-query count, outer iterations and interface-solve count.

This is a fixed-geometry, spatially uniform common wall with explicitly
authored conductances. Transport/source splitting remains first order.
Initial wall temperatures must lie in the original gas caloric interval.
Spatial wall gradients, moving geometry, multiple shared groups on one gas
cell, inferred material conductivities and empirical heat-transfer qualification
remain separate capabilities. The scalar solve does not establish a general
discrete entropy theorem.

Seven focused physical regressions cover an independent three-body
backward-Euler solution with drag heating, nonlinear-Cv equilibration and stiff
passive response, large independent gas energy datums, exact hydrostatic and
thermal rest at a large common potential datum, one original latent inventory,
heat-driven field motion, source/work/overlap refusal, and final cancellation
with deterministic history replay.

### Shared walls with evolving mixture composition

With `experimental-euler`, the multicomponent network's
`wall_exchange::NetworkWallOperator::with_shared_walls` binds the same original
finite inventories and `SharedWallPort` values to evolving mixture cells and
plenums. The frozen and mixture adapters share one scalar energy-balance
kernel; every mixture query borrows the current constituent densities and
original gas caloric owner. Chemistry or transport between source halves
therefore changes the actual gas heat capacity and inversion used by the next
half, without freezing composition or replacing the source energy reference.

Local drag precedes shared heat. Every gas target must have an adiabatic local
contact without competing storage and may occur only once across the shared
groups. Any bound gas transport law also admits the heated endpoint. A moving
plenum uses the composing owner's current trial volume; only that owner can
change the accepted geometry, and its pressure work remains separate.

Candidate wall storage carries through both source halves and publishes with
the complete network transaction. Accepted accessors expose only the last
complete state. `SharedWallControl` counts original source queries across both
halves and resets on retry; `last_shared_wall_report` combines their storage
residuals, contact residuals and work. A failed second half cannot publish a
first-half wall state. Original h/u, latent energy, pressure-work conventions,
uniform-wall applicability, authored conductances and first-order composition
retain the limits described above.

## Original material conductivity on shared gas/wall ports

With `experimental-euler,material`, both the frozen and evolving-mixture
network wall owners expose `with_shared_wall_conduction`. Rows follow the
existing shared-wall order and entries follow each wall's gas-port order.
Each optional `MaterialWallConduction` retains its original temperature or
phase-dependent conductivity source and planar length/area geometry. A bound
port's authored conductance remains the gas-side film; solid resistance is
added once in series. Unbound ports retain their complete authored coefficient.
The resulting conductance is never written back as the next step's gas film,
and no material density, heat capacity or latent energy becomes extra storage.

`InitialState` samples original k at the beginning of each source half and
checks its final state. `ImplicitEndpoint` resolves k(T(h),phi(h)) at each
proposed wall inventory inside the existing simultaneous gas/wall balance.
The phase form requires the same original equilibrium chart as the wall's
caloric owner, so conductivity can change on an isothermal latent plateau.
The existing temperature-coordinate form retains its original source and
admitted geometry. Implicit evaluation requires original caloric storage.

The shared energy-root kernel is unchanged in ownership and physical balance.
An implicit conductivity binding suppresses the incomplete capacity-only
Newton slope; the bounded scalar solver uses bisection without inventing dk/dh
or asserting nonlinear uniqueness. Every actual conductivity query spends the
existing source-attempt budget, including failed proposals and final checks.
If wall-to-wall links subsequently move the inventory, all bound conductivity
sources admit that final candidate before publication. Domain, work or
cancellation refusal preserves the previously accepted gas and wall states.

This is a uniform-wall resistance model, not a resolved through-wall thermal
profile. The material supplies conductivity and the caller supplies the actual
resistance geometry and gas film. Source support at sampled states does not
establish a continuous-path certificate or empirical heat-transfer accuracy.

## Connected finite wall inventories

`NetworkWallOperator::with_shared_wall_links` connects the existing shared-wall
inventories through ordered `SharedWallLink` thermal ports. Each has two
distinct original wall indices and one nonnegative authored conductance in
W/K, complete unless material resistance is explicitly bound below. Cycles and
parallel interfaces retain authored order; a zero-conductance
link transports no heat. No gas mass or species passes through these ports.

After the original simultaneous gas/shared-wall solves, each link solves
`Q = dt G (Ta(-Q) - Tb(Q))` with the existing safeguarded scalar owner. Original
capacity integrals or T(h/u) charts supply both endpoint states, including
latent plateaus and their conserved phase energy. Exact available tangents
accelerate the same equation; absent local tangents use its bounded bisection.
For fixed conductance the heat interval intersects both original source domains
and the passive initial-flux bound. No source extrapolation, fitted capacity, duplicate wall
mass or alternate nonlinear solver is introduced.

Acceptance checks both actual represented storage changes and the stage heat
equation. Their signed sum enters the existing whole-network wall-energy
account; wall-to-wall transfer is internal and never counted as external heat
or gas heating. Source queries share `SharedWallControl` across gas groups and
all link stages, including failed evaluations and final checks. Prepared
receipt arrays are reused. Failure or final cancellation leaves gas fields,
wall inventories, clocks and all previously accepted reports unchanged.

`shared_wall_link_reports` retains every link's initial and endpoint original
states and temperatures, signed heat received by the second wall, both actual
storage changes, energy/flux residuals and scalar iteration count. These are
stage endpoints: a later link can change a wall again. The final wall state is
not substituted into an earlier gas or wall heat equation. All links use the
actual interval accepted by the original network transport/retry owner.

This is first-order composition of finite thermal ports, with observable
link-order error that must be controlled by timestep refinement. It is not a
simultaneous backward-Euler solve of the whole thermal graph, nor a spatial
PDE discretization. Geometric/material spatial conduction remains owned by
`fs-thermal`; no conductivity or geometric resolution is inferred from the
authored W/K ports. Existing fixed-geometry and lumped-cell validity conditions
remain. The scalar solves alone do not establish a general entropy theorem.

G0/G1 regressions cover independent staged gas/wall calorimetry; an analytical
three-wall chain with first-order convergence and vanishing link-order bias;
latent inventory transfer and exit into sensible cooling; exact thermal rest
with cyclic/parallel links and isolated storage; exact whole-attempt source
budget boundaries and late source-domain refusal; and atomic cancellation with
bitwise checkpoint replay. Source-specific calorimetry, energy closure and
temperature-hull checks do not claim experimental heat-transfer accuracy.

### Connected walls in evolving-mixture transactions

The multicomponent network wall owner also exposes `with_shared_wall_links`,
`shared_wall_links` and `shared_wall_link_reports`, using the same port types
and scalar exchange kernel as the frozen network. Each source half first
solves its current gas/shared-wall groups, then the authored wall links. The
next half starts from those candidate wall states and the newly transported
or reacted gas composition. Wall-to-wall heat changes only the two original
wall stores; it adds neither external heat nor another gas heating term.

The complete source-attempt budget covers both gas and wall links across both
halves. Report buffers are allocated at binding for two full link sequences.
The accepted report slice contains one chunk in authored link order for each
source half: one chunk for standalone `advance_for`, two for a composed
network step. Every receipt retains its actual half duration and original
stage endpoints. A later stage cannot replace the temperatures of an earlier
flux equation. Accepted states and the full report sequence publish together;
failure or cancellation preserves the previous complete transaction, and a
retry resets every candidate state, receipt and source-work counter.

Focused tests check independent three-wall staged calorimetry, transfer out of
an original latent plateau into sensible cooling, both-half receipts and exact
source-work limits, final cancellation with replay, and source-domain refusal
after gas heating. The same first-order, uniform-wall, authored-conductance and
source pressure-work limits apply. This extension adds no spatial gradient,
empirical heat-transfer qualification or simultaneous graph-endpoint claim.

### Original material resistance between finite walls

Both network wall facades expose `with_shared_wall_link_conduction` and
`shared_wall_link_conduction`. Rows follow authored link order; each row's
two optional `MaterialWallConduction` bindings belong to that link's first
and second wall. Each side retains its own original caloric chart, phase
state, conductivity source and length/area. Distinct materials never share an
averaged temperature, enthalpy datum or inferred conductivity law.

On a bound link, the original authored G represents finite interface/contact
conductance. Each supplied bulk resistance L/(k A) is added once in series
with it using scaled conductance arithmetic. The authored G remains unchanged;
zero still means insulation, and two absent layers retain the complete-G path.
The finite interface model does not infer perfect contact. Geometry describes
each bulk segment between its retained cell state and the common interface;
it creates no additional mass or latent-energy account.

Initial-state layers sample at that link's actual entering state, including
changes made by earlier links. Implicit layers instead query each proposed
endpoint; the two sides may use different evaluation modes. When conductivity
is implicit, the heat bracket uses original passive storage limits rather than
the initial-flux cap, because increasing k can make the accepted heat exceed
that initial estimate. Safeguarded bisection avoids inventing dk/dh or claiming
nonlinear uniqueness. Fixed and lagged links retain their available tangents.

Every runtime conductivity query shares the complete gas/wall source budget.
Lagged endpoints are admitted, and all bound link sources are checked again
after the full ordered sequence: a later link cannot leave an earlier source's
wall outside its support. Both facades retain the same original storage,
physical heat receipts and atomic publication. The first-order composition,
fixed-geometry, finite-cell and material-domain limits above still apply.

Focused G1/G3 cases compare independent two-material calorimetry in both
orientations and an increasing-conductivity latent root beyond the initial-flux
bound. G4 cases check source-work limits, cancellation/replay and late-link
material-domain refusal in both facades. Actual-source native compilation and
28 focused gas tests pass; workspace/DSR verification remains separate.

## Finite-time thermal exchange with multicomponent phase equilibrium

With `material,experimental-euler`, `material::phase_thermal::exchange_mixture_parcel`
couples the original `IsochoricMixtureEquilibriumParcel` to a held-temperature
reservoir, an authored constant-capacity finite wall, or one retained
`MaterialWallCaloric`. The parcel's immutable feed, mass and vessel volume
remain unchanged; pressure, temperature and phase amounts follow the original
source-bound U,V equilibrium. Every trial uses that same caloric owner and its
finite phase-screen evidence. No equilibrium Cp/Cv, dilute transport, extra
latent heat, chemical reaction or phase-change kinetics is inferred.

The existing bounded wall scalar owner solves
`Q = dt G (Twall(-Q) - Tgas(Q))` in signed heat coordinates. `G` is the complete
authored conductance in W/K, fixed during the supplied physical interval.
The caller declares a heat bracket; it is intersected with the initial passive
flux bound and original finite-wall storage limits. Every remaining gas endpoint
must satisfy its original pressure/temperature/source brackets. An unsupported
trial, missing sign bracket or unresolved original phase search refuses rather
than extrapolating or switching to a homogeneous approximation. Numerical
secants accelerate the safeguarded scalar root; they carry no heat-capacity or
physical-tangent authority. The accepted source-bearing trial is retained and
published without another U,V solve.

Actual represented gas energy, original finite-wall h/u changes and endpoint
heat flux must independently close in joules. Original wall enthalpy keeps its
declared pressure-work convention; it is not relabelled rigid-volume internal
energy. Reservoir heat is external, finite-wall exchange is internal, and
latent storage already belongs to the original thermodynamic potentials.
The U,V residual allowance reserves part of the thermal budget; a nonzero
initial flux tightens the absolute root allowance so a small exchange does not
vanish into a loose configured tolerance. Source increments cancel species
formation and entropy datums before rounding.

`MixtureThermalControl` bounds heat queries and all original PT selections for
the complete interval, in addition to the existing per-trial pressure/phase
limits. Remaining PT work is passed into every original U,V solve. Reports
retain actual aggregate heat/PT/pressure work and the sum of original EOS work
ceilings, separately from the accepted trial's caloric report. Final
cancellation, source/balance failures and work exhaustion leave both supplied
inventories unchanged. Zero duration, zero conductance and equal initial
temperatures retain both states exactly with zero phase evaluations. A cloned
parcel and retained original wall state supply deterministic retry state.

This is a first-order backward-Euler, spatially uniform thermal interface.
Uniform gas/wall assumptions, authored conductance applicability and time
refinement remain caller obligations. It does not establish spatial conduction,
evaporation rate, slip, moving-boundary work, global Gibbs stability, calibrated
heat-transfer correlations or experimental fluid accuracy. The gas source
entropy increment is retained; an arbitrary finite wall may lack its own
entropy law, so the API makes no combined entropy-production theorem.

The public `phase_equilibrium_thermal` integration target adds six focused
G0/G1/G3/G4 tests: independent binary isochoric calorics with a held bath;
finite-wall heating across the phase boundary; original latent-wall cooling;
independent species-datum and extensive-scaling invariance; exact no-work rest
and finite source-domain refusal; and whole-call work, final cancellation and
deterministic retry. The methane/ethane fixture uses synthetic NASA capacities;
its original cubic/fugacity references test constitutive consistency, not a
measured-mixture validation claim.

Native validation compiled the production path and passed five of these six
cases. The finite-wall phase-crossing fixture now declares a supported heat
bracket and a 256-iteration original flash budget, retaining the original
strict tolerances. Its final rerun was blocked by shared disk exhaustion and
remains unverified; neither that case nor the complete suite is claimed green.

## Adiabatic equilibrium mixture compressors and expanders

With `material`, `material::phase_machinery::adiabatic_mixture_machine`
resolves a steady adiabatic compressor or expander over one immutable
`ThermochemicalPengRobinsonModel`. It first solves the original pressure/entropy
increment constraint `Delta s = 0` at the outlet pressure. For an authored
isentropic efficiency `0 < eta <= 1`, the actual enthalpy target is
`Delta h_s / eta` for compression and `eta * Delta h_s` for expansion. The
existing pressure/enthalpy increment owner resolves the actual outlet using
exactly the same original phase-screen policy. Source intervals may differ
between the ideal and actual temperature solves; phase policies may not.

Both solves retain the complete feed, original NASA laws, EOS departures,
phase mixing entropy and finite original phase-screen evidence. Liquid/vapor
amounts may change. Species energy and entropy reference datums cancel through
the original integrated increments before rounding. No mixed heat capacity,
extra latent energy, dilute transport extrapolation or replacement EOS enters
the machine model.

The signed shaft work **into** the fluid is the actual source-integrated
enthalpy rise. This steady-flow balance assumes no external heat and negligible
inlet/outlet kinetic and potential energy differences. Stream enthalpy already
contains pressure-flow work; another pressure-volume debit would double count
it. The result retains the ideal endpoint, actual endpoint, source work,
actual entropy production, enthalpy residual and ideal entropy residual.
The caller may convert work to power with an independently supplied nonnegative
typed mass flow. Pressure ratio and efficiency do not predict that flow.

Admission checks the pressure direction and representable work sign. Actual
source entropy change must be nonnegative within a declared numerical
allowance; the reported value is never clipped. A single bounded PT allowance
covers both original solves, including repeated inlet queries and all bracket
states. Each nested owner receives only the remaining allowance. Cancellation,
source-domain refusal, unresolved phase work or an exhausted budget returns no
partial outlet. Efficiency one reuses the accepted PS endpoint; equal pressure
returns exact zero work and entropy change after one original inlet selection
when the reference is inside its admitted bracket.

This is an equilibrium endpoint model with a caller-authored efficiency, not a
validated machine map, rotor/transient model, reaction/phase-kinetic model,
resolved flow geometry, choking model, or global phase-stability certificate.
The declared tolerances bound numerical residuals, not model-form error.

Seven focused regressions cover an independent two-phase isentropic expansion,
reversible round trip and both efficiency conventions, actual nonnegative
entropy production, large independent species reference shifts, exact rest,
whole-call budget boundaries, malformed policies and final cancellation/replay.
Rustfmt parsed the new source and tests. Native compilation and execution remain
unverified because the shared build disk is exhausted; no passing native result
is claimed for this new module.

## Simultaneous reaction and finite reservoir heat exchange

The prepared constant-volume and constant-pressure reactors accept an explicit
`HeatBath`: nonnegative total conductance G [W/K] and a held positive bath
temperature Tb. The original extent integrator advances Qdot=G(Tb-T) together
with the reaction laws. Both backward Euler and SDIRK2 eliminate stage heat
inside the original caloric root; constrained temperature derivatives use
original reaction delta-u/(Cv+hG) or delta-h/(Cp+hG), with reference-volume
scaling. Physical heat-capacity reports remain Cv/Cp. There is no additional
integrator, split temperature update, reaction heat source or heap workspace.

Original initial U/H targets retain their meanings. Accepted energy balances
subtract compensated external heat once; fixed-pressure results retain actual
volume and p delta-V work. Full, half and rejected stages carry independent
heat state. Step doubling checks heat as well as species and temperature using
a formation-datum-independent scale. Caloric residuals and local heat
differences retain the compensation terms.

Parcel entropy may decrease on cooling. Admission uses parcel-plus-reservoir
production delta-S-Q/Tb; reports retain parcel delta-S, reservoir transfer,
total production and minimum accepted-step production. Original element,
charge, mass, nonnegative inventory, source support and energy gates remain.
Zero conductance preserves adiabatic results. Explicit-energy zero-duration
calls resolve initial U/H before any exchange.

Six new native tests passed: frozen heating/cooling, simultaneous reaction and
heat against an independent closed form including resonance, fixed-pressure
first law/work, large elemental formation-reference shifts, insulated parity,
cancellation/replay, explicit targets and constrained-Jacobian differences.
The bounded native fs-gas library run passed 68 tests, with two lengthy existing
hydrogen/methane reference cases excluded; the SDIRK hydrogen reference passed.
These are focused G0/G1/G3/G4 checks, not a full repository Gauntlet result.

Scope remains a sealed homogeneous ideal-gas parcel with caller-declared
constant conductance and an infinite fixed-temperature reservoir. No finite
wall capacity, evolving wall source, radiation law, streams, spatial transport
or empirical mechanism validation is implied. Time-error estimates are local.

## Finite material bodies in adaptive reacting parcels

With the `material` feature, `chemistry::heat::material::FiniteHeatBody`
retains an exact original equilibrium enthalpy chart, invariant body mass,
declared constant body pressure and total thermal conductance. Both prepared
gas reactors expose `advance_with_finite_body`; the returned gas and original-
source body endpoints publish together. Source, cancellation, work and energy
refusals publish neither endpoint.

The existing backward-Euler/SDIRK2 chemical stages eliminate the finite body
on its original piecewise-linear caloric segments. Binary segment lookup and
analytic elimination introduce no second time integrator or nested iterative
temperature root. The constrained gas Jacobian includes the body's changing
temperature response. Latent plateaus retain their actual enthalpy and phase
fraction; the adaptive comparison includes both body temperature and liquid
fraction. Caloric bracketing respects both original source domains.

Gas heat and body heat cancel once. The original material isobaric heat owner
reports density-derived volume, pressure work, enthalpy and internal-energy
changes. Fixed-volume gas closes its U change with body U and body p*dV;
fixed-pressure gas also reports its own outward p*dV. Formation and latent
energy remain in their respective original potentials. Entropy admission uses
the actual gas change plus the original body integral of dh/T. Unresolved
enthalpy-datum heat increments refuse rather than disappear.

This is a homogeneous finite-body model with constant conductance, held
pressures and no gas/body mass exchange. The caller must justify negligible
internal body gradients and the pressure applicability of the original chart.
It does not resolve wall geometry, contact mechanics, moving interfaces or
spatial flame transport. Source traversal is bounded to 4096 original knots.

Five focused G1/G3/G4 tests cover independent sensible/reacting two-body
solutions, sharp latent heating/cooling and onset, both gas constraints and
time methods, enthalpy datum changes, both pressure-work accounts, finite-body
constrained Jacobian differences, source refusal, cancellation, replay and
insulated agreement with the original adiabatic owner. The final actual-source
material-feature heat target passed all 11 tests (five new, six existing).
A broader material-feature unit run passed 78 tests, excluding the two known
long hydrogen/methane Cantera cases, before the final finite-denominator guard
and additional Jacobian test; that final focused run includes both changes.
Default-feature compilation and all six existing heat tests also passed.
These are native fallback builds against compatible cached dependencies;
DSR and required sibling checkouts were unavailable. They do not establish
a full workspace verdict or measured-material validation.

### Finite support inertia in the original network wall transaction

`NetworkWallOperator::with_finite_supports` binds optional independent
`SlidingSupport` masses/velocities in the original branch/cell order. Selected
supports replace prescribed wall motion, retain authored or source-refreshed
resistance, and use the original exact reduced-mass recoil mechanics. All slip
loss thermalizes once in gas sensible energy. The existing local,
material-caloric, shared-thermal and conduction owner then advances heat with
drag disabled; no thermal inventory is copied or reconstructed. Both material
transport endpoint checks use actual support-relative slip, and the accepted
contact retains the returned support velocity.

Support candidates follow the same begin-trial, source-half, retry and publish
lifecycle as gas and walls. Workspaces allocate at attachment; stepping adds no
allocations or clock. Read-only support states and per-cell `SupportExchange`
receipts retain actual complete support endpoint impulse/kinetic change and
accumulated friction-only gas impulse/energy. Independent dimensional momentum
and energy budgets reject before publication.
`NetworkWallReport.mechanical_work_to_gas_j` includes finite-support work as
internal exchange with that support; the support receipt closes this additional
kinetic storage. Branch axes are not summed.

Each selected cell has its own mechanical support; identical descriptors do not
share inertia. Shared thermal storage remains separate. The feature does not
move geometry, evolve plenum bulk momentum, infer friction coefficients or
supply multi-cell rigid support dynamics. The all-to-gas slip heat partition is
explicit; other heat partitions require a jointly coupled thermal owner.
Five selected-source native tests cover independent analytic recoil/heat
balance, both halves and endpoint accounting, cancellation/retry,
material-relative domain checks/Galilean boosts, and unselected identity and
admission. Two additional full network/shared-thermal source tests were authored
and parsed but not executed in the delivery environment; no full Cargo,
material or advection pass is claimed.

### One finite support exchanging momentum with multiple gas amounts

`SharedRecoilWorkspace` advances fixed gas amounts against one common
`SlidingSupport` along an authored scalar axis. Each `SharedRecoilCell`
retains its own density, volume, momentum and nonnegative frozen linear and
quadratic resistance coefficients. For unequal drag laws or heterogeneous nonlinear states, the original exact
reduced-mass pair mechanics is composed in a palindromic Strang sweep. This general path
approximates simultaneous shared-support dynamics with one common inertia.
Single-active-cell intervals reuse the exact original pair map. Multiple
active cells with one common linear rate and no quadratic drag use the exact
collective/deviation modes described below.

Bounded step doubling estimates local errors in each velocity, physical slip
heat and Galilean-invariant work. The support velocity is included. Dimensional
velocity/energy tolerances and a relative tolerance are separate from the
strict endpoint momentum and energy budgets in `RecoilControl`. On the general
refined path, each candidate interval also resolves the resistance scale:
`h*(max(k_i) + sum(m_i/M*k_i)) <= 0.5`, with
`k_i=a_i+2*b_i*R` and R the initial velocity-hull width of the support and active cells. This prevents
saturated split pair maps from masquerading as resolved heat allocation.
Rejected attempts consume the explicit step budget; under-resolved,
non-advancing or exhausted intervals refuse. These are local time-error
estimates, not a certified global or stiffness-independent accuracy bound.

Every exact pair flow exchanges opposite impulse and dissipates nonnegative
physical slip energy. Accepted endpoints audit actual represented momenta,
support kinetic storage, and per-cell kinetic change/work/heat. Numerical
residuals are reported or refused; they are not deposited as heat. Work and
heat remain separate per-cell outputs. The caller's original caloric owner
must apply total-energy work and retain the associated slip heat once.
Thermal storage, friction calibration, geometry motion, arbitrary 3D support
topology and a new flow clock are outside this mechanical interface.

Workspace allocation occurs at preparation; advancing uses retained buffers
and bounded cancellation checkpoints. A refusal returns no accepted output,
and retry reinitializes all candidates from the supplied original inputs.

Seven focused G0/G1/G3/G4 native tests check unequal-mass/rate analytic
shared motion and individual heat/work, tolerance refinement, symmetric
nonlinear decay, initially comoving cells, Galilean changes and cell order,
clock subdivision, exact single-active parity and workspace reuse, stiff-rate
refusal, and cancellation at every exposed checkpoint with exact retry.
The run compiled original mechanics/drag/controls and fs-math/fs-qty; the
surrounding parent error type was a test fixture. No full crate or workspace
build is claimed.

### Shared branch inertia on the original network source transaction

`NetworkWallOperator::with_shared_branch_supports` binds one optional
`SlidingSupport` per original branch. Every cell on a selected branch reacts
against that one mass and velocity through `SharedRecoilWorkspace`.
Independent per-cell supports may coexist on other branches; overlap is
rejected in either binding order. `None` retains the existing path. Distinct
branch axes are never summed or treated as one spatial body.

Each original source half resolves all selected coefficients at its actual
entering common support velocity. Shared mechanics then supplies per-cell
momentum and mechanical work before the original thermal owner runs with
drag disabled. Slip loss thermalizes once through total energy minus kinetic
energy. Local wall capacity, original material calorics/conduction, and shared
thermal storage retain their existing ownership. Original transport admission
checks actual relative-slip endpoints, including after shared thermal exchange.
No thermal inventory, mass, mesh volume or flow clock is added.

Accepted support states and `SupportExchange` receipts publish in the same
transaction as gas and thermal walls. Both source halves accumulate actual
friction-only gas impulse/energy; the common support's impulse and kinetic
change use the complete accepted endpoints. Cancellation, exhausted mechanical
work, original source-domain failures and balance refusals preserve accepted
gas, support, thermal states and reports. Preparation allocates the work arrays;
stepping does not allocate. The fixed-axis/frozen-geometry and local estimated
mechanical accuracy limits of the shared solver remain explicit.

Six focused selected-source native tests passed with the final solver and
network source: nonuniform two-species volumes/masses with linear, Darcy and
combined friction plus original local heat; an independent two-rate shared
continuum with individual heat; both source halves and endpoint receipts;
prior-state cancellation/replay and work refusal; disjoint independent
supports, overlap rejection and None-path identity; and support-relative
transport/Galilean endpoint admission. The native run used the actual network
wall owner, heat, friction, transport and mesh source with explicit constant-Cv
gas/network/transport fixtures. Material and shared-thermal paths were trapped,
not exercised. Full Cargo/material/advection verification remains separate.

## Two reacting chambers sharing one finite piston

`ReactingEuler1d::step_with_shared_piston` composes both original chamber
mechanisms with the existing two-field, moving-volume pressure/piston solve.
Each source bank retains its species order, caloric references, original
reactor and chemistry controls. The second source halves use the actual
moved cell volumes. A shortened pressure interval restarts both complete
source/pressure/source paths from their accepted endpoints.

Whole-system admission uses both actual chemical-energy changes, the native
piston kinetic change and applied-force displacement work. The two fixed
outer support impulses remain explicit in momentum closure. Independent
source energy gauges are never mixed, and no extra heat-release term is
introduced. Existing trial fields/reactors are reused; both fields, meshes,
clocks and the common piston publish only after the final cancellation poll.

The native regression cases cover independent exothermic mechanisms and
source references, cumulative physical balances, shared-duration replay,
postcompression source refusal and cancellation/retry after prior motion.
They parsed with the pinned formatter but were not compiled or executed in
the constrained editing environment. This entry is adiabatic chemistry plus
Euler transport; molecular transport, finite-wall heat and leakage retain
their separate owners. Symmetric splitting alone establishes no accuracy
order, entropy-stability theorem or combustion-mechanism validity.

### Exact collective modes for common linear shared-support drag

When every active gas amount has the same positive linear acceleration
coefficient and zero quadratic coefficient, the shared recoil workspace uses
the exact collective/deviation solution. Active gas mean/support slip decays
at `a*(1+sum(m_i)/M)`; differences between gas velocities decay at `a`.
Arbitrary gas masses and unequal initial velocities are retained. The mode
uses O(cells) work and one accepted attempt, including stiff/light-support
cases that would otherwise exhaust the generic resolved-step budget.
Disconnected zero-coefficient cells keep their original state and do not
enter the collective mass.

Per-cell physical slip heat and support work use integrated exponential
moments with explicit cancellation-safe weak-coupling evaluation. In
particular, an initially comoving cell retains its induced second-order
impulse and third-order heat. Work remains separate from heat; no kinetic
roundoff residual is relabelled dissipation. Existing actual endpoint
momentum/energy acceptance, compensated accumulators, cancellation and
private-candidate publication remain authoritative. Other shared systems use
the original bounded refined solver, except for the uniform nonlinear groups
described below. The one-active-cell path retains the original exact pair
calculation.

The original network support owner receives this exact mode through its
existing workspace call. No network clock, source splitting, caloric storage
or material-coefficient ownership changes. Source coefficients remain frozen
over the original supplied interval; an exact frozen mechanical solution does
not establish accuracy of a changing material law or the complete multiphysics
composition.

Twenty-one focused actual-source native tests passed: four new collective
regressions, seven shared-support regressions, and ten original drag/pair
regressions. Independent high-precision references cover unequal masses and
initial velocities, individual heat/work, boosts and cell order, both sides
of the moment-evaluation switch, initially comoving cubic heat, sub-ULP mass
contrast, overflowed but fully decayed exponents, disconnected cells, and
cancellation at every exposed checkpoint with exact retry and workspace reuse.
The run compiled the actual controls, mechanics, shared solver and drag with
fs-math/fs-qty; unrelated surrounding error variants were explicit fixtures.
The added stiff common-rate network regression was source-reviewed and
formatter-parsed but not compiled or executed. Full Cargo/network/material
verification is not claimed by this selected-source run.

### Exact shared recoil for uniform nonlinear gas groups

If all active gas cells enter with exactly the same actual velocity and
exactly the same frozen linear/quadratic acceleration coefficients, a
positive quadratic law preserves their common velocity. The shared solver
now advances this invariant group through one original exact reduced-mass
pair map using the compensated sum of active gas masses and the one actual
support inertia. Pure quadratic and combined linear/quadratic laws both
use one accepted interval, independent of the collective relaxation rate.
The existing general common-linear mode still handles unequal velocities.

Each gas cell receives its own mass share of physical slip heat and
integrated work in the initial support frame. The original outer owner
restores work using each cell's actual represented impulse and checks actual
gas/support endpoint momentum and energy. No conservation residual becomes
heat. Disconnected cells do not enter the group and remain unchanged.
Group detection uses exact equality, with no tolerance merging different
states or laws; heterogeneous systems retain the original bounded solver.

The implementation uses the existing workspace with no per-advance vector
allocation. The aggregate unit-density coordinate is private mechanics,
not a new gas cell, volume, caloric state, transport clock or material model.
Original network source owners receive the mode through their unchanged
shared-workspace call. Coefficients remain frozen over the supplied source
interval; material evolution and complete source/flow splitting accuracy
retain their original ownership.

Four new actual-source native regression groups check independent Riccati
slip decay and individual heat/work for both slip signs and pure/combined
laws; unequal density/volume and mass partition; boosts, ordering and clock
partition; stiff/light-support one-interval success and exact heterogeneity
fallback; disconnected cells; and every cancellation checkpoint with exact
retry and retained workspace storage. All 25 selected-source native tests passed, including
the existing common-linear, shared, drag and pair regressions. Full Cargo
and network/material integration execution remain outside this native run.
