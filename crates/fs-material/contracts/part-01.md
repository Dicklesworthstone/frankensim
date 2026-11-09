<!-- fs-material contract, part 1 of 4; a higher number is newer. Index: ../CONTRACT.md. Keep this file under 900 KB: the GitHub contents API returns empty content above 1 MiB, and a writer that sees an empty file overwrites it. -->
## Physical viscous stress and shear traction

`fluid_transport::viscous` evaluates the three-dimensional isotropic law
`tau = 2 mu dev(sym(grad u)) + zeta div(u) I` in physical SI units and native
tensor-Voigt xy/yz/zx order. `NewtonianViscosity::from_transport` consumes the
original gas/liquid snapshot; independent nonnegative bulk viscosity is always
explicit. Spin produces no stress and isotropic expansion produces no shear.
Thermodynamic pressure stays with its original owner.

The same response supplies oriented normal/tangential traction, nonnegative
physical `tau:D` heating, and an allocation-free exact tangent action. The
optional 6-by-6 tangent uses tensor double-contraction weights. Tangential
traction is evaluated from deviatoric stress before adding bulk stress.
`ThermalCarreauYasuda::viscous_response` evaluates the original bounded
temperature/shear law at `sqrt(2 dev(D):dev(D))`, retains its original response,
and returns full strain-rate and temperature derivatives. The positive original
differential stress avoids subtractive cancellation in the thinning tangent;
the tensor tangent remains finite at rest even when the scalar viscosity slope
is singular. Invalid states or unrepresentable outputs return no response.

This is a local constitutive evaluation, not a velocity-gradient reconstruction,
momentum advance or heat insertion. Shear calibration is not an extensional or
compressible-flow validation; independent bulk viscosity has no inferred
temperature law. No elastic memory, yield stress, turbulence or numerical
loss is added. G0/G1/G3 regressions cover analytic Couette/extension/dilation,
all shear orientations, rotation/spin invariance, tensor and thermal derivative
finite differences, work conjugacy, bulk/shear separation and source support.

## Independent heat and work in the driven-mixture first law

`DrivenMixtureEquilibriumParcel::advance` checks the independently supplied
heat, outward boundary work and represented internal-energy change with the
existing error-free `TwoSum` primitives. Its reported first-law residual is
`Delta U + W_by - Q`; it retains the low part of `Q - W_by` instead of comparing
only against that rounded target. If target rounding alone exceeds the actual
extensive energy allowance, admission refuses before any phase query. The
accepted source increment receives the same compensated physical audit.

This prevents a small heat input from disappearing beside much larger
mechanical work. The original equilibrium solver, source bounds, tolerances,
mass/species inventory and transactional publication remain the owners. The
allocation-free correction changes accounting precision, not the constitutive
model or physical heat/work convention. Focused numerical regressions cover a
false-zero binary64 balance, exact closed/insulated accounts, nonfinite inputs,
and early refusal with the parcel unchanged.

### Causal rigid-frame porous diffusion

`porous::diffusive::PorousDiffusive` retains the original `PorousFluidLaw` and one coherent thermoviscous fluid snapshot. High-frequency storage is `alpha_inf*rho/phi` and `K_s/phi`, with the original static Darcy drag and positive viscous/thermal relaxation arms. Analytic positive measures are discretized by bounded midpoint quadrature; thermal relaxation includes its isolated real pole. Cp=Cv removes only thermal memory. Fluid substitution keeps intrinsic permeability and actual adiabatic modulus. Admission checks both complex response and each nonzero real loss against the original JCA/JCAL law on 129 logarithmic frequencies inside the card's support. Reported maxima are sampled discrepancies, not continuous-band bounds. This frozen rigid-frame realization supplies coefficients rather than an integrator; it does not qualify elastic-frame waves, mean flow, phase change, or source accuracy outside its declared band.

## Free, clamped and passive-support elastic surface reflection

`elastic::acoustic::interface::surface::IsotropicElasticSurface` retains one
original explicitly isotropic elastic state. `TractionFree` imposes zero
traction; `Clamped` imposes zero velocity at an ideal stationary support. The
solid occupies z < 0, with the bonded-interface owner's modal velocity, outward
traction and peak phasor conventions. There is no artificial second material,
finite-impedance approximation to a clamp, or guessed loss.

Preparation reuses the original tensor-based traces and critical-slowness
arithmetic. One row-scaled 2x2 native complex LU supplies coupled P/SV
reflection; SH has its independent free/clamped sign. The fixed-size prepared
operator applies coherent amplitudes without allocation or another solve.
Converted evanescent P fields retain phase and traction while carrying zero
individual normal power. Finite drives enforce the original half-space strain
bound and independently check reflected power and surface work. Incoming
nonpropagating modes, anisotropic sources, singular systems and unrepresentable
inputs refuse; preparation is cancellable.

G0/G1/G3/G4 regressions compare free/clamped analytical coefficients, normal
reflection phase and clamp reaction, coherent power unitarity, supercritical
conversion, common material scaling, original domain limits and cancellation.
This is an infinite bulk-wave surface model, not a flexural plate-edge law,
Rayleigh-wave source, nonlocal support, rough surface, bulk attenuation model or
experimental qualification of the source material.

`ElasticSurfaceCondition::Impedance` adds a passive local support with separate
normal and isotropic tangential specific impedances [Pa s/m], evaluated by the
caller at the original frequency. In the same exp(-i omega t) convention,
`traction = -Z velocity`; nonnegative resistance absorbs power and either
reactance sign retains phase. The same P/SV solve and an independent SH solve
preserve oblique mode conversion and outgoing evanescence. Bounded boundary
coordinates reconstruct weak wall motion beside stiff supports from both
coherent traces. Absorption is positive local `sum Re(Z_i)|v_i|²/2`, separate
from the independently evaluated surface flux and reflected-power deficit.
G0/G1/G3 cases cover matched normal termination, oblique SH coefficients,
coherent supercritical reflection, purely reactive supports, weak resistance,
stiff-support absorption below reflection roundoff and invalid active loads.
This extension supplies no causal time-domain realization, spatially nonlocal
support or measured support-law qualification.

`ElasticSurfaceCondition::TensorImpedance` admits a finite complex-symmetric
3x3 specific impedance in the same Cartesian surface frame. Its original real
resistance uses the shared passive-work Gram admission, including zero and
rank-deficient loss without a positive floor; resistance and reactance need
not share principal axes. `from_resistance_gram` instead retains three original
dissipative channels directly, preserving their null directions and weak loss
below Cartesian coefficient roundoff; it only forms their traction matrix.
Normal/tangential coupling converts all P/SV/SH
polarizations while retaining the same original isotropic bulk material. One
row-scaled native 6x6 solve prepares reflected amplitudes and actual wall
velocities together, preserving stiff-support motion without cancellation or
normal equations. Coherent drives reuse fixed-size maps and evaluate positive
loss as `sum |L^T v|²/2`, separately from the original wave-flux check. The
axisymmetric path keeps its smaller solve. Analytical rotated shear and
normal/shear conversion, noncommuting resistance/reactance, evanescence,
weak/rank-deficient loss, stiff motion and bounded cancellation cover the
tensor path. This remains a linear local reciprocal boundary law at one
frequency; it adds neither bulk anisotropy nor nonreciprocal support physics.

## Nonideal passive plant residence

`recovery::residence::nonideal::NonidealTankCascade` advances one named
conserved constituent through up to 64 unequal active/stagnant tank pairs.
Declared total flow splits into an instantaneous bypass and the train flow;
active residence uses the latter. Explicit counterflow exchanges mass with
stagnant regions, including during an external-flow stop. Zero exchange keeps
a dead region exactly isolated. Geometry, initial inventories, original
provenance and each interval's flows remain explicit.

Positive exponential propagation independently integrates the train outlet and
both exchange directions. Active, stagnant and overall balances must close;
there is no balancing plug or normalization. Work is bounded by 4096 transition
exposures per interval and 256 series terms per substep. Invalid, nonfinite,
unrepresentable or cancelled calculations leave the original state unchanged.
This is a passive forward transport model, not a recovery law, fitted residence
distribution, tracer deconvolution or two-phase slip model. Binary64 results are
not interval enclosures. Analytical regressions cover ideal-tank reduction,
bypass, unequal tanks, exchange-only shutdown, retained dead inventory, resumed
flow, partitioned intervals and transactional refusal/cancellation. The
`nonideal_residence` example and
`data/matdb/seed-v1/MINING_NONIDEAL_RESIDENCE.md` demonstrate a complete authored
pulse/stop/restart/flush inventory account without implying plant calibration.

## Experimental Austin grinding population balances

`comminution::austin::resolve_austin_grinding` resolves seven scalar Austin
parameters and an exact class-size curve from one original observation with
identical support. Typed envelope queries admit every source class and the
complete requested contact interval. The existing `SelectionBreakageModel`
retains original claims, observation, uncertainty and verified receipts through
prepared passes. Published B_ij indexing generates conservative daughter columns
without clipping or normalization; the terminal class collects unresolved fines
with zero selection. Unsupported operating conditions and mixed sources refuse.

The `quartz-assiros-petrakis-2017-austin` pack supplies three Table 2 laboratory
calibrations and 24 claims. Its exact six-node sieve grid and 0–240 s declared
calculation window preserve the source's operating conditions and its reported
mismatch between original fitted and measured product size distributions.
The incomplete Table 3 amplitude is not reconstructed. This fixed-grid batch
model does not infer plant scale-up, liberation, power or media wear. The
`austin_grinding_reference` example and `MINING_AUSTIN_GRINDING.md` demonstrate
source-to-population-balance use; analytical and compiler-backed regressions
check source selection rates, first/second populations, daughter mass closure,
element transport, retained source identity and domain refusal.

## Bounded grind-to-spec inversion and declared-power batch energy

`comminution::kinetics::target` adds `SelectionBreakageModel::time_to_target`
from an actual feed PSD to a cumulative fraction in a named original source
class and every no-larger class. It verifies that active target classes cannot
transfer mass outside the target set and calls the original `prepare` and
`transform` owners for every trial. Successful nonzero results retain actual
below-target and meeting-target points within the declared time tolerance,
the original product and source evidence, and the full forward-call count.
Comparisons include twice the original Poisson-tail bound and an explicit
caller rounding allowance; ambiguous crossings, unreachable targets, exhausted
budgets, cancellation and original source refusals return no accepted time.
These nominal binary64 allowances do not cover experimental or model uncertainty.

`GrindTargetSolution::energy_at_constant_net_power` accounts for `E = P_net*t`
at both time endpoints and `E/m_dry` at the meeting endpoint. Net power and
batch mass are explicit inputs; an admitted `dry-charge-mass` must match.
The PSD flow scale supplies weights only. No power law, idle loss, arbitrary
sieve interpolation or plant scale-up is inferred. Analytical mixed-feed and
secondary-breakage tests check crossing brackets and conservation; the genuine
source compiler integration independently checks all three original Austin
calibrations against `ln(2)/S_0`. The reference example accepts optional
`--net-power-w DECLARED_WATTS`; see `MINING_AUSTIN_GRINDING.md` for scope.

### Partial-enthalpy species exchange between real-gas volumes

`ThermochemicalPengRobinsonModel::partial_calorics_at_density` combines each
original NASA species enthalpy with its original mixture-dependent fixed-T,p PR
partial departure and converts using that species' molar mass. The borrowed
result retains the exact material and EOS source identities, canonical species
order and fixed-capacity storage. Its query allocates nothing. The original
homogeneous mechanical, NASA-support and complete-Cv admission checks apply;
actual endpoint resolution separately checks all original transport sources.

`IsochoricRealGasParcel::exchange_species_with` advances simultaneous signed
species streams from one immutable pair snapshot. Positive rows move left to
right, negative rows move right to left. Each stream carries its actual upstream
partial enthalpy, including flow work once. Opposing species may exchange energy
and composition with exactly zero net mass. Final volumes and extra source axes
remain fixed; unchanged total mass retains the original accepted density bits.
The returned candidate publishes both reconstructed endpoints or neither.

The consumer reuses the original ideal-source catalog/inventory join and the
original real-gas critical/pair binding owner. Shared NASA, transport, pack,
critical and known interaction data must agree exactly. Every new PR species
pair must be supplied by an original source; absent cross-catalog data is never
filled with zero. Rows are canonicalized by SpeciesId, duplicate/unknown rows
refuse, and each withdrawal is bounded against its original species inventory.
Intended and actual local/global species masses retain absolute roundoff budgets.

For each side, its final composition and density at its old temperature own the
caloric anchor. Before arithmetic, ideal formation datums cancel by conserved
species inventory. The remaining extensive target is
`m0*uR0 - m1*uRanchor + sum(delta_n_i *
[Delta_h_ideal_i(T0,Tup) + hbarR_i(up) + R*T0])`.
Molar counterflows are summed before multiplying by RT. The existing finite PR
caloric owner then reconstructs temperature and the complete pressure, acoustic
and transport state on that final-composition homogeneous branch. Local and
paired energy residuals use these small original-law increments; they never
subtract enormous displayed formation-energy totals. The unchanged temperature
is admitted for a sub-ULP caloric target only when it is the closest represented
original-law endpoint within the unchanged caller tolerance, with the complete
unmodified residual still reported and checked.

This is a first-order prescribed stream primitive. It supplies no rate law,
membrane/separation work, diffusion coefficients, phase flash or entropy claim
for arbitrary authored fluxes. Those remain responsibilities of the physical
transport owner. Source-domain, representation, work-budget, conservation and
cancellation failures publish neither endpoint.

Focused regressions check dense zero-net-mass counterflow against original
absolute EOS energy, Euler partial-enthalpy closure, recovery of the existing
bulk transfer for proportional species streams, independent species energy
reference shifts of order 1e24, source/catalog refusals, permutation/reversal,
identity, exact retained density, and late failure/cancellation. The existing
bulk-transfer and caloric regressions remain active.

## Source-backed Ohmic leakage with dielectric memory

`field_storage::dielectric::leaky::ResolvedLeakyDielectric` combines fixed
isotropic multipole Debye polarization with independently selected nonnegative
DC conductivity from the same complete source bundle. Explicit sigma=0 remains
a valid sourced insulating limit; negative, nonfinite or unowned dynamic
conductivity refuses. Consuming decomposition retains the full original
material once, together with the exact conductivity quantity/name.

Midpoint total-current elimination adds dt*sigma/2 to the effective field
coefficient and dt*sigma*E_before/2 to its history. Recovery reports retained
polarization, physical dielectric storage, separate polarization/Joule heat,
and stable total-current work. Conductivity adds no storage. The numerical
coefficient is not physical permittivity. This fixed-condition scalar lane
does not supply anisotropic leakage, temperature feedback, induction,
full-wave validity, breakdown or monotone decay for unresolved time steps.

Six focused cases in tests/dielectric/leaky.rs exercise source retention,
zero conductivity, independent midpoint/Maxwell relaxation laws, second-order
convergence, DC/transient heat separation, biased increments, representable
extreme coefficients, and invalid-source/state refusal. Native integration
validation remains required; these fixtures are not measured qualification.

### Simultaneous partial-enthalpy transport on real-gas networks

`nonideal::energy::parcel::transfer::species::network` advances prescribed
species streams on branching and cyclic maps of original PR parcels. All edges
read the same unchanged donor snapshot; each affected endpoint is reconstructed
once after every species and heat contribution is combined. The pair-selective
API delegates to this same borrowed core, retaining its signed-row admission and
canonical SpeciesId reports without cloning an input map of complete materials.

The original shared stream, heat, control and error types are reused. Work is
bounded by 256 parcels and 4096 submitted species/internal-heat/external-heat
records before allocation. Parallel streams are sorted and combined with
compensated summation. GROSS withdrawal is admitted against each ORIGINAL
species inventory before source joins or caloric evaluation. A receipt cannot
fund an outflow or forward an initially absent species during this interval.

Each active donor's original partial calorics is evaluated once. A receiving
catalog joins only that node and its original incoming donors, preserving all
NASA, transport, pack, critical and known PR pair checks. Disconnected sources
do not impose fictitious compatibility requirements. Idle states remain exact;
heat-only nodes retain their original material and accepted density. Exact
zero-net-mass transport likewise preserves the retained density bits.

Each final composition at final density and old node temperature owns the
caloric anchor. Signed species amounts are grouped before multiplying formation
references or RT; original NASA enthalpy increments and PR partial departures
supply the relative stream terms. The original finite PR caloric solver resolves
the resulting extensive energy once. Local/global mass, individual species and
energy defects retain their absolute budgets, independent of arbitrary energy
reference magnitudes. Nonzero net increments that disappear in floating-point
inventory arithmetic refuse; exactly balanced throughflow remains admissible.

Internal antisymmetric heat and uniquely identified external heat enter that
same reconstruction, separately from transported species enthalpy. The external
thermal owner retains opposite storage and any complete entropy audit. These
prescribed records supply no rate law, separation work, phase flash or general
entropy-production guarantee. Source refusal, budget exhaustion or cancellation
returns no proposed network. Regressions cover dense branches/cycles, joint heat,
independent absolute EOS energy, canonical parallel rows, large species reference
shifts, balanced throughflow, node-local catalogs, original-inventory depletion,
work bounds and late transactional failures.

### Pressure-driven real-gas networks with shared enthalpy storage

`nonideal::energy::parcel::transfer::nozzle::network` samples every authored
converging restriction from one original PR parcel snapshot. Actual pressure
selects each upstream owner independently of endpoint array order. The original
`nozzle_flow_to` owner resolves the same density isentrope, finite enthalpy drop,
back-pressure matching and sonic choking, retaining complete throat source data.
Distinct path IDs denote distinct physical restrictions, including parallel paths.

Bulk nozzle outflow is expanded into its original species amounts and submitted
to the single simultaneous real-gas partial-enthalpy network update. Optional
prescribed selective streams, internal heat and explicit external heat share that
same endpoint reconstruction. Additional species amounts are integrated kilograms
and heat records are integrated joules over the interval. The expansion itself is
bounded by the existing 4096-record limit. No receipt funds another outflow; the
strictest active nozzle donor-fraction bound applies to all outgoing nozzle and
selective mass from that original owner, alongside species-specific gross bounds.

The receiving volume stores the original stagnation enthalpy once. No additional
jet kinetic-energy heat is deposited. Full final pressure, temperature, composition,
acoustic and transport views come from the original caloric/EOS/material owners.
An active path whose reconstructed donor pressure falls below its receiver refuses
the entire trial; transferred mass is not clipped. Closed/equal-pressure paths
retain explicit zero-flow reports. Source failures, malformed paths, inventory or
work exhaustion and late cancellation publish no endpoints.

This is a first-order explicit quasi-steady physical interval on admitted
homogeneous branches. The caller refines the interval and resamples the rates;
there is no automatic time-error estimate, phase flash, valve momentum, pipe wave,
shock or arbitrary selective-stream entropy guarantee. Node energy controls own
the final collective reconstruction; each original nozzle control retains its
existing flow-query admission. Focused tests compare the original single-nozzle
oracle and finite-pair result, audit branched/throughflow storage independently,
combine selective counterflow and heat, verify aggregate donor limits and pressure
crossing refusal, preserve large-reference invariance, and check source/cancellation
and lowered-record boundaries.

## Dynamic closed grinding circuit with retained inventory

DESIGN: `population_balance::circuit::transient` extends the original
`ContinuousGrindingCircuit` with `advance_inventory`. The mill retains finite
size-class dry and elemental masses plus a separate water stock; classifier
and return line have zero inventory and instantaneous return. Constant fresh
feed, residence, selection/daughter parameters and classifier partitions apply
during each interval. For each dry or elemental class,
`dH_i/dt=fresh_i-(S_i+escape_i/tau)H_i+sum_j(b_ij*S_j*H_j)`, with
`escape_i=purge+returned*overflow_i`. Water has its own external escape.

The original positive coarse-to-fine backward-Euler step is shared with the
open mill owner. Gross discharge is H/tau. Every initial/numerical endpoint
passes through the original classifier, and source support is checked using
actual fresh plus returned mill feed, residence and feed solids fraction.
A declared returned/fresh dry-flow ceiling is enforced at every point.
Positive fresh dry flow and an admitted initial wet inventory are required;
this closed-circuit API does not model empty startup or shutdown.

The immutable response retains original model/input/query/settings, initial
and final streams, final inventory, integrated fresh feed, actual mill feed,
gross discharge, underflow, recycle, overflow and purge. Independent mill
and whole-circuit audits check dry, water and every element. Internal recycle
never enters external product; no remainder repairs a stream. Cancellation,
work exhaustion, source failure or failed conservation returns no partial state.
The first-order numerical trajectory and endpoint source checks supply neither
a continuous-time enclosure nor a trajectory-error estimate.

FACT (constructed tests): seven circuit regressions plus 33 existing
population-balance tests pass in the original production-module closure.
Independent two-class one-second inventories are 1000/19 kg coarse,
172000/4009 kg fine and 1000/31 kg water. Other checks cover the original steady
fixed point, zero-return agreement with the open owner, ODE refinement, assay
changes, finite accumulation, actual source/circulation refusals and atomic
cancellation/work limits. The `transient_grinding_circuit` example carries a
copper feed change across three intervals and prints retained/exported Cu and
both balance residuals. These are constructed model checks; no rate-power,
liberation, classifier delay or laboratory-to-plant calibration is supplied.

### Adaptive physical-time evolution of real-gas nozzle networks

`nonideal::energy::parcel::transfer::nozzle::network::adaptive::advance_adaptive`
advances a closed, unforced network using the original pressure/nozzle and
partial-enthalpy snapshot owners. Each comparison uses one full interval and
two half intervals, resampling original pressure, composition, nozzle throat
and donor enthalpy at every interval. Only the conservative fine pair is
accepted. No Richardson-extrapolated material, clipped mass or comparison
stream enters the physical account.

Local error controls act on total and individual species inventories,
temperature, and the extensive coordinate
`U - sum_i m_i*u_i_ideal(T_start)`. The latter is evaluated from original finite
NASA internal-energy increments and the accepted PR departure; arbitrary
species formation datums cancel before arithmetic even as mixtures change.
Relative scales use changes from the accepted state, not absolute caloric
references or initial stock. The first-order step-doubling correction accounts
for unequal representable half spans. This is a local temporal estimate, not
a global enclosure, phase-stability proof or experimental accuracy claim.

Original gross-inventory and pressure-crossing refusals cause bounded interval
reduction. Malformed inputs, unsupported sources, root/caloric failures and
cancellation propagate. A final interval below the floor may be attempted
once and must still meet accuracy. Unresolved pressure equalization or reversal
at the floor returns the original named refusal; it is not passed by clipping.
A failed window publishes no partial trajectory and leaves all original
parcels unchanged.

Explicit trial and accepted-step ceilings bound work. The implementation keeps
three candidate maps and fixed per-path compensated accounts, without storing
an unbounded step history. Canonical path reports retain separate gross mass
and stagnation-enthalpy integrals for both directions, including zero paths;
rejected trials contribute only work counters. A final total/species inventory
audit and a whole-window sum of absolute accepted local caloric residuals have
independent caller budgets. The next suggested interval is explicit.

Focused G1/G3 cases compare the exact original half-step trajectory, audit
absolute storage independently, check the adiabatic donor entropy identity and
temporal refinement, exercise changing mixtures with independent huge formation
shifts, and verify gross-limit subdivision, canonical paths, short final spans,
bounded work, source refusal and late cancellation. This owner supplies no
external heat/species schedule, phase flash, pipe momentum or wave evolution.

### External species and phase energy on real-gas networks

The real-gas selective network now admits explicit external species/phase
boundaries together with internal species streams and internal/external heat.
`exchange_candidate_with_boundaries` and
`exchange_candidate_with_boundaries_and_external_heat` share the original
snapshot, gross-inventory admission and single finite-PR reconstruction.
Boundary owners and species consume the same 4096-record work bound. External
GROSS withdrawal joins every internal outgoing stream before caloric evaluation
or catalog mixing; incoming external material cannot fund same-interval outflow.

`RealGasSpeciesNetworkBoundary` retains its original PR catalog and the neutral
`BoundarySpeciesTransfer` chart:
`E_into = net_mass*h_ideal_species(T_original) + relative_energy`.
The external owner authors its complete phase or uptake energy. Negative mass
does not silently acquire a gas departure or latent heat. The source-bound
`RealGasBoundaryStream::inflow` and `outflow` conveniences instead supply actual
original gas partial enthalpy: finite NASA enthalpy differences plus the donor's
analytic PR partial departure, with formation references cancelled before
arithmetic. Pressure-flow work is included once. The inlet helper's receiving
original snapshot defines its reference temperature and must be the same node
subsequently supplied to the network.

Every affected node retains its original source and its actual active incoming
or external catalogs, including NASA, transport, critical-constant and known
binary-interaction compatibility checks. Unrelated idle nodes impose no source
requirements. External net mass/species amounts are subtracted from open-system
conservation defects; their relative energy enters the same grouped-reference
caloric target as all internal streams. The neutral boundary report explicitly
retains the ideal-species reference enthalpy and the complete authored energy.
Original-density preservation and unrepresentable-inventory refusal remain.

External donor/phase storage, its aggregate inventory allowance, entropy and
any interfacial kinetic law belong to the enclosing physical owner. The gas
boundary does not create a reservoir EOS, infer a phase transition or deposit
latent heat again. Zero rows preserve states; source/domain/work failure and
late cancellation return no proposed network. Seven focused cases compare
actual-gas inlet/outlet against the internal PR transfer owner and independent
stored energy, combine phase uptake with both heat channels, shift independent
species references, introduce source-bound new species, combine gross depletion,
and exercise canonical order, work bounds and transactional refusals.

## Finite classifier-feed sump and dilution water

DESIGN: `population_balance::circuit::transient::sump` adds
`ContinuousGrindingCircuit::advance_with_sump`. A `ClassifierFeedSump`
retains an explicitly sourced or authored positive common drain time. Mill
and sump carry separate size-class dry, elemental and water inventories;
independently declared makeup water enters the sump. The original
instantaneous-return API remains a separate zero-hold-up model.

Each backward-Euler step uses a positive per-class Schur complement to reuse
the existing coarse-to-fine mill step. The sump stock is reconstructed
independently. Actual sump discharge feeds the original classifier; actual
fresh-plus-return feed enters the original mill's source query. Eight
integrated stream accounts and independent mill, sump and whole-circuit
conservation audits retain both stocks, makeup water and internal return.
No dense circuit solver or balancing stream is introduced. Inputs remain
immutable on source/domain, arithmetic, budget or cancellation refusal.

FACT (constructed G0/G1/G3/G4 checks): seven dedicated regressions and all 47
population-balance tests pass in the original production-module closure.
The checks include an independently solved rational six-by-six system,
exact open-mill final-stock agreement at zero return, the original steady
fixed point, refinement against a continuous two-eigenvalue solution,
actual source/circulation restrictions and atomic work/cancellation refusal.
The `grinding_sump` example demonstrates simultaneous grade and dilution
changes; unsupported dilution refuses before successful stdout.

The declared drain law applies the same residence to solids, elements and
water, with fixed controls during each interval. It supplies no settling/slip,
pump or level model, tank-volume inference, empty-sump operating policy,
plant calibration or continuous-time error enclosure. Numerical conservation
does not certify the first-order trajectory's discretization error.


## Irreversible ductile initiation from native plastic history

`ductile` owns local damage initiation as
`omega += delta_alpha / epsilon_f(eta, delta_alpha / dt, T)`, with native
accumulated equivalent PLASTIC strain, tensile-positive physical stress
triaxiality and an explicit physical duration. Endpoint-held conditions make
this first-order quadrature for changing stress/rate/temperature histories.
Elastic unloading preserves utilization; reversed plastic flow adds to it.
History is immutable until its caller accepts a complete candidate. Overshoot
is retained and omega >= 1 means initiation, not zero carrying capacity.

Explicit Johnson--Cook parameters use `exp(-d3 * eta)`, whose d3 sign is the
opposite of the original publication's convention. Positive finite factors
are admitted over the complete rectangular domain. Original tabulated loci
retain every knot and use binary-search piecewise-linear evaluation with no
extrapolation. Rate and triaxiality are checked when plastic flow occurs;
temperature support applies on every step. Physical Cauchy or Kirchhoff stress
is required, never first Piola components or stress relative to backstress.
Scaled invariant evaluation avoids unnecessary overflow; undefined
hydrostatic stress with positive plastic flow refuses.

`ductile::source::resolve_ductile_locus` requires the exact property
`ductile-fracture-equivalent-plastic-strain` against
`stress-triaxiality`, dimensionless typed coordinates, original positive
`equivalent-plastic-strain-rate` bounds in 1/s, and a semantic absolute
temperature. It retains original card/state identity, selected claim, fixed
context, every original knot and verified query receipt. It authors numerical
interpolation without changing the source's interpolation policy or promoting
intermediate values to measurements. Engineering elongation, reduction of area
and fatigue ductility are not aliases. Unsupported evolving coordinates refuse.

A resolved one-dimensional locus holds the selected temperature fixed and
inherits the original rate applicability interval; it does not fit dependence
on either coordinate. Source/native binding requires the same original card,
named material state and complete fixed context in both directions, excluding
only the actual runtime triaxiality and plastic rate. Native-only environmental
coordinates cannot silently acquire fracture calibration.

Eight core tests and seven source regression definitions cover independent
formula/analytic initiation, signed-factor admission, path irreversibility,
rotation/scaling, restart, source meaning, typed domains, ambiguity pins and
original material/context binding. The original initiation and source tests
now pass alongside the Lode additions: 20 native core tests and 19 material-source
integration tests execute against the original numerical/source bodies and
genuine dependencies. These focused builds are not a full-crate or workspace
gate. Initiation alone adds no post-initiation softening, fracture energy,
characteristic length, crack topology, mesh-objective failure or experimentally
validated component failure. Progressive small-
strain evolution is now provided by the separate owner below; broader fracture
tasks remain open.

## Hyperelastic models with explicit Cauchy tensile limits

`strength::ResolvedTensileHyperelastic::from_cauchy_limit` binds the original
`ResolvedHyperelastic` and an explicit isotropic tensile limit from the same
immutable material card and complete typed query. The exact constitutive model
must be a member of that card; source identity covers model selection, strength
property and conditions. Energy, Piola stress and tangent delegate to the native
owner unchanged. The constructor explicitly selects a current-area Cauchy-strength
interpretation; it does not convert engineering UTS or supply a damage law.
Changing conditions requires readmitting both original owners. Tests in
`tests/tensile_hyperelastic.rs` cover unchanged native response, model/card/query
mismatches, law-domain refusal and distinct source bindings.

## Pressure-sensitive porous-solid compaction

`plastic::porous::PorousPlasticity` adds the bounded infinitesimal-strain,
associated Deshpande--Fleck self-similar ellipse
`phi = sqrt((q^2 + a^2 p^2)/(1 + a^2/9))`. Mean stress is tension-positive;
the authored initial yield stress is the uniaxial value. Elastic moduli,
nonnegative energetic hardening `H`, and dimensionless pressure sensitivity
`a` stay fixed through a load path. Positive sensitivity permits hydrostatic
yielding and permanent volumetric plastic strain. Explicit `a = 0` recovers
isotropic-hardening J2, including its isochoric plastic history.

The existing `SmallStrainLaw` owner supplies native solid assembly with a
fixed-array scalar safeguarded return and its analytic discrete tangent in
tensor-Voigt coordinates. `response` returns typed refusals without mutating
input history. It validates publicly mutable coefficients, finite elastic and
total strain norms within the smaller of the elastic limit and 0.05, and the
necessary dual-norm bound between retained plastic offset and accumulated
work-conjugate strain. Unrepresentable or unconverged returns refuse; the
infallible stress/tangent interface marks a refusal nonfinite and preserves
history.

`advance_with_energy` validates the old committed elastic stress and the new
yield/reconstruction residuals against authored allowances. It separates
elastic and energetic-hardening storage, physical yield heat
`yield_stress * delta_alpha`, and backward-Euler numerical dissipation.
Integration loss is never relabeled as physical heat. No parameter-change
work or thermal feedback is inferred.

This law is symmetric in hydrostatic tension/compression. It does not model
evolving porosity, large-strain densification, pore-fluid pressure, soil caps,
asymmetric tensile damage or measured foam calibration. Its reference is
Deshpande and Fleck (2000), JMPS 48, 1253--1283, equations 2, 3, 11, 13 and 21
(<https://doi.org/10.1016/S0022-5096(99)00082-4>).

Focused numerical regressions cover hydrostatic and pure-shear closed forms,
all six tangent columns after mixed history, weighted tangent symmetry, the
native J2 limit and reversal, rotations, stress-unit scaling, refusal paths,
and the physical/numerical work split. These are numerical checks, not
experimental validation.

## Source-resolved porous compaction

`plastic::porous::source` resolves Young's modulus, Poisson ratio, initial
uniaxial yield stress, explicit porous energetic hardening and dimensionless
Deshpande--Fleck pressure sensitivity at one fixed material state.
`ResolvedPorousPlasticity` preserves the complete source/query/usage bundle
and delegates `SmallStrainLaw` and checked responses to the native return map.
The explicit source keys are `porous_hardening_modulus` [Pa] and
`porous_pressure_sensitivity` [dimensionless]. Hardening is the slope against
this law's work-conjugate accumulated plastic strain; density, porosity, J2
coefficients and total-strain tangents provide no fallback. Perfect plasticity
and the J2 limit each require an explicit zero source.

Fixed external-condition interpolation is admitted by the common material
resolver, while evolving strain, stress, history, rate and frequency sources
refuse. Exact SI quantity dimensions and original native parameter limits are
checked even when adopting a separately resolved shared material bundle.
The native small-strain and symmetric-yield limits remain binding; source
re-resolution neither transfers history nor accounts for parameter work.
Additional physics properties and their original receipts remain available.

Four focused regressions cover source-driven hydrostatic compaction and mixed
loading, fixed-temperature interpolation with extra physics and receipts,
missing/wrong/other-law calibration refusal, explicit zero limits, ambiguous
source pinning, and refusal to freeze evolving mechanical conditions.


## Progressive ductile J2 softening with explicit work and length

`ductile::evolution::DuctileJ2` composes the ORIGINAL native J2 radial return,
energetic isotropic/Prager hardening, native work audit and algorithmic tangent
with irreversible exponential attenuation. After initiation,
`g = exp[-h/Gp * integral (sigma_y + H_iso*alpha) d(alpha)]`. The receiving
geometry supplies the original fixed characteristic length h. Gp is explicitly
POST-INITIATION EQUIVALENT PLASTIC WORK PER AREA; it is not total thermodynamic
fracture dissipation, cohesive Gc or toughness. The scalar degraded-work
diagnostic is `Gp*(1-g)`; it is not added again to mechanical dissipation.

The tensile volumetric and deviatoric response degrade while compressive bulk
stiffness remains intact. The stored potential is
`g*psi_native + (1-g)*K/2*min(trace(strain),0)^2`. Physical stress and tangent
differentiate this split consistently. The tangent includes native accumulated-
plastic-strain sensitivity and moving onset, including triaxiality, Lode parameter
and rate effects during the crossing increment. Initiation remains endpoint
quadrature;
after its first accepted crossing the initiation observation freezes, while
the original native plastic clock continues. Supported temperature is required
even during elastic unloading. The tangent may be nonsymmetric and indefinite.

Endpoint external work closes into stored-energy change, native plastic
dissipation, released stored energy from damage, and backward-Euler numerical
loss. Numerical loss is not plastic heat or fracture energy. This discrete
identity does not certify an unresolved physical trajectory. An explicit
maximum damage increment refuses excessive trials without clipping damage or
the tangent; native continuation must refine the physical loading interval.
There is no residual stiffness floor. Underflowed integrity and nonfinite
work refuse. Trial updates are fixed-size and allocation-free, preserve the
old state, and bind history to the exact immutable material/length owner.

`ductile::evolution::source::ResolvedDuctileJ2` retains the complete original
`ResolvedJ2Plasticity` source bundle and the same-card initiation calibration.
It requires exactly `ductile_post_initiation_equivalent_plastic_work_per_area`
with scalar J/m2 dimensions, positive original data and supported fixed
external conditions. Engineering elongation, toughness and other fracture-work
properties are not aliases. Card/state identity and the complete fixed context
must match in both directions. Geometry supplies length explicitly; calibration
cannot be overridden through numerical settings. Prepared runtime updates
make no material-database queries or source-card copies.

FACT (bounded actual-source execution): all 20 core tests and all 19 material-source
integration tests pass. They exercise intact/native agreement, independent
hardening-work and Lode-dependent onset, all 36 tangent entries through coupled
initiation and compressive reclosure, length/work scaling, irreversibility,
plastic reversal work balance, repeatability, transactional refusals and
numerical-loss decay on sufficiently refined paths. The exact held-strain
restart preserves accepted state with zero incremental work; nonzero loading
still obeys the original rate support. Source tests exercise original receipts,
complete typed context, physical attenuation and source/domain refusals.

The build uses the original native J2, elasticity, tensor, energy, ductile and
source bodies with genuine numerical and material-data dependencies. It selects
those modules without unrelated exports; this is not a full-crate, workspace or
integrated structural gate. Structural execution is reported by fs-solid.

This owner advances the post-initiation part of MR63
(`frankensim-material-reality-ui2if.7.11`) without closing that broader bead.
Characteristic-length scaling does not establish arbitrary-mesh crack-path
objectivity. Finite deformation, crack topology, new-surface contact, erosion,
fragments, history transfer, heat partition and experimental
qualification remain outside this model.

### Lode-dependent initiation and progressive softening

The shared `ductile::DuctileStressState` carries tensile-positive triaxiality
and `xi = (27/2) det(s) / q^3 = cos(3 theta)`. Both are invariant under rigid
rotation and positive stress scaling. This convention is distinct from an
angle in radians or `1 - 6 theta/pi`. Separate total-stress/deviator scaling
avoids unnecessary overflow and loss of small deviators. Exactly hydrostatic
stress supplies no invariant pair; positive plastic flow there still refuses.

`DuctileInitiationLaw::fracture_strain_at` propagates the complete stress state
through the original history owner and every consumer of that owner. Existing
Johnson--Cook and one-dimensional tables keep their triaxiality-only response.
Lode-dependent laws refuse the scalar entry point instead of inventing xi.
The progressive J2 tangent includes the exact polynomial invariant derivative
`d xi = [27/2 dev(U^2) - 9/2 xi U] : d sigma / q`, where `U = s/q`.
It introduces no Lode-angle singularity at the tensile/compressive meridians,
and it retains the original native return map and endpoint onset quadrature.

`ductile::source::resolve_ductile_stress_locus` prepares the exact named
fracture-plastic-strain property from an original `Surface2D` claim on
`stress-triaxiality` and `stress-lode-parameter`. Either source-axis order is
accepted. Original `LinearInside` triangles, samples, full typed query points
and verified receipts remain attached. The existing `PreparedSurfaceResponse`
owns values, exact support predicates and cached analytical cell gradients;
holes and isolated samples gain no solver support from their bounding box.
Shared edges use its explicit first-cell semismooth branch. Runtime evaluation
performs no database query or source selection. Temperature remains fixed at
the original observation, and the positive plastic-rate interval remains an
applicability statement, not fitted rate dependence.

The original `ResolvedDuctileLocus`, source-bound `DuctileHistory`, and
`ResolvedDuctileJ2` bindings accept both curve and surface preparations.
Card/state identity and fixed external context still match in both directions;
only actual runtime stress invariants and plastic rate are excluded from that
fixed-context comparison. Synthetic regressions cover equal-triaxiality stress
states, rotation/scaling, all tangent entries, source holes and semantics,
native softening, and assembled tetrahedral utilization. The 20 native core
and 19 material-source tests execute successfully against the original bodies
and genuine dependencies. Two assembled source-bound tetrahedron regressions
also pass through the original J2, history and source-binding owners. These
tests do not qualify an alloy or supply experimentally measured Lode calibration.

## Finite ductile softening with original multiplicative plastic history

`plastic::finite::ductile::FiniteDuctileJ2` retains the ORIGINAL `FiniteJ2`
inverse plastic metric, Kirchhoff return, linear/Voce hardening, logarithmic
derivatives and explicit native domains. Shared initiation receives the actual
undegraded Kirchhoff triaxiality, Lode parameter and plastic rate. After onset,
`g = exp[-h/Gp * integral R(alpha) d(alpha)]`, where
`R(alpha) = Y + H*alpha + Q*(1-exp(-b*alpha))`. The Voce term uses the original
stable translated hardening primitive. The reference band width h stays fixed;
Gp is equivalent KIRCHHOFF plastic work per REFERENCE area after initiation.
`Gp*(1-g)` is a diagnostic, not an additional physical dissipation term.

With `psi_c = K/2 * min(log(J),0)^2`, storage is
`psi = g*psi_native + (1-g)*psi_c` and physical first Piola stress is
`P = g*P_native + (1-g)*P_c`. Compressive bulk strength remains intact.
The full analytical 9-by-9 `dP/dF` includes the native return, moving initiation
and the complete stress-invariant/rate damage chain. Its generally nonsymmetric
tangent is NOT a scalar-potential Hessian: this law implements no
`FiniteStrainPotential`. The stress-only residual path omits tangent assembly.

All energy accounts refer to original volume. Physical plastic loss is
`g_new*D_native`; damage release is
`(g_old-g_new)*(psi_native_old-psi_c_old)`, evaluated from positive active
storage. Both stay nonnegative. Independently assembled endpoint work closes as
`P:delta(F) = delta(psi) + plastic_loss + damage_release + endpoint_defect`.
The endpoint defect is SIGNED, including finite rotation of a prestressed body;
it is not heat, fracture dissipation or an assumed nonnegative numerical loss.
Consumers must bound its absolute magnitude and refine the physical path.
Exact holds preserve accepted history and all work accounts bit-for-bit.
Foreign owners, unsupported domains, excessive damage increments, nonfinite
results and unresolved integrity/history increments refuse without mutation or
clipping. Initiation freezes after its first accepted crossing; native plastic
history continues. No stiffness floor or element deletion is introduced.

`plastic::finite::ductile::source::ResolvedFiniteDuctileJ2` retains the original
`ResolvedFiniteJ2`, work bundle and `ResolvedDuctileLocus`. The finite model must
be an actual member of the same immutable `MaterialCard`; equal coefficients
or source names do not prove membership. Complete fixed typed contexts agree.
The exact positive scalar J/m2 property is
`finite_ductile_post_initiation_equivalent_kirchhoff_plastic_work_per_reference_area`.
Small-strain Gp, current-area work, cohesive Gc and total fracture energy cannot
alias this ordinate. Original claims, observations and usage receipts remain
attached; prepared trials perform no source query or calibration substitution.

FACT (bounded actual-source execution): seven new numerical tests and four new
source tests pass in `tests/plastic_finite_ductile.rs` and
`tests/plastic_finite_ductile_source.rs`. They cover intact native equivalence,
independent linear/Voce work, all 81 tangent entries through Lode/rate onset and
active compressed flow, objectivity, signed rotation work, retained bulk
strength, irreversibility, exact holds, refinement and transactional/source
refusals. Eight original finite tests and 20 shared ductile core tests also
pass with the reused derivative owners. Execution selects original source
modules with genuine dependencies; it is not a full Cargo crate/workspace gate
or experimental material qualification. Structural execution has its own
fs-solid contract. The native backbone is fixed-temperature rate-independent
plasticity; initiation may depend on actual rate and temperature. Native creep
backbones, finite thermomechanical feedback, crack topology, separated-surface
contact, fragments and history transfer are outside this owner.

## Source-bound tensile initiation and elastic handoff

`brittle::source::initiation` binds the native tensile detector and onset to the
complete original `ResolvedBrittleCrackBand`. The detected plane may be resolved
at a receiving geometry's declared crack-band width without requerying or
replacing source coefficients. Handoff requires the same original source bundle,
all fracture parameters, strain domain and material-work gate, then delegates
native stress/storage continuity and undamaged-history checks. Additional
multiphysics receipts remain attached; an equal-valued foreign material or state
cannot substitute. This bridge supplies no independent crack geometry,
existing-history reset, all-plane mixed-mode criterion or experimental
calibration. Focused regression definitions cover physical plane discovery,
source/context mismatch, temperature-dependent onset and bounded refusals; Rust
parsing/formatting passed, while compiled execution remains outstanding.

## Finite porous compaction and native potential coupling

`plastic::finite::porous::FinitePorous` adds multiplicative
pressure-sensitive plasticity with quadratic Hencky elastic storage and
linear or explicitly selected Voce energetic hardening. Its coefficients act on Kirchhoff
stress; energy and yield heat are per original reference volume. An authored
finite law is required: a small-strain or Cauchy-yield calibration is not
silently reinterpreted.

The accepted SPD inverse plastic metric retains its determinant. Trials use
`b_trial = F Cp^-1 F^T`, the original finite spectral/logarithm owner, and
the same scalar pressure-sensitive return kernel as `plastic::porous`.
The corrected metric is `F^-1 exp(2 e_returned) F^-T`, without determinant
normalization. The law reports both total and plastic volume ratios.
History admission includes finite SPD resolution, explicit elastic/plastic
log-strain and accumulation bounds, positive total Jacobian, and the
necessary weighted SPD-distance reachability bound. Hydrostatic compression
can retain a finite permanent volume loss; the zero-sensitivity limit
retains native finite J2 behavior.

The analytic full Piola tangent differentiates the same return through the
Frechet matrix logarithm, including the inverse-deformation geometric term.
The reduced incremental potential is stored energy plus physical yield
dissipation. `FrozenPorousTrial` implements the existing
`FiniteStrainPotential` interface, borrowing one accepted history through
all Newton and line-search evaluations. A trial never commits history.
Callers admit all material points and their full coupled step before
publishing returned candidates. The endpoint Piola work defect is distinct
from physical heat; this law does not assert that the finite-step defect is
nonnegative.

This is isothermal, rate independent, symmetric in mean Kirchhoff stress,
and bounded by explicitly supplied applicability limits. It does not infer
evolving porosity, cell-wall densification, pore-fluid pressure, tensile
damage, or experimental foam accuracy. No mixed-pressure or locking-free
element claim is added. Numerical G0/G3 fixtures cover independent
hydrostatic volume/storage/heat, all nine tangent columns and potential
gradient after mixed history, objectivity and exact holds, native J2,
logarithmic small-strain correspondence, unit scaling and refusals.

`FinitePorous::with_voce_hardening(Q,b)` adds the same native finite-J2
saturation `Q*(1-exp(-b*alpha))`, with Q an added strength in Pa and b per unit
work-conjugate accumulation. Its stable primitive is recoverable storage;
physical yield heat remains `Y0*delta_alpha`. The existing pressure/deviatoric
return uses a bracketed scalar solve and its exact nonlinear derivative,
including the stable `Y-lambda*Y'` curvature term. The zero-Q case retains the
linear path. No new history variable, tensor Newton solver, source query or
allocation is added to the constitutive return. G0/G1 cases in
`plastic_finite_porous_voce` check independent hydrostatic compression and
unloading, all nine tangent columns, incremental-energy gradients, native
Voce-J2 equivalence at zero pressure sensitivity, scaling and refusals. This
is isotropic saturation, not evolving pore fraction or cell-wall densification.

### Source-defined finite porous compression

`plastic::finite::porous::source::ResolvedFinitePorous` binds the explicit
`finite-porous-hencky-kirchhoff` model card (law/state version 1) to this same
native owner. All four pressure-valued coefficients, dimensionless ellipse
sensitivity, and four native deformation/history limits are mandatory. Source
identity, validity and the complete typed fixed context remain available.
Virgin history means the identity inverse plastic metric, including its volume,
and zero accumulated work-conjugate strain. Unknown parameters, nonvirgin
initialization, evolving mechanical/rate coordinates, or a Cauchy/small-strain
law selector refuse. No measured stress-strain table is silently converted.

The separate `finite-porous-voce-hencky-kirchhoff` selector additionally requires
`voce_saturation_pa` [Pa] and dimensionless `voce_rate_per_strain`, retaining the
same law/state version and virgin-history convention. The linear selector
rejects these extra coefficients. Source resolution invokes the native
saturation constructor; structural consumers retain the complete selected law.

`FinitePorousMaterial` lets structural consumers borrow either this source or
the authored native law. Source admission occurs once; existing frozen-history
potentials, return mapping, analytic tangents and energy accounting remain the
hot path. G0/G1 tests in `plastic_finite_porous_source` cover independent
hydrostatic crush, permanent volume loss and elastic reload, the Kirchhoff/
Cauchy distinction, physical yield heat, explicit J2/perfect-plastic limits,
source domains and propagation of native bounds. This adds no calibration,
thermal evolution, porosity evolution or cell-wall densification model.

## Original density for brittle reference-mass dynamics

`brittle::source::dynamics` resolves original fracture coefficients and positive
density together through the existing material-state owner. The density getter
admits only the exact scalar SI bulk-density property at fixed external
conditions and independently enforces positivity on adopted multiphysics
bundles. Dynamic consumers bind that snapshot to original reference geometry
and retain its mass inventory; later density queries do not authorize matter or
momentum changes. No current-to-reference conversion, porosity evolution or
mass-transfer law is inferred. Quasistatic brittle resolution keeps its existing
requirements. Four regression definitions cover complete source ownership,
quantity/positivity refusals, thermal snapshots and mechanically restricted
density. New Rust parsed/formatted; compiled execution remains outstanding.

## Source-integrated Fourier profiles at gray conducting gaps

`thermal::conduction::GrayLayerGap::equal_area_full_view` explicitly declares
two equal-area, opaque diffuse-gray faces with unit mutual view factor and a
fixed nonradiative gap HTC. `TwoLayerConduction::evaluate_gray_gap` resolves
their two material-face temperatures using the original k(T) Kirchhoff
integrals and inverses. The ordinary contact and gray-gap paths share the
existing supported scalar flux bracket and safeguarded iteration. No extra
temperature profile, thermal timestepper, wall storage or source query is added.

The interface relation is `j = h_g*(Ua-Ub) + epsilon_eff*sigma*(Ua^4-Ub^4)`;
each half cell satisfies its own integrated Fourier equation at that same j.
The full analytic bulk-temperature derivatives include the ratios of bulk
to surface conductivity. Individual source spans need cover only their own
material profiles, including when the gap bridges disjoint temperature spans.
Disabling radiation delegates to the original contact-resistance owner and
retains its exact result and near-perfect-contact representation gate.

`ImplicitLayerLink::with_gray_gap` binds this complete response to the original
chain/graph enthalpy solve, with the original caloric dT/dh chain rule. One
equal/opposite heat transfer enters its energy and entropy accounts; reported
conductive and radiative components are descriptions of that same transfer.
The attachment requires zero authored contact resistance, no separate contact
law, equal original covered areas and planar geometry. Source layers retain
their explicit pure-phase restriction. Radial radiation geometry, nonunit view
factors, phase fronts inside a half cell, evolving emissivity and gap storage
remain outside this closure.

Regression definitions compare independent nonlinear material integrals and
surface balances, the original constant-k gray link, exact radiation-disabled
contact recovery, reversed ordered materials, full flux tangents, and a
manufactured original-source chain/graph endpoint with conserved heat. Source,
ownership, nonlinear work and cancellation refusals retain the original state.


## Source-integrated thermal eigenstrain in native tensile fracture

`brittle::source::expansion::ResolvedThermalBrittle` binds the original brittle source bundle to the original integrated isotropic thermal-expansion source. Its resolver selects the fracture, elastic and instantaneous alpha(T) claims together at one fixed reference, then retains those claims and integrates the selected expansion through every supported temperature cell. Missing path support, different evidence, coordinate conventions, source context or reference refuses. Negative and zero expansion coefficients are admissible; alpha is never replaced by an endpoint secant approximation.

`ResolvedThermalBrittle::with_resolved_band` prepares a discovered crack plane
and receiving geometric width while retaining the exact integrated expansion
and immutable reference. The complete material bundle, including any original
density receipt, elastic/fracture coefficients and strain/work bounds must
agree. Equal numerical coefficients from a different source or temperature
context do not authorize replacement. Geometry preparation transfers no
history: the original tensile owner must separately admit an elastic seed and
its stress/storage continuity, and existing damaged states keep their native
ownership checks. Three native regression definitions in
`tests/brittle_source/expansion/initiation.rs` cover source/reference retention,
mechanical-strain fracture after handoff, unchanged later thermal refresh,
source/context/control refusal and preservation of additional density receipts.

The borrowed total-strain law supplies mechanical strain equal to total tensor strain minus the integrated free strain times the identity. Both total and mechanical Frobenius bounds are checked. Tensor shear values, native tangent conventions, fixed plane/width, cohesive history and original fracture-energy matching are retained. No independent damage model or thermal history copy is introduced.

A coefficient-only preparation transfers the original mechanical history to the admitted successor source. A solid consumer then solves total displacement and the new eigenstrain together; it must not precommit a fully restrained damage trial before free expansion equilibrates. A distinct optional held-total-strain operation represents actual restraint: it orders coefficient work before the native mechanical-strain change and uses the original exact mean stress for eigenstrain work. Signed coefficient/eigenstrain work, coefficient-induced and strain-induced fracture, storage change and checked work defects remain separate. Cooling, reheating and increased toughness do not erase original damage or cumulative fracture work.

Six focused regression definitions cover exact multi-knot expansion integration, free versus restrained cooling and reheating, separate signed coefficient/eigenstrain work and native stress/tangent/shear equivalence, negative expansion, source/reference/context/strain-domain refusals, and rejection of numerically identical expansion from different selected evidence. Pinned Rustfmt parsing and bounded API/physics review completed; native compilation and runtime execution remain unperformed in the full-filesystem environment. This is an explicitly resolved, coefficient-first small-strain model: it does not supply finite thermal stretches, evolving crack planes, thermal-cycle integration, a caloric equation, inferred fracture heat or entropy production.


## Original-source mass and mean forces for thermal tensile motion

The thermal brittle source now exposes `source::expansion::dynamics::resolve_thermal_brittle_dynamic_material`. It extends the existing shared reference/current resolver with strictly positive density in the same original selected fracture, elastic and instantaneous-expansion bundle. Static thermal resolution retains its existing requirements. Current-temperature density queries are explicit material values: a receiving motion owner binds its chosen initial snapshot to reference geometry once and must retain that inventory. The expansion reference does not silently infer a current/reference volume conversion or a mass/momentum source.

`ResolvedThermalBrittle::mean_small_strain_law` supplies the original analytic path-mean stress after total-to-mechanical strain conversion. The checked trial retains the original optional mean derivative; the trait preserves a missing derivative for an explicitly admitted native solver fallback. Force-only assembly does not compute an unused derivative. State updates retain the original unaveraged endpoint history and physical fracture law. The additive thermal strain has zero displacement derivative at fixed new temperature, so the displacement tangent remains the original native mechanical mean tangent with unchanged tensor-shear conventions.

The material force is work-conjugate to both total displacement strain and free thermal strain. Their separate work ports sum to original mechanical-strain path work; a coupled consumer must retain the old thermal source to evaluate that split. This adapter supplies neither its own time integrator nor a heat/entropy model. Coefficients are transferred first, then total/eigenstrain follow one admitted straight increment under the new coefficients; no second-order thermal trajectory or unseen-reversal claim follows.

Four regression definitions cover complete original density/expansion/fracture source binding; native mean work with independent full total-strain tangent columns; matched free thermal motion versus restrained damage and reheating; and density, source-owner, strain-domain and nonsmooth-derivative refusals. Changed Rust parsed/formatted with pinned Rustfmt and received bounded source/API review. Native compilation and runtime execution remain unperformed in the full-filesystem environment.

## Original IF97 equilibrium water density and high-pressure acoustics

`phase::water::acoustics::HomogeneousEquilibriumWater` retains an original
saturation pair and an explicit rapid-phase-equilibrium closure. The original
Regions 1/2 entry remains available; the regional entry also consumes prepared
Region 3 liquid/vapor endpoints below the original critical temperature. Density
comes from phase-mass-weighted specific volume. Compressibility and sound speed
follow the original Region 4 saturation slope and the reversible-work tangent
`du + p dv = 0`, with no fitted mixture equation or resampling of the EOS.

Dense-phase pressure inversion residuals remain visible. The phase caloric
tangents include the difference between each achieved phase pressure and the
common coexistence pressure, and latent work uses the original
`delta_u + p delta_v`. Original entropy, Gibbs, Clapeyron, saturation-pressure
and latent-enthalpy consistency defects are reported without shifting source
properties or treating numerical inconsistency as physical entropy production.
Both phase masses must be positive; exact endpoints refuse an ambiguous
equilibrium tangent. The caller owns rapid thermal/phase equilibration,
negligible slip and the validity of this acoustic-frequency limit.

Evaluation is fixed-size and allocation-free after the caller prepares the
source pair. The focused regressions define original low-temperature parity,
dense-water reversible-work finite perturbations, phase-inventory density,
source/pressure-defect retention, energy closure and invalid-quality controls.
Independent source-equation calculations checked the three IAPWS Table 33
reference states and nine high-pressure coexistence states at 630, 640 and
646 K. Native Rust compilation and execution of the new regression definitions
remain unperformed because the execution environment became unavailable.

### Acoustic views of solved regional heating and boiling states

`phase::water::acoustics::thermal` retains an original `RegionalWaterState`
and an explicit frozen-phase-mass or homogeneous-equilibrium closure for
coexistence. Single-phase states use their original IF97 response; exact
saturation endpoints admit the explicitly frozen response and refuse an
ambiguous equilibrium tangent. The adapter preserves the thermal state's
density exactly and the selected closure's bulk modulus, retaining the
underlying response and original thermal/source diagnostics. It implements
the existing `FluidAcoustics` interface without another source query,
pressure solve, interpolation, phase-kinetic model or implicit closure choice.

Preparing another acoustic view from the next accepted thermal state updates
density and compressibility together. The view does not modify the thermal
state, phase inventory or conservation account. The added regression
definitions exercise heated-state density, source retention, distinct
coexistence responses, high-pressure and single-phase dispatch, exact phase
limits, and use through the shared acoustic-consumer interface. These Rust
regressions remain unexecuted in the unavailable execution environment.

### Thermal longitudinal reduction of complete anisotropic piezoelectric sources

`ThermalLongitudinalPiezoelectric::from_tensor` binds an original affine
strain-charge tensor source to the existing longitudinal rod/circuit thermal
owner. It retains all original compliance, signed coupling, permittivity and
density roles, their source selections, fixed frame transform and oriented
specimen/electrode direction. Complete tensor contractions include the lateral,
shear and off-diagonal temperature derivatives. Both full coupled endpoint
blocks must pass the native stability checks, including directions not excited
by the longitudinal projection; affine source blocks then remain positive
throughout the admitted interval.

Every temperature trial reconstructs the original source tensor and its
longitudinal reduction. No scalar source claims are synthesized. A matching
original query preserves the original resolved snapshot. Caloric binding checks
the reduction receipt as well as source values, so a changed frame or electrode
orientation cannot silently reuse another thermal law. Signed coupling may
cross zero; projection roundoff is admitted against original source magnitudes
without modifying or clamping actual material values.

The existing free-energy, entropy, internal-energy and field heat-capacity
responses serve the rod/circuit and finite-caloric consumers. The scalar
constructor continues to require scalar source roles. This model retains
small-signal homogeneous stress and longitudinal-field constraints; it does not
supply a finite-specimen insulating lateral-boundary solve, thermal expansion,
pyroelectric offsets, hysteresis, depoling or arbitrary non-affine source paths.
Affine stress-charge sources are refused because their inverted strain-charge
coefficients are generally non-affine.

Four synthetic regression definitions exercise independent full-tensor
contractions and thermodynamic derivatives, signed coupling through zero,
corotation and electrode reversal, incompatible reduction receipts, and
unstable or non-affine unexcited source directions. Source/API review and
independent analytic calculations were performed. Native Rust compilation and
execution remain unverified because the execution environment became
unavailable.

### Finite-rate density and compression response of regional water

The source-retaining
phase::water::acoustics::relaxation::RelaxingRegionalWater binds the original
RegionalWaterState to its existing frozen and homogeneous-equilibrium acoustic
responses. Its real density is exactly the thermal density. One declared
positive fixed-compression relaxation time supplies a positive Maxwell arm
with delta_K = K_frozen - K_equilibrium. The shared BulkRelaxationArm response
also serves the original molecular bulk owner; no second Debye evaluator,
new EOS query or source-property adjustment is introduced.

With s=-i*omega, K=K_equilibrium+delta_K*s*tau/(1+s*tau). Compressibility,
density susceptibility, forward wavenumber and pressure/velocity impedance
all follow this same K and rho. The corresponding fixed-pressure compliance
time is tau*K_frozen/K_equilibrium. The material evaluates DC explicitly,
retains both endpoint owners and refuses unavailable coexistence tangents,
invalid times/frequencies and nonpositive relaxation strength. Evaluation
uses fixed-size scalar arithmetic without successful-path heap allocation.

This is an explicitly selected single-time constitutive model. IF97 supplies
the density and endpoint moduli, not the relaxation time or its validity band.
The caller owns the rate data, small-amplitude coexistence range and common
pressure/velocity homogenization. The model does not evolve phase-specific
temperatures, vapor quality, bubble inertia, slip or noncondensable gas, and
adds no molecular, wall or viscothermal loss. Source consistency defects stay
with their original owners and do not become physical dissipation.

Focused native validation compiled the exact changed source against genuine
cached dependencies. All four new finite-rate response tests, ten existing
water acoustic/thermal tests and three original molecular bulk-memory tests
passed. The new tests check original low/high-pressure density and endpoints,
complex-response and energy-flux coherence, independent native ramp histories
with refinement, time/frequency scaling and input refusals. Rustfmt passed.
This is a native source-module validation, not a whole-workspace Cargo/DSR pass.

### Exact finite-rate compression histories with density and heat

WaterCompressionHistory uses the same retained RelaxingRegionalWater and the
existing visco::PronyRamp to advance prescribed linear volume-compression
ramps, including frozen jumps, holds, reversals and recovery. Its Maxwell
memory follows qdot=edot-q/tau, pressure perturbation is K_eq*e+delta_K*q,
and independent arm heat is the integral of delta_K*q^2/tau. Every accepted
step reports work, stored-energy change, heat and the floating-point balance.

For fixed parcel mass, V/V0=1-e and rho=rho0/(1-e). A separately evaluated
density increment preserves small signals lost in subtraction of absolute
densities. The explicit compression bound is positive and less than one;
it guarantees positive parcel volume, while the caller remains responsible
for the stricter phase and linearization regime. This geometric conservation
identity is not a nonlinear IF97 constitutive claim. Coefficients, phase
inventory and the original thermal state remain fixed, and all energy
densities are perturbation quantities per original volume.

PreparedWaterCompressionRamp shares the original exponential coefficients
among histories with equal relaxation time. Accepted application needs no
allocation, exponentials or EOS evaluation. Invalid state, amplitude,
duration, time-ratio representation or incompatible preparation leaves
history unchanged. The hosting system owns thermal deposition and background
pressure work; the history introduces no spatial transport or split solver.

Four focused native Rust regressions passed against the exact source and
genuine existing PronyRamp: frozen jumps and analytic relaxation heat,
fixed-mass geometry, independent Simpson pressure-work/heat integration,
exact-ramp subdivision, periodic loss versus the harmonic response, atomic
refusals and retained sub-ulp density increments. Rustfmt passed. Combined
with the source response and prior water/molecular regressions, 21 focused
material tests passed. This is not a whole-workspace Cargo/DSR pass.

### Exponential diffusion on finite binary gas networks

`species::diffusion::network::binary::advance_candidate` advances a closed,
fixed-volume binary gas network through a complete requested duration. Each
symmetric sweep applies the original `diffuse_binary_exponential_with` owner
in canonical physical-edge order: forward half updates, one central full
update, then reverse half updates. Every invocation refreshes actual face
T/p, source-selected diffusivity, upstream species enthalpies and both original
caloric endpoints. Source, withdrawal, local T/p, mass, energy and pair-entropy
checks remain authoritative. Input order, reporting IDs and edge orientation
do not change the endpoint arithmetic.

Branches, cycles, disconnected components and unequal finite volumes are
supported. Each open edge must join exactly two original species and supply
their binary material. At most one edge connects an unordered endpoint pair;
distinct parallel source laws are not combined. Closed faces retain the
original source-free identity.

Caller controls select symmetric substeps and bound actual pair evaluations.
For E>0 links, each sweep costs 2*E-1 pair calls, including closed faces.
Integrated species and partial-enthalpy transfers and pair energy residuals
use compensated accumulation; detailed receipts retain only the last sweep.
Independent original-to-final species, total-mass, source-relative energy and
entropy audits and the energy-residual sum must fit fixed whole-window
budgets, without multiplying tolerances by the step count. Late failure or
cancellation returns no partial network.

A single isothermal/isobaric constant-D edge retains its analytic exponential
endpoint without an explicit diffusion CFL restriction. Symmetric composition
of exact constant-state pair flows is second order under refinement.
State-dependent or nonisothermal frozen-coefficient pair updates remain first
order; symmetry alone does not raise their order. Large stiff sweeps can retain
substantial network splitting error. This API requires substep refinement and
supplies no adaptive time-error estimate or unconditional nonlinear stability
claim. Bulk flow, pressure diffusion, reactions and moving volumes retain
their separate owners. G1/G3/G4 fixtures cover analytic stiff relaxation,
chain refinement and conservation, variable-source cycles, canonical ordering,
closed faces and transactional work/source/cancellation refusal.


### Affine stress-charge piezoelectric thermal sources

`ThermalLongitudinalPiezoelectric::from_stress_charge_tensor` admits complete affine C^E(T), e(T), epsilon^S(T) source tensors over a positive absolute-temperature interval with constant reference density. It preserves source selections, fixed frame/direction, and native transverse relaxation, differentiates the resulting non-affine longitudinal inverse exactly through second order, and supplies the existing rod/circuit thermal owner with consistent Helmholtz energy, entropy, internal energy, and fixed-coordinate field capacity. Complete source endpoint SPD establishes source-interval stability; a whole-interval inverse bound supplies conservative field-capacity admission. This is analytic model admission with roundoff padding, not an interval-arithmetic certificate. Synthetic G0/G1/G3 checks do not establish measured-material accuracy. No pyroelectric offset, thermal expansion, phase transition, depoling, or hysteresis is inferred.

## Regional fixed-volume liquid-vapor density and energy

`IsochoricWater::new_regional_coexistence` reuses the original fixed-specific-volume
internal-energy inverse and heat transaction with a retained
`If97FluidSaturationPair`, including Region 3. The domain is
273.16 <= T < 647.096 K. Each bounded temperature window must have either
maximum <= 623.15 K or minimum > 623.15 K; this entry does not stitch the
separately fitted Gibbs and Helmholtz source equations.

Every evaluated original pair must bracket the fixed specific volume and lie
inside the declared pressure bounds. Exact phase endpoints retain their pair;
departure from coexistence refuses. The original dense pressure owner retains
its controls, work and cancellation. `rebind_specific_volume` preserves the
source mode and dense controls when a finite inventory changes V/M.

The analytic fixed-volume heat capacity differentiates each original phase's
u using that phase's actual EOS pressure. Declared saturation pressure and the
two root-pressure residuals remain separate. Raw h, u and s are unchanged:
`pressure_closure_defect_j_kg` and the heat report's endpoint joule defects
expose h-u-p_pair*v without crediting it as heat or repairing a source value.
Density is 1/v_source; the requested fixed volume and recovery residual remain
independently visible. Existing original heat/enthalpy balance budgets still
apply to every accepted step.

No single-phase crossover, finite-rate phase conversion, critical-state
extrapolation or global heat-capacity enclosure is claimed. The existing
finite-cell owner remains responsible for conserved mass, energy, stream
enthalpy and reciprocal boundary work.

Five added native tests check legacy Gibbs bit parity, original Region 3 pair
retention and density, an independent raw-u finite-difference Cv oracle at
630/640/646 K, energy inversion and finite-inventory heat/cool balances,
source-preserving rebinding, exact endpoints, bounded-domain refusal and
cancellation. The exact source module passed all 11 native tests: these five
plus six unchanged original isochoric regressions. This is a focused native
module check against genuine dependencies, not a full workspace Cargo/DSR gate.

## Adaptive exponential diffusion across finite binary gas networks

`thermochemical_gas::energy::parcel::transfer::species::diffusion::network::binary::advance_adaptive_candidate`
completes a requested duration through the existing molecular temporal controller.
Each comparison uses one full exponential sweep and two half sweeps; only the
accepted fine path contributes state and integrated transport. Error scales
compare actual species inventories and endpoint temperatures with conservative
first-order step-doubling scaling. Every attempted edge refreshes its source
state and upstream partial enthalpy. The constant-state exact-flow special case
retains its original order; no higher nonlinear order is inferred from symmetry.

Separate complete-trial and actual pair-call budgets include rejected and
partially failed sweeps. Only original-owner withdrawal/finite-step entropy
refusals and typed caloric bracket/convergence failures permit cutback.
Source/domain and malformed-input failures remain direct. Original whole-window
species, mass, source-relative energy and entropy budgets apply unchanged;
late errors and cancellation publish no partial network. Gross transfers include
only accepted fine intervals, and detailed receipts retain only the last one.

This extends the declared binary molar-frame closure. It adds no multicomponent,
pressure-frame, Soret/Dufour or fully coupled implicit constitutive authority.
G1/G4 regression definitions cover analytic unequal-volume relaxation, network
mode refinement, temperature-dependent transport against an independent ODE,
partial-sweep work charging, source failures and atomic work/cancellation gates.


## Directional metal yielding and permanent strain (Hill48)

`plastic::hill::Hill48` implements small-strain associated Hill48 plasticity
with isotropic elasticity and explicitly energetic linear isotropic hardening.
The six authored dimensionless coefficients define
`q² = F(syy-szz)² + G(szz-sxx)² + H(sxx-syy)² + 2L syz² + 2M szx² + 2N sxy²`.
The quadratic must be positive definite on all five deviatoric directions;
valid negative individual F/G/H are admitted. Hydrostatic pressure does not
yield this dense-metal law, and retained plastic strain is isochoric.

The distinct `HillState::alpha` is accumulated work-conjugate flow, not a
silently reused von Mises equivalent strain. Five orthonormal modes are
prepared at construction. The backward-Euler update uses fixed-size arrays,
one bounded monotone scalar solve, and its analytical consistent tangent.
No numerical derivative, per-point matrix eigensolve or heap allocation enters
the checked return. Failed updates leave accepted history untouched; the
infallible material trait returns nonfinite stresses/tangents on refusal.

`stored_energy` inspects elastic and hardening storage without evolving state.
`advance_with_energy` separates yield-offset physical dissipation from stored
hardening and backward-Euler numerical loss. Neither of the latter is heat.
`oriented::OrientedHill48` precomputes proper tensor transformations for fixed
material axes, returns spatial stress/tangents and retains local-axis history.
Identity orientation preserves the native arithmetic. It rejects reflections
and nonorthogonal frames rather than normalizing them.

`source::ResolvedHill48` binds the complete original material bundle and its
ten scalar requirements at one fixed physical point. Reference yield radius,
Hill work-conjugate hardening and all six directional coefficients have
explicit property names; generic tensile yield, J2 hardening, hardness and
elongation do not fill missing data. Original claims, conditions, quantities,
uncertainty, extra physics receipts and explicit orientation survive. Evolving
strain/rate/history-dependent coefficients cannot be frozen into this model.

G0/G1/G3 regressions in `plastic_hill`, `plastic_hill_orientation` and
`plastic_hill_source` compare J2 limits, independent directional yield/flow,
nonproportional tangent finite differences, pressure and coordinate changes,
work/heat balance and original source admission. Native 2-D/3-D consumers are
covered by `fs-solid/tests/plastic_hill.rs`. These are numerical constitutive
checks, not measured-metal calibration, finite-strain anisotropic forming,
texture evolution, necking, fracture, fatigue or forming-limit predictions.

## Homogeneous multicomponent liquid sources

`liquid::solution::multicomponent::MulticomponentLiquidSolutionLaw` binds the
original N-component regular-solution owner to one pinned liquid-property
bundle. Explicit species bindings carry molar masses and N or N-1 distinct typed
mole-fraction axes; any dependent species is explicitly declared. Sorting keeps
each species, mass and source coordinate together. Inventory inputs and all
returned arrays use the thermodynamic owner's canonical species order.

`homogeneous` retains the actual finite component amounts, evaluates the same
homogeneous free energy and refreshes the original density, viscosity, heat
capacity and conductivity at that composition and temperature. Pressure and
other source coordinates stay fixed. Mass, mass fractions, source-derived
volume and molar concentration feed existing fluid/thermal transport consumers.
The retained thermodynamic state supplies cross-component Maxwell-Stefan forces
and simplex curvature; negative curvature does not trigger an implicit flash.
Positive traces lost in normalization or inaccurately represented by the
dependent source coordinate refuse. An absent dependent component requires
exactly zero chart remainder; a positive component requires relative agreement
within 64 machine epsilons. N explicit coordinates avoid subtractive loss when
no species remains abundant across all physical phases. Every explicit source
coordinate is refreshed from the retained composition, including positive
traces. Source failures and cancellation publish no partial cell.

The caller qualifies pair parameters, ordered chemistry and mixture-property
data. The homogeneous evaluation performs no phase search and introduces no
pure-fluid averaging, inferred expansion, spatial evolution or
excess-enthalpy/Cp consistency claim. G0/G1/
G3/G4 fixtures in `tests/liquid_multicomponent.rs` exercise authored ternary
surfaces, independent component/mass/volume accounts, exact equal-diffusivity
fluxes, unstable homogeneous compositions, source/species permutation, traces,
typed-coordinate refusal and atomic source/cancellation failures.

## Source-bound multicomponent liquid calorics

`liquid::solution::multicomponent::caloric::MulticomponentLiquidCaloricLaw`
binds the original homogeneous mixture law to one explicitly qualified
pure-liquid `IsobaricHeatCapacityLaw` per species. Canonicalization keeps
species, source receipts and original molar masses together. Every original
Cp knot survives the mass-to-molar conversion. The original typed mixture
pressure must match every pure source, and the numerical owner admits only
their common datum and supported temperature interval.

Preparation, unequal-temperature mixing and signed heating/cooling delegate
all caloric work to `MulticomponentLiquidCaloricModelV1`. The same nonideal
potential supplies mixing heat and entropy once; total-mixture transport Cp
is never substituted for pure caloric sources. Only the accepted thermal
endpoint refreshes actual density, viscosity, conductivity and transport Cp
through the original mixture source binding. The returned state retains the
complete caloric inventory, actual resolved liquid cell and pure source laws.
Source-relative heat and entropy increments also require the same original
component inventory and complete source owner.

Source/domain refusal or cancellation returns no partially published endpoint.
The mixture sources must support the accepted temperature/composition even
when the pure Cp sources admit a larger interval. An isolated zero-heat mix
inherits the numerical owner's entropy-decrease refusal; with external heat,
system entropy needs a separate reservoir exchange before an entropy-production
claim. The caller qualifies pure-liquid source roles and interaction data.
No equilibrium phase split, pressure evolution or consistency certification
between independently authored transport Cp and caloric Cp is introduced.
Focused G1/G3 integration tests cover ternary nonideal mixing, original receipt
retention, actual material refresh, variable source Cp heating/cooling, pressure
and datum mismatch, foreign source provenance and atomic final-source refusal.

### Homogeneous multicomponent thermal storage, conduction and streams

`homogeneous_thermal_storage` binds those original calorics to `ThermalLaw`:
mass-normalized h, original safeguarded inversion, exact Cp/T entropy and a
conservative source-knot CFL bound. Closed storage uses a cold-endpoint datum.
`MulticomponentLiquidConductionLaw` pins the original mixture conductivity
role on the actual composition/pressure path, retaining every source knot and
receipt through the existing scalar conductivity owner. One caloric inverse
supplies exact dT/dh=1/Cp and dk/dh=(dk/dT)/Cp to the native implicit solvers.

`homogeneous_stream_storage_from_mass_fractions` additionally retains the
common pure-species enthalpy datum and original species-reduced entropy.
Different compositions and thermal brackets share compatible source references;
the native stationary-species, advection and entropy owners retain transport
and accepted-stage accounting. Original nonideal mixing heat and entropy enter
once. Authored mass fractions, positive traces and complete source identities
are preserved. The homogeneous constraint, fixed pressure and reference
geometry remain explicit; no phase conductivity, pressure evolution or
transient composition storage is inferred. Numeric caloric trials do not
refresh unrelated mixture transport properties.

`tests/liquid_multicomponent/caloric/storage` checks nonlinear Cp, reference and
inventory scaling, source conductivity/tangents, native two-cell implicit heat
exchange and ternary BE/SDIRK mixing against independent energy/entropy
expressions, with original-source and cancellation refusals.

## Source-resolved multicomponent equilibrium phases

`liquid::solution::multicomponent::equilibrium` resolves every positive phase
of an accepted `MulticomponentEquilibriumCandidateV1` through the original
`MulticomponentLiquidSolutionLaw`. Complete model equality, common temperature,
original feed and phase inventories, accepted phase compositions and species/
mass conservation are checked with relative floating-point bounds and exact
absences. Each phase retains its own source-resolved density, viscosity, heat
capacity, conductivity, mass, volume and component inventory. Only extensive
mass and actual phase volumes are aggregated; effective transport properties
and morphology are not inferred.

The original seeded local candidate and its numerical report are retained.
Material resolution performs no additional phase search and provides no global
Gibbs-minimum or phase-count certificate. Source support is required at the
actual phases rather than a hypothetical homogeneous feed. Any phase source,
domain, trace-coordinate or cancellation failure rejects the whole unpublished
result. Complete explicit composition charts permit trace-sensitive source
evaluation across strongly separated phases without subtractive reconstruction.

Focused integration tests check distinct phase properties, finite-feed scaling,
one-phase parity and closed supports, original-owner and source refusals,
dependent-coordinate trace refusal, an all-explicit trace-sensitive source,
and cancellation at every observed polling boundary.

### Source-bound multicomponent equilibrium heat capacity

`liquid::solution::multicomponent::caloric::MulticomponentLiquidCaloricLaw`
binds that original liquid owner to one explicitly qualified pure-liquid
`IsobaricHeatCapacityLaw` per species. Named inputs are canonicalized together;
each original molar mass converts its own complete Cp knot sequence to molar
units without resampling. Original source laws and receipts remain retained.
The common datum, positive static pressure and intersected original temperature
support must agree. A species label alone does not establish pure-liquid
source qualification.

`isobaric_caloric_response` accepts only a resolved equilibrium belonging to
the SAME complete material and thermodynamic owner, including original source
bindings, component molar masses and coordinate roles. It borrows the exact
accepted phases and delegates once to the original thermochemical candidate's
analytic response. Total equilibrium Cp divided by the original finite mass
is exposed in J/(kg K), alongside frozen/redistribution capacities and
species-conserving phase derivatives. A singular branch, invalid source,
exhausted response budget or cancellation publishes no result.

This is a local differential response qualified by the retained seeded
candidate and its numerical report. Phase transport Cp is not replaced or
declared consistent with the independently supplied free-energy law. No
effective transport property, global-minimum certificate, phase-onset
continuation or enthalpy inversion is inferred. Native integration cases
check complete source retention and ordered mass conversion, an independent
three-phase enthalpy derivative, extensive scaling, complete foreign-owner
refusal and cancellation.

### Source-bound finite multicomponent equilibrium heating

The existing `MulticomponentLiquidCaloricLaw::equilibrium_caloric_state` binds
an already resolved equilibrium to its complete original material and caloric
owner. The original phases are retained without another source query or phase
search. `add_heat_equilibrium` delegates signed heat and bounded fixed-phase
continuation to that SAME thermochemical owner, including temperature-dependent
pure Cp and thermal redistribution of phase inventories.

Only the accepted final candidate refreshes each phase's original density,
viscosity, conductivity, mass and volume. The returned caloric/material state
retains that exact candidate, source receipts, source-normalized specific heat
and finite heat/entropy balance. A foreign owner, unsupported final phase,
exhausted continuation budget or cancellation returns no partial endpoint and
leaves the original state unchanged. Numeric thermal trials perform no
material-source queries. Phase transport Cp remains its independently selected
property rather than an inferred effective equilibrium transport coefficient.

This operation preserves the accepted phase count and support. It supplies
seeded local numerical continuation with explicit branch/domain refusal, not
automatic phase appearance/disappearance or global-equilibrium certification.
Native integration tests heat and cool three phases across an original Cp
knot against independent scalar common-tangent and analytic H/S integrals,
check actual composition-dependent material refresh and conserved mass, and
reject a final phase leaving its original source chart, foreign caloric
references and last-poll cancellation atomically.

### Spatial heat transport with multicomponent liquid phase redistribution

`MulticomponentLiquidCaloricLaw::equilibrium_thermal_storage` prepares one
resolved local equilibrium branch on the existing native `ThermalLaw` owner.
Total composition, pressure and reference geometry stay fixed while actual
phase amounts and compositions respond to heat. The original pure Cp/CpT
integrals, phase mixing storage, continuation and safeguarded enthalpy inverse
retain the full redistribution heat capacity and entropy. No sampled chart or
independent apparent Cp replaces those sources. Closed storage uses the cold
endpoint as its enthalpy datum; numerical heat increments cancel original
reference offsets before arithmetic. The original downward-rounded pure-Cp
minimum bounds dT/dh on the admitted stable local branch.

Both endpoint phase bundles resolve through their original material sources.
Preparation uses two separately bounded continuations and later queries retain
the original per-query phase-solve limit. The complete temperature interval
must admit the same phase count and support; disappearing or singular phases
refuse. This preserves the native continuous storage contract and does not
promote finite phase discovery into a global Gibbs or branch-identity proof.
Phase appearance and genuine first-order latent intervals need their own
coexistence relation; they are not smoothed into this branch.

`MulticomponentLiquidEquilibriumConductionLaw` binds this exact storage to an
explicitly declared bulk `ConductivityCurve`. Individual phase conductivities
do not determine morphology or interphase resistance and are never averaged
implicitly. One caloric inverse supplies T and dT/dh=mass/Cp; dk/dh uses the
same temperature and original conductivity derivative. All original knots and
domain refusals remain. Native G1 integration compares a two-cell nonlinear
heat exchange against independent common-tangent, enthalpy and entropy
calculations while phase inventories redistribute across Cp/conductivity knots;
it also checks the analytic tangent, conservation and cancellation. No
phase-selective mass flow or separate interphase temperature is inferred.

### Stationary flow into redistributing liquid phases

`equilibrium_stream_storage_from_mass_fractions` binds an actual resolved
equilibrium inventory to the existing stationary-species, implicit advection
and accepted-stage entropy owners. Exact-unit-sum canonical mass fractions
must match the original species inventory to relative arithmetic roundoff;
positive traces and exact absences remain intact. Homogeneous and separated
streams share the same complete pure-source enthalpy and reduced-entropy
reference, independent of composition, prepared bracket and phase topology.
A separate local caloric identity retains the actual cold phase allocations.

Each receiving cell uses its original equilibrium enthalpy inverse and full
redistribution Cp. Actual incoming mixing enthalpy and each retained phase's
mixing entropy enter once through the existing flow balance. Equal-temperature
inlets can therefore change the receiver temperature when they mix. There is
no extra heat-of-mixing source, averaged phase conductivity or homogeneous
substitute for the separated receiver. Closed storage keeps its zero datum;
only stream storage restores the common source datum, and unrepresentable
interior enthalpies refuse.

The prepared local phase count/support, stationary overall composition, fixed
pressure and reference geometry remain explicit. Phase-selective flow,
transient species storage, changing morphology and finite-rate interphase
heat transfer are outside this owner. Original continuation and source bounds
apply to every caloric query; cancellation publishes no partial stream law.

Two added native integration tests compare actual phase entropy and native
backward-Euler/SDIRK mixing against independent ternary common-tangent H/S
calculations, including source-Cp-knot crossing and equal-temperature heat of
mixing. They check common references across brackets and compositions, source
and inventory refusal, accepted-stage entropy balances and bit-identical
thermal evolution with the optional entropy audit enabled or disabled. The
whole equilibrium and homogeneous integration roots pass 28 and 15 tests
respectively against rebuilt original material/thermochemical sources and
real native dependencies. These are focused native checks, not a full
workspace build or measured-material calibration.

### Native liquid heat transport through a genuine latent interval

`MulticomponentLiquidEquilibriumThermalStorage::join_coexistence` joins adjacent
same-source, same-feed original branches only through the L1 common-tangent
owner. `MulticomponentLiquidCoexistenceThermalStorage` supplies one continuous
nondecreasing T(h) path with both sensible branches and the actual isothermal
latent interval. Its h=0 datum is the lower-enthalpy coexistence endpoint,
keeping a small latent heat out of an unrelated cold-reference offset.
Original branch Cp, continuation limits and direct H/S increments remain;
the latent interval has exactly zero dT/dh, including the selected one-sided
corner tangent. A temperature-only enthalpy or Cp query at coexistence refuses
the set-valued relation instead of choosing an unrequested phase amount.

The native `ThermalLaw` entropy is its defining integral dh/T. The separate
`specific_source_entropy_change_j_kg_k` query retains actual original phase
entropy; the two differ only by the admitted portion of the visible original
latent closure residual and arithmetic normalization. No source value is
changed to satisfy the thermal account. Both sensible branches use the same
conserved normalization mass and retain original source/domain refusals.

`MulticomponentLiquidCoexistenceConductionLaw` binds this path to an explicitly
authored bulk k(T) curve. One original inverse supplies T, dT/dh and dk/dh;
both derivatives are exactly zero on the plateau. Native adaptive observables
retain the upper-ensemble mass fraction while temperature is constant. Both
ensembles are liquid: this is the named higher-enthalpy ensemble fraction,
not vapor quality or a solid/liquid fraction. The bulk conductivity closure
does not infer phase morphology, interphase resistance or effective transport
from individual source conductivities.

Independent ternary common-tangent and original Cp/CpT integrals exercise
source-preserving sensible/latent inversion, actual three/four/one-phase
inventories, source entropy closure and native two-cell implicit heat transfer.
The join retains finite offered branch/common-tangent qualification; it does
not turn bounded phase discovery into a global phase-selection certificate.

### Source-bound conserved phase relaxation and mixing

`MulticomponentLiquidCaloricLaw::relax_homogeneous` and `mix_equilibrium` extend
the SAME source-bound phase-caloric owner to an actual homogeneous parcel and
to pairs of homogeneous or separated parcels at different temperatures. The
explicit admitted branch anchor must retain the original material, pure Cp
and pressure binding and exactly the combined physical species inventory.
It is numerical phase evidence; its stored energy is not another physical input.
Original pure Cp/CpT integrals and all actual phase mixing storage determine
the energy-conserving endpoint, including the temperature change from zero-heat
demixing and the total isolated-entropy check.

The existing bounded phase continuation and caloric root perform every thermal
trial. The original final-state resolver then queries each actual phase's
composition and temperature, retaining the exact accepted candidate, source
receipts, mass, density, volume and transport properties. Only complete source
resolution publishes the balance. A final phase outside its source chart,
foreign input binding, branch/inventory mismatch or cancellation leaves all
input states unchanged. No averaged phase transport, isothermal heat sink,
extra latent storage, global phase selection or kinetic interpretation is added.
Integration regressions cover adiabatic demixing across an original Cp knot,
unequal-temperature phase/homogeneous mixing against independent H/S oracles,
actual per-phase properties and volumes, final-source refusal and cancellation.

### Source-bound thermal endpoints with liquid phase discovery

`MulticomponentLiquidCaloricLaw::add_heat_discovered` and `mix_discovered`
delegate to the original bounded thermochemical discovery/caloric owner, so
accepted thermal endpoints may contain a different number of liquid phases.
Every input retains its complete original material, pure-Cp and fixed-pressure
binding. All numerical temperature samples use the same declared discovery
policy and actual phase H/S; only the accepted endpoint refreshes material
sources. The original phase resolver queries each actual final composition
and temperature and retains its receipts, density, viscosity, conductivity,
mass and volume. No averaged phase transport law is introduced.

The material balance retains the exact numerical candidate, original final
discovery report, bounded work and thermal/entropy accounts. Source refusal
for a newly appearing phase, foreign input binding, exhausted numerical work
or last-poll cancellation publishes no partial result. Finite discovery remains
qualified local evidence, not global stability, a continuous phase trajectory,
a latent-heat plateau or nucleation kinetics. Focused integration regressions
cover three-to-one heating and one-to-three cooling against independent H/S
integrals, unequal-temperature mixing of actual separated/homogeneous inputs,
source refresh and conservation, newly discovered phases outside the original
material chart, foreign references and final publication cancellation.

### Full three-dimensional piezoelectric thermal storage

`piezoelectric::thermal::ThermalTensorPiezoelectric` binds complete affine original stress-charge stiffness, signed coupling and clamped-permittivity source tensors over a positive absolute-temperature interval with constant reference density. It preserves full source selections, nonthermal context and the fixed frame transformation. At fixed full Mandel strain and Cartesian electric displacement, the native tensor/factor owner supplies Helmholtz energy, entropy, internal energy and field heat capacity, without scalar projection or fresh response-time factorization. Complete endpoint tensor admission establishes affine-source stability. A rotation-invariant bound supplies a conservative whole-interval capacity floor under explicit strain/displacement envelopes; consuming caloric owners add their actual capacity and enforce those bounds. Reversible thermal exchange remains separate from material loss.

`with_thermal_bias` optionally binds six explicit Mandel thermal-expansion
coefficients and three explicit **clamped** pyroelectric coefficients from the
same original card, query and source frame. Every component, including zeros,
must be sourced and constant throughout the admitted temperature interval.
The caller declares the in-range reference temperature at which thermal strain
and polarization vanish. Source refreshes retain the selections and reference;
second-rank thermal strain and vector polarization co-rotate with the tensor.
For `A=alpha*(T-Tref)` and `P=p*(T-Tref)`, the shared potential is
`psi/V=(S-A):C:(S-A)/2+(D-e*S-P).epsilonS^-1.(D-e*S-P)/2`.
The constitutive stress, field and displacement, thermal entropy and curvature
all derive from that potential. Coupling acts on total strain: the free shorted
pyroelectric response includes `p+e*alpha`, so a measured total coefficient must
be converted before entering a clamped source role.

The caloric chart is the stress-free, zero-field reference inventory;
an unconverted clamped heat capacity is not interchangeable with it. With
thermal expansion and temperature-dependent stiffness, field heat capacity
can have either sign. Directional endpoint bounds include expansion curvature,
stiffness-slope cross terms and pyroelectric inverse-dielectric response, while
the consuming thermal owner enforces positive combined caloric/field capacity.
No spontaneous polarization, hysteresis, phase transition, depoling or
measured-material validation is inferred. Synthetic regression definitions
cover full potentials and derivatives, source/refusal boundaries and frame
rotation; native execution remains a separate verification obligation.

## Directional tensile strength in orthotropic material axes

`strength::OrthotropicTensileLimit` resolves three distinct, caller-named
tensile-strength claims in Pa from one original material card and complete
physical query. Each source must be positive, finite and static at the retained
external conditions. Missing transverse strength, unsupported rate/frequency
dependence and foreign-card refreshes refuse; no longitudinal-to-transverse
inference, safety factor or implicit strength selection is supplied. Refreshing
temperature or other supported conditions preserves the original claim choices.

`assess_material_stress` evaluates the tensile branches of the orthotropic
maximum-stress criterion, `max_i(max(0, sigma_ii) / X_ti)`, in fixed material
axes. These axes are not replaced by the principal directions of each load.
The fixed-size calculation allocates nothing, retains all six signed
tensor-shear stresses, reports each directional utilization and resolves exact
ties by the first material axis. The spatial consumer owns the actual frame
transform. Nonfinite stresses and unrepresentable utilization refuse.

This is a tensile-limit observation, not a complete failure surface:
compression, material-axis shear failure, stress interaction, damage evolution,
fracture energy and load redistribution are unassessed. Zero tensile demand
under pure material-axis shear does not imply safety. Engineering UTS is not
converted to a finite-strain Cauchy failure law. The three tensile branches
follow the maximum-stress definition in the Abaqus documentation, “Generalized
Damage Initiation Criteria for Anisotropic Materials”:
<https://docs.software.vt.edu/abaqusv2025/English/SIMACAEMATRefMap/simamat-c-damageinitgeneral.htm>.

G0/G3 regressions in `tests/orthotropic_tensile.rs` cover each directional
threshold, transverse versus longitudinal loading, mixed signed loading,
retained unassessed shear, deterministic ties, stress-unit scaling, original
source/temperature refresh, missing data, units and arithmetic refusal. All
coefficients are manufactured; no measured-material validation is asserted.

### Oriented magnetic permeability and reluctivity

`field_storage::spatial::SpatialMagneticLaw` admits either the original isotropic
absolute permeability or all three `magnetic_permeability_1/2/3` source values in
H/m with an explicit fixed orthonormal material frame. Every principal source
retains the shared fixed-condition admission; relative or porous permeability,
missing components, and field/frequency/strain/stress-dependent sources refuse.
The law supplies the complete symmetric positive tensor, including off-diagonal
flux, and retains its principal values and axes. `SpatialPermeability` applies
`B=mu H`, resolves `H=mu^-1 B` through those original principal directions, and
prepares inverse weights by rotating reciprocal principal values instead of
inverting a rounded matrix. Rounded tensors that lose positive definiteness
refuse; field assembly retains its own interval-SPD check.

The scalar path keeps its original multiply/divide arithmetic. Coefficients are
prepared once, and evaluation allocates nothing. The fixed-frame
`reluctivity_change_action` computes `(nu_next-nu_old)B` directly from principal
changes and returns exactly zero for held values. It does not itself authorize
source refresh or account for thermal entropy; those remain consumer duties.
No anisotropic saturation, history, gyrotropy, frequency loss, mechanical stress
or measured-source qualification is inferred. Focused G0/G3 regressions compare
full fields, inverse response, energy and parameter changes with independent
principal-frame calculations, source substitutions and source/domain refusals.

### Source-bound isotropic vector-speed Forchheimer coefficients

`PermeabilityTemperatureLaw::with_isotropic_forchheimer` binds one exact scalar
beta key from the original resolved material bundle to `f = -rho beta |u| u`,
with Euclidean superficial speed and no factor one half. Beta is finite,
nonnegative and dimensionally inverse metres. This explicit convention is
independent of the intrinsic-permeability frame and is not inferred from equal
separable coefficients or isotropic permeability. One permeability law admits
only one inertia convention.

`ResolvedPermeability::isotropic_forchheimer_beta_per_m` identifies the scalar
vector law; the original directional accessor continues to identify separable
inertia. The retained positive finite speed bound applies to the full Euclidean
norm for the scalar law and to each principal component for the separable law.
Temperature continuation preserves that convention, all selected source claims
and all other original material conditions. Zero beta remains exactly zero.

Consumers enforce the declared speed norm on the same velocities used by their
force and dissipation quadrature. The native MAC coupling uses its original
eight face tuples for force, physical heat and support checks. No pore-speed
conversion, second thermal store, moving skeleton, saturation, turbulence or
compressibility closure is inferred. Focused synthetic regressions cover source
continuation, distinct conventions, unsupported calibration and source-switching
refusal; they do not establish experimental validity.


## Hill48 stress and history evaluation

The local `Hill48` and `OrientedHill48` implementations share the same checked
constitutive return across full responses, stress-only residual evaluations,
and accepted-history updates. Stress and history requests do not construct a
6x6 tangent; oriented residual evaluations rotate only stress, while the
accepted history remains in the original material coordinates. Full responses
retain the analytic discrete tangent and tensor-shear transformation.

Public effective-stress evaluation first preserves finite deviatoric differences
in the original stress scale. It normalizes finite stress components only when
the direct evaluation is nonfinite, recovering representable answers despite
intermediate overflow without discarding small represented differences on a
large hydrostatic offset. It still refuses nonfinite inputs or an unrepresentable
effective stress. Independent numerical oracles and exact path-parity cases live
in `tests/plastic_hill_robustness.rs`; near-hydrostatic scale regressions are in
`tests/plastic_hill.rs`.


## Retained layered plane-stress plastic sections

`plastic::section::LayeredPlaneStressSection` integrates the original native
plane-stress reduction through an ordered stack in one declared material frame.
Every thickness point retains the complete three-dimensional material history
and transverse strains. Generalized engineering strains
`[ex,ey,gamma_xy,kx,ky,kxy]` produce work-conjugate membrane forces and bending
moments with the full consistent 6-by-6 tangent. Asymmetric layers, a shifted
reference surface, partial yielding and native nonsymmetric tangents retain
their coupling; the original constitutive owner performs every return map.

Geometry and two- or three-point Gauss rules are prepared once, with explicit
layer subdivisions and allocation bounds. Evaluations can reuse the original
plane-stress factorization workspace. Each trial borrows the accepted history,
and timed laws advance their physical interval only in the returned accepted
candidate. In-plane face strains are bounded in addition to sampled full
tensor strains. Elastic A/B/D integration is exact; moving yield fronts require
thickness refinement and no quadrature-error bound is inferred.

`plastic::section::energy` audits one exact accepted native increment for J2,
Voce, linear/power/Voce overstress and Chaboche laws. It integrates the original
`advance_with_energy` account, retaining elastic, isotropic-hardening and
kinematic storage separately from physical dissipation and numerical loss.
Chaboche recovery is counted once. Valid kinetic overstress is preserved;
rate-independent incoming histories must satisfy their native yield allowance.
Complete returned-history correspondence is required. Signed transverse
constraint-residual work closes generalized membrane/bending work without
being converted into another heat source.

This remains small-strain plane stress in fixed material axes. It supplies no
transverse-shear law, finite-rotation forming, layer separation, heat transport,
temperature evolution or measured-material calibration. Constitutive energy
interpretations and coefficients must remain unchanged over an audited step.
Six section regressions and three energy regressions define elastic/resultant
work identities, reference shifts, partial yielding and springback,
nonproportional Chaboche tangents, timed holds, domain/history refusals, native
heat and numerical-loss separation. These are added test definitions; native
compilation and execution remain outstanding in the resource-exhausted runner.


## Hill48 calibration from planar plastic strain ratios

`Hill48Coefficients::from_planar_r_values` converts explicit Lankford plastic
width/thickness ratios at 0, 45 and 90 degrees to the rolling direction into
the existing associated Hill surface. The rolling-direction uniaxial yield
is the reference radius: `G+H=1`, `H/G=r0`, `H/F=r90`, and
`N=(r45+1/2)(F+G)`. Independently supplied transverse shear yield ratios in
planes 23 and 31, normalized by `Y_rolling/sqrt(3)`, determine L and M.
The three planar tests alone do not identify those shear coefficients.

`source::planar::resolve_hill48_from_planar_r_values` and
`ResolvedHill48::try_from_planar_r_material` retain the original material
bundle, fixed query, selected source claims, uncertainty and usage receipts.
Their nine scalar requirements include the ordinary elastic coefficients,
explicit rolling-reference yield and Hill hardening modulus,
`hill_lankford_r_0`, `hill_lankford_r_45`, `hill_lankford_r_90`,
`hill_shear_yield_ratio_23` and `hill_shear_yield_ratio_31`. Derived F/G/H/L/M/N
are not passed off as original observations. Extra physics receipts survive.
Missing shear data, wrong quantities, evolving mechanical conditions,
nonpositive ratios, and numerically unrepresentable calibrations refuse.

This calibration fits plastic flow ratios within the fixed Hill48 surface;
it does not independently fit every directional yield observation. The input
ratios are plastic strain increments, not engineering elongation, total
strain or a failure measure. The original small-strain, small-rotation,
tension/compression-symmetric applicability and explicit material-frame
requirements remain unchanged. `tests/plastic_hill_calibration.rs` checks
independent 0/45/90-degree tension and compression endpoints, rolling/shear
yield normalization and explicit refusals; `tests/plastic_hill_planar_source.rs`
checks original receipts, fixed conditions, orientation and missing data.


`NativeSectionEnergyLaw::physical_duration_s` exposes each timed native
return's original interval; rate-independent laws return no interval. Linear,
power and Voce overstress and borrowed-law adapters retain their actual native
duration. The audited dynamic shell consumer refuses a different inertia
interval rather than advancing plastic history on an inconsistent physical
clock. This metadata performs no material retiming or temperature evolution.


## Source-bound thermoplastic plane-stress sections

`plastic::section::thermal::ThermoplasticSectionStep` freezes each original physical layer's source, old/new caloric enthalpy and actual duration, then resolves its native law once per layer. The original plane-stress reduction supplies total tensor strains, complete J2/Voce/Perzyna/Chaboche histories and consistent tangents. The exact native thermal increment is audited at every retained point; coefficient-dependent free energy and the change in mechanical entropy times temperature remain separate from physical/reversible heat and numerical loss. Original thickness weights produce per-layer and total reference-area accounts. Generalized work plus signed transverse constraint-residual work must close the first law, and full returned histories must match the accepted section. Consumers compare the reported resultants against their equilibrium endpoint.

Temperature is uniform within each physical layer; original geometry/density remain fixed and the source family's fully-solid and small-strain bounds apply. Thermal transport and atomic coupled publication stay with their existing owners. No melting, finite forming, evolving thickness or layer separation is inferred. Four focused synthetic regressions cover analytical thermal moments and reversal, yielded thermal histories with coefficient/entropy accounting, timed multi-branch relaxation and invalid bindings. These checks do not establish experimental material calibration or quadrature/time accuracy.


### Vector dielectric switching and temperature source continuity

The default-off `experimental-dielectric-hysteresis` feature exposes
`dielectric::hysteresis::vector::{VectorDielectricHysteresis,
ResolvedVectorDielectricHysteresis}`. Explicit isotropic quadratic branch fits
reuse the native coupled Euclidean play/spherical saturation solver and shared
immutable parameters. Electric-field and displacement-controlled candidates
retain three-dimensional polarization, consistent `dD/dE`, its native SPD
inverse, physical SI residuals, bounded work and cancellation. Original scalar
source receipts survive; they do not become rotating-field calibration evidence.

Physical switching heat is `sum(Ec |delta P|)`. Endpoint electrical work closes
with stored-energy change, this heat and a separately retained nonnegative
numerical remainder. Frozen-history `Phi(E;P0) = E.D - W - Q` has gradient `D`
and Hessian `dD/dE`; its displacement-controlled dual is `W + Q`. These are
incremental Newton objectives, not extra stored energy. Allocation-free,
factored pairwise differences retain the SAME accepted `P0`, preserving small
trial changes beneath large unrelated field energy. Active-set corners retain
the native selected tangent; turns and reversals require path refinement.
The fixed isotropic fit supplies no anisotropy, domain kinetics, breakdown,
strain or temperature-dependent switching law.

The scalar source adapter admits an owned absolute-temperature interval only
when every original constant fitted coefficient explicitly covers that complete
interval. `UniformHystereticCapacitor::at_source_query` rechecks the same selected
fit and fixed geometry without creating or resetting polarization history.
Temperature may change source-resolved transport and the resulting trajectory;
it does not change these fitted coefficients or imply electrocaloric heating.

`tests/dielectric_vector_hysteresis.rs` exercises scalar and linear limits,
coupled spherical stationarity, rotation covariance, inverse and variational
tangents, tiny trial changes about a large bias, cycling/refinement and source
continuation/refusal. `tests/dielectric_hysteresis.rs` also exercises explicit
temperature coverage and history-preserving source refresh.

The shared magnetic/electric vector saturation return refines its original
safeguarded multiplier until the normalized polarization also satisfies the
original full stationarity allowance. Radius convergence alone is insufficient
when radial normalization changes the direction of a small switching increment.
The same iteration budget, pinning norm, spherical bound, constitutive residual
and physical work gates remain active. An independent tangential-equilibrium
regression covers a finite rotating slip that previously passed the radial
criterion and then failed the final stationarity check.

## Three-dimensional orthotropic strength interaction

`strength::OrthotropicTsaiWu` admits the original nine positive strength
magnitudes `[X1t,X1c,X2t,X2c,X3t,X3c,S12,S23,S31]` in Pa and three explicitly
selected normalized interactions `[f12,f23,f31]`. All twelve named claims
retain the same card, complete typed conditions and selection receipts. Missing
compression, transverse, shear or biaxial data never receives a default. Source
refresh preserves these exact choices and rechecks the complete surface.

The full normal interaction matrix must be positive definite. Three separate
pair bounds do not establish that condition; open or singular surfaces refuse.
The normalized determinant must exceed `64*epsilon` times the sum of the
absolute determinant terms; positivity unresolved at roundoff also refuses,
without clipping source coefficients or regularizing the strength surface.
The fixed material-frame tensor-shear stress drives the original 3-D Tsai--Wu
polynomial. A cached Cholesky factor, scaled norm and stable quadratic root
return a nonnegative homogeneous utilization without squared dimensional
stresses, eigensolves, allocation or source queries during assessment. Its
reciprocal scales every component of the current stress, including thermal or
residual stress, and is not a load factor with those components held fixed.

`assess_affine_load` instead finds the first positive multiplier on
`sigma = residual + multiplier * load`, holding thermal/manufacturing stress
fixed. A strictly interior residual state and a nonzero linear unit-load
response are required. The cached quadratic form and a scaled, rationalized
root resolve the first crossing without iteration or allocation, including
initial unloading and load reversal. The returned multiplier retains a
reassessment of the actual boundary stress. Material axes, conditions and the
strength surface remain fixed; this is not nonlinear redistribution or an
interval bound. Analytical normal/shear intercepts, an independent multiaxial
polynomial, load rescaling and near-boundary residuals have focused regressions.

`tests/orthotropic_tsai_wu.rs` checks every normal/shear intercept, an independent
dimensional polynomial, combined tensile/shear crossings, three-axis convexity,
stress scaling over 560 orders of magnitude, compressive asymmetry, source
refresh and unsupported data. These are manufactured G0/G1/G3 checks, not
measured-material qualification. The criterion supplies no failure-mode
decomposition, automatic stiffness loss, fracture energy, progressive damage,
finite-strain stress conversion or experimental calibration.


## Native creep in layered mechanical and thermal sections

`NativeSectionEnergyLaw` admits `NortonCreep`, `NortonJ2` and their original source-resolved forms. The adapter invokes the native creep increment at its actual duration and preserves original physical and numerical losses. Pure primary resistance stores no energy. Joint creep/plasticity reconstructs elastic and isotropic-hardening storage through the original J2 energy owner using the shared elastic strain; creep accumulation is never substituted for plastic hardening history. The existing section quadrature, transverse reduction and work gates remain owners.

`plastic::section::thermal::source::ThermalSectionSourceFamily` carries complete native histories and per-volume work without a plastic-history projection. Existing thermoplastic families retain their original resolution and energy arithmetic. `CaloricCreepSource` binds any native `CreepSourceFamily`, including explicitly expanded standalone or joint primary creep, to the original same-card caloric chart and reference density. Admission retains fixed conditions, pinned Cp/density claims, a fully solid supported temperature interval and positive combined caloric/mechanical heat capacity. Transport selections must be instantaneous and depend only on owned temperature and fixed external conditions. A retained thermal creep temperature must equal the actual accepted caloric temperature; inspection never advances creep.

The thermal section integrates physical creep work as a subaccount of total physical heat. Native reciprocal expansion heat, free-energy change, mechanical entropy and numerical loss retain their separate meanings. One prepared interval resolves native laws per physical layer; each audited point replays the same old history and exact physical duration. Small-strain/rotation, fixed-reference geometry, thermal-domain and transverse-equilibrium limits remain in force. No tertiary damage, rupture, empirical calibration or finite-shell-forming claim is added.

Two new section regressions and native shell/coupled consumers exercise independent primary/secondary axial returns, distinct creep/plastic storage, sustained loads, relaxation and thermal history. These tests were added but not executed in the implementation environment; formatting and source checks are not compilation or numerical execution.

## Primary creep with independent static resistance recovery

`creep::recovery::RecoveringNortonCreep` keeps permanent `CreepState` and
nonnegative `resistance_strain` in separate history fields. The explicitly
selected rate is `alpha_dot = rate_Norton(q,T) h(r)`, with the existing primary
multiplier `h(r) = f + (1-f)(1+r/a)^(-m)` and
`r_dot = alpha_dot - k(T) r`. Independently supplied Arrhenius recovery
coefficients determine `k(T)`; they are never inferred from the flow activation
energy, yield behavior or melting point. This phenomenological resistance
model requires its own calibration and carries no stored energy in `r`.

For `c = delta_alpha` and `g = 1 + k(T) dt`, backward Euler gives
`r1 = (r0+c)/g`. The existing monotone scalar return and analytic tangent are
reused with private prepared scale `a*g` and kinetic memory `r0`. Permanent
strain validation, nondecreasing accumulated path and native stress/work
checks remain unchanged. The actual solved increment is retained even when
adding it to accumulated alpha rounds away. Unloaded recovery changes only
resistance, without erasing permanent deformation or generating heat. Physical
creep heat remains `q*c`. With zero recovery and `r=alpha`, stress, tangent,
history and energy reduce exactly to the original primary law.

`RecoveringNortonSource` selects and pins all thirteen original claims: the ten
primary Norton coefficients and three independent recovery coefficients. Every
thermal knot must retain the same reference coefficients and fixed external
conditions. Endpoint resolution refreshes the complete receipt bundle once;
constitutive trials perform no source queries or additional nonlinear solve.
The complete associated history passes through `CreepLaw`, sourced expansion,
native tetrahedra, `CaloricCreepSource` and thermal layered shells. Explicit
isothermal section adapters also retain its original duration and native work.

Analytical quadratic returns, unloaded exponential refinement, renewed creep
after a dwell, nonproportional tangents, heat closure, source admission and
native section history are covered by authored regressions. Compilation and
native execution remain pending the prescribed build runner. This increment
does not add recovery to joint plasticity, finite-deformation laws or frame
fiber dispatch, nor does it establish tertiary damage, rupture or service life.

## Joint recovering creep/plasticity and its exact axial restriction

`creep::recovery::plastic::RecoveringNortonJ2` extends the standalone recovery
model to simultaneous Norton flow and energetic isotropic J2 plasticity. The
original active-set return, including yield crossing, perfect plasticity and
positive hardening, uses recovering resistance as its kinetic memory. Permanent
creep and plastic tensors and their accumulated paths remain independent. The
shared `RecoveryStep` prepares the original monotone return once per physical
interval; no additional nonlinear solve, source query or history projection is
introduced in the stress/tangent path. Resistance evolves from the actual solved
creep increment even when it rounds away in accumulated permanent strain.

The joint work account retains elastic plus `H alpha_p^2/2` storage and separate
`q delta_alpha_c` creep and `yield_stress delta_alpha_p` plastic heat. Unloaded
resistance recovery neither erases permanent strain nor supplies heat or stored
energy. `NativeSectionEnergyLaw` queries the same original J2 storage owner.
`RecoveringNortonJ2Source` selects all fifteen original claims together, admits
the full thermal span through the existing source owners and refreshes the
complete pinned bundle once per interval. Its two permanent strain fields pass
unchanged through generic sourced expansion, solids and caloric sections.

`creep::recovery::uniaxial::UniaxialRecoveringNorton` is the exact free-sided
restriction of standalone recovering creep. It reuses the original axial root
with Young's modulus, retains signed strain, accumulated path and resistance,
and returns the compatible transverse strain and consistent axial tangent.
The original tensor recovery owner audits its work and heat. Source-resolved
axial laws retain all thirteen receipts and expose original physical coefficients,
never the privately timestep-scaled primary strain.

Authored regressions cover independent quadratic joint/axial returns, both J2
active branches, recovery-induced branch changes, nonproportional and condensed
tangents, zero-recovery equivalence, tiny solved increments, source admission,
native section history and heated solid dwell/reload. Scoped Rust parsing and
format checks do not constitute compilation or test execution. Finite-strain
recovery, kinematic hardening, tertiary damage and rupture remain outside these
bounded phenomenological laws.

## Hill plasticity in energetic layered sections

`plastic::section::energy::NativeSectionEnergyLaw` admits the original `Hill48`,
`OrientedHill48` and `ResolvedHill48` owners. The same checked native return
retains local-axis plastic strain and its Hill work-conjugate multiplier,
recoverable elastic/isotropic/kinematic storage, physical yield heat and separate
backward-Euler numerical loss. Layered plane-stress quadrature and the existing
native shell energy/dynamics paths now consume these laws without a J2 history
projection or inferred heat fraction. Coefficients and material axes remain fixed
over an isothermal account; invalid incoming or substituted endpoint histories
refuse through the original owners.

`tests/plastic_section/hill.rs` supplies authored directional coupon, rotation,
partial-thickness yield/springback, source-parity and history-refusal regressions.
These are numerical conformance cases, not measured rolled-alloy calibration,
finite-strain forming, texture evolution or fracture predictions.

## Source-bound thermal Hill plasticity

`plastic::hill::thermal::{ThermalHillSource, ThermalHill}` bind the original
oriented Hill return to the same fully solid caloric card and retained source
claims. The fixed Hill metric and frame preserve local plastic strain and its
work-conjugate accumulated multiplier. E(T) and energetic H(T) are affine, nu
and reference density are fixed, yield resistance may vary with temperature,
and selected constant isotropic expansion supplies stress-free strain. Every
source-line cell is admitted; unsupported axes, nonaffine stored coefficients,
metric drift, mismatched context/claims and nonpositive combined heat capacity
refuse. The capacity gate uses the caloric owner's actual sensible-capacity
bounds, not an endpoint enthalpy secant.

Old history is validated at its actual old temperature before the original
endpoint return. A valid state may therefore yield under heating at fixed total
strain. The original physical yield work, frozen-history coefficient storage,
mechanical entropy and signed reversible heat close the discrete first law;
backward-Euler loss remains separate. Native point evaluations retain the
original analytical tangent and fixed-array scalar return. Layer preparation
resolves endpoint sources once before native equilibrium evaluations.

The direct `ThermalSectionSourceFamily` implementation preserves complete Hill
history in existing layered sections and thermoplastic shells, including their
motion/contact/thermal owners. `tests/plastic_hill_thermal.rs` covers constrained
softening, an independent rotated tensor free-energy/entropy account, unequal
thermal layers and source/capacity refusals. Scope is fixed-axis small-strain
solid plasticity, not finite-rotation texture, melting, fracture or an empirical
calibration. No J2 property/history substitute or heat fraction is inferred.

### Complete free-sided recovering creep/yield and thermal histories

`creep::recovery::plastic::uniaxial::UniaxialRecoveringNortonJ2` restricts the
same simultaneous recovering tensor law to free-sided axial stress. The original
joint axial return uses Young's modulus in its balance and retains the actual
solved creep increment. Signed creep/plastic strains, both accumulated paths
and recovering resistance stay independent. The native recovering tensor law
checks the compatible lateral strain, stress, histories and full work account;
there is no fitted axial replacement or extra stress/tangent root. All fifteen
original source receipts remain attached to source-resolved fibers.

`creep::thermal::uniaxial::recovery` supplies standalone and joint thermal
recovering fibers through one shared implementation. Each accepted state retains
its actual temperature alongside the complete mechanical history. Old compatible
strain is reconstructed at that retained temperature; trial free strain and both
Arrhenius kinetics use the explicit endpoint temperature. The original thermal
tensor owner supplies elastic/hardening storage, physical creep/plastic heat,
signed reciprocal heat, `T` times mechanical entropy and numerical loss. The
source constructor resolves this tensor once and builds its axial reduction
from the same resolved mechanical law and selected expansion claim.

Manufactured regressions specify independent axial return/tangent equations,
source continuation, thermal loading and softening, stress-free hot dwell,
reloading, all retained histories, the exact native tensor energy account and
condensed tangent. Scoped parsing, formatting and scalar reference checks do
not replace native compilation or execution. These are bounded small-strain
material laws with prescribed endpoint temperature; a thermal fiber alone does
not solve heat transport, add damage or establish calibrated service life.

## Source-bound mixed and thermal Hill hardening

`hill48_mixed_requirements`, `resolve_hill48_mixed` and the explicit
`ResolvedHill48::with_linear_kinematic_hardening` opt-in select the independent
`hill_kinematic_hardening_modulus` C [Pa]. Native Prager backstress is X=C eps_p
and recoverable storage is C eps_p:eps_p/2, including tensor shear twice. J2
uniaxial H_k, scalar tensile slopes and missing directional data cannot stand in
for this source. Original source claims, fixed material axes, scalar return,
analytical tangent and local plastic history remain owned by the native law.
The ordinary constructor retains C=0 until the explicit opt-in.

`ThermalHillSource::with_linear_kinematic_hardening` pins a nonnegative affine
C(T) over every source-line cell and refuses repeated selection or aliasing.
Each endpoint refreshes that same claim. Mechanical entropy includes
-C' eps_p:eps_p/2; frozen-history coefficient work includes
(C1-C0) eps_p0:eps_p0/2. The complete native yield heat and numerical loss remain
separate. C may reach zero without dividing by it, rescaling a backstress or
replacing retained plastic strain. Affine C has no additional frozen-history
second-temperature derivative in the existing capacity admission.

`tests/plastic_hill_source_mixed.rs` covers source/native cyclic parity, exact
zero-C equivalence, directional Bauschinger thresholds/tangents and source-role
refusals. `tests/plastic_hill_thermal/mixed.rs` covers native section reversals,
independent tensor storage/entropy through C=0, original claim pins and corrupt
history refusal. These authored conformance cases do not establish nonlinear
recovery, texture evolution, rate dependence or measured-alloy calibration.

### Source-bound thermal Hill saturation

`hill48_voce_requirements`, `resolve_hill48_voce` and
`ResolvedHill48::with_voce_isotropic_hardening` select independent
`hill_voce_saturation_increment` Q [Pa] and `hill_voce_saturation_rate` b.
The rate is per Hill work-conjugate alpha; J2-equivalent-strain coefficients
cannot substitute for directional data. Native material axes, optional Prager
hardening and the complete original source bundle remain intact.

`ThermalHillSource::with_voce_isotropic_hardening` admits nonnegative globally
affine Q(T) and positive constant b over every original source-line cell.
Pinned claims survive endpoint refresh and either Prager/Voce selection order.
The original Hill return, consistent tangent and full local history serve
thermal plane-stress sections and their existing shell consumers.

With g(alpha)=alpha-(1-exp(-b*alpha))/b, extra free energy is Q*g,
mechanical entropy is -Q'*g, and frozen-history coefficient work is delta_Q*g.
The shared stable primitive evaluates signed coefficients directly, preserving
representable thermal work even when the unit-Q shape underflows. A selected
Q=0 endpoint or reference retains nonzero Q' entropy. Affine Q and fixed b add
no frozen-history heat-capacity derivative; physical yield heat and numerical
loss keep their native, separate accounts. Unselected/identically zero-Q paths
retain the original arithmetic. Missing or aliased roles, negative/kinked Q,
varying b and incomplete source support refuse. Focused regressions in
`tests/plastic_hill_thermal/voce.rs` cover source/native parity, thermal section
history, rotated first-law accounts and these limits. This is fixed-axis
small-strain constitutive verification, not empirical alloy calibration,
anisotropic texture evolution, melting or thermal annealing.

## Source-bound phosphorescence kinetics

`phosphorescence::PhosphorescenceStatePoint` resolves the original singlet and triplet excitation energies and six independent first-order rates into the existing `fs_thermal::luminescence::transfer::SingletTripletDecay`. Energies use joules per excitation relative to the same ground state, with `E_s >= E_t > 0`; radiative, nonradiative and forward/reverse crossing rates use inverse seconds. An explicitly zero rate is admitted, including a stable store. Missing rates never receive defaults, and lifetimes or branching yields are not independently invented. The native kinetics owner retains positive evolution, integrated events, crossing heat and reverse-crossing thermal withdrawal.

The adapter preserves the complete original material bundle and its receipts, including unrelated optical/thermal properties. Exact scalar quantities and fixed external source contexts are rechecked on ingress. Original temperature and explicitly typed concentration/composition curves or surfaces may supply changing rates, including quencher-dependent nonradiative loss. Refresh pins the same card, coordinates, requirements and selected claims. An excited-history owner must retain fixed state energies or account for their parameter work. Source authors own the physical validity of reverse rates. Nonradiative channels mean locally thermalized loss; excited-quencher storage and chemistry are absent.

This is a dilute two-state constitutive snapshot, with no trap distribution, population saturation, annihilation, bleaching or measured-material claim. Photon absorption, spectra and transport remain their existing owners. Source resolution occurs at preparation/refresh; reading the native law allocates nothing and performs no source query. Focused manufactured tests in `tests/phosphorescence.rs` check exact roles, complete-bundle retention, pinned thermal refresh, Stern–Volmer-type quenching against independent exponential decay, explicit stable stores and wrong-unit/history/source-domain refusals.

`phosphorescence::pump::GroundStatePumpStatePoint` resolves two separately authored ground-to-singlet/triplet absorption cross sections in square metres per emitter. Vacuum wavelength is an explicitly typed SI-metre coordinate; an optional owned absolute-temperature coordinate refreshes at the thermal state. Original native queries retain the complete material bundle, exact selected claims, coordinate conventions and source support, with no fitted spectrum or inferred molecular count. Both target cross sections may be nonzero at one wavelength: a transport owner must use their sum for attenuation and their separate shares for excitation. They are not mass absorption, molar extinction or a bulk-index conversion. Missing/negative roles, wrong units, unowned field/history dependence and unsupported refreshes refuse; explicit zero is retained. The convenience wavelength/temperature refresh holds every other external coordinate fixed and requires temperature optionality to match its declared owner. This source snapshot supplies coefficients only; photon depletion and finite-population/energy evolution remain their existing transport and kinetic owners. Manufactured tests in `tests/phosphorescence/pump.rs` cover native wavelength/temperature surfaces, overlapping-source claim pinning, whole-bundle retention, coordinate/domain refusal and permissive-bundle revalidation.

### Triplet-pair annihilation source coefficients

`phosphorescence::annihilation::TripletAnnihilationStatePoint` resolves the
original `luminescence_triplet_pair_event_rate` in m^3/s and
`luminescence_triplet_pair_singlet_yield` in [0, 1]. The convention is explicitly
`R = k N_t^2 / V` pair events per second: each event removes two triplets. A
triplet-loss coefficient must be converted by its author, and a molar rate
coefficient is not interchangeable with this count-based convention. Neither
diffusion rates nor spin branching factors are inferred from material names.

The snapshot preserves the complete original material bundle, exact selected
claims and fixed external conditions. An explicitly owned absolute-temperature
axis refreshes both coefficients at the native thermal state. Wrong units,
negative rates, probabilities outside [0, 1], unsupported source domains and
unowned field/history response refuse. Zero remains an explicit physical value.
The caller supplies the active volume to `law(volume_m3)`; the native
`TripletAnnihilation` owner checks that each productive pair can supply its
singlet energy and owns pair counts, populations and heat. Coefficient refresh
does not evolve or replace those histories.

Four focused G0/G1/G3 tests in `tests/phosphorescence/annihilation.rs` cover
original thermal interpolation, the physical factor of two and exact isolated
pair kinetics, active-volume admission, whole-bundle and claim retention,
explicit zeros, permissive-bundle revalidation and unsupported refreshes. All
14 phosphorescence source tests passed in the existing narrow actual-source
validation build. These manufactured claims establish implementation behavior,
not empirical coefficients for a particular phosphor.

## Bounded tabulated plastic hardening

`plastic::tabulated::TabulatedHardening` retains explicitly supplied absolute
yield stress versus accumulated equivalent plastic strain. The origin is zero
plastic strain, strains increase strictly, and positive strengths are
nondecreasing; plateaus are allowed. Slopes and exact hardening primitives are
prepared once. Returns and lookups allocate no material-point memory. Support
ends at the final knot: no constant extension, slope extrapolation, smoothing,
or automatic conversion from engineering total-strain measurements occurs.
At an interior kink the right-hand derivative is selected; the final knot uses
the last measured segment. No two-sided differentiability is claimed there.

`TabulatedJ2` reuses the native small-strain J2 tensor flow, plastic history,
Prager translation and algorithmic tangent. `UniaxialTabulated` restricts the
same scalar consistency relation to free-sided axial stress, including lateral
plastic flow and Poisson contraction. Native plane-stress and layered-section
consumers retain the full original history and energy account.

`FiniteJ2::with_tabulated_hardening` explicitly selects a Kirchhoff-stress
table for the original multiplicative Hencky model. Its initial yield must
match, its original linear modulus must be zero, and table support must cover
the declared history domain. Linear, Voce and table selections cannot silently
replace one another. The original inverse plastic metric, analytic matrix-log
tangent, finite-strain limits, plane-stress thickness solve and volume/area
conventions remain authoritative. Native Perzyna, joint Norton/J2 and ductile
work paths use the same table rather than reverting to a linear backbone.

Both stress conventions select energetic isotropic hardening explicitly:
stored hardening is the integral of `Y(alpha)-Y(0)`; yield dissipation is
`Y(0)*delta_alpha`. Native viscous/creep work remains separate. A measured flow
curve alone does not identify this partition experimentally. Numerical work
defects are not material heat. These fixed-condition tables do not supply
temperature evolution, annealing, anisotropic texture, a forming-limit curve,
necking, or experimental calibration. G0/G1/G3 regressions target knot and
plateau returns, analytical scalar relations, independent tensor derivatives,
history under reversal, energy partition, and domain refusals in
`plastic_tabulated`, `plastic_tabulated_reduction`, and
`plastic_finite_tabulated`.

### Bonded common-Poisson circular transverse shear

`section::transverse_stress::layered::LayeredCircularTransverseShearRecovery`
retains distinct original E/G and material identities in perfectly bonded
concentric isotropic layers. All layers must explicitly share one Poisson
ratio. Exact radial warping and shear-traction transmission use the shared
pivoted tridiagonal owner in O(layers) preparation; each material-point query
is allocation-free O(1). Original bore/exterior radial tractions vanish, the
integrated stress recovers the supplied transverse resultant, and each side
of an interface retains its own strain and stress. The whole selected annulus
engineering-shear bound includes analytic interior extrema of both stress-map
eigenvalues. Unsupported or numerically unresolved input refuses without
clipping, homogenizing moduli or changing original shear correction factors.

The original homogeneous owner retains its admission boundary. The layered
owner is an explicit alternative for prismatic small-strain observations; it
does not replace reduced mechanical storage or supply varying-Poisson-ratio,
interface slip, plasticity, thermal mismatch, end effects or cross-section
inertia corrections. Kernel tests compare closed-form solid/tube limits,
heterogeneous interface continuity, exact resultants and material-subdivision
invariance. Source/native regression definitions add the independent two-layer
center solution, source-role refusal and actual rod fatigue consumption.

## Original sensible temperatures for ordered contact sources

`thermal::interface::ThermalInterfacePoint::at_temperatures` and
`resolve_resistance_at_temperatures` admit actual temperatures supplied by
original sensible-storage owners without manufacturing an equilibrium phase
state. Both temperatures must be finite positive absolute temperatures on
distinct existing axes. Original surface A/B identities, optional declared
mean temperature, typed fixed conditions, static-source admission and exact
claim selection retain the existing interface owner. Any declared liquid-mass
fraction axis refuses this path; the original phase-state methods keep their
full phase/card checks.

`CaloricConductionLayer::supports_material_card` checks exact evidence
membership for source contacts. Cp, Cv and fixed-volume layers retain their
original material card. An assembled body admits an original constituent card
or its explicitly selected effective-conductivity card. The caller declares
which card describes the exposed surface; membership alone does not infer
surface topology, constituent exposure or a homogenization law. The native
contact adapter in `fs-conduction` also binds the actual original caloric
objects and validates endpoint mass/volume before using this seam.

## Calibrated Walker stress-life mean stress

`fatigue::walker` adds an explicitly supplied Walker exponent and a closed
load-ratio interval to the shared S-N owner. Every counted full or endpoint half
cycle maps to the original curve's retained reference ratio using
`a_ref = a * ((1 - R_ref) / (1 - R))^(1 - gamma)`. The supported domain is
`0 <= gamma <= 1` and `-1 <= R_min < R_max < 1`, with R_ref inside that same
interval. Positive-amplitude cycles require a tensile maximum. Zero-amplitude
observations add no fatigue cycle; the correction supplies no static-strength
assessment. The original finite-life and endurance declarations still own the
corrected amplitude's admissibility.

`MeanStressCorrection::Walker` works with fully reversed and non-fully-reversed
curves without changing their knots. Normalized geometric means avoid powers of
dimensional stresses, with direct gamma=0, 1/2 and 1 limits. Domain checks use
stable `1-R`; reference-equivalence tolerance is relative to amplitude so
ratios near unity cannot hide a large amplitude correction. No default gamma,
compression credit below R=-1, extrapolation or constitutive yield adjustment
is inferred. The reference normalization follows NASA/TM-2005-213886/REV2
section 2.3 equations (23)-(25).

`resolve_walker_sn` retains the separate `fatigue-walker-gamma` claim, selected
answer and receipt alongside all original S-N observations. Gamma is a fixed
scalar or an original temperature-only curve; provenance, observation population
and complete typed support must match the S-N claim except for the explicitly
declared R interval. The returned curve requires its exact retained Walker
policy. A different authored gamma or support cannot silently inherit that
source binding.

Six focused regressions in `tests/fatigue_walker_stress_life.rs` cover analytic
limits, actual reference ratios, variable means, source receipts/population,
transactional refusals, scale invariance and near-unity ratio conditioning.
These are synthetic implementation checks, not measured component-life evidence.

### Fugacity-driven real-gas molecular diffusion

`nonideal::energy::parcel::transfer::species::diffusion` binds the original
Peng–Robinson/NASA catalog to the existing Maxwell–Stefan friction owner and
explicit selected binary material claims. Both physical density/composition
traces are evaluated at the same declared central temperature; no PT root is
reselected and no ideal p/(RT) density enters. Original critical constants,
interaction pairs, calorics, transport and source packs must match. Common
positive species support is required; one-sided logarithmic zeros refuse.

With g_i=Delta ln(x_i p phi_i)/distance, the mass-frame force is
d_i=x_i g_i−Y_i sum_j(x_j g_j). Finite differences are not assumed to obey
Gibbs–Duhem. Removing a mass-proportional force leaves d unchanged. Scaled
mass ratios, mean mass and the mechanical-gauge cancellation remain
double-double, including the independent force-power reduction. The reported
raw/reduced forces, original row and mass-frame residuals, positive pair
friction and reconstructed constrained chemical power remain separate.
Friction entropy is not deposited as additional heat.

Every selected binary MS claim is pinned and refreshed at actual endpoint and
central T/p. This rung supports T/p-only coefficient laws; original data must
qualify their composition-independent use over the intended mixtures.
Composition/humidity-dependent coefficient claims refuse. PR chemical
potentials do not supply dense-fluid diffusivity, infer ideal pD scaling or
promote pure-gas measurements into mixed-gas authority.

The network samples every face from one unchanged inventory, shares dense
scratch among equal-size catalogs, and invokes the original simultaneous
partial-enthalpy/PR reconstruction once per affected parcel. Actual final source
support and authored T/p limits are rechecked. The original closed-species
caloric/entropy audit is reused through borrowed physical inventory views, with
the original PR density-chart departures added once. Whole-system mass,
species, source-relative U and actual S must pass before publishing any endpoint.
The ideal-gas audit retains its original formulas and reduction order.

This is first-order explicit transport with a caller-bounded local isothermal
approximation. It supplies no Soret/Fourier/thermal-transpiration model, bulk
momentum, flash, phase-stability theorem, implicit stability or global/spatial
error bound. Homogeneous mechanical/caloric admission alone does not establish
stable single-phase equilibrium.

Six exact production constitutive tests passed natively against the genuine
current fs-thermochem/fs-math owners: binary mass-frame and ternary friction
identities, non-dyadic large mechanical-gauge invariance, trace forces, an
independent dense PR/binary oracle, and refusal/cancellation. Five material-face
and finite-network regression groups cover original sources, conservation,
reference/order invariance and atomic refusal; they remain unexecuted because
the full dependency build is blocked before fs-material. All changed Rust
sources parsed; this is not a full crate type-check.

### Adaptive source-bound real-gas diffusion windows

`nonideal::energy::parcel::transfer::species::diffusion::network::adaptive`
advances complete finite PR networks through the original gas step-doubling
controller. A narrow inventory/candidate type parameter preserves the ideal,
binary and source-conjugate thermal paths' arithmetic and accepted ledger.
Real-gas probes retain the original PR endpoint reconstruction, source
qualification and finite-interval audit; no second integrator or EOS is used.

Each original species' actual extensive mass and the original PR endpoint
temperature are compared between one full step and two half steps. Only the
two-half solution is accepted, with first-order local error divisor one and
no extrapolated state or trace floor. Temperature control avoids error scales
dominated by arbitrary formation-energy references. The estimate is local,
not a global temporal or spatial accuracy certificate.

Accepted gross species/partial-enthalpy transfers are accumulated in canonical
physical-edge order; reverse flow owns a separate row. Rejected probes change
no accepted inventory. Original source data and fugacity are refreshed on
every probe. Matrix scratch is reused without retaining a physical state.
Both complete candidate evaluations and actual real-gas face calls have
cumulative finite budgets, including rejected, failed and partially sampled
probes. Geometrically closed faces consume no constitutive calls.
Cancellation and any terminal refusal return no partially advanced network.

Only typed gross withdrawal, endpoint caloric-bracket and finite-step entropy
refusals may shorten a trial. Source/domain, mechanical/phase, conservation and
nonconvergence refusals remain fatal. The accepted source-relative energy
ledger AND independent original-to-final species/U/S audit must satisfy the
same fixed budgets, never allowances multiplied by step count. The original
PR density-chart energy/entropy departures and species-reference cancellation
are retained. The local isothermal-face approximation, positive common
species support and T/p-only binary coefficient qualification still apply;
no implicit stability, Soret/Fourier law or phase-stability claim is added.

Five G1/G3/G4 integration groups cover native fine-state/enthalpy equivalence,
reference invariance, refined species/temperature convergence, actual rejected
work, fatal-source separation and late cancellation. They are authored but
unexecuted while the full dependency graph remains blocked before fs-material.
All six changed Rust sources passed pinned-toolchain format/syntax checks;
this is not a full crate type-check.

### Nonisothermal real-gas molecular and reciprocal heat transport

`nonideal::energy::parcel::transfer::species::diffusion::thermal` adds
source-bound Soret/Dufour and Fourier transport on actual PR states. Original
NASA increments and actual ln(x_i p phi_i) form Delta(mu_i/T); no ideal density,
held-temperature surrogate or new PT root enters. The chosen h*_i is the
existing NASA divided caloric plus symmetric actual PR partial residual
enthalpy. That residual term is a finite-face approximation, not an exact
Gibbs–Helmholtz integral across a temperature jump.

The mass-frame effort is [Delta(mu_i/T)−h*_i Delta(1/T)−M_i H_i Delta(1/T)]/R.
The same h* converts original donor enthalpy carriage through additional heat
q_Fourier + sum(H_i M_i N_i) + sum(N_i(h*_i−h_i,donor)). Formation references
and the common H datum cancel before arithmetic; mechanical gauge removal
reuses the original double-double projection. Independent actual donor-energy/
chemical power must match original MS pair friction plus Fourier entropy.
At equal T and H=k=0, species transport matches the original PR diffusion law;
centered residual enthalpy can still differ from its upwind energy carriage.

Every binary coefficient, signed H_i and zero-species-flux conductivity is an
explicit pinned T/p-only source claim qualified over the intended mixtures.
PR thermodynamics supplies no missing kinetic data. The unchanged positive
common-support, original-source and homogeneous PR admission boundaries apply.
Local PSD/reciprocity does not establish phase stability or finite-step safety.

One shared real-gas network owner combines species and additional heat before
its original native endpoint reconstruction and actual species/U/S audit.
The original temporal controller also checks each thermal face's signed and
both gross-sign mass/heat integrals. Rejected probes spend actual cumulative
face work but no inventory; fixed whole-window budgets and cancellation remain.

Six new native constitutive tests passed with genuine original NASA/PR/MS
owners, including independent source power, equal-T parity, Onsager cross
response, temperature-only refinement and large-reference invariance.
Three material/network integration groups are authored but unexecuted because
the full graph remains blocked before fs-material. All ten changed Rust files
parsed/formatted; no full material crate type-check is claimed.

### Source-conjugate gas thermodiffusion

`thermochemical_gas::energy::parcel::transfer::species::diffusion::thermal` prepares an original NASA species catalog and complete pinned binary-diffusivity, mass-specific heat-of-transport and zero-diffusive-flux conductivity sources. Actual traces and the central face requery those same claims. Only owned canonical temperature/pressure dependence is admitted. A reusable Maxwell–Stefan workspace supplies the mass-average molecular flux; centered heat-of-transport coefficients couple Soret forcing to reciprocal heat. Original NASA enthalpy/entropy increments form divided caloric enthalpies, so heat added to existing donor-enthalpy transport includes the finite-distance caloric correction without subtracting arbitrary formation references. Common exact zeros stay absent, positive trace species are retained, and one-sided logarithmic zeros refuse.

The thermal network samples every edge from one unchanged snapshot, then combines molecular masses and additional heat with all existing authored transport paths before one inventory/caloric reconstruction per parcel. Gross original-species withdrawals, mass/energy budgets, pinned source support on the combined endpoint and actual total entropy are checked before publication. This is explicit first-order finite-volume transport, with authored local T/p limits; it supplies no implicit stability guarantee, momentum closure, thermal transpiration, source uncertainty bound or empirical coefficient inference. Six analytic/source/native-network tests cover Soret counterflow, reciprocal heat, variable-Cp divided-enthalpy entropy pairing, common-heat and formation-reference invariance, exact-zero/trace behavior, original endpoint source support, cancellation and coupled conservation.

### Adaptive source-conjugate gas thermodiffusion

`species::diffusion::thermal::network::advance_candidate` completes a closed
nonisothermal gas-network window through the existing molecular temporal owner.
Every full/half probe runs the original thermal snapshot transaction: source
T/p/composition, Soret/Fourier/reciprocal heat, original donor species enthalpy,
other declared transport paths and actual endpoint source admission remain
coupled. Only accepted half steps contribute inventory and integrated accounts.

Independent local criteria compare species inventories, temperatures and each
original thermal link's signed and both gross-sign constituent masses and
additional heat. The latter is separate from species enthalpy, not a second
charge for total thermal carriage. Compensated coarse/fine comparisons retain
low parts and prevent opposing links or within-interval reversals from hiding
transport error. The optional report criterion is absent for original endpoint-
only controllers, whose acceptance arithmetic is unchanged. Other authored
paths share endpoint control and accounting but have no per-path integral
criterion through this API. These first-order local estimates are not certified
whole-window, spatial or source-uncertainty bounds.

Rejected probes spend the original evaluation budget but no accepted inventory.
Only native withdrawal, endpoint entropy and caloric-bracket refusals permit
cutback; source/domain failures remain fatal. Original whole-window species,
mass, source-relative energy and entropy budgets apply once. Last-step source
and flux receipts describe the start of the final accepted half step; actual
final material support is admitted separately. This is explicit transport,
without an implicit stability or bulk-momentum claim.

The exact production integral-comparison helper passed three native tests with
the genuine fs-math double-double owner. Five analytic thermal-relaxation,
transfer-accuracy, caloric-reference, equilibrium and transactional-refusal
integration groups are authored but unexecuted because the dependency build
remains blocked before fs-material. Source parsing is not crate type checking.

## Original measured hardening, rate dependence and thermal source transport

The tabulated J2 family retains the original accumulated-equivalent-plastic-strain
knots and bounded nondecreasing absolute yield ordinates. Its isotropic free
energy is the exact integral of resistance above the initial yield stress.
Isothermal source preparation distinguishes the exact Cauchy and Kirchhoff
property names; total or engineering tensile strain is not silently converted.
Original material identity, claim, policy, uncertainty, fixed context and each
queried observation remain inspectable. Prepared constitutive iterations do not
query the material database or allocate observation receipts.

`TabulatedOverstressJ2` selects linear or power-law associative Perzyna flow at an
explicit physical duration. Linear flow uses the shared piecewise-linear scalar
return; nonlinear flow locates one measured segment and uses a bounded monotone
log-increment solve. The original J2/Prager tensor return, consistent tangent and
history remain authoritative. Layered-section energy receives the same duration.
Recoverable table storage, initial-yield loss, viscous loss and integration loss
remain separate; only physical yield/viscous loss is heat. Machine-scale
last-knot consistency roundoff is bounded inside the measured domain without
admitting a physical overshoot.

Finite thermal tables retain fixed strain knots with affine temperature-dependent
Kirchhoff ordinates. Endpoint admission proves positive nondecreasing resistance
over the declared interval. The exact signed temperature derivative of table
storage contributes to entropy, internal energy, coefficient work and reversible
caloric exchange, including an endpoint curve that becomes a plateau. It adds no
fixed-history heat-capacity correction because this hardening free energy is
affine in temperature. Perzyna and simultaneous primary/secondary Norton flow
retain the same table and native plastic/creep histories.

`resolve_thermal_tabulated_hardening` requires original complete strain/absolute-
temperature grids, measured family endpoints and affine columns across every
intervening source row. All original sample receipts survive, including samples
outside a narrower requested thermal interval. TabulatedOnly observations remain
discrete evidence for explicitly selected constitutive interpolation. LinearInside
triangles additionally require a compatible affine plane and original continuous
rectangle support; bilinear and triangular interpolation are not interchangeable.
Holes, thermal derivative kinks, single-temperature derivatives and extrapolation
refuse.

The finite thermal source has separate tabulated affine, Perzyna, Norton/J2 and
primary Norton/J2 selectors. Thirteen mechanical/thermal parameters remain
explicit; competing scalar yield/hardening values and their slopes are rejected.
The same-card surface supplies those coefficients and derivatives. Law/state
versions retain the existing one/two distinction for simultaneous creep. The
original source result is shared across family clones. Added material, source and
structural regression fixtures are synthetic numerical cases, not measured
calibration. Rust syntax and scoped whitespace were checked; native compilation
and regression execution remain unverified in the resource-constrained session.


## Analytic frequency dispersion of original Maxwell spectra

`GeneralizedMaxwell::modulus_frequency_derivative` returns the original
storage and loss slopes with respect to angular frequency [Pa s]. Scalar
frequency retains signed parity: storage is even, loss is odd, storage slope
is odd and loss slope is even. `TensorProny::complex_modulus_frequency_derivative`
returns those full Mandel tensor slopes plus `d(C''/omega)/domega` [Pa s²],
including signed anisotropic coupling entries. The original nonnegative
applicability band is enforced, including its endpoints; the algebraic law
supplies the admitted one-sided endpoint slopes. At admitted DC, storage and
effective-viscosity slopes are zero and loss slope is `C_kv + sum(C_j tau_j)`.
Kelvin viscosity contributes its original constant loss slope and exactly
zero effective-viscosity slope.

The derivatives use the authored branch stiffnesses and relaxation times,
with binary-scaled products and compensated accumulation in branch order.
No spectral refit, perturbation solve or frequency differencing is performed.
Reciprocal evaluation and delayed physical scaling retain representable weak
storage, loss and derivative entries when `omega*tau` or its square cannot
be represented. Nonfinite inputs or required unrepresentable derivatives
return `ViscoError::Parameters`; tensor band violations remain
`ViscoError::OutOfBand`. A finite primal query does not require a finite
tangent, and the scalar modulus slope does not require the additional
effective-viscosity slope to fit. The small fixed tensor queries allocate no
workspace and introduce no worker-dependent reduction order.

`tests/visco_frequency.rs` contains ten synthetic G0/G1/G3 regressions for
closed-form branch slopes, full anisotropy and Kelvin terms, independent
native-primal finite differences, signed scalar parity, tensor band admission,
DC and representability boundaries, weak mutual entries over extreme scales,
subnormal relaxation times and cancellation between relaxation slopes.
These integration tests are authored but have not run in the implementation
environment: the real dependency graph is blocked while compiling
`asupersync`. Six inline tests of the actual shared spectrum kernel passed
with the pinned native compiler. Neither result establishes measured-material
qualification or replaces the original hereditary states during time
integration; `C''/omega` here is a harmonic constitutive quantity.

### Binary-liquid diffusivity at distinct endpoint temperatures

`MaterialBinaryDiffusivity::resolve_nonisothermal_face` refreshes the same pinned original claim at each endpoint's own temperature/composition and at their arithmetic mean face. It preserves held external coordinates and original joint interpolation/support. Equal compositions do not reuse a receipt across unequal temperatures. The isothermal `resolve_face` delegates to this path with identical temperatures, retaining its existing semantics. Two source tests exercise distinct-temperature receipts and triangular joint support; two spatial integration tests verify the actual source face coefficient and atomic refusal when a local endpoint leaves the original source support. This adapter supplies molecular kinetics only; heat carriage and explicitly authored Soret/Dufour coefficients remain the thermochemical/spatial owner's responsibility.


## Complex plane-stress harmonic frequency derivatives

`TensorProny::plane_stress_complex_frequency_derivative` differentiates the
original complex transverse-stress constraint in all four ordered material
planes. It returns signed storage/loss slopes in the same engineering-shear
coordinates as `plane_stress_complex`. The primal material evaluation,
transverse solve and reduced passivity admission are shared once; the public
primal query retains its original pair result and admission semantics.

For `X=C_oo^-1 C_oi` and `T=[I;-X]`, the exact complex-symmetric Schur derivative
is `Q_omega=T^T C_omega T`, with ordinary transpose and no conjugation. This
retains the changing transverse-strain coupling between storage and loss
without a second factorization or a derivative solve. The shared material
aggregation requests only modulus slopes: an unrepresentable derivative of
effective viscosity does not invalidate a finite plane-stress slope.

Each reduced entry has its own factored, compensated complex contraction.
A large eliminated normal derivative therefore cannot globally rescale away
an unrelated weak in-plane derivative. Physical units are restored after
signed accumulation; finite-precision cancellation is not a relative-error
certificate. The result preserves exact symmetry and the original frequency
band. Derivatives need not themselves be passive or positive matrices.

Six integration regressions cover independent full stress constraints with
noncommuting anisotropic branches in every plane, native frequency differences,
complex isotropic formulas, independent viscosity-slope admission, an in-plane
`1e-100` slope beside an eliminated `1e300` slope, and typed source refusals.
They remain unexecuted with the full dependency graph. Four native tests of the
actual production projection module passed against the real `fs-math` library,
including complex product signs, overflowed products that cancel, weak retained
responses and minimum-subnormal coefficient scaling. This is a material query
for plate/shell consumers; it does not add an assembled plate/shell frequency
solver or replace the original hereditary material history.

## Inertial dielectric resonance at fixed material conditions

`field_storage::dielectric::lorentz::ResolvedLorentzMaterial` retains a complete
original material snapshot and explicitly named high-frequency permittivity,
positive oscillator increments, positive natural angular frequencies, and
nonnegative damping rates. Natural frequencies require the semantic rad/s
quantity; Hz and dimension-only inverse seconds refuse. No oscillator or loss
is inferred from missing source data. The empty bank recovers the instantaneous
lossless law.

Each pole obeys `P'' + gamma P' + omega0^2 P = delta_epsilon omega0^2 E`.
Its `exp(+i omega t)` susceptibility uses scaled factored detuning, retains
negative real permittivity and passive loss, and refuses the actual undamped
periodic pole without artificial broadening. Original polarization and current
retain positive restoring plus kinetic storage. Normalization maps this same
law to the existing exact `fs-modal` forcing owner; there is no second integrator.

Prepared linear-field and continuous-tone transitions advance one to three
polarization components without warm-step allocation or transcendental calls.
The harmonic convention is `E(t)=Ec cos(Omega t)+Es sin(Omega t)`, equivalent to
phasor `Ec-i Es`. Exact lossless resonant growth remains finite over finite time.
Independent supplied work and frictional heat accompany storage and balance
residuals. Failed immutable trials publish no partial state.

Ten native tests of the actual production oscillator modules pass, covering
analytic response, high-Q and true lossless poles, optical-scale ringdown,
linear and harmonic subdivision, resonant secular growth and independent cycle
heat. Five material-source integration tests are authored; their source-card
fixtures compile against real dependencies, but the full material adapter tests
remain unexecuted without the complete constellation dependency environment.
This is a linear isotropic fixed-condition law, not measured qualification,
free-carrier Drude dynamics, nonlinear saturation or a transient Maxwell solver.

## Transport sources across ionization population boundaries

`ionized_gas::IonizedGasTransportBank` retains momentum laws independently of
the currently present carriers. Every charged species in the original
thermodynamic bank has exactly one authored law. Kinetic electron drag retains
one original momentum-transfer cross section for every neutral source row,
including initially absent targets. Canonical carrier and neutral indices are
prepared once; endpoint evaluation only integrates populated neutral targets.

Material-backed laws retain their original complete resolved requirements and
immutable cards. The same source-refresh helper used by `at_frozen_state` and
`at_reacting_state` requeries actual Th, Te and pressure without replacing
selected claims or widening a common-temperature schema. Neutral drag uses
actual target densities and the weak-drift Maxwellian momentum moment. Optional
Coulomb drag rebinds the original screening, cutoff and limits to the actual
created carriers and reevaluates every charged pair.

`evaluate` accepts population support changes at the original fixed volume and
thermodynamic bank. An exactly neutral endpoint returns `None`; its modeled
free-carrier current and heating vanish. Positive trace populations are never
floored away: they produce the native `IonizedGasConductivity` or an explicit
source/domain/representability refusal. The resulting owner feeds existing DC,
Joule, harmonic and magnetized consumers. No positive PDE coefficient is
fabricated for an insulating endpoint, and no neutral dielectric response is
inferred. Species changes are prescribed by the caller's chemistry or transport
owner; this constitutive query does not itself evolve or certify reactions.

Four focused G0/G1/G3 regressions cover exact neutral-to-carrier-to-neutral
transitions and trace carriers, independent constant-cross-section DC and
Drude response, conserved H recombination that activates a new neutral target
in an inert buffer, original Coulomb rebinding, incomplete/foreign/duplicate
source refusal and atomic cancellation. These material integration tests remain
unexecuted without the complete native dependency chain. Source/bank admission,
format checks and static API review do not constitute that missing test pass.

## Prepared full plane-stress harmonic strain recovery

`TensorProny::prepare_plane_stress_harmonic` reuses the original full complex
constraint solve to retain engineering storage/loss and all six Mandel strain
components in every supported ordered material plane. Physical peak phasors
use `exp(-i omega t)`: the transverse constraint is conjugated consistently
with `C' - i C''`. At DC strains must be real. Preparation retains the original
frequency admission and strain ceiling; each response remains immutable.

`PlaneStressHarmonic::recover_strain` uses compensated, factored contractions
without per-point allocation or another material solve. The common
`TensorProny::validate_harmonic_mandel_strain` admits the maximum physical norm
over the complete cycle through the scaled two-column Gram eigenvalue.
It includes transverse strain, admits rotating strains that a complex-norm
bound would reject, and catches maxima missed by checking only quadrature
phases. The original full-coordinate solid harmonic owner uses this same
admission. This is a pointwise source-domain check, not a spatial error bound.

Nine tests of the actual numeric strain/projection modules passed against
the real `fs-math` library, covering phase-invariant extrema, large and
subnormal physical scales, finite refusal, complex signs and cancellation.
Five material integration regressions cover full anisotropic stress constraints
in every ordered plane, independently computed complex Poisson strain, complete
cycle limits, transverse-domain refusal and original frequency/DC admission.
Those integration tests remain unexecuted because the real runtime dependency
build is killed before reaching the material crate in this environment.

## Original-unit factors for weak passive constitutive work

`plane_stress_viscosity::factor_passive_work_matrix` factors six work-conjugate
coordinates with fixed-storage compensated LDL updates in the original units.
The three-component resultant and passive shell A/B/D owners share this kernel.
Square roots occur only in output channels, so subtracting rounded square-root
products or normalizing a source matrix cannot distort a resolved weak pivot.
No global positive rank threshold removes an independent small loss channel.
Per-row admission and reconstruction allowances remain `64*EPSILON`; negative
diagonals, coupled zero diagonals and unresolved active work refuse.

Source stress/power and actual CST/DKT strain-work regressions cover nonbinary
weak couplings, both signs of membrane/bending coupling, exact nulls and loss
scales separated by hundreds of decades. The original strict elastic factor
and successful full-rank source stress arithmetic remain unchanged. This
matrix admission cannot recover accuracy absent from the supplied coefficients.


## Source-preserving tabulated flow, axial reduction and small-strain thermal families

`plastic::tabulated::rate::source` binds original Cauchy hardening observations and the explicit canonical Perzyna coefficients from one fixed material context. Linear and power-law kinetics are explicit selections; missing coefficients are not inferred from unrelated viscosity roles. `ResolvedTabulatedOverstressJ2::for_duration` retains the original source and changes only the actual positive physical interval. Immutable knot, slope and stored-energy arrays are shared, and maximum table slope is cached, so duration rebinding neither copies nor scans the measured table.

`UniaxialTabulatedOverstress` reduces the same native bounded return to free-sided axial response with elastic correction E, independent Prager history and exact lateral contraction. It retains `UniaxialJ2State`, including valid accepted viscous overstress. Native tensor storage, incremental potential and energy remain authoritative; yield work and physical overstress work are separated from recoverable hardening and numerical loss. Original source evidence and curve storage survive changing durations. Unsupported measured-strain continuation refuses without projecting or resetting history.

`plastic::tabulated::thermal::AffineTabulatedHardening` shares fixed strain knots, affine absolute-yield ordinates and the exact signed derivative of hardening storage between thermal owners. Its derivative bounds include interior segment extrema. The existing finite Kirchhoff wrapper and resolver keep their public interpretation. A distinct `resolve_thermal_cauchy_tabulated_hardening` entry point admits only original Cauchy observations; no stress relabeling or unit-Jacobian assumption joins these sources.

`ThermoplasticSource::with_tabulated_hardening` supplies `ThermalTabulatedSource` on the existing small-strain caloric/source-family boundary. It retains the original complete source grid and receipts, fixed context and caloric span, requires explicit zero linear H, excludes competing Voce/total-strain tangent laws, and checks agreement of the scalar initial-yield claim across its complete original temperature line. Native expansion, temperature-dependent elasticity and Prager storage, and explicitly selected linear or power-law kinetics use their original sources at the actual interval.

Thermal softening can return a loaded old history outside the new yield surface after that history is admitted at its own old temperature. A crate-private checked-energy path retains native old-history/domain checks and every new-state yield/work gate; public isothermal admission remains strict. Exact table storage derivatives enter mechanical entropy and coefficient work through the existing thermal account. Only physical dissipation and signed reversible exchange feed caloric heat. The existing uniform caloric owner and layered-section family adapter retain the same accepted history through their trials.

Focused regressions cover source retention, tensor/axial parity, cyclic reversal, duration-dependent relaxation, analytical rate returns, measured support endpoints, heated loaded histories, exact thermal storage/entropy, source refusals and coupled section behavior. These are manufactured constitutive fixtures, not experimental alloy or forming calibration. Rust syntax and formatting, bounded source review and independent arithmetic checks were performed; native Cargo compilation and regression execution remain unverified because the shared build environment has no available storage.

## Owned thermal intervals for axial malleability

`ThermalTabulatedSource::at_uniaxial_temperatures` prepares one owned physical interval for the existing free-sided tabulated J2 or Perzyna return. Its immutable endpoint laws, original source choices and exact thermal work preparation are shared across fibers. Accepted state retains permanent strain, accumulated plastic flow and actual temperature. Temperature continuity, original mechanical/kinetic/expansion selections and the source-bound caloric identity must agree before a successor interval is admitted; that caloric identity includes the original Cp/density receipts, not only equal numerical knots.

Thermal expansion shifts all normal strains and the original uniaxial reduction retains transverse plastic flow and Poisson contraction. Both total and mechanical tensor strains remain inside the original small-strain envelope. Heating can lower strength and produce further permanent extension under fixed force; a held member can relax over the actual interval. Unsupported source temperatures, table continuation, chronology or source substitution refuse without advancing history. Preparing, cloning or inspecting an interval does not evolve material state.

The original thermal expansion owner now prepares an owned coefficient account for reuse during constitutive trials. Its shared finish operation performs no source query, allocation or return map. Existing borrowed thermal families use the same arithmetic for elastic, Prager, Voce and table entropy, coefficient work, exact expansion increments and the change in temperature times mechanical entropy. Native physical and numerical dissipation stay with the original material law; numerical loss is never heated. The scalar interval checks its response, lateral freedom and retained state against that native tensor account.

Four G0/G1 regression sources cover heated constant-force knot crossing and retained extension with independent storage/entropy equations; actual-duration linear and power-law holds; original caloric/source identity and atomic refusals; and stress-free expansion with full tensor-domain admission. They also compare the owned account to the existing borrowed thermal family. Rustfmt syntax/format, scoped whitespace and independent scalar arithmetic checks passed. DSR and RCH are unavailable and shared build storage is full; native Cargo compilation and regression execution remain unverified. The declared small-strain and empirical source domains remain binding.

## Original harmonic stress, spring storage and physical heat

`TensorProny::prepare_harmonic_material` borrows the original complete tensor law
and prepares its storage/loss spectrum once at a positive finite frequency inside
the original source band. Allocation-free point observations admit the complete
six-component physical strain cycle and report the full stress phasor, original
equilibrium and Maxwell spring energy, separate Kelvin and Maxwell heat, and an
independent signed stress-work balance. Peak phasors use exp(-i omega t), with
quarter-weight mean spring energy and half-weight mean dashpot heat. DC uses the
original static/time-domain owner and is refused by this sinusoidal observer.

Original Gram factors determine positive storage and heat. Binary-factored
Maxwell weights remain inside the original Gram contractions, their squares and
division by the original relaxation time. A spring strain too small to represent
does not erase representable branch stress or heat; an overflowing intermediate
Gram strain does not reject a representable scaled energy. Both energy and heat
reuse each original Gram contraction. `branch_response` retains source branch
order and reports the original memory-equation spring strain and branch stress;
its rounded strain is not used as the stress/heat intermediate. Work is evaluated
independently from complete stress and strain rate, without forcing closure or
inferring a temperature evolution. Original source conditions and histories are
unchanged. Outputs remain f64 and do not promise relative accuracy at cancellation.

Twelve tests of the actual strain/projection arithmetic passed against the real
`fs-math`, including spring-strain underflow with representable stress/heat and
the opposite intermediate-overflow case. Five native material regressions cover
full anisotropic stress/energy, original ramp-history convergence, transverse
stress constraints, source/domain refusal and extreme relaxation scales. Their
integration execution remains pending because the original runtime dependency
build is blocked by this environment's memory limit.

`TensorPronyHarmonic::prepare_window` extends that same original spectrum to
physical finite intervals. Its allocation-free prepared observer reuses one
analytic temporal Gram and source-clock rotations across all material points.
`observe_window` is the single-point convenience path. Only uniform windows
represent unweighted physical work, heat and endpoint conservation; shaped
detector profiles refuse. Positive frequency, original source band and the
complete-cycle strain envelope remain binding even for a short interval near
a strain node.

The interval account separately evaluates original equilibrium/Maxwell spring
energy at both endpoints, Kelvin and Maxwell heat, and complete stress times
strain rate. Its signed balance is work minus recoverable-energy change minus
physical heat. The storage change is independently integrated from the original
spring energy rate, so a short interval need not lose its change when both
rounded endpoint energies coincide. Centered sine/cosine weights retain the
twice-carrier terms that cycle averaging discards. Source time and duration
remain separate through phase rotation. Original Gram channels, complex Prony
weights, rates, duration and relaxation time retain their binary scales through
the contractions; a rounded-away spring strain does not erase representable
heat. Neither work nor heat is synthesized from a closure residual.
Interval stress work contracts the original equilibrium, Kelvin and Maxwell
coefficients with their binary spectrum, rate and duration before rounding.
A loss entry too small for the cached complex modulus therefore does not erase
representable interval work after a long physical exposure.

This remains one settled sinusoidal particular response at fixed source
conditions. It does not infer startup history, temperature evolution or
cross-frequency interference by adding independent carrier accounts. The
native 3-D solid observer consumes the interval account at every original
constitutive quadrature point with its original geometry and source limits.

Seven interval regressions passed against the selected actual production
constructors, spectrum, source admission, Gram contractions and complete
harmonic response/window owners, alongside eighteen existing numerical tests.
They cover independent real-field Maxwell/Simpson oracles, phase and interval
partition covariance, equal rounded endpoint energies with a retained short
change, and both Kelvin and Maxwell work rescued from underflowed loss entries.
The full Cargo/workspace build remains blocked by the runtime dependency's
memory failure; the separate solid integration tests require that native build.

## Viscoelastic solid acoustic absorption

`elastic::acoustic::viscoelastic::ViscoelasticSlab` evaluates normal-incidence
solid absorption from borrowed original bulk and shear Maxwell laws. The
causal longitudinal modulus `M*=K*+4G*/3` supplies complex propagation and both
physical boundary conditions. Coherent incident pressures on either face
produce independent outgoing powers and positive spatial bulk/shear
strain-work integrals.

Inward-decaying fields, a factored round-trip determinant, and logarithmic
even/odd wave integrals retain opaque layers and weak thin-layer heating
without subtracting nearly equal boundary powers. Preparation bounds the
combined branch count; queries allocate nothing, retain source identity and
frequency support, enforce a conservative all-depth strain ceiling, and poll
cancellation before publication.

Focused numerical regressions cover independent boundary equations and
constitutive dissipation, lossless and opaque limits, thin standing waves,
coherent subdivision and source/support refusal. The model is frozen,
linear, isotropic, planar and normal-incidence. It does not infer measured
authority, oblique P/S coupling, poroelastic frame motion, finite-aperture
scattering or thermal evolution.

`elastic::acoustic::viscoelastic::interface::FluidViscoelasticInterface`
extends these original bulk/shear Maxwell sources to oblique fluid/solid
half-space reflection. Conserved real tangential slowness and outgoing complex
P/SV roots retain refraction, conversion, evanescence and reflection phase.
The prepared three-by-three contact solve enforces normal velocity/traction
and vanishing shear traction. Actual drives enforce the original solid strain
ceiling and explicit fluid small-signal limits; depth queries allocate nothing.

Power and bulk/shear dissipation use the full coherent fields, including the
P/S interference required in an absorbing solid. Volumetric strain is retained
analytically when evanescent normal strains nearly cancel. Entering solid flux
is distinguished from interface heating and from subsequent bulk absorption.
Focused G1/G3 tests compare complex normal impedance, the independent scalar
surface-impedance formula, the original lossless elastic contact across critical
angles, and depth-integrated constitutive work. G0/G4 tests cover source support,
amplitude limits and cancellation. This planar small-strain harmonic model does
not claim roughness, anisotropy, finite-thickness layers, an incident solid wave,
thermal evolution or experimental validation.

## Finite fluid-loaded viscoelastic plates at oblique incidence

`elastic::acoustic::viscoelastic::interface::layer::FluidViscoelasticLayer`
retains the original homogeneous bulk/shear Maxwell laws, density, identity,
frequency/strain support and both exterior fluid sources. A bounded six-by-six
solve prepares coherent two-face pressure scattering with P/SV conversion.
Waves anchored at their respective faces remain bounded in opaque/evanescent
plates; an algebraic cosine/sinc basis preserves the independent solution at
exact internal P or SV critical incidence. A zero-thickness layer reduces to
the original two-fluid contact. Outgoing evanescent or cutoff exit-fluid fields
are retained, while nonpropagating incoming ports are refused.

Drives and depth queries allocate nothing. Original fluid small-signal limits
and a conservative all-depth Mandel strain bound admit the complete coherent
field. Reciprocal power-normalized scattering is exposed when both fluid ports
propagate. Physical face flux, reflection/transmission and the flux deficit are
distinct from local positive bulk/shear constitutive heat; interfaces receive
no invented heat. Optional allocation-free depth work uses phase-resolved
positive Gauss quadrature, explicit tolerance/work limits and 8/16-point error
estimates, refusing unresolved ranges. This retains weak standing-wave heat
below boundary-power subtraction roundoff. The estimate is not an interval
certificate, and quadrature cost grows with resolved complex phase.

Eight focused regressions cover independent complex normal Fabry–Pérot
scattering, source reuse, distinct-fluid reciprocity, lossless resonances, exact
P/S critical roots and basis crossover, evanescent tunnelling, opaque half-space
limits, coherent bulk/shear work and depth subdivision, weak thin odd-wave
heating, cutoff/evanescent exterior admission, source limits and cancellation.
This is
a fixed, isotropic, small-strain, laterally infinite single plate. Exact undriven
guided-mode singularities may refuse. No anisotropy, roughness, compliant contact,
arbitrary oblique stacks, finite aperture, thermal evolution or experimental
material validation is inferred.

## Calibrated temperature snapshots of source-resolved plane-stress viscosity

`TemperatureShift::shift_plane_stress_viscosity` applies the existing calibrated
WLF/Arrhenius kinetics to the original reduced Kelvin source: `V_T = a_T V_ref`.
The source's exact reference coordinate must match the declared reference
temperature and carry absolute/static semantics, or the legacy `T` convention.
All selected coefficients must retain only admitted fixed external conditions;
an existing dependence on that temperature axis or a second temperature owner
refuses. Other held moisture, pressure or composition context remains intact.

The original source receipts, selected properties and material plane stay
unchanged. Original matrix entries scale by `a_T`; cached passive Gram factors
scale by its square root without refactorization or loss of resolved weak
channels. Exact zeros/signs survive; nonzero underflow and overflow refuse.
The original tensor-shift frequency policy transports support to
`[omega_lo/a_T, omega_hi/a_T]`, intersecting its upper endpoint with finite
representable frequencies. Reference-temperature evaluation clones identically.
The returned `ShiftedViscoelastic` retains the actual temperature/calibration.

Each low-level call requires the original reference law, as for `shift_tensor`.
This is a prepared fixed-temperature source interpretation, not live thermal
history, an elastic modulus shift, a caloric balance or independent authority
for time-temperature superposition. Focused tests use genuine resolved scalar
sources to check stress/heat, weak passive channels, held moisture, source-axis
refusals, immutable reference data and support/representability bounds.

## Passive Mandel tensor work factors

Full three-dimensional `TensorProny` branches and optional Kelvin--Voigt
viscosity use `plane_stress_viscosity::factor_passive_work_matrix` for their
passive Mandel work factors. Admission uses the original dimensional
coefficients and the shared compensated Schur arithmetic; it does not square a
global Frobenius norm, discard positive loss relative to the strongest axis, or
introduce a damping floor. Resolved rank-deficient null directions remain
undamped. Nonfinite, asymmetric, active, or unrepresentable data are refused
under the shared factor's local admission and reconstruction bounds.

The original constitutive tensors, relaxation times, branch order, and
physical Mandel internal strains remain unchanged. Exact affine-ramp stress
and mean stress still use the source tensors; branch storage and independently
integrated nonnegative heat use their admitted work factors. The equilibrium
positive-definite factor and elastic admission path are unchanged.

Five native regressions demonstrate and now pass against the original defects,
including threefold Kelvin heat, phantom null-direction work, passive scale
overflow/underflow and hidden active cross terms. Their independent oracle
contracts original tensors and integrates physical heat by quadrature. The
shared replacement factor also has passing source and shell coverage. The
full-tensor native gate covers all five original regressions unchanged.

## Edited Prony branch source coherence

`PronyBranch::c_pa` remains editable before law admission. Both `TensorProny::new`
and `new_mandel` reuse cached Gram rows only when the existing factor-defined
source reconstructs exactly to the current coefficients. An unchanged explicit
Gram source retains its original row coordinates and weak channels. Edited
coefficients instead pass through the shared passive work factor; an old cache
cannot restore removed damping or conceal active, asymmetric or nonfinite data.

Each admitted branch supplies the same rows and transposed factor to the
physical-strain history, factor-coordinate history and strain projection.
Rotation transports those admitted channels. Source consistency is checked
once during construction, with no additional factorization or coefficient
comparison in stepping. Relaxation times, branch order and zero-branch history
slots remain intact.

Five focused regressions fail with the old cache behavior and pass after this
repair, including original-source stress, stored energy, integrated heat,
projection, rotation, harmonic response, zeroing and invalid-source refusal.
The existing 256-by-512 random-history gate now supplies the original generated
Gram rows: rounding their rank-deficient product into a matrix can introduce
an active direction. Its random draws and physical assertions are unchanged;
matrix-only passivity admission is not weakened to repair synthetic input.

## Complete cyclic history in owned thermal axial intervals

`ThermalChabocheSource::at_uniaxial_temperatures` prepares source-bound,
shared axial intervals retaining one through eight independent Chaboche
branches, Prager and Voce history, thermal expansion and actual physical time.
The original scalar return also supplies linear or power-law Perzyna response
through `UniaxialChabocheOverstress`; its tangent differentiates the same joint
branch/flow equation. No effective single backstress replaces the original
ordered branches.

Temperature changes preserve each strain-like branch coordinate X/C before
the native mechanical update. Original source and caloric identities, accepted
temperature and branch roles govern successor admission. Source queries,
thermal integration and immutable preparation occur outside constitutive
trials. The original tensor owner independently checks scalar stress, transverse
freedom and every history component in stress units. The audited native state
and shared thermal account retain branch storage/entropy, coefficient work,
physical yield/recovery/viscous heat and numerical loss separately. Numerical
loss never heats the material. Native storage inspection performs no return or
relaxation, including on accepted viscous states.

The borrowed and owned thermal paths share branch thermodynamics. A restricted
old-history admission seam permits thermal softening only after the old history
is checked at its own temperature, including explicitly absent branches; public
isothermal admission retains its prior behavior. Five focused regression
definitions cover thermal cycling/reversal, independent rate/tangent and heat
equations, unloaded branch transport, original-source/clock/domain refusal, and
the zero-branch limit. Independent scalar fixture arithmetic passed. The
execution host became unavailable before edits; no Rust formatting, compilation
or native regression execution was performed. Small strain, globally affine
positive-or-absent branch moduli, fixed recovery coefficients and original
source domains remain binding; finite cyclic forming, static recovery and
experimental calibration are not supplied by this extension.

## Coherent absorption in bonded viscoelastic solid layers

`elastic::acoustic::viscoelastic::layers::ViscoelasticLayers` retains each
layer's original bulk/shear Maxwell laws, density, thickness, source identity
and frequency/strain support. Balanced real reference impedances are only wave
coordinates at the exact bonded interfaces; they add no fluid gap, material,
propagation or thermal receiver. A compensated suffix scattering recursion
and cached internal input maps enforce coherent traction/velocity continuity
with O(n) preparation and drive work. Repeated drives and local observations
allocate nothing.

Each physical layer uses the original positive bulk/shear strain-work integrals,
with independent exterior power accounting. The shared slab kernel now cancels
common impedance before applying transmission phase and evaluates even/odd wave
amplitudes directly, preserving weak thin-layer heating under subdivision.
Material sources remain immutable; each distinct exterior source is queried
once. Layer/branch budgets, conservative all-depth strain limits and
cancellation apply before publication.

Focused numerical checks cover an independent global physical boundary system,
constitutive heat quadrature, unequal-material interface continuity,
reference-coordinate invariance, subdivision, opaque shielding, weak standing-
wave heat, reciprocity and coherent drive scaling. The scope remains frozen
linear isotropic planar normal-incidence solids with perfect bonds. Oblique
conversion, interface compliance/slip, poroelastic frames, thermal evolution
and measured-material qualification are not inferred; arithmetic residuals
are not physical error certificates.

## Original normal-stress-free harmonic material with retained shear

`TensorProny::prepare_normal_stress_harmonic` borrows the complete original
anisotropic Kelvin/Maxwell law and eliminates only the transverse normal stress.
Its five engineering coordinates retain both transverse shear tractions and
every in-plane/shear cross term. One original complex constraint supplies the
modulus, physical exp(-i omega t) full Mandel strain recovery and analytic
frequency tangent. Storage/loss use constrained original positive work; the
tangent uses the same cached constraint in the complex-symmetric
`T^T C_omega T` contraction. Shared scaled complex division, binary-factored
contractions and passive-work admission retain their numerical ownership.

The source frequency band and complete six-component whole-cycle strain ceiling
remain binding. Prepared recovery allocates nothing and rejects imaginary DC
strains. This is a local constitutive reduction: it supplies no shear correction,
through-thickness trial profile, finite-element kinematics, laminate traction
equilibration or hereditary time integration. Existing full plane stress keeps
its three retained coordinates and original 3x3 elimination arithmetic.

Six native integration regressions and one scalar-owner regression are authored
for independent one-normal Schur/traction/work identities, all ordered planes,
physical harmonic phase, all shear/cross terms, analytic and finite-difference
tangents, source-domain admission, isotropic shear and extreme unrelated loss
scales. Native compilation and execution remain unverified because the shared
workspace is unavailable; no substitute runtime or dependency stub was used.

### Passive midpoint coupling of inertial polarization

`field_storage::dielectric::lorentz` now supplies source-owned midpoint
condensation and recovery for original polarization/current vectors. The
endpoint field coefficient is `epsilon_infinity + sum delta*(h*omega0)^2 /
(1+h*gamma+(h*omega0)^2)`, with `h=dt/2`; the history retains both P0 and J0.
Recovery keeps each original pole in its authored order and separately reports
instantaneous-field, restoring and kinetic storage, actual endpoint increments,
midpoint electric work, and nonnegative frictional heat. Heat is evaluated from
the original midpoint current and is never a balance residual. Rounded P/J
publication is reflected in the reported energy residual. One-, two- and
three-dimensional states are supported; an empty bank is the original lossless
instantaneous material. Invalid interval, vector shape, nonfinite input, and
unrepresentable recovered quantities refuse without mutating the inputs.

`PreparedLorentzMaterialMidpoint` caches one fixed interval's original pole
coefficients and positive condensed permittivity for repeated spatial cells.
Warm histories allocate nothing and neither histories nor recovery re-prepare
the coefficients or query the material. Only the recovered bounded state vector
allocates. The prepared cache does not clone or replace the complete source
record; the field owner retains that record and its identity. Convenience methods
reuse the same recovery kernel. Dimensional products retain binary exponents
and use the existing quad-double math owner; no pole fit or second constitutive
model is introduced. This is second-order passive implicit midpoint, not an
exact propagator or an L-stable method. Timestep refinement must resolve
resonance phase, damping and field variation. Spatial convergence, global
source/boundary work admission, thermal parameter work and full-wave Maxwell
regime validity remain with their existing owners.

Twenty-one native tests of the actual oscillator implementation pass, including
six new midpoint tests of the original equations/condensation, lossless phase
and energy, damped passivity, continuous-solution convergence, representable
extreme scales and refusal/replay. Three source-bank integration regressions
are authored but unexecuted without the complete material dependency closure;
they cover cached/convenience equality, original multipole work and energy,
source retention, the empty bank, and refusals.

## Complete harmonic laminate work in an explicit physical shear trial

`laminate::harmonic::HarmonicLaminate` binds immutable actual full Elastic/Prony
layers to a separately supplied original `LaminateShearRecovery`. Layer geometry,
reference faces and material angles must agree; actual I0/I1/I2 and source
receipts remain with the actual stack. No Prony surrogate, new shear correction
or reference profile is manufactured. Reduced-only actual layers refuse.

The prepared section integrates all 8x8 membrane/bending/transverse-shear terms,
including anisotropic cross blocks and signed offset coupling, using each
source's complete five-component normal-stress-only reduction. Original cached
shear Bernstein controls and their original material rotation supply one fixed
physical trial field. Positive three-point thickness Gauss integrates its
degree-four energy exactly in real arithmetic. Existing binary-factored work
contractions own accumulation; operators retain complete engineering work
coordinates and exp(-i omega t) physical phase. The final two resultants are
generalized forces conjugate to the declared shear trial coordinates:
`Qs = integral(P_s^T tau_actual dz)`. Outside the reference equilibrated source,
these need not equal the unweighted actual shear-force integral; consumers must
retain the declared trial work map. Prepared exact frequency slopes
reuse original normal constraints, branch times and Kelvin loss with the
declared geometry/profile held fixed.

Allocation-free resultant/recovery observations preserve every original
six-component strain ceiling. Three Bernstein controls per layer provide a
conservative complete-thickness and whole-cycle admission, including shear peaks
hidden by the exterior faces. Caller budgets bound retained source/frequency
payload and branch-dependent work; cancellation returns no partial result.
The original profile accessor also exposes material-axis Bernstein controls,
using the existing recovery's cached physical rotation.

This is a fixed-profile Galerkin section. Only the reference elastic source owns
its two-probe traction equilibration. General actual anisotropy or frequency
dispersion may produce different exterior/interface tractions, which remain
visible rather than being discarded. Frequency-dependent warping, native
enriched MITC assembly, its original physical mass/gauge coupling, time
integration and thermal evolution remain separate work.

Six authored native regressions compare original elastic A/B/D/H, actual density
moments, independent five-point full-tensor stress integration, retained cross
work, phase, source rotations, exact frequency slopes, interior strain admission,
source incompatibility, budgets and cancellation. These tests have not been
compiled or executed because the workspace is unavailable; source review and
exact remote blob readback do not substitute for native execution.

### Source-owned conduction with inertial dielectric resonance

`field_storage::dielectric::lorentz::leaky::ResolvedLeakyLorentzMaterial`
binds the original Lorentz material to an independently selected scalar DC
conductivity. Adoption moves the existing complete source snapshot without
re-resolution, cloning it, or dropping unrelated multiphysics roles. The
one-shot resolver uses the same bounded Lorentz requirement builder and adds
the separate conductivity role. Angular natural frequencies retain their
original quantity semantics; pole damping and conductivity independently admit
exact zero. Both Debye and Lorentz composition use the same DC source gate.
Missing or aliased roles, wrong quantities, negative/nonfinite values, and
unowned frequency/field/strain/stress restrictions refuse.

The consuming decomposition supplies the unchanged original dielectric, exact
conductivity property name and source value to the existing spatial owner.
Polarization current remains dP/dt, with separate restoring/kinetic storage and
physical pole-friction heat. Conduction contributes sigma E and sigma |E|^2,
without storage; the selected DC source must exclude the loss already assigned
to resonant poles. The field owner enforces total-current continuity and owns
conduction integration, boundary work and heat deposition. No duplicate local
oscillator, integrator, coefficient fit or material-refresh path is introduced.
This remains fixed-condition isotropic Ohmic transport, without carrier inertia,
anisotropic conductivity, breakdown, induction or measured qualification.

Three focused G0/G1 source regressions are authored for complete identity and
role retention, independent zero/finite-conductivity and zero/finite-pole-loss
limits, unchanged original polarization/current recovery, and source refusals.
Formatting and a focused source review pass. These integration tests are not
claimed executed without the complete material dependency closure.


### Timed static recovery of cyclic backstress

`ChabocheStaticRecoveryJ2<N>` explicitly adds `dX_i/dt = -r_i X_i`
with independent rates in s^-1 and the actual positive physical interval,
alone or with the original linear/power-law Perzyna kinetics. Fixed preparation
divides C, gamma and incoming X by `1+r*dt`, then reuses the original scalar
consistency, direction drift and exact tensor/axial tangents for simultaneous
recovery and plastic flow. Original C remains authoritative for stored energy;
the original work account separates physical dynamic/static recovery from
backward-Euler numerical loss. No response performs a new source query.

Incoming history and accepted candidates retain original bounds. Rate/duration
overflow and unresolved active coefficients refuse. Branches remain independent
and proposals never mutate accepted history. Three authored G0/G1 regressions
cover independent dwell recurrence and continuous-time refinement, simultaneous
nonproportional flow/tangent/energy, zero-rate and absent-branch reductions, and
refusal. Rust parsing/formatting passed; native compilation and tests remain
unverified. This is linear static recovery, not nonlinear recovery, isotropic
annealing, recrystallization, finite kinematic plasticity or empirical calibration.

## Depth-resolved viscoelastic acoustic absorption

`DrivenViscoelasticLayers::absorption_between` integrates the original positive
bulk/shear strain work over a physical depth interval without adding a material
interface or subtracting cumulative heats. Stable even/odd restriction preserves
weak thin standing waves; independently scaled inward waves and logarithmic
attenuation preserve deep opaque heat. Full intervals retain the original
whole-layer arithmetic, zero-width intervals yield zero, and invalid depths
or cancellation refuse. Each query takes constant work and allocates nothing.

Four native regressions compare independent constitutive quadrature, arbitrary
partition closure, thin standing waves and square-underflow opaque tails.
The query retains the source identity and the original conservative whole-layer
strain bound. It does not itself advance a thermal history or infer changes to
the frozen acoustic material laws.


### Original-source thermal dwell recovery across cyclic consumers

`ThermalChabocheSource::with_static_recovery` explicitly selects one complete
`-r X` coefficient in s^-1 for each ordered branch. Original same-card receipts
and nonnegative full caloric-span support are pinned: no aliased roles, source
extrapolation, inferred rates or duration-free evolution. Temperature-dependent
kinetic rates add no stored energy or entropy. Native thermal history transport
precedes the simultaneous implicit recovery/flow return; original branch moduli
still own storage, coefficient work and entropy.

Tensor, owned axial and uniform caloric paths preserve actual duration, all
material history and physical recovery heat exactly once. Existing thermal
fiber/frame and continuum enthalpy consumers inherit this source family through
the same APIs. Family identity also retains the ordered static-rate names.
Selection on a uniform caloric point must precede loading and physical time.

Three synthetic G0/G1 regression definitions cover hot dwell/reversal with
independent branch/rate/storage equations and native tensor/axial parity, full
source/quantity/family/time refusals, and uniform adiabatic recovery heat with
unchanged plastic accumulation and cancellation/retry. Rustfmt checks passed;
compilation and native tests remain unverified. This is linear static recovery;
nonlinear recovery, isotropic annealing, recrystallization, finite cyclic
deformation and empirical calibration remain outside this implementation.

## Direct reciprocal exchange between finite phase bodies

`phase::exchange::pair::exchange` and `exchange_with_checkpoint` integrate an
isolated fixed-conductance trajectory between two `FinitePhaseBody` values.
Each retains its own original chart, state identity, fixed mass, material
context and enthalpy datum. One signed heat Q debits one body and credits the
other; requested exposure is G times duration in J/K. The defining integral
uses the actual evolving hot/cold temperature difference on both original
enthalpy paths.

Authored affine T(h), constant original Cp and isothermal latent intervals
advance analytically. Varying affine Cp(T) retains the original sensible
inverse and uses source-split Simpson quadrature with analytic truncation
bounds and explicit floating-point allowances, followed by bounded scalar
inversion. Both heat and exposure accuracy gates apply. Incomplete integral
lower bounds can reject a trial but can never accept a solution. Original
source intervals are prepared once; total prepared intervals, root iterations
and total nonlinear panel visits have explicit caps. Cancellation is polled
through preparation, integration and final entropy evaluation.

The proposed states retain hot/cold ordering. A source endpoint reached at
positive temperature gap has finite available exposure; larger requests
refuse. Thermal equilibrium is approached asymptotically and is never crossed.
Original latent heat is consumed once, without an apparent heat capacity,
additional latent source or prescribed phase-front time. Equal initial
temperatures and zero exposure preserve both accepted states.

`FinitePhaseExchange` reports retained relative heat, both represented
absolute-enthalpy changes and their individual residuals. Trajectory error
does not absorb these storage residuals: callers must enforce their physical
energy budgets. The signed entropy account integrates the original charts
between the represented states, including latent entropy, and exposes its
numerical allowance rather than clipping a negative result. The numerical
bounds are not interval-arithmetic certificates or source uncertainties.

This owner conserves the charts' declared enthalpy account under their original
pressure/work convention. It introduces no mass transfer, external source,
changing conductance, spatial approximation, pressure-work source or claim of
constant-volume internal-energy conservation. Native conduction admission
and geometry remain with the lumped adapter. Focused regressions exercise
analytic and nonlinear calorimetry, phase boundaries, reversal, subdivision,
source support, representation budgets, bounded work and cancellation.
Pinned Rust parsing was available; native Cargo compilation and test execution
were unavailable in this resource-exhausted session.

## Spatial acoustic absorption on native thermal histories

`elastic::acoustic::viscoelastic::thermal::ViscoelasticHeating` binds a complete
nonoverlapping depth partition of the original physical solid layers to active
native thermal cells. Every bin retains its original material identity and
positive bulk/shear work integral. Cell volume and mass must agree with the
authored area, depth and original density. No acoustic interfaces, uniform
layer-heating approximation, midpoint source samples or fitted loss law are
introduced to obtain thermal resolution.

The original enthalpy owner performs storage, latent heat, conduction, boundary
events and acceptance. Exact prepared cycle-mean powers enter its interval source
hook once; the accepted native duration owns each thermal credit. Independent
incident, outgoing and constitutive absorbed energies remain separate from the
actual cell heat credit. Positive depth powers are never renormalized to force
boundary closure. Cumulative arithmetic discrepancies use an explicit joule
allowance. A failed endpoint, changed mass/geometry/caloric binding, exhausted
work allowance or cancellation publishes no acoustic or thermal state.

A bounded optional affine power envelope scales a fixed coherent two-face
carrier. Its exact mean intensity, jumps and dark gaps use native transition
clipping; no acoustic wave solve or acoustic material/exterior query is repeated
during thermal advances. Preparation sorts bins once and evaluates each physical
layer/bin. Repeated source work makes bounded passes over bins with the native
ordered-cell and source-map lookups. Immutable prepared data is
shared across candidate clones. Restart preserves cumulative acoustic energy.

The laws, density, geometry and carrier frequency remain frozen within explicitly
declared temperature support. This is quasi-steady cycle-mean heating: envelope
changes must be slow relative to carrier, propagation and ring-up scales.
Temperature/phase-dependent acoustic response, carrier transients, changing
mechanical stored energy, streaming, deformation and measured material authority
are not inferred. Spatial and temporal thermal refinement retain the native
cell-uniformity and supplied-conductance assumptions.

Focused native regressions cover heterogeneous depth heating, scalar pulse
fluence, boundary/contact events and restart, original latent storage, a
three-rung conduction convergence ladder against an independent two-cell solution,
and partition/inventory/source-domain/work/cancellation refusals.


## Original NASA CHARMS Corning 7980 refraction

`silica_refraction::Corning7980Silica` resolves the original Leviton–Frey Table-3 temperature-dependent Sellmeier fit for the measured Corning 7980 synthetic fused silica specimens. It admits only 30–300 K and 400–2600 vacuum nm, retains source and frozen-temperature identities, and exposes the existing three-term oscillator convention with squared resonance wavelengths. The material owner has no renderer dependency. Preparation and queries use fixed storage without allocation.

Phase index, group index, wavelength derivative, and first/second temperature derivatives all come from the same original law, including the temperature dependence of resonance wavelengths. These are derivatives of the global index fit, not the separately processed Table-6 thermo-optic estimates; their discrepancy near the temperature endpoints is explicit. Source applicability and condition-dependent measurement uncertainty remain unchanged. Absorption, scattering, pressure response, stress optics, thermal expansion, and other silica grades are not inferred.

The NASA primary source and public-use record are retained in the module. Four focused G0/G1/G3 tests compare 21 independently printed indices, independent high-precision numerical derivative controls, temperature/wavelength finite differences, frozen coefficients, identities, and domain refusals. These tests and strict Clippy passed in a narrow Cargo build of the exact source and original fs-math/fs-blake3 dependencies; no full-workspace or experimental validation is claimed.

## Finite primary creep with independent resistance recovery

`plastic::finite::creep::recovery::FiniteRecoveringNortonCreep` extends the
original authored finite Norton interval with an explicit nonnegative recovery
rate at the same fixed temperature. Its complete `FiniteRecoveringCreepState`
retains the original isochoric creep metric and cumulative `native.alpha`, plus
independent `resistance_strain`. Backward Euler solves
`r_new = (r_old + delta_alpha)/(1+k*dt)` alongside the same primary flow law.
The private prepared scale `a*(1+k*dt)` and old resistance enter the shared
logarithmic primary root and inverse-rate potential; permanent-history bounds,
finite metric reconstruction and analytic `dP/dF` stay with the native owner.
Each trial uses one scalar constitutive return, and line searches omit tangent
assembly. Original coefficients remain inspectable without exposing the
duration-transformed scale as a material parameter.

Recovery neither resets finite permanent deformation nor reduces cumulative
creep. Resistance stores no energy in this phenomenological model; physical
creep heat remains `q*delta_alpha`, with q the equivalent Kirchhoff stress.
The original standalone finite `plastic_*` field names retain their creep
meaning. `FiniteJ2Response` accepts a complete typed history while defaulting
to its original state, so recovery candidates carry one authoritative history.
Inspection performs no flow or recovery. Invalid, unreachable or unrepresentable
states refuse without modifying accepted history. Zero recovery with r=alpha
reduces exactly to the original primary law, and disabled primary resistance
preserves the original secondary mechanical return.

`tests/plastic_finite_creep_recovery.rs` covers G0 positive heat, history and
tangent invariants; G1 independent quadratic returns, cubic kinetic potentials,
stationary constant-stress flow and unloaded dwell; and G3 rotations, hydrostatic
stretch and exact disabled reductions. Reload after a dwell resumes faster creep
without recovering the finite permanent set. These are authored numerical
verification cases. Automatic temperature shifting, source-card selection,
simultaneous yielding, thermal feedback, tertiary damage, rupture and measured
service-life qualification are not supplied by this native material-point owner.

## Original density through the shared thermal interface

`PcmComposite`, `CrystallizationCaloricBranch`, `FreezingWaterThermalLaw`,
`IsobaricSuperfluidLaw` and `GasEnthalpyLaw` expose their original density at the
requested specific enthalpy through `ThermalLaw::bulk_density_kg_m3`. These
adapters retain the original domain admission, phase mass fractions, additive
phase/constituent volumes, integrated thermal expansion and the gas mixture's
retained pressure and composition. They do not substitute reference mesh
density, average densities by mass fraction or fit another equation of state.
The existing opt-in material-density conduction layers can therefore derive
their normal distances from conserved mass for these material laws. Density
evaluation itself changes neither inventory nor stored geometry.

Focused G1 checks cover the additive-volume PCM and ice/water plateaus,
crystallization-branch expansion, the analytic superfluid density integral,
the isobaric mixture gas density and bracket, and two explicit PCM conduction
steps with updated resistance and conserved mass and energy. Conductivity,
finite-step accuracy, pressure work and the geometric constraints remain the
responsibilities of their existing consumers.

## Sourced temperature-dependent thermo-creep elasticity

`ExpandingCreepSource::from_norton` and `from_norton_with_strain_hardening` admit independently affine original `E(T)`, `nu(T)` and `alpha(T)` through the existing thermoelastic potential. `from_norton_j2` and `from_norton_j2_with_strain_hardening` additionally retain simultaneous J2 yielding and independently sourced yield softening. Norton reference/primary coefficients remain constant over the complete selected temperature span. J2 energetic hardening must be a constant, directly sourced plastic-strain modulus; a post-yield tangent basis is refused because its converted energetic modulus changes with E. Unsupported interior knots, missing expansion data and changed source roles are refused.

The native mechanical return uses actual endpoint stiffness and retains all accepted irreversible histories. The thermal owner corrects accepted free-energy storage by `delta(Psi_T) = e_0:[C(T_1)-C(T_0)]:e_0/2` and transfers reversible heat `sigma_1:I delta(f) - delta(Psi_T) - delta(Ts)`. Physical creep/plastic heat, signed reciprocal heat, entropy storage and numerical dissipation remain separate. The original potential owns the full-span thermal-capacity bound. Free-sided reductions reconstruct accepted transverse strain with its retained temperature's Poisson ratio and share the same resolved source with the tensor audit.

Existing bare source constructors retain their constant-storage requirement. A source prepared for varying elasticity refuses bare `at_temperature` continuation; its `at_thermoelastic_temperature` decomposition is for the complete thermal owner. These laws remain bounded small-strain models. The separate recovering-source families retain their existing fixed-elasticity construction. `tests/creep_thermal_softening.rs` checks analytic thermal cycles, mixed-shear relaxation and tangent, nonvirgin axial/J2 continuation and source refusals. These manufactured cases do not constitute measured material calibration.

## Typed finite inelastic points for spatial continuation

`plastic::finite::FiniteInelasticMaterial` retains the complete associated
history of native J2, Perzyna, Norton, simultaneous Norton/J2 and independent
primary-recovery laws. Existing source-admitted `FiniteJ2Material` wrappers keep
their original native owner through a blanket adapter. The protocol never
unwraps a kinetic law to a rate-independent inspection backbone. Its optional
duration is the original physical constitutive interval; spatial owners check
their common material clock and, with inertia, agreement with the physical step.
Frozen-history line searches and inspections consume no additional interval.

Stress, analytic tangent, original-volume storage, integrated flow potential
and physical inelastic heat remain the native response, including the original
sole-creep field convention. The flow potential is used for variational audits;
only physical dissipation belongs in a heat account. This seam introduces no
thermal shift or source calibration. `tests/plastic_finite_inelastic.rs` checks
generic-consumer relaxation, stress-free recovery and analytic reload, plus
preserved equilibrium/overstress checks and distinct kinetic heat/potential.

## Recovering thermo-creep with sourced elastic softening

`ExpandingCreepSource::from_recovering_norton` and `from_recovering_norton_j2` bind independently affine original E(T), nu(T) and alpha(T) to the existing primary-resistance recovery laws and thermal potential. The thirteen-/fifteen-claim bundles retain independently sourced Norton and recovery Arrhenius parameters, their complete temperature support and selected claim identities. The original constructors retain fixed elastic storage; a source prepared for varying elasticity refuses bare continuation, including access through the joint source's recovery/Norton children. Joint energetic hardening must remain a constant, directly sourced plastic-strain modulus.

Endpoint construction rebuilds the original native law with actual E/nu while retaining permanent creep, plastic strain, both accumulated paths and independent resistance. Recovery still obeys `r1=(r0+delta_alpha_c)/(1+k(T1) dt)` and owns no recoverable energy or heat. The original thermal owner supplies coefficient-work correction, entropy storage, signed reciprocal heat and the full-span capacity bound; physical creep/plastic heating and numerical loss stay separate. Free-sided reductions reconstruct accepted lateral strain with nu at its retained temperature and share the complete resolved law with their tensor energy audit.

The existing thermal-fiber test seam specifies G1 heating, loaded relaxation, stress-free hot dwell, cooling and reloading against independent axial flow/tangent and Lame-energy equations for both recovery families. G0 cases specify exact zero-recovery reduction under nonproportional loading and refusal of changing recovery reference claims, nonsmooth elasticity, unsupported hardening, bare storage bypass and cancellation. The original fs-couple tetrahedral/enthalpy fixture additionally exercises loaded elastic softening, independent temperature-dependent recovery, hot dwell and reload through actual nodal equilibrium and caloric feedback. Its independent reference retains the original elastic trace/deviator, coefficient work and fixed-history entropy. These are bounded small-strain phenomenological source models, not experimental calibration, tertiary damage or rupture predictions. Native execution remains unverified after dependency compilation was killed.

## Sourced liquid-water absorption at 25 degrees Celsius

`fs-material::water_absorption` retains all 169 original Hale–Querry (1973) Table I imaginary-index knots, using the CC0 refractiveindex.info transcription at Git blob `114173ed00461bc5b4910ffd6711aa6759ec86c5`. Source support is 200–200000 vacuum nm at the single nominal temperature 298.15 K. The historical compilation is not a modern visible-accuracy claim or a temperature/pressure-dependent absorption model. The allocation-free evaluator shares the original table, interpolates dimensionless k, and then computes intensity absorption `4*pi*k/lambda_vacuum_m`; no RGB conversion, knot reduction or extrapolation is applied. Its 32-byte identity binds the source, every knot, units, temperature, band and interpolation. Phase/group dispersion, scattering and strongly absorbing interface optics require their separate owners.

Validation: all 169 k values and exact micrometre-to-nanometre conversions were compared with the pinned source transcription. Three focused tests of the original source module passed via the existing small source facade, including independent SI controls, interpolation order and invalid/support-boundary refusals; strict Clippy passed. This records source-module execution, not a full workspace build.

### Original ordinary-water wavelength curvature

The IAPWS R9-97 refraction owner also reports
`WaterRefractiveSample::wavelength_second_derivative_per_nm2`, obtained by
analytically differentiating the original rational wavelength terms and the
same nonlinear Lorentz-Lorenz phase law. Preparation retains the wavelength
coefficients for allocation-free T/rho sampling. This derivative holds the
original temperature and density fixed, preserves the exact vacuum limit and
the original 200--1100 nm domain, and introduces no EOS or absorption model.
Native phase and group wavelength differences verify the added curvature;
the release's original uncertainty and applicability limits remain unchanged.

### Freezing-water signed volume response and density-dependent slab layers

`FreezingWaterState::specific_volume_enthalpy_derivative_m3_j` returns the
original isobaric dv/dh in m³/J. Sensible ice and liquid use the retained source
Gibbs derivative (dv/dT)/Cp. Fusion uses (v_liquid-v_ice)/L, including both exact
latent endpoints under the existing inclusive plateau convention. Melting
contraction and liquid water's negative expansion below its density maximum
retain their signs even when dT/dh is zero. The state query allocates nothing,
performs no inverse or source reevaluation, and refuses nonrepresentable
tangents. The original implicit conductivity binding supplies this signed
tangent to the shared nonlinear geometry owner.

`FreezingSlab::with_material_density_layers` opts the existing explicit,
implicit and adaptive slab paths into their original material-density layer
geometry. Each cell keeps its original mass and fixed footprint area; current
specific volume determines its normal conductive half-distances. The original
reference coordinates remain inventory labels. Binding requires the native
initial mass/density/reference-volume admission and must occur while those
quantities agree. Accepted continuation retains the selected geometry, and
refusal or cancellation leaves the original slab unchanged.

Conductivity, caloric enthalpy, source domains, imposed pressure, fixed face
areas and contact resistances retain their original owners. The independent
H and U+pV accounts use the actual signed material volume change. This is
declared normal layer expansion at fixed footprint, not lateral mechanics,
confinement stress, advection or an additional pressure law.

Eight complete native freezing-core tests passed, including source-derived
sensible signs, latent contraction and endpoint tangents. The two new native
slab regressions also passed: explicit geometry refresh/cancellation and an
independent two-cell series-mixture backward-Euler phase-volume quadratic,
including fixed-geometry discrimination, energy closure and continuation.
All 32 freezing test bodies compiled in the bounded integration build; only
the two slab regressions ran in that build. It used complete fs-thermal source
at cdb00fe5f4eb2db4eeb21ab9a0cae934747916b7 with native linear/nonlinear numerical
modules and nightly-2026-08-31, excluding only unrelated runtime wrappers.
No thermal solver or material algorithm was substituted. Full workspace Cargo
validation remains unrun.

## Explicit-time anisotropic Hill overstress

`plastic::hill::rate::LinearOverstressHill48` supplies the associative n=1 kinetic law `alpha_dot=max(q(sigma-X)-Y-R(alpha),0)/eta`, with independently selected positive viscosity [Pa s] and actual step duration [s]. Alpha remains the original Hill work-conjugate multiplier; coefficients fitted against a differently normalized J2 strain cannot be substituted implicitly. Aligned and fixed-oriented constructors retain the complete original linear/Voce/Prager law and local material history.

The original five-mode scalar return and exact tensor tangent include `eta/dt` only in the incremental radius. Original yield and energetic hardening coefficients remain unchanged. Accepted viscous histories may have positive static overstress; duration-free history admission checks their finite, isochoric, reachable small-strain state without a static-yield projection. Repeated proposals borrow the same accepted history and consume no additional interval. A further accepted hold interval can relax stress and advance plastic history.

Physical dissipation is `Y*delta_alpha + (eta/dt)*delta_alpha^2`. Original elastic, isotropic and Prager storage and backward-Euler losses retain their existing owners. The endpoint-frozen thermal seam preserves old-coefficient history admission and original material axes; coefficient work and entropy belong to the thermal owner. Native tests specify independent anisotropic shear relaxation and endpoint slopes, J2 correspondence, rotated nonproportional tangents, plane-stress holds and stress-unit covariance. This bounded law does not supply sub-yield creep, a nonlinear rate exponent, dynamic backstress recovery, evolving texture, finite-strain forming or measured-alloy calibration.

### Source-bound linear overstress for thermal Hill solids and sections

`ThermalHillSource::with_linear_overstress` explicitly selects
`hill_overstress_viscosity` in Pa s for Hill's work-conjugate alpha. It pins
the original claim and admits every source-line crossing over the original
fully solid caloric chart. The viscosity must remain positive and finite;
hidden rate/frequency axes, partial temperature support, substitute fluid
viscosity, and implicit J2 normalizations are refused. Dissipative viscosity
may vary nonlinearly with temperature because it contributes no stored
energy, entropy, or frozen-history heat capacity.

Each `resolve_section_step` prepares the native oriented linear-overstress
law once at its actual positive physical duration. The same old enthalpy,
original material-frame history and interval feed the mechanical solve and
the independently checked work account. Accepted viscous histories may
remain outside the static yield surface; validation inspects their native
strain/reachability without a second return or an artificial duration.
Source refresh retains every selected Prager and Voce parameter regardless
of builder order. The unselected rate-independent path remains unchanged.

The endpoint heat is `(Y + eta*delta_alpha/dt)*delta_alpha`. Elastic,
energetic isotropic/Prager/Voce storage, coefficient work, mechanical entropy,
and backward-Euler numerical loss keep their original separate accounts.
Existing thermal sections, shells and generic thermal-solid consumers use
the native Hill state directly. This is small-strain associative threshold
flow with fixed anisotropy and material axes; it is not sub-yield creep,
finite-strain texture evolution, damage, or measured material qualification.

`tests/plastic_hill_thermal/rate.rs` adds an independent exact shear
loading/hold/reversal recurrence with varying temperature and duration, a
rotated Voce/Prager thermal first-law check, actual section stress relaxation,
and source-path/time/history refusals. These authored numerical fixtures do
not establish empirical material accuracy.

## Finite combined and thermal creep resistance recovery

`FiniteRecoveringNortonJ2` uses the original simultaneous finite return with
permanent plastic alpha, permanent creep_alpha, and independent primary resistance.
Backward Euler evolves resistance as `r1=(r0+dc)/(1+k*dt)`; its bound uses the
creep path. Recovery preserves the permanent metric and energetic hardening,
creates no heat, and shares preparation with standalone finite recovery. The
native analytic tangent, tangent-free line search, physical mechanism losses and
integrated inverse-rate potential remain authoritative.

`AffineThermalRecoveringCreep` adds independently activated Arrhenius recovery
to the existing affine Kirchhoff/logarithmic thermal family. Complete associated
histories pass through `FiniteThermalInelasticMaterial` and its frozen-point
protocol. Existing J2 families and the default thermal increment state retain
their original APIs. One shared thermal account owns total-F mapping, old-state
inspection at old temperature, coefficient storage, entropy/internal energy,
reversible heat, physical dissipation and the signed endpoint-work defect.
Inspection consumes no duration; recovery adds no stored energy or heat capacity.

`thermal::source::recovery` explicitly binds standalone, simultaneous linear-J2,
and simultaneous Voce-J2 recovering Norton selectors. Law v1/state v3 requires
primary coefficients, an independent recovery reference rate and activation,
and the complete original source context/domain. Original cards and hashes remain
intact; no small-strain or Cauchy-stress calibration is silently reinterpreted.
The caloric consumer still must admit positive total heat capacity. Invalid,
unrepresentable or unsupported histories/rates refuse transactionally.

Focused native and thermal recovery tests retain independent manufactured flow,
potential, entropy, tangent, objectivity, dwell/reload and source checks. These
are bounded authored constitutive models; measured calibration, tertiary damage,
rupture, and service-life qualification remain separate claims.

## Original caloric radial boundary conduction

`IntegratedRadialCaloricBoundary` binds the original caloric material layer to an explicitly cylindrical or spherical reservoir/ambient face. `RadialCaloricBoundaryGeometry` retains physical face/trace radii, complete cylinder length where applicable, and the supplied circumferential or solid-angle fraction. Its exact geometric multiplier is `C = 1/integral(dr/A(r))`; the native distance is the resistance coordinate `A_face/C`, never physical thickness or cell inventory. Both the actual face area and multiplier scale with the declared angular fraction, with no inferred lateral spreading.

The original material surface solve integrates the same k(T), evaluates the actual material-face temperature, and chains the same endpoint derivative through original Cp/Cv storage. The areal film becomes total resistance `R''/A_face` exactly once. The native chain/graph retains the original reference volume, mass, heat, bath entropy and accepted time. Geometry-pinned admission rejects altered native measures, and density-driven planar thickness tangents remain unsupported. The model assumes positive radial traces, source-free quasi-steady legs, fixed geometry/source context and insulated cylinder ends. It supplies no extra thermal inventory, moving-boundary work or phase-front discretization.

Regression definitions cover independent thick-shell logarithmic/inverse-radius fluxes (distinguishable from a planar-thickness approximation), both radial faces, angular fractions, original Cp/Cv derivatives, a supported material surface with an out-of-chart bath, and native implicit heat/entropy evolution using independently computed physical shell mass.

## Current radial caloric conduction

`CurrentRadialCaloricLink` prepares explicit current cylindrical or spherical `LayerGeometry` on the same original integrated caloric/contact owner. Physical radii and explicit cylinder length remain accessible; `LayerGeometry::cylinder_length_m` never infers an extent from resistance coordinates. Preparation checks the original contact and per-side reference coverage, retains the radial family and ordered material sides, and rejects planar/radial or cylinder/sphere replacement. Exact logarithmic or inverse-radius metrics are lowered as `A_contact * integral(dr/A(r))` only for native transport. Every native query requires those prepared metrics plus unchanged reference coverage, and reuses the original k(T) root, original caloric inverse, and full original energy-dependent contact derivatives.

`CurrentRadialCaloricBoundary` similarly retains the original radial boundary's source, material side, areal film and angular fraction while accepting explicit new radii/axial length. It prepares the current exact radial multiplier and `R''/A_current` together, retaining the original material-face support rules and complete energy tangent. Reference mass, cell geometry, covered areas and energy/entropy ownership remain with the native heat owner. These immutable adapters hold the chosen mechanical configuration during thermal trials; they do not infer density-driven distance tangents, physical inventory or mechanical work from an equivalent distance.

Regression definitions independently construct actual shell volumes and volume-preserving motion, compare total-power/contact and film roots against physical cylindrical/spherical integrals, exercise both heat directions and original Cp/Cv contact tangents, reject stale/class/side/coverage changes, and drive native graph and reservoir BE heat/conservation checks. Pinned Rustfmt parsed all changed Rust sources, and independent numerical evaluation checked oracle brackets, supported temperatures and film-area discrimination. Native Rust integration tests were not executed in this constrained environment.

## Original conductivity interval means

`ConductivityCurve::mean_conductivity_w_m_k` evaluates the mean of its original bounded piecewise-affine conductivity over a supported interval, retaining every crossed source knot. Equal endpoints return the point law; reversing an interval traverses the same ordered knots and gives identical arithmetic. Positive convex segment-mean accumulation avoids first forming a dimensional integral that can overflow or underflow even when the mean is representable. No table is copied and the original finite knot bound controls work. This canonical method supplies the directional source adapter in the native tensor conduction owner. Two tests compiled from the actual production curve and original error type passed for crossed knots, support refusal, symmetry and extreme representability.

### Nonlinear anisotropic Hill rate dependence

`plastic::hill::rate::PowerLawOverstressHill48` selects
`alpha_dot=(s0/eta)*max((q(sigma-X)-Y-R(alpha))/s0,0)^n`, with independent
positive eta [Pa s], s0 [Pa], physical dt [s] and finite n>=1. Alpha remains
the native Hill work-conjugate multiplier. Fixed axes, linear/Voce isotropic
hardening and energetic Prager translation retain their existing owners.
The n=1 case delegates to the existing linear law exactly.

The original five-mode return solves one monotone scalar equation, with
analytic kinetic bounds and dimensionless logarithms. Its exact tangent uses
the nonlinear resistance derivative and the positive remainder
`f_rate*(1-1/n)`; neither a secant tangent nor a second material update is
used. Stress-only line searches avoid tangent assembly. Endpoints may retain
static overstress, and subsequent positive-duration holds relax that history.

Physical dissipation is `(Y+f_rate)*delta_alpha`. Kinetic coefficients store
no energy; original hardening storage and material integration loss remain
separate. A bounded scalar solve or unrepresentable state refuses the whole
candidate. Five authored regressions cover prescribed anisotropic multiplier
and slope, variable-duration reversal/relaxation, rotated full tangents,
exact n=1 and isotropic J2 correspondence, and stress-unit covariance.
Independent arithmetic checks are implementation evidence, not measured
material calibration. No sub-yield creep, backstress recovery, evolving
texture, finite-strain forming or damage is supplied by this law.

Both timed Hill laws implement the existing `NativeSectionEnergyLaw` with their
actual physical duration. Isothermal layered plane-stress sections and shell
energy/motion consumers retain the original rate history, yield/viscous heat,
and hardening storage; native mechanical and inertia clocks must agree.
Two independent shear-section load/hold/reversal regressions verify the
linear and nonlinear accounts without a second relaxation during auditing.

### Source-bound nonlinear Hill overstress

`ThermalHillSource::with_power_law_overstress` selects independent original
`hill_overstress_viscosity` [Pa s], `hill_overstress_stress_scale` [Pa], and
`hill_overstress_exponent` [1] for
`alpha_dot=(s0/eta)*max((q(sigma-X)-Y-R(alpha))/s0,0)^n`.
All three claim identities remain pinned. The complete original caloric span,
including every source-line knot/crossing, must keep eta and s0 positive and
finite and n finite and at least one. No hidden frequency/rate axes or
extrapolation are admitted. Selection is exclusive with linear kinetics and
commutes with the existing independent Voce/Prager selections.

Kinetic coefficients may vary nonlinearly with temperature: they have no
stored potential, mechanical entropy or heat-capacity term. Endpoint laws
prepare the actual physical interval once; original-temperature history is
admitted without static projection or an artificial extra interval. The
existing thermal section, shell and Hex8 consumers receive the same Hill
history and exact tangent through the original source-family interface.
Yield plus nonlinear viscous work is physical heat; hardening/coefficient
storage, reversible caloric exchange and integration loss remain distinct.
The n=1 case follows the existing linear arithmetic exactly.

Four focused regressions specify independent scalar load/hold/reversal with
varying kinetics and duration, rotated Voce/Prager thermodynamics and source
pins, n=1 parity, and original-source/duplicate-selection refusal. Coefficients
are authored fixtures, not empirical material validation. Fixed axes and
small strain remain the native model boundaries.

### Source-bound density for exact-Cp kinetic crystallization

`IsobaricCrystallizationPair::try_new_with_free_expansion` binds each original
exact-Cp phase branch to an explicitly selected `PreparedFreeExpansionDensity`.
Existing `IsobaricFreeExpansion` admission requires the same immutable material
card, typed context and temperature endpoints, with explicit density/strain
source support for the original caloric pressure. Phase roles and
matched-specimen applicability remain caller declarations. Reference density
and engineering-length sources are retained directly; no sampled density or
enthalpy chart is introduced.

Common-reference enthalpy and original Cp integrals remain unchanged.
`CrystallizationCaloric` and `ThermalLaw` expose the bound density;
`ThermalLaw` first inverts the rebased phase enthalpy. Branch/pair identities
include expansion binding. The pair's `isobaric_parcel` verifies phase identities
and supplies retained pressure to the existing conserved-mass kinetic parcel.
At crystalline mass fraction X,
`V=m[(1-X)/rho_parent(T)+X/rho_crystal(T)]`. Mechanical reference volume does not
replace current volume. `Delta H=Q` and `Delta U+p Delta V=Q`; latent enthalpy
and pressure work occur once.

Unbound caloric pairs retain their existing thermal behavior and still refuse
isobaric volume construction without phase density. Source-domain failures,
incompatible pressure/history and cancellation publish no replacement parcel.
This is declared free phase geometry at maintained pressure, with no inferred
confinement, pressure-dependent EOS, Gibbs stability, melting support or melt
flow. Regressions cover different original caloric datums, unequal density and
expansion, adiabatic phase shrinkage, subsequent heat, resumed kinetics,
explicit pressure support, phase identity and late cancellation.

## Corning 7980 reference stress birefringence

`silica_stress_optic` retains Priestley (Corning), SPIE 4346 (2001), equation (2) and Table 1, DOI `10.1117/12.435666`: the nominal 193/248/633 nm lines have mean stress-optic coefficients 51/42/35 nm/cm/MPa and between-specimen one-sigma spreads 3/1/1 in those units. They convert to inverse Pa with factor `1e-13`. The source used uniaxial compression, did not specify temperature dependence, and plotted different fit ranges by wavelength; the 1.24 MPa load-cell capacity is not an admitted multi-axial stress envelope. These remain three discrete source lines, without a borrowed continuous dispersion fit.

The evaluator borrows the original tension-positive tensor-shear `Voigt` stress and projects its transverse anisotropy under an explicitly adopted isotropic linear constitutive extension. It returns differential indices and relative phase with the higher-index slow axis; hydrostatic stress cancels. Common index/group shifts, splitting and walk-off are outside this law. Applicability to the caller's stress, temperature and specimen remains the caller's responsibility, distinct from scalar CHARMS provenance.

Validation: the three original source assertion bodies executed successfully in a small dynamic-std executable using the exact production module, exact original Voigt alias and original cached fs-math/fs-blake3. They cover the published SI/one-sigma controls, signed tensor shear and principal stresses, oblique projection/frame covariance, hydrostatic null and refusals. Rustfmt passed. Cargo libtest and Clippy could not write their outputs because the environment was full; no workspace lint or native-renderer test pass is claimed.

### Paired finite equilibrium tanks with integrated enthalpy transfer

`nonideal::flash::parcel::transfer::finite` transfers a prescribed bulk portion
between two existing rigid equilibrium tanks. The original donor S,V withdrawal
retains overall composition and specific entropy while temperature, pressure and
phase amounts evolve. Both phases are entrained in equilibrium proportions.
Its integrated outgoing H drives the original receiving source join and U,V
closure; freezing initial donor h over the finite portion is not permitted.
The receiving inventory is credited with the donor's actual mass debit.

The native integrated-H correction enters the receiving relative-energy target
before original species datums are added. The accepted donor phase state is reused
without another inlet flash. Original NASA/PR identities, known cross-interactions
and common phase policy remain authoritative. Pair mass, each species, reference-
free energy and total entropy production are checked independently; pressure-flow
work and latent energy already belong to the original H/U owners. Both native work
caps remain finite and their combined receipt counters use checked addition.

Both inventories publish together after final cancellation and balance checks.
Source, budget, entropy or pressure refusal leaves both unchanged. Zero transfer
is exact idle. Initial and final pressure-head checks are necessary endpoint
conditions, not a guarantee of the intervening path. This owner supplies no
rate law, empty-tank filling, selective vapor vent, jet energy, heat or shaft work.
Focused regressions specify two-phase discharge and actual endpoint U closure,
dilute analytic blowdown/filling with changing composition, portion subdivision,
large energy/entropy datums and atomic late refusal/cancellation.


## Shared air-relative catalog-glass refraction

`catalog_refraction::AirRelativeSellmeier` owns the original relative Sellmeier coefficients, full fitted band, index-reference Ciddor atmosphere and independently declared wavelength-coordinate atmosphere. It returns absolute vacuum phase/group index and analytic scaled ordinary wavelength derivatives, including both atmospheric product and coordinate-chain curvature terms. Complete source-band admission rejects active poles before intersecting the mapped air support. Boundary inversion occurs once during preparation; bounded, allocation-free queries do not fit coefficients or solve a wavelength root. Additional catalog curvature is evaluated only by the requested jet.

Three original-source tests passed against independent phase/group/curvature controls, distinct atmospheric conventions, source-band and constant-relative limits. The shared owner supports both LUMEN and native Maxwell consumers without reversing layer dependencies. It supplies no absorption measurement, glass temperature evolution, source authentication or causal broadband pulse certificate.


## Finite material walls during molecular vacuum blowdown

`thermochemical_gas::energy::parcel::transfer::species::network::transport::effusion::channel::thermal::vacuum::step_candidate`
couples the original molecular vacuum outlets to finite `FixedVolumeChannelWall`
inventories. All channels sample the original gas and wall states; accommodation
heat is pooled by stable wall identity before one original isochoric storage
inversion. Subsequent calls use the accepted wall temperatures. Shared original
species withdrawal and escaped molecular energy stay with the native vacuum
owner. The final audit replaces the thermostat term with actual wall storage:
`Delta U_gas + E_escaped + Delta U_wall`, in the existing source-relative chart.
Original gas and final wall source domains are admitted before either state is
published. The closed-channel owner shares the same wall inversion and retains
its double-double energy and entropy accumulation.

This is an explicit first-order, fixed-geometry stage for stationary,
collisionless circular channels with diffuse, fully accommodating walls.
Reported wall entropy is a material storage change; vacuum entropy and total
open-system entropy production are not assigned. Optional source-backed wall
conduction is supplied by the extension below. Phase change, adsorption and
in-flight molecular storage remain outside this owner.

Six focused regression groups cover independent monatomic beam/nonlinear-wall
balances, thermal feedback, shared-wall pooling and input ordering, caloric
reference shifts, inactive walls, shared inventory, overshoot and cancellation.
Native Cargo execution did not reach these modules: the dependency build ran
out of workspace space and the executor disconnected. The recovered source and tests match the original implementation after Rustfmt.
Formatting and whitespace checks pass; full native crate execution remains pending.

## Coherent multitone tensor-Prony material energy

`TensorProny::prepare_harmonic_spectrum` retains up to 64 authored positive
angular frequencies, their original immutable tensor law and original material
band. `TensorPronySpectrumWindow` observes the settled real strain field over an
unweighted physical interval using the shared `fs_fft::harmonic_window` pair
moments. The complete Kelvin and Maxwell heat includes every carrier and beat
cross term. Original spring factors determine endpoint storage and directly
integrated storage change; original constitutive tensors independently determine
stress work. Neither a diagonal cycle average nor a work-balance residual is
substituted for physical heat. Binary material, relaxation and strain scales
remain factored until physical results; equal-frequency lines combine with
retained compensation before any product or source-domain check.

The sum of the grouped strain ellipses' peak norms bounds the complete waveform
for source admission, including outside the requested window. It is a
conservative all-time bound and may refuse an admissible waveform; it is not a
certified supremum. Frequency preparation and cached temporal pairs are bounded;
point observations allocate no heap and the controlled evaluator propagates the
caller's cancellation/error channel between bounded material groups. No material
history is mutated. DC/prestrain, startup transients, evolving temperature,
spatial error or empirical source validation are not inferred by this observer.
The physical regression definitions compare independent real Maxwell ODE fields
and quadrature, beat heat, coherent cancellation, one-tone parity, source-clock
translations, original domain refusal and cancellation/reuse.


## Refined diffuse cylindrical-channel transmission

`species::effusion::channel::CylindricalMolecularGeometry::resolve_transmission`
solves the diffuse-wall Clausing equation on nested axial rings. Exact integrated
ring-to-ring and ring-to-end view factors retain reciprocity without angular
quadrature or random sampling. Positive factorizations preserve thin-ring and
distant-ring probabilities. Compensated conjugate gradients solve the
antisymmetric half system with O(N) storage and O(N^2) work per iteration.
Geometry is resolved once; ordinary species and thermal sampling reuses the
retained probabilities and performs no additional geometry solve.

Explicit `ClausingControls` bound ring count (at most 2048), iterations and
absolute/relative probability tolerances. At least three nested grids are used,
except for the exact zero-length aperture. The original operator recomputes its
final residual before acceptance. `ClausingResolution` separates the observed
grid change from both discrete solve estimates. The latter excludes floating-
point roundoff; the combined refinement estimate is not a rigorous continuum
error enclosure. Exhausted budgets and cancellation refuse without publishing
partial geometry. Probability normalization is admitted only within rounding
allowance. No measured performance or empirical wall-model validation is claimed.

The original diameter, length, collision qualification, provenance, direct view
factor and total wall-hit probability remain unchanged. Separate positive
returned and wall-transmitted probabilities preserve wall heat even when total
transmission rounds to one. `DiffuseMolecularChannel` and its existing closed,
vacuum and finite-wall consumers use these same probabilities. The default
constructor retains the documented Santeler estimate; resolution is explicit.
Physical assumptions remain stationary Maxwellian entrance reservoirs,
collisionless flight, diffuse nonadsorbing walls and the consumer's declared
thermal accommodation. No outgoing momentum distribution is inferred.

Equations and diameter-normalized reference values follow Matuda (2015),
DOI 10.3131/jvsj2.58.311, and Gomez-Goni/Lobo, DOI 10.1116/1.1568746.
Twelve native tests of the unmodified production geometry modules, linked to
original `fs-math` with a minimal error/cancellation envelope, passed using the
pinned nightly. Seven new cases cover independent reference transmission
(absolute discrepancy below 1e-5 at L/D=1,2,5,10), view-factor conservation,
refinement, tiny and zero lengths, scaling, determinism, budgets and cancellation;
five existing geometry tests also passed. A source-backed species/accommodation
consumer regression is added but awaits the full material crate build. These
focused tests do not establish full `fs-material` or `fs-couple` integration.


### Adaptive finite-wall molecular vacuum blowdown

`species::network::transport::effusion::channel::thermal::vacuum::trajectory`
advances finite gases and finite material channel walls over one complete
physical window. It shares the existing closed-channel time driver and original
single-step owners: compare a full step against two resampled halves, retain
only accepted halves, and publish both inventories or neither. Error estimates
cover gas temperature/pressure/species, source-relative wall storage, per-outlet
escaped and returned species, relative escaped energy and accommodation heat.
Wall accuracy uses the original caloric integral anchored to its source domain,
so an arbitrary storage datum cannot relax the error tolerance.

The original native vacuum audit consumes borrowed accepted outlet reports;
no source catalogs are copied into an additional flattened history. One fixed
mass/species/energy budget covers the entire original-to-final window. Each
finite wall's actual storage is checked against its accumulated sampled heat,
then replaces the thermostat account. Output row bounds, attempted/accepted
step limits, clock representability and cancellation are explicit. Only typed
step-size/domain refusals retry. Coarse and rejected proposals never enter
conservation accounts. Vacuum entropy is unavailable; the reported wall entropy
is a material storage change only. This remains first-order local temporal
error estimation with fixed geometry/areas, not a global error certificate.

Four focused G1/G3/G4 test groups specify analytic monatomic blowdown, nonlinear
finite-wall feedback versus small native steps, gas/wall energy-datum invariance,
and work/minimum-step/late-cancellation atomicity. An independent monatomic
calculation supports the authored discretization tolerances. Pinned Rustfmt
and whitespace checks pass; native material/thermal integration remains
unexecuted because the shared build filesystem lacks capacity and DSR/RCH are
unavailable. This does not close the broader porous/capillary transport bead.


## Coherent plane-stress harmonic spectra

`TensorProny::prepare_plane_stress_harmonic_spectrum` accepts bounded settled
engineering strain lines in an explicitly ordered material plane. Every distinct
frequency retains the full original complex transverse constraint, including
Kelvin viscosity and all Maxwell branches. No constant Poisson ratio or
independently condensed branch replaces that frequency-dependent recovery.
Exact equal-frequency coordinates combine before recovery and admission; both
compensated words survive the full six-component Mandel map and original energy
contractions. The conservative whole-waveform strain envelope includes eliminated
normal and shear strains and preserves the original material ceiling.

`PlaneStressSpectrumWindow` reuses the original full-tensor finite-window owner
for independent stress work, stored energy change, Kelvin heat and Maxwell heat,
including carrier and beat cross terms. Results are energy densities [J/m³];
plate/shell callers supply their actual thickness quadrature volumes. Up to 64
positive source-admitted authored lines are supported; empty input is zero AC.
Uniform windows, explicit numeric payload accounting, fixed point scratch and
caller cancellation preserve the tensor spectrum contract. DC/prestrain, causal
startup, spatial convergence and temperature feedback are not inferred.

`visco_plane_stress_spectrum` checks independent real-field integration with
frequency-dependent transverse relaxation, all ordered anisotropic planes,
finite-beat heating, equal-frequency cancellation, source-clock invariance,
original full-strain admission, payload bounds and cancellation reuse.


### Partial energy accommodation over repeated diffuse molecular collisions

`DiffuseMolecularChannel::with_thermal_accommodation` binds an explicitly
qualified common per-encounter total-energy coefficient to its original mouth,
source species and circular geometry. Diffuse angles remain unchanged. An N-hit
path thermalizes `1-(1-alpha)^N`: separate transmitted and returned heat weights
come from a discounted Clausing resolvent, not an alpha multiplier on fully
accommodating heat. Geometry and both weights refine on the same rings. One
shared compensated CG implementation solves symmetry-reduced systems with
O(N) storage, bounded work and cancellation. Small alpha uses positive heat;
near-unit alpha uses retained-energy complements. Exact zero/unit and zero-length
limits preserve supplied geometry. Reports expose residuals, refinement changes
and estimates excluding roundoff, not a certified continuum enclosure.

The existing beam-wall helper uses these retained weights with original NASA
relative energies. Held/finite walls, reciprocal gases and vacuum escape share
that implementation. Gross returned molecules still spend the original
withdrawal inventory at alpha=0; an elastic diffuse wall performs no wall caloric
query. The builder binds weights to their matching geometry. The coefficient
requires caller qualification over every catalog species and physical support;
it is not inferred from a bulk material or momentum coefficient. Stationary,
velocity-independent diffuse walls, collisionless flight and negligible
residence storage remain assumptions.

The collision energy definition follows Chiggiato, arXiv:1501.07162, section
2.2.2; repeated-history weights compose it with the Clausing kernel. Six new
native geometry tests and all twelve existing geometry tests pass using actual
production geometry and real fs-math through a minimal error/cancellation
driver. Tests include independent collision-count enumeration, E[N]=2L/D at
small alpha, tiny nonzero heat, exact limits, refinement/scaling, work refusals
and cancellation. Seven additional material/finite-wall consumer test groups
are authored but unexecuted under the shared build-capacity limit. Formatting
and whitespace checks pass; full crate integration and DSR remain outstanding.


## Conservative gas-energy coordinates below a temperature quantum

`GasCaloricBranch::sample_density_internal_energy_coordinate` is an explicit
numeric-stage view for an ODE owner retaining its own conservative energy
coordinate. An otherwise refused nonzero request can sample the unchanged
reference temperature only if the exact original Cv integral to the directional
adjacent float proves it is the closest represented temperature, the neighbor
is inside the admitted branch and the caller's absolute energy tolerance is
met. The unresolved energy is returned in the sample residual and MUST remain
in that separate coordinate. Ordinary converged samples preserve their bits.

This does not relax standalone gas heating or accepted parcel reconstruction:
those strict APIs still refuse a nonzero increment that would disappear into
an unchanged state. Caloric bubble and hydraulic pocket stages retain the
energy coordinate; their original-source endpoint and closed-window energy
checks remain strict. Source failures, bounds and cancellation propagate.

All nine focused original prepared/floating-point-floor regressions passed,
including the retained-coordinate test and the existing strict heat-increment
refusals. The native cavitation source subset passed 82 tests. These are
focused numerical checks, not full workspace/DSR or experimental validation.


## Species-specific molecular channel thermal accommodation

A `DiffuseMolecularChannel` may bind a complete canonical species-to-coefficient/provenance table to its immutable original gas catalog, within the existing 32-species effusion bound. Missing, duplicate, foreign, invalid-probability and unqualified rows are refused. Original caloric, transport and source-pack admission remains authoritative. Common and species-specific coefficient modes are explicit and mutually exclusive; changing the mode replaces the prior qualification.

All distinct per-collision coefficients resolve through one shared Clausing kernel and mass solution on each refinement grid. Duplicate coefficients reuse thermal work. Refinement accepts only when the common mass transmission and both heat probabilities for every distinct interior coefficient meet their estimates. Responses preserve input order; exact zero remains elastic, and coefficient one uses the final shared geometry. Numerical reports distinguish authored batch size, distinct interior thermal solves, shared work, and each response's heat estimates. Storage remains linear in ring count plus bounded coefficient count.

Species-specific weights change only thermal exchange. Transmission and gross reflected species inventories remain geometric. Held-wall, finite-wall, vacuum and trajectory consumers share the existing original-source caloric owner. Absent species and present exactly elastic species do not query inactive wall calorics; thermally active species and accepted endpoints retain their original support checks. No thermodynamic storage, material source or species coefficient is inferred. Coefficients remain authored constant total-energy accommodation laws for stationary diffuse walls within their qualified domains.


## Relaxed DC strain in coherent harmonic observations

The full-tensor and plane-stress spectrum windows expose `evaluate_with_bias`
and `evaluate_with_bias_and_checkpoint`; source `validate_strains_with_bias`
admits the complete prescribed waveform. A finite nonzero DC strain requires
zero frequency in the original material band. Exactly zero bias preserves the
existing AC-only applicability and numerical result. The full six-component
strain envelope includes the constant norm and every coherent AC ellipse.

Only the equilibrium spring carries the explicitly relaxed static strain.
Endpoint storage combines DC and AC before squaring original Gram channels;
its change and independent stress work retain the cross terms through original
zero-frequency FFT moments, without subtracting rounded endpoint energies.
Maxwell and Kelvin heat remains the independently computed AC heat. The
plane-stress bias uses its original zero-frequency transverse constraint,
separate from every frequency-dependent AC constraint; the prepared source
caches that map or its original refusal without adding per-point solves.

This is settled linear constitutive accounting, not a causal startup, retained
initial Maxwell memory, finite-strain response, or geometric prestress model.
The existing line, payload, physical-window and cancellation bounds remain.
The spectrum tests check independent real-field biased energy, unchanged heat,
pure relaxed static storage, a positive 1e-20-second storage/work increment
below endpoint-energy resolution, original DC applicability, duplicate-tone
cancellation, complete strain admission and dispersive plane-stress recovery.


## Exactly evacuated receivers behind finite thermal molecular channels

`species::network::transport::effusion::channel::filling` composes the original evacuated-volume, finite-wall vacuum and reciprocal gas/wall owners. The first stage captures the accepted transmitted species and their wall-adjusted molecular energy inside the authored receiver volume. Its original caloric inversion determines temperature and pressure without seed gas or an assigned vacuum temperature. The source-relative energy coordinate adds only transmitted wall heat to the aperture's RT/(2M) excess; reflected mass and returned heat stay with the original donor. Zero capture preserves the evacuated owner. Subsequent stages use the actual filled receiver's return beam through the ordinary reciprocal finite-wall transaction.

Donor, receiver and material wall are immutable until one candidate passes the original gross-withdrawal, caloric, transport and active species wall-support gates. Common and species-specific accommodation use the same channel helper. Actual source-relative gas storage plus finite material storage must close under the fixed species/mass/energy budgets. Original gas entropy, including newly occupied volume and mixing, plus the material's authored sensible entropy must satisfy the complete closed-system entropy budget. No entropy or thermal state is assigned to the initial vacuum.

This is an explicitly first-order quasistatic step for well-mixed ideal gases in fixed rigid vessels, a stationary diffuse channel and negligible flight time. Fixed vessel supports absorb momentum without boundary work. The conducting-assembly extension below supplies original material heat paths between finite walls; moving/rotating receivers and spatial filling transients still require their own owners. Six focused tests cover wall-controlled first-fill temperature and pressure, reciprocal continuation, closed/elastic limits, species fractionation and accommodation, caloric reference invariance, and atomic late source/inventory/cancellation refusal. Full material-owner test execution remains pending the complete dependency build.

## Original gas candidates in the adaptive network owner

`species::network::advance_candidate_adaptive` accepts original
`SpeciesNetworkCandidate` probes through the existing bounded full-step versus
two-half-step temporal owner. The public controls, report and candidate aliases
retain that owner; no second integrator or caloric solver is introduced.
Callbacks must use the supplied current snapshot, physical duration and original
endpoint/physical controls. Closed-window admission refuses changed parcel
identities, external species paths and external heat. The original source audit
retains fixed volumes, canonical species and caloric/transport authorities.

Only accepted fine candidates enter the compensated mass/partial-enthalpy
ledger. `DiffusionNetworkAdvanceReport::heat` separately integrates canonical
paired heat, including its sign; existing molecular/binary kernels report an
empty heat list. Species and heat rows share the original 4096-row window bound.
Probe energy residuals and actual original-to-final source-relative internal
energy remain independent checks. Each species, total mass, energy and actual
endpoint entropy uses the same fixed physical allowance for the entire window.

Only this callback path adds typed `GasEnergyError::OutsideBracket` recovery,
alongside the existing typed withdrawal, local-error and entropy refinement.
It cannot convert arbitrary source/domain failures into retries. Existing
molecular and binary retry behavior is unchanged. Rejected probes consume their
work budget but cannot contribute heat, enthalpy or accepted state. No
extrapolated source state or second-order/global-error guarantee is supplied.

The native turbulent gas regressions exercise this adapter through source-owned
variable-Cp heat/species mixing, independent time-refinement references,
formation-reference invariance, accepted-only heat accounting, typed
thermal/withdrawal recovery and transactional failure. Scoped pinned formatting
and static review passed; native integration execution is pending the complete
pinned dependency build.


## Saturated-water vapor calorics for nonlinear phase exchange

`cavitation::vapor::SaturatedWaterVapor` retains one original IF97 saturated-liquid/vapor pair and its analytic temperature jet. Pressure follows the original Region 4 derivative; vapor density and internal energy use the original Region 2 Gibbs derivatives along that curve. Raw Region 1/2 caloric references, latent enthalpy and Gibbs mismatch remain visible. No fitted equation of state, ideal-vapor replacement or Clapeyron reconciliation is introduced.

For a uniform saturated vapor volume, `m = rho_v(T) V` and `U = m u_v(T)`. With liquid transfer enthalpy `h_l`, the explicit open-system first law is `dU = dQ + h_l dm - p dV`. The source exposes its temperature coefficient `rho_v du_v/dT + (u_v-h_l) d rho_v/dT` and the original latent enthalpy. This coefficient belongs to an exchanging phase inventory, not the isochoric heat capacity of a closed gas. Fixed-temperature growth consequently requires `rho_v (h_v-h_l) dV` of interfacial heat.

The source query is allocation-free on its successful path, has no iteration, and retains the existing 273.15..=623.15 K Gibbs-pair domain. It supplies no phase geometry, kinetic evaporation coefficient, liquid heat-transfer closure, moving interface, nucleation or critical-region extrapolation. Consumers own those physical assumptions and their admissibility. Five focused tests execute the original source and check independent saturation differences, open-system energy, latent heat, retained endpoint data and domain refusal.

## Shared projected Maxwell spectra

`visco::PronySpectrumWeights` exposes the existing binary-scaled original
Maxwell storage/loss multipliers to coupled frequency-domain consumers.
`new(omega, tau)` retains finite signed-frequency/DC semantics and requires a
positive finite relaxation time; `primal(coefficient)` preserves signed modal
cross terms while refusing nonfinite or overflowing results. Equilibrium
stiffness is neither inserted nor subtracted. Positive loss follows the
material owner's exp(+i omega t) convention, so exp(-i omega t) consumers
negate the imaginary part. Source validity and full-form passivity remain on
the original material/projection owners. Existing internal tangents and
arithmetic are unchanged. All six original spectrum tests passed in a native
standalone module build, including extreme-rate/time and subnormal controls.


## Intrinsic rate-dependent felt in strain and stress units

`fiber::felt_rate::RateFelt` composes the original `WoolFelt` nonlinear,
unilateral compaction envelope with up to eight positive series Kelvin creep
arms. Parameters use Pa and Pa s; state uses engineering strain and energy
density. Geometry, instrument identity, fitted defaults and measured-data
authority remain with the consumer. With no arms, original strain-driven
stress, tangent and hysteresis remain bit-for-bit unchanged.

A prepared duration caches the original exponential branch response and loss
factors. The strain-driven response solves one monotone scalar equation and
returns its consistent tangent, every original creep strain and a separate
trial checkpoint. Stress control follows the original spring inverse and
exact held-stress creep, including free recovery after stress removal. Exact
repeated stress reuses the committed spring state and cannot create spurious
new compaction. Invalid parameters, unresolved physical increments, failed
work accounts and domain refusals leave the caller's state unchanged.

Work integrates the declared quasistatic spring path plus the common
held-stress Kelvin path. Recoverable spring/arm storage, direct positive
Kelvin heat and independent original compaction heat are separate accounts;
heat is never defined by subtracting stored energy from work. Under a varying
strain command the endpoint held-stress discretization is first order and
requires time refinement. Zero-force control chooses the residual spring
strain; it does not determine the width of an open physical gap.

Focused native source tests cover 10,000 histories of unchanged WoolFelt
behavior, analytical creep/recovery and heat, consistent tangents, independent
continuous-ODE refinement, rate-dependent stiffness, stationary stress and
atomic numerical/domain refusals. These synthetic checks do not supply
coupon calibration, G2 experimental evidence or consumer geometry. The shared
owner does not yet replace the example and impact-specific adapters.

## Independently sourced dielectric intensity absorption

`SpectralDielectricAbsorptionStatePoint` projects the existing nine SI vacuum
wavelength and nine typed alpha properties without requiring phase-index claims.
Its resolver selects those 18 canonical roles; adoption retains every additional
original requirement and receipt. Refresh re-evaluates exact selected claims,
including wavelength coordinates, on the same immutable card and coordinate
schema. It does not fit absorption(T), infer scattering or manufacture phase data.

The original tabulated and Sellmeier dielectric sources share the same bounded,
allocation-free piecewise-linear alpha interpolation with this projection. No
resampling, alpha/k round trip, positive-value clipping or transparent-band
extension is introduced. Consumers admit their own interpolation gaps and the
intersection with independently supported refraction. Two focused regressions
specify phase-free hot refresh, exact existing-law parity, original source pins
and typed/nonnegative/ordered ingress. Native execution is not claimed here.


### Original Prony cards for the shared felt owner

`RateFelt::from_prony` algebraically realizes an authored
`GeneralizedMaxwell` compliance as passive series Kelvin arms, then replaces
only its instantaneous linear spring with the supplied WoolFelt envelope.
It introduces no fitted relaxation time, damping ratio or material default.
The eight-active-source-arm cap applies before exactly equal source times
combine. This is not a conversion of a pre-existing finite Maxwell history.

The single-arm path uses analytic positive products. Multi-arm matrix
diagonals use positive complements, avoiding loss of the equilibrium spring
by subtraction. Binary-factored source products retain representable moduli,
viscosities and weak weights across separated scales. In addition to the
native generalized eigenproblem residual and equilibrium compliance, every
original constitutive row must satisfy its unweighted relative equation.
Missing weak eigenvector rotations and unresolved coupling refuse explicitly;
a small global residual cannot authorize an inaccurate weak damping arm.

Four native source test groups cover analytic single-arm limits including a
1e20 modulus ratio, original complex-compliance parity within 1e-12 for normal,
weak and repeated-time spectra, zero/invalid/cap admission, and a weak-arm
counterexample invisible to a total-compliance check. Together with the seven
stepping groups these pass in an isolated actual-source harness using original
WoolFelt, Maxwell, Cholesky, Jacobi and generalized-eigen functions and real
fs-math/fs-evidence dependencies. Full fs-material integration remains unrun;
the broader dependency build was killed while compiling asupersync.


## Explicit Sherwood transport for original spherical evaporation

`evaporation::gas::droplet` exposes phase and reciprocal-heat variants with
separate molecular diffusivity and diameter-based Sherwood inputs. The common
geometry owner converts them to `k_m = Sh D_AB / diameter`; the original phase,
kinetic, inventory and gas/liquid caloric owners still perform the exchange.
Reports retain true molecular D, explicit Sh and the actually used coefficient.
No effective diffusivity, fitted correlation or inferred sensible heat is added.

The existing stagnant APIs call the same owner with Sh=2 and retain the exact
original D/r arithmetic. A finite positive Sh is required; D=0 seals mass
transfer while independently authored sensible heat can still change density
and radius. Area, coefficient and transport inputs remain frozen through one
first-order interval, and final size follows actual mass and original material
density. Transport correlation validity remains with the caller. Focused
regressions specify both phase directions, original-surface equivalence, exact
stagnant parity, independent sealed heat and invalid-input atomicity. The
session could review these sources but could not compile or execute the tests.

## Original Maxwell history under coherent harmonic forcing

`TensorPronyHarmonicSpectrum::prepare_causal_window` and its DD-clock variant
drive an existing `TensorPronyState` over a uniform physical interval. The
candidate preserves the original raw-Mandel Maxwell spring-strain inventory
and cumulative ledgers, returning endpoint strain, stress and memory together
with independently integrated original-tensor stress work, Gram-channel heat
and direct stored-energy change. A nonzero constant strain does not erase or
artificially relax supplied memory.

Start strain must be continuous and cached initial storage must agree with
the original constitutive energy within scale-aware roundoff. The original
AC/DC applicability and complete coherent strain-envelope admission apply.
Law, branch ordering, temperature and reference configuration remain fixed.
The public material state does not establish source identity; its owner must
retain the same law. This is prescribed full-tensor forcing, with no constrained
plane-stress startup or band-limited transient claim.

Shared causal moments retain normalized DD clocks, short-interval ODE
coefficients and positive heat accumulation; unresolved positive energy is
explicitly refused. The existing equilibrium/Kelvin temporal cache retains the
same source clock. Prepared temporal payload and per-point work have explicit
bounds, including fixed point scratch; endpoint history is transactional.
Seven G1/G3/G4 regression tests specify anisotropic retained memory against
independent real ODE integration, native time-domain continuation, tiny startup
heat, free relaxation, settled parity, chained clocks, large constitutive
scales and refusal atomicity. Original branch-state queries also contract
stress work against original tensor columns, avoiding six redundant queries
and an unused stress-square overflow per branch.
The new causal source/tests were not native-compiled or executed because the
execution environment became unavailable; earlier settled tests do not prove
this new path.


## Adaptive evacuated-receiver filling with the original finite channel wall

`species::network::transport::effusion::channel::filling::trajectory::advance_candidate`
advances an exactly evacuated finite receiver, its original finite donor, and
one finite material channel wall over an authored physical window. It returns
an atomic candidate for all three endpoints together with the actual accepted
first-fill, first reciprocal, and continuation reports. The original empty
volume has no invented gas seed, temperature, pressure, or entropy.

Startup compares one full original capture with one half capture followed
immediately by a real reciprocal finite-wall half step. Thus the newly formed
receiver's return beam participates in the startup error estimate; only the
refined pair is accepted. The existing finite-wall
`thermal::trajectory::advance_history` owns the continuation. Seven helper/type
visibility changes expose its existing implementation to this sibling; its
continuing adaptive algorithm is unchanged.

Error estimates include the original gas temperature, pressure and species
inventories; chamber energies relative to the original species' internal
energies at the initial donor temperature; nonlinear wall temperature and
storage; separate directional transmitted and gross reflected species masses;
and positive/negative accommodation heat into each gas and the wall. Zero-mass
species do not query inactive caloric sources. These are local temporal
estimates, not certified global error bounds.

A conservative initial pressure ceiling uses the original transmitted molar
flux, actual receiver volume, receiver solve temperature bracket and admitted
mouth/donor pressure. This protects the original first-fill step without
giving vacuum a thermodynamic state. A ceiling below the authored minimum
step returns a typed initial-pressure-limit refusal. Original source,
composition, caloric and material-domain errors remain hard failures; only
the original native physical/adaptive retry categories may reduce a step.

Startup and continuation share attempt, accepted-interval and retained-report
budgets. Coarse and rejected candidates count as attempted work; their
physical reports are not accepted history. Cancellation remains typed through
the pressure prelude, every physical probe, and final publication. The clock
guard rejects a startup interval swallowed by the represented remaining
duration and requires startup plus continuation to equal the authored window.
For this single-wall entry point, a closed path retains the exactly empty
receiver and unchanged donor/wall for the complete window with one inactive
evaluation. The conducting-assembly extension below instead evolves its
active material heat paths while preserving the evacuated receiver.

The whole accepted window reuses the parent's original gas-plus-wall mass,
source-relative energy, and entropy audit with one fixed physical budget.
Wall mass, volume and source identity must remain unchanged. The original
sensible-entropy material law owns finite wall entropy; the newly formed
receiver contributes its actual gas entropy. All three final states publish
together only after this audit succeeds.

These remain quasistatic well-mixed gases in rigid vessels and a stationary
diffuse channel. The initial non-Maxwellian filling transient and channel
flight-time storage are outside this model.

Validation evidence: the initial production trajectory and audit/helper
changes were parsed/formatted by the pinned Rust formatter before the executor
disconnected. Recovered final edits and four new native test groups were
statically reviewed against immutable original APIs, but were not compiled
or executed. Tests cover independent continuous monatomic mass/internal-energy
balances with nonlinear wall storage, tighter-step convergence, exact accepted
history replay, species accommodation/fractionation, formation-datum invariance,
closed/elastic/pressure-limited startup, work/minimum-time bounds, original
source refusals and early/final cancellation with unchanged original endpoints.
Independent arithmetic probes support the fixtures, but do not constitute
execution of the Rust physical owners.


## Finite diffuse-channel axial momentum and accepted-transfer impulses

The original circular Clausing geometry now resolves incident-normalized axial
cosine moments for direct, wall-transmitted, and returned molecules through
`CylindricalMolecularGeometry::resolve_axial_momentum`. The returned geometry
and all three moments share one final mass grid and escape-probability solve.
The circular cosine-law angular kernels have an analytic axial ring average
and positive eight-point Gauss angular quadrature. Spatial and angular
refinement are separate, cancellable, bounded by the authored Clausing
controls, and use O(number of rings) storage. Reports retain angular work,
spatial changes, linear-solve estimates, and the equilibrium reciprocity
residual. These are numerical estimates, not certified error bounds.
The exact zero-length aperture retains direct moment 2/3 with no wall paths.

`species::effusion::channel::momentum::DiffuseChannelMomentumInterface` binds
those moments to the original channel catalog, mouth qualification, original
species caloric owners, and explicit speed-distribution provenance. Construction
rebinds the existing common or species-specific heat accommodation to that same
final geometry without another geometry solve. It qualifies the original independent per-collision Maxwellian reset or
speed-preserving diffuse event for every species, including intermediate
probabilities. The joint angular/collision-history extension below supplies
matching heat and speed partitions on the same final geometry. An aperture has no active wall and
does not acquire a wall-temperature or accommodation dependency.

The new rate APIs preserve the original full vacuum or reciprocal mass/energy
flux and return per-species incident, transmitted, returned, and stationary-wall
forces in N. Flux-weighted molecular speed is 3/2 times the original thin-mouth
mean normal speed. Each outgoing speed uses its conditional angular moment,
with the authored original gas or wall temperature as appropriate; using the
thin-aperture speed times channel conductance is not the channel momentum law.
`vacuum_outlet`, `reciprocal_step`, and `finite_wall_step` enrich original
accepted native receipts with impulses in N s. They use the actual accepted
transmitted and returned masses, retain original mouth/source identity and
exact geometry/surface qualification, and reject noncanonical, duplicated,
unknown, negative, nonfinite, or omitted active species inventory. These receipt
queries do not resample the caloric owners.

For each directed beam the wall account is incident plus returned minus
transmitted momentum. Exterior forces/impulses use one left-to-right axis;
their sum with the wall contribution closes the molecular momentum account.
A vessel containing its source gas and channel wall recoils against the
escaping momentum. These are stationary-wall local-rest-frame molecular
boundary data; the module adds neither other container tractions nor a
continuum nozzle pressure-thrust term. A moving mechanics consumer must retain
the accompanying changing mass, complete molecular energy, and frame transport.

Validation evidence for this increment: the narrow native geometry executable
using the real fs-math owner passed all 23 existing geometry tests and five
of six new moment tests. The remaining test stopped in its existing mass-only
L/D=10 comparison because its 1024-ring work limit was insufficient. Its limit
was corrected to 2048 without changing tolerances; that native rerun was
interrupted by loss of the executor and is not claimed as a pass. Passed new
checks include an independent 160,000-history diffuse-ray reference, an
independent direct-path integral, aperture/scaling behavior, cancellation and
work bounds. A bounded JavaScript replay confirms the corrected comparison's
numerical work estimate only. Four additional native consumer test groups cover
equilibrium pressure/momentum, beaming and wall-speed scaling, endpoint species
accommodation, and original accepted vacuum/reciprocal receipt identity and
conservation; these were authored and statically reviewed but not compiled or
executed after the workspace became unavailable. Recovered source has no new
native formatting or compilation claim.


### Shared Kelvin creep for simultaneous contact solves

`fiber::felt_rate::{CreepSpectrum, PreparedCreep}` exposes the same intrinsic
series Kelvin owner used by `RateFelt`, without requiring an invented nonlinear
spring merely to convert a Prony card. Cold conversion still retains every
resolved observable source arm and checks the original source equations.
The consumer supplies area/thickness and its own simultaneous mechanical
force solve. No contact geometry, default envelope or material fit enters
the intrinsic owner.

Fixed-clock preparation caches the decay, complementary response and direct
positive loss factor. A weak exposure uses a compensated increment; a strong
exposure uses equilibrium plus the surviving exponential tail. Thus a
representable free recovery at exp(-100) survives even when 1-exp(-100)
rounds to one. Successful warm operations use bounded eight-arm arithmetic,
fixed-size checkpoints, no allocation and no branch exponential.

Held-stress trials report accepted-state work and storage change, independent
Kelvin heat, the actual balance residual and a visible numerical estimate of
state-rounding error. An equilibrated history can remain at the same float
when its physical increment is below checkpoint resolution; that defect is
reported and must fit the estimate. This is not an interval certificate or a
time/model-error allowance. The existing full RateFelt owner retains its
stricter material-state acceptance. Neither API mutates caller history.

Six new focused native source groups pass: analytic creep/recovery and heat,
surviving long-time tails, 25,000 held-stress steps, common-kernel parity with
RateFelt, empty/zero-clock/invalid-state behavior, and direct Prony realization.
All eleven existing felt/Prony source groups also pass. This source harness
uses original numerical functions and real fs-math/fs-evidence dependencies;
full fs-material/runtime integration remains unverified.

## Collective spherical droplets in one finite atmosphere

`evaporation::gas::droplet::chamber` derives each physical droplet's area and
diameter from its positive mass, original enthalpy chart and liquid density.
Explicit molecular D and diameter-based Sh supply `k_m = Sh D / diameter`.
Stable IDs canonically order heterogeneous droplets; duplicate IDs and rosters
above the caller's positive work bound refuse before phase work.

One original chamber owns the shared vapor donor, species update and gas
caloric endpoint. Gross condensation cannot spend another droplet's same-step
evaporation. Optional reciprocal sensible heat uses the original simultaneous
shared-gas heat owner, followed by the same joint phase interval. Original
finite kinetics, caloric reference compatibility, endpoint equilibrium and
inventory/temperature limits remain active. Latent and carried energy are not
added twice. Receipts retain individual phase/heat accounts plus the sole gas
account and noncancelling extensive mass/energy defects.

Geometry and transfer inputs freeze for this first-order interval. Accepted
radii use actual final mass and source density, including heat-only density
changes at D=0. Any source, geometry, conservation or cancellation refusal
returns no partial candidate. Gas volume stays externally prescribed; no
displacement work, statistical parcel weights, droplet shielding, trajectories,
breakup, coalescence, nucleation or complete dry-out is inferred.

Five focused G0/G3/G4 regression groups specify independent spherical rates,
mixed directions/sizes, original scalar/collective equivalence, permutation,
shared donor limits, reciprocal heat and atomic retry. Pinned rustfmt parsing
and formatting checks passed; compilation and runtime execution were unavailable.

## Coherent kinematic envelopes for held thermal expansion

`TensorPronyHarmonicSpectrum::validate_kinematic_envelope_with_bias` admits
both the original material norm and an explicit finite positive absolute
Mandel-component ceiling over the complete waveform. Original equal-frequency
grouping occurs before bounds, so large cancelling authored fields retain their
physical residual. Distinct-frequency envelopes remain conservative all-time
bounds, not sampled peaks or certified suprema.

This is amplitude admission for total or mechanical kinematics. A constant
free-expansion offset does not request a fictitious DC material response.
The actual constitutive evaluation independently admits its mechanical bias
and original frequency band. Both rules are retained when held thermal strains
enter the causal solid history.

The two original-source admission regressions pass natively, including
coherent cancellation, future-phase ceiling refusal, independent norm/component
limits and retained DC refusal for actual mechanical bias. Follow-up native
execution also passed all seven causal material regressions, eleven settled
spectrum regressions and seven one-tone window regressions (27 total). These
used the original material modules and original fs-math/harmonic-window owners;
the latter passed all 90 original window tests including the nine causal tests. Full Cargo
and Cx-backed solid/coupling integration are separate unexecuted gates.


### Prepared original spectral absorption over temperature

`state_point::SpectralDielectricAbsorptionTemperaturePath` prepares the original
nine wavelength and nine alpha claims over a typed absolute-temperature interval,
holding every other source coordinate fixed. Each requirement uses its exact
selected claim through the existing matdb source-line owner; auxiliary roles,
sampled-surface cell crossings, source holes and original interpolation policy
remain effective. The caller bounds the complete retained breakpoint count and
can cancel preparation without publishing a partial path.

Prepared paths share compact original piecewise-affine values. Native
temperature/wavelength queries allocate nothing and perform no material queries.
Wavelength knots may move independently with temperature; wavelength ordering,
positivity, nonnegative absorption and every auxiliary constitutive requirement
are admitted at the complete source breakpoints. Runtime arithmetic still refuses
an unrepresentable or collapsed spectrum. Transparent sources retain their finite
bands. The common wavelength intersection over the whole temperature interval
may be absent even though every local spectrum remains usable.

This adds no absorption-temperature fit, fixed-grid resampling, extrapolation,
phase/absorption relation, thermal evolution or uncertainty enclosure. Four
focused source regressions cover original pins and interlaced knots, fixed-pressure
surface slices and holes, interior crossings and auxiliary-domain exits, bounded
cancellation, and transparent moving support. These regressions are authored and
source-reviewed; native compilation/execution remains unverified in this session.


## Joint angular and collision-history accommodation for channel momentum

`CylindricalMolecularGeometry::resolve_axial_accommodation_batch` resolves
mass transmission, total axial angular moments, native thermal accommodation,
and retained-speed/reset-speed angular partitions on one final Clausing grid.
Coefficients retain their authored order; duplicate coefficients reuse their
solves. The law is the channel's existing explicit independent per-collision
Bernoulli reset: with probability alpha, speed and internal state equilibrate
with the stationary isothermal wall; otherwise both are retained while angles
remain diffuse. A thermal accommodation coefficient alone does not establish
this stronger distribution law; it must be covered by the original qualification.

With q=1-alpha, the discounted visit field solves (I-qW)y=eL. Retained wall
angular moments use q*y dotted with the appropriate end angular kernel; reset
moments use alpha*y dotted with the common angular escape field (I-W)^-1*m.
These retain the correlation between collision count and eventual exit angle.
A scalar heat fraction times a total angular moment is not this joint law.

Each refinement shares one mass solve and one coefficient-independent angular
escape parity pair. Each distinct interior coefficient needs one discounted
parity pair, reused for both heat and angular weights. Only one coefficient's
ring fields are active at a time; scratch storage is O(rings + species).
Tiny reset terms and tiny retained complements are evaluated positively.
Their resolvent partition is checked against numerical residual bounds before
obtaining the dominant component by subtraction. Exact alpha=0, alpha=1 and
zero-length limits remain algebraic. All mass, total angular, heat and split
angular estimates must satisfy the same authored controls before publication.
Reports retain work, residual and refinement estimates; they are not certified
continuum or floating-point enclosures.

`DiffuseChannelMomentumInterface` now admits intermediate common and
species-specific reset probabilities. It binds the returned native heat
responses directly into the original channel, preserving source and surface
qualifications without a second geometry solve. Immutable split responses
are shared between interface clones. Rate and accepted-impulse calculations
combine the original gas and wall Maxwellian mean speeds with their respective
angular partitions. They do not take the square root of an averaged energy
temperature. Source and stationary-wall momentum accounts, original accepted
species masses, exact aperture behavior, and inactive elastic wall semantics
remain unchanged. The finite-wall vessel adapter consumes this same interface.

Focused native geometry validation passes all 35 tests: the original 29 plus
six new groups covering independent collision-generation enumeration and
160,000 three-dimensional diffuse ray histories with counted collisions,
different angular versus thermal reset fractions, exact endpoints/aperture,
shared geometry and duplicate work reuse, alpha=1e-18 with relative-only
control, alpha=1-epsilon, thin passages, scale invariance, exhausted budgets,
and cancellation. The launcher uses actual production geometry modules and a
real fs-math library with minimal enclosing error/poll context; it is not a
full fs-material or workspace build. Pinned formatting passes all nine changed
Rust files. The material integration tests exercise partial-reset wall-speed
scaling, pressure/momentum closure, matching heat weights, inactive species
sources and actual vacuum/reciprocal accepted impulses; those integration tests
were statically reviewed but not compiled or executed in this workspace.

## Causal harmonic plane stress on the original mixed-control history

`PreparedTensorMixedControl::prepare_plane_stress_harmonic_window` drives ordered
engineering strain carriers through the existing source-bound, full-tensor
Maxwell realization. The original coupled memory modes enforce zero transverse
stress throughout the interval. No independent branchwise plane-stress fit or
settled-state initialization replaces the accepted history.

Equal carriers combine before admission. The common exp(-i omega t) clock keeps
its double-double start and separate duration; point evaluation checks strain
continuity, original source identity, the matching plane partition, and the
complete imposed/transient strain envelope. Preparation and evaluation retain
explicit payload/work ceilings and cancellation. Refusal returns no candidate
and never changes the borrowed accepted history.

Physical Maxwell heat integrates positive modal squares. Stress work uses the
original constrained reaction ports independently of direct spring-storage
increments; the original work/storage/heat balance must resolve before a
candidate can be accepted. Endpoint strain, stress, all branch coordinates and
cumulative work/heat remain in the original `TensorPronyHistory`, allowing
continued mixed-control ramps and further harmonic intervals.

This lane retains the mixed controller's fixed-material, fixed-temperature and
zero-Kelvin domain. It does not certify a structural force equation, evolve
temperature, or substitute a sampled maximum for its conservative strain bound.

Focused physical tests: `tests/visco_plane_stress_causal.rs` compares all four
ordered material planes with an independently integrated original-tensor ODE,
then resumes the accepted history through the existing ramp controller. It also
covers retained tiny-interval heat, large-clock subdivision, duplicate carrier
cancellation with memory, transient transverse strain admission, and atomic
source/domain/budget/cancellation refusals.


### Reacting rigid-cell internal-energy storage

`thermochemical_gas::energy::equilibrium::volume::storage::GasIsochoricEquilibriumThermalStorage` connects the original canonical TV/UV equilibrium path to the existing thermal law with `ThermalStorageBasis::InternalEnergy`. Conserved feed mass and rigid volume define density; pressure and molecular composition evolve. Cells may scale mass and volume together, but registration rejects a density mismatch beyond a 64-EPS representation allowance. The existing thermal timeline owns finite heat transfer. Stored U already contains chemical formation energy, and rigid cells add neither a separate reaction-heat source nor pressure-volume boundary work. No open-flow/advection identity is inferred from the cold-relative coordinate.

The source interval is continuous and bounded. Original per-species Cv coefficient bounds, conservatively adjusted for admitted inventory error, bound the reacting inverse slope globally; every queried reacting Cv is checked independently. Complete material endpoints retain original cards, auxiliary axes and transport admissions. The exact advertised upper stored-energy endpoint passes the same relative-U gate as later full queries. Interior caloric support does not promise interior transport support. Canonical scratch restoration makes immutable queries repeatable; refusal/cancellation never publishes a partial material endpoint.

Original material freezing re-closes pressure at exact mass/volume. Its amount-only U correction preserves represented mole differences in double-double arithmetic before multiplication by formation-sensitive species energy. Local entropy keeps reacting composition/mixing and evolving pressure, checks source U against the stored interval and accepted endpoint residuals, and admits its inverse-temperature mean only inside the actual endpoint interval with explicit roundoff allowance. Coincident accepted temperatures use du/T with accuracy limited by the original UV tolerance; unresolved source intervals otherwise refuse without clamping.

`storage::conduction::GasIsochoricEquilibriumConductionLaw` shares this exact caloric law with the original implicit thermal operator. An explicitly authored bulk k(T) curve covers the complete fixed-density branch; dT/du uses reacting Cv, dk/du uses the curve chain rule, and the rigid specific-volume derivative is exactly zero. Nonzero derivative underflow refuses. Regression cases cover reacting Cv/P/U/S, timeline heat/cool conservation, dissociating conduction derivatives, rigid-density and transport-domain refusal, failed/late-cancelled query repeatability, and formation-sensitive freeze corrections.


## Matched liquid calorics for resolved vapor-cavitation heat fields

`cavitation::vapor::liquid::SaturatedWaterLiquidThermal` retains the original IF97 liquid enthalpy used by the vapor mass-transfer port, its analytic saturation-path derivative, and original industrial water conductivity and viscosity. Its storage derivative is `dh_l,sat/dT = cp + (v - T dv/dT|p) dp_sat/dT`. The fixed-pressure Gibbs heat capacity and the relative pressure-path correction remain separately visible; a consumer cannot silently label the saturation-path derivative as cp.

This is an immutable near-saturation incompressible caloric reference for matching liquid storage and advective enthalpy to the existing vapor interface. It does not assign the actual liquid pressure to saturation, invent a metastable equation of state, or include pressure-motion heating. Density, thermal expansion, isothermal compressibility, sound speed and the unmodified source pair are exposed for consumer applicability checks. The source retains the intersection of its thermodynamic and transport domains, 273.16..=623.15 K; the 273.15 K IF97 caloric extension is refused by the original transport source.

Successful queries allocate no memory and require no root solve. Four new tests compare the analytic capacity with independent original-source differences, establish exact liquid/vapor transfer-enthalpy identity, retain actual source transport and the pressure-path defect, and exercise source-domain refusal. All 29 selected material-source and prior caloric/phase regressions passed.

## Adaptive collective spherical evaporation

`evaporation::gas::droplet::chamber::adaptive` evolves a fixed roster of positive
spherical liquid inventories in one original fixed-volume gas parcel. Each
comparison runs one original coarse interval and two original fine intervals;
only the conservative fine pair is accepted. Actual interval-start mass sets
radius, area and transfer length again, with shared vapor depletion and original
caloric closure retained. The common callback also permits original source-backed
intervals without inventing film coefficients in the material controller.

Local first-order error includes each droplet's mass, enthalpy, temperature and
radius, the shared gas temperature, and every original gas species inventory.
Relative scales use coordinate changes from the comparison start, keeping
enthalpy datum shifts out of accuracy decisions. Original method-order arithmetic
handles unequal representable half spans. The accepted endpoint is not
Richardson-extrapolated; this estimate is not a global trajectory certificate.

Only original positive step-size limits permit physical retries. Source,
admissibility, conservation and cancellation refusals retain their identity.
Zero duration still performs original physical admission. Explicit comparison,
reduction and span bounds constrain work, with atomic refusal and reusable row
buffers; no trajectory history is retained.

Accepted fine intervals alone enter the compensated physical account and final
payload. Noncancelling mass and energy defects share one complete-window
allowance. Final per-droplet, per-species and total gas inventory changes are
checked independently, including carrier species. Latent energy is a subaccount
of the original phase exchange, not an extra gas energy source.

Four regression groups cover refinement against 512 original intervals,
reciprocal sensible heat and enthalpy-datum invariance, canonical permutation and
resume, and physical step reductions/budgets/cancellation/zero-time admission.
Formatting and bounded static review passed; these tests were not compiled or
executed because this workspace has no runnable Cargo toolchain. Complete
dry-out, nucleation, drag, shielding and moving-volume work remain outside this
positive-seed, well-mixed chamber path.

## Viscoelastic coatings on a mechanical backing

`elastic::acoustic::viscoelastic::layers::SupportedViscoelasticLayers` prepares
normal longitudinal reflection from one original fluid into bonded
viscoelastic layers terminated by
`ViscoelasticBacking::{TractionFree, Clamped, Impedance}`. Each coating retains
its original bulk/shear Maxwell sources, identities, density, frequency band,
branch budgets and strain admission. Empty coatings retain the direct
fluid/backing boundary; zero-thickness layers still enforce material domains.

The native convention is `exp(+i omega t)`, compression `p = -sigma_zz`,
velocity directed into the backing, and `p = Z v`. A finite backing requires
finite `Z` with nonnegative real part. An adapter whose original mechanical
law uses `exp(-i omega t)` conjugates that law's impedance explicitly. This
normal scalar model supplies neither oblique mode conversion nor a causal
realization of an arbitrary impedance snapshot.

Driving reuses the stable scattering stack and original positive bulk/shear
strain-work kernels. `SupportedViscoelasticPower` reports total incident,
reflected, intrinsic bulk, intrinsic shear and backing powers in watts over
the requested area. Backing work is independently
`area * Re(Z) * |v_rear|^2 / 2`; no reflected-power deficit supplies heat.

Actual rear compression and velocity, every layer's interior field and
per-layer/subinterval intrinsic absorption remain available. Rear motion is
retained directly through the scaled termination, including a finite stiff
backing whose reflection rounds to the clamped limit; its exact rear field
uses that trace. Free and clamped boundaries retain exactly zero traction or
velocity. Cancellation, unrepresentable solves and original strain limits
fail before response publication. Warm drives and observations allocate
nothing and do not query the constitutive sources again.

Five focused regression definitions cover an independent dissipative
two-layer reference and positive heat partition, empty/zero/lossless phases,
free/clamped/stiff supports, sub-ulp intrinsic heat, opaque shielding, scaling
and original source-domain/cancellation admission.

## Causal Kelvin–Maxwell harmonic plane stress

`kelvin_control::PreparedTensorKelvinPlaneStress` binds the original immutable
tensor law and source-bound history. Engineering in-plane carriers and fixed
bias drive all original Maxwell memories plus the transverse elastic state;
zero transverse traction determines the actual transverse strain rate.
Preparation uses the native positive QR/SVD factor and each harmonic interval
reuses the original causal decay moments. No material branch is fitted,
discarded or replaced by a settled-state initialization.

The full original Kelvin tensor may be positive semidefinite. Its transverse
block may be resolved positive definite, or have rank zero, one or two when the
nonzero transverse projections of the original viscosity Gram columns are
independent and numerically resolved. Redundant or unresolved projected factors
retain an explicit refusal. Material coefficients and temperature are fixed.
The existing no-Kelvin mixed controller and full-stress Kelvin/DAE owner retain
their domains and APIs.

Point evaluation returns an unpublished original history, endpoint strain rate,
direct strain increment and conservative whole-interval strain bounds. Original
source identity, AC/DC coverage, coherent carrier grouping, strain continuity,
material limits, resource ceilings and cancellation remain admission gates.
Original reaction work, positive combined Kelvin/Maxwell heat and direct spring
storage increments are calculated independently. The account does not infer a
separate Kelvin-versus-Maxwell heat partition.

`tests/visco_kelvin_plane_stress_harmonic.rs` contains independent original-tensor
ODE references for all four ordered planes, retained-memory and tiny-interval
heat, pure-Kelvin cubic startup, DD subdivision, duplicate cancellation and
transactional/domain refusals. These tests require native execution; source
review or formatting alone is not numerical validation.

## Original tabulated thermal-refraction cells

`SpectralDielectricTemperaturePath` prepares the existing nine original
wavelength/index/absorption claim roles over one typed absolute-temperature
interval. It reuses the pinned absorption-line owner on the complete original
resolved bundle, retaining every auxiliary requirement, fixed physical
coordinate, source receipt and actual interpolation breakpoint. Source-knot
budget and cancellation remain atomic. All phase values must satisfy the native
`n >= 1`, `n² <= f64::MAX/64` domain throughout their original lines;
positive values alone do not establish a usable geometric-optics index.
The cached maximum original wavelength gap includes every source-breakpoint
extremum, so consumers can retain their interpolation-gap ceiling throughout
thermal motion of the measured grid.

A monochromatic `PreparedSpectralDielectricTemperatureCell` admits the entire
declared temperature interval only when its active wavelength, index and loss
lines are exactly affine and the original moving wavelength bracket stays
fixed. Exact equal-slope continuation is admitted; genuine temperature cusps
and interior moving-knot crossings refuse derivative preparation. Phase-only
values remain available at source cusps. Thermal endpoints use the derivative
from inside the admitted cell. Stationary interior spectral knots retain
thermal jets, with group index absent when original spectral slopes disagree.

Hot samples evaluate original phase, independently sourced intensity absorption,
`n_g = n - lambda dn/dlambda`, and analytic first/second temperature
derivatives of the rational moving-grid interpolant without material queries,
allocation, finite differences or refitting. Finite nonpositive group indices
remain source derivatives for consumers to admit or refuse for their clocks.
Transparent spectra retain finite source support; a disjoint global wavelength
intersection does not erase valid individual-temperature value queries. This
owner supplies source optics, not a temperature solver, birefringence model,
complex-interface law, uncertainty enclosure or invented smooth continuation.

## Normal reflection from original anisotropic elasticity

`elastic::acoustic::interface::normal::NormalElasticCharacteristic` retains the
original full elastic state and a checked proper surface-to-material frame.
The native Christoffel owner supplies ascending-speed quasi-modes, their
polarizations, impedances and energy-velocity walkoff. Original tensor strain
and traction, Cartesian characteristic/admittance tensors and mode mappings
remain available without an isotropic projection or another eigensolver.
Velocity uses the same native polarization in both directions: free reflection
is +1 and clamped reflection -1 for every mode. Degenerate mode labels identify
only a subspace. Velocity and traction use surface axes; Mandel strains retain
the original material axes.

One native 3x3 factorization prepares full passive tensor reflection and a
separate wall-motion map, retaining tiny reflection near matching and weak
motion at stiff supports. Original Gram-channel work stays positive and
separate from independently reconstructed wave flux. The welded two-material
owner uses a native 6x6 solve for all velocity/traction continuities, coherent
incidence from both sides, and each original source's whole-halfspace strain
bound. Prepared drives allocate nothing; preparation is cancellable.

Six authored regression groups cover known full-tensor eigenpairs, an
independent noncommuting-support oracle, frame covariance and degenerate
isotropic parity, stiff/matched limits, source bounds, and an independent
two-sided welded-interface oracle. This is lossless normal-wave bulk physics;
it adds no oblique anisotropic branches, attenuation, source calibration or
experimental qualification.

### Resolved normal anisotropic coatings

`elastic::acoustic::interface::normal::coating` binds original incident and
coating characteristics to a nonnegative thickness and the existing passive
backing. Each checked frame maps the same physical surface axes into its own
original material axes. Native Christoffel modes propagate with their original
speeds and full tensor strain/traction; the existing double-double phase owner
retains phase-product and rear-distance low bits. One row-scaled native 12x12
solve caches reflection, all six interior amplitudes and independent rear
velocity. Weak stiff-support motion and original positive Gram-channel work
remain observable separately from cancelling wave traces.

Prepared coherent drives and depth queries allocate nothing. Whole-halfspace
and whole-layer Mandel norm bounds enforce both original source domains,
including internal standing-wave strain peaks between safe face samples.
Front continuity, rear support fields and independent normal fluxes must close.
Exactly zero thickness delegates to the original supported incident half-space
and has no invented coating field. Preparation remains bounded and cancellable.

Six authored test groups cover an independent full-tensor Cartesian continuum
reference with interior strain, the quarter-wave transformer and original
isotropic owner, independent material-frame rotations, exact absent-layer
behavior, stiff-support work, and internal strain/domain/cancellation limits.
The coating is lossless; only the supplied backing absorbs work. This adds no
oblique anisotropic propagation, multilayer recursion, bulk damping, finite
lateral geometry or experimental material qualification.

## Source-backed conduction through finite molecular-channel walls

The closed and vacuum finite-wall owners expose
`step_candidate_with_conduction`; their existing `step_candidate` functions
delegate with no conductive links. `thermal::conduction::ChannelWallConduction`
retains an original `FixedVolumeConductionLink` or an original prepared
`FixedVolumeLayerConduction` with its explicit nonlinear work controls.
Endpoint IDs, caloric/source identity, reference mass basis, fixed distances,
area, aligned scalar conductivity and area-specific contact resistance retain
their original meaning. No deposition geometry, thermal capacity, source
curve or thermostat is invented to join the state owners.

The volumetric source owner supplies read-only binding access and stateless
temperature entry points to its existing conductivity kernels. Its ordinary
timeline APIs keep the original cell/source/geometry admission before using
those same kernels. Cellwise conduction evaluates each original endpoint k(T)
and the native series/contact resistance. Integrated conduction uses the
already prepared original curves and native source-free one-dimensional layer
solve, retaining interface temperatures, endpoint derivatives and residuals.
No material chart is copied or prepared inside a time step. Every conductive
link must match the actual finite wall laws; a changed mass basis, datum,
entropy interpretation or source is refused. An integrated path admits both
original source spans at the initial and final wall temperatures, without a
second nonlinear solve. Missing material support remains a hard refusal.

Molecular and conductive heat sample one original gas/wall snapshot. Each
conductive transfer enters the two original wall pools with opposite signs,
then all molecular and conductive contributions pool before one original
material inversion per affected wall. Walls without a channel can conduct
heat through the same finite inventory. Per-wall receipts separate molecular,
conductive and net heat; per-link receipts retain direction and the native
source response. Closed paths query no inactive conductivity. Canonical link
ordering, duplicate/unknown endpoint refusal, existing work bounds and
cancellation apply before any endpoint is published.

Actual nonlinear material storage and entropy remain authoritative. The closed
system checks gas plus all finite walls against its original energy and
entropy budgets. The vacuum audit replaces only the native molecular bath
term with actual wall storage; internal conduction cancels across finite
walls. Vacuum retains no invented temperature or total entropy production.
Refusal leaves all original gas and wall states unchanged.

Both native adaptive owners expose `advance_candidate_with_conduction` through
the existing, unchanged `advance_history` time controller. Each trial resamples
the original conductors. Error estimates control molecular and conductive
wall heat separately and each conductive link's two heat directions, so
opposing mechanisms or a reversal cannot hide transfer error. Conductive
receipts count toward the existing retained-row bound. Accepted histories
alone enter one fixed original-to-final budget; every wall's actual storage
must agree with its accumulated net heat. These are explicit first-order local
temporal estimates, not global error certificates.

Six focused G1/G3/G4 test groups cover independent nonlinear volumetric energy
and entropy integrals, reference-density/datum invariance, the analytic
Kirchhoff flux for linear k(T), reversed links, simultaneous molecular and
conductive pooling in closed and vacuum paths, original-law and final cross-leg
support refusal, analytic two-wall relaxation with tighter-step convergence,
accepted native replay, conductive-heat-only adaptivity, closed paths,
reordering, work bounds and late cancellation. The working implementation and
tests were parsed/formatted by the pinned Rust formatter before the executor
disconnected. The final GitHub blobs recover those focused source edits from
their original published bases and pass whitespace checks, but recovered
formatting has not been rerun natively. These integration tests were authored
and statically reviewed but not compiled or executed; no full material,
Cargo or DSR validation is claimed. This does not close the broader
porous/capillary transport bead.

Each material wall is still a uniform-temperature finite inventory. Links
describe declared material heat paths between inventories, not spatially
resolved temperature variation along one molecular bore. Geometry and wall
mass remain fixed; expansion work, phase change, adsorption and finite flight
time remain outside this channel model.

### Shared reacting reservoir outflow quadrature and caloric reference

The existing reacting `EquilibriumGasReservoir` retains a native opaque caloric
anchor at each accepted physical state. Its finite withdrawal energy uses the
original direct current-minus-anchor source integral; the immutable birth
isentrope and elemental ratios remain unchanged. The finite exported stagnation
enthalpy still accounts for flow work and nozzle kinetic energy exactly once.

`blowdown::midpoint` is the common bounded implicit mass equation
`dm = dt*mdot(m_old-dm/2)`. The original rocket caller delegates to it and keeps
its original initial/midpoint/final flow and reservoir admission. Safeguarded
Newton proposals use an optional actual rate derivative; the measured mass
residual accepts the result. The outgoing fraction is at most one half and
excess duration refuses with a midpoint-rule retry bound, without clipping.
Exact zero rate preserves mass without a trial; final cancellation refuses.
`blowdown::prepared` shares the original birth-SP workspace, density inversion
and finite-energy candidate with other flow consumers. An intervening throat
query cannot be published as the remaining tank: final inventory is resampled.


### Pure-fluid PR isobaric latent thermal storage

`thermochemical_gas::nonideal::equilibrium::isobaric` binds the original pure
PR/NASA source and accepted saturation to the existing ThermalLaw enthalpy
coordinate. Saturated-liquid h is zero; sensible branches use original stable
PR responses. The inclusive latent interval selects phase amount at constant T,
with zero dT/dh and nonzero specific-volume derivative. Temperature alone at
coexistence cannot select the amount. The global inverse-temperature bound
uses the original analytic NASA Cv floor and native nonnegative alpha-curvature
proof; critical, multicomponent and unproved domains refuse.

Complete phase states publish through the original material owner with source
receipts and its distinct isochoric response. Thermal entropy integrates dh/T;
a separate original-source entropy query and visible H-T Delta S residual retain
accepted coexistence error. Direct local H and reciprocal-temperature mean
gates refuse unresolved intervals without clamping. The implicit adapter
borrows this exact caloric law, exposes phase fraction for time control and
requires authored bulk k(T), explicitly independent of quality on the plateau.
It does not infer morphology or average phase transport. Fixed-pressure
expansion belongs to enthalpy; reference geometry is not a rigid-volume lock.

One-sided latent endpoints remain usable for finite stepping. The explicit
native harmonic admission refuses both exact derivative corners, while
interior latent storage keeps zero temperature response and finite energy
response. No second heat integrator, phase kinetics or open-flow source datum
is added. Five original-material regressions cover phase publication, latent
heating/cooling, derivatives, datum invariance and refusals. They and the added
native harmonic-corner consumer regression are authored but unexecuted behind
the runtime build limit; final source files parse and received independent
source review. The underlying five PR caloric tests executed successfully.

## Fully relaxed source-bound material initial conditions

`TensorPronyHistory::initial_relaxed_strain` and
`initial_relaxed_plane_stress` create explicit initial conditions from an
original virgin history template. They retain the identical immutable model
binding and every original Maxwell row, with zero branch coordinates and zero
cumulative work/heat. Nonvirgin templates refuse; initialization never resets
an accepted history. Initial stored energy is supplied by the caller as part
of the initial condition, not manufactured as work or thermal deposition.

Full-strain input uses tensor-shear Voigt coordinates. The plane method accepts
ordered engineering in-plane strain and reuses the original DC harmonic
plane-stress elimination to recover the complete equilibrium strain. Initial
stress and storage come from the original equilibrium Gram, with transverse
traction checked before publishing a plane state. Kelvin viscosity has zero
initial rate and supplies no initial stress or loss. Nonzero static strain
requires original DC coverage; the zero relaxed state remains admissible for
an AC-only source. Full-strain validity, representability and cancellation
remain mandatory.

These initial conditions support independently admitted held thermal free
strains and restrained thermal preloads. They do not solve assembled force
equilibrium, infer a preceding loading path or relaxation duration, evolve
temperature, or authorize a changing material source. Four focused regressions
cover analytic plane elasticity in all ordered planes, original source binding,
tensor-shear energy, and atomic domain/reset/cancellation refusals.

### Original thermochemical conductivity intervals

`thermochemical_gas::conductivity::PreparedThermochemicalConductivity`
prepares selected original pure-species viscosity and conductivity paths
through `ClaimSet::query_line_typed` on a complete continuous temperature
bracket. It owns their receipts, canonical composition, material and pack
identities, fixed pressure/auxiliary context and union of all source knots.
Static source support, unchanged selection and positive values must hold
through every interval; supported endpoints alone do not establish admission.

Point queries reuse the original Wilke interaction and the same WMS sum as
`ThermochemicalGasModel::evaluate`. No NASA inverse, EOS call, material-card
selection or mixture fit is introduced per point. The normal original
O(N^2) interaction and small evaluation receipts remain. Interval means
integrate that actual mixture on each crossed source segment with bounded
Gauss-7/Kronrod-15 quadrature, normalized to avoid overflowing an unnecessary
dimensional integral. Reversed bounds use identical arithmetic; equal bounds
query the original point.

Absolute/relative mean tolerances, a complete-request subdivision allowance
and fixed depth/stack cap admit the numerical error ESTIMATE. This is not a
rigorous enclosure, empirical source uncertainty or physical heat deposit.
The point checkpoint runs once per original mixture sample, allowing a caller
to share cancellation/work authority across cells and solver stages.
Preparation checkpoints surround the existing bounded synchronous source-line
operations separately. Refusal returns no partial integral or updated state.

Five tests compiled from the exact production quadrature leaf passed for
crossed knots/reversal/equality, a steep endpoint layer, representable normalized
means, controls/work/depth limits and cancellation. Four original-material
regressions cover point parity, an independent nonlinear two-species integral,
source holes/nonpositive paths/context, and query refusal. The latter tests
were parsed and reviewed but remain unexecuted under the current dependency
and storage constraints.

## Conducting finite-wall assemblies during molecular receiver filling

`channel::filling::step_candidate_with_conduction` extends the original
exact-vacuum receiver owner to a map of finite material walls and the original
`ChannelWallConduction` links. The existing single-wall API delegates with an
empty conduction graph. First capture uses the native finite-wall vacuum
transaction; a filled receiver uses the native reciprocal transaction. Both
sample molecular and conductive heat together before one inversion per wall.

Captured energy contains only the original transmitted molecules and their
accommodation heat. Conductive heat into an auxiliary wall is not added to the
receiver. Every original material inventory, including walls without a molecular
passage, contributes its actual storage and sensible entropy to the closed
donor/receiver/wall audit. Source identity, mass basis and final conductivity
support remain those of the original material owners.

`filling::trajectory::advance_candidate_with_conduction` carries the complete
wall map through the existing positive-capture startup comparison and native
adaptive history. It controls gas state, source-relative chamber energy,
individual wall storage, gross molecular transfers and separate signed
molecular/conductive heat. Accepted histories alone enter fixed whole-window
mass, energy and entropy checks; each wall's stored change is also compared
with its accepted net heat under one original-to-final allowance.

A closed molecular passage may still have active wall conduction. That path
uses the existing vacuum thermal history with zero molecular escape and leaves
the receiver exactly evacuated. Its distinct `ClosedPath` report contains no
invented first-fill event, vacuum temperature or seed inventory. All gas and
wall endpoints publish atomically after the final audit and cancellation check.

The existing well-mixed, fixed-volume, quasistatic and uniform-wall assumptions
remain binding. Conduction connects separate inventories; it does not resolve
temperature variation along a bore or the spatial collisionless filling
transient. Local step-doubling estimates do not certify global model accuracy.

Five focused G1/G3/G4 test groups specify independent nonlinear storage and
entropy, continuous gas/two-wall refinement, accepted-history replay, analytic
closed-port conductive relaxation, original-source/work/cancellation refusal,
and unchanged single-wall delegation. These are numerical implementation
fixtures, not empirical validation of the declared material or flow regime.

### Supported viscoelastic coating heat on the original thermal timeline

`viscoelastic::thermal::SupportedViscoelasticHeating` binds the original
`SupportedViscoelasticLayers` to complete physical depth bins on the existing
enthalpy simulation. Exact native interval work supplies separate bulk and shear
heating. It shares the existing two-fluid heating owner's partition checks,
source identities, caloric/mass/geometry guards, power envelope, cancellation
and atomic accepted-time observer. The original two-fluid constructor,
source/input/power getter types and behavior remain available.

Backing resistance work follows the same accepted exposure independently.
Without a receiver it is reported as `exported_backing_j`.
`with_backing_receiver` can instead bind one explicitly authored active thermal
cell distinct from every coating bin before the first advance. Its original
mass, geometry, caloric identity and storage basis remain authoritative; no
backing thermal mass or constitutive law is inferred from impedance.
`backing_thermal_input_j` is a subset of total native `thermal_input_j`.
Backing elastic and kinetic storage never becomes heat.

The physical acoustic balance subtracts all backing dissipation; the
thermal-credit balance includes only explicitly received backing work.
Empty or zero-thickness supported coatings retain their direct boundary and
can heat only the explicit backing receiver. Coefficients remain frozen under
the caller's declared coating and support temperature bands: no temperature
refresh, carrier transient, resonator ring-up or mechanical-energy release is
inferred. The native thermal owner continues to control conduction, boundary
events, material storage and accepted durations.

The supported-heating regression fixtures cover tiny backing work and intrinsic
absorption below reflected-power roundoff, native enthalpy retention when the
reported temperature rounds unchanged, receiver and source-identity refusal,
zero-thickness support, and compatibility of the original two-fluid API.


### Original singular-transverse Kelvin plane-stress constraints

The original `PreparedTensorKelvinPlaneStress` now admits directional Kelvin
viscosity without adding a transverse dashpot floor. Positive QR resolves the
original null-traction constraints; a restricted positive factor evolves the
remaining elastic coordinates and every original Maxwell memory. Exactly zero
transverse Gram projections remain direct in-plane viscosity. Nonzero projected
columns must be independent and resolved; arbitrary redundant PSD factorizations
are not admitted by this extension. The full-rank path retains its original
preparation and evolution arithmetic.

`transverse_viscosity_rank()` reports the admitted rank. Public
`state_coordinates()` and the existing resource ceiling retain the conservative
original coordinate count, while dynamic modal loops use only the constrained
coordinates. Pure in-plane Kelvin loss with no Maxwell branches has no dynamic
modes: original elastic traction closure and direct positive heat still apply.
No original Maxwell coordinate is removed or refitted.

An accepted original history must already satisfy the null-traction constraint;
preparation does not project or jump its initial strain or memory. The endpoint
null strain rate follows the differentiated original spring constraint. Combined
Kelvin/Maxwell heat, reaction work and spring storage remain independently
accounted, with the original source, strain, frequency, cancellation and atomic
candidate gates. A separate heat-channel partition is not inferred.

Four focused material tests specify independent original DAE trajectories in all
four ordered planes, ranks zero through two, retained memory and changed/held
banks, an analytic zero-mode heat case, incompatible histories and unresolved
rank refusals. The downstream original laminate regression now exercises the
admitted directional law. These sources passed pinned Rust parsing; the new
numerical fixtures have not yet been executed in a native build.

## Conserved original gas inventory for mixed cavitation

`cavitation::vapor::gas::ConservedBubbleGas` retains the exact mass,
composition, original NASA-9 calorics, transport sources and typed auxiliary
conditions of a `ClosedGasParcel`. Its enclosing bubble supplies the current
temperature and volume; this material view has no independently accepted gas
clock or mutable second inventory. Preparation admits the existing continuous
caloric branch. Successful trial samples evaluate original Cp/Cv and the
fixed-mass ideal EOS without allocation or caloric inversion.

The original `m_g integral Cv dT` supplies finite energy changes without
subtracting formation energies. Full endpoint qualification resolves the
original material at the dry component's actual partial pressure, retaining
the original auxiliary query conditions and checking fixed-density coherence.
It never reconstructs the conserved mass from rounded endpoint pressure.
Explicitly declared water-bearing gas is refused; gas-phase stability and
negligible dissolution of the retained species remain consumer conditions.

The returned conductivity/acoustics of a qualified dry component are not
properties of a water-vapor mixture. A mixed cavitation consumer must
separately admit total inertia, partial versus total pressure, interfacial
species momentum, thermal equilibration and species diffusion using physical
wet-mixture data or explicitly authored bounds over its declared domain.

Four new test groups exercise nonlinear original capacities and EOS,
analytical caloric increments and formation-reference invariance, complete
source/auxiliary-coordinate retention, and source-domain/water/transport
refusals. Source and API review was completed against the current original
owners. The executor became unavailable before this new component could be
compiled, formatted or tested; no new runtime result or full-workspace pass
is asserted for this checkpoint.

## Complete-network droplet time accuracy and physical windows

The material droplet adaptive network module supplies shared accuracy and
accepted-transfer accounting for one droplet cloud coupled to an entire original
finite gas network. Every parcel participates, including idle and droplet-free
volumes. Numerical controls and explicit parcel/species work ceilings are admitted
without querying liquid or transport sources.

The local coarse/fine norm retains the original first-order component scaling.
It covers droplet enthalpy, temperature, mass and radius, and each gas parcel's
temperature, total mass and union of species inventories. Relative scales use
changes from the original comparison state. Initial liquid temperature and
geometry come from the original coarse interval's receipts; an external source
is not sampled before the enclosing transport-then-phase owner admits it.

GasDropletNetworkWindow records only consecutive accepted original intervals.
It retains sparse per-parcel species accounts and the final liquid states,
with no trajectory history. Each record consumes the previous account, so a
refused partially accumulated account cannot be reused. Original transport,
joint sensible heat and phase caloric defects share one noncancelling mass and
energy allowance for the complete window. Formation-reference totals are not
subtracted to manufacture a new energy residual, and latent energy remains a
subaccount of the original phase exchange.

Final actual inventories are checked independently for every droplet, gas
parcel and species, plus global gas/liquid mass and species closure. New arrivals
consume the same explicit retained-species work allowance. Closed internal gas
transport and its original molecular energy corrections are supported; external
gas streams or external heat need their opposite physical owner and are refused.
A gas-only entropy receipt is not a coupled liquid/gas passivity certificate.

Focused regressions exercise remote-parcel and newly arriving species error,
manual original-interval composition, and independent whole-window balances.
Compilation, runtime tests and rustfmt were unavailable after the execution
workspace disconnected; no executable validation pass is claimed.

## Positive finite-interval tensor-Prony interference

The shared tensor-Prony/plane-stress spectrum window now retains bounded real
temporal field modes for eligible Uniform intervals. It combines original
binary-scaled constitutive channels before Kelvin/Maxwell heat squares and
independent stress:strain-rate work products. Centered coefficients are exposed
read-only by both shared owned/borrowed temporal factor APIs so constitutive
scales and original relaxation times need not first fit ordinary DD amplitudes.
Bounds and profile/parity requirements remain explicit.

For admitted 3..32 actual-real frequencies with maximum absolute half-width
at most0.5, at most17 modes replace the quadratic pair Gram cache. Preparation
and retained storage are O(F*M); material point work is O(F*M*(branches+1)*6),
with fixed point scratch and no query heap allocation. Original wider/count
fallbacks remain. The original spring increment is assembled before
DeltaU=Deltaq*(q_start+Deltaq/2), including equilibrium bias, instead of
subtracting bright stored endpoints. Endpoint reporting remains physical.
These are unweighted physical heat/work intervals; shaped detector means do
not redefine constitutive energy balance.

Validation: the original Kelvin three-tone rate-null regression returned
5.254869241218064e-46 for a centered1e-10s interval whose physical leading heat
is2.53125e-51. The new retained-field result and independently assembled
original-viscosity work agree with81*T^5/320. Eighteen focused actual-source
tests pass (eleven existing and seven new), including original extreme
viscosity/finite work, Maxwell memory below ordinary f64 range with its
original tau, anisotropic real ODE/quadrature closure with equilibrium bias,
precise source clocks, and shared owned/borrowed coefficient parity/bounds.
The shared temporal library and material projection/window test graph compile
with warnings denied. The full native workspace and public TensorProny
integration suite were not run. Finite-mode tail bounds exclude roundoff and
do not certify relative error at arbitrary cancellation nulls.

### Sealed pure-fluid PR coexistence thermal storage

`thermochemical_gas::nonideal::equilibrium::isochoric::PureFluidIsochoricThermalStorage`
connects one original pure PR/NASA source to the existing thermal timeline
using INTERNAL ENERGY at fixed density. u=0 is the prepared cold endpoint.
The original UV solver selects temperature, saturation pressure and phase
amount with a private prepared source-relative target. Every trial retains
original coexistence, volume, source and whole-call work admission, including
four homogeneous caloric evaluations per nonidentical source increment.

Only an interior two-phase subcritical interval is admitted. Original
nonnegative kappa below Tc gives nonpositive attraction derivative, decreasing
p_T(v), nested saturation volumes and nonnegative residual/moving-phase Cv.
The analytic NASA Cv floor supplies the global inverse-temperature bound.
Frozen quality and an isothermal plateau do not replace this response.
Complete original phase/material receipts and source entropy remain available.
Thermal entropy additionally checks direct local U closure and the actual
endpoint reciprocal-temperature interval; unresolved intervals refuse.

Cell admission requires m=rho*V within the original 128-epsilon relative volume
allowance. Fixed density and its identically zero energy derivative use only
coordinate-domain admission, avoiding unnecessary UV inversions; actual
temperature, pressure and phase queries retain the native solve. The implicit
adapter requires explicit bulk k(T), reports phase fraction, and refuses known
conductivity derivative corners for harmonic use.

No second heat integrator, boundary-work debit, phase kinetics,
morphology/transport averaging or open-stream datum is inferred. Five
original-material regressions cover UV responses, sealed heating/cooling,
implicit derivatives and entropy, large formation references, and
source/phase/work/cancellation refusal. They are authored but uncompiled and
unexecuted during the development environment outage.

### Shared finite outflow accounts for receiving-volume consumers

The original reacting reservoir's `step_at_mass` is the single scalar owner of
finite fixed-volume `U=m*h-p*V` export. Both ordinary endpoint publication and
`EquilibriumReservoirDepletion::sample_withdrawal` consume that same account;
scalar coupling trials allocate no product receipts and cannot rebase the
reservoir's immutable birth source or policy. `with_conserved_mass` consumes an
EOS-reconstructed reservoir and admits a prescribed conserved inventory only
within an explicit absolute `rho*V-m` allowance, leaving its T,p and original
source state unchanged. A receiving consumer remains responsible for admitting
the actual returned energy/entropy balances. These operations introduce no
alternative EOS, finite-energy law, species source or material transport model.

## Two-source physical constitutive intervals

The positive tensor-Prony interval path now includes two actual-real source
frequencies through shared admission. Original one-tone and broad-band paths
remain unchanged. The retained numerical payload bound now covers both the
17-mode two-source corner and the pair cache. Nineteen actual material
projection/window tests pass with warnings denied; the added centered1ps
cos(t)-cos(2t) rate verifies physical Kelvin heat9*T^5/320 and independent
original-viscosity stress work. Full public material integration remains unrun.

## Admitted-mass energy coupling during exact-vacuum filling

`EvacuatedGasVolume::fill_candidate_with_energy` retains the original first-fill inventory, composition and caloric owners while accepting one additional internal-energy callback at the exact admitted positive mass. The callback runs once before the sole inversion; its contribution and the supplied source-relative coordinate remain separate receipts. The original `fill_candidate` delegates with zero additional energy.

A composing moving-gas owner can remove retained bulk kinetic energy without assigning a temperature, velocity or seed mass to vacuum, duplicating inventory selection, or performing a second thermodynamic solve. Original source support, extensive mass/energy budgets, iteration limits and cancellation remain binding. The returned material parcel stores thermal internal energy; momentum, opposite storage and work remain the composing owner's responsibility. Callback side effects cannot be rolled back.

Focused regressions specify exact admitted-mass evaluation, independent caloric response, original zero-additional-energy parity, callback refusal and cancellation immediately after callback evaluation.


### Relaxed DC and AC interference in finite material work

Short real material windows now admit one through 32 authored AC frequencies within the existing max(omega)*duration/2 <= 0.5 bound. A singleton reuses the shared real factor with an exact zero-frequency preparation column; only the authored AC coefficients are retained. Authored relaxed DC strain enters the constant temporal mode before the independently supplied constitutive stress is multiplied by strain rate. Thus cancellation between DC and AC stress is retained before physical work, including one-tone and multitone fields. The existing factored spring increment now also handles the singleton short-window case. Original material/source admission, Maxwell zero-bias policy, unweighted physical interval, wide-window pair fallback and retained payload bounds remain in force; point evaluation adds no heap allocation.

Focused actual numerical-source tests reproduce cos(t)-1 over [0, 1e-16]: the previous storage increment was zero and stress work was one-third of the physical 1.25e-65; both now retain that value. Regressions cover signed reverse intervals, multiple AC tones plus DC, an unrepresentable individual stress with finite physical work, singleton payload admission, existing wide-window fallback and prior anisotropic Maxwell/work checks. This focused validation executes the original projection/window/spectrum graph with actual fs-math/fs-fft; the full TensorProny/plane-stress owner integration and workspace were not built in this storage-constrained environment. The finite basis tail bound excludes roundoff and does not certify relative accuracy at arbitrary nulls.

## Original equilibrium factors for mechanical projection consumers

ProjectedTensorProny retains the six original volume-scaled equilibrium Gram
rows already computed by TensorProny::project_strain and exposes them through
equilibrium_gram_rows. They define the same full work-conjugate stiffness as
the existing projected matrix. Original column order, material factor basis,
cross terms and null directions remain intact; no consumer-side factorization
or equilibrium-modulus reconstruction is required.

Projection arithmetic, its scalar-work budget, source law, compatible strain
map, Maxwell factors/times and Kelvin factors are unchanged. The additional
retained rows use six times the admitted coordinate count of f64 storage.
Existing clones and explicitly bound eigenstrain projections retain those
same stiffness factors; affine eigenstrain terms remain separate original
potential terms and are not inferred by this getter.

Verification: The complete spatial-material bridge's four actual-source
groups exercise this getter against original factor rows, full coupled
equilibrium/Maxwell/Kelvin spectra, SI work and native memory continuation.
The bounded run retained the original material projection/source owners and
canonical numerical dependencies; the full material/workspace suite was not run.
