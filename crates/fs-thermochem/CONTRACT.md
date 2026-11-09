# CONTRACT: fs-thermochem

## Channel-resolved fusion-neutron emission

`fusion::neutrons::FusionNeutronEmissionV1` retains the actual integrated
neutron amounts and kinetic energies of D-T and D-D→He3+n separately, with
their accepted interval duration. Each channel converts its amount using the
exact SI Avogadro constant and exposes its mean particle energy only when it
contains neutrons. Empty channels create no fictitious zero-energy particles.
The fixed-size representation checks finite/nonnegative data, consistent
number/energy pairs, positive duration for emission, and representable totals.

`FusionNeutronIntervalV1` retains one complete emission with its original
accepted start/end times. Both endpoints must be finite and nonnegative,
and end must be strictly after start. The emission duration must exactly
match the represented endpoint difference. The immutable, fixed-size record
never moves an endpoint, changes a duration, or rescales number/energy
moments to accommodate clock roundoff. Empty emission can occupy a positive
interval. This supplies accepted source history without selecting an
emission profile inside a step.

This is a common L1 handoff for already accepted source accounts. It creates
neither a neutron spectrum nor a temporal/angular emission profile, and it
does not turn emitted kinetic energy into heat. Consumers must explicitly
choose their projection, transport, deposition and validity assumptions.
Focused G0/G3 regressions cover channel identity, SI count/energy conversion,
empty intervals, invalid domains, arithmetic extremes and extensive scaling.
Interval regressions also cover exact endpoint retention, empty steps, invalid
windows and duration mismatch at ordinary and large absolute clock offsets.

## Original-source isotropic test-ion energy drift and diffusion

`coulomb::test_particle::TestParticleCoulombV1` binds a positive heavy projectile
by its original source index and an explicit kinetic/rest-energy limit. Prepared
backgrounds refresh actual populations, temperatures, density and screening in
fixed-size storage. Source identity, mass, charge, weak coupling, nondegeneracy,
cutoff and thermal support retain the original Coulomb owner's admission.
The raw-point caller owns total hybrid charge and energy closure.

The finite-mass Landau energy moments include drift `A=<dE>/dt`, diffusion
`D=<dE*dE>/(2*dt)`, and velocity drag separately for heavy/electron baths.
Velocity slowing is not mistaken for loss of mean particle kinetic energy.
Stable small-speed series retain finite heating at E=0, where diffusion
vanishes; high-speed energy loss retains the finite target-mass factor.
Each target freezes its Coulomb logarithm at the parent's declared thermal
relative-velocity convention, so `A = D' + D*(1/(2E)-1/(kT_b))` and the
equal-temperature stationary energy density is proportional to
`sqrt(E)*exp(-E/kT)`. Same-species background collisions are included.

Six focused source checks cover independent SI/incomplete-gamma coefficient
anchors, per-bath detailed balance, thermal-tail signs and limits, source refresh,
admission and cancellation. The coefficients follow *NRL Plasma Formulary 2023*,
Collisions and Transport, pp. 31–35:
<https://homepage.physics.uiowa.edu/~ghowes/teach/phys4731/resources/NRL_Plasma_Formulary_2023.pdf>.
This is a dilute, isotropic, nonrelativistic leading-log source law; it owns no
distribution evolution and supplies no neutral stopping, velocity-dependent
screening, quantum collision integral or collective instability model.

## Light-ion nuclear source binding and reaction-weighted product energy

`fusion::plasma::LightIonFusionPlasmaV1` binds original p+, D+, T+, He3++,
He4++ and electron sources in arbitrary canonical order. Isotope associations
are explicit. Bare-nucleus masses, H/He axes, charge, fixed monatomic sensible
calorics and zero formation constants are checked before admission. Both
helium isotopes contribute charge-square four. All potential secondary fuels
and products require original temperature support even when their amount is
zero. Shared source and physical-scale kernels retain the four-source D-T API;
neither binding creates replacement material cards or projects populations.

`FusionReactionV1::thermal_event` is the common allocation-free per-event
energy owner for all four Bosch-Hale channels. It returns the original
reactivity, reaction-weighted incoming ion energy and the two product kinetic
energies. The derivative moment is `Erel=kB*Ti*(3/2+dln<sigma v>/dlnTi)`;
center-of-mass translation contributes `3*kB*Ti/2`. Product energies sum to
Q plus the incoming debit using the existing leading nonrelativistic
inverse-mass relative-energy partition. Common Maxwellian ion temperature,
no relative drift and neglected mass-defect boost corrections remain explicit;
the rate-fit residual is not an error bound for its derivative or energy moment.

Twenty-three native fusion tests pass against the actual source crate,
including six new six-species admission checks and two new all-channel energy
and independent-derivative checks. The binding evolves no state and chooses
neither reaction channels nor charged-product deposition fractions.

## Electron-ion free-free emissive power on original plasma sources

`plasma_radiation::ElectronIonBremsstrahlungV1` evaluates frequency-integrated
electron-ion thermal bremsstrahlung from the existing canonical species bank:
`P/V = 1.426e-40 * ne * sum(ni Zi^2 gi(Te)) * sqrt(Te)` in SI. The coefficient
is the direct CGS-to-SI conversion of [Itoh et al., Eq.(11)](https://arxiv.org/abs/astro-ph/0011044).
Newly populated helium ash contributes charge-square weight four.

`BremsstrahlungGauntFactorV1` explicitly selects `Constant(f64)` or
`ItohNonrelativistic`. The latter uses the eleven literal Table 2 coefficients
in a degree-ten Horner evaluation and its analytic temperature derivative,
separately for each populated positive ion. The published variable is
`gamma^2 = Zi^2 * 1.579e5 K / Te`, with inclusive support `1e-3..=100`.
Thus D/T require 1,579–157,900,000 K and He++ requires 6,316–631,600,000 K,
intersected with the caller's temperature interval and all original source
and physical limits. Endpoint admission precedes conversion to gamma-squared;
there is no clipping, extrapolation or automatic constant-factor fallback.
The published fit residual is not a physical uncertainty or derivative-error
certificate.

Allocation-free point evaluation checks the original populated caloric sources,
canonical amounts, volume, electron identity and mass, electron degeneracy and
relativity, and the neglected ion/electron thermal-variance ratio. These
admission limits are explicit policies. The plasma owner still admits weak
coupling and negligible collective spectral effects.

Results distinguish emissive power, its full temperature derivative and the
instantaneous `P/sqrt(Te)` coefficient, which is temperature independent only
for the constant Gaunt law at fixed populations and volume. Escaped radiation
requires an explicit optically thin consumer closure. This law supplies no
absorption, opacity, electron-electron, line, recombination or synchrotron
radiation. Nine focused native tests pass: the six constant-law checks plus
published-coefficient/derivative anchors, charge-temperature scaling, mixed-ion
fit evaluation, exact per-charge endpoints and newly populated ash refusals.

## Quantum-diffraction Coulomb exchange and changing populations

`CoulombModelV1::new_with_cutoff` retains an explicit classical closest-approach
or quantum-diffraction prescription. The original constructor retains its
classical-scattering refusal. Quantum diffraction uses
`ln(lambda_D/max(b90,hbar/(2*mu*sqrt(S))))` in the same reciprocal Maxwellian
pair moments. The selected fixed thermal convention and static Debye screening
remain explicit; this is a nondegenerate, nonrelativistic leading-log
approximation with the original screening, coupling, occupation, source and
drift gates. It is not a velocity-resolved quantum collision integral.
The cutoff follows the conventional prescription after Eq.3.119 in
[Fitzpatrick's Coulomb Logarithm](https://farside.ph.utexas.edu/teaching/plasma/lectures/node39.html).

`evaluate_point` refreshes canonical trial amounts, volume and separate
temperatures without allocating a frozen thermodynamic state. The retained
original bank supplies every populated row, including newly formed ash.
The caller owns charge balance, caloric state admission and final publication;
existing population/material receipts retain their original binding.

Twelve focused native Coulomb tests pass, including independently reconstructed
hot-plasma SI coefficients, retained classical behavior and physical gates,
fresh ash/density/volume evaluation, and cancellation. These are model and
implementation checks, not a physical-uncertainty certificate.

## Maxwellian nuclear fusion and original plasma thermodynamics

`fusion` evaluates the separate D-T, D-D neutron, D-D proton and D-He3
channels from Bosch and Hale (1992), equation (14), table VII
([primary source](https://doi.org/10.1088/0029-5515/32/4/I07)).
Inputs are common Maxwellian ion temperatures in kelvin; outputs are
particle-pair reactivities in m³/s and analytic temperature derivatives.
The published intervals are inclusive; extrapolation is refused. D-D event
counting includes its one-half identical-pair factor and consumes two deuterons.
Integer nuclear charge and baryon inventories are distinct from chemical
element conservation. Q values use [CODATA 2022 bare-nucleus masses](https://physics.nist.gov/cuu/Constants/Table/allascii.txt).

The product energy split is leading nonrelativistic two-body kinematics.
The fitted-rate derivative supplies the reaction-weighted incoming relative
kinetic-energy moment. Neither fit residuals nor the regression tolerances
certify that moment's physical uncertainty. Mass-defect boost and relativistic
product corrections, screening, beam distributions and transport are outside
this law's scope.

`fusion::plasma::DtFusionPlasmaV1` explicitly binds D+, T+ and He4++ to the
existing source-bearing ionized-gas bank and its electron. Admission requires
the exact four-species charge/element layout, nominal bare-particle masses,
original monatomic calorics with zero formation energy, and source support
including initially absent ash. The caller supplies weak-coupling,
nondegeneracy and nonrelativistic limits. Point admission uses fixed-size
storage; it does not replace a source bank, insert neutron chemistry, or infer
isotopes from material names. Burn evolution belongs to L3.

The focused regression suite contains published table-VIII anchors, SI
conversion, independent derivative differences, source boundaries,
stoichiometry, Q/energy moments, deterministic repetition, and source/domain
binding refusals. These assertions test the declared model, not ignition or a
complete fusion reactor.

## Shared half-life and first-order removal laws

`decay::DecayLawV1` admits a finite positive half-life or a finite nonnegative
inverse-second hazard. Stability is explicit: rate zero and no finite
half-life. The mean lifetime is `1/lambda`, while the half-life is
`ln(2)/lambda`. Independent competing removal channels add their hazards;
this operation does not model consecutive daughter decays or transfer one
channel's heat/particle attribution to another.

The allocation-free fraction closure consumes the canonical `fs-math::decay`
kernel. Survival and loss are evaluated separately, preserving tiny positive
losses when rounded survival is still one, and tiny positive survivors when
rounded loss is already one. Finite nonnegative elapsed time is required,
zero time is exact, and enormous positive exposures saturate to full decay.
Positive rate/exposure underflow and unrepresentable finite lifetimes refuse;
the eventual survivor may round to zero. These are binary64 evaluations,
not interval certificates or independently validated isotope data. All
half-lives, rates and their applicable physical conditions remain caller data.
This L1 closure owns no evolving inventory, decay-chain solver or heat transport.
The decay and activity conventions follow IAEA TECDOC 1741:
<https://www-pub.iaea.org/MTCD/Publications/PDF/TE-1741_web.pdf>.

`ElementaryReactionV1::from_half_life` binds the same hazard into the existing
irreversible elementary ideal-gas path with one first-order reactant, fixed
temperature-independent Arrhenius coefficient and no pressure modifier. The
reactant cannot be regenerated among the products. This is chemical removal:
existing mechanism admission still enforces atoms, charge, source molecular
masses and species thermodynamics. It does not permit nuclear transmutation.
Another reaction feeding that species, or a subsequently attached pressure
modifier, changes the actual inventory's observed halving time. G0/G3 tests
cover one-half survival, mean lifetime, time partition/rescaling, separate
small tails, independent competing removal and invalid/range refusals.

## Silver-sulfide leaching and sulfur cyanide demand

Three matched NBS packs add fourteen literal reference facts for historical
alpha-orthorhombic Ag2S, aqueous thiosulfate and thiocyanate. Conserved
silver-sulfide-to-thiosulfate, thiosulfate-to-thiocyanate and elemental-sulfur-to-
thiocyanate endpoints share the existing gold reaction's finite cyanide and
dissolved-oxygen inventory. Species identities join the earlier pyrrhotite,
sulfite and native-gold owners without a second atom/charge implementation.

The standalone `sulfur_cyanide_demand` example uses explicit stages, preserves
free versus metal-bound versus thiocyanate CN, and prints source masses and
both direct-G and H-minus-T-times-S diagnostics. All stages and source queries
must succeed before a result is printed. Shared reagent overdraw refuses.
Stage ordering is the declared gross-inventory accounting policy; no kinetic
ordering, product distribution, recovery or pH is inferred.

The thiosulfate erratum marker and negative conventional ionic SCN Cp remain
in their original provenance and property family. All claims require 298.15 K
and 100000 Pa with unstated uncertainty. Reference Cp points are not temperature
laws or bulk-solution heat capacities. `mining_sulfur_cyanide_sources` checks
literal data, phase/reference admission, independent H/G/mass calculations,
both stages, all tracked atoms/charge and finite CN/O2 refusal through the
original pack compiler and reaction owners. See
`data/matdb/seed-v1/MINING_SULFUR_CYANIDE.md` for sources, equations and usage.

## Finite neutral-solute saturation

`solution::solubility::IdealSolubilityV1` binds explicit solvent and neutral
solute identities, reference concentration in mol/m³ of incompressible solvent,
reference temperature, reference dissolution enthalpy and an inclusive
temperature domain. Without additional caloric data the enthalpy is constant.
`with_phase_heat_capacities` binds explicit positive solid and dissolved Cp(T)
curves for the same component across the complete original temperature domain.
Their Kirchhoff integral supplies `delta_H(T) = delta_H_ref + integral(delta_Cp)`;
integrating `d ln(c_sat)/dT = delta_H(T)/(R*T*T)` supplies saturation from those
same curves. The ideal-dilute law increases saturation wherever dissolution is
endothermic, and admits a turning point where its enthalpy changes sign.
Logarithmic evaluation preserves
representable results even when a separately formed exponential would overflow;
reference-temperature values and the original zero-enthalpy/equal-Cp limit.
Original phase curves and their merged integration knots share immutable
storage. Evaluations allocate nothing, use stable local analytic integrals and
never infer caloric measurements or phase evidence from a species name.

`equilibrate` retains a finite total component stock and returns separate
dissolved and pure-solid amounts. Dissolved amount is `min(total, c_sat*V)`;
the absent-solid branch may be undersaturated and the present-solid branch
has exactly the supplied saturation capacity. A positive dissolved trace is
retained independently of a much larger solid stock. The compensated component
balance reports floating-point residuals without claiming a certified enclosure.
`partition_at_saturation` is allocation-free and noniterative, allowing an
isothermal material consumer to cache the original pair-law evaluation while
repartitioning at each solvent-volume trial. G0/G3 tests cover finite phase
complementarity, drying/redissolution, temperature direction and slope,
extensive scaling, trace inventories and domain/range refusal.

This is an explicitly declared, independent neutral-solute approximation, not
an electrolyte solubility-product, pH, complexation or common-ion model. It
does not infer measured material parameters, solid-phase selection, solvent
density corrections, nucleation barriers, finite dissolution rates, solid
volume or heat evolution. Concentration is on the stated solvent-volume basis;
molality or total-solution-volume data require an explicit compatible conversion.
The reversible finite-phase inequality follows USGS PHREEQC:
<https://water.usgs.gov/water-resources/software/PHREEQC/documentation/phreeqc3-html/phreeqc3-13.htm>;
its constant-enthalpy temperature treatment is documented in
<https://water.usgs.gov/water-resources/software/PHREEQC/documentation/phreeqc3-html/phreeqc3-36.htm>.
For ionic solids, IUPAC defines the solubility product in terms of ion
activities: <https://goldbook.iupac.org/terms/view/S05742>.

### Chemical work shared by phase exchange and molecular transport

The same `IdealSolubilityV1` exposes the aqueous-minus-solid chemical potential
`RT ln(c/c_sat)` and the solid-relative free energy
`RT n_aqueous [ln(n_aqueous/(V*c_sat))-1]`. The absent-aqueous chemical potential
is explicitly `None` (the negative-infinite ideal limit); the free energy has
its finite zero-inventory limit. No concentration floor or fictitious nuclei
are introduced. Fixed-volume changes retain the saturation reference term:
dissolution changes aqueous moles, so cancelling that term as if only dissolved
inventory were conserved gives the wrong energy for phase-bearing transport.

`solid_relative_free_energy_change_j` integrates the chemical potential using
relative logarithms and a compensated midpoint series for nearby inventories.
It avoids subtracting large absolute energies, retains the quadratic minimum
at saturation, and performs bounded allocation-free work. The common solid
datum cancels only for a closed, isothermal system with identical original pair
laws and separately conserved dissolved-plus-solid totals. These functions do
not supply absolute Gibbs energy, open-system chemical reference fluxes,
mechanical volume work, nonideal activities, or another heat-of-solution source.
G0/G1 tests check the potential derivative, phase-equilibrium minimum, tiny
saturation departures, fixed-area/particle dissipation, Gibbs--Helmholtz
consistency, extensive scaling and domain refusals. No interval enclosure or
universal relative floating-point error bound is claimed.

`solid_relative_inventory_free_energy_change_j` additionally retains actual
solvent-volume changes. A stable fixed-volume inventory difference plus
`RT*n_after*ln(V_before/V_after)` gives the complete isothermal chemical change,
including adjacent-volume dilution. The original solid reference is retained;
the common solvent reference cancels only under closed solvent-volume balance.
Pressure-response work remains a separate mechanical account. Focused checks
recover `dF/dV=-Pi`, path independence, absent-phase limits and the exact
same-volume reduction, without inventing a thermal or pressure-dependent law.

## Finite dissolution and precipitation in one thermal bath

`solution::solubility::caloric::SolubilityCaloricModelV1` closes an isobaric
heat input against one incompressible solvent and 1–64 independently soluble
neutral components. It retains separate dissolved and solid stocks, complete
component totals, exact original solvent volume, and one shared temperature.
An initial state may explicitly contain unequilibrated phase inventories;
zero applied heat then gives adiabatic dissolution or precipitation.

Every component binds its original `IdealSolubilityV1`. Constant-enthalpy
components explicitly declare one common dissolved/solid molar heat-capacity
curve. `SolubilityCaloricComponentV1::from_phase_heat_capacities` instead uses
the exact distinct phase curves already held by the saturation law, preventing
a second caloric source from changing its thermodynamics. Solvent and actual
phase stocks each contribute sensible heat once; changing dissolved amount
contributes the original phase enthalpy at the trial temperature once.
The shared `MolarHeatCapacityCurveV1` extracts
the existing exact piecewise-linear Cp integration without implying a phase
qualification. `PureLiquidHeatCapacityV1` preserves its liquid-specific meaning.

The closed first-law increment integrates Cp directly from the initial to trial
temperature, weights the solid/dissolved sensible terms by the initial phase
stocks, and uses `delta_n_dissolved * delta_H(T_final)` for phase transfer.
It never subtracts two
large reference enthalpies. Equilibrium enthalpy is continuous at solid
exhaustion and strictly increases with temperature: positive sensible capacity
is augmented by n_dissolved*(delta_H(T))^2/(R*T^2) for each saturated component.
One bounded scalar bisection therefore resolves the common bath without an
inner chemical solve. Trial evaluations allocate nothing and construct accepted
phase arrays only once. The original pair and Cp domains, explicit temperature
bracket, energy/temperature tolerances and iteration budget remain binding;
out-of-domain heat and late cancellation refuse without changing input state.

`prepare_at_solvent_volume` preserves a material consumer's exact volume
representation while deriving its solvent amount once from the declared molar
volume. Amount-based preparation remains available. Original source identity,
ordered totals, phase partitions, solvent amount/volume and direct energy
residuals remain inspectable. Independent analytic cases cover both signs of
dissolution enthalpy, precipitation, opposed phase transfers into one bath,
solid exhaustion/recrystallization, distinct phase Cp, the independent equilibrium
enthalpy derivative, variable Cp, small heat increments, volume roundoff,
work limits and immutable refusal.

This is the thermal closure of the stated ideal neutral-solute approximation.
It supplies no electrolyte activities, common-ion chemistry, pressure-dependent
solubility, variable-pressure work, separate phase temperatures, nucleation or
finite dissolution kinetics. Numeric caller coefficients do not authenticate
measured material data.

## Finite-rate neutral dissolution and crystal growth

`solution::solubility::kinetics::SolubilityKineticsV1` uses the original pair law
and an explicitly declared fixed-area conductance `G=kA [m³/s]`. At fixed
temperature and solvent volume its active equation is
`dn_dissolved/dt = G*(c_sat - n_dissolved/V)`. The analytical exponential update
retains actual dissolved and solid stocks, their conserved component total,
signed transfer and a separate equilibrium target. A finite solid stock is
exhausted at its analytical event; the remaining interval has no dissolution.
It does not clip an Euler step or substitute the equilibrium target for the
accepted kinetic state.

`ExistingSolidOnly` permits precipitation only while actual solid is present.
`SeededSurface` explicitly declares a persistent growth substrate with no added
solute moles; it supplies neither a nucleation barrier nor an induction-time
law. Dilution or temperature changes can therefore leave a finite lag, and an
unseeded supersaturated solution may remain supersaturated. The fixed-area
well-mixed approximation is derived in the module from the Noyes–Whitney law:
<https://doi.org/10.1021/ja02086a003>.

Scalar updates allocate nothing. Stable exponential increments retain short
steps; explicit phase storage retains small remaining donors. An unrepresentable
asymptotic tail may round to zero but is never labeled a finite exhaustion event.
G0/G1/G3 cases check independent exponential/event solutions, semigroup and
extensive scaling, zero contact/time, seeded/unseeded growth, trace stocks and
source/domain refusal. Its conductance and experimental applicability remain
explicit caller inputs; reaction/speciation and thermal evolution remain
separate owners.

`with_particle_surface(reference_solid)` selects a fixed-number, self-similar
particle population with constant density and film coefficient. The original
conductance applies at that reference amount; the actual exposed area gives
`G(s)=G_ref*(s/s_ref)^(2/3)`. The phase step integrates the closed finite-bath
radius equation `r'=-G_ref*(r^3+V*c_sat-n_total)/(3*V*s_ref^(2/3))`, where
`r=s^(1/3)`, with bounded allocation-free scalar work. Area is not frozen at
the beginning of the step. At `V*c_sat=n_total` the exact algebraic tail replaces
the fixed-area exponential. Supersaturation grows existing particles; a finite
exhaustion event removes their area and no later step silently renucleates them.
`SeededSurface` is refused for this geometry because an inert persistent seed
requires a different area/volume law. `saturation_deficit_half_life` refuses
active particle kinetics, which have no constant exponential half-life.

The original reference stock survives every state copy and transport/thermal
handoff. `conductance_at` exposes the actual area dependence. G0/G1/G3 particle
tests compare independent radius-ODE and exhaustion-integral oracles, exact
algebraic decay, split-step composition, amount/volume rescaling, trace phases,
growth, no-nucleation and invalid geometry. The Hixson–Crowell self-similar
premise is explicit: <https://doi.org/10.1021/ie50260a018>. Breakup, aggregation,
polydispersity, changing film thickness and solute-volume displacement are not
modeled. No performance or experimental calibration claim is inferred from
the numeric reference coefficients.

`advance_at_saturation` lets an explicit coupled thermodynamic owner supply a
conditional surface concentration while reusing these same exact fixed-area
and bounded evolving-particle trajectories. It preserves the original source,
temperature support, rate coefficient, particle reference and growth policy.
Its scalar stage holds the supplied concentration fixed; the enclosing owner
must qualify that activity/source calculation and control changes in other
species. Returned target/exhaustion receipts use the supplied conditional
saturation. This numerical entry point supplies no nonideal calibration,
excess calorics or new nuclei. Focused checks establish identity with the same
native trajectories at each supplied target and retain source/domain refusals.

## Coupled finite-rate dissolution and thermal evolution

`solution::solubility::caloric::kinetics::SolubilityThermalKineticsV1` binds
the original caloric model and the exact ordered kinetic laws. Actual dissolved
and solid stocks, the original solvent volume, temperature and source owner
travel together. Different conductances, particle reference stocks, growth
policies or phase heat-capacity sources cannot silently replace the declared
physical model.

Without an interaction extension, each stage evaluates the original analytical fixed-temperature kinetic law at
the unknown endpoint temperature. One bounded scalar heat balance solves that
temperature using the shared sensible-plus-dissolution enthalpy equation.
On an active branch the temperature derivative is positive phase-weighted Cp
plus `w*n_sat*delta_H(T)^2/(R*T^2)`. Fixed area has `w=1-exp(-G*h/V)`;
particle area instead uses the nonnegative sensitivity of final dissolved
amount to saturation capacity. Its scalar sensitivity ODE has positive forcing
`G(s)/V` and zero initial value. Finite-stock and growth-surface branches preserve
the same monotone energy closure without pretending particle area stays fixed.

A full-stage/two-half-stage comparison controls changing-temperature evolution
and accepts the two-half-stage state. This is a first-order approximation to
continuous thermal kinetics, not an exact continuous-event solution or a
certified global error bound. Dissolved and solid errors are scaled separately.
For `advance`, the original window heat is delivered uniformly in time; each stage is anchored
to its cumulative heat target and the final energy residual is independently
recomputed from original/final actual phases. Stage-root tolerance does not
silently accumulate into an enlarged heat allowance.

Trial storage is bounded by the existing 64-component limit and allocates no
heap memory. Only returned public states allocate. Input states remain immutable
on cancellation, domain refusal, unbracketed heat, temporal tolerance failure or
work exhaustion. Independent cases cover closed-surface sensible heating,
zero-enthalpy kinetics, both heat signs, finite-rate/equilibrium limits,
first-order convergence against an independent ODE, unequal phase Cp, shared
bath balances and transactional refusals. The model retains fixed solvent
volume, reference conductance and particle geometry law over each window;
actual selected particle area evolves in each phase stage. Simultaneous spatial
transport, evolving film coefficients and nucleation physics remain separate.

### Temperature-qualified interactions on the same thermal owner

`SolubilityThermalKineticsV1::with_interactions` binds one explicit immutable
`SolubilityThermalInteractionV1` source to the original caloric/kinetic owner.
The source supplies conditional saturation and the additional enthalpy increment
from its SAME thermodynamic potential, including temperature-dependent pure-solid
calibration. Core thermochemistry still calls each original `advance_at_saturation`
and owns all phase states, particle laws, heat roots and adaptive stepping.
Material interactions cannot replace kinetic laws or publish another phase state.
The exact source Arc travels with states and outcomes; source replacement and
cross-owner energy/advancement refuse. Temperature support is the intersection
of the original caloric and explicitly supplied interaction domains.

For interacting phases, one declared-order phase stage uses the current
temperature, then the existing scalar heat solve projects temperature with
the resulting stocks frozen. Source binding must qualify positive frozen-stock
heat capacity throughout its domain. Composition convexity alone does not
establish this thermal condition, and the independent ideal endpoint-phase
monotonicity argument above is not asserted for interacting species. The same
full/two-half controller measures the complete first-order splitting error.
Native trace stocks, particle reference geometry, cumulative heat targets and
direct original-to-final energy checks remain unchanged. Trial evaluations
allocate nothing and use the original bounded work/cancellation limits.

`MolarHeatCapacityCurveV1::heat_capacity_bounds` evaluates exact extrema of the
original piecewise-linear curve on a supported interval, including interior
source knots. `IdealSolubilityV1::log_saturation_ratio` and
`dissolution_enthalpy_change_j_mol` use the original integrated van't Hoff and
Kirchhoff kernels directly, preserving small thermal changes that subtraction
of rounded concentrations or large enthalpy datums would erase. Focused native
cases cover that neighboring-temperature limit, original-pair phase/heat
splitting, exact Cp extrema, source identity and thermal-domain refusal.

### Finite-conductance heating and cooling of soluble phases

`SolubilityThermalKineticsV1::advance_with_heat_transfer` couples the SAME
thermal/phase owner to a prescribed-temperature reservoir through a complete
contact conductance `G [W/K]`. The boundary law is Newton heat transfer,
`power_into_solution=G*(T_reservoir-T_solution)`; geometry and transfer
coefficients must already be included in G. This is distinct from the molecular
film conductance `kA [m³/s]`. The standard boundary law is described in the
[NASA Passive Thermal Control Engineering Guidebook](https://ntrs.nasa.gov/citations/20220006584).

Each stage solves the implicit endpoint heat `Q_stage=G*h*(T_reservoir-T_end)`
inside the original bounded scalar root. Its residual derivative adds `G*h`
to the previously qualified positive caloric derivative. The original adaptive
full/two-half comparison controls the phase and boundary-feedback error together.
Only the accepted two-half-stage boundary heats enter the compensated heat sum;
coarse or rejected trials contribute nothing. Direct original-to-final enthalpy
is checked against this independently accumulated supplied heat, so the reported
energy residual is not made zero by defining heat from the solution enthalpy.

Zero G uses the original adiabatic path exactly. Reservoir temperature must be
positive and finite but need not lie inside the solution's source interval;
every solution endpoint still obeys its original temperature bracket. Domain,
work-limit, cancellation and publication rules remain unchanged. The reservoir
is prescribed and has no evolving finite heat capacity in this entry point.
Focused cases compare sensible cooling against the independent exponential
solution, verify first-order refinement and actual boundary heat, and exercise
particle growth, nonideal source retention and transactional refusals.

### Finite sensible heat reservoirs on the original thermal timeline

`SensibleHeatReservoirV1` declares a positive constant complete heat capacity
`C [J/K]` and an original positive reference temperature. It retains a
compensated energy inventory relative to that reference. Published temperature
is derived from this inventory; continuing a window must pass the returned
reservoir, not reconstruct its energy from a rounded endpoint temperature.
The constant-C constitutive law is explicit and supplies no material-dependent
heat-capacity curve, reservoir phase change or mechanical work.

`SolubilityThermalKineticsV1::advance_with_finite_heat_reservoir` uses the SAME
thermal stages and adaptive timeline. With `kappa=G*h`, elimination of the
reservoir endpoint gives the implicit heat into the solution as
`Q_stage=kappa/(1+kappa/C)*(T_reservoir,start-T_solution,end)`.
The reservoir receives exactly the opposite accepted heat. Both capacities
participate in the transient; finite storage is not approximated by holding
the starting reservoir temperature fixed over the window. The effective
contact term is nonnegative, preserving the monotone scalar caloric root.

The original full/two-half error includes the reservoir temperature as well
as solution temperature and both phase inventories. Rejected/coarse exchanges
never enter either accepted energy store. The returned
`SolubilityThermalReservoirOutcomeV1` retains the original thermal balance,
final reservoir, its direct energy increment and the independent closed-system
energy residual. Repeated windows preserve the original reservoir datum and
small heat increments. All original solution/source/work/cancellation limits
remain binding; finite reservoir storage introduces no second phase solver.

## Binary liquid miscibility and source-owned thermodynamic forces

`solution::BinaryRegularSolutionV1` retains two explicit ordered species and a
caller-declared temperature domain for the default `chi(T)=A+B/T`. Its single symmetric
molar free energy supplies activity coefficients, chemical potentials relative
to pure-liquid standards, mixing Gibbs energy/enthalpy/entropy and the
thermodynamic factor `Gamma=1-2 chi x(1-x)`. Pure states have zero mixing energy;
absent-species chemical potentials are `None`, while their Gibbs–Duhem-balanced
force products use the finite one-sided limit. Trace composition keeps the
majority species' `log(1-x)` contribution even when `1-x` rounds to one.

For `chi>2`, a bounded allocation-free log-odds bisection excludes the false
central root and resolves the common-tangent binodal, with a small-argument
series near criticality. It distinguishes that global coexistence gap from
the local spinodal. A finite closed inventory uses the molar lever rule and
reports one or two phases, plus stable/binodal/metastable/spinodal/critical feed
classification. Species mass or phase volume must use original molecular masses
and density laws in the consumer. The final root bracket is numerical evidence,
not a rigorous interval certificate. Unrepresentable rich/dilute coexistence
and exhausted work budgets refuse; cancellation never publishes a partial result.

`BinarySolutionStateV1::driving_forces` supplies isothermal/isobaric
`[Gamma grad(x_A), -Gamma grad(x_A)]` to both existing Maxwell–Stefan frames.
Negative `Gamma` identifies demixing; it is not a molecular diffusivity for a
Fick evolution scheme. Tests verify independent energy/chemical/caloric
derivatives, pure and trace limits, equal chemical potentials, component
conservation, critical limits, temperature-driven mixing/demixing, species
exchange and original Maxwell–Stefan flux/drag-entropy equations.

`evaluate_composition` and `equilibrate_composition` accept both independently
represented mole fractions, with an explicit `8*f64::EPSILON` sum-roundoff
budget and common normalization. A dilute component is retained even when its
dominant partner rounds to one. Binodal comparisons and phase amounts use the
actual minority distances `x_B-x_binodal` and `x_A-x_binodal`, preserving a
finite minority phase that subtraction of rich-side coordinates would erase.
The scalar APIs delegate and retain their original pure/binodal endpoint
conventions. Paired-composition regressions cover the unequal-inventory trace
case, both species orders, common normalization, critical/pure/binodal states,
finite dilute phase amounts, and component conservation.

The constitutive model does not infer empirical pair parameters, pressure
dependence, electrolyte or unequal-size corrections, interfacial tension,
nucleation kinetics or a spatial phase-field evolution. Its free-energy form
is the molar equivalent of the INL MOOSE regular-solution implementation:
<https://mooseframework.inl.gov/docs/doxygen/modules/classRegularSolutionFreeEnergy.html>.

`with_excess_heat_capacity(C,T0)` adds an explicitly authored, signed constant
molar excess-Cp coefficient. With `q=x_A*x_B`, the same interaction supplies
`H_ex=q*[R*B+C*(T-T0)]`, `S_ex=q*[-R*A+C*ln(T/T0)]` and `G_ex=H_ex-T*S_ex`.
Thus homogeneous excess Cp is `C*q`; activity, phase equilibrium and spatial
composition heat/entropy all follow that same free energy. The temperature
remainder and local entropy increments retain nearby temperatures without
subtracting large reference values. Zero C preserves the original law and
its owner identity. Signed C may produce closed/reentrant miscibility windows;
no phase solver assumes a single upper critical temperature.

### Multicomponent homogeneous liquid miscibility and diffusion

`solution::multicomponent::MulticomponentRegularSolutionV1` extends the liquid
constitutive path to 1..=128 explicitly identified components. Two complete,
symmetric, zero-diagonal source matrices define `chi_ij(T)=A_ij+B_ij/T`.
Canonicalization moves both matrices with their original species; missing
interactions are refused, while explicit zero coefficients admit ideal pairs.
The same molar free energy `RT sum(x_i ln x_i + sum_{j>i} chi_ij x_i x_j)`
supplies log activities, chemical potentials relative to pure-liquid standards,
mixing enthalpy/entropy and partial molar excess enthalpies. Exact absences
have no finite chemical potential. Independent positive trace fractions are
retained, including the carrier's small logarithm and activity correction.

The source-owned state applies the FULL coupled thermodynamic-factor operator
to a zero-sum composition gradient. Antisymmetric pair accumulation preserves
Gibbs-Duhem balance to roundoff; a species can respond even when its own mole
fraction gradient is zero. `mass_diffusion` and `molar_diffusion` pass these
forces to the existing Maxwell-Stefan solve only after matching original
species identities and canonical order. Liquid molar concentration and positive
pair diffusivities remain independently supplied sources. Neither is inferred
from an ideal-gas law or from a negative thermodynamic curvature.

Directional simplex curvature reports local homogeneous demixing directions;
positive curvature in one direction is not a global stability certificate.
The response retains its immutable original owner and has no post-construction
heap allocation beyond shared-owner reference counting. Energy/force evaluation
is O(n²), with fixed-size bounded scratch and cancellation at pair-row boundaries.
Tests cover extensive Gibbs/enthalpy derivatives, Gibbs-Duhem and Euler
identities, binary reduction, ternary cross diffusion and cosolvent effects,
both transport frames and friction entropy, pure/trace/maximum-size limits,
canonical permutations, source refusals and atomic cancellation.

The homogeneous response does not perform a phase search; the separate
`equilibrium` module supplies the seeded local candidates described below.
Neither path infers equilibrium phase count, polymer/electrolyte entropy,
measured interaction parameters, pressure/Soret forces or an N-component
spatial phase field. Pure sensible enthalpies must accompany
the partial excess enthalpies in a full diffusive-energy account. The model is
the symmetric-pair specialization of the multicomponent Margules formulation:
<https://cantera.org/stable/cxx/d7/dfe/classCantera_1_1MargulesVPSSTP.html>.

### Complete local miscibility stability and tangent planes

`solution::multicomponent::stability` examines the entire present-species
composition simplex. The entropy congruence
`I + diag(sqrt(x)) chi diag(sqrt(x))`, restricted perpendicular to `sqrt(x)`,
preserves curvature signs while avoiding an overflowing `1/x` matrix.
A deterministic Householder basis and bounded cyclic Jacobi solve report the
sorted spectrum, a smallest-mode composition tangent, negative-mode count,
completed work and numerical uncertainty. Eigenvalues use the ideal-entropy
metric; the returned tangent has unit `sum(d_i^2/x_i)`. Explicit scratch and
sweep budgets admit O(N²) storage and bounded O(sweeps*N³) work. Cancellation
or exhausted work returns no partial assessment. Pure supports have no modes.

Positive local curvature does not rule out metastability. The separate
`tangent_plane_distance` evaluates a finite trial using the original model,
common temperature and closed species support. A negative distance is an
instability witness; finitely many nonnegative trials are not a global
miscibility certificate. The local numerical uncertainty is not an interval
proof. No concentration floor, absent-species inventory, interaction fit or
phase-count inference is introduced. Directional curvature also forms its
physical ideal-entropy norm before interaction rescaling, retaining finite
curvatures for subnormal trace fractions.

G0/G1/G3/G4 tests cover binary analytic spinodals, collective ternary instability
despite positive pairwise exchanges, a locally stable but globally demixing
ternary witness, all 128 components with a subnormal trace, pure/absent supports,
relative-entropy tangent planes, species relabeling and atomic budget/cancellation
refusal. The thermodynamic tangent-plane criterion follows Michelsen (1982):
<https://doi.org/10.1016/0378-3812(82)85001-2>.

### Finite-inventory multicomponent liquid phase candidates

`solution::multicomponent::equilibrium` searches an explicitly supplied phase
count and positive seed on the feed's exact species support. Species-column
allocations are normalized to the original finite inventory; each largest
allocation is reconstructed so dilute allocations are retained independently.
A reusable, entropy-scaled analytic extensive-Gibbs Hessian drives damped
Newton steps with positive-definite regularization, positivity-limited steps
and bounded energy-decreasing backtracking. With d=(phases-1)*present_species,
dense work is O(d³) and numerical scratch is admitted before allocation.
Iteration, regularization, line-search, phase-count and memory limits are explicit.

An accepted candidate retains the original finite feed, complete model, common
temperature, every phase's actual component moles and thermodynamic state, and
its controls and numerical report. Final independent checks require relative
per-species balance, agreement of every present-species chemical potential,
distinct compositions, positive phase amounts, and locally stable phase
curvature. The total Gibbs energy cannot exceed either the homogeneous feed
or feasible initial seed beyond ONE original final allowance. The report
records actual residuals, extensive feed/seed/final mixing energies and work.
Cancellation or any failed check returns no partially optimized inventory.

This is a seeded local candidate, not automatic phase-count discovery or a
global Gibbs-minimum certificate. Coincident, disappearing, unstable and
numerically unresolved phases refuse; none are silently merged, removed or
assigned a concentration floor. No interfacial energy, nucleation kinetics,
spatial evolution, volume law or averaged transport property is inferred.
Pure-liquid standard terms cancel because each species is conserved at the
same temperature and caller-qualified pressure.

G0/G1 tests include the existing independent binary coexistence oracle,
symmetric three-phase outer roots and asymmetric lever fractions, a metastable
ternary feed with positive local curvature but lower three-phase free energy,
unequal pair interactions, finite traces, exact absent species, pure feeds,
spurious stationary/coincident roots, bounded work and atomic cancellation.

### Shared safeguarded caloric temperature inverse

`solution::caloric::temperature` exposes the existing bracketed inverse through
`TemperatureControl`, `TemperatureSample`, `TemperatureRoot` and `solve`.
Original binary, multicomponent and finite-solid callers retain their own
energy evaluator, source admission, cancellation and error type. The shared
owner validates finite positive increasing temperature bounds, positive
temperature tolerance, nonnegative energy tolerance and 1..256 iterations.
An admitted positive capacity supplies a safeguarded Newton proposal; an
unusable optional tangent falls back to bisection without replacing the
physical energy residual. This numerical inverse does not itself certify a
material branch, positive capacity between samples or global uniqueness.

### Homogeneous multicomponent liquid calorics

`solution::multicomponent::caloric` binds every original species to one
`PureLiquidHeatCapacityV1`, preserving the source knots, common sensible datum
and common finite temperature domain. Explicit homogeneous inventories retain
all component amounts, original mixing state, total enthalpy and exact isobaric
capacity. Unequal-temperature mixing and signed heating/cooling solve the
conserved total enthalpy with the same safeguarded temperature inverse used by
the binary caloric owner. Original source Cp and Cp/T integrals handle every
crossed knot; arbitrary sensible-reference offsets cancel before evaluation.

The original A+B/T potential supplies nonideal enthalpy and entropy exactly
once. At fixed homogeneous composition these terms have no temperature
derivative, so each thermal trial evaluates only bounded pure-source integrals
and the positive exact capacity tangent. Composition and pair interactions are
prepared outside the inverse; thermal trials allocate no phase or composition
state. Ideal mixing entropy uses a compensated relative-entropy formulation
that retains trace and dominant-species contributions without species floors.

The caller supplies finite root, energy, entropy and iteration budgets.
Source mismatch, unsupported temperature, unrepresentable component addition,
unbracketed heat, exhausted work and cancellation return no endpoint. A closed
zero-heat transaction also refuses if its homogeneous endpoint decreases
entropy beyond the explicit numerical tolerance. With prescribed nonzero
heat, reported system entropy change requires the actual reservoir exchange
before it can be interpreted as entropy production. This owner neither
asserts equilibrium miscibility nor evolves phase splits, vapor, solids,
nucleation or pressure changes; pressure and pure-liquid source qualification
remain caller responsibilities.

G1 regressions cover ternary unequal-temperature mixing with all pair H/S
terms, invariance to large sensible-reference offsets, source Cp-knot heating
and cooling, reduction to the original binary caloric owner, all 128 species
with exact absences and a retained dilute inventory, an entropy-decreasing
constrained homogeneous endpoint, source identity and atomic refusal.

`MulticomponentCaloricStateV1::reduced_entropy_j_per_k` retains the original
pure Cp/T integrals and homogeneous mixing entropy while analytically removing
one shared set of linear pure-species constants at the original datum. This
coordinate supports different stationary stream compositions; it does not
infer absolute entropy or pressure dependence. Direct same-inventory entropy
increments remain the stable owner for nearby thermal states.

### Multicomponent local equilibrium heat capacity

`solution::multicomponent::equilibrium::caloric` differentiates the SAME
accepted phase candidate at fixed pressure and conserved species inventory.
Explicitly ordered original `PureLiquidHeatCapacityV1` curves must share one
datum and support the actual temperature. Their sum `sum(n_i Cp_i)` supplies
the frozen capacity. For the retained `chi_ij=A_ij+B_ij/T` law, fixed-composition
excess enthalpy has zero temperature derivative; additional capacity comes
from changing phase amounts and compositions.

The existing entropy-scaled extensive-Gibbs Hessian and original partial
molar excess enthalpies supply that redistribution analytically. One
unregularized Cholesky solve yields the local thermal response, its scaled
linear residual and a nonnegative quadratic-form heat-capacity contribution.
Log-temperature scaling avoids forming T squared; a scaled Cholesky norm
avoids cancellation between opposing species heat contributions. No perturbed
temperature flash or independently fitted mixture Cp is used. Each species'
anchor derivative is reconstructed from its paired transfers, preserving
the conserved inventory and exact absent species.

The borrowed response retains the original candidate, ordered pure sources,
common temperature support, total/frozen/redistribution capacities, and every
phase's component and total mole derivatives. This remains a smooth, seeded
LOCAL fixed-phase-count response with the original candidate's numerical
accuracy. Singular or unresolved joint curvature refuses even if every
individual phase is locally stable; regularization is never substituted for
physical Cp. It does not certify global equilibrium, cross phase onset or
supply an enthalpy inverse. Scratch is admitted before allocation, work is
O(d cubed) for d=(phases-1)*present_species, and cancellation is atomic.

Focused native tests compare binary closed-form and independent symmetric
ternary common-tangent responses, independently differenced equilibrium
enthalpy/entropy and changing phase inventories, extensive scaling, species
permutations, exact absences and positive traces. Source, budget, cancellation
and genuinely unresolved near-critical joint-response cases refuse explicitly.

### Finite multicomponent equilibrium heat exchange

`solution::multicomponent::caloric::equilibrium` binds the SAME accepted phase
candidate to the existing `MulticomponentLiquidCaloricModelV1` and original
pure Cp sources. `add_heat_equilibrium` applies signed heat at fixed pressure
and conserved species inventory while phase amounts and compositions adjust.
The existing safeguarded temperature inverse uses analytic equilibrium Cp as
its Newton tangent and accepts the actual source-relative enthalpy residual.
Both the temperature-width and energy budgets remain explicit.

Each bounded temperature continuation step predicts the original phase
component amounts from their analytic derivatives, reconstructs the largest
species allocation for conservation, and corrects with the original equilibrium
solver. Every corrected state must retain the phase count/support, remain in
the local predictor neighborhood, and admit the unregularized physical Cp.
The caller bounds temperature steps and total phase solves; each phase solve
retains the original candidate's inner controls. Disappearing, coincident or
unresolved phases refuse rather than changing the physical model implicitly.
For nonzero heat, the complete offered thermal bracket must be reachable;
zero heat within the bracket re-admits the original state without phase solves.
This is seeded
numerical continuation, not a mathematical branch-identity proof, global Gibbs
envelope or a phase-appearance/latent-transition solver.

Direct H/S increments integrate original pure Cp and Cp/T and difference the
actual phase amounts and pair-composition products locally before extensive
scaling. Large pure-reference totals are never subtracted to recover heat.
The result retains the exact final candidate, total/frozen/redistribution Cp,
signed heat, entropy change, caloric residual, bracket and phase-solve work.
Finite entropy must satisfy the endpoint Clausius bounds within the declared
entropy, energy and roundoff allowance; external reservoir entropy remains a
separate account. Source/branch/work failure and cancellation publish no endpoint.

`equilibrium_at_temperature` exposes that same bounded branch continuation as
a caloric query for native spatial heat storage. It retains the original feed,
phase count, support, predictor, inner candidate controls and unregularized Cp.
Its actual source H/S increments must agree with the temperature direction and
satisfy the same finite Clausius bounds. Querying a temperature does not infer
an external heat source, a globally stable Gibbs envelope or a latent interval.
Material thermal consumers retain the original pure-source capacity minimum
and analytic redistribution response rather than sampling another caloric law.

`MulticomponentEquilibriumCaloricStateV1::reduced_entropy_j_per_k` exposes
the same pure-species entropy convention as the homogeneous owner. Conserved
pure Cp/T integrals are combined with each actual phase's mixing entropy and
actual phase amount. A homogeneous feed's mixing entropy never substitutes for
that separated inventory. Stationary species balance cancels the omitted
common linear species constants, permitting native open-flow entropy accounts
across homogeneous and separated states without an absolute-entropy claim.
Small same-inventory changes continue to use the direct original H/S increment
methods, avoiding subtraction of anchored totals.

### Original multicomponent latent storage at an offered common tangent

`solution::multicomponent::caloric::equilibrium::coexistence` admits two
original local phase ensembles from the SAME complete caloric owner, exactly
shared temperature and conserved species feed. Every present component must
have the same original reduced chemical potential in every retained phase,
within the declared numerical allowance. The caller bounds retained phases
and funds separate chemical-potential and entropy-closure tolerances; the
actual maximum chemical-potential residual remains visible.

Original phase mixing changes through the same homogeneous feed recover
positive latent H and S before extensive reference datums enter arithmetic.
The actual residual DeltaS-DeltaH/T is retained and must fit the stated entropy
allowance; no source entropy is fitted to the heat. `state_at_enthalpy_increment`
uses the actual latent lever rule and returns a borrowed inventory of both
ensembles' original phases, with conserved species amounts and immutable
original compositions and chemical potentials. Endpoint-zero ensembles are
absent; a positive amount that cannot be represented refuses. The represented
inventory's heat residual remains explicit. Plateau queries need no nonlinear
phase solve, numerical derivative or artificial heat capacity.

This is finite numerical admission of an offered common tangent, not automatic
transition-temperature discovery, a global Gibbs-envelope proof or a kinetic
or morphology model. Foreign owners, unmatched potentials, nonpositive latent
storage, exhausted phase limits and cancellation publish no partial state.
Focused regressions use the independent symmetric ternary three-phase/one-phase
common tangent, its four-phase interior inventory, exact source latent H/S,
large pure-reference cancellation, control limits and transactional refusals.

### Conserved mixing and relaxation into local liquid phase branches

`solution::multicomponent::caloric::equilibrium::transactions` extends the SAME
finite-heat owner to actual homogeneous or separated input parcels.
`relax_homogeneous` conserves a homogeneous parcel's full stored enthalpy while
it enters an explicitly admitted phase branch; at zero supplied heat, demixing
can therefore warm or cool the parcel. `mix_equilibrium` combines two actual
input inventories at possibly different temperatures, retaining each original
phase's amounts, composition and mixing storage. A caller-supplied branch
anchor must belong to the same complete caloric owner and contain exactly the
combined species inventory. It supplies numerical continuation evidence, never
additional matter, heat or a global phase-selection claim.

The original analytic predictor/corrector, fixed phase count, whole-call work
limits, unregularized physical Cp and safeguarded temperature inverse remain
the only solve path. Actual input-to-final Cp integrals enter the heat residual
directly; requested heat is not first rounded into a large artificial branch
heat offset. Original pair-product differences and relative ideal mixing
entropy are evaluated through the same conserved homogeneous feed, which
cancels between input and final phase inventories. Every actual phase amount
weights its original excess H/S once. Pure standard-state datums cancel before
arithmetic, and no apparent heat capacity or latent-heat source is added.

The accepted balance reports total input-to-final enthalpy and system entropy.
An isolated transaction with decreasing entropy beyond the caller's allowance
refuses, including an attempt to force a higher-Gibbs local branch. For nonzero
heat, the actual reservoir entropy remains external. Source/owner or inventory
mismatch, unsupported phase continuation, exhausted work and cancellation
publish no partial result. These are bulk, fixed-pressure local branch
transactions; phase onset/disappearance, interfacial energy and kinetics remain
separate. Focused regressions cover the binary adiabatic-demixing oracle,
arbitrary reference shifts, unequal-temperature ternary mixing with homogeneous
and separated inputs, distinct enthalpies and unequal weights of asymmetric
solvent phases, the isolated entropy guard and transactional refusals.

### Conserved thermal endpoints with bounded liquid phase discovery

`solution::multicomponent::caloric::equilibrium::discovery` adds
`add_heat_discovered` and `mix_discovered` for actual homogeneous or separated
parcels. Every temperature sample invokes the original deterministic phase
discovery with exactly the same configured seed policy, then binds its selected
candidate through the original caloric owner and unregularized physical Cp.
The original shared temperature inverse accepts the full source-relative
enthalpy residual; the final phase count can differ from every input's count.
Pure Cp/CpT integrals and the existing common-feed phase H/S differences count
every actual initial and final phase once, including unequal input temperatures.

The result retains the original final discovery report and controls, accepted
caloric candidate, finite thermal/entropy balance, temperature bracket, actual
search count and total candidate attempts. Each new sample reserves its full
configured search allowance before work begins; completed searches debit only
actual attempts. Exhausting whole-call work never silently shortens the seed
policy. Isolated endpoints must pass the original entropy-increase admission.
Nonzero supplied heat still needs the caller's actual reservoir entropy account.

This is a bounded lowest-accepted-Gibbs search, not a proved global envelope or
a continuous phase trajectory. A failed search never falls back to an unstable
homogeneous feed. An unresolved physical response or enthalpy jump refuses;
no interpolation manufactures a latent plateau, extra coexistence phases,
interfacial storage or nucleation kinetics. Numerical and source failures and
cancellation leave all inputs intact. Focused physical regressions cover binary
2-to-1 heating and 1-to-2 cooling against original common-tangent calorics,
large reference shifts, adiabatic demixing and unequal-temperature mixing,
three-phase disappearance with an absent species, metastable entropy refusal
and the symmetric ternary first-order enthalpy gap.

### Closed isobaric liquid calorics

`solution::caloric` combines that SAME pair law with explicitly ordered
pure-liquid molar Cp(T) curves. `PureLiquidHeatCapacityV1` admits at most 4096
positive, strictly ordered source points and exactly integrates their
piecewise-linear interpolation, including crossed knots. Cooling negates the
same ascending integration. No gas standard-state substitution, total-mixture
Cp average, Cp(T) times a finite temperature change, or fabricated property
fit is used. The caller qualifies the pure-liquid sources and their common
pressure; the numerical model preserves species order, common datum and the
intersection of the source temperature domains.

`BinaryLiquidCaloricModelV1` prepares finite equilibrated or explicitly
homogeneous parcels, mixes two parcels at possibly different temperatures,
and adds/removes prescribed heat. Each bounded temperature trial recomputes
the original common-tangent phase equilibrium. Pure sensible enthalpy plus
the pair's original temperature-dependent excess closes the isobaric energy balance, including
the heat capacity contributed by changing coexistence compositions. A
homogeneous unstable feed can demix adiabatically and change temperature;
heating/cooling can remove/create the miscibility gap without changing either
component inventory. Positive heat is INTO the liquid inventory. This is an
enthalpy balance with constant-pressure boundary work, not a constant-volume
internal-energy closure.

`BinaryCaloricStateV1::isobaric_heat_capacity_j_per_k` differentiates that same
enthalpy at fixed component inventory. A homogeneous state sums the original
species Cp(T) and `n*C*x_A*x_B`; an equilibrated two-liquid state also includes
`n [H_interaction(T)/T]^2/R * (1-2x)^2 x(1-x) / [1-2 chi x(1-x)]` from changing coexistence
compositions. The retained positive log-odds root supplies a cancellation-free
response near criticality, including its finite split-side limit. No new
phase solve, temperature differencing, fitted total Cp or allocation occurs.
An explicitly homogeneous nonequilibrium parcel retains its fixed-composition
constraint. Temperature-independent interactions have zero redistribution
capacity, including at critical/binodal compositions. The redistribution
term also vanishes wherever the interaction enthalpy coefficient is zero.

A negative excess Cp is admitted only when total fixed-composition Cp stays
positive for every composition throughout the common source domain. Checking
`-C < (sqrt(Cp_A)+sqrt(Cp_B))^2` at both endpoints and all original source
knots proves this without temperature/composition sampling. A shared total-Cp
lower bound includes negative mixing capacity at the feed composition and
the nonnegative redistribution response; it supports conservative thermal
storage bounds. Cropped endpoints use outward-rounded source interpolation.

The strict heat-capacity query refuses distinct one-sided critical/binodal
limits and feeds unresolved by the retained numerical coexistence bracket.
`isobaric_heat_capacity_response` exposes the same smooth derivative, but
explicitly qualifies a selected homogeneous-side tangent at either kind of
boundary for an implicit thermal solve. That tangent changes no enthalpy or
accepted energy residual and is not reported as a unique physical Cp. Source
and arithmetic errors still refuse; there is no general fallback capacity.
Focused tests compare the response with independent original enthalpy and
entropy differences for upper/lower critical solution laws, variable pure Cp,
inventory/species scaling, the near-critical jump and boundary qualification.

The thermal residual integrates sensible increments directly from each
original parcel and differences the mixing products before inventory scaling.
It never subtracts large absolute reference enthalpies to recover a small
heat input. Compensated sums retain cancellation between sensible and excess
terms. Stored reference enthalpy is a rounded constitutive diagnostic; the
reported incremental heat residual is the acceptance account. Both explicit
temperature-width and absolute energy budgets must pass, unless a sampled
endpoint is an exact numerical energy root. A supplied bracket stays inside
all source domains; unsupported phase roots, out-of-bracket heat, exhausted
work, unrepresentable inventories and cancellation return no completed state.
Original parcels remain immutable and cannot be reused with another pair or
caloric law. Evaluation is allocation-free after constructing bounded source
curves, apart from shared-owner reference counting; no source querying occurs
inside the root loop.


Relative entropy uses the SAME source curves and accepted phase states.
`PureLiquidHeatCapacityV1::entropy_change_j_per_mol_k` analytically integrates
each original affine Cp segment divided by temperature. Positive endpoint
weights, local series and scaled sums retain small temperature increments;
underflow of a nonzero contribution refuses. Reversing the endpoints negates
the same ascending integration.

`BinaryLiquidCaloricModelV1::entropy_change_j_per_k` requires identical
component inventories and the original caloric owner. It combines direct
pure-component Cp/T increments with the actual homogeneous or equilibrium
mixing-entropy change, before extensive scaling. For fixed-composition
isobaric equilibrium endpoints this includes phase redistribution in dH/T.
Accepted source/phase/root budgets still limit precision. Prescribed heat
alone supplies no entropy account for an unspecified external reservoir.
Focused tests cover affine Cp, tiny increments, subnormal arithmetic, source
and inventory refusal, UCST/LCST Clausius bounds and adiabatic demixing.

`BinaryCaloricStateV1::reduced_entropy_j_per_k` supplies a source-relative
entropy coordinate for transport between compositions. It integrates each
original pure Cp/T from the shared sensible datum and retains the actual
accepted homogeneous or coexistence mixing entropy, including signed excess
Cp. Only common linear pure-species entropy constants are removed; those
cancel under conserved species. No absolute species entropy or pressure
dependence is inferred. Small fixed-inventory changes still use the direct
increment above. Native tests compare signed excess terms, pure-datum limits,
local increments and a two-liquid binodal against independent formulas.

G0/G1/G3 tests in `src/solution/caloric/tests.rs` check exact source integrals,
ideal unequal-capacity temperature mixing, both heat-of-mixing signs, an
independent high-precision two-phase thermal root, the additional equilibrium
heat capacity, reversible heating/cooling across coexistence, adiabatic
demixing, small heating against large reference energy, source mismatch,
domain/work refusal, cancellation and the phase-weighted Gibbs identity.
Fixtures validate numerical behavior, not measured properties of a named
chemical pair. The model still excludes vapor/solid phases, variable pressure,
spatial diffusion, nucleation rates, capillarity and hydrodynamic evolution.

## Finite calcium-sulfate phases and lime neutralization

`MineralReaction` provides gypsum precipitation with both hydrate waters,
anhydrite precipitation, anhydrite hydration, and the net hydrated-lime/sulfate
neutralization endpoint. Original exact atom/charge admission governs all
four. `CaSO4(s)` and `CaSO4.2H2O(s)` remain distinct participants; the latter
contains four hydrogens and six oxygens per formula unit.

The matched NBS anhydrite and gypsum-selenite packs retain eight literal
molar-mass, formation-enthalpy, independent formation-Gibbs and entropy facts
at 298.15 K and 100000 Pa. The insoluble anhydrite row is distinct from the
source's soluble alpha/beta rows and hemihydrates. The selenite variety and its
two waters of crystallization remain explicit. Uncertainty stays unstated;
these source values do not replace the existing USGS evaluations.

The `calcium_sulfate` example binds the actual seven compatible NBS packs to
the existing single-reaction `equilibrate_ideal` owner. A caller selects gypsum
or anhydrite and supplies free calcium, free sulfate, finite solid inventory,
explicit charge-balancing Na+/Cl- spectators and a fixed ideal activity basis.
Gypsum also tracks finite water and conserves water plus crystal-hydrate water.
The activity-reference amount remains independently fixed; the model does not
recompute solvent molality, solution volume or water activity from hydration.

The result distinguishes an interior saturation equilibrium from exhaustion
of the supplied solid, reports signed extent, remaining solid/source mass,
calcium/sulfur/charge/water balances and signed standard reaction heat.
A separate source-mass-driven lime example admits a prescribed finite extent
against explicit lime, acid and sulfate stocks. Source receipts, independent
direct/caloric Gibbs routes and finite-reactant refusal remain visible.

The model follows one selected free-ion reaction coordinate. It does not
jointly solve sulfate protonation, calcium-sulfate ion pairing, pH or competing
solid-phase selection, and supplies no nucleation/scale-growth kinetics,
electrolyte activity correction, operating heat duty or plant removal forecast.
Actual-source tests check the eight facts and scope refusals, all four reaction
energies and mass accounts, finite reagent limits, independently derived ideal
saturation/common-ion results, complete solid dissolution and basis scaling.

## Immutable admitted species definitions

`ConservedReaction::species_charge` and `species_elements` expose the exact
charge and nonzero elemental counts retained from the original admitted
matrices for an active participant. Absent and zero-coefficient species
return `None`; neither method parses a label or infers a formula.
The immutable accessors let coupled process owners verify specific transfer
models, including the neutral monovalent-counterion/gold-cyanide pair, while
retaining the original reaction conservation certificate.

## Sulfite reagents and the free-cyanide SO2/air endpoint

Seven `MineralReaction` equations retain crystalline sodium sulfite/metabisulfite
preparation, HSO3-/SO3(2-) dissociation, productive free-CN oxidation and separate
bisulfite/sulfite side oxidation. Four matched NBS packs supply 15 literal claims
at 298.15 K and 1 bar. Negative conventional ionic entropy stays conventional;
missing metabisulfite Gibbs energy/entropy cannot be replaced by its enthalpy.
The existing acid/base owner admits the explicitly limited HSO3-/SO3(2-) family.

`cyanide_treatment sulfite` resolves actual source packs and completes finite
preparation and oxidation stages before printing accepted output. Productive
and side reactions consume the same original dissolved oxygen and reagent
stocks. Sodium, sulfate, water and proton equivalents stay in the inventory;
nitrogen stays in cyanide/cyanate. Source mass, enthalpy and source identities
are reported. Exact conservation, source facts/refusals, independently summed
heats, staged demand and shared-pH checks are G0/G1/G3 evidence in
`xtask/tests/mining_sulfite_sources.rs`. No gas-transfer, copper-catalyst kinetics,
WAD-CN conversion, automatic pH, commercial dose or treatment-completion claim
is made. See `data/matdb/seed-v1/MINING_SULFITE_TREATMENT.md`.

## Finite peroxide treatment and the ammonia acid/base family

`MineralReaction` adds free-CN oxidation by aqueous peroxide, cyanate hydrolysis
to ammonium/carbonate, ammonium dissociation, and separate peroxide-decomposition
endpoints for gas and dissolved oxygen. All five equations pass the original
exact atom/charge admission. Cyanate's `OCN[-](aq)` identity stays distinct from
the fulminate isomer; aqueous peroxide/ammonia are not pure-liquid/gas aliases.
Four matched NBS packs supply sixteen literal M/H/G/S facts, checked against
original table images and the corrected 1989 CSV. Source-state selectors,
ionic entropy conventions and unstated uncertainty remain explicit.

The `cyanide_treatment` example admits oxidation and decomposition against one
finite peroxide stock, then hydrolysis against completed first-stage cyanate
and water. Both stages and all source diagnostics must succeed before accepted
output is printed. Original network rules reject reagent overdraw and same-stage
product borrowing. Independent direct-Gibbs and H−TS results remain visible.
Nitrogen stays in the named species; prescribed extents establish no rate,
WAD/total-cyanide removal, catalyst demand, water-quality result or operating
heat duty.

Its separate `speciate` mode binds source-backed ammonium, carbonate and
uncomplexed-cyanide families to the same aqueous pH owner. Caller analytical
totals, explicit spectators and unit-activity amount are required. NBS unit
molality uses solvent mass, not solution volume. This ideal fixed-solvent
result is not a finite-solvent H/O ledger and does not silently equilibrate
residual cyanate/peroxide, metal complexes or gas loss. Actual-pack integration
tests cover source rows, five reaction sums, finite-stage limits, mass closure
and coupled family/charge/mass-action balances.

## Coupled finite metal/ligand competition and aqueous pH

`acid_base::complexation::IdealLigandSystem` couples one original neutral
`HL -> H+ + L-` acid family to explicitly declared, disjoint aqueous metal
families. Each family has one to eight original exactly conserved association
steps `parent + q L- -> complex`, with q between one and eight. The first step
may consume two ligands and subsequent steps one each, as in the declared
Cu(I)/dicyano/tricyano/tetracyano family. Exact shared elemental/charge
definitions, common standard references and identical reused selected
formation-Gibbs sources are required. The original queries, answers, receipts,
reactions and aqueous model remain reachable from the result.

The caller supplies finite analytical ligand and metal totals, every other
acid-family total, explicit spectator-ion stocks and the common ideal activity
amount `n_ref`. No counterion or missing complex is invented. At trial free
ligand activity l, metal-member weights are `beta_j*l^n_j`, normalized in
logarithmic coordinates. A reduced copy of the original aqueous owner treats
free L- and charged metal members as spectators and solves electroneutrality;
neutral complexes retain their own positive amounts. It then infers
`n_HL=n_L*a_H/Ka`. A bounded logarithmic root closes the total of free,
protonated and metal-bound ligand. Other acid families respond to the same pH.

The retained inner pH bracket contributes its required-ligand variation to
root signs and convergence; an ambiguous sign cannot discard an endpoint.
The original full aqueous owner restores the complete HL/L family at the
accepted uncomplexed pool. Normally that pool is the independently positive
sum `n_L+n_HL`. Near floating logarithmic resolution, strong ion cancellation
can make inferred HL poorly conditioned; a well-resolved analytical component
difference `ligand_total-bound_ligand` may instead set the pool BEFORE that
full solve. This is an alternate component equation, not normalization or
clipping of previously solved species. An unresolved trace difference keeps
the independent positive-sum route. In both cases, independently reduced final
species inventories, original aqueous charge balance, all metal association
mass actions and the bracket containing the restored free ligand must meet
the requested controls. The reported ligand numerical spread is zero for an
analytically fixed final pool; final pH-bracket uncertainty still participates
in the association mass-action acceptance.

Results expose every metal member's original identity, charge, bound-ligand
count, amount, activity and log fraction; finite component totals, actual free
and bound ligand, the complete original aqueous report, independent component
and ligand residuals, log-equilibrium residual, numerical free-ligand bracket
and nested work counts. Exact empty pools stay empty; representable positive
trace stocks are retained even if their fraction underflows. Unrepresentable
positive stocks, overlapping or inconsistent sources, inadequate inner or
outer resolution, out-of-bracket roots, exhausted work and cancellation return
no partial result. Ligand and charge tolerances control those residuals; they
do not imply a universal pH error bound in strongly conditioned solutions.

This fixed-temperature ideal model supplies no omitted metal redox states,
precipitation, nonideal activity corrections, finite-solvent balance,
kinetic constants or stability claim against undeclared complexes. Seven
analytical test groups cover manufactured simultaneous states, two/one/one
ligand chains, competing metals, buffering, extensive scaling, conditioned
charge cancellation, empty and trace inventories, source identity and bounded
work/refusal. Actual matched NBS Cu(I)/CN/HCN/water packs separately exercise
source-backed competition and its coupling to the gold-leaching runtime.

## Native silver and zinc cementation of precious-metal cyanides

`MineralReaction` includes gas- and dissolved-oxygen native-silver cyanidation
and separate gold/silver zinc-cementation endpoints. Every equation is
admitted by the original exact atom/charge owner. The cyanide complexes retain
their full ligand inventory; Zn(CN)4(2-) is distinct from Zn(II), free cyanide
and solid zinc cyanide. Native silver does not represent sulfide ore or an
electrum activity model.

The original `ReactionNetwork` admits both cementation extents against one
initial zinc stock, with no clipping, inferred priority or same-stage product
borrowing. Actual NBS packs supply point-state standard enthalpy and independent
Gibbs diagnostics. Literal 1982-row values were independently checked against
the corrected 1989 CSV; both source references remain in the four profiles.
Ionic entropy remains conventional; source disagreements and unstated
uncertainty are not reconciled. No plant recovery, rate, side-reaction demand,
deaeration, filtration or operating heat duty follows from these endpoints.

The `precious_metal_recovery` example supplies explicit finite leach or mixed
cementation stages. `xtask/tests/mining_precious_metal_sources.rs` compiles the
actual packs, checks all 16 new facts and four H/G/S balances, shared-zinc
refusal, source-mass closure, stage separation and source applicability.

## Source-bound ideal gas/aqueous partition

`stoichiometry::thermodynamics::gas_liquid::IdealGasLiquidPartition` binds an
original exactly conserved neutral equation `species(g) -> species(aq)` with
coefficients -1/+1 to matched formation-Gibbs sources. Selected claims must
own the declared gas/pure-species or aqueous/infinite-dilution phase and the
same temperature, reference pressure and elemental convention. Original
queries, selected answers and receipts remain in the standard reaction owner.

The explicit ideal model uses `a_aq=n_aq/n_ref`, `a_g=p_partial/p_ref`, and
`a_aq=K*a_g`, with `ln(K)=-Delta_r G/(R*T)`. The caller supplies `n_ref` as an
amount at unit aqueous activity; unit molality is not a volume concentration.
For a finite ideal headspace, `p_partial=n_g*R*T/V_g` closes an algebraic
partition of a supplied gas-plus-aqueous total. Headspace volume and aqueous
activity basis stay separate. Results expose both stocks, partial pressure,
independently reduced component and log-equilibrium residuals, and the original
model. Empty totals have zero stocks and undefined log-equilibrium residual.

Logarithmic ratios avoid overflow and retain representable trace stocks even
when the fractional partition underflows. Neither phase is calculated by
subtracting the other from the total. A dominant stock can round to the total
while a positive trace remains separately represented; its balance discrepancy
is retained. Unrepresentable positive stocks or pressures refuse. This bounded
algebra has no iterative work or persistent runtime state; numerical residuals
are diagnostics, not interval certificates or propagated source uncertainty.

The source temperature remains fixed. An ideal partial-pressure evaluation is
not an experimentally validated pressure range. This owner supplies no measured
Henry coefficient, solvent-density conversion, salinity/nonideal correction,
dissolved chemical reaction, mass-transfer coefficient or transfer kinetics.
Six focused test groups cover independent finite-partition and activity laws,
extensive scaling, determinism, zero inventories, logarithmic trace stocks,
source/phase/reference refusals and numerical range loss.

## Finite-headspace solubility coupled to aqueous pH

`acid_base::gas_liquid::VolatileAcidBaseSystem` binds original
`IdealGasLiquidPartition` owners to neutral members of the original
`IdealAcidBaseSystem`. A neutral member may be the acid, the base or an
intermediate protonation state. Each selected family has one gas identity;
multiple volatile families share the same headspace and the same aqueous
proton activity. Admission checks exact shared atom/charge definitions and
identical original selected source answers/receipts, temperature, standard
pressure and elemental reference convention. A gas identity cannot also be
an aqueous spectator ion.

Inputs explicitly distinguish a volatile family's **gas plus all aqueous
forms** total from an unbound family's aqueous-only total. Headspace volume
and the common amount at unit aqueous activity remain independent; no volume
to molality conversion or counterion is inferred. At each trial pH, the
original aqueous member fractions give the neutral fraction `alpha`, and
`n_g/n_aqueous = alpha*V_g*p_ref/(n_ref*R*T*K)` closes the finite family.
Ionization therefore changes dissolved uptake and gas partial pressure while
the same original charge balance determines pH. The extra gas weight is
neutral, preserving the monotonic charge residual. The implementation reuses
the original single pH bisection, without nested equilibrium solves or
per-iteration allocations. Original nonvolatile-family semantics are retained.

The original complete aqueous result reports dissolved-only family totals.
Per volatile family, the enclosing result retains gas and total aqueous
amounts, partial pressure, the independently summed component balance, and
the recomputed logarithmic gas/neutral mass-action residual. Both phase tails
are represented independently, including a trace smaller than a subtraction
from the dominant total can resolve; unrepresentable positive stocks refuse.
Empty pools have zero gas and aqueous stocks with undefined gas mass action.
Source, bracket, numerical budget and cancellation refusals publish no result.

G0/G1/G3/G4 fixtures cover independent manufactured volatile-acid and
volatile-base endpoints, base-promoted acid uptake with depletion of the same
finite gas, multiple gases sharing pH and a nonvolatile buffer, extensive
scaling, deterministic declaration order, positive trace stocks, empty pools,
source/identity mismatch, invalid inputs, exhausted work and cancellation.
These are ideal-model arithmetic checks, not measured-solubility validation.
The common source temperature, ideal gas law and `a=n/n_ref` approximation
remain fixed. This path adds no salting-out or other nonideal activity law,
fixed-pressure bubble, thermal balance, transfer kinetics, finite solvent
inventory, redox, complexation or precipitation. The phase distinctions and
fixed-volume mass-action formulation follow the
[USGS PHREEQC gas-component equations](https://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/html/final-10.html).

## Coupled pH and finite mineral precipitation

`acid_base::precipitation::HydrolysisPrecipitation` composes the original
`IdealAcidBaseSystem` with one explicitly selected, exactly conserved reaction
`M(z+) + w H2O -> solid + z H+`. The cation and solid coefficients are -1/+1;
water consumption is a positive integer and z is the original charge number
between 1 and 32. The original reaction-network owner checks shared elemental
and charge definitions. Aqueous, pure-liquid and pure-solid source phases must
match their declared roles. All reused species must retain identical selected
source answers/receipts and a common reference state. The chosen solid does not
automatically establish stability against omitted competing phases.

The input is the finite analytical total of this dissolved cation plus this
solid, alongside all acid/base-family totals and fixed spectator ions. The
cation cannot also occur in a family or caller spectator list. Each trial
partition calls the original complete aqueous charge/family solve. Thus
precipitation changes proton activity through electroneutrality; pH is not
frozen at a pre-precipitation value. At coexistence, `a_M = a_H^z / K_precip`.
If the all-dissolved state is undersaturated, the solid remains absent and its
mass-action equality is not forced. An empty component has undefined saturation
and a separately reported `EmptyComponent` status.

The outer residual decreases monotonically with dissolved amount: increasing
cation raises pH under the inner charge balance and decreases proton activity.
Bounded bisection in log dissolved amount retains positive trace metal that a
reaction-extent subtraction could erase. The caller declares a positive search
floor, logarithmic partition/saturation tolerances and iteration bound; both
inner and outer solves poll cancellation. Unsupported brackets, inadequate
numerical resolution, inconsistent sources, cancellation and exhausted work
return no partial result. Present-solid success requires both the logarithmic
bracket and saturation residual to meet the requested tolerances.
The residual includes its remaining variation over the inner pH bracket;
an approximate pH cannot silently satisfy a tighter mineral tolerance.
Phase absence also requires its inequality sign to hold across that bracket.

Results retain the original phase/family models, total/dissolved/solid amounts,
full aqueous speciation, independently reduced component residual, natural-log
IAP/K for the dissolution direction, numerical bracket and nested work counts.
Solid amount is evaluated from the total minus the solved dissolved amount;
binary64 rounding can make a large solid amount equal the total while a tiny
positive dissolved amount remains separately represented. The exposed balance
residual preserves that rounding discrepancy; no stock is clipped or silently
discarded. Numerical brackets and residuals are estimates, not certificates or
experimental uncertainty bounds.

The model keeps the original common `a=n/n_ref` ideal activity approximation,
fixed solvent basis and unit water activity. It does not close finite solvent
H/O stocks, include omitted metal complexes/nonideal activities, choose between
solids, or predict nucleation, rate, reagent purity or commercial lime dose.
Analytical tests cover a manufactured charge/saturation state, a buffering
family at the same pH, absent/empty phases, trace dissolved stocks, scaling and
source/work/cancellation refusal. Actual USGS iron/goethite/H+/OH-/water packs
exercise acid suppression and declared-base promotion of precipitation. The
`mineral_precipitation` example exposes these ideal-solution calculations.

## Finite salt precipitation coupled to an analytical acid family

`acid_base::precipitation::salt::SaltPrecipitation` binds one exactly conserved
integer-stoichiometry precipitation equation p M + q A -> solid, optionally
consuming solvent water. M is an independent cation; A is the terminal anion of
an existing acid/base family. Positive original coefficients p and q satisfy
exact charge balance, and the source equation has one solid formula unit.
The original 1:1 construction and results remain supported.

Every family, water autoprotolysis and the selected phase retains its original
compatible formation-Gibbs sources. Cation and analytical family totals include
p and q component units per mole of selected solid. At each trial the same
aqueous owner solves all protonation states and charge balance. Saturation
uses the original ion activities to powers p and q. Precipitation cannot borrow
a missing component or retain a dissolved family total while consuming it.

The outer logarithmic variable is the limiting dissolved formula-equivalent
amount (cation/p or family/q), retaining positive traces independently of a
large solid stock. The excess component remains explicit. Compensated
stoichiometric comparisons preserve nearly tied inventories and avoid overflow
of a representable partition. The all-dissolved endpoint must satisfy the
absence inequality; a present phase must close saturation and both component
balances with their original coefficients. The anion coefficient also scales
the retained inner-pH uncertainty bound. Work limits, caller bounds, source
refusals and cancellation publish no partial partition.

The original aqueous owner prepares canonical identities and inventory buffers
once per phase solve. Nested salt trials change only the selected family and
cation amounts, reuse the original scalar charge equation, and retain compact
pH results without constructing complete species reports. One final source-
backed aqueous report is materialized after restoring the accepted amounts.
Public aqueous and finite-gas entrypoints reuse this same prepared path.
Focused native tests retain the original 1:1/gas/complex/hydrolysis coverage
and add 3:2 salts, both limiting inventories, common ions, finite exhaustion,
positive traces and difficult stoichiometric scaling.

`ScoroditePrecipitation` supplies the Fe/As/O/H-conserved dihydrate equation.
The `scorodite_partition` example and `mining_scorodite_sources` integration use
the NMK2014 arsenate/iron/scorodite evaluation with its CODATA-based water and
ionic auxiliaries. The nine profiles preserve nine Gibbs values and eight
reported uncertainty magnitudes; a reported uncertainty with unstated coverage
does not acquire a confidence level or covariance. These profiles remain
separate from the earlier NBS arsenic and iron evaluation.

This is one selected pure phase in an ideal, fixed-solvent solution. It does
not choose the stable solid assemblage, include ferric hydrolysis or dissolved
iron-arsenate complexes, model nonideal activities or adsorption, predict
nucleation/rates, or establish a mine-water treatment yield. The source is
limited to 298.15 K and 100000 Pa; no hot pressure-oxidation extrapolation is
inferred. Fe/As component and charge closure are distinct from finite-solvent
H/O conservation.

## Finite-rate ionic dissolution and crystal growth

`acid_base::precipitation::kinetics::MineralKineticsV1` binds one original
`SaltPrecipitation` or `HydrolysisPrecipitation` source to an explicitly supplied
`MineralSurfaceRateV1`. Actual dissolved component stocks and actual solid
formula units are declared separately. Preparation and basis changes equilibrate
only water ions and analytical acid/base families; they do not move mineral
mass or force a present crystal to saturation. The complete analytical stocks,
original spectator identities and source model remain attached to the state.

The selected fixed-surface constitutive law is `du/dt = K_A*(1-Omega)`, where
u is dissolved solid formula equivalents, `Omega = IAP/K_dissolution`, and
K_A has units mol formula units/s. For a salt, the original integer coefficients
supply `ln(Omega) = p*ln(a_M) + q*ln(a_A) + ln(K_precipitation)`; hydrolysis
uses `ln(a_M) - z*ln(a_H) + ln(K_precipitation)`. The original aqueous owner
solves charge balance and every acid-family fraction at each candidate
inventory. Thus pH, common-ion excess and protonation affect the rate through
the actual free-ion activities, with no invented neutral saturation fit.

`ExistingSolidOnly` blocks spontaneous crystal creation when the actual solid
stock is zero. `SeededSurface` explicitly declares a persistent growth surface
without supplying any component moles. Neither option infers a nucleation law
or an empirical kinetic coefficient from equilibrium formation energies.
This is the unit-affinity-exponent form of the empirical mineral rate family
described in the USGS PHREEQC theory, [Rate Expressions](https://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/html/final-20.html).
Its thermodynamic sign is consistent with the original reaction affinity;
material applicability and the surface coefficient remain caller declarations.

One bounded implicit scalar phase stage reuses the prepared original aqueous
solver. Full-step/two-half-step comparison controls the accepted time evolution;
only the accepted fine trajectory is published. Separate root, time and aqueous
work limits, original pH uncertainty and cancellation remain explicit. Reports
retain actual aqueous speciation, phase/component balances, saturation and work
counts. Solid exhaustion is a bounded endpoint, not a claimed exact continuous
event time. A failed or cancelled advance leaves its input state unchanged.

Each stage bounds its possible transfer using the initial rate and aqueous
pH spread before evaluating a trial inventory. A short step within the declared
pH domain does not require speciation of the distant, fully depleted donor.
A nonzero accepted transfer must change at least one retained phase stock;
unresolved increments refuse instead of silently accumulating lost mass.
A final time remainder below the minimum step may be omitted only when its
whole-motion estimate fits both phase-stock tolerances and the absolute extent
tolerance. Its normalized bound remains in the reported maximum accepted error.

Strong precipitation solves the implicit relation in the logarithm of the
remaining dissolved formula stock when an extent subtraction would lose that
stock. The residual is `ln(Omega(u))-ln(1+(u_start-u)/(h*K_A))`; both terms
are evaluated without overflowing the rate product or saturation ratio. The
original salt stoichiometry and terminal-anion charge balance, or the original
hydrolysis charge balance, supply its monotone lower slope. Aqueous pH spread
and floating-point logarithm margins contribute to numerical root-location
estimates, with both relative trace accuracy and the caller's extent/amount
tolerances enforced. These are numerical estimates, not interval certificates.

The positive remaining donor is retained directly alongside the solid stock.
A state close to saturation may remain unchanged when either its estimated
maximum possible later relaxation or the motion over the whole remaining
window fits the same tolerances. The time-limited allowance uses the remaining
window rather than accumulating independently omitted small stages. Root and
stationary allowances contribute to the existing temporal error comparison,
so a small positive dissolved trace remains usable on subsequent advances.

The `mineral_precipitation --kinetics` mode runs the same original USGS
iron/goethite source selection with explicit initial dissolved/solid stocks and
a separately supplied surface rate. Its equilibrium comparison is diagnostic;
it never replaces a finite-rate state.

This capability holds source temperature, pressure and ideal activity conventions
fixed, and omits competing minerals, additional dissolved complexes, nonideal
ionic activities, reaction heat, finite-solvent H/O balance and evolving particle
geometry.

## Coupled source-backed ideal aqueous acid/base speciation

`stoichiometry::thermodynamics::acid_base::IdealAcidBaseSystem` couples
water autoprotolysis and up to 32 analytical protonation families. Each family
contains one to eight original, exactly conserved one-proton dissociations.
Original selected formation-Gibbs sources supply Ka and Kw at a common T, p
and elemental reference. Reused species must retain identical selected
source answers/receipts; the existing reaction-network owner checks exact
shared elemental and charge definitions. Overlapping families and spectators
that duplicate reactive species refuse.

For a shared proton activity h, each family's member weights are
`w_j=product(K_1..K_j)/h^j`, normalized in logarithmic coordinates. Water
supplies `aOH=Kw/h`. One monotone charge-balance root includes all family
amounts and explicitly declared spectator-ion stocks. The solver never
sequentially equilibrates independent reactions, invents counterions, or
substitutes total analytical cyanide for free CN-. It retains the complete
original model, analytical-member amounts/fractions/activities, spectator
stocks, pH bracket and recomputed charge/mass-action residuals.

Every aqueous species uses the SAME explicit `a=n/n_ref` approximation.
For the NBS unit-molality references, n_ref corresponds to solvent mass
times 1 mol/kg; solution volume does not determine it. Water activity is one
and solvent basis is fixed: no closed finite-water H/O balance, nonideal
activity coefficients, metal complexation, redox, precipitation, volatilization
or leaching rate is claimed. The HCN+CN- pool excludes metal-bound and WAD
cyanide. Numerical brackets/residuals are estimates, not uncertainty or
interval-arithmetic certificates.

The caller bounds pH, charge accuracy and iteration work. Out-of-bracket
solutions, numerical range loss, inconsistent sources or cancellation return
no partial result. Eight independent analytical test groups cover mixed
mono/polyprotic families, water/strong-ion quadratics, counterions, dilution,
extensive scaling and tiny positive members. All 63 stoichiometry Cargo tests
pass. Three actual-source integration tests compile the five NBS packs and
check coupled chemistry independently, including the water contribution to
HCN/cyanide-salt roots. The Cargo-built `mineral_speciation --materials
PACK_DIR` example runs five declared-stock scenarios and prints all residuals.

## Source-backed arsenic protonation and peroxide oxidation

Five `MineralReaction` endpoints retain exact H/As/O and charge conservation:
the three successive As(V) proton losses, the first As(III) proton loss, and
`H3AsO3(aq) + H2O2(aq) -> H3AsO4(aq) + H2O(l)`. Six compatible NBS aqueous
profiles supply literal molar masses and formation H/G/S at 298.15 K and
100000 Pa. All claims require the source oxidation-state flag; the As(III)
claims additionally require their original component-plus-water convention.
Unstated uncertainty and unavailable heat capacities remain explicit.

`mineral_speciation --family arsenate --materials PACK_DIR` constructs one
three-stage family in the existing shared-charge owner, retaining all four
As(V) forms, water autoprotolysis and explicit spectator inventories. Original
Gibbs rows determine each Ka. No later evaluation is substituted to tune the
source values. `mineral_processing --reaction
arsenous-acid-peroxide-oxidation` exposes source-backed amounts and H/G/S
diagnostics; the conserved-inventory API bounds an explicit extent by the
actual peroxide stock. A favorable standard Gibbs value does not set an
extent, rate or peroxide selectivity.

`mining_arsenic_sources` compiles all ten participating source packs through
the ordinary compiler. It checks the 24 new claims and receipts, four direct
reaction H/G/pKa values, independently evaluated triprotic fractions over
pH 2-12 and extensive scaling, peroxide-limited inventory advancement, and
missing-context, extrapolation, missing-property and cancellation refusals.
The As(III) family contains only its first proton step. These sources supply
neither strongly alkaline As(III) speciation nor mineral dissolution,
nonideal electrolyte behavior, sorption or arsenic fixation.

## Source-backed aqueous cyanide, protons and dissolved oxygen

Three matched NBS/1989-errata packs add twelve literal molar-mass and H/G/S
claims for molecular HCN(aq), H+(aq) and O2(aq). Neutral HCN uses the original
`ao` entry, not the adjacent already-ionized `ai` component. Dissolved oxygen
retains its own nonzero formation properties and aqueous unit-activity
reference; gas oxygen is a distinct participant. Hydrogen retains the
individual-ion H+-zero convention. All claims remain at 298.15 K / 100000 Pa
under `source-elemental-reference-nbs1989`.

The original conserved mineral owner adds HCN dissociation, water dissociation
and an Elsner endpoint consuming dissolved oxygen. All three are runnable in
`mineral_processing`. The existing single-reaction ideal equilibrium owner
can now consume the complete aqueous HCN/water source sets with an explicit
amount at unit activity. It does not infer that amount from volume. The direct
Gibbs data give pKa=9.232640551245206 and pKw=13.995246497841048 on their
retained reference basis. A dissolved-O2 Elsner extent has direct
DeltaG=-407118 J, distinct from the gas-O2 endpoint's -390718 J.

Four actual-source tests compile the nine compatible packs and check all
literal rows, source conventions, the finite HCN quadratic, water equilibrium,
solvent-basis scaling and dissolved-oxygen limits. The three actual CLI cases
run with both thermodynamic routes. These source additions alone do not
provide coupled speciation, gas transfer, cyanide-complex competition or
dissolution kinetics. Complete sources and commands are documented in
`data/matdb/seed-v1/MINING_AQUEOUS_CYANIDE.md`.

### Source-complete gold cyanidation

The existing conserved Elsner endpoint now has a complete six-pack NBS source
set, exposed as `mineral_processing --reaction gold-cyanidation`. Gold,
cyanide, oxygen, liquid water, dicyanoaurate and hydroxide retain their own
molar mass, H/G/S claims, original source receipts and 298.15 K / 100000 Pa
domain under `source-elemental-reference-nbs1989`. Individual-ion entropy
remains conventional under the H+(aq)-zero reference. No aqueous cyanide
species is substituted for solid commercial NaCN or total cyanide.

The actual source compiler and original reaction/thermodynamic owners give
−583916 J/mol equation extent for H, −642.358 J/(mol extent K) for S and
−390718 J/mol extent from directly tabulated G. The independent H−TΔS route
differs by 1678.9623 J/mol; both are retained. NBS's historical Au(CN)2− value
is not silently combined with the later USBM evaluation. Four source-backed
integration tests cover all 24 claims, receipts, both thermodynamic routes,
activity response, cyanide/oxygen inventory limits and reference-state refusal.
The endpoint does not supply dissolution rate, pH/speciation, oxygen transfer,
preg-robbing, adsorption or actual gold recovery. The complete source guide
and runnable commands are in `data/matdb/seed-v1/MINING_GOLD_CYANIDATION.md`.



## Mixed-phase stoichiometric and standard-enthalpy balances

`stoichiometry::ConservedReaction` admits one equation through the existing
`fs-qty` exact elemental and charge matrices/certificate. It computes consumed
and formed molar amounts, converts a caller-specified reacted mass to extent,
and bounds extent from explicit reactant inventories. Canonical species order,
bounded layouts and exactly convertible coefficients make the calculation
deterministic; missing/duplicate species, nonfinite inputs, imbalance and
unrepresentable nonzero arithmetic refuse. It owns no evolving reactor state.

`limiting_extent` rounds a ratio down when an FMA residual reveals that its
coefficient product would exceed available stock. `advance_inventory` applies
an explicit forward extent to a complete inventory of active species, including
zero for absent products. It returns canonical finite nonnegative inventories
without mutating inputs; excessive extent, overflow and nonzero changes below
inventory resolution refuse rather than clipping or disappearing. The exact
certificate proves the integer equation; binary64 final amounts are separately
evaluated, with no claim of exact numerical mass conservation.

`stoichiometry::minerals::MineralReaction` provides explicit reference
endpoints: lime hydration, calcite/dolomite acid neutralization to dissolved
CO2 or retained bicarbonate, hydrated-lime neutralization, malachite acid
leaching to dissolved copper and CO2, separate calcite/malachite pathways to
the USGS source-defined un-ionized H2CO3(aq), the net ionic Elsner gold-cyanidation
balance, pyrite oxidation to either ferrous ions or ideal ferric hydroxide,
ferrous oxidation to ferric ions, ferric attack on pyrite,
magnesite calcination to crystalline MgO and gaseous CO2, magnetite
oxidation to hematite, cuprous-oxide oxidation, and carbon-monoxide reduction
of either cupric or cuprous oxide to crystalline copper and gaseous CO2,
dry pyrite roasting to hematite/SO2, SO2 oxidation to gaseous SO3, SO3 hydration
to pure liquid H2SO4, the net pyrite-to-hematite/pure-acid endpoint, and copper
electrowinning with an oxygen-evolving anode.
The bicarbonate pathways consume half as many protons per mole of carbonate
mineral as the CO2 pathways. The equations do not assign a universal measured
acid consumption to an ore, nor infer gas transfer or commercial reagent purity.
Pyrite and magnetite equations use integer multiples; extent is per equation
as written. Endpoint choice is declared, never inferred from pH, Eh, time or
thermal phase stability. The actual JANAF source packs close carbonate and
iron-oxide standard reaction heats at 298.15 K / 1 bar; they do not determine
kiln duty, reaction rate or decomposition temperature.
The copper equations separately retain oxide oxidation state, reductant and
oxygen demand, and metal production. Their JANAF formation sources give
−282838 J/mol for `2 Cu2O + O2 -> 4 CuO`, −126932 J/mol for
`CuO + CO -> Cu + CO2`, and −112288 J/mol for `Cu2O + CO -> 2 Cu + CO2`.
These are per written equation at the admitted standard state. The runnable
`mineral_processing` example accepts explicit compatible source sets,
prints each consumed/formed mass, retains independently rounded molar masses,
and exposes their mass-sum residual instead of adjusting a source value.
The sulfur endpoints reuse the existing pyrite JANAF source with its explicit
`source-phase-pyrite=1` condition. Their standard heats are −3339566 J/mol for
`4 FeS2 + 11 O2 -> 2 Fe2O3 + 8 SO2`, −197846 J/mol for
`2 SO2 + O2 -> 2 SO3`, and −132394 J/mol for
`SO3(g) + H2O(l) -> H2SO4(l)`. The separately balanced five-source net equation
`4 FeS2 + 15 O2 + 8 H2O -> 2 Fe2O3 + 8 H2SO4(l)` gives −5190102 J/mol,
equal to the sum of the three stages at their matching extents. It explicitly
assigns all sulfur to pure acid; it does not predict sulfur capture, catalytic
conversion, acid concentration, dilution heat, or plant heat duty. The example
accepts this five-pack source set as `--reaction pyrite-to-acid`.
`faraday_yield` converts integrated charge, positive electron count, supplied
molar mass and explicit current efficiency to product amount and mass.

`standard_enthalpy_change` consumes a source pack/query for every active species,
using `standard-molar-formation-enthalpy` in J/mol. It requires the selected
claims themselves to declare a common typed temperature and pressure, an
explicit elemental-reference flag, standard-state identity and phase. Caller
query flags alone cannot supply missing source authority. Selection and receipt
verification stay in `fs-matdb`; the result retains every original evidence
sample/uncertainty and source receipt. Products-minus-reactants enthalpy is
negative for exothermic constant-pressure conversion. No derived uncertainty
band is asserted. Multiplying by extent must not double-count chemical energy
already present in a total-energy formulation.

G0/G1/G3 tests cover exact balances, source-backed CaO hydration through the real
pack compiler, distinct pyrite and carbonate endpoints, acid-limited malachite
inventory advancement, Elsner reagent ratios, source-backed copper oxide
redox heat/metal/reductant balances with limiting inventories, sulfur-route
oxygen/water/acid amounts and a source-backed Hess cycle, Faraday copper yield,
unit/extent scaling and invalid/mislabelled reference states. This slice does
not extend the ideal-gas EOS/species model, select thermodynamic pathways,
authenticate caller species-to-pack associations, solve aqueous speciation,
predict reaction rates/recovery/current efficiency, or supply heat of dilution.
Independently rounded source molar masses need not sum bit-exactly even when
elemental conservation is exact; they are not a new mass-conservation proof.

### Ferric copper-sulfide endpoints

Four additional `MineralReaction` endpoints name chalcopyrite, chalcocite,
covellite and bornite explicitly. They consume respectively 4, 4, 2 and 12 mol
Fe3+ per mole of mineral, form 1, 2, 1 and 5 mol Cu2+, and assign all sulfur to
orthorhombic S on an atomic molar basis. They preserve exact elemental/charge
certification and the original finite-inventory transaction. These are chosen
net reaction endpoints, not automatic mineral-pathway or passivation selection.

Five new crystalline USGS B2131 packs pair with the existing conventional
aqueous Cu2+, Fe3+ and Fe2+ packs at 298.15 K / 100000 Pa. The
`mineral_processing` example accepts the four `ferric-*-leach` cases, with
explicit mineral/Fe3+/Cu2+/Fe2+/S pack order and the existing H, G, entropy and
activity consumers. All participants retain their original source receipts.
The source tests check literal molar quantities, phase/reference refusal,
conserved reagent-limited extent and sequential ferric regeneration. Bornite's
independently tabulated G and H-TΔS differ by 1.98631 kJ/mol mineral; both remain
visible. No source reconciliation, kinetic rate, extrapolated phase stability,
heat of dilution or plant recovery is inferred. Details and runnable commands
are in `data/matdb/seed-v1/MINERAL_THERMOCHEMISTRY.md`.

### Mixed-phase formation Gibbs energy and explicit activities

`stoichiometry::thermodynamics` evaluates an already conserved equation from
complete `standard-molar-formation-gibbs-energy` claims in J/mol. It reuses the
enthalpy source selection, typed reference-state checks and receipt verifier;
all active species must share temperature, reference pressure and elemental
reference. Aqueous claims must themselves fix the hypothetical unit-activity,
infinite-dilution state. Charged participants additionally require the H+(aq)
zero convention, identified from the exact charge matrix rather than a parsed
species label. Non-aqueous sources require neutral pure-species references.
Formation enthalpy minus absolute species entropy is never substituted for a
missing formation-Gibbs claim.

`standard_gibbs_change` returns the canonical source sum Delta_r G°, finite
ln(K°) = -Delta_r G°/(R T), and direct K° only when positive and representable.
The original sources, their uncertainty and their receipts remain attached.
`at_activities` requires one finite positive dimensionless activity for every
active species, including solvent and solid phases, relative to those source
states. Activities greater than one are valid inputs; absent activities and
zero-activity logarithmic limits refuse. Fixed-order compensated logarithmic
sums evaluate ln(Q) without forming an overflowing quotient. The result
Delta_r G = R T ln(Q/K°) is a local dG/dxi diagnostic; forward affinity is its
negative. It is not a finite-extent total free-energy change, an equilibrium
composition solve, a phase-selection certificate or a kinetic prediction.
No concentrations, activity coefficients, pH or ideal-gas EOS are inferred.

The USGS B2131 source set covers pyrite oxidation, ferrous-oxidant regeneration,
ferric pyrite leaching and calcite bicarbonate neutralization. Two ferric-leach
extents plus seven ferrous-oxidation extents reproduce one integer net pyrite
oxidation equation, including cancellation of intermediate iron and acid used
by regeneration. Exact stoichiometry and source-backed H/G Hess sums are checked
independently. Each finite stage still needs its reactants in inventory.
The example exposes standard diagnostics through `--gibbs` and accepts repeated
`--activity 'SPECIES=VALUE'` inputs for a complete explicit activity state.
Pyrite oxidation retains finite ln(K°) even though direct K° overflows; the
calcite source regression verifies reversal of driving force on lowering proton
activity while holding product activities fixed.

G0/G3 regressions cover independent Decimal known answers, reaction reversal
and scaling, conserved reference shifts, source/activity permutation,
log-only range limits, retained receipts and missing/mismatched reference
conditions. The formulas follow IUPAC S05915 (standard equilibrium constant)
and A00115 (relative activity); numerical checks do not calibrate a plant.

### Finite-inventory equilibrium along one reaction

`ConservedReaction::equilibrate_ideal` solves the signed reaction extent from
complete finite stocks, direct formation-Gibbs sources and explicit activity
models. Each aqueous participant requires a positive caller-declared amount
at unit activity; its ideal law is `a = n / n_ref` on a fixed solvent basis.
Each pure solid/liquid requires an explicit `PureCondensed` law with unit
activity. The selected original source claim must own the matching phase;
caller labels cannot assign an aqueous ion to the condensed branch. Gases
and solid solutions require other laws. Standard-state source data do not
authenticate the chosen activity approximations or their mixture consistency.

The exact conserved equation supplies `n_i = n_i0 + nu_i * xi`. Both directions
are bounded by available stocks, using inward-rounded ratios and FMA residuals.
The solver does not clip inventories to zero, and a representable inventory
limit can retain a tiny positive remainder. Interior `ln(Q/K)` is monotone
because its derivative is the sum of `nu_i^2 / n_i` over ideal-aqueous species.
Bounded deterministic bisection requires BOTH an acceptable finite logarithmic
residual and an extent bracket within the caller's tolerances. Separate logs
avoid forming overflowing or underflowing activity ratios. A limiting aqueous
stock uses its analytic logarithmic endpoint; an unresolvable interior root
refuses rather than being called complete conversion. Budget exhaustion,
unrepresentable amounts and nonzero changes below stock resolution also refuse.

Pure-phase constraints can terminate at either inventory boundary with the
appropriate one-sided Gibbs inequality. A missing mineral may precipitate
when this chosen reaction favors it, or remain absent when its formation is
unfavorable. A zero-width feasible interval returns `InventoryBlocked` without
evaluating undefined activities. All-condensed reactions with exactly zero
standard Gibbs change return `NonUnique`, retaining the full feasible interval.
The result retains canonical initial/final inventories, signed extent, numerical
bracket, status, work count, declared activity laws and every standard source
receipt. The bracket is not an interval-arithmetic or experimental uncertainty
certificate. The finite result is equilibrium only along this reaction and
under these assumptions; it does not solve competing reactions, charge-balance
speciation, concentrated electrolytes, solvent depletion effects, stable-phase
selection, nucleation barriers or kinetic recovery.

Analytical tests cover ideal isomerization, equation reversal/scaling, mineral
saturation, finite dissolution, precipitation from an initially absent phase,
all-condensed degeneracy, zero stocks, source/model mismatch and bounded-work
refusal. Actual USGS goethite/iron/water/proton packs exercise source-backed
precipitation, dissolution and retained source receipts. The
`mineral_equilibrium` example requires explicit inventories and an aqueous
amount scale; it prints the assumptions and constrained result.

### Source-bound reversible electrochemical potential

`stoichiometry::electrolysis::ElectrolysisReaction::standard_potential` evaluates
formation-Gibbs sources on its own original conserved reaction, then computes
the signed typed voltage `E° = -Delta_r G°/(n F)`. The electron count is the
original productive-species association times its written-equation coefficient;
counts that cannot convert exactly to binary64 refuse. The result retains the
original electrolysis model, standard Gibbs evidence and receipts. Its
`at_activities` method reuses the same complete positive dimensionless activity
admission and reports local `E = -Delta_r G/(n F)` without changing the source
or electron convention.

For the USGS copper/water/oxygen/proton set, the four-electron written copper
electrowinning equation gives `Delta_r G° = 344000 J/mol` and
`E° = -0.891327190438547 V`. With proton activity 0.01 and all other activities
one, the potential is `-0.7730084910689824 V`. The sign follows the written
direction; the opposing magnitude is a reversible thermodynamic requirement,
not an operating cell voltage. Current efficiency, resistance, overpotential,
side reactions and finite-extent work are not inferred. The original
charge-driven inventory/yield API remains separate. See
`data/matdb/seed-v1/MINING_ELECTROWINNING.md` for the source set and runnable case.

### Shared finite reagent inventories and source-bound heat

`stoichiometry::network::ReactionNetwork` joins admitted equations while
retaining each active species' exact nonzero element counts and charge. A
shared identifier with conflicting chemical definitions refuses. Reactions
and the union of active species have canonical order; inventory and extent
inputs must cover their complete axes, including explicit zeros.

`advance_inventory` checks every reaction's gross withdrawals against one
initial inventory before crediting products. Calcite gangue and malachite
therefore cannot independently spend the same protons. Products are available
to a later caller-declared stage; they never fund simultaneous withdrawals or
an otherwise empty reaction cycle. Results retain per-reaction extents and
amounts, complete gross consumed/formed accounts and final species inventories.
The private floating-point reduction retains product residuals when checking
shared stock, so an overdraw that rounds back onto stock still refuses.
Inventory, extent and numerical-range failures return no partial result.

`advance_with_standard_enthalpy` also requires exactly one formation source
for every species in the network. A shared species uses the same original
source in every reaction. All selected sources must declare one common typed
temperature, pressure and elemental reference, including disconnected or
zero-extent reactions. The result contains each admitted extent, each original
`ReactionEnthalpy` with its evidence/receipts, its extent-scaled energy and the
compensated sum in joules. This is standard chemical enthalpy, not heat of
mixing/dilution, sensible heat, aqueous activities or a reactor temperature.
No material source or exact equation supplies a reaction rate or acid-allocation
policy; a caller must declare those separately.

## Prepared fixed-pressure reaction calorics

`IdealGasKineticsV1::reaction_enthalpies_at` evaluates canonical reaction
Delta h from the same coefficient-combined original NASA sources as Delta u
and Delta g. All original species, including spectators, retain source and Cv
admission; the successful prepared query allocates nothing.
`KineticWorkspaceV1::reaction_enthalpies` and `species_enthalpies` retain the
corresponding values from the actual temperature evaluation, with reusable
storage included in the reported workspace size.

These values support fixed-pressure enthalpy evolution in fs-gas. They do not
add reaction heat to a total-energy balance or replace the existing
fixed-volume Delta u derivatives. The original source offsets, continuity
policies and elemental reference cancellation remain authoritative. Focused
G3 caloric tests compare the typed standalone query, prepared arrays and
stoichiometrically conserved reference shifts.

`reaction_caloric_change_at` also returns the directly reduced standard-state
reaction entropy. `reaction_entropies_at` fills its canonical reaction array
without successful-call allocation. The original entropy coefficients and
explicit corrections are combined before evaluation: conserved elemental
datum shifts cancel, but unequal reaction entropy differences remain. Entropy
is never recovered by subtracting rounded Delta g from Delta h. These standard
state values exclude the evolving parcel's ideal mixing and expansion entropy.


> Status: ACTIVE, CODE-FIRST SLICES. This contract covers typed,
> provenance-bound NASA-9 ideal-gas standard-state evaluation and bounded
> frozen-composition ideal-gas mixture evaluation, plus a positive-state
> mechanical ideal-gas EOS rung and a bounded Wilke dilute-gas viscosity
> mixing closure, plus bounded standard-state ideal-gas reaction equilibrium
> and a caller-declared single-reaction stoichiometric mass-action closure,
> plus bounded ideal-gas TP/HP/TV/UV composition equilibrium with explicit
> elemental and electrical-neutrality constraints.
> Central batch compilation and test execution for the newest slice are
> pending; the parent thermochemistry bead remains in progress.

## Purpose and layer

`fs-thermochem` is the L1 thermochemical law-data and standard-state evaluation
layer. It reuses the exact species, element, reaction, stoichiometric, charge,
and conservation artifacts owned by `fs-qty`; it does not create a competing
chemistry identity or conservation system.

The first slice admits immutable NASA-9 coefficient cards from `fs-matdb`,
evaluates typed molar `cp`, `h`, and `s`, and derives `u` and `g` only under an
explicit gas/ideal-gas/reference-pressure/elemental-reference convention. The
second combines up to 128 strictly positive, canonically ordered components
into frozen ideal-gas molar and mass-specific properties. The third supplies
the typed mechanical closure from positive finite temperature, pressure, and
molar mass to density, molar volume, specific gas constant, and `Z = 1`. The
fourth combines up to 128 positive caller-supplied pure-species viscosities by
Wilke's dilute-gas rule at one explicitly declared common temperature. The
fifth combines one exactly conserved `fs-qty` stoichiometric column with all
and only its active NASA-9 models to evaluate standard reaction Gibbs energy,
`ln(K_p)`, and `K_p` when directly representable. The sixth evaluates
forward, reverse, and optional net progress for that one reaction only after
the caller explicitly declares stoichiometric coefficients as kinetic orders
and supplies positive dimensionless ideal-gas activities.

Direct runtime dependencies are L1 or lower: `fs-qty`, `fs-matdb`, and
`fs-math`. The direct `fs-evidence` edge is development-only so conformance
fixtures construct the same validity type that reaches runtime transitively
inside `fs-matdb` cards. The crate has no L2-L6 dependency. Its bounded
composition-equilibrium solver is an immutable local thermodynamic closure;
it owns no transport solver/protocol, evolving state, session, persistence,
or orchestration behavior.

### Trial-state caloric evaluation

Both NASA-9 species and frozen-mixture models expose an
allocation-free-on-success evaluate_properties path for nonlinear trial
states. It shares the full evaluator's operation tree, canonical reduction
order, domain checks, formation reference, mixing entropy and pressure
correction; it is not a constant-capacity approximation or a second law.
The immutable model retains its original source cards. Trial properties
do not contain a new evaluation receipt; accepted/publication states use
the full evaluate path. G5 tests require bit-identical properties across
NASA region boundaries, compositions and pressures, and identical physical
refusals. Existing receipt fields and evaluator versions are unchanged.

### Relative entropy

NASA-9 species and frozen ideal-gas mixtures expose direct entropy changes
between admitted endpoints. Factored polynomial differences and a stable
logarithmic ratio cancel common entropy references before rounding; distinct
region reference offsets are retained. Frozen mixtures include the pressure
term and cancel their unchanged ideal mixing term. Successful calls allocate
nothing. These are relative state functions, not reaction or mixing entropy,
and endpoint admission does not establish continuous source-path coverage.

### Relative enthalpy

NASA-9 species and frozen ideal-gas mixtures also expose direct molar and
mass-specific enthalpy changes. Factored differences of the original NASA-9
energy powers cancel common formation references before rounding, preserving
small sensible-energy changes for nozzle acceleration. Distinct region
coefficients and reference jumps remain observable. Ideal-gas enthalpy is
pressure-independent; these APIs keep composition fixed and do not replace
reaction, mixing, or open-volume chemical-energy balances. Successful calls
allocate nothing and reversing endpoints negates the result exactly. G0/G1
tests cover every caloric term, reference shifts, region jumps, variable heat
capacity, mixture mass basis and physical/numerical refusal boundaries.

### Reaction caloric differences

Kinetic evaluator version 5 and reaction-equilibrium evaluator version 2 share
one reduction of the original NASA-9 coefficients under each exact
stoichiometric column before evaluating reaction enthalpy, internal and Gibbs
energies. Conserved elemental enthalpy and entropy references cancel before
absolute species properties are rounded; input differences already lost when
source coefficients were authored cannot be recovered. In that reduction,
distinct formation energies, molecular-count pressure-work differences, raw
regional jumps and explicitly declared continuity offsets remain. One shared
temperature basis serves all reactions without additional per-reaction logs.
The same reaction energy drives reverse-rate temperature derivatives,
adiabatic temperature rates and the gas reactor's caloric Jacobian. Successful
evaluations reuse prepared storage. `reaction_internal_energies_at` exposes
the original reaction energies for conservative evolution; they are not an
additional heat source in a total-chemical-energy balance. G0/G3 regressions
cover elemental-reference shifts, changed binding energies, molecular-count
changes and raw versus explicitly reconciled source seams.

`reaction_caloric_change_at` exposes the same source-coefficient reduction for
one canonical reaction as typed `ReactionCaloricChangeV1` enthalpy, internal
energy and Gibbs-energy changes. It evaluates Delta h and Delta u independently,
retaining Delta h - Delta u = Delta nu R T without subtracting absolute species
energies. Full original-source support and positive-Cv admission include
spectators. Successful calls allocate no workspace; cancellation and an invalid
reaction index refuse. Compartment burning uses this owner for enthalpy power
conversion and the original internal-energy reference correction.

### Thermodynamic mass-action entropy

The declared stoichiometric mass-action evaluator uses A/(R T) =
log(P_forward) - log(P_reverse) + log(K_p), with compensated fixed-order
summation. Net progress uses the dominant directional rate and deterministic
expm1, retaining near-equilibrium differences without subtracting rounded
gross rates. The original finite-range directional statuses remain authoritative;
unrepresentable nonzero net rates remain unavailable. dimensionless_affinity
exposes that same force; entropy_production_w_m3_k returns R*r*A/(R T).
One-sided exact species absence retains the positive-infinite ideal-mixture
limit. Both absent directions and exact kinetic shutdown give zero; numerical
range loss gives None. These are instantaneous local constitutive quantities,
not a finite-step entropy or global passivity certificate. Common kinetic
multipliers preserve affinity. G0/G1/G3 regressions cover detailed balance,
tiny departures in both directions over 400 decades of rate scale, multiplier
semantics, exact absence and numerical underflow.

## Public types and semantics

- `SpeciesId`, `ElementId`, `ReactionId`, `ElementalMatrix`,
  `StoichiometricMatrix`, `ChargeVector`, `ConservationCertificate`,
  `MassAmountBasis`, and `verify_conservation` are direct `fs-qty` re-exports.
  Exact `A N = 0` and `z^T N = 0` authority remains in that one owner.
- `MolarHeatCapacityV1`, `MolarEnthalpyV1`, `MolarEntropyV1`,
  `MolarInternalEnergyV1`, and `MolarGibbsEnergyV1` are distinct semantic
  wrappers around coherent-SI quantities. Equal dimensions do not make two
  thermodynamic meanings interchangeable.
- `ElementalReferenceIdV1` is a validated opaque convention name. It does not
  authenticate a reference-element table or establish formation-property
  authority.
- `StandardStateConventionV1` retains the phase, reference EOS, positive finite
  reference pressure, and elemental-reference id required by every model.
  Version 1 supports exactly `Gas` plus `IdealGas`.
- `Nasa9RegionV1` is admitted only from a complete, immutable
  `ConstitutiveModelCard` with law id `nasa9-standard-state`, law version 1,
  state-schema version 0, zero internal state, exactly one positive finite `T`
  validity axis, and exactly `a0..a8` plus `reference_pressure`.
- `Nasa9StandardStateModelV1` binds one `SpeciesId`, positive finite
  `MolarMass`, explicit standard-state convention, and 1 through 16 ordered
  regions. Regions may have explicit gaps but may not overlap. At a shared
  boundary, the upper region wins deterministically. Because the current
  `fs-matdb` law-card schema carries no species, molar mass, phase, EOS, or
  elemental-reference identity, those fields are caller declarations: they
  are receipt-bound but not source-authenticated. Only reference pressure is
  checked directly against every retained card.
- `Nasa9EvaluationV1` contains typed properties and an immutable
  `Nasa9EvaluationReceiptV1`. The receipt binds evaluator and deterministic
  math versions, gas-constant bits, species and molar-mass bits, input
  temperature bits, every convention field, selected region and bound bits,
  all coefficient bits, and the full source-card content identity.
- `Composition`, `CompositionBasis`, and `SemanticError` are direct `fs-qty`
  re-exports. `FrozenIdealGasMixtureModelV1` accepts only mass or mole fractions,
  requires an exact canonical sum of one, refuses zero entries and duplicates,
  converts mass fractions through the shared typed molar-mass path, and stores
  both declared and evaluated mole compositions.
- `FrozenIdealGasMixturePropertiesV1` carries molar `cp`, `cv`, `h`, `u`, `s`,
  and `g`, mixture molar mass, mass-specific `cp` and `cv`, and `R/M_mix` as
  distinct semantic wrappers. `FrozenIdealGasMixtureReceiptV1` binds both
  composition bases, exact sums and `sum(x ln x)`, common conventions,
  evaluator/math/quantity versions, `T`, `p`, `p0`, mixture molar mass, and
  every nested species standard-state receipt.
- `Nasa9StandardStateModelV1::evaluate_heat_capacity` and
  `FrozenIdealGasMixtureModelV1::evaluate_heat_capacities` evaluate only the
  original source Cp polynomial and frozen Cp/Cv mass/molar conversions. They
  share the complete evaluator's operation tree, region selection and canonical
  reduction, including at reconciled source seams, with no allocation or
  logarithms. Their narrower admission checks temperature/source support and
  finite requested capacities; it does not evaluate H/S/U/G, pressure, or
  their reference arithmetic, and does not prove positive Cv over a path.
  Complete accepted-state receipts and solver caloric-stability admission
  remain separate. `tests/heat_capacity.rs` checks bitwise agreement, variable
  capacities, both composition bases, source limits and arithmetic refusals.
- `IdealGasEosV1` retains one positive finite caller-declared molar mass.
  `evaluate_pt` returns semantic wrappers for positive finite mass density,
  molar volume, the same canonical `MassSpecificGasConstantV1` used by frozen
  mixtures, and compressibility factor. Its `IdealGasEosReceiptV1` binds the
  evaluator/quantity versions, dedicated
  `MechanicalEquationOfStateV1::IdealGas` rung, gas-constant bits, all input
  bits, and all output bits. The receipt does not authenticate species, phase,
  molar mass, or material identity.
- `WilkeGasViscosityComponentV1` binds one canonical `SpeciesId`, positive
  finite typed molar mass, and positive finite typed dynamic viscosity. The
  component is an immutable caller declaration; it does not authenticate a
  pure-species transport law or the state at which that law was evaluated.
- `WilkeGasViscosityMixtureV1` admits at most 128 active components, moves
  fractions with components into canonical species order, requires exact-one
  active composition, converts mass fractions to mole fractions through typed
  molar masses, and refuses volume fractions. Its one retained positive finite
  temperature is the caller's common-state declaration.
- `WilkeGasViscosityEvaluationV1` contains one typed mixed dynamic viscosity
  and an immutable `WilkeGasViscosityReceiptV1`. The receipt binds the
  operation-tree/math/quantity versions, explicit `Wilke` rule, temperature,
  both composition bases and exact sums, canonical species, molar masses,
  supplied viscosities, fractions, every outer-component denominator, and the
  final result.
- `IdealGasReactionEquilibriumV1` retains one bounded
  `StoichiometricMatrix`, an existing matching `ConservationCertificate`, one
  selected `ReactionId`, every nonzero exact coefficient, and all and only the
  corresponding NASA-9 models in canonical species order. It requires at
  least one reactant and product, coefficients exactly representable in
  binary64, and one exact phase/EOS/reference-pressure/elemental-reference
  convention across all active species.
- `IdealGasReactionEquilibriumEvaluationV1` carries typed standard reaction
  Gibbs energy, finite dimensionless `ln(K_p)`, a positive finite direct `K_p`
  when exponentiation stays representable, and an exact-field receipt. A
  finite log result survives direct-constant overflow/underflow under an
  explicit `EquilibriumConstantStatusV1`; zero or infinity is never published
  as a physical equilibrium constant.
- `IdealGasReactionEquilibriumReceiptV1` binds the evaluator/math/quantity
  versions, exact gas constant, reaction and canonical column, stoichiometric
  and conservation-certificate identities, temperature, common convention,
  every coefficient/Gibbs/contribution bit pattern, every nested NASA-9
  receipt, standard reaction Gibbs energy, log constant, and direct-constant
  representation state.
- `ThermodynamicReverseRateClosureV1` binds one positive finite forward
  reaction-progress-rate scale to that equilibrium authority. Forward and
  reverse scales have units mol m^-3 s^-1 because the assumed mass-action
  activities are dimensionless `p_i / p0` ratios. Its evaluation retains
  `ln(k_reverse / k_forward)`, an optional direct ratio and reverse scale, an
  explicit ratio/scale range status, and a receipt containing the complete
  nested equilibrium receipt. It does not choose kinetic orders or evaluate a
  net reaction rate.
- `DeclaredStoichiometricMassActionV1` binds that reverse-rate closure to the
  explicit
  `MassActionKineticOrderConventionV1::CallerDeclaredStoichiometricCoefficients`
  choice. Its input activities are nonnegative finite dimensionless `p_i / p0`
  ratios covering all and only the active reaction species. Its evaluation
  retains canonical forward/reverse log activity products, optional direct
  activity products and directional progress rates, an optional signed net
  progress rate, per-direction range statuses, and the complete nested
  thermodynamic receipt. The type records a caller's law declaration; it does
  not infer or certify elementary-reaction kinetics.
  Operation-tree version 3 admits exact absence: an absent participant gives
  its direction a `ZeroActivity` status, exact zero product/rate and a retained
  negative-infinite log. This is distinct from underflow of a positive product,
  which still supplies no direct rate. Products need no artificial trace seed.

## NASA-9 operation tree

For absolute temperature `T` in kelvin and coefficients `a0..a8`, version 1
uses this fixed scalar operation tree:

```text
cp/R  = a0/T^2 + a1/T + a2 + a3*T + a4*T^2 + a5*T^3 + a6*T^4
h/RT  = -a0/T^2 + a1*ln(T)/T + a2 + a3*T/2 + a4*T^2/3
        + a5*T^3/4 + a6*T^4/5 + a7/T
s/R   = -a0/(2*T^2) - a1/T + a2*ln(T) + a3*T + a4*T^2/2
        + a5*T^3/3 + a6*T^4/4 + a8
```

`R` is pinned to the binary64 encoding of
`8.31446261815324 J mol^-1 K^-1`, the decimal product fixed by the post-2019 SI
Boltzmann and Avogadro constants. Its exact evaluator bits are retained in the
receipt. `fs_math::det::ln` is used; platform libm is not part of the evaluator.
Under the only admitted v1 convention, `u = h - R T` and `g = h - T s`.

The dimensioned logarithm is interpreted as `ln(T / 1 K)`; the evaluator's
coherent-SI `Temperature` scalar is therefore the numerical kelvin value.

NASA/TP-2002-211556 labels its seven heat-capacity coefficients `a1..a7` and
its two integration constants `b1,b2`. This crate follows Cantera's zero-based
nine-slot layout:

```text
local/Cantera [a0, a1, a2, a3, a4, a5, a6, a7, a8]
NASA report   [a1, a2, a3, a4, a5, a6, a7, b1, b2]
```

The NASA report's tabulation convention uses
`R = 8.314510 J mol^-1 K^-1`; this evaluator deliberately uses the current-SI
value above, matching current Cantera policy. Direct physical-unit values from
the report consequently differ by about 5.7 parts per million even when the
dimensionless polynomials match. External oracle tests must compare `cp/R`,
`h/(R T)`, and `s/R`, or explicitly normalize both sides to one declared gas
constant. A raw SI comparison across those dialects is not admissible evidence.

The formula source is NASA/TP-2002-211556, *NASA Glenn Coefficients for
Calculating Thermodynamic Properties of Individual Species*:
<https://ntrs.nasa.gov/api/citations/20020085330/downloads/20020085330.pdf>.
Cantera's independent NASA-9 documentation is retained as a development
cross-reference, never a runtime dependency:
<https://cantera.org/dev/reference/thermo/species-thermo.html>.

`Nasa9StandardStateModelV1::evaluate_heat_capacity_temperature_derivative`
analytically differentiates that same retained Cp polynomial and returns a
dimensioned J/(mol K²) quantity. It allocates nothing and uses the original
temperature/region and finite-Cp admission. At a shared source seam, both Cp
and its derivative must agree between the two original regions to rounding
tolerance; raw jumps or derivative corners refuse. Source-domain endpoints
use inward derivatives. Caloric reference offsets do not alter Cp or its
tangent. This does not extrapolate, certify positive Cp/Cv along a path, or
turn source uncertainty into a derivative error bound. Focused regressions
cover original polynomial channels, independent differences, reference shifts,
source gaps, endpoints and shared-seam refusal.

## Frozen ideal-gas mixture operation tree

For exact canonical mole fractions `x_i`, common reference pressure `p0`, and
species standard states evaluated at `T`, version 1 uses canonical species
order and defines:

```text
M_mix = sum(x_i M_i)
cp    = sum(x_i cp_i)
cv    = cp - R
h     = sum(x_i h_i)
Qx    = sum(x_i ln(x_i))
Lp    = ln(p) - ln(p0)
s     = sum(x_i s_i) - R Qx - R Lp
u     = h - R T
g     = h - T s
cp_m  = cp / M_mix
cv_m  = cv / M_mix
R_mix = R / M_mix
```

`ln(p) - ln(p0)` avoids overflow or underflow from first forming `p/p0`.
Every listed fraction is strictly positive, so `ln(x_i)` is total; absent
species must be omitted rather than represented by zero. These are
frozen-composition derivatives only. Reacting or equilibrium derivatives do
not inherit `cv = cp - R` without their own composition-response terms.
Cantera's ideal-gas implementation is a development cross-reference for the
mixing and pressure terms, not a runtime dependency:
<https://www.cantera.org/3.0/doxygen/html/d7/dd4/IdealGasPhase_8cpp_source.html>.

## Local response from density derivatives

`eos::density::DensityThermodynamicResponseV1` derives signed isobaric
expansion, isothermal/adiabatic compressibility, Cv and sound speed from one
homogeneous fixed-composition phase's original `rho`, `rho_T|p`, `rho_p|T`,
absolute temperature and isobaric Cp. It uses the local Gibbs identities,
including `Cp-Cv = T beta²/rho_p`, without an ideal-gas or incompressible
assumption. All original scalars and the directly evaluated Cp-Cv correction
remain available. Fixed-size evaluation allocates nothing. Nonpositive
mechanical/caloric response, a caloric residual within 16 rounding units of Cp,
and unrepresentable arithmetic refuse. Zero and negative expansion are valid.

G0/G1/G3 regressions compare ideal gases, an independently differentiated Gibbs
potential, original IF97 liquid/steam response, density scaling, a density
maximum and small slopes. This local algebra grants no source authenticity,
phase selection, uncertainty bound or global integrability of independently
supplied density and Cp functions. Source ownership remains with its consumer.

## Mechanical ideal-gas EOS operation tree

For positive finite molar mass `M`, absolute temperature `T`, and pressure `p`,
version 1 uses ordinary binary64 arithmetic in this fixed order:

```text
R_specific = R / M
RT         = R * T
V_m        = RT / p
rho        = M / V_m
Z          = 1
```

Every exposed derived scalar must be strictly positive and finite. A zero or
non-finite intermediate/result refuses without a partial evaluation. This is a
mechanical ideal-gas model-domain check, not a physical phase-stability test.

## Wilke dilute-gas viscosity operation tree

For positive mole fractions `x_i`, pure-species dynamic viscosities `mu_i`,
and molar masses `M_i`, version 1 evaluates canonical `(i, j)` pairs in this
fixed binary64 tree:

```text
r_mu       = mu_i / mu_j
r_M_rev    = M_j / M_i
q_M        = sqrt(sqrt(r_M_rev))
root_prod  = sqrt(r_mu) * q_M
n_base     = 1 + root_prod
numerator  = n_base * n_base
r_M_fwd    = M_i / M_j
d_rad      = 8 * (1 + r_M_fwd)
Phi_ij     = numerator / sqrt(d_rad)
D_i        = sum_j(x_j * Phi_ij)
mu_mix     = sum_i(x_i * (mu_i / D_i))
```

This is the usual Wilke expression in an explicitly pinned, algebraically
equivalent operation order. `fs_math::det::sqrt` supplies every square root;
the evaluator neither calls platform libm nor retains an `N x N` interaction
matrix. Every intermediate and fixed-order partial sum must remain strictly
positive and finite. Cantera's `GasTransport` documentation is retained as a
development cross-reference for the rule and units, never as a runtime
dependency:
<https://cantera.org/3.1/cxx/d8/d58/classCantera_1_1GasTransport.html>.

`WilkeGasViscosityEvaluationV1::temperature_tangent` consumes pure dmu/dT in
the evaluation's canonical species order and differentiates the same retained
values, fractions, masses and interaction operation tree. It returns dmu_mix/dT
and every dD_i/dT for reuse by the WMS conductivity owner. Signed/zero slopes
are admitted, the pure-species limit is exact, and invalid lengths or
unrepresentable arithmetic refuse. The original value evaluation and receipt
are unchanged. Work remains bounded O(N²), with stack logarithmic derivatives
and one O(N) result vector, without a pair matrix or perturbed value solve.
The caller must bind the supplied slopes to the original pure source laws;
this algebra establishes neither their physical validity nor derivative
uncertainty. Regressions compare complete mixture and denominator tangents
with differences of original evaluations, common scaling and signed pure limits.

## Standard-state reaction-equilibrium operation tree

For one exactly conserved reaction column with negative reactant and positive
product coefficients `nu_i`, common ideal-gas standard-state models, and
absolute temperature `T`, version 2 admits active species once in canonical
`SpeciesId` order and retains their complete original evaluation receipts.
Reaction calorics use the same compensated source-coefficient reduction as
kinetics, followed by the original NASA-9 enthalpy and entropy basis:

```text
delta_a_j = compensated_sum_i(binary64_exact(nu_i) * a_i,j)
h_offset  = compensated_sum_i(nu_i * declared_h_offset_i)
s_offset  = compensated_sum_i(nu_i * declared_s_offset_i)
delta_h   = R * compensated_sum_j(delta_a_j * h_basis_j(T)) + h_offset
delta_s   = R * compensated_sum_j(delta_a_j * s_basis_j(T)) + s_offset
delta_g   = delta_h - T * delta_s
R_T       = R * T
ln(K_p)   = -delta_g / R_T
K_p       = deterministic_exp(ln(K_p)) when positive finite
```

Inactive polynomial powers are skipped. Conserved elemental reference terms
cancel before species-property rounding, while physical formation differences,
raw source-region jumps and explicitly declared constant caloric corrections
remain active. The original per-species `nu_i * g_i^0(T)` values remain in receipts
as source diagnostics; their rounded sum is not the version 2 reaction total.
Every source evaluation and diagnostic contribution must still be finite.
Nested reverse-rate and mass-action receipts carry this version 2 equilibrium
receipt while preserving their original source NASA-9 receipts.

`K_p` is dimensionless because each ideal-gas activity is `p_i / p0`; all
active models must bind the same exact `p0`. Coefficients outside
`[-2^53, 2^53]` refuse because conversion would lose their exact integer
identity. Finite `ln(K_p)` remains the successful representation when direct
exponentiation overflows or underflows. This operation tree evaluates a
standard-state constant only; it does not solve for equilibrium composition.

For one positive finite forward progress-rate scale `k_forward`, version 1
then closes the reverse scale under the same dimensionless ideal-gas activity
convention:

```text
log_ratio = ln(k_reverse / k_forward) = -ln(K_p)
ratio     = deterministic_exp(log_ratio) when positive finite
k_reverse = k_forward * ratio when positive finite
```

The finite log ratio remains authoritative when either exponentiation or the
final scale multiplication leaves binary64 range. No zero or infinite direct
ratio/rate is published. This is only the thermodynamic coefficient relation
for a mass-action formulation using the retained standard-state activities;
it does not itself select kinetic orders.

## Declared stoichiometric mass-action operation tree

After the caller explicitly selects stoichiometric coefficients as kinetic
orders and supplies all nonnegative dimensionless activities `a_i = p_i / p0`,
operation-tree version 3 preserves the positive-state operation order:

```text
ell_i       = deterministic_ln(a_i)
ell_f       = fixed_order_sum_(nu_i < 0)((-nu_i) * ell_i)
ell_r       = fixed_order_sum_(nu_i > 0)(nu_i * ell_i)
A_forward   = deterministic_exp(ell_f) when positive finite
A_reverse   = deterministic_exp(ell_r) when positive finite
q_forward   = k_forward * A_forward when positive finite
q_reverse   = k_reverse * A_reverse when positive finite
q_net       = q_forward - q_reverse when both directions are direct
```

An exact zero activity sets its direction's log product to negative infinity
and publishes a `ZeroActivity` product/rate of exactly zero. No logarithm of
zero enters the finite positive-state reduction. Finite `ell_f` and `ell_r`
remain authoritative if a positive direct product or directional multiplication
leaves binary64 range. Every missing direct value has a per-direction status;
numerical underflow is never published as physical absence. `q_net` is absent unless both directional
rates are representable. This bounded tree is a declared algebraic closure,
not evidence that stoichiometric coefficients are physically valid kinetic
orders for the named reaction.

## Invariants

- ONE CHEMISTRY AUTHORITY: exact bookkeeping and conservation come from
  `fs-qty`; this crate only re-exports them.
- DIMS AT ADMISSION: `a0..a8` dimensions exactly cancel their documented
  temperature powers. Reference pressure has `Pressure::DIMS`. Foreign,
  missing, or dimensionally wrong parameters refuse.
- PROVENANCE IS RETAINED: the complete source card remains available and its
  canonical `fs-matdb` content identity is copied into every evaluation
  receipt. A receipt records provenance; it does not upgrade evidence color or
  certify source truth.
- NO IMPLICIT EXTRAPOLATION: evaluation outside all admitted regions, including
  an explicit inter-region gap, refuses.
- EXACT PRESSURE CONVENTION: every region's reference-pressure IEEE bits must
  equal the model convention's bits. Numerically close is not identical.
- FIXED REGION PRECEDENCE: exactly shared boundaries select the upper region;
  all other admitted endpoints are inclusive.
- TOTAL SUCCESS VALUES: no successful property is NaN or infinite. The first
  non-finite derived property causes a typed refusal and no partial result.
- DERIVED POTENTIALS ARE CONDITIONAL: `u` and `g` are exposed only by the
  version whose type-level alternatives admit gas plus ideal gas. Adding a
  phase or EOS requires a new explicit derivation and versioned semantics.
- BOUNDED METADATA: a model contains at most 16 regions, making selection and
  receipt construction bounded independently of caller input size.
- CANONICAL MIXTURE ORDER: component/fraction pairs sort together by
  `SpeciesId`; duplicates refuse. Caller permutation therefore cannot change
  reduction order, properties, or receipt.
- EXACT ACTIVE COMPOSITION: all declared and converted mole fractions must be
  strictly positive and sum to bit-exact `1.0` in canonical order. The wider
  `fs-qty` composition tolerance is not treated as exact thermodynamic
  normalization at this boundary.
- COMMON MIXTURE CONVENTION: every component must match phase, EOS,
  reference-pressure bits, and elemental-reference id. This checks internal
  consistency, not source authenticity or phase truth.
- MIXTURE RECEIPTS ARE PROVENANCE, NOT CERTIFICATES: nested card identities and
  exact arithmetic inputs permit replay but do not establish coefficient
  accuracy, stability, or evidence color.
- POSITIVE IDEAL EOS DOMAIN: `M`, `T`, and `p`, plus every exposed mechanical
  EOS output, are strictly positive and finite. `Z` is exactly one by model
  definition; this does not prove a real fluid is stable or ideal.
- BOUNDED QUADRATIC TRANSPORT WORK: Wilke count and composition-length gates
  run before component scanning, sorting, or `O(N^2)` pair work. `N <= 128`,
  canonical sorting uses no auxiliary interaction matrix, and every evaluator
  vector reservation is fallible before publication.
- ACTIVE WILKE COMPOSITION: zero fractions refuse instead of retaining inert
  receipt entries. Canonical declared and converted mole fractions must each
  sum to bit-exact one. Species duplicates refuse after canonical sorting.
- WILKE VALUES ARE CONDITIONAL: positive finite component values and a common
  caller-declared temperature are structural prerequisites, not proof that the
  values share a state or that a gas is dilute, single-phase, or described by
  Wilke's empirical approximation.
- ONE CONSERVATION PROOF OWNER: reaction equilibrium accepts only an existing
  `fs-qty` certificate whose exact stoichiometric identity matches. It neither
  reconstructs nor weakens `A N = 0` or `z^T N = 0` authority.
- BOUNDED REACTION WORK: matrix axes are each capped at 128, the admitted
  matrix identity at 16,384 cells, and one selected reaction at 128 active
  terms. All retained vectors reserve fallibly before publication.
- CANONICAL REACTION ORDER: supplied NASA-9 model order is irrelevant; active
  terms and fixed-order sums follow the canonical stoichiometric species axis.
  Missing, duplicate, inactive, or convention-incompatible models refuse.
- LOG AUTHORITY SURVIVES RANGE LOSS: every successful `ln(K_p)` is finite.
  Direct `K_p` is optional and appears only when deterministic exponentiation
  is strictly positive and finite; range loss is explicit, not clamped.
- REVERSE-RATE CONVENTION IS NOT INFERRED: reverse-rate closure consumes the
  exact equilibrium convention and admits only a positive finite typed forward
  progress-rate scale. It uses `k_reverse / k_forward = 1 / K_p` only for
  dimensionless standard-state activities; it does not silently crosswalk a
  concentration-based or dimensionful kinetic law.
- REVERSE RANGE LOSS IS EXPLICIT: the finite log ratio and nested equilibrium
  receipt survive ratio or final-scale overflow/underflow. A direct reverse
  scale appears only when both operations produce positive finite values.
- MASS-ACTION ORDERS ARE EXPLICIT: constructing the mass-action closure binds
  the caller-declared stoichiometric-order convention in the model and every
  receipt. Exact conservation bookkeeping is never treated as kinetic-law
  evidence.
- EXACT ACTIVITY COVERAGE: nonnegative finite dimensionless activities cover all
  and only the active reaction terms. Input order is irrelevant; duplicates,
  omissions, and inactive species refuse before thermodynamic evaluation.
- MASS-ACTION RANGE LOSS IS EXPLICIT: finite canonical forward and reverse
  positive-state log products survive direct-product and rate multiplication range loss.
  Per-direction statuses distinguish product overflow/underflow, unavailable
  reverse scale, and final rate overflow/underflow. Net progress appears only
  when both directional rates are direct and finite.
  Exact zero is a separate admitted physical state with a negative-infinite
  log sentinel; it is never used to hide positive-state underflow. The existing
  equilibrium owner still evaluates every bound reaction species at the actual
  temperature, including species with zero activity.

## Error model

Library paths are total and contain no `unwrap`, `expect`, or intentional
panic. `ThermochemErrorV1` reports source-card validation failures; law,
version, and state mismatches; missing, foreign, or wrongly dimensioned
parameters; invalid validity and convention data; invalid molar mass; empty,
excess, overlapping, or pressure-inconsistent regions; invalid or out-of-range
temperatures; and non-finite arithmetic. Float-bearing refusals preserve exact
IEEE-754 bits where those bits identify the rejected value.

`FrozenIdealGasMixtureErrorV1` adds count/length/zero/exact-sum/duplicate/basis
refusals, mass-to-mole underflow, typed expected/found cross-component
convention mismatch, invalid pressure or derived mixture molar mass, contextual
species-evaluation failure, and non-finite mixture fields.

`IdealGasEosErrorV1` separately names invalid molar mass, temperature, pressure,
and the first zero, negative, or non-finite EOS output/intermediate in fixed
operation-tree order. Float-bearing refusals retain exact IEEE-754 bits.

`WilkeGasViscosityErrorV1` separately names count and length gates, shared
composition failures, invalid component scalars and temperature, duplicate or
zero active species, non-exact declared or converted composition, mass-to-mole
underflow, bounded allocation stage, and the first failed pair/outer component
plus operation-tree intermediate. Float-bearing refusals retain exact IEEE-754
bits and no partial viscosity or receipt escapes.

`ReactionEquilibriumErrorV1` separately names matrix/cell/model bounds,
certificate drift, unknown/zero/one-sided reactions, unavailable matrix
entries, inexact integer-to-binary64 coefficients, duplicate/missing/inactive
models, convention drift, bounded allocation stage, nested species evaluation,
and the first non-finite term/sum/`R T`/log value. No partial reaction result or
receipt escapes.

`ReverseRateClosureErrorV1` names an invalid forward progress-rate scale, the
exact nested equilibrium refusal, or an otherwise impossible NaN reverse
scale. Positive finite ratio/final-scale overflow or underflow is a successful
typed log-only status rather than an error or clamp.

`StoichiometricMassActionErrorV1` names invalid activity values, count and
canonical species-coverage drift, duplicate activity identities, bounded
receipt allocation, the exact nested reverse-rate refusal, or an otherwise
impossible non-finite fixed-order intermediate. Expected activity-product or
directional-rate range loss is a successful typed status, not a clamp.

## Determinism class

Version 1 is fixed-order deterministic for identical inputs under the same
compiled arithmetic target. Region traversal is caller-significant order,
ties have one rule, parameter maps are canonical `BTreeMap`s, and logarithms
use the versioned deterministic `fs-math` implementation. Repeated evaluations
on one target are expected to be bit-identical and are covered by a G5 replay
test.

Mixture component/fraction pairs are canonicalized by species before every
fixed-order sum. The 128-component cap and nested 16-region cap make this
operation tree bounded. Declared-input permutation is expected to produce the
same model, properties, and exact-field receipt.

The mechanical EOS is one fixed allocation-free scalar tree. Repeated
same-target evaluation is expected to return bit-identical outputs and receipt.

The Wilke evaluator canonicalizes by `SpeciesId`, then executes fixed nested
`(i, j)` loops and fixed-order sums. Caller permutation cannot change the model,
result, denominator trace, or receipt. `fs_math::det::sqrt` and the explicit
evaluator version bind the elementary-math dialect; repeated evaluation on one
target is expected to be bit-identical.

Reaction equilibrium canonicalizes models to the exact stoichiometric species
axis, evaluates each active species once, and uses one fixed compensated
coefficient reduction shared with kinetics, then `fs_math::det::exp`.
Caller model permutation cannot change the admitted model,
standard reaction Gibbs energy, log/direct constant, or receipt. Forward and
exactly reversed columns retain opposite Gibbs/log values under the same
models.

Reverse-rate closure negates the retained equilibrium log bits, calls the same
deterministic exponential once, and performs one fixed multiplication. Replay
is exact on one target; exact reaction reversal negates the log ratio and makes
representable direct ratios reciprocal to the documented elementary-math
tolerance.

Declared mass action canonicalizes activities by `SpeciesId`, calls
`fs_math::det::ln` once per active term, accumulates separate directional log
products in fixed species order, and uses the fixed direct-product and scale
multiplication tree above. Caller permutation therefore cannot change the
evaluation or receipt; repeated same-target evaluation is expected to be
bit-identical.

Cross-ISA bit identity is not claimed until the central Gauntlet runs retain
evidence for both reference ISA families. NASA/mixture receipts record their
evaluator and `fs-math` versions; the allocation-free EOS receipt records its
evaluator and `fs-qty` versions because that tree calls no elementary math.

## Cancellation behavior

Model construction after region admission scans at most 16 regions, selection
scans at most 16 regions, and evaluation is one fixed scalar expression. There
is no useful tile boundary at which to poll a `Cx` in those paths.

`Nasa9RegionV1::from_card` is different: the generic input card may contain an
unbounded parameter/source collection or provenance strings, and upstream card
validation plus content hashing are input-linear and currently non-cancellable.
This slice makes no bounded-admission claim for hostile cards. A follow-up must
add explicit card byte/count budgets before this boundary is exposed to
untrusted bulk ingestion.

Frozen-mixture construction sorts at most 128 components and evaluation makes
one canonical pass whose nested species work is bounded by 16 regions each;
there is no useful `Cx` tile boundary. Future composition-equilibrium,
mechanism-wide kinetics, or database operations require explicit work budgets
and cancellation/drain semantics before landing.

Mechanical ideal-gas construction and evaluation are fixed-size scalar work
with no useful cancellation tile. They neither consume a `Cx` budget nor own
request-drain-finalize behavior.

Wilke admission performs bounded `N <= 128` sorting and basis conversion;
evaluation performs at most 16,384 pair interactions and fixed-size receipt
work. This first L1 closure therefore has a hard defensive work ceiling and no
useful asynchronous `Cx` tile. It does not claim caller-budget consumption,
deadline preemption inside scalar arithmetic, or request-drain-finalize
ownership. Larger or iterative transport closures require explicit budgets and
cancellation semantics.

Standard-state reaction admission hashes at most 16,384 already-constructed
stoichiometric cells and binds at most 128 active models; evaluation performs
at most 128 bounded NASA-9 evaluations and one canonical reduction. This first
closure therefore has a hard work ceiling and no useful asynchronous `Cx` tile.
Composition-equilibrium iterations, mechanism-wide kinetics, or larger sparse
reaction systems require explicit budgets and request-drain-finalize semantics.

Reverse-rate closure adds one deterministic exponential and multiplication to
that same bounded equilibrium evaluation. It introduces no mechanism loop,
activity product, evolving state, callback, or additional cancellation tile.

Declared mass action adds at most 128 deterministic logarithms, two bounded
fixed-order reductions, two exponentials, and three scalar rate operations.
It introduces no mechanism loop, evolving state, callback, or useful
asynchronous cancellation tile.

## Unsafe boundary

The crate uses `#![forbid(unsafe_code)]`. There is no unsafe boundary.

## Feature flags

None. No frontier or moonshot capability is promoted by this slice.

## Conformance tests

Inline tests in `src/lib.rs` define the current executable contract:

- G0 independently pins the hard-coded coefficient-dimension table, then basis
  fixtures exercise every NASA-9 coefficient channel against the published
  operation tree, plus a constant-`cp/R` closed form and both derived
  potentials.
- G0 exact chemistry verifies the shared `fs-qty` elemental and charge
  certificate for `2 H2 + O2 -> 2 H2O`.
- G3 falsifiers cover wrong law/version/state, missing/foreign/wrong-dimension
  parameters, wrong validity axis, invalid convention ids and pressure, invalid
  molar mass, empty/excess/overlapping/pressure-inconsistent regions, explicit
  region gaps and both adjacent endpoints, invalid/out-of-range temperatures,
  and finite inputs whose evaluation overflows.
- G5 replay pins bit-identical repeated evaluation and upper-region selection
  at a shared boundary. Receipt assertions bind coefficient and source-card
  identities.

Inline tests in `src/mixture.rs` add:

- G0 pure-species reduction; binary weighted `cp/h`, ideal mixing and pressure
  entropy, derived `cv/u/g`, molar mass, direct mass-specific `cp/cv/R_mix`
  oracles, and molar/mass-basis gas-constant identities; plus equivalent
  mass-to-mole basis conversion.
- G3 permutation invariance, pressure scaling, alternative Gibbs formulation,
  and refusals for count/length/zero/nonexact/duplicate/volume compositions,
  convention drift, invalid pressure, component-domain failure, positive mass
  fractions that underflow during mole conversion, and a derived molar mass
  that underflows to zero.
- G5 repeated evaluation with exact nested receipts, plus exact assertions and
  controlled mutations for version, convention, state, composition, and nested
  component receipt fields.

Inline tests in `src/eos.rs` add:

- G0 typed `p V_m = R T`, `p = rho R_specific T`, `rho V_m = M`, and exact
  `Z = 1` identities;
- G3 pressure/temperature rescaling plus zero, negative, NaN, infinity,
  intermediate/output overflow, and output underflow refusals;
- G5 bit-identical replay, exact assertions for every receipt field, and
  controlled state/output mutation coverage.

Inline tests in `src/transport.rs` add:

- G0 exact pure-species and equal-property limits, a binary independent Wilke
  expression, typed mass-to-mole basis conversion, and the exact 128-component
  admission boundary;
- G3 canonical permutation invariance, common viscosity scaling, count/length,
  invalid scalar/temperature, duplicate/zero/nonexact/volume composition, and
  mass-to-mole underflow and extreme-ratio arithmetic refusals;
- G5 bit-identical replay, exact assertions over every receipt field, and a
  supplied-viscosity mutation that changes the receipt.

Inline tests in `src/equilibrium.rs` add:

- G0 exact conserved-water reaction binding, canonical active terms, an
  independent analytic constant-Cp reaction Gibbs/log result, evaluator version
  2, and complete nested stoichiometric/conservation receipt identities;
- G0/G3 elemental enthalpy and entropy reference invariance, changed molecular
  binding energy, exact reaction reversal and reverse-rate source receipts;
  kinetic/equilibrium crosswalks also cover raw region jumps and explicit
  caloric continuity corrections;
- G3 unknown, zero, and one-sided columns; certificate drift; missing species
  model; reference-pressure mismatch; maximum exact coefficient with both
  log-only overflow and underflow retention; and first inexact coefficient
  refusal;
- G5 model-order replay plus forward/reverse standard Gibbs and log-constant
  polarity and reciprocal direct constants.
- G0 exact reverse/forward log identity, finite direct closure, typed units,
  nested-receipt equality, and direct-scale reconstruction;
- G3 invalid forward scales, nested temperature refusal, ratio overflow and
  underflow, and final reverse-scale overflow and underflow with retained log
  authority;
- G5 exact replay plus reaction-reversal log polarity and reciprocal direct
  ratios.
- G0 canonical declared stoichiometric activity exponents, forward/reverse
  products and typed directional/net progress, plus complete nested receipt
  binding;
- G0 exact absence gives zero directional progress, permits reverse reaction
  and remains distinct from numerical underflow of a positive activity product;
- G3 invalid activities, count/duplicate/missing/inactive coverage drift,
  nested thermodynamic refusal, log-only reverse-scale unavailability,
  activity-product overflow/underflow, and final progress-rate
  overflow/underflow with retained log authority;
- G5 caller-permutation invariance and bit-identical evaluation/receipt replay.

The newest tests are code-first and batch-verification pending. A sourced
external NASA/Cantera numerical oracle battery, an independently tabulated
Wilke battery, adversarial source-card mutation battery, and retained cross-ISA
evidence are required follow-ups; the synthetic algebraic fixtures do not
substitute for them.

## No-claim boundaries

This slice does **not** claim:

- authenticity, licensing, or accuracy of caller-supplied NASA coefficients;
- authenticity of the caller-declared species/molar-mass association to a
  source card whose schema does not itself carry those fields;
- authenticity of the caller-declared gas/EOS/elemental-reference association;
  in particular, a condensed-species card mislabeled as gas is not detected by
  this schema and must not be treated as authority for `u = h - R T`;
- an external numerical-oracle match or any evidence-color promotion;
- reacting/equilibrium composition derivatives beyond the dilute ideal-gas
  thermodynamic response described below,
  activities or chemical potentials beyond the retained ideal-gas `p_i/p0`
  standard-state convention, fugacity, departure functions, cubic/tabular
  real-gas or multiphase EOS behavior;
- a physical phase-stability verdict from the positive ideal-gas model-domain
  check, flash calculations, nonideal or multiphase equilibrium-composition
  solves, source-authenticated reaction rates, mechanism evolution or kinetics
  integration, thermal conductivity, diffusion coefficients, or transport
  solves;
- chemical meaningfulness, reversibility, kinetic accessibility, detailed
  balance, or reverse-rate authority from exact bookkeeping and a computed
  standard-state `K_p` alone;
- a validated kinetic order, elementary-reaction classification,
  concentration-unit conversion, forward-law authenticity, mechanism
  consistency, stiffness claim, or time integration from the caller-declared
  single-reaction mass-action result;
- authenticity or same-state consistency of caller-supplied pure-species
  viscosities, viscosity-law validity outside their source domains, dense-gas
  corrections, or physical applicability/accuracy of the Wilke approximation;
- uncertainty propagation, interval enclosure, rounding-error certification,
  or validity outside the exact admitted temperature regions;
- cross-ISA bit stability before retained G5 evidence exists;
- any L3 gas-state protocol, time evolution, L6 session state, ledger storage,
  planner authority, or admission policy.

Those capabilities require separate typed laws, contracts, tests, budgets, and
proof-bearing beads. These initial evaluators must not be used as evidence that
the broader `fs-thermochem` roadmap is complete.




## Explicit common kinetic multiplier

Mass-action operation-tree version 3 adds
`DeclaredStoichiometricMassActionV1::evaluate_with_kinetic_multiplier`.
The caller supplies one finite nonnegative dimensionless factor that applies
to both directional rates after the original thermodynamic evaluation.
The original `evaluate` path is the unit-factor case and preserves its
arithmetic; each receipt now retains the exact multiplier bits.

A physical zero factor produces `ZeroKineticMultiplier` and zero progress
without replacing real activities, activity products, term receipts or the
thermodynamic reverse closure. A direction already stopped by `ZeroActivity`
retains that more specific status. A positive-factor product that underflows
or overflows has no direct rate and retains an explicit range status. Original
log-only base rates remain unavailable; this extension does not rescale them
back into range. Invalid coefficients or source domains still refuse even
when the multiplier is zero.

The thermochemistry owner does not infer collider concentrations, efficiencies
or falloff models. The material gas coupling supplies those laws using its
existing gas inventory and source resolver. Focused native regressions check
unit-factor identity, common forward/reverse scaling, real zero with retained
activities, range loss and invalid-factor/source refusal.

### Explicit continuous caloric references

`Nasa9StandardStateModelV1::with_caloric_continuity` is an explicit model adjustment for contiguous NASA-9 regions. Raw construction never applies it automatically. Original source cards, card identities, coefficients, and Cp curves remain unchanged. The caller declares an anchor temperature and maximum absolute enthalpy and entropy corrections in molar units. The selected anchor region retains its original h/s reference throughout that region. For formation-reference preservation, choose a region containing the source's formation/reference temperature. Choosing a different region can change low-temperature formation values within the declared correction bounds.

Construction walks outward from the anchor and computes seam H/S differences by reducing the original source coefficients before evaluating their temperature basis. Common formation and entropy datums therefore cannot hide a correction, change its budget admission, or create a reaction heat jump. Every cumulative per-region correction must fit its declared bound. Source gaps, nonfinite corrections, exceeded bounds, and unrepresentable corrected properties refuse. Continuity applies to the coefficient-defined laws within floating-point rounding; independently rounded absolute endpoints are not forced to be bit-identical. These bounds describe the permitted reference adjustment; they are not estimates of source-data uncertainty. Cp may retain a source-fit jump at a region boundary. Continuity admission does not establish positive Cp or Cv.

Corrected evaluation derives u = h − RT and g = h − Ts from the same h/s owner. Corrected relative enthalpy and entropy integrate the retained Cp and Cp/T curves over each crossed region, avoiding subtraction of large absolute references or independently rounded seam corrections. Absolute and relative evaluations agree within floating-point rounding. Raw models retain their original region references, source jumps, and evaluator version 1. Explicitly corrected models use evaluator version 3 and bind the anchor, bounds, original card identities, and exact offset bits in immutable policy metadata shared by the model and accepted receipts. Construction and cross-region changes are bounded by the existing 16-region admission limit; successful property-only evaluation adds no allocation. G0/G3 regressions cover positive/negative common datums, all anchor regions, cumulative correction budgets, and reaction enthalpy/entropy/Gibbs consistency across the source seam.


## Ideal-gas composition, ionization and adiabatic equilibrium

`equilibrium::composition::IdealGasEquilibriumV1` binds the original typed
`ElementalMatrix`, signed `ChargeVector` and canonically ordered NASA-9 models.
It minimizes ideal-gas Gibbs energy at fixed temperature and pressure subject
to the actual feed elemental inventory and zero net electrical charge. Its HP
closure solves the same composition problem inside an explicit temperature bracket and closes total
enthalpy. A caller may supply a common inlet temperature or the sum of the
original enthalpies of separate inlet streams. Source formation enthalpy enters
once; there is no additional heat-release term.

Damped log-amount Newton steps retain positive species and distinguish a
structurally absent feed element from a representable logarithmic trace whose
direct fraction underflows. The ordinary constraint-space system has at most 18
unknowns. A bounded full KKT fallback, allocated lazily once, avoids squaring
trace-species conditioning when the reduced system loses rank. Scratch buffers
are reused across iterations and warm-started HP temperature trials. Species,
elements, Newton/backtracking steps and temperature work all have explicit
bounds; cancellation is polled before result publication.

Results retain accepted original-source receipts, element, charge and chemical-potential
residuals, independent mass closure, amounts, logarithmic mole fractions,
enthalpy and iteration counts. The optional equilibrium heat capacity includes
composition response through a nonnegative Gibbs-Hessian quadratic form; it is
not the frozen-composition heat capacity and does not imply `cv = cp - R` for
reacting products. The material freeze boundary must handle its own exact-one
floating-point composition convention and must not silently erase an
underflowed positive trace.

Positive and negative ions, multiple ionization stages, and one nucleus-free
charge -1 electron may share the same source bank. Feed must be neutral to
floating-point summation precision; net-charged feeds require a field model
and refuse. Charge coefficients are exact representable signed integers.
Charge balance and its differentiated constraint use logarithmically scaled
carrier amounts, independently of the neutral inventory. Thus even a charge
population whose direct amounts underflow cannot pass with an uncontrolled
electron/ion ratio. No trace cutoff or electron floor is introduced. A charge
sign without an eligible countercharge is structurally excluded.

Electrons contribute particle pressure, mixing entropy and energy through
their original NASA-9 sources. The formation-energy difference supplies
ionization heat once. Fixed-volume pressure brackets account for the extra
electrons using the nuclear inventory and positive charge bound. Equilibrium
Cp, Cv, compressibility and sound speed differentiate these same elemental and
charge constraints. Accepted endpoints share immutable original source models;
`is_bound_to` checks the complete model bank and canonical element/charge
identities before a consumer interprets the species as carriers.

Charged equilibrium assumes a dilute, nondegenerate, common-temperature local
thermodynamic equilibrium. Source applicability remains explicit; no Coulomb
correction, continuum lowering, pressure ionization, space-charge field,
two-temperature kinetics or radiation population balance is inferred.
Neutral Wilke/WMS transport is not a consequence of charged equilibrium.
Focused G1/G3 tests use independent Sackur-Tetrode/Saha hydrogen formulas for
TP/TV populations, electron pressure, ionization energy, HP/UV cooling and
equilibrium caloric/acoustic response, plus negative/multiply charged species,
trace charge balance, feed scaling and source-binding refusals.

`equilibrium::composition::FrozenIonizedGasStateV1` is the separate
actual-population boundary for transport and energy consumers. It keeps dense
canonical amounts, volume, the complete immutable source bank, heavy temperature
Th and electron temperature Te. It never solves equilibrium or creates a seed.
The nucleus-free charge -1 row identifies the electron independently of its
name. Quasineutrality is checked against actual charged amounts, including
subnormal traces, independently of the neutral inventory. Net charge requires
a different field-owning state. Positive populated-source Cv and the exact
source temperature domain are checked before publication.

Pressure is `R/V * (sum_heavy n_i Th + n_e Te)`; internal energy is the sum of
the original species `n_i u_i(T_i)`, with formation references counted once.
The two pool energies and frozen heat capacities are retained. Heavy internal
modes are at Th and the electron distribution is Maxwellian at Te; this does
not establish a two-temperature reaction law or chemical equilibrium.
Absent species impose no source-temperature domain on an actual frozen state;
a consumer creating products must admit their sources separately. Original
element, charge, source-model and coefficient arrays are shared immutably, as
are cloned frozen populations. The state contains matter energy only.

This closure does not infer a chemical mechanism,
authenticate a source-to-species declaration, determine phase stability, or
implement condensed phases, finite-rate combustion or a real-gas EOS.
Residual tolerances bound the equations, not relative errors of arbitrarily
tiny trace species. Manufactured analytic fixtures are not a CEA oracle match.

Native verification with the repository-pinned compiler ran the complete
`fs-thermochem` library: 89 tests passed, including 12 new composition/HP tests.
They cover analytic pressure-dependent dissociation, minor-element balance,
strong-product KKT fallback, trace underflow, independent adiabatic energy
closure, distinct inlet temperatures, equilibrium heat-capacity derivatives,
feed/reference invariance, source and work-budget refusals, and cancellation.

### Fixed-volume and adiabatic internal-energy equilibrium

`equilibrium::composition::volume` adds TV and UV endpoints through the same
original NASA-9 and element/charge-constrained Gibbs owner. `equilibrate_tv`
determines the pressure satisfying `p V = N(T,p) R T`; it does not hold the feed pressure
when reactions change molecular count. Conserved elemental amounts bound the
pressure search, including electron production by ionization. Safeguarded
Newton steps use the existing analytic reacting compressibility, with a bounded log-pressure bracket and an explicit relative
volume residual. One workspace, including its lazily allocated trace-species
KKT fallback, is reused through every trial.

`equilibrate_uv` derives total feed internal energy from the original species
`u = h - R T`. `equilibrate_uv_energy` accepts that target explicitly, including
separate-temperature feeds or declared heat/work contributions. Both conserve
volume, elements and charge while solving temperature and pressure; they neither reuse
constant-pressure enthalpy nor add formation energy a second time. The UV
Newton slope is the equilibrium constant-volume heat capacity, including
chemical composition response. Every eligible source must support the complete
temperature bracket; traversing NASA regions requires the existing explicit
caloric-continuity policy. Absolute energy and volume tolerances remain separate.

Only accepted endpoints allocate source receipts. Results retain the full
chemical state and thermodynamic response, original source identity, total
internal energy, volume/energy residuals and actual pressure/composition work.
Pressure and temperature trial counts are explicitly bounded, cancellation
prevents publication, and no supplied chemistry or feed is mutated. G1/G3/G4
regressions compare dissociation pressure to an independent closed form, rigid
exothermic energy balance to a separate scalar root, reacting Cv to fixed-volume
energy differences, and feed/reference scaling, source continuity, work budgets
and cancellation. These are dilute ideal-gas numerical closures, not a
finite-rate combustion, ignition, detonation, real-gas or phase-stability model,
and the manufactured fixtures are not an external CEA validation.

## Elementary finite-rate kinetics algebra

`kinetics` binds up to 512 caller-authored elementary reactions, with an
explicit `ForwardRateV1::Arrhenius` or `ForwardRateV1::ElectronImpact`,
to the original NASA-9 bank and the exact `fs-qty`
conservation certificate. Separate integer reactant/product orders retain
explicit catalysts. Molecular-mass consistency is independently checked against
the caller's bounded tolerance; species inventories are never repaired.

Rates use concentrations normalized by 1 mol/m3. Reversible-channel equilibrium
constants include the original `p0/(R T)` concentration standard, and analytic
concentration/temperature Jacobians include changing molecular count and exact
zero species. The Gibbs-derived affinity is independent of a disabled rate or
absent catalyst; an indeterminate both-sides-absent affinity is explicit.
Species sources are evaluated once per trial into reusable bounded scratch;
accepted rate evaluations allocate no heap storage. Near-equilibrium progress
uses affinity and `expm1` for reversible channels, with their ideal-mixture
entropy-production boundary limits retained. Original source thermodynamics supplies the adiabatic
instantaneous temperature rate, with formation energy included once.

`KineticWorkspaceV1::workspace_bytes` exposes actual allocated numerical array
capacities, including optional affinities, so higher-level coupled memory
budgets include this original owner. Immutable source storage is excluded.

Seven new native regressions cover Arrhenius scaling and heat release,
dissociation equilibrium against the Gibbs owner, zero species and catalysts,
400-decade rate scales near equilibrium, the complete constrained-UV Jacobian,
source/conservation/range/cancellation refusals and deterministic scratch reuse.
All 111 `fs-thermochem` library tests passed through the existing scoped native
Cargo manifest using the full original crate source. This is not a workspace
or cross-ISA qualification. The primary rate convention follows the
[Cantera 3.2 elementary-reaction and modified-Arrhenius reference](https://cantera.org/3.2/reference/kinetics/reaction-rates.html).

This L1 module owns immutable laws and algebra. It supplies no time integrator,
empirical mechanism inference, unsupported pressure laws, fractional-order
global reactions, or real-gas kinetics. A supplied
mechanism's physical accuracy remains a property of its original source data.

## Electron collision cross sections and kinetic rates

`electron_collision` owns source-retained incident-energy cross sections in
SI, with distinct ionization and elastic momentum-transfer meanings. A table
declares piecewise-linear interpolation, its complete threshold-to-maximum
energy coverage, and a nonnegative upper envelope for its unresolved tail.
Missing high-energy data are never implicitly zero. Immutable tables are
shared across evaluated sources. Independent implementations of NIST SRD 107's
generalized BEB Eq. (1) provide neutral H/He single-ionization sources with the
published orbital constants and dipole Q, explicit 10 keV model limit, citation
and original compilation terms. They are theoretical models, not experimental
validation or a license to redistribute NIST's compiled database.

The bounded, allocation-free Maxwellian evaluator integrates event
`<sigma v>` or displaced-Maxwellian fluid drag `(2/3)<sigma_mt v E/(kB Te)>`.
These are different moments: a constant momentum cross section gives
`(4/3) sigma <v>`. The latter is not the exact energy-dependent
Lorentz-Boltzmann mobility. Analytic temperature derivatives use the same
integral. Threshold shifting preserves finite log rates after direct
underflow; unrepresentable direct coefficients are `None`. The shift is retained
inside table interpolation and the BEB logarithm, so a positive threshold excess
is not rounded away by reconstructing absolute electron energy. Each admitted
table interval is reused without a new search at every quadrature node. Analytic
linear-threshold moments and H/He BEB threshold asymptotics test the factored
integral as well as its rate and temperature derivative. Positive Gauss
rules estimate quadrature error for both required moments, while an explicit
source envelope bounds the omitted tail. Exhausted work or unresolved tails
refuse. Numerical estimates and conditional model-tail bounds do not certify
experimental uncertainty or non-Maxwellian physics.

`ElectronImpactRateV1` feeds the existing exact stoichiometric kinetics owner
with `NA <sigma v> c_e c_target`, for a binary electron/neutral single-ionization
channel. It verifies the source target, structural electron, physical electron
mass, target-speed approximation, exact reaction orders and absence of a
pressure multiplier. No arbitrary Arrhenius surrogate is constructed. The
original common-temperature Gibbs source still sets the reverse coefficient;
this does not evaluate molecular detailed balance at a separate Te. Both
zero-electron absence and the conservative stiff gas reactor are preserved.

Primary model equations and parameter sources:
https://physics.nist.gov/PhysRefData/Ionization/Eqs/latex.html,
https://physics.nist.gov/cgi-bin/Ionization/bun.php?Id=HeI,
and Banks (1966), https://doi.org/10.1016/0032-0633(66)90024-9.
Analytic and independently evaluated numerical moment tests check source
thresholds, derivatives, molar conversion, detailed balance, coverage and
refusals; they do not validate an observed plasma history.

### Shared caloric-path admission for evolving gas

`kinetics::admit_caloric_interval` is the common L1 admission used by the
finite-rate reactor and frozen compressible Euler owner in `fs-gas`. It checks
original source support and explicit continuity across NASA regions, and uses
bounded outward-rounded polynomial subdivision to establish strictly positive
species Cv throughout the requested interval. Positive endpoint samples alone
do not admit a path. Unresolved bounds, gaps, discontinuity or cancellation
refuse without changing source coefficients or formation/entropy references.
The subdivision uses a fixed-size stack with at most 4095 work items per region.
The native regression includes positive endpoint heat capacities with a
negative interior Cv. This addition brings the full-source native library
suite to 111 passing tests. Time integration remains in L3.

## Pressure-dependent detailed-balance reaction rates

Kinetic operation-tree version 2 adds explicitly authored
`ThirdBodyEfficienciesV1`, `PressureDependenceV1`, `FalloffBroadeningV1` and
`TroeFalloffV1`. Elementary reactions remain pressure independent unless the
caller attaches a pressure law. Species-specific collider efficiencies are
canonical, finite, nonnegative and bound to the exact original species bank;
unknown or duplicate species are refused. Their explicit default applies to
unlisted species.

The effective collider concentration is `M = sum(epsilon_i c_i)`. Third-body
rates multiply the elementary coefficient by `M/c*`, with `c*=1 mol/m3`.
Normal falloff uses `k = kinf Pr/(1+Pr) F`, `Pr=k0(M/c*)/kinf`, with Lindemann
or three-/four-parameter Troe broadening. The original forward Arrhenius law is
the high-pressure law and the additional low-pressure law includes one extra
reference concentration. Coherent-SI conversion remains the mechanism author's
responsibility.

Both directions use the same effective coefficient, preserving the original
NASA-9 equilibrium constant, affinity and exact element/charge accounting.
Analytic concentration and temperature Jacobians include collider dependence,
both Arrhenius limits and Troe derivatives; existing stiff reactors consume
these derivatives directly. Logarithmic evaluation avoids artificial overflow
in `k0*M` and underflow in reduced pressure. Exact zero collider density has its
correct inward derivative, including Troe's actual limiting broadening. The
scaled Troe rational form remains finite at the conventional apparent pole.
Zero authored rate scales or all-zero efficiencies bypass irrelevant rate
arithmetic and retain exact kinetic shutdown. Successful evaluations allocate
nothing; collider bindings are prepared once.

Seven added algebra regressions cover rates and all Jacobian columns, pressure
limits and extreme intermediate values, Troe values/derivatives, exact zero
colliders, the apparent pole, near-equilibrium entropy, declaration refusals and
shutdown. The complete native thermochemistry library passed 118 tests; the
expanded shutdown regression passed again after its final guard correction.
An additional prepared-reactor regression matches the analytic density-dependent
third-body isomerization relaxation at two physical volumes.

This admits constant collider efficiencies and normal Lindemann/Troe falloff.
Chemically activated reactions, SRI, P-log, temperature-dependent efficiencies
and mechanism calibration remain outside this closure. Troe alpha is admitted
in [0,1] with positive finite characteristic temperatures. References:
<https://cantera.org/3.2/reference/kinetics/reaction-rates.html> and
<https://cantera.org/3.2/reference/kinetics/rate-constants.html>.

## Ordinary-water saturation source (IAPWS SR1-86(1992))

`water::SaturatedWaterState` evaluates equations (1)–(9) of the IAPWS
Revised Supplementary Release on Saturation Properties of Ordinary Water
Substance (September 1992, ITS-90):
<https://iapws.org/public/documents/6dGkr/Supp-sat.pdf>.
Its inclusive source range is 273.16–647.096 K. The immutable state carries
saturation pressure, analytic dp/dT, both phase densities, and a common-source
enthalpy, entropy, and internal-energy datum. It is a saturation source, not
an off-saturation equation of state or a phase-change time integrator.

The implementation shares deterministic fractional powers and evaluates the
phase volume gap before subtractive cancellation. Latent quantities obey
Delta h = T(dp/dT)Delta v, Delta s = (dp/dT)Delta v, and
Delta u = (T dp/dT - p)Delta v. The exact critical state has zero phase gaps;
consumers requiring distinct phase inventories must reject that degeneracy.
Published coefficient precision, including the small triple-point reference
residuals, is retained without an extra datum shift.

Pressure inversion uses the original correlation with safeguarded Newton and
bisection, a caller-declared absolute pressure residual budget, bounded work,
and cancellation. The returned source pressure is never replaced by the target.
Admission uses evaluated endpoint pressures rather than rounded table entries.
Numerical solve tolerances are not claims about the release's material uncertainty.
Tests cover every published reference-table row at its stated precision,
independent pressure differentiation, Clapeyron/Gibbs/pressure-work identities,
critical coalescence, inversion, domain limits, work bounds, and cancellation.

### Ordinary-water/pure-vapor surface tension

`water::surface::WaterSurfaceTensionState` evaluates the unchanged IAPWS
R1-76(2014) ITS-90 correlation, attributed to IAPWS:
<https://iapws.org/technical-guidance/release/Surf-H2O.download>.
It returns typed N/m surface tension and its analytic N/(m K) derivative along
saturation for 273.16--647.096 K. One deterministic fractional power supplies
both values without allocation. The critical value and limiting derivative
are exactly zero; this does not establish a distinct interface there.
`SaturatedWaterState` and `If97SaturationPair` expose the same source at their
existing temperature, without changing saturation pressure or caloric datums.
IF97's sub-triple extension does not widen this interfacial source domain.
The source is ordinary water in equilibrium with pure water vapor; arbitrary
carrier gases, solutes, surfactants, wetting, curvature corrections and surface
heat/storage closures remain separate. G0/G1 tests compare literal Table 1
values, independent temperature differences, exact domain boundaries and
adjacent-to-critical temperatures. Table-rounding tolerances do not describe
the source's experimental uncertainty.

## Maxwell-Stefan molecular diffusion with explicit reference frames

`transport::maxwell_stefan` supplies a bounded L1 constitutive friction solve,
without evolving any spatial field. `MaxwellStefanV1` retains canonical original
species IDs and molar masses. Its reusable workspace accepts actual mole
fractions, their spatial gradients, T, p, grad(log p), and a complete symmetric
positive binary-diffusivity matrix in coherent SI.

For `z_i=j_i/M_i`, it solves
`sum_j (x_j z_i-x_i z_j)/D_ij = -c [grad(x_i)+(x_i-Y_i) grad(log p)]`
with the barycentric constraint `sum_i M_i z_i=0`. `evaluate` returns mass fluxes
relative to mass-averaged velocity. `evaluate_molar` instead imposes `sum_i z_i=0`
in that same friction solve, with distinct molar controls and result types;
its mass fluxes generally do not sum to zero. The existing material parcel
diffusion path consumes this molar-frame owner while retaining its original
source queries, isothermal/isobaric admission and finite-volume exchange checks.
Unequal molecular masses, coupled
crossdiffusion and barodiffusion are retained. Exact zero mole fractions are
admitted without division by x, seeding or flux clipping.

`evaluate_driving_forces` and `evaluate_molar_driving_forces` admit explicit
source-owned dimensionless thermodynamic forces `d_i [1/m]` and typed total
molar concentration, using those same friction rows and frame constraints.
The force sum must close the caller's Gibbs–Duhem budget; incompatible data
are refused rather than projected. Chemical-potential derivatives, activity
coefficients, EOS concentration, phase stability, thermal coupling and reciprocal
energy terms remain the source owner's responsibility. The API does not infer
any of them. The positive symmetric pair model retains local nonnegative drag
entropy; it does not establish a complete nonideal or thermodiffusive entropy
balance without those source closures. Independent fixtures use a stable binary
regular-solution chemical potential, both reference frames, exact ideal
composition/pressure reduction, invalid-force refusal and cancellation/reuse.

The maximum is 128 species. Scaled pivot elimination reuses fixed-dimension
matrix and output buffers. Every original friction row, including the gauge
row replaced by the chosen frame constraint, and the actual published primary
flux sum must close their explicit absolute/relative budgets before return.
Mass flux is primary in the barycentric frame; molar flux is primary in the
molar frame and does not disappear merely because conversion to a much smaller
mass flux underflows. Cancellation
is polled at bounded assembly/elimination/check boundaries. Successful
evaluations allocate no storage. The optional local pair-friction entropy
production is absent at an exact-zero composition carrying nonzero flux; no
finite chemical potential is invented there.

Seven focused native tests pass: unequal-mass binary analytic flux; independent
ternary crossdiffusion (including flux at zero own-species gradient);
manufactured unequal-mass ternary original equations; exact-zero and structural
species; analytic pressure diffusion; canonicalization/scaling/buffer reuse;
and maximum-species, invalid-data and late-cancellation behavior. Built-in thermal
(Soret) coefficients, external species forces, ions, nonideal activities and a
fully discrete entropy theorem remain outside this algebra. Pair coefficients
and their thermodynamic-state support remain caller-owned source data.
Additional focused regressions cover the unequal-mass binary molar Fick limit,
the independently reconstructed mass-to-molar reference-velocity transformation,
frame-invariant pair-friction entropy, subnormal primary molar flux and reuse
after cancellation. The frame change does not add physical transport, thermal
diffusion, an empirical coefficient, or a fully discrete entropy claim.
Primary formulation:
<https://doc.comsol.com/6.4/doc/com.comsol.help.chem/chem_ug_chemsptrans.08.257.html>.



### Exact composition reduction for sparse species populations

Evaluator version 3 solves an original two-species system analytically. For a
larger species bank, each exact-zero-composition row is diagonal and is solved
first. Its potentially nonzero entering flux remains in the positive-composition
rows and in the chosen mass- or molar-frame constraint. Only the k species with
positive fractions enter the coupled elimination. A single carrier closes the
frame directly. There is no concentration floor, deleted species, zero-force
assumption for absent species, or approximation to their effect on the carriers.

Every supplied pair still must be finite, positive and exactly symmetric.
Absent/absent pairs multiply identically zero friction terms, so they no longer
set the internal diffusivity scale. Their source data remain required and become
active immediately when the corresponding composition changes. This also avoids
an irrelevant tiny pair coefficient forcing an otherwise representable response
outside the numerical scaling range. Version 3 retains version 2's balanced
binary-exponent flux reconstruction and all original-row/frame residual checks.

Assembly, full-data validation and final equation checks remain O(n²); coupled
elimination becomes O(k³). Existing O(n²) heap workspaces are reused, with bounded
fixed-size stack scratch and no successful-call allocation. These are operation
and storage bounds, not a measured speedup. Existing gas, material and liquid
consumers inherit this exact algebra through the original constitutive owner.

Six added G0/G1/G3 groups cover independently manufactured one-, two- and
many-carrier mixtures through 128 source species in both frames; nonzero entering
flux for absent species; irrelevant pair scales down to the smallest subnormal;
positive traces whose friction remains finite; pure/binary extreme mass ratios;
composition changes, source-data refusals and cancellation/workspace reuse.
The exact production module and its genuine fs-qty/fs-blake3 dependencies passed
all 22 native leaf tests. This focused result is not a workspace/DSR result or
experimental transport qualification.

### Equal-pair diffusion without a cubic solve

Evaluator version 4 specializes exactly equal active pair diffusivities.
Unequal pairs retain the original coupled solve; absent/absent pairs retain
their complete source admission without setting this physical coefficient.
For the scaled response, the common friction matrix is S*I-x*1^T, S=sum(x).
Its frame-constrained solution uses the linear formula
z=(f-x*(w.f)/(w.x))/S, with the original omitted force row and frame weights.

Near subtractive cancellation, that same numerator is evaluated through
off-diagonal pairs w_j*(x_j*f_i-x_i*f_j). The self-pair cancels algebraically;
balanced products retain representable trace responses under extreme molar-mass
contrasts. The original frame closes the final component directly. No species,
force or flux is clipped, and all published original-row/frame checks remain.

Ordinary solution algebra is O(n), and the bounded pair evaluation is O(n^2)
in the worst case. Complete admission and original-equation checks remain
O(n^2); reusable heap storage and bounded stack scratch require no successful
evaluation allocation. These are operation bounds, not a measured speedup.
Two new G1 groups recover dense 3–128-species fluxes in both frames across
240 orders of pair diffusivity, and retain extreme-mass and subnormal trace
responses. All 24 native tests of the real module and its genuine fs-qty and
fs-blake3 dependencies passed. This is a focused module result, not a full
workspace/DSR result or experimental transport qualification.


### Reusing unequal-pair friction across spatial directions

Evaluator version 5 retains a factorization of the exact assembled active
friction and reference-frame matrix in the existing workspace. A changed
dimension or any changed coefficient bit refactors it; neither source IDs nor
approximate fingerprints can select stale factors. Each direction still admits
the complete current source table and fractions, forms fresh absent-species
responses and driving terms, rescales the published physical flux, and checks
every original friction equation and the chosen frame.

The first unequal-pair solve costs O(k^3) for k positive-composition species;
an unchanged matrix then costs O(k^2) per driving direction. Row scaling,
pivot choices and chronological right-hand-side arithmetic retain the original
scalar elimination order. The two additional dense buffers and pivot indices
are included in `workspace_bytes`; no successful evaluation allocates.
Interrupted refactorization invalidates the cache before mutation, while an
interrupted right-hand-side solve leaves complete factors reusable.

Six focused native groups cover dense, exact-absent and positive-trace mixtures
through 128 species in both frames; bitwise equality with fresh workspaces;
composition, coefficient, frame and molecular-mass changes; current source
admission and diffusivity rescaling; late pivots; cancellation/reuse; and
capacity accounting. All 30 tests of the real constitutive module and genuine
fs-qty/fs-blake3 dependencies passed. This is a focused native module result,
not a full workspace gate or experimental transport qualification.

### Analytic directional diffusion response

The mass- and molar-frame driving-force APIs also expose analytic tangent
evaluations for one dimensionless parameter. The caller supplies derivatives
of composition, thermodynamic driving forces, total concentration and every
binary diffusivity in the same original canonical order. Species and molecular
masses stay fixed. Derivative composition closes the simplex tangent plane;
force derivatives close Gibbs–Duhem and pair derivatives are finite symmetric
with zero diagonal. Physical source differentiability, validity and one-sided
paths at exact-zero composition remain explicit caller responsibilities.

The tangent differentiates the original equations,
`A z' = -c d' - c' d - A' z`, with the same reference-frame constraint.
It evaluates both right-hand sides through the original scalar owner and
reuses its exact matrix factors. A separate check reconstructs every
differentiated physical row from the published value and derivative fluxes.
The caller's original row/frame controls apply independently to both solves.
No finite difference, concentration floor, source extrapolation or entropy
derivative is inferred.

The dedicated workspace adds two O(N) retained primal arrays; ordinary scalar
storage and its operation tree are unchanged. All allocated capacities are
reported, cancellation publishes no partial result, and successful evaluation
does not allocate. The tangent operation has its own evaluator version 1.

Six focused native groups cover a nonideal binary analytical derivative,
manufactured original-row responses through 128 species in both frames with
exact absence and positive traces, an independent interior finite-difference
convergence ladder, zero/trace limits, derivative/source refusals, linearity,
cancellation and storage reuse. All 36 real-module tests passed with genuine
fs-qty/fs-blake3 dependencies. This is constitutive verification, not a coupled
nonlinear-solver or physical-source qualification.

## Source-bound dilute-gas collision transport

`transport::kinetic::KineticGasTransportV1` retains the exact original NASA
species bank, explicit model ID/version, canonical Lennard-Jones diameters and
well depths, and an authored T/p domain. Molecular masses come exclusively from
the original species sources. First Chapman-Enskog viscosity and binary
diffusion use exact SI Boltzmann/Avogadro constants, Lorentz-Berthelot pair
combining rules, and the complete Neufeld-Janzen-Aziz collision-integral fits.
The declared domain must remain within 0.3 <= T* <= 100 for every pure/pair
collision law. No dense-gas or polar correction is silently inferred.

Pair molecular factors are precomputed. A reusable workspace provides the
symmetric canonical binary-diffusion matrix (zero diagonal), pure viscosities,
pure conductivities, mixed viscosity and thermal conductivity without
successful-evaluation allocation. Only Omega11 is evaluated in the quadratic
binary-diffusion loop. Mixture viscosity follows Wilke and is checked against
the existing shared Wilke owner.

Conductivity is explicitly the Eucken estimate
`lambda_i = mu_i [Cp_i(T)+5R/4]/M_i`, with Cp from the original NASA law, followed
by the arithmetic/harmonic mixture rule. It recovers monatomic Pr=2/3 without
imposing a constant polyatomic Prandtl number. Polyatomic internal relaxation,
Soret/Dufour transport, polar interactions and experimental validation are not
claimed by this estimate.

The default classical translational admission remains `Cp >= 5R/2`, with only
its original 32-epsilon arithmetic tolerance. An explicit consuming
`with_translational_cp_deficit_allowance(MolarThermalQuantityV1)` declaration
can admit a small source-fit deficit in J/(mol K). The retained accessor and
model equality distinguish that declaration from strict mode. The caller must
justify it against the original fit and physical domain; every actual Cp is
still checked, including possible products with zero current fraction.
Allowances must be finite, nonnegative and no larger than one part per million
of `5R/2`. This hard policy ceiling prevents this small-fit control from
waiving a substantial physical deficit; it is not an empirical uncertainty
bound. Zero restores exactly the original strict admission. The allowance
changes neither the source bank nor its Cp: Eucken conductivity continues to
consume the actual source value, without a floor or any coefficient edits.

The retained Cantera hydrogen source supplies a concrete case: its original H
high-temperature NASA fit has a deficit of about `1.4445e-8 J/(mol K)` near
1456.31 K. Two focused regressions exercise explicit `2e-8 J/(mol K)` admission
against that actual source, strict and insufficient-budget refusal, unchanged
source identity and exact original-Cp Eucken evaluation, invalid/overlarge
allowance refusal, and continued rejection of a physical `Cp=2R` source even
at the maximum permitted allowance. This opt-in does not validate the Eucken
approximation, prove a source-fit error bound, or extend the collision domain.

Nine native tests pass on actual source: the two Cp-admission cases above;
independent SI prefactors; complete
collision-correlation values and admitted domain; monatomic Prandtl/Schmidt
relations; pure/trace/Wilke limits; pressure, diameter and temperature scaling;
original variable-Cp conductivity; and source identity, canonicalization,
cancellation and workspace reuse.
References: Neufeld, Janzen and Aziz, J. Chem. Phys.57,1100-1102(1972),
<https://doi.org/10.1063/1.1678363>;
<https://cantera.org/3.2/cxx/d5/de4/GasTransport_8cpp_source.html>;
<https://cantera.org/3.2/cxx/dd/d16/MixTransport_8cpp_source.html>.

## Original IF97 compressed-liquid and steam source

`water::if97` implements the original IAPWS industrial formulation in native
Rust: Region 1 liquid, Region 2 vapor, Region 5 high-temperature steam,
Region 4 saturation and the B23 boundary.
The source is the [IAPWS R7-97(2012) revised release](https://iapws.org/technical-guidance/release/IF97-Rev.download).
The released coefficients and equations are attributed to IAPWS; this is a
separate source from the existing SR1-86(1992) saturation formulation. They are
never silently combined into one caloric chart.

States retain typed pressure, temperature, density, specific volume, enthalpy,
internal energy, entropy and Gibbs energy, together with Cp, Cv, sound speed and
analytic thermodynamic derivatives. Admission follows the original regional
bounds through 100 MPa and 1073.15 K, followed by the revised Region 5 domain
through 50 MPa and 2273.15 K. The Region 5 equation also admits its common
1073.15 K boundary when selected explicitly. Its original six ideal and six
residual terms (equations 32–34, Tables 37/38) feed the same Gibbs derivative
algebra with smaller fixed stack power tables; no Region 2 extrapolation or
constant-Cp substitution is used. The source assumes pure **undissociated**
water, so chemical dissociation requires a separate model. The separate
metastable-vapor formulation remains absent. Region 3 is supplied by the
distinct density/temperature source documented below. The source's 273.15–273.16 K extension
includes the metastability explicitly described in the release; it does not
establish stable liquid water against ice. Saturation pressure and temperature
extend to the critical point, while paired Region 1/2 saturated calorics stop
at 623.15 K. Original regional Gibbs mismatch is exposed, not repaired by
shifting source caloric values. Canonical closed saturation and B23 endpoints
are admitted without widening a regional domain by a tolerance band.

`bounds::isobaric_heat_capacity_bound` proves a positive Cp enclosure over a
caller-selected regional temperature interval using outward-rounded arithmetic
and bounded adaptive subdivision. The result bounds the retained binary64
coefficient model, not experimental uncertainty or source-fit error. It exposes
the inverse-temperature slope bound and work count; unsupported paths, exhausted
budgets and cancellation refuse explicitly.

`increments::isobaric_difference` evaluates changes of h, s, v and u on one
original-source isobar. Cached integer-power difference recurrences preserve the
input temperature increment before subtracting reduced temperatures or large
absolute caloric values. The first-law relation is delta_u=delta_h-p*delta_v.
Both endpoints and the complete regional path are admitted; cancellation is
polled and no foreign numerical runtime is used.

Region 5 participates in these original-source heat-capacity enclosures and
finite caloric/volume increments, and in both density-jet APIs, including the
complete pressure/temperature Hessian. Its pressure curvature uses the original
residual terms to preserve the dilute-gas limit. The extended focused cases
cover all three published Table 42 verification points, independent Maxwell
and constant-volume Cv derivatives, finite Cp/Cp-over-T integrals, adjacent
binary64 temperature increments, full-domain Cp enclosures and source limits.

The 17 focused tests cover published verification tables, independent
thermodynamic derivative identities, density anomaly and Cv behavior, raw
saturation Gibbs mismatch, closed source boundaries, whole-interval Cp bounds,
adjacent-representable-temperature increments, independent Cp/Cp-over-T
integrals, reversal, work refusal and cancellation. The merged native
thermochem source suite passes all 156 tests. This does not claim validation of
unsupported regions or a full-workspace Cargo/DSR build.

## Maxwell--Stefan scaling across trace gradients

Fixed operation-tree evaluator version 2 introduced forming each mass flux with its solved
amplitude before rounding a scale. Normal intermediates retain an ordinary
arithmetic path; otherwise finite products and quotients balance exact binary
exponents and round only the final normal/subnormal magnitude. This avoids
rejecting a representable species flux because the common
`concentration * minimum_diffusivity * driving_scale` underflowed first.

Published mass flux reconstructs the normalized coefficients for the original
friction-row and mass-frame checks. A naturally zero-rounded flux is admitted
only if those unchanged budgets pass. There is no gradient cutoff, species
clipping or tolerance increase, and true overflow still refuses. The fixed
operation-tree version changes because operation ordering can change rounding.

All ten exact-production-module tests pass. Three new independent cases cover
normal/subnormal Fick flux below the absolute residual tolerance, genuinely
unrepresentable flux with both admitted and stricter-refused residual budgets,
and an inactive `D12=1e-300` pair whose common normalization underflows while
the actual `1e-30` mass flux remains representable. The correction addresses
the concrete refined N2/O2 molecular-manifold refusal, not a model calibration
or new entropy/stability certificate.

## Original IF97 Region 3 dense water

`water::if97::region3` extends the original IF97 source with the Region 3
density/temperature Helmholtz equation, using equation 28 and Tables 30–33 of
the attributed [IAPWS R7-97(2012) release](https://iapws.org/technical-guidance/release/IF97-Rev.download).
`If97DenseState` returns typed p, rho, T, v, h, u, s, Helmholtz/Gibbs energies,
Cp, Cv, sound speed and analytic pressure/caloric/volume derivatives. Cached
integer powers and one coefficient pass evaluate the source without
allocation or a nonlinear density solve. Existing Region 1/2 APIs retain their
original operation order and references.

Admission uses 623.15–863.15 K, the original B23/100 MPa pressure limits,
positive Cv and isothermal stiffness, and the stable Region 4 side below Tc.
A conservative 805 kg/m3 density envelope excludes remote extrapolated
polynomial roots; independent outward interval arithmetic verifies that its
upper boundary lies beyond 100 MPa over the complete source temperature
interval. The exact critical state refuses its singular response functions.
No critical neighborhood is smoothed, and no pressure, entropy or energy is
shifted to join independently fitted regions. The direct density/temperature API does not itself supply a mixture lever rule,
ice or Region 5; pressure inversion and high-temperature saturation endpoints
are supplied by the separate module documented below.

`region3::increments::DenseIsochoricInterval` prepares an immutable fixed-density
temperature path. Bounded outward interval subdivision establishes its complete
source-domain membership, stable derivatives and positive Cv enclosure, including
B23, pressure and subcritical saturation-side constraints. Coefficients are
collected once at the fixed density and evaluated with interval Horner arithmetic.
Explicit work exhaustion, unresolved boundary inclusion and cancellation refuse;
the proof does not widen physical bounds by a tolerance band. Enclosures concern
the retained binary64 coefficient model, not source-fit or experimental uncertainty.

Repeated `difference` queries reuse this proof and factor the original
integer-power changes before subtracting large absolute properties. Original
temperature differences survive adjacent-representable endpoints; delta_h equals
delta_u+delta_p/rho. Changes of entropy, pressure and both caloric potentials
remain on the original datum, with exact reversal and zero-change behavior.
The convenience free function has a documented 4095-interval setup ceiling;
repeated consumers should retain the prepared interval.

All 169 tests pass in the current native thermochemistry source suite. The ten
new tests include all published Table 33 states at their printed precision;
independent thermodynamic derivatives, isobaric Cp and isentropic acoustic paths;
stable/metastable/critical/domain refusals; the density-envelope proof; whole-path
Cv bounds; adjacent-temperature increments; independent Cv and Cv/T integrals;
additivity, reversal, bounded work and cancellation. The existing 37 unchanged
source files match the inspected main-tree blobs. Full workspace DSR remains
unavailable in this environment.

## Pressure-driven dense water and stable caloric changes

`region3::pressure` adds original-equation pressure/temperature inversion and
Region 3/4 saturation endpoints for 623.15 < T < 647.096 K. Stable outer roots
are selected using the original saturation equation; temperature powers are
cached once during safeguarded analytic Newton/bisection. The returned source
pressure, density, h, u and s remain unmodified. Paired endpoints expose both
pressure residuals, the requested-pressure/Region 4 inverse residual, latent
enthalpy and entropy, and the original Gibbs mismatch. A pair is not a bulk
mixture model or a new caloric reference.

Acceptance requires an explicit pressure tolerance and a numerical density
bracket width. An outward evaluation verifies positive Cv and isothermal
stiffness throughout that final bracket. The pressure-sign bracket is a
floating-point solver diagnostic, not a certified real-arithmetic root-error
enclosure or an experimental uncertainty. Per-endpoint work is limited to
1..=256 branch/root iterations, with bounded cancellation checks and no
allocation on success. Exact critical states, ambiguous saturation, unresolved
branches/stability, insufficient resolution, exhausted work and unsupported
source domains refuse. No critical neighborhood is smoothed or admitted by
widening a source boundary.

`region3::difference::difference` provides stable changes of original h, u, s,
v, p and p/rho between validated source states with changing density AND
temperature. Factored integer-power changes and the existing stable entropy
log-ratio primitive preserve small input increments. Compensated accumulation
and canonical endpoint ordering give fixed work, exact reversal and exact zero;
delta_h equals delta_u plus delta(p/rho). These are endpoint differences, not
evidence that an interpolated thermodynamic path is admissible. Consumers of
absolute storage must still audit represented endpoint energy and pressure work.

The native original-source thermochemistry suite passes 177 tests. Eight new
tests cover all Table 33 pressure inversions with printed-value precision;
liquid/vapor/supercritical round trips and derivatives; distinct high-temperature
saturation pairs; original pressure/Gibbs residuals; closed 100 MPa and B23
boundaries; work, domain, resolution and early/late cancellation refusals;
adjacent-representable density/temperature changes; independent first-law and
entropy integrals; and agreement with the existing prepared isochore.
There is no Region 3 mixture transport, automatic Region 1/2/3 caloric stitching,
ice or Region 5 claim. Full-workspace DSR/Cargo validation remains separate.

## Original water regional routing and source junctions

`water::if97::regions` selects the original stable Region 1, 2, 3 or 5 equation
from pressure and temperature. `If97FluidState` retains its original source
state and typed h, u, s, v, density, Cp, Cv, sound speed and derivatives.
Requested pressure, achieved dense-source pressure residual and density work
remain distinct. Gibbs evaluation introduces no density solve. Both canonical
Region 4 anchors, p(T) and T(p), are recognized as ambiguous saturation rather
than selecting a phase from inverse-equation roundoff.

`If97FluidSaturationPair` routes low-temperature coexistence to the original
Region 1/2 pair and high-temperature coexistence to the Region 3 pair. Original
latent h/s, Gibbs mismatch and saturation-pressure residual remain visible.
The exact minimum-temperature and 623.15 K pressure anchors use their original
forward-temperature constructors, preserving the requested pressure bit-for-bit.
Noncanonical inverse results below the stable ingress minimum refuse; no
temperature clipping or caloric adjustment is applied.

`If97RegionalJunction` prepares the two original source endpoints at a
single-phase regional boundary. Region 1/3 uses the 623.15 K isotherm;
Region 3/2 uses bounded bisection of the original forward B23 equation to obtain
admitted one-sided temperatures. Their finite width, work and pressure residuals
are exposed. The temperature bracket is a numerical admission diagnostic, not
a certified real-arithmetic root enclosure. A 100 MPa B23 endpoint pair is
admitted within the tested 1e-8 K caller budget without widening source domains.
Region 2/5 uses the exact 1073.15 K boundary through 50 MPa, with zero root
iterations. Automatic routing retains Region 2 at that boundary and selects
Region 5 above it; explicit source selection and the junction retain both
unadjusted source endpoints and their finite matching errors.

Signed original h/u/s/v/Gibbs differences at these junctions are retained.
Section 10/Table 43 of the attributed IAPWS release describes finite regional
matching errors; a junction is not a phase transition and its caloric mismatch
is not latent heat. No interpolation, entropy/enthalpy reference correction or
common shifted equation is introduced by this router. The engineering ingress
is 273.16–2273.15 K with original regional pressure limits through 100 MPa
up to 1073.15 K and 50 MPa above it.
Ice, metastable continuation and the exact critical singularity remain
unsupported; unresolved near-critical solves retain their explicit refusal.

Five focused tests verify regional source identity, canonical saturation
ambiguity and endpoints, positive/negative junction mismatches at 20/25/50/100
MPa, independent common pressure derivatives, source-domain and work limits,
and early/late cancellation. All 182 native thermochemistry tests pass.
Production routing is allocation-free and retains bounded source/boundary work.

## Pinned hydrogen/oxygen finite-rate mechanism

`kinetics::hydrogen::cantera_h2o2` constructs the complete ten-species,
29-reaction ideal-gas `ohmech` from Cantera 3.2.0, commit
`4a8358eb80cfeb50474386b5f9ec0b3a83519889`, mechanism blob
`65e6321a60bafb1feb4f3b1be5347e6d1e89c937`. Atomic weights retain source
`Elements.cpp` blob `ebe9da8b1a3e7944d1de070851b36fdb858ad8d4`.
The original typed species bank, NASA owner and kinetic laws perform all work;
no external runtime dependency or second thermodynamic implementation is added.
All duplicate reactions, explicit collider channels, five generic third-body
laws and the Troe falloff law remain. Source cm/mol/s and cal/mol rate data
are converted to SI using their actual reaction order.

The canonical axis is AR, H, H2, H2O, H2O2, HO2, N2, O, O2, OH. NASA7 is
embedded exactly in NASA9 with zero inverse-temperature coefficients, on the
common 300–3500 K interval. Callers explicitly choose original source h/s or
`ContinuousAt300K`: the latter preserves the low-temperature region and Cp,
using the existing continuity owner to bound upper-region h/s shifts by
0.01 J/mol and 2e-5 J/(mol K). These are numerical reconciliation budgets,
not experimental uncertainty. Source cards and corrections remain inspectable.
Reverse rates and equilibrium constants follow the same chosen thermodynamics.

Five focused source tests cover every species Cp/h/s and every reaction net
rate at 700 and 1800 K, source identities, pressure/collider channels, equilibrium,
continuity and cancellation. The native thermochemistry library passes 164 tests.
An independent development-only Cantera 3.2.0 calculation agrees with all
29 fixture rates to 2.45e-13 relative and thermodynamics to 3.51e-16.
This establishes source conformance, not experimental validation.
The loader supplies the hydrogen/oxygen submechanism, not full GRI methane
chemistry, a time integrator, transport, soot, real-gas phases or a fire model.

## Original ice Ih thermodynamics and melting

`water::ice` evaluates the original IAPWS R10-06(2009) Gibbs equation using
the released IAPWS-95-compatible caloric reference. It retains the original
complex poles and pressure-polynomial coefficients. A small native complex
helper and analytic Gibbs derivatives supply h, u, s, g, Helmholtz energy,
volume, density, Cp, Cv, expansion, pressure coefficient and compressibility.
No enthalpy or entropy datum is shifted to force agreement with another EOS.

The numerical chart is deliberately limited to 50..=273.16 K and
0..=210 MPa, explicitly assuming hexagonal ice Ih. This source evaluator does
not select the globally stable phase against fluid water or other ice
polymorphs. Consumers must admit their phase path and retain the release's
experimental limitations separately from arithmetic residuals.

The independent original R14-08(2011) Ih melting correlation supplies p(T),
its analytic derivative, and a bounded safeguarded T(p) inverse over
251.165..=273.16 K. Residual and temperature-width budgets must both pass.
Exact forward anchors at 251.165 K, 273.15 K and 273.16 K are retained, including
the lower-temperature overlap needed by the original IF97 liquid equation.
Reported temperature brackets concern evaluated floating-point correlation
values; they are not real-arithmetic root enclosures. Expm1/log1p expressions
retain the small pressure changes near the triple point.

Stable isobaric differences evaluate changes in h, u, s and v without
subtracting large reference quantities. The analytic increments preserve
small temperature steps, exact zero, reversal and subdivision, and maintain
Delta h=Delta u+p*Delta v.

`water::ice::bounds` prepares a whole-temperature-interval positive Cp
enclosure and a conservative dT/dh bound. It reduces the Gibbs second
temperature derivative to real rational expressions and rounds every
interval operation outward, including pressure-dependent source coefficients.
Adaptive subdivision uses explicit work/depth limits and fixed stack storage.
Samples never certify the bound; unresolved work or cancellation returns no
partial result. These are bounds on the retained source equation, not
experimental uncertainty or a global ice phase-stability certificate.

All 188 native thermochemistry tests pass: the 182 previous tests, four ice
source tests, and two whole-interval Cp-bound tests. Verification covers all
16 rows and three columns of the release's Table 6, independent thermodynamic
derivative identities, stable tiny/reversed/subdivided increments, the released
melting value, pressure/temperature boundaries and inverse budgets, and Cp
enclosures across the admitted assumed-Ih chart. Original source attribution:
https://iapws.org/technical-guidance/release/Ice-2009.download and
https://iapws.org/technical-guidance/release/MeltSub.download.

## Complete GRI-Mech 3.0 methane and nitrogen kinetics

`kinetics::gri30::cantera_gri30` loads all 53 species and 325 channels from
Cantera 3.2.0 commit `4a8358eb80cfeb50474386b5f9ec0b3a83519889`,
`data/gri30.yaml` Git blob `9a62004ec17bbec8ac11e97785932d426dd58610`.
Original GRI-Mech 3.0 attribution and the Cantera BSD redistribution notice
are retained. Atomic weights come from the same pinned `Elements.cpp`
source as the hydrogen mechanism. Source tables convert cm/mol/s and cal/mol
once, preserving species and channel identities, duplicate contributions,
nine explicit collision partners, twelve generic third-body laws, three
Lindemann and 26 Troe falloff channels. No runtime dependency is added.

`Gri30Calorics::Original` retains source NASA7 Cp/h/s and original seams.
`ContinuousAt300K` preserves low-region h/s and all Cp values, using the
existing continuity owner to bound high-region offsets by 0.2 J/mol and
2e-4 J/(mol K). The largest offsets are approximately +0.156299 J/mol for
C3H7 and +0.000148985 J/(mol K) for C3H8. These are numerical reconciliation
budgets, not empirical uncertainty. Original species intervals and seams
remain, including HCNO at 1382 K, HOCN at 1368 K and HNCO at 1478 K.
The common admitted interval is 300–3000 K, limited above by CH3O.

The elementary owner admits at most 512 channels.
`ElementaryReactionV1::irreversible()` explicitly disables reverse progress
and reverse derivatives; reversible construction remains the default.
All 16 originally irreversible GRI channels retain that directionality.
Their `R*q*affinity` contribution is signed, exposing physical incompatibility
rather than changing source rates or manufacturing reverse chemistry. The
closure operation version is 3. This does not relax a consumer's finite-step
entropy gate or establish detailed balance or nonnegative entropy production
for the complete mechanism.

Five source tests compare all 325 net rates and all 53 Cp/h/s triples at
700 and 1800 K with independent Cantera 3.2.0 output, check source identities,
conservation, continuity/domain, derivatives across irreversible, duplicate,
collider and falloff channels, and cancellation/reuse/pure-inert states.
Two additional elementary regressions check forward-only progress,
product-only states, signed entropy and collider derivatives. All 27 focused
native kinetics tests pass. This establishes source and numerical conformance,
not experimentally validated ignition, emissions or fire predictions.
The source supplies methane/natural-gas chemistry with nitrogen species;
it does not supply a soot model, PMMA/wood volatile mechanism, transport
closure, sprays, real-gas phases or a time integrator.


## Near-atmospheric water and freshwater-ice conductivity

`water::conductivity` supplies the original IAPWS SR6-08(2011) equation (8)
liquid conductivity and the explicitly empirical Yen (1981), CRREL Report
81-10 equation (33) ice conductivity. Both return typed W/(m K) values and
analytic W/(m K²) temperature derivatives without changing caloric data.
Liquid admission is 273.15–383.15 K and the original IF97 saturation pressure
through 0.3 MPa; the temperature-only relation does not select stability
against ice. The selected ice window is 243.15–273.16 K for ordinary
polycrystalline freshwater ice near atmospheric pressure. The ice source
supplies no pressure correction, crystal-axis tensor, porosity/salinity law
or quantified experimental uncertainty bound. Original coefficients are not
retuned to the report's rounded discussion of conductivity near melting.

Source domains refuse rather than clip or extrapolate. These allocation-free
safe Rust correlations use deterministic repository math and add neither heat
capacity nor latent heat. G0/G1 tests reproduce the admitted original liquid
verification values, check analytic derivatives and the unchanged ice
regression, and exercise exact domain and pressure boundaries. Numerical
agreement is implementation evidence, not validation of an arbitrary frozen
specimen or an effective mixed-phase microstructure.


## Peng--Robinson real-gas mechanical EOS

### Isothermal phase equilibrium and tangent-plane searches

`eos::peng_robinson::stability` searches the original component fugacities
from liquid-like and vapor-like Wilson seeds, using damped logarithmic
successive substitution and the minimum-Gibbs root at each trial composition.
A negative dimensionless tangent-plane distance is a concrete instability
witness for the original feed phase. Two nonnegative converged stationary
searches are not a global stability certificate. Iteration exhaustion refuses;
neither a failed search nor a missing Wilson bracket declares a stable phase.
Exact feed absences stay absent, and underflow of present species refuses.

The existing `two_phase_flash_pt` permits incipient single-phase iterates and
can initialize from a tangent-plane witness if Wilson iteration fails. It
publishes only a nontrivial split with Rachford--Rice material balance,
component fugacity equality and no Gibbs increase relative to the best
homogeneous feed root. The Gibbs change is returned explicitly; it does not
prove stability against a third phase. Pure-fluid or unresolved azeotropic
coexistence needs a supplied phase amount or another extensive constraint and
remains outside this two-phase amount solver. Both APIs have cancellation-aware variants, bounded
evaluation counts, immutable inputs and no partial publication. Phase trial
models reuse admitted component/interaction storage instead of cloning and
revalidating the quadratic interaction matrix every iteration.

G0/G1/G3/G4 regressions cover an independently solved methane/ethane split
outside Wilson's dew bracket, negative tangent-plane witnesses, a homogeneous
gas, component permutation and absence, Gibbs/material/fugacity closure,
iteration limits and cancellation. The physics is the supplied cubic EOS and
mixing rule; this is not experimental validation of its parameters. Reference:
Michelsen (1982), <https://doi.org/10.1016/0378-3812(82)85001-2>.

### Fixed-composition mechanical evaluations

`eos::peng_robinson` adds a bounded fixed-composition Peng--Robinson mechanical
closure. Each component carries caller-supplied critical temperature, critical
pressure, Pitzer acentric factor and molar mass; an optional explicit symmetric
binary-interaction matrix supplies classical quadratic mixing. No critical
property, interaction coefficient or phase identity is inferred from a species
name. Up to 64 components are retained.

A PT evaluation keeps every real cubic root satisfying the physical covolume
constraint `Z>B`, with molar volume, density, residual Gibbs energy and one
log fugacity coefficient per component. Callers may request the smallest root,
largest root, or the minimum residual-Gibbs root at the SAME fixed composition.
The last choice is a mechanical branch selection, not a multiphase flash or
global tangent-plane stability proof. No root clipping, van-der-Waals loop
smoothing, phase lever rule or critical-region regularization is added.

Focused tests pin an independently evaluated supercritical methane root and
density, retain all three subcritical methane roots, recover the ideal-gas
dilute limit, verify binary-component permutation/fugacity invariance, and
exercise malformed source/composition/state refusals. This creates the reusable
nonideal mechanical/fugacity primitive needed by later dense-gas flow and flash
owners; it does not by itself change the existing ideal-gas solvers.

### Conservative density coordinates

`PengRobinsonMixtureV1::evaluate_density_temperature` evaluates the same mixing
law directly at `v=M/rho`, retaining the caller's finite positive mass density.
It shares mixture parameters, pressure derivatives, residual internal energy,
residual Cv and fugacity evaluation with the PT owner, but does not solve a
cubic or choose a different density. Covolume violations, nonpositive pressure,
unrepresentable arithmetic and mechanical instability `(dp/dv)_T >= 0` refuse.
The returned homogeneous root is not an inventory of all PT roots. Mechanical
stability does not establish equilibrium against phase separation.

`evaluate_density_caloric` supplies the same homogeneous pressure, residual
U/H/entropy/Cv and pressure derivatives in a fixed-size result. Its shared
mixing reduction uses four stack arrays bounded by the admitted 64 components;
no successful scalar evaluation allocates or calculates component fugacities.
The original pair arithmetic, including `sqrt(a_i*a_j)` on either side of an
alpha-base sign change, is retained. Full evaluations record the same pair
values for their fugacity owner. Both paths share the residual and derivative
arithmetic and the pressure/covolume/mechanical-stability checks. Scalar trials
produce no Gibbs or phase-equilibrium result; a consumer needing fugacity or
phase admission must still request the full state. Focused tests require bitwise
agreement in dilute, dense, metastable, alpha-sign and 64-component cases, and
matching physical refusals. The original fixed operation-tree version remains 4.

G1/G3 tests compare PT and density-coordinate states, the independent pressure
equation, residual-energy temperature derivatives and the dilute limit. The
source adapter supplies ideal-reference Cv before complete caloric and acoustic
stability can be admitted. This coordinate extension enables conservative
energy reconstruction; it does not implement a spatial compressible solver.
The cubic residual-thermodynamics convention is described in
<https://pages.nist.gov/teqp-docs/en/latest/models/cubics.html>.

Residual entropy exposes its reference coordinates explicitly. The density
form is relative to ideal gas at the same temperature and molar volume,
`sR = R ln(1-b/v) + da/dT I(v)`. The pressure form adds `R ln Z`, giving
`R ln(Z-B) + da/dT I(v)` relative to ideal gas at the same temperature and
pressure. Both use the original PR attractive integral and mixing derivatives;
they introduce no new caloric fit. Their finite differences recover `CvR/T`
and the original pressure/temperature Maxwell derivative, while the pressure
form retains the Gibbs identity with the existing residual enthalpy. These
primitives support source-bound entropy inversion in `fs-material` without
coupling this mechanical owner to transport cards or mutable parcel state.

### Stable finite homogeneous calorics

`prepare_density_caloric_change` borrows the original mixture, admits its initial
T/rho state and fixes the final density. Its `caloric::PengRobinsonCaloricReferenceV1`
prepares density differences once and retains at most 64 component records.
Each `evaluate` returns the original scalar endpoint and finite residual U/H
and density-reference entropy changes in one original mixing traversal, with
no allocation, second quadratic pair pass or independently fitted EOS. Pair
differences factor Delta sqrt(T) on their actual absolute-alpha branch; a
changed alpha-product sign retains the original endpoint pair values and
derivative jump. Density integrals, excluded-volume work and attraction-volume
work use factored differences and symmetric product identities. The absolute
state API and its fixed operation-tree version 4 are unchanged.

Ideal-reference source increments remain separate. For entropy, an ideal
source evaluated at unchanged pressure receives the explicit correction
`-R [ln(T1/T0)+ln(rho1/rho0)]` and the PR density departure. This avoids
cancellation between separately rounded real pressure and `R ln Z` terms.
Preparation/trials retain original density, covolume, pressure and local
mechanical-stability refusals; complete Cv and material/transport domains
remain consumer obligations. These are fixed-composition homogeneous endpoint
differences, not a heat/work path, mixing law or phase-equilibrium latent heat.
Focused references use independent 85-digit dimensional PR calculations for
thermal/density changes thousands of coordinate ULPs above resolution, together
with exact endpoint reversal, alpha-sign, full-component-bound and refusal cases.

## Pure-fluid coexistence and conserved phase amount

`eos::peng_robinson::saturation` closes pure-fluid liquid/vapor coexistence by
equal logarithmic fugacities from the original EOS. Its bounded solve uses
`d(ln(phi_l)-ln(phi_v))/d(ln p)=Z_l-Z_v` with bisection safeguarding and accepts
pressure brackets whose endpoints are on distinct single-phase branches. The
published outer roots must both have negative isothermal pressure-volume
derivatives, a resolved volume gap, and closed fugacity and pressure-correction
tests. Numerical critical degeneracy refuses; these tests are not interval
enclosures or evidence of stability against a solid phase.

The maximum work is a whole-call EOS evaluation cap, including both endpoints.
Cancellation is checked before every EOS evaluation and result publication.
PT leaves quality unset. An explicit mean molar volume fixes vapor amount by
the lever rule; density uses conserved mass over additive phase volume. Latent
enthalpy/internal energy, Clapeyron slope and amount-weighted residual h/u/s
come from the original departure functions. Their matching ideal reference is
added once by the caller; there is no extra latent source or ideal-mixing
entropy between identical phases. Two-phase transport, interfaces, metastable
kinetics and chemical reactions are separate closures. The explicit equilibrium
pressure/entropy response below adds its own acoustic assumptions.

Six focused regressions compare independently computed Maxwell-area methane
and ethane saturation states, differentiated Clapeyron response and latent
energy, conservative volume/quality closure, single-phase bracket endpoints,
source-parameter changes, critical-domain refusal, cancellation and work caps.
The original constitutive law is Peng and Robinson (1976),
<https://doi.org/10.1021/i160057a011>; these verify that law's implementation,
not experimental accuracy across all fluids or critical conditions.

## Pure-fluid equilibrium isochoric caloric response

The equilibrium isochoric residual-energy tangent differentiates the same
moving saturation roots and volume lever rule using original pressure
derivatives and the Clapeyron slope. The material owner adds its ideal molar Cv
once. This is not a frozen-quality average; at either saturated endpoint the
tangent is one-sided from inside the dome. Interior fixed-volume finite
differences independently check the derivative at several qualities.

## Maxwellian Coulomb momentum and thermal moments

`coulomb::CoulombModelV1` retains original species thermodynamics, masses,
charges, populations, volume, caller-selected Debye screening and explicit
classical weak-coupling limits. One reciprocal coefficient supplies each
electron-ion or ion-ion pair's momentum transfer and its matching thermal
moment. Original source support, screening occupation, Coulomb logarithm,
coupling, phase-space density, quantum length, relativistic temperature and
actual relative drift are checked without a logarithm floor or a quantum fit.
Pair drift heating separates positive friction shares from the signed
unequal-temperature correction; their sum is relative-velocity dissipation.

`evaluate_at_state` explicitly refreshes the same law at a new source-bound
state with the exact original thermodynamic bank. Newly populated and vanished
charged species, density, volume and separate temperatures enter the same
admission and pair evaluator. Bounded stack carrier metadata makes successful
refresh allocation-free. This does not relax `is_bound_to_populations`:
material transport still needs a new receipt for changed populations or volume.
L1 returns local moments and does not own time evolution or a plasma solver.
Focused tests compare refreshed summaries with separately constructed models,
including both screening choices, source-domain rejection and cancellation.

## Homogeneous-equilibrium pressure and entropy response

`PengRobinsonSaturatedMixtureV1::equilibrium_response` differentiates the
original pure-fluid coexistence under instantaneous thermal and phase
equilibration. In the resolved interior of the dome, `p = p_sat(T)`,
`ds = (Cv_eq/T) dT + p_sat' dv`, and `(dp/dv)_s = -T (p_sat')^2/Cv_eq`.
The supplied original ideal molar Cv is added exactly once to the existing
moving-phase residual tangent. Both complete single-phase capacities must be
positive; no mixed constant-pressure heat capacity is invented.

Common phase velocity and negligible interfacial storage give
`K_s = -v (dp/dv)_s` and `c^2 = K_s/rho`. Exact or numerically unresolved
phase endpoints refuse this bidirectional small-signal response because the
interior and homogeneous branches have distinct one-sided derivatives.
Original saturation admissibility, including critical coalescence refusal,
remains authoritative. This equilibrium limit supplies no phase relaxation
times, valid frequency band, slip, nucleation, bubble resonance or attenuation.

The frozen and equilibrated acoustic closures are physically distinct; see
Lund and Flåtten (2010), *Equilibrium conditions and sound velocities in
two-phase flows*,
<https://www.sintef.no/globalassets/project/co2-dynamics/publications/lund_an10.pdf>.
Core regressions check the Maxwell relation, pressure-volume/energy closure,
original ideal-Cv dependence and endpoint/reference refusal. Material-level
regressions independently perturb full entropy-preserving equilibrium states.

### Screened homogeneous and separated PT admission

`eos::peng_robinson::phase_screen` composes the existing bounded tangent-plane searches and original two-phase flash. Both Wilson-seeded stationary searches must converge. A negative tangent-plane witness must resolve an admitted fugacity/material/Gibbs-checked split; failure does not fall back to the feed phase. Without a negative witness, the retained minimum-Gibbs feed root must have negative original `dp/dv|T`. Distinct mechanically stable feed roots within the declared Gibbs tolerance refuse because pressure and temperature do not determine their phase amount, including pure-fluid or azeotropic coexistence.

Each result retains both search receipts, all original feed roots, the exact numerical options, and the selected homogeneous root or separated phases. Two converged nonnegative stationary searches are finite screening evidence and do not certify a global minimum over all compositions or phases. No homogeneous vapor amount or liquid/vapor label is inferred.

The stability receipt reports its actual cubic EOS count (one feed plus both trial counts). A phase-selection receipt reports that observed screening work plus the original flash's conservative work ceiling when a split was attempted. The whole-call ceiling is `1 + 2*stability.maximum_iterations + 6*flash.maximum_iterations + 2`; it is not an observed EOS iteration count. Original per-iteration cancellation and a final publication checkpoint are retained.

## Stable Peng–Robinson pressure root inventory

Peng–Robinson evaluator operation version 3 retains the version 2 correction
in the single-real-root Cardano branch. For the depressed polynomial `y^3+p*y+q=0`,
it takes the larger-magnitude cube-root argument and recovers the other factor
from `u*v=-p/3`. This avoids subtracting nearly equal `q/2` and the square root
of the discriminant before taking a cube root. The source then admits at most
two finite Newton corrections, each required to decrease the residual of the
original cubic.

The near-degenerate depressed-cubic branch now polishes the best-separated
trigonometric candidate and recovers the other two roots through the original
Vieta sum/product equations. Its scaled quadratic discriminant distinguishes
true small-pressure liquid roots from a cancellation-created phantom branch.
Quadratic division avoids subtracting nearly equal terms; root deduplication
uses the roots' own relative scale. Dense molar volumes therefore remain
resolved as pressure decreases, even when their compressibility factors are
tiny beside the gas root. The existing fixed work limit is retained.

Mechanical admission, phase/stability tolerances and fugacity/material/Gibbs
acceptance criteria are unchanged. This is an original EOS arithmetic
correction, not an alternative phase solver or relaxed equilibrium acceptance.
No global stability, near-critical root enclosure, experimental EOS accuracy or
general convergence guarantee follows from the correction.

Focused G1 tests retain the original binary state at
`T=217.1522724700456 K`, `p=2500001.5210795477 Pa`: the independent liquid
composition makes the small direct Cardano argument round to zero, although
its unique real root is well conditioned. Tests check the original polynomial
residual and independent compressibility factor, a small positive/negative
depressed-p family, and the complete screened binary flash against independent
simultaneous cubic/fugacity phase compositions and lever-rule amount. The
phase and stability tolerances stay at `1e-12`; iteration ceilings are retained.

Additional G1 regressions compare the 150 K methane liquid molar volume at
100 through 1e-8 Pa against independent 70-digit dimensional-polynomial roots,
check every returned original cubic residual on its own scale, and retain the
unstable middle branch without promoting it into a material. Above-critical
200–400 K methane has only its real gas root through the corresponding dilute
pressure range. Small, zero and repeated algebraic roots exercise deflation's
distinct conditioning cases. These tests establish numerical evaluation of the
declared EOS, not the physical persistence of metastable liquid in a vacuum.

### Pressure-consistent real-gas departures

Complete caloric/mechanical admission requires `Cp >= Cv > 0` and negative
`(dp/dv)_T`. The original identity `Cp-Cv = -T (dp/dT)_v²/(dp/dv)_T`
permits equality at zero thermal expansion; equality may also result when a
positive correction rounds away beside Cv. Both cases retain the original
isentropic sound-speed response and all mechanical/caloric stability checks.
G1 regressions compare a zero-expansion PR state with independent U, H and
pressure derivatives and admit a positive correction below Cv resolution.

Operation version 4 preserves the original cubic root inventory and evaluates
fugacity, internal energy, enthalpy, entropy and Cv from common immutable
residual geometry. Deterministic log1p arithmetic forms the covolume offsets
before rounding volume ratios to one. Near Z=1, the original cubic identity
`Z-1 = B-A*(Z-B)/(Z²+2BZ-B²)` retains the pressure-work departure even when
the published compressibility factor rounds to one. Dense roots retain their
original positive Z inside the logarithm. The shared attraction integral and
volume logarithms are evaluated once per root, so caloric queries introduce
no additional transcendental evaluations or allocation.

The same original PR parameters and source-temperature derivatives continue
to determine all departures. Phase selection, stability and flash tolerances
are unchanged. Focused regressions check the independent second-virial limits
of pure and partial-mixture fugacity and all caloric departures down to
1e-12 Pa, the Gibbs identity, and independently differentiated dilute-gas
energy/entropy Cv identities. These numerical limits do not extend the
physical validity of the source EOS or establish a continuum-flow regime.

## Homogeneous binary square-gradient constitutive owner

`solution::gradient::BinarySquareGradientV1` retains the exact original binary
regular-solution pair at one admitted temperature, an explicit positive site
concentration `c` [mol/m³], and nonnegative square-gradient coefficient
`kappa` [J/m]. It supplies the homogeneous exchange chemical potential,
convex-split potential and positive implicit curvature used by spatial
miscibility evolution. Local cells are never independently replaced by bulk
coexistence states. In the spinodal-capable `chi > 2` regime, zero `kappa`
refuses; the convex/ideal regime admits zero for the diffusion limit.

The spatial energy density is `c*g_mix(x,T) + kappa*|grad x|²/2`.
`kappa` is an independently qualified interfacial constitutive parameter;
it is not a liquid-liquid surface tension [J/m²]. Constant `c` assumes equal,
fixed partial molar volumes. Unequal species molar masses are permitted, but a
spatial consumer must declare the molar/volume reference frame.

Both composition coordinates remain independently represented. Pure-state
mixing energy is defined, while requesting a finite logarithmic chemical
potential at a pure endpoint refuses. Near the critical composition, the
exchange-potential series retains the small cubic term after cancellation;
direct free-energy increments retain the critical fourth-order contribution,
trace changes and exact endpoint reversal without subtracting two large
absolute energies. The split treats logarithmic entropy, negative-chi
interactions and gradient energy implicitly; positive-chi concavity uses the
previous composition. The owner supplies constitutive quantities only, with no
time integrator or claim of nucleation kinetics.

Five focused native tests cover original-energy derivatives, reversible and
trace increments, critical higher-order terms, the convex-split energy bound
for both signs of chi, and endpoint/regularization admission.

The same owner supplies direct `mixing_enthalpy_change_j_per_mol` and
`mixing_entropy_change_j_per_mol_k` for spatial caloric consumers.
For `chi=A+B/T`, these use the original `R*B*delta(x_A*x_B)` and original
ideal/excess entropy with coefficient A. The enthalpy product increment is
formed before multiplying by B and R; a finite increment need not form an
overflowing reference enthalpy. Entropy shares the existing reduced-energy
kernel, retaining trace changes and critical fourth-order terms. Exact endpoint
reversal and the existing free-energy evaluation order are preserved. These
are homogeneous composition increments only: pure Cp/CpT, gradient energy,
heat balance and temporal integration remain with their corresponding owners.

Four additional native regressions cover direct H-minus-T-S closure, pure and
one-ULP trace changes, quadratic enthalpy/critical quartic entropy, and finite
increments beside overflowing reference energies.

## Planar interface tension from original binary coexistence

`solution::gradient::interface` computes the equilibrium planar liquid-liquid
surface tension and a characteristic interface width from the SAME
`BinarySquareGradientV1`. For the symmetric regular solution, the common
tangent between original binodals is horizontal. With
`Delta_g = g_mix(x)-g_mix(x_alpha)`, the model evaluates
`sigma = integral sqrt(2*kappa*c*Delta_g) dx` between the two binodals,
consistent with the existing `kappa*|grad x|^2/2` energy convention.
The width is the composition gap divided by the central profile slope; the
exact equilibrium profile has infinite tails, not a finite total thickness.

The calculation retains the original model, coexistence receipt, requested
controls, quadrature estimate and sample count. Near-critical conditioning
uses direct cancellation-resistant free-energy increments and normalization
by the center barrier. Adaptive Simpson quadrature has one global
absolute-plus-relative accuracy budget, bounded work (at most 65,537 samples)
and cancellation. Its reported error is a truncation estimate, not a rigorous
certificate; floating-point roundoff and the separately reported binodal
uncertainty are not included. Negative or unresolved interior excess energy
refuses rather than being clipped. Miscible and critical states return None.

This is a model prediction from an explicitly qualified kappa, not a measured
material surface-tension claim, wall contact angle or automatic curved-interface
force law. Five focused native tests cover independent high-precision regular-
solution quadrature, the near-critical sigma~(chi-2)^(3/2) and width~(chi-2)^(-1/2)
limits with their prefactors, kappa/concentration/component-order scaling,
physical absence/refusal cases and bounded work/cancellation.

## Original IF97 isobaric density curvature

`water::if97::If97State::isobaric_density_jet` supplies density and its first/second temperature derivatives at the original fixed pressure. It retains the existing Region 1/2/5 state value and first Gibbs derivative, evaluating only the additional mixed third Gibbs partial from the original released coefficients on demand. Fixed-size power tables and the native compensated accumulator allocate nothing; ordinary state queries retain their previous operation tree and cost. The exact relation is rho_TT = rho * (2 * (v_T/v)^2 - v_TT/v).

Focused actual-source controls cover fourth-order differences of the original scalar EOS and retained first derivatives at liquid, vapor and high-temperature steam states, the dilute ideal-gas limit, water's density-maximum sign change and the separate saturation-branch endpoints. The jet retains the existing source domain and uncertainty. It is a single-phase branch derivative, not a two-phase mixture tangent, automatic regional continuation, or a pressure-derivative API.

`If97State::density_jet` additionally supplies the complete pressure/temperature density Hessian, including the mixed partial. The value and first responses reuse the same admitted state, and its temperature derivatives agree exactly with `isobaric_density_jet`. The extra original Gibbs partials share fixed-size power tables and compensated sums; the isobaric path omits pressure work. Vapor pressure curvature differentiates the residual of `pi*gamma_pi` directly, retaining the finite virial limit without cancellation of ideal-gas terms. The original scalar EOS operation tree, domain and source uncertainty remain unchanged.

Focused controls compare the pressure and both mixed-derivative directions with fourth-order differences of the original state across liquid, vapor and high-temperature steam conditions, preserve the separate saturation tangents, and recover the independently evaluated dilute virial curvature down to 1e-140 Pa. These are single-phase local derivatives; the API supplies no phase-mixture tangent or regional crossing.

`If97DenseState::density_jet` extends that same typed response to Region 3 by implicit differentiation of the original Helmholtz pressure law. Three compensated coefficient sums evaluate its pressure Hessian directly; first responses retain the admitted state's stiffness and thermal pressure slope. This allocation-free query requires no pressure inversion and leaves ordinary state cost unchanged. Its isobaric projection uses the same derivatives. `If97FluidState` routes both jet queries to the retained Region 1, 2, 3 or 5 source, preserving achieved pressure and solve diagnostics. Saturation endpoints keep their single-phase tangent and near-critical amplification remains visible; no finite mixture derivative, smoothing, extrapolation or reduced source uncertainty is implied. Focused regressions compare native pressure and density derivatives in both mixed directions with independent fourth-order differences, and cover dense coexistence, critical amplification and exact regional delegation.

## Source-bound pyrrhotite oxidant and acid balances

Three matched NBS/1989 packs add 15 point claims for crystalline alpha FeS,
Fe7S8 and orthorhombic sulfur at 298.15 K and 1 bar. Four conserved
`MineralReaction` endpoints separate finite dissolved-oxygen oxidation to
sulfate from ferric attack to elemental sulfur. Formula-unit and mineral-phase
identity remain explicit: Fe7S8 is not silently rescaled to Fe0.875S, and the
sulfur reference is per mole of S atoms. The four `mineral_processing` modes
compute source-driven mass, standard reaction heat and direct Gibbs sums.
Compiler-backed regressions also exercise both pyrrhotite phases competing for
one finite oxygen stock and refuse an overdraw before publishing inventories.

The source rows do not support continuous Fe1-xS interpolation or kinetics.
Endpoint selection and reacted extents are caller inputs; no oxygen gas supply,
passivation, ferric regeneration or downstream iron hydrolysis is inferred.
The FeS sulfate endpoint's zero direct proton term is not a declaration of zero
long-term acid generation. Source provenance, standard states, historical
molar masses and unstated uncertainty are retained in
`data/matdb/seed-v1/MINING_PYRRHOTITE.md` and its referenced packs.

## Clay-gangue acid demand and staged aluminum hydrolysis

`MineralReaction::{KaoliniteAcidLeach, AluminumHydrolysisGibbsite}` add conserved
six-proton clay dissolution and the reverse acid-generation endpoint for selected
crystalline gibbsite. Five matched NBS packs provide 23 point claims for kaolinite,
free Al3+, the formal H4SiO4 aqueous component, doubled-unit gibbsite and CuO.
The existing CuO acid-leach endpoint now has a compatible NBS source route.
Required selectors retain the source's H2SiO3(aq)+H2O(l) component convention and
its Al2O3·3H2O gibbsite formula unit; phase or formula-basis substitutions refuse.
All claims stay at 298.15 K and 100000 Pa with unstated source uncertainty.

`mineral_processing --clay-acid-inventory` uses the existing `ReactionNetwork`
for three explicit stages: competing CuO/clay leaching from one finite acid
stock, selected aluminum hydrolysis, then later CuO leaching. Every stage
consumes only its initial inventory; later or simultaneous acid/aluminum
production cannot fund an earlier consumption. Original source mass, atom,
charge and standard enthalpy accounts remain visible. Extents are caller inputs,
not predicted dissolution, precipitation kinetics, pH, stable phases or plant
doses. Separate CLI modes resolve source formation-H and direct formation-G.
`mining_clay_acid_sources` compiles all five new packs and checks all 23 facts,
identity/domain refusals, three H/G sums, finite reactants and staged acid return.
See `data/matdb/seed-v1/MINING_CLAY_ACID_DEMAND.md` for equations and commands.

### Partial molar calorics for selective real-gas streams

`PengRobinsonMixtureV1::evaluate_density_partial_calorics` supplies the original
homogeneous scalar state and each component's fixed-temperature, fixed-pressure
partial enthalpy departure and partial molar volume. The result borrows the
original mixture and retains its exact component order, including zero-fraction
entries. At most 64 records are stored; evaluation allocates nothing and shares
one original mixing traversal. Analytic composition and pressure derivatives
retain the original alpha-sign behavior. A converged analytic series evaluates
the attractive-integral derivative in its dilute cancellation regime; it is not
a replacement caloric law.

The material owner adds each original NASA species enthalpy to its partial
departure and converts with that species' molar mass. This already includes flow
work. Partial molar volumes remain signed. The Euler identities recover the
mixture enthalpy and volume without redistributing rounding residuals. Admission
retains the original positive-pressure, covolume and local mechanical-stability
checks; complete Cv, material support and phase equilibrium remain separate
original-source responsibilities. The API does not select phases, supply a
diffusion rate or certify a two-phase bulk partial enthalpy.

Focused checks cover independent 90-digit extensive derivatives at fixed T,p,
original fugacity/Gibbs–Helmholtz derivatives, Euler closure, pure and dilute
limits, dense/alpha-sign states, negative partial volumes, permutation, exact
absence, the 64-component bound and physical refusals. Existing scalar-state and
finite-caloric regressions remain unchanged.

## Source-bound copper sulfide precipitation and cyanide recovery

`acid_base::complexation::precipitation::CopperSulfidePrecipitation` couples
one selected pure Cu2S phase to the original `IdealLigandSystem`. Construction
binds the original CN ligand owner, Cu(I) and H2S/HS- families, conserved
`2 Cu+ + HS- -> Cu2S + H+` equation and identical shared formation-Gibbs
receipts. Each phase-partition trial reruns the original finite metal/ligand
and charge solver; free cyanide and pH are not frozen from an earlier state.

Three matched NBS profiles supply 13 original point facts for H2S(aq), HS-(aq)
and crystalline source-specific Cu2S(alpha), at 298.15 K and 100000 Pa.
Molecular H2S, the individual-ion convention and the full Cu2S formula unit
remain explicit. Missing heat capacities and all source uncertainties remain
unavailable/Unstated. Existing differently evaluated chalcocite is not silently
substituted. Source links, checksum, exact states and commands are retained in
`data/matdb/seed-v1/MINING_CYANIDE_RECOVERY.md`.

Cu and sulfide totals include solid at two Cu atoms and one S per formula
unit; CN stays in the aqueous/complex inventory. A directly retained dissolved
limiting-component coordinate preserves tiny Cu/S traces even when subtracting
the precipitate from the total would round to zero. Present solid must satisfy
the declared log-saturation tolerance; absent solid requires resolved
undersaturation. Finite Cu/S balances, original ligand/charge checks and bounded
work/cancellation remain mandatory. Returned phase brackets and inner spreads
describe numerical resolution, not experimental uncertainty or certificates.

The complete final pH bracket must lie in [2,10], an explicit selected-model
approximation domain. The assemblage omits S(2-), polysulfides, copper sulfide
complexes, CuCN/other solids, redox, gas escape and nonideal activities.
Fixed solvent reference does not constitute a full H/O solvent ledger.
No reaction rate, filtration efficiency, heat duty or calibrated recovery
follows from this selected equilibrium calculation.

FACT: direct source sums give phase Delta-r G=-198240 J/mol equation and
ln(K_precip)=79.9691150886. The net tricyanocuprate recovery endpoint gives
Delta-r G=-187680 J/mol equation. Its heat refuses because the selected
copper-cyanide sources lack formation enthalpies. All 29 exact-stoichiometry,
eight coupled-phase and five actual-source compiler tests pass in the original
source closure. Independent tests cover both limiting components, absent/empty
phases, incompatible sources, scaling, finite supply, cancellation and retained
copper traces near 1e-25 mol.

The `copper_cyanide_recovery` consumer completes feed speciation, declared
NaHS/HCl additions, explicit primary-solid removal and NaOH addition before
publishing a result. Its authored example removes 0.000999290419 mol Cu2S;
free CN rises from 0.00277691 to 0.00732330 mol after re-alkalization while
retaining residual copper binding and both dissolved Cu/S stocks. Secondary
solid, if present, stays in the final inventory. Independent sequence balances
include the removed primary solid. Unsupported final pH refuses with empty
successful stdout. These constructed checks do not establish plant operation.

## Finite salt assemblages with shared calcium

`acid_base::precipitation::salt::assemblage::CommonCationSaltSystem` composes
one through eight original `SaltPrecipitation` owners with one common cation,
distinct acid families and identical complete aqueous/source definitions.
Every coordinate reuses the original single-salt and acid/base solvers while
holding other solid stocks fixed. Its available cation includes dissolved
stock and only its own current solid. No coordinate can spend calcium already
held by another phase.

The final acceptance condition recomputes every phase's dissolution saturation
residual at one common aqueous state and pH. Present solids require resolved
saturation; absent solids require resolved undersaturation. Empty original
components retain undefined saturation diagnostics. Inner pH resolution enters
the phase decision. Independently summed cation and complete aqueous-family
accounts include all stoichiometric solid holdings and must close within the
declared allowances. Dissolved stocks are carried directly, retaining positive
traces even when a solid amount rounds to its analytical total.

Canonical phase order, at most 256 sweeps, at most 65,536 total aqueous solves,
the original inner solver bounds and cancellation govern the entire operation.
Remaining global work is reserved before each coordinate call. Exhausted work,
unresolvable complementarity or numerical range loss refuses without a partial
result or mutation of caller stocks. The immutable result retains original
phase/source owners, controls, all component accounts, shared aqueous state,
final phase diagnostics and actual work counts. It is a nominal numerical
result under a selected ideal model, not an interval or source-uncertainty
certificate. Shared-family/polymorph competition is explicitly unsupported.

The matched `calcite-solid-nbs1989-thermochemistry` pack adds five original
reference facts for pure crystalline calcite, with explicit phase, NBS reference,
298.15 K and 100000 Pa selectors and unstated source uncertainty. The independent
USGS calcite evaluation remains separate. `MineralReaction::CalcitePrecipitation`
adds the conserved Ca2+ + CO3(2-) -> CaCO3 equation. The original direct-Gibbs
column is retained independently of H-minus-T-times-S; source discrepancies are
not silently reconciled.

`limestone_sulfate_equilibrium` loads eleven authentic compiled source packs
and couples finite calcite/gypsum partition to the original carbonic-component,
sulfate and water speciation. Both selected phases consume one original calcium
inventory. In the authored default cases, adding 0.005 mol HCl changes pH from
6.6436764 to 6.5367336, calcite from 0.0151054091 to 0.0121458109 mol and gypsum
from 0.0324344778 to 0.0349775905 mol. Free calcium also increases. Every requested
case must succeed before output. These values demonstrate the selected model;
they do not calibrate a plant acid dose.

Seven independent assemblage tests and fourteen existing salt tests pass,
including a closed-form two-phase state, bitwise one-phase reduction, absent
and empty phases, scale/order invariance, trace preservation, incompatible
sources, budgets and cancellation. The original reference-reaction conservation
test and five actual-compiler/source integration tests also pass. The latter
construct an independent pH-8 state from literal formation energies and
protonation polynomials, and check source facts, complete accounts, final
saturation and acid addition. Negative calcium or invalid activity basis
refuses with empty successful stdout.

The model has fixed solvent activity and ideal aqueous activities. It does not
close a finite H/O solvent ledger or include calcium ion pairing, other phases,
CO2 exchange, kinetics, passivation, thermal evolution or specimen qualification.
Sources, original values, model equations and normal workspace commands are in
`data/matdb/seed-v1/MINING_LIMESTONE_SULFATE.md`.

## Atomic photon-ionization source

`kinetics::atomic_ionization::photoionization::AtomicPhotoionizationV1` binds
a single-photon event cross section to the original validated ground-state
atom/ion/electron channel. The first original table sample declares its physical
onset, which may have nonzero cross section. The effective onset cannot lie
below the original chemical energy gap; the original samples and provenance
remain unchanged. No upper-energy extrapolation or missing-source zero is used.

Discrete photon groups carry photon energy in joules and bandwidth-integrated
number flux in m^-2 s^-1. An allocation-free, bounded evaluation returns the
event rate and absorbed, chemical-storage and photoelectron powers per neutral.
Chemical storage uses the original source gap once. Local angular illumination
uses angle-integrated number intensity, not net vector flux. Compensated sums
retain input order and cancellation is checked every 32 groups and before return.

This L1 law performs no population update, attenuation, thermalization or
radiative-recombination evolution. A finite-volume consumer must conserve
photons and neutral inventory, and an electron-temperature consumer must account
for the changing electron population. The source assumes one ground-state
ionization event per absorbed photon and stationary heavy targets; it does not
infer excited-state, multielectron, recoil, secondary-ionization or nonthermal
photoelectron behavior. Focused tests specify edge/support refusal, independent
event/energy moments and in-loop cancellation.

### Published H I photoionization fit

`AtomicPhotoionizationV1::verner_hydrogen` supplies the Verner, Ferland,
Korista and Yakovlev (1996, ApJ 465, 487) ground-state H I analytic fit using
the published H I parameters and megabarn-to-SI conversion. Its published
13.60 eV onset and 50000 eV upper support remain distinct from the original
thermochemical gap used in every energy account. Structural hydrogen,
ground-state weights, mass and gap compatibility are checked before admission;
species display names do not select the law.

The fit evaluates without allocating or tabulating an invented interpolation
grid. Original caller tables remain available, and `points()` distinguishes
actual interpolation samples from the analytic source. This single bound-free
event fit does not supply total opacity, excited levels, secondary electron
cascades or relativistic/thermal validity for a gas consumer. Published source:
<https://www.pa.uky.edu/~verner/photo.html> and
<https://arxiv.org/abs/astro-ph/9601009>.

### Published He I and He II ground-state photon sources

`AtomicPhotoionizationV1::verner_helium` and `verner_helium_ion` evaluate the
literal He I and He II rows of Verner et al. (1996), Table 1, with the complete
shifted Eq.(1). The neutral fit retains its nonzero yw, y0 and y1 parameters;
it is not obtained by hydrogenic scaling. Source admission checks the He
elemental row, adjacent charges, original level statistics, He-4 masses and
the compatibility of the thermochemical gap with the published rounded onset.
Arbitrary species display names are allowed. Cross-section evaluation is
allocation-free and uses the existing deterministic elementary functions.

He II retains the published 54.42..=50000 eV domain. For He I the equation's
published 24.59 eV onset is retained, but this ground-state-only owner limits
the upper energy to 65.39 eV. The first He+ excited level is 329179.2940 cm^-1
above its ground state (about 40.81303 eV), so an incident photon can leave
excited He+ above about 65.40042 eV. The total fit cannot then assign the upper
ion's excitation and photoelectron heat to this single ground-state channel.
The conservative 65.39 eV limit remains below the opening across the admitted
0.005 eV source-gap rounding interval. It is an explicit restriction on this
consumer, not a claim that the published fit ends there. Source levels:
<https://physics.nist.gov/PhysRefData/Handbook/Tables/heliumtable6.htm> and
<https://physics.nist.gov/cgi-bin/ASD/ie.pl?spectra=He&units=1>.

All three selected fit laws decrease above their admitted onsets. The He I
derivative proof is retained next to `maximum_cross_section_in_domain`; it
uses the literal shifted-fit parameters and a strictly positive quadratic,
not a hydrogen-only assumption or sampled maximum. This supplies the exact
in-domain maximum for the existing Milne envelope admission. A global
envelope beyond support still requires the caller's separate assertion.
Low-temperature He I reciprocal capture can resolve its finite supported
tail; hotter states refuse when omitted contributions cannot be bounded.
The fit never supplies missing channel branches or a total recombination rate.

G0/G1 regressions cover independently evaluated cross-section anchors,
structural refusals, physical onset and upper support, monotonicity and Milne
envelope admission, original-gap photoelectron heat, and independently
integrated reciprocal event/kinetic-energy moments. These are numerical source
and reciprocity checks, not experimental validation of a plasma trajectory.

## Reciprocal spontaneous radiative recombination

`kinetics::atomic_ionization::radiative_recombination::AtomicRadiativeRecombinationV1`
retains one original photon-ionization channel and inverts its cross section
through the Milne relation. Electron spin and photon polarization cancel in
this statistical-weight convention. A shared Maxwellian quadrature returns
the event coefficient alpha [m^3/s], captured-electron kinetic-energy
coefficient L [J m^3/s], and emitted-photon coefficient chi/NA*alpha+L.
These multiply actual electron and ion number densities; molar concentration
consumers must supply the corresponding Avogadro factors.

The caller supplies a global nonnegative cross-section upper envelope and
its separate provenance. Construction checks all known source values; beyond
the source endpoint the bound remains an explicit caller assertion. Both
event and kinetic-energy tails must resolve within the requested allowance.
Unknown tails are not added as guessed rates or silently set to zero.
Quadrature estimates describe integration of the declared source law, not
uncertainty in that law. Source knots, thermal scales and bounded cancellation
are retained; original electron caloric support and a declared nonrelativistic
thermal-energy admission policy are enforced.

The original photon onset remains distinct from the chemical gap. In
particular, the rounded Verner H I onset retains the declared small
zero-cross-section interval above the thermochemical threshold; numerical
quadrature does not correct this source approximation. The law is spontaneous,
source-resolved capture, not total Case A/B recombination or a Planck/Saha
radiation equilibrium model. Stimulated recombination, excited-state cascades,
photon trapping, dielectronic completion and recoil are not inferred.
References: Castor, Radiation Hydrodynamics (2003), Eq. 8.58; Mao, Kaastra &
Badnell (2017), A&A 599 A10, <https://arxiv.org/abs/1612.01017>.

### Original-source nonisothermal binary-liquid transport

`solution::caloric::transport` supplies a molar-frame face with equal fixed partial molar volumes, original-source Cp/interaction calorics, and an explicitly authored A-minus-B heat-of-transport contrast. The two-point chemical/enthalpy split preserves the original chemical-potential-over-temperature identity. Total energy carries source partial enthalpy once; reduced heat includes reciprocal Dufour transport. Conductivity is defined at zero diffusive species flux. Unequal pure heat capacities do not fabricate Soret flow when the authored heat contrast is zero. Exact absent support is refused rather than floored.

The face entropy form is nonnegative. A spatial owner must separately enforce kappa=0, convex admission over its complete temperature interval (including the possible interior interaction maximum), source D/k/Q domains, conservative finite transfers, and finite-step entropy acceptance. Changing-inventory caloric recovery stays homogeneous and uses original-source enthalpy/entropy increments without an equilibrium phase flash. Nine source-native tests cover the conjugate energy/entropy identity, reciprocal heat, variable and excess Cp, critical and trace limits, changing inventories, bounded caloric inversion, interior convexity and cancellation.


### Original-source nonisothermal multicomponent transport calorics

`solution::multicomponent::caloric::transport` supplies original homogeneous
enthalpy and entropy changes for independently changed positive component
amounts. Signed exchanges retain each pure-liquid datum, original Cp and Cp/T
integrals, and the original A+B/T pair potential. Small changes use direct
increments rather than subtraction of total cell energies. The bounded
temperature inverse prepares composition once outside its scalar iteration;
each trial uses only original source integrals and the exact positive Cp.
No individual-cell entropy gate is imposed on an open exchange transaction.

The finite-face source record retains every chemical effort, divided partial
enthalpy and arithmetic mean composition in canonical species order. Enthalpy
contrasts cancel the common species gauge before absolute-value rounding,
retaining finite mixing heat under a large common sensible datum. In the
zero-net-molar frame it obeys Delta(mu_i/T)=effort_i+h_star_i*Delta(1/T), up to
one common gauge. The pure enthalpy division uses the original Gibbs–Helmholtz
integrals, including variable Cp. Uniform composition has exactly zero
chemical effort at unequal temperature; Soret/Dufour data are not inferred.
The complete original positive support and abundant-carrier complement are
retained without a concentration floor or dependent-species projection.

The optional temperature-interval admission applies a sufficient global
composition-curvature lower bound to the original full tangent-space
eigensolve at both bracket endpoints. Its validity throughout the bracket
uses the original interaction's affine dependence on reciprocal temperature.
This bounded numerical test is sufficient, not necessary, and is not an
interval certificate or phase-equilibrium calculation. These constitutive
primitives own no mesh, transport coefficients, flow, or evolution. A spatial
owner must enforce its molar frame, fixed-volume assumptions, closed finite
species/energy balances and actual total entropy production.

Seven source-native regressions cover signed exchange with variable Cp,
large-datum and quadratic increments, independent affine-Cp divided
enthalpies and the chemical-potential identity, retained trace carriers,
whole-bracket/full-tangent admission, common-datum contrast invariance,
unsupported targets and cancellation.


## Adjacent atomic charge stages on original sources

`ElectronImpactRateV1` and `AtomicIonizationV1` admit the single-electron
channel `e + A(z) <-> A(z+1) + 2e` for a nonnegative integer source charge z.
The original exact element/charge conservation, physical electron mass,
reaction mass compatibility, positive source energy gap, ground-state
statistics, massive-target approximation and cross-section support remain
mandatory. Negative-ion detachment, skipped charge stages and overflowing
charge increments refuse. No cross section or level data are inferred from
the element name or charge number.

The lower/upper stage indices and lower charge are directly inspectable.
The existing `neutral_index` and statistical-weight `neutral` names denote
the lower stage for compatibility with the original neutral-channel callers;
they no longer imply that its actual source charge is zero. The source-bound
Saha factor and three-body reverse use the two original level weights and
masses at the actual electron temperature. Generic tabulated photoionization
and its Milne radiative reverse use the same admitted stages and chemical gap.
Each published Verner H I/He I/He II convenience source checks its own declared
element and charge stage; the neutral H law is not reused for ionic targets.

Three focused G0/G1 regressions compare higher-stage impact and Saha rates,
their common-temperature derivatives, source charge/nucleus balance and
generic photon/Milne energy moments against independent analytic formulas.
These are manufactured ground-state laws, not experimental validation or a
multilevel, continuum-lowered, dense-plasma or total-recombination model.


## Differential reciprocal recombination emission

`AtomicRadiativeRecombinationV1::spectrum` prepares the original Milne source
at an admitted Maxwellian electron temperature without allocation.
`RecombinationSpectrumV1::evaluate` returns the photon-energy differential
event coefficient, captured-electron kinetic-energy density and the original
photoionization cross section at that energy. It retains the original onset,
ground-state statistical weights, finite support and nonrelativistic temperature
admission. The density is angle-integrated, per joule of photon energy; neither
a frequency Jacobian nor an angular/escape factor is inferred.

`tail_bound` analytically integrates the retained global cross-section envelope
above a supplied energy, including beyond original source coverage. It returns
absolute event, captured-electron and emitted-photon energy bounds. This is
conditional on the original envelope assertion, not source extrapolation or an
interval/physical-uncertainty certificate. The outward arithmetic allowance
scales with log-term magnitudes; positive omitted bounds retain a least-positive
representable value on underflow. Differential positive underflow refuses.
A consumer may deliberately omit an insignificant far tail only by retaining
and checking its bound while continuing to evaluate any applicable absorption law.

This closure owns no populations, angular field, photon transport or evolution.
Focused tests independently evaluate constant-cross-section differential and
incomplete-gamma moments, recover shifted-onset totals through composite Simpson
refinement, and exercise original-domain, zero-law and underflow behavior.
Manufactured source agreement is numerical verification, not experimental
recombination validation.


### Repeatable reacting isobaric caloric storage

`composition::enthalpy::isobaric::EquilibriumIsobaricCaloricV1` prepares the original neutral ideal-gas Gibbs/HP owner at fixed pressure and conserved element ratios. Its closed-storage coordinate is `(H-H_cold)/feed_mass`; zero is the admitted lower-temperature equilibrium and is not an independently compatible flow-stream enthalpy datum. The complete original source bracket, hot equilibrium, and original bounded inverse at the exact advertised inclusive upper specific-enthalpy endpoint are admitted before construction succeeds. No tolerance is widened to accommodate coordinate rounding.

Every forward or inverse query restores one bounded original numerical anchor, including caloric source evaluations, so accepted constitutive values and publication do not depend on query order or failed trials. Original reacting Cp, density, composition, source receipts, work limits, and cancellation remain authoritative. Failed queries clear publication. Direct changes between two storage states use a fresh local original-source reference, retaining small enthalpy and entropy increments far from the cold datum; requested stored heat and actual bounded-solve enthalpy changes remain separately observable.

This native law supplies reacting calorics only. The material/thermal adapter owns global slope certification, transport admission, finite heat transfers, and the fixed-pressure enthalpy/pressure-work convention; no sampled caloric surrogate, frozen-Cp substitute, or second heat integrator is introduced. Five added regressions cover analytical isomer/dissociation calorics, mixing entropy, molecular-count response, canonical repeatability, small offset-invariant heat increments, bounds, work limits, and cancellation. They are authored but unexecuted after executor/resource failure; prior native runs predate this extension.


### Repeatable reacting isochoric caloric storage

`composition::volume::energy_path::isochoric::EquilibriumIsochoricCaloricV1` binds the original TV/UV Gibbs owner to fixed conserved element ratios and feed mass/volume. Its storage coordinate is `(U-U_cold)/feed_mass`, with zero at the admitted lower-temperature equilibrium. Pressure and molecular count respond to heat; the Newton tangent is the original reacting Cv, never frozen Cp-R. The complete continuous source bracket, hot TV state and exact advertised upper specific-energy endpoint pass the original bounded inverse before construction succeeds.

Every query restores one bounded original numerical/source anchor. Numerical history and failed trials do not change accepted constitutive values; refusals clear publication. Direct local U/S changes retain small increments far from the cold datum, while actual and requested energy changes and original volume/energy residuals remain separately observable. Fixed conserved density is distinguished from the native response density carrying its admitted residuals. This law supplies no time integration, transport model, finite-rate chemistry, pressure controller or flow-enthalpy datum. Native regressions cover reacting U/Cv/pressure/mixing entropy, dissociation and extensive scaling, repeatability, formation offsets, source bounds, work ceilings and cancellation.


## Pure-fluid isobaric sensible and latent calorics

`eos::peng_robinson::saturation::isobaric::PengRobinsonIsobaricCaloricV1`
binds an accepted pure-fluid coexistence point to the original PR component
and ideal-reference NASA source. Pressure remains the accepted saturation
pressure. Saturated-liquid enthalpy is zero; the inclusive latent interval
selects quality at constant temperature, with zero dT/dh and nonzero latent
volume response. Temperature alone at coexistence refuses an undetermined
phase amount. Sensible branches retain original stable roots and EOS Cp.

Continuous source admission, positive ideal Cv, nonnegative original PR alpha
curvature, and separated subcritical coexistence bound this fluid-only path.
Source-relative NASA and finite PR departure increments cancel formation and
entropy datums before rounding. Local changes retain source entropy separately
from the caloric latent integral Delta h/T and expose the original H-T Delta S
closure residual. No source data are fitted or normalized.

Every inverse starts from immutable endpoints with bounded PT work, simultaneous
energy/temperature convergence and final cancellation admission. Unsupported
sources, critical/multicomponent domains and unrepresentable nonzero phase
amounts refuse. Five exact-production-source native G1/G3/G4 regressions passed
for independent Maxwell latent data, sensible inversion and derivatives, datum
invariance and small heat, cross-phase H/S/reversal, and domain/work/cancellation
boundaries. This is a constitutive closure; heat integration, phase kinetics and
a global multicomponent flash certificate remain outside its claim.

### Direct caloric anchors on an immutable chemical isentrope

`composition::isentrope::EquilibriumCaloricAnchorV1` captures an accepted
current workspace's normalized species amounts and original source enthalpies.
`caloric_anchor` and `specific_enthalpy_change_from` are bounded, cancellable
operations on the existing source-integral owner. A direct difference rejects
changed species sources, elemental/charge inventory or feed normalization and
never rebases the path's birth entropy. Failed current samples cannot supply
an anchor or a direct difference. This resolves weak current energy changes
after a large historical change without subtracting rounded absolute or
birth-relative enthalpies. It is a caloric operation, not a flow integrator.

### Pure-fluid coexistence caloric increments at fixed density

`eos::peng_robinson::saturation::isochoric::PengRobinsonIsochoricCaloricReferenceV1`
borrows one original pure PR/NASA source, a resolved interior coexistence
state, a fixed density and a continuous subcritical NASA interval. Each query
retains original component, temperature, phase and volume admission, including
the original volume lever rule's 128-epsilon relative allowance.

Direct final-minus-reference internal energy and source entropy combine the
existing prepared homogeneous density-caloric increments with the symmetric,
compensated phase-amount product identity. Common ideal, formation and entropy
references enter once. Identical species add no mixing entropy, and rigid
volume adds no pressure work. The original moving-phase Cv, pressure, entropy
and acoustic response remains authoritative. Original nonnegative alpha
curvature and mechanically stable phase roots establish the residual-Cv sign;
a few sampled capacities are not used as a global bound.

Nonidentical queries report four bounded phase caloric algebra evaluations;
identical rest reports zero. Queries poll cancellation and publish no partial
result. This module introduces no saturation/UV solver, time integrator,
kinetics, interfacial storage, phase-departure continuation or critical
regularization. Source entropy retains the original accepted Gibbs defect.

The focused G1/G3/G4 source regressions cover original absolute U/S oracles,
moving-phase derivatives, weak heat under large reference datums, composition
and reversal, source/density/phase refusals and late cancellation. These new
regressions were authored but could not be compiled or executed while the
development executor was unavailable.


### Source-relative receiving UV equilibrium from an accepted SP state

The prepared `composition::volume::energy_path` owner accepts an opaque
`EquilibriumCaloricAnchorV1` through `sample_from_caloric_anchor`. Exact original
chemistry, elemental/charge inventory and feed normalization must agree. The
original TV/Gibbs closure and safeguarded UV root solve the requested extensive
`U-U_anchor`; per-species source integrals and amount-reference terms precede
summation. The shared internal-energy integral retains double-double amount,
formation-energy and RT arithmetic, matching the existing enthalpy treatment.
Neither the prepared UV birth reference nor the anchor's SP entropy reference
changes. The returned anchored balances are separate from the ordinary sample's
immutable UV-reference balances. Failure or cancellation invalidates current
publication. Direct SP entropy differences use the same source/mixing-entropy
owner and exact anchor binding.

This is a caloric/source integration seam, not a mass-flow law or a different
reaction solver. A receiving-volume consumer owns its mass scaling, finite
exported enthalpy, EOS reconciliation and total energy/entropy acceptance.
Focused frozen-limit, large-formation-datum and source/cancellation regressions
are in `fs-gas::equilibrium_flow::tests::nozzle_pair_anchor`.

Three focused native integration regressions are authored but uncompiled and
unexecuted while the development executor is disconnected.
