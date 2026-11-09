<!-- fs-material contract, part 4 of 4; a higher number is newer. Index: ../CONTRACT.md. Keep this file under 900 KB: the GitHub contents API returns empty content above 1 MiB, and a writer that sees an empty file overwrites it. -->

## Original strain domains in checked solid consumers

`SmallStrainLaw::admits_strain` retains each original law's nonallocating
strain admission. Native plane-stress reduction and P1/plane/HEX8
force, tangent and endpoint consumers enforce it before constitutive work.
Isotropic/orthotropic elastic laws retain their calibrated tensor Frobenius
limits; source tensor and frozen thermal laws retain their original Mandel
and total/elastic bounds. B-bar checks compatible and assumed strains.

Numerical allowances cannot extend a material domain. Retained plane-stress
strains remain numerical warm starts; the material admits the current
constitutive strains, not an inferred history. Raw constitutive and evidence
extrapolation semantics remain available. This adds no finite-rotation or
general material-history admission. Focused material and native solid
regressions cover calibration, transverse relaxation, heterogeneous points,
source bounds and B-bar projection.

## Native progressive axial ductility

`ductile::uniaxial::UniaxialDuctileJ2` shares immutable native J2,
initiation and receiving-band parameters. Owned nongeneric states retain their
exact model identity across frame clones and physical-interval preparation.
The explicitly one-dimensional model attenuates all axial stored energy and
stress in both signs; it does not inherit the three-dimensional intact
compressive-bulk split or a crack-closure law.

Its tangent differentiates native plasticity and moving rate-dependent onset.
Original native energy closes endpoint work into stored energy, plastic loss,
damage release and numerical loss, with separate accounts. The existing
`DuctileWorkEvolution` supplies exponential and exact finite-support profiles;
`Gp` keeps its degraded equivalent-plastic-work capacity meaning, not a new
total fracture-surface-energy claim. Last accepted incremental energy is
retained without another return map; unchanged-strain holds record zero work.
Negative and zero tangents remain physical outputs for the structural owner.

The six `ductile_uniaxial` regressions cover original native parity, independent
softening/work integrals and tangents, rate sensitivity, reversals and energy,
band scaling/terminal failure, and clone/owner/refusal boundaries. No stiffness
floor, erosion, heat partition, finite-strain necking or crack contact is added.

## Mixed multiplicative J2 plasticity

`plastic::finite::mixed::FrozenFiniteJ2Remainder` borrows the original
rate-independent `FiniteJ2` and accepted plastic metric. Linear, Voce and
tabulated hardening use the same native radial return, history, physical yield
heat and consistent logarithmic tangent. The bulk-free stress, tangent and
stored/incremental energy omit the original bulk contribution at assembly;
they do not subtract a large penalty or replace its modulus.

The generic `MixedFiniteStrainPotential` boundary requires native parameter/
history validation and an original pressure dual, without requiring an unrelated
`Hyperelastic` owner. Existing hyperelastic, viscous and thermal-reference
implementations preserve their prior arithmetic. The J2 adapter uses the
original log-squared volume law on its locally convex `log(J)<1` branch, or
explicit exact incompressibility without evaluating the unused bulk penalty.
Both actual deformation and pressure-reconstructed volume retain the original
Jacobian and elastic-strain domain checks. Spatial pressure equations impose
volume weakly; physical bulk storage is `U(J(p))`, not a second `U(det F)`.

`remainder_response`, `stored_remainder_energy` and
`validate_remainder_committed` expose candidate history, physical remainder
storage and independent accepted-state admission to native mixed solid owners.
Newton trials never mutate accepted history. This baseline adapter adds no
thermal, porous, ductile, spatial locking-convergence or experimental-validation claim.
Focused G0/G1/G3 cases in `tests/plastic_finite_mixed.rs` compare the full native
return and pressure condensation, independently difference every mixed tangent
block, continue reverse plastic loading, and test an unused overflowing bulk
penalty plus unchanged-history refusals.

## Mixed finite creep and overstress

`plastic::finite::mixed::MixedFiniteInelasticMaterial` extends the original
complete-history protocol to the same mixed pressure formulation.
`FrozenFiniteInelasticRemainder` borrows the original source owner and typed
accepted state. Perzyna, Norton, simultaneous Norton/J2, primary hardening and
independent resistance recovery reuse their original implicit returns and
physical interval. The original rate-independent API and source adapters remain
available. Bulk contributions are omitted directly; no kinetic law is replaced
by its inspection backbone or a different modulus.

Accepted Perzyna overstress and creep stress remain physical; combined J2 still
enforces equilibrium plastic yield. Recovering resistance stays distinct from
permanent plastic/creep accumulation and carries no invented storage or heat.
Flow potential, physical plastic/creep losses and pressure-constraint work remain
separate. Fixed-history energy/domain checks do not flow or recover again, and
exact incompressibility never evaluates unused penalty storage or stress.
Spatial consumers must use one common kinetic interval and publish complete
history only after equilibrium and work acceptance. Evolving temperature uses
the original thermal family's mixed extension below.

`tests/plastic_finite_mixed_kinetic.rs` covers full native/condensed parity for
all admitted kinetic families, independent value/stress derivatives, an unused
overflowing penalty, and analytical unloaded recovery with unchanged-history
refusals. These synthetic regressions make no material-calibration claim.

## Mixed finite thermal inelasticity

`plastic::finite::thermal::mixed` extends the original thermal point protocol
with a borrowed mixed potential and fixed-history mixed thermodynamics.
Native and source-bound affine J2, Perzyna, Norton, simultaneous creep/J2,
linear/Voce/tabulated hardening, primary creep and recovering-resistance
families preserve their complete original history and actual physical interval.
No new return map, thermal coefficient convention or caloric material is inferred.

These native thermal laws store energy per ORIGINAL volume, unlike a relaxed-
volume hyperelastic thermal reference. With `q=theta(T)^3`, physical pressure
`p` maps to native mechanical pressure `q*p`. The original-volume pressure dual
is `C(q*p)+p*(q-1)`, reconstructed volume is `q*Jm(q*p)`, and compliance is
`q^2*C''(q*p)`. Bulk storage remains the original `U(Jm)`, without a factor q.
Exact elastic incompressibility imposes `J=q` weakly and evaluates no unused
bulk penalty. The original actual and reconstructed strain/volume bounds remain.

Bulk-free thermal inspection and incremental work reuse the native owner with
bulk omitted at its original assembly sites. The mixed bulk account adds only
pressure-reconstructed storage. Its entropy retains
`s_bulk=3*a(T)*p*J_reconstructed-(K'(T)/K(T))*U_bulk`, including the nonzero
pressure/expansion entropy in exact incompressibility. Mixed thermodynamic
results deliberately do not relabel the displacement-only heat capacity or
temperature-stress derivative as a fixed-pressure derivative.

`checked_mixed_increment` reuses already frozen old/new material points;
`checked_mixed_increment_for_step` supplies the family convenience path.
Accepted history is checked at OLD temperature before a new-temperature return,
so genuine thermal softening can produce additional yielding. Recovering
resistance remains independent of permanent accumulated strain. Inspection
never flows or recovers and kinetic inspection cannot substitute for a step.

Physical inelastic heat, reversible caloric exchange, coefficient work,
`delta(T*s)`, and signed endpoint defects remain separate. Thermal bulk-volume
work uses `p_new*(q_new-q_old)*Jm_old`; the remaining actual endpoint-volume
work error belongs to the numerical defect. Neither that error nor flow
potential is deposited as heat. The point identity closes mechanical work into
mechanical internal-energy change, physical caloric transfer and signed defect.
Supplying temperatures does not establish a coupled enthalpy capacity, heat
transport solution, temporal accuracy, or material calibration.

`tests/plastic_finite_mixed_thermal.rs` specifies native pressure-condensation
parity, all mixed tangent blocks, exact thermal volume and pressure entropy,
temperature-induced yielding, activated creep/recovery and independent work
accounts. These bounded manufactured cases retain the original source domains.

## Source-bound Johnson--Cook ductile initiation

`ductile::source::resolve_johnson_cook_initiation` admits one explicitly pinned
`ConstitutiveModelCard` carried by the original immutable material card. Law
`ductile-johnson-cook-minus-triaxiality`, law/state version 1, fixes the
`exp(-d3*eta)` convention and virgin accumulated plastic strain/utilization.
All eight scalar parameters are mandatory and unit-checked: dimensionless
`d1`–`d5`, `reference_strain_rate_s_inv`, and absolute-kelvin
`transition_temperature_k`/`melting_temperature_k`. No original plus-exponent
convention, independent parameter fit, tabulated strain observation, or coupon
elongation is silently converted. The whole fit and its calibration artifact
hashes, provenance, material identity and complete typed context are retained.

Original typed triaxiality, equivalent-plastic-strain-rate and absolute-temperature
bounds stay live, including temperatures beyond the original admission point.
The numeric constructor checks the complete finite calibration rectangle;
runtime queries enforce it without source queries or allocations. Extra evolving
coordinates and a second temperature axis refuse. Other environmental coordinates
stay fixed and source binding compares them in both directions. Matching an
existing material bundle uses its retained coordinate views without allocation.

`ResolvedJohnsonCookInitiation` delegates the native response and analytical
triaxiality/rate derivatives, so existing small- and finite-strain progressive
ductility laws use the same onset chain rule. This source adapter supplies no
temperature-dependent mechanical coefficients, caloric feedback, parameter fit,
post-fracture topology, or independent experimental validation. Focused G0/G1/G3
regressions in `tests/ductile_johnson_cook_source.rs` cover the analytic response,
rate/temperature-sensitive irreversible initiation, native progressive parity,
derivatives, exact fit selection, source/context binding, and calibration refusals.

## Exact finite-support ductile completion

`DuctileWorkEvolution` selects the post-initiation profile explicitly; all existing numeric and source constructors retain `Exponential`. `FiniteSupportQuadratic` uses `x=(h/Gp) integral R d(alpha)`, `r=max(1-x/3,0)`, `g=r^2`, and `dg/dx=-2r/3` below its support, zero thereafter. Exact failure occurs at `x=3`, after undegraded equivalent work `3Gp`; `integral g dx=1`, so `Gp` remains integrated degraded equivalent plastic-work capacity, not total thermodynamic fracture energy.

Cumulative and incremental attenuation-work diagnostics use `Phi=1-r^3` and a stable factored increment. They differ from damage and its increment. Physical plastic loss, release of old active storage, signed finite endpoint-work defects, numerical loss, and heat destinations retain their existing separate meanings. Perzyna attenuation retains the equilibrium resistance `R`; physical viscous loss remains in its original account.

Exact terminal predicates cannot be triggered by rounded exponential damage or a stiffness cutoff. Full three-dimensional terminal states retain their original compressive bulk support and native histories. The finite tangent uses the exact zero profile and derivative without evaluating active spectral derivatives that no longer contribute. Explicit source builders retain the original capacity receipts while declaring an authored shape assumption; a scalar capacity claim does not establish measured support for that shape.

Terminal finite plane-stress candidates return `PlaneStressTerminalFailure` before dividing by integrity or inverting the transverse block. Small-strain reduction also refuses its singular block. No stiffness floor, arbitrary terminal thickness, automatic crack direction, deletion, erased inventory, or experimental qualification is implied. Component fracture/contact transitions require their own geometry, traction, inventory, and work admission.

Focused G0/G3 regressions in `tests/ductile_completion.rs` cover analytic capacity and stable terminal tails, full native tangents, retained compression and irreversible histories, physical-rate hold completion, exponential-default parity, and the plane-stress boundary. In this implementation session, pinned rustfmt parsed the source; Cargo typechecking and test execution remain unverified because the available workspace was full and the real dependency build had been killed by the memory limit.

## Progressive ductile plane-stress section work

`DuctileJ2Law` audits accepted increments through the original
`LayeredPlaneStressSection` reduction. Every thickness point retains its full
tensor strain and exact immutable material owner. The audit replays the single
old-to-candidate native update and verifies accepted history and transverse
equilibrium without advancing candidate history again.

Stored energy, plastic dissipation, damage release, endpoint numerical loss
and signed transverse residual work remain separate. The integrated
`Gp*delta_damage/h` diagnostic is nonadditive and is per original section area.
Every layer represents the same actual interval. Characteristic length remains
with the original material; shell thickness never replaces an in-plane band
width. Four focused tests in `ductile_section` cover analytical twisting-bend
states/work, unloading/reversal, nonsymmetric condensed tangent finite
differences and owner/time/refusal semantics.

## Progressive ductile plane stress

`plastic::finite::ductile::plane_stress` solves all three transverse Piola
tractions with positive thickness and the complete original progressive/Perzyna
return. Fixed accepted history, temperature and physical interval survive every
local retry. A normalized pivoted 3-by-3 solve produces the unsymmetrized Schur
tangent; physical traction is normalized by current integrity so attenuation
cannot conceal an unequilibrated root. All local work uses fixed-size storage.

The original candidate retains separate reference-volume storage, physical
plastic dissipation, released damage storage and signed endpoint-work defect.
The method admits a regular local transverse equilibrium with explicit
iteration/domain refusals. It does not certify branch stability, uniqueness,
necking, shell bending or fracture topology.

G1/G3 cases in `plastic_finite_ductile_plane_stress` cover independent biaxial
softening/thickness/heat/negative slope, intact native parity, unloading, all
condensed tangent components with moving onset/transverse shear/compression,
held rate relaxation and transactional failures. Synthetic coefficients are
numerical checks, without measured-alloy qualification.

## Rate-dependent progressive finite ductility

`FiniteDuctileJ2::with_overstress` explicitly selects the original fixed-temperature
finite Perzyna flow. Every trial uses its supplied physical duration for both
implicit native flow and initiation; solver iterations never consume history.
Linear, Voce and tabulated hardening retain their original metric and return.
Held overstressed deformation can relax and continue accumulating damage. The
full Piola tangent includes the inverse-rate modulus in plastic-accumulation
and moving-onset derivatives.

Attenuation retains the calibrated equilibrium-resistance work integral;
viscous overstress joins physical plastic dissipation, without changing the
supplied damage-work calibration. Storage, damage release and signed
endpoint-work defects remain separate. `validate_committed` and accepted stress
getters inspect history without another rate step or equilibrium-yield
projection. Coefficients are explicit authored, fixed-temperature selections;
no source-backed viscosity, evolving thermal coefficients, creep or fracture
topology is inferred.

G1/G3 regressions in `plastic_finite_ductile_rate` cover independent closed-form
work and dissipation, physical-time relaxation/refinement, native hardening
family parity, all 81 tangent entries and transactional refusals. These are
numerical constitutive checks, not empirical qualification of a material.

## Timed thermal strain in original tensile fracture

`brittle::source::expansion::rate` composes original integrated thermal expansion
with the native physical fracture-plane viscosity. Endpoint and path-mean force
adapters subtract free strain before the original timed material update; the
native mechanical history, derivatives and irreversible fracture remain the
same. Total path-mean stress supplies both total-displacement work and the signed
eigenstrain port. Their sum includes intrinsic fracture and separately reported
bond-viscous loss. Coefficient work remains a preceding held-history operation.

Explicit source selection retains all original fracture, Pa s/m viscosity,
instantaneous alpha(T), and optional density receipts. Thermal successors retain
the same reference, claims, geometry and jump-rate bound. Actual temperature
rate can be converted with current sourced alpha. The affine free-strain rate
of a declared interval is a distinct quantity, particularly for varying alpha;
neither is guessed from endpoint temperatures. Actual endpoint stress accepts
the declared free-strain rate separately from total deformation rate.

`transfer_temperature_over_time` supplies an explicitly restrained total-strain
experiment. It reuses the original coefficient transfer and timed path, keeping
eigenstrain work, recoverable storage, intrinsic fracture and viscous dissipation
separate through cooling, reheating and complete failure. Untimed evolution
still refuses surviving viscous bonds. These are prescribed source trajectories,
not a caloric heat/entropy model or a second-order temperature integration claim.

## Physical rate response of original bulk fracture bonds

`brittle::rate` projects the original timed cohesive law through the existing
fixed-plane crack band. Explicit fracture-plane viscosities in Pa s/m supply
surviving-bond traction, consistent endpoint and path-mean tangents, and the
analytic viscous work integral. Intrinsic fracture energy, elastic storage and
viscous loss remain separate; no second damage history or numerical relaxation
law is introduced. Complete failure removes bond viscosity while retaining the
independent compression response. `stress_at_rate` observes actual accepted
strain rate without advancing history; midpoint secant stress is not silently
reported as an instantaneous endpoint stress.

`brittle::source::rate` retains the complete original strength, fracture,
viscosity and optional reference-density receipts. The two viscosity properties
resolve from the same card and external conditions with explicit interface
units and a declared rate bound. Coefficient-temperature continuation retains
that selection, original claims and intrinsic history. Bulk Pa s coefficients,
rate-conditioned frozen values and elastic-only initiation of a viscous history
refuse. Thermal eigenstrain uses the explicit total-strain timed adapters in
`brittle::source::expansion::rate`.

Coefficients are not inferred by dividing bulk damping by band width. Intrinsic
Gc excludes the additional viscous loss, and arbitrary-width objectivity of
total viscous fracture work is not claimed. Focused regressions cover exact
work, subdivision, endpoint and path derivatives, rotation, rate observations,
release after rupture, source units, provenance and temperature continuation.

## Resistive deposition with a finite nonisothermal gas

`sublimation::gas::chamber::coating::resistance` couples existing two-temperature
coatings to the original joint molecular/gas owner. Each surface retains its
nonvolatile support mass and enthalpy, actual deposit mass and enthalpy, two
original inverse caloric charts, and an explicit finite conductance in W/K.
The deposit receives molecular/reflected energy at its own temperature. The
original resistive-storage solve shares that energy with the support through
equal/opposite endpoint interface heat; substrate mass never enters vapor stock.
Only one original gas endpoint is reconstructed for the canonical surface set.

Mass-transfer film coefficients and thermal conductance are independent inputs.
The final represented film mass participates in thermal storage, and fixed
formation offsets cancel before sensible heat is solved. Separate support and
deposit residuals must each fit the original energy tolerance. Gas species mass,
whole condensed/gas mass and energy retain the original joint admission. An
initial gas/deposit equilibrium cannot bypass gross inventory limits when an
unequal-temperature support will heat or cool the deposit. Zero thermal
conductance preserves the genuine insulated limit.

Both condensed masses must remain positive. Bare nucleation and complete film
extinction require an explicit transition and refuse here. The existing
isothermal owner continues to admit bare substrates. Frozen gas/molecular rates
and endpoint thermal exchange are first-order approximations; time refinement
is required. Conductance is supplied, not inferred from frost porosity or
thickness. No spatial temperature profile, moving interface or automatic
morphology is claimed. Focused tests compare independent two-body and ideal-gas
balances, datum changes, canonical shared-gas exchange, reflection, inventory
limits and cancellation.

### Bare-film initialization with finite thermal contact

`sublimation::gas::chamber::coating::resistance::advance_resistive_chamber_with_births`
explicitly initializes a thermally resistive film from a zero volatile inventory.
The authored zero-mass deposit enthalpy must invert to the original substrate
temperature within floating-point roundoff. It is a caloric reference, not a
mass or energy seed. The original bare-interface snapshot emits no nonexistent
deposit and accepts incoming molecules under the same shared gas inventory
ceiling. Their energy enters the actual positive final deposit mass before the
original finite-conductance kernel shares heat with the substrate. Gas-film
mass resistance and reflected-carrier heat retain their original owners.

Every bare row must receive positive representable mass in that interval.
Absent vapor, zero exposure, or sealed mass transfer refuses the birth; heating
a still-bare support requires its separate thermal owner. Existing positive-film
entry points retain their positive-inventory admission and extinction refusal.
The returned state enters those existing APIs without a mass floor or special
born-film representation.

This is a first-order frozen bare-exposure interval selected by the caller. It
does not locate a birth/nucleation event or model re-emission during that first
interval. Equal-half-step adaptivity from zero is not claimed: the newborn's
gross withdrawal fraction need not decrease when its age and trial step are
halved together. Complete extinction, nucleation barriers and general lifecycle
transitions remain separate models.

Three G0/G1/G3 integration regressions cover independent molecular, two-body and
ideal-gas balances, independent caloric datums, positive-film re-entry, reflected
heat, finite gas films, inventory refusal and transactional cancellation.


### Adaptive two-temperature deposition against a finite gas

`sublimation::gas::chamber::coating::resistance::advance` evolves the original
resistive coating owner over a complete duration. Every trial compares one full
interval with two half intervals and accepts only the conservative fine history;
it does not extrapolate states or clip species withdrawals. Both condensed
inventories must initially be positive. Bare-film birth requires its explicit
single-interval owner before this history can be admitted.

Local error includes deposit mass independently of support mass, every retained
gas species, support/deposit/gas temperatures, separate original support/deposit
specific-enthalpy increments, integrated interface heat, and gas sensible energy
from the original active-mixture caloric integral. Relative caloric scales use
changes from the current state, so a large enthalpy datum cannot relax thermal
accuracy. Gas formation constants cancel inside the original caloric owner.

One whole-call trial allowance includes coarse, rejected and inventory-refused
work, with at most three bounded chamber invocations per trial. The audit counts
attempted invocations; only accepted half steps contribute to physical streams,
reflection, interface heat and separate condensed storage. Canonical surface
order and compensated original-owner sums retain deterministic reductions.
Whole-window total/species mass and gas/support/deposit energy defects must pass
the original physical tolerances. Late cancellation or any refusal returns no
partial candidate. Only the original inventory-limit error is retryable.

The fixed-geometry entry points retain authored areas, thermal conductances and
independent mass-film coefficients. This is first-order local temporal control,
not a global error bound, source-fit accuracy claim, nucleation/extinction model,
or spatial temperature solution. G1 regressions define independent two-body
exponential and finite-gas ODE limits; G3 covers accepted integrals, datums and
surface order; G4 covers inventory retry, work/minimum-step limits and cancellation.


### Growing resistive coatings on the original finite-gas history

`sublimation::gas::chamber::coating::resistance::growing::advance_growing_resistive_chamber`
composes the shared two-temperature adaptive history with the original
mass-owned coating geometry and borrowed `ConductivityLaw`. Planar, spherical
and cylindrical layers derive current thickness and exposed area from the
actual deposit mass. Conductivity is queried at the current deposit-specific
enthalpy on its original material chart. The declared contact resistance uses
the fixed inner contact area. Source errors retain `DepositionError`; no source
law, caloric chart, effective heat capacity or molecular stream is replaced.

Each full or half interval refreshes area and thermal conductance before the
original chamber runs. Its accepted endpoint refreshes outer area and exact
mass-derived recession; the next half step therefore uses the newly deposited
mass, altered area and layer resistance. The single adaptive owner accepts only
the conservative two-half-step history. Thickness and conductance join its
independent support/deposit enthalpy increments, temperatures, deposit mass,
gas-species inventories, gas sensible energy and interface-heat error controls.
Canonical model references are prepared once and borrowed across linear walks;
material evidence is not reselected inside trial loops.

The gas volume is the caller's fixed accessible volume after excluding fixed
supports. Total deposited volume must remain below the authored fraction of
that volume, at most one percent. This is a dilute-displacement approximation,
not an error bound for omitted pressure work. Constant phase density is checked
against the original phase law at every sampled state. Conductivity is frozen
at the start of each physical interval; quasi-steady full-thickness resistance
and lumped deposit storage are a two-node approximation. No porosity, thermal
expansion, spatial profile, gas-volume work, bare birth or extinction transition
is inferred by this history owner.

G0/G1/G3 regressions in `resistance/growing/tests.rs` compare two original
physical half intervals with an independently calculated planar circuit, check
curved outer area against direct mass/volume formulas while inner contact area
stays fixed, and preserve original conductivity-domain errors, density/volume
admission and transactional cancellation.

### Explicit heating and cooling of finite-gas, two-temperature deposits

The resistive chamber accepts separate prescribed support and deposit heat through
`GasResistiveDrive` and `advance_resistive_chamber_driven`. Positive joules heat
the named condensed body and negative joules cool it. These sources enter the
original two-body caloric solve, alongside the existing molecular energy target;
they never enter the gas molecular stream or add another latent source. Every
surface and chamber receipt exposes supplied heat. The two external body inputs
remain separate in thermal receipts, even for an equal/opposite pair. Energy
defects mean condensed storage plus gas storage minus supplied heat, with
original caloric/formation datums cancelled before arithmetic. A retained
error-free addition remainder prevents an unrepresentable small deposit source
from silently vanishing into a larger molecular target.

An optional borrowed `ResistiveInterfaceLaw` replaces the authored scalar
contact inside that same endpoint solve. Gas-film mass transfer remains an
independent input. Zero duration permits an explicitly prescribed energy impulse
but no mass transfer or contact heat. Positive coating inventories remain
required by the ordinary driven entry point. Either nonzero external source
disables the stationary-chamber inventory shortcut, regardless of their sum.

`advance_resistive_chamber_adaptive_driven` takes `GasResistivePower` in watts.
Every actual full or half trial receives its own `P * dt`; rejected intervals
cannot contribute supplied heat. One reused drive map and the existing adaptive
controller carry both original caloric states and the finite gas. Both source
integrals enter the physical error norm, and the accepted whole-window source
for each body is checked independently against `P * duration`. Unknown IDs,
nonfinite/unrepresentable work, exhausted work, source failures and cancellation
publish no candidate. This is prescribed-power evolution, not an imposed
temperature thermostat, feedback controller or unbounded gas reservoir.

### Original conductivity integrals across the coating and its contact

`ConductiveCoatingInterface` binds an actual positive coating mass, the original
positive conductivity curve and the fixed inner contact area. The original
material/film owner integrates conductivity across the actual coating temperature
drop and solves the coating-side contact-face temperature. Positive power leaves
the support and enters the deposit; the contact jump and coating flux share that
same power. The two-temperature deposition kernel uses this law inside its
existing caloric balance instead of freezing conductivity at a bulk trace.

The adapter reuses the original planar/radial geometric multiplier, including
mass-derived sub-ULP shell thicknesses. Areal contact resistance is divided by
the fixed inner area once. The deposit and solved coating face must remain in
the conductivity source domain. A finite contact jump may admit a support
temperature outside that domain without treating support material as coating.
Perfect contact requires both material traces to be supported. The caller owns
the identity of the borrowed conductivity curve and its consistency with the
deposit; borrowing knots alone does not establish source/material identity.

Binding and successful flux queries allocate no storage and perform no
material-card resolution. Original knot limits and the explicit contact-solve
budget bound work; a borrowed checkpoint can interrupt original inner trials.
The richer response retains heat continuity, actual face temperature and
endpoint derivatives, while the low-level thermal trait preserves cancellation
and budget refusal. This is source-free, isotropic, quasi-steady conduction
through a layer with two original lumped caloric inventories. It does not add
another heat store or resolve a spatial transient temperature field.

### Growing layers with nonlinear original conductivity and prescribed power

`GasIntegratedGrowingCoating` and `advance_growing_integrated_chamber` combine
the original mass-owned planar/spherical/cylindrical geometry, conductivity
integrals and finite inner contact with the same two-temperature finite-gas
history. Only geometry is frozen within a physical interval. Every thermal
endpoint evaluates the original conductivity curve over the actual coating
temperature drop and solves its contact face. Actual deposited mass refreshes
area and thickness between half intervals. The reported equivalent circuit uses
the original average conductivity across only the material segment; it does not
subtract nearly equal total/contact resistances or freeze a trace conductivity.

The existing adaptive owner retains molecular streams, both independent caloric
inventories, one gas reconstruction and one equal/opposite interface transfer.
Thickness and source-integrated secant conductance participate in its temporal
error estimate. Optional signed support/deposit powers are integrated over each
actual trial interval, and only accepted half intervals contribute source heat.
Each accepted body source is independently checked against prescribed power
times the full requested duration. The layer owns its interface law; a competing
interface in the power map is refused. Piecewise power schedules use successive
explicit time windows.

Original conductivity knots, finite inner and outer solve budgets, and the same
caller cancellation gate bound work. Source failures or late cancellation return
no partial geometry, heat history or gas candidate. The earlier frozen-law
growing entry point retains its original path without allocating nonlinear
interface objects or source maps. Both paths retain positive-film, constant-density,
shape-preserving and fixed-accessible-gas-volume assumptions. The deposit-volume
ceiling bounds omitted displacement, not pressure-work error. This extension
does not resolve nucleation/extinction, porosity, densification or a transient
temperature field through the coating.

### Live gas-film transport in the original deposition history

`resistance::film::advance_resistive_chamber_adaptive_resolved` and
`resistance::growing::advance_growing_integrated_chamber_resolved` accept a
bounded source resolver for mass-transfer coefficients. Every physical full or
half interval supplies the current original gas, canonical surface identities
and both caloric states. Growing layers first rebind the actual mass-owned
exposed area. The second half interval therefore resolves against its actual
updated gas and deposit, rather than a coefficient frozen over the full window.

The existing adaptive owner, molecular heat/mass streams, nonlinear coating
contact, prescribed-power accounting and single gas reconstruction are retained.
Only accepted half intervals contribute physical history. The fixed-film APIs
use the same interval seam without allocating a source map. A resolved history
reuses map nodes; each retained value is invalidated before source evaluation
and must be overwritten or explicitly removed. Omitted IDs use the original
free-molecular limit. Unknown IDs, negative/nonfinite coefficients and stale
outputs refuse through the original chamber admission.

Resolution also runs at zero-duration admission. The resolver must be
deterministic for physical inputs, bounded, cancellation-aware and free of
trial-committed physical side effects. Original source errors remain typed in
`GasResistiveFilmError::Source`; core evolution errors remain separate. Either
failure or late cancellation returns no partial gas or coating candidate.
The coefficient is still frozen within each physical interval; the existing
full/half comparison controls coupled temporal error. No extra spatial gas
film, source extrapolation or continuum correlation is inferred by this seam.

Focused G0/G3/G4 tests compare live rebinding with two original physical
intervals, exercise changed gas inventory and mass-owned area, and cover stale
outputs, original source errors, identity refusal and cancellation.

### Closed support/coating/gas phase-equilibrium endpoints

The resistive chamber's equilibrium module supplies an isolated thermodynamic
endpoint for one nonvolatile support, one volatile coating and the original
finite multicomponent gas. Positive exposed area, scalar support/coating
conductance and molecular accommodation are required exchange permissions.
The endpoint imposes a common temperature and partitions the conserved
volatile mass by the original saturation law: positive deposit coexists with
saturated vapor; insufficient volatile inventory leaves an exactly bare,
undersaturated support. No mass tolerance deletes a meaningful residual film.

A caller-declared support-specific-enthalpy bracket selects one increasing
branch. Every nonstationary trial uses the original phase source and gas
partial enthalpy, derives coating enthalpy with the fixed datum conversion,
and verifies that value through the original coating inverse. The residual
uses original species enthalpy increments so formation constants cancel before
energy admission. Stored caloric encoding has its own energy gate. Both body
datums remain independent. An exact stationary endpoint is preserved without
probing unused root endpoints.

The original gas inventory/known-temperature owner reconstructs the accepted
composition, pressure, transport and acoustic state at most once. Actual species
inventories, total mass, source-relative energy, common temperature and phase
complementarity are admitted before publishing the candidate. Exact depletion
of one gas species retains its source catalog; complete gas vacuum refuses.
Explicit evaluation limits and the existing caller cancellation gate bound work.

The fixed accessible gas volume excludes supports. Actual initial and trial
deposit volumes must remain below the explicit fraction ceiling in (0, 0.01].
This bounds omitted displacement, not pressure-work error. The endpoint is not
a transient molecular step, a nucleation/extinction clock, a global free-energy
minimizer or a morphology update. It adds no external heat, pore model or
kinetic equilibration time.

## Source-bound finite membrane entropy and thermal stress

`HyperelasticMembrane::thermoelastic_state` prepares the original expansion
path's exact local thermal jet, reduces the original finite-compliance potential
through transverse equilibrium, and applies the existing engineering Green-strain
chain rule. `frozen_visco_thermoelastic_state` adds the complete explicit native
Maxwell history while retaining the exact source equilibrium and natural map.
Neither query advances time, changes accepted memory or creates heat.

Mechanical entropy, internal energy, heat capacity and thermal resultants use
the original reference thickness and area. The entropy-constrained tangent adds
the thermal outer-product correction using an explicit positive caloric capacity
plus the mechanical contribution exactly once. This numeric caloric argument
does not grant source admission or define its thermodynamic reference.

`SheetExpansionReference::thermal_reference_jet` reuses the source owner's
`at_temperature_with_smooth_tangent`: one binary search checks at most two
adjacent original segments. True slope changes and source endpoints refuse a
two-sided derivative; redundant datum knots remain valid. Fixed directional
source permutations retain their exact axes. Existing three-dimensional
thermodynamic consumers share this source derivative admission, including their
explicit one-sided finite-heating endpoint convention.

Exact-volume thermal potentials and reciprocal transient heat exchange remain
separate owners. Four source-backed regressions compare relaxed energies and
thermal resultants, entropy-constrained tangents, source support and retained
nonproportional viscous memory. Source fixture model domains now retain the same
typed Kelvin descriptor as their original expansion queries.

## Finite-strain thermodynamics with retained viscous memory

`ViscoHyperelastic::frozen_thermoelastic_state` derives mechanical Helmholtz
energy, entropy, internal energy, fixed-deformation capacity and thermal Piola
stress from the original zero-time physical potential at fixed accepted viscous
metrics. The shared natural-map chain rule retains noncommuting reference
derivatives and the existing elastic-source coefficient derivatives. Current
law/reference/history admission remains native. Positive-duration incremental
potentials and relaxation loss are not interpreted as reversible storage.

`FiniteThermoelasticPlaneStressState` uses the original full-compliance
transverse solver and its stable Schur block. Eliminating all transverse
tractions changes thermal stress and adds the required positive Schur correction
to the mechanical heat capacity. Elastic, source-derived and fixed-history
viscoelastic snapshots share this reduction without a second transverse solve.
Every density remains per original volume; no relaxed thickness becomes mass.

The complete caloric capacity must be admitted separately. Exact-volume pressure
thermodynamics and thermal/history evolution are outside these local snapshots.
Native stable finite energy increments remain available for elastic snapshots;
fixed-visco endpoints explicitly refuse that separate increment operation until
its complete history-dependent owner exists.

`hyper_thermoelastic_plane_stress` and `visco_thermoelastic` specify independent
thermal derivatives, Maxwell reciprocity, objectivity, nonproportional retained
history, relaxation-time independence, native plane reduction and refusals.

## Supported-film lifecycle with finite nonisothermal gas

`sublimation::gas::chamber::coating::lifecycle::advance_supported_lifecycle`
advances one original `GasSupportedSurface` through zero-film formation and
complete extinction. The thermal molecular-turnover owner integrates the
finite reservoir's changing backpressure at frozen initial support/gas
temperatures, area and optional dilute-film coefficient. Both integrated
opposing streams then enter the original supported caloric storage and the
same shared source-relative gas transaction. Composition, pressure, transport
and gas caloric state are reconstructed once; support mass never enters vapor.

The interface explicitly admits barrier-free heterogeneous deposition and
instantaneous re-emission. A bare undersaturated support re-emits every sticking
incoming molecule at its own temperature, giving equal nonzero streams and
sensible heat exchange without negative film inventory. Supersaturation forms
a film with concurrent emission. After an actual zero-film hit, the remaining
interval follows the same bare-interface rule. An initially absent vapor may
be emitted, mixed and reabsorbed within one interval using its retained original
caloric source at gas temperature. No first-step emission suppression, arbitrary
trace-mass seed or latent-energy double count is introduced.

`GasSupportedLifecycleSolve` carries the existing gas and supported-storage
solvers, source-compatibility tolerance and an explicit deposit-volume ceiling
in `(0,0.01]`. Gross streams may recycle inventories many times, so this owner
uses a private single-surface integrated-inventory admission instead of the
ordinary chamber's old-inventory gross bounds. It independently admits net
inventories, stream/mass representation residuals, source identity and complete
gas/support energy balances. The ordinary multi-surface chamber retains its
original gross admission. More than one independently integrated surface or
additional reflected-heat accommodation refuses rather than sharing the same
finite vapor twice or freezing its changing pressure in a second heat model.

The candidate retains the original supported-chamber endpoint plus the exact
frozen-temperature turnover receipt and both endpoint deposit-volume fractions.
Its turnover reservoir is an intermediate kinetic state; `gas()` is the final
thermochemical state. Positive asymptotic inventory tails that cannot survive
the original net-transfer representation refuse, even below mass tolerance;
they are never relabeled as an extinction event. Zero duration still admits
sources and volume and preserves state. Cancellation or any unsupported source,
unrepresentable inventory, volume or caloric balance publishes no partial state.

The coupled step is first order in evolving temperatures; its kinetic extinction
time is a local event estimate, not an exact nonisothermal clock. This model
has one internally isothermal supported coating, fixed exposed area and fixed
accessible gas volume, with specular nonsticking/carrier reflection. The volume
cap bounds omitted deposit displacement, not pressure-work error. Nucleation,
adsorption residence, porosity, gas-film heat storage, multiple competing fronts,
finite support/coating thermal resistance and resolved morphology remain outside
this opt-in model. Focused G0/G1 tests specify analytical finite-vapor growth,
extinction and continuing bare sensible exchange, same-interval species reentry,
the original small-step limit, datum invariance and transactional refusals.

## Adaptive supported deposition with finite gas

`sublimation::gas::timeline::coating::adaptive::advance_coatings_adaptive_with_gas`
advances the original supported-coating, heat-transport and closed-gas owners
on one event clock. It compares one complete interval with two half intervals
and accepts only the conservative fine history, without extrapolating physical
inventories. The shared `GasTimelineAdaptiveSolve` controls bound local errors
in gas and substrate temperatures, every gas-species inventory, coating mass,
and substrate/coating specific enthalpy. Deposit mass has its own error scale;
heavy nonvolatile supports cannot conceal small-film errors. Relative caloric
scales use increments from the interval's initial state, not arbitrary datums.

Trials stop at birth/contact events; new surfaces acquire no preceding exposure,
and empty pre-birth spans jump to their next event. Accepted heat, incoming
material and molecular streams are accumulated once in canonical cell order.
Every completed original thermal substep, including rejected/coarse probes,
spends one whole-call work allowance. Only gross-inventory refusals are retried
with a shorter interval; source-domain failures, exhausted budgets and
cancellation preserve their original refusal and publish no partial state.
Whole-window acceptance separately checks each deposit's integrated mass,
all gas-species inventories, unchanged support mass, and gas/solid energy defects.

This is a first-order local temporal estimator, not a global accuracy bound or
experimental calibration. Original fixed reference areas, internally isothermal
thin coatings, explicitly supplied molecular/film resistance and reflected heat
remain the physical model. Spatial frost structure, nucleation and moving
geometry are not inferred. The focused coating regressions include an independent
simultaneous-rate reference, birth shielding, datum/mass-scale invariance, and
shared-budget/cancellation behavior.

### Integrated coating lifecycle on the original timeline

`timeline::coating::advance_coating_lifecycle_with_gas` and
`adaptive::advance_coating_lifecycle_adaptive_with_gas` select the original
supported-film lifecycle interval inside this SAME heat/event/controller owner.
Exactly one coating is registered; other nonvolatile thermal cells and contacts
remain available. The registered original mass-transfer film, caloric charts,
source identity and event-clipped exposed area are retained. Bare growth,
prompt molecular return and complete coating extinction no longer depend on
withdrawing only a fraction of the old film. Its nonvolatile substrate, heat
sources and contacts continue after its coating reaches exactly zero.

The same whole-window species and energy balances, accepted-fine histories,
local state norm, cancellation and total attempted-heat-work budget apply.
Scheduled births receive no earlier exposure; endpoint births and zero-duration
calls still admit original sources and the declared deposit-volume fraction.
An integrated phase step reconstructs the original gas once. Coarse and rejected
probes never enter accepted molecular or external-source totals.

The opt-in zero-residence heterogeneous boundary and specular nonsticking/carrier
assumptions of the lifecycle interval remain explicit. Temperatures evolve
through the original first-order heat/molecular split; an analytical event of
the frozen-temperature mass subproblem is not an exact thermally coupled event.
No multi-surface analytical reservoir, nucleation barrier, adsorbate storage,
finite coating resistance, moving morphology or global error bound is claimed.
The five G0/G1 regressions specify bare undersaturation, growth followed by a
scheduled heat pulse and extinction with continuing contact, accepted manual
half-step equivalence, delayed birth and transactional source/cancellation
admission. Their independent storage check reads actual thermal and gas states.

## Thermal finite porous compression

`plastic::finite::porous::thermal::AffineThermalFinitePorous` extends the
original finite pressure-sensitive return with explicit affine temperature
dependence of shear/bulk moduli, yield and energetic hardening, including
optional Voce saturation. The ellipse parameter and saturation rate stay
fixed. Isotropic logarithmic expansion gives `Fm=F/theta(T)`; the native
non-isochoric plastic metric remains intact. Elastic volume is
`Je=J/(theta^3*Jp)`. The exact total-Piola tangent scales the native tangent
by `1/theta^2`, retaining its geometric and nonlinear hardening terms.

Original-volume free energy supplies fixed-history entropy, internal energy,
the mechanical heat-capacity contribution and the Piola/temperature derivative.
The complete temperature interval and native strain/history limits are admitted
at construction. Caloric owners must combine their own capacity with the
declared mechanical lower bound; that bound is not a uniqueness theorem for
the coupled plastic/thermal return. Old history is checked at its old temperature
before an endpoint-temperature return, so heating can cause further compaction
at unchanged total deformation. Physical yield heat, signed reversible caloric
transfer and signed mechanical endpoint-work defect remain distinct.

`ResolvedAffineThermalFinitePorous` binds exact affine or affine-Voce thermal
model selectors, law/state v1, with all eighteen base parameters mandatory and
three additional Voce parameters when selected. One typed absolute-temperature
source axis covers the full family interval; original card, reference context
and fixed external conditions persist. `FiniteThermalPorousMaterial` exposes
the complete family and cannot bypass thermal mapping through an isothermal
material adapter. The seven `plastic_finite_porous_thermal` regressions specify
independent hot hydrostatic compaction/unloading, noncoaxial tangent and thermal
derivatives, heat/storage bookkeeping, isothermal and zero-saturation limits,
temperature/history admission and source retention/refusal.

This is a dry, symmetric, rate-independent thermal extension. It infers no
evolving pore fraction, cell-wall densification, pore-fluid pressure, phase
transition, separate Taylor--Quinney fraction or measured calibration.

## Source-prepared natural strain for finite membranes

`sheet::hyperelastic::expansion` binds a reusable `PreparedFreeExpansion`
engineering-length path to the original sheet. Every temperature endpoint is
evaluated relative to the retained original datum; expansion ratios never
compound and the original mass, thickness and property receipts never refresh.
The exact selected hyperelastic model is admitted at the endpoint temperature.
The source card, typed temperature coordinate and fixed external conditions
remain bound to the prepared path.

Isotropic sources and explicitly aligned fixed-principal sources are supported.
Directional mappings must be proper signed axis permutations into the fixed
warp/weft/normal frame. Unsupported thermal shear, inconsistent frame identities,
foreign cards, ambiguous temperatures and out-of-domain endpoints refuse.
Preparation runs outside native mechanical iterations and shares the immutable
source path across endpoints.

`with_hyperelastic_free_expansion` selects the path once;
`at_hyperelastic_temperature` prepares later endpoints. Accepted viscous states
must continue through their native reference-change transaction and its separate
signed work accounts. This is prescribed thermal strain preparation, not a
caloric or entropy closure; native viscoelastic reference events still require
unchanged equilibrium and ordered branch energies.
`tests/sheet_hyperelastic/expansion.rs` covers retained datum/inventory,
directional and exact-volume free shapes, source/frame refusals and continuation
of the original native Maxwell memory.

## Prescribed membrane strain-reference changes with retained memory

`ViscoHyperelastic::change_membrane_reference_from` preserves the installed
in-plane Green strain and original material-frame Maxwell metrics while the
native plane-stress reference transaction changes the prescribed natural
configuration. Both finite compliance and exact volume reuse the accepted
transverse endpoint for consistent tangent condensation, without a second
transverse Newton solve. No time passes and no relaxation heat is generated.

Original-area storage, fixed-full-deformation reference work, pressure-volume
constraint input, native transverse loading and physical endpoint transverse
work remain separate. Their signed finite-event integration defect must be
gated by the spatial owner. Ordered branch energies, equilibrium law and
elastic domain retain the native compatibility rules; changing a material
reference does not rename coordinates or supply a thermal evolution law.
`tests/visco_membrane_reference.rs` checks retained nonproportional history,
work closure, all tangent columns, kinetic-only changes, unused-bulk
independence, continued relaxation and refusal atomicity.

## Finite-strain membrane relaxation with original material history

`visco_hyper::membrane::FrozenViscoMembrane` reduces the original finite-strain
Maxwell potential to engineering Green strain through the same exact Cholesky
chain rule as elastic sheets. Every trial borrows one accepted history and one
physical interval. Finite compliance relaxes all transverse tractions; the
explicit exact-volume variant retains native remainder history and its physical
pressure multiplier. Material-frame viscous metrics survive nonproportional
loading and spatial rotation without being reset to a new reference.

`step_response` proposes complete native memory, physical stress, reference-area
storage, loading work and nonnegative relaxation loss. These remain separate
from the incremental sheet potential and any spatial endpoint-work error.
Original thickness and natural-reference volume determine work; no current-area
or current-thickness substitution creates mass or heat. The adapter supplies no
new source qualification, bending, wrinkling, damage or evolving temperature.
`tests/visco_membrane.rs` covers retained relaxation, consistent tangent,
independent exact-volume energy, natural-reference volume and refusal atomicity.

## Exact source heat capacity through solid and phase transport

`ResolvedHeatCapacityCurve::linear_heat_exchange` uses the existing exact
affine-Cp exchange kernel for fixed-temperature reservoirs and steady channels.
Original source Cp knots are retained once at admission and across caloric
joins. Queries allocate nothing, visit only crossed original intervals and
ignore optional enthalpy-chart subdivisions. Relative heat, the exchange mean
and inlet/wall derivatives share that same source integral. Source pressure,
composition and transfer-model assumptions remain explicit; the adapter adds
no phase change, hydraulic feedback or independent pressure-work account.

Source heat-capacity ingress retains the original affine `Cp(T)` on every
generated sensible interval. Its enthalpy integral, quadratic temperature
inverse and `dT/dh = 1/Cp(T)` share the same immutable curve. The public boundary
chart still honors the existing subdivision allowance and knot budget; runtime
calorics no longer replace the source law by those boundary chords. Explicitly
authored h/T charts keep their authored interpolation, and a latent plateau
retains exactly zero temperature slope and its single original latent account.

`EquilibriumEnthalpyPhaseCurve::try_from_heat_capacity_exact` and
`ResolvedHeatCapacityCurve::try_new_exact` retain only the original Cp/density
source-knot union. They use the same analytical calorics and source admission,
without generating enthalpy-chord subdivisions or evaluating their error
diagnostics. The supplied source-knot budget remains explicit. Existing
constructors still admit an optional finer public boundary chart; consumers
of the compact path use caloric methods rather than interpolate boundary
chords. `material_surface::caloric` checks crossed source cells, original
receipts, exact heat/entropy/density, support refusals and finite storage that
an unnecessary chord-error intermediate could not represent.

`sensible_heat_capacities_on_interval` identifies retained source intervals.
`specific_enthalpy_on_sensible_interval` supplies their original temperature
inverse, and `sensible_heat_capacity_bounds_j_kg_k` supplies positive minimum
and maximum capacities for coupled stability/admission bounds. None of these
methods extrapolates or supplies missing material data. Source context,
receipts, reference datum, and Cp metadata remain part of curve identity.

Adjacent joins, source melting charts and fixed-mass composites retain the
original capacities and mass-weighted sensible energy. Temperature-coordinate
relative heat integrates Cp directly; entropy integrates Cp/T with stable
small-step kernels and the deterministic logarithm owner. The finite entropy
mean and both analytic endpoint derivatives use that same nonlinear T(h).
Density remains its original rho(T); the mixing-volume remainder includes the
Cp-dependent curvature of h(T). Crystallization admission checks interior
quadratic enthalpy-gap minima as well as the merged temperature boundaries.

Prepared phase conductivity places original k(T) breakpoints at their true
h(T), evaluates k at the original chart temperature and uses its analytic chain
rule. A clipped upper endpoint retains the prepared interval's left derivative,
including next to a latent plateau. Required unrepresentable source coordinates
or tangents refuse. Evaluation performs no new source selection, hashing or
allocation. Work remains bounded by original source/caloric intervals; no
adaptive table resampling is added.

Focused regressions cover coarse boundary charts, both Cp slopes, energy-datum
shifts, entropy means, original conductivity knots, latent joins, composites and
crystallization minima. Numerical-leaf execution does not establish full native
consumer compilation, measured-material validity, pressure evolution, or
cross-module integration-test success.

## Linear kinematic hardening of anisotropic Hill shear

`plastic::hill::Hill48::with_linear_kinematic_hardening` explicitly adds the
Prager backstress X=C*eps_p in fixed material axes to the existing Hill owner.
C is a nonnegative Pa coefficient on tensor plastic strain; the J2 uniaxial
hardening convention would use C=2*H_k/3. Existing constructors and material
source adapters retain C=0. Callers select independent kinematic data on the
native law before retaining history. The original five-mode return and
HillState are reused; no independent backstress fit or second solver is added.

Yielding uses q(sigma-X), while observations return physical Cauchy stress.
The consistent tangent retains the physical elastic driving modulus and the
same original orientation transforms. The backstress shifts reverse shear
yielding in all three material directions. Recoverable Prager storage is
C*eps_p:eps_p/2; physical dissipation remains Y*delta_alpha, separate from
both hardening storage and backward-Euler numerical loss. This is an energetic
model interpretation, not a calibrated heat fraction, nonlinear cyclic
hardening, saturation, ratcheting law or finite-strain texture evolution.

Seven focused tests in `tests/hill_kinematic.rs` check directional reverse-yield
thresholds and closed-form tangents, nonzero pure-kinematic shear stiffness,
all 36 mixed-loading tangent entries, cyclic energy closure, exact C=0
compatibility, the original oriented owner and invalid-modulus refusal.
The convention follows LAME's elastic-plastic model, equation 4.27:
<https://www.sandia.gov/files/sierra/SM_LAME_5_28/main/models/elastic_plastic.html>.

## Bounded circular torsion and exact layer resultants

`ConcentricShearSection::torsional_stress_at_twist_bounded` extends the
original `torsional_stress_at_twist` observation with an explicit small-shear
limit at the selected layer's outer radius. It returns the same
`section::torsion::TorsionalStressPoint`, delegating point stress, tensor
conventions and representability checks to that canonical law. Original G
remains cached at source admission; neither observation allocates source keys.
Both sides of a material interface retain their separate G and stress, with
no averaged material or transverse-shear correction factor.

`torsional_layer_response` gives the analytical outer-radius shear maximum,
signed layer torque GJk and reversible energy per length GJk^2/2. Their layer
sums recover the original circular rod stiffness and storage. An explicit
small-engineering-shear bound applies even when querying near the axis.
Invalid layers, bore/exterior points, nonfinite inputs, and unrepresentable
outputs (including positive layer storage underflow) refuse without changing
source state. Layer GJ uses the same multiplication order as native admission.

`tests/section_shear.rs` checks analytical solid/tubular and layered resultants,
torque/energy sums, reversal, explicit zero twist, exact delegation and hidden
surface overstrain. The original torsional point-law tests remain unchanged. This is
Saint-Venant circular torsion with perfect bonding;
it does not recover transverse-shear or bending stress, end effects, noncircular
warping, plastic redistribution, shear strength or fracture.

## Native implicit conduction from volumetric heat capacity

`thermal::transport::volumetric::timeline::FixedVolumeThermalLaw` implements the
existing `ConductivityLaw` and `ImplicitConductionLaw` directly. Its original
source-integrated storage is `e = integral C_vol(T) dT / rho0`; the analytic
tangents are `dT/de = rho0/C_vol` and `dk/de = (dk/dT)*rho0/C_vol`. Conductivity
derivatives use original conductivity intervals, so a nearby capacity-only
breakpoint cannot create a spurious slope. Knot derivatives are right-sided,
except at the upper support endpoint. Runtime queries allocate and hash nothing.

The native implicit chain, graph, SDIRK2 and adaptive owners retain time stepping,
sources, signed pair heat, residual gates and cancellation. The same caloric
object retains mass/volume admission, source receipts, datum and optional
sensible entropy. `DeclaredSensibleStorage` remains distinct from Cp, Cv, stream
enthalpy and phase-change storage. Fixed geometry and the explicitly selected
scalar conduction direction remain required.

`tests/volumetric_implicit.rs` covers analytic source tangents, independent
nonlinear two-body heat balance on chain/graph owners, reference-mass/datum
changes, crossed capacity knots, source boundaries, mass mismatch and
cancellation. These are synthetic G0/G1/G3 model checks, not experimental
material validation.

## Globally convex finite-compression bulk storage

`hyper::Volumetric::SimoTaylor` explicitly selects
`U(J)=kappa/4*(J^2-1-2*ln(J))`. Its positive curvature
`kappa/2*(1+1/J^2)` holds for every positive volume, and its energy diverges
under volume collapse. It shares the existing Gent, Ogden and HGO energy,
Piola, tangent, plane-stress, reference-deformation and mixed consumers.
Material cards select it with `volumetric_kind=2`; existing selectors retain
their original meanings. No experimental calibration or default replacement
is implied. Reference: Hartmann & Neff (2003), Table 4,
doi:10.1016/S0020-7683(03)00086-6.

The original-material pressure dual uses an analytic scaled, rationalized
inverse under compression. Its signed volume increment and stored energy are
resolved near the reference even when the volume rounds to one. The scalar
discrete-gradient owner uses the same original bulk energy, a cancellation-safe
secant and its exact endpoint derivative on the caller's positive bracket;
it has no log-squared `J<e` restriction. Interval work pressure remains distinct
from endpoint pressure. Existing stretch, work and iteration limits still apply.
G0/G1/G3 tests in `hyper_simo_taylor`, `hyper_mixed`, `hyper_discrete_bulk`,
`hyper_discrete_bulk_reference` and `hyperelastic_gent` check analytical
stress/curvature, pressure inversion, mixed condensation, reversible finite
compression work, held references, source-card selection and rotated tangents.


## Local shear stress in layered circular torsion

`ConcentricShearSection::torsional_stress_at_twist` and
`torsional_stress_at_torque` recover the elastic Saint-Venant torsion field
from the original per-layer G and summed GJ. In the local x-axial frame,
`tau_xy=-G*kappa*z`, `tau_xz=G*kappa*y`; tensor strain halves the engineering
shear coordinates and point storage is their work-conjugate elastic energy.
The section's transverse shear correction never changes this torsional field.
See MIT 1.050, chapter 8, “Stresses/Deflections Shafts in Torsion”:
<https://ocw.mit.edu/courses/1-050-solid-mechanics-fall-2004/8c30fdf15d9d50c02f4215903c847d7f_emech8_04.pdf>.

Explicit layer selection admits either material side of an interface and refuses
the bore, exterior, wrong indices, nonfinite inputs and unrepresentable output.
Evaluation is bounded and allocation-free, changes no history, and retains the
original frozen source state. G0/G1/G3 tests in `section/torsion/stress.rs` and
the existing section tests cover signed stress, tensor work, annular integration,
material jumps, subdivision and independence from transverse shear correction.
This is small-strain elastic torsion of perfectly bonded circular layers;
axial/bending stress, transverse-force shear, viscous history, end effects,
slip and strength/failure remain separate models.

## Circular transverse-force stress with original isotropic sources

`section::transverse_stress::CircularTransverseShearRecovery` prepares the
Saint-Venant flexure field for an elastically homogeneous solid circle or
concentric tube. Original explicit E/nu, E/G or K/G roles supply Poisson ratio
and uncorrected G. Physical layer subdivision retains its original materials
and density; unequal elastic coefficients or an undeclared isotropic law refuse.
The analytic field satisfies section equilibrium, isotropic compatibility,
the applied transverse resultant and traction-free bore/exterior contours.
See Di Paola, Pirrotta and Santoro, Acta Mechanica 217 (2011), sections 2 and
5.2, Eq. (28), <https://doi.org/10.1007/s00707-010-0376-8>.

Preparation is linear in the bounded layer count. Each point query is constant
work, allocation-free and returns physical tensor shear, engineering strain,
local elastic storage and the exact whole-selected-annulus strain maximum.
An explicit small-strain limit therefore cannot be bypassed by querying a
stress-free point. The native authored shear correction is neither changed nor
used to infer point stress; recovered storage must not be added to or silently
identified with that reduced model's energy. This quasistatic prismatic-section
recovery excludes load-introduction effects, cross-sectional inertia, elastic
heterogeneity, plasticity and a failure criterion.

`tests/circular_transverse_stress.rs` defines independent Cartesian/Laurent,
equilibrium/compatibility, integrated force/torque/energy, auxiliary Poisson
ratio, subdivision, source-refusal and strain-envelope regressions. Native
rod-fatigue consumers use this same prepared owner. These synthetic checks do
not establish measured component strength or fatigue life.

## Frozen brittle and cohesive stiffness

`CohesivePoint::frozen_tangent_pa_m` differentiates the original recoverable
storage while holding irreversible damage and its onset mode mix fixed.
`BrittleCrackBand::frozen_tangent` projects that same stiffness through the
original crack-band modes and parallel elastic remainder. Native remaining
integrity is retained directly, including values lost by subtracting rounded
damage from one. Neither operation advances history, dissipated work or time.

Opening and shear lose transfer exactly at complete failure; normal compression
and in-plane bulk elasticity remain. The bulk tensor is positive semidefinite
under tensor-shear work weights, with physical null directions preserved and
no residual stiffness. Unequal opening/compression stiffness at zero normal
jump returns `CohesiveError::NormalClosureSwitch`; an arbitrary active branch
is never presented as a unique derivative. Owner and strain admission remain
the original brittle checks. Both fixed-size observations allocate nothing.

This is a conditional frozen-history stiffness for vibration/wave consumers,
not the tangent of ongoing fracture loading. Finite oscillations must avoid
further damage and opening/closure switches. It supplies neither geometric
prestress stiffness, crack-face contact, damping nor experimental authority.
Focused G0/G1/G3 regressions check unloading/storage, mixed-mode anisotropy,
exact failed null directions, compression, closure refusal, retained history,
rotation and work symmetry.

## Completed brittle-plane release admission

`BrittleCrackBand::release_failed_plane` admits an explicitly selected physical
face only when the original native history has completely failed, the face
normal aligns with that fixed material plane, and the actual full-stress
traction meets its declared Pa gate. Active normal compression refuses even
with a loose traction tolerance. Surviving bonds cannot be released by a
numerical tolerance.

The returned observation preserves full stress, in-plane elastic storage and
all cumulative fracture dissipation, including thermal history. This bounded,
allocation-free query changes no state and creates no new cohesive energy or
impulse. Structural owners must retain those accounts once through topology
transfer. A generalized crack-band opening is not a geometric face gap; this
entry does not admit partial-damage surface insertion or contact response.

## Native brittle mean-stress assembly adapter

`BrittleCrackBand::mean_small_strain_law` borrows the original crack-band owner
for native small-strain assembly. Its force response and algorithmic tangent
are the analytic path-mean stress and its derivative at fixed accepted history;
its state update remains the original physical endpoint `BrittleState`. Mean
response and endpoint stress must remain distinct in receiving observations.
A valid path without a classical derivative retains `None` in the checked API
and a nonfinite tangent in the infallible trait. This includes virgin held
zero paths and nonsmooth retained maxima. The adapter never substitutes an
endpoint or numerical tangent; an explicitly admitted nonlinear-solver fallback
owns such cases. Original source identity, strain domains and work gates remain.

## Exact brittle path-mean stress and derivative

`brittle::path` projects the existing analytic cohesive path mean into the
original tensor-strain crack band, including its positive in-plane remainder.
`trial_with_mean_stress` performs one native cohesive evolution and returns
the unchanged endpoint plus the exact affine-path stress average.
The average's double-contracted work must match original material work within
the declared volumetric tolerance. Rupture overshoot, unloading, normal
reclosure and nonproportional first initiation retain their original histories.

`trial_with_mean_stress_tangent` additionally projects the native derivative
with respect to endpoint strain: half the parallel elastic stiffness plus the
cohesive mean derivative and original tensor-shear column weights. It retains
first-onset BK sensitivities and returns None at the native unresolved kinks;
the endpoint tangent is not substituted. Force-only calls avoid derivative work.
Both operations remain fixed-size and allocation-free and share original
owner, Frobenius-strain and energy admission with the ordinary endpoint trial.

Eight G0/G1/G3 tests in `tests/brittle_path.rs` cover exact elasticity,
fracture work through overshoot, unload/reclosure, mixed-mode derivative finite
differences, rotation, native endpoint equality and transactional refusals.
This provides work-conjugate forces for time-integration consumers; it does not
supply mass, advance time, add damping or claim an FE time-stepping method.

## Gent finite-chain elasticity through native material consumers

`hyper::HyperelasticModel::Gent` retains the explicitly supplied shear
modulus, dimensionless limiting invariant increment `Jm`, and original
quadratic or log-squared bulk potential. Its decoupled energy is
`-mu*Jm/2*ln(1-(I1bar-3)/Jm) + U(J)`. This is the declared compressible
extension of Gent's incompressible law (1996, doi:10.5254/1.3538357), not
an inferred bulk response or a fitted default. The existing exact-volume
membrane and mixed owners select incompressibility explicitly.

The original objective small-strain invariant and a bounded logarithm-ratio
series retain small recoverable energies and the large-`Jm` limit without
forming `mu*Jm`. Stress and tangent differentiate this same potential. Every
checked evaluation rejects `I1bar-3 >= Jm` before differentiation; the
existing principal-stretch restriction also remains active. No stiffness
cap, damage, rupture, hysteresis or material-specific calibration is inferred.
The native local plane-stress solver retains its stable-seed and iteration
limits; this law does not establish global branch uniqueness.

The `gent-decoupled` v1 source card requires dimensioned `mu` and `kappa`,
dimensionless `locking_invariant` and `stretch_limit`, and an explicit
`volumetric_kind`. `ResolvedHyperelastic` retains the original card and fixed
source context for native solid, mixed and sheet use; graph registration
retains the original law, version and validity declaration.
`tests/hyperelastic_gent.rs` checks independent shear
and uniaxial formulas, chain refusal, frame objectivity, tangent differences,
the Neo-Hookean limit, native mixed condensation, compressible plane stress
and a source-bound exactly incompressible sheet. These are G0/G1/G3 model
checks with synthetic coefficients, not experimental validation.
The generic native model card therefore carries an unbounded empirical
discrepancy band until a caller supplies separate calibration evidence.


## Tensile activation at changed-coefficient endpoints

`TensileInitiator::activation_at` and `activation_at_controlled` observe
the declared Rankine criterion at one explicitly supplied elastic strain.
Thermal consumers pass held MECHANICAL strain under the new coefficients.
The result retains actual strain/stress, original elastic owner and strength,
principal plane/stresses, multiplicity and spectral residual. It carries no
strain fraction, onset time, bisection receipt or history-seeding operation.
The existing affine `between` detector still refuses an activated start.

Ordinary subcritical observations use the original allocation-free guarded
threshold screen; loaded or unresolved observations use one native 3x3
spectrum. Both APIs share the same objective principal-plane selection.
Repeated principal tension still requires a physical direction hint, and
conflicting screen/spectrum decisions refuse.

`ResolvedTensileActivation` preserves the complete evaluated source bundle.
Its `resolve_band` prepares only the observed plane and receiving width,
retaining every coefficient, control and additional receipt. The consuming
solid must preserve the OLD accepted source/state, admit conservative handoff,
and replay coefficient transfer through the original fracture owner. No
unseen temperature path, current damaged history or experimental calibration
is inferred. Four native regression definitions cover fixed-strain weakening,
threshold equality, rotation/scaling, physical hints, source continuity,
geometry preparation, cancellation and unchanged affine-path refusal.
Focused formatting and source/API/physics review completed; native build and
test execution remain unavailable without the required runners and pinned
constellation checkouts.

## Tensile crack-plane initiation and elastic handoff

`brittle::initiation::TensileInitiator` detects a maximum-principal-stress
tensile crossing on an explicitly affine, previously elastic strain increment.
It borrows the original isotropic elastic owner and uses the explicit tensile
strength. Fixed-size Gershgorin/positive-definiteness screens keep ordinary
elastic queries allocation-free; threshold uncertainty and the event plane use
the existing checked native Jacobi eigensolver. Strain domain, stress/fraction
resolution, bisection count and cancellation are explicit.

Repeated principal tension has no unique plane. An explicit material-direction
hint is projected into its principal eigenspace; missing or perpendicular hints
refuse. The event reports its numerical fraction bracket and spectral residual.
These are numerical resolution observations, not interval certificates.

`TensileOnset::seed` hands off at the original increment start, retaining its
elastic storage, then callers replay the actual increment through the original
fixed-plane crack band. Exact material/strength agreement, matching plane,
both endpoint strain domains, zero damage/dissipation, no prior frozen mixed
mode, and stress/storage continuity are required. Virtual preload work is not
charged again. No existing damaged history can be supplied to this initializer.

G0/G1/G3/G4 coverage in `tests/brittle_initiation.rs` includes analytic and
nonproportional onset, native fracture after conservative handoff, compression,
rotation/scaling, repeated principal stresses and bounded refusal/cancellation.
This tensile criterion does not select all-plane mixed-mode or crushing
failure, evolve crack paths, integrate FE fracture events or change topology.
It adds no experimental strength calibration.

## Brittle thermal weakening with retained fracture history

`brittle::temperature` adds an explicit fixed-strain coefficient/history
transfer to the original crack-band law. Young's modulus, normal/shear
strengths, fracture energies and the BK exponent may change. The native
cohesive temperature owner retains the physical onset mode, irreversible
damage and already dissipated fracture energy; cooling cannot rebond a crack.
The in-plane elastic remainder participates in the signed parameter-work
account. Publication requires parameter work minus total storage change minus
new fracture release to meet the caller's volumetric energy tolerance.

The plane, width and Poisson ratio stay fixed because they define the modal
kinematics. Source refresh preserves the prepared projection bit for bit,
including oblique normals. Strain must fit both owners' domains. Mechanical
work is zero during this endpoint transfer; parameter work is not heat by
definition. Thermal expansion, entropy, a temperature-path integration,
crushing and thermal fatigue remain outside this model.

`ResolvedBrittleCrackBand::at_temperature` refreshes the complete original
material bundle at one positive absolute-temperature coordinate using the same
card and selected claims. All other coordinates, source conventions, geometry
and controls stay fixed; additional properties keep their source identities.
`transfer_temperature` checks both original state ownership and the thermal
successor before moving history. A raw law substitution still refuses.
Refreshing a density observation does not modify mass or inertia.

G0/G1/G3 coverage in `tests/brittle_temperature.rs` exercises full elastic
parameter work, thermal failure, nonhealing on cooling/toughening, retained
fracture loss, compressive reclosure, rotated mixed-mode histories, source
refresh and oblique-plane identity. Numerical fixtures are not measured
material qualification.

## Source indentation hardness for contact and wear

`ResolvedDepthHardness::hardness_integral_between` integrates the original
piecewise-affine profile over exposed material depth in Pa m.
`advance_by_integral` inverts that monotone integral for a prescribed Archard
measure `k * p * ds`, crossing original source knots without freezing hardness.
It returns the retained depth increment and the original profile's interval
mean for the existing wear-volume owner. The mean is a quadrature value, not a
new material claim. Local coordinates and compensated accumulation preserve
small increments late in a wear history; the inverse scales hardness before
solving its quadratic. Runtime work allocates nothing and visits only crossed
source segments after binary search. Unsupported removal and lost arithmetic
or depth resolution refuse; source limits never extrapolate. Pressure feedback
and mechanical work remain with contact. Synthetic affine, cross-knot,
partition, extreme-scale and late-small-increment references exercise this
bounded integration; material calibration remains separate.

`hardness::plastic_history::ResolvedJ2Hardness::hardness_pa_at_accepted` observes
an accepted native J2 point without a second return map. Original mechanical
source identity, fixed conditions, hardening convention and law must match;
native small-strain/history/yield admission and prepared hardness support remain
binding. Scalar observation allocates nothing. `observe_accepted` retains a
small numeric endpoint for the existing typed-source export; it copies no source
receipts and adds no work, storage or heat. The spatial solid observer composes
this same route for each original element. Numerical admission is not evidence
that an arbitrary supplied tensor came from a converged physical trajectory.

`hardness::cyclic_history::ResolvedCyclicHardness` extends the same prepared
original H(alpha) source to native source-bound Voce, single-branch Chaboche and
ordered multi-branch Chaboche laws at fixed external conditions. The complete
mechanical source, nonlinear coefficients and named branch order remain binding,
even when two branches happen to have equal coefficients. Every branch passes
the native history bounds and the accepted endpoint passes its actual nonlinear
yield surface before hardness is sampled. A cyclic state is never projected
onto a bare linear-J2 yield surface. Scalar observation allocates nothing;
optional typed export also enforces the originally prepared support interval.
Dynamic backstress recovery does not erase accumulated strain or supply
annealing kinetics, and an independent hardness change adds no work or heat.
`tests/hardness_cyclic_history.rs` exercises native loading, unloading, reversal,
nonproportional recovery, full typed-source retention and atomic source/history
refusals. These synthetic checks do not establish a fatigue-life relation.

`hardness::history_field::AcceptedHardnessField` exports a compact ordered
alpha/H snapshot from actual accepted J2, Voce or complete Chaboche tuples.
One shared original mechanical/hardness descriptor retains fixed conditions,
the exact test, prepared support and native family/branch order. Every actual
point law must match that source. Construction uses the existing accepted
checks once per point, without a return map, database query or per-point source
cloning; late failure or a caller checkpoint publishes no partial field.
Source identity is independent of point values, while snapshot identity also
binds their declared ordering. Checked point lookup supplies explicit consumer
maps without interpolation. This frozen external/preconditioning property
transfer does not carry residual stress, plastic displacement, geometry or
elastic compliance, and does not attest a mesh or loading trajectory. Contact
consumers retain their own homogeneous half-space geometry and mapping; contact
depth/recession cannot replace alpha or duplicate bulk indentation work.

`hardness::depth::ResolvedDepthHardness::at_temperature` re-prepares the exact
original H(depth) test at a new prescribed canonical absolute T. The immutable
card is required even for a no-op; the selected claim, complete test, all other
typed coordinates and maximum original-material depth remain pinned. The same
preparation owner rechecks continuous support and positive hardness across the
entire depth interval, including source-cell boundaries, holes and nonpositive
interiors despite positive endpoints. Runtime depth evaluation and exact depth
integration retain the existing prepared numerical kernels. This supports
original H(T, depth) tables without annealing kinetics, heat, expansion or any
change to removed material or indentation history. It admits the new depth
line only, not the unobserved intervening temperature path. The existing
`hardness_depth_integral` target covers independent hot/cold affine integrals,
typed-coordinate retention, competing claims, pinned cooling, no-op identity and
whole-profile/card/temperature refusals using synthetic data.

`hardness::thermal_history::ThermoplasticHardnessSource` observes an independent
original H(T, accumulated-equivalent-plastic-strain) test at accepted native
`CaloricPoint` states. J2, Voce/viscous and complete N-branch Chaboche histories
retain their original return maps, thermal coefficient transport and work/heat
accounts. The adapter returns the unchanged native feedback report and typed
hardness together, refusing foreign mechanical/caloric sources, changed tests,
missing endpoint coverage, unowned conditions and nonpositive hardness without
advancing the input point. Hardness is sampled at old/accepted endpoints; this
does not certify support or integrate hardness between those states, and supplies
no annealing kinetics, yield-to-hardness conversion or additional heat. Focused
`hardness_thermal_history` tests compare the complete native histories and energy
reports under heating, reversal and viscous cyclic loading, and exercise source,
coverage and cancellation refusals. These are synthetic regressions.

`try_from_native_state` / `at_native_state` share the same source admission with
spatial physical owners. The actual caloric chart, original mechanical card,
complete selected mechanical claims, hardening convention and resolved thermal
context must agree before observing the complete native history. This borrowed
ingress repeats no constitutive solve or mechanical coefficient refresh; the
calling physical owner remains responsible for accepted equilibrium and history.

`hardness::IndentationHardnessStatePoint` selects a positive static
pressure-valued `IndentationHardness` claim from an original material bundle.
It retains the exact apparatus, load/dwell, protocol, specimen and source
selection through the existing complete property key. Refresh uses the same
card, query conventions and pinned claims, preserving any additional coupled
properties and refusing absent coverage or nonpositive hardness atomically.

The shared borrowed `indentation_hardness_property` admission is used by dry
sliding and Mindlin fretting wear. Existing explicitly selected dimension-only
Pa values remain caller-declared hardness estimates; they do not become typed
test-qualified measurements. Typed fluid pressure/stress, empirical hardness
numbers, and tensor components refuse. Consumer owners still admit the source
coordinates and fixed/thermal operating conditions. This adapter evolves no
plasticity or wear and infers no HRC/Shore/Vickers conversion, yield stress,
fracture toughness, hardening law, or empirical validation.

Focused G0/G3 tests exercise distinct tests under the same property name,
temperature refresh of original claims, full-bundle retention, exact context,
changed-card/coverage refusal, and the declared-versus-typed admission boundary.

`hardness::depth::ResolvedDepthHardness` prepares an exact original indentation
test along the dimensional `exposed-material-depth` coordinate, from the
reference surface to an explicit supported depth. Every other coordinate is a
fixed external condition. Native source-line preparation checks continuous
support and positive interior breakpoints, refusing holes, discrete-only data
and extrapolation. The complete original support receipts remain available;
ordinary positive piecewise-affine evaluation uses binary search without source
queries or allocations. Typed endpoint export still requires the original card
and pinned claim. Depth means original material exposed by actual removal; it is
not indenter penetration, relative plastic indentation, inferred case depth,
layered elasticity or a coating-failure law. Synthetic contact/wear tests check
source holes, exact tests, graded-profile refinement and zero-depth endpoints.

### Joint original depth and accepted plastic history

`ResolvedDepthHardness::{from_j2_accepted,from_cyclic_accepted}` prepare an
original H(z, alpha) profile after the native J2, Voce or Chaboche owner validates
the complete accepted tensor history and instantaneous yield surface. The exact
indentation test, original material card, fixed external context, native law and
complete cyclic branch identity remain bound. The private
`AcceptedDepthHardnessHistory` receipt exposes provenance and actual
dimensionless alpha without cloning constitutive histories into wear proposals.
Ordinary construction still cannot admit an arbitrary caller-supplied alpha.

Native history supplies alpha and actual removal supplies dimensional exposed
material depth. Continuous source coverage, positive interior knots, prepared
binary-search evaluation and the existing exact piecewise-affine integral and
stable inverse are reused. Preparation queries the source once; wear stepping
does not rerun plasticity or material resolution. Exact no-op temperature
refresh is admitted, but changed temperature requires newly validated native
history and a newly prepared profile.

This is a frozen, externally preconditioned material column: the caller must
establish that its admitted alpha represents the declared depth interval.
Surface indentation does not advance alpha, and no residual stress, bulk
plastic work, thermal history or hardness-to-yield conversion is transferred.
Native load/unload/reversal, full cyclic-branch admission, exact integral
inversion, selected-test identity and support-hole/interior-positivity refusal
regressions specify this numerical contract. They do not establish empirical
material qualification.

## Bulk brittle crack bands through the original cohesive history

`brittle::BrittleCrackBand` composes isotropic small-strain elasticity with the
original mixed-mode `BilinearCohesive` owner on an explicitly supplied fixed
plane and band width `h`. A positive plane-stress remainder retains in-plane
elasticity. Three work-conjugate band coordinates carry normal opening and
the two shears; their intact stiffnesses are `(lambda+2*mu)/h` and `mu/h`.
The normal coordinate includes the elastic Poisson correction. Intact stress
and tangent recover the original three-dimensional elastic law, including
auxetic backbones. No spectral tensile-strain split or inferred crack direction
is used. The normal coordinate is not a measured geometric face gap.

The original cohesive law owns quadratic initiation, frozen-onset BK mode
mix, irreversible damage, unloading, softening tangents and analytic path
work. Complete failure removes tensile and shear transfer across the chosen
plane while retaining in-plane elasticity and compressive normal response.
Its fracture work per volume is the original work per area divided by `h`;
pure-mode complete failure spends exactly `Gc/h`. Construction refuses bands
whose elastic initiation energy leaves no positive softening interval. Width
scaling regularizes this declared band experiment; it is not a general mesh-
or orientation-independent crack-growth theorem. Positive local softening
does not rule out structural snap-back; the original solver retains its
nonconvergence/refusal and displacement-control boundaries.

`small_strain_law` borrows the immutable prepared owner for existing generic
solid solvers. Each state retains its original law, strain and cohesive history;
rejected trials never publish history. The checked `trial` returns stress,
the classical tangent where it exists, a selected branch tangent, storage,
fracture work and the energy defect. The legacy infallible trait signals
refusals with nonfinite responses. Mixed-mode initiation can have a
nonsymmetric tangent and advertises that fact to the native linear solver.
Prepared mappings and point updates use fixed-size arrays without per-point
law/table clones or heap allocation.

`brittle::source` retains a complete `ResolvedMaterialStatePoint`, including
other physics properties, while requiring explicit brittle tensile/shear
strength, mode-I/mode-II fracture energy and BK exponent alongside E and nu.
Only fixed supported external conditions are admitted. No modulus-of-rupture,
hardness, yield stress, surface energy or material-name default supplies an
absent fracture datum. A newly resolved temperature cannot replace the law
under an existing history. Source uncertainty and selection receipts survive;
they do not establish a statistical breaking probability. Without calibration,
the model card retains an unbounded model-form discrepancy.

This is a fixed-plane, small-strain constitutive estimate, without automatic
crack nucleation, multiple crack families, evolving crack geometry, fragment
creation, finite rotations, frictional recontact, thermal-history transfer or
experimental calibration. Existing separated-interface/contact owners retain
those responsibilities. G0/G1/G3 regressions specify intact elasticity, pure-
mode fracture energy, width scaling, unloading/reclosure, rotation, consistent
tangents and source/owner/domain refusals. The broader MR63 fragment-transfer
acceptance remains open.

## Source-bound isotropic tensile limit

`strength::IsotropicTensileLimit` binds one explicitly named pressure-valued
strength claim at its original material/environment state. The caller selects
the isotropic maximum-principal-stress criterion; it is never inferred from a
material name or substituted for a ductile yield law. Assessment consumes
tensile-positive tensor-shear Cauchy stress `[xx,yy,zz,xy,yz,zx]` and reports
`max(0, sigma_max)/strength`, using the shared scaled, allocation-free 3x3
Jacobi kernel. Thus hydrostatic tension is retained and compression does not
consume tensile capacity. Missing, nonpositive, mistyped or dynamically
conditioned strength refuses; fixed external conditions and source-pinned
state refresh use the original resolver. Nonfinite stress and unrepresentable
results refuse without partially updating any state.

This is a bounded local criterion observation, with floating-point eigenvalue
estimates. It is not an allowable-stress conversion, uncertainty bound,
anisotropic failure surface or fracture evolution. Engineering ultimate
strength is not automatically a true-stress limit after finite strain/necking.
Toughness, damage, crack energy, plastic redistribution, fatigue and
post-failure load transfer require their own laws. G0/G3 regressions in
`tests/tensile_strength.rs` cover tensile/compressive states, hydrostatic
tension, tensor shear, rotations/repeated roots, extreme scaling, source
admission and source-preserving temperature changes.

## Paris-region fatigue crack growth

`fatigue::crack_growth` adds an immutable, explicitly scoped Paris calibration
and an evolving `ParisCrackGrowth` state. The normalized law is
`da/dN = reference_rate * (delta_K/reference_delta_K)^m`, with length in metres,
full cycles, and the complete stress-intensity range in Pa sqrt(m). The
reference rate avoids the exponent-dependent units of the conventional C.
The original `FatigueScope` retains material, process, surface/geometry,
environment, load basis, temperature and caller-declared source reference.
These declarations do not confer material-database or measured-data authority.

`ParisLaw::new_walker` explicitly selects calibrated tensile load-ratio
dependence. Its reference rate belongs to both the reference full Delta-K and
`WalkerParameters::reference_load_ratio`. The normalized multiplier is
`((1-R_ref)/(1-R))^[m*(1-gamma)]`, following AFGROW DTD Handbook equation
5.1.3; gamma is the inner Walker exponent and m is the existing growth
exponent. The admitted gamma branch is `[0,1]`: gamma=1 retains plain Paris
arithmetic, and gamma=0 gives Kmax dependence for supported positive ranges.
The reference R and reference Kmax must themselves satisfy the original fit.
Every runtime check still uses actual R, actual full Delta-K and actual Kmax.
Near-equal ratios use their direct difference before logarithmic evaluation.
At fixed stress extrema this correction changes the rate multiplier while the
existing constant/table geometry power and analytical integrations remain
unchanged. No compression, closure, overload retardation, near-threshold law
or high-K fracture acceleration is inferred by selecting Walker.

`crack_growth::source::resolve_paris` resolves four original material-card
claims: squared reference full Delta-K, growth in metres per full cycle,
the Delta-K exponent, and squared exclusive maximum-K validity cutoff.
Squared intensities use Pa² m (`Dims([-1, 2, -4, 0, 0, 0])`), with an explicit
square-root conversion at the adapter; K is never mislabeled as pressure.
All four claims retain one provenance, observation population and identical
complete typed support domain, including semantic absolute temperature,
load ratio and squared Delta-K. Original scalars and temperature-only
interpolation are admitted; range/ratio-dependent fits and 2D surface boxes
cannot supply unsupported joint support. Original claims, uncertainty,
selections, answers and receipts remain in the immutable law. Native consumers
bind the actual original material card/state and every fixed physical query
coordinate; only calibration R and squared Delta-K are excluded from operating
context matching. No missing growth data, toughness, threshold or measured
authority is inferred from unrelated material properties.

`crack_growth::source::resolve_walker` resolves all six original parameters
from distinct `walker-*` properties, including the reference load ratio and
gamma. All six retain one provenance, observation population and identical
complete typed validity domain; only scalar or temperature-only parameter
interpolation is admitted. The four common normalized parameter roles and
the two extra original claims, answers and selections remain inspectable.
Neither resolver substitutes the other model's property keys, and querying
a plain Paris fit cannot silently enable Walker. Reference R is an original
fit parameter, not the operating R chosen in the source query. Native binding
continues to require the actual original material card/state and fixed typed
operating context for either explicitly selected model.

`CrackGeometry` retains either constant Y or an immutable, bounded original
table of Y versus crack length. Both use `K=Y(a)*sigma*sqrt(pi*a)`. The constant
case retains its analytical O(1) path, including the exponential m=2 limit.
`TabulatedCrackGeometry` retains the original source reference, crack-length
convention, nominal-stress basis and all points without sorting or repairing
them. Within each log-log interval, `Y=Y0*(a/a0)^s`, so the Paris growth rate is
proportional to `a^[m*(s+1/2)]` and integrates analytically. Increasing,
constant and decreasing intensities are supported. Work is bounded by crossed
intervals, independent of cycle count; point caps bound allocation and interval
work, and cancellation is polled at construction and interval boundaries.
There is no extrapolation. Geometry-table interpolation and adequacy of its
resolution remain declared model approximations, not certified error bounds.
Ordered spectra preserve their input order, and clones retain the complete
accepted state with shared immutable geometry. This memoryless Paris law
supplies no overload/retardation sequence interaction.
It must not consume rainflow completion order as a chronological spectrum.

`FiniteWidthCentreCrack` supplies the elastic secant geometry
`Y=1/sqrt(cos(pi*a/W))` for a centred through crack under remote uniform
gross-section tension. W is the full plate width and a is the physical
half-crack length. The supported formula range is `0<a<=0.4W`, following
AFGROW DTD Handbook Table 11.3.8 and NASA TN D-5140 p. 6; this geometry does not
calibrate the Paris law or establish plate aspect-ratio/grip applicability.
The helper retains width, references and stress/length conventions in tables
sampled at exactly the caller's ordered lengths. Deterministic direct factor
evaluation, bounded construction and cancellation share the existing owner.
No knots are inserted, and table resolution remains explicit: refinement tests
compare life against independent integration of the unsampled secant formula.
There is no finite-height, bending, plastic-zone or 3D crack-front correction.

Alternatively, `tabulate_to_tolerance` generates a table over a declared
half-length interval with a caller-supplied relative-Y interpolation tolerance
and point cap. In logarithmic length, the secant model's positive increasing
curvature bounds the exact linear-interpolation error on each interval.
Bounded logarithmic refinement uses that analytical criterion, concentrates
points where geometry changes faster, and preserves the exact requested
endpoints. This criterion controls geometry interpolation in exact arithmetic;
ordinary floating-point evaluation is not an interval-arithmetic certificate,
and neither the handbook approximation nor material uncertainty is reduced by
refining the table. Cancellation, unrepresentable refinement and exhausted
point budgets refuse without publishing a partial geometry.

Positive-range cycles require the supplied tensile R interval and Paris-region
range bounds. Every sample can be checked against the current K limit before
cycle closure. Growth checks every crossed interval's endpoints, including
interior intensity peaks or valleys of the complete geometry table, against
the declared fitted range and exclusive upper-K limits. Advancing beyond the
last geometry point reports the supported length and unconsumed exposure;
it publishes no partial extension. Forecasting and advancement use the same
compensated cumulative interval counts, so forecasts reach original knots
without tolerance-based clipping. Nonrepresentable interval counts refuse.
Zero-range exposure adds no fatigue growth but still checks static K.
The lower fit bound is not a threshold: unsupported low loads refuse. Blocks,
spectra, cancellation and failed endpoint admission publish no partial growth.

`ParisCrackGrowth::cycles_to_length` provides a read-only analytical forecast
to a caller-declared length strictly above the accepted crack length, under
fixed tensile stress extrema. It returns fractional additional full cycles and
both endpoint stress intensities, with constant work for constant Y and bounded
interval traversal for tabulated Y. `cycles_to_length_with_cancel` exposes the
same cancellation-aware calculation to native observers. The inverse
retains the m=2 limit and tiny relative extensions without subtracting nearly
equal powers. Geometry, R, Delta-K and exclusive maximum-K support must hold
over the complete interval. One block and the predicted cycles must fit the
remaining budgets, with a distinct finite cumulative exposure; no budget or
state is consumed. Zero range, unsupported targets and nonrepresentable inverse
arithmetic refuse, including rounding to the finite-cycle singularity. The
selected target is not an inferred critical crack size, service-life estimate
or failure certificate.

G0/G1/G3 regressions compare independent m=1/2/4 integrals, the m=2 limiting
case, range/full-cycle conventions, partition/restart behavior, source-scope
retention, static and evolving-domain refusals, budgets and rollback. The
formula follows NASA-HDBK-5010 Volume 2 Revision A §7.3 and AFGROW DTD Handbook
§5.1.2. These synthetic checks verify implementation, not component life.
This lane supplies no inferred crack initiation, closure/compression correction,
mixed-mode growth, thermal transfer, resolved mesh
crack, stiffness degradation, heat, toughness measurement or rupture guarantee.

## Objective small-strain magnitude

`tensor::frobenius_norm` evaluates the physical symmetric-tensor norm in the
native `[xx, yy, zz, xy, yz, zx]` tensor-shear convention. It counts shear twice,
scales before squaring, preserves NaN and distinguishes total strain from a
deviatoric J2 measure. `evidence_stress` uses this rotation-invariant magnitude
for the generic `strain-magnitude` coordinate; a rotation cannot move the same
strain into a scalar calibration domain by spreading it across components.
Nonfinite strain or stress always gives out-of-domain model evidence.

Native Norton creep uses the same norm for its existing physical strain and
stress bounds. Full-tensor, frozen thermoelastic and crystallization adapters
retain their original componentwise Mandel admission in checked kernels. Their
generic evidence card uses the inscribed Frobenius ball, so it is conservative
and does not replace the source's anisotropic domain or validate a material.
G0/G3 checks in `tensor` and `tests/conformance.rs` exercise physical shear,
hydrostatic strain, rotations, extreme representable scales and nonfinite input.

Isotropic and orthotropic construction require finite positive strain bounds
and representable stiffness. A positive isotropic Young's modulus cannot
silently become zero shear stiffness through underflow. The existing compensated
spherical/deviatoric isotropic stress owner, including its auxetic response,
is retained; the constructor gate does not replace it with a dense matrix
evaluation. Focused G0/G1 regressions cover domain and stiffness admission,
nearly incompressible response and all six columns of the original consistent
tangent. These remain small-strain laws; callers own strain admission, and
no finite-rotation response or material-accuracy claim is implied.

## Integrated bulk conduction with an evolving source contact

`ImplicitLayerLink::with_contact` composes the original
`ImplicitContactLaw` with both source-integrated material half-cells. It
requires zero authored contact resistance, one contact owner, finite source
bounds wholly inside the supported pure-phase layer domains and the same
ordered caloric objects. Every trial checks those bounds, evaluates the original
contact once, and uses the same
Kirchhoff solve for power and its full derivative:
`dQ/dh = (dQ/dT)(dT/dh) + (dQ/dR'')(dR''/dh)`.
The analytic areal-resistance derivative uses the accepted interface flux and
differential conductance without a second solve. Unrepresentable nonzero
tangents refuse. Contact-state admission is shared with native thermal links.

The contact is prepared on the shared admissible branch; wider contact domains
refuse instead of heuristically resampling source tangents. At a contact upper
bound inside the bulk chart, the bulk tangent uses its exact left caloric
segment, keeping both chain-rule contributions on the same admitted side.
The original thermal chain/graph retains energy, entropy, bounded work and
transactional publication. Source response coordinates and fixed external
pressure/composition keep their original meaning; no interface-temperature law
or evolving pressure is inferred. This optional composition does not alter the
fixed-resistance `PreparedLayerConduction::resolve` adapter. Pure-phase,
quasi-steady half-cell assumptions and independent time/spatial refinement
remain. Regressions check complete-flux finite differences, domain/ownership
refusals and rollback, and actual `PreparedMaterialContact` source selection
through a manufactured 30 J exchange in both native implicit solvers.


## Source-backed ionic solubility and osmotic particles

`osmosis::solubility::electrolyte::ElectrolyteSolubilitySolution` retains one
original `HydrolysisPrecipitation` or `SaltPrecipitation` model, its selected
formation-Gibbs sources, complete analytical cation/family stocks, declared
spectator stocks and accepted phase/pH result. Its shared ideal aqueous view
counts the actual H+, OH-, acid-family members and spectators exactly once.
Pure solids and a second nominal dissociation multiplier contribute no particles.
The original analytical totals remain distinct from dissolved-only reports.

`SolventActivityBasis` explicitly binds solvent molar volume and reference
amount per physical volume. Each solvent-volume trial scales the common
amount at unit activity from that declared density and re-equilibrates the
original chemistry. Unit-molality sources require the caller's solvent-mass
conversion; no density, standard-state conversion or charge is inferred from
a name. Volume trials share immutable analytical-stock and species-axis
storage. The original bounded chemistry solver still owns source validation,
component/charge balance, work limits and cancellation. This adapter does not
claim allocation-free ionic volume trials.

`with_analytical_totals` permits finite stock updates on the same source,
family and spectator layout, including material held in the selected solid.
Failed construction, volume changes, dosing and cancellation publish no partial
state. The native tests use source-qualified manufactured 1:2 salt and
hydrolysis equilibria to check particle counts, reversible phase changes,
preserved sources and inventory/layout refusals.

The original fixed-temperature, fixed-pressure ideal chemistry uses unit
water activity. Osmotic pressure and water activity are the leading dilute
particle approximation, not a concentrated-electrolyte activity correction.
Actual H+/OH- inventories retain the small pure-water autoionization reference
term (2*sqrt(Kw) times the declared reference-amount density); no reference
osmolarity subtraction is made. Empty analytical stocks therefore do not
supply an exactly normalized pure-water solvent activity. Chemical reaction
water, crystal/solute partial volumes, phase-change heat,
extra aqueous complexes, competing solids and kinetic rates are not supplied
by this equilibrium adapter.

### Finite-rate ionic material states

`electrolyte::kinetics::KineticElectrolyteSolubilitySolution` retains one original
`MineralKineticsV1` and its actual mineral/aqueous state. Constructors require
explicit dissolved cation, dissolved analytical-family and solid inventories.
The original source-bound chemistry owns fast aqueous speciation and finite
mineral exchange. The material adapter never initializes an undeclared
equilibrium partition or substitutes total family concentration for free anion
activity.

Solvent-volume changes preserve the actual crystal and dissolved component
stocks, scale the common activity basis with the original `SolventActivityBasis`,
and re-solve only aqueous pH/protonation. Subsequent time advancement then
dissolves or grows the solid under the explicitly declared surface rate and
growth policy. This makes the delay after dilution or solvent removal visible
instead of immediately imposing mineral equilibrium. Complete analytical totals
remain distinct from the current dissolved-only osmotic particle view.

`with_dissolved_analytical_totals` replaces explicitly declared aqueous cation,
family and spectator amounts on their existing identity/charge layout. The
actual solid remains fixed while fast aqueous chemistry responds to dosing or
removal. The resulting complete totals include that unchanged solid with its
original stoichiometric coefficients. Subsequent mineral transfer uses the same
rate source; missing counterions or undeclared added species are not inferred.

Both equilibrium and kinetic material paths use one shared publisher for H+,
OH-, each analytical-family member and each spectator exactly once. The solid
and solvent do not enter that dissolved particle axis. Accepted updates retain
the source and species identities; failures or cancellation publish no partial
material endpoint. The native kinetic outcome remains available for inspection
of actual saturation, balances and numerical work. The same fixed-temperature,
ideal-dilute and fixed-solvent chemical limitations above apply; the adapter does
not add ionic heat transport, complexes or nonideal activity corrections.

## Absolute thermo-optic response for refraction

`thermo_optic::SchottThermoOpticLaw` retains explicitly supplied SCHOTT TIE-19
coefficients, an absolute reference temperature, a nonzero source association,
and positive finite temperature/vacuum-wavelength domains. Active resonances
anywhere in the inclusive source band refuse. The source association does not
establish measurement authority or qualify coefficients for a material.

`prepare` evaluates the spectral coefficients once from the original absolute
phase index and optional original group index. The group correction retains
both the wavelength derivative of the reference-index prefactor and the
resonance term. Unknown reference group index remains unknown. Each subsequent
allocation-free cubic evaluation returns phase/group indices and exact first
and second temperature derivatives, always relative to the original reference
state. The reference-temperature indices retain their bits; unsupported states,
nonpositive indices and unrepresentable responses refuse without clamping.

The law supplies no absorption, thermal expansion, stress birefringence, air
correction or inferred material data. Controls include independent
high-precision phase/group values, wavelength and temperature derivatives,
reference-state recovery, signed coefficients, pole/domain refusal and identity.

## Source-bound hyperelastic membranes

`ResolvedSheet::with_hyperelastic_membrane(card, model_index, settings)` explicitly
selects a native Neo-Hookean, Mooney--Rivlin, Ogden or HGO model from the exact
original material card at the same retained typed query. Independently sourced
reference areal mass and thickness, extra selected properties, receipts and
model provenance survive the shared `MembraneLaw::Hyperelastic` carrier.
Mechanical/frequency-conditioned reference properties, competing selected
membrane capabilities, missing model/context and native domain failures refuse.
The supplied thickness and areal mass belong to the original mesh reference.
Without an explicit natural-reference selection that reference is stress-free;
bending and evolving temperature/history are not inferred.

For engineering Green strain `[E11,E22,2E12]`, an objective Cholesky
representative supplies the original native plane-stress owner, including all
three transverse traction conditions and its Schur-complement tangent.
Reference thickness scales energy, work-conjugate second-Piola resultants and
the full chain-rule tangent; the stress times second kinematic derivative is
retained. Existing `fs-solid::membrane::MembranePatch` therefore uses these
models through its original 3-D kinematics, geometric stiffness, dead/pressure
equilibrium and reference mass, without another spatial solver. HGO activation
and unresolved principal-stretch boundaries refuse two-sided vibration
admission. Native local stable-branch limits remain; no wrinkling or global
stability certificate is added.

`with_incompressible_hyperelastic_membrane` explicitly selects elastic `Je=1`
through the original `hyper::incompressible_plane_stress` owner. Reference
thickness is unchanged. With the original natural reference, current normal
stretch is the reciprocal in-plane area ratio.
Both transverse shears relax against the native energy remainder and the full
constrained tangent flows through the same Green-strain chain rule. No bulk
penalty is evaluated or replaced. The original model card, source bundle,
reference mass and stretch bounds remain authoritative; exact incompressibility
is an imposed modeling choice, not new experimental source evidence. The
retained `volumetric_treatment` distinguishes this choice from the original
finite-compliance selection, including in membrane equality. The existing
`MembranePatch` equilibrium, pressure and vibration consumers use either choice
through the same `MembraneLaw` without another spatial solver.

`with_hyperelastic_reference_stretches` explicitly supplies positive natural
length ratios in warp/weft/normal order. The fixed diagonal G maps the mesh
reference to the native model's stress-free frame: Fe=F G^-1 and
W_ref=det(G) W_native(Fe). Exact Piola and tangent pullbacks precede the same
Green-strain chain rule; no second local Newton solve is introduced. The
original thickness and areal mass remain fixed. Native stretch and fibre
activation bounds apply to Fe. Exact incompressibility therefore constrains
elastic volume while total det(F)=det(G), permitting a changed natural volume.
The default reference and explicit identity G retain original response bits.

`with_hyperelastic_mean_thermal_expansion` binds three finite mean-expansion
endpoint sources and explicit reference/endpoint temperature-axis roles. Each
length ratio is 1+the source's engineering strain. The exact card, typed fixed
context and both temperatures must agree across sources; original mass,
thickness and initial model belong to the reconstructed reference context.
The same selected model is admitted at the endpoint context. All three source
objects and axis roles survive in the material carrier. Small-strain alpha
integrals are not reinterpreted as finite reference maps. Reference selection
is explicit and occurs once; no temperature evolution, reciprocal heat or
work balance for changing G is supplied by this fixed-endpoint preparation.

`sheet_hyperelastic` checks independent compressible Neo-Hookean reduction,
all native model tangents, tilted-fibre transverse relaxation and source/reference
refusals. Its exact-incompressibility controls add an independent analytic
Neo-Hookean Green-strain potential, all four laws' work/tangent derivatives,
bitwise independence of unused bulk parameters, tilted-fibre traction closure,
and retained source/mass/domain controls. `hyperelastic_membrane` checks 3-D
force/tangent objectivity, independent mass, and finite-deformation
dead-load/pressure equilibrium against an independent pyramid oracle.

Natural-reference regressions add independent finite/exact Neo-Hookean
energy/resultant/tangent formulas, all-model derivatives, elastic-domain and
HGO-activation admission, typed thermal sources and unchanged mass. Bounded
original-source numerical checks execute the pullback and constitutive chain;
full typed-source and native consumer tests still require the complete workspace.

## Fixed-mass PCM storage in a sensible solid matrix

`phase::composite::PcmComposite` borrows one sensible all-solid matrix branch
and one complete solid/liquid PCM chart. Its declared PCM mass fraction remains
fixed. The constituents share temperature, retain their original source receipts
and density interpolation laws, and have additive enthalpy and specific volume:
`1/rho = (1-w)/rho_matrix + w/rho_pcm`. Actual composite liquid mass fraction is
`w*f_pcm`; complete PCM melting leaves the matrix solid. PCM volume fraction is
reported separately. A supplied total mass remains fixed as the returned free
volume changes; the law does not update reference geometry or thermal inventory.

Preparation merges every original caloric breakpoint, preserving the admitted
piecewise-affine constituent `T(h)` exactly in real arithmetic, including latent
plateaus. The explicit composite lower enthalpy datum is independent of both
constituent datums. Cached slope and logarithmic, allocation-free state queries
feed the existing `ThermalLaw`, thermal timeline and backward-Euler residual;
entropy remains the mass-weighted original constituent integrals. Density is
evaluated on the original constituents, never fitted to the compiled mixture
knots. The combined original-knot preparation ceiling is 8192.

The caller supplies compatible source pressure/work conditions and a nonreacting,
local-equilibrium, additive-volume model with no heat of mixing. Unsupported
temperatures, nonpositive fractions and unresolved caloric coordinates refuse.
A matrix with an ambiguous solid plateau is excluded; two melting constituents
would require a separate phase-allocation closure. No confinement, evolving
geometry or extra pressure-work solver is supplied. Manufactured controls cover
latent energy, mass/volume fractions, original density laws, merged sensible
breakpoints, entropy, independent datum shifts, and explicit/implicit native
heating across the latent plateau.

`thermal::conduction::composite::PreparedPcmConductivity` binds those exact
original constituent caloric objects and their prepared material conductivity
sources to the existing implicit conduction owner. The caller explicitly selects
parallel or series layers; actual constituent volume fractions determine the
directional conductivity, while mass fractions continue to determine storage.
Analytic tangents include both original conductivity and density changes,
including nonzero conductivity changes on an isothermal fusion plateau. Source
knots retain one-sided derivatives inside the admitted composite domain. The
shared adaptive observable reports total composite liquid mass fraction `w*f_pcm`,
not a normalized PCM-only fraction. State evaluation requires no allocation,
source lookup, new caloric chart, finite-difference query or extra latent source.
These layered homogenizations do not claim isotropic foam morphology, interface
resistance, percolation, melt convection or constituent thermal nonequilibrium.
Reference geometry and expansion work remain with their original owners.

## Source calorics, kinetics and shrinkage for crystallization

`phase::crystallization::CrystallizationCaloricPair` borrows two original pure
material charts over a declared common temperature interval. It rejects
coexistence plateaus, uncovered metastable temperatures, mismatched material
identities and nonpositive crystallization enthalpy gaps. The union of both
charts' knots establishes positivity of the represented phase enthalpy gap;
each original heat-capacity/density interpolation and source receipt survives.
Source-backed and explicitly authored inputs retain their respective authority.
The caller must supply compatible pressure/composition conditions and one
physical enthalpy reference. Temperature-dependent latent heat is the original
branch difference, not a new heater or duplicated apparent heat capacity.

`kinetics::PreparedCrystallizationRate` prepares the canonical typed material
source line for Nakamura K(T), in inverse seconds. It preserves every source
breakpoint and original receipt, rejects invalid units/rates and unsupported
conditions, and performs allocation-free binary-search interpolation during
stepping. Source Cp lower bounds accelerate the shared caloric inversion without
replacing the original charts. All parameters and caloric/kinetic source
identities enter the resulting immutable state identity.

Actual kinetic MASS crystallinity supplies additive constituent specific
volumes, bulk density and a separately reported crystalline VOLUME fraction.
`mechanics::CrystallizationLock` stores an explicitly selected stress-free
temperature, fraction and total strain. Its frozen `SmallStrainLaw` view uses
the same density-derived free-volume change and caller-supplied original
elasticity in native continuum solvers. Separate total, elastic and free-strain
bounds apply; prestress, mismatched sources and unsupported melting/unlocking
refuse. No universal solid-support threshold or stiffness/fraction law is
inferred. The mechanical response is one-way: elastic/parameter work is not
silently charged to the thermal crystallization account. Locking and shrinkage
carry no calibrated model-discrepancy bound.

Focused regressions exercise original source receipts/rates, distinct Cp and
latent feedback, absolute datum shifts, mass/volume conversion, cooling through
the shared thermal timeline and native free/constrained shrinkage. The runnable
`crystallization_cooling` example exchanges heat with a finite mold and prints
temperature, both crystal fractions, free volume and the energy balance. Its
inputs are explicitly manufactured; it is not a measured material prediction.

### Original properties on the accepted crystallization history

`phase::crystallization::properties::CrystallizationPropertySource` binds the
original material card, caloric pair, kinetic law and exact selected scalar
property claims. It updates the accepted temperature and explicitly named
crystalline mass and/or volume fraction coordinates. Absolute crystal mass
fraction is distinct from normalized conversion and from liquid fraction;
the parent may itself be an amorphous solid. Volume fraction uses the original
constituent densities. Other admitted pressure/composition conditions remain
fixed, and unowned evolving coordinates refuse.

Every evaluation validates the original kinetic state and refreshes the same
source claims as one atomic material bundle. Unsupported crystallinity or
temperature cannot select a different source, extrapolate or clamp. Thermal,
electrical and mechanical consumers can request their own typed scalar
requirements on that shared state. No conductivity or stiffness mixing rule is
invented. A separately queried density/Cp does not replace the original mass,
caloric or pressure-work owner, and changing elastic parameters does not by
itself close thermomechanical energy feedback.

Regressions cover distinct mass/volume coordinates, fixed conditions, original
claims, foreign-history refusal and finite-body contact against an independent
analytic heat balance with kinetic conductivity. Native triangle/quad mechanics
consume the same source bundle through the existing elastic and shrinkage laws.
The cooling example now preserves incoming crystal history and evaluates its
declared conductivity source within the existing implicit network residual.
These manufactured cases test numerical integration, not material calibration,
source support between endpoints, morphology or coupled mechanical work.

`phase::crystallization::source::CrystallizationSourcePair` lets these same
property, locking, shrinkage and unlocking adapters consume the original
`IsobaricCrystallizationPair`. Its exact Cp integrals and phase-specific free
expansion are borrowed directly; no sampled replacement chart is introduced.
The common mass-to-volume conversion retains the pair and both density-source
identities. Exact-Cp property queries must preserve the original temperature
axis and fixed pressure/composition coordinates. A calorics-only pair refuses
volume-dependent operations until both original phase densities are supplied.
Isobaric elastic-release transfers additionally check the actual source pressure
and retain the existing separate heat, internal-energy and pressure-work account.
The `crystallization_exact_cp` regressions exercise analytic recalescence, live
mass/volume properties, free/constrained shrinkage and conservative elastic
release on manufactured affine-Cp and free-expansion sources. Continuous
thermomechanical feedback and measured material validation remain separate work.

### Material properties of retained crystal populations

The original thermal `CrystallizationPopulation` equilibrates incoming masses
at conserved enthalpy while retaining their individual kinetic exposures. Its
parallel-cohort closure uses one common temperature and exposure increment;
averaging ages is unavailable. The original caloric inverse also serves this
population, and immutable cohort storage is shared across trial steps.
`CrystallizationPropertySource::for_population` and `evaluate_population`
resolve the actual aggregate mass and density-derived volume fractions on the
same pinned claims, without manufacturing a single-history phase state.
An explicitly admitted aggregate property law remains necessary: this API does
not infer conductivity/stiffness homogenization, inter-cohort impingement,
spatial mass flux, remelting or calibrated molecular mixing.

### Bulk conduction on kinetic interface exchange

`thermal::interface::crystallization::CrystallizationBulkConduction` prepares
both original conductivity claims with actual kinetic coordinates and fixed
positive layer distances. Supplied pressure coordinates must agree with the
parcel pressure. `CrystallizationInterfaceContact::new_with_bulk_conduction`
binds that preparation to the same ordered caloric/kinetic laws and evaluates
both conductivities and the original interface resistance inside each existing
implicit heat candidate. The shared series-resistance owner computes total
conductance. Zero interface resistance is admitted with bulk resistance;
supported zero conductivity insulates without freezing phase history.

Accepted responses retain both bulk source answers and the interface receipt.
Unsupported domains cannot select another claim or publish either parcel.
These stationary isotropic layer laws do not infer moving geometry, spatial
temperature profiles or material calibration. Focused regressions specify
analytic heat/latent balances, distinct crystal coordinates, both resistances,
insulation and transactional source refusal.

`CrystallizationInterfaceContact::new_with_crystallinity` additionally refreshes
separate optional crystalline mass and volume coordinates for both ordered
interface sides before every pinned resistance query. Growth and partial
remelting therefore change conductance at their actual candidate state, including
on an isothermal melting plateau. Volume fractions use the original phase
densities; neither relative Avrami conversion nor liquid fraction substitutes
for crystalline mass fraction. Original chart and exact-Cp pairs share this
path, including optional bulk conduction and the existing atomic heat/history/
pressure-work transaction. Unsupported joint source domains refuse without
switching claims or publishing either parcel. Applying each uniform bulk phase
inventory to its interface is an explicitly admitted source assumption, not a
prediction of surface enrichment or morphology. Manufactured regressions check
independent remelting and growth heat balances, density-dependent coordinates,
accepted receipts, unit/axis admission and transactional source refusal.

## Original mixed pressure dual without unused remainder tangents

`MixedFiniteStrainPotential::mixed_pressure_potential` exposes the original
scalar pressure/volume dual for prescribed-displacement pressure solves. It uses
the same equilibrium law, deformation/stretch checks, finite-pressure checks,
positive reconstructed volume and bulk-branch restrictions as the full mixed
response, including the original log-squared convex branch. It does not evaluate
the remainder Hessian or substitute a new bulk model. The native residual owner
must separately admit remainder/history; this scalar response is not a complete
material-state admission. Full mixed-state and tangent APIs retain their existing
semantics. The native solid regressions compare both original bulk laws against
the full pressure block and independent residual derivatives; execution remains
unverified under the workspace resource blocker.

## Frozen source thermoelasticity in native small-strain mechanics

`ResolvedThermoelastic::try_mechanical_law` prepares a borrowed
`FrozenThermoelasticLaw` for existing plane, tetrahedral and hexahedral
`SmallStrainLaw` consumers. It retains the original temperature, complete
elastic/expansion source responses, full oriented stiffness and independently
integrated thermal strain. Per-element views admit nonuniform temperatures.
Native plane stress relaxes all three transverse strains, including thermal
shear; total and elastic Mandel strain retain the original separate bounds.

Preparation converts the work-conjugate tangent once. Constitutive evaluations
use fixed arrays without source queries, matrix inversion or thermal derivative
work. The same shared Mandel/tensor-Voigt helpers serve frozen acoustic elasticity.
The borrowed source's `state` remains the free-energy, entropy, internal-energy
and thermal-derivative owner. A temperature update requires a fresh source
response; this fixed-temperature view supplies no heat step, parameter-work
account, finite-strain or geometric-prestress model. Focused regressions cover
potential derivatives, anisotropic free expansion through P1/Q1 plane solves,
and a nonuniformly heated, loaded body with independent segmentwise response.

## Species exchange with separately owned internal energy

`ClosedGasParcel::exchange_species_boundary_relative_with_energy` uses the
original species-boundary inventory and caloric reconstruction with one extra
internal-energy coordinate. A callback receives the admitted final total mass
exactly once before temperature reconstruction. This permits a coupled momentum
owner to remove its actual bulk kinetic change from total stream energy without
predicting a different final mass, inventing a species, or reconstructing gas
twice. Explicit external heat/work may use the same coordinate when its owner
closes the opposite account.

The extra increment is reported separately and included in the original total
energy budget. Formation constants still cancel before caloric inversion, and
all source, mass, species, temperature, cancellation and representability gates
remain active. Existing entry points supply exactly zero extra energy and retain
their original operation tree. Focused regressions specify zero-extra parity,
one callback at actual final mass, heat-only use without species placeholders,
reference shifts and early/late atomic refusals.

## Moving molecular interface mass, momentum and energy

`sublimation::gas::kinetic::evaluate_molecular_interface` resolves one local
solid/gas interface from the original finite gas snapshot and retained species
catalog. The original phase law supplies emission and sticking; the shared
drifting-Maxwellian kernel supplies complete incident moments. An explicitly
declared Maxwell model splits nonsticking molecules into specular and diffuse
branches. Carrier pressure, phase recoil, sticking momentum and reflected
normal/shear loads therefore share the same molecular account. These complete
loads must not be added to a second pressure load for the same molecules.

Original species enthalpy differences cancel formation constants before energy
evaluation. Each row reports surface-relative and laboratory momentum, energy
beyond the original gas partial enthalpy, mechanical wall power, and transported
bulk kinetic energy. Laboratory transformation is P_lab=P_rel+w*j and
q_lab=q_rel+w·P_rel+|w|²*j/2. The condensed body's variable-mass force is
−P_rel times actual exposed area; it excludes transported bulk momentum. Patch
torque uses the explicitly declared centroid and moment origin. The near-zero
normal-drift net phase rate retains the original saturation difference and the
separately computed drift increment instead of subtracting rounded gross beams.

Pure solids use the original h/T callback. An explicit sharp fusion chart uses
its original solid constituent enthalpy and density; bulk mixed-phase enthalpy
is not substituted for the solid/vapor latent anchor. All-liquid states refuse,
and exposed solid area is never inferred from bulk phase fraction. Catalog
source identity, complete scattering coverage, temperatures, geometry, finite
arithmetic and cancellation are admitted before returning the rate.

This operation is an instantaneous boundary constitutive query. It does not
evolve finite inventories, bulk gas/body momentum or moving gas volume. Such a
consumer must bound gross withdrawals, reconstruct original thermal stores and
subtract actual bulk kinetic change from laboratory stream energy. No extra
pressure work, fabricated caloric law or continuum flow solution is supplied.
Focused regressions specify stationary-owner agreement, independent drifting
mixture quadrature, diffuse/specular slip, inertial covariance, energy-reference
invariance, sub-ulp saturation drift, fusion and atomic source/frame refusal.

## Explicit Landau two-fluid material response

`superfluid::TwoFluidState` admits immutable, explicitly supplied SI density,
temperature, superfluid fraction, specific entropy, Cp, adiabatic sound speed,
signed thermal expansion and normal viscosity. Normal and superfluid fractions
must both be positive. Entropy/Cp use the total mass basis. Thermodynamic
compressibility and Cv follow the same inputs; the two finite-expansion Landau
sound roots use a scaled discriminant and a product-derived small root.
`SoundModes` reports fast/slow speeds without assuming mode polarization.

Component velocities retain separate mass and entropy currents. Heat current
is relative to the barycentric frame, `rho_s*s*T*(v_n-v_s)`, with no advected
enthalpy or pressure-work duplication. Kinetic storage splits into barycentric
and relative motion; the relative inertia is `rho_n*rho_s/rho`.
Thermomechanical acceleration follows `dmu=dp/rho-s*dT` and its reciprocal
normal momentum equation. Normal momentum diffusivity uses `mu_n/rho_n`.

Explicit linear drag and isotropic Gorter-Mellink drag have exact homogeneous
frozen-state updates, opposite impulses and one nonnegative dissipated-energy
density. That heat is returned for a thermal owner to accept once; this update
does not evolve temperature or silently discard/re-add kinetic energy.
Deterministic `fs_math::det::expm1` resolves weak linear drag without cancellation.
Strong linear and Gorter-Mellink updates retain their independently evaluated
exponential or algebraic counterflow tail, reconstructing component velocities
from conserved barycentric motion instead of subtracting nearly equal values.

The full state implements neither ordinary Newtonian transport nor scalar
acoustics. `first_sound_projection` explicitly retains the original state and
implements `FluidAcoustics` only for zero-expansion, in-phase density waves.
This approximation excludes the entropy mode and classical thermoviscous walls.
No phase diagram, zero-temperature endpoint, vortex topology/nucleation,
critical velocity, turbulent onset, nonlinear spatial evolution, or generic
helium property fit is inferred. Focused G0/G1/G3 tests cover thermodynamic
identities, both modes, fountain balance, currents, Galilean boosts, opposite
friction impulses, analytic decay, semigroup/heat closure and admission limits.

`VortexMutualFriction` retains caller-supplied circulation magnitude, vortex line
density and Hall–Vinen B/B′ coefficients. `relax_hall_vinen` exactly integrates
the frozen straight aligned population: parallel counterflow remains unchanged,
transverse motion damps and rotates, impulses are opposite, and only transverse
damping releases heat. The B/B′ convention is `alpha=rho_n*B/(2*rho)`; the
reactive term is signed and does no work. Scaled direction normalization and
`expm1`/small-angle increments retain weak responses without an explicit-step
stiffness limit. The isotropic unpolarized force projection is
`gamma=(rho_n*rho_s/rho)*B*kappa*L/3`; it is consumed by the full two-fluid
channel and its existing heat account. This is not averaging independently
evolved orientations. Vortex geometry, line density, coefficients and material
state are frozen: no nucleation, reconnection, tension, critical velocity,
frequency-dependent fit, or evolving tangle is inferred. Direct regressions
cover analytic anisotropic decay/rotation, reactive zero heat, weak/strong
limits, momentum/heat closure, semigroup and frame covariance, the angular
force average, and source refusal.

`NormalViscousWall` supplies positive authored permeability on the fluid-volume
basis, with steady circular-capillary `R²/8` and parallel-plate full-gap²/12
helpers. Normal shear viscosity drives `F_n=-eta*v_n/k`, with normal inertia
`rho_n`; superfluid velocity is unchanged and the stationary wall receives
the opposite external impulse. `prepare_normal_wall_relaxation` evaluates
the exact frozen exponential source once for an allocation-free spatial batch.
Its heat is lost normal kinetic energy, separately owned from mutual friction;
the caller deposits it once. Weak impulse/heat, strong decay tails and scaled
rate arithmetic are retained. This reduced law assumes laminar, fully developed,
continuum no-slip flow at low Womersley number. It does not resolve transverse
startup profiles, oscillatory boundary layers, entrance effects, porous
tortuosity, moving walls, or second viscosity. A spatial composition must
separately control its time error. `tests/superfluid_viscosity.rs` checks
Poiseuille resistance, component and wall impulses, heat balance, coflow versus
superflow, semigroup/rotation, weak/stiff limits and finite-domain refusal.

## Longitudinal Landau–Khalatnikov transport

`superfluid::transport::LongitudinalTransport` retains the original
`TwoFluidState`, its normal shear viscosity, explicit Fourier conductivity and
immutable `SecondViscosities`. The signed cross coefficient satisfies
`zeta4=zeta1`; exact binary64 product comparison admits only `zeta2,zeta3>=0`
and `zeta1²<=zeta2*zeta3`, including semidefinite limits without a tolerance.
The constitutive form follows Nikuni–Griffin, Phys. Rev. A 63, 033608 (2001),
equations (86)–(87), and Escobedo et al., Phys. Rev. A 79, 063623 (2009),
equations (8)–(15). No measured transport values are inferred.

With `q0=rho_s(v_s-v_n)`, the scalar longitudinal stress is
`(4 eta/3+zeta2)*div(v_n)+zeta1*div(q0)` and the specific chemical-potential
correction is `-zeta1*div(v_n)-zeta3*div(q0)`. Cached symmetric resistance and
kinetic-energy-normalized diffusion matrices preserve the original total and
reduced component densities. Factored local viscous dissipation stays
nonnegative near rank-deficient coefficients. Conductive temperature diffusion
uses `kappa/(rho Cv)`; Fourier entropy flux and availability destruction are
separate from physical viscous heat, so a thermal consumer must not deposit
conductive availability loss as additional Joules.

Evaluation allocates nothing, uses bounded constant work and refuses nonfinite
or unrepresentable responses. The law has no spatial boundary conditions,
frequency-dependent relaxation, microscopic fit, wall/vortex double counting,
phase transition or nonlinear temperature feedback. G0/G1/G3 regressions in
`tests/superfluid_transport.rs` cover exact PSD admission at overflow/underflow
scales, signed reciprocal coupling, independent stress/chemical-potential work,
the correct density bases, null directions and Fourier entropy accounting.

## Authored Kapitza solid/helium thermal interfaces

`superfluid::kapitza::KapitzaInterface` evaluates the explicit surface law
`Q = area * sigma * (T_from^n - T_to^n)`, with positive heat leaving the first
side. Its immutable `KapitzaParameters` retain actual area, the source's
power-law coefficient in W/(m² K^n), exponent `n >= 1`, and inclusive positive
temperature bounds. A coefficient is never inferred from a material name.
The common quartic relation and its temperature-dependent secant conductance
are described in Bottura, *Cable Stability*, CERN-2014-005, equations (28)-(29),
<https://cds.cern.ch/record/1974064/files/CERN-2014-005-p401.pdf>.

Both sides must remain inside the declared source range. Tangents are positive,
and the lower-endpoint tangent bounds them throughout that interval. Integer
power differences are factored; general exponents use the original deterministic
math owner with a corrected logarithm and `expm1`. Nearly equal temperatures
therefore retain their small heat transfer without subtracting absolute powers.
The evaluation allocates no memory and exposes source parameters, tangents and
the validity range to implicit thermal consumers. The reference-plus-offset
entry point preserves temperature differences below one ULP of the absolute
reference; an error-free sum checks source bounds even when an invalid offset
would otherwise round onto an accepted boundary.

`entropy_production_w_k` reports `Q*(1/T_to - 1/T_from) >= 0`. This entropy
account is not a second heat source: a finite thermal composition must debit
and credit the original heat once. The law has no interface storage and does
not supply bulk helium conductivity, material calibration, boiling or lambda
transition physics. `tests/superfluid_kapitza.rs` checks independent integer
and fractional-power solutions, one-ULP temperature differences, sign and area
covariance, monotonicity, entropy production and domain refusals (G0/G1/G3).

## Continuous gold leaching with oxygen uptake and stripping

`cyanide_leach::flow::aeration::AeratedCyanidationCstr` attaches the original
prescribed-pressure `OxygenTransfer` to an admitted `CyanidationCstr`, including
a tank already advanced without gas. The original charge, stream inventories,
local kinetic card, aqueous chemistry, finite competing-metal chemistry,
elapsed clock and prior liquid/reaction accounts remain intact.

For fixed temperature, solution volume and activity basis, dissolved molecular
oxygen follows `dn/dt = Q*C_in - Q*n/V + kLa*(n_sat-n) - r_Au/4`.
The retained gas/liquid Gibbs sources determine
`n_sat = n_ref*K*p_O2/p_ref`, or the original measured-Henry source determines
`n_sat = V*Hcp*p_O2` using the physical liquid volume. The independent original
transfer card supplies kLa and its validity. Feed, flow and prescribed oxygen partial pressure may
change between advances. Stopping liquid flow permits gas exchange and
dissolution to continue. Neither solubility nor transfer data calibrate the
separate local gold-dissolution coefficient.

The shared transfer-source evaluator checks actual DO, pO2, saturation,
temperature, volume, absolute process clock and additional declared conditions
at evaluated liquid, gas and reaction stages. It cannot overwrite owned
chemical coordinates or turn the local leach card into a batch-age law.
Original source owners, receipts and final admitted conditions remain
accessible. The inner tank is exposed read-only so external advancement cannot
silently bypass its gas accounts.

Symmetric exact positive liquid transport and gas relaxation surround the
existing RK4 reaction step. Gross gas uptake and stripping are accumulated
independently of liquid inlet/outlet and Au consumption. Raw full/two-half
differences control an L1 amount estimate over stocks, exports, Au progress and
both gas totals; no smooth-order divisor is assumed at flux reversals.
The relative allowance scales with all admitted amounts and a gas-exchange
scale, so trace-O2 accuracy requires appropriate explicit tolerances.
These estimates exclude source uncertainty, equilibrium tolerance and model
error; evaluated-stage validity is not a continuous-path enclosure.

Interval and cumulative component/O2 balances include both explicit gas
accounts without modifying the original charge or liquid streams. All checks
and the final cancellation checkpoint must pass before publishing an advance.
No clipping, balance repair, inferred gas-supply capacity, aerator power,
finite headspace, temperature feedback, carbon adsorption or tank-train claim
is made. Focused cases include the independent unsplit flow/gas solution with
a within-interval flux reversal, the zero-transfer/zero-flow limits, gas-enabled
reaction from zero liquid oxygen, finite copper competition and atomic refusal.

## Continuous speciated gold leaching

`cyanide_leach::flow::CyanidationCstr` advances one full, fixed-volume,
isothermal mixed tank. Incoming and outgoing solution volumes are equal.
Accessible solid Au, dissolved Au complex, analytical cyanide, dissolved O2,
background acid families, spectator ions and competing-metal totals all share
the declared hydraulic residence. Feed concentrations and Q use the same
physical solution-volume basis; solids displacement and ore tonnage are not
inferred. The caller's original initial charge remains retained separately
from the canonical internal component layout.

The distinct `LocalGoldCyanidation` source declares
`local-speciated-gold-cyanidation` and the original empirical form
`r_Au=k*n_Au*c_free_CN*c_O2`. Its card is retained unchanged. Batch-age and
residence-time coordinates cannot provide local kinetic authority. Required
typed source coordinates use actual pH, free CN, DO, accessible Au, temperature
and volume. Original water, cyanide and optional metal-family sources own
the shared equilibrium calculation at every evaluated stage.

Transport acts on analytical totals: `dn/dt=Q*C_in-Q*n/V+nu*r_Au`.
The cyanide component includes ligand bound to competing metals and excludes
ligands in the separately tracked Au complex. Metal and background inventories
can enter and leave; their actual partition is re-solved after transport and
reaction. Every initial family, spectator identity/charge and metal total
must have exactly one feed entry, including explicit zeros.

Exact positive transport, bounded RK4 reaction and exact transport form a
second-order split. Full/two-half differences divided by three estimate local
error in the L1 norm over inventories, independently accumulated exports and
Au reaction progress. Time-weighted allowances bound the sum of accepted
estimates in each call. These are numerical estimates, not certified global
bounds or equilibrium/source uncertainty. Source admission covers numerical
split endpoints and reaction stages, not an enclosure of the unsplit path.

Interval and cumulative inlet/outlet accounts remain separate from reaction
progress and final tank inventory. Independent component, Au, cyanide, O2 and
Na balances retain the original ligand-equilibrium residual. A flow stop
preserves inlet/outlet accounts while chemical dissolution can continue.
Invalid feeds, source refusal, depletion, unresolved numerical updates,
work exhaustion and cancellation publish no partial state or stream account.

Oxygen is supplied only by initial liquid and explicit feed. Gas transfer,
carbon adsorption, solids slip, mineral liberation, competing-metal redox,
precipitation, variable volume, heat feedback and finite solvent H/O closure
remain outside this owner. The source-backed example and input conventions
are in `data/matdb/seed-v1/MINING_CONTINUOUS_CYANIDE_LEACH.md`.

## Simultaneous carbon-in-leach and explicit preg-robbing inventories

`cyanide_leach::cil::CarbonInLeachBatch` advances accessible-gold dissolution
and reversible capture by independently declared retained carbon populations
during the same finite-bath time interval. Each `CarbonPopulation` retains
its original analytical-Au Langmuir isotherm, local linear-driving-force
source, dry mass, initial loading, identity and extra source conditions.
`ProductCarbon` and `PregRobbing` are explicit recovery destinations; neither
label supplies a capacity, affinity or kinetic coefficient. Reports retain
both destinations separately and use signed net interval capture, including
desorption from preloaded carbon.

Attachment takes a fresh original `GoldCyanidationBatch` with zero clock and
zero gold progress. Accepted initial NaCN/NaOH doses and optional competing
metal families survive attachment. Original accessible gold, solution volume,
temperature, acid/base sources, activity-reference amount and finite dissolved
oxygen remain authoritative. Every stage reuses the original live pH/free-CN
and finite metal/ligand solve. The current leach law has no product-inhibition
factor: neutral ion-pair capture must not invent a pH benefit or faster gold
dissolution by returning bound cyanide to the free ligand pool.

`IonPairTransfer` explicitly declares
`Au(CN)2- + monovalent countercation -> neutral adsorbed pair`.
The sodium constructor selects `Na[+](aq)` and `NaAu(CN)2(ads)`.
Generic declarations require exactly three distinct participants with written
coefficients -1, -1 and +1; their admitted charges must be -1, +1 and zero.
The countercation cannot contain gold, and the gold-complex formula must
agree with the original conserved Elsner equation. Element and charge
definitions are retained from exact reaction artifacts, never inferred
from labels. The CIL owner also refuses overlapping solution-family
ownership and requires an explicit aqueous countercation spectator.

Adsorption moves one gold atom, two intact cyanide ligands and one
countercation together; desorption returns the same pair. Preloaded carbon
brings these inventories with it. The initial aqueous countercation must
cover all initially dissolved and potentially dissolvable gold. This finite
reserve restriction avoids inventing a counterion-limited kinetic law.
The supported neutral-pair equation is an explicit bookkeeping/model choice,
not a universal microscopic adsorption mechanism or an adsorbed formation
thermochemistry claim.

At every rate stage and endpoint the actual dissolved-gold concentration,
loading, solution volume and semantic absolute temperature are admitted by
the original isotherm/local-rate owners. CIL additionally supplies actual
`pH`, `free-cyanide-concentration`, elapsed time, dissolved O2 and the
original leach inventory/activity-reference coordinates, plus
`ion-pair-counterion-concentration`, to both sources. Caller overrides of
any coordinate owned by the original leach or adsorption models are refused. Source validity is checked at the actual
state; no extrapolation or clipping is permitted.

The bounded full/two-half RK4 controller retains separate dissolution
progress, dissolved complex, aqueous countercation, dissolved oxygen,
adsorbed amounts and vacant capacities. It recomputes Au, analytical cyanide,
countercation, oxygen and per-population site balances independently.
The cyanide report retains the original explicitly bounded metal/ligand
equilibrium residual. Numerical full/two-half differences divided by 15
estimate local inventory error; the time-weighted allowance and accumulated
estimate are not certified bounds or experimental uncertainty. Nested
speciation shares the same cancellation/work budget. A source exit,
cancellation, exhausted budget or unrepresentable inventory leaves the
complete accepted state and clock unchanged.

This owner is isothermal with fixed liquid volume and retained dry masses.
It has finite initial dissolved O2, with no gas replenishment, further dosing,
flow-through residence-time model, carbon movement or regeneration. Those
require explicit additional owners. The runnable
`simultaneous_gold_carbon_in_leach` example and
`data/matdb/seed-v1/MINING_CARBON_IN_LEACH.md` combine original NBS aqueous
chemistry with clearly marked authored adsorption/leach coefficients; they
do not claim calibrated ore or carbon performance.

## Source-admitted instantaneous aqueous adsorption

`AqueousLangmuir::loading_at` evaluates the original isotherm at actual
dissolved analytical-element concentration and absolute temperature. The
report retains its original model and admitted source point, together with
independently evaluated occupied and vacant equilibrium loading.
`LocalLangmuirTransfer::rate_at` evaluates the signed local loading
derivative at actual concentration, loading, volume and temperature.
Equilibrium and kinetic sources remain independently admitted. Near
saturation the derivative uses the difference of vacancies to preserve
representable trace desorption. Caller-supplied extra conditions cannot
override coordinates owned by these methods. These operations expose the
original constitutive responses for coupled solvers; they do not impose
finite-batch equilibrium or add measured kinetic coefficients.

## Dissolved cyanide and caustic additions during leaching

The closed, finite-headspace and prescribed-pressure cyanidation owners expose
`dose_reagent` for finite, already-dissolved `SodiumCyanide` and
`SodiumHydroxide` equivalents. NaCN adds analytical cyanide and sodium; NaOH
adds sodium/alkalinity. The original shared acid/base owner re-solves pH,
HCN, free CN and OH. Hydroxide is not installed as a second fixed spectator,
and a dose never assigns a target pH or assumes all cyanide remains free.
Optional competing-metal families retain their original totals and sources;
their ligand partition is re-solved against the increased analytical pool.

The immutable initial charge remains distinct from cumulative `ReagentInputs`.
At cumulative gold dissolution `x`, the available cyanide before partition is
`T0+n_NaCN-2*x`, and sodium is `Na0+n_NaCN+n_NaOH`. Dose and advance receipts
retain the original model and independent cyanide/sodium residuals. Adding a
reagent preserves elapsed time, gold progress, actual liquid/gas oxygen and
the existing oxygen input/uptake/stripping accounts. The finite chemical
ceiling expands with cyanide input; a closed bath still cannot consume more
than its original dissolved oxygen supplies.

Each dose prepares and admits a complete candidate state before publication.
Negative, nonfinite or numerically unresolved additions, conflicting sodium
ownership, source-domain refusal, exhausted work and cancellation leave the
accepted state and cumulative inputs unchanged. The finite-headspace owner
also re-admits its retained transfer source at the changed chemistry. The
prescribed-pressure owner does not retain a boundary between advance calls;
its next advance admits that call's pressure and transfer source against the
updated solution.

The boundary is fixed temperature and volume with caller-declared dissolved
equivalents. Solid dissolution completion, commercial purity, carrier liquid,
heat of mixing, cyanide destruction and an optimal operating dose are not
inferred. The NBS source-mass conversion and runnable `--dose-reagents`
example are documented in
`data/matdb/seed-v1/MINING_REAGENT_DISSOLUTION.md`.

## Competing dissolved metals during gold cyanidation

`GoldCyanidationBatch::new_with_competing_metals` couples explicit dissolved
metal association families to the existing finite leach owner. Its
`CompetingMetals` input supplies original source-backed families, every finite
metal total, and bounded ligand/pH solve controls. The coupled system uses the
original HCN/CN and water owners; shared chemical definitions, Gibbs references
and receipts must agree. Metal counterions are explicit charge inputs. Neither
an ore copper assay nor a Cu(II) source becomes a dissolved Cu(I) inventory
without a separately declared physical and chemical model.

The supplied charge's `uncomplexed_cyanide` is HCN+CN before optional metal
equilibration, excluding ligands in any initial gold complex. At cumulative
gold progress `x`, the remaining available ligand total is `T0+n_NaCN-2*x`.
Every rate stage jointly partitions it among HCN, free CN and the declared metal
complexes while solving electroneutrality with the current gold-complex stock.
The actual free CN and pH enter the original source-admitted empirical rate;
the owned uncomplexed-cyanide source coordinate is the actual HCN+CN pool.
Metal distributions can change as gold consumes ligand. Fixed metal totals
are conserved; binding does not destroy cyanide or supply a fitted consumption
factor. Zero metal totals remain explicit and recover the original chemistry
within the declared equilibrium tolerances.

Each state retains the complete ligand equilibrium, its source owner, metal
stocks, ligand counts, charge and component residuals, numerical bracket,
and work counts alongside the original aqueous report. The
independently recomputed cyanide account includes metal-bound ligands and
retains the signed, tolerance-checked ligand-equilibrium residual. Subtracting
that reported numerical residual must leave only roundoff accounting error;
stocks are never clipped or normalized to close the account. Kinetic local
error estimates remain separate from equilibrium/source/model uncertainty.

Both the finite-headspace and prescribed-pressure oxygen adapters carry this
same competition through attachment and continued advance. The latter's
`from_batch` preserves an already admitted fresh bath and its metal families.
Leach and oxygen-transfer source conditions use actual competing-solution
coordinates. All original atomic source/work/cancellation refusals still
apply. This coupling does not infer redox, mineral dissolution, CuCN solid
equilibrium, irreversible cyanide oxidation, carbon uptake, nonideal activity,
or an ore-specific kinetic calibration.

## Prescribed-pressure oxygen contact during cyanidation

`cyanide_leach::aeration::AeratedGoldCyanidationBatch` couples the original
speciated gold law to a maintained oxygen partial pressure. The existing
`IdealGasLiquidPartition` supplies the equilibrium dissolved amount from
original gas/aqueous formation-Gibbs receipts, explicit solvent activity basis,
species partial pressure and the retained standard reference pressure.
Alternatively, `OxygenTransfer::from_henry_source` retains the same
`gas_transfer::OxygenHenrySource` used by the finite-headspace owner. Its
explicit concentration/pressure convention gives `O_sat = Hcp*p_O2*V`; the
independent solvent activity amount never substitutes for physical volume.
Batch, continuous CSTR and simultaneous carbon-in-leach routes share this
evaluation and preserve the original pinned Hcp receipt in `henry_answer`.
The explicit-pressure CIL reagent-dose path uses that same post-dose query,
retains its receipt, and leaves the gas and carbon inventories unchanged.
The original Gibbs accessor is `Some` only for its own route; measured data do
not synthesize formation energies or a Gibbs reference convention. A
separate version-one `prescribed-pressure-oxygen-transfer` source card supplies
nonnegative `volumetric-transfer-coefficient` kLa in inverse seconds; no rate is
inferred from equilibrium data.

At fixed volume, `J=kLa*(O_sat-O)` and the coupled balances are `dx/dt=r_Au`,
`dO/dt=J-r_Au/4`. Actual DO is independently integrated, so an initially
oxygen-free charge can start reacting after uptake. The finite chemical
ceiling is `min(Au0,(CN0+n_NaCN)/2)`; closed-batch initial DO is not an
open-boundary ceiling. The original chemical evaluator, simultaneous acid/base speciation,
conserved equation and empirical source admission remain authoritative.

The owner separately integrates gross uptake and gross stripping and preserves
interval and cumulative accounts over calls. Original stoichiometric oxygen
consumption independently closes `O+x/4+out=O0+in`. The reported accumulated
arithmetic allowance never repairs a stock. Partial pressure can change between
advances; zero pressure strips under active contact, whereas zero kLa stops
exchange. Finite headspace and metered dosing belong to the existing
`gas_transfer` owner, not this maintained-reservoir boundary.

Both kinetic cards admit actual conditions at every stage and endpoint.
Caller coordinates cannot override chemistry, gas pressure or saturation.
Gibbs gas and chemical sources share temperature, reference pressure and elemental
convention; the exact network joins original species definitions. Measured Hcp
requeries its retained selected claim at the actual aqueous temperature, with
the original explicit pure-water, molecular-O2, infinite-dilution and unqualified
engineering-reference flags required in advance conditions. The source's
uncertainty, observation references and temperature support remain unchanged;
no saline/slurry correction or empirical ore calibration is inferred. The external
gas identity is excluded from liquid families and spectators, even for empty
inventories. Original gas-equilibrium errors are retained in the public refusal.

Focused G1 fixtures cover uptake and stripping against exponential solutions
with physical volume deliberately different from the activity amount, an
independent unsplit flow/gas target, and measured-Hcp versus equivalent-Gibbs
simultaneous leaching/carbon capture. G4 fixtures retain the entire batch and
oxygen accounts on missing/contradictory source flags and late cancellation.

Bounded full/two-half RK4 trials evolve gold, DO and both gross transfer accounts.
Raw componentwise discrepancies control local numerical estimates without
assuming smooth error scaling across uptake/stripping reversal. Each call starts
fresh interval transfer coordinates, so large historical totals cannot hide
small increments. Negative or unrepresentable positive stocks, unresolved clock
stages, unsupported source points, exhausted budgets and cancellation refuse
without publishing any candidate state. Numerical estimates are not certified
global bounds or calibration uncertainty.

The prescribed boundary asserts no finite gas-feed capacity, gas power, changing
volume, heat feedback, nonideal activity correction or measured ore calibration.

## Continuous population-balance grinding with classifier recycle

`comminution::population_balance::circuit::ContinuousGrindingCircuit` closes
the original continuous perfect-mixer selection/breakage law around the
original classifier and an explicit common underflow return/purge partition.
At prescribed residence, a coarse-to-fine triangular solve retains dry mass
and every original elemental channel. Water has its own recycle closure.
The positive escape coefficient uses independently retained overflow/purge
fractions. Fed traps refuse; unforced trapped circulation is zero under the
existing circuit convention.

The model reconstructs actual circulating feed and reevaluates both original
physical owners. Their source applicability is checked at actual residence,
throughput, diameter and solids fraction. Original circuit tolerances govern
recomputed local and external balances; excess circulation or any failed
source/closure check refuses the whole result. No source values are normalized
or caller coordinates substituted to force admission.

`ContinuousCircuitSolution` retains the complete original model, existing
`CircuitSolution`, and original `ContinuousMillSolution` with class/water
hold-up and source query. The additional `at_dry_holdup` solve requires a
positive target dry inventory, explicit residence bracket, mass tolerances
and bounded iteration budget. It evaluates the original admitted circuit at
every trial and returns an actual converged circuit with residual diagnostics.
Unbracketed targets, exhausted work, unresolved numerical intervals and
original source/circuit refusals produce no accepted result; the solver does
not widen brackets or relax source support. No selection-rate calibration,
power draw, capacity, liberation or plant recovery is inferred.
Seven focused regressions and `continuous_grinding_circuit` exercise analytical
answers, assays, water, purge, traps, applicability, limits and residence.
Seven additional hold-up regressions independently check the inferred
residence and product distribution, water exclusion, purge/no-recycle limits,
source admission, infeasible targets, work exhaustion and numerical refusal.

`population_balance::circuit::power::ContinuousBondCircuit` composes this
original circuit with an original `BondGrindingCard`. Prescribed-residence
and fixed-dry-inventory calls independently solve the circuit, then recompute
fresh-feed F80 and the combined external overflow-plus-purge P80. The original
Bond owner admits that actual size reduction; its specific energy multiplies
fresh dry throughput, excluding circulating solids and carrier water. Results
retain both original models, complete circuit streams and all available
hold-up diagnostics. Caller-edited public results cannot stand in for admission.

Six focused tests cover analytical energy/power, an inventory-controlled P80
change under doubled feed, purge, fresh-tonnage accounting, original size/source
refusals and preservation of original circuit/hold-up errors. Discrete P80 can
jump as throughput changes, so this adapter makes no monotone power/capacity
inverse or selection-rate/power calibration claim. Its net pinion estimate
supplies no equipment power-draw law, drive efficiency or auxiliary demand.

## Finite gas/liquid oxygen supply during cyanidation

`cyanide_leach::gas_transfer::AeratedGoldCyanidation` attaches a finite,
fixed-volume ideal oxygen headspace to the original admitted cyanidation bath.
The original charge, source-backed acid/base families, empirical dissolution
law, reaction progress and absolute clock survive attachment. For the Gibbs
partition route, the original reaction-network owner joins exact O2 species
definitions; matching labels alone cannot admit a different molecular formula.
Those thermochemical references must share temperature, reference pressure and
elemental convention.

`OxygenTransfer` binds either the original `IdealGasLiquidPartition` or a
measured `OxygenHenrySource` to a separate `local-oxygen-gas-liquid-transfer`
source card declaring kLa in inverse seconds.
At each stage, current gas stock supplies ideal partial pressure and hence
equilibrium dissolved amount D*. The explicit local law is `J=kLa*(D*-D)`.
With original gold rate r, the coupled balances are `dx/dt=r`,
`dD/dt=J-r/4`, `dG/dt=-J`; another coordinate independently integrates signed
gas-to-liquid transfer. Thus oxygen supply changes while reaction consumes it.
An explicit contact switch stops exchange without stopping dissolved-oxygen
consumption. Metered `dose_oxygen` calls add actual O2 to the gas stock between
advances and retain a separate cumulative input account.

`OxygenHenrySource::from_pack` resolves an observation-backed Hcp claim in
mol/(m3 Pa), retaining the original pack, content hash, typed query and receipt.
The selected claim itself must declare semantic absolute temperature and the
five molecular-O2/pure-water/infinite-dilution/engineering-reference/unqualified-
specimen flags. Unsupported additional source axes are refused. Each stage
requeries the pinned claim using actual temperature and all declared flags;
its equilibrium amount is `Hcp*p_O2*V_liquid`, with physical liquid volume
independent of the activity-reference amount. The accepted state retains
`henry_answer`. This route supplies no formation Gibbs data or elemental
reference convention: `equilibrium()` returns `None`, while `henry_source()`
returns the original source. Applying the pure-water datum to a cyanide liquor
is an explicit engineering approximation and supplies no salinity correction,
kLa or ore kinetics.

Bounded RK4 full/two-half trials retain the positive two-half solution without
clipping or correcting balances. Difference/15 is measured in the stated
Au-equivalent L1 norm `|dx|+4*(|dD|+|dG|+|dTransfer|)`; time-weighted allowances
control the sum of accepted estimates within each call. Accumulated estimates
are reported since attachment, separate from earlier closed-batch error. They
are not certified global bounds or source/speciation uncertainty. Every stage
and final point must satisfy both original kinetic sources. Transfer-source
pH, free CN and other known changing chemistry coordinates use the actual
solution owner; caller-fixed extras cannot replace them. Source, range, work
or cancellation refusal publishes no partial state, clock or dose account.

Results retain the original models, source conditions, gas/liquid O2, actual
free CN/pH, Au inventories, signed interval/cumulative transfer, metered doses,
original equation amounts and independent Au/CN/liquid-O2/gas-O2/total-O2
residuals. Original-equation OH formation remains distinct from equilibrium OH.
Fixed temperature, liquid/gas volumes, solvent activity basis and ideal phase
partition are explicit. The model does not infer continuous air flow, venting,
gas purity, bubble hydrodynamics, salinity correction, slurry kLa calibration,
finite solvent H/O closure, energy balance or simultaneous carbon adsorption.
The actual-source `aerated_gold_leach` example and source compilation are
documented in `data/matdb/seed-v1/MINING_OXYGEN_TRANSFER.md`.

## Finite gold cyanidation with coupled aqueous speciation

`cyanide_leach::GoldCyanidationBatch` advances accessible gold dissolution
against finite uncomplexed HCN+CN and dissolved-O2 inventories. Its separate
kinetic card declares the empirical law `r_Au=k*n_Au*c_free_CN*c_O2`, with
`k` in m6/(mol2 s). These first-order factors are explicit model choices;
equilibrium data, chemical stoichiometry and ore names supply no rate constant.
Accessible gold is a declared reactive inventory, not inferred liberation.

Cumulative progress is mol Au dissolved. The original conserved dissolved-O2
Elsner equation supplies consumption/formation at extent `x/4`. The original
`IdealAcidBaseSystem` jointly solves all admitted acid families and water at
every rate stage. In the absence of optional metal competition, the remaining
HCN+CN total is `T0+n_NaCN-2*x`; dicyanoaurate is an
internally updated spectator with amount `G0+x`. Thus HCN remains available
through dissociation, and free CN is never confused with total cyanide.
Original Gibbs-source receipts, initial counterions and other family totals
remain attached to the model and accepted state. Source temperature must
match the actual fixed bath temperature. The common activity amount and
liquid volume are separate declarations.

Required source axes cover actual pH, free CN, DO, accessible gold, absolute
time, volume and temperature; optional owned inventory/activity axes cannot
be overridden by caller coordinates. Bounded RK4 full/two-half trials retain
the positive two-half result without extrapolation or stock clipping.
Time-weighted local difference/15 estimates control accepted numerical work;
they are not certified global bounds and exclude speciation/model uncertainty.
Every stage is source admitted. Cancellation, exhausted work, unresolved
progress or source refusal leaves the full accepted state and clock unchanged.
Reports separate interval gold dissolution, cumulative equation extent,
stoichiometric amounts, actual species speciation and analytical Au/CN/DO
balance residuals.

Water retains fixed unit activity. Stoichiometric OH production is distinct
from equilibrium OH, and no finite solvent H/O or thermal balance is claimed.
The scope excludes gas transfer, cyanide destruction, other
cyanicide minerals, passivation, through-flow and calibrated ore/plant recovery.
The `gold_leach_to_carbon` consumer hands analytical dissolved Au to the
existing finite-time carbon owner and checks the combined solid/liquid/carbon
Au balance. It specifies carbon-stage pH externally and supplies no ligand or
counterion adsorption chemistry or simultaneous CIL model. See
`data/matdb/seed-v1/MINING_CYANIDE_LEACH_KINETICS.md` for source compilation
and the complete coupled example.

## Continuous solution flow through retained sorbent populations

`aqueous_adsorption::flow::AdsorptionCstr` advances one full, ideally mixed
solution tank containing one to 128 retained dry sorbent populations. All
populations share the actual dissolved inventory of one analytical element.
The original `AqueousLangmuir` owner supplies equilibrium loading; a distinct
`local-aqueous-langmuir-linear-driving-force` source supplies the inverse-second
rate in `dA_i/dt=k_i*(M_i*q_eq(D/V)-A_i)`. The solution balance is
`dD/dt=Q*(c_in-D/V)-sum(dA_i/dt)`. Each population retains its own original
isotherm, rate, dry mass, loading and additional chemistry conditions. Batch
age and residence-time correlations cannot silently supply a local rate.

Piecewise-constant feed intervals retain the accepted dissolved, occupied and
vacant stocks, elapsed time and cumulative inlet/outlet accounts. Loaded
sorbents can desorb while another population adsorbs. With the pump stopped,
sorption continues while inlet and outlet are exactly zero. Carbon remains
retained; the solution outlet contains only the dissolved analytical element.

Backward Euler reduces each implicit step to one monotone common-concentration
balance. Occupied and vacant stocks use bounded convex updates without
clipping. Outlet is independently integrated as `h*Q*D_new/V`, never supplied
by closing the balance. Accepted states use two half steps; their comparison
with a full step uses divisor one for this first-order method. The tolerance
applies to each accepted local L1 stock/outlet comparison. Its accumulated
estimate is reported separately and is not a certified global-error bound.
Original source support is checked at every full- and half-step endpoint;
this does not enclose the exact continuous trajectory. Step, scalar-root and
population work are bounded. Any refusal or cancellation publishes no partial
inventory, clock or cumulative account.

Results expose original models and source points, signed interval/cumulative
transfer, retained dissolved and sorbent stocks, independent outlet integration,
interval/cumulative balance residuals and numerical work. Fixed liquid volume,
temperature, dry masses and externally declared chemistry are explicit. The
owner does not infer charged-complex removal, counterion transfer, simultaneous
leaching, carbon movement, regeneration chemistry or plant calibration.

Eight focused tests cover independent washout and finite-bath solutions,
first-order refinement, stationary loading, heterogeneous retained populations,
stops/restart and desorption, scaling/subdivision/order/continuation, trace
stocks and atomic source/work/cancellation refusal. The
`mineral_adsorption_flow` example runs four continued feed/stop/wash intervals
and prints the complete metal accounts. Its coefficients are explicitly
synthetic. Run with
`cargo run -p fs-material --example mineral_adsorption_flow`.

## Finite-time aqueous adsorption with competing sorbent populations

`aqueous_adsorption::kinetics::AdsorptionBatch` adds persistent reversible
uptake/desorption to the original `AqueousLangmuir` equilibrium owner.
One to 128 independently declared sorbents share one finite dissolved
inventory of the same analytical element. The original equilibrium source
supplies capacity and affinity; a separate kinetic card explicitly supplies
the constant inverse-second rate in `dq/dt=k*(q_eq(c)-q)`. Every population
retains its original dry mass, loaded stock, source cards and conditions.
No rate or preg-robbing response is inferred from a carbon name.

An implicit step reduces to one monotone scalar balance for the common
dissolved stock. Occupied and vacant sites use bounded convex interpolation
from the nearer endpoint, preserving identical endpoints exactly. Actual
source conditions are checked at intermediate and final accepted states.
Full-step/two-half-step L1 comparison selects the positive two-half result
without extrapolating, clipping stocks or normalizing conservation. All root,
step and nested work is bounded; failure or cancellation leaves the complete
accepted bath and absolute clock unchanged.

Reports retain original models, signed interval/cumulative transfer, complete
dissolved/site inventories, source points, conservation residuals and numerical
work. Loaded sorbents may release metal while another population takes it up.
Local and accumulated error estimates are separate diagnostics, not certified
global bounds or source/model uncertainty. Fixed solution volume, dry masses,
temperature and external chemistry are explicit. This owner does not supply
aqueous speciation, competitive co-adsorbate chemistry, simultaneous dissolution,
flow, heat, pore-diffusion identification or calibrated CIP/CIL recovery.

Eight new focused groups check an independent finite-bath analytical curve,
equilibrium approach, multiple sorbents, preloading/desorption, time/rate and
extensive scaling, subdivision/order behavior, continuation, trace inventories
and atomic source/work/cancellation refusals. All 124 tests in the original
separation, comminution, reaction_batch and aqueous_adsorption native module
scope pass, including the six original aqueous-equilibrium tests. The actual
`mineral_adsorption_kinetics` example compiles and runs both synthetic
activated/natural-carbon scenarios through three continued intervals; its
initially unloaded natural carbon later desorbs as the shared bath depletes.
The model and commands are documented in
`data/matdb/seed-v1/MINING_SOLUTION_RECOVERY.md`.

## Persistent size-dependent leaching with one finite reagent bath

`reaction_batch::shrinking_core::population::ShrinkingCoreBatch` advances
1..=128 explicitly assigned mineral/size populations in a single closed,
well-mixed bath. Each population keeps its original conserved equation,
source model, reactive-grain radius, original amount and accepted core.
Film, surface-reaction and product-layer diffusion resistances remain in
series. All local laws must consume the same reagent; products may not feed
any modeled reactant in this irreversible closure.

For penetration `u=1-r/R` and conversion `X=1-(1-u)^3`, the original local
law integrates to
`F(u)=(a/b)*rho*R*[r_f*X/3+r_s*u+R*r_d*u^2*(3-2u)/6]`.
Every unfinished population receives the same incremental exposure
`H=integral(C dt)`, so `F(u_new)-F(u_old)=H`. A scalar bath-clock integral
closes reagent depletion against all original mineral stocks. Size populations
therefore compete through actual reagent consumption, without processing order
priority or clipping independent excess-reagent conversions. Positive
incremental/tail polynomials and bounded inversions retain remaining cores;
quadrature splits at exact population-completion exposure events.

Repeated calls retain absolute elapsed time, original mineral-total loading
(summed over its size bins), shrinking cores and the depleted bath. Dividing
one identical bin into two does not change the source loading coordinate.
Repeated reaction IDs require the same original terms and conservation
certificate; per-bin extents aggregate into the existing exact reaction
network. Reports expose **cumulative** network consumption/formation from the
original complete inventory, per-population remaining amounts, geometry/molar
residuals and original source points. Exposure continues accumulating during
a post-completion hold. Invalid stock assignments, unsupported source domains,
nonrepresentable positive stocks, work exhaustion or cancellation leave the
accepted state untouched.

The caller bounds all work, including nested core inversions, and chooses
absolute/relative time tolerances. Adaptive quadrature and root residuals are
estimates, not interval or global-error certificates. Fixed volume,
temperature and coefficients, spherical fixed outer radii, additive local
rates and one shared reagent are explicit assumptions. This owner does not
infer mixture calibration, mineral liberation, activity/speciation, product
inhibition, gas transfer, pore evolution, fragmentation, energy coupling or
metal recovery.

Seven public-API test groups use independent analytical film, diffusion and
mixed-resistance solutions, different reagent stoichiometries, completion
events, permutation/subdivision/continuation equivalence and atomic refusals;
two additional intrinsic-coordinate tests cover tiny progress and positive
tails. The complete original separation/comminution/reaction-batch module
scope passes 110 native tests. The actual `grinding_to_leaching` example
connects the existing selection/breakage mill, source-backed charge adapter
and persistent shared bath. It consumes actual USGS CuO/calcite molar-mass
packs and compares two synthetic grind durations through two consecutive
leach advances while retaining dry-solids, carrier-water, copper and acid
accounts. Its rates, mineral fractions, liberated-sphere interpretation and
bath conditions are explicitly authored illustrations, not measured ore or
plant forecasts.

Run with compiled source packs:
`cargo run -p fs-material --example grinding_to_leaching -- --materials PACK_DIR`.
The two required pack names and assumptions are in the example's module docs.

## Grinding-stream collection for explicit mineral leach charges

`comminution::leach_charge::CollectedLeachCharge` collects an original
`ParticleStream` over an explicit duration and converts caller-assigned
mineral fractions into typed amounts. It retains size-class and population IDs,
the original stream, dry/element/water masses and unallocated mass per class.
Collection duration is distinct from milling residence and leach contact time.

Each allocation names its solid species, fraction with artifact-backed
provenance, and a `MolarMassSource` resolved through the actual typed material
query and receipt verifier. The source pack identity, query, answer and
uncertainty remain inspectable. Neither elemental assays nor a source-pack
label authenticates mineralogy. Reactive radius needs an independently sourced
value or an explicit liberated-spherical-particle declaration allowing
class diameter/2; bulk ore PSD does not silently become mineral-grain geometry.

Canonical population IDs and exact fraction budgets prevent double allocation.
The bounded owner handles at most 128 source classes and 1024 declared mineral
populations; a downstream solver may impose a smaller bound. Positive numerical
underflow, overflow, duplicate IDs, absent classes and source/schema failures
refuse without repair. Unallocated dry mass has unknown mineralogy and is not
certified inert. Original elemental accounts describe the whole charge and
must not be assigned again to either the mineral or residual inventory.

Six focused tests exercise real typed molar-mass evidence/receipts, independent
kg/s-to-kg-to-mol arithmetic, duration scaling, original assays and water,
zero/unallocated controls, radius/fraction evidence and numerical extremes.
This adapter prepares populations for reaction owners; it does not choose
kinetics, infer liberation, certify mineral/assay closure or supply a solution.

## Continuous stirred-tank chemistry with finite shared reagents

`reaction_flow::ReactiveCstr` evolves complete molar inventories in one
isothermal, fixed-full-volume ideal stirred tank. It joins 1..=32 original
`LocalBilinearReaction` laws consuming one explicitly named shared reagent
and distinct other reactants. Products may be shared, but cannot feed a
reactant in this slice. Each conserved equation and species definition remains
owned by `fs-thermochem::ReactionNetwork`; chemistry uses simultaneous gross
stock admission so competing pathways cannot independently spend the same acid.

Each local law declares `r_j = k_j * n_Aj * n_Bj / V` in mol/s with explicit
species roles and `k` in m³/(mol s). The distinct version-1
`local-bilinear-reaction` card requires zero hidden state, original source
hashes/provenance and finite typed volume, temperature and both concentration
domains. Additional source conditions remain mandatory. Batch age, batch
volume and residence-time domains cannot be promoted into this local law.
Mineral inventory divided by volume is an empirical availability variable;
neither reaction orders nor rate constants are inferred from stoichiometry,
mineral identity, pure-solid activity or particle surface area.

`advance` takes complete initial stocks, complete inlet concentrations, actual
volumetric flow and one piecewise-constant duration/temperature/condition set.
Every species, including solids and dissolved products, shares the declared
ideal residence. At zero flow, reaction continues with no inlet or outlet.
The returned full inventory is the state for a later feed change or pump
restart. Transport integrates its own positive exponential outlet amounts;
outlet is never reconstructed as a mass-balance plug or counted as recovered
metal. Results expose initial/final/incoming/outgoing amounts, cumulative
reaction extents and gross species consumption/formation, source endpoint
queries and independently evaluated complete-species balance residuals.
`result.final_inventory` is the final tank state; the nested
`result.reactions.final_inventory` is the cumulative reaction account funded
by initial stocks plus all inlet, before subtracting the separately integrated
outlet. The two inventories have intentionally different accounting scopes.

The integrator composes exact transport half-steps with simultaneous RK4
reaction steps. One full step versus two half steps estimates the generally
second-order splitting error, divided by three. Error controls include stocks,
integrated outlets and reaction extents, with time-weighted per-step allowances,
accumulated estimates and bounded accepted/rejected attempts. These are local
numerical estimates, not certified global errors. Original source support is
checked at each split endpoint and reaction stage; rectangular admission
covers each monotone split path, not an enclosure of the unsplit trajectory.
Stock overdraw, erased positive finite-time reactants, range/resolution loss,
out-of-domain stages or exhausted work refuse without publishing partial state.
Neither negative stock clipping nor renormalized reaction progress is used.

Seven focused regressions cover an independent stopped-flow bilinear solution,
exact nonreactive washout and outlet integrals, shared acid with unequal proton
demands, flow interruption, atom/charge/species balances, steady-state step
refinement, source refusals and atomic budget failure. The runnable
`mineral_leach_train` example carries copper-oxide and calcite stocks through
acid-feed, flushing, stopped-flow and restart intervals with explicitly
synthetic rate declarations. It reports retained/exported copper and acid
accounts. This slice does not supply ore calibration, selective solids
residence, speciation/equilibrium feedback, reagent recycling through reaction
products, heat evolution, gas transfer or changing tank volume.

## Particle-size-dependent shrinking-core leaching

`reaction_batch::shrinking_core::ShrinkingCoreReaction` binds one exactly
conserved equation, explicit solid/reagent species roles and an original
version-1 `spherical-shrinking-core-reaction` constitutive card. For
`a A(aq) + b B(s) -> products`, a monodisperse population has fixed outer
particle radius R and an initially unreacted core. With
`u=1-r/R` and `X=u*(3-3*u+u*u)`, the local progress law is

`dt/du = (a/b)*rho_B*R * [r_f*(1-u)^2 + r_s + R*r_d*u*(1-u)] / C_A(u)`.

The finite, well-mixed bath supplies
`C_A(u) = [n_A0 - (a/b)*n_B0*X(u)] / V`.
The declared resistances are `film-transfer-resistance` and
`surface-reaction-resistance` in s/m, plus
`product-layer-diffusion-resistance` in s/m². They are inverses of the film
coefficient, surface coefficient and effective diffusivity. Each is nonnegative,
at least one is positive, and zero explicitly neglects that resistance.
`solid-molar-density` is positive mol/m³ of the chosen solid formula unit.
No rate, surface coefficient or diffusivity is inferred from mineral identity.

Source admission requires original hashes/provenance, zero hidden state and
six finite typed ranges: elapsed time from zero, temperature, batch volume,
initial particle **radius**, reagent concentration and initial solid loading.
Additional declared conditions remain mandatory and fixed; caller coordinates
cannot replace these owned axes. Initial loading is a retained batch descriptor.
Reagent concentration decreases monotonically, so initial/final admission
covers the path through the rectangular source domain.

`advance` inverts the integrated resistance law with bounded adaptive Simpson
quadrature and bisection. Its reported time error is an estimate, not a rigorous
interval certificate. A shared evaluation budget and cancellation polling bound
work. Nearly stoichiometric short-time requests do not integrate an unnecessary
completion tail; an analytical tail upper bound permits complete solid
dissolution only with sufficient excess reagent and elapsed time. Positive
finite-time stocks must remain representable. Reagent depletion is integrated
into the rate, never imposed by clipping a separate recovery curve.

The original `ConservedReaction::advance_inventory` publishes the final
extent, all species, stoichiometric amounts and explicit reactant accounts.
The result retains the original model, initial/final source points, elapsed
time, conversion, residual core radius and numerical-work/error report.
Failures publish no partial result or mutation. The separate
`constant_concentration_time` method evaluates the analytic maintained-bath
limit; it does not silently replenish a closed batch.

Eight G0/G1/G3 regression groups cover the independent finite-reagent film
solution, matched-capacity surface/diffusion/mixed integrals, R versus R²
scaling, series resistance, complete/excess/zero-stock cases, stoichiometric and
batch scaling, source-domain refusal, work limits/cancellation, and early time
with nearly stoichiometric excess. The runnable
`cargo run -p fs-material --example mineral_shrinking_core` uses the new
chalcopyrite/ferric endpoint with explicitly synthetic coefficients and prints
size-sensitive conversion, copper production and remaining ferric stock.

This is one irreversible, initially virgin, spherical solid population at
fixed volume, temperature and transport coefficients. It does not infer
activity/speciation, equilibrium, evolving porosity or film properties,
non-spherical geometry, reagent competition between particle populations,
temperature evolution or plant calibration. An already partially reacted
particle requires a history-aware model; repeated calls with a smaller stock
are not a continuation of its existing core geometry. The resistance form
follows Shiea, Crema and Macchi (2024), equation 3,
<https://doi.org/10.1016/j.ces.2024.119987>.

## Time-dependent grinding population balance

`comminution::kinetics::SelectionBreakageModel` binds the original
`SizeRedistributionCard` daughter mass fractions to one explicitly supplied
inverse-second selection rate per class. Its version-1 constitutive law is
`constant-selection-breakage`; parameter names are
`selection-rate:<class-id>`. A finite typed `grinding-contact-time` window
must start at zero. Source hashes, provenance, specimen identity and all
operating conditions remain attached. Positive-rate classes require strictly
finer daughters; a zero-rate class is an absorbing sink. Missing fines classes,
self daughters for an active class, invalid units and unsupported conditions
refuse. Selection rates are not inferred from mineral identity, Mohs hardness,
Bond work index or rotation speed.

The discrete mass balance is
`dm_i/dt = -S_i m_i + sum_j b_ij S_j m_j`.
For constant rates, `prepare` evaluates a nonnegative Poisson uniformization
with `q=max(S)` and `P=I+Q/q`. The explicit approximation caps the event
count at K and assigns the remaining probability to `P^K`. It does not drop
tail mass or renormalize the source columns. The reported nominal L1 bound is
twice the Poisson-tail bound; binary64 rounding and source uncertainty are
separate. The original redistribution owner must accept the resulting columns
without repair. Preparation refuses `q*t>64`, exhausted term/product budgets,
unrepresentable positive arithmetic or cancellation, without a partial result.

`KineticGrindingPass` retains both original cards, its exact typed preparation
point, elapsed contact time and numerical-work report. Its transition is
directly usable by the existing grinding-circuit API at that point.
`transform` uses the original `ParticleStream` owner to transport dry solids
and every assayed element congruently; carrier water is unchanged. This
constant-rate model supplies particle-size evolution, not mineral liberation,
a residence-time distribution, power draw, media wear or plant calibration.

Seven G0/G1/G3 regression groups cover independent equal-rate Erlang and
unequal-rate solutions, explicit truncation error, dry/element/water closure,
zero-time/zero-rate identity, split-interval consistency, source-domain
admission, cancellation and work bounds. Run the synthetic tin-bearing example
with `cargo run -p fs-material --example mineral_grinding_kinetics` to compare
30 and 120 seconds and inspect all balance residuals and tail bounds.

## Finite-batch reaction with reagent depletion

`reaction_batch::BilinearBatchReaction` binds an existing exactly conserved
equation, two explicit consumed-species roles and the complete original
constitutive card. For `a A + b B -> products`, the declared empirical law is
`dξ/dt = (k/V) (nA0-aξ) (nB0-bξ)`, with `k` in m³/(mol s), constant positive
batch volume and absolute temperature, and no inflow/outflow. Both stocks
affect progress continuously. The two kinetic orders are declared as one;
stoichiometric coefficients do not determine them. Mineral inventory divided
by batch volume is an empirical availability variable, not pure-solid activity
or exposed surface area. The source card must supply the coefficient and its
support; admission does not establish an ore calibration.

The card requires version 1, zero hidden state, original source hashes,
provenance and finite typed ranges for elapsed time, volume, temperature and
both actual concentrations. Time zero and the final time/concentrations must
lie inside the original rectangular domain. Monotone concentration depletion
makes endpoint admission sufficient for this law's complete path. All other
declared conditions remain fixed and mandatory. Caller conditions cannot
replace the owned coordinates; zero progress still requires admission.

The analytic equal/unequal-capacity solutions advance the complete inventory
through `ConservedReaction::advance_inventory`. Results retain the original
model/reaction, complete initial/final species amounts, gross stoichiometric
flows, extent, availability ceiling and consumed/remaining reactant accounts.
The equation conserves atoms and charge exactly; floating-point amounts are
not a numerical certificate. Unrepresentable progress or loss of a positive
finite-time reactant remainder refuses without mutation or clipping.

`mineral_leaching` demonstrates an explicitly synthetic dolomite/acid rate and
optionally consumes five real USGS material packs for standard heat at the
computed extent. It does not predict solution activities, equilibrium,
surface evolution, multi-reaction competition or temperature evolution.

## Recovery contact-time curves and conservative residence

`recovery::RecoveryKinetics` consumes the complete original constitutive card
for either a single first-order rate or a uniform population on `[0,kmax]`.
It requires the exact law/version, dimensioned ultimate-recovery/rate parameters,
zero internal state, original source hashes, provenance and a finite typed
elapsed-time validity window. Every additional declared condition remains
required. Actual contact time is separate from a hydraulic residence-time
distribution; an ore name never supplies a kinetic parameter. Stable
small-time evaluation retains positive representable response. The original
card and any fit evidence remain inspectable without a new calibration claim.

`recovery::residence::StirredTankCascade` advances a conserved constituent
through equal, full, ideal stirred tanks under piecewise-constant feed and
actual train flow. Concentration is constituent mass per slurry volume. The
caller supplies active mixed volume and initial inventories; bypass/dead space
is external and cannot silently alter the meaning of train flow. State persists
across flow changes and stops. Each interval independently integrates outlet
mass, reports its flow-weighted concentration, and checks inventory/inlet/outlet
closure without a plug. Failed steps publish no new state.

Positive source-scaled Poisson propagation has explicit 64-tank and 4096-exposure
limits. Longer intervals must be split. Premature probability underflow cannot
erase a representable physical response; truly unrepresentable positive paths
refuse. G0/G1/G3 tests cover analytic single-tank endpoint/mean, independent
Erlang pulse integration, stopped/variable flow, split intervals, numerical
extremes, work bounds and transactionality. `mineral_residence` demonstrates
the eight-tank transient bias with a separately declared downstream extraction
and independently displays both contact-time curves. No deconvolution, field
calibration, aqueous chemistry, segregation, adsorption, coupled reaction rate
or actual process recovery is inferred by this forward hydraulic model.

## Transient dissolved and undissolved constituent inventories

`recovery::reactive::ReactiveTankCascade` advances both forms of the same
constituent through 1..=64 equal, full, ideal tanks. Both remain in the tanks
and advect with the declared carrier flow. Local irreversible first-order
conversion transfers undissolved constituent to dissolved constituent;
dissolved outlet is a feed to downstream recovery, not already refined metal.
The original `LocalTransferCard` must explicitly authorize the transient local
law, supply its inverse-second rate and retain source identity, provenance and
typed validity conditions. Batch-contact and stationary-recovery declarations
are not promoted into transient laws. Rates may differ between tanks, while
each rate and its admitted external conditions remain fixed during a step.

Each piecewise-constant feed/flow interval independently propagates per-tank
inventories, both final-tank outlet integrals and the internal conversion
integral. The result exposes initial/final/inlet/outlet constituent accounts
and separate phase balance residuals. A stopped flow admits continued reaction
with zero inlet and outlet. Restart exports retained dissolved and undissolved
material according to their common hydraulic residence. Failed admission,
numerical range or conservation checks leave the original state unchanged.

Positive uniformization limits total combined exposure to 4096 and splits it
into substeps of at most eight; longer intervals must be split by the caller.
The floating-point residual is a check, not an interval certificate. Reagent
depletion, passivation, phase slip, precipitation, temperature evolution and
changes in carrier volume remain outside this dilute-constituent model.
Independent one-tank solutions, stopped-flow conversion, hydraulic zero-rate
equivalence, per-tank steady limits and pulse/restart split-interval checks
exercise the new slice. `mineral_residence` prints the two retained forms,
their effluents and internal conversion separately for an authored scenario.

## Continuous recovery from actual train flow

`recovery::continuous::ContinuousRecoveryCard` admits a separately sourced
stationary local first-order law, with a single rate or a uniform population
on `[0,kmax]`. Its law IDs differ from finite-time batch cards; it rejects an
`elapsed-time` domain. The source must bound typed positive mean residence
`V/Q`, integer tank count, and every additional declared condition. Active
mixed volume and actual train flow remain separate typed inputs. Zero flow
refuses a steady through-flow calculation rather than inventing recovery
during a stop.

For N equal full ideal tanks, single-rate survival is `(1+k*V/(N*Q))^(-N)`.
The uniform law analytically integrates that survival over the original rate
population. Rate classes keep their identities through the complete train;
the remaining population is never reset after a tank. This does not evaluate
batch recovery at the mean or silently extend finite test support into the
infinite residence-time tail. No automatic batch-card conversion exists.

`ContinuousTreatment` retains the original complete card, admitted conditions,
actual volume/flow and mean residence, recovered and surviving fractions,
contained-constituent feed/recovered/tail flows, and an independent conservation
residual. Its dilute-constituent assumption leaves carrier flow unchanged.
Stable small-response formulas and independent surviving-tail evaluation
refuse unrepresentable positive results. Eight focused tests cover independent
single-tank balances and ratewise train integration, source and geometry
admission, batch-card refusal, changed flow, stops and numerical extremes.

`continuous_recovery` compares one/eight tanks and two distributions at equal
mean rate. `mineral_slurry` directly supplies its real source-backed volume
flow to both continuous laws, with separately authored 17000 m³ active volume
and a 3 g/t dry-feed contained-metal assay. Increasing water changes residence
and recovery; unsupported residence prints `OutsideValidity` while retaining
the valid upstream slurry result. These examples establish arithmetic and
source flow integration, not measured calibration, transient reaction,
adsorption, reagent limitation, bulk concentrate grade or volume change.

## Distributed rates and conservative particle separation

`recovery::distributed::RecoveryCard` extends the source-qualified kinetic
families to gamma and explicit discrete populations. Batch and ideal plug
flow evaluate the population expectation of `1-exp(-kt)`. An equal-volume
steady CSTR cascade evaluates each discrete rate with remaining fraction
`(1+k*t/N)^(-N)` using TOTAL train residence. A mean rate never substitutes
for its distribution. Continuous-distribution CSTR requests explicitly refuse.
The immutable card preserves rate units, original weights, specimen/method,
provenance and typed source conditions. Its response separately retains
recovered, recoverable remaining and intrinsically unrecoverable fractions.
Scaled compensated sums preserve small representable rate populations.
These steady fractions do not advance the separate transient tank inventories.

`separation::ParticleStream` carries dry particle-class flows, representative
diameters, explicit solid-phase elemental flows and water. All classes declare
the same element axes; unknown assay is never zero. `split` requires exactly
one product/reject partition per homogeneous class and a separate water split.
Both products retain class and element identities, including zero-flow classes.
Reported grade uses actual product dry solids and is absent for zero solids.
Small representable product/reject tails are evaluated independently; arithmetic
that erases a positive flow refuses without changing the feed.

`ClassifierCard` consumes an original constitutive card for the corrected
Rosin–Rammler partition `c=1-exp(-ln(2)*(d/d50)^m)` and explicit water bypass
`b+(1-b)*c`. Version, parameters, dimensions, source hashes and source conditions
are required. Classification supplies actual particle diameter and stream solids
mass fraction; caller metadata cannot override them. Near-cut evaluation retains
adjacent representable sizes even with a sharp curve. Cut size and sharpness
are declared test parameters, not predictions from cyclone geometry.

Seven separation and nine distributed-recovery tests cover independent analytic
solutions, element/water/class conservation, reactor distinctions, small tails,
zero products, source refusals and ordering/scaling behavior. The runnable
`mineral_separation` example connects Bond dry throughput, authored size/assay
classes, classifier overflow and separate valuable/gangue kinetic populations.
It prints concentrate/tailings grades and complete one-pass balance residuals.
Underflow remains a separate outlet: recycle convergence, breakage, liberation,
aqueous inventories, entrainment and reagent response are separate capabilities.

## Source-bound mineral slurry streams

`slurry::SlurryStream` converts explicit dry constituent and carrier-liquid
mass flows into additive volume flow, pulp density, solids mass fraction and
solids volume fraction. Each `DensitySource` is a typed query against an
actual normalized material pack; the complete pack identity, original query,
evidence, uncertainty and verified usage receipt survive all flow operations.
Reported density ranges require an explicitly selected endpoint. Grain-density
references do not become bulk-ore measurements or common-temperature evidence.

The caller explicitly selects additive constituent volumes with no entrained
gas. `dilute_to_mass_fraction` adds the original carrier while preserving dry
constituent flows. A concentration target reports required liquid removal;
it never silently clips the request. Nonselective splitting and canonical
mixing conserve each component to binary64 rounding. Identical component names
with different source/query contexts refuse to merge. Zero-flow streams retain
zero extensive rates and `None` for undefined density/composition. Distinct
mass- and volume-fraction types prevent accidental basis substitution.

There are at most 128 solids or mixing inlets. Fixed sorted/compensated sums,
finite positive density checks and explicit overflow/underflow refusals keep
these synchronous operations bounded. They own no evolving plant inventory.
G0/G1/G3 tests include a real pack/compiler journey from Bond dry throughput
through magnetite/pyrite/dolomite/water densities, dilution and split/recombine.
The runnable `mineral_slurry` example exposes the complete source-backed flow
calculation. No dissolution, reaction, thermal equilibration, volume contraction,
settling, rheology, liberation, separation efficiency or field calibration is
inferred; source uncertainty is retained without inventing a mixture interval.

## Selective separation of source-bound slurry constituents

`slurry::separation::SelectiveSeparator` applies one explicit concentrate
partition to each exact solid identity and a separate carrier-liquid partition.
The existing `SlurryStream` outlets retain every original density source,
typed query, evidence and usage receipt. Complementary concentrate/tailings
flows preserve each constituent; no new material, dissolution or liberation
is inferred. Provenance records the prescribed fractions without promoting
them to calibrated flotation, gravity or cyclone performance.

The result reports per-component flows, computed recovery and numerical closure,
dry-component fractions in all three streams, carrier accounting, and total
dry-mass yield. Fractions refer to the supplied identity set, which may include
mineral-by-size classes; they are not elemental assays. Missing, extra or
duplicate partitions refuse, including missing declarations for zero-flow IDs.
Zero-feed recoveries and zero-solids grades are absent. Every fraction lies in
[0,1]; positive unrepresentable outlets or metrics refuse atomically.

Seven focused tests cover independent balances, two-stage rougher/cleaner
accounting, density-source preservation, nonselective equivalence, distinct
size-class routing, canonical ordering and undefined/numerical boundaries.
`mineral_partition_circuit` runs that two-stage example using actual compiled
magnetite, pyrite, dolomite and water packs, combines both tailings streams,
and prints dry grade, recovery, volume flow and component/water residuals.
It complements the separate size/element `ParticleStream` classifier; no
automatic conversion between their distinct composition models is implied.

## Yield-stress slurry rheology and pipe pressure

`slurry::rheology::HerschelBulkleyCard` retains an explicit yield stress,
reference excess stress at a reference shear rate, flow index and original
specimen/test provenance. This reference-rate form keeps stress in Pa even
for noninteger exponents. Signed stress and inverse shear rate use the same
finite temperature, solids MASS fraction and absolute-rate domain. Zero rate
returns the complete static yield interval; positive shearing reports
nonnegative dissipation. `SlurryConditions::from_stream` transfers the existing
stream's mass fraction and requires an explicit temperature; density metadata
cannot silently establish thermal equilibrium or a rheological fit.

`LaminarCircularPipe` integrates fully developed no-slip Herschel–Bulkley flow
analytically, retaining the unsheared plug and stationary yield threshold.
It returns signed flow, pressure drop, wall stress/rate, plug radius and
hydraulic dissipation. Pressure inversion uses bounded bisection and an explicit
volume-flow tolerance inside both the source domain and caller-declared
laminar wall-rate ceiling. The ceiling is an applicability declaration, not
a predicted transition. Overflow, positive underflow and unsupported conditions
refuse; no turbulent, settling, entrance, elevation or pump-efficiency correction
is inferred.

`profile_at_pressure` retains the same accepted flow and borrowed source to
evaluate axial velocity, radial traction, velocity gradient and local irreversible
heat without allocation or quadrature. A moving or stationary plug can carry
stress with zero local heat. Complementary radial coordinates preserve near-axis
shear and thin yielded annuli; the no-slip wall has zero velocity and finite
viscous heating. This frozen-state profile adds no heat/advection evolution.

The rheology tests check Poiseuille and Buckingham–Reiner limits,
independent radial quadrature for a noninteger exponent, signed reversal,
inverse residuals, yield behavior, source limits and numerical refusal. Radial
velocity and heat independently integrate to the original flow and pressure work,
including thinning and thickening laws, with near-axis and wall regressions. The
stream adapter has a separate real-pack mass/volume-basis regression.
`cargo run -p fs-material --example mineral_pipe -- 0.005 0.6` exposes an
authored pressure-sizing example; its parameters are not a plant calibration.

## Conserved size redistribution and particle-stream mixing

`comminution::breakage::SizeRedistributionCard` retains a complete source-declared
per-pass mass-transfer matrix over at most 128 named size classes. Every column
closes to unit mass within its documented binary64 allowance without
normalization. Positive off-diagonal transfers must enter strictly smaller
diameters; diagonal entries retain the original class. Daughter material inherits
each parent's explicit elemental composition. This is congruent redistribution,
not inferred mineral liberation, selection-rate kinetics or chemical reaction.
Water passes unchanged. Original context, typed validity, columns and independent
retained/fragmented fractions remain available. Exact feed layouts and original
source coordinates are required; positive transfer underflow refuses atomically.

`ParticleStream::mix` conserves the same explicit class/element/water layout,
rejecting incompatible diameters or missing assays. Magnitude-ordered compensated
sums make inlet order irrelevant. `scaled` preserves layout even at zero flow.
`discrete_passing_size` returns an existing representative diameter from the
dry-solids CDF, grouping equal sizes without inventing within-bin interpolation.
Six redistribution and six stream tests cover known Cu/Fe/water balances,
identity, scaling/order invariance, tiny retained tails, numerical refusal,
typed applicability and discrete quantile boundaries.

## Closed grinding and classifier recycle

`comminution::circuit::ClosedGrindingCircuit` solves `x = f + diag(r) B x`
for actual mill feed, with fresh feed f, the original per-pass redistribution B
and the classifier's retained fractions r. Strictly smaller daughter classes
allow a direct decreasing-size triangular solve: work is quadratic in at most
128 declared classes and linear in the explicit element accounts. Positive
product and fragmentation exits avoid cancellation near closed loops, while
fused coefficient arithmetic does not erase representable physical transfers.
There is no fitted circulating load or iterative convergence threshold.

Water has its own recycle balance. The mill's optional solids-mass-fraction and
mill-dry-solids-flow validity axes are bound from the solved mill feed; callers
cannot override them. Its original transform and the classifier are then
evaluated on actual streams and conditions. Admission therefore never treats a
guessed initial feed as the final operating point. Recomputed class/element/water
balances and external product-versus-fresh accounts must all satisfy the
caller's absolute/relative allowance. A separate maximum circulating load is
required. A fed class with no represented exit refuses; an unforced closed class
retains the explicitly selected zero-circulation state. Failed solves publish
no streams and never change the fresh input or source cards.

Seven G0/G1/G3 tests cover independent rational balances, geometric-pass limits,
very large circulation with retained small product tails, trapped classes,
actual-domain/coordinate refusals, scaling and positive-flow numerical limits.
The `mineral_circuit` example reports all streams, discrete P80, circulating
load and residuals, with an independent Bond estimate on FRESH circuit throughput.
It does not multiply circuit-specific Bond energy by recirculating mill tonnage.
Fixed per-pass coefficients remain a source-declared approximation, without
liberation, selection-rate kinetics, transient hold-up or a process calibration.

## Source-specific comminution references

`comminution::source::resolve_bond_ball_mill` resolves four original scalar
material claims for one specimen's work index, closing screen, F80 and P80.
The claims must share one observation and the same source support. The
metric-tonne and short-ton property names remain distinct after normalization
to J/kg; missing/wrong units, unbacked values, mixed specimens, interpolation
and unsupported source coordinates refuse. `BondGrindingCard` retains the
original pack hash, query, claims, observation, evidence and verified receipts
through ordinary cloning and its existing physical consumers. Reported test
sizes remain separate from the explicitly supplied operational size domain;
neither scope is expanded by resolution. The `bond_grinding_reference` example
uses the five historical Joanna composites from the real material pack.
Source-backed test values do not make its energy/throughput outputs measured
plant performance. Unit regressions and `xtask`'s `mining_bond_sources` exercise
source preservation, independent energy/throughput answers and refusals.

`comminution::drop_weight::resolve_jk_drop_weight` resolves the independent
A-fraction and inverse-specific-energy b parameters of the JKDWT impact law
from one scalar MatDB observation with identical source support. Original
claims, evidence, uncertainty, observation, query, pack identity and verified
receipts survive in the returned card. A must be in (0,1] and b positive;
percent-valued A must be normalized before resolution. Missing or mixed
parameters and unsupported source conditions refuse. `response` evaluates
`t10 = A*(1-exp(-b*Ecs))` in coherent SI and preserves the complementary
retained fraction; `energy_for_t10` inverts reachable targets below A. Both
honor the explicitly supplied finite nonnegative energy interval, which does
not establish experimental support. No PSD, abrasion ta, selection rate,
liberation or SAG power is inferred. The `bueno-jkdwt-2021` material pack supplies
five published fitted specimens; `drop_weight_reference` and the ordinary
`mining_drop_weight_sources` compiler integration exercise those source rows.

`comminution` supplies a typed, uncorrected Bond grinding-energy closure,
dry-solids throughput from net pinion power, and ideal single-ball critical
speed. Work-index construction distinguishes kWh per dry metric tonne from
kWh per dry US short ton and retains the original denominator. The immutable
card records specimen identity, test/preparation procedure, closing-screen
aperture and source provenance. Its explicit F80/P80 validity window is a
caller declaration, not proof of site calibration. This baseline refuses
P80 below 70 micrometres and requires actual size reduction.

Energy outputs are J/kg dry solids; throughput is kg/s dry solids with explicit
utilization; critical speed is revolutions/second from the inside-liner and
ball diameters and caller-supplied gravity. Motor nameplate power cannot stand
in for net pinion power. Finite/positive admission, out-of-domain inputs and
unrepresentable arithmetic return `ComminutionError`; nothing is clamped.
Fixed scalar operation order has no mutable state or IO. Blend work is linear
in the number of explicitly declared feeds.

`AdditiveGrindingBlend` explicitly selects the approximation
`E_blend = sum(m_i * E_i) / sum(m_i)` for one common circuit and target P80.
Every `GrindingFeed` retains its own dry mass, F80 and original
`BondGrindingCard`; every source admits the requested product independently.
Admission requires positive masses and a nonempty common product window,
with strict size reduction for each feed. The weighted energy is converted
to throughput only after summation, rather than averaging work indices,
F80 values or individual throughputs. Scaling the extensive factors avoids
overflow of a finite mean; unrepresentable nonzero contributions refuse.
This opt-in rule does not establish physical additivity of an actual ore blend
or estimate interaction, liberation or classification effects.

`BondGrindingCard::product_size_for_specific_energy` inverts the original
Bond relation for P80; `product_size_for_throughput` first obtains energy
from explicit pinion power, utilization and dry-solids flow. Both retain the
same source/test context and require the resulting P80 to lie in the original
size window, with actual reduction. Positive flow with zero utilization
refuses. Exact forward endpoint energies round-trip to their declared sizes;
values outside that energy window refuse rather than clamping. Rationalized
small reductions and scaled power/flow arithmetic avoid loss through
intermediate cancellation, overflow or underflow. No recovery or circuit
efficiency is inferred by inversion.

G0/G1/G3 tests cover hand-calculated energy/tonnage, the published GMG 2021
stage-energy example, metric/short-ton equivalence, monotonicity, one-ULP
size reductions, validity refusal and critical-speed force balance/scaling.
Blend tests independently reproduce 8 and 10 kWh/t feeds yielding 9.5 kWh/t
at 25%/75% dry mass, reject incompatible source windows, and preserve energy
under dry-mass scaling, feed permutation and lot splitting.
Inverse tests recover 100 micrometres from 43200 J/kg or 12 MW at 0.9
utilization and 250 kg/s, cover domain endpoints and one-ULP reductions,
and reject unsupported grind targets, invalid inputs and downtime.
These establish reference-formula behavior, not site predictive validity.
There is no SAG/AG, circuit correction, slurry rheology, liberation, recovery,
wear, DEM trajectory or recommended operating-speed claim. Source and
calibration uncertainty are not converted into a certified prediction band.

## Explicit phase coordinates in spatial conductivity

`SpatialConductorMaterial` retains an optional original liquid-mass-fraction
coordinate and a distinct stationary-Ohmic source admission.
`ResolvedOrthotropicConductor3::try_from_phase_material` admits three exact
principal sources with that same existing typed coordinate in [0,1]. Refresh
preserves all keys, independent conductivity/resistivity basis, complete source
bundle, phase role and fixed principal axes. No mixture law or missing principal
coefficient is inferred. The original static-source admission used by induction
is unchanged; a phase-dependent snapshot does not silently widen that owner.
The tetrahedral current/finite-enthalpy regression supplies the actual original
caloric fraction on every endpoint trial, including latent-temperature plateaus.

## Original phase response in the Ohmic constitutive node

`conductor::graph::MaterialOhmicNode::try_new_with_liquid_mass_fraction`
extends the same native uniform conductor with an explicit original fraction
axis. Its three graph inputs are temperature change, signed current and liquid
MASS fraction change, referenced to the original complete material query. The
two-input constructor and its sensible/free-expansion arithmetic remain intact;
the three-input declaration uses law version 2 and fixed geometry.

The same pinned source supplies resistance and nonnegative Joule heat. Constant
scalars, supported linear T/f curves and canonical triangulated T/f surfaces
provide both analytic derivatives. Conductivity derivatives pass through the
native reciprocal relation rather than a fitted resistivity curve. Surfaces
use `fs-matdb::PreparedSurfaceResponse`, retaining original source cells and
their declared one-sided tangents. Every runtime value still resolves the
original complete bundle and actual applicability coordinates. Independent
T/f bounds never imply that an unmeasured source-surface corner or hole is
supported. Unsupported states return no physical candidate.

This stationary scalar DC law supplies no phase mixing rule, density/geometry
evolution, phase kinetics, interfacial electrical impedance or additional
latent heat. The thermal owner supplies actual fraction and caloric storage.
Native graph/enthalpy regressions exercise scalar/curve/surface source swaps,
both source bases, both axis orders, latent current/voltage duality, continued
heating, exact constant-source compatibility and transactional refusal. Source
and parser checks do not establish native execution or measured qualification.

## Variable-amplitude critical-plane strain-life accumulation

`fatigue::critical_plane::miner::CriticalPlaneMiner` projects each accepted
stress/total-mechanical-strain tensor onto every explicitly declared fixed plane
and delegates each signed normal history to the existing bounded
rainflow/Palmgren-Miner strain-life owner. The selected plane maximizes total
accumulated Miner damage over the complete observation; exact damage ties use
maximum strain amplitude and then declared plane order. No full tensor history
is retained.

Uncorrected counting accepts arbitrary finite simultaneous normal stress while
using the projected strain history. Morrow and Smith-Watson-Topper retain the
existing paired-history requirement: stress/strain reversals must remain
coincident, and phase lag or relaxation that would hide stress extrema refuses
instead of inventing a scalar cycle correction. Pushes and finalization are
transactional across the whole plane bank. Sample count and samples-times-planes
projection work are bounded before publication.

This is a finite-plane, variable-amplitude multiaxial fatigue estimate. It does
not optimize over the continuum sphere, infer extrema between observations,
model sequence interaction beyond linear Miner summation, establish cyclic
stabilization, evolve stiffness/cracks, or qualify component life. Calibration
scope, temperature, process, surface/environment and source evidence retain
their existing strain-life meaning.

## Fatemi–Socie shear-strain life on native tensor cycles

`fatigue::critical_plane::fatemi_socie` evaluates the original
maximum-shear-strain-plane Fatemi–Socie criterion from paired physical stress
and total mechanical strain tensors. It selects the greatest engineering shear
strain amplitude, then the greatest normal stress for exact amplitude ties.
The normal-stress maximum is taken over the full cycle; it need not occur at
the shear extremum. Selecting the maximum combined parameter would be a
different model and is not this plane-selection convention.

Engineering shear is twice the tangential projection of the tensor strain.
Its two-dimensional path uses the same bounded enclosing-circle kernel as
Findley, retaining rotation and phase rather than reducing the history to an
unsigned equivalent strain. The positive calibrated relation is
`gamma_a * (1 + k * sigma_n,max / sigma_y0)
= tau_f'/G * (2N)^b_gamma + gamma_f' * (2N)^c_gamma`.
The source declares initial yield normalization and a named normal-stress
convention: Signed refuses a nonpositive multiplier; TensionOnly uses the
tensile part of the maximum normal stress. Neither imports Findley clipping.

The shear calibration has distinct coefficients and a strict finite cycle
domain. Only the private safeguarded scalar inverse is shared with axial CMB;
no axial coefficient or Young's modulus is relabeled as shear data.
`source::resolve_isotropic_fatemi_socie` retains eight original claims:
E, nu, initial yield, four torsion coefficients and k. It explicitly derives
G from isotropic E/nu, preserves original claim selections and uncertainties,
normalizes an explicitly declared N or 2N power basis, and intersects the
torsional fit's life support with the separate k fit. Original intervals must
be valid before intersection. Convention-specific k keys prevent an implicit
signed/tension-only source switch. Every parameter supports the actual typed
absolute temperature; fitted coefficients require the original reversed
torsional baseline and declared finite life support.

All six stress and strain components must close within explicit tolerances.
A closed representative cycle is caller-declared; closure alone does not prove
cyclic stabilization, resolve internal variable-amplitude cycles or count
repetitions. Results expose full-cycle life, selected plane, shear amplitude,
normal extrema, correction, closure defects and numerical work. Sample,
plane, projection and circle limits and cancellation bound computation.
The retained history supports inexpensive append rollback. Accumulation is
available only through the explicit fixed-plane mode below. No stiffness loss,
crack growth, inferred endurance limit or validation of component service life
is claimed.

The plane and strain conventions follow the
[original Fatemi–Socie formulation](https://doi.org/10.1111/j.1460-2695.1988.tb01169.x)
and [strain-based fatigue model reference](https://doc.comsol.com/6.4/doc/com.comsol.help.fatigue/fatigue_ug_sme.4.23.html).

### Explicit fixed-plane accumulation across resolved shear cycles

`FatemiSocieCycle::enable_miner(maximum_cycles)` opts into bounded FS/Miner
accumulation before the first loading excursion, with zero or one observation.
It can be enabled only once, and the positive complete-cycle limit survives
window restarts. Every successful finalization represents exactly ONE
caller-declared resolved complete cycle. A window containing several reversals
is not automatically decomposed into multiaxial cycles.

For every original material plane p and accepted cycle i, the same vector-path
amplitude and whole-window normal-stress maximum produce its calibrated life
N(p,i). The retained rule is `D(p) = sum_i 1/N(p,i)`; the summary reports
`max_p D(p)`. It never adds the winners from different planes into one damage
history. This maximum-accumulated-damage approximation is a separate model
choice from the original maximum-shear-plane representative-cycle result,
which remains unchanged and separately inspectable.

`miner_summary()` exposes completed and maximum cycle counts, plane count,
selected original normal and accumulated damage. `miner_plane_damage()`
exposes every plane's total in the original normal order. Exact damage ties
select the first plane. Sums use compensation and refuse nonfinite arithmetic.
Zero shear contributes zero; EVERY positive candidate-plane parameter must
remain inside the original finite life calibration. Thus a weak off-axis
plane below the fit can refuse accumulation even when the representative
maximum-shear plane is supported. No endurance limit or extrapolated life is
silently substituted.

The projection and enclosing-circle work is shared with the existing result;
only the opted-in mode inverts each positive plane's life. Additional retained
state is O(number of planes), independent of previous cycle sample counts.
Finalization stages all totals and publishes after the last domain, work and
cancellation check. Repeated finalization is idempotent. `fresh()` carries
completed totals while discarding an unfinished window without credit;
`discard_after()` retracts the current window's provisional contribution
without removing previously carried cycles. Neither operation resets native
constitutive history or supplies missing loading exposure.

This empirical linear damage sum does not model load-order interactions,
nonproportional hardening, creep, changing material conditions, or fatigue
feedback into stiffness and heat. Cycle segmentation, sampling resolution,
conditioning and calibration applicability remain explicit caller obligations.
Six focused numerical regressions cover changing damaging planes, varying
amplitudes and normal peaks, rotating shear, chunking, prefix rollback,
late-plane refusal, cancellation and cycle limits. They validate implementation
identities rather than component service life.

## Shear-path Findley fatigue on declared elastic stress windows

`fatigue::critical_plane::findley` evaluates the original Findley parameter
`P = tau_a + k * sigma_n,max` on explicitly declared candidate planes. The
shear amplitude is the radius of the minimum enclosing circle of the sampled
two-dimensional shear-vector trajectory on each plane. Retaining both shear
components distinguishes rotating shear from constant shear, including a
rotating path whose unsigned shear magnitude never changes. Plane selection
maximizes the combined Findley parameter, with deterministic tie handling.

The result is a nonnegative usage factor relative to an explicitly supplied
Findley limit factor. It is an empirical limit comparison, not a cycles-to-failure
estimate, Miner damage, a crack-growth rate or a certified endurance claim.
The full tensor stress window must close within declared stress tolerances and
contain actual cyclic variation. Candidate-plane and time resolution remain
explicit numerical limitations; finite plane sampling cannot certify the
continuum maximum. Sample, plane, projection-work and circle-work limits bound
the calculation, and cancellation/refusal leaves accepted observations intact.

`findley::source::resolve_findley` resolves the distinct original sensitivity
and limit-factor properties from one immutable material card. Their selected
claims retain source receipts, uncertainty, observations and validity domains;
both parameters must belong to the same declared fit and observation population.
Each must admit the actual semantic absolute temperature in kelvin. Fixed
scalars and source-supported temperature curves are admitted without inferring
fatigue parameters from static strength or a single axial S-N datum. The
calibration exposes its original card, material state and typed query for native
consumer binding. Scope labels declare applicability; they do not prove transfer
to another geometry, process, environment or load history.

The shear-path convention follows the
[COMSOL critical-plane methods reference](https://doc.comsol.com/6.3/doc/com.comsol.help.fatigue/fatigue_ug_sme.4.18.html).
The numerical tests use analytical synthetic histories and do not constitute
experimental calibration or validation of a material.

The constitutive-law and authority-separated identifiability kernel: materials
as mathematical objects — calibration domains, CONSISTENT tangent operators,
thermodynamic guardrails, hysteresis, uncertainty — plus source-closed,
multi-case questions whose execution and assessment identities stay distinct.

Ambition tags: elastic/hyperelastic/J2/RC-fiber laws + calibration [S];
multi-case identifiability authority schema [F]. Ogden and tension-only
HGO/GOH are [S] constitutive laws; specimen/tissue calibration is separate.

## One magnetic history shared by multiple windings

`field_storage::hysteresis::shared::SharedHystereticCore` couples one original
`UniformHystereticCore` to one through sixteen physical windings. Each winding
declares its signed turns ratio to the reference core and positive physical
leakage inductance. Total linkage is `lambda_i = D_i i_i + r_i lambda_ref`;
the common magnetomotive balance is `I_ref = sum(r_i i_i)`. Eliminating these
currents gives a positive linear load line on the original scalar material
inverse. No independent winding histories, guessed coupling coefficient,
numerical leakage, or fitted core-loss resistance is introduced.

The source-owned polarizations, reversible core energy and irreversible
pinning heat occur once. A candidate integrates the exact piecewise-linear
material path along a straight chord in all total linkages. Its work-conjugate
mean winding currents use integration by parts, retaining pinning and saturation
crossings where endpoint current averages would give incorrect work. Independently
checked winding work equals the stable change of leakage and core energy plus
the original pinning heat. A zero common drive advances only the leakage modes.
The full row-major active-set inductance matrix is `D + L_ref r r^T`, including
signed mutual entries and the original one-sided tangent at thresholds.

All candidates retain the original winding geometry and accepted core memory.
Prepared remanence must come from an explicit core history with compatible total
linkages; a source query can refresh only the original admitted constant fit.
Incomplete/nonfinite states, incompatible magnetic balance, unsupported source
conditions and unrepresentable work refuse without changing accepted state.
Geometry and constitutive coefficients remain fixed. Core reversals require
separate increments; the result is a rate-independent path law, not a complex
AC hysteresis model or a thermal/motion integrator. Copper and eddy-current heat
belong to their existing circuit/field owners and must not duplicate pinning heat.

Focused authored-fixture tests independently integrate exact active-set spans,
checking multiple winding polarities/orderings, threshold crossings, reversal,
prepared remanence, differential leakage modes, full inductance and unchanged
state on refusal. These verify the declared constitutive construction, not
experimental material qualification.

### Source-backed series gaps and reciprocal magnetic work

`SharedHystereticCore::with_gap` extends the fixed-geometry path above with one
physical series gap. Its original `ResolvedInductor` supplies constant static
permeability, area, absolute reference turns and preparation separation; each
accepted state retains the actual gap within caller-declared positive open
limits. The original medium source/query stays fixed. The common balance is
`sum(r_i i_i) = I_core + R(g) lambda_ref`, with
`R(g) = g/(mu A N_ref^2)`. The full mutual tangent includes this reluctance;
core memory, physical leakage and gap energy remain distinct original stores.

`candidate_linkages_and_gap` advances the total linkage/gap chord without
advancing any accepted memory. The common flux stays monotone on that chord.
Original play thresholds partition the path into affine material spans; exact
linear-over-linear flux moments provide both mean winding currents and mean
attractive force `-R' mean(lambda_ref^2)/2`. Positive balanced susceptibility
sums and compensated moment reductions bound work by the original branch count,
without field sampling or a fitted midpoint force. Electrical work independently
equals total stored change, one original pinning heat, and work delivered to
mechanics. Constant-flux chords retain the original polarizations exactly.

Source/geometry mismatches, unresolved thresholds and unrepresentable work
refuse without publication. Core-temperature queries do not silently change
the separately retained gap medium. Native mechanics owns inertia, supports
and trajectory closure; fringing, contact, changing area and thermal expansion
are outside this uniform-gap reduction. G0/G1 fixtures in
`tests/shared_hysteretic_gap.rs` compare full linear mutual inductance, fixed-
linkage virtual work, negative reference polarity, pinning/saturation/reversal
and exact constant-flux history against independent bisection/quadrature.


## Layer-resolved thermal expansion of elastic laminates

`laminate::thermal::PreparedLaminateThermal` shares the original immutable
`LaminateSection`, its ordered source layers, physical faces, reference surface,
rotated stiffness and inertia. Each layer receives its own integrated engineering
free strain, uniform or affine between its faces in its principal 1/2 plane.
Adapters lower the existing isotropic and orthotropic temperature-path owners;
callers retain those original paths and source receipts. No CTE averaging,
endpoint-CTE shortcut or inferred thickness restraint enters this plane-stress law.

Two-point thickness integration prepares thermal membrane force and signed
bending moment, and evaluates physical layer stress and elastic energy without
allocation, source lookup or trigonometry. Energy is integrated from the actual
elastic strain before squaring; a freely bending incompatible stack retains its
residual self-stress energy even when net force and moment vanish. Cancellation
is polled per layer before publication. Invalid profiles, overflow and Prony
layers refuse: frequency-reference stiffness is not a thermal-history law.

The numerical fixtures in `tests/laminate_thermal.rs` cover equal and unequal
bimetals, compatible thickness gradients, reference-surface shifts, rotated
orthotropy, work derivatives, source integration and refusal/cancellation.
This is prescribed-temperature infinitesimal-strain elasticity with fixed
section properties; heat evolution, viscoelastic relaxation, interlaminar
traction, debonding and measured laminate validation remain outside this lane.

### Temperature-resolved elastic laminate sources

`laminate::thermal::source::LaminateThermalSource` binds each original immutable
card, full reference context, source selection, layer geometry and reference
density once. Each current layer temperature refreshes the pinned elastic
claims and integrates the pinned CTE claims from that layer's reference point.
The returned section contains the current stiffness and must be used in native
shell assembly. Reference thickness, orientation, density and all mass moments
remain unchanged; density is not queried at the new temperature.

Isotropic E/nu, E/G and K/G schemas and measured orthotropic plane-12 elasticity
are supported independently of isotropic or principal-direction expansion.
The orthotropic plane requires both transverse shear moduli without inventing
E3 or out-of-plane Poisson ratios. Original reference/current receipts and
expansion paths remain accessible. An explicit thermal-cell subdivision keeps
one temperature per sublamina; no average temperature or CTE is substituted.
This is held-temperature elastic preparation, with no implied mechanical
history, heat/work feedback, or solved-strain calibration admission.

`tests/laminate_thermal_source.rs` checks independent stiffness/resultant
formulas, nonlinear CTE integration, fixed mass, heating/cooling, claim pinning,
directional response and source/domain/cancellation refusals.

### Source-pinned plane-stress thermoelastic potentials

`laminate::thermal::thermoelastic::AffineLaminaThermoelastic` differentiates the
original affine E/nu, E/G, K/G or orthotropic plane-12 sources before evaluating
their rational reduced stiffness. Current source receipts and the original
integrated expansion path remain authoritative. Exact first/second temperature
derivatives supply work-conjugate stress, entropy `S = -Psi_T`, internal energy
`U = Psi + T*S`, fixed-total-strain capacity and the mixed stress/temperature
derivative. No unmeasured thickness-normal elastic constants are invented.

The declared principal engineering-component envelope covers both total and
elastic strain. Outward interval arithmetic bounds the elastic capacity over
the complete source interval and refuses an unresolved rational pole; a thermal
consumer must add its original caloric capacity and admit a positive total.
Exact zero source slopes and equality of every original pinned support-knot
value identify temperature-independent elasticity, without mistaking a tiny
interior source variation for a constant factorization.
In-plane evaluation allocates nothing and performs no source query or rotation.
Local transverse shear uses RAW G13/G23: native whole-stack Whitney corrections
and their temperature derivatives are separate section physics, not a universal
per-layer 5/6 factor. These potentials neither advance heat nor reset memory.

`tests/laminate_thermoelastic.rs` uses independent compliance/eigenchannel
oracles for all four schemas, work/Maxwell and thermal derivatives, capacity
bounds, source pinning and domain/kink refusals. Targeted Cargo metadata checks
pass; numerical execution and experimental validation are separate evidence.

### Coupled elastic shear from declared cylindrical-bending probes

`LaminateSection::try_new_equilibrated_shear` prepares a full transverse shear
stiffness from two explicit orthogonal cylindrical-bending probes. Complete
rotated ply stiffness, membrane-bending coupling and transverse compliance
enter the equilibrium traction fields and their cross-energy. Resultants are
normalized before positive three-point Gauss integration of the quartic energy.
Preparation is linear in layer count; existing shell consumers use the cached
2x2 stiffness with no extra constitutive integration in their element loops.

`shear_model()` retains the selected closure and physical bending-axis angle.
Coordinate covariance requires co-rotating those axes with the laminate. This
is an explicit elastic reduction, not a load-independent 3-D warping or free-edge
solution. Prony layers refuse. The original directional Whitney constructor
and its thermal derivatives remain unchanged; derivatives of this new closure
are not supplied. Focused regressions include homogeneous 5Gh/6, an independent
exact-rational off-axis cross-energy result, frame/reference invariance and
actual MITC shell work for positive and negative mixed shear.

`LaminateShearRecovery` retains the exact section through `Arc` and caches the
same kernel's normalized quadratic traction maps once per original ply. Local
queries recover continuous bonded-interface shear traction, zero exterior-face
traction to roundoff, and each ply's actual engineering strain through its
uncorrected transverse compliance. Full off-axis cross terms remain present.
Integrated local complementary energy equals the original coupled H energy;
the section shear coordinate is not mislabeled as a constant physical ply
strain. Preparation is O(layers); a point query allocates nothing and is O(1).
Directional Whitney sections refuse this distinct recovery model. Reference
and material/probe rotations retain their existing meaning. No strength limit,
delamination, transverse normal stress or general free-edge solution is added.
`tests/laminate_shear_recovery.rs` checks the homogeneous parabola, layered
resultants/work, interface continuity, compliance jumps and frame invariance.

### Nonlocal thermal derivatives of native laminate shear

`laminate::thermal::shear::ResolvedLaminateShearThermoelastic` differentiates
the SAME directional Whitney shear operator with original source-pinned Q/G
derivatives, fixed reference geometry and `fs-ad`. Neutral-axis motion,
through-thickness energy integrals and off-axis congruence remain in the chain
rule. Native coefficient agreement is checked; its numerical fallback refuses
instead of acquiring invented derivatives. The one stack potential supplies
per-cell entropy conjugates, `U=Psi+sum(T_i*S_i)`, directional stiffness
derivatives and the entropy Hessian action, including mixed layer derivatives.
There is no arbitrary allocation of shear internal energy among cells.

Preparation and Hessian actions use O(n²) work and O(n) storage; directional
derivatives use O(n) work, and point evaluation uses cached gradients. Its
directional native Whitney approximation does not establish general off-axis
multilayer covariance or a fully coupled 3-D warping model.
`tests/laminate_thermal_shear.rs` covers the homogeneous closed form, unequal
temperature derivatives, mixed symmetry, reference/axis transformations and
source/domain/cancellation/fallback refusals. Numerical execution and physical
validation require their own evidence.

`shear::bounds::LaminateShearCapacityBound` extends this potential to the complete
Cartesian product of the original source-temperature intervals. Outward source
jets and analytic derivative-norm propagation bound each absolute row sum of
the shear entropy Hessian, per reference area, over the declared principal
strain envelope. Positive square integrals avoid cancellation in the Whitney
energy enclosure; interval poles, unresolved positivity and unrepresentable
bounds refuse. The one-layer homogeneous reduction retains the exact 5Gh/6
identity. O(n²) preparation and O(n) storage avoid a dense thermal Hessian.
The caller must combine this bound with its actual caloric and in-plane entropy
derivatives; it is not itself a heat step or a discrete passivity certificate.

`temperature_increment` gives the finite change of the same source operator
one original temperature coordinate at a time. It averages forward and reverse
spatial-order sweeps, using cached endpoint coefficient records without new
source queries or AD preparations. The coordinate matrices telescope to H1-H0;
reversing the endpoints reverses every contribution, and unchanged coordinates
receive exact zero. Endpoints must retain the exact same immutable source
owner, geometry, original affine intervals and strain envelopes. Roundoff is
bounded and reported, never distributed into a physical heat source.

For endpoint shear strains g0/g1 and reference area A, the per-coordinate
temperature work is `A*g0^T*delta_H_i*g1/2`, and reversible heat is its negative
minus `delta(T_i*S_i)`. Summing these transfers gives mechanical shear work
minus the change in the ONE stack internal energy. This symmetric coordinate
path is a declared finite-step approximation requiring refinement; it is not
the exact straight-path entropy integral or invariant under arbitrary layer
reordering. The material carrier does not advance or replace caloric storage.
`tests/laminate_thermal_shear_bounds.rs` and
`tests/laminate_thermal_shear_increment.rs` cover closed forms, mixed-temperature
enclosures, work telescoping/reversal, consistency and atomic refusals. Focused
compiler checks and independent analytical references are distinct from Rust
numerical execution and measured material validation.

### Thermal-strain histories of compatible viscoelastic laminates

`laminate::thermal::history::LaminateThermalHistory` retains original full-3D
TensorProny histories at two thickness points per layer. It admits one common
plane-stress strain extension shared by the equilibrium tensor, every Maxwell
branch and optional Kelvin viscosity. Fixed-Poisson isotropic spectra and
compatible anisotropic spectra are supported; differing relaxation maps refuse
explicitly. No independently condensed branch law or frequency-reference
stiffness is substituted for the original material memory.

Prescribed membrane/curvature and integrated layer free-strain profiles advance
as exact affine ramps through the native history owner. Mean stress separates
actual mechanical work from the eigenstrain-work port; equilibrium/branch
storage and independently integrated nonnegative dissipation remain distinct.
Layer face stresses are recovered from the exact affine field. Complete-state
proposal/commit semantics preserve accepted history on cancellation, source
domain exit or a refused Kelvin jump. Geometry and reference inertia stay fixed.

`tests/laminate_thermal_history.rs` checks restrained ramps and holds, moving
thermal/mechanical work, compatible expansion, bimetal residual memory, rotated
anisotropic shear, unsupported spectra and transactional failures. Tensors and
relaxation times are held fixed; varying-temperature relaxation parameters,
mechanical equilibrium and reciprocal thermal evolution are separate owners.

## Ionized-gas DC transport from actual populations

`ionized_gas::IonizedGasConductivity` consumes the exact source-bound
`FrozenIonizedGasStateV1`, including actual finite-rate amounts, volume and
separate electron/heavy temperatures. Its original LTE constructor retains
the accepted `IdealGasEquilibriumV1` solution and delegates to the same
transport calculation. `equilibrium()` is optional; a frozen state does not
acquire an equilibrium receipt. Original charge numbers, particle masses and
actual number densities supply the carriers. There is no Saha projection of a
kinetic endpoint and no invented electron seed.
Each present carrier requires a separately selected positive material property
in reciprocal seconds for its effective momentum-relaxation rate. The complete
original source bundle, claim/key, query and usage receipts remain attached;
temperature and pressure must match the actual endpoint. Common-temperature
sources use typed `T,p`; separate-temperature sources explicitly use typed
`T_h,T_e,p` and cannot reuse a stale common-temperature rate. No neutral-gas
viscosity rule, species-name constant or fabricated electrical-property claim
is substituted for these sources.

Electron drag can also consume the original source-retained Maxwellian
momentum-transfer evaluations. Every present neutral requires exactly one
matching target source; actual densities supply `nu_e=sum(n_target*k_m)`.
The displaced-Maxwellian moment is `k_m=(2/3)<sigma_mt*v*E/(kTe)>`, not the
event-count rate `<sigma*v>`. Ionization/event coefficients refuse as drag.
Original cross sections, evaluation controls, target densities and the exact
frozen state remain attached. The electron source mass must agree with the
collision velocity model. A caller-declared bound on
both `sqrt(me/m_target)` and `sqrt((Th/m_target)/(Te/me))` admits the
stationary-heavy approximation. Cold targets cannot bypass the recoil-mass
requirement by having a small thermal speed.
This is a Maxwellian fluid momentum closure, not an exact solution of the
energy-dependent Lorentz-Boltzmann electrical mobility problem.

With independent background drag, the scalar law is
`sigma=sum(n_s*(z_s*e)^2/(m_s*nu_s))`. The optional consuming
`with_coulomb_drag` builder retains the original `CoulombModelV1` and evaluates
every unordered charged pair at the snapshot's actual two temperatures.
Electron-ion and ion-ion forces use the same symmetric pair coefficient K.
The full momentum system is
`a_s*u_s + sum_t K_st*(u_s-u_t) = n_s*q_s*E`, where
`a_s=n_s*m_s*nu_background_s`. Each original positive background anchor remains;
there is no independent-carrier substitution of an electron-ion collision
frequency and no arbitrary gauge for a pure-Coulomb common velocity.
The existing diagonally scaled Cholesky owner solves at most 128 carrier
velocities once per snapshot. Original force balances and
`sum_s E.J_s = sum_s a_s*|u_s|² + sum_pairs K_st*|u_s-u_t|²`
must close to a relative numerical tolerance of 1e-10. Signed carrier mobility,
current contributions and complete pair receipts remain inspectable.

Both static paths require a prescribed massive neutral background, zero magnetic
field and zero frequency. They are local Maxwellian fluid closures; background
dynamics, nonlocal transport, Spitzer-Harm kinetic electrical conductivity,
Hall effects, electrical breakdown, non-Maxwellian distributions, excited-state
kinetics and radiation are outside their scope. Explicit
caller bounds admit the charged mole fraction, classical phase-space density,
Coulomb coupling and carrier drift relative to thermal speed. These are domain
limits, not certified error bounds or empirical validation. Each carrier's
actual temperature supplies its classical phase-space and drift screens;
the coldest present carrier temperature supplies a conservative Coulomb screen.
The actual field
must satisfy the weak-drift limit before current and heat can be published.
`joule_heating_w_m3` checks that same actual field, including every pair's
relative-drift limit. With Coulomb coupling, it separates nonnegative
background/pair friction heating from the original source's signed O(drift²)
correction to unequal-temperature thermal exchange. Each pair's mass-weighted
friction shares sum to `K*|u_s-u_t|²`; the thermal correction has equal and
opposite pool contributions. A pool's combined thermal power may be negative
and is never clipped or replaced by that carrier's electrical work. Three
scalar thermal coefficients are prepared once, so actual cell heating remains
O(1) regardless of pair count. Electron and heavy powers feed the original
thermal budgets; the unsplit `E.J` and the correction must not be charged a
second time. Zero-drift Coulomb/neutral thermal exchange remains a separate
conservative process in the existing energy owner. State changes require
newly resolved background transport and refreshed original Coulomb moments.

The sealed `conductor::planar::PlanarOhmicResponse` trait exposes the common
static field response. Existing source-material conductor implementations use
their unchanged electrical-property admission; gas states retain their
thermochemical and collision sources instead. Material-specific contact and
refresh APIs still require `PlanarConductorMaterial`.

G0/G1/G3 tests in `tests/ionized_gas.rs` check actual carrier number density and
mass conversion, conductivity and collision-rate scaling, temperature-driven
ionization, frozen carrier abundance, separate temperatures, complete source
coverage, momentum/event distinction, stale-state and stationary-target
refusals. Further cases compare the coupled two-fluid result with an independent
analytic inverse, prove that ion-ion drag changes relative motion, and retain
signed thermal corrections and exact population binding. Cancellation is
polled per species/pair before publication. There
is no unsafe code, evolving state, default source, or new runtime dependency.

### Settled harmonic current from the same carrier sources

`IonizedGasConductivity::prepare_harmonic` retains the exact original state,
background anchors and Coulomb pair moments. For peak phasors with
`exp(+i omega t)`, it solves `(K+i omega M)u=bE`, where `M_ss=n_s m_s`
and `b_s=n_s q_s`. Independent carriers keep an O(N) diagonal path;
coupled transport diagonalizes the mass-normalized symmetric drag once and
reuses its real spectrum for O(N²) frequency responses, capped at 128 carriers.
Collision rates remain reciprocal-second rates without a spurious 2*pi factor.
Original per-carrier momentum and power balances must close within 1e-10.

The response reports three-dimensional peak velocities/current, cycle-mean
background and pair friction heat, the original signed drift-dependent thermal
correction, and separate mean drift kinetic storage U. The conventional
reactive power `Im(sum(E*conj(J)))/2` equals `2*omega*U`; it is not heat.
No zero-drift thermal exchange is attributed to the electrical drive.
Actual instantaneous peak speed, including elliptical polarization and every
pair-relative motion, is checked against the original carrier/source bounds.
Only positive finite angular frequencies are accepted. Exactly DC uses the
static owner; the low-frequency sinusoidal mean approaches one half of the
DC heat at the same peak field.

Populations, temperatures and Markovian rates remain frozen. This supplies no
Hall response, kinetic/nonlocal conductivity, transient initialization,
dielectric background, Maxwell propagation or empirical broadband validation.
Cancellation is polled through species/pair work and around the bounded dense
preparation. `tests/ionized_gas/harmonic.rs` targets independent Drude and
two-fluid complex inverses, DC/high-frequency limits, heat/storage separation,
polarized peak admission, source retention and invalid/cancelled responses.

### Magnetized harmonic carriers and Hall current

`IonizedGasHarmonic::prepare_magnetized` adds an explicitly prescribed uniform
static B to the same frozen carrier populations, signed charges, masses,
background relaxation sources and complete Coulomb pair friction. Under
`exp(+i omega t)` the momentum equation is
`(K+i omega M)u + nq*(B cross u) = nq*E`. A caller-declared ceiling in `(0,1)`
on every carrier's `hbar*|qB/m|/(kT)` screens classical thermal applicability;
it is not a measured material qualification or an error bound.

The parallel response reuses the original relaxation eigensystem. Two circular
polarizations retain the signed cyclotron frequencies and all mutual forces.
Independent carriers keep O(N) preparation; coupled carriers use two bounded,
diagonally scaled N-by-N solves through the existing cancellable complex LU,
with the original 128-carrier limit. Complete original momentum balances must
close within 1e-10. The prepared fixed-frequency/fixed-B owner reuses its
mobilities for O(N + pair-count) arbitrary polarized field queries. Original
circular factors are retained instead of first subtracting large Hall and
transverse coefficients near resonance. The full nonsymmetric conductivity
tensor and circular conductivities remain inspectable.

Actual instantaneous peak carrier and pair-relative drifts enforce the same
source domains, including elliptical polarization and resonant amplification.
Heating uses only the original background/pair friction, including its signed
drift-dependent two-temperature correction. Magnetic force does no real work.
Conventional reactive electrical power includes a signed gyroscopic term:
`Im(E.conj(J))/2 = 2*omega*U + Im(sum(nq*(B cross u).conj(u)))/2`.
This term is not heat; at cyclotron resonance kinetic storage can remain
nonzero even when electrical reactance vanishes. Heat, real magnetic work and
reactive closure are checked separately, with resonant cancellation normalized
by the gross reactive terms. Cancellation publishes no partial preparation or
response, and no original source state is mutated.

`tests/ionized_gas/magnetized.rs` covers the original B=0 response, an independent
Cartesian Drude-Hall inverse, both charge signs' circular resonances, a complete
two-fluid Coulomb inverse and thermal partition, proper-rotation covariance,
Onsager-Casimir field reversal, source-domain refusals and cancellation. This
is a positive-frequency classical local fluid response with prescribed B and
fixed isotropic Markovian collision coefficients. It does not implement
magnetized collision kinetics, quantum or finite-Larmor-radius response,
transient initialization, background momentum evolution, Maxwell propagation,
or empirical broadband validation. Its tensor is deliberately not admitted by
the existing scalar unmagnetized plane-wave or symmetric-Ohmic field owners.

### Original-source refresh after frozen-population heating

`IonizedGasConductivity::at_frozen_state` rebuilds transport at the actual new
heavy/electron temperatures and ideal-gas pressure of the same fixed-volume
parcel. It requires exact original thermodynamics, species amounts and volume.
Every retained electron-neutral momentum cross section is reintegrated at the
new Te using its original weighting, numerical controls and target-speed bound.
Material-backed rates require each original immutable card exactly once, retain
their complete requirements/selected claims and coordinate conventions, and
requery the actual T or T_h/T_e and p. Common-temperature evidence cannot be
widened to unequal temperatures; ambiguous schemas, missing/extra/duplicate
cards and source-domain failures refuse the entire refresh.

The original optional Coulomb model is retained and reevaluated at the actual
new temperatures, with all mutual drag and original transport admission limits.
The new owner carries its frozen-state receipt, without a stale LTE equilibrium
claim. Cancellation or failure leaves the old snapshot intact. This evaluates
an endpoint; it neither evolves thermal energy nor claims source support along
an intervening path. The focused refresh tests exercise source/card retention,
the constant-cross-section sqrt(Te) law, Coulomb reevaluation, fixed-population
and volume binding, common-temperature refusal and transactional cancellation.

`IonizedGasConductivity::at_reacting_state` extends the same source refresh to
changed species amounts after finite-rate chemistry. Original thermodynamics,
volume and present-species support remain exact. Electron drag is reintegrated
at the actual Te and multiplied by the new neutral densities; all carrier
densities and resulting conductivity are rebuilt from the accepted endpoint.
An optional Coulomb model is rebound to these populations using its original
screening, short-range cutoff and physical limits, then reevaluated together
with the original background anchors. Material cards retain their exact
identity, requirements and supported temperature/pressure schema.

No LTE reset, extra ionization heat or stale collision frequency enters this
handoff. Newly appearing or exhausted species require a separately admitted
transport source bank; weak-ionization, source and classical-domain limits
continue to apply. Publication is atomic and the frozen-population entry point
still rejects changing amounts. Focused regressions independently reconstruct
carrier/neutral density scaling, constant-cross-section temperature dependence,
Coulomb population rebinding, Joule heat, reverse population evolution, source
support and cancellation. This is an endpoint material evaluation, not a
reaction integrator or a continuously coupled field/chemistry solve.

## Permanent secondary creep

`creep::NortonCreep` adds sub-yield permanent deformation on the existing
isotropic small-strain elasticity and tensor-Voigt convention. Its normalized
Norton rate uses explicit reference stress [Pa], reference rate [1/s], exponent
at least one, molar activation energy [J/mol], absolute temperature [K] and
physical duration [s]. Arrhenius kinetics are admitted only in the declared
temperature interval. Signed traceless creep strain and its nondecreasing
equivalent path length remain specimen history, separate from plastic strain.

The implicit radial solve is monotone and bounded; a logarithmic safeguarded
Newton solve avoids forming overflowing stress powers. The algorithmic tangent
includes transverse and radial derivatives, including the zero-stress linear
creep limit. Trial evaluations never mutate history. The checked candidate
admission audits its flow residual, consistency with stored strain, and the
independent identity `endpoint work = delta elastic storage + q delta_alpha +
backward-Euler loss`. Numerical loss is never deposited as heat. Elastic
coefficients are fixed during this work account; temperature changes kinetics.

Total, elastic and signed creep Frobenius strain norms must remain within the
same ceiling, at most 0.05; invalid, nonfinite or unrepresentable histories
refuse. This is an authored secondary-creep model until source-bound. It does
not infer creep from melting temperature or yield, and makes no primary/tertiary
creep, anisotropy, finite-strain, damage, rupture or lifetime claim. Focused G0/G1
tests cover stress-controlled creep, strain-held relaxation and timestep
convergence, shear, temperature, reversal, objectivity, tangent consistency,
hydrostatic invariance, work/heat and unchanged inputs on refusal.

`creep::uniaxial::UniaxialNorton` eliminates transverse stress through
`q + E dt rate(q) = |E (epsilon - epsilon_c_old)|`, using the same bounded
logarithmic kinetic solve as the tensor law. The axial tangent is condensed
at zero transverse stress; it is not the tensor C_xxxx entry. Signed axial
history retains Poisson contraction, traceless lateral creep and accumulated
path length. Checked work/heat admission delegates to the original tensor
owner and independently compares axial stress, zero lateral stress and history.
The existing `Uniaxial` interface evaluates immutable physical-step candidates;
consumers commit history once per actual interval. Source-bound laws retain
the original complete material bundle through shared immutable storage.
`tests/creep_uniaxial.rs` supplies analytical constant-stress and held-strain
references, refinement, finite-difference/Schur-complement tangents, reversal,
source-temperature continuation and refusal checks.

`creep::thermal::ZeroExpansionCaloricPoint` couples the source-bound native
history to the existing uniform enthalpy-transport protocol. Its type explicitly
selects zero thermal expansion, fixed reference density/volume and constant
elastic storage; absent expansion data is not a measured zero. One immutable
card, every fixed source condition, paired pinned Cp/density receipts and a
bounded wholly solid caloric interval must agree before execution. Every thermal
proposal starts from the same accepted state and physical clock; final mechanics
are reevaluated at the returned temperature. Physical creep heat alone enters
caloric storage. Separate transport, temperature, heat and combined work gates
precede immutable candidate publication, with bounded iterations and cancellation.
The adiabatic convenience path uses the existing caloric energy-step owner.
This reduction provides no phase change, pressure work, expansion coupling,
spatial-temperature resolution or unconditional thermal stability guarantee.

`creep::source::NortonCreepSource` binds seven explicitly selected claims from
one immutable material card: E, Poisson ratio, reference stress/rate, stress
exponent, molar activation energy and reference temperature. It checks complete
pinned temperature lines, including every original interior knot, and refuses
unsupported spans or changes in those reference coefficients. The Arrhenius law
owns kinetic temperature dependence; elastic storage stays fixed. No yield data
or material-name defaults are used. Each physical step refreshes the original
receipts at its actual absolute temperature and duration, preserving the
original source axis schema. Source-bound
laws implement the same `SmallStrainLaw`, preserving native stresses, tangents
and history. The source checks and manufactured temperature/material swaps in
`tests/creep_source.rs` establish implementation behavior, not measured creep
or service-life qualification.

## Elastic wave and hydrostatic response

`elastic::acoustic::ElasticAcousticState` derives wave and pressure properties
from the original positive-definite Mandel stiffness, inertial density and
material-state identity. Isotropic and orthotropic constructors use the existing
elastic laws and their strain ceilings. Full tensors, including normal/shear
coupling, use the same orthonormal basis and transformations as native mechanics;
no material-name table or assumed isotropy is introduced.

The same immutable state implements `SmallStrainLaw`, so the native plane,
tetrahedral and hexahedral material consumers can use the complete tensor,
including off-axis normal/shear coupling. Its tensor-Voigt tangent is
`D^-1 C_M D`, with `D=diag(1,1,1,sqrt(2),sqrt(2),sqrt(2))`; symmetry is in the
work-conjugate tensor metric. The checked `strain_response` returns Mandel
stress and `epsilon:sigma/2` with the original identity and componentwise
strain ceiling. Infallible trait calls return nonfinite values outside that
ceiling so checked consumers refuse; they cannot silently widen the source
domain. This is one frozen, state-free elastic law, with no plastic history
or thermal eigenstrain. G1/G3 tests compare native P1/Q1 plane stress and
strain against an independent three-dimensional block solve, retained
transverse shear relaxation, tangent finite differences and frame covariance.

Hydrostatic compliance solves `C s = I`, giving compressibility `I:s`, bulk
modulus `1/(I:s)`, pressure strain `-p s` and stored energy `p^2(I:s)/2`.
Directional modes solve the Christoffel problem `Gamma a = rho c^2 a` through
the existing `fs-la` spectral owner and return phase speed, polarization and
wave-normal traction/velocity impedance `rho c`. The original identity survives.
The same eigenpair and full tensor give energy velocity
`v_j = a_i C_ijkl a_k n_l/(rho*c)` without angular differencing or another
eigensolve. It is energy flux divided by total mean wave energy and equals
`grad_k omega` on a smooth separated branch; `n.v=c`, but v can have a
transverse component. At a repeated root it describes the selected polarization,
not a unique differentiable branch. Charge-free piezoelectric waves also retain
the electrical correction `-E_n D/(rho*c)`, with `n.D=0`: differentiating only
a frozen direction-dependent effective stiffness omits this transverse flux.
The fixed-field branch keeps the original mechanical result.
Exact symmetry, positive definiteness, finite representability, solve residuals
and the declared hydrostatic strain ceiling are checked; invalid states refuse
without repairing source data. Degenerate polarizations identify a subspace.

These are frozen, unstressed, small-strain model predictions. Static or
isothermal source moduli do not establish measured adiabatic or high-frequency
ultrasonic response. There is no attenuation, dispersion,
prestress, finite-strain or thermoelastic-relaxation claim. Focused tests cover
analytic isotropic and orthotropic modes, oblique anisotropy, full-tensor frame
covariance, hydrostatic strain/work, density scaling and invalid inputs.
G0/G3 energy-transport regressions compare an independent orthotropic shear
ray, analytic piezoelectric dispersion, full frequency gradients and co-rotation.

## Source-bound axial modulus references

`axial_reference::ResolvedAxialModulusReference` resolves density and an
explicitly named tensile/Young modulus for tension or compression modulus for
compression from one immutable material
card and one complete typed source point. It retains both original usage
receipts and uncertainties. Tension and compression are separate modes;
dynamic or flexural modulus, strength and pointwise stress cannot substitute
for an axial modulus. A source-reported Young modulus does not establish a
tensile-test protocol, directional equivalence or a measured linear range. No Poisson ratio or three-dimensional isotropy is inferred.

This is an authored, zero-origin **linear modulus reference**. A reported
compression modulus fitted over a finite strain interval is not reclassified
as a measured initial tangent. The caller declares a strain ceiling, original
length and area. The existing generalized-Maxwell equilibrium-spring limit
supplies nominal stress and ramp work; the coupon returns `E*A/L`, mass,
force, displacement and **model loading work**.

Loading must increase monotonically in the selected direction. Holds,
unloading, reversal, excess strain and nonrepresentable values refuse before
history changes. When the selected claim declares a crosshead-speed axis, the
caller must identify it and the displacement/time ramp must match. When source
speed is unstated, ramp duration remains an authored model choice. A query's
unconstrained extra axis cannot acquire source-rate authority.

The model does not predict relaxation, creep, cyclic hysteresis, heat generation,
recoverable material energy, lateral contraction or failure. Sampled strain,
rate and history laws cannot be silently frozen into this reference.

Focused checks cover independent axial stiffness/force/mass/work references,
separate loading modes, source/protocol refusals and unchanged state after
rejected steps:

```sh
cargo test -p fs-material axial_reference
```

The `material_axial_reference` example reads a canonical source pack and derives
axis types from its selected source properties. It accepts explicit geometry,
loading mode, source coordinates and up to 64 monotonic displacement ramps.
It validates the entire requested path before emitting a successful response.
The material seed README includes complete PEEK tension and ROHACELL compression
recipes with independent force, mass and model-work expectations.


## Source-bound compression observation references

`compression_reference::ResolvedCompressionReference` binds a selected nominal
engineering-compression stress curve and a strain-independent reference density
from one immutable material card. Original typed queries, claim selections,
uncertainties and exact knot receipts survive. Source interpolation is the
default requirement. The explicit `PiecewiseLinearBetweenObservations` choice
authors a bounded linear reference between discrete observations without
changing their `TabulatedOnly` source policy or issuing an interpolated source
receipt. All original knots must remain supported; no origin or missing segment
is fabricated. Scalars, strengths, moduli, true-stress keys, sampled density and
unowned rate dependence refuse.

A uniform specimen supplies original length and area. It starts at the first
source compression and keeps its original mass. Strictly increasing shortening
produces nominal force, one-sided segment stiffness and piecewise-integrated
loading work **relative to that first state**. Known crosshead speed is enforced;
unknown rate leaves ramp duration authored. Rejected steps preserve history.
No unloading, holds, creep, recovery, cyclic hysteresis, lateral contraction,
failure, dissipated heat or recoverable-energy claim is made. The model does
not identify a source test's geometry with the caller's geometry.

`material_compression_reference` loads canonical packs and validates the entire
requested path before printing successful output. Its explicit
`--interpolate-observations` flag exposes the model-authoring choice.
`tests/compression_reference.rs` checks the source and numerical boundaries;
the catalog's compression-reference tests exercise four actual foam grades
and reinforced PEEK through the production compiler and this owner.

## Accepted-state fatigue observation for cyclic metal laws

`fatigue::cyclic::SourceBoundFatigueLaw` admits the original source-bound
J2/Prager, Voce and Chaboche laws, including one through eight independently
resolved recovery branches and borrowed native laws. Static dispatch retains
the actual associated state and complete immutable material receipts. The
common elastic/initial-yield backbone serves calibration binding and accepted
stress recovery; it never replaces the native nonlinear constitutive update.

Accepted physical stress is the original elastic map of total mechanical strain
minus accepted plastic strain. Backstresses translate the yield surface and
are not subtracted from that physical stress. Consumers that already retain
native stress use that tensor directly. Observation performs no return map,
state update, source query or per-sample source clone. Native history validation
retains finite/strain/plastic bounds and every Chaboche branch's deviatoric,
reachable and saturation bounds. Equilibrium and yield convergence remain
owned by the native solve, without an invented observation tolerance.

The interface is sealed to these fixed-condition laws. Thermal eigenstrain,
creep, viscous histories and changing-temperature calibrations require their
own admitted mechanical-strain decomposition and are not silently included.
Focused regressions in `plastic_chaboche_source` and
`plastic_chaboche_branches` compare accepted physical stress to native updates,
exercise nonproportional/reversed strain and borrowed sources, and reject
unreachable histories including an invalid last branch in an eight-branch law.
Synthetic numerical checks do not establish measured fatigue life.

## Source-resolved strain-life calibration

`fatigue::strain_life::source::resolve_cmb` resolves the original elastic modulus
and four fitted fatigue coefficients from one immutable material card. Every
property retains its selected claim, full material answer, uncertainty and query
receipt. The four fatigue coefficients must share one declared fit provenance,
observation population, fully reversed stress baseline and finite life interval.
Missing cyclic coefficients are never inferred from monotonic material properties.

The source power convention is explicit: a coefficient multiplying `N^p` becomes
`A * 2^(-p)` in the native `(2N)^p` law. The original coefficient is unchanged in
its receipt. Source life bounds are converted to full cycles and become the
curve's strict model domain. Coefficients may be scalar or vary with the declared
fixed temperature; a life-dependent coefficient is not a CMB fit and refuses.
Source interpolation remains owned by the material database.

The retained calibration can check the actual structural material card, named
state, external query coordinates and quantity conventions, and original modulus
claim. Calibration-only life and stress-ratio coordinates are distinguished from
operating conditions. This source binding establishes traceability and supported
numerical use, not measured component life or a validated population transfer.
Eight focused tests in `tests/fatigue_strain_source.rs` exercise analytical count
conversions, original receipts, source selection, context and fit-domain refusals.

## Paired stress and mechanical-strain fatigue corrections

`fatigue::strain_life::PairedStrainLifeMiner` evaluates named Morrow and
Smith-Watson-Topper policies from simultaneous signed normal stress and total
mechanical strain on one fixed plane. Morrow uses the actual signed mean stress;
SWT uses the actual maximum tensile stress times strain amplitude. The original
finite cycle interval remains binding, with each correction's amplitude support
recomputed at those bounds. The same safeguarded log-life inverse handles all
three strain-life policies. A scalar strain-only counter refuses corrected
policies because it has no stress observation.

The existing rainflow owner is generalized with a per-turning-point context:
scalar paths have zero-sized context and paired paths retain the observed
stress. There is one counting algorithm and one compensated Miner sum. Samples,
blocks and finalization retain their transactional cancellation/refusal behavior.
Stress and strain increments must have matching directions; stress changes at
unchanged strain and nested reductions that would hide an unmatched stress
extreme refuse. These are observed coincident-extrema checks, not proof of
unsampled peaks, stabilized material behavior or general multiaxial validity.

Morrow's compressive-mean credit is part of that explicitly selected empirical
model. A nonpositive Morrow elastic fatigue coefficient refuses; it is not a
static-strength failure criterion. Nonzero-strain SWT cycles without a tensile
maximum also refuse rather than imply infinite fatigue life. Applying these
constant-cycle equations to evolving cycles remains an empirical approximation;
no hardening, relaxation, sequence, heat or constitutive degradation law is added.
Seven focused analytical and history tests cover both equations, cycle/strain
conventions, real stress pairing, nested histories, chunking and rollback.

## Closed-cycle tensor SWT critical-plane fatigue

`fatigue::critical_plane::SwtCycle` evaluates one resolved representative
stress/total-mechanical-strain cycle on an explicitly declared set of material
planes. The tensor convention is `[xx, yy, zz, xy, yz, zx]` with tensor shear
in one fixed frame. Thermal eigenstrain must already have been removed by its
physical owner. Every observation contributes independently to normal-strain
and normal-stress extrema on every plane; stress peaks need not coincide with
strain turning points.

The selection convention is maximum normal strain range, followed by maximum
normal stress for exactly equal ranges, then stable input order. The selected
plane alone enters the existing source-bound CMB/SWT life inverse. This is not
a search for the maximum SWT product, a shear-failure model or a multiaxial
rainflow rule. The existing finite life interval and compression-only refusal
remain binding. Zero strain variation produces no fatigue contribution, not
a validated infinite-life statement.

Finalization requires the first and last full tensors to agree within explicit
closure tolerances and every nonconstant projected strain trace to contain one
cyclic excursion. Open paths and nested/repeated strain subcycles refuse.
Closing a sampled stress/strain loop does not prove constitutive stabilization,
unsampled extrema or applicability to a particular material population.
The caller must justify the representative cycle and fixed-temperature CMB
calibration. No evolving hardening is replaced by a repeated native state.

Storage is proportional to the declared plane count, independently of the
number of observations. Sample, plane and projection-work budgets are explicit;
samples, blocks and finalization publish atomically under cancellation/refusal.
The deterministic hemisphere helper supplies a finite angular search whose
resolution can be refined. Results report the selected sampled plane and
observed closure defects; no continuous-sphere optimality or certified search
error is claimed. Scalar rainflow and the existing paired counter are unchanged.

## Purpose and layer

### Shared stress-life and total mechanical strain-life fatigue

`fatigue` owns one bounded streaming rainflow reducer for both S-N/Miner
stress-life and `fatigue::strain_life` Coffin-Manson-Basquin/Miner strain-life.
Signed physical normal stress and signed total mechanical normal strain are
different inputs. Tensor shear uses `[xx, yy, zz, xy, yz, zx]`; unsigned von
Mises stress and accumulated equivalent plastic strain cannot substitute for a
reversing history. Fixed-plane projection is not a multiaxial failure criterion.
Counters retain unresolved reversals across blocks, compensate the damage sum,
check explicit sample/residue limits, and publish samples, blocks and final
endpoint half cycles transactionally. Finalization is idempotent. The original
source or calibration remains shared through the immutable curve.

`preview_finish` evaluates the accepted prefix as if observation ended now,
including the pending endpoint and all residual half cycles. Its returned
summary describes that hypothetical finalized prefix (`finished = true`);
the live counter's summary, residue, sample budget and continuation are unchanged.
The cancellable form uses the same cycle admission and life law as `finish`.
Unsupported provisional amplitudes still return the original domain error.
The query neither commits provisional cycles nor estimates unobserved future
loading. Stress-life, scalar and paired strain-life, and variable-amplitude
critical-plane counters share this meaning. Critical-plane previews compare
the complete provisional damage of every declared plane without closing any
plane's live history.

The total-strain law is `epsilon_a = sigma_f'/E*(2*N_f)^b + epsilon_f'*(2*N_f)^c`:
amplitude is half the signed range, `N_f` is full cycles and `2*N_f` is reversals.
Both empirical contributions use the supplied calibration and finite life domain;
a safeguarded deterministic log-life inverse performs bounded work. Zero
amplitude contributes no damage; positive amplitudes outside the declared
domain refuse rather than inventing endurance or extrapolation. The explicitly
selected `Uncorrected` strain-life policy neglects mean stress; strain mean does
not supply the stress extrema needed by Morrow or Smith-Watson-Topper models.
The formulation follows the combined relation documented in COMSOL's
[Strain Life Models](https://doc.comsol.com/6.4/doc/com.comsol.help.fatigue/fatigue_ug_sme.4.21.html).

`fatigue::source::resolve_log_log_sn` admits an original material-card
`fatigue-stress-amplitude` or `fatigue-stress-range` curve against explicitly
named cycles-to-failure or reversals-to-failure. The caller deliberately authors
log-log interpolation between the original observations. Every original knot,
claim, uncertainty, validity condition, typed query and selection receipt stays
attached to `SnCurve::source`; a tabulated-only source remains tabulated-only.
No point is dropped to repair monotonicity, avoid a source conflict or escape a
validity boundary. The original resolver requires R = -1;
`resolve_log_log_sn_at_stress_ratio` retains one exact source-declared fixed
ratio in -1 <= R < 1, including pulsating and tension-tension calibrations.
Both require the actual absolute temperature. Strength, censored runouts and an
unsupported endurance limit cannot replace finite-life observations.
Applicability labels remain caller declarations; material/process/environment
equivalence is not inferred.

`SnCurve::new_at_stress_ratio` likewise retains an explicitly authored fixed
calibration ratio. `MeanStressCorrection::MatchCalibration` checks every
counted full cycle and endpoint half cycle against that ratio, allowing only
floating-point roundoff. The cycle maximum must be tensile. Amplitude remains
half the signed range; the measured mean-stress effect is already in the curve.
Non-R=-1 curves cannot silently use uncorrected or Goodman policies. A mismatch
refuses the requested block or finalization transactionally, preserving the
accepted prefix for continuation. Normalizing amplitude and mean before the
ratio check avoids intermediate stress overflow; each cycle adds constant work
and no allocation. No load-ratio interpolation or mean-stress law is inferred.

`SnCalibrationSource::matches_material_state` binds the original card, named
material state and every fixed typed physical query coordinate to a native
material bundle. Only the calibration life axis and calibration R axis are
excluded from physical-state matching. Custom temperature axes and additional
fixed environmental coordinates remain exact. A mean-stress correction is a
separately declared empirical model; it does not confer new source support.

Focused G0/G1/G3 regressions check analytic strain-life inversion and single-term
limits, cycle/reversal and range/amplitude conversion, hand-countable damage,
shared rainflow behavior, chunk/restart equivalence, source conflicts and
unsupported-domain rollback. These are implementation checks using synthetic
fixtures. Miner exhaustion is an empirical estimate, not a rupture certificate,
stiffness degradation, heat source, crack-growth law, validated component life
or a claim of transferable calibration across temperature or load populations.

Layer **L3** (FLUX support). Runtime deps: `std`, fs-ad (dual-number
energy derivatives), fs-blake3 (typed content identities), fs-evidence
(model cards, Evidence, and V&V artifacts), fs-matdb (immutable material and
constitutive-model cards), fs-qty, fs-math.
Consumers: fs-solid elasticity (tfz.13), fiber-section beams, lattice
homogenization, the P2 milestone.

## Public types and semantics

- `state_point::VisibleConductorStatePoint::try_from_material` adopts all
  eighteen canonical dimensionless n/k knots from a complete resolved bundle,
  retaining additional selected properties and the full source identity. This
  lets spectral wave consumers retain optical constants and permeability from
  one snapshot. Adoption and pinned refresh recheck positive real index,
  nonnegative extinction and fixed external source conditions; missing knots,
  substituted quantities and unowned dynamic coordinates refuse atomically.

- `state_point::SpectralComplexIndexStatePoint` adopts 2–16 explicitly supplied
  vacuum-wavelength/full-property-key pairs from that same complete material
  bundle. It retains density, auxiliary properties, exact quantity conventions
  and all original receipts. Real n requires a dimension-only scalar; k also
  admits the existing static optical-extinction quantity. Tensor, hardness and
  unowned dynamic/spectral contexts cannot fill those roles. Every source key
  occupies one declared knot role, with finite positive increasing wavelengths,
  positive n and nonnegative k. No visible table is implicitly extended to IR.
  The projection identity binds the original bundle plus wavelength/key roles.
  `at_wavelength` returns exact source knots or nominal piecewise-linear n/k
  within the closed support, without allocations or source queries. It refuses
  extrapolation; it does not certify unresolved spectral features or propagate
  marginal source uncertainty. `at_query_point` refreshes every retained optical
  and auxiliary property under its original claim pin before publishing the
  candidate. Synthetic G0/G1/G3 regressions cover IR interpolation, exact keys,
  full-bundle thermal refresh, competing sources, passivity and atomic refusal.
  The wavelength attribution is an explicit caller declaration; this adapter
  supplies no additional measured-data authority or thermal-emissivity model.

- `state_point::VisibleOpticalStatePoint::at_query_point` and the conductor /
  dielectric variants continue the original visible optical law at a new
  complete physical condition. They reuse the shared pinned-claim refresh,
  retaining the immutable card, every exact property key, axis convention,
  consumer bound, optical family and spectral-versus-RGB absorption schema.
  The original decoders also reapply dielectric passivity before publication.
  Competing claims cannot replace an unsupported original source. This is an
  atomic endpoint observation, not phase or surface evolution, nor a guarantee
  of source support along the intervening path. Focused G0 regressions cover
  overlapping replacement claims, unsupported endpoints, passive RGB bounds
  and continued spectral absorption on a card that also declares RGB data.

- `wood_moisture::ResolvedWoodMoisture` binds the original 12% dry-basis
  density and independent green-to-ovendry radial, tangential and volumetric
  shrinkages, plus the source-declared fiber-saturation approximation. Exact
  keys, typed source queries, card identity and all five receipts survive.
  Equation (4-9) gives dimension ratios relative to 12% moisture; independent
  volume data and water mass then give density while conserving ovendry mass.
  The admitted moisture interval intersects the coefficient domains with the
  below-fiber-saturation range. No longitudinal dimension or deformation
  tensor is inferred from independent population averages. Moisture is
  prescribed; this does not solve drying, diffusion, hysteresis or stress.

- `section_temperature::thermoelastic::ConcentricThermoelasticLaw` derives
  layer force, Helmholtz storage, entropy and thermal curvature from the same
  source-bound axial potential. Each layer retains its own temperature, card,
  original annular geometry and integrated thermal-expansion path. Scalar
  reference density and expansion, plus scalar or globally affine `E(T)`, are
  admitted; every original modulus knot is checked. Continuous material support
  and all original property selections survive reevaluation, including auxiliary
  electrical/thermal properties. For elastic strain `e`, the per-reference-length
  potential is `psi=E*A*e^2/2`, entropy is
  `s=A*(E*alpha*e-E'*e^2/2)` and curvature is
  `psi_TT=A*(E*alpha^2-2*E'*alpha*e)`. Consumers add `psi+T*s` to the original
  stress-free caloric energy and must admit positive combined heat capacity,
  `C_caloric-T*psi_TT`, on their actual displacement and temperature domain.
  This shared law supplies constitutive derivatives; continuing motion, heat
  transfer and atomic publication belong to the native coupling owners.
  It does not infer zero expansion, average materials, evolve mass, model slip,
  or certify a radial/finite-strain/phase-changing constitutive law. The focused
  G0/G1 tests cover unequal layer temperatures, opposing stresses, source pins,
  negative expansion and independent free-energy derivative identities. Those
  tests were added but not executed natively in the connector-only session;
  syntax parsing and independent derivative arithmetic were checked.
- `section::winding::CloseWoundSection` retains a sourced circular core and a
  distinct round-wire winding. Exact helical length gives winding line mass;
  the close-coiled spring law gives its axial rigidity with core and winding
  constrained to share end displacement while slipping along the span.
  `WoundCoreSection` retains the core's sourced Kelvin–Voigt bending viscosity,
  independent per-layer Prony histories, or original `ConcentricModalLossSection`
  loss-angle/thermoelastic law, with common-state and intersected core/wrap
  source-band admission. The consumer projects core bending arms and evaluates
  modal loss using total winding-plus-core dynamic mass and actual frequency;
  the wrapping material never inherits a bending spectrum or becomes an extra
  bonded bending layer. Circular thermoelastic core loss explicitly retains the
  insulated homogeneous-core thermal boundary and original core radius. It
  supplies no core/wrap heat transfer or thermal bending law for the winding.
  Geometric closest-turn separation rejects interpenetration. Thin wire,
  small helix angle and long wavelengths are explicit approximation domains.
  The consumer propagates these coefficients through the existing string/body
  mechanics, gas drag and physical pressure. G1 compares mass against resolved
  helix quadrature and stiffness against spring strain energy; G3 changes only
  winding material data. This sliding/core-dominated approximation supplies no
  winding friction/loss, torsion, thermal mismatch or resolved contact, and
  does not claim acoustic equivalence to a measured or resolved specimen.

- `section::wound::HelicalWrapping` is the mass-only interim for wraps outside
  the close-coiled admission domain (thick wire/coil ratios, several wraps, as
  on piano bass strings). Each wrap lies on the diameter beneath it; its line
  mass is the wire density times area times the exact helix length per axial
  length. Neighbouring turns must clear one wire diameter at their exact
  closest approach (any helix angle); `close_wound_pitch_m` solves the
  tightest admissible pitch. Wraps resolve at the core's query coordinates.
  They carry no tension and add no EA, EI, loss, friction or wrap prestress:
  the core section owns those. G0 tests check the helix integral, the second
  wrap's radius, brute-force turn clearance and every refusal.
- `section::transverse::ConcentricTransverseViscosity` retains each original
  shear section, material state and explicit shear-viscosity selector. Its
  corrected rigidity is `C_s = kappa * sum(eta_s A)` in each transverse
  direction, using the same declared uniform-shear correction as elasticity.
  Fixed-state constant coefficients, original source bands and layer order
  survive; absent, negative, wrong-quantity and frequency-local transient
  substitutes are refused. This is a small-strain section reduction, without
  resolved cross-sectional warping or interlayer slip.
- `section::strength::AxialTensileCapacity` binds a caller-named tensile-limit
  claim [Pa] (yield, proof or ultimate strength; never chosen, factored or
  converted to an allowable here) from every bonded layer. A common axial
  strain gives layer stress `E_i N / EA`; the capacity is
  `min_i(limit_i EA / E_i)` and that layer governs. A layer without the claim
  refuses as missing data, never as unbreakable. The capacity retains the
  original section and the complete strength answers, including the explicit
  caller-owned diameter-resolved path. Both constructors require strength
  sources static at held external conditions; moisture/pressure/composition
  dependence is retained, while unowned strain, stress, rate and frequency
  applicability refuses. The external-source path additionally admits a
  positive, explicitly length-typed `wire_diameter` held as a fixed drawing
  coordinate; it cannot freeze a live response coordinate into a constant.
  `assess_axial_bending` reports each annulus's normal-stress extrema
  `E_i N/EA +/- E_i r_outer_i |M|/EI`, comparing only positive maximum stress
  with that layer's own limit. Signed compression and biaxial bending are
  admitted. `assess_thermoelastic_axial_bending` requires the exact original
  section/material states and reuses its integrated thermal eigenstrain and
  current elastic layer stresses, including self-equilibrated tension at zero
  resultant force. Total and elastic strain bounds apply to the bending
  extrema as well. Nonfinite/unrepresentable arithmetic refuses. G0/G1/G3
  analytical regressions cover combined loads, layer-specific governing
  limits, thermal self-stress, rotation, subdivision and source/state refusals.
  These are small-strain, perfectly bonded elastic stress assessments;
  residual drawing/winding stress, shear/torsional failure, yielding
  redistribution, crack evolution, buckling and fatigue are absent.
- `section::strength::KelvinVoigtTensileCapacity` prepares both original axial
  and bending viscosity owners against the exact tensile-capacity section,
  geometry, layer order and material sources. Each layer retains its own E,
  axial eta and bending eta: mean stress is `E_i epsilon + eta_ai epsilon_dot`,
  and bending amplitude is `r_i |E_i kappa + eta_bi kappa_dot|`. Signed elastic
  and viscous vectors combine before their magnitude; section-average eta/E
  cannot assign another layer's loss stress. Original viscosity selectors,
  receipts, integrated rigidities and intersected source bands survive; actual
  history/frequency admission remains with the dynamics owner. Cached layer
  coefficients and single-layer assessment allocate nothing during evaluation.
  `AxialTensileCapacity::assess_layer_normal_stress` also admits actual signed
  mean stress and nonnegative bending amplitude from a consumer that owns the
  complete layer history, without an inferred equivalent stiffness. Both paths
  compare only positive total normal stress against the unchanged sourced
  limit. These observations supply no rate-dependent strength, thermal
  eigenstrain, hereditary history, plasticity or fracture evolution. Manufactured
  regressions cover distinct layer viscosities, signed cancellation, rotation,
  subdivision, lossless limits, source/state mismatch and numerical refusals.
- `section::loss_angle::ConcentricLossAngleSection` binds a named dimensionless
  loss tangent `Im(E)/Re(E)` from every bonded layer: the complex bending
  stiffness is `EI + iH`, `H = sum(E_i I_i tan(delta_i))`, with no averaged
  layer loss. A quality factor or decrement is never converted into it. A
  modal consumer uses the equal-energy-per-cycle rate `H k^4/(mu omega)` at
  each actual mode, inside the intersected source band. Missing, negative,
  wrong-quantity, state-varying and band-less claims refuse; explicit zero is
  lossless. Thermoelastic relaxation from thermal data is in `section::ted`.
- `section::ted::ConcentricThermoelasticBending` computes Zener thermoelastic
  bending loss of a homogeneous solid circular member from its one material
  state: `Delta = E alpha^2 T/(rho c_p)` with T from the state's single typed
  absolute-temperature axis, and the complete insulated-cylinder thermal
  diffusion response in `visco::circular_ted`. An allocation-free Bessel
  continued fraction and scaled asymptotics retain the infinite-spectrum
  loss, its `7u/24` low-frequency and `4/u` high-frequency limits, and signed
  frequency symmetry, where `u = omega rho c_p R^2/k`. Finite extreme scales
  do not turn physical loss into zero through squared-frequency overflow;
  invalid frequency or thermal time returns NaN. The scalar path does not
  initialize thermal roots. The explicit 32-mode table remains available:
  `tau_n = rho c_p R^2/(k xi_n^2)`, `J1'(xi_n) = 0`, and weights
  `8/(xi_n^2 (xi_n^2 - 1))`. Its omitted static weight below 2e-6 does not
  bound relative loss in the high-frequency tail. Independent high-precision
  Bessel fixtures, spectral moments and extreme-scale checks cover the scalar
  response; certified J1 and projection integrals check the explicit table.
  This remains first-order weak thermoelastic coupling for a homogeneous
  isotropic section with an insulated boundary. Layered, hollow,
  untyped-temperature and incomplete
  inputs refuse. It is a lower bound on bending loss: intrinsic anelastic
  excess is the consumer's measured addend or a reported unknown.

- `section::relaxation::ConcentricPronySection`,
  `section::relaxation::axial::ConcentricPronyAxialSection`,
  `section::torsion::ConcentricPronyTorsionalSection` and
  `section::relaxation::transverse::ConcentricPronyTransverseSection` compose
  independent generalized-Maxwell bending, axial, torsion and transverse shear
  spectra on the existing bonded concentric geometry. All four use one
  `LayerRelaxationProperties` coefficient resolver;
  the bending-specific spelling names the same source-role type. Each layer
  selects its equilibrium modulus and zero or more modulus/time pairs.
  An explicit empty spectrum is elastic; an omitted layer or absent claim
  refuses. Equilibrium must match the actual channel's Young/shear modulus.
  Each arm retains `delta_E * I` for bending, `delta_E * A` [N] for axial
  extension or `delta_G * J`, `J=2I`, for circular torsion. An axial member
  supplies its reference length to obtain spring stiffness `delta_E*A/L`;
  each layer's original Young-modulus spectrum and area remain distinct.
  Transverse arms retain `kappa*A*delta_G` [N] in their
  distinct corrected-area record, never a polar-moment rigidity. The original
  shear roles, correction and equilibrium stiffness remain authoritative.
  A receiving member supplies its own spatial strain basis and memory states;
  no G(t) is inferred from E(t). Time constants and layer identities are never averaged.
  The original torsional Kelvin viscosity remains a separate parallel dashpot;
  it is not counted again in the Prony spectrum. Geometry, mass, equilibrium
  stiffness and all source receipts survive without authority promotion.
  Density and every coefficient must be constant at the resolved temperature.
  Each relaxation time (or equilibrium modulus for an empty spectrum) supplies
  frequency support, intersected with every other active coefficient domain.
  At most 128 arms per selected section channel are admitted. Missing, ambiguous,
  wrong-quantity, incompatible-band, negative-increment, nonpositive-time and
  unrepresentable-stiffness inputs refuse atomically. G1/G3 tests check independent
  complex axial/bending/torsional stiffness, unchanged parallel viscous loss and
  homogeneous subdivision; G0 covers empty layers and invalid source roles.
  This is frozen-state small-strain bonded section relaxation, not helical
  winding, interlayer slip, heating or experimental validation. Axial creep
  and relaxation histories belong to the receiving dynamics owner; the section
  binding supplies the original material law without introducing a time stepper.
- `state_point` is the fail-closed bridge from immutable `fs-matdb`
  `MaterialCard`s to executable scalar property bundles. Consumers declare
  property names, dimensions, numerical domains, the complete physical query
  point, and an explicit claim-selection rule. Resolution is atomic, retains
  each `PropertyUsageReceipt`, canonicalizes requirement order, and binds the
  card, state point, selected values, consumer domains, and receipts into one
  identity. It never infers properties from chemistry names or substitutes a
  representative material.
  `resolve_visible_tabulated_dielectric_state_point` resolves eighteen explicit
  phase-index and intensity-absorption values at 380, 430, ..., 780 vacuum nm.
  `optical_dielectric_eta_*` is dimensionless with values in `[1, 3]`;
  `optical_absorption_per_m_*` has inverse-length dimensions and is nonnegative.
  These are separate from conductor complex-index `k` and from Cauchy/RGB
  properties. All knots must be supported at one complete material state;
  missing knots, wrong units, unsupported states and invalid values refuse.
  The generic visible optical resolver recognizes this family and refuses
  simultaneous complete alternative families. Explicit family selection remains
  available. The renderer interpolates index and absorption between knots without
  fitting or RGB conversion. Source receipts retain their original authority;
  the real-index boundary with interior absorption is a weakly absorbing,
  unpolarized geometric-optics model, not a complex-interface or thermal law.
  G0/G1 material tests cover state interpolation, all retained source receipts,
  sample completeness, dimensions, phase-index domain and family ambiguity.
  `ScalarPropertyRequirement::try_with_quantity` carries an exact `QuantitySpec`
  through both bulk and interface queries, including caller-pinned selection.
  The answer retains the same quantity kind and value form. Dimension-only
  requirements remain explicit and cannot consume semantic claims. Bundles
  containing semantic requirements use identity domain v2 and bind every
  requirement's canonical quantity token; all-dimensional bundles retain the
  exact v1 preimage. Named convenience resolvers still declare their existing
  dimensional property requirements; axis descriptors travel exactly as
  declared by the caller, without inference from names.
  `ScalarPropertyRequirement::try_with_key` also retains the existing matdb
  hardness test context. Both material and interface resolution forward it
  under all three selection policies; a claim pin cannot fill in an omitted
  apparatus/protocol/specimen requirement. The resolved requirement exposes
  that context, and selected claim/usage-receipt identities bind it into the
  existing bundle identity. The focused synthetic G0 regression proves this
  adapter behavior, not measured hardness accuracy or standards compliance.
  Distinct complete keys may share a name. Duplicate keys refuse, and same-name
  pins must identify distinct claims matching the required contexts. Bundles
  order same-name answers by selected claim identity; reversed requirement/pin
  order cannot alter the result. `property(name)` returns `None` for an ambiguous
  name; `property_by_key` requires the complete quantity/test/tensor context.
  The elastic resolver uses that exact lookup for all 36 matrix coefficients.
  `ResolvedScalarProperty::is_constant_at_fixed_temperature` retains the exact
  selected claim's independent axes and admits either validity-wide scalar constants
  or T-only curves evaluated at positive absolute temperature. Frequency curves,
  multiaxis surfaces, other state axes and temperature-difference coordinates refuse. Legacy T
  remains kelvin; typed coordinates must be dimensional temperature or static
  absolute temperature. The original interpolation decision, evidence and
  identities remain unchanged. This classifies a declared source model;
  consumers must still check every operating frequency against its validity.
  G0 tests cover axis semantics, identical values from differently pinned
  curves, and unchanged source receipts. It does not establish measured accuracy
  or admit a changing-temperature constitutive history.
  `ResolvedScalarProperty::angular_frequency_band_rad_s` reads one unambiguous
  frequency coordinate from the exact query and converts its selected-claim
  applicability to rad/s. Static cyclic frequency is multiplied by 2*pi;
  static angular frequency and legacy `omega` retain their values. Semantic
  descriptors take precedence over axis names. Multiple frequency coordinates,
  incompatible `omega` semantics, negative or overflowing bounds refuse;
  inverse-time dimensions alone do not identify an unrelated axis as frequency.
  Absent source bounds remain unconstrained. This nominal unit conversion does
  not modify receipts, extend source support or freeze a sampled frequency curve.
  G0 tests cover conversions, source preservation, ambiguity and invalid bounds.
- `ResolvedMaterialStatePoint::query` reconstructs its complete physical query;
  `at_query_point` refreshes every retained property from the exact original
  card and selected claims. Temperature, pressure, moisture, shear rate and
  other declared coordinate values can change together. All coordinate names
  and quantity conventions, complete property keys (including same-name tensor
  and test contexts), and consumer admissibility remain fixed. A dependent
  property leaving source support refuses the whole candidate instead of
  switching sources or leaving a stale value. `ResolvedLiquidState` uses this
  same refresh and recalculates all derived transport coefficients with its
  original dynamic/kinematic viscosity basis; section and liquid temperature
  laws delegate to it. Focused G0 regressions cover simultaneous multiaxis
  updates, complete-key preservation, typed-coordinate refusals, source pins
  and atomic failure. These are endpoint evaluations, not history integration,
  support along an intervening path, or measured multiaxis validation.
  `ResolvedInterfaceStatePoint::{query, at_query_point}` uses the same refresh
  implementation for ordered surface systems. Every property, complete key,
  selected claim, coordinate convention and consumer domain is retained; any
  unsupported property refuses the whole candidate. A different surface card,
  preparation/history identity or axis schema cannot be substituted by refresh.
  `ResolvedElectricalInterface` and `ResolvedLiquidWettingState` can adopt an
  already resolved multiphysics interface bundle without querying again. Their
  electrical-resistance and equilibrium-angle admission remains active even
  when the generic bundle requested only finite values. Both adapters and
  `ResolvedLinearInterface` rebuild their derived response through this shared
  endpoint refresh, retaining all other interface properties and original source
  receipts. G0/G3 `tests/interface_refresh.rs` exercises simultaneous temperature
  and pressure changes across electrical, mechanical and wetting projections,
  physical coefficient reductions, deterministic refresh, source pinning and
  atomic refusal. This supplies no wear/hysteresis integration, contact geometry
  update, path-support proof or empirical validation; the caller still owns
  evolution and the accepted-state transaction.
  The shared `Surface2D` payload evaluates both source coordinates together
  through the original matdb evaluator. `ResolvedScalarProperty::source_axes`
  retains every declared dependence even at an exact vertex; a surface cannot
  masquerade as a constant or temperature-only coefficient. Sparse source holes
  remain refusals during complete-bundle refresh. Fixed-context temperature
  paths use the canonical `ClaimSet::query_line_typed` seam to retain every
  crossed source-cell breakpoint and prove continuous support. Caloric,
  thermal-expansion and Seebeck consumers integrate those same source slices
  with their existing physics owners. Endpoint support alone never authorizes
  integration across a source hole or nonzero discrete-only span.
- `conductor::resolve_uniform_conductor` binds exactly one sourced
  `electrical_resistivity` claim in ohm-metre dimensions to caller-owned
  positive uniform length and cross-sectional area, retaining the resolved
  material bundle and receipt with `R = rho_e L / A`. It refuses invalid
  geometry, source/unit/domain/selection failures, and unrepresentable derived
  resistance. `joule_power_w` evaluates finite signed-current `I^2 R` only;
  it is not a circuit, contact-resistance, current-distribution, temperature,
  or energy-balance model. The adapter is synchronous, fixed-cost arithmetic
  and has no cancellation claim.
  `OhmicDrive` distinguishes ideal current and voltage control;
  `current_a_for_drive` resolves `I` or `V/R`, refusing nonfinite inputs and
  unrepresentable nonzero currents. The same `joule_power_w` law then supplies
  nonnegative source dissipation. No dynamic circuit or source impedance is inferred.
  `ElectricalPropertyBasis` also admits explicit scalar DC conductivity through
  `resolve_conductor_state_with_basis` and `resolve_uniform_conductor_with_basis`.
  The source is evaluated at the query point before taking its reciprocal;
  reciprocal knot interpolation is not substituted. The original property,
  uncertainty and receipts remain intact under `source_property()`; derived
  resistivity is never inserted as a measured claim. Existing resistivity query
  entrypoints keep their requirement unchanged. `try_from_material` can adopt a
  complete multiphysics bundle with either one source basis; two bases require
  an explicit independent-variable choice and are not averaged. Consequently,
  concentric sections can mix conductivity- and resistivity-sourced layers in
  their existing parallel-current and layer-specific Joule-heat calculation.
  This scalar reciprocal law does not invert tensors or admit AC response.
  `resolve_conductor_state_with_property` and
  `resolve_uniform_conductor_with_property` bind a caller-selected exact
  `PropertyKey` to that same explicit resistivity/conductivity basis. They do
  not normalize names or infer aliases. The key must retain the basis's exact
  dimensional quantity; the complete original bundle, conditions and receipts
  remain authoritative. `source_key()` survives bulk-to-uniform lowering and
  material refresh. Planar source-name iterators borrow the actual bound names,
  and coupled refresh guards retain the selected key even when two properties
  in the same bundle have equal values. The canonical entrypoints retain their
  original requirements and refusal behavior.
- `field_storage::{ResolvedCapacitor, ResolvedInductor}` adopt complete material
  state bundles and derive `C = epsilon A / d` and `L = mu N^2 A / l` from
  explicit uniform geometry. Absolute SI electric permittivity [F/m] and
  magnetic permeability [H/m] are required; relative coefficients cannot be
  silently substituted. Original source claims and state receipts are retained.
  Charge and flux-linkage energies use `q^2/(2C)` and `lambda^2/(2L)`.
  Only linear isotropic, nondispersive, lossless uniform-field reductions are
  admitted. Curves and surfaces may be resolved at fixed temperature and typed
  pressure, composition, moisture or concentration. The shared fixed-external
  source policy checks every value axis and applicability bound; dimensional
  equality does not turn stress into pressure or strain into composition.
  Frequency-, field-, strain- or stress-dependent/restricted sources require
  another storage model. Changing a resolved coefficient is a new snapshot;
  live evolution must account its parameter work and any field entropy.
  Invalid geometry, ambiguous sources and unrepresentable coefficients/energy
  refuse. These reductions introduce no circuit solver, hysteresis or loss law.
- `thermoelectric::ResolvedSeebeckSpan` resolves one source-owned scalar, linear-inside curve
  or fixed-context two-coordinate surface slice of
  `absolute_seebeck_coefficient` [V/K] at both ends of a fixed-context temperature
  span and at its actual bulk temperature, preserving original receipts and the
  same selected claim. Joule-heated bulk states cannot escape source support
  merely because the junctions remain cooler. Its
  zero-magnetic-field convention is `J = sigma(-grad(V)-S grad(T))`,
  `q = -k grad(T)+S T J`. Exact piecewise-linear integration supplies `S_left`
  and `S_right`, the integral means on a-to-bulk and bulk-to-b. The voltage drop
  is `integral(S dT)`; heat into the three stores is `[-I S_left T_a,
  I (S_left-S_right) T_bulk, I S_right T_b]`. These cross terms exchange equal
  and opposite electrical/thermal work and zero net entropy analytically.
  The bulk term is an interior-node Thomson allocation. Endpoint edge averages
  contain the remaining spatial allocation, so these are discrete nodal heat
  transfers, not exact continuum endpoint Peltier fluxes or a whole-leg Thomson
  source. Spatial refinement is required to resolve that distribution. A
  constant scalar retains the original zero-bulk-heat Kelvin law. Isothermal
  edges use the actual sourced S; signed integral cancellation is permitted.
  Source-cell breakpoints are retained before integration; curves on a fixed
  non-temperature axis contribute their evaluated constant coefficient.
  Relative thermopower, nonzero tabulated-only spans, source holes and
  extrapolation refuse. The entire temperature path, including the actual
  bulk state, must remain supported by the SAME selected claim.
  Original conductor, Fourier and caloric owners supply all
  dissipative transport and storage. This module supplies no independent circuit
  or thermal integrator and no tensor, magnetic-field or phase-change model.
- `liquid::resolve_liquid_state` resolves positive density, specific heat,
  conductivity and dynamic viscosity at one caller-supplied source condition.
  It retains the complete four-property bundle and derives kinematic viscosity,
  thermal diffusivity and Prandtl number, refusing non-finite or non-positive
  results. Liquid phase, formulation, concentration and operating support come
  from the selected card and query; the adapter adds no EOS, concentration
  conversion, pressure-independence or phase-stability evidence. Arithmetic
  after bounded resolution is synchronous and fixed-cost. Source-driven heat
  and flow comparisons remain numerical checks of the declared properties,
  not experimental qualification. The separate
  `liquid::resolve_liquid_state_from_kinematic` entry point requires a
  kinematic-viscosity claim instead and derives dynamic viscosity as `rho*nu`
  at the same query point. It retains the original kinematic receipt and
  refuses unrepresentable conversion; neither API silently substitutes the
  other viscosity basis or certifies consistency of independent source tables.
- `liquid::density::ResolvedLiquidDensityState` binds that liquid's exact
  selected `rho(T,p)` Surface2D claim to its same-state Cp through the shared
  `fs-thermochem::eos::density` identities. It supplies expansion, isothermal
  and adiabatic compressibility, Cv and sound speed through the existing
  `FluidAcoustics`, `FluidTransport`, `ThermalFluidTransport` and
  `ThermoviscousFluid` traits, retaining one original density everywhere.
  Original LinearInside cells supply analytic two-sided derivatives; exposed
  edges, holes, unequal incident gradients and declared T/p validity endpoints
  refuse. Typed coordinate meanings, card/claim identity and every non-T/p
  coordinate remain fixed. Refresh updates all original transport properties;
  immutable prepared source cells are shared across states. Pre-bound expansion
  or independent beta refuses, and no derived property is inserted as a source
  claim. Positive compressibility and Cv are mandatory; signed or zero expansion
  is retained. The source caps bound work. This local single-phase equilibrium
  closure supplies no derivative-uncertainty, phase-stability or finite-path
  thermodynamic-integrability certificate.
- `resolve_isotropic_elastic_state_point` requests only density, Young's
  modulus, and Poisson ratio for linear vibration/elasticity. The distinct
  `resolve_isotropic_solid_state_point` additionally requires yield stress for
  contact/plastic admission; consumers do not need to fabricate an unrelated
  yield datum merely to compute modes.
  `resolve_isotropic_thermoelastic_state_point` atomically resolves density,
  Young's modulus, Poisson ratio, linear expansion, specific heat, and thermal
  conductivity at an explicit positive `T` in kelvin. It retains all six
  property-use receipts and their original uncertainty and supplies the existing
  Zener law. Missing properties and out-of-domain states refuse; no density or
  chemistry-name preset is selected. This is the homogeneous isotropic,
  single thermal relaxation mode approximation. Thermal conductivity may be
  explicitly dimensional (legacy data) or the exact static ThermalConductivity
  schema; a semantic claim of any other kind/form refuses. The NASA TN D-6448
  seed's v3 manifest declares its source-to-runtime names, thermal kind and
  absolute-temperature axis. Typed T coordinates must be absolute temperature;
  a temperature difference cannot supply the Zener law's reference temperature.
  The compiler-to-state integration test uses that actual source at 300 K,
  retaining its Unstated uncertainty and point-only temperature validity.
  Resolved material/interface bundles expose `axis_quantities()` from their
  member receipts. Existing bundle hashes transitively bind every descriptor
  through these receipt identities without duplicating the quantity system.
  This is not an anisotropic law or a
  resolved thermomechanical energy balance. Its synthetic temperature curves
  and independent Zener checks establish numerical behavior, not experimental
  validation.
- `resolve_isotropic_thermoelastic_state_point_with_requirements` resolves the
  same mandatory six-property bundle together with explicitly requested auxiliary
  properties. One card, state point, exact pin plan and resolved identity cover
  the complete bundle; auxiliary values cannot replace or ambiguate required
  thermal/elastic names. The existing 64-property cap is checked before cloning
  auxiliary requirements. The original convenience resolver uses an empty
  auxiliary set and retains its prior requirements and identities.
  The material-driven plate compiler consumes this API to bind thermal loss and
  independently selected bending viscosity to the same physical state. G0/G1/G3
  integration tests in fs-couple cover required-name shadowing, count limits,
  exact claim pins in the presence of conflicting heat-capacity claims, typed
  conductivity, temperature-dependent input curves and actual emitted pressure.
- `resolve_isotropic_thermal_expansion_state_point` independently resolves the
  instantaneous isotropic linear expansion coefficient, including its exact
  state-point validity and property-use receipt. Negative expansion is
  admissible; material names imply neither its sign nor its value. This result
  is a coefficient, not a total strain: a coupled driver must integrate it over
  the admitted temperature path before solid mechanics use.
  `integrate_isotropic_thermal_expansion` performs that integration exactly for
  scalar, piecewise-linear curve and fixed-context surface claims. It requires
  absolute temperature (including legacy kelvin `T`), keeps every other state
  axis fixed, proves continuous support through the canonical source line,
  and pins the complete path to one selected claim. Crossed source-cell
  breakpoints preserve changes in interpolation slope. The complete path and
  signed `integral(alpha(T) dT)` result are identity-bound.
  `integrate_orthotropic_thermal_expansion` uses that same integration owner
  for the three explicit `linear_thermal_expansion_coefficient_1/2/3`
  properties in inverse-kelvin dimensions. Each principal coefficient retains
  its own reference claim through the complete fixed-context path. The result
  retains all three source records and signed principal strains; negative
  coefficients, cooling and fixed-context surface slices use the same rules.
  Missing directions, ambiguous claims, wrong dimensions and support gaps
  refuse the whole result. These are additive small strains in one fixed
  principal frame. The consumer must supply its proper rotation into the mesh;
  neither a general six-component coefficient law nor changing principal axes
  are inferred. The focused regression source is
  `tests/orthotropic_thermal_expansion.rs`.
  `resolve_mean_thermal_expansion_strain` separately resolves source-declared
  interval means as engineering strain `mean_alpha * (T_end - T_reference)`.
  Both positive absolute-temperature coordinates must be typed and constrained
  by the selected source validity. The derived result retains the parent
  evidence and binds its reference/endpoint roles in its identity; it refuses
  nonfinite results and nonpositive length ratios. Callers must explicitly
  select a mean coefficient: equal dimensions do not establish that convention.
  This endpoint consequence supplies no continuous alpha law, density law,
  thermal stress, or propagated statistical uncertainty. The stored steel
  source integration in `xtask/tests/matdb_pack_cli.rs` covers the five observed
  intervals, receipt preservation, replay and unsupported-state refusals.
  `resolve_orthotropic_elastic_state_point`
  atomically resolves density plus the nine principal-axis engineering
  constants and refuses a non-positive-definite compliance. Spatial material
  orientation remains an independently identity-bound field owned by the
  geometry/solid consumer; no chemistry name implies isotropy or a grain axis.
  `resolve_orthotropic_elastic_state_point_from_reference_ratios` derives E/G
  as the exact binary64 products of receipted E_L and the source's ratios.
  L/R/T map to axes1/2/3, with shears in G_LR,G_RT,G_TL order. The caller
  declares Full3d or an ordered material plane and the MajorRatios Poisson
  rule. Full3d returns the existing complete law; a plane result has only its
  two Young moduli, major Poisson value and shear modulus, never fabricated
  out-of-plane constants. Missing required ratios/density refuse by name.
  Independently measured minor Poisson ratios for requested pairs, when
  present, retain receipts and signed relative reciprocity residuals; they
  are never averaged into the major set. Unrepresentable residuals refuse.
  `ReferenceRatioDerivation` uses domain
  `org.frankensim.fs-material.orthotropic-reference-ratio.v1` to bind source
  bundle, set, axes, rule, strain limit and residuals. Canonical-key states
  retain their existing bundle identity and have no ratio derivation.
  Derived products are not inserted as fake source claims. FPL E_L remains
  the quasi-static reference modulus with its original moisture validity and
  Unstated uncertainty; no dynamic correction, moisture law, spatial grain
  orientation or propagated covariance is inferred. This is bounded exact
  algebra after ordinary source resolution, with no iterative cancellation
  claim. G0 tests cover compliance roundtrip2e-15 (8ulp), product/receipt bits,
  reciprocity, plane support and refusals; the source-to-consumer lane checks
  all23 FPL packs (12 Full3d,18 Plane12) and bitwise plate section/mode controls.
  Native qualification is recorded separately; these tests are not specimen
  validation. No material pack schema or source row is changed.
  `resolve_elastic_tensor_state_point` resolves density and all 36 explicitly
  addressed coefficients from one card at one supported point. The requested
  keys must share a source tensor identity, basis/frame, notation, order and
  symmetry, and address each cell exactly; missing entries and implicit zeros
  refuse. The existing scalar requirements carry the complete descriptor through
  all bulk/interface selection policies and retain each usage receipt.
  `ElasticTensorStatePoint` exposes the nominal source matrix and its evidence
  bundle. Full-law symmetry/stability and spatial rotation are admitted by
  `fs-solid::linear3::TetElasticMaterial::from_resolved_elastic_tensor`.
  Marginal uncertainties remain upstream; this adapter does not assume
  independence or produce a certified rotated uncertainty bound.
  `resolve_strain_tensor_state_point` resolves all six explicit symmetric strain
  coordinates through the same scalar requirement/selection path. Exact keys
  must share notation, order, frame and source tensor identity and address each
  position once. Negative and zero strain components are admitted; missing
  components never imply zero. `StrainTensorStatePoint` retains the six nominal
  values and all selected evidence at the query point. The solid thermal adapter
  can consume this bundle with a caller-declared integrated stress-free strain
  interpretation. Neither point support nor a source tensor name proves that
  temperature integration was performed or that a complete path is supported.
  `resolve_stress_tensor_state_point` likewise resolves six explicit symmetric
  Cauchy stress coordinates from one material card and supported query point.
  Complete scalar requirements retain the stress descriptor through all bulk
  and interface selection policies; equal pressure dimensions do not alias
  stiffness. `StressTensorStatePoint` exposes pascal values in the declared
  notation/order/frame, source tensor identity and the selected evidence bundle.
  Missing entries, mixed tensors/bases and context-free requirements refuse.
  Signed stress is admitted; engineering stress shear is physical sigma_ij.
  Source compiler v8 and normalized-pack v7 feed this path. The synthetic
  Pa/kPa import test compares it with the existing rotated isotropic law's
  constrained thermal response. It adds no solved field, initial-stress load,
  inferred covariance, finite-deformation stress or experimental validation.
  `resolve_joint_stress_tensor_state_point` resolves the same six nominal
  stress coordinates and their joint receipt from one immutable material-card
  pack, point and selection policy. Each nominal member receipt must match the
  joint answer. Private `JointStressTensorStatePoint` fields bind both under
  `joint-stress-state-point.v1`; a caller cannot substitute unrelated covariance.
  The source-coordinate covariance or exact unknown-correlation reason is
  retained. Nominal values are not asserted to be distribution means, marginal
  widths are not converted to covariance, and no orientation uncertainty is
  inferred. The existing two joint-query policies are supported. The fixed-frame
  numerical transform is owned by fs-solid, with no new pack schema or dependency.
  `resolve_joint_strain_tensor_state_point` binds that nominal state to the
  existing joint query on the same `NormalizedMaterialCardPack`, component
  requirements, physical point and selection policy. It compares all member
  receipts and retains the source joint receipt, including covariance in source
  component order or the exact unknown-correlation reason. Its private
  `JointStrainTensorStatePoint` fields prevent injection of an unrelated matrix.
  Domain `joint-strain-state-point.v1` binds nominal state and joint receipt.
  The existing single-claim and observation-backed joint policies are supported;
  this adapter adds no pinned joint policy or marginal-width-to-variance inference.
  The two-material reduced modulus on the contact-capable isotropic rung uses
  the standard ordered Hertz half-space compliance and binds both resolved
  state identities. Out-of-range temperature or phase-dependent data refuse
  through the material query; this module does not silently extend a solid
  card through melting or select a new phase/model.
- `phase::EquilibriumEnthalpyPhaseCurve` is the generic phase-state ingress:
  source-provided monotone specific-enthalpy knots resolve temperature, bulk
  density, and solid/liquid mass fractions, including an isothermal
  latent-heat plateau. It contains no material-name presets and never
  extrapolates. `try_single_phase` admits only explicitly all-solid or
  all-liquid curves with strictly increasing temperature; it refuses a
  latent plateau or mixed fractions disguised as a one-phase curve. It is not
  a heat-transfer, deformation, remeshing, free-surface, acoustic, or optical
  solver.
  `temperature_derivative_at_specific_enthalpy` returns the admitted chart's
  piecewise-linear dT/dh, exactly zero on a latent plateau. Interior knots use
  the right segment; the upper endpoint uses the left. This is an explicit
  generalized derivative for semismooth transport, not a smoothed heat capacity
  or a classical derivative at a corner. Nonfinite/out-of-domain inputs and
  unrepresentable positive slopes refuse. Signed-zero enthalpy endpoints resolve
  identically in the state and derivative queries. The focused public derivative
  tests include independent segment differences, latent plateaus, knot sides,
  single-phase charts, domain errors and representability boundaries.
  `try_from_heat_capacity` derives that same one-phase chart from positive,
  strictly increasing source `Cp(T)` knots, taking its enthalpy reference at
  the first source temperature. Linear `Cp` segments are integrated exactly
  and subdivided under a caller-provided real-arithmetic chord-error bound;
  all source endpoints remain knots. Sensible samples retain the source
  piecewise-linear rho(T) at the evaluated temperature, including constant-Cp
  spans that require no caloric subdivision. The density interpolation basis
  is bound into chart identity. Generic supplied h/rho charts retain their
  original specific-volume interpolation and identity. The enthalpy bound
  excludes source uncertainty and floating-point error.
  `density_law_on_interval` exposes that original interpolation through
  `PhaseIntervalDensityLaw`: affine specific volume for generic/coexistence
  intervals, or affine density for sourced sensible branches. The query uses
  one adjacent original interval and refuses invalid indices; it neither fits
  a replacement EOS nor changes chart identity. Conservative thermal-volume
  consumers use this distinction to account for nonlinear mixing volume.
  `mixing_specific_volume_m3_kg` evaluates the original-chart remainder
  `v(hr) + v'(hr)*(hd-hr) - v(hd)` across sensible and latent boundaries.
  Each crossed affine-density span uses a stable reciprocal quadratic
  remainder; affine-volume spans contribute their original slope changes.
  Signed knot terms and local enthalpy differences retain the physical mixing
  correction without subtracting large caloric datums or refitting density.
  Interior knots use their right derivative; the upper endpoint uses its left.
  The query allocates/hashes nothing, preflights its crossed-interval budget,
  returns the actual count for whole-solve accounting, and polls caller
  cancellation before/after evaluation and every 64 spans. Equal enthalpies
  consume zero intervals. Original domain and representability refusals apply;
  this constitutive query owns neither material transport nor a heat source.
  Generated charts cap at 4096 knots before allocation. This is bounded
  numerical construction, not experimental validation or a mechanical model.
  `ResolvedHeatCapacityCurve::try_new` admits that chart directly from one
  immutable material card using explicit Cp/density keys and complete typed
  endpoint queries. Both properties require continuous support over the same
  fixed-context temperature span through `ClaimSet::query_line_typed`.
  The union of their curve knots and crossed surface-cell breakpoints is
  retained (at most 4096 distinct temperatures), including density-only knots. Ordered
  Cp/density receipts remain accessible and are bound into the curve identity.
  Wrong dimensions, changing non-temperature conditions, unsupported spans and
  discrete-only data refuse. The reference enthalpy and single-phase choice
  remain caller assumptions; this ingress does not infer phase stability,
  phase-transition endpoints, source accuracy or material compatibility.
  `UniformEnthalpyStepInput` and `UniformEnthalpyStep<R>` carry a proposed
  uniform-body heat-transfer step between physical owners without a reverse
  solver dependency. They include specimen mass/geometry, initial chart/state,
  duration, internal heat, actual boundary heat, a distinct solve-residual
  budget and the concrete solver report. They are data, not an admission or
  transport solver. The transport must justify lumping; the coupled owner must
  verify chart identity, energy balance and regime before publishing, and the
  callback must not debit an external source before that commit.
- `resolve_interface_state_point` applies the same atomic query contract to an
  ordered `InterfaceSystemCard`. Its identity therefore includes both surface
  material states and texture frames, medium, third body, environment,
  preparation/history, complete query point, claim-selection policy, selected
  values, and property-use receipts. Interface properties never migrate onto
  either bulk material merely because the two chemistry names match.
- `ResolvedLinearInterface` lowers explicit stiffness/viscosity-per-area claims
  from that same resolved interface into passive bilateral `k` [N/m] and `c`
  [N s/m] over a declared positive area. It retains the complete ordered source,
  selected property names and area; exact receiving material states and actual
  angular frequencies can be checked without replacing receipts. Both claims
  require explicit intersected frequency support and fixed-temperature constant
  coefficients. Missing loss, frequency-dependent curves, two zero coefficients,
  negative values, overflow and disappearing nonzero products refuse. Geometry,
  work-conjugate force maps, rest displacement and stepping remain the consumer's
  existing owners; this is shared material lowering, not another joint solver.

- `SmallStrainLaw` trait (Voigt 6-space, TENSOR shear components):
  `stress` / `tangent` / `update_state` / `admissibility` / `card`.
  **The tangent contract**: `tangent` is the exact derivative of the
  ALGORITHMIC stress update at the same committed state — FD-gated for
  every law in conformance (merge-gate discipline, same as adjoints).
- `IsotropicElastic` (E, ν → Lamé), `OrthotropicElastic` (engineering
  constants; construction REFUSES thermodynamically inadmissible Poisson
  sets via compliance positive-definiteness minors).
- `J2Plasticity`: radial-return mapping with linear isotropic or explicitly
  split isotropic/Prager kinematic hardening
  and the Simo–Hughes algorithmic moduli
  `C = κ I⊗I + 2μθ I_dev − 2μθ̄ n̂⊗n̂`.
- `Hyperelastic` (`NeoHookean`, `MooneyRivlin`, decoupled `Ogden`,
  `HolzapfelGasserOgden`): stored energies written
  ONCE generic over the fs-ad `Real` scalar; `piola` is the exact dual
  gradient, `tangent` the exact nested-dual Hessian (9×9). `det F ≤ 0`
  and nonfinite states refuse structurally. Ogden uses the original
  `mu/alpha` convention: infinitesimal shear `mu = sum(mu_p alpha_p)/2`.
  At most eight finite terms with `mu_p alpha_p > 0` are admitted, including
  negative-mu/negative-alpha terms. The shared principal-stretch helper uses
  the existing symmetric Jacobi owner; a spectral trace lift carries the
  exact value/gradient/Hessian through repeated roots, using the analytic
  divided-difference limit without perturbing stretches or differentiating
  eigenvectors. Derivatives above second order are not claimed.
  HGO admits at most four explicitly unit reference directions, `c10 > 0`,
  `k1 >= 0`, `k2 > 0`, dispersion in `[0,1/3]`, and an explicit full or
  isochoric structural invariant. Its energy per family is
  `k1/(2 k2) (exp(k2 <E>_+²)-1)`. At zero dispersion `I4 <= 1` contributes
  exactly zero; with dispersion, tension refers to the distributed
  structural strain `E`, preserving C1 continuity and direction independence
  at dispersion 1/3. The inactive-side Hessian is returned at E=0; a
  two-sided Hessian there is not claimed. A single aligned family supplies
  `2 k1 (lambda²-1) exp(k2 (lambda²-1)²) lambda²` incompressible Cauchy
  stress; two identical families supply twice this value.
  For `k2 E² <= 1/8`, the fibre potential factors as
  `k1 E²/2 * exprel(k2 E²)`, with the analytic `exprel(0)=1` extension
  evaluated by the original scalar-generic exponential series. This retains
  finite quadratic fibre storage, stress and tangent as positive `k2` tends
  to zero, including subnormal exponents, without constructing `k1/k2`.
  Zero-strength fibres return zero before the exponential. Focused regressions
  in `hyperelastic_ogden_hgo.rs` compare the quadratic energy, all stress
  components and the full tangent against an independent polynomial oracle
  at finite and acoustic-scale strains, exercise both invariant formulations
  and dispersion, and check tangents by finite differences.
  All four variants enforce all principal stretches in `[1/limit,limit]` in
  stress and tangent evaluation, with eight-EPSILON eigensolver boundary
  roundoff; limit*(1+1e-12) refuses. The limit must be finite and at least one,
  including after mutation of the public field. Generic-scalar energy
  evaluation retains its documented admissible-F precondition. The Ogden and
  HGO variants offer quadratic or log-squared
  bulk energies and separate distortional/bulk accessors. Full-invariant
  HGO's fibre accessor intentionally retains volume dependence.
  Compressible Neo-Hookean energy retains tiny elastic storage near any rigid
  rotation through the objective increment `A=F^T F-I`. Its small-increment
  branch evaluates `trace(A)-log(det(I+A))` without subtracting unit-size
  invariants; a bounded atanh series keeps its value, gradient and Hessian
  consistent. The ordinary invariant expression remains outside that numerical
  branch. Full and split Neo-Hookean energy use the same helper. This addresses
  the zero-energy cancellation of tiny simple shear without changing the
  physical potential or the original stretch domain.
  Ogden also retains small elastic storage through `L=dev(log(C))`, evaluated
  from the objective increment with bounded matrix-log and exponential series.
  Cayley-Hamilton scalar recurrences avoid an eigensolve and heap allocation
  on this branch. Every term starts at quadratic strain order; factoring the
  original `mu*alpha` preserves positive/negative and near-zero exponent
  limits without cancellation. The series use the same `Real` scalar for
  energy, stress and tangent. Admission requires `C-I` entries at most `1/64`,
  a normalized traceless strain infinity norm at most `1/16`, and an exponent
  bound guaranteeing `||alpha L/2||_inf <= 1/8`; the original spectral owner
  remains outside those numerical bounds. This preserves the finite-strain
  law and stretch envelope, including pure dilation and finite rotations.
  Four focused G0/G1 regressions in `hyperelastic_ogden_hgo.rs` check independent
  simple-shear work/stress/tangent down to `1e-12`, signed and near-zero
  exponents, dilation/rotation, branch continuity and all tangent columns.
- `hyper::acoustic` prepares the original full finite-strain tangent once and
  evaluates current-direction wave modes without repeating constitutive AD.
  It uses the 3-D compressible law, not a plane-stress/incompressibility
  projection or a Maxwell algorithmic tangent.
  With `q=F^T n`, it uses `Q_ik=A_iJkL q_J q_L/J` and conserved current density
  `rho=rho0/J`. The original Hessian already includes initial-stress effects;
  no small-strain replacement or extra geometric term is added. The shared
  Jacobi eigensolver returns signed squared phase speeds and polarizations;
  negative or zero eigenvalues retain their sign and have no positive real
  phase speed. The physical tensor, current density and eigen residual remain
  available. Positive modes also return spatial energy velocity
  `v_r = a_i A_iJkL a_k q_L F_rJ/(rho0*c)` from that same cached Hessian,
  with no additional AD or eigenproblems. It equals the frequency gradient on
  a smooth separated branch; at a repeated root it retains the chosen
  polarization's flux. Nonpositive modes have no real energy velocity.
  Stable modes in one direction do not certify global ellipticity or an
  attenuation/dispersion model.
  `hyper_acoustic.rs` checks independent finite-NH formulas, the small-strain
  limit, all four laws against rank-one work derivatives, anisotropic rotation,
  density/modulus scaling and an admitted expanded NH state with a negative
  longitudinal eigenvalue. Energy velocities are checked against independent
  finite-NH dispersion, all four frequency gradients, current-frame rotation,
  normal projection and density/direction scaling. The formulation follows
  [Abiza, Destrade & Ogden (2012)](https://doi.org/10.1016/j.wavemoti.2011.12.002).
  NH's former unconditional polyconvex declaration is restricted to the
  sufficient domain condition `lambda=0` or
  `3 ln(limit) <= 1+mu/lambda`; the original log-volume energy can lose
  rank-one convexity beyond it. MR retains a sufficient certificate when
  `c01=0`; other MR parameters return unknown, because the existing
  `|cof F|²/J^(4/3)` term lacks that convex-perspective certificate. Unknown
  does not assert non-polyconvexity, and near-reference random samples do not
  establish a domain-wide certificate.
- `hyper::plane_stress` relaxes all three entries of the transverse deformation
  column until the original law's transverse nominal tractions vanish. Thin
  planar hyperelastic consumers therefore retain thickness contraction and,
  for tilted anisotropic fibres, transverse shear relaxation. The original
  reference-volume energy and in-plane Piola stress remain work conjugate;
  the consistent tangent eliminates the transverse block by its Schur
  complement. Bounded Newton/backtracking uses the original principal-stretch
  domain and refuses nonpositive transverse stiffness or exhausted budgets.
  If unit thickness itself violates the material domain, initialization first
  tries `||F_in_plane||_F/sqrt(2)`, the positive thickness minimizing `I1bar`
  at zero transverse shear. This retains admissible compressible Gent states
  whose chain limit excludes the unit-thickness guess; the full traction solve
  still determines thickness and shear. This candidate consumes one of the
  same `max_backtracks` seed budget, and every original material bound remains
  active. Independent scalar Gent traction roots under biaxial extension and
  compression, all three bulk laws, and reduced work/tangent differences cover
  this recovery in `hyper_plane_stress.rs`.
  If an admitted initial thickness has nonpositive transverse stiffness,
  initialization tries a volume-preserving thickness. An out-of-domain target is approached by
  bounded admissibility bisection until a stable seed is found; the same
  `max_backtracks` cap bounds this search. No accepted deformation or material
  stiffness is clipped. Large biaxial extension is covered for both admissible
  and inadmissible isochoric seeds, plus compression and decoupled bulk laws.
  Traction normalization uses inverse transverse-compliance trace, so a large
  bulk penalty cannot hide a shear-scale residual. The energy/stress-only path
  forms only a three-direction local AD Hessian; the full path forms the native
  9x9 Hessian once after convergence. The tests in `hyper_plane_stress.rs`
  cover analytic isotropic plane stress, all four original energy models,
  tilted fibres, work/tangent finite differences, frame and modulus scaling,
  near-incompressible traction closure, and domain/budget refusals. This is a
  local stable plane-stress branch; no unique global minimum, shell bending,
  wrinkling, exactly incompressible law, or locking-free volume formulation is
  claimed. The formulation follows the quadrature-local elimination of
  [Ahmadi et al. (2025)](https://doi.org/10.1007/s00707-025-04310-z), extended
  here to the complete transverse column by differentiating `P[:,2]=0`.
- `hyper::incompressible_plane_stress` is the explicit exactly incompressible
  alternative: `F22=1/det(F_in_plane)` enforces `J=1`, and the two transverse
  shears relax against the original `isochoric_energy` remainder. No bulk
  penalty or large bulk-modulus arithmetic is evaluated. The full physical
  Piola stress is `dW_remainder/dF - p cof(F)`, with the compression-positive
  multiplier `p=(dW_remainder/dF22)/det(F_in_plane)`. This multiplier follows
  the native remainder's convention; for NH it is `mu*(F22²-1)`, and it is
  not generally the negative mean Cauchy stress. The returned tangent
  differentiates the six-variable constrained energy, including the second
  derivatives of the exact thickness, then condenses the two shear variables.
  The state-only path uses two-direction local AD. Original stretch admission,
  positive shear stiffness, unit-covariant traction tolerances and bounded
  Newton/backtracking apply without clipping. `hyper_incompressible_plane_stress.rs`
  tests independent NH finite-deformation formulas, all four laws and both HGO
  invariant choices, physical transverse traction closure, energy/tangent finite
  differences, rotation, modulus scaling and bitwise independence of bulk
  parameters through 1e300. This locally stable planar branch does not claim
  3-D mixed displacement/pressure elements, shell bending, wrinkling or a
  globally unique equilibrium.
- `hyper::cards` instantiates strict, dimensioned `ConstitutiveModelCard`s
  for `neo-hookean-log`, `mooney-rivlin-decoupled`, `ogden-decoupled` and
  `hgo-dispersed` v1 and registers the real
  stress/tangent/energy factories in `graph::LawRegistry`. Graph inputs are
  H=F-I, outputs P, so the zero-input admission probe is the reference state.
  Moduli require Pa, exponents/directions/counts/selectors are dimensionless;
  missing/extra parameters, fractional counts, versions and undeclared
  invariant selectors refuse. FSMODPK v1 remains unchanged.
  Neo-Hookean supplies `mu` and Lamé `lambda` in Pa for the native
  `mu/2 (I1-3) - mu ln(J) + lambda/2 ln(J)^2` energy. Mooney–Rivlin
  supplies `c10`, `c01` and `kappa` in Pa for its native isochoric-invariant
  energy plus `kappa/2 (J-1)^2`. Both require dimensionless `stretch_limit`;
  neither accepts an unused bulk selector. Native admission requires finite
  coefficients and representable reference shear/longitudinal stiffness.
  All four laws traverse the existing model-pack decode and registry test;
  the isotropic laws and Ogden declare smooth nodes, HGO piecewise smooth.
  `hyper::cards::ResolvedHyperelastic` retains the exact source model card and
  admitted `QueryPoint` beside its native law. `try_from_model_card_at` requires
  the typed source domain to contain that fixed query; bounded conditions
  support static temperature, pressure, composition, moisture and concentration
  quantities, plus legacy kelvin `T`. The context-free constructor admits only
  unconstrained domains. Evolving kinematic calibration axes refuse until an
  actual coordinate owner exists. `HyperelasticMaterial` lets spatial consumers
  borrow this same energy/Piola/tangent owner. Native parameter and actual
  principal-stretch checks remain unchanged; the wrapper adds no thermal
  evolution or uncertainty propagation.
  The offline `xtask matdb-pack` `constitutive-v1` profile accepts explicit
  model/law/state/initialization/validity declarations and dimensioned
  parameter/validity rows, normalizes through fs-qty, retains each transform
  receipt, and self-decodes the existing runtime model pack. Its compiler
  tests transport the Ogden/HGO laws into the real registry. A generic profile does
  not certify any law's physical parameters or infer a tissue identity.
  Polyconvexity is parameter-derived: Ogden sufficient exponents are
  alpha>=3/2 (F power perspective) or alpha<=-3 (cof F perspective), with
  positive coefficients and convex bulk. Full-invariant HGO has convex
  positive-part structural energy. Other sets return `None`; no numerical
  spot check upgrades an unknown result to a proof. Log-squared bulk is
  certified only on declared J<=e, where its second derivative is nonnegative.
- `active_tissue::ActiveTissue` composes a passive `Hyperelastic` with a
  prescribed `ActiveFibreStress` along a reference fibre `a0` (unit; for an
  HGO passive card it must be one of its fibre families). The active Cauchy
  stress is `T_max a f_L(lambda_f) m (x) m` with `lambda_f = |F a0|` and the
  activation `a` in `[0, 1]`; its first Piola image is
  `(J / lambda_f^2) T_max a f_L (F a0) (x) a0`. The factor `f_L` is a Gaussian,
  a parabola or a tabulated piecewise-linear family (at most 64 knots). The
  law refuses a fibre stretch outside its declared inclusive domain (eight ULPs
  of roundoff at the bound) and an activation outside `[0, 1]`
  (`MaterialError::State`), and a non-positive or non-finite `T_max`, a non-unit
  fibre or an invalid family (`MaterialError::Parameters`). Zero activation
  returns the passive stress and tangent bit for bit and does not consult the
  domain. The active stress is a PRESCRIBED, non-conservative input: the tangent
  is the passive (symmetric) tangent plus the exact forward-mode derivative of
  the active stress, which is not symmetric (`tangent_symmetric: false`), and
  its work `P_a : Fdot` is booked as an external power input, never as stored
  energy. Activation is a performance control that no law input, card or pack
  sets here.
- `WoolFelt::response` evaluates engineering compression strain with a finite,
  inclusive `eps_densify` validity bound and returns stress, the consistent
  tangent and proposed compaction history together. It readmits public
  parameters and verifies the retained maximum against this law's power-law
  envelope, including during unloading and tensile-gap opening. A trial beyond
  densification, nonfinite input, or inconsistent prior excursion refuses;
  the law does not extend its domain by clipping strain or resetting history.
  The infallible `Uniaxial` stress/tangent methods return NaN on refusal and
  its state update preserves the input. The analytic `envelope` remains
  available without history/domain admission. `wool_felt_strain.rs` checks
  inclusive bounds, transactional refusal, residual-gap reloading, consistent
  derivatives and the analytic dissipated-work/conservative-cycle limits.
- `Uniaxial` trait + the RC flagship pair: `MenegottoPintoSteel`
  (R0/a1/a2 curvature degradation, Bauschinger via branch-state
  asymptote intersection) and `ManderConcrete` (confined envelope
  `f = f′cc·x·r/(r−1+xʳ)`, elastic unload/reload lines with residual
  strain, zero tension). `ManderConcrete::response` readmits finite positive
  parameters with representable `Ec > Esec` and `r > 1`, and admits the
  committed maximum against that same envelope before producing stress,
  tangent and a proposed history. The ultimate compression strain `eps_cu`
  is inclusive; exceeding it or providing a corrupt/nonfinite history
  refuses without mutation. At zero virgin strain and at residual-gap closure
  the tangent is the compression-side derivative, while an open tensile gap
  has zero stress and tangent. This unilateral convention supplies no tensile
  capacity. The infallible trait stress/tangent methods return NaN on refusal;
  its history update retains the accepted input. Consumers needing a diagnosis
  use the checked response. `envelope` remains the analytic curve evaluator,
  without history or domain admission. Focused G0 tests in `conformance.rs`
  verify the rational `r=2` envelope, one-sided startup, gap reloading,
  the inclusive limit and invalid history/parameter refusal. The formula is
  Mander, Priestley and Park (1988), Eqs. (3), (4), (6), (8),
  doi:10.1061/(ASCE)0733-9445(1988)114:8(1804); the elastic cyclic lines remain
  the declared simplified rule. No post-crushing damage evolution is supplied.
  Menegotto–Pinto's virgin state has direction zero. Its first nonzero
  excursion selects the signed yield asymptote while retaining `R0` and zero
  reversal history; starting in compression is not a prior plastic cycle.
  Zero-strain commits preserve that virgin state. Focused G0/G3 conformance
  tests check the analytic yield ordinate, zero/subnormal first increments,
  mirrored cyclic stresses and consistent tangents. The shared steel-section
  regression checks reflected axial and bending histories. This adds neither
  a stored-energy potential, fatigue law nor material calibration.
  Reversals retain separate tensile/compressive extrema, seeded at ±yield
  strain. Curvature uses the new asymptote intersection and the prior extremum
  in the new loading direction; `xi_max` is only a diagnostic maximum.
  Subyield unloading therefore retains the near-elastic knee, while nested
  loops recover their directional curvature. Independent no-isotropic-shift
  Steel02 values check stress, tangent and the native section response.
- `calibrate_bilinear`: segmented least squares recovering (E, σ_y, H)
  from monotonic data with standard-error envelopes and RMS residual.
- `evidence_stress`: wraps any law's stress in `Evidence` with the card
  attached and `in_domain` FLAGGING calibration-domain exit.
- `tensor`: Voigt helpers (deviator, contraction with shear doubling,
  von Mises, Rodrigues rotations) used by the objectivity gates.
  Shared `StressTensorBasis` / `StressTensorNotation` descriptors distinguish
  physical/engineering stress shear (both `sigma_ij`) from Mandel stress shear.
  They are re-exported from fs-matdb for fs-solid's recovered thermal-stress
  outputs; they do not reuse engineering-strain factor two or imply Piola stress
  or material/frame calibration. Source-pack stress transport uses the separate
  StressTensorComponent context described above.

### `gas` — first-principles ambient gas state (zdmm1 slice 1)

`thermochemical_gas::packs::from_packs` joins explicit verified species,
NASA-9 and pure-transport material-card packs. Names supply no properties.
The model and each resolved state retain the three pack hashes per canonical
component, alongside existing thermodynamic/transport receipts. Missing NASA-9
or transport data refuse. `gas_species::input::read_gas_record` requires an
already-loaded registry for thermochemical input and returns the full resolved
state; the legacy two-argument reader preserves its old bits and refuses that
record without a registry. A named `mixture` record keeps its species and
mole fractions beside its state (`GasRecord::Mixture(state, composition)`;
a mass record is converted by the species' own molar masses), so a consumer
that needs the composition reads it and never infers it from the state.
Phase admission is explicit: frozen gas makes no
condensation claim; the liquid-water option applies the existing Buck limit
only inside its stated temperature window. No real-gas or molecular-relaxation
accuracy is supplied by this adapter. Literal NASA/JANAF disagreements remain
test failures after accounting separately for each printed half-unit. Relative
bands are unchanged; cp also satisfies its relative-only gate. Rows absent from
the source are recorded as absent, never interpolated into alleged facts.

`gas_species::input::temperature::IsobaricGasLaw` retains the original pressure,
composition and property laws for local temperature feedback. Legacy humid air
preserves water mole fraction rather than relative humidity and refuses later
supersaturation. Named mixtures keep their original source windows. Canonical
thermochemical packs enter through `read_source_temperature_law`, retaining the
NASA-9 and transport models, typed coordinates and pack identities at every
evaluation; this evolving path requires explicit frozen-gas phase admission.
Density, heat capacity, viscosity and conductivity come from that same local
state. `GasRecord` implements `FluidAcoustics` and `ThermoviscousFluid` alongside
its existing flow and thermal transport interfaces. Generic ducts, radiation,
FE cavities and aeroacoustic consumers can borrow the complete parsed or
temperature-refreshed record directly, retaining composition and NASA-9/transport
receipts. Each view delegates to the original owner and its complete coherence
admission; no new gas law or default atmosphere is introduced. These are frozen
local snapshots, with explicit temperature updates owned by `IsobaricGasLaw`.
No pressure evolution, species transport,
compressible inertia, phase equilibrium or experimental validation is implied.
G0/G3 unit tests cover fixed-composition EOS scaling, viscosity/conductivity
response, humid cooling refusal and source/domain preservation.

`IsobaricGasLaw::transport_temperature_tangent` retains that same evaluated
`GasRecord` and returns analytic isobaric derivatives of density, Cp, conductivity
and dynamic viscosity for the original engineering and thermochemical gas models. The derivative
holds pressure and mole fractions fixed, including water. It differentiates
Sutherland viscosity, the original USSA/Eucken conductivity, and the complete
Wilke/WMS quotient, including each temperature-dependent interaction denominator.
Ideal-gas density gives `-rho/T`; calorically perfect models' frozen Cp gives zero. Values and
phase/domain admission remain with `at_temperature`, and mixture interaction
values use the original native function. Additional derivative storage is a
fixed stack array for at most six engineering species. No finite-difference material query,
replacement transport fit or experimental accuracy claim is introduced.
At a validity endpoint the tangent belongs only to the admitted branch.
Thermochemical Cp derivatives come from the original NASA-9 region owner,
weighted by the retained mole fractions and converted with the same mixture
molar mass. Pure transport tangents differentiate the actual selected constant,
linear-temperature or triangulated temperature/pressure claims; constants require explicit nonzero thermal
validity. Source-domain endpoints use their inward derivatives. Interior curve
corners and NASA-9 Cp jumps/corners refuse unless adjacent values/slopes agree
to rounding tolerance. Tabulated-only points and unsupported source forms
retain their value admission without gaining derivative support.

Original `Surface2D` transport holds pressure exactly fixed and accepts either
source axis order. The source owner's exact cell predicates admit the local
thermal slice and its inward endpoints, including pressure-boundary states;
claim validity clips the permitted directions. All contributing thermal slopes
must agree to relative rounding tolerance. Genuine corners, holes and isolated
samples refuse without promoting value-only or tabulated data. Neither a
pressure derivative nor an unowned changing coordinate is inferred.
`IsobaricGasLaw::prepare_transport_temperature_tangent` returns a borrowed
`PreparedIsobaricGasTransport` retaining the actual selected surface geometry
and gradients for repeated evaluations. Each value still comes from the same
complete gas/source query, and the selected claim identity is checked against
the retained preparation. Scalar/curve and engineering paths need no cache
allocation. The convenience tangent API prepares surfaces on demand; iterative
consumers retain the view to avoid repeated geometry validation and copying.

The already evaluated native Wilke mixture supplies its analytic viscosity and
interaction-denominator derivatives. WMS conductivity differentiates its own
selected pure conductivity and reuses those same denominators and derivatives.
Additional work is bounded by the original 128-species admission, uses fixed
stack component arrays plus one O(N) denominator vector, and builds no N-by-N
matrix. Derivative work performs no extra material queries, source reselection,
model clone or altered caloric reference. The complete source-bearing gas record
and original pressure/composition/context remain unchanged.
Focused G0/G1/G3 tests compare original values, independently differentiated
pure-gas equations, high-order differences of the original mixture evaluator,
pressure scaling, component permutation, variable NASA-9 Cp, unequal pure
transport slopes, source knots, explicit constant support, saturation and
source-boundary refusal. These derivatives add no model accuracy or phase claim.

`IsobaricGasLaw::specific_enthalpy_change` supplies the retained gas's relative
isobaric h increment without subtracting large formation enthalpies. Legacy
and named engineering gases use their original constant mixture Cp; canonical
thermochemical gases use the original NASA-9 difference. Both endpoint states
retain source/phase admission, and crossed NASA source gaps refuse. Contiguous
raw-region reference jumps remain authored; only an explicitly reconciled
source integrates a continuous piecewise Cp branch. The operation preserves
pressure/composition and refuses nonfinite or nonpositive-Cp-directed changes.

Frozen gas parcels and nozzles use the thermochemical owner's direct entropy
differences instead of subtracting absolute reference entropies. The relative
density/entropy closure targets a requested change from the retained initial
state and preserves the original bracket, cancellation and source checks.
Reversible volume changes target zero entropy change. Midpoint advancement
reports its actual numerical entropy change without asserting exact adiabatic
reversibility. Iterative entropy and nozzle trial states use allocation-free
caloric properties; accepted states retain complete source receipts.

The converging `nozzle_flow_to` closure accepts vacuum back pressure without
querying an EOS at zero pressure. Sonic and kinetic-energy residuals use the
thermochemical owner's direct enthalpy differences; common formation offsets
cannot erase available nozzle energy.

`nozzle::rocket` adds ideal converging-diverging propulsion using the same
frozen NASA-9 composition. A bounded sonic solve determines mass flow from
chamber state and throat area, then a bounded supersonic temperature solve
closes exit mass flux using source enthalpy and entropy differences. Trial
caloric evaluations allocate no material receipts; accepted throat and exit
states retain full source and transport resolution. Explicit mass, energy and
entropy residuals are checked before publication. Outputs include exit Mach,
momentum and pressure thrust, total thrust, effective exhaust velocity,
characteristic velocity and specific impulse. The latter uses standard
gravity only as a unit convention. The off-design extension below admits
source-supported unchoked discharge, ideal internal normal shocks and attached
external compression in addition to matched/underexpanded and vacuum flow.
It does not infer viscous separation, combustion equilibrium, real-gas
validity, nozzle losses, adjoints or resolved multidimensional flow.

`nozzle::blowdown` removes finite mass from a frozen, well-mixed, adiabatic
fixed-volume reservoir on the original density/entropy path: `d(mu)=h dm`
implies `ds=0` for the remaining gas even with temperature-dependent heat
capacity. It exports actual `U0-U1` once, retaining formation conventions and
every species. `rocket_blowdown` closes `dm=dt*mdot(m0-dm/2)` by bounded
implicit midpoint iteration and reports midpoint thrust impulse. Initial,
midpoint and endpoint nozzle/source domains and the donor-fraction limit are
checked transactionally. Smooth time accuracy is second order, as checked
against analytic perfect-gas mass, pressure and impulse refinement; a mass
fraction limit does not certify time error. No wall heat, moving tank,
combustion, resolved chamber/jet dynamics or double-counted exhaust kinetic
energy is supplied by this lane.

Source-backed gas temperature and pressure support comes from the selected
NASA-9 regions and every pure-transport claim, not from the legacy USSA/Sutherland
air constructor's 50..2000 K and 100 bar bounds. `GasState::validate_coherence`
checks positive finite coordinates and the common EOS/caloric/transport/acoustic
identities; it does not assert source coverage. Energy and entropy inversion
accept positive finite temperature brackets and retain their original source
checks, residual tolerances, cancellation and atomic updates. Consequently
source-supported hot gases can heat, compress and supply the same density,
viscosity, conductivity and sound speed to existing consumers without changing
their physics code. Legacy `GasSpec` constructors keep their original limits.
G1/G3 manufactured tests cover heating through 2000 K, variable-cp reversible
compression, data-only composition changes and unchanged state after leaving
transport support. The NASA-9 family covers species-specific high-temperature
regions (NASA/TP-2002-211556, <https://ntrs.nasa.gov/citations/20020085330>);
the active source bounds still govern each calculation. This does not admit
condensed phases, real-gas behavior, equilibrium chemistry or plasma physics,
and source support alone is not an experimental-accuracy claim.

`gas_species::input::read_isochoric_parcel` retains the same parsed
thermochemical model, typed conditions and canonical source associations in
the existing `IsochoricGasParcel` owner. Its mass is rho*V from the declared
volume and initial state; later signed heat evolves the original NASA-9 energy,
EOS and transport rather than a frozen initial cp. `input::packs::read_source_parcel`
uses the same bounded canonical pack loader and complete-record checks as the
snapshot reader, and the parcel owns its sources after caller bytes are dropped.
Borrowed prefix input is unchanged on every refusal and retains a caller-owned
suffix on success; the complete pack reader rejects trailing fields.

Finite-parcel ingress explicitly requires `FrozenGas`: an initial liquid-water
saturation check cannot certify later cooling, so `LiquidWaterLimit` is refused
instead of silently weakened. The caller remains responsible for a valid
frozen-gas phase domain during evolution. No condensation, reaction, spatial
flow or experimental accuracy is supplied by this input adapter. The existing
snapshot grammar, source selection and phase checks are shared without changing
their numerical evaluation. Focused G1/G3/G4 regressions in
`fs-couple/tests/gas_source_packet.rs` check source retention after input drop,
nonlinear caloric heating and pressure in both composition bases, exact snapshot
parity, and transactional source/volume/phase/suffix refusals.

`gas_field::AxialGasField` stores 2..=4096 supplied axial stations with
temperature and canonical mole composition at one explicit common pressure.
It interpolates only primitive coordinates; NASA-9 thermodynamics, EOS, Wilke
viscosity and WMS conductivity are reevaluated by the existing owner at every
sample. Uniform fields and original station values preserve direct-evaluation
bits. No extrapolation, invented transport data or independent interpolation of
rho/c/cp/mu is allowed. `with_composition` retains the active species' original
laws, transport selections and pack identities; explicit zeros omit inactive
species without fabricating trace concentrations. The full source model remains
on the field. Only generated interpolation roundoff receives a bounded simplex
correction; authored fractions remain subject to the existing exact-sum rule.
The existing liquid-water phase policy is explicit; supersaturation refuses
instead of silently removing water. The field is prescribed, not a heat/species
transport or condensation solver, and source-station measurement uncertainty
remains with its caller. Construction and sampling are synchronous bounded
operations; callers can chunk sampling. Focused tests cover uniform bit identity,
fresh primitive evaluation, zero-concentration ingress, source retention and
station/sample/domain/phase refusals. Native verification is recorded separately.

The five-species acoustic fixture binds NASA/TP-2002-211556 Appendix D to
primary zero-density transport correlations: Lemmon/Jacobsen (N2/O2/Ar),
Laesecke/Muzny and Huber (CO2), and IAPWS R12-08/R15-11 (H2O). The stored
linear curves use 256 uniform intervals plus declared query knots; tests
measure 64/128/256 rungs and require relative representation error <= 1e-4.
Source uncertainty remains separate, without invented confidence intervals.
Water's 250 K dilute-vapor continuation has no quantified sub-triple accuracy
certificate. These curves do not certify finite-density or mixture errors.
The SRM1720 fixture uses the pack's information values, not certified complete
air. Cramer1993 remains source-unverified: the 108-state humid grid computes
NASA-9 sound speeds but emits NO-DATA and fails until the primary comparison
and its stated uncertainty can be verified. No secondary transcription is
used to qualify the requested 5e-4 gate.

`relaxation::RelaxationAbsorption` (M0.1-W2) is the molecular-relaxation law
for any composition. Each species declares its relaxing internal reservoirs:
vibrational modes with Planck-Einstein heat capacity, and classical rotation.
They enter one exact complex `c_v(omega)`, and `k = omega/c(omega) + i alpha_SK`.
Stokes-Kirchhoff is counted once, there is no separate bulk viscosity, and the
strengths are not linearized.

For propagated histories, `RelaxationState::causal_wavenumber` instead solves
the coupled homogeneous linear continuity, viscous-momentum and energy
equations with the original complex `cv(omega)` and ideal-gas EOS. The same
caloric response governs compression and Fourier heat conduction; an
equilibrium-loss factor is no longer multiplied onto the molecular wave.
The shared thermoviscous dispersion owner selects the acoustic root with a
rationalized quadratic, retaining weak thermal loss and refusing unresolved
acoustic/thermal roots. No extra bulk viscosity, wall attenuation or thermal
mode is added. Zero reservoir strength recovers the full original classical
dispersion, including outside the leading Stokes--Kirchhoff approximation.
For complex calorics, the acoustic root continues from zero transport at the
requested frequency. An analytic positive-real discriminant minimum retains the
inexpensive separated-mode fast path. Otherwise a quadratic Bezier enclosure
bounds each complete continuation interval away from zero and tracks the
square-root sheet across its artificial principal cut. Constant storage and at
most 256 attempted intervals bound this work; coalescing, numerically unresolved
or exhausted paths refuse. A stable equivalent quadratic formula avoids a
cancelling denominator after a sheet change while the original rationalized
formula retains weak thermal loss. This is pointwise transport continuation,
not a proof of global frequency analyticity. Independent fine-step inverse-root
continuation checks the positive-reservoir Eucken crossing that previously
refused, with original polynomial and continuity/momentum/energy residuals,
exact cut endpoints, coefficient scaling, coalescence refusals and dimensional
similarity. The shared dispersion tests and `tests/relaxation_absorption.rs`
also retain the old-multiplier falsifier and full zero-reservoir limit.
Continuum, linear-response, constant-transport and original rate-data support
remain required; no experimental or finite-wavefront claim is inferred.

`RelaxationState::caloric_response` exposes that same bounded, allocation-free
reservoir sum as complex mass-specific `cv`, `cp = cv + R`, and bulk modulus
`K = p cp/cv`, including the equilibrium DC limit. Under `exp(-i omega t)`,
`Im(K) <= 0`; the cycle-mean compressional work is
`-omega Im(K) |strain|^2 / 2` per volume for peak phasors. Duct wall thermodynamics
and coupled cavity pressure use this constitutive owner directly. It adds no
shear/thermal propagation attenuation, boundary-layer loss or second bulk
viscosity. Single-reservoir real-arithmetic rational response, DC/frozen limits,
the original wavenumber relation and the zero-reservoir lossless limit are
checked in `tests/relaxation_absorption.rs`.

`RelaxationState::bulk_relaxation` realizes that same finite caloric law as
equilibrium bulk storage plus positive Maxwell compression arms. It uses the
existing symmetric Jacobi owner on the diagonal-plus-rank-one caloric matrix,
checks eigenpairs and static/frozen limits, and retains every original source
reservoir. Equal molecular times combine exactly; no fitting, pole truncation,
propagation loss or new gas-specific constitutive law is introduced. The work
is bounded by the original 32-reservoir cap. Single-reservoir analytic
elimination, multi-reservoir complex response, equal-time/zero-strength limits
and pressure/rate scaling are checked in `relaxation::bulk::tests`.

Relaxation times are data: Landau-Teller collision-partner rates with a linear
mole-fraction rule, a partner-independent power law, or ISO 9613-1's `f_rO`
and `f_rN` (`iso9613_relaxation_frequencies_hz`, bit-identical to the ISO
formula's). The ISO laws are a declared terrestrial-air variant: dry-basis O2
and N2 must be within 1 % of 0.209 and 0.781, with CO2 at most 1e-3. Any other
mixture is refused, naming the species (the same admission as
`iso9613::TerrestrialAir`). `RelaxationAbsorption::iso9613_terrestrial_air()`
is that data set: the O2 and N2 vibrations at ISO's theta and rates, with Ar,
CO2 and H2O declared without reservoirs (ISO's own assumption). Rotation is
not a reservoir in it, so the rotational bulk loss inside ISO's classical
constant is left out.
`RelaxationAbsorption::evaluate_iso9613_rates_held(medium)` evaluates the
same reservoirs for N2 and O2 in any proportion (with Ar, H2O, and CO2 at or
below ISO's ceiling) with ISO's relaxation frequencies at the medium's own T, p
and water fraction HELD: a declared assumption, since the partner rates of any
other proportion are not sourced. At ISO's own fractions it is ISO's rate
law. Any other species present refuses by name.

The law also refuses:
- a present species without data;
- a present collision partner without a rate;
- T outside a rate's validity;
- more than 32 reservoirs;
- relaxing heat capacity above translation.

G0 results (`tests/relaxation_absorption.rs`):
- zero reservoirs equal Stokes-Kirchhoff bitwise;
- the single-process `alpha lambda` peak equals its closed form (2.6e-15);
- `alpha >= 0` and `c0 <= c(omega) <= c_inf` at 10^4 random states;
- the Kramers-Kronig transform of `alpha` matches the dispersion (1.6e-9).

Against ISO 9613-1's own vibrational terms, with ISO's theta and rates:
- dry air differs by 2.0e-3 (253 K) to 1.0e-2 (323 K). This is the standard's
  small-strength, equilibrium-heat-capacity linearization.
- humid air additionally differs by about -x_H2O, the O2/N2 dilution that the
  standard omits.

The 1e-3 G1/G2 gates against ISO are therefore not met and remain open.
NO-CLAIM: no sourced H2O, CO2 or rotational rate data is bundled, so the law
refuses those species rather than defaulting them.

The universal ambient-medium primitive of the musical-acoustics
doctrine: `GasState::try_new(&GasSpec, T, p)` derives the COMPLETE
tuple (density, sound speed, dynamic viscosity, thermal conductivity,
gamma, specific gas constant, cp, Prandtl, characteristic impedance)
bottoms-up — ideal gas + Sutherland + a declared conductivity model
(USSA-1976 empirical air fit or the kinetic-theory Eucken relation).
`GasSpec::dry_air_ussa1976()` carries the public-domain U.S. Standard
Atmosphere 1976 constants and reproduces its printed sea-level values
(rho 1.2250, c 340.294, mu 1.7894e-5) to printed precision; evaluated
at 20 C it IS the pair the acoustic crates' `air()` conveniences
hardcode (asserted in the fs-couple casebook).
`GasState::stokes_kirchhoff_absorption(ω)` is the classical
viscothermal `α ~ ω²` law from `(μ, κ, γ, c_p)`.
`GasState::iso9613_absorption(air, h_r, ω)` is ISO 9613-1
(molecular + classical) with humidity as an explicit `[0, 1]`
argument; do not add Stokes on top. It needs `air`, an
`iso9613::TerrestrialAir` admission bound to this state. The admission is
made from the gas's explicit composition (`TerrestrialAir::admit`, the
relaxation law's own policy) or with the legacy moist-air state
(`TerrestrialAir::moist_air`, U.S. Standard Atmosphere 1976 dry air). Any
other composition is refused by `Iso9613Refusal`, which names the species:
a foreign species (He, SF6, ...), O2 or N2 off ISO's dry fraction, or CO2
above 1e-3. `iso9613::iso9613_terms` exposes ISO's own decomposition, its
classical coefficient and each relaxing species' strength and `f_r`.
Consumers that need ISO's dispersion (its exact Debye pairs) build it from
there, not from a second transcription. `iso9613_relative_humidity` inverts
ISO's humidity-to-water conversion. Pr is DERIVED, never an
input — the emergent ~0.71 for air is the built-in kinetic-theory
falsifier, and the Eucken-vs-USSA divergence (measured 0.1%..12% over
200..600 K, monotone-worsening) is asserted as the relation's honest
boundary. Validity window T in [50, 2000] K, p in (0, 1e7] Pa refuses
outside; the calorically-perfect gamma degrades above ~600 K
(vibrational excitation, order 1-2% there), densities near the 100 bar
ceiling carry a 1-3% real-gas error at 300 K (worse when cold),
and the USSA/Sutherland fits are EXTRAPOLATED beyond their published
bands (~187-288 K atmosphere tables; ~170-1900 K Sutherland-air) —
all documented as accuracy boundaries so extreme-regime studies
evaluate with honest error language. NO-CLAIM: phase validity is not
checked — an ideal-gas model cannot know about condensation, so a
returned state is the gas-phase extrapolation and never evidence the
gas phase exists at that (T, p); callers own the phase check until
vapor-pressure/melting taxonomy rows land. The Eucken 5/4 coefficient
is pinned ABSOLUTELY by the monatomic identity Pr = 2/3 (exact,
mu-independent), alongside the relative Eucken-vs-USSA divergence
envelope.

`resolve_sutherland_gas_state` connects a material card to this same evaluator.
It resolves five exact parameter keys at a shared, explicitly typed temperature
and pressure, retaining the complete source bundle and receipts. Reference
viscosity is converted back to Sutherland beta; it is not mistaken for viscosity
at the operating temperature. Reference temperature is absolute, Sutherland S
is an interval, and conductivity model selection is explicit. The returned
`ResolvedGasState` exposes the selected model separately; its parameter-bundle
identity does not claim to bind that model choice. Missing, ambiguous, mistyped,
out-of-source-domain or invalid parameters refuse atomically. The
`air-dry-ussa1976` seed defines a narrower 273.15–313.15 K / 80–110 kPa dry
engineering application range, with its 1976 composition and model limits.
This is source-backed model evaluation, not real-gas or humidity qualification.
The G1/G3 `g1_g3_sourced_dry_air_reaches_acoustic_loss` integration case in
`xtask/tests/matdb_pack_cli.rs` exercises real source compilation, persistent
store reopening, discovery, this adapter and the existing cylinder air-loss
law at four T/P states, with independent references, receipts, replay and
source-sensitivity/refusal controls. This evidence is numerical, not empirical.

`moist_air::resolve_moist_air_state` resolves the carrier and dilute water-vapor
cards separately at the same typed temperature and total mixture pressure,
then uses the shared moist-gas evaluator. Both parameter bundles and all ten
usage receipts remain accessible; relative humidity, vapor mole fraction and
vapor mass fraction are distinct outputs. Explicit component flags and normal
source applicability prevent silently substituting another composition.
Neither component identity alone identifies the mixture. USSA carrier
conductivity, Eucken vapor conductivity and Wilke/Wassiljewa–Mason–Saxena
mixing are explicit fixed model choices. The component evaluation does not
assert stable pure steam at the mixture's total pressure. This is bounded
source-backed engineering model evaluation, not condensation or a real-gas EOS.

### `visco` — viscoelastic damping tiers (bead ybc75)

Four tiers, one runtime path. `RayleighDamping` (tier 0: αM + βK
coefficients and the two-point fit, which must reproduce both finite damping
pins with nonnegative coefficients or refuse; materially negative coefficients
are never clipped into a different resonant response, and dimensionless
frequency ratios preserve the fit under time-unit rescaling); `FractionalZener` (fitting
canonical, 4 params, `Real`-generic so exact AD parameter tangents share
the value's code path); `GeneralizedMaxwell` (runtime Prony series with
the unconditionally stable exact-exponential recursion for
piecewise-linear strain and a per-step work/stored/dissipated energy
ledger; `matching_loss` builds one branch that hits an authored
`η(ω₀)`); `lower_to_prony` (deterministic relaxation-time ladder + NNLS,
emitting a `LoweredModel` whose certificate is the MEASURED supremum
relative modulus error over a dense in-band grid and whose evaluation
REFUSES outside the certified band, `FS-MAT-VISCO-OUT-OF-BAND`);
`ThermoelasticZener` (single-relaxation uniaxial thermoelastic loss,
  plus handbook `aluminum` / `structural_steel` constructors);
`fit_fractional_zener` (deterministic Levenberg–Marquardt over an
admissibility-preserving reparameterization, exact AD Jacobian, and a
caller-AUTHORED residual budget `tol_rel`: a stationary point above the
budget refuses with `FS-MAT-VISCO-FIT-DIVERGED` instead of returning a
model the data cannot support). The e2e evidence lane is the
`visco_casebook` binary: synthesize (deterministic FNV-keyed ±0.1%
perturbation, no RNG state) → fit → lower → single-DOF velocity-Verlet
ring-down → re-measure η from the decay envelope, one JSON line per
stage, bit-identical across reruns except `elapsed_ms`.

Scalar fractional-Zener, Maxwell and thermoelastic
spectra evaluate bounded reciprocal forms so weak physical loss and relaxation
tails are not erased by complex-product cancellation or squared-frequency overflow.
The fractional value and AD tangents share that arithmetic. These numerical
improvements do not extend any source applicability band or replace the light-loss
assumption behind modal `eta/2`.

`PronyRampSensitivity` prepares opt-in exact transient derivatives of the original
Maxwell ramp. `PronyRampLinearization` retains the original primal tuple and the
Jacobian of endpoint branch strain, interval mean and independently integrated
physical heat with respect to initial branch strain, imposed strain increment and
log relaxation time at fixed duration. Reusable preparation and tangent/pullback
application allocate nothing; forward-only consumers retain their original cost.
Cancellation-free coefficient derivatives preserve weak physical loss. Elastic
zero-duration and fully relaxed limits are explicit. Heat derivatives follow the
original positive integral, never subtraction of nearly equal work and storage.

The initial history is held fixed in each local parameter derivative. A complete
experiment must propagate it through the accepted steps; the API supplies both
push-forward and transpose pullback for that purpose. Focused tests independently
integrate the sensitivity equations, differentiate the energy identity, reverse
multibranch loading/reversal histories, and recover an authored relaxation time
from synthetic transient stresses. This is a derivative of the prescribed-strain
constitutive law, not an automatic derivative of a coupled dynamics solve, a
stress-control eigenbasis, a source selection or fractional-to-Prony fit, and it
does not establish measured-material identifiability or a new applicability band.

`visco::history_fit::fit_prony_history[_with_cancel]` fits a passive runtime
material directly to original, nonuniform strain/stress observations on a
caller-supplied relaxation-time ladder. The first observation explicitly has
fully relaxed branch memory, including at nonzero total strain. Exact original
`PronyRamp` histories form the design columns; the existing nonnegative solver
fits positive equilibrium and nonnegative branch moduli. Each observation's
authored absolute stress error weights its row. Factored column equilibration
preserves stress-unit scaling, and an independent replay through the original
material `step` must satisfy EVERY observation bound. Sample/branch capacities,
unexcited columns, duplicate times, unrepresentable arithmetic and cancellation
have explicit refusals. The small dense path supports at most 63 branches plus
equilibrium and polls before/after its bounded NNLS solve, not inside it.
G1/G3 tests recover a multibranch law through loading/reversal, reproduce heat
on a different loading history, exercise unequal observation errors and preserve
stress/clock unit changes. These are synthetic numerical checks; a successful
sampled fit does not identify unknown relaxation times, establish unique physical
parameters, validate unsampled histories or invent a source applicability band.

`visco::slab_ted` evaluates the complete insulated-slab thermal diffusion
spectrum in constant work, using a positive low-frequency series and scaled
closed form. `slab_loss_factor` is the uniaxial bending loss;
`plate_bending_relaxation_matrix` supplies the plane-stress rank-one tensor
`T (E alpha/(1-nu))^2 h^3/(12 rho c_p) [1,1,0]^T [1,1,0]`.
Thermal loss therefore acts on mean curvature and vanishes for pure twist and
trace-free saddle bending. This is first-order weak thermoelastic coupling in
a homogeneous isotropic slab with insulated faces, without lateral diffusion
or a transient heat solver. Focused tests compare the closed form with an
independent thermal-mode sum and recover uniaxial anticlastic bending.

`SlabThermoelasticBending::new_coupled` extends the explicit anisotropic
plane-stress law to coupled linear thermoelastic memory. It treats the supplied
elastic tensor as isothermal and converts the supplied constant-stress heat
capacity to `C_epsilon = C_sigma - T alpha^T Q alpha`, refusing a nonpositive
or numerically indeterminate result. The same fixed-strain capacity supplies
both the positive rank-one relaxation tensor and diffusion times. The original
`new` constructor and diagonal weak-loss consumers retain their constant-stress
capacity convention. Both models assume a homogeneous frozen state,
through-thickness Fourier conduction and insulated faces; neither supplies
lateral heat transport or evolving mean temperature.

`slab_relaxation_response` and `complex_bending_increment` expose the complete
complex thermal increment, including storage dispersion, DC and signed-frequency
conjugacy. `lower_to_prony` supplies a bounded positive odd-mode prefix to the
existing Maxwell state owners. Retained weights are never renormalized.
`SlabPronySpectrum` records the declared frequency band, analytic omitted-tail
bound, retained coefficient-rounding bound and total normalized complex-error
bound. The tail uses the convex midpoint integral bound on the infinite thermal
series; admission does not rely on a frequency grid. An exhausted term cap
returns `TruncationTolerance`. This scalar constitutive bound excludes source
tensor, FE projection, modal truncation and runtime roundoff, and is not an
arbitrary-force trajectory guarantee. It does not create material provenance.

`fit_fractional_zener_components` fits storage and loss with separate relative
residuals and AD Jacobian rows. Every measured component must meet its own
samplewise budget, so weak damping cannot disappear behind elastic stiffness.
Relative loss fitting requires positive measured loss; a zero or censored value
needs an experimental detection limit this API does not invent.
`lower_to_prony_components` applies independent storage/loss gates to the shared
positive-Prony solver and its 512-point verification grid, retaining the original
complex-error API and coefficients. At those samples the loss-factor relative
error is bounded by `(loss_rel + storage_rel)/(1 - storage_rel)` when
`storage_rel < 1`. The `visco_casebook` consumer uses both component-aware APIs
and reports both lowering errors before its ring-down. These remain numerical
sample/grid checks, without an inter-sample or experimental validation claim.

The fractional-Zener spectrum accepts signed angular frequencies with the
`exp(iωt)` convention: storage is even, loss is odd, and material-parameter
tangents obey the same conjugate symmetry needed by real Fourier histories.
At exactly zero frequency, modulus is `(E0, 0)` and the fixed-frequency
gradients with respect to `[E0, Einf, alpha, tau]` are exactly `[1, 0, 0, 0]`
and `[0, 0, 0, 0]`. This relaxed limit avoids undefined logarithmic AD products;
the DC branch does not claim a derivative with respect to frequency, whose
slope is not finite for fractional order below one.

The public `visco::prony_ramp` is the shared allocation-free scalar branch owner:
initial elastic strain plus total strain increment and `dt/tau` produce endpoint,
interval mean and independent physical dashpot heat. Consumers retain their
original spatial force assembly, branch stiffness, energy units and source
validation. `PronyRamp` prepares the same exponential and small-rate coefficients
once per branch and step; immutable `apply` calls reuse them across strain
coordinates without allocation or further transcendental evaluations.
Zero time is an elastic jump and infinite `dt/tau` the fully relaxed
limit. G0/G1 regressions cover passive-fit refusal, frequency scaling, weak-loss
AD tangents, and scalar spectral tails; they are numerical-law checks, not
experimental material validation.

`visco::tensor` (`TensorProny`, `PronyBranch`, `TensorPronyState`, `Mandel6`, IRP M0.2-W2a) extends the runtime Prony model to passive anisotropic 6×6 Mandel tensors (`[xx, yy, zz, √2 xy, √2 yz, √2 zx]`). Admission requires a symmetric positive-definite equilibrium stiffness `C_eq ≻ 0` (diagonally scaled Cholesky + deterministic Jacobi eigendecomposition via `fs_la::eigen::jacobi_eigh`), an optional positive-semidefinite Kelvin–Voigt tensor `kv_pa_s ⪰ 0`, at most 64 relaxation branches `(C_j, τ_j)` with `C_j ⪰ 0` (eigenvalues verified non-negative without clipping) and `τ_j > 0`, and a finite applicability band `band_rad_s`. Each branch is factored once at admission as `C_j = L_j L_jᵀ`. Constructors include general `TensorProny::new`, shared-eigenbasis `TensorProny::kelvin(c_eq, per_mode)`, and `TensorProny::isotropic(k, g)` (whose pure shear and pure volumetric histories reduce bit-for-bit to scalar `GeneralizedMaxwell::step`). Frame rotation `TensorProny::rotated` uses the single crate-wide Mandel Bond rotation owner `tensor::rotate_mandel` (`C → Q_6 C Q_6ᵀ`, shared with `fs-solid::linear3`). `TensorProny::plane_stress_complex` computes the frequency-dependent complex Schur complement `(Q′(ω), Q″(ω))` in any ordered `MaterialPlane` and checks passivity (`Re Q ≻ 0`, `Im Q ⪰ 0`) at `ω`; it does not claim a frequency-independent branch-wise plane-stress Prony series (lowering a plane-stress Schur target to reduced-plane Prony branches is M0.2-W2b's responsibility). Time stepping `TensorProny::step` advances `q_j ∈ ℝ⁶` componentwise via `prony_ramp` and integrates dashpot heat by polarization as a sum of squares `A_j ‖L_jᵀ (q_{j,0} + h_j δ)‖² + R_j ‖L_jᵀ δ‖² ≥ 0` with no clamp, closing `work = Δψ + dissipated` to `≤ 1e-12` relative. Intrinsic material loss only is represented here; thermoelastic, contact, support, and acoustic radiation losses are separate owners and are not folded into `TensorProny`.

The `visco_tensor` harmonic APIs use the same `TensorProny::complex_modulus`
evaluation: `modulus_mandel`, `plane_stress` and its plate sections retain the
original frequency band, including after fractional lowering. A finite input
that produces unrepresentable storage or loss refuses. G0 regressions admit
both band endpoints, refuse a one-ULP excursion, and compare the work-conjugate
plane-stress spectra through both public entry points.
Both tensor history APIs retain Kelvin stress and heat on positive-duration
ramps and refuse an instantaneous strain change when nonzero Kelvin viscosity
is present. `TensorProny::step` publishes its candidate memory and work/heat
ledger together only after finite-response checks. G0/G1 tests cover atomic
overflow refusal and retry, explicit zero viscosity, and stress/storage/heat/
work parity between Mandel and tensor-shear histories.

`TensorPronyHistory::advance_temperature_ramp` retains the original reference
law, its tensor factors and every Maxwell memory coordinate through prescribed
linear temperature and physical-time strain ramps. The existing WLF/Arrhenius
calibration scales branch times by `a_T` and Kelvin viscosity by `a_T`; all
elastic tensors remain fixed. A caller-bounded midpoint subdivision reuses the
original exact affine-strain/positive-heat transaction without rebuilding or
factoring a material law. Its logarithmic-shift variation limit is a coefficient
resolution control, not a certified trajectory-error bound. Both thermal and
strain histories require refinement.

The returned mean stress is weighted in physical time. Local signed work and
independent positive heat are accumulated before touching cumulative ledgers;
the complete history's endpoint stress uses actual endpoint-temperature Kelvin
viscosity and the prescribed physical strain rate. Constant reference temperature
retains original advance arithmetic; existing isothermal adjoints are unchanged.
This temperature-ramp extension supplies no temperature/path adjoint. Zero-time
temperature changes, out-of-calibration paths, unrepresentable shifted parameters,
exhausted work limits and cancellation refuse transactionally, including at the
final poll. No thermal expansion, temperature-dependent elastic modulus, caloric
storage, entropy or self-heating closure is inferred. Six focused regressions
cover fixed-temperature parity, loaded-memory thermal cycles, endpoint/mean
Kelvin stress, independent full-tensor ODE refinement of stress/work/heat, weak
local increments after large cumulative histories and atomic refusals.

`TensorPronyHistory::advance_linearized` retains exact matrix-free forward and
reverse derivatives of the original anisotropic Maxwell/Kelvin ramp. Both
endpoints' tensor-shear strains, every original six-coordinate branch memory and
independent log relaxation times participate; endpoint/mean stress, recoverable
storage, work and independently integrated heat are differentiated. Original
tensor factors and their directions remain fixed, including rank-deficient
channels. A shared primal transaction preserves ordinary results; immutable
tensors are shared by the linear-size tape, and Jacobian/transpose actions use
caller buffers with atomic cancellation/refusal and no allocations or
branch-count-squared matrices. Complete-history gradients require memory
propagation and accumulation of endpoint-strain/log-time cotangents. Nonzero
Kelvin viscosity requires positive duration for the full endpoint-gradient map;
pure Maxwell/elastic zero-time jumps retain their exact limit. No factor,
source, temperature or surrounding coupled-solver derivative is claimed. Six
focused checks cover tensor-shear duality and actual-history finite differences,
original-primal equality, independent heat integrals and weak/stiff limits,
differentiated energy balance, and transactional refusals.

`visco_tensor::history_fit::fit_tensor_prony_histories[_with_cancel]` jointly
fits a passive Maxwell/Kelvin law to independent multiaxial experiments. Up to
64 nonnegative strengths multiply authored Gram tensor directions for
equilibrium, fixed-time Maxwell arms and Kelvin viscosity. Signed cross terms,
equal-time directions, original order and zero strengths remain explicit.
Each experiment starts with relaxed Maxwell memory and zero strain rate at its
supplied first strain. Original linear ramps generate exact history columns;
each stress component carries its own absolute error bound. Factored column
equilibration feeds the existing NNLS owner. The equilibrium sum must pass SPD
admission, and independent original-history replay must satisfy every supplied
observation bound before publication. Seven focused actual-source checks cover
coupled recovery, independently integrated held-out heat, complementary
experiments, frame/unit/factor scaling, error bounds and transactional refusals.
This is sampled calibration on an authored tensor cone: unknown directions,
times, prehistory, uniqueness and unsampled applicability are not inferred.

`visco_tensor::stress_control::PreparedTensorStressControl` binds an original
history for exact full-tensor creep, recovery and piecewise-linear stress ramps.
Existing strain and branch memory survive switches between strain/stress
control. Original QR/SVD owners factor the positive compliance transform,
retaining weak equilibrium springs without subtracting nearly equal tensors;
rotated nontriangular Gram factors and rank-deficient branch directions remain
physical inputs. Existing Prony ramps integrate modal memory and strain rate.
Work and nonnegative heat are independent integrals, with heat weights applied
before squaring to retain extremely weak representable dissipation. Local
quadratic storage increments, original-source stress closure and conservative
whole-path strain bounds gate the atomic history update. Preparation has
explicit dimension/storage/work ceilings and bounded native factors; each ramp
uses quadratic work in original memory coordinates. Five focused actual-source
checks cover closed-form creep/recovery, an independent noncommuting tensor ODE
ladder with work/heat, original-history switching and public material rotation,
weak/stiff limits down to representable heat near 4e-304, and atomic refusals.
Nonzero Kelvin viscosity and certified forward-error bounds are outside this
stress-control API. Its prescribed temperature-ramp extension is specified below.

`visco_tensor::stress_control::PreparedTensorMixedControl` shares that same
positive compliance preparation and exact ramp transaction for any partition
of prescribed stress and strain components. The stress-controlled principal
block is eliminated from the original full tensor, retaining all branch
coordinates and noncommuting couplings. Prescribed stresses and reaction
stresses contribute independently to port work; heat remains the integral of
the original dashpot squares. Plane-stress histories prescribe the three
out-of-plane stresses and retain their evolving strain and memory, rather than
condensing each relaxation branch independently. Controllers with different
partitions can alternate on the same history; full stress and full strain are
included as limiting partitions. Original source closure, local energy balance,
whole-path strain limits and transactional refusals apply to every partition.
Five focused checks cover closed-form plane-stress creep/recovery and reaction
work, an independent noncommuting original-equation ODE refinement ladder,
all-stress/all-strain parity, partition switching and atomic refusals.
The component partition and elastic tensors stay fixed within each prepared
controller; no Kelvin or finite-strain extension is inferred.

Both controllers also expose `advance_temperature_ramp` for prescribed temperature
and stress/mixed loading that are affine in physical time. The shared
`TemperatureRamp` and `TemperatureRampControls` admit the original reference law
and a bounded midpoint-temperature partition. Each physical interval h reuses the
same prepared controller for reduced duration h/a_T; every original branch
coordinate and full-tensor constraint remains present, without refactorization,
material replacement or history resets. Mechanical work and independently
integrated dashpot heat are invariant under this time substitution and receive
no extra shift-factor multiplier. Elastic tensors stay fixed; this is the same
common kinetic-shift law used by the strain-controlled thermal owner.

`TensorStressTemperatureStep` carries the admitted endpoint shifts, actual
substep count and original response. Local strain increments, quadratic storage
changes, work and heat use compensated sums, while the reported strain bounds
include the initial state and every accepted subinterval. Each native whole-path
strain, original-source and energy admission remains binding. Positive physical
intervals whose reduced durations underflow to zero or become nonfinite refuse.
The changing-temperature transaction publishes only after all substeps and its
final cancellation poll succeed. Constant reference temperature delegates to the
original controller with unchanged response bits; zero physical time at fixed
temperature retains the existing elastic jump, and a zero-time temperature
change refuses.

The coefficient-variation bound controls midpoint kinetic resolution, not
trajectory error. Thermal/loading refinement is still required. No Kelvin
extension, thermal expansion, temperature-dependent elastic modulus, caloric
feedback or certified forward-error bound is inferred. The existing
`ramp_tangent` remains the reference-law tangent, not a derivative of this
variable-temperature path. Five focused regression sources cover isothermal
Zener creep/recovery with independent work/heat, second-order thermal creep
convergence against an independently integrated exposure, reference-temperature
and deterministic-replay parity, coupled plane-stress thermal reversals against
original physical-time equations with tensor-shear work, and atomic budget,
applicability, whole-path strain, late cancellation and reduced-time refusals.

Both tensor stress and mixed controllers expose `ramp_tangent` for exact
endpoint strain/stress derivatives with respect to the six prescribed endpoint
controls. The initial history, material and duration stay fixed. Analytical
Prony responses and the retained compliance factors supply this six-column map
without perturbed advances or numerical differentiation. Reconstructing stress
derivatives from original equilibrium and branch coordinates retains weak
equilibrium stiffness after stiff relaxation; no cancellation of full stiffness
and nearly equal relaxed corrections is required. Public matrices use physical
tensor-shear input/output conventions, including the mixed reaction rows.
Four checks cover preloaded-history finite differences, native all-strain and
analytical creep limits, independent zero-time constraint elimination, weak
equilibrium springs and cancellation/refusal. This is a load tangent for
implicit coupling, not an incoming-memory, material, time or complete-history
adjoint; a finite load still requires the primal domain and balance admission.

`visco_tensor::kelvin_control::PreparedTensorKelvinStressControl` advances
original tensor histories under stress ramps when the Kelvin viscosity is
nonzero and positive semidefinite. Its complete passive operator includes all equilibrium and
Maxwell energy coordinates; the existing QR/SVD owners prepare coupled loss
coordinates without subtracting a large relaxed-equilibrium state. Exact Prony
responses integrate strain, original memory, stress work and independent total
Kelvin/Maxwell heat. Amplitudes are weighted before unresolved primitive
integrals so representable physical increments survive weak and stiff limits.
With positive-definite viscosity, stress jumps keep strain and memory continuous
and do zero work. With directional viscosity, the exact original zero columns
identify its elastic nullspace. Native orthogonal factors eliminate its algebraic
stress constraints before preparing the same loss modes, with no positive
viscosity floor or removal of weak positive directions. The original state retains
six equilibrium coordinates plus every branch coordinate; its dynamic loss
dimension is `rank(V) + 6*branches`. A stress jump moves only the elastic null
directions and corresponding original branch spring strains, retains viscous
memory, and does reversible storage work with exactly zero heat. Endpoint strain
rate is reconstructed from the original constitutive equation, including the
physical stress slope in elastic null directions; a jump reports its subsequent
hold rate. Original source,
local work/storage/heat and whole-path strain checks precede the atomic update.
Preparation has explicit dimension/storage/work ceilings, with no dense
factorization during a ramp. Six focused checks cover Kelvin and combined
Maxwell/Kelvin closed forms, a noncommuting original-equation ODE ladder,
rotation and loading switches, extreme rates and transactional refusals.
Additional directional/shear-only checks cover exact jump and ramp work, weak
positive directions beside strong directions, noncommuting original-equation
DAE refinement, rotation, subdivision, history continuation and atomic refusals.
Zero Kelvin viscosity remains with the existing Maxwell controller. Unresolved
rank, factors or source balances refuse; this API does not infer mixed control,
derivatives or certified error bounds.

`PreparedTensorKelvinStressControl::advance_temperature_ramp` consumes the
existing calibrated temperature schedule for physical-time stress loading.
Both original Maxwell times and the complete Kelvin viscosity scale by the
same `a_T`; fixed equilibrium/branch stiffness and original memory remain
unchanged. Each midpoint-frozen subinterval reuses the reference controller for
`h/a_T`, with no factorization or fit. Stress interpolates in physical time;
local work, storage changes and independent Kelvin/Maxwell heat accumulate
without rescaling. The reported endpoint relaxation rate divides the reference
hold rate by the actual final shift, rather than the last midpoint shift. The
elastic-nullspace contribution retains the physical stress slope without scaling.
The whole history publishes only after final-rate representability and final
cancellation admission, including constant-temperature and zero-duration calls.
Constant reference temperature retains the original response bits.

The shared coefficient-variation and substep limits control this second-order
midpoint kinetic approximation; they are not trajectory-error certificates.
Elastic moduli stay fixed, temperature is prescribed, and no thermal expansion,
caloric feedback or thermal-state evolution is inferred. Focused regressions
cover shifted-law/reference parity, endpoint-rate physics, noncommuting thermal
reversals against an independent physical-coordinate ODE and refinement,
original-memory continuation into strain control, and atomic late cancellation,
applicability, budget, binding and physical-rate overflow refusals. Semidefinite
checks also cover separately shifted-law parity, the independent elastic
endpoint slope, noncommuting thermal DAE refinement and atomic jump-rate overflow.

`TensorProny::step_with_mean` exposes that same Mandel update's endpoint
stress, exact interval-mean stress, signed step work and independently integrated
step heat through `TensorPronyMandelStep`. The increments are returned before
addition to cumulative counters, retaining weak loss after long histories.
Legacy `step` delegates to this owner and retains its endpoint-stress API.
`mean_tangent_mandel(dt)` differentiates the work-mean stress with respect to
the final strain using the original stable `PronyRamp` coefficients; it is not
an endpoint-stress trapezoid. The shared solid midpoint dynamics uses this
pair to close kinetic/material energy against physical heat and applied work.

`visco::tensor::lower` (`TensorComplexModulus`, `KelvinFractional`, `ComplexEngineering`, `PlaneStressTarget`, `TensorLoweringTolerance`, `TensorCertificate`, `LoweredTensorModel`, `lower_tensor_to_prony`, IRP M0.2-W2b) lowers causal frequency-domain anisotropic complex-modulus targets `C*(ω) = C′(ω) + i C″(ω)` to passive `TensorProny` laws over `[f_lo_hz, f_hi_hz]`: 1×1 scalar `FractionalZener` targets delegate to `lower_to_prony` bit-identically; shared-eigenbasis `KelvinFractional` targets solve Lawson–Hanson NNLS (`visco::nnls`) per Kelvin mode on the shared log-spaced `τ` ladder; general non-commuting targets (`ComplexEngineering`, `PlaneStressTarget`) solve PSD-constrained weighted least squares in a simultaneously diagonalized congruence basis at the band midpoint with `jacobi_eigh` PSD cone projection inside the iteration loop (`ViscoError::FitDiverged` on cap exhaustion, no post-fit clipping). `TensorCertificate` certifies `grid_sup_storage + lipschitz_term` and `grid_sup_loss + lipschitz_term` on a 512-point log-spaced verification grid using generalized eigenvalue relative quadratic-form errors `|λ(C_fit − C, C)|` plus the inter-sample Lipschitz enclosure `L · Δln(ω) / 2`, refusing with `ViscoError::FitTolerance` when either bound exceeds `TensorLoweringTolerance`.

`ComplexEngineering` normalizes its complex normal compliance before cofactor
inversion, so singularity admission has no dimensional determinant cutoff that
would reject ordinary metal stiffnesses or change under uniform modulus scaling.
Nonfinite compliance and unrepresentable stiffness are refused; equilibrium,
storage and loss definiteness checks are retained. Regressions independently
check 100 GPa isotropic storage/loss against the Lamé identities, rotated
anisotropic spectra across modulus scales, and singular/indefinite refusals.

## Invariants

1. **Tangent consistency (the merge gate)**: every law's tangent matches
   central FD of its own stress to ≤1e−5 relative across elastic branch,
   plastic branch, cyclic states, and 9×9 hyperelastic components.
2. **Frame indifference**: isotropic small-strain σ(QεQᵀ) = Qσ(ε)Qᵀ;
   hyperelastic P(QF) = Q·P(F) — randomized rotation battery.
3. **Return-map consistency**: after every J2 update, the yield function
   at the returned stress satisfies f ≤ tolerance; dissipation increments
   σ:Δεₚ ≥ 0 (associative flow), total cycle dissipation > 0.
4. **Hysteresis fixture behavior**: M-P virgin curve approaches the b·E₀
   asymptote, tangents E₀/b·E₀ at the extremes, reverse branches soften
   below the elastic line (Bauschinger), symmetric cycles dissipate;
   Mander peaks exactly at (ε_cc, f′cc) with slope 0, softens post-peak,
   unloads to the residual strain, reloads rejoining the envelope.
5. **Calibration round-trip**: synthetic bilinear data recovers E within
   1%, H within 5%, σ_y within 2%, truth inside the fitted envelope.
6. **Inadmissible parameters refuse at construction** (ν bounds,
   compliance definiteness, Ec > Esec, b ∈ [0,1), positive yield).

## Error model

`MaterialError`: `Parameters`, `State` (e.g. det F ≤ 0), `Calibration`.
The generic `evidence_stress` adapter flags calibration-domain exit through
`Evidence.model.in_domain`. Checked constitutive paths retain their own explicit
domain refusals, including the Mander ultimate-compression bound above.

## Determinism class

**D0**: pure f64 arithmetic with fs_math::det transcendentals; no
iteration counts depend on ambient state (radial return is closed-form).

## Cancellation behavior

All updates are closed-form, allocation-free; P7 by boundedness.

## Unsafe boundary

Zero `unsafe`.

## Feature flags

None.

- Moist air (music bead 3ez8g.3.5): `GasState::try_new_moist_air`
  derives the ideal dry-air + water-vapor mixture from relative
  humidity — breath vs hall as physics, not a detune knob. Both
  mixture rules are EXACT for an ideal mixture (molar mass linear in
  `x_w`; the isochoric-heat identity `1/(γ−1) = Σ x_i/(γ_i−1)`); the
  RH → vapor-fraction step is the Buck 1996 liquid-water fit
  (`e_s = 611.21 exp((18.678 − t/234.5)t/(257.14 + t))` Pa, quoted
  −20..+50 °C, refused outside when RH > 0). Water-vapor spec
  provenance: M = 18.01528e-3 kg/mol, γ = 33.590/(33.590−8.31432)
  from NIST-JANAF H-064 cp(H2O g, 298.15 K) = 33.590 J/(mol K).
  The COMSOL/White steam reference μ = 1.12e-5 Pa s at 350 K and
  S = 1064 K gives β = 2.4186073348106757e-6. These correct the
  earlier rounded/misattributed Cp and β. The dry limit RH = 0 is
  the SAME CODE PATH — bitwise, so no golden anywhere moves.

- Cane packs (music bead 3ez8g.3.6): the reed card's material
  population, licensing-first — three CC-BY primaries in
  data/matdb/seed-v1 (reed-BLANK density + E_L with the radial cortex
  gradient named, Karczewski 2025 MDPI; reed-HEEL E_L/E_T dry AND
  water-soaked, Ukshini & Dirckx 2020 MDPI; whole-wall DMA storage
  modulus + logarithmic decrement, Weidenfeller 2018 SciELO). Blank
  and heel measure DIFFERENT zones of one plant (cortex-weighted vs
  inner) — recorded separately, never averaged; moisture conditions
  recorded per row (unstated conditions carry a
  moisture_content_known 0 flag). NO licensable shear modulus or E_R
  exists (Obataya/Taillard/Kawasaki are publisher-copyrighted) — the
  absence is asserted, not papered over.

- Tissue packs (music bead 3ez8g.3.7): the vocal-fold/lip population
  for the reduce lab and glottis islands — SIX CC-BY packs where a
  refusal was expected: human vertical stiffness gradient
  (3.2/5.5/7.0 kPa indentation with full protocol), porcine dynamic
  modulus across the phonation band (postmortem time typed per row),
  the felid shear rheometry loss-factor find (tan delta 0.33/0.15 —
  the only licensable native-tissue damping), the lone in-vivo lip
  stiffness (45 kPa WITH the hard-teeth-in-series caveat), the
  canonical two-mass MODEL CARD behind a register firewall
  (model_-prefixed names, declared not-constitutive), and the
  STATED-not-measured soft-tissue density register (1.04 g/cm3 as
  stated; the measured Perlman-Titze chain is unlicensable and no
  pack claims it). Method-scale disjointness is a first-class fact:
  Pa-scale matrix shear vs kPa-scale indentation registers are
  recorded separately, never averaged.

- Porous-absorber facility (music bead 3ez8g.5.4): `porous` module
  — JCA(L) equivalent fluid (five transport parameters + optional
  Lafarge thermal permeability; `None` degrades to Champoux–Allard
  with the degeneracy EXPLICIT: k0' = φΛ'²/8) yielding effective
  density/bulk, characteristic impedance/wavenumber, rigid-backed
  surface impedance, and normal-incidence absorption. Cards are
  licensing-first: Nguyen 2024 measured melamine characterization
  (CC-BY, uncertainties quoted in provenance) and the FOAM 02
  Basotect/Pinta identified rows (CC-BY, thick-sample fits only —
  thin-sample rows pin parameters at fit bounds and are refused as
  card material). ORACLE: the committed FOAM 02 measured mean
  spectra (machine-readable CC-BY data; no figure digitization) —
  the model reproduces measured absorption to worst 0.041 / mean
  0.018 (Basotect 60 mm) and 0.029 / 0.006 (Pinta 50 mm) inside
  authored 0.06/0.03 bands. `porous::PorousFluidLaw` (also available through
  the existing `porous::gas::PorousGasLaw` spelling) retains the
  original card, provenance, frequency support and explicitly supplied
  characterization-fluid viscosity. Intrinsic permeability
  `k0 = mu_ref / sigma_ref` stays fixed; current viscosity supplies
  `sigma = mu / k0`. `effective_fluid` and `surface_impedance_rigid_fluid`
  consume one original `ThermoviscousFluid` state: density, adiabatic bulk
  modulus, viscosity, Cp/Cv and thermal diffusion remain coherent. The shared
  JCA(L) kernel uses the fluid's actual `K_s`, including `rho*c²` for a sourced
  liquid or nonideal gas; liquid pressure is never substituted for stiffness.
  The existing gas entry preserves its `gamma*p` arithmetic in that same
  kernel. `Cp=Cv` gives exactly `K_s/porosity` without evaluating inactive
  thermal relaxation; viscous drag remains. Missing transport or Cv refuses
  through the original material owner. `PorousFluidResponse` retains the
  numerical state, dynamic density/modulus, characteristic response and
  rigid-backed surface response; original sources retain their authority.
  `fs-duct::network::NetworkLoad::PorousRigid`
  consumes the rigid-backed response as a normal-incidence boundary,
  with the exact phasor and exposed-area conversion. Source/geometry
  refusals remain explicit. No porous constitutive model is duplicated.
  The skeleton remains rigid with pores completely filled by one fluid;
  model transfer does not experimentally qualify an air-characterized
  specimen in a liquid. Original frequency support is still enforced.
  Substitution does not establish partial saturation, condensation, swelling, moving-frame waves,
  distributed lining or transient deposition of absorbed acoustic heat.
  General-fluid G0/G1 regressions cover gas parity, Darcy and high-frequency
  inertia, isothermal/adiabatic liquid stiffness, thin-layer compliance,
  exact zero thermal coupling, material-only changes and original refusals.
  Focused gas-source and duct-cylinder tests exercise data-only changes of
  temperature, pressure, gas mixture and pore permeability, including
  independent acoustic power balance through both network load paths.

## Conformance tests

`tests/active_tissue.rs` (suite `irp-s0-materials`, bead M0.3-W3): zero
activation bitwise equal to the passive law (1000 random F, stress and
tangent), the isometric Cauchy stress `T_max a f_L` against closed forms of
three force-length fixtures on a 50-point stretch grid (1e-12, orthogonal
components exactly zero), the traction-free free-contraction root at prescribed
J = 1 against an analytic balance (1e-10, three activations), the active-work
ledger `W_total = Delta Psi_passive + W_active` over a time-varying activation
(1e-10), objectivity, the finite-difference tangent of the active stress and of
the total, and the refusals above at the bound and just beyond it.

`tests/conformance.rs` (JSON verdicts, suite `fs-material/conformance`):
mt-001 FD tangent gate (every law, cyclic states); mt-002 objectivity;
mt-003 return-map consistency + dissipation; mt-004 hysteresis fixtures;
mt-005 calibration round-trip (+ degenerate refusals); mt-006 rank-one
convexity sampling for NH/MR, constructor refusals, card completeness,
Evidence domain flagging.

Visco (inline in `src/visco.rs`, 14 tests): closed forms over six
decades, α→1 classical-Zener reduction, Prony recursion vs direct
convolution at 1e-12, dissipation-ledger closure with sign-flip
mutation, ≤8-term plateau lowering with out-of-band refusal and
dropped-term mutation, modal injection vs exact quadratic roots, AD
gradients vs finite differences, bitwise-deterministic lowering,
dt→0 elastic-jump/small-dt expm1 regression, LM fit ground-truth
recovery from an offset start with negated-loss mutation refusal,
degenerate-input refusals, and the η ≥ 0 / monotone-relaxation property
sweep. `visco_casebook` is the end-to-end lane (fit noise floor,
certified lowering, ring-down η within an authored, measured gate).

`tests/prony_tensor.rs` (9 tests, suite `irp-s0` / `prony_tensor`):
PD/PSD/τ admission refusals, 64 random multiaxial energy-closure and
non-negative dissipation histories, bitwise isotropic reduction to
scalar `GeneralizedMaxwell`, Kelvin constructor vs general assembly,
SO(3) Bond rotation spectrum/dissipation invariance, Kramers–Kronig
consistency, Fourier-quadrature relaxation vs complex modulus,
plane-stress complex Schur complement passivity and elastic limit,
and out-of-band refusal.

`tests/prony_tensor_lowering.rs` (8 tests, suite `irp-s0` / `prony_tensor_lowering`):
1×1 scalar bitwise reproduction of `lower_to_prony`, KelvinFractional
shared-eigenbasis per-mode lowering within storage/loss certificates,
orthotropic `ComplexEngineering` non-commuting lowering with PSD
branches, `PlaneStressTarget` Schur-complement lowering, 65,536-point
independent dense-grid enclosure check, branch-removal mutation breaking
the certificate, term-count/band/tolerance/divergence refusals, and
bitwise determinism across repeated runs.

`src/gas.rs` moist battery — the two-source saturation rule (Buck vs
the independent Alduchov–Eskridge Magnus constants, ≤0.5% across the
window; IAPWS absolute pins 611.657 Pa at the triple point and
2339.2 Pa at 20 °C, measured 611.65/2338.3); the OIML R 111 / CIPM
conventional-density anchor (ρ = 1.2 kg/m³ DEFINED at 20 °C, 101325
Pa, 50% RH; measured 1.19885 — an independent metrology-community
check on the whole mixture chain, and humid air comes out LIGHTER);
the bitwise dry limit; mixture rules re-derived in the test from raw
literals (the anti-tautology discipline: never compare two quantities
off the same GasState); RH-monotone c/ρ with the saturated-20 °C
uplift in the published 0.8–1.8 m/s class (measured 1.263, the
Wong-1986/Cramer-1993 ~0.35% class); refusals (RH ∉ [0,1], Buck
window, the x_w ≤ 0.15 transport ceiling) and bitwise repeats.

`tests/cane_pack.rs` g2 gate — positivity + schema hygiene (every
scalar carries uncertainty + validity rows), fiber anisotropy in the
source's own band (E_L/E_T = 10 dry, wider wet), the moisture law
(water softens both axes, transverse more), cross-source RANGE OVERLAP
between the heel's static band and the whole-wall DMA range, the
blank-vs-heel DISJOINT-zones assertion with the gradient named in both
observations, the derived loss factor delta/pi in the woody class with
conversion ownership stated, and the named-absence assertion that no
shear row exists without a source.

`tests/tissue_pack.rs` g2 gate — positivity + schema hygiene, the
human vertical gradient (the source's own finding), method-scale
disjointness with named method dependence in every observation,
postmortem-time separation, loss tangents in the soft-tissue class
with tiger above lion, the lip teeth-caveat presence, the model-card
firewall (model_ prefix + declared register), and the
stated-density assertions (no measured tissue-density row may exist
without a licensable source).

`src/porous.rs` `porous_tests` (pa-001..pa-006): both measured-
spectrum oracles; analytic degeneracies (thin layer -> 0, rigid-
backed compliance limit Zs -> K/(jωd) at 7e-5, air-like card ->
free-medium Zc within the disclosed 0.07% constant mismatch);
refusals by name (parameters, out-of-band frequency, negative
thickness); the composed lined-tube demonstration (a simulated
ISO-10534-2 two-microphone rig recovers the direct absorption to
2e-16, and the melamine-lined tube's fundamental decays at
T60 = 0.13 s where the rigid tube never decays); card validation +
bitwise determinism.

## No-claim boundaries

- The porous facility claims LININGS (surface impedance/absorption
  of a rigid-framed layer), never the hammer's compression stack;
  no elastic-frame Biot waves (equivalent fluid only — frame-wave
  escalation is a named successor); no room acoustics; no oblique
  incidence in v1 (normal incidence + locally-reacting consumption).
  Pinta Plano Polar is polyester-fleece FELT-CLASS nonwoven, not
  wool piano felt — the wool-felt characterization is a recorded
  vein.


- **No intrinsic tissue calibration is inferred from the law family.**
  `tests/hyperelastic_ogden_hgo.rs` gates Ogden-1972's published three-term
  analytic nominal curves under prescribed J=1 with lateral-pressure
  elimination (uniaxial 1..7.5, equibiaxial 1..4.5, pure shear 1..5, 1e-3)
  and HGO one-/two-family incompressible Cauchy oracles (1e-6). Separate
  traction-free finite-bulk roots measure the O(1/kappa) modelling residual:
  finite kappa/mu=1e4 is not identically incompressible. The literal 1e-3
  finite-bulk Ogden endpoint band and 1e-6 HGO band at kappa/c10=1e5 are
  contradicted by the exact specified laws; tests retain those counterexamples
  and their derived residuals rather than changing the law or widening gates.
  Treloar observations retrieved/digitized by Bucchi & Hearn (2013, Appendix A)
  are diagnostic only: source uncertainty, redistribution licence and fit band
  are unstated. The M8-W1 ledger rows `ogden-1972-treloar-three-term-parameters`,
  `treloar-1944-points-bucchi-hearn-2013` (fixture-only) and
  `hgo-2000-fig14-media-parameters` record the provenance. The original Ogden
  (1972) paper was not opened: its coefficients were read from an open reproduction
  (Twizell and Ogden 1981, Table I, no unit printed). Literature coefficients are
  not a calibrated tissue card.
- **The active fibre stress is isometric, prescribed and uncalibrated.** There is
  no force-velocity factor, no fatigue or length history, and no activation
  dynamics. The force-length families in `tests/active_tissue.rs` are test
  fixtures, not muscle data; no vocal-fold or lip active curve is gated, and no
  `T_max` is sourced here (none was found for orbicularis oris). The G2 gate
  that reproduces a published force-length or prestress curve with its source's
  own parameters is not shipped yet.
- **Mander cyclic rules are the simplified elastic-unload variant**
  (declared on the card): no cyclic stiffness degradation beyond the
  residual-strain rule, no tension stiffening. Full Mander cyclic rules
  are follow-up work for the fiber-beam bead.
- **J2 cyclic response is bounded**: the opt-in linear Prager component
  translates the yield surface, but supplies no nonlinear backstress saturation,
  calibrated ratcheting, damage accumulation or fatigue-life prediction.
- **No damage/softening 3D laws yet**; failure envelopes are declarative.
- **Rank-one convexity is SAMPLED**, a necessary-condition spot test —
  not a polyconvexity proof (interval-certified convexity belongs to the
  certifier tier).
- **Calibration v0 is bilinear segmented LSQ**; nonlinear (Mander/M-P
  parameter) fitting and posterior envelopes land with fs-io CSV
  catalogs and the UQ stack.
- **Homogenization hooks**: homogenized laws register as ordinary
  `ModelCard`-carrying laws; the unit-cell pipeline itself is the
  lattice bead's scope.
- **Visco lowering certificates are measured, not analytic**: the
  `LoweredModel` bound is a dense-grid supremum on the stated band —
  no global proof, and nothing is claimed outside the band (evaluation
  there refuses). The casebook's ring-down comparison is a NUMERICAL
  consistency check of the fit/lower/step chain against its own
  declared physics; it calibrates no real material, and its η gate is
  authored from an observed run with headroom, not derived from an
  integrator error bound. The fractional→Prony path serves 3D FEM via
  fs-time's damping seam; a spatially assembled viscoelastic FEM
  operator is NOT part of this crate's claim.

