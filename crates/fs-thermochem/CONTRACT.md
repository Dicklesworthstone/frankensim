# CONTRACT: fs-thermochem

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

## Standard-state reaction-equilibrium operation tree

For one exactly conserved reaction column with negative reactant and positive
product coefficients `nu_i`, common ideal-gas standard-state models, and
absolute temperature `T`, version 1 evaluates active species once in canonical
`SpeciesId` order:

```text
term_i   = binary64_exact(nu_i) * g_i^0(T)
delta_g  = fixed_order_sum_i(term_i)
R_T      = R * T
ln(K_p)  = -delta_g / R_T
K_p      = deterministic_exp(ln(K_p)) when positive finite
```

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
axis, evaluates each active species once, and uses one fixed reduction and
`fs_math::det::exp`. Caller model permutation cannot change the admitted model,
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
  independently assembled fixed-order Gibbs/log result, and complete nested
  stoichiometric/conservation receipt identities;
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

Construction walks outward from the anchor and admits constant h/s offsets only when every cumulative per-region correction fits its declared bound and adjacent endpoint evaluations agree. Source gaps, nonfinite corrections, exceeded bounds, and unrepresentable endpoint agreement refuse. These bounds describe the permitted reference adjustment; they are not estimates of source-data uncertainty. Cp may retain a source-fit jump at a region boundary. Continuity admission does not establish positive Cp or Cv.

Corrected evaluation derives u = h − RT and g = h − Ts from the same h/s owner. Corrected relative enthalpy and entropy integrate the retained Cp and Cp/T curves over each crossed region, avoiding subtraction of large absolute references or independently rounded seam corrections. Absolute and relative evaluations agree within floating-point rounding. Raw models retain their original region references, source jumps, and evaluator version 1. Explicitly corrected models use evaluator version 2 and bind the anchor, bounds, original card identities, and exact offset bits in immutable policy metadata shared by the model and accepted receipts. Construction and cross-region changes are bounded by the existing 16-region admission limit; successful property-only evaluation adds no allocation.


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

`kinetics` binds up to 512 caller-authored elementary Arrhenius
reactions to the original NASA-9 bank and the exact `fs-qty`
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
global reactions, ions, or real-gas kinetics. A supplied
mechanism's physical accuracy remains a property of its original source data.

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
Rust: Region 1 liquid, Region 2 vapor, Region 4 saturation and the B23 boundary.
The source is the [IAPWS R7-97(2012) revised release](https://iapws.org/technical-guidance/release/IF97-Rev.download).
The released coefficients and equations are attributed to IAPWS; this is a
separate source from the existing SR1-86(1992) saturation formulation. They are
never silently combined into one caloric chart.

States retain typed pressure, temperature, density, specific volume, enthalpy,
internal energy, entropy and Gibbs energy, together with Cp, Cv, sound speed and
analytic thermodynamic derivatives. Admission follows the original regional
bounds through 100 MPa and 1073.15 K. Region 5 and the separate metastable-vapor formulation are absent. Region 3
is supplied by the distinct density/temperature source documented below. The source's 273.15–273.16 K extension
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

The 17 focused tests cover published verification tables, independent
thermodynamic derivative identities, density anomaly and Cv behavior, raw
saturation Gibbs mismatch, closed source boundaries, whole-interval Cp bounds,
adjacent-representable-temperature increments, independent Cp/Cp-over-T
integrals, reversal, work refusal and cancellation. The merged native
thermochem source suite passes all 156 tests. This does not claim validation of
unsupported regions or a full-workspace Cargo/DSR build.

## Maxwell--Stefan scaling across trace gradients

Fixed operation-tree evaluator version 2 forms each mass flux with its solved
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

`water::if97::regions` selects the original stable Region 1, 2 or 3 equation
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

Signed original h/u/s/v/Gibbs differences at these junctions are retained.
Section 10/Table 43 of the attributed IAPWS release describes finite regional
matching errors; a junction is not a phase transition and its caloric mismatch
is not latent heat. No interpolation, entropy/enthalpy reference correction or
common shifted equation is introduced by this router. The engineering ingress
is 273.16–1073.15 K with original regional pressure limits through 100 MPa.
Ice, Region 5, metastable continuation and the exact critical singularity remain
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

G1/G3 tests compare PT and density-coordinate states, the independent pressure
equation, residual-energy temperature derivatives and the dilute limit. The
source adapter supplies ideal-reference Cv before complete caloric and acoustic
stability can be admitted. This coordinate extension enables conservative
energy reconstruction; it does not implement a spatial compressible solver.
The cubic residual-thermodynamics convention is described in
<https://pages.nist.gov/teqp-docs/en/latest/models/cubics.html>.
