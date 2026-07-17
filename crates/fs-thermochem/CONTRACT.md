# CONTRACT: fs-thermochem

> Status: ACTIVE, CODE-FIRST SLICES. This contract covers typed,
> provenance-bound NASA-9 ideal-gas standard-state evaluation and bounded
> frozen-composition ideal-gas mixture evaluation, plus a positive-state
> mechanical ideal-gas EOS rung and a bounded Wilke dilute-gas viscosity
> mixing closure, plus bounded standard-state ideal-gas reaction equilibrium.
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
`ln(K_p)`, and `K_p` when directly representable.

Direct runtime dependencies are L1 or lower: `fs-qty`, `fs-matdb`, and
`fs-math`. The direct `fs-evidence` edge is development-only so conformance
fixtures construct the same validity type that reaches runtime transitively
inside `fs-matdb` cards. The crate has no L2-L6 dependency and owns no solver,
transport protocol, evolving state, session, persistence, or orchestration
behavior.

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
for a mass-action formulation using the retained stoichiometric exponents and
activities; it is not a rate-expression, activity, or kinetics-integrator
implementation.

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
- reacting/equilibrium composition derivatives, chemical potentials,
  activities beyond the retained ideal-gas `p_i/p0` standard-state convention,
  fugacity, departure functions, cubic/tabular real-gas or multiphase EOS
  behavior;
- a physical phase-stability verdict from the positive ideal-gas model-domain
  check, flash calculations, equilibrium-composition solves, reaction rates,
  kinetics integration, thermal conductivity, diffusion coefficients, or
  transport solves;
- chemical meaningfulness, reversibility, kinetic accessibility, detailed
  balance, or reverse-rate authority from exact bookkeeping and a computed
  standard-state `K_p` alone;
- a validated kinetic order, elementary-reaction classification, activity
  product, concentration-unit conversion, forward-law authenticity, net
  progress rate, mechanism consistency, stiffness claim, or time integration
  from the bounded reverse-progress-scale relation;
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
