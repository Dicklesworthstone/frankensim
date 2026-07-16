# CONTRACT: fs-thermochem

> Status: ACTIVE, CODE-FIRST SLICES. This contract covers typed,
> provenance-bound NASA-9 ideal-gas standard-state evaluation and bounded
> frozen-composition ideal-gas mixture evaluation, plus a positive-state
> mechanical ideal-gas EOS rung. Central batch compilation and test execution
> are pending; the parent thermochemistry bead remains in progress.

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
molar mass to density, molar volume, specific gas constant, and `Z = 1`.

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
there is no useful `Cx` tile boundary. Future equilibrium, kinetics, or
database operations require explicit work budgets and cancellation/drain
semantics before landing.

Mechanical ideal-gas construction and evaluation are fixed-size scalar work
with no useful cancellation tile. They neither consume a `Cx` budget nor own
request-drain-finalize behavior.

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

These tests are code-first and batch-verification pending. A sourced external
NASA/Cantera numerical oracle battery, adversarial source-card mutation battery,
and retained cross-ISA evidence are required follow-ups; the synthetic
algebraic fixtures do not substitute for them.

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
  activities beyond the aggregate ideal-mixture law, fugacity, departure
  functions, cubic/tabular real-gas or multiphase EOS behavior;
- a physical phase-stability verdict from the positive ideal-gas model-domain
  check, flash calculations, chemical equilibrium, reaction rates,
  kinetics integration, transport coefficients, or transport solves;
- uncertainty propagation, interval enclosure, rounding-error certification,
  or validity outside the exact admitted temperature regions;
- cross-ISA bit stability before retained G5 evidence exists;
- any L3 gas-state protocol, time evolution, L6 session state, ledger storage,
  planner authority, or admission policy.

Those capabilities require separate typed laws, contracts, tests, budgets, and
proof-bearing beads. These initial evaluators must not be used as evidence that
the broader `fs-thermochem` roadmap is complete.
