# CONTRACT: fs-qty

## Purpose and layer
Compile-time dimensional analysis (`Qty`), runtime-checked `QtyAny`, SI
unit-expression parsing, opt-in semantic quantity kinds, and exact chemical
conservation metadata — the "units" pillar of the Five Explicits (plan §11.5,
Appendix B). Layer: UTIL; its only production dependency is the Franken-only
`fs-blake3` content-identity primitive.

## Public types and semantics
- `Qty<const M: i8, const KG: i8, const S: i8, const K: i8, const A: i8,
  const MOL: i8>(f64)`
  — value in coherent SI base units; dimension carried in the type.
  Same-dimension `Add/Sub/Neg/PartialOrd` compile; mixed-dimension ones do
  not. `Mul`/`Div` perform dimension arithmetic (`generic_const_exprs`).
  Scalar `*`/`/` by `f64`. `total_cmp` provides deterministic NaN-last
  ordering. `erase()` produces a `QtyAny`.
- Aliases: `Dimensionless, Length, Area, Volume, Time, Frequency, Velocity,
  Acceleration, Mass, Density, Force, Stress/Pressure, Energy, Power,
  DynViscosity, KinViscosity, SurfaceTension, Temperature, Current,
  ElectricCharge, Voltage, MagneticFlux, Inductance, Resistance, Conductance,
  Capacitance, MagneticFluxDensity, Amount, MolarMass, AmountConcentration,
  MassFlowRate, VolumetricFlowRate, AngularVelocity, Angle`.
- `units::*` constructors (meters, millimeters, seconds, hours, kilograms,
  moles, kilograms_per_mole, moles_per_cubic_meter, kelvin, celsius, newtons,
  pascals, pascal_seconds, joules, watts, coulombs, volts, webers, henries,
  ohms, siemens, farads, teslas, hertz, liters, liters_per_second,
  meters_per_second, newtons_per_meter, radians, degrees).
- `Dims([i8; 6])` — runtime dimension vector `[m, kg, s, K, A, mol]` with
  SATURATING `plus/minus/times` (agent-facing paths must not panic on
  adversarial exponent chains; consumers reject long before saturation).
- `QtyAny { value, dims }` — checked `try_add/try_sub` (returning
  `DimensionMismatch`), checked `powi` (returning `DimensionOverflow`
  before metadata can saturate), `Mul/Div`, and `to_typed::<...>()`
  downcasts.
- `parse::parse_qty(&str) -> Result<QtyAny, ParseError>` — the FrankenScript
  literal grammar (`0.12Pa*s`, `0.5L/s`, `65deg`, `0.03m2/s3`, `9.81m/s^2`,
  `20degC`, `15%`, `2mol`, `12V`, `3Wb`, `4H`, `5Ohm`, `6S`, `7F`, `8T`);
  strict left-to-right `*`/`·`/`/`; prefixes p n u µ m c d k M G T;
  whole-symbol match beats prefix+symbol (`min` is minutes, `T` is tesla when
  it is the complete symbol).
- `json::to_json/from_json` — canonical v2
  `{"schema_version":2,"value":V,"dims":[m,kg,s,K,A,mol]}` with bit-exact
  finite-value round-trip. `decode_json` also accepts the exact historical
  implicit or explicit v1 five-vector wire, appends `mol=0`, and returns an
  immutable BLAKE3 `old_hash -> new_hash` semantic-crosswalk receipt. The
  source bytes must equal their version's canonical writer shape; whitespace,
  field-order, and numeric-spelling mutations refuse before a receipt can be
  issued. The convenience `from_json` refuses v1 so callers cannot discard
  that evidence.
- `semantic` — privately validated `SemanticQty`/`SemanticType` carriers keep
  dimensionally identical meanings distinct: absolute/delta temperature,
  mechanical/electrical angle and angular velocity, torque/energy,
  pressure/stress, tensor/engineering strain, composition basis,
  mass/amount/concentration, entropy/heat capacity, and acoustic
  pressure/power. Named operations make Celsius offsets, revolutions,
  rad/s versus rpm, pole-pair phase offsets, shear factors, molar-mass basis
  changes, and sinusoidal peak/RMS factors explicit. Whole-vector
  `Composition`, paired `PhasorQty`, and referenced `AcousticLevel` prevent
  convention fragments from drifting independently.
- `chemistry` — validated opaque `SpeciesId`, `ElementId`, and `ReactionId`;
  immutable `ElementalMatrix`, `StoichiometricMatrix`, and `ChargeVector` with
  domain-separated content identities; exact checked-`i128` verification of
  `A N = 0` and `z^T N = 0`; and an immutable `ConservationCertificate` that
  binds all three inputs. `MassAmountBasis` records whether a coherent mass/
  amount pair was derived from mass or amount and requires a positive finite
  molar mass.

## Invariants
- All stored values are coherent SI base units; unit conversion happens ONLY
  at parse/constructor boundaries.
- Typed and erased algebra agree bit-for-bit (same f64 operations).
- `parse_qty` never panics on any input (garbage-battery-tested); every
  failure is a structured `ParseError` with position, kind, and help.
- Angles are dimensionless radians; `deg` converts numerically. `degC` is
  affine and legal only as a lone unit; compounds are rejected with guidance.
- Accumulated unit exponents beyond ±60 are rejected as unphysical.
- Semantic carriers validate dimensions, finiteness, amplitude/range domains,
  and exact source/target kinds before any named conversion. An exhaustive,
  closed kind/form matrix admits instantaneous, peak, RMS, and paired phasor
  forms only for temperature differences, angles/angular velocities, torque,
  pressure, stress/strain, and acoustic pressure. Absolute temperature,
  energy, composition, mass/amount/molar mass/concentrations,
  entropy/heat-capacity, and acoustic power are static-only and fail with a
  typed form-policy error at carrier construction. Pole-pair counts are
  positive and their electrical phase offsets are finite; angle maps apply the
  offset while angular-velocity maps do not. Offset-bearing operations accept
  only static/instantaneous point values, while linear waveform amplitudes may
  cross domains without applying an offset.
- Positive scalar mass/amount and concentration-basis conversions never retain
  a rounded zero. If the exact positive result is below the representable f64
  domain, the conversion returns a typed representability refusal; a true zero
  source remains a legal exact zero. Non-finite overflow remains a structured
  finite-value refusal.
- Composition values are immutable, nonempty, finite unit fractions whose
  deterministic input-order sum is one within the fixed documented tolerance;
  mass/mole conversion is whole-vector and requires one positive finite molar
  mass per component. Conversion rescales by the minimum/maximum active molar
  mass before normalization so finite extreme inputs do not overflow or erase
  every active weight.
- Chemical identifiers, axis order, matrix values, and charge values are
  immutable after validation. Conservation certificates are issued only after
  exact checked arithmetic proves both matrix equalities and bind the exact
  content identities presented to the verifier.

## Error model
`DimensionMismatch { op, left, right }`, `DimensionOverflow { op, dims, factor }`,
`ParseError { input, at, kind, help }`, `JsonError { at, message }`,
`SemanticError`, and `ChemistryError` are structured values with contextual
operations, kinds, axes, indices, or arithmetic laws as applicable (P10
errors-as-guidance); no panics across the crate boundary.

## Determinism class
Deterministic: pure functions of inputs; no RNG, no time, no I/O, no
platform-dependent elementary math. Semantic conversions use only fixed
constants and basic arithmetic; acoustic levels deliberately do not calculate
logarithms. JSON writing uses Rust's shortest-round-trip float formatting;
migration receipts hash the exact supplied v1 bytes and exact canonical v2
bytes with `fs-blake3`; noncanonical source spellings cannot mint receipts.
Chemical content identities use deterministic canonical encodings and
domain-separated `fs-blake3` hashes.

## Cancellation behavior
Scalar operations are O(1), parsing and composition conversion are O(input
length), and exact conservation verification is O(elements * species *
reactions). These are metadata-boundary operations rather than tile kernels;
no Cx is required.

## Unsafe boundary
None. `unsafe_code` denied.

## Feature flags
None. Nightly liability: `generic_const_exprs` for Mul/Div dimension
arithmetic (documented fallback: macro-generated products over the alias set;
public API unchanged).

## Conformance tests
`tests/conformance.rs`: Appendix C literal battery (qty-001), typed/erased
bit-agreement (qty-002), v1/v2 JSON and pinned crosswalk receipts (qty-003),
dimension safety (qty-004), parser totality over garbage (qty-005), plus 1,200
deterministic fs-propcheck cases for six-component dimension-vector laws with
shrinking enabled and fixed cases retained. The exhaustive qty-007 Cartesian
kind/form/phasor matrix emits one JSONL admission/refusal record per case, and
qty-008 proves parser-to-semantic unit-rescaling invariance. Unit tests include
the 20k-case seeded garbage battery, new-token collision cases, nonzero-mole
round trips, strict canonical v1/v2 mutation and arity refusals, every semantic
distinction and named conversion (including exact-bit subnormal boundaries),
immutable chemistry identities, axis/order mismatches, exact elemental/charge
conservation, and checked multiply/add overflow. Compile-fail doctests prove
the type-level rejections.

## No-claim boundaries
- Luminous-intensity (`cd`) dimensions; candela stays out until photometry is
  real.
- Raw `Qty` aliases and `QtyAny` enforce dimensions only; aliases such as
  `Stress`/`Pressure` remain the same Rust type. Boundaries that need semantic
  distinctions must opt into `semantic::SemanticQty` rather than retag raw
  values by convention.
- Information/monetary units (refused with a pointer to fs-ir budgets).
- Dimensioned roots (sqrt only on `Dimensionless`).
- Unit RECONSTRUCTION in display (`kg·m^-1·s^-2` exponent form only — no
  derived-unit naming like "Pa"); format→parse round-trip is guaranteed for
  dimensionless only.
- General complex arithmetic or signal processing; `PhasorQty` only binds a
  real/imaginary pair to one waveform-capable kind and peak/RMS convention.
- Logarithmic acoustic conversion; `AcousticLevel` validates and retains an
  explicit positive physical reference but delegates logarithms to the owning
  acoustics/math layer.
- Periodic-table lookup, chemical-formula parsing, chemistry/kinetics/
  thermodynamic validity, or reconciliation of opaque species labels with
  caller-supplied elemental/charge tables. A conservation certificate proves
  only the exact `A N = 0` and `z^T N = 0` bookkeeping laws for its bound
  artifacts; zero reaction columns remain legal and carry no meaning claim.
