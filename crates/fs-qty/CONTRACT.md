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
  AUTHORITATIVE `checked_plus/checked_minus/checked_times -> Option<Dims>`
  (bead sj31i.11: `None` on any i8 exponent overflow — saturation could
  cancel back into the valid range carrying wrong physics, so overflow
  surfaces typed at the operation). The renamed
  `saturating_plus/saturating_minus/saturating_times` survive ONLY as
  explicitly non-authoritative diagnostics (display/logging/pre-admission
  estimates that a checked chokepoint re-derives); nothing a semantic
  contract consumes may use them. The const-generic `Qty` product/quotient
  path refuses exponent overflow at COMPILE TIME (const-eval error;
  compile-fail doctest pinned).
- `QtyAny { value, dims }` — checked `try_add/try_sub` (returning
  `DimensionMismatch`), checked `try_mul/try_div/powi` (returning
  `DimensionOverflow`, which now also carries the right-hand operand for
  binary composition), `Mul/Div` operator conveniences that PANIC on
  exponent overflow (developer contract; agent-facing paths use the
  `try_*` forms), and `to_typed::<...>()` downcasts.
- `parse::parse_qty_with_budget(&str, ParseBudget) -> Result<QtyAny,
  ParseError>` — the explicitly bounded FrankenScript literal grammar;
  `parse_qty` is the compatibility entry point and always applies the public
  4,096-byte / 256-factor / 64-token-byte / 256-diagnostic-byte default.
  Grammar examples include `0.12Pa*s`, `0.5L/s`, `65deg`, `0.03m2/s3`,
  `9.81m/s^2`, `20degC`, `15%`, `2mol`, `12V`, `3Wb`, `4H`, `5Ohm`, `6S`,
  `7F`, and `8T`; evaluation is strict left-to-right `*`/`·`/`/`; prefixes are
  p n u µ m c d k M G T; whole-symbol match beats prefix+symbol (`min` is
  minutes, `T` is tesla when it is the complete symbol).
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
- `QuantitySpec` (also available as `semantic::QuantitySpec`) is the reusable
  scalar-schema identity descriptor:
  either explicit six-base dimensions with no kind claim, or one exact
  `SemanticType`. Its version-1 canonical token is exactly 12 bytes: version,
  six raw signed exponent bytes, variant, kind, two kind parameters, and value
  form. Strict decoding rejects wrong length/version/tags, nonzero reserved
  payload, and semantic dimensions that disagree with the kind; accepted
  bytes re-encode exactly. The descriptor does not attach or admit a value.
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
- `parse_qty` never panics on any input (garbage-battery-tested). Total byte
  admission precedes trimming, hashing, number scanning, or token allocation;
  factor and token caps are checked before evaluating or retaining the next
  factor/token. Every failure has an exact byte offset, kind, help, total source
  length, and UTF-8-safe bounded excerpt. Errors for byte-admitted inputs carry
  a domain-separated full-input hash; oversized inputs are deliberately not
  scanned merely to manufacture a hash after refusal.
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
- Quantity-spec encoding is injective over the current descriptor domain:
  dimension-only and semantic variants have disjoint tags, parameter slots
  are exact, unused bytes are zero, and semantic dimensions are canonicalized
  from the kind. Unknown future versions and tags refuse instead of falling
  back to dimension-only semantics.
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
`ParseError { preview, preview_start, input_bytes, source_hash, at, kind, help }`,
`JsonError { at, message }`, `SemanticError`, `QuantitySpecDecodeError`, and
`ChemistryError` are
structured values with contextual operations, kinds, axes, indices, or
arithmetic laws as applicable (P10 errors-as-guidance); no panics across the
crate boundary. `ParseError::verifies_source` checks an admitted source against
the retained full-input identity. `BudgetExceeded` identifies input bytes,
token bytes, or factor work with its limit and exact/proven-lower-bound
observation.

## Determinism class
Deterministic: pure functions of inputs; no RNG, no time, no I/O, no
platform-dependent elementary math. Semantic conversions use only fixed
constants and basic arithmetic; acoustic levels deliberately do not calculate
logarithms. JSON writing uses Rust's shortest-round-trip float formatting;
migration receipts hash the exact supplied v1 bytes and exact canonical v2
bytes with `fs-blake3`; noncanonical source spellings cannot mint receipts.
Chemical content identities use deterministic canonical encodings and
domain-separated `fs-blake3` hashes. Quantity-spec tokens use a closed,
fixed-width integer encoding with no allocation, hashing, or host byte order.

## Cancellation behavior
Scalar operations are O(1). Quantity parsing is O(admitted input bytes) with
explicit byte/token/factor bounds and bounded diagnostic retention; rejected
oversized input takes O(diagnostic bytes), not O(input bytes). Composition
conversion is O(input length), and exact conservation verification is
O(elements * species * reactions). These are metadata-boundary operations
rather than tile kernels; no Cx is required.

## Unsafe boundary
Production code has none; workspace `unsafe_code` remains denied. The isolated
`tests/parse_allocation.rs` binary uses a test-only global allocator whose
unsafe implementation delegates every allocation operation unchanged to
`System` and adds only non-allocating relaxed-atomic measurement. Its safety
invariant is documented at the implementation site and it is never linked into
the library.

## Feature flags
None. Nightly liability: `generic_const_exprs` for Mul/Div dimension
arithmetic (documented fallback: macro-generated products over the alias set;
public API unchanged).

## Conformance tests
`tests/conformance.rs`: Appendix C literal battery (qty-001), typed/erased
bit-agreement (qty-002), v1/v2 JSON and pinned crosswalk receipts (qty-003),
dimension safety (qty-004), parser totality over garbage (qty-005), exact
byte/factor/token budget boundaries and deterministic bounded diagnostics
(qty-005b), plus 1,200 deterministic fs-propcheck cases for six-component
dimension-vector laws with shrinking enabled and fixed cases retained. The
exhaustive qty-007 Cartesian
kind/form/phasor matrix emits one JSONL admission/refusal record per case, and
qty-008 proves parser-to-semantic unit-rescaling invariance. qty-009 exhausts
the quantity-spec kind/form codec, pins representative bytes, and mutation-
tests strict canonical decoding while proving value admission remains separate.
Unit tests include
the 20k-case seeded garbage battery, default/custom/zero budget boundaries,
UTF-8 excerpt and exact-offset cases, long unknown-unit/exponent refusals,
new-token collision cases, nonzero-mole round trips, strict canonical v1/v2
mutation and arity refusals, every semantic
distinction and named conversion (including exact-bit subnormal boundaries),
immutable chemistry identities, axis/order mismatches, exact elemental/charge
conservation, and checked multiply/add overflow. Compile-fail doctests prove
the type-level rejections. The separate `tests/parse_allocation.rs` binary
wraps `System` with a peak-request probe and proves a one-megabyte byte-refused
literal never causes a source-sized transient allocation.

## No-claim boundaries
- Luminous-intensity (`cd`) dimensions; candela stays out until photometry is
  real.
- Raw `Qty` aliases and `QtyAny` enforce dimensions only; aliases such as
  `Stress`/`Pressure` remain the same Rust type. Boundaries that need semantic
  distinctions must opt into `semantic::SemanticQty` rather than retag raw
  values by convention.
- `QuantitySpec` is an identity descriptor, not a value carrier, conversion
  receipt, or physics certificate. It does not retain source unit spelling,
  coordinate frame, vector component order, covariance/observation role,
  support geometry, instrument calibration, or acquisition-cost meaning;
  actual scalar values still require `SemanticQty` admission where semantics
  matter.
- Information/monetary units (refused with a pointer to fs-ir budgets).
- A full content hash for source rejected by the byte-admission gate. Such an
  error retains exact source length plus a bounded position-aware excerpt, but
  hashing the entire already-unadmitted source would violate the work bound;
  callers needing a pre-admission identity must compute it in their own
  separately budgeted ingestion layer.
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
