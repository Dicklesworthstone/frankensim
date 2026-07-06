# CONTRACT: fs-qty

## Purpose and layer
Compile-time dimensional analysis (`Qty`), runtime-checked `QtyAny`, and SI
unit-expression parsing — the "units" pillar of the Five Explicits (plan
§11.5, Appendix B). Layer: UTIL (dependency-free; usable by every layer).

## Public types and semantics
- `Qty<const M: i8, const KG: i8, const S: i8, const K: i8, const A: i8>(f64)`
  — value in coherent SI base units; dimension carried in the type.
  Same-dimension `Add/Sub/Neg/PartialOrd` compile; mixed-dimension ones do
  not. `Mul`/`Div` perform dimension arithmetic (`generic_const_exprs`).
  Scalar `*`/`/` by `f64`. `total_cmp` provides deterministic NaN-last
  ordering. `erase()` produces a `QtyAny`.
- Aliases: `Dimensionless, Length, Area, Volume, Time, Frequency, Velocity,
  Acceleration, Mass, Density, Force, Stress/Pressure, Energy, Power,
  DynViscosity, KinViscosity, SurfaceTension, Temperature, Current,
  MassFlowRate, VolumetricFlowRate, AngularVelocity, Angle`.
- `units::*` constructors (meters, millimeters, seconds, hours, kilograms,
  kelvin, celsius, newtons, pascals, pascal_seconds, joules, watts, hertz,
  liters, liters_per_second, meters_per_second, newtons_per_meter, radians,
  degrees).
- `Dims([i8; 5])` — runtime dimension vector `[m, kg, s, K, A]` with
  SATURATING `plus/minus/times` (agent-facing paths must not panic on
  adversarial exponent chains; consumers reject long before saturation).
- `QtyAny { value, dims }` — checked `try_add/try_sub` (returning
  `DimensionMismatch`), unchecked-but-correct `Mul/Div/powi`, and
  `to_typed::<...>()` downcasts.
- `parse::parse_qty(&str) -> Result<QtyAny, ParseError>` — the FrankenScript
  literal grammar (`0.12Pa*s`, `0.5L/s`, `65deg`, `0.03m2/s3`, `9.81m/s^2`,
  `20degC`, `15%`); strict left-to-right `*`/`·`/`/`; prefixes p n u µ m c d
  k M G T; whole-symbol match beats prefix+symbol (`min` is minutes).
- `json::to_json/from_json` — canonical `{"value":V,"dims":[m,kg,s,K,A]}`;
  bit-exact round-trip for finite values; strict parser (the writer is ours,
  deviation = corruption).

## Invariants
- All stored values are coherent SI base units; unit conversion happens ONLY
  at parse/constructor boundaries.
- Typed and erased algebra agree bit-for-bit (same f64 operations).
- `parse_qty` never panics on any input (garbage-battery-tested); every
  failure is a structured `ParseError` with position, kind, and help.
- Angles are dimensionless radians; `deg` converts numerically. `degC` is
  affine and legal only as a lone unit; compounds are rejected with guidance.
- Accumulated unit exponents beyond ±60 are rejected as unphysical.

## Error model
`DimensionMismatch { op, left, right }`, `ParseError { input, at, kind, help }`,
`JsonError { at, message }` — all structured values with teaching messages
(P10 errors-as-guidance); no panics across the crate boundary.

## Determinism class
Deterministic: pure functions of inputs; no RNG, no time, no I/O, no
platform-dependent math (arithmetic only — elementary functions live in
fs-math). JSON writing uses Rust's shortest-round-trip float formatting.

## Cancellation behavior
No long-running paths; all operations are O(input length). No Cx required.

## Unsafe boundary
None. `unsafe_code` denied.

## Feature flags
None. Nightly liability: `generic_const_exprs` for Mul/Div dimension
arithmetic (documented fallback: macro-generated products over the alias set;
public API unchanged).

## Conformance tests
`tests/conformance.rs`: Appendix C literal battery (qty-001), typed/erased
bit-agreement (qty-002), JSON round-trip (qty-003), dimension safety
(qty-004), parser totality over garbage (qty-005). Unit tests: 29 cases
including the 20k-case seeded garbage battery and boundary battery.
Compile-fail doctests prove the type-level rejections.

## No-claim boundaries
- `mol`/`cd` dimensions (outside the 5-vector).
- Information/monetary units (refused with a pointer to fs-ir budgets).
- Dimensioned roots (sqrt only on `Dimensionless`).
- Unit RECONSTRUCTION in display (`kg·m^-1·s^-2` exponent form only — no
  derived-unit naming like "Pa"); format→parse round-trip is guaranteed for
  dimensionless only.
- Affine temperature arithmetic beyond parse-time conversion (differences of
  Celsius are the caller's responsibility to express in kelvin).
