//! fs-qty — compile-time dimensional analysis (`Qty`), runtime-checked
//! [`QtyAny`], and SI unit-expression parsing.
//!
//! This crate is the "units" pillar of the Five Explicits (plan §11.5,
//! Appendix B): a pressure cannot be added to a time IN THE TYPE SYSTEM,
//! and runtime-loaded data carries its dimensions as checked values.
//!
//! Dimension vector: `(M, KG, S, K, A, MOL)` — metre, kilogram, second,
//! kelvin, ampere, and amount-of-substance exponents as `i8`. Angles are
//! dimensionless (radians); `deg` parses with numeric conversion. Luminous
//! intensity (cd) remains out of scope until photometry is real, as do
//! information/monetary units (`GiB` budgets belong to fs-ir).
//!
//! Nightly note: multiplication/division dimension arithmetic uses
//! `generic_const_exprs` (a documented nightly liability, see CONTRACT.md;
//! addition/subtraction/comparison are stable-compatible). If that feature
//! ever regresses, the fallback is macro-generated products over the alias
//! set — the public API is designed so that swap is not a breaking change.

#![feature(generic_const_exprs)]
#![allow(incomplete_features)]

pub mod chemistry;
pub mod json;
pub mod parse;
pub mod semantic;

use core::cmp::Ordering;
use core::fmt;
use core::ops::{Add, Div, Mul, Neg, Sub};

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Number of admitted SI base dimensions: `[m, kg, s, K, A, mol]`.
pub const DIMENSION_COUNT: usize = 6;

// ---------------------------------------------------------------------------
// Dims: the runtime dimension vector shared by Qty (const) and QtyAny (value).
// ---------------------------------------------------------------------------

/// A dimension vector `[m, kg, s, K, A, mol]` of SI base-unit exponents.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct Dims(pub [i8; DIMENSION_COUNT]);

impl Dims {
    /// The dimensionless vector.
    pub const NONE: Dims = Dims([0; DIMENSION_COUNT]);

    /// AUTHORITATIVE component-wise sum (dimension of a product):
    /// `None` when any exponent leaves `i8` (bead sj31i.11). Saturation
    /// is NOT admissible for anything a semantic contract consumes —
    /// a saturate-then-cancel chain re-enters the valid exponent range
    /// carrying wrong physics (e.g. `(100 + 50) − 30` saturates to `97`
    /// instead of the true `120`), so overflow must surface typed at
    /// the operation, not be clamped.
    #[must_use]
    pub const fn checked_plus(self, other: Dims) -> Option<Dims> {
        let a = self.0;
        let b = other.0;
        let mut out = [0i8; DIMENSION_COUNT];
        let mut i = 0;
        while i < DIMENSION_COUNT {
            match a[i].checked_add(b[i]) {
                Some(e) => out[i] = e,
                None => return None,
            }
            i += 1;
        }
        Some(Dims(out))
    }

    /// AUTHORITATIVE component-wise difference (dimension of a
    /// quotient); `None` on `i8` overflow. See [`Dims::checked_plus`].
    #[must_use]
    pub const fn checked_minus(self, other: Dims) -> Option<Dims> {
        let a = self.0;
        let b = other.0;
        let mut out = [0i8; DIMENSION_COUNT];
        let mut i = 0;
        while i < DIMENSION_COUNT {
            match a[i].checked_sub(b[i]) {
                Some(e) => out[i] = e,
                None => return None,
            }
            i += 1;
        }
        Some(Dims(out))
    }

    /// AUTHORITATIVE exponent scaling (dimension of an integer power);
    /// `None` on `i8` overflow (including `-128 * -1`). See
    /// [`Dims::checked_plus`].
    #[must_use]
    pub const fn checked_times(self, n: i8) -> Option<Dims> {
        let a = self.0;
        let mut out = [0i8; DIMENSION_COUNT];
        let mut i = 0;
        while i < DIMENSION_COUNT {
            match a[i].checked_mul(n) {
                Some(e) => out[i] = e,
                None => return None,
            }
            i += 1;
        }
        Some(Dims(out))
    }

    /// Non-authoritative DIAGNOSTIC sum, saturating at the `i8` bounds
    /// (bead sj31i.11: the renamed former `plus`). Only for display,
    /// logging, or pre-admission estimates that a checked chokepoint
    /// re-derives — a saturated vector can alias valid physics after
    /// later cancellation and MUST NOT feed a semantic contract.
    #[must_use]
    pub const fn saturating_plus(self, other: Dims) -> Dims {
        let a = self.0;
        let b = other.0;
        Dims([
            a[0].saturating_add(b[0]),
            a[1].saturating_add(b[1]),
            a[2].saturating_add(b[2]),
            a[3].saturating_add(b[3]),
            a[4].saturating_add(b[4]),
            a[5].saturating_add(b[5]),
        ])
    }

    /// Non-authoritative DIAGNOSTIC difference; see
    /// [`Dims::saturating_plus`].
    #[must_use]
    pub const fn saturating_minus(self, other: Dims) -> Dims {
        let a = self.0;
        let b = other.0;
        Dims([
            a[0].saturating_sub(b[0]),
            a[1].saturating_sub(b[1]),
            a[2].saturating_sub(b[2]),
            a[3].saturating_sub(b[3]),
            a[4].saturating_sub(b[4]),
            a[5].saturating_sub(b[5]),
        ])
    }

    /// Non-authoritative DIAGNOSTIC scaling; see
    /// [`Dims::saturating_plus`].
    #[must_use]
    pub const fn saturating_times(self, n: i8) -> Dims {
        let a = self.0;
        Dims([
            a[0].saturating_mul(n),
            a[1].saturating_mul(n),
            a[2].saturating_mul(n),
            a[3].saturating_mul(n),
            a[4].saturating_mul(n),
            a[5].saturating_mul(n),
        ])
    }

    /// True if all exponents are zero.
    #[must_use]
    pub const fn is_none(self) -> bool {
        matches!(self.0, [0, 0, 0, 0, 0, 0])
    }

    /// Canonical unit string, e.g. `kg·m^-1·s^-2`; `1` for dimensionless.
    /// Order follows SI custom for mechanics: kg, m, s, K, A, mol.
    #[must_use]
    pub fn unit_string(self) -> String {
        let [m, kg, s, k, a, mol] = self.0;
        let mut parts: Vec<String> = Vec::new();
        for (sym, e) in [
            ("kg", kg),
            ("m", m),
            ("s", s),
            ("K", k),
            ("A", a),
            ("mol", mol),
        ] {
            match e {
                0 => {}
                1 => parts.push(sym.to_string()),
                e => parts.push(format!("{sym}^{e}")),
            }
        }
        if parts.is_empty() {
            "1".to_string()
        } else {
            parts.join("·")
        }
    }
}

impl fmt::Debug for Dims {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Dims({})", self.unit_string())
    }
}

// ---------------------------------------------------------------------------
// Qty: compile-time dimensioned scalar.
// ---------------------------------------------------------------------------

/// A dimensioned `f64`: `Qty<M, KG, S, K, A, MOL>` carries SI base-unit exponents
/// in its type. Same-dimension addition/subtraction/comparison compile;
/// mixed-dimension ones do not (see the compile-fail doctests below).
///
/// ```
/// use fs_qty::{Length, Time, Velocity};
/// let d = Length::new(6.0);
/// let t = Time::new(2.0);
/// let v: Velocity = d / t;
/// assert!((v.value() - 3.0).abs() < 1e-12);
/// ```
///
/// Adding quantities of different dimensions is a compile error:
///
/// ```compile_fail
/// use fs_qty::{Length, Time};
/// let _ = Length::new(1.0) + Time::new(1.0); // ERROR: mismatched types
/// ```
///
/// So is assigning a product to the wrong dimension:
///
/// ```compile_fail
/// use fs_qty::{Length, Volume};
/// let _: Volume = Length::new(2.0) * Length::new(3.0); // ERROR: Area, not Volume
/// ```
///
/// Exponent overflow on the const-generic path refuses at COMPILE TIME
/// (bead sj31i.11 — const/runtime parity with [`Dims::checked_plus`]):
///
/// ```compile_fail
/// use fs_qty::Qty;
/// let ceiling: Qty<127, 0, 0, 0, 0, 0> = Qty(1.0);
/// let meter: Qty<1, 0, 0, 0, 0, 0> = Qty(1.0);
/// let _ = ceiling * meter; // ERROR: const-eval `127 + 1` overflows i8
/// ```
#[derive(Clone, Copy, PartialEq, PartialOrd, Default)]
pub struct Qty<const M: i8, const KG: i8, const S: i8, const K: i8, const A: i8, const MOL: i8>(
    pub f64,
);

impl<const M: i8, const KG: i8, const S: i8, const K: i8, const A: i8, const MOL: i8>
    Qty<M, KG, S, K, A, MOL>
{
    /// This type's dimension vector as a value.
    pub const DIMS: Dims = Dims([M, KG, S, K, A, MOL]);

    /// Wrap a raw value already expressed in coherent SI base units.
    #[must_use]
    pub const fn new(value: f64) -> Self {
        Qty(value)
    }

    /// The raw value in coherent SI base units.
    #[must_use]
    pub const fn value(self) -> f64 {
        self.0
    }

    /// Absolute value.
    #[must_use]
    pub fn abs(self) -> Self {
        Qty(self.0.abs())
    }

    /// Component-wise minimum.
    #[must_use]
    pub fn min(self, other: Self) -> Self {
        Qty(self.0.min(other.0))
    }

    /// Component-wise maximum.
    #[must_use]
    pub fn max(self, other: Self) -> Self {
        Qty(self.0.max(other.0))
    }

    /// Total order for sorting (NaN sorts last, deterministically — the
    /// deterministic tie-breaking discipline reaches down to comparators).
    #[must_use]
    pub fn total_cmp(&self, other: &Self) -> Ordering {
        self.0.total_cmp(&other.0)
    }

    /// Erase the compile-time dimensions into a runtime-checked [`QtyAny`].
    #[must_use]
    pub const fn erase(self) -> QtyAny {
        QtyAny {
            value: self.0,
            dims: Self::DIMS,
        }
    }
}

impl<const M: i8, const KG: i8, const S: i8, const K: i8, const A: i8, const MOL: i8> fmt::Debug
    for Qty<M, KG, S, K, A, MOL>
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} {}", self.0, Self::DIMS.unit_string())
    }
}

impl<const M: i8, const KG: i8, const S: i8, const K: i8, const A: i8, const MOL: i8> fmt::Display
    for Qty<M, KG, S, K, A, MOL>
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Debug::fmt(self, f)
    }
}

impl<const M: i8, const KG: i8, const S: i8, const K: i8, const A: i8, const MOL: i8> Add
    for Qty<M, KG, S, K, A, MOL>
{
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Qty(self.0 + rhs.0)
    }
}

impl<const M: i8, const KG: i8, const S: i8, const K: i8, const A: i8, const MOL: i8> Sub
    for Qty<M, KG, S, K, A, MOL>
{
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Qty(self.0 - rhs.0)
    }
}

impl<const M: i8, const KG: i8, const S: i8, const K: i8, const A: i8, const MOL: i8> Neg
    for Qty<M, KG, S, K, A, MOL>
{
    type Output = Self;
    fn neg(self) -> Self {
        Qty(-self.0)
    }
}

impl<const M: i8, const KG: i8, const S: i8, const K: i8, const A: i8, const MOL: i8> Mul<f64>
    for Qty<M, KG, S, K, A, MOL>
{
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        Qty(self.0 * rhs)
    }
}

impl<const M: i8, const KG: i8, const S: i8, const K: i8, const A: i8, const MOL: i8>
    Mul<Qty<M, KG, S, K, A, MOL>> for f64
{
    type Output = Qty<M, KG, S, K, A, MOL>;
    fn mul(self, rhs: Qty<M, KG, S, K, A, MOL>) -> Qty<M, KG, S, K, A, MOL> {
        Qty(self * rhs.0)
    }
}

impl<const M: i8, const KG: i8, const S: i8, const K: i8, const A: i8, const MOL: i8> Div<f64>
    for Qty<M, KG, S, K, A, MOL>
{
    type Output = Self;
    fn div(self, rhs: f64) -> Self {
        Qty(self.0 / rhs)
    }
}

/// Dimension-arithmetic multiplication: `Qty<a> * Qty<b> = Qty<a+b>`.
impl<
    const M1: i8,
    const KG1: i8,
    const S1: i8,
    const K1: i8,
    const A1: i8,
    const MOL1: i8,
    const M2: i8,
    const KG2: i8,
    const S2: i8,
    const K2: i8,
    const A2: i8,
    const MOL2: i8,
> Mul<Qty<M2, KG2, S2, K2, A2, MOL2>> for Qty<M1, KG1, S1, K1, A1, MOL1>
where
    Qty<{ M1 + M2 }, { KG1 + KG2 }, { S1 + S2 }, { K1 + K2 }, { A1 + A2 }, { MOL1 + MOL2 }>: Sized,
{
    type Output =
        Qty<{ M1 + M2 }, { KG1 + KG2 }, { S1 + S2 }, { K1 + K2 }, { A1 + A2 }, { MOL1 + MOL2 }>;
    fn mul(self, rhs: Qty<M2, KG2, S2, K2, A2, MOL2>) -> Self::Output {
        Qty(self.0 * rhs.0)
    }
}

/// Dimension-arithmetic division: `Qty<a> / Qty<b> = Qty<a-b>`.
impl<
    const M1: i8,
    const KG1: i8,
    const S1: i8,
    const K1: i8,
    const A1: i8,
    const MOL1: i8,
    const M2: i8,
    const KG2: i8,
    const S2: i8,
    const K2: i8,
    const A2: i8,
    const MOL2: i8,
> Div<Qty<M2, KG2, S2, K2, A2, MOL2>> for Qty<M1, KG1, S1, K1, A1, MOL1>
where
    Qty<{ M1 - M2 }, { KG1 - KG2 }, { S1 - S2 }, { K1 - K2 }, { A1 - A2 }, { MOL1 - MOL2 }>: Sized,
{
    type Output =
        Qty<{ M1 - M2 }, { KG1 - KG2 }, { S1 - S2 }, { K1 - K2 }, { A1 - A2 }, { MOL1 - MOL2 }>;
    fn div(self, rhs: Qty<M2, KG2, S2, K2, A2, MOL2>) -> Self::Output {
        Qty(self.0 / rhs.0)
    }
}

impl Dimensionless {
    /// Square root (dimensionless only; dimensioned roots are a no-claim
    /// until a certified even-exponent path exists).
    #[must_use]
    pub fn sqrt(self) -> Dimensionless {
        Qty(self.0.sqrt())
    }
}

impl From<f64> for Dimensionless {
    fn from(v: f64) -> Self {
        Qty(v)
    }
}

impl From<Dimensionless> for f64 {
    fn from(q: Dimensionless) -> f64 {
        q.0
    }
}

// ---------------------------------------------------------------------------
// The working alias set (plan §11.5 / Appendix B).
// ---------------------------------------------------------------------------

/// Pure number.
pub type Dimensionless = Qty<0, 0, 0, 0, 0, 0>;
/// Metres.
pub type Length = Qty<1, 0, 0, 0, 0, 0>;
/// Square metres.
pub type Area = Qty<2, 0, 0, 0, 0, 0>;
/// Cubic metres.
pub type Volume = Qty<3, 0, 0, 0, 0, 0>;
/// Seconds.
pub type Time = Qty<0, 0, 1, 0, 0, 0>;
/// Hertz.
pub type Frequency = Qty<0, 0, -1, 0, 0, 0>;
/// Metres per second.
pub type Velocity = Qty<1, 0, -1, 0, 0, 0>;
/// Metres per second squared.
pub type Acceleration = Qty<1, 0, -2, 0, 0, 0>;
/// Kilograms.
pub type Mass = Qty<0, 1, 0, 0, 0, 0>;
/// Kilograms per cubic metre.
pub type Density = Qty<-3, 1, 0, 0, 0, 0>;
/// Newtons (kg·m·s⁻²).
pub type Force = Qty<1, 1, -2, 0, 0, 0>;
/// Pascals (kg·m⁻¹·s⁻²). Stress and pressure share a dimension.
pub type Stress = Qty<-1, 1, -2, 0, 0, 0>;
/// Alias of [`Stress`].
pub type Pressure = Stress;
/// Joules (kg·m²·s⁻²).
pub type Energy = Qty<2, 1, -2, 0, 0, 0>;
/// Watts (kg·m²·s⁻³).
pub type Power = Qty<2, 1, -3, 0, 0, 0>;
/// Dynamic viscosity, Pa·s (kg·m⁻¹·s⁻¹).
pub type DynViscosity = Qty<-1, 1, -1, 0, 0, 0>;
/// Kinematic viscosity, m²/s.
pub type KinViscosity = Qty<2, 0, -1, 0, 0, 0>;
/// Surface tension, N/m (kg·s⁻²).
pub type SurfaceTension = Qty<0, 1, -2, 0, 0, 0>;
/// Kelvin.
pub type Temperature = Qty<0, 0, 0, 1, 0, 0>;
/// Amperes.
pub type Current = Qty<0, 0, 0, 0, 1, 0>;
/// Coulombs (A·s).
pub type ElectricCharge = Qty<0, 0, 1, 0, 1, 0>;
/// Volts (kg·m²·s⁻³·A⁻¹).
pub type Voltage = Qty<2, 1, -3, 0, -1, 0>;
/// Webers (kg·m²·s⁻²·A⁻¹).
pub type MagneticFlux = Qty<2, 1, -2, 0, -1, 0>;
/// Henries (kg·m²·s⁻²·A⁻²).
pub type Inductance = Qty<2, 1, -2, 0, -2, 0>;
/// Ohms (kg·m²·s⁻³·A⁻²).
pub type Resistance = Qty<2, 1, -3, 0, -2, 0>;
/// Siemens (kg⁻¹·m⁻²·s³·A²).
pub type Conductance = Qty<-2, -1, 3, 0, 2, 0>;
/// Farads (kg⁻¹·m⁻²·s⁴·A²).
pub type Capacitance = Qty<-2, -1, 4, 0, 2, 0>;
/// Teslas (kg·s⁻²·A⁻¹).
pub type MagneticFluxDensity = Qty<0, 1, -2, 0, -1, 0>;
/// Moles.
pub type Amount = Qty<0, 0, 0, 0, 0, 1>;
/// Kilograms per mole.
pub type MolarMass = Qty<0, 1, 0, 0, 0, -1>;
/// Moles per cubic metre.
pub type AmountConcentration = Qty<-3, 0, 0, 0, 0, 1>;
/// Kilograms per second.
pub type MassFlowRate = Qty<0, 1, -1, 0, 0, 0>;
/// Cubic metres per second.
pub type VolumetricFlowRate = Qty<3, 0, -1, 0, 0, 0>;
/// Radians per second (angle is dimensionless).
pub type AngularVelocity = Frequency;
/// Angle in radians (dimensionless by SI convention; `deg` parses with
/// numeric conversion — see [`parse`]).
pub type Angle = Dimensionless;

/// Unit-bearing constructors, all returning coherent-SI values.
pub mod units {
    use super::{
        Amount, AmountConcentration, Angle, Capacitance, Conductance, DynViscosity, ElectricCharge,
        Energy, Force, Frequency, Inductance, Length, MagneticFlux, MagneticFluxDensity, Mass,
        MolarMass, Power, Pressure, Qty, Resistance, SurfaceTension, Temperature, Time, Velocity,
        Voltage, Volume, VolumetricFlowRate,
    };

    /// Metres.
    #[must_use]
    pub const fn meters(v: f64) -> Length {
        Qty(v)
    }
    /// Millimetres.
    #[must_use]
    pub const fn millimeters(v: f64) -> Length {
        Qty(v * 1e-3)
    }
    /// Seconds.
    #[must_use]
    pub const fn seconds(v: f64) -> Time {
        Qty(v)
    }
    /// Hours (3600 s).
    #[must_use]
    pub const fn hours(v: f64) -> Time {
        Qty(v * 3600.0)
    }
    /// Kilograms.
    #[must_use]
    pub const fn kilograms(v: f64) -> Mass {
        Qty(v)
    }
    /// Moles.
    #[must_use]
    pub const fn moles(v: f64) -> Amount {
        Qty(v)
    }
    /// Kilograms per mole.
    #[must_use]
    pub const fn kilograms_per_mole(v: f64) -> MolarMass {
        Qty(v)
    }
    /// Moles per cubic metre.
    #[must_use]
    pub const fn moles_per_cubic_meter(v: f64) -> AmountConcentration {
        Qty(v)
    }
    /// Kelvin.
    #[must_use]
    pub const fn kelvin(v: f64) -> Temperature {
        Qty(v)
    }
    /// Degrees Celsius, converted to kelvin (affine: K = °C + 273.15).
    #[must_use]
    pub const fn celsius(v: f64) -> Temperature {
        Qty(v + 273.15)
    }
    /// Newtons.
    #[must_use]
    pub const fn newtons(v: f64) -> Force {
        Qty(v)
    }
    /// Pascals.
    #[must_use]
    pub const fn pascals(v: f64) -> Pressure {
        Qty(v)
    }
    /// Pascal-seconds (dynamic viscosity).
    #[must_use]
    pub const fn pascal_seconds(v: f64) -> DynViscosity {
        Qty(v)
    }
    /// Joules.
    #[must_use]
    pub const fn joules(v: f64) -> Energy {
        Qty(v)
    }
    /// Watts.
    #[must_use]
    pub const fn watts(v: f64) -> Power {
        Qty(v)
    }
    /// Coulombs.
    #[must_use]
    pub const fn coulombs(v: f64) -> ElectricCharge {
        Qty(v)
    }
    /// Volts.
    #[must_use]
    pub const fn volts(v: f64) -> Voltage {
        Qty(v)
    }
    /// Webers.
    #[must_use]
    pub const fn webers(v: f64) -> MagneticFlux {
        Qty(v)
    }
    /// Henries.
    #[must_use]
    pub const fn henries(v: f64) -> Inductance {
        Qty(v)
    }
    /// Ohms.
    #[must_use]
    pub const fn ohms(v: f64) -> Resistance {
        Qty(v)
    }
    /// Siemens.
    #[must_use]
    pub const fn siemens(v: f64) -> Conductance {
        Qty(v)
    }
    /// Farads.
    #[must_use]
    pub const fn farads(v: f64) -> Capacitance {
        Qty(v)
    }
    /// Teslas.
    #[must_use]
    pub const fn teslas(v: f64) -> MagneticFluxDensity {
        Qty(v)
    }
    /// Hertz.
    #[must_use]
    pub const fn hertz(v: f64) -> Frequency {
        Qty(v)
    }
    /// Litres (1e-3 m³).
    #[must_use]
    pub const fn liters(v: f64) -> Volume {
        Qty(v * 1e-3)
    }
    /// Litres per second.
    #[must_use]
    pub const fn liters_per_second(v: f64) -> VolumetricFlowRate {
        Qty(v * 1e-3)
    }
    /// Metres per second.
    #[must_use]
    pub const fn meters_per_second(v: f64) -> Velocity {
        Qty(v)
    }
    /// Newtons per metre (surface tension).
    #[must_use]
    pub const fn newtons_per_meter(v: f64) -> SurfaceTension {
        Qty(v)
    }
    /// Radians (dimensionless).
    #[must_use]
    pub const fn radians(v: f64) -> Angle {
        Qty(v)
    }
    /// Degrees, converted to radians (dimensionless).
    #[must_use]
    pub fn degrees(v: f64) -> Angle {
        Qty(v * (core::f64::consts::PI / 180.0))
    }
}

// ---------------------------------------------------------------------------
// QtyAny: runtime-checked dimensioned value.
// ---------------------------------------------------------------------------

/// Error type for runtime dimension violations — a structured value carrying
/// both dimension vectors, per the errors-as-guidance doctrine (P10).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DimensionMismatch {
    /// Operation that failed (`"add"`, `"sub"`, `"convert"`).
    pub op: &'static str,
    /// Left/actual dimension.
    pub left: Dims,
    /// Right/expected dimension.
    pub right: Dims,
}

impl fmt::Display for DimensionMismatch {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "dimension mismatch in {}: [{}] vs [{}] — quantities must share a dimension; \
             check the unit annotations on the inputs",
            self.op,
            self.left.unit_string(),
            self.right.unit_string()
        )
    }
}

impl core::error::Error for DimensionMismatch {}

/// Runtime dimension arithmetic exceeded the representable exponent domain.
/// Refusing is safer than attaching saturated metadata to an unsaturated value.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DimensionOverflow {
    /// Operation that failed.
    pub op: &'static str,
    /// Dimensions of the left/base operand.
    pub dims: Dims,
    /// Dimensions of the right operand for binary composition
    /// (`None` for exponent scaling).
    pub rhs: Option<Dims>,
    /// Requested scale factor (`0` for binary composition).
    pub factor: i32,
}

impl fmt::Display for DimensionOverflow {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.rhs {
            Some(rhs) => write!(
                f,
                "dimension overflow in {}: [{}] composed with [{}] leaves the supported i8 \
                 exponent domain",
                self.op,
                self.dims.unit_string(),
                rhs.unit_string()
            ),
            None => write!(
                f,
                "dimension overflow in {}: [{}] scaled by {} leaves the supported i8 exponent \
                 domain",
                self.op,
                self.dims.unit_string(),
                self.factor
            ),
        }
    }
}

impl core::error::Error for DimensionOverflow {}

/// A runtime-dimensioned value: the dimension vector travels as data. Used
/// for IR values, data files, and anywhere compile-time `Qty` cannot reach
/// (G3 demands unit-rescaling invariance through RUNTIME paths too).
#[derive(Clone, Copy, PartialEq)]
pub struct QtyAny {
    /// Value in coherent SI base units.
    pub value: f64,
    /// SI base-unit exponents.
    pub dims: Dims,
}

impl QtyAny {
    /// A dimensionless value.
    #[must_use]
    pub const fn dimensionless(value: f64) -> Self {
        QtyAny {
            value,
            dims: Dims::NONE,
        }
    }

    /// Construct from value and dimension vector.
    #[must_use]
    pub const fn new(value: f64, dims: Dims) -> Self {
        QtyAny { value, dims }
    }

    /// Checked addition: dimensions must match.
    ///
    /// # Errors
    /// Returns [`DimensionMismatch`] when the dimensions differ.
    pub fn try_add(self, rhs: QtyAny) -> Result<QtyAny, DimensionMismatch> {
        if self.dims == rhs.dims {
            Ok(QtyAny {
                value: self.value + rhs.value,
                dims: self.dims,
            })
        } else {
            Err(DimensionMismatch {
                op: "add",
                left: self.dims,
                right: rhs.dims,
            })
        }
    }

    /// Checked subtraction: dimensions must match.
    ///
    /// # Errors
    /// Returns [`DimensionMismatch`] when the dimensions differ.
    pub fn try_sub(self, rhs: QtyAny) -> Result<QtyAny, DimensionMismatch> {
        if self.dims == rhs.dims {
            Ok(QtyAny {
                value: self.value - rhs.value,
                dims: self.dims,
            })
        } else {
            Err(DimensionMismatch {
                op: "sub",
                left: self.dims,
                right: rhs.dims,
            })
        }
    }

    /// Integer power with checked dimension scaling.
    ///
    /// # Errors
    /// [`DimensionOverflow`] when any scaled exponent would leave `i8`.
    pub fn powi(self, n: i8) -> Result<QtyAny, DimensionOverflow> {
        let scaled = self.dims.checked_times(n).ok_or(DimensionOverflow {
            op: "powi",
            dims: self.dims,
            rhs: None,
            factor: i32::from(n),
        })?;
        Ok(QtyAny {
            value: powi_pinned(self.value, i32::from(n)),
            dims: scaled,
        })
    }

    /// AUTHORITATIVE multiplication (bead sj31i.11): dimensions add
    /// with checked exponent arithmetic.
    ///
    /// # Errors
    /// [`DimensionOverflow`] when any summed exponent leaves `i8`.
    pub fn try_mul(self, rhs: QtyAny) -> Result<QtyAny, DimensionOverflow> {
        let dims = self.dims.checked_plus(rhs.dims).ok_or(DimensionOverflow {
            op: "mul",
            dims: self.dims,
            rhs: Some(rhs.dims),
            factor: 0,
        })?;
        Ok(QtyAny {
            value: self.value * rhs.value,
            dims,
        })
    }

    /// AUTHORITATIVE division (bead sj31i.11): dimensions subtract with
    /// checked exponent arithmetic.
    ///
    /// # Errors
    /// [`DimensionOverflow`] when any subtracted exponent leaves `i8`.
    pub fn try_div(self, rhs: QtyAny) -> Result<QtyAny, DimensionOverflow> {
        let dims = self.dims.checked_minus(rhs.dims).ok_or(DimensionOverflow {
            op: "div",
            dims: self.dims,
            rhs: Some(rhs.dims),
            factor: 0,
        })?;
        Ok(QtyAny {
            value: self.value / rhs.value,
            dims,
        })
    }

    /// Downcast to a compile-time `Qty`, checking the dimension.
    ///
    /// # Errors
    /// Returns [`DimensionMismatch`] when this value's dimension differs
    /// from the target type's.
    pub fn to_typed<
        const M: i8,
        const KG: i8,
        const S: i8,
        const K: i8,
        const A: i8,
        const MOL: i8,
    >(
        self,
    ) -> Result<Qty<M, KG, S, K, A, MOL>, DimensionMismatch> {
        let want = Qty::<M, KG, S, K, A, MOL>::DIMS;
        if self.dims == want {
            Ok(Qty(self.value))
        } else {
            Err(DimensionMismatch {
                op: "convert",
                left: self.dims,
                right: want,
            })
        }
    }
}

/// Multiplication: dimensions add with CHECKED exponent arithmetic
/// (bead sj31i.11). The operator form is a developer convenience that
/// PANICS on exponent overflow — agent-facing paths (parser, IR,
/// admission) must use [`QtyAny::try_mul`] for the typed refusal.
/// Silent saturation is gone: it could alias wrong physics back into
/// the valid exponent range after cancellation.
impl Mul for QtyAny {
    type Output = QtyAny;
    fn mul(self, rhs: QtyAny) -> QtyAny {
        self.try_mul(rhs)
            .expect("QtyAny `*` overflowed the i8 exponent domain; use try_mul for typed refusal")
    }
}

/// Division: dimensions subtract with CHECKED exponent arithmetic; the
/// operator form panics on overflow — see [`QtyAny::try_div`].
impl Div for QtyAny {
    type Output = QtyAny;
    fn div(self, rhs: QtyAny) -> QtyAny {
        self.try_div(rhs)
            .expect("QtyAny `/` overflowed the i8 exponent domain; use try_div for typed refusal")
    }
}

impl fmt::Debug for QtyAny {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} {}", self.value, self.dims.unit_string())
    }
}

impl fmt::Display for QtyAny {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Debug::fmt(self, f)
    }
}

/// Pinned-order integer power (LSB-first square-and-multiply, one final
/// reciprocal for ordinary negative results). A reciprocal-base recovery
/// preserves representable subnormals when the positive intermediate
/// overflows.
///
/// `f64::powi` is FORBIDDEN in value-producing paths: its rounding is
/// optimization-level-dependent (bead 4xnt), which would make parsed
/// unit scales and `QtyAny::powi` build-mode-dependent. This is a
/// deliberate private copy of `fs_math::det::powi`: fs-qty is UTIL and
/// may not depend on L0 fs-math (xtask layer rule), so the three-line
/// loop is duplicated rather than the layer boundary broken. Keep the
/// two implementations textually in sync.
pub(crate) fn powi_pinned(x: f64, n: i32) -> f64 {
    fn unsigned(mut a: f64, mut b: u64) -> f64 {
        let mut r = 1.0f64;
        loop {
            if b & 1 == 1 {
                r *= a;
            }
            b /= 2;
            if b == 0 {
                break;
            }
            a *= a;
        }
        r
    }

    let b = i64::from(n).unsigned_abs();
    let positive = unsigned(x, b);
    if n >= 0 {
        return positive;
    }
    let reciprocal = 1.0 / positive;
    if reciprocal == 0.0 && x.is_finite() && x != 0.0 {
        return unsigned(1.0 / x, b);
    }
    reciprocal
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn version_is_stamped() {
        assert!(!VERSION.is_empty());
    }

    #[test]
    fn typed_algebra_produces_correct_dimensions() {
        let rho = Density::new(1000.0);
        let v = Velocity::new(2.0);
        let l = Length::new(0.05);
        let mu = DynViscosity::new(1e-3);
        // Reynolds number: rho v L / mu must be dimensionless.
        let re: Dimensionless = rho * v * l / mu;
        assert!((re.value() - 1.0e5).abs() < 1e-6);
    }

    #[test]
    fn pressure_is_force_over_area() {
        let f = Force::new(10.0);
        let a = Area::new(2.0);
        let p: Pressure = f / a;
        assert!((p.value() - 5.0).abs() < 1e-12);
        assert_eq!(Pressure::DIMS, Dims([-1, 1, -2, 0, 0, 0]));
    }

    #[test]
    fn qty_any_checked_add_rejects_mismatch() {
        let p = Pressure::new(1.0).erase();
        let t = Time::new(1.0).erase();
        let err = p.try_add(t).unwrap_err();
        assert_eq!(err.op, "add");
        let msg = err.to_string();
        assert!(
            msg.contains("dimension mismatch"),
            "teaching message expected, got: {msg}"
        );
    }

    #[test]
    fn qty_any_power_refuses_dimension_saturation() {
        let q = QtyAny::new(2.0, Dims([2, 0, 0, 0, 0, 0]));
        // det-ok: QtyAny::powi — routes to powi_pinned (4xnt)
        let err = q.powi(64).expect_err("m^128 cannot fit runtime dims");
        assert_eq!(err.factor, 64);

        let admitted = QtyAny::new(2.0, Dims([1, 0, 0, 0, 0, 0]))
            // det-ok: QtyAny::powi — routes to powi_pinned (4xnt)
            .powi(101)
            .expect("m^101 is representable");
        assert_eq!(admitted.dims, Dims([101, 0, 0, 0, 0, 0]));
        assert_eq!(powi_pinned(2.0, -1024).to_bits(), 1_u64 << 50);
    }

    #[test]
    fn qty_any_algebra_matches_typed_algebra() {
        // Property loop over a deterministic value grid: erased algebra must
        // agree with typed algebra bit-for-bit (same f64 operations).
        for i in 1..50u32 {
            let x = f64::from(i) * 0.37;
            let y = f64::from(i).mul_add(1.11, 0.5);
            let typed = (Length::new(x) / Time::new(y)).value();
            let erased = Length::new(x).erase() / Time::new(y).erase();
            assert_eq!(erased.value.to_bits(), typed.to_bits());
            assert_eq!(erased.dims, Velocity::DIMS);
        }
    }

    #[test]
    fn round_trip_typed_erase_downcast() {
        let mu = DynViscosity::new(0.12);
        let any = mu.erase();
        let back: DynViscosity = any.to_typed().expect("dims match");
        assert_eq!(back.value().to_bits(), mu.value().to_bits());
        let wrong: Result<Pressure, _> = any.to_typed();
        assert!(wrong.is_err());
    }

    #[test]
    fn unit_strings_are_canonical() {
        assert_eq!(Pressure::DIMS.unit_string(), "kg·m^-1·s^-2");
        assert_eq!(Dimensionless::DIMS.unit_string(), "1");
        assert_eq!(Velocity::DIMS.unit_string(), "m·s^-1");
    }

    #[test]
    fn total_cmp_sorts_nan_deterministically() {
        let mut v = [Length::new(f64::NAN), Length::new(1.0), Length::new(-1.0)];
        v.sort_by(Qty::total_cmp);
        assert_eq!(v[0].value().to_bits(), (-1.0f64).to_bits());
        assert_eq!(v[1].value().to_bits(), 1.0f64.to_bits());
        assert!(v[2].value().is_nan());
    }

    #[test]
    fn units_module_constructors() {
        assert!((units::millimeters(3.0).value() - 0.003).abs() < 1e-15);
        assert!((units::hours(2.0).value() - 7200.0).abs() < 1e-9);
        assert_eq!(units::moles(2.0).erase().dims, Amount::DIMS);
        assert!((units::celsius(20.0).value() - 293.15).abs() < 1e-9);
        assert!((units::degrees(180.0).value() - core::f64::consts::PI).abs() < 1e-12);
        assert!((units::liters_per_second(0.5).value() - 5e-4).abs() < 1e-15);
    }
}
