//! Original IAPWS ice Ih Gibbs equation and pressure-dependent melting curve.
//!
//! Ice properties use R10-06(2009), equation (1), Tables 1--4, with its
//! IAPWS-95-compatible entropy reference. The melting correlation is the
//! independent R14-08(2011), equation (1). Attribution: International
//! Association for the Properties of Water and Steam; see the source URLs.
//! No coefficient, enthalpy or entropy datum is adjusted to join another EOS.
//!
//! This evaluator assumes hexagonal ice Ih. Its deliberately restricted
//! numerical chart is 50..=273.16 K and 0..=210 MPa. It does not select a stable
//! phase against liquid/vapor or other ice polymorphs; consumers must admit
//! their own equilibrium path. The release's low-temperature experimental
//! limitations and physical uncertainties are distinct from solve residuals.

use super::MassSpecificVolume;
use super::if97::{
    IsothermalCompressibility, SpecificVolumePressureDerivative,
    SpecificVolumeTemperatureDerivative, ThermalExpansion,
};
use crate::{MassSpecificEnergyQuantityV1, MassSpecificThermalQuantityV1};
use core::fmt;
use core::ops::{Add, Div, Mul, Neg, Sub};
use fs_math::det;
use fs_qty::{Density, Pressure, Qty, Temperature};

/// Outward-rounded whole-path Cp bounds for conservative thermal integration.
pub mod bounds;

/// Original coefficient release and temperature scale.
pub const SOURCE_RELEASE: &str = "IAPWS R10-06(2009), Ice Ih, ITS-90";
/// Original source equation, coefficients, derivative and verification tables.
pub const SOURCE_URL: &str = "https://iapws.org/technical-guidance/release/Ice-2009.download";
/// Original melting-pressure correlation source.
pub const MELTING_SOURCE_RELEASE: &str = "IAPWS R14-08(2011), ice Ih melting, ITS-90";
/// Original melting-pressure equation and verification table.
pub const MELTING_SOURCE_URL: &str =
    "https://iapws.org/technical-guidance/release/MeltSub.download";
/// Deterministic original-source operation tree.
pub const EVALUATOR_VERSION: u32 = 1;
/// Lower admitted numerical temperature [K]; no zero-temperature extrapolation.
pub const MINIMUM_TEMPERATURE_K: f64 = 50.0;
/// Triple-point temperature and upper admitted ice chart temperature [K].
pub const TRIPLE_POINT_TEMPERATURE_K: f64 = 273.16;
/// Experimental triple-point pressure used by both released equations [Pa].
pub const TRIPLE_POINT_PRESSURE_PA: f64 = 611.657;
/// Normal pressure used by the pressure polynomials [Pa].
pub const NORMAL_PRESSURE_PA: f64 = 101_325.0;
/// Upper pressure of this explicitly assumed-Ih source chart [Pa].
pub const MAXIMUM_PRESSURE_PA: f64 = 210e6;
/// Ice Ih/ice III/liquid endpoint of the melting correlation [K].
pub const MELTING_MINIMUM_TEMPERATURE_K: f64 = 251.165;
/// Hard work limit in addition to the caller's numerical solve allowance.
pub const MAX_MELTING_ITERATIONS: usize = 128;

// Table 2, retaining the IAPWS-95 caloric convention exactly as released.
const G_COEFFICIENTS: [f64; 5] = [
    -0.632_020_233_335_886e6,
    0.655_022_213_658_955,
    -0.189_369_929_326_131e-7,
    0.339_746_123_271_053e-14,
    -0.556_464_869_058_991e-21,
];
const S0: f64 = -0.332_733_756_492_168e4;
const T_POLES: [Complex; 2] = [
    Complex::new(0.368_017_112_855_051e-1, 0.510_878_114_959_572e-1),
    Complex::new(0.337_315_741_065_416, 0.335_449_415_919_309),
];
const R1: Complex = Complex::new(0.447_050_716_285_388e2, 0.656_876_847_463_481e2);
const R2_COEFFICIENTS: [Complex; 3] = [
    Complex::new(-0.725_974_574_329_220e2, -0.781_008_427_112_870e2),
    Complex::new(-0.557_107_698_030_123e-4, 0.464_578_634_580_806e-4),
    Complex::new(0.234_801_409_215_913e-10, -0.285_651_142_904_972e-10),
];
const MELTING_TERMS: [(f64, f64); 3] = [
    (0.119_539_337e7, 0.300_000e1),
    (0.808_183_159e5, 0.257_500e2),
    (0.333_826_860e4, 0.103_750e3),
];

/// Source-domain, numerical representation, work or cancellation refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum IceIhError {
    Invalid(&'static str),
    OutsideDomain {
        pressure_pa: f64,
        temperature_k: f64,
    },
    Unrepresentable(&'static str),
    NotConverged {
        iterations: usize,
        pressure_residual_pa: f64,
    },
    Cancelled,
}
impl fmt::Display for IceIhError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "IAPWS ice Ih: {self:?}")
    }
}
impl std::error::Error for IceIhError {}

/// Second temperature derivative of specific Gibbs energy [J/(kg K²)].
pub type SpecificGibbsTemperatureSecondDerivative = Qty<2, 0, -2, -2, 0, 0>;
/// Ice pressure coefficient at fixed specific volume [Pa/K].
pub type PressureTemperatureDerivative = Qty<-1, 1, -2, -1, 0, 0>;

/// Immutable original Gibbs source state, explicitly assuming the ice Ih phase.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IceIhState {
    pressure_pa: f64,
    temperature_k: f64,
    gibbs_j_kg: f64,
    enthalpy_j_kg: f64,
    internal_energy_j_kg: f64,
    entropy_j_kg_k: f64,
    volume_m3_kg: f64,
    cp_j_kg_k: f64,
    cv_j_kg_k: f64,
    volume_temperature_derivative: f64,
    volume_pressure_derivative: f64,
}
impl IceIhState {
    /// Original Gibbs evaluation without allocation or inferred phase selection.
    pub fn at_pressure_temperature(
        pressure: Pressure,
        temperature: Temperature,
    ) -> Result<Self, IceIhError> {
        let (p, t) = (pressure.value(), temperature.value());
        admit(p, t)?;
        let tau = t / TRIPLE_POINT_TEMPERATURE_K;
        let (g0, g0p, g0pp, r2, r2p, r2pp) = pressure_coefficients(p);
        let mut h_term = 0.0;
        let mut derivative_term = 0.0;
        let mut second_term = 0.0;
        let mut f2 = Complex::ZERO;
        let mut fp2 = Complex::ZERO;
        for (index, (pole, r)) in T_POLES.into_iter().zip([R1, r2]).enumerate() {
            let (h, fp, fpp) = thermal_terms(pole, tau);
            h_term += (r * h).re;
            derivative_term += (r * fp).re;
            second_term += (r * fpp).re;
            if index == 1 {
                f2 = h + fp * tau;
                fp2 = fp;
            }
        }
        let enthalpy = g0 + TRIPLE_POINT_TEMPERATURE_K * h_term;
        let entropy = S0 - derivative_term;
        let gibbs = enthalpy - t * entropy;
        let volume = g0p + TRIPLE_POINT_TEMPERATURE_K * (r2p * f2).re;
        let vt = (r2p * fp2).re;
        let vp = g0pp + TRIPLE_POINT_TEMPERATURE_K * (r2pp * f2).re;
        let cp = -tau * second_term;
        let cv = cp + t * vt * vt / vp;
        let internal = enthalpy - p * volume;
        if ![gibbs, enthalpy, entropy, volume, vt, vp, cp, cv, internal]
            .into_iter()
            .all(f64::is_finite)
            || volume <= 0.0
            || vp >= 0.0
            || cp <= 0.0
            || cv <= 0.0
        {
            return Err(IceIhError::Unrepresentable(
                "ice Gibbs state or stable response",
            ));
        }
        Ok(Self {
            pressure_pa: p,
            temperature_k: t,
            gibbs_j_kg: gibbs,
            enthalpy_j_kg: enthalpy,
            internal_energy_j_kg: internal,
            entropy_j_kg_k: entropy,
            volume_m3_kg: volume,
            cp_j_kg_k: cp,
            cv_j_kg_k: cv,
            volume_temperature_derivative: vt,
            volume_pressure_derivative: vp,
        })
    }
    pub const fn pressure(self) -> Pressure {
        Pressure::new(self.pressure_pa)
    }
    pub const fn temperature(self) -> Temperature {
        Temperature::new(self.temperature_k)
    }
    pub const fn specific_gibbs_energy(self) -> MassSpecificEnergyQuantityV1 {
        MassSpecificEnergyQuantityV1::new(self.gibbs_j_kg)
    }
    pub const fn specific_enthalpy(self) -> MassSpecificEnergyQuantityV1 {
        MassSpecificEnergyQuantityV1::new(self.enthalpy_j_kg)
    }
    pub const fn specific_internal_energy(self) -> MassSpecificEnergyQuantityV1 {
        MassSpecificEnergyQuantityV1::new(self.internal_energy_j_kg)
    }
    pub const fn specific_entropy(self) -> MassSpecificThermalQuantityV1 {
        MassSpecificThermalQuantityV1::new(self.entropy_j_kg_k)
    }
    pub const fn specific_volume(self) -> MassSpecificVolume {
        MassSpecificVolume::new(self.volume_m3_kg)
    }
    pub fn density(self) -> Density {
        Density::new(1.0 / self.volume_m3_kg)
    }
    pub const fn cp(self) -> MassSpecificThermalQuantityV1 {
        MassSpecificThermalQuantityV1::new(self.cp_j_kg_k)
    }
    pub const fn cv(self) -> MassSpecificThermalQuantityV1 {
        MassSpecificThermalQuantityV1::new(self.cv_j_kg_k)
    }
    pub fn specific_helmholtz_energy(self) -> MassSpecificEnergyQuantityV1 {
        MassSpecificEnergyQuantityV1::new(self.gibbs_j_kg - self.pressure_pa * self.volume_m3_kg)
    }
    pub const fn volume_temperature_derivative(self) -> SpecificVolumeTemperatureDerivative {
        SpecificVolumeTemperatureDerivative::new(self.volume_temperature_derivative)
    }
    pub const fn volume_pressure_derivative(self) -> SpecificVolumePressureDerivative {
        SpecificVolumePressureDerivative::new(self.volume_pressure_derivative)
    }
    pub fn gibbs_temperature_second_derivative(self) -> SpecificGibbsTemperatureSecondDerivative {
        SpecificGibbsTemperatureSecondDerivative::new(-self.cp_j_kg_k / self.temperature_k)
    }
    pub fn thermal_expansion(self) -> ThermalExpansion {
        ThermalExpansion::new(self.volume_temperature_derivative / self.volume_m3_kg)
    }
    pub fn pressure_coefficient(self) -> PressureTemperatureDerivative {
        PressureTemperatureDerivative::new(
            -self.volume_temperature_derivative / self.volume_pressure_derivative,
        )
    }
    pub fn isothermal_compressibility(self) -> IsothermalCompressibility {
        IsothermalCompressibility::new(-self.volume_pressure_derivative / self.volume_m3_kg)
    }
    pub fn isentropic_compressibility(self) -> IsothermalCompressibility {
        IsothermalCompressibility::new(
            (-self.volume_pressure_derivative
                - self.temperature_k * self.volume_temperature_derivative.powi(2) / self.cp_j_kg_k)
                / self.volume_m3_kg,
        )
    }
}

fn admit(p: f64, t: f64) -> Result<(), IceIhError> {
    if !p.is_finite() || !t.is_finite() {
        return Err(IceIhError::Invalid(
            "finite pressure and temperature required",
        ));
    }
    if !(0.0..=MAXIMUM_PRESSURE_PA).contains(&p)
        || !(MINIMUM_TEMPERATURE_K..=TRIPLE_POINT_TEMPERATURE_K).contains(&t)
    {
        return Err(IceIhError::OutsideDomain {
            pressure_pa: p,
            temperature_k: t,
        });
    }
    Ok(())
}

fn thermal_coefficients(p: f64) -> [Complex; 2] {
    let x = (p - NORMAL_PRESSURE_PA) / TRIPLE_POINT_PRESSURE_PA;
    [
        R1,
        (R2_COEFFICIENTS[2] * x + R2_COEFFICIENTS[1]) * x + R2_COEFFICIENTS[0],
    ]
}
fn pressure_coefficients(p: f64) -> (f64, f64, f64, Complex, Complex, Complex) {
    let x = (p - NORMAL_PRESSURE_PA) / TRIPLE_POINT_PRESSURE_PA;
    let c = G_COEFFICIENTS;
    let g = (((c[4] * x + c[3]) * x + c[2]) * x + c[1]) * x + c[0];
    let gp = ((4.0 * c[4] * x + 3.0 * c[3]) * x + 2.0 * c[2]) * x + c[1];
    let gpp = (12.0 * c[4] * x + 6.0 * c[3]) * x + 2.0 * c[2];
    let r2 = thermal_coefficients(p)[1];
    let r2p = (R2_COEFFICIENTS[2] * (2.0 * x) + R2_COEFFICIENTS[1]) / TRIPLE_POINT_PRESSURE_PA;
    let r2pp = R2_COEFFICIENTS[2] * (2.0 / (TRIPLE_POINT_PRESSURE_PA * TRIPLE_POINT_PRESSURE_PA));
    (
        g,
        gp / TRIPLE_POINT_PRESSURE_PA,
        gpp / (TRIPLE_POINT_PRESSURE_PA * TRIPLE_POINT_PRESSURE_PA),
        r2,
        r2p,
        r2pp,
    )
}

// Return F-tau F', F', F''. Algebraically equivalent to release Table 4.
// The reciprocal identity avoids cancellation in F'' at colder temperatures.
fn thermal_terms(pole: Complex, tau: f64) -> (Complex, Complex, Complex) {
    let ratio = Complex::new(tau, 0.0) / pole;
    let minus = (-ratio).ln_1p();
    let plus = ratio.ln_1p();
    let h = pole * (minus + plus) + ratio * tau;
    let fp = plus - minus - ratio * 2.0;
    let fpp =
        Complex::new(2.0 * tau * tau, 0.0) / (pole * (pole * pole - Complex::new(tau * tau, 0.0)));
    (h, fp, fpp)
}

/// Stable changes on one original ice Ih isobar; no large reference subtraction.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IceIhDifference {
    enthalpy_j_kg: f64,
    internal_energy_j_kg: f64,
    entropy_j_kg_k: f64,
    volume_m3_kg: f64,
}
impl IceIhDifference {
    pub const fn specific_enthalpy_change(self) -> MassSpecificEnergyQuantityV1 {
        MassSpecificEnergyQuantityV1::new(self.enthalpy_j_kg)
    }
    pub const fn specific_internal_energy_change(self) -> MassSpecificEnergyQuantityV1 {
        MassSpecificEnergyQuantityV1::new(self.internal_energy_j_kg)
    }
    pub const fn specific_entropy_change(self) -> MassSpecificThermalQuantityV1 {
        MassSpecificThermalQuantityV1::new(self.entropy_j_kg_k)
    }
    pub const fn specific_volume_change(self) -> MassSpecificVolume {
        MassSpecificVolume::new(self.volume_m3_kg)
    }
}

/// Integrate original Cp, Cp/T and volume response between two admitted temperatures.
pub fn isobaric_difference(
    pressure: Pressure,
    initial: Temperature,
    final_temperature: Temperature,
) -> Result<IceIhDifference, IceIhError> {
    let (p, ta, tb) = (pressure.value(), initial.value(), final_temperature.value());
    admit(p, ta)?;
    admit(p, tb)?;
    if ta == tb {
        return Ok(IceIhDifference {
            enthalpy_j_kg: 0.0,
            internal_energy_j_kg: 0.0,
            entropy_j_kg_k: 0.0,
            volume_m3_kg: 0.0,
        });
    }
    let (ta, tb, sign) = if ta < tb {
        (ta, tb, 1.0)
    } else {
        (tb, ta, -1.0)
    };
    let a = ta / TRIPLE_POINT_TEMPERATURE_K;
    let dtau = (tb - ta) / TRIPLE_POINT_TEMPERATURE_K;
    let square_delta =
        ((tb - ta) / TRIPLE_POINT_TEMPERATURE_K) * ((tb + ta) / TRIPLE_POINT_TEMPERATURE_K);
    let mut dh = 0.0;
    let mut ds = 0.0;
    let mut df2 = Complex::ZERO;
    for (index, (pole, r)) in T_POLES.into_iter().zip(thermal_coefficients(p)).enumerate() {
        let h_delta = pole
            * (-Complex::new(square_delta, 0.0) / (pole * pole - Complex::new(a * a, 0.0))).ln_1p()
            + Complex::new(square_delta, 0.0) / pole;
        let fp_delta = (Complex::new(dtau, 0.0) / (pole + Complex::new(a, 0.0))).ln_1p()
            - (-Complex::new(dtau, 0.0) / (pole - Complex::new(a, 0.0))).ln_1p()
            - Complex::new(2.0 * dtau, 0.0) / pole;
        dh += (r * h_delta).re;
        ds -= (r * fp_delta).re;
        if index == 1 {
            let fp_b = thermal_terms(pole, tb / TRIPLE_POINT_TEMPERATURE_K).1;
            df2 = h_delta + fp_delta * a + fp_b * dtau;
        }
    }
    let r2p = pressure_coefficients(p).4;
    let dh = sign * TRIPLE_POINT_TEMPERATURE_K * dh;
    let ds = sign * ds;
    let dv = sign * TRIPLE_POINT_TEMPERATURE_K * (r2p * df2).re;
    let du = dh - p * dv;
    if ![dh, du, ds, dv].into_iter().all(f64::is_finite)
        || dh == 0.0
        || ds == 0.0
        || (dh > 0.0) != (sign > 0.0)
    {
        return Err(IceIhError::Unrepresentable("isobaric ice difference"));
    }
    Ok(IceIhDifference {
        enthalpy_j_kg: dh,
        internal_energy_j_kg: du,
        entropy_j_kg_k: ds,
        volume_m3_kg: dv,
    })
}

/// Original R14-08(2011) ice Ih melting pressure, on 251.165..=273.16 K.
pub fn melting_pressure(temperature: Temperature) -> Result<Pressure, IceIhError> {
    melting_value_slope(temperature.value()).map(|(p, _)| Pressure::new(p))
}
/// Analytic derivative of the original melting pressure correlation [Pa/K].
pub fn melting_pressure_derivative(
    temperature: Temperature,
) -> Result<PressureTemperatureDerivative, IceIhError> {
    melting_value_slope(temperature.value()).map(|(_, dp)| PressureTemperatureDerivative::new(dp))
}
fn melting_value_slope(t: f64) -> Result<(f64, f64), IceIhError> {
    if !t.is_finite() || !(MELTING_MINIMUM_TEMPERATURE_K..=TRIPLE_POINT_TEMPERATURE_K).contains(&t)
    {
        return Err(IceIhError::Invalid(
            "ice Ih melting temperature outside 251.165..=273.16 K",
        ));
    }
    let x = (t - TRIPLE_POINT_TEMPERATURE_K) / TRIPLE_POINT_TEMPERATURE_K;
    let log_theta = ln_1p(x);
    let mut reduced = 1.0;
    let mut derivative = 0.0;
    for (a, b) in MELTING_TERMS {
        reduced -= a * det::expm1(b * log_theta);
        derivative -= a * b * det::exp((b - 1.0) * log_theta);
    }
    let p = TRIPLE_POINT_PRESSURE_PA * reduced;
    let dp = (TRIPLE_POINT_PRESSURE_PA / TRIPLE_POINT_TEMPERATURE_K) * derivative;
    if !p.is_finite() || !dp.is_finite() || p <= 0.0 || dp >= 0.0 {
        return Err(IceIhError::Unrepresentable("melting curve"));
    }
    Ok((p, dp))
}

/// Explicit residual, bracket-width and work budgets for the original melting inverse.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IceMeltingSolve {
    pub pressure_tolerance_pa: f64,
    pub temperature_tolerance_k: f64,
    pub maximum_iterations: usize,
}
/// Accepted original correlation and its retained inverse diagnostics.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IceMeltingSolution {
    temperature_k: f64,
    actual_pressure_pa: f64,
    requested_pressure_pa: f64,
    temperature_bounds_k: [f64; 2],
    iterations: usize,
}
impl IceMeltingSolution {
    pub const fn temperature(self) -> Temperature {
        Temperature::new(self.temperature_k)
    }
    pub const fn pressure(self) -> Pressure {
        Pressure::new(self.actual_pressure_pa)
    }
    pub const fn requested_pressure(self) -> Pressure {
        Pressure::new(self.requested_pressure_pa)
    }
    pub fn pressure_residual_pa(self) -> f64 {
        self.actual_pressure_pa - self.requested_pressure_pa
    }
    /// Bracket of floating-point original correlation values, not a real-arithmetic enclosure.
    pub const fn temperature_bounds_k(self) -> [f64; 2] {
        self.temperature_bounds_k
    }
    pub const fn iterations(self) -> usize {
        self.iterations
    }
}

/// Bounded safeguarded Newton/bisection on the decreasing melting-pressure curve.
/// Both residual and temperature-width budgets must pass. Exact forward anchors
/// at the source endpoints and 273.15 K are retained without clipping.
pub fn melting_temperature(
    pressure: Pressure,
    controls: IceMeltingSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<IceMeltingSolution, IceIhError> {
    poll(&mut cancelled)?;
    let p = pressure.value();
    if !p.is_finite()
        || !controls.pressure_tolerance_pa.is_finite()
        || controls.pressure_tolerance_pa < 0.0
        || !controls.temperature_tolerance_k.is_finite()
        || controls.temperature_tolerance_k < 0.0
        || !(1..=MAX_MELTING_ITERATIONS).contains(&controls.maximum_iterations)
    {
        return Err(IceIhError::Invalid(
            "finite nonnegative inverse budgets and 1..=128 iterations required",
        ));
    }
    let mut lo = MELTING_MINIMUM_TEMPERATURE_K;
    let mut hi = TRIPLE_POINT_TEMPERATURE_K;
    let maximum_p = melting_value_slope(lo)?.0;
    if p < TRIPLE_POINT_PRESSURE_PA || p > maximum_p {
        return Err(IceIhError::OutsideDomain {
            pressure_pa: p,
            temperature_k: f64::NAN,
        });
    }
    for t in [lo, 273.15, hi] {
        let value = melting_value_slope(t)?.0;
        if value == p {
            poll(&mut cancelled)?;
            return Ok(IceMeltingSolution {
                temperature_k: t,
                actual_pressure_pa: value,
                requested_pressure_pa: p,
                temperature_bounds_k: [t, t],
                iterations: 0,
            });
        }
    }
    let mut t = hi + (p - TRIPLE_POINT_PRESSURE_PA) / melting_value_slope(hi)?.1;
    let mut best = (f64::INFINITY, t, f64::NAN);
    let mut used = 0;
    for iteration in 1..=controls.maximum_iterations {
        poll(&mut cancelled)?;
        if !(t > lo && t < hi) {
            t = lo + 0.5 * (hi - lo);
        }
        if !(t > lo && t < hi) {
            break;
        }
        let (value, slope) = melting_value_slope(t)?;
        let residual = value - p;
        used = iteration;
        if residual.abs() < best.0 {
            best = (residual.abs(), t, value);
        }
        if residual == 0.0 {
            lo = t;
            hi = t;
        } else if residual > 0.0 {
            lo = t;
        } else {
            hi = t;
        }
        if residual.abs() <= controls.pressure_tolerance_pa
            && hi - lo <= controls.temperature_tolerance_k
        {
            poll(&mut cancelled)?;
            return Ok(IceMeltingSolution {
                temperature_k: t,
                actual_pressure_pa: value,
                requested_pressure_pa: p,
                temperature_bounds_k: [lo, hi],
                iterations: iteration,
            });
        }
        let newton = t - residual / slope;
        t = if newton > lo && newton < hi && newton != t {
            newton
        } else {
            lo + 0.5 * (hi - lo)
        };
    }
    if best.0 <= controls.pressure_tolerance_pa && hi - lo <= controls.temperature_tolerance_k {
        poll(&mut cancelled)?;
        return Ok(IceMeltingSolution {
            temperature_k: best.1,
            actual_pressure_pa: best.2,
            requested_pressure_pa: p,
            temperature_bounds_k: [lo, hi],
            iterations: used,
        });
    }
    Err(IceIhError::NotConverged {
        iterations: used,
        pressure_residual_pa: best.2 - p,
    })
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), IceIhError> {
    if cancelled() {
        Err(IceIhError::Cancelled)
    } else {
        Ok(())
    }
}

// A fixed-size native complex helper is sufficient for the two released poles.
// All logarithms are principal values using the repository's deterministic math.
#[derive(Clone, Copy, Debug, PartialEq)]
struct Complex {
    re: f64,
    im: f64,
}
impl Complex {
    const ZERO: Self = Self::new(0.0, 0.0);
    const fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }
    fn ln_1p(self) -> Self {
        let norm_delta = self.re * (2.0 + self.re) + self.im * self.im;
        Self::new(0.5 * ln_1p(norm_delta), det::atan2(self.im, 1.0 + self.re))
    }
}
impl Add for Complex {
    type Output = Self;
    fn add(self, r: Self) -> Self {
        Self::new(self.re + r.re, self.im + r.im)
    }
}
impl Sub for Complex {
    type Output = Self;
    fn sub(self, r: Self) -> Self {
        Self::new(self.re - r.re, self.im - r.im)
    }
}
impl Neg for Complex {
    type Output = Self;
    fn neg(self) -> Self {
        Self::new(-self.re, -self.im)
    }
}
impl Mul for Complex {
    type Output = Self;
    fn mul(self, r: Self) -> Self {
        Self::new(
            self.re * r.re - self.im * r.im,
            self.re * r.im + self.im * r.re,
        )
    }
}
impl Div for Complex {
    type Output = Self;
    fn div(self, r: Self) -> Self {
        let d = r.re * r.re + r.im * r.im;
        Self::new(
            (self.re * r.re + self.im * r.im) / d,
            (self.im * r.re - self.re * r.im) / d,
        )
    }
}
impl Mul<f64> for Complex {
    type Output = Self;
    fn mul(self, r: f64) -> Self {
        Self::new(self.re * r, self.im * r)
    }
}
impl Div<f64> for Complex {
    type Output = Self;
    fn div(self, r: f64) -> Self {
        Self::new(self.re / r, self.im / r)
    }
}

fn ln_1p(x: f64) -> f64 {
    if x.abs() >= 0.5 {
        return det::ln(1.0 + x);
    }
    let y = x / (2.0 + x);
    let y2 = y * y;
    let mut term = y;
    let mut sum = y;
    for k in 1..=20 {
        term *= y2;
        sum += term / (2 * k + 1) as f64;
    }
    2.0 * sum
}

#[cfg(test)]
mod tests;
