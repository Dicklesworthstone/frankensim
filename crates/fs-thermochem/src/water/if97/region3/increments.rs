//! Stable caloric differences on a prepared original Region 3 isochore.
//!
//! A bounded outward-interval proof admits the complete temperature path and
//! positive Cv once. Repeated differences reuse that immutable result, then
//! factor the original integer powers before subtracting absolute energies.
//! No source coefficients, phase references or caloric states are adjusted.

use super::{
    CRITICAL_DENSITY_KG_M3, CRITICAL_TEMPERATURE_K, If97DenseError, If97DenseState,
    LOG_COEFFICIENT, MAXIMUM_PRESSURE_PA, R_SPECIFIC_J_KG_K, TERMS,
};
use crate::{MassSpecificEnergyQuantityV1, MassSpecificThermalQuantityV1};
use core::fmt;
use fs_qty::{Density, Pressure, Temperature};

/// Hard limit in addition to the caller's explicit interval work allowance.
pub const MAX_INTERVAL_EVALUATIONS: usize = 65_535;
const MAX_DEPTH: u8 = 24;

/// Work allowance for preparing a whole fixed-density temperature interval.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DenseIsochoricIntervalSolve {
    pub maximum_intervals: usize,
}

/// Failed admission never publishes a partial path or an incomplete Cv bound.
#[derive(Clone, Debug, PartialEq)]
pub enum DenseIsochoricIntervalError {
    Source(If97DenseError),
    Invalid(&'static str),
    NotProven {
        interval_evaluations: usize,
        lower_temperature_k: f64,
        upper_temperature_k: f64,
    },
    Cancelled,
}
impl From<If97DenseError> for DenseIsochoricIntervalError {
    fn from(error: If97DenseError) -> Self {
        Self::Source(error)
    }
}
impl fmt::Display for DenseIsochoricIntervalError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "IF97 dense isochoric interval: {self:?}")
    }
}
impl std::error::Error for DenseIsochoricIntervalError {}

/// Proven original-source interval and positive Cv bounds at one density.
///
/// The enclosure concerns the retained binary64 coefficient model, not source
/// fit or experimental uncertainty. Boundary intervals can refuse when their
/// exact inclusion cannot be proved with the permitted work; no tolerance band
/// silently widens B23, Region 4, stability or the 100 MPa limit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DenseIsochoricInterval {
    density_kg_m3: f64,
    lower_temperature_k: f64,
    upper_temperature_k: f64,
    minimum_cv_j_kg_k: f64,
    maximum_cv_j_kg_k: f64,
    interval_evaluations: usize,
}
impl DenseIsochoricInterval {
    pub fn new(
        density: Density,
        lower: Temperature,
        upper: Temperature,
        controls: DenseIsochoricIntervalSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DenseIsochoricIntervalError> {
        use DenseIsochoricIntervalError as Error;
        if cancelled() {
            return Err(Error::Cancelled);
        }
        let (rho, lo, hi) = (density.value(), lower.value(), upper.value());
        if !lo.is_finite() || !hi.is_finite() || lo > hi {
            return Err(Error::Invalid("finite ordered temperature bounds required"));
        }
        if !(1..=MAX_INTERVAL_EVALUATIONS).contains(&controls.maximum_intervals) {
            return Err(Error::Invalid("invalid interval evaluation allowance"));
        }
        If97DenseState::at_density_temperature(density, lower)?;
        If97DenseState::at_density_temperature(density, upper)?;
        if rho == CRITICAL_DENSITY_KG_M3
            && lo <= CRITICAL_TEMPERATURE_K
            && hi >= CRITICAL_TEMPERATURE_K
        {
            return Err(If97DenseError::CriticalSingularity.into());
        }
        // Split the source's distinct subcritical saturation-side condition
        // before adaptive subdivision. The fixed-size depth-first stack bounds
        // storage independently of the caller's total evaluation allowance.
        let mut stack = [(0.0, 0.0, 0_u8); MAX_DEPTH as usize + 3];
        let mut pending = if lo < CRITICAL_TEMPERATURE_K && hi > CRITICAL_TEMPERATURE_K {
            stack[0] = (CRITICAL_TEMPERATURE_K, hi, 0);
            stack[1] = (lo, CRITICAL_TEMPERATURE_K, 0);
            2
        } else {
            stack[0] = (lo, hi, 0);
            1
        };
        let mut work = 0;
        let mut minimum = f64::INFINITY;
        let mut maximum = 0.0_f64;
        let polynomial = IsochoricPolynomial::new(rho);
        while pending != 0 {
            if cancelled() {
                return Err(Error::Cancelled);
            }
            pending -= 1;
            let (a, b, depth) = stack[pending];
            let unresolved = |work| Error::NotProven {
                interval_evaluations: work,
                lower_temperature_k: a,
                upper_temperature_k: b,
            };
            if work == controls.maximum_intervals {
                return Err(unresolved(work));
            }
            work += 1;
            let bounds = path_bounds(rho, a, b, &polynomial);
            if bounds.admitted && bounds.cv.hi <= 2.0 * bounds.cv.lo {
                minimum = minimum.min(bounds.cv.lo);
                maximum = maximum.max(bounds.cv.hi);
                continue;
            }
            let midpoint = a + 0.5 * (b - a);
            if depth == MAX_DEPTH || midpoint <= a || midpoint >= b {
                return Err(unresolved(work));
            }
            stack[pending] = (midpoint, b, depth + 1);
            stack[pending + 1] = (a, midpoint, depth + 1);
            pending += 2;
        }
        if cancelled() {
            return Err(Error::Cancelled);
        }
        Ok(Self {
            density_kg_m3: rho,
            lower_temperature_k: lo,
            upper_temperature_k: hi,
            minimum_cv_j_kg_k: minimum,
            maximum_cv_j_kg_k: maximum,
            interval_evaluations: work,
        })
    }
    pub const fn density(&self) -> Density {
        Density::new(self.density_kg_m3)
    }
    pub const fn temperature_bounds_k(&self) -> [f64; 2] {
        [self.lower_temperature_k, self.upper_temperature_k]
    }
    pub const fn minimum_heat_capacity_j_kg_k(&self) -> f64 {
        self.minimum_cv_j_kg_k
    }
    pub const fn maximum_heat_capacity_j_kg_k(&self) -> f64 {
        self.maximum_cv_j_kg_k
    }
    /// Proven upper bound on dT/du over this entire original-source interval.
    pub fn maximum_temperature_slope(&self) -> f64 {
        (1.0 / self.minimum_cv_j_kg_k).next_up()
    }
    pub const fn interval_evaluations(&self) -> usize {
        self.interval_evaluations
    }
    /// Original-source changes, with no repeat of interval preparation.
    pub fn difference(
        &self,
        from: Temperature,
        to: Temperature,
    ) -> Result<DenseIsochoricDifference, If97DenseError> {
        for t in [from, to] {
            if !t.value().is_finite()
                || t.value() < self.lower_temperature_k
                || t.value() > self.upper_temperature_k
            {
                return Err(If97DenseError::InvalidInput(
                    "temperature outside prepared isochore",
                ));
            }
            If97DenseState::at_density_temperature(self.density(), t)?;
        }
        difference(self.density_kg_m3, from.value(), to.value())
    }
}

/// Stable changes along a single-phase, fixed-density original source path.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DenseIsochoricDifference {
    internal_energy_j_kg: f64,
    enthalpy_j_kg: f64,
    entropy_j_kg_k: f64,
    pressure_pa: f64,
}
impl DenseIsochoricDifference {
    pub const fn specific_internal_energy_change(&self) -> MassSpecificEnergyQuantityV1 {
        MassSpecificEnergyQuantityV1::new(self.internal_energy_j_kg)
    }
    pub const fn specific_enthalpy_change(&self) -> MassSpecificEnergyQuantityV1 {
        MassSpecificEnergyQuantityV1::new(self.enthalpy_j_kg)
    }
    pub const fn specific_entropy_change(&self) -> MassSpecificThermalQuantityV1 {
        MassSpecificThermalQuantityV1::new(self.entropy_j_kg_k)
    }
    pub const fn pressure_change(&self) -> Pressure {
        Pressure::new(self.pressure_pa)
    }
}

/// Convenience one-shot query with a fixed 4095-interval preparation ceiling.
///
/// Repeated or cancellation-aware consumers should prepare an interval once
/// with explicit controls and call its `difference` method. Unresolved path
/// admission refuses; this helper never joins valid endpoints across a gap.
pub fn isochoric_difference(
    density: Density,
    from: Temperature,
    to: Temperature,
) -> Result<DenseIsochoricDifference, If97DenseError> {
    let (lo, hi) = if from.value() <= to.value() {
        (from, to)
    } else {
        (to, from)
    };
    let interval = DenseIsochoricInterval::new(
        density,
        lo,
        hi,
        DenseIsochoricIntervalSolve {
            maximum_intervals: 4095,
        },
        || false,
    )
    .map_err(|error| match error {
        DenseIsochoricIntervalError::Source(source) => source,
        _ => If97DenseError::Unrepresentable("fixed-density source path not proven"),
    })?;
    interval.difference(from, to)
}

fn difference(rho: f64, from: f64, to: f64) -> Result<DenseIsochoricDifference, If97DenseError> {
    if from == to {
        return Ok(DenseIsochoricDifference {
            internal_energy_j_kg: 0.0,
            enthalpy_j_kg: 0.0,
            entropy_j_kg_k: 0.0,
            pressure_pa: 0.0,
        });
    }
    let (a, b, sign) = if from < to {
        (from, to, 1.0)
    } else {
        (to, from, -1.0)
    };
    let dt = b - a;
    let tau0 = CRITICAL_TEMPERATURE_K / a;
    let tau1 = CRITICAL_TEMPERATURE_K / b;
    let delta_tau = -tau0 * (dt / b);
    let mut delta_power = [0.0; 27];
    let mut original_power = 1.0;
    for j in 1..delta_power.len() {
        delta_power[j] = tau1 * delta_power[j - 1] + delta_tau * original_power;
        original_power *= tau0;
    }
    let delta = rho / CRITICAL_DENSITY_KG_M3;
    let mut density_power = [1.0; 12];
    for i in 1..density_power.len() {
        density_power[i] = density_power[i - 1] * delta;
    }
    let (mut u, mut s, mut p) = (Sum::default(), Sum::default(), Sum::default());
    p.add(LOG_COEFFICIENT * dt);
    for &(i, j, n) in &TERMS {
        let coefficient = n * density_power[i];
        if j != 0 {
            u.add(coefficient * j as f64 * delta_power[j - 1]);
        }
        s.add(coefficient * (j as f64 - 1.0) * delta_power[j]);
        if i != 0 {
            let product_change = if j == 0 {
                dt
            } else {
                CRITICAL_TEMPERATURE_K * delta_power[j - 1]
            };
            p.add(coefficient * i as f64 * product_change);
        }
    }
    let du = sign * R_SPECIFIC_J_KG_K * CRITICAL_TEMPERATURE_K * u.value;
    let ds = sign * R_SPECIFIC_J_KG_K * s.value;
    let dp = sign * rho * R_SPECIFIC_J_KG_K * p.value;
    let dh = du + dp / rho;
    if ![du, ds, dp, dh].into_iter().all(f64::is_finite) {
        return Err(If97DenseError::Unrepresentable(
            "fixed-density source increment",
        ));
    }
    Ok(DenseIsochoricDifference {
        internal_energy_j_kg: du,
        enthalpy_j_kg: dh,
        entropy_j_kg_k: ds,
        pressure_pa: dp,
    })
}

#[derive(Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) {
        crate::entropy::accumulate(&mut self.value, &mut self.correction, value);
    }
}

#[derive(Clone, Copy)]
struct Interval {
    lo: f64,
    hi: f64,
}
impl Interval {
    const fn point(x: f64) -> Self {
        Self { lo: x, hi: x }
    }
    fn add(self, b: Self) -> Self {
        Self {
            lo: (self.lo + b.lo).next_down(),
            hi: (self.hi + b.hi).next_up(),
        }
    }
    fn neg(self) -> Self {
        Self {
            lo: -self.hi,
            hi: -self.lo,
        }
    }
    fn sub(self, b: Self) -> Self {
        self.add(b.neg())
    }
    fn mul(self, b: Self) -> Self {
        let products = [
            self.lo * b.lo,
            self.lo * b.hi,
            self.hi * b.lo,
            self.hi * b.hi,
        ];
        if products.into_iter().any(f64::is_nan) {
            return Self::whole();
        }
        Self {
            lo: products
                .into_iter()
                .fold(f64::INFINITY, f64::min)
                .next_down(),
            hi: products
                .into_iter()
                .fold(f64::NEG_INFINITY, f64::max)
                .next_up(),
        }
    }
    fn div(self, b: Self) -> Self {
        if b.lo <= 0.0 && b.hi >= 0.0 {
            return Self::whole();
        }
        self.mul(Self {
            lo: (1.0 / b.hi).next_down(),
            hi: (1.0 / b.lo).next_up(),
        })
    }
    fn sqrt(self) -> Self {
        if self.lo < 0.0 {
            return Self::whole();
        }
        Self {
            lo: fs_math::det::sqrt(self.lo).next_down(),
            hi: fs_math::det::sqrt(self.hi).next_up(),
        }
    }
    fn scale(self, factor: f64) -> Self {
        self.mul(Self::point(factor))
    }
    const fn whole() -> Self {
        Self {
            lo: f64::NEG_INFINITY,
            hi: f64::INFINITY,
        }
    }
}
struct PathBounds {
    admitted: bool,
    cv: Interval,
}
// Density is fixed throughout preparation. Group equal temperature powers
// once, with outward rounding, before Horner evaluation on each subinterval.
// This resolves coefficient cancellation before introducing a temperature
// interval and avoids rebuilding all density powers for every subdivision.
struct IsochoricPolynomial {
    pressure: [Interval; 27],
    stiffness: [Interval; 27],
    cv: [Interval; 27],
}
impl IsochoricPolynomial {
    fn new(rho: f64) -> Self {
        let mut result = Self {
            pressure: [Interval::point(0.0); 27],
            stiffness: [Interval::point(0.0); 27],
            cv: [Interval::point(0.0); 27],
        };
        result.pressure[0] = Interval::point(LOG_COEFFICIENT);
        result.stiffness[0] = Interval::point(LOG_COEFFICIENT);
        let delta = Interval::point(rho).div(Interval::point(CRITICAL_DENSITY_KG_M3));
        let mut rho_powers = [Interval::point(1.0); 12];
        for i in 1..rho_powers.len() {
            rho_powers[i] = rho_powers[i - 1].mul(delta);
        }
        for &(i, j, n) in &TERMS {
            let coefficient = rho_powers[i].scale(n);
            if i != 0 {
                result.pressure[j] = result.pressure[j].add(coefficient.scale(i as f64));
                result.stiffness[j] =
                    result.stiffness[j].add(coefficient.scale((i * (i + 1)) as f64));
            }
            if j > 1 {
                result.cv[j] = result.cv[j].add(coefficient.scale(-((j * (j - 1)) as f64)));
            }
        }
        result
    }
    fn evaluate(coefficients: &[Interval; 27], tau: Interval) -> Interval {
        coefficients
            .iter()
            .rev()
            .fold(Interval::point(0.0), |sum, coefficient| {
                sum.mul(tau).add(*coefficient)
            })
    }
}
fn path_bounds(rho: f64, lo: f64, hi: f64, polynomial: &IsochoricPolynomial) -> PathBounds {
    let t = Interval { lo, hi };
    let tau = Interval::point(CRITICAL_TEMPERATURE_K).div(t);
    let p = t
        .mul(IsochoricPolynomial::evaluate(&polynomial.pressure, tau))
        .scale(R_SPECIFIC_J_KG_K)
        .scale(rho);
    let stiffness = t
        .mul(IsochoricPolynomial::evaluate(&polynomial.stiffness, tau))
        .scale(R_SPECIFIC_J_KG_K);
    let cv = IsochoricPolynomial::evaluate(&polynomial.cv, tau).scale(R_SPECIFIC_J_KG_K);
    let b23 = super::super::B23_COEFFICIENTS;
    let p23 = Interval::point(b23[0])
        .sub(t.scale(b23[1]))
        .add(t.scale(b23[2]).mul(t))
        .scale(1e6);
    let mut admitted = cv.lo > 0.0
        && cv.hi.is_finite()
        && stiffness.lo > 0.0
        && stiffness.hi.is_finite()
        && p.lo >= p23.hi
        && p.hi <= MAXIMUM_PRESSURE_PA;
    if lo < CRITICAL_TEMPERATURE_K {
        let ps = saturation_interval(t);
        admitted &= if rho > CRITICAL_DENSITY_KG_M3 {
            p.lo >= ps.hi
        } else if rho < CRITICAL_DENSITY_KG_M3 {
            p.hi <= ps.lo
        } else {
            false
        };
    }
    PathBounds { admitted, cv }
}
fn saturation_interval(t: Interval) -> Interval {
    let n = super::super::R4;
    let theta = t.add(Interval::point(n[8]).div(t.sub(Interval::point(n[9]))));
    let theta2 = theta.mul(theta);
    let a = theta2.add(theta.scale(n[0])).add(Interval::point(n[1]));
    let b = theta2
        .scale(n[2])
        .add(theta.scale(n[3]))
        .add(Interval::point(n[4]));
    let c = theta2
        .scale(n[5])
        .add(theta.scale(n[6]))
        .add(Interval::point(n[7]));
    let beta = c
        .scale(2.0)
        .div(b.neg().add(b.mul(b).sub(a.mul(c).scale(4.0)).sqrt()));
    let square = beta.mul(beta);
    square.mul(square).scale(1e6)
}

#[cfg(test)]
mod tests;
