//! Conservative isobaric heat-capacity bounds for the original IF97 equations.
//!
//! Every arithmetic operation in the Gibbs-monomial enclosure is rounded
//! outward. Adaptive subdivision resolves dependency between temperature
//! powers; positive samples are never used as proof. These bounds concern the
//! equations with the retained binary64 source coefficients, not uncertainty
//! in the empirical water formulation.

use super::{
    If97Error, If97Region, R_SPECIFIC_J_KG_K, R1_PRESSURE_STAR_PA, R1_TEMPERATURE_STAR_K, R1_TERMS,
    R2_IDEAL, R2_PRESSURE_STAR_PA, R2_RESIDUAL, R2_TEMPERATURE_STAR_K,
};
use core::fmt;
use fs_qty::{Pressure, Temperature};

/// Hard ceiling in addition to the caller's explicit work budget.
pub const MAX_BOUND_INTERVALS: usize = 65_535;
const MAX_DEPTH: u8 = 24;

/// Bounded work for a proof over a complete, immutable isobaric source interval.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct If97HeatCapacityBoundSolve {
    /// Maximum Gibbs-interval evaluations, including intervals later subdivided.
    pub maximum_intervals: usize,
}

/// Source admission, unresolved enclosure, or transactional cancellation.
#[derive(Clone, Debug, PartialEq)]
pub enum If97HeatCapacityBoundError {
    /// The requested source path is not admitted.
    Source(If97Error),
    /// Invalid interval order or numerical work allowance.
    Invalid(&'static str),
    /// An interval remains unresolved when the permitted work is exhausted.
    NotProven {
        interval_evaluations: usize,
        lower_temperature_k: f64,
        upper_temperature_k: f64,
    },
    /// No partial bound is published after cancellation.
    Cancelled,
}
impl fmt::Display for If97HeatCapacityBoundError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "IF97 heat-capacity bound: {self:?}")
    }
}
impl std::error::Error for If97HeatCapacityBoundError {}

/// Proven bounds over every temperature in the requested source interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct If97HeatCapacityBound {
    minimum_j_kg_k: f64,
    maximum_j_kg_k: f64,
    interval_evaluations: usize,
}
impl If97HeatCapacityBound {
    /// Strictly positive lower bound for Cp [J/(kg K)].
    pub const fn minimum_heat_capacity_j_kg_k(&self) -> f64 {
        self.minimum_j_kg_k
    }
    /// Finite upper bound for Cp [J/(kg K)].
    pub const fn maximum_heat_capacity_j_kg_k(&self) -> f64 {
        self.maximum_j_kg_k
    }
    /// Conservative upper bound on dT/dh [kg K/J].
    pub fn maximum_temperature_slope(&self) -> f64 {
        (1.0 / self.minimum_j_kg_k).next_up()
    }
    /// Actual evaluations, including every rejected coarse enclosure.
    pub const fn interval_evaluations(&self) -> usize {
        self.interval_evaluations
    }
}

/// Prove Cp bounds on one original Region 1 or Region 2 isobaric path.
///
/// The returned lower bound also bounds the inverse caloric law globally.
/// A local enclosure must have upper/lower <= 2 before acceptance, so a nearly
/// zero positive coarse bound cannot make an otherwise useful thermal CFL
/// allowance arbitrarily pessimistic. No heat-capacity floor, interpolation,
/// source extrapolation, or modification of the original equation is applied.
pub fn isobaric_heat_capacity_bound(
    pressure: Pressure,
    region: If97Region,
    lower: Temperature,
    upper: Temperature,
    controls: If97HeatCapacityBoundSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<If97HeatCapacityBound, If97HeatCapacityBoundError> {
    use If97HeatCapacityBoundError as Error;
    if cancelled() {
        return Err(Error::Cancelled);
    }
    let lo = lower.value();
    let hi = upper.value();
    if !lo.is_finite() || !hi.is_finite() || lo > hi {
        return Err(Error::Invalid("temperatures must be finite and ordered"));
    }
    if controls.maximum_intervals == 0 || controls.maximum_intervals > MAX_BOUND_INTERVALS {
        return Err(Error::Invalid("invalid interval evaluation allowance"));
    }
    super::increments::admit_isobaric_interval(pressure.value(), lo, hi, region)
        .map_err(Error::Source)?;
    let mut stack = [(0.0, 0.0, 0_u8); MAX_DEPTH as usize + 2];
    stack[0] = (lo, hi, 0);
    let (mut pending, mut work) = (1, 0);
    let (mut minimum, mut maximum) = (f64::INFINITY, 0.0_f64);
    while pending != 0 {
        if cancelled() {
            return Err(Error::Cancelled);
        }
        pending -= 1;
        let (lo, hi, depth) = stack[pending];
        let unresolved = |work| Error::NotProven {
            interval_evaluations: work,
            lower_temperature_k: lo,
            upper_temperature_k: hi,
        };
        if work == controls.maximum_intervals {
            return Err(unresolved(work));
        }
        work += 1;
        let cp = heat_capacity_interval(pressure.value(), region, lo, hi);
        if cp.lo > 0.0 && cp.hi.is_finite() && cp.hi <= 2.0 * cp.lo {
            minimum = minimum.min(cp.lo);
            maximum = maximum.max(cp.hi);
            continue;
        }
        let mid = lo + 0.5 * (hi - lo);
        if depth == MAX_DEPTH || mid <= lo || mid >= hi {
            return Err(unresolved(work));
        }
        stack[pending] = (mid, hi, depth + 1);
        stack[pending + 1] = (lo, mid, depth + 1);
        pending += 2;
    }
    if cancelled() {
        return Err(Error::Cancelled);
    }
    Ok(If97HeatCapacityBound {
        minimum_j_kg_k: minimum,
        maximum_j_kg_k: maximum,
        interval_evaluations: work,
    })
}

// Enclosures treat the finite input/source constants as exact binary64 values.
// No hardware rounding-mode change or platform transcendental is involved.
#[derive(Clone, Copy, Debug)]
struct Interval {
    lo: f64,
    hi: f64,
}
impl Interval {
    const fn point(x: f64) -> Self {
        Self { lo: x, hi: x }
    }
    fn add(self, other: Self) -> Self {
        Self {
            lo: (self.lo + other.lo).next_down(),
            hi: (self.hi + other.hi).next_up(),
        }
    }
    fn negate(self) -> Self {
        Self {
            lo: -self.hi,
            hi: -self.lo,
        }
    }
    fn subtract(self, other: Self) -> Self {
        self.add(other.negate())
    }
    fn multiply(self, other: Self) -> Self {
        let products = [
            self.lo * other.lo,
            self.lo * other.hi,
            self.hi * other.lo,
            self.hi * other.hi,
        ];
        if products.into_iter().any(f64::is_nan) {
            return Self {
                lo: f64::NEG_INFINITY,
                hi: f64::INFINITY,
            };
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
    fn reciprocal_positive(self) -> Self {
        Self {
            lo: (1.0 / self.hi).next_down().max(0.0),
            hi: (1.0 / self.lo).next_up(),
        }
    }
    fn positive_power(mut self, power: i32) -> Self {
        if power < 0 {
            self = self.reciprocal_positive();
        }
        let mut result = Self::point(1.0);
        let mut n = power.unsigned_abs();
        while n != 0 {
            if n & 1 != 0 {
                result = result.multiply(self);
                result.lo = result.lo.max(0.0);
            }
            n >>= 1;
            if n != 0 {
                self = self.multiply(self);
                self.lo = self.lo.max(0.0);
            }
        }
        result
    }
}

fn heat_capacity_interval(p: f64, region: If97Region, lo: f64, hi: f64) -> Interval {
    let (p_star, t_star, shift, terms) = match region {
        If97Region::Region1 => (
            R1_PRESSURE_STAR_PA,
            R1_TEMPERATURE_STAR_K,
            1.222,
            &R1_TERMS[..],
        ),
        If97Region::Region2 => (
            R2_PRESSURE_STAR_PA,
            R2_TEMPERATURE_STAR_K,
            0.5,
            &R2_RESIDUAL[..],
        ),
    };
    let tau = Interval { lo, hi }
        .reciprocal_positive()
        .multiply(Interval::point(t_star));
    let y = tau.subtract(Interval::point(shift));
    let pi = Interval::point(p).multiply(Interval::point(p_star).reciprocal_positive());
    let x = match region {
        If97Region::Region1 => Interval::point(7.1).subtract(pi),
        If97Region::Region2 => pi,
    };
    let minus_r = Interval::point(-R_SPECIFIC_J_KG_K);
    let tau_squared = tau.positive_power(2);
    let mut cp = Interval::point(0.0);
    for term in terms {
        if term.j == 0 || term.j == 1 {
            continue;
        }
        let value = Interval::point(term.n)
            .multiply(Interval::point(f64::from(term.j)))
            .multiply(Interval::point(f64::from(term.j - 1)))
            .multiply(x.positive_power(term.i))
            .multiply(y.positive_power(term.j - 2))
            .multiply(tau_squared)
            .multiply(minus_r);
        cp = cp.add(value);
    }
    if region == If97Region::Region2 {
        for term in &R2_IDEAL {
            if term.j == 0 || term.j == 1 {
                continue;
            }
            // Here tau^2 cancels algebraically before enclosing the original
            // ideal-gas Gibbs second derivative, improving the global bound.
            let value = Interval::point(term.n)
                .multiply(Interval::point(f64::from(term.j)))
                .multiply(Interval::point(f64::from(term.j - 1)))
                .multiply(tau.positive_power(term.j))
                .multiply(minus_r);
            cp = cp.add(value);
        }
    }
    cp
}

#[cfg(test)]
mod tests;
