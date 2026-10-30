//! Whole-interval isobaric heat-capacity bounds for the original ice Ih EOS.
//!
//! The Gibbs second temperature derivative reduces to rational functions of
//! real temperature. Outward-rounded arithmetic encloses those functions and
//! the original pressure-dependent coefficients; no sampled Cp minimum or
//! fitted heat-capacity floor is used. The result bounds the retained source
//! equation, not its experimental uncertainty or the stability of other ice
//! polymorphs. Phase admission remains the caller's responsibility.

use super::{
    IceIhError, IceIhState, NORMAL_PRESSURE_PA, R1, R2_COEFFICIENTS, T_POLES,
    TRIPLE_POINT_PRESSURE_PA, TRIPLE_POINT_TEMPERATURE_K,
};
use core::fmt;
use fs_qty::{Pressure, Temperature};

/// Hard ceiling in addition to the caller's work allowance.
pub const MAX_BOUND_INTERVALS: usize = 65_535;
const MAX_DEPTH: u8 = 24;

/// Work allowed to enclose an immutable pressure/temperature interval.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct IceIhHeatCapacityBoundSolve {
    /// Includes coarse intervals that require further subdivision.
    pub maximum_intervals: usize,
}

/// Source admission, bounded numerical work, or cancellation refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum IceIhHeatCapacityBoundError {
    Source(IceIhError),
    Invalid(&'static str),
    NotProven {
        interval_evaluations: usize,
        lower_temperature_k: f64,
        upper_temperature_k: f64,
    },
    Cancelled,
}
impl fmt::Display for IceIhHeatCapacityBoundError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "ice Ih heat-capacity bound: {self:?}")
    }
}
impl std::error::Error for IceIhHeatCapacityBoundError {}

/// Positive Cp enclosure over every temperature in the prepared interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IceIhHeatCapacityBound {
    minimum_j_kg_k: f64,
    maximum_j_kg_k: f64,
    interval_evaluations: usize,
}
impl IceIhHeatCapacityBound {
    pub const fn minimum_heat_capacity_j_kg_k(&self) -> f64 {
        self.minimum_j_kg_k
    }
    pub const fn maximum_heat_capacity_j_kg_k(&self) -> f64 {
        self.maximum_j_kg_k
    }
    /// Conservative upper bound on dT/dh [kg K/J].
    pub fn maximum_temperature_slope(&self) -> f64 {
        (1.0 / self.minimum_j_kg_k).next_up()
    }
    pub const fn interval_evaluations(&self) -> usize {
        self.interval_evaluations
    }
}

/// Bound original ice Ih Cp and the inverse caloric slope on one isobar.
///
/// All arithmetic in the enclosure is rounded outward. A local bound must
/// have upper/lower <= 2 before acceptance, preventing an almost-zero coarse
/// lower bound from imposing an unnecessarily tiny thermal timestep. Failure
/// to prove a bound within the caller's allowance returns no partial result.
pub fn isobaric_heat_capacity_bound(
    pressure: Pressure,
    lower: Temperature,
    upper: Temperature,
    controls: IceIhHeatCapacityBoundSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<IceIhHeatCapacityBound, IceIhHeatCapacityBoundError> {
    use IceIhHeatCapacityBoundError as Error;
    if cancelled() {
        return Err(Error::Cancelled);
    }
    let (lo, hi) = (lower.value(), upper.value());
    if !lo.is_finite() || !hi.is_finite() || lo > hi {
        return Err(Error::Invalid("finite ordered temperatures required"));
    }
    if !(1..=MAX_BOUND_INTERVALS).contains(&controls.maximum_intervals) {
        return Err(Error::Invalid("invalid interval evaluation allowance"));
    }
    IceIhState::at_pressure_temperature(pressure, lower).map_err(Error::Source)?;
    IceIhState::at_pressure_temperature(pressure, upper).map_err(Error::Source)?;
    let coefficients = coefficient_intervals(pressure.value());
    let mut stack = [(0.0, 0.0, 0_u8); MAX_DEPTH as usize + 2];
    stack[0] = (lo, hi, 0);
    let (mut pending, mut evaluations) = (1, 0);
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
        if evaluations == controls.maximum_intervals {
            return Err(unresolved(evaluations));
        }
        evaluations += 1;
        let cp = heat_capacity_interval(coefficients, lo, hi);
        if cp.lo > 0.0 && cp.hi.is_finite() && cp.hi <= 2.0 * cp.lo {
            minimum = minimum.min(cp.lo);
            maximum = maximum.max(cp.hi);
            continue;
        }
        let mid = lo + 0.5 * (hi - lo);
        if depth == MAX_DEPTH || mid <= lo || mid >= hi {
            return Err(unresolved(evaluations));
        }
        stack[pending] = (mid, hi, depth + 1);
        stack[pending + 1] = (lo, mid, depth + 1);
        pending += 2;
    }
    if cancelled() {
        return Err(Error::Cancelled);
    }
    Ok(IceIhHeatCapacityBound {
        minimum_j_kg_k: minimum,
        maximum_j_kg_k: maximum,
        interval_evaluations: evaluations,
    })
}

#[derive(Clone, Copy, Debug)]
struct Interval {
    lo: f64,
    hi: f64,
}
impl Interval {
    const fn point(value: f64) -> Self {
        Self {
            lo: value,
            hi: value,
        }
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
    fn square(self) -> Self {
        let (a, b) = (self.lo * self.lo, self.hi * self.hi);
        Self {
            lo: if self.lo <= 0.0 && self.hi >= 0.0 {
                0.0
            } else {
                a.min(b).next_down().max(0.0)
            },
            hi: a.max(b).next_up(),
        }
    }
    fn reciprocal_positive(self) -> Self {
        if self.lo <= 0.0 {
            return Self {
                lo: f64::NEG_INFINITY,
                hi: f64::INFINITY,
            };
        }
        Self {
            lo: (1.0 / self.hi).next_down().max(0.0),
            hi: (1.0 / self.lo).next_up(),
        }
    }
}

fn coefficient_intervals(p: f64) -> [[Interval; 2]; 2] {
    let q = Interval::point(p)
        .subtract(Interval::point(NORMAL_PRESSURE_PA))
        .multiply(Interval::point(TRIPLE_POINT_PRESSURE_PA).reciprocal_positive());
    let component = |imaginary: bool| {
        let get = |i: usize| {
            Interval::point(if imaginary {
                R2_COEFFICIENTS[i].im
            } else {
                R2_COEFFICIENTS[i].re
            })
        };
        get(2).multiply(q).add(get(1)).multiply(q).add(get(0))
    };
    [
        [Interval::point(R1.re), Interval::point(R1.im)],
        [component(false), component(true)],
    ]
}

fn heat_capacity_interval(coefficients: [[Interval; 2]; 2], lo: f64, hi: f64) -> Interval {
    let tau = Interval { lo, hi }
        .multiply(Interval::point(TRIPLE_POINT_TEMPERATURE_K).reciprocal_positive());
    let tau2 = tau.square();
    let minus_twice_tau3 = tau2.multiply(tau).multiply(Interval::point(-2.0));
    let mut cp = Interval::point(0.0);
    for (t, [real, imaginary]) in T_POLES.into_iter().zip(coefficients) {
        let (a, b) = (Interval::point(t.re), Interval::point(t.im));
        let (aa, bb) = (a.square(), b.square());
        // Cp = -2*tau^3 Re[r / (t*(t^2-tau^2))]. Keeping the
        // denominator factored gives an explicitly positive constant term.
        let real_denominator = a.multiply(
            aa.subtract(bb.multiply(Interval::point(3.0)))
                .subtract(tau2),
        );
        let imaginary_denominator = b.multiply(
            aa.multiply(Interval::point(3.0))
                .subtract(bb)
                .subtract(tau2),
        );
        let denominator = aa.add(bb).multiply(
            aa.subtract(bb)
                .subtract(tau2)
                .square()
                .add(aa.multiply(bb).multiply(Interval::point(4.0))),
        );
        let numerator = real
            .multiply(real_denominator)
            .add(imaginary.multiply(imaginary_denominator));
        cp = cp.add(
            minus_twice_tau3
                .multiply(numerator)
                .multiply(denominator.reciprocal_positive()),
        );
    }
    cp
}

#[cfg(test)]
mod tests;
