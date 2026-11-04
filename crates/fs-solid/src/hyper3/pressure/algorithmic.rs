//! Discrete gradient of the native clipped tetrahedral pressure mean.
//!
//! Averaging coordinate increments over all 24 orders is permutation symmetric
//! and time reversible. Its nonnegative weights close the exact endpoint mean
//! change. Exact divided-difference arithmetic avoids subtracting nearly equal
//! potentials, and has the analytic derivative as its coincident-endpoint limit.
//! Only partially immersed elements need this bounded four-coordinate work.
//! If an endpoint wet edge collapses onto an already coplanar vertex, its
//! individual rational factors have a removable 0/0. Use the complete moment
//! quotient there; only a cancellation-limited numerator selects the analytic
//! midpoint derivative. The native time owner still gates measured work.

use core::ops::{Add, Div, Mul, Neg, Sub};
use fs_ad::{Dual, Dual64, Real};

use super::free_surface::{MomentScalar, tet_mean, tet_mean_scalar};
use super::{HyperTetError, finite};

pub(in crate::hyper3) struct PressureMean {
    pub before_pa: f64,
    pub after_pa: f64,
    pub after_derivative: [f64; 4],
    pub weights: [f64; 4],
    pub tangent_after: [[f64; 4]; 4],
}

/// Endpoint values and an exact divided difference with respect to one nodal
/// pressure. It can itself carry native Duals for the endpoint force tangent.
#[derive(Clone, Copy)]
struct Secant<R> {
    before: R,
    after: R,
    slope: R,
}

impl<R: Real> Secant<R> {
    fn fixed(value: R) -> Self {
        Self {
            before: value,
            after: value,
            slope: R::zero(),
        }
    }
    fn mean(self) -> R {
        self.before * R::from_f64(0.5) + self.after * R::from_f64(0.5)
    }
}

impl<R: Real> MomentScalar for Secant<R> {
    fn constant(value: f64) -> Self {
        Self::fixed(R::from_f64(value))
    }
    fn value(self) -> f64 {
        self.mean().value()
    }
}
impl<R: Real> Add for Secant<R> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self {
            before: self.before + rhs.before,
            after: self.after + rhs.after,
            slope: self.slope + rhs.slope,
        }
    }
}
impl<R: Real> Neg for Secant<R> {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            before: -self.before,
            after: -self.after,
            slope: -self.slope,
        }
    }
}
impl<R: Real> Sub for Secant<R> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        self + -rhs
    }
}
impl<R: Real> Mul for Secant<R> {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        Self {
            before: self.before * rhs.before,
            after: self.after * rhs.after,
            slope: self.mean() * rhs.slope + rhs.mean() * self.slope,
        }
    }
}
impl<R: Real> Div for Secant<R> {
    type Output = Self;
    fn div(self, rhs: Self) -> Self {
        let before = self.before / rhs.before;
        let after = self.after / rhs.after;
        let quotient_mean = before * R::from_f64(0.5) + after * R::from_f64(0.5);
        Self {
            before,
            after,
            slope: (self.slope - quotient_mean * rhs.slope) / rhs.mean(),
        }
    }
}

fn segment<R: Real>(fixed: [R; 4], node: usize, before: R, after: R) -> R {
    if (before.value() == 0.0 || after.value() == 0.0)
        && before.value().max(after.value()) > 0.0
        && fixed
            .iter()
            .enumerate()
            .any(|(i, p)| i != node && p.value() == 0.0)
    {
        // A vanishing wet edge can have undefined individual fractions while
        // its COMPLETE moment and derivatives are finite. Keep the derivative
        // of the other coplanar coordinate; making that edge constant loses it.
        let mut initial = fixed;
        let mut final_values = fixed;
        initial[node] = before;
        final_values[node] = after;
        let first = tet_mean_scalar(initial);
        let last = tet_mean_scalar(final_values);
        let change = last - first;
        if change.value().abs() > 64.0 * f64::EPSILON * first.value().abs().max(last.value().abs())
        {
            return change / (after - before);
        }
        // Decide from cancellation in the moment, NOT the magnitude of dry
        // nodal pressures: arbitrarily thin wet caps retain resolved quotients.
        let mut midpoint = fixed.map(Dual::<R, 1>::constant);
        midpoint[node] = Dual::variable(before * R::from_f64(0.5) + after * R::from_f64(0.5), 0);
        return tet_mean_scalar(midpoint).eps[0];
    }
    let mut values = fixed.map(Secant::fixed);
    values[node] = Secant {
        before,
        after,
        slope: R::one(),
    };
    tet_mean_scalar(values).slope
}

fn coordinate_secant<R: Real>(fixed: [R; 4], node: usize, before: R, after: R) -> R {
    let (a, b) = (before.value(), after.value());
    if (a < 0.0 && b > 0.0) || (a > 0.0 && b < 0.0) {
        // The same positive decomposition applies on each side of the zero
        // crossing. Combine their divided differences by interval lengths,
        // scaled before adding opposite magnitudes (no overflowing b-a).
        let scale = R::from_f64(a.abs().max(b.abs()));
        let left = before.abs() / scale;
        let right = after.abs() / scale;
        let sum = left + right;
        (left / sum) * segment(fixed, node, before, R::zero())
            + (right / sum) * segment(fixed, node, R::zero(), after)
    } else {
        segment(fixed, node, before, after)
    }
}

fn discrete_gradient<R: Real>(before: [R; 4], after: [R; 4]) -> [R; 4] {
    if before.iter().chain(&after).all(|p| p.value() <= 0.0) {
        return [R::zero(); 4];
    }
    if before.iter().chain(&after).all(|p| p.value() >= 0.0) {
        return [R::from_f64(0.25); 4];
    }
    core::array::from_fn(|node| {
        let mut weight = R::zero();
        for subset in 0_u32..16 {
            if subset & (1 << node) != 0 {
                continue;
            }
            let fixed = core::array::from_fn(|i| {
                if subset & (1 << i) == 0 {
                    before[i]
                } else {
                    after[i]
                }
            });
            // |S|! (3-|S|)! / 4! counts all paths whose earlier coordinates
            // are exactly S. Thirty-two edges replace 24 full permutations.
            let multiplicity = match subset.count_ones() {
                0 | 3 => 0.25,
                _ => 1.0 / 12.0,
            };
            weight = weight
                + R::from_f64(multiplicity)
                    * coordinate_secant(fixed, node, before[node], after[node]);
        }
        weight
    })
}

pub(in crate::hyper3) fn evaluate(
    before: [f64; 4],
    after: [f64; 4],
    clipped: bool,
    tangent: bool,
) -> Result<PressureMean, HyperTetError> {
    if !clipped {
        return Ok(PressureMean {
            before_pa: finite(before.into_iter().map(|p| p / 4.0).sum())?,
            after_pa: finite(after.into_iter().map(|p| p / 4.0).sum())?,
            after_derivative: [0.25; 4],
            weights: [0.25; 4],
            tangent_after: [[0.0; 4]; 4],
        });
    }
    let mut result = PressureMean {
        before_pa: tet_mean(before)?,
        after_pa: tet_mean(after)?,
        after_derivative: [0.0; 4],
        weights: [0.0; 4],
        tangent_after: [[0.0; 4]; 4],
    };
    if tangent {
        let seeded = core::array::from_fn(|i| Dual64::<4>::variable(after[i], i));
        result.after_derivative = tet_mean_scalar(seeded).eps;
        let weights = discrete_gradient(before.map(Dual64::constant), seeded);
        for (i, value) in weights.into_iter().enumerate() {
            result.weights[i] = value.re;
            result.tangent_after[i] = value.eps;
        }
    } else {
        result.weights = discrete_gradient(before, after);
    }
    for value in result
        .weights
        .iter()
        .chain(&result.after_derivative)
        .chain(result.tangent_after.iter().flatten())
    {
        finite(*value)?;
    }
    Ok(result)
}

#[cfg(test)]
#[path = "algorithmic/tests.rs"]
mod tests;
