//! Pull each scalar contact discrete gradient back through its own collocation.
//! A global energy correction along the full state increment can add net
//! internal force and alter kinetic efforts even for an equal/opposite row.
//! This owner keeps the contact effort in the span of the original rows.

use super::OpeningContactStep;
use crate::{ContactStorage, Obstacle};
use fs_math::det;
use fs_phs::ConstitutiveDiscreteGradient;

impl ContactStorage {
    pub(crate) fn has_contact_discrete_gradient(&self) -> bool {
        if self.inner.discrete_gradient_rule().is_some() {
            return true;
        }
        let Some(q) = self.inner.exact_quadratic_matrix() else {
            return false;
        };
        let n = q.len().isqrt();
        // A merely tolerated asymmetric quadratic needs its original rule.
        // Inspect the current declaration: an inner material owner may refresh.
        n > 0
            && n.checked_mul(n) == Some(q.len())
            && self
                .n_modes
                .checked_mul(2)
                .is_some_and(|prefix| prefix <= n)
            && q.iter()
                .enumerate()
                .all(|(index, value)| value.is_finite() && *value == q[(index % n) * n + index / n])
    }

    fn discrete_gradient_shape(&self, n: usize) -> bool {
        self.n_modes
            .checked_mul(2)
            .is_some_and(|prefix| prefix <= n)
            && self.obstacles.iter().all(|ob| {
                ob.n_points.checked_mul(self.n_modes) == Some(ob.collocation.len())
                    && ob.gaps.len() == ob.n_points
                    && ob.weights.len() == ob.n_points
            })
    }
}

impl ConstitutiveDiscreteGradient for ContactStorage {
    fn evaluate(&self, a: &[f64], b: &[f64], midpoint: &[f64], out: &mut [f64]) {
        let n = a.len();
        if b.len() != n
            || midpoint.len() != n
            || out.len() != n
            || !self.discrete_gradient_shape(n)
            || a.iter().chain(b).chain(midpoint).any(|v| !v.is_finite())
        {
            out.fill(f64::NAN);
            return;
        }
        if let Some(rule) = self.inner.discrete_gradient_rule() {
            rule.evaluate(a, b, midpoint, out);
        } else if self.inner.exact_quadratic_matrix().is_some_and(|q| {
            n.checked_mul(n) == Some(q.len()) && self.has_contact_discrete_gradient()
        }) {
            self.inner.gradient(midpoint, out);
        } else {
            out.fill(f64::NAN);
            return;
        }
        for ob in &self.obstacles {
            for i in 0..ob.n_points {
                let (point, after) = point_step(ob, i, self.n_modes, a, b);
                let Ok((elastic, _)) = point.coefficients(after) else {
                    out.fill(f64::NAN);
                    return;
                };
                for k in 0..self.n_modes {
                    out[2 * k] += elastic * ob.collocation[i * self.n_modes + k];
                }
            }
        }
    }

    fn end_jacobian_action(
        &self,
        a: &[f64],
        b: &[f64],
        direction: &[f64],
        out: &mut [f64],
    ) -> bool {
        let n = a.len();
        if b.len() != n
            || direction.len() != n
            || out.len() != n
            || !self.discrete_gradient_shape(n)
            || a.iter().chain(b).chain(direction).any(|v| !v.is_finite())
        {
            return false;
        }
        if let Some(rule) = self.inner.discrete_gradient_rule() {
            if !rule.end_jacobian_action(a, b, direction, out) {
                return false;
            }
        } else {
            let Some(q) = self.inner.exact_quadratic_matrix() else {
                return false;
            };
            if n.checked_mul(n) != Some(q.len()) || !self.has_contact_discrete_gradient() {
                return false;
            }
            for i in 0..n {
                out[i] = (0..n).map(|j| 0.5 * q[i * n + j] * direction[j]).sum();
            }
        }
        for ob in &self.obstacles {
            for i in 0..ob.n_points {
                let row = &ob.collocation[i * self.n_modes..(i + 1) * self.n_modes];
                let rate: f64 = row
                    .iter()
                    .enumerate()
                    .map(|(k, v)| v * direction[2 * k])
                    .sum();
                if rate == 0.0 {
                    continue;
                }
                let (point, after) = point_step(ob, i, self.n_modes, a, b);
                let Some(slope) = endpoint_closure_slope(&point, after) else {
                    return false;
                };
                for (k, coefficient) in row.iter().enumerate() {
                    out[2 * k] += slope * rate * coefficient;
                }
            }
        }
        out.iter().all(|value| value.is_finite())
    }
}

fn point_step(
    ob: &Obstacle,
    i: usize,
    n_modes: usize,
    a: &[f64],
    b: &[f64],
) -> (OpeningContactStep, f64) {
    let row = &ob.collocation[i * n_modes..(i + 1) * n_modes];
    let mut before = 0.0;
    let mut after = 0.0;
    for (k, coefficient) in row.iter().enumerate() {
        before -= coefficient * a[2 * k];
        after -= coefficient * b[2 * k];
    }
    (
        OpeningContactStep {
            before,
            gap: ob.gaps[i],
            stiffness: ob.stiffness,
            weight: ob.weights[i],
            alpha: ob.alpha,
            // Loss is a separate port; it must not alter or refuse the
            // conservative effort through an overflowing damping coefficient.
            loss: 0.0,
        },
        after,
    )
}

// d[ (V(b)-V(a))/(b-a) ]/db in the signed CLOSURE coordinate.
// At coincident zero penetration alpha=1 has a kink; choose its inactive
// one-sided Newton representative, as the original contact Hessian does.
fn endpoint_closure_slope(point: &OpeningContactStep, after: f64) -> Option<f64> {
    let (elastic, _) = point.coefficients(after).ok()?;
    let a = -point.before - point.gap;
    let b = -after - point.gap;
    let delta = point.before - after;
    if (a <= 0.0 && b <= 0.0) || point.weight == 0.0 || point.stiffness == 0.0 {
        return Some(0.0);
    }
    let center = f64::midpoint(a, b);
    let z = (0.5 * delta) / center;
    let slope = if a > 0.0 && b > 0.0 && z.abs() <= 0.25 && (point.alpha * z).abs() <= 0.5 {
        // Integral_0^1 t*V''(a+t*(b-a)) dt. This centered binomial series
        // avoids subtracting two nearly equal forces near coincidence.
        // Under this branch condition terms after the first contract by <=1/4;
        // 33 terms leave the truncation below f64 rounding.
        let mut term = 1.0;
        let mut sum = 0.5;
        for j in 1..=32 {
            term *= ((point.alpha - j as f64) / j as f64) * z;
            let denominator = (j + 1 + j % 2) as f64;
            sum += 0.5 * term / denominator;
        }
        point.weight * point.stiffness * point.alpha * det::pow(center, point.alpha - 1.0) * sum
    } else if delta != 0.0 {
        let force_after = point.weight * point.stiffness * det::pow(b.max(0.0), point.alpha);
        (force_after - elastic) / delta
    } else {
        0.5 * point.weight * point.stiffness * point.alpha * det::pow(b.max(0.0), point.alpha - 1.0)
    };
    (slope.is_finite() && slope >= 0.0).then_some(slope)
}
