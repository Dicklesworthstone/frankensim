//! Original logarithmic extension, including its scalar component work law.

use super::{PhsError, RelaxationStrainDiscreteGradient, RelaxationStrainMap};
use crate::interconnection::accumulation::Sum;
use fs_math::det;

/// Scalar logarithmic strain `epsilon(x)=ln(reference_stretch + projection.x)`.
///
/// The stretch must be positive. The projection and reference are original
/// geometry, in units making stretch dimensionless; no coordinate, material,
/// equilibrium stress or applicable strain range is inferred. Zero projection
/// entries remain exactly uncoupled, including independent kinetic coordinates.
/// This is a scalar source map, not a tensor finite-strain material model.
/// Maxwell branch stiffness, clocks and history still belong to the composing
/// owner. All source actions are deterministic and allocate no warm scratch.
#[derive(Clone, Debug)]
pub struct LogarithmicStrainMap {
    reference_stretch: f64,
    projection: Box<[f64]>,
}

impl LogarithmicStrainMap {
    pub fn new(reference_stretch: f64, projection: Vec<f64>) -> Result<Self, PhsError> {
        if !reference_stretch.is_finite()
            || reference_stretch <= 0.0
            || projection.is_empty()
            || projection.iter().any(|v| !v.is_finite())
        {
            return Err(PhsError::Dimension {
                what: "finite positive reference stretch and complete projection required",
            });
        }
        Ok(Self {
            reference_stretch,
            projection: projection.into_boxed_slice(),
        })
    }

    fn stretch(&self, state: &[f64]) -> Option<f64> {
        if state.len() != self.projection.len() || state.iter().any(|v| !v.is_finite()) {
            return None;
        }
        let mut sum = Sum::dot(&self.projection, state)?;
        sum.product(1.0, self.reference_stretch)?;
        let stretch = sum.value()?;
        (stretch > 0.0).then_some(stretch)
    }

    fn scaled_increment(&self, a: &[f64], b: &[f64], scale: f64) -> Option<f64> {
        if a.len() != self.projection.len()
            || b.len() != a.len()
            || a.iter().chain(b).any(|v| !v.is_finite())
        {
            return None;
        }
        let mut sum = Sum::default();
        for i in 0..a.len() {
            sum.product(self.projection[i] / scale, b[i] - a[i])?;
        }
        sum.value()
    }
}

// L(a,b)=(ln b-ln a)/(b-a)=A/m. The normalized second-endpoint
// derivative is B=m^2*dL/db. For nearby stretches use the analytic even
// series atanh(r)/r, r=(b-a)/(2m), instead of subtracting close logarithms.
// At |r|<=1/8 the omitted A term is below 3e-18; the derivative term is
// below 4e-16. Outside this range the divided differences are well separated.
fn logarithmic_mean(a: f64, b: f64) -> Option<(f64, f64, f64)> {
    let m = f64::midpoint(a, b);
    let r = 0.5 * ((b - a) / m);
    let (mean, derivative) = if r.abs() <= 0.125 {
        let r2 = r * r;
        let (mut power, mut value, mut slope) = (1.0, 1.0, 0.0);
        for k in 1..=9 {
            slope += (2 * k) as f64 / (2 * k + 1) as f64 * r * power;
            power *= r2;
            value += power / (2 * k + 1) as f64;
        }
        (value, 0.5 * ((1.0 - r) * slope - value))
    } else {
        let ratio = b / a;
        let change = if ratio.is_finite() && ratio > 0.0 {
            det::ln(ratio)
        } else {
            det::ln(b) - det::ln(a)
        };
        let value = (m / (b - a)) * change;
        (value, (m / (b - a)) * (m / b - value))
    };
    [m, mean, derivative]
        .iter()
        .all(|v| v.is_finite())
        .then_some((m, mean, derivative))
}

// ln(1+t)-t is the curvature part of the original strain increment. Retain
// its leading quadratic/cubic work when |t| is small; never recover it from
// a subtraction of total stored energies. Outside the series interval the
// two terms are separated. This helper is specific to this constitutive law.
fn log_curvature(t: f64) -> Option<f64> {
    if !t.is_finite() || t <= -1.0 {
        return None;
    }
    let result = if t.abs() <= 0.125 {
        let mut power = t * t;
        let mut sum = -0.5 * power;
        for k in 3..=20 {
            power *= -t;
            sum -= power / k as f64;
        }
        sum
    } else {
        det::ln(1.0 + t) - t
    };
    result.is_finite().then_some(result)
}

impl RelaxationStrainMap for LogarithmicStrainMap {
    fn state_dimension(&self) -> usize {
        self.projection.len()
    }
    fn strain(&self, state: &[f64]) -> f64 {
        self.stretch(state).map_or(f64::NAN, det::ln)
    }
    fn gradient(&self, state: &[f64], out: &mut [f64]) {
        out.fill(f64::NAN);
        if out.len() != self.projection.len() {
            return;
        }
        let Some(stretch) = self.stretch(state) else {
            return;
        };
        for (o, &weight) in out.iter_mut().zip(&self.projection) {
            *o = weight / stretch;
        }
    }
    fn hessian_vector(&self, state: &[f64], direction: &[f64], out: &mut [f64]) {
        out.fill(f64::NAN);
        if out.len() != self.projection.len()
            || direction.len() != out.len()
            || direction.iter().any(|v| !v.is_finite())
        {
            return;
        }
        let Some(stretch) = self.stretch(state) else {
            return;
        };
        let mut change = Sum::default();
        for (&weight, &v) in self.projection.iter().zip(direction) {
            if change.product(weight / stretch, v).is_none() {
                return;
            }
        }
        let Some(change) = change.value() else {
            return;
        };
        for (o, &weight) in out.iter_mut().zip(&self.projection) {
            *o = -(weight / stretch) * change;
        }
    }
    fn discrete_gradient_rule(&self) -> Option<&dyn RelaxationStrainDiscreteGradient> {
        Some(self)
    }
}

impl RelaxationStrainDiscreteGradient for LogarithmicStrainMap {
    fn evaluate(&self, a: &[f64], b: &[f64], mid: &[f64], out: &mut [f64]) -> bool {
        out.fill(f64::NAN);
        if out.len() != self.projection.len() || self.stretch(mid).is_none() {
            return false;
        }
        let (Some(sa), Some(sb)) = (self.stretch(a), self.stretch(b)) else {
            return false;
        };
        let Some((m, mean, _)) = logarithmic_mean(sa, sb) else {
            return false;
        };
        for (o, &weight) in out.iter_mut().zip(&self.projection) {
            *o = (weight / m) * mean;
        }
        out.iter().all(|v| v.is_finite())
    }

    fn end_jacobian_action(&self, a: &[f64], b: &[f64], v: &[f64], out: &mut [f64]) -> bool {
        out.fill(f64::NAN);
        if out.len() != self.projection.len()
            || v.len() != out.len()
            || v.iter().any(|v| !v.is_finite())
        {
            return false;
        }
        let (Some(sa), Some(sb)) = (self.stretch(a), self.stretch(b)) else {
            return false;
        };
        let Some((m, _, derivative)) = logarithmic_mean(sa, sb) else {
            return false;
        };
        let mut change = Sum::default();
        for (&weight, &v) in self.projection.iter().zip(v) {
            if change.product(weight / m, v).is_none() {
                return false;
            }
        }
        let Some(change) = change.value() else {
            return false;
        };
        for (o, &weight) in out.iter_mut().zip(&self.projection) {
            *o = (weight / m) * (derivative * change);
        }
        out.iter().all(|v| v.is_finite())
    }

    fn midpoint_energy_remainder(
        &self,
        a: &[f64],
        b: &[f64],
        mid: &[f64],
        scale: f64,
        offsets: [f64; 3],
    ) -> Option<f64> {
        if !scale.is_finite() || scale <= 0.0 || offsets.iter().any(|v| !v.is_finite()) {
            return None;
        }
        // Independent original strain algebra, including represented-midpoint
        // rounding and memory work. It does not call the discrete-gradient rule.
        self.stretch(a)?;
        self.stretch(b)?;
        let sm = self.stretch(mid)?;
        let left = self.scaled_increment(mid, a, sm)?;
        let right = self.scaled_increment(mid, b, sm)?;
        let ca = log_curvature(left)?;
        let cb = log_curvature(right)?;
        let mut remainder = Sum::default();
        remainder.product(scale, cb - ca)?;
        let mut change = remainder;
        change.product(scale, right - left)?;
        change.product(-1.0, offsets[1] - offsets[0])?;
        let mut offset = Sum::default();
        offset.product(scale, left + right)?;
        offset.product(scale, ca)?;
        offset.product(scale, cb)?;
        offset.product(-1.0, (offsets[0] - offsets[2]) + (offsets[1] - offsets[2]))?;
        let mut work = Sum::default();
        work.product(0.5 * offset.value()?, change.value()?)?;
        work.product(scale * det::ln(sm) - offsets[2], remainder.value()?)?;
        work.value()
    }
}
