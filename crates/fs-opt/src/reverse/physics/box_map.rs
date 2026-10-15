//! Smooth finite-interval coordinates for any explicitly bound physics model.
//!
//! `y_i = lower_i + (upper_i-lower_i) sigmoid(x_i)`. Values and total
//! derivatives still come from the original provider. The transformed Hessian
//! is `a * H_y(a*v) + b * g_y * v`, with `a=dy/dx` and `b=d2y/dx2`.
//! The diagonal coordinate-curvature term is not optional.

use super::{PhysicsError, PhysicsHessian, PhysicsModel, PhysicsSample, PhysicsSignature};
use fs_exec::Cx;
use fs_qty::Dims;

fn poll(cx: Option<&Cx<'_>>) -> Result<(), PhysicsError> {
    if let Some(cx) = cx { cx.checkpoint().map_err(|_| PhysicsError::Cancelled)?; }
    Ok(())
}
fn finite(value: f64, quantity: &'static str) -> Result<f64, PhysicsError> {
    if value.is_finite() { Ok(value) }
    else { Err(PhysicsError::NonFinite { quantity, bits: value.to_bits() }) }
}

/// An immutable open box, with limits expressed in `physical_dims`.
/// New optimization coordinates are dimensionless. Finite inputs that round
/// onto an endpoint or lose a representable positive slope are refused; this
/// avoids reporting a spurious stationary point from numerical saturation.
///
/// This is a smooth coordinate chart, NOT a closed-box active-set algorithm.
/// Its numerical KKT tests are in chart coordinates. Boundary stationarity,
/// global optimality, and the underlying physics domain are separate obligations.
#[derive(Debug, Clone)]
pub struct LogisticBox {
    lower: Vec<f64>,
    upper: Vec<f64>,
    widths: Vec<f64>,
    physical_dims: Dims,
}
struct Linearization {
    physical: Vec<f64>,
    first: Vec<f64>,
    second: Vec<f64>,
}
impl LogisticBox {
    /// Validate every finite, nonempty interval, including a representable
    /// interior midpoint and width. No provider call or solve is performed.
    pub fn new(lower: &[f64], upper: &[f64], physical_dims: Dims) -> Result<Self, PhysicsError> {
        if lower.is_empty() || lower.len() != upper.len() {
            return Err(PhysicsError::InvalidInput("nonempty matching box limits required"));
        }
        let mut widths = Vec::with_capacity(lower.len());
        for (&lo, &hi) in lower.iter().zip(upper) {
            if !(lo.is_finite() && hi.is_finite() && lo < hi) {
                return Err(PhysicsError::InvalidInput("box intervals must be finite and strictly ordered"));
            }
            let width = finite(hi-lo, "box width")?;
            let mid = lo + width*0.5;
            if !(mid > lo && mid < hi && width*0.25 > 0.0) {
                return Err(PhysicsError::InvalidInput("box has no resolved interior and positive slope"));
            }
            widths.push(width);
        }
        Ok(Self { lower: lower.to_vec(), upper: upper.to_vec(), widths, physical_dims })
    }
    pub fn dimensions(&self) -> usize { self.lower.len() }
    pub fn physical_dims(&self) -> Dims { self.physical_dims }
    pub fn lower(&self) -> &[f64] { &self.lower }
    pub fn upper(&self) -> &[f64] { &self.upper }

    /// Encode a strictly interior physical point without an overflowing ratio.
    /// The returned point is checked against the same representability gate as
    /// decoding. Rounding means this is a numerical, not bitwise, inverse.
    pub fn encode(&self, physical: &[f64], cx: Option<&Cx<'_>>) -> Result<Vec<f64>, PhysicsError> {
        poll(cx)?;
        if physical.len() != self.dimensions() {
            return Err(PhysicsError::InvalidInput("box physical point shape mismatch"));
        }
        let mut point = Vec::with_capacity(physical.len());
        for (i, &y) in physical.iter().enumerate() {
            if i % 256 == 0 { poll(cx)?; }
            finite(y, "box physical coordinate")?;
            if !(y > self.lower[i] && y < self.upper[i]) {
                return Err(PhysicsError::Domain("physical point must be strictly inside the box"));
            }
            point.push(finite(fs_math::det::ln(y-self.lower[i])
                - fs_math::det::ln(self.upper[i]-y), "box logit")?);
        }
        self.linearize(&point, cx)?;
        Ok(point)
    }
    pub fn decode(&self, point: &[f64], cx: Option<&Cx<'_>>) -> Result<Vec<f64>, PhysicsError> {
        Ok(self.linearize(point, cx)?.physical)
    }
    /// Borrow the original provider; cloning this wrapper cannot duplicate or
    /// reset its physics work counters. Signature dimensions are checked now.
    pub fn bind<'a>(&'a self, model: &'a dyn PhysicsModel) -> Result<LogisticPhysicsModel<'a>, PhysicsError> {
        let signature = model.signature();
        if signature.parameter_count != self.dimensions() || signature.parameter_dims != self.physical_dims {
            return Err(PhysicsError::InvalidInput("box shape or physical dimensions do not match provider"));
        }
        Ok(LogisticPhysicsModel { coordinates: self, model })
    }
    fn linearize(&self, point: &[f64], cx: Option<&Cx<'_>>) -> Result<Linearization, PhysicsError> {
        poll(cx)?;
        if point.len() != self.dimensions() {
            return Err(PhysicsError::InvalidInput("box logit point shape mismatch"));
        }
        let mut physical = Vec::with_capacity(point.len());
        let mut first = Vec::with_capacity(point.len());
        let mut second = Vec::with_capacity(point.len());
        for (i, &x) in point.iter().enumerate() {
            if i % 256 == 0 { poll(cx)?; }
            finite(x, "box logit coordinate")?;
            // exp is evaluated only at a nonpositive argument. Compute small
            // tail probabilities directly rather than subtracting from one.
            let tail = fs_math::det::exp(-x.abs());
            let denominator = 1.0+tail;
            let distance = self.widths[i]*(tail/denominator);
            let y = if x >= 0.0 { self.upper[i]-distance } else { self.lower[i]+distance };
            let a = distance/denominator;
            let b = a * if x >= 0.0 { (tail-1.0)/denominator } else { (1.0-tail)/denominator };
            if !(y > self.lower[i] && y < self.upper[i] && a > 0.0) {
                return Err(PhysicsError::Domain("box coordinate numerically saturated; move toward the interior"));
            }
            physical.push(finite(y, "decoded box coordinate")?);
            first.push(finite(a, "box slope")?);
            second.push(finite(b, "box curvature")?);
        }
        poll(cx)?;
        Ok(Linearization { physical, first, second })
    }
}

/// A coordinate adapter, not a new physical solver or a second cache. Scientific
/// output and study name are unchanged; input coordinates become dimensionless.
/// Second-order capability is exposed only when the underlying provider has it.
/// Temporary storage is O(n); there is no dense Hessian or finite-difference path.
#[derive(Debug, Clone, Copy)]
pub struct LogisticPhysicsModel<'a> {
    coordinates: &'a LogisticBox,
    model: &'a dyn PhysicsModel,
}
impl LogisticPhysicsModel<'_> {
    fn sample(&self, physical: &[f64], cx: Option<&Cx<'_>>) -> Result<PhysicsSample, PhysicsError> {
        let sample = self.model.value_gradient(physical, cx)?;
        poll(cx)?;
        if sample.gradient.len() != self.coordinates.dimensions() {
            return Err(PhysicsError::InvalidInput("physical gradient shape mismatch"));
        }
        finite(sample.value, "physical value")?;
        for (i, &g) in sample.gradient.iter().enumerate() {
            if i % 256 == 0 { poll(cx)?; } finite(g, "physical gradient")?;
        }
        Ok(sample)
    }
}
impl PhysicsModel for LogisticPhysicsModel<'_> {
    fn signature(&self) -> PhysicsSignature<'_> {
        PhysicsSignature { parameter_dims: Dims::NONE, ..self.model.signature() }
    }
    fn value_gradient(&self, point: &[f64], cx: Option<&Cx<'_>>) -> Result<PhysicsSample, PhysicsError> {
        let map = self.coordinates.linearize(point, cx)?;
        let sample = self.sample(&map.physical, cx)?;
        let mut gradient = Vec::with_capacity(point.len());
        for (i, (&a, &g)) in map.first.iter().zip(&sample.gradient).enumerate() {
            if i % 256 == 0 { poll(cx)?; }
            gradient.push(finite(a*g, "box gradient")?);
        }
        poll(cx)?;
        Ok(PhysicsSample { value: sample.value, gradient })
    }
    fn second_order(&self) -> Option<&dyn PhysicsHessian> {
        self.model.second_order().map(|_| self as &dyn PhysicsHessian)
    }
}
impl PhysicsHessian for LogisticPhysicsModel<'_> {
    fn hessian_vector_product(&self, point: &[f64], direction: &[f64], cx: Option<&Cx<'_>>)
        -> Result<Vec<f64>, PhysicsError>
    {
        poll(cx)?;
        if direction.len() != self.coordinates.dimensions() {
            return Err(PhysicsError::InvalidInput("box Hessian direction shape mismatch"));
        }
        for (i, &v) in direction.iter().enumerate() {
            if i % 256 == 0 { poll(cx)?; } finite(v, "box Hessian direction")?;
        }
        let hessian = self.model.second_order().ok_or(PhysicsError::InvalidInput("provider has no second-order capability"))?;
        let map = self.coordinates.linearize(point, cx)?;
        // Even a zero direction must validate the underlying physical domain.
        let sample = self.sample(&map.physical, cx)?;
        if direction.iter().all(|&v| v == 0.0) { poll(cx)?; return Ok(vec![0.0; direction.len()]); }
        let mut tangent = Vec::with_capacity(direction.len());
        for (i, (&a, &v)) in map.first.iter().zip(direction).enumerate() {
            if i % 256 == 0 { poll(cx)?; }
            tangent.push(finite(a*v, "physical tangent")?);
        }
        let physical = hessian.hessian_vector_product(&map.physical, &tangent, cx)?;
        poll(cx)?;
        if physical.len() != direction.len() {
            return Err(PhysicsError::InvalidInput("physical Hessian product shape mismatch"));
        }
        let mut product = Vec::with_capacity(direction.len());
        for i in 0..direction.len() {
            if i % 256 == 0 { poll(cx)?; }
            finite(physical[i], "physical Hessian product")?;
            product.push(finite(map.first[i]*physical[i]+map.second[i]*sample.gradient[i]*direction[i],
                "box Hessian product")?);
        }
        poll(cx)?;
        Ok(product)
    }
}

#[cfg(test)]
mod tests;
