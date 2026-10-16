//! Reuse the fixed part of a joint GP posterior across candidate probes.
//!
//! With V = L_train^-1 K_train,points, the posterior covariance is K - V^T V.
//! Factor its baseline block once; append only candidate rows of its Cholesky
//! factor on each probe. The baseline projections and factor are never rebuilt.
//! This is block Gaussian conditioning, not a sparse or independent-output
//! approximation. See BoTorch's cached-Cholesky acquisition implementation.

use super::{Gp, forward_sub};
use fs_la::factor::{Cholesky, cholesky};

/// Preparation/extension refusals; no partial factor or changed cache escapes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PreparedError {
    /// Empty/wrong-dimensional/non-finite points, invalid jitter or model shape.
    InvalidInput,
    /// A coordinate repeats the baseline or another point, including signed zero.
    DuplicatePoint,
    /// A required row-major length cannot be represented.
    ShapeOverflow,
    /// The declared regularization does not admit a positive definite covariance.
    NotPositiveDefinite,
    /// Posterior or factor arithmetic is non-finite.
    NonFinite,
    /// The caller paused preparation or extension.
    Cancelled,
}
impl std::fmt::Display for PreparedError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "prepared GP posterior: {self:?}")
    }
}
impl std::error::Error for PreparedError {}

/// Candidate means and ONLY the new rows of a full joint lower factor.
#[derive(Debug, Clone, PartialEq)]
pub struct PosteriorExtension {
    /// Latent posterior means, in candidate order.
    pub mean: Vec<f64>,
    /// Row-major `mean.len() x (baseline_len + mean.len())`. Row i is row
    /// baseline_len+i of the joint factor; its unused upper entries are zero.
    pub lower_rows: Vec<f64>,
    /// Number of baseline columns, whose factor belongs to the preparation.
    pub baseline_len: usize,
}

/// Borrow one fitted GP and cache its baseline projections and factor.
///
/// `diagonal_jitter` is an explicit, nonnegative POSTERIOR covariance variance
/// added once to every distinct baseline/candidate diagonal. It is not fitted
/// observation noise. A fixed value makes the baseline independent of probes;
/// unlike `Gp::predict_joint`, this path does not escalate jitter or discard
/// correlations on failure. It may therefore refuse where that API falls back.
/// Different regularization/arithmetic can change acquisition and design bits.
///
/// Memory is O(n*b + b*b), with n training rows and b distinct baseline points.
/// Each q-point extension needs O(n*q + b*q + q*q) additional memory, not another
/// b*b allocation. The GP cannot be mutated while this preparation is borrowed.
pub struct PreparedPosterior<'a> {
    gp: &'a Gp,
    points: Vec<Vec<f64>>,
    mean: Vec<f64>,
    projections: Vec<Vec<f64>>,
    lower: Cholesky,
    jitter: f64,
}

fn check(keep_going: &mut dyn FnMut() -> bool) -> Result<(), PreparedError> {
    if keep_going() { Ok(()) } else { Err(PreparedError::Cancelled) }
}
fn entries(rows: usize, cols: usize) -> Result<usize, PreparedError> {
    rows.checked_mul(cols).filter(|v| *v <= isize::MAX as usize / std::mem::size_of::<f64>())
        .ok_or(PreparedError::ShapeOverflow)
}
fn points_valid(points: &[Vec<f64>], dim: usize) -> Result<(), PreparedError> {
    if points.is_empty() || points.iter().any(|p|
        p.len() != dim || p.iter().any(|v| !v.is_finite())) {
        return Err(PreparedError::InvalidInput);
    }
    for (i, point) in points.iter().enumerate() {
        if points[..i].contains(point) { return Err(PreparedError::DuplicatePoint); }
    }
    Ok(())
}
fn project(gp: &Gp, point: &[f64]) -> Result<(f64, Vec<f64>), PreparedError> {
    let mut v: Vec<f64> = gp.x.iter().map(|x| gp.kernel.eval(x, point)).collect();
    let mean = v.iter().zip(&gp.alpha).map(|(a, b)| a * b).sum();
    forward_sub(&gp.chol, &mut v, gp.x.len());
    if !f64::is_finite(mean) || v.iter().any(|v| !v.is_finite()) {
        return Err(PreparedError::NonFinite);
    }
    Ok((mean, v))
}
fn covariance(gp: &Gp, a: &[f64], va: &[f64], b: &[f64], vb: &[f64],
    diagonal: bool, jitter: f64) -> Result<f64, PreparedError>
{
    let dot: f64 = va.iter().zip(vb).map(|(a, b)| a * b).sum();
    let value = gp.kernel.eval(a, b) - dot;
    if !value.is_finite() { return Err(PreparedError::NonFinite); }
    // Same marginal roundoff convention as Gp::predict_joint, but no adaptive
    // jitter. Check finiteness BEFORE max, which would otherwise hide NaN.
    let value = if diagonal { value.max(0.0) + jitter } else { value };
    if value.is_finite() { Ok(value) } else { Err(PreparedError::NonFinite) }
}

impl<'a> PreparedPosterior<'a> {
    /// Prepare distinct baseline points with the supplied covariance jitter.
    pub fn new(gp: &'a Gp, baseline: &[Vec<f64>], diagonal_jitter: f64)
        -> Result<Self, PreparedError>
    {
        Self::new_controlled(gp, baseline, diagonal_jitter, &mut || true)
    }

    /// Checks continuation before projections, covariance rows and factorization.
    /// One training triangular solve or dense Cholesky remains non-preemptible.
    pub fn new_controlled(gp: &'a Gp, baseline: &[Vec<f64>], diagonal_jitter: f64,
        keep_going: &mut dyn FnMut() -> bool) -> Result<Self, PreparedError>
    {
        let dim = gp.kernel.lengthscales.len();
        if dim == 0 || gp.x.is_empty() || gp.alpha.len() != gp.x.len()
            || gp.x.iter().any(|x| x.len() != dim || x.iter().any(|v| !v.is_finite()))
            || !gp.kernel.signal.is_finite() || gp.kernel.signal <= 0.0
            || gp.kernel.lengthscales.iter().any(|v| !v.is_finite() || *v <= 0.0)
            || !diagonal_jitter.is_finite() || diagonal_jitter < 0.0 {
            return Err(PreparedError::InvalidInput);
        }
        points_valid(baseline, dim)?;
        let b = baseline.len();
        let length = entries(b, b)?;
        entries(b, gp.x.len())?;
        let mut mean = Vec::with_capacity(b);
        let mut projections = Vec::with_capacity(b);
        for point in baseline {
            check(keep_going)?;
            let (m, v) = project(gp, point)?;
            mean.push(m);
            projections.push(v);
        }
        let mut sigma = vec![0.0; length];
        for i in 0..b {
            check(keep_going)?;
            for j in 0..=i {
                sigma[i * b + j] = covariance(gp, &baseline[i], &projections[i],
                    &baseline[j], &projections[j], i == j, diagonal_jitter)?;
            }
        }
        check(keep_going)?;
        let lower = cholesky(&sigma, b).map_err(|_| PreparedError::NotPositiveDefinite)?;
        check(keep_going)?;
        Ok(Self { gp, points: baseline.to_vec(), mean, projections, lower, jitter: diagonal_jitter })
    }

    /// Distinct baseline points in factor order.
    #[must_use]
    pub fn baseline(&self) -> &[Vec<f64>] { &self.points }
    /// Baseline posterior means.
    #[must_use]
    pub fn mean(&self) -> &[f64] { &self.mean }
    /// The retained baseline factor; strict upper queries return zero.
    #[must_use]
    pub fn factor(&self) -> &Cholesky { &self.lower }
    /// Applied posterior variance regularization, unchanged across probes.
    #[must_use]
    pub fn diagonal_jitter(&self) -> f64 { self.jitter }

    /// Compute only new joint factor rows for distinct, non-baseline candidates.
    pub fn extend(&self, candidates: &[Vec<f64>]) -> Result<PosteriorExtension, PreparedError> {
        self.extend_controlled(candidates, &mut || true)
    }

    /// Cancellable row extension. A failed/cancelled call cannot change the
    /// baseline or a future probe. No training or baseline refactorization is
    /// performed; Schur pivots must pass under the originally declared jitter.
    pub fn extend_controlled(&self, candidates: &[Vec<f64>], keep_going: &mut dyn FnMut() -> bool)
        -> Result<PosteriorExtension, PreparedError>
    {
        points_valid(candidates, self.gp.kernel.lengthscales.len())?;
        if candidates.iter().any(|p| self.points.contains(p)) {
            return Err(PreparedError::DuplicatePoint);
        }
        let b = self.points.len();
        let q = candidates.len();
        let width = b.checked_add(q).ok_or(PreparedError::ShapeOverflow)?;
        let length = entries(q, width)?;
        entries(q, self.gp.x.len())?;
        let mut mean = Vec::with_capacity(q);
        let mut projections = Vec::with_capacity(q);
        for point in candidates {
            check(keep_going)?;
            let (m, v) = project(self.gp, point)?;
            mean.push(m);
            projections.push(v);
        }
        let mut rows = vec![0.0; length];
        for i in 0..q {
            check(keep_going)?;
            for j in 0..=b + i {
                let (point, projection) = if j < b { (&self.points[j], &self.projections[j]) }
                    else { (&candidates[j - b], &projections[j - b]) };
                let mut value = covariance(self.gp, &candidates[i], &projections[i], point,
                    projection, j == b + i, self.jitter)?;
                for k in 0..j {
                    let previous = if j < b { self.lower.l(j, k) } else { rows[(j - b) * width + k] };
                    value = (-rows[i * width + k]).mul_add(previous, value);
                }
                if !value.is_finite() { return Err(PreparedError::NonFinite); }
                if j == b + i {
                    if value <= 0.0 { return Err(PreparedError::NotPositiveDefinite); }
                    rows[i * width + j] = fs_math::det::sqrt(value);
                } else {
                    let pivot = if j < b { self.lower.l(j, j) } else { rows[(j - b) * width + j] };
                    rows[i * width + j] = value / pivot;
                }
                if !rows[i * width + j].is_finite() { return Err(PreparedError::NonFinite); }
            }
        }
        check(keep_going)?;
        Ok(PosteriorExtension { mean, lower_rows: rows, baseline_len: b })
    }
}

#[cfg(test)]
mod tests;
