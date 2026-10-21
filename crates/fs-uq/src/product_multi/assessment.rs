//! Cross-response dependence and confidence for one fixed joint event.
use core::fmt;
use fs_evidence::Color;
use crate::product_sensitivity::Sum;
use crate::{AnytimeEstimate, UqComplianceError, UqStatus};
use super::{MultiQoiExecution, QoiDefinition};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MultiQoiAssessmentError {
    RefusedExecution,
    Budget,
    Cancelled,
    NumericalRange,
    Compliance(UqComplianceError),
}
impl fmt::Display for MultiQoiAssessmentError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::RefusedExecution => f.write_str("refused paired execution cannot supply an assessment"),
            Self::Budget => f.write_str("paired covariance exceeds the declared product allowance"),
            Self::Cancelled => f.write_str("paired assessment cancelled before publication"),
            Self::NumericalRange => f.write_str("paired assessment exceeds finite numerical range"),
            Self::Compliance(error) => fmt::Display::fmt(error, f),
        }
    }
}
impl std::error::Error for MultiQoiAssessmentError {}

/// Empirical paired moments, NOT a joint probability law or PSD certificate.
#[derive(Debug, Clone, PartialEq)]
pub struct MultiQoiCovariance {
    pub samples: usize,
    /// Matrix indices and units follow the original output declaration order.
    pub quantities: Vec<QoiDefinition>,
    /// Sample covariance with divisor n-1. Entry (i,j) has units unit_i*unit_j.
    /// None means nonzero underflow/overflow, not zero covariance.
    pub covariance: Vec<Vec<Option<f64>>>,
    /// Estimated covariance of the vector sample mean, using divisor n*(n-1).
    /// Computed separately so it can survive overflow of sample covariance.
    pub mean_covariance: Vec<Vec<Option<f64>>>,
    /// Dimensionless empirical Pearson correlations. Constant columns have
    /// undefined correlations, including on their diagonal. Rounding overshoot
    /// alone is clipped to [-1,1]; no PSD repair or population claim is made.
    pub correlation: Vec<Vec<Option<f64>>>,
    pub evidence_color: Color,
}

impl MultiQoiExecution {
    /// Paired sample covariance, mean covariance and correlation without physics.
    /// At most max_products scalar products, with an additional hard cap of
    /// 67,108,864. The count n*d*(d+1)/2 is admitted BEFORE scratch allocation.
    /// Normalization uses O(n*d) scratch, already bounded by retained storage.
    /// Checks cancellation every 512 entries and between columns/pairs.
    ///
    /// Center before scaling whenever possible, preserving small changes on a
    /// large common offset. Opposite finite extremes use scaling first. Mixed
    /// physical units never enter a shared normalization. Partial-run matrices
    /// are descriptive, not optional-stopping inference. Fewer than two rows
    /// return None. Refusal/cancellation changes no retained observation.
    pub fn assess_covariance(&self, max_products: usize, mut cancelled: impl FnMut() -> bool)
        -> Result<Option<MultiQoiCovariance>, MultiQoiAssessmentError>
    {
        self.assessment_admission(&mut cancelled)?;
        let n = self.samples_accepted();
        if n < 2 { return Ok(None); }
        let d = self.quantities.len();
        let work = n.checked_mul(d).and_then(|v| v.checked_mul(d + 1)).map(|v| v / 2)
            .ok_or(MultiQoiAssessmentError::Budget)?;
        if work > max_products.min(67_108_864) { return Err(MultiQoiAssessmentError::Budget); }
        let mut columns = Vec::with_capacity(d);
        let mut scales = Vec::with_capacity(d);
        let mut variances = Vec::with_capacity(d);
        for j in 0..d {
            poll(&mut cancelled)?;
            let (column, scale) = centered(self.column(j), &mut cancelled)?;
            variances.push(dot(&column, &column, &mut cancelled)? / (n - 1) as f64);
            columns.push(column);
            scales.push(scale);
        }
        let mut covariance = vec![vec![None; d]; d];
        let mut mean_covariance = covariance.clone();
        let mut correlation = covariance.clone();
        for i in 0..d {
            for j in 0..=i {
                poll(&mut cancelled)?;
                let unit = if i == j { variances[i] }
                    else { dot(&columns[i], &columns[j], &mut cancelled)? / (n - 1) as f64 };
                let cov = rescale(unit, scales[i], scales[j]);
                let mean_unit = unit / n as f64;
                let mean_cov = if unit != 0.0 && mean_unit == 0.0 { None }
                    else { rescale(mean_unit, scales[i], scales[j]) };
                let corr = if variances[i] == 0.0 || variances[j] == 0.0 { None }
                    else if i == j { Some(1.0) }
                    else {
                        let value = unit / fs_math::det::sqrt(variances[i])
                            / fs_math::det::sqrt(variances[j]);
                        if !value.is_finite() { return Err(MultiQoiAssessmentError::NumericalRange); }
                        Some(value.clamp(-1.0, 1.0))
                    };
                covariance[i][j] = cov; covariance[j][i] = cov;
                mean_covariance[i][j] = mean_cov; mean_covariance[j][i] = mean_cov;
                correlation[i][j] = corr; correlation[j][i] = corr;
            }
        }
        poll(&mut cancelled)?;
        Ok(Some(MultiQoiCovariance {
            samples: n, quantities: self.quantities.clone(), covariance, mean_covariance, correlation,
            evidence_color: Color::Estimated {
                estimator: "empirical paired covariance; no population, PSD or model-error certificate".into(),
                dispersion: f64::INFINITY,
            },
        }))
    }

    /// Bernoulli-mixture confidence sequence for satisfying ALL declared limits.
    /// Uses the existing fs-eproc owner on each complete paired response. This is
    /// one fixed event; arbitrary output dependence is retained. No multiplication
    /// of marginal probabilities or covariance-based distribution is substituted.
    ///
    /// Fix alpha, the physical model and ALL limits before inspecting outcomes.
    /// Under the owner's fixed-conditional-probability assumptions, the interval
    /// is time-uniform for this event; selecting models/seeds/limits or several
    /// events requires separate multiplicity control. No finite-grid, rounding,
    /// solver/model-error bound or physical qualification is supplied.
    ///
    /// None means no rows. Convergence requires BOTH asymmetric distances from
    /// the empirical mean to meet target_half_width. This read-only assessment
    /// never changes Complete's meaning (the original physical budget finished).
    /// Polls each row and brackets the bounded interval inversion; O(n*d) work,
    /// O(1) extra statistical storage, no physical calls or partial publication.
    pub fn assess_joint_compliance(&self, alpha: f64, target_half_width: f64,
        mut cancelled: impl FnMut() -> bool)
        -> Result<Option<AnytimeEstimate>, MultiQoiAssessmentError>
    {
        use MultiQoiAssessmentError::Compliance;
        if !alpha.is_finite() || alpha <= 0.0 || alpha >= 1.0 || !alpha.recip().is_finite() {
            return Err(Compliance(UqComplianceError::InvalidAlpha));
        }
        if !target_half_width.is_finite() || target_half_width < 0.0 {
            return Err(Compliance(UqComplianceError::InvalidHalfWidth));
        }
        if !self.has_limits() { return Err(Compliance(UqComplianceError::MissingThreshold)); }
        self.assessment_admission(&mut cancelled)?;
        let mut cs = fs_eproc::bernoulli::BernoulliMixtureCs::new(alpha)
            .map_err(|_| MultiQoiAssessmentError::NumericalRange)?;
        for i in 0..self.samples_accepted() {
            poll(&mut cancelled)?;
            cs.observe(self.complies(i)).map_err(|_| MultiQoiAssessmentError::NumericalRange)?;
        }
        poll(&mut cancelled)?;
        let interval = cs.interval().map_err(|_| MultiQoiAssessmentError::NumericalRange)?;
        poll(&mut cancelled)?;
        Ok(interval.map(|i| AnytimeEstimate {
            mean: i.mean, lo: i.lo, hi: i.hi, n: i.n,
            converged: (i.mean - i.lo).max(i.hi - i.mean) <= target_half_width,
        }))
    }

    fn assessment_admission(&self, cancelled: &mut impl FnMut() -> bool)
        -> Result<(), MultiQoiAssessmentError>
    {
        poll(cancelled)?;
        if self.primary.status == UqStatus::Refused || self.primary.failure.is_some() {
            Err(MultiQoiAssessmentError::RefusedExecution)
        } else { Ok(()) }
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), MultiQoiAssessmentError> {
    if cancelled() { Err(MultiQoiAssessmentError::Cancelled) } else { Ok(()) }
}

// Compensated, independently normalized columns. No global output-unit scale.
fn centered(values: &[f64], cancelled: &mut impl FnMut() -> bool)
    -> Result<(Vec<f64>, f64), MultiQoiAssessmentError>
{
    let origin = values[0];
    let mut scale = 0.0_f64;
    let mut absolute = 0.0_f64;
    for (i, &x) in values.iter().enumerate() {
        if i % 512 == 0 { poll(cancelled)?; }
        scale = scale.max((x - origin).abs());
        absolute = absolute.max(x.abs());
    }
    let scale_first = !scale.is_finite();
    if scale_first { scale = absolute; }
    if scale == 0.0 { scale = 1.0; }
    let mut column = Vec::with_capacity(values.len());
    let mut total = Sum::default();
    for (i, &x) in values.iter().enumerate() {
        if i % 512 == 0 { poll(cancelled)?; }
        let u = if scale_first { x / scale - origin / scale } else { (x - origin) / scale };
        total.add(u); column.push(u);
    }
    let mean = total.value / values.len() as f64;
    for (i, value) in column.iter_mut().enumerate() {
        if i % 512 == 0 { poll(cancelled)?; }
        *value -= mean;
    }
    Ok((column, scale))
}
fn dot(a: &[f64], b: &[f64], cancelled: &mut impl FnMut() -> bool)
    -> Result<f64, MultiQoiAssessmentError>
{
    let mut sum = Sum::default();
    for (i, (&x, &y)) in a.iter().zip(b).enumerate() {
        if i % 512 == 0 { poll(cancelled)?; }
        sum.add(x * y);
    }
    if sum.value.is_finite() { Ok(sum.value) }
    else { Err(MultiQoiAssessmentError::NumericalRange) }
}
// Multiply smallest and largest magnitudes first: mixed units cannot overflow
// an intermediate product when the final nonzero result is representable.
fn rescale(value: f64, a: f64, b: f64) -> Option<f64> {
    if value == 0.0 { return Some(0.0); }
    let mut factors = [value.abs(), a, b];
    factors.sort_by(f64::total_cmp);
    let result = (factors[0] * factors[2]) * factors[1];
    (result.is_finite() && result > 0.0).then_some(result.copysign(value))
}
