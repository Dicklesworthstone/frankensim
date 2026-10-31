//! Physical-stress acceptance of freshly solved two-phase response endpoints.
//! No gray qp-relaxation, threshold gradient, or cached-equilibrium admission.
use super::*;
use crate::sdf3::response::limits::strength::ResponseStressError3;
use crate::sdf3::stress::{von_mises, StressError3};
use std::ops::ControlFlow;

mod scan;
pub use scan::{BinaryStressCandidate3, BinaryStressScan3, BinaryStressScanOptions3};

/// Fixed-grid endpoint requirements, independent of the gray design's stress
/// relaxation. Physical stress is tested in BOTH solid and positive-ersatz cells.
#[derive(Debug, Clone, Copy)]
pub struct BinaryResponseStressOptions3 {
    pub endpoint: BinaryEndpointOptions3,
    /// Positive physical stress cap, in the reference material's declared units.
    pub stress_limit: f64,
    /// Unweighted p-norm over all retained bulk samples and experiments.
    /// No volume/case averaging can hide a failed sample. Sampling changes this
    /// conservative finite-set aggregate; it is NOT a continuum bound.
    pub aggregation_power: f64,
    /// Whole-family sampling allowance; each case's point vector is released
    /// before sampling the next. Existing cumulative solve budgets also apply.
    pub max_points: usize,
}
impl Default for BinaryResponseStressOptions3 {
    fn default() -> Self {
        Self { endpoint: BinaryEndpointOptions3::default(), stress_limit: 1.0,
            aggregation_power: 8.0, max_points: 1_000_000 }
    }
}
impl BinaryResponseStressOptions3 {
    fn validate(self) -> Result<(), ResponseStressError3> {
        self.endpoint.validate()?;
        if !self.stress_limit.is_finite() || self.stress_limit <= 0.0
            || !(1.0/self.stress_limit).is_finite() || !self.aggregation_power.is_finite()
            || !(2.0..=64.0).contains(&self.aggregation_power) || self.max_points == 0 {
            return Err(ResponseStressError3::Invalid("invalid binary physical-stress policy"));
        }
        Ok(())
    }
}

/// Original bulk quadrature, including weak cells. Zero means no nonzero stress
/// in that phase; it does not establish the absence of a geometric phase.
#[derive(Debug, Clone, PartialEq)]
pub struct BinaryStressCase3 {
    pub points: usize,
    pub sampled_physical_max: f64,
    pub solid_max: f64,
    pub void_max: f64,
}
/// Re-solved response evidence plus actual physical stress at exactly its scales.
/// No differentiable density or stress sensitivity is attached to a threshold.
#[derive(Debug, Clone)]
pub struct BinaryResponseStressEvaluation3 {
    pub response: BinaryLimitEvaluation3,
    pub cases: Vec<BinaryStressCase3>,
    pub point_count: usize,
    pub aggregation_power: f64,
    pub stress_limit: f64,
    pub stress_aggregate: f64,
    pub sampled_physical_max: f64,
    pub stress_residual: f64,
    pub sampled_physical_residual: f64,
    /// Same upper log-sum-exp of response residuals plus the stress residual.
    pub constraint: f64,
    /// Worst original response or actual physical-stress sample residual.
    pub worst_residual: f64,
    pub work: SolveWork,
}
impl BinaryResponseStressEvaluation3 {
    /// Joint aggregate AND every original response AND physical sample must pass.
    /// This does not certify continuum stress, zero-void strength or manufacturability.
    #[must_use]
    pub fn feasible(&self, tolerance: f64) -> bool {
        self.response.feasible(tolerance)
            && [self.constraint, self.worst_residual, self.stress_residual, self.sampled_physical_residual]
                .iter().all(|v| v.is_finite() && *v <= tolerance)
    }
}

// Rescaled streaming p-norm. Retain just the maximum and normalized power sum,
// not every point from every case. At least one summand equals one when m>0.
// Changing m is an algebraic rescaling, not clipping or a changed functional.
#[derive(Default)]
struct StressNorm { maximum: f64, sum: f64 }
impl StressNorm {
    fn add(&mut self, value: f64, p: f64) -> Result<(), ResponseStressError3> {
        if !value.is_finite() || value < 0.0 {
            return Err(ResponseStressError3::Invalid("invalid physical von Mises stress"));
        }
        if value > self.maximum {
            self.sum = if self.maximum == 0.0 { 1.0 }
                else { self.sum*fs_math::det::pow(self.maximum/value, p)+1.0 };
            self.maximum = value;
        } else if value > 0.0 {
            self.sum += fs_math::det::pow(value/self.maximum, p);
        }
        if !self.sum.is_finite() { return Err(ResponseStressError3::Invalid("binary stress sum overflow")); }
        Ok(())
    }
    fn value(&self, p: f64) -> Result<f64, ResponseStressError3> {
        let value = if self.maximum == 0.0 { 0.0 }
            else { self.maximum*fs_math::det::pow(self.sum, 1.0/p) };
        if !value.is_finite() || value < self.maximum {
            return Err(ResponseStressError3::Invalid("binary stress aggregate failed its finite-sample bound"));
        }
        Ok(value)
    }
}
fn stress_poll(control: &mut SolveControl<'_>) -> ControlFlow<()> {
    if control.checkpoint("response-binary-stress-physics").is_ok() { ControlFlow::Continue(()) }
    else { ControlFlow::Break(()) }
}

impl<O: AdaptiveSdf3Elasticity> CutDensityStudy3<O> {
    /// Threshold and solve the original response experiment family, then sample
    /// its actual physical stress. Every field is fresh for this exact binary
    /// material; no gray field, relaxed stress, adjoint or optimizer cache enters.
    /// The original endpoint owner handles thresholds, ties, protected regions,
    /// current-material motion lifting and true equilibrium-residual admission.
    ///
    /// Stress sampling reuses those fields with their exact stiffness installed;
    /// it adds NO preparation, primal solve or geometry integration. Both phases
    /// restore incoming scales on every returned result. Sampling failures retain
    /// consumed solve work but do not publish a partial joint endpoint. All original
    /// dead loads and support patches remain, including on weak cells; material
    /// body force follows the binary physical mask. Pure callback
    /// panics are not caught. The cap covers retained samples, not the continuum.
    pub fn evaluate_binary_response_stress_limits(&mut self, rho: &[f64],
        cases: &[ResponseLimitCase3<'_>], options: BinaryResponseStressOptions3,
        control: &mut SolveControl<'_>) -> Result<BinaryResponseStressEvaluation3, ResponseStressError3> {
        control.checkpoint("response-binary-stress-start")?;
        options.validate()?;
        let response = self.evaluate_binary_response_limits(rho, cases, options.endpoint, control)?;
        self.binary_physical_stress(response, options, control)
    }

    // Private: the only callers just obtained this endpoint from the original
    // solver in the same synchronous evaluation. Not an API for cached fields.
    fn binary_physical_stress(&mut self, response: BinaryLimitEvaluation3,
        options: BinaryResponseStressOptions3, control: &mut SolveControl<'_>)
        -> Result<BinaryResponseStressEvaluation3, ResponseStressError3> {
        control.checkpoint("response-binary-stress-sample")?;
        let previous = self.operator.scales().to_vec();
        self.operator.set_scales(&response.stiffness_scales()).map_err(ResponseError3::from)?;
        let result = (|| {
            let mut norm = StressNorm::default();
            let mut point_count = 0usize;
            let mut cases = Vec::with_capacity(response.displacements.len());
            for u in &response.displacements {
                control.checkpoint("response-binary-stress-case")?;
                let points = self.operator.bulk_stress(u, options.max_points-point_count,
                    || stress_poll(control)).map_err(StressError3::from)?;
                if points.is_empty() { return Err(ResponseStressError3::Invalid("empty binary stress quadrature")); }
                point_count = point_count.checked_add(points.len())
                    .filter(|n| *n <= options.max_points).ok_or(StressError3::PointBudget)?;
                let mut case = BinaryStressCase3 { points: points.len(), sampled_physical_max: 0.0,
                    solid_max: 0.0, void_max: 0.0 };
                for (i, point) in points.iter().enumerate() {
                    if i % 256 == 0 { control.checkpoint("response-binary-stress-measure")?; }
                    let value = von_mises(&point.stress).0;
                    norm.add(value, options.aggregation_power)?;
                    case.sampled_physical_max = case.sampled_physical_max.max(value);
                    if response.solid[point.cell] { case.solid_max = case.solid_max.max(value); }
                    else { case.void_max = case.void_max.max(value); }
                }
                cases.push(case);
            }
            let stress_aggregate = norm.value(options.aggregation_power)?;
            let stress_residual = finite(stress_aggregate/options.stress_limit-1.0)?;
            let sampled_physical_residual = finite(norm.maximum/options.stress_limit-1.0)?;
            let constraint = aggregate(&[response.constraint, stress_residual], options.endpoint.limits.sharpness)?.0;
            let worst_residual = response.worst_residual.max(sampled_physical_residual);
            control.checkpoint("response-binary-stress-publish")?;
            Ok(BinaryResponseStressEvaluation3 { response, cases, point_count,
                aggregation_power: options.aggregation_power, stress_limit: options.stress_limit,
                stress_aggregate, sampled_physical_max: norm.maximum, stress_residual,
                sampled_physical_residual, constraint, worst_residual, work: control.work() })
        })();
        self.operator.set_scales(&previous).expect("incoming binary stress material remains valid");
        result
    }
}

#[cfg(test)]
mod tests;
