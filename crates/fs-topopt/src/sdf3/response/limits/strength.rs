//! Joint response inequalities and a sampled-stress cap on the SAME design.
//!
//! Let R=logsumexp(beta*c_i)/beta be the existing response upper aggregate,
//! and S=stress/cap-1. logsumexp(beta*[R,S])/beta is exactly the upper aggregate
//! of ALL original response residuals plus S in exact arithmetic. No averaging,
//! priority weights or removal of inactive limits can hide an original failure.
//! The derivative is the softmax-weighted sum of the existing exact gradients.
use super::*;
use crate::pipeline::LoadCase;
use crate::sdf3::stress::{self as stress_owner, StressError3, StressEvaluation3,
    StressMeasure3, StressOptions3};

pub mod design;

/// Explicit reference-model policy. Stress always uses SampledPeakBound,
/// including every declared experiment; a normalized average is NOT a cap on
/// each sample. The stress p-norm and outer log-sum-exp both add conservatism.
#[derive(Debug, Clone, Copy)]
pub struct ResponseStressOptions3 {
    pub response: ResponseLimitOptions3,
    pub stress: StressOptions3,
    /// Positive cap in the declared reference-material stress units.
    pub stress_limit: f64,
}
impl Default for ResponseStressOptions3 {
    fn default() -> Self {
        Self { response: ResponseLimitOptions3::default(), stress: StressOptions3::default(),
            stress_limit: 1.0 }
    }
}
#[derive(Debug, Clone, PartialEq)]
pub enum ResponseStressError3 {
    Invalid(&'static str),
    Response(ResponseError3),
    Stress(StressError3),
    Evaluation(EvaluationStop),
}
impl From<ResponseError3> for ResponseStressError3 {
    fn from(e: ResponseError3) -> Self {
        match e { ResponseError3::Evaluation(e) => Self::Evaluation(e), other => Self::Response(other) }
    }
}
impl From<StressError3> for ResponseStressError3 {
    fn from(e: StressError3) -> Self {
        match e { StressError3::Evaluation(e) => Self::Evaluation(e), other => Self::Stress(other) }
    }
}
impl From<EvaluationStop> for ResponseStressError3 {
    fn from(e: EvaluationStop) -> Self { Self::Evaluation(e) }
}
impl std::fmt::Display for ResponseStressError3 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "joint response/stress evaluation refused: {self:?}")
    }
}
impl std::error::Error for ResponseStressError3 {}

/// Complete evidence at one raw design. Both components use identical force,
/// prescribed-motion, geometry and material laws; their distinct adjoints are
/// retained, not relabeled as an adjoint of the joint aggregate.
#[derive(Debug, Clone)]
pub struct ResponseStressEvaluation3 {
    pub response: ResponseLimitEvaluation3,
    pub stress: StressEvaluation3,
    /// Joint smooth upper aggregate in dimensionless constraint coordinates.
    pub constraint: f64,
    pub constraint_gradient: Vec<f64>,
    /// [response-aggregate, stress-aggregate] chain-rule coefficients.
    pub component_weights: [f64; 2],
    pub stress_limit: f64,
    /// Stress p-norm/cap-1; distinct from an original pointwise sample residual.
    pub stress_residual: f64,
    pub sampled_relaxed_residual: f64,
    pub sampled_physical_residual: f64,
    /// Worst ORIGINAL response or retained stress-sample residual.
    pub worst_residual: f64,
    /// Cumulative work after BOTH component evaluations and final admission.
    pub work: SolveWork,
}
impl ResponseStressEvaluation3 {
    /// Requires the joint aggregate, stress aggregate, original response limits
    /// AND both separately evaluated stress maxima. This is finite-grid numerical
    /// feasibility, not a continuum stress, actuator-energy or safety certificate.
    #[must_use]
    pub fn feasible(&self, tolerance: f64) -> bool {
        self.response.feasible(tolerance)
            && [self.constraint, self.stress_residual, self.sampled_relaxed_residual,
                self.sampled_physical_residual, self.worst_residual]
                .iter().all(|v| v.is_finite() && *v <= tolerance)
    }
}

fn admit<O: AdaptiveSdf3Elasticity>(study: &CutDensityStudy3<O>, rho: &[f64],
    cases: &[ResponseLimitCase3<'_>], options: ResponseStressOptions3)
    -> Result<(), ResponseStressError3> {
    // Admit the ENTIRE policy before calling a physical law or doing a solve.
    if !options.stress_limit.is_finite() || options.stress_limit <= 0.0
        || !(1.0/options.stress_limit).is_finite() {
        return Err(ResponseStressError3::Invalid("positive invertible stress cap required"));
    }
    admit_limits(study, cases, options.response)?;
    if cases.len() > options.stress.max_cases {
        return Err(ResponseStressError3::Invalid("joint stress case allowance exhausted"));
    }
    let loads: Vec<_> = cases.iter().map(|c| LoadCase { force: c.force, weight: 1.0 }).collect();
    let motions: Vec<_> = cases.iter().map(|c| c.prescribed).collect();
    stress_owner::admit(study, rho, &loads, options.stress)?;
    stress_owner::admit_motions(study, &loads, Some(&motions))?;
    Ok(())
}
impl<O: AdaptiveSdf3Elasticity> CutDensityStudy3<O> {
    /// Compose the original reaction/response and motion-aware stress evaluators.
    /// Each independently solves the SAME experiment family, with every case
    /// present in the stress cap. Pure reference laws are mandatory. The two
    /// evaluators deliberately retain their own preparations and primal solves:
    /// this reference integration is not a fused or reduced-work implementation.
    /// All work, including a successful first phase before a failed second phase,
    /// remains charged to the caller's ONE SolveControl. Both phases restore
    /// incoming scales on every returned result; no partial joint result escapes.
    ///
    /// The exact material derivative contains reaction direct terms, stress
    /// direct terms, changes of stiffness AND prescribed-motion lifting, then
    /// the original filter/projection/SIMP pullback. No finite-difference tangent
    /// or differentiation through Krylov iterations is introduced.
    pub fn evaluate_response_stress_limits(&mut self, rho: &[f64],
        cases: &[ResponseLimitCase3<'_>], options: ResponseStressOptions3,
        control: &mut SolveControl<'_>) -> Result<ResponseStressEvaluation3, ResponseStressError3> {
        control.checkpoint("response-stress-start")?;
        admit(self, rho, cases, options)?;
        let response = self.evaluate_response_limits(rho, cases, options.response, control)?;
        control.checkpoint("response-stress-between")?;
        let loads: Vec<_> = cases.iter().map(|c| LoadCase { force: c.force, weight: 1.0 }).collect();
        let motions: Vec<_> = cases.iter().map(|c| c.prescribed).collect();
        let stress = self.evaluate_stress_with_prescribed_motion(rho, &loads, &motions,
            options.stress, StressMeasure3::SampledPeakBound, control)?;
        if response.rho != stress.rho || response.projected_rho != stress.projected_rho
            || response.scales != stress.scales || response.volume_fraction != stress.volume_fraction
            || response.volume_gradient != stress.volume_gradient {
            return Err(ResponseStressError3::Invalid("joint evaluations disagree on the material map"));
        }
        let stress_residual = finite(stress.aggregate/options.stress_limit-1.0)?;
        let sampled_relaxed_residual = finite(stress.sampled_relaxed_max/options.stress_limit-1.0)?;
        let sampled_physical_residual = finite(stress.sampled_physical_max/options.stress_limit-1.0)?;
        let (constraint, _, weights) = aggregate(&[response.constraint, stress_residual], options.response.sharpness)?;
        let mut constraint_gradient = Vec::with_capacity(self.cells());
        for (i, (&r, &s)) in response.constraint_gradient.iter().zip(&stress.gradient).enumerate() {
            if i % 192 == 0 { control.checkpoint("response-stress-combine")?; }
            constraint_gradient.push(finite(weights[0]*r+weights[1]*(s/options.stress_limit))?);
        }
        let worst_residual = response.worst_residual.max(sampled_relaxed_residual).max(sampled_physical_residual);
        control.checkpoint("response-stress-publish")?;
        Ok(ResponseStressEvaluation3 { response, stress, constraint, constraint_gradient,
            component_weights: [weights[0], weights[1]], stress_limit: options.stress_limit,
            stress_residual, sampled_relaxed_residual, sampled_physical_residual, worst_residual,
            work: control.work() })
    }
}

#[cfg(test)]
mod tests;
