//! Smooth worst-load compliance on the real three-field CutFEM design map.
//! References are performance scales, NOT probabilities. Every declared load
//! is solved independently. The existing weighted-sum APIs remain unchanged.
pub mod optimize;
pub mod limits;
pub mod load_space;
use super::{RobustEvaluation3, RobustOptions3};
use super::super::{CutDensityStudy3, Sdf3Elasticity, failure};
use crate::gradient_check::{GradientCheckOptions, GradientSample, MultiLoadGradientCheck, audit_design_map};
use crate::pipeline::LoadCase;
use crate::{EvaluationStop, SolveControl, SolveWork};

/// Immutable dimensionless maximum policy. In exact arithmetic,
/// max(C_i/ref_i) <= value <= max(C_i/ref_i) + temperature*ln(case_count).
/// This is an algebraic smoothing bound, not a floating-point enclosure or a
/// continuum certificate. Positive references have the units of compliance.
#[derive(Debug, Clone, PartialEq)]
pub struct LoadEnvelope3 {
    references: Vec<f64>,
    temperature: f64,
    bias_bound: f64,
}
impl LoadEnvelope3 {
    pub fn new(references: Vec<f64>, temperature: f64, max_load_cases: usize)
        -> Result<Self, &'static str> {
        if references.is_empty() || references.len() > max_load_cases
            || references.iter().any(|r| !r.is_finite() || *r <= 0.0)
            || !temperature.is_finite() || temperature <= 0.0 {
            return Err("invalid load-envelope references, temperature or case budget");
        }
        let bias_bound = temperature * fs_math::det::ln(references.len() as f64);
        if !bias_bound.is_finite() { return Err("unrepresentable load-envelope smoothing bound"); }
        Ok(Self { references, temperature, bias_bound })
    }
    #[must_use] pub fn references(&self) -> &[f64] { &self.references }
    #[must_use] pub fn temperature(&self) -> f64 { self.temperature }
    #[must_use] pub fn bias_bound(&self) -> f64 { self.bias_bound }

    /// Stable shifted log-sum-exp; exact ties select the first declared case.
    /// Derivatives are with respect to PHYSICAL compliances, including 1/ref_i.
    /// A failed range check never becomes a zero objective or missing scenario.
    pub fn aggregate(&self, compliances: &[f64]) -> Result<EnvelopeValue3, EvaluationStop> {
        if compliances.len() != self.references.len() {
            return Err(failure("sdf3-envelope-case-count"));
        }
        let mut ratios = Vec::with_capacity(compliances.len());
        let mut maximum = 0.0;
        let mut worst_case = 0;
        for (i, (&c, &r)) in compliances.iter().zip(&self.references).enumerate() {
            let ratio = c / r;
            if !c.is_finite() || c < 0.0 || !ratio.is_finite() {
                return Err(failure("sdf3-envelope-range"));
            }
            if ratio > maximum { maximum = ratio; worst_case = i; }
            ratios.push(ratio);
        }
        let mut weights: Vec<_> = ratios.iter()
            .map(|r| fs_math::det::exp((r - maximum) / self.temperature)).collect();
        // At least one exponent is exactly zero; no all-underflow normalization.
        let total: f64 = weights.iter().sum();
        let value = maximum + self.temperature * fs_math::det::ln(total);
        if !total.is_finite() || total < 1.0 || !value.is_finite() {
            return Err(failure("sdf3-envelope-range"));
        }
        for (w, r) in weights.iter_mut().zip(&self.references) {
            *w = (*w / total) / r;
            if !w.is_finite() { return Err(failure("sdf3-envelope-derivative")); }
        }
        Ok(EnvelopeValue3 { value, maximum, worst_case, ratios, compliance_weights: weights })
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct EnvelopeValue3 {
    /// Smooth dimensionless upper approximation to the maximum normalized case.
    pub value: f64,
    pub maximum: f64,
    /// Zero-based original case index; deterministic first-case tie breaking.
    pub worst_case: usize,
    pub ratios: Vec<f64>,
    /// d(value)/d(C_i), NOT probabilities or constant optimization weights.
    pub compliance_weights: Vec<f64>,
}
#[derive(Debug, Clone)]
pub struct EnvelopeEvaluation3 {
    pub policy: LoadEnvelope3,
    /// Actual robust fields. The embedded compliance/gradient is the UNIFORM
    /// MEAN diagnostic, NOT the envelope objective or its derivative below.
    pub physical: RobustEvaluation3,
    pub envelope: EnvelopeValue3,
    pub gradient: Vec<f64>,
    /// Cumulative work including the extra exact design-map pullback.
    pub work: SolveWork,
}

impl<O: Sdf3Elasticity> CutDensityStudy3<O> {
    /// Evaluate a normalized worst-load objective and exact discrete derivative.
    /// Uses one robust load-family solve, never finite-difference optimization
    /// or one state solve per design coordinate. Reuses eroded displacements to
    /// combine self-adjoint stiffness contractions before ONE filter pullback.
    ///
    /// Authored force vectors are density-independent, with homogeneous supports.
    /// An attached material body load uses each eroded occupancy and contributes
    /// its exact force derivative. References, case order and temperature must
    /// stay fixed throughout an optimization.
    /// Every solve shares control; incoming nominal scales remain unchanged on
    /// success and every returned stop. No partial family is published.
    pub fn evaluate_robust_envelope(&mut self, rho: &[f64], forces: &[&[f64]],
        robustness: RobustOptions3, policy: &LoadEnvelope3, control: &mut SolveControl<'_>)
        -> Result<EnvelopeEvaluation3, EvaluationStop> {
        if forces.len() != policy.references.len() {
            return Err(failure("sdf3-envelope-case-count"));
        }
        control.checkpoint("sdf3-envelope-start")?;
        let weight = 1.0 / forces.len() as f64;
        let loads: Vec<_> = forces.iter().map(|&force| LoadCase { force, weight }).collect();
        let physical = self.evaluate_robust(rho, &loads, robustness, control)?;
        let envelope = policy.aggregate(&physical.eroded.objective.case_compliances)?;
        // Reapply the SAME filter to retain the exact slope including physical
        // overrides and clamping. Do not infer it from projected densities.
        let filtered = self.filter_apply(rho, control)?;
        let design = self.design_from_filtered(&filtered, self.params.eta + robustness.threshold_offset, control)?;
        let mut local = vec![0.0; self.cells()];
        let mut occupancy = vec![0.0; self.cells()];
        for (u, &w) in physical.eroded.objective.displacements.iter().zip(&envelope.compliance_weights) {
            control.checkpoint("sdf3-envelope-adjoint")?;
            // Sealed backends implement dK/dscale, independent of the incoming
            // nominal scales, including Nitsche/ghost and hanging-node terms.
            let energies = self.operator.scale_quadratic_forms(u).map_err(|_| failure("sdf3-envelope-energy"))?;
            for (g, e) in local.iter_mut().zip(energies) { *g -= w * e; }
            // Compliance differentiates both the equilibrium load and the
            // external-work functional: d(f^T u) = 2 u^T df - u^T dK u.
            self.add_material_body_gradient(u, 2.0 * w, &mut occupancy, control)?;
        }
        let gradient = self.physical_gradient(&design, &local, &occupancy, control)?;
        control.checkpoint("sdf3-envelope-publish")?;
        Ok(EnvelopeEvaluation3 { policy: policy.clone(), physical, envelope, gradient, work: control.work() })
    }
}

/// Binds the non-nominal objective inputs to the existing numerical experiment.
#[derive(Debug, Clone)]
pub struct EnvelopeGradientCheck3 {
    pub policy: LoadEnvelope3,
    pub robustness: RobustOptions3,
    /// Its compliance quantities mean the dimensionless ENVELOPE in this check.
    pub numerical: MultiLoadGradientCheck,
}
pub fn controlled_envelope_gradient_check<O: Sdf3Elasticity>(
    study: &mut CutDensityStudy3<O>, rho: &[f64], forces: &[&[f64]],
    robustness: RobustOptions3, policy: &LoadEnvelope3, options: GradientCheckOptions,
    control: &mut SolveControl<'_>,
) -> Result<EnvelopeGradientCheck3, EvaluationStop> {
    options.assert_valid();
    robustness.assert_valid(study.params());
    let numerical = audit_design_map(study.params(), rho, options, control, |point, _, control| {
        let e = study.evaluate_robust_envelope(point, forces, robustness, policy, control)?;
        Ok(GradientSample { compliance: e.envelope.value, compliance_gradient: e.gradient,
            volume: e.physical.dilated_volume_fraction, volume_gradient: e.physical.dilated_volume_gradient })
    })?;
    Ok(EnvelopeGradientCheck3 { policy: policy.clone(), robustness, numerical })
}

#[cfg(test)]
#[path = "envelope/tests.rs"]
mod tests;
