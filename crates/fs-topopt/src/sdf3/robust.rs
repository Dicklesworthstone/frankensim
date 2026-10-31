//! Three-field erosion/dilation robustness on real 3-D implicit-cut elasticity.
//!
//! One graph-filtered raw design produces thresholds eta+delta (eroded), eta
//! (nominal) and eta-delta (dilated). Independent external loads are solved on
//! the ERODED material; a configured material body force follows that same
//! eroded occupancy. A material budget applies to the DILATED volume. Both
//! derivatives include the exact projection and filter transpose. Prescribed
//! solid/void regions are unchanged in all realizations.
//!
//! This is a discrete manufacturing-tolerance surrogate, not a geometric
//! offset, guaranteed minimum feature size, continuum certificate, or safety
//! proof. Authored dead loads, homogeneous supports, geometry and the
//! SIMP/filter/projection model stay fixed. Neither evaluations nor numerical
//! gradient checks change the caller's nominal model or incoming scales.
//! With occupancy-dependent loading, eroded compliance need not be the largest
//! across manufacturing realizations; this owner evaluates the declared eroded
//! functional and does not assert that extra worst-case ordering.

pub mod continuation;
pub mod goal;
pub mod refit;
pub mod adaptive;
pub mod envelope;
mod optimize;
pub use optimize::{RobustDesignIteration3, RobustDesignOptions3, RobustDesignPoint3, RobustDesignStudy3};

use super::{CutDensityEvaluation3, CutDensityStudy3, Sdf3Elasticity, assert_loads, failure};
use crate::gradient_check::{GradientCheckOptions, GradientSample, MultiLoadGradientCheck, audit_design_map};
use crate::pipeline::LoadCase;
use crate::{EvaluationStop, SimpParams, SolveControl, SolveWork};

/// Symmetric uncertainty in the projection threshold, NOT a length tolerance.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RobustOptions3 {
    pub threshold_offset: f64,
}
impl Default for RobustOptions3 {
    fn default() -> Self { Self { threshold_offset: 0.1 } }
}
impl RobustOptions3 {
    fn valid(self, params: SimpParams) -> bool {
        self.threshold_offset.is_finite() && self.threshold_offset >= 0.0
            && params.eta - self.threshold_offset > 0.0
            && params.eta + self.threshold_offset < 1.0
    }
    fn assert_valid(self, params: SimpParams) {
        params.assert_valid();
        assert!(self.valid(params),
            "robust thresholds must lie strictly inside (0,1)");
    }
}

/// All three fields at ONE raw design and ONE declared uncertainty model.
#[derive(Debug, Clone)]
pub struct RobustEvaluation3 {
    pub robustness: RobustOptions3,
    /// Objective, displacements and compliance gradient belong to this field.
    /// Its volume and volume gradient are eroded, not the constraint volume.
    pub eroded: CutDensityEvaluation3,
    pub nominal_rho: Vec<f64>,
    pub nominal_volume_fraction: f64,
    pub nominal_volume_gradient: Vec<f64>,
    pub dilated_rho: Vec<f64>,
    /// Use this quantity and its gradient for the robust material constraint.
    pub dilated_volume_fraction: f64,
    pub dilated_volume_gradient: Vec<f64>,
    /// Cumulative work including every filter, independent load and pullback.
    pub work: SolveWork,
}

impl<O: Sdf3Elasticity> CutDensityStudy3<O> {
    /// Evaluate three realizations with ONE forward filter and one eroded
    /// equilibrium per independent load. Nominal/dilated fields need no extra
    /// elasticity solve. Weights are not normalized and loads are never summed.
    ///
    /// Invalid model inputs follow the existing panic-before-mutation contract.
    /// Every returned failure AND success restores incoming stiffness scales;
    /// all spent work remains charged. No partial realization family escapes.
    pub fn evaluate_robust(&mut self, rho: &[f64], loads: &[LoadCase<'_>],
        options: RobustOptions3, control: &mut SolveControl<'_>)
        -> Result<RobustEvaluation3, EvaluationStop> {
        options.assert_valid(self.params);
        assert_loads(self, loads);
        assert_eq!(rho.len(), self.cells(), "one raw density per active cut cell required");
        assert!(rho.iter().all(|r| r.is_finite() && (0.0..=1.0).contains(r)), "invalid raw density");
        control.checkpoint("sdf3-robust-evaluation")?;
        let filtered = self.filter_apply(rho, control)?;
        let nominal = self.design_from_filtered(&filtered, self.params.eta, control)?;
        let dilated = self.design_from_filtered(&filtered, self.params.eta - options.threshold_offset, control)?;
        let eroded = self.design_from_filtered(&filtered, self.params.eta + options.threshold_offset, control)?;
        // Refuse a failed ordering, rather than reorder fields and thereby
        // detach their derivatives from the physical maps that produced them.
        if eroded.projected.iter().zip(&nominal.projected).zip(&dilated.projected)
            .any(|((&e, &n), &d)| e > n + 32.0 * f64::EPSILON || n > d + 32.0 * f64::EPSILON) {
            return Err(failure("sdf3-robust-ordering"));
        }
        let nominal_local: Vec<f64> = nominal.slope.iter().zip(&self.mass).map(|(s, m)| s*m).collect();
        let dilated_local: Vec<f64> = dilated.slope.iter().zip(&self.mass).map(|(s, m)| s*m).collect();
        let nominal_volume_gradient = self.pullback(&nominal_local, control)?;
        let dilated_volume_gradient = self.pullback(&dilated_local, control)?;
        let previous = self.operator.scales().to_vec();
        let result = (|| {
            let eroded = self.evaluate_design(eroded, loads, control)?;
            control.checkpoint("sdf3-robust-publish")?;
            Ok(RobustEvaluation3 {
                robustness: options, eroded,
                nominal_rho: nominal.projected, nominal_volume_fraction: nominal.volume,
                nominal_volume_gradient,
                dilated_rho: dilated.projected, dilated_volume_fraction: dilated.volume,
                dilated_volume_gradient, work: control.work(),
            })
        })();
        self.operator.set_scales(&previous).expect("previous scales remain valid");
        result
    }
}

/// Numerical evidence for eroded compliance AND dilated material volume.
#[derive(Debug, Clone)]
pub struct RobustGradientCheck3 {
    /// Binds the extra projection-threshold input absent from SimpParams.
    pub robustness: RobustOptions3,
    pub numerical: MultiLoadGradientCheck,
}

/// Reuses the existing bound-aware heterogeneous-direction experiment, with
/// independently solved eroded load families at every finite-difference probe.
/// The returned check's `params` are the NOMINAL model; `robustness` supplies the
/// two shifted thresholds. No state is installed and no continuum claim follows.
pub fn controlled_sdf3_robust_gradient_check<O: Sdf3Elasticity>(
    study: &mut CutDensityStudy3<O>, loads: &[LoadCase<'_>], rho: &[f64],
    robustness: RobustOptions3, options: GradientCheckOptions, control: &mut SolveControl<'_>,
) -> Result<RobustGradientCheck3, EvaluationStop> {
    robustness.assert_valid(study.params);
    options.assert_valid();
    assert_loads(study, loads);
    assert_eq!(rho.len(), study.cells(), "one raw density per active cut cell required");
    assert!(!rho.is_empty() && rho.iter().all(|r| r.is_finite() && (0.0..=1.0).contains(r)),
        "gradient-check densities must be nonempty, finite and in [0,1]");
    let numerical = audit_design_map(study.params, rho, options, control, |point, _, control| {
        let evaluated = study.evaluate_robust(point, loads, robustness, control)?;
        Ok(GradientSample {
            compliance: evaluated.eroded.objective.compliance,
            compliance_gradient: evaluated.eroded.objective.gradient,
            volume: evaluated.dilated_volume_fraction,
            volume_gradient: evaluated.dilated_volume_gradient,
        })
    })?;
    Ok(RobustGradientCheck3 { robustness, numerical })
}
