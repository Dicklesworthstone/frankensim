//! Refined-grid assessment of the ERODED robust design, not nominal stiffness.
use std::ops::ControlFlow;
use super::RobustOptions3;
use super::super::{AdaptiveSdf3Elasticity, CutDensityEvaluation3, CutDensityStudy3};
use crate::pipeline::LoadCase;
use crate::{SimpParams, SolveControl};
use crate::sdf3_goal::{ComplianceRefinement3, GoalReferenceLoad3,
    GoalRefinementError3, GoalRefinementOptions3};
use fs_cutfem::elastic3::adaptive::AdaptiveElasticity3;

/// One freshly re-solved eroded design and its complete two-grid experiment.
#[derive(Debug, Clone)]
pub struct RobustComplianceRefinement3 {
    /// Nominal SIMP/projection parameters; erosion shifts eta by the offset.
    pub params: SimpParams,
    pub robustness: RobustOptions3,
    pub rho: Vec<f64>,
    /// All fields, gradients and volume here are ERODED, not nominal/dilated.
    pub eroded: CutDensityEvaluation3,
    /// Inherits eroded physical stiffness, reintegrates reference loads and
    /// retains consistency/algebraic terms. This is NOT a new fine-grid design
    /// problem, a continuum error bound, or a dilated-volume error estimate.
    pub estimate: ComplianceRefinement3,
}

impl<O: AdaptiveSdf3Elasticity> CutDensityStudy3<O> {
    /// Re-solve a raw design at an explicitly supplied robust material model,
    /// then run the EXISTING mixed-reference-load compliance DWR experiment.
    /// The coarse eroded scales are installed only during that experiment and
    /// inherited by the fine operator; fine RAW controls are never refiltered.
    /// This prevents a nominal fine solve from assessing eroded displacements.
    ///
    /// Explicit params allow assessing a continuation report while the study
    /// retains its original nominal model. Cached input fields are not trusted:
    /// the raw design is freshly solved with the exact supplied load laws.
    /// Weights remain independent and unnormalized. Reference laws must be
    /// pure, density-independent and identical on both grids. Enriched geometry,
    /// supports and base material must represent the same physical problem.
    /// Nonzero prescribed motion is not part of this homogeneous-support path.
    ///
    /// All state/filter/enriched solves share control. Success and every returned
    /// error restore the caller's original parameters and scales. No partial
    /// load family escapes. Model/input panics occur before state mutation.
    #[allow(clippy::too_many_arguments)]
    pub fn estimate_robust_reference_compliance_enrichment(
        &mut self, enriched: AdaptiveElasticity3, loads: &[GoalReferenceLoad3<'_>],
        rho: &[f64], params: SimpParams, robustness: RobustOptions3,
        options: GoalRefinementOptions3, control: &mut SolveControl<'_>,
    ) -> Result<RobustComplianceRefinement3, GoalRefinementError3> {
        if self.has_material_body_load() {
            return Err(GoalRefinementError3::Invalid("material body loads are not supported by robust goal refinement"));
        }
        robustness.assert_valid(params);
        if rho.len() != self.cells() || rho.iter().any(|r| !r.is_finite() || !(0.0..=1.0).contains(r))
            || loads.is_empty() || loads.len() > options.max_load_cases
            || !loads.iter().any(|l| l.weight > 0.0)
            || loads.iter().any(|l| !l.weight.is_finite() || l.weight < 0.0)
            || options.max_transfer_terms == 0
            || ![options.numerical.residual_tolerance, options.numerical.identity_tolerance]
                .iter().all(|v| v.is_finite() && *v > 0.0 && *v < 1.0) {
            return Err(GoalRefinementError3::Invalid("invalid robust goal design, loads or resource policy"));
        }
        control.checkpoint("sdf3-robust-goal-start")?;
        let incoming_params = self.params;
        let incoming_scales = self.operator.scales().to_vec();
        let result = (|| {
            let mut forces = Vec::with_capacity(loads.len());
            for load in loads {
                forces.push(self.operator.adaptive().reference_load(load.load, || {
                    if control.checkpoint("sdf3-robust-goal-load").is_ok() { ControlFlow::Continue(()) }
                    else { ControlFlow::Break(()) }
                })?);
            }
            let cases: Vec<_> = forces.iter().zip(loads).map(|(force, load)|
                LoadCase { force, weight: load.weight }).collect();
            self.params = params;
            self.params.eta += robustness.threshold_offset;
            // This is exactly the robust evaluator's eroded map, sharing its
            // filter/projection/SIMP implementation and physical-region rules.
            let eroded = self.evaluate(rho, &cases, control)?;
            let estimate = self.estimate_reference_compliance_enrichment(enriched,
                loads, &eroded.objective.displacements, options, control)?;
            control.checkpoint("sdf3-robust-goal-publish")?;
            Ok(RobustComplianceRefinement3 { params, robustness, rho: rho.to_vec(), eroded, estimate })
        })();
        self.params = incoming_params;
        self.operator.set_scales(&incoming_scales).expect("incoming scales remain valid");
        result
    }
}

#[cfg(test)]
#[path = "goal_tests.rs"]
mod tests;
