//! Same-material compliance descent with a sampled-stress constraint at EVERY
//! accepted update, not only at the end of a complete optimization trajectory.
//!
//! The feasible baseline and each improving projected candidate are independently
//! solved by the existing stress evaluator. A stress refusal contracts the next
//! proposal from the last accepted state; it never publishes the refused field.
//! This is a bounded feasible-design search, not a stress-adjoint/KKT method.
//! The stress bound applies only to the evaluator's deterministic sample set.

use std::convert::Infallible;
use std::ops::ControlFlow;

use fs_cutfem::CutFemError;

use crate::checkpoint::OptimizeCheckpoint;
use crate::stress::evaluate_sampled_stress_controlled;
use crate::projected::{ProjectedOptimizer, ProjectedProgress, ProjectedSettings, ProjectedStage};
use crate::volume::VolumeProjectionSettings;
use crate::{EvaluatedFinalState, SampledStressEvaluation, SampledStressLimit, evaluate_sampled_stress};

mod lifecycle;
pub use lifecycle::ProjectedStressSetupStage;

/// Stress evaluation for one otherwise improving, same-area candidate.
#[derive(Debug, Clone)]
pub struct ProjectedStressCheck {
    /// Candidate index in the shared bounded proposal search.
    pub index: usize,
    /// Independent mechanics and sample-scoped stress, when the solve succeeded.
    pub evaluation: Option<SampledStressEvaluation>,
    /// Exact failed constraint or numerical refusal; `None` means admitted.
    pub refusal: Option<String>,
}

/// The shared descent outcome plus the stress checks performed for this update.
#[derive(Debug, Clone)]
pub struct ProjectedStressUpdate {
    /// Accepted, iteration-limited, or no improving feasible candidate.
    pub progress: ProjectedProgress,
    /// Only candidates that passed area/compliance gates require a stress solve.
    pub stress_checks: Vec<ProjectedStressCheck>,
}

/// An optimizer whose retained geometry satisfies both declared constraints.
///
/// The baseline is the feasible design at construction/resume, not an earlier
/// overfilled design. Stress and geometry cannot be mutated through this API.
#[derive(Debug, Clone)]
pub struct ProjectedStressOptimizer {
    optimizer: ProjectedOptimizer,
    limit: SampledStressLimit,
    baseline: SampledStressEvaluation,
    current: SampledStressEvaluation,
}

fn refused(what: impl Into<String>) -> CutFemError {
    CutFemError::InvalidElasticityInput { what: what.into() }
}

fn admission_reason(
    evaluation: &SampledStressEvaluation,
    state: EvaluatedFinalState,
    limit: SampledStressLimit,
) -> Option<String> {
    if evaluation.snapshot != state.snapshot {
        return Some("stress evaluation does not describe the projected geometry".to_string());
    }
    if evaluation.sample_count == 0
        || !evaluation.sampled_max_von_mises.is_finite()
        || evaluation.sampled_max_von_mises < 0.0
    {
        return Some("stress evaluation has no finite nonnegative sampled maximum".to_string());
    }
    if evaluation.sampled_max_von_mises > limit.admitted_max() {
        return Some(format!(
            "sampled von Mises stress {:.17e} exceeds admitted limit {:.17e}",
            evaluation.sampled_max_von_mises, limit.admitted_max(),
        ));
    }
    None
}

impl ProjectedStressOptimizer {
    /// Add a stress constraint to an already area-feasible optimizer.
    ///
    /// Independently evaluates and admits its CURRENT geometry before any
    /// update. An overstressed baseline is refused, not silently called feasible.
    ///
    /// # Errors
    /// Refuses invalid public limit fields, a failed solve, or baseline stress.
    pub fn new(optimizer: ProjectedOptimizer, limit: SampledStressLimit) -> Result<Self, CutFemError> {
        let limit = SampledStressLimit::new(limit.max_von_mises, limit.absolute_tolerance)?;
        let checkpoint = optimizer.checkpoint();
        let current = evaluate_sampled_stress(
            checkpoint.geometry(), checkpoint.fixture(), checkpoint.settings(),
        )?;
        if let Some(reason) = admission_reason(&current, optimizer.current(), limit) {
            return Err(refused(format!("stress-constrained baseline refused: {reason}")));
        }
        Ok(Self { optimizer, limit, baseline: current.clone(), current })
    }

    /// Resume an exact feasible checkpoint under explicitly supplied constraints.
    ///
    /// The geometry, ordinal and AL multiplier are preserved. Constraint/search
    /// settings must be retained by the caller alongside the checkpoint. The
    /// resumed baseline denotes this segment; historical evidence is not minted.
    ///
    /// # Errors
    /// Refuses invalid constraints, changed fixed nodes, area, or sampled stress.
    pub fn from_checkpoint(
        checkpoint: OptimizeCheckpoint,
        fixed: Vec<(usize, f64)>,
        projection: VolumeProjectionSettings,
        controls: ProjectedSettings,
        limit: SampledStressLimit,
    ) -> Result<Self, CutFemError> {
        // Refuse a malformed limit before running projection or any PDE solve.
        let limit = SampledStressLimit::new(limit.max_von_mises, limit.absolute_tolerance)?;
        let optimizer = ProjectedOptimizer::from_checkpoint(checkpoint, fixed, projection, controls)?;
        Self::new(optimizer, limit)
    }

    /// Exact accepted geometry, global iteration and proposal multiplier.
    #[must_use]
    pub fn checkpoint(&self) -> &OptimizeCheckpoint { self.optimizer.checkpoint() }

    /// Independently evaluated same-area baseline for this segment.
    #[must_use]
    pub fn baseline(&self) -> &SampledStressEvaluation { &self.baseline }

    /// Independently evaluated mechanics/stress of the retained geometry.
    #[must_use]
    pub fn current(&self) -> &SampledStressEvaluation { &self.current }

    /// Admitted fixed stress policy, including its explicit numerical allowance.
    #[must_use]
    pub const fn limit(&self) -> SampledStressLimit { self.limit }

    /// Attempt one stress-limited, same-area accepted update.
    ///
    /// # Errors
    /// Propagates checkpoint admission errors. Candidate solve/constraint failures
    /// remain visible in the bounded attempt history without changing state.
    pub fn advance_one(&mut self) -> Result<ProjectedStressUpdate, CutFemError> {
        match self.advance_one_controlled(|_| ControlFlow::<Infallible>::Continue(()))? {
            ControlFlow::Continue(update) => Ok(update),
            ControlFlow::Break(never) => match never {},
        }
    }

    /// Transactional controlled update using the shared descent checkpoints.
    ///
    /// Cancellation reaches both final CG solves and each cell's stress probes.
    /// No accepted geometry, multiplier, ordinal or stress changes on a stop,
    /// even after the full candidate stress evaluation. Assembly, area quadrature
    /// and one cell's sampling remain indivisible; this is not a hard deadline.
    ///
    /// # Errors
    /// Propagates checkpoint admission errors without publishing a partial state.
    pub fn advance_one_controlled<B>(
        &mut self,
        control: impl FnMut(ProjectedStage) -> ControlFlow<B>,
    ) -> Result<ControlFlow<B, ProjectedStressUpdate>, CutFemError> {
        let fixture = self.checkpoint().fixture();
        let settings = self.checkpoint().settings();
        let limit = self.limit;
        let poll_iters = self.optimizer.controls().poll_iters;
        let mut next = self.current.clone();
        let mut stress_checks = Vec::new();
        let progress = self.optimizer.advance_one_admitted_controlled(|index, geometry, state, control| {
            let (evaluation, refusal) = match evaluate_sampled_stress_controlled(
                geometry, fixture, settings, poll_iters,
                |stage| control(ProjectedStage::Stress(index, stage)),
            ) {
                Ok(ControlFlow::Break(reason)) => return Ok(ControlFlow::Break(reason)),
                Ok(ControlFlow::Continue(evaluation)) => {
                    let refusal = admission_reason(&evaluation, state, limit);
                    if refusal.is_none() { next = evaluation.clone(); }
                    (Some(evaluation), refusal)
                }
                Err(error) => (None, Some(format!("stress solve refused: {error}"))),
            };
            stress_checks.push(ProjectedStressCheck { index, evaluation, refusal: refusal.clone() });
            Ok(ControlFlow::Continue(refusal))
        }, control)?;
        match progress {
            ControlFlow::Break(reason) => Ok(ControlFlow::Break(reason)),
            ControlFlow::Continue(progress) => {
                if matches!(&progress, ProjectedProgress::Accepted(_)) {
                    self.current = next;
                }
                Ok(ControlFlow::Continue(ProjectedStressUpdate { progress, stress_checks }))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Cantilever, GridSdf, OptimizeSettings};

    fn projection() -> VolumeProjectionSettings {
        VolumeProjectionSettings { target: 0.6, tolerance: 1e-4, max_shift: 2.0, max_evaluations: 64 }
    }

    fn controls() -> ProjectedSettings {
        ProjectedSettings { max_candidates: 8, ..ProjectedSettings::default() }
    }

    fn fixed(phi: &GridSdf) -> Vec<(usize, f64)> {
        phi.nodes().iter().copied().enumerate().filter(|(i, _)| i % 9 == 0 || i % 9 == 8).collect()
    }

    fn optimizer() -> ProjectedOptimizer {
        let phi = GridSdf::from_fn(8, &|_, y| (y - 0.5).abs() - 0.35);
        let fixed = fixed(&phi);
        ProjectedOptimizer::new(phi, Cantilever { load: 1.0, band: 0.125 },
            OptimizeSettings { level: 3, iterations: 2, volfrac: 0.6, move_cells: 0.1,
                nucleation_period: 0, ..OptimizeSettings::default() },
            fixed, projection(), controls()).unwrap()
    }

    fn loose_limit() -> SampledStressLimit { SampledStressLimit::new(1e12, 0.0).unwrap() }

    #[test]
    fn overstressed_baseline_and_forged_limit_refuse() {
        let source = optimizer();
        let stress = evaluate_sampled_stress(source.checkpoint().geometry(),
            source.checkpoint().fixture(), source.checkpoint().settings()).unwrap();
        assert!(stress.sampled_max_von_mises > 0.0);
        let limit = SampledStressLimit::new(0.5 * stress.sampled_max_von_mises, 0.0).unwrap();
        assert!(ProjectedStressOptimizer::new(source.clone(), limit).is_err());
        assert!(ProjectedStressOptimizer::new(source,
            SampledStressLimit { max_von_mises: f64::NAN, absolute_tolerance: 0.0 }).is_err());
    }

    #[test]
    fn accepted_update_has_same_area_lower_compliance_and_replayed_stress() {
        let mut optimizer = ProjectedStressOptimizer::new(optimizer(), loose_limit()).unwrap();
        let before = optimizer.current().clone();
        let update = optimizer.advance_one().unwrap();
        let ProjectedProgress::Accepted(step) = update.progress else { panic!("informative beam step") };
        assert!(optimizer.current().compliance < before.compliance);
        assert!((optimizer.current().volume - 0.6).abs() <= 1e-4);
        assert_eq!(step.state.snapshot, optimizer.current().snapshot);
        assert!(update.stress_checks.iter().any(|check| check.refusal.is_none()));
        let replay = evaluate_sampled_stress(optimizer.checkpoint().geometry(),
            optimizer.checkpoint().fixture(), optimizer.checkpoint().settings()).unwrap();
        assert_eq!(&replay, optimizer.current());
        assert!(replay.sampled_max_von_mises <= optimizer.limit().admitted_max());
    }

    #[test]
    fn cancel_after_stress_solve_preserves_geometry_multiplier_and_stress() {
        let mut optimizer = ProjectedStressOptimizer::new(optimizer(), loose_limit()).unwrap();
        let before = optimizer.clone();
        let result = optimizer.advance_one_controlled(|stage| {
            if matches!(stage, ProjectedStage::Publish(_)) { ControlFlow::Break("cancel") }
            else { ControlFlow::Continue(()) }
        }).unwrap();
        assert!(matches!(result, ControlFlow::Break("cancel")));
        assert_eq!(optimizer.current(), before.current());
        assert_eq!(optimizer.checkpoint().geometry().nodes(), before.checkpoint().geometry().nodes());
        assert_eq!(optimizer.checkpoint().next_iteration(), 0);
        assert_eq!(optimizer.checkpoint().ell().to_bits(), before.checkpoint().ell().to_bits());
    }

    #[test]
    fn resumed_segment_reproduces_the_uninterrupted_accepted_endpoint() {
        let mut uninterrupted = ProjectedStressOptimizer::new(optimizer(), loose_limit()).unwrap();
        assert!(matches!(uninterrupted.advance_one().unwrap().progress, ProjectedProgress::Accepted(_)));
        let checkpoint = uninterrupted.checkpoint().clone();
        let mut resumed = ProjectedStressOptimizer::from_checkpoint(checkpoint.clone(),
            fixed(checkpoint.geometry()), projection(), controls(), loose_limit()).unwrap();
        let a = uninterrupted.advance_one().unwrap();
        let b = resumed.advance_one().unwrap();
        assert_eq!(std::mem::discriminant(&a.progress), std::mem::discriminant(&b.progress));
        assert_eq!(uninterrupted.current(), resumed.current());
        assert_eq!(uninterrupted.checkpoint().geometry().nodes(), resumed.checkpoint().geometry().nodes());
        assert_eq!(uninterrupted.checkpoint().ell().to_bits(), resumed.checkpoint().ell().to_bits());
        assert_eq!(uninterrupted.checkpoint().next_iteration(), resumed.checkpoint().next_iteration());
    }
}
