//! Cancellable stress admission and two-solve checkpoint recovery.
use super::*;
use crate::evaluated::DesignEvaluationStage;
use crate::projected::ProjectedSetupStage;

/// Recovery distinguishes area-feasible staging from full stress admission.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProjectedStressSetupStage {
    /// Area/fixed-node verification and independent compliance replay.
    Area(ProjectedSetupStage),
    /// Independent stress solve and deterministic cell sampling.
    Stress(DesignEvaluationStage),
    /// Both constraints passed; no optimizer has yet been returned.
    Publish,
}

impl ProjectedStressOptimizer {
    /// Admit a stress limit without consuming the caller's area-feasible state.
    ///
    /// Cancellation reaches CG and every cell's sampling. A partial sampled
    /// maximum is never returned, and the input remains available for retry.
    ///
    /// # Errors
    /// Refuses invalid limits, failed mechanics, or an overstressed baseline.
    pub fn new_controlled<B>(
        optimizer: &ProjectedOptimizer,
        limit: SampledStressLimit,
        mut control: impl FnMut(ProjectedStressSetupStage) -> ControlFlow<B>,
    ) -> Result<ControlFlow<B, Self>, CutFemError> {
        let limit = SampledStressLimit::new(limit.max_von_mises, limit.absolute_tolerance)?;
        let checkpoint = optimizer.checkpoint();
        let current = match evaluate_sampled_stress_controlled(
            checkpoint.geometry(), checkpoint.fixture(), checkpoint.settings(),
            optimizer.controls().poll_iters,
            |stage| control(ProjectedStressSetupStage::Stress(stage)),
        )? {
            ControlFlow::Continue(current) => current,
            ControlFlow::Break(reason) => return Ok(ControlFlow::Break(reason)),
        };
        if let Some(reason) = admission_reason(&current, optimizer.current(), limit) {
            return Err(refused(format!("stress-constrained baseline refused: {reason}")));
        }
        if let ControlFlow::Break(reason) = control(ProjectedStressSetupStage::Publish) {
            return Ok(ControlFlow::Break(reason));
        }
        Ok(ControlFlow::Continue(Self {
            optimizer: optimizer.clone(), limit, baseline: current.clone(), current,
        }))
    }

    /// Reconstruct a saved stress-constrained state with cancellable PDE replay.
    ///
    /// Preserves the exact geometry, original iteration target and AL multiplier.
    /// Area and stress are independently re-solved; no previous update is run.
    /// `Break` is a stop reason, not a feasible baseline or a constraint refusal.
    ///
    /// # Errors
    /// Refuses invalid policy, changed fixed nodes, infeasible area/stress or PDE failure.
    pub fn from_checkpoint_controlled<B>(
        checkpoint: &OptimizeCheckpoint,
        fixed: Vec<(usize, f64)>,
        projection: VolumeProjectionSettings,
        controls: ProjectedSettings,
        limit: SampledStressLimit,
        mut control: impl FnMut(ProjectedStressSetupStage) -> ControlFlow<B>,
    ) -> Result<ControlFlow<B, Self>, CutFemError> {
        let limit = SampledStressLimit::new(limit.max_von_mises, limit.absolute_tolerance)?;
        let optimizer = match ProjectedOptimizer::from_checkpoint_controlled(
            checkpoint, fixed, projection, controls,
            |stage| control(ProjectedStressSetupStage::Area(stage)),
        )? {
            ControlFlow::Continue(optimizer) => optimizer,
            ControlFlow::Break(reason) => return Ok(ControlFlow::Break(reason)),
        };
        Self::new_controlled(&optimizer, limit, control)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Cantilever, GridSdf, OptimizeSettings};

    fn input() -> GridSdf { GridSdf::from_fn(8, &|_, y| (y - 0.5).abs() - 0.35) }
    fn fixture() -> Cantilever { Cantilever { load: 1.0, band: 0.125 } }
    fn settings() -> OptimizeSettings {
        OptimizeSettings { level: 3, iterations: 2, volfrac: 0.6, move_cells: 0.1,
            nucleation_period: 0, ..OptimizeSettings::default() }
    }
    fn controls() -> ProjectedSettings {
        ProjectedSettings { max_candidates: 8, poll_iters: 1, ..ProjectedSettings::default() }
    }
    fn projection() -> VolumeProjectionSettings {
        VolumeProjectionSettings { target: 0.6, tolerance: 1e-4, max_shift: 2.0, max_evaluations: 64 }
    }
    fn fixed(phi: &GridSdf) -> Vec<(usize, f64)> {
        phi.nodes().iter().copied().enumerate().filter(|(i, _)| i % 9 == 0 || i % 9 == 8).collect()
    }
    fn bits(phi: &GridSdf) -> Vec<u64> { phi.nodes().iter().map(|v| v.to_bits()).collect() }
    fn limit() -> SampledStressLimit { SampledStressLimit::new(1e12, 0.0).unwrap() }
    fn source() -> ProjectedOptimizer {
        let phi = input();
        ProjectedOptimizer::new(phi.clone(), fixture(), settings(), fixed(&phi), projection(), controls()).unwrap()
    }
    fn same(a: &ProjectedStressOptimizer, b: &ProjectedStressOptimizer) {
        assert_eq!(bits(a.checkpoint().geometry()), bits(b.checkpoint().geometry()));
        assert_eq!(a.checkpoint().ell().to_bits(), b.checkpoint().ell().to_bits());
        assert_eq!(a.checkpoint().next_iteration(), b.checkpoint().next_iteration());
        assert_eq!(a.current(), b.current());
    }

    #[test]
    fn final_solve_and_stress_cancellation_preserve_accepted_state_and_retry_history() {
        let before = ProjectedStressOptimizer::new(source(), limit()).unwrap();
        let mut uninterrupted = before.clone();
        let expected = uninterrupted.advance_one().unwrap();
        assert!(matches!(&expected.progress, ProjectedProgress::Accepted(_)));
        for stop in 0..5 {
            let mut optimizer = before.clone();
            let result = optimizer.advance_one_controlled(|stage| {
                let stop_here = match stop {
                    0 => matches!(stage, ProjectedStage::Evaluation(_, DesignEvaluationStage::Solve(n)) if n > 0),
                    1 => matches!(stage, ProjectedStage::Stress(_, DesignEvaluationStage::Solve(n)) if n > 0),
                    2 => matches!(stage, ProjectedStage::Stress(_, DesignEvaluationStage::StressCell(32))),
                    3 => matches!(stage, ProjectedStage::Stress(_, DesignEvaluationStage::Publish)),
                    _ => matches!(stage, ProjectedStage::Publish(_)),
                };
                if stop_here { ControlFlow::Break(format!("stop {stop}")) }
                else { ControlFlow::Continue(()) }
            }).unwrap();
            assert!(matches!(result, ControlFlow::Break(ref why) if *why == format!("stop {stop}")));
            same(&optimizer, &before);
            let retry = optimizer.advance_one().unwrap();
            same(&optimizer, &uninterrupted);
            assert_eq!(format!("{retry:?}"), format!("{expected:?}"));
        }
    }

    #[test]
    fn baseline_setup_is_interruptible_and_matches_ordinary_feasible_geometry() {
        let phi = input();
        let original = bits(&phi);
        for stop in 0..4 {
            let result = ProjectedOptimizer::new_controlled(&phi, fixture(), settings(),
                fixed(&phi), projection(), controls(), |stage| {
                    let stop_here = match stop {
                        0 => matches!(stage, ProjectedSetupStage::Prepare),
                        1 => matches!(stage, ProjectedSetupStage::Projection(_)),
                        2 => matches!(stage, ProjectedSetupStage::Evaluation(DesignEvaluationStage::Solve(n)) if n > 0),
                        _ => matches!(stage, ProjectedSetupStage::Publish),
                    };
                    if stop_here { ControlFlow::Break(stop) } else { ControlFlow::Continue(()) }
                }).unwrap();
            assert!(matches!(result, ControlFlow::Break(reason) if reason == stop));
            assert_eq!(bits(&phi), original);
        }
        let ControlFlow::Continue(area) = ProjectedOptimizer::new_controlled(&phi, fixture(), settings(),
            fixed(&phi), projection(), controls(), |_| ControlFlow::<()>::Continue(())).unwrap()
        else { panic!("setup interrupted") };
        let ordinary = source();
        assert_eq!(bits(area.checkpoint().geometry()), bits(ordinary.checkpoint().geometry()));
        assert_eq!(area.current().compliance.to_bits(), ordinary.current().compliance.to_bits());
        let ControlFlow::Continue(stress) = ProjectedStressOptimizer::new_controlled(&area, limit(),
            |_| ControlFlow::<()>::Continue(())).unwrap() else { panic!("stress setup interrupted") };
        same(&stress, &ProjectedStressOptimizer::new(ordinary, limit()).unwrap());
    }

    #[test]
    fn recovery_can_stop_in_either_solve_or_sampling_without_losing_the_source() {
        let mut original = ProjectedStressOptimizer::new(source(), limit()).unwrap();
        assert!(matches!(original.advance_one().unwrap().progress, ProjectedProgress::Accepted(_)));
        let checkpoint = original.checkpoint();
        let original_bits = bits(checkpoint.geometry());
        for stop in 0..5 {
            let result = ProjectedStressOptimizer::from_checkpoint_controlled(checkpoint,
                fixed(checkpoint.geometry()), projection(), controls(), limit(), |stage| {
                    let stop_here = match stop {
                        0 => matches!(stage, ProjectedStressSetupStage::Area(ProjectedSetupStage::Prepare)),
                        1 => matches!(stage, ProjectedStressSetupStage::Area(ProjectedSetupStage::Evaluation(DesignEvaluationStage::Solve(n))) if n > 0),
                        2 => matches!(stage, ProjectedStressSetupStage::Stress(DesignEvaluationStage::Solve(n)) if n > 0),
                        3 => matches!(stage, ProjectedStressSetupStage::Stress(DesignEvaluationStage::StressCell(32))),
                        _ => matches!(stage, ProjectedStressSetupStage::Publish),
                    };
                    if stop_here { ControlFlow::Break(stop) } else { ControlFlow::Continue(()) }
                }).unwrap();
            assert!(matches!(result, ControlFlow::Break(reason) if reason == stop));
            assert_eq!(bits(checkpoint.geometry()), original_bits);
            assert_eq!(checkpoint.next_iteration(), 1);
        }
        let ControlFlow::Continue(mut resumed) = ProjectedStressOptimizer::from_checkpoint_controlled(
            checkpoint, fixed(checkpoint.geometry()), projection(), controls(), limit(),
            |_| ControlFlow::<()>::Continue(())).unwrap() else { panic!("recovery interrupted") };
        same(&resumed, &original);
        let a = resumed.advance_one().unwrap();
        let b = original.advance_one().unwrap();
        same(&resumed, &original);
        assert_eq!(format!("{a:?}"), format!("{b:?}"));
    }

    #[test]
    fn invalid_limits_and_constraint_cancellation_cannot_start_another_candidate() {
        let mut area = source();
        let before = bits(area.checkpoint().geometry());
        let invalid = SampledStressLimit { max_von_mises: f64::NAN, absolute_tolerance: 0.0 };
        assert!(ProjectedStressOptimizer::from_checkpoint_controlled(area.checkpoint(),
            fixed(area.checkpoint().geometry()), projection(), controls(), invalid,
            |_| -> ControlFlow<()> { panic!("invalid limit must refuse before replay") }).is_err());
        let tiny = SampledStressLimit::new(1e-30, 0.0).unwrap();
        assert!(ProjectedStressOptimizer::new_controlled(&area, tiny,
            |_| ControlFlow::<()>::Continue(())).is_err());
        let mut checked = 0;
        let result = area.advance_one_admitted_controlled(|_, _, _, _| {
            checked += 1;
            Ok(ControlFlow::Break("constraint stop"))
        }, |_| ControlFlow::Continue(())).unwrap();
        assert!(matches!(result, ControlFlow::Break("constraint stop")));
        assert_eq!(checked, 1);
        assert_eq!(bits(area.checkpoint().geometry()), before);
        assert_eq!(area.checkpoint().next_iteration(), 0);
    }
}
