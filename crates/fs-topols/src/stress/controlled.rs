//! Interruptible mechanics and cell-by-cell stress sampling.
use super::*;
use crate::evaluated::DesignEvaluationStage;
use crate::evaluated::controlled::solve_fields_controlled;
use std::ops::ControlFlow;

// Both public evaluators use this exact sampler, preserving probe order and
// first-maximum tie handling. Only a completed residual-admitted field enters.
pub(super) fn sample<B>(
    phi: &GridSdf,
    grid: &Quadtree,
    nodal: &NodalField,
    lambda: f64,
    mu: f64,
    control: &mut impl FnMut(DesignEvaluationStage) -> ControlFlow<B>,
) -> Result<ControlFlow<B, (f64, [f64; 2], usize)>, CutFemError> {
    let mut sampled_max = 0.0_f64;
    let mut max_location = [0.0, 0.0];
    let mut sample_count = 0usize;
    for (index, cell) in grid.leaves().enumerate() {
        if let ControlFlow::Break(reason) = control(DesignEvaluationStage::StressCell(index)) {
            return Ok(ControlFlow::Break(reason));
        }
        let (lo, hi) = grid.rect(cell);
        let enclosure = phi.enclose(lo, hi);
        if enclosure.lo() > 0.0 {
            continue;
        }
        if enclosure.hi() < 0.0 {
            for point in full_cell_points(lo, hi) {
                observe_stress(
                    grid, nodal, lambda, mu, point,
                    &mut sampled_max, &mut max_location, &mut sample_count,
                )?;
            }
        } else {
            let rules = cut_cell_rules(phi, lo, hi, 2);
            for &(point, weight) in &rules.bulk {
                if weight > 0.0 {
                    observe_stress(
                        grid, nodal, lambda, mu, point,
                        &mut sampled_max, &mut max_location, &mut sample_count,
                    )?;
                }
            }
            for &(point, weight, _) in &rules.iface {
                if weight > 0.0 {
                    observe_stress(
                        grid, nodal, lambda, mu, point,
                        &mut sampled_max, &mut max_location, &mut sample_count,
                    )?;
                }
            }
        }
    }
    if sample_count == 0 {
        return Err(invalid("stress evaluation found no material stress samples"));
    }
    Ok(ControlFlow::Continue((sampled_max, max_location, sample_count)))
}

/// Independently solve and sample stress with cooperative interruption.
///
/// The sample set and maximum are identical to [`evaluate_sampled_stress`].
/// Polling occurs inside the existing CG correction solver and before every
/// cell's probes. Cancellation even at `Publish` discards the complete local
/// result; no partial maximum can masquerade as an admitted stress constraint.
/// Assembly, area quadrature and an individual cell's probes are indivisible.
///
/// # Errors
/// Refuses zero polling intervals, malformed inputs or failed mechanics/stress.
pub fn evaluate_sampled_stress_controlled<B>(
    phi: &GridSdf,
    fixture: Cantilever,
    settings: OptimizeSettings,
    poll_iters: usize,
    mut control: impl FnMut(DesignEvaluationStage) -> ControlFlow<B>,
) -> Result<ControlFlow<B, SampledStressEvaluation>, CutFemError> {
    let fields = match solve_fields_controlled(phi, fixture, settings, poll_iters, &mut control)? {
        ControlFlow::Continue(fields) => fields,
        ControlFlow::Break(reason) => return Ok(ControlFlow::Break(reason)),
    };
    let (sampled_max, max_location, sample_count) = match sample(
        phi, &fields.grid, fields.solution.nodal(), fields.lambda, fields.mu, &mut control,
    )? {
        ControlFlow::Continue(samples) => samples,
        ControlFlow::Break(reason) => return Ok(ControlFlow::Break(reason)),
    };
    let result = SampledStressEvaluation {
        compliance: fields.state.compliance,
        volume: fields.state.volume,
        sampled_max_von_mises: sampled_max,
        max_location, sample_count, snapshot: fields.state.snapshot,
    };
    if let ControlFlow::Break(reason) = control(DesignEvaluationStage::Publish) {
        return Ok(ControlFlow::Break(reason));
    }
    Ok(ControlFlow::Continue(result))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::evaluated::{evaluate_compliance_design, evaluate_compliance_design_controlled};

    fn fixture() -> Cantilever { Cantilever { load: 1.0, band: 0.125 } }
    fn settings() -> OptimizeSettings { OptimizeSettings { level: 3, ..OptimizeSettings::default() } }
    fn beam() -> GridSdf { GridSdf::from_fn(8, &|_, y| (y - 0.5).abs() - 0.35) }

    #[test]
    fn projected_final_evaluation_batches_match_uncontrolled_mechanics_and_stress_bits() {
        for phi in [beam(), GridSdf::from_fn(8, &|_, _| -1.0)] {
            let ordinary = evaluate_sampled_stress(&phi, fixture(), settings()).unwrap();
            let compliance = evaluate_compliance_design(&phi, fixture(), settings()).unwrap();
            for poll_iters in [1, 7, 32, usize::MAX] {
                let mut counts = Vec::new();
                let result = evaluate_sampled_stress_controlled(&phi, fixture(), settings(), poll_iters, |stage| {
                    if let DesignEvaluationStage::Solve(count) = stage { counts.push(count); }
                    ControlFlow::<()>::Continue(())
                }).unwrap();
                let ControlFlow::Continue(result) = result else { panic!("unexpected interruption") };
                assert_eq!(result.snapshot, ordinary.snapshot);
                assert_eq!(result.sample_count, ordinary.sample_count);
                for (a, b) in [result.compliance, result.volume, result.sampled_max_von_mises,
                    result.max_location[0], result.max_location[1]].into_iter().zip([
                    ordinary.compliance, ordinary.volume, ordinary.sampled_max_von_mises,
                    ordinary.max_location[0], ordinary.max_location[1]])
                { assert_eq!(a.to_bits(), b.to_bits()); }
                assert!(counts.last().copied().unwrap() > 0);
                assert!(counts.windows(2).all(|pair| pair[0] <= pair[1]
                    && pair[1] - pair[0] <= poll_iters));
                let result = evaluate_compliance_design_controlled(&phi, fixture(), settings(), poll_iters,
                    |_| ControlFlow::<()>::Continue(())).unwrap();
                let ControlFlow::Continue(result) = result else { panic!("unexpected interruption") };
                assert_eq!(result.compliance.to_bits(), compliance.compliance.to_bits());
                assert_eq!(result.volume.to_bits(), compliance.volume.to_bits());
                assert_eq!(result.snapshot, compliance.snapshot);
            }
        }
    }

    #[test]
    fn projected_stress_interruptions_never_return_partial_evidence_and_retry_replays() {
        let phi = beam();
        let before: Vec<_> = phi.nodes().iter().map(|v| v.to_bits()).collect();
        let expected = evaluate_sampled_stress(&phi, fixture(), settings()).unwrap();
        for stop in [DesignEvaluationStage::Prepare, DesignEvaluationStage::Assemble,
            DesignEvaluationStage::Solve(1), DesignEvaluationStage::Area,
            DesignEvaluationStage::StressCell(32), DesignEvaluationStage::Publish]
        {
            let result = evaluate_sampled_stress_controlled(&phi, fixture(), settings(), 1, |stage| {
                if stage == stop { ControlFlow::Break(format!("stopped at {stop:?}")) }
                else { ControlFlow::Continue(()) }
            }).unwrap();
            assert!(matches!(result, ControlFlow::Break(ref why) if *why == format!("stopped at {stop:?}")));
            assert_eq!(phi.nodes().iter().map(|v| v.to_bits()).collect::<Vec<_>>(), before);
            let retry = evaluate_sampled_stress_controlled(&phi, fixture(), settings(), 7,
                |_| ControlFlow::<()>::Continue(())).unwrap();
            let ControlFlow::Continue(retry) = retry else { panic!("retry interrupted") };
            assert_eq!(retry, expected);
        }
        for stop in [DesignEvaluationStage::Solve(1), DesignEvaluationStage::Publish] {
            let result = evaluate_compliance_design_controlled(&phi, fixture(), settings(), 1, |stage| {
                if stage == stop { ControlFlow::Break("cancel") } else { ControlFlow::Continue(()) }
            }).unwrap();
            assert!(matches!(result, ControlFlow::Break("cancel")));
        }
    }

    #[test]
    fn projected_evaluation_invalid_controls_start_no_work_and_errors_are_not_cancellation() {
        let phi = beam();
        assert!(evaluate_sampled_stress_controlled(&phi, fixture(), settings(), 0,
            |_| -> ControlFlow<()> { panic!("zero poll interval must fail before work") }).is_err());
        assert!(evaluate_compliance_design_controlled(&phi, fixture(), settings(), 0,
            |_| -> ControlFlow<()> { panic!("zero poll interval must fail before work") }).is_err());
        let bad = OptimizeSettings { youngs: f64::NAN, ..settings() };
        assert!(matches!(evaluate_sampled_stress_controlled(&phi, fixture(), bad, 1,
            |_| ControlFlow::<()>::Continue(())), Err(CutFemError::InvalidElasticityInput { .. })));
    }
}
