//! Budget, cancellation, failure, and deterministic-resumption regressions.
//! Bead: frankensim-extreal-program-f85xj.6.7 (execution prerequisite).

use std::cell::Cell;

use fs_evidence::Color;
use fs_uq::{
    CorrelationModel, ParameterUncertainty, PropagationMethod, UqExecution, UqPlan,
    UqPropagator, UqStatus,
};

fn plan(samples: usize) -> UqPlan {
    UqPlan::new("junction_maximum", PropagationMethod::MonteCarlo, samples)
        .with_parameter(ParameterUncertainty::gaussian("ambient", 300.0, 5.0, "K"))
        .with_parameter(ParameterUncertainty::gaussian("power", 50.0, 2.0, "W"))
        .with_correlation(CorrelationModel::JointGaussian {
            matrix: vec![vec![1.0, 0.6], vec![0.6, 1.0]],
        })
        .with_compliance_threshold(345.0)
}

#[allow(clippy::unnecessary_wraps)]
fn model(parameters: &[f64]) -> Result<f64, &'static str> {
    Ok(parameters[0] + 0.8 * parameters[1])
}

#[test]
fn every_chunking_has_identical_observations_and_final_evidence() {
    let plan = plan(79);
    let reference = UqPropagator::run(&plan, |x| model(x).unwrap());
    for chunk in [1, 2, 7, 31, 79, usize::MAX] {
        let mut execution = UqExecution::new(&plan).unwrap();
        let mut calls = 0;
        while execution.report().status != UqStatus::Complete {
            let before = calls;
            let result = execution.advance(chunk, || false, |x| {
                calls += 1;
                model(x)
            });
            assert!(calls - before <= chunk);
            assert!(matches!(result.evidence_color, Color::Estimated { .. }));
            assert!(matches!(result.status, UqStatus::Complete | UqStatus::BudgetTruncated));
        }
        assert_eq!(calls, plan.budget_max_samples);
        assert_eq!(execution.evaluations_attempted(), calls);
        assert_eq!(execution.report(), reference);
        assert_eq!(execution.report().content_hash(), reference.content_hash());
        // Independently reconstruct the original Philox sample addressing and
        // 2x2 Cholesky transform, rather than comparing two wrappers alone.
        for (ordinal, &observed) in execution.observations().iter().enumerate() {
            let mut stream = fs_rand::StreamKey {
                seed: plan.seed,
                kernel: 0x0517,
                tile: ordinal as u32,
            }.stream();
            let a = stream.next_normal();
            let b = stream.next_normal();
            let expected = (300.0 + 5.0 * a)
                + 0.8 * (50.0 + 2.0 * (0.6 * a + (1.0_f64 - 0.6 * 0.6).sqrt() * b));
            assert_eq!(observed.to_bits(), expected.to_bits());
        }
    }
}

#[test]
fn cancellation_before_first_sample_returns_no_invented_statistics() {
    let mut execution = UqExecution::new(&plan(10)).unwrap();
    let result = execution.advance(10, || true, |_| -> Result<f64, &str> {
        panic!("cancelled execution invoked model")
    });
    assert_eq!(result.status, UqStatus::Cancelled);
    assert_eq!(result.samples_evaluated, 0);
    assert_eq!(result.mean, None);
    assert_eq!(result.std_dev, None);
    assert_eq!(result.percentiles, None);
    assert_eq!(result.probability_of_compliance, None);
    assert_eq!(execution.advance(10, || false, model).status, UqStatus::Complete);
}

#[test]
fn cancellation_retains_paid_work_and_resumes_the_next_ordinal() {
    let plan = plan(10);
    let calls = Cell::new(0);
    let mut execution = UqExecution::new(&plan).unwrap();
    let paused = execution.advance(10, || calls.get() == 3, |x| {
        calls.set(calls.get() + 1);
        model(x)
    });
    assert_eq!(paused.status, UqStatus::Cancelled);
    assert_eq!(paused.samples_evaluated, 3);
    let prefix = execution.observations().to_vec();
    let completed = execution.advance(usize::MAX, || false, |x| {
        calls.set(calls.get() + 1);
        model(x)
    });
    assert_eq!(calls.get(), 10);
    assert_eq!(&execution.observations()[..3], prefix.as_slice());
    assert_eq!(completed, UqPropagator::run(&plan, |x| model(x).unwrap()));
}

#[test]
fn zero_work_and_terminal_calls_do_not_invoke_callbacks() {
    let mut execution = UqExecution::new(&plan(4)).unwrap();
    let zero = execution.advance(0, || panic!("zero work polled cancellation"), model);
    assert_eq!(zero.status, UqStatus::BudgetTruncated);
    assert_eq!(zero.samples_evaluated, 0);
    let completed = execution.advance(usize::MAX, || false, model);
    assert_eq!(completed.status, UqStatus::Complete);
    assert_eq!(execution.advance(1, || panic!("completed run polled cancellation"), |_| -> Result<f64, &str> {
        panic!("completed run invoked model")
    }), completed);
}

#[test]
fn original_lifetime_budget_is_not_reset_or_mutated_by_the_caller() {
    let mut original = plan(5);
    let mut execution = UqExecution::new(&original).unwrap();
    original.budget_max_samples = 1_000_000;
    original.seed += 1;
    for expected in [2, 4, 5, 5] {
        execution.advance(2, || false, model);
        assert_eq!(execution.evaluations_attempted(), expected);
    }
    assert_eq!(execution.plan().budget_max_samples, 5);
    assert_eq!(execution.plan().seed + 1, original.seed);
}

#[test]
fn solver_failure_is_terminal_and_never_drops_an_unfavorable_sample() {
    let mut execution = UqExecution::new(&plan(20)).unwrap();
    let mut calls = 0;
    let failure = execution.advance(20, || false, |x| {
        calls += 1;
        if calls == 3 { Err("conduction residual did not converge") } else { model(x) }
    });
    assert_eq!(failure.status, UqStatus::Refused);
    assert_eq!(failure.samples_evaluated, 3);
    assert_eq!(execution.observations().len(), 2);
    assert_eq!(failure.mean, None);
    assert_eq!(failure.probability_of_compliance, None);
    assert!(failure.rejection_reason.as_ref().unwrap().contains("conduction residual"));
    assert_eq!(execution.advance(20, || false, |_| -> Result<f64, &str> {
        panic!("failed sample must not be silently skipped or retried")
    }), failure);
}

#[test]
fn nonfinite_model_outputs_fail_closed() {
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let mut execution = UqExecution::new(&plan(4)).unwrap();
        execution.advance(1, || false, model);
        let failure = execution.advance(3, || false, |_| Ok::<_, &str>(bad));
        assert_eq!(failure.status, UqStatus::Refused);
        assert_eq!(failure.samples_evaluated, 2);
        assert_eq!(execution.observations().len(), 1);
        assert_eq!(failure.std_dev, None);
        assert_eq!(execution.advance(3, || false, model), failure);
    }
}

#[test]
fn one_observation_has_no_fabricated_variance() {
    let mut execution = UqExecution::new(&plan(4)).unwrap();
    let result = execution.advance(1, || false, |_| Ok::<_, &str>(340.0));
    assert_eq!(result.status, UqStatus::BudgetTruncated);
    assert_eq!(result.mean, Some(340.0));
    assert_eq!(result.std_dev, None);
    assert_eq!(result.percentiles, Some([340.0; 3]));
    assert_eq!(result.interval_bounds, [340.0; 2]);
}

#[test]
fn unrepresentable_partial_dispersion_does_not_change_final_execution() {
    let mut execution = UqExecution::new(&plan(3)).unwrap();
    let mut outputs = [f64::MAX, -f64::MAX, 0.0].into_iter();
    let partial = execution.advance(2, || false, |_| Ok::<_, &str>(outputs.next().unwrap()));
    assert_eq!(partial.status, UqStatus::BudgetTruncated);
    assert_eq!(partial.std_dev, None);
    let completed = execution.advance(1, || false, |_| Ok::<_, &str>(outputs.next().unwrap()));
    assert_eq!(completed.status, UqStatus::Complete);
    assert_eq!(completed.mean, Some(0.0));
    assert_eq!(completed.std_dev, Some(f64::MAX));
}

#[test]
fn admission_and_terminal_statistical_overflow_are_preserved() {
    let mut invalid = plan(4);
    invalid.method = PropagationMethod::QuasiMonteCarlo;
    assert!(UqExecution::new(&invalid).is_err());
    let mut execution = UqExecution::new(&plan(2)).unwrap();
    let mut values = [f64::MAX, -f64::MAX].into_iter();
    let result = execution.advance(2, || false, |_| Ok::<_, &str>(values.next().unwrap()));
    assert_eq!(result.status, UqStatus::Refused);
    assert_eq!(execution.advance(10, || false, model), result);
}
