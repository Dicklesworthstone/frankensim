use fs_uq::{
    CorrelationModel, MlmcExecution, ParameterUncertainty, PropagationMethod, UqExecution, UqPlan,
    UqStatus,
};

fn plan(allocation: &[usize]) -> UqPlan {
    let calls = allocation[0] + 2 * allocation[1..].iter().sum::<usize>();
    let mut plan = UqPlan::new(
        "temperature ratio",
        PropagationMethod::MultilevelMonteCarlo,
        calls,
    )
    .with_parameter(ParameterUncertainty::uniform(
        "cooling rate",
        1.0,
        2.0,
        "1/s",
    ));
    plan.seed = 83;
    plan
}

// Backward Euler for dT/dt = -k*T, T(0)=1, over one second. Every model level
// sees the same physical uncertain k; refinement changes only the time grid.
fn cooling(level: usize, parameters: &[f64]) -> Result<f64, &'static str> {
    let steps = 1usize << level;
    let factor = 1.0 + parameters[0] / steps as f64;
    let mut temperature = 1.0;
    for _ in 0..steps {
        temperature /= factor;
    }
    Ok(temperature)
}

fn variance_of_mean(values: &[f64]) -> f64 {
    let n = values.len() as f64;
    let mean = values.iter().sum::<f64>() / n;
    values.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (n * (n - 1.0))
}

#[test]
fn physical_refinement_ladder_owns_pairing_and_telescoping_statistics() {
    let allocation = [4096, 1024, 256, 128];
    let p = plan(&allocation);
    let mut execution = MlmcExecution::new(&p, &allocation).unwrap();
    let mut calls = Vec::new();
    let report = execution.advance(
        usize::MAX,
        || false,
        |level, parameters| {
            let value = cooling(level, parameters)?;
            calls.push((level, parameters.to_vec(), value));
            Ok::<_, &str>(value)
        },
    );
    assert_eq!(report.status, UqStatus::Complete);
    assert_eq!(report.model_calls_evaluated, p.budget_max_samples);
    assert_eq!(report.samples_accepted, allocation.iter().sum::<usize>());
    assert!(!report.pending_pair);
    let mut cursor = 0usize;
    let mut reconstructed_mean = 0.0;
    let mut reconstructed_variance = 0.0;
    for (level, &count) in allocation.iter().enumerate() {
        let mut corrections = Vec::new();
        for _ in 0..count {
            let fine = &calls[cursor];
            assert_eq!(fine.0, level);
            cursor += 1;
            let correction = if level == 0 {
                fine.2
            } else {
                let coarse = &calls[cursor];
                assert_eq!(coarse.0, level - 1);
                assert_eq!(fine.1, coarse.1, "coupling uses identical physical inputs");
                cursor += 1;
                fine.2 - coarse.2
            };
            corrections.push(correction);
        }
        assert_eq!(
            execution.corrections(level).unwrap(),
            corrections.as_slice()
        );
        let mean = corrections.iter().sum::<f64>() / count as f64;
        let variance = variance_of_mean(&corrections);
        let observed = report.levels[level].correction.as_ref().unwrap();
        assert!((observed.mean - mean).abs() < 2e-14);
        assert!((observed.estimator_variance.unwrap() - variance).abs() < 1e-16);
        reconstructed_mean += mean;
        reconstructed_variance += variance;
    }
    let estimate = report.estimate.unwrap();
    assert!((estimate.mean - reconstructed_mean).abs() < 2e-14);
    assert!((estimate.estimator_variance.unwrap() - reconstructed_variance).abs() < 1e-16);
    // Exact expectation of the finest discrete solution for uniform k in [1,2].
    let expected = 8.0 / 7.0 * ((1.0_f64 + 1.0 / 8.0).powi(-7) - (1.0_f64 + 2.0 / 8.0).powi(-7));
    assert!((estimate.mean - expected).abs() < 6.0 * estimate.standard_error.unwrap());
    assert!(matches!(
        report.evidence_color,
        fs_evidence::Color::Estimated { .. }
    ));
    assert_eq!(
        execution.plan().method,
        PropagationMethod::MultilevelMonteCarlo
    );
    assert_eq!(execution.allocation(), allocation);
}

#[test]
fn one_level_reuses_mc_samples_and_singular_joint_gaussian_admission() {
    let allocation = [256];
    let p = UqPlan::new("sum", PropagationMethod::MultilevelMonteCarlo, 256)
        .with_parameter(ParameterUncertainty::gaussian("a", 2.0, 3.0, "m"))
        .with_parameter(ParameterUncertainty::gaussian("b", -1.0, 2.0, "m"))
        .with_correlation(CorrelationModel::JointGaussian {
            matrix: vec![vec![1.0, 1.0], vec![1.0, 1.0]],
        })
        .with_compliance_threshold(1.0);
    let mut ordinary = p.clone();
    ordinary.method = PropagationMethod::MonteCarlo;
    let mut mc = UqExecution::new(&ordinary).unwrap();
    let reference = mc.advance(256, || false, |x| Ok::<_, &str>(x[0] + x[1]));
    let mut mlmc = MlmcExecution::new(&p, &allocation).unwrap();
    let report = mlmc.advance(
        256,
        || false,
        |level, x| {
            assert_eq!(level, 0);
            assert!(((x[0] - 2.0) / 3.0 - (x[1] + 1.0) / 2.0).abs() < 1e-13);
            Ok::<_, &str>(x[0] + x[1])
        },
    );
    assert_eq!(mlmc.corrections(0).unwrap(), mc.observations());
    let estimate = report.estimate.unwrap();
    assert_eq!(estimate.mean, reference.mean.unwrap());
    assert!((estimate.standard_error.unwrap() - reference.sampling_error).abs() < 1e-14);
    assert_eq!(
        report.compliance.unwrap().mean,
        reference.probability_of_compliance.unwrap()
    );
}

#[test]
fn a_paid_fine_half_survives_chunks_cancellation_interruption_and_clone() {
    let allocation = [2, 2];
    let p = plan(&allocation);
    let mut full = MlmcExecution::new(&p, &allocation).unwrap();
    let mut full_calls = Vec::new();
    let expected = full.advance(
        99,
        || false,
        |level, x| {
            full_calls.push((level, x.to_vec()));
            cooling(level, x)
        },
    );
    let mut split = MlmcExecution::new(&p, &allocation).unwrap();
    let mut paid_calls = Vec::new();
    let first = split.advance(
        1,
        || false,
        |level, x| {
            paid_calls.push((level, x.to_vec()));
            cooling(level, x)
        },
    );
    assert!(
        first.levels[0]
            .correction
            .as_ref()
            .unwrap()
            .standard_error
            .is_none()
    );
    assert!(
        first.levels[0]
            .correction
            .as_ref()
            .unwrap()
            .estimator_variance
            .is_none()
    );
    assert!(first.estimate.is_none());
    assert!(matches!(first.evidence_color,
        fs_evidence::Color::Estimated { dispersion, .. } if dispersion == f64::INFINITY));
    let half = split.advance(
        2,
        || false,
        |level, x| {
            paid_calls.push((level, x.to_vec()));
            cooling(level, x)
        },
    );
    assert_eq!(half.model_calls_evaluated, 3);
    assert_eq!(half.samples_accepted, 2);
    assert!(half.pending_pair && half.levels[1].correction.is_none());
    assert_eq!(
        split.advance(
            0,
            || panic!("zero allowance"),
            |_, _| -> Result<f64, &str> { panic!("zero allowance") }
        ),
        half
    );
    let paused = split.advance(
        1,
        || true,
        |_, _| -> Result<f64, &str> { panic!("cancelled") },
    );
    assert_eq!(paused.status, UqStatus::Cancelled);
    assert_eq!(paused.model_calls_evaluated, 3);
    let interrupted = split.advance_interruptible(
        1,
        || false,
        |level, x| {
            assert_eq!(level, 0);
            assert_eq!(x, full_calls[3].1.as_slice());
            Ok::<_, &str>(None)
        },
    );
    assert_eq!(interrupted.model_calls_evaluated, 3);
    assert!(interrupted.pending_pair);
    let mut resumed = split.clone();
    for _ in 0..3 {
        resumed.advance(
            1,
            || false,
            |level, x| {
                paid_calls.push((level, x.to_vec()));
                cooling(level, x)
            },
        );
    }
    assert_eq!(
        paid_calls, full_calls,
        "no completed half is evaluated twice"
    );
    assert_eq!(resumed.report(), expected);
    assert!(
        split.report().pending_pair,
        "the checkpoint remains unchanged"
    );
    assert_eq!(
        resumed.advance(
            1,
            || panic!("complete is terminal"),
            |_, _| -> Result<f64, &str> { panic!("complete is terminal") }
        ),
        expected
    );
}

#[test]
fn compliance_telescopes_indicators_of_physical_levels() {
    let allocation = [2, 2];
    let p = UqPlan::new("response", PropagationMethod::MultilevelMonteCarlo, 6)
        .with_parameter(ParameterUncertainty::uniform("fixed", 0.0, 0.0, "1"))
        .with_compliance_threshold(0.5);
    let mut run = MlmcExecution::new(&p, &allocation).unwrap();
    let report = run.advance(6, || false, |level, _| Ok::<_, &str>(level as f64));
    assert_eq!(run.corrections(0).unwrap(), [0.0, 0.0]);
    assert_eq!(run.corrections(1).unwrap(), [1.0, 1.0]);
    assert_eq!(run.compliance_corrections(0).unwrap(), [1.0, 1.0]);
    assert_eq!(run.compliance_corrections(1).unwrap(), [-1.0, -1.0]);
    assert_eq!(report.estimate.unwrap().mean, 1.0);
    let probability = report.compliance.unwrap();
    assert_eq!(probability.mean, 0.0);
    assert_eq!(probability.standard_error, Some(0.0));
}

#[test]
fn invalid_allocations_and_probability_laws_refuse_before_work() {
    let p = plan(&[2, 2]);
    for allocation in [
        vec![],
        vec![1, 2],
        vec![2, 1],
        vec![2],
        vec![2; 65],
        vec![usize::MAX, 2],
    ] {
        assert!(MlmcExecution::new(&p, &allocation).is_err());
    }
    let mut wrong = p.clone();
    wrong.method = PropagationMethod::MonteCarlo;
    assert!(MlmcExecution::new(&wrong, &[2, 2]).is_err());
    let unknown = p
        .clone()
        .with_parameter(ParameterUncertainty::uniform("other", 0.0, 1.0, "1"));
    assert!(MlmcExecution::new(&unknown, &[2, 2]).is_err());
    let invalid = p
        .clone()
        .with_parameter(ParameterUncertainty::interval("interval", 0.0, 1.0, "m"));
    assert!(MlmcExecution::new(&invalid, &[2, 2]).is_err());
    let mut unnamed = p.clone();
    unnamed.parameters[0].unit.clear();
    assert!(MlmcExecution::new(&unnamed, &[2, 2]).is_err());
    let mut excessive = p;
    excessive.budget_max_samples = 1_000_001;
    assert!(MlmcExecution::new(&excessive, &[1_000_001]).is_err());
    let maximum_levels = [2; 64];
    assert!(MlmcExecution::new(&plan(&maximum_levels), &maximum_levels).is_ok());
}

#[test]
fn coarse_refusals_are_terminal_and_finite_extremes_keep_honest_dispersion() {
    let allocation = [2, 2];
    for output in [Err("physical refusal"), Ok(f64::NAN), Ok(f64::INFINITY)] {
        let mut run = MlmcExecution::new(&plan(&allocation), &allocation).unwrap();
        run.advance(3, || false, cooling);
        let report = run.advance(1, || false, |_, _| output);
        assert_eq!(report.status, UqStatus::Refused);
        assert_eq!(report.model_calls_evaluated, 4);
        assert_eq!(report.samples_accepted, 2);
        assert!(report.estimate.is_none() && report.compliance.is_none());
        assert!(report.levels.iter().all(|l| l.correction.is_none()));
        assert_eq!(
            run.advance(
                8,
                || panic!("terminal"),
                |_, _| -> Result<f64, &str> { panic!("terminal") }
            ),
            report
        );
    }
    let mut run = MlmcExecution::new(&plan(&allocation), &allocation).unwrap();
    let report = run.advance(6, || false, |_, _| Ok::<_, &str>(f64::MAX));
    assert_eq!(report.status, UqStatus::Complete);
    let estimate = report.estimate.unwrap();
    assert_eq!(estimate.mean, f64::MAX);
    assert_eq!(estimate.standard_error, Some(0.0));
    assert_eq!(estimate.estimator_variance, Some(0.0));
    let mut wide = MlmcExecution::new(&plan(&[2]), &[2]).unwrap();
    let mut ordinal = 0;
    let report = wide.advance(
        2,
        || false,
        |_, _| {
            ordinal += 1;
            Ok::<_, &str>(if ordinal == 1 { -f64::MAX } else { f64::MAX })
        },
    );
    let estimate = report.estimate.unwrap();
    assert_eq!(estimate.mean, 0.0);
    assert_eq!(estimate.standard_error, Some(f64::MAX));
    assert!(estimate.estimator_variance.is_none());
    let mut overflow = MlmcExecution::new(&plan(&allocation), &allocation).unwrap();
    let report = overflow.advance(
        6,
        || false,
        |level, _| Ok::<_, &str>(if level == 0 { -f64::MAX } else { f64::MAX }),
    );
    assert_eq!(report.status, UqStatus::Refused);
    assert_eq!(report.model_calls_evaluated, 4);
    assert!(report.estimate.is_none());
}
