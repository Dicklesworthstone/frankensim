//! Physical rare-event probability, likelihood accounting and resumable work.
use fs_eproc::GaussianMixtureCs;
use fs_uq::{
    CorrelationModel, ImportanceConfig, ImportanceExecution, ParameterUncertainty,
    PropagationMethod, UncertaintyKind, UqComplianceError, UqExecution, UqPlan, UqStatus,
};

fn plan(n: usize) -> UqPlan {
    let mut plan = UqPlan::new("response", PropagationMethod::MonteCarlo, n)
        .with_parameter(ParameterUncertainty::gaussian("load", 0.0, 1.0, "N"))
        .with_correlation(CorrelationModel::Independent)
        .with_compliance_threshold(5.0);
    plan.seed = 83;
    plan
}

fn config(germ_shift: Vec<f64>, defensive_fraction: f64) -> ImportanceConfig {
    ImportanceConfig {
        germ_shift,
        defensive_fraction,
    }
}

fn response(x: &[f64]) -> Result<f64, &'static str> {
    Ok(x[0])
}

#[test]
fn five_sigma_failure_probability_is_resolved_by_physical_proposal_calls() {
    let mut run = ImportanceExecution::new(&plan(8192), config(vec![5.0], 0.1)).unwrap();
    let report = run.advance(usize::MAX, || false, response);
    assert_eq!(report.status, UqStatus::Complete);
    assert_eq!(report.evaluations_attempted, 8192);
    assert!(report.proposal_failures > 3000);
    let estimate = report.estimate.unwrap();
    // Analytic standard-normal upper tail at five, not a fitted surrogate.
    let target = 2.866515718791933e-7;
    assert!((estimate.failure_probability / target - 1.0).abs() < 0.15);
    assert!(estimate.standard_error.unwrap() < target * 0.06);
    assert!(estimate.failure_effective_samples.unwrap() > 500.0);
    let confidence = run.assess_failure_probability(0.05, 1e-6).unwrap().unwrap();
    assert!(!confidence.converged);
    assert!(confidence.hi > estimate.failure_probability * 1000.0);
    assert!(matches!(
        report.evidence_color,
        fs_evidence::Color::Estimated { .. }
    ));
}

#[test]
fn physical_units_correlations_and_singular_laws_preserve_the_target_event() {
    for rho in [0.6_f64, 1.0] {
        let mut physical = plan(8192);
        physical.parameters = vec![
            ParameterUncertainty::gaussian("load", 100.0, 3.0, "N"),
            ParameterUncertainty::gaussian("modulus", 20.0, 2.0, "Pa"),
        ];
        physical.correlation = CorrelationModel::JointGaussian {
            matrix: vec![vec![1.0, rho], vec![rho, 1.0]],
        };
        let denominator = (2.0 + 2.0 * rho).sqrt();
        let shift = vec![
            5.0 * (1.0 + rho) / denominator,
            5.0 * (1.0 - rho * rho).sqrt() / denominator,
        ];
        let mut run = ImportanceExecution::new(&physical, config(shift, 0.1)).unwrap();
        let report = run.advance(
            usize::MAX,
            || false,
            |x| {
                let a = (x[0] - 100.0) / 3.0;
                let b = (x[1] - 20.0) / 2.0;
                if rho == 1.0 {
                    assert!((a - b).abs() < 1e-13);
                }
                Ok::<_, &str>((a + b) / denominator)
            },
        );
        assert_eq!(report.status, UqStatus::Complete);
        assert!(
            (report.estimate.unwrap().failure_probability / 2.866515718791933e-7 - 1.0).abs()
                < 0.15
        );
        assert_eq!(run.plan(), &physical);
    }
}

#[test]
fn likelihood_moments_and_confidence_match_independent_reconstruction() {
    let physical = plan(256).with_compliance_threshold(1.0);
    let mut run = ImportanceExecution::new(&physical, config(vec![2.0], 0.25)).unwrap();
    let estimate = run.advance(256, || false, response).estimate.unwrap();
    let mut weighted = Vec::new();
    let mut sum_weights = 0.0;
    let mut cs = GaussianMixtureCs::new(0.5, 1.0, 0.05);
    for observation in run.observations() {
        let reference = 1.0 / (0.25 + 0.75 * (2.0 * observation.qoi - 2.0).exp());
        assert!((observation.likelihood_weight / reference - 1.0).abs() < 1e-13);
        assert_eq!(
            observation.scaled_weight / 0.25,
            observation.likelihood_weight
        );
        assert!(observation.scaled_weight > 0.0 && observation.scaled_weight <= 1.0);
        let value = if observation.qoi > 1.0 {
            reference
        } else {
            0.0
        };
        weighted.push(value);
        sum_weights += reference;
        cs.observe(if observation.failed {
            observation.scaled_weight
        } else {
            0.0
        });
    }
    let n = weighted.len() as f64;
    let sum: f64 = weighted.iter().sum();
    let mean = sum / n;
    let variance = weighted.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (n - 1.0);
    let ess = sum * sum / weighted.iter().map(|x| x * x).sum::<f64>();
    assert!((estimate.failure_probability - mean).abs() < 1e-14);
    assert!((estimate.standard_error.unwrap() - (variance / n).sqrt()).abs() < 1e-14);
    assert!((estimate.failure_effective_samples.unwrap() - ess).abs() < 1e-10);
    assert!((estimate.mean_likelihood_weight - sum_weights / n).abs() < 1e-14);
    let (center, radius) = cs.interval().unwrap();
    let assessed = run.assess_failure_probability(0.05, 0.02).unwrap().unwrap();
    assert_eq!(assessed.mean, center / 0.25);
    assert_eq!(assessed.lo, ((center - radius) / 0.25).max(0.0));
    assert_eq!(assessed.hi, ((center + radius) / 0.25).min(1.0));
    assert_eq!(assessed.converged, radius / 0.25 <= 0.02);
}

#[test]
fn zero_shift_and_full_defensive_mass_preserve_ordinary_mc_bits() {
    let mut physical = plan(32).with_compliance_threshold(0.0);
    physical
        .parameters
        .push(ParameterUncertainty::gaussian("other", 2.0, 3.0, "Pa"));
    physical.correlation = CorrelationModel::JointGaussian {
        matrix: vec![vec![1.0, -0.6], vec![-0.6, 1.0]],
    };
    let mut mc = UqExecution::new(&physical).unwrap();
    let mut expected = Vec::new();
    mc.advance(
        32,
        || false,
        |x| {
            expected.push(x.to_vec());
            response(x)
        },
    );
    for proposal in [config(vec![0.0, 0.0], 0.1), config(vec![5.0, -2.0], 1.0)] {
        let mut run = ImportanceExecution::new(&physical, proposal.clone()).unwrap();
        let mut actual = Vec::new();
        run.advance(
            32,
            || false,
            |x| {
                actual.push(x.to_vec());
                response(x)
            },
        );
        assert_eq!(actual, expected);
        assert_eq!(run.config(), &proposal);
        let mut ordinary_cs = GaussianMixtureCs::new(0.5, 1.0, 0.05);
        for (observation, &qoi) in run.observations().iter().zip(mc.observations()) {
            assert_eq!(observation.qoi.to_bits(), qoi.to_bits());
            assert_eq!(observation.likelihood_weight, 1.0);
            assert_eq!(observation.scaled_weight, proposal.defensive_fraction);
            ordinary_cs.observe(if observation.failed { 1.0 } else { 0.0 });
        }
        let assessed = run.assess_failure_probability(0.05, 0.01).unwrap().unwrap();
        let (center, radius) = ordinary_cs.interval().unwrap();
        assert_eq!(assessed.mean, center);
        assert_eq!(assessed.lo, (center - radius).max(0.0));
        assert_eq!(assessed.hi, (center + radius).min(1.0));
    }
}

#[test]
fn finite_sample_weighted_probability_is_not_self_normalized_or_clipped() {
    let mut saw_above_one = false;
    for seed in 0..32 {
        let mut physical = plan(2);
        physical.seed = seed;
        let mut run = ImportanceExecution::new(&physical, config(vec![2.0], 0.5)).unwrap();
        let estimate = run
            .advance(2, || false, |_| Ok::<_, &str>(6.0))
            .estimate
            .unwrap();
        let raw = run
            .observations()
            .iter()
            .map(|v| v.likelihood_weight)
            .sum::<f64>()
            / 2.0;
        assert!((estimate.failure_probability - raw).abs() < 1e-14);
        saw_above_one |= estimate.failure_probability > 1.0;
    }
    assert!(saw_above_one);
}

#[test]
fn chunks_cancellation_unfinished_calls_and_clone_retain_paid_work() {
    let physical = plan(32);
    let proposal = config(vec![5.0], 0.1);
    let mut full = ImportanceExecution::new(&physical, proposal.clone()).unwrap();
    let mut expected_inputs = Vec::new();
    let expected = full.advance(
        99,
        || false,
        |x| {
            expected_inputs.push(x.to_vec());
            response(x)
        },
    );
    let mut split = ImportanceExecution::new(&physical, proposal).unwrap();
    assert!(
        split
            .assess_failure_probability(0.05, 0.01)
            .unwrap()
            .is_none()
    );
    let mut paid_inputs = Vec::new();
    let partial = split.advance(
        3,
        || false,
        |x| {
            paid_inputs.push(x.to_vec());
            response(x)
        },
    );
    assert_eq!(partial.evaluations_accepted, 3);
    assert!(partial.estimate.is_some());
    assert_eq!(
        split.advance(0, || panic!("zero allowance"), response),
        partial
    );
    assert_eq!(
        split.advance(1, || true, response).status,
        UqStatus::Cancelled
    );
    let paused = split.advance_interruptible(
        1,
        || false,
        |x| {
            assert_eq!(x, expected_inputs[3].as_slice());
            Ok::<_, &str>(None)
        },
    );
    assert_eq!(paused.evaluations_attempted, 3);
    let mut resumed = split.clone();
    for count in [1, 4, 100] {
        resumed.advance(
            count,
            || false,
            |x| {
                paid_inputs.push(x.to_vec());
                response(x)
            },
        );
    }
    assert_eq!(paid_inputs, expected_inputs);
    assert_eq!(resumed.observations(), full.observations());
    assert_eq!(resumed.report(), expected);
    assert_eq!(
        resumed.assess_failure_probability(0.05, 0.01),
        full.assess_failure_probability(0.05, 0.01)
    );
    assert_eq!(
        resumed.advance(1, || panic!("complete"), response),
        expected
    );
}

#[test]
fn no_observed_failures_does_not_prove_zero_target_risk() {
    let mut run = ImportanceExecution::new(&plan(16), config(vec![2.0], 0.1)).unwrap();
    let estimate = run
        .advance(16, || false, |_| Ok::<_, &str>(5.0))
        .estimate
        .unwrap();
    assert_eq!(estimate.failure_probability, 0.0); // Equality is compliant.
    assert_eq!(estimate.standard_error, Some(0.0));
    assert_eq!(estimate.failure_effective_samples, None);
    let assessed = run.assess_failure_probability(0.05, 0.01).unwrap().unwrap();
    assert!(assessed.hi > 0.0);
    assert!(!assessed.converged);
}

#[test]
fn invalid_laws_thresholds_and_unrepresentable_proposals_refuse_at_admission() {
    let original = plan(16);
    for epsilon in [0.0, -0.1, 1.1, f64::NAN, f64::INFINITY, f64::from_bits(1)] {
        assert!(ImportanceExecution::new(&original, config(vec![5.0], epsilon)).is_err());
    }
    for shift in [vec![], vec![1.0, 2.0], vec![f64::NAN], vec![f64::MAX]] {
        assert!(ImportanceExecution::new(&original, config(shift, 0.1)).is_err());
    }
    let mut invalid = original.clone();
    invalid.compliance_threshold = None;
    assert!(ImportanceExecution::new(&invalid, config(vec![5.0], 0.1)).is_err());
    invalid = original.clone();
    invalid.correlation = CorrelationModel::Unknown;
    assert!(ImportanceExecution::new(&invalid, config(vec![5.0], 0.1)).is_err());
    invalid = original;
    invalid.parameters[0].kind = UncertaintyKind::AleatoryUniform { lo: 0.0, hi: 1.0 };
    assert!(ImportanceExecution::new(&invalid, config(vec![5.0], 0.1)).is_err());
}

#[test]
fn numerical_and_model_failures_are_terminal_and_never_publish_selected_prefixes() {
    // A representable configuration whose shifted draw has unrepresentable
    // weight: it must refuse BEFORE spending a physical-model call.
    let mut underflow = ImportanceExecution::new(&plan(16), config(vec![100.0], 1e-300)).unwrap();
    let failed = underflow.advance(
        1,
        || false,
        |_| -> Result<f64, &str> { panic!("invalid likelihood reached physical evaluator") },
    );
    assert_eq!(failed.status, UqStatus::Refused);
    assert_eq!(failed.evaluations_attempted, 0);
    assert!(failed.rejection_reason.unwrap().contains("weight"));
    for output in [Err("solver failed"), Ok(f64::NAN)] {
        let mut run = ImportanceExecution::new(&plan(16), config(vec![5.0], 0.1)).unwrap();
        run.advance(2, || false, response);
        let failed = run.advance(1, || false, |_| output);
        assert_eq!(failed.status, UqStatus::Refused);
        assert_eq!(failed.evaluations_attempted, 3);
        assert_eq!(failed.evaluations_accepted, 2);
        assert!(failed.estimate.is_none());
        assert_eq!(
            run.assess_failure_probability(0.05, 0.01),
            Err(UqComplianceError::RefusedExecution)
        );
        assert_eq!(run.advance(1, || panic!("terminal"), response), failed);
    }
    let run = ImportanceExecution::new(&plan(16), config(vec![5.0], 0.1)).unwrap();
    assert_eq!(
        run.assess_failure_probability(0.0, 0.01),
        Err(UqComplianceError::InvalidAlpha)
    );
    assert_eq!(
        run.assess_failure_probability(0.05, -1.0),
        Err(UqComplianceError::InvalidHalfWidth)
    );
}
