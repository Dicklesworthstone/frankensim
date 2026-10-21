//! Physical Gaussian propagation, independent regression diagnostics and resume.
use fs_uq::{
    CorrelationModel, ParameterUncertainty, PceConfig, PceExecution, PropagationMethod,
    UncertaintyKind, UqExecution, UqPlan, UqStatus,
};

fn thermal_plan(samples: usize) -> UqPlan {
    let mut plan = UqPlan::new(
        "junction temperature",
        PropagationMethod::PolynomialChaos,
        samples,
    )
    .with_parameter(ParameterUncertainty::gaussian("power", 5.0, 0.3, "W"))
    .with_parameter(ParameterUncertainty::gaussian(
        "resistance",
        2.0,
        0.1,
        "K/W",
    ))
    .with_correlation(CorrelationModel::Independent)
    .with_compliance_threshold(310.0);
    plan.seed = 83;
    plan
}

fn temperature(x: &[f64]) -> Result<f64, &'static str> {
    Ok(300.0 + x[0] * x[1])
}

#[test]
fn thermal_polynomial_recovers_physical_moments_and_usable_surrogate() {
    let plan = thermal_plan(128);
    let config = PceConfig::new(2, 96);
    let mut run = PceExecution::new(&plan, config).unwrap();
    let mut inputs = Vec::new();
    let report = run.advance(
        usize::MAX,
        || false,
        |x| {
            inputs.push(x.to_vec());
            temperature(x)
        },
    );
    assert_eq!(report.status, UqStatus::Complete);
    assert_eq!(report.samples_evaluated, 128);
    assert_eq!(report.training_samples_accepted, 96);
    assert_eq!(report.validation_samples_accepted, 32);
    assert_eq!(report.basis_terms, 6);
    assert!(report.r_condition_estimate < config.max_condition_estimate);
    assert_eq!(run.plan(), &plan);
    for (input, germ) in inputs.iter().zip(run.germs()) {
        assert_eq!(input, &run.physical_parameters(germ).unwrap());
    }
    let summary = report.summary.unwrap();
    let variance = 4.0 * 0.09 + 25.0 * 0.01 + 0.09 * 0.01;
    assert!((summary.surrogate_mean - 310.0).abs() < 1e-10);
    assert!((summary.surrogate_variance.unwrap() - variance).abs() < 1e-10);
    assert!(summary.training_rmse < 1e-10 && summary.validation_rmse < 1e-10);
    let accepted = run.observations()[96..]
        .iter()
        .filter(|&&t| t <= 310.0)
        .count();
    let probability = accepted as f64 / 32.0;
    assert_eq!(summary.probability_of_compliance, Some(probability));
    let error = fs_math::det::sqrt(probability * (1.0 - probability) / 31.0);
    assert!((summary.compliance_standard_error.unwrap() - error).abs() < 1e-14);
    let model = run.model().unwrap();
    let germ = [0.3, -0.5];
    let x = run.physical_parameters(&germ).unwrap();
    let (prediction, derivative) = model.eval_with_gradient(&germ);
    assert!((prediction - temperature(&x).unwrap()).abs() < 1e-10);
    assert!((derivative[0] - 0.3 * x[1]).abs() < 1e-10);
    assert!((derivative[1] - 0.1 * x[0]).abs() < 1e-10);
    let sensitivity = model.sobol_indices().unwrap();
    assert!((sensitivity.first_order[0] - 0.36 / variance).abs() < 1e-10);
    assert!((sensitivity.total_order[0] - 0.3609 / variance).abs() < 1e-10);
    assert!(matches!(report.evidence_color,
        fs_evidence::Color::Estimated { dispersion, .. } if dispersion == f64::INFINITY));
}

#[test]
fn joint_gaussian_and_singular_laws_retain_mc_physical_sample_identity() {
    for rho in [0.75, 1.0] {
        let plan = UqPlan::new("sum", PropagationMethod::PolynomialChaos, 32)
            .with_parameter(ParameterUncertainty::gaussian("a", 2.0, 3.0, "m"))
            .with_parameter(ParameterUncertainty::gaussian("b", -1.0, 2.0, "m"))
            .with_correlation(CorrelationModel::JointGaussian {
                matrix: vec![vec![1.0, rho], vec![rho, 1.0]],
            });
        let mut ordinary = plan.clone();
        ordinary.method = PropagationMethod::MonteCarlo;
        let mut mc = UqExecution::new(&ordinary).unwrap();
        let mut reference = Vec::new();
        mc.advance(
            32,
            || false,
            |x| {
                reference.push(x.to_vec());
                Ok::<_, &str>(x[0] + x[1])
            },
        );
        let mut pce = PceExecution::new(&plan, PceConfig::new(1, 24)).unwrap();
        let mut actual = Vec::new();
        let report = pce.advance(
            32,
            || false,
            |x| {
                actual.push(x.to_vec());
                Ok::<_, &str>(x[0] + x[1])
            },
        );
        assert_eq!(report.status, UqStatus::Complete);
        assert_eq!(actual, reference);
        assert_eq!(pce.observations(), mc.observations());
        let summary = report.summary.unwrap();
        assert!((summary.surrogate_mean - 1.0).abs() < 1e-12);
        assert!((summary.surrogate_variance.unwrap() - (13.0 + 12.0 * rho)).abs() < 1e-10);
        assert!(summary.validation_rmse < 1e-12);
    }
}

#[test]
fn held_out_physics_changes_diagnostics_without_leaking_into_coefficients() {
    let plan = thermal_plan(32);
    let config = PceConfig::new(2, 24);
    let mut original = PceExecution::new(&plan, config).unwrap();
    original.advance(32, || false, temperature);
    let mut shifted = PceExecution::new(&plan, config).unwrap();
    let training_inputs: Vec<_> = shifted.germs()[..24]
        .iter()
        .map(|germ| shifted.physical_parameters(germ).unwrap())
        .collect();
    // A deliberately discontinuous, but fixed deterministic physical function:
    // it agrees at every training input and differs by 10 K everywhere else.
    let report = shifted.advance(
        32,
        || false,
        |x| {
            let offset = if training_inputs.iter().any(|p| p.as_slice() == x) {
                0.0
            } else {
                10.0
            };
            Ok::<_, &str>(temperature(x)? + offset)
        },
    );
    assert_eq!(original.model(), shifted.model());
    let summary = report.summary.unwrap();
    assert!(summary.training_rmse < 1e-10);
    assert!((summary.validation_rmse - 10.0).abs() < 1e-10);
    assert!((summary.validation_mean_error - 10.0).abs() < 1e-10);
    assert_eq!(summary.probability_of_compliance, Some(0.0));
}

#[test]
fn smooth_cooling_ode_has_analytic_gaussian_moments_and_measured_holdout_error() {
    // Exact cooling solution T(1)/T(0)=exp(-k), k~N(1, 0.15^2).
    let plan =
        UqPlan::new("temperature ratio", PropagationMethod::PolynomialChaos, 640).with_parameter(
            ParameterUncertainty::gaussian("cooling rate", 1.0, 0.15, "1/s"),
        );
    let mut run = PceExecution::new(&plan, PceConfig::new(3, 512)).unwrap();
    let report = run.advance(640, || false, |x| Ok::<_, &str>(fs_math::det::exp(-x[0])));
    assert_eq!(report.status, UqStatus::Complete);
    let summary = report.summary.unwrap();
    let mean = fs_math::det::exp(-1.0 + 0.0225 / 2.0);
    let variance = fs_math::det::exp(-2.0 + 0.0225) * (fs_math::det::exp(0.0225) - 1.0);
    assert!((summary.surrogate_mean - mean).abs() < 3e-4);
    assert!((summary.surrogate_variance.unwrap() - variance).abs() < 3e-5);
    assert!(summary.validation_rmse > 0.0 && summary.validation_rmse < 3e-4);
    assert!(summary.probability_of_compliance.is_none());
}

#[test]
fn cancellation_interruption_and_clone_preserve_paid_work_and_fit_identity() {
    let plan = thermal_plan(32);
    let config = PceConfig::new(2, 24);
    let mut full = PceExecution::new(&plan, config).unwrap();
    let mut expected_inputs = Vec::new();
    let expected = full.advance(
        32,
        || false,
        |x| {
            expected_inputs.push(x.to_vec());
            temperature(x)
        },
    );
    let mut split = PceExecution::new(&plan, config).unwrap();
    let mut inputs = Vec::new();
    let first = split.advance(
        5,
        || false,
        |x| {
            inputs.push(x.to_vec());
            temperature(x)
        },
    );
    assert!(first.summary.is_none() && split.model().is_none());
    assert_eq!(
        split.advance(
            0,
            || panic!("zero allowance"),
            |_| -> Result<f64, &str> { panic!("zero allowance") }
        ),
        first
    );
    let paused = split.advance(1, || true, |_| -> Result<f64, &str> { panic!("cancelled") });
    assert_eq!(paused.status, UqStatus::Cancelled);
    let interrupted = split.advance_interruptible(
        1,
        || false,
        |x| {
            assert_eq!(x, expected_inputs[5].as_slice());
            Ok::<_, &str>(None)
        },
    );
    assert_eq!(interrupted.samples_evaluated, 5);
    let mut resumed = split.clone();
    for allowance in [3, 7, 99] {
        resumed.advance(
            allowance,
            || false,
            |x| {
                inputs.push(x.to_vec());
                temperature(x)
            },
        );
    }
    assert_eq!(inputs, expected_inputs);
    assert_eq!(resumed.report(), expected);
    assert_eq!(resumed.model(), full.model());
    assert_eq!(split.observations().len(), 5);
    // Cancellation immediately before or after the fit cannot repay a model call.
    for stop_after in [32, 33] {
        let mut run = PceExecution::new(&plan, config).unwrap();
        let mut polls = 0;
        let paused = run.advance(
            32,
            || {
                polls += 1;
                polls > stop_after
            },
            temperature,
        );
        assert_eq!(paused.status, UqStatus::Cancelled);
        assert_eq!(paused.samples_evaluated, 32);
        assert!(paused.summary.is_none() && run.model().is_none());
        let resumed = run.advance(
            1,
            || false,
            |_| -> Result<f64, &str> { panic!("already paid") },
        );
        assert_eq!(resumed, expected);
        assert_eq!(run.model(), full.model());
        assert_eq!(
            run.advance(
                1,
                || panic!("terminal"),
                |_| -> Result<f64, &str> { panic!("terminal") }
            ),
            expected
        );
    }
}

#[test]
fn unsupported_laws_oversampling_storage_and_conditioning_refuse_before_work() {
    let plan = thermal_plan(32);
    for config in [
        PceConfig::new(2, 11),
        PceConfig::new(2, 31),
        PceConfig::new(9, 24),
        PceConfig::new(2, usize::MAX),
    ] {
        assert!(PceExecution::new(&plan, config).is_err());
    }
    for limit in [f64::NAN, f64::INFINITY, 0.0, 1e9, 1.0] {
        let config = PceConfig {
            max_condition_estimate: limit,
            ..PceConfig::new(2, 24)
        };
        assert!(PceExecution::new(&plan, config).is_err());
    }
    let mut invalid = plan.clone();
    invalid.correlation = CorrelationModel::Unknown;
    assert!(PceExecution::new(&invalid, PceConfig::new(2, 24)).is_err());
    invalid.correlation = CorrelationModel::Independent;
    for kind in [
        UncertaintyKind::AleatoryUniform { lo: 2.0, hi: 1.0 },
        UncertaintyKind::EpistemicInterval { lo: 1.0, hi: 2.0 },
        UncertaintyKind::AleatoryGaussian {
            mean: 0.0,
            std_dev: -1.0,
        },
    ] {
        invalid.parameters[0].kind = kind;
        assert!(PceExecution::new(&invalid, PceConfig::new(2, 24)).is_err());
    }
    let mut wrong = plan.clone();
    wrong.method = PropagationMethod::MonteCarlo;
    assert!(PceExecution::new(&wrong, PceConfig::new(2, 24)).is_err());
    let mut huge = plan.clone();
    huge.budget_max_samples = 65_537;
    assert!(PceExecution::new(&huge, PceConfig::new(2, 24)).is_err());
    let mut wide = plan;
    for i in 2..17 {
        wide.parameters.push(ParameterUncertainty::gaussian(
            format!("p{i}"),
            0.0,
            1.0,
            "1",
        ));
    }
    assert!(PceExecution::new(&wide, PceConfig::new(0, 24)).is_err());
    wide.parameters.pop();
    wide.budget_max_samples = 1024;
    assert!(PceExecution::new(&wide, PceConfig::new(2, 512)).is_err()); // 153 terms.
    let mut expensive = thermal_plan(65_536);
    expensive
        .parameters
        .push(ParameterUncertainty::gaussian("third", 0.0, 1.0, "1"));
    assert!(PceExecution::new(&expensive, PceConfig::new(4, 65_534)).is_err());
}

#[test]
fn zero_uncertainty_is_constant_and_failed_physics_is_terminal() {
    let plan = UqPlan::new("constant", PropagationMethod::PolynomialChaos, 4)
        .with_parameter(ParameterUncertainty::gaussian("fixed", 7.0, 0.0, "K"));
    let mut run = PceExecution::new(&plan, PceConfig::new(0, 2)).unwrap();
    let report = run.advance(
        4,
        || false,
        |x| {
            assert_eq!(x, [7.0]);
            Ok::<_, &str>(x[0])
        },
    );
    let summary = report.summary.unwrap();
    assert!((summary.surrogate_mean - 7.0).abs() < 1e-13);
    assert_eq!(summary.surrogate_variance, Some(0.0));
    assert_eq!(summary.surrogate_std_dev, Some(0.0));
    assert!(run.physical_parameters(&[]).is_err());
    assert!(run.physical_parameters(&[f64::NAN]).is_err());
    for output in [Err("solver failed"), Ok(f64::NAN), Ok(f64::INFINITY)] {
        let mut run = PceExecution::new(&thermal_plan(32), PceConfig::new(2, 24)).unwrap();
        run.advance(2, || false, temperature);
        let failed = run.advance(1, || false, |_| output);
        assert_eq!(failed.status, UqStatus::Refused);
        assert_eq!(failed.samples_evaluated, 3);
        assert_eq!(run.observations().len(), 2);
        assert!(failed.summary.is_none() && run.model().is_none());
        assert_eq!(
            run.advance(
                1,
                || panic!("terminal"),
                |_| -> Result<f64, &str> { panic!("terminal") }
            ),
            failed
        );
    }
}
