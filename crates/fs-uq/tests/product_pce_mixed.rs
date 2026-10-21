//! Uniform tolerances and Gaussian scatter use one product execution owner.
use fs_uq::{
    CorrelationModel, ParameterUncertainty, PceBasis, PceConfig, PceExecution,
    PropagationMethod, UqExecution, UqPlan, UqStatus,
};

fn plan(samples: usize) -> UqPlan {
    UqPlan::new("junction temperature", PropagationMethod::PolynomialChaos, samples)
        .with_parameter(ParameterUncertainty::uniform("power", 4.0, 6.0, "W"))
        .with_parameter(ParameterUncertainty::gaussian("resistance", 2.0, 0.1, "K/W"))
        .with_correlation(CorrelationModel::Independent)
        .with_compliance_threshold(310.0)
}
fn thermal(x: &[f64]) -> Result<f64, &'static str> { Ok(300.0 + x[0]*x[1]) }
fn close(a: f64, b: f64) {
    assert!((a-b).abs() < 1e-10*(1.0+a.abs().max(b.abs())), "{a} != {b}");
}

#[test]
fn mixed_thermal_model_recovers_moments_and_physical_input_sensitivity() {
    let mut run = PceExecution::new(&plan(128), PceConfig::new(2,96)).unwrap();
    assert_eq!(run.basis(), &[PceBasis::Uniform01, PceBasis::StandardNormal]);
    let report = run.advance(128, || false, thermal);
    assert_eq!(report.status, UqStatus::Complete);
    let summary = report.summary.unwrap();
    let vp = 1.0/3.0;
    let vr = 0.01;
    let variance = 4.0*vp + 25.0*vr + vp*vr;
    close(summary.surrogate_mean, 310.0);
    close(summary.surrogate_variance.unwrap(), variance);
    assert!(summary.training_rmse < 1e-10 && summary.validation_rmse < 1e-10);
    let model = run.model().unwrap();
    assert_eq!(model.basis(), run.basis());
    let sobol = model.sobol_indices().unwrap();
    close(sobol.first_order[0], 4.0*vp/variance);
    close(sobol.total_order[0], (4.0+vr)*vp/variance);
    close(sobol.first_order[1], 25.0*vr/variance);
    close(sobol.total_order[1], (25.0+vp)*vr/variance);
    for germ in [[0.0,0.0], [0.5,-0.5], [1.0,2.0]] {
        let physical = run.physical_parameters(&germ).unwrap();
        let (value, gradient) = model.eval_with_gradient(&germ);
        close(value, thermal(&physical).unwrap());
        close(gradient[0], 2.0*physical[1]);
        close(gradient[1], 0.1*physical[0]);
    }
    let passed = run.observations()[96..].iter().filter(|&&x| x <= 310.0).count();
    assert_eq!(summary.probability_of_compliance, Some(passed as f64/32.0));
    assert!(matches!(report.evidence_color,
        fs_evidence::Color::Estimated { dispersion, .. } if dispersion == f64::INFINITY));
}

#[test]
fn uniform_quadratic_retains_interior_extremum_and_analytic_moments() {
    let plan = UqPlan::new("response", PropagationMethod::PolynomialChaos, 48)
        .with_parameter(ParameterUncertainty::uniform("position", -1.0, 1.0, "m"));
    let mut run = PceExecution::new(&plan, PceConfig::new(2,32)).unwrap();
    let report = run.advance(48, || false, |x| Ok::<_, &str>(1.0-x[0]*x[0]));
    assert_eq!(report.status, UqStatus::Complete);
    let summary = report.summary.unwrap();
    close(summary.surrogate_mean, 2.0/3.0);
    close(summary.surrogate_variance.unwrap(), 4.0/45.0);
    let model = run.model().unwrap();
    for (u, value, derivative) in [(0.0,0.0,4.0), (0.5,1.0,0.0), (1.0,0.0,-4.0)] {
        let (actual, gradient) = model.eval_with_gradient(&[u]);
        close(actual, value);
        close(gradient[0], derivative);
    }
}

#[test]
fn mixed_draws_match_canonical_mc_in_both_orders_with_constant_parameters() {
    for reverse in [false,true] {
        let mut plan = plan(32).with_parameter(ParameterUncertainty::uniform("fixed", 7.0, 7.0, "K"));
        if reverse { plan.parameters.reverse(); }
        let mut ordinary = plan.clone();
        ordinary.method = PropagationMethod::MonteCarlo;
        let mut mc = UqExecution::new(&ordinary).unwrap();
        let mut reference = Vec::new();
        mc.advance(32, || false, |x| {
            reference.push(x.to_vec());
            Ok::<_, &str>(x.iter().sum::<f64>())
        });
        let mut run = PceExecution::new(&plan, PceConfig::new(1,24)).unwrap();
        let mut actual = Vec::new();
        let report = run.advance(32, || false, |x| {
            actual.push(x.to_vec());
            Ok::<_, &str>(x.iter().sum::<f64>())
        });
        assert_eq!(report.status, UqStatus::Complete);
        assert_eq!(actual, reference);
        assert_eq!(run.observations(), mc.observations());
        for (x, germ) in actual.iter().zip(run.germs()) {
            assert_eq!(x, &run.physical_parameters(germ).unwrap());
        }
    }
}

#[test]
fn mixed_paid_work_and_fit_survive_interruption_cancellation_and_clone() {
    let plan = plan(32);
    let config = PceConfig::new(2,24);
    let mut full = PceExecution::new(&plan, config).unwrap();
    let mut expected_inputs = Vec::new();
    let expected = full.advance(32, || false, |x| {
        expected_inputs.push(x.to_vec());
        thermal(x)
    });
    let mut split = PceExecution::new(&plan, config).unwrap();
    let mut paid = Vec::new();
    let first = split.advance(7, || false, |x| { paid.push(x.to_vec()); thermal(x) });
    assert!(first.summary.is_none() && split.model().is_none());
    let paused = split.advance_interruptible(1, || false, |x| {
        assert_eq!(x, expected_inputs[7].as_slice());
        Ok::<_, &str>(None)
    });
    assert_eq!(paused.samples_evaluated, 7);
    assert_eq!(split.advance(1, || true, |_| -> Result<f64,&str> { panic!("cancelled") }).samples_evaluated, 7);
    let mut resumed = split.clone();
    for allowance in [1,8,64] {
        resumed.advance(allowance, || false, |x| { paid.push(x.to_vec()); thermal(x) });
    }
    assert_eq!(paid, expected_inputs);
    assert_eq!(resumed.report(), expected);
    assert_eq!(resumed.model(), full.model());
    assert_eq!(split.observations().len(), 7);
    for boundary in [32,33] {
        let mut run = PceExecution::new(&plan, config).unwrap();
        let mut polls = 0;
        let paused = run.advance(32, || { polls += 1; polls > boundary }, thermal);
        assert_eq!(paused.status, UqStatus::Cancelled);
        assert!(run.model().is_none());
        assert_eq!(run.advance(1, || false, |_| -> Result<f64,&str> { panic!("already paid") }), expected);
        assert_eq!(run.model(), full.model());
    }
}

#[test]
fn mixed_holdout_data_do_not_change_the_fitted_polynomial() {
    let config = PceConfig::new(2,24);
    let mut baseline = PceExecution::new(&plan(32), config).unwrap();
    baseline.advance(32, || false, thermal);
    let mut shifted = PceExecution::new(&plan(32), config).unwrap();
    let training: Vec<_> = shifted.germs()[..24].iter()
        .map(|g| shifted.physical_parameters(g).unwrap()).collect();
    let report = shifted.advance(32, || false, |x| {
        let offset = if training.iter().any(|p| p.as_slice() == x) { 0.0 } else { 100.0 };
        Ok::<_, &str>(thermal(x)?+offset)
    });
    assert_eq!(report.status, UqStatus::Complete);
    assert_eq!(baseline.model(), shifted.model());
    let summary = report.summary.unwrap();
    assert!(summary.training_rmse < 1e-10);
    close(summary.validation_rmse, 100.0);
    close(summary.validation_mean_error, 100.0);
    assert_eq!(summary.probability_of_compliance, Some(0.0));
}

#[test]
fn uniforms_do_not_erase_dependence_refusals_or_germ_domain_checks() {
    let mut input = plan(32);
    for dependence in [
        CorrelationModel::Unknown,
        CorrelationModel::Correlated { matrix: vec![vec![1.0,0.5],vec![0.5,1.0]] },
        CorrelationModel::JointGaussian { matrix: vec![vec![1.0,0.5],vec![0.5,1.0]] },
    ] {
        input.correlation = dependence;
        assert!(PceExecution::new(&input, PceConfig::new(2,24)).is_err());
    }
    let input = UqPlan::new("extreme", PropagationMethod::PolynomialChaos, 4)
        .with_parameter(ParameterUncertainty::uniform("span", -f64::MAX, f64::MAX, "m"));
    let run = PceExecution::new(&input, PceConfig::new(0,2)).unwrap();
    assert_eq!(run.physical_parameters(&[0.0]).unwrap(), [-f64::MAX]);
    assert_eq!(run.physical_parameters(&[0.5]).unwrap(), [0.0]);
    assert_eq!(run.physical_parameters(&[1.0]).unwrap(), [f64::MAX]);
    for u in [-0.001,1.001,f64::NAN,f64::INFINITY] {
        assert!(run.physical_parameters(&[u]).is_err());
    }
    let fixed = UqPlan::new("fixed", PropagationMethod::PolynomialChaos, 4)
        .with_parameter(ParameterUncertainty::uniform("fixed", 7.0, 7.0, "K"));
    let mut run = PceExecution::new(&fixed, PceConfig::new(0,2)).unwrap();
    let report = run.advance(4, || false, |x| { assert_eq!(x,[7.0]); Ok::<_, &str>(x[0]) });
    assert_eq!(report.status, UqStatus::Complete);
    close(report.summary.unwrap().surrogate_mean, 7.0);
    assert!(run.physical_parameters(&[-1.0]).is_err());
}
