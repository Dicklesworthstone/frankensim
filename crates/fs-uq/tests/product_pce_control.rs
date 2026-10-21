use fs_uq::{CorrelationModel, ParameterUncertainty, PceConfig, PceControlError, PceExecution,
    PceMeanControl, PropagationMethod, UqExecution, UqPlan, UqStatus};

fn plan() -> UqPlan {
    UqPlan::new("temperature", PropagationMethod::PolynomialChaos, 128)
        .with_parameter(ParameterUncertainty::uniform("power", 4.0, 6.0, "W"))
        .with_parameter(ParameterUncertainty::gaussian("resistance", 2.0, 0.1, "K/W"))
        .with_correlation(CorrelationModel::Independent)
        .with_compliance_threshold(310.0)
}
fn temperature(x: &[f64]) -> Result<f64, &'static str> { Ok(300.0 + x[0] * x[1]) }
fn fitted(plan: &UqPlan, f: impl FnMut(&[f64]) -> Result<f64, &'static str>) -> PceExecution {
    let mut fit = PceExecution::new(plan, PceConfig::new(2, 96)).unwrap();
    assert_eq!(fit.advance(128, || false, f).status, UqStatus::Complete);
    fit
}
fn production(plan: &UqPlan) -> UqExecution {
    let mut plan = plan.clone();
    plan.method = PropagationMethod::MonteCarlo;
    plan.seed ^= 0x5043_454d_4541_4e31;
    plan.budget_max_samples = 256;
    UqExecution::new(&plan).unwrap()
}

#[test]
fn mixed_thermal_control_estimates_physics_without_changing_raw_compliance() {
    let plan = plan();
    let fit = fitted(&plan, temperature);
    let mut raw = production(&plan);
    let control = PceMeanControl::new(&raw, &fit, 1.0).unwrap();
    assert!(control.assess(&raw).unwrap().is_none());
    let report = raw.advance(256, || false, temperature);
    let before = raw.observations().to_vec();
    let estimate = control.assess(&raw).unwrap().unwrap();
    assert_eq!(estimate.samples, 256);
    assert!((estimate.mean - 310.0).abs() < 1e-10);
    assert!(estimate.standard_error.unwrap() < 1e-11);
    assert!(estimate.variance_ratio.unwrap() < 1e-18);
    assert!(estimate.raw_standard_error.unwrap() > 0.01);
    assert_eq!(raw.observations(), before);
    assert_eq!(raw.report(), report);
    assert_eq!(control.plan(), raw.plan());
    assert_eq!(control.fit_seed(), fit.plan().seed);
    assert_eq!(control.weight(), 1.0);
}

#[test]
fn imperfect_and_variance_increasing_controls_do_not_replace_the_physical_mean() {
    let plan = plan();
    let biased = fitted(&plan, |x| Ok(400.0 + 0.5 * x[0] * x[1]));
    let exact = fitted(&plan, temperature);
    let mut raw = production(&plan);
    let imperfect = PceMeanControl::new(&raw, &biased, 1.0).unwrap();
    let harmful = PceMeanControl::new(&raw, &exact, -3.0).unwrap();
    let zero = PceMeanControl::new(&raw, &exact, 0.0).unwrap();
    raw.advance(256, || false, temperature);
    let estimate = imperfect.assess(&raw).unwrap().unwrap();
    assert!((estimate.mean - (155.0 + 0.5 * estimate.raw_mean)).abs() < 1e-10);
    assert!((estimate.variance_ratio.unwrap() - 0.25).abs() < 1e-9);
    assert!((estimate.mean - biased.model().unwrap().mean()).abs() > 90.0);
    let worse = harmful.assess(&raw).unwrap().unwrap();
    assert!((worse.variance_ratio.unwrap() - 16.0).abs() < 1e-8);
    let unchanged = zero.assess(&raw).unwrap().unwrap();
    assert_eq!(unchanged.mean, unchanged.raw_mean);
    assert_eq!(unchanged.standard_error, unchanged.raw_standard_error);
    assert_eq!(unchanged.variance_ratio, Some(1.0));
}

#[test]
fn singular_gaussian_controls_replay_original_germs_with_constant_parameters() {
    for rho in [-1.0, 1.0] {
        let plan = UqPlan::new("reaction", PropagationMethod::PolynomialChaos, 128)
            .with_parameter(ParameterUncertainty::gaussian("a", 1.0, 2.0, "1"))
            .with_parameter(ParameterUncertainty::gaussian("b", 2.0, 3.0, "1"))
            .with_parameter(ParameterUncertainty::gaussian("constant", 7.0, 0.0, "1"))
            .with_correlation(CorrelationModel::JointGaussian {
                matrix: vec![vec![1.0, rho, 0.0], vec![rho, 1.0, 0.0], vec![0.0, 0.0, 1.0]],
            });
        let physical = |x: &[f64]| Ok::<_, &'static str>(x[0] * x[1] + x[2]);
        let fit = fitted(&plan, physical);
        let mut raw = production(&plan);
        let control = PceMeanControl::new(&raw, &fit, 1.0).unwrap();
        raw.advance(256, || false, physical);
        let estimate = control.assess(&raw).unwrap().unwrap();
        assert!((estimate.mean - (9.0 + 6.0 * rho)).abs() < 1e-10);
        assert!(estimate.standard_error.unwrap() < 1e-10);
    }
}

#[test]
fn binding_and_freeze_rules_refuse_without_altering_physics() {
    let plan = plan();
    let fit = fitted(&plan, temperature);
    let mut raw = production(&plan);
    let incomplete = PceExecution::new(&plan, PceConfig::new(2, 96)).unwrap();
    assert!(matches!(PceMeanControl::new(&raw, &incomplete, 1.0), Err(PceControlError::IncompleteFit)));
    for weight in [f64::NAN, f64::INFINITY, -f64::INFINITY] {
        assert!(matches!(PceMeanControl::new(&raw, &fit, weight), Err(PceControlError::InvalidWeight)));
    }
    let mut reuse = raw.plan().clone();
    reuse.seed = fit.plan().seed;
    assert!(matches!(PceMeanControl::new(&UqExecution::new(&reuse).unwrap(), &fit, 1.0), Err(PceControlError::ReusedSeed)));
    let mut changed = raw.plan().clone();
    changed.parameters[0].unit = "kW".into();
    assert!(matches!(PceMeanControl::new(&UqExecution::new(&changed).unwrap(), &fit, 1.0), Err(PceControlError::LawMismatch)));
    let control = PceMeanControl::new(&raw, &fit, 1.0).unwrap();
    let mut changed = raw.plan().clone();
    changed.budget_max_samples -= 1;
    assert_eq!(control.assess(&UqExecution::new(&changed).unwrap()), Err(PceControlError::PlanMismatch));
    raw.advance(1, || false, temperature);
    assert!(matches!(PceMeanControl::new(&raw, &fit, 1.0), Err(PceControlError::AlreadySampled)));
    let singleton = control.assess(&raw).unwrap().unwrap();
    assert_eq!(singleton.samples, 1);
    assert_eq!(singleton.standard_error, None);
    assert_eq!(singleton.variance_ratio, None);
    raw.advance(1, || false, |_| Err::<f64, _>("physical solve failed"));
    assert_eq!(control.assess(&raw), Err(PceControlError::RefusedExecution));
    assert_eq!(raw.observations().len(), 1);
}

#[test]
fn assessment_cancellation_and_chunked_physics_preserve_paid_observations() {
    let plan = plan();
    let fit = fitted(&plan, temperature);
    let mut raw = production(&plan);
    let control = PceMeanControl::new(&raw, &fit, 1.0).unwrap();
    let mut chunks = raw.clone();
    raw.advance(256, || false, temperature);
    for count in [1, 7, 31, 300] { chunks.advance(count, || false, temperature); }
    assert_eq!(chunks.observations(), raw.observations());
    let before = raw.report();
    let mut polls = 0;
    let expected = control.assess_interruptible(&raw, || { polls += 1; false }).unwrap();
    for stop in [1, 2, polls / 2, polls - 1, polls] {
        let mut seen = 0;
        assert_eq!(control.assess_interruptible(&raw, || { seen += 1; seen == stop }), Err(PceControlError::Cancelled));
        assert_eq!(raw.report(), before);
    }
    assert_eq!(control.clone().assess(&chunks).unwrap(), expected);
    assert_eq!(control.assess(&raw).unwrap(), expected);
}
