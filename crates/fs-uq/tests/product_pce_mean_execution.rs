use fs_uq::{CorrelationModel, ParameterUncertainty, PceConfig, PceMeanConfig, PceMeanControl,
    PceMeanExecution, PceMeanStage, PropagationMethod, UqExecution, UqPlan, UqStatus};

fn plan() -> UqPlan {
    UqPlan::new("temperature", PropagationMethod::MonteCarlo, 384)
        .with_parameter(ParameterUncertainty::uniform("power", 4.0, 6.0, "W"))
        .with_parameter(ParameterUncertainty::gaussian("resistance", 2.0, 0.1, "K/W"))
        .with_correlation(CorrelationModel::Independent)
        .with_compliance_threshold(310.0)
}
fn config() -> PceMeanConfig {
    PceMeanConfig { fit_evaluations: 128, regression: PceConfig::new(2, 96), control_weight: 1.0 }
}
fn physics(x: &[f64]) -> Result<f64, &'static str> { Ok(300.0 + x[0] * x[1]) }

#[test]
fn one_budget_separates_fit_data_and_reproduces_standalone_physical_production() {
    let plan = plan();
    let mut run = PceMeanExecution::new(&plan, config()).unwrap();
    let report = run.advance(usize::MAX, || false, physics);
    assert_eq!(report.status, UqStatus::Complete);
    assert_eq!(report.stage, PceMeanStage::Complete);
    assert_eq!(report.model_calls_evaluated, 384);
    assert_eq!(report.model_calls_planned, 384);
    assert_eq!(report.fit.samples_evaluated, 128);
    assert_eq!(report.fit.validation_samples_accepted, 32);
    assert_eq!(run.production().observations().len(), 256);
    assert_eq!(run.fit().observations().len(), 128);
    assert_eq!(report.fit_seed ^ report.production_seed, 0x5043_454d_4541_4e31);
    let mut reference = UqExecution::new(run.production().plan()).unwrap();
    let control = PceMeanControl::new(&reference, run.fit(), 1.0).unwrap();
    let raw = reference.advance(256, || false, physics);
    assert_eq!(run.production().observations(), reference.observations());
    assert_eq!(report.production, Some(raw));
    assert_eq!(report.estimate, control.assess(&reference).unwrap());
    let estimate = report.estimate.unwrap();
    assert!((estimate.mean - 310.0).abs() < 1e-10);
    assert!(estimate.variance_ratio.unwrap() < 1e-18);
    assert_eq!(run.plan(), &plan);
    assert_eq!(run.config(), config());
}

#[test]
fn chunks_interruptions_and_clone_never_repeat_paid_fitting_or_production() {
    let mut full = PceMeanExecution::new(&plan(), config()).unwrap();
    let mut expected_inputs = Vec::new();
    let expected = full.advance(usize::MAX, || false, |x| {
        expected_inputs.push(x.to_vec()); physics(x)
    });
    let mut split = PceMeanExecution::new(&plan(), config()).unwrap();
    let mut inputs = Vec::new();
    for n in [1, 17, 109, 2] {
        let before = split.evaluations_attempted();
        let report = split.advance(n, || false, |x| { inputs.push(x.to_vec()); physics(x) });
        assert!(report.model_calls_evaluated - before <= n);
        assert!(report.estimate.is_none());
    }
    assert_eq!(split.evaluations_attempted(), 129);
    let before = split.report();
    assert_eq!(split.advance(0, || panic!("zero allowance"), |_| -> Result<f64, &str> { panic!("zero allowance") }), before);
    let mut interrupted = None;
    let paused = split.advance_interruptible(1, || false, |x| {
        interrupted = Some(x.to_vec()); Ok::<_, &str>(None)
    });
    assert_eq!(paused.status, UqStatus::Cancelled);
    assert_eq!(paused.model_calls_evaluated, 129);
    let mut resumed = split.clone();
    let report = resumed.advance(usize::MAX, || false, |x| {
        if inputs.len() == 129 { assert_eq!(Some(x.to_vec()), interrupted); }
        inputs.push(x.to_vec()); physics(x)
    });
    assert_eq!(inputs, expected_inputs);
    assert_eq!(report, expected);
    assert_eq!(split.evaluations_attempted(), 129);
    assert_eq!(resumed.advance(1, || panic!("complete"), |_| -> Result<f64, &str> { panic!("complete") }), expected);
}

#[test]
fn cancellation_at_fit_transition_and_final_assessment_resumes_exactly() {
    let mut full = PceMeanExecution::new(&plan(), config()).unwrap();
    let mut expected_inputs = Vec::new();
    let mut polls = 0;
    let expected = full.advance(usize::MAX, || { polls += 1; false }, |x| {
        expected_inputs.push(x.to_vec()); physics(x)
    });
    // The first poll is owned by the study. PCE then polls before every fit
    // call and around fitting; production and final control replay follow.
    for boundary in [1, 2, 128, 129, 130, 131, 132, 133, polls - 1, polls] {
        let mut run = PceMeanExecution::new(&plan(), config()).unwrap();
        let mut seen = 0;
        let mut inputs = Vec::new();
        let paused = run.advance(usize::MAX, || { seen += 1; seen == boundary }, |x| {
            inputs.push(x.to_vec()); physics(x)
        });
        assert_eq!(paused.status, UqStatus::Cancelled, "checkpoint {boundary}");
        assert!(paused.estimate.is_none());
        assert_eq!(paused.model_calls_evaluated, inputs.len());
        if boundary == polls {
            assert_eq!(paused.stage, PceMeanStage::Assessment);
            assert_eq!(paused.model_calls_evaluated, 384);
        }
        let resumed = run.advance(usize::MAX, || false, |x| { inputs.push(x.to_vec()); physics(x) });
        assert_eq!(inputs, expected_inputs, "checkpoint {boundary}");
        assert_eq!(resumed, expected, "checkpoint {boundary}");
    }
}

#[test]
fn physical_refusals_in_either_phase_are_terminal_and_never_censor_an_observation() {
    for fail_at in [10, 129] {
        let mut run = PceMeanExecution::new(&plan(), config()).unwrap();
        let mut calls = 0;
        let refused = run.advance(usize::MAX, || false, |x| {
            calls += 1;
            if calls == fail_at { Err("physical solve refused") } else { physics(x) }
        });
        assert_eq!(refused.status, UqStatus::Refused);
        assert_eq!(refused.model_calls_evaluated, fail_at);
        assert!(refused.estimate.is_none());
        assert!(refused.rejection_reason.as_deref().unwrap().contains("physical solve refused"));
        assert_eq!(run.advance(1, || panic!("refused"), |_| -> Result<f64, &str> { panic!("refused") }), refused);
    }
}

#[test]
fn whole_study_admission_reserves_both_phases_and_checks_control_work_before_physics() {
    let mut p = plan();
    for total in [0, 128, 129, 1_000_001] {
        p.budget_max_samples = total;
        assert!(PceMeanExecution::new(&p, config()).is_err());
    }
    p = plan();
    p.method = PropagationMethod::PolynomialChaos;
    assert!(PceMeanExecution::new(&p, config()).is_err());
    p = plan();
    for weight in [f64::NAN, f64::INFINITY] {
        assert!(PceMeanExecution::new(&p, PceMeanConfig { control_weight: weight, ..config() }).is_err());
    }
    assert!(PceMeanExecution::new(&p, PceMeanConfig { fit_evaluations: 3, ..config() }).is_err());
    assert!(PceMeanExecution::new(&p, PceMeanConfig { regression: PceConfig::new(8, 96), ..config() }).is_err());
    p.budget_max_samples = 1_000_000;
    for i in 2..6 {
        p.parameters.push(ParameterUncertainty::gaussian(format!("x{i}"), 0.0, 1.0, "1"));
    }
    // 28 terms * 6 coordinates * (1,000,000 - 128) exceeds replay allowance.
    assert!(PceMeanExecution::new(&p, config()).is_err());
}
