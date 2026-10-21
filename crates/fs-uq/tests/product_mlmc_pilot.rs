use fs_uq::{CorrelationModel, MlmcExecution, ParameterUncertainty, PilotMlmcConfig,
    PilotMlmcExecution, PilotMlmcStage, PropagationMethod, UqPlan, UqStatus};

fn plan() -> UqPlan {
    UqPlan::new("cooling", PropagationMethod::MultilevelMonteCarlo, 8192)
        .with_parameter(ParameterUncertainty::uniform("initial_temperature", 0.0, 1.0, "K"))
        .with_correlation(CorrelationModel::Independent)
}
fn config() -> PilotMlmcConfig {
    PilotMlmcConfig { pilot_samples_per_level: 16, physical_model_costs: vec![1.0, 2.0, 4.0],
        target_standard_error: 0.025, minimum_production_samples_per_level: 4,
        correction_std_dev_floor: 1e-4 }
}
fn cooling(level: usize, x: &[f64]) -> Result<f64, &'static str> {
    let steps = 1usize << level;
    Ok(x[0] * (1.0 / (1.0 + 1.0 / steps as f64)).powi(steps as i32))
}

#[test]
fn pilot_selects_an_affordable_production_design_without_reusing_its_observations() {
    let mut plan = plan();
    plan.compliance_threshold = Some(0.2);
    let mut run = PilotMlmcExecution::new(&plan, config()).unwrap();
    let mut inputs = Vec::new();
    let report = run.advance(usize::MAX, || false, |level, x| {
        inputs.push((level, x[0].to_bits())); cooling(level, x)
    });
    assert_eq!(report.status, UqStatus::Complete);
    assert_eq!(report.stage, PilotMlmcStage::Production);
    assert_eq!(run.plan(), &plan);
    assert_eq!(report.pilot.model_calls_evaluated, 80);
    assert_eq!(report.model_calls_evaluated, report.model_calls_planned.unwrap());
    assert!(report.model_calls_evaluated <= plan.budget_max_samples);
    assert_ne!(report.pilot_seed, report.production_seed);
    let allocation = report.allocation.as_ref().unwrap();
    assert!(allocation.predicted_standard_error <= config().target_standard_error);
    let production = run.production().unwrap();
    assert_eq!(production.plan().seed, report.production_seed);
    let mut replay = MlmcExecution::new(production.plan(), &allocation.samples_per_level).unwrap();
    let mut replay_inputs = Vec::new();
    let reference = replay.advance(usize::MAX, || false, |level, x| {
        replay_inputs.push((level, x[0].to_bits())); cooling(level, x)
    });
    assert_eq!(inputs[80..], replay_inputs);
    assert_eq!(report.estimate, reference.estimate);
    assert_eq!(report.compliance, reference.compliance);
    assert_eq!(report.production.as_ref().unwrap(), &reference);
    for level in 0..3 { assert_eq!(production.corrections(level), replay.corrections(level)); }
    assert_ne!(inputs[0].1, inputs[80].1);
    let estimate = report.estimate.unwrap();
    let expected = 0.5 * (1.0_f64 / 1.25).powi(4);
    assert!((estimate.mean - expected).abs() < 8.0 * estimate.standard_error.unwrap() + 0.01);
    assert_eq!(report.observed_sampling_target_met,
        estimate.standard_error.map(|se| se <= config().target_standard_error));
}

#[test]
fn half_pairs_chunks_interrupted_calls_and_clone_preserve_both_streams_bitwise() {
    let mut full = PilotMlmcExecution::new(&plan(), config()).unwrap();
    let mut expected_inputs = Vec::new();
    let expected = full.advance(usize::MAX, || false, |level, x| {
        expected_inputs.push((level, x[0].to_bits())); cooling(level, x)
    });
    let mut split = PilotMlmcExecution::new(&plan(), config()).unwrap();
    let mut paid = Vec::new();
    split.advance(17, || false, |level, x| { paid.push((level, x[0].to_bits())); cooling(level, x) });
    assert!(split.report().pilot.pending_pair);
    let before = split.report();
    assert_eq!(split.advance(0, || panic!("zero allowance"), |_, _| -> Result<f64, &str> { panic!("zero allowance") }), before);
    let stopped = split.advance(1, || true, |_, _| -> Result<f64, &str> { panic!("cancelled") });
    assert_eq!(stopped.status, UqStatus::Cancelled);
    let paused = split.advance_interruptible(1, || false, |level, x| {
        assert_eq!((level, x[0].to_bits()), expected_inputs[17]); Ok::<_, &str>(None)
    });
    assert_eq!(paused.model_calls_evaluated, 17);
    let mut resumed = split.clone();
    for allowance in [1, 62, 1, 3, usize::MAX] {
        resumed.advance(allowance, || false, |level, x| { paid.push((level, x[0].to_bits())); cooling(level, x) });
    }
    assert_eq!(paid, expected_inputs);
    assert_eq!(resumed.report(), expected);
    assert_eq!(split.evaluations_attempted(), 17);
    assert_eq!(resumed.advance(10, || panic!("terminal"), |_, _| -> Result<f64, &str> { panic!("terminal") }), expected);
}

#[test]
fn cancellation_at_the_phase_boundary_never_repeats_paid_pilot_calls() {
    let mut run = PilotMlmcExecution::new(&plan(), config()).unwrap();
    let mut polls = 0;
    let report = run.advance(80, || { polls += 1; polls == 81 }, cooling);
    assert_eq!(report.pilot.status, UqStatus::Complete);
    assert_eq!(report.status, UqStatus::Cancelled);
    assert_eq!(report.stage, PilotMlmcStage::Pilot);
    assert!(run.production().is_none());
    assert_eq!(report.model_calls_evaluated, 80);
    let report = run.advance(1, || false, cooling);
    assert_eq!(report.stage, PilotMlmcStage::Production);
    assert_eq!(report.model_calls_evaluated, 81);
    assert_eq!(report.pilot.model_calls_evaluated, 80);
}

#[test]
fn impossible_target_retains_pilot_diagnostics_and_performs_no_production_calls() {
    let mut plan = plan();
    plan.budget_max_samples = 12;
    let mut cfg = config();
    cfg.physical_model_costs = vec![1.0, 2.0];
    cfg.pilot_samples_per_level = 2;
    cfg.minimum_production_samples_per_level = 2;
    cfg.target_standard_error = 1e-12;
    let mut run = PilotMlmcExecution::new(&plan, cfg).unwrap();
    let mut calls = 0;
    let report = run.advance(usize::MAX, || false, |_, _| { calls += 1; Ok::<_, &str>(0.0) });
    assert_eq!(calls, 6);
    assert_eq!(report.status, UqStatus::Refused);
    assert_eq!(report.pilot.status, UqStatus::Complete);
    assert!(report.allocation_refusal.is_some() && report.rejection_reason.is_some());
    assert!(report.estimate.is_none() && run.production().is_none());
    assert_eq!(run.advance(1, || panic!("terminal"), |_, _| -> Result<f64, &str> { panic!("terminal") }), report);
}

#[test]
fn model_failure_in_either_phase_is_terminal_without_censoring() {
    for production in [false, true] {
        for value in [Err("physical solver failed"), Ok(f64::NAN), Ok(f64::INFINITY)] {
            let mut run = PilotMlmcExecution::new(&plan(), config()).unwrap();
            if production { run.advance(80, || false, cooling); }
            let before = run.evaluations_attempted();
            let report = run.advance(1, || false, |_, _| value);
            assert_eq!(report.status, UqStatus::Refused);
            assert_eq!(report.model_calls_evaluated, before + 1);
            assert!(report.estimate.is_none() && report.rejection_reason.is_some());
            assert_eq!(run.advance(1, || panic!("terminal"), |_, _| -> Result<f64, &str> { panic!("terminal") }), report);
        }
    }
}

#[test]
fn constructor_rejects_incomplete_designs_and_preserves_joint_gaussian_admission() {
    let mut invalid = plan();
    invalid.budget_max_samples = 80; // Pilot consumes everything; no minimum production.
    assert!(PilotMlmcExecution::new(&invalid, config()).is_err());
    invalid = plan(); invalid.method = PropagationMethod::MonteCarlo;
    assert!(PilotMlmcExecution::new(&invalid, config()).is_err());
    for cfg in [PilotMlmcConfig { pilot_samples_per_level: 1, ..config() },
        PilotMlmcConfig { physical_model_costs: vec![], ..config() },
        PilotMlmcConfig { physical_model_costs: vec![f64::NAN], ..config() },
        PilotMlmcConfig { target_standard_error: 0.0, ..config() }] {
        assert!(PilotMlmcExecution::new(&plan(), cfg).is_err());
    }
    let joint = UqPlan::new("singular", PropagationMethod::MultilevelMonteCarlo, 8192)
        .with_parameter(ParameterUncertainty::gaussian("a", 0.0, 1.0, "1"))
        .with_parameter(ParameterUncertainty::gaussian("b", 0.0, 2.0, "1"))
        .with_correlation(CorrelationModel::JointGaussian { matrix: vec![vec![1.0, 1.0], vec![1.0, 1.0]] });
    let mut cfg = config(); cfg.target_standard_error = 0.1;
    let mut run = PilotMlmcExecution::new(&joint, cfg).unwrap();
    let report = run.advance(usize::MAX, || false, |level, x| {
        assert_eq!(x[1], 2.0 * x[0]); Ok::<_, &str>(x[0] / (level + 1) as f64)
    });
    assert_eq!(report.status, UqStatus::Complete);
    assert!(report.model_calls_evaluated <= joint.budget_max_samples);
}

#[test]
fn completed_production_reports_a_missed_sampling_target_without_resampling() {
    let mut cfg = config();
    cfg.physical_model_costs = vec![1.0]; cfg.pilot_samples_per_level = 2;
    cfg.minimum_production_samples_per_level = 2; cfg.target_standard_error = 0.1;
    let mut run = PilotMlmcExecution::new(&plan(), cfg).unwrap();
    run.advance(2, || false, |_, _| Ok::<_, &str>(0.0));
    // Deliberately adversarial callback responses test reporting, not physical
    // model validity. A pilot target must never manufacture a success claim.
    let mut i = 0;
    let report = run.advance(usize::MAX, || false, |_, _| { i += 1; Ok::<_, &str>(if i == 1 { -1.0 } else { 1.0 }) });
    assert_eq!(i, 2);
    assert_eq!(report.status, UqStatus::Complete);
    assert_eq!(report.observed_sampling_target_met, Some(false));
    assert_eq!(report.model_calls_evaluated, 4);
}
