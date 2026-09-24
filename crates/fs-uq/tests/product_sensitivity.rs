use fs_uq::{CorrelationModel, ParameterUncertainty, PropagationMethod, SobolExecution, UqPlan,
    UqExecution, UqStatus};

fn plan(rows: usize, dimensions: usize) -> UqPlan {
    let mut plan = UqPlan::new("response", PropagationMethod::MonteCarlo, rows * (dimensions + 2))
        .with_correlation(CorrelationModel::Independent);
    plan.seed = 81231;
    for i in 0..dimensions {
        plan = plan.with_parameter(ParameterUncertainty::uniform(&format!("x{i}"), -1.0, 1.0, "1"));
    }
    plan
}
fn run(plan: &UqPlan, model: impl Fn(&[f64]) -> f64) -> fs_uq::SobolReport {
    SobolExecution::new(plan).unwrap().advance(usize::MAX, || false, |x| Ok::<_, &str>(model(x)))
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance, "{actual} != {expected}");
}

#[test]
fn additive_and_interaction_models_recover_analytical_variance_shares() {
    let plan = plan(16_384, 3);
    let additive = run(&plan, |x| x[0] + 2.0*x[1]);
    assert_eq!(additive.status, UqStatus::Complete);
    let estimate = additive.estimate.unwrap();
    for (effect, expected) in estimate.effects.iter().zip([0.2, 0.8, 0.0]) {
        near(effect.first_order, expected, 0.04);
        near(effect.total_order, expected, 0.04);
    }
    assert_eq!(estimate.effects[2].total_order, 0.0);
    let interaction = run(&plan, |x| x[0]*x[1]);
    let effects = interaction.estimate.unwrap().effects;
    for effect in &effects[..2] {
        near(effect.first_order, 0.0, 0.04);
        near(effect.total_order, 1.0, 0.04);
    }
    near(effects[2].first_order, 0.0, 0.04);
    assert_eq!(effects[2].total_order, 0.0);
    assert!(effects.iter().map(|x| x.total_order).sum::<f64>() > 1.8);
}

#[test]
fn gaussian_variance_shares_use_the_existing_marginal_sampler() {
    let mut plan = plan(16_384, 2);
    plan.parameters = vec![
        ParameterUncertainty::gaussian("ambient", 290.0, 2.0, "K"),
        ParameterUncertainty::gaussian("load", 30.0, 3.0, "W"),
    ];
    // Var(X0) = 4; Var(2*X1) = 36, so shares are 0.1 and 0.9.
    let estimate = run(&plan, |x| x[0] + 2.0*x[1]).estimate.unwrap();
    for (effect, expected) in estimate.effects.iter().zip([0.1, 0.9]) {
        near(effect.first_order, expected, 0.04);
        near(effect.total_order, expected, 0.04);
    }
    assert_eq!(estimate.effects[0].parameter, "ambient");
    assert_eq!(estimate.effects[1].unit, "W");
}

#[test]
fn rows_reuse_counter_addressed_bases_and_replace_only_the_named_coordinate() {
    let plan = plan(3, 3);
    let mut base_plan = plan.clone(); base_plan.budget_max_samples = 6;
    let mut bases = Vec::new();
    UqExecution::new(&base_plan).unwrap().advance(6, || false, |x| {
        bases.push(x.to_vec()); Ok::<_, &str>(0.0)
    });
    let mut inputs = Vec::new();
    let report = SobolExecution::new(&plan).unwrap().advance(usize::MAX, || false, |x| {
        inputs.push(x.to_vec()); Ok::<_, &str>(x[0] + x[1])
    });
    assert_eq!(report.evaluations_attempted, 15);
    assert_eq!(report.evaluations_accepted, 15);
    for row in 0..3 {
        let group = &inputs[row*5..row*5+5];
        assert_eq!(group[0], bases[2*row]); assert_eq!(group[1], bases[2*row+1]);
        for coordinate in 0..3 {
            let mut expected = bases[2*row].clone();
            expected[coordinate] = bases[2*row+1][coordinate];
            assert_eq!(group[coordinate+2], expected);
        }
    }
}

#[test]
fn every_split_and_interrupted_hybrid_resume_without_repeating_paid_calls() {
    let plan = plan(4, 2);
    let full = run(&plan, |x| x[0] + x[1]*x[1]);
    for split in 0..=plan.budget_max_samples {
        let mut execution = SobolExecution::new(&plan).unwrap();
        let prefix = execution.advance(split, || false, |x| Ok::<_, &str>(x[0] + x[1]*x[1]));
        if split < plan.budget_max_samples { assert!(prefix.estimate.is_none()); }
        let mut execution = execution.clone();
        let mut calls = 0;
        let result = execution.advance(usize::MAX, || false, |x| {
            calls += 1; Ok::<_, &str>(x[0] + x[1]*x[1])
        });
        assert_eq!(calls, plan.budget_max_samples - split);
        assert_eq!(result, full);
    }
    let mut execution = SobolExecution::new(&plan).unwrap();
    execution.advance(3, || false, |x| Ok::<_, &str>(x[0] + x[1]*x[1]));
    let paid = execution.observations().to_vec();
    let mut interrupted = Vec::new();
    let report = execution.advance_interruptible(9, || false, |x| {
        interrupted = x.to_vec(); Ok::<_, &str>(None)
    });
    assert_eq!(report.status, UqStatus::Cancelled);
    assert_eq!(report.evaluations_attempted, 3);
    assert_eq!(execution.observations(), paid);
    let mut first = true;
    let result = execution.advance(usize::MAX, || false, |x| {
        if first { assert_eq!(x, interrupted); first = false; }
        Ok::<_, &str>(x[0] + x[1]*x[1])
    });
    assert_eq!(result, full);
    let terminal = execution.advance(1, || panic!("terminal poll"), |_| -> Result<f64,&str> {
        panic!("terminal callback")
    });
    assert_eq!(terminal, full);
}

#[test]
fn pre_callback_cancellation_and_model_failures_do_not_filter_designs() {
    let plan = plan(4, 2);
    let mut execution = SobolExecution::new(&plan).unwrap();
    let mut polls = 0;
    let partial = execution.advance(16, || { polls += 1; polls == 4 }, |_| Ok::<_, &str>(1.0));
    assert_eq!(partial.status, UqStatus::Cancelled);
    assert_eq!(partial.evaluations_accepted, 3);
    assert_eq!(partial.completed_rows, 0);
    assert!(partial.estimate.is_none());
    for result in [Err("solver refused"), Ok(f64::NAN), Ok(f64::INFINITY)] {
        let mut execution = execution.clone();
        let refusal = execution.advance(1, || false, |_| result);
        assert_eq!(refusal.status, UqStatus::Refused);
        assert_eq!(refusal.evaluations_attempted, 4);
        assert_eq!(refusal.evaluations_accepted, 3);
        assert!(refusal.estimate.is_none());
        assert!(refusal.rejection_reason.is_some());
        assert_eq!(execution.advance(16, || panic!("refused"), |_| -> Result<f64,&str> {
            panic!("refused")
        }), refusal);
    }
}

#[test]
fn dependence_layout_fixed_input_and_unsupported_method_refuse_at_admission() {
    let valid = plan(4, 2);
    let mut bad = valid.clone(); bad.correlation = CorrelationModel::Unknown;
    assert!(SobolExecution::new(&bad).is_err());
    bad.correlation = CorrelationModel::JointGaussian { matrix: vec![vec![1.0,0.0],vec![0.0,1.0]] };
    bad.parameters = vec![ParameterUncertainty::gaussian("a",0.0,1.0,"1"),
        ParameterUncertainty::gaussian("b",0.0,1.0,"1")];
    assert!(SobolExecution::new(&bad).is_err(), "zero correlation is not this adapter's independence declaration");
    for count in [0, 1, 4, 7, 15, 17, usize::MAX] {
        let mut bad = valid.clone(); bad.budget_max_samples = count;
        assert!(SobolExecution::new(&bad).is_err());
    }
    let mut bad = valid.clone(); bad.method = PropagationMethod::QuasiMonteCarlo;
    assert!(SobolExecution::new(&bad).is_err());
    let mut bad = valid.clone();
    bad.parameters[0] = ParameterUncertainty::uniform("fixed", 1.0, 1.0, "1");
    assert!(SobolExecution::new(&bad).is_err());
}

#[test]
fn ratios_survive_output_rescaling_extremes_and_common_offsets() {
    let plan = plan(128, 2);
    let baseline = run(&plan, |x| x[0] + 2.0*x[1]).estimate.unwrap();
    for scale in [1e200, -1e200, 1e-200] {
        let scaled = run(&plan, |x| scale*(x[0] + 2.0*x[1])).estimate.unwrap();
        for (a,b) in baseline.effects.iter().zip(scaled.effects) {
            near(a.first_order, b.first_order, 1e-13);
            near(a.total_order, b.total_order, 1e-13);
        }
    }
    // Exact integer increments remain distinguishable near 2^40.
    let small = run(&plan, |x| (x[0]*32.0).round() + 2.0*(x[1]*32.0).round()).estimate.unwrap();
    let offset = run(&plan, |x| 2_f64.powi(40) + (x[0]*32.0).round() + 2.0*(x[1]*32.0).round()).estimate.unwrap();
    assert_eq!(small.effects, offset.effects);
    let extreme = run(&plan, |x| if x[0] < 0.0 { -f64::MAX } else { f64::MAX });
    let estimate = extreme.estimate.unwrap();
    assert_eq!(estimate.effects[0].first_order, 1.0);
    assert!(estimate.effects.iter().all(|e| e.first_order.is_finite() && e.total_order.is_finite()));
    let constant = run(&plan, |_| f64::MAX);
    assert_eq!(constant.status, UqStatus::Complete);
    assert!(constant.estimate.is_none());
    assert!(constant.unavailable_reason.unwrap().contains("variance"));
}

#[test]
fn finite_sample_indices_are_not_clipped_into_plausible_ranges() {
    // Manufacture a valid deterministic response on the eight visited inputs;
    // these values make the second total index far greater than one.
    let plan = plan(2, 2);
    let mut ordinal = 0;
    let report = SobolExecution::new(&plan).unwrap().advance(8, || false, |_| {
        let values = [0.0,1.0,1.0,10.0, 1.0,0.0,0.0,-10.0];
        let value = values[ordinal]; ordinal += 1; Ok::<_, &str>(value)
    });
    let effects = report.estimate.unwrap().effects;
    assert!(effects[1].total_order > 1.0);
    assert!(effects[1].first_order < 0.0);
}
