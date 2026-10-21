use fs_uq::{CorrelationModel, MultiQoiExecution, ParameterUncertainty, PropagationMethod,
    QoiDefinition, UqExecution, UqPlan, UqStatus};

fn plan(n: usize) -> UqPlan {
    UqPlan::new("thermal study", PropagationMethod::MonteCarlo, n)
        .with_parameter(ParameterUncertainty::uniform("power", 4.0, 6.0, "W"))
        .with_parameter(ParameterUncertainty::gaussian("resistance", 2.0, 0.1, "K/W"))
        .with_correlation(CorrelationModel::Independent)
}
fn quantities() -> Vec<QoiDefinition> {
    vec![QoiDefinition::new("temperature", "K").with_upper_limit(310.0),
         QoiDefinition::new("power", "W").with_upper_limit(5.4),
         QoiDefinition::new("temperature rise", "K")]
}
fn physical(x: &[f64]) -> Result<Vec<f64>, &'static str> {
    Ok(vec![300.0 + x[0] * x[1], x[0], x[0] * x[1]])
}

#[test]
fn one_physical_solve_supplies_every_canonical_marginal() {
    let plan = plan(64);
    let definitions = quantities();
    let mut run = MultiQoiExecution::new(&plan, &definitions).unwrap();
    assert!(run.report().outputs.is_empty());
    let mut actual_inputs = Vec::new();
    let report = run.advance(64, || false, |x| {
        actual_inputs.push(x.to_vec());
        physical(x)
    });
    assert_eq!(report.status, UqStatus::Complete);
    assert_eq!(actual_inputs.len(), 64);
    assert_eq!(report.samples_evaluated, 64);
    assert_eq!(report.samples_accepted, 64);
    assert_eq!(run.plan(), &plan);
    assert_eq!(run.quantities(), definitions);
    for (j, q) in definitions.iter().enumerate() {
        let mut scalar_plan = plan.clone();
        scalar_plan.target_qoi = q.name.clone();
        scalar_plan.compliance_threshold = q.upper_limit;
        let mut scalar = UqExecution::new(&scalar_plan).unwrap();
        let mut inputs = Vec::new();
        let expected = scalar.advance(64, || false, |x| {
            inputs.push(x.to_vec());
            Ok::<_, &str>(physical(x)?[j])
        });
        assert_eq!(inputs, actual_inputs);
        assert_eq!(report.outputs[j].quantity, *q);
        assert_eq!(report.outputs[j].statistics, expected);
        assert_eq!(run.observations(j).unwrap(), scalar.observations());
    }
    for (i, x) in actual_inputs.iter().enumerate() {
        assert_eq!(run.observation(i).unwrap(), physical(x).unwrap());
    }
    assert!(run.observations(3).is_none());
    assert!(run.observation(64).is_none());
    let successes = actual_inputs.iter().filter(|x| {
        let y = physical(x).unwrap(); y[0] <= 310.0 && y[1] <= 5.4
    }).count();
    assert_eq!(report.joint_compliant_samples, Some(successes));
    assert_eq!(report.joint_compliance.unwrap().mean, successes as f64 / 64.0);
}

#[test]
fn mutually_exclusive_limits_do_not_multiply_marginal_probabilities() {
    let plan = UqPlan::new("opposed loads", PropagationMethod::MonteCarlo, 256)
        .with_parameter(ParameterUncertainty::uniform("x", -1.0, 1.0, "N"));
    let definitions = [QoiDefinition::new("left", "N").with_upper_limit(-0.1),
        QoiDefinition::new("right", "N").with_upper_limit(-0.1)];
    let mut run = MultiQoiExecution::new(&plan, &definitions).unwrap();
    let report = run.advance(256, || false, |x| Ok::<_, &str>(vec![x[0], -x[0]]));
    assert_eq!(report.status, UqStatus::Complete);
    for output in &report.outputs {
        let p = output.statistics.probability_of_compliance.unwrap();
        assert!(p > 0.2 && p < 0.7);
    }
    assert_eq!(report.joint_compliant_samples, Some(0));
    assert_eq!(report.joint_compliance.unwrap().mean, 0.0);
    assert!(matches!(report.evidence_color,
        fs_evidence::Color::Estimated { dispersion, .. } if dispersion == f64::INFINITY));
}

#[test]
fn chunks_interruption_and_clone_keep_every_paid_row_aligned() {
    let plan = plan(16);
    let mut full = MultiQoiExecution::new(&plan, &quantities()).unwrap();
    let mut expected_inputs = Vec::new();
    let expected = full.advance(99, || false, |x| {
        expected_inputs.push(x.to_vec()); physical(x)
    });
    let mut run = MultiQoiExecution::new(&plan, &quantities()).unwrap();
    let first = run.advance(1, || false, physical);
    assert!(first.outputs.iter().all(|r| r.statistics.std_dev.is_none()));
    assert!(first.joint_compliance.unwrap().standard_error.is_none());
    assert_eq!(run.advance(0, || panic!("zero"), |_| -> Result<Vec<f64>, &str> {
        panic!("zero")
    }), run.report());
    let cancelled = run.advance(1, || true, |_| -> Result<Vec<f64>, &str> { panic!("cancelled") });
    assert_eq!(cancelled.status, UqStatus::Cancelled);
    assert_eq!(cancelled.samples_accepted, 1);
    let interrupted = run.advance_interruptible(99, || false, |x| {
        assert_eq!(x, expected_inputs[1]); Ok::<_, &str>(None)
    });
    assert_eq!(interrupted.samples_evaluated, 1);
    let mut clone = run.clone();
    let mut inputs = vec![expected_inputs[0].clone()];
    for allowance in [1, 4, 99] {
        clone.advance(allowance, || false, |x| { inputs.push(x.to_vec()); physical(x) });
    }
    assert_eq!(inputs, expected_inputs);
    assert_eq!(clone.report(), expected);
    for j in 0..3 { assert_eq!(clone.observations(j), full.observations(j)); }
    assert_eq!(run.samples_accepted(), 1);
    assert_eq!(clone.advance(1, || panic!("complete"), |_| -> Result<Vec<f64>, &str> {
        panic!("complete")
    }), expected);
}

#[test]
fn malformed_vectors_refuse_the_whole_draw_without_partial_column_updates() {
    for bad in [vec![], vec![1.0], vec![1.0, 2.0], vec![1.0, f64::NAN, 3.0],
        vec![1.0, 2.0, f64::INFINITY], vec![1.0, 2.0, 3.0, 4.0]] {
        let mut run = MultiQoiExecution::new(&plan(8), &quantities()).unwrap();
        run.advance(3, || false, physical);
        let before: Vec<_> = (0..3).map(|i| run.observation(i).unwrap()).collect();
        let refused = run.advance(1, || false, |_| Ok::<_, &str>(bad.clone()));
        assert_eq!(refused.status, UqStatus::Refused);
        assert_eq!(refused.samples_evaluated, 4);
        assert_eq!(refused.samples_accepted, 3);
        assert!(refused.outputs.is_empty() && refused.joint_compliance.is_none());
        for j in 0..3 { assert_eq!(run.observations(j).unwrap().len(), 3); }
        assert_eq!((0..3).map(|i| run.observation(i).unwrap()).collect::<Vec<_>>(), before);
        assert_eq!(run.advance(1, || panic!("terminal"), |_| -> Result<Vec<f64>, &str> {
            panic!("terminal")
        }), refused);
    }
    let mut run = MultiQoiExecution::new(&plan(8), &quantities()).unwrap();
    let refused = run.advance(1, || false, |_| Err::<Vec<f64>, _>("physical solver failed"));
    assert_eq!(refused.samples_evaluated, 1);
    assert_eq!(refused.samples_accepted, 0);
    assert!(refused.rejection_reason.unwrap().contains("physical solver failed"));
}

#[test]
fn numerical_refusal_of_any_completed_marginal_hides_all_estimates() {
    let quantities = [QoiDefinition::new("small", "K"), QoiDefinition::new("large", "Pa")];
    for chunk in [1, 2] {
        let mut run = MultiQoiExecution::new(&plan(2), &quantities).unwrap();
        let mut calls = 0;
        for _ in 0..2/chunk {
            run.advance(chunk, || false, |_| {
                calls += 1;
                Ok::<_, &str>(vec![0.0, if calls == 1 { f64::MAX } else { -f64::MAX }])
            });
        }
        let report = run.report();
        assert_eq!(report.status, UqStatus::Refused);
        assert_eq!(report.samples_evaluated, 2);
        assert_eq!(report.samples_accepted, 2);
        assert!(report.outputs.is_empty() && report.joint_compliance.is_none());
        assert!(report.rejection_reason.unwrap().contains("large"));
    }
}

#[test]
fn output_admission_is_explicit_and_unconstrained_outputs_do_not_create_events() {
    let p = plan(8);
    assert!(MultiQoiExecution::new(&p, &[]).is_err());
    for definitions in [vec![QoiDefinition::new("", "K")], vec![QoiDefinition::new("q", "")],
        vec![QoiDefinition::new("q", "K").with_upper_limit(f64::NAN)],
        vec![QoiDefinition::new("q", "K"), QoiDefinition::new("q", "Pa")],
        (0..33).map(|i| QoiDefinition::new(format!("q{i}"), "1")).collect()] {
        assert!(MultiQoiExecution::new(&p, &definitions).is_err());
    }
    assert!(MultiQoiExecution::new(&p.clone().with_compliance_threshold(1.0), &quantities()).is_err());
    let mut huge = p.clone(); huge.budget_max_samples = 1_000_000;
    let many: Vec<_> = (0..5).map(|i| QoiDefinition::new(format!("q{i}"), "1")).collect();
    assert!(MultiQoiExecution::new(&huge, &many).is_err());
    let mut unknown = p.clone(); unknown.correlation = CorrelationModel::Unknown;
    assert!(MultiQoiExecution::new(&unknown, &quantities()).is_err());
    let mut run = MultiQoiExecution::new(&p, &[QoiDefinition::new("value", "W")]).unwrap();
    let report = run.advance(8, || false, |x| Ok::<_, &str>(vec![x[0]]));
    assert_eq!(report.status, UqStatus::Complete);
    assert!(report.joint_compliance.is_none() && report.joint_compliant_samples.is_none());
}

#[test]
fn joint_gaussian_singular_inputs_keep_original_physical_sampling() {
    let p = UqPlan::new("linked loads", PropagationMethod::MonteCarlo, 16)
        .with_parameter(ParameterUncertainty::gaussian("a", 2.0, 3.0, "N"))
        .with_parameter(ParameterUncertainty::gaussian("b", -1.0, 2.0, "N"))
        .with_correlation(CorrelationModel::JointGaussian { matrix: vec![vec![1.0, 1.0], vec![1.0, 1.0]] });
    let mut run = MultiQoiExecution::new(&p, &[QoiDefinition::new("sum", "N"),
        QoiDefinition::new("difference", "N")]).unwrap();
    let mut inputs = Vec::new();
    run.advance(16, || false, |x| { inputs.push(x.to_vec()); Ok::<_, &str>(vec![x[0]+x[1], x[0]-x[1]]) });
    let mut scalar = UqExecution::new(&p).unwrap();
    let mut expected = Vec::new();
    scalar.advance(16, || false, |x| { expected.push(x.to_vec()); Ok::<_, &str>(x[0]+x[1]) });
    assert_eq!(inputs, expected);
    assert_eq!(run.observations(0).unwrap(), scalar.observations());
}
