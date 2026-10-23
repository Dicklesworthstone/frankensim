//! Physical mixed-law conditional pairs, not Gaussian approximations to tolerances.
use fs_blake3::ContentHash;
use fs_uq::{CorrelationModel, ParameterUncertainty, PropagationMethod, ShapleyExecution,
    UqExecution, UqPlan, UqStatus};
use fs_uq::product_shapley::ShapleyBatchExecution;

fn mixed(rows: usize) -> UqPlan {
    let mut plan = UqPlan::new("response", PropagationMethod::MonteCarlo, 14 * rows)
        .with_parameter(ParameterUncertainty::gaussian("scatter", 0.0, 1.0, "1"))
        .with_parameter(ParameterUncertainty::uniform("width", -1.0, 1.0, "m"))
        .with_parameter(ParameterUncertainty::uniform("load", -1.0, 1.0, "N"))
        .with_correlation(CorrelationModel::Independent);
    plan.seed = 83;
    plan
}
fn model(x: &[f64]) -> Result<f64, &'static str> { Ok(2.0 * x[0] + 3.0 * x[1] * x[2]) }
fn bits(x: &[f64]) -> Vec<u64> { x.iter().map(|v| v.to_bits()).collect() }

#[test]
fn every_mixed_coalition_freezes_the_exact_declared_physical_coordinates() {
    let plan = mixed(4);
    let mut run = ShapleyExecution::new(&plan).unwrap();
    let mut inputs = Vec::new();
    run.advance(56, || false, |x| { inputs.push(x.to_vec()); model(x) });
    let mut base = plan.clone(); base.budget_max_samples = 8;
    let mut reference = Vec::new();
    UqExecution::new(&base).unwrap().advance(8, || false, |x| {
        reference.push(x.to_vec()); model(x)
    });
    for (row, values) in inputs.chunks_exact(14).enumerate() {
        let a = &reference[2 * row]; let b = &reference[2 * row + 1];
        for subset in 0..7 {
            assert_eq!(bits(&values[2 * subset]), bits(a));
            for j in 0..3 {
                let expected = if subset & (1 << j) == 0 { b[j] } else { a[j] };
                assert_eq!(values[2 * subset + 1][j].to_bits(), expected.to_bits());
            }
        }
        assert!(values.iter().all(|x| (-1.0..=1.0).contains(&x[1])
            && (-1.0..=1.0).contains(&x[2])));
    }
    let effects = run.report().estimate.unwrap().effects;
    for (effect, parameter) in effects.iter().zip(&plan.parameters) {
        assert_eq!(effect.parameter, parameter.name);
        assert_eq!(effect.unit, parameter.unit);
    }
}

#[test]
fn mixed_nonlinear_interaction_is_shared_without_changing_uniform_variance() {
    let mut run = ShapleyExecution::new(&mixed(8192)).unwrap();
    let estimate = run.advance(usize::MAX, || false, model).estimate.unwrap();
    // Gaussian term has variance 4. Independent U[-1,1] products have variance
    // 1/9, so 3*U*V has variance 1, shared equally between the two uniforms.
    for (effect, expected) in estimate.effects.iter().zip([0.8, 0.1, 0.1]) {
        assert!((effect.index - expected).abs() < 0.04, "{effect:?}");
    }
    assert!((estimate.base_output_std_dev.unwrap() - fs_math::det::sqrt(5.0)).abs() < 0.1);
    assert!((estimate.effects.iter().map(|e| e.index).sum::<f64>() - 1.0).abs() < 1e-14);
}

#[test]
fn all_uniform_three_way_interaction_has_the_expected_physical_effects() {
    let mut plan = mixed(8192);
    plan.parameters[0] = ParameterUncertainty::uniform("scatter", -1.0, 1.0, "1");
    let mut run = ShapleyExecution::new(&plan).unwrap();
    let estimate = run.advance(usize::MAX, || false,
        |x| Ok::<_, &str>(x[0] * x[1] * x[2])).estimate.unwrap();
    for effect in &estimate.effects { assert!((effect.index - 1.0 / 3.0).abs() < 0.04); }
    assert!((estimate.base_output_std_dev.unwrap() - fs_math::det::sqrt(1.0 / 27.0)).abs() < 0.01);
}

#[test]
fn dummy_coordinates_stay_dummy_in_every_marginal_order() {
    for active in 0..3 {
        let mut plan = mixed(8);
        plan.parameters.swap(0, active);
        let mut run = ShapleyExecution::new(&plan).unwrap();
        let estimate = run.advance(usize::MAX, || false,
            |x| Ok::<_, &str>(x[active])).estimate.unwrap();
        for (j, effect) in estimate.effects.iter().enumerate() {
            let expected = if j == active { 1.0 } else { 0.0 };
            assert!((effect.index - expected).abs() < 1e-14);
        }
    }
}

#[test]
fn every_mixed_prefix_and_interrupted_half_pair_restore_identically() {
    let plan = mixed(2); let identity = ContentHash([37; 32]);
    let mut full = ShapleyExecution::new(&plan).unwrap();
    let mut inputs = Vec::new();
    let expected = full.advance(28, || false, |x| { inputs.push(bits(x)); model(x) });
    for split in 0..=28 {
        let mut run = ShapleyExecution::new(&plan).unwrap();
        run.advance(split, || false, model);
        if split < 28 {
            run.advance_interruptible(1, || false, |x| {
                assert_eq!(bits(x), inputs[split]); Ok::<_, &str>(None)
            });
        }
        let bytes = run.checkpoint(identity).unwrap();
        let mut restored = ShapleyExecution::restore(&plan, identity, &bytes).unwrap();
        assert_eq!(restored.report(), run.report());
        let mut ordinal = split;
        let actual = restored.advance(usize::MAX, || false, |x| {
            assert_eq!(bits(x), inputs[ordinal]); ordinal += 1; model(x)
        });
        assert_eq!(ordinal, 28); assert_eq!(actual, expected);
        assert_eq!(restored.checkpoint(identity).unwrap(), full.checkpoint(identity).unwrap());
    }
}

#[test]
fn mixed_parallel_recovery_retains_later_successes_and_the_original_report() {
    let plan = mixed(3); let identity = ContentHash([38; 32]);
    let mut serial = ShapleyExecution::new(&plan).unwrap();
    let mut inputs = Vec::new();
    let expected = serial.advance(42, || false, |x| { inputs.push(bits(x)); model(x) });
    let mut batch = ShapleyBatchExecution::new(&plan).unwrap();
    let mut calls = 0;
    let partial = batch.advance_batch::<&str>(7, || false, |jobs, slots| {
        for j in [1, 4] { slots[j] = Some(model(&jobs[j].parameters)); calls += 1; }
    }).unwrap();
    assert_eq!(partial.retained_out_of_order, 2);
    assert_eq!(partial.contiguous_observations, 0);
    assert!(partial.result.is_none());
    let bytes = batch.checkpoint(identity).unwrap();
    let mut restored = ShapleyBatchExecution::restore(&plan, identity, &bytes).unwrap();
    while restored.report().status != UqStatus::Complete {
        restored.advance_batch::<&str>(3, || false, |jobs, slots| {
            for (job, slot) in jobs.iter().zip(slots) {
                assert_ne!(job.ordinal, 1); assert_ne!(job.ordinal, 4);
                assert_eq!(bits(&job.parameters), inputs[job.ordinal]);
                *slot = Some(model(&job.parameters)); calls += 1;
            }
        }).unwrap();
    }
    assert_eq!(calls, 42);
    assert_eq!(restored.report().result.unwrap(), expected);
    assert_eq!(bits(restored.observations()), bits(serial.observations()));
}

#[test]
fn uniform_supports_with_unrepresentable_width_use_original_convex_interpolation() {
    let plan = UqPlan::new("normalized", PropagationMethod::MonteCarlo, 32)
        .with_parameter(ParameterUncertainty::uniform("x", -f64::MAX, f64::MAX, "m"))
        .with_correlation(CorrelationModel::Independent);
    let mut run = ShapleyExecution::new(&plan).unwrap();
    let report = run.advance(32, || false, |x| {
        assert!(x[0].is_finite()); Ok::<_, &str>(x[0] / f64::MAX)
    });
    assert_eq!(report.status, UqStatus::Complete);
    assert_eq!(report.estimate.unwrap().effects[0].index, 1.0);
}

#[test]
fn mixed_laws_and_saved_physical_units_cannot_be_silently_reinterpreted() {
    let plan = mixed(2); let identity = ContentHash([39; 32]);
    let mut run = ShapleyExecution::new(&plan).unwrap();
    run.advance(3, || false, model);
    let saved = run.checkpoint(identity).unwrap();
    for correlation in [CorrelationModel::Unknown,
        CorrelationModel::JointGaussian { matrix: vec![vec![1.0,0.0,0.0],vec![0.0,1.0,0.0],vec![0.0,0.0,1.0]] }]
    {
        let mut bad = plan.clone(); bad.correlation = correlation;
        assert!(ShapleyExecution::new(&bad).is_err());
    }
    for marginal in [ParameterUncertainty::uniform("width", 1.0, 1.0, "m"),
        ParameterUncertainty::interval("width", -1.0, 1.0, "m")]
    {
        let mut bad = plan.clone(); bad.parameters[1] = marginal;
        assert!(ShapleyExecution::new(&bad).is_err());
    }
    for marginal in [ParameterUncertainty::uniform("width", -2.0, 1.0, "m"),
        ParameterUncertainty::uniform("width", -1.0, 1.0, "mm"),
        ParameterUncertainty::gaussian("width", 0.0, 1.0, "m")]
    {
        let mut changed = plan.clone(); changed.parameters[1] = marginal;
        assert!(ShapleyExecution::restore(&changed, identity, &saved).is_err());
    }
    let failed = run.advance(1, || false, |_| Err::<f64, _>("physical refusal"));
    assert_eq!(failed.status, UqStatus::Refused);
    assert!(failed.estimate.is_none()); assert!(run.checkpoint(identity).is_err());
}
