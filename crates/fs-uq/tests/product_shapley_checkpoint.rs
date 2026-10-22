use fs_blake3::{ContentHash, hash_domain};
use fs_uq::{CorrelationModel, ParameterUncertainty, PropagationMethod,
    ShapleyExecution, UqCheckpointError, UqExecution, UqPlan, UqStatus};

fn plan(correlation: CorrelationModel) -> UqPlan {
    UqPlan::new("thermal response", PropagationMethod::MonteCarlo, 42)
        .with_parameter(ParameterUncertainty::gaussian("ambient", 300.0, 1.0, "K"))
        .with_parameter(ParameterUncertainty::gaussian("power", 5.0, 0.2, "W"))
        .with_parameter(ParameterUncertainty::gaussian("resistance", 2.0, 0.1, "K/W"))
        .with_correlation(correlation)
}
fn model() -> ContentHash { hash_domain("test-shapley-model", b"T=ambient+power*resistance") }
fn physics(x: &[f64]) -> Result<f64, &'static str> { Ok(x[0] + x[1] * x[2]) }
fn bits(values: &[f64]) -> Vec<u64> { values.iter().map(|x| x.to_bits()).collect() }
fn reseal(bytes: &mut [u8]) {
    let end = bytes.len() - 32;
    let digest = hash_domain("org.frankensim.uq.shapley.checkpoint.v1", &bytes[..end]);
    bytes[end..].copy_from_slice(digest.as_bytes());
}

#[test]
fn every_half_pair_and_row_boundary_restores_without_repaying_physics() {
    for law in [CorrelationModel::Independent, CorrelationModel::JointGaussian {
        matrix: vec![vec![1.0, 0.5, 0.5], vec![0.5, 1.0, 0.5], vec![0.5, 0.5, 1.0]],
    }, CorrelationModel::JointGaussian { matrix: vec![vec![1.0; 3]; 3] }] {
        let plan = plan(law);
        let mut full = ShapleyExecution::new(&plan).unwrap();
        let mut inputs = Vec::new();
        let expected = full.advance(usize::MAX, || false, |x| {
            inputs.push(bits(x)); physics(x)
        });
        let final_bytes = full.checkpoint(model()).unwrap();
        for prefix in 0..=plan.budget_max_samples {
            let mut run = ShapleyExecution::new(&plan).unwrap();
            run.advance(prefix, || false, physics);
            if prefix < plan.budget_max_samples {
                run.advance(1, || true, |_| -> Result<f64, &str> { panic!("cancelled") });
            }
            let saved = run.checkpoint(model()).unwrap();
            let mut restored = ShapleyExecution::restore(&plan, model(), &saved).unwrap();
            assert_eq!(restored.report(), run.report());
            assert_eq!(restored.checkpoint(model()).unwrap(), saved);
            let mut next = prefix;
            let actual = restored.advance(usize::MAX, || false, |x| {
                assert_eq!(bits(x), inputs[next]); next += 1; physics(x)
            });
            assert_eq!(next, plan.budget_max_samples);
            assert_eq!(actual, expected);
            assert_eq!(bits(restored.observations()), bits(full.observations()));
            assert_eq!(restored.checkpoint(model()).unwrap(), final_bytes);
        }
    }
}

#[test]
fn unfinished_conditional_partner_retries_its_exact_input_after_restore() {
    let plan = plan(CorrelationModel::JointGaussian { matrix: vec![vec![1.0; 3]; 3] });
    let mut run = ShapleyExecution::new(&plan).unwrap();
    run.advance(3, || false, physics);
    let mut interrupted = Vec::new();
    run.advance_interruptible(1, || false, |x| {
        interrupted = bits(x); Ok::<_, &str>(None)
    });
    let mut restored = ShapleyExecution::restore(&plan, model(), &run.checkpoint(model()).unwrap()).unwrap();
    assert_eq!(restored.observations().len(), 3);
    let report = restored.advance(1, || false, |x| {
        assert_eq!(bits(x), interrupted); physics(x)
    });
    assert_eq!(report.evaluations_attempted, 4);
    assert!(report.estimate.is_none());
}

#[test]
fn changed_semantics_and_other_sampling_codecs_cannot_resume_a_shapley_prefix() {
    let plan = plan(CorrelationModel::Independent);
    let run = ShapleyExecution::new(&plan).unwrap();
    let bytes = run.checkpoint(model()).unwrap();
    assert!(matches!(ShapleyExecution::restore(&plan, ContentHash([42; 32]), &bytes),
        Err(UqCheckpointError::IdentityMismatch)));
    for change in 0..6 {
        let mut other = plan.clone();
        match change {
            0 => other.seed ^= 1,
            1 => other.target_qoi.push_str("changed"),
            2 => other.parameters[0].unit = "degC".into(),
            3 => other.parameters.swap(0, 1),
            4 => other.budget_max_samples += 14,
            _ => other.compliance_threshold = Some(320.0),
        }
        assert!(ShapleyExecution::restore(&other, model(), &bytes).is_err());
    }
    let mc = UqExecution::new(&plan).unwrap().checkpoint(model()).unwrap();
    assert!(ShapleyExecution::restore(&plan, model(), &mc).is_err());
    assert!(UqExecution::restore(&plan, model(), &bytes).is_err());
}

#[test]
fn malformed_count_status_values_and_corruption_are_rejected() {
    let plan = plan(CorrelationModel::Independent);
    let mut run = ShapleyExecution::new(&plan).unwrap();
    run.advance(3, || false, physics);
    let saved = run.checkpoint(model()).unwrap();
    for length in 0..saved.len() {
        assert!(ShapleyExecution::restore(&plan, model(), &saved[..length]).is_err());
    }
    let mut trailing = saved.clone(); trailing.push(0);
    assert!(ShapleyExecution::restore(&plan, model(), &trailing).is_err());
    let mut corrupt = saved.clone(); corrupt[49] ^= 1;
    assert!(matches!(ShapleyExecution::restore(&plan, model(), &corrupt),
        Err(UqCheckpointError::IntegrityMismatch)));
    for status in [2, 3, 255] {
        let mut bad = saved.clone(); bad[40] = status; reseal(&mut bad);
        assert!(ShapleyExecution::restore(&plan, model(), &bad).is_err());
    }
    for count in [0, 4, u64::MAX] {
        let mut bad = saved.clone(); bad[41..49].copy_from_slice(&count.to_le_bytes()); reseal(&mut bad);
        assert!(ShapleyExecution::restore(&plan, model(), &bad).is_err());
    }
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let mut bad = saved.clone(); bad[49..57].copy_from_slice(&value.to_bits().to_le_bytes()); reseal(&mut bad);
        assert!(ShapleyExecution::restore(&plan, model(), &bad).is_err());
    }
}

#[test]
fn undefined_complete_effects_and_signed_zero_bits_remain_terminal() {
    let plan = plan(CorrelationModel::Independent);
    let mut run = ShapleyExecution::new(&plan).unwrap();
    let mut ordinal = 0;
    let expected = run.advance(usize::MAX, || false, |_| {
        ordinal += 1; Ok::<_, &str>(if ordinal % 2 == 0 { 0.0 } else { -0.0 })
    });
    assert_eq!(expected.status, UqStatus::Complete);
    assert!(expected.estimate.is_none() && expected.unavailable_reason.is_some());
    let bytes = run.checkpoint(model()).unwrap();
    let mut restored = ShapleyExecution::restore(&plan, model(), &bytes).unwrap();
    assert_eq!(bits(restored.observations()), bits(run.observations()));
    assert_eq!(restored.advance(1, || panic!("complete"),
        |_| -> Result<f64, &str> { panic!("complete") }), expected);
}

#[test]
fn failed_physics_cannot_be_serialized_as_a_resumable_selected_prefix() {
    for result in [Err("model failed"), Ok(f64::NAN), Ok(f64::INFINITY)] {
        let mut run = ShapleyExecution::new(&plan(CorrelationModel::Independent)).unwrap();
        run.advance(3, || false, physics);
        assert_eq!(run.advance(1, || false, |_| result).status, UqStatus::Refused);
        assert_eq!(run.checkpoint(model()), Err(UqCheckpointError::NotResumable));
    }
}
