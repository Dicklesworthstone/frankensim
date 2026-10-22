//! Recovery retains paid physics and its original, unnormalized likelihoods.
use fs_blake3::{ContentHash, hash_domain};
use fs_uq::{CorrelationModel, ImportanceConfig, ImportanceExecution, ParameterUncertainty,
    PropagationMethod, UqCheckpointError, UqPlan, UqStatus};

fn plan() -> UqPlan {
    let mut plan = UqPlan::new("force", PropagationMethod::MonteCarlo, 17)
        .with_parameter(ParameterUncertainty::gaussian("load", 1.0, 0.2, "N"))
        .with_parameter(ParameterUncertainty::gaussian("other", 3.0, 0.5, "N"))
        .with_correlation(CorrelationModel::JointGaussian { matrix: vec![vec![1.0, 1.0], vec![1.0, 1.0]] })
        .with_compliance_threshold(4.5);
    plan.seed = 83;
    plan
}
fn proposal() -> ImportanceConfig {
    ImportanceConfig { germ_shift: vec![2.0, -1.0], defensive_fraction: 0.2 }
}
fn model() -> ContentHash { hash_domain("importance-test-model", b"sum of two forces") }
fn response(x: &[f64]) -> Result<f64, &'static str> { Ok(x[0] + x[1]) }
fn reseal(bytes: &mut [u8]) {
    let end = bytes.len() - 32;
    let sum = hash_domain("org.frankensim.uq.importance.checkpoint.v1", &bytes[..end]);
    bytes[end..].copy_from_slice(sum.as_bytes());
}

#[test]
fn every_prefix_restores_identical_samples_weights_and_final_probability() {
    let plan = plan(); let proposal = proposal();
    let mut full = ImportanceExecution::new(&plan, proposal.clone()).unwrap();
    let mut inputs = Vec::new();
    full.advance(99, || false, |x| { inputs.push(x.to_vec()); response(x) });
    let expected = full.checkpoint(model()).unwrap();
    for prefix in 0..=plan.budget_max_samples {
        let mut run = ImportanceExecution::new(&plan, proposal.clone()).unwrap();
        run.advance(prefix, || false, response);
        let bytes = run.checkpoint(model()).unwrap();
        assert_eq!(&bytes[..8], b"FSIMP001");
        assert_eq!(bytes.len(), 81 + 16 * prefix);
        let mut restored = ImportanceExecution::restore(&plan, proposal.clone(), model(), &bytes).unwrap();
        assert_eq!(restored.report(), run.report());
        assert_eq!(restored.observations(), run.observations());
        assert_eq!(restored.checkpoint(model()).unwrap(), bytes);
        let mut calls = 0;
        restored.advance(99, || false, |x| {
            assert_eq!(x, inputs[prefix + calls]); calls += 1; response(x)
        });
        assert_eq!(calls, plan.budget_max_samples - prefix);
        assert_eq!(restored.report(), full.report());
        assert_eq!(restored.checkpoint(model()).unwrap(), expected);
        assert_eq!(restored.assess_failure_probability(0.05, 0.01), full.assess_failure_probability(0.05, 0.01));
    }
}

#[test]
fn interrupted_callback_retries_the_same_unpaid_physical_input() {
    let plan = plan(); let proposal = proposal();
    let mut run = ImportanceExecution::new(&plan, proposal.clone()).unwrap();
    run.advance(3, || false, response);
    let mut pending = Vec::new();
    run.advance_interruptible(5, || false, |x| {
        pending = x.to_vec(); Ok::<_, &str>(None)
    });
    let bytes = run.checkpoint(model()).unwrap();
    let mut restored = ImportanceExecution::restore(&plan, proposal, model(), &bytes).unwrap();
    assert_eq!(restored.report().status, UqStatus::Cancelled);
    assert_eq!(restored.report().evaluations_attempted, 3);
    let mut calls = 0;
    restored.advance(1, || false, |x| { assert_eq!(x, pending); calls += 1; response(x) });
    assert_eq!(calls, 1);
    assert_eq!(restored.report().evaluations_accepted, 4);
    assert_eq!(&restored.observations()[..3], run.observations());
}

#[test]
fn target_shortcuts_signed_zero_and_threshold_equality_retain_exact_bits() {
    let mut plan = plan(); plan.budget_max_samples = 4; plan.compliance_threshold = Some(0.0);
    let values = [-0.0, 0.0, f64::MIN_POSITIVE, -f64::MAX];
    for proposal in [ImportanceConfig { germ_shift: vec![0.0, -0.0], defensive_fraction: 0.125 },
        ImportanceConfig { germ_shift: vec![4.0, -2.0], defensive_fraction: 1.0 }]
    {
        let mut run = ImportanceExecution::new(&plan, proposal.clone()).unwrap(); let mut i = 0;
        run.advance(4, || false, |_| { let value = values[i]; i += 1; Ok::<_, &str>(value) });
        let bytes = run.checkpoint(model()).unwrap();
        let mut restored = ImportanceExecution::restore(&plan, proposal, model(), &bytes).unwrap();
        for (observation, value) in restored.observations().iter().zip(values) {
            assert_eq!(observation.qoi.to_bits(), value.to_bits());
            assert_eq!(observation.likelihood_weight, 1.0);
            assert_eq!(observation.failed, value > 0.0);
        }
        let report = restored.advance(99, || panic!("complete"), |_| -> Result<f64, &str> { panic!("complete") });
        assert_eq!(report.status, UqStatus::Complete);
        assert_eq!(restored.checkpoint(model()).unwrap(), bytes);
    }
}

#[test]
fn proposal_event_law_budget_and_model_identity_cannot_change_on_restore() {
    let original = plan(); let proposal = proposal();
    let mut run = ImportanceExecution::new(&original, proposal.clone()).unwrap();
    run.advance(3, || false, response); let bytes = run.checkpoint(model()).unwrap();
    for field in 0..5 {
        let mut changed = original.clone();
        match field {
            0 => changed.seed += 1,
            1 => changed.budget_max_samples += 1,
            2 => changed.compliance_threshold = Some(5.0),
            3 => changed.parameters[0].unit = "Pa".into(),
            _ => changed.correlation = CorrelationModel::Independent,
        }
        assert!(matches!(ImportanceExecution::restore(&changed, proposal.clone(), model(), &bytes), Err(UqCheckpointError::IdentityMismatch)));
    }
    for changed in [ImportanceConfig { defensive_fraction: 0.3, ..proposal.clone() },
        ImportanceConfig { germ_shift: vec![2.0, 1.0], ..proposal.clone() }]
    {
        assert!(matches!(ImportanceExecution::restore(&original, changed, model(), &bytes), Err(UqCheckpointError::IdentityMismatch)));
    }
    assert!(matches!(ImportanceExecution::restore(&original, proposal, ContentHash([9; 32]), &bytes), Err(UqCheckpointError::IdentityMismatch)));
}

#[test]
fn malformed_framing_counts_status_and_weights_refuse_even_after_resealing() {
    let plan = plan(); let proposal = proposal();
    let mut run = ImportanceExecution::new(&plan, proposal.clone()).unwrap();
    run.advance(3, || false, response); let original = run.checkpoint(model()).unwrap();
    let restore = |bytes: &[u8]| ImportanceExecution::restore(&plan, proposal.clone(), model(), bytes);
    for n in [0, 7, 40, 80, original.len() - 1] { assert!(restore(&original[..n]).is_err()); }
    let mut bytes = original.clone(); bytes.push(0); assert!(restore(&bytes).is_err());
    let mut bytes = original.clone(); bytes[50] ^= 1;
    assert!(matches!(restore(&bytes), Err(UqCheckpointError::IntegrityMismatch)));
    for status in [2, 3, 255] {
        let mut bytes = original.clone(); bytes[40] = status; reseal(&mut bytes); assert!(restore(&bytes).is_err());
    }
    for count in [0_u64, 18, u64::MAX] {
        let mut bytes = original.clone(); bytes[41..49].copy_from_slice(&count.to_le_bytes()); reseal(&mut bytes);
        assert!(restore(&bytes).is_err());
    }
    for value in [0.0_f64, -0.0, -1.0, 1.01, f64::NAN, f64::INFINITY] {
        let mut bytes = original.clone(); bytes[57..65].copy_from_slice(&value.to_bits().to_le_bytes()); reseal(&mut bytes);
        assert!(restore(&bytes).is_err());
    }
    let mut bytes = original.clone(); bytes[49..57].copy_from_slice(&f64::NAN.to_bits().to_le_bytes()); reseal(&mut bytes);
    assert!(restore(&bytes).is_err());
    let mut zero = ImportanceExecution::new(&plan, ImportanceConfig { germ_shift: vec![0.0; 2], defensive_fraction: 0.2 }).unwrap();
    zero.advance(1, || false, response); let mut bytes = zero.checkpoint(model()).unwrap();
    bytes[57..65].copy_from_slice(&0.3_f64.to_bits().to_le_bytes()); reseal(&mut bytes);
    assert!(ImportanceExecution::restore(&plan, zero.config().clone(), model(), &bytes).is_err());
}

#[test]
fn terminal_failures_cannot_be_saved_as_selected_successful_prefixes() {
    for output in [Err("physical refusal"), Ok(f64::NAN), Ok(f64::INFINITY)] {
        let mut run = ImportanceExecution::new(&plan(), proposal()).unwrap();
        run.advance(3, || false, response);
        let report = run.advance(1, || false, |_| output);
        assert_eq!(report.status, UqStatus::Refused);
        assert_eq!(report.evaluations_attempted, 4);
        assert_eq!(run.checkpoint(model()), Err(UqCheckpointError::NotResumable));
        assert_eq!(run.advance(9, || panic!("refused"), |_| -> Result<f64, &str> { panic!("refused") }), report);
    }
}
