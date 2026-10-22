//! Recovery of the original learning experiment, including phase boundaries.
use fs_blake3::{ContentHash, hash_domain};
use fs_uq::{CorrelationModel, ImportanceConfig, ImportanceFitConfig, ParameterUncertainty,
    PilotImportanceConfig, PilotImportanceExecution, UqCheckpointError, UqPlan,
    PropagationMethod, UqStatus};

fn plan(singular: bool) -> UqPlan {
    let mut plan = UqPlan::new("response", PropagationMethod::MonteCarlo, 24)
        .with_parameter(ParameterUncertainty::gaussian("load", 0.0, 1.0, "N"))
        .with_parameter(ParameterUncertainty::gaussian("temperature", 300.0, 0.5, "K"))
        .with_correlation(if singular { CorrelationModel::JointGaussian {
            matrix: vec![vec![1.0, 1.0], vec![1.0, 1.0]],
        }} else { CorrelationModel::Independent })
        .with_compliance_threshold(2.0);
    plan.seed = u64::MAX - 3;
    plan
}
fn config() -> PilotImportanceConfig {
    PilotImportanceConfig {
        pilot_rounds: 2, samples_per_pilot: 6,
        initial_proposal: ImportanceConfig { germ_shift: vec![0.25, -0.0], defensive_fraction: 0.25 },
        fit: ImportanceFitConfig { elite_fraction: 0.5, relaxation: 0.7, maximum_shift_norm: 4.0 },
    }
}
fn model() -> ContentHash { hash_domain("test-pilot-recovery", b"load+0.5*(temperature-300)") }
fn response(x: &[f64]) -> Result<f64, &'static str> { Ok(x[0] + 0.5 * (x[1] - 300.0)) }
fn new(singular: bool) -> PilotImportanceExecution {
    PilotImportanceExecution::new(&plan(singular), config()).unwrap()
}
fn restore(run: &PilotImportanceExecution) -> PilotImportanceExecution {
    PilotImportanceExecution::restore(run.plan(), run.config().clone(), model(),
        &run.checkpoint(model()).unwrap()).unwrap()
}

#[test]
fn every_paid_prefix_restores_exact_proposals_weights_and_final_report() {
    for singular in [false, true] {
        let mut full = new(singular);
        let mut inputs = Vec::new();
        let expected = full.advance(usize::MAX, || false, |x| {
            inputs.push(x.to_vec()); response(x)
        });
        assert_eq!(expected.status, UqStatus::Complete);
        for prefix in 0..=24 {
            let mut run = new(singular);
            run.advance(prefix, || false, response);
            let bytes = run.checkpoint(model()).unwrap();
            let mut resumed = restore(&run);
            assert_eq!(resumed.report(), run.report());
            assert_eq!(resumed.checkpoint(model()).unwrap(), bytes);
            let mut next = prefix;
            let actual = resumed.advance(usize::MAX, || false, |x| {
                assert_eq!(x, inputs[next]); next += 1; response(x)
            });
            assert_eq!(next, 24);
            assert_eq!(actual, expected);
            assert_eq!(resumed.checkpoint(model()).unwrap(), full.checkpoint(model()).unwrap());
        }
    }
}

#[test]
fn cancellation_during_fits_and_after_last_call_never_repeats_physics() {
    let mut full = new(false);
    let mut checkpoints = 0;
    let expected = full.advance(usize::MAX, || { checkpoints += 1; false }, response);
    let mut saw_pending_fit = false;
    let mut saw_pending_publication = false;
    for stop in 1..=checkpoints {
        let mut run = new(false);
        let mut polls = 0;
        let mut calls = 0;
        let report = run.advance(usize::MAX, || { polls += 1; polls == stop }, |x| {
            calls += 1; response(x)
        });
        assert_eq!(report.status, UqStatus::Cancelled);
        saw_pending_fit |= run.production().is_none() && run.current().report().status == UqStatus::Complete;
        saw_pending_publication |= calls == 24;
        let mut resumed = restore(&run);
        assert_eq!(resumed.report(), report);
        let actual = resumed.advance(usize::MAX, || false, |x| { calls += 1; response(x) });
        assert_eq!(calls, 24, "repeated paid work at checkpoint {stop}");
        assert_eq!(actual, expected);
    }
    assert!(saw_pending_fit && saw_pending_publication);
}

#[test]
fn unfinished_call_retries_identical_input_after_disk_restore() {
    let mut run = new(true);
    run.advance(7, || false, response);
    let mut unfinished = Vec::new();
    run.advance_interruptible(1, || false, |x| {
        unfinished = x.to_vec(); Ok::<_, &str>(None)
    });
    let mut resumed = restore(&run);
    assert_eq!(resumed.evaluations_attempted(), 7);
    let report = resumed.advance(1, || false, |x| {
        assert_eq!(x, unfinished); response(x)
    });
    assert_eq!(report.model_calls_evaluated, 8);
    assert!(report.estimate.is_none());
}

#[test]
fn changed_model_plan_and_every_fitting_choice_refuse() {
    let run = new(false);
    let bytes = run.checkpoint(model()).unwrap();
    let mut variants = Vec::new();
    let mut changed = config(); changed.pilot_rounds = 1; variants.push(changed);
    let mut changed = config(); changed.samples_per_pilot = 5; variants.push(changed);
    let mut changed = config(); changed.initial_proposal.defensive_fraction = 0.5; variants.push(changed);
    let mut changed = config(); changed.initial_proposal.germ_shift[0] = -0.25; variants.push(changed);
    let mut changed = config(); changed.fit.elite_fraction = 0.25; variants.push(changed);
    let mut changed = config(); changed.fit.relaxation = 1.0; variants.push(changed);
    let mut changed = config(); changed.fit.maximum_shift_norm = 5.0; variants.push(changed);
    for policy in variants {
        assert!(matches!(PilotImportanceExecution::restore(&plan(false), policy, model(), &bytes),
            Err(UqCheckpointError::IdentityMismatch)));
    }
    let mut changed = plan(false); changed.seed = 7;
    assert!(PilotImportanceExecution::restore(&changed, config(), model(), &bytes).is_err());
    assert!(PilotImportanceExecution::restore(&plan(true), config(), model(), &bytes).is_err());
    assert!(PilotImportanceExecution::restore(&plan(false), config(), ContentHash([0; 32]), &bytes).is_err());
    for prefix in 0..bytes.len() {
        assert!(PilotImportanceExecution::restore(&plan(false), config(), model(), &bytes[..prefix]).is_err());
    }
}

fn reseal(bytes: &mut [u8]) {
    let end = bytes.len() - 32;
    let sum = hash_domain("org.frankensim.uq.pilot-importance.checkpoint.v1", &bytes[..end]);
    bytes[end..].copy_from_slice(sum.as_bytes());
}

#[test]
fn malformed_history_and_phase_cannot_change_work_accounting() {
    let mut run = new(false); run.advance(6, || false, response);
    assert_eq!(run.pilot_history().len(), 1);
    let original = run.checkpoint(model()).unwrap();
    for (offset, replacement) in [
        (40, vec![2]), (41, u64::MAX.to_le_bytes().to_vec()),
        (49, 0_u64.to_le_bytes().to_vec()), (73, vec![2]),
        (74, f64::NAN.to_bits().to_le_bytes().to_vec()),
        (91, 1e100_f64.to_bits().to_le_bytes().to_vec()),
        (107, u64::MAX.to_le_bytes().to_vec()),
    ] {
        let mut bytes = original.clone();
        bytes[offset..offset + replacement.len()].copy_from_slice(&replacement);
        reseal(&mut bytes);
        assert!(PilotImportanceExecution::restore(run.plan(), config(), model(), &bytes).is_err());
    }
    let mut corrupt = original.clone(); corrupt[50] ^= 1;
    assert!(matches!(PilotImportanceExecution::restore(run.plan(), config(), model(), &corrupt),
        Err(UqCheckpointError::IntegrityMismatch)));
    let mut trailing = original; trailing.push(0);
    assert!(PilotImportanceExecution::restore(run.plan(), config(), model(), &trailing).is_err());
}

#[test]
fn refused_runs_stay_terminal_and_constant_pilot_ties_round_trip() {
    let mut failed = new(false);
    failed.advance(3, || false, |_| Err::<f64, _>("physical refusal"));
    assert_eq!(failed.checkpoint(model()), Err(UqCheckpointError::NotResumable));
    let mut policy = config(); policy.initial_proposal.defensive_fraction = 1.0;
    let mut run = PilotImportanceExecution::new(&plan(false), policy).unwrap();
    run.advance(12, || false, |_| Ok::<_, &str>(2.0));
    assert!(run.pilot_history().iter().all(|p| p.fit.selected == 6 && p.fit.observed_failures == 0));
    let mut resumed = restore(&run);
    let expected = run.advance(12, || false, |_| Ok::<_, &str>(2.0));
    assert_eq!(resumed.advance(12, || false, |_| Ok::<_, &str>(2.0)), expected);
    assert_eq!(resumed.checkpoint(model()), run.checkpoint(model()));
}
