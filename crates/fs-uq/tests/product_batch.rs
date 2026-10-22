use fs_blake3::{ContentHash, hash_domain};
use fs_uq::{CorrelationModel, ParameterUncertainty, PropagationMethod, UqBatchExecution,
    UqCheckpointError, UqExecution, UqPlan, UqStatus};

fn plan(joint: bool) -> UqPlan {
    let mut p = UqPlan::new("temperature", PropagationMethod::MonteCarlo, 19)
        .with_parameter(ParameterUncertainty::gaussian("load", 3.0, 0.5, "W"))
        .with_parameter(if joint { ParameterUncertainty::gaussian("inlet", 300.0, 0.2, "K") }
            else { ParameterUncertainty::uniform("inlet", 299.0, 301.0, "K") })
        .with_correlation(if joint { CorrelationModel::JointGaussian { matrix: vec![vec![1.0,1.0],vec![1.0,1.0]] } }
            else { CorrelationModel::Independent }).with_compliance_threshold(307.0);
    p.seed = 73;
    p
}
fn model(x: &[f64]) -> f64 { x[1] + 2.0 * x[0] }
fn identity() -> ContentHash { hash_domain("test-batch-model", b"T=inlet+2*load") }
fn finish(run: &mut UqBatchExecution, width: usize) {
    while run.report().status != UqStatus::Complete {
        let report = run.advance_batch::<&str>(width, || false, |jobs, slots| {
            for (job, slot) in jobs.iter().zip(slots).rev() { *slot = Some(Ok(model(&job.parameters))); }
        }).unwrap();
        assert_ne!(report.status, UqStatus::Refused, "{report:?}");
    }
}

#[test]
fn widths_and_actual_scoped_workers_preserve_scalar_sample_and_report_bits() {
    for joint in [false, true] {
        let p = plan(joint);
        let mut scalar = UqExecution::new(&p).unwrap();
        scalar.advance(19, || false, |x| Ok::<_, &str>(model(x)));
        for width in [1, 2, 3, 8, 256] {
            let mut run = UqBatchExecution::new(&p).unwrap();
            while run.report().status != UqStatus::Complete {
                run.advance_batch::<&str>(width, || false, |jobs, slots| {
                    std::thread::scope(|scope| {
                        for (job, slot) in jobs.iter().zip(slots) {
                            scope.spawn(move || { *slot = Some(Ok(model(&job.parameters))); });
                        }
                    });
                }).unwrap();
            }
            assert_eq!(run.observations(), scalar.observations());
            assert_eq!(run.report().result, Some(scalar.report()));
            assert_eq!(run.evaluations_attempted(), 19);
        }
    }
}

#[test]
fn later_successes_survive_process_recovery_and_only_holes_are_retried() {
    let p = plan(true);
    let mut run = UqBatchExecution::new(&p).unwrap();
    let mut expected = Vec::new();
    let r = run.advance_batch::<&str>(8, || false, |jobs, slots| {
        expected = jobs.iter().map(|j| j.parameters.clone()).collect();
        for (job, slot) in jobs.iter().zip(slots).skip(1).rev() { *slot = Some(Ok(model(&job.parameters))); }
    }).unwrap();
    assert_eq!(r.status, UqStatus::Cancelled);
    assert_eq!(r.evaluations_accepted, 7);
    assert_eq!(r.contiguous_observations, 0);
    assert_eq!(r.retained_out_of_order, 7);
    assert!(r.result.is_none());
    let saved = run.checkpoint(identity()).unwrap();
    let mut restored = UqBatchExecution::restore(&p, identity(), &saved).unwrap();
    assert_eq!(restored.report(), run.report());
    assert_eq!(restored.checkpoint(identity()).unwrap(), saved);
    restored.advance_batch::<&str>(2, || false, |jobs, slots| {
        assert_eq!(jobs.len(), 1);
        assert_eq!(jobs[0].ordinal, 0);
        assert_eq!(jobs[0].parameters, expected[0]);
        slots[0] = Some(Ok(model(&jobs[0].parameters)));
    }).unwrap();
    assert_eq!(restored.evaluations_attempted(), 8);
    finish(&mut restored, 3);
    finish(&mut run, 7);
    assert_eq!(restored.report(), run.report());
    assert_eq!(restored.checkpoint(identity()).unwrap(), run.checkpoint(identity()).unwrap());
}

#[test]
fn every_saved_prefix_and_cancellation_checkpoint_keep_paid_observations() {
    let p = plan(false);
    let mut full = UqBatchExecution::new(&p).unwrap(); finish(&mut full, 5);
    for stop in 0..=19 {
        let mut run = UqBatchExecution::new(&p).unwrap();
        run.advance_batch::<&str>(stop, || false, |jobs, slots| {
            for (job, slot) in jobs.iter().zip(slots) { *slot = Some(Ok(model(&job.parameters))); }
        }).unwrap();
        let mut copy = UqBatchExecution::restore(&p, identity(), &run.checkpoint(identity()).unwrap()).unwrap();
        finish(&mut copy, 4);
        assert_eq!(copy.report(), full.report());
    }
    // A four-input wave has polls before setup, each input, dispatch, and after drain.
    for cancel_at in 0..7 {
        let mut run = UqBatchExecution::new(&p).unwrap();
        let mut polls = 0;
        let mut called = 0;
        let report = run.advance_batch::<&str>(4, || {
            let now = polls; polls += 1; now == cancel_at
        }, |jobs, slots| {
            for (job, slot) in jobs.iter().zip(slots) {
                called += 1; *slot = Some(Ok(model(&job.parameters)));
            }
        }).unwrap();
        assert_eq!(report.evaluations_accepted, called);
        let mut copy = UqBatchExecution::restore(&p, identity(), &run.checkpoint(identity()).unwrap()).unwrap();
        finish(&mut copy, 2);
        assert_eq!(copy.report(), full.report());
    }
}

#[test]
fn genuine_failures_are_terminal_even_with_successful_later_siblings() {
    let mut run = UqBatchExecution::new(&plan(false)).unwrap();
    let report = run.advance_batch(4, || false, |_, slots| {
        slots[3] = Some(Ok(99.0));
        slots[2] = Some(Err("later failure"));
        slots[0] = Some(Err("first failure"));
    }).unwrap();
    assert_eq!(report.status, UqStatus::Refused);
    assert_eq!(report.evaluations_attempted, 3);
    assert_eq!(report.evaluations_accepted, 1);
    assert_eq!(run.observation(3), Some(99.0));
    assert!(report.rejection_reason.unwrap().contains("evaluation 0 refused: first failure"));
    assert!(report.result.is_none());
    assert_eq!(run.checkpoint(identity()), Err(UqCheckpointError::NotResumable));
    run.advance_batch::<&str>(4, || false, |_, _| panic!("refused run dispatched")).unwrap();
    let mut nonfinite = UqBatchExecution::new(&plan(false)).unwrap();
    let r = nonfinite.advance_batch::<&str>(1, || false, |_, slots| slots[0] = Some(Ok(f64::NAN))).unwrap();
    assert_eq!(r.status, UqStatus::Refused);
}

#[test]
fn zero_terminal_and_invalid_allowances_do_not_dispatch() {
    let mut run = UqBatchExecution::new(&plan(false)).unwrap();
    assert!(run.advance_batch::<&str>(257, || false, |_, _| panic!("oversized wave")).is_err());
    run.advance_batch::<&str>(0, || true, |_, _| panic!("zero wave")).unwrap();
    assert_eq!(run.evaluations_attempted(), 0);
    finish(&mut run, 4);
    run.advance_batch::<&str>(4, || true, |_, _| panic!("complete wave")).unwrap();
}

#[test]
fn checkpoint_identity_corruption_and_noncanonical_slots_refuse() {
    let p = plan(false);
    let mut run = UqBatchExecution::new(&p).unwrap();
    run.advance_batch::<&str>(3, || false, |_, slots| slots[2] = Some(Ok(-0.0))).unwrap();
    let bytes = run.checkpoint(identity()).unwrap();
    let restored = UqBatchExecution::restore(&p, identity(), &bytes).unwrap();
    assert_eq!(restored.observation(2).unwrap().to_bits(), (-0.0_f64).to_bits());
    let mut changed = p.clone(); changed.seed += 1;
    assert!(matches!(UqBatchExecution::restore(&changed, identity(), &bytes), Err(UqCheckpointError::IdentityMismatch)));
    assert!(matches!(UqBatchExecution::restore(&p, ContentHash([1;32]), &bytes), Err(UqCheckpointError::IdentityMismatch)));
    for offset in [0, 8, 40, 41, 49, 57, bytes.len()-1] {
        let mut corrupt = bytes.clone(); corrupt[offset] ^= 0x80;
        assert!(UqBatchExecution::restore(&p, identity(), &corrupt).is_err());
    }
    // A resealed full-prefix first pending value must still fail canonical framing.
    let mut malformed = bytes.clone(); malformed[57] = 1;
    let end = malformed.len()-32;
    let digest = hash_domain("org.frankensim.uq.batch.checkpoint.v1", &malformed[..end]);
    malformed[end..].copy_from_slice(digest.as_bytes());
    assert!(UqBatchExecution::restore(&p, identity(), &malformed).is_err());
}
