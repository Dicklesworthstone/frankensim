//! G3/G4/G5: original sampling, ordered reductions and interrupted physical work.
use fs_blake3::{ContentHash, hash_domain};
use fs_uq::{CorrelationModel, ParameterUncertainty, PropagationMethod, SobolExecution,
    UqCheckpointError, UqPlan, UqStatus};
use fs_uq::product_sensitivity::SobolBatchExecution;

fn plan() -> UqPlan {
    let mut p = UqPlan::new("response", PropagationMethod::MonteCarlo, 32)
        .with_parameter(ParameterUncertainty::uniform("ambient", -2.0, 3.0, "K"))
        .with_parameter(ParameterUncertainty::gaussian("load", 4.0, 0.5, "W"))
        .with_correlation(CorrelationModel::Independent);
    p.seed = 147;
    p
}
fn identity() -> ContentHash { ContentHash([37; 32]) }
fn physics(x: &[f64]) -> f64 { x[0] + 0.5 * x[1] * x[1] + x[0] * x[1] }
fn reference() -> (SobolExecution, Vec<Vec<f64>>) {
    let mut run = SobolExecution::new(&plan()).unwrap();
    let mut inputs = Vec::new();
    run.advance(32, || false, |x| { inputs.push(x.to_vec()); Ok::<_, &str>(physics(x)) });
    (run, inputs)
}
fn finish(run: &mut SobolBatchExecution, width: usize, inputs: &[Vec<f64>], paid: &mut [usize]) {
    while run.report().status != UqStatus::Complete {
        let report = run.advance_batch(width, || false, |jobs, slots| {
            for (job, slot) in jobs.iter().zip(slots).rev() {
                assert_eq!(job.parameters, inputs[job.ordinal]);
                paid[job.ordinal] += 1;
                *slot = Some(Ok::<_, &str>(physics(&job.parameters)));
            }
        }).unwrap();
        assert_ne!(report.status, UqStatus::Refused);
        if report.status != UqStatus::Complete { assert!(report.result.is_none()); }
    }
}

#[test]
fn every_sample_and_final_index_match_serial_across_row_and_wave_boundaries() {
    let (serial, inputs) = reference();
    for width in [1, 3, 4, 7, 31, 64, 256] {
        let mut run = SobolBatchExecution::new(&plan()).unwrap();
        let mut paid = vec![0; 32];
        finish(&mut run, width, &inputs, &mut paid);
        assert_eq!(paid, vec![1; 32]);
        assert_eq!(run.observations(), serial.observations());
        assert_eq!(run.report().result.as_ref(), Some(&serial.report()));
        assert_eq!(run.report().retained_out_of_order, 0);
        assert_eq!(run.report().completed_rows, 8);
    }
}

#[test]
fn all_gap_patterns_restore_without_repeating_successful_hybrids() {
    let (serial, inputs) = reference();
    for mask in 0_u32..128 {
        let mut run = SobolBatchExecution::new(&plan()).unwrap();
        let mut paid = vec![0; 32];
        run.advance_batch(7, || false, |jobs, slots| {
            for (job, slot) in jobs.iter().zip(slots) {
                assert_eq!(job.parameters, inputs[job.ordinal]);
                if mask & (1 << job.ordinal) != 0 {
                    paid[job.ordinal] += 1;
                    *slot = Some(Ok::<_, &str>(physics(&job.parameters)));
                }
            }
        }).unwrap();
        assert_eq!(run.evaluations_attempted(), mask.count_ones() as usize);
        assert!(run.report().result.is_none());
        let bytes = run.checkpoint(identity()).unwrap();
        let mut restored = SobolBatchExecution::restore(&plan(), identity(), &bytes).unwrap();
        assert_eq!(restored.report(), run.report());
        assert_eq!(restored.checkpoint(identity()).unwrap(), bytes);
        for i in 0..7 { assert_eq!(restored.observation(i), run.observation(i)); }
        let mut cloned = restored.clone();
        let mut cloned_paid = paid.clone();
        finish(&mut cloned, 5, &inputs, &mut cloned_paid);
        finish(&mut restored, 2, &inputs, &mut paid);
        assert_eq!(paid, vec![1; 32]);
        assert_eq!(cloned_paid, vec![1; 32]);
        assert_eq!(restored.report().result.as_ref(), Some(&serial.report()));
        assert_eq!(restored.checkpoint(identity()), cloned.checkpoint(identity()));
    }
}

#[test]
fn scoped_workers_overlap_and_cancellation_never_discards_drained_results() {
    let (serial, inputs) = reference();
    for stop in 1..=9 {
        let mut run = SobolBatchExecution::new(&plan()).unwrap();
        let mut polls = 0;
        let mut paid = vec![0; 32];
        run.advance_batch(3, || { polls += 1; polls == stop }, |jobs, slots| {
            let barrier = std::sync::Barrier::new(jobs.len());
            std::thread::scope(|scope| {
                let handles: Vec<_> = jobs.iter().map(|job| {
                    let barrier = &barrier;
                    scope.spawn(move || { barrier.wait(); physics(&job.parameters) })
                }).collect();
                for ((job, slot), handle) in jobs.iter().zip(slots).zip(handles) {
                    paid[job.ordinal] += 1;
                    *slot = Some(Ok::<_, &str>(handle.join().unwrap()));
                }
            });
        }).unwrap();
        assert_eq!(run.evaluations_attempted(), paid.iter().sum::<usize>());
        let mut restored = SobolBatchExecution::restore(&plan(), identity(), &run.checkpoint(identity()).unwrap()).unwrap();
        finish(&mut restored, 7, &inputs, &mut paid);
        assert_eq!(paid, vec![1; 32]);
        assert_eq!(restored.report().result.as_ref(), Some(&serial.report()));
    }
}

#[test]
fn actual_failures_are_terminal_and_lowest_ordinal_wins_after_all_successes_are_retained() {
    let mut run = SobolBatchExecution::new(&plan()).unwrap();
    let report = run.advance_batch(4, || false, |_, slots| {
        slots[3] = Some(Ok::<_, &str>(9.0));
        slots[2] = Some(Err("later failure"));
        slots[1] = Some(Err("first failure"));
        slots[0] = Some(Ok(7.0));
    }).unwrap();
    assert_eq!(report.status, UqStatus::Refused);
    assert!(report.rejection_reason.unwrap().contains("evaluation 1 refused: first failure"));
    assert_eq!(report.evaluations_attempted, 4);
    assert_eq!(report.evaluations_accepted, 2);
    assert_eq!(run.observation(3), Some(9.0));
    assert!(report.result.is_none());
    assert_eq!(run.checkpoint(identity()), Err(UqCheckpointError::NotResumable));
    run.advance_batch::<&str>(4, || panic!("terminal poll"), |_, _| panic!("terminal dispatch")).unwrap();
    let mut nonfinite = SobolBatchExecution::new(&plan()).unwrap();
    assert_eq!(nonfinite.advance_batch::<&str>(1, || false, |_, slots| {
        slots[0] = Some(Ok(f64::NAN));
    }).unwrap().status, UqStatus::Refused);
}

#[test]
fn checkpoint_binds_the_original_design_and_rejects_noncanonical_pending_values() {
    let mut run = SobolBatchExecution::new(&plan()).unwrap();
    run.advance_batch::<&str>(3, || false, |_, slots| { slots[1] = Some(Ok(-0.0)); }).unwrap();
    let bytes = run.checkpoint(identity()).unwrap();
    assert_eq!(&bytes[..8], b"FSSOBB01");
    let restored = SobolBatchExecution::restore(&plan(), identity(), &bytes).unwrap();
    assert_eq!(restored.observation(1).unwrap().to_bits(), (-0.0_f64).to_bits());
    let mut changed = plan(); changed.seed += 1;
    assert!(SobolBatchExecution::restore(&changed, identity(), &bytes).is_err());
    assert!(SobolBatchExecution::restore(&plan(), ContentHash([0; 32]), &bytes).is_err());
    let serial = SobolExecution::new(&plan()).unwrap().checkpoint(identity()).unwrap();
    assert!(SobolBatchExecution::restore(&plan(), identity(), &serial).is_err());
    assert!(SobolExecution::restore(&plan(), identity(), &bytes).is_err());
    for end in [0, 8, 23, bytes.len() - 1] {
        assert!(SobolBatchExecution::restore(&plan(), identity(), &bytes[..end]).is_err());
    }
    let pending = 24 + u64::from_le_bytes(bytes[8..16].try_into().unwrap()) as usize;
    for malformed in [0, 1] {
        let mut bad = bytes.clone();
        if malformed == 0 { bad[pending] = 1; }
        else { bad[pending + 10..pending + 18].copy_from_slice(&f64::NAN.to_bits().to_le_bytes()); }
        let end = bad.len() - 32;
        let checksum = hash_domain("org.frankensim.uq.sobol-batch.checkpoint.v1", &bad[..end]);
        bad[end..].copy_from_slice(checksum.as_bytes());
        assert!(SobolBatchExecution::restore(&plan(), identity(), &bad).is_err());
    }
}

#[test]
fn unchanged_admission_and_undefined_variance_do_not_invent_sensitivity() {
    let mut correlated = plan(); correlated.correlation = CorrelationModel::Unknown;
    assert!(SobolBatchExecution::new(&correlated).is_err());
    let mut wrong_count = plan(); wrong_count.budget_max_samples = 31;
    assert!(SobolBatchExecution::new(&wrong_count).is_err());
    let mut run = SobolBatchExecution::new(&plan()).unwrap();
    run.advance_batch::<&str>(0, || panic!("zero poll"), |_, _| panic!("zero dispatch")).unwrap();
    assert!(run.advance_batch::<&str>(257, || false, |_, _| panic!("invalid dispatch")).is_err());
    let report = run.advance_batch::<&str>(32, || false, |_, slots| {
        for slot in slots { *slot = Some(Ok(f64::MAX)); }
    }).unwrap().result.unwrap();
    assert_eq!(report.status, UqStatus::Complete);
    assert!(report.estimate.is_none());
    assert!(report.unavailable_reason.is_some());
    let mut restored = SobolBatchExecution::restore(&plan(), identity(), &run.checkpoint(identity()).unwrap()).unwrap();
    assert_eq!(restored.report().result, Some(report));
    restored.advance_batch::<&str>(1, || panic!("complete poll"), |_, _| panic!("complete dispatch")).unwrap();
}
