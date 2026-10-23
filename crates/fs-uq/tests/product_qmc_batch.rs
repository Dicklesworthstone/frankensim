//! G0/G4/G5: original nets, original reductions, drained workers and exact recovery.
use fs_blake3::{ContentHash, hash_domain};
use fs_uq::{CorrelationModel, ParameterUncertainty, PropagationMethod, QmcConfig,
    QmcExecution, QmcReport, UqCheckpointError, UqPlan, UqStatus};
use fs_uq::product_qmc::{QmcBatchExecution, QmcBatchReport};

fn plan() -> UqPlan {
    UqPlan::new("temperature", PropagationMethod::QuasiMonteCarlo, 32)
        .with_parameter(ParameterUncertainty::uniform("ambient", 290.0, 310.0, "K"))
        .with_parameter(ParameterUncertainty::gaussian("power", 30.0, 2.0, "W"))
        .with_correlation(CorrelationModel::Independent).with_compliance_threshold(330.0)
}
fn layout() -> QmcConfig { QmcConfig { replicates: 4, samples_per_replicate: 8 } }
fn model() -> ContentHash { ContentHash([73; 32]) }
fn qoi(x: &[f64]) -> f64 { x[0] + 0.8 * x[1] + 0.001 * x[1] * x[1] }
fn reference(plan: &UqPlan) -> (Vec<Vec<f64>>, QmcReport) {
    let mut inputs = Vec::new();
    let report = QmcExecution::new(plan, layout()).unwrap().advance(32, || false, |x| {
        inputs.push(x.to_vec()); Ok::<_, &str>(qoi(x))
    });
    (inputs, report)
}
fn finish(run: &mut QmcBatchExecution, width: usize) -> QmcBatchReport {
    loop {
        let report = run.advance_batch::<&str>(width, || false, |jobs, slots| {
            for (job, slot) in jobs.iter().zip(slots) { *slot = Some(Ok(qoi(&job.parameters))); }
        }).unwrap();
        if report.status == UqStatus::Complete { return report; }
        assert_eq!(report.status, UqStatus::BudgetTruncated);
        assert!(report.result.is_none());
    }
}

#[test]
fn scoped_workers_and_cross_replicate_waves_match_the_original_scalar_bits() {
    let mut singular = plan();
    singular.parameters = vec![ParameterUncertainty::gaussian("a", 300.0, 2.0, "K"),
        ParameterUncertainty::gaussian("b", 30.0, 1.0, "W")];
    singular.correlation = CorrelationModel::JointGaussian { matrix: vec![vec![1.0, 1.0], vec![1.0, 1.0]] };
    for p in [plan(), singular] {
        let (inputs, expected) = reference(&p);
        for width in [1, 3, 9, 64] {
            let mut run = QmcBatchExecution::new(&p, layout()).unwrap();
            loop {
                let report = run.advance_batch::<&str>(width, || false, |jobs, slots| {
                    for job in jobs { assert_eq!(job.parameters, inputs[job.ordinal]); }
                    std::thread::scope(|scope| {
                        let handles: Vec<_> = jobs.iter().enumerate().map(|(i, job)|
                            (i, scope.spawn(move || qoi(&job.parameters)))).collect();
                        for (i, handle) in handles.into_iter().rev() { slots[i] = Some(Ok(handle.join().unwrap())); }
                    });
                }).unwrap();
                if report.status == UqStatus::Complete {
                    assert_eq!(report.result, Some(expected.clone())); break;
                }
                assert!(report.result.is_none());
            }
        }
    }
}

#[test]
fn later_nets_survive_a_hole_and_only_the_missing_point_is_retried() {
    let p = plan();
    let (inputs, expected) = reference(&p);
    let mut run = QmcBatchExecution::new(&p, layout()).unwrap();
    let report = run.advance_batch::<&str>(12, || false, |jobs, slots| {
        for (job, slot) in jobs.iter().zip(slots) {
            if job.ordinal != 0 { *slot = Some(Ok(qoi(&job.parameters))); }
        }
    }).unwrap();
    assert_eq!(report.status, UqStatus::Cancelled);
    assert_eq!(report.evaluations_attempted, 11);
    assert_eq!(report.contiguous_observations, 0);
    assert_eq!(report.retained_out_of_order, 11);
    assert!(report.result.is_none());
    assert_eq!(run.observation(10), Some(qoi(&inputs[10])));
    let bytes = run.checkpoint(model()).unwrap();
    let mut run = QmcBatchExecution::restore(&p, layout(), model(), &bytes).unwrap();
    assert_eq!(run.report(), report);
    let report = run.advance_batch::<&str>(1, || false, |jobs, slots| {
        assert_eq!(jobs.len(), 1);
        assert_eq!(jobs[0].ordinal, 0);
        assert_eq!(jobs[0].parameters, inputs[0]);
        slots[0] = Some(Ok(qoi(&jobs[0].parameters)));
    }).unwrap();
    assert_eq!(report.contiguous_observations, 12);
    assert_eq!(report.evaluations_attempted, 12);
    assert!(report.result.is_none(), "complete leading nets do not license selected-design inference");
    assert_eq!(finish(&mut run, 3).result, Some(expected));
}

#[test]
fn every_prefix_and_cancellation_checkpoint_preserves_the_remaining_experiment() {
    let p = plan();
    let expected = reference(&p).1;
    for split in 0..=32 {
        let mut run = QmcBatchExecution::new(&p, layout()).unwrap();
        run.advance_batch::<&str>(split, || false, |jobs, slots| {
            for (job, slot) in jobs.iter().zip(slots) { *slot = Some(Ok(qoi(&job.parameters))); }
        }).unwrap();
        let saved = run.checkpoint(model()).unwrap();
        let mut restored = QmcBatchExecution::restore(&p, layout(), model(), &saved).unwrap();
        assert_eq!(restored.report(), run.report());
        assert_eq!(finish(&mut restored, 7).result, Some(expected.clone()));
    }
    for stop in 0..=10 {
        let mut run = QmcBatchExecution::new(&p, layout()).unwrap();
        let mut polls = 0;
        run.advance_batch::<&str>(8, || { let cancel = polls == stop; polls += 1; cancel }, |jobs, slots| {
            for (job, slot) in jobs.iter().zip(slots) { *slot = Some(Ok(qoi(&job.parameters))); }
        }).unwrap();
        let mut restored = QmcBatchExecution::restore(&p, layout(), model(), &run.checkpoint(model()).unwrap()).unwrap();
        assert_eq!(finish(&mut restored, 5).result, Some(expected.clone()));
    }
}

#[test]
fn failed_jobs_are_terminal_but_all_finite_siblings_are_retained() {
    let mut run = QmcBatchExecution::new(&plan(), layout()).unwrap();
    let report = run.advance_batch::<&str>(5, || false, |jobs, slots| {
        for (job, slot) in jobs.iter().zip(slots).rev() {
            *slot = Some(if job.ordinal == 1 { Err("first failure") }
                else if job.ordinal == 3 { Ok(f64::NAN) } else { Ok(qoi(&job.parameters)) });
        }
    }).unwrap();
    assert_eq!(report.status, UqStatus::Refused);
    assert_eq!(report.evaluations_attempted, 5);
    assert_eq!(report.evaluations_accepted, 3);
    assert!(report.rejection_reason.as_ref().unwrap().contains("evaluation 1 refused: first failure"));
    assert!(run.observation(4).is_some());
    assert!(report.result.is_none());
    assert_eq!(run.checkpoint(model()), Err(UqCheckpointError::NotResumable));
    assert_eq!(run.advance_batch::<&str>(1, || panic!("terminal"), |_, _| panic!("terminal")).unwrap(), report);
}

fn reseal(bytes: &mut [u8]) {
    let end = bytes.len() - 32;
    let checksum = hash_domain("org.frankensim.uq.qmc-batch.checkpoint.v1", &bytes[..end]);
    bytes[end..].copy_from_slice(checksum.as_bytes());
}

#[test]
fn changed_layout_identity_corruption_and_forged_pending_values_refuse() {
    let p = plan();
    let mut run = QmcBatchExecution::new(&p, layout()).unwrap();
    run.advance_batch::<&str>(4, || false, |_, slots| { slots[2] = Some(Ok(-0.0)); }).unwrap();
    let bytes = run.checkpoint(model()).unwrap();
    assert_eq!(QmcBatchExecution::restore(&p, QmcConfig { replicates: 8, samples_per_replicate: 4 }, model(), &bytes).unwrap_err(), UqCheckpointError::IdentityMismatch);
    assert_eq!(QmcBatchExecution::restore(&p, layout(), ContentHash([1; 32]), &bytes).unwrap_err(), UqCheckpointError::IdentityMismatch);
    let mut changed = p.clone(); changed.seed += 1;
    assert_eq!(QmcBatchExecution::restore(&changed, layout(), model(), &bytes).unwrap_err(), UqCheckpointError::IdentityMismatch);
    for end in 0..bytes.len() { assert!(QmcBatchExecution::restore(&p, layout(), model(), &bytes[..end]).is_err()); }
    let mut extra = bytes.clone(); extra.push(0);
    assert!(QmcBatchExecution::restore(&p, layout(), model(), &extra).is_err());
    let scalar_len = u64::from_le_bytes(bytes[8..16].try_into().unwrap()) as usize;
    let start = 24 + scalar_len;
    for (offset, data) in [(16, u64::MAX.to_le_bytes().to_vec()),
        (start, vec![1]), (start + 1, (-0.0_f64).to_bits().to_le_bytes().to_vec()),
        (start + 19, f64::NAN.to_bits().to_le_bytes().to_vec())]
    {
        let mut bad = bytes.clone(); bad[offset..offset + data.len()].copy_from_slice(&data); reseal(&mut bad);
        assert!(QmcBatchExecution::restore(&p, layout(), model(), &bad).is_err());
    }
    let restored = QmcBatchExecution::restore(&p, layout(), model(), &bytes).unwrap();
    assert_eq!(restored.observation(2).unwrap().to_bits(), (-0.0_f64).to_bits());
    let scalar = QmcExecution::new(&p, layout()).unwrap().checkpoint(model()).unwrap();
    assert!(QmcBatchExecution::restore(&p, layout(), model(), &scalar).is_err());
}

#[test]
fn bounded_admission_zero_work_and_extreme_values_keep_qmc_semantics() {
    let p = plan();
    let mut run = QmcBatchExecution::new(&p, layout()).unwrap();
    assert!(run.advance_batch::<&str>(257, || panic!("invalid"), |_, _| panic!("invalid")).is_err());
    assert_eq!(run.advance_batch::<&str>(0, || panic!("zero"), |_, _| panic!("zero")).unwrap(), run.report());
    let mut invalid = p.clone(); invalid.method = PropagationMethod::MonteCarlo;
    assert!(QmcBatchExecution::new(&invalid, layout()).is_err());
    for value in [-0.0_f64, f64::MAX, -f64::MAX] {
        let mut run = QmcBatchExecution::new(&p, layout()).unwrap();
        let result = run.advance_batch::<&str>(32, || false, |_, slots| {
            for slot in slots { *slot = Some(Ok(value)); }
        }).unwrap();
        let expected = QmcExecution::new(&p, layout()).unwrap().advance(32, || false, |_| Ok::<_, &str>(value));
        assert_eq!(result.result, Some(expected));
        let restored = QmcBatchExecution::restore(&p, layout(), model(), &run.checkpoint(model()).unwrap()).unwrap();
        assert!(restored.observations().iter().all(|x| x.to_bits() == value.to_bits()));
        assert_eq!(restored.report(), result);
    }
}
