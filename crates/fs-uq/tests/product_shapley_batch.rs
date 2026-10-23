//! G0/G4/G5: preserve the physical conditional experiment across scheduling.
use fs_blake3::{ContentHash, hash_domain};
use fs_uq::{CorrelationModel, ParameterUncertainty, PropagationMethod, ShapleyExecution,
    UqCheckpointError, UqPlan, UqStatus};
use fs_uq::product_shapley::ShapleyBatchExecution;
use std::cell::Cell;

const MODEL: ContentHash = ContentHash([29; 32]);
fn plan(dim: usize, rho: Option<f64>) -> UqPlan {
    let mut p = UqPlan::new("response", PropagationMethod::MonteCarlo, 4 * ((1 << dim) - 1));
    for i in 0..dim {
        p = p.with_parameter(ParameterUncertainty::gaussian(&format!("x{i}"), i as f64, 0.5, "K"));
    }
    p.correlation = rho.map_or(CorrelationModel::Independent, |rho| CorrelationModel::JointGaussian {
        matrix: (0..dim).map(|i| (0..dim).map(|j| if i == j { 1.0 } else { rho }).collect()).collect(),
    });
    p.seed = 73;
    p
}
fn qoi(x: &[f64]) -> f64 {
    x.iter().enumerate().map(|(i, x)| (i + 1) as f64 * x).sum::<f64>() + 0.125 * x[0] * x[0]
}
fn reference(p: &UqPlan) -> (ShapleyExecution, Vec<Vec<u64>>) {
    let mut run = ShapleyExecution::new(p).unwrap();
    let mut inputs = Vec::new();
    run.advance(usize::MAX, || false, |x| {
        inputs.push(x.iter().map(|v| v.to_bits()).collect()); Ok::<_, &str>(qoi(x))
    });
    (run, inputs)
}
fn finish(run: &mut ShapleyBatchExecution, width: usize) {
    while run.report().status != UqStatus::Complete {
        let report = run.advance_batch::<&str>(width, || false, |jobs, slots| {
            for i in (0..jobs.len()).rev() { slots[i] = Some(Ok(qoi(&jobs[i].parameters))); }
        }).unwrap();
        assert_ne!(report.status, UqStatus::Refused);
    }
}

#[test]
fn all_original_conditional_inputs_and_effects_survive_wave_and_row_boundaries() {
    for dim in [1, 2, 3, 8] {
        for rho in [None, Some(0.6), Some(1.0)] {
            let p = plan(dim, rho);
            let (serial, inputs) = reference(&p);
            for width in [1, 3, 7, 256] {
                let mut run = ShapleyBatchExecution::new(&p).unwrap();
                while run.report().status != UqStatus::Complete {
                    let report = run.advance_batch::<&str>(width, || false, |jobs, slots| {
                        for i in (0..jobs.len()).rev() {
                            let job = &jobs[i];
                            assert_eq!(job.parameters.iter().map(|v| v.to_bits()).collect::<Vec<_>>(), inputs[job.ordinal]);
                            slots[i] = Some(Ok(qoi(&job.parameters)));
                        }
                    }).unwrap();
                    assert_ne!(report.status, UqStatus::Refused);
                    if report.status != UqStatus::Complete { assert!(report.result.is_none()); }
                }
                assert_eq!(run.report().result.unwrap(), serial.report());
                assert_eq!(run.observations().iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
                    serial.observations().iter().map(|x| x.to_bits()).collect::<Vec<_>>());
                assert_eq!(run.evaluations_attempted(), p.budget_max_samples);
            }
        }
    }
}

#[test]
fn every_hole_pattern_retains_paid_halves_and_retries_only_missing_physics() {
    let p = plan(3, Some(0.6));
    let (serial, inputs) = reference(&p);
    for mask in 0_u32..128 {
        let mut run = ShapleyBatchExecution::new(&p).unwrap();
        let mut paid = vec![0; p.budget_max_samples];
        let report = run.advance_batch::<&str>(7, || false, |jobs, slots| {
            for (i, job) in jobs.iter().enumerate() {
                if mask & (1 << i) != 0 {
                    slots[i] = Some(Ok(qoi(&job.parameters))); paid[job.ordinal] += 1;
                }
            }
        }).unwrap();
        assert_eq!(report.evaluations_accepted, mask.count_ones() as usize);
        assert!(report.result.is_none());
        let bytes = run.checkpoint(MODEL).unwrap();
        let mut resumed = ShapleyBatchExecution::restore(&p, MODEL, &bytes).unwrap();
        assert_eq!(resumed.report(), report);
        assert_eq!(resumed.checkpoint(MODEL).unwrap(), bytes);
        for width in [1, 2, 5, 256].into_iter().cycle().take(p.budget_max_samples) {
            if resumed.report().status == UqStatus::Complete { break; }
            resumed.advance_batch::<&str>(width, || false, |jobs, slots| {
                for (job, slot) in jobs.iter().zip(slots) {
                    assert_eq!(paid[job.ordinal], 0, "repeated paid ordinal {}", job.ordinal);
                    assert_eq!(job.parameters.iter().map(|x| x.to_bits()).collect::<Vec<_>>(), inputs[job.ordinal]);
                    *slot = Some(Ok(qoi(&job.parameters))); paid[job.ordinal] += 1;
                }
            }).unwrap();
        }
        assert!(paid.iter().all(|&n| n == 1));
        assert_eq!(resumed.report().result.unwrap(), serial.report());
    }
}

#[test]
fn scoped_workers_overlap_and_publish_original_slots() {
    let p = plan(2, Some(0.6));
    let mut run = ShapleyBatchExecution::new(&p).unwrap();
    let barrier = std::sync::Barrier::new(4);
    run.advance_batch::<&str>(4, || false, |jobs, slots| {
        std::thread::scope(|scope| {
            let barrier = &barrier;
            let handles: Vec<_> = jobs.iter().map(|job| scope.spawn(move || {
                barrier.wait(); qoi(&job.parameters)
            })).collect();
            for (slot, handle) in slots.iter_mut().zip(handles) { *slot = Some(Ok(handle.join().unwrap())); }
        });
    }).unwrap();
    finish(&mut run, 3);
    assert_eq!(run.report().result.unwrap(), reference(&p).0.report());
}

#[test]
fn cancellation_at_each_wave_checkpoint_retains_every_completed_response() {
    let p = plan(3, Some(1.0));
    let expected = reference(&p).0.report();
    for stop in 1..=10 {
        let mut run = ShapleyBatchExecution::new(&p).unwrap();
        let polls = Cell::new(0);
        let calls = Cell::new(0);
        run.advance_batch::<&str>(7, || {
            polls.set(polls.get() + 1); polls.get() == stop
        }, |jobs, slots| {
            for (job, slot) in jobs.iter().zip(slots) {
                calls.set(calls.get() + 1); *slot = Some(Ok(qoi(&job.parameters)));
            }
        }).unwrap();
        assert_eq!(run.evaluations_attempted(), calls.get());
        assert!(run.report().result.is_none());
        let mut cloned = run.clone();
        let mut restored = ShapleyBatchExecution::restore(&p, MODEL, &run.checkpoint(MODEL).unwrap()).unwrap();
        finish(&mut cloned, 3); finish(&mut restored, 5);
        assert_eq!(cloned.report().result.unwrap(), expected);
        assert_eq!(restored.report().result.unwrap(), expected);
        assert_eq!(cloned.checkpoint(MODEL), restored.checkpoint(MODEL));
    }
}

#[test]
fn actual_failures_win_over_cancellation_and_cannot_yield_resumable_prefixes() {
    let p = plan(2, Some(0.6));
    for nonfinite in [None, Some(f64::NAN), Some(f64::INFINITY)] {
        let mut run = ShapleyBatchExecution::new(&p).unwrap();
        let cancel = Cell::new(false);
        let report = run.advance_batch::<&str>(5, || cancel.get(), |jobs, slots| {
            for i in (0..jobs.len()).rev() {
                slots[i] = Some(if i == 1 { nonfinite.map_or(Err("first failure"), Ok) }
                    else if i == 3 { Err("later failure") } else { Ok(qoi(&jobs[i].parameters)) });
            }
            cancel.set(true);
        }).unwrap();
        assert_eq!(report.status, UqStatus::Refused);
        assert_eq!(report.evaluations_attempted, 5);
        assert_eq!(report.evaluations_accepted, 3);
        assert!(run.observation(4).is_some());
        assert!(report.rejection_reason.unwrap().contains("evaluation 1"));
        assert!(report.result.is_none());
        assert_eq!(run.checkpoint(MODEL), Err(UqCheckpointError::NotResumable));
        run.advance_batch::<&str>(1, || panic!("terminal poll"), |_, _| panic!("terminal driver")).unwrap();
    }
}

#[test]
fn checkpoint_binds_the_physical_law_and_rejects_malformed_pending_responses() {
    let p = plan(2, Some(0.6));
    let mut run = ShapleyBatchExecution::new(&p).unwrap();
    run.advance_batch::<&str>(3, || false, |jobs, slots| {
        slots[1] = Some(Ok(qoi(&jobs[1].parameters)));
        slots[2] = Some(Ok(qoi(&jobs[2].parameters)));
    }).unwrap();
    let bytes = run.checkpoint(MODEL).unwrap();
    assert_eq!(&bytes[..8], b"FSSHPB01");
    assert!(ShapleyBatchExecution::restore(&p, ContentHash([30; 32]), &bytes).is_err());
    for changed in [plan(2, None), plan(2, Some(1.0))] {
        assert!(ShapleyBatchExecution::restore(&changed, MODEL, &bytes).is_err());
    }
    let mut changed = p.clone(); changed.seed += 1;
    assert!(ShapleyBatchExecution::restore(&changed, MODEL, &bytes).is_err());
    for end in 0..bytes.len() { assert!(ShapleyBatchExecution::restore(&p, MODEL, &bytes[..end]).is_err()); }
    let mut extra = bytes.clone(); extra.push(0);
    assert!(ShapleyBatchExecution::restore(&p, MODEL, &extra).is_err());
    for at in [0, 8, 16, 30, bytes.len() - 1] {
        let mut bad = bytes.clone(); bad[at] ^= 0x80;
        assert!(ShapleyBatchExecution::restore(&p, MODEL, &bad).is_err());
    }
    let inner_len = u64::from_le_bytes(bytes[8..16].try_into().unwrap()) as usize;
    let first = 24 + inner_len;
    for (at, value) in [(first, vec![1]), (first + 1, (-0.0_f64).to_bits().to_le_bytes().to_vec()),
        (first + 10, f64::INFINITY.to_bits().to_le_bytes().to_vec())] {
        let mut bad = bytes.clone(); bad[at..at + value.len()].copy_from_slice(&value);
        let end = bad.len() - 32;
        let hash = hash_domain("org.frankensim.uq.shapley-batch.checkpoint.v1", &bad[..end]);
        bad[end..].copy_from_slice(hash.as_bytes());
        assert!(ShapleyBatchExecution::restore(&p, MODEL, &bad).is_err());
    }
    let serial = ShapleyExecution::new(&p).unwrap().checkpoint(MODEL).unwrap();
    assert!(ShapleyBatchExecution::restore(&p, MODEL, &serial).is_err());
    assert!(ShapleyExecution::restore(&p, MODEL, &bytes).is_err());
}

#[test]
fn admission_and_undefined_effects_keep_original_semantics() {
    assert!(ShapleyBatchExecution::new(&plan(9, None)).is_err());
    let mut p = plan(2, None); p.budget_max_samples -= 1;
    assert!(ShapleyBatchExecution::new(&p).is_err());
    let p = plan(2, Some(1.0));
    let mut run = ShapleyBatchExecution::new(&p).unwrap();
    assert!(run.advance_batch::<&str>(257, || panic!("invalid poll"), |_, _| panic!("invalid driver")).is_err());
    run.advance_batch::<&str>(0, || panic!("zero poll"), |_, _| panic!("zero driver")).unwrap();
    for value in [-0.0_f64, f64::MAX] {
        let mut run = run.clone();
        run.advance_batch::<&str>(256, || false, |_, slots| slots.fill(Some(Ok(value)))).unwrap();
        let restored = ShapleyBatchExecution::restore(&p, MODEL, &run.checkpoint(MODEL).unwrap()).unwrap();
        assert!(restored.observations().iter().all(|v| v.to_bits() == value.to_bits()));
        let report = restored.report().result.unwrap();
        assert_eq!(report.status, UqStatus::Complete);
        assert!(report.estimate.is_none()); assert!(report.unavailable_reason.is_some());
        run.advance_batch::<&str>(1, || panic!("complete poll"), |_, _| panic!("complete driver")).unwrap();
    }
}
