use super::*;
use crate::{Kernel, Matern};
use crate::noisy::minimize_noisy;

fn config() -> NoisyBoConfig {
    NoisyBoConfig { bounds: (0.0, 1.0),
        kernel: Kernel { family: Matern::FiveHalves, signal: 1.0, lengthscales: vec![0.3] },
        prior_mean: 0.5, q: 2, mc_samples: 16, acq_starts: 1, acq_evals: 16, seed: 97 }
}

fn observation(x: &[f64]) -> NoisyObservation {
    NoisyObservation { value: (x[0] - 0.37).powi(2), noise_variance: 0.01 + 0.03 * x[0] }
}

fn finish(study: &mut NoisyStudy, reverse: bool) -> usize {
    let mut physical_calls = 0;
    while !study.is_complete() {
        let mut batch = study.ask();
        if reverse { batch.reverse(); }
        for request in batch {
            let result = observation(&request.x);
            physical_calls += 1;
            assert_eq!(study.tell(&request, result), Ok(TellDisposition::Recorded));
        }
        study.advance().unwrap();
    }
    physical_calls
}

#[test]
fn g3_out_of_order_results_match_synchronous_sequential_and_batched_bo() {
    for q in [1, 2] {
        let mut c = config(); c.q = q;
        let fixed = minimize_noisy(&mut observation, 1, 3, 2, &c);
        let mut study = NoisyStudy::new(42, 1, 3, 2, &c);
        assert_eq!(finish(&mut study, true), 3 + 2 * q);
        assert_eq!(study.report(), &fixed);
        assert_eq!(study.work(), NoisyStudyWork { model_fit_attempts: 3, acquisition_batch_attempts: 2 });
        assert!(study.pending().is_empty());
        assert!(study.ask().is_empty());
        let work = study.work();
        study.advance().unwrap();
        assert_eq!(study.work(), work);
    }
}

#[test]
fn g5_mid_batch_checkpoint_retains_results_without_repeating_simulations() {
    let c = config();
    let mut study = NoisyStudy::new(91, 1, 3, 2, &c);
    let jobs = study.ask();
    study.tell(&jobs[1], observation(&jobs[1].x)).unwrap();
    assert_eq!(study.ask(), vec![jobs[0].clone(), jobs[2].clone()]);
    assert_eq!(study.advance(), Err(NoisyStudyError::WaitingForObservations));
    assert_eq!(study.work(), NoisyStudyWork::default());
    let mut checkpoint = study.clone();
    assert_eq!(finish(&mut study, false), 6);
    assert_eq!(finish(&mut checkpoint, true), 6);
    assert_eq!(study.report(), checkpoint.report());
    assert_eq!(study.report(), &minimize_noisy(&mut observation, 1, 3, 2, &c));
}

#[test]
fn g0_delivery_retries_are_idempotent_and_conflicts_never_overwrite() {
    let mut study = NoisyStudy::new(7, 1, 2, 1, &config());
    let jobs = study.ask();
    let job = &jobs[0];
    let value = observation(&job.x);
    assert_eq!(study.tell(job, value), Ok(TellDisposition::Recorded));
    assert_eq!(study.tell(job, value), Ok(TellDisposition::Replayed));
    let pending = study.pending().to_vec();
    let changed = NoisyObservation { value: value.value + 1.0, ..value };
    assert_eq!(study.tell(job, changed), Err(NoisyStudyError::ConflictingObservation));
    assert_eq!(study.pending(), pending.as_slice());
    finish(&mut study, true);
    assert_eq!(study.tell(job, value), Ok(TellDisposition::Replayed));
    assert_eq!(study.tell(job, changed), Err(NoisyStudyError::ConflictingObservation));
    assert_eq!(study.report().observations[0], value);
}

#[test]
fn g0_foreign_unknown_modified_and_invalid_responses_do_not_consume_jobs() {
    let mut study = NoisyStudy::new(3, 1, 2, 0, &config());
    let jobs = study.ask();
    let mut job = jobs[0].clone();
    let value = observation(&job.x);
    job.study_id += 1;
    assert_eq!(study.tell(&job, value), Err(NoisyStudyError::ForeignStudy));
    job = jobs[0].clone(); job.evaluation = usize::MAX;
    assert_eq!(study.tell(&job, value), Err(NoisyStudyError::UnknownEvaluation));
    job = jobs[0].clone(); job.x[0] += 0.1;
    assert_eq!(study.tell(&job, value), Err(NoisyStudyError::RequestMismatch));
    for invalid in [NoisyObservation { value: f64::NAN, ..value },
        NoisyObservation { value: f64::INFINITY, ..value },
        NoisyObservation { noise_variance: -1.0, ..value },
        NoisyObservation { noise_variance: f64::NAN, ..value }] {
        assert_eq!(study.tell(&jobs[0], invalid), Err(NoisyStudyError::InvalidObservation));
    }
    assert_eq!(study.ask(), jobs);
    assert_eq!(study.work(), NoisyStudyWork::default());
    assert_eq!(finish(&mut study, true), 2);
    assert!(study.is_complete());
}

#[test]
fn g4_cancelled_acquisition_keeps_received_results_and_records_spent_attempts() {
    let c = config();
    let mut study = NoisyStudy::new(8, 1, 3, 1, &c);
    for job in study.ask() { study.tell(&job, observation(&job.x)).unwrap(); }
    let pending = study.pending().to_vec();
    let mut resumed = study.clone();
    assert_eq!(study.advance_controlled(&mut || false), Err(NoisyStudyError::Cancelled));
    assert_eq!(study.work(), NoisyStudyWork::default());
    let mut checks = 0;
    assert_eq!(study.advance_controlled(&mut || { checks += 1; checks < 5 }),
        Err(NoisyStudyError::Cancelled));
    assert_eq!(checks, 5); // after completing the first greedy slot
    assert_eq!(study.pending(), pending.as_slice());
    assert!(study.report().x.is_empty());
    assert!(study.ask().is_empty()); // all physical results are already retained
    assert_eq!(study.work(), NoisyStudyWork { model_fit_attempts: 1, acquisition_batch_attempts: 1 });
    assert_eq!(finish(&mut study, true), 2);
    assert_eq!(finish(&mut resumed, false), 2);
    assert_eq!(study.report(), resumed.report());
    assert_eq!(study.work().model_fit_attempts, resumed.work().model_fit_attempts + 1);
    assert_eq!(study.work().acquisition_batch_attempts, resumed.work().acquisition_batch_attempts + 1);
}

#[test]
fn g0_inadmissible_completed_data_are_retained_not_rescheduled_or_noise_repaired() {
    let mut c = config();
    // Three initial points, but only two representable coordinates in this
    // interval: at least two exact noiseless constraints must coincide.
    c.bounds = (1.0, 1.0 + f64::EPSILON);
    let mut study = NoisyStudy::new(19, 1, 3, 1, &c);
    for job in study.ask() {
        study.tell(&job, NoisyObservation { value: 1.0, noise_variance: 0.0 }).unwrap();
    }
    let pending = study.pending().to_vec();
    assert_eq!(study.advance(), Err(NoisyStudyError::InvalidModel));
    assert_eq!(study.advance(), Err(NoisyStudyError::InvalidModel));
    assert_eq!(study.pending(), pending.as_slice());
    assert!(study.ask().is_empty());
    assert!(!study.is_complete());
    assert!(study.report().observations.is_empty());
    assert_eq!(study.work().model_fit_attempts, 2);
    assert_eq!(study.work().acquisition_batch_attempts, 0);
}

fn learning() -> NoisyLearningConfig {
    NoisyLearningConfig { refit_every: 2, fit: crate::hyper::HeteroFitConfig {
        lengthscale_bounds: vec![(0.05, 2.0)], signal_bounds: (0.01, 5.0),
        starts: 2, max_iterations: 12, max_evaluations: 20, gradient_tolerance: 1e-7, seed: 13,
    } }
}

fn assert_same_fits(a: &[NoisyFitRecord], b: &[NoisyFitRecord]) {
    assert_eq!(a.len(), b.len());
    for (a, b) in a.iter().zip(b) {
        assert_eq!((a.after_batches, a.observation_count, a.evaluations, a.seed),
            (b.after_batches, b.observation_count, b.evaluations, b.seed));
        assert_eq!(a.kernel.family, b.kernel.family);
        assert_eq!(a.kernel.lengthscales, b.kernel.lengthscales);
        assert_eq!(a.kernel.signal.to_bits(), b.kernel.signal.to_bits());
        assert_eq!(a.initial_lml.to_bits(), b.initial_lml.to_bits());
        assert_eq!(a.lml.to_bits(), b.lml.to_bits());
        assert_eq!(a.projected_gradient_inf.to_bits(), b.projected_gradient_inf.to_bits());
        assert_eq!(a.evaluation_limit_reached, b.evaluation_limit_reached);
    }
}

#[test]
fn learned_external_studies_match_synchronous_refits_and_all_outputs() {
    for q in [1, 2] {
        let mut c = config(); c.q = q;
        let policy = learning();
        let expected = crate::learning::minimize_noisy_with_learning(
            &mut observation, 1, 3, 3, &c, &policy).unwrap();
        let mut study = NoisyStudy::new_with_learning(45, 1, 3, 3, &c, &policy).unwrap();
        assert_eq!(finish(&mut study, true), 3 + 3 * q);
        assert_eq!(study.report(), &expected.report);
        assert_same_fits(study.fit_records(), &expected.fits);
        assert_eq!(study.fit_records().iter().map(|r| r.after_batches).collect::<Vec<_>>(), vec![0, 2]);
        assert!(study.fit_records()[0].lml > study.fit_records()[0].initial_lml);
        assert_eq!(study.kernel().signal, expected.fits.last().unwrap().kernel.signal);
    }
}

#[test]
fn learned_partial_batch_clone_preserves_the_next_warm_start() {
    let c = config();
    let mut study = NoisyStudy::new_with_learning(46, 1, 3, 2, &c, &learning()).unwrap();
    for job in study.ask() { study.tell(&job, observation(&job.x)).unwrap(); }
    study.advance().unwrap();
    let job = study.ask()[1].clone();
    study.tell(&job, observation(&job.x)).unwrap();
    let mut clone = study.clone();
    assert_eq!(finish(&mut study, true), 3);
    assert_eq!(finish(&mut clone, false), 3);
    assert_eq!(study.report(), clone.report());
    assert_same_fits(study.fit_records(), clone.fit_records());
}

#[test]
fn cancellation_after_kernel_learning_does_not_commit_the_kernel_or_fit_record() {
    let c = config();
    let policy = learning();
    let mut study = NoisyStudy::new_with_learning(47, 1, 3, 1, &c, &policy).unwrap();
    for job in study.ask() { study.tell(&job, observation(&job.x)).unwrap(); }
    let before = study.clone();
    // Measure the successful path's checkpoints, then refuse the final publish
    // check. No dependence on the number of backtracks chosen by the learner.
    let mut probe = before.clone();
    let mut count = 0;
    probe.advance_controlled(&mut || { count += 1; true }).unwrap();
    let mut checks = 0;
    assert_eq!(study.advance_controlled(&mut || { checks += 1; checks < count }),
        Err(NoisyStudyError::Cancelled));
    assert_eq!(study.pending(), before.pending());
    assert_eq!(study.report(), before.report());
    assert!(study.fit_records().is_empty());
    assert_eq!(study.kernel().signal.to_bits(), c.kernel.signal.to_bits());
    assert_eq!(study.work(), NoisyStudyWork { model_fit_attempts: 1, acquisition_batch_attempts: 1 });
    study.advance().unwrap();
    assert_eq!(study.ask(), probe.ask());
    assert_eq!(study.report(), probe.report());
    assert_same_fits(study.fit_records(), probe.fit_records());
}

#[test]
fn cancellation_inside_learner_retains_received_jobs_for_retry() {
    let mut study = NoisyStudy::new_with_learning(48, 1, 3, 0, &config(), &learning()).unwrap();
    for job in study.ask() { study.tell(&job, observation(&job.x)).unwrap(); }
    let pending = study.pending().to_vec();
    let mut checks = 0;
    assert_eq!(study.advance_controlled(&mut || { checks += 1; checks < 4 }),
        Err(NoisyStudyError::Cancelled));
    assert_eq!(study.pending(), pending);
    assert!(study.report().x.is_empty());
    assert!(study.fit_records().is_empty());
    assert_eq!(finish(&mut study, false), 0);
    assert_eq!(study.fit_records().len(), 1);
    assert_eq!(study.work().model_fit_attempts, 2);
}

#[test]
fn learned_constructor_refuses_invalid_work_policy_and_one_probe_recovers_fixed_study() {
    let c = config();
    let mut bad = learning(); bad.refit_every = 0;
    assert!(matches!(NoisyStudy::new_with_learning(49, 1, 3, 2, &c, &bad),
        Err(HeteroFitError::InvalidConfig)));
    bad = learning(); bad.fit.max_evaluations = usize::MAX;
    assert!(matches!(NoisyStudy::new_with_learning(49, 1, 3, 2, &c, &bad),
        Err(HeteroFitError::InvalidConfig)));
    let mut policy = learning(); policy.fit.max_evaluations = 1;
    let mut learned = NoisyStudy::new_with_learning(50, 1, 3, 2, &c, &policy).unwrap();
    let mut fixed = NoisyStudy::new(50, 1, 3, 2, &c);
    finish(&mut learned, true); finish(&mut fixed, false);
    assert_eq!(learned.report(), fixed.report());
    assert_eq!(learned.work(), fixed.work());
    assert!(learned.fit_records().iter().all(|fit| fit.evaluations == 1));
}

fn checkpoint_limits() -> CheckpointLimits {
    CheckpointLimits { max_bytes: 1 << 20, max_evaluations: 1000 }
}

fn restore(study: &NoisyStudy) -> NoisyStudy {
    let bytes = study.checkpoint_bytes();
    let restored = NoisyStudy::from_checkpoint(&bytes, checkpoint_limits()).unwrap();
    assert_eq!(restored.checkpoint_bytes(), bytes, "canonical exact-bit roundtrip");
    assert_eq!(restored.report(), study.report());
    assert_eq!(restored.pending(), study.pending());
    assert_eq!(restored.work(), study.work());
    assert_eq!(restored.is_complete(), study.is_complete());
    assert_same_fits(restored.fit_records(), study.fit_records());
    restored
}

#[test]
fn checkpoint_restarts_fixed_and_learned_studies_after_every_response_and_advance() {
    for learned in [false, true] {
        let c = config(); let policy = learning();
        let mut live = if learned {
            NoisyStudy::new_with_learning(51, 1, 3, 2, &c, &policy).unwrap()
        } else { NoisyStudy::new(51, 1, 3, 2, &c) };
        live = restore(&live);
        let mut calls = 0;
        while !live.is_complete() {
            for job in live.ask().into_iter().rev() {
                let response = observation(&job.x);
                live.tell(&job, response).unwrap(); calls += 1;
                live = restore(&live);
                assert_eq!(live.tell(&job, response), Ok(TellDisposition::Replayed));
            }
            live.advance().unwrap(); live = restore(&live);
        }
        assert_eq!(calls, 7);
        if learned {
            let expected = crate::learning::minimize_noisy_with_learning(
                &mut observation, 1, 3, 2, &c, &policy).unwrap();
            assert_eq!(live.report(), &expected.report);
            assert_same_fits(live.fit_records(), &expected.fits);
        } else {
            assert_eq!(live.report(), &minimize_noisy(&mut observation, 1, 3, 2, &c));
        }
    }
}

#[test]
fn checkpoint_keeps_cancelled_fit_attempts_and_does_not_reschedule_received_physics() {
    let mut live = NoisyStudy::new_with_learning(52, 1, 3, 1, &config(), &learning()).unwrap();
    for job in live.ask() { live.tell(&job, observation(&job.x)).unwrap(); }
    let mut calls = 0;
    assert_eq!(live.advance_controlled(&mut || { calls += 1; calls < 4 }), Err(NoisyStudyError::Cancelled));
    let mut resumed = restore(&live);
    assert!(resumed.ask().is_empty());
    assert_eq!(resumed.work().model_fit_attempts, 1);
    assert_eq!(finish(&mut resumed, false), 2);
    assert_eq!(finish(&mut live, true), 2);
    assert_eq!(resumed.report(), live.report());
    assert_same_fits(resumed.fit_records(), live.fit_records());
}

#[test]
fn checkpoint_retains_invalid_model_pending_results_for_inspection() {
    let mut c = config(); c.bounds = (1.0, 1.0 + f64::EPSILON);
    let mut live = NoisyStudy::new(53, 1, 3, 0, &c);
    for job in live.ask() {
        live.tell(&job, NoisyObservation { value: 1.0, noise_variance: 0.0 }).unwrap();
    }
    assert_eq!(live.advance(), Err(NoisyStudyError::InvalidModel));
    let mut resumed = restore(&live);
    assert!(resumed.ask().is_empty());
    assert_eq!(resumed.advance(), Err(NoisyStudyError::InvalidModel));
    assert_eq!(resumed.pending(), live.pending());
    assert_eq!(resumed.work().model_fit_attempts, 2);
}

#[test]
fn checkpoint_preserves_signed_zero_response_bits_and_large_study_ids() {
    let mut live = NoisyStudy::new(u64::MAX, 1, 2, 0, &config());
    let job = live.ask()[0].clone();
    let response = NoisyObservation { value: -0.0, noise_variance: -0.0 };
    live.tell(&job, response).unwrap();
    let mut resumed = restore(&live);
    assert_eq!(resumed.tell(&job, response), Ok(TellDisposition::Replayed));
    assert_eq!(resumed.tell(&job, NoisyObservation { value: 0.0, ..response }),
        Err(NoisyStudyError::ConflictingObservation));
}

fn reseal(bytes: &mut [u8]) {
    let end = bytes.len() - 8;
    let mut sum = 0xcbf29ce484222325u64;
    for byte in &bytes[..end] { sum = (sum ^ u64::from(*byte)).wrapping_mul(0x100000001b3); }
    bytes[end..].copy_from_slice(&sum.to_le_bytes());
}

#[test]
fn checkpoint_refuses_truncation_corruption_versions_and_caller_limits() {
    let frame = NoisyStudy::new(54, 1, 3, 1, &config()).checkpoint_bytes();
    for end in 0..frame.len() {
        assert!(NoisyStudy::from_checkpoint(&frame[..end], checkpoint_limits()).is_err());
    }
    let mut bad = frame.clone(); bad[6] = 2;
    assert!(matches!(NoisyStudy::from_checkpoint(&bad, checkpoint_limits()), Err(CheckpointError::IncompatibleVersion)));
    bad = frame.clone(); bad[40] ^= 1;
    assert!(matches!(NoisyStudy::from_checkpoint(&bad, checkpoint_limits()), Err(CheckpointError::ChecksumMismatch)));
    let mut limits = checkpoint_limits(); limits.max_bytes = frame.len() - 1;
    assert!(matches!(NoisyStudy::from_checkpoint(&frame, limits), Err(CheckpointError::LimitExceeded)));
    limits = checkpoint_limits(); limits.max_evaluations = 4; // full plan is 3 + 1*2
    assert!(matches!(NoisyStudy::from_checkpoint(&frame, limits), Err(CheckpointError::LimitExceeded)));
}

#[test]
fn checkpoint_structural_admission_is_not_only_a_checksum_check() {
    let frame = NoisyStudy::new(55, 1, 3, 1, &config()).checkpoint_bytes();
    // v1, dimension one, fixed kernel: these offsets are frozen grammar fields,
    // not random byte fuzz. Recompute checksum to reach structural admission.
    for (offset, word) in [(16, 0), (24, 0), (64, f64::NAN.to_bits()),
        (88, 0), (128, 2), (152, 4), (160, u64::MAX), (168, 0)] {
        let mut bad = frame.clone(); bad[offset..offset + 8].copy_from_slice(&word.to_le_bytes());
        reseal(&mut bad);
        assert!(NoisyStudy::from_checkpoint(&bad, checkpoint_limits()).is_err(), "offset {offset}");
    }
    let mut trailing = frame[..frame.len() - 8].to_vec();
    trailing.extend_from_slice(&[0; 16]); reseal(&mut trailing);
    assert!(matches!(NoisyStudy::from_checkpoint(&trailing, checkpoint_limits()), Err(CheckpointError::InvalidState)));
    let mut learned = NoisyStudy::new_with_learning(56, 1, 3, 0, &config(), &learning()).unwrap();
    finish(&mut learned, false);
    // Private test-only mutation simulates structurally well-framed stored data
    // with a wrong refit cadence, zero allowance, or non-finite recommendation.
    let mut wrong = learned.clone(); wrong.fits[0].after_batches = 1;
    assert!(NoisyStudy::from_checkpoint(&wrong.checkpoint_bytes(), checkpoint_limits()).is_err());
    wrong = learned.clone(); wrong.fits[0].evaluations = 0;
    assert!(NoisyStudy::from_checkpoint(&wrong.checkpoint_bytes(), checkpoint_limits()).is_err());
    wrong = learned; wrong.report.incumbent_trace[0].mean = f64::INFINITY;
    assert!(NoisyStudy::from_checkpoint(&wrong.checkpoint_bytes(), checkpoint_limits()).is_err());
}
