use super::*;

fn limits() -> CheckpointLimits {
    CheckpointLimits { max_bytes: 1 << 20, max_evaluations: 64 }
}
fn config(learned: bool) -> ConstrainedBoConfig {
    let k = Kernel { family: Matern::FiveHalves, signal: 1.0, lengthscales: vec![0.3] };
    let learning = |cadence| NoisyLearningConfig { refit_every: cadence, fit: HeteroFitConfig {
        lengthscale_bounds: vec![(0.1, 1.0)], signal_bounds: (0.1, 3.0),
        starts: 1, max_iterations: 3, max_evaluations: 5, gradient_tolerance: 1e-7, seed: 31,
    }};
    ConstrainedBoConfig {
        search: NoisyBoConfig { bounds: (0.0, 1.0), kernel: k.clone(), prior_mean: 0.5,
            q: 1, mc_samples: 32, acq_starts: 2, acq_evals: 48, seed: 97 },
        objective_learning: learned.then(|| learning(2)),
        constraints: vec![
            OutcomeConstraint { kernel: k.clone(), prior_mean: 0.0, upper_bound: 0.8,
                learning: learned.then(|| learning(3)) },
            OutcomeConstraint { kernel: k, prior_mean: 0.0, upper_bound: 1.0, learning: None },
        ],
        reference: 4.0, recommendation_probability: 0.5,
    }
}
fn observation(x: &[f64]) -> ConstrainedObservation {
    ConstrainedObservation {
        objective: NoisyObservation { value: (x[0] - 0.37).powi(2), noise_variance: 0.02 },
        constraints: vec![
            NoisyObservation { value: x[0], noise_variance: 0.01 },
            NoisyObservation { value: 0.3 - x[0], noise_variance: 0.03 },
        ],
    }
}
fn restore(s: &AsyncConstrainedStudy) -> AsyncConstrainedStudy {
    let bytes = s.checkpoint_bytes();
    let r = AsyncConstrainedStudy::from_checkpoint(&bytes, limits()).unwrap();
    assert_eq!(r.checkpoint_bytes(), bytes, "exact bits, kernels, counters and model membership");
    assert_eq!(r.ask(), s.ask());
    assert_eq!(r.model_report(), s.model_report());
    assert_eq!(r.is_complete(), s.is_complete());
    r
}

#[test]
fn restart_after_every_arrival_and_refill_preserves_fixed_and_learned_campaigns() {
    for learned in [false, true] {
        let mut live = AsyncConstrainedStudy::new_with_learning(42, 1, 3, 6, 3, &config(learned)).unwrap();
        let mut resumed = restore(&live);
        let mut evaluations = Vec::new();
        while !live.is_complete() {
            // Keep the oldest jobs running while replacements finish first.
            let request = live.ask().pop().unwrap();
            evaluations.push(request.evaluation);
            let o = observation(&request.x);
            live.tell(&request, o.clone()).unwrap(); resumed.tell(&request, o.clone()).unwrap();
            resumed = restore(&resumed);
            assert_eq!(resumed.tell(&request, o).unwrap(), TellDisposition::Replayed);
            assert_eq!(live.refill().unwrap(), resumed.refill().unwrap());
            resumed = restore(&resumed);
            assert_eq!(live.checkpoint_bytes(), resumed.checkpoint_bytes());
        }
        evaluations.sort_unstable();
        assert_eq!(evaluations, (0..6).collect::<Vec<_>>());
        let work = resumed.work();
        assert_eq!(resumed.refill().unwrap(), 0);
        assert_eq!(resumed.work(), work);
    }
}

#[test]
fn initial_waves_and_late_unmodeled_results_do_not_become_new_jobs_or_modeled_rows() {
    let c = config(false);
    let mut wave = AsyncConstrainedStudy::new(10, 1, 5, 5, 2, &c).unwrap();
    while !wave.is_complete() {
        wave = restore(&wave);
        for job in wave.ask() { wave.tell(&job, observation(&job.x)).unwrap(); }
        wave = restore(&wave); wave.refill().unwrap();
        assert!(wave.in_flight() <= 2);
    }
    let mut s = AsyncConstrainedStudy::new(11, 1, 3, 3, 3, &c).unwrap();
    let jobs = s.ask();
    for j in &jobs[1..] { s.tell(j, observation(&j.x)).unwrap(); }
    s.refill().unwrap();
    s.tell(&jobs[0], observation(&jobs[0].x)).unwrap();
    let mut r = restore(&s);
    assert!(r.ask().is_empty() && !r.is_complete());
    assert_eq!(r.model_report().unwrap().observed_evaluations, vec![1, 2]);
    assert_eq!(r.refill().unwrap(), 0);
    assert!(r.is_complete());
    assert_eq!(r.model_report().unwrap().observed_evaluations, vec![0, 1, 2]);
}

#[test]
fn cancelled_learning_and_failed_models_keep_results_and_attempts_across_restart() {
    let mut s = AsyncConstrainedStudy::new_with_learning(12, 1, 2, 2, 2, &config(true)).unwrap();
    let j = s.ask()[1].clone(); s.tell(&j, observation(&j.x)).unwrap();
    let before = s.jobs().to_vec(); let mut checks = 0;
    assert!(s.refill_controlled(&mut || { checks += 1; checks < 4 }).is_err());
    assert_eq!(checks, 4);
    assert_eq!(s.jobs(), before.as_slice());
    let mut r = restore(&s);
    assert_eq!(r.work().model_fit_attempts, 1);
    assert_eq!(s.refill().unwrap(), r.refill().unwrap());
    assert_eq!(s.checkpoint_bytes(), r.checkpoint_bytes());

    let mut c = config(false); c.search.bounds = (1.0, 1.0 + f64::EPSILON);
    let mut bad = AsyncConstrainedStudy::new(13, 1, 3, 3, 3, &c).unwrap();
    for j in bad.ask() {
        let mut o = observation(&j.x); o.objective.noise_variance = 0.0;
        bad.tell(&j, o).unwrap();
    }
    assert!(bad.refill().is_err());
    let mut bad = restore(&bad);
    assert!(bad.ask().is_empty() && !bad.is_complete());
    assert!(bad.refill().is_err());
    assert_eq!(bad.work().model_fit_attempts, 2);
}

fn reseal(bytes: &mut [u8]) {
    let end = bytes.len() - 8;
    let h = checksum(&bytes[..end]); bytes[end..].copy_from_slice(&h.to_le_bytes());
}
#[test]
fn bounded_decoder_rejects_truncation_unknown_versions_counts_and_trailing_words() {
    let s = AsyncConstrainedStudy::new_with_learning(14, 1, 3, 9, 3, &config(true)).unwrap();
    let bytes = s.checkpoint_bytes();
    for end in 0..bytes.len() { assert!(AsyncConstrainedStudy::from_checkpoint(&bytes[..end], limits()).is_err()); }
    let mut bad = bytes.clone(); bad[6] ^= 1;
    assert!(matches!(AsyncConstrainedStudy::from_checkpoint(&bad, limits()), Err(CheckpointError::IncompatibleVersion)));
    let mut bad = bytes.clone(); bad[60] ^= 1;
    assert!(matches!(AsyncConstrainedStudy::from_checkpoint(&bad, limits()), Err(CheckpointError::ChecksumMismatch)));
    let mut small = limits(); small.max_bytes = bytes.len() - 1;
    assert!(matches!(AsyncConstrainedStudy::from_checkpoint(&bytes, small), Err(CheckpointError::LimitExceeded)));
    small = limits(); small.max_evaluations = 8; // total budget, not just three issued jobs
    assert!(matches!(AsyncConstrainedStudy::from_checkpoint(&bytes, small), Err(CheckpointError::LimitExceeded)));
    let mut bad = bytes.clone(); bad[40..48].copy_from_slice(&u64::MAX.to_le_bytes()); reseal(&mut bad);
    assert!(AsyncConstrainedStudy::from_checkpoint(&bad, limits()).is_err());
    let mut bad = bytes[..bytes.len() - 8].to_vec(); bad.extend_from_slice(&[0; 16]); reseal(&mut bad);
    assert!(matches!(AsyncConstrainedStudy::from_checkpoint(&bad, limits()), Err(CheckpointError::InvalidState)));
}

#[test]
fn structural_admission_rejects_inconsistent_models_even_with_valid_checksums() {
    let mut s = AsyncConstrainedStudy::new(15, 1, 3, 3, 3, &config(false)).unwrap();
    for j in s.ask() { s.tell(&j, observation(&j.x)).unwrap(); }
    s.refill().unwrap();
    assert!(s.model_report().unwrap().recommendation.is_some());
    for mutation in 0..10 {
        let mut b = s.clone();
        match mutation {
            0 => b.model_report.as_mut().unwrap().observed_evaluations = vec![1, 1],
            1 => b.model_report.as_mut().unwrap().observed_evaluations[0] = 99,
            2 => b.model_updates = 4,
            3 => b.jobs[0].observation = None,
            4 => b.kernels[0].signal *= 2.0, // fixed covariance may not change
            5 => b.initial[0][0] = 2.0,
            6 => b.max_in_flight = 0,
            7 => b.likelihood_evaluations = 1, // no output learns
            8 => b.config.search.q = 0,
            _ => b.model_report.as_mut().unwrap().recommendation.as_mut().unwrap().observation_index = 3,
        }
        assert!(AsyncConstrainedStudy::from_checkpoint(&b.checkpoint_bytes(), limits()).is_err(), "mutation {mutation}");
    }
}

#[test]
fn signed_zero_and_complete_constraint_vectors_roundtrip_without_refitting() {
    let mut c = config(false); c.search.prior_mean = -0.0;
    let mut s = AsyncConstrainedStudy::new(16, 1, 2, 2, 2, &c).unwrap();
    let j = s.ask()[0].clone();
    let zero = NoisyObservation { value: -0.0, noise_variance: -0.0 };
    s.tell(&j, ConstrainedObservation { objective: zero, constraints: vec![zero; 2] }).unwrap();
    let r = restore(&s);
    let received = r.jobs()[0].observation.as_ref().unwrap();
    assert_eq!(received.objective.value.to_bits(), (-0.0f64).to_bits());
    assert_eq!(received.constraints[1].noise_variance.to_bits(), (-0.0f64).to_bits());
    assert_eq!(r.work(), NoisyStudyWork::default());
}
