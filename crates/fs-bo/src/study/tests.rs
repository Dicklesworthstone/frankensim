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
