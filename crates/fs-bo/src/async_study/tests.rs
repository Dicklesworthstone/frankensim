use super::*;
use crate::gp::{Kernel, Matern};
use crate::noisy::minimize_noisy;

fn config() -> NoisyBoConfig {
    NoisyBoConfig { bounds: (0.0, 1.0),
        kernel: Kernel { family: Matern::FiveHalves, signal: 1.0, lengthscales: vec![0.3] },
        prior_mean: 0.5, q: 2, mc_samples: 32, acq_starts: 2, acq_evals: 32, seed: 97 }
}
fn observation(x: &[f64]) -> NoisyObservation {
    NoisyObservation { value: (x[0] - 0.37).powi(2), noise_variance: 0.01 + 0.03 * x[0] }
}
fn receive(study: &mut AsyncNoisyStudy, request: &NoisyRequest) {
    assert_eq!(study.tell(request, observation(&request.x)), Ok(TellDisposition::Recorded));
}

#[test]
fn replacement_uses_received_data_and_all_pending_jobs_without_a_batch_barrier() {
    let c = config();
    let mut s = AsyncNoisyStudy::new(41, 1, 3, 7, 3, &c).unwrap();
    let initial = s.ask();
    receive(&mut s, &initial[2]);
    let x = vec![initial[2].x.clone()];
    let o = observation(&x[0]);
    let gp = Gp::try_fit_diag(&x, &[o.value - c.prior_mean], c.kernel.clone(), &[o.noise_variance]).unwrap();
    let pending = vec![initial[0].x.clone(), initial[1].x.clone()];
    let expected = c.propose_pending(&gp, &x, &pending, 1, 0, &mut || true).unwrap();
    assert_eq!(s.refill(), Ok(1));
    assert_eq!(s.jobs().len(), 4);
    assert_eq!(s.in_flight(), 3);
    assert_eq!(s.jobs()[3].request.x, expected[0]);
    assert!(s.jobs()[0].observation.is_none() && s.jobs()[1].observation.is_none());
    let r = s.model_report().unwrap();
    assert_eq!(r.observed_evaluations, vec![2]);
    assert_eq!(r.recommended_evaluation, 2);
    let (mean, variance) = gp.predict(&x[0]);
    assert_eq!((r.mean, r.variance), (mean + c.prior_mean, variance));
}

#[test]
fn straggler_does_not_prevent_exhausting_the_other_workers_job_budget() {
    let mut c = config(); c.q = 1;
    let mut s = AsyncNoisyStudy::new(7, 1, 2, 6, 2, &c).unwrap();
    let slow = s.ask()[0].clone();
    for next in 1..6 {
        let fast = s.ask().into_iter().find(|r| r.evaluation != slow.evaluation).unwrap();
        assert_eq!(fast.evaluation, next);
        receive(&mut s, &fast);
        assert_eq!(s.refill().unwrap(), usize::from(next < 5));
        assert!(s.jobs()[0].observation.is_none());
        assert!(s.jobs().len() <= 6 && s.in_flight() <= 2);
    }
    assert_eq!(s.ask(), vec![slow.clone()]);
    assert!(!s.is_complete());
    receive(&mut s, &slow);
    assert!(!s.is_complete(), "the final posterior must still be updated");
    assert_eq!(s.refill(), Ok(0));
    assert!(s.is_complete());
    assert_eq!(s.model_report().unwrap().observed_evaluations, (0..6).collect::<Vec<_>>());
    let work = s.work();
    assert_eq!(s.refill(), Ok(0));
    assert_eq!(s.work(), work);
}

#[test]
fn initial_design_and_nonmultiple_final_refill_obey_both_caps() {
    let mut c = config(); c.q = 4;
    let mut s = AsyncNoisyStudy::new(9, 1, 5, 8, 2, &c).unwrap();
    assert_eq!(s.jobs().len(), 2);
    let initial = c.initial_design(1, 5);
    let mut physical_calls = 0;
    while !s.is_complete() {
        for r in s.ask().into_iter().rev() {
            receive(&mut s, &r);
            physical_calls += 1;
        }
        s.refill().unwrap();
        assert!(s.in_flight() <= 2 && s.jobs().len() <= 8);
    }
    assert_eq!(physical_calls, 8);
    for (job, expected) in s.jobs()[..5].iter().zip(initial) { assert_eq!(job.request.x, expected); }
}

#[test]
fn equal_availability_sets_replay_independently_of_delivery_order() {
    let c = config();
    let mut a = AsyncNoisyStudy::new(13, 1, 3, 7, 3, &c).unwrap();
    let jobs = a.ask();
    receive(&mut a, &jobs[2]);
    let mut b = a.clone(); // partial in-memory checkpoint
    receive(&mut a, &jobs[1]);
    receive(&mut b, &jobs[1]);
    a.refill().unwrap(); b.refill().unwrap();
    assert_eq!(a.jobs(), b.jobs());
    assert_eq!(a.model_report(), b.model_report());
    assert_eq!(a.model_report().unwrap().observed_evaluations, vec![1, 2]);
    let mut other = AsyncNoisyStudy::new(13, 1, 3, 7, 3, &c).unwrap();
    receive(&mut other, &jobs[1]); receive(&mut other, &jobs[2]); other.refill().unwrap();
    assert_eq!(a.jobs(), other.jobs());
    assert_eq!(a.model_report(), other.model_report());
}

#[test]
fn cancellation_discards_partial_proposals_not_received_results_or_attempts() {
    let mut s = AsyncNoisyStudy::new(19, 1, 3, 7, 3, &config()).unwrap();
    let initial = s.ask();
    receive(&mut s, &initial[1]); receive(&mut s, &initial[2]);
    let mut clean = s.clone();
    let before = s.jobs().to_vec();
    assert_eq!(s.refill_controlled(&mut || false), Err(AsyncStudyError::Cancelled));
    assert_eq!(s.work(), NoisyStudyWork::default());
    let mut checks = 0;
    assert_eq!(s.refill_controlled(&mut || { checks += 1; checks < 6 }), Err(AsyncStudyError::Cancelled));
    assert_eq!(checks, 6); // second slot boundary, after the first search
    assert_eq!(s.jobs(), before.as_slice());
    assert!(s.model_report().is_none());
    assert_eq!(s.work(), NoisyStudyWork { model_fit_attempts: 1, acquisition_batch_attempts: 1 });
    s.refill().unwrap(); clean.refill().unwrap();
    assert_eq!(s.jobs(), clean.jobs());
    assert_eq!(s.model_report(), clean.model_report());
    assert_eq!(s.work().acquisition_batch_attempts, clean.work().acquisition_batch_attempts + 1);
}

#[test]
fn response_retries_do_not_overwrite_data_or_release_an_extra_slot() {
    let mut s = AsyncNoisyStudy::new(23, 1, 2, 2, 2, &config()).unwrap();
    let r = s.ask()[1].clone(); let o = observation(&r.x);
    receive(&mut s, &r);
    let before = s.jobs().to_vec();
    assert_eq!(s.tell(&r, o), Ok(TellDisposition::Replayed));
    assert_eq!(s.tell(&r, NoisyObservation { value: o.value + 1.0, ..o }), Err(NoisyStudyError::ConflictingObservation));
    let mut bad = r.clone(); bad.study_id += 1;
    assert_eq!(s.tell(&bad, o), Err(NoisyStudyError::ForeignStudy));
    bad = r.clone(); bad.evaluation = usize::MAX;
    assert_eq!(s.tell(&bad, o), Err(NoisyStudyError::UnknownEvaluation));
    bad = r.clone(); bad.x[0] += 0.1;
    assert_eq!(s.tell(&bad, o), Err(NoisyStudyError::RequestMismatch));
    assert_eq!(s.tell(&r, NoisyObservation { noise_variance: -1.0, ..o }), Err(NoisyStudyError::InvalidObservation));
    assert_eq!(s.jobs(), before.as_slice());
    assert_eq!(s.in_flight(), 1);
    let outstanding = s.ask()[0].clone(); receive(&mut s, &outstanding); s.refill().unwrap();
    assert!(s.is_complete());
    assert_eq!(s.tell(&r, o), Ok(TellDisposition::Replayed));
}

#[test]
fn final_model_refusal_retains_responses_without_minting_measurement_noise() {
    let mut c = config(); c.bounds = (1.0, 1.0 + f64::EPSILON);
    let mut s = AsyncNoisyStudy::new(31, 1, 3, 3, 3, &c).unwrap();
    for r in s.ask() { s.tell(&r, NoisyObservation { value: 1.0, noise_variance: 0.0 }).unwrap(); }
    let before = s.jobs().to_vec();
    assert_eq!(s.refill(), Err(AsyncStudyError::InvalidModel));
    assert_eq!(s.refill(), Err(AsyncStudyError::InvalidModel));
    assert_eq!(s.jobs(), before.as_slice());
    assert!(s.ask().is_empty() && !s.is_complete());
    assert!(s.model_report().is_none());
    assert_eq!(s.work(), NoisyStudyWork { model_fit_attempts: 2, acquisition_batch_attempts: 0 });
}

#[test]
fn one_worker_q1_recovers_the_synchronous_noisy_trajectory() {
    let mut c = config(); c.q = 1;
    let reference = minimize_noisy(&mut observation, 1, 3, 2, &c);
    let mut s = AsyncNoisyStudy::new(37, 1, 3, 5, 1, &c).unwrap();
    while !s.is_complete() {
        let r = s.ask()[0].clone(); receive(&mut s, &r); s.refill().unwrap();
    }
    assert_eq!(s.jobs().iter().map(|j| j.request.x.clone()).collect::<Vec<_>>(), reference.x);
    assert_eq!(s.jobs().iter().map(|j| j.observation.unwrap()).collect::<Vec<_>>(), reference.observations);
    let expected = reference.incumbent_trace.last().unwrap();
    let actual = s.model_report().unwrap();
    assert_eq!(actual.recommended_evaluation, expected.observation_index);
    assert_eq!((actual.mean, actual.variance), (expected.mean, expected.variance));
}

#[test]
fn empty_model_and_policy_admission_do_not_issue_unbudgeted_work() {
    let c = config();
    for (n_init, budget, cap) in [(0, 4, 2), (3, 2, 2), (2, 4, 0), (2, usize::MAX, 3)] {
        assert!(matches!(AsyncNoisyStudy::new(1, 1, n_init, budget, cap, &c), Err(AsyncStudyError::InvalidPolicy)));
    }
    let mut s = AsyncNoisyStudy::new(1, 1, 1, 4, 3, &c).unwrap();
    let before = s.jobs().to_vec();
    assert_eq!(s.refill(), Err(AsyncStudyError::WaitingForObservations));
    assert_eq!(s.jobs(), before.as_slice());
    assert_eq!(s.work(), NoisyStudyWork::default());
}
