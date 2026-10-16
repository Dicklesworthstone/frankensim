use super::*;
use crate::{Kernel, Matern};

fn config() -> NoisyBoConfig {
    NoisyBoConfig { bounds: (0.0, 1.0),
        kernel: Kernel { family: Matern::FiveHalves, signal: 1.0, lengthscales: vec![0.3] },
        prior_mean: 0.5, q: 1, mc_samples: 64, acq_starts: 2, acq_evals: 32, seed: 97 }
}
fn receive(s: &mut AsyncNoisyStudy, r: &NoisyRequest) {
    s.tell(r, NoisyObservation { value: (r.x[0] - 0.37).powi(2),
        noise_variance: 0.01 + 0.03 * r.x[0] }).unwrap();
}

#[test]
fn cached_refills_replay_keep_stragglers_and_use_the_direct_cached_proposal() {
    let c = config();
    let mut a = AsyncNoisyStudy::new(7, 1, 2, 6, 2, &c).unwrap();
    let initial = a.ask(); receive(&mut a, &initial[1]);
    let mut b = a.clone();
    let observed = a.jobs()[1].observation.unwrap();
    let x = vec![initial[1].x.clone()];
    let gp = Gp::try_fit_diag(&x, &[observed.value - c.prior_mean], c.kernel.clone(),
        &[observed.noise_variance]).unwrap();
    let expected = c.propose_pending_cached(&gp, &x, &[initial[0].x.clone()], 1, 0, 1e-6, &mut || true).unwrap();
    assert_eq!(a.refill_cached(1e-6), Ok(1));
    assert_eq!(b.refill_cached_controlled(1e-6, &mut || true), Ok(1));
    assert_eq!(a.jobs()[2].request.x, expected[0]);
    for id in 2..6 {
        let next = a.ask().into_iter().find(|r| r.evaluation == id).unwrap();
        receive(&mut a, &next); receive(&mut b, &next);
        assert_eq!(a.refill_cached(1e-6).unwrap(), usize::from(id < 5));
        b.refill_cached(1e-6).unwrap();
        assert_eq!(a.jobs(), b.jobs()); assert_eq!(a.model_report(), b.model_report());
        assert!(a.in_flight() <= 2 && a.jobs().len() <= 6);
        assert!(a.jobs()[0].observation.is_none(), "slow job remains reserved");
    }
    receive(&mut a, &initial[0]); receive(&mut b, &initial[0]);
    a.refill_cached(1e-6).unwrap(); b.refill_cached(1e-6).unwrap();
    assert!(a.is_complete() && b.is_complete());
    assert_eq!(a.model_report(), b.model_report());
    assert_eq!(a.model_report().unwrap().observed_evaluations.len(), 6);
}

#[test]
fn cached_preparation_failure_keeps_results_and_can_retry_with_explicit_jitter() {
    let c = config();
    let mut s = AsyncNoisyStudy::new(11, 1, 1, 2, 1, &c).unwrap();
    let r = s.ask()[0].clone();
    s.tell(&r, NoisyObservation { value: 0.0, noise_variance: 0.0 }).unwrap();
    let before = s.jobs().to_vec();
    assert_eq!(s.refill_cached(0.0), Err(AsyncStudyError::CachedPosterior(PreparedError::NotPositiveDefinite)));
    assert_eq!(s.jobs(), before.as_slice()); assert!(s.model_report().is_none());
    assert_eq!(s.work().acquisition_batch_attempts, 1);
    assert_eq!(s.refill_cached(1e-6), Ok(1));
    assert_eq!(s.jobs()[0].observation.unwrap().noise_variance, 0.0);
    assert_eq!(s.work().acquisition_batch_attempts, 2);
}

#[test]
fn cancelled_cached_work_does_not_publish_a_model_or_jobs() {
    let c = config();
    let mut s = AsyncNoisyStudy::new(13, 1, 2, 4, 2, &c).unwrap();
    let first = s.ask()[1].clone(); receive(&mut s, &first);
    let original = s.jobs().to_vec(); let mut clean = s.clone();
    let mut calls = 0;
    assert_eq!(s.refill_cached_controlled(1e-6, &mut || { calls += 1; calls < 100 }),
        Err(AsyncStudyError::Cancelled));
    assert_eq!(calls, 100); assert_eq!(s.jobs(), original.as_slice());
    assert!(s.model_report().is_none());
    s.refill_cached(1e-6).unwrap(); clean.refill_cached(1e-6).unwrap();
    assert_eq!(s.jobs(), clean.jobs()); assert_eq!(s.model_report(), clean.model_report());
    assert_eq!(s.work().model_fit_attempts, clean.work().model_fit_attempts + 1);
}

#[test]
fn invalid_jitter_does_not_release_initial_jobs_or_spend_work() {
    let mut s = AsyncNoisyStudy::new(17, 1, 3, 5, 1, &config()).unwrap();
    let r = s.ask()[0].clone(); receive(&mut s, &r);
    let before = s.jobs().to_vec();
    for invalid in [f64::NAN, f64::INFINITY, -1.0] {
        assert_eq!(s.refill_cached(invalid), Err(AsyncStudyError::InvalidPolicy));
        assert_eq!(s.jobs(), before.as_slice()); assert_eq!(s.work(), NoisyStudyWork::default());
    }
    assert_eq!(s.refill_cached(1e-6), Ok(1));
    assert_eq!(s.work(), NoisyStudyWork::default()); // initial jobs need no model
}
