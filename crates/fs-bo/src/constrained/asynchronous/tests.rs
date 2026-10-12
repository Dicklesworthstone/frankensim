use super::*;
use crate::{Kernel, Matern};
use crate::constrained::{OutcomeConstraint, minimize_constrained};
use crate::noisy::{NoisyBoConfig, NoisyObservation};
use crate::hyper::HeteroFitError;

fn config() -> ConstrainedBoConfig {
    let k = Kernel { family: Matern::FiveHalves, signal: 1.0, lengthscales: vec![0.3] };
    ConstrainedBoConfig {
        search: NoisyBoConfig { bounds: (0.0, 1.0), kernel: k.clone(), prior_mean: 0.0,
            q: 2, mc_samples: 64, acq_starts: 2, acq_evals: 32, seed: 91 },
        objective_learning: None,
        constraints: vec![OutcomeConstraint { kernel: Kernel { signal: 0.1, ..k },
            prior_mean: -1.0, upper_bound: 0.0, learning: None }],
        reference: 2.0, recommendation_probability: 0.9,
    }
}
fn obs(value: f64, noise_variance: f64) -> NoisyObservation { NoisyObservation { value, noise_variance } }
fn evaluate(x: &[f64]) -> ConstrainedObservation {
    ConstrainedObservation { objective: obs((x[0] - 0.37).powi(2), 0.001),
        constraints: vec![obs(-1.0 + 0.1 * x[0], 0.001 + 0.003 * x[0])] }
}
fn receive(s: &mut AsyncConstrainedStudy, request: &NoisyRequest) {
    assert_eq!(s.tell(request, evaluate(&request.x)), Ok(TellDisposition::Recorded));
}

#[test]
fn g3_one_worker_matches_synchronous_constrained_objectives_and_models() {
    let mut c = config(); c.search.q = 1;
    let reference = minimize_constrained(&mut evaluate, 1, 3, 2, &c).unwrap();
    let mut s = AsyncConstrainedStudy::new(7, 1, 3, 5, 1, &c).unwrap();
    while !s.is_complete() {
        let r = s.ask()[0].clone(); receive(&mut s, &r); s.refill().unwrap();
    }
    assert_eq!(s.jobs().iter().map(|j| j.request.x.clone()).collect::<Vec<_>>(), reference.x);
    assert_eq!(s.jobs().iter().map(|j| j.observation.clone().unwrap()).collect::<Vec<_>>(), reference.observations);
    assert_eq!(&s.model_report().unwrap().recommendation, reference.recommendations.last().unwrap());
}

#[test]
fn g3_straggler_remains_pending_while_feasible_replacements_exhaust_the_budget() {
    let mut c = config(); c.search.q = 1;
    let mut s = AsyncConstrainedStudy::new(9, 1, 2, 6, 2, &c).unwrap();
    let slow = s.ask()[0].clone();
    for id in 1..6 {
        let fast = s.ask().into_iter().find(|r| r.evaluation != slow.evaluation).unwrap();
        assert_eq!(fast.evaluation, id);
        receive(&mut s, &fast);
        assert_eq!(s.refill().unwrap(), usize::from(id < 5));
        assert!(s.jobs()[0].observation.is_none());
        assert!(s.in_flight() <= 2 && s.jobs().len() <= 6);
        assert_eq!(s.model_report().unwrap().observed_evaluations, (1..=id).collect::<Vec<_>>());
    }
    receive(&mut s, &slow);
    assert!(!s.is_complete());
    assert_eq!(s.refill(), Ok(0));
    assert!(s.is_complete());
    assert_eq!(s.model_report().unwrap().observed_evaluations, (0..6).collect::<Vec<_>>());
    let work = s.work(); assert_eq!(s.refill(), Ok(0)); assert_eq!(s.work(), work);
}

#[test]
fn g3_late_constraint_results_can_remove_a_recommendation_without_losing_history() {
    let mut c = config(); c.search.kernel.lengthscales = vec![10.0];
    c.constraints[0].kernel.lengthscales = vec![10.0];
    let mut s = AsyncConstrainedStudy::new(11, 1, 2, 2, 2, &c).unwrap();
    let jobs = s.ask();
    s.tell(&jobs[1], ConstrainedObservation { objective: obs(0.0, 0.001),
        constraints: vec![obs(-1.0, 0.1)] }).unwrap();
    assert_eq!(s.refill(), Ok(0));
    assert!(s.model_report().unwrap().recommendation.is_some());
    s.tell(&jobs[0], ConstrainedObservation { objective: obs(0.5, 0.001),
        constraints: vec![obs(10.0, 1e-12)] }).unwrap();
    s.refill().unwrap();
    assert!(s.is_complete());
    assert!(s.model_report().unwrap().recommendation.is_none());
    assert_eq!(s.jobs()[1].observation.as_ref().unwrap().constraints[0], obs(-1.0, 0.1));
}

#[test]
fn g0_two_constraints_reject_either_violation_and_map_model_rows_to_job_ids() {
    let mut c = config(); c.search.kernel.lengthscales = vec![1e-8];
    c.constraints[0].kernel.lengthscales = vec![1e-8];
    c.constraints.push(c.constraints[0].clone());
    let mut s = AsyncConstrainedStudy::new(13, 1, 4, 4, 4, &c).unwrap();
    let r = s.ask();
    for (id, f, g, h) in [(3, -100.0, 10.0, -1.0), (1, -50.0, -1.0, 10.0), (2, 1.0, -1.0, -1.0)] {
        s.tell(&r[id], ConstrainedObservation { objective: obs(f, 1e-6),
            constraints: vec![obs(g, 1e-6), obs(h, 1e-6)] }).unwrap();
    }
    s.refill().unwrap();
    let report = s.model_report().unwrap();
    assert_eq!(report.observed_evaluations, vec![1, 2, 3]);
    let best = report.recommendation.as_ref().unwrap();
    assert_eq!(best.observation_index, 1);
    assert_eq!(report.observed_evaluations[best.observation_index], 2);
    assert!(best.constraint_means.iter().all(|v| *v < 0.0));
}

#[test]
fn g0_incomplete_or_invalid_constraint_vectors_never_consume_jobs() {
    let mut s = AsyncConstrainedStudy::new(17, 1, 2, 4, 2, &config()).unwrap();
    let job = s.ask()[1].clone();
    let original = s.jobs().to_vec();
    let mut bad = evaluate(&job.x); bad.constraints.clear();
    assert!(matches!(s.tell(&job, bad), Err(AsyncConstrainedError::Model(ConstrainedBoError::OutputCount { .. }))));
    for value in [obs(f64::NAN, 1.0), obs(1.0, -1.0), obs(1.0, f64::INFINITY)] {
        let mut bad = evaluate(&job.x); bad.constraints[0] = value;
        assert_eq!(s.tell(&job, bad), Err(AsyncConstrainedError::Model(
            ConstrainedBoError::Observation { evaluation: 1, output: 1 })));
    }
    let mut foreign = job.clone(); foreign.study_id += 1;
    assert_eq!(s.tell(&foreign, evaluate(&job.x)), Err(AsyncConstrainedError::Delivery(NoisyStudyError::ForeignStudy)));
    let mut unknown = job.clone(); unknown.evaluation = usize::MAX;
    assert_eq!(s.tell(&unknown, evaluate(&job.x)), Err(AsyncConstrainedError::Delivery(NoisyStudyError::UnknownEvaluation)));
    let mut changed = job.clone(); changed.x[0] += 0.1;
    assert_eq!(s.tell(&changed, evaluate(&job.x)), Err(AsyncConstrainedError::Delivery(NoisyStudyError::RequestMismatch)));
    assert_eq!(s.jobs(), original.as_slice());
    receive(&mut s, &job);
    assert_eq!(s.tell(&job, evaluate(&job.x)), Ok(TellDisposition::Replayed));
    let mut changed = evaluate(&job.x); changed.constraints[0].noise_variance *= 2.0;
    assert_eq!(s.tell(&job, changed), Err(AsyncConstrainedError::Delivery(NoisyStudyError::ConflictingObservation)));
    assert_eq!(s.in_flight(), 1);
}

#[test]
fn g4_pause_between_output_fits_or_greedy_slots_preserves_results_and_attempts() {
    let mut s = AsyncConstrainedStudy::new(19, 1, 3, 7, 3, &config()).unwrap();
    for r in s.ask().into_iter().skip(1) { receive(&mut s, &r); }
    let mut clean = s.clone();
    let received = s.jobs().to_vec();
    // hook 1 = entry; 2 = objective fit; 3 = constraint fit.
    let mut checks = 0;
    assert_eq!(s.refill_controlled(&mut || { checks += 1; checks < 3 }), Err(AsyncConstrainedError::Cancelled));
    assert_eq!(checks, 3);
    assert_eq!(s.jobs(), received.as_slice());
    assert!(s.model_report().is_none());
    assert_eq!(s.work(), NoisyStudyWork { model_fit_attempts: 1, acquisition_batch_attempts: 0 });
    // entry + two output checks + fit completion + proposal + selection entry
    // + first slot + second slot. At 8 one greedy search has been discarded.
    checks = 0;
    assert_eq!(s.refill_controlled(&mut || { checks += 1; checks < 8 }), Err(AsyncConstrainedError::Cancelled));
    assert_eq!(checks, 8);
    assert_eq!(s.jobs(), received.as_slice());
    assert!(s.model_report().is_none());
    assert_eq!(s.work(), NoisyStudyWork { model_fit_attempts: 2, acquisition_batch_attempts: 1 });
    s.refill().unwrap(); clean.refill().unwrap();
    assert_eq!(s.jobs(), clean.jobs());
    assert_eq!(s.model_report(), clean.model_report());
}

#[test]
fn g4_a_failed_constraint_model_does_not_publish_a_partial_objective_update() {
    let mut c = config(); c.search.bounds = (1.0, 1.0 + f64::EPSILON);
    let mut s = AsyncConstrainedStudy::new(23, 1, 3, 3, 3, &c).unwrap();
    for r in s.ask() {
        // Objective covariance is admissible, constraint covariance is singular.
        s.tell(&r, ConstrainedObservation { objective: obs(1.0, 0.1),
            constraints: vec![obs(-1.0, 0.0)] }).unwrap();
    }
    let original = s.jobs().to_vec();
    assert_eq!(s.refill(), Err(AsyncConstrainedError::Model(ConstrainedBoError::Model {
        output: 1, source: HeteroFitError::InvalidInitialModel })));
    assert!(s.model_report().is_none() && !s.is_complete());
    assert_eq!(s.jobs(), original.as_slice());
    assert!(s.ask().is_empty());
}

#[test]
fn g5_availability_order_clone_replay_and_nonmultiple_budget_caps() {
    let mut c = config(); c.search.q = 4;
    let mut a = AsyncConstrainedStudy::new(29, 1, 5, 8, 2, &c).unwrap();
    let r = a.ask(); receive(&mut a, &r[1]);
    let mut b = a.clone();
    while !a.is_complete() {
        for r in a.ask() { receive(&mut a, &r); }
        for r in b.ask().into_iter().rev() { receive(&mut b, &r); }
        assert_eq!(a.refill().unwrap(), b.refill().unwrap());
        assert!(a.in_flight() <= 2 && a.jobs().len() <= 8);
        assert_eq!(a.jobs(), b.jobs());
        assert_eq!(a.model_report(), b.model_report());
    }
    assert!(b.is_complete());
    assert_eq!(a.jobs().len(), 8);
    assert_eq!(a.jobs().iter().filter(|j| j.observation.is_some()).count(), 8);
}

#[test]
fn g0_async_learning_and_bad_budgets_refuse_before_any_jobs_are_created() {
    let c = config();
    for (n, total, cap) in [(0, 4, 2), (3, 2, 2), (2, 4, 0), (2, usize::MAX, 3)] {
        assert!(matches!(AsyncConstrainedStudy::new(1, 1, n, total, cap, &c), Err(AsyncConstrainedError::InvalidPolicy)));
    }
    let mut c = c;
    let learning = crate::learning::NoisyLearningConfig { refit_every: 1,
        fit: crate::hyper::HeteroFitConfig { lengthscale_bounds: vec![(0.05, 1.0)], signal_bounds: (0.05, 2.0),
            starts: 1, max_iterations: 2, max_evaluations: 3, gradient_tolerance: 1e-5, seed: 7 } };
    c.objective_learning = Some(learning.clone());
    assert!(matches!(AsyncConstrainedStudy::new(1, 1, 2, 4, 2, &c), Err(AsyncConstrainedError::LearningUnsupported)));
    c.objective_learning = None; c.constraints[0].learning = Some(learning);
    assert!(matches!(AsyncConstrainedStudy::new(1, 1, 2, 4, 2, &c), Err(AsyncConstrainedError::LearningUnsupported)));
}
