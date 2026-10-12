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

fn learned_config() -> ConstrainedBoConfig {
    let mut c = config();
    let policy = crate::learning::NoisyLearningConfig { refit_every: 2,
        fit: crate::hyper::HeteroFitConfig { lengthscale_bounds: vec![(0.05, 2.0)],
            signal_bounds: (0.01, 2.0), starts: 1, max_iterations: 5,
            max_evaluations: 8, gradient_tolerance: 1e-7, seed: 713 } };
    c.objective_learning = Some(policy.clone());
    c.constraints[0].learning = Some(crate::learning::NoisyLearningConfig {
        refit_every: 3, ..policy
    });
    c
}

fn kernel_bits(kernels: &[Kernel]) -> Vec<u64> {
    kernels.iter().flat_map(|k| std::iter::once(k.signal.to_bits())
        .chain(k.lengthscales.iter().map(|x| x.to_bits()))).collect()
}

#[test]
fn g3_async_output_learning_follows_data_revisions_and_exact_warm_starts() {
    let c = learned_config();
    let mut s = AsyncConstrainedStudy::new_with_learning(41, 1, 4, 4, 4, &c).unwrap();
    let requests = s.ask();
    let mut kernels = c.initial_kernels();
    let mut reference = ConstrainedBoReport { x: Vec::new(), observations: Vec::new(),
        recommendations: Vec::new(), fits: Vec::new(), likelihood_evaluations: 0, posterior_only_fits: 0 };
    for (stage, id) in [2, 0, 3, 1].into_iter().enumerate() {
        receive(&mut s, &requests[id]);
        reference.x = s.jobs().iter().filter(|j| j.observation.is_some())
            .map(|j| j.request.x.clone()).collect();
        reference.observations = s.jobs().iter().filter_map(|j| j.observation.clone()).collect();
        let models = condition_controlled(&mut reference, &c, &mut kernels, stage,
            &mut || true).unwrap().unwrap();
        assert_eq!(s.refill(), Ok(0));
        assert_eq!(s.model_updates(), stage + 1);
        assert_eq!(kernel_bits(s.kernels()), kernel_bits(&kernels));
        assert_eq!(s.model_report().unwrap().recommendation, recommend(&models, &reference.x, &c).unwrap());
        assert_eq!(s.likelihood_evaluations(), reference.likelihood_evaluations);
        assert_eq!(s.posterior_only_fits(), reference.posterior_only_fits);
    }
    assert_eq!(reference.fits.iter().map(|r| (r.fit.after_batches, r.output)).collect::<Vec<_>>(),
        vec![(0, 0), (0, 1), (2, 0), (3, 1)]);
    assert!(reference.fits.iter().any(|r| r.fit.lml > r.fit.initial_lml), "learning is active");
    assert!(s.is_complete());
    assert!(s.likelihood_evaluations() <= 4 * 8);
}

#[test]
fn g3_one_worker_learned_async_trajectory_matches_synchronous_learning() {
    let mut c = learned_config(); c.search.q = 1;
    let reference = minimize_constrained(&mut evaluate, 1, 3, 2, &c).unwrap();
    let mut s = AsyncConstrainedStudy::new_with_learning(43, 1, 3, 5, 1, &c).unwrap();
    while !s.is_complete() {
        let request = s.ask()[0].clone(); receive(&mut s, &request); s.refill().unwrap();
    }
    assert_eq!(s.jobs().iter().map(|j| j.request.x.clone()).collect::<Vec<_>>(), reference.x);
    assert_eq!(s.jobs().iter().map(|j| j.observation.clone().unwrap()).collect::<Vec<_>>(), reference.observations);
    assert_eq!(&s.model_report().unwrap().recommendation, reference.recommendations.last().unwrap());
    assert_eq!(s.likelihood_evaluations(), reference.likelihood_evaluations);
    assert_eq!(s.posterior_only_fits(), reference.posterior_only_fits);
    assert_eq!(s.model_updates(), 3);
}

#[test]
fn g0_same_data_refills_do_not_retrain_or_advance_the_learning_schedule() {
    let mut c = learned_config(); c.search.q = 1;
    c.objective_learning.as_mut().unwrap().refit_every = 1;
    c.constraints[0].learning.as_mut().unwrap().refit_every = 1;
    let mut s = AsyncConstrainedStudy::new_with_learning(47, 1, 2, 6, 4, &c).unwrap();
    for request in s.ask() { receive(&mut s, &request); }
    assert_eq!(s.refill(), Ok(1));
    let kernels = kernel_bits(s.kernels());
    let learning_work = s.likelihood_evaluations();
    let conditioning = s.posterior_only_fits();
    assert!(learning_work > 0);
    assert_eq!(s.refill(), Ok(1));
    assert_eq!(s.model_updates(), 1);
    assert_eq!(s.likelihood_evaluations(), learning_work);
    assert_eq!(s.posterior_only_fits(), conditioning + 2);
    assert_eq!(kernel_bits(s.kernels()), kernels);
    assert_eq!(s.model_report().unwrap().observed_evaluations, vec![0, 1]);
    assert_eq!(s.jobs().len(), 4);
}

#[test]
fn g4_cancelled_learning_or_selection_preserves_accepted_kernels_and_cadence() {
    let c = learned_config();
    let mut s = AsyncConstrainedStudy::new_with_learning(53, 1, 3, 7, 3, &c).unwrap();
    for request in s.ask().into_iter().skip(1) { receive(&mut s, &request); }
    let mut clean = s.clone();
    let pending = s.jobs().to_vec();
    let kernels = kernel_bits(s.kernels());
    let mut checks = 0;
    assert_eq!(s.refill_controlled(&mut || { checks += 1; checks < 4 }),
        Err(AsyncConstrainedError::Cancelled)); // inside first learner's gradient
    assert_eq!(checks, 4);
    assert_eq!(kernel_bits(s.kernels()), kernels);
    assert_eq!(s.model_updates(), 0);
    assert_eq!(s.likelihood_evaluations(), 0);
    assert_eq!(s.jobs(), pending.as_slice());

    // Count the real fit boundaries, then refuse after one candidate slot.
    let mut trial = ConstrainedBoReport {
        x: pending.iter().filter(|j| j.observation.is_some()).map(|j| j.request.x.clone()).collect(),
        observations: pending.iter().filter_map(|j| j.observation.clone()).collect(),
        recommendations: Vec::new(), fits: Vec::new(), likelihood_evaluations: 0, posterior_only_fits: 0,
    };
    let mut fit_checks = 0;
    let _models = condition_controlled(&mut trial, &c, &mut c.initial_kernels(), 0,
        &mut || { fit_checks += 1; true }).unwrap().unwrap();
    checks = 0;
    assert_eq!(s.refill_controlled(&mut || { checks += 1; checks < fit_checks + 5 }),
        Err(AsyncConstrainedError::Cancelled));
    assert_eq!(s.jobs(), pending.as_slice());
    assert_eq!(kernel_bits(s.kernels()), kernels);
    assert_eq!(s.model_updates(), 0);
    assert!(s.model_report().is_none());
    assert_eq!(s.work().model_fit_attempts, 2);
    assert_eq!(s.work().acquisition_batch_attempts, 1);
    s.refill().unwrap(); clean.refill().unwrap();
    assert_eq!(s.jobs(), clean.jobs());
    assert_eq!(s.model_report(), clean.model_report());
    assert_eq!(kernel_bits(s.kernels()), kernel_bits(clean.kernels()));
    assert_eq!(s.model_updates(), clean.model_updates());
    assert_eq!(s.likelihood_evaluations(), clean.likelihood_evaluations());
}

#[test]
fn g0_async_learning_work_is_admitted_over_the_whole_received_data_budget() {
    let mut c = learned_config();
    c.constraints[0].learning = None;
    c.objective_learning.as_mut().unwrap().fit.max_evaluations = usize::MAX;
    assert!(matches!(AsyncConstrainedStudy::new_with_learning(59, 1, 2, 4, 2, &c),
        Err(AsyncConstrainedError::InvalidPolicy)));
    c.objective_learning.as_mut().unwrap().fit.max_evaluations = 8;
    c.objective_learning.as_mut().unwrap().refit_every = 0;
    assert!(matches!(AsyncConstrainedStudy::new_with_learning(59, 1, 2, 4, 2, &c),
        Err(AsyncConstrainedError::Model(ConstrainedBoError::InvalidConfig))));
}
