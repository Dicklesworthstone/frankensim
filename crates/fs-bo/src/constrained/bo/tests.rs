use super::*;
use crate::{Matern, hyper::HeteroFitConfig};

fn kernel() -> Kernel {
    Kernel { family: Matern::FiveHalves, signal: 1.0, lengthscales: vec![0.3] }
}
fn config() -> ConstrainedBoConfig {
    ConstrainedBoConfig {
        search: NoisyBoConfig { bounds: (0.0, 1.0), kernel: kernel(), prior_mean: 0.0,
            q: 2, mc_samples: 32, acq_starts: 1, acq_evals: 24, seed: 91 },
        objective_learning: None,
        constraints: vec![OutcomeConstraint { kernel: kernel(), prior_mean: 0.0, upper_bound: 0.0, learning: None }],
        reference: 2.0, recommendation_probability: 0.9,
    }
}
fn obs(value: f64, noise_variance: f64) -> NoisyObservation { NoisyObservation { value, noise_variance } }
fn evaluation(x: &[f64]) -> ConstrainedObservation {
    ConstrainedObservation { objective: obs((x[0] - 0.2).powi(2), 0.001),
        constraints: vec![obs(0.6 - x[0], 0.001 + 0.003 * x[0])] }
}
fn learning(refit_every: usize, max_evaluations: usize) -> NoisyLearningConfig {
    NoisyLearningConfig { refit_every, fit: HeteroFitConfig {
        lengthscale_bounds: vec![(0.05, 2.0)], signal_bounds: (0.01, 5.0), starts: 1,
        max_iterations: 8, max_evaluations, gradient_tolerance: 1e-7, seed: 93,
    }}
}

#[test]
fn g3_every_output_and_recommendation_uses_complete_observation_history() {
    let c = config();
    let mut calls = Vec::new();
    let r = minimize_constrained(&mut |x| { calls.push(x.to_vec()); evaluation(x) }, 1, 3, 2, &c).unwrap();
    assert_eq!(calls, r.x);
    assert_eq!(r.x.len(), 7);
    assert_eq!(r.recommendations.len(), 3);
    assert_eq!(r.posterior_only_fits, 6); // two outputs, three boundaries
    assert_eq!(r.likelihood_evaluations, 0);
    for (x, observed) in r.x.iter().zip(&r.observations) { assert_eq!(*observed, evaluation(x)); }
    for (stage, recommendation) in r.recommendations.iter().enumerate() {
        let n = 3 + stage * c.search.q;
        let mut models = Vec::new();
        for (output, policy) in c.outputs().into_iter().enumerate() {
            let observations: Vec<_> = r.observations[..n].iter().map(|o|
                if output == 0 { o.objective } else { o.constraints[output - 1] }).collect();
            let y: Vec<_> = observations.iter().map(|o| o.value - policy.prior).collect();
            let noise: Vec<_> = observations.iter().map(|o| o.noise_variance).collect();
            models.push(Gp::try_fit_diag(&r.x[..n], &y, policy.initial.clone(), &noise).unwrap());
        }
        // Independent marginal reconstruction: do not use the recommendation
        // helper to decide eligibility or choose an index in this regression.
        let mut best: Option<(usize, f64, f64, f64)> = None;
        for i in 0..n {
            let (mean, variance) = models[0].predict(&r.x[i]);
            let (g, gv) = models[1].predict(&r.x[i]);
            let p = crate::phi_cdf(-g / fs_math::det::sqrt(gv)).clamp(0.0, 1.0);
            if p >= c.recommendation_probability && best.is_none_or(|(_, old, _, _)| mean < old) {
                best = Some((i, mean, variance, p));
            }
        }
        match (recommendation, best) {
            (None, None) => (),
            (Some(actual), Some((i, mean, variance, p))) => {
                assert_eq!(actual.observation_index, i);
                assert_eq!(actual.objective_mean, mean);
                assert_eq!(actual.objective_variance, variance);
                assert_eq!(actual.joint_feasibility, p);
            }
            pair => panic!("recommendation mismatch: {pair:?}"),
        }
    }
}

#[test]
fn g0_recommendation_rejects_either_violated_constraint_not_just_bad_objectives() {
    let mut c = config();
    c.search.kernel.lengthscales = vec![1e-8];
    c.constraints[0].kernel.lengthscales = vec![1e-8];
    c.constraints.push(c.constraints[0].clone());
    let mut call = 0;
    let r = minimize_constrained(&mut |_| {
        let (f, g, h) = match call { 0 => (-100.0, 10.0, -1.0), 1 => (1.0, -1.0, -1.0), _ => (-50.0, -1.0, 10.0) };
        call += 1;
        ConstrainedObservation { objective: obs(f, 1e-6), constraints: vec![obs(g, 1e-6), obs(h, 1e-6)] }
    }, 1, 3, 0, &c).unwrap();
    let recommended = r.recommendations[0].as_ref().unwrap();
    assert_eq!(recommended.observation_index, 1);
    assert!(recommended.constraint_means.iter().all(|v| *v < 0.0));
    assert!(recommended.joint_feasibility > 0.99);
    assert_eq!(r.observations[0].objective.value, -100.0);
}

#[test]
fn g0_unresolved_feasibility_returns_none_even_with_a_negative_mean() {
    let c = config();
    for (value, noise) in [(10.0, 1e-6), (-0.01, 1.0)] {
        let r = minimize_constrained(&mut |_| ConstrainedObservation {
            objective: obs(-100.0, 0.01), constraints: vec![obs(value, noise)],
        }, 1, 1, 0, &c).unwrap();
        assert!(r.recommendations[0].is_none());
    }
}

#[test]
fn g3_original_units_are_centered_consistently_in_constraints_and_reference() {
    let c = config();
    let mut shifted = c.clone();
    shifted.search.prior_mean = 128.0;
    shifted.reference += 128.0;
    shifted.constraints[0].prior_mean = 512.0;
    shifted.constraints[0].upper_bound += 512.0;
    let a = minimize_constrained(&mut |_| ConstrainedObservation {
        objective: obs(0.25, 0.001), constraints: vec![obs(-0.5, 0.001)],
    }, 1, 3, 1, &c).unwrap();
    let b = minimize_constrained(&mut |_| ConstrainedObservation {
        objective: obs(128.25, 0.001), constraints: vec![obs(511.5, 0.001)],
    }, 1, 3, 1, &shifted).unwrap();
    assert_eq!(a.x, b.x, "adding physical offsets cannot change the centered acquisition");
    assert_eq!(a.observations.len(), b.observations.len());
    for (old, new) in a.observations.iter().zip(&b.observations) {
        assert_eq!(old.objective.noise_variance, new.objective.noise_variance);
        assert_eq!(old.constraints[0].noise_variance, new.constraints[0].noise_variance);
    }
}

#[test]
fn g3_independent_refit_cadences_and_one_probe_fits_preserve_fixed_trajectory() {
    let c = config();
    let mut learned = c.clone();
    learned.objective_learning = Some(learning(2, 1));
    learned.constraints[0].learning = Some(learning(1, 1));
    let fixed = minimize_constrained(&mut evaluation, 1, 3, 2, &c).unwrap();
    let r = minimize_constrained(&mut evaluation, 1, 3, 2, &learned).unwrap();
    assert_eq!(r.x, fixed.x);
    assert_eq!(r.observations, fixed.observations);
    assert_eq!(r.recommendations, fixed.recommendations);
    assert_eq!(r.fits.iter().map(|f| (f.output, f.fit.after_batches)).collect::<Vec<_>>(),
        vec![(0, 0), (1, 0), (1, 1), (0, 2), (1, 2)]);
    assert_eq!(r.likelihood_evaluations, 5);
    assert_eq!(r.posterior_only_fits, 1);
    assert!(r.fits.iter().all(|f| f.fit.evaluation_limit_reached));
    assert_ne!(r.fits[0].fit.seed, r.fits[1].fit.seed);
}

#[test]
fn g5_active_objective_and_constraint_learning_replay_without_replacing_noise() {
    let mut c = config();
    c.objective_learning = Some(learning(1, 12));
    c.constraints[0].learning = Some(learning(1, 12));
    let a = minimize_constrained(&mut evaluation, 1, 3, 1, &c).unwrap();
    let b = minimize_constrained(&mut evaluation, 1, 3, 1, &c).unwrap();
    assert_eq!(a.x, b.x); assert_eq!(a.observations, b.observations); assert_eq!(a.recommendations, b.recommendations);
    assert_eq!(a.likelihood_evaluations, b.likelihood_evaluations);
    assert!(a.likelihood_evaluations <= 48);
    assert!(a.fits.iter().any(|f| f.fit.lml > f.fit.initial_lml), "learning must be active");
    for (fa, fb) in a.fits.iter().zip(&b.fits) {
        assert_eq!(fa.fit.kernel.signal.to_bits(), fb.fit.kernel.signal.to_bits());
        assert_eq!(fa.fit.kernel.lengthscales, fb.fit.kernel.lengthscales);
        assert!(fa.fit.lml >= fa.fit.initial_lml);
        let n = fa.fit.observation_count;
        let observed: Vec<_> = a.observations[..n].iter().map(|o|
            if fa.output == 0 { o.objective } else { o.constraints[fa.output - 1] }).collect();
        let y: Vec<_> = observed.iter().map(|o| o.value).collect();
        let noise: Vec<_> = observed.iter().map(|o| o.noise_variance).collect();
        let direct = Gp::try_fit_diag(&a.x[..n], &y, fa.fit.kernel.clone(), &noise).unwrap();
        assert_eq!(direct.lml.to_bits(), fa.fit.lml.to_bits());
    }
}

#[test]
fn g0_invalid_policies_refuse_before_callbacks_and_bad_outputs_stop_immediately() {
    for variant in 0..7 {
        let mut c = config();
        match variant {
            0 => c.constraints[0].upper_bound = f64::NAN,
            1 => c.constraints[0].kernel.lengthscales.clear(),
            2 => c.recommendation_probability = 0.0,
            3 => c.reference = f64::INFINITY,
            4 => c.constraints[0].learning = Some(learning(0, 12)),
            5 => c.search.q = usize::MAX,
            _ => c.objective_learning = Some(learning(1, usize::MAX)),
        }
        let mut calls = 0;
        let error = minimize_constrained(&mut |x| { calls += 1; evaluation(x) }, 1, 3, 2, &c).unwrap_err();
        assert_eq!(error, ConstrainedBoError::InvalidConfig);
        assert_eq!(calls, 0);
    }
    for bad in [ConstrainedObservation { objective: obs(0.0, 0.1), constraints: vec![] },
        ConstrainedObservation { objective: obs(0.0, 0.1), constraints: vec![obs(0.0, -1.0)] }] {
        let mut calls = 0;
        let error = minimize_constrained(&mut |_| { calls += 1; bad.clone() }, 1, 3, 2, &config()).unwrap_err();
        assert_eq!(calls, 1);
        assert!(matches!(error, ConstrainedBoError::OutputCount { evaluation: 0, .. }
            | ConstrainedBoError::Observation { evaluation: 0, output: 1 }));
    }
}

#[test]
fn g0_singular_constraint_fit_is_not_repaired_by_inventing_measurement_noise() {
    let c = config();
    let mut r = ConstrainedBoReport { x: vec![vec![0.0], vec![-0.0]],
        observations: vec![ConstrainedObservation { objective: obs(1.0, 0.1), constraints: vec![obs(0.0, 0.0)] }; 2],
        recommendations: Vec::new(), fits: Vec::new(), likelihood_evaluations: 0, posterior_only_fits: 0 };
    let mut kernels = vec![kernel(), kernel()];
    let failure = condition(&mut r, &c, &mut kernels, 0);
    assert!(matches!(failure, Err(ConstrainedBoError::Model { output: 1, source: HeteroFitError::InvalidInitialModel })));
    for observed in &mut r.observations { observed.constraints[0].noise_variance = 0.1; }
    assert!(condition(&mut r, &c, &mut kernels, 0).is_ok(), "noisy replicate constraints remain observations");
}
