use super::*;

fn kernel(family: Matern) -> Kernel {
    Kernel { family, signal: 1.3, lengthscales: vec![0.4, 0.8] }
}

fn policy() -> HeteroFitConfig {
    HeteroFitConfig {
        lengthscale_bounds: vec![(0.05, 3.0), (0.1, 4.0)],
        signal_bounds: (0.05, 5.0), starts: 2, max_iterations: 40,
        max_evaluations: 100, gradient_tolerance: 1e-7, seed: 93,
    }
}

fn fixture() -> (Vec<Vec<f64>>, Vec<f64>, Vec<f64>) {
    // Coincident points have different uncertainties; neither is deduplicated
    // in the observation likelihood. The noisy outlier must stay downweighted.
    (vec![vec![0.0, 0.1], vec![0.3, 0.7], vec![0.3, 0.7], vec![0.9, 0.2]],
        vec![0.2, 0.8, -5.0, -0.3], vec![0.01, 0.04, 20.0, 0.09])
}

#[test]
fn g0_all_matern_log_gradients_match_independent_refits() {
    let (x, y, noise) = fixture();
    for family in [Matern::Half, Matern::ThreeHalves, Matern::FiveHalves] {
        let k = kernel(family);
        let evaluation = evaluate(&x, &y, &noise, k.clone(), &mut || true).unwrap().unwrap();
        for coordinate in 0..3 {
            let mut plus = k.clone();
            let mut minus = k.clone();
            let eps = 1e-5f64;
            if coordinate < 2 {
                plus.lengthscales[coordinate] *= eps.exp();
                minus.lengthscales[coordinate] *= (-eps).exp();
            } else {
                plus.signal *= eps.exp();
                minus.signal *= (-eps).exp();
            }
            let fp = Gp::try_fit_diag(&x, &y, plus, &noise).unwrap().lml;
            let fm = Gp::try_fit_diag(&x, &y, minus, &noise).unwrap().lml;
            let fd = (fp - fm) / (2.0 * eps);
            assert!((fd - evaluation.gradient[coordinate]).abs() < 2e-5 * (1.0 + fd.abs()),
                "{family:?} coordinate {coordinate}: analytic={} fd={fd}", evaluation.gradient[coordinate]);
        }
        assert!(evaluation.gradient.iter().any(|g| g.abs() > 0.01), "nonvacuous derivatives");
    }
}

#[test]
fn g1_one_point_recovers_closed_form_signal_variance() {
    // C = signal + noise. LML' = (y^2-C)/(2 C^2), so signal = y^2-noise = 2.
    let k = Kernel { family: Matern::FiveHalves, signal: 0.4, lengthscales: vec![0.7] };
    let config = HeteroFitConfig { lengthscale_bounds: vec![(0.7, 0.7)],
        signal_bounds: (0.1, 4.0), starts: 1, max_iterations: 200,
        max_evaluations: 400, gradient_tolerance: 1e-8, seed: 7 };
    let report = fit_heteroscedastic(&[vec![0.0]], &[1.5], &[0.25], &k, &config).unwrap();
    assert!((report.model.kernel.signal - 2.0).abs() < 1e-5);
    assert_eq!(report.model.kernel.lengthscales, vec![0.7]);
    assert!(report.model.lml > report.initial_lml + 0.5);
    assert!(report.projected_gradient_inf < 1e-6);
    assert!(report.evaluations <= config.max_evaluations);
}

#[test]
fn g3_training_preserves_every_observation_variance() {
    let (x, y, noise) = fixture();
    let original_noise = noise.clone();
    let report = fit_heteroscedastic(&x, &y, &noise, &kernel(Matern::FiveHalves), &policy()).unwrap();
    let direct = Gp::try_fit_diag(&x, &y, report.model.kernel.clone(), &noise).unwrap();
    assert_eq!(noise, original_noise);
    assert_eq!(report.model.lml.to_bits(), direct.lml.to_bits());
    for point in [vec![0.1, 0.4], vec![0.3, 0.7], vec![0.8, 0.9]] {
        assert_eq!(report.model.predict(&point), direct.predict(&point));
    }
    assert!(report.model.lml >= report.initial_lml);
    assert!(report.model.predict(&[0.3, 0.7]).0 > 0.0, "noisy outlier must not dominate");
}

#[test]
fn g0_budget_and_frozen_coordinates_preserve_initial_model() {
    let (x, y, noise) = fixture();
    let k = kernel(Matern::ThreeHalves);
    let mut c = policy();
    c.max_evaluations = 1;
    let report = fit_heteroscedastic(&x, &y, &noise, &k, &c).unwrap();
    assert_eq!(report.evaluations, 1);
    assert_eq!(report.starts_attempted, 1);
    assert_eq!(report.accepted_steps, 0);
    assert_eq!(report.model.kernel.signal.to_bits(), k.signal.to_bits());
    assert_eq!(report.model.kernel.lengthscales, k.lengthscales);
    assert_eq!(report.model.lml.to_bits(), report.initial_lml.to_bits());
    assert!(report.evaluation_limit_reached);
    assert!(report.projected_gradient_inf > 1e-5, "budget exhaustion is not convergence");

    c.max_evaluations = 100;
    c.starts = 1;
    c.lengthscale_bounds = vec![(0.4, 0.4), (0.8, 0.8)];
    c.signal_bounds = (1.3, 1.3);
    let frozen = fit_heteroscedastic(&x, &y, &noise, &k, &c).unwrap();
    assert_eq!(frozen.evaluations, 1);
    assert_eq!(frozen.projected_gradient_inf, 0.0);
}

#[test]
fn g0_best_restart_is_retained_not_just_the_last_model() {
    let (x, y, noise) = fixture();
    let k = kernel(Matern::Half);
    let mut c = policy();
    c.starts = 5;
    c.max_iterations = 0;
    c.max_evaluations = 5;
    let result = fit_heteroscedastic(&x, &y, &noise, &k, &c).unwrap();
    let mut best = Gp::try_fit_diag(&x, &y, k.clone(), &noise).unwrap().lml;
    let bounds = c.log_bounds();
    for start in 1..5 {
        let mut rng = fs_rand::StreamKey { seed: c.seed, kernel: 0x4845_5446, tile: start }.stream();
        let p: Vec<f64> = bounds.iter().map(|(lo, hi)| (hi - lo).mul_add(rng.next_f64(), *lo)).collect();
        best = best.max(Gp::try_fit_diag(&x, &y, c.kernel(k.family, &p), &noise).unwrap().lml);
    }
    assert_eq!(result.evaluations, 5);
    assert_eq!(result.starts_attempted, 5);
    assert_eq!(result.model.lml.to_bits(), best.to_bits());
}

#[test]
fn g0_admission_and_singular_covariances_do_not_invent_noise() {
    let (x, y, noise) = fixture();
    let k = kernel(Matern::Half);
    let mut checks = 0;
    let mut stop = || { checks += 1; true };
    assert!(matches!(fit_heteroscedastic_controlled(&x, &y, &[-1.0; 4], &k, &policy(), &mut stop),
        Err(HeteroFitError::InvalidData)));
    assert_eq!(checks, 0);
    let mut c = policy();
    c.signal_bounds = (0.1, 0.2);
    assert!(matches!(fit_heteroscedastic(&x, &y, &noise, &k, &c), Err(HeteroFitError::InvalidConfig)));
    assert!(matches!(fit_heteroscedastic(&x, &y, &[0.0; 4], &k, &policy()),
        Err(HeteroFitError::InvalidInitialModel)), "duplicate zero-noise rows must refuse");
    let distinct = vec![x[0].clone(), x[1].clone(), x[3].clone()];
    let exact = fit_heteroscedastic(&distinct, &[0.2, 0.8, -0.3], &[0.0; 3], &k, &policy()).unwrap();
    assert_eq!(exact.model.noise, 0.0);
}

#[test]
fn g4_stop_hook_refuses_inside_gradient_work() {
    let (x, y, noise) = fixture();
    let mut checks = 0;
    let result = fit_heteroscedastic_controlled(&x, &y, &noise, &kernel(Matern::FiveHalves),
        &policy(), &mut || { checks += 1; checks < 3 });
    assert!(matches!(result, Err(HeteroFitError::Cancelled)));
    assert_eq!(checks, 3);
}

#[test]
fn g5_entire_fit_replays_and_reports_returned_model_stationarity() {
    let (x, y, noise) = fixture();
    let k = kernel(Matern::FiveHalves);
    let c = policy();
    let a = fit_heteroscedastic(&x, &y, &noise, &k, &c).unwrap();
    let b = fit_heteroscedastic(&x, &y, &noise, &k, &c).unwrap();
    assert_eq!(a.model.kernel.signal.to_bits(), b.model.kernel.signal.to_bits());
    assert_eq!(a.model.kernel.lengthscales, b.model.kernel.lengthscales);
    assert_eq!(a.model.lml.to_bits(), b.model.lml.to_bits());
    assert_eq!((a.evaluations, a.accepted_steps, a.starts_attempted),
        (b.evaluations, b.accepted_steps, b.starts_attempted));
    assert_eq!(a.projected_gradient_inf.to_bits(), b.projected_gradient_inf.to_bits());
    let scored = evaluate(&x, &y, &noise, a.model.kernel.clone(), &mut || true).unwrap().unwrap();
    assert_eq!(a.projected_gradient_inf,
        projected_norm(&parameters(&a.model.kernel), &scored.gradient, &c.log_bounds()));
}

#[test]
fn g0_outward_boundary_gradient_is_not_a_false_stationarity_failure() {
    let k = Kernel { family: Matern::Half, signal: 1.0, lengthscales: vec![0.5] };
    let c = HeteroFitConfig { lengthscale_bounds: vec![(0.5, 0.5)], signal_bounds: (0.1, 1.0),
        starts: 1, max_iterations: 10, max_evaluations: 20, gradient_tolerance: 1e-8, seed: 0 };
    let report = fit_heteroscedastic(&[vec![0.0]], &[10.0], &[0.1], &k, &c).unwrap();
    assert_eq!(report.model.kernel.signal, 1.0);
    assert_eq!(report.projected_gradient_inf, 0.0);
    assert_eq!(report.evaluations, 1);
}
