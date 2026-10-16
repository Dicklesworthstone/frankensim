use fs_eproc::bernoulli::{BernoulliCsError, BernoulliMixtureCs};
use fs_eproc::GaussianMixtureCs;
use fs_math::det;

fn state(successes: usize, failures: usize, alpha: f64) -> BernoulliMixtureCs {
    let mut cs = BernoulliMixtureCs::new(alpha).unwrap();
    for _ in 0..successes { cs.observe(true).unwrap(); }
    for _ in 0..failures { cs.observe(false).unwrap(); }
    cs
}
fn factorial(n: usize) -> f64 { (1..=n).map(|v| v as f64).product() }

#[test]
fn likelihood_matches_independent_half_integer_beta_identity() {
    for n in 1..=12 {
        for s in 0..=n {
            let f = n - s;
            let marginal = factorial(2*s)*factorial(2*f)
                / (4.0_f64.powi(n as i32)*factorial(s)*factorial(f)*factorial(n));
            let cs = state(s, f, 0.05);
            for p in [0.01_f64, 0.3, 0.5, 0.9, 0.99] {
                let expected = marginal / (p.powi(s as i32)*(1.0-p).powi(f as i32));
                let actual = det::exp(cs.log_e_value(p).unwrap());
                assert!((actual/expected - 1.0).abs() < 2e-12, "{n}/{s} at {p}");
            }
        }
    }
}

#[test]
fn asymmetric_bounds_retain_probability_endpoints_and_nonzero_uncertainty() {
    let alpha = 0.05;
    let yes = state(1024, 0, alpha);
    let no = state(0, 1024, alpha);
    let a = yes.interval().unwrap().unwrap();
    let b = no.interval().unwrap().unwrap();
    let expected_lower = det::exp((yes.log_e_value(1.0).unwrap() + det::ln(alpha))/1024.0);
    assert!((a.lo - expected_lower).abs() < 2e-14);
    assert_eq!(a.hi, 1.0); assert_eq!(b.lo, 0.0);
    assert!((a.lo + b.hi - 1.0).abs() < 2e-14);
    assert!(a.lo > 0.99 && a.lo < 1.0);
    assert!(b.hi > 0.0 && b.hi < 0.01);
    let mut generic = GaussianMixtureCs::new(0.5, 1.0, alpha);
    for _ in 0..1024 { generic.observe(1.0); }
    let (mean, radius) = generic.interval().unwrap();
    assert!(mean-radius < 0.99);
}

#[test]
fn mixed_bounds_bracket_the_likelihood_threshold_on_both_sides() {
    for (s, f) in [(1, 1), (1, 99), (73, 27), (99, 1)] {
        let cs = state(s, f, 0.05);
        let interval = cs.interval().unwrap().unwrap();
        assert!(interval.lo < interval.mean && interval.mean < interval.hi);
        for p in [interval.lo, interval.hi] {
            assert!((cs.log_e_value(p).unwrap() + det::ln(0.05)).abs() < 1e-10);
        }
        assert!(cs.log_e_value(interval.mean).unwrap() < 0.0);
    }
}

#[test]
fn finite_horizon_crossing_probability_sums_all_binary_paths() {
    fn crossing(cs: BernoulliMixtureCs, p: f64, mass: f64, left: usize, threshold: f64) -> f64 {
        if cs.log_e_value(p).unwrap() >= threshold { return mass; }
        if left == 0 { return 0.0; }
        let mut yes = cs.clone(); yes.observe(true).unwrap();
        let mut no = cs; no.observe(false).unwrap();
        crossing(yes, p, mass*p, left-1, threshold) + crossing(no, p, mass*(1.0-p), left-1, threshold)
    }
    for p in [0.01, 0.1, 0.3, 0.5, 0.9, 0.99] {
        for alpha in [0.05, 0.2, 0.8] {
            let total = crossing(BernoulliMixtureCs::new(alpha).unwrap(), p, 1.0, 12, -det::ln(alpha));
            assert!(total <= alpha + 1e-12, "p={p}, alpha={alpha}: {total}");
        }
    }
}

#[test]
fn clone_resume_and_endpoint_nulls_preserve_the_original_process() {
    let mut full = BernoulliMixtureCs::new(0.05).unwrap();
    assert!(full.is_empty()); assert_eq!(full.interval().unwrap(), None);
    assert_eq!(full.log_e_value(0.0).unwrap(), 0.0);
    let mut split = full.clone();
    for i in 0..70 { full.observe(i%7 != 0).unwrap(); }
    for i in 0..17 { split.observe(i%7 != 0).unwrap(); }
    let mut resumed = split.clone();
    for i in 17..70 { resumed.observe(i%7 != 0).unwrap(); }
    assert_eq!(full, resumed); assert_eq!(full.interval(), resumed.interval());
    assert_eq!(full.log_e_value(0.0).unwrap(), f64::INFINITY);
    assert_eq!(full.log_e_value(1.0).unwrap(), f64::INFINITY);
    assert_eq!(full.len(), 70); assert_eq!(full.successes(), 60);
}

#[test]
fn invalid_inputs_refuse_and_tiny_alpha_does_not_overflow_a_reciprocal() {
    for alpha in [0.0, -0.1, 1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(BernoulliMixtureCs::new(alpha), Err(BernoulliCsError::InvalidAlpha));
    }
    let cs = state(1, 1, 0.05);
    for p in [-0.1, 1.1, f64::NAN, f64::INFINITY] {
        assert_eq!(cs.log_e_value(p), Err(BernoulliCsError::InvalidProbability));
    }
    let narrow_alpha = state(1, 0, 1e-300).interval().unwrap().unwrap();
    assert_eq!(narrow_alpha.hi, 1.0);
    assert!(narrow_alpha.lo < 1e-20); // bounded inversion may conservatively keep zero
}
