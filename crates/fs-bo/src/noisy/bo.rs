//! A declared-noise BO consumer of joint-posterior q-NEI.

use super::{joint_normal_bank, q_noisy_expected_improvement};
use crate::gp::{Gp, Kernel};

/// One stochastic objective evaluation in the objective's original units.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NoisyObservation {
    /// Observed objective value (or an estimated mean).
    pub value: f64,
    /// Variance of THIS observation, not its standard deviation. For an
    /// estimated mean, supply the variance of the mean, not of raw samples.
    pub noise_variance: f64,
}

/// Explicit fixed-prior policy for noisy Bayesian optimization.
///
/// Signal and noise variances use squared objective units; ARD lengthscales
/// use the input units. No standardization, learned noise floor, or automatic
/// hyperparameter refit changes these declared quantities during the run.
#[derive(Debug, Clone)]
pub struct NoisyBoConfig {
    /// Common finite search interval for each coordinate.
    pub bounds: (f64, f64),
    /// Fixed GP covariance, including its positive signal variance.
    pub kernel: Kernel,
    /// Fixed prior mean in the original objective units.
    pub prior_mean: f64,
    /// Number of evaluations selected together in each adaptive batch.
    pub q: usize,
    /// Fixed normal-bank samples per joint acquisition evaluation.
    pub mc_samples: usize,
    /// CMA-ES restarts per greedy batch position.
    pub acq_starts: usize,
    /// CMA-ES objective-evaluation budget per restart.
    pub acq_evals: usize,
    /// Root seed for initialization, normal banks and acquisition search.
    pub seed: u64,
}

/// Posterior recommendation among inputs evaluated so far.
///
/// These are model estimates, NOT the smallest observation or a certificate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NoisyIncumbent {
    /// Index into the report's aligned input/observation history.
    pub observation_index: usize,
    /// Latent posterior mean in the original objective units.
    pub mean: f64,
    /// Latent posterior variance (measurement noise is not re-added).
    pub variance: f64,
}

/// Complete observed history and posterior recommendations.
#[derive(Debug, Clone, PartialEq)]
pub struct NoisyBoReport {
    /// All callback inputs in evaluation order.
    pub x: Vec<Vec<f64>>,
    /// Unmodified, aligned observations and their declared variances.
    pub observations: Vec<NoisyObservation>,
    /// One recommendation after initialization and after each batch. This is
    /// not a monotone best-observed trace: learning may revise an incumbent.
    pub incumbent_trace: Vec<NoisyIncumbent>,
}

fn validate(dim: usize, n_init: usize, iters: usize, c: &NoisyBoConfig) {
    assert!((1..=fs_rand::qmc::MAX_SOBOL_DIM).contains(&dim), "noisy BO dimension outside Sobol table");
    assert!(n_init > 0 && u32::try_from(n_init).is_ok(), "invalid noisy BO initial count");
    let (lo, hi) = c.bounds;
    let span = hi - lo;
    assert!(lo.is_finite() && hi.is_finite() && span.is_finite() && span > 0.0,
        "noisy BO bounds must have a positive finite width");
    assert!(0.2 * span > 0.0, "noisy BO search scale underflow");
    assert!(c.prior_mean.is_finite(), "noisy BO prior mean must be finite");
    assert!(c.kernel.signal.is_finite() && c.kernel.signal > 0.0, "invalid noisy BO signal variance");
    assert_eq!(c.kernel.lengthscales.len(), dim, "noisy BO kernel dimension mismatch");
    assert!(c.kernel.lengthscales.iter().all(|l| l.is_finite() && *l > 0.0), "invalid noisy BO lengthscale");
    assert!(c.q > 0, "noisy BO batch must be nonempty");
    assert!(c.mc_samples > 0 && u32::try_from(c.mc_samples).is_ok(), "invalid noisy BO sample count");
    assert!(c.acq_starts > 0 && u32::try_from(c.acq_starts).is_ok(), "invalid noisy BO acquisition starts");
    assert!(c.acq_evals > 0, "noisy BO acquisition budget must be positive");
    let total = iters.checked_mul(c.q).and_then(|n| n.checked_add(n_init))
        .expect("noisy BO observation count overflow");
    total.checked_mul(dim).expect("noisy BO input length overflow");
    total.checked_mul(total).expect("noisy BO covariance length overflow");
    total.checked_mul(c.mc_samples).expect("noisy BO sample-bank length overflow");
}

fn evaluate(f: &mut dyn FnMut(&[f64]) -> NoisyObservation, x: &[f64], prior: f64) -> NoisyObservation {
    let observation = f(x);
    assert!(observation.value.is_finite() && (observation.value - prior).is_finite(),
        "noisy BO observed/centered value must be finite");
    assert!(observation.noise_variance.is_finite() && observation.noise_variance >= 0.0,
        "noisy BO observation variance must be finite and nonnegative");
    observation
}

fn fit(x: &[Vec<f64>], observations: &[NoisyObservation], c: &NoisyBoConfig) -> Gp {
    let values: Vec<f64> = observations.iter().map(|o| o.value - c.prior_mean).collect();
    let variances: Vec<f64> = observations.iter().map(|o| o.noise_variance).collect();
    let gp = Gp::try_fit_diag(x, &values, c.kernel.clone(), &variances)
        .expect("noisy BO covariance must be SPD; no undeclared noise floor is added");
    assert!(gp.lml.is_finite(), "noisy BO fit must be finite");
    gp
}

fn recommend(gp: &Gp, x: &[Vec<f64>], prior: f64) -> NoisyIncumbent {
    let mut best: Option<NoisyIncumbent> = None;
    for (index, point) in x.iter().enumerate() {
        let (centered, variance) = gp.predict(point);
        let mean = centered + prior;
        assert!(mean.is_finite() && variance.is_finite(), "noisy BO posterior must be finite");
        if best.as_ref().is_none_or(|previous| mean < previous.mean) {
            best = Some(NoisyIncumbent { observation_index: index, mean, variance });
        }
    }
    best.expect("validated nonempty noisy BO history")
}

fn argmax(acquisition: &dyn Fn(&[f64]) -> f64, dim: usize, c: &NoisyBoConfig, seed: u64) -> Vec<f64> {
    let (lo, hi) = c.bounds;
    let sobol = fs_rand::qmc::Sobol::scrambled(dim, seed);
    let mut point = vec![0.0; dim];
    let mut best: Option<(f64, Vec<f64>)> = None;
    for start in 0..c.acq_starts {
        sobol.point(u32::try_from(start + 1).expect("validated starts"), &mut point);
        let initial: Vec<f64> = point.iter().map(|u| (hi - lo).mul_add(*u, lo)).collect();
        let mut objective = |x: &[f64]| {
            let bounded: Vec<f64> = x.iter().map(|v| v.clamp(lo, hi)).collect();
            let value = acquisition(&bounded);
            assert!(value.is_finite() && value >= 0.0, "invalid noisy BO acquisition value");
            -value
        };
        let params = fs_dfo::CmaParams::standard(dim, 0.2 * (hi - lo), c.acq_evals, f64::NEG_INFINITY);
        let report = fs_dfo::cmaes(&mut objective, &initial, &params,
            seed ^ (start as u64).wrapping_mul(0x9e37_79b9_7f4a_7c15));
        let value = -report.f_best;
        assert!(value.is_finite(), "noisy BO acquisition search must return a finite value");
        if best.as_ref().is_none_or(|(previous, _)| value > *previous) {
            let bounded = report.x_best.iter().map(|v| v.clamp(lo, hi)).collect();
            best = Some((value, bounded));
        }
    }
    best.expect("validated nonzero acquisition starts").1
}

/// Minimize a stochastic objective with declared per-observation variances.
///
/// Runs exactly `n_init + iters * config.q` objective callbacks on successful
/// completion. Each adaptive batch is selected before any callback from that
/// batch: sequential greedy q-NEI uses one fixed joint normal bank throughout
/// its CMA-ES searches. All evaluated inputs remain in the incumbent set.
///
/// The kernel and prior mean are fixed caller declarations in physical units.
/// Zero observation variance is accepted, but a singular training covariance
/// refuses rather than adding invented measurement noise. q-NEI measures new
/// latent improvement, not the value of information from replicate sampling.
///
/// Replay requires the same configuration AND the same callback outcomes.
/// This synchronous driver inherits the existing CMA-ES/GP work boundaries;
/// it adds no within-factorization cancellation or statistical stopping rule.
/// Dense joint posterior work scales with the complete observation history.
///
/// # Panics
/// Invalid configuration refuses before the first callback. Invalid observed
/// values/variances refuse at that callback; inadmissible GP fits or numerical
/// acquisition failures refuse without returning a successful report.
pub fn minimize_noisy(
    f: &mut dyn FnMut(&[f64]) -> NoisyObservation,
    dim: usize,
    n_init: usize,
    iters: usize,
    config: &NoisyBoConfig,
) -> NoisyBoReport {
    config.run_with_model(f, dim, n_init, iters, &mut |x, observations, _| {
        Ok::<_, std::convert::Infallible>(fit(x, observations, config))
    }).unwrap_or_else(|never| match never {})
}

impl NoisyBoConfig {
    pub(crate) fn validate_study(&self, dim: usize, n_init: usize, iters: usize) {
        validate(dim, n_init, iters, self);
    }

    pub(crate) fn initial_design(&self, dim: usize, n_init: usize) -> Vec<Vec<f64>> {
        let (lo, hi) = self.bounds;
        let sobol = fs_rand::qmc::Sobol::scrambled(dim, self.seed);
        let mut point = vec![0.0; dim];
        (0..n_init).map(|index| {
            sobol.point(u32::try_from(index + 1).expect("validated initial count"), &mut point);
            point.iter().map(|u| (hi - lo).mul_add(*u, lo)).collect()
        }).collect()
    }

    pub(crate) fn study_incumbent(&self, gp: &Gp, x: &[Vec<f64>]) -> NoisyIncumbent {
        recommend(gp, x, self.prior_mean)
    }

    // One seed policy and acquisition engine for callback and ask/tell studies.
    // Cancellation is checked between greedy slots, not inside a CMA-ES search.
    pub(crate) fn select_batch(&self, gp: &Gp, x: &[Vec<f64>], dim: usize,
        iteration: usize, keep_going: &mut dyn FnMut() -> bool) -> Option<Vec<Vec<f64>>>
    {
        if !keep_going() { return None; }
        let seed = self.seed ^ 0x4E45_4942 ^ (iteration as u64).wrapping_mul(0x9e37_79b9_7f4a_7c15);
        let width = x.len() + self.q; // admitted by total-count validation
        let bank = joint_normal_bank(self.mc_samples, width, seed);
        let mut batch: Vec<Vec<f64>> = Vec::new();
        for slot in 0..self.q {
            if !keep_going() { return None; }
            let active_width = x.len() + slot + 1;
            let mut prefix = Vec::with_capacity(self.mc_samples * active_width);
            for row in bank.chunks_exact(width) {
                prefix.extend_from_slice(&row[..active_width]);
            }
            let candidate = argmax(&|input: &[f64]| {
                let mut trial = batch.clone();
                trial.push(input.to_vec());
                q_noisy_expected_improvement(gp, x, &trial, &prefix)
            }, dim, self, seed ^ (slot as u64).wrapping_mul(0xD1B5_4A32_D192_ED03));
            batch.push(candidate);
        }
        if !keep_going() { return None; }
        Some(batch)
    }

    // Shared acquisition/history engine. Model fitting happens only at complete
    // batch boundaries; a failed fit cannot launch another objective callback.
    pub(crate) fn run_with_model<E>(
        &self,
        f: &mut dyn FnMut(&[f64]) -> NoisyObservation,
        dim: usize,
        n_init: usize,
        iters: usize,
        fit_model: &mut dyn FnMut(&[Vec<f64>], &[NoisyObservation], usize) -> Result<Gp, E>,
    ) -> Result<NoisyBoReport, E> {
    let config = self;
    validate(dim, n_init, iters, config);
    let mut x = Vec::new();
    let mut observations = Vec::new();
    for input in config.initial_design(dim, n_init) {
        observations.push(evaluate(f, &input, config.prior_mean));
        x.push(input);
    }
    let mut gp = fit_model(&x, &observations, 0)?;
    let mut incumbent_trace = vec![recommend(&gp, &x, config.prior_mean)];
    for iteration in 0..iters {
        let batch = config.select_batch(&gp, &x, dim, iteration, &mut || true)
            .expect("unconditional continuation");
        for input in batch {
            observations.push(evaluate(f, &input, config.prior_mean));
            x.push(input);
        }
        gp = fit_model(&x, &observations, iteration + 1)?;
        incumbent_trace.push(recommend(&gp, &x, config.prior_mean));
    }
    Ok(NoisyBoReport { x, observations, incumbent_trace })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gp::Matern;
    use std::panic::{AssertUnwindSafe, catch_unwind};

    fn config() -> NoisyBoConfig {
        NoisyBoConfig {
            bounds: (0.0, 1.0),
            kernel: Kernel { family: Matern::FiveHalves, signal: 1.0, lengthscales: vec![0.3] },
            prior_mean: 0.5, q: 2, mc_samples: 32, acq_starts: 1, acq_evals: 24, seed: 97,
        }
    }

    fn observation(x: &[f64]) -> NoisyObservation {
        NoisyObservation { value: (x[0] - 0.37).powi(2), noise_variance: 0.01 + 0.03 * x[0] }
    }

    #[test]
    fn batched_driver_retains_noise_and_replays_every_output() {
        let c = config();
        let mut calls = Vec::new();
        let result = minimize_noisy(&mut |x| { calls.push(x.to_vec()); observation(x) }, 1, 3, 2, &c);
        assert_eq!(calls, result.x);
        assert_eq!(result.x.len(), 7);
        assert_eq!(result.observations.len(), 7);
        assert_eq!(result.incumbent_trace.len(), 3);
        assert!(result.x.iter().flatten().all(|v| (0.0..=1.0).contains(v)));
        for (point, observed) in result.x.iter().zip(&result.observations) {
            assert_eq!(*observed, observation(point));
        }
        let replay = minimize_noisy(&mut observation, 1, 3, 2, &c);
        assert_eq!(result, replay);
        for (iteration, recorded) in result.incumbent_trace.iter().enumerate() {
            let n = 3 + iteration * c.q;
            let model = fit(&result.x[..n], &result.observations[..n], &c);
            let expected = recommend(&model, &result.x[..n], c.prior_mean);
            assert_eq!(*recorded, expected);
            assert!(recorded.observation_index < n);
        }
    }

    #[test]
    fn recommendation_does_not_select_the_noisiest_raw_minimum() {
        let mut c = config();
        c.prior_mean = 2.0;
        c.kernel.lengthscales = vec![1e-6];
        let mut count = 0;
        let result = minimize_noisy(&mut |_| {
            count += 1;
            if count == 1 {
                NoisyObservation { value: -20.0, noise_variance: 10000.0 }
            } else {
                NoisyObservation { value: 1.0, noise_variance: 1e-6 }
            }
        }, 1, 2, 0, &c);
        let best = result.incumbent_trace[0];
        assert_eq!(best.observation_index, 1);
        assert!((best.mean - (2.0 - 1.0 / 1.000001)).abs() < 1e-10);
        assert!((best.variance - (1.0 - 1.0 / 1.000001)).abs() < 1e-10);
        assert_eq!(result.observations[0].value, -20.0);
    }

    #[test]
    fn zero_variance_is_retained_without_an_invented_noise_floor() {
        let c = config();
        let result = minimize_noisy(&mut |x| NoisyObservation { value: x[0], noise_variance: 0.0 }, 1, 2, 0, &c);
        assert!(result.observations.iter().all(|o| o.noise_variance == 0.0));
        assert!(result.incumbent_trace[0].variance < 1e-12);
    }

    #[test]
    fn invalid_configuration_refuses_before_physics_callbacks() {
        let mut cases = Vec::new();
        let c = config();
        cases.push((0, 2, 1, c.clone()));
        cases.push((1, 0, 1, c.clone()));
        let mut bad = c.clone(); bad.q = 0; cases.push((1, 2, 1, bad));
        let mut bad = c.clone(); bad.mc_samples = 0; cases.push((1, 2, 1, bad));
        let mut bad = c.clone(); bad.acq_starts = 0; cases.push((1, 2, 1, bad));
        let mut bad = c.clone(); bad.acq_evals = 0; cases.push((1, 2, 1, bad));
        let mut bad = c.clone(); bad.prior_mean = f64::NAN; cases.push((1, 2, 1, bad));
        let mut bad = c.clone(); bad.bounds = (-f64::MAX, f64::MAX); cases.push((1, 2, 1, bad));
        let mut bad = c.clone(); bad.kernel.signal = -1.0; cases.push((1, 2, 1, bad));
        let mut bad = c.clone(); bad.kernel.lengthscales = vec![0.0]; cases.push((1, 2, 1, bad));
        cases.push((1, 2, usize::MAX, c));
        for (dim, initial, iters, policy) in cases {
            let mut calls = 0;
            let refused = catch_unwind(AssertUnwindSafe(|| minimize_noisy(
                &mut |x| { calls += 1; observation(x) }, dim, initial, iters, &policy)));
            assert!(refused.is_err());
            assert_eq!(calls, 0);
        }
    }

    #[test]
    fn invalid_observation_refuses_at_the_offending_callback() {
        for bad in [NoisyObservation { value: f64::NAN, noise_variance: 1.0 },
            NoisyObservation { value: 0.0, noise_variance: -0.1 },
            NoisyObservation { value: 0.0, noise_variance: f64::INFINITY }] {
            let mut calls = 0;
            let refused = catch_unwind(AssertUnwindSafe(|| minimize_noisy(&mut |_| {
                calls += 1;
                bad
            }, 1, 2, 1, &config())));
            assert!(refused.is_err());
            assert_eq!(calls, 1);
        }
    }

    #[test]
    fn single_candidate_path_uses_noisy_not_plugin_improvement() {
        let mut c = config(); c.q = 1;
        let report = minimize_noisy(&mut observation, 1, 2, 1, &c);
        let initial_model = fit(&report.x[..2], &report.observations[..2], &c);
        let bank = joint_normal_bank(c.mc_samples, 3, c.seed ^ 0x4E45_4942);
        let gain = q_noisy_expected_improvement(&initial_model, &report.x[..2], &report.x[2..], &bank);
        assert!(gain > 0.0, "the selected candidate must improve the joint acquisition");
        assert_eq!(report.x.len(), 3);
        assert_eq!(report.incumbent_trace.len(), 2);
    }
}
