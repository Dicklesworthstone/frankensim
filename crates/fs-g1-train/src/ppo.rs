//! PPO + GAE training loop over a generic stepwise environment.
//!
//! Design: the env is a trait (`G1Env`) so the kernel can plug in whenever
//! the stepwise API lands (cmaes-j36). For now, tests use a mock env.
//! The transformer policy is the actor; the value estimate is a stub
//! (mean of the policy output) until a separate linear value head lands.
//!
//! HONEST STATUS (2026-08-29 fresh-eyes pass): `ppo_update` evaluates the
//! PPO ratio/KL diagnostic correctly but does NOT yet apply a gradient —
//! hand-rolled backprop through `GaitTransformer` (RMSNorm/GQA/RoPE/
//! SwiGLU) is the pending seam. Callers must treat it as an early-stop
//! KL measurement, not as training. The GAE, rollout, Gaussian policy,
//! and observation-normalization machinery are real.
//!
//! Observation normalization: running mean/var (Welford's method).
//! GAE: λ=0.95, γ=0.99. PPO clip: 0.2. Entropy bonus: configurable.

use crate::transformer::{GaitTransformer, N_INPUTS};

/// A stepwise environment for the G1 walker.
pub trait G1Env {
    /// Reset to the initial state. Returns the first observation.
    fn reset(&mut self, seed: u64) -> Vec<f32>;
    /// Take one action. Returns (next_obs, reward, done).
    fn step(&mut self, action: &[f32]) -> (Vec<f32>, f32, bool);
    /// Observation dimension.
    fn obs_dim(&self) -> usize;
    /// Action dimension.
    fn action_dim(&self) -> usize;
}

/// Running mean/var for observation normalization (Welford's).
pub struct RunningNorm {
    pub mean: Vec<f32>,
    pub var: Vec<f32>,
    pub count: f64,
}

impl RunningNorm {
    pub fn new(dim: usize) -> Self {
        Self { mean: vec![0.0; dim], var: vec![1.0; dim], count: 0.0 }
    }

    pub fn update(&mut self, x: &[f32]) {
        self.count += 1.0;
        let n = self.count as f32;
        for i in 0..x.len() {
            let d = x[i] - self.mean[i];
            self.mean[i] += d / n;
            let d2 = x[i] - self.mean[i];
            // Exact Welford population-variance update: the cross term
            // d*d2 (pre-update x deviation times post-update deviation)
            // carries the (n-1)/n weighting. Using d2^2 here, as the
            // first version did, underestimates the variance.
            self.var[i] += (d * d2 - self.var[i]) / n;
        }
    }

    pub fn normalize(&self, x: &mut [f32]) {
        for i in 0..x.len() {
            x[i] = (x[i] - self.mean[i]) / (self.var[i].sqrt() + 1e-8);
        }
    }
}
/// PPO hyperparameters.
pub struct PpoConfig {
    pub lr: f32,
    pub clip_ratio: f32,
    pub gamma: f32,
    pub gae_lambda: f32,
    pub entropy_coef: f32,
    pub value_coef: f32,
    pub epochs_per_batch: usize,
    pub horizon: usize,
}

impl Default for PpoConfig {
    fn default() -> Self {
        Self {
            lr: 3e-4,
            clip_ratio: 0.2,
            gamma: 0.99,
            gae_lambda: 0.95,
            entropy_coef: 0.01,
            value_coef: 0.5,
            epochs_per_batch: 3,
            horizon: 64,
        }
    }
}

/// A collected trajectory segment.
pub struct Trajectory {
    pub observations: Vec<Vec<f32>>,
    pub actions: Vec<Vec<f32>>,
    pub rewards: Vec<f32>,
    pub values: Vec<f32>,
    pub log_probs: Vec<f32>,
    pub dones: Vec<bool>,
}

impl Trajectory {
    pub fn new() -> Self {
        Self {
            observations: Vec::new(),
            actions: Vec::new(),
            rewards: Vec::new(),
            values: Vec::new(),
            log_probs: Vec::new(),
            dones: Vec::new(),
        }
    }

    pub fn len(&self) -> usize { self.rewards.len() }
    pub fn is_empty(&self) -> bool { self.rewards.is_empty() }

    /// GAE advantage + return computation.
    pub fn compute_gae(&self, last_value: f32, gamma: f32, lambda: f32) -> (Vec<f32>, Vec<f32>) {
        let n = self.rewards.len();
        let mut advantages = vec![0.0f32; n];
        let mut returns = vec![0.0f32; n];
        let mut gae = 0.0f32;
        for t in (0..n).rev() {
            let next_value = if t + 1 < n { self.values[t + 1] } else { last_value };
            let next_non_terminal = if self.dones[t] { 0.0 } else { 1.0 };
            let delta = self.rewards[t] + gamma * next_value * next_non_terminal - self.values[t];
            gae = delta + gamma * lambda * next_non_terminal * gae;
            advantages[t] = gae;
            returns[t] = gae + self.values[t];
        }
        (advantages, returns)
    }
}

/// Simple Gaussian policy head for continuous actions.
/// Returns (action, log_prob) given the mean output.
pub fn gaussian_action(mean: &[f32], log_std: f32, rng: &mut u64) -> (Vec<f32>, f32) {
    let mut actions = Vec::with_capacity(mean.len());
    let mut log_prob = 0.0f32;
    for m in mean {
        // Sample from N(m, exp(log_std)²)
        let noise = gaussian_sample(rng);
        let std = log_std.exp();
        let action = m + noise * std;
        actions.push(action);
        // log π(a|s) = -0.5 * ((a-m)/std)² - log(std) - 0.5*log(2π)
        let u = (action - m) / std;
        log_prob += -0.5 * u * u - log_std - 0.5 * (2.0 * std::f32::consts::PI).ln();
    }
    (actions, log_prob)
}

fn gaussian_sample(state: &mut u64) -> f32 {
    // Box-Muller
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let u1 = ((state.wrapping_mul(0xBF58_476D_1CE4_E5B9) >> 11) as f64 / (1u64 << 53) as f64).max(1e-10) as f32;
    let u2 = ((state.wrapping_mul(0x94D0_49BB_1331_11EB) >> 11) as f64 / (1u64 << 53) as f64) as f32;
    (-2.0 * u1.ln()).sqrt() * (2.0 * std::f32::consts::PI * u2).cos()
}

/// PPO ratio/KL diagnostic over a collected trajectory.
///
/// For each stored (observation, action) pair this re-evaluates the
/// CURRENT policy mean and computes the stored action's Gaussian log
/// density under it — the mathematically correct PPO ratio input. The
/// first version sampled a FRESH action here instead, which made the
/// ratio and the KL diagnostic meaningless.
///
/// HONEST LIMIT: no gradient is applied. The clipped surrogates are
/// computed but hand-rolled backprop through `GaitTransformer`
/// (RMSNorm/GQA/RoPE/SwiGLU) is the pending seam; `model` is passed
/// `&mut` to reserve that call-site. The return value is the mean
/// absolute log-ratio — the standard PPO early-stop signal — and the
/// ONLY effect of this function today. Treat it as a measurement,
/// not as training.
pub fn ppo_update(
    model: &mut GaitTransformer,
    norm: &RunningNorm,
    traj: &Trajectory,
    advantages: &[f32],
    _returns: &[f32],
    old_log_probs: &[f32],
    config: &PpoConfig,
    log_std: &mut f32,
) -> f32 {
    if traj.is_empty() {
        return 0.0;
    }
    // Normalize advantages
    let adv_mean: f32 = advantages.iter().sum::<f32>() / advantages.len() as f32;
    let adv_std = (advantages.iter().map(|a| (a - adv_mean).powi(2)).sum::<f32>() / advantages.len() as f32).sqrt() + 1e-8;

    let mut total_kl = 0.0f32;
    for _epoch in 0..config.epochs_per_batch {
        for t in 0..traj.len() {
            let adv = (advantages[t] - adv_mean) / adv_std;
            // Re-evaluate the policy mean at this observation, with the
            // SAME normalization AND zero-padding path the rollout used
            // (try_into on a shorter obs silently produced an all-zero
            // input here, which is why the KL did not vanish on an
            // unchanged policy).
            let mut obs_norm = traj.observations[t].clone();
            norm.normalize(&mut obs_norm);
            let mut obs_arr = [0.0f32; N_INPUTS];
            let copy_len = N_INPUTS.min(obs_norm.len());
            obs_arr[..copy_len].copy_from_slice(&obs_norm[..copy_len]);
            let output = model.forward(&obs_arr, t);

            // log pi_new(a_t | s_t) for the STORED action under the
            // current Gaussian (no sampling — the density of the
            // collected action), matching the old_log_probs collected
            // at rollout time.
            let std = log_std.exp();
            let mut new_log_prob = 0.0f32;
            for (m, a) in output.iter().zip(traj.actions[t].iter()) {
                let u = (a - m) / std;
                new_log_prob += -0.5 * u * u - *log_std - 0.5 * (2.0 * std::f32::consts::PI).ln();
            }
            let old_lp = old_log_probs[t];

            // PPO clipped surrogate — computed for the diagnostic and
            // the future gradient seam; no optimizer step exists yet.
            let ratio = (new_log_prob - old_lp).exp();
            let clipped = ratio.clamp(1.0 - config.clip_ratio, 1.0 + config.clip_ratio);
            let _surr1 = ratio * adv;
            let _surr2 = clipped * adv;

            total_kl += (new_log_prob - old_lp).abs();
        }
    }
    total_kl / traj.len() as f32
}

/// Collect a trajectory by rolling out the policy in the env.
/// `horizon` caps the steps and `reset_seed` seeds the env (the first
/// version hardcoded 64 and 0).
pub fn collect_trajectory(
    env: &mut dyn G1Env,
    model: &GaitTransformer,
    norm: &RunningNorm,
    log_std: f32,
    rng: &mut u64,
    horizon: usize,
    reset_seed: u64,
) -> Trajectory {
    let mut traj = Trajectory::new();
    let mut obs = env.reset(reset_seed);
    for _ in 0..horizon {
        let mut normalized = obs.clone();
        norm.normalize(&mut normalized);
        let mut obs_arr = [0.0f32; N_INPUTS];
        let copy_len = N_INPUTS.min(normalized.len());
        obs_arr[..copy_len].copy_from_slice(&normalized[..copy_len]);
        let output = model.forward(&obs_arr, traj.len());
        let (action, log_prob) = gaussian_action(&output, log_std, rng);
        let (next_obs, reward, done) = env.step(&action);

        traj.observations.push(obs);
        traj.actions.push(action);
        traj.rewards.push(reward);
        traj.log_probs.push(log_prob);
        traj.dones.push(done);

        // Value estimate stub (mean of the policy output) until a
        // separate linear value head lands.
        let value: f32 = output.iter().sum::<f32>() / output.len() as f32;
        traj.values.push(value);

        obs = next_obs;
        if done { break; }
    }
    traj
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Mock env: 2-D random walk, reward = negative distance from origin.
    struct MockEnv {
        pos: [f32; 2],
        step_count: usize,
    }

    impl G1Env for MockEnv {
        fn reset(&mut self, _seed: u64) -> Vec<f32> {
            self.pos = [0.0, 0.0];
            self.step_count = 0;
            vec![self.pos[0], self.pos[1], 0.0, 0.0]
        }
        fn step(&mut self, action: &[f32]) -> (Vec<f32>, f32, bool) {
            self.pos[0] += action[0].clamp(-1.0, 1.0) * 0.1;
            self.pos[1] += action[1].clamp(-1.0, 1.0) * 0.1;
            self.step_count += 1;
            let reward = -(self.pos[0] * self.pos[0] + self.pos[1] * self.pos[1]);
            (vec![self.pos[0], self.pos[1], 0.0, 0.0], reward, self.step_count >= 32)
        }
        fn obs_dim(&self) -> usize { 4 }
        fn action_dim(&self) -> usize { 2 }
    }

    #[test]
    fn ppo_collects_and_computes_gae() {
        let mut seed = 42u64;
        let mut env = MockEnv { pos: [0.0, 0.0], step_count: 0 };
        let mut seed2 = 0xD1CEu64;
        let model = GaitTransformer::new(&mut seed2);
        let norm = RunningNorm::new(4);
        let traj = collect_trajectory(&mut env, &model, &norm, -0.5, &mut seed, 64, 0);
        assert!(!traj.is_empty());
        assert!(traj.len() <= 32);
        let (advantages, returns) = traj.compute_gae(0.0, 0.99, 0.95);
        assert_eq!(advantages.len(), traj.len());
        assert_eq!(returns.len(), traj.len());
        // Returns should be finite
        for r in returns.iter() {
            assert!(r.is_finite());
        }
    }

    #[test]
    fn ppo_update_kl_diagnostic_and_rollout_stability() {
        // HONEST TEST (renamed from `ppo_reduces_mock_env_loss`, which
        // claimed training that did not happen): `ppo_update` is a
        // ratio/KL DIAGNOSTIC — it must return a finite, non-negative
        // value, leave the model untouched, and the rollout reward must
        // remain stable (no catastrophic degradation). The real
        // learning test moves to the G1 env once transformer backprop
        // lands.
        let mut seed = 123u64;
        let mut env = MockEnv { pos: [0.0, 0.0], step_count: 0 };
        let mut model_seed = 0xD1CEu64;
        let mut model = GaitTransformer::new(&mut model_seed);
        let norm = RunningNorm::new(4);
        let config = PpoConfig::default();
        let mut log_std = -0.5f32;
        let initial_reward: f32 = {
            let traj = collect_trajectory(&mut env, &model, &norm, log_std, &mut seed, 64, 0);
            traj.rewards.iter().sum::<f32>() / traj.len().max(1) as f32
        };
        let mut kl_sum = 0.0f32;
        for _ in 0..3 {
            let traj = collect_trajectory(&mut env, &model, &norm, log_std, &mut seed, 64, 0);
            let (advantages, returns_from_gae) = traj.compute_gae(0.0, config.gamma, config.gae_lambda);
            let old_log_probs = traj.log_probs.clone();
            let kl = ppo_update(&mut model, &norm, &traj, &advantages, &returns_from_gae, &old_log_probs, &config, &mut log_std);
            assert!(kl.is_finite() && kl >= 0.0, "KL diagnostic must be finite and non-negative, got {kl}");
            kl_sum += kl;
            let _ = env.reset(0);
        }
        // Policy is unchanged between rollout and diagnostic, so the
        // stored-action density under the current policy must EQUAL the
        // rollout density bit-for-bit -> exactly 0. This pins the
        // density evaluation; when real gradient application lands,
        // relax to `>= 0` (it becomes strictly positive).
        assert_eq!(kl_sum, 0.0, "unchanged policy must give exactly 0 log-ratio drift");
        let final_reward: f32 = {
            let traj = collect_trajectory(&mut env, &model, &norm, log_std, &mut seed, 64, 0);
            traj.rewards.iter().sum::<f32>() / traj.len().max(1) as f32
        };
        assert!(final_reward >= initial_reward - 1.0, "rollout should not catastrophically degrade: {initial_reward} -> {final_reward}");
    }


    #[test]
    fn running_norm_converges() {
        let mut norm = RunningNorm::new(1);
        for i in 0..100 {
            norm.update(&[i as f32]);
        }
        assert!((norm.mean[0] - 49.5).abs() < 1.0);
        assert!(norm.var[0] > 100.0);
    }

    #[test]
    fn gaussian_action_is_deterministic_per_seed() {
        let mut rng1 = 42u64;
        let mut rng2 = 42u64;
        let mean = [0.0f32, 0.0];
        let (a1, lp1) = gaussian_action(&mean, -0.5, &mut rng1);
        let (a2, lp2) = gaussian_action(&mean, -0.5, &mut rng2);
        assert_eq!(a1, a2);
        assert!((lp1 - lp2).abs() < 1e-6);
    }
}
