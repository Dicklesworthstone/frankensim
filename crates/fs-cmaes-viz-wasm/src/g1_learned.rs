//! Transformer adapter for the learned-G1-policy hook (feature `g1-learned`).
//! The trait, trace type, and observation flattening live in `g1_walking.rs`.

// ─── fs-g1-train adapter (feature-gated dependency) ───

use crate::g1_walking::{flatten_observation, G1_LEARNED_OBS_DIMS, EpisodeTrace, LearnedG1Policy};
use fs_mbd::robot_models::{G1PolicyObservation, G1ResidualPolicy};
use fs_g1_train::ppo::{gaussian_action, log_gaussian_prob, PolicyLogStd};
use fs_mbd::robot_models::G1_POLICY_ACTUATORS;
use fs_g1_train::transformer::GaitTransformer;

/// Transformer policy adapter: stepwise inference over the G1 rollout with
/// Gaussian exploration, recording (obs, action, log-prob, value) for PPO.
pub struct TransformerG1Policy {
    pub model: GaitTransformer,
    pub log_std: PolicyLogStd,
    pub rng: u64,
    /// The tuned composed policy (5040-D catalog parameter vector) whose
    /// per-step residual output is the base the learned delta refines.
    /// Zero-init of the policy head makes the adapter START at exactly the
    /// tuned controller's behavior — so PPO refines rather than replaces.
    pub curriculum_params: Vec<f64>,
    collected_obs: Vec<Vec<f32>>,
    collected_actions: Vec<Vec<f32>>,
    collected_log_probs: Vec<f32>,
    collected_values: Vec<f32>,
}

impl TransformerG1Policy {
    pub fn new(
        model: GaitTransformer,
        log_std: PolicyLogStd,
        rng: u64,
        curriculum_params: Vec<f64>,
    ) -> Self {
        let mut model = model;
        // Residual-on-base formulation: delta head starts at zero, so the
        // executed residual equals `base_residual` before any training.
        for w in model.policy_head.params.iter_mut() {
            *w = 0.0;
        }
        Self {
            model,
            log_std,
            rng,
            curriculum_params,
            collected_obs: Vec::new(),
            collected_actions: Vec::new(),
            collected_log_probs: Vec::new(),
            collected_values: Vec::new(),
        }
    }

    pub fn model_mut(&mut self) -> &mut GaitTransformer {
        &mut self.model
    }

    /// Split borrow: (model, log_std) for the PPO update.
    pub fn training_parts(&mut self) -> (&mut GaitTransformer, &mut PolicyLogStd) {
        (&mut self.model, &mut self.log_std)
    }

    /// RMS of the exploration std (training-health diagnostic).
    #[must_use]
    pub fn log_std_std(&self) -> f32 {
        let n = self.log_std.log_std.len() as f32;
        let mean = self.log_std.log_std.iter().sum::<f32>() / n;
        (self.log_std.log_std.iter().map(|l| (l - mean) * (l - mean)).sum::<f32>() / n).sqrt()
    }



    pub fn log_std_mut(&mut self) -> &mut PolicyLogStd {
        &mut self.log_std
    }

    /// Reset the inference cache and hand the collected PPO transitions
    /// back to the trainer. Call once per rollout.
    pub fn begin_episode(&mut self) {
        self.model.reset_cache();
        self.collected_obs.clear();
        self.collected_actions.clear();
        self.collected_log_probs.clear();
        self.collected_values.clear();
    }

    /// (observations, actions, log-probs, values) — the PPO batch.
    #[must_use]
    pub fn take_collected(
        &mut self,
    ) -> (Vec<Vec<f32>>, Vec<Vec<f32>>, Vec<f32>, Vec<f32>) {
        (
            std::mem::take(&mut self.collected_obs),
            std::mem::take(&mut self.collected_actions),
            std::mem::take(&mut self.collected_log_probs),
            std::mem::take(&mut self.collected_values),
        )
    }
}

impl LearnedG1Policy for TransformerG1Policy {
    fn act(&mut self, obs: &G1PolicyObservation, step: usize) -> [f64; G1_POLICY_ACTUATORS] {
        // Base: the tuned composed residual for THIS observation.
        let base = G1ResidualPolicy::new(&self.curriculum_params)
            .expect("curriculum params validate")
            .evaluate(obs)
            .expect("composed residual evaluates");
        // The transformer outputs the DELTA mean; the Gaussian explores
        // deltas; the executed residual is base + delta, clamped to the
        // same [-1, 1] envelope the composed tanh residual lives in. The
        // recorded PPO action is the DELTA (its density under the Gaussian
        // is what ppo_update re-evaluates).
        let flat = flatten_observation(obs);
        let (delta_mean, value) = self.model.forward_step(&flat, step);
        let delta = gaussian_action(&delta_mean, &self.log_std.log_std, &mut self.rng);
        let (log_prob, _, _) = log_gaussian_prob(&delta_mean, &self.log_std.log_std, &delta);
        self.collected_obs.push(flat.to_vec());
        self.collected_actions.push(delta.clone());
        self.collected_log_probs.push(log_prob);
        self.collected_values.push(value);
        let mut out = [0.0f64; G1_POLICY_ACTUATORS];
        for (i, o) in out.iter_mut().enumerate() {
            let executed = base[i] + delta[i] as f64;
            *o = executed.clamp(-1.0, 1.0);
        }
        out
    }
}
