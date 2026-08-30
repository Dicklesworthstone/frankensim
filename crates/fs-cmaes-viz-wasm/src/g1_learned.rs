//! Transformer adapter for the learned-G1-policy hook (feature `g1-learned`).
//! The trait, trace type, and observation flattening live in `g1_walking.rs`.

// ─── fs-g1-train adapter (feature-gated dependency) ───

use crate::g1_walking::{G1_LEARNED_OBS_DIMS, EpisodeTrace, LearnedG1Policy};
use fs_g1_train::ppo::{gaussian_action, log_gaussian_prob, PolicyLogStd};
use fs_mbd::robot_models::G1_POLICY_ACTUATORS;
use fs_g1_train::transformer::GaitTransformer;

/// Transformer policy adapter: stepwise inference over the G1 rollout with
/// Gaussian exploration, recording (obs, action, log-prob, value) for PPO.
pub struct TransformerG1Policy {
    pub model: GaitTransformer,
    pub log_std: PolicyLogStd,
    pub rng: u64,
    collected_obs: Vec<Vec<f32>>,
    collected_actions: Vec<Vec<f32>>,
    collected_log_probs: Vec<f32>,
    collected_values: Vec<f32>,
}

impl TransformerG1Policy {
    pub fn new(model: GaitTransformer, log_std: PolicyLogStd, rng: u64) -> Self {
        Self {
            model,
            log_std,
            rng,
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
    fn act(&mut self, obs: &[f32; G1_LEARNED_OBS_DIMS], step: usize) -> [f64; G1_POLICY_ACTUATORS] {
        let (mean, value) = self.model.forward_step(obs, step);
        let action = gaussian_action(&mean, &self.log_std.log_std, &mut self.rng);
        let (log_prob, _, _) = log_gaussian_prob(&mean, &self.log_std.log_std, &action);
        self.collected_obs.push(obs.to_vec());
        self.collected_actions.push(action.clone());
        self.collected_log_probs.push(log_prob);
        self.collected_values.push(value);
        let mut out = [0.0f64; G1_POLICY_ACTUATORS];
        for (o, a) in out.iter_mut().zip(action.iter()) {
            *o = *a as f64;
        }
        out
    }
}
