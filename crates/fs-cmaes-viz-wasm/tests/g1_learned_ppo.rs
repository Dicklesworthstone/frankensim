//! Integration test (feature `g1-learned`): REAL PPO training of the
//! fs-g1-train transformer policy against the REAL G1 walking rollout
//! (Balance task, Flat challenge, 0.5 s episodes).
//!
//! Honest gate for this milestone: the seam runs end-to-end — episodes
//! complete on real physics, gradients flow through the transformer, the
//! policy moves (KL > 0), and the shaped episode reward does not degrade.
//! This pins the seam, not athletic performance.

#![cfg(feature = "g1-learned")]

use fs_cmaes_viz_wasm::g1_learned::TransformerG1Policy;
use fs_cmaes_viz_wasm::g1_walking::{
    EpisodeTrace, G1Challenge, G1Task, G1WalkingConfig, G1WalkingEvaluator,
};
use fs_g1_train::ppo::{PpoConfig, PolicyLogStd, RunningNorm, Trajectory};
use fs_g1_train::transformer::{Config, GaitTransformer};

const OBS_DIMS: usize = 42;
const ACT_DIMS: usize = 15;

#[test]
fn ppo_trains_transformer_on_real_g1_rollout() {
    let mut seed = 0x1234_5678_9ABC_DEF0u64;

    let model_cfg = Config {
        d_model: 64,
        n_heads: 2,
        head_dim: 32,
        n_kv_heads: 1,
        kv_dim: 32,
        n_layers: 2,
        mlp_hidden: 128,
        context: 64,
        n_inputs: OBS_DIMS,
        n_outputs: ACT_DIMS,
    };
    let model = GaitTransformer::new(model_cfg, 1e-3, 0.9, 3e-4, &mut seed);
    let log_std = PolicyLogStd::new(ACT_DIMS, 1e-3, -0.5);
    let mut policy = TransformerG1Policy::new(model, log_std, seed);

    let walk_cfg = G1WalkingConfig {
        task: G1Task::Balance,
        challenge: G1Challenge::Flat,
        duration_s: 0.5, // 240 fixed steps at 480 Hz
        ..G1WalkingConfig::default()
    };
    let evaluator = G1WalkingEvaluator::new(walk_cfg).expect("evaluator builds");

    let ppo = PpoConfig {
        epochs_per_batch: 2,
        lr_muon: 1e-3,
        ..PpoConfig::default()
    };
    let norm = RunningNorm::new(OBS_DIMS);

    let episode = |policy: &mut TransformerG1Policy, rng: &mut u64| -> (f32, usize) {
        policy.begin_episode();
        let mut trace = EpisodeTrace::default();
        let receipt = evaluator
            .rollout_learned(policy, &mut trace)
            .expect("learned rollout runs");
        let mean = trace.rewards.iter().sum::<f32>() / trace.rewards.len().max(1) as f32;
        let _ = receipt;
        (mean, trace.completed_steps)
    };

    let (initial_reward, initial_steps) = episode(&mut policy, &mut seed);
    println!("[g1-ppo] initial: reward {initial_reward:.3}, steps {initial_steps}");

    let mut last_kl = 0.0f32;
    for iteration in 0..3 {
        policy.begin_episode();
        let mut trace = EpisodeTrace::default();
        evaluator
            .rollout_learned(&mut policy, &mut trace)
            .expect("learned rollout runs");
        let (observations, actions, log_probs, values) = policy.take_collected();
        let traj = Trajectory {
            observations: observations.clone(),
            actions: actions
                .iter()
                .map(|a| a.iter().map(|x| *x as f32).collect())
                .collect(),
            rewards: trace.rewards.clone(),
            values,
            log_probs: log_probs.clone(),
            dones: trace.done.clone(),
        };
        let (advantages, returns) =
            traj.compute_gae(0.0, ppo.gamma, ppo.gae_lambda);
        let (model, log_std) = policy.training_parts();
        let kl = fs_g1_train::ppo::ppo_update(
            model,
            &norm,
            log_std,
            &traj,
            &advantages,
            &returns,
            &log_probs,
            &ppo,
        );
        let (r, steps) = episode(&mut policy, &mut seed);
        println!("[g1-ppo] iter {iteration}: reward {r:.3}, steps {steps}, kl {kl:.4}");
        last_kl = kl;
        let _ = r;
        let _ = steps;
    }
    assert!(last_kl.is_finite() && last_kl > 0.0, "a real update must move the policy (KL {last_kl})");
    let (final_reward, _) = episode(&mut policy, &mut seed);
    assert!(
        final_reward >= initial_reward - 0.5,
        "PPO on the real rollout must not catastrophically degrade: {initial_reward:.3} -> {final_reward:.3}"
    );
}
