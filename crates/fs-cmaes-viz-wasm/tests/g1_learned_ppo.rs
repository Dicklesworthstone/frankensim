//! Integration test (feature `g1-learned`): REAL PPO training of the
//! fs-g1-train transformer policy against the REAL G1 walking rollout
//! (Balance task, Flat challenge, 0.5 s episodes).
//!
//! Honest gate for this milestone: the seam runs end-to-end — episodes
//! complete on real physics, gradients flow through the transformer, the
//! policy moves (KL > 0), deterministic learned deltas causally change owner
//! receipts, and the shaped episode reward avoids catastrophic regression.
//! This pins the seam and its counterfactuals, not locomotion performance.

#![cfg(feature = "g1-learned")]

use fs_cmaes_viz_wasm::g1_learned::TransformerG1Policy;
use fs_cmaes_viz_wasm::g1_walking::{
    EpisodeTrace, G1Challenge, G1Task, G1WalkingConfig, G1WalkingEvaluator, G1WalkingReceipt,
};
use fs_g1_train::ppo::{PolicyLogStd, PpoConfig, RunningNorm, Trajectory};
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
    // Sigma ~0.1: refine rather than swamp the composed base.
    let log_std = PolicyLogStd::new(ACT_DIMS, 1e-3, -2.3);
    // Balance task -> the stabilizing curriculum is the right composed base.
    let curriculum_params = fs_cmaes_viz_wasm::g1_walking::g1_stabilizing_policy_mean().to_vec();
    let mut policy = TransformerG1Policy::new(model, log_std, seed, curriculum_params);

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

    let walk_cfg_duration_steps = 240usize; // 0.5 s at 480 Hz
    let episode = |policy: &mut TransformerG1Policy| -> (f32, G1WalkingReceipt) {
        policy.begin_episode();
        let mut trace = EpisodeTrace::default();
        let receipt = evaluator
            .rollout_learned(policy, &mut trace)
            .expect("learned rollout runs");
        // Horizon-normalized: a fall at step 30 must score below a full
        // survival, even though the mean over surviving steps looks fine.
        let horizon = walk_cfg_duration_steps;
        let normalized = trace.rewards.iter().sum::<f32>() / horizon as f32;
        (normalized, receipt)
    };

    policy.set_exploration_enabled(false);
    let (initial_reward, initial_receipt) = episode(&mut policy);
    let (_, initial_actions, _, _) = policy.take_collected();
    assert!(
        !initial_actions.is_empty(),
        "the deterministic baseline must exercise the learned-policy adapter"
    );
    assert!(
        initial_actions
            .iter()
            .flatten()
            .all(|action| *action == 0.0),
        "the zero-initialized deterministic delta head must exactly reproduce the composed base policy"
    );
    println!(
        "[g1-ppo] deterministic baseline: reward {initial_reward:.3}, steps {}, objective {:.6}",
        initial_receipt.completed_steps, initial_receipt.objective
    );

    policy.set_exploration_enabled(true);

    let mut last_kl = 0.0f32;
    for iteration in 0..40 {
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
                .map(|a| a.to_vec())
                .collect(),
            rewards: trace.rewards.clone(),
            values,
            log_probs: log_probs.clone(),
            dones: trace.done.clone(),
        };
        let (advantages, returns) = traj.compute_gae(0.0, ppo.gamma, ppo.gae_lambda);
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
        let (r, receipt) = episode(&mut policy);
        let sigma = policy.log_std_std();
        println!(
            "[g1-ppo] iter {iteration}: reward {r:.3}, steps {}, kl {kl:.4}, sigma {sigma:.3}",
            receipt.completed_steps
        );
        last_kl = kl;
    }
    assert!(
        last_kl.is_finite() && last_kl > 0.0,
        "a real update must move the policy (KL {last_kl})"
    );
    let head_l2 = policy.policy_head_l2_norm();
    assert!(
        head_l2.is_finite() && head_l2 > 0.0,
        "the real-owner PPO update must move the zero-initialized policy head (L2 {head_l2})"
    );

    policy.set_exploration_enabled(false);
    let (final_reward, final_receipt) = episode(&mut policy);
    let (_, final_actions, _, _) = policy.take_collected();
    assert!(
        final_actions
            .iter()
            .flatten()
            .any(|action| action.abs() > 1e-8),
        "deterministic learned deltas must reach the owner engine after training"
    );
    let receipt_delta = (final_receipt.objective - initial_receipt.objective).abs()
        + (final_receipt.distance_m - initial_receipt.distance_m).abs()
        + (final_receipt.actuator_work_j - initial_receipt.actuator_work_j).abs();
    assert!(
        receipt_delta > 1e-12,
        "nonzero deterministic learned actions must causally change owner-engine receipts"
    );
    println!(
        "[g1-ppo] deterministic learned: reward {final_reward:.3}, steps {}, objective {:.6}, distance {:.6} m, work {:.6} J, head L2 {head_l2:.6}",
        final_receipt.completed_steps,
        final_receipt.objective,
        final_receipt.distance_m,
        final_receipt.actuator_work_j
    );
    assert!(
        final_reward >= initial_reward - 0.5,
        "PPO on the real rollout must not catastrophically degrade: {initial_reward:.3} -> {final_reward:.3}"
    );
    let sigma_final = policy.log_std_std();
    assert!(
        sigma_final > 1e-3,
        "exploration std collapsed: {sigma_final}"
    );
}
