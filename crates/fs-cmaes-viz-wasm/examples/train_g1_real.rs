//! PPO training of the transformer against the REAL G1 walking rollout.
//!
//! The shipped artifact was trained against a kinematic stand-in whose forward
//! speed is a formula, not physics, and which caps distance at the commanded
//! speed. Everything measured there is measured against a toy. This trains the
//! same architecture on the owner the flagship actually uses: articulated-body
//! dynamics, real contact, the same receipts the site reports.
//!
//! The transformer is a RESIDUAL on the tuned 5,040-D curriculum policy, and
//! its head starts at zero, so training begins at exactly the tuned
//! controller's behaviour and PPO refines from there rather than relearning
//! locomotion from scratch. That is what the artifact was built for; the
//! browser had been running it standalone, where a zero head means zero action.
//!
//! The existing seam test uses one 0.5 s episode per PPO update, which is a
//! batch of one and cannot learn — its reward is flat across every iteration.
//! This collects a real batch of full-length episodes per update.
//!
//!   cargo run --release --features g1-learned --example train_g1_real
//!
//! Environment overrides: ITERATIONS, EPISODES_PER_ITER, DURATION_S, LR_MUON.

use std::time::Instant;

use fs_cmaes_viz_wasm::g1_learned::TransformerG1Policy;
use fs_cmaes_viz_wasm::g1_walking::{
    EpisodeTrace, G1Challenge, G1Task, G1WalkingConfig, G1WalkingEvaluator,
};
use fs_g1_train::ppo::{ppo_update, PolicyLogStd, PpoConfig, RunningNorm, Trajectory};
use fs_g1_train::transformer::{Config, GaitTransformer};

const OBS_DIMS: usize = 42;
const ACT_DIMS: usize = 15;

fn env_usize(key: &str, fallback: usize) -> usize {
    std::env::var(key).ok().and_then(|v| v.parse().ok()).unwrap_or(fallback)
}

fn env_f64(key: &str, fallback: f64) -> f64 {
    std::env::var(key).ok().and_then(|v| v.parse().ok()).unwrap_or(fallback)
}

fn main() {
    let iterations = env_usize("ITERATIONS", 60);
    let episodes_per_iter = env_usize("EPISODES_PER_ITER", 8);
    let duration_s = env_f64("DURATION_S", 1.5);
    let lr_muon = env_f64("LR_MUON", 1e-3) as f32;

    let mut seed = 0x51D0_2026_0905_A11Cu64;
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
    let log_std = PolicyLogStd::new(ACT_DIMS, 1e-3, -2.3);
    let curriculum = fs_cmaes_viz_wasm::g1_walking::g1_walking_curriculum_mean().to_vec();
    let mut policy = TransformerG1Policy::new(model, log_std, seed, curriculum);

    let walk_cfg = G1WalkingConfig {
        task: G1Task::Walking,
        challenge: G1Challenge::Flat,
        duration_s,
        ..G1WalkingConfig::default()
    };
    let steps_per_episode = (duration_s * 480.0).round() as usize;
    let evaluator = G1WalkingEvaluator::new(walk_cfg).expect("evaluator builds");

    let ppo = PpoConfig {
        epochs_per_batch: 4,
        lr_muon,
        ..PpoConfig::default()
    };
    let norm = RunningNorm::new(OBS_DIMS);

    // Deterministic measurement: exploration off, so the number reported is the
    // policy's own behaviour rather than a lucky sample.
    let mut measure = |policy: &mut TransformerG1Policy| {
        policy.set_exploration_enabled(false);
        policy.begin_episode();
        let mut trace = EpisodeTrace::default();
        let receipt = evaluator.rollout_learned(policy, &mut trace).expect("rollout");
        let _ = policy.take_collected();
        policy.set_exploration_enabled(true);
        receipt
    };

    let baseline = measure(&mut policy);
    println!(
        "baseline (tuned controller, zero head): distance {:.4} m  objective {:.3}  steps {}",
        baseline.distance_m, baseline.objective, baseline.completed_steps
    );

    // Selection tracks the OWNER's multi-factor walking objective — the same
    // number CMA-ES minimises on the flagship — not PPO's shaping reward. They
    // are different quantities, and improving the shaping reward while the
    // owner's verdict worsens is exactly the confusion this project keeps
    // having to unpick.
    let mut best_objective = baseline.objective;
    let mut best_distance = baseline.distance_m;
    let started = Instant::now();

    for iteration in 1..=iterations {
        // Collect a real batch: one PPO update over many full episodes, rather
        // than the single 0.5 s episode the seam test uses.
        let mut traj = Trajectory::new();
        let mut rewards: Vec<f32> = Vec::new();
        let mut episode_reward_sum = 0.0f32;
        for _ in 0..episodes_per_iter {
            policy.begin_episode();
            let mut trace = EpisodeTrace::default();
            let _ = evaluator.rollout_learned(&mut policy, &mut trace).expect("rollout");
            let (obs, actions, log_probs, values) = policy.take_collected();
            let taken = obs.len().min(trace.rewards.len());
            for i in 0..taken {
                traj.observations.push(obs[i].clone());
                traj.actions.push(actions[i].clone());
                traj.log_probs.push(log_probs[i]);
                traj.values.push(values[i]);
                rewards.push(trace.rewards[i]);
                // Terminal on the episode's last transition, so GAE does not
                // bootstrap a return across the boundary into the next episode.
                traj.dones.push(i + 1 == taken);
            }
            episode_reward_sum +=
                trace.rewards.iter().sum::<f32>() / steps_per_episode as f32;
        }
        traj.rewards = rewards;

        let (advantages, returns) = traj.compute_gae(0.0, ppo.gamma, ppo.gae_lambda);
        let old_log_probs = traj.log_probs.clone();
        let (model, log_std) = policy.training_parts();
        let kl = ppo_update(
            model,
            &norm,
            log_std,
            &traj,
            &advantages,
            &returns,
            &old_log_probs,
            &ppo,
        );

        if iteration % 5 == 0 || iteration == iterations {
            let receipt = measure(&mut policy);
            let improved = receipt.objective < best_objective;
            if improved {
                best_objective = receipt.objective;
                best_distance = receipt.distance_m;
            }
            println!(
                "iter {iteration:>4}  distance {:.4} m  objective {:.3}  steps {:>4}  kl {kl:.4}  batch-reward {:.4}  {}  ({:.0}s)",
                receipt.distance_m,
                receipt.objective,
                receipt.completed_steps,
                episode_reward_sum / episodes_per_iter as f32,
                if improved { "*" } else { " " },
                started.elapsed().as_secs_f64(),
            );
        }
    }

    println!();
    println!(
        "baseline objective {:.3} (distance {:.4} m)",
        baseline.objective, baseline.distance_m
    );
    println!(
        "best     objective {:.3} (distance {:.4} m) after {} episodes",
        best_objective,
        best_distance,
        iterations * episodes_per_iter
    );
}
