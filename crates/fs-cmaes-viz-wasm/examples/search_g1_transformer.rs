//! Search the transformer's residual head against the REAL G1 objective.
//!
//! Two earlier attempts and what each taught:
//!
//! 1. The browser trained a head in a kinematic stand-in whose forward speed is
//!    a formula with a hard 7.80 m cap. Every number there is about a toy.
//!
//! 2. PPO on the real rollout (examples/train_g1_real.rs) optimises a shaping
//!    reward, and that reward is NOT the owner's verdict. Measured: its batch
//!    reward stayed around 0.13-0.23 while the owner's objective went from
//!    -59 to over 500,000 and the robot fell at step 202 instead of surviving
//!    720. Optimising the wrong scalar destroyed a working controller.
//!
//! So this searches the one thing that matters: the owner's own multi-factor
//! walking objective, the same scalar CMA-ES minimises on the flagship, on the
//! same articulated-body physics. The transformer contributes a residual delta
//! on top of the tuned 5,040-D controller and its head starts at zero, so the
//! search begins at exact parity with that controller and every accepted step
//! is a measured improvement over it.
//!
//! Only the policy head is searched (n_outputs x d_model). The trunk is frozen,
//! which keeps the search small enough to run in minutes on a CPU — the model
//! is tiny, and native Rust needs no GPU for this.
//!
//!   cargo run --release --features g1-learned --example search_g1_transformer
//!
//! Overrides: ITERATIONS, PAIRS, SIGMA, LR, DURATION_S, CHALLENGE=flat|terrain.

use std::time::Instant;

const WEIGHTS_MAGIC: &[u8; 4] = b"FSGT";
const LAYOUT_VERSION: u32 = 1;

/// Serialise the searched model in the same FSGT layout the browser loader
/// reads, so the result of a run is a file rather than a console number that
/// disappears when the process exits.
fn dump_weights(model: &GaitTransformer) -> Vec<u8> {
    let cfg = &model.cfg;
    let mut out = Vec::new();
    out.extend_from_slice(WEIGHTS_MAGIC);
    out.extend_from_slice(&LAYOUT_VERSION.to_le_bytes());
    for dim in [
        cfg.d_model,
        cfg.n_heads,
        cfg.head_dim,
        cfg.n_kv_heads,
        cfg.kv_dim,
        cfg.n_layers,
        cfg.mlp_hidden,
        cfg.context,
        cfg.n_inputs,
        cfg.n_outputs,
    ] {
        out.extend_from_slice(&(dim as u32).to_le_bytes());
    }
    let mut arrays: Vec<&[f32]> = Vec::new();
    arrays.push(&model.embed.params[..]);
    for layer in &model.layers {
        arrays.push(&layer.wq.weights);
        arrays.push(&layer.wk.weights);
        arrays.push(&layer.wv.weights);
        arrays.push(&layer.wo.weights);
        arrays.push(&layer.w_gate.weights);
        arrays.push(&layer.w_up.weights);
        arrays.push(&layer.w_down.weights);
        arrays.push(&layer.norm1.params);
        arrays.push(&layer.norm2.params);
    }
    arrays.push(&model.final_norm.params);
    arrays.push(&model.policy_head.params);
    arrays.push(&model.value_w.params);
    arrays.push(&model.value_b.params);
    // Observations reach the policy already flattened by the walking owner, so
    // this model carries an identity normalisation rather than fitted stats.
    let identity_mean = vec![0.0f32; cfg.n_inputs];
    let identity_var = vec![1.0f32; cfg.n_inputs];
    arrays.push(&identity_mean);
    arrays.push(&identity_var);
    out.extend_from_slice(&(arrays.len() as u32).to_le_bytes());
    for arr in arrays {
        out.extend_from_slice(&(arr.len() as u32).to_le_bytes());
        for value in arr {
            out.extend_from_slice(&value.to_le_bytes());
        }
    }
    out
}

use fs_cmaes_viz_wasm::g1_learned::TransformerG1Policy;
use fs_cmaes_viz_wasm::g1_walking::{
    EpisodeTrace, G1Challenge, G1Task, G1WalkingConfig, G1WalkingEvaluator,
};
use fs_g1_train::ppo::PolicyLogStd;
use fs_g1_train::transformer::{Config, GaitTransformer};

const OBS_DIMS: usize = 42;
const ACT_DIMS: usize = 15;

fn env_usize(key: &str, fallback: usize) -> usize {
    std::env::var(key).ok().and_then(|v| v.parse().ok()).unwrap_or(fallback)
}
fn env_f64(key: &str, fallback: f64) -> f64 {
    std::env::var(key).ok().and_then(|v| v.parse().ok()).unwrap_or(fallback)
}

/// Deterministic standard normal, so a run is reproducible from its seed.
fn normal(state: &mut u64) -> f32 {
    let mut next = || {
        *state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((*state >> 33) as f64 / (1u64 << 31) as f64) as f64
    };
    let u = next().max(1e-12);
    let v = next();
    (( -2.0 * u.ln()).sqrt() * (std::f64::consts::TAU * v).cos()) as f32
}

fn main() {
    let iterations = env_usize("ITERATIONS", 200);
    let pairs = env_usize("PAIRS", 8);
    let sigma0 = env_f64("SIGMA", 0.02) as f32;
    let lr = env_f64("LR", 0.5) as f32;
    let duration_s = env_f64("DURATION_S", 1.5);
    let eval_every = env_usize("EVAL_EVERY", 5);
    let challenge = match std::env::var("CHALLENGE").as_deref() {
        Ok("terrain") => G1Challenge::TerrainAndPush,
        _ => G1Challenge::Flat,
    };

    let mut seed = 0x5EA4_2026_0905_C1A5u64;
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
    policy.set_exploration_enabled(false);

    let evaluator = G1WalkingEvaluator::new(G1WalkingConfig {
        task: G1Task::Walking,
        challenge,
        duration_s,
        ..G1WalkingConfig::default()
    })
    .expect("evaluator builds");

    let dim = policy.model.policy_head.params.len();

    // Score = the owner's objective (lower is better), from the same receipt
    // the site publishes. No shaping, no proxy.
    let mut score = |policy: &mut TransformerG1Policy, head: &[f32]| -> (f64, f64, usize) {
        policy.model.policy_head.params.copy_from_slice(head);
        policy.begin_episode();
        let mut trace = EpisodeTrace::default();
        let receipt = evaluator.rollout_learned(policy, &mut trace).expect("rollout");
        let _ = policy.take_collected();
        (receipt.objective, receipt.distance_m, receipt.completed_steps)
    };

    let mut theta = vec![0.0f32; dim];
    let (base_obj, base_dist, base_steps) = score(&mut policy, &theta);
    println!(
        "baseline (tuned controller, zero head): objective {base_obj:.4}  distance {base_dist:.4} m  steps {base_steps}  head dim {dim}"
    );

    let mut best = theta.clone();
    let mut best_obj = base_obj;
    let mut best_dist = base_dist;
    let mut sigma = sigma0;
    let mut rng = 0x9E37_79B9_7F4A_7C15u64;
    let mut plus = vec![0.0f32; dim];
    let mut minus = vec![0.0f32; dim];
    let mut noise = vec![0.0f32; dim];
    let mut grad = vec![0.0f32; dim];
    let started = Instant::now();
    let mut episodes = 1usize;

    for iteration in 1..=iterations {
        grad.iter_mut().for_each(|g| *g = 0.0);
        let mut scored: Vec<(f64, f32, Vec<f32>)> = Vec::with_capacity(pairs * 2);
        for _ in 0..pairs {
            for value in noise.iter_mut() {
                *value = normal(&mut rng);
            }
            for i in 0..dim {
                plus[i] = theta[i] + sigma * noise[i];
                minus[i] = theta[i] - sigma * noise[i];
            }
            let (fp, _, _) = score(&mut policy, &plus);
            let (fm, _, _) = score(&mut policy, &minus);
            episodes += 2;
            scored.push((fp, 1.0, noise.clone()));
            scored.push((fm, -1.0, noise.clone()));
        }
        // Lower objective is better, so sort ascending and weight the best most.
        scored.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        let n = scored.len() as f32;
        for (rank, (_, sign, draw)) in scored.iter().enumerate() {
            let weight = (n - 1.0 - 2.0 * rank as f32) / (n - 1.0);
            let scale = weight * sign;
            for i in 0..dim {
                grad[i] += scale * draw[i];
            }
        }
        let step = lr * sigma / n;
        for i in 0..dim {
            theta[i] += step * grad[i];
        }

        if iteration % eval_every == 0 || iteration == iterations {
            let (obj, dist, steps) = score(&mut policy, &theta);
            episodes += 1;
            let improved = obj < best_obj;
            if improved {
                best_obj = obj;
                best_dist = dist;
                best.copy_from_slice(&theta);
                sigma = (sigma * 1.05).min(sigma0 * 4.0);
            } else {
                // Elitist restart. The tuned controller sits at a sharp optimum:
                // without this the mean wanders into regions where every sample
                // falls, and nothing pulls it back. Return to the best policy
                // found and search a smaller neighbourhood around it.
                theta.copy_from_slice(&best);
                sigma = (sigma * 0.8).max(sigma0 * 0.02);
            }
            println!(
                "iter {iteration:>4}  objective {obj:>10.4}  distance {dist:.4} m  steps {steps:>4}  best {best_obj:.4}  sigma {sigma:.4}  ({episodes} eps, {:.0}s)",
                started.elapsed().as_secs_f64()
            );
        }
    }

    let (final_obj, final_dist, final_steps) = score(&mut policy, &best);
    println!();
    println!("tuned controller : objective {base_obj:.4}   distance {base_dist:.4} m");
    println!(
        "transformer residual: objective {final_obj:.4}   distance {final_dist:.4} m   steps {final_steps}"
    );
    println!(
        "improvement: {:.4} objective ({:.1}%), {:+.4} m distance, over {episodes} episodes",
        base_obj - final_obj,
        100.0 * (base_obj - final_obj) / base_obj.abs(),
        final_dist - base_dist,
    );
    let l2 = policy.policy_head_l2_norm();
    println!("residual head L2 norm: {l2:.6} (zero would mean it never moved)");

    // Only write when the search actually beat the tuned controller. A file
    // named for a trained policy that is worse than the thing it started from
    // is worse than no file.
    if final_obj < base_obj {
        let bytes = dump_weights(&policy.model);
        let out_dir = std::env::var("OUT_DIR_G1")
            .unwrap_or_else(|_| "target/g1-residual".to_string());
        std::fs::create_dir_all(&out_dir).expect("create output dir");
        let weights_path = format!("{out_dir}/g1-residual-head.bin");
        std::fs::write(&weights_path, &bytes).expect("write weights");
        let receipt = format!(
            "{{\n  \"searchedParams\": {dim},\n  \"episodes\": {episodes},\n  \"challenge\": \"{challenge:?}\",\n  \"durationSeconds\": {duration_s},\n  \"tunedControllerObjective\": {base_obj},\n  \"tunedControllerDistanceMeters\": {base_dist},\n  \"residualObjective\": {final_obj},\n  \"residualDistanceMeters\": {final_dist},\n  \"residualCompletedSteps\": {final_steps},\n  \"headL2Norm\": {l2},\n  \"wallclockSeconds\": {:.1}\n}}\n",
            started.elapsed().as_secs_f64()
        );
        std::fs::write(format!("{out_dir}/search-receipt.json"), receipt)
            .expect("write receipt");
        println!("wrote {weights_path} ({} bytes) and search-receipt.json", bytes.len());
    } else {
        println!("no improvement over the tuned controller; nothing written");
    }
}
