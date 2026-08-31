//! Train the real GaitTransformer policy on the disclosed action-causal
//! stand-in env and
//! export everything the browser ablation needs:
//!   1. `g1-ablation-weights-v1.bin` — length-prefixed little-endian f32 dump
//!      (layout documented in `dump_weights`; mirrored in
//!      cmaes_explainer `app/lib/gaitTransformer.ts`).
//!   2. `g1-ablation-train-receipt.json` — REAL training provenance: config,
//!      parameter count, samples consumed, wallclock, learning-curve summary,
//!      greedy-eval metrics at the training horizon (240 steps) and the
//!      flagship horizon (720 steps).
//!   3. `g1-ablation-golden.json` — deterministic forward outputs from the
//!      RELOADED weight file, consumed by the TS parity test.
//!
//! Run: cargo run --release --example train_ablation -- <out_dir> [iters]
//! Everything is deterministic given the fixed seed chain.

use fs_g1_train::ppo::{
    collect_trajectory, ppo_update, G1Env, PolicyLogStd, PpoConfig, RunningNorm,
};
use fs_g1_train::standin_env::{StandinEnv, STANDIN_CONTRACT_ID};
use fs_g1_train::transformer::{Config, GaitTransformer};
use fs_g1_train::MuonParam;

use std::fs;
use std::time::Instant;

const WEIGHTS_MAGIC: &[u8; 4] = b"FSGT";
const LAYOUT_VERSION: u32 = 1;

fn f32_le_bytes(v: &[f32]) -> Vec<u8> {
    let mut out = Vec::with_capacity(v.len() * 4);
    for x in v {
        out.extend_from_slice(&x.to_le_bytes());
    }
    out
}

/// Dump the full parameter set in the documented order. Every array is
/// length-prefixed (u32 LE) so the TS loader can fail closed on any mismatch.
fn dump_weights(model: &GaitTransformer, norm: &RunningNorm) -> Vec<u8> {
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
    let mut arrays: Vec<&[f32]> = Vec::with_capacity(2 + cfg.n_layers * 9 + 5);
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
    arrays.push(&norm.mean);
    arrays.push(&norm.var);
    out.extend_from_slice(&(arrays.len() as u32).to_le_bytes());
    for arr in arrays {
        out.extend_from_slice(&(arr.len() as u32).to_le_bytes());
        out.extend_from_slice(&f32_le_bytes(arr));
    }
    out
}

/// Rebuild a model's weights from a dump (proves the file round-trips; the
/// golden vectors and evaluations then come from the FILE, not live memory).
fn load_weights(model: &mut GaitTransformer, norm: &mut RunningNorm, bytes: &[u8]) {
    assert_eq!(&bytes[0..4], WEIGHTS_MAGIC, "bad magic");
    let mut pos = 4usize;
    let read_u32 = |pos: &mut usize| -> u32 {
        let v = u32::from_le_bytes(bytes[*pos..*pos + 4].try_into().unwrap());
        *pos += 4;
        v
    };
    let read_f32s = |pos: &mut usize| -> Vec<f32> {
        let len = read_u32(pos) as usize;
        let mut v = Vec::with_capacity(len);
        for i in 0..len {
            let off = *pos + i * 4;
            v.push(f32::from_le_bytes(bytes[off..off + 4].try_into().unwrap()));
        }
        *pos += len * 4;
        v
    };

    assert_eq!(read_u32(&mut pos), LAYOUT_VERSION, "layout version");
    let cfg = &model.cfg;
    let dims: Vec<usize> = (0..10).map(|_| read_u32(&mut pos) as usize).collect();
    assert_eq!(
        dims,
        vec![
            cfg.d_model,
            cfg.n_heads,
            cfg.head_dim,
            cfg.n_kv_heads,
            cfg.kv_dim,
            cfg.n_layers,
            cfg.mlp_hidden,
            cfg.context,
            cfg.n_inputs,
            cfg.n_outputs
        ],
        "config dims must match the live model"
    );
    let count = read_u32(&mut pos) as usize;
    let mut arrays: Vec<Vec<f32>> = Vec::with_capacity(count);
    for _ in 0..count {
        arrays.push(read_f32s(&mut pos));
    }
    assert_eq!(pos, bytes.len(), "trailing bytes in weights file");

    let mut it = arrays.into_iter();
    model.embed.params.copy_from_slice(&it.next().unwrap());
    for layer in model.layers.iter_mut() {
        let put = |dst: &mut MuonParam, src: Vec<f32>| {
            assert_eq!(dst.weights.len(), src.len());
            dst.weights.copy_from_slice(&src);
        };
        put(&mut layer.wq, it.next().unwrap());
        put(&mut layer.wk, it.next().unwrap());
        put(&mut layer.wv, it.next().unwrap());
        put(&mut layer.wo, it.next().unwrap());
        put(&mut layer.w_gate, it.next().unwrap());
        put(&mut layer.w_up, it.next().unwrap());
        put(&mut layer.w_down, it.next().unwrap());
        layer.norm1.params.copy_from_slice(&it.next().unwrap());
        layer.norm2.params.copy_from_slice(&it.next().unwrap());
    }
    model.final_norm.params.copy_from_slice(&it.next().unwrap());
    model
        .policy_head
        .params
        .copy_from_slice(&it.next().unwrap());
    model.value_w.params.copy_from_slice(&it.next().unwrap());
    model.value_b.params.copy_from_slice(&it.next().unwrap());
    norm.mean.copy_from_slice(&it.next().unwrap());
    norm.var.copy_from_slice(&it.next().unwrap());
    assert!(it.next().is_none(), "unconsumed arrays");
    model.reset_cache();
}

/// Greedy (mean-action) rollout — the honest deployment-mode evaluation.
fn greedy_eval(
    env: &mut StandinEnv,
    model: &mut GaitTransformer,
    norm: &RunningNorm,
) -> serde_json::Value {
    let horizon = env.max_steps();
    let mut obs = env.reset(0);
    model.reset_cache();
    let mut total_reward = 0.0f64;
    let mut steps = 0usize;
    let mut normalized = vec![0.0f32; model.cfg.n_inputs];
    loop {
        normalized.copy_from_slice(&obs[..model.cfg.n_inputs]);
        norm.normalize(&mut normalized);
        let (mean, _value) = model.forward_step(&normalized, steps);
        let (next_obs, reward, done) = env.step(&mean);
        total_reward += reward as f64;
        steps += 1;
        obs = next_obs;
        if done {
            break;
        }
    }
    let distance = env.cumulative_distance();
    let work_j = env.cumulative_work_joules();
    serde_json::json!({
        "completedSteps": steps,
        "fell": steps < horizon,
        "distanceMeters": distance,
        "averageSpeedMps": distance / (steps as f64 / 60.0),
        "totalReward": total_reward,
        "actuatorWorkJoules": work_j,
        "costOfTransport": if distance > 1e-9 {
            work_j / (45.0 * 9.81 * distance)
        } else {
            0.0
        },
        "survivalRatePercent": (steps as f64 / horizon as f64) * 100.0,
    })
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let out_dir = args
        .get(1)
        .cloned()
        .unwrap_or_else(|| "artifacts".to_string());
    let iterations: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(400);
    let episode_steps = 240usize; // full episodes: no truncated-GAE bias
    fs::create_dir_all(&out_dir).expect("create out dir");

    let mut rng: u64 = 0x5EED_C0FF_EE50_40u64;
    let mut model_seed: u64 = 0xD1CE_5040_2026_0830u64;

    let cfg = Config::default();
    let mut model = GaitTransformer::new(cfg, 7e-4, 0.9, 1e-4, &mut model_seed);
    let norm = RunningNorm::new(cfg.n_inputs);
    // Fall-gated env: sampled exploration noise across 15 joints must keep
    // Σ|action| below ~3.6 or the episode ends in an instant fall, drowning
    // the learning signal in −50 penalties. Start near-deterministic
    // (σ = e^-1.5 ≈ 0.22 per dim) and let entropy shrink/grow from there.
    let mut log_std = PolicyLogStd::new(model.cfg.n_outputs, 3e-3, -1.5);
    let mut ppo_cfg = PpoConfig::default();
    ppo_cfg.entropy_coef = 0.001;
    let param_count = model.param_count();

    // Zero-init the policy output layer (standard PPO/SAC practice): the
    // untrained head would otherwise emit large random actions, and this
    // env's fall gate trips on ANY sizeable Σ|action| — exploration then
    // starts from near-zero motion and the gradient signal stays clean.
    for v in model.policy_head.params.iter_mut() {
        *v = 0.0;
    }
    println!(
        "[train] params {} iters {} horizon {}",
        param_count, iterations, episode_steps
    );
    let t0 = Instant::now();
    let mut samples: usize = 0;
    let mut mean_rewards: Vec<f64> = Vec::with_capacity(iterations);
    let mut best_mean = f64::NEG_INFINITY;
    let mut best_weights: Option<Vec<u8>> = None;
    let mut best_iter = 0usize;

    for it in 0..iterations {
        let mut env = StandinEnv::new(episode_steps);
        let traj = collect_trajectory(
            &mut env,
            &mut model,
            &norm,
            &log_std,
            &mut rng,
            episode_steps,
            it as u64,
        );
        samples += traj.len();
        let mean_r = traj.rewards.iter().sum::<f32>() as f64 / traj.len().max(1) as f64;
        mean_rewards.push(mean_r);
        let (advantages, returns) = traj.compute_gae(0.0, ppo_cfg.gamma, ppo_cfg.gae_lambda);
        let _kl = ppo_update(
            &mut model,
            &norm,
            &mut log_std,
            &traj,
            &advantages,
            &returns,
            &traj.log_probs,
            &ppo_cfg,
        );

        // Publication selection is based on the post-update greedy policy
        // under the action-causal contract. The previous exporter selected
        // iteration zero from the noisy training rollout before performing
        // its update, allowing an all-zero head to become the "best" walker.
        let mut selection_env = StandinEnv::new(episode_steps);
        let selection = greedy_eval(&mut selection_env, &mut model, &norm);
        let greedy_total = selection["totalReward"]
            .as_f64()
            .expect("greedy total reward is finite");
        let greedy_mean = greedy_total / episode_steps as f64;
        let greedy_distance = selection["distanceMeters"]
            .as_f64()
            .expect("greedy distance is finite");
        let greedy_work = selection["actuatorWorkJoules"]
            .as_f64()
            .expect("greedy work is finite");
        let policy_head_l2 = model
            .policy_head
            .params
            .iter()
            .map(|weight| f64::from(*weight).powi(2))
            .sum::<f64>()
            .sqrt();
        let release_candidate =
            policy_head_l2 > 1e-8 && greedy_distance > 1e-6 && greedy_work > 1e-9;
        if release_candidate && greedy_mean > best_mean {
            best_mean = greedy_mean;
            best_iter = it;
            best_weights = Some(dump_weights(&model, &norm));
            // Continuous checkpointing: the operator may stop the run at any
            // time (compute budgets are real). Whatever is on disk is always
            // the best-so-far checkpoint with its receipt — honestly labeled
            // with the completed iteration count.
            fs::write(
                format!("{out_dir}/g1-ablation-weights-v1.bin"),
                best_weights.as_deref().expect("just assigned"),
            )
            .expect("write checkpoint weights");
            fs::write(
                format!("{out_dir}/g1-ablation-checkpoint.json"),
                serde_json::to_string_pretty(&serde_json::json!({
                    "note": "written only for a post-update, nonzero-head greedy policy that causes motion under the action-causal stand-in; final full receipt lands on normal completion",
                    "environmentContract": STANDIN_CONTRACT_ID,
                    "completedIterations": it + 1,
                    "iterationsPlanned": iterations,
                    "samplesConsumed": samples,
                    "bestMeanReward": best_mean,
                    "bestIteration": best_iter,
                    "bestGreedyEvaluation": selection,
                    "policyHeadL2": policy_head_l2,
                    "wallclockSeconds": t0.elapsed().as_secs_f64(),
                    "stoppedEarly": true,
                }))
                .expect("serialize checkpoint"),
            )
            .expect("write checkpoint receipt");
        }
        if it % 10 == 0 || it == iterations - 1 {
            println!(
                "[train] iter {}/{} mean_reward {:.3} (best {:.3} @ {}) samples {} elapsed {:.1}s",
                it + 1,
                iterations,
                mean_r,
                best_mean,
                best_iter,
                samples,
                t0.elapsed().as_secs_f64()
            );
        }
    }
    let wallclock = t0.elapsed().as_secs_f64();

    // Export weights, then RELOAD them into fresh models for evaluation and
    // golden vectors — the shipped file is the single source of truth.
    let bytes = best_weights.expect(
        "no publishable checkpoint: training must produce a nonzero policy head, nonzero greedy action-caused distance, and nonzero actuator work",
    );
    fs::write(format!("{out_dir}/g1-ablation-weights-v1.bin"), &bytes).expect("write weights");

    let rebuild = || {
        let mut s: u64 = 0xD1CE_5040_2026_0830u64;
        let m = GaitTransformer::new(Config::default(), 2e-3, 0.9, 3e-4, &mut s);
        let n = RunningNorm::new(Config::default().n_inputs);
        (m, n)
    };

    let eval_at = |steps: usize| -> serde_json::Value {
        let (mut m, mut n) = rebuild();
        load_weights(&mut m, &mut n, &bytes);
        let mut env = StandinEnv::new(steps);
        greedy_eval(&mut env, &mut m, &n)
    };
    let eval_240 = eval_at(episode_steps);
    let eval_720 = eval_at(720);

    // Golden forward vectors from the reloaded file (deterministic obs chain).
    let (mut golden_model, mut golden_norm) = rebuild();
    load_weights(&mut golden_model, &mut golden_norm, &bytes);
    let seq_len = 70usize;
    let mut state = 0xA5A5_5A5A_1234_ABCDu64;
    let mut next_f32 = move || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (((state >> 33) as f64 / (1u64 << 31) as f64) * 2.0 - 1.0) as f32
    };
    let obs_rows: Vec<Vec<f32>> = (0..seq_len)
        .map(|_| (0..cfg.n_inputs).map(|_| next_f32()).collect())
        .collect();
    let record_at = [0usize, 1, 63, 64, 69];
    let mut cases = Vec::new();
    for (t, row) in obs_rows.iter().enumerate() {
        let (mean, value) = golden_model.forward_step(row, t);
        if record_at.contains(&t) {
            cases.push(serde_json::json!({
                "position": t,
                "obs": row,
                "meanAction": mean,
                "value": value,
            }));
        }
    }
    let golden = serde_json::json!({
        "layoutVersion": LAYOUT_VERSION,
        "note": "forward outputs of the RELOADED g1-ablation-weights-v1.bin; TS must match within 1e-3 abs (f32-native vs f64-accumulated JS)",
        "cases": cases,
    });
    fs::write(
        format!("{out_dir}/g1-ablation-golden.json"),
        serde_json::to_string_pretty(&golden).unwrap(),
    )
    .expect("write golden");

    let avg = |xs: &[f64]| xs.iter().sum::<f64>() / xs.len().max(1) as f64;
    let receipt = serde_json::json!({
        "crate": "fs-g1-train",
        "example": "train_ablation",
        "layoutVersion": LAYOUT_VERSION,
        "architecture": {
            "dModel": cfg.d_model, "nHeads": cfg.n_heads, "headDim": cfg.head_dim,
            "nKvHeads": cfg.n_kv_heads, "nLayers": cfg.n_layers,
            "mlpHidden": cfg.mlp_hidden, "context": cfg.context,
            "nInputs": cfg.n_inputs, "nOutputs": cfg.n_outputs,
            "parameterCount": param_count,
        },
        "training": {
            "algorithm": "PPO (clipped) + GAE; Muon on hidden 2-D weights; Adam on embed/heads/norms/log_std",
            "environment": "StandinEnv — action-causal explanatory port of cmaes_explainer app/lib/g1StepwiseEnv.ts, full 240-step episodes; not the owner SE(3) rollout",
            "environmentContract": STANDIN_CONTRACT_ID,
            "iterations": iterations,
            "episodeSteps": episode_steps,
            "samplesConsumed": samples,
            "wallclockSeconds": wallclock,
            "meanRewardFirst10": avg(&mean_rewards[..10.min(mean_rewards.len())]),
            "meanRewardLast10": avg(&mean_rewards[mean_rewards.len().saturating_sub(10)..]),
            "bestMeanReward": best_mean,
            "bestIteration": best_iter,
            "shippedCheckpoint": "best post-update greedy mean reward among candidates with nonzero policy head, action-caused distance, and actuator work (snapshot reloaded for eval/golden)",
            "seedChain": { "rolloutRng": "0x5EEDC0FFEE5040", "modelInit": "0xD1CE504020260830" },
        },
        "evaluation": { "greedy240": eval_240, "greedy720": eval_720 },
        "hostNote": "Trained natively (Apple M4, single thread); wallclock is host-specific and disclosed as such.",
    });
    fs::write(
        format!("{out_dir}/g1-ablation-train-receipt.json"),
        serde_json::to_string_pretty(&receipt).unwrap(),
    )
    .expect("write receipt");

    println!(
        "[train] done: {} samples in {:.1}s → {}",
        samples, wallclock, out_dir
    );
}
