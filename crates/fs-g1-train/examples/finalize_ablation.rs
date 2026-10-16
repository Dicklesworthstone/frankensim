//! Finalize the ablation artifact set from an existing (possibly
//! early-stopped) checkpoint: reload `g1-ablation-weights-v1.bin`, run the
//! greedy evaluations, regenerate the golden vectors, and emit the full
//! training receipt. Run this after stopping `train_ablation` early.
//!
//! Run: cargo run --release --example finalize_ablation -- <out_dir>

use fs_g1_train::ppo::{G1Env, RunningNorm};
use fs_g1_train::standin_env::{StandinEnv, STANDIN_CONTRACT_ID};
use fs_g1_train::transformer::{Config, GaitTransformer};
use fs_g1_train::MuonParam;

use std::fs;

const WEIGHTS_MAGIC: &[u8; 4] = b"FSGT";
const LAYOUT_VERSION: u32 = 1;

/// Rebuild a model's weights from a dump (same contract as the trainer).
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

fn rebuild() -> (GaitTransformer, RunningNorm) {
    let mut s: u64 = 0xD1CE_5040_2026_0830u64;
    let m = GaitTransformer::new(Config::default(), 7e-4, 0.9, 1e-4, &mut s);
    let n = RunningNorm::new(Config::default().n_inputs);
    (m, n)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let out_dir = args
        .get(1)
        .cloned()
        .unwrap_or_else(|| "artifacts".to_string());

    let bytes =
        fs::read(format!("{out_dir}/g1-ablation-weights-v1.bin")).expect("read weights bin");
    let checkpoint: serde_json::Value = serde_json::from_str(
        &fs::read_to_string(format!("{out_dir}/g1-ablation-checkpoint.json"))
            .expect("read checkpoint"),
    )
    .expect("parse checkpoint json");

    let cfg = Config::default();
    let (mut probe, mut probe_norm) = rebuild();
    load_weights(&mut probe, &mut probe_norm, &bytes);
    let param_count = probe.param_count();

    let eval_at = |steps: usize| -> serde_json::Value {
        let (mut m, mut n) = rebuild();
        load_weights(&mut m, &mut n, &bytes);
        let mut env = StandinEnv::new(steps);
        greedy_eval(&mut env, &mut m, &n)
    };
    let episode_steps = 240usize;
    let eval_240 = eval_at(episode_steps);
    let eval_720 = eval_at(720);

    // Refuse to launder the legacy iteration-zero, all-zero-head artifact
    // through the finalizer. Only checkpoints selected by the current
    // action-causal release gate may produce a full receipt and golden set.
    assert_eq!(
        checkpoint["environmentContract"].as_str(),
        Some(STANDIN_CONTRACT_ID),
        "checkpoint predates the action-causal stand-in contract"
    );
    let policy_head_l2 = probe
        .policy_head
        .params
        .iter()
        .map(|weight| f64::from(*weight).powi(2))
        .sum::<f64>()
        .sqrt();
    assert!(policy_head_l2 > 1e-8, "refusing an all-zero policy head");
    assert!(
        eval_240["distanceMeters"].as_f64().unwrap_or(0.0) > 1e-6,
        "refusing a policy with no action-caused greedy distance"
    );
    assert!(
        eval_240["actuatorWorkJoules"].as_f64().unwrap_or(0.0) > 1e-9,
        "refusing a policy with no greedy actuator work"
    );

    // Golden forward vectors from the reloaded file (deterministic obs chain,
    // identical LCG to train_ablation.rs).
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

    let completed = checkpoint["completedIterations"].as_u64().unwrap_or(0) as usize;
    let planned = checkpoint["iterationsPlanned"].as_u64().unwrap_or(400) as usize;
    let receipt = serde_json::json!({
        "crate": "fs-g1-train",
        "example": "train_ablation (finalized by finalize_ablation)",
        "layoutVersion": LAYOUT_VERSION,
        "architecture": {
            "dModel": cfg.d_model, "nHeads": cfg.n_heads, "headDim": cfg.head_dim,
            "nKvHeads": cfg.n_kv_heads, "nLayers": cfg.n_layers,
            "mlpHidden": cfg.mlp_hidden, "context": cfg.context,
            "nInputs": cfg.n_inputs, "nOutputs": cfg.n_outputs,
            "parameterCount": param_count,
        },
        "training": {
            "algorithm": "PPO (clipped) + GAE; Muon on hidden 2-D weights; Adam on embed/heads/norms/log_std; zero-initialized policy head",
            "environment": "StandinEnv — action-causal explanatory port of cmaes_explainer app/lib/g1StepwiseEnv.ts, full 240-step episodes; not the owner SE(3) rollout",
            "environmentContract": STANDIN_CONTRACT_ID,
            "stoppedEarly": true,
            "completedIterations": completed,
            "iterationsPlanned": planned,
            "episodeSteps": episode_steps,
            "samplesConsumed": checkpoint["samplesConsumed"],
            "wallclockSeconds": checkpoint["wallclockSeconds"],
            "bestMeanReward": checkpoint["bestMeanReward"],
            "bestIteration": checkpoint["bestIteration"],
            "shippedCheckpoint": "best post-update greedy mean reward among candidates with nonzero policy head, action-caused distance, and actuator work (snapshot reloaded for eval/golden)",
            "policyHeadL2": policy_head_l2,
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

    println!("[finalize] artifacts completed in {out_dir}");
}
