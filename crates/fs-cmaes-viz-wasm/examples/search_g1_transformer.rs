//! Search the transformer's parameters against the REAL G1 objective, using
//! this project's own CMA-ES, on every core.
//!
//! Three earlier attempts and what each taught:
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
//! 3. A hand-rolled rank-shaped ES on the 960-parameter head, single threaded,
//!    reached -86.60 against the tuned controller's -59.41. That worked, but it
//!    needed an elitist restart bolted on to stop the mean wandering off a
//!    sharp optimum, it used one core out of ten, and it left 99% of the
//!    network frozen.
//!
//! Backprop is not what was missing, and no GPU is either. You cannot
//! backpropagate through the articulated-body solver, which is exactly why PPO
//! had to invent a differentiable shaping reward and why that reward pointed
//! somewhere else. What was missing was a real optimiser, real parallelism, and
//! permission to move more than the last layer.
//!
//! So this drives `fs_dfo::CmaOptimizer` — the same LM-CMA the flagship runs,
//! built to scale to 5,040-D — over the transformer's parameters, evaluating a
//! whole generation in parallel across the machine's cores. The transformer
//! contributes a residual delta on top of the tuned 5,040-D controller and its
//! head starts at zero, so the search begins at exact parity with that
//! controller and every accepted step is a measured improvement over it.
//!
//!   cargo run --release --features g1-learned --example search_g1_transformer
//!
//! Overrides:
//!   SCOPE=head|block|all   which parameters move (default head)
//!   CHALLENGES=flat|terrain|both   conditions each candidate must satisfy
//!   GENERATIONS, SIGMA, POPULATION, THREADS, DURATION_S, OUT_DIR_G1

use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

use fs_cmaes_viz_wasm::g1_learned::TransformerG1Policy;
use fs_cmaes_viz_wasm::g1_walking::{
    EpisodeTrace, G1Challenge, G1Task, G1WalkingConfig, G1WalkingEvaluator,
};
use fs_dfo::{CmaConfig, CmaFamily, CmaOptimizer};
use fs_g1_train::ppo::PolicyLogStd;
use fs_g1_train::transformer::{Config, GaitTransformer};

const OBS_DIMS: usize = 42;
const ACT_DIMS: usize = 15;
const WEIGHTS_MAGIC: &[u8; 4] = b"FSGT";
const LAYOUT_VERSION: u32 = 1;
const MODEL_SEED: u64 = 0x5EA4_2026_0905_C1A5;

/// Which parameters the search is allowed to move.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Scope {
    /// Output layer only: the readout of features the frozen trunk computes.
    Head,
    /// Final transformer block, its norm, and the head.
    Block,
    /// Every parameter that affects the action.
    All,
}

fn env_usize(key: &str, fallback: usize) -> usize {
    std::env::var(key).ok().and_then(|v| v.parse().ok()).unwrap_or(fallback)
}
fn env_f64(key: &str, fallback: f64) -> f64 {
    std::env::var(key).ok().and_then(|v| v.parse().ok()).unwrap_or(fallback)
}

fn model_config() -> Config {
    Config {
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
    }
}

/// Build the policy. Deterministic in `MODEL_SEED`, so every worker thread
/// constructs a bit-identical starting model without sharing one.
fn build_policy() -> TransformerG1Policy {
    let mut seed = MODEL_SEED;
    let model = GaitTransformer::new(model_config(), 1e-3, 0.9, 3e-4, &mut seed);
    let log_std = PolicyLogStd::new(ACT_DIMS, 1e-3, -2.3);
    let curriculum = fs_cmaes_viz_wasm::g1_walking::g1_walking_curriculum_mean().to_vec();
    let mut policy = TransformerG1Policy::new(model, log_std, seed, curriculum);
    policy.set_exploration_enabled(false);
    policy
}

/// The parameter slices the given scope may move, in a fixed order so that a
/// flat search vector means the same thing on every thread and across runs.
fn searchable(model: &mut GaitTransformer, scope: Scope) -> Vec<&mut [f32]> {
    let n_layers = model.cfg.n_layers;
    let mut slices: Vec<&mut [f32]> = Vec::new();
    let first_layer = match scope {
        Scope::Head => n_layers,
        Scope::Block => n_layers - 1,
        Scope::All => 0,
    };
    if scope == Scope::All {
        slices.push(&mut model.embed.params[..]);
    }
    for (index, layer) in model.layers.iter_mut().enumerate() {
        if index < first_layer {
            continue;
        }
        slices.push(&mut layer.wq.weights[..]);
        slices.push(&mut layer.wk.weights[..]);
        slices.push(&mut layer.wv.weights[..]);
        slices.push(&mut layer.wo.weights[..]);
        slices.push(&mut layer.w_gate.weights[..]);
        slices.push(&mut layer.w_up.weights[..]);
        slices.push(&mut layer.w_down.weights[..]);
        slices.push(&mut layer.norm1.params[..]);
        slices.push(&mut layer.norm2.params[..]);
    }
    if scope != Scope::Head {
        slices.push(&mut model.final_norm.params[..]);
    }
    slices.push(&mut model.policy_head.params[..]);
    slices
}

fn searchable_len(scope: Scope) -> usize {
    let mut model = build_policy().model;
    searchable(&mut model, scope).iter().map(|s| s.len()).sum()
}

/// Read the searchable parameters into a flat vector.
fn gather(model: &mut GaitTransformer, scope: Scope) -> Vec<f64> {
    let mut out = Vec::new();
    for slice in searchable(model, scope) {
        out.extend(slice.iter().map(|v| f64::from(*v)));
    }
    out
}

/// Write a flat vector back into the searchable parameters.
fn scatter(model: &mut GaitTransformer, scope: Scope, flat: &[f64]) {
    let mut cursor = 0usize;
    for slice in searchable(model, scope) {
        for value in slice.iter_mut() {
            *value = flat[cursor] as f32;
            cursor += 1;
        }
    }
    debug_assert_eq!(cursor, flat.len(), "flat vector length matches the scope");
}

/// One candidate's verdict: the mean of the owner's objective over every
/// condition in the evaluation set. Averaging rather than taking the best
/// stops a candidate buying a flat-ground win with a terrain collapse.
fn score(
    policy: &mut TransformerG1Policy,
    evaluators: &[G1WalkingEvaluator],
    scope: Scope,
    flat: &[f64],
) -> (f64, f64, usize) {
    scatter(&mut policy.model, scope, flat);
    let mut objective = 0.0;
    let mut distance = 0.0;
    let mut steps = usize::MAX;
    for evaluator in evaluators {
        policy.begin_episode();
        let mut trace = EpisodeTrace::default();
        let receipt = evaluator.rollout_learned(policy, &mut trace).expect("rollout");
        let _ = policy.take_collected();
        objective += receipt.objective;
        distance += receipt.distance_m;
        steps = steps.min(receipt.completed_steps);
    }
    let n = evaluators.len() as f64;
    (objective / n, distance / n, steps)
}

fn build_evaluators(challenges: &[G1Challenge], duration_s: f64) -> Vec<G1WalkingEvaluator> {
    challenges
        .iter()
        .map(|challenge| {
            G1WalkingEvaluator::new(G1WalkingConfig {
                task: G1Task::Walking,
                challenge: *challenge,
                duration_s,
                ..G1WalkingConfig::default()
            })
            .expect("evaluator builds")
        })
        .collect()
}

/// Serialise the searched model in the FSGT layout the browser loader reads, so
/// a run leaves a file rather than a console number that dies with the process.
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
    // Observations reach the policy already flattened by the walking owner, so
    // this model carries an identity normalisation rather than fitted stats.
    let identity_mean = vec![0.0f32; cfg.n_inputs];
    let identity_var = vec![1.0f32; cfg.n_inputs];
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

fn main() {
    let scope = match std::env::var("SCOPE").as_deref() {
        Ok("block") => Scope::Block,
        Ok("all") => Scope::All,
        _ => Scope::Head,
    };
    let challenges: Vec<G1Challenge> = match std::env::var("CHALLENGES").as_deref() {
        Ok("terrain") => vec![G1Challenge::TerrainAndPush],
        Ok("both") => vec![G1Challenge::Flat, G1Challenge::TerrainAndPush],
        _ => vec![G1Challenge::Flat],
    };
    let generations = env_usize("GENERATIONS", 300);
    let duration_s = env_f64("DURATION_S", 1.5);
    let sigma = env_f64("SIGMA", 0.002);
    let threads = env_usize(
        "THREADS",
        std::thread::available_parallelism().map_or(1, std::num::NonZeroUsize::get),
    );

    let dim = searchable_len(scope);
    let population = env_usize("POPULATION", 4 + (3.0 * (dim as f64).ln()) as usize);

    // The zero head is exactly the tuned controller, so this baseline is not a
    // separate policy — it is this policy before the search moves it.
    let mut probe = build_policy();
    let probe_evaluators = build_evaluators(&challenges, duration_s);
    let start_point = gather(&mut probe.model, scope);
    let (base_obj, base_dist, base_steps) =
        score(&mut probe, &probe_evaluators, scope, &start_point);
    println!(
        "baseline (tuned controller, zero head): objective {base_obj:.4}  distance {base_dist:.4} m  steps {base_steps}"
    );
    println!(
        "scope {scope:?}  dim {dim}  population {population}  threads {threads}  challenges {}",
        challenges.len()
    );

    let mut optimizer = CmaOptimizer::new(CmaConfig {
        family: CmaFamily::LmCma,
        mean: start_point.clone(),
        sigma,
        max_evaluations: generations * population,
        seed: 0x0905_2026_C1A5_5EA4,
        population_size: Some(population),
        memory: None,
    })
    .expect("optimizer builds");

    let mut best = start_point.clone();
    let mut best_obj = base_obj;
    let mut best_dist = base_dist;
    let mut best_steps = base_steps;
    let mut evaluations = 0usize;
    let started = Instant::now();

    for generation in 1..=generations {
        let Ok(batch) = optimizer.ask() else { break };
        let candidates = batch.candidates().to_vec();

        // Every candidate in a generation is independent, so the whole
        // generation is one parallel batch. Each worker builds its own policy
        // and evaluators once and then pulls candidates until they run out,
        // which keeps every core busy even though falls finish early.
        let cursor = AtomicUsize::new(0);
        let results: Vec<Vec<(usize, f64, f64, usize)>> = std::thread::scope(|s| {
            let cursor = &cursor;
            let candidates = &candidates;
            let challenges = &challenges;
            let handles: Vec<_> = (0..threads)
                .map(|_| {
                    s.spawn(move || {
                        let mut policy = build_policy();
                        let evaluators = build_evaluators(challenges, duration_s);
                        let mut local = Vec::new();
                        loop {
                            let index = cursor.fetch_add(1, Ordering::Relaxed);
                            if index >= candidates.len() {
                                break;
                            }
                            let (objective, distance, steps) =
                                score(&mut policy, &evaluators, scope, &candidates[index]);
                            local.push((index, objective, distance, steps));
                        }
                        local
                    })
                })
                .collect();
            handles.into_iter().map(|h| h.join().expect("worker")).collect()
        });

        let mut objectives = vec![f64::INFINITY; candidates.len()];
        for (index, objective, distance, steps) in results.into_iter().flatten() {
            // A non-finite objective is a candidate that broke the rollout, not
            // a good one; rank it last rather than letting it poison the update.
            objectives[index] = if objective.is_finite() { objective } else { f64::MAX };
            if objectives[index] < best_obj {
                best_obj = objectives[index];
                best_dist = distance;
                best_steps = steps;
                best.copy_from_slice(&candidates[index]);
            }
        }
        evaluations += candidates.len();
        if optimizer.tell(&batch, &objectives).is_err() {
            break;
        }

        if generation % 10 == 0 || generation == generations {
            let elapsed = started.elapsed().as_secs_f64();
            println!(
                "gen {generation:>4}  best {best_obj:>12.4}  distance {best_dist:.4} m  steps {best_steps:>4}  ({evaluations} evals, {elapsed:.0}s, {:.1} evals/s)",
                evaluations as f64 / elapsed.max(1e-9),
            );
        }
    }

    let mut final_policy = build_policy();
    let final_evaluators = build_evaluators(&challenges, duration_s);
    let (final_obj, final_dist, final_steps) =
        score(&mut final_policy, &final_evaluators, scope, &best);
    println!();
    println!("tuned controller    : objective {base_obj:.4}   distance {base_dist:.4} m");
    println!(
        "transformer residual: objective {final_obj:.4}   distance {final_dist:.4} m   steps {final_steps}"
    );
    println!(
        "improvement: {:.4} objective ({:.1}%), {:+.4} m distance, over {evaluations} evaluations in {:.0}s",
        base_obj - final_obj,
        100.0 * (base_obj - final_obj) / base_obj.abs(),
        final_dist - base_dist,
        started.elapsed().as_secs_f64(),
    );

    // Only write when the search actually beat the tuned controller. A file
    // named for a trained policy that is worse than the thing it started from
    // is worse than no file.
    if final_obj < base_obj {
        let bytes = dump_weights(&final_policy.model);
        let out_dir =
            std::env::var("OUT_DIR_G1").unwrap_or_else(|_| "target/g1-residual".to_string());
        std::fs::create_dir_all(&out_dir).expect("create output dir");
        let weights_path = format!("{out_dir}/g1-residual-head.bin");
        std::fs::write(&weights_path, &bytes).expect("write weights");
        let receipt = format!(
            "{{\n  \"optimizer\": \"LM-CMA (fs-dfo)\",\n  \"scope\": \"{scope:?}\",\n  \"searchedParams\": {dim},\n  \"population\": {population},\n  \"evaluations\": {evaluations},\n  \"challenges\": {},\n  \"durationSeconds\": {duration_s},\n  \"tunedControllerObjective\": {base_obj},\n  \"tunedControllerDistanceMeters\": {base_dist},\n  \"residualObjective\": {final_obj},\n  \"residualDistanceMeters\": {final_dist},\n  \"residualCompletedSteps\": {final_steps},\n  \"wallclockSeconds\": {:.1}\n}}\n",
            challenges.len(),
            started.elapsed().as_secs_f64()
        );
        std::fs::write(format!("{out_dir}/search-receipt.json"), receipt).expect("write receipt");
        println!("wrote {weights_path} ({} bytes) and search-receipt.json", bytes.len());
    } else {
        println!("no improvement over the tuned controller; nothing written");
    }
}
