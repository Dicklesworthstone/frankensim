//! Transformer adapter for the learned-G1-policy hook (feature `g1-learned`).
//! The trait, trace type, and observation flattening live in `g1_walking.rs`.

// ─── fs-g1-train adapter (feature-gated dependency) ───

use crate::g1_walking::{flatten_observation, LearnedG1Policy};
use fs_g1_train::ppo::{gaussian_action, log_gaussian_prob, PolicyLogStd};
use fs_g1_train::transformer::GaitTransformer;
use fs_mbd::robot_models::{G1PolicyObservation, G1ResidualPolicy, G1_POLICY_ACTUATORS};

/// Collected observations, actions, log probabilities, and values for one PPO
/// rollout, in that order.
pub type G1PpoRolloutBatch = (Vec<Vec<f32>>, Vec<Vec<f32>>, Vec<f32>, Vec<f32>);

/// Transformer policy adapter: stepwise inference over the G1 rollout with
/// Gaussian exploration, recording (obs, action, log-prob, value) for PPO.
pub struct TransformerG1Policy {
    pub model: GaitTransformer,
    pub log_std: PolicyLogStd,
    pub rng: u64,
    exploration_enabled: bool,
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
            exploration_enabled: true,
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
        (self
            .log_std
            .log_std
            .iter()
            .map(|l| (l - mean) * (l - mean))
            .sum::<f32>()
            / n)
            .sqrt()
    }

    /// Enable Gaussian sampling for PPO collection or disable it for a
    /// deterministic owner-engine evaluation of the policy-head mean.
    pub fn set_exploration_enabled(&mut self, enabled: bool) {
        self.exploration_enabled = enabled;
    }

    /// L2 norm of the learned residual head. A nonzero value proves the
    /// zero-initialized policy head received an optimizer update.
    #[must_use]
    pub fn policy_head_l2_norm(&self) -> f32 {
        self.model
            .policy_head
            .params
            .iter()
            .map(|weight| weight * weight)
            .sum::<f32>()
            .sqrt()
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
    pub fn take_collected(&mut self) -> G1PpoRolloutBatch {
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
        let delta = if self.exploration_enabled {
            gaussian_action(&delta_mean, &self.log_std.log_std, &mut self.rng)
        } else {
            delta_mean.clone()
        };
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

// ─── Browser-side training against the real walking owner ───

use crate::g1_walking::{
    g1_walking_curriculum_mean, EpisodeTrace, G1Challenge, G1Task, G1WalkingConfig,
    G1WalkingEvaluator,
};
use fs_dfo::{CmaAsk, CmaConfig, CmaFamily, CmaOptimizer};
use fs_g1_train::transformer::Config as TransformerConfig;

const OBS_DIMS: usize = 42;
const ACT_DIMS: usize = 15;
const WEIGHTS_MAGIC: &[u8; 4] = b"FSGT";
const LAYOUT_VERSION: u32 = 1;
/// Fixed so a browser run and the native `search_g1_transformer` example start
/// from a bit-identical frozen trunk and their numbers can be compared.
const MODEL_SEED: u64 = 0x5EA4_2026_0905_C1A5;
/// Generations one LM-CMA run may spend before the trainer restarts it. The
/// owner requires a finite budget, so an open-ended browser session is a chain
/// of finite runs rather than one unbounded one.
const RUN_GENERATION_ALLOWANCE: usize = 4096;

/// The small causal transformer the residual search trains. Deliberately tiny:
/// the whole point is that this needs no GPU.
#[must_use]
pub const fn browser_transformer_config() -> TransformerConfig {
    TransformerConfig {
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

/// Serialise a model in the FSGT layout the browser's TypeScript loader reads.
#[must_use]
pub fn dump_gait_transformer(model: &GaitTransformer) -> Vec<u8> {
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

/// Words in a progress packet, in order: status, generation, evaluations,
/// best objective, best distance, baseline objective, baseline distance,
/// restarts, population, last objective, last distance, last completed steps.
pub const TRAINER_PACKET_WORDS: usize = 12;
/// A candidate was evaluated and the run continues.
pub const TRAINER_STATUS_RUNNING: f64 = 0.0;
/// The generation just closed and the distribution was updated.
pub const TRAINER_STATUS_GENERATION: f64 = 1.0;
/// The optimizer refused to continue; `best_head` still holds the result.
pub const TRAINER_STATUS_STOPPED: f64 = 2.0;

/// LM-CMA over the transformer's policy head against the real G1 walking
/// objective, advanced ONE candidate at a time.
///
/// The native `search_g1_transformer` example evaluates a whole generation in
/// parallel across cores. A browser tab has no cores to spare and must stay
/// responsive, so this exposes the same search at single-candidate
/// granularity: each `pump` runs one 1.5 s rollout, which is a unit small
/// enough to keep a progress bar honest and to interrupt.
///
/// The policy head starts at zero, which on a residual policy is exactly the
/// tuned controller's behaviour, so the baseline this reports is not a
/// different policy — it is this policy before the search moves it, and every
/// improvement is measured against the controller the flagship already ships.
pub struct G1TransformerTrainer {
    policy: TransformerG1Policy,
    evaluators: Vec<G1WalkingEvaluator>,
    optimizer: Option<CmaOptimizer>,
    pending: Option<CmaAsk>,
    cursor: usize,
    objectives: Vec<f64>,
    best: Vec<f64>,
    best_objective: f64,
    best_distance: f64,
    baseline_objective: f64,
    baseline_distance: f64,
    last: [f64; 3],
    evaluations: usize,
    restarts: usize,
    population: usize,
    sigma: f64,
    seed: u64,
    since_improvement: usize,
    /// Whether any candidate in the generation currently being evaluated beat
    /// the incumbent. Kept separate from `since_improvement`: deriving
    /// "improved" from that counter being zero is true on the first generation
    /// by construction, which silently pinned it at zero and meant the stall
    /// limit — and therefore every IPOP restart — could never be reached.
    improved_this_generation: bool,
    stall_limit: usize,
    generation: u64,
}

impl G1TransformerTrainer {
    /// `challenge`: 0 flat, 1 terrain-with-push, anything else both averaged.
    #[must_use]
    pub fn new(challenge: u32, duration_s: f64, sigma: f64, seed: u64) -> Self {
        let challenges: Vec<G1Challenge> = match challenge {
            0 => vec![G1Challenge::Flat],
            1 => vec![G1Challenge::TerrainAndPush],
            _ => vec![G1Challenge::Flat, G1Challenge::TerrainAndPush],
        };
        let evaluators = challenges
            .into_iter()
            .map(|challenge| {
                G1WalkingEvaluator::new(G1WalkingConfig {
                    task: G1Task::Walking,
                    challenge,
                    duration_s,
                    ..G1WalkingConfig::default()
                })
                .expect("walking evaluator builds")
            })
            .collect();

        let mut model_seed = MODEL_SEED;
        let model = GaitTransformer::new(
            browser_transformer_config(),
            1e-3,
            0.9,
            3e-4,
            &mut model_seed,
        );
        let log_std = PolicyLogStd::new(ACT_DIMS, 1e-3, -2.3);
        let mut policy = TransformerG1Policy::new(
            model,
            log_std,
            model_seed,
            g1_walking_curriculum_mean().to_vec(),
        );
        policy.set_exploration_enabled(false);

        let dim = policy.model.policy_head.params.len();
        let population = 4 + (3.0 * (dim as f64).ln()) as usize;
        let mut trainer = Self {
            policy,
            evaluators,
            optimizer: None,
            pending: None,
            cursor: 0,
            objectives: Vec::new(),
            best: vec![0.0; dim],
            best_objective: f64::INFINITY,
            best_distance: 0.0,
            baseline_objective: f64::INFINITY,
            baseline_distance: 0.0,
            last: [0.0; 3],
            evaluations: 0,
            restarts: 0,
            population,
            sigma,
            seed,
            since_improvement: 0,
            improved_this_generation: false,
            stall_limit: 25,
            generation: 0,
        };
        let zero = vec![0.0f64; dim];
        let (objective, distance, _) = trainer.score(&zero);
        trainer.baseline_objective = objective;
        trainer.baseline_distance = distance;
        trainer.best_objective = objective;
        trainer.best_distance = distance;
        trainer.start_run();
        trainer
    }

    fn start_run(&mut self) {
        self.optimizer = CmaOptimizer::new(CmaConfig {
            family: CmaFamily::LmCma,
            mean: self.best.clone(),
            sigma: self.sigma,
            // The browser run is open-ended — the user stops it, not a budget —
            // but the owner admits a finite evaluation budget and rejects a
            // saturating one, so each run gets a large finite allowance and is
            // simply restarted when it is spent.
            max_evaluations: self.population.saturating_mul(RUN_GENERATION_ALLOWANCE),
            seed: self
                .seed
                .wrapping_add((self.restarts as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15)),
            population_size: Some(self.population),
            memory: None,
        })
        .ok();
        self.pending = None;
        self.since_improvement = 0;
        self.improved_this_generation = false;
    }

    /// Generations since the incumbent last improved. Test-only: the stall
    /// counter drives every restart, and a counter pinned at zero looks
    /// identical from outside to a search that simply keeps improving.
    #[cfg(test)]
    #[must_use]
    pub const fn since_improvement(&self) -> usize {
        self.since_improvement
    }

    /// Mean objective, mean distance, and worst completed-step count over every
    /// configured challenge. Averaging rather than taking the best stops a
    /// candidate buying a flat-ground win with a terrain collapse.
    fn score(&mut self, head: &[f64]) -> (f64, f64, usize) {
        for (slot, value) in self.policy.model.policy_head.params.iter_mut().zip(head) {
            *slot = *value as f32;
        }
        let mut objective = 0.0;
        let mut distance = 0.0;
        let mut steps = usize::MAX;
        for evaluator in &self.evaluators {
            self.policy.begin_episode();
            let mut trace = EpisodeTrace::default();
            match evaluator.rollout_learned(&mut self.policy, &mut trace) {
                Ok(receipt) => {
                    objective += receipt.objective;
                    distance += receipt.distance_m;
                    steps = steps.min(receipt.completed_steps);
                }
                // A rollout the owner refuses is not a good candidate; rank it
                // last rather than letting it poison the distribution update.
                Err(_) => return (f64::MAX, 0.0, 0),
            }
            let _ = self.policy.take_collected();
        }
        let n = self.evaluators.len() as f64;
        (objective / n, distance / n, steps)
    }

    /// Advance the search by exactly one rollout.
    pub fn pump(&mut self) -> Vec<f64> {
        if self.optimizer.is_none() {
            return self.packet(TRAINER_STATUS_STOPPED);
        }
        if self.pending.is_none() {
            let Some(optimizer) = self.optimizer.as_mut() else {
                return self.packet(TRAINER_STATUS_STOPPED);
            };
            match optimizer.ask() {
                Ok(batch) => {
                    self.objectives = vec![f64::MAX; batch.len()];
                    self.cursor = 0;
                    self.pending = Some(batch);
                }
                Err(_) => {
                    // The run spent its allowance. Restart from the best point
                    // rather than ending the user's session.
                    self.restarts += 1;
                    self.start_run();
                    return self.packet(TRAINER_STATUS_GENERATION);
                }
            }
        }

        let candidate = {
            let Some(batch) = self.pending.as_ref() else {
                return self.packet(TRAINER_STATUS_STOPPED);
            };
            batch.candidates()[self.cursor].clone()
        };
        let (objective, distance, steps) = self.score(&candidate);
        let objective = if objective.is_finite() { objective } else { f64::MAX };
        self.objectives[self.cursor] = objective;
        self.last = [objective, distance, steps as f64];
        self.evaluations += 1;
        self.cursor += 1;
        if objective < self.best_objective {
            self.best_objective = objective;
            self.best_distance = distance;
            self.best.copy_from_slice(&candidate);
            self.improved_this_generation = true;
        }

        let complete = self
            .pending
            .as_ref()
            .is_some_and(|batch| self.cursor >= batch.len());
        if !complete {
            return self.packet(TRAINER_STATUS_RUNNING);
        }

        let batch = self.pending.take().expect("pending generation");
        let improved = std::mem::take(&mut self.improved_this_generation);
        let told = self
            .optimizer
            .as_mut()
            .map(|optimizer| optimizer.tell(&batch, &self.objectives));
        self.generation += 1;
        if !matches!(told, Some(Ok(_))) {
            self.optimizer = None;
            return self.packet(TRAINER_STATUS_STOPPED);
        }
        self.since_improvement = if improved { 0 } else { self.since_improvement + 1 };
        if self.since_improvement >= self.stall_limit {
            // IPOP: a converged run stands still, so restart from the best
            // point with a doubled population instead of burning the user's
            // time on a distribution that has stopped moving.
            self.restarts += 1;
            self.population = (self.population * 2).min(1024);
            self.start_run();
        }
        self.packet(TRAINER_STATUS_GENERATION)
    }

    fn packet(&self, status: f64) -> Vec<f64> {
        vec![
            status,
            self.generation as f64,
            self.evaluations as f64,
            self.best_objective,
            self.best_distance,
            self.baseline_objective,
            self.baseline_distance,
            self.restarts as f64,
            self.population as f64,
            self.last[0],
            self.last[1],
            self.last[2],
        ]
    }

    /// Current progress without advancing the search.
    #[must_use]
    pub fn progress(&self) -> Vec<f64> {
        self.packet(TRAINER_STATUS_RUNNING)
    }

    /// The best policy head found so far.
    #[must_use]
    pub fn best_head(&self) -> Vec<f64> {
        self.best.clone()
    }

    /// The best policy found so far, in the FSGT layout, ready to download.
    #[must_use]
    pub fn export_weights(&mut self) -> Vec<u8> {
        let restore: Vec<f32> = self.policy.model.policy_head.params.clone();
        for (slot, value) in self
            .policy
            .model
            .policy_head
            .params
            .iter_mut()
            .zip(&self.best)
        {
            *slot = *value as f32;
        }
        let bytes = dump_gait_transformer(&self.policy.model);
        self.policy.model.policy_head.params.copy_from_slice(&restore);
        bytes
    }
}

#[cfg(test)]
mod trainer_tests {
    use super::*;

    /// The trainer's baseline must BE the tuned controller, not merely resemble
    /// it: a zero policy head on a residual policy executes the composed
    /// controller exactly. If this drifts, every improvement the browser
    /// reports is measured against the wrong thing.
    #[test]
    fn baseline_is_the_tuned_controller() {
        let trainer = G1TransformerTrainer::new(0, 1.5, 0.002, 7);
        let packet = trainer.progress();
        assert_eq!(packet.len(), TRAINER_PACKET_WORDS);
        let baseline_objective = packet[5];
        let baseline_distance = packet[6];
        assert!(
            (baseline_objective - -59.413_444_990_206_53).abs() < 1e-6,
            "baseline objective drifted: {baseline_objective}"
        );
        assert!(
            (baseline_distance - 0.308_372_115_535_312_35).abs() < 1e-9,
            "baseline distance drifted: {baseline_distance}"
        );
    }

    /// The stall counter must actually count. It drives every IPOP restart,
    /// and the first version derived "did this generation improve?" from
    /// `since_improvement == 0`, which is true on the first generation by
    /// construction: the counter reset itself forever, no restart could ever
    /// fire, and from outside that is indistinguishable from a search that
    /// simply keeps improving.
    #[test]
    fn stall_counter_advances_when_a_generation_does_not_improve() {
        let mut trainer = G1TransformerTrainer::new(0, 1.5, 0.002, 7);
        let population = trainer.progress()[8] as usize;
        assert!(population > 1, "population must cover a real generation");

        let mut generations_checked = 0;
        for _ in 0..4 {
            let before = trainer.progress()[3];
            let mut closed = false;
            for _ in 0..population {
                if trainer.pump()[0] == TRAINER_STATUS_GENERATION {
                    closed = true;
                    break;
                }
            }
            assert!(closed, "a full generation should close within its population");
            let after = trainer.progress()[3];
            if after < before {
                assert_eq!(
                    trainer.since_improvement(),
                    0,
                    "an improving generation must reset the stall counter"
                );
            } else {
                assert!(
                    trainer.since_improvement() > 0,
                    "a generation that did not improve must advance the stall counter"
                );
                generations_checked += 1;
            }
        }
        assert!(
            generations_checked > 0,
            "no non-improving generation occurred, so the counter was never exercised"
        );
    }

    /// Pumping must actually search: the best objective has to fall below the
    /// baseline, and the exported artifact has to carry a head that moved.
    #[test]
    fn pumping_improves_on_the_tuned_controller() {
        let mut trainer = G1TransformerTrainer::new(0, 1.5, 0.002, 7);
        let baseline = trainer.progress()[5];
        for _ in 0..120 {
            let packet = trainer.pump();
            assert_eq!(packet.len(), TRAINER_PACKET_WORDS);
            assert!(packet[0] != TRAINER_STATUS_STOPPED, "search stopped early");
        }
        let progress = trainer.progress();
        assert!(
            progress[3] < baseline,
            "no improvement after 120 rollouts: {} vs {baseline}",
            progress[3]
        );
        assert!(progress[2] >= 120.0, "evaluations not counted");

        let head = trainer.best_head();
        let moved: f64 = head.iter().map(|w| w * w).sum();
        assert!(moved > 0.0, "best head is still all zeros");

        // The export must round-trip as a real FSGT artifact, not just bytes.
        let bytes = trainer.export_weights();
        assert_eq!(&bytes[0..4], WEIGHTS_MAGIC, "export lost its magic");
        assert_eq!(
            u32::from_le_bytes(bytes[4..8].try_into().unwrap()),
            LAYOUT_VERSION
        );
    }
}
