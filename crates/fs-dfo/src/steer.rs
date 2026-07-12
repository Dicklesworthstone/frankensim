//! MOO INTEGRATION LANES (bead qlvf): the two lanes gated on other
//! crates, now landed.
//!
//! (a) INTERACTIVE STEERING via WORLD-FORKING (P9): an agent
//! re-weights objectives mid-campaign and the study FORKS rather than
//! mutates — the parent branch is untouched, the steering event is a
//! ledger-ready record, and BOTH branches replay bitwise from their
//! lineage (state = the deterministic (population, stream-index,
//! weights) triple).
//!
//! (b) CHANCE CONSTRAINTS live in `fs_uq::chance` (fs-uq sits above
//! fs-dfo through fs-bo, so the integration points that way — the
//! dependency cycle the first draft hit is the layer diagram talking).

use crate::moo::Individual;

/// The forkable study state — plain data, resumable by construction.
#[derive(Debug, Clone)]
pub struct StudyState {
    /// Current population.
    pub population: Vec<Individual>,
    /// The deterministic stream cursor (advances per generation).
    pub stream_index: u64,
    /// Objective weights (the steering surface).
    pub weights: Vec<f64>,
}

/// One ledger-ready steering record.
#[derive(Debug, Clone, PartialEq)]
pub struct SteerEvent {
    /// Logical time (generation count at the fork).
    pub at_generation: u64,
    /// Weights before.
    pub from: Vec<f64>,
    /// Weights after.
    pub to: Vec<f64>,
}

impl SteerEvent {
    /// Ledger payload.
    #[must_use]
    pub fn to_json(&self) -> String {
        format!(
            "{{\"at\":{},\"from\":{:?},\"to\":{:?}}}",
            self.at_generation, self.from, self.to
        )
    }
}

/// A steered study branch: state + the lineage that reproduces it.
#[derive(Debug, Clone)]
pub struct SteeredStudy {
    /// The branch's current state.
    pub state: StudyState,
    /// Base seed (shared by all branches of the world tree).
    pub seed: u64,
    /// The steering lineage (ledgered ops, replayable).
    pub lineage: Vec<SteerEvent>,
}

fn unit(seed: u64, k: u64) -> f64 {
    let mut z = seed ^ 0x9e37_79b9_7f4a_7c15u64.wrapping_mul(k.wrapping_add(1));
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^= z >> 31;
    (z >> 11) as f64 / (1u64 << 53) as f64
}

impl SteeredStudy {
    /// Start a study: uniform-random population in the box, uniform
    /// weights, empty lineage.
    #[must_use]
    pub fn start(
        objectives: &mut dyn FnMut(&[f64]) -> Vec<f64>,
        dim: usize,
        bounds: (f64, f64),
        pop: usize,
        n_obj: usize,
        seed: u64,
    ) -> SteeredStudy {
        // `advance` selects parents via `... % population.len()`; an empty
        // population would divide by zero. Fail closed at construction.
        assert!(pop > 0, "steered study needs a positive population");
        let (lo, hi) = bounds;
        let population: Vec<Individual> = (0..pop)
            .map(|i| {
                let x: Vec<f64> = (0..dim)
                    .map(|j| lo + (hi - lo) * unit(seed, (i * dim + j) as u64))
                    .collect();
                let f = objectives(&x);
                Individual { x, f }
            })
            .collect();
        SteeredStudy {
            state: StudyState {
                population,
                stream_index: 0,
                weights: vec![1.0 / n_obj as f64; n_obj],
            },
            seed,
            lineage: Vec::new(),
        }
    }

    /// FORK: a new branch with re-weighted objectives. The PARENT is
    /// untouched (world-forking, not mutation); the steering event is
    /// recorded in the child's lineage.
    #[must_use]
    pub fn fork(&self, new_weights: Vec<f64>) -> SteeredStudy {
        let mut child = self.clone();
        child.lineage.push(SteerEvent {
            at_generation: self.state.stream_index,
            from: self.state.weights.clone(),
            to: new_weights.clone(),
        });
        child.state.weights = new_weights;
        child
    }

    /// Advance `generations` of deterministic weighted evolution
    /// (tournament on the weighted sum + box-clamped Gaussian-ish
    /// mutation, all counter-based from (seed, stream_index)).
    pub fn advance(
        &mut self,
        objectives: &mut dyn FnMut(&[f64]) -> Vec<f64>,
        bounds: (f64, f64),
        generations: u64,
    ) {
        let (lo, hi) = bounds;
        let sigma = 0.1 * (hi - lo);
        for _ in 0..generations {
            let g = self.state.stream_index;
            let w = self.state.weights.clone();
            let score =
                |ind: &Individual| -> f64 { ind.f.iter().zip(&w).map(|(f, wi)| f * wi).sum() };
            let n = self.state.population.len();
            let mut next = Vec::with_capacity(n);
            for i in 0..n {
                // Deterministic binary tournament.
                let k = self.seed ^ (g.wrapping_mul(0x9e37) ^ i as u64);
                #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
                let a = (unit(k, 0) * n as f64) as usize % n;
                #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
                let b = (unit(k, 1) * n as f64) as usize % n;
                let parent = if score(&self.state.population[a]) <= score(&self.state.population[b])
                {
                    &self.state.population[a]
                } else {
                    &self.state.population[b]
                };
                // Counter-based mutation.
                let x: Vec<f64> = parent
                    .x
                    .iter()
                    .enumerate()
                    .map(|(j, v)| {
                        let z: f64 =
                            (0..6).map(|q| unit(k, 2 + j as u64 * 6 + q)).sum::<f64>() - 3.0;
                        (v + sigma * z / 3.0).clamp(lo, hi)
                    })
                    .collect();
                let f = objectives(&x);
                next.push(Individual { x, f });
            }
            // Elitist merge on the weighted score.
            let mut all = std::mem::take(&mut self.state.population);
            all.extend(next);
            all.sort_by(|p, q| {
                score(p)
                    .total_cmp(&score(q))
                    .then_with(|| p.x[0].total_cmp(&q.x[0]))
            });
            all.truncate(n);
            self.state.population = all;
            self.state.stream_index += 1;
        }
    }

    /// The branch's best weighted score (diagnostics).
    #[must_use]
    pub fn best_score(&self) -> f64 {
        self.state
            .population
            .iter()
            .map(|ind| {
                ind.f
                    .iter()
                    .zip(&self.state.weights)
                    .map(|(f, w)| f * w)
                    .sum::<f64>()
            })
            .fold(f64::INFINITY, f64::min)
    }

    /// A deterministic fingerprint of the state (replay witness).
    /// Canonical replay identity encoding (gp3.14): the former bare
    /// concatenation of variable-length x/f/weight streams was
    /// non-injective — a value could migrate across a section boundary
    /// (ind.x tail vs ind.f head) without moving the hash. Sections
    /// now carry typed length prefixes.
    #[must_use]
    pub fn fingerprint(&self) -> u64 {
        let mut b = fs_obs::ident::IdentityBuilder::new("dfo-steered-study")
            .u64("stream_index", self.state.stream_index)
            .u64("population", self.state.population.len() as u64);
        for ind in &self.state.population {
            b = b.u64("x_len", ind.x.len() as u64);
            for &v in &ind.x {
                b = b.f64_bits("x", v);
            }
            b = b.u64("f_len", ind.f.len() as u64);
            for &v in &ind.f {
                b = b.f64_bits("f", v);
            }
        }
        b = b.u64("weights_len", self.state.weights.len() as u64);
        for &v in &self.state.weights {
            b = b.f64_bits("w", v);
        }
        b.finish().root()
    }
}
