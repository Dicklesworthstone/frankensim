//! CMA-ES outer hyperparameter optimization: a (1+λ)-ES over the 8-D
//! hyperparameter vector, where each fitness evaluation runs one truncated
//! inner training. Reuses fs-dfo's natural-gradient IGO framework for the
//! sampling + recombination (the hyperparameter space is R^8, so the same
//! Euclidean CMA-ES machinery applies).
//!
//! This is the "demotion" from weight-space to HPO-space described in
//! RESEARCH_G1_LEARNING.md section 3.3: CMA-ES is genuinely strongest for
//! small-D expensive black-box problems, which is exactly what outer-loop
//! hyperparameter search is.
//!
//! # Review fixes (2026-08-29, RoseBasin fresh-eyes pass)
//!
//! The first version of this module had four compounding defects that let
//! its tests pass vacuously: (1) `incumbent_fitness` initialized to +∞
//! with higher-is-better fitness, so no evaluation could ever improve the
//! record; (2) `tell` updated the record but never anchored the incumbent
//! VECTOR on the best candidate; (3) the sigma adaptation's `improved`
//! flag was computed after the record update, so sigma only ever widened;
//! (4) log-scale sampling took `ln(linear_range)` (negative for ranges
//! < 1, so log-space candidates could only drift DOWN) and used an
//! uncentered [0,1) variate. All four are fixed here, and `sigma` now
//! actually scales the ask() step size.

pub const HPO_DIM: usize = 8;

/// Names of the 8 hyperparameters being optimized.
pub const HPO_NAMES: [&str; HPO_DIM] = [
    "inner_lr",
    "muon_momentum",
    "entropy_coef",
    "reward_weight_upright",
    "reward_weight_forward",
    "reward_weight_contact",
    "gae_lambda",
    "value_coef",
];

/// Default hyperparameter vector (the current production values).
pub const HPO_DEFAULT: [f32; HPO_DIM] = [
    3e-4, // inner_lr (Adam baseline for embed/head)
    0.9,  // muon_momentum
    0.01, // entropy_coef
    1.0,  // reward_weight_upright
    1.0,  // reward_weight_forward
    0.5,  // reward_weight_contact
    0.95, // gae_lambda
    0.5,  // value_coef
];

/// Per-hyperparameter search bounds (log-scale for LR-like params).
pub const HPO_BOUNDS: [(f32, f32); HPO_DIM] = [
    (1e-5, 3e-3), // inner_lr (log-uniform)
    (0.5, 0.999), // muon_momentum
    (0.001, 0.1), // entropy_coef (log-uniform)
    (0.1, 5.0),   // reward_weight_upright
    (0.1, 5.0),   // reward_weight_forward
    (0.05, 2.0),  // reward_weight_contact
    (0.8, 0.99),  // gae_lambda
    (0.1, 2.0),   // value_coef
];

/// Whether each hyperparameter is searched in log space.
pub const HPO_LOG_SCALE: [bool; HPO_DIM] = [true, false, true, false, false, false, false, false];

/// Clip a hyperparameter vector to the search bounds.
pub fn clip_to_bounds(hp: &[f32; HPO_DIM]) -> [f32; HPO_DIM] {
    let mut clipped = *hp;
    for i in 0..HPO_DIM {
        let (lo, hi) = HPO_BOUNDS[i];
        clipped[i] = clipped[i].max(lo).min(hi);
    }
    clipped
}

/// Simple (1+λ)-ES over the hyperparameter space. Each generation:
/// `ask` returns the incumbent plus λ−1 sigma-scaled perturbations of it;
/// `tell` anchors the incumbent on the best-evaluated candidate when it
/// beats the record and adapts sigma by the 1/5th rule (widen on a
/// generation with no record improvement, tighten on improvement).
pub struct HyperparameterSearch {
    pub incumbent: [f32; HPO_DIM],
    pub incumbent_fitness: f32,
    pub sigma: f32,
    pub generation: usize,
    /// (generation, incumbent fitness at that generation) — the running
    /// record, not the raw per-generation best.
    pub best_history: Vec<(usize, f32)>,
    state: u64,
}

impl HyperparameterSearch {
    pub fn new(sigma: f32, state: u64) -> Self {
        Self {
            incumbent: HPO_DEFAULT,
            // Higher-is-better fitness: the record starts at −∞ so the
            // first finite evaluation always anchors.
            incumbent_fitness: f32::NEG_INFINITY,
            sigma,
            generation: 0,
            best_history: Vec::new(),
            state,
        }
    }

    /// Ask for λ candidate hyperparameter vectors: the incumbent first,
    /// then λ−1 perturbations. Step size is `0.2 * sigma` of the
    /// parameter's range (linear) or log-range (log-scale params),
    /// symmetric around the incumbent.
    pub fn ask(&mut self, lambda: usize) -> Vec<[f32; HPO_DIM]> {
        let mut candidates = Vec::with_capacity(lambda);
        candidates.push(self.incumbent); // always include the incumbent
        for _ in 1..lambda {
            let mut candidate = self.incumbent;
            for i in 0..HPO_DIM {
                self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
                let z = ((self.state >> 11) as f64 / (1u64 << 53) as f64) as f32;
                let (lo, hi) = HPO_BOUNDS[i];
                let range = hi - lo;
                if HPO_LOG_SCALE[i] {
                    // Log-space width is ln(hi/lo) (always positive);
                    // symmetric centered step, sigma-scaled.
                    let log_width = (hi / lo).ln() as f32;
                    let log_cur = (candidate[i] / lo).ln();
                    let step = (z * 2.0 - 1.0) * log_width * 0.2 * self.sigma;
                    let new_log = (log_cur + step).max(0.0);
                    candidate[i] = lo * new_log.exp();
                } else {
                    candidate[i] += (z * 2.0 - 1.0) * range * 0.2 * self.sigma;
                }
            }
            candidates.push(clip_to_bounds(&candidate));
        }
        candidates
    }

    /// Tell: record the evaluated candidates, anchor the incumbent on
    /// the best one when it beats the record, and adapt sigma.
    /// `fitnesses[i]` = the inner training's final mean reward for
    /// `candidates[i]` (higher = better).
    pub fn tell(&mut self, candidates: &[[f32; HPO_DIM]], fitnesses: &[f32]) {
        assert_eq!(
            candidates.len(),
            fitnesses.len(),
            "hpo tell: candidates and fitnesses must pair up",
        );
        let mut best_idx = 0usize;
        let mut best = f32::NEG_INFINITY;
        for (idx, f) in fitnesses.iter().enumerate() {
            if *f > best {
                best = *f;
                best_idx = idx;
            }
        }
        // 1/5th rule, computed BEFORE the record update: improvement
        // means strictly beating the previous record.
        let improved = best > self.incumbent_fitness;
        if improved {
            self.incumbent_fitness = best;
            self.incumbent = candidates[best_idx];
            self.sigma *= 0.95;
        } else {
            self.sigma *= 1.15;
        }
        self.sigma = self.sigma.clamp(1e-6, 10.0);
        self.generation += 1;
        self.best_history
            .push((self.generation, self.incumbent_fitness));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hpo_candidates_are_bounded() {
        let mut hpo = HyperparameterSearch::new(0.5, 42);
        let candidates = hpo.ask(8);
        assert_eq!(candidates.len(), 8);
        for c in candidates.iter() {
            for i in 0..HPO_DIM {
                let (lo, hi) = HPO_BOUNDS[i];
                assert!(
                    c[i] >= lo && c[i] <= hi,
                    "param {i} = {} out of [{lo}, {hi}]",
                    c[i]
                );
            }
        }
    }

    #[test]
    fn hpo_first_tell_anchors_finite_record() {
        // The record starts at −∞ (higher-is-better); the first tell
        // must replace it with a finite value (regression for the +∞
        // initialization bug that made every evaluation vacuous).
        let mut hpo = HyperparameterSearch::new(0.5, 42);
        let candidates = hpo.ask(8);
        let fitnesses: Vec<f32> = candidates.iter().map(|_| 0.1).collect();
        hpo.tell(&candidates, &fitnesses);
        assert_eq!(hpo.incumbent_fitness, 0.1);
        assert!(hpo.incumbent_fitness.is_finite());
    }

    #[test]
    fn hpo_improves_mock_fitness() {
        let mut hpo = HyperparameterSearch::new(1.0, 42);
        // Mock fitness: higher inner_lr and higher entropy_coef = better.
        let fitness = |hp: &[f32; HPO_DIM]| -> f32 { hp[0] / 3e-3 + hp[2] / 0.1 };
        for _ in 0..20 {
            let candidates = hpo.ask(8);
            let fitnesses: Vec<f32> = candidates.iter().map(|c| fitness(c)).collect();
            hpo.tell(&candidates, &fitnesses);
        }
        // The default's fitness is 0.1 + 0.1 = 0.2. With anchoring, the
        // search must climb well past it (the max is 2.0 at the bounds).
        assert!(
            hpo.incumbent_fitness > 0.5,
            "HPO should improve on the default's 0.2, got {}",
            hpo.incumbent_fitness
        );
        // Anchoring must have MOVED the incumbent vector.
        assert!(
            hpo.incumbent != HPO_DEFAULT,
            "incumbent vector must move when improvements anchor"
        );
    }

    #[test]
    fn hpo_worse_batch_keeps_incumbent_vector() {
        let mut hpo = HyperparameterSearch::new(0.5, 42);
        // First tell anchors the −∞ record on the batch best.
        let candidates = hpo.ask(8);
        let fitnesses: Vec<f32> = candidates.iter().map(|_| 0.1).collect();
        hpo.tell(&candidates, &fitnesses);
        let anchored = hpo.incumbent;
        // A strictly worse batch must not move the anchored incumbent.
        let candidates2 = hpo.ask(8);
        let fitnesses2: Vec<f32> = candidates2.iter().map(|_| 0.05).collect();
        hpo.tell(&candidates2, &fitnesses2);
        assert_eq!(hpo.incumbent, anchored);
        assert_eq!(hpo.incumbent_fitness, 0.1);
    }

    #[test]
    fn hpo_sigma_adapts_both_ways() {
        // Regression for the dead sigma logic: widening on a no-record
        // improvement, tightening on an improvement. Checks are RELATIVE
        // (an absolute comparison across multiple adaptations is
        // order-dependent; widen 1.15 then tighten 0.95 nets above 1).
        let mut hpo = HyperparameterSearch::new(1.0, 7);
        // First batch anchors the −inf record -> tighten.
        let candidates = hpo.ask(8);
        let fitnesses: Vec<f32> = candidates.iter().map(|_| 0.1).collect();
        hpo.tell(&candidates, &fitnesses);
        let before_widen = hpo.sigma;
        // Batch that does not beat the record -> widen.
        let candidates = hpo.ask(8);
        let fitnesses: Vec<f32> = candidates.iter().map(|_| 0.05).collect();
        hpo.tell(&candidates, &fitnesses);
        let after_widen = hpo.sigma;
        assert!(
            after_widen > before_widen,
            "sigma must widen when the record does not improve"
        );
        // Batch that beats the record -> tighten relative to post-widen.
        let candidates = hpo.ask(8);
        let mut fitnesses: Vec<f32> = candidates.iter().map(|_| 0.05).collect();
        fitnesses[3] = 1.0; // a clear record-beater
        hpo.tell(&candidates, &fitnesses);
        assert!(
            hpo.sigma < after_widen,
            "a record improvement must tighten sigma"
        );
    }

    #[test]
    fn hpo_log_scale_samples_both_directions() {
        // Regression for the ln(linear_range) bias: log-scale candidates
        // must move BOTH up and down from the incumbent in log space.
        let mut hpo = HyperparameterSearch::new(1.0, 99);
        let lo = HPO_BOUNDS[0].0;
        let log_incumbent = (HPO_DEFAULT[0] / lo).ln();
        let mut ups = 0usize;
        let mut downs = 0usize;
        for _ in 0..20 {
            let candidates = hpo.ask(8);
            for c in candidates.iter().skip(1) {
                let log_c = (c[0] / lo).ln();
                if log_c > log_incumbent + 1e-6 {
                    ups += 1;
                } else if log_c < log_incumbent - 1e-6 {
                    downs += 1;
                }
            }
        }
        assert!(
            ups > 0 && downs > 0,
            "log-scale sampling must go both ways (ups={ups}, downs={downs})"
        );
    }
}
