//! ACCEPT/REJECT ECONOMICS + telemetry (bead lmp4.3): the layer that
//! turns "propose then verify" into a self-improving loop. Accept
//! OUTRIGHT when the certified bound meets tolerance; otherwise the
//! candidate becomes a WARM START whose iteration savings are
//! MEASURED — a rejected candidate is never wasted. Telemetry is
//! simultaneously the surrogate training signal, the planner's cost
//! model, and the drift detector: an accept-rate collapse in a regime
//! IS the distribution-shift alarm, and it demotes the offending
//! proposer there with hysteresis (no flapping).

use crate::fem1d::solve_nonlinear;
use crate::zoo::{Outcome, Registry, SpeculationQuery, speculate};
use std::collections::BTreeMap;
use std::fmt::Write as _;

/// One speculation's economic outcome.
#[derive(Debug, Clone)]
pub enum EconDecision {
    /// The certified bound met tolerance: the answer ships with its
    /// verified color, no solve at all.
    AcceptedOutright {
        /// Which proposer won.
        proposer: &'static str,
        /// The certified bound.
        bound: f64,
    },
    /// Rejected candidates warm-start the true solve; savings are
    /// measured and RECORDED CLAMPED at ≥ 0 (a worse-than-cold start
    /// is never a win) with the raw delta logged.
    WarmStarted {
        /// Iterations from cold.
        cold: u32,
        /// Iterations from the candidate.
        warm: u32,
        /// Recorded savings (`max(cold − warm, 0)`).
        saved: u32,
        /// Raw delta (negative when the warm start was WORSE).
        raw_delta: i64,
    },
    /// Nothing to try: the full cold solve.
    ColdSolve {
        /// Iterations spent.
        iterations: u32,
    },
}

/// Hysteresis-guarded drift detection state per (proposer, regime).
#[derive(Debug, Default)]
pub struct DriftGuard {
    counts: BTreeMap<(String, String), (u64, u64)>, // (accepts, tries)
    demoted: BTreeMap<(String, String), u32>,       // failed probations
    savings: BTreeMap<(String, String), Vec<u32>>,
}

impl DriftGuard {
    /// Record one try.
    pub fn record(&mut self, proposer: &str, regime: &str, accepted: bool) {
        let e = self
            .counts
            .entry((proposer.to_string(), regime.to_string()))
            .or_insert((0, 0));
        e.1 += 1;
        if accepted {
            e.0 += 1;
        }
    }

    /// Record measured warm-start savings for the regime's ledger.
    pub fn record_savings(&mut self, proposer: &str, regime: &str, saved: u32) {
        self.savings
            .entry((proposer.to_string(), regime.to_string()))
            .or_default()
            .push(saved);
    }

    /// Accept rate, or a CONSERVATIVE prior (0.0) for a regime with
    /// zero telemetry — never a divide-by-zero, never optimism.
    #[must_use]
    pub fn accept_rate_or_prior(&self, proposer: &str, regime: &str) -> f64 {
        self.counts
            .get(&(proposer.to_string(), regime.to_string()))
            .map_or(0.0, |&(a, t)| if t == 0 { 0.0 } else { a as f64 / t as f64 })
    }

    /// Drift check: demote proposers in regimes where the accept rate
    /// collapsed below `threshold` after ≥ `min_tries` (a single
    /// unlucky reject can never demote). Returns new demotions.
    pub fn check_drift(&mut self, threshold: f64, min_tries: u64) -> Vec<(String, String)> {
        let mut out = Vec::new();
        for ((p, r), &(a, t)) in &self.counts {
            let key = (p.clone(), r.clone());
            if t >= min_tries
                && (a as f64 / t as f64) < threshold
                && !self.demoted.contains_key(&key)
            {
                self.demoted.insert(key.clone(), 0);
                out.push(key);
            }
        }
        out
    }

    /// Is (proposer, regime) demoted?
    #[must_use]
    pub fn is_demoted(&self, proposer: &str, regime: &str) -> bool {
        self.demoted
            .contains_key(&(proposer.to_string(), regime.to_string()))
    }

    /// Probation: re-admit a demoted proposer for a probe window; it
    /// re-promotes ONLY if the probe rate clears `promote_threshold`
    /// (strictly above the demotion threshold — hysteresis), else it
    /// stays demoted with the failure counted. No flapping: each
    /// failed probation doubles the evidence needed for the next.
    pub fn probation(
        &mut self,
        proposer: &str,
        regime: &str,
        probe_accepts: u64,
        probe_tries: u64,
        promote_threshold: f64,
    ) -> bool {
        let key = (proposer.to_string(), regime.to_string());
        let Some(failures) = self.demoted.get(&key).copied() else {
            return true; // not demoted
        };
        let needed_tries = 5u64 << failures.min(10); // doubles per failure
        if probe_tries >= needed_tries
            && probe_tries > 0
            && (probe_accepts as f64 / probe_tries as f64) >= promote_threshold
        {
            self.demoted.remove(&key);
            // Reset the regime's window so old collapse data cannot
            // immediately re-demote.
            self.counts.insert(key, (probe_accepts, probe_tries));
            true
        } else {
            self.demoted.insert(key, failures + 1);
            false
        }
    }

    /// The kernel × regime × proposer dashboard rows.
    #[must_use]
    pub fn dashboard(&self, kernel: &str) -> Vec<String> {
        self.counts
            .iter()
            .map(|((p, r), &(a, t))| {
                let med = self
                    .savings
                    .get(&(p.clone(), r.clone()))
                    .map_or(0, |v| {
                        let mut s = v.clone();
                        s.sort_unstable();
                        s.get(s.len() / 2).copied().unwrap_or(0)
                    });
                let mut row = String::new();
                let _ = write!(
                    row,
                    "{{\"kernel\":\"{kernel}\",\"proposer\":\"{p}\",\"regime\":\"{r}\",\
                     \"accepts\":{a},\"tries\":{t},\"rate\":{:.4},\
                     \"median_savings\":{med},\"demoted\":{}}}",
                    a as f64 / t.max(1) as f64,
                    self.is_demoted(p, r)
                );
                row
            })
            .collect()
    }
}

/// The four solve-node ledger fields (the schema amendment, stored as
/// a `speculation` extension record in fs-ledger).
#[must_use]
pub fn solve_node_record(
    proposer_id: &str,
    accepted: bool,
    bound: f64,
    iterations_saved: u32,
) -> String {
    format!(
        "{{\"proposer_id\":\"{proposer_id}\",\"accepted\":{accepted},\
         \"bound\":{bound:.6e},\"iterations_saved\":{iterations_saved}}}"
    )
}

/// The control loop: speculate; accept outright on a certified pass;
/// otherwise warm-start the true (nonlinear-class) solve from the best
/// rejected candidate and MEASURE the savings. Deterministic given the
/// query and registry state.
pub fn run_speculative(
    query: &SpeculationQuery,
    registry: &Registry,
    zoo_telemetry: &mut crate::zoo::ZooTelemetry,
    guard: &mut DriftGuard,
    max_iter: u32,
) -> EconDecision {
    match speculate(query, registry, zoo_telemetry) {
        Outcome::Accepted(ans) => {
            guard.record(ans.proposer(), &query.regime, true);
            EconDecision::AcceptedOutright {
                proposer: ans.proposer(),
                bound: ans.report().bound.hi,
            }
        }
        Outcome::AllRejected { .. } => {
            // A rejected candidate is never wasted: gather the best
            // (smallest certified bound) as the warm start.
            let mut best: Option<(&'static str, Vec<f64>, f64)> = None;
            for p in registry_iter(registry) {
                if zoo_telemetry.is_demoted(p.name(), &query.regime) {
                    continue;
                }
                if let Some(prop) = p.propose(query) {
                    let rep = crate::estimator::verify(&query.problem, &prop.candidate, query.tolerance);
                    guard.record(p.name(), &query.regime, false);
                    let better = best
                        .as_ref()
                        .is_none_or(|(_, _, b)| rep.bound.hi < *b);
                    if rep.bound.hi.is_finite() && better {
                        best = Some((p.name(), prop.candidate, rep.bound.hi));
                    }
                }
            }
            let zero = vec![0.0; query.problem.mesh.len()];
            let (_, cold) = solve_nonlinear(&query.problem, &zero, max_iter);
            match best {
                Some((name, candidate, _)) => {
                    let (_, warm) = solve_nonlinear(&query.problem, &candidate, max_iter);
                    let raw = i64::from(cold) - i64::from(warm);
                    let saved = u32::try_from(raw.max(0)).unwrap_or(0);
                    guard.record_savings(name, &query.regime, saved);
                    EconDecision::WarmStarted {
                        cold,
                        warm,
                        saved,
                        raw_delta: raw,
                    }
                }
                None => EconDecision::ColdSolve { iterations: cold },
            }
        }
        Outcome::NoCandidates => {
            let zero = vec![0.0; query.problem.mesh.len()];
            let (_, cold) = solve_nonlinear(&query.problem, &zero, max_iter);
            EconDecision::ColdSolve { iterations: cold }
        }
    }
}

fn registry_iter(registry: &Registry) -> impl Iterator<Item = &(dyn crate::zoo::Proposer)> {
    registry.proposers_dyn()
}
