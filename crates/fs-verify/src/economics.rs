//! ACCEPT/REJECT ECONOMICS + telemetry (bead lmp4.3): the layer that
//! turns "propose then verify" into a self-improving loop. Accept
//! OUTRIGHT when the certified bound meets tolerance; otherwise the
//! candidate becomes a WARM START whose iteration savings are
//! MEASURED — a rejected candidate is never wasted. Telemetry is
//! simultaneously the surrogate training signal, the planner's cost
//! model, and the drift detector: an accept-rate collapse in a regime
//! IS the distribution-shift alarm, and it demotes the offending
//! proposer there with hysteresis (no flapping).

use crate::fem1d::{
    Fem1dError, MAX_FEM1D_NEWTON_ITERATIONS, require_converged, solve_nonlinear, try_zeroed,
    validate_identity,
};
use crate::zoo::{Outcome, Registry, SpeculationQuery, speculate};
use std::collections::{BTreeMap, VecDeque};
use std::fmt::Write as _;

fn json_escape(value: &str) -> String {
    let mut escaped = String::with_capacity(value.len());
    for character in value.chars() {
        match character {
            '"' => escaped.push_str("\\\""),
            '\\' => escaped.push_str("\\\\"),
            '\u{08}' => escaped.push_str("\\b"),
            '\u{0c}' => escaped.push_str("\\f"),
            '\n' => escaped.push_str("\\n"),
            '\r' => escaped.push_str("\\r"),
            '\t' => escaped.push_str("\\t"),
            control if control <= '\u{1f}' => {
                let _ = write!(escaped, "\\u{:04x}", u32::from(control));
            }
            other => escaped.push(other),
        }
    }
    escaped
}

fn finite_json(value: f64) -> String {
    if value.is_finite() {
        format!("{value:.6e}")
    } else {
        "null".to_string()
    }
}

const MAX_DRIFT_KEYS: usize = 4_096;
const MAX_RETAINED_SAVINGS_PER_KEY: usize = 1_024;

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
    savings: BTreeMap<(String, String), VecDeque<u32>>,
}

impl DriftGuard {
    /// Record one try.
    pub fn record(
        &mut self,
        proposer: &str,
        regime: &str,
        accepted: bool,
    ) -> Result<(), Fem1dError> {
        validate_identity(proposer, "drift proposer")?;
        validate_identity(regime, "drift regime")?;
        if let Some((_, counts)) = self
            .counts
            .iter_mut()
            .find(|((stored_p, stored_r), _)| stored_p == proposer && stored_r == regime)
        {
            return increment_counts(counts, accepted);
        }
        if self.counts.len() >= MAX_DRIFT_KEYS {
            return Err(Fem1dError::ResourceLimit {
                resource: "drift telemetry keys",
                requested: self.counts.len().saturating_add(1),
                limit: MAX_DRIFT_KEYS,
            });
        }
        let mut counts = (0, 0);
        increment_counts(&mut counts, accepted)?;
        self.counts
            .insert((proposer.to_string(), regime.to_string()), counts);
        Ok(())
    }

    /// Record measured warm-start savings in a bounded most-recent window.
    pub fn record_savings(
        &mut self,
        proposer: &str,
        regime: &str,
        saved: u32,
    ) -> Result<(), Fem1dError> {
        validate_identity(proposer, "drift proposer")?;
        validate_identity(regime, "drift regime")?;
        if !self.counts.keys().any(|(stored_proposer, stored_regime)| {
            stored_proposer == proposer && stored_regime == regime
        }) {
            return Err(Fem1dError::InvalidScalar {
                field: "drift savings key",
                reason: "requires a previously recorded attempt",
            });
        }
        if let Some((_, samples)) = self
            .savings
            .iter_mut()
            .find(|((stored_p, stored_r), _)| stored_p == proposer && stored_r == regime)
        {
            push_savings_sample(samples, saved)?;
        } else {
            let mut samples = VecDeque::new();
            push_savings_sample(&mut samples, saved)?;
            self.savings
                .insert((proposer.to_string(), regime.to_string()), samples);
        }
        Ok(())
    }

    /// Accept rate, or a CONSERVATIVE prior (0.0) for a regime with
    /// zero telemetry — never a divide-by-zero, never optimism.
    #[must_use]
    pub fn accept_rate_or_prior(&self, proposer: &str, regime: &str) -> f64 {
        self.counts
            .iter()
            .find(|((stored_proposer, stored_regime), _)| {
                stored_proposer == proposer && stored_regime == regime
            })
            .map(|(_, counts)| counts)
            .map_or(
                0.0,
                |&(a, t)| if t == 0 { 0.0 } else { a as f64 / t as f64 },
            )
    }

    /// Drift check: demote proposers in regimes where the accept rate
    /// collapsed below `threshold` after ≥ `min_tries` (a single
    /// unlucky reject can never demote). Returns new demotions.
    pub fn check_drift(
        &mut self,
        threshold: f64,
        min_tries: u64,
    ) -> Result<Vec<(String, String)>, Fem1dError> {
        if !threshold.is_finite() || !(0.0..=1.0).contains(&threshold) {
            return Err(Fem1dError::InvalidScalar {
                field: "drift threshold",
                reason: "must be finite and within [0, 1]",
            });
        }
        if min_tries == 0 {
            return Err(Fem1dError::InvalidScalar {
                field: "drift min_tries",
                reason: "must be positive",
            });
        }
        let mut out = Vec::new();
        out.try_reserve_exact(self.counts.len())
            .map_err(|_| Fem1dError::AllocationFailed {
                stage: "drift demotions",
                requested: self.counts.len(),
            })?;
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
        Ok(out)
    }

    /// Is (proposer, regime) demoted?
    #[must_use]
    pub fn is_demoted(&self, proposer: &str, regime: &str) -> bool {
        self.demoted.keys().any(|(stored_proposer, stored_regime)| {
            stored_proposer == proposer && stored_regime == regime
        })
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
    ) -> Result<bool, Fem1dError> {
        validate_identity(proposer, "drift proposer")?;
        validate_identity(regime, "drift regime")?;
        if !promote_threshold.is_finite() || !(0.0..=1.0).contains(&promote_threshold) {
            return Err(Fem1dError::InvalidScalar {
                field: "probation promote threshold",
                reason: "must be finite and within [0, 1]",
            });
        }
        if probe_tries == 0 || probe_accepts > probe_tries {
            return Err(Fem1dError::InvalidScalar {
                field: "probation evidence",
                reason: "tries must be positive and accepts must not exceed tries",
            });
        }
        let key = (proposer.to_string(), regime.to_string());
        let Some(failures) = self.demoted.get(&key).copied() else {
            return Ok(true); // not demoted
        };
        let needed_tries = 5u64 << failures.min(10); // doubles per failure
        if probe_tries >= needed_tries
            && probe_tries > 0
            && (probe_accepts as f64 / probe_tries as f64) >= promote_threshold
        {
            self.demoted.remove(&key);
            // Reset the regime's window so old collapse data cannot
            // immediately re-demote or bias the recovered cost model.
            self.savings.remove(&key);
            self.counts.insert(key, (probe_accepts, probe_tries));
            Ok(true)
        } else {
            let failures = failures.checked_add(1).ok_or(Fem1dError::CounterOverflow {
                counter: "failed probations",
            })?;
            self.demoted.insert(key, failures);
            Ok(false)
        }
    }

    /// The kernel × regime × proposer dashboard rows.
    #[must_use]
    pub fn dashboard(&self, kernel: &str) -> Vec<String> {
        let kernel = json_escape(kernel);
        self.counts
            .iter()
            .map(|((p, r), &(a, t))| {
                let proposer = json_escape(p);
                let regime = json_escape(r);
                let med = self
                    .savings
                    .iter()
                    .find(|((stored_proposer, stored_regime), _)| {
                        stored_proposer == p && stored_regime == r
                    })
                    .map_or(0, |(_, values)| {
                        let mut s: Vec<u32> = values.iter().copied().collect();
                        s.sort_unstable();
                        s.get(s.len() / 2).copied().unwrap_or(0)
                    });
                let mut row = String::new();
                let _ = write!(
                    row,
                    "{{\"kernel\":\"{kernel}\",\"proposer\":\"{proposer}\",\"regime\":\"{regime}\",\
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

fn increment_counts(counts: &mut (u64, u64), accepted: bool) -> Result<(), Fem1dError> {
    let tries = counts.1.checked_add(1).ok_or(Fem1dError::CounterOverflow {
        counter: "drift tries",
    })?;
    let accepts = if accepted {
        counts.0.checked_add(1).ok_or(Fem1dError::CounterOverflow {
            counter: "drift accepts",
        })?
    } else {
        counts.0
    };
    *counts = (accepts, tries);
    Ok(())
}

fn push_savings_sample(samples: &mut VecDeque<u32>, saved: u32) -> Result<(), Fem1dError> {
    if samples.len() == MAX_RETAINED_SAVINGS_PER_KEY {
        samples.pop_front();
    } else {
        samples
            .try_reserve(1)
            .map_err(|_| Fem1dError::AllocationFailed {
                stage: "drift savings window",
                requested: samples.len().saturating_add(1),
            })?;
    }
    samples.push_back(saved);
    Ok(())
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
    let proposer_id = json_escape(proposer_id);
    let bound = finite_json(bound);
    format!(
        "{{\"proposer_id\":\"{proposer_id}\",\"accepted\":{accepted},\
         \"bound\":{bound},\"iterations_saved\":{iterations_saved}}}"
    )
}

/// The control loop: speculate; accept outright on a certified pass;
/// otherwise warm-start the true (nonlinear-class) solve from the best
/// rejected candidate and MEASURE the savings. Deterministic given the
/// query and registry state.
///
/// # Errors
/// Returns [`Fem1dError`] for invalid/bounded-work inputs, proposal failures,
/// numerical failures, or finite nonconvergence. A failed solve never becomes
/// an iteration-savings observation.
pub fn run_speculative(
    query: &SpeculationQuery,
    registry: &Registry,
    zoo_telemetry: &mut crate::zoo::ZooTelemetry,
    guard: &mut DriftGuard,
    max_iter: u32,
) -> Result<EconDecision, Fem1dError> {
    if max_iter > MAX_FEM1D_NEWTON_ITERATIONS {
        return Err(Fem1dError::ResourceLimit {
            resource: "Newton iterations",
            requested: max_iter as usize,
            limit: MAX_FEM1D_NEWTON_ITERATIONS as usize,
        });
    }
    match speculate(query, registry, zoo_telemetry)? {
        Outcome::Accepted(ans) => {
            for proposer in ans.rejected_before() {
                guard.record(proposer, &query.regime, false)?;
            }
            guard.record(ans.proposer(), &query.regime, true)?;
            Ok(EconDecision::AcceptedOutright {
                proposer: ans.proposer(),
                bound: ans.report().bound.hi,
            })
        }
        Outcome::AllRejected(rejected) => {
            // Consume first-pass verified rejects. Proposers are never rerun:
            // stateful proposal side effects cannot change the warm-start choice.
            let (attempted, best) = rejected.into_parts();
            let zero = try_zeroed("economic cold start", query.problem.mesh.len())?;
            let cold = solve_nonlinear(&query.problem, &zero, max_iter)?;
            require_converged(&cold, "economic cold solve")?;
            if let Some(rejected) = best {
                let warm = solve_nonlinear(&query.problem, &rejected.candidate, max_iter)?;
                require_converged(&warm, "economic warm solve")?;
                let raw = i64::from(cold.iterations) - i64::from(warm.iterations);
                let saved = u32::try_from(raw.max(0)).unwrap_or(0);
                for proposer in attempted {
                    guard.record(proposer, &query.regime, false)?;
                }
                guard.record_savings(rejected.proposer, &query.regime, saved)?;
                Ok(EconDecision::WarmStarted {
                    cold: cold.iterations,
                    warm: warm.iterations,
                    saved,
                    raw_delta: raw,
                })
            } else {
                for proposer in attempted {
                    guard.record(proposer, &query.regime, false)?;
                }
                Ok(EconDecision::ColdSolve {
                    iterations: cold.iterations,
                })
            }
        }
        Outcome::NoCandidates => {
            let zero = try_zeroed("economic cold start", query.problem.mesh.len())?;
            let cold = solve_nonlinear(&query.problem, &zero, max_iter)?;
            require_converged(&cold, "economic cold solve")?;
            Ok(EconDecision::ColdSolve {
                iterations: cold.iterations,
            })
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn drift_state_is_bounded_and_counter_updates_are_atomic() {
        let mut guard = DriftGuard::default();
        for index in 0..MAX_DRIFT_KEYS {
            guard
                .record("p", &format!("r{index}"), false)
                .expect("key within drift cap");
        }
        assert!(matches!(
            guard.record("p", "overflow", false),
            Err(Fem1dError::ResourceLimit {
                resource: "drift telemetry keys",
                ..
            })
        ));

        let mut overflow = DriftGuard::default();
        overflow
            .counts
            .insert(("p".to_string(), "r".to_string()), (u64::MAX, 7));
        assert!(matches!(
            overflow.record("p", "r", true),
            Err(Fem1dError::CounterOverflow {
                counter: "drift accepts"
            })
        ));
        assert_eq!(
            overflow.counts.get(&("p".to_string(), "r".to_string())),
            Some(&(u64::MAX, 7)),
            "a failed counter update must not partially increment tries"
        );
    }

    #[test]
    fn savings_are_a_bounded_recent_window_and_policy_inputs_refuse() {
        let mut guard = DriftGuard::default();
        guard.record("p", "r", false).expect("seed attempt");
        for saved in 0..1_100u32 {
            guard
                .record_savings("p", "r", saved)
                .expect("bounded savings sample");
        }
        let samples = guard
            .savings
            .get(&("p".to_string(), "r".to_string()))
            .expect("savings window");
        assert_eq!(samples.len(), MAX_RETAINED_SAVINGS_PER_KEY);
        assert_eq!(samples.front(), Some(&76));
        assert_eq!(samples.back(), Some(&1_099));
        assert!(guard.dashboard("kernel")[0].contains("\"median_savings\":588"));

        assert!(matches!(
            guard.record_savings("missing", "r", 1),
            Err(Fem1dError::InvalidScalar {
                field: "drift savings key",
                ..
            })
        ));
        assert!(matches!(
            guard.check_drift(f64::NAN, 1),
            Err(Fem1dError::InvalidScalar {
                field: "drift threshold",
                ..
            })
        ));
        assert!(matches!(
            guard.check_drift(0.5, 0),
            Err(Fem1dError::InvalidScalar {
                field: "drift min_tries",
                ..
            })
        ));
        assert!(matches!(
            guard.probation("p", "r", 2, 1, 0.5),
            Err(Fem1dError::InvalidScalar {
                field: "probation evidence",
                ..
            })
        ));
    }
}
