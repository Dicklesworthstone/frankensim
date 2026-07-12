//! The PROPOSER ZOO (bead lmp4.2): untrusted, hot-swappable fast
//! models behind ONE `propose()` interface, feeding the certified
//! verifier. The speculative-decoding pattern transplanted to
//! numerics, licensed by the check/produce asymmetry: checking is
//! cheap, so proposers may be as reckless as they like.
//!
//! THE SAFETY INVARIANT LIVES IN THE TYPES: a [`CertifiedAnswer`] has
//! no public constructor — the only way one comes into existence is
//! [`speculate`] passing a candidate through [`crate::estimator::verify`]
//! and receiving an accept. A bad proposer can waste a check; it can
//! never corrupt a result.
//!
//! Self-reported confidence is ADVISORY ONLY: it orders which proposer
//! gets tried first (the economics), and it NEVER enters the accept
//! decision — a NaN-confidence candidate that verifies is accepted; a
//! confidence-1.0 garbage candidate is rejected.

use crate::estimator::{VerifierReport, verify};
use crate::fem1d::{
    Fem1dError, MAX_FEM1D_MESH_NODES, MmsProblem, solve_p1, validate_candidate,
    validate_finite_scalar, validate_identity, validate_problem, validate_tolerance,
};
use std::collections::BTreeMap;
use std::fmt::Write as _;

const MAX_ZOO_PROPOSERS: usize = 4_096;
const MAX_NEIGHBOR_CACHE_ENTRIES: usize = 4_096;
const MAX_ZOO_TELEMETRY_KEYS: usize = 4_096;

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

/// A speculation query: the problem, where in design space it sits,
/// and the tolerance the answer must certify against.
#[derive(Debug, Clone)]
pub struct SpeculationQuery {
    /// The target problem.
    pub problem: MmsProblem,
    /// The design-space coordinate (v0: one parameter).
    pub theta: f64,
    /// The certification tolerance.
    pub tolerance: f64,
    /// The regime key (telemetry / demotion granularity).
    pub regime: String,
}

/// One proposer's candidate.
#[derive(Debug, Clone)]
pub struct Proposal {
    /// Nodal candidate values.
    pub candidate: Vec<f64>,
    /// Self-reported confidence — ADVISORY ONLY (ordering hint;
    /// never enters any certificate or accept decision).
    pub confidence: f64,
}

/// The uniform proposer interface. None is trusted; all are useful;
/// each is independently retirable.
pub trait Proposer: Send + Sync {
    /// Stable name (registry, telemetry, ledger rows).
    fn name(&self) -> &'static str;
    /// Produce a candidate, or decline (`None` = nothing to offer).
    ///
    /// # Errors
    /// Returns [`Fem1dError`] when bounded proposal construction cannot proceed.
    fn propose(&self, query: &SpeculationQuery) -> Result<Option<Proposal>, Fem1dError>;
}

/// A CERTIFIED answer: candidate + the verifier's report (verified
/// color included). NO public constructor — the type is the proof
/// that the verifier said yes.
#[derive(Debug, Clone)]
pub struct CertifiedAnswer {
    candidate: Vec<f64>,
    report: VerifierReport,
    proposer: &'static str,
    rejected_before: Vec<&'static str>,
}

impl CertifiedAnswer {
    /// The accepted candidate.
    #[must_use]
    pub fn candidate(&self) -> &[f64] {
        &self.candidate
    }

    /// The verifier's report (bound, color, tolerance).
    #[must_use]
    pub fn report(&self) -> &VerifierReport {
        &self.report
    }

    /// Which proposer won.
    #[must_use]
    pub fn proposer(&self) -> &'static str {
        self.proposer
    }

    pub(crate) fn rejected_before(&self) -> &[&'static str] {
        &self.rejected_before
    }
}

#[derive(Debug, Clone)]
pub(crate) struct RejectedCandidate {
    pub(crate) proposer: &'static str,
    pub(crate) candidate: Vec<f64>,
    pub(crate) report: VerifierReport,
}

/// Sealed details from one first-pass all-rejected speculation.
///
/// Callers can inspect the attempt count, while the integrated economics loop
/// consumes the retained verified rejects without invoking proposers again.
#[derive(Debug, Clone)]
pub struct AllRejected {
    tried: u32,
    attempted: Vec<&'static str>,
    best: Option<RejectedCandidate>,
}

impl AllRejected {
    /// Number of first-pass proposals checked.
    #[must_use]
    pub fn tried(&self) -> u32 {
        self.tried
    }

    pub(crate) fn into_parts(self) -> (Vec<&'static str>, Option<RejectedCandidate>) {
        (self.attempted, self.best)
    }
}

/// The speculation outcome.
#[derive(Debug, Clone)]
pub enum Outcome {
    /// A proposal passed the verifier.
    Accepted(Box<CertifiedAnswer>),
    /// Proposals were tried; every one was rejected (fall back to the
    /// full solve — nothing was corrupted, only checks were spent).
    AllRejected(AllRejected),
    /// No enabled proposer had a candidate to offer.
    NoCandidates,
}

/// Per-proposer, per-regime accept telemetry with the auto-demotion
/// hook (an accept-rate collapse disables a proposer in that regime).
#[derive(Debug, Default)]
pub struct ZooTelemetry {
    counts: BTreeMap<(String, String), (u64, u64)>, // (accepts, tries)
    demoted: BTreeMap<(String, String), bool>,
}

impl ZooTelemetry {
    fn record(&mut self, proposer: &str, regime: &str, accepted: bool) -> Result<(), Fem1dError> {
        validate_identity(proposer, "proposer name")?;
        validate_identity(regime, "regime")?;
        if let Some((_, counts)) = self
            .counts
            .iter_mut()
            .find(|((stored_p, stored_r), _)| stored_p == proposer && stored_r == regime)
        {
            return increment_counts(counts, accepted, "zoo tries", "zoo accepts");
        }
        if self.counts.len() >= MAX_ZOO_TELEMETRY_KEYS {
            return Err(Fem1dError::ResourceLimit {
                resource: "zoo telemetry keys",
                requested: self.counts.len().saturating_add(1),
                limit: MAX_ZOO_TELEMETRY_KEYS,
            });
        }
        let mut counts = (0, 0);
        increment_counts(&mut counts, accepted, "zoo tries", "zoo accepts")?;
        self.counts
            .insert((proposer.to_string(), regime.to_string()), counts);
        Ok(())
    }

    /// Accept rate for (proposer, regime).
    #[must_use]
    pub fn accept_rate(&self, proposer: &str, regime: &str) -> Option<f64> {
        self.counts
            .iter()
            .find(|((stored_proposer, stored_regime), _)| {
                stored_proposer == proposer && stored_regime == regime
            })
            .map(|(_, counts)| counts)
            .map(|&(a, t)| a as f64 / t.max(1) as f64)
    }

    /// Demote proposers whose accept rate in a regime collapsed below
    /// `threshold` after at least `min_tries` attempts. Returns the
    /// demotions performed.
    pub fn demote_collapsed(
        &mut self,
        threshold: f64,
        min_tries: u64,
    ) -> Result<Vec<(String, String)>, Fem1dError> {
        if !threshold.is_finite() || !(0.0..=1.0).contains(&threshold) {
            return Err(Fem1dError::InvalidScalar {
                field: "zoo demotion threshold",
                reason: "must be finite and within [0, 1]",
            });
        }
        if min_tries == 0 {
            return Err(Fem1dError::InvalidScalar {
                field: "zoo demotion min_tries",
                reason: "must be positive",
            });
        }
        let mut out = Vec::new();
        out.try_reserve_exact(self.counts.len())
            .map_err(|_| Fem1dError::AllocationFailed {
                stage: "zoo demotions",
                requested: self.counts.len(),
            })?;
        for ((p, r), &(a, t)) in &self.counts {
            if t >= min_tries && (a as f64 / t as f64) < threshold {
                let key = (p.clone(), r.clone());
                if !self.demoted.get(&key).copied().unwrap_or(false) {
                    self.demoted.insert(key.clone(), true);
                    out.push(key);
                }
            }
        }
        Ok(out)
    }

    /// Is a proposer demoted in a regime?
    #[must_use]
    pub fn is_demoted(&self, proposer: &str, regime: &str) -> bool {
        self.demoted
            .iter()
            .any(|((stored_proposer, stored_regime), demoted)| {
                stored_proposer == proposer && stored_regime == regime && *demoted
            })
    }

    /// Ledger rows (per proposer × regime).
    #[must_use]
    pub fn rows(&self) -> Vec<String> {
        self.counts
            .iter()
            .map(|((p, r), &(a, t))| {
                let proposer = json_escape(p);
                let regime = json_escape(r);
                let mut s = String::new();
                let _ = write!(
                    s,
                    "{{\"proposer\":\"{proposer}\",\"regime\":\"{regime}\",\"accepts\":{a},\
                     \"tries\":{t},\"rate\":{:.4},\"demoted\":{}}}",
                    a as f64 / t.max(1) as f64,
                    self.is_demoted(p, r)
                );
                s
            })
            .collect()
    }
}

fn increment_counts(
    counts: &mut (u64, u64),
    accepted: bool,
    tries_counter: &'static str,
    accepts_counter: &'static str,
) -> Result<(), Fem1dError> {
    let tries = counts.1.checked_add(1).ok_or(Fem1dError::CounterOverflow {
        counter: tries_counter,
    })?;
    let accepts = if accepted {
        counts.0.checked_add(1).ok_or(Fem1dError::CounterOverflow {
            counter: accepts_counter,
        })?
    } else {
        counts.0
    };
    *counts = (accepts, tries);
    Ok(())
}

/// The hot-swap registry.
#[derive(Default)]
pub struct Registry {
    proposers: Vec<Box<dyn Proposer>>,
}

impl Registry {
    /// Empty registry.
    #[must_use]
    pub fn new() -> Self {
        Registry::default()
    }

    /// Register a proposer (consumers never change).
    ///
    /// # Errors
    /// Returns [`Fem1dError`] for an invalid/duplicate identity, a full
    /// synchronous registry, or allocation failure.
    pub fn register(&mut self, p: Box<dyn Proposer>) -> Result<(), Fem1dError> {
        if self.proposers.len() >= MAX_ZOO_PROPOSERS {
            return Err(Fem1dError::ResourceLimit {
                resource: "registered proposers",
                requested: self.proposers.len().saturating_add(1),
                limit: MAX_ZOO_PROPOSERS,
            });
        }
        let name = p.name();
        validate_identity(name, "proposer name")?;
        if self
            .proposers
            .iter()
            .any(|existing| existing.name() == name)
        {
            return Err(Fem1dError::InvalidScalar {
                field: "proposer name",
                reason: "must be unique within a registry",
            });
        }
        self.proposers
            .try_reserve_exact(1)
            .map_err(|_| Fem1dError::AllocationFailed {
                stage: "proposer registry",
                requested: self.proposers.len().saturating_add(1),
            })?;
        self.proposers.push(p);
        Ok(())
    }

    /// Deregister by name (independently retirable).
    pub fn deregister(&mut self, name: &str) {
        self.proposers.retain(|p| p.name() != name);
    }

    /// Registered names in order.
    #[must_use]
    pub fn names(&self) -> Vec<&'static str> {
        self.proposers.iter().map(|p| p.name()).collect()
    }

    /// Iterate registered proposers (the economics loop's view).
    pub fn proposers_dyn(&self) -> impl Iterator<Item = &dyn Proposer> {
        self.proposers.iter().map(AsRef::as_ref)
    }
}

/// Drive one speculation: gather proposals from enabled proposers,
/// order by ADVISORY confidence (descending; NaN sorts last;
/// deterministic name tie-break), and verify until one is accepted.
/// The ONLY path to a [`CertifiedAnswer`] is through the verifier.
///
/// # Errors
/// Returns [`Fem1dError`] before proposal work for an invalid query/resource
/// envelope, or when a built-in proposer cannot construct its candidate.
pub fn speculate(
    query: &SpeculationQuery,
    registry: &Registry,
    telemetry: &mut ZooTelemetry,
) -> Result<Outcome, Fem1dError> {
    validate_query(query)?;
    if registry.proposers.len() > MAX_ZOO_PROPOSERS {
        return Err(Fem1dError::ResourceLimit {
            resource: "registered proposers",
            requested: registry.proposers.len(),
            limit: MAX_ZOO_PROPOSERS,
        });
    }
    let mut proposals: Vec<(&'static str, Proposal)> = Vec::new();
    proposals
        .try_reserve_exact(registry.proposers.len())
        .map_err(|_| Fem1dError::AllocationFailed {
            stage: "proposal ordering",
            requested: registry.proposers.len(),
        })?;
    for p in &registry.proposers {
        let name = p.name();
        validate_identity(name, "proposer name")?;
        if telemetry.is_demoted(name, &query.regime) {
            continue;
        }
        if let Some(prop) = p.propose(query)? {
            if prop.candidate.len() > MAX_FEM1D_MESH_NODES {
                return Err(Fem1dError::ResourceLimit {
                    resource: "proposal candidate values",
                    requested: prop.candidate.len(),
                    limit: MAX_FEM1D_MESH_NODES,
                });
            }
            proposals.push((name, prop));
        }
    }
    if proposals.is_empty() {
        return Ok(Outcome::NoCandidates);
    }
    // Advisory ordering: confidence desc, NaN last, name tie-break.
    proposals.sort_by(|a, b| {
        let ca = if a.1.confidence.is_nan() {
            f64::NEG_INFINITY
        } else {
            a.1.confidence
        };
        let cb = if b.1.confidence.is_nan() {
            f64::NEG_INFINITY
        } else {
            b.1.confidence
        };
        cb.partial_cmp(&ca)
            .expect("NaN normalized")
            .then(a.0.cmp(b.0))
    });
    let mut attempted = Vec::new();
    attempted
        .try_reserve_exact(proposals.len())
        .map_err(|_| Fem1dError::AllocationFailed {
            stage: "rejected proposer identities",
            requested: proposals.len(),
        })?;
    let mut tried = 0u32;
    let mut best = None;
    for (name, prop) in proposals {
        tried = tried.checked_add(1).ok_or(Fem1dError::ResourceLimit {
            resource: "proposals tried",
            requested: usize::MAX,
            limit: u32::MAX as usize,
        })?;
        let report = verify(&query.problem, &prop.candidate, query.tolerance);
        let accepted = report.accept;
        telemetry.record(name, &query.regime, accepted)?;
        if accepted {
            return Ok(Outcome::Accepted(Box::new(CertifiedAnswer {
                candidate: prop.candidate,
                report,
                proposer: name,
                rejected_before: attempted,
            })));
        }
        attempted.push(name);
        let better = report.refusal.is_none()
            && report.bound.hi.is_finite()
            && best.as_ref().is_none_or(|best: &RejectedCandidate| {
                report.bound.hi < best.report.bound.hi
                    || (report.bound.hi.to_bits() == best.report.bound.hi.to_bits()
                        && name < best.proposer)
            });
        if better {
            best = Some(RejectedCandidate {
                proposer: name,
                candidate: prop.candidate,
                report,
            });
        }
    }
    Ok(Outcome::AllRejected(AllRejected {
        tried,
        attempted,
        best,
    }))
}

fn validate_query(query: &SpeculationQuery) -> Result<(), Fem1dError> {
    validate_problem(&query.problem)?;
    validate_finite_scalar(query.theta, "theta")?;
    validate_tolerance(query.tolerance)?;
    validate_identity(&query.regime, "regime")?;
    Ok(())
}

/// Emulate fp16 storage (10-bit mantissa truncation): the precision
/// discipline demo — speculate LOW, verify HIGH; the proposer's
/// precision is nobody's business.
#[must_use]
pub fn quantize_f16(x: f64) -> f64 {
    if !x.is_finite() {
        return x;
    }
    let bits = x.to_bits();
    // Keep sign + exponent + top 10 mantissa bits of the f64.
    let mask: u64 = !((1u64 << 42) - 1);
    f64::from_bits(bits & mask)
}

/// Proposer 1 — NEIGHBOR EXTRAPOLATION: retrieve the nearest CERTIFIED
/// run in design space; apply a first-order Taylor correction when a
/// cached sensitivity is available, degrade gracefully to zeroth-order
/// otherwise. Equidistant neighbors tie-break to the SMALLER θ
/// (deterministic).
pub struct NeighborExtrapolation {
    /// Certified prior runs: (θ, nodal solution, optional dU/dθ).
    pub cache: Vec<(f64, Vec<f64>, Option<Vec<f64>>)>,
}

impl Proposer for NeighborExtrapolation {
    fn name(&self) -> &'static str {
        "neighbor-extrapolation"
    }

    fn propose(&self, query: &SpeculationQuery) -> Result<Option<Proposal>, Fem1dError> {
        validate_query(query)?;
        if self.cache.len() > MAX_NEIGHBOR_CACHE_ENTRIES {
            return Err(Fem1dError::ResourceLimit {
                resource: "neighbor cache entries",
                requested: self.cache.len(),
                limit: MAX_NEIGHBOR_CACHE_ENTRIES,
            });
        }
        for (theta, values, sensitivity) in &self.cache {
            validate_finite_scalar(*theta, "neighbor cache theta")?;
            let distance = *theta - query.theta;
            if !distance.is_finite() {
                return Err(Fem1dError::NonFiniteIntermediate {
                    stage: "neighbor distance",
                    index: None,
                });
            }
            validate_candidate(&query.problem, values, "neighbor values")?;
            if let Some(sensitivity) = sensitivity {
                validate_candidate(&query.problem, sensitivity, "neighbor sensitivity")?;
            }
        }
        // Nearest by |θ − θ_i|; ties to the smaller θ.
        let best = self.cache.iter().min_by(|a, b| {
            let da = (a.0 - query.theta).abs();
            let db = (b.0 - query.theta).abs();
            da.total_cmp(&db).then(a.0.total_cmp(&b.0))
        });
        let Some(best) = best else {
            return Ok(None);
        };
        let (theta0, u0, sens) = best;
        let dt = query.theta - theta0;
        let mut candidate = Vec::new();
        let candidate_len = sens.as_ref().map_or(u0.len(), |du| u0.len().min(du.len()));
        candidate
            .try_reserve_exact(candidate_len)
            .map_err(|_| Fem1dError::AllocationFailed {
                stage: "neighbor candidate",
                requested: candidate_len,
            })?;
        match sens {
            Some(du) => candidate.extend(u0.iter().zip(du).map(|(u, d)| d.mul_add(dt, *u))),
            None => candidate.extend_from_slice(u0),
        }
        validate_candidate(&query.problem, &candidate, "neighbor candidate")?;
        Ok(Some(Proposal {
            candidate,
            // Advisory: nearer neighbors report higher confidence.
            confidence: 1.0 / (1.0 + dt.abs() * 10.0),
        }))
    }
}

/// Proposer 2 — COARSE-RUNG PROLONGATION: solve on the halved mesh
/// (rung k−1), prolongate linearly to the target mesh. Classical and
/// reliable; the asupersync speculative-race form (loser drained at a
/// tile boundary) is the CONTRACT no-claim.
pub struct CoarseRungProlongation;

impl Proposer for CoarseRungProlongation {
    fn name(&self) -> &'static str {
        "coarse-rung-prolongation"
    }

    fn propose(&self, query: &SpeculationQuery) -> Result<Option<Proposal>, Fem1dError> {
        validate_query(query)?;
        let mesh = query.problem.mesh();
        if mesh.len() < 5 {
            return Ok(None); // no coarser rung exists
        }
        // Coarse mesh: every other node (keeping the endpoints).
        let coarse_capacity = mesh.len() / 2 + 1;
        let mut coarse_mesh = Vec::new();
        coarse_mesh
            .try_reserve_exact(coarse_capacity)
            .map_err(|_| Fem1dError::AllocationFailed {
                stage: "coarse-rung mesh",
                requested: coarse_capacity,
            })?;
        coarse_mesh.extend(mesh.iter().step_by(2).copied());
        if mesh.len().is_multiple_of(2) {
            coarse_mesh.push(mesh[mesh.len() - 1]);
        }
        let coarse = query.problem.with_mesh(coarse_mesh)?;
        let cu = solve_p1(&coarse)?;
        // Linear prolongation onto the fine mesh.
        let mut candidate = Vec::new();
        candidate
            .try_reserve_exact(mesh.len())
            .map_err(|_| Fem1dError::AllocationFailed {
                stage: "coarse-rung prolongation",
                requested: mesh.len(),
            })?;
        let mut segment = 0usize;
        for (index, &x) in mesh.iter().enumerate() {
            while segment + 1 < coarse.mesh().len() - 1 && x > coarse.mesh()[segment + 1] {
                segment += 1;
            }
            let (x0, x1) = (coarse.mesh()[segment], coarse.mesh()[segment + 1]);
            if x < x0 || x > x1 {
                return Err(Fem1dError::NonFiniteIntermediate {
                    stage: "coarse-rung segment routing",
                    index: Some(index),
                });
            }
            let t = (x - x0) / (x1 - x0);
            let value = cu[segment] * (1.0 - t) + cu[segment + 1] * t;
            if !t.is_finite() || !value.is_finite() {
                return Err(Fem1dError::NonFiniteIntermediate {
                    stage: "coarse-rung prolongation",
                    index: Some(index),
                });
            }
            candidate.push(value);
        }
        Ok(Some(Proposal {
            candidate,
            confidence: 0.7,
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn telemetry_keys_and_counters_fail_closed() {
        let mut telemetry = ZooTelemetry::default();
        for index in 0..MAX_ZOO_TELEMETRY_KEYS {
            telemetry
                .record("p", &format!("r{index}"), false)
                .expect("key within telemetry cap");
        }
        assert!(matches!(
            telemetry.record("p", "overflow", false),
            Err(Fem1dError::ResourceLimit {
                resource: "zoo telemetry keys",
                ..
            })
        ));

        let mut overflow = ZooTelemetry::default();
        overflow
            .counts
            .insert(("p".to_string(), "r".to_string()), (0, u64::MAX));
        assert!(matches!(
            overflow.record("p", "r", false),
            Err(Fem1dError::CounterOverflow {
                counter: "zoo tries"
            })
        ));
        let mut accept_overflow = ZooTelemetry::default();
        accept_overflow
            .counts
            .insert(("p".to_string(), "r".to_string()), (u64::MAX, 7));
        assert!(matches!(
            accept_overflow.record("p", "r", true),
            Err(Fem1dError::CounterOverflow {
                counter: "zoo accepts"
            })
        ));
        assert_eq!(
            accept_overflow
                .counts
                .get(&("p".to_string(), "r".to_string())),
            Some(&(u64::MAX, 7)),
            "a failed counter update must not partially increment tries"
        );
        assert!(matches!(
            overflow.demote_collapsed(f64::NAN, 1),
            Err(Fem1dError::InvalidScalar {
                field: "zoo demotion threshold",
                ..
            })
        ));
    }
}
