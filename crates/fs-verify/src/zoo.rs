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
use crate::fem1d::{MmsProblem, solve_p1};
use std::collections::BTreeMap;
use std::fmt::Write as _;

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
    fn propose(&self, query: &SpeculationQuery) -> Option<Proposal>;
}

/// A CERTIFIED answer: candidate + the verifier's report (verified
/// color included). NO public constructor — the type is the proof
/// that the verifier said yes.
#[derive(Debug, Clone)]
pub struct CertifiedAnswer {
    candidate: Vec<f64>,
    report: VerifierReport,
    proposer: &'static str,
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
}

/// The speculation outcome.
#[derive(Debug, Clone)]
pub enum Outcome {
    /// A proposal passed the verifier.
    Accepted(Box<CertifiedAnswer>),
    /// Proposals were tried; every one was rejected (fall back to the
    /// full solve — nothing was corrupted, only checks were spent).
    AllRejected {
        /// How many candidates were checked.
        tried: u32,
    },
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
    fn record(&mut self, proposer: &str, regime: &str, accepted: bool) {
        let e = self
            .counts
            .entry((proposer.to_string(), regime.to_string()))
            .or_insert((0, 0));
        e.1 += 1;
        if accepted {
            e.0 += 1;
        }
    }

    /// Accept rate for (proposer, regime).
    #[must_use]
    pub fn accept_rate(&self, proposer: &str, regime: &str) -> Option<f64> {
        self.counts
            .get(&(proposer.to_string(), regime.to_string()))
            .map(|&(a, t)| a as f64 / t.max(1) as f64)
    }

    /// Demote proposers whose accept rate in a regime collapsed below
    /// `threshold` after at least `min_tries` attempts. Returns the
    /// demotions performed.
    pub fn demote_collapsed(&mut self, threshold: f64, min_tries: u64) -> Vec<(String, String)> {
        let mut out = Vec::new();
        for ((p, r), &(a, t)) in &self.counts {
            if t >= min_tries && (a as f64 / t as f64) < threshold {
                let key = (p.clone(), r.clone());
                if !self.demoted.get(&key).copied().unwrap_or(false) {
                    self.demoted.insert(key.clone(), true);
                    out.push(key);
                }
            }
        }
        out
    }

    /// Is a proposer demoted in a regime?
    #[must_use]
    pub fn is_demoted(&self, proposer: &str, regime: &str) -> bool {
        self.demoted
            .get(&(proposer.to_string(), regime.to_string()))
            .copied()
            .unwrap_or(false)
    }

    /// Ledger rows (per proposer × regime).
    #[must_use]
    pub fn rows(&self) -> Vec<String> {
        self.counts
            .iter()
            .map(|((p, r), &(a, t))| {
                let mut s = String::new();
                let _ = write!(
                    s,
                    "{{\"proposer\":\"{p}\",\"regime\":\"{r}\",\"accepts\":{a},\
                     \"tries\":{t},\"rate\":{:.4},\"demoted\":{}}}",
                    a as f64 / t.max(1) as f64,
                    self.is_demoted(p, r)
                );
                s
            })
            .collect()
    }
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
    pub fn register(&mut self, p: Box<dyn Proposer>) {
        self.proposers.push(p);
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
#[must_use]
pub fn speculate(
    query: &SpeculationQuery,
    registry: &Registry,
    telemetry: &mut ZooTelemetry,
) -> Outcome {
    let mut proposals: Vec<(&'static str, Proposal)> = Vec::new();
    for p in &registry.proposers {
        if telemetry.is_demoted(p.name(), &query.regime) {
            continue;
        }
        if let Some(prop) = p.propose(query) {
            proposals.push((p.name(), prop));
        }
    }
    if proposals.is_empty() {
        return Outcome::NoCandidates;
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
    let mut tried = 0;
    for (name, prop) in proposals {
        tried += 1;
        let report = verify(&query.problem, &prop.candidate, query.tolerance);
        let accepted = report.accept;
        telemetry.record(name, &query.regime, accepted);
        if accepted {
            return Outcome::Accepted(Box::new(CertifiedAnswer {
                candidate: prop.candidate,
                report,
                proposer: name,
            }));
        }
    }
    Outcome::AllRejected { tried }
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

    fn propose(&self, query: &SpeculationQuery) -> Option<Proposal> {
        // Nearest by |θ − θ_i|; ties to the smaller θ.
        let best = self.cache.iter().min_by(|a, b| {
            let da = (a.0 - query.theta).abs();
            let db = (b.0 - query.theta).abs();
            da.partial_cmp(&db)
                .expect("finite thetas")
                .then(a.0.partial_cmp(&b.0).expect("finite thetas"))
        })?;
        let (theta0, u0, sens) = best;
        let dt = query.theta - theta0;
        let candidate: Vec<f64> = match sens {
            Some(du) => u0.iter().zip(du).map(|(u, d)| d.mul_add(dt, *u)).collect(),
            None => u0.clone(),
        };
        Some(Proposal {
            candidate,
            // Advisory: nearer neighbors report higher confidence.
            confidence: 1.0 / (1.0 + dt.abs() * 10.0),
        })
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

    fn propose(&self, query: &SpeculationQuery) -> Option<Proposal> {
        let mesh = &query.problem.mesh;
        if mesh.len() < 5 {
            return None; // no coarser rung exists
        }
        // Coarse mesh: every other node (keeping the endpoints).
        let coarse_mesh: Vec<f64> = mesh
            .iter()
            .step_by(2)
            .copied()
            .chain(if mesh.len().is_multiple_of(2) {
                Some(*mesh.last().expect("nonempty"))
            } else {
                None
            })
            .collect();
        let coarse = MmsProblem::new(&query.problem.name, query.problem.u.clone(), coarse_mesh);
        let cu = solve_p1(&coarse);
        // Linear prolongation onto the fine mesh.
        let mut candidate = Vec::with_capacity(mesh.len());
        for &x in mesh {
            let seg = coarse
                .mesh
                .windows(2)
                .position(|w| x >= w[0] && x <= w[1])
                .unwrap_or(coarse.mesh.len() - 2);
            let (x0, x1) = (coarse.mesh[seg], coarse.mesh[seg + 1]);
            let t = if x1 > x0 { (x - x0) / (x1 - x0) } else { 0.0 };
            candidate.push(cu[seg] * (1.0 - t) + cu[seg + 1] * t);
        }
        Some(Proposal {
            candidate,
            confidence: 0.7,
        })
    }
}
