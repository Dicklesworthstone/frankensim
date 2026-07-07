//! fs-constraint — the constraint CALCULUS (plan §9.1, patch Rev F).
//! Layer: L4.
//!
//! Constraints with SEMANTICS, not anonymous `g(x) ≤ 0`:
//!
//! - [`ConstraintKind`] taxonomy — Hard (never traded), Soft (penalty
//!   law), Chance (estimator + validity machinery: a chance constraint
//!   reports satisfied only when its CONFIDENCE BOUND clears the level,
//!   not when the raw empirical rate does), Robust (uncertainty box,
//!   proven by interval evaluation), Certification (REFUSES to report
//!   satisfied without a proof artifact; the in-house interval prover
//!   is a real one), Fabrication and Code (domain semantics carried to
//!   the ledger). Each kind maps to an optimizer [`Treatment`].
//! - [`ConstraintEvidence`] per evaluation: status, EVIDENCE-TYPED
//!   violation magnitude (fs-evidence certificates), active-set role,
//!   and ranked repair suggestions.
//! - Infeasibility DIAGNOSIS ([`crate::diagnose`]): elastic-relaxation
//!   solves find where a design space fights itself; deletion filtering
//!   extracts a MINIMAL unsat core (dropping any member restores
//!   feasibility); repairs come RANKED with feasibility estimates
//!   calibrated against enumeration — optimizer failures become design
//!   conversations.
//!
//! fs-opt hosts the expression graphs; this crate owns what the
//! constraints MEAN.

mod diagnose;
mod ival;

pub use diagnose::{
    Diagnosis, DomainBox, ElasticReport, RepairAction, RepairKind, diagnose_infeasibility,
    elastic_solve,
};
pub use ival::{Iv, IvalError, interval_eval};

use fs_evidence::{NumericalCertificate, StatisticalCertificate};
use fs_opt::{NodeId, Problem};

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// How a chance constraint's probability is estimated.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ChanceEstimator {
    /// Deterministic-stream Monte Carlo with a Hoeffding confidence
    /// bound at level `delta` (the validity machinery: SATISFIED means
    /// the LOWER bound clears the target, not the raw rate).
    MonteCarlo {
        /// Sample count.
        samples: u32,
        /// Bound failure probability (e.g. 0.05).
        delta: f64,
    },
}

/// Which proof a certification constraint demands.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProofKind {
    /// Rigorous interval evaluation over the stated domain (in-house,
    /// available NOW via [`interval_eval`]).
    Interval,
    /// Sum-of-squares certificate (fs-sos; represented, not yet
    /// executable — a CONTRACT no-claim).
    Sos,
}

/// A proof artifact attached to a certification constraint.
#[derive(Debug, Clone, PartialEq)]
pub enum ProofArtifact {
    /// Interval proof: `hi ≤ 0` over the domain box (carried bound).
    IntervalBound {
        /// The proven upper bound of `g` over the domain.
        hi: f64,
    },
    /// An external SOS certificate reference (opaque until fs-sos).
    SosReference {
        /// Artifact identifier.
        id: String,
    },
}

/// Soft-constraint penalty laws.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PenaltyLaw {
    /// `w · max(g, 0)²`.
    Quadratic {
        /// Weight.
        weight: f64,
    },
    /// `w · max(g, 0)`.
    Hinge {
        /// Weight.
        weight: f64,
    },
}

/// The constraint taxonomy (semantics, not just inequalities).
#[derive(Debug, Clone, PartialEq)]
pub enum ConstraintKind {
    /// Physics/safety: never traded away; violations demand
    /// feasibility restoration before anything else.
    Hard,
    /// Preference: violations price into the objective by law.
    Soft(PenaltyLaw),
    /// `P(g ≤ 0) ≥ level` under a declared noise model, estimated by a
    /// declared estimator WITH validity machinery.
    Chance {
        /// Required probability level.
        level: f64,
        /// The estimator (and its validity parameters).
        estimator: ChanceEstimator,
    },
    /// `g ≤ 0` for ALL parameter draws in an uncertainty box around
    /// the design point (proven conservatively by interval evaluation).
    Robust {
        /// Half-widths of the uncertainty box per variable component.
        half_widths: Vec<f64>,
    },
    /// Requires a PROOF artifact; refuses "satisfied" without one.
    Certification {
        /// The proof kind demanded.
        proof: ProofKind,
    },
    /// Manufacturability semantics (process named for the ledger;
    /// fs-fab supplies process models).
    Fabrication {
        /// Process name (e.g. "cnc-3axis").
        process: String,
    },
    /// Design-code semantics (standard named for the ledger; the frame
    /// flagship's AISC-class checks).
    Code {
        /// Standard name (e.g. "aisc-360").
        standard: String,
    },
}

/// How an optimizer must treat a kind (routing metadata).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Treatment {
    /// Restore feasibility with priority; never trade.
    FeasibilityRestoration,
    /// Fold into the objective via the penalty law.
    PenaltyTerm,
    /// Estimate with the declared estimator, act on the BOUND.
    EstimateThenBound,
    /// Prove over the set (interval/SOS) or escalate.
    ProveOrEscalate,
    /// Evaluate the domain rule; report to the ledger.
    DomainCheck,
}

impl ConstraintKind {
    /// The optimizer treatment this kind demands.
    #[must_use]
    pub fn treatment(&self) -> Treatment {
        match self {
            ConstraintKind::Hard => Treatment::FeasibilityRestoration,
            ConstraintKind::Soft(_) => Treatment::PenaltyTerm,
            ConstraintKind::Chance { .. } => Treatment::EstimateThenBound,
            ConstraintKind::Robust { .. } | ConstraintKind::Certification { .. } => {
                Treatment::ProveOrEscalate
            }
            ConstraintKind::Fabrication { .. } | ConstraintKind::Code { .. } => {
                Treatment::DomainCheck
            }
        }
    }

    /// Stable kind name for ledger rows.
    #[must_use]
    pub fn kind_name(&self) -> &'static str {
        match self {
            ConstraintKind::Hard => "hard",
            ConstraintKind::Soft(_) => "soft",
            ConstraintKind::Chance { .. } => "chance",
            ConstraintKind::Robust { .. } => "robust",
            ConstraintKind::Certification { .. } => "certification",
            ConstraintKind::Fabrication { .. } => "fabrication",
            ConstraintKind::Code { .. } => "code",
        }
    }
}

/// One typed constraint over an fs-opt graph node (`g ≤ 0` semantics;
/// the node must be scalar in the host problem).
#[derive(Debug, Clone, PartialEq)]
pub struct ConstraintSpec {
    /// Human name (diagnostics, ledger).
    pub name: String,
    /// The `g` node in the host problem.
    pub node: NodeId,
    /// Semantics.
    pub kind: ConstraintKind,
    /// Active-set tolerance (|g| ≤ tol counts as active).
    pub active_tol: f64,
}

/// Evaluation status of one constraint.
#[derive(Debug, Clone, PartialEq)]
pub enum Status {
    /// `g ≤ −tol` (strictly inside).
    Satisfied,
    /// `|g| ≤ tol` (on the boundary — active).
    Active,
    /// `g > tol`.
    Violated,
    /// Certification kind without its proof artifact: NOT allowed to
    /// claim satisfied (the refusal is the feature).
    NeedsProof {
        /// What proof is missing.
        proof: ProofKind,
    },
    /// Proven over the stated set (interval/SOS artifact attached).
    Proven,
    /// Chance kind: the confidence BOUND does not clear the level even
    /// if the raw estimate does (validity machinery speaking).
    BoundNotCleared {
        /// Empirical satisfaction rate.
        empirical: f64,
        /// The lower confidence bound that failed to clear.
        lower_bound: f64,
    },
}

/// Active-set role for optimizer consumption.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ActiveRole {
    /// Slack: locally ignorable.
    Inactive,
    /// On the boundary: shapes the local geometry.
    Active,
    /// Currently violated: drives restoration.
    Violating,
}

/// The per-evaluation artifact: status + evidence-typed magnitude +
/// role (+ repairs, attached by the diagnosis pass when violated).
#[derive(Debug, Clone, PartialEq)]
pub struct ConstraintEvidence {
    /// Which constraint.
    pub name: String,
    /// Kind name (ledger).
    pub kind: &'static str,
    /// Status.
    pub status: Status,
    /// Violation magnitude `max(g, 0)` with its numerical certificate
    /// (exact for algebraic graphs).
    pub violation: f64,
    /// Certificate for the violation value.
    pub certificate: NumericalCertificate,
    /// Statistical certificate (chance kinds).
    pub statistical: StatisticalCertificate,
    /// Active-set role.
    pub role: ActiveRole,
    /// Soft-kind penalty contribution (0 otherwise).
    pub penalty: f64,
}

impl ConstraintEvidence {
    /// Canonical ledger row (Rev S table shape).
    #[must_use]
    pub fn to_ledger_row(&self) -> String {
        let status = match &self.status {
            Status::Satisfied => "satisfied".to_string(),
            Status::Active => "active".to_string(),
            Status::Violated => "violated".to_string(),
            Status::NeedsProof { proof } => format!("needs-proof:{proof:?}"),
            Status::Proven => "proven".to_string(),
            Status::BoundNotCleared { .. } => "bound-not-cleared".to_string(),
        };
        format!(
            "{{\"constraint\":\"{}\",\"kind\":\"{}\",\"status\":\"{}\",\
             \"violation\":{:.6e},\"penalty\":{:.6e}}}",
            self.name, self.kind, status, self.violation, self.penalty
        )
    }
}

/// Errors this crate teaches with.
#[derive(Debug, Clone, PartialEq)]
pub enum ConError {
    /// The referenced node is not scalar in the host problem.
    NotScalar {
        /// Node id.
        node: u32,
    },
    /// Underlying fs-opt evaluation failed (carried through).
    Eval(fs_opt::OptError),
    /// Interval proof attempt failed (with the reason).
    NotProvable {
        /// Why the interval engine refused.
        why: String,
    },
    /// A parameter left its valid range.
    BadParam {
        /// What.
        what: &'static str,
        /// Value.
        value: f64,
    },
    /// Serialized text failed to parse.
    Parse {
        /// 1-based line.
        line: usize,
        /// What went wrong.
        what: String,
    },
}

impl core::fmt::Display for ConError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            ConError::NotScalar { node } => write!(
                f,
                "constraint node {node} is not scalar; reduce it (dot/norm_sq/component) first"
            ),
            ConError::Eval(e) => write!(f, "evaluation failed: {e}"),
            ConError::NotProvable { why } => write!(
                f,
                "interval proof refused: {why}; tighten the domain, rewrite the \
                 expression, or escalate to an SOS certificate"
            ),
            ConError::BadParam { what, value } => {
                write!(f, "`{what}` = {value} is outside its valid range")
            }
            ConError::Parse { line, what } => write!(f, "parse error at line {line}: {what}"),
        }
    }
}

impl std::error::Error for ConError {}

impl From<fs_opt::OptError> for ConError {
    fn from(e: fs_opt::OptError) -> Self {
        ConError::Eval(e)
    }
}

pub(crate) fn scalar_at(problem: &Problem, node: NodeId, x: &[f64]) -> Result<f64, ConError> {
    let v = fs_opt::eval(problem, node, std::slice::from_ref(&x.to_vec()))?;
    v.scalar().ok_or(ConError::NotScalar { node: node.0 })
}

/// Evaluate one constraint at a design point (single Rn-variable
/// problems — the v1 host shape). `noise` supplies parameter draws for
/// chance kinds (deterministic streams in tests).
///
/// # Errors
/// Teaching errors ([`ConError`]).
#[allow(clippy::too_many_lines)] // one arm per kind: the calculus IS the dispatch
pub fn evaluate(
    problem: &Problem,
    spec: &ConstraintSpec,
    x: &[f64],
    noise: Option<&dyn Fn(u64) -> Vec<f64>>,
) -> Result<ConstraintEvidence, ConError> {
    let g = scalar_at(problem, spec.node, x)?;
    let violation = g.max(0.0);
    let base_status = if g > spec.active_tol {
        Status::Violated
    } else if g >= -spec.active_tol {
        Status::Active
    } else {
        Status::Satisfied
    };
    let role = match base_status {
        Status::Violated => ActiveRole::Violating,
        Status::Active => ActiveRole::Active,
        _ => ActiveRole::Inactive,
    };
    let mut ev = ConstraintEvidence {
        name: spec.name.clone(),
        kind: spec.kind.kind_name(),
        status: base_status,
        violation,
        certificate: NumericalCertificate::exact(violation),
        statistical: StatisticalCertificate::None,
        role,
        penalty: 0.0,
    };
    match &spec.kind {
        ConstraintKind::Hard | ConstraintKind::Fabrication { .. } | ConstraintKind::Code { .. } => {
        }
        ConstraintKind::Soft(law) => {
            ev.penalty = match *law {
                PenaltyLaw::Quadratic { weight } => weight * violation * violation,
                PenaltyLaw::Hinge { weight } => weight * violation,
            };
        }
        ConstraintKind::Chance { level, estimator } => {
            let ChanceEstimator::MonteCarlo { samples, delta } = *estimator;
            if !(*level > 0.0 && *level < 1.0) {
                return Err(ConError::BadParam {
                    what: "chance level",
                    value: *level,
                });
            }
            let noise = noise.ok_or(ConError::BadParam {
                what: "chance noise model (required)",
                value: f64::NAN,
            })?;
            let mut hits = 0u32;
            for s in 0..samples {
                let draw = noise(u64::from(s));
                let shifted: Vec<f64> = x.iter().zip(&draw).map(|(a, b)| a + b).collect();
                if scalar_at(problem, spec.node, &shifted)? <= 0.0 {
                    hits += 1;
                }
            }
            let empirical = f64::from(hits) / f64::from(samples);
            // Hoeffding lower confidence bound at failure prob delta.
            let half_width = ((1.0 / delta).ln() / (2.0 * f64::from(samples))).sqrt();
            let lower = empirical - half_width;
            ev.statistical = StatisticalCertificate::HalfWidth {
                half_width,
                confidence: 1.0 - delta,
            };
            ev.status = if lower >= *level {
                Status::Satisfied
            } else if empirical >= *level {
                // The raw rate clears but the BOUND does not: refuse —
                // this is the validity machinery earning its keep.
                Status::BoundNotCleared {
                    empirical,
                    lower_bound: lower,
                }
            } else {
                Status::Violated
            };
            ev.violation = (*level - lower).max(0.0);
            ev.certificate = NumericalCertificate::estimate(ev.violation, ev.violation);
            ev.role = if matches!(ev.status, Status::Satisfied) {
                ActiveRole::Inactive
            } else {
                ActiveRole::Violating
            };
        }
        ConstraintKind::Robust { half_widths } => {
            // Prove sup over the uncertainty box via interval eval.
            let boxes: Vec<(f64, f64)> = x
                .iter()
                .zip(half_widths)
                .map(|(c, h)| (c - h, c + h))
                .collect();
            match interval_eval(problem, spec.node, &boxes) {
                Ok(iv) => {
                    ev.violation = iv.hi.max(0.0);
                    ev.certificate = NumericalCertificate::enclosure(iv.lo, iv.hi);
                    ev.status = if iv.hi <= 0.0 {
                        Status::Proven
                    } else {
                        Status::Violated
                    };
                    ev.role = if iv.hi <= 0.0 {
                        ActiveRole::Inactive
                    } else {
                        ActiveRole::Violating
                    };
                }
                Err(e) => {
                    return Err(ConError::NotProvable { why: e.to_string() });
                }
            }
        }
        ConstraintKind::Certification { proof } => {
            // Without an artifact the status is NeedsProof — REGARDLESS
            // of how good g(x) looks pointwise.
            ev.status = Status::NeedsProof { proof: *proof };
            ev.role = ActiveRole::Violating;
        }
    }
    Ok(ev)
}

/// Attempt the interval proof for a certification constraint over a
/// stated domain box; success attaches the artifact and the PROVEN
/// status.
///
/// # Errors
/// [`ConError::NotProvable`] with the engine's reason (an honest gap,
/// not a failure).
pub fn prove_interval(
    problem: &Problem,
    spec: &ConstraintSpec,
    domain: &[(f64, f64)],
) -> Result<(ConstraintEvidence, ProofArtifact), ConError> {
    let iv = interval_eval(problem, spec.node, domain)
        .map_err(|e| ConError::NotProvable { why: e.to_string() })?;
    if iv.hi <= 0.0 {
        Ok((
            ConstraintEvidence {
                name: spec.name.clone(),
                kind: spec.kind.kind_name(),
                status: Status::Proven,
                violation: 0.0,
                certificate: NumericalCertificate::enclosure(iv.lo, iv.hi),
                statistical: StatisticalCertificate::None,
                role: ActiveRole::Inactive,
                penalty: 0.0,
            },
            ProofArtifact::IntervalBound { hi: iv.hi },
        ))
    } else {
        Err(ConError::NotProvable {
            why: format!(
                "interval bound over the domain is [{:.3e}, {:.3e}]; the upper end \
                 exceeds 0",
                iv.lo, iv.hi
            ),
        })
    }
}

/// Serialize a constraint set (canonical line form; floats as bits).
#[must_use]
pub fn serialize_specs(specs: &[ConstraintSpec]) -> String {
    use std::fmt::Write as _;
    let hex = |v: f64| format!("{:016X}", v.to_bits());
    let mut s = String::from("fscon v1\n");
    for c in specs {
        let kind = match &c.kind {
            ConstraintKind::Hard => "hard".to_string(),
            ConstraintKind::Soft(PenaltyLaw::Quadratic { weight }) => {
                format!("soft quadratic {}", hex(*weight))
            }
            ConstraintKind::Soft(PenaltyLaw::Hinge { weight }) => {
                format!("soft hinge {}", hex(*weight))
            }
            ConstraintKind::Chance {
                level,
                estimator: ChanceEstimator::MonteCarlo { samples, delta },
            } => format!("chance {} mc {samples} {}", hex(*level), hex(*delta)),
            ConstraintKind::Robust { half_widths } => {
                let ws: Vec<String> = half_widths.iter().map(|w| hex(*w)).collect();
                format!("robust {}", ws.join(","))
            }
            ConstraintKind::Certification { proof } => match proof {
                ProofKind::Interval => "certification interval".to_string(),
                ProofKind::Sos => "certification sos".to_string(),
            },
            ConstraintKind::Fabrication { process } => {
                format!("fabrication {}", process.replace(' ', "%20"))
            }
            ConstraintKind::Code { standard } => {
                format!("code {}", standard.replace(' ', "%20"))
            }
        };
        let _ = writeln!(
            s,
            "constraint {} {} {} {kind}",
            c.name.replace(' ', "%20"),
            c.node.0,
            hex(c.active_tol)
        );
    }
    s
}

/// Parse [`serialize_specs`] output (round-trip identity).
///
/// # Errors
/// [`ConError::Parse`] with line numbers.
#[allow(clippy::too_many_lines)] // one grammar rule per kind
pub fn parse_specs(text: &str) -> Result<Vec<ConstraintSpec>, ConError> {
    let unhex = |s: &str| -> Option<f64> { u64::from_str_radix(s, 16).ok().map(f64::from_bits) };
    let perr = |line: usize, what: &str| ConError::Parse {
        line,
        what: what.to_string(),
    };
    let mut out = Vec::new();
    for (ln0, line) in text.lines().enumerate() {
        let ln = ln0 + 1;
        let toks: Vec<&str> = line.split(' ').collect();
        match toks.first().copied() {
            Some("fscon") => {
                if toks.get(1) != Some(&"v1") {
                    return Err(perr(ln, "unsupported version"));
                }
            }
            Some("constraint") => {
                let name = toks
                    .get(1)
                    .ok_or_else(|| perr(ln, "missing name"))?
                    .replace("%20", " ");
                let node: u32 = toks
                    .get(2)
                    .and_then(|t| t.parse().ok())
                    .ok_or_else(|| perr(ln, "bad node"))?;
                let active_tol = toks
                    .get(3)
                    .and_then(|t| unhex(t))
                    .ok_or_else(|| perr(ln, "bad tol"))?;
                let kind = match toks.get(4).copied() {
                    Some("hard") => ConstraintKind::Hard,
                    Some("soft") => {
                        let w = toks
                            .get(6)
                            .and_then(|t| unhex(t))
                            .ok_or_else(|| perr(ln, "bad weight"))?;
                        match toks.get(5).copied() {
                            Some("quadratic") => {
                                ConstraintKind::Soft(PenaltyLaw::Quadratic { weight: w })
                            }
                            Some("hinge") => ConstraintKind::Soft(PenaltyLaw::Hinge { weight: w }),
                            _ => return Err(perr(ln, "unknown penalty law")),
                        }
                    }
                    Some("chance") => {
                        let level = toks
                            .get(5)
                            .and_then(|t| unhex(t))
                            .ok_or_else(|| perr(ln, "bad level"))?;
                        if toks.get(6) != Some(&"mc") {
                            return Err(perr(ln, "unknown estimator"));
                        }
                        let samples = toks
                            .get(7)
                            .and_then(|t| t.parse().ok())
                            .ok_or_else(|| perr(ln, "bad samples"))?;
                        let delta = toks
                            .get(8)
                            .and_then(|t| unhex(t))
                            .ok_or_else(|| perr(ln, "bad delta"))?;
                        ConstraintKind::Chance {
                            level,
                            estimator: ChanceEstimator::MonteCarlo { samples, delta },
                        }
                    }
                    Some("robust") => {
                        let ws = toks.get(5).ok_or_else(|| perr(ln, "missing widths"))?;
                        let half_widths: Option<Vec<f64>> = ws.split(',').map(unhex).collect();
                        ConstraintKind::Robust {
                            half_widths: half_widths.ok_or_else(|| perr(ln, "bad widths"))?,
                        }
                    }
                    Some("certification") => match toks.get(5).copied() {
                        Some("interval") => ConstraintKind::Certification {
                            proof: ProofKind::Interval,
                        },
                        Some("sos") => ConstraintKind::Certification {
                            proof: ProofKind::Sos,
                        },
                        _ => return Err(perr(ln, "unknown proof kind")),
                    },
                    Some("fabrication") => ConstraintKind::Fabrication {
                        process: toks
                            .get(5)
                            .ok_or_else(|| perr(ln, "missing process"))?
                            .replace("%20", " "),
                    },
                    Some("code") => ConstraintKind::Code {
                        standard: toks
                            .get(5)
                            .ok_or_else(|| perr(ln, "missing standard"))?
                            .replace("%20", " "),
                    },
                    _ => return Err(perr(ln, "unknown kind")),
                };
                out.push(ConstraintSpec {
                    name,
                    node: NodeId(node),
                    kind,
                    active_tol,
                });
            }
            Some("") | None => {}
            Some(other) => {
                return Err(ConError::Parse {
                    line: ln,
                    what: format!("unknown directive `{other}`"),
                });
            }
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn version_is_stamped() {
        assert!(!VERSION.is_empty());
    }

    #[test]
    fn treatments_map_per_kind() {
        assert_eq!(
            ConstraintKind::Hard.treatment(),
            Treatment::FeasibilityRestoration
        );
        assert_eq!(
            ConstraintKind::Soft(PenaltyLaw::Hinge { weight: 1.0 }).treatment(),
            Treatment::PenaltyTerm
        );
        assert_eq!(
            ConstraintKind::Chance {
                level: 0.9,
                estimator: ChanceEstimator::MonteCarlo {
                    samples: 64,
                    delta: 0.05
                }
            }
            .treatment(),
            Treatment::EstimateThenBound
        );
        assert_eq!(
            ConstraintKind::Certification {
                proof: ProofKind::Interval
            }
            .treatment(),
            Treatment::ProveOrEscalate
        );
        assert_eq!(
            ConstraintKind::Code {
                standard: "aisc-360".into()
            }
            .treatment(),
            Treatment::DomainCheck
        );
    }
}
