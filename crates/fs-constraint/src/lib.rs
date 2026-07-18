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
pub use ival::{Iv, IvalError};

use fs_evidence::{NumericalCertificate, NumericalKind, StatisticalCertificate};
use fs_opt::{Manifold, NodeId, Problem, ProblemSemanticId};

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

pub(crate) fn push_json_string(out: &mut String, value: &str) {
    use core::fmt::Write as _;

    out.push('"');
    for ch in value.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\u{0008}' => out.push_str("\\b"),
            '\u{000c}' => out.push_str("\\f"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            control if control <= '\u{001f}' => {
                let _ = write!(out, "\\u{:04x}", u32::from(control));
            }
            printable => out.push(printable),
        }
    }
    out.push('"');
}

fn wire_token_is_safe(byte: u8) -> bool {
    byte.is_ascii_alphanumeric() || matches!(byte, b'-' | b'.' | b'_' | b'~')
}

fn encode_wire_token(value: &str) -> String {
    use core::fmt::Write as _;

    // A lone percent is the empty-token sentinel. A literal percent is always
    // `%25`, so the sentinel cannot collide with caller text.
    if value.is_empty() {
        return "%".to_string();
    }
    let mut encoded = String::with_capacity(value.len());
    for byte in value.bytes() {
        if wire_token_is_safe(byte) {
            encoded.push(char::from(byte));
        } else {
            let _ = write!(encoded, "%{byte:02X}");
        }
    }
    encoded
}

fn decode_wire_nibble(byte: u8) -> Option<u8> {
    match byte {
        b'0'..=b'9' => Some(byte - b'0'),
        b'A'..=b'F' => Some(byte - b'A' + 10),
        _ => None,
    }
}

fn decode_wire_token(token: &str, line: usize, field: &'static str) -> Result<String, ConError> {
    if token == "%" {
        return Ok(String::new());
    }
    let bytes = token.as_bytes();
    let mut decoded = Vec::with_capacity(bytes.len());
    let mut index = 0;
    while index < bytes.len() {
        let byte = bytes[index];
        if byte == b'%' {
            let Some((&hi, &lo)) = bytes.get(index + 1).zip(bytes.get(index + 2)) else {
                return Err(ConError::Parse {
                    line,
                    what: format!("truncated percent escape in {field}"),
                });
            };
            let value = decode_wire_nibble(hi)
                .zip(decode_wire_nibble(lo))
                .map(|(hi, lo)| (hi << 4) | lo)
                .ok_or_else(|| ConError::Parse {
                    line,
                    what: format!("noncanonical percent escape in {field}"),
                })?;
            decoded.push(value);
            index += 3;
        } else if wire_token_is_safe(byte) {
            decoded.push(byte);
            index += 1;
        } else {
            return Err(ConError::Parse {
                line,
                what: format!("unescaped byte in {field}"),
            });
        }
    }
    let decoded = String::from_utf8(decoded).map_err(|_| ConError::Parse {
        line,
        what: format!("{field} is not valid UTF-8"),
    })?;
    if encode_wire_token(&decoded) != token {
        return Err(ConError::Parse {
            line,
            what: format!("noncanonical encoding of {field}"),
        });
    }
    Ok(decoded)
}

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

/// Exact subject identity for a proof artifact.
///
/// The problem component is the admitted, full-width semantic identity from
/// `fs-opt`; the node and every domain endpoint remain exact rather than being
/// collapsed into a local short hash. Fields are sealed so callers cannot mint
/// an apparently admitted subject by struct literal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProofSubject {
    problem: ProblemSemanticId,
    node: NodeId,
    domain_bits: Vec<(u64, u64)>,
}

impl ProofSubject {
    /// Admitted problem semantic identity.
    #[must_use]
    pub fn problem(&self) -> ProblemSemanticId {
        self.problem
    }

    /// Exact scalar expression node proved.
    #[must_use]
    pub fn node(&self) -> NodeId {
        self.node
    }

    /// Exact `(lo_bits, hi_bits)` sequence of the proved domain.
    #[must_use]
    pub fn domain_bits(&self) -> &[(u64, u64)] {
        &self.domain_bits
    }

    fn push_json(&self, out: &mut String) {
        use core::fmt::Write as _;

        out.push_str("{\"problem\":");
        push_json_string(out, &self.problem.to_hex());
        let _ = write!(out, ",\"node\":{},\"domain_bits\":[", self.node.0);
        for (index, &(lo, hi)) in self.domain_bits.iter().enumerate() {
            if index > 0 {
                out.push(',');
            }
            let _ = write!(out, "[\"{lo:016X}\",\"{hi:016X}\"]");
        }
        out.push_str("]}");
    }
}

/// A sealed interval-proof artifact attached to a certification constraint.
///
/// Only [`prove_interval`] mints this type. It binds both ends of the computed
/// enclosure to the admitted problem, exact node, and exact domain bytes; an
/// artifact from another constraint or box therefore cannot be confused with
/// this subject by a consumer that checks [`Self::is_bound_to`].
#[derive(Debug, Clone, PartialEq)]
pub struct ProofArtifact {
    subject: ProofSubject,
    bound: Iv,
}

impl ProofArtifact {
    /// The exact proof subject.
    #[must_use]
    pub fn subject(&self) -> &ProofSubject {
        &self.subject
    }

    /// The carried interval enclosure.
    #[must_use]
    pub fn interval_bound(&self) -> Iv {
        self.bound
    }

    /// Whether this artifact is bound to this admitted problem, node, and
    /// exact domain spelling. Invalid domains or inadmissible problems return
    /// `false`; this check never upgrades a malformed subject.
    #[must_use]
    pub fn is_bound_to(&self, problem: &Problem, node: NodeId, domain: &[(f64, f64)]) -> bool {
        proof_subject(problem, node, domain).is_ok_and(|subject| subject == self.subject)
    }

    /// Whether this sealed artifact is the exact proof retained by a valid
    /// `Proven` evidence value. This compares the full subject and both
    /// interval endpoint bits; a proof for a neighboring box cannot authorize
    /// the evidence even when its human-readable bounds happen to round alike.
    #[must_use]
    pub fn verifies_evidence(&self, evidence: &ConstraintEvidence) -> bool {
        evidence.invalid_reason().is_none()
            && matches!(evidence.status, Status::Proven)
            && evidence.proof_subject.as_ref() == Some(&self.subject)
            && evidence.proof_bound.is_some_and(|bound| {
                bound.lo.to_bits() == self.bound.lo.to_bits()
                    && bound.hi.to_bits() == self.bound.hi.to_bits()
            })
    }
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
    /// Sealed exact proof-subject binding for `Proven` evidence.
    proof_subject: Option<ProofSubject>,
    /// Exact interval whose upper endpoint cleared zero for `Proven` evidence.
    proof_bound: Option<Iv>,
}

impl ConstraintEvidence {
    fn invalid_reason(&self) -> Option<&'static str> {
        if !matches!(
            self.kind,
            "hard" | "soft" | "chance" | "robust" | "certification" | "fabrication" | "code"
        ) {
            return Some("unknown-constraint-kind");
        }
        if !self.violation.is_finite() {
            return Some("nonfinite-violation");
        }
        if self.violation < 0.0 {
            return Some("negative-violation");
        }
        if !self.penalty.is_finite() {
            return Some("nonfinite-penalty");
        }
        if self.penalty < 0.0 {
            return Some("negative-penalty");
        }
        if self.kind != "soft" && self.penalty != 0.0 {
            return Some("penalty-on-nonsoft-constraint");
        }
        if matches!(self.status, Status::Proven) {
            if self.proof_subject.is_none() {
                return Some("proven-status-missing-proof-subject");
            }
            let Some(proof_bound) = self.proof_bound else {
                return Some("proven-status-missing-proof-bound");
            };
            if !proof_bound.lo.is_finite() || !proof_bound.hi.is_finite() {
                return Some("nonfinite-proof-bound");
            }
            if proof_bound.lo > proof_bound.hi {
                return Some("reversed-proof-bound");
            }
            if proof_bound.hi > 0.0 {
                return Some("proof-bound-does-not-clear-zero");
            }
        } else if self.proof_subject.is_some() || self.proof_bound.is_some() {
            return Some("unexpected-proof-evidence");
        }
        match &self.status {
            Status::Satisfied => {
                if self.violation != 0.0 || self.role != ActiveRole::Inactive {
                    return Some("inconsistent-satisfied-status");
                }
            }
            Status::Active => {
                if self.role != ActiveRole::Active {
                    return Some("inconsistent-active-status");
                }
            }
            Status::Violated => {
                if self.violation <= 0.0 || self.role != ActiveRole::Violating {
                    return Some("inconsistent-violated-status");
                }
            }
            Status::NeedsProof { .. } => {
                if self.kind != "certification" || self.role != ActiveRole::Violating {
                    return Some("inconsistent-needs-proof-status");
                }
            }
            Status::Proven => {
                if !matches!(self.kind, "robust" | "certification")
                    || self.violation != 0.0
                    || self.role != ActiveRole::Inactive
                {
                    return Some("inconsistent-proven-status");
                }
            }
            Status::BoundNotCleared { .. } => {
                if self.kind != "chance"
                    || self.violation <= 0.0
                    || self.role != ActiveRole::Violating
                {
                    return Some("inconsistent-bound-status");
                }
            }
        }
        let status_matches_kind = match self.kind {
            "chance" => matches!(
                self.status,
                Status::Satisfied | Status::Violated | Status::BoundNotCleared { .. }
            ),
            "robust" => matches!(self.status, Status::Proven | Status::Violated),
            "certification" => {
                matches!(self.status, Status::NeedsProof { .. } | Status::Proven)
            }
            "hard" | "soft" | "fabrication" | "code" => matches!(
                self.status,
                Status::Satisfied | Status::Active | Status::Violated
            ),
            _ => false,
        };
        if !status_matches_kind {
            return Some("status-kind-mismatch");
        }
        if self.certificate.kind == NumericalKind::NoClaim {
            return Some("numerical-certificate-no-claim");
        }
        let expected_numerical_kind = match (self.kind, &self.status) {
            ("chance", _) => NumericalKind::Estimate,
            ("robust", _) | ("certification", Status::Proven) => NumericalKind::Enclosure,
            _ => NumericalKind::Exact,
        };
        if self.certificate.kind != expected_numerical_kind {
            return Some("numerical-certificate-kind-mismatch");
        }
        if !self.certificate.lo.is_finite() || !self.certificate.hi.is_finite() {
            return Some("nonfinite-numerical-certificate");
        }
        if self.certificate.lo > self.certificate.hi {
            return Some("reversed-numerical-certificate");
        }
        if self.certificate.kind == NumericalKind::Exact
            && self.certificate.lo != self.certificate.hi
        {
            return Some("nonpoint-exact-certificate");
        }
        if self.violation < self.certificate.lo || self.violation > self.certificate.hi {
            return Some("violation-outside-numerical-certificate");
        }
        match self.statistical {
            StatisticalCertificate::None => {}
            StatisticalCertificate::EValue { e, alpha } => {
                if !e.is_finite() || e < 0.0 {
                    return Some("invalid-e-value");
                }
                if !(alpha.is_finite() && alpha > 0.0 && alpha < 1.0) {
                    return Some("invalid-e-value-level");
                }
            }
            StatisticalCertificate::HalfWidth {
                half_width,
                confidence,
            } => {
                if !half_width.is_finite() || half_width < 0.0 {
                    return Some("invalid-statistical-half-width");
                }
                if !(confidence.is_finite() && confidence > 0.0 && confidence < 1.0) {
                    return Some("invalid-statistical-confidence");
                }
            }
        }
        if self.kind == "chance"
            && !matches!(self.statistical, StatisticalCertificate::HalfWidth { .. })
        {
            return Some("chance-statistical-certificate-kind-mismatch");
        }
        if self.kind != "chance" && !matches!(self.statistical, StatisticalCertificate::None) {
            return Some("statistical-certificate-on-nonchance-constraint");
        }
        if let Status::BoundNotCleared {
            empirical,
            lower_bound,
        } = &self.status
        {
            if !(empirical.is_finite() && (0.0..=1.0).contains(empirical)) {
                return Some("invalid-empirical-rate");
            }
            if !lower_bound.is_finite() {
                return Some("nonfinite-lower-bound");
            }
            if lower_bound > empirical {
                return Some("lower-bound-exceeds-empirical-rate");
            }
        }
        None
    }

    /// Canonical ledger row (Rev S table shape). Dynamic text is escaped.
    /// Publicly forged malformed evidence emits a deterministic invalid,
    /// `no-claim` row rather than retaining a positive status while replacing
    /// its required numbers with JSON `null`.
    #[must_use]
    pub fn to_ledger_row(&self) -> String {
        use core::fmt::Write as _;

        if let Some(reason) = self.invalid_reason() {
            let mut row = "{\"valid\":false,\"reason\":".to_string();
            push_json_string(&mut row, reason);
            row.push_str(",\"constraint\":");
            push_json_string(&mut row, &self.name);
            row.push_str(",\"kind\":");
            push_json_string(&mut row, self.kind);
            row.push_str(",\"status\":\"no-claim\",\"violation\":null,\"penalty\":null}");
            return row;
        }

        let status = match &self.status {
            Status::Satisfied => "satisfied".to_string(),
            Status::Active => "active".to_string(),
            Status::Violated => "violated".to_string(),
            Status::NeedsProof { proof } => format!("needs-proof:{proof:?}"),
            Status::Proven => "proven".to_string(),
            Status::BoundNotCleared { .. } => "bound-not-cleared".to_string(),
        };
        let mut row = "{\"constraint\":".to_string();
        push_json_string(&mut row, &self.name);
        row.push_str(",\"kind\":");
        push_json_string(&mut row, self.kind);
        row.push_str(",\"status\":");
        push_json_string(&mut row, &status);
        row.push_str(",\"violation\":");
        let _ = write!(row, "{:.6e}", self.violation);
        row.push_str(",\"penalty\":");
        let _ = write!(row, "{:.6e}", self.penalty);
        if let Some(subject) = &self.proof_subject {
            row.push_str(",\"proof_subject\":");
            subject.push_json(&mut row);
            let proof_bound = self
                .proof_bound
                .expect("validated proven evidence has its proof bound");
            let _ = write!(
                row,
                ",\"proof_bound_bits\":[\"{:016X}\",\"{:016X}\"]",
                proof_bound.lo.to_bits(),
                proof_bound.hi.to_bits()
            );
        }
        row.push('}');
        row
    }

    /// Exact subject binding for proven evidence, when present.
    #[must_use]
    pub fn proof_subject(&self) -> Option<&ProofSubject> {
        self.proof_subject.as_ref()
    }

    /// Exact interval whose upper endpoint established the proven claim.
    #[must_use]
    pub fn proof_bound(&self) -> Option<Iv> {
        self.proof_bound
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
    /// A proof engine was asked to mint the wrong proof kind, or to prove a
    /// non-certification constraint.
    ProofKindMismatch {
        /// Proof kind declared by the constraint, if it is a certification
        /// constraint at all.
        requested: Option<ProofKind>,
        /// Proof engine the caller invoked.
        attempted: ProofKind,
    },
    /// The v1 constraint host/domain contract was not admitted.
    InvalidDomain(DomainError),
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
            ConError::ProofKindMismatch {
                requested,
                attempted,
            } => write!(
                f,
                "cannot satisfy requested proof {requested:?} with the {attempted:?} prover"
            ),
            ConError::InvalidDomain(error) => write!(f, "invalid constraint domain: {error}"),
            ConError::Parse { line, what } => write!(f, "parse error at line {line}: {what}"),
        }
    }
}

/// Allocation-free admission failures for a v1 constraint host/domain.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum DomainError {
    /// The v1 evaluator/elastic solver requires exactly one host variable.
    HostVariableCount {
        /// Number of variables declared by the host problem.
        got: usize,
    },
    /// The sole host variable is not Euclidean `Rn`.
    HostVariableManifold {
        /// Supplied manifold descriptor.
        got: Manifold,
    },
    /// The declared `Rn` dimension cannot be represented on this target.
    PointDimensionUnrepresentable {
        /// Raw dimension from the manifold descriptor.
        declared: u32,
    },
    /// The caller supplied the wrong number of point/domain components.
    DimensionMismatch {
        /// Point dimension declared by the sole `Rn` variable.
        expected: usize,
        /// Number of point/domain components supplied by the caller.
        got: usize,
    },
    /// One component range failed interval admission.
    InvalidRange {
        /// Zero-based component index.
        axis: usize,
        /// Supplied lower endpoint.
        lo: f64,
        /// Supplied upper endpoint.
        hi: f64,
        /// Exact range rule that failed.
        reason: DomainRangeError,
    },
}

/// Why one elastic-domain component range was refused.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DomainRangeError {
    /// At least one endpoint is NaN or infinite.
    NonFiniteEndpoint,
    /// The lower endpoint exceeds the upper endpoint.
    Reversed,
    /// Finite ordered endpoints have an unrepresentable difference.
    UnrepresentableSpan,
}

impl core::fmt::Display for DomainError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            DomainError::HostVariableCount { got } => {
                write!(f, "host declares {got} variables; expected exactly one")
            }
            DomainError::HostVariableManifold { got } => {
                write!(f, "host variable uses {got:?}; expected Rn")
            }
            DomainError::PointDimensionUnrepresentable { declared } => write!(
                f,
                "Rn point dimension {declared} is not representable on this target"
            ),
            DomainError::DimensionMismatch { expected, got } => write!(
                f,
                "received {got} point/domain components; the Rn variable needs {expected}"
            ),
            DomainError::InvalidRange {
                axis,
                lo,
                hi,
                reason,
            } => write!(f, "axis {axis} has {reason} (lo={lo}, hi={hi})"),
        }
    }
}

impl core::fmt::Display for DomainRangeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            DomainRangeError::NonFiniteEndpoint => f.write_str("a non-finite endpoint"),
            DomainRangeError::Reversed => f.write_str("lower endpoint above upper endpoint"),
            DomainRangeError::UnrepresentableSpan => {
                f.write_str("finite endpoints with an unrepresentable span")
            }
        }
    }
}

impl std::error::Error for DomainError {}

impl std::error::Error for ConError {}

impl From<fs_opt::OptError> for ConError {
    fn from(e: fs_opt::OptError) -> Self {
        ConError::Eval(e)
    }
}

fn require_finite_nonnegative(value: f64, what: &'static str) -> Result<(), ConError> {
    if value.is_finite() && value >= 0.0 {
        Ok(())
    } else {
        Err(ConError::BadParam { what, value })
    }
}

fn validate_spec_policy(spec: &ConstraintSpec, point_dim: Option<usize>) -> Result<(), ConError> {
    require_finite_nonnegative(spec.active_tol, "active-set tolerance")?;
    match &spec.kind {
        ConstraintKind::Soft(PenaltyLaw::Quadratic { weight })
        | ConstraintKind::Soft(PenaltyLaw::Hinge { weight }) => {
            require_finite_nonnegative(*weight, "soft-constraint weight")?;
        }
        ConstraintKind::Chance {
            level,
            estimator: ChanceEstimator::MonteCarlo { samples, delta },
        } => {
            if !(level.is_finite() && *level > 0.0 && *level < 1.0) {
                return Err(ConError::BadParam {
                    what: "chance level",
                    value: *level,
                });
            }
            if !(delta.is_finite() && *delta > 0.0 && *delta < 1.0) {
                return Err(ConError::BadParam {
                    what: "chance failure probability (delta)",
                    value: *delta,
                });
            }
            let confidence = 1.0 - *delta;
            if !(confidence > 0.0 && confidence < 1.0) {
                return Err(ConError::BadParam {
                    what: "chance confidence representation",
                    value: confidence,
                });
            }
            if *samples == 0 {
                return Err(ConError::BadParam {
                    what: "chance sample count",
                    value: f64::from(*samples),
                });
            }
        }
        ConstraintKind::Robust { half_widths } => {
            if half_widths.is_empty() {
                return Err(ConError::BadParam {
                    what: "robust half-width count",
                    value: 0.0,
                });
            }
            if let Some(point_dim) = point_dim
                && half_widths.len() != point_dim
            {
                return Err(ConError::BadParam {
                    what: "robust half-width count",
                    value: half_widths.len() as f64,
                });
            }
            for &half_width in half_widths {
                require_finite_nonnegative(half_width, "robust half-width")?;
            }
        }
        ConstraintKind::Hard
        | ConstraintKind::Certification { .. }
        | ConstraintKind::Fabrication { .. }
        | ConstraintKind::Code { .. } => {}
    }
    Ok(())
}

fn validate_interval_values(boxes: &[(f64, f64)]) -> Result<(), IvalError> {
    if boxes
        .iter()
        .any(|&(lo, hi)| !lo.is_finite() || !hi.is_finite() || lo > hi || !(hi - lo).is_finite())
    {
        Err(IvalError::BadBindings)
    } else {
        Ok(())
    }
}

fn validate_interval_bindings(problem: &Problem, boxes: &[(f64, f64)]) -> Result<(), IvalError> {
    let expected = match problem.vars() {
        [] => 0,
        [variable] => {
            let Manifold::Rn { dim } = variable.manifold else {
                return Err(IvalError::BadBindings);
            };
            usize::try_from(dim).map_err(|_| IvalError::BadBindings)?
        }
        _ => return Err(IvalError::BadBindings),
    };
    if boxes.len() != expected {
        return Err(IvalError::BadBindings);
    }
    validate_interval_values(boxes)
}

/// Evaluate a scalar node over a direct interval box after fail-closed domain
/// admission. A constant zero-variable problem takes an empty box; otherwise
/// the problem must have exactly one `Rn` host and one box per component.
/// Every endpoint must be finite and ordered, and every span must be
/// representable; malformed bindings refuse before the interval evaluator can
/// inspect the expression graph.
///
/// # Errors
/// [`IvalError::BadBindings`] for malformed boxes, otherwise the interval
/// engine's typed refusal.
pub fn interval_eval(
    problem: &Problem,
    node: NodeId,
    boxes: &[(f64, f64)],
) -> Result<Iv, IvalError> {
    validate_interval_bindings(problem, boxes)?;
    let interval = ival::interval_eval(problem, node, boxes)?;
    if !interval.lo.is_finite() || !interval.hi.is_finite() || interval.lo > interval.hi {
        return Err(IvalError::BadBindings);
    }
    Ok(interval)
}

fn proof_subject(
    problem: &Problem,
    node: NodeId,
    domain: &[(f64, f64)],
) -> Result<ProofSubject, ConError> {
    validate_interval_values(domain).map_err(|error| ConError::NotProvable {
        why: error.to_string(),
    })?;
    let expected = match problem.vars() {
        [] => 0,
        [variable] => {
            let Manifold::Rn { dim } = variable.manifold else {
                return Err(ConError::NotProvable {
                    why: format!(
                        "interval proof subjects require an Rn variable, got {:?}",
                        variable.manifold
                    ),
                });
            };
            usize::try_from(dim).map_err(|_| ConError::NotProvable {
                why: format!("interval proof dimension {dim} is not representable on this target"),
            })?
        }
        variables => {
            return Err(ConError::NotProvable {
                why: format!(
                    "interval proof subjects require a zero-variable constant problem or one Rn \
                     variable, got {} variables",
                    variables.len()
                ),
            });
        }
    };
    if domain.len() != expected {
        return Err(ConError::NotProvable {
            why: format!(
                "interval proof domain has {} axes but the Rn variable has {expected}",
                domain.len()
            ),
        });
    }
    let admission = problem.admit().map_err(|report| ConError::NotProvable {
        why: format!("proof subject problem did not admit: {report}"),
    })?;
    Ok(ProofSubject {
        problem: admission.semantic_id(),
        node,
        domain_bits: domain
            .iter()
            .map(|&(lo, hi)| (lo.to_bits(), hi.to_bits()))
            .collect(),
    })
}

fn robust_boxes(x: &[f64], half_widths: &[f64]) -> Result<Vec<(f64, f64)>, ConError> {
    if half_widths.len() != x.len() {
        return Err(ConError::BadParam {
            what: "robust half-width count",
            value: half_widths.len() as f64,
        });
    }
    for (&center, &half_width) in x.iter().zip(half_widths) {
        require_finite_nonnegative(half_width, "robust half-width")?;
        if !center.is_finite() {
            return Err(ConError::BadParam {
                what: "robust interval center",
                value: center,
            });
        }
        let lo = center - half_width;
        let hi = center + half_width;
        if !lo.is_finite() || !hi.is_finite() || lo > hi || !(hi - lo).is_finite() {
            return Err(ConError::BadParam {
                what: "robust interval range",
                value: half_width,
            });
        }
    }
    Ok(x.iter()
        .zip(half_widths)
        .map(|(&center, &half_width)| (center - half_width, center + half_width))
        .collect())
}

pub(crate) fn scalar_at(problem: &Problem, node: NodeId, x: &[f64]) -> Result<f64, ConError> {
    let v = fs_opt::eval(problem, node, std::slice::from_ref(&x.to_vec()))?;
    v.scalar().ok_or(ConError::NotScalar { node: node.0 })
}

fn validate_evaluate_host(problem: &Problem, x: &[f64]) -> Result<usize, ConError> {
    if problem.vars().len() != 1 {
        return Err(ConError::InvalidDomain(DomainError::HostVariableCount {
            got: problem.vars().len(),
        }));
    }
    let variable = &problem.vars()[0];
    let Manifold::Rn { dim } = variable.manifold else {
        return Err(ConError::InvalidDomain(DomainError::HostVariableManifold {
            got: variable.manifold,
        }));
    };
    let expected = usize::try_from(dim).map_err(|_| {
        ConError::InvalidDomain(DomainError::PointDimensionUnrepresentable { declared: dim })
    })?;
    if x.len() != expected {
        return Err(ConError::InvalidDomain(DomainError::DimensionMismatch {
            expected,
            got: x.len(),
        }));
    }
    Ok(expected)
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
    // The v1 evaluator's uncertainty model is raw Euclidean addition. Admit
    // that exact host shape before graph evaluation or a chance-noise callback
    // can observe work; manifold-aware perturbations belong to a successor.
    let point_dim = validate_evaluate_host(problem, x)?;
    // Policy values are public fields and therefore untrusted. Admit them
    // before evaluating the expression so NaN tolerances and negative weights
    // cannot turn a real violation into a satisfied or rewarded result.
    validate_spec_policy(spec, Some(point_dim))?;
    let g = scalar_at(problem, spec.node, x)?;
    // fs-opt currently returns a typed EvalNonFinite refusal before this point.
    // Retain the finite guard as defense in depth for any future evaluator that
    // carries a non-finite scalar value: IEEE comparisons and `max` must never
    // turn an undefined constraint into an exact-zero satisfied claim.
    let finite = g.is_finite();
    let violation = if finite { g.max(0.0) } else { f64::INFINITY };
    let base_status = if !finite || g > spec.active_tol {
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
        proof_subject: None,
        proof_bound: None,
    };
    match &spec.kind {
        ConstraintKind::Hard | ConstraintKind::Fabrication { .. } | ConstraintKind::Code { .. } => {
        }
        ConstraintKind::Soft(law) => {
            ev.penalty = match *law {
                PenaltyLaw::Quadratic { weight } => weight * violation * violation,
                PenaltyLaw::Hinge { weight } => weight * violation,
            };
            if !ev.penalty.is_finite() || ev.penalty < 0.0 {
                return Err(ConError::BadParam {
                    what: "soft-constraint penalty result",
                    value: ev.penalty,
                });
            }
        }
        ConstraintKind::Chance { level, estimator } => {
            let ChanceEstimator::MonteCarlo { samples, delta } = *estimator;
            let noise = noise.ok_or(ConError::BadParam {
                what: "chance noise model (required)",
                value: f64::NAN,
            })?;
            let mut hits = 0u32;
            for s in 0..samples {
                let draw = noise(u64::from(s));
                if draw.len() != x.len() {
                    return Err(ConError::BadParam {
                        what: "chance noise draw dimension",
                        value: draw.len() as f64,
                    });
                }
                let shifted: Vec<f64> = x.iter().zip(&draw).map(|(a, b)| a + b).collect();
                if scalar_at(problem, spec.node, &shifted)? <= 0.0 {
                    hits += 1;
                }
            }
            let empirical = f64::from(hits) / f64::from(samples);
            // Hoeffding lower confidence bound at failure prob delta.
            // `-ln(delta)` is algebraically `ln(1/delta)` but remains finite
            // for every positive finite binary64 delta; forming `1/delta`
            // first would overflow for valid subnormal policy values.
            let half_width = (-delta.ln() / (2.0 * f64::from(samples))).sqrt();
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
            let boxes = robust_boxes(x, half_widths)?;
            match interval_eval(problem, spec.node, &boxes) {
                Ok(iv) => {
                    ev.violation = iv.hi.max(0.0);
                    ev.certificate =
                        NumericalCertificate::enclosure(iv.lo.max(0.0), iv.hi.max(0.0));
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
                    if matches!(ev.status, Status::Proven) {
                        ev.proof_subject = Some(proof_subject(problem, spec.node, &boxes)?);
                        ev.proof_bound = Some(iv);
                    }
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
    let requested = match &spec.kind {
        ConstraintKind::Certification { proof } => Some(*proof),
        _ => None,
    };
    if requested != Some(ProofKind::Interval) {
        return Err(ConError::ProofKindMismatch {
            requested,
            attempted: ProofKind::Interval,
        });
    }
    validate_spec_policy(spec, None)?;
    let subject = proof_subject(problem, spec.node, domain)?;
    let iv = interval_eval(problem, spec.node, domain)
        .map_err(|e| ConError::NotProvable { why: e.to_string() })?;
    if iv.hi <= 0.0 {
        Ok((
            ConstraintEvidence {
                name: spec.name.clone(),
                kind: spec.kind.kind_name(),
                status: Status::Proven,
                violation: 0.0,
                certificate: NumericalCertificate::enclosure(0.0, 0.0),
                statistical: StatisticalCertificate::None,
                role: ActiveRole::Inactive,
                penalty: 0.0,
                proof_subject: Some(subject.clone()),
                proof_bound: Some(iv),
            },
            ProofArtifact { subject, bound: iv },
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

/// Encode a constraint set in canonical line form (floats as bits).
///
/// This infallible writer preserves even untrusted public field bits. Use
/// [`parse_specs`] as the admission boundary before treating those bytes as a
/// valid policy; writer/parser identity is guaranteed for admitted specs.
#[must_use]
pub fn serialize_specs(specs: &[ConstraintSpec]) -> String {
    use std::fmt::Write as _;
    let hex = |v: f64| format!("{:016X}", v.to_bits());
    let mut s = String::from("fscon v2\n");
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
                format!("fabrication {}", encode_wire_token(process))
            }
            ConstraintKind::Code { standard } => {
                format!("code {}", encode_wire_token(standard))
            }
        };
        let _ = writeln!(
            s,
            "constraint {} {} {} {kind}",
            encode_wire_token(&c.name),
            c.node.0,
            hex(c.active_tol)
        );
    }
    s
}

/// Parse and admit [`serialize_specs`] output (admitted round-trip identity).
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
    let mut saw_header = false;
    for (ln0, line) in text.lines().enumerate() {
        let ln = ln0 + 1;
        let toks: Vec<&str> = line.split(' ').collect();
        match toks.first().copied() {
            Some("fscon") => {
                if ln != 1 || saw_header || toks.as_slice() != ["fscon", "v2"] {
                    return Err(perr(ln, "expected the single canonical `fscon v2` header"));
                }
                saw_header = true;
            }
            Some("constraint") => {
                if !saw_header {
                    return Err(perr(ln, "constraint appears before the fscon header"));
                }
                let name = decode_wire_token(
                    toks.get(1).ok_or_else(|| perr(ln, "missing name"))?,
                    ln,
                    "constraint name",
                )?;
                let node: u32 = toks
                    .get(2)
                    .and_then(|t| t.parse().ok())
                    .ok_or_else(|| perr(ln, "bad node"))?;
                let active_tol = toks
                    .get(3)
                    .and_then(|t| unhex(t))
                    .ok_or_else(|| perr(ln, "bad tol"))?;
                let kind = match toks.get(4).copied() {
                    Some("hard") => {
                        if toks.len() != 5 {
                            return Err(perr(ln, "hard constraint has trailing fields"));
                        }
                        ConstraintKind::Hard
                    }
                    Some("soft") => {
                        if toks.len() != 7 {
                            return Err(perr(ln, "soft constraint has the wrong field count"));
                        }
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
                        if toks.len() != 9 {
                            return Err(perr(ln, "chance constraint has the wrong field count"));
                        }
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
                        if toks.len() != 6 {
                            return Err(perr(ln, "robust constraint has the wrong field count"));
                        }
                        let ws = toks.get(5).ok_or_else(|| perr(ln, "missing widths"))?;
                        let half_widths: Option<Vec<f64>> = ws.split(',').map(unhex).collect();
                        ConstraintKind::Robust {
                            half_widths: half_widths.ok_or_else(|| perr(ln, "bad widths"))?,
                        }
                    }
                    Some("certification") => {
                        if toks.len() != 6 {
                            return Err(perr(
                                ln,
                                "certification constraint has the wrong field count",
                            ));
                        }
                        match toks.get(5).copied() {
                            Some("interval") => ConstraintKind::Certification {
                                proof: ProofKind::Interval,
                            },
                            Some("sos") => ConstraintKind::Certification {
                                proof: ProofKind::Sos,
                            },
                            _ => return Err(perr(ln, "unknown proof kind")),
                        }
                    }
                    Some("fabrication") => {
                        if toks.len() != 6 {
                            return Err(perr(
                                ln,
                                "fabrication constraint has the wrong field count",
                            ));
                        }
                        ConstraintKind::Fabrication {
                            process: decode_wire_token(
                                toks.get(5).ok_or_else(|| perr(ln, "missing process"))?,
                                ln,
                                "fabrication process",
                            )?,
                        }
                    }
                    Some("code") => {
                        if toks.len() != 6 {
                            return Err(perr(ln, "code constraint has the wrong field count"));
                        }
                        ConstraintKind::Code {
                            standard: decode_wire_token(
                                toks.get(5).ok_or_else(|| perr(ln, "missing standard"))?,
                                ln,
                                "code standard",
                            )?,
                        }
                    }
                    _ => return Err(perr(ln, "unknown kind")),
                };
                let spec = ConstraintSpec {
                    name,
                    node: NodeId(node),
                    kind,
                    active_tol,
                };
                validate_spec_policy(&spec, None)
                    .map_err(|error| perr(ln, &format!("invalid policy: {error}")))?;
                out.push(spec);
            }
            Some("") | None => return Err(perr(ln, "blank lines are not canonical")),
            Some(other) => {
                return Err(ConError::Parse {
                    line: ln,
                    what: format!("unknown directive `{other}`"),
                });
            }
        }
    }
    if !saw_header {
        return Err(perr(1, "missing `fscon v2` header"));
    }
    let canonical = serialize_specs(&out);
    if canonical != text {
        let mismatch = canonical
            .bytes()
            .zip(text.bytes())
            .take_while(|(expected, actual)| expected == actual)
            .count();
        let line = text.as_bytes()[..mismatch.min(text.len())]
            .iter()
            .filter(|&&byte| byte == b'\n')
            .count()
            + 1;
        return Err(perr(
            line,
            "input is not the canonical fscon v2 byte spelling",
        ));
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

    #[test]
    fn internally_corrupted_proof_bounds_lose_claim_authority() {
        let mut builder = fs_opt::ProblemBuilder::new();
        let variable = builder
            .var("x", Manifold::Rn { dim: 1 }, fs_qty::Dims::NONE)
            .expect("variable");
        let variable_ref = builder.var_ref(variable).expect("variable ref");
        let x = builder.component(variable_ref, 0).expect("component");
        let one = builder.konst(1.0, fs_qty::Dims::NONE).expect("constant");
        let constraint = builder.sub(x, one).expect("constraint");
        let objective = builder.norm_sq(variable_ref).expect("objective");
        builder
            .objective(objective, fs_opt::Sense::Minimize, 1.0)
            .expect("objective entry");
        let problem = builder.finish();
        let spec = ConstraintSpec {
            name: "proof-bound-corruption".into(),
            node: constraint,
            kind: ConstraintKind::Certification {
                proof: ProofKind::Interval,
            },
            active_tol: 0.0,
        };
        let domain = [(0.0, 0.5)];
        let (mut evidence, artifact) =
            prove_interval(&problem, &spec, &domain).expect("valid interval proof");
        assert!(artifact.verifies_evidence(&evidence));

        evidence.proof_bound = Some(Iv { lo: -1.0, hi: 1.0 });
        assert!(
            evidence
                .to_ledger_row()
                .contains("\"reason\":\"proof-bound-does-not-clear-zero\"")
        );
        assert!(!artifact.verifies_evidence(&evidence));

        evidence.proof_bound = Some(Iv {
            lo: f64::NAN,
            hi: -0.5,
        });
        assert!(
            evidence
                .to_ledger_row()
                .contains("\"reason\":\"nonfinite-proof-bound\"")
        );
        assert!(!artifact.verifies_evidence(&evidence));
    }
}
