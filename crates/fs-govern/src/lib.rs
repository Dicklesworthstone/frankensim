//! fs-govern — the addendum's machine-readable RISK REGISTER (plan addendum,
//! Part V). Layer: UTIL (pure data + audit; `fs-blake3` identity dependency).
//!
//! Design principle P8 says the plan itself must be falsifiable, and
//! Governance Rule 2 says a risk (or kill criterion) whose measurement was
//! never instrumented counts as unmanaged — "unmeasured survival is not
//! survival". This crate makes that enforceable: it encodes the ten named
//! risks R1–R10, each with its mitigation, its EARLY-WARNING METRIC, the
//! THRESHOLD at which to act, and the OWNING bead — and provides an [`audit`]
//! that a CI gate can run to fail the build if any risk is missing a metric or
//! an owner, and to report how many early-warning metrics are actually
//! instrumented.
//!
//! The register is the canonical record; [`to_json`] emits it for dashboards.
//! Instrumentation is EVIDENCE, not a flag (bead xpck.9): a risk counts as
//! operationally managed only when it carries a fresh, content-identified
//! [`InstrumentationReceipt`]; every receipt is `None` today
//! (the honest baseline), so the register is schema-complete but
//! operationally RED — and the audits say both things separately.
//!
//! Sibling modules encode the rest of the addendum's governance as data: the
//! design principles + governance rules ([`doctrine`]) and the nineteen
//! proposals with their kill metrics + owning beads + a completeness audit
//! ([`proposals`]).

pub mod crates;
pub mod doctrine;
pub mod proposals;

pub use crates::{AddendumCrate, CrateAudit, addendum_crates, crate_audit, crates_json};
pub use doctrine::{GovernanceRule, PRINCIPLES, Principle, RULES, principles, rules};
pub use fs_blake3::ContentHash;
pub use proposals::{GovernanceAudit, Proposal, governance_audit, proposals, proposals_json};

/// The ten addendum risks (Part V).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RiskId {
    /// Estimator constants (Proposal 9's hard core).
    R1,
    /// Determinism at scale.
    R2,
    /// Stable entity identity.
    R3,
    /// Loose sensitivity bounds.
    R4,
    /// Spectral-gap fragility.
    R5,
    /// Restriction-map quality.
    R6,
    /// Model-form composition.
    R7,
    /// Registration well-posedness.
    R8,
    /// Standards-body latency.
    R9,
    /// Breadth death.
    R10,
}

impl RiskId {
    /// Every risk id, in order.
    pub const ALL: [RiskId; 10] = [
        RiskId::R1,
        RiskId::R2,
        RiskId::R3,
        RiskId::R4,
        RiskId::R5,
        RiskId::R6,
        RiskId::R7,
        RiskId::R8,
        RiskId::R9,
        RiskId::R10,
    ];

    /// The stable code (`"R1"` … `"R10"`).
    #[must_use]
    pub fn code(self) -> &'static str {
        match self {
            RiskId::R1 => "R1",
            RiskId::R2 => "R2",
            RiskId::R3 => "R3",
            RiskId::R4 => "R4",
            RiskId::R5 => "R5",
            RiskId::R6 => "R6",
            RiskId::R7 => "R7",
            RiskId::R8 => "R8",
            RiskId::R9 => "R9",
            RiskId::R10 => "R10",
        }
    }
}

/// A content-identified instrumentation receipt (bead xpck.9): the claim
/// "this metric is live on a dashboard" carried with its supporting evidence,
/// not as a mutable boolean.
///
/// The fields are private so a caller cannot accidentally mutate a receipt
/// without invalidating its identity. The identity binds the subject,
/// dashboard, verifier, evidence artifact, and verification day. It is an
/// unkeyed content identity, not a signature: issuer authorization and the
/// scientific adequacy of the referenced evidence remain external policy.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct InstrumentationReceipt {
    /// Dashboard locator (non-empty).
    dashboard: &'static str,
    /// Day the feed was last verified live (days since 2026-01-01).
    verified_day: u32,
    /// Identity of the person, service, or policy that performed the check.
    verifier: &'static str,
    /// Content address of the evidence supporting the live-feed check.
    evidence_artifact: ContentHash,
    /// Domain-separated identity of every preceding field plus the subject.
    identity: ContentHash,
}

/// Receipts older than this demote to [`InstrumentationStatus::Stale`].
pub const MAX_RECEIPT_AGE_DAYS: u32 = 45;

/// Domain for canonical instrumentation-receipt identities.
pub const INSTRUMENTATION_RECEIPT_IDENTITY_DOMAIN: &str =
    "frankensim.fs-govern.instrumentation-receipt.v1";

/// Why an instrumentation receipt could not be constructed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReceiptError {
    /// The governed risk or proposal id was empty.
    EmptySubject,
    /// The dashboard locator was empty.
    EmptyDashboard,
    /// The verifier identity was empty.
    EmptyVerifier,
    /// The supporting evidence artifact used the all-zero missing-value
    /// sentinel rather than a content address.
    EmptyEvidenceArtifact,
}

impl core::fmt::Display for ReceiptError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str(match self {
            ReceiptError::EmptySubject => "instrumentation receipt subject is empty",
            ReceiptError::EmptyDashboard => "instrumentation receipt dashboard is empty",
            ReceiptError::EmptyVerifier => "instrumentation receipt verifier is empty",
            ReceiptError::EmptyEvidenceArtifact => {
                "instrumentation receipt evidence artifact is the all-zero missing-value sentinel"
            }
        })
    }
}

fn push_identity_field(out: &mut Vec<u8>, tag: u8, bytes: &[u8]) {
    out.push(tag);
    let len = u64::try_from(bytes.len()).expect("field length fits in u64");
    out.extend_from_slice(&len.to_le_bytes());
    out.extend_from_slice(bytes);
}

/// Derive the canonical content identity for an instrumentation receipt.
///
/// Every variable-length field is tagged and length-prefixed, and BLAKE3's
/// derive-key mode separates this identity from all other artifact types.
#[must_use]
pub fn receipt_identity(
    subject: &str,
    dashboard: &str,
    verifier: &str,
    evidence_artifact: ContentHash,
    verified_day: u32,
) -> ContentHash {
    let mut canonical = Vec::new();
    push_identity_field(&mut canonical, 1, subject.as_bytes());
    push_identity_field(&mut canonical, 2, dashboard.as_bytes());
    push_identity_field(&mut canonical, 3, verifier.as_bytes());
    push_identity_field(&mut canonical, 4, evidence_artifact.as_bytes());
    push_identity_field(&mut canonical, 5, &verified_day.to_le_bytes());
    fs_blake3::hash_domain(INSTRUMENTATION_RECEIPT_IDENTITY_DOMAIN, &canonical)
}

impl InstrumentationReceipt {
    /// Construct a sealed receipt whose identity binds all semantic fields.
    pub fn new(
        subject: &str,
        dashboard: &'static str,
        verifier: &'static str,
        evidence_artifact: ContentHash,
        verified_day: u32,
    ) -> Result<Self, ReceiptError> {
        if subject.trim().is_empty() {
            return Err(ReceiptError::EmptySubject);
        }
        if dashboard.trim().is_empty() {
            return Err(ReceiptError::EmptyDashboard);
        }
        if verifier.trim().is_empty() {
            return Err(ReceiptError::EmptyVerifier);
        }
        if evidence_artifact.as_bytes().iter().all(|byte| *byte == 0) {
            return Err(ReceiptError::EmptyEvidenceArtifact);
        }
        Ok(Self {
            dashboard,
            verified_day,
            verifier,
            evidence_artifact,
            identity: receipt_identity(
                subject,
                dashboard,
                verifier,
                evidence_artifact,
                verified_day,
            ),
        })
    }

    /// Dashboard locator asserted by this receipt.
    #[must_use]
    pub fn dashboard(self) -> &'static str {
        self.dashboard
    }

    /// Day on which the dashboard feed was checked.
    #[must_use]
    pub fn verified_day(self) -> u32 {
        self.verified_day
    }

    /// Identity of the verifier that asserted the live-feed check.
    #[must_use]
    pub fn verifier(self) -> &'static str {
        self.verifier
    }

    /// Content address of the supporting evidence artifact.
    #[must_use]
    pub fn evidence_artifact(self) -> ContentHash {
        self.evidence_artifact
    }

    /// Canonical content identity of this receipt.
    #[must_use]
    pub fn identity(self) -> ContentHash {
        self.identity
    }

    /// Whether this receipt's sealed identity is valid for `subject`.
    #[must_use]
    pub fn is_consistent_for(self, subject: &str) -> bool {
        self.identity
            == receipt_identity(
                subject,
                self.dashboard,
                self.verifier,
                self.evidence_artifact,
                self.verified_day,
            )
    }

    /// Deterministic receipt provenance for dashboards and ledger records.
    #[must_use]
    pub fn to_json(self) -> String {
        use core::fmt::Write as _;
        let mut out = String::new();
        write!(
            out,
            "{{\"dashboard\":\"{}\",\"verified_day\":{},\"verifier\":\"{}\",\"evidence_artifact\":\"{}\",\"identity\":\"{}\"}}",
            json_escape(self.dashboard),
            self.verified_day,
            json_escape(self.verifier),
            self.evidence_artifact,
            self.identity,
        )
        .expect("writing to a String is infallible");
        out
    }
}

/// Operational status of one subject's instrumentation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InstrumentationStatus {
    /// Receipt present, content-identity-consistent, and fresh.
    Verified {
        /// Days since the feed was verified live.
        age_days: u32,
    },
    /// Receipt content-identity-consistent but older than
    /// [`MAX_RECEIPT_AGE_DAYS`] — the dashboard may be dead.
    Stale {
        /// Days since the feed was verified live.
        age_days: u32,
    },
    /// Receipt present but its content identity is inconsistent with the
    /// governed subject, or its verification date is in the future.
    BadReceipt,
    /// No receipt at all.
    Uninstrumented,
}

impl InstrumentationStatus {
    /// Stable lowercase name for JSON/ledger rows.
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            InstrumentationStatus::Verified { .. } => "verified",
            InstrumentationStatus::Stale { .. } => "stale",
            InstrumentationStatus::BadReceipt => "bad-receipt",
            InstrumentationStatus::Uninstrumented => "uninstrumented",
        }
    }
}

/// Judge one subject's receipt as of `today_day` (days since
/// 2026-01-01). Fail-closed: only a content-identity-consistent, non-empty,
/// fresh receipt is [`InstrumentationStatus::Verified`].
#[must_use]
pub fn instrumentation_status(
    subject: &str,
    receipt: Option<&InstrumentationReceipt>,
    today_day: u32,
) -> InstrumentationStatus {
    let Some(r) = receipt else {
        return InstrumentationStatus::Uninstrumented;
    };
    if r.dashboard.trim().is_empty()
        || r.verifier.trim().is_empty()
        || !r.is_consistent_for(subject)
        || r.verified_day > today_day
    {
        return InstrumentationStatus::BadReceipt;
    }
    let age_days = today_day - r.verified_day;
    if age_days > MAX_RECEIPT_AGE_DAYS {
        InstrumentationStatus::Stale { age_days }
    } else {
        InstrumentationStatus::Verified { age_days }
    }
}

/// One risk-register entry.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Risk {
    /// The risk id.
    pub id: RiskId,
    /// Short name.
    pub name: &'static str,
    /// What can go wrong.
    pub description: &'static str,
    /// The mitigation already embedded in the proposals.
    pub mitigation: &'static str,
    /// The EARLY-WARNING METRIC that makes the risk visible before it is fatal.
    pub early_warning: &'static str,
    /// The threshold / condition at which to act.
    pub threshold: &'static str,
    /// The owning bead (or governance owner).
    pub owner: &'static str,
    /// Evidence that the early-warning metric is live on a dashboard
    /// (`None` = uninstrumented; a bare flag cannot claim coverage).
    pub receipt: Option<InstrumentationReceipt>,
}

/// The canonical R1–R10 register.
const REGISTER: [Risk; 10] = [
    Risk {
        id: RiskId::R1,
        name: "Estimator constants",
        description: "Guaranteed a-posteriori bounds beyond the elliptic/FEEC class need genuinely hard reliability constants.",
        mitigation: "Constant-free equilibrated-flux estimators for the elliptic class first; DWR accepts carry estimated color; nonlinear/transient candidates are warm-starts only.",
        early_warning: "accept-rate stratified by problem class",
        threshold: "a problem class stuck at warm-start-only (accept-rate near zero)",
        owner: "frankensim-epic-flywheel-lmp4.1",
        receipt: None,
    },
    Risk {
        id: RiskId::R2,
        name: "Determinism at scale",
        description: "GPU reductions and dynamic scheduling threaten the bitwise reproducibility Proposals 2 and 10 stand on.",
        mitigation: "Fixed-order reductions and deterministic scheduling as certified contracts with Gauntlet tests; the measured perf tax is accepted explicitly.",
        early_warning: "reproducibility test failures per release",
        threshold: "determinism tax exceeding ~15% on dominant kernels",
        owner: "frankensim-epic-flywheel-lmp4.6",
        receipt: None,
    },
    Risk {
        id: RiskId::R3,
        name: "Stable entity identity",
        description: "Semantic diff and merge require persistent IDs through topology-changing edits (the CAD topological-naming problem).",
        mitigation: "FrankenSim owns its kernel; IDs are first-class from day one; edits are ledgered ops that transform IDs explicitly.",
        early_warning: "fraction of diffs falling back to unattributed geometric comparison",
        threshold: "a material fraction of diffs cannot attribute to a causal edit",
        owner: "frankensim-epic-flywheel-lmp4.10",
        receipt: None,
    },
    Risk {
        id: RiskId::R4,
        name: "Loose sensitivity bounds",
        description: "Interval derivatives on nonlinear ops may be so pessimistic that the recompute frontier balloons to the whole DAG.",
        mitigation: "Graceful degradation to hash memoization; adjoint-sharpened bounds where Proposal 1 is live; per-op skip-yield telemetry.",
        early_warning: "skip-yield (fraction of DAG certifiably skipped)",
        threshold: "skip-yield delivers <2x median wall-clock vs hash memoization",
        owner: "frankensim-epic-flywheel-lmp4.7",
        receipt: None,
    },
    Risk {
        id: RiskId::R5,
        name: "Spectral-gap fragility",
        description: "The fixable/structural triage (and merge adjudication) degrades in ill-conditioned regions.",
        mitigation: "Spectral-health monitoring with mandatory low-confidence propagation into merge outputs; the Gauntlet refusal-is-the-pass suite.",
        early_warning: "gap-collapse incidence per assembly class",
        threshold: "gap collapse observed outside synthetic cases at volume",
        owner: "frankensim-epic-selfknow-knh1.3",
        receipt: None,
    },
    Risk {
        id: RiskId::R6,
        name: "Restriction-map quality",
        description: "The sheaf propagates garbage with certificates attached if trace/conversion operators are inaccurate.",
        mitigation: "The Proposal 7 conformance suite (functoriality, adjoint consistency, MMS tolerance honesty) applied to first-party converters with third-party severity.",
        early_warning: "conformance-tier distribution of converters on the hot path",
        threshold: "hot-path converters cluster in low conformance tiers",
        owner: "frankensim-epic-gtm-jwq8.2",
        receipt: None,
    },
    Risk {
        id: RiskId::R7,
        name: "Model-form composition",
        description: "Estimated-color quantities do not compose with the clean algebra of verified bounds.",
        mitigation: "The type system's laundering refusal; discrepancy-probe maps as the empirical substitute; the weakest-input rule on headlines.",
        early_warning: "audit rate of estimated-color claims in decision-critical positions",
        threshold: "estimated-color claims silently drive decisions without probe maps",
        owner: "frankensim-epic-epistype-qmao.2",
        receipt: None,
    },
    Risk {
        id: RiskId::R8,
        name: "Registration well-posedness",
        description: "Scan-to-design alignment error can exceed the deviations being certified.",
        mitigation: "Design-for-verification fiducials pushed upstream; point-sensor assimilation as the registration-free fallback.",
        early_warning: "registration-uncertainty-to-signal ratio per part class",
        threshold: "registration uncertainty exceeds the geometric deviations being certified",
        owner: "frankensim-epic-coupling-bk0o.4",
        receipt: None,
    },
    Risk {
        id: RiskId::R9,
        name: "Standards-body latency",
        description: "Machine-checkable evidence may sit ahead of what auditors will engage with.",
        mitigation: "The vocabulary crosswalk speaks their language; the package doubles as internal QA and B2B diligence collateral, so investment is not stranded.",
        early_warning: "auditor engagement rate in the first regulated-vertical cycle",
        threshold: "no auditor engages the machine-checkable format even as supplementary evidence",
        owner: "frankensim-epic-epistype-qmao.9",
        receipt: None,
    },
    Risk {
        id: RiskId::R10,
        name: "Breadth death",
        description: "Nineteen proposals executed as nineteen parallel programs is an obituary.",
        mitigation: "The governance rules: one research bet at a time, quarterly kill enforcement, phase gating.",
        early_warning: "headcount-weighted work-in-progress outside the current phase",
        threshold: "WIP outside the current phase grows materially",
        owner: "frankensim-epic-addendum-xpck.1",
        receipt: None,
    },
];

/// The canonical risk register (R1–R10, in order).
#[must_use]
pub fn register() -> &'static [Risk] {
    &REGISTER
}

/// Look up a single risk.
#[must_use]
pub fn risk(id: RiskId) -> &'static Risk {
    // ALL and REGISTER share the same order, so the index is the position.
    let idx = RiskId::ALL
        .iter()
        .position(|r| *r == id)
        .expect("every RiskId is registered");
    &REGISTER[idx]
}

/// The result of auditing the register: DECLARATION (schema) and LIVE
/// OPERATION (receipts) are distinct verdicts — collapsing them was
/// the false-green this bead removed (xpck.9).
#[derive(Debug, Clone, PartialEq)]
pub struct RiskAudit {
    /// Total risks.
    pub total: usize,
    /// Risks that DECLARE both a non-empty early-warning metric and an owner.
    pub declared: usize,
    /// Risks whose early-warning metric carries a fresh, identity-consistent
    /// receipt with verifier and evidence provenance.
    pub verified_instrumented: usize,
    /// `(risk, reason)` for every declaration gap.
    pub schema_gaps: Vec<(RiskId, &'static str)>,
    /// `(risk, status)` for every risk NOT verified live — the exact
    /// operational gaps Governance Rule 2 demands be visible.
    pub operational_gaps: Vec<(RiskId, InstrumentationStatus)>,
}

impl RiskAudit {
    /// Does every risk DECLARE a metric and an owner? (Schema only —
    /// says nothing about whether anything is actually measured.)
    #[must_use]
    pub fn declared_schema_ok(&self) -> bool {
        self.schema_gaps.is_empty()
    }

    /// Is every risk OPERATIONALLY managed — declared AND its metric carries
    /// a fresh, identity-consistent receipt? Fails closed on
    /// any uninstrumented, stale, or bad-receipt entry. "Unmeasured
    /// survival is not survival."
    #[must_use]
    pub fn operationally_managed(&self) -> bool {
        self.declared_schema_ok() && self.operational_gaps.is_empty()
    }
}

/// Audit the canonical register as of `today_day` (see [`audit_slice`]).
#[must_use]
pub fn audit(today_day: u32) -> RiskAudit {
    audit_slice(&REGISTER, today_day)
}

/// Audit an arbitrary risk slice: every risk must carry an early-warning
/// metric and an owner (declaration), and Governance Rule 2 demands the
/// metric be VERIFIED LIVE (a risk with no working measurement is
/// unmanaged) — the two verdicts are reported separately and the
/// operational one fails closed.
#[must_use]
pub fn audit_slice(risks: &[Risk], today_day: u32) -> RiskAudit {
    let mut schema_gaps = Vec::new();
    let mut operational_gaps = Vec::new();
    let mut declared = 0usize;
    let mut verified_instrumented = 0usize;
    for r in risks {
        let mut ok = true;
        if r.early_warning.trim().is_empty() {
            schema_gaps.push((r.id, "missing early-warning metric"));
            ok = false;
        }
        if r.owner.trim().is_empty() {
            schema_gaps.push((r.id, "missing owner"));
            ok = false;
        }
        if ok {
            declared += 1;
        }
        match instrumentation_status(r.id.code(), r.receipt.as_ref(), today_day) {
            InstrumentationStatus::Verified { .. } => verified_instrumented += 1,
            other => operational_gaps.push((r.id, other)),
        }
    }
    RiskAudit {
        total: risks.len(),
        declared,
        verified_instrumented,
        schema_gaps,
        operational_gaps,
    }
}

/// Escape a string for embedding in JSON.
pub(crate) fn json_escape(s: &str) -> String {
    use core::fmt::Write as _;

    let mut out = String::with_capacity(s.len() + 2);
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\t' => out.push_str("\\t"),
            '\r' => out.push_str("\\r"),
            '\u{0000}'..='\u{001f}' => {
                write!(out, "\\u{:04x}", u32::from(c)).expect("writing to a String is infallible");
            }
            _ => out.push(c),
        }
    }
    out
}

/// Emit the register as a machine-readable JSON array (one object per
/// risk) as of `today_day`, for dashboards and CI gates. Deterministic
/// (risks in order). Each entry carries its instrumentation STATUS
/// (verified/stale/bad-receipt/uninstrumented) — never an ambiguous
/// "complete" flag.
#[must_use]
pub fn to_json(today_day: u32) -> String {
    use core::fmt::Write as _;
    let mut out = String::from("[");
    for (i, r) in REGISTER.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        let status = instrumentation_status(r.id.code(), r.receipt.as_ref(), today_day);
        write!(
            out,
            "{{\"id\":\"{}\",\"name\":\"{}\",\"early_warning\":\"{}\",\"threshold\":\"{}\",\"owner\":\"{}\",\"instrumentation\":\"{}\",\"receipt\":{},\"mitigation\":\"{}\"}}",
            r.id.code(),
            json_escape(r.name),
            json_escape(r.early_warning),
            json_escape(r.threshold),
            json_escape(r.owner),
            status.name(),
            r.receipt
                .map_or_else(|| "null".to_owned(), InstrumentationReceipt::to_json),
            json_escape(r.mitigation),
        )
        .expect("writing to a String is infallible");
    }
    out.push(']');
    out
}
