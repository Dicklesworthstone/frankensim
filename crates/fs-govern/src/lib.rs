//! fs-govern — the addendum's machine-readable RISK REGISTER (plan addendum,
//! Part V). Layer: UTIL (pure data + audit; no dependencies).
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
//! `instrumented` defaults to `false` for every risk (nothing is wired yet —
//! the honest baseline); a risk flips to instrumented only when its
//! early-warning metric is live on a dashboard.
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
    /// Is the early-warning metric actually live on a dashboard yet?
    pub instrumented: bool,
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
        instrumented: false,
    },
    Risk {
        id: RiskId::R2,
        name: "Determinism at scale",
        description: "GPU reductions and dynamic scheduling threaten the bitwise reproducibility Proposals 2 and 10 stand on.",
        mitigation: "Fixed-order reductions and deterministic scheduling as certified contracts with Gauntlet tests; the measured perf tax is accepted explicitly.",
        early_warning: "reproducibility test failures per release",
        threshold: "determinism tax exceeding ~15% on dominant kernels",
        owner: "frankensim-epic-flywheel-lmp4.6",
        instrumented: false,
    },
    Risk {
        id: RiskId::R3,
        name: "Stable entity identity",
        description: "Semantic diff and merge require persistent IDs through topology-changing edits (the CAD topological-naming problem).",
        mitigation: "FrankenSim owns its kernel; IDs are first-class from day one; edits are ledgered ops that transform IDs explicitly.",
        early_warning: "fraction of diffs falling back to unattributed geometric comparison",
        threshold: "a material fraction of diffs cannot attribute to a causal edit",
        owner: "frankensim-epic-flywheel-lmp4.10",
        instrumented: false,
    },
    Risk {
        id: RiskId::R4,
        name: "Loose sensitivity bounds",
        description: "Interval derivatives on nonlinear ops may be so pessimistic that the recompute frontier balloons to the whole DAG.",
        mitigation: "Graceful degradation to hash memoization; adjoint-sharpened bounds where Proposal 1 is live; per-op skip-yield telemetry.",
        early_warning: "skip-yield (fraction of DAG certifiably skipped)",
        threshold: "skip-yield delivers <2x median wall-clock vs hash memoization",
        owner: "frankensim-epic-flywheel-lmp4.7",
        instrumented: false,
    },
    Risk {
        id: RiskId::R5,
        name: "Spectral-gap fragility",
        description: "The fixable/structural triage (and merge adjudication) degrades in ill-conditioned regions.",
        mitigation: "Spectral-health monitoring with mandatory low-confidence propagation into merge outputs; the Gauntlet refusal-is-the-pass suite.",
        early_warning: "gap-collapse incidence per assembly class",
        threshold: "gap collapse observed outside synthetic cases at volume",
        owner: "frankensim-epic-selfknow-knh1.3",
        instrumented: false,
    },
    Risk {
        id: RiskId::R6,
        name: "Restriction-map quality",
        description: "The sheaf propagates garbage with certificates attached if trace/conversion operators are inaccurate.",
        mitigation: "The Proposal 7 conformance suite (functoriality, adjoint consistency, MMS tolerance honesty) applied to first-party converters with third-party severity.",
        early_warning: "conformance-tier distribution of converters on the hot path",
        threshold: "hot-path converters cluster in low conformance tiers",
        owner: "frankensim-epic-gtm-jwq8.2",
        instrumented: false,
    },
    Risk {
        id: RiskId::R7,
        name: "Model-form composition",
        description: "Estimated-color quantities do not compose with the clean algebra of verified bounds.",
        mitigation: "The type system's laundering refusal; discrepancy-probe maps as the empirical substitute; the weakest-input rule on headlines.",
        early_warning: "audit rate of estimated-color claims in decision-critical positions",
        threshold: "estimated-color claims silently drive decisions without probe maps",
        owner: "frankensim-epic-epistype-qmao.2",
        instrumented: false,
    },
    Risk {
        id: RiskId::R8,
        name: "Registration well-posedness",
        description: "Scan-to-design alignment error can exceed the deviations being certified.",
        mitigation: "Design-for-verification fiducials pushed upstream; point-sensor assimilation as the registration-free fallback.",
        early_warning: "registration-uncertainty-to-signal ratio per part class",
        threshold: "registration uncertainty exceeds the geometric deviations being certified",
        owner: "frankensim-epic-coupling-bk0o.4",
        instrumented: false,
    },
    Risk {
        id: RiskId::R9,
        name: "Standards-body latency",
        description: "Machine-checkable evidence may sit ahead of what auditors will engage with.",
        mitigation: "The vocabulary crosswalk speaks their language; the package doubles as internal QA and B2B diligence collateral, so investment is not stranded.",
        early_warning: "auditor engagement rate in the first regulated-vertical cycle",
        threshold: "no auditor engages the machine-checkable format even as supplementary evidence",
        owner: "frankensim-epic-epistype-qmao.9",
        instrumented: false,
    },
    Risk {
        id: RiskId::R10,
        name: "Breadth death",
        description: "Nineteen proposals executed as nineteen parallel programs is an obituary.",
        mitigation: "The governance rules: one research bet at a time, quarterly kill enforcement, phase gating.",
        early_warning: "headcount-weighted work-in-progress outside the current phase",
        threshold: "WIP outside the current phase grows materially",
        owner: "frankensim-epic-addendum-xpck.1",
        instrumented: false,
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

/// The result of auditing the register for governance completeness.
#[derive(Debug, Clone, PartialEq)]
pub struct RiskAudit {
    /// Total risks.
    pub total: usize,
    /// Risks that have BOTH a non-empty early-warning metric and an owner.
    pub complete: usize,
    /// Risks whose early-warning metric is actually instrumented.
    pub instrumented: usize,
    /// `(risk, reason)` for every incomplete entry.
    pub gaps: Vec<(RiskId, &'static str)>,
}

impl RiskAudit {
    /// Is the register complete (every risk has a metric AND an owner)?
    #[must_use]
    pub fn ok(&self) -> bool {
        self.gaps.is_empty()
    }
}

/// Audit the canonical register (see [`audit_slice`]).
#[must_use]
pub fn audit() -> RiskAudit {
    audit_slice(&REGISTER)
}

/// Audit an arbitrary risk slice: every risk must carry an early-warning
/// metric and an owner (Governance Rule 2 — a risk with no measurement is
/// unmanaged). Also counts how many early-warning metrics are instrumented
/// (the standing gap between "declared" and "watched").
#[must_use]
pub fn audit_slice(risks: &[Risk]) -> RiskAudit {
    let mut gaps = Vec::new();
    let mut complete = 0usize;
    let mut instrumented = 0usize;
    for r in risks {
        let mut ok = true;
        if r.early_warning.trim().is_empty() {
            gaps.push((r.id, "missing early-warning metric"));
            ok = false;
        }
        if r.owner.trim().is_empty() {
            gaps.push((r.id, "missing owner"));
            ok = false;
        }
        if ok {
            complete += 1;
        }
        if r.instrumented {
            instrumented += 1;
        }
    }
    RiskAudit {
        total: risks.len(),
        complete,
        instrumented,
        gaps,
    }
}

/// Escape a string for embedding in JSON.
pub(crate) fn json_escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len() + 2);
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\t' => out.push_str("\\t"),
            '\r' => out.push_str("\\r"),
            _ => out.push(c),
        }
    }
    out
}

/// Emit the register as a machine-readable JSON array (one object per risk),
/// for dashboards and CI gates. Deterministic (risks in order).
#[must_use]
pub fn to_json() -> String {
    use core::fmt::Write as _;
    let mut out = String::from("[");
    for (i, r) in REGISTER.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        write!(
            out,
            "{{\"id\":\"{}\",\"name\":\"{}\",\"early_warning\":\"{}\",\"threshold\":\"{}\",\"owner\":\"{}\",\"instrumented\":{},\"mitigation\":\"{}\"}}",
            r.id.code(),
            json_escape(r.name),
            json_escape(r.early_warning),
            json_escape(r.threshold),
            json_escape(r.owner),
            r.instrumented,
            json_escape(r.mitigation),
        )
        .expect("writing to a String is infallible");
    }
    out.push(']');
    out
}
