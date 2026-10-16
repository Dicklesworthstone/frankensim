//! The external-adapter policy ratification (bead f85xj.11.1): the explicit
//! ruling that resolves the tension between "no foreign kernels in
//! production" and real-world CAD/CAE adoption — recorded the way this
//! project records claims, with the considered options, the chosen ruling,
//! the verbatim mission amendment, the non-negotiable trust invariants, named
//! mechanically-evaluable falsifiers, and a review date.
//!
//! The ruling is OFFICIAL QUARANTINED ADAPTERS: FrankenSim may ship optional
//! adapter binaries wrapping foreign tools as SEPARATE out-of-process
//! executables with SEPARATE distribution, never linked into the production
//! dependency graph; adapter OUTPUT enters the workspace only through the
//! fs-io quarantine boundary with Estimate-only authority and the adapter
//! identity + version recorded in receipts. In-process foreign kernels (FFI
//! plugins) remain forbidden in production paths, and `check-deps`
//! enforcement is unchanged — the ruling legalizes a distribution channel,
//! not a dependency edge.
//!
//! Authority note: the bead flagged this decision as requiring human
//! sign-off because it amends mission text. The project owner explicitly
//! delegated that sign-off on 2026-07-23 ("use your best judgment on all
//! that stuff and proceed"), and the record retains that delegation verbatim
//! so the authority chain is auditable rather than implied.
//!
//! The record is fail-closed like the vertical ratification: the accessor
//! re-validates structural completeness before handing the record out, and a
//! conformance test (`tests/adapter_policy.rs`) enforces the three-way
//! agreement between this record, the AGENTS.md mission text, and the xtask
//! `check-deps` policy language, so the three cannot drift apart silently.

use crate::json_escape;
use crate::ratification::Falsifier;
use core::fmt::Write as _;

/// Stable decision-record id cited by AGENTS.md, xtask, and downstream beads.
pub const ADAPTER_POLICY_ID: &str = "ADPT-2026-07";

/// One considered (and possibly rejected) policy option.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AdapterPolicyOption {
    /// Stable option slug.
    pub id: &'static str,
    /// Human-facing label.
    pub label: &'static str,
    /// Why the option was rejected (empty only for the chosen option).
    pub rejection_reason: &'static str,
}

impl AdapterPolicyOption {
    /// Is every structural field populated?
    #[must_use]
    pub fn is_complete(self) -> bool {
        !self.id.trim().is_empty() && !self.label.trim().is_empty()
    }
}

/// One trust invariant that holds under EVERY option, chosen or not.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TrustInvariant {
    /// Stable invariant id.
    pub id: &'static str,
    /// The non-negotiable statement.
    pub statement: &'static str,
}

impl TrustInvariant {
    /// Is every field populated?
    #[must_use]
    pub fn is_complete(self) -> bool {
        !self.id.trim().is_empty() && !self.statement.trim().is_empty()
    }
}

/// The external-adapter policy decision record.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AdapterPolicyRatification {
    /// Stable record id.
    pub id: &'static str,
    /// Decision date.
    pub decided_on: &'static str,
    /// The authority chain, including the recorded human delegation.
    pub authority: &'static str,
    /// Chosen option slug.
    pub chosen_option: &'static str,
    /// Every considered option, chosen first.
    pub options: &'static [AdapterPolicyOption],
    /// The operative clauses of the ruling.
    pub ruling: &'static [&'static str],
    /// The verbatim mission-amendment heading now present in AGENTS.md; the
    /// conformance test enforces its literal presence.
    pub mission_amendment_heading: &'static str,
    /// Key amendment clauses that must appear (whitespace-normalized) in
    /// AGENTS.md; the conformance test enforces each.
    pub mission_amendment_clauses: &'static [&'static str],
    /// Trust invariants that remain non-negotiable under any option.
    pub invariants: &'static [TrustInvariant],
    /// Named mechanically-evaluable falsifiers / review triggers.
    pub falsifiers: &'static [Falsifier],
    /// Next scheduled review.
    pub review_due: &'static str,
    /// Downstream beads gated on this record.
    pub downstream_gates: &'static [&'static str],
}

/// Why the adapter-policy record refuses to stand.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdapterPolicyError {
    /// A required text field is empty.
    EmptyField {
        /// Which field.
        field: &'static str,
    },
    /// Fewer than two options were considered: a ruling with no rejected
    /// alternative is not a decision.
    MissingAlternatives,
    /// An option entry is structurally incomplete.
    IncompleteOption {
        /// The offending option id (or its ordinal label when empty).
        id: &'static str,
    },
    /// The chosen option is not among the considered options.
    ChosenOptionUnlisted,
    /// A rejected option carries no rejection reason.
    UnreasonedRejection {
        /// The offending option id.
        id: &'static str,
    },
    /// The record names no ruling clause.
    MissingRuling,
    /// The record names no amendment clause for the mission text.
    MissingAmendment,
    /// The record names no trust invariant.
    MissingInvariants,
    /// A trust invariant is structurally incomplete.
    IncompleteInvariant {
        /// The offending invariant id.
        id: &'static str,
    },
    /// The record names no falsifier at all.
    MissingFalsifiers,
    /// A named falsifier is structurally incomplete.
    IncompleteFalsifier {
        /// The offending falsifier id.
        id: &'static str,
    },
    /// The record gates no downstream work: a ruling nobody consumes is
    /// governance theater.
    MissingDownstreamGates,
}

impl core::fmt::Display for AdapterPolicyError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::EmptyField { field } => {
                write!(formatter, "adapter policy field {field} is empty")
            }
            Self::MissingAlternatives => {
                write!(formatter, "adapter policy considers fewer than two options")
            }
            Self::IncompleteOption { id } => {
                write!(formatter, "adapter policy option {id:?} is incomplete")
            }
            Self::ChosenOptionUnlisted => {
                write!(
                    formatter,
                    "chosen option is not among the considered options"
                )
            }
            Self::UnreasonedRejection { id } => {
                write!(
                    formatter,
                    "rejected option {id:?} carries no rejection reason"
                )
            }
            Self::MissingRuling => write!(formatter, "adapter policy has no ruling clauses"),
            Self::MissingAmendment => {
                write!(
                    formatter,
                    "adapter policy names no mission-amendment clause"
                )
            }
            Self::MissingInvariants => {
                write!(formatter, "adapter policy names no trust invariant")
            }
            Self::IncompleteInvariant { id } => {
                write!(formatter, "trust invariant {id:?} is incomplete")
            }
            Self::MissingFalsifiers => write!(formatter, "adapter policy names no falsifier"),
            Self::IncompleteFalsifier { id } => {
                write!(formatter, "falsifier {id:?} is incomplete")
            }
            Self::MissingDownstreamGates => {
                write!(formatter, "adapter policy gates no downstream work")
            }
        }
    }
}

impl core::error::Error for AdapterPolicyError {}

const OPTIONS: &[AdapterPolicyOption] = &[
    AdapterPolicyOption {
        id: "official-quarantined-adapters",
        label: "Official quarantined adapters (CHOSEN)",
        rejection_reason: "",
    },
    AdapterPolicyOption {
        id: "pure-core-only",
        label: "Pure core only; users bring their own conversion tooling",
        rejection_reason: "Production CAD semantics (AP242 assemblies, trimmed NURBS, PMI) \
                           cannot realistically be rebuilt from scratch quickly, so refusing to \
                           distribute any adapter makes the missing tooling every adopter's \
                           first experience: an adoption death sentence that silently defeats \
                           the wedge's cycle-time mission while pretending to be principled.",
    },
    AdapterPolicyOption {
        id: "ffi-plugin-ecosystem",
        label: "In-process FFI plugin ecosystem with backend trust contracts",
        rejection_reason: "Highest adoption, highest trust erosion: an in-process foreign \
                           kernel shares our address space, so its memory-unsafety, panics, and \
                           nondeterminism become indistinguishable from ours, the workspace's \
                           #![deny(unsafe_code)] and determinism-class claims stop being \
                           auditable at the crate boundary, and the no-laundering rule loses \
                           its enforcement point. A future re-evaluation requires a full \
                           backend trust-contract design, not ad-hoc FFI.",
    },
];

const RULING: &[&str] = &[
    "FrankenSim may ship optional official adapter binaries that wrap foreign tools \
     (OpenCASCADE/gmsh-class kernels) as separate out-of-process executables with separate \
     distribution.",
    "Adapters are never build-time or link-time dependencies of any workspace crate; the \
     production runtime dependency graph stays Franken-only and check-deps enforcement is \
     unchanged.",
    "Adapter output enters the workspace only through the existing fs-io quarantine boundary: \
     evidence-bearing receipts, Estimate-only authority, and the adapter identity plus version \
     recorded in the receipt.",
    "An adapter can never mint authority by itself; promotion beyond Estimate requires the \
     independent verifier machinery, exactly as for any other external input.",
    "In-process foreign kernels (FFI plugins) remain forbidden in production paths.",
];

const MISSION_AMENDMENT_HEADING: &str = "### External-Adapter Ruling (ADPT-2026-07)";

const MISSION_AMENDMENT_CLAUSES: &[&str] = &[
    "separate out-of-process executables with separate distribution",
    "never build-time or link-time dependencies of any workspace crate",
    "only through the fs-io quarantine boundary",
    "Estimate-only authority",
    "adapter identity and version recorded in the receipt",
    "In-process foreign kernels (FFI plugins) remain forbidden in production paths",
];

const INVARIANTS: &[TrustInvariant] = &[
    TrustInvariant {
        id: "quarantine-receipts",
        statement: "Every external artifact enters through evidence-bearing quarantine \
                    receipts; there is no receiptless side door under any option.",
    },
    TrustInvariant {
        id: "no-authority-minting",
        statement: "No adapter, importer, or conversion step can raise an artifact above \
                    Estimate authority by itself; promotion requires the independent verifier \
                    machinery.",
    },
    TrustInvariant {
        id: "content-addressing",
        statement: "Adapter inputs and outputs are content-addressed so receipts bind exact \
                    bytes, not filenames or timestamps.",
    },
    TrustInvariant {
        id: "process-isolation",
        statement: "Foreign code runs in a separate OS process with a typed file/stream \
                    boundary; it never shares the workspace address space in production.",
    },
    TrustInvariant {
        id: "recorded-provenance",
        statement: "The adapter identity and version travel in the receipt, so a defect in one \
                    adapter release is traceable to exactly the artifacts it produced.",
    },
];

const FALSIFIERS: &[Falsifier] = &[
    Falsifier {
        id: "adapter-authority-leak",
        statement: "An adapter-origin artifact reaches Verified or Validated authority without \
                    an independent verifier receipt.",
        measurement: "Quarterly audit of fs-io receipts and checker admissions for \
                      adapter-origin artifacts.",
        threshold: "Any single occurrence forces an immediate policy review.",
    },
    Falsifier {
        id: "quarantine-bypass",
        statement: "A workspace crate links a foreign kernel or embeds an adapter binary \
                    inside the production dependency graph.",
        measurement: "check-deps plus the source-manifest external-tool surface review.",
        threshold: "Any runtime dependency edge or embedded foreign binary.",
    },
    Falsifier {
        id: "adoption-pressure",
        statement: "The out-of-process boundary demonstrably blocks a committed design partner \
                    on cycle time or capability.",
        measurement: "Wedge cycle-time metrics and partner-feedback review.",
        threshold: "Two consecutive quarterly reviews citing the adapter boundary as the \
                    binding constraint triggers a re-evaluation of the FFI option with a full \
                    backend trust-contract design.",
    },
    Falsifier {
        id: "maintenance-starvation",
        statement: "Adapter maintenance starves core program lanes.",
        measurement: "Share of program velocity spent on adapter deliverables per quarter.",
        threshold: "More than twenty percent for two consecutive quarters forces re-scoping \
                    toward pure-core.",
    },
];

const ADAPTER_POLICY: AdapterPolicyRatification = AdapterPolicyRatification {
    id: ADAPTER_POLICY_ID,
    decided_on: "2026-07-23",
    authority: "Bead f85xj.11.1 flagged this ruling as requiring human sign-off because it \
                amends mission text. The project owner explicitly delegated the decision on \
                2026-07-23 (\"use your best judgment on all that stuff and proceed\"), and the \
                executing session ratified the bead's own recommended option under that \
                delegation.",
    chosen_option: "official-quarantined-adapters",
    options: OPTIONS,
    ruling: RULING,
    mission_amendment_heading: MISSION_AMENDMENT_HEADING,
    mission_amendment_clauses: MISSION_AMENDMENT_CLAUSES,
    invariants: INVARIANTS,
    falsifiers: FALSIFIERS,
    review_due: "2026-10-23",
    downstream_gates: &[
        "frankensim-extreal-program-f85xj.6.13",
        "frankensim-extreal-program-f85xj.11.5",
    ],
};

fn validate(record: &AdapterPolicyRatification) -> Result<(), AdapterPolicyError> {
    for (field, value) in [
        ("id", record.id),
        ("decided_on", record.decided_on),
        ("authority", record.authority),
        ("chosen_option", record.chosen_option),
        (
            "mission_amendment_heading",
            record.mission_amendment_heading,
        ),
        ("review_due", record.review_due),
    ] {
        if value.trim().is_empty() {
            return Err(AdapterPolicyError::EmptyField { field });
        }
    }
    if record.options.len() < 2 {
        return Err(AdapterPolicyError::MissingAlternatives);
    }
    let mut chosen_listed = false;
    for option in record.options {
        if !option.is_complete() {
            return Err(AdapterPolicyError::IncompleteOption {
                id: if option.id.trim().is_empty() {
                    "(unnamed)"
                } else {
                    option.id
                },
            });
        }
        if option.id == record.chosen_option {
            chosen_listed = true;
        } else if option.rejection_reason.trim().is_empty() {
            return Err(AdapterPolicyError::UnreasonedRejection { id: option.id });
        }
    }
    if !chosen_listed {
        return Err(AdapterPolicyError::ChosenOptionUnlisted);
    }
    if record.ruling.is_empty() || record.ruling.iter().any(|clause| clause.trim().is_empty()) {
        return Err(AdapterPolicyError::MissingRuling);
    }
    if record.mission_amendment_clauses.is_empty()
        || record
            .mission_amendment_clauses
            .iter()
            .any(|clause| clause.trim().is_empty())
    {
        return Err(AdapterPolicyError::MissingAmendment);
    }
    if record.invariants.is_empty() {
        return Err(AdapterPolicyError::MissingInvariants);
    }
    for invariant in record.invariants {
        if !invariant.is_complete() {
            return Err(AdapterPolicyError::IncompleteInvariant { id: invariant.id });
        }
    }
    if record.falsifiers.is_empty() {
        return Err(AdapterPolicyError::MissingFalsifiers);
    }
    for falsifier in record.falsifiers {
        if !falsifier.is_complete() {
            return Err(AdapterPolicyError::IncompleteFalsifier { id: falsifier.id });
        }
    }
    if record.downstream_gates.is_empty()
        || record
            .downstream_gates
            .iter()
            .any(|gate| gate.trim().is_empty())
    {
        return Err(AdapterPolicyError::MissingDownstreamGates);
    }
    Ok(())
}

/// The ratified external-adapter policy, re-validated fail-closed on every
/// access so an incomplete edit cannot stand on stale authority.
///
/// # Errors
/// A structural completeness refusal naming the offending field.
pub fn adapter_policy() -> Result<&'static AdapterPolicyRatification, AdapterPolicyError> {
    validate(&ADAPTER_POLICY)?;
    Ok(&ADAPTER_POLICY)
}

/// Deterministic JSON rendering of the validated record for reports and
/// dashboards.
///
/// # Errors
/// The same structural refusals as [`adapter_policy`].
pub fn adapter_policy_json() -> Result<String, AdapterPolicyError> {
    let record = adapter_policy()?;
    let mut out = String::from("{\"adapter_policy\":{");
    let _ = write!(
        out,
        "\"id\":\"{}\",\"decided_on\":\"{}\",\"chosen_option\":\"{}\",\"authority\":\"{}\",\"review_due\":\"{}\"",
        json_escape(record.id),
        json_escape(record.decided_on),
        json_escape(record.chosen_option),
        json_escape(record.authority),
        json_escape(record.review_due),
    );
    out.push_str(",\"options\":[");
    for (index, option) in record.options.iter().enumerate() {
        if index > 0 {
            out.push(',');
        }
        let _ = write!(
            out,
            "{{\"id\":\"{}\",\"label\":\"{}\",\"rejection_reason\":\"{}\"}}",
            json_escape(option.id),
            json_escape(option.label),
            json_escape(option.rejection_reason),
        );
    }
    out.push_str("],\"ruling\":[");
    for (index, clause) in record.ruling.iter().enumerate() {
        if index > 0 {
            out.push(',');
        }
        let _ = write!(out, "\"{}\"", json_escape(clause));
    }
    out.push_str("],\"invariants\":[");
    for (index, invariant) in record.invariants.iter().enumerate() {
        if index > 0 {
            out.push(',');
        }
        let _ = write!(
            out,
            "{{\"id\":\"{}\",\"statement\":\"{}\"}}",
            json_escape(invariant.id),
            json_escape(invariant.statement),
        );
    }
    out.push_str("],\"falsifiers\":[");
    for (index, falsifier) in record.falsifiers.iter().enumerate() {
        if index > 0 {
            out.push(',');
        }
        let _ = write!(
            out,
            "{{\"id\":\"{}\",\"statement\":\"{}\",\"measurement\":\"{}\",\"threshold\":\"{}\"}}",
            json_escape(falsifier.id),
            json_escape(falsifier.statement),
            json_escape(falsifier.measurement),
            json_escape(falsifier.threshold),
        );
    }
    out.push_str("],\"downstream_gates\":[");
    for (index, gate) in record.downstream_gates.iter().enumerate() {
        if index > 0 {
            out.push(',');
        }
        let _ = write!(out, "\"{}\"", json_escape(gate));
    }
    out.push_str("]}}");
    Ok(out)
}
