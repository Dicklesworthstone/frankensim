//! Generated requirement-to-evidence traceability ledger.
//!
//! The extension charter's Appendix E is a seed, not a status dashboard. This
//! module turns that seed into typed, reviewable registry data and derives the
//! machine-readable ledger from it. Callers that ingest equivalent rows from
//! Beads or contracts use the same fail-closed validator and renderer.

use core::fmt::Write as _;
use std::collections::BTreeSet;

use crate::json_escape;

/// Stable schema tag for generated traceability artifacts.
pub const TRACEABILITY_SCHEMA: &str = "frankensim-requirement-traceability-v1";

/// Authority label emitted by this pure declaration registry.
///
/// Scientific promotion requires a separate source-bound admission receipt;
/// this generator intentionally cannot mint one.
pub const TRACEABILITY_AUTHORITY: &str = "declaration-only";

/// Number of executable proof obligations in the extension charter.
pub const PROOF_OBLIGATION_COUNT: usize = 25;

/// Number of B1-B14 and RQ-* rows in the extension charter seed.
pub const REQUIREMENT_COUNT: usize = 30;

/// Maximum requirement rows accepted from one source snapshot.
pub const MAX_REQUIREMENT_ROWS: usize = 256;

/// Maximum UTF-8 bytes accepted in one scalar source field.
pub const MAX_TRACEABILITY_FIELD_BYTES: usize = 16 * 1024;

/// Maximum PO links on one requirement row.
pub const MAX_REQUIREMENT_PO_LINKS: usize = PROOF_OBLIGATION_COUNT;

/// Maximum owning Beads on one proof-obligation definition.
pub const MAX_PROOF_OBLIGATION_OWNERS: usize = 16;

const REQUIRED_REQUIREMENT_IDS: [&str; REQUIREMENT_COUNT] = [
    "B1",
    "B2",
    "B3",
    "B4",
    "B5",
    "B6",
    "B7",
    "B8",
    "B9",
    "B10",
    "B11",
    "B12",
    "B13",
    "B14",
    "RQ-ROLL",
    "RQ-GEAR",
    "RQ-FRICTION",
    "RQ-CONSTITUTIVE",
    "RQ-DENSITY",
    "RQ-MECHMAT",
    "RQ-ELEC",
    "RQ-MAG",
    "RQ-PHASE",
    "RQ-FLUID",
    "RQ-PERMEATE",
    "RQ-WET",
    "RQ-MOTORGEN",
    "RQ-ICE",
    "RQ-ACOUSTIC",
    "RQ-ACTIVE",
];

const VALID_STATUSES: [&str; 6] = [
    "proposed",
    "in_progress",
    "implemented",
    "proof-pending",
    "refused",
    "retired",
];

/// One source row for the generated requirement-to-evidence ledger.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RequirementRow<'a> {
    /// Stable requirement identifier (`B1`..`B14` or `RQ-*`).
    pub requirement_id: &'a str,
    /// Capability or scientific property the requirement exists to deliver.
    pub capability_property: &'a str,
    /// Concrete current blocker, never an aspirational restatement.
    pub blocker: &'a str,
    /// Owning crate, type, Bead, or generated artifact.
    pub owner_artifact: &'a str,
    /// Earlier phase capability that must exist before this row can execute.
    pub prerequisite_phase: &'a str,
    /// Phase or milestone whose gate consumes this row.
    pub milestone: &'a str,
    /// Flagship or cross-cutting vertical slice forced by the row.
    pub flagship: &'a str,
    /// Benchmark, dataset, or executable battery that supplies evidence.
    pub benchmark_data: &'a str,
    /// IDs in the complete [`proof_obligations`] registry.
    pub proof_obligations: &'a [&'a str],
    /// Honest limit beyond which the row grants no claim authority.
    pub claim_boundary: &'a str,
    /// Declarative lifecycle state supplied by the source registry.
    pub status: &'a str,
}

/// One entry in the complete PO-1 through PO-25 index.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ProofObligation<'a> {
    /// Stable `PO-n` identifier.
    pub id: &'a str,
    /// Executable obligation summary. The owner Beads retain full detail.
    pub summary: &'a str,
    /// Beads that own the executable evidence for this obligation.
    pub owner_beads: &'a [&'a str],
}

/// A field that can make a traceability source row orphaned or unrenderable.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum TraceabilityField {
    /// The ledger scope itself.
    Ledger,
    /// `requirement_id`.
    RequirementId,
    /// `capability_property`.
    CapabilityProperty,
    /// `blocker`.
    Blocker,
    /// `owner_artifact`.
    OwnerArtifact,
    /// `prerequisite_phase`.
    PrerequisitePhase,
    /// `milestone`.
    Milestone,
    /// `flagship`.
    Flagship,
    /// `benchmark_data`.
    BenchmarkData,
    /// `proof_obligations` or the PO index.
    ProofObligation,
    /// `claim_boundary`.
    ClaimBoundary,
    /// `status`.
    Status,
}

impl TraceabilityField {
    /// Stable field name used in structured diagnostics.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Ledger => "ledger",
            Self::RequirementId => "requirement_id",
            Self::CapabilityProperty => "capability_property",
            Self::Blocker => "blocker",
            Self::OwnerArtifact => "owner_artifact",
            Self::PrerequisitePhase => "prerequisite_phase",
            Self::Milestone => "milestone",
            Self::Flagship => "flagship",
            Self::BenchmarkData => "benchmark_data",
            Self::ProofObligation => "proof_obligations",
            Self::ClaimBoundary => "claim_boundary",
            Self::Status => "status",
        }
    }
}

/// Structured fail-closed diagnostic emitted by traceability generation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TraceabilityDiagnostic {
    /// Requirement that failed, or a synthetic registry scope such as
    /// `<proof-obligation-index>`.
    pub requirement_id: String,
    /// Exact missing or inconsistent field.
    pub field: TraceabilityField,
    /// Actionable explanation of the failure.
    pub reason: String,
}

/// Completeness report for a proposed traceability ledger.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TraceabilityAudit {
    /// Number of supplied requirement rows.
    pub total: usize,
    /// Rows with no row-local diagnostic.
    pub complete: usize,
    /// Every structured refusal, in deterministic source order.
    pub diagnostics: Vec<TraceabilityDiagnostic>,
}

impl TraceabilityAudit {
    /// Whether the non-empty source and complete PO index can be rendered.
    #[must_use]
    pub fn ok(&self) -> bool {
        self.total > 0 && self.complete == self.total && self.diagnostics.is_empty()
    }
}

fn diagnostic(
    requirement_id: impl Into<String>,
    field: TraceabilityField,
    reason: impl Into<String>,
) -> TraceabilityDiagnostic {
    TraceabilityDiagnostic {
        requirement_id: requirement_id.into(),
        field,
        reason: reason.into(),
    }
}

fn valid_po_id(id: &str) -> Option<u8> {
    let digits = id.strip_prefix("PO-")?;
    if digits.is_empty()
        || digits.starts_with('0')
        || !digits.bytes().all(|byte| byte.is_ascii_digit())
    {
        return None;
    }
    let value = digits.parse::<u8>().ok()?;
    (1..=PROOF_OBLIGATION_COUNT as u8)
        .contains(&value)
        .then_some(value)
}

fn field_diagnostics(
    diagnostics: &mut Vec<TraceabilityDiagnostic>,
    scope: &str,
    field: TraceabilityField,
    value: &str,
    missing_reason: &'static str,
) {
    if value.len() > MAX_TRACEABILITY_FIELD_BYTES {
        diagnostics.push(diagnostic(
            scope,
            field,
            format!(
                "field {} is {} bytes; maximum is {}",
                field.name(),
                value.len(),
                MAX_TRACEABILITY_FIELD_BYTES
            ),
        ));
    } else if value.trim().is_empty() {
        diagnostics.push(diagnostic(scope, field, missing_reason));
    }
}

fn requirement_scope(id: &str) -> &str {
    if id.len() > MAX_TRACEABILITY_FIELD_BYTES {
        "<oversized-requirement-id>"
    } else if id.trim().is_empty() {
        "<orphaned-requirement>"
    } else {
        id
    }
}

/// Audit arbitrary rows from code, Beads, or contract adapters before output.
///
/// Every declared field is mandatory. The audit also requires an exact,
/// unique PO-1..PO-25 index and refuses duplicate requirements or dangling PO
/// links. No partial ledger is returned by [`generate_traceability_ledger`].
#[must_use]
#[allow(clippy::too_many_lines)] // exhaustive field and registry refusals remain one audit
pub fn audit_traceability(
    rows: &[RequirementRow<'_>],
    obligations: &[ProofObligation<'_>],
) -> TraceabilityAudit {
    let mut diagnostics = Vec::new();
    if rows.len() > MAX_REQUIREMENT_ROWS {
        diagnostics.push(diagnostic(
            "<ledger>",
            TraceabilityField::Ledger,
            format!(
                "source contains {} requirements; maximum is {MAX_REQUIREMENT_ROWS}",
                rows.len()
            ),
        ));
        return TraceabilityAudit {
            total: rows.len(),
            complete: 0,
            diagnostics,
        };
    }
    if obligations.len() > PROOF_OBLIGATION_COUNT {
        diagnostics.push(diagnostic(
            "<proof-obligation-index>",
            TraceabilityField::ProofObligation,
            format!(
                "index contains {} definitions; exact maximum is {PROOF_OBLIGATION_COUNT}",
                obligations.len()
            ),
        ));
        return TraceabilityAudit {
            total: rows.len(),
            complete: 0,
            diagnostics,
        };
    }
    let mut known_obligations = BTreeSet::new();

    for obligation in obligations {
        let scope = if obligation.id.len() > MAX_TRACEABILITY_FIELD_BYTES {
            "<oversized-proof-obligation-id>"
        } else if obligation.id.trim().is_empty() {
            "<proof-obligation-index>"
        } else {
            obligation.id
        };
        if obligation.id.len() > MAX_TRACEABILITY_FIELD_BYTES {
            diagnostics.push(diagnostic(
                scope,
                TraceabilityField::ProofObligation,
                format!(
                    "proof-obligation id is {} bytes; maximum is {MAX_TRACEABILITY_FIELD_BYTES}",
                    obligation.id.len()
                ),
            ));
            continue;
        }
        if obligation.id.trim().is_empty() {
            diagnostics.push(diagnostic(
                scope,
                TraceabilityField::ProofObligation,
                "proof-obligation definition is missing its id",
            ));
            continue;
        }
        if valid_po_id(obligation.id).is_none() {
            diagnostics.push(diagnostic(
                scope,
                TraceabilityField::ProofObligation,
                format!(
                    "proof-obligation id {:?} is outside the closed PO-1..PO-25 registry",
                    obligation.id
                ),
            ));
        }
        if !known_obligations.insert(obligation.id) {
            diagnostics.push(diagnostic(
                scope,
                TraceabilityField::ProofObligation,
                format!("duplicate proof-obligation definition {:?}", obligation.id),
            ));
        }
        if obligation.summary.len() > MAX_TRACEABILITY_FIELD_BYTES {
            diagnostics.push(diagnostic(
                scope,
                TraceabilityField::ProofObligation,
                format!(
                    "proof-obligation summary is {} bytes; maximum is {MAX_TRACEABILITY_FIELD_BYTES}",
                    obligation.summary.len()
                ),
            ));
        } else if obligation.summary.trim().is_empty() {
            diagnostics.push(diagnostic(
                scope,
                TraceabilityField::ProofObligation,
                "proof-obligation definition is missing its executable summary",
            ));
        }
        if obligation.owner_beads.is_empty() {
            diagnostics.push(diagnostic(
                scope,
                TraceabilityField::OwnerArtifact,
                "proof-obligation definition has no owning Bead",
            ));
        } else if obligation.owner_beads.len() > MAX_PROOF_OBLIGATION_OWNERS {
            diagnostics.push(diagnostic(
                scope,
                TraceabilityField::OwnerArtifact,
                format!(
                    "proof-obligation has {} owners; maximum is {MAX_PROOF_OBLIGATION_OWNERS}",
                    obligation.owner_beads.len()
                ),
            ));
        } else {
            let mut owners = BTreeSet::new();
            for owner in obligation.owner_beads {
                if owner.len() > MAX_TRACEABILITY_FIELD_BYTES {
                    diagnostics.push(diagnostic(
                        scope,
                        TraceabilityField::OwnerArtifact,
                        format!(
                            "owner Bead is {} bytes; maximum is {MAX_TRACEABILITY_FIELD_BYTES}",
                            owner.len()
                        ),
                    ));
                } else if owner.trim().is_empty() {
                    diagnostics.push(diagnostic(
                        scope,
                        TraceabilityField::OwnerArtifact,
                        "proof-obligation definition has no owning Bead",
                    ));
                } else if !owners.insert(*owner) {
                    diagnostics.push(diagnostic(
                        scope,
                        TraceabilityField::OwnerArtifact,
                        format!("duplicate proof-obligation owner {owner:?}"),
                    ));
                }
            }
        }
    }

    for number in 1..=PROOF_OBLIGATION_COUNT {
        let id = format!("PO-{number}");
        if !known_obligations.contains(id.as_str()) {
            diagnostics.push(diagnostic(
                "<proof-obligation-index>",
                TraceabilityField::ProofObligation,
                format!("complete index is missing {id}"),
            ));
        }
    }

    if rows.is_empty() {
        diagnostics.push(diagnostic(
            "<ledger>",
            TraceabilityField::Ledger,
            "traceability generation refuses an empty requirement scope",
        ));
    }

    let mut requirement_ids = BTreeSet::new();
    for row in rows {
        let scope = requirement_scope(row.requirement_id);
        let required = [
            (
                TraceabilityField::RequirementId,
                row.requirement_id,
                "requirement has no stable id",
            ),
            (
                TraceabilityField::CapabilityProperty,
                row.capability_property,
                "requirement has no capability/property",
            ),
            (
                TraceabilityField::Blocker,
                row.blocker,
                "requirement has no concrete blocker",
            ),
            (
                TraceabilityField::OwnerArtifact,
                row.owner_artifact,
                "requirement has no owner or artifact route",
            ),
            (
                TraceabilityField::PrerequisitePhase,
                row.prerequisite_phase,
                "requirement has no prerequisite phase",
            ),
            (
                TraceabilityField::Milestone,
                row.milestone,
                "requirement has no executable milestone gate",
            ),
            (
                TraceabilityField::Flagship,
                row.flagship,
                "requirement has no forcing flagship or vertical slice",
            ),
            (
                TraceabilityField::BenchmarkData,
                row.benchmark_data,
                "requirement has no benchmark/data evidence route",
            ),
            (
                TraceabilityField::ClaimBoundary,
                row.claim_boundary,
                "requirement has no honest claim boundary",
            ),
            (
                TraceabilityField::Status,
                row.status,
                "requirement has no lifecycle status",
            ),
        ];
        for (field, value, reason) in required {
            field_diagnostics(&mut diagnostics, scope, field, value, reason);
        }
        if row.status.len() <= MAX_TRACEABILITY_FIELD_BYTES
            && !row.status.trim().is_empty()
            && !VALID_STATUSES.contains(&row.status)
        {
            diagnostics.push(diagnostic(
                scope,
                TraceabilityField::Status,
                format!(
                    "unsupported declaration status {:?}; this unbound registry cannot mint scientific proof status from tracker state",
                    row.status
                ),
            ));
        }

        if row.requirement_id.len() <= MAX_TRACEABILITY_FIELD_BYTES
            && !row.requirement_id.trim().is_empty()
            && !requirement_ids.insert(row.requirement_id)
        {
            diagnostics.push(diagnostic(
                scope,
                TraceabilityField::RequirementId,
                format!("duplicate requirement id {:?}", row.requirement_id),
            ));
        }
        if row.requirement_id.len() <= MAX_TRACEABILITY_FIELD_BYTES
            && !row.requirement_id.trim().is_empty()
            && !REQUIRED_REQUIREMENT_IDS.contains(&row.requirement_id)
        {
            diagnostics.push(diagnostic(
                scope,
                TraceabilityField::RequirementId,
                format!(
                    "requirement id {:?} is outside the closed B1-B14/RQ-* registry",
                    row.requirement_id
                ),
            ));
        }
        if row.proof_obligations.is_empty() {
            diagnostics.push(diagnostic(
                scope,
                TraceabilityField::ProofObligation,
                "requirement has no proof-obligation gate",
            ));
        } else if row.proof_obligations.len() > MAX_REQUIREMENT_PO_LINKS {
            diagnostics.push(diagnostic(
                scope,
                TraceabilityField::ProofObligation,
                format!(
                    "requirement links {} proof obligations; maximum is {MAX_REQUIREMENT_PO_LINKS}",
                    row.proof_obligations.len()
                ),
            ));
        } else {
            let mut row_obligations = BTreeSet::new();
            for id in row.proof_obligations {
                if id.len() > MAX_TRACEABILITY_FIELD_BYTES {
                    diagnostics.push(diagnostic(
                        scope,
                        TraceabilityField::ProofObligation,
                        format!(
                            "proof-obligation link is {} bytes; maximum is {MAX_TRACEABILITY_FIELD_BYTES}",
                            id.len()
                        ),
                    ));
                } else if id.trim().is_empty() {
                    diagnostics.push(diagnostic(
                        scope,
                        TraceabilityField::ProofObligation,
                        "requirement contains an empty proof-obligation link",
                    ));
                } else if !known_obligations.contains(id) {
                    diagnostics.push(diagnostic(
                        scope,
                        TraceabilityField::ProofObligation,
                        format!("requirement links unknown proof obligation {id:?}"),
                    ));
                } else if !row_obligations.insert(*id) {
                    diagnostics.push(diagnostic(
                        scope,
                        TraceabilityField::ProofObligation,
                        format!("requirement links proof obligation {id:?} more than once"),
                    ));
                }
            }
        }
    }

    for required_id in REQUIRED_REQUIREMENT_IDS {
        if !requirement_ids.contains(required_id) {
            diagnostics.push(diagnostic(
                required_id,
                TraceabilityField::RequirementId,
                format!("closed registry is missing requirement {required_id}"),
            ));
        }
    }

    diagnostics.sort_by(|left, right| {
        left.requirement_id
            .cmp(&right.requirement_id)
            .then_with(|| left.field.cmp(&right.field))
            .then_with(|| left.reason.cmp(&right.reason))
    });
    let incomplete_scopes: BTreeSet<&str> = diagnostics
        .iter()
        .map(|diagnostic| diagnostic.requirement_id.as_str())
        .collect();
    let complete = rows
        .iter()
        .filter(|row| !incomplete_scopes.contains(requirement_scope(row.requirement_id)))
        .count();

    TraceabilityAudit {
        total: rows.len(),
        complete,
        diagnostics,
    }
}

/// Validate and render a deterministic JSON ledger.
///
/// Input row order and PO-index order do not affect the artifact. Generation
/// returns the complete structured audit on any failure and never emits a
/// plausible-looking partial dashboard.
pub fn generate_traceability_ledger(
    rows: &[RequirementRow<'_>],
    obligations: &[ProofObligation<'_>],
) -> Result<String, TraceabilityAudit> {
    let audit = audit_traceability(rows, obligations);
    if !audit.ok() {
        return Err(audit);
    }

    let mut ordered_rows: Vec<&RequirementRow<'_>> = rows.iter().collect();
    ordered_rows.sort_unstable_by_key(|row| {
        REQUIRED_REQUIREMENT_IDS
            .iter()
            .position(|required| *required == row.requirement_id)
            .unwrap_or(usize::MAX)
    });
    let mut ordered_obligations: Vec<&ProofObligation<'_>> = obligations.iter().collect();
    ordered_obligations
        .sort_unstable_by_key(|obligation| valid_po_id(obligation.id).unwrap_or(u8::MAX));

    let mut out = String::new();
    write!(
        out,
        "{{\"schema\":\"{}\",\"authority\":\"{}\",\"source_snapshot\":null,\"requirements\":[",
        TRACEABILITY_SCHEMA, TRACEABILITY_AUTHORITY,
    )
    .expect("writing to a String is infallible");
    for (index, row) in ordered_rows.iter().enumerate() {
        if index > 0 {
            out.push(',');
        }
        let mut links = row.proof_obligations.to_vec();
        links.sort_unstable_by_key(|id| valid_po_id(id).unwrap_or(u8::MAX));
        write!(
            out,
            concat!(
                "{{\"requirement_id\":\"{}\",",
                "\"capability_property\":\"{}\",",
                "\"blocker\":\"{}\",",
                "\"owner_artifact\":\"{}\",",
                "\"prerequisite_phase\":\"{}\",",
                "\"milestone\":\"{}\",",
                "\"flagship\":\"{}\",",
                "\"benchmark_data\":\"{}\",",
                "\"proof_obligations\":["
            ),
            json_escape(row.requirement_id),
            json_escape(row.capability_property),
            json_escape(row.blocker),
            json_escape(row.owner_artifact),
            json_escape(row.prerequisite_phase),
            json_escape(row.milestone),
            json_escape(row.flagship),
            json_escape(row.benchmark_data),
        )
        .expect("writing to a String is infallible");
        for (link_index, id) in links.iter().enumerate() {
            if link_index > 0 {
                out.push(',');
            }
            write!(out, "\"{}\"", json_escape(id)).expect("writing to a String is infallible");
        }
        write!(
            out,
            "],\"claim_boundary\":\"{}\",\"status\":\"{}\"}}",
            json_escape(row.claim_boundary),
            json_escape(row.status),
        )
        .expect("writing to a String is infallible");
    }
    out.push_str("],\"proof_obligation_index\":[");
    for (index, obligation) in ordered_obligations.iter().enumerate() {
        if index > 0 {
            out.push(',');
        }
        let mut owners = obligation.owner_beads.to_vec();
        owners.sort_unstable();
        write!(
            out,
            "{{\"id\":\"{}\",\"summary\":\"{}\",\"owner_beads\":[",
            json_escape(obligation.id),
            json_escape(obligation.summary),
        )
        .expect("writing to a String is infallible");
        for (owner_index, owner) in owners.iter().enumerate() {
            if owner_index > 0 {
                out.push(',');
            }
            write!(out, "\"{}\"", json_escape(owner)).expect("writing to a String is infallible");
        }
        out.push_str("]}");
    }
    out.push_str("]}");
    Ok(out)
}

const PROOF_OBLIGATIONS: [ProofObligation<'static>; PROOF_OBLIGATION_COUNT] = [
    ProofObligation {
        id: "PO-1",
        summary: "incidence remains material-independent while constitutive blocks prove units, objectivity, stability, dissipation, reversible coupling, tangent consistency, quadrature, and coefficient robustness",
        owner_beads: &[
            "frankensim-ext-feec-weighted-operators-vxth",
            "frankensim-ext-constitutive-graph-kagp",
        ],
    },
    ProofObligation {
        id: "PO-2",
        summary: "junctions have zero power defect and storage, dissipation, sources, streams, losses, and thermal credits close in one accounting chart",
        owner_beads: &[
            "frankensim-ext-couple-port-schema-3feh",
            "frankensim-ext-thermal-domain-je8y",
        ],
    },
    ProofObligation {
        id: "PO-3",
        summary: "variational multibody claims remain conditional and distinguish residuals, state distance, projection defects, virtual power, and friction feasibility",
        owner_beads: &["frankensim-ext-mbd-core-vqqt"],
    },
    ProofObligation {
        id: "PO-4",
        summary: "complete event coverage proves compact-domain, regular-mode, true-flow, isolable-guard, reset-closure, and Zeno obligations; scans only falsify",
        owner_beads: &[
            "frankensim-ext-events-validated-prescribed-6b8h",
            "frankensim-ext-time-validated-step-ow2o",
        ],
    },
    ProofObligation {
        id: "PO-5",
        summary: "swept and envelope enclosures are sound and branch-complete, with conjugacy additionally proving contact and non-interference conditions",
        owner_beads: &["frankensim-ext-motion-swept-envelope-c58q"],
    },
    ProofObligation {
        id: "PO-6",
        summary: "reactive moving-mesh updates preserve EOS-specific entropy and admissibility while closing GCL, sources, boundaries, and nonlinear-solve defect",
        owner_beads: &["frankensim-ext-gas-entropy-lane-60np"],
    },
    ProofObligation {
        id: "PO-7",
        summary: "evidence closure contains every load-bearing receipt, excludes irrelevant ingredients, and keeps verification and validation orthogonal",
        owner_beads: &[
            "frankensim-ext-matdb-core-5hmy",
            "frankensim-ext-vv-artifact-schemas-x68z",
        ],
    },
    ProofObligation {
        id: "PO-8",
        summary: "weighted adjoint identities declare their metric and complex convention, regular discrete sensitivities are independently checked, and nonsmooth exceptions refuse",
        owner_beads: &[
            "frankensim-ext-adjoint-composition-easb",
            "frankensim-ext-sliding-interface-auh5",
        ],
    },
    ProofObligation {
        id: "PO-9",
        summary: "sliding and nonmatching interfaces prove compatibility, signed adjointness, preservation, stability, power balance, and moving-interface GCL",
        owner_beads: &["frankensim-ext-sliding-interface-auh5"],
    },
    ProofObligation {
        id: "PO-10",
        summary: "mobility, self-stress, and index reports agree with independent rank and duality checks or return Unknown; finite mobility needs higher-order evidence",
        owner_beads: &["frankensim-ext-kinematics-tangent-rigidity-e2iz"],
    },
    ProofObligation {
        id: "PO-11",
        summary: "reaction and diffusion laws conserve elements and charge, declare frame and nullspace, and produce nonnegative entropy on the independent-flux subspace",
        owner_beads: &[
            "frankensim-ext-thermochem-core-5fkv",
            "frankensim-ext-porous-capillary-biz5",
        ],
    },
    ProofObligation {
        id: "PO-12",
        summary: "circuit transitions preserve admitted state relations or solve distributional MNA with explicit energy-defect and initialization receipts, otherwise refusing",
        owner_beads: &["frankensim-ext-circuit-descriptor-mna-htfy"],
    },
    ProofObligation {
        id: "PO-13",
        summary: "electromagnetics passes gauge, source, winding, frame, force-convention, held-variable, and moving-conductor power-closure checks",
        owner_beads: &[
            "frankensim-ext-em-moving-conductor-jven",
            "frankensim-ext-em-forces-losses-o9im",
        ],
    },
    ProofObligation {
        id: "PO-14",
        summary: "moving meshes and fresh cells close GCL and remap balances plus equal-and-opposite fluid/body impulse, torque, work, and power",
        owner_beads: &["frankensim-ext-flux-lbm-moving-4v7w"],
    },
    ProofObligation {
        id: "PO-15",
        summary: "material and interface queries preserve units, frames, domains, definitions, covariance, provenance, and licenses without invalid permeability-family crosswalks",
        owner_beads: &[
            "frankensim-ext-matdb-core-5hmy",
            "frankensim-ext-porous-capillary-biz5",
        ],
    },
    ProofObligation {
        id: "PO-16",
        summary: "performance gates bind model, resolution, QoI error, machine fingerprint, and baseline; a naked throughput number has no claim authority",
        owner_beads: &["frankensim-ext-scale-qualification-0h2j"],
    },
    ProofObligation {
        id: "PO-17",
        summary: "theorem promotion requires machine-readable statements and assumptions, admitted nonvacuous instances, reproducible checkers, and an explicit TCB",
        owner_beads: &[
            "frankensim-ext-theorem-foundry-infra-zxob",
            "frankensim-ext-e8-summit-tek2",
        ],
    },
    ProofObligation {
        id: "PO-18",
        summary: "electromechanical and gluing cards prove gauge and cut invariance plus whole-interface work and GCL; force bounds enclose every named contribution",
        owner_beads: &["frankensim-ext-theorem-lane-electromech-u7gi"],
    },
    ProofObligation {
        id: "PO-19",
        summary: "whole-machine thermodynamics composes open-system first and second laws across multirate windows without hiding unresolved defects in components",
        owner_beads: &[
            "frankensim-ext-theorem-lane-thermo-tbss",
            "frankensim-ext-couple-cosim-lanes-pelj",
        ],
    },
    ProofObligation {
        id: "PO-20",
        summary: "fidelity descent and conjugate geometry preserve evidence monotonicity, bound naturality defects, and retain rather than erase admitted counterexamples",
        owner_beads: &[
            "frankensim-ext-theorem-lane-thermo-tbss",
            "frankensim-ext-theorem-lane-mechanics-ysu2",
        ],
    },
    ProofObligation {
        id: "PO-21",
        summary: "multicontact action/reaction, stored energy, dissipation, impact work, frictional heat, and wear updates close globally with model discrepancy exposed",
        owner_beads: &["frankensim-ext-contact-nonsmooth-lane-oh0i"],
    },
    ProofObligation {
        id: "PO-22",
        summary: "acoustic lanes prove applicable radiation, complex-power, dispersion, and interface reciprocity while keeping source validity separate from propagation",
        owner_beads: &["frankensim-ext-acoustics-core-0ja4"],
    },
    ProofObligation {
        id: "PO-23",
        summary: "sampled, periodic, and hybrid controls declare conventions and discharge stability, reachability, observability, implementation-error, timing, and fault obligations",
        owner_beads: &["frankensim-ext-control-runtime-0a1k"],
    },
    ProofObligation {
        id: "PO-24",
        summary: "life and reliability claims bind population, process, load spectrum, failure mode, dependence, sampling, and blind held-out validation",
        owner_beads: &[
            "frankensim-ext-solid-life-ladder-gahl",
            "frankensim-ext-uq-reliability-upgrades-pwui",
        ],
    },
    ProofObligation {
        id: "PO-25",
        summary: "stable identity, lineage, transactional snapshots, and safety cases reject ambiguous rebinding and half-committed state while keeping conformance orthogonal",
        owner_beads: &[
            "frankensim-ext-machine-ir-6iq3",
            "frankensim-ext-safety-emc-assurance-te0w",
        ],
    },
];

/// Complete executable-obligation index in canonical numeric order.
#[must_use]
pub fn proof_obligations() -> &'static [ProofObligation<'static>] {
    &PROOF_OBLIGATIONS
}

const REQUIREMENTS: [RequirementRow<'static>; REQUIREMENT_COUNT] = [
    RequirementRow {
        requirement_id: "B1",
        capability_property: "moving geometry, certified sweeps, envelopes, and ALE",
        blocker: "motion is not yet bound to geometry and only the ball case of Minkowski expansion exists",
        owner_artifact: "fs-motion::CertifiedMotorTube; Rep Router; fs-flux",
        prerequisite_phase: "E0a-E0d typed identity, operators, evidence, and replay spine",
        milestone: "E1 prescribed motion and E5 moving reactive flow",
        flagship: "prescribed-motion mechanisms and the reactive-flow machine stack",
        benchmark_data: "trochoid oracle plus moving-wall and geometric-conservation-law decks",
        proof_obligations: &["PO-5", "PO-14"],
        claim_boundary: "rigid motion ships first; deformable motion remains a separately gated successor",
        status: "proposed",
    },
    RequirementRow {
        requirement_id: "B2",
        capability_property: "multibody joints, DAEs, and nonholonomic dynamics",
        blocker: "there is no production joint graph, constrained DAE, nonholonomic solver, or complete event path",
        owner_artifact: "fs-kinematics; fs-mbd; fs-time",
        prerequisite_phase: "E0a-E0d typed machine graph and replay spine",
        milestone: "E1 nonlinear kinematics through E2 dynamic interaction",
        flagship: "fs-geneva-e2e and constant-width rolling mechanisms",
        benchmark_data: "IFToMM mechanisms, rolling disk, and Chaplygin sleigh decks",
        proof_obligations: &["PO-3", "PO-4", "PO-10"],
        claim_boundary: "regular and inclusion-valued event classes remain explicit and are never silently conflated",
        status: "proposed",
    },
    RequirementRow {
        requirement_id: "B3",
        capability_property: "body contact, continuous collision detection, and penetration handling",
        blocker: "deformable-to-fixed-obstacle contact exists but general body-to-body CCD and penetration do not",
        owner_artifact: "capability-routed fs-query; fs-contact; fs-tribo",
        prerequisite_phase: "E1 certified motion and geometry queries",
        milestone: "E2 dynamic interaction",
        flagship: "fs-geneva-e2e and multicontact machine-element slices",
        benchmark_data: "Hertz, Painleve, and adversarial multicontact decks",
        proof_obligations: &["PO-21"],
        claim_boundary: "nonintersection is conditional on the admitted motion and contact capability class",
        status: "proposed",
    },
    RequirementRow {
        requirement_id: "B4",
        capability_property: "gears, cams, Geneva drives, and general linkage mechanisms",
        blocker: "the production architecture has no complete mechanism vocabulary or solver surface",
        owner_artifact: "fs-kinematics; fs-machine",
        prerequisite_phase: "E0a-E0d machine identity plus E1 certified motion",
        milestone: "E1 prescribed kinematics through E3 machine elements",
        flagship: "fs-geneva-e2e and fs-gear-e2e",
        benchmark_data: "law-of-gearing, conjugacy, Geneva, loaded transmission-error, and continuation decks",
        proof_obligations: &["PO-5", "PO-10", "PO-20"],
        claim_boundary: "each mechanism family retains its own contact and standards scope",
        status: "proposed",
    },
    RequirementRow {
        requirement_id: "B5",
        capability_property: "electromagnetic formulations and force coupling",
        blocker: "there is no curl-curl assembly, gauge path, EM material field, complex Krylov lane, or audited EM force surface",
        owner_artifact: "fs-rep-mesh; fs-feec; fs-em",
        prerequisite_phase: "E0a typed identity and E0b weighted operators",
        milestone: "E0 operator spine through E4 electromechanical drive",
        flagship: "fs-motor-e2e and generator validation slices",
        benchmark_data: "exactness, orientation, manufactured solutions, TEAM problems, and motor/generator decks",
        proof_obligations: &["PO-1", "PO-8", "PO-9", "PO-13"],
        claim_boundary: "magnetoquasistatics is the first admitted rung; full-wave behavior is separate",
        status: "proposed",
    },
    RequirementRow {
        requirement_id: "B6",
        capability_property: "typed transient thermal transport, radiation, contact, and phase change",
        blocker: "generic scalar diffusion fixtures do not constitute a thermal-domain solver stack",
        owner_artifact: "fs-thermal; fs-xform; fs-couple",
        prerequisite_phase: "E0a typed quantities and E0b coupled operators",
        milestone: "E0 thermal law spine through E5 machine thermal closure",
        flagship: "fs-motor-e2e and fs-wankel-e2e thermal ladders",
        benchmark_data: "Stefan, NAFEMS, laser-flash, and calorimetry decks",
        proof_obligations: &["PO-1", "PO-2"],
        claim_boundary: "data, model, radiation, contact, and phase-change rungs retain separate evidence colors",
        status: "proposed",
    },
    RequirementRow {
        requirement_id: "B7",
        capability_property: "compressible, species, thermochemical, and combustion transport",
        blocker: "no production EOS, species, reactive-gas transport, or combustion API exists",
        owner_artifact: "fs-thermochem; fs-gas; fs-flux",
        prerequisite_phase: "E0a typed amount units plus E0b operators and E0c evidence",
        milestone: "E5 reactive-flow machine stack through E6 high fidelity",
        flagship: "fs-wankel-e2e and fs-genset-e2e",
        benchmark_data: "Riemann, nozzle, Engine Combustion Network, and engine decks",
        proof_obligations: &["PO-6", "PO-11"],
        claim_boundary: "closure validity and experimental validity are named separately for every QoI",
        status: "proposed",
    },
    RequirementRow {
        requirement_id: "B8",
        capability_property: "cross-physics quantities, materials, evidence, and interface-history discipline",
        blocker: "materials are mechanics-first and key property paths still bypass typed cross-domain queries",
        owner_artifact: "fs-qty; fs-evidence; fs-matdb; fs-material",
        prerequisite_phase: "constitutional typed-units and provenance substrate",
        milestone: "E0a units and identity plus E0c data and V&V registry",
        flagship: "all coupled machine vertical slices",
        benchmark_data: "offline-pack, typed-query, covariance, license, and cross-version batteries",
        proof_obligations: &["PO-7", "PO-11", "PO-15"],
        claim_boundary: "immutable scientific evidence remains distinct from mutable runtime state",
        status: "proposed",
    },
    RequirementRow {
        requirement_id: "B9",
        capability_property: "typed multiphysics ports, transfer operators, and multirate coupling",
        blocker: "the coupling surface is a scalar fixture without field transfer, mortar, vector acceleration, or multirate semantics",
        owner_artifact: "fs-couple::PortSchema; domain transfer adapters",
        prerequisite_phase: "E0a typed identity and E0b neutral ports",
        milestone: "E0 coupling spine through E7 whole-machine synthesis",
        flagship: "motor, Wankel, genset, pump-bearing, and active-material vertical slices",
        benchmark_data: "power-pairing, transfer-adjoint, interface-GCL, and coupled vertical-slice batteries",
        proof_obligations: &["PO-2", "PO-8", "PO-9", "PO-19"],
        claim_boundary: "graph topology alone never proves passivity or whole-machine thermodynamics",
        status: "proposed",
    },
    RequirementRow {
        requirement_id: "B10",
        capability_property: "sparse nonlinear solvers, validated time integration, and composed adjoints",
        blocker: "implicit integrators are dense-only and production adjoints cover only narrow seed paths",
        owner_artifact: "fs-solver; fs-time; fs-ad; fs-adjoint",
        prerequisite_phase: "E0a execution identity and E0b operator protocols",
        milestone: "E0 solver spine through E6 high-fidelity escalation",
        flagship: "event-driven, reactive-flow, and optimization vertical slices",
        benchmark_data: "manufactured solutions, event, cancellation, restart, and independent gradient gates",
        proof_obligations: &["PO-3", "PO-4", "PO-6", "PO-8"],
        claim_boundary: "nonsmooth and irregular paths carry explicit no-gradient or generalized-sensitivity boundaries",
        status: "proposed",
    },
    RequirementRow {
        requirement_id: "B11",
        capability_property: "sparse mechanism spectral analysis and continuation",
        blocker: "the generic spectral crate does not yet own the sparse kinematic and sheaf complex required by machine mechanisms",
        owner_artifact: "fs-spectral; fs-kinematics; fs-time",
        prerequisite_phase: "E0b sparse solver service and E1 mechanism complex",
        milestone: "E1 nonlinear kinematics through E3 rotordynamics",
        flagship: "mechanism continuation, gear, and rotating-machine slices",
        benchmark_data: "Bennett, Bricard, eigen-gap, continuation, and Campbell-diagram decks",
        proof_obligations: &["PO-10"],
        claim_boundary: "scaling, gauge, multiplicity, and unresolved spectral gaps are serialized rather than guessed",
        status: "proposed",
    },
    RequirementRow {
        requirement_id: "B12",
        capability_property: "oriented two-dimensional complexes and FEEC families for planar fields",
        blocker: "fs-feec is three-dimensional simplicial and the existing surface mesh is not a production 2-D FEEC complex",
        owner_artifact: "fs-rep-mesh::TriComplex2; fs-feec 2-D families",
        prerequisite_phase: "E0a stable topology identity and E0b weighted operators",
        milestone: "E0 operator spine through E4 planar electromagnetics",
        flagship: "planar EM and fs-motor-e2e",
        benchmark_data: "incidence exactness, orientation, commuting projection, manufactured solutions, and planar EM decks",
        proof_obligations: &["PO-1", "PO-13"],
        claim_boundary: "a surface half-edge mesh alone does not satisfy the complex or FEEC claim",
        status: "proposed",
    },
    RequirementRow {
        requirement_id: "B13",
        capability_property: "machine flexibility, life, NVH, power electronics, and control",
        blocker: "the production stack lacks flexible machines, rotor-bearing-seal dynamics, degradation, acoustics, and closed-loop validation",
        owner_artifact: "fs-mbd; fs-solid; fs-acoustics; fs-control; fs-uq upgrades",
        prerequisite_phase: "E0-E2 typed machine, operators, motion, and interaction",
        milestone: "E3 machine elements through E7 whole-machine synthesis",
        flagship: "fs-gear-e2e, fs-motor-e2e, fs-wankel-e2e, and fs-genset-e2e",
        benchmark_data: "Campbell, NVH, fault, coupon-to-component, and held-out population/QoI decks",
        proof_obligations: &["PO-22", "PO-23", "PO-24"],
        claim_boundary: "population, process, failure mode, source validity, and control regime remain explicit",
        status: "proposed",
    },
    RequirementRow {
        requirement_id: "B14",
        capability_property: "machine V&V, safety, EMC, standards, workflow, and scale qualification",
        blocker: "no operational assurance case, conformance axis, engineer-interchange workflow, or competitive qualification program exists",
        owner_artifact: "Machine IR; V&V artifacts; safety case; ScaleQualification and CompetitiveCapability ledgers",
        prerequisite_phase: "E0a-E0d identity, schema, evidence, and replay spine",
        milestone: "E0 assurance spine through E8 theorem-foundry summit",
        flagship: "every whole-machine and theorem-foundry flagship",
        benchmark_data: "blind tests, injected faults, FMI/SSP quarantine, scale suites, and independent reproduction packs",
        proof_obligations: &["PO-16", "PO-17", "PO-25"],
        claim_boundary: "scientific receipts, standards conformance, regulatory approval, and certification remain orthogonal",
        status: "proposed",
    },
    RequirementRow {
        requirement_id: "RQ-ROLL",
        capability_property: "truly rolling constant-width and Reuleaux mechanisms",
        blocker: "support and centroid motion, nonholonomic rolling, branch switches, and tolerance effects are not integrated",
        owner_artifact: "nonholonomic fs-mbd; fs-contact; fs-constant-width-e2e",
        prerequisite_phase: "E1 certified geometry and E2 nonholonomic contact",
        milestone: "E2 dynamic interaction and E7 whole-machine synthesis",
        flagship: "fs-constant-width-e2e",
        benchmark_data: "support-function, centroid, branch-switch, rolling-without-slip, and tolerance decks",
        proof_obligations: &["PO-3", "PO-4", "PO-21"],
        claim_boundary: "rolling contact is distinct from sliding Wankel-seal contact",
        status: "proposed",
    },
    RequirementRow {
        requirement_id: "RQ-GEAR",
        capability_property: "broad real gear families and loaded transmission behavior",
        blocker: "conjugate, spatial, compliant, and standards-scoped gear mechanisms are absent",
        owner_artifact: "conjugate and spatial fs-kinematics; fs-machine",
        prerequisite_phase: "E1 certified conjugacy and E2 contact",
        milestone: "E1 prescribed kinematics through E3 machine elements",
        flagship: "fs-gear-e2e",
        benchmark_data: "law-of-gearing, conjugacy, loaded transmission error, contact, and exact ISO-scope decks",
        proof_obligations: &["PO-5", "PO-10", "PO-20"],
        claim_boundary: "each family retains its own contact class, validity domain, and standards edition",
        status: "proposed",
    },
    RequirementRow {
        requirement_id: "RQ-FRICTION",
        capability_property: "history-dependent friction, contact, lubrication, and wear",
        blocker: "friction is not yet represented as a stateful interface-system law with thermal and wear history",
        owner_artifact: "InterfaceSystemCard; fs-tribo runtime state",
        prerequisite_phase: "E1 material/interface identity and E2 contact",
        milestone: "E2 dry interaction through E3 lubrication and wear",
        flagship: "fs-geneva-e2e and fs-gear-e2e contact systems",
        benchmark_data: "incline, Hertz, EHL film and traction, flash-temperature, and wear holdouts",
        proof_obligations: &["PO-15", "PO-21"],
        claim_boundary: "friction is system-, state-, process-, and history-specific rather than a universal pair constant",
        status: "proposed",
    },
    RequirementRow {
        requirement_id: "RQ-CONSTITUTIVE",
        capability_property: "weighted, nonlinear, coupled, and history-dependent matter laws",
        blocker: "the common constitutive graph and its coupled thermodynamic admission checks are not complete",
        owner_artifact: "fs-material::ConstitutiveGraph; domain constitutive nodes",
        prerequisite_phase: "E0a typed quantities and E0b weighted operators",
        milestone: "E0 constitutive spine through E8 theorem integration",
        flagship: "all coupled domain and whole-machine flagships",
        benchmark_data: "coupled-law manufactured solutions, tangent checks, entropy audits, and held-out material tests",
        proof_obligations: &["PO-1", "PO-7"],
        claim_boundary: "reversible skew coupling, dissipation, stability, and empirical validity remain separate claims",
        status: "proposed",
    },
    RequirementRow {
        requirement_id: "RQ-DENSITY",
        capability_property: "density to mass, center-of-mass, and inertia propagation",
        blocker: "mass-property evidence is not yet joined to material provenance and cross-representation geometry",
        owner_artifact: "fs-matdb property receipt; fs-query::GeometricMoments",
        prerequisite_phase: "E0a typed quantities and E0c material evidence",
        milestone: "E1 prescribed motion and geometric moments",
        flagship: "all rotating and translating machine flagships",
        benchmark_data: "analytic solids, transformed bodies, and cross-representation mass-property fixtures",
        proof_obligations: &["PO-7", "PO-15"],
        claim_boundary: "nonuniform density escalates to an admitted integration and uncertainty path",
        status: "proposed",
    },
    RequirementRow {
        requirement_id: "RQ-MECHMAT",
        capability_property: "ductility, hardness, toughness, fatigue, creep, and corrosion",
        blocker: "mechanical material cards do not yet close the coupon-to-component-to-machine life ladder",
        owner_artifact: "fs-material; fs-solid life ladder",
        prerequisite_phase: "E0c material/V&V registry and E2 load histories",
        milestone: "E3 machine life through E7 reliability",
        flagship: "gear, bearing, motor, Wankel, and genset life cases",
        benchmark_data: "coupon, crack-growth, creep, corrosion, component, and machine holdouts",
        proof_obligations: &["PO-15", "PO-24"],
        claim_boundary: "no coupon result or generic property card grants a universal service-life claim",
        status: "proposed",
    },
    RequirementRow {
        requirement_id: "RQ-ELEC",
        capability_property: "electrical and thermal conductivity plus dielectric and insulation response",
        blocker: "typed property cards are not yet consumed by coupled electrical, EM, thermal, and insulation operators",
        owner_artifact: "weighted operators; fs-em; fs-thermal; fs-power",
        prerequisite_phase: "E0a typed properties and E0b coupled operators",
        milestone: "E0 operator spine through E4 electromechanical drive",
        flagship: "fs-motor-e2e and power-conversion slices",
        benchmark_data: "resistance, Joule-heat closure, dielectric, insulation, and thermal-conduction fixtures",
        proof_obligations: &["PO-1", "PO-2", "PO-13"],
        claim_boundary: "frequency, temperature, field, process, and insulation regime are named for every claim",
        status: "proposed",
    },
    RequirementRow {
        requirement_id: "RQ-MAG",
        capability_property: "permeability, BH response, remanence, hysteresis, and body-integrated magnetic moment",
        blocker: "magnetic state laws and force/moment conventions are not integrated into the EM evidence stack",
        owner_artifact: "energy and state cards; fs-em; body-integral moment artifact",
        prerequisite_phase: "E0 constitutive spine and E4 planar EM",
        milestone: "E4 electromechanical drive stack",
        flagship: "fs-motor-e2e and generator slices",
        benchmark_data: "TEAM problems, BH and loss fixtures, force checks, and motor dynamometer data",
        proof_obligations: &["PO-1", "PO-13"],
        claim_boundary: "there is no geometry-free total moment or universal scalar-spline magnetic theorem",
        status: "proposed",
    },
    RequirementRow {
        requirement_id: "RQ-PHASE",
        capability_property: "latent heat, phase-dependent density, and phase-front remap",
        blocker: "total-enthalpy and phase-state updates are not yet joined to thermal transport and conservative remap",
        owner_artifact: "fs-thermal total enthalpy; phase state; fs-xform remap",
        prerequisite_phase: "E0 thermal laws and typed state identity",
        milestone: "E0 thermal spine through E5 reactive-flow machines",
        flagship: "fs-wankel-e2e thermal and phase-change slices",
        benchmark_data: "Stefan, laser-flash, phase-density, and remap-conservation decks",
        proof_obligations: &["PO-1", "PO-2"],
        claim_boundary: "apparent-heat regularization is labeled and never presented as exact sharp-interface dynamics",
        status: "proposed",
    },
    RequirementRow {
        requirement_id: "RQ-FLUID",
        capability_property: "viscosity, bulk modulus, vapor pressure, hydraulics, and EHL",
        blocker: "fluid property definitions, phase regimes, hydraulic transients, and tribological consumption are not unified",
        owner_artifact: "fs-matdb; fs-flux; fs-gas; fs-tribo",
        prerequisite_phase: "E0 typed properties and E2 contact/interface state",
        milestone: "E3 machine lubrication through E5 reactive flow",
        flagship: "pump-bearing, gear, Wankel, and genset slices",
        benchmark_data: "Couette, bearing, cavitation, water-hammer, EHL, and leakage decks",
        proof_obligations: &["PO-2", "PO-14", "PO-15", "PO-21"],
        claim_boundary: "rheology, compressibility, cavitation, thermal, and phase regimes remain explicit",
        status: "proposed",
    },
    RequirementRow {
        requirement_id: "RQ-PERMEATE",
        capability_property: "intrinsic permeability, gas transport, porous media, and membranes",
        blocker: "porous and membrane transport adapters do not yet preserve their distinct definitions and evidence",
        owner_artifact: "transport adapters from charter section 5.13; fs-porous-capillary",
        prerequisite_phase: "E0 typed quantities and E5 thermochemical transport",
        milestone: "E5 reactive flow through E6 high fidelity",
        flagship: "Wankel leakage and turbo-e-fuel transport slices",
        benchmark_data: "Darcy, layered-media, membrane, diffusion-frame, and definition-crosswalk decks",
        proof_obligations: &["PO-11", "PO-15"],
        claim_boundary: "permeability, hydraulic conductivity, and permeance are distinct typed properties",
        status: "proposed",
    },
    RequirementRow {
        requirement_id: "RQ-WET",
        capability_property: "wetting, dynamic contact angles, and capillarity",
        blocker: "wetting is not yet represented as a dynamic interface system coupled to free surfaces and porous transport",
        owner_artifact: "dynamic InterfaceSystemCard; free-surface and porous-capillary lanes",
        prerequisite_phase: "E2 interface state and E5 multiphase transport",
        milestone: "E5 reactive flow through E6 high fidelity",
        flagship: "seal-leakage, spray, porous, and thermal-management slices",
        benchmark_data: "Young-Laplace, capillary rise, and advancing/receding contact-angle tests",
        proof_obligations: &["PO-2", "PO-15"],
        claim_boundary: "rate, roughness, hysteresis, contamination, and validity domain are always bound",
        status: "proposed",
    },
    RequirementRow {
        requirement_id: "RQ-MOTORGEN",
        capability_property: "motors and generators with electronics, control, and thermal behavior",
        blocker: "EM, circuit, converter, control, thermal, fault, and assurance artifacts are not yet composed",
        owner_artifact: "fs-motor-e2e",
        prerequisite_phase: "E0-E3 typed machine, operators, interaction, and machine elements",
        milestone: "E4 electromechanical drive and E7 whole-machine synthesis",
        flagship: "fs-motor-e2e",
        benchmark_data: "TEAM, converter, closed-loop fault, thermal-derating, dynamometer, and EMC decks",
        proof_obligations: &["PO-2", "PO-8", "PO-9", "PO-13", "PO-23", "PO-25"],
        claim_boundary: "fidelity rung, held electrical variables, control regime, and validation population are explicit",
        status: "proposed",
    },
    RequirementRow {
        requirement_id: "RQ-ICE",
        capability_property: "piston and Wankel internal-combustion engines plus gensets",
        blocker: "moving geometry, contact, thermochemistry, reactive flow, heat transfer, emissions, and reliability are not integrated",
        owner_artifact: "fs-wankel-e2e; fs-genset-e2e",
        prerequisite_phase: "E0-E4 typed machine through electromechanical drive",
        milestone: "E5 reactive-flow machines through E7 whole-machine synthesis",
        flagship: "fs-wankel-e2e and fs-genset-e2e",
        benchmark_data: "pressure, Engine Combustion Network, calorimetry, emissions, seal, blow-by, and held-out engine decks",
        proof_obligations: &["PO-2", "PO-4", "PO-6", "PO-11", "PO-14", "PO-24"],
        claim_boundary: "correlation and validation colors remain configuration- and QoI-specific",
        status: "proposed",
    },
    RequirementRow {
        requirement_id: "RQ-ACOUSTIC",
        capability_property: "NVH, aeroacoustics, combustion acoustics, and structural-acoustic coupling",
        blocker: "there is no production acoustic formulation ladder or validated source-transfer path",
        owner_artifact: "fs-acoustics",
        prerequisite_phase: "E2 moving interfaces and E3 flexible machine elements",
        milestone: "E3 machine NVH through E7 whole-machine validation",
        flagship: "gear, motor, Wankel, and genset NVH slices",
        benchmark_data: "NAFEMS R0083, NASA CAA, analytic radiation, dispersion, and measured acoustic datasets",
        proof_obligations: &["PO-22"],
        claim_boundary: "source-model validity is separate from propagation and radiation correctness",
        status: "proposed",
    },
    RequirementRow {
        requirement_id: "RQ-ACTIVE",
        capability_property: "piezoelectric, magnetostrictive, thermoelectric, and electrochemical coupling",
        blocker: "active-material cross-couplings and their thermodynamic admission checks are not yet implemented",
        owner_artifact: "charter section 5.13 coupled-law adapters; fs-material ConstitutiveGraph nodes",
        prerequisite_phase: "E0 coupled constitutive spine and relevant E4-E6 domain solvers",
        milestone: "E4 electromechanics through E8 theorem integration",
        flagship: "active-material sensor, actuator, harvester, and electrochemical expansion slices",
        benchmark_data: "direct/converse reciprocity, Kelvin relation, energy/entropy closure, parity, and held-out material decks",
        proof_obligations: &["PO-1", "PO-2"],
        claim_boundary: "frontier and moonshot paths remain guarded until independent validation and stability evidence exist",
        status: "proposed",
    },
];

/// Canonical B1-B14 plus RQ-* requirement registry in charter order.
#[must_use]
pub fn requirements() -> &'static [RequirementRow<'static>] {
    &REQUIREMENTS
}

/// Generate the canonical machine-readable ledger from the code registry.
///
/// The result remains fallible so a future incomplete edit cannot silently
/// publish a partial catalog.
pub fn traceability_ledger_json() -> Result<String, TraceabilityAudit> {
    generate_traceability_ledger(requirements(), proof_obligations())
}
