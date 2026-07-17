//! Semantic-presence join for loaded traceability sources.
//!
//! The filesystem adapter binds exact source bytes and extracts only canonical
//! top-level Beads ids. This module performs the next deliberately narrow
//! operation: every proof-obligation owner named by the declaration registry
//! must be present in that bound id index. Presence is not closure, status,
//! assignment, dependency satisfaction, scientific evidence, or proof.

use crate::traceability::{
    BoundTraceabilityLedger, ProofObligation, RequirementRow, TraceabilityAudit,
    audit_traceability, generate_traceability_ledger_from_snapshot,
};
use crate::traceability_fs::LoadedTraceabilitySourceSnapshot;
use fs_blake3::ContentHash;
use std::collections::BTreeSet;

/// Versioned meaning of the lexical Beads-owner presence join.
pub const TRACEABILITY_OWNER_JOIN_VERSION: &str = "frankensim-traceability-owner-presence-join-v1";

/// Field named by one owner-presence diagnostic.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum TraceabilityOwnerJoinField {
    /// Aggregate loaded-source join shape.
    Join,
    /// One proof-obligation owner Bead.
    OwnerBead,
}

impl TraceabilityOwnerJoinField {
    /// Stable diagnostic field name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Join => "join",
            Self::OwnerBead => "owner_bead",
        }
    }
}

/// One deterministic owner-presence refusal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TraceabilityOwnerJoinDiagnostic {
    /// Proof-obligation id, or `<join>` for an aggregate refusal.
    pub proof_obligation_id: String,
    /// Missing owner id, or `<none>` for an aggregate refusal.
    pub owner_bead: String,
    /// Failed join field.
    pub field: TraceabilityOwnerJoinField,
    /// Actionable refusal reason.
    pub reason: String,
}

/// Complete lexical owner-presence audit over one loaded snapshot.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TraceabilityOwnerJoinAudit {
    /// Number of owner references in the supplied obligation index.
    pub total_owner_references: usize,
    /// Number of distinct declared owner ids.
    pub distinct_owner_beads: usize,
    /// Number of bound Beads source files.
    pub beads_source_count: usize,
    /// Number of distinct indexed Beads records across those files.
    pub indexed_beads: usize,
    /// Every deterministic join refusal.
    pub diagnostics: Vec<TraceabilityOwnerJoinDiagnostic>,
}

impl TraceabilityOwnerJoinAudit {
    /// Whether a nonempty owner set is completely present in the loaded index.
    #[must_use]
    pub fn ok(&self) -> bool {
        self.total_owner_references > 0
            && self.distinct_owner_beads > 0
            && self.beads_source_count > 0
            && self.indexed_beads > 0
            && self.diagnostics.is_empty()
    }
}

/// Green receipt for the lexical owner-presence join.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TraceabilityOwnerJoinReceipt {
    version: &'static str,
    total_owner_references: usize,
    distinct_owner_beads: usize,
    beads_source_count: usize,
    indexed_beads: usize,
    source_snapshot_identity: ContentHash,
}

impl TraceabilityOwnerJoinReceipt {
    /// Owner-join semantics version.
    #[must_use]
    pub const fn version(&self) -> &'static str {
        self.version
    }

    /// Total owner references checked.
    #[must_use]
    pub const fn total_owner_references(&self) -> usize {
        self.total_owner_references
    }

    /// Distinct declared owner ids checked.
    #[must_use]
    pub const fn distinct_owner_beads(&self) -> usize {
        self.distinct_owner_beads
    }

    /// Bound Beads source files contributing ids.
    #[must_use]
    pub const fn beads_source_count(&self) -> usize {
        self.beads_source_count
    }

    /// Distinct indexed Beads records.
    #[must_use]
    pub const fn indexed_beads(&self) -> usize {
        self.indexed_beads
    }

    /// Exact source-snapshot identity whose id index was audited.
    #[must_use]
    pub const fn source_snapshot_identity(&self) -> ContentHash {
        self.source_snapshot_identity
    }
}

/// Bound declaration ledger paired with its lexical owner-presence receipt.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LoadedTraceabilityLedger {
    ledger: BoundTraceabilityLedger,
    owner_join: TraceabilityOwnerJoinReceipt,
}

impl LoadedTraceabilityLedger {
    /// Existing declaration-only ledger bound to the exact source snapshot.
    #[must_use]
    pub const fn ledger(&self) -> &BoundTraceabilityLedger {
        &self.ledger
    }

    /// Lexical owner-presence receipt for that same source snapshot.
    #[must_use]
    pub const fn owner_join(&self) -> &TraceabilityOwnerJoinReceipt {
        &self.owner_join
    }

    /// Consume into the bound declaration ledger and owner-presence receipt.
    #[must_use]
    pub fn into_parts(self) -> (BoundTraceabilityLedger, TraceabilityOwnerJoinReceipt) {
        (self.ledger, self.owner_join)
    }
}

/// Failure stage for loaded-source ledger generation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LoadedTraceabilityGenerationError {
    /// The declaration rows or proof-obligation index are structurally invalid.
    Declaration(TraceabilityAudit),
    /// At least one declared owner is absent from the bound Beads id index.
    OwnerJoin(TraceabilityOwnerJoinAudit),
}

fn diagnostic(
    proof_obligation_id: impl Into<String>,
    owner_bead: impl Into<String>,
    field: TraceabilityOwnerJoinField,
    reason: impl Into<String>,
) -> TraceabilityOwnerJoinDiagnostic {
    TraceabilityOwnerJoinDiagnostic {
        proof_obligation_id: proof_obligation_id.into(),
        owner_bead: owner_bead.into(),
        field,
        reason: reason.into(),
    }
}

/// Audit declared proof-obligation owners against exact loaded Beads bytes.
///
/// This checks canonical id presence only. Fields such as `status`,
/// `dependencies`, `assignee`, and `closed_at` are intentionally ignored.
#[must_use]
pub fn audit_traceability_owner_join(
    obligations: &[ProofObligation<'_>],
    loaded: &LoadedTraceabilitySourceSnapshot,
) -> TraceabilityOwnerJoinAudit {
    let mut diagnostics = Vec::new();
    let mut indexed_ids = BTreeSet::new();
    let mut beads_source_count = 0usize;
    for source in loaded.receipt().sources() {
        if source.kind() != crate::traceability::TraceabilitySourceKind::Beads {
            continue;
        }
        beads_source_count += 1;
        indexed_ids.extend(source.beads_ids().iter().map(String::as_str));
    }
    if beads_source_count == 0 || indexed_ids.is_empty() {
        diagnostics.push(diagnostic(
            "<join>",
            "<none>",
            TraceabilityOwnerJoinField::Join,
            "loaded source snapshot has no indexed Beads ids",
        ));
    }

    let mut total_owner_references = 0usize;
    let mut distinct_owners = BTreeSet::new();
    for obligation in obligations {
        for owner in obligation.owner_beads {
            total_owner_references = match total_owner_references.checked_add(1) {
                Some(total) => total,
                None => {
                    diagnostics.push(diagnostic(
                        "<join>",
                        "<none>",
                        TraceabilityOwnerJoinField::Join,
                        "owner reference count overflowed",
                    ));
                    continue;
                }
            };
            distinct_owners.insert(*owner);
            if !indexed_ids.contains(owner) {
                diagnostics.push(diagnostic(
                    obligation.id,
                    *owner,
                    TraceabilityOwnerJoinField::OwnerBead,
                    "declared proof-obligation owner is absent from the bound Beads id index",
                ));
            }
        }
    }
    if total_owner_references == 0 {
        diagnostics.push(diagnostic(
            "<join>",
            "<none>",
            TraceabilityOwnerJoinField::Join,
            "proof-obligation owner index is empty",
        ));
    }
    diagnostics.sort_by(|left, right| {
        left.proof_obligation_id
            .cmp(&right.proof_obligation_id)
            .then_with(|| left.owner_bead.cmp(&right.owner_bead))
            .then_with(|| left.field.cmp(&right.field))
            .then_with(|| left.reason.cmp(&right.reason))
    });
    diagnostics.dedup();
    TraceabilityOwnerJoinAudit {
        total_owner_references,
        distinct_owner_beads: distinct_owners.len(),
        beads_source_count,
        indexed_beads: indexed_ids.len(),
        diagnostics,
    }
}

/// Generate a declaration-only ledger from structurally valid declarations and
/// a concrete source snapshot whose Beads bytes contain every declared owner.
///
/// # Errors
/// Returns [`LoadedTraceabilityGenerationError::Declaration`] before joining
/// when the row/obligation registry is invalid, or
/// [`LoadedTraceabilityGenerationError::OwnerJoin`] when an owner id is absent.
/// No partial ledger is returned.
pub fn generate_traceability_ledger_from_loaded_sources(
    rows: &[RequirementRow<'_>],
    obligations: &[ProofObligation<'_>],
    loaded: &LoadedTraceabilitySourceSnapshot,
) -> Result<LoadedTraceabilityLedger, LoadedTraceabilityGenerationError> {
    let declaration_audit = audit_traceability(rows, obligations);
    if !declaration_audit.ok() {
        return Err(LoadedTraceabilityGenerationError::Declaration(
            declaration_audit,
        ));
    }
    let join_audit = audit_traceability_owner_join(obligations, loaded);
    if !join_audit.ok() {
        return Err(LoadedTraceabilityGenerationError::OwnerJoin(join_audit));
    }
    let ledger = generate_traceability_ledger_from_snapshot(rows, obligations, loaded.snapshot())
        .map_err(LoadedTraceabilityGenerationError::Declaration)?;
    let owner_join = TraceabilityOwnerJoinReceipt {
        version: TRACEABILITY_OWNER_JOIN_VERSION,
        total_owner_references: join_audit.total_owner_references,
        distinct_owner_beads: join_audit.distinct_owner_beads,
        beads_source_count: join_audit.beads_source_count,
        indexed_beads: join_audit.indexed_beads,
        source_snapshot_identity: loaded.snapshot().identity(),
    };
    Ok(LoadedTraceabilityLedger { ledger, owner_join })
}
