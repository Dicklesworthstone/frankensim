//! VALUE-OF-INFORMATION QUERIES (addendum Proposal C, bead knh1.6;
//! \[F\] — behind the `voi-queries` feature): THE IGNORANCE MARKET, v0
//! as a RANKED LIST. Across everything the ledger is uncertain about,
//! where does one dollar of evidence most change the downstream
//! decision? Decision sensitivity (does the decision FLIP inside the
//! node's interval?) is computed from CACHED surrogate sweeps, crossed
//! with a PRICED PROBE MENU that unifies computational and physical
//! experiments, ranked by SAMPLED FLIP-FRACTION REDUCTION PER DOLLAR.
//!
//! MYOPIC one-step VoI ONLY (the proposal's own discipline: full
//! sequential VoI is intractable and myopic captures most of the
//! value). The output surfaces as (i) the hint on query results — the
//! upgrade the fs-ir anytime module's CONTRACT reserved for Proposal C
//! — and (ii) the scheduler for discrepancy probes.
//!
//! THE KILL CRITERION AS CODE: [`VoiScheduler`] owns one append-only
//! chronological pairwise e-process, the remaining spend budget, and consumed
//! decision snapshots. Scheduling consults the CURRENT audit verdict under
//! exclusive mutation; stale reports carry no spending authority. The raw
//! path remains advisory; authenticated scheduling additionally reverifies
//! sealed ledger receipts for the ranked catalog, decision/model identities,
//! prospective outcomes, branch snapshot, policy versions, and expiry.
//!
//! Decision evaluation is cooperatively cancellation-aware: each call is
//! bracketed by [`Cx`] checkpoints and receives a private declared-work permit
//! under an explicit [`DecisionBudget`]. [`LiveDecision`] is only the bounded-
//! invocation adapter for an already-cached callback; arbitrary callback time
//! and memory remain outside the library's enforceable claim.

use std::collections::BTreeSet;

pub use asupersync::Cx;
use fs_blake3::{ContentHash, hash_bytes, hash_domain};
use fs_eproc::{LossSpan, PairwiseRace};
use fs_ledger::vcs::Vcs;
use fs_ledger::{EdgeRole, ExecMode, FiveExplicits, Ledger, OpArtifactEdge, OpOutcome};

/// Maximum uncertainty nodes in one myopic VoI request.
pub const MAX_VOI_NODES: usize = 256;
/// Maximum probe menu entries (and scheduled ranked entries).
pub const MAX_VOI_PROBES: usize = 1024;
/// Maximum visible-ASCII byte length of node, probe, target, and audit names.
pub const MAX_VOI_NAME_BYTES: usize = 128;
/// Maximum interval-sweep grid size.
pub const MAX_VOI_GRID: usize = 1024;
/// Maximum surrogate evaluations admitted by one public VoI call.
pub const MAX_VOI_EVALUATIONS: usize = 4096;
/// Maximum abstract oracle work units admitted by one public VoI call.
pub const MAX_VOI_WORK_UNITS: u64 = 1_000_000_000;
/// Maximum matched-cost observations admitted by one prospective audit.
pub const MAX_VOI_AUDIT_RECORDS: usize = 4096;
/// Maximum distinct ranking snapshots retained by one live scheduler.
pub const MAX_VOI_SCHEDULED_CONTEXTS: usize = 4096;
/// Fixed anytime-valid false-activation level for VoI scheduling authority.
pub const VOI_AUDIT_ALPHA: f64 = 0.05;
/// Canonical ledger-receipt format shared by ranked menus and audit outcomes.
pub const VOI_LEDGER_RECEIPT_VERSION: u32 = 1;
/// Maximum canonical bytes admitted for one VoI ledger receipt.
pub const MAX_VOI_LEDGER_RECEIPT_BYTES: usize = 64 * 1024;
/// Maximum canonical bytes admitted for one complete probe catalog artifact.
pub const MAX_VOI_PROBE_CATALOG_BYTES: usize = 512 * 1024;
/// Maximum bytes read to authenticate one decision/model identity artifact.
pub const MAX_VOI_AUTHORITY_ARTIFACT_BYTES: usize = 1024 * 1024;
/// Ledger artifact kind for the canonical decision/callback identity manifest.
pub const VOI_DECISION_ARTIFACT_KIND: &str = "fs-plan-voi-decision";
/// Ledger artifact kind for the canonical model/surrogate identity manifest.
pub const VOI_MODEL_ARTIFACT_KIND: &str = "fs-plan-voi-model";
/// Ledger artifact kind for [`voi_probe_catalog_bytes`].
pub const VOI_PROBE_CATALOG_ARTIFACT_KIND: &str = "fs-plan-voi-probe-catalog";
/// Exact metadata for a complete canonical probe catalog artifact.
pub const VOI_PROBE_CATALOG_ARTIFACT_METADATA: &str =
    "{\"schema\":\"fs-plan-voi-probe-catalog-v1\",\"complete\":true}";
/// Ledger artifact kind for [`voi_audit_outcome_bytes`].
pub const VOI_AUDIT_OUTCOME_ARTIFACT_KIND: &str = "fs-plan-voi-outcome";
/// Exact metadata for a canonical prospective outcome artifact.
pub const VOI_AUDIT_OUTCOME_ARTIFACT_METADATA: &str = "{\"schema\":\"fs-plan-voi-outcome-v1\"}";

/// Canonical semantics of the validated inputs to a ranked VoI menu.
pub const VOI_RANKED_SOURCE_IDENTITY_VERSION: u32 = 2;
/// Canonical semantics of the ranked output rows.
pub const VOI_RANKED_MENU_IDENTITY_VERSION: u32 = 2;
/// Canonical semantics of a chronological matched-audit prefix.
pub const VOI_AUDIT_CONTEXT_IDENTITY_VERSION: u32 = 2;

const RANKED_MENU_SOURCE_DOMAIN: &str = "frankensim.fs-plan.voi-ranked-source.v2";
const RANKED_MENU_CONTEXT_DOMAIN: &str = "frankensim.fs-plan.voi-ranked-menu.v2";
const AUDIT_CONTEXT_DOMAIN: &str = "frankensim.fs-plan.voi-audit.v2";
const RANKED_LEDGER_RECEIPT_MAGIC: &[u8] = b"fs-plan-voi-ranked-ledger-receipt-v1\0";
const AUDIT_LEDGER_RECEIPT_MAGIC: &[u8] = b"fs-plan-voi-audit-ledger-receipt-v1\0";
const PROBE_CATALOG_MAGIC: &[u8] = b"fs-plan-voi-probe-catalog-v1\0";
const AUDIT_OUTCOME_MAGIC: &[u8] = b"fs-plan-voi-audit-outcome-v1\0";
const RANKED_RECEIPT_KIND: &str = "fs-plan-voi-ranked-receipt";
const AUDIT_RECEIPT_KIND: &str = "fs-plan-voi-audit-receipt";
const RANKED_RECEIPT_META: &str =
    "{\"schema\":\"fs-plan-voi-ranked-ledger-receipt-v1\",\"trust\":\"ledger-lineage\"}";
const AUDIT_RECEIPT_META: &str =
    "{\"schema\":\"fs-plan-voi-audit-ledger-receipt-v1\",\"trust\":\"ledger-lineage\"}";
const RECEIPT_SEED: &[u8] = b"fs-plan-voi-ledger-receipt-v1";
const RECEIPT_VERSIONS: &str = "{\"fs-plan-voi-ledger-receipt\":1,\"ranked_source_identity\":2,\"ranked_menu_identity\":2,\"audit_context_identity\":2}";
const RECEIPT_BUDGET: &str = "{\"max_receipt_bytes\":65536,\"max_catalog_bytes\":524288,\"max_authority_artifact_bytes\":1048576}";
const RECEIPT_CAPABILITY: &str = "{\"ops\":[\"voi.receipt\"]}";
const NANOSECONDS_PER_DAY: i64 = 86_400_000_000_000;

/// Owner-local ranked-source declaration consumed by `xtask check-identities`.
pub const VOI_RANKED_SOURCE_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-plan:voi-ranked-source",
    "version_const=VOI_RANKED_SOURCE_IDENTITY_VERSION",
    "version=2",
    "domain=frankensim.fs-plan.voi-ranked-source.v2",
    "domain_const=RANKED_MENU_SOURCE_DOMAIN",
    "encoder=ranked_source_context",
    "encoder_helpers=ranked_source_context_with_schema,ranked_source_context_with_declared_counts,push_u32,push_text,ProbeKind::identity_tag",
    "schema_constants=RANKED_MENU_SOURCE_DOMAIN,VOI_RANKED_SOURCE_IDENTITY_VERSION,crates/fs-blake3/src/lib.rs#IV,crates/fs-blake3/src/lib.rs#MSG_PERMUTATION,crates/fs-blake3/src/lib.rs#BLOCK_LEN,crates/fs-blake3/src/lib.rs#CHUNK_LEN,crates/fs-blake3/src/lib.rs#CHUNK_START,crates/fs-blake3/src/lib.rs#CHUNK_END,crates/fs-blake3/src/lib.rs#PARENT,crates/fs-blake3/src/lib.rs#ROOT,crates/fs-blake3/src/lib.rs#DERIVE_KEY_CONTEXT,crates/fs-blake3/src/lib.rs#DERIVE_KEY_MATERIAL,crates/fs-blake3/src/lib.rs#MAX_DEPTH",
    "schema_functions=check_identity_version,check_ranked_source_identity_version,RankedMenu::source_identity_version,crates/fs-blake3/src/lib.rs#hash_domain,crates/fs-blake3/src/lib.rs#ContentHash::as_bytes,crates/fs-blake3/src/lib.rs#g,crates/fs-blake3/src/lib.rs#round,crates/fs-blake3/src/lib.rs#permute,crates/fs-blake3/src/lib.rs#compress,crates/fs-blake3/src/lib.rs#words_from_block,crates/fs-blake3/src/lib.rs#first_8_words,crates/fs-blake3/src/lib.rs#Output::chaining_value,crates/fs-blake3/src/lib.rs#Output::root_hash,crates/fs-blake3/src/lib.rs#parent_output,crates/fs-blake3/src/lib.rs#ChunkState::new,crates/fs-blake3/src/lib.rs#ChunkState::len,crates/fs-blake3/src/lib.rs#ChunkState::start_flag,crates/fs-blake3/src/lib.rs#ChunkState::update,crates/fs-blake3/src/lib.rs#ChunkState::output,crates/fs-blake3/src/lib.rs#Blake3::new_internal,crates/fs-blake3/src/lib.rs#Blake3::push_stack,crates/fs-blake3/src/lib.rs#Blake3::pop_stack,crates/fs-blake3/src/lib.rs#Blake3::add_chunk_chaining_value,crates/fs-blake3/src/lib.rs#Blake3::update,crates/fs-blake3/src/lib.rs#Blake3::finalize",
    "schema_dependencies=none",
    "digest=fs-blake3",
    "encoding=typed-binary",
    "sources=UncertaintyNode,Probe,DecisionOracleMetadata,DecisionComputationReceipt,DecisionBudget",
    "source_fields=UncertaintyNode.name:semantic,UncertaintyNode.lo:semantic,UncertaintyNode.hi:semantic,UncertaintyNode.nominal:semantic,Probe.name:semantic,Probe.target:semantic,Probe.cost:semantic,Probe.shrink:semantic,Probe.kind:semantic,DecisionOracleMetadata.arity:semantic,DecisionOracleMetadata.work_units_per_evaluation:semantic,DecisionComputationReceipt.evaluations:semantic,DecisionComputationReceipt.work_units:semantic,DecisionComputationReceipt.budget:derived:nested-budget-fields-encoded-separately,DecisionBudget.max_evaluations:semantic,DecisionBudget.max_work_units:semantic",
    "source_bindings=UncertaintyNode.name>node-name,UncertaintyNode.lo>node-lo,UncertaintyNode.hi>node-hi,UncertaintyNode.nominal>node-nominal,Probe.name>probe-name,Probe.target>probe-target,Probe.cost>probe-cost,Probe.shrink>probe-shrink,Probe.kind>probe-kind,DecisionOracleMetadata.arity>oracle-arity,DecisionOracleMetadata.work_units_per_evaluation>work-units-per-evaluation,DecisionComputationReceipt.evaluations>decision-evaluations,DecisionComputationReceipt.work_units>decision-work-units,DecisionBudget.max_evaluations>decision-evaluation-budget,DecisionBudget.max_work_units>decision-work-budget",
    "external_semantic_fields=artifact-domain,identity-version,policy-scope,snapshot-id,grid,node-count,node-order,probe-count",
    "semantic_fields=artifact-domain,identity-version,policy-scope,snapshot-id,grid,oracle-arity,work-units-per-evaluation,decision-evaluations,decision-work-units,decision-evaluation-budget,decision-work-budget,node-count,node-order,node-name,node-lo,node-nominal,node-hi,probe-count,probe-name,probe-target,probe-cost,probe-shrink,probe-kind",
    "excluded_fields=source-menu-input-order:canonicalized-by-validated-probe-name,allocation-capacity:representation-only",
    "consumers=rank_purchases,RankedMenu::source_context_id,RankedMenu::context_id,record_ranked_menu_receipt,verify_ranked_menu_receipt,VoiScheduler::schedule,VoiScheduler::schedule_authenticated",
    "mutations=artifact-domain:crates/fs-plan/src/voi.rs#voi_ranked_source_identity_mutation_battery,identity-version:crates/fs-plan/src/voi.rs#voi_ranked_source_identity_mutation_battery,policy-scope:crates/fs-plan/src/voi.rs#voi_ranked_source_identity_mutation_battery,snapshot-id:crates/fs-plan/src/voi.rs#voi_ranked_source_identity_mutation_battery,grid:crates/fs-plan/src/voi.rs#voi_ranked_source_identity_mutation_battery,oracle-arity:crates/fs-plan/src/voi.rs#voi_ranked_source_identity_mutation_battery,work-units-per-evaluation:crates/fs-plan/src/voi.rs#voi_ranked_source_identity_mutation_battery,decision-evaluations:crates/fs-plan/src/voi.rs#voi_ranked_source_identity_mutation_battery,decision-work-units:crates/fs-plan/src/voi.rs#voi_ranked_source_identity_mutation_battery,decision-evaluation-budget:crates/fs-plan/src/voi.rs#voi_ranked_source_identity_mutation_battery,decision-work-budget:crates/fs-plan/src/voi.rs#voi_ranked_source_identity_mutation_battery,node-count:crates/fs-plan/src/voi.rs#voi_ranked_source_identity_mutation_battery,node-order:crates/fs-plan/src/voi.rs#voi_ranked_source_identity_mutation_battery,node-name:crates/fs-plan/src/voi.rs#voi_ranked_source_identity_mutation_battery,node-lo:crates/fs-plan/src/voi.rs#voi_ranked_source_identity_mutation_battery,node-nominal:crates/fs-plan/src/voi.rs#voi_ranked_source_identity_mutation_battery,node-hi:crates/fs-plan/src/voi.rs#voi_ranked_source_identity_mutation_battery,probe-count:crates/fs-plan/src/voi.rs#voi_ranked_source_identity_mutation_battery,probe-name:crates/fs-plan/src/voi.rs#voi_ranked_source_identity_mutation_battery,probe-target:crates/fs-plan/src/voi.rs#voi_ranked_source_identity_mutation_battery,probe-cost:crates/fs-plan/src/voi.rs#voi_ranked_source_identity_mutation_battery,probe-shrink:crates/fs-plan/src/voi.rs#voi_ranked_source_identity_mutation_battery,probe-kind:crates/fs-plan/src/voi.rs#voi_ranked_source_identity_mutation_battery",
    "nonsemantic_mutations=source-menu-input-order:crates/fs-plan/src/voi.rs#voi_ranked_source_menu_input_order_is_nonsemantic,allocation-capacity:crates/fs-plan/src/voi.rs#voi_identity_allocation_capacity_is_nonsemantic",
    "field_guard=classify_voi_ranked_source_identity_fields",
    "transport_guard=ranked_source_context",
    "version_guard=crates/fs-plan/tests/voi.rs#voi_identity_versions_fail_closed",
    "coupling_surface=fs-plan:voi-ranked-source",
];

/// Owner-local ranked-menu declaration consumed by `xtask check-identities`.
pub const VOI_RANKED_MENU_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-plan:voi-ranked-menu",
    "version_const=VOI_RANKED_MENU_IDENTITY_VERSION",
    "version=2",
    "domain=frankensim.fs-plan.voi-ranked-menu.v2",
    "domain_const=RANKED_MENU_CONTEXT_DOMAIN",
    "encoder=ranked_output_context",
    "encoder_helpers=ranked_output_context_with_schema,ranked_output_context_with_declared_count,push_u32,push_text",
    "schema_constants=RANKED_MENU_CONTEXT_DOMAIN,VOI_RANKED_MENU_IDENTITY_VERSION,crates/fs-blake3/src/lib.rs#IV,crates/fs-blake3/src/lib.rs#MSG_PERMUTATION,crates/fs-blake3/src/lib.rs#BLOCK_LEN,crates/fs-blake3/src/lib.rs#CHUNK_LEN,crates/fs-blake3/src/lib.rs#CHUNK_START,crates/fs-blake3/src/lib.rs#CHUNK_END,crates/fs-blake3/src/lib.rs#PARENT,crates/fs-blake3/src/lib.rs#ROOT,crates/fs-blake3/src/lib.rs#DERIVE_KEY_CONTEXT,crates/fs-blake3/src/lib.rs#DERIVE_KEY_MATERIAL,crates/fs-blake3/src/lib.rs#MAX_DEPTH",
    "schema_functions=check_identity_version,check_ranked_source_identity_version,check_ranked_menu_identity_version,RankedMenu::source_identity_version,RankedMenu::identity_version,RankedMenu::admit_retained_identity_versions,crates/fs-blake3/src/lib.rs#hash_domain,crates/fs-blake3/src/lib.rs#ContentHash::as_bytes,crates/fs-blake3/src/lib.rs#g,crates/fs-blake3/src/lib.rs#round,crates/fs-blake3/src/lib.rs#permute,crates/fs-blake3/src/lib.rs#compress,crates/fs-blake3/src/lib.rs#words_from_block,crates/fs-blake3/src/lib.rs#first_8_words,crates/fs-blake3/src/lib.rs#Output::chaining_value,crates/fs-blake3/src/lib.rs#Output::root_hash,crates/fs-blake3/src/lib.rs#parent_output,crates/fs-blake3/src/lib.rs#ChunkState::new,crates/fs-blake3/src/lib.rs#ChunkState::len,crates/fs-blake3/src/lib.rs#ChunkState::start_flag,crates/fs-blake3/src/lib.rs#ChunkState::update,crates/fs-blake3/src/lib.rs#ChunkState::output,crates/fs-blake3/src/lib.rs#Blake3::new_internal,crates/fs-blake3/src/lib.rs#Blake3::push_stack,crates/fs-blake3/src/lib.rs#Blake3::pop_stack,crates/fs-blake3/src/lib.rs#Blake3::add_chunk_chaining_value,crates/fs-blake3/src/lib.rs#Blake3::update,crates/fs-blake3/src/lib.rs#Blake3::finalize",
    "schema_dependencies=fs-plan:voi-ranked-source",
    "digest=fs-blake3",
    "encoding=typed-binary",
    "sources=RankedPurchase,Probe",
    "source_fields=RankedPurchase.probe:derived:nested-probe-fields-classified-separately,RankedPurchase.flip_before:semantic,RankedPurchase.flip_after:semantic,RankedPurchase.score:semantic,Probe.name:semantic,Probe.target:derived:already-bound-by-source-context-root,Probe.cost:derived:already-bound-by-source-context-root,Probe.shrink:derived:already-bound-by-source-context-root,Probe.kind:derived:already-bound-by-source-context-root",
    "source_bindings=RankedPurchase.flip_before>flip-before,RankedPurchase.flip_after>flip-after,RankedPurchase.score>score,Probe.name>ranked-probe-name",
    "external_semantic_fields=artifact-domain,identity-version,source-context-root,row-count,row-order",
    "semantic_fields=artifact-domain,identity-version,source-context-root,row-count,row-order,ranked-probe-name,flip-before,flip-after,score",
    "excluded_fields=ranked-probe-payload:already-bound-by-source-context-root,allocation-capacity:representation-only",
    "consumers=RankedMenu::context_id,QueryHint,record_ranked_menu_receipt,verify_ranked_menu_receipt,VoiScheduler::schedule,VoiScheduler::schedule_authenticated",
    "mutations=artifact-domain:crates/fs-plan/src/voi.rs#voi_ranked_menu_identity_mutation_battery,identity-version:crates/fs-plan/src/voi.rs#voi_ranked_menu_identity_mutation_battery,source-context-root:crates/fs-plan/src/voi.rs#voi_ranked_menu_identity_mutation_battery,row-count:crates/fs-plan/src/voi.rs#voi_ranked_menu_identity_mutation_battery,row-order:crates/fs-plan/src/voi.rs#voi_ranked_menu_identity_mutation_battery,ranked-probe-name:crates/fs-plan/src/voi.rs#voi_ranked_menu_identity_mutation_battery,flip-before:crates/fs-plan/src/voi.rs#voi_ranked_menu_identity_mutation_battery,flip-after:crates/fs-plan/src/voi.rs#voi_ranked_menu_identity_mutation_battery,score:crates/fs-plan/src/voi.rs#voi_ranked_menu_identity_mutation_battery",
    "nonsemantic_mutations=ranked-probe-payload:crates/fs-plan/src/voi.rs#voi_ranked_menu_probe_payload_is_bound_by_source_context,allocation-capacity:crates/fs-plan/src/voi.rs#voi_identity_allocation_capacity_is_nonsemantic",
    "field_guard=classify_voi_ranked_menu_identity_fields",
    "transport_guard=ranked_output_context",
    "version_guard=crates/fs-plan/tests/voi.rs#voi_identity_versions_fail_closed",
    "coupling_surface=fs-plan:voi-ranked-menu",
];

/// Owner-local audit-context declaration consumed by `xtask check-identities`.
pub const VOI_AUDIT_CONTEXT_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-plan:voi-audit-context",
    "version_const=VOI_AUDIT_CONTEXT_IDENTITY_VERSION",
    "version=2",
    "domain=frankensim.fs-plan.voi-audit.v2",
    "domain_const=AUDIT_CONTEXT_DOMAIN",
    "encoder=audit_context",
    "encoder_helpers=audit_context_with_schema,audit_context_with_declared_count,push_u32,push_text",
    "schema_constants=AUDIT_CONTEXT_DOMAIN,VOI_AUDIT_ALPHA,VOI_AUDIT_CONTEXT_IDENTITY_VERSION,MAX_VOI_AUDIT_RECORDS,crates/fs-blake3/src/lib.rs#IV,crates/fs-blake3/src/lib.rs#MSG_PERMUTATION,crates/fs-blake3/src/lib.rs#BLOCK_LEN,crates/fs-blake3/src/lib.rs#CHUNK_LEN,crates/fs-blake3/src/lib.rs#CHUNK_START,crates/fs-blake3/src/lib.rs#CHUNK_END,crates/fs-blake3/src/lib.rs#PARENT,crates/fs-blake3/src/lib.rs#ROOT,crates/fs-blake3/src/lib.rs#DERIVE_KEY_CONTEXT,crates/fs-blake3/src/lib.rs#DERIVE_KEY_MATERIAL,crates/fs-blake3/src/lib.rs#MAX_DEPTH",
    "schema_functions=check_identity_version,check_audit_context_identity_version,AuditReport::identity_version,AuditReport::admit_retained_identity_version,crates/fs-blake3/src/lib.rs#hash_domain,crates/fs-blake3/src/lib.rs#ContentHash::as_bytes,crates/fs-blake3/src/lib.rs#g,crates/fs-blake3/src/lib.rs#round,crates/fs-blake3/src/lib.rs#permute,crates/fs-blake3/src/lib.rs#compress,crates/fs-blake3/src/lib.rs#words_from_block,crates/fs-blake3/src/lib.rs#first_8_words,crates/fs-blake3/src/lib.rs#Output::chaining_value,crates/fs-blake3/src/lib.rs#Output::root_hash,crates/fs-blake3/src/lib.rs#parent_output,crates/fs-blake3/src/lib.rs#ChunkState::new,crates/fs-blake3/src/lib.rs#ChunkState::len,crates/fs-blake3/src/lib.rs#ChunkState::start_flag,crates/fs-blake3/src/lib.rs#ChunkState::update,crates/fs-blake3/src/lib.rs#ChunkState::output,crates/fs-blake3/src/lib.rs#Blake3::new_internal,crates/fs-blake3/src/lib.rs#Blake3::push_stack,crates/fs-blake3/src/lib.rs#Blake3::pop_stack,crates/fs-blake3/src/lib.rs#Blake3::add_chunk_chaining_value,crates/fs-blake3/src/lib.rs#Blake3::update,crates/fs-blake3/src/lib.rs#Blake3::finalize",
    "schema_dependencies=none",
    "digest=fs-blake3",
    "encoding=typed-binary",
    "sources=MatchedAuditRecord",
    "source_fields=MatchedAuditRecord.observation_id:semantic,MatchedAuditRecord.recommended_id:semantic,MatchedAuditRecord.alternative_id:semantic,MatchedAuditRecord.provenance:semantic,MatchedAuditRecord.matched_cost:semantic,MatchedAuditRecord.recommended_changed_decision:semantic,MatchedAuditRecord.alternative_changed_decision:semantic",
    "source_bindings=MatchedAuditRecord.observation_id>observation-id,MatchedAuditRecord.recommended_id>recommended-id,MatchedAuditRecord.alternative_id>alternative-id,MatchedAuditRecord.provenance>provenance,MatchedAuditRecord.matched_cost>matched-cost,MatchedAuditRecord.recommended_changed_decision>recommended-changed-decision,MatchedAuditRecord.alternative_changed_decision>alternative-changed-decision",
    "external_semantic_fields=artifact-domain,identity-version,policy-scope,audit-alpha,max-audit-records,record-count,record-order",
    "semantic_fields=artifact-domain,identity-version,policy-scope,audit-alpha,max-audit-records,record-count,record-order,observation-id,recommended-id,alternative-id,provenance,matched-cost,recommended-changed-decision,alternative-changed-decision",
    "excluded_fields=allocation-capacity:representation-only,report-verdict:derived-from-chronological-e-process",
    "consumers=VoiScheduler::audit_report,VoiScheduler::schedule,VoiScheduler::schedule_authenticated,record_prospective_audit_receipt,verify_prospective_audit_receipt,audit_scheduling,AuditReport::audit_context_id",
    "mutations=artifact-domain:crates/fs-plan/src/voi.rs#voi_audit_context_identity_mutation_battery,identity-version:crates/fs-plan/src/voi.rs#voi_audit_context_identity_mutation_battery,policy-scope:crates/fs-plan/src/voi.rs#voi_audit_context_identity_mutation_battery,audit-alpha:crates/fs-plan/src/voi.rs#voi_audit_context_identity_mutation_battery,max-audit-records:crates/fs-plan/src/voi.rs#voi_audit_context_identity_mutation_battery,record-count:crates/fs-plan/src/voi.rs#voi_audit_context_identity_mutation_battery,record-order:crates/fs-plan/src/voi.rs#voi_audit_context_identity_mutation_battery,observation-id:crates/fs-plan/src/voi.rs#voi_audit_context_identity_mutation_battery,recommended-id:crates/fs-plan/src/voi.rs#voi_audit_context_identity_mutation_battery,alternative-id:crates/fs-plan/src/voi.rs#voi_audit_context_identity_mutation_battery,provenance:crates/fs-plan/src/voi.rs#voi_audit_context_identity_mutation_battery,matched-cost:crates/fs-plan/src/voi.rs#voi_audit_context_identity_mutation_battery,recommended-changed-decision:crates/fs-plan/src/voi.rs#voi_audit_context_identity_mutation_battery,alternative-changed-decision:crates/fs-plan/src/voi.rs#voi_audit_context_identity_mutation_battery",
    "nonsemantic_mutations=allocation-capacity:crates/fs-plan/src/voi.rs#voi_identity_allocation_capacity_is_nonsemantic,report-verdict:crates/fs-plan/src/voi.rs#voi_audit_report_verdict_is_derived_and_nonsemantic",
    "field_guard=classify_voi_audit_context_identity_fields",
    "transport_guard=audit_context",
    "version_guard=crates/fs-plan/tests/voi.rs#voi_identity_versions_fail_closed",
    "coupling_surface=fs-plan:voi-audit-context",
];

/// Why a VoI query, audit, or schedule was refused.
#[derive(Debug, Clone, PartialEq)]
pub enum VoiError {
    /// A bounded collection falls outside its admitted size range.
    SizeLimit {
        /// Collection being validated.
        collection: &'static str,
        /// Supplied element count.
        count: usize,
        /// Inclusive lower bound.
        min: usize,
        /// Inclusive upper bound.
        max: usize,
    },
    /// The surrogate's declared arity differs from the node vector.
    ArityMismatch {
        /// Declared surrogate arity.
        arity: usize,
        /// Supplied node count.
        node_count: usize,
    },
    /// A node/probe/target/audit identity is not bounded visible ASCII.
    InvalidName {
        /// Name category.
        kind: &'static str,
        /// Position in its collection.
        index: usize,
        /// Supplied UTF-8 byte length.
        bytes: usize,
        /// Inclusive byte limit.
        max_bytes: usize,
    },
    /// A supposedly unique name occurs more than once.
    DuplicateName {
        /// Name category.
        kind: &'static str,
        /// Duplicate value (already bounded by [`MAX_VOI_NAME_BYTES`]).
        name: String,
    },
    /// An interval is nonfinite, unordered, too wide for finite
    /// arithmetic, or excludes its nominal value.
    InvalidInterval {
        /// Node name.
        node: String,
        /// Lower endpoint.
        lo: f64,
        /// Nominal value.
        nominal: f64,
        /// Upper endpoint.
        hi: f64,
    },
    /// A surrogate returned a nonfinite decision margin.
    NonFiniteMargin {
        /// Returned margin.
        value: f64,
    },
    /// A sensitivity request names a missing node.
    NodeIndexOutOfRange {
        /// Supplied node index.
        node_idx: usize,
        /// Supplied node count.
        node_count: usize,
    },
    /// The sweep grid is zero or exceeds the declared cap.
    InvalidGrid {
        /// Supplied grid size.
        grid: usize,
        /// Inclusive upper bound.
        max: usize,
    },
    /// A request would exceed the surrogate-evaluation budget.
    EvaluationLimitExceeded {
        /// Required evaluations.
        requested: usize,
        /// Inclusive limit.
        max: usize,
    },
    /// A caller-supplied evaluation budget is zero or above the public cap.
    InvalidEvaluationBudget {
        /// Supplied evaluation limit.
        supplied: usize,
        /// Inclusive public cap.
        max: usize,
    },
    /// A caller-supplied work budget is zero or above the public cap.
    InvalidWorkBudget {
        /// Supplied abstract work-unit limit.
        supplied: u64,
        /// Inclusive public cap.
        max: u64,
    },
    /// The oracle's declared cost per evaluation is zero or outside the cap.
    InvalidOracleWorkUnits {
        /// Supplied abstract work units per evaluation.
        work_units_per_evaluation: u64,
        /// Inclusive upper bound.
        max: u64,
    },
    /// A request would exceed the explicit oracle work budget.
    WorkLimitExceeded {
        /// Required abstract work units.
        requested: u64,
        /// Inclusive caller-supplied limit.
        max: u64,
    },
    /// The asupersync context refused at an evaluation boundary.
    DecisionEvaluationCancelled,
    /// A probe has an invalid numeric field.
    InvalidProbeValue {
        /// Probe name.
        probe: String,
        /// Invalid field (`cost` or `shrink`).
        field: &'static str,
        /// Supplied value.
        value: f64,
    },
    /// A probe target resolves to zero or multiple nodes.
    TargetResolution {
        /// Probe name.
        probe: String,
        /// Requested target.
        target: String,
        /// Number of matching nodes.
        matches: usize,
    },
    /// The scheduling budget is nonfinite or negative.
    InvalidBudget {
        /// Supplied budget.
        budget: f64,
    },
    /// An audit observation has a malformed finite matched-cost pair.
    InvalidAuditCost {
        /// Observation identity.
        observation: String,
        /// Recommended-purchase cost.
        recommended_cost: f64,
        /// Alternative-purchase cost.
        alternative_cost: f64,
    },
    /// An audit compares a purchase with itself.
    InvalidAuditPair {
        /// Observation identity.
        observation: String,
    },
    /// An audit repeats an observation identity and could double-count evidence.
    DuplicateAuditObservation {
        /// Repeated observation identity.
        observation: String,
    },
    /// A ranked menu belongs to a different scheduling policy.
    PolicyScopeMismatch {
        /// Policy scope owned by the scheduler.
        expected: String,
        /// Policy scope bound into the ranked menu.
        actual: String,
    },
    /// One decision snapshot was already evaluated by this scheduler.
    RankingSnapshotAlreadyConsumed {
        /// Duplicate caller-declared decision/ledger snapshot.
        snapshot_id: String,
    },
    /// Scheduling was requested before the live audit crossed its threshold.
    MissingSchedulingAuthority,
    /// A receipt context token, date interval, or snapshot binding is invalid.
    InvalidReceiptContext {
        /// Field that failed validation.
        field: &'static str,
        /// Bounded diagnostic.
        detail: String,
    },
    /// A retained receipt or its ledger envelope differs from the expected
    /// canonical binding.
    ReceiptMismatch {
        /// Receipt component that disagreed.
        field: &'static str,
    },
    /// The receipt is presented before its declared issuance day.
    ReceiptNotYetValid {
        /// Declared issuance day.
        issued_day: u32,
        /// Caller-supplied verification day.
        current_day: u32,
    },
    /// The receipt is presented after its inclusive expiry day.
    ReceiptExpired {
        /// Inclusive last valid day.
        expires_day: u32,
        /// Caller-supplied verification day.
        current_day: u32,
    },
    /// The receipt's ledger branch no longer has the bound snapshot at its
    /// current head.
    StaleLedgerSnapshot {
        /// Bound snapshot root.
        expected: ContentHash,
        /// Current branch head, or `None` for an unknown/uncommitted branch.
        actual: Option<ContentHash>,
    },
    /// A receipt input is not a finished output artifact in its bound snapshot.
    SnapshotArtifactMissing {
        /// Semantic role of the missing artifact.
        role: &'static str,
        /// Missing content identity.
        artifact: ContentHash,
    },
    /// A prospective receipt is not the next chronological audit record.
    AuditSequenceMismatch {
        /// Next sequence required by the scheduler.
        expected: usize,
        /// Sequence bound into the receipt.
        actual: usize,
    },
    /// A ledger read/write or VCS lookup failed while issuing or verifying a
    /// receipt.
    LedgerEvidence {
        /// Operation being performed.
        operation: &'static str,
        /// Underlying bounded diagnostic.
        detail: String,
    },
    /// Receipt issuance requires ownership of its own transaction.
    ReceiptTransactionActive,
    /// Retained identity bytes declare semantics unknown to this build.
    UnsupportedIdentityVersion {
        /// Identity surface being admitted.
        identity: &'static str,
        /// Version declared by the retained producer.
        declared: u32,
        /// Exact version supported by this build.
        supported: u32,
    },
    /// Finite inputs could not produce a finite, monotone result.
    ArithmeticRefusal {
        /// Operation that failed.
        operation: &'static str,
        /// Bounded subject name.
        subject: String,
    },
}

impl core::fmt::Display for VoiError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::SizeLimit {
                collection,
                count,
                min,
                max,
            } => write!(
                f,
                "{collection} has {count} entries; expected an inclusive range of {min}..={max}"
            ),
            Self::ArityMismatch { arity, node_count } => write!(
                f,
                "surrogate declares arity {arity}, but {node_count} uncertainty node(s) were supplied"
            ),
            Self::InvalidName {
                kind,
                index,
                bytes,
                max_bytes,
            } => write!(
                f,
                "{kind} name at index {index} is not nonempty visible ASCII or is {bytes} bytes long (limit {max_bytes})"
            ),
            Self::DuplicateName { kind, name } => {
                write!(f, "duplicate {kind} name {name:?}")
            }
            Self::InvalidInterval {
                node,
                lo,
                nominal,
                hi,
            } => write!(
                f,
                "node {node:?} has invalid interval [{lo:?}, {hi:?}] with nominal {nominal:?}"
            ),
            Self::NonFiniteMargin { value } => {
                write!(f, "surrogate returned nonfinite margin {value:?}")
            }
            Self::NodeIndexOutOfRange {
                node_idx,
                node_count,
            } => write!(
                f,
                "node index {node_idx} is out of range for {node_count} node(s)"
            ),
            Self::InvalidGrid { grid, max } => {
                write!(f, "sweep grid {grid} is outside 1..={max}")
            }
            Self::EvaluationLimitExceeded { requested, max } => write!(
                f,
                "VoI request needs {requested} surrogate evaluations; the limit is {max}"
            ),
            Self::InvalidEvaluationBudget { supplied, max } => write!(
                f,
                "decision evaluation budget is {supplied}; expected 1..={max}"
            ),
            Self::InvalidWorkBudget { supplied, max } => {
                write!(f, "decision work budget is {supplied}; expected 1..={max}")
            }
            Self::InvalidOracleWorkUnits {
                work_units_per_evaluation,
                max,
            } => write!(
                f,
                "decision oracle declares {work_units_per_evaluation} work units per evaluation; expected 1..={max}"
            ),
            Self::WorkLimitExceeded { requested, max } => write!(
                f,
                "VoI request needs {requested} oracle work units; the explicit limit is {max}"
            ),
            Self::DecisionEvaluationCancelled => write!(
                f,
                "decision evaluation was cancelled or its asupersync budget was exhausted"
            ),
            Self::InvalidProbeValue {
                probe,
                field,
                value,
            } => write!(f, "probe {probe:?} has invalid {field} {value:?}"),
            Self::TargetResolution {
                probe,
                target,
                matches,
            } => write!(
                f,
                "probe {probe:?} target {target:?} resolves to {matches} uncertainty node(s), expected exactly one"
            ),
            Self::InvalidBudget { budget } => {
                write!(
                    f,
                    "probe budget must be finite and non-negative, got {budget:?}"
                )
            }
            Self::InvalidAuditCost {
                observation,
                recommended_cost,
                alternative_cost,
            } => write!(
                f,
                "audit observation {observation:?} requires equal finite positive matched costs, got {recommended_cost:?} and {alternative_cost:?}"
            ),
            Self::InvalidAuditPair { observation } => write!(
                f,
                "audit observation {observation:?} compares a purchase with itself"
            ),
            Self::DuplicateAuditObservation { observation } => write!(
                f,
                "audit observation {observation:?} appears more than once"
            ),
            Self::PolicyScopeMismatch { expected, actual } => write!(
                f,
                "ranked menu policy scope {actual:?} does not match scheduler scope {expected:?}"
            ),
            Self::RankingSnapshotAlreadyConsumed { snapshot_id } => write!(
                f,
                "decision snapshot {snapshot_id:?} was already consumed by this scheduler"
            ),
            Self::MissingSchedulingAuthority => write!(
                f,
                "the live VoI audit has not crossed the anytime-valid scheduling threshold"
            ),
            Self::InvalidReceiptContext { field, detail } => {
                write!(
                    f,
                    "VoI receipt context field {field:?} is invalid: {detail}"
                )
            }
            Self::ReceiptMismatch { field } => {
                write!(f, "VoI ledger receipt does not match its expected {field}")
            }
            Self::ReceiptNotYetValid {
                issued_day,
                current_day,
            } => write!(
                f,
                "VoI ledger receipt is not valid before day {issued_day}; verification used day {current_day}"
            ),
            Self::ReceiptExpired {
                expires_day,
                current_day,
            } => write!(
                f,
                "VoI ledger receipt expired after day {expires_day}; verification used day {current_day}"
            ),
            Self::StaleLedgerSnapshot { expected, actual } => write!(
                f,
                "VoI receipt binds ledger snapshot {} but the current branch head is {}",
                expected,
                actual.map_or_else(|| "<none>".to_string(), |root| root.to_string())
            ),
            Self::SnapshotArtifactMissing { role, artifact } => write!(
                f,
                "VoI receipt {role} artifact {artifact} is not a finished output in its bound snapshot"
            ),
            Self::AuditSequenceMismatch { expected, actual } => write!(
                f,
                "VoI audit receipt sequence is {actual}; the live prospective stream requires {expected}"
            ),
            Self::LedgerEvidence { operation, detail } => {
                write!(f, "VoI ledger evidence {operation} failed: {detail}")
            }
            Self::ReceiptTransactionActive => write!(
                f,
                "VoI receipt issuance owns an atomic ledger transaction; finish the caller transaction first"
            ),
            Self::UnsupportedIdentityVersion {
                identity,
                declared,
                supported,
            } => write!(
                f,
                "retained {identity} identity v{declared} is unsupported; this build accepts exactly v{supported}"
            ),
            Self::ArithmeticRefusal { operation, subject } => {
                write!(
                    f,
                    "{operation} for {subject:?} did not remain finite and monotone"
                )
            }
        }
    }
}

impl std::error::Error for VoiError {}

fn check_identity_version(
    identity: &'static str,
    declared: u32,
    supported: u32,
) -> Result<(), VoiError> {
    if declared == supported {
        Ok(())
    } else {
        Err(VoiError::UnsupportedIdentityVersion {
            identity,
            declared,
            supported,
        })
    }
}

/// Refuse retained ranked-source identities from stale or future schemas.
///
/// # Errors
/// [`VoiError::UnsupportedIdentityVersion`] unless `declared` is exactly the
/// current ranked-source identity version.
pub fn check_ranked_source_identity_version(declared: u32) -> Result<(), VoiError> {
    check_identity_version(
        "VoI ranked source",
        declared,
        VOI_RANKED_SOURCE_IDENTITY_VERSION,
    )
}

/// Refuse retained ranked-menu identities from stale or future schemas.
///
/// # Errors
/// [`VoiError::UnsupportedIdentityVersion`] unless `declared` is exactly the
/// current ranked-menu identity version.
pub fn check_ranked_menu_identity_version(declared: u32) -> Result<(), VoiError> {
    check_identity_version(
        "VoI ranked menu",
        declared,
        VOI_RANKED_MENU_IDENTITY_VERSION,
    )
}

/// Refuse retained audit-context identities from stale or future schemas.
///
/// # Errors
/// [`VoiError::UnsupportedIdentityVersion`] unless `declared` is exactly the
/// current audit-context identity version.
pub fn check_audit_context_identity_version(declared: u32) -> Result<(), VoiError> {
    check_identity_version(
        "VoI audit context",
        declared,
        VOI_AUDIT_CONTEXT_IDENTITY_VERSION,
    )
}

fn validate_receipt_token(field: &'static str, value: &str) -> Result<(), VoiError> {
    if value.is_empty()
        || value.len() > MAX_VOI_NAME_BYTES
        || !value.bytes().all(|byte| {
            byte.is_ascii_alphanumeric() || matches!(byte, b'-' | b'_' | b'.' | b':' | b'/')
        })
    {
        return Err(VoiError::InvalidReceiptContext {
            field,
            detail: format!(
                "expected 1..={MAX_VOI_NAME_BYTES} bytes from [A-Za-z0-9._:/-], got {} bytes",
                value.len()
            ),
        });
    }
    Ok(())
}

fn receipt_time_ns(day: u32) -> Result<i64, VoiError> {
    i64::from(day)
        .checked_mul(NANOSECONDS_PER_DAY)
        .ok_or_else(|| VoiError::InvalidReceiptContext {
            field: "issued_day",
            detail: format!("day {day} cannot be represented as a ledger nanosecond timestamp"),
        })
}

/// Immutable authority scope bound into every ledger-backed VoI receipt.
///
/// The snapshot root is a real [`Vcs`] commit on `branch`; issuance and ranked
/// scheduling additionally require it to remain that branch's current head.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VoiLedgerContext {
    session_id: String,
    decision_version: String,
    model_version: String,
    policy_version: String,
    branch: i64,
    snapshot_root: ContentHash,
    issued_day: u32,
    expires_day: u32,
}

impl VoiLedgerContext {
    /// Construct one bounded receipt scope.
    ///
    /// # Errors
    /// [`VoiError::InvalidReceiptContext`] for malformed tokens, a nonpositive
    /// branch, an all-zero snapshot root, an inverted validity interval, or an
    /// issuance day that cannot map to the ledger nanosecond domain.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        session_id: impl Into<String>,
        decision_version: impl Into<String>,
        model_version: impl Into<String>,
        policy_version: impl Into<String>,
        branch: i64,
        snapshot_root: ContentHash,
        issued_day: u32,
        expires_day: u32,
    ) -> Result<Self, VoiError> {
        let session_id = session_id.into();
        let decision_version = decision_version.into();
        let model_version = model_version.into();
        let policy_version = policy_version.into();
        for (field, value) in [
            ("session_id", session_id.as_str()),
            ("decision_version", decision_version.as_str()),
            ("model_version", model_version.as_str()),
            ("policy_version", policy_version.as_str()),
        ] {
            validate_receipt_token(field, value)?;
        }
        if branch <= 0 {
            return Err(VoiError::InvalidReceiptContext {
                field: "branch",
                detail: format!("expected a positive branch id, got {branch}"),
            });
        }
        if snapshot_root.as_bytes().iter().all(|byte| *byte == 0) {
            return Err(VoiError::InvalidReceiptContext {
                field: "snapshot_root",
                detail: "all-zero is not a committed ledger snapshot".to_string(),
            });
        }
        if issued_day > expires_day {
            return Err(VoiError::InvalidReceiptContext {
                field: "validity",
                detail: format!("issued day {issued_day} exceeds expiry day {expires_day}"),
            });
        }
        receipt_time_ns(issued_day)?;
        Ok(Self {
            session_id: session_id.as_str().to_owned(),
            decision_version: decision_version.as_str().to_owned(),
            model_version: model_version.as_str().to_owned(),
            policy_version: policy_version.as_str().to_owned(),
            branch,
            snapshot_root,
            issued_day,
            expires_day,
        })
    }

    /// Session whose decision/audit stream owns the receipt.
    #[must_use]
    pub fn session_id(&self) -> &str {
        &self.session_id
    }

    /// Exact decision implementation/version identity.
    #[must_use]
    pub fn decision_version(&self) -> &str {
        &self.decision_version
    }

    /// Exact model/surrogate version identity.
    #[must_use]
    pub fn model_version(&self) -> &str {
        &self.model_version
    }

    /// Scheduling/audit policy version.
    #[must_use]
    pub fn policy_version(&self) -> &str {
        &self.policy_version
    }

    /// Ledger branch containing the snapshot.
    #[must_use]
    pub fn branch(&self) -> i64 {
        self.branch
    }

    /// Bound semantic VCS root.
    #[must_use]
    pub fn snapshot_root(&self) -> ContentHash {
        self.snapshot_root
    }

    /// First valid logical day.
    #[must_use]
    pub fn issued_day(&self) -> u32 {
        self.issued_day
    }

    /// Inclusive last valid logical day.
    #[must_use]
    pub fn expires_day(&self) -> u32 {
        self.expires_day
    }
}

/// Durable locator for one canonical receipt artifact and its sole producer.
///
/// Locators are intentionally constructible for replay/import; they carry no
/// authority until an independent verifier checks the ledger artifact,
/// producer operation, lineage, snapshot, and expiry.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct VoiLedgerReceipt {
    operation: i64,
    artifact: ContentHash,
}

impl VoiLedgerReceipt {
    /// Reconstruct a retained locator. This does not verify or mint authority.
    #[must_use]
    pub const fn from_parts(operation: i64, artifact: ContentHash) -> Self {
        Self {
            operation,
            artifact,
        }
    }

    /// Sole ledger operation claimed to have produced the receipt.
    #[must_use]
    pub const fn operation(&self) -> i64 {
        self.operation
    }

    /// Content-addressed receipt artifact.
    #[must_use]
    pub const fn artifact(&self) -> ContentHash {
        self.artifact
    }
}

/// One uncertainty node touching a live decision: a named quantity the
/// ledger only knows to an interval.
#[derive(Debug, Clone, PartialEq)]
pub struct UncertaintyNode {
    /// Ledger name.
    pub name: String,
    /// Current uncertainty interval.
    pub lo: f64,
    /// Upper end.
    pub hi: f64,
    /// Nominal (decision-time) value.
    pub nominal: f64,
}

/// Explicit resource envelope for one decision-oracle query.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DecisionBudget {
    max_evaluations: usize,
    max_work_units: u64,
}

impl DecisionBudget {
    /// Construct a bounded decision-oracle envelope.
    ///
    /// # Errors
    /// [`VoiError`] when either limit is zero or exceeds the public cap.
    pub fn new(max_evaluations: usize, max_work_units: u64) -> Result<Self, VoiError> {
        if !(1..=MAX_VOI_EVALUATIONS).contains(&max_evaluations) {
            return Err(VoiError::InvalidEvaluationBudget {
                supplied: max_evaluations,
                max: MAX_VOI_EVALUATIONS,
            });
        }
        if !(1..=MAX_VOI_WORK_UNITS).contains(&max_work_units) {
            return Err(VoiError::InvalidWorkBudget {
                supplied: max_work_units,
                max: MAX_VOI_WORK_UNITS,
            });
        }
        Ok(Self {
            max_evaluations,
            max_work_units,
        })
    }

    /// Maximum oracle calls authorized by this envelope.
    #[must_use]
    pub const fn max_evaluations(self) -> usize {
        self.max_evaluations
    }

    /// Maximum abstract oracle work authorized by this envelope.
    #[must_use]
    pub const fn max_work_units(self) -> u64 {
        self.max_work_units
    }
}

/// Exact evaluation count and declared-work charge retained by a ranking.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DecisionComputationReceipt {
    evaluations: usize,
    work_units: u64,
    budget: DecisionBudget,
}

impl DecisionComputationReceipt {
    /// Number of completed oracle evaluations.
    #[must_use]
    pub const fn evaluations(self) -> usize {
        self.evaluations
    }

    /// Deterministic declared work charged for the evaluations.
    #[must_use]
    pub const fn work_units(self) -> u64 {
        self.work_units
    }

    /// Caller-supplied resource envelope.
    #[must_use]
    pub const fn budget(self) -> DecisionBudget {
        self.budget
    }
}

/// One library-issued permit for a deterministic oracle evaluation.
///
/// The private fields prevent callers from inventing permits. Charged work is
/// the oracle's declaration, not a wall-clock, allocation, or instruction
/// measurement; cooperative oracles must use [`Cx`] checkpoints internally.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DecisionEvaluationPermit {
    ordinal: usize,
    total_evaluations: usize,
    charged_work_units: u64,
    remaining_evaluations: usize,
    remaining_work_units: u64,
    envelope: DecisionBudget,
}

impl DecisionEvaluationPermit {
    /// Zero-based position in the canonical evaluation sequence.
    #[must_use]
    pub const fn ordinal(self) -> usize {
        self.ordinal
    }

    /// Exact number of evaluations admitted for the whole query.
    #[must_use]
    pub const fn total_evaluations(self) -> usize {
        self.total_evaluations
    }

    /// Declared work charged for this evaluation.
    #[must_use]
    pub const fn charged_work_units(self) -> u64 {
        self.charged_work_units
    }

    /// Admitted evaluations remaining after this evaluation.
    #[must_use]
    pub const fn remaining_evaluations(self) -> usize {
        self.remaining_evaluations
    }

    /// Caller-authorized work units remaining after this evaluation's charge.
    #[must_use]
    pub const fn remaining_work_units(self) -> u64 {
        self.remaining_work_units
    }

    /// Caller-supplied envelope enclosing the complete query.
    #[must_use]
    pub const fn envelope(self) -> DecisionBudget {
        self.envelope
    }
}

/// Fallible, cooperatively cancellation-aware decision surface consumed by
/// VoI ranking.
///
/// Implementations must checkpoint `cx` inside long-running work and must
/// refuse a permit whose charged work is insufficient. The library cannot
/// preempt an implementation that blocks or ignores this protocol.
pub trait DecisionOracle {
    /// Number of scalar inputs expected by the oracle. This metadata accessor
    /// must be bounded, deterministic, and side-effect free.
    fn arity(&self) -> usize;

    /// Deterministic abstract cost charged before every oracle call. This
    /// metadata accessor must be bounded and side-effect free.
    fn work_units_per_evaluation(&self) -> u64;

    /// Evaluate a decision margin under the supplied cancellation context and
    /// library-issued declared-work permit.
    fn evaluate(
        &self,
        cx: &Cx,
        permit: DecisionEvaluationPermit,
        values: &[f64],
    ) -> Result<f64, VoiError>;
}

/// Adapter for an already-cached, caller-budgeted synchronous margin.
///
/// This adapter charges one declared work unit and only checks cancellation
/// before and after the callback. It makes no time or memory bound for an
/// arbitrary closure. Long-running or fallible work must implement
/// [`DecisionOracle`] directly and checkpoint its [`Cx`] internally.
pub struct LiveDecision<'a> {
    /// The already-cached surrogate margin.
    pub margin: &'a dyn Fn(&[f64]) -> f64,
    /// Node count.
    pub arity: usize,
}

impl DecisionOracle for LiveDecision<'_> {
    fn arity(&self) -> usize {
        self.arity
    }

    fn work_units_per_evaluation(&self) -> u64 {
        1
    }

    fn evaluate(
        &self,
        _cx: &Cx,
        _permit: DecisionEvaluationPermit,
        values: &[f64],
    ) -> Result<f64, VoiError> {
        Ok((self.margin)(values))
    }
}

#[derive(Debug, Clone, Copy)]
struct DecisionOracleMetadata {
    arity: usize,
    work_units_per_evaluation: u64,
}

fn decision_oracle_metadata(
    decision: &dyn DecisionOracle,
) -> Result<DecisionOracleMetadata, VoiError> {
    let metadata = DecisionOracleMetadata {
        arity: decision.arity(),
        work_units_per_evaluation: decision.work_units_per_evaluation(),
    };
    if !(1..=MAX_VOI_WORK_UNITS).contains(&metadata.work_units_per_evaluation) {
        return Err(VoiError::InvalidOracleWorkUnits {
            work_units_per_evaluation: metadata.work_units_per_evaluation,
            max: MAX_VOI_WORK_UNITS,
        });
    }
    Ok(metadata)
}

fn validate_size(
    collection: &'static str,
    count: usize,
    min: usize,
    max: usize,
) -> Result<(), VoiError> {
    if (min..=max).contains(&count) {
        Ok(())
    } else {
        Err(VoiError::SizeLimit {
            collection,
            count,
            min,
            max,
        })
    }
}

fn validate_name(kind: &'static str, index: usize, name: &str) -> Result<(), VoiError> {
    if name.is_empty()
        || name.len() > MAX_VOI_NAME_BYTES
        || !name.bytes().all(|byte| byte.is_ascii_graphic())
    {
        Err(VoiError::InvalidName {
            kind,
            index,
            bytes: name.len(),
            max_bytes: MAX_VOI_NAME_BYTES,
        })
    } else {
        Ok(())
    }
}

fn validate_nodes(arity: usize, nodes: &[UncertaintyNode]) -> Result<(), VoiError> {
    validate_size("uncertainty nodes", nodes.len(), 1, MAX_VOI_NODES)?;
    if arity != nodes.len() {
        return Err(VoiError::ArityMismatch {
            arity,
            node_count: nodes.len(),
        });
    }
    let mut names = BTreeSet::new();
    for (index, node) in nodes.iter().enumerate() {
        validate_name("uncertainty node", index, &node.name)?;
        if !names.insert(node.name.as_str()) {
            return Err(VoiError::DuplicateName {
                kind: "uncertainty node",
                name: node.name.clone(),
            });
        }
        let width = node.hi - node.lo;
        if !node.lo.is_finite()
            || !node.hi.is_finite()
            || !node.nominal.is_finite()
            || node.lo > node.hi
            || node.nominal < node.lo
            || node.nominal > node.hi
            || !width.is_finite()
        {
            return Err(VoiError::InvalidInterval {
                node: node.name.clone(),
                lo: node.lo,
                nominal: node.nominal,
                hi: node.hi,
            });
        }
    }
    Ok(())
}

fn validate_grid(grid: usize) -> Result<(), VoiError> {
    if (1..=MAX_VOI_GRID).contains(&grid) {
        Ok(())
    } else {
        Err(VoiError::InvalidGrid {
            grid,
            max: MAX_VOI_GRID,
        })
    }
}

fn validate_evaluations(requested: usize) -> Result<(), VoiError> {
    if requested <= MAX_VOI_EVALUATIONS {
        Ok(())
    } else {
        Err(VoiError::EvaluationLimitExceeded {
            requested,
            max: MAX_VOI_EVALUATIONS,
        })
    }
}

fn admit_decision_computation(
    metadata: DecisionOracleMetadata,
    evaluations: usize,
    budget: DecisionBudget,
) -> Result<DecisionComputationReceipt, VoiError> {
    validate_evaluations(evaluations)?;
    if evaluations > budget.max_evaluations {
        return Err(VoiError::EvaluationLimitExceeded {
            requested: evaluations,
            max: budget.max_evaluations,
        });
    }
    let evaluations_u64 = u64::try_from(evaluations).map_err(|_| VoiError::ArithmeticRefusal {
        operation: "decision evaluation count conversion",
        subject: "decision oracle".to_string(),
    })?;
    let work_units = evaluations_u64
        .checked_mul(metadata.work_units_per_evaluation)
        .ok_or_else(|| VoiError::ArithmeticRefusal {
            operation: "decision oracle work accounting",
            subject: "decision oracle".to_string(),
        })?;
    if work_units > budget.max_work_units {
        return Err(VoiError::WorkLimitExceeded {
            requested: work_units,
            max: budget.max_work_units,
        });
    }
    Ok(DecisionComputationReceipt {
        evaluations,
        work_units,
        budget,
    })
}

struct DecisionComputationMeter {
    receipt: DecisionComputationReceipt,
    work_units_per_evaluation: u64,
    completed_evaluations: usize,
    charged_work_units: u64,
}

impl DecisionComputationMeter {
    const fn new(metadata: DecisionOracleMetadata, receipt: DecisionComputationReceipt) -> Self {
        Self {
            receipt,
            work_units_per_evaluation: metadata.work_units_per_evaluation,
            completed_evaluations: 0,
            charged_work_units: 0,
        }
    }

    fn next_permit(&mut self) -> Result<DecisionEvaluationPermit, VoiError> {
        if self.completed_evaluations >= self.receipt.evaluations {
            return Err(VoiError::ArithmeticRefusal {
                operation: "decision evaluation permit overrun",
                subject: "decision oracle".to_string(),
            });
        }
        let ordinal = self.completed_evaluations;
        let completed_evaluations =
            ordinal
                .checked_add(1)
                .ok_or_else(|| VoiError::ArithmeticRefusal {
                    operation: "decision evaluation permit accounting",
                    subject: "decision oracle".to_string(),
                })?;
        let charged_work_units = self
            .charged_work_units
            .checked_add(self.work_units_per_evaluation)
            .ok_or_else(|| VoiError::ArithmeticRefusal {
                operation: "decision evaluation permit work accounting",
                subject: "decision oracle".to_string(),
            })?;
        let remaining_evaluations = self
            .receipt
            .evaluations
            .checked_sub(completed_evaluations)
            .ok_or_else(|| VoiError::ArithmeticRefusal {
                operation: "decision evaluation permit remainder",
                subject: "decision oracle".to_string(),
            })?;
        let remaining_work_units = self
            .receipt
            .budget
            .max_work_units
            .checked_sub(charged_work_units)
            .ok_or_else(|| VoiError::ArithmeticRefusal {
                operation: "decision evaluation permit work remainder",
                subject: "decision oracle".to_string(),
            })?;
        self.completed_evaluations = completed_evaluations;
        self.charged_work_units = charged_work_units;
        Ok(DecisionEvaluationPermit {
            ordinal,
            total_evaluations: self.receipt.evaluations,
            charged_work_units: self.work_units_per_evaluation,
            remaining_evaluations,
            remaining_work_units,
            envelope: self.receipt.budget,
        })
    }

    fn finish(&self) -> Result<(), VoiError> {
        if self.completed_evaluations == self.receipt.evaluations
            && self.charged_work_units == self.receipt.work_units
        {
            Ok(())
        } else {
            Err(VoiError::ArithmeticRefusal {
                operation: "decision evaluation permit underrun",
                subject: "decision oracle".to_string(),
            })
        }
    }
}

fn evaluate_margin(
    cx: &Cx,
    decision: &dyn DecisionOracle,
    meter: &mut DecisionComputationMeter,
    values: &[f64],
) -> Result<f64, VoiError> {
    cx.checkpoint()
        .map_err(|_| VoiError::DecisionEvaluationCancelled)?;
    let permit = meter.next_permit()?;
    let margin = decision.evaluate(cx, permit, values);
    cx.checkpoint()
        .map_err(|_| VoiError::DecisionEvaluationCancelled)?;
    let margin = margin?;
    if margin.is_finite() {
        Ok(margin)
    } else {
        Err(VoiError::NonFiniteMargin { value: margin })
    }
}

fn nominal_values(nodes: &[UncertaintyNode]) -> Vec<f64> {
    nodes.iter().map(|node| node.nominal).collect()
}

fn node_at(nodes: &[UncertaintyNode], node_idx: usize) -> Result<&UncertaintyNode, VoiError> {
    nodes.get(node_idx).ok_or(VoiError::NodeIndexOutOfRange {
        node_idx,
        node_count: nodes.len(),
    })
}

fn nominal_verdict_validated(
    cx: &Cx,
    decision: &dyn DecisionOracle,
    meter: &mut DecisionComputationMeter,
    nodes: &[UncertaintyNode],
) -> Result<bool, VoiError> {
    Ok(evaluate_margin(cx, decision, meter, &nominal_values(nodes))? > 0.0)
}

#[allow(clippy::too_many_arguments)] // keeps the validated sweep inputs explicit
fn flip_probability_validated(
    cx: &Cx,
    decision: &dyn DecisionOracle,
    meter: &mut DecisionComputationMeter,
    nodes: &[UncertaintyNode],
    base: bool,
    node_idx: usize,
    lo: f64,
    hi: f64,
    grid: usize,
) -> Result<f64, VoiError> {
    let mut values = nominal_values(nodes);
    let node = node_at(nodes, node_idx)?;
    let width = hi - lo;
    let mut flips = 0usize;
    for k in 0..grid {
        #[allow(clippy::cast_precision_loss)]
        let t = (k as f64 + 0.5) / grid as f64;
        let sample = lo + t * width;
        if !sample.is_finite() {
            return Err(VoiError::ArithmeticRefusal {
                operation: "interval sweep",
                subject: node.name.clone(),
            });
        }
        *values
            .get_mut(node_idx)
            .ok_or(VoiError::NodeIndexOutOfRange {
                node_idx,
                node_count: nodes.len(),
            })? = sample;
        if (evaluate_margin(cx, decision, meter, &values)? > 0.0) != base {
            flips += 1;
        }
    }
    #[allow(clippy::cast_precision_loss)]
    let probability = flips as f64 / grid as f64;
    Ok(probability)
}

impl LiveDecision<'_> {
    /// The nominal verdict.
    ///
    /// # Errors
    /// [`VoiError`] when node/arity/interval invariants fail or the
    /// cached surrogate returns a nonfinite margin.
    pub fn nominal_verdict(
        &self,
        cx: &Cx,
        nodes: &[UncertaintyNode],
        budget: DecisionBudget,
    ) -> Result<bool, VoiError> {
        let metadata = decision_oracle_metadata(self)?;
        validate_nodes(metadata.arity, nodes)?;
        let computation = admit_decision_computation(metadata, 1, budget)?;
        let mut meter = DecisionComputationMeter::new(metadata, computation);
        let verdict = nominal_verdict_validated(cx, self, &mut meter, nodes)?;
        meter.finish()?;
        Ok(verdict)
    }

    /// DECISION SENSITIVITY of one node: sweep the node's interval on
    /// the cached surrogate (others at nominal, `grid` points) and
    /// return the fraction of MIDPOINT GRID SAMPLES where the verdict differs
    /// from nominal. This is a myopic estimate under the uniform interval
    /// measure, not a certified probability.
    ///
    /// # Errors
    /// [`VoiError`] when the request is malformed, exceeds the declared
    /// sweep/evaluation limits, or the surrogate returns a nonfinite
    /// margin.
    pub fn flip_probability(
        &self,
        cx: &Cx,
        nodes: &[UncertaintyNode],
        node_idx: usize,
        grid: usize,
        budget: DecisionBudget,
    ) -> Result<f64, VoiError> {
        let metadata = decision_oracle_metadata(self)?;
        validate_nodes(metadata.arity, nodes)?;
        let node = node_at(nodes, node_idx)?;
        validate_grid(grid)?;
        let evaluations = grid
            .checked_add(1)
            .ok_or_else(|| VoiError::ArithmeticRefusal {
                operation: "sweep evaluation count",
                subject: node.name.clone(),
            })?;
        let computation = admit_decision_computation(metadata, evaluations, budget)?;
        let mut meter = DecisionComputationMeter::new(metadata, computation);
        let base = nominal_verdict_validated(cx, self, &mut meter, nodes)?;
        let probability = flip_probability_validated(
            cx, self, &mut meter, nodes, base, node_idx, node.lo, node.hi, grid,
        )?;
        meter.finish()?;
        Ok(probability)
    }
}

/// The kind of evidence purchase — the menu UNIFIES computational and
/// physical experiments (the epistemic-engine identity made concrete).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProbeKind {
    /// Climb a fidelity rung / refine / add solver accuracy.
    Computational,
    /// Wind-tunnel anchor, CT scan, strain gauge — reality as evidence.
    Physical,
}

impl ProbeKind {
    const fn identity_tag(self) -> u8 {
        match self {
            Self::Computational => 0,
            Self::Physical => 1,
        }
    }
}

/// One priced probe: buying it SHRINKS a node's interval around its
/// nominal by `shrink` (0 < shrink < 1).
#[derive(Debug, Clone, PartialEq)]
pub struct Probe {
    /// Menu name ("climb-to-rung-96", "wind-tunnel-anchor", ...).
    pub name: String,
    /// Which node it tightens.
    pub target: String,
    /// Price in dollars.
    pub cost: f64,
    /// Post-probe interval width as a fraction of the current width.
    pub shrink: f64,
    /// Computational or physical.
    pub kind: ProbeKind,
}

/// One ranked purchase: the myopic VoI score.
#[derive(Debug, Clone, PartialEq)]
pub struct RankedPurchase {
    /// The probe.
    probe: Probe,
    /// Grid-sampled flip fraction before the purchase.
    flip_before: f64,
    /// Grid-sampled flip fraction after the declared contraction.
    flip_after: f64,
    /// THE SCORE: sampled flip-fraction reduction per dollar.
    score: f64,
}

impl RankedPurchase {
    /// The validated probe purchase.
    #[must_use]
    pub fn probe(&self) -> &Probe {
        &self.probe
    }

    /// Grid-sampled flip fraction before the purchase.
    #[must_use]
    pub fn flip_before(&self) -> f64 {
        self.flip_before
    }

    /// Grid-sampled flip fraction after the declared contraction.
    #[must_use]
    pub fn flip_after(&self) -> f64 {
        self.flip_after
    }

    /// Grid-sampled flip-fraction reduction per dollar.
    #[must_use]
    pub fn score(&self) -> f64 {
        self.score
    }
}

#[allow(dead_code)]
fn classify_voi_ranked_source_identity_fields(
    node: &UncertaintyNode,
    probe: &Probe,
    metadata: &DecisionOracleMetadata,
    computation: &DecisionComputationReceipt,
    budget: &DecisionBudget,
) {
    let UncertaintyNode {
        name: _,
        lo: _,
        hi: _,
        nominal: _,
    } = node;
    let Probe {
        name: _,
        target: _,
        cost: _,
        shrink: _,
        kind: _,
    } = probe;
    let DecisionOracleMetadata {
        arity: _,
        work_units_per_evaluation: _,
    } = metadata;
    let DecisionComputationReceipt {
        evaluations: _,
        work_units: _,
        budget: _,
    } = computation;
    let DecisionBudget {
        max_evaluations: _,
        max_work_units: _,
    } = budget;
}

#[allow(dead_code)]
fn classify_voi_ranked_menu_identity_fields(row: &RankedPurchase, probe: &Probe) {
    let RankedPurchase {
        probe: _,
        flip_before: _,
        flip_after: _,
        score: _,
    } = row;
    let Probe {
        name: _,
        target: _,
        cost: _,
        shrink: _,
        kind: _,
    } = probe;
}

/// A complete, canonical ranking for one validated supplied
/// uncertainty/menu/grid snapshot. Rows and context are private so safe callers
/// cannot omit, splice, or reorder rows after ranking. This value alone does
/// not authenticate the caller-declared provenance.
#[derive(Debug, PartialEq)]
pub struct RankedMenu {
    rows: Vec<RankedPurchase>,
    source_context_id: ContentHash,
    context_id: ContentHash,
    policy_scope: String,
    snapshot_id: String,
    grid: usize,
    computation: DecisionComputationReceipt,
}

impl RankedMenu {
    /// Number of ranked purchases in the complete menu.
    #[must_use]
    pub fn len(&self) -> usize {
        self.rows.len()
    }

    /// A ranked menu produced by [`rank_purchases`] is never empty.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.rows.is_empty()
    }

    /// Inspect one canonical row without exposing mutable membership.
    #[must_use]
    pub fn get(&self, index: usize) -> Option<&RankedPurchase> {
        self.rows.get(index)
    }

    /// Inspect the highest-ranked purchase.
    #[must_use]
    pub fn top(&self) -> Option<&RankedPurchase> {
        self.rows.first()
    }

    /// Iterate over canonical rows for reporting only.
    #[must_use]
    pub fn iter(&self) -> impl ExactSizeIterator<Item = &RankedPurchase> {
        self.rows.iter()
    }

    /// Midpoint grid used for every sampled flip estimate in this menu.
    #[must_use]
    pub fn grid(&self) -> usize {
        self.grid
    }

    /// BLAKE3 identity of the validated node/menu/grid snapshot.
    ///
    /// This binds supplied content but does not identify callback code, prove
    /// catalog completeness, or prove that the snapshot remains current;
    /// callers must compare it with their ledger/session snapshot before use.
    #[must_use]
    pub fn context_id(&self) -> ContentHash {
        self.context_id
    }

    /// Identity of policy, snapshot, nodes, source menu, and grid before
    /// evaluating the decision surrogate.
    #[must_use]
    pub fn source_context_id(&self) -> ContentHash {
        self.source_context_id
    }

    /// Ranked-source identity semantics used for [`Self::source_context_id`].
    #[must_use]
    pub const fn source_identity_version(&self) -> u32 {
        VOI_RANKED_SOURCE_IDENTITY_VERSION
    }

    /// Ranked-output identity semantics used for [`Self::context_id`].
    #[must_use]
    pub const fn identity_version(&self) -> u32 {
        VOI_RANKED_MENU_IDENTITY_VERSION
    }

    /// Admit versions carried beside retained ranked-menu roots.
    ///
    /// Safe `RankedMenu` values are sealed by [`rank_purchases`]; this gate is
    /// for callers that retain the two roots with their producer-declared
    /// versions. Root comparison must follow this fail-closed version check.
    ///
    /// # Errors
    /// [`VoiError::UnsupportedIdentityVersion`] for any stale or future
    /// ranked-source or ranked-output version.
    pub fn admit_retained_identity_versions(
        &self,
        declared_source_version: u32,
        declared_menu_version: u32,
    ) -> Result<(), VoiError> {
        check_identity_version(
            "VoI ranked source",
            declared_source_version,
            self.source_identity_version(),
        )?;
        check_identity_version(
            "VoI ranked menu",
            declared_menu_version,
            self.identity_version(),
        )
    }

    /// Caller-declared policy/version scope bound into this ranking.
    #[must_use]
    pub fn policy_scope(&self) -> &str {
        &self.policy_scope
    }

    /// Caller-declared decision/ledger snapshot bound into this ranking.
    #[must_use]
    pub fn snapshot_id(&self) -> &str {
        &self.snapshot_id
    }

    /// Exact admitted decision-oracle resource use and enclosing budget.
    #[must_use]
    pub const fn computation(&self) -> DecisionComputationReceipt {
        self.computation
    }
}

/// Structured, grid-qualified query hint. Its private optional purchase keeps
/// the no-sampled-change state distinct from an authoritative zero claim.
#[derive(Debug, Clone, PartialEq)]
pub struct QueryHint {
    context_id: ContentHash,
    grid: usize,
    purchase: Option<RankedPurchase>,
}

impl QueryHint {
    /// Ranked snapshot identity supporting this estimate.
    #[must_use]
    pub fn context_id(&self) -> ContentHash {
        self.context_id
    }

    /// Midpoint grid supporting this estimate.
    #[must_use]
    pub fn grid(&self) -> usize {
        self.grid
    }

    /// Estimated top purchase, or `None` when no sampled row changed the
    /// decision on this grid. `None` is not a proof that no purchase can help.
    #[must_use]
    pub fn purchase(&self) -> Option<&RankedPurchase> {
        self.purchase.as_ref()
    }

    /// Safe deterministic text. Identifiers are escaped and every finite
    /// scalar uses Rust's shortest round-tripping representation.
    #[must_use]
    pub fn render_text(&self) -> String {
        match &self.purchase {
            Some(top) => format!(
                "estimated top evidence on the supplied menu from a {}-point midpoint sweep: {} (${}) - sampled flip fraction {} -> {} on {} ({}/$)",
                self.grid,
                escape_text(&top.probe.name),
                top.probe.cost,
                top.flip_before,
                top.flip_after,
                escape_text(&top.probe.target),
                top.score,
            ),
            None => format!(
                "no sampled purchase changed the decision on the {}-point midpoint sweep; this estimate does not prove that further evidence has zero value",
                self.grid
            ),
        }
    }

    /// Strict JSON rendering for logs and evidence payloads.
    #[must_use]
    pub fn to_json(&self) -> String {
        let context = self.context_id.to_hex();
        match &self.purchase {
            Some(top) => format!(
                "{{\"schema\":\"fs-plan.voi-hint.v1\",\"kind\":\"estimated_purchase\",\"context\":\"{context}\",\"grid\":{},\"probe\":{},\"target\":{},\"cost_dollars\":{},\"sampled_flip_before\":{},\"sampled_flip_after\":{},\"score_per_dollar\":{}}}",
                self.grid,
                json_string(&top.probe.name),
                json_string(&top.probe.target),
                top.probe.cost,
                top.flip_before,
                top.flip_after,
                top.score,
            ),
            None => format!(
                "{{\"schema\":\"fs-plan.voi-hint.v1\",\"kind\":\"no_sampled_change\",\"context\":\"{context}\",\"grid\":{},\"authoritative_zero\":false}}",
                self.grid
            ),
        }
    }
}

impl core::fmt::Display for QueryHint {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str(&self.render_text())
    }
}

fn escape_text(value: &str) -> String {
    value.chars().flat_map(char::escape_default).collect()
}

fn json_string(value: &str) -> String {
    let mut out = String::with_capacity(value.len() + 2);
    out.push('"');
    for byte in value.bytes() {
        match byte {
            b'"' => out.push_str("\\\""),
            b'\\' => out.push_str("\\\\"),
            _ => out.push(char::from(byte)),
        }
    }
    out.push('"');
    out
}

fn compare_ranked(a: &RankedPurchase, b: &RankedPurchase) -> core::cmp::Ordering {
    b.score
        .total_cmp(&a.score)
        .then(a.probe.cost.total_cmp(&b.probe.cost))
        .then(a.probe.name.cmp(&b.probe.name))
}

fn validate_probe(index: usize, probe: &Probe) -> Result<(), VoiError> {
    validate_name("probe", index, &probe.name)?;
    validate_name("probe target", index, &probe.target)?;
    if !probe.cost.is_finite() || probe.cost <= 0.0 {
        return Err(VoiError::InvalidProbeValue {
            probe: probe.name.clone(),
            field: "cost",
            value: probe.cost,
        });
    }
    if !probe.shrink.is_finite() || probe.shrink <= 0.0 || probe.shrink >= 1.0 {
        return Err(VoiError::InvalidProbeValue {
            probe: probe.name.clone(),
            field: "shrink",
            value: probe.shrink,
        });
    }
    Ok(())
}

fn validate_menu(nodes: &[UncertaintyNode], menu: &[Probe]) -> Result<Vec<usize>, VoiError> {
    validate_size("probe menu", menu.len(), 1, MAX_VOI_PROBES)?;
    let mut names = BTreeSet::new();
    let mut targets = Vec::with_capacity(menu.len());
    for (index, probe) in menu.iter().enumerate() {
        validate_probe(index, probe)?;
        if !names.insert(probe.name.as_str()) {
            return Err(VoiError::DuplicateName {
                kind: "probe",
                name: probe.name.clone(),
            });
        }
        let mut matched = None;
        let mut matches = 0usize;
        for (node_idx, node) in nodes.iter().enumerate() {
            if node.name == probe.target {
                matches += 1;
                matched = Some(node_idx);
            }
        }
        let Some(node_idx) = matched.filter(|_| matches == 1) else {
            return Err(VoiError::TargetResolution {
                probe: probe.name.clone(),
                target: probe.target.clone(),
                matches,
            });
        };
        targets.push(node_idx);
    }
    Ok(targets)
}

fn push_u32(out: &mut Vec<u8>, value: usize, subject: &'static str) -> Result<(), VoiError> {
    let value = u32::try_from(value).map_err(|_| VoiError::ArithmeticRefusal {
        operation: "VoI context length",
        subject: subject.to_string(),
    })?;
    out.extend_from_slice(&value.to_le_bytes());
    Ok(())
}

fn push_text(out: &mut Vec<u8>, value: &str, subject: &'static str) -> Result<(), VoiError> {
    push_u32(out, value.len(), subject)?;
    out.extend_from_slice(value.as_bytes());
    Ok(())
}

/// Canonical complete probe-catalog artifact bytes used by ledger receipts.
///
/// Input order is nonsemantic: unique validated probes are sorted by name and
/// encoded with an explicit count and exact floating-point bits.
///
/// # Errors
/// [`VoiError`] for an empty/oversized catalog, malformed probe, duplicate
/// identity, or an impossible receipt-size bound.
pub fn voi_probe_catalog_bytes(menu: &[Probe]) -> Result<Vec<u8>, VoiError> {
    validate_size("VoI probe menu", menu.len(), 1, MAX_VOI_PROBES)?;
    let mut names = BTreeSet::new();
    for (index, probe) in menu.iter().enumerate() {
        validate_probe(index, probe)?;
        validate_name("probe target", index, &probe.target)?;
        if !names.insert(probe.name.as_str()) {
            return Err(VoiError::DuplicateName {
                kind: "probe",
                name: probe.name.clone(),
            });
        }
    }
    let mut canonical_menu = menu.iter().collect::<Vec<_>>();
    canonical_menu.sort_by(|left, right| left.name.cmp(&right.name));
    let mut canonical = Vec::new();
    canonical.extend_from_slice(PROBE_CATALOG_MAGIC);
    canonical.extend_from_slice(&VOI_LEDGER_RECEIPT_VERSION.to_le_bytes());
    push_u32(&mut canonical, canonical_menu.len(), "probe catalog")?;
    for probe in canonical_menu {
        push_text(&mut canonical, &probe.name, "probe name")?;
        push_text(&mut canonical, &probe.target, "probe target")?;
        canonical.extend_from_slice(&probe.cost.to_bits().to_le_bytes());
        canonical.extend_from_slice(&probe.shrink.to_bits().to_le_bytes());
        canonical.push(probe.kind.identity_tag());
    }
    if canonical.len() > MAX_VOI_PROBE_CATALOG_BYTES {
        return Err(VoiError::SizeLimit {
            collection: "VoI probe catalog bytes",
            count: canonical.len(),
            min: 1,
            max: MAX_VOI_PROBE_CATALOG_BYTES,
        });
    }
    Ok(canonical)
}

fn ranked_source_context(
    nodes: &[UncertaintyNode],
    menu: &[Probe],
    grid: usize,
    policy_scope: &str,
    snapshot_id: &str,
    metadata: DecisionOracleMetadata,
    computation: DecisionComputationReceipt,
) -> Result<ContentHash, VoiError> {
    ranked_source_context_with_schema(
        RANKED_MENU_SOURCE_DOMAIN,
        VOI_RANKED_SOURCE_IDENTITY_VERSION,
        nodes,
        menu,
        grid,
        policy_scope,
        snapshot_id,
        metadata,
        computation,
    )
}

#[allow(clippy::too_many_arguments)]
fn ranked_source_context_with_schema(
    domain: &str,
    producer_version: u32,
    nodes: &[UncertaintyNode],
    menu: &[Probe],
    grid: usize,
    policy_scope: &str,
    snapshot_id: &str,
    metadata: DecisionOracleMetadata,
    computation: DecisionComputationReceipt,
) -> Result<ContentHash, VoiError> {
    ranked_source_context_with_declared_counts(
        domain,
        producer_version,
        nodes,
        nodes.len(),
        menu,
        menu.len(),
        grid,
        policy_scope,
        snapshot_id,
        metadata,
        computation,
    )
}

#[allow(clippy::too_many_arguments)]
fn ranked_source_context_with_declared_counts(
    domain: &str,
    producer_version: u32,
    nodes: &[UncertaintyNode],
    declared_node_count: usize,
    menu: &[Probe],
    declared_probe_count: usize,
    grid: usize,
    policy_scope: &str,
    snapshot_id: &str,
    metadata: DecisionOracleMetadata,
    computation: DecisionComputationReceipt,
) -> Result<ContentHash, VoiError> {
    let mut canonical = Vec::new();
    canonical.extend_from_slice(&producer_version.to_le_bytes());
    push_text(&mut canonical, policy_scope, "VoI policy scope")?;
    push_text(&mut canonical, snapshot_id, "VoI snapshot identity")?;
    push_u32(&mut canonical, grid, "grid")?;
    push_u32(&mut canonical, metadata.arity, "decision oracle arity")?;
    canonical.extend_from_slice(&metadata.work_units_per_evaluation.to_le_bytes());
    push_u32(
        &mut canonical,
        computation.evaluations,
        "decision evaluations",
    )?;
    canonical.extend_from_slice(&computation.work_units.to_le_bytes());
    push_u32(
        &mut canonical,
        computation.budget.max_evaluations,
        "decision evaluation budget",
    )?;
    canonical.extend_from_slice(&computation.budget.max_work_units.to_le_bytes());
    push_u32(&mut canonical, declared_node_count, "uncertainty nodes")?;
    for node in nodes {
        push_text(&mut canonical, &node.name, "uncertainty node name")?;
        canonical.extend_from_slice(&node.lo.to_bits().to_le_bytes());
        canonical.extend_from_slice(&node.nominal.to_bits().to_le_bytes());
        canonical.extend_from_slice(&node.hi.to_bits().to_le_bytes());
    }
    let mut canonical_menu: Vec<&Probe> = menu.iter().collect();
    canonical_menu.sort_by(|left, right| left.name.cmp(&right.name));
    push_u32(&mut canonical, declared_probe_count, "probe menu")?;
    for probe in canonical_menu {
        push_text(&mut canonical, &probe.name, "probe name")?;
        push_text(&mut canonical, &probe.target, "probe target")?;
        canonical.extend_from_slice(&probe.cost.to_bits().to_le_bytes());
        canonical.extend_from_slice(&probe.shrink.to_bits().to_le_bytes());
        canonical.push(probe.kind.identity_tag());
    }
    Ok(hash_domain(domain, &canonical))
}

fn ranked_output_context(
    source_context_id: ContentHash,
    rows: &[RankedPurchase],
) -> Result<ContentHash, VoiError> {
    ranked_output_context_with_schema(
        RANKED_MENU_CONTEXT_DOMAIN,
        VOI_RANKED_MENU_IDENTITY_VERSION,
        source_context_id,
        rows,
    )
}

fn ranked_output_context_with_schema(
    domain: &str,
    producer_version: u32,
    source_context_id: ContentHash,
    rows: &[RankedPurchase],
) -> Result<ContentHash, VoiError> {
    ranked_output_context_with_declared_count(
        domain,
        producer_version,
        source_context_id,
        rows,
        rows.len(),
    )
}

fn ranked_output_context_with_declared_count(
    domain: &str,
    producer_version: u32,
    source_context_id: ContentHash,
    rows: &[RankedPurchase],
    declared_row_count: usize,
) -> Result<ContentHash, VoiError> {
    let mut canonical = Vec::new();
    canonical.extend_from_slice(&producer_version.to_le_bytes());
    canonical.extend_from_slice(source_context_id.as_bytes());
    push_u32(&mut canonical, declared_row_count, "ranked output rows")?;
    for row in rows {
        push_text(&mut canonical, &row.probe.name, "ranked probe name")?;
        canonical.extend_from_slice(&row.flip_before.to_bits().to_le_bytes());
        canonical.extend_from_slice(&row.flip_after.to_bits().to_le_bytes());
        canonical.extend_from_slice(&row.score.to_bits().to_le_bytes());
    }
    Ok(hash_domain(domain, &canonical))
}

#[derive(Debug, Clone, Copy)]
struct PreparedProbe {
    node_idx: usize,
    post_lo: f64,
    post_hi: f64,
}

fn prepare_probes(
    nodes: &[UncertaintyNode],
    menu: &[Probe],
    targets: &[usize],
) -> Result<Vec<PreparedProbe>, VoiError> {
    let mut prepared = Vec::with_capacity(menu.len());
    for (probe, &node_idx) in menu.iter().zip(targets) {
        let node = node_at(nodes, node_idx)?;
        let contracted_left = (node.nominal - node.lo) * probe.shrink;
        let contracted_right = (node.hi - node.nominal) * probe.shrink;
        let post_lo = node.nominal - contracted_left;
        let post_hi = node.nominal + contracted_right;
        let post_width = post_hi - post_lo;
        let expected_width = (node.hi - node.lo) * probe.shrink;
        if !contracted_left.is_finite()
            || !contracted_right.is_finite()
            || !post_lo.is_finite()
            || !post_hi.is_finite()
            || !post_width.is_finite()
            || !expected_width.is_finite()
            || (node.nominal > node.lo && contracted_left == 0.0)
            || (node.hi > node.nominal && contracted_right == 0.0)
            || post_lo < node.lo
            || post_lo > node.nominal
            || post_hi < node.nominal
            || post_hi > node.hi
            || (node.hi > node.lo && post_width == 0.0)
        {
            return Err(VoiError::ArithmeticRefusal {
                operation: "post-probe interval contraction",
                subject: probe.name.clone(),
            });
        }
        prepared.push(PreparedProbe {
            node_idx,
            post_lo,
            post_hi,
        });
    }
    Ok(prepared)
}

/// Rank the probe menu by sampled flip-fraction reduction per dollar for the live
/// decision — MYOPIC one-step VoI (each probe is evaluated against the
/// CURRENT state only; no sequential tree).
///
/// The complete evaluation count and declared-work charge are admitted before
/// the first callback. Calls then receive canonical ordinal permits and are
/// bracketed by [`Cx`] checkpoints. A refusal returns no [`RankedMenu`]; oracle
/// implementations remain responsible for internal checkpoints and truthful
/// declared work.
///
/// # Errors
/// [`VoiError`] when the decision, node set, menu, targets, grid, probe
/// economics, callback margins, or derived arithmetic are invalid.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)] // explicit query/proof pipeline
pub fn rank_purchases(
    cx: &Cx,
    decision: &dyn DecisionOracle,
    nodes: &[UncertaintyNode],
    menu: &[Probe],
    grid: usize,
    budget: DecisionBudget,
    policy_scope: &str,
    snapshot_id: &str,
) -> Result<RankedMenu, VoiError> {
    validate_name("VoI policy scope", 0, policy_scope)?;
    validate_name("VoI snapshot identity", 0, snapshot_id)?;
    let metadata = decision_oracle_metadata(decision)?;
    validate_nodes(metadata.arity, nodes)?;
    validate_grid(grid)?;
    let targets = validate_menu(nodes, menu)?;
    let unique_targets = targets.iter().copied().collect::<BTreeSet<_>>();
    let sweep_count = unique_targets
        .len()
        .checked_add(menu.len())
        .ok_or_else(|| VoiError::ArithmeticRefusal {
            operation: "ranking sweep count",
            subject: "probe menu".to_string(),
        })?;
    let evaluations = grid
        .checked_mul(sweep_count)
        .and_then(|sweeps| sweeps.checked_add(1))
        .ok_or_else(|| VoiError::ArithmeticRefusal {
            operation: "ranking evaluation count",
            subject: "probe menu".to_string(),
        })?;
    let computation = admit_decision_computation(metadata, evaluations, budget)?;
    // All input-derived intervals are prepared before the first callback, so
    // a malformed later probe cannot leave observable partial callback work.
    let prepared = prepare_probes(nodes, menu, &targets)?;
    let source_context_id = ranked_source_context(
        nodes,
        menu,
        grid,
        policy_scope,
        snapshot_id,
        metadata,
        computation,
    )?;
    let mut meter = DecisionComputationMeter::new(metadata, computation);
    let base = nominal_verdict_validated(cx, decision, &mut meter, nodes)?;
    let mut flip_before = vec![None; nodes.len()];
    for node_idx in unique_targets {
        let node = node_at(nodes, node_idx)?;
        flip_before[node_idx] = Some(flip_probability_validated(
            cx, decision, &mut meter, nodes, base, node_idx, node.lo, node.hi, grid,
        )?);
    }

    let mut ranked = Vec::with_capacity(menu.len());
    let mut evaluation_order = (0..menu.len()).collect::<Vec<_>>();
    evaluation_order.sort_by(|left, right| menu[*left].name.cmp(&menu[*right].name));
    for index in evaluation_order {
        let probe = &menu[index];
        let prepared = prepared[index];
        let before = flip_before[prepared.node_idx].ok_or_else(|| VoiError::ArithmeticRefusal {
            operation: "pre-probe sweep lookup",
            subject: probe.name.clone(),
        })?;
        let flip_after = flip_probability_validated(
            cx,
            decision,
            &mut meter,
            nodes,
            base,
            prepared.node_idx,
            prepared.post_lo,
            prepared.post_hi,
            grid,
        )?;
        let score = (before - flip_after).max(0.0) / probe.cost;
        if !score.is_finite() || score < 0.0 {
            return Err(VoiError::ArithmeticRefusal {
                operation: "sampled flip-fraction-per-dollar score",
                subject: probe.name.clone(),
            });
        }
        ranked.push(RankedPurchase {
            probe: probe.clone(),
            flip_before: before,
            flip_after,
            score,
        });
    }
    meter.finish()?;
    ranked.sort_by(compare_ranked);
    let context_id = ranked_output_context(source_context_id, &ranked)?;
    Ok(RankedMenu {
        rows: ranked,
        source_context_id,
        context_id,
        policy_scope: policy_scope.to_string(),
        snapshot_id: snapshot_id.to_string(),
        grid,
        computation,
    })
}

/// Surface a structured QUERY-RESULT HINT. Every scalar is explicitly a
/// grid-sampled estimate; a sampled zero is never rendered as proof that no
/// evidence could change the decision.
#[must_use]
pub fn hint_for_query(ranked: &RankedMenu) -> QueryHint {
    QueryHint {
        context_id: ranked.context_id,
        grid: ranked.grid,
        purchase: ranked.rows.iter().find(|row| row.score > 0.0).cloned(),
    }
}

/// The audit verdict for reporting. This enum is not scheduling authority;
/// only the live [`VoiScheduler`] owns executable state.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AuditVerdict {
    /// Anytime-valid evidence crossed the fixed activation threshold.
    KeepScheduling,
    /// Evidence is absent, insufficient, or has not crossed the threshold.
    DemoteToReporting,
}

/// One validated matched-cost prospective-audit observation.
///
/// Fields are private so raw booleans and unmatched prices cannot enter the
/// e-process without identity, provenance, and economic validation.
#[derive(Debug, Clone, PartialEq)]
pub struct MatchedAuditRecord {
    observation_id: String,
    recommended_id: String,
    alternative_id: String,
    provenance: String,
    matched_cost: f64,
    recommended_changed_decision: bool,
    alternative_changed_decision: bool,
}

impl MatchedAuditRecord {
    /// Construct one matched-cost comparison.
    ///
    /// # Errors
    /// [`VoiError`] unless identities/provenance are bounded visible ASCII,
    /// candidates differ, and both finite positive costs are bit-identical.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        observation_id: impl Into<String>,
        recommended_id: impl Into<String>,
        alternative_id: impl Into<String>,
        provenance: impl Into<String>,
        recommended_cost: f64,
        alternative_cost: f64,
        recommended_changed_decision: bool,
        alternative_changed_decision: bool,
    ) -> Result<Self, VoiError> {
        let observation_id = observation_id.into();
        let recommended_id = recommended_id.into();
        let alternative_id = alternative_id.into();
        let provenance = provenance.into();
        for (kind, value) in [
            ("audit observation", observation_id.as_str()),
            ("recommended purchase", recommended_id.as_str()),
            ("alternative purchase", alternative_id.as_str()),
            ("audit provenance", provenance.as_str()),
        ] {
            validate_name(kind, 0, value)?;
        }
        if recommended_id == alternative_id {
            return Err(VoiError::InvalidAuditPair {
                observation: observation_id,
            });
        }
        if !recommended_cost.is_finite()
            || recommended_cost <= 0.0
            || recommended_cost.to_bits() != alternative_cost.to_bits()
        {
            return Err(VoiError::InvalidAuditCost {
                observation: observation_id,
                recommended_cost,
                alternative_cost,
            });
        }
        // Caller-owned Strings can carry arbitrary spare capacity despite
        // bounded content. Rebuild every retained identity from its validated
        // slice so the audit-record cap is also an operational memory bound.
        let observation_id = observation_id.as_str().to_owned();
        let recommended_id = recommended_id.as_str().to_owned();
        let alternative_id = alternative_id.as_str().to_owned();
        let provenance = provenance.as_str().to_owned();
        Ok(Self {
            observation_id,
            recommended_id,
            alternative_id,
            provenance,
            matched_cost: recommended_cost,
            recommended_changed_decision,
            alternative_changed_decision,
        })
    }

    /// Stable observation identity used to prevent duplicate evidence.
    #[must_use]
    pub fn observation_id(&self) -> &str {
        &self.observation_id
    }

    /// Recommended-purchase identity.
    #[must_use]
    pub fn recommended_id(&self) -> &str {
        &self.recommended_id
    }

    /// Matched alternative-purchase identity.
    #[must_use]
    pub fn alternative_id(&self) -> &str {
        &self.alternative_id
    }

    /// Caller-supplied provenance identity.
    #[must_use]
    pub fn provenance(&self) -> &str {
        &self.provenance
    }

    /// Exact matched cost.
    #[must_use]
    pub fn matched_cost(&self) -> f64 {
        self.matched_cost
    }

    /// Whether the recommended purchase changed the realized decision.
    #[must_use]
    pub fn recommended_changed_decision(&self) -> bool {
        self.recommended_changed_decision
    }

    /// Whether the matched alternative changed the realized decision.
    #[must_use]
    pub fn alternative_changed_decision(&self) -> bool {
        self.alternative_changed_decision
    }
}

#[allow(dead_code)]
fn classify_voi_audit_context_identity_fields(record: &MatchedAuditRecord) {
    let MatchedAuditRecord {
        observation_id: _,
        recommended_id: _,
        alternative_id: _,
        provenance: _,
        matched_cost: _,
        recommended_changed_decision: _,
        alternative_changed_decision: _,
    } = record;
}

/// One caller-scoped, single-epoch scheduling suggestion. It retains the ranked
/// snapshot, audit evidence root, and exact budget transition, but remains
/// caller-supplied and is not ledger authority. Use
/// [`VoiScheduler::schedule_authenticated`] for an authenticated result.
#[derive(Debug, PartialEq)]
pub struct ScheduledPurchase {
    purchase: RankedPurchase,
    ranked_context_id: ContentHash,
    ranked_source_context_id: ContentHash,
    policy_scope: String,
    snapshot_id: String,
    ranked_grid: usize,
    audit_context_id: ContentHash,
    audit_observations: usize,
    audit_log_e_value: f64,
    budget_dollars: f64,
    remaining_budget_dollars: f64,
}

impl ScheduledPurchase {
    /// Suggested purchase (authenticated only when wrapped by
    /// [`AuthenticatedScheduledPurchase`]).
    #[must_use]
    pub fn purchase(&self) -> &RankedPurchase {
        &self.purchase
    }

    /// Ranked node/menu/grid snapshot identity.
    #[must_use]
    pub fn ranked_context_id(&self) -> ContentHash {
        self.ranked_context_id
    }

    /// Source policy/snapshot/node/menu/grid identity before evaluation.
    #[must_use]
    pub fn ranked_source_context_id(&self) -> ContentHash {
        self.ranked_source_context_id
    }

    /// Scheduling policy/version scope shared by the audit and ranking.
    #[must_use]
    pub fn policy_scope(&self) -> &str {
        &self.policy_scope
    }

    /// Decision/ledger snapshot whose ranking was consumed.
    #[must_use]
    pub fn snapshot_id(&self) -> &str {
        &self.snapshot_id
    }

    /// Midpoint grid supporting the sampled purchase score.
    #[must_use]
    pub fn ranked_grid(&self) -> usize {
        self.ranked_grid
    }

    /// Anytime-valid matched-audit evidence identity.
    #[must_use]
    pub fn audit_context_id(&self) -> ContentHash {
        self.audit_context_id
    }

    /// Matched-cost observation count supporting the audit verdict.
    #[must_use]
    pub fn audit_observations(&self) -> usize {
        self.audit_observations
    }

    /// Final log e-value supporting the audit verdict.
    #[must_use]
    pub fn audit_log_e_value(&self) -> f64 {
        self.audit_log_e_value
    }

    /// Admitted scheduling budget in dollars.
    #[must_use]
    pub fn budget_dollars(&self) -> f64 {
        self.budget_dollars
    }

    /// Exact remaining budget in dollars after this one purchase.
    #[must_use]
    pub fn remaining_budget_dollars(&self) -> f64 {
        self.remaining_budget_dollars
    }
}

/// A scheduling result whose ranked menu and every prospective audit outcome
/// were independently reverified against immutable ledger receipts.
#[derive(Debug)]
pub struct AuthenticatedScheduledPurchase {
    scheduled: ScheduledPurchase,
    ranked_receipt: VoiLedgerReceipt,
    audit_receipts: Vec<VoiLedgerReceipt>,
    valid_through_day: u32,
}

impl AuthenticatedScheduledPurchase {
    /// The fully retained scheduling decision and budget transition.
    #[must_use]
    pub fn scheduled(&self) -> &ScheduledPurchase {
        &self.scheduled
    }

    /// Receipt authenticating the complete ranked menu and current snapshot.
    #[must_use]
    pub const fn ranked_receipt(&self) -> VoiLedgerReceipt {
        self.ranked_receipt
    }

    /// Prospectively ordered receipts authenticating every audit observation.
    #[must_use]
    pub fn audit_receipts(&self) -> &[VoiLedgerReceipt] {
        &self.audit_receipts
    }

    /// Inclusive minimum expiry across ranking and audit receipts.
    #[must_use]
    pub const fn valid_through_day(&self) -> u32 {
        self.valid_through_day
    }
}

/// Immutable reporting snapshot of one live prospective-audit prefix. This is
/// never scheduling authority; only the owning [`VoiScheduler`] can spend.
#[derive(Debug, Clone, PartialEq)]
pub struct AuditReport {
    policy_scope: String,
    audit_context_id: ContentHash,
    observations: usize,
    log_e_value: f64,
    verdict: AuditVerdict,
}

impl AuditReport {
    /// Reporting verdict.
    #[must_use]
    pub fn verdict(&self) -> AuditVerdict {
        self.verdict
    }

    /// Caller-declared audit policy/version scope.
    #[must_use]
    pub fn policy_scope(&self) -> &str {
        &self.policy_scope
    }

    /// Content identity of the canonical evidence prefix.
    #[must_use]
    pub fn audit_context_id(&self) -> ContentHash {
        self.audit_context_id
    }

    /// Number of matched-cost observations evaluated.
    #[must_use]
    pub fn observations(&self) -> usize {
        self.observations
    }

    /// Final log e-value, useful for reporting progress before activation.
    #[must_use]
    pub fn log_e_value(&self) -> f64 {
        self.log_e_value
    }

    /// Audit-context identity semantics used for [`Self::audit_context_id`].
    #[must_use]
    pub const fn identity_version(&self) -> u32 {
        VOI_AUDIT_CONTEXT_IDENTITY_VERSION
    }

    /// Admit a version carried beside a retained audit-context root.
    ///
    /// Root comparison must follow this fail-closed version check.
    ///
    /// # Errors
    /// [`VoiError::UnsupportedIdentityVersion`] for a stale or future audit
    /// identity version.
    pub fn admit_retained_identity_version(&self, declared: u32) -> Result<(), VoiError> {
        check_identity_version("VoI audit context", declared, self.identity_version())
    }
}

fn audit_context(
    policy_scope: &str,
    records: &[MatchedAuditRecord],
) -> Result<ContentHash, VoiError> {
    audit_context_with_schema(
        AUDIT_CONTEXT_DOMAIN,
        VOI_AUDIT_CONTEXT_IDENTITY_VERSION,
        policy_scope,
        records,
        VOI_AUDIT_ALPHA,
        MAX_VOI_AUDIT_RECORDS,
    )
}

fn audit_context_with_schema(
    domain: &str,
    producer_version: u32,
    policy_scope: &str,
    records: &[MatchedAuditRecord],
    audit_alpha: f64,
    max_audit_records: usize,
) -> Result<ContentHash, VoiError> {
    audit_context_with_declared_count(
        domain,
        producer_version,
        policy_scope,
        records,
        records.len(),
        audit_alpha,
        max_audit_records,
    )
}

fn audit_context_with_declared_count(
    domain: &str,
    producer_version: u32,
    policy_scope: &str,
    records: &[MatchedAuditRecord],
    declared_record_count: usize,
    audit_alpha: f64,
    max_audit_records: usize,
) -> Result<ContentHash, VoiError> {
    let mut canonical = Vec::new();
    canonical.extend_from_slice(&producer_version.to_le_bytes());
    push_text(&mut canonical, policy_scope, "VoI audit policy scope")?;
    canonical.extend_from_slice(&audit_alpha.to_bits().to_le_bytes());
    push_u32(&mut canonical, max_audit_records, "maximum audit records")?;
    push_u32(&mut canonical, declared_record_count, "audit records")?;
    for record in records {
        push_text(&mut canonical, &record.observation_id, "audit observation")?;
        push_text(
            &mut canonical,
            &record.recommended_id,
            "recommended purchase",
        )?;
        push_text(
            &mut canonical,
            &record.alternative_id,
            "alternative purchase",
        )?;
        push_text(&mut canonical, &record.provenance, "audit provenance")?;
        canonical.extend_from_slice(&record.matched_cost.to_bits().to_le_bytes());
        canonical.push(u8::from(record.recommended_changed_decision));
        canonical.push(u8::from(record.alternative_changed_decision));
    }
    Ok(hash_domain(domain, &canonical))
}

fn ledger_evidence_error(operation: &'static str, error: fs_ledger::LedgerError) -> VoiError {
    VoiError::LedgerEvidence {
        operation,
        detail: error.to_string(),
    }
}

fn verify_receipt_day(context: &VoiLedgerContext, current_day: u32) -> Result<(), VoiError> {
    if current_day < context.issued_day {
        return Err(VoiError::ReceiptNotYetValid {
            issued_day: context.issued_day,
            current_day,
        });
    }
    if current_day > context.expires_day {
        return Err(VoiError::ReceiptExpired {
            expires_day: context.expires_day,
            current_day,
        });
    }
    Ok(())
}

fn push_receipt_context(
    canonical: &mut Vec<u8>,
    ledger_id: ContentHash,
    context: &VoiLedgerContext,
) -> Result<(), VoiError> {
    canonical.extend_from_slice(&VOI_LEDGER_RECEIPT_VERSION.to_le_bytes());
    canonical.extend_from_slice(ledger_id.as_bytes());
    canonical.extend_from_slice(&context.branch.to_le_bytes());
    canonical.extend_from_slice(context.snapshot_root.as_bytes());
    push_text(canonical, &context.session_id, "receipt session")?;
    push_text(
        canonical,
        &context.decision_version,
        "receipt decision version",
    )?;
    push_text(canonical, &context.model_version, "receipt model version")?;
    push_text(canonical, &context.policy_version, "receipt policy version")?;
    canonical.extend_from_slice(&context.issued_day.to_le_bytes());
    canonical.extend_from_slice(&context.expires_day.to_le_bytes());
    Ok(())
}

fn validate_receipt_len(canonical: Vec<u8>) -> Result<Vec<u8>, VoiError> {
    if canonical.is_empty() || canonical.len() > MAX_VOI_LEDGER_RECEIPT_BYTES {
        return Err(VoiError::SizeLimit {
            collection: "VoI ledger receipt bytes",
            count: canonical.len(),
            min: 1,
            max: MAX_VOI_LEDGER_RECEIPT_BYTES,
        });
    }
    Ok(canonical)
}

/// Exact metadata for a decision identity artifact.
///
/// # Errors
/// [`VoiError::InvalidReceiptContext`] when `decision_version` is not a
/// bounded receipt token.
pub fn voi_decision_artifact_metadata(decision_version: &str) -> Result<String, VoiError> {
    validate_receipt_token("decision_version", decision_version)?;
    Ok(format!(
        "{{\"schema\":\"fs-plan-voi-decision-v1\",\"version\":\"{}\"}}",
        decision_version
    ))
}

/// Exact metadata for a model identity artifact.
///
/// # Errors
/// [`VoiError::InvalidReceiptContext`] when `model_version` is not a bounded
/// receipt token.
pub fn voi_model_artifact_metadata(model_version: &str) -> Result<String, VoiError> {
    validate_receipt_token("model_version", model_version)?;
    Ok(format!(
        "{{\"schema\":\"fs-plan-voi-model-v1\",\"version\":\"{}\"}}",
        model_version
    ))
}

fn verify_snapshot(
    ledger: &Ledger,
    vcs: &Vcs,
    context: &VoiLedgerContext,
    current_day: u32,
    require_current_head: bool,
) -> Result<(ContentHash, fs_ledger::ViewSnapshot), VoiError> {
    verify_receipt_day(context, current_day)?;
    ledger
        .checked_instance_id()
        .map_err(|error| ledger_evidence_error("check physical ledger identity", error))?;
    let ledger_id = ledger
        .vcs_identity()
        .map_err(|error| ledger_evidence_error("read VCS ledger identity", error))?;
    if require_current_head {
        let head = vcs
            .head(ledger, context.branch)
            .map_err(|error| ledger_evidence_error("read VCS branch head", error))?;
        if head != Some(context.snapshot_root) {
            return Err(VoiError::StaleLedgerSnapshot {
                expected: context.snapshot_root,
                actual: head,
            });
        }
    }
    let snapshot = vcs
        .checkout(ledger, context.branch, &context.snapshot_root)
        .map_err(|error| ledger_evidence_error("checkout bound VCS snapshot", error))?;
    Ok((ledger_id, snapshot))
}

fn require_snapshot_artifact(
    snapshot: &fs_ledger::ViewSnapshot,
    role: &'static str,
    artifact: ContentHash,
) -> Result<(), VoiError> {
    if snapshot.artifacts.contains(&artifact) {
        Ok(())
    } else {
        Err(VoiError::SnapshotArtifactMissing { role, artifact })
    }
}

fn verify_artifact_envelope(
    ledger: &Ledger,
    artifact: ContentHash,
    expected_kind: &str,
    expected_meta: &str,
    expected_bytes: Option<&[u8]>,
    max_bytes: u64,
) -> Result<(), VoiError> {
    let info = ledger
        .artifact_info(&artifact)
        .map_err(|error| ledger_evidence_error("read receipt input artifact", error))?
        .ok_or(VoiError::ReceiptMismatch {
            field: "artifact presence",
        })?;
    if info.kind != expected_kind || info.meta.as_deref() != Some(expected_meta) {
        return Err(VoiError::ReceiptMismatch {
            field: "artifact envelope",
        });
    }
    if info.len == 0 {
        return Err(VoiError::ReceiptMismatch {
            field: "nonempty artifact bytes",
        });
    }
    let stored = ledger
        .get_artifact_bounded(&artifact, max_bytes)
        .map_err(|error| ledger_evidence_error("read bounded receipt input", error))?
        .ok_or(VoiError::ReceiptMismatch {
            field: "artifact bytes",
        })?;
    if u64::try_from(stored.len()).ok() != Some(info.len) || hash_bytes(&stored) != artifact {
        return Err(VoiError::ReceiptMismatch {
            field: "artifact content identity",
        });
    }
    if let Some(expected_bytes) = expected_bytes
        && stored.as_slice() != expected_bytes
    {
        return Err(VoiError::ReceiptMismatch {
            field: "artifact bytes",
        });
    }
    Ok(())
}

fn ranked_catalog_bytes(ranked: &RankedMenu) -> Result<Vec<u8>, VoiError> {
    let probes = ranked
        .rows
        .iter()
        .map(|row| row.probe.clone())
        .collect::<Vec<_>>();
    voi_probe_catalog_bytes(&probes)
}

fn validate_ranked_receipt_inputs(
    ledger: &Ledger,
    snapshot: &fs_ledger::ViewSnapshot,
    ranked: &RankedMenu,
    context: &VoiLedgerContext,
    decision_artifact: ContentHash,
    model_artifact: ContentHash,
    probe_catalog_artifact: ContentHash,
) -> Result<Vec<u8>, VoiError> {
    if ranked.policy_scope != context.policy_version {
        return Err(VoiError::ReceiptMismatch {
            field: "ranked policy version",
        });
    }
    if ranked.snapshot_id != context.snapshot_root.to_hex() {
        return Err(VoiError::ReceiptMismatch {
            field: "ranked snapshot root",
        });
    }
    let catalog = ranked_catalog_bytes(ranked)?;
    for (role, artifact) in [
        ("decision", decision_artifact),
        ("model", model_artifact),
        ("probe catalog", probe_catalog_artifact),
    ] {
        require_snapshot_artifact(snapshot, role, artifact)?;
    }
    verify_artifact_envelope(
        ledger,
        decision_artifact,
        VOI_DECISION_ARTIFACT_KIND,
        &voi_decision_artifact_metadata(&context.decision_version)?,
        None,
        u64::try_from(MAX_VOI_AUTHORITY_ARTIFACT_BYTES).unwrap_or(u64::MAX),
    )?;
    verify_artifact_envelope(
        ledger,
        model_artifact,
        VOI_MODEL_ARTIFACT_KIND,
        &voi_model_artifact_metadata(&context.model_version)?,
        None,
        u64::try_from(MAX_VOI_AUTHORITY_ARTIFACT_BYTES).unwrap_or(u64::MAX),
    )?;
    verify_artifact_envelope(
        ledger,
        probe_catalog_artifact,
        VOI_PROBE_CATALOG_ARTIFACT_KIND,
        VOI_PROBE_CATALOG_ARTIFACT_METADATA,
        Some(&catalog),
        u64::try_from(MAX_VOI_PROBE_CATALOG_BYTES).unwrap_or(u64::MAX),
    )?;
    Ok(catalog)
}

#[allow(clippy::too_many_arguments)]
fn ranked_receipt_bytes(
    ledger_id: ContentHash,
    ranked: &RankedMenu,
    context: &VoiLedgerContext,
    decision_artifact: ContentHash,
    model_artifact: ContentHash,
    probe_catalog_artifact: ContentHash,
) -> Result<Vec<u8>, VoiError> {
    let mut canonical = Vec::new();
    canonical.extend_from_slice(RANKED_LEDGER_RECEIPT_MAGIC);
    push_receipt_context(&mut canonical, ledger_id, context)?;
    canonical.extend_from_slice(&ranked.source_identity_version().to_le_bytes());
    canonical.extend_from_slice(&ranked.identity_version().to_le_bytes());
    canonical.extend_from_slice(ranked.source_context_id.as_bytes());
    canonical.extend_from_slice(ranked.context_id.as_bytes());
    push_u32(&mut canonical, ranked.rows.len(), "ranked receipt rows")?;
    canonical.extend_from_slice(decision_artifact.as_bytes());
    canonical.extend_from_slice(model_artifact.as_bytes());
    canonical.extend_from_slice(probe_catalog_artifact.as_bytes());
    validate_receipt_len(canonical)
}

/// Canonical bytes for one side of a prospective matched-cost outcome.
///
/// Retaining these bytes as an `fs-plan-voi-outcome` artifact lets receipt
/// verification authenticate the exact candidate, cost, provenance, and
/// realized decision-change bit instead of trusting a caller boolean.
///
/// # Errors
/// [`VoiError`] if a bounded field cannot be encoded.
pub fn voi_audit_outcome_bytes(
    record: &MatchedAuditRecord,
    recommended: bool,
) -> Result<Vec<u8>, VoiError> {
    let mut canonical = Vec::new();
    canonical.extend_from_slice(AUDIT_OUTCOME_MAGIC);
    canonical.extend_from_slice(&VOI_LEDGER_RECEIPT_VERSION.to_le_bytes());
    push_text(&mut canonical, &record.observation_id, "audit observation")?;
    push_text(
        &mut canonical,
        if recommended {
            &record.recommended_id
        } else {
            &record.alternative_id
        },
        "audit candidate",
    )?;
    push_text(&mut canonical, &record.provenance, "audit provenance")?;
    canonical.extend_from_slice(&record.matched_cost.to_bits().to_le_bytes());
    canonical.push(u8::from(if recommended {
        record.recommended_changed_decision
    } else {
        record.alternative_changed_decision
    }));
    validate_receipt_len(canonical)
}

#[allow(clippy::too_many_arguments)]
fn audit_receipt_bytes(
    ledger_id: ContentHash,
    record: &MatchedAuditRecord,
    sequence: usize,
    context: &VoiLedgerContext,
    decision_artifact: ContentHash,
    model_artifact: ContentHash,
    recommended_outcome_artifact: ContentHash,
    alternative_outcome_artifact: ContentHash,
) -> Result<Vec<u8>, VoiError> {
    let mut canonical = Vec::new();
    canonical.extend_from_slice(AUDIT_LEDGER_RECEIPT_MAGIC);
    push_receipt_context(&mut canonical, ledger_id, context)?;
    canonical.extend_from_slice(&VOI_AUDIT_CONTEXT_IDENTITY_VERSION.to_le_bytes());
    push_u32(&mut canonical, sequence, "audit receipt sequence")?;
    push_text(&mut canonical, &record.observation_id, "audit observation")?;
    push_text(
        &mut canonical,
        &record.recommended_id,
        "recommended purchase",
    )?;
    push_text(
        &mut canonical,
        &record.alternative_id,
        "alternative purchase",
    )?;
    push_text(&mut canonical, &record.provenance, "audit provenance")?;
    canonical.extend_from_slice(&record.matched_cost.to_bits().to_le_bytes());
    canonical.push(u8::from(record.recommended_changed_decision));
    canonical.push(u8::from(record.alternative_changed_decision));
    canonical.extend_from_slice(decision_artifact.as_bytes());
    canonical.extend_from_slice(model_artifact.as_bytes());
    canonical.extend_from_slice(recommended_outcome_artifact.as_bytes());
    canonical.extend_from_slice(alternative_outcome_artifact.as_bytes());
    validate_receipt_len(canonical)
}

fn receipt_ir(kind: &str, artifact: ContentHash, snapshot: ContentHash) -> String {
    format!(
        "{{\"op\":\"voi.receipt\",\"kind\":\"{kind}\",\"artifact\":\"{artifact}\",\"snapshot\":\"{snapshot}\"}}"
    )
}

fn persist_receipt(
    ledger: &Ledger,
    context: &VoiLedgerContext,
    kind: &'static str,
    meta: &'static str,
    canonical: &[u8],
    inputs: &[ContentHash],
) -> Result<VoiLedgerReceipt, VoiError> {
    if ledger.in_transaction() {
        return Err(VoiError::ReceiptTransactionActive);
    }
    let artifact = hash_bytes(canonical);
    if let Some(operation) = ledger
        .artifact_output_seal(&artifact)
        .map_err(|error| ledger_evidence_error("look up idempotent receipt", error))?
    {
        let receipt = VoiLedgerReceipt {
            operation,
            artifact,
        };
        verify_receipt_envelope(ledger, context, receipt, kind, meta, canonical, inputs)?;
        return Ok(receipt);
    }
    let ir = receipt_ir(kind, artifact, context.snapshot_root);
    let issued_at_ns = receipt_time_ns(context.issued_day)?;
    let explicits = FiveExplicits {
        seed: RECEIPT_SEED,
        versions: RECEIPT_VERSIONS,
        budget: RECEIPT_BUDGET,
        capability: RECEIPT_CAPABILITY,
    };
    ledger
        .begin()
        .map_err(|error| ledger_evidence_error("begin receipt transaction", error))?;
    let result = (|| {
        let operation = ledger
            .begin_op_on(
                context.branch,
                ExecMode::Deterministic,
                Some(context.session_id.as_bytes()),
                &ir,
                &explicits,
                issued_at_ns,
            )
            .map_err(|error| ledger_evidence_error("begin receipt operation", error))?;
        let mut unique_inputs = BTreeSet::new();
        for input in inputs {
            if unique_inputs.insert(*input) {
                ledger
                    .link(operation, input, EdgeRole::In)
                    .map_err(|error| ledger_evidence_error("link receipt input", error))?;
            }
        }
        let stored = ledger
            .put_artifact(kind, canonical, Some(meta))
            .map_err(|error| ledger_evidence_error("store receipt artifact", error))?;
        if stored.hash != artifact
            || stored.len != u64::try_from(canonical.len()).unwrap_or(u64::MAX)
        {
            return Err(VoiError::ReceiptMismatch {
                field: "stored receipt",
            });
        }
        ledger
            .link(operation, &artifact, EdgeRole::Out)
            .map_err(|error| ledger_evidence_error("link receipt output", error))?;
        ledger
            .seal_artifact_output(&artifact, operation)
            .map_err(|error| ledger_evidence_error("seal receipt producer", error))?;
        ledger
            .finish_op(operation, OpOutcome::Ok, None, issued_at_ns)
            .map_err(|error| ledger_evidence_error("finish receipt operation", error))?;
        Ok(VoiLedgerReceipt {
            operation,
            artifact,
        })
    })();
    match result {
        Ok(receipt) => {
            if let Err(error) = ledger.commit() {
                let _ = ledger.rollback();
                return Err(ledger_evidence_error("commit receipt transaction", error));
            }
            Ok(receipt)
        }
        Err(error) => {
            let _ = ledger.rollback();
            Err(error)
        }
    }
}

fn verify_receipt_envelope(
    ledger: &Ledger,
    context: &VoiLedgerContext,
    receipt: VoiLedgerReceipt,
    kind: &'static str,
    meta: &'static str,
    canonical: &[u8],
    inputs: &[ContentHash],
) -> Result<(), VoiError> {
    let expected_artifact = hash_bytes(canonical);
    if receipt.operation <= 0 || receipt.artifact != expected_artifact {
        return Err(VoiError::ReceiptMismatch {
            field: "receipt locator",
        });
    }
    verify_artifact_envelope(
        ledger,
        receipt.artifact,
        kind,
        meta,
        Some(canonical),
        u64::try_from(MAX_VOI_LEDGER_RECEIPT_BYTES).unwrap_or(u64::MAX),
    )?;
    let operation = ledger
        .op(receipt.operation)
        .map_err(|error| ledger_evidence_error("read receipt operation", error))?
        .ok_or(VoiError::ReceiptMismatch {
            field: "receipt operation",
        })?;
    let expected_ir = receipt_ir(kind, receipt.artifact, context.snapshot_root);
    let issued_at_ns = receipt_time_ns(context.issued_day)?;
    if operation.session.as_deref() != Some(context.session_id.as_bytes())
        || operation.ir != expected_ir
        || operation.seed.as_slice() != RECEIPT_SEED
        || operation.versions != RECEIPT_VERSIONS
        || operation.budget != RECEIPT_BUDGET
        || operation.capability != RECEIPT_CAPABILITY
        || operation.t_start != issued_at_ns
        || operation.t_end != Some(issued_at_ns)
        || operation.outcome.as_deref() != Some("ok")
        || operation.diag.is_some()
    {
        return Err(VoiError::ReceiptMismatch {
            field: "receipt operation envelope",
        });
    }
    let execution = ledger
        .op_execution_context(receipt.operation)
        .map_err(|error| ledger_evidence_error("read receipt execution context", error))?
        .ok_or(VoiError::ReceiptMismatch {
            field: "receipt execution context",
        })?;
    if execution.branch != context.branch || execution.exec_mode != ExecMode::Deterministic {
        return Err(VoiError::ReceiptMismatch {
            field: "receipt execution context",
        });
    }
    let unique_inputs = inputs.iter().copied().collect::<BTreeSet<_>>();
    for input in &unique_inputs {
        if !ledger
            .edge_exists(receipt.operation, input, EdgeRole::In)
            .map_err(|error| ledger_evidence_error("verify receipt input edge", error))?
        {
            return Err(VoiError::ReceiptMismatch {
                field: "receipt input lineage",
            });
        }
    }
    if !ledger
        .edge_exists(receipt.operation, &receipt.artifact, EdgeRole::Out)
        .map_err(|error| ledger_evidence_error("verify receipt output edge", error))?
        || ledger
            .artifact_output_seal(&receipt.artifact)
            .map_err(|error| ledger_evidence_error("verify receipt producer seal", error))?
            != Some(receipt.operation)
    {
        return Err(VoiError::ReceiptMismatch {
            field: "receipt output lineage",
        });
    }
    let mut expected_edges = unique_inputs
        .into_iter()
        .map(|artifact| OpArtifactEdge {
            role: EdgeRole::In,
            artifact,
        })
        .collect::<Vec<_>>();
    expected_edges.push(OpArtifactEdge {
        role: EdgeRole::Out,
        artifact: receipt.artifact,
    });
    let actual_edges = ledger
        .op_artifact_edges_bounded(receipt.operation, expected_edges.len())
        .map_err(|error| ledger_evidence_error("verify exact receipt lineage", error))?;
    if actual_edges.truncated || actual_edges.edges != expected_edges {
        return Err(VoiError::ReceiptMismatch {
            field: "exact receipt lineage",
        });
    }
    Ok(())
}

/// A ranked menu whose ledger receipt was independently verified.
#[derive(Debug)]
pub struct AuthenticatedRankedMenu {
    ranked: RankedMenu,
    context: VoiLedgerContext,
    receipt: VoiLedgerReceipt,
    decision_artifact: ContentHash,
    model_artifact: ContentHash,
    probe_catalog_artifact: ContentHash,
}

impl AuthenticatedRankedMenu {
    /// Verified ranked menu for reporting before it is consumed.
    #[must_use]
    pub fn ranked(&self) -> &RankedMenu {
        &self.ranked
    }

    /// Verified receipt locator.
    #[must_use]
    pub const fn receipt(&self) -> VoiLedgerReceipt {
        self.receipt
    }

    /// Bound ledger authority scope.
    #[must_use]
    pub fn context(&self) -> &VoiLedgerContext {
        &self.context
    }

    /// Exact decision identity artifact admitted by the verifier.
    #[must_use]
    pub const fn decision_artifact(&self) -> ContentHash {
        self.decision_artifact
    }

    /// Exact model identity artifact admitted by the verifier.
    #[must_use]
    pub const fn model_artifact(&self) -> ContentHash {
        self.model_artifact
    }

    /// Complete canonical probe catalog admitted by the verifier.
    #[must_use]
    pub const fn probe_catalog_artifact(&self) -> ContentHash {
        self.probe_catalog_artifact
    }
}

/// Mint or exactly replay one ranked-menu receipt through an atomic,
/// deterministic ledger operation whose input lineage is the decision, model,
/// and complete catalog.
///
/// # Errors
/// [`VoiError`] when the context is stale/expired, the ranked menu is not bound
/// to the exact current root and policy, any input artifact is absent or has a
/// wrong envelope, or the atomic ledger write fails.
#[allow(clippy::too_many_arguments)]
pub fn record_ranked_menu_receipt(
    ledger: &Ledger,
    vcs: &Vcs,
    ranked: &RankedMenu,
    context: &VoiLedgerContext,
    decision_artifact: ContentHash,
    model_artifact: ContentHash,
    probe_catalog_artifact: ContentHash,
) -> Result<VoiLedgerReceipt, VoiError> {
    let (ledger_id, snapshot) = verify_snapshot(ledger, vcs, context, context.issued_day, true)?;
    let _catalog = validate_ranked_receipt_inputs(
        ledger,
        &snapshot,
        ranked,
        context,
        decision_artifact,
        model_artifact,
        probe_catalog_artifact,
    )?;
    let canonical = ranked_receipt_bytes(
        ledger_id,
        ranked,
        context,
        decision_artifact,
        model_artifact,
        probe_catalog_artifact,
    )?;
    persist_receipt(
        ledger,
        context,
        RANKED_RECEIPT_KIND,
        RANKED_RECEIPT_META,
        &canonical,
        &[decision_artifact, model_artifact, probe_catalog_artifact],
    )
}

/// Independently verify a ranked receipt and mint the opaque capability used
/// by authenticated scheduling.
///
/// # Errors
/// [`VoiError`] for any root/policy/catalog/artifact/operation/lineage/expiry
/// mismatch. The ranked menu is returned only inside the verified wrapper.
#[allow(clippy::too_many_arguments)]
pub fn verify_ranked_menu_receipt(
    ledger: &Ledger,
    vcs: &Vcs,
    ranked: RankedMenu,
    context: VoiLedgerContext,
    decision_artifact: ContentHash,
    model_artifact: ContentHash,
    probe_catalog_artifact: ContentHash,
    receipt: VoiLedgerReceipt,
    current_day: u32,
) -> Result<AuthenticatedRankedMenu, VoiError> {
    let (ledger_id, snapshot) = verify_snapshot(ledger, vcs, &context, current_day, true)?;
    let _catalog = validate_ranked_receipt_inputs(
        ledger,
        &snapshot,
        &ranked,
        &context,
        decision_artifact,
        model_artifact,
        probe_catalog_artifact,
    )?;
    let canonical = ranked_receipt_bytes(
        ledger_id,
        &ranked,
        &context,
        decision_artifact,
        model_artifact,
        probe_catalog_artifact,
    )?;
    verify_receipt_envelope(
        ledger,
        &context,
        receipt,
        RANKED_RECEIPT_KIND,
        RANKED_RECEIPT_META,
        &canonical,
        &[decision_artifact, model_artifact, probe_catalog_artifact],
    )?;
    Ok(AuthenticatedRankedMenu {
        ranked,
        context,
        receipt,
        decision_artifact,
        model_artifact,
        probe_catalog_artifact,
    })
}

/// One audit observation whose outcome artifacts and receipt were verified.
#[derive(Debug, Clone)]
pub struct AuthenticatedAuditRecord {
    record: MatchedAuditRecord,
    sequence: usize,
    context: VoiLedgerContext,
    receipt: VoiLedgerReceipt,
    decision_artifact: ContentHash,
    model_artifact: ContentHash,
    recommended_outcome_artifact: ContentHash,
    alternative_outcome_artifact: ContentHash,
}

impl AuthenticatedAuditRecord {
    /// Exact prospectively ordered record.
    #[must_use]
    pub fn record(&self) -> &MatchedAuditRecord {
        &self.record
    }

    /// Zero-based prospective sequence bound into the receipt.
    #[must_use]
    pub const fn sequence(&self) -> usize {
        self.sequence
    }

    /// Verified receipt locator.
    #[must_use]
    pub const fn receipt(&self) -> VoiLedgerReceipt {
        self.receipt
    }

    /// Bound ledger authority scope.
    #[must_use]
    pub fn context(&self) -> &VoiLedgerContext {
        &self.context
    }

    /// Exact decision identity artifact admitted by the verifier.
    #[must_use]
    pub const fn decision_artifact(&self) -> ContentHash {
        self.decision_artifact
    }

    /// Exact model identity artifact admitted by the verifier.
    #[must_use]
    pub const fn model_artifact(&self) -> ContentHash {
        self.model_artifact
    }

    /// Recommended-side outcome artifact admitted by the verifier.
    #[must_use]
    pub const fn recommended_outcome_artifact(&self) -> ContentHash {
        self.recommended_outcome_artifact
    }

    /// Alternative-side outcome artifact admitted by the verifier.
    #[must_use]
    pub const fn alternative_outcome_artifact(&self) -> ContentHash {
        self.alternative_outcome_artifact
    }
}

#[allow(clippy::too_many_arguments)]
fn validate_audit_receipt_inputs(
    ledger: &Ledger,
    snapshot: &fs_ledger::ViewSnapshot,
    record: &MatchedAuditRecord,
    context: &VoiLedgerContext,
    decision_artifact: ContentHash,
    model_artifact: ContentHash,
    recommended_outcome_artifact: ContentHash,
    alternative_outcome_artifact: ContentHash,
) -> Result<(Vec<u8>, Vec<u8>), VoiError> {
    if recommended_outcome_artifact == alternative_outcome_artifact {
        return Err(VoiError::ReceiptMismatch {
            field: "distinct matched outcome artifacts",
        });
    }
    let recommended = voi_audit_outcome_bytes(record, true)?;
    let alternative = voi_audit_outcome_bytes(record, false)?;
    for (role, artifact) in [
        ("decision", decision_artifact),
        ("model", model_artifact),
        ("recommended outcome", recommended_outcome_artifact),
        ("alternative outcome", alternative_outcome_artifact),
    ] {
        require_snapshot_artifact(snapshot, role, artifact)?;
    }
    verify_artifact_envelope(
        ledger,
        decision_artifact,
        VOI_DECISION_ARTIFACT_KIND,
        &voi_decision_artifact_metadata(&context.decision_version)?,
        None,
        u64::try_from(MAX_VOI_AUTHORITY_ARTIFACT_BYTES).unwrap_or(u64::MAX),
    )?;
    verify_artifact_envelope(
        ledger,
        model_artifact,
        VOI_MODEL_ARTIFACT_KIND,
        &voi_model_artifact_metadata(&context.model_version)?,
        None,
        u64::try_from(MAX_VOI_AUTHORITY_ARTIFACT_BYTES).unwrap_or(u64::MAX),
    )?;
    verify_artifact_envelope(
        ledger,
        recommended_outcome_artifact,
        VOI_AUDIT_OUTCOME_ARTIFACT_KIND,
        VOI_AUDIT_OUTCOME_ARTIFACT_METADATA,
        Some(&recommended),
        u64::try_from(MAX_VOI_LEDGER_RECEIPT_BYTES).unwrap_or(u64::MAX),
    )?;
    verify_artifact_envelope(
        ledger,
        alternative_outcome_artifact,
        VOI_AUDIT_OUTCOME_ARTIFACT_KIND,
        VOI_AUDIT_OUTCOME_ARTIFACT_METADATA,
        Some(&alternative),
        u64::try_from(MAX_VOI_LEDGER_RECEIPT_BYTES).unwrap_or(u64::MAX),
    )?;
    Ok((recommended, alternative))
}

/// Mint or exactly replay one prospective matched-outcome receipt through a
/// sealed ledger op.
///
/// # Errors
/// [`VoiError`] when the snapshot is stale, an outcome artifact does not
/// authenticate the exact record side, or ledger persistence fails.
#[allow(clippy::too_many_arguments)]
pub fn record_prospective_audit_receipt(
    ledger: &Ledger,
    vcs: &Vcs,
    record: &MatchedAuditRecord,
    sequence: usize,
    context: &VoiLedgerContext,
    decision_artifact: ContentHash,
    model_artifact: ContentHash,
    recommended_outcome_artifact: ContentHash,
    alternative_outcome_artifact: ContentHash,
) -> Result<VoiLedgerReceipt, VoiError> {
    if sequence >= MAX_VOI_AUDIT_RECORDS {
        return Err(VoiError::SizeLimit {
            collection: "VoI audit receipt sequence",
            count: sequence.saturating_add(1),
            min: 1,
            max: MAX_VOI_AUDIT_RECORDS,
        });
    }
    let (ledger_id, snapshot) = verify_snapshot(ledger, vcs, context, context.issued_day, true)?;
    let _outcomes = validate_audit_receipt_inputs(
        ledger,
        &snapshot,
        record,
        context,
        decision_artifact,
        model_artifact,
        recommended_outcome_artifact,
        alternative_outcome_artifact,
    )?;
    let canonical = audit_receipt_bytes(
        ledger_id,
        record,
        sequence,
        context,
        decision_artifact,
        model_artifact,
        recommended_outcome_artifact,
        alternative_outcome_artifact,
    )?;
    persist_receipt(
        ledger,
        context,
        AUDIT_RECEIPT_KIND,
        AUDIT_RECEIPT_META,
        &canonical,
        &[
            decision_artifact,
            model_artifact,
            recommended_outcome_artifact,
            alternative_outcome_artifact,
        ],
    )
}

#[allow(clippy::too_many_arguments)]
fn verify_audit_receipt_inner(
    ledger: &Ledger,
    vcs: &Vcs,
    record: &MatchedAuditRecord,
    sequence: usize,
    context: &VoiLedgerContext,
    decision_artifact: ContentHash,
    model_artifact: ContentHash,
    recommended_outcome_artifact: ContentHash,
    alternative_outcome_artifact: ContentHash,
    receipt: VoiLedgerReceipt,
    current_day: u32,
    require_current_head: bool,
) -> Result<(), VoiError> {
    if sequence >= MAX_VOI_AUDIT_RECORDS {
        return Err(VoiError::SizeLimit {
            collection: "VoI audit receipt sequence",
            count: sequence.saturating_add(1),
            min: 1,
            max: MAX_VOI_AUDIT_RECORDS,
        });
    }
    let (ledger_id, snapshot) =
        verify_snapshot(ledger, vcs, context, current_day, require_current_head)?;
    let _outcomes = validate_audit_receipt_inputs(
        ledger,
        &snapshot,
        record,
        context,
        decision_artifact,
        model_artifact,
        recommended_outcome_artifact,
        alternative_outcome_artifact,
    )?;
    let canonical = audit_receipt_bytes(
        ledger_id,
        record,
        sequence,
        context,
        decision_artifact,
        model_artifact,
        recommended_outcome_artifact,
        alternative_outcome_artifact,
    )?;
    verify_receipt_envelope(
        ledger,
        context,
        receipt,
        AUDIT_RECEIPT_KIND,
        AUDIT_RECEIPT_META,
        &canonical,
        &[
            decision_artifact,
            model_artifact,
            recommended_outcome_artifact,
            alternative_outcome_artifact,
        ],
    )
}

/// Independently verify one prospective outcome receipt before it can enter
/// the scheduler's authenticated e-process.
///
/// # Errors
/// [`VoiError`] for any snapshot, expiry, outcome byte, receipt, operation, or
/// lineage mismatch.
#[allow(clippy::too_many_arguments)]
pub fn verify_prospective_audit_receipt(
    ledger: &Ledger,
    vcs: &Vcs,
    record: MatchedAuditRecord,
    sequence: usize,
    context: VoiLedgerContext,
    decision_artifact: ContentHash,
    model_artifact: ContentHash,
    recommended_outcome_artifact: ContentHash,
    alternative_outcome_artifact: ContentHash,
    receipt: VoiLedgerReceipt,
    current_day: u32,
) -> Result<AuthenticatedAuditRecord, VoiError> {
    verify_audit_receipt_inner(
        ledger,
        vcs,
        &record,
        sequence,
        &context,
        decision_artifact,
        model_artifact,
        recommended_outcome_artifact,
        alternative_outcome_artifact,
        receipt,
        current_day,
        true,
    )?;
    Ok(AuthenticatedAuditRecord {
        record,
        sequence,
        context,
        receipt,
        decision_artifact,
        model_artifact,
        recommended_outcome_artifact,
        alternative_outcome_artifact,
    })
}

/// Single-owner live VoI audit and purchase scheduler.
///
/// Audit observations enter one append-only [`PairwiseRace`] in prospective
/// order. The scheduler owns the remaining budget and the bounded set of
/// decision snapshots it has already evaluated, so one process cannot reuse a
/// stale ranking or reset the budget through the safe API. The original
/// [`VoiScheduler::observe_audit`] / [`VoiScheduler::schedule`] path remains an
/// explicit caller-supplied reporting/advisory surface. Ledger authority is
/// minted only by [`VoiScheduler::observe_authenticated_audit`] followed by
/// [`VoiScheduler::schedule_authenticated`].
#[derive(Debug)]
pub struct VoiScheduler {
    policy_scope: String,
    remaining_budget_dollars: f64,
    audit_records: Vec<MatchedAuditRecord>,
    observation_ids: BTreeSet<String>,
    race: PairwiseRace,
    consumed_snapshots: BTreeSet<String>,
    authenticated_audits: Vec<AuthenticatedAuditRecord>,
}

impl VoiScheduler {
    /// Create one live scheduler for a fixed policy/version and total budget.
    ///
    /// # Errors
    /// [`VoiError`] when the policy identity or budget is malformed.
    pub fn new(policy_scope: impl Into<String>, budget: f64) -> Result<Self, VoiError> {
        let policy_scope = policy_scope.into();
        validate_name("VoI policy scope", 0, &policy_scope)?;
        if !budget.is_finite() || budget < 0.0 {
            return Err(VoiError::InvalidBudget { budget });
        }
        Ok(Self {
            policy_scope: policy_scope.as_str().to_owned(),
            remaining_budget_dollars: budget,
            audit_records: Vec::new(),
            observation_ids: BTreeSet::new(),
            race: PairwiseRace::new(LossSpan::ONE),
            consumed_snapshots: BTreeSet::new(),
            authenticated_audits: Vec::new(),
        })
    }

    /// Append one prospectively ordered matched-cost result to the one live
    /// e-process. Duplicate identities and limit+1 refuse before wealth or
    /// retained state changes.
    ///
    /// # Errors
    /// [`VoiError`] on a duplicate observation, the audit cap, or invalid
    /// e-process arithmetic.
    pub fn observe_audit(&mut self, record: MatchedAuditRecord) -> Result<(), VoiError> {
        if self.audit_records.len() >= MAX_VOI_AUDIT_RECORDS {
            return Err(VoiError::SizeLimit {
                collection: "VoI audit records",
                count: self.audit_records.len().saturating_add(1),
                min: 0,
                max: MAX_VOI_AUDIT_RECORDS,
            });
        }
        if self.observation_ids.contains(&record.observation_id) {
            return Err(VoiError::DuplicateAuditObservation {
                observation: record.observation_id.clone(),
            });
        }
        let recommended_loss = f64::from(u8::from(!record.recommended_changed_decision));
        let alternative_loss = f64::from(u8::from(!record.alternative_changed_decision));
        let mut next_race = self.race.clone();
        next_race
            .observe(recommended_loss, alternative_loss)
            .map_err(|_| VoiError::ArithmeticRefusal {
                operation: "VoI matched-cost e-process",
                subject: record.observation_id.clone(),
            })?;
        self.observation_ids.insert(record.observation_id.clone());
        self.audit_records.push(record);
        self.race = next_race;
        Ok(())
    }

    /// Append one independently verified prospective outcome to the only
    /// audit stream eligible for ledger-backed scheduling authority.
    ///
    /// The receipt sequence must equal the current chronological length.
    /// Mixing caller-supplied records into an authenticated scheduler fails
    /// closed, as do session/decision/model/policy, branch, or exact
    /// decision/model artifact changes within one stream.
    ///
    /// # Errors
    /// [`VoiError::AuditSequenceMismatch`] for reordered/skipped receipts;
    /// [`VoiError::ReceiptMismatch`] for mixed or cross-scope evidence; other
    /// audit validation errors are inherited from [`Self::observe_audit`].
    pub fn observe_authenticated_audit(
        &mut self,
        evidence: AuthenticatedAuditRecord,
    ) -> Result<(), VoiError> {
        let expected = self.audit_records.len();
        if evidence.sequence != expected {
            return Err(VoiError::AuditSequenceMismatch {
                expected,
                actual: evidence.sequence,
            });
        }
        if self.authenticated_audits.len() != self.audit_records.len() {
            return Err(VoiError::ReceiptMismatch {
                field: "mixed authenticated and caller-supplied audit stream",
            });
        }
        if evidence.context.policy_version != self.policy_scope {
            return Err(VoiError::PolicyScopeMismatch {
                expected: self.policy_scope.clone(),
                actual: evidence.context.policy_version.clone(),
            });
        }
        if let Some(first) = self.authenticated_audits.first()
            && (first.context.session_id != evidence.context.session_id
                || first.context.decision_version != evidence.context.decision_version
                || first.context.model_version != evidence.context.model_version
                || first.context.policy_version != evidence.context.policy_version
                || first.context.branch != evidence.context.branch
                || first.decision_artifact != evidence.decision_artifact
                || first.model_artifact != evidence.model_artifact)
        {
            return Err(VoiError::ReceiptMismatch {
                field: "authenticated audit stream scope",
            });
        }
        self.observe_audit(evidence.record.clone())?;
        self.authenticated_audits.push(evidence);
        Ok(())
    }

    /// Immutable reporting snapshot of the current chronological audit prefix.
    /// The report carries no scheduling capability.
    ///
    /// # Errors
    /// [`VoiError`] if the e-process or content identity leaves its finite
    /// bounded domain.
    pub fn audit_report(&self) -> Result<AuditReport, VoiError> {
        let audit_context_id = audit_context(&self.policy_scope, &self.audit_records)?;
        let log_e_value = self.race.log_e_value();
        if !log_e_value.is_finite() {
            return Err(VoiError::ArithmeticRefusal {
                operation: "VoI audit log e-value",
                subject: audit_context_id.to_hex(),
            });
        }
        Ok(AuditReport {
            policy_scope: self.policy_scope.clone(),
            audit_context_id,
            observations: self.audit_records.len(),
            log_e_value,
            verdict: if self.race.a_beats_b(VOI_AUDIT_ALPHA) {
                AuditVerdict::KeepScheduling
            } else {
                AuditVerdict::DemoteToReporting
            },
        })
    }

    /// Fixed policy/version scope owned by this scheduler.
    #[must_use]
    pub fn policy_scope(&self) -> &str {
        &self.policy_scope
    }

    /// Current unspent scheduler budget.
    #[must_use]
    pub fn remaining_budget_dollars(&self) -> f64 {
        self.remaining_budget_dollars
    }

    /// Number of decision snapshots already evaluated by this scheduler.
    #[must_use]
    pub fn consumed_snapshots(&self) -> usize {
        self.consumed_snapshots.len()
    }

    /// Evaluate at most one advisory purchase from one previously unseen
    /// caller-declared decision snapshot. The current live audit must still be
    /// above threshold. This raw path does not verify ledger provenance; use
    /// [`Self::schedule_authenticated`] before treating a result as authority.
    /// All validation and arithmetic precede mutation; success (including a
    /// no-affordable-purchase result) consumes the snapshot atomically.
    ///
    /// # Errors
    /// [`VoiError`] when the audit threshold is currently unmet, policy scopes
    /// differ, the snapshot was already consumed, retained snapshot capacity
    /// is exhausted, or budget arithmetic cannot decrease monotonically.
    pub fn schedule(&mut self, ranked: RankedMenu) -> Result<Option<ScheduledPurchase>, VoiError> {
        if ranked.policy_scope != self.policy_scope {
            return Err(VoiError::PolicyScopeMismatch {
                expected: self.policy_scope.clone(),
                actual: ranked.policy_scope,
            });
        }
        if self.consumed_snapshots.contains(&ranked.snapshot_id) {
            return Err(VoiError::RankingSnapshotAlreadyConsumed {
                snapshot_id: ranked.snapshot_id,
            });
        }
        if self.consumed_snapshots.len() >= MAX_VOI_SCHEDULED_CONTEXTS {
            return Err(VoiError::SizeLimit {
                collection: "VoI consumed ranking snapshots",
                count: self.consumed_snapshots.len().saturating_add(1),
                min: 0,
                max: MAX_VOI_SCHEDULED_CONTEXTS,
            });
        }
        if !self.race.a_beats_b(VOI_AUDIT_ALPHA) {
            return Err(VoiError::MissingSchedulingAuthority);
        }
        let audit = self.audit_report()?;
        let budget = self.remaining_budget_dollars;
        let purchase = ranked
            .rows
            .iter()
            .find(|row| row.score > 0.0 && row.probe.cost <= budget)
            .cloned();
        let remaining = if let Some(purchase) = &purchase {
            let remaining = budget - purchase.probe.cost;
            if !remaining.is_finite() || remaining < 0.0 || remaining >= budget {
                return Err(VoiError::ArithmeticRefusal {
                    operation: "remaining-budget subtraction",
                    subject: purchase.probe.name.clone(),
                });
            }
            remaining
        } else {
            budget
        };

        self.consumed_snapshots.insert(ranked.snapshot_id.clone());
        self.remaining_budget_dollars = remaining;
        Ok(purchase.map(|purchase| ScheduledPurchase {
            purchase,
            ranked_context_id: ranked.context_id,
            ranked_source_context_id: ranked.source_context_id,
            policy_scope: ranked.policy_scope,
            snapshot_id: ranked.snapshot_id,
            ranked_grid: ranked.grid,
            audit_context_id: audit.audit_context_id,
            audit_observations: audit.observations,
            audit_log_e_value: audit.log_e_value,
            budget_dollars: budget,
            remaining_budget_dollars: remaining,
        }))
    }

    /// Reverify the current ranked snapshot and every prospective audit
    /// receipt, then atomically consume at most one affordable purchase.
    ///
    /// The ranked receipt must still bind the current branch head. Historical
    /// audit roots may be older, but must remain known immutable snapshots on
    /// the same branch and all receipts must share one session, decision,
    /// model, and policy scope plus the exact decision/model artifact hashes.
    /// No scheduler state changes until every receipt and expiry check has
    /// passed.
    ///
    /// # Errors
    /// [`VoiError`] for incomplete/mixed audit authority, a cross-scope
    /// receipt, stale or expired evidence, any independently detected ledger
    /// mismatch, or the validation failures documented by [`Self::schedule`].
    pub fn schedule_authenticated(
        &mut self,
        ledger: &Ledger,
        vcs: &Vcs,
        evidence: AuthenticatedRankedMenu,
        current_day: u32,
    ) -> Result<Option<AuthenticatedScheduledPurchase>, VoiError> {
        if self.authenticated_audits.is_empty()
            || self.authenticated_audits.len() != self.audit_records.len()
        {
            return Err(VoiError::MissingSchedulingAuthority);
        }

        let first_audit = &self.authenticated_audits[0];
        if evidence.context.session_id != first_audit.context.session_id
            || evidence.context.decision_version != first_audit.context.decision_version
            || evidence.context.model_version != first_audit.context.model_version
            || evidence.context.policy_version != first_audit.context.policy_version
            || evidence.context.branch != first_audit.context.branch
            || evidence.decision_artifact != first_audit.decision_artifact
            || evidence.model_artifact != first_audit.model_artifact
        {
            return Err(VoiError::ReceiptMismatch {
                field: "ranked and audit authority scope",
            });
        }

        let (ledger_id, snapshot) =
            verify_snapshot(ledger, vcs, &evidence.context, current_day, true)?;
        let _catalog = validate_ranked_receipt_inputs(
            ledger,
            &snapshot,
            &evidence.ranked,
            &evidence.context,
            evidence.decision_artifact,
            evidence.model_artifact,
            evidence.probe_catalog_artifact,
        )?;
        let canonical = ranked_receipt_bytes(
            ledger_id,
            &evidence.ranked,
            &evidence.context,
            evidence.decision_artifact,
            evidence.model_artifact,
            evidence.probe_catalog_artifact,
        )?;
        verify_receipt_envelope(
            ledger,
            &evidence.context,
            evidence.receipt,
            RANKED_RECEIPT_KIND,
            RANKED_RECEIPT_META,
            &canonical,
            &[
                evidence.decision_artifact,
                evidence.model_artifact,
                evidence.probe_catalog_artifact,
            ],
        )?;

        let mut valid_through_day = evidence.context.expires_day;
        let mut audit_receipts = Vec::with_capacity(self.authenticated_audits.len());
        for audit in &self.authenticated_audits {
            if audit.context.session_id != evidence.context.session_id
                || audit.context.decision_version != evidence.context.decision_version
                || audit.context.model_version != evidence.context.model_version
                || audit.context.policy_version != evidence.context.policy_version
                || audit.context.branch != evidence.context.branch
                || audit.decision_artifact != evidence.decision_artifact
                || audit.model_artifact != evidence.model_artifact
            {
                return Err(VoiError::ReceiptMismatch {
                    field: "authenticated audit authority scope",
                });
            }
            verify_audit_receipt_inner(
                ledger,
                vcs,
                &audit.record,
                audit.sequence,
                &audit.context,
                audit.decision_artifact,
                audit.model_artifact,
                audit.recommended_outcome_artifact,
                audit.alternative_outcome_artifact,
                audit.receipt,
                current_day,
                false,
            )?;
            valid_through_day = valid_through_day.min(audit.context.expires_day);
            audit_receipts.push(audit.receipt);
        }

        let ranked_receipt = evidence.receipt;
        self.schedule(evidence.ranked).map(|scheduled| {
            scheduled.map(|scheduled| AuthenticatedScheduledPurchase {
                scheduled,
                ranked_receipt,
                audit_receipts,
                valid_through_day,
            })
        })
    }
}

/// Recompute a bounded reporting-only audit from a supplied chronological
/// prefix. This helper never returns scheduling authority; executable callers
/// must retain one live [`VoiScheduler`] and append observations to it.
///
/// # Errors
/// [`VoiError`] for malformed policy identity, duplicate/oversized records, or
/// invalid e-process arithmetic.
pub fn audit_scheduling(
    policy_scope: &str,
    records: &[MatchedAuditRecord],
) -> Result<AuditReport, VoiError> {
    let mut scheduler = VoiScheduler::new(policy_scope, 0.0)?;
    for record in records {
        scheduler.observe_audit(record.clone())?;
    }
    scheduler.audit_report()
}

#[cfg(test)]
mod retained_capacity_tests {
    use super::*;

    fn spare_capacity(value: &str) -> String {
        let mut out = String::with_capacity(4096);
        out.push_str(value);
        out
    }

    fn next_up(value: f64) -> f64 {
        assert!(value.is_finite() && value >= 0.0);
        let bits = value.to_bits();
        let rendered = format!("{value:.12}");
        let next = f64::from_bits(bits + 1);
        assert_ne!(
            next.to_bits(),
            bits,
            "one-ULP mutation must change exact bits"
        );
        assert_eq!(next.to_bits(), bits + 1, "mutation must be exactly one ULP");
        assert_eq!(
            format!("{next:.12}"),
            rendered,
            "chosen display rendering must hide the exact-bit mutation"
        );
        next
    }

    #[derive(Clone)]
    struct SourceIdentityFixture {
        domain: String,
        producer_version: u32,
        nodes: Vec<UncertaintyNode>,
        probes: Vec<Probe>,
        grid: usize,
        policy_scope: String,
        snapshot_id: String,
        metadata: DecisionOracleMetadata,
        computation: DecisionComputationReceipt,
    }

    impl SourceIdentityFixture {
        fn root(&self) -> ContentHash {
            self.root_with_declared_counts(self.nodes.len(), self.probes.len())
        }

        fn root_with_declared_counts(
            &self,
            declared_node_count: usize,
            declared_probe_count: usize,
        ) -> ContentHash {
            ranked_source_context_with_declared_counts(
                &self.domain,
                self.producer_version,
                &self.nodes,
                declared_node_count,
                &self.probes,
                declared_probe_count,
                self.grid,
                &self.policy_scope,
                &self.snapshot_id,
                self.metadata,
                self.computation,
            )
            .expect("bounded ranked-source identity fixture")
        }
    }

    fn source_identity_fixture() -> SourceIdentityFixture {
        SourceIdentityFixture {
            domain: RANKED_MENU_SOURCE_DOMAIN.to_string(),
            producer_version: VOI_RANKED_SOURCE_IDENTITY_VERSION,
            nodes: vec![
                UncertaintyNode {
                    name: "a".to_string(),
                    lo: 0.25,
                    hi: 2.0,
                    nominal: 1.0,
                },
                UncertaintyNode {
                    name: "b".to_string(),
                    lo: 2.0,
                    hi: 4.0,
                    nominal: 3.0,
                },
            ],
            probes: vec![
                Probe {
                    name: "alpha".to_string(),
                    target: "a".to_string(),
                    cost: 2.0,
                    shrink: 0.5,
                    kind: ProbeKind::Computational,
                },
                Probe {
                    name: "beta".to_string(),
                    target: "b".to_string(),
                    cost: 4.0,
                    shrink: 0.25,
                    kind: ProbeKind::Physical,
                },
            ],
            grid: 8,
            policy_scope: "policy-v1".to_string(),
            snapshot_id: "snapshot-v1".to_string(),
            metadata: DecisionOracleMetadata {
                arity: 2,
                work_units_per_evaluation: 3,
            },
            computation: DecisionComputationReceipt {
                evaluations: 9,
                work_units: 27,
                budget: DecisionBudget {
                    max_evaluations: 16,
                    max_work_units: 64,
                },
            },
        }
    }

    fn assert_source_identity_moves(mutate: impl FnOnce(&mut SourceIdentityFixture)) {
        let fixture = source_identity_fixture();
        let baseline = fixture.root();
        let mut changed = fixture;
        mutate(&mut changed);
        assert_ne!(baseline, changed.root());
    }

    fn ranked_row(name: &str, flip_before: f64, flip_after: f64, score: f64) -> RankedPurchase {
        RankedPurchase {
            probe: Probe {
                name: name.to_string(),
                target: "a".to_string(),
                cost: 2.0,
                shrink: 0.5,
                kind: ProbeKind::Computational,
            },
            flip_before,
            flip_after,
            score,
        }
    }

    #[derive(Clone)]
    struct RankedIdentityFixture {
        domain: String,
        producer_version: u32,
        source_context_id: ContentHash,
        rows: Vec<RankedPurchase>,
    }

    impl RankedIdentityFixture {
        fn root(&self) -> ContentHash {
            self.root_with_declared_count(self.rows.len())
        }

        fn root_with_declared_count(&self, declared_row_count: usize) -> ContentHash {
            ranked_output_context_with_declared_count(
                &self.domain,
                self.producer_version,
                self.source_context_id,
                &self.rows,
                declared_row_count,
            )
            .expect("bounded ranked-menu identity fixture")
        }
    }

    fn ranked_identity_fixture() -> RankedIdentityFixture {
        RankedIdentityFixture {
            domain: RANKED_MENU_CONTEXT_DOMAIN.to_string(),
            producer_version: VOI_RANKED_MENU_IDENTITY_VERSION,
            source_context_id: hash_domain("identity-test-source", b"source-a"),
            rows: vec![
                ranked_row("alpha", 0.75, 0.25, 0.125),
                ranked_row("beta", 0.5, 0.125, 0.0625),
            ],
        }
    }

    fn assert_ranked_identity_moves(mutate: impl FnOnce(&mut RankedIdentityFixture)) {
        let fixture = ranked_identity_fixture();
        let baseline = fixture.root();
        let mut changed = fixture;
        mutate(&mut changed);
        assert_ne!(baseline, changed.root());
    }

    fn audit_record_fixture(id: &str, cost: f64, recommended_wins: bool) -> MatchedAuditRecord {
        MatchedAuditRecord::new(
            id,
            format!("recommended-{id}"),
            format!("alternative-{id}"),
            format!("provenance-{id}"),
            cost,
            cost,
            recommended_wins,
            !recommended_wins,
        )
        .expect("valid audit identity fixture")
    }

    #[derive(Clone)]
    struct AuditIdentityFixture {
        domain: String,
        producer_version: u32,
        policy_scope: String,
        records: Vec<MatchedAuditRecord>,
        audit_alpha: f64,
        max_audit_records: usize,
    }

    impl AuditIdentityFixture {
        fn root(&self) -> ContentHash {
            self.root_with_declared_count(self.records.len())
        }

        fn root_with_declared_count(&self, declared_record_count: usize) -> ContentHash {
            audit_context_with_declared_count(
                &self.domain,
                self.producer_version,
                &self.policy_scope,
                &self.records,
                declared_record_count,
                self.audit_alpha,
                self.max_audit_records,
            )
            .expect("bounded audit identity fixture")
        }
    }

    fn audit_identity_fixture() -> AuditIdentityFixture {
        AuditIdentityFixture {
            domain: AUDIT_CONTEXT_DOMAIN.to_string(),
            producer_version: VOI_AUDIT_CONTEXT_IDENTITY_VERSION,
            policy_scope: "policy-v1".to_string(),
            records: vec![
                audit_record_fixture("a", 2.0, true),
                audit_record_fixture("b", 4.0, false),
            ],
            audit_alpha: VOI_AUDIT_ALPHA,
            max_audit_records: MAX_VOI_AUDIT_RECORDS,
        }
    }

    fn assert_audit_identity_moves(mutate: impl FnOnce(&mut AuditIdentityFixture)) {
        let fixture = audit_identity_fixture();
        let baseline = fixture.root();
        let mut changed = fixture;
        mutate(&mut changed);
        assert_ne!(baseline, changed.root());
    }

    #[test]
    fn voi_ranked_source_identity_mutation_battery() {
        assert_source_identity_moves(|fixture| fixture.domain.push_str(".alternate"));
        assert_source_identity_moves(|fixture| fixture.producer_version += 1);
        assert_source_identity_moves(|fixture| fixture.policy_scope.push('x'));
        assert_source_identity_moves(|fixture| fixture.snapshot_id.push('x'));
        assert_source_identity_moves(|fixture| fixture.grid += 1);
        assert_source_identity_moves(|fixture| fixture.metadata.arity += 1);
        assert_source_identity_moves(|fixture| fixture.metadata.work_units_per_evaluation += 1);
        assert_source_identity_moves(|fixture| fixture.computation.evaluations += 1);
        assert_source_identity_moves(|fixture| fixture.computation.work_units += 1);
        assert_source_identity_moves(|fixture| fixture.computation.budget.max_evaluations += 1);
        assert_source_identity_moves(|fixture| fixture.computation.budget.max_work_units += 1);

        let fixture = source_identity_fixture();
        assert_ne!(
            fixture.root(),
            fixture.root_with_declared_counts(fixture.nodes.len() + 1, fixture.probes.len()),
            "the node-count frame must move independently of node bytes",
        );
        assert_source_identity_moves(|fixture| fixture.nodes.reverse());
        assert_source_identity_moves(|fixture| fixture.nodes[0].name.push('x'));
        assert_source_identity_moves(|fixture| fixture.nodes[0].lo = next_up(fixture.nodes[0].lo));
        assert_source_identity_moves(|fixture| {
            fixture.nodes[0].nominal = next_up(fixture.nodes[0].nominal);
        });
        assert_source_identity_moves(|fixture| fixture.nodes[0].hi = next_up(fixture.nodes[0].hi));

        let fixture = source_identity_fixture();
        assert_ne!(
            fixture.root(),
            fixture.root_with_declared_counts(fixture.nodes.len(), fixture.probes.len() + 1),
            "the probe-count frame must move independently of probe bytes",
        );
        assert_source_identity_moves(|fixture| fixture.probes[0].name.push('x'));
        assert_source_identity_moves(|fixture| fixture.probes[0].target = "b".to_string());
        assert_source_identity_moves(|fixture| {
            fixture.probes[0].cost = next_up(fixture.probes[0].cost);
        });
        assert_source_identity_moves(|fixture| {
            fixture.probes[0].shrink = next_up(fixture.probes[0].shrink);
        });
        assert_source_identity_moves(|fixture| fixture.probes[0].kind = ProbeKind::Physical);
    }

    #[test]
    fn voi_ranked_source_menu_input_order_is_nonsemantic() {
        let fixture = source_identity_fixture();
        let baseline = fixture.root();
        let mut reordered = fixture;
        reordered.probes.reverse();
        assert_eq!(baseline, reordered.root());
    }

    #[test]
    fn voi_identity_allocation_capacity_is_nonsemantic() {
        let source = source_identity_fixture();
        let mut source_with_spare_capacity = source.clone();
        source_with_spare_capacity.nodes = Vec::with_capacity(64);
        source_with_spare_capacity
            .nodes
            .extend(source.nodes.iter().cloned());
        source_with_spare_capacity.probes = Vec::with_capacity(64);
        source_with_spare_capacity
            .probes
            .extend(source.probes.iter().cloned());
        assert_eq!(source.root(), source_with_spare_capacity.root());

        let ranked = ranked_identity_fixture();
        let mut ranked_with_spare_capacity = ranked.clone();
        ranked_with_spare_capacity.rows = Vec::with_capacity(64);
        ranked_with_spare_capacity
            .rows
            .extend(ranked.rows.iter().cloned());
        assert_eq!(ranked.root(), ranked_with_spare_capacity.root());

        let audit = audit_identity_fixture();
        let mut audit_with_spare_capacity = audit.clone();
        audit_with_spare_capacity.records = Vec::with_capacity(64);
        audit_with_spare_capacity
            .records
            .extend(audit.records.iter().cloned());
        assert_eq!(audit.root(), audit_with_spare_capacity.root());
    }

    #[test]
    fn voi_ranked_menu_probe_payload_is_bound_by_source_context() {
        let fixture = ranked_identity_fixture();
        let baseline = fixture.root();
        let mut changed = fixture;
        changed.rows[0].probe.target = "different-source-target".to_string();
        changed.rows[0].probe.cost = changed.rows[0].probe.cost.next_up();
        changed.rows[0].probe.shrink = changed.rows[0].probe.shrink.next_up();
        changed.rows[0].probe.kind = ProbeKind::Physical;
        assert_eq!(baseline, changed.root());
    }

    #[test]
    fn voi_ranked_menu_identity_mutation_battery() {
        assert_ranked_identity_moves(|fixture| fixture.domain.push_str(".alternate"));
        assert_ranked_identity_moves(|fixture| fixture.producer_version += 1);
        assert_ranked_identity_moves(|fixture| {
            fixture.source_context_id = hash_domain("identity-test-source", b"source-b");
        });
        let fixture = ranked_identity_fixture();
        assert_ne!(
            fixture.root(),
            fixture.root_with_declared_count(fixture.rows.len() + 1),
            "the ranked-row count frame must move independently of row bytes",
        );
        assert_ranked_identity_moves(|fixture| fixture.rows.reverse());
        assert_ranked_identity_moves(|fixture| fixture.rows[0].probe.name.push('x'));
        assert_ranked_identity_moves(|fixture| {
            fixture.rows[0].flip_before = next_up(fixture.rows[0].flip_before);
        });
        assert_ranked_identity_moves(|fixture| {
            fixture.rows[0].flip_after = next_up(fixture.rows[0].flip_after);
        });
        assert_ranked_identity_moves(|fixture| {
            fixture.rows[0].score = next_up(fixture.rows[0].score);
        });
    }

    #[test]
    fn voi_audit_context_identity_mutation_battery() {
        assert_audit_identity_moves(|fixture| fixture.domain.push_str(".alternate"));
        assert_audit_identity_moves(|fixture| fixture.producer_version += 1);
        assert_audit_identity_moves(|fixture| fixture.policy_scope.push('x'));
        assert_audit_identity_moves(|fixture| {
            fixture.audit_alpha = next_up(fixture.audit_alpha);
        });
        assert_audit_identity_moves(|fixture| fixture.max_audit_records += 1);
        let fixture = audit_identity_fixture();
        assert_ne!(
            fixture.root(),
            fixture.root_with_declared_count(fixture.records.len() + 1),
            "the audit-record count frame must move independently of record bytes",
        );
        assert_audit_identity_moves(|fixture| fixture.records.reverse());
        assert_audit_identity_moves(|fixture| fixture.records[0].observation_id.push('x'));
        assert_audit_identity_moves(|fixture| fixture.records[0].recommended_id.push('x'));
        assert_audit_identity_moves(|fixture| fixture.records[0].alternative_id.push('x'));
        assert_audit_identity_moves(|fixture| fixture.records[0].provenance.push('x'));
        assert_audit_identity_moves(|fixture| {
            fixture.records[0].matched_cost = next_up(fixture.records[0].matched_cost);
        });
        assert_audit_identity_moves(|fixture| {
            fixture.records[0].recommended_changed_decision =
                !fixture.records[0].recommended_changed_decision;
        });
        assert_audit_identity_moves(|fixture| {
            fixture.records[0].alternative_changed_decision =
                !fixture.records[0].alternative_changed_decision;
        });
    }

    #[test]
    fn voi_audit_report_verdict_is_derived_and_nonsemantic() {
        let fixture = audit_identity_fixture();
        let root = fixture.root();
        let report = AuditReport {
            policy_scope: fixture.policy_scope,
            audit_context_id: root,
            observations: fixture.records.len(),
            log_e_value: 0.0,
            verdict: AuditVerdict::DemoteToReporting,
        };
        let mut changed = report.clone();
        changed.verdict = AuditVerdict::KeepScheduling;
        assert_eq!(report.audit_context_id(), changed.audit_context_id());
        assert_ne!(report.verdict(), changed.verdict());
    }

    #[test]
    fn audit_authority_rebuilds_caller_owned_string_capacity() {
        let record = MatchedAuditRecord::new(
            spare_capacity("obs-1"),
            spare_capacity("recommended"),
            spare_capacity("alternative"),
            spare_capacity("ledger-row-1"),
            1.0,
            1.0,
            true,
            false,
        )
        .expect("bounded matched audit record");
        for value in [
            &record.observation_id,
            &record.recommended_id,
            &record.alternative_id,
            &record.provenance,
        ] {
            assert!(value.capacity() <= MAX_VOI_NAME_BYTES);
        }

        let scheduler =
            VoiScheduler::new(spare_capacity("policy-v1"), 1.0).expect("bounded scheduler policy");
        assert!(scheduler.policy_scope.capacity() <= MAX_VOI_NAME_BYTES);
    }
}
