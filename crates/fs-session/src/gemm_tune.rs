//! The production GEMM autotune loop (bead yqug): measure → cache →
//! model → dispatch, closed end-to-end.
//!
//! [`gemm_f64_session`] is the production consumer the tuner was built
//! for: it resolves an MC/NC [`GemmBlockPlan`] for the caller's shape
//! class (pins beat cached rows beat the documented cold-start default),
//! runs a BOUNDED candidate sweep when the machine is cold, records the
//! ranked wall-time evidence as a tune row, applies the caller's explicit
//! cache policy, and dispatches
//! `fs_la::gemm_f64_parallel_with_pool` with the selected plan.
//!
//! Honesty boundaries, in the fs-exec tuner's division:
//! - The KERNEL KEY embeds fs-la's `GEMM_BIT_SEMANTICS_VERSION`, so rows
//!   measured under a different accumulation contract can never match a
//!   lookup (semantic filtering by construction). Rows are additionally
//!   bound to the exact probe dims, requested/normalized thread budget,
//!   resolved ISA tier, placement policy, and implementation version, then
//!   the exact compiler/profile/codegen build fingerprint, then machine-
//!   fingerprint-keyed. The ledger read path refuses stale,
//!   differently scoped, non-canonical, and params/body-disagreeing rows.
//! - MC/NC are BIT-NEUTRAL by fs-la's determinism contract, and the
//!   sweep ENFORCES that: every repeat of every effective candidate is
//!   compared word-for-word with the first output, else the loop fails
//!   closed with [`GemmTuneError::BitDrift`] and records nothing. KC is part
//!   of the bit contract and is NOT in this loop. The resolved SIMD tier is
//!   bit-neutral but remains performance identity.
//! - The "cost model" is declared and minimal: argmin of the per-
//!   candidate MINIMUM wall time, ties to the earlier candidate in
//!   lattice order — a recorded selection rule, never a statistical
//!   confidence claim.
//!
//! Determinism class: dispatch results are bit-identical to serial
//! `gemm_f64` for every plan the loop can select (enforced by the sweep
//! and gated in tests); WHICH plan wins is wall-clock-dependent by
//! nature and travels as evidence + a pinnable decision, never inside
//! numeric results.

use fs_exec::{
    CancelGate, GEMM_KERNEL_PREFIX, GemmBlockPlan, GemmExecutionIdentity, GemmTuneKey,
    PreparedGemmDecision, PreparedGemmRow, TilePool, TuneError, TuneEvidence, TuneObservation,
    TuneRow, TuneSource, Tuner,
};
use fs_ledger::Ledger;

pub use fs_la::{GEMM_DEPGRAPH_RECEIPT_DOMAIN, GemmGraphEvidenceClass};

/// The bounded sweep lattice: up to 4 × 2 candidates, lattice order
/// (mc-major ascending). Candidates that clamp to an identical effective
/// `(mc, nc)` pair are deduplicated before measurement. Chosen around the
/// measured xlvx s5 landscape: thin bands won both reference machines; the
/// extremes document the neighborhood.
const SWEEP_MC: [usize; 4] = [16, 32, 64, 128];
const SWEEP_NC_CAP: [usize; 2] = [512, 2048];

/// Probe M/K dims are capped so a cold-start sweep stays bounded (seconds,
/// not minutes) even when the caller's problem is huge. N has a separate
/// cap: it must extend beyond the smaller NC candidate or that axis is never
/// measured at all.
const PROBE_MK_DIM_CAP: usize = 512;
const PROBE_N_DIM_CAP: usize = 2048;

/// Wall-time samples per candidate (min-of ranking, all survive in the
/// evidence row).
const SWEEP_SAMPLES: usize = 3;

/// Session tune-metadata plan (bead wf9.15.1): the checked logical byte
/// bound for every session-owned autotune metadata structure of one sweep —
/// candidate/ranking/observation collections, the reused sample buffer,
/// canonical plan labels, and the sealed-row strings. All components derive
/// from the const sweep lattice and documented schema caps, so the plan is a
/// pure constant: a freshly measured row and the same row adopted later
/// derive the identical plan, keeping `receipt_identity` stable across both
/// paths. Numeric probe buffers are wf9.15's domain and generic TilePool
/// internals are wf9.16's; neither is claimed here.
mod tune_metadata_plan {
    use super::{SWEEP_MC, SWEEP_NC_CAP, SWEEP_SAMPLES, SweepCandidate};
    use fs_exec::{GemmBlockPlan, TuneObservation};

    /// Schema tag bound into the sealed-row receipt.
    pub const SCHEMA: &str = "fs-session-tune-metadata-plan-v1";

    /// Maximum distinct sweep candidates (the const lattice).
    pub const CANDIDATE_CAP: usize = SWEEP_MC.len() * SWEEP_NC_CAP.len();

    /// Documented cap on one canonical plan label (`mc=NNN,nc-cap=NNNN`);
    /// enforced at observation time, not assumed.
    pub const LABEL_BYTES_CAP: usize = 64;

    /// Documented cap on the sealed row's `params` JSON; enforced at seal
    /// time. The params string is the canonical plan spelling.
    pub const PARAMS_BYTES_CAP: usize = 256;

    /// Documented cap on the sealed row's `measured` evidence JSON
    /// (observations for the full lattice at `SWEEP_SAMPLES` samples plus
    /// fixed fields); enforced at seal time.
    pub const MEASURED_BYTES_CAP: usize = 16 * 1024;

    /// Total logical bytes the plan charges against the caller envelope
    /// before any probe allocation.
    #[must_use]
    pub fn requested_bytes() -> u128 {
        let candidates = CANDIDATE_CAP * size_of::<SweepCandidate>();
        let ranked = CANDIDATE_CAP * size_of::<(u64, usize, GemmBlockPlan)>();
        let observations = CANDIDATE_CAP * size_of::<TuneObservation>();
        let samples_reused = SWEEP_SAMPLES * size_of::<u64>();
        let labels = CANDIDATE_CAP * LABEL_BYTES_CAP;
        let sample_storage = CANDIDATE_CAP * SWEEP_SAMPLES * size_of::<u64>();
        let sealed_row = PARAMS_BYTES_CAP + MEASURED_BYTES_CAP;
        candidates as u128
            + ranked as u128
            + observations as u128
            + samples_reused as u128
            + labels as u128
            + sample_storage as u128
            + sealed_row as u128
    }

    /// Canonical receipt fragment for the sealed-row JSON.
    #[must_use]
    pub fn receipt_fragment() -> String {
        format!(
            "{{\"schema\":\"{SCHEMA}\",\"requested_bytes\":{}}}",
            requested_bytes()
        )
    }
}

/// Public schema tag of the session tune-metadata plan bound into every
/// sealed-row receipt (bead wf9.15.1).
pub const GEMM_TUNE_METADATA_PLAN_SCHEMA: &str = tune_metadata_plan::SCHEMA;

/// Logical bytes the session tune-metadata plan charges against the caller
/// envelope before any probe allocation. A pure constant of the sweep
/// lattice and documented schema caps — the same value fresh measurement
/// and later adoption bind into `receipt_identity`.
#[must_use]
pub fn gemm_tune_metadata_plan_bytes() -> u128 {
    tune_metadata_plan::requested_bytes()
}

/// Durable identity of the autotune producer algorithm: candidate lattice,
/// probe dimensions/sample policy, ranking, and plan-to-dispatch mapping. Any
/// semantic change to those choices must bump this value before old rows may be
/// considered compatible.
pub const GEMM_TUNER_SCHEMA_VERSION: u32 = 1;

/// Logical stream seed for the compatibility session pool. GEMM itself is
/// non-stochastic; caller-owned pools retain their study seed for Cx identity.
const SESSION_GEMM_POOL_SEED: u64 = 0x4653_2D53_4553_534E;

const GEMM_SWEEP_RUN_DOMAIN: &str = "org.frankensim.fs-session.gemm-sweep-run.v1";

/// Globally unique BLAKE3 derive-key context for canonical GEMM tune-row
/// receipts.
pub const GEMM_TUNE_ROW_RECEIPT_DOMAIN: &str = "org.frankensim.fs-session.gemm-tune-row-receipt.v2";

/// Current schema carried by retained GEMM tune-row receipt metadata.
pub const GEMM_TUNE_ROW_RECEIPT_IDENTITY_VERSION: u32 = 2;

/// Maximum exact retained tune-row receipt admitted by the replay parser.
const MAX_GEMM_TUNE_ROW_RECEIPT_BYTES: usize = 128 * 1024;
/// Maximum one outer receipt string admitted by the replay parser.
const MAX_GEMM_TUNE_ROW_RECEIPT_STRING_BYTES: usize = 64 * 1024;
/// Maximum nesting admitted while validating the embedded canonical row JSON.
const MAX_GEMM_TUNE_ROW_RECEIPT_JSON_DEPTH: usize = 64;

/// Schema of the tagged binary transport for deterministic GEMM execution facts.
pub const GEMM_EXECUTION_RECEIPT_IDENTITY_VERSION: u32 = 1;
/// BLAKE3 derive-key context for the deterministic GEMM execution receipt.
pub const GEMM_EXECUTION_RECEIPT_DOMAIN: &str =
    "org.frankensim.fs-session.gemm-execution-receipt.v1";
const GEMM_EXECUTION_RECEIPT_MAGIC: &[u8] = b"frankensim-gemm-execution-receipt";
const MAX_GEMM_EXECUTION_RECEIPT_BYTES: usize = 64 * 1024 * 1024;
const MAX_GEMM_EXECUTION_RECEIPT_PANELS: usize = 1 << 20;
const MAX_GEMM_EXECUTION_RECEIPT_STRING_BYTES: usize = 64 * 1024;

const EXEC_TAG_DOMAIN: u8 = 0x01;
const EXEC_TAG_VERSION: u8 = 0x02;
const EXEC_TAG_DECLARED_RUN: u8 = 0x10;
const EXEC_TAG_COMPLETED_TILES: u8 = 0x11;
const EXEC_TAG_TOTAL_TILES: u8 = 0x12;
const EXEC_TAG_MEMORY: u8 = 0x20;
const EXEC_TAG_MEMORY_LIMIT: u8 = 0x21;
const EXEC_TAG_MEMORY_STAGING: u8 = 0x22;
const EXEC_TAG_MEMORY_B_PACK: u8 = 0x23;
const EXEC_TAG_MEMORY_BAND_METADATA: u8 = 0x24;
const EXEC_TAG_MEMORY_POOL_RUN: u8 = 0x25;
const EXEC_TAG_MEMORY_ARENA_PER_WORKER: u8 = 0x26;
const EXEC_TAG_MEMORY_ACTIVE_WORKERS: u8 = 0x27;
const EXEC_TAG_MEMORY_ARENA: u8 = 0x28;
const EXEC_TAG_MEMORY_REQUESTED: u8 = 0x29;
const EXEC_TAG_PANELS: u8 = 0x30;
const EXEC_TAG_PANEL: u8 = 0x31;
const EXEC_TAG_PANEL_KERNEL: u8 = 0x32;
const EXEC_TAG_PANEL_MODE: u8 = 0x33;
const EXEC_TAG_PANEL_DECLARED_RUN: u8 = 0x34;
const EXEC_TAG_PANEL_COMPLETED: u8 = 0x35;
const EXEC_TAG_PANEL_TOTAL: u8 = 0x36;
const EXEC_TAG_END: u8 = 0xff;

/// Owner-local tune-row receipt declaration consumed by `xtask check-identities`.
#[allow(dead_code)]
pub const GEMM_TUNE_ROW_RECEIPT_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-session:gemm-tune-row-receipt",
    "version_const=GEMM_TUNE_ROW_RECEIPT_IDENTITY_VERSION",
    "version=2",
    "domain=org.frankensim.fs-session.gemm-tune-row-receipt.v2",
    "domain_const=GEMM_TUNE_ROW_RECEIPT_DOMAIN",
    "encoder=ValidatedGemmTuneRow::receipt_identity",
    "encoder_helpers=ValidatedGemmTuneRow::receipt_json,push_json_string,tune_metadata_plan::receipt_fragment,tune_metadata_plan::requested_bytes,parse_validated_gemm_tune_row_receipt,ExactJsonCursor::take,ExactJsonCursor::is_finished,ExactJsonCursor::canonical_string,ExactJsonCursor::take_hex_quad,ExactJsonCursor::canonical_u64,ExactJsonCursor::canonical_u128,ExactJsonCursor::canonical_value,ExactJsonCursor::parse_value,ExactJsonCursor::parse_number",
    "schema_constants=GEMM_TUNE_ROW_RECEIPT_IDENTITY_VERSION,GEMM_TUNE_ROW_RECEIPT_DOMAIN,GEMM_TUNE_METADATA_PLAN_SCHEMA,GEMM_TUNER_SCHEMA_VERSION,SWEEP_MC,SWEEP_NC_CAP,SWEEP_SAMPLES,PROBE_MK_DIM_CAP,PROBE_N_DIM_CAP,MAX_GEMM_TUNE_ROW_RECEIPT_BYTES,MAX_GEMM_TUNE_ROW_RECEIPT_STRING_BYTES,MAX_GEMM_TUNE_ROW_RECEIPT_JSON_DEPTH,crates/fs-exec/src/tune.rs#GEMM_TUNE_KEY_IDENTITY_VERSION,crates/fs-exec/src/tune.rs#TUNE_ROW_IDENTITY_VERSION,crates/fs-exec/src/tune.rs#TUNING_DECISION_IDENTITY_VERSION,crates/fs-exec/src/tune.rs#GEMM_TUNE_KEY_IDENTITY_DOMAIN,crates/fs-exec/src/tune.rs#TUNE_ROW_IDENTITY_DOMAIN,crates/fs-exec/src/tune.rs#TUNING_DECISION_IDENTITY_DOMAIN,crates/fs-blake3/src/lib.rs#IV,crates/fs-blake3/src/lib.rs#MSG_PERMUTATION,crates/fs-blake3/src/lib.rs#BLOCK_LEN,crates/fs-blake3/src/lib.rs#CHUNK_LEN,crates/fs-blake3/src/lib.rs#CHUNK_START,crates/fs-blake3/src/lib.rs#CHUNK_END,crates/fs-blake3/src/lib.rs#PARENT,crates/fs-blake3/src/lib.rs#ROOT,crates/fs-blake3/src/lib.rs#DERIVE_KEY_CONTEXT,crates/fs-blake3/src/lib.rs#DERIVE_KEY_MATERIAL,crates/fs-blake3/src/lib.rs#MAX_DEPTH",
    "schema_functions=ValidatedGemmTuneRow::from_prepared,ValidatedGemmTuneRow::matches_decision,ValidatedGemmTuneRow::matches_ledger_row,ValidatedGemmTuneRow::admit_receipt_json,parse_validated_gemm_tune_row_receipt,ExactJsonCursor::take,ExactJsonCursor::canonical_string,ExactJsonCursor::canonical_u64,ExactJsonCursor::canonical_u128,ExactJsonCursor::canonical_value,ExactJsonCursor::parse_value,ExactJsonCursor::parse_number,probe_buffer_bytes_for_dims,crates/fs-exec/src/tune.rs#PreparedGemmRow::key,crates/fs-exec/src/tune.rs#PreparedGemmRow::params_json,crates/fs-exec/src/tune.rs#PreparedGemmRow::row_json,crates/fs-exec/src/tune.rs#PreparedGemmDecision::key,crates/fs-exec/src/tune.rs#PreparedGemmDecision::plan,crates/fs-exec/src/tune.rs#PreparedGemmDecision::source,crates/fs-exec/src/tune.rs#GemmTuneKey::kernel,crates/fs-exec/src/tune.rs#GemmTuneKey::shape_class,crates/fs-exec/src/tune.rs#GemmTuneKey::execution,crates/fs-exec/src/tune.rs#TuneRow::from_canonical_json,crates/fs-exec/src/tune.rs#TuneRow::to_canonical_json,crates/fs-exec/src/tune.rs#TuneRow::kernel,crates/fs-exec/src/tune.rs#TuneRow::shape_class,crates/fs-exec/src/tune.rs#TuneRow::machine,crates/fs-exec/src/tune.rs#TuneRow::params,crates/fs-blake3/src/lib.rs#hash_domain,crates/fs-blake3/src/lib.rs#g,crates/fs-blake3/src/lib.rs#round,crates/fs-blake3/src/lib.rs#permute,crates/fs-blake3/src/lib.rs#compress,crates/fs-blake3/src/lib.rs#words_from_block,crates/fs-blake3/src/lib.rs#first_8_words,crates/fs-blake3/src/lib.rs#Output::chaining_value,crates/fs-blake3/src/lib.rs#Output::root_hash,crates/fs-blake3/src/lib.rs#parent_output,crates/fs-blake3/src/lib.rs#ChunkState::new,crates/fs-blake3/src/lib.rs#ChunkState::len,crates/fs-blake3/src/lib.rs#ChunkState::start_flag,crates/fs-blake3/src/lib.rs#ChunkState::update,crates/fs-blake3/src/lib.rs#ChunkState::output,crates/fs-blake3/src/lib.rs#Blake3::new_internal,crates/fs-blake3/src/lib.rs#Blake3::push_stack,crates/fs-blake3/src/lib.rs#Blake3::pop_stack,crates/fs-blake3/src/lib.rs#Blake3::add_chunk_chaining_value,crates/fs-blake3/src/lib.rs#Blake3::update,crates/fs-blake3/src/lib.rs#Blake3::finalize",
    "schema_dependencies=fs-exec:gemm-tune-key,fs-exec:tune-row,fs-exec:tuning-decision",
    "digest=fs-blake3",
    "encoding=canonical-transport-exact-bits",
    "sources=ValidatedGemmTuneRow",
    "source_fields=ValidatedGemmTuneRow.kernel:semantic,ValidatedGemmTuneRow.shape_class:semantic,ValidatedGemmTuneRow.machine:semantic,ValidatedGemmTuneRow.params:semantic,ValidatedGemmTuneRow.measured:semantic,ValidatedGemmTuneRow.memory_limit_bytes:semantic,ValidatedGemmTuneRow.probe_buffer_bytes:semantic",
    "source_bindings=ValidatedGemmTuneRow.kernel>kernel,ValidatedGemmTuneRow.shape_class>shape-class,ValidatedGemmTuneRow.machine>machine-fingerprint,ValidatedGemmTuneRow.params>selected-params,ValidatedGemmTuneRow.measured>measured-row,ValidatedGemmTuneRow.memory_limit_bytes>memory-limit-bytes,ValidatedGemmTuneRow.probe_buffer_bytes>probe-buffer-bytes",
    "external_semantic_fields=artifact-domain,identity-version,canonical-field-order,machine-hex-width,metadata-plan-schema,metadata-plan-requested-bytes",
    "semantic_fields=artifact-domain,identity-version,canonical-field-order,kernel,shape-class,machine-fingerprint,machine-hex-width,selected-params,measured-row,memory-limit-bytes,probe-buffer-bytes,metadata-plan-schema,metadata-plan-requested-bytes",
    "excluded_fields=none",
    "consumers=ValidatedGemmTuneRow::matches_decision,ValidatedGemmTuneRow::publish_to_ledger,ValidatedGemmTuneRow::publish_if_absent_or_identical,ValidatedGemmTuneRow::persist,ValidatedGemmTuneRow::replace_cache_row,install_sweep_row,adopt_cached_row,fs-roofline::KernelExecutionBinding",
    "mutations=artifact-domain:crates/fs-session/src/gemm_tune.rs#gemm_tune_row_receipt_identity_fields_move_independently,identity-version:crates/fs-session/src/gemm_tune.rs#gemm_tune_row_receipt_identity_fields_move_independently,canonical-field-order:crates/fs-session/src/gemm_tune.rs#gemm_tune_row_receipt_identity_fields_move_independently,kernel:crates/fs-session/src/gemm_tune.rs#gemm_tune_row_receipt_identity_fields_move_independently,shape-class:crates/fs-session/src/gemm_tune.rs#gemm_tune_row_receipt_identity_fields_move_independently,machine-fingerprint:crates/fs-session/src/gemm_tune.rs#gemm_tune_row_receipt_identity_fields_move_independently,machine-hex-width:crates/fs-session/src/gemm_tune.rs#gemm_tune_row_receipt_identity_fields_move_independently,selected-params:crates/fs-session/src/gemm_tune.rs#gemm_tune_row_receipt_identity_fields_move_independently,measured-row:crates/fs-session/src/gemm_tune.rs#gemm_tune_row_receipt_identity_fields_move_independently,memory-limit-bytes:crates/fs-session/src/gemm_tune.rs#gemm_tune_row_receipt_identity_fields_move_independently,probe-buffer-bytes:crates/fs-session/src/gemm_tune.rs#gemm_tune_row_receipt_identity_fields_move_independently,metadata-plan-schema:crates/fs-session/src/gemm_tune.rs#gemm_tune_row_receipt_identity_fields_move_independently,metadata-plan-requested-bytes:crates/fs-session/src/gemm_tune.rs#gemm_tune_row_receipt_identity_fields_move_independently",
    "nonsemantic_mutations=none",
    "field_guard=classify_validated_gemm_tune_row_identity_fields",
    "transport_guard=ValidatedGemmTuneRow::admit_receipt_json",
    "version_guard=crates/fs-session/src/gemm_tune.rs#gemm_tune_row_receipt_versions_fail_closed",
    "coupling_surface=fs-session:gemm-tune-row-receipt",
];

/// Owner-local execution-receipt declaration consumed by `xtask check-identities`.
#[allow(dead_code)]
pub const GEMM_EXECUTION_RECEIPT_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-session:gemm-execution-receipt",
    "version_const=GEMM_EXECUTION_RECEIPT_IDENTITY_VERSION",
    "version=1",
    "domain=org.frankensim.fs-session.gemm-execution-receipt.v1",
    "domain_const=GEMM_EXECUTION_RECEIPT_DOMAIN",
    "encoder=GemmExecutionReceipt::receipt_identity",
    "encoder_helpers=GemmDispatch::execution_receipt,GemmMemoryReceipt::from,GemmExecutionReceipt::from_report,GemmExecutionReceipt::canonical_bytes,GemmExecutionReceipt::canonical_bytes_with_schema,GemmExecutionReceipt::is_complete,GemmExecutionReceiptCodecError::fmt,execution_receipt_codec_error,checked_receipt_len_add,execution_receipt_encoded_len,push_execution_text,push_execution_u32,push_execution_u64,push_execution_u128,ExecutionReceiptCursor::take_exact,ExecutionReceiptCursor::take_tag,ExecutionReceiptCursor::fixed,ExecutionReceiptCursor::u32,ExecutionReceiptCursor::u64,ExecutionReceiptCursor::u128,ExecutionReceiptCursor::text,ExecutionReceiptCursor::is_finished",
    "schema_constants=GEMM_EXECUTION_RECEIPT_IDENTITY_VERSION,GEMM_EXECUTION_RECEIPT_DOMAIN,GEMM_EXECUTION_RECEIPT_MAGIC,MAX_GEMM_EXECUTION_RECEIPT_BYTES,MAX_GEMM_EXECUTION_RECEIPT_PANELS,MAX_GEMM_EXECUTION_RECEIPT_STRING_BYTES,EXEC_TAG_DOMAIN,EXEC_TAG_VERSION,EXEC_TAG_DECLARED_RUN,EXEC_TAG_COMPLETED_TILES,EXEC_TAG_TOTAL_TILES,EXEC_TAG_MEMORY,EXEC_TAG_MEMORY_LIMIT,EXEC_TAG_MEMORY_STAGING,EXEC_TAG_MEMORY_B_PACK,EXEC_TAG_MEMORY_BAND_METADATA,EXEC_TAG_MEMORY_POOL_RUN,EXEC_TAG_MEMORY_ARENA_PER_WORKER,EXEC_TAG_MEMORY_ACTIVE_WORKERS,EXEC_TAG_MEMORY_ARENA,EXEC_TAG_MEMORY_REQUESTED,EXEC_TAG_PANELS,EXEC_TAG_PANEL,EXEC_TAG_PANEL_KERNEL,EXEC_TAG_PANEL_MODE,EXEC_TAG_PANEL_DECLARED_RUN,EXEC_TAG_PANEL_COMPLETED,EXEC_TAG_PANEL_TOTAL,EXEC_TAG_END,crates/fs-la/src/gemm.rs#GEMM_PANEL_RUN_DOMAIN,crates/fs-blake3/src/lib.rs#IV,crates/fs-blake3/src/lib.rs#MSG_PERMUTATION,crates/fs-blake3/src/lib.rs#BLOCK_LEN,crates/fs-blake3/src/lib.rs#CHUNK_LEN,crates/fs-blake3/src/lib.rs#CHUNK_START,crates/fs-blake3/src/lib.rs#CHUNK_END,crates/fs-blake3/src/lib.rs#PARENT,crates/fs-blake3/src/lib.rs#ROOT,crates/fs-blake3/src/lib.rs#DERIVE_KEY_CONTEXT,crates/fs-blake3/src/lib.rs#DERIVE_KEY_MATERIAL,crates/fs-blake3/src/lib.rs#MAX_DEPTH",
    "schema_functions=GemmExecutionReceipt::from_report,GemmExecutionReceipt::from_canonical_bytes,GemmExecutionReceipt::is_complete,GemmMemoryReceipt::from,execution_receipt_encoded_len,ExecutionReceiptCursor::take_tag,ExecutionReceiptCursor::take_exact,ExecutionReceiptCursor::u32,ExecutionReceiptCursor::u64,ExecutionReceiptCursor::u128,ExecutionReceiptCursor::text,crates/fs-la/src/gemm.rs#gemm_panel_run_id,crates/fs-blake3/src/lib.rs#hash_domain,crates/fs-blake3/src/lib.rs#g,crates/fs-blake3/src/lib.rs#round,crates/fs-blake3/src/lib.rs#permute,crates/fs-blake3/src/lib.rs#compress,crates/fs-blake3/src/lib.rs#words_from_block,crates/fs-blake3/src/lib.rs#first_8_words,crates/fs-blake3/src/lib.rs#Output::chaining_value,crates/fs-blake3/src/lib.rs#Output::root_hash,crates/fs-blake3/src/lib.rs#parent_output,crates/fs-blake3/src/lib.rs#ChunkState::new,crates/fs-blake3/src/lib.rs#ChunkState::len,crates/fs-blake3/src/lib.rs#ChunkState::start_flag,crates/fs-blake3/src/lib.rs#ChunkState::update,crates/fs-blake3/src/lib.rs#ChunkState::output,crates/fs-blake3/src/lib.rs#Blake3::new_internal,crates/fs-blake3/src/lib.rs#Blake3::push_stack,crates/fs-blake3/src/lib.rs#Blake3::pop_stack,crates/fs-blake3/src/lib.rs#Blake3::add_chunk_chaining_value,crates/fs-blake3/src/lib.rs#Blake3::update,crates/fs-blake3/src/lib.rs#Blake3::finalize",
    "schema_dependencies=none",
    "digest=fs-blake3",
    "encoding=typed-binary",
    "sources=GemmExecutionReceipt,GemmMemoryReceipt,GemmPanelReceipt",
    "source_fields=GemmExecutionReceipt.declared_run:semantic,GemmExecutionReceipt.completed_tiles:semantic,GemmExecutionReceipt.total_tiles:semantic,GemmExecutionReceipt.memory:derived:nested-memory-fields-classified-separately,GemmExecutionReceipt.panels:semantic,GemmMemoryReceipt.limit_bytes:semantic,GemmMemoryReceipt.staging_bytes:semantic,GemmMemoryReceipt.b_pack_bytes:semantic,GemmMemoryReceipt.band_metadata_bytes:semantic,GemmMemoryReceipt.pool_run_bytes:semantic,GemmMemoryReceipt.arena_bytes_per_worker:semantic,GemmMemoryReceipt.active_arena_workers:semantic,GemmMemoryReceipt.arena_bytes:semantic,GemmMemoryReceipt.requested_bytes:semantic,GemmPanelReceipt.kernel:semantic,GemmPanelReceipt.mode:semantic,GemmPanelReceipt.declared_run:semantic,GemmPanelReceipt.completed:semantic,GemmPanelReceipt.total:semantic",
    "source_bindings=GemmExecutionReceipt.declared_run>declared-run,GemmExecutionReceipt.completed_tiles>completed-tiles,GemmExecutionReceipt.total_tiles>total-tiles,GemmExecutionReceipt.panels>panel-count+panel-order,GemmMemoryReceipt.limit_bytes>memory-limit-bytes,GemmMemoryReceipt.staging_bytes>memory-staging-bytes,GemmMemoryReceipt.b_pack_bytes>memory-b-pack-bytes,GemmMemoryReceipt.band_metadata_bytes>memory-band-metadata-bytes,GemmMemoryReceipt.pool_run_bytes>memory-pool-run-bytes,GemmMemoryReceipt.arena_bytes_per_worker>memory-arena-bytes-per-worker,GemmMemoryReceipt.active_arena_workers>memory-active-arena-workers,GemmMemoryReceipt.arena_bytes>memory-arena-bytes,GemmMemoryReceipt.requested_bytes>memory-requested-bytes,GemmPanelReceipt.kernel>panel-kernel,GemmPanelReceipt.mode>panel-mode,GemmPanelReceipt.declared_run>panel-declared-run,GemmPanelReceipt.completed>panel-completed,GemmPanelReceipt.total>panel-total",
    "external_semantic_fields=artifact-domain,identity-version,canonical-field-order",
    "semantic_fields=artifact-domain,identity-version,canonical-field-order,declared-run,completed-tiles,total-tiles,memory-limit-bytes,memory-staging-bytes,memory-b-pack-bytes,memory-band-metadata-bytes,memory-pool-run-bytes,memory-arena-bytes-per-worker,memory-active-arena-workers,memory-arena-bytes,memory-requested-bytes,panel-count,panel-order,panel-kernel,panel-mode,panel-declared-run,panel-completed,panel-total",
    "excluded_fields=schedule-steals:observed-scheduling-only,schedule-cross-ccd-steals:observed-scheduling-only,schedule-cancel-latencies:observed-scheduling-only,schedule-worker-distribution:observed-scheduling-only,schedule-memory-peak:observed-allocation-only,schedule-memory-refused:failure-path-only",
    "consumers=GemmDispatch::execution_receipt,GemmExecutionReceipt::is_complete,fs-roofline::KernelExecutionBinding,fs-roofline::execution_path_shape_eq",
    "mutations=artifact-domain:crates/fs-session/src/gemm_tune.rs#gemm_execution_receipt_identity_fields_move_independently,identity-version:crates/fs-session/src/gemm_tune.rs#gemm_execution_receipt_identity_fields_move_independently,canonical-field-order:crates/fs-session/src/gemm_tune.rs#gemm_execution_receipt_identity_fields_move_independently,declared-run:crates/fs-session/src/gemm_tune.rs#gemm_execution_receipt_identity_fields_move_independently,completed-tiles:crates/fs-session/src/gemm_tune.rs#gemm_execution_receipt_identity_fields_move_independently,total-tiles:crates/fs-session/src/gemm_tune.rs#gemm_execution_receipt_identity_fields_move_independently,memory-limit-bytes:crates/fs-session/src/gemm_tune.rs#gemm_execution_receipt_identity_fields_move_independently,memory-staging-bytes:crates/fs-session/src/gemm_tune.rs#gemm_execution_receipt_identity_fields_move_independently,memory-b-pack-bytes:crates/fs-session/src/gemm_tune.rs#gemm_execution_receipt_identity_fields_move_independently,memory-band-metadata-bytes:crates/fs-session/src/gemm_tune.rs#gemm_execution_receipt_identity_fields_move_independently,memory-pool-run-bytes:crates/fs-session/src/gemm_tune.rs#gemm_execution_receipt_identity_fields_move_independently,memory-arena-bytes-per-worker:crates/fs-session/src/gemm_tune.rs#gemm_execution_receipt_identity_fields_move_independently,memory-active-arena-workers:crates/fs-session/src/gemm_tune.rs#gemm_execution_receipt_identity_fields_move_independently,memory-arena-bytes:crates/fs-session/src/gemm_tune.rs#gemm_execution_receipt_identity_fields_move_independently,memory-requested-bytes:crates/fs-session/src/gemm_tune.rs#gemm_execution_receipt_identity_fields_move_independently,panel-count:crates/fs-session/src/gemm_tune.rs#gemm_execution_receipt_identity_fields_move_independently,panel-order:crates/fs-session/src/gemm_tune.rs#gemm_execution_receipt_identity_fields_move_independently,panel-kernel:crates/fs-session/src/gemm_tune.rs#gemm_execution_receipt_identity_fields_move_independently,panel-mode:crates/fs-session/src/gemm_tune.rs#gemm_execution_receipt_identity_fields_move_independently,panel-declared-run:crates/fs-session/src/gemm_tune.rs#gemm_execution_receipt_identity_fields_move_independently,panel-completed:crates/fs-session/src/gemm_tune.rs#gemm_execution_receipt_identity_fields_move_independently,panel-total:crates/fs-session/src/gemm_tune.rs#gemm_execution_receipt_identity_fields_move_independently",
    "nonsemantic_mutations=schedule-steals:crates/fs-session/src/gemm_tune.rs#execution_receipt_excludes_schedule_measurements,schedule-cross-ccd-steals:crates/fs-session/src/gemm_tune.rs#execution_receipt_excludes_schedule_measurements,schedule-cancel-latencies:crates/fs-session/src/gemm_tune.rs#execution_receipt_excludes_schedule_measurements,schedule-worker-distribution:crates/fs-session/src/gemm_tune.rs#execution_receipt_excludes_schedule_measurements,schedule-memory-peak:crates/fs-session/src/gemm_tune.rs#execution_receipt_excludes_schedule_measurements,schedule-memory-refused:crates/fs-session/src/gemm_tune.rs#execution_receipt_excludes_schedule_measurements",
    "field_guard=classify_gemm_execution_receipt_identity_fields",
    "transport_guard=GemmExecutionReceipt::from_canonical_bytes",
    "version_guard=crates/fs-session/src/gemm_tune.rs#gemm_execution_receipt_versions_fail_closed",
    "coupling_surface=fs-session:gemm-execution-receipt",
];

/// Build/dependency identity available to a root before it admits or publishes
/// GEMM tune evidence.
///
/// An operator-observed receipt can be retained as an exact artifact. It is not
/// promoted to independently verified graph evidence by this projection.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GemmTuneBuildEvidence {
    /// Exact generated-code fingerprint already bound into every tune key.
    pub build_fingerprint: &'static str,
    /// Dependency evidence trust class.
    pub graph_class: GemmGraphEvidenceClass,
    /// Fingerprint-bound graph class identity.
    pub graph_class_identity: &'static str,
    /// Exact canonical receipt suitable for artifact retention, when supplied.
    pub dependency_receipt: Option<&'static str>,
    /// Domain-separated digest of `dependency_receipt`, when supplied.
    pub dependency_receipt_digest: Option<&'static str>,
}

/// Evidence a root must inspect before citing or retaining this binary's GEMM
/// tuning results.
#[must_use]
pub const fn gemm_tune_build_evidence() -> GemmTuneBuildEvidence {
    let graph = fs_la::gemm_graph_evidence();
    GemmTuneBuildEvidence {
        build_fingerprint: fs_la::GEMM_BUILD_FINGERPRINT,
        graph_class: graph.class,
        graph_class_identity: graph.class_identity,
        dependency_receipt: graph.receipt,
        dependency_receipt_digest: graph.receipt_digest,
    }
}

fn push_json_string(out: &mut String, value: &str) {
    use core::fmt::Write as _;

    out.push('"');
    for ch in value.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if c.is_control() => {
                let _ = write!(out, "\\u{:04x}", u32::from(c));
            }
            c => out.push(c),
        }
    }
    out.push('"');
}

/// Cursor for the exact, whitespace-free JSON subset emitted by the tune-row
/// receipt writer. It validates embedded JSON without normalizing any retained
/// bytes, so a successful parse can be checked as a writer/parser fixed point.
struct ExactJsonCursor<'a> {
    input: &'a str,
    offset: usize,
}

impl ExactJsonCursor<'_> {
    fn take(&mut self, expected: &str) -> Option<()> {
        self.input
            .get(self.offset..)?
            .starts_with(expected)
            .then(|| self.offset += expected.len())
    }

    fn is_finished(&self) -> bool {
        self.offset == self.input.len()
    }

    fn canonical_string(&mut self) -> Option<String> {
        let start = self.offset;
        self.take("\"")?;
        let mut value = String::new();
        loop {
            let rest = self.input.get(self.offset..)?;
            let ch = rest.chars().next()?;
            match ch {
                '"' => {
                    self.offset += 1;
                    break;
                }
                '\\' => {
                    self.offset += 1;
                    let escape = self.input.get(self.offset..)?.chars().next()?;
                    self.offset += escape.len_utf8();
                    match escape {
                        '"' => value.push('"'),
                        '\\' => value.push('\\'),
                        'n' => value.push('\n'),
                        'r' => value.push('\r'),
                        't' => value.push('\t'),
                        'u' => {
                            let first = self.take_hex_quad()?;
                            let scalar = if (0xd800..=0xdbff).contains(&first) {
                                self.take("\\u")?;
                                let second = self.take_hex_quad()?;
                                if !(0xdc00..=0xdfff).contains(&second) {
                                    return None;
                                }
                                0x1_0000
                                    + ((u32::from(first) - 0xd800) << 10)
                                    + (u32::from(second) - 0xdc00)
                            } else if (0xdc00..=0xdfff).contains(&first) {
                                return None;
                            } else {
                                u32::from(first)
                            };
                            value.push(char::from_u32(scalar)?);
                        }
                        _ => return None,
                    }
                }
                c if c.is_control() => return None,
                c => {
                    value.push(c);
                    self.offset += c.len_utf8();
                }
            }
            if value.len() > MAX_GEMM_TUNE_ROW_RECEIPT_STRING_BYTES {
                return None;
            }
        }
        let consumed = self.input.get(start..self.offset)?;
        let mut canonical = String::new();
        push_json_string(&mut canonical, &value);
        (canonical == consumed).then_some(value)
    }

    fn take_hex_quad(&mut self) -> Option<u16> {
        let hex = self.input.get(self.offset..self.offset.checked_add(4)?)?;
        if !hex
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
        {
            return None;
        }
        self.offset += 4;
        u16::from_str_radix(hex, 16).ok()
    }

    fn canonical_u64(&mut self) -> Option<u64> {
        let value = self.canonical_u128()?;
        u64::try_from(value).ok()
    }

    fn canonical_u128(&mut self) -> Option<u128> {
        let start = self.offset;
        while self
            .input
            .as_bytes()
            .get(self.offset)
            .is_some_and(u8::is_ascii_digit)
        {
            self.offset += 1;
        }
        let digits = self.input.get(start..self.offset)?;
        if digits.is_empty() || (digits.len() > 1 && digits.starts_with('0')) {
            return None;
        }
        let value = digits.parse::<u128>().ok()?;
        (value.to_string() == digits).then_some(value)
    }

    fn canonical_value(&mut self) -> Option<&str> {
        let start = self.offset;
        self.parse_value(0)?;
        self.input.get(start..self.offset)
    }

    fn parse_value(&mut self, depth: usize) -> Option<()> {
        if depth > MAX_GEMM_TUNE_ROW_RECEIPT_JSON_DEPTH {
            return None;
        }
        match *self.input.as_bytes().get(self.offset)? {
            b'"' => self.canonical_string().map(|_| ()),
            b'{' => {
                self.offset += 1;
                if self.input.as_bytes().get(self.offset) == Some(&b'}') {
                    self.offset += 1;
                    return Some(());
                }
                let mut keys = std::collections::BTreeSet::new();
                loop {
                    if !keys.insert(self.canonical_string()?) {
                        return None;
                    }
                    self.take(":")?;
                    self.parse_value(depth + 1)?;
                    match self.input.as_bytes().get(self.offset)? {
                        b',' => self.offset += 1,
                        b'}' => {
                            self.offset += 1;
                            return Some(());
                        }
                        _ => return None,
                    }
                }
            }
            b'[' => {
                self.offset += 1;
                if self.input.as_bytes().get(self.offset) == Some(&b']') {
                    self.offset += 1;
                    return Some(());
                }
                loop {
                    self.parse_value(depth + 1)?;
                    match self.input.as_bytes().get(self.offset)? {
                        b',' => self.offset += 1,
                        b']' => {
                            self.offset += 1;
                            return Some(());
                        }
                        _ => return None,
                    }
                }
            }
            b't' => self.take("true"),
            b'f' => self.take("false"),
            b'n' => self.take("null"),
            b'-' | b'0'..=b'9' => self.parse_number(),
            _ => None,
        }
    }

    fn parse_number(&mut self) -> Option<()> {
        let start = self.offset;
        let negative = self.input.as_bytes().get(self.offset) == Some(&b'-');
        if negative {
            self.offset += 1;
        }
        match self.input.as_bytes().get(self.offset)? {
            b'0' => {
                self.offset += 1;
                if negative {
                    return None;
                }
                if self
                    .input
                    .as_bytes()
                    .get(self.offset)
                    .is_some_and(u8::is_ascii_digit)
                {
                    return None;
                }
            }
            b'1'..=b'9' => {
                self.offset += 1;
                while self
                    .input
                    .as_bytes()
                    .get(self.offset)
                    .is_some_and(u8::is_ascii_digit)
                {
                    self.offset += 1;
                }
            }
            _ => return None,
        }
        (self.offset > start).then_some(())
    }
}

/// A structured autotune-loop failure. Every variant fails closed: sweep
/// failures record no row and nothing dispatches under unvalidated blocking.
/// A cancellation during the final dispatch may retain the already validated
/// measured row, but records no successful decision and does not commit `C`.
#[derive(Debug)]
pub enum GemmTuneError {
    /// The cancel gate was requested. Compute may have completed in private
    /// staging, but the caller's output was not committed.
    Cancelled {
        /// The caller-visible envelope in force.
        limit_bytes: u64,
        /// Largest session-owned logical reservation concurrency reached.
        peak_used_bytes: u128,
        /// Drained numerical-run report when cancellation was returned by
        /// fs-la. `None` means the gate was observed between dispatch calls;
        /// earlier completed probes may still contribute to the peak.
        report: Option<Box<fs_la::GemmRunReport>>,
    },
    /// Tuner-side refusal (invalid pin, evidence, or adoption).
    Tune(TuneError),
    /// Ledger cache I/O failed (the loop does not guess around storage).
    Ledger(String),
    /// Two sweep candidates produced different output bits: the
    /// bit-neutrality contract is broken and NO plan may be selected.
    BitDrift {
        /// Canonical params of the candidate that diverged.
        candidate: String,
        /// One-based repeat whose exact output bits diverged.
        repeat: usize,
    },
    /// The GEMM path refused at its memory boundary (wf9.15): the plan
    /// exceeded the caller envelope or an allocator declined a reservation.
    /// Output is not committed; the retained report may contain drained
    /// private progress from panels that completed before refusal.
    MemoryRefused {
        /// Which reservation was refused.
        what: &'static str,
        /// Bytes the refused reservation asked for.
        requested_bytes: u128,
        /// The envelope in force.
        limit_bytes: u64,
        /// Largest session-owned logical reservation concurrency reached.
        peak_used_bytes: u128,
        /// Drained numerical-run report when refusal occurred inside fs-la.
        report: Option<Box<fs_la::GemmRunReport>>,
    },
    /// Checked arithmetic could not represent the session or fs-la memory plan.
    MemoryPlanOverflow {
        /// Component whose arithmetic overflowed.
        what: &'static str,
        /// The envelope in force.
        limit_bytes: u64,
    },
    /// The TilePool contained a tile panic or executor invariant refusal.
    Executor {
        /// Structured fs-exec outcome with logical tile provenance.
        error: fs_exec::RunError,
        /// The caller-visible envelope in force.
        limit_bytes: u64,
        /// Largest session-owned logical reservation concurrency reached.
        peak_used_bytes: u128,
        /// Drained numerical-run report, including memory and tile progress.
        report: Box<fs_la::GemmRunReport>,
    },
}

impl core::fmt::Display for GemmTuneError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Cancelled {
                limit_bytes,
                peak_used_bytes,
                report,
            } => write!(
                f,
                "gemm work cancelled after {}/{} compute tiles under a {limit_bytes}-byte envelope \
                 after reaching {peak_used_bytes} logical bytes; output not committed",
                report.as_ref().map_or(0, |run| run.completed_tiles),
                report.as_ref().map_or(0, |run| run.total_tiles),
            ),
            Self::Tune(e) => write!(f, "gemm autotune: {e}"),
            Self::Ledger(detail) => write!(f, "gemm autotune ledger cache: {detail}"),
            Self::MemoryRefused {
                what,
                requested_bytes,
                limit_bytes,
                peak_used_bytes,
                report,
            } => write!(
                f,
                "gemm autotune memory refused at {what}: {requested_bytes} bytes against a \
                 {limit_bytes}-byte envelope after reaching {peak_used_bytes} logical bytes and \
                 {}/{} compute tiles; output not committed",
                report.as_ref().map_or(0, |run| run.completed_tiles),
                report.as_ref().map_or(0, |run| run.total_tiles),
            ),
            Self::MemoryPlanOverflow { what, limit_bytes } => write!(
                f,
                "gemm autotune memory-plan arithmetic overflowed at {what} under the \
                 {limit_bytes}-byte envelope; output not touched"
            ),
            Self::BitDrift { candidate, repeat } => write!(
                f,
                "gemm autotune: candidate {candidate} repeat {repeat} broke the MC/NC bit-neutrality contract"
            ),
            Self::Executor {
                error,
                limit_bytes,
                peak_used_bytes,
                report,
            } => write!(
                f,
                "gemm executor failed after {}/{} compute tiles under a {limit_bytes}-byte \
                 envelope after reaching {peak_used_bytes} logical bytes: {error}",
                report.completed_tiles, report.total_tiles,
            ),
        }
    }
}

impl core::error::Error for GemmTuneError {}

impl From<TuneError> for GemmTuneError {
    fn from(e: TuneError) -> Self {
        Self::Tune(e)
    }
}

impl From<fs_la::GemmCancelled> for GemmTuneError {
    fn from(cancelled: fs_la::GemmCancelled) -> Self {
        let limit_bytes = cancelled.report.memory.limit_bytes;
        let peak_used_bytes = cancelled.report.memory.peak_used_bytes;
        Self::Cancelled {
            limit_bytes,
            peak_used_bytes,
            report: Some(cancelled.report),
        }
    }
}

impl From<fs_la::GemmRunError> for GemmTuneError {
    fn from(error: fs_la::GemmRunError) -> Self {
        match error {
            fs_la::GemmRunError::Cancelled(cancelled) => Self::from(cancelled),
            fs_la::GemmRunError::Executor { error, report } => {
                let limit_bytes = report.memory.limit_bytes;
                let peak_used_bytes = report.memory.peak_used_bytes;
                Self::Executor {
                    error,
                    limit_bytes,
                    peak_used_bytes,
                    report,
                }
            }
            fs_la::GemmRunError::MemoryRefused {
                what,
                requested_bytes,
                limit_bytes,
                report,
            } => Self::MemoryRefused {
                what,
                requested_bytes,
                limit_bytes,
                peak_used_bytes: report.memory.peak_used_bytes,
                report: Some(report),
            },
            fs_la::GemmRunError::MemoryPlanOverflow { what, limit_bytes } => {
                Self::MemoryPlanOverflow { what, limit_bytes }
            }
        }
    }
}

fn gemm_error_with_session_memory(
    error: fs_la::GemmRunError,
    envelope: fs_la::GemmMemoryEnvelope,
    session_bytes: u128,
) -> GemmTuneError {
    match GemmTuneError::from(error) {
        GemmTuneError::MemoryRefused {
            what,
            requested_bytes,
            peak_used_bytes,
            report,
            ..
        } => match session_bytes.checked_add(peak_used_bytes) {
            Some(peak_used_bytes) => GemmTuneError::MemoryRefused {
                what,
                requested_bytes,
                limit_bytes: envelope.limit_bytes,
                peak_used_bytes,
                report,
            },
            None => GemmTuneError::MemoryPlanOverflow {
                what: "session-plus-gemm-peak",
                limit_bytes: envelope.limit_bytes,
            },
        },
        GemmTuneError::Cancelled {
            peak_used_bytes,
            report,
            ..
        } => match session_bytes.checked_add(peak_used_bytes) {
            Some(peak_used_bytes) => GemmTuneError::Cancelled {
                limit_bytes: envelope.limit_bytes,
                peak_used_bytes,
                report,
            },
            None => GemmTuneError::MemoryPlanOverflow {
                what: "session-plus-gemm-peak",
                limit_bytes: envelope.limit_bytes,
            },
        },
        GemmTuneError::Executor {
            error,
            peak_used_bytes,
            report,
            ..
        } => match session_bytes.checked_add(peak_used_bytes) {
            Some(peak_used_bytes) => GemmTuneError::Executor {
                error,
                limit_bytes: envelope.limit_bytes,
                peak_used_bytes,
                report,
            },
            None => GemmTuneError::MemoryPlanOverflow {
                what: "session-plus-gemm-peak",
                limit_bytes: envelope.limit_bytes,
            },
        },
        other => other,
    }
}

fn cancelled_before_compute(envelope: fs_la::GemmMemoryEnvelope) -> GemmTuneError {
    GemmTuneError::Cancelled {
        limit_bytes: envelope.limit_bytes,
        peak_used_bytes: 0,
        report: None,
    }
}

fn cancelled_with_live_probe_memory(
    envelope: fs_la::GemmMemoryEnvelope,
    session_bytes: u128,
    numerical_peak: u128,
) -> GemmTuneError {
    let Some(peak_used_bytes) = session_bytes.checked_add(numerical_peak) else {
        return GemmTuneError::MemoryPlanOverflow {
            what: "session-plus-gemm-peak",
            limit_bytes: envelope.limit_bytes,
        };
    };
    GemmTuneError::Cancelled {
        limit_bytes: envelope.limit_bytes,
        peak_used_bytes,
        report: None,
    }
}

/// The receipt for one autotuned dispatch: what ran, under which plan,
/// and where the plan came from. A study records this; replay pins it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GemmDispatch {
    /// Exact scoped kernel key: numerical version plus complete execution
    /// identity.
    pub kernel: String,
    /// Shape class the plan was resolved for.
    pub shape_class: String,
    /// The MC/NC plan that dispatched.
    pub plan: GemmBlockPlan,
    /// Plan provenance (pinned / tuned / cold-start).
    pub source: TuneSource,
    /// True when this call ran the measurement sweep (cold cache).
    pub swept: bool,
    /// Sealed newly measured row. Read-only cache users can retain this
    /// process-locally and publish it only after their enclosing run passes
    /// admission. `None` when no sweep ran.
    pub new_tune_row: Option<ValidatedGemmTuneRow>,
    /// Sealed row adopted or measured during this call. Callers that need a
    /// citable execution receipt retain this identity across later warm-cache
    /// dispatches. `None` when this call reused an already local row or bypassed
    /// tuning.
    pub validated_tune_row: Option<ValidatedGemmTuneRow>,
    /// Final production execution receipt. Its `pool_runs` prove the selected
    /// plan traversed the caller's TilePool rather than a detached thread path.
    pub run: fs_la::GemmRunReport,
}

impl GemmDispatch {
    /// Deterministic execution facts suitable for replay and evidence binding.
    /// Scheduling measurements (steals, worker distribution, and cancellation
    /// latency) deliberately remain outside this identity.
    #[must_use]
    pub fn execution_receipt(&self) -> GemmExecutionReceipt {
        GemmExecutionReceipt::from_report(&self.run)
    }
}

/// Stable facts for one TilePool traversal of an NC/KC panel.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GemmPanelReceipt {
    /// Stable TileKernel name.
    pub kernel: String,
    /// Deterministic/fast execution mode recorded by the pool.
    pub mode: String,
    /// Deterministic NC/KC panel ordinal used as the fs-exec declared run.
    pub declared_run: u64,
    /// Logical M-band tiles completed.
    pub completed: u64,
    /// Logical M-band tiles planned.
    pub total: u64,
}

/// Deterministic production-path receipt for a GEMM dispatch.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GemmExecutionReceipt {
    /// Caller-ledgered identity of the final production dispatch.
    pub declared_run: u64,
    /// Completed bounded GEMM microtiles.
    pub completed_tiles: usize,
    /// Total bounded GEMM microtiles.
    pub total_tiles: usize,
    /// Declared logical-memory plan. Schedule-observed peak and refusal bytes
    /// remain in the full run report and are excluded from replay identity.
    pub memory: GemmMemoryReceipt,
    /// Ordered NC/KC panel traversals through the caller's TilePool.
    pub panels: Vec<GemmPanelReceipt>,
}

/// Identity-stable projection of [`fs_la::GemmMemoryReport`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GemmMemoryReceipt {
    /// Caller-declared memory ceiling.
    pub limit_bytes: u64,
    /// Transactional C staging bytes.
    pub staging_bytes: u128,
    /// Shared B-pack bytes.
    pub b_pack_bytes: u128,
    /// M-band metadata bytes.
    pub band_metadata_bytes: u128,
    /// fs-la panel-receipt vector bytes.
    pub pool_run_bytes: u128,
    /// Fresh arena reservation per active worker.
    pub arena_bytes_per_worker: u64,
    /// Planned maximum active arena workers.
    pub active_arena_workers: usize,
    /// Planned active arena bytes.
    pub arena_bytes: u128,
    /// Checked fs-la-owned plan total.
    pub requested_bytes: u128,
}

/// Refusal from the bounded execution-receipt codec.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GemmExecutionReceiptCodecError {
    detail: &'static str,
}

impl core::fmt::Display for GemmExecutionReceiptCodecError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str(self.detail)
    }
}

impl std::error::Error for GemmExecutionReceiptCodecError {}

fn execution_receipt_codec_error(detail: &'static str) -> GemmExecutionReceiptCodecError {
    GemmExecutionReceiptCodecError { detail }
}

fn checked_receipt_len_add(total: &mut usize, add: usize) -> Option<()> {
    *total = total.checked_add(add)?;
    (*total <= MAX_GEMM_EXECUTION_RECEIPT_BYTES).then_some(())
}

fn execution_receipt_encoded_len(receipt: &GemmExecutionReceipt, domain: &str) -> Option<usize> {
    if receipt.panels.len() > MAX_GEMM_EXECUTION_RECEIPT_PANELS
        || domain.len() > MAX_GEMM_EXECUTION_RECEIPT_STRING_BYTES
    {
        return None;
    }
    let mut total = GEMM_EXECUTION_RECEIPT_MAGIC.len();
    checked_receipt_len_add(&mut total, 1 + 4 + domain.len())?;
    checked_receipt_len_add(&mut total, 1 + 4)?;
    checked_receipt_len_add(&mut total, 3 * (1 + 8))?;
    checked_receipt_len_add(&mut total, 1)?;
    checked_receipt_len_add(&mut total, 1 + 8)?;
    checked_receipt_len_add(&mut total, 4 * (1 + 16))?;
    checked_receipt_len_add(&mut total, 1 + 8)?;
    checked_receipt_len_add(&mut total, 1 + 8)?;
    checked_receipt_len_add(&mut total, 2 * (1 + 16))?;
    checked_receipt_len_add(&mut total, 1 + 8)?;
    for panel in &receipt.panels {
        if panel.kernel.len() > MAX_GEMM_EXECUTION_RECEIPT_STRING_BYTES
            || panel.mode.len() > MAX_GEMM_EXECUTION_RECEIPT_STRING_BYTES
        {
            return None;
        }
        checked_receipt_len_add(&mut total, 1)?;
        checked_receipt_len_add(&mut total, 1 + 4 + panel.kernel.len())?;
        checked_receipt_len_add(&mut total, 1 + 4 + panel.mode.len())?;
        checked_receipt_len_add(&mut total, 3 * (1 + 8))?;
    }
    checked_receipt_len_add(&mut total, 1)?;
    Some(total)
}

fn push_execution_u32(out: &mut Vec<u8>, tag: u8, value: u32) {
    out.push(tag);
    out.extend_from_slice(&value.to_le_bytes());
}

fn push_execution_u64(out: &mut Vec<u8>, tag: u8, value: u64) {
    out.push(tag);
    out.extend_from_slice(&value.to_le_bytes());
}

fn push_execution_u128(out: &mut Vec<u8>, tag: u8, value: u128) {
    out.push(tag);
    out.extend_from_slice(&value.to_le_bytes());
}

fn push_execution_text(
    out: &mut Vec<u8>,
    tag: u8,
    value: &str,
) -> Result<(), GemmExecutionReceiptCodecError> {
    let len = u32::try_from(value.len())
        .ok()
        .filter(|&len| {
            usize::try_from(len).is_ok_and(|len| len <= MAX_GEMM_EXECUTION_RECEIPT_STRING_BYTES)
        })
        .ok_or_else(|| execution_receipt_codec_error("execution receipt string exceeds its cap"))?;
    out.push(tag);
    out.extend_from_slice(&len.to_le_bytes());
    out.extend_from_slice(value.as_bytes());
    Ok(())
}

struct ExecutionReceiptCursor<'a> {
    bytes: &'a [u8],
    offset: usize,
}

impl ExecutionReceiptCursor<'_> {
    fn take_exact(&mut self, expected: &[u8]) -> Option<()> {
        let end = self.offset.checked_add(expected.len())?;
        (self.bytes.get(self.offset..end)? == expected).then(|| self.offset = end)
    }

    fn take_tag(&mut self, expected: u8) -> Option<()> {
        (self.bytes.get(self.offset).copied()? == expected).then(|| self.offset += 1)
    }

    fn fixed<const N: usize>(&mut self) -> Option<[u8; N]> {
        let end = self.offset.checked_add(N)?;
        let value = self.bytes.get(self.offset..end)?.try_into().ok()?;
        self.offset = end;
        Some(value)
    }

    fn u32(&mut self, tag: u8) -> Option<u32> {
        self.take_tag(tag)?;
        Some(u32::from_le_bytes(self.fixed()?))
    }

    fn u64(&mut self, tag: u8) -> Option<u64> {
        self.take_tag(tag)?;
        Some(u64::from_le_bytes(self.fixed()?))
    }

    fn u128(&mut self, tag: u8) -> Option<u128> {
        self.take_tag(tag)?;
        Some(u128::from_le_bytes(self.fixed()?))
    }

    fn text(&mut self, tag: u8) -> Option<String> {
        self.take_tag(tag)?;
        let len = usize::try_from(u32::from_le_bytes(self.fixed()?)).ok()?;
        if len > MAX_GEMM_EXECUTION_RECEIPT_STRING_BYTES {
            return None;
        }
        let end = self.offset.checked_add(len)?;
        let value = core::str::from_utf8(self.bytes.get(self.offset..end)?)
            .ok()?
            .to_string();
        self.offset = end;
        Some(value)
    }

    fn is_finished(&self) -> bool {
        self.offset == self.bytes.len()
    }
}

impl From<fs_la::GemmMemoryReport> for GemmMemoryReceipt {
    fn from(report: fs_la::GemmMemoryReport) -> Self {
        Self {
            limit_bytes: report.limit_bytes,
            staging_bytes: report.staging_bytes,
            b_pack_bytes: report.b_pack_bytes,
            band_metadata_bytes: report.band_metadata_bytes,
            pool_run_bytes: report.pool_run_bytes,
            arena_bytes_per_worker: report.arena_bytes_per_worker,
            active_arena_workers: report.active_arena_workers,
            arena_bytes: report.arena_bytes,
            requested_bytes: report.requested_bytes,
        }
    }
}

impl GemmExecutionReceipt {
    /// Project a successful numerical run report onto identity-stable fields.
    /// Error paths retain their full report in [`GemmTuneError`] instead.
    #[must_use]
    pub fn from_report(report: &fs_la::GemmRunReport) -> Self {
        Self {
            declared_run: report.declared_run.0,
            completed_tiles: report.completed_tiles,
            total_tiles: report.total_tiles,
            memory: report.memory.into(),
            panels: report
                .pool_runs
                .iter()
                .map(|panel| GemmPanelReceipt {
                    kernel: panel.kernel.to_string(),
                    mode: panel.mode.to_string(),
                    declared_run: panel.declared_run.0,
                    completed: panel.completed,
                    total: panel.total,
                })
                .collect(),
        }
    }

    fn canonical_bytes_with_schema(
        &self,
        domain: &str,
        identity_version: u32,
    ) -> Result<Vec<u8>, GemmExecutionReceiptCodecError> {
        let encoded_len = execution_receipt_encoded_len(self, domain).ok_or_else(|| {
            execution_receipt_codec_error("execution receipt exceeds its canonical transport cap")
        })?;
        let declared_run = self.declared_run;
        let completed_tiles = u64::try_from(self.completed_tiles).map_err(|_| {
            execution_receipt_codec_error("completed tile count does not fit canonical u64")
        })?;
        let total_tiles = u64::try_from(self.total_tiles).map_err(|_| {
            execution_receipt_codec_error("total tile count does not fit canonical u64")
        })?;
        let active_arena_workers =
            u64::try_from(self.memory.active_arena_workers).map_err(|_| {
                execution_receipt_codec_error(
                    "active arena worker count does not fit canonical u64",
                )
            })?;
        let panel_count = u64::try_from(self.panels.len())
            .map_err(|_| execution_receipt_codec_error("panel count does not fit canonical u64"))?;

        let mut out = Vec::new();
        out.try_reserve_exact(encoded_len).map_err(|_| {
            execution_receipt_codec_error("execution receipt allocation was refused")
        })?;
        out.extend_from_slice(GEMM_EXECUTION_RECEIPT_MAGIC);
        push_execution_text(&mut out, EXEC_TAG_DOMAIN, domain)?;
        push_execution_u32(&mut out, EXEC_TAG_VERSION, identity_version);
        push_execution_u64(&mut out, EXEC_TAG_DECLARED_RUN, declared_run);
        push_execution_u64(&mut out, EXEC_TAG_COMPLETED_TILES, completed_tiles);
        push_execution_u64(&mut out, EXEC_TAG_TOTAL_TILES, total_tiles);
        out.push(EXEC_TAG_MEMORY);
        push_execution_u64(&mut out, EXEC_TAG_MEMORY_LIMIT, self.memory.limit_bytes);
        push_execution_u128(&mut out, EXEC_TAG_MEMORY_STAGING, self.memory.staging_bytes);
        push_execution_u128(&mut out, EXEC_TAG_MEMORY_B_PACK, self.memory.b_pack_bytes);
        push_execution_u128(
            &mut out,
            EXEC_TAG_MEMORY_BAND_METADATA,
            self.memory.band_metadata_bytes,
        );
        push_execution_u128(
            &mut out,
            EXEC_TAG_MEMORY_POOL_RUN,
            self.memory.pool_run_bytes,
        );
        push_execution_u64(
            &mut out,
            EXEC_TAG_MEMORY_ARENA_PER_WORKER,
            self.memory.arena_bytes_per_worker,
        );
        push_execution_u64(
            &mut out,
            EXEC_TAG_MEMORY_ACTIVE_WORKERS,
            active_arena_workers,
        );
        push_execution_u128(&mut out, EXEC_TAG_MEMORY_ARENA, self.memory.arena_bytes);
        push_execution_u128(
            &mut out,
            EXEC_TAG_MEMORY_REQUESTED,
            self.memory.requested_bytes,
        );
        push_execution_u64(&mut out, EXEC_TAG_PANELS, panel_count);
        for panel in &self.panels {
            out.push(EXEC_TAG_PANEL);
            push_execution_text(&mut out, EXEC_TAG_PANEL_KERNEL, &panel.kernel)?;
            push_execution_text(&mut out, EXEC_TAG_PANEL_MODE, &panel.mode)?;
            push_execution_u64(&mut out, EXEC_TAG_PANEL_DECLARED_RUN, panel.declared_run);
            push_execution_u64(&mut out, EXEC_TAG_PANEL_COMPLETED, panel.completed);
            push_execution_u64(&mut out, EXEC_TAG_PANEL_TOTAL, panel.total);
        }
        out.push(EXEC_TAG_END);
        if out.len() != encoded_len {
            return Err(execution_receipt_codec_error(
                "execution receipt encoded length disagrees with its checked plan",
            ));
        }
        Ok(out)
    }

    /// Exact tagged binary transport for deterministic execution facts.
    ///
    /// The frame carries its domain and version, uses fixed-width little-endian
    /// integers, length-prefixes every string and collection, and preserves
    /// panel order. It is identical across ISAs for equal receipt values.
    ///
    /// # Errors
    ///
    /// Returns an error if a platform-sized count cannot fit the fixed-width
    /// schema, a string or collection exceeds its cap, or allocation fails.
    #[must_use]
    pub fn canonical_bytes(&self) -> Result<Vec<u8>, GemmExecutionReceiptCodecError> {
        self.canonical_bytes_with_schema(
            GEMM_EXECUTION_RECEIPT_DOMAIN,
            GEMM_EXECUTION_RECEIPT_IDENTITY_VERSION,
        )
    }

    /// Parse one exact current tagged binary transport.
    ///
    /// Stale domains/versions, missing or reordered tags, oversized lengths,
    /// invalid UTF-8, trailing bytes, and any non-fixed-point spelling fail
    /// closed.
    ///
    /// # Errors
    ///
    /// Returns an error unless `bytes` are the exact current bounded transport.
    #[must_use]
    pub fn from_canonical_bytes(bytes: &[u8]) -> Result<Self, GemmExecutionReceiptCodecError> {
        let fail = || {
            execution_receipt_codec_error(
                "execution receipt is not an exact current canonical transport",
            )
        };
        if bytes.len() > MAX_GEMM_EXECUTION_RECEIPT_BYTES {
            return Err(fail());
        }
        let mut parser = ExecutionReceiptCursor { bytes, offset: 0 };
        parser
            .take_exact(GEMM_EXECUTION_RECEIPT_MAGIC)
            .ok_or_else(fail)?;
        if parser.text(EXEC_TAG_DOMAIN).ok_or_else(fail)? != GEMM_EXECUTION_RECEIPT_DOMAIN {
            return Err(fail());
        }
        if parser.u32(EXEC_TAG_VERSION).ok_or_else(fail)? != GEMM_EXECUTION_RECEIPT_IDENTITY_VERSION
        {
            return Err(fail());
        }
        let declared_run = parser.u64(EXEC_TAG_DECLARED_RUN).ok_or_else(fail)?;
        let completed_tiles =
            usize::try_from(parser.u64(EXEC_TAG_COMPLETED_TILES).ok_or_else(fail)?)
                .map_err(|_| fail())?;
        let total_tiles = usize::try_from(parser.u64(EXEC_TAG_TOTAL_TILES).ok_or_else(fail)?)
            .map_err(|_| fail())?;
        parser.take_tag(EXEC_TAG_MEMORY).ok_or_else(fail)?;
        let memory = GemmMemoryReceipt {
            limit_bytes: parser.u64(EXEC_TAG_MEMORY_LIMIT).ok_or_else(fail)?,
            staging_bytes: parser.u128(EXEC_TAG_MEMORY_STAGING).ok_or_else(fail)?,
            b_pack_bytes: parser.u128(EXEC_TAG_MEMORY_B_PACK).ok_or_else(fail)?,
            band_metadata_bytes: parser
                .u128(EXEC_TAG_MEMORY_BAND_METADATA)
                .ok_or_else(fail)?,
            pool_run_bytes: parser.u128(EXEC_TAG_MEMORY_POOL_RUN).ok_or_else(fail)?,
            arena_bytes_per_worker: parser
                .u64(EXEC_TAG_MEMORY_ARENA_PER_WORKER)
                .ok_or_else(fail)?,
            active_arena_workers: usize::try_from(
                parser
                    .u64(EXEC_TAG_MEMORY_ACTIVE_WORKERS)
                    .ok_or_else(fail)?,
            )
            .map_err(|_| fail())?,
            arena_bytes: parser.u128(EXEC_TAG_MEMORY_ARENA).ok_or_else(fail)?,
            requested_bytes: parser.u128(EXEC_TAG_MEMORY_REQUESTED).ok_or_else(fail)?,
        };
        let panel_count = usize::try_from(parser.u64(EXEC_TAG_PANELS).ok_or_else(fail)?)
            .ok()
            .filter(|&count| count <= MAX_GEMM_EXECUTION_RECEIPT_PANELS)
            .ok_or_else(fail)?;
        let mut panels = Vec::new();
        panels.try_reserve_exact(panel_count).map_err(|_| fail())?;
        for _ in 0..panel_count {
            parser.take_tag(EXEC_TAG_PANEL).ok_or_else(fail)?;
            panels.push(GemmPanelReceipt {
                kernel: parser.text(EXEC_TAG_PANEL_KERNEL).ok_or_else(fail)?,
                mode: parser.text(EXEC_TAG_PANEL_MODE).ok_or_else(fail)?,
                declared_run: parser.u64(EXEC_TAG_PANEL_DECLARED_RUN).ok_or_else(fail)?,
                completed: parser.u64(EXEC_TAG_PANEL_COMPLETED).ok_or_else(fail)?,
                total: parser.u64(EXEC_TAG_PANEL_TOTAL).ok_or_else(fail)?,
            });
        }
        parser.take_tag(EXEC_TAG_END).ok_or_else(fail)?;
        if !parser.is_finished() {
            return Err(fail());
        }
        let receipt = Self {
            declared_run,
            completed_tiles,
            total_tiles,
            memory,
            panels,
        };
        if receipt.canonical_bytes().map_err(|_| fail())? != bytes {
            return Err(fail());
        }
        Ok(receipt)
    }

    /// Domain-separated digest of [`Self::canonical_bytes`].
    ///
    /// # Errors
    ///
    /// Returns the same transport-construction errors as
    /// [`Self::canonical_bytes`].
    #[must_use]
    pub fn receipt_identity(
        &self,
    ) -> Result<fs_ledger::ContentHash, GemmExecutionReceiptCodecError> {
        let bytes = self.canonical_bytes()?;
        Ok(fs_blake3::hash_domain(
            GEMM_EXECUTION_RECEIPT_DOMAIN,
            &bytes,
        ))
    }

    /// Whether every panel completed and carries the exact child RunId derived
    /// from this receipt's declared operation identity and panel ordinal.
    #[must_use]
    pub fn is_complete(&self) -> bool {
        if self.completed_tiles != self.total_tiles {
            return false;
        }
        if self.total_tiles == 0 {
            return self.panels.is_empty();
        }
        !self.panels.is_empty()
            && self.panels.iter().enumerate().all(|(ordinal, panel)| {
                let Ok(ordinal) = u64::try_from(ordinal) else {
                    return false;
                };
                panel.declared_run
                    == fs_la::gemm_panel_run_id(fs_exec::RunId(self.declared_run), ordinal).0
                    && !panel.kernel.is_empty()
                    && !panel.mode.is_empty()
                    && panel.total > 0
                    && panel.completed == panel.total
            })
    }
}

#[allow(dead_code)]
fn classify_gemm_execution_receipt_identity_fields(
    receipt: &GemmExecutionReceipt,
    memory_source: &GemmMemoryReceipt,
    panel_source: &GemmPanelReceipt,
) {
    let GemmExecutionReceipt {
        declared_run,
        completed_tiles,
        total_tiles,
        memory,
        panels,
    } = receipt;
    let GemmMemoryReceipt {
        limit_bytes,
        staging_bytes,
        b_pack_bytes,
        band_metadata_bytes,
        pool_run_bytes,
        arena_bytes_per_worker,
        active_arena_workers,
        arena_bytes,
        requested_bytes,
    } = memory_source;
    let GemmPanelReceipt {
        kernel,
        mode,
        declared_run: panel_declared_run,
        completed,
        total,
    } = panel_source;
    let _ = (
        declared_run,
        completed_tiles,
        total_tiles,
        memory,
        panels,
        limit_bytes,
        staging_bytes,
        b_pack_bytes,
        band_metadata_bytes,
        pool_run_bytes,
        arena_bytes_per_worker,
        active_arena_workers,
        arena_bytes,
        requested_bytes,
        kernel,
        mode,
        panel_declared_run,
        completed,
        total,
    );
}

/// Explicit access policy for the durable GEMM tune cache.
///
/// A read-only caller may adopt a previously admitted row but cannot publish a
/// new measurement during speculative or not-yet-admitted work. Newly measured
/// rows remain available through [`GemmDispatch::new_tune_row`].
#[derive(Clone, Copy)]
pub enum GemmTuneCache<'a> {
    /// Do not read or write durable tuning state.
    Disabled,
    /// Adopt validated rows, but never write the ledger.
    ReadOnly(&'a Ledger),
    /// Adopt validated rows and persist a newly measured row before local
    /// installation.
    ReadWrite(&'a Ledger),
}

impl<'a> GemmTuneCache<'a> {
    fn reader(self) -> Option<&'a Ledger> {
        match self {
            Self::Disabled => None,
            Self::ReadOnly(ledger) | Self::ReadWrite(ledger) => Some(ledger),
        }
    }
}

/// A validated tune row that can be published after a wider admission gate.
///
/// Fields are private: callers can neither forge nor alter the scoped kernel,
/// shape, machine fingerprint, selected parameters, or measured evidence.
/// Instances are created only after fs-exec validates the complete row.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ValidatedGemmTuneRow {
    kernel: String,
    shape_class: String,
    machine: [u8; 8],
    params: String,
    measured: String,
    memory_limit_bytes: u64,
    probe_buffer_bytes: u128,
}

fn parse_validated_gemm_tune_row_receipt(receipt_json: &str) -> Option<ValidatedGemmTuneRow> {
    if receipt_json.len() > MAX_GEMM_TUNE_ROW_RECEIPT_BYTES {
        return None;
    }
    let mut parser = ExactJsonCursor {
        input: receipt_json,
        offset: 0,
    };
    parser.take("{\"kernel\":")?;
    let kernel = parser.canonical_string()?;
    parser.take(",\"shape_class\":")?;
    let shape_class = parser.canonical_string()?;
    parser.take(",\"machine\":")?;
    let machine_hex = parser.canonical_string()?;
    if machine_hex.len() != 16
        || !machine_hex
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
    {
        return None;
    }
    let machine_value = u64::from_str_radix(&machine_hex, 16).ok()?;
    let machine = machine_value.to_le_bytes();
    parser.take(",\"params\":")?;
    let params_value = parser.canonical_string()?;
    let mut params = String::new();
    push_json_string(&mut params, &params_value);
    parser.take(",\"measured\":")?;
    let measured = parser.canonical_value()?.to_string();
    if !measured.starts_with('{')
        || measured.len() > tune_metadata_plan::MEASURED_BYTES_CAP
        || params.len() > tune_metadata_plan::PARAMS_BYTES_CAP
        || kernel.len() > MAX_GEMM_TUNE_ROW_RECEIPT_STRING_BYTES
        || shape_class.len() > MAX_GEMM_TUNE_ROW_RECEIPT_STRING_BYTES
    {
        return None;
    }
    let measured_row = TuneRow::from_canonical_json(&measured).ok()?;
    if measured_row.kernel() != kernel.as_str()
        || measured_row.shape_class() != shape_class.as_str()
        || measured_row.machine() != machine_value
        || measured_row.params() != params_value.as_str()
    {
        return None;
    }
    parser.take(",\"memory_limit_bytes\":")?;
    let memory_limit_bytes = parser.canonical_u64()?;
    parser.take(",\"probe_buffer_bytes\":")?;
    let probe_buffer_bytes = parser.canonical_u128()?;
    parser.take(",\"metadata_plan\":{\"schema\":")?;
    if parser.canonical_string()? != GEMM_TUNE_METADATA_PLAN_SCHEMA {
        return None;
    }
    parser.take(",\"requested_bytes\":")?;
    if parser.canonical_u128()? != gemm_tune_metadata_plan_bytes() {
        return None;
    }
    parser.take("}}")?;
    if !parser.is_finished() {
        return None;
    }
    let row = ValidatedGemmTuneRow {
        kernel,
        shape_class,
        machine,
        params,
        measured,
        memory_limit_bytes,
        probe_buffer_bytes,
    };
    (row.receipt_json() == receipt_json).then_some(row)
}

impl ValidatedGemmTuneRow {
    fn from_prepared(prepared: &PreparedGemmRow, machine: u64) -> Result<Self, GemmTuneError> {
        let execution = prepared.key().execution();
        let probe_buffer_bytes = probe_buffer_bytes_for_dims(execution.probe_dims()).ok_or(
            GemmTuneError::MemoryPlanOverflow {
                what: "tune-probe-buffers",
                limit_bytes: execution.memory_limit_bytes(),
            },
        )?;
        let params = prepared.params_json();
        let measured = prepared.row_json();
        // Session tune-metadata plan caps are ENFORCED at seal time (bead
        // wf9.15.1): a sealed-row string beyond its documented cap means the
        // schema grew past the plan, which must fail loudly here rather
        // than silently undercount.
        assert!(
            params.len() <= tune_metadata_plan::PARAMS_BYTES_CAP,
            "sealed-row params JSON exceeds the tune-metadata plan cap (schema drift)"
        );
        assert!(
            measured.len() <= tune_metadata_plan::MEASURED_BYTES_CAP,
            "sealed-row measured JSON exceeds the tune-metadata plan cap (schema drift)"
        );
        Ok(Self {
            kernel: prepared.key().kernel().to_string(),
            shape_class: prepared.key().shape_class().to_string(),
            machine: machine.to_le_bytes(),
            params,
            measured,
            memory_limit_bytes: execution.memory_limit_bytes(),
            probe_buffer_bytes,
        })
    }

    /// Canonical JSON preimage for this exact tune-table tuple.
    #[must_use]
    pub fn receipt_json(&self) -> String {
        use core::fmt::Write as _;

        let mut out = String::new();
        out.push_str("{\"kernel\":");
        push_json_string(&mut out, &self.kernel);
        out.push_str(",\"shape_class\":");
        push_json_string(&mut out, &self.shape_class);
        let _ = write!(
            out,
            ",\"machine\":\"{:016x}\",\"params\":",
            u64::from_le_bytes(self.machine)
        );
        out.push_str(&self.params);
        out.push_str(",\"measured\":");
        out.push_str(&self.measured);
        let _ = write!(
            out,
            ",\"memory_limit_bytes\":{},\"probe_buffer_bytes\":{}",
            self.memory_limit_bytes, self.probe_buffer_bytes
        );
        // Bind the session tune-metadata plan into the row receipt (bead
        // wf9.15.1). The plan is a pure constant of the sweep lattice and
        // schema caps, so a freshly measured row and the same row adopted
        // later derive the identical fragment — receipt identity stays
        // stable across both paths. (This field's introduction rotates
        // pre-plan row identities once; see CONTRACT.)
        out.push_str(",\"metadata_plan\":");
        out.push_str(&tune_metadata_plan::receipt_fragment());
        out.push('}');
        out
    }

    /// Domain-separated identity of [`Self::receipt_json`]. It is stable for a
    /// freshly measured row and the same row adopted later.
    #[must_use]
    pub fn receipt_identity(&self) -> fs_ledger::ContentHash {
        fs_blake3::hash_domain(GEMM_TUNE_ROW_RECEIPT_DOMAIN, self.receipt_json().as_bytes())
    }

    /// Admit retained tune-row receipt bytes under explicit schema metadata.
    ///
    /// The v2 receipt preimage stays byte-for-byte compatible with the shipped
    /// writer. Replay carries the domain/version beside those bytes; this path
    /// refuses stale metadata, re-adopts the embedded fs-exec tune row through
    /// its exact semantic parser, cross-checks every duplicated key field,
    /// requires the current metadata-plan fragment, rebuilds the sealed private
    /// tuple, and accepts only a byte-identical writer/parser fixed point.
    #[must_use]
    pub fn admit_receipt_json(
        identity_domain: &str,
        identity_version: u32,
        receipt_json: &str,
    ) -> Option<fs_ledger::ContentHash> {
        if identity_domain != GEMM_TUNE_ROW_RECEIPT_DOMAIN
            || identity_version != GEMM_TUNE_ROW_RECEIPT_IDENTITY_VERSION
        {
            return None;
        }
        let row = parse_validated_gemm_tune_row_receipt(receipt_json)?;
        Some(row.receipt_identity())
    }

    /// Whether this sealed row is the exact evidence behind one dispatched
    /// decision.
    #[must_use]
    pub fn matches_decision(
        &self,
        scoped_kernel: &str,
        shape_class: &str,
        machine: u64,
        canonical_plan: &str,
    ) -> bool {
        self.kernel == scoped_kernel
            && self.shape_class == shape_class
            && self.machine == machine.to_le_bytes()
            && self.params == format!("\"{canonical_plan}\"")
    }

    /// Whether a ledger query returned this exact sealed tuple.
    #[must_use]
    pub fn matches_ledger_row(&self, row: &fs_ledger::TuneRow) -> bool {
        self.kernel == row.kernel
            && self.shape_class == row.shape_class
            && self.machine.as_slice() == row.machine
            && self.params == row.params
            && self.measured == row.measured
    }

    /// Publish this already validated row without replacing a different row,
    /// preserving the caller's ledger transaction when one is active.
    ///
    /// # Errors
    /// Propagates the original ledger diagnostic.
    pub fn publish_to_ledger(&self, ledger: &Ledger) -> Result<(), fs_ledger::LedgerError> {
        self.publish_if_absent_or_identical(ledger)
    }

    /// Insert this row into an evidence ledger without replacing a different
    /// tune decision already stored under the same key. An identical row is an
    /// idempotent success; a conflict fails closed.
    ///
    /// # Errors
    /// Propagates ledger failures or returns [`fs_ledger::LedgerError::Invalid`]
    /// when the destination key contains a different tuple.
    pub fn publish_if_absent_or_identical(
        &self,
        ledger: &Ledger,
    ) -> Result<(), fs_ledger::LedgerError> {
        ledger.tune_put_if_absent(
            &self.kernel,
            &self.shape_class,
            &self.machine,
            &self.params,
            &self.measured,
        )?;
        let stored = ledger
            .tune_get(&self.kernel, &self.shape_class, &self.machine)?
            .ok_or_else(|| fs_ledger::LedgerError::Invalid {
                field: "tune".to_string(),
                problem: "insert-if-absent returned without a stored tune row".to_string(),
            })?;
        if self.matches_ledger_row(&stored) {
            Ok(())
        } else {
            Err(fs_ledger::LedgerError::Invalid {
                field: "tune".to_string(),
                problem: format!(
                    "refusing to replace a conflicting tune row for kernel {:?}, shape {:?}",
                    self.kernel, self.shape_class
                ),
            })
        }
    }

    /// Persist this already validated row to a durable tune ledger.
    ///
    /// # Errors
    /// Returns [`GemmTuneError::Ledger`] when the ledger refuses the write.
    pub fn persist(&self, ledger: &Ledger) -> Result<(), GemmTuneError> {
        self.publish_to_ledger(ledger)
            .map_err(|error| GemmTuneError::Ledger(error.to_string()))
    }

    /// Install this sealed row as the current mutable cache decision.
    ///
    /// Unlike [`ValidatedGemmTuneRow::publish_if_absent_or_identical`], this
    /// method deliberately replaces a stale or malformed row under the same
    /// cache key. The `tune` table is a dispatch cache, not an append-only
    /// evidence history; citable benchmark publication uses the insert-only
    /// method and content-addressed ledger artifacts instead.
    ///
    /// # Errors
    /// Returns [`GemmTuneError::Ledger`] when the upsert or exact read-back
    /// verification fails.
    pub fn replace_cache_row(&self, ledger: &Ledger) -> Result<(), GemmTuneError> {
        ledger
            .tune_put(
                &self.kernel,
                &self.shape_class,
                &self.machine,
                &self.params,
                &self.measured,
            )
            .map_err(|error| GemmTuneError::Ledger(error.to_string()))?;
        let stored = ledger
            .tune_get(&self.kernel, &self.shape_class, &self.machine)
            .map_err(|error| GemmTuneError::Ledger(error.to_string()))?
            .ok_or_else(|| {
                GemmTuneError::Ledger(
                    "cache upsert returned without a stored GEMM tune row".to_string(),
                )
            })?;
        if self.matches_ledger_row(&stored) {
            Ok(())
        } else {
            Err(GemmTuneError::Ledger(format!(
                "cache read-back disagrees with the sealed GEMM row for kernel {:?}, shape {:?}",
                self.kernel, self.shape_class
            )))
        }
    }
}

#[allow(dead_code)]
fn classify_validated_gemm_tune_row_identity_fields(row: &ValidatedGemmTuneRow) {
    let ValidatedGemmTuneRow {
        kernel,
        shape_class,
        machine,
        params,
        measured,
        memory_limit_bytes,
        probe_buffer_bytes,
    } = row;
    let _ = (
        kernel,
        shape_class,
        machine,
        params,
        measured,
        memory_limit_bytes,
        probe_buffer_bytes,
    );
}

/// The kernel key for this build's GEMM accumulation contract.
#[must_use]
pub fn gemm_kernel_key() -> String {
    format!(
        "{GEMM_KERNEL_PREFIX}{}",
        fs_la::gemm::GEMM_BIT_SEMANTICS_VERSION
    )
}

/// Bucket one extent to its shape-class quantum (next power of two,
/// clamped to [8, 65536]).
fn bucket(extent: usize) -> usize {
    extent.clamp(8, 65_536).next_power_of_two()
}

/// The shape class for an (m, n, k) problem: power-of-two buckets. Exact
/// measured probe dims remain in [`GemmTuneKey`], so a bucket never erases
/// the context that produced a row.
#[must_use]
pub fn gemm_shape_class(m: usize, n: usize, k: usize) -> String {
    format!("m{}-n{}-k{}", bucket(m), bucket(n), bucket(k))
}

fn probe_dims(m: usize, n: usize, k: usize) -> [usize; 3] {
    [
        m.clamp(1, PROBE_MK_DIM_CAP),
        n.clamp(1, PROBE_N_DIM_CAP),
        k.clamp(1, PROBE_MK_DIM_CAP),
    ]
}

fn probe_buffer_bytes_for_dims([m, n, k]: [u64; 3]) -> Option<u128> {
    let m = u128::from(m);
    let n = u128::from(n);
    let k = u128::from(k);
    let elements = m
        .checked_mul(k)?
        .checked_add(k.checked_mul(n)?)?
        .checked_add(m.checked_mul(n)?.checked_mul(2)?)?;
    elements.checked_mul(core::mem::size_of::<u64>() as u128)
}

/// Construct the exact persistent tuning identity for this invocation.
/// Studies normally replay the recorded decision key directly; exposing this
/// constructor also lets admission and diagnostics explain why two calls do
/// or do not share evidence.
///
/// # Errors
/// [`GemmTuneError::Tune`] if a dimension or implementation identity cannot
/// be represented canonically.
pub fn gemm_tune_key(
    threads: usize,
    m: usize,
    n: usize,
    k: usize,
) -> Result<GemmTuneKey, GemmTuneError> {
    gemm_tune_key_budgeted(threads, m, n, k, fs_la::GemmMemoryEnvelope::UNBOUNDED)
}

/// Construct the persistent tuning identity under an explicit memory envelope.
/// Otherwise-identical calls with different ceilings cannot share rows or pins.
///
/// # Errors
/// As [`gemm_tune_key`].
pub fn gemm_tune_key_budgeted(
    threads: usize,
    m: usize,
    n: usize,
    k: usize,
    envelope: fs_la::GemmMemoryEnvelope,
) -> Result<GemmTuneKey, GemmTuneError> {
    let pool = TilePool::for_host(threads, SESSION_GEMM_POOL_SEED);
    gemm_tune_key_with_pool_budgeted(&pool, m, n, k, envelope)
}

/// Construct the persistent tuning identity from the TilePool that will
/// actually execute the measured and selected plans.
///
/// # Errors
/// [`GemmTuneError::Tune`] if a pool dimension or identity component cannot
/// be represented canonically.
pub fn gemm_tune_key_with_pool(
    pool: &TilePool,
    m: usize,
    n: usize,
    k: usize,
) -> Result<GemmTuneKey, GemmTuneError> {
    gemm_tune_key_with_pool_budgeted(pool, m, n, k, fs_la::GemmMemoryEnvelope::UNBOUNDED)
}

/// Construct the persistent tuning identity from the executing pool and an
/// explicit memory envelope.
///
/// # Errors
/// As [`gemm_tune_key_with_pool`].
pub fn gemm_tune_key_with_pool_budgeted(
    pool: &TilePool,
    m: usize,
    n: usize,
    k: usize,
    envelope: fs_la::GemmMemoryEnvelope,
) -> Result<GemmTuneKey, GemmTuneError> {
    gemm_tune_key_for_execution(
        pool.workers(),
        pool.workers(),
        envelope.limit_bytes,
        &pool.placement_identity(),
        m,
        n,
        k,
    )
}

fn gemm_tune_key_for_execution(
    requested_threads: usize,
    thread_budget: usize,
    memory_limit_bytes: u64,
    placement: &str,
    m: usize,
    n: usize,
    k: usize,
) -> Result<GemmTuneKey, GemmTuneError> {
    gemm_tune_key_for_execution_schema(
        requested_threads,
        thread_budget,
        memory_limit_bytes,
        placement,
        m,
        n,
        k,
        GEMM_TUNER_SCHEMA_VERSION,
    )
}

#[allow(clippy::too_many_arguments)]
fn gemm_tune_key_for_execution_schema(
    requested_threads: usize,
    thread_budget: usize,
    memory_limit_bytes: u64,
    placement: &str,
    m: usize,
    n: usize,
    k: usize,
    tuner_schema: u32,
) -> Result<GemmTuneKey, GemmTuneError> {
    debug_assert!(tuner_schema > 0);
    let implementation = format!(
        "fs-la-{}-gemm-v{}-fs-session-tuner-v{tuner_schema}",
        fs_la::VERSION,
        fs_la::GEMM_IMPLEMENTATION_VERSION
    );
    let execution = GemmExecutionIdentity::new(
        requested_threads,
        thread_budget,
        memory_limit_bytes,
        probe_dims(m, n, k),
        fs_la::gemm_execution_tier(),
        placement,
        implementation,
        fs_la::gemm_build_identity(),
    )?;
    Ok(GemmTuneKey::new(
        gemm_kernel_key(),
        gemm_shape_class(m, n, k),
        execution,
    )?)
}

#[track_caller]
fn checked_product(label: &str, lhs: usize, rhs: usize) -> usize {
    lhs.checked_mul(rhs)
        .unwrap_or_else(|| panic!("{label} extent overflow: {lhs} * {rhs}"))
}

fn try_filled_buffer<T: Copy>(
    len: usize,
    value: T,
    what: &'static str,
    envelope: fs_la::GemmMemoryEnvelope,
    peak_used_bytes: u128,
) -> Result<Vec<T>, GemmTuneError> {
    let requested_bytes = (len as u128)
        .checked_mul(core::mem::size_of::<T>() as u128)
        .ok_or(GemmTuneError::MemoryPlanOverflow {
            what,
            limit_bytes: envelope.limit_bytes,
        })?;
    let mut values = Vec::new();
    values
        .try_reserve_exact(len)
        .map_err(|_| GemmTuneError::MemoryRefused {
            what,
            requested_bytes,
            limit_bytes: envelope.limit_bytes,
            peak_used_bytes,
            report: None,
        })?;
    values.resize(len, value);
    Ok(values)
}

/// Mirror fs-la's public contiguous-slice precondition before consulting or
/// mutating tuning state. fs-la validates again at the execution boundary;
/// this ordering is the session-level no-phantom-row guarantee.
#[track_caller]
fn assert_contiguous_shapes(m: usize, n: usize, k: usize, a: &[f64], b: &[f64], c: &[f64]) {
    let a_len = checked_product("a", m, k);
    let b_len = checked_product("b", k, n);
    let c_len = checked_product("c", m, n);
    assert_eq!(a.len(), a_len, "a must be m*k = {a_len}");
    assert_eq!(b.len(), b_len, "b must be k*n = {b_len}");
    assert_eq!(c.len(), c_len, "c must be m*n = {c_len}");
}

/// Deterministic probe fill (splitmix64 bits folded to [-0.5, 0.5)):
/// integer-only, so probe inputs are bit-identical on every ISA.
fn probe_fill(buf: &mut [f64], salt: u64) {
    for (i, slot) in buf.iter_mut().enumerate() {
        let mut z = (i as u64)
            .wrapping_add(salt)
            .wrapping_add(0x9E37_79B9_7F4A_7C15);
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^= z >> 31;
        // 53 mantissa bits → [0, 1), then center.
        *slot = (z >> 11) as f64 / 9_007_199_254_740_992.0 - 0.5;
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct SweepCandidate {
    plan: GemmBlockPlan,
    effective_mc: usize,
    effective_nc: usize,
}

#[derive(Debug)]
struct SweepResult {
    winner: GemmBlockPlan,
    evidence: TuneEvidence,
}

/// Build the lattice the kernel will ACTUALLY execute. Nominal plans that
/// collapse to the same clamped `(mc, nc)` pair are measured only once.
fn effective_sweep_candidates(
    pm: usize,
    pn: usize,
    envelope: fs_la::GemmMemoryEnvelope,
) -> Result<Vec<SweepCandidate>, GemmTuneError> {
    // Dedup by linear scan (bead wf9.15.1): the lattice is at most
    // CANDIDATE_CAP entries, and BTreeSet node overhead is not honestly
    // accountable under the session tune-memory plan.
    let mut candidates = try_reserved_vec::<SweepCandidate>(
        tune_metadata_plan::CANDIDATE_CAP,
        "tune-metadata-candidates",
        envelope,
    )?;
    for (mc, nc_cap) in SWEEP_MC
        .iter()
        .flat_map(|&mc| SWEEP_NC_CAP.iter().map(move |&nc| (mc, nc)))
    {
        let plan = GemmBlockPlan::new(mc, nc_cap)?;
        // These are the clamps applied by fs-la's packed parallel engine.
        let effective_mc = plan.mc.max(8).min(pm.max(8));
        let effective_nc = pn.min(plan.nc_cap).max(4);
        if !candidates.iter().any(|c: &SweepCandidate| {
            c.effective_mc == effective_mc && c.effective_nc == effective_nc
        }) {
            candidates.push(SweepCandidate {
                plan,
                effective_mc,
                effective_nc,
            });
        }
    }
    Ok(candidates)
}

/// Checked logical peaks for the probe output and exact-bit reference
/// buffers: (first-output peak, both-live peak).
fn probe_output_peaks(
    output_len: usize,
    envelope: fs_la::GemmMemoryEnvelope,
    base_used_bytes: u128,
) -> Result<(u128, u128), GemmTuneError> {
    let output_bytes = (output_len as u128)
        .checked_mul(core::mem::size_of::<u64>() as u128)
        .ok_or(GemmTuneError::MemoryPlanOverflow {
            what: "tune-probe-output",
            limit_bytes: envelope.limit_bytes,
        })?;
    let first_output_peak =
        base_used_bytes
            .checked_add(output_bytes)
            .ok_or(GemmTuneError::MemoryPlanOverflow {
                what: "tune-probe-output-peak",
                limit_bytes: envelope.limit_bytes,
            })?;
    let live_probe_bytes =
        first_output_peak
            .checked_add(output_bytes)
            .ok_or(GemmTuneError::MemoryPlanOverflow {
                what: "tune-probe-reference-peak",
                limit_bytes: envelope.limit_bytes,
            })?;
    Ok((first_output_peak, live_probe_bytes))
}

/// One candidate's evidence observation: an exact-length sample copy
/// (accounted by the plan's per-candidate sample-storage component) under a
/// label whose schema cap is ENFORCED, not assumed (bead wf9.15.1).
fn candidate_observation(
    candidate: &SweepCandidate,
    sample_scratch: &[u64],
    envelope: fs_la::GemmMemoryEnvelope,
) -> Result<TuneObservation, GemmTuneError> {
    let label = candidate.plan.canonical();
    assert!(
        label.len() <= tune_metadata_plan::LABEL_BYTES_CAP,
        "canonical plan label exceeds the tune-metadata plan cap (schema drift)"
    );
    let mut samples_owned =
        try_reserved_vec::<u64>(sample_scratch.len(), "tune-metadata-sample-copy", envelope)?;
    samples_owned.extend_from_slice(sample_scratch);
    Ok(TuneObservation::wall_time(label, samples_owned)?)
}

/// Fallibly pre-reserve one bounded metadata collection (bead wf9.15.1);
/// allocator refusal is a typed diagnostic, never an abort.
fn try_reserved_vec<T>(
    cap: usize,
    what: &'static str,
    envelope: fs_la::GemmMemoryEnvelope,
) -> Result<Vec<T>, GemmTuneError> {
    let mut values = Vec::new();
    values.try_reserve_exact(cap).map_err(|_| {
        let requested_bytes = (cap as u128).saturating_mul(size_of::<T>() as u128);
        GemmTuneError::MemoryRefused {
            what,
            requested_bytes,
            limit_bytes: envelope.limit_bytes,
            peak_used_bytes: 0,
            report: None,
        }
    })?;
    Ok(values)
}

/// Measure candidate executions supplied by `run`. Keeping this core
/// injectable lets the Gauntlet force drift in each repeat and cache faults
/// without adding test behavior to the production GEMM implementation.
fn measure_candidates<R>(
    gate: &CancelGate,
    candidates: &[SweepCandidate],
    output_len: usize,
    envelope: fs_la::GemmMemoryEnvelope,
    base_used_bytes: u128,
    mut run: R,
) -> Result<SweepResult, GemmTuneError>
where
    R: FnMut(&SweepCandidate, &mut [f64]) -> Result<u128, GemmTuneError>,
{
    let (first_output_peak, live_probe_bytes) =
        probe_output_peaks(output_len, envelope, base_used_bytes)?;
    let mut c = try_filled_buffer(
        output_len,
        0.0_f64,
        "tune-probe-c",
        envelope,
        base_used_bytes,
    )?;
    let mut reference_bits = try_filled_buffer(
        output_len,
        0_u64,
        "tune-probe-reference",
        envelope,
        first_output_peak,
    )?;
    let mut observations = try_reserved_vec::<TuneObservation>(
        candidates.len(),
        "tune-metadata-observations",
        envelope,
    )?;
    let mut ranked = try_reserved_vec::<(u64, usize, GemmBlockPlan)>(
        candidates.len(),
        "tune-metadata-ranked",
        envelope,
    )?;
    // One reused sample buffer for the whole sweep (bead wf9.15.1): the
    // per-candidate observation still receives its own exact-length copy,
    // but scratch growth happens once, fallibly.
    let mut sample_scratch =
        try_reserved_vec::<u64>(SWEEP_SAMPLES, "tune-metadata-samples", envelope)?;
    let mut reference_initialized = false;
    let mut numerical_peak = 0_u128;
    for (index, candidate) in candidates.iter().enumerate() {
        if gate.is_requested() {
            return Err(cancelled_with_live_probe_memory(
                envelope,
                live_probe_bytes,
                numerical_peak,
            ));
        }
        sample_scratch.clear();
        for repeat in 1..=SWEEP_SAMPLES {
            if gate.is_requested() {
                return Err(cancelled_with_live_probe_memory(
                    envelope,
                    live_probe_bytes,
                    numerical_peak,
                ));
            }
            c.fill(0.0);
            let t0 = std::time::Instant::now();
            numerical_peak = numerical_peak.max(run(candidate, &mut c)?);
            let ns = u64::try_from(t0.elapsed().as_nanos()).unwrap_or(u64::MAX);
            sample_scratch.push(ns.max(1));
            if gate.is_requested() {
                return Err(cancelled_with_live_probe_memory(
                    envelope,
                    live_probe_bytes,
                    numerical_peak,
                ));
            }

            // Compare every output word directly. A fixed-width digest is
            // not a proof of bit-neutrality and would also hide which repeat
            // drifted. `to_bits` intentionally distinguishes signed zero and
            // every NaN payload.
            if !reference_initialized {
                for (dst, value) in reference_bits.iter_mut().zip(&c) {
                    *dst = value.to_bits();
                }
                reference_initialized = true;
            } else if !reference_bits
                .iter()
                .zip(&c)
                .all(|(&expected, value)| expected == value.to_bits())
            {
                return Err(GemmTuneError::BitDrift {
                    candidate: candidate.plan.canonical(),
                    repeat,
                });
            }
        }
        let best = sample_scratch.iter().copied().min().unwrap_or(u64::MAX);
        ranked.push((best, index, candidate.plan));
        observations.push(candidate_observation(candidate, &sample_scratch, envelope)?);
    }
    ranked.sort_unstable_by_key(|&(ns, index, _)| (ns, index));
    let winner = ranked
        .first()
        .map(|entry| entry.2)
        .ok_or_else(|| TuneError {
            detail: "the effective GEMM candidate lattice is empty".to_string(),
        })?;
    let evidence = TuneEvidence::ranked_wall_times(observations)?;
    Ok(SweepResult { winner, evidence })
}

/// Run the bounded candidate sweep for one exact probe. This function only
/// measures and validates; its caller persists first and commits the tuner
/// row second so a cache failure cannot leave a phantom in-memory success.
fn run_sweep(
    gate: &CancelGate,
    pool: &TilePool,
    declared_run: fs_exec::RunId,
    m: usize,
    n: usize,
    k: usize,
    envelope: fs_la::GemmMemoryEnvelope,
) -> Result<SweepResult, GemmTuneError> {
    // Probe at the CALLER's dims (capped): the oracle lane showed that
    // probing at the class's power-of-two bucket flips winners — at
    // m = 320 the band count under each mc differs from m = 512, and
    // band balance decides the ranking. The row retains the bucketed shape
    // class, but the exact capped probe is also part of the scoped key so a
    // neighboring caller cannot silently inherit different evidence.
    let [pm, pn, pk] = probe_dims(m, n, k);
    let probe_dims_u64 = [pm, pn, pk]
        .map(|extent| u64::try_from(extent).expect("capped GEMM probe dimensions fit u64"));
    let probe_buffer_bytes =
        probe_buffer_bytes_for_dims(probe_dims_u64).ok_or(GemmTuneError::MemoryPlanOverflow {
            what: "tune-probe-buffers",
            limit_bytes: envelope.limit_bytes,
        })?;
    if probe_buffer_bytes > u128::from(envelope.limit_bytes) {
        return Err(GemmTuneError::MemoryRefused {
            what: "tune-probe-envelope",
            requested_bytes: probe_buffer_bytes,
            limit_bytes: envelope.limit_bytes,
            peak_used_bytes: 0,
            report: None,
        });
    }
    // Session tune-metadata plan (bead wf9.15.1): after the probe buffers
    // clear the envelope on their own, charge the bounded metadata byte
    // plan ON TOP, still before any allocation. A refusal here loses
    // nothing — no fs-la report exists yet — and the plan constant is the
    // same one the sealed row binds into its receipt.
    let metadata_bytes = tune_metadata_plan::requested_bytes();
    let planned_bytes = probe_buffer_bytes.checked_add(metadata_bytes).ok_or(
        GemmTuneError::MemoryPlanOverflow {
            what: "tune-metadata-plan",
            limit_bytes: envelope.limit_bytes,
        },
    )?;
    if planned_bytes > u128::from(envelope.limit_bytes) {
        return Err(GemmTuneError::MemoryRefused {
            what: "tune-metadata-plan",
            requested_bytes: planned_bytes,
            limit_bytes: envelope.limit_bytes,
            peak_used_bytes: 0,
            report: None,
        });
    }
    let child_limit_bytes = if envelope == fs_la::GemmMemoryEnvelope::UNBOUNDED {
        u64::MAX
    } else {
        // The metadata plan's bytes are live for the whole sweep, so the
        // child envelope excludes them alongside the probe buffers.
        u64::try_from(u128::from(envelope.limit_bytes) - planned_bytes)
            .expect("bounded probe preflight leaves a u64 child envelope")
    };
    let child_envelope = fs_la::GemmMemoryEnvelope {
        limit_bytes: child_limit_bytes,
    };

    let a_len = checked_product("tune probe A", pm, pk);
    let b_len = checked_product("tune probe B", pk, pn);
    let a_bytes = (a_len as u128)
        .checked_mul(core::mem::size_of::<f64>() as u128)
        .ok_or(GemmTuneError::MemoryPlanOverflow {
            what: "tune-probe-a",
            limit_bytes: envelope.limit_bytes,
        })?;
    let b_bytes = (b_len as u128)
        .checked_mul(core::mem::size_of::<f64>() as u128)
        .ok_or(GemmTuneError::MemoryPlanOverflow {
            what: "tune-probe-b",
            limit_bytes: envelope.limit_bytes,
        })?;
    let ab_bytes = a_bytes
        .checked_add(b_bytes)
        .ok_or(GemmTuneError::MemoryPlanOverflow {
            what: "tune-probe-a-plus-b",
            limit_bytes: envelope.limit_bytes,
        })?;
    let mut a = try_filled_buffer(a_len, 0.0_f64, "tune-probe-a", envelope, 0)?;
    let mut b = try_filled_buffer(b_len, 0.0_f64, "tune-probe-b", envelope, a_bytes)?;
    probe_fill(&mut a, 0xA);
    probe_fill(&mut b, 0xB);
    let candidates = effective_sweep_candidates(pm, pn, envelope)?;
    let mut sweep_ordinal = 0_u64;
    measure_candidates(
        gate,
        &candidates,
        checked_product("tune probe C", pm, pn),
        envelope,
        ab_bytes,
        |candidate, c| {
            let sweep_run = declared_run.derive(GEMM_SWEEP_RUN_DOMAIN, sweep_ordinal);
            sweep_ordinal = sweep_ordinal.checked_add(1).ok_or_else(|| TuneError {
                detail: "GEMM sweep run ordinal exhausted".to_string(),
            })?;
            fs_la::gemm_f64_parallel_with_pool_budgeted(
                pm,
                pn,
                pk,
                1.0,
                &a,
                &b,
                0.0,
                c,
                pool,
                candidate.effective_mc,
                candidate.effective_nc,
                gate,
                sweep_run,
                child_envelope,
            )
            .map(|report| report.memory.peak_used_bytes)
            .map_err(|error| gemm_error_with_session_memory(error, envelope, probe_buffer_bytes))
        },
    )
}

/// Persist a validated measured row before installing it in the process-local
/// tuner. `persist` is injectable so the failure-atomic boundary is directly
/// testable without corrupting a real ledger connection.
fn install_sweep_row<P>(
    tuner: &mut Tuner,
    key: &GemmTuneKey,
    sweep: SweepResult,
    persist: P,
) -> Result<(GemmBlockPlan, ValidatedGemmTuneRow), GemmTuneError>
where
    P: FnOnce(&ValidatedGemmTuneRow) -> Result<(), GemmTuneError>,
{
    let prepared = tuner.prepare_gemm_row(key, sweep.winner, sweep.evidence)?;
    let validated = ValidatedGemmTuneRow::from_prepared(&prepared, tuner.machine())?;
    persist(&validated)?;
    let winner = sweep.winner;
    tuner.commit_gemm_row(prepared)?;
    Ok((winner, validated))
}

fn adopt_cached_row(
    tuner: &mut Tuner,
    key: &GemmTuneKey,
    params: &str,
    measured: &str,
) -> Result<Option<ValidatedGemmTuneRow>, GemmTuneError> {
    let Ok(prepared) = tuner.prepare_adopt_gemm_row_json(key, measured) else {
        return Ok(None);
    };
    if params != prepared.params_json() {
        return Ok(None);
    }
    let validated = ValidatedGemmTuneRow::from_prepared(&prepared, tuner.machine())?;
    tuner.commit_gemm_row(prepared)?;
    Ok(Some(validated))
}

fn execute_prepared_decision<R, F>(
    tuner: &mut Tuner,
    decision: PreparedGemmDecision,
    run: F,
) -> Result<(GemmBlockPlan, TuneSource, R), GemmTuneError>
where
    F: FnOnce(GemmBlockPlan) -> Result<R, GemmTuneError>,
{
    let plan = decision.plan();
    let source = decision.source();
    let output = run(plan)?;
    // Exclusive access to `tuner` spans prepare -> run -> commit, so no
    // applicable pin/row can change and make this prepared decision stale.
    tuner
        .commit_gemm_decision(decision)
        .expect("exclusive tuner borrow preserves a prepared GEMM decision");
    Ok((plan, source, output))
}

/// The production autotuned f64 GEMM: `c = alpha·a·b + beta·c` through
/// the measure → cache → model → dispatch loop.
///
/// Resolution order after shape and cancellation preflight: a pinned plan
/// dispatches without measurement; else an exact cached row (in the tuner,
/// seeded from the cache when permitted); else the bounded sweep measures,
/// applies the explicit write policy, commits it locally, and dispatches.
/// Serial, small-M, and no-product calls bypass tuning entirely.
///
/// # Errors
/// [`GemmTuneError`] — cancellation, tuner refusals, ledger I/O, or a
/// bit-neutrality violation. On every returned error, `c` retains its exact
/// original bits. Cancellable GEMM computes in private staging, drains its
/// workers, and commits only after its final poll.
///
/// # Panics
/// Inherits fs-la's structured shape panics for mismatched slice
/// lengths.
#[allow(clippy::too_many_arguments)] // BLAS-shape signature + orchestration handles
pub fn gemm_f64_session(
    tuner: &mut Tuner,
    cache: GemmTuneCache<'_>,
    gate: &CancelGate,
    threads: usize,
    m: usize,
    n: usize,
    k: usize,
    alpha: f64,
    a: &[f64],
    b: &[f64],
    beta: f64,
    c: &mut [f64],
) -> Result<GemmDispatch, GemmTuneError> {
    gemm_f64_session_budgeted(
        tuner,
        cache,
        gate,
        threads,
        m,
        n,
        k,
        alpha,
        a,
        b,
        beta,
        c,
        fs_la::GemmMemoryEnvelope::UNBOUNDED,
    )
}

/// As [`gemm_f64_session`], under an explicit memory envelope bound into tune
/// identity and every numerical dispatch.
///
/// # Errors
/// As [`gemm_f64_session`], plus structured memory refusal.
///
/// # Panics
/// As [`gemm_f64_session`].
#[allow(clippy::too_many_arguments)]
pub fn gemm_f64_session_budgeted(
    tuner: &mut Tuner,
    cache: GemmTuneCache<'_>,
    gate: &CancelGate,
    threads: usize,
    m: usize,
    n: usize,
    k: usize,
    alpha: f64,
    a: &[f64],
    b: &[f64],
    beta: f64,
    c: &mut [f64],
    envelope: fs_la::GemmMemoryEnvelope,
) -> Result<GemmDispatch, GemmTuneError> {
    let pool = TilePool::for_host(threads, SESSION_GEMM_POOL_SEED);
    gemm_f64_session_with_pool_budgeted(
        tuner, cache, &pool, gate, m, n, k, alpha, a, b, beta, c, envelope,
    )
}

/// The production autotuned f64 GEMM on a caller-owned, reusable TilePool.
/// The same pool executes every sweep candidate and the selected plan; its
/// normalized worker budget and placement policy are bound into the tune key.
///
/// # Errors
/// As [`gemm_f64_session`], plus a structured executor failure if TilePool
/// contains a tile panic or detects an incomplete traversal.
///
/// # Panics
/// Inherits fs-la's structured shape panics for mismatched slice lengths.
#[allow(clippy::too_many_arguments)]
pub fn gemm_f64_session_with_pool(
    tuner: &mut Tuner,
    cache: GemmTuneCache<'_>,
    pool: &TilePool,
    gate: &CancelGate,
    m: usize,
    n: usize,
    k: usize,
    alpha: f64,
    a: &[f64],
    b: &[f64],
    beta: f64,
    c: &mut [f64],
) -> Result<GemmDispatch, GemmTuneError> {
    gemm_f64_session_with_pool_budgeted(
        tuner,
        cache,
        pool,
        gate,
        m,
        n,
        k,
        alpha,
        a,
        b,
        beta,
        c,
        fs_la::GemmMemoryEnvelope::UNBOUNDED,
    )
}

/// As [`gemm_f64_session_with_pool`], under an explicit memory envelope.
///
/// # Errors
/// As [`gemm_f64_session_with_pool`], plus structured memory refusal.
///
/// # Panics
/// As [`gemm_f64_session_with_pool`].
#[allow(clippy::too_many_arguments)]
pub fn gemm_f64_session_with_pool_budgeted(
    tuner: &mut Tuner,
    cache: GemmTuneCache<'_>,
    pool: &TilePool,
    gate: &CancelGate,
    m: usize,
    n: usize,
    k: usize,
    alpha: f64,
    a: &[f64],
    b: &[f64],
    beta: f64,
    c: &mut [f64],
    envelope: fs_la::GemmMemoryEnvelope,
) -> Result<GemmDispatch, GemmTuneError> {
    gemm_f64_session_with_pool_declared_budgeted(
        tuner,
        cache,
        pool,
        gate,
        fs_exec::RunId::default(),
        m,
        n,
        k,
        alpha,
        a,
        b,
        beta,
        c,
        envelope,
    )
}

/// As [`gemm_f64_session_with_pool`], with the caller-ledgered identity of the
/// final production dispatch. Sweep repetitions receive separate
/// domain-derived children and cannot collide with the final run's tile
/// streams.
///
/// # Errors
/// As [`gemm_f64_session_with_pool`].
///
/// # Panics
/// As [`gemm_f64_session_with_pool`].
#[allow(clippy::too_many_arguments)]
pub fn gemm_f64_session_with_pool_declared(
    tuner: &mut Tuner,
    cache: GemmTuneCache<'_>,
    pool: &TilePool,
    gate: &CancelGate,
    declared_run: fs_exec::RunId,
    m: usize,
    n: usize,
    k: usize,
    alpha: f64,
    a: &[f64],
    b: &[f64],
    beta: f64,
    c: &mut [f64],
) -> Result<GemmDispatch, GemmTuneError> {
    gemm_f64_session_with_pool_declared_budgeted(
        tuner,
        cache,
        pool,
        gate,
        declared_run,
        m,
        n,
        k,
        alpha,
        a,
        b,
        beta,
        c,
        fs_la::GemmMemoryEnvelope::UNBOUNDED,
    )
}

/// As [`gemm_f64_session_with_pool_declared`], under an explicit memory
/// envelope bound into tune identity, sweep admission, and final dispatch.
///
/// # Errors
/// As [`gemm_f64_session_with_pool_declared`], plus structured memory refusal.
///
/// # Panics
/// As [`gemm_f64_session_with_pool_declared`].
#[allow(clippy::too_many_arguments)]
pub fn gemm_f64_session_with_pool_declared_budgeted(
    tuner: &mut Tuner,
    cache: GemmTuneCache<'_>,
    pool: &TilePool,
    gate: &CancelGate,
    declared_run: fs_exec::RunId,
    m: usize,
    n: usize,
    k: usize,
    alpha: f64,
    a: &[f64],
    b: &[f64],
    beta: f64,
    c: &mut [f64],
    envelope: fs_la::GemmMemoryEnvelope,
) -> Result<GemmDispatch, GemmTuneError> {
    // Public slice/extent preconditions are checked before tier resolution,
    // cache reads, sweeps, rows, or decisions. Invalid work cannot poison the
    // tuning state and `c` is still untouched when this panics.
    assert_contiguous_shapes(m, n, k, a, b, c);
    if gate.is_requested() {
        return Err(cancelled_before_compute(envelope));
    }

    let key = gemm_tune_key_with_pool_budgeted(pool, m, n, k, envelope)?;
    let kernel = key.kernel().to_string();
    let shape_class = gemm_shape_class(m, n, k);
    let mut swept = false;
    let mut new_tune_row = None;
    let mut validated_tune_row = None;

    // No product, one-thread, and small-M routes do not have a meaningful
    // production MC/NC choice. Dispatch them cancellation-correctly under the
    // documented cold plan without reading or mutating tune state.
    if !fs_la::gemm_tuning_is_effective(m, n, k, alpha, pool.workers()) {
        let plan = GemmBlockPlan::COLD_START;
        let run = fs_la::gemm_f64_parallel_with_pool_budgeted(
            m,
            n,
            k,
            alpha,
            a,
            b,
            beta,
            c,
            pool,
            plan.mc,
            n.min(plan.nc_cap).max(1),
            gate,
            declared_run,
            envelope,
        )
        .map_err(GemmTuneError::from)?;
        return Ok(GemmDispatch {
            kernel,
            shape_class,
            plan,
            source: TuneSource::ColdStart,
            swept,
            new_tune_row,
            validated_tune_row,
            run,
        });
    }

    if !tuner.has_gemm_pin(&key) && !tuner.has_gemm_row(&key) {
        // Cache tier: try the ledger before measuring. Stale
        // (other-machine) or non-canonical rows are refused by
        // prepare_adopt_gemm_row_json and we fall through to a fresh sweep.
        // The ledger's separate params column must agree byte-for-byte with
        // the validated row body before either is allowed into the tuner.
        if let Some(ledger) = cache.reader() {
            let cached = ledger
                .tune_get(
                    key.kernel(),
                    key.shape_class(),
                    &tuner.machine().to_le_bytes(),
                )
                .map_err(|e| GemmTuneError::Ledger(e.to_string()))?;
            if let Some(row) = cached {
                validated_tune_row = adopt_cached_row(tuner, &key, &row.params, &row.measured)?;
            }
        }
        if validated_tune_row.is_none() {
            let sweep = run_sweep(gate, pool, declared_run, m, n, k, envelope)?;
            swept = true;
            let (_, validated) = install_sweep_row(tuner, &key, sweep, |row| match cache {
                GemmTuneCache::ReadWrite(ledger) => row.replace_cache_row(ledger),
                GemmTuneCache::Disabled | GemmTuneCache::ReadOnly(_) => Ok(()),
            })?;
            new_tune_row = Some(validated.clone());
            validated_tune_row = Some(validated);
        }
    }

    if gate.is_requested() {
        return Err(cancelled_before_compute(envelope));
    }
    let decision = tuner.prepare_gemm_decision(&key);
    let (plan, source, run) = execute_prepared_decision(tuner, decision, |plan| {
        fs_la::gemm_f64_parallel_with_pool_budgeted(
            m,
            n,
            k,
            alpha,
            a,
            b,
            beta,
            c,
            pool,
            plan.mc,
            n.min(plan.nc_cap).max(1),
            gate,
            declared_run,
            envelope,
        )
        .map_err(GemmTuneError::from)
    })?;
    Ok(GemmDispatch {
        kernel,
        shape_class,
        plan,
        source,
        swept,
        new_tune_row,
        validated_tune_row,
        run,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn root_build_evidence_projects_exact_fs_la_material() {
        let evidence = gemm_tune_build_evidence();
        let graph = fs_la::gemm_graph_evidence();
        assert_eq!(evidence.build_fingerprint, fs_la::gemm_build_identity());
        assert_eq!(evidence.graph_class, graph.class);
        assert_eq!(evidence.graph_class_identity, graph.class_identity);
        assert_eq!(evidence.dependency_receipt, graph.receipt);
        assert_eq!(evidence.dependency_receipt_digest, graph.receipt_digest);
        if let (Some(receipt), Some(digest)) = (
            evidence.dependency_receipt,
            evidence.dependency_receipt_digest,
        ) {
            assert_eq!(
                fs_blake3::hash_domain(GEMM_DEPGRAPH_RECEIPT_DOMAIN, receipt.as_bytes())
                    .to_string(),
                digest
            );
        }
        match evidence.graph_class {
            GemmGraphEvidenceClass::OperatorObservedReceipt => {
                assert!(evidence.dependency_receipt.is_some());
                assert!(evidence.dependency_receipt_digest.is_some());
            }
            GemmGraphEvidenceClass::DevelopmentEquivalenceSalt => {
                assert!(evidence.dependency_receipt.is_none());
                assert!(evidence.dependency_receipt_digest.is_none());
            }
        }
    }

    #[test]
    fn tuner_schema_bump_separates_durable_keys() {
        let v1 =
            gemm_tune_key_for_execution_schema(4, 4, u64::MAX, "test-placement", 512, 640, 512, 1)
                .expect("schema-v1 key");
        let v2 =
            gemm_tune_key_for_execution_schema(4, 4, u64::MAX, "test-placement", 512, 640, 512, 2)
                .expect("schema-v2 key");

        assert_ne!(v1.kernel(), v2.kernel());
        assert_eq!(v1.shape_class(), v2.shape_class());
        assert!(v1.kernel().contains("fs-session-tuner-v1"));
        assert!(v2.kernel().contains("fs-session-tuner-v2"));
    }

    fn identity_test_measured_tune_row(
        kernel: &str,
        shape_class: &str,
        machine: u64,
        refresh: u32,
    ) -> String {
        let winner = GemmBlockPlan::new(16, 512).expect("winner plan");
        let runner_up = GemmBlockPlan::new(32, 512).expect("runner-up plan");
        let evidence = TuneEvidence::ranked_wall_times(vec![
            TuneObservation::wall_time(winner.canonical(), vec![10, 11, 12])
                .expect("winner evidence"),
            TuneObservation::wall_time(runner_up.canonical(), vec![20, 21, 22])
                .expect("runner-up evidence"),
        ])
        .expect("ranked evidence");
        TuneRow::new(
            kernel,
            shape_class,
            machine,
            winner.canonical(),
            evidence,
            refresh,
        )
        .expect("valid measured GEMM tune row")
        .to_canonical_json()
        .expect("canonical measured GEMM tune row")
    }

    fn identity_test_validated_tune_row() -> ValidatedGemmTuneRow {
        let memory_limit_bytes = 1 << 30;
        let key = gemm_tune_key_budgeted(
            4,
            320,
            288,
            300,
            fs_la::GemmMemoryEnvelope {
                limit_bytes: memory_limit_bytes,
            },
        )
        .expect("identity-test key");
        let machine = 0x0102_0304_0506_0708_u64;
        let selected = GemmBlockPlan::new(16, 512)
            .expect("selected plan")
            .canonical();
        let mut params = String::new();
        push_json_string(&mut params, &selected);
        let probe_buffer_bytes = probe_buffer_bytes_for_dims(key.execution().probe_dims())
            .expect("identity-test probe plan");
        ValidatedGemmTuneRow {
            kernel: key.kernel().to_string(),
            shape_class: key.shape_class().to_string(),
            machine: machine.to_le_bytes(),
            params,
            measured: identity_test_measured_tune_row(key.kernel(), key.shape_class(), machine, 1),
            memory_limit_bytes,
            probe_buffer_bytes,
        }
    }

    #[test]
    fn gemm_tune_row_receipt_identity_fields_move_independently() {
        let base_row = identity_test_validated_tune_row();
        let base_json = base_row.receipt_json();
        let base_identity = base_row.receipt_identity();
        assert_eq!(
            ValidatedGemmTuneRow::admit_receipt_json(
                GEMM_TUNE_ROW_RECEIPT_DOMAIN,
                GEMM_TUNE_ROW_RECEIPT_IDENTITY_VERSION,
                &base_json,
            ),
            Some(base_identity),
            "the current writer must be an exact parser fixed point"
        );

        let mut variants = Vec::new();
        let mut kernel = base_row.clone();
        kernel.kernel.push_str("-v2");
        variants.push(("kernel", kernel, false));
        let mut shape = base_row.clone();
        shape.shape_class.push_str("-v2");
        variants.push(("shape-class", shape, false));
        let mut machine = base_row.clone();
        machine.machine[0] ^= 1;
        variants.push(("machine-fingerprint", machine, false));
        let mut params = base_row.clone();
        params.params = "\"mc=32,nc-cap=512\"".to_string();
        variants.push(("selected-params", params, false));
        let mut measured = base_row.clone();
        measured.measured = identity_test_measured_tune_row(
            &base_row.kernel,
            &base_row.shape_class,
            u64::from_le_bytes(base_row.machine),
            2,
        );
        variants.push(("measured-row", measured, true));
        let mut memory_limit = base_row.clone();
        memory_limit.memory_limit_bytes += 1;
        variants.push(("memory-limit-bytes", memory_limit, true));
        let mut probe_buffers = base_row.clone();
        probe_buffers.probe_buffer_bytes += 1;
        variants.push(("probe-buffer-bytes", probe_buffers, true));

        let mut identities = std::collections::BTreeSet::new();
        for (field, variant, should_admit) in variants {
            let json = variant.receipt_json();
            let identity = variant.receipt_identity();
            assert_ne!(identity, base_identity, "{field} did not move identity");
            assert!(
                identities.insert(identity.to_string()),
                "{field} collided with another independent field mutation"
            );
            let admitted = ValidatedGemmTuneRow::admit_receipt_json(
                GEMM_TUNE_ROW_RECEIPT_DOMAIN,
                GEMM_TUNE_ROW_RECEIPT_IDENTITY_VERSION,
                &json,
            );
            if should_admit {
                assert_eq!(
                    admitted,
                    Some(identity),
                    "{field} did not survive the exact current transport"
                );
            } else {
                assert!(
                    admitted.is_none(),
                    "outer-only {field} mutation disagrees with the embedded fs-exec row"
                );
            }
        }

        let stale_domain = "org.frankensim.fs-session.gemm-tune-row-receipt.v3";
        assert!(
            ValidatedGemmTuneRow::admit_receipt_json(
                stale_domain,
                GEMM_TUNE_ROW_RECEIPT_IDENTITY_VERSION,
                &base_json,
            )
            .is_none(),
            "artifact-domain mutation must be refused"
        );
        assert!(
            ValidatedGemmTuneRow::admit_receipt_json(
                GEMM_TUNE_ROW_RECEIPT_DOMAIN,
                GEMM_TUNE_ROW_RECEIPT_IDENTITY_VERSION + 1,
                &base_json,
            )
            .is_none(),
            "identity-version mutation must be refused"
        );

        let mut kernel_json = String::new();
        push_json_string(&mut kernel_json, &base_row.kernel);
        let mut shape_json = String::new();
        push_json_string(&mut shape_json, &base_row.shape_class);
        let kernel_field = format!("\"kernel\":{kernel_json}");
        let shape_field = format!("\"shape_class\":{shape_json}");
        let reordered = base_json.replacen(
            &format!("{kernel_field},{shape_field}"),
            &format!("{shape_field},{kernel_field}"),
            1,
        );
        let machine_hex = format!("{:016x}", u64::from_le_bytes(base_row.machine));
        let narrow_machine = base_json.replacen(
            &format!("\"machine\":\"{machine_hex}\""),
            &format!("\"machine\":\"{}\"", machine_hex.trim_start_matches('0')),
            1,
        );
        let moved_plan_schema = base_json.replacen(
            GEMM_TUNE_METADATA_PLAN_SCHEMA,
            "fs-session-tune-metadata-plan-v2",
            1,
        );
        let requested = gemm_tune_metadata_plan_bytes();
        let moved_plan_bytes = base_json.replacen(
            &format!("\"requested_bytes\":{requested}"),
            &format!("\"requested_bytes\":{}", requested + 1),
            1,
        );
        for (field, hostile) in [
            ("canonical-field-order", reordered),
            ("machine-hex-width", narrow_machine),
            ("metadata-plan-schema", moved_plan_schema),
            ("metadata-plan-requested-bytes", moved_plan_bytes),
        ] {
            assert_ne!(hostile, base_json, "{field} mutation did not change bytes");
            assert!(
                ValidatedGemmTuneRow::admit_receipt_json(
                    GEMM_TUNE_ROW_RECEIPT_DOMAIN,
                    GEMM_TUNE_ROW_RECEIPT_IDENTITY_VERSION,
                    &hostile,
                )
                .is_none(),
                "{field} mutation must fail exact replay admission"
            );
        }
    }

    #[test]
    fn gemm_tune_row_receipt_versions_fail_closed() {
        let row = identity_test_validated_tune_row();
        let json = row.receipt_json();
        for (domain, version) in [
            (
                "org.frankensim.fs-session.gemm-tune-row-receipt.v1",
                GEMM_TUNE_ROW_RECEIPT_IDENTITY_VERSION,
            ),
            (
                GEMM_TUNE_ROW_RECEIPT_DOMAIN,
                GEMM_TUNE_ROW_RECEIPT_IDENTITY_VERSION - 1,
            ),
            ("", GEMM_TUNE_ROW_RECEIPT_IDENTITY_VERSION),
        ] {
            assert!(
                ValidatedGemmTuneRow::admit_receipt_json(domain, version, &json).is_none(),
                "stale receipt metadata {domain:?} v{version} must fail closed"
            );
        }
        assert!(
            ValidatedGemmTuneRow::admit_receipt_json(
                GEMM_TUNE_ROW_RECEIPT_DOMAIN,
                GEMM_TUNE_ROW_RECEIPT_IDENTITY_VERSION,
                &format!("{json} "),
            )
            .is_none(),
            "trailing display content is not an exact retained transport"
        );

        let stale_nested_domain = row.receipt_json().replacen(
            "org.frankensim.fs-exec.tune-row.v2",
            "org.frankensim.fs-exec.tune-row.v1",
            1,
        );
        let stale_nested_version =
            row.receipt_json()
                .replacen("\"identity_version\":2", "\"identity_version\":1", 1);
        for stale in [stale_nested_domain, stale_nested_version] {
            assert!(
                ValidatedGemmTuneRow::admit_receipt_json(
                    GEMM_TUNE_ROW_RECEIPT_DOMAIN,
                    GEMM_TUNE_ROW_RECEIPT_IDENTITY_VERSION,
                    &stale,
                )
                .is_none(),
                "stale embedded fs-exec tune-row metadata must fail closed"
            );
        }
    }

    fn identity_test_execution_receipt() -> GemmExecutionReceipt {
        let operation = fs_exec::RunId(17);
        GemmExecutionReceipt {
            declared_run: operation.0,
            completed_tiles: 23,
            total_tiles: 23,
            memory: GemmMemoryReceipt {
                limit_bytes: 1 << 30,
                staging_bytes: 101,
                b_pack_bytes: 102,
                band_metadata_bytes: 103,
                pool_run_bytes: 104,
                arena_bytes_per_worker: 105,
                active_arena_workers: 7,
                arena_bytes: 106,
                requested_bytes: 107,
            },
            panels: vec![
                GemmPanelReceipt {
                    kernel: "fs-la/gemm-panel-a".to_string(),
                    mode: "deterministic".to_string(),
                    declared_run: fs_la::gemm_panel_run_id(operation, 0).0,
                    completed: 11,
                    total: 11,
                },
                GemmPanelReceipt {
                    kernel: "fs-la/gemm-panel-b".to_string(),
                    mode: "deterministic".to_string(),
                    declared_run: fs_la::gemm_panel_run_id(operation, 1).0,
                    completed: 12,
                    total: 12,
                },
            ],
        }
    }

    #[test]
    #[allow(clippy::too_many_lines)] // one table proves every independently framed source field
    fn gemm_execution_receipt_identity_fields_move_independently() {
        let base = identity_test_execution_receipt();
        let base_bytes = base.canonical_bytes().expect("canonical execution receipt");
        let base_identity = base.receipt_identity().expect("execution identity");
        assert_eq!(
            GemmExecutionReceipt::from_canonical_bytes(&base_bytes),
            Ok(base.clone()),
            "the writer and parser must be an exact field-preserving fixed point"
        );

        let mut variants = Vec::new();
        let mut declared_run = base.clone();
        declared_run.declared_run += 1;
        variants.push(("declared-run", declared_run));
        let mut completed_tiles = base.clone();
        completed_tiles.completed_tiles += 1;
        variants.push(("completed-tiles", completed_tiles));
        let mut total_tiles = base.clone();
        total_tiles.total_tiles += 1;
        variants.push(("total-tiles", total_tiles));
        let mut limit = base.clone();
        limit.memory.limit_bytes += 1;
        variants.push(("memory-limit-bytes", limit));
        let mut staging = base.clone();
        staging.memory.staging_bytes += 1;
        variants.push(("memory-staging-bytes", staging));
        let mut b_pack = base.clone();
        b_pack.memory.b_pack_bytes += 1;
        variants.push(("memory-b-pack-bytes", b_pack));
        let mut band_metadata = base.clone();
        band_metadata.memory.band_metadata_bytes += 1;
        variants.push(("memory-band-metadata-bytes", band_metadata));
        let mut pool_run = base.clone();
        pool_run.memory.pool_run_bytes += 1;
        variants.push(("memory-pool-run-bytes", pool_run));
        let mut arena_per_worker = base.clone();
        arena_per_worker.memory.arena_bytes_per_worker += 1;
        variants.push(("memory-arena-bytes-per-worker", arena_per_worker));
        let mut active_workers = base.clone();
        active_workers.memory.active_arena_workers += 1;
        variants.push(("memory-active-arena-workers", active_workers));
        let mut arena = base.clone();
        arena.memory.arena_bytes += 1;
        variants.push(("memory-arena-bytes", arena));
        let mut requested = base.clone();
        requested.memory.requested_bytes += 1;
        variants.push(("memory-requested-bytes", requested));
        let mut panel_count = base.clone();
        panel_count.panels.push(GemmPanelReceipt {
            kernel: "fs-la/gemm-panel-c".to_string(),
            mode: "deterministic".to_string(),
            declared_run: 99,
            completed: 1,
            total: 1,
        });
        variants.push(("panel-count", panel_count));
        let mut panel_order = base.clone();
        panel_order.panels.swap(0, 1);
        variants.push(("panel-order", panel_order));
        let mut panel_kernel = base.clone();
        panel_kernel.panels[0].kernel.push_str("-v2");
        variants.push(("panel-kernel", panel_kernel));
        let mut panel_mode = base.clone();
        panel_mode.panels[0].mode = "fast".to_string();
        variants.push(("panel-mode", panel_mode));
        let mut panel_run = base.clone();
        panel_run.panels[0].declared_run += 1;
        variants.push(("panel-declared-run", panel_run));
        let mut panel_completed = base.clone();
        panel_completed.panels[0].completed += 1;
        variants.push(("panel-completed", panel_completed));
        let mut panel_total = base.clone();
        panel_total.panels[0].total += 1;
        variants.push(("panel-total", panel_total));

        let mut identities = std::collections::BTreeSet::new();
        for (field, variant) in variants {
            let bytes = variant
                .canonical_bytes()
                .unwrap_or_else(|error| panic!("{field} canonical bytes: {error}"));
            let identity = variant
                .receipt_identity()
                .unwrap_or_else(|error| panic!("{field} identity: {error}"));
            assert_ne!(identity, base_identity, "{field} did not move identity");
            assert!(
                identities.insert(identity.to_string()),
                "{field} collided with another field mutation"
            );
            assert_eq!(
                GemmExecutionReceipt::from_canonical_bytes(&bytes),
                Ok(variant),
                "{field} did not preserve its exact typed transport"
            );
        }

        let stale_domain = base
            .canonical_bytes_with_schema(
                "org.frankensim.fs-session.gemm-execution-receipt.v2",
                GEMM_EXECUTION_RECEIPT_IDENTITY_VERSION,
            )
            .expect("stale-domain fixture");
        assert!(
            GemmExecutionReceipt::from_canonical_bytes(&stale_domain).is_err(),
            "artifact-domain mutation must be refused"
        );
        let stale_version = base
            .canonical_bytes_with_schema(
                GEMM_EXECUTION_RECEIPT_DOMAIN,
                GEMM_EXECUTION_RECEIPT_IDENTITY_VERSION + 1,
            )
            .expect("stale-version fixture");
        assert!(
            GemmExecutionReceipt::from_canonical_bytes(&stale_version).is_err(),
            "identity-version mutation must be refused"
        );

        let prefix_len = GEMM_EXECUTION_RECEIPT_MAGIC.len()
            + 1
            + 4
            + GEMM_EXECUTION_RECEIPT_DOMAIN.len()
            + 1
            + 4;
        let mut reordered = base_bytes;
        assert_eq!(reordered[prefix_len], EXEC_TAG_DECLARED_RUN);
        assert_eq!(reordered[prefix_len + 9], EXEC_TAG_COMPLETED_TILES);
        reordered[prefix_len..prefix_len + 18].rotate_left(9);
        assert!(
            GemmExecutionReceipt::from_canonical_bytes(&reordered).is_err(),
            "canonical-field-order mutation must be refused by exact tags"
        );
    }

    #[test]
    fn gemm_execution_receipt_versions_fail_closed() {
        let receipt = identity_test_execution_receipt();
        let stale_domain = receipt
            .canonical_bytes_with_schema(
                "org.frankensim.fs-session.gemm-execution-receipt.v0",
                GEMM_EXECUTION_RECEIPT_IDENTITY_VERSION,
            )
            .expect("stale domain fixture");
        let stale_version = receipt
            .canonical_bytes_with_schema(
                GEMM_EXECUTION_RECEIPT_DOMAIN,
                GEMM_EXECUTION_RECEIPT_IDENTITY_VERSION + 1,
            )
            .expect("stale version fixture");
        let mut trailing = receipt.canonical_bytes().expect("current fixture");
        trailing.push(0);
        for stale in [stale_domain, stale_version, trailing] {
            assert!(GemmExecutionReceipt::from_canonical_bytes(&stale).is_err());
        }
    }

    #[test]
    fn execution_receipt_excludes_schedule_measurements() {
        let operation_run = fs_exec::RunId(7);
        let no_product = fs_la::GemmRunReport {
            declared_run: operation_run,
            completed_tiles: 0,
            total_tiles: 0,
            pool_runs: Vec::new(),
            memory: fs_la::GemmMemoryReport::default(),
        };
        assert!(
            GemmExecutionReceipt::from_report(&no_product).is_complete(),
            "a successful no-product dispatch is complete without panel traversals"
        );
        let base_panel = fs_exec::RunReport {
            kernel: "fs-la/gemm-f64-m-band-v1",
            mode: "deterministic",
            declared_run: fs_la::gemm_panel_run_id(operation_run, 0),
            completed: 4,
            total: 4,
            steals: 0,
            cross_ccd_steals: 0,
            cancel_latencies_ns: Vec::new(),
            tiles_by_worker: vec![2, 2],
        };
        let first = fs_la::GemmRunReport {
            declared_run: operation_run,
            completed_tiles: 32,
            total_tiles: 32,
            pool_runs: vec![base_panel.clone()],
            memory: fs_la::GemmMemoryReport::default(),
        };
        let mut noisy_panel = base_panel;
        noisy_panel.steals = 99;
        noisy_panel.cross_ccd_steals = 17;
        noisy_panel.cancel_latencies_ns = vec![3, 5, 8];
        noisy_panel.tiles_by_worker = vec![4, 0];
        let mut second = fs_la::GemmRunReport {
            declared_run: operation_run,
            completed_tiles: 32,
            total_tiles: 32,
            pool_runs: vec![noisy_panel],
            memory: fs_la::GemmMemoryReport::default(),
        };
        second.memory.peak_used_bytes = 999;
        second.memory.refused_bytes = 17;
        assert_eq!(
            GemmExecutionReceipt::from_report(&first),
            GemmExecutionReceipt::from_report(&second),
            "steal, latency, and worker-distribution envelopes are not replay identity"
        );
        assert_eq!(
            GemmExecutionReceipt::from_report(&first)
                .receipt_identity()
                .expect("first identity"),
            GemmExecutionReceipt::from_report(&second)
                .receipt_identity()
                .expect("second identity"),
            "schedule-only and observed-memory mutations must not move the canonical identity"
        );
        assert!(GemmExecutionReceipt::from_report(&first).is_complete());
        let mut different_memory_plan = second.clone();
        different_memory_plan.memory.limit_bytes = 1 << 20;
        assert_ne!(
            GemmExecutionReceipt::from_report(&first),
            GemmExecutionReceipt::from_report(&different_memory_plan),
            "the declared memory plan is replay identity"
        );
        let mut different_run = first;
        different_run.pool_runs[0].declared_run = fs_exec::RunId(1);
        assert_ne!(
            GemmExecutionReceipt::from_report(&different_run),
            GemmExecutionReceipt::from_report(&second),
            "declared logical run is part of replay identity"
        );
    }

    fn synthetic_sweep() -> SweepResult {
        let winner = GemmBlockPlan::new(16, 512).expect("winner plan");
        let runner_up = GemmBlockPlan::new(32, 512).expect("runner-up plan");
        let evidence = TuneEvidence::ranked_wall_times(vec![
            TuneObservation::wall_time(winner.canonical(), vec![10, 11, 12])
                .expect("winner evidence"),
            TuneObservation::wall_time(runner_up.canonical(), vec![20, 21, 22])
                .expect("runner-up evidence"),
        ])
        .expect("ranked evidence");
        SweepResult { winner, evidence }
    }

    #[test]
    fn exact_bits_gate_catches_drift_in_every_repeat() {
        let candidates =
            effective_sweep_candidates(320, 2048, fs_la::GemmMemoryEnvelope::UNBOUNDED)
                .expect("candidate lattice");
        assert!(candidates.len() >= 2);
        for drift_repeat in 1..=SWEEP_SAMPLES {
            let mut call = 0usize;
            let error = measure_candidates(
                &CancelGate::new(),
                &candidates,
                2,
                fs_la::GemmMemoryEnvelope::UNBOUNDED,
                0,
                |_, c| {
                    let candidate = call / SWEEP_SAMPLES;
                    let repeat = call % SWEEP_SAMPLES + 1;
                    call += 1;
                    c[0] = 0.0;
                    c[1] = f64::from_bits(0x7ff8_0000_0000_0001);
                    if candidate == 1 && repeat == drift_repeat {
                        // Both changes are invisible to ordinary floating-point
                        // equality: signed zero compares equal and NaNs compare
                        // unequal regardless of payload. The contract is bits.
                        c[0] = -0.0;
                        c[1] = f64::from_bits(0x7ff8_0000_0000_0002);
                    }
                    Ok(0)
                },
            )
            .expect_err("the injected repeat must fail closed");
            assert!(
                matches!(
                    error,
                    GemmTuneError::BitDrift {
                        repeat,
                        ..
                    } if repeat == drift_repeat
                ),
                "repeat {drift_repeat}: {error}"
            );
        }
    }

    #[test]
    fn effective_candidate_lattice_is_unique_and_exercises_nc() {
        let narrow = effective_sweep_candidates(320, 288, fs_la::GemmMemoryEnvelope::UNBOUNDED)
            .expect("narrow lattice");
        assert_eq!(narrow.len(), SWEEP_MC.len());
        let narrow_pairs: std::collections::BTreeSet<_> = narrow
            .iter()
            .map(|candidate| (candidate.effective_mc, candidate.effective_nc))
            .collect();
        assert_eq!(narrow_pairs.len(), narrow.len());

        let wide = effective_sweep_candidates(320, 2048, fs_la::GemmMemoryEnvelope::UNBOUNDED)
            .expect("wide lattice");
        let wide_pairs: std::collections::BTreeSet<_> = wide
            .iter()
            .map(|candidate| (candidate.effective_mc, candidate.effective_nc))
            .collect();
        assert_eq!(wide_pairs.len(), wide.len());
        assert_eq!(
            wide.iter()
                .map(|candidate| candidate.effective_nc)
                .collect::<std::collections::BTreeSet<_>>(),
            std::collections::BTreeSet::from([512, 2048]),
            "n > 512 must measure both a multi-panel NC=512 execution and the wider panel"
        );

        let mut executed = Vec::new();
        measure_candidates(
            &CancelGate::new(),
            &wide,
            1,
            fs_la::GemmMemoryEnvelope::UNBOUNDED,
            0,
            |candidate, c| {
                executed.push((candidate.effective_mc, candidate.effective_nc));
                c[0] = 1.0;
                Ok(0)
            },
        )
        .expect("synthetic sweep");
        for pair in wide_pairs {
            assert_eq!(
                executed
                    .iter()
                    .filter(|&&observed| observed == pair)
                    .count(),
                SWEEP_SAMPLES,
                "each unique effective pair runs every repeat"
            );
        }
    }

    #[test]
    fn cancellation_between_repeats_returns_no_partial_evidence() {
        let candidates =
            effective_sweep_candidates(320, 2048, fs_la::GemmMemoryEnvelope::UNBOUNDED)
                .expect("candidate lattice");
        let gate = CancelGate::new();
        let error = measure_candidates(
            &gate,
            &candidates,
            1,
            fs_la::GemmMemoryEnvelope::UNBOUNDED,
            0,
            |_, c| {
                c[0] = 1.0;
                gate.request();
                Ok(0)
            },
        )
        .expect_err("the post-repeat poll must observe cancellation");
        assert!(matches!(
            error,
            GemmTuneError::Cancelled {
                peak_used_bytes: 16,
                report: None,
                ..
            }
        ));
    }

    #[test]
    fn cache_persistence_failure_is_atomic_and_retryable() {
        let key = gemm_tune_key(4, 320, 288, 300).expect("key");
        let mut tuner = Tuner::cold(0xAA55);
        let error = install_sweep_row(&mut tuner, &key, synthetic_sweep(), |_| {
            Err(GemmTuneError::Ledger("injected write failure".to_string()))
        })
        .expect_err("faulted cache write");
        assert!(matches!(error, GemmTuneError::Ledger(_)));
        assert!(!tuner.has_gemm_row(&key));
        assert!(tuner.decisions().is_empty());

        let (winner, _) = install_sweep_row(&mut tuner, &key, synthetic_sweep(), |_| Ok(()))
            .expect("retry installs the row");
        assert_eq!(winner, GemmBlockPlan::new(16, 512).expect("plan"));
        assert!(tuner.has_gemm_row(&key));
    }

    #[test]
    fn cached_params_and_body_must_agree_before_adoption() {
        let key = gemm_tune_key(4, 320, 288, 300).expect("key");
        let mut producer = Tuner::cold(0xAA55);
        let mut params = String::new();
        let mut measured = String::new();
        let mut sealed = None;
        install_sweep_row(&mut producer, &key, synthetic_sweep(), |validated| {
            params.clone_from(&validated.params);
            measured.clone_from(&validated.measured);
            sealed = Some(validated.clone());
            Ok(())
        })
        .expect("produce cached row");
        let sealed = sealed.expect("sealed row");

        let mut consumer = Tuner::cold(0xAA55);
        assert!(
            adopt_cached_row(&mut consumer, &key, "\"mc=32,nc-cap=512\"", &measured)
                .expect("mismatch is a cache miss")
                .is_none()
        );
        assert!(!consumer.has_gemm_row(&key));
        let adopted = adopt_cached_row(&mut consumer, &key, &params, &measured)
            .expect("adopt")
            .expect("validated adopted row");
        assert_eq!(adopted.receipt_identity(), sealed.receipt_identity());
        assert!(adopted.matches_decision(
            key.kernel(),
            key.shape_class(),
            0xAA55,
            "mc=16,nc-cap=512"
        ));
        assert!(consumer.has_gemm_row(&key));

        let original_identity = sealed.receipt_identity();
        let mut field_tampers = Vec::new();
        let mut tampered = sealed.clone();
        tampered.kernel.push('x');
        field_tampers.push(tampered);
        let mut tampered = sealed.clone();
        tampered.shape_class.push('x');
        field_tampers.push(tampered);
        let mut tampered = sealed.clone();
        tampered.machine[0] ^= 1;
        field_tampers.push(tampered);
        let mut tampered = sealed.clone();
        tampered.params.push(' ');
        field_tampers.push(tampered);
        let mut tampered = sealed.clone();
        tampered.measured.push(' ');
        field_tampers.push(tampered);
        let mut tampered = sealed.clone();
        tampered.memory_limit_bytes ^= 1;
        field_tampers.push(tampered);
        let mut tampered = sealed.clone();
        tampered.probe_buffer_bytes ^= 1;
        field_tampers.push(tampered);
        assert!(
            field_tampers
                .iter()
                .all(|tampered| tampered.receipt_identity() != original_identity),
            "every ledger tuple field must participate in the derive-key identity"
        );

        let other_probe = gemm_tune_key(4, 320, 289, 300).expect("other key");
        let mut wrong_context = Tuner::cold(0xAA55);
        assert!(
            adopt_cached_row(&mut wrong_context, &other_probe, &params, &measured)
                .expect("wrong context is a cache miss")
                .is_none()
        );
        assert!(!wrong_context.has_gemm_row(&other_probe));
    }

    #[test]
    fn cancelled_dispatch_preserves_progress_but_records_no_success_decision() {
        let key = gemm_tune_key(4, 320, 288, 300).expect("key");
        let mut tuner = Tuner::cold(0xAA55);
        tuner
            .pin_gemm_blocking(&key, GemmBlockPlan::COLD_START)
            .expect("pin");
        let decision = tuner.prepare_gemm_decision(&key);
        let error = execute_prepared_decision(&mut tuner, decision, |_| {
            Err::<(), _>(GemmTuneError::from(fs_la::GemmCancelled {
                report: Box::new(fs_la::GemmRunReport {
                    declared_run: fs_exec::RunId(9),
                    completed_tiles: 7,
                    total_tiles: 19,
                    pool_runs: Vec::new(),
                    memory: fs_la::GemmMemoryReport {
                        limit_bytes: 1_024,
                        requested_bytes: 512,
                        peak_used_bytes: 384,
                        ..fs_la::GemmMemoryReport::default()
                    },
                }),
            }))
        })
        .expect_err("cancelled producer");
        let GemmTuneError::Cancelled {
            limit_bytes,
            peak_used_bytes,
            report: Some(report),
        } = error
        else {
            panic!("expected retained cancellation report");
        };
        assert_eq!(limit_bytes, 1_024);
        assert_eq!(peak_used_bytes, 384);
        assert_eq!(report.completed_tiles, 7);
        assert_eq!(report.total_tiles, 19);
        assert_eq!(report.memory.requested_bytes, 512);
        assert_eq!(report.memory.peak_used_bytes, 384);
        assert!(tuner.decisions().is_empty());
        assert!(tuner.has_gemm_pin(&key));
    }

    #[test]
    fn executor_failure_retains_full_memory_and_progress_report() {
        let report = fs_la::GemmRunReport {
            declared_run: fs_exec::RunId(12),
            completed_tiles: 3,
            total_tiles: 11,
            pool_runs: Vec::new(),
            memory: fs_la::GemmMemoryReport {
                limit_bytes: 2_048,
                requested_bytes: 1_024,
                peak_used_bytes: 768,
                ..fs_la::GemmMemoryReport::default()
            },
        };
        let error = GemmTuneError::from(fs_la::GemmRunError::Executor {
            error: fs_exec::RunError::Incomplete {
                kernel: "fixture",
                tile: 4,
            },
            report: Box::new(report),
        });
        let GemmTuneError::Executor {
            limit_bytes,
            peak_used_bytes,
            report,
            ..
        } = error
        else {
            panic!("expected retained executor report");
        };
        assert_eq!(limit_bytes, 2_048);
        assert_eq!(peak_used_bytes, 768);
        assert_eq!(report.completed_tiles, 3);
        assert_eq!(report.total_tiles, 11);
        assert_eq!(report.memory.requested_bytes, 1_024);
        assert_eq!(report.memory.peak_used_bytes, 768);
    }

    #[test]
    fn error_reports_are_indirected_off_the_success_path() {
        assert!(
            core::mem::size_of::<GemmTuneError>() <= 128,
            "GemmTuneError grew to {} bytes; keep full reports behind indirection",
            core::mem::size_of::<GemmTuneError>()
        );
    }
}
