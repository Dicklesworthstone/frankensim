//! fs-ledger: the Design Ledger v0 (plan §11.2, Bet 10; Decalogue P4/P9).
//!
//! FrankenSQLite-backed system of record: content-addressed artifacts,
//! event-sourced ops with the frozen Five Explicits, lineage edges, metric
//! time series, the autotuner cache, and the fine-grained event stream —
//! one file per project, WAL mode, snapshot-isolated readers.
//!
//! Layer: L6 (HELM). Runtime dependencies: `std` + `fsqlite` (constellation).
//!
//! Concurrency contract (from FrankenSQLite's documented model): open one
//! [`Ledger`] (one connection) per thread within one process. Readers get
//! snapshot isolation and never block the appending writer; contention
//! surfaces as a retryable [`LedgerError::Busy`], never a hang or a silent
//! lost write. Multi-process multi-writer use is a no-claim boundary
//! (CONTRACT.md).
//!
//! Time travel, forkable worlds, `explain()`, the replay audit, and GC
//! live in the [`travel`] module (schema v2).

pub mod colors;
pub mod hash;
pub mod schema;
pub mod tombstone;
pub mod travel;
pub mod vcs;

pub use colors::{
    ColorDemotion, ColorGraph, ColorNode, ColorReplayError, ColorStructureRejection,
    ColorWriteError, NoSourceOriginVerifier, NoWaiverVerifier, PolicyDecision, SourceOrigin,
    SourceOriginRejection, SourceOriginRequest, SourceOriginVerifier, WAIVER_SCOPE_COLOR_UPGRADE,
    WAIVER_SCOPE_SOURCE_COLOR, Waiver, WaiverDependency, WaiverGrant, WaiverRejection,
    WaiverVerifier,
};
pub use hash::{Blake3, ContentHash, hash_bytes};
pub use schema::{ALL_TABLES, SCHEMA_VERSION, STORAGE_CHUNK_LEN, V1_TABLES};
pub use travel::{
    BranchDiff, BranchInfo, ExecMode, ExplainNode, ExplainOp, GcReport, MAIN_BRANCH,
    ReplayMismatch, ReplayVerdict, ViewSnapshot,
};

use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{SystemTime, UNIX_EPOCH};

use fsqlite::{Connection, FrankenError, SqliteValue};

/// Crate version (compile-time stamp).
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

pub(crate) const VCS_IDENTITY_EVENT_KIND: &str = "vcs-identity";

// ---------------------------------------------------------------------------
// Error model
// ---------------------------------------------------------------------------

/// Structured, machine-actionable ledger errors (Decalogue P10). Never a
/// panic across the crate boundary.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LedgerError {
    /// The database could not be opened or configured.
    Open {
        /// Path passed to [`Ledger::open`].
        path: String,
        /// Underlying engine detail.
        detail: String,
    },
    /// The file was written by a newer fs-ledger than this build supports.
    FutureSchema {
        /// `PRAGMA user_version` found in the file.
        found: i64,
        /// Highest version this build can read/write.
        supported: i64,
    },
    /// A SQL-level failure that is not retryable contention.
    Sql {
        /// What the ledger was doing.
        context: String,
        /// Underlying engine detail.
        detail: String,
    },
    /// Retryable contention (busy/locked/write-conflict). Retry the call,
    /// ideally with backoff; no state was silently lost.
    Busy {
        /// What the ledger was doing.
        context: String,
        /// Underlying engine detail.
        detail: String,
    },
    /// A required Five Explicits field is missing or malformed (P4/P10).
    MissingExplicit {
        /// Field name: `seed`, `versions`, `budget`, or `capability`.
        field: String,
        /// What is wrong and how to fix it.
        problem: String,
    },
    /// An input failed validation (structured; names the field).
    Invalid {
        /// Field name.
        field: String,
        /// What is wrong and how to fix it.
        problem: String,
    },
    /// Stored bytes no longer match their recorded content hash.
    Corrupt {
        /// Hex hash of the first corrupted artifact.
        hash_hex: String,
        /// Diagnosis.
        detail: String,
    },
    /// A referenced row does not exist.
    NotFound {
        /// Description of the missing row.
        what: String,
    },
    /// `finish_op` was called on an op that already has an outcome.
    DoubleFinish {
        /// The op id.
        op: i64,
    },
    /// An [`ArtifactWriter`] cannot be opened inside an explicit transaction
    /// (it owns its own transaction for crash atomicity).
    WriterInTransaction,
    /// The on-disk objects do not match the schema its `user_version`
    /// claims (bead gp3.18): pre-existing incompatible tables, a
    /// partially initialized file, wrong columns/affinities, missing
    /// indexes, or foreign objects. Refused BEFORE any migration
    /// advances the version — `CREATE TABLE IF NOT EXISTS` must never
    /// launder an alien schema into a labeled one.
    SchemaMismatch {
        /// The `PRAGMA user_version` the file claims.
        claimed_version: i64,
        /// Every attestation violation found (object-level diffs).
        violations: Vec<String>,
    },
    /// Identical bytes offered with a DIFFERENT envelope (kind or
    /// metadata) than the stored artifact carries (bead gp3.19).
    /// Byte deduplication stays content-addressed, but an envelope
    /// disagreement refuses instead of silently retaining whichever
    /// arrived first — provenance must not depend on insertion order.
    ArtifactEnvelopeConflict {
        /// Hex content hash of the artifact.
        hash_hex: String,
        /// Which envelope field conflicts: `kind` or `meta`.
        field: &'static str,
        /// The envelope value already stored (`<none>` for NULL meta).
        stored: String,
        /// The envelope value this call offered.
        offered: String,
    },
}

impl LedgerError {
    /// Stable machine-readable error code.
    #[must_use]
    pub fn code(&self) -> &'static str {
        match self {
            LedgerError::Open { .. } => "LedgerOpen",
            LedgerError::FutureSchema { .. } => "LedgerFutureSchema",
            LedgerError::Sql { .. } => "LedgerSql",
            LedgerError::Busy { .. } => "LedgerBusy",
            LedgerError::MissingExplicit { .. } => "LedgerMissingExplicit",
            LedgerError::Invalid { .. } => "LedgerInvalid",
            LedgerError::Corrupt { .. } => "LedgerCorruption",
            LedgerError::NotFound { .. } => "LedgerNotFound",
            LedgerError::DoubleFinish { .. } => "LedgerDoubleFinish",
            LedgerError::WriterInTransaction => "LedgerWriterInTransaction",
            LedgerError::SchemaMismatch { .. } => "LedgerSchemaMismatch",
            LedgerError::ArtifactEnvelopeConflict { .. } => "LedgerArtifactEnvelopeConflict",
        }
    }
}

impl std::fmt::Display for LedgerError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            LedgerError::Open { path, detail } => {
                write!(f, "LedgerOpen: cannot open ledger at {path}: {detail}")
            }
            LedgerError::FutureSchema { found, supported } => write!(
                f,
                "LedgerFutureSchema: file has schema v{found} but this build supports up to \
                 v{supported} — upgrade fs-ledger; refusing to touch a newer file"
            ),
            LedgerError::Sql { context, detail } => {
                write!(f, "LedgerSql: {context}: {detail}")
            }
            LedgerError::Busy { context, detail } => write!(
                f,
                "LedgerBusy: {context}: {detail} — retryable contention; retry with backoff"
            ),
            LedgerError::MissingExplicit { field, problem } => write!(
                f,
                "LedgerMissingExplicit: Five Explicits field '{field}' rejected: {problem} \
                 (units travel inside the typed IR; seeds/versions/budget/capability are \
                 mandatory columns — Decalogue P4)"
            ),
            LedgerError::Invalid { field, problem } => {
                write!(f, "LedgerInvalid: field '{field}' rejected: {problem}")
            }
            LedgerError::Corrupt { hash_hex, detail } => write!(
                f,
                "LedgerCorruption: artifact {hash_hex}: {detail} — refuse to trust or replay \
                 from a tampered ledger"
            ),
            LedgerError::NotFound { what } => write!(f, "LedgerNotFound: {what}"),
            LedgerError::DoubleFinish { op } => write!(
                f,
                "LedgerDoubleFinish: op {op} already has an outcome; ops are event-sourced \
                 facts and cannot be finished twice"
            ),
            LedgerError::WriterInTransaction => write!(
                f,
                "LedgerWriterInTransaction: ArtifactWriter owns its own transaction; commit \
                 or roll back the explicit transaction first"
            ),
            LedgerError::SchemaMismatch {
                claimed_version,
                violations,
            } => write!(
                f,
                "LedgerSchemaMismatch: the file claims schema v{claimed_version} but its \
                 objects do not attest ({} violation(s)): {} — refusing to advance \
                 user_version over an alien or partially initialized schema; migrate the \
                 data out manually or delete the file if it is disposable",
                violations.len(),
                violations.join("; ")
            ),
            LedgerError::ArtifactEnvelopeConflict {
                hash_hex,
                field,
                stored,
                offered,
            } => write!(
                f,
                "LedgerArtifactEnvelopeConflict: artifact {hash_hex} already stores \
                 {field}={stored} but this call offered {field}={offered}; identical \
                 bytes dedupe only under an AGREEING envelope — match the stored \
                 envelope (or offer no metadata to accept it) instead of relying on \
                 insertion order",
            ),
        }
    }
}

impl std::error::Error for LedgerError {}

fn is_retryable(e: &FrankenError) -> bool {
    matches!(
        e,
        FrankenError::Busy
            | FrankenError::BusyRecovery
            | FrankenError::BusySnapshot { .. }
            | FrankenError::DatabaseLocked { .. }
            | FrankenError::WriteConflict { .. }
            | FrankenError::SerializationFailure { .. }
    )
}

fn is_duplicate_key(e: &FrankenError) -> bool {
    matches!(
        e,
        FrankenError::PrimaryKeyViolation | FrankenError::UniqueViolation { .. }
    )
}

fn sql_err(context: &str, e: &FrankenError) -> LedgerError {
    if is_retryable(e) {
        LedgerError::Busy {
            context: context.to_string(),
            detail: e.to_string(),
        }
    } else {
        LedgerError::Sql {
            context: context.to_string(),
            detail: e.to_string(),
        }
    }
}

// ---------------------------------------------------------------------------
// Public row/receipt types
// ---------------------------------------------------------------------------

/// Receipt returned by artifact writes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PutReceipt {
    /// Content identity.
    pub hash: ContentHash,
    /// Total byte length.
    pub len: u64,
    /// True if identical bytes were already stored (no new row).
    pub deduped: bool,
    /// True if stored as chunk rows rather than inline bytes.
    pub chunked: bool,
}

/// Metadata for one stored artifact.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ArtifactInfo {
    /// Content identity.
    pub hash: ContentHash,
    /// Caller-declared kind (e.g. "field", "mesh", "study-ir").
    pub kind: String,
    /// Total byte length.
    pub len: u64,
    /// Number of chunk rows (0 = stored inline).
    pub chunk_count: u64,
    /// Optional JSON metadata.
    pub meta: Option<String>,
    /// Wall-clock nanoseconds at first insertion (provenance envelope; not
    /// part of the content identity).
    pub created_at: i64,
}

/// Outcome of a finished op.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OpOutcome {
    /// Completed successfully.
    Ok,
    /// Failed with a structured diagnostic.
    Error,
    /// Cancelled (request → drain → finalize; P7).
    Cancelled,
}

impl OpOutcome {
    fn as_str(self) -> &'static str {
        match self {
            OpOutcome::Ok => "ok",
            OpOutcome::Error => "error",
            OpOutcome::Cancelled => "cancelled",
        }
    }
}

/// Lineage edge direction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EdgeRole {
    /// The op consumed the artifact.
    In,
    /// The op produced the artifact.
    Out,
}

impl EdgeRole {
    fn as_str(self) -> &'static str {
        match self {
            EdgeRole::In => "in",
            EdgeRole::Out => "out",
        }
    }
}

/// The frozen Five Explicits of an op (P4). Units are the fifth explicit and
/// travel inside the typed IR itself (fs-qty dimensions), so they have no
/// separate column; the other four are mandatory here.
#[derive(Debug, Clone, Copy)]
pub struct FiveExplicits<'a> {
    /// RNG seed bytes (non-empty).
    pub seed: &'a [u8],
    /// Constellation/crate versions, JSON (e.g. the lock hash).
    pub versions: &'a str,
    /// Budget grant, JSON (accuracy/time/memory).
    pub budget: &'a str,
    /// Capability grant, JSON (ops/cores/mem/wall).
    pub capability: &'a str,
}

/// One recorded op row.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OpRow {
    /// Rowid.
    pub id: i64,
    /// Session identity, if any.
    pub session: Option<Vec<u8>>,
    /// Frozen IR (JSON).
    pub ir: String,
    /// Frozen seed bytes.
    pub seed: Vec<u8>,
    /// Frozen versions (JSON).
    pub versions: String,
    /// Frozen budget (JSON).
    pub budget: String,
    /// Frozen capability (JSON).
    pub capability: String,
    /// Start wall-clock ns.
    pub t_start: i64,
    /// End wall-clock ns (None while in flight).
    pub t_end: Option<i64>,
    /// Outcome text (None while in flight).
    pub outcome: Option<String>,
    /// Structured diagnostic (JSON), if any.
    pub diag: Option<String>,
}

/// One event-stream row to append.
#[derive(Debug, Clone, Copy)]
pub struct EventRow<'a> {
    /// Session identity, if any.
    pub session: Option<&'a [u8]>,
    /// Logical or wall time (caller-controlled so deterministic replays can
    /// use logical time).
    pub t: i64,
    /// Event kind (fs-obs kind names recommended).
    pub kind: &'a str,
    /// JSON payload, if any.
    pub payload: Option<&'a str>,
}

/// One autotuner cache row.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TuneRow {
    /// Kernel identity.
    pub kernel: String,
    /// Shape class.
    pub shape_class: String,
    /// Machine fingerprint bytes (fs-substrate probe hash).
    pub machine: Vec<u8>,
    /// Chosen parameters (JSON).
    pub params: String,
    /// Measured results (JSON).
    pub measured: String,
}

/// Rev S extension tables (sparse in v0; uniform `(name, body JSON)` shape).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExtensionTable {
    /// Requirement records.
    Requirements,
    /// Model cards for physical/surrogate models.
    ModelCards,
    /// Evidence records (Gauntlet artifacts, certificates).
    Evidence,
    /// Scenario definitions.
    Scenarios,
    /// Constraint definitions.
    Constraints,
    /// Hardware capability probe snapshots.
    CapabilityProbes,
    /// External import receipts.
    Imports,
    /// Unsafe-capsule registry mirror.
    UnsafeCapsules,
    /// Speculation telemetry: solve-node records carrying
    /// `(proposer_id, accepted, bound, iterations_saved)` (v3,
    /// bead lmp4.3).
    Speculation,
}

impl ExtensionTable {
    /// The underlying table name.
    #[must_use]
    pub fn table_name(self) -> &'static str {
        match self {
            ExtensionTable::Requirements => "requirements",
            ExtensionTable::Speculation => "speculation",
            ExtensionTable::ModelCards => "model_cards",
            ExtensionTable::Evidence => "evidence",
            ExtensionTable::Scenarios => "scenarios",
            ExtensionTable::Constraints => "constraints",
            ExtensionTable::CapabilityProbes => "capability_probes",
            ExtensionTable::Imports => "imports",
            ExtensionTable::UnsafeCapsules => "unsafe_capsules",
        }
    }
}

/// Referential/shape hygiene report. All-zero means clean.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct LintReport {
    /// Edges whose op id does not exist.
    pub orphan_edge_ops: u64,
    /// Edges whose artifact hash does not exist.
    pub orphan_edge_artifacts: u64,
    /// Metrics whose op id does not exist.
    pub orphan_metric_ops: u64,
    /// Artifacts violating the inline-XOR-chunked storage invariant.
    pub malformed_artifacts: u64,
    /// Chunked artifacts whose actual chunk-row count differs.
    pub chunk_count_mismatches: u64,
    /// Artifacts whose stored byte length differs from `len`.
    pub len_mismatches: u64,
    /// Chunk rows without a parent artifact row (e.g. abandoned staging).
    pub orphan_chunks: u64,
    /// Ops with exactly one of (t_end, outcome) set.
    pub half_finished_ops: u64,
    /// Ops whose branch id does not exist (v2).
    pub orphan_op_branches: u64,
    /// Branches whose parent id does not exist (v2).
    pub orphan_branch_parents: u64,
}

impl LintReport {
    /// True when every counter is zero.
    #[must_use]
    pub fn is_clean(&self) -> bool {
        *self == LintReport::default()
    }
}

/// Result of a full integrity re-hash.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct IntegrityReport {
    /// Artifacts checked.
    pub checked: u64,
    /// Hex hashes of artifacts whose bytes no longer match their identity.
    pub corrupted: Vec<String>,
}

impl IntegrityReport {
    /// True when nothing is corrupted.
    #[must_use]
    pub fn is_clean(&self) -> bool {
        self.corrupted.is_empty()
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Wall-clock nanoseconds since the Unix epoch (provenance envelope only;
/// never part of content identity — P2).
#[must_use]
pub fn now_wall_ns() -> i64 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0, |d| i64::try_from(d.as_nanos()).unwrap_or(i64::MAX))
}

/// Total i64 from a length (lengths are capped far below `i64::MAX` by the
/// engine's 1 GB value limit; saturate rather than wrap on absurd inputs).
fn int_from_usize(n: usize) -> i64 {
    i64::try_from(n).unwrap_or(i64::MAX)
}

fn int_from_u64(n: u64) -> i64 {
    i64::try_from(n).unwrap_or(i64::MAX)
}

/// Every user object in `sqlite_master` as `(type, name) -> normalized
/// SQL` (whitespace-collapsed; internal `sqlite_*` entries and DDL-less
/// auto-indexes excluded). The stored SQL text carries STRICT, CHECKs,
/// and FOREIGN KEY clauses verbatim, so equality attests them all.
fn schema_objects(
    conn: &Connection,
) -> Result<std::collections::BTreeMap<(String, String), String>, LedgerError> {
    let rows = conn
        .query("SELECT type, name, COALESCE(sql, '') FROM sqlite_master ORDER BY type, name")
        .map_err(|e| sql_err("attest: read sqlite_master", &e))?;
    let mut objects = std::collections::BTreeMap::new();
    for row in &rows {
        let kind = row_text(row, 0, "attest: sqlite_master type")?;
        let name = row_text(row, 1, "attest: sqlite_master name")?;
        // SQL LIKE treats `_` as a wildcard, so filtering with
        // `NOT LIKE 'sqlite_%'` would also hide legal user objects such as
        // `sqlitex_hidden`. Only the literal reserved prefix is internal.
        if name.starts_with("sqlite_") {
            continue;
        }
        let sql = row_text(row, 2, "attest: sqlite_master sql")?;
        let normalized = sql.split_whitespace().collect::<Vec<_>>().join(" ");
        objects.insert((kind, name), normalized);
    }
    Ok(objects)
}

/// One table's columns as stable signature strings
/// (`name:type:notnull:default:pk`) from `PRAGMA table_info`.
fn table_columns(conn: &Connection, table: &str) -> Result<Vec<String>, LedgerError> {
    let rows = conn
        .query(&format!("PRAGMA table_info({table})"))
        .map_err(|e| sql_err("attest: table_info", &e))?;
    let mut columns = Vec::with_capacity(rows.len());
    for row in &rows {
        let name = row_text(row, 1, "attest: column name")?;
        let declared = row_text(row, 2, "attest: column type")?;
        let not_null = row_i64(row, 3, "attest: column not-null")? != 0;
        let default_sql = match row.get(4) {
            Some(SqliteValue::Text(value)) => value.as_str().to_string(),
            Some(SqliteValue::Null) | None => "<none>".to_string(),
            Some(other) => format!("{other:?}"),
        };
        let pk = row_i64(row, 5, "attest: column pk")? != 0;
        columns.push(format!("{name}:{declared}:{not_null}:{default_sql}:{pk}"));
    }
    columns.sort();
    Ok(columns)
}

/// A fresh in-memory database with the first `steps` shipped migration
/// batches applied — the attestation reference.
fn reference_connection(steps: usize) -> Result<Connection, LedgerError> {
    let reference = Connection::open(":memory:").map_err(|e| LedgerError::Sql {
        context: "attest: open reference".to_string(),
        detail: e.to_string(),
    })?;
    for batch in schema::MIGRATIONS.iter().take(steps) {
        for ddl in *batch {
            reference.execute(ddl).map_err(|e| LedgerError::Sql {
                context: "attest: build reference".to_string(),
                detail: format!("{e} while executing: {}", ddl.get(..60).unwrap_or(ddl)),
            })?;
        }
    }
    Ok(reference)
}

fn row_text(row: &fsqlite::Row, idx: usize, context: &str) -> Result<String, LedgerError> {
    match row.get(idx) {
        Some(SqliteValue::Text(value)) => Ok(value.as_str().to_string()),
        other => Err(LedgerError::Sql {
            context: context.to_string(),
            detail: format!("expected TEXT at column {idx}, got {other:?}"),
        }),
    }
}

fn text_param(s: &str) -> SqliteValue {
    SqliteValue::Text(s.into())
}

fn blob_param(b: &[u8]) -> SqliteValue {
    SqliteValue::Blob(b.to_vec().into())
}

fn opt_text_param(s: Option<&str>) -> SqliteValue {
    s.map_or(SqliteValue::Null, text_param)
}

fn opt_blob_param(b: Option<&[u8]>) -> SqliteValue {
    b.map_or(SqliteValue::Null, blob_param)
}

fn row_i64(row: &fsqlite::Row, idx: usize, context: &str) -> Result<i64, LedgerError> {
    match row.get(idx) {
        Some(SqliteValue::Integer(v)) => Ok(*v),
        other => Err(LedgerError::Sql {
            context: context.to_string(),
            detail: format!("column {idx}: expected INTEGER, got {other:?}"),
        }),
    }
}

fn nonnegative_u64(value: i64, context: &str) -> Result<u64, LedgerError> {
    u64::try_from(value).map_err(|_| LedgerError::Sql {
        context: context.to_string(),
        detail: format!("expected non-negative INTEGER, got {value}"),
    })
}

fn row_u64(row: &fsqlite::Row, idx: usize, context: &str) -> Result<u64, LedgerError> {
    nonnegative_u64(row_i64(row, idx, context)?, context)
}

fn inline_artifact_bytes(row: &fsqlite::Row, expected_len: u64) -> Result<&[u8], String> {
    let Some(SqliteValue::Blob(bytes)) = row.get(0) else {
        return Err(format!("inline bytes: expected BLOB, got {:?}", row.get(0)));
    };
    let actual_len = u64::try_from(bytes.len())
        .map_err(|_| "inline byte length does not fit the ledger length domain".to_string())?;
    if actual_len != expected_len {
        return Err(format!(
            "inline length mismatch: recorded {expected_len}, found {actual_len}"
        ));
    }
    Ok(bytes.as_ref())
}

struct ArtifactChunkValidator {
    expected_len: u64,
    expected_count: u64,
    streamed: u64,
    actual_count: u64,
}

impl ArtifactChunkValidator {
    const fn new(info: &ArtifactInfo) -> Self {
        Self {
            expected_len: info.len,
            expected_count: info.chunk_count,
            streamed: 0,
            actual_count: 0,
        }
    }

    fn accept<'row>(&mut self, row: &'row fsqlite::Row) -> Result<&'row [u8], String> {
        let actual_seq = match row.get(0) {
            Some(SqliteValue::Integer(seq)) => u64::try_from(*seq)
                .map_err(|_| format!("chunk sequence must be non-negative, got {seq}"))?,
            other => {
                return Err(format!("chunk sequence: expected INTEGER, got {other:?}"));
            }
        };
        if actual_seq != self.actual_count {
            return Err(format!(
                "chunk sequence mismatch: expected {}, found {actual_seq}",
                self.actual_count
            ));
        }

        let Some(SqliteValue::Blob(bytes)) = row.get(1) else {
            return Err(format!("chunk bytes: expected BLOB, got {:?}", row.get(1)));
        };
        if bytes.len() > STORAGE_CHUNK_LEN {
            return Err(format!(
                "chunk {actual_seq} exceeds the {STORAGE_CHUNK_LEN}-byte storage bound: found {} \
                 bytes",
                bytes.len()
            ));
        }
        let chunk_len = u64::try_from(bytes.len()).map_err(|_| {
            format!("chunk {actual_seq} length does not fit the ledger length domain")
        })?;
        let next_count = self
            .actual_count
            .checked_add(1)
            .ok_or_else(|| "actual chunk count overflowed u64".to_string())?;
        if next_count > self.expected_count {
            return Err(format!(
                "chunk count exceeds recorded {} at sequence {actual_seq}",
                self.expected_count
            ));
        }
        let next_streamed = self
            .streamed
            .checked_add(chunk_len)
            .ok_or_else(|| format!("byte total overflowed u64 at chunk sequence {actual_seq}"))?;
        if next_streamed > self.expected_len {
            return Err(format!(
                "chunk bytes exceed recorded length {} at sequence {actual_seq}",
                self.expected_len
            ));
        }

        self.actual_count = next_count;
        self.streamed = next_streamed;
        Ok(bytes.as_ref())
    }

    fn finish(self) -> Result<u64, String> {
        if self.actual_count != self.expected_count {
            return Err(format!(
                "chunk count mismatch: recorded {}, found {}",
                self.expected_count, self.actual_count
            ));
        }
        if self.streamed != self.expected_len {
            return Err(format!(
                "chunk length mismatch: recorded {}, streamed {}",
                self.expected_len, self.streamed
            ));
        }
        Ok(self.streamed)
    }
}

// ---------------------------------------------------------------------------
// The Ledger
// ---------------------------------------------------------------------------

/// One handle on the Design Ledger: a single fsqlite connection plus the
/// schema/pragma contract. Not `Send`: open one per thread (the engine's
/// documented model); snapshot isolation happens below this API.
pub struct Ledger {
    conn: Connection,
    path: String,
}

impl Ledger {
    /// Open (creating and migrating as needed) the ledger at `path`.
    ///
    /// Applies the pragma contract: WAL journal, `synchronous=FULL`
    /// (fsync-before-publish durability), `busy_timeout`, and enforced
    /// foreign keys. `":memory:"` is supported for tests.
    ///
    /// # Errors
    /// [`LedgerError::Open`] on engine failure; [`LedgerError::FutureSchema`]
    /// if the file was written by a newer fs-ledger.
    pub fn open(path: &str) -> Result<Ledger, LedgerError> {
        let conn = Connection::open(path).map_err(|e| LedgerError::Open {
            path: path.to_string(),
            detail: e.to_string(),
        })?;
        for pragma in [
            "PRAGMA journal_mode=WAL",
            "PRAGMA synchronous=FULL",
            "PRAGMA busy_timeout=5000",
            "PRAGMA foreign_keys=ON",
        ] {
            conn.query(pragma).map_err(|e| LedgerError::Open {
                path: path.to_string(),
                detail: format!("{pragma}: {e}"),
            })?;
        }
        let ledger = Ledger {
            conn,
            path: path.to_string(),
        };
        ledger.migrate()?;
        Ok(ledger)
    }

    /// The path this ledger was opened at.
    #[must_use]
    pub fn path(&self) -> &str {
        &self.path
    }

    /// The schema version recorded in the file.
    ///
    /// # Errors
    /// [`LedgerError::Sql`] on engine failure.
    pub fn schema_version(&self) -> Result<i64, LedgerError> {
        let row = self
            .conn
            .query_row("PRAGMA user_version")
            .map_err(|e| sql_err("read user_version", &e))?;
        row_i64(&row, 0, "read user_version")
    }

    fn migrate(&self) -> Result<(), LedgerError> {
        let found = self.schema_version()?;
        if found > SCHEMA_VERSION {
            return Err(LedgerError::FutureSchema {
                found,
                supported: SCHEMA_VERSION,
            });
        }
        // ATTESTATION BEFORE ADVANCEMENT (bead gp3.18): `CREATE TABLE IF
        // NOT EXISTS` silently tolerates pre-existing objects, so without
        // this check an alien or partially initialized file would be
        // stamped as the current schema. A v0 file must be EMPTY; a file
        // claiming v>0 must attest object-for-object against a reference
        // built from the shipped DDL for that version.
        if found == 0 {
            // ATOMIC INITIALIZATION: the whole ladder plus the final
            // version stamp in ONE transaction. The emptiness attestation is
            // inside that same transaction so another connection cannot add
            // an object between the check and the first CREATE statement.
            self.conn
                .begin_transaction()
                .map_err(|e| sql_err("init: begin", &e))?;
            let init = (|| {
                let objects = schema_objects(&self.conn)?;
                if !objects.is_empty() {
                    let violations = objects
                        .keys()
                        .map(|(kind, name)| {
                            format!("pre-existing {kind} `{name}` in an unversioned file")
                        })
                        .collect();
                    return Err(LedgerError::SchemaMismatch {
                        claimed_version: 0,
                        violations,
                    });
                }
                for batch in schema::MIGRATIONS {
                    for ddl in *batch {
                        self.conn.execute(ddl).map_err(|error| LedgerError::Sql {
                            context: "initialize schema".to_string(),
                            detail: format!(
                                "{error} while executing: {}",
                                ddl.get(..60).unwrap_or(ddl)
                            ),
                        })?;
                    }
                }
                self.conn
                    .execute(&format!("PRAGMA user_version = {SCHEMA_VERSION}"))
                    .map_err(|error| sql_err("init: set user_version", &error))?;
                self.conn
                    .commit_transaction()
                    .map_err(|error| sql_err("init: commit", &error))
            })();
            if let Err(error) = init {
                let _ = self.conn.rollback_transaction();
                return Err(error);
            }
            return Ok(());
        }
        self.attest_schema(found)?;
        let start = usize::try_from(found).unwrap_or(usize::MAX);
        for (step, batch) in schema::MIGRATIONS.iter().enumerate().skip(start) {
            let target = step + 1;
            self.conn
                .begin_transaction()
                .map_err(|e| sql_err("migrate: begin", &e))?;
            let migration = (|| {
                for ddl in *batch {
                    let already_applied = schema::RECOVERABLE_ADDED_COLUMNS
                        .iter()
                        .find(|column| column.ddl == *ddl)
                        .map_or(Ok(false), |column| {
                            self.recoverable_column_is_present(target, column)
                        })?;
                    if already_applied {
                        continue;
                    }
                    self.conn.execute(ddl).map_err(|error| LedgerError::Sql {
                        context: format!("migrate to v{target}"),
                        detail: format!(
                            "{error} while executing: {}",
                            ddl.get(..60).unwrap_or(ddl)
                        ),
                    })?;
                }
                self.conn
                    .execute(&format!("PRAGMA user_version = {target}"))
                    .map_err(|error| sql_err("migrate: set user_version", &error))?;
                self.conn
                    .commit_transaction()
                    .map_err(|error| sql_err("migrate: commit", &error))
            })();
            if let Err(error) = migration {
                let _ = self.conn.rollback_transaction();
                return Err(error);
            }
        }
        Ok(())
    }

    fn recoverable_column_is_present(
        &self,
        target: usize,
        expected: &schema::RecoverableAddedColumn,
    ) -> Result<bool, LedgerError> {
        let rows = self
            .conn
            .query(&format!("PRAGMA table_info({})", expected.table))
            .map_err(|error| sql_err("migrate: inspect recoverable column", &error))?;
        for row in &rows {
            let name = match row.get(1) {
                Some(SqliteValue::Text(value)) => value.as_str(),
                other => {
                    return Err(LedgerError::Sql {
                        context: format!("migrate to v{target}"),
                        detail: format!(
                            "PRAGMA table_info({}) returned a non-TEXT column name: {other:?}",
                            expected.table
                        ),
                    });
                }
            };
            if name != expected.name {
                continue;
            }
            let declared_type = match row.get(2) {
                Some(SqliteValue::Text(value)) => value.as_str(),
                other => {
                    return Err(LedgerError::Sql {
                        context: format!("migrate to v{target}"),
                        detail: format!(
                            "recoverable column {}.{} has non-TEXT declared type {other:?}",
                            expected.table, expected.name
                        ),
                    });
                }
            };
            let not_null = row_i64(row, 3, "migrate: inspect column not-null")? != 0;
            let default_sql = match row.get(4) {
                Some(SqliteValue::Null) => None,
                Some(SqliteValue::Text(value)) => Some(value.as_str()),
                other => {
                    return Err(LedgerError::Sql {
                        context: format!("migrate to v{target}"),
                        detail: format!(
                            "recoverable column {}.{} has invalid default metadata {other:?}",
                            expected.table, expected.name
                        ),
                    });
                }
            };
            let primary_key = row_i64(row, 5, "migrate: inspect column primary-key")? != 0;
            if declared_type.eq_ignore_ascii_case(expected.declared_type)
                && not_null == expected.not_null
                && default_sql == expected.default_sql
                && primary_key == expected.primary_key
            {
                return Ok(true);
            }
            return Err(LedgerError::Sql {
                context: format!("migrate to v{target}"),
                detail: format!(
                    "existing column {}.{} does not match the recoverable migration definition: \
                     type={declared_type:?}, not_null={not_null}, default={default_sql:?}, \
                     primary_key={primary_key}",
                    expected.table, expected.name
                ),
            });
        }
        Ok(false)
    }

    /// ATTEST the on-disk objects against a REFERENCE database built by
    /// replaying the shipped DDL up to `claimed` (bead gp3.18). Compares
    /// sqlite_master object-for-object (tables, indexes, triggers, views —
    /// the stored SQL text covers foreign keys, CHECKs, and STRICT) and
    /// every table column-for-column via PRAGMA table_info (name, declared
    /// type/affinity, not-null, default, primary key). RECOVERY TOLERANCE
    /// (the tt_001b crash-window contract, generalizing
    /// RECOVERABLE_ADDED_COLUMNS): objects or columns beyond `claimed` are
    /// tolerated IFF they match the CURRENT shipped schema exactly — a
    /// committed-DDL-but-stale-marker file heals, while any divergent
    /// early object fails closed.
    fn attest_schema(&self, claimed: i64) -> Result<(), LedgerError> {
        let steps = usize::try_from(claimed).unwrap_or(usize::MAX);
        let at_claimed = reference_connection(steps)?;
        let at_current = reference_connection(schema::MIGRATIONS.len())?;
        let mut violations = Vec::new();
        let expected = schema_objects(&at_claimed)?;
        let current = schema_objects(&at_current)?;
        let actual = schema_objects(&self.conn)?;
        for (key, sql) in &expected {
            match actual.get(key) {
                None => violations.push(format!("missing {} `{}`", key.0, key.1)),
                Some(found) if found != sql && Some(found) != current.get(key) => violations.push(
                    format!("{} `{}` differs from the shipped definition", key.0, key.1),
                ),
                Some(_) => {}
            }
        }
        for (key, sql) in &actual {
            if !expected.contains_key(key) && current.get(key) != Some(sql) {
                violations.push(format!(
                    "unexpected {} `{}` (conflicting or foreign object)",
                    key.0, key.1
                ));
            }
        }
        // Column-level attestation for every expected table: COLUMN-level
        // diagnostics, with the same future-form tolerance (a column from a
        // committed-but-unmarked later batch must match its shipped
        // definition exactly).
        for (kind, table) in expected.keys() {
            if kind != "table" || !actual.contains_key(&(kind.clone(), table.clone())) {
                continue;
            }
            let want = table_columns(&at_claimed, table)?;
            let full = table_columns(&at_current, table)?;
            let have = table_columns(&self.conn, table)?;
            for col in &want {
                if !have.contains(col) {
                    violations.push(format!("table `{table}`: missing or altered column {col}"));
                }
            }
            for col in &have {
                if !want.contains(col) && !full.contains(col) {
                    violations.push(format!("table `{table}`: unexpected column {col}"));
                }
            }
        }
        if violations.is_empty() {
            Ok(())
        } else {
            violations.sort();
            Err(LedgerError::SchemaMismatch {
                claimed_version: claimed,
                violations,
            })
        }
    }

    /// `json_valid` check through the same engine that enforces the schema
    /// CHECKs, so pre-validation and enforcement can never disagree.
    fn json_valid(&self, s: &str) -> Result<bool, LedgerError> {
        let rows = self
            .conn
            .query_with_params("SELECT json_valid(?1)", &[text_param(s)])
            .map_err(|e| sql_err("json_valid", &e))?;
        let row = rows.first().ok_or_else(|| LedgerError::Sql {
            context: "json_valid".to_string(),
            detail: "no row returned".to_string(),
        })?;
        Ok(row_i64(row, 0, "json_valid")? == 1)
    }

    fn require_json(&self, field: &str, value: &str, explicit: bool) -> Result<(), LedgerError> {
        let problem = if value.trim().is_empty() {
            Some("empty string; supply a JSON value".to_string())
        } else if !self.json_valid(value)? {
            Some(format!(
                "not valid JSON: {:?}",
                value.get(..40).unwrap_or(value)
            ))
        } else {
            None
        };
        match problem {
            None => Ok(()),
            Some(problem) if explicit => Err(LedgerError::MissingExplicit {
                field: field.to_string(),
                problem,
            }),
            Some(problem) => Err(LedgerError::Invalid {
                field: field.to_string(),
                problem,
            }),
        }
    }

    // -- transactions -------------------------------------------------------

    /// Begin an explicit transaction (for atomic op+edges+metrics groups).
    ///
    /// # Errors
    /// [`LedgerError::Busy`] on contention; [`LedgerError::Sql`] otherwise.
    pub fn begin(&self) -> Result<(), LedgerError> {
        self.conn
            .begin_transaction()
            .map_err(|e| sql_err("begin", &e))
    }

    /// Commit the explicit transaction.
    ///
    /// # Errors
    /// [`LedgerError::Busy`] on commit-time conflict (retry the whole group).
    pub fn commit(&self) -> Result<(), LedgerError> {
        self.conn
            .commit_transaction()
            .map_err(|e| sql_err("commit", &e))
    }

    /// Roll back the explicit transaction.
    ///
    /// # Errors
    /// [`LedgerError::Sql`] on engine failure.
    pub fn rollback(&self) -> Result<(), LedgerError> {
        self.conn
            .rollback_transaction()
            .map_err(|e| sql_err("rollback", &e))
    }

    /// True while an explicit transaction is open.
    #[must_use]
    pub fn in_transaction(&self) -> bool {
        self.conn.in_transaction()
    }

    // -- artifacts ----------------------------------------------------------

    /// Store `bytes` content-addressed. Identical bytes dedupe to one row,
    /// across forks, forever; the receipt says whether this call created a
    /// new row. Artifacts larger than [`STORAGE_CHUNK_LEN`] are stored as
    /// chunk rows (fsqlite has no incremental-blob API; RAM per row stays
    /// bounded).
    ///
    /// # Errors
    /// [`LedgerError::Invalid`] for empty `kind`/malformed `meta`;
    /// [`LedgerError::Busy`]/[`LedgerError::Sql`] from the engine.
    pub fn put_artifact(
        &self,
        kind: &str,
        bytes: &[u8],
        meta: Option<&str>,
    ) -> Result<PutReceipt, LedgerError> {
        self.validate_artifact_inputs(kind, meta)?;
        let h = hash_bytes(bytes);
        if let Some(info) = self.artifact_info(&h)? {
            self.attest_artifact_envelope(&h, &info, kind, meta)?;
            return Ok(PutReceipt {
                hash: h,
                len: info.len,
                deduped: true,
                chunked: info.chunk_count > 0,
            });
        }
        let len = bytes.len() as u64;
        if bytes.len() <= STORAGE_CHUNK_LEN {
            return self.insert_inline_artifact(&h, kind, bytes, meta);
        }
        // Chunked path, atomically.
        let owns_txn = !self.conn.in_transaction();
        if owns_txn {
            self.begin()?;
        }
        let result = self.insert_chunked_artifact(&h, kind, bytes, meta);
        match (&result, owns_txn) {
            (Ok(()), true) => {
                if let Err(e) = self.commit() {
                    let _ = self.rollback();
                    return Err(e);
                }
            }
            (Err(_), true) => {
                let _ = self.rollback();
            }
            _ => {}
        }
        result?;
        Ok(PutReceipt {
            hash: h,
            len,
            deduped: false,
            chunked: true,
        })
    }

    fn validate_artifact_inputs(&self, kind: &str, meta: Option<&str>) -> Result<(), LedgerError> {
        if kind.is_empty() {
            return Err(LedgerError::Invalid {
                field: "kind".to_string(),
                problem: "empty; name the artifact kind (e.g. \"field\", \"mesh\")".to_string(),
            });
        }
        if let Some(m) = meta {
            self.require_json("meta", m, false)?;
        }
        Ok(())
    }

    /// Envelope agreement gate at every dedupe site (bead gp3.19): a
    /// byte-identical artifact dedupes only when the offered `kind`
    /// matches exactly and the offered metadata (when a claim is made)
    /// canonically equals the stored metadata. Offering `meta: None`
    /// makes NO claim and accepts the stored envelope (the streaming
    /// dedupe contract); offering metadata against a row stored without
    /// any is a conflict — silent claim-dropping is the bug this gate
    /// removes.
    fn attest_artifact_envelope(
        &self,
        h: &ContentHash,
        info: &ArtifactInfo,
        kind: &str,
        meta: Option<&str>,
    ) -> Result<(), LedgerError> {
        if kind != info.kind {
            return Err(LedgerError::ArtifactEnvelopeConflict {
                hash_hex: h.to_hex(),
                field: "kind",
                stored: info.kind.clone(),
                offered: kind.to_string(),
            });
        }
        let Some(offered) = meta else {
            return Ok(());
        };
        let conflict = |stored: String| LedgerError::ArtifactEnvelopeConflict {
            hash_hex: h.to_hex(),
            field: "meta",
            stored,
            offered: offered.to_string(),
        };
        let Some(stored) = info.meta.as_deref() else {
            return Err(conflict("<none>".to_string()));
        };
        // Canonical comparison through the engine (whitespace-insensitive;
        // key order remains significant, as documented in CONTRACT.md).
        let rows = self
            .conn
            .query_with_params(
                "SELECT json(?1) = json(?2)",
                &[text_param(stored), text_param(offered)],
            )
            .map_err(|e| sql_err("artifact meta compare", &e))?;
        let equal = rows
            .first()
            .map_or(Ok(0), |row| row_i64(row, 0, "artifact meta compare"))?;
        if equal == 1 {
            Ok(())
        } else {
            Err(conflict(stored.to_string()))
        }
    }

    fn insert_inline_artifact(
        &self,
        h: &ContentHash,
        kind: &str,
        bytes: &[u8],
        meta: Option<&str>,
    ) -> Result<PutReceipt, LedgerError> {
        let insert = self
            .conn
            .prepare(
                "INSERT INTO artifacts(hash, kind, bytes, len, chunk_count, meta, created_at) \
                 VALUES (?1, ?2, ?3, ?4, 0, ?5, ?6)",
            )
            .map_err(|e| sql_err("artifact insert prepare", &e))?
            .execute_with_params(&[
                blob_param(h.as_bytes()),
                text_param(kind),
                blob_param(bytes),
                SqliteValue::Integer(int_from_usize(bytes.len())),
                opt_text_param(meta),
                SqliteValue::Integer(now_wall_ns()),
            ]);
        match insert {
            Ok(_) => Ok(PutReceipt {
                hash: *h,
                len: bytes.len() as u64,
                deduped: false,
                chunked: false,
            }),
            // A concurrent writer stored the same content first: that IS the
            // dedupe contract, not an error — but only under an AGREEING
            // envelope (gp3.19), so re-read theirs and attest.
            Err(e) if is_duplicate_key(&e) => {
                if let Some(info) = self.artifact_info(h)? {
                    self.attest_artifact_envelope(h, &info, kind, meta)?;
                }
                Ok(PutReceipt {
                    hash: *h,
                    len: bytes.len() as u64,
                    deduped: true,
                    chunked: false,
                })
            }
            Err(e) => Err(sql_err("artifact insert", &e)),
        }
    }

    fn insert_chunked_artifact(
        &self,
        h: &ContentHash,
        kind: &str,
        bytes: &[u8],
        meta: Option<&str>,
    ) -> Result<(), LedgerError> {
        let chunks = bytes.chunks(STORAGE_CHUNK_LEN);
        let chunk_count = int_from_usize(chunks.len());
        for (seq, chunk) in chunks.enumerate() {
            let insert = self
                .conn
                .prepare("INSERT INTO artifact_chunks(hash, seq, bytes) VALUES (?1, ?2, ?3)")
                .map_err(|e| sql_err("chunk insert prepare", &e))?
                .execute_with_params(&[
                    blob_param(h.as_bytes()),
                    SqliteValue::Integer(int_from_usize(seq)),
                    blob_param(chunk),
                ]);
            match insert {
                Ok(_) => {}
                Err(e) if is_duplicate_key(&e) => {
                    // Concurrent identical store; the other writer wins —
                    // if their committed envelope agrees (gp3.19).
                    if let Some(info) = self.artifact_info(h)? {
                        self.attest_artifact_envelope(h, &info, kind, meta)?;
                    }
                    return Ok(());
                }
                Err(e) => return Err(sql_err("chunk insert", &e)),
            }
        }
        self.conn
            .prepare(
                "INSERT INTO artifacts(hash, kind, bytes, len, chunk_count, meta, created_at) \
                 VALUES (?1, ?2, NULL, ?3, ?4, ?5, ?6)",
            )
            .map_err(|e| sql_err("artifact insert prepare", &e))?
            .execute_with_params(&[
                blob_param(h.as_bytes()),
                text_param(kind),
                SqliteValue::Integer(int_from_usize(bytes.len())),
                SqliteValue::Integer(chunk_count),
                opt_text_param(meta),
                SqliteValue::Integer(now_wall_ns()),
            ])
            .map_err(|e| sql_err("artifact insert", &e))?;
        Ok(())
    }

    /// Begin a streaming artifact write (for fields too large to hold in
    /// memory). The writer owns a transaction: a crash or drop before
    /// `finish` leaves zero residue.
    ///
    /// # Errors
    /// [`LedgerError::WriterInTransaction`] if an explicit transaction is
    /// open; engine errors otherwise.
    pub fn artifact_writer(&self, kind: &str) -> Result<ArtifactWriter<'_>, LedgerError> {
        self.validate_artifact_inputs(kind, None)?;
        if self.conn.in_transaction() {
            return Err(LedgerError::WriterInTransaction);
        }
        self.begin()?;
        Ok(ArtifactWriter {
            ledger: self,
            kind: kind.to_string(),
            hasher: Blake3::new(),
            provisional: provisional_key(),
            next_seq: 0,
            buf: Vec::new(),
            len: 0,
            finished: false,
        })
    }

    /// Metadata for one artifact, if present.
    ///
    /// # Errors
    /// Engine errors only; absence is `Ok(None)`.
    pub fn artifact_info(&self, h: &ContentHash) -> Result<Option<ArtifactInfo>, LedgerError> {
        let rows = self
            .conn
            .query_with_params(
                "SELECT kind, len, chunk_count, meta, created_at FROM artifacts WHERE hash = ?1",
                &[blob_param(h.as_bytes())],
            )
            .map_err(|e| sql_err("artifact_info", &e))?;
        let Some(row) = rows.first() else {
            return Ok(None);
        };
        let kind = match row.get(0) {
            Some(SqliteValue::Text(t)) => t.as_str().to_string(),
            other => {
                return Err(LedgerError::Sql {
                    context: "artifact_info".to_string(),
                    detail: format!("kind: expected TEXT, got {other:?}"),
                });
            }
        };
        let meta = match row.get(3) {
            Some(SqliteValue::Text(t)) => Some(t.as_str().to_string()),
            Some(SqliteValue::Null) | None => None,
            other => {
                return Err(LedgerError::Sql {
                    context: "artifact_info".to_string(),
                    detail: format!("meta: expected TEXT or NULL, got {other:?}"),
                });
            }
        };
        Ok(Some(ArtifactInfo {
            hash: *h,
            kind,
            len: row_u64(row, 1, "artifact_info.len")?,
            chunk_count: row_u64(row, 2, "artifact_info.chunk_count")?,
            meta,
            created_at: row_i64(row, 4, "artifact_info.created_at")?,
        }))
    }

    /// Fetch an artifact's full bytes (assembles chunked storage in memory —
    /// prefer [`Ledger::read_artifact_chunks`] for very large fields). The
    /// materializer reserves fallibly from bytes actually read; recorded
    /// length metadata is never used as an allocation request.
    ///
    /// # Errors
    /// Engine errors, [`LedgerError::Corrupt`] when storage does not match its
    /// recorded shape, or [`LedgerError::Invalid`] when memory cannot be
    /// reserved for the result; absence is `Ok(None)`.
    pub fn get_artifact(&self, h: &ContentHash) -> Result<Option<Vec<u8>>, LedgerError> {
        let mut out = Vec::new();
        let mut allocation_error = None;
        let streamed = self.read_artifact_chunks(h, &mut |chunk| {
            if allocation_error.is_some() {
                return;
            }
            if let Err(error) = out.try_reserve(chunk.len()) {
                allocation_error = Some(LedgerError::Invalid {
                    field: "artifact_bytes".to_string(),
                    problem: format!(
                        "could not reserve {} additional bytes ({error}); use \
                         Ledger::read_artifact_chunks to process the artifact incrementally",
                        chunk.len()
                    ),
                });
                return;
            }
            // `try_reserve` above established sufficient capacity, so this
            // copy cannot trigger another allocation.
            out.extend_from_slice(chunk);
        })?;
        if streamed.is_none() {
            return Ok(None);
        }
        if let Some(error) = allocation_error {
            return Err(error);
        }
        Ok(Some(out))
    }

    /// Stream an artifact's bytes chunk-by-chunk without materializing the
    /// whole value. Returns the total length streamed, or `Ok(None)` if the
    /// artifact does not exist.
    ///
    /// Each callback receives only a prefix whose row sequence, count, and
    /// cumulative length are valid at that point. A later row or the final
    /// exact-count/exact-length check can still return `Err`; callback side
    /// effects for an already-delivered prefix are not rolled back. The first
    /// row that violates a recorded bound is not delivered. The callback has
    /// no error channel of its own, so callers that perform fallible work must
    /// record that failure and make subsequent invocations no-ops.
    ///
    /// # Errors
    /// Engine errors; [`LedgerError::Corrupt`] on malformed rows, non-dense
    /// sequences, oversized chunks, arithmetic overflow, or disagreement
    /// with the recorded chunk count or byte length.
    pub fn read_artifact_chunks(
        &self,
        h: &ContentHash,
        f: &mut dyn FnMut(&[u8]),
    ) -> Result<Option<u64>, LedgerError> {
        let Some(info) = self.artifact_info(h)? else {
            return Ok(None);
        };
        if info.chunk_count == 0 {
            let rows = self
                .conn
                .query_with_params(
                    "SELECT bytes FROM artifacts WHERE hash = ?1",
                    &[blob_param(h.as_bytes())],
                )
                .map_err(|e| sql_err("read_artifact_chunks", &e))?;
            let row = rows.first().ok_or_else(|| LedgerError::Corrupt {
                hash_hex: h.to_hex(),
                detail: "artifact metadata disappeared before its inline bytes were read"
                    .to_string(),
            })?;
            let bytes =
                inline_artifact_bytes(row, info.len).map_err(|detail| LedgerError::Corrupt {
                    hash_hex: h.to_hex(),
                    detail,
                })?;
            f(bytes);
            return Ok(Some(info.len));
        }

        let mut validator = ArtifactChunkValidator::new(&info);
        let mut corrupt_detail = None;
        self.conn
            .query_with_params_for_each(
                "SELECT seq, bytes FROM artifact_chunks WHERE hash = ?1 ORDER BY seq",
                &[blob_param(h.as_bytes())],
                |row| {
                    if corrupt_detail.is_some() {
                        return Ok(());
                    }
                    match validator.accept(row) {
                        Ok(bytes) => f(bytes),
                        Err(detail) => corrupt_detail = Some(detail),
                    }
                    Ok(())
                },
            )
            .map_err(|e| sql_err("read_artifact_chunks", &e))?;
        if let Some(detail) = corrupt_detail {
            return Err(LedgerError::Corrupt {
                hash_hex: h.to_hex(),
                detail,
            });
        }
        validator
            .finish()
            .map(Some)
            .map_err(|detail| LedgerError::Corrupt {
                hash_hex: h.to_hex(),
                detail,
            })
    }

    /// Re-hash every stored artifact against its recorded identity.
    /// Byte-level corruption fails LOUDLY here (Decalogue P9).
    ///
    /// # Errors
    /// Engine errors; corruption is reported in the result, not an `Err`.
    pub fn verify_artifact_integrity(&self) -> Result<IntegrityReport, LedgerError> {
        let rows = self
            .conn
            .query("SELECT hash FROM artifacts")
            .map_err(|e| sql_err("integrity scan", &e))?;
        let mut report = IntegrityReport::default();
        for row in &rows {
            let stored = match row.get(0) {
                Some(SqliteValue::Blob(b)) => ContentHash::from_slice(b),
                _ => None,
            };
            let Some(stored) = stored else {
                report.corrupted.push("<malformed hash column>".to_string());
                continue;
            };
            let mut hasher = Blake3::new();
            report.checked += 1;
            match self.read_artifact_chunks(&stored, &mut |chunk| hasher.update(chunk)) {
                Ok(Some(_)) => {}
                Ok(None) | Err(LedgerError::Corrupt { .. }) => {
                    report.corrupted.push(stored.to_hex());
                    continue;
                }
                Err(error) => return Err(error),
            }
            if hasher.finalize() != stored {
                report.corrupted.push(stored.to_hex());
            }
        }
        Ok(report)
    }

    /// Test hook: overwrite one artifact's stored bytes so integrity checks
    /// can prove they fail loudly. Never call outside tests.
    ///
    /// # Errors
    /// Engine errors; [`LedgerError::NotFound`] if absent.
    pub fn corrupt_artifact_for_test(&self, h: &ContentHash) -> Result<(), LedgerError> {
        let Some(info) = self.artifact_info(h)? else {
            return Err(LedgerError::NotFound {
                what: format!("artifact {}", h.to_hex()),
            });
        };
        let (sql, params): (&str, Vec<SqliteValue>) = if info.chunk_count == 0 {
            (
                "UPDATE artifacts SET bytes = X'DEADBEEF', len = 4 WHERE hash = ?1",
                vec![blob_param(h.as_bytes())],
            )
        } else {
            (
                "UPDATE artifact_chunks SET bytes = X'DEADBEEF' WHERE hash = ?1 AND seq = 0",
                vec![blob_param(h.as_bytes())],
            )
        };
        self.conn
            .prepare(sql)
            .map_err(|e| sql_err("corrupt hook prepare", &e))?
            .execute_with_params(&params)
            .map_err(|e| sql_err("corrupt hook", &e))?;
        Ok(())
    }

    // -- ops and lineage ----------------------------------------------------

    /// Record the start of an op with its frozen Five Explicits (P4) on the
    /// main branch in deterministic mode (see [`Ledger::begin_op_on`] for
    /// forks and fast mode). The caller controls `t_start_ns` so
    /// deterministic replays can use logical time; [`now_wall_ns`] is the
    /// conventional wall-clock source.
    ///
    /// # Errors
    /// [`LedgerError::MissingExplicit`] naming the offending field;
    /// engine errors otherwise.
    pub fn begin_op(
        &self,
        session: Option<&[u8]>,
        ir: &str,
        explicits: &FiveExplicits<'_>,
        t_start_ns: i64,
    ) -> Result<i64, LedgerError> {
        self.begin_op_on(
            MAIN_BRANCH,
            ExecMode::Deterministic,
            session,
            ir,
            explicits,
            t_start_ns,
        )
    }

    /// Record an op's outcome. Each op finishes exactly once.
    ///
    /// # Errors
    /// [`LedgerError::DoubleFinish`] on a second finish;
    /// [`LedgerError::NotFound`] for an unknown op.
    pub fn finish_op(
        &self,
        op: i64,
        outcome: OpOutcome,
        diag: Option<&str>,
        t_end_ns: i64,
    ) -> Result<(), LedgerError> {
        if let Some(d) = diag {
            self.require_json("diag", d, false)?;
        }
        let affected = self
            .conn
            .prepare(
                "UPDATE ops SET t_end = ?1, outcome = ?2, diag = ?3 \
                 WHERE id = ?4 AND outcome IS NULL",
            )
            .map_err(|e| sql_err("op finish prepare", &e))?
            .execute_with_params(&[
                SqliteValue::Integer(t_end_ns),
                text_param(outcome.as_str()),
                opt_text_param(diag),
                SqliteValue::Integer(op),
            ])
            .map_err(|e| sql_err("op finish", &e))?;
        if affected == 1 {
            return Ok(());
        }
        match self.op(op)? {
            Some(_) => Err(LedgerError::DoubleFinish { op }),
            None => Err(LedgerError::NotFound {
                what: format!("op {op}"),
            }),
        }
    }

    /// Fetch one op row, if present.
    ///
    /// # Errors
    /// Engine errors; absence is `Ok(None)`.
    pub fn op(&self, id: i64) -> Result<Option<OpRow>, LedgerError> {
        let rows = self
            .conn
            .query_with_params(
                "SELECT id, session, ir, seed, versions, budget, capability, t_start, t_end, \
                 outcome, diag FROM ops WHERE id = ?1",
                &[SqliteValue::Integer(id)],
            )
            .map_err(|e| sql_err("op fetch", &e))?;
        let Some(row) = rows.first() else {
            return Ok(None);
        };
        let text_at = |idx: usize| -> Result<String, LedgerError> {
            match row.get(idx) {
                Some(SqliteValue::Text(t)) => Ok(t.as_str().to_string()),
                other => Err(LedgerError::Sql {
                    context: "op fetch".to_string(),
                    detail: format!("column {idx}: expected TEXT, got {other:?}"),
                }),
            }
        };
        let opt_text_at = |idx: usize| -> Option<String> {
            match row.get(idx) {
                Some(SqliteValue::Text(t)) => Some(t.as_str().to_string()),
                _ => None,
            }
        };
        let session = match row.get(1) {
            Some(SqliteValue::Blob(b)) => Some(b.to_vec()),
            _ => None,
        };
        let seed = match row.get(3) {
            Some(SqliteValue::Blob(b)) => b.to_vec(),
            other => {
                return Err(LedgerError::Sql {
                    context: "op fetch".to_string(),
                    detail: format!("seed: expected BLOB, got {other:?}"),
                });
            }
        };
        let t_end = match row.get(8) {
            Some(SqliteValue::Integer(v)) => Some(*v),
            _ => None,
        };
        Ok(Some(OpRow {
            id: row_i64(row, 0, "op.id")?,
            session,
            ir: text_at(2)?,
            seed,
            versions: text_at(4)?,
            budget: text_at(5)?,
            capability: text_at(6)?,
            t_start: row_i64(row, 7, "op.t_start")?,
            t_end,
            outcome: opt_text_at(9),
            diag: opt_text_at(10),
        }))
    }

    /// Link an op to an artifact in the lineage DAG. Foreign keys are
    /// enforced: both rows must exist.
    ///
    /// # Errors
    /// [`LedgerError::Invalid`] on a dangling reference (names which side).
    pub fn link(&self, op: i64, artifact: &ContentHash, role: EdgeRole) -> Result<(), LedgerError> {
        let insert = self
            .conn
            .prepare("INSERT INTO edges(op, artifact, role) VALUES (?1, ?2, ?3)")
            .map_err(|e| sql_err("edge insert prepare", &e))?
            .execute_with_params(&[
                SqliteValue::Integer(op),
                blob_param(artifact.as_bytes()),
                text_param(role.as_str()),
            ]);
        match insert {
            Ok(_) => Ok(()),
            Err(FrankenError::ForeignKeyViolation) => Err(LedgerError::Invalid {
                field: "edge".to_string(),
                problem: format!(
                    "op {op} or artifact {} does not exist; record both before linking",
                    artifact.to_hex()
                ),
            }),
            Err(e) => Err(sql_err("edge insert", &e)),
        }
    }

    /// Whether the lineage DAG contains this exact role-qualified edge.
    ///
    /// This is the verifier-side companion to [`Ledger::link`]: callers can
    /// prove that a content-addressed artifact was an input or output of the
    /// claimed operation without scanning or reconstructing the whole DAG.
    ///
    /// # Errors
    /// Engine errors only. Missing operations, artifacts, or edges return
    /// `Ok(false)`.
    pub fn edge_exists(
        &self,
        op: i64,
        artifact: &ContentHash,
        role: EdgeRole,
    ) -> Result<bool, LedgerError> {
        let rows = self
            .conn
            .query_with_params(
                "SELECT 1 FROM edges WHERE op = ?1 AND artifact = ?2 AND role = ?3 LIMIT 1",
                &[
                    SqliteValue::Integer(op),
                    blob_param(artifact.as_bytes()),
                    text_param(role.as_str()),
                ],
            )
            .map_err(|error| sql_err("edge existence query", &error))?;
        Ok(!rows.is_empty())
    }

    // -- metrics, events, tune ---------------------------------------------

    /// Append one metric sample. `value` must be finite (REAL NOT NULL).
    ///
    /// # Errors
    /// [`LedgerError::Invalid`] for non-finite values; engine errors
    /// otherwise (including PK conflicts on duplicate `(op, t, name)`).
    pub fn record_metric(
        &self,
        op: i64,
        t: i64,
        name: &str,
        value: f64,
    ) -> Result<(), LedgerError> {
        if !value.is_finite() {
            return Err(LedgerError::Invalid {
                field: "value".to_string(),
                problem: format!("{value} is not finite; metrics are REAL NOT NULL"),
            });
        }
        self.conn
            .prepare("INSERT INTO metrics(op, t, name, value) VALUES (?1, ?2, ?3, ?4)")
            .map_err(|e| sql_err("metric insert prepare", &e))?
            .execute_with_params(&[
                SqliteValue::Integer(op),
                SqliteValue::Integer(t),
                text_param(name),
                SqliteValue::Float(value),
            ])
            .map_err(|e| sql_err("metric insert", &e))?;
        Ok(())
    }

    /// Append one event-stream row; returns its rowid.
    ///
    /// # Errors
    /// [`LedgerError::Invalid`] for malformed payload JSON or a reserved
    /// internal event kind; engine errors.
    pub fn append_event(&self, event: &EventRow<'_>) -> Result<i64, LedgerError> {
        if event.kind == VCS_IDENTITY_EVENT_KIND {
            return Err(LedgerError::Invalid {
                field: "kind".to_string(),
                problem: format!(
                    "event kind {VCS_IDENTITY_EVENT_KIND:?} is reserved for ledger identity"
                ),
            });
        }
        self.append_event_unchecked(event)
    }

    pub(crate) fn append_vcs_identity_event(
        &self,
        event: &EventRow<'_>,
    ) -> Result<i64, LedgerError> {
        if event.kind != VCS_IDENTITY_EVENT_KIND {
            return Err(LedgerError::Invalid {
                field: "kind".to_string(),
                problem: "the VCS identity insertion path accepts only its reserved event kind"
                    .to_string(),
            });
        }
        self.append_event_unchecked(event)
    }

    fn append_event_unchecked(&self, event: &EventRow<'_>) -> Result<i64, LedgerError> {
        if let Some(p) = event.payload {
            self.require_json("payload", p, false)?;
        }
        self.conn
            .prepare("INSERT INTO events(session, t, kind, payload) VALUES (?1, ?2, ?3, ?4)")
            .map_err(|e| sql_err("event insert prepare", &e))?
            .execute_with_params(&[
                opt_blob_param(event.session),
                SqliteValue::Integer(event.t),
                text_param(event.kind),
                opt_text_param(event.payload),
            ])
            .map_err(|e| sql_err("event insert", &e))?;
        Ok(self.conn.last_insert_rowid())
    }

    /// Append a batch of events in one transaction (the append-heavy write
    /// path of plan §11.2).
    ///
    /// # Errors
    /// On any failure the whole batch rolls back (when this call owns the
    /// transaction).
    pub fn append_events(&self, batch: &[EventRow<'_>]) -> Result<(), LedgerError> {
        let owns_txn = !self.conn.in_transaction();
        if owns_txn {
            self.begin()?;
        }
        for event in batch {
            if let Err(e) = self.append_event(event) {
                if owns_txn {
                    let _ = self.rollback();
                }
                return Err(e);
            }
        }
        if owns_txn && let Err(e) = self.commit() {
            let _ = self.rollback();
            return Err(e);
        }
        Ok(())
    }

    /// Number of stored events (stress/verification helper).
    ///
    /// # Errors
    /// Engine errors.
    pub fn table_count(&self, table: &str) -> Result<u64, LedgerError> {
        if !ALL_TABLES.contains(&table) {
            return Err(LedgerError::Invalid {
                field: "table".to_string(),
                problem: format!("unknown table {table:?}; see fs_ledger::ALL_TABLES"),
            });
        }
        let row = self
            .conn
            .query_row(&format!("SELECT COUNT(*) FROM {table}"))
            .map_err(|e| sql_err("table_count", &e))?;
        row_u64(&row, 0, "table_count")
    }

    /// Upsert one autotuner cache row (`kernel` × `shape_class` × machine
    /// fingerprint).
    ///
    /// Implemented as UPDATE-then-INSERT rather than
    /// `INSERT .. ON CONFLICT .. DO UPDATE`: the upsert form corrupts the
    /// database at fsqlite HEAD when the conflict seek lands after a leaf
    /// split (bead u8og; Dicklesworthstone/frankensqlite#123, regression in
    /// upstream d1a543e). Revert to the single-statement upsert once #123 is
    /// fixed. Equivalent under this connection model: the engine is
    /// single-writer per connection, so no row can appear between the two
    /// statements.
    ///
    /// # Errors
    /// [`LedgerError::Invalid`] for malformed JSON; engine errors.
    pub fn tune_put(
        &self,
        kernel: &str,
        shape_class: &str,
        machine: &[u8],
        params: &str,
        measured: &str,
    ) -> Result<(), LedgerError> {
        self.require_json("params", params, false)?;
        self.require_json("measured", measured, false)?;
        let updated = self
            .conn
            .prepare(
                "UPDATE tune SET params = ?4, measured = ?5 \
                 WHERE kernel = ?1 AND shape_class = ?2 AND machine = ?3",
            )
            .map_err(|e| sql_err("tune update prepare", &e))?
            .execute_with_params(&[
                text_param(kernel),
                text_param(shape_class),
                blob_param(machine),
                text_param(params),
                text_param(measured),
            ])
            .map_err(|e| sql_err("tune update", &e))?;
        if updated == 0 {
            self.conn
                .prepare(
                    "INSERT INTO tune(kernel, shape_class, machine, params, measured) \
                     VALUES (?1, ?2, ?3, ?4, ?5)",
                )
                .map_err(|e| sql_err("tune insert prepare", &e))?
                .execute_with_params(&[
                    text_param(kernel),
                    text_param(shape_class),
                    blob_param(machine),
                    text_param(params),
                    text_param(measured),
                ])
                .map_err(|e| sql_err("tune insert", &e))?;
        }
        Ok(())
    }

    /// Insert one autotuner row only when its exact storage key is absent.
    /// Existing rows are never modified. Callers that require idempotent exact
    /// identity should fetch and compare after this call.
    ///
    /// # Errors
    /// [`LedgerError::Invalid`] for malformed JSON; engine errors otherwise.
    pub fn tune_put_if_absent(
        &self,
        kernel: &str,
        shape_class: &str,
        machine: &[u8],
        params: &str,
        measured: &str,
    ) -> Result<(), LedgerError> {
        self.require_json("params", params, false)?;
        self.require_json("measured", measured, false)?;
        self.conn
            .prepare(
                "INSERT INTO tune(kernel, shape_class, machine, params, measured) \
                 VALUES (?1, ?2, ?3, ?4, ?5) \
                 ON CONFLICT(kernel, shape_class, machine) DO NOTHING",
            )
            .map_err(|e| sql_err("tune insert-if-absent prepare", &e))?
            .execute_with_params(&[
                text_param(kernel),
                text_param(shape_class),
                blob_param(machine),
                text_param(params),
                text_param(measured),
            ])
            .map_err(|e| sql_err("tune insert-if-absent", &e))?;
        Ok(())
    }

    /// Fetch one autotuner cache row, if present.
    ///
    /// # Errors
    /// Engine errors; absence is `Ok(None)`.
    pub fn tune_get(
        &self,
        kernel: &str,
        shape_class: &str,
        machine: &[u8],
    ) -> Result<Option<TuneRow>, LedgerError> {
        let rows = self
            .conn
            .query_with_params(
                "SELECT params, measured FROM tune \
                 WHERE kernel = ?1 AND shape_class = ?2 AND machine = ?3",
                &[
                    text_param(kernel),
                    text_param(shape_class),
                    blob_param(machine),
                ],
            )
            .map_err(|e| sql_err("tune fetch", &e))?;
        let Some(row) = rows.first() else {
            return Ok(None);
        };
        let text_at = |idx: usize| -> Result<String, LedgerError> {
            match row.get(idx) {
                Some(SqliteValue::Text(t)) => Ok(t.as_str().to_string()),
                other => Err(LedgerError::Sql {
                    context: "tune fetch".to_string(),
                    detail: format!("column {idx}: expected TEXT, got {other:?}"),
                }),
            }
        };
        Ok(Some(TuneRow {
            kernel: kernel.to_string(),
            shape_class: shape_class.to_string(),
            machine: machine.to_vec(),
            params: text_at(0)?,
            measured: text_at(1)?,
        }))
    }

    /// All autotuner cache rows for one kernel, across shape classes and
    /// machine fingerprints (staleness scans: "a target that was never
    /// re-measured is a lie waiting to happen", plan §14.1).
    ///
    /// # Errors
    /// Engine errors; an unknown kernel is an empty vec.
    pub fn tune_rows(&self, kernel: &str) -> Result<Vec<TuneRow>, LedgerError> {
        let rows = self
            .conn
            .query_with_params(
                "SELECT shape_class, machine, params, measured FROM tune \
                 WHERE kernel = ?1 ORDER BY shape_class",
                &[text_param(kernel)],
            )
            .map_err(|e| sql_err("tune scan", &e))?;
        let mut out = Vec::with_capacity(rows.len());
        for row in &rows {
            let text_at = |idx: usize| -> Result<String, LedgerError> {
                match row.get(idx) {
                    Some(SqliteValue::Text(t)) => Ok(t.as_str().to_string()),
                    other => Err(LedgerError::Sql {
                        context: "tune scan".to_string(),
                        detail: format!("column {idx}: expected TEXT, got {other:?}"),
                    }),
                }
            };
            let machine = match row.get(1) {
                Some(SqliteValue::Blob(b)) => b.to_vec(),
                other => {
                    return Err(LedgerError::Sql {
                        context: "tune scan".to_string(),
                        detail: format!("machine: expected BLOB, got {other:?}"),
                    });
                }
            };
            out.push(TuneRow {
                kernel: kernel.to_string(),
                shape_class: text_at(0)?,
                machine,
                params: text_at(2)?,
                measured: text_at(3)?,
            });
        }
        Ok(out)
    }

    // -- Rev S extension tables ----------------------------------------------

    /// Upsert a named record in one of the Rev S extension tables.
    ///
    /// # Errors
    /// [`LedgerError::Invalid`] for empty name / malformed JSON; engine
    /// errors otherwise.
    pub fn put_extension(
        &self,
        table: ExtensionTable,
        name: &str,
        body_json: &str,
    ) -> Result<(), LedgerError> {
        if name.is_empty() {
            return Err(LedgerError::Invalid {
                field: "name".to_string(),
                problem: "empty; extension records are keyed by name".to_string(),
            });
        }
        self.require_json("body", body_json, false)?;
        self.conn
            .prepare(&format!(
                "INSERT INTO {}(name, body, created_at) VALUES (?1, ?2, ?3) \
                 ON CONFLICT(name) DO UPDATE SET body = excluded.body",
                table.table_name()
            ))
            .map_err(|e| sql_err("extension upsert prepare", &e))?
            .execute_with_params(&[
                text_param(name),
                text_param(body_json),
                SqliteValue::Integer(now_wall_ns()),
            ])
            .map_err(|e| sql_err("extension upsert", &e))?;
        Ok(())
    }

    /// Fetch a named extension record's JSON body, if present.
    ///
    /// # Errors
    /// Engine errors; absence is `Ok(None)`.
    pub fn get_extension(
        &self,
        table: ExtensionTable,
        name: &str,
    ) -> Result<Option<String>, LedgerError> {
        let rows = self
            .conn
            .query_with_params(
                &format!("SELECT body FROM {} WHERE name = ?1", table.table_name()),
                &[text_param(name)],
            )
            .map_err(|e| sql_err("extension fetch", &e))?;
        match rows.first().and_then(|r| r.get(0)) {
            Some(SqliteValue::Text(t)) => Ok(Some(t.as_str().to_string())),
            _ => Ok(None),
        }
    }

    // -- hygiene --------------------------------------------------------------

    /// Referential/shape hygiene scan. A crash-recovered ledger must lint
    /// clean: orphan edges are the acceptance criterion of the kill-storm
    /// battery.
    ///
    /// # Errors
    /// Engine errors.
    pub fn lint(&self) -> Result<LintReport, LedgerError> {
        let count = |sql: &str, context: &str| -> Result<u64, LedgerError> {
            let row = self.conn.query_row(sql).map_err(|e| sql_err(context, &e))?;
            row_u64(&row, 0, context)
        };
        Ok(LintReport {
            orphan_edge_ops: count(
                "SELECT COUNT(*) FROM edges e LEFT JOIN ops o ON e.op = o.id WHERE o.id IS NULL",
                "lint orphan_edge_ops",
            )?,
            orphan_edge_artifacts: count(
                "SELECT COUNT(*) FROM edges e LEFT JOIN artifacts a ON e.artifact = a.hash \
                 WHERE a.hash IS NULL",
                "lint orphan_edge_artifacts",
            )?,
            orphan_metric_ops: count(
                "SELECT COUNT(*) FROM metrics m LEFT JOIN ops o ON m.op = o.id \
                 WHERE o.id IS NULL",
                "lint orphan_metric_ops",
            )?,
            malformed_artifacts: count(
                "SELECT COUNT(*) FROM artifacts WHERE NOT \
                 ((bytes IS NOT NULL AND chunk_count = 0) OR \
                  (bytes IS NULL AND chunk_count > 0))",
                "lint malformed_artifacts",
            )?,
            chunk_count_mismatches: count(
                "SELECT COUNT(*) FROM artifacts a WHERE a.chunk_count > 0 AND \
                 (SELECT COUNT(*) FROM artifact_chunks c WHERE c.hash = a.hash) \
                 != a.chunk_count",
                "lint chunk_count_mismatches",
            )?,
            len_mismatches: count(
                "SELECT COUNT(*) FROM artifacts a WHERE \
                 (a.bytes IS NOT NULL AND length(a.bytes) != a.len) OR \
                 (a.chunk_count > 0 AND \
                  (SELECT COALESCE(SUM(length(c.bytes)), 0) FROM artifact_chunks c \
                   WHERE c.hash = a.hash) != a.len)",
                "lint len_mismatches",
            )?,
            orphan_chunks: count(
                "SELECT COUNT(*) FROM artifact_chunks c LEFT JOIN artifacts a \
                 ON c.hash = a.hash WHERE a.hash IS NULL",
                "lint orphan_chunks",
            )?,
            // AND/OR form rather than `(a IS NULL) != (b IS NULL)`: fsqlite
            // mis-associates postfix IS NULL against comparison operators
            // when re-parsing stored CHECK text (upstream bug; see bead).
            half_finished_ops: count(
                "SELECT COUNT(*) FROM ops WHERE \
                 (t_end IS NULL AND outcome IS NOT NULL) OR \
                 (t_end IS NOT NULL AND outcome IS NULL)",
                "lint half_finished_ops",
            )?,
            orphan_op_branches: count(
                "SELECT COUNT(*) FROM ops o LEFT JOIN branches b ON o.branch = b.id \
                 WHERE b.id IS NULL",
                "lint orphan_op_branches",
            )?,
            orphan_branch_parents: count(
                "SELECT COUNT(*) FROM branches c LEFT JOIN branches p ON c.parent = p.id \
                 WHERE c.parent IS NOT NULL AND p.id IS NULL",
                "lint orphan_branch_parents",
            )?,
        })
    }
}

// ---------------------------------------------------------------------------
// Streaming artifact writer
// ---------------------------------------------------------------------------

static WRITER_NONCE: AtomicU64 = AtomicU64::new(0);

/// A provisional (non-content) chunk key for staging streamed chunks inside
/// the writer's transaction. Collision with a real BLAKE3 content hash would
/// require finding a preimage; treated as impossible (CONTRACT.md).
fn provisional_key() -> [u8; 32] {
    let mut h = Blake3::new();
    h.update(b"fs-ledger provisional chunk key v1");
    h.update(&std::process::id().to_le_bytes());
    h.update(&WRITER_NONCE.fetch_add(1, Ordering::Relaxed).to_le_bytes());
    h.update(&now_wall_ns().to_le_bytes());
    h.finalize().0
}

/// Streaming content-addressed artifact writer (see
/// [`Ledger::artifact_writer`]). Bytes are hashed incrementally and staged
/// as chunk rows under a provisional key inside a writer-owned transaction;
/// `finish` resolves dedupe, rewrites the key to the final content hash, and
/// commits. Dropping without `finish` rolls everything back.
pub struct ArtifactWriter<'a> {
    ledger: &'a Ledger,
    kind: String,
    hasher: Blake3,
    provisional: [u8; 32],
    next_seq: i64,
    buf: Vec<u8>,
    len: u64,
    finished: bool,
}

impl ArtifactWriter<'_> {
    /// Absorb more bytes.
    ///
    /// # Errors
    /// Engine errors while flushing full chunks.
    pub fn write(&mut self, data: &[u8]) -> Result<(), LedgerError> {
        self.hasher.update(data);
        self.len += data.len() as u64;
        self.buf.extend_from_slice(data);
        while self.buf.len() > STORAGE_CHUNK_LEN {
            let rest = self.buf.split_off(STORAGE_CHUNK_LEN);
            let full = std::mem::replace(&mut self.buf, rest);
            self.flush_chunk(&full)?;
        }
        Ok(())
    }

    fn flush_chunk(&mut self, chunk: &[u8]) -> Result<(), LedgerError> {
        self.ledger
            .conn
            .prepare("INSERT INTO artifact_chunks(hash, seq, bytes) VALUES (?1, ?2, ?3)")
            .map_err(|e| sql_err("stream chunk prepare", &e))?
            .execute_with_params(&[
                blob_param(&self.provisional),
                SqliteValue::Integer(self.next_seq),
                blob_param(chunk),
            ])
            .map_err(|e| sql_err("stream chunk insert", &e))?;
        self.next_seq += 1;
        Ok(())
    }

    /// Finalize: dedupe, promote staged chunks to the content hash, commit.
    ///
    /// # Errors
    /// Engine errors; on error the transaction is rolled back and nothing
    /// is stored.
    pub fn finish(mut self, meta: Option<&str>) -> Result<PutReceipt, LedgerError> {
        let result = self.finish_inner(meta);
        self.finished = true;
        if result.is_err() {
            let _ = self.ledger.rollback();
        }
        result
    }

    fn finish_inner(&mut self, meta: Option<&str>) -> Result<PutReceipt, LedgerError> {
        if let Some(m) = meta {
            self.ledger.require_json("meta", m, false)?;
        }
        let h = self.hasher.finalize();
        if let Some(info) = self.ledger.artifact_info(&h)? {
            // Identical content already stored: keep theirs — if the
            // envelope agrees (gp3.19) — and discard staging.
            self.ledger
                .attest_artifact_envelope(&h, &info, &self.kind, meta)?;
            self.discard_staging()?;
            self.ledger.commit()?;
            return Ok(PutReceipt {
                hash: h,
                len: info.len,
                deduped: true,
                chunked: info.chunk_count > 0,
            });
        }
        if self.next_seq == 0 {
            // Everything fit in the buffer: store inline.
            let buf = std::mem::take(&mut self.buf);
            let receipt = self
                .ledger
                .insert_inline_artifact(&h, &self.kind, &buf, meta)?;
            self.ledger.commit()?;
            return Ok(receipt);
        }
        if !self.buf.is_empty() {
            let tail = std::mem::take(&mut self.buf);
            self.flush_chunk(&tail)?;
        }
        self.ledger
            .conn
            .prepare("UPDATE artifact_chunks SET hash = ?1 WHERE hash = ?2")
            .map_err(|e| sql_err("stream promote prepare", &e))?
            .execute_with_params(&[blob_param(h.as_bytes()), blob_param(&self.provisional)])
            .map_err(|e| sql_err("stream promote", &e))?;
        self.ledger
            .conn
            .prepare(
                "INSERT INTO artifacts(hash, kind, bytes, len, chunk_count, meta, created_at) \
                 VALUES (?1, ?2, NULL, ?3, ?4, ?5, ?6)",
            )
            .map_err(|e| sql_err("stream artifact prepare", &e))?
            .execute_with_params(&[
                blob_param(h.as_bytes()),
                text_param(&self.kind),
                SqliteValue::Integer(int_from_u64(self.len)),
                SqliteValue::Integer(self.next_seq),
                opt_text_param(meta),
                SqliteValue::Integer(now_wall_ns()),
            ])
            .map_err(|e| sql_err("stream artifact insert", &e))?;
        self.ledger.commit()?;
        Ok(PutReceipt {
            hash: h,
            len: self.len,
            deduped: false,
            chunked: true,
        })
    }

    fn discard_staging(&self) -> Result<(), LedgerError> {
        self.ledger
            .conn
            .prepare("DELETE FROM artifact_chunks WHERE hash = ?1")
            .map_err(|e| sql_err("stream discard prepare", &e))?
            .execute_with_params(&[blob_param(&self.provisional)])
            .map_err(|e| sql_err("stream discard", &e))?;
        Ok(())
    }
}

impl Drop for ArtifactWriter<'_> {
    fn drop(&mut self) {
        if !self.finished {
            let _ = self.ledger.rollback();
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn mem() -> Ledger {
        Ledger::open(":memory:").expect("open :memory:")
    }

    const FX: FiveExplicits<'static> = FiveExplicits {
        seed: &[0x5E, 0xED, 0x00, 0x01],
        versions: r#"{"constellation":"f92683cc4572a198"}"#,
        budget: r#"{"wall_s":10}"#,
        capability: r#"{"ops":["test.*"]}"#,
    };

    #[test]
    fn open_migrates_to_current_version() {
        let l = mem();
        assert_eq!(l.schema_version().unwrap(), SCHEMA_VERSION);
        for table in ALL_TABLES {
            // A fresh ledger is empty except the seeded main branch row.
            let expected = u64::from(*table == "branches");
            assert_eq!(
                l.table_count(table).unwrap(),
                expected,
                "{table} fresh count"
            );
        }
    }

    #[test]
    fn signed_database_counts_refuse_negative_values() {
        assert_eq!(nonnegative_u64(0, "fixture").unwrap(), 0);
        assert_eq!(
            nonnegative_u64(i64::MAX, "fixture").unwrap(),
            9_223_372_036_854_775_807u64
        );
        let error = nonnegative_u64(-1, "artifact_info.len").unwrap_err();
        assert!(matches!(
            error,
            LedgerError::Sql { context, detail }
                if context == "artifact_info.len"
                    && detail == "expected non-negative INTEGER, got -1"
        ));
    }

    fn insert_raw_chunked_artifact(
        ledger: &Ledger,
        tag: &[u8],
        len: i64,
        chunk_count: i64,
        chunks: &[(i64, &[u8])],
    ) -> ContentHash {
        let hash = hash_bytes(tag);
        for (seq, bytes) in chunks {
            ledger
                .conn
                .prepare("INSERT INTO artifact_chunks(hash, seq, bytes) VALUES (?1, ?2, ?3)")
                .unwrap()
                .execute_with_params(&[
                    blob_param(hash.as_bytes()),
                    SqliteValue::Integer(*seq),
                    blob_param(bytes),
                ])
                .unwrap();
        }
        ledger
            .conn
            .prepare(
                "INSERT INTO artifacts(hash, kind, bytes, len, chunk_count, meta, created_at) \
                 VALUES (?1, 'corrupt-fixture', NULL, ?2, ?3, NULL, 0)",
            )
            .unwrap()
            .execute_with_params(&[
                blob_param(hash.as_bytes()),
                SqliteValue::Integer(len),
                SqliteValue::Integer(chunk_count),
            ])
            .unwrap();
        hash
    }

    #[test]
    fn inline_artifact_length_mismatch_is_structured_corruption() {
        let ledger = mem();
        let receipt = ledger.put_artifact("blob", b"abc", None).unwrap();
        ledger
            .conn
            .prepare("UPDATE artifacts SET len = 4 WHERE hash = ?1")
            .unwrap()
            .execute_with_params(&[blob_param(receipt.hash.as_bytes())])
            .unwrap();

        let error = ledger.get_artifact(&receipt.hash).unwrap_err();
        assert!(matches!(
            error,
            LedgerError::Corrupt { hash_hex, detail }
                if hash_hex == receipt.hash.to_hex()
                    && detail == "inline length mismatch: recorded 4, found 3"
        ));
    }

    #[test]
    fn chunked_i64_max_length_is_rejected_without_metadata_allocation() {
        let ledger = mem();
        let hash =
            insert_raw_chunked_artifact(&ledger, b"i64-max-length", i64::MAX, 1, &[(0, b"x")]);

        // This metadata used to feed Vec::with_capacity and could panic or
        // attempt an enormous allocation before inspecting the one-byte row.
        let error = ledger.get_artifact(&hash).unwrap_err();
        assert!(matches!(
            error,
            LedgerError::Corrupt { hash_hex, detail }
                if hash_hex == hash.to_hex()
                    && detail
                        == "chunk length mismatch: recorded 9223372036854775807, streamed 1"
        ));
    }

    #[test]
    fn chunked_count_mismatch_returns_corruption_after_valid_prefix() {
        let ledger = mem();
        let hash = insert_raw_chunked_artifact(&ledger, b"count-mismatch", 1, 2, &[(0, b"x")]);
        let mut prefix = Vec::new();

        let error = ledger
            .read_artifact_chunks(&hash, &mut |chunk| prefix.extend_from_slice(chunk))
            .unwrap_err();
        assert_eq!(prefix, b"x", "a previously validated prefix is retained");
        assert!(matches!(
            error,
            LedgerError::Corrupt { hash_hex, detail }
                if hash_hex == hash.to_hex()
                    && detail == "chunk count mismatch: recorded 2, found 1"
        ));
    }

    #[test]
    fn chunked_sequence_gap_is_rejected_before_callback() {
        let ledger = mem();
        let hash = insert_raw_chunked_artifact(&ledger, b"sequence-gap", 1, 1, &[(1, b"x")]);
        let mut callback_count = 0;

        let error = ledger
            .read_artifact_chunks(&hash, &mut |_| callback_count += 1)
            .unwrap_err();
        assert_eq!(callback_count, 0, "the invalid row must not be delivered");
        assert!(matches!(
            error,
            LedgerError::Corrupt { hash_hex, detail }
                if hash_hex == hash.to_hex()
                    && detail == "chunk sequence mismatch: expected 0, found 1"
        ));
    }

    #[test]
    fn chunked_dense_rows_stream_and_materialize_exactly() {
        let ledger = mem();
        let hash =
            insert_raw_chunked_artifact(&ledger, b"dense-valid", 5, 2, &[(0, b"abc"), (1, b"de")]);
        let mut chunks = Vec::new();
        let streamed = ledger
            .read_artifact_chunks(&hash, &mut |chunk| chunks.push(chunk.to_vec()))
            .unwrap();

        assert_eq!(streamed, Some(5));
        assert_eq!(chunks, [b"abc".as_slice(), b"de".as_slice()]);
        assert_eq!(ledger.get_artifact(&hash).unwrap().unwrap(), b"abcde");
    }

    #[test]
    fn artifact_dedupe_inline() {
        let l = mem();
        let a = l.put_artifact("blob", b"same bytes", None).unwrap();
        let b = l.put_artifact("blob", b"same bytes", None).unwrap();
        assert!(!a.deduped);
        assert!(b.deduped);
        assert_eq!(a.hash, b.hash);
        assert_eq!(l.table_count("artifacts").unwrap(), 1);
        assert_eq!(l.get_artifact(&a.hash).unwrap().unwrap(), b"same bytes");
    }

    #[test]
    fn artifact_meta_must_be_json() {
        let l = mem();
        let err = l.put_artifact("blob", b"x", Some("not json")).unwrap_err();
        assert_eq!(err.code(), "LedgerInvalid");
    }

    #[test]
    fn five_explicits_are_enforced_field_by_field() {
        let l = mem();
        let empty_seed = FiveExplicits { seed: &[], ..FX };
        let err = l.begin_op(None, "{}", &empty_seed, 1).unwrap_err();
        assert!(matches!(err, LedgerError::MissingExplicit { ref field, .. } if field == "seed"));

        let bad_budget = FiveExplicits {
            budget: "not json",
            ..FX
        };
        let err = l.begin_op(None, "{}", &bad_budget, 1).unwrap_err();
        assert!(matches!(err, LedgerError::MissingExplicit { ref field, .. } if field == "budget"));

        let empty_versions = FiveExplicits { versions: "", ..FX };
        let err = l.begin_op(None, "{}", &empty_versions, 1).unwrap_err();
        assert!(
            matches!(err, LedgerError::MissingExplicit { ref field, .. } if field == "versions")
        );
    }

    #[test]
    fn op_lifecycle_and_double_finish() {
        let l = mem();
        let op = l
            .begin_op(Some(b"sess".as_slice()), r#"{"op":"test"}"#, &FX, 100)
            .unwrap();
        let row = l.op(op).unwrap().unwrap();
        assert_eq!(row.outcome, None);
        assert_eq!(row.seed, FX.seed);
        l.finish_op(op, OpOutcome::Ok, None, 200).unwrap();
        let row = l.op(op).unwrap().unwrap();
        assert_eq!(row.outcome.as_deref(), Some("ok"));
        assert_eq!(row.t_end, Some(200));
        let err = l.finish_op(op, OpOutcome::Error, None, 300).unwrap_err();
        assert_eq!(err, LedgerError::DoubleFinish { op });
        let err = l.finish_op(9999, OpOutcome::Ok, None, 1).unwrap_err();
        assert_eq!(err.code(), "LedgerNotFound");
    }

    #[test]
    fn edges_reject_dangling_references() {
        let l = mem();
        let ghost = hash_bytes(b"never stored");
        let op = l.begin_op(None, "{}", &FX, 1).unwrap();
        let err = l.link(op, &ghost, EdgeRole::Out).unwrap_err();
        assert_eq!(err.code(), "LedgerInvalid");
        let real = l.put_artifact("blob", b"real", None).unwrap();
        l.link(op, &real.hash, EdgeRole::Out).unwrap();
        assert!(l.edge_exists(op, &real.hash, EdgeRole::Out).unwrap());
        assert!(!l.edge_exists(op, &real.hash, EdgeRole::In).unwrap());
        assert!(!l.edge_exists(op + 1, &real.hash, EdgeRole::Out).unwrap());
        assert_eq!(l.table_count("edges").unwrap(), 1);
    }

    #[test]
    fn metrics_reject_non_finite() {
        let l = mem();
        let op = l.begin_op(None, "{}", &FX, 1).unwrap();
        assert!(l.record_metric(op, 0, "residual", 1.0e-9).is_ok());
        let err = l.record_metric(op, 1, "residual", f64::NAN).unwrap_err();
        assert_eq!(err.code(), "LedgerInvalid");
    }

    #[test]
    fn events_batch_is_atomic() {
        let l = mem();
        let good = EventRow {
            session: None,
            t: 1,
            kind: "tile_complete",
            payload: None,
        };
        let bad = EventRow {
            session: None,
            t: 2,
            kind: "tile_complete",
            payload: Some("nope"),
        };
        let err = l.append_events(&[good, bad]).unwrap_err();
        assert_eq!(err.code(), "LedgerInvalid");
        assert_eq!(l.table_count("events").unwrap(), 0, "batch rolled back");
        l.append_events(&[good, good]).unwrap();
        assert_eq!(l.table_count("events").unwrap(), 2);
    }

    #[test]
    fn tune_upserts() {
        let l = mem();
        l.tune_put(
            "gemm",
            "f64-512",
            b"m1",
            r#"{"mc":256}"#,
            r#"{"gflops":100}"#,
        )
        .unwrap();
        l.tune_put(
            "gemm",
            "f64-512",
            b"m1",
            r#"{"mc":384}"#,
            r#"{"gflops":120}"#,
        )
        .unwrap();
        assert_eq!(l.table_count("tune").unwrap(), 1);
        let row = l.tune_get("gemm", "f64-512", b"m1").unwrap().unwrap();
        assert_eq!(row.params, r#"{"mc":384}"#);
        assert!(l.tune_get("gemm", "f64-512", b"m2").unwrap().is_none());
    }

    #[test]
    fn tune_insert_if_absent_preserves_conflicts_and_transactions() {
        let l = mem();
        let original_params = r#"{"mc":256}"#;
        let original_measured = r#"{"gflops":100}"#;
        l.tune_put_if_absent("gemm", "f64-512", b"m1", original_params, original_measured)
            .expect("insert absent row");
        l.tune_put_if_absent("gemm", "f64-512", b"m1", original_params, original_measured)
            .expect("identical insert is an idempotent no-op");
        l.tune_put_if_absent(
            "gemm",
            "f64-512",
            b"m1",
            r#"{"mc":384}"#,
            r#"{"gflops":120}"#,
        )
        .expect("conflicting insert is a non-overwriting no-op");
        let retained = l
            .tune_get("gemm", "f64-512", b"m1")
            .expect("query")
            .expect("retained row");
        assert_eq!(retained.params, original_params);
        assert_eq!(retained.measured, original_measured);
        assert_eq!(l.table_count("tune").expect("count"), 1);

        let malformed = l
            .tune_put_if_absent("gemm", "bad-json", b"m1", "not-json", "{}")
            .expect_err("malformed params must be refused");
        assert_eq!(malformed.code(), "LedgerInvalid");
        let malformed = l
            .tune_put_if_absent("gemm", "bad-json", b"m1", "{}", "not-json")
            .expect_err("malformed measured evidence must be refused");
        assert_eq!(malformed.code(), "LedgerInvalid");
        assert!(l.tune_get("gemm", "bad-json", b"m1").unwrap().is_none());

        l.begin().expect("begin transaction");
        l.tune_put_if_absent("gemm", "f64-1024", b"m1", "{}", "{}")
            .expect("transactional insert");
        assert!(
            l.tune_get("gemm", "f64-1024", b"m1")
                .expect("query in transaction")
                .is_some()
        );
        l.rollback().expect("rollback transaction");
        assert!(
            l.tune_get("gemm", "f64-1024", b"m1")
                .expect("query after rollback")
                .is_none()
        );
    }

    #[test]
    fn tune_rows_scans_across_machines() {
        let l = mem();
        l.tune_put("axpy", "roofline-v1", b"m1", "{}", r#"{"gbs":50}"#)
            .unwrap();
        l.tune_put("axpy", "roofline-v1", b"m2", "{}", r#"{"gbs":60}"#)
            .unwrap();
        l.tune_put("gemm", "f64-512", b"m1", "{}", "{}").unwrap();
        let rows = l.tune_rows("axpy").unwrap();
        assert_eq!(rows.len(), 2);
        assert!(rows.iter().all(|r| r.kernel == "axpy"));
        assert!(rows.iter().any(|r| r.machine == b"m1"));
        assert!(rows.iter().any(|r| r.machine == b"m2"));
        assert!(l.tune_rows("nonexistent").unwrap().is_empty());
    }

    #[test]
    fn extension_tables_upsert_and_fetch() {
        let l = mem();
        l.put_extension(
            ExtensionTable::UnsafeCapsules,
            "fs-simd/neon",
            r#"{"lines":134}"#,
        )
        .unwrap();
        l.put_extension(
            ExtensionTable::UnsafeCapsules,
            "fs-simd/neon",
            r#"{"lines":140}"#,
        )
        .unwrap();
        assert_eq!(
            l.get_extension(ExtensionTable::UnsafeCapsules, "fs-simd/neon")
                .unwrap()
                .as_deref(),
            Some(r#"{"lines":140}"#)
        );
        assert!(
            l.get_extension(ExtensionTable::Evidence, "missing")
                .unwrap()
                .is_none()
        );
    }

    #[test]
    fn lint_clean_on_fresh_and_populated_ledger() {
        let l = mem();
        assert!(l.lint().unwrap().is_clean());
        let op = l.begin_op(None, "{}", &FX, 1).unwrap();
        let art = l.put_artifact("blob", b"bytes", None).unwrap();
        l.link(op, &art.hash, EdgeRole::Out).unwrap();
        l.record_metric(op, 0, "m", 1.0).unwrap();
        l.finish_op(op, OpOutcome::Ok, None, 2).unwrap();
        assert!(l.lint().unwrap().is_clean());
    }

    #[test]
    fn integrity_detects_corruption() {
        let l = mem();
        let a = l.put_artifact("blob", b"precious bytes", None).unwrap();
        assert!(l.verify_artifact_integrity().unwrap().is_clean());
        l.corrupt_artifact_for_test(&a.hash).unwrap();
        let report = l.verify_artifact_integrity().unwrap();
        assert_eq!(report.corrupted, vec![a.hash.to_hex()]);
    }

    #[test]
    fn writer_streams_and_dedupes() {
        let l = mem();
        let data: Vec<u8> = (0..100_000u32).map(|i| (i % 251) as u8).collect();
        let mut w = l.artifact_writer("field").unwrap();
        for piece in data.chunks(7919) {
            w.write(piece).unwrap();
        }
        let r1 = w.finish(None).unwrap();
        assert_eq!(r1.hash, hash_bytes(&data));
        assert!(!r1.deduped);
        assert_eq!(l.get_artifact(&r1.hash).unwrap().unwrap(), data);
        // Same content again → dedupe, still one artifact row.
        let mut w = l.artifact_writer("field").unwrap();
        w.write(&data).unwrap();
        let r2 = w.finish(None).unwrap();
        assert!(r2.deduped);
        assert_eq!(l.table_count("artifacts").unwrap(), 1);
        assert!(l.lint().unwrap().is_clean());
    }

    #[test]
    fn dropped_writer_leaves_zero_residue() {
        let l = mem();
        let mut w = l.artifact_writer("field").unwrap();
        w.write(&[7u8; 10_000]).unwrap();
        drop(w);
        assert_eq!(l.table_count("artifacts").unwrap(), 0);
        assert_eq!(l.table_count("artifact_chunks").unwrap(), 0);
        assert!(!l.in_transaction());
        assert!(l.lint().unwrap().is_clean());
    }

    #[test]
    fn writer_rejected_inside_transaction() {
        let l = mem();
        l.begin().unwrap();
        let err = l
            .artifact_writer("field")
            .err()
            .map(|e| e.code().to_string());
        assert_eq!(err.as_deref(), Some("LedgerWriterInTransaction"));
        l.rollback().unwrap();
    }
}
