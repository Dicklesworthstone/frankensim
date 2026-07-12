//! Design Ledger schema v0 (plan §11.2 + Appendix D, patch Rev S extensions).
//!
//! All tables are STRICT. One deliberate divergence from Appendix D as
//! written: SQLite STRICT tables only admit INT/INTEGER/REAL/TEXT/BLOB/ANY
//! column types, so every `JSON` column in the appendix is declared `TEXT`
//! with a `json_valid(...)` CHECK constraint — same semantics, actually
//! enforceable. A second divergence: `artifacts` gains `len`/`chunk_count`
//! and the `artifact_chunks` sibling table, because fsqlite has no
//! incremental-blob API and multi-GiB fields must be stored as bounded-size
//! chunk rows (CONTRACT.md documents the storage invariant).
//!
//! Migrations are versioned through `PRAGMA user_version`; each version marker
//! is committed in the same transaction as its DDL. The v2 recovery metadata
//! also recognizes the exact columns an older build could commit before
//! crashing ahead of its formerly separate version bump. Schema v4 adds one
//! immutable database-instance identity row. The row is seeded by Rust inside
//! the same migration transaction as the table and version marker because the
//! identity bytes must come from the in-tree domain-separated generator rather
//! than from engine-specific SQL randomness.

/// The schema version this crate writes and reads.
pub const SCHEMA_VERSION: i64 = 4;

/// Storage chunk length for large artifacts (bytes). Artifacts strictly
/// larger than this are stored as `artifact_chunks` rows of at most this
/// size; smaller ones live inline in `artifacts.bytes`.
pub const STORAGE_CHUNK_LEN: usize = 4 * 1024 * 1024;

/// Migration ladder: `MIGRATIONS[i]` migrates a database at `user_version`
/// `i` to `i + 1`. Append-only; never edit a shipped batch.
pub(crate) const MIGRATIONS: &[&[&str]] = &[V1, V2, V3, V4];

/// v1: the six core tables (Appendix D), chunk storage, and the Rev S
/// extension tables (sparse in v0 but present EARLY so downstream crates can
/// rely on them existing). Public so migration tests can construct genuine
/// v1 databases and prove the upgrade path.
pub const V1: &[&str] = &[
    // -- core six ---------------------------------------------------------
    "CREATE TABLE IF NOT EXISTS artifacts(
        hash BLOB PRIMARY KEY CHECK(length(hash) = 32),
        kind TEXT NOT NULL CHECK(length(kind) > 0),
        bytes BLOB,
        len INTEGER NOT NULL CHECK(len >= 0),
        chunk_count INTEGER NOT NULL DEFAULT 0 CHECK(chunk_count >= 0),
        meta TEXT CHECK(meta IS NULL OR json_valid(meta)),
        created_at INTEGER NOT NULL,
        CHECK((bytes IS NOT NULL AND chunk_count = 0) OR (bytes IS NULL AND chunk_count > 0))
    ) STRICT",
    "CREATE TABLE IF NOT EXISTS artifact_chunks(
        hash BLOB NOT NULL,
        seq INTEGER NOT NULL CHECK(seq >= 0),
        bytes BLOB NOT NULL,
        PRIMARY KEY(hash, seq)
    ) STRICT",
    "CREATE TABLE IF NOT EXISTS ops(
        id INTEGER PRIMARY KEY,
        session BLOB,
        ir TEXT NOT NULL CHECK(json_valid(ir)),
        seed BLOB NOT NULL CHECK(length(seed) > 0),
        versions TEXT NOT NULL CHECK(json_valid(versions)),
        budget TEXT NOT NULL CHECK(json_valid(budget)),
        capability TEXT NOT NULL CHECK(json_valid(capability)),
        t_start INTEGER NOT NULL,
        t_end INTEGER,
        outcome TEXT CHECK(outcome IN ('ok','error','cancelled')),
        diag TEXT CHECK(diag IS NULL OR json_valid(diag)),
        CHECK((t_end IS NULL AND outcome IS NULL) OR
              (t_end IS NOT NULL AND outcome IS NOT NULL))
    ) STRICT",
    "CREATE TABLE IF NOT EXISTS edges(
        op INTEGER NOT NULL REFERENCES ops(id),
        artifact BLOB NOT NULL REFERENCES artifacts(hash),
        role TEXT NOT NULL CHECK(role IN ('in','out')),
        PRIMARY KEY(op, artifact, role)
    ) STRICT",
    "CREATE TABLE IF NOT EXISTS metrics(
        op INTEGER NOT NULL,
        t INTEGER NOT NULL,
        name TEXT NOT NULL CHECK(length(name) > 0),
        value REAL NOT NULL,
        PRIMARY KEY(op, t, name)
    ) STRICT",
    "CREATE TABLE IF NOT EXISTS tune(
        kernel TEXT NOT NULL,
        shape_class TEXT NOT NULL,
        machine BLOB NOT NULL,
        params TEXT NOT NULL CHECK(json_valid(params)),
        measured TEXT NOT NULL CHECK(json_valid(measured)),
        PRIMARY KEY(kernel, shape_class, machine)
    ) STRICT",
    "CREATE TABLE IF NOT EXISTS events(
        id INTEGER PRIMARY KEY,
        session BLOB,
        t INTEGER NOT NULL,
        kind TEXT NOT NULL CHECK(length(kind) > 0),
        payload TEXT CHECK(payload IS NULL OR json_valid(payload))
    ) STRICT",
    // -- indexes for the query shapes the plan names ----------------------
    "CREATE INDEX IF NOT EXISTS idx_edges_artifact ON edges(artifact)",
    "CREATE INDEX IF NOT EXISTS idx_events_session_t ON events(session, t)",
    "CREATE INDEX IF NOT EXISTS idx_ops_session ON ops(session)",
    // -- Rev S extension tables (sparse in v0, uniform shape) --------------
    "CREATE TABLE IF NOT EXISTS requirements(
        id INTEGER PRIMARY KEY,
        name TEXT NOT NULL UNIQUE CHECK(length(name) > 0),
        body TEXT NOT NULL CHECK(json_valid(body)),
        created_at INTEGER NOT NULL
    ) STRICT",
    "CREATE TABLE IF NOT EXISTS model_cards(
        id INTEGER PRIMARY KEY,
        name TEXT NOT NULL UNIQUE CHECK(length(name) > 0),
        body TEXT NOT NULL CHECK(json_valid(body)),
        created_at INTEGER NOT NULL
    ) STRICT",
    "CREATE TABLE IF NOT EXISTS evidence(
        id INTEGER PRIMARY KEY,
        name TEXT NOT NULL UNIQUE CHECK(length(name) > 0),
        body TEXT NOT NULL CHECK(json_valid(body)),
        created_at INTEGER NOT NULL
    ) STRICT",
    "CREATE TABLE IF NOT EXISTS scenarios(
        id INTEGER PRIMARY KEY,
        name TEXT NOT NULL UNIQUE CHECK(length(name) > 0),
        body TEXT NOT NULL CHECK(json_valid(body)),
        created_at INTEGER NOT NULL
    ) STRICT",
    "CREATE TABLE IF NOT EXISTS constraints(
        id INTEGER PRIMARY KEY,
        name TEXT NOT NULL UNIQUE CHECK(length(name) > 0),
        body TEXT NOT NULL CHECK(json_valid(body)),
        created_at INTEGER NOT NULL
    ) STRICT",
    "CREATE TABLE IF NOT EXISTS capability_probes(
        id INTEGER PRIMARY KEY,
        name TEXT NOT NULL UNIQUE CHECK(length(name) > 0),
        body TEXT NOT NULL CHECK(json_valid(body)),
        created_at INTEGER NOT NULL
    ) STRICT",
    "CREATE TABLE IF NOT EXISTS imports(
        id INTEGER PRIMARY KEY,
        name TEXT NOT NULL UNIQUE CHECK(length(name) > 0),
        body TEXT NOT NULL CHECK(json_valid(body)),
        created_at INTEGER NOT NULL
    ) STRICT",
    "CREATE TABLE IF NOT EXISTS unsafe_capsules(
        id INTEGER PRIMARY KEY,
        name TEXT NOT NULL UNIQUE CHECK(length(name) > 0),
        body TEXT NOT NULL CHECK(json_valid(body)),
        created_at INTEGER NOT NULL
    ) STRICT",
];

/// v2: forkable worlds and replay provenance (plan §11.2 time travel).
/// `branches` models the op-log branch tree (main = row 1, created here);
/// `ops` gains its branch and the recorded execution mode (replays of
/// `deterministic` ops must reproduce artifact hashes exactly; `fast` ops
/// may diverge and the replay audit reports them separately).
///
/// `ADD COLUMN` keeps the defaults NON-NULL so every pre-v2 op lands on the
/// main branch as a deterministic op — the correct reading of v1 history.
/// The `INSERT ... WHERE NOT EXISTS` seed is idempotent. The two `ADD COLUMN`
/// statements predate atomic version markers, so their exact definitions are
/// also registered in [`RECOVERABLE_ADDED_COLUMNS`] for crash-window healing.
pub(crate) const V2_ADD_BRANCH_COLUMN: &str =
    "ALTER TABLE ops ADD COLUMN branch INTEGER NOT NULL DEFAULT 1";
pub(crate) const V2_ADD_EXEC_MODE_COLUMN: &str =
    "ALTER TABLE ops ADD COLUMN exec_mode TEXT NOT NULL DEFAULT 'deterministic'";

/// Ordered v2 DDL batch.
pub const V2: &[&str] = &[
    "CREATE TABLE IF NOT EXISTS branches(
        id INTEGER PRIMARY KEY,
        name TEXT NOT NULL UNIQUE CHECK(length(name) > 0),
        parent INTEGER,
        fork_op INTEGER,
        created_at INTEGER NOT NULL
    ) STRICT",
    "INSERT INTO branches(id, name, parent, fork_op, created_at)
     SELECT 1, 'main', NULL, NULL, 0
     WHERE NOT EXISTS (SELECT 1 FROM branches WHERE id = 1)",
    V2_ADD_BRANCH_COLUMN,
    V2_ADD_EXEC_MODE_COLUMN,
    "CREATE INDEX IF NOT EXISTS idx_ops_branch ON ops(branch)",
];

/// Exact metadata for a non-idempotent `ADD COLUMN` shipped before migration
/// version markers became transactionally atomic.
#[derive(Debug, Clone, Copy)]
pub(crate) struct RecoverableAddedColumn {
    pub ddl: &'static str,
    pub table: &'static str,
    pub name: &'static str,
    pub declared_type: &'static str,
    pub not_null: bool,
    pub default_sql: Option<&'static str>,
    pub primary_key: bool,
}

/// Columns that may already exist while `user_version` still names the prior
/// schema. Recovery skips an `ALTER` only after every declared property agrees.
pub(crate) const RECOVERABLE_ADDED_COLUMNS: &[RecoverableAddedColumn] = &[
    RecoverableAddedColumn {
        ddl: V2_ADD_BRANCH_COLUMN,
        table: "ops",
        name: "branch",
        declared_type: "INTEGER",
        not_null: true,
        default_sql: Some("1"),
        primary_key: false,
    },
    RecoverableAddedColumn {
        ddl: V2_ADD_EXEC_MODE_COLUMN,
        table: "ops",
        name: "exec_mode",
        declared_type: "TEXT",
        not_null: true,
        default_sql: Some("'deterministic'"),
        primary_key: false,
    },
];

/// Names of every table the v1 schema owns (used by lint and tests).
pub const V1_TABLES: &[&str] = &[
    "artifacts",
    "artifact_chunks",
    "ops",
    "edges",
    "metrics",
    "tune",
    "events",
    "requirements",
    "model_cards",
    "evidence",
    "scenarios",
    "constraints",
    "capability_probes",
    "imports",
    "unsafe_capsules",
];

/// v3 (bead lmp4.3): speculation telemetry — solve nodes gain
/// `(proposer_id, accepted, bound, iterations_saved)` as speculation
/// records keyed by solve-op identity. Additive: every existing query
/// is untouched (the migration regression test proves it).
pub const V3: &[&str] = &["CREATE TABLE IF NOT EXISTS speculation(
        id INTEGER PRIMARY KEY,
        name TEXT NOT NULL UNIQUE CHECK(length(name) > 0),
        body TEXT NOT NULL CHECK(json_valid(body)),
        created_at INTEGER NOT NULL
    ) STRICT"];

/// v4 (bead pifg): one immutable, move-stable identity for the physical
/// ledger instance. File-backed ledgers retain this row across path aliases and
/// reopenings; a replacement database at the same path receives a new value.
/// Independent in-memory handles likewise receive distinct values that live in
/// the handle rather than depending on a movable Rust address.
pub const V4: &[&str] = &["CREATE TABLE IF NOT EXISTS ledger_identity(
    singleton INTEGER PRIMARY KEY CHECK(singleton = 1),
    instance_id BLOB NOT NULL CHECK(length(instance_id) = 16)
) STRICT"];

/// Every table the CURRENT schema owns (v1 set + v2/v3/v4 additions); the
/// `table_count`/lint whitelist.
pub const ALL_TABLES: &[&str] = &[
    "artifacts",
    "artifact_chunks",
    "ops",
    "edges",
    "metrics",
    "tune",
    "events",
    "requirements",
    "model_cards",
    "evidence",
    "scenarios",
    "constraints",
    "capability_probes",
    "imports",
    "unsafe_capsules",
    "branches",
    "speculation",
    "ledger_identity",
];
