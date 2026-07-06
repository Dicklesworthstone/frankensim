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
//! Migrations are versioned through `PRAGMA user_version`; every DDL batch is
//! idempotent (`IF NOT EXISTS`) so a crash between DDL and the version bump
//! re-applies harmlessly on the next open.

/// The schema version this crate writes and reads.
pub const SCHEMA_VERSION: i64 = 2;

/// Storage chunk length for large artifacts (bytes). Artifacts strictly
/// larger than this are stored as `artifact_chunks` rows of at most this
/// size; smaller ones live inline in `artifacts.bytes`.
pub const STORAGE_CHUNK_LEN: usize = 4 * 1024 * 1024;

/// Migration ladder: `MIGRATIONS[i]` migrates a database at `user_version`
/// `i` to `i + 1`. Append-only; never edit a shipped batch.
pub(crate) const MIGRATIONS: &[&[&str]] = &[V1, V2];

/// v1: the six core tables (Appendix D), chunk storage, and the Rev S
/// extension tables (sparse in v0 but present EARLY so downstream crates can
/// rely on them existing).
const V1: &[&str] = &[
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
/// The `INSERT ... WHERE NOT EXISTS` seed is idempotent (crash between DDL
/// and the version bump re-applies harmlessly, like every batch here).
const V2: &[&str] = &[
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
    "ALTER TABLE ops ADD COLUMN branch INTEGER NOT NULL DEFAULT 1",
    "ALTER TABLE ops ADD COLUMN exec_mode TEXT NOT NULL DEFAULT 'deterministic'",
    "CREATE INDEX IF NOT EXISTS idx_ops_branch ON ops(branch)",
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

/// Every table the CURRENT schema owns (v1 set + v2 additions); the
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
];
