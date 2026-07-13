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
//! database-instance identity row, seeded by Rust inside the same migration
//! transaction as the table and version marker. Schema v5 makes that row
//! immutable through attested update/delete refusal triggers. Schema v6 adds
//! immutable terminal-session receipts and deterministic flush-batch markers
//! so retry after a database commit cannot append the same audit events twice.

/// The schema version this crate writes and reads.
pub const SCHEMA_VERSION: i64 = 6;

/// Storage chunk length for large artifacts (bytes). Artifacts strictly
/// larger than this are stored as `artifact_chunks` rows of at most this
/// size; smaller ones live inline in `artifacts.bytes`.
pub const STORAGE_CHUNK_LEN: usize = 4 * 1024 * 1024;

/// Migration ladder: `MIGRATIONS[i]` migrates a database at `user_version`
/// `i` to `i + 1`. Append-only; never edit a shipped batch.
pub(crate) const MIGRATIONS: &[&[&str]] = &[V1, V2, V3, V4, V5, V6];

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

/// v4 (bead pifg): one persisted, move-stable identity for the physical
/// ledger instance. File-backed ledgers retain this row across path aliases and
/// reopenings; a replacement database at the same path receives a new value.
/// Independent in-memory handles likewise receive distinct values that live in
/// the handle rather than depending on a movable Rust address.
pub const V4: &[&str] = &["CREATE TABLE IF NOT EXISTS ledger_identity(
    singleton INTEGER PRIMARY KEY CHECK(singleton = 1),
    instance_id BLOB NOT NULL CHECK(length(instance_id) = 16)
) STRICT"];

/// v5: make the physical-ledger identity immutable under ordinary SQL writes.
/// All guards are schema-attested like every other shipped object. The
/// update guard fires even for a no-op assignment: ledger identity is creation
/// metadata, never mutable application state. The insert guard permits the
/// initial seed only while the singleton row is absent, closing `OR REPLACE`.
pub const V5: &[&str] = &[
    "CREATE TRIGGER IF NOT EXISTS trg_ledger_identity_immutable_update
     BEFORE UPDATE ON ledger_identity
     BEGIN
       SELECT RAISE(ABORT, 'ledger_identity is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_ledger_identity_immutable_delete
     BEFORE DELETE ON ledger_identity
     BEGIN
       SELECT RAISE(ABORT, 'ledger_identity is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_ledger_identity_immutable_reinsert
     BEFORE INSERT ON ledger_identity
     WHEN EXISTS(SELECT 1 FROM ledger_identity WHERE singleton = 1)
     BEGIN
       SELECT RAISE(ABORT, 'ledger_identity is immutable');
     END",
];

/// v6: immutable session-mutation claims, authenticated terminal receipts,
/// explicit ownership links to global audit events, and exact flush-batch
/// membership witnesses. The public API verifies every bounded BLOB/hash and
/// commits each terminal plus its owned event group in one transaction.
pub const V6: &[&str] = &[
    "CREATE TABLE IF NOT EXISTS session_claims(
        authority BLOB NOT NULL PRIMARY KEY CHECK(length(authority) = 32),
        ledger_instance_id BLOB NOT NULL CHECK(length(ledger_instance_id) = 16),
        governor_hash BLOB NOT NULL CHECK(length(governor_hash) = 32),
        session_open_hash BLOB NOT NULL CHECK(length(session_open_hash) = 32),
        registry_schema_version INTEGER NOT NULL CHECK(registry_schema_version = 1),
        kind TEXT NOT NULL CHECK(
            length(CAST(kind AS BLOB)) BETWEEN 1 AND 64 AND
            length(CAST(kind AS BLOB)) = length(kind) AND
            kind NOT GLOB '*[^!-~]*'
        ),
        session BLOB NOT NULL CHECK(length(session) = 8),
        ledger_scope TEXT NOT NULL CHECK(
            length(CAST(ledger_scope AS BLOB)) BETWEEN 1 AND 128 AND
            length(CAST(ledger_scope AS BLOB)) = length(ledger_scope) AND
            ledger_scope NOT GLOB '*[^!-~]*'
        ),
        generation BLOB NOT NULL CHECK(length(generation) = 8),
        causal_ordinal BLOB CHECK(
            causal_ordinal IS NULL OR
            (typeof(causal_ordinal) = 'blob' AND length(causal_ordinal) = 8)
        ),
        payload BLOB NOT NULL CHECK(length(payload) BETWEEN 0 AND 1048576),
        payload_hash BLOB NOT NULL CHECK(length(payload_hash) = 32),
        claim_hash BLOB NOT NULL CHECK(length(claim_hash) = 32),
        created_at INTEGER NOT NULL
    ) STRICT",
    "CREATE TABLE IF NOT EXISTS session_terminals(
        authority BLOB NOT NULL PRIMARY KEY CHECK(length(authority) = 32)
            REFERENCES session_claims(authority),
        receipt BLOB NOT NULL CHECK(length(receipt) BETWEEN 1 AND 1048576),
        receipt_hash BLOB NOT NULL CHECK(length(receipt_hash) = 32),
        event_count INTEGER NOT NULL CHECK(event_count BETWEEN 0 AND 1024),
        events_hash BLOB NOT NULL CHECK(length(events_hash) = 32),
        encoded_bytes INTEGER NOT NULL CHECK(encoded_bytes BETWEEN 1 AND 4194304),
        created_at INTEGER NOT NULL
    ) STRICT",
    "CREATE TABLE IF NOT EXISTS session_terminal_events(
        authority BLOB NOT NULL CHECK(length(authority) = 32)
            REFERENCES session_terminals(authority),
        seq INTEGER NOT NULL CHECK(seq BETWEEN 0 AND 1023),
        event_id INTEGER NOT NULL UNIQUE CHECK(event_id > 0) REFERENCES events(id),
        PRIMARY KEY(authority, seq)
    ) STRICT",
    "CREATE TABLE IF NOT EXISTS session_flush_batches(
        batch_id BLOB NOT NULL PRIMARY KEY CHECK(length(batch_id) = 32),
        ledger_instance_id BLOB NOT NULL CHECK(length(ledger_instance_id) = 16),
        registry_schema_version INTEGER NOT NULL CHECK(registry_schema_version = 1),
        terminal_count INTEGER NOT NULL CHECK(terminal_count BETWEEN 1 AND 1024),
        event_count INTEGER NOT NULL CHECK(event_count BETWEEN 0 AND 1024),
        encoded_bytes INTEGER NOT NULL CHECK(encoded_bytes BETWEEN 1 AND 4194304),
        created_at INTEGER NOT NULL
    ) STRICT",
    "CREATE TABLE IF NOT EXISTS session_flush_batch_members(
        batch_id BLOB NOT NULL CHECK(length(batch_id) = 32)
            REFERENCES session_flush_batches(batch_id),
        seq INTEGER NOT NULL CHECK(seq BETWEEN 0 AND 1023),
        authority BLOB NOT NULL CHECK(length(authority) = 32)
            REFERENCES session_terminals(authority),
        PRIMARY KEY(batch_id, seq),
        UNIQUE(batch_id, authority)
    ) STRICT",
    "CREATE INDEX IF NOT EXISTS idx_session_flush_batch_members_authority
     ON session_flush_batch_members(authority)",
    "CREATE TRIGGER IF NOT EXISTS trg_session_claims_immutable_update
     BEFORE UPDATE ON session_claims
     BEGIN
       SELECT RAISE(ABORT, 'session claim is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_session_claims_immutable_delete
     BEFORE DELETE ON session_claims
     BEGIN
       SELECT RAISE(ABORT, 'session claim is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_session_claims_immutable_reinsert
     BEFORE INSERT ON session_claims
     WHEN EXISTS(SELECT 1 FROM session_claims WHERE authority = NEW.authority)
     BEGIN
       SELECT RAISE(ABORT, 'session claim is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_session_terminals_immutable_update
     BEFORE UPDATE ON session_terminals
     BEGIN
       SELECT RAISE(ABORT, 'session terminal is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_session_terminals_immutable_delete
     BEFORE DELETE ON session_terminals
     BEGIN
       SELECT RAISE(ABORT, 'session terminal is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_session_terminals_immutable_reinsert
     BEFORE INSERT ON session_terminals
     WHEN EXISTS(
         SELECT 1 FROM session_terminals WHERE authority = NEW.authority
     )
     BEGIN
       SELECT RAISE(ABORT, 'session terminal is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_session_terminal_events_immutable_update
     BEFORE UPDATE ON session_terminal_events
     BEGIN
       SELECT RAISE(ABORT, 'session terminal event link is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_session_terminal_events_immutable_delete
     BEFORE DELETE ON session_terminal_events
     BEGIN
       SELECT RAISE(ABORT, 'session terminal event link is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_session_terminal_events_immutable_reinsert
     BEFORE INSERT ON session_terminal_events
     WHEN EXISTS(
         SELECT 1 FROM session_terminal_events
         WHERE (authority = NEW.authority AND seq = NEW.seq)
            OR event_id = NEW.event_id
     )
     BEGIN
       SELECT RAISE(ABORT, 'session terminal event link is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_session_flush_batches_immutable_update
     BEFORE UPDATE ON session_flush_batches
     BEGIN
       SELECT RAISE(ABORT, 'session flush batch is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_session_flush_batches_immutable_delete
     BEFORE DELETE ON session_flush_batches
     BEGIN
       SELECT RAISE(ABORT, 'session flush batch is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_session_flush_batches_immutable_reinsert
     BEFORE INSERT ON session_flush_batches
     WHEN EXISTS(
         SELECT 1 FROM session_flush_batches WHERE batch_id = NEW.batch_id
     )
     BEGIN
       SELECT RAISE(ABORT, 'session flush batch is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_session_flush_batch_members_immutable_update
     BEFORE UPDATE ON session_flush_batch_members
     BEGIN
       SELECT RAISE(ABORT, 'session flush batch member is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_session_flush_batch_members_immutable_delete
     BEFORE DELETE ON session_flush_batch_members
     BEGIN
       SELECT RAISE(ABORT, 'session flush batch member is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_session_flush_batch_members_immutable_reinsert
     BEFORE INSERT ON session_flush_batch_members
     WHEN EXISTS(
         SELECT 1 FROM session_flush_batch_members
         WHERE (batch_id = NEW.batch_id AND seq = NEW.seq)
            OR (batch_id = NEW.batch_id AND authority = NEW.authority)
     )
     BEGIN
       SELECT RAISE(ABORT, 'session flush batch member is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_owned_session_events_immutable_update
     BEFORE UPDATE ON events
     WHEN EXISTS(
         SELECT 1 FROM session_terminal_events WHERE event_id = OLD.id
     )
     BEGIN
       SELECT RAISE(ABORT, 'owned session event is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_owned_session_events_immutable_delete
     BEFORE DELETE ON events
     WHEN EXISTS(
         SELECT 1 FROM session_terminal_events WHERE event_id = OLD.id
     )
     BEGIN
       SELECT RAISE(ABORT, 'owned session event is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_owned_session_events_immutable_reinsert
     BEFORE INSERT ON events
     WHEN EXISTS(
         SELECT 1 FROM session_terminal_events WHERE event_id = NEW.id
     )
     BEGIN
       SELECT RAISE(ABORT, 'owned session event is immutable');
     END",
];

/// Every table the CURRENT schema owns (v1 set + v2 through v6 additions); the
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
    "session_claims",
    "session_terminals",
    "session_terminal_events",
    "session_flush_batches",
    "session_flush_batch_members",
];
