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
//! Schema v7 adds causal-ordinal ownership and insert guards without rewriting
//! the shipped v6 tables. Schema v8 adds an immutable, independently indexed
//! discovery witness for session claims and splits OR-based reinsert guards so
//! each refusal probe follows one existing unique index. Schema v9 adds the
//! two covering lineage indexes used by capped verifier reads and immutable
//! per-artifact output seals for consumers that require exactly one producer,
//! plus immutable operation-edge-set seals for exact-lineage consumers.
//! Schema v10 adds immutable solver-checkpoint receipts that bind one pause
//! authority to an existing solver-state artifact and executor-originated
//! drain/finalize report. Schema v11 adds immutable five-to-six quantity
//! crosswalk rows whose source and target hashes both reference retained
//! artifacts. Schema v12 adds portable semantic state-checkpoint receipts that
//! bind runtime bytes to exact law, parameter, schema, and implementation
//! identities. Schema v13 adds immutable strong-identity migration receipts:
//! exact legacy and canonical bytes remain distinct from raw content IDs,
//! typed semantic IDs, legacy FNV replay tokens, and authority state. Schema
//! v14 gives every artifact an exact typed content-identity companion row;
//! migration backfills existing hashes and a trigger dual-writes later rows
//! without guessing a semantic schema. Schema v15 carries that typed raw-byte
//! identity across every role-qualified operation/artifact lineage edge.
//! Schema v16 adds explicit, receipt-backed artifact semantic bindings while
//! leaving all unknown historical meanings unbound.
//! Schema v17 adds the same explicit receipt boundary for exact retained
//! evidence JSON bytes and makes a bound evidence body immutable.

/// The schema version this crate writes and reads.
pub const SCHEMA_VERSION: i64 = 17;

/// Storage chunk length for large artifacts (bytes). Artifacts strictly
/// larger than this are stored as `artifact_chunks` rows of at most this
/// size; smaller ones live inline in `artifacts.bytes`.
pub const STORAGE_CHUNK_LEN: usize = 4 * 1024 * 1024;

/// Migration ladder: `MIGRATIONS[i]` migrates a database at `user_version`
/// `i` to `i + 1`. Append-only; never edit a shipped batch.
pub(crate) const MIGRATIONS: &[&[&str]] = &[
    V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13, V14, V15, V16, V17,
];

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
     WHEN EXISTS (SELECT 1 FROM ledger_identity WHERE singleton = 1)
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
    "CREATE INDEX IF NOT EXISTS idx_session_claims_recovery_pending
     ON session_claims(
         governor_hash, session_open_hash, kind, session, ledger_scope, generation, authority
     )",
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

/// v7: strengthen the shipped v6 session registry without rebuilding tables.
/// New causal ordinals are canonical positive signed-ledger values, one
/// governor/kind ordinal has one owner. Legacy v6 submission claims may have a NULL
/// admission ordinal; typed recovery reads the ordinal from their authenticated
/// receipt, while the insert guard requires every new submission to bind it.
pub const V7: &[&str] = &[
    "CREATE INDEX IF NOT EXISTS idx_session_claims_governor_authority
     ON session_claims(governor_hash, authority)",
    "CREATE UNIQUE INDEX IF NOT EXISTS idx_session_claims_unique_causal_ordinal
     ON session_claims(governor_hash, kind, causal_ordinal)
     WHERE causal_ordinal IS NOT NULL",
    "CREATE TRIGGER IF NOT EXISTS trg_session_claims_causal_ordinal_range
     BEFORE INSERT ON session_claims
     WHEN NEW.causal_ordinal IS NOT NULL AND (
          typeof(NEW.causal_ordinal) != 'blob' OR
          length(NEW.causal_ordinal) != 8 OR
          NEW.causal_ordinal <= X'0000000000000000' OR
          NEW.causal_ordinal > X'7FFFFFFFFFFFFFFF'
     )
     BEGIN
       SELECT RAISE(ABORT, 'session causal ordinal is outside 1..=i64::MAX');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_session_submission_requires_admission_ordinal
     BEFORE INSERT ON session_claims
     WHEN NEW.kind = 'submission' AND NEW.causal_ordinal IS NULL
     BEGIN
       SELECT RAISE(ABORT, 'session submission requires admission ordinal');
    END",
];

/// v8: add a compact second copy of the authenticated claim-discovery fields.
/// Recovery can scan this witness without repeatedly decoding the potentially
/// large claim payload, then compare every field with `session_claims` before
/// trusting it. Existing claims are copied before the new table becomes
/// immutable. The two v6 reinsert guards that used `OR` are also replaced by
/// one indexable trigger per uniqueness constraint.
pub const V8: &[&str] = &[
    "CREATE TABLE IF NOT EXISTS session_claim_discovery(
        authority BLOB NOT NULL PRIMARY KEY CHECK(length(authority) = 32)
            REFERENCES session_claims(authority),
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
        payload_hash BLOB NOT NULL CHECK(length(payload_hash) = 32),
        claim_hash BLOB NOT NULL CHECK(length(claim_hash) = 32)
    ) STRICT",
    "INSERT INTO session_claim_discovery(
         authority, ledger_instance_id, governor_hash, session_open_hash,
         registry_schema_version, kind, session, ledger_scope, generation,
         causal_ordinal, payload_hash, claim_hash
     )
     SELECT authority, ledger_instance_id, governor_hash, session_open_hash,
            registry_schema_version, kind, session, ledger_scope, generation,
            causal_ordinal, payload_hash, claim_hash
     FROM session_claims AS claim
     WHERE NOT EXISTS(
         SELECT 1 FROM session_claim_discovery AS discovery
         WHERE discovery.authority = claim.authority
     )
     ORDER BY authority",
    "CREATE INDEX IF NOT EXISTS idx_session_claim_discovery_recovery_pending
     ON session_claim_discovery(
         governor_hash, session_open_hash, kind, session, ledger_scope, generation, authority
     )",
    "CREATE INDEX IF NOT EXISTS idx_session_claim_discovery_governor_authority
     ON session_claim_discovery(governor_hash, authority)",
    "CREATE UNIQUE INDEX IF NOT EXISTS idx_session_claim_discovery_governor_kind_ordinal
     ON session_claim_discovery(governor_hash, kind, causal_ordinal)
     WHERE causal_ordinal IS NOT NULL",
    "CREATE TRIGGER IF NOT EXISTS trg_session_claim_discovery_causal_ordinal_range
     BEFORE INSERT ON session_claim_discovery
     WHEN NEW.causal_ordinal IS NOT NULL AND (
          typeof(NEW.causal_ordinal) != 'blob' OR
          length(NEW.causal_ordinal) != 8 OR
          NEW.causal_ordinal <= X'0000000000000000' OR
          NEW.causal_ordinal > X'7FFFFFFFFFFFFFFF'
     )
     BEGIN
       SELECT RAISE(ABORT, 'session causal ordinal is outside 1..=i64::MAX');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_session_claim_discovery_submission_requires_admission_ordinal
     BEFORE INSERT ON session_claim_discovery
     WHEN NEW.kind = 'submission' AND NEW.causal_ordinal IS NULL
     BEGIN
       SELECT RAISE(ABORT, 'session submission requires admission ordinal');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_session_claim_discovery_immutable_update
     BEFORE UPDATE ON session_claim_discovery
     BEGIN
       SELECT RAISE(ABORT, 'session claim discovery witness is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_session_claim_discovery_immutable_delete
     BEFORE DELETE ON session_claim_discovery
     BEGIN
       SELECT RAISE(ABORT, 'session claim discovery witness is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_session_claim_discovery_immutable_reinsert
     BEFORE INSERT ON session_claim_discovery
     WHEN EXISTS(
         SELECT 1 FROM session_claim_discovery WHERE authority = NEW.authority
     )
     BEGIN
       SELECT RAISE(ABORT, 'session claim discovery witness is immutable');
     END",
    "DROP TRIGGER IF EXISTS trg_session_terminal_events_immutable_reinsert",
    "CREATE TRIGGER IF NOT EXISTS trg_session_terminal_events_immutable_reinsert
     BEFORE INSERT ON session_terminal_events
     WHEN EXISTS(
         SELECT 1 FROM session_terminal_events
         WHERE authority = NEW.authority AND seq = NEW.seq
     )
     BEGIN
       SELECT RAISE(ABORT, 'session terminal event link is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_session_terminal_events_immutable_event_id_reinsert
     BEFORE INSERT ON session_terminal_events
     WHEN EXISTS(
         SELECT 1 FROM session_terminal_events WHERE event_id = NEW.event_id
     )
     BEGIN
       SELECT RAISE(ABORT, 'session terminal event link is immutable');
     END",
    "DROP TRIGGER IF EXISTS trg_session_flush_batch_members_immutable_reinsert",
    "CREATE TRIGGER IF NOT EXISTS trg_session_flush_batch_members_immutable_reinsert
     BEFORE INSERT ON session_flush_batch_members
     WHEN EXISTS(
         SELECT 1 FROM session_flush_batch_members
         WHERE batch_id = NEW.batch_id AND seq = NEW.seq
     )
     BEGIN
       SELECT RAISE(ABORT, 'session flush batch member is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_session_flush_batch_members_immutable_authority_reinsert
     BEFORE INSERT ON session_flush_batch_members
     WHEN EXISTS(
         SELECT 1 FROM session_flush_batch_members
         WHERE batch_id = NEW.batch_id AND authority = NEW.authority
     )
     BEGIN
       SELECT RAISE(ABORT, 'session flush batch member is immutable');
     END",
];

/// v9: make verifier-side lineage caps bound SQLite work as well as returned
/// rows, and provide a durable generic seal for artifacts whose output
/// provenance must remain single-producer. Operation-edge seals independently
/// freeze one bounded, already-validated edge set. Once present, the attested
/// guards reject every conflicting producer or edge-set mutation and every
/// attempt to rewrite or remove either seal.
pub const V9: &[&str] = &[
    "CREATE INDEX IF NOT EXISTS idx_edges_artifact_role_op
     ON edges(artifact, role, op)",
    "CREATE INDEX IF NOT EXISTS idx_edges_op_role_artifact
     ON edges(op, role, artifact)",
    "CREATE TABLE IF NOT EXISTS artifact_output_seals(
        artifact BLOB NOT NULL PRIMARY KEY CHECK(length(artifact) = 32),
        op INTEGER NOT NULL,
        role TEXT NOT NULL CHECK(role = 'out'),
        FOREIGN KEY(op, artifact, role) REFERENCES edges(op, artifact, role)
    ) STRICT",
    "CREATE TABLE IF NOT EXISTS op_artifact_edge_seals(
        op INTEGER NOT NULL PRIMARY KEY REFERENCES ops(id),
        edge_count INTEGER NOT NULL CHECK(edge_count BETWEEN 0 AND 1024)
    ) STRICT",
    "CREATE TRIGGER IF NOT EXISTS trg_artifact_output_seals_exact_producer
     BEFORE INSERT ON artifact_output_seals
     WHEN NOT EXISTS(
              SELECT 1 FROM edges INDEXED BY idx_edges_artifact_role_op
              WHERE artifact = NEW.artifact AND role = 'out' AND op = NEW.op
              LIMIT 1
          ) OR EXISTS(
              SELECT 1 FROM edges INDEXED BY idx_edges_artifact_role_op
              WHERE artifact = NEW.artifact AND role = 'out' AND op != NEW.op
              LIMIT 1
          )
     BEGIN
       SELECT RAISE(ABORT, 'artifact output seal requires one exact producer');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_artifact_output_seals_immutable_update
     BEFORE UPDATE ON artifact_output_seals
     BEGIN
       SELECT RAISE(ABORT, 'artifact output seal is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_artifact_output_seals_immutable_delete
     BEFORE DELETE ON artifact_output_seals
     BEGIN
       SELECT RAISE(ABORT, 'artifact output seal is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_artifact_output_seals_immutable_reinsert
     BEFORE INSERT ON artifact_output_seals
     WHEN EXISTS(
         SELECT 1 FROM artifact_output_seals WHERE artifact = NEW.artifact
     )
     BEGIN
       SELECT RAISE(ABORT, 'artifact output seal is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_op_artifact_edge_seals_exact_count
     BEFORE INSERT ON op_artifact_edge_seals
     WHEN NOT EXISTS(SELECT 1 FROM ops WHERE id = NEW.op LIMIT 1) OR
          NEW.edge_count != (
              SELECT COUNT(*) FROM (
                  SELECT 1 FROM edges INDEXED BY idx_edges_op_role_artifact
                  WHERE op = NEW.op LIMIT 1025
              )
          )
     BEGIN
       SELECT RAISE(ABORT, 'op artifact-edge seal requires the exact bounded edge count');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_op_artifact_edge_seals_immutable_update
     BEFORE UPDATE ON op_artifact_edge_seals
     BEGIN
       SELECT RAISE(ABORT, 'op artifact-edge seal is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_op_artifact_edge_seals_immutable_delete
     BEFORE DELETE ON op_artifact_edge_seals
     BEGIN
       SELECT RAISE(ABORT, 'op artifact-edge seal is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_op_artifact_edge_seals_immutable_reinsert
     BEFORE INSERT ON op_artifact_edge_seals
     WHEN EXISTS(SELECT 1 FROM op_artifact_edge_seals WHERE op = NEW.op)
     BEGIN
       SELECT RAISE(ABORT, 'op artifact-edge seal is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_edges_sealed_output_insert
     BEFORE INSERT ON edges
     WHEN NEW.role = 'out' AND EXISTS(
         SELECT 1 FROM artifact_output_seals
         WHERE artifact = NEW.artifact AND op != NEW.op
     )
     BEGIN
       SELECT RAISE(ABORT, 'sealed artifact rejects a different output producer');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_edges_sealed_output_update
     BEFORE UPDATE ON edges
     WHEN NEW.role = 'out' AND EXISTS(
         SELECT 1 FROM artifact_output_seals
         WHERE artifact = NEW.artifact AND op != NEW.op
     )
     BEGIN
       SELECT RAISE(ABORT, 'sealed artifact rejects a different output producer');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_edges_sealed_output_update_existing
     BEFORE UPDATE ON edges
     WHEN OLD.role = 'out' AND EXISTS(
         SELECT 1 FROM artifact_output_seals
         WHERE artifact = OLD.artifact AND op = OLD.op
     )
     BEGIN
       SELECT RAISE(ABORT, 'sealed artifact output edge is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_edges_sealed_output_delete
     BEFORE DELETE ON edges
     WHEN OLD.role = 'out' AND EXISTS(
         SELECT 1 FROM artifact_output_seals
         WHERE artifact = OLD.artifact AND op = OLD.op
     )
     BEGIN
       SELECT RAISE(ABORT, 'sealed artifact output edge is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_edges_sealed_op_insert
     BEFORE INSERT ON edges
     WHEN EXISTS(SELECT 1 FROM op_artifact_edge_seals WHERE op = NEW.op)
     BEGIN
       SELECT RAISE(ABORT, 'sealed operation artifact-edge set is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_edges_sealed_op_update
     BEFORE UPDATE ON edges
     WHEN EXISTS(SELECT 1 FROM op_artifact_edge_seals WHERE op = OLD.op) OR
          EXISTS(SELECT 1 FROM op_artifact_edge_seals WHERE op = NEW.op)
     BEGIN
       SELECT RAISE(ABORT, 'sealed operation artifact-edge set is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_edges_sealed_op_delete
     BEFORE DELETE ON edges
     WHEN EXISTS(SELECT 1 FROM op_artifact_edge_seals WHERE op = OLD.op)
     BEGIN
       SELECT RAISE(ABORT, 'sealed operation artifact-edge set is immutable');
     END",
];

/// v10: one immutable, idempotent solver-checkpoint receipt per pause
/// authority. Every scalar uses a fixed-width BLOB so the full `u64` domain is
/// represented exactly; the artifact foreign key makes persisted solver state
/// a prerequisite rather than a caller-authored hash claim.
pub const V10: &[&str] = &[
    "CREATE TABLE IF NOT EXISTS session_checkpoint_receipts(
        receipt_hash BLOB NOT NULL PRIMARY KEY CHECK(length(receipt_hash) = 32),
        ledger_instance_id BLOB NOT NULL CHECK(length(ledger_instance_id) = 16),
        session BLOB NOT NULL CHECK(length(session) = 8),
        run BLOB NOT NULL CHECK(length(run) = 8),
        pause_authority BLOB NOT NULL UNIQUE CHECK(length(pause_authority) = 32),
        gate_generation BLOB NOT NULL CHECK(length(gate_generation) = 8),
        solver_state_artifact BLOB NOT NULL CHECK(length(solver_state_artifact) = 32)
            REFERENCES artifacts(hash),
        drain_report_hash BLOB NOT NULL CHECK(length(drain_report_hash) = 32),
        registered_workers BLOB NOT NULL CHECK(length(registered_workers) = 8),
        drained_workers BLOB NOT NULL CHECK(length(drained_workers) = 8),
        created_at INTEGER NOT NULL,
        CHECK(registered_workers = drained_workers),
        CHECK(registered_workers != X'0000000000000000')
    ) STRICT",
    "CREATE INDEX IF NOT EXISTS idx_session_checkpoint_artifact
     ON session_checkpoint_receipts(solver_state_artifact)",
    "CREATE TRIGGER IF NOT EXISTS trg_session_checkpoint_receipts_immutable_update
     BEFORE UPDATE ON session_checkpoint_receipts
     BEGIN
       SELECT RAISE(ABORT, 'session checkpoint receipt is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_session_checkpoint_receipts_immutable_delete
     BEFORE DELETE ON session_checkpoint_receipts
     BEGIN
       SELECT RAISE(ABORT, 'session checkpoint receipt is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_session_checkpoint_receipts_immutable_receipt_reinsert
     BEFORE INSERT ON session_checkpoint_receipts
     WHEN EXISTS(
         SELECT 1 FROM session_checkpoint_receipts
         WHERE receipt_hash = NEW.receipt_hash
     )
     BEGIN
       SELECT RAISE(ABORT, 'session checkpoint receipt is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_session_checkpoint_receipts_immutable_pause_reinsert
     BEFORE INSERT ON session_checkpoint_receipts
     WHEN EXISTS(
         SELECT 1 FROM session_checkpoint_receipts
         WHERE pause_authority = NEW.pause_authority
     )
     BEGIN
       SELECT RAISE(ABORT, 'session checkpoint pause authority is immutable');
     END",
];

/// v11: immutable exact-byte crosswalks from the historical five-base
/// quantity JSON schema to the canonical six-base schema. Migration infers no
/// rows: only the typed fs-qty receipt plus both retained artifacts may create
/// one through the public ledger API.
pub const V11: &[&str] = &[
    "CREATE TABLE IF NOT EXISTS qty_dimension_crosswalks(
        old_hash BLOB NOT NULL PRIMARY KEY CHECK(length(old_hash) = 32)
            REFERENCES artifacts(hash),
        new_hash BLOB NOT NULL CHECK(length(new_hash) = 32)
            REFERENCES artifacts(hash),
        source_version INTEGER NOT NULL CHECK(source_version = 1),
        target_version INTEGER NOT NULL CHECK(target_version = 2),
        rule TEXT NOT NULL CHECK(rule = 'append-mole-zero'),
        created_at INTEGER NOT NULL,
        CHECK(old_hash != new_hash)
    ) STRICT",
    "CREATE INDEX IF NOT EXISTS idx_qty_dimension_crosswalks_new_hash
     ON qty_dimension_crosswalks(new_hash)",
    "CREATE TRIGGER IF NOT EXISTS trg_qty_dimension_crosswalks_immutable_update
     BEFORE UPDATE ON qty_dimension_crosswalks
     BEGIN
       SELECT RAISE(ABORT, 'quantity dimension crosswalk is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_qty_dimension_crosswalks_immutable_delete
     BEFORE DELETE ON qty_dimension_crosswalks
     BEGIN
       SELECT RAISE(ABORT, 'quantity dimension crosswalk is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_qty_dimension_crosswalks_immutable_reinsert
     BEFORE INSERT ON qty_dimension_crosswalks
     WHEN EXISTS(
         SELECT 1 FROM qty_dimension_crosswalks
         WHERE old_hash = NEW.old_hash
     )
     BEGIN
       SELECT RAISE(ABORT, 'quantity dimension crosswalk source is immutable');
     END",
];

/// v12: portable semantic state-checkpoint receipts. The runtime-state bytes
/// are retained as an artifact; every other field is an immutable semantic
/// identity that a replayer must know exactly before decoding those bytes.
/// Migration infers no rows from generic solver snapshots because they do not
/// carry this complete binding.
pub const V12: &[&str] = &[
    "CREATE TABLE IF NOT EXISTS semantic_state_checkpoint_receipts(
        receipt_hash BLOB NOT NULL PRIMARY KEY CHECK(length(receipt_hash) = 32),
        state_slot BLOB NOT NULL
            CHECK(length(state_slot) = 32 AND state_slot != X'0000000000000000000000000000000000000000000000000000000000000000'),
        law_id BLOB NOT NULL CHECK(length(law_id) BETWEEN 1 AND 256),
        law_version INTEGER NOT NULL CHECK(law_version BETWEEN 0 AND 4294967295),
        state_schema_version INTEGER NOT NULL
            CHECK(state_schema_version BETWEEN 0 AND 4294967295),
        runtime_state_artifact BLOB NOT NULL CHECK(length(runtime_state_artifact) = 32)
            REFERENCES artifacts(hash),
        canonical_parameters_hash BLOB NOT NULL
            CHECK(length(canonical_parameters_hash) = 32
                AND canonical_parameters_hash != X'0000000000000000000000000000000000000000000000000000000000000000'),
        contract_and_code_hash BLOB NOT NULL
            CHECK(length(contract_and_code_hash) = 32
                AND contract_and_code_hash != X'0000000000000000000000000000000000000000000000000000000000000000'),
        created_at INTEGER NOT NULL
    ) STRICT",
    "CREATE INDEX IF NOT EXISTS idx_semantic_state_checkpoint_runtime_artifact
     ON semantic_state_checkpoint_receipts(runtime_state_artifact)",
    "CREATE INDEX IF NOT EXISTS idx_semantic_state_checkpoint_slot
     ON semantic_state_checkpoint_receipts(state_slot)",
    "CREATE INDEX IF NOT EXISTS idx_semantic_state_checkpoint_law
     ON semantic_state_checkpoint_receipts(law_id, law_version, state_schema_version)",
    "CREATE TRIGGER IF NOT EXISTS trg_semantic_state_checkpoint_immutable_update
     BEFORE UPDATE ON semantic_state_checkpoint_receipts
     BEGIN
       SELECT RAISE(ABORT, 'semantic state checkpoint receipt is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_semantic_state_checkpoint_immutable_delete
     BEFORE DELETE ON semantic_state_checkpoint_receipts
     BEGIN
       SELECT RAISE(ABORT, 'semantic state checkpoint receipt is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_semantic_state_checkpoint_immutable_receipt_reinsert
     BEFORE INSERT ON semantic_state_checkpoint_receipts
     WHEN EXISTS(
         SELECT 1 FROM semantic_state_checkpoint_receipts
         WHERE receipt_hash = NEW.receipt_hash
     )
     BEGIN
       SELECT RAISE(ABORT, 'semantic state checkpoint receipt is immutable');
     END",
];

/// v13: immutable crosswalk receipts for the strong-identity migration.
/// Exact legacy bytes and their quarantined FNV token remain separate from
/// canonical bytes, plain content IDs, typed semantic identity, schema
/// identity, and authority state. Multiple receipts may name one legacy
/// content ID; callers must inspect the bounded candidate set and explicitly
/// select a typed schema instead of treating row presence as authority.
pub const V13: &[&str] = &[
    "CREATE TABLE IF NOT EXISTS identity_migration_receipts(
        receipt_id BLOB NOT NULL PRIMARY KEY CHECK(length(receipt_id) = 32),
        legacy_bytes BLOB NOT NULL CHECK(length(legacy_bytes) <= 1048576),
        legacy_content_id BLOB NOT NULL CHECK(length(legacy_content_id) = 32),
        legacy_fnv BLOB NOT NULL CHECK(length(legacy_fnv) = 8),
        canonical_bytes BLOB NOT NULL CHECK(length(canonical_bytes) <= 1048576),
        canonical_content_id BLOB NOT NULL CHECK(length(canonical_content_id) = 32),
        semantic_rule BLOB NOT NULL CHECK(length(semantic_rule) BETWEEN 1 AND 256),
        semantic_id BLOB NOT NULL CHECK(length(semantic_id) = 32),
        identity_role INTEGER NOT NULL CHECK(identity_role BETWEEN 1 AND 12),
        identity_domain BLOB NOT NULL CHECK(length(identity_domain) BETWEEN 1 AND 256),
        identity_schema_name BLOB NOT NULL
            CHECK(length(identity_schema_name) BETWEEN 1 AND 256),
        identity_schema_id BLOB NOT NULL CHECK(length(identity_schema_id) = 32),
        identity_schema_version INTEGER NOT NULL
            CHECK(identity_schema_version BETWEEN 1 AND 4294967295),
        identity_context BLOB NOT NULL CHECK(length(identity_context) BETWEEN 1 AND 4096),
        canonical_preimage_id BLOB NOT NULL CHECK(length(canonical_preimage_id) = 32),
        canonical_frame_bytes BLOB NOT NULL CHECK(length(canonical_frame_bytes) = 8),
        field_count INTEGER NOT NULL CHECK(field_count BETWEEN 0 AND 4294967295),
        collection_items BLOB NOT NULL CHECK(length(collection_items) = 8),
        max_canonical_bytes BLOB NOT NULL CHECK(length(max_canonical_bytes) = 8),
        max_field_bytes BLOB NOT NULL CHECK(length(max_field_bytes) = 8),
        max_fields INTEGER NOT NULL CHECK(max_fields BETWEEN 1 AND 4294967295),
        max_collection_items BLOB NOT NULL CHECK(length(max_collection_items) = 8),
        cancellation_poll_bytes INTEGER NOT NULL
            CHECK(cancellation_poll_bytes BETWEEN 1 AND 4294967295),
        trust_state INTEGER NOT NULL CHECK(trust_state BETWEEN 0 AND 3),
        anchor_content_id BLOB CHECK(anchor_content_id IS NULL OR length(anchor_content_id) = 32),
        verifier_id BLOB CHECK(verifier_id IS NULL OR length(verifier_id) = 32),
        key_policy_id BLOB CHECK(key_policy_id IS NULL OR length(key_policy_id) = 32),
        no_claim_state INTEGER NOT NULL CHECK(no_claim_state BETWEEN 0 AND 1),
        created_at INTEGER NOT NULL,
        CHECK(
            (trust_state = 0 AND anchor_content_id IS NULL AND verifier_id IS NULL
                AND key_policy_id IS NULL AND no_claim_state = 0)
            OR
            (trust_state IN (1, 2) AND anchor_content_id IS NOT NULL
                AND verifier_id IS NOT NULL AND key_policy_id IS NOT NULL
                AND no_claim_state = 0)
            OR
            (trust_state = 3 AND anchor_content_id IS NOT NULL
                AND verifier_id IS NOT NULL AND key_policy_id IS NOT NULL
                AND no_claim_state = 1)
        )
    ) STRICT",
    "CREATE INDEX IF NOT EXISTS idx_identity_migration_legacy
     ON identity_migration_receipts(legacy_content_id, receipt_id)",
    "CREATE INDEX IF NOT EXISTS idx_identity_migration_canonical
     ON identity_migration_receipts(canonical_content_id, receipt_id)",
    "CREATE INDEX IF NOT EXISTS idx_identity_migration_semantic
     ON identity_migration_receipts(
         identity_role, identity_domain, identity_schema_version, semantic_id, receipt_id
     )",
    "CREATE TRIGGER IF NOT EXISTS trg_identity_migration_receipts_immutable_update
     BEFORE UPDATE ON identity_migration_receipts
     BEGIN
       SELECT RAISE(ABORT, 'identity migration receipt is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_identity_migration_receipts_immutable_delete
     BEFORE DELETE ON identity_migration_receipts
     BEGIN
       SELECT RAISE(ABORT, 'identity migration receipt is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_identity_migration_receipts_immutable_reinsert
     BEFORE INSERT ON identity_migration_receipts
     WHEN EXISTS(
         SELECT 1 FROM identity_migration_receipts WHERE receipt_id = NEW.receipt_id
     )
     BEGIN
       SELECT RAISE(ABORT, 'identity migration receipt is immutable');
     END",
];

/// v14: typed content identity for every artifact compatibility hash.
///
/// `artifacts.hash` remains the compatibility key used by schema-v1 readers.
/// The companion table names the same exact digest as a raw-byte `ContentId`;
/// equality is enforced in SQL, not inferred by readers. The backfill copies
/// no semantic identity or authority. An `AFTER INSERT` trigger makes writes
/// from already-open old handles dual-write once this migration commits.
pub const V14: &[&str] = &[
    "CREATE TABLE IF NOT EXISTS artifact_content_identities(
        artifact_hash BLOB NOT NULL PRIMARY KEY
            REFERENCES artifacts(hash) ON DELETE CASCADE
            CHECK(length(artifact_hash) = 32),
        content_id BLOB NOT NULL CHECK(length(content_id) = 32),
        row_schema_version INTEGER NOT NULL CHECK(row_schema_version = 1),
        CHECK(content_id = artifact_hash)
    ) STRICT",
    "CREATE UNIQUE INDEX IF NOT EXISTS idx_artifact_content_identity_content
     ON artifact_content_identities(content_id)",
    "INSERT OR IGNORE INTO artifact_content_identities(
         artifact_hash, content_id, row_schema_version
     )
     SELECT hash, hash, 1 FROM artifacts",
    "CREATE TRIGGER IF NOT EXISTS trg_artifact_content_identity_dual_write
     AFTER INSERT ON artifacts
     BEGIN
       INSERT INTO artifact_content_identities(
           artifact_hash, content_id, row_schema_version
       ) VALUES (NEW.hash, NEW.hash, 1);
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_artifact_content_identity_immutable_update
     BEFORE UPDATE ON artifact_content_identities
     BEGIN
       SELECT RAISE(ABORT, 'artifact content identity is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_artifact_content_identity_guard_delete
     BEFORE DELETE ON artifact_content_identities
     WHEN EXISTS(
         SELECT 1 FROM artifacts WHERE hash = OLD.artifact_hash
     )
     BEGIN
       SELECT RAISE(ABORT, 'artifact content identity is retained with its artifact');
     END",
];

/// v15: typed artifact content identity for every role-qualified lineage edge.
///
/// The v1 `edges.artifact` BLOB remains the compatibility key. This companion
/// repeats it only as a typed raw-byte `ContentId`, never as a semantic ID or
/// authority. Backfill joins through the authenticated v14 artifact sidecar;
/// a trigger dual-writes every later edge, including old-handle inserts.
pub const V15: &[&str] = &[
    "CREATE TABLE IF NOT EXISTS edge_content_identities(
        op INTEGER NOT NULL,
        artifact_hash BLOB NOT NULL CHECK(length(artifact_hash) = 32),
        role TEXT NOT NULL CHECK(role IN ('in','out')),
        content_id BLOB NOT NULL CHECK(length(content_id) = 32),
        row_schema_version INTEGER NOT NULL CHECK(row_schema_version = 1),
        PRIMARY KEY(op, artifact_hash, role),
        FOREIGN KEY(op, artifact_hash, role)
            REFERENCES edges(op, artifact, role) ON DELETE CASCADE,
        FOREIGN KEY(content_id)
            REFERENCES artifact_content_identities(content_id),
        CHECK(content_id = artifact_hash)
    ) STRICT",
    "CREATE INDEX IF NOT EXISTS idx_edge_content_identity_content
     ON edge_content_identities(content_id, role, op)",
    "INSERT OR IGNORE INTO edge_content_identities(
         op, artifact_hash, role, content_id, row_schema_version
     )
     SELECT e.op, e.artifact, e.role, i.content_id, 1
     FROM edges AS e
     JOIN artifact_content_identities AS i ON i.artifact_hash = e.artifact",
    "CREATE TRIGGER IF NOT EXISTS trg_edge_content_identity_dual_write
     AFTER INSERT ON edges
     BEGIN
       INSERT INTO edge_content_identities(
           op, artifact_hash, role, content_id, row_schema_version
       ) VALUES (NEW.op, NEW.artifact, NEW.role, NEW.artifact, 1);
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_edge_content_identity_immutable_update
     BEFORE UPDATE ON edge_content_identities
     BEGIN
       SELECT RAISE(ABORT, 'edge content identity is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_edge_content_identity_guard_delete
     BEFORE DELETE ON edge_content_identities
     WHEN EXISTS(
         SELECT 1 FROM edges
         WHERE op = OLD.op AND artifact = OLD.artifact_hash AND role = OLD.role
     )
     BEGIN
       SELECT RAISE(ABORT, 'edge content identity is retained with its edge');
     END",
];

/// v16: explicit semantic bindings between retained artifacts and independently
/// verified v13 identity-migration receipts.
///
/// Multiple receipts may bind one artifact and one semantic identity may name
/// multiple artifacts. Neither storage nor lookup chooses a winner. Migration
/// creates an empty table and deliberately infers nothing from artifact kind,
/// metadata, edges, FNV values, or generic extension rows.
pub const V16: &[&str] = &[
    "CREATE TABLE IF NOT EXISTS artifact_semantic_bindings(
        artifact_hash BLOB NOT NULL CHECK(length(artifact_hash) = 32),
        receipt_id BLOB NOT NULL CHECK(length(receipt_id) = 32),
        semantic_id BLOB NOT NULL CHECK(length(semantic_id) = 32),
        identity_role INTEGER NOT NULL CHECK(identity_role BETWEEN 1 AND 12),
        identity_schema_id BLOB NOT NULL CHECK(length(identity_schema_id) = 32),
        identity_schema_version INTEGER NOT NULL
            CHECK(identity_schema_version BETWEEN 1 AND 4294967295),
        trust_state INTEGER NOT NULL CHECK(trust_state BETWEEN 0 AND 3),
        no_claim_state INTEGER NOT NULL CHECK(no_claim_state BETWEEN 0 AND 1),
        created_at INTEGER NOT NULL,
        PRIMARY KEY(artifact_hash, receipt_id),
        FOREIGN KEY(artifact_hash)
            REFERENCES artifact_content_identities(artifact_hash),
        FOREIGN KEY(receipt_id)
            REFERENCES identity_migration_receipts(receipt_id)
    ) STRICT",
    "CREATE INDEX IF NOT EXISTS idx_artifact_semantic_binding_receipt
     ON artifact_semantic_bindings(receipt_id, artifact_hash)",
    "CREATE INDEX IF NOT EXISTS idx_artifact_semantic_binding_semantic
     ON artifact_semantic_bindings(
         identity_role, identity_schema_id, identity_schema_version,
         semantic_id, artifact_hash, receipt_id
     )",
    "CREATE TRIGGER IF NOT EXISTS trg_artifact_semantic_binding_immutable_update
     BEFORE UPDATE ON artifact_semantic_bindings
     BEGIN
       SELECT RAISE(ABORT, 'artifact semantic binding is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_artifact_semantic_binding_immutable_delete
     BEFORE DELETE ON artifact_semantic_bindings
     BEGIN
       SELECT RAISE(ABORT, 'artifact semantic binding is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_artifact_semantic_binding_immutable_reinsert
     BEFORE INSERT ON artifact_semantic_bindings
     WHEN EXISTS(
         SELECT 1 FROM artifact_semantic_bindings
         WHERE artifact_hash = NEW.artifact_hash AND receipt_id = NEW.receipt_id
     )
     BEGIN
       SELECT RAISE(ABORT, 'artifact semantic binding is immutable');
     END",
];

/// v17: explicit semantic bindings for exact retained evidence JSON bytes.
///
/// Binding is opt-in and receipt-backed. Migration creates no rows because an
/// evidence name or JSON shape cannot establish a nominal schema or authority.
/// Once a binding exists, the exact evidence name/body pair is retained; exact
/// response-loss retries remain legal while reinterpretation is refused.
pub const V17: &[&str] = &[
    "CREATE TABLE IF NOT EXISTS evidence_semantic_bindings(
        evidence_name TEXT NOT NULL
            CHECK(length(CAST(evidence_name AS BLOB)) BETWEEN 1 AND 65536),
        receipt_id BLOB NOT NULL CHECK(length(receipt_id) = 32),
        content_id BLOB NOT NULL CHECK(length(content_id) = 32),
        semantic_id BLOB NOT NULL CHECK(length(semantic_id) = 32),
        identity_role INTEGER NOT NULL CHECK(identity_role BETWEEN 1 AND 12),
        identity_schema_id BLOB NOT NULL CHECK(length(identity_schema_id) = 32),
        identity_schema_version INTEGER NOT NULL
            CHECK(identity_schema_version BETWEEN 1 AND 4294967295),
        trust_state INTEGER NOT NULL CHECK(trust_state BETWEEN 0 AND 3),
        no_claim_state INTEGER NOT NULL CHECK(no_claim_state BETWEEN 0 AND 1),
        created_at INTEGER NOT NULL,
        PRIMARY KEY(evidence_name, receipt_id),
        FOREIGN KEY(evidence_name) REFERENCES evidence(name),
        FOREIGN KEY(receipt_id)
            REFERENCES identity_migration_receipts(receipt_id)
    ) STRICT",
    "CREATE INDEX IF NOT EXISTS idx_evidence_semantic_binding_content
     ON evidence_semantic_bindings(content_id, evidence_name, receipt_id)",
    "CREATE INDEX IF NOT EXISTS idx_evidence_semantic_binding_receipt
     ON evidence_semantic_bindings(receipt_id, evidence_name)",
    "CREATE INDEX IF NOT EXISTS idx_evidence_semantic_binding_semantic
     ON evidence_semantic_bindings(
         identity_role, identity_schema_id, identity_schema_version,
         semantic_id, evidence_name, receipt_id
     )",
    "CREATE TRIGGER IF NOT EXISTS trg_evidence_semantic_binding_immutable_update
     BEFORE UPDATE ON evidence_semantic_bindings
     BEGIN
       SELECT RAISE(ABORT, 'evidence semantic binding is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_evidence_semantic_binding_immutable_delete
     BEFORE DELETE ON evidence_semantic_bindings
     BEGIN
       SELECT RAISE(ABORT, 'evidence semantic binding is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_evidence_semantic_binding_immutable_reinsert
     BEFORE INSERT ON evidence_semantic_bindings
     WHEN EXISTS(
         SELECT 1 FROM evidence_semantic_bindings
         WHERE evidence_name = NEW.evidence_name AND receipt_id = NEW.receipt_id
     )
     BEGIN
       SELECT RAISE(ABORT, 'evidence semantic binding is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_evidence_semantic_binding_guard_source_update
     BEFORE UPDATE ON evidence
     WHEN EXISTS(
         SELECT 1 FROM evidence_semantic_bindings
         WHERE evidence_name = OLD.name
     ) AND (NEW.name != OLD.name OR CAST(NEW.body AS BLOB) != CAST(OLD.body AS BLOB))
     BEGIN
       SELECT RAISE(ABORT, 'bound evidence identity is immutable');
     END",
    "CREATE TRIGGER IF NOT EXISTS trg_evidence_semantic_binding_guard_source_delete
     BEFORE DELETE ON evidence
     WHEN EXISTS(
         SELECT 1 FROM evidence_semantic_bindings
         WHERE evidence_name = OLD.name
     )
     BEGIN
       SELECT RAISE(ABORT, 'bound evidence identity is retained');
     END",
];

/// Every table the CURRENT schema owns (v1 set + v2 through v17 additions); the
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
    "session_claim_discovery",
    "artifact_output_seals",
    "op_artifact_edge_seals",
    "session_checkpoint_receipts",
    "qty_dimension_crosswalks",
    "semantic_state_checkpoint_receipts",
    "identity_migration_receipts",
    "artifact_content_identities",
    "edge_content_identities",
    "artifact_semantic_bindings",
    "evidence_semantic_bindings",
];
