//! Mini Design Ledger on FrankenSQLite (plan §11.2 in embryo): ops +
//! content-addressed artifacts + lineage edges. Replay = re-execute the
//! recorded study and compare hashes; stored-byte corruption fails loudly.
//!
//! Hashing is domain-separated BLAKE3 (bead frankensim-ynsl): this crate
//! patterns the flagship runner binaries and sits inside the flywheel
//! loop, where content-address integrity is the merge/skip soundness
//! basis — a 64-bit FNV placeholder was collision-cheap there. Ledger
//! files carry a FORMAT VERSION; pre-BLAKE3 (v1/FNV) files are
//! version-refused with a teaching error, never silently misread.

use fsqlite::{AsyncConnection, SqliteValue};
use std::marker::PhantomData;
use std::rc::Rc;

/// The ledger format this crate writes and reads. v1 was the FNV-1a era
/// (16-hex hashes, no meta table); v2 is domain-separated BLAKE3.
pub const LEDGER_FORMAT_VERSION: &str = "2";

/// Semantic version of the mini-ledger artifact content address.
pub const ARTIFACT_CONTENT_IDENTITY_VERSION: u32 = 2;

/// Domain separating mini-ledger artifact bytes from every other BLAKE3 use.
pub const ARTIFACT_CONTENT_IDENTITY_DOMAIN: &str = "frankensim.fs-vskeleton.artifact.v2";

/// Owner-local declaration consumed by `xtask check-identities`.
pub const ARTIFACT_CONTENT_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-vskeleton:artifact-content",
    "version_const=ARTIFACT_CONTENT_IDENTITY_VERSION",
    "version=2",
    "domain=frankensim.fs-vskeleton.artifact.v2",
    "domain_const=ARTIFACT_CONTENT_IDENTITY_DOMAIN",
    "encoder=content_hash",
    "encoder_helpers=content_hash_with_domain",
    "schema_functions=artifact_content_identity_version_is_supported,MiniLedger::open,MiniLedger::put_artifact,MiniLedger::verify_artifact_integrity,crates/fs-blake3/src/lib.rs#ContentHash::to_hex,crates/fs-blake3/src/lib.rs#hash_domain",
    "schema_constants=ARTIFACT_CONTENT_IDENTITY_VERSION,ARTIFACT_CONTENT_IDENTITY_DOMAIN,LEDGER_FORMAT_VERSION",
    "schema_dependencies=none",
    "digest=fs-blake3",
    "encoding=canonical-transport-exact-bits",
    "sources=ArtifactContentIdentityInput",
    "source_fields=ArtifactContentIdentityInput.bytes:semantic",
    "source_bindings=ArtifactContentIdentityInput.bytes>artifact-bytes",
    "external_semantic_fields=artifact-domain",
    "semantic_fields=artifact-domain,artifact-bytes",
    "excluded_fields=artifact-kind:metadata-only-not-content-address",
    "consumers=MiniLedger::put_artifact,MiniLedger::verify_artifact_integrity,crate::run_study,crate::replay",
    "mutations=artifact-domain:crates/fs-vskeleton/src/ledger.rs#artifact_content_domain_moves_identity,artifact-bytes:crates/fs-vskeleton/src/ledger.rs#artifact_content_bytes_move_identity",
    "nonsemantic_mutations=artifact-kind:crates/fs-vskeleton/src/ledger.rs#artifact_kind_does_not_move_content_identity",
    "field_guard=classify_artifact_content_identity_fields",
    "transport_guard=content_hash",
    "version_guard=crates/fs-vskeleton/src/ledger.rs#artifact_content_identity_version_fails_closed",
    "coupling_surface=fs-vskeleton:artifact-content",
];

struct ArtifactContentIdentityInput<'a> {
    bytes: &'a [u8],
}

/// Whether a retained mini-ledger content address uses the only identity
/// semantics accepted by this build.
#[must_use]
pub const fn artifact_content_identity_version_is_supported(declared: u32) -> bool {
    declared == ARTIFACT_CONTENT_IDENTITY_VERSION
}

// The hash domain MUST carry the format version, so a future version bump that
// forgets to re-tag the domain cannot let two ledger formats hash identical
// bytes to the SAME content address (the cross-version replay this format gate
// exists to prevent). The `.v2` above and `LEDGER_FORMAT_VERSION` were two
// independent literals; pin them together at COMPILE TIME — the domain must end
// with `v{LEDGER_FORMAT_VERSION}`, so bumping the version without updating the
// domain fails the build rather than shipping a silent cross-version collision.
const _: () = {
    let d = ARTIFACT_CONTENT_IDENTITY_DOMAIN.as_bytes();
    let v = LEDGER_FORMAT_VERSION.as_bytes();
    assert!(
        d.len() > v.len(),
        "hash domain must be longer than the version tag"
    );
    assert!(
        d[d.len() - v.len() - 1] == b'v',
        "hash domain must end with 'v' + LEDGER_FORMAT_VERSION"
    );
    let mut i = 0;
    while i < v.len() {
        assert!(
            d[d.len() - v.len() + i] == v[i],
            "hash domain version tag must equal LEDGER_FORMAT_VERSION"
        );
        i += 1;
    }
};

/// Content hash rendered as fixed-width hex (64 chars, BLAKE3).
#[must_use]
pub fn content_hash(bytes: &[u8]) -> String {
    content_hash_with_domain(
        ARTIFACT_CONTENT_IDENTITY_DOMAIN,
        &ArtifactContentIdentityInput { bytes },
    )
}

fn content_hash_with_domain(domain: &str, input: &ArtifactContentIdentityInput<'_>) -> String {
    fs_blake3::hash_domain(domain, input.bytes).to_hex()
}

#[allow(dead_code)] // exhaustive source-shape guard consumed by xtask
fn classify_artifact_content_identity_fields(input: &ArtifactContentIdentityInput<'_>) {
    let ArtifactContentIdentityInput { bytes: _ } = input;
}

/// A thin ledger over one fsqlite database file.
pub struct MiniLedger {
    conn: AsyncConnection,
    _single_thread: PhantomData<Rc<()>>,
}

impl MiniLedger {
    /// Open (creating the schema if needed).
    ///
    /// # Errors
    /// Returns a message if the database cannot be opened or migrated.
    pub fn open(path: &str) -> Result<MiniLedger, String> {
        let conn =
            AsyncConnection::open_sync(path).map_err(|e| format!("ledger open {path}: {e}"))?;
        for ddl in [
            "CREATE TABLE IF NOT EXISTS artifacts(hash TEXT PRIMARY KEY, kind TEXT, bytes BLOB)",
            "CREATE TABLE IF NOT EXISTS ops(id INTEGER PRIMARY KEY, kind TEXT, ir TEXT, seed TEXT)",
            "CREATE TABLE IF NOT EXISTS edges(op INTEGER, artifact TEXT, role TEXT)",
            "CREATE TABLE IF NOT EXISTS vskeleton_meta(key TEXT PRIMARY KEY, value TEXT)",
        ] {
            conn.execute_sync(ddl)
                .map_err(|e| format!("ledger DDL: {e}"))?;
        }
        // FORMAT ATTESTATION (the gp3.18 doctrine at embryo scale): a
        // version row is stamped on first use of an EMPTY ledger; a file
        // with artifacts but no (or a different) version is another
        // format's data and is refused with the migration named — never
        // silently misread under a new hash function.
        let version = conn
            .query_sync("SELECT value FROM vskeleton_meta WHERE key = 'format_version'")
            .map_err(|e| format!("ledger version read: {e}"))?;
        match version.first().and_then(|row| row.get(0)) {
            Some(SqliteValue::Text(v)) if v.as_str() == LEDGER_FORMAT_VERSION => {}
            Some(SqliteValue::Text(v)) => {
                return Err(format!(
                    "LedgerFormatMismatch: {path} is format v{} but this build reads/writes \
                     v{LEDGER_FORMAT_VERSION}; replay it with a matching build or re-run the \
                     study into a fresh ledger — hashes are not comparable across formats",
                    v.as_str()
                ));
            }
            _ => {
                let artifacts = conn
                    .query_sync("SELECT hash FROM artifacts LIMIT 1")
                    .map_err(|e| format!("ledger census: {e}"))?;
                if !artifacts.is_empty() {
                    return Err(format!(
                        "LedgerFormatMismatch: {path} holds artifacts but no format version — \
                         a pre-v2 (FNV-era) ledger; its 16-hex hashes are not comparable to \
                         v{LEDGER_FORMAT_VERSION} BLAKE3 addresses; replay the original study \
                         into a fresh ledger instead of migrating hashes in place"
                    ));
                }
                let sql = "INSERT INTO vskeleton_meta(key, value) VALUES ('format_version', ?1)";
                conn.prepare_sync(sql)
                    .map_err(|e| format!("ledger version stamp prepare: {e}"))?;
                conn.execute_with_params_sync(
                    sql,
                    &[SqliteValue::Text(LEDGER_FORMAT_VERSION.into())],
                )
                .map_err(|e| format!("ledger version stamp: {e}"))?;
            }
        }
        Ok(MiniLedger {
            conn,
            _single_thread: PhantomData,
        })
    }

    /// Store an artifact (content-addressed; identical bytes dedupe).
    ///
    /// # Errors
    /// Returns a message on write failure.
    pub fn put_artifact(&self, kind: &str, bytes: &[u8]) -> Result<String, String> {
        let hash = content_hash(bytes);
        let existing = self
            .conn
            .query_with_params_sync(
                "SELECT hash FROM artifacts WHERE hash = ?1",
                &[SqliteValue::Text(hash.clone().into())],
            )
            .map_err(|e| format!("artifact lookup: {e}"))?;
        if existing.is_empty() {
            let sql = "INSERT INTO artifacts(hash, kind, bytes) VALUES (?1, ?2, ?3)";
            self.conn
                .prepare_sync(sql)
                .map_err(|e| format!("artifact insert prepare: {e}"))?;
            self.conn
                .execute_with_params_sync(
                    sql,
                    &[
                        SqliteValue::Text(hash.clone().into()),
                        SqliteValue::Text(kind.into()),
                        SqliteValue::Blob(bytes.to_vec().into()),
                    ],
                )
                .map_err(|e| format!("artifact insert: {e}"))?;
        }
        Ok(hash)
    }

    /// Record an executed op with its (frozen) IR and seed.
    ///
    /// # Errors
    /// Returns a message on write failure.
    pub fn record_op(&self, kind: &str, ir: &str, seed_hex: &str) -> Result<i64, String> {
        let sql = "INSERT INTO ops(kind, ir, seed) VALUES (?1, ?2, ?3)";
        self.conn
            .prepare_sync(sql)
            .map_err(|e| format!("op prepare: {e}"))?;
        self.conn
            .execute_with_params_sync(
                sql,
                &[
                    SqliteValue::Text(kind.into()),
                    SqliteValue::Text(ir.into()),
                    SqliteValue::Text(seed_hex.into()),
                ],
            )
            .map_err(|e| format!("op insert: {e}"))?;
        let row = self
            .conn
            .query_row_sync("SELECT MAX(id) FROM ops")
            .map_err(|e| format!("op id: {e}"))?;
        match row.get(0) {
            Some(SqliteValue::Integer(id)) => Ok(*id),
            other => Err(format!("op id: unexpected value {other:?}")),
        }
    }

    /// Link an op to an artifact with a role ("in"/"out").
    ///
    /// # Errors
    /// Returns a message on write failure.
    pub fn link(&self, op: i64, artifact: &str, role: &str) -> Result<(), String> {
        let sql = "INSERT INTO edges(op, artifact, role) VALUES (?1, ?2, ?3)";
        self.conn
            .prepare_sync(sql)
            .map_err(|e| format!("edge prepare: {e}"))?;
        self.conn
            .execute_with_params_sync(
                sql,
                &[
                    SqliteValue::Integer(op),
                    SqliteValue::Text(artifact.into()),
                    SqliteValue::Text(role.into()),
                ],
            )
            .map_err(|e| format!("edge insert: {e}"))?;
        Ok(())
    }

    /// The recorded study IR (the first `study-ir` artifact).
    ///
    /// # Errors
    /// Errors if the ledger holds no study.
    pub fn get_study_ir(&self) -> Result<String, String> {
        let rows = self
            .conn
            .query_sync("SELECT bytes FROM artifacts WHERE kind = 'study-ir'")
            .map_err(|e| format!("study lookup: {e}"))?;
        let row = rows
            .first()
            .ok_or("ledger holds no study-ir artifact — nothing to replay")?;
        match row.get(0) {
            Some(SqliteValue::Blob(b)) => String::from_utf8(b.to_vec())
                .map_err(|e| format!("study-ir is not UTF-8 (corruption?): {e}")),
            other => Err(format!("study-ir: unexpected value {other:?}")),
        }
    }

    /// Every stored artifact's bytes must still hash to its recorded key —
    /// byte-level corruption fails LOUDLY here.
    ///
    /// # Errors
    /// Names the first corrupted artifact.
    pub fn verify_artifact_integrity(&self) -> Result<(), String> {
        let rows = self
            .conn
            .query_sync("SELECT hash, bytes FROM artifacts")
            .map_err(|e| format!("integrity scan: {e}"))?;
        for row in &rows {
            let (Some(SqliteValue::Text(h)), Some(SqliteValue::Blob(b))) = (row.get(0), row.get(1))
            else {
                return Err("integrity scan: malformed artifact row".to_string());
            };
            let actual = content_hash(b);
            if actual != h.as_str() {
                return Err(format!(
                    "LedgerCorruption: artifact recorded as {h} hashes to {actual} — bytes were \
                     modified after recording; refuse to replay from a tampered ledger"
                ));
            }
        }
        Ok(())
    }

    /// Every non-study artifact hash (row order is UNSPECIFIED — the query
    /// carries no `ORDER BY`; callers comparing against a recomputed list must
    /// treat the result as a multiset and sort, as `replay` does).
    ///
    /// # Errors
    /// Returns a message on read failure.
    pub fn artifact_hashes_excluding_study(&self) -> Result<Vec<String>, String> {
        let rows = self
            .conn
            .query_sync("SELECT hash FROM artifacts WHERE kind != 'study-ir'")
            .map_err(|e| format!("hash scan: {e}"))?;
        let mut out = Vec::with_capacity(rows.len());
        for row in &rows {
            match row.get(0) {
                Some(SqliteValue::Text(h)) => out.push(h.as_str().to_string()),
                other => return Err(format!("hash scan: unexpected value {other:?}")),
            }
        }
        Ok(out)
    }

    /// Deliberately corrupt one artifact's bytes (test hook for the
    /// corruption-fails-loudly exit criterion).
    ///
    /// # Errors
    /// Returns a message on write failure.
    pub fn corrupt_first_artifact_for_test(&self) -> Result<(), String> {
        self.conn
            .execute_sync("UPDATE artifacts SET bytes = X'DEADBEEF' WHERE kind != 'study-ir'")
            .map_err(|e| format!("corruption hook: {e}"))?;
        Ok(())
    }
}

#[cfg(test)]
mod identity_tests {
    use super::*;

    #[test]
    fn artifact_content_domain_moves_identity() {
        let input = ArtifactContentIdentityInput { bytes: b"payload" };
        assert_ne!(
            content_hash_with_domain(ARTIFACT_CONTENT_IDENTITY_DOMAIN, &input),
            content_hash_with_domain("frankensim.fs-vskeleton.artifact.v2.alternate", &input)
        );
    }

    #[test]
    fn artifact_content_bytes_move_identity() {
        assert_ne!(content_hash(b"payload-a"), content_hash(b"payload-b"));
    }

    #[test]
    fn artifact_kind_does_not_move_content_identity() {
        let identity_for_kind = |_kind: &str| content_hash(b"shared payload");
        assert_eq!(identity_for_kind("study-ir"), identity_for_kind("field"));
    }

    #[test]
    fn artifact_content_identity_version_fails_closed() {
        assert_eq!(ARTIFACT_CONTENT_IDENTITY_VERSION, 2);
        assert_eq!(LEDGER_FORMAT_VERSION, "2");
        assert!(ARTIFACT_CONTENT_IDENTITY_DOMAIN.ends_with(".v2"));
        assert!(artifact_content_identity_version_is_supported(2));
        assert!(!artifact_content_identity_version_is_supported(1));
        assert!(!artifact_content_identity_version_is_supported(3));
    }
}
