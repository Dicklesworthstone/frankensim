//! Mini Design Ledger on FrankenSQLite (plan §11.2 in embryo): ops +
//! content-addressed artifacts + lineage edges. Replay = re-execute the
//! recorded study and compare hashes; stored-byte corruption fails loudly.
//!
//! Hashing is FNV-1a 64 (fs-obs) — a documented placeholder for the
//! BLAKE3-class tree hash that fs-ledger-core owns.

use fsqlite::{Connection, SqliteValue};

/// Content hash rendered as fixed-width hex.
#[must_use]
pub fn content_hash(bytes: &[u8]) -> String {
    format!("{:016x}", fs_obs::fnv1a64(bytes))
}

/// A thin ledger over one fsqlite database file.
pub struct MiniLedger {
    conn: Connection,
}

impl MiniLedger {
    /// Open (creating the schema if needed).
    ///
    /// # Errors
    /// Returns a message if the database cannot be opened or migrated.
    pub fn open(path: &str) -> Result<MiniLedger, String> {
        let conn = Connection::open(path).map_err(|e| format!("ledger open {path}: {e}"))?;
        for ddl in [
            "CREATE TABLE IF NOT EXISTS artifacts(hash TEXT PRIMARY KEY, kind TEXT, bytes BLOB)",
            "CREATE TABLE IF NOT EXISTS ops(id INTEGER PRIMARY KEY, kind TEXT, ir TEXT, seed TEXT)",
            "CREATE TABLE IF NOT EXISTS edges(op INTEGER, artifact TEXT, role TEXT)",
        ] {
            conn.execute(ddl).map_err(|e| format!("ledger DDL: {e}"))?;
        }
        Ok(MiniLedger { conn })
    }

    /// Store an artifact (content-addressed; identical bytes dedupe).
    ///
    /// # Errors
    /// Returns a message on write failure.
    pub fn put_artifact(&self, kind: &str, bytes: &[u8]) -> Result<String, String> {
        let hash = content_hash(bytes);
        let existing = self
            .conn
            .query_with_params(
                "SELECT hash FROM artifacts WHERE hash = ?1",
                &[SqliteValue::Text(hash.clone().into())],
            )
            .map_err(|e| format!("artifact lookup: {e}"))?;
        if existing.is_empty() {
            self.conn
                .prepare("INSERT INTO artifacts(hash, kind, bytes) VALUES (?1, ?2, ?3)")
                .map_err(|e| format!("artifact insert prepare: {e}"))?
                .execute_with_params(&[
                    SqliteValue::Text(hash.clone().into()),
                    SqliteValue::Text(kind.into()),
                    SqliteValue::Blob(bytes.to_vec().into()),
                ])
                .map_err(|e| format!("artifact insert: {e}"))?;
        }
        Ok(hash)
    }

    /// Record an executed op with its (frozen) IR and seed.
    ///
    /// # Errors
    /// Returns a message on write failure.
    pub fn record_op(&self, kind: &str, ir: &str, seed_hex: &str) -> Result<i64, String> {
        self.conn
            .prepare("INSERT INTO ops(kind, ir, seed) VALUES (?1, ?2, ?3)")
            .map_err(|e| format!("op prepare: {e}"))?
            .execute_with_params(&[
                SqliteValue::Text(kind.into()),
                SqliteValue::Text(ir.into()),
                SqliteValue::Text(seed_hex.into()),
            ])
            .map_err(|e| format!("op insert: {e}"))?;
        let row = self
            .conn
            .query_row("SELECT MAX(id) FROM ops")
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
        self.conn
            .prepare("INSERT INTO edges(op, artifact, role) VALUES (?1, ?2, ?3)")
            .map_err(|e| format!("edge prepare: {e}"))?
            .execute_with_params(&[
                SqliteValue::Integer(op),
                SqliteValue::Text(artifact.into()),
                SqliteValue::Text(role.into()),
            ])
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
            .query("SELECT bytes FROM artifacts WHERE kind = 'study-ir'")
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
            .query("SELECT hash, bytes FROM artifacts")
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
            .query("SELECT hash FROM artifacts WHERE kind != 'study-ir'")
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
            .execute("UPDATE artifacts SET bytes = X'DEADBEEF' WHERE kind != 'study-ir'")
            .map_err(|e| format!("corruption hook: {e}"))?;
        Ok(())
    }
}
