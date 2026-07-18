//! Immutable persistence for versioned semantic crosswalk receipts.
//!
//! The first supported family is fs-qty's exact-byte five-to-six dimension
//! migration. A row is admitted only after both content-addressed artifacts
//! are re-read under a fixed budget and the producer receipt is independently
//! replayed. Lookup repeats that verification; row presence alone is never
//! evidence that a migration remains reproducible.

use std::str;

use fs_qty::json::{
    DimensionCrosswalkReceipt, FiveToSixRule, QtyWireVersion, decode_json, to_json,
};
use fsqlite::SqliteValue;

use crate::{ContentHash, Ledger, LedgerError, blob_param, now_wall_ns, sql_err, text_param};

/// Maximum exact source or target JSON bytes retained for one quantity
/// crosswalk. Canonical fs-qty values are tiny; this explicit envelope keeps a
/// hostile artifact from forcing an unbounded allocation during admission or
/// replay.
pub const MAX_QTY_CROSSWALK_JSON_BYTES: u64 = 4 * 1024;

const APPEND_MOLE_ZERO_RULE: &str = "append-mole-zero";

/// Result of recording one immutable five-to-six quantity crosswalk.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct QtyDimensionCrosswalkWrite {
    old_hash: ContentHash,
    new_hash: ContentHash,
    deduped: bool,
}

impl QtyDimensionCrosswalkWrite {
    /// Exact historical JSON content identity.
    #[must_use]
    pub const fn old_hash(self) -> ContentHash {
        self.old_hash
    }

    /// Exact canonical six-base JSON content identity.
    #[must_use]
    pub const fn new_hash(self) -> ContentHash {
        self.new_hash
    }

    /// Whether an identical verified row already existed.
    #[must_use]
    pub const fn deduped(self) -> bool {
        self.deduped
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct StoredQtyDimensionCrosswalk {
    old_hash: ContentHash,
    new_hash: ContentHash,
}

fn invalid(field: &str, problem: impl Into<String>) -> LedgerError {
    LedgerError::Invalid {
        field: field.to_string(),
        problem: problem.into(),
    }
}

fn stored_corrupt(old_hash: ContentHash, detail: impl Into<String>) -> LedgerError {
    LedgerError::Corrupt {
        hash_hex: old_hash.to_hex(),
        detail: format!(
            "quantity dimension crosswalk row is corrupt: {}",
            detail.into()
        ),
    }
}

fn replay_refusal(
    row: StoredQtyDimensionCrosswalk,
    stored_row: bool,
    field: &'static str,
    detail: impl Into<String>,
) -> LedgerError {
    let detail = detail.into();
    if stored_row {
        stored_corrupt(row.old_hash, format!("{field}: {detail}"))
    } else {
        invalid(field, detail)
    }
}

fn hash_from_value(
    value: Option<&SqliteValue>,
    old_hash: ContentHash,
    field: &'static str,
) -> Result<ContentHash, LedgerError> {
    let Some(SqliteValue::Blob(bytes)) = value else {
        return Err(stored_corrupt(old_hash, format!("{field} is not a BLOB")));
    };
    ContentHash::from_slice(bytes).ok_or_else(|| {
        stored_corrupt(
            old_hash,
            format!(
                "{field} must contain exactly 32 bytes, found {}",
                bytes.len()
            ),
        )
    })
}

impl Ledger {
    fn stored_qty_dimension_crosswalk(
        &self,
        old_hash: ContentHash,
    ) -> Result<Option<StoredQtyDimensionCrosswalk>, LedgerError> {
        let preflight = self
            .conn
            .query_with_params(
                "SELECT typeof(new_hash), length(new_hash), \
                        typeof(source_version), typeof(target_version), \
                        typeof(rule), length(rule), typeof(created_at) \
                 FROM qty_dimension_crosswalks WHERE old_hash = ?1 LIMIT 2",
                &[blob_param(old_hash.as_bytes())],
            )
            .map_err(|error| sql_err("quantity dimension crosswalk preflight", &error))?;
        if preflight.is_empty() {
            return Ok(None);
        }
        if preflight.len() != 1 {
            return Err(stored_corrupt(
                old_hash,
                "one source hash names multiple rows",
            ));
        }
        let Some(metadata) = preflight.first() else {
            return Err(stored_corrupt(
                old_hash,
                "row disappeared before bounded metadata validation",
            ));
        };
        let type_is = |index, expected: &str| matches!(metadata.get(index), Some(SqliteValue::Text(value)) if value.as_str() == expected);
        let integer_is = |index, expected: i64| matches!(metadata.get(index), Some(SqliteValue::Integer(value)) if *value == expected);
        if !type_is(0, "blob")
            || !integer_is(1, 32)
            || !type_is(2, "integer")
            || !type_is(3, "integer")
            || !type_is(4, "text")
            || !integer_is(
                5,
                i64::try_from(APPEND_MOLE_ZERO_RULE.len()).unwrap_or(i64::MAX),
            )
            || !type_is(6, "integer")
        {
            return Err(stored_corrupt(
                old_hash,
                "stored types or bounded lengths violate the v11 row envelope",
            ));
        }

        // Repeat every type, length, version, and rule predicate in the
        // payload query. A row changed after preflight therefore disappears
        // instead of materializing hostile variable-size storage.
        let rows = self
            .conn
            .query_with_params(
                "SELECT old_hash, \
                        CASE WHEN typeof(new_hash) = 'blob' AND length(new_hash) = 32 \
                                   AND typeof(source_version) = 'integer' \
                                   AND source_version = 1 \
                                   AND typeof(target_version) = 'integer' \
                                   AND target_version = 2 \
                                   AND typeof(rule) = 'text' \
                                   AND length(rule) = ?2 AND rule = ?3 \
                                   AND typeof(created_at) = 'integer' \
                             THEN new_hash ELSE NULL END \
                 FROM qty_dimension_crosswalks \
                 WHERE old_hash = ?1 \
                   AND typeof(old_hash) = 'blob' AND length(old_hash) = 32 \
                 LIMIT 2",
                &[
                    blob_param(old_hash.as_bytes()),
                    SqliteValue::Integer(
                        i64::try_from(APPEND_MOLE_ZERO_RULE.len()).unwrap_or(i64::MAX),
                    ),
                    text_param(APPEND_MOLE_ZERO_RULE),
                ],
            )
            .map_err(|error| sql_err("quantity dimension crosswalk guarded read", &error))?;
        if rows.len() != 1 {
            return Err(stored_corrupt(
                old_hash,
                "row disappeared after bounded metadata preflight",
            ));
        }
        let Some(row) = rows.first() else {
            return Err(stored_corrupt(
                old_hash,
                "row disappeared after bounded metadata preflight",
            ));
        };
        let stored_old = hash_from_value(row.get(0), old_hash, "old_hash")?;
        if stored_old != old_hash {
            return Err(stored_corrupt(
                old_hash,
                "indexed old_hash disagrees with the selected row",
            ));
        }
        let new_hash = hash_from_value(row.get(1), old_hash, "new_hash")?;
        if new_hash == old_hash {
            return Err(stored_corrupt(
                old_hash,
                "old_hash and new_hash must remain distinct",
            ));
        }
        Ok(Some(StoredQtyDimensionCrosswalk { old_hash, new_hash }))
    }

    fn replay_qty_dimension_crosswalk(
        &self,
        row: StoredQtyDimensionCrosswalk,
        stored_row: bool,
    ) -> Result<DimensionCrosswalkReceipt, LedgerError> {
        let old_bytes = self
            .get_artifact_bounded(&row.old_hash, MAX_QTY_CROSSWALK_JSON_BYTES)?
            .ok_or_else(|| {
                replay_refusal(
                    row,
                    stored_row,
                    "qty_dimension_crosswalk.old_hash",
                    format!("source artifact {} does not exist", row.old_hash),
                )
            })?;
        let new_bytes = self
            .get_artifact_bounded(&row.new_hash, MAX_QTY_CROSSWALK_JSON_BYTES)?
            .ok_or_else(|| {
                replay_refusal(
                    row,
                    stored_row,
                    "qty_dimension_crosswalk.new_hash",
                    format!("target artifact {} does not exist", row.new_hash),
                )
            })?;
        let old_json = str::from_utf8(&old_bytes).map_err(|error| {
            replay_refusal(
                row,
                stored_row,
                "qty_dimension_crosswalk.old_hash",
                format!(
                    "source artifact {} is not UTF-8 JSON: {error}",
                    row.old_hash
                ),
            )
        })?;
        let decoded = decode_json(old_json).map_err(|error| {
            replay_refusal(
                row,
                stored_row,
                "qty_dimension_crosswalk.old_hash",
                format!(
                    "source artifact {} is not canonical legacy quantity JSON: {error}",
                    row.old_hash
                ),
            )
        })?;
        let receipt = decoded.migration().cloned().ok_or_else(|| {
            replay_refusal(
                row,
                stored_row,
                "qty_dimension_crosswalk.source_version",
                format!(
                    "source artifact {} decoded without the required five-to-six receipt",
                    row.old_hash
                ),
            )
        })?;
        let canonical = to_json(decoded.qty()).map_err(|error| {
            replay_refusal(
                row,
                stored_row,
                "qty_dimension_crosswalk.new_hash",
                format!(
                    "source artifact {} could not reproduce canonical six-base JSON: {error}",
                    row.old_hash
                ),
            )
        })?;
        if canonical.as_bytes() != new_bytes {
            return Err(replay_refusal(
                row,
                stored_row,
                "qty_dimension_crosswalk.new_hash",
                format!(
                    "target artifact {} is not the exact canonical six-base output for source {}",
                    row.new_hash, row.old_hash
                ),
            ));
        }
        if receipt.old_hash() != row.old_hash || receipt.new_hash() != row.new_hash {
            return Err(stored_corrupt(
                row.old_hash,
                "replayed receipt hashes disagree with the immutable row",
            ));
        }
        if !receipt.verifies(&old_bytes, &new_bytes) {
            return Err(stored_corrupt(
                row.old_hash,
                "fs-qty refused the retained source-to-target byte mapping",
            ));
        }
        Ok(receipt)
    }

    /// Persist or exactly replay one fs-qty five-to-six crosswalk.
    ///
    /// Both source and target artifacts must already exist. The call owns one
    /// transaction, re-reads both artifacts under a 4 KiB-per-side budget,
    /// reproduces the canonical target with fs-qty, and compares every
    /// producer-receipt field before inserting. An identical retry reports
    /// `deduped`; a changed target for the same historical source refuses.
    ///
    /// # Errors
    /// Open caller transactions, unsupported receipt versions/rules, missing,
    /// oversized, corrupt, or non-canonical artifacts, conflicting rows, and
    /// database failures.
    pub fn record_qty_dimension_crosswalk(
        &self,
        receipt: &DimensionCrosswalkReceipt,
    ) -> Result<QtyDimensionCrosswalkWrite, LedgerError> {
        if self.in_transaction() {
            return Err(invalid(
                "qty_dimension_crosswalk.transaction",
                "crosswalk recording must own its transaction; commit or roll back first",
            ));
        }
        if receipt.source_version() != QtyWireVersion::LegacyFive
            || receipt.target_version() != QtyWireVersion::SixBase
            || receipt.rule() != FiveToSixRule::AppendMoleZero
        {
            return Err(invalid(
                "qty_dimension_crosswalk.receipt",
                "only the exact legacy-five to six-base append-mole-zero receipt is supported",
            ));
        }
        let offered = StoredQtyDimensionCrosswalk {
            old_hash: receipt.old_hash(),
            new_hash: receipt.new_hash(),
        };

        self.begin()?;
        let write = (|| {
            let replayed = self.replay_qty_dimension_crosswalk(offered, false)?;
            if &replayed != receipt {
                return Err(invalid(
                    "qty_dimension_crosswalk.receipt",
                    "offered producer receipt differs from independent replay",
                ));
            }
            match self.stored_qty_dimension_crosswalk(offered.old_hash)? {
                Some(stored) if stored == offered => {
                    self.replay_qty_dimension_crosswalk(stored, true)?;
                    Ok(QtyDimensionCrosswalkWrite {
                        old_hash: offered.old_hash,
                        new_hash: offered.new_hash,
                        deduped: true,
                    })
                }
                Some(stored) => Err(invalid(
                    "qty_dimension_crosswalk.old_hash",
                    format!(
                        "source {} is already immutably mapped to {}, not {}",
                        offered.old_hash, stored.new_hash, offered.new_hash
                    ),
                )),
                None => {
                    self.conn
                        .prepare(
                            "INSERT INTO qty_dimension_crosswalks(\
                                old_hash, new_hash, source_version, target_version, rule, created_at\
                             ) VALUES (?1, ?2, ?3, ?4, ?5, ?6)",
                        )
                        .map_err(|error| {
                            sql_err("quantity dimension crosswalk insert prepare", &error)
                        })?
                        .execute_with_params(&[
                            blob_param(offered.old_hash.as_bytes()),
                            blob_param(offered.new_hash.as_bytes()),
                            SqliteValue::Integer(i64::from(fs_qty::json::LEGACY_WIRE_VERSION)),
                            SqliteValue::Integer(i64::from(fs_qty::json::WIRE_VERSION)),
                            text_param(APPEND_MOLE_ZERO_RULE),
                            SqliteValue::Integer(now_wall_ns()),
                        ])
                        .map_err(|error| {
                            sql_err("quantity dimension crosswalk insert", &error)
                        })?;
                    Ok(QtyDimensionCrosswalkWrite {
                        old_hash: offered.old_hash,
                        new_hash: offered.new_hash,
                        deduped: false,
                    })
                }
            }
        })();
        match write {
            Ok(write) => {
                if let Err(error) = self.commit() {
                    let _ = self.rollback();
                    return Err(error);
                }
                Ok(write)
            }
            Err(error) => {
                let _ = self.rollback();
                Err(error)
            }
        }
    }

    /// Recover and independently re-verify one persisted five-to-six
    /// crosswalk by its exact historical content hash.
    ///
    /// Absence is `Ok(None)`. A present row returns only after both retained
    /// artifacts re-hash correctly and fs-qty reproduces the same typed
    /// receipt and exact canonical target bytes.
    ///
    /// # Errors
    /// Malformed stored rows, missing, oversized, corrupt, or non-canonical
    /// artifacts, and database failures.
    pub fn qty_dimension_crosswalk(
        &self,
        old_hash: ContentHash,
    ) -> Result<Option<DimensionCrosswalkReceipt>, LedgerError> {
        let Some(row) = self.stored_qty_dimension_crosswalk(old_hash)? else {
            return Ok(None);
        };
        self.replay_qty_dimension_crosswalk(row, true).map(Some)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{SCHEMA_VERSION, schema};

    fn drop_v11_objects(ledger: &Ledger) {
        for ddl in [
            "DROP TRIGGER IF EXISTS trg_qty_dimension_crosswalks_immutable_reinsert",
            "DROP TRIGGER IF EXISTS trg_qty_dimension_crosswalks_immutable_delete",
            "DROP TRIGGER IF EXISTS trg_qty_dimension_crosswalks_immutable_update",
            "DROP INDEX IF EXISTS idx_qty_dimension_crosswalks_new_hash",
            "DROP TABLE IF EXISTS qty_dimension_crosswalks",
        ] {
            ledger.conn.execute(ddl).expect("remove v11 fixture object");
        }
    }

    fn retained_receipt(ledger: &Ledger) -> DimensionCrosswalkReceipt {
        const OLD: &str = r#"{"value":0.25,"dims":[-1,1,-1,0,0]}"#;
        let decoded = decode_json(OLD).expect("canonical legacy fixture");
        let receipt = decoded.migration().cloned().expect("migration receipt");
        let new = to_json(decoded.qty()).expect("canonical target fixture");
        ledger
            .put_artifact("quantity-json-v1", OLD.as_bytes(), None)
            .expect("retain source fixture");
        ledger
            .put_artifact("quantity-json-v2", new.as_bytes(), None)
            .expect("retain target fixture");
        receipt
    }

    #[test]
    fn genuine_v10_and_stale_v11_markers_migrate_to_v11() {
        let ledger = Ledger::open(":memory:").expect("fresh v11 ledger");
        drop_v11_objects(&ledger);
        ledger
            .conn
            .execute("PRAGMA user_version = 10")
            .expect("mark genuine v10 fixture");
        ledger
            .migrate_from_observed_version(10)
            .expect("migrate genuine v10 fixture");
        assert_eq!(ledger.schema_version().unwrap(), SCHEMA_VERSION);
        assert_eq!(ledger.table_count("qty_dimension_crosswalks").unwrap(), 0);

        ledger
            .conn
            .execute("PRAGMA user_version = 10")
            .expect("install stale v10 marker over exact v11 objects");
        ledger
            .migrate_from_observed_version(10)
            .expect("heal exact pre-applied v11 objects");
        assert_eq!(ledger.schema_version().unwrap(), SCHEMA_VERSION);
        assert_eq!(ledger.table_count("qty_dimension_crosswalks").unwrap(), 0);
    }

    #[test]
    fn divergent_early_v11_object_refuses_before_marker_advances() {
        let ledger = Ledger::open(":memory:").expect("fresh v11 ledger");
        drop_v11_objects(&ledger);
        ledger
            .conn
            .execute("CREATE TABLE qty_dimension_crosswalks(alien INTEGER) STRICT")
            .expect("install divergent early object");
        ledger
            .conn
            .execute("PRAGMA user_version = 10")
            .expect("mark v10 fixture");
        assert!(matches!(
            ledger.migrate_from_observed_version(10),
            Err(LedgerError::SchemaMismatch {
                claimed_version: 10,
                ..
            })
        ));
        assert_eq!(ledger.schema_version().unwrap(), 10);
    }

    #[test]
    fn v11_constraints_and_immutability_guards_refuse_raw_bypasses() {
        let ledger = Ledger::open(":memory:").expect("fresh v11 ledger");
        let old = ledger
            .put_artifact("fixture", b"old", None)
            .expect("old artifact")
            .hash;
        let new = ledger
            .put_artifact("fixture", b"new", None)
            .expect("new artifact")
            .hash;
        let insert = |source_version: i64,
                      target_version: i64,
                      rule: &str,
                      old_hash: ContentHash,
                      new_hash: ContentHash| {
            ledger
                .conn
                .prepare(
                    "INSERT INTO qty_dimension_crosswalks(\
                        old_hash, new_hash, source_version, target_version, rule, created_at\
                     ) VALUES (?1, ?2, ?3, ?4, ?5, 0)",
                )
                .expect("prepare raw constraint fixture")
                .execute_with_params(&[
                    blob_param(old_hash.as_bytes()),
                    blob_param(new_hash.as_bytes()),
                    SqliteValue::Integer(source_version),
                    SqliteValue::Integer(target_version),
                    text_param(rule),
                ])
        };
        assert!(insert(9, 2, APPEND_MOLE_ZERO_RULE, old, new).is_err());
        assert!(insert(1, 9, APPEND_MOLE_ZERO_RULE, old, new).is_err());
        assert!(insert(1, 2, "foreign-rule", old, new).is_err());
        assert!(insert(1, 2, APPEND_MOLE_ZERO_RULE, old, old).is_err());
        assert!(
            insert(
                1,
                2,
                APPEND_MOLE_ZERO_RULE,
                ContentHash([0x71; 32]),
                ContentHash([0x72; 32]),
            )
            .is_err()
        );
        assert_eq!(ledger.table_count("qty_dimension_crosswalks").unwrap(), 0);

        let receipt = retained_receipt(&ledger);
        ledger
            .record_qty_dimension_crosswalk(&receipt)
            .expect("public verified insert");
        let update = ledger
            .conn
            .execute("UPDATE qty_dimension_crosswalks SET created_at = created_at + 1")
            .expect_err("immutable update trigger");
        assert!(
            update
                .to_string()
                .contains("quantity dimension crosswalk is immutable")
        );
        let delete = ledger
            .conn
            .execute("DELETE FROM qty_dimension_crosswalks")
            .expect_err("immutable delete trigger");
        assert!(
            delete
                .to_string()
                .contains("quantity dimension crosswalk is immutable")
        );
        let reinsert = ledger
            .conn
            .execute(
                "INSERT INTO qty_dimension_crosswalks(\
                    old_hash, new_hash, source_version, target_version, rule, created_at\
                 ) SELECT old_hash, new_hash, source_version, target_version, rule, created_at \
                   FROM qty_dimension_crosswalks LIMIT 1",
            )
            .expect_err("immutable source reinsert trigger");
        assert!(
            reinsert
                .to_string()
                .contains("quantity dimension crosswalk source is immutable")
        );
        assert_eq!(ledger.table_count("qty_dimension_crosswalks").unwrap(), 1);
        assert!(ledger.lint().unwrap().is_clean());
    }

    #[test]
    fn intact_artifacts_with_a_rebound_target_fail_closed_on_lookup() {
        let ledger = Ledger::open(":memory:").expect("fresh v11 ledger");
        let receipt = retained_receipt(&ledger);
        ledger
            .record_qty_dimension_crosswalk(&receipt)
            .expect("public verified insert");
        let foreign_target = ledger
            .put_artifact("quantity-json-v2", b"intact-but-foreign", None)
            .expect("foreign target artifact")
            .hash;

        ledger
            .conn
            .execute("DROP TRIGGER trg_qty_dimension_crosswalks_immutable_update")
            .expect("open raw-tamper fixture boundary");
        ledger
            .conn
            .prepare("UPDATE qty_dimension_crosswalks SET new_hash = ?1 WHERE old_hash = ?2")
            .expect("prepare raw target rebind")
            .execute_with_params(&[
                blob_param(foreign_target.as_bytes()),
                blob_param(receipt.old_hash().as_bytes()),
            ])
            .expect("inject raw target rebind");
        ledger
            .conn
            .execute(schema::V11.get(2).expect("immutable-update trigger DDL"))
            .expect("restore exact immutable-update guard");

        let Err(LedgerError::Corrupt { hash_hex, detail }) =
            ledger.qty_dimension_crosswalk(receipt.old_hash())
        else {
            panic!("rebound intact target must fail as stored crosswalk corruption");
        };
        assert_eq!(hash_hex, receipt.old_hash().to_hex());
        assert!(detail.contains("qty_dimension_crosswalk.new_hash"));
        assert!(detail.contains("exact canonical six-base output"));
    }

    #[test]
    fn migration_ladder_preserves_v11_through_v14_before_the_v15_batch() {
        assert_eq!(
            schema::MIGRATIONS.len(),
            usize::try_from(SCHEMA_VERSION).unwrap()
        );
        assert_eq!(schema::MIGRATIONS.get(10), Some(&schema::V11));
        assert_eq!(schema::MIGRATIONS.get(11), Some(&schema::V12));
        assert_eq!(schema::MIGRATIONS.get(12), Some(&schema::V13));
        assert_eq!(schema::MIGRATIONS.get(13), Some(&schema::V14));
        assert_eq!(schema::MIGRATIONS.last(), Some(&schema::V15));
    }
}
