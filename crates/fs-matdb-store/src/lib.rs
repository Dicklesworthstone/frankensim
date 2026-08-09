//! FrankenSQLite-backed queryable material-property store
//! (bead frankensim-oecdy, owner directive 2026-08-08): a "little
//! database of material properties of all kinds for all common
//! materials", built over the fs-matdb pack corpus without ever
//! becoming a second source of truth.
//!
//! ARCHITECTURE — index for search, canonical bytes for evaluation:
//! SQL tables answer the discovery questions (which materials carry a
//! property, which properties a material carries, what is valid at a
//! given ambient point), while every actual property EVALUATION
//! decodes the stored pack bytes through
//! `NormalizedPack::from_bytes_verified` (hash checked first) and
//! delegates to the in-memory `ClaimSet::query` — the SAME evaluator,
//! receipts, and refusals as direct pack use. Parity with the
//! in-memory layer therefore holds by construction: there is exactly
//! one evaluator, and index tampering cannot poison an answer (it can
//! only misdirect discovery, which `verify_index` detects).
//!
//! STALENESS FAILS CLOSED: `seal_corpus` records a domain-separated
//! BLAKE3 digest folded over every pack's content hash in canonical
//! pack-id order. Every discovery/evaluation entry point recomputes
//! the digest from the stored packs and REFUSES
//! (`StoreError::CorpusChanged`) when it disagrees with the seal —
//! a store that drifted from its sealed corpus answers nothing.
//!
//! The pack corpus stays the provenance-bearing source of truth: this
//! crate never compiles TSV sources (that is `xtask matdb-pack`'s
//! fail-closed job) and never fabricates rows; absent data is a named
//! refusal, which is the demand signal of the population strategy in
//! `docs/MATERIAL_PROPERTY_TAXONOMY.md`.
//!
//! Determinism: canonical ingest order, deterministic DDL, and a
//! bitwise-identical rebuild (G5) are asserted in the tests.

use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::{MatDbError, MaterialAnswer, NormalizedPack, QueryPoint, SelectionPolicy};
use fsqlite::{AsyncConnection, FrankenError, Row, SqliteValue};

/// Domain string for the corpus-staleness digest.
const CORPUS_DIGEST_DOMAIN: &str = "org.frankensim.fs-matdb-store.corpus.v1";

/// Store schema version recorded in `PRAGMA user_version`.
pub const STORE_SCHEMA_VERSION: i64 = 1;

/// Typed store errors with stable `FS-MATDB-STORE-*` codes.
#[derive(Debug)]
pub enum StoreError {
    /// Underlying FrankenSQLite failure.
    Sql {
        /// What was being done.
        context: &'static str,
        /// Driver error.
        error: FrankenError,
    },
    /// The sealed corpus digest disagrees with the stored packs.
    CorpusChanged {
        /// Digest recorded by `seal_corpus`.
        sealed: String,
        /// Digest recomputed from the stored packs.
        live: String,
    },
    /// The store has packs but was never sealed (or vice versa).
    NotSealed,
    /// A stored pack failed hash-verified decoding.
    PackCorrupt {
        /// Which pack.
        pack_id: String,
    },
    /// No pack with this id.
    UnknownPack {
        /// Requested id.
        pack_id: String,
    },
    /// A duplicate pack id was ingested.
    DuplicatePack {
        /// Offending id.
        pack_id: String,
    },
    /// An index row disagrees with the decoded pack (tampering or a
    /// store bug) — found by [`MaterialStore::verify_index`].
    IndexMismatch {
        /// Which pack.
        pack_id: String,
        /// What disagreed.
        what: &'static str,
    },
    /// Pack-level admission failure (empty license/redistribution).
    Inadmissible {
        /// Which pack.
        pack_id: String,
        /// Why.
        what: &'static str,
    },
    /// The delegated fs-matdb evaluation refused (extrapolation,
    /// unknown property, ambiguity...) — passed through UNCHANGED.
    MatDb(MatDbError),
    /// No claim in the whole corpus carries this property name — a
    /// typo is not an empty result set.
    UnknownProperty {
        /// Requested property.
        property: String,
    },
    /// A stored value could not be decoded as the expected SQL type.
    Malformed {
        /// Where.
        context: &'static str,
    },
}

impl core::fmt::Display for StoreError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            StoreError::Sql { context, error } => {
                write!(f, "FS-MATDB-STORE-SQL: {context}: {error:?}")
            }
            StoreError::CorpusChanged { sealed, live } => write!(
                f,
                "FS-MATDB-STORE-CORPUS-CHANGED: sealed {sealed} but stored packs hash to {live}; \
                 re-seal deliberately or restore the corpus"
            ),
            StoreError::NotSealed => write!(
                f,
                "FS-MATDB-STORE-NOT-SEALED: seal_corpus must run before queries"
            ),
            StoreError::PackCorrupt { pack_id } => {
                write!(f, "FS-MATDB-STORE-PACK-CORRUPT: {pack_id}")
            }
            StoreError::UnknownPack { pack_id } => {
                write!(f, "FS-MATDB-STORE-UNKNOWN-PACK: {pack_id}")
            }
            StoreError::DuplicatePack { pack_id } => {
                write!(f, "FS-MATDB-STORE-DUPLICATE-PACK: {pack_id}")
            }
            StoreError::IndexMismatch { pack_id, what } => {
                write!(f, "FS-MATDB-STORE-INDEX-MISMATCH: {pack_id}: {what}")
            }
            StoreError::Inadmissible { pack_id, what } => {
                write!(f, "FS-MATDB-STORE-INADMISSIBLE: {pack_id}: {what}")
            }
            StoreError::MatDb(error) => write!(f, "{error}"),
            StoreError::UnknownProperty { property } => {
                write!(f, "FS-MATDB-STORE-UNKNOWN-PROPERTY: {property}")
            }
            StoreError::Malformed { context } => {
                write!(f, "FS-MATDB-STORE-MALFORMED: {context}")
            }
        }
    }
}

impl std::error::Error for StoreError {}

impl From<MatDbError> for StoreError {
    fn from(error: MatDbError) -> Self {
        StoreError::MatDb(error)
    }
}

fn sql_err(context: &'static str) -> impl FnOnce(FrankenError) -> StoreError {
    move |error| StoreError::Sql { context, error }
}

/// One discovery row: a property carried by a material pack.
#[derive(Debug, Clone, PartialEq)]
pub struct PropertyRow {
    /// Owning pack id.
    pub pack_id: String,
    /// Property name.
    pub property: String,
    /// Scalar value (SI) for scalar claims; `None` for curves.
    pub scalar_value: Option<f64>,
    /// `"scalar"` or `"curve"`.
    pub kind: String,
    /// Whether the claim carries at least one observation.
    pub observation_backed: bool,
}

/// One validity axis of a claim.
#[derive(Debug, Clone, PartialEq)]
pub struct ValidityRow {
    /// Axis name (SI post-normalization).
    pub axis: String,
    /// Lower bound.
    pub lo: f64,
    /// Upper bound.
    pub hi: f64,
}

/// The queryable store.
pub struct MaterialStore {
    conn: AsyncConnection,
}

const DDL: &[&str] = &[
    "CREATE TABLE IF NOT EXISTS packs(
        pack_id TEXT PRIMARY KEY,
        content_hash BLOB NOT NULL,
        bytes BLOB NOT NULL,
        compiler TEXT NOT NULL,
        redistribution TEXT NOT NULL
    )",
    "CREATE TABLE IF NOT EXISTS claims(
        pack_id TEXT NOT NULL,
        claim_hash BLOB NOT NULL,
        property TEXT NOT NULL,
        kind TEXT NOT NULL,
        scalar_value REAL,
        observation_backed INTEGER NOT NULL,
        license TEXT NOT NULL,
        source TEXT NOT NULL
    )",
    "CREATE TABLE IF NOT EXISTS validity(
        pack_id TEXT NOT NULL,
        claim_hash BLOB NOT NULL,
        axis TEXT NOT NULL,
        lo REAL NOT NULL,
        hi REAL NOT NULL
    )",
    "CREATE TABLE IF NOT EXISTS corpus_seal(
        singleton INTEGER PRIMARY KEY CHECK (singleton = 1),
        digest BLOB NOT NULL,
        pack_count INTEGER NOT NULL
    )",
    "CREATE INDEX IF NOT EXISTS idx_claims_property ON claims(property)",
    "CREATE INDEX IF NOT EXISTS idx_claims_pack ON claims(pack_id)",
    "CREATE INDEX IF NOT EXISTS idx_validity_claim ON validity(pack_id, claim_hash)",
];

impl MaterialStore {
    /// Open (or create) a store at `path` (`":memory:"` supported) and
    /// apply the schema.
    ///
    /// # Errors
    /// [`StoreError::Sql`] on driver failures.
    pub fn open(path: &str) -> Result<Self, StoreError> {
        let conn = AsyncConnection::open_sync(path).map_err(sql_err("open"))?;
        for ddl in DDL {
            conn.execute_sync(ddl).map_err(sql_err("schema ddl"))?;
        }
        conn.execute_sync(&format!("PRAGMA user_version = {STORE_SCHEMA_VERSION}"))
            .map_err(sql_err("schema version"))?;
        Ok(MaterialStore { conn })
    }

    /// Ingest one verified pack: canonical bytes plus derived index
    /// rows. Refuses duplicates and packs with empty license or
    /// redistribution terms (the license gate survives the store).
    ///
    /// # Errors
    /// [`StoreError`] on duplicates, inadmissible metadata, or driver
    /// failures.
    pub fn ingest_pack(&self, pack: &NormalizedPack) -> Result<(), StoreError> {
        let pack_id = pack.pack_id().to_string();
        if pack.redistribution_terms().trim().is_empty() {
            return Err(StoreError::Inadmissible {
                pack_id,
                what: "empty redistribution terms",
            });
        }
        // Admission PRE-PASS before any row is written (review finding:
        // a mid-ingest refusal must not leave a committed partial pack).
        for (_claim_id, claim) in pack.claims().claims_ordered() {
            if claim.provenance.license.trim().is_empty() {
                return Err(StoreError::Inadmissible {
                    pack_id,
                    what: "claim with empty license",
                });
            }
        }
        let existing = self
            .conn
            .query_with_params_sync(
                "SELECT 1 FROM packs WHERE pack_id = ?1",
                &[text_param(&pack_id)],
            )
            .map_err(sql_err("duplicate probe"))?;
        if !existing.is_empty() {
            return Err(StoreError::DuplicatePack { pack_id });
        }
        let bytes = pack.to_bytes();
        let hash = pack.content_hash();
        // Transactional ingest: any failure below rolls the whole pack
        // back, so retries never hit DuplicatePack on a half-ingested
        // id and the index can never be partially populated.
        self.conn
            .execute_sync("BEGIN IMMEDIATE")
            .map_err(sql_err("begin ingest"))?;
        let result = self.ingest_rows(&pack_id, pack, &bytes, hash);
        match &result {
            Ok(()) => {
                self.conn
                    .execute_sync("COMMIT")
                    .map_err(sql_err("commit ingest"))?;
            }
            Err(_) => {
                let _ = self.conn.execute_sync("ROLLBACK");
            }
        }
        result
    }

    fn ingest_rows(
        &self,
        pack_id: &str,
        pack: &NormalizedPack,
        bytes: &[u8],
        hash: ContentHash,
    ) -> Result<(), StoreError> {
        self.conn
            .execute_with_params_sync(
                "INSERT INTO packs(pack_id, content_hash, bytes, compiler, redistribution) \
                 VALUES (?1, ?2, ?3, ?4, ?5)",
                &[
                    text_param(pack_id),
                    blob_param(hash.as_bytes()),
                    blob_param(bytes),
                    text_param(pack.compiler()),
                    text_param(pack.redistribution_terms()),
                ],
            )
            .map_err(sql_err("insert pack"))?;
        // Derived index rows, in the ClaimSet's canonical order.
        for (claim_id, claim) in pack.claims().claims_ordered() {
            let (kind, scalar_value) = match &claim.value {
                fs_matdb::PropertyValue::Scalar { value, .. } => ("scalar", Some(*value)),
                fs_matdb::PropertyValue::Curve { .. } => ("curve", None),
            };
            self.conn
                .execute_with_params_sync(
                    "INSERT INTO claims(pack_id, claim_hash, property, kind, scalar_value, \
                     observation_backed, license, source) VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8)",
                    &[
                        text_param(pack_id),
                        blob_param(claim_id.0.as_bytes()),
                        text_param(claim.key.name()),
                        text_param(kind),
                        scalar_value.map_or(SqliteValue::Null, SqliteValue::Float),
                        SqliteValue::Integer(i64::from(!claim.observations.is_empty())),
                        text_param(&claim.provenance.license),
                        text_param(&claim.provenance.source),
                    ],
                )
                .map_err(sql_err("insert claim"))?;
            for (axis, (lo, hi)) in claim.validity.bounds() {
                self.conn
                    .execute_with_params_sync(
                        "INSERT INTO validity(pack_id, claim_hash, axis, lo, hi) \
                         VALUES (?1, ?2, ?3, ?4, ?5)",
                        &[
                            text_param(pack_id),
                            blob_param(claim_id.0.as_bytes()),
                            text_param(axis),
                            SqliteValue::Float(*lo),
                            SqliteValue::Float(*hi),
                        ],
                    )
                    .map_err(sql_err("insert validity"))?;
            }
        }
        Ok(())
    }

    fn live_digest(&self) -> Result<(ContentHash, i64), StoreError> {
        let rows = self
            .conn
            .query_sync("SELECT pack_id, content_hash FROM packs ORDER BY pack_id")
            .map_err(sql_err("digest scan"))?;
        let mut hasher = DomainHasher::new(CORPUS_DIGEST_DOMAIN);
        for row in &rows {
            let id = row_text(row, 0, "packs.pack_id")?;
            let hash = row_blob(row, 1, "packs.content_hash")?;
            hasher.update(&(id.len() as u64).to_le_bytes());
            hasher.update(id.as_bytes());
            hasher.update(&(hash.len() as u64).to_le_bytes());
            hasher.update(&hash);
        }
        Ok((
            hasher.finalize(),
            i64::try_from(rows.len()).expect("pack count fits"),
        ))
    }

    /// Record the corpus digest. Queries refuse until this runs, and
    /// refuse again if the stored packs later drift from the seal.
    ///
    /// # Errors
    /// [`StoreError::Sql`] on driver failures.
    pub fn seal_corpus(&self) -> Result<ContentHash, StoreError> {
        let (digest, count) = self.live_digest()?;
        self.conn
            .execute_with_params_sync(
                "INSERT INTO corpus_seal(singleton, digest, pack_count) VALUES (1, ?1, ?2) \
                 ON CONFLICT(singleton) DO UPDATE SET digest = ?1, pack_count = ?2",
                &[blob_param(digest.as_bytes()), SqliteValue::Integer(count)],
            )
            .map_err(sql_err("seal"))?;
        Ok(digest)
    }

    /// The staleness gate every entry point calls: recompute the live
    /// digest and compare with the seal.
    ///
    /// # Errors
    /// [`StoreError::NotSealed`] / [`StoreError::CorpusChanged`].
    pub fn require_sealed(&self) -> Result<(), StoreError> {
        let rows = self
            .conn
            .query_sync("SELECT digest FROM corpus_seal WHERE singleton = 1")
            .map_err(sql_err("read seal"))?;
        let Some(row) = rows.first() else {
            return Err(StoreError::NotSealed);
        };
        let sealed = row_blob(row, 0, "corpus_seal.digest")?;
        let (live, _) = self.live_digest()?;
        if sealed != live.as_bytes() {
            return Err(StoreError::CorpusChanged {
                sealed: hex(&sealed),
                live: live.to_hex(),
            });
        }
        Ok(())
    }

    /// Discovery: every property row of one material pack.
    ///
    /// # Errors
    /// Staleness refusals plus driver failures.
    pub fn properties_of(&self, pack_id: &str) -> Result<Vec<PropertyRow>, StoreError> {
        self.require_sealed()?;
        let rows = self
            .conn
            .query_with_params_sync(
                "SELECT pack_id, property, scalar_value, kind, observation_backed \
                 FROM claims WHERE pack_id = ?1 ORDER BY property",
                &[text_param(pack_id)],
            )
            .map_err(sql_err("properties_of"))?;
        if rows.is_empty() {
            let known = self
                .conn
                .query_with_params_sync(
                    "SELECT 1 FROM packs WHERE pack_id = ?1",
                    &[text_param(pack_id)],
                )
                .map_err(sql_err("pack probe"))?;
            if known.is_empty() {
                return Err(StoreError::UnknownPack {
                    pack_id: pack_id.to_string(),
                });
            }
        }
        rows.iter().map(property_row).collect()
    }

    /// Discovery: every material carrying `property`, optionally
    /// filtered to claims whose scalar value lies in `[lo, hi]`.
    ///
    /// # Errors
    /// Staleness refusals plus driver failures.
    pub fn materials_with(
        &self,
        property: &str,
        scalar_range: Option<(f64, f64)>,
    ) -> Result<Vec<PropertyRow>, StoreError> {
        self.require_sealed()?;
        let rows = match scalar_range {
            None => self
                .conn
                .query_with_params_sync(
                    "SELECT pack_id, property, scalar_value, kind, observation_backed \
                     FROM claims WHERE property = ?1 ORDER BY pack_id",
                    &[text_param(property)],
                )
                .map_err(sql_err("materials_with"))?,
            Some((lo, hi)) => self
                .conn
                .query_with_params_sync(
                    "SELECT pack_id, property, scalar_value, kind, observation_backed \
                     FROM claims WHERE property = ?1 AND scalar_value IS NOT NULL \
                     AND scalar_value >= ?2 AND scalar_value <= ?3 ORDER BY pack_id",
                    &[
                        text_param(property),
                        SqliteValue::Float(lo),
                        SqliteValue::Float(hi),
                    ],
                )
                .map_err(sql_err("materials_with range"))?,
        };
        if rows.is_empty() {
            self.require_known_property(property)?;
        }
        rows.iter().map(property_row).collect()
    }

    fn require_known_property(&self, property: &str) -> Result<(), StoreError> {
        let any = self
            .conn
            .query_with_params_sync(
                "SELECT 1 FROM claims WHERE property = ?1 LIMIT 1",
                &[text_param(property)],
            )
            .map_err(sql_err("property probe"))?;
        if any.is_empty() {
            return Err(StoreError::UnknownProperty {
                property: property.to_string(),
            });
        }
        Ok(())
    }

    /// Discovery: materials whose `property` claim is VALID at the
    /// given axis point (axis present with `lo <= value <= hi`, or no
    /// bound on that axis at all — matching `ValidityDomain`
    /// semantics where a missing axis is unconstrained).
    ///
    /// # Errors
    /// Staleness refusals plus driver failures.
    pub fn valid_at(
        &self,
        property: &str,
        axis: &str,
        value: f64,
    ) -> Result<Vec<PropertyRow>, StoreError> {
        self.require_sealed()?;
        let rows = self
            .conn
            .query_with_params_sync(
                "SELECT c.pack_id, c.property, c.scalar_value, c.kind, c.observation_backed \
                 FROM claims c WHERE c.property = ?1 AND NOT EXISTS (\
                    SELECT 1 FROM validity v WHERE v.pack_id = c.pack_id \
                    AND v.claim_hash = c.claim_hash AND v.axis = ?2 \
                    AND (v.lo > ?3 OR v.hi < ?3)\
                 ) ORDER BY c.pack_id",
                &[
                    text_param(property),
                    text_param(axis),
                    SqliteValue::Float(value),
                ],
            )
            .map_err(sql_err("valid_at"))?;
        if rows.is_empty() {
            self.require_known_property(property)?;
        }
        rows.iter().map(property_row).collect()
    }

    /// Load a pack's canonical bytes and decode them HASH-VERIFIED,
    /// behind the staleness gate like every other surface.
    ///
    /// # Errors
    /// Staleness refusals plus [`StoreError::UnknownPack`] /
    /// [`StoreError::PackCorrupt`].
    pub fn load_pack(&self, pack_id: &str) -> Result<NormalizedPack, StoreError> {
        self.require_sealed()?;
        self.load_pack_unchecked(pack_id)
    }

    /// The gate-free decode used internally after a caller has already
    /// paid `require_sealed` (and by `verify_index`, which must work on
    /// a DRIFTED store precisely because its job is investigating one).
    fn load_pack_unchecked(&self, pack_id: &str) -> Result<NormalizedPack, StoreError> {
        let rows = self
            .conn
            .query_with_params_sync(
                "SELECT content_hash, bytes FROM packs WHERE pack_id = ?1",
                &[text_param(pack_id)],
            )
            .map_err(sql_err("load pack"))?;
        let Some(row) = rows.first() else {
            return Err(StoreError::UnknownPack {
                pack_id: pack_id.to_string(),
            });
        };
        let hash_bytes = row_blob(row, 0, "packs.content_hash")?;
        let bytes = row_blob(row, 1, "packs.bytes")?;
        let expected = ContentHash::from_slice(&hash_bytes).ok_or(StoreError::Malformed {
            context: "packs.content_hash length",
        })?;
        NormalizedPack::from_bytes_verified(expected, &bytes).map_err(|_| StoreError::PackCorrupt {
            pack_id: pack_id.to_string(),
        })
    }

    /// AUTHORITATIVE evaluation: decode the hash-verified pack and
    /// delegate to the in-memory receipted query — the same evaluator,
    /// receipt, and refusal semantics as direct pack use. Index
    /// tampering cannot reach this path.
    ///
    /// # Errors
    /// Staleness refusals; [`StoreError::MatDb`] passing through every
    /// fs-matdb refusal (extrapolation, unknown property, ambiguity)
    /// unchanged.
    pub fn evaluate(
        &self,
        pack_id: &str,
        property: &str,
        point: &QueryPoint,
        policy: SelectionPolicy,
    ) -> Result<MaterialAnswer, StoreError> {
        self.require_sealed()?;
        let pack = self.load_pack_unchecked(pack_id)?;
        let answer = pack.claims().query(property, point, policy)?;
        pack.claims().verify_receipt(&answer.receipt)?;
        Ok(answer)
    }

    /// Cross-check every index row of a pack against its decoded
    /// canonical bytes — the tamper/bug detector for the derived
    /// tables.
    ///
    /// # Errors
    /// [`StoreError::IndexMismatch`] naming the first disagreement.
    pub fn verify_index(&self, pack_id: &str) -> Result<(), StoreError> {
        let pack = self.load_pack_unchecked(pack_id)?;
        let index_rows = self
            .conn
            .query_with_params_sync(
                "SELECT property, kind, scalar_value, claim_hash FROM claims \
                 WHERE pack_id = ?1 ORDER BY claim_hash",
                &[text_param(pack_id)],
            )
            .map_err(sql_err("verify_index scan"))?;
        let mut derived: Vec<(String, String, Option<f64>, Vec<u8>)> = pack
            .claims()
            .claims_ordered()
            .map(|(claim_id, claim)| {
                let (kind, value) = match &claim.value {
                    fs_matdb::PropertyValue::Scalar { value, .. } => {
                        ("scalar".to_string(), Some(*value))
                    }
                    fs_matdb::PropertyValue::Curve { .. } => ("curve".to_string(), None),
                };
                (
                    claim.key.name().to_string(),
                    kind,
                    value,
                    claim_id.0.as_bytes().to_vec(),
                )
            })
            .collect();
        derived.sort_by(|a, b| a.3.cmp(&b.3));
        if index_rows.len() != derived.len() {
            return Err(StoreError::IndexMismatch {
                pack_id: pack_id.to_string(),
                what: "claim count",
            });
        }
        for (row, expect) in index_rows.iter().zip(derived.iter()) {
            if row_text(row, 0, "claims.property")? != expect.0 {
                return Err(StoreError::IndexMismatch {
                    pack_id: pack_id.to_string(),
                    what: "property name",
                });
            }
            if row_text(row, 1, "claims.kind")? != expect.1 {
                return Err(StoreError::IndexMismatch {
                    pack_id: pack_id.to_string(),
                    what: "claim kind",
                });
            }
            let stored_value = match row.get(2) {
                Some(SqliteValue::Float(v)) => Some(*v),
                Some(SqliteValue::Null) | None => None,
                _ => {
                    return Err(StoreError::Malformed {
                        context: "claims.scalar_value type",
                    });
                }
            };
            let matches = match (stored_value, expect.2) {
                (Some(a), Some(b)) => a.to_bits() == b.to_bits(),
                (None, None) => true,
                _ => false,
            };
            if !matches {
                return Err(StoreError::IndexMismatch {
                    pack_id: pack_id.to_string(),
                    what: "scalar value",
                });
            }
            if row_blob(row, 3, "claims.claim_hash")? != expect.3 {
                return Err(StoreError::IndexMismatch {
                    pack_id: pack_id.to_string(),
                    what: "claim hash",
                });
            }
        }
        Ok(())
    }

    /// Canonical dump of the derived tables (for the bitwise-rebuild
    /// determinism proof): every claims + validity row rendered in a
    /// fixed order and format.
    ///
    /// # Errors
    /// Driver failures.
    pub fn canonical_dump(&self) -> Result<String, StoreError> {
        let mut out = String::new();
        let rows = self
            .conn
            .query_sync(
                "SELECT pack_id, property, kind, ifnull(scalar_value, 'NULL'), \
                 observation_backed, lower(hex(claim_hash)) FROM claims \
                 ORDER BY pack_id, claim_hash",
            )
            .map_err(sql_err("dump claims"))?;
        for row in &rows {
            for i in 0..6 {
                out.push_str(&render_value(row.get(i)));
                out.push('|');
            }
            out.push('\n');
        }
        let rows = self
            .conn
            .query_sync(
                "SELECT pack_id, lower(hex(claim_hash)), axis, lo, hi FROM validity \
                 ORDER BY pack_id, claim_hash, axis",
            )
            .map_err(sql_err("dump validity"))?;
        for row in &rows {
            for i in 0..5 {
                out.push_str(&render_value(row.get(i)));
                out.push('|');
            }
            out.push('\n');
        }
        Ok(out)
    }
}

fn render_value(value: Option<&SqliteValue>) -> String {
    match value {
        Some(SqliteValue::Text(s)) => s.to_string(),
        Some(SqliteValue::Integer(v)) => v.to_string(),
        Some(SqliteValue::Float(v)) => format!("{:016x}", v.to_bits()),
        Some(SqliteValue::Blob(b)) => hex(b),
        Some(SqliteValue::Null) | None => "NULL".to_string(),
    }
}

fn text_param(s: &str) -> SqliteValue {
    SqliteValue::Text(s.into())
}

fn blob_param(b: &[u8]) -> SqliteValue {
    SqliteValue::Blob(b.to_vec().into())
}

fn hex(bytes: &[u8]) -> String {
    use core::fmt::Write as _;
    let mut out = String::with_capacity(bytes.len() * 2);
    for b in bytes {
        write!(out, "{b:02x}").expect("writing to String cannot fail");
    }
    out
}

fn row_text(row: &Row, idx: usize, context: &'static str) -> Result<String, StoreError> {
    match row.get(idx) {
        Some(SqliteValue::Text(s)) => Ok(s.to_string()),
        _ => Err(StoreError::Malformed { context }),
    }
}

fn row_blob(row: &Row, idx: usize, context: &'static str) -> Result<Vec<u8>, StoreError> {
    match row.get(idx) {
        Some(SqliteValue::Blob(b)) => Ok(b.to_vec()),
        _ => Err(StoreError::Malformed { context }),
    }
}

fn property_row(row: &Row) -> Result<PropertyRow, StoreError> {
    let scalar_value = match row.get(2) {
        Some(SqliteValue::Float(v)) => Some(*v),
        Some(SqliteValue::Null) | None => None,
        _ => {
            return Err(StoreError::Malformed {
                context: "claims.scalar_value",
            });
        }
    };
    let observation_backed = match row.get(4) {
        Some(SqliteValue::Integer(v)) => *v != 0,
        _ => {
            return Err(StoreError::Malformed {
                context: "claims.observation_backed",
            });
        }
    };
    Ok(PropertyRow {
        pack_id: row_text(row, 0, "claims.pack_id")?,
        property: row_text(row, 1, "claims.property")?,
        scalar_value,
        kind: row_text(row, 3, "claims.kind")?,
        observation_backed,
    })
}
