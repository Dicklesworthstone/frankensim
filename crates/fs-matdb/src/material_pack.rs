//! Canonical runtime artifact for one revision-0 material-state card.
//!
//! The ordinary [`NormalizedPack`] already owns bounded claim, observation,
//! joint-statistics, and normalization transport, but it carries no material
//! identity: FSMATPK's `pack_id` is a free-text name, not a
//! [`MaterialStateId`]. This wrapper binds a caller-declared named material
//! state around an already-admitted claim pack so a runtime consumer can
//! reconstruct the exact [`MaterialCard`] a project binding references by
//! content hash — without inventing a second claim codec. V1 deliberately
//! carries no constitutive model cards; a material law needs a separately
//! versioned model-pack binding rather than an opaque executable payload.
//!
//! Claim selection is NOT performed here. A card transports its complete
//! claim set; requirement-driven selection (which claim answers
//! thermal-conductivity at which query point, under which policy) happens at
//! binding time, where it leaves a replayable usage receipt.

use fs_blake3::{ContentHash, hash_domain};

use crate::{MATDB_PACK_TARGET_BASIS, MaterialCard, MaterialStateId, NormalizedPack, PackError};

/// Current normalized material-card-pack wire schema.
pub const MATERIAL_CARD_PACK_SCHEMA_VERSION: u32 = 1;
/// Coherent numeric basis inherited by the nested claim pack.
pub const MATERIAL_CARD_PACK_TARGET_BASIS: &str = MATDB_PACK_TARGET_BASIS;

const MAGIC: &[u8; 8] = b"FSMCDPK\0";
const MATERIAL_CARD_PACK_HASH_DOMAIN: &str =
    "org.frankensim.fs-matdb.normalized-material-card-pack.v1";
const MAX_MATERIAL_CARD_PACK_BYTES: usize = 256 * 1024 * 1024;
const MAX_STRING_BYTES: usize = 1_048_576;

/// Runtime-loadable result of an admitted offline material-card compilation.
///
/// `claims_pack` remains the sole owner of property values, provenance,
/// uncertainty, joint statistics, and normalization receipts. `card` binds
/// those claims to the caller-declared named material state at revision 0.
#[derive(Debug, Clone, PartialEq)]
pub struct NormalizedMaterialCardPack {
    card: MaterialCard,
    claims_pack: NormalizedPack,
}

impl NormalizedMaterialCardPack {
    /// Admit a revision-0 material card around an already-admitted claim pack.
    ///
    /// V1 has no model-card argument by design. The nested claim-pack
    /// identity and the reconstructed material-card identity are both stored
    /// and checked during decode. The declared state must sit at revision 0:
    /// lineage advances only through [`MaterialCard::supersede`], which this
    /// transport wrapper deliberately cannot express.
    ///
    /// # Errors
    /// Refuses a blank chemistry/phase/process, an over-limit identity
    /// string, a nonzero revision (via [`MaterialCard::assemble`]), and an
    /// encoded artifact beyond the byte cap.
    pub fn new(state: MaterialStateId, claims_pack: NormalizedPack) -> Result<Self, PackError> {
        validate_state(&state)?;
        let card = MaterialCard::assemble(state, claims_pack.claims().clone(), Vec::new())?;
        let pack = Self { card, claims_pack };
        let encoded_bytes = pack.to_bytes().len();
        if encoded_bytes > MAX_MATERIAL_CARD_PACK_BYTES {
            return Err(limit(
                "material_card_pack_bytes",
                MAX_MATERIAL_CARD_PACK_BYTES,
                encoded_bytes,
            ));
        }
        Ok(pack)
    }

    /// Reconstructed immutable material card.
    #[must_use]
    pub fn card(&self) -> &MaterialCard {
        &self.card
    }

    /// Nested canonical claim artifact.
    #[must_use]
    pub fn claims_pack(&self) -> &NormalizedPack {
        &self.claims_pack
    }

    /// Stable pack name supplied by the source manifest.
    #[must_use]
    pub fn pack_id(&self) -> &str {
        self.claims_pack.pack_id()
    }

    /// Compiler/version identity that made the admission decisions.
    #[must_use]
    pub fn compiler(&self) -> &str {
        self.claims_pack.compiler()
    }

    /// Hash of the exact raw source envelope.
    #[must_use]
    pub fn source_artifact(&self) -> ContentHash {
        self.claims_pack.source_artifact()
    }

    /// Retained redistribution decision/terms.
    #[must_use]
    pub fn redistribution_terms(&self) -> &str {
        self.claims_pack.redistribution_terms()
    }

    /// Canonical binary representation consumed by L1.
    #[must_use]
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut writer = Writer::default();
        writer.bytes.extend_from_slice(MAGIC);
        writer.u32(MATERIAL_CARD_PACK_SCHEMA_VERSION);
        encode_state(&mut writer, self.card.id());
        writer.hash(self.card.content_hash());
        writer.hash(self.claims_pack.content_hash());
        writer.blob(&self.claims_pack.to_bytes());
        writer.bytes
    }

    /// Domain-separated identity of the canonical material-card-pack bytes.
    #[must_use]
    pub fn content_hash(&self) -> ContentHash {
        hash_domain(MATERIAL_CARD_PACK_HASH_DOMAIN, &self.to_bytes())
    }

    /// Verify an externally pinned whole-artifact identity before decoding.
    pub fn from_bytes_verified(expected: ContentHash, bytes: &[u8]) -> Result<Self, PackError> {
        if bytes.len() > MAX_MATERIAL_CARD_PACK_BYTES {
            return Err(limit(
                "material_card_pack_bytes",
                MAX_MATERIAL_CARD_PACK_BYTES,
                bytes.len(),
            ));
        }
        let actual = hash_domain(MATERIAL_CARD_PACK_HASH_DOMAIN, bytes);
        if actual != expected {
            return Err(PackError::IdentityMismatch {
                kind: "material_card_pack",
                expected,
                actual,
            });
        }
        Self::from_bytes(bytes)
    }

    /// Decode and semantically re-admit a canonical material-card pack.
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, PackError> {
        if bytes.len() > MAX_MATERIAL_CARD_PACK_BYTES {
            return Err(limit(
                "material_card_pack_bytes",
                MAX_MATERIAL_CARD_PACK_BYTES,
                bytes.len(),
            ));
        }
        let mut reader = Reader::new(bytes);
        reader.expect(MAGIC, "normalized material-card-pack magic")?;
        let version = reader.u32()?;
        if version != MATERIAL_CARD_PACK_SCHEMA_VERSION {
            return Err(reader.malformed(format!(
                "unsupported schema version {version}; expected {MATERIAL_CARD_PACK_SCHEMA_VERSION}"
            )));
        }
        let state = decode_state(&mut reader)?;
        let expected_card = reader.hash()?;
        let expected_claims_pack = reader.hash()?;
        let claims_bytes = reader.blob("nested_claims_pack", MAX_MATERIAL_CARD_PACK_BYTES)?;
        reader.finish()?;

        let claims_pack = NormalizedPack::from_bytes_verified(expected_claims_pack, claims_bytes)?;
        let pack = Self::new(state, claims_pack)?;
        let actual_card = pack.card.content_hash();
        if actual_card != expected_card {
            return Err(PackError::IdentityMismatch {
                kind: "material_card",
                expected: expected_card,
                actual: actual_card,
            });
        }
        if pack.to_bytes() != bytes {
            return Err(PackError::Malformed {
                at: 0,
                detail: "material-card pack is semantically valid but not canonically encoded"
                    .to_string(),
            });
        }
        Ok(pack)
    }
}

fn validate_state(state: &MaterialStateId) -> Result<(), PackError> {
    require_text("material_state", "state.chemistry", &state.chemistry)?;
    require_text("material_state", "state.phase", &state.phase)?;
    require_text("material_state", "state.process", &state.process)
}

fn require_text(field: &'static str, label: &str, value: &str) -> Result<(), PackError> {
    if value.trim().is_empty() {
        return Err(invalid(field, format!("{label} must not be blank")));
    }
    if value.len() > MAX_STRING_BYTES {
        return Err(PackError::ResourceLimit {
            resource: field,
            limit: MAX_STRING_BYTES,
            observed: value.len(),
        });
    }
    Ok(())
}

fn invalid(field: &'static str, detail: impl Into<String>) -> PackError {
    PackError::InvalidField {
        field,
        detail: detail.into(),
    }
}

fn limit(resource: &'static str, maximum: usize, observed: usize) -> PackError {
    PackError::ResourceLimit {
        resource,
        limit: maximum,
        observed,
    }
}

fn encode_state(writer: &mut Writer, state: &MaterialStateId) {
    writer.string(&state.chemistry);
    writer.string(&state.phase);
    writer.string(&state.process);
    writer.u32(state.revision);
}

fn decode_state(reader: &mut Reader<'_>) -> Result<MaterialStateId, PackError> {
    Ok(MaterialStateId {
        chemistry: reader.string()?,
        phase: reader.string()?,
        process: reader.string()?,
        revision: reader.u32()?,
    })
}

#[derive(Default)]
struct Writer {
    bytes: Vec<u8>,
}

impl Writer {
    fn u32(&mut self, value: u32) {
        self.bytes.extend_from_slice(&value.to_le_bytes());
    }

    fn string(&mut self, value: &str) {
        self.u32(u32::try_from(value.len()).unwrap_or(u32::MAX));
        self.bytes.extend_from_slice(value.as_bytes());
    }

    fn hash(&mut self, value: ContentHash) {
        self.bytes.extend_from_slice(&value.0);
    }

    fn blob(&mut self, value: &[u8]) {
        self.u32(u32::try_from(value.len()).unwrap_or(u32::MAX));
        self.bytes.extend_from_slice(value);
    }
}

struct Reader<'a> {
    bytes: &'a [u8],
    cursor: usize,
}

impl<'a> Reader<'a> {
    fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, cursor: 0 }
    }

    fn malformed(&self, detail: impl Into<String>) -> PackError {
        PackError::Malformed {
            at: self.cursor,
            detail: detail.into(),
        }
    }

    fn take(&mut self, length: usize) -> Result<&'a [u8], PackError> {
        let end = self
            .cursor
            .checked_add(length)
            .ok_or_else(|| self.malformed("byte offset overflow"))?;
        let slice = self
            .bytes
            .get(self.cursor..end)
            .ok_or_else(|| self.malformed(format!("truncated field needs {length} bytes")))?;
        self.cursor = end;
        Ok(slice)
    }

    fn expect(&mut self, expected: &[u8], name: &str) -> Result<(), PackError> {
        let actual = self.take(expected.len())?;
        if actual == expected {
            Ok(())
        } else {
            Err(self.malformed(format!("invalid {name}")))
        }
    }

    fn u32(&mut self) -> Result<u32, PackError> {
        let bytes: [u8; 4] = self
            .take(4)?
            .try_into()
            .map_err(|_| self.malformed("u32 width"))?;
        Ok(u32::from_le_bytes(bytes))
    }

    fn string(&mut self) -> Result<String, PackError> {
        let raw = self.u32()?;
        let length =
            usize::try_from(raw).map_err(|_| self.malformed("string length does not fit usize"))?;
        if length > MAX_STRING_BYTES {
            return Err(limit(
                "material_card_string_bytes",
                MAX_STRING_BYTES,
                length,
            ));
        }
        let start = self.cursor;
        let bytes = self.take(length)?;
        std::str::from_utf8(bytes)
            .map(str::to_string)
            .map_err(|error| PackError::Malformed {
                at: start + error.valid_up_to(),
                detail: "string field is not UTF-8".to_string(),
            })
    }

    fn hash(&mut self) -> Result<ContentHash, PackError> {
        let bytes: [u8; 32] = self
            .take(32)?
            .try_into()
            .map_err(|_| self.malformed("content-hash width"))?;
        Ok(ContentHash(bytes))
    }

    fn blob(&mut self, resource: &'static str, maximum: usize) -> Result<&'a [u8], PackError> {
        let raw = self.u32()?;
        let length = usize::try_from(raw)
            .map_err(|_| self.malformed(format!("{resource} length does not fit usize")))?;
        if length > maximum {
            return Err(limit(resource, maximum, length));
        }
        self.take(length)
    }

    fn finish(self) -> Result<(), PackError> {
        if self.cursor == self.bytes.len() {
            Ok(())
        } else {
            Err(PackError::Malformed {
                at: self.cursor,
                detail: format!(
                    "{} trailing bytes after canonical material-card pack",
                    self.bytes.len() - self.cursor
                ),
            })
        }
    }
}
