//! Canonical runtime artifact for one ordered interface-system card.
//!
//! The ordinary [`NormalizedPack`] already owns bounded claim, observation,
//! joint-statistics, and normalization transport. This wrapper adds the
//! identity-bearing ordered surfaces and system context without inventing a
//! second claim codec. V1 deliberately carries no constitutive model cards;
//! an interface law needs a separately versioned model-pack binding rather
//! than an opaque executable payload.

use fs_blake3::{ContentHash, hash_domain};

use crate::{
    InterfaceSystemCard, MATDB_PACK_TARGET_BASIS, MaterialStateId, NormalizedPack, PackError,
    SurfaceSpec, SystemContext,
};

/// Current normalized interface-pack wire schema.
pub const INTERFACE_PACK_SCHEMA_VERSION: u32 = 1;
/// Coherent numeric basis inherited by the nested claim pack.
pub const INTERFACE_PACK_TARGET_BASIS: &str = MATDB_PACK_TARGET_BASIS;

const MAGIC: &[u8; 8] = b"FSINTPK\0";
const INTERFACE_PACK_HASH_DOMAIN: &str = "org.frankensim.fs-matdb.normalized-interface-pack.v1";
const MAX_INTERFACE_PACK_BYTES: usize = 256 * 1024 * 1024;
const MAX_STRING_BYTES: usize = 1_048_576;

/// Runtime-loadable result of an admitted offline interface-pack compilation.
///
/// `claims_pack` remains the sole owner of property values, provenance,
/// uncertainty, joint statistics, and normalization receipts. `card` binds
/// those claims to the ordered surfaces, medium, optional third body,
/// environment, and named history state.
#[derive(Debug, Clone, PartialEq)]
pub struct NormalizedInterfacePack {
    card: InterfaceSystemCard,
    claims_pack: NormalizedPack,
}

impl NormalizedInterfacePack {
    /// Admit an ordered interface card around an already-admitted claim pack.
    ///
    /// V1 has no model-card argument by design. The nested material-pack
    /// identity and the reconstructed interface-card identity are both stored
    /// and checked during decode.
    pub fn new(
        surface_a: SurfaceSpec,
        surface_b: SurfaceSpec,
        context: SystemContext,
        claims_pack: NormalizedPack,
    ) -> Result<Self, PackError> {
        validate_surface("surface_a", &surface_a)?;
        validate_surface("surface_b", &surface_b)?;
        validate_context(&context)?;
        let card = InterfaceSystemCard::assemble(
            surface_a,
            surface_b,
            context,
            claims_pack.claims().clone(),
            Vec::new(),
        )?;
        let pack = Self { card, claims_pack };
        let encoded_bytes = pack.to_bytes().len();
        if encoded_bytes > MAX_INTERFACE_PACK_BYTES {
            return Err(limit(
                "interface_pack_bytes",
                MAX_INTERFACE_PACK_BYTES,
                encoded_bytes,
            ));
        }
        Ok(pack)
    }

    /// Reconstructed immutable interface-system card.
    #[must_use]
    pub fn card(&self) -> &InterfaceSystemCard {
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
        writer.u32(INTERFACE_PACK_SCHEMA_VERSION);
        encode_surface(&mut writer, self.card.surface_a());
        encode_surface(&mut writer, self.card.surface_b());
        encode_context(&mut writer, self.card.context());
        writer.hash(self.card.content_hash());
        writer.hash(self.claims_pack.content_hash());
        writer.blob(&self.claims_pack.to_bytes());
        writer.bytes
    }

    /// Domain-separated identity of the canonical interface-pack bytes.
    #[must_use]
    pub fn content_hash(&self) -> ContentHash {
        hash_domain(INTERFACE_PACK_HASH_DOMAIN, &self.to_bytes())
    }

    /// Verify an externally pinned whole-artifact identity before decoding.
    pub fn from_bytes_verified(expected: ContentHash, bytes: &[u8]) -> Result<Self, PackError> {
        if bytes.len() > MAX_INTERFACE_PACK_BYTES {
            return Err(limit(
                "interface_pack_bytes",
                MAX_INTERFACE_PACK_BYTES,
                bytes.len(),
            ));
        }
        let actual = hash_domain(INTERFACE_PACK_HASH_DOMAIN, bytes);
        if actual != expected {
            return Err(PackError::IdentityMismatch {
                kind: "interface_pack",
                expected,
                actual,
            });
        }
        Self::from_bytes(bytes)
    }

    /// Decode and semantically re-admit a canonical interface pack.
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, PackError> {
        if bytes.len() > MAX_INTERFACE_PACK_BYTES {
            return Err(limit(
                "interface_pack_bytes",
                MAX_INTERFACE_PACK_BYTES,
                bytes.len(),
            ));
        }
        let mut reader = Reader::new(bytes);
        reader.expect(MAGIC, "normalized interface-pack magic")?;
        let version = reader.u32()?;
        if version != INTERFACE_PACK_SCHEMA_VERSION {
            return Err(reader.malformed(format!(
                "unsupported schema version {version}; expected {INTERFACE_PACK_SCHEMA_VERSION}"
            )));
        }
        let surface_a = decode_surface(&mut reader)?;
        let surface_b = decode_surface(&mut reader)?;
        let context = decode_context(&mut reader)?;
        let expected_card = reader.hash()?;
        let expected_claims_pack = reader.hash()?;
        let claims_bytes = reader.blob("nested_claims_pack", MAX_INTERFACE_PACK_BYTES)?;
        reader.finish()?;

        let claims_pack = NormalizedPack::from_bytes_verified(expected_claims_pack, claims_bytes)?;
        let pack = Self::new(surface_a, surface_b, context, claims_pack)?;
        let actual_card = pack.card.content_hash();
        if actual_card != expected_card {
            return Err(PackError::IdentityMismatch {
                kind: "interface_card",
                expected: expected_card,
                actual: actual_card,
            });
        }
        if pack.to_bytes() != bytes {
            return Err(PackError::Malformed {
                at: 0,
                detail: "interface pack is semantically valid but not canonically encoded"
                    .to_string(),
            });
        }
        Ok(pack)
    }
}

fn validate_surface(label: &str, surface: &SurfaceSpec) -> Result<(), PackError> {
    require_text(
        "interface_surface_material",
        &format!("{label}.material.chemistry"),
        &surface.material.chemistry,
    )?;
    require_text(
        "interface_surface_material",
        &format!("{label}.material.phase"),
        &surface.material.phase,
    )?;
    require_text(
        "interface_surface_material",
        &format!("{label}.material.process"),
        &surface.material.process,
    )?;
    require_text(
        "interface_texture_frame",
        &format!("{label}.texture_frame"),
        &surface.texture_frame,
    )
}

fn validate_context(context: &SystemContext) -> Result<(), PackError> {
    require_text("interface_context", "context.medium", &context.medium)?;
    if let Some(third_body) = &context.third_body {
        require_text("interface_context", "context.third_body", third_body)?;
    }
    require_text(
        "interface_context",
        "context.environment",
        &context.environment,
    )?;
    require_text("interface_context", "context.history", &context.history)
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

fn encode_surface(writer: &mut Writer, surface: &SurfaceSpec) {
    writer.string(&surface.material.chemistry);
    writer.string(&surface.material.phase);
    writer.string(&surface.material.process);
    writer.u32(surface.material.revision);
    writer.string(&surface.texture_frame);
}

fn decode_surface(reader: &mut Reader<'_>) -> Result<SurfaceSpec, PackError> {
    Ok(SurfaceSpec {
        material: MaterialStateId {
            chemistry: reader.string()?,
            phase: reader.string()?,
            process: reader.string()?,
            revision: reader.u32()?,
        },
        texture_frame: reader.string()?,
    })
}

fn encode_context(writer: &mut Writer, context: &SystemContext) {
    writer.string(&context.medium);
    writer.optional_string(context.third_body.as_deref());
    writer.string(&context.environment);
    writer.string(&context.history);
}

fn decode_context(reader: &mut Reader<'_>) -> Result<SystemContext, PackError> {
    Ok(SystemContext {
        medium: reader.string()?,
        third_body: reader.optional_string()?,
        environment: reader.string()?,
        history: reader.string()?,
    })
}

#[derive(Default)]
struct Writer {
    bytes: Vec<u8>,
}

impl Writer {
    fn u8(&mut self, value: u8) {
        self.bytes.push(value);
    }

    fn u32(&mut self, value: u32) {
        self.bytes.extend_from_slice(&value.to_le_bytes());
    }

    fn string(&mut self, value: &str) {
        self.u32(u32::try_from(value.len()).unwrap_or(u32::MAX));
        self.bytes.extend_from_slice(value.as_bytes());
    }

    fn optional_string(&mut self, value: Option<&str>) {
        match value {
            None => self.u8(0),
            Some(value) => {
                self.u8(1);
                self.string(value);
            }
        }
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

    fn u8(&mut self) -> Result<u8, PackError> {
        Ok(self.take(1)?[0])
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
            return Err(limit("interface_string_bytes", MAX_STRING_BYTES, length));
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

    fn optional_string(&mut self) -> Result<Option<String>, PackError> {
        match self.u8()? {
            0 => Ok(None),
            1 => self.string().map(Some),
            tag => Err(self.malformed(format!("unknown optional-string tag {tag}"))),
        }
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
                    "{} trailing bytes after canonical interface pack",
                    self.bytes.len() - self.cursor
                ),
            })
        }
    }
}
