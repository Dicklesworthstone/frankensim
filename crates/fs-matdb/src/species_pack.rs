//! Canonical runtime transport for one source-declared species association.
//!
//! [`NormalizedSpeciesPack`] carries source-declared thermochemical context
//! alongside separately validated model cards without pretending that the
//! generic constitutive-card schema contains those fields or directly links the
//! artifacts. The pack is data only: it retains the source declaration and
//! normalization evidence but does not authenticate chemistry or evaluate a
//! model.

use fs_blake3::{ContentHash, hash_domain};
use fs_qty::Dims;
use fs_qty::chemistry::SpeciesId;

use crate::{PackError, Provenance};

/// Current normalized species-association wire schema.
pub const SPECIES_PACK_SCHEMA_VERSION: u32 = 1;
/// Canonical coherent-SI target basis and explicit dimension order.
pub const SPECIES_PACK_TARGET_BASIS: &str = "SI-six-base[m,kg,s,K,A,mol]";
/// Coherent-SI dimensions of molar mass, kilograms per mole.
pub const SPECIES_MOLAR_MASS_DIMS: Dims = Dims([0, 1, 0, 0, 0, -1]);
/// Coherent-SI dimensions of reference pressure.
pub const SPECIES_REFERENCE_PRESSURE_DIMS: Dims = Dims([-1, 1, -2, 0, 0, 0]);

const MAGIC: &[u8; 8] = b"FSSPCPK\0";
const SPECIES_PACK_HASH_DOMAIN: &str = "org.frankensim.fs-matdb.normalized-species-pack.v1";
const MAX_PACK_BYTES: usize = 16 * 1024 * 1024;
const MAX_STRING_BYTES: usize = 1_048_576;
const MAX_SOURCES: usize = 4_096;
const MAX_ELEMENTAL_REFERENCE_BYTES: usize = 128;
const NORMALIZATION_COUNT: usize = 2;

/// Numeric field normalized by the offline species compiler.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum SpeciesNormalizationTarget {
    /// Positive molar mass in coherent kilograms per mole.
    MolarMass,
    /// Positive standard-state reference pressure in pascals.
    ReferencePressure,
}

/// Exact unit/basis transform for one species-association numeric field.
#[derive(Debug, Clone, PartialEq)]
pub struct SpeciesNormalizationReceipt {
    target: SpeciesNormalizationTarget,
    source_literal: ContentHash,
    dims: Dims,
    scale: f64,
    offset: f64,
    source_basis: String,
    target_basis: String,
}

impl SpeciesNormalizationReceipt {
    /// Construct an immutable normalization receipt.
    #[must_use]
    pub fn new(
        target: SpeciesNormalizationTarget,
        source_literal: ContentHash,
        dims: Dims,
        scale: f64,
        offset: f64,
        source_basis: impl Into<String>,
        target_basis: impl Into<String>,
    ) -> Self {
        Self {
            target,
            source_literal,
            dims,
            scale,
            offset,
            source_basis: source_basis.into(),
            target_basis: target_basis.into(),
        }
    }

    /// Structurally linked numeric field.
    #[must_use]
    pub const fn target(&self) -> SpeciesNormalizationTarget {
        self.target
    }

    /// Hash of the exact source literal and unit token.
    #[must_use]
    pub const fn source_literal(&self) -> ContentHash {
        self.source_literal
    }

    /// Six-base dimensions of the normalized field.
    #[must_use]
    pub const fn dims(&self) -> Dims {
        self.dims
    }

    /// Multiplicative term in `si = source * scale + offset`.
    #[must_use]
    pub const fn scale(&self) -> f64 {
        self.scale
    }

    /// Additive term in `si = source * scale + offset`.
    #[must_use]
    pub const fn offset(&self) -> f64 {
        self.offset
    }

    /// Explicit source unit/basis expression.
    #[must_use]
    pub fn source_basis(&self) -> &str {
        &self.source_basis
    }

    /// Canonical normalized target basis.
    #[must_use]
    pub fn target_basis(&self) -> &str {
        &self.target_basis
    }
}

/// Immutable source-declared thermochemical association for one species.
///
/// Version 1 deliberately admits only the exact gas/ideal-gas convention used
/// by the first downstream standard-state evaluator. The elemental-reference
/// id remains an opaque source-declared convention name.
#[derive(Debug, Clone, PartialEq)]
pub struct SpeciesAssociation {
    species: SpeciesId,
    molar_mass: f64,
    standard_state_phase: String,
    reference_eos: String,
    reference_pressure: f64,
    elemental_reference: String,
    sources: Vec<ContentHash>,
    provenance: Provenance,
}

impl SpeciesAssociation {
    /// Admit one bounded, provenance-complete source association.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        species: SpeciesId,
        molar_mass: f64,
        standard_state_phase: impl Into<String>,
        reference_eos: impl Into<String>,
        reference_pressure: f64,
        elemental_reference: impl Into<String>,
        mut sources: Vec<ContentHash>,
        provenance: Provenance,
    ) -> Result<Self, PackError> {
        require_positive("species.molar_mass", molar_mass)?;
        require_positive("species.reference_pressure", reference_pressure)?;
        let standard_state_phase = standard_state_phase.into();
        let reference_eos = reference_eos.into();
        let elemental_reference = elemental_reference.into();
        if standard_state_phase != "gas" {
            return Err(invalid(
                "species.standard_state_phase",
                "species pack v1 admits exactly the explicit 'gas' phase",
            ));
        }
        if reference_eos != "ideal-gas" {
            return Err(invalid(
                "species.reference_eos",
                "species pack v1 admits exactly the explicit 'ideal-gas' reference EOS",
            ));
        }
        require_elemental_reference(&elemental_reference)?;
        if sources.is_empty() {
            return Err(invalid(
                "species.sources",
                "a portable species association requires at least one source artifact",
            ));
        }
        if sources.len() > MAX_SOURCES {
            return Err(limit("species_sources", MAX_SOURCES, sources.len()));
        }
        sources.sort_unstable();
        if sources.windows(2).any(|pair| pair[0] == pair[1]) {
            return Err(invalid(
                "species.sources",
                "source artifacts must be unique",
            ));
        }
        require_text("species.provenance.source", &provenance.source)?;
        require_text("species.provenance.license", &provenance.license)?;
        let artifact = provenance.artifact.ok_or_else(|| {
            invalid(
                "species.provenance.artifact",
                "a portable species association requires an exact provenance artifact",
            )
        })?;
        if sources.binary_search(&artifact).is_err() {
            return Err(invalid(
                "species.provenance.artifact",
                "the provenance artifact must be one of the retained source artifacts",
            ));
        }
        Ok(Self {
            species,
            molar_mass,
            standard_state_phase,
            reference_eos,
            reference_pressure,
            elemental_reference,
            sources,
            provenance,
        })
    }

    /// Exact canonical species identifier.
    #[must_use]
    pub const fn species(&self) -> &SpeciesId {
        &self.species
    }

    /// Positive coherent-SI molar mass in kilograms per mole.
    #[must_use]
    pub const fn molar_mass(&self) -> f64 {
        self.molar_mass
    }

    /// Explicit source-declared standard-state phase (`gas` in v1).
    #[must_use]
    pub fn standard_state_phase(&self) -> &str {
        &self.standard_state_phase
    }

    /// Explicit source-declared reference EOS (`ideal-gas` in v1).
    #[must_use]
    pub fn reference_eos(&self) -> &str {
        &self.reference_eos
    }

    /// Positive coherent-SI standard-state pressure in pascals.
    #[must_use]
    pub const fn reference_pressure(&self) -> f64 {
        self.reference_pressure
    }

    /// Opaque elemental-reference convention identifier.
    #[must_use]
    pub fn elemental_reference(&self) -> &str {
        &self.elemental_reference
    }

    /// Canonically ordered source and record artifacts.
    #[must_use]
    pub fn sources(&self) -> &[ContentHash] {
        &self.sources
    }

    /// Exact source citation, license, and artifact provenance.
    #[must_use]
    pub const fn provenance(&self) -> &Provenance {
        &self.provenance
    }
}

/// Runtime-loadable result of an admitted offline species compilation.
#[derive(Debug, Clone, PartialEq)]
pub struct NormalizedSpeciesPack {
    pack_id: String,
    compiler: String,
    source_artifact: ContentHash,
    redistribution_terms: String,
    association: SpeciesAssociation,
    normalizations: Vec<SpeciesNormalizationReceipt>,
}

impl NormalizedSpeciesPack {
    /// Admit one association and complete numeric normalization evidence.
    pub fn new(
        pack_id: impl Into<String>,
        compiler: impl Into<String>,
        source_artifact: ContentHash,
        redistribution_terms: impl Into<String>,
        association: SpeciesAssociation,
        mut normalizations: Vec<SpeciesNormalizationReceipt>,
    ) -> Result<Self, PackError> {
        let pack_id = pack_id.into();
        let compiler = compiler.into();
        let redistribution_terms = redistribution_terms.into();
        require_text("pack_id", &pack_id)?;
        require_text("compiler", &compiler)?;
        require_text("redistribution_terms", &redistribution_terms)?;
        if pack_id != association.species.as_str() {
            return Err(invalid(
                "pack_id",
                format!(
                    "species pack id {pack_id:?} must equal canonical SpeciesId {:?}",
                    association.species.as_str()
                ),
            ));
        }
        normalizations.sort_by_key(SpeciesNormalizationReceipt::target);
        if normalizations.len() != NORMALIZATION_COUNT {
            return Err(invalid(
                "species_normalizations",
                "species pack v1 requires exactly molar-mass and reference-pressure receipts",
            ));
        }
        if normalizations[0].target != SpeciesNormalizationTarget::MolarMass
            || normalizations[1].target != SpeciesNormalizationTarget::ReferencePressure
        {
            return Err(invalid(
                "species_normalizations",
                "normalization targets must cover molar mass and reference pressure exactly once",
            ));
        }
        for receipt in &normalizations {
            validate_receipt(receipt)?;
        }
        let pack = Self {
            pack_id,
            compiler,
            source_artifact,
            redistribution_terms,
            association,
            normalizations,
        };
        let encoded_len = pack.to_bytes().len();
        if encoded_len > MAX_PACK_BYTES {
            return Err(limit("species_pack_bytes", MAX_PACK_BYTES, encoded_len));
        }
        Ok(pack)
    }

    /// Stable pack name, exactly equal to the retained species id.
    #[must_use]
    pub fn pack_id(&self) -> &str {
        &self.pack_id
    }

    /// Compiler/version identity that made the admission decisions.
    #[must_use]
    pub fn compiler(&self) -> &str {
        &self.compiler
    }

    /// Hash of the exact raw source envelope.
    #[must_use]
    pub const fn source_artifact(&self) -> ContentHash {
        self.source_artifact
    }

    /// Retained redistribution decision and terms.
    #[must_use]
    pub fn redistribution_terms(&self) -> &str {
        &self.redistribution_terms
    }

    /// One immutable source-declared species association.
    #[must_use]
    pub const fn association(&self) -> &SpeciesAssociation {
        &self.association
    }

    /// Canonically ordered, complete numeric normalization receipts.
    #[must_use]
    pub fn normalizations(&self) -> &[SpeciesNormalizationReceipt] {
        &self.normalizations
    }

    /// Canonical binary representation consumed by L1.
    #[must_use]
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut writer = Writer::default();
        writer.bytes.extend_from_slice(MAGIC);
        writer.u32(SPECIES_PACK_SCHEMA_VERSION);
        writer.string(&self.pack_id);
        writer.string(&self.compiler);
        writer.hash(self.source_artifact);
        writer.string(&self.redistribution_terms);
        writer.string(self.association.species.as_str());
        writer.f64(self.association.molar_mass);
        writer.string(&self.association.standard_state_phase);
        writer.string(&self.association.reference_eos);
        writer.f64(self.association.reference_pressure);
        writer.string(&self.association.elemental_reference);
        writer.count(self.association.sources.len());
        for source in &self.association.sources {
            writer.hash(*source);
        }
        writer.string(&self.association.provenance.source);
        writer.string(&self.association.provenance.license);
        writer.optional_hash(self.association.provenance.artifact);
        writer.count(self.normalizations.len());
        for receipt in &self.normalizations {
            writer.u8(match receipt.target {
                SpeciesNormalizationTarget::MolarMass => 0,
                SpeciesNormalizationTarget::ReferencePressure => 1,
            });
            writer.hash(receipt.source_literal);
            writer.dims(receipt.dims);
            writer.f64(receipt.scale);
            writer.f64(receipt.offset);
            writer.string(&receipt.source_basis);
            writer.string(&receipt.target_basis);
        }
        writer.bytes
    }

    /// Domain-separated identity of the canonical pack bytes.
    #[must_use]
    pub fn content_hash(&self) -> ContentHash {
        hash_domain(SPECIES_PACK_HASH_DOMAIN, &self.to_bytes())
    }

    /// Verify an externally pinned whole-pack identity before decoding.
    pub fn from_bytes_verified(expected: ContentHash, bytes: &[u8]) -> Result<Self, PackError> {
        if bytes.len() > MAX_PACK_BYTES {
            return Err(limit("species_pack_bytes", MAX_PACK_BYTES, bytes.len()));
        }
        let actual = hash_domain(SPECIES_PACK_HASH_DOMAIN, bytes);
        if actual != expected {
            return Err(PackError::IdentityMismatch {
                kind: "species pack",
                expected,
                actual,
            });
        }
        Self::from_bytes(bytes)
    }

    /// Decode and semantically re-admit one canonical species pack.
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, PackError> {
        if bytes.len() > MAX_PACK_BYTES {
            return Err(limit("species_pack_bytes", MAX_PACK_BYTES, bytes.len()));
        }
        let mut reader = Reader::new(bytes);
        reader.expect(MAGIC, "normalized species-pack magic")?;
        let version = reader.u32()?;
        if version != SPECIES_PACK_SCHEMA_VERSION {
            return Err(reader.malformed(format!(
                "unsupported schema version {version}; expected {SPECIES_PACK_SCHEMA_VERSION}"
            )));
        }
        let pack_id = reader.string()?;
        let compiler = reader.string()?;
        let source_artifact = reader.hash()?;
        let redistribution_terms = reader.string()?;
        let species_text = reader.string()?;
        let species = SpeciesId::new(species_text).map_err(|error| {
            invalid(
                "species.id",
                format!("fs-qty refused decoded SpeciesId: {error}"),
            )
        })?;
        let molar_mass = reader.f64()?;
        let standard_state_phase = reader.string()?;
        let reference_eos = reader.string()?;
        let reference_pressure = reader.f64()?;
        let elemental_reference = reader.string()?;
        let source_count = reader.count("species_sources", MAX_SOURCES)?;
        reader.require_items(source_count, 32, "species source hashes")?;
        let mut sources = Vec::with_capacity(source_count);
        for _ in 0..source_count {
            sources.push(reader.hash()?);
        }
        let provenance = Provenance {
            source: reader.string()?,
            license: reader.string()?,
            artifact: reader.optional_hash()?,
        };
        let normalization_count = reader.count("species_normalizations", NORMALIZATION_COUNT)?;
        let mut normalizations = Vec::with_capacity(normalization_count);
        for _ in 0..normalization_count {
            let target = match reader.u8()? {
                0 => SpeciesNormalizationTarget::MolarMass,
                1 => SpeciesNormalizationTarget::ReferencePressure,
                tag => return Err(reader.malformed(format!("unknown species target tag {tag}"))),
            };
            normalizations.push(SpeciesNormalizationReceipt::new(
                target,
                reader.hash()?,
                reader.dims()?,
                reader.f64()?,
                reader.f64()?,
                reader.string()?,
                reader.string()?,
            ));
        }
        reader.finish()?;
        let association = SpeciesAssociation::new(
            species,
            molar_mass,
            standard_state_phase,
            reference_eos,
            reference_pressure,
            elemental_reference,
            sources,
            provenance,
        )?;
        let pack = Self::new(
            pack_id,
            compiler,
            source_artifact,
            redistribution_terms,
            association,
            normalizations,
        )?;
        if pack.to_bytes() != bytes {
            return Err(PackError::Malformed {
                at: bytes.len(),
                detail: "decoded fields do not reproduce the canonical species-pack byte stream"
                    .to_string(),
            });
        }
        Ok(pack)
    }
}

fn validate_receipt(receipt: &SpeciesNormalizationReceipt) -> Result<(), PackError> {
    require_positive("species_normalization.scale", receipt.scale)?;
    if receipt.offset.to_bits() != 0 {
        return Err(invalid(
            "species_normalization.offset",
            "molar mass and pressure require a linear transform with canonical positive-zero offset",
        ));
    }
    require_text("species_normalization.source_basis", &receipt.source_basis)?;
    if receipt.target_basis != SPECIES_PACK_TARGET_BASIS {
        return Err(invalid(
            "species_normalization.target_basis",
            format!("target basis must be {SPECIES_PACK_TARGET_BASIS:?}"),
        ));
    }
    let expected_dims = match receipt.target {
        SpeciesNormalizationTarget::MolarMass => SPECIES_MOLAR_MASS_DIMS,
        SpeciesNormalizationTarget::ReferencePressure => SPECIES_REFERENCE_PRESSURE_DIMS,
    };
    if receipt.dims != expected_dims {
        return Err(invalid(
            "species_normalization.dims",
            format!(
                "target {:?} requires dimensions {expected_dims:?}, found {:?}",
                receipt.target, receipt.dims
            ),
        ));
    }
    Ok(())
}

fn require_elemental_reference(value: &str) -> Result<(), PackError> {
    if value.is_empty() {
        return Err(invalid(
            "species.elemental_reference",
            "elemental-reference id must not be empty",
        ));
    }
    if value.len() > MAX_ELEMENTAL_REFERENCE_BYTES {
        return Err(limit(
            "elemental_reference_bytes",
            MAX_ELEMENTAL_REFERENCE_BYTES,
            value.len(),
        ));
    }
    let bytes = value.as_bytes();
    if !bytes[0].is_ascii_alphanumeric()
        || !bytes.iter().all(|byte| {
            byte.is_ascii_alphanumeric() || matches!(byte, b'-' | b'_' | b'.' | b':' | b'/')
        })
    {
        return Err(invalid(
            "species.elemental_reference",
            "elemental-reference id must use compact ASCII without whitespace",
        ));
    }
    Ok(())
}

fn require_positive(field: &'static str, value: f64) -> Result<(), PackError> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(invalid(field, "value must be positive and finite"))
    }
}

fn require_text(field: &'static str, value: &str) -> Result<(), PackError> {
    if value.is_empty() {
        return Err(invalid(field, "value must not be empty"));
    }
    if value.len() > MAX_STRING_BYTES {
        return Err(limit("string_bytes", MAX_STRING_BYTES, value.len()));
    }
    if value.chars().any(char::is_control) {
        return Err(invalid(field, "value contains a control character"));
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

    fn u64(&mut self, value: u64) {
        self.bytes.extend_from_slice(&value.to_le_bytes());
    }

    fn count(&mut self, value: usize) {
        self.u32(u32::try_from(value).unwrap_or(u32::MAX));
    }

    fn string(&mut self, value: &str) {
        self.count(value.len());
        self.bytes.extend_from_slice(value.as_bytes());
    }

    fn hash(&mut self, value: ContentHash) {
        self.bytes.extend_from_slice(&value.0);
    }

    fn optional_hash(&mut self, value: Option<ContentHash>) {
        match value {
            None => self.u8(0),
            Some(value) => {
                self.u8(1);
                self.hash(value);
            }
        }
    }

    fn dims(&mut self, value: Dims) {
        for exponent in value.0 {
            self.u8(exponent.cast_unsigned());
        }
    }

    fn f64(&mut self, value: f64) {
        self.u64(value.to_bits());
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

    fn require_items(
        &self,
        count: usize,
        minimum_width: usize,
        field: &str,
    ) -> Result<(), PackError> {
        let minimum_bytes = count
            .checked_mul(minimum_width)
            .ok_or_else(|| self.malformed(format!("{field} byte length overflow")))?;
        if self.bytes.len().saturating_sub(self.cursor) < minimum_bytes {
            Err(self.malformed(format!(
                "truncated {field} needs at least {minimum_bytes} bytes"
            )))
        } else {
            Ok(())
        }
    }

    fn expect(&mut self, expected: &[u8], name: &str) -> Result<(), PackError> {
        if self.take(expected.len())? == expected {
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

    fn u64(&mut self) -> Result<u64, PackError> {
        let bytes: [u8; 8] = self
            .take(8)?
            .try_into()
            .map_err(|_| self.malformed("u64 width"))?;
        Ok(u64::from_le_bytes(bytes))
    }

    fn count(&mut self, resource: &'static str, maximum: usize) -> Result<usize, PackError> {
        let count = usize::try_from(self.u32()?)
            .map_err(|_| self.malformed(format!("{resource} count does not fit usize")))?;
        if count > maximum {
            Err(limit(resource, maximum, count))
        } else {
            Ok(count)
        }
    }

    fn string(&mut self) -> Result<String, PackError> {
        let length = self.count("string_bytes", MAX_STRING_BYTES)?;
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

    fn optional_hash(&mut self) -> Result<Option<ContentHash>, PackError> {
        match self.u8()? {
            0 => Ok(None),
            1 => self.hash().map(Some),
            tag => Err(self.malformed(format!("unknown optional-hash tag {tag}"))),
        }
    }

    fn dims(&mut self) -> Result<Dims, PackError> {
        let mut dims = [0i8; 6];
        for exponent in &mut dims {
            *exponent = i8::from_ne_bytes([self.u8()?]);
        }
        Ok(Dims(dims))
    }

    fn f64(&mut self) -> Result<f64, PackError> {
        Ok(f64::from_bits(self.u64()?))
    }

    fn finish(self) -> Result<(), PackError> {
        if self.cursor == self.bytes.len() {
            Ok(())
        } else {
            Err(PackError::Malformed {
                at: self.cursor,
                detail: format!(
                    "{} trailing bytes after canonical species pack",
                    self.bytes.len() - self.cursor
                ),
            })
        }
    }
}
