//! Bounded canonical transport primitives for V&V artifacts.

use core::fmt;

use fs_blake3::ContentHash;

use super::model::*;

const MAGIC: &[u8; 4] = b"FSVV";
const CANONICAL_RULE: &str = "vv-canonical-identity";
const ROOT_ARTIFACT: u8 = 0;
const ROOT_CASE: u8 = 1;

/// Maximum accepted size of one canonical V&V artifact or case transport.
pub const MAX_VV_CANONICAL_BYTES: usize = 4 * 1024 * 1024;
/// Maximum UTF-8 bytes accepted for one schema string.
pub(crate) const MAX_VV_STRING_BYTES: usize = MAX_VV_TEXT_BYTES;
/// Maximum entries accepted for one schema collection.
pub(crate) const MAX_VV_COLLECTION_ITEMS: usize = MAX_VV_MATRIX_DIMENSION * MAX_VV_MATRIX_DIMENSION;
/// Aggregate decoded entries across nested collections in one case.
pub(crate) const MAX_VV_TOTAL_COLLECTION_ITEMS: usize = 256 * 1024;

/// A bounded canonical-transport refusal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VvCodecError {
    offset: usize,
    detail: String,
}

impl VvCodecError {
    pub(crate) fn at(offset: usize, detail: impl Into<String>) -> Self {
        Self {
            offset,
            detail: detail.into(),
        }
    }

    /// Stable rule identifier for every wire-level refusal.
    #[must_use]
    pub const fn rule_name(&self) -> &'static str {
        CANONICAL_RULE
    }

    /// Byte offset at which decoding or encoding refused.
    #[must_use]
    pub const fn offset(&self) -> usize {
        self.offset
    }

    /// Human-readable detail that does not participate in rule matching.
    #[must_use]
    pub fn detail(&self) -> &str {
        &self.detail
    }
}

impl fmt::Display for VvCodecError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{CANONICAL_RULE} at byte {}: {}",
            self.offset, self.detail
        )
    }
}

impl std::error::Error for VvCodecError {}

/// Checked fixed-order encoder used by the schema model.
pub(crate) struct Encoder {
    bytes: Vec<u8>,
    collection_items: usize,
}

impl Encoder {
    pub(crate) fn new() -> Result<Self, VvCodecError> {
        let mut this = Self {
            bytes: Vec::new(),
            collection_items: 0,
        };
        this.raw(MAGIC)?;
        this.u32(VV_SCHEMA_VERSION)?;
        this.u32(VV_RULESET_VERSION)?;
        Ok(this)
    }

    fn reserve(&mut self, additional: usize) -> Result<(), VvCodecError> {
        let requested = self
            .bytes
            .len()
            .checked_add(additional)
            .ok_or_else(|| VvCodecError::at(self.bytes.len(), "encoded length overflow"))?;
        if requested > MAX_VV_CANONICAL_BYTES {
            return Err(VvCodecError::at(
                self.bytes.len(),
                format!(
                    "canonical transport would require {requested} bytes; limit is {MAX_VV_CANONICAL_BYTES}"
                ),
            ));
        }
        self.bytes
            .try_reserve(additional)
            .map_err(|_| VvCodecError::at(self.bytes.len(), "allocation refused"))
    }

    pub(crate) fn raw(&mut self, value: &[u8]) -> Result<(), VvCodecError> {
        self.reserve(value.len())?;
        self.bytes.extend_from_slice(value);
        Ok(())
    }

    pub(crate) fn u8(&mut self, value: u8) -> Result<(), VvCodecError> {
        self.raw(&[value])
    }

    pub(crate) fn bool(&mut self, value: bool) -> Result<(), VvCodecError> {
        self.u8(u8::from(value))
    }

    pub(crate) fn u32(&mut self, value: u32) -> Result<(), VvCodecError> {
        self.raw(&value.to_le_bytes())
    }

    pub(crate) fn u64(&mut self, value: u64) -> Result<(), VvCodecError> {
        self.raw(&value.to_le_bytes())
    }

    pub(crate) fn usize(&mut self, value: usize) -> Result<(), VvCodecError> {
        self.u64(u64::try_from(value).map_err(|_| {
            VvCodecError::at(self.bytes.len(), "collection length does not fit u64")
        })?)
    }

    pub(crate) fn f64(&mut self, value: f64) -> Result<(), VvCodecError> {
        self.u64(value.to_bits())
    }

    pub(crate) fn hash(&mut self, value: ContentHash) -> Result<(), VvCodecError> {
        self.raw(value.as_bytes())
    }

    pub(crate) fn string(&mut self, value: &str) -> Result<(), VvCodecError> {
        if value.len() > MAX_VV_STRING_BYTES {
            return Err(VvCodecError::at(
                self.bytes.len(),
                format!(
                    "string has {} bytes; per-string limit is {MAX_VV_STRING_BYTES}",
                    value.len()
                ),
            ));
        }
        self.usize(value.len())?;
        self.raw(value.as_bytes())
    }

    pub(crate) fn count(&mut self, count: usize) -> Result<(), VvCodecError> {
        if count > MAX_VV_COLLECTION_ITEMS {
            return Err(VvCodecError::at(
                self.bytes.len(),
                format!("collection has {count} entries; limit is {MAX_VV_COLLECTION_ITEMS}"),
            ));
        }
        self.collection_items = self
            .collection_items
            .checked_add(count)
            .ok_or_else(|| VvCodecError::at(self.bytes.len(), "aggregate item count overflow"))?;
        if self.collection_items > MAX_VV_TOTAL_COLLECTION_ITEMS {
            return Err(VvCodecError::at(
                self.bytes.len(),
                format!("aggregate collection entries exceed {MAX_VV_TOTAL_COLLECTION_ITEMS}"),
            ));
        }
        self.usize(count)
    }

    pub(crate) fn finish(self) -> Vec<u8> {
        self.bytes
    }
}

/// Checked fixed-order decoder used by the schema model.
pub(crate) struct Decoder<'a> {
    bytes: &'a [u8],
    offset: usize,
    collection_items: usize,
}

impl<'a> Decoder<'a> {
    pub(crate) fn new(bytes: &'a [u8]) -> Result<Self, VvCodecError> {
        if bytes.len() > MAX_VV_CANONICAL_BYTES {
            return Err(VvCodecError::at(
                0,
                format!(
                    "canonical transport has {} bytes; limit is {MAX_VV_CANONICAL_BYTES}",
                    bytes.len()
                ),
            ));
        }
        let mut this = Self {
            bytes,
            offset: 0,
            collection_items: 0,
        };
        this.exact(MAGIC)?;
        let schema = this.u32()?;
        if schema != VV_SCHEMA_VERSION {
            return Err(VvCodecError::at(
                4,
                format!("unsupported V&V schema version {schema}"),
            ));
        }
        let ruleset = this.u32()?;
        if ruleset != VV_RULESET_VERSION {
            return Err(VvCodecError::at(
                8,
                format!("unsupported V&V ruleset version {ruleset}"),
            ));
        }
        Ok(this)
    }

    pub(crate) fn position(&self) -> usize {
        self.offset
    }

    fn take(&mut self, count: usize) -> Result<&'a [u8], VvCodecError> {
        let end = self
            .offset
            .checked_add(count)
            .ok_or_else(|| VvCodecError::at(self.offset, "decode offset overflow"))?;
        let value = self
            .bytes
            .get(self.offset..end)
            .ok_or_else(|| VvCodecError::at(self.offset, "truncated canonical transport"))?;
        self.offset = end;
        Ok(value)
    }

    fn exact(&mut self, expected: &[u8]) -> Result<(), VvCodecError> {
        let start = self.offset;
        if self.take(expected.len())? == expected {
            Ok(())
        } else {
            Err(VvCodecError::at(start, "unexpected transport domain"))
        }
    }

    pub(crate) fn u8(&mut self) -> Result<u8, VvCodecError> {
        Ok(self.take(1)?[0])
    }

    pub(crate) fn bool(&mut self) -> Result<bool, VvCodecError> {
        let offset = self.offset;
        match self.u8()? {
            0 => Ok(false),
            1 => Ok(true),
            _ => Err(VvCodecError::at(offset, "invalid Boolean tag")),
        }
    }

    pub(crate) fn u32(&mut self) -> Result<u32, VvCodecError> {
        let bytes: [u8; 4] = self
            .take(4)?
            .try_into()
            .map_err(|_| VvCodecError::at(self.offset, "invalid u32 field"))?;
        Ok(u32::from_le_bytes(bytes))
    }

    pub(crate) fn u64(&mut self) -> Result<u64, VvCodecError> {
        let bytes: [u8; 8] = self
            .take(8)?
            .try_into()
            .map_err(|_| VvCodecError::at(self.offset, "invalid u64 field"))?;
        Ok(u64::from_le_bytes(bytes))
    }

    pub(crate) fn usize(&mut self) -> Result<usize, VvCodecError> {
        usize::try_from(self.u64()?)
            .map_err(|_| VvCodecError::at(self.offset, "length does not fit usize"))
    }

    pub(crate) fn f64(&mut self) -> Result<f64, VvCodecError> {
        Ok(f64::from_bits(self.u64()?))
    }

    pub(crate) fn hash(&mut self) -> Result<ContentHash, VvCodecError> {
        let offset = self.offset;
        ContentHash::from_slice(self.take(32)?)
            .ok_or_else(|| VvCodecError::at(offset, "invalid content hash"))
    }

    pub(crate) fn string(&mut self) -> Result<String, VvCodecError> {
        let offset = self.offset;
        let len = self.usize()?;
        if len > MAX_VV_STRING_BYTES {
            return Err(VvCodecError::at(
                offset,
                format!("string length {len} exceeds {MAX_VV_STRING_BYTES}"),
            ));
        }
        let value = self.take(len)?;
        let value = std::str::from_utf8(value)
            .map_err(|_| VvCodecError::at(offset, "string is not UTF-8"))?;
        let mut owned = String::new();
        owned
            .try_reserve_exact(value.len())
            .map_err(|_| VvCodecError::at(offset, "string allocation refused"))?;
        owned.push_str(value);
        Ok(owned)
    }

    pub(crate) fn count(&mut self) -> Result<usize, VvCodecError> {
        let offset = self.offset;
        let count = self.usize()?;
        if count > MAX_VV_COLLECTION_ITEMS {
            return Err(VvCodecError::at(
                offset,
                format!("collection length {count} exceeds {MAX_VV_COLLECTION_ITEMS}"),
            ));
        }
        self.collection_items = self
            .collection_items
            .checked_add(count)
            .ok_or_else(|| VvCodecError::at(offset, "aggregate item count overflow"))?;
        if self.collection_items > MAX_VV_TOTAL_COLLECTION_ITEMS {
            return Err(VvCodecError::at(
                offset,
                format!("aggregate collection entries exceed {MAX_VV_TOTAL_COLLECTION_ITEMS}"),
            ));
        }
        Ok(count)
    }

    pub(crate) fn finish(self) -> Result<(), VvCodecError> {
        if self.offset == self.bytes.len() {
            Ok(())
        } else {
            Err(VvCodecError::at(
                self.offset,
                format!(
                    "{} trailing byte(s) after the canonical value",
                    self.bytes.len() - self.offset
                ),
            ))
        }
    }
}

fn model_error(offset: usize, context: &str, error: VvErrors) -> VvCodecError {
    VvCodecError::at(offset, format!("{context} refused by the model: {error}"))
}

fn decode_model<T>(
    decoder: &Decoder<'_>,
    context: &str,
    result: Result<T, VvErrors>,
) -> Result<T, VvCodecError> {
    result.map_err(|error| model_error(decoder.position(), context, error))
}

fn bounded_vec<T>(
    decoder: &Decoder<'_>,
    count: usize,
    context: &str,
) -> Result<Vec<T>, VvCodecError> {
    let mut values = Vec::new();
    values.try_reserve_exact(count).map_err(|_| {
        VvCodecError::at(
            decoder.position(),
            format!("allocation refused for {context} with {count} entries"),
        )
    })?;
    Ok(values)
}

fn ensure_strictly_increasing<T: Ord>(
    previous: Option<&T>,
    current: &T,
    offset: usize,
    context: &str,
) -> Result<(), VvCodecError> {
    if previous.is_some_and(|value| value >= current) {
        Err(VvCodecError::at(
            offset,
            format!("{context} keys are duplicated or out of canonical order"),
        ))
    } else {
        Ok(())
    }
}

fn encode_id(encoder: &mut Encoder, value: &str) -> Result<(), VvCodecError> {
    encoder.string(value)
}

fn decode_id<T>(
    decoder: &mut Decoder<'_>,
    context: &str,
    constructor: impl FnOnce(String) -> Result<T, VvErrors>,
) -> Result<T, VvCodecError> {
    let offset = decoder.position();
    let value = decoder.string()?;
    constructor(value).map_err(|error| model_error(offset, context, error))
}

macro_rules! id_codec {
    ($encode:ident, $decode:ident, $type:ty, $context:literal) => {
        fn $encode(encoder: &mut Encoder, value: &$type) -> Result<(), VvCodecError> {
            encode_id(encoder, value.as_str())
        }

        fn $decode(decoder: &mut Decoder<'_>) -> Result<$type, VvCodecError> {
            decode_id(decoder, $context, <$type>::from_canonical)
        }
    };
}

id_codec!(
    encode_artifact_id,
    decode_artifact_id,
    ArtifactId,
    "artifact id"
);
id_codec!(encode_qoi_id, decode_qoi_id, QoiId, "QoI id");
id_codec!(
    encode_observation_id,
    decode_observation_id,
    ObservationId,
    "observation id"
);
id_codec!(
    encode_assumption_id,
    decode_assumption_id,
    AssumptionId,
    "assumption id"
);
id_codec!(encode_axis_id, decode_axis_id, AxisId, "axis id");
id_codec!(encode_unit_id, decode_unit_id, UnitId, "unit id");

fn encode_artifact_kind(encoder: &mut Encoder, value: ArtifactKind) -> Result<(), VvCodecError> {
    encoder.u8(value.canonical_wire_tag())
}

fn decode_artifact_kind(decoder: &mut Decoder<'_>) -> Result<ArtifactKind, VvCodecError> {
    let offset = decoder.position();
    match decoder.u8()? {
        0 => Ok(ArtifactKind::ContextOfUse),
        1 => Ok(ArtifactKind::ValidationPlan),
        2 => Ok(ArtifactKind::ExperimentArtifact),
        3 => Ok(ArtifactKind::CalibrationSplit),
        4 => Ok(ArtifactKind::SolutionVerificationReceipt),
        5 => Ok(ArtifactKind::PredictionAssessment),
        6 => Ok(ArtifactKind::AssumptionsLedger),
        tag => Err(VvCodecError::at(
            offset,
            format!("unknown artifact-kind tag {tag}"),
        )),
    }
}

fn encode_rule(encoder: &mut Encoder, value: VvRule) -> Result<(), VvCodecError> {
    encoder.u8(match value {
        VvRule::SchemaIdentity => 0,
        VvRule::SchemaCardinality => 1,
        VvRule::SplitPartitionsDisjoint => 2,
        VvRule::SplitBlindHoldoutSealed => 3,
        VvRule::ColorCategoricalOnly => 4,
        VvRule::ValidationRequiresPhysicalReferent => 5,
        VvRule::QoiDependencyClosed => 6,
        VvRule::QoiDependencyIsolated => 7,
        VvRule::WaterfallModeDeclared => 8,
        VvRule::WaterfallArithmetic => 9,
        VvRule::WaterfallDependenceDeclared => 10,
        VvRule::ExperimentInstrumentCalibration => 11,
        VvRule::ExperimentClockSynchronization => 12,
        VvRule::ExperimentRepeatabilityCovariance => 13,
        VvRule::ExperimentDataAuthenticity => 14,
        VvRule::DiagnosticObservability => 15,
        VvRule::DiagnosticIdentifiability => 16,
        VvRule::DiagnosticConfounding => 17,
        VvRule::DiagnosticInverseCrime => 18,
        VvRule::ValidationMetricUncertainty => 19,
        VvRule::SolutionVerificationComplete => 20,
        VvRule::ApplicabilityDecision => 21,
        VvRule::ApplicabilityPolicy => 22,
        VvRule::ProcessConformanceSeparate => 23,
        VvRule::AssumptionRowComplete => 24,
        VvRule::AssumptionDomainEnforced => 25,
        VvRule::AssumptionA001 => 26,
        VvRule::AssumptionA002 => 27,
        VvRule::AssumptionA003 => 28,
        VvRule::AssumptionA004 => 29,
        VvRule::AssumptionA005 => 30,
        VvRule::AssumptionA006 => 31,
        VvRule::AssumptionA007 => 32,
        VvRule::AssumptionA008 => 33,
        VvRule::ReceiptBinding => 34,
    })
}

fn decode_rule(decoder: &mut Decoder<'_>) -> Result<VvRule, VvCodecError> {
    let offset = decoder.position();
    match decoder.u8()? {
        0 => Ok(VvRule::SchemaIdentity),
        1 => Ok(VvRule::SchemaCardinality),
        2 => Ok(VvRule::SplitPartitionsDisjoint),
        3 => Ok(VvRule::SplitBlindHoldoutSealed),
        4 => Ok(VvRule::ColorCategoricalOnly),
        5 => Ok(VvRule::ValidationRequiresPhysicalReferent),
        6 => Ok(VvRule::QoiDependencyClosed),
        7 => Ok(VvRule::QoiDependencyIsolated),
        8 => Ok(VvRule::WaterfallModeDeclared),
        9 => Ok(VvRule::WaterfallArithmetic),
        10 => Ok(VvRule::WaterfallDependenceDeclared),
        11 => Ok(VvRule::ExperimentInstrumentCalibration),
        12 => Ok(VvRule::ExperimentClockSynchronization),
        13 => Ok(VvRule::ExperimentRepeatabilityCovariance),
        14 => Ok(VvRule::ExperimentDataAuthenticity),
        15 => Ok(VvRule::DiagnosticObservability),
        16 => Ok(VvRule::DiagnosticIdentifiability),
        17 => Ok(VvRule::DiagnosticConfounding),
        18 => Ok(VvRule::DiagnosticInverseCrime),
        19 => Ok(VvRule::ValidationMetricUncertainty),
        20 => Ok(VvRule::SolutionVerificationComplete),
        21 => Ok(VvRule::ApplicabilityDecision),
        22 => Ok(VvRule::ApplicabilityPolicy),
        23 => Ok(VvRule::ProcessConformanceSeparate),
        24 => Ok(VvRule::AssumptionRowComplete),
        25 => Ok(VvRule::AssumptionDomainEnforced),
        26 => Ok(VvRule::AssumptionA001),
        27 => Ok(VvRule::AssumptionA002),
        28 => Ok(VvRule::AssumptionA003),
        29 => Ok(VvRule::AssumptionA004),
        30 => Ok(VvRule::AssumptionA005),
        31 => Ok(VvRule::AssumptionA006),
        32 => Ok(VvRule::AssumptionA007),
        33 => Ok(VvRule::AssumptionA008),
        34 => Ok(VvRule::ReceiptBinding),
        tag => Err(VvCodecError::at(
            offset,
            format!("unknown V&V-rule tag {tag}"),
        )),
    }
}

fn encode_artifact_ref(encoder: &mut Encoder, value: &ArtifactRef) -> Result<(), VvCodecError> {
    encode_artifact_kind(encoder, value.kind())?;
    encode_artifact_id(encoder, value.id())?;
    encoder.hash(value.hash())
}

fn decode_artifact_ref(decoder: &mut Decoder<'_>) -> Result<ArtifactRef, VvCodecError> {
    Ok(ArtifactRef::new(
        decode_artifact_kind(decoder)?,
        decode_artifact_id(decoder)?,
        decoder.hash()?,
    ))
}

fn encode_seed(encoder: &mut Encoder, value: &SeedDeclaration) -> Result<(), VvCodecError> {
    match value {
        SeedDeclaration::Fixed(seed) => {
            encoder.u8(0)?;
            encoder.u64(*seed)
        }
        SeedDeclaration::NotApplicable { reason } => {
            encoder.u8(1)?;
            encoder.string(reason)
        }
    }
}

fn decode_seed(decoder: &mut Decoder<'_>) -> Result<SeedDeclaration, VvCodecError> {
    let offset = decoder.position();
    match decoder.u8()? {
        0 => Ok(SeedDeclaration::Fixed(decoder.u64()?)),
        1 => Ok(SeedDeclaration::NotApplicable {
            reason: decoder.string()?,
        }),
        tag => Err(VvCodecError::at(
            offset,
            format!("unknown seed-declaration tag {tag}"),
        )),
    }
}

fn encode_f64_budget(
    encoder: &mut Encoder,
    value: &DeclaredBudget<f64>,
) -> Result<(), VvCodecError> {
    match value {
        DeclaredBudget::Limit(limit) => {
            encoder.u8(0)?;
            encoder.f64(*limit)
        }
        DeclaredBudget::NotApplicable { reason } => {
            encoder.u8(1)?;
            encoder.string(reason)
        }
    }
}

fn decode_f64_budget(decoder: &mut Decoder<'_>) -> Result<DeclaredBudget<f64>, VvCodecError> {
    let offset = decoder.position();
    match decoder.u8()? {
        0 => Ok(DeclaredBudget::Limit(decoder.f64()?)),
        1 => Ok(DeclaredBudget::NotApplicable {
            reason: decoder.string()?,
        }),
        tag => Err(VvCodecError::at(
            offset,
            format!("unknown floating-point budget tag {tag}"),
        )),
    }
}

fn encode_u64_budget(
    encoder: &mut Encoder,
    value: &DeclaredBudget<u64>,
) -> Result<(), VvCodecError> {
    match value {
        DeclaredBudget::Limit(limit) => {
            encoder.u8(0)?;
            encoder.u64(*limit)
        }
        DeclaredBudget::NotApplicable { reason } => {
            encoder.u8(1)?;
            encoder.string(reason)
        }
    }
}

fn decode_u64_budget(decoder: &mut Decoder<'_>) -> Result<DeclaredBudget<u64>, VvCodecError> {
    let offset = decoder.position();
    match decoder.u8()? {
        0 => Ok(DeclaredBudget::Limit(decoder.u64()?)),
        1 => Ok(DeclaredBudget::NotApplicable {
            reason: decoder.string()?,
        }),
        tag => Err(VvCodecError::at(
            offset,
            format!("unknown integer budget tag {tag}"),
        )),
    }
}

fn encode_header(encoder: &mut Encoder, value: &ArtifactHeader) -> Result<(), VvCodecError> {
    encode_artifact_id(encoder, value.id())?;
    encoder.count(value.units().len())?;
    for unit in value.units() {
        encode_unit_id(encoder, unit)?;
    }
    encode_seed(encoder, value.seed())?;
    encode_f64_budget(encoder, value.accuracy())?;
    encode_u64_budget(encoder, value.time_ms())?;
    encode_u64_budget(encoder, value.memory_bytes())?;
    encoder.count(value.versions().len())?;
    for (component, version) in value.versions() {
        encoder.string(component)?;
        encoder.string(version)?;
    }
    encoder.count(value.capabilities().len())?;
    for capability in value.capabilities() {
        encoder.string(capability)?;
    }
    Ok(())
}

fn decode_header(decoder: &mut Decoder<'_>) -> Result<ArtifactHeader, VvCodecError> {
    let id = decode_artifact_id(decoder)?;

    let count = decoder.count()?;
    let mut units = bounded_vec(decoder, count, "header units")?;
    for _ in 0..count {
        let offset = decoder.position();
        let unit = decode_unit_id(decoder)?;
        ensure_strictly_increasing(units.last(), &unit, offset, "header unit")?;
        units.push(unit);
    }

    let seed = decode_seed(decoder)?;
    let accuracy = decode_f64_budget(decoder)?;
    let time_ms = decode_u64_budget(decoder)?;
    let memory_bytes = decode_u64_budget(decoder)?;

    let count = decoder.count()?;
    let mut versions = bounded_vec(decoder, count, "header versions")?;
    for _ in 0..count {
        let offset = decoder.position();
        let component = decoder.string()?;
        ensure_strictly_increasing(
            versions.last().map(|(key, _)| key),
            &component,
            offset,
            "header version",
        )?;
        versions.push((component, decoder.string()?));
    }

    let count = decoder.count()?;
    let mut capabilities = bounded_vec(decoder, count, "header capabilities")?;
    for _ in 0..count {
        let offset = decoder.position();
        let capability = decoder.string()?;
        ensure_strictly_increasing(
            capabilities.last(),
            &capability,
            offset,
            "header capability",
        )?;
        capabilities.push(capability);
    }

    decode_model(
        decoder,
        "artifact header",
        ArtifactHeader::try_new(
            id,
            units,
            seed,
            accuracy,
            time_ms,
            memory_bytes,
            versions,
            capabilities,
        ),
    )
}

fn encode_numeric_domain_axis(
    encoder: &mut Encoder,
    value: &NumericDomainAxis,
) -> Result<(), VvCodecError> {
    encode_axis_id(encoder, value.axis())?;
    encode_unit_id(encoder, value.unit())?;
    let (lo, hi) = value.bounds();
    encoder.f64(lo)?;
    encoder.f64(hi)
}

fn decode_numeric_domain_axis(
    decoder: &mut Decoder<'_>,
) -> Result<NumericDomainAxis, VvCodecError> {
    let axis = decode_axis_id(decoder)?;
    let unit = decode_unit_id(decoder)?;
    let lo = decoder.f64()?;
    let hi = decoder.f64()?;
    decode_model(
        decoder,
        "numeric applicability axis",
        NumericDomainAxis::try_new(axis, unit, lo, hi),
    )
}

fn encode_categorical_domain_axis(
    encoder: &mut Encoder,
    value: &CategoricalDomainAxis,
) -> Result<(), VvCodecError> {
    encode_axis_id(encoder, value.axis())?;
    encoder.count(value.allowed().len())?;
    for allowed in value.allowed() {
        encoder.string(allowed)?;
    }
    Ok(())
}

fn decode_categorical_domain_axis(
    decoder: &mut Decoder<'_>,
) -> Result<CategoricalDomainAxis, VvCodecError> {
    let axis = decode_axis_id(decoder)?;
    let count = decoder.count()?;
    let mut allowed = bounded_vec(decoder, count, "categorical-domain values")?;
    for _ in 0..count {
        let offset = decoder.position();
        let value = decoder.string()?;
        ensure_strictly_increasing(allowed.last(), &value, offset, "categorical-domain value")?;
        allowed.push(value);
    }
    decode_model(
        decoder,
        "categorical applicability axis",
        CategoricalDomainAxis::try_new(axis, allowed),
    )
}

fn encode_applicability_domain(
    encoder: &mut Encoder,
    value: &ApplicabilityDomain,
) -> Result<(), VvCodecError> {
    encoder.count(value.numeric().len())?;
    for row in value.numeric().values() {
        encode_numeric_domain_axis(encoder, row)?;
    }
    encoder.count(value.categorical().len())?;
    for row in value.categorical().values() {
        encode_categorical_domain_axis(encoder, row)?;
    }
    Ok(())
}

fn decode_applicability_domain(
    decoder: &mut Decoder<'_>,
) -> Result<ApplicabilityDomain, VvCodecError> {
    let count = decoder.count()?;
    let mut numeric = bounded_vec(decoder, count, "numeric applicability axes")?;
    for _ in 0..count {
        let offset = decoder.position();
        let row = decode_numeric_domain_axis(decoder)?;
        ensure_strictly_increasing(
            numeric.last().map(NumericDomainAxis::axis),
            row.axis(),
            offset,
            "numeric applicability axis",
        )?;
        numeric.push(row);
    }
    let count = decoder.count()?;
    let mut categorical = bounded_vec(decoder, count, "categorical applicability axes")?;
    for _ in 0..count {
        let offset = decoder.position();
        let row = decode_categorical_domain_axis(decoder)?;
        ensure_strictly_increasing(
            categorical.last().map(CategoricalDomainAxis::axis),
            row.axis(),
            offset,
            "categorical applicability axis",
        )?;
        categorical.push(row);
    }
    decode_model(
        decoder,
        "applicability domain",
        ApplicabilityDomain::try_new(numeric, categorical),
    )
}

fn encode_applicability_point(
    encoder: &mut Encoder,
    value: &ApplicabilityPoint,
) -> Result<(), VvCodecError> {
    encoder.count(value.numeric().len())?;
    for (axis, coordinate) in value.numeric() {
        encode_axis_id(encoder, axis)?;
        encoder.f64(*coordinate)?;
    }
    encoder.count(value.categorical().len())?;
    for (axis, coordinate) in value.categorical() {
        encode_axis_id(encoder, axis)?;
        encoder.string(coordinate)?;
    }
    Ok(())
}

fn decode_applicability_point(
    decoder: &mut Decoder<'_>,
) -> Result<ApplicabilityPoint, VvCodecError> {
    let count = decoder.count()?;
    let mut numeric = bounded_vec(decoder, count, "numeric applicability point")?;
    for _ in 0..count {
        let offset = decoder.position();
        let axis = decode_axis_id(decoder)?;
        ensure_strictly_increasing(
            numeric.last().map(|(axis, _)| axis),
            &axis,
            offset,
            "numeric applicability-point axis",
        )?;
        numeric.push((axis, decoder.f64()?));
    }
    let count = decoder.count()?;
    let mut categorical = bounded_vec(decoder, count, "categorical applicability point")?;
    for _ in 0..count {
        let offset = decoder.position();
        let axis = decode_axis_id(decoder)?;
        ensure_strictly_increasing(
            categorical.last().map(|(axis, _)| axis),
            &axis,
            offset,
            "categorical applicability-point axis",
        )?;
        categorical.push((axis, decoder.string()?));
    }
    decode_model(
        decoder,
        "applicability point",
        ApplicabilityPoint::try_new(numeric, categorical),
    )
}

fn encode_domain_violation(
    encoder: &mut Encoder,
    value: &DomainViolation,
) -> Result<(), VvCodecError> {
    match value {
        DomainViolation::Missing { axis } => {
            encoder.u8(0)?;
            encode_axis_id(encoder, axis)
        }
        DomainViolation::Numeric {
            axis,
            value,
            lo,
            hi,
        } => {
            encoder.u8(1)?;
            encode_axis_id(encoder, axis)?;
            encoder.f64(*value)?;
            encoder.f64(*lo)?;
            encoder.f64(*hi)
        }
        DomainViolation::Categorical { axis, value } => {
            encoder.u8(2)?;
            encode_axis_id(encoder, axis)?;
            encoder.string(value)
        }
        DomainViolation::Assumption { id } => {
            encoder.u8(3)?;
            encode_assumption_id(encoder, id)
        }
    }
}

fn decode_domain_violation(decoder: &mut Decoder<'_>) -> Result<DomainViolation, VvCodecError> {
    let offset = decoder.position();
    match decoder.u8()? {
        0 => Ok(DomainViolation::Missing {
            axis: decode_axis_id(decoder)?,
        }),
        1 => Ok(DomainViolation::Numeric {
            axis: decode_axis_id(decoder)?,
            value: decoder.f64()?,
            lo: decoder.f64()?,
            hi: decoder.f64()?,
        }),
        2 => Ok(DomainViolation::Categorical {
            axis: decode_axis_id(decoder)?,
            value: decoder.string()?,
        }),
        3 => Ok(DomainViolation::Assumption {
            id: decode_assumption_id(decoder)?,
        }),
        tag => Err(VvCodecError::at(
            offset,
            format!("unknown domain-violation tag {tag}"),
        )),
    }
}

fn encode_applicability_policy(
    encoder: &mut Encoder,
    value: ApplicabilityPolicy,
) -> Result<(), VvCodecError> {
    encoder.u8(match value {
        ApplicabilityPolicy::Demote => 0,
        ApplicabilityPolicy::Refuse => 1,
    })
}

fn decode_applicability_policy(
    decoder: &mut Decoder<'_>,
) -> Result<ApplicabilityPolicy, VvCodecError> {
    let offset = decoder.position();
    match decoder.u8()? {
        0 => Ok(ApplicabilityPolicy::Demote),
        1 => Ok(ApplicabilityPolicy::Refuse),
        tag => Err(VvCodecError::at(
            offset,
            format!("unknown applicability-policy tag {tag}"),
        )),
    }
}

fn encode_violations(
    encoder: &mut Encoder,
    violations: &[DomainViolation],
) -> Result<(), VvCodecError> {
    encoder.count(violations.len())?;
    for violation in violations {
        encode_domain_violation(encoder, violation)?;
    }
    Ok(())
}

fn decode_violations(decoder: &mut Decoder<'_>) -> Result<Vec<DomainViolation>, VvCodecError> {
    let count = decoder.count()?;
    let mut violations = bounded_vec(decoder, count, "domain violations")?;
    for _ in 0..count {
        violations.push(decode_domain_violation(decoder)?);
    }
    Ok(violations)
}

fn encode_applicability_decision(
    encoder: &mut Encoder,
    value: &ApplicabilityDecision,
) -> Result<(), VvCodecError> {
    match value {
        ApplicabilityDecision::InDomain => encoder.u8(0),
        ApplicabilityDecision::Demoted { violations } => {
            encoder.u8(1)?;
            encode_violations(encoder, violations)
        }
        ApplicabilityDecision::Refused { violations } => {
            encoder.u8(2)?;
            encode_violations(encoder, violations)
        }
    }
}

fn decode_applicability_decision(
    decoder: &mut Decoder<'_>,
) -> Result<ApplicabilityDecision, VvCodecError> {
    let offset = decoder.position();
    match decoder.u8()? {
        0 => Ok(ApplicabilityDecision::InDomain),
        1 => Ok(ApplicabilityDecision::Demoted {
            violations: decode_violations(decoder)?,
        }),
        2 => Ok(ApplicabilityDecision::Refused {
            violations: decode_violations(decoder)?,
        }),
        tag => Err(VvCodecError::at(
            offset,
            format!("unknown applicability-decision tag {tag}"),
        )),
    }
}

fn encode_acceptance(
    encoder: &mut Encoder,
    value: &AcceptanceCriterion,
) -> Result<(), VvCodecError> {
    match value {
        AcceptanceCriterion::ClosedRange { lo, hi } => {
            encoder.u8(0)?;
            encoder.f64(*lo)?;
            encoder.f64(*hi)
        }
        AcceptanceCriterion::AbsoluteErrorAtMost { limit } => {
            encoder.u8(1)?;
            encoder.f64(*limit)
        }
        AcceptanceCriterion::RelativeErrorAtMost { limit } => {
            encoder.u8(2)?;
            encoder.f64(*limit)
        }
        AcceptanceCriterion::CategoryEquals { expected } => {
            encoder.u8(3)?;
            encoder.string(expected)
        }
    }
}

fn decode_acceptance(decoder: &mut Decoder<'_>) -> Result<AcceptanceCriterion, VvCodecError> {
    let offset = decoder.position();
    match decoder.u8()? {
        0 => Ok(AcceptanceCriterion::ClosedRange {
            lo: decoder.f64()?,
            hi: decoder.f64()?,
        }),
        1 => Ok(AcceptanceCriterion::AbsoluteErrorAtMost {
            limit: decoder.f64()?,
        }),
        2 => Ok(AcceptanceCriterion::RelativeErrorAtMost {
            limit: decoder.f64()?,
        }),
        3 => Ok(AcceptanceCriterion::CategoryEquals {
            expected: decoder.string()?,
        }),
        tag => Err(VvCodecError::at(
            offset,
            format!("unknown acceptance-criterion tag {tag}"),
        )),
    }
}

fn encode_qoi_spec(encoder: &mut Encoder, value: &QoiSpec) -> Result<(), VvCodecError> {
    encode_qoi_id(encoder, value.id())?;
    encoder.string(value.name())?;
    encode_unit_id(encoder, value.unit())?;
    encode_acceptance(encoder, value.acceptance())
}

fn decode_qoi_spec(decoder: &mut Decoder<'_>) -> Result<QoiSpec, VvCodecError> {
    let id = decode_qoi_id(decoder)?;
    let name = decoder.string()?;
    let unit = decode_unit_id(decoder)?;
    let acceptance = decode_acceptance(decoder)?;
    decode_model(
        decoder,
        "QoI specification",
        QoiSpec::try_new(id, name, unit, acceptance),
    )
}

fn encode_context(encoder: &mut Encoder, value: &ContextOfUse) -> Result<(), VvCodecError> {
    encode_header(encoder, value.header())?;
    encoder.string(value.decision())?;
    encoder.count(value.qois().len())?;
    for qoi in value.qois().values() {
        encode_qoi_spec(encoder, qoi)?;
    }
    encode_applicability_domain(encoder, value.applicability())?;
    encode_applicability_policy(encoder, value.applicability_policy())
}

fn decode_context(decoder: &mut Decoder<'_>) -> Result<ContextOfUse, VvCodecError> {
    let header = decode_header(decoder)?;
    let decision = decoder.string()?;
    let count = decoder.count()?;
    let mut qois = bounded_vec(decoder, count, "context QoIs")?;
    for _ in 0..count {
        let offset = decoder.position();
        let qoi = decode_qoi_spec(decoder)?;
        ensure_strictly_increasing(
            qois.last().map(QoiSpec::id),
            qoi.id(),
            offset,
            "context QoI",
        )?;
        qois.push(qoi);
    }
    let applicability = decode_applicability_domain(decoder)?;
    let policy = decode_applicability_policy(decoder)?;
    decode_model(
        decoder,
        "context of use",
        ContextOfUse::try_new(header, decision, qois, applicability, policy),
    )
}

fn encode_diagnostic_record(
    encoder: &mut Encoder,
    value: &DiagnosticRecord,
) -> Result<(), VvCodecError> {
    encoder.bool(value.passed())?;
    encoder.hash(value.artifact_hash())?;
    encoder.string(value.detail())
}

fn decode_diagnostic_record(decoder: &mut Decoder<'_>) -> Result<DiagnosticRecord, VvCodecError> {
    let passed = decoder.bool()?;
    let artifact_hash = decoder.hash()?;
    let detail = decoder.string()?;
    decode_model(
        decoder,
        "diagnostic record",
        DiagnosticRecord::try_new(passed, artifact_hash, detail),
    )
}

fn encode_diagnostic_plan(
    encoder: &mut Encoder,
    value: &DiagnosticPlan,
) -> Result<(), VvCodecError> {
    encode_diagnostic_record(encoder, value.observability())?;
    encode_diagnostic_record(encoder, value.identifiability())?;
    encode_diagnostic_record(encoder, value.confounding())?;
    encode_diagnostic_record(encoder, value.inverse_crime())
}

fn decode_diagnostic_plan(decoder: &mut Decoder<'_>) -> Result<DiagnosticPlan, VvCodecError> {
    Ok(DiagnosticPlan::new(
        decode_diagnostic_record(decoder)?,
        decode_diagnostic_record(decoder)?,
        decode_diagnostic_record(decoder)?,
        decode_diagnostic_record(decoder)?,
    ))
}

fn encode_validation_metric_spec(
    encoder: &mut Encoder,
    value: &ValidationMetricSpec,
) -> Result<(), VvCodecError> {
    match value {
        ValidationMetricSpec::IntervalAgreement => encoder.u8(0),
        ValidationMetricSpec::NormalizedDiscrepancy { maximum } => {
            encoder.u8(1)?;
            encoder.f64(*maximum)
        }
        ValidationMetricSpec::PosteriorPredictive {
            minimum_tail_probability,
        } => {
            encoder.u8(2)?;
            encoder.f64(*minimum_tail_probability)
        }
    }
}

fn decode_validation_metric_spec(
    decoder: &mut Decoder<'_>,
) -> Result<ValidationMetricSpec, VvCodecError> {
    let offset = decoder.position();
    match decoder.u8()? {
        0 => Ok(ValidationMetricSpec::IntervalAgreement),
        1 => Ok(ValidationMetricSpec::NormalizedDiscrepancy {
            maximum: decoder.f64()?,
        }),
        2 => Ok(ValidationMetricSpec::PosteriorPredictive {
            minimum_tail_probability: decoder.f64()?,
        }),
        tag => Err(VvCodecError::at(
            offset,
            format!("unknown validation-metric specification tag {tag}"),
        )),
    }
}

fn encode_qoi_validation_plan(
    encoder: &mut Encoder,
    value: &QoiValidationPlan,
) -> Result<(), VvCodecError> {
    encode_qoi_id(encoder, value.qoi())?;
    encoder.count(value.experiments().len())?;
    for experiment in value.experiments() {
        encode_artifact_ref(encoder, experiment)?;
    }
    encode_artifact_ref(encoder, value.split())?;
    encoder.count(value.metrics().len())?;
    for metric in value.metrics() {
        encode_validation_metric_spec(encoder, metric)?;
    }
    encode_diagnostic_plan(encoder, value.diagnostics())
}

fn decode_qoi_validation_plan(
    decoder: &mut Decoder<'_>,
) -> Result<QoiValidationPlan, VvCodecError> {
    let qoi = decode_qoi_id(decoder)?;
    let count = decoder.count()?;
    let mut experiments = bounded_vec(decoder, count, "validation experiment references")?;
    for _ in 0..count {
        let offset = decoder.position();
        let experiment = decode_artifact_ref(decoder)?;
        ensure_strictly_increasing(
            experiments.last(),
            &experiment,
            offset,
            "validation experiment reference",
        )?;
        experiments.push(experiment);
    }
    let split = decode_artifact_ref(decoder)?;
    let count = decoder.count()?;
    let mut metrics: Vec<ValidationMetricSpec> =
        bounded_vec(decoder, count, "validation metric specifications")?;
    for _ in 0..count {
        let offset = decoder.position();
        let metric = decode_validation_metric_spec(decoder)?;
        if metrics
            .last()
            .is_some_and(|previous| previous.canonical_key() >= metric.canonical_key())
        {
            return Err(VvCodecError::at(
                offset,
                "validation metric specifications are duplicated or out of canonical order",
            ));
        }
        metrics.push(metric);
    }
    let diagnostics = decode_diagnostic_plan(decoder)?;
    decode_model(
        decoder,
        "QoI validation-plan row",
        QoiValidationPlan::try_new(qoi, experiments, split, metrics, diagnostics),
    )
}

fn encode_validation_plan(
    encoder: &mut Encoder,
    value: &ValidationPlan,
) -> Result<(), VvCodecError> {
    encode_header(encoder, value.header())?;
    encode_artifact_ref(encoder, value.context())?;
    encoder.count(value.by_qoi().len())?;
    for row in value.by_qoi().values() {
        encode_qoi_validation_plan(encoder, row)?;
    }
    Ok(())
}

fn decode_validation_plan(decoder: &mut Decoder<'_>) -> Result<ValidationPlan, VvCodecError> {
    let header = decode_header(decoder)?;
    let context = decode_artifact_ref(decoder)?;
    let count = decoder.count()?;
    let mut rows = bounded_vec(decoder, count, "validation-plan QoIs")?;
    for _ in 0..count {
        let offset = decoder.position();
        let row = decode_qoi_validation_plan(decoder)?;
        ensure_strictly_increasing(
            rows.last().map(QoiValidationPlan::qoi),
            row.qoi(),
            offset,
            "validation-plan QoI",
        )?;
        rows.push(row);
    }
    decode_model(
        decoder,
        "validation plan",
        ValidationPlan::try_new(header, context, rows),
    )
}

fn encode_experiment_origin(
    encoder: &mut Encoder,
    value: &ExperimentOrigin,
) -> Result<(), VvCodecError> {
    match value {
        ExperimentOrigin::Physical {
            apparatus_id,
            facility_id,
        } => {
            encoder.u8(0)?;
            encode_artifact_id(encoder, apparatus_id)?;
            encode_artifact_id(encoder, facility_id)
        }
        ExperimentOrigin::SyntheticHighFidelity { producer } => {
            encoder.u8(1)?;
            encode_artifact_id(encoder, producer)
        }
        ExperimentOrigin::SecondImplementation { producer } => {
            encoder.u8(2)?;
            encode_artifact_id(encoder, producer)
        }
    }
}

fn decode_experiment_origin(decoder: &mut Decoder<'_>) -> Result<ExperimentOrigin, VvCodecError> {
    let offset = decoder.position();
    match decoder.u8()? {
        0 => Ok(ExperimentOrigin::Physical {
            apparatus_id: decode_artifact_id(decoder)?,
            facility_id: decode_artifact_id(decoder)?,
        }),
        1 => Ok(ExperimentOrigin::SyntheticHighFidelity {
            producer: decode_artifact_id(decoder)?,
        }),
        2 => Ok(ExperimentOrigin::SecondImplementation {
            producer: decode_artifact_id(decoder)?,
        }),
        tag => Err(VvCodecError::at(
            offset,
            format!("unknown experiment-origin tag {tag}"),
        )),
    }
}

fn encode_instrument_calibration(
    encoder: &mut Encoder,
    value: &InstrumentCalibration,
) -> Result<(), VvCodecError> {
    encode_artifact_id(encoder, value.instrument_id())?;
    encoder.hash(value.certificate_hash())?;
    encoder.bool(value.current())
}

fn decode_instrument_calibration(
    decoder: &mut Decoder<'_>,
) -> Result<InstrumentCalibration, VvCodecError> {
    Ok(InstrumentCalibration::new(
        decode_artifact_id(decoder)?,
        decoder.hash()?,
        decoder.bool()?,
    ))
}

fn encode_clock_synchronization(
    encoder: &mut Encoder,
    value: &ClockSynchronization,
) -> Result<(), VvCodecError> {
    match value {
        ClockSynchronization::SingleClock { clock_id } => {
            encoder.u8(0)?;
            encode_artifact_id(encoder, clock_id)
        }
        ClockSynchronization::Synchronized {
            clock_ids,
            method,
            max_skew_seconds,
            evidence_hash,
        } => {
            encoder.u8(1)?;
            encoder.count(clock_ids.len())?;
            for clock_id in clock_ids {
                encode_artifact_id(encoder, clock_id)?;
            }
            encoder.string(method)?;
            encoder.f64(*max_skew_seconds)?;
            encoder.hash(*evidence_hash)
        }
    }
}

fn decode_clock_synchronization(
    decoder: &mut Decoder<'_>,
) -> Result<ClockSynchronization, VvCodecError> {
    let offset = decoder.position();
    match decoder.u8()? {
        0 => Ok(ClockSynchronization::SingleClock {
            clock_id: decode_artifact_id(decoder)?,
        }),
        1 => {
            let count = decoder.count()?;
            let mut clock_ids = bounded_vec(decoder, count, "synchronized clock ids")?;
            for _ in 0..count {
                let offset = decoder.position();
                let clock_id = decode_artifact_id(decoder)?;
                ensure_strictly_increasing(
                    clock_ids.last(),
                    &clock_id,
                    offset,
                    "synchronized clock id",
                )?;
                clock_ids.push(clock_id);
            }
            let method = decoder.string()?;
            let max_skew_seconds = decoder.f64()?;
            let evidence_hash = decoder.hash()?;
            decode_model(
                decoder,
                "clock synchronization",
                ClockSynchronization::synchronized(
                    clock_ids,
                    method,
                    max_skew_seconds,
                    evidence_hash,
                ),
            )
        }
        tag => Err(VvCodecError::at(
            offset,
            format!("unknown clock-synchronization tag {tag}"),
        )),
    }
}

fn encode_covariance_matrix(
    encoder: &mut Encoder,
    value: &CovarianceMatrix,
) -> Result<(), VvCodecError> {
    encoder.usize(value.dimension())?;
    encoder.count(value.lower_triangle().len())?;
    for entry in value.lower_triangle() {
        encoder.f64(*entry)?;
    }
    Ok(())
}

fn decode_covariance_matrix(decoder: &mut Decoder<'_>) -> Result<CovarianceMatrix, VvCodecError> {
    let dimension = decoder.usize()?;
    let count = decoder.count()?;
    let mut lower_triangle = bounded_vec(decoder, count, "covariance entries")?;
    for _ in 0..count {
        lower_triangle.push(decoder.f64()?);
    }
    decode_model(
        decoder,
        "covariance matrix",
        CovarianceMatrix::try_new(dimension, lower_triangle),
    )
}

fn encode_repeatability(
    encoder: &mut Encoder,
    value: &RepeatabilitySummary,
) -> Result<(), VvCodecError> {
    encoder.count(value.qoi_order().len())?;
    for qoi in value.qoi_order() {
        encode_qoi_id(encoder, qoi)?;
    }
    encoder.u32(value.replicates())?;
    encode_covariance_matrix(encoder, value.covariance())
}

fn decode_repeatability(decoder: &mut Decoder<'_>) -> Result<RepeatabilitySummary, VvCodecError> {
    let count = decoder.count()?;
    let mut qoi_order = bounded_vec(decoder, count, "repeatability covariance QoI axes")?;
    for _ in 0..count {
        qoi_order.push(decode_qoi_id(decoder)?);
    }
    let replicates = decoder.u32()?;
    let covariance = decode_covariance_matrix(decoder)?;
    decode_model(
        decoder,
        "repeatability summary",
        RepeatabilitySummary::try_new_for_qois(replicates, qoi_order, covariance),
    )
}

fn encode_data_authenticity(
    encoder: &mut Encoder,
    value: &DataAuthenticity,
) -> Result<(), VvCodecError> {
    encoder.hash(value.source_bytes_hash())?;
    encoder.hash(value.custody_receipt_hash())?;
    encoder.bool(value.authenticated())
}

fn decode_data_authenticity(decoder: &mut Decoder<'_>) -> Result<DataAuthenticity, VvCodecError> {
    Ok(DataAuthenticity::new(
        decoder.hash()?,
        decoder.hash()?,
        decoder.bool()?,
    ))
}

fn encode_experiment(
    encoder: &mut Encoder,
    value: &ExperimentArtifact,
) -> Result<(), VvCodecError> {
    encode_header(encoder, value.header())?;
    encode_artifact_id(encoder, value.dataset_id())?;
    encode_experiment_origin(encoder, value.origin())?;
    encoder.count(value.qois().len())?;
    for qoi in value.qois() {
        encode_qoi_id(encoder, qoi)?;
    }
    encode_observation_manifest(encoder, value.manifest())?;
    encoder.count(value.instruments().len())?;
    for instrument in value.instruments() {
        encode_instrument_calibration(encoder, instrument)?;
    }
    encode_clock_synchronization(encoder, value.clocks())?;
    encode_repeatability(encoder, value.repeatability())?;
    encode_data_authenticity(encoder, value.authenticity())
}

fn decode_experiment(decoder: &mut Decoder<'_>) -> Result<ExperimentArtifact, VvCodecError> {
    let header = decode_header(decoder)?;
    let dataset_id = decode_artifact_id(decoder)?;
    let origin = decode_experiment_origin(decoder)?;

    let count = decoder.count()?;
    let mut qois = bounded_vec(decoder, count, "experiment QoIs")?;
    for _ in 0..count {
        let offset = decoder.position();
        let qoi = decode_qoi_id(decoder)?;
        ensure_strictly_increasing(qois.last(), &qoi, offset, "experiment QoI")?;
        qois.push(qoi);
    }

    let manifest = decode_observation_manifest(decoder)?;
    let count = decoder.count()?;
    let mut instruments = bounded_vec(decoder, count, "instrument calibrations")?;
    for _ in 0..count {
        let offset = decoder.position();
        let instrument = decode_instrument_calibration(decoder)?;
        ensure_strictly_increasing(
            instruments.last().map(InstrumentCalibration::instrument_id),
            instrument.instrument_id(),
            offset,
            "instrument calibration id",
        )?;
        instruments.push(instrument);
    }
    let clocks = decode_clock_synchronization(decoder)?;
    let repeatability = decode_repeatability(decoder)?;
    let authenticity = decode_data_authenticity(decoder)?;
    decode_model(
        decoder,
        "experiment artifact",
        ExperimentArtifact::try_new(
            header,
            dataset_id,
            origin,
            qois,
            manifest,
            instruments,
            clocks,
            repeatability,
            authenticity,
        ),
    )
}

fn encode_blind_release(
    encoder: &mut Encoder,
    value: &BlindReleaseReceipt,
) -> Result<(), VvCodecError> {
    encode_artifact_ref(encoder, value.split())?;
    encoder.hash(value.blind_commitment())?;
    encoder.hash(value.authority_receipt_hash())
}

fn decode_blind_release(decoder: &mut Decoder<'_>) -> Result<BlindReleaseReceipt, VvCodecError> {
    let split = decode_artifact_ref(decoder)?;
    let blind_commitment = decoder.hash()?;
    let authority_receipt_hash = decoder.hash()?;
    decode_model(
        decoder,
        "blind-release receipt",
        BlindReleaseReceipt::new(split, blind_commitment, authority_receipt_hash),
    )
}

fn encode_partition(encoder: &mut Encoder, value: &EvidencePartition) -> Result<(), VvCodecError> {
    match value {
        EvidencePartition::Validation => encoder.u8(0),
        EvidencePartition::BlindHoldout { release } => {
            encoder.u8(1)?;
            encode_blind_release(encoder, release)
        }
    }
}

fn decode_partition(decoder: &mut Decoder<'_>) -> Result<EvidencePartition, VvCodecError> {
    let offset = decoder.position();
    match decoder.u8()? {
        0 => Ok(EvidencePartition::Validation),
        1 => Ok(EvidencePartition::BlindHoldout {
            release: decode_blind_release(decoder)?,
        }),
        tag => Err(VvCodecError::at(
            offset,
            format!("unknown evidence-partition tag {tag}"),
        )),
    }
}

fn encode_observation_selection(
    encoder: &mut Encoder,
    value: &ObservationSelection,
) -> Result<(), VvCodecError> {
    encode_artifact_ref(encoder, value.split())?;
    encoder.count(value.ids().len())?;
    for id in value.ids() {
        encode_observation_id(encoder, id)?;
    }
    encode_partition(encoder, value.partition())
}

fn decode_observation_selection(
    decoder: &mut Decoder<'_>,
) -> Result<ObservationSelection, VvCodecError> {
    let split = decode_artifact_ref(decoder)?;
    let count = decoder.count()?;
    let mut ids = bounded_vec(decoder, count, "observation selection")?;
    for _ in 0..count {
        let offset = decoder.position();
        let id = decode_observation_id(decoder)?;
        ensure_strictly_increasing(ids.last(), &id, offset, "observation-selection id")?;
        ids.push(id);
    }
    let partition = decode_partition(decoder)?;
    decode_model(
        decoder,
        "observation selection",
        ObservationSelection::from_canonical(split, ids, partition),
    )
}

// Schema v3 experiment rows bind scientific/metrology meaning as well as the
// complete dataset/locator/extraction source authority. This codec is
// intentionally distinct from the hash-only blind-partition codec below.
fn encode_observation_source_ref(
    encoder: &mut Encoder,
    source: &ObservationSourceRef,
) -> Result<(), VvCodecError> {
    encoder.hash(source.dataset_source_bytes_hash())?;
    encoder.string(source.locator_domain())?;
    encoder.u32(source.locator_contract_version())?;
    encoder.hash(source.locator_hash())?;
    encoder.hash(source.extraction_receipt_hash())
}

fn decode_observation_source_ref(
    decoder: &mut Decoder<'_>,
) -> Result<ObservationSourceRef, VvCodecError> {
    let offset = decoder.position();
    let dataset_source_bytes_hash = decoder.hash()?;
    let locator_domain = decoder.string()?;
    let locator_contract_version = decoder.u32()?;
    let locator_hash = decoder.hash()?;
    let extraction_receipt_hash = decoder.hash()?;
    ObservationSourceRef::try_new(
        dataset_source_bytes_hash,
        locator_domain,
        locator_contract_version,
        locator_hash,
        extraction_receipt_hash,
    )
    .map_err(|error| model_error(offset, "observation source reference", error))
}

fn encode_observation_manifest(
    encoder: &mut Encoder,
    manifest: &ObservationManifest,
) -> Result<(), VvCodecError> {
    encoder.count(manifest.rows().len())?;
    for (id, row) in manifest.rows() {
        encode_observation_id(encoder, id)?;
        encode_observation_source_ref(encoder, row.source_ref())?;
        encode_qoi_id(encoder, row.qoi())?;
        encode_artifact_id(encoder, row.instrument())?;
        encode_artifact_id(encoder, row.acquisition_channel())?;
        encode_artifact_id(encoder, row.clock())?;
    }
    Ok(())
}

fn decode_observation_manifest(
    decoder: &mut Decoder<'_>,
) -> Result<ObservationManifest, VvCodecError> {
    let manifest_offset = decoder.position();
    let count = decoder.count()?;
    let mut rows: Vec<(ObservationId, ObservationManifestRow)> =
        bounded_vec(decoder, count, "experiment manifest rows")?;
    for _ in 0..count {
        let row_offset = decoder.position();
        let id = decode_observation_id(decoder)?;
        ensure_strictly_increasing(
            rows.last().map(|(last, _)| last),
            &id,
            row_offset,
            "experiment manifest row",
        )?;
        let source = decode_observation_source_ref(decoder)?;
        let qoi = decode_qoi_id(decoder)?;
        let instrument = decode_artifact_id(decoder)?;
        let acquisition_channel = decode_artifact_id(decoder)?;
        let clock = decode_artifact_id(decoder)?;
        let row = decode_model(
            decoder,
            "experiment manifest row",
            ObservationManifestRow::try_new(source, qoi, instrument, acquisition_channel, clock),
        )?;
        rows.push((id, row));
    }
    ObservationManifest::try_new(rows)
        .map_err(|error| model_error(manifest_offset, "experiment manifest", error))
}

// bead xl3yi: blind partitions deliberately retain the length-framed
// hash-only `(id, source)` wire form in strict canonical id order. Their v2
// locator commitment is checked against the richer v3 experiment manifest at
// case admission.
fn encode_hash_manifest_rows(
    encoder: &mut Encoder,
    rows: &std::collections::BTreeMap<ObservationId, ContentHash>,
) -> Result<(), VvCodecError> {
    encoder.count(rows.len())?;
    for (id, source) in rows {
        encode_observation_id(encoder, id)?;
        encoder.hash(*source)?;
    }
    Ok(())
}

fn decode_hash_manifest_rows(
    decoder: &mut Decoder<'_>,
    context: &str,
) -> Result<Vec<(ObservationId, ContentHash)>, VvCodecError> {
    let count = decoder.count()?;
    let mut rows: Vec<(ObservationId, ContentHash)> = bounded_vec(decoder, count, context)?;
    for _ in 0..count {
        let offset = decoder.position();
        let id = decode_observation_id(decoder)?;
        ensure_strictly_increasing(rows.last().map(|(last, _)| last), &id, offset, context)?;
        let source = decoder.hash()?;
        rows.push((id, source));
    }
    Ok(rows)
}

fn encode_observation_set(
    encoder: &mut Encoder,
    values: &std::collections::BTreeSet<ObservationId>,
) -> Result<(), VvCodecError> {
    encoder.count(values.len())?;
    for value in values {
        encode_observation_id(encoder, value)?;
    }
    Ok(())
}

fn decode_observation_set(
    decoder: &mut Decoder<'_>,
    context: &str,
) -> Result<Vec<ObservationId>, VvCodecError> {
    let count = decoder.count()?;
    let mut values = bounded_vec(decoder, count, context)?;
    for _ in 0..count {
        let offset = decoder.position();
        let value = decode_observation_id(decoder)?;
        ensure_strictly_increasing(values.last(), &value, offset, context)?;
        values.push(value);
    }
    Ok(values)
}

fn encode_calibration_split(
    encoder: &mut Encoder,
    value: &CalibrationSplit,
) -> Result<(), VvCodecError> {
    encode_header(encoder, value.header())?;
    encode_artifact_ref(encoder, value.experiment())?;
    encoder.hash(value.preregistration_hash())?;
    encode_observation_set(encoder, value.calibration_ids())?;
    encode_observation_set(encoder, value.validation_ids())?;
    encode_hash_manifest_rows(encoder, value.blind_sources())?;
    encoder.hash(value.blind_commitment())
}

fn decode_calibration_split(decoder: &mut Decoder<'_>) -> Result<CalibrationSplit, VvCodecError> {
    let header = decode_header(decoder)?;
    let experiment = decode_artifact_ref(decoder)?;
    let preregistration_hash = decoder.hash()?;
    let calibration = decode_observation_set(decoder, "calibration partition")?;
    let validation = decode_observation_set(decoder, "validation partition")?;
    let blind_holdout = decode_hash_manifest_rows(decoder, "blind-holdout row")?;
    let encoded_commitment = decoder.hash()?;
    let split = decode_model(
        decoder,
        "calibration split",
        CalibrationSplit::try_new(
            header,
            experiment,
            preregistration_hash,
            calibration,
            validation,
            blind_holdout,
        ),
    )?;
    if split.blind_commitment() != encoded_commitment {
        return Err(VvCodecError::at(
            decoder.position(),
            "blind-holdout commitment does not match the canonical partition",
        ));
    }
    Ok(split)
}

fn encode_numerical_uncertainty(
    encoder: &mut Encoder,
    value: &NumericalUncertainty,
) -> Result<(), VvCodecError> {
    encoder.f64(value.half_width())?;
    encoder.hash(value.evidence_hash())
}

fn decode_numerical_uncertainty(
    decoder: &mut Decoder<'_>,
) -> Result<NumericalUncertainty, VvCodecError> {
    let half_width = decoder.f64()?;
    let evidence_hash = decoder.hash()?;
    decode_model(
        decoder,
        "numerical uncertainty",
        NumericalUncertainty::try_new(half_width, evidence_hash),
    )
}

fn encode_solution_verification(
    encoder: &mut Encoder,
    value: &SolutionVerificationReceipt,
) -> Result<(), VvCodecError> {
    encode_header(encoder, value.header())?;
    encode_artifact_id(encoder, value.solve_id())?;
    encode_qoi_id(encoder, value.qoi())?;
    encode_unit_id(encoder, value.unit())?;
    encode_numerical_uncertainty(encoder, value.mesh())?;
    encode_numerical_uncertainty(encoder, value.time())?;
    encode_numerical_uncertainty(encoder, value.nonlinear())?;
    encode_numerical_uncertainty(encoder, value.iterative())?;
    encoder.f64(value.combined_half_width())
}

fn decode_solution_verification(
    decoder: &mut Decoder<'_>,
) -> Result<SolutionVerificationReceipt, VvCodecError> {
    let header = decode_header(decoder)?;
    let solve_id = decode_artifact_id(decoder)?;
    let qoi = decode_qoi_id(decoder)?;
    let unit = decode_unit_id(decoder)?;
    let mesh = decode_numerical_uncertainty(decoder)?;
    let time = decode_numerical_uncertainty(decoder)?;
    let nonlinear = decode_numerical_uncertainty(decoder)?;
    let iterative = decode_numerical_uncertainty(decoder)?;
    let encoded_combined = decoder.f64()?;
    let receipt = decode_model(
        decoder,
        "solution-verification receipt",
        SolutionVerificationReceipt::try_new(
            header, solve_id, qoi, unit, mesh, time, nonlinear, iterative,
        ),
    )?;
    if receipt.combined_half_width().to_bits() != encoded_combined.to_bits() {
        return Err(VvCodecError::at(
            decoder.position(),
            "combined numerical uncertainty is not the canonical derived value",
        ));
    }
    Ok(receipt)
}

fn encode_evidence_target(
    encoder: &mut Encoder,
    value: &EvidenceTarget,
) -> Result<(), VvCodecError> {
    match value {
        EvidenceTarget::VvArtifact(reference) => {
            encoder.u8(0)?;
            encode_artifact_ref(encoder, reference)
        }
        EvidenceTarget::External { family, id, hash } => {
            encoder.u8(1)?;
            encode_artifact_id(encoder, family)?;
            encode_artifact_id(encoder, id)?;
            encoder.hash(*hash)
        }
    }
}

fn decode_evidence_target(decoder: &mut Decoder<'_>) -> Result<EvidenceTarget, VvCodecError> {
    let offset = decoder.position();
    match decoder.u8()? {
        0 => Ok(EvidenceTarget::VvArtifact(decode_artifact_ref(decoder)?)),
        1 => Ok(EvidenceTarget::External {
            family: decode_artifact_id(decoder)?,
            id: decode_artifact_id(decoder)?,
            hash: decoder.hash()?,
        }),
        tag => Err(VvCodecError::at(
            offset,
            format!("unknown evidence-target tag {tag}"),
        )),
    }
}

fn encode_dependency_role(
    encoder: &mut Encoder,
    value: DependencyRole,
) -> Result<(), VvCodecError> {
    encoder.u8(match value {
        DependencyRole::CodeVerification => 0,
        DependencyRole::SolutionVerification => 1,
        DependencyRole::PhysicalValidation => 2,
        DependencyRole::ModelDiscrepancy => 3,
        DependencyRole::ParameterData => 4,
        DependencyRole::PosteriorPredictive => 5,
        DependencyRole::ProcessConformance => 6,
    })
}

fn decode_dependency_role(decoder: &mut Decoder<'_>) -> Result<DependencyRole, VvCodecError> {
    let offset = decoder.position();
    match decoder.u8()? {
        0 => Ok(DependencyRole::CodeVerification),
        1 => Ok(DependencyRole::SolutionVerification),
        2 => Ok(DependencyRole::PhysicalValidation),
        3 => Ok(DependencyRole::ModelDiscrepancy),
        4 => Ok(DependencyRole::ParameterData),
        5 => Ok(DependencyRole::PosteriorPredictive),
        6 => Ok(DependencyRole::ProcessConformance),
        tag => Err(VvCodecError::at(
            offset,
            format!("unknown dependency-role tag {tag}"),
        )),
    }
}

fn encode_evidence_dependency(
    encoder: &mut Encoder,
    value: &EvidenceDependency,
) -> Result<(), VvCodecError> {
    encode_qoi_id(encoder, value.qoi())?;
    encode_dependency_role(encoder, value.role())?;
    encode_evidence_target(encoder, value.target())?;
    encoder.bool(value.observations().is_some())?;
    if let Some(observations) = value.observations() {
        encode_observation_selection(encoder, observations)?;
    }
    Ok(())
}

fn decode_evidence_dependency(
    decoder: &mut Decoder<'_>,
) -> Result<EvidenceDependency, VvCodecError> {
    let qoi = decode_qoi_id(decoder)?;
    let role = decode_dependency_role(decoder)?;
    let target = decode_evidence_target(decoder)?;
    let observations = if decoder.bool()? {
        Some(decode_observation_selection(decoder)?)
    } else {
        None
    };
    match (role, target, observations) {
        (
            DependencyRole::PhysicalValidation,
            EvidenceTarget::VvArtifact(experiment),
            Some(observations),
        ) => Ok(EvidenceDependency::physical_validation(
            qoi,
            experiment,
            observations,
        )),
        (role, target, None) => Ok(EvidenceDependency::new(qoi, role, target)),
        _ => Err(VvCodecError::at(
            decoder.position(),
            "observation selection is only canonical on a physical V&V-artifact dependency",
        )),
    }
}

fn encode_prediction_uncertainty_kind(
    encoder: &mut Encoder,
    value: PredictionUncertaintyKind,
) -> Result<(), VvCodecError> {
    encoder.u8(match value {
        PredictionUncertaintyKind::ModelForm => 0,
        PredictionUncertaintyKind::Parameter => 1,
        PredictionUncertaintyKind::Numerical => 2,
        PredictionUncertaintyKind::Data => 3,
        PredictionUncertaintyKind::Aleatory => 4,
        PredictionUncertaintyKind::Epistemic => 5,
    })
}

fn decode_prediction_uncertainty_kind(
    decoder: &mut Decoder<'_>,
) -> Result<PredictionUncertaintyKind, VvCodecError> {
    let offset = decoder.position();
    match decoder.u8()? {
        0 => Ok(PredictionUncertaintyKind::ModelForm),
        1 => Ok(PredictionUncertaintyKind::Parameter),
        2 => Ok(PredictionUncertaintyKind::Numerical),
        3 => Ok(PredictionUncertaintyKind::Data),
        4 => Ok(PredictionUncertaintyKind::Aleatory),
        5 => Ok(PredictionUncertaintyKind::Epistemic),
        tag => Err(VvCodecError::at(
            offset,
            format!("unknown prediction-uncertainty-kind tag {tag}"),
        )),
    }
}

fn encode_uncertainty_term(
    encoder: &mut Encoder,
    value: &UncertaintyTerm,
) -> Result<(), VvCodecError> {
    encode_prediction_uncertainty_kind(encoder, value.kind())?;
    encoder.f64(value.magnitude())?;
    encode_evidence_target(encoder, value.source())
}

fn decode_uncertainty_term(decoder: &mut Decoder<'_>) -> Result<UncertaintyTerm, VvCodecError> {
    let kind = decode_prediction_uncertainty_kind(decoder)?;
    let magnitude = decoder.f64()?;
    let source = decode_evidence_target(decoder)?;
    decode_model(
        decoder,
        "uncertainty term",
        UncertaintyTerm::try_new(kind, magnitude, source),
    )
}

fn encode_correlation_matrix(
    encoder: &mut Encoder,
    value: &CorrelationMatrix,
) -> Result<(), VvCodecError> {
    encoder.usize(value.dimension())?;
    encoder.count(value.values().len())?;
    for entry in value.values() {
        encoder.f64(*entry)?;
    }
    Ok(())
}

fn decode_correlation_matrix(decoder: &mut Decoder<'_>) -> Result<CorrelationMatrix, VvCodecError> {
    let dimension = decoder.usize()?;
    let count = decoder.count()?;
    let mut values = bounded_vec(decoder, count, "correlation entries")?;
    for _ in 0..count {
        values.push(decoder.f64()?);
    }
    decode_model(
        decoder,
        "correlation matrix",
        CorrelationMatrix::try_new(dimension, values),
    )
}

fn encode_waterfall_mode(encoder: &mut Encoder, value: &WaterfallMode) -> Result<(), VvCodecError> {
    match value {
        WaterfallMode::GuaranteedBound => encoder.u8(0),
        WaterfallMode::Probabilistic {
            confidence,
            dependence,
        } => {
            encoder.u8(1)?;
            encoder.f64(*confidence)?;
            encode_correlation_matrix(encoder, dependence)
        }
    }
}

fn decode_waterfall_mode(decoder: &mut Decoder<'_>) -> Result<WaterfallMode, VvCodecError> {
    let offset = decoder.position();
    match decoder.u8()? {
        0 => Ok(WaterfallMode::GuaranteedBound),
        1 => Ok(WaterfallMode::Probabilistic {
            confidence: decoder.f64()?,
            dependence: decode_correlation_matrix(decoder)?,
        }),
        tag => Err(VvCodecError::at(
            offset,
            format!("unknown waterfall-mode tag {tag}"),
        )),
    }
}

fn encode_waterfall(
    encoder: &mut Encoder,
    value: &UncertaintyWaterfall,
) -> Result<(), VvCodecError> {
    encode_qoi_id(encoder, value.qoi())?;
    encode_unit_id(encoder, value.unit())?;
    encode_waterfall_mode(encoder, value.mode())?;
    encoder.count(value.terms().len())?;
    for term in value.terms() {
        encode_uncertainty_term(encoder, term)?;
    }
    encoder.f64(value.total())
}

fn decode_waterfall(decoder: &mut Decoder<'_>) -> Result<UncertaintyWaterfall, VvCodecError> {
    let qoi = decode_qoi_id(decoder)?;
    let unit = decode_unit_id(decoder)?;
    let mode = decode_waterfall_mode(decoder)?;
    let count = decoder.count()?;
    let mut terms: Vec<UncertaintyTerm> = bounded_vec(decoder, count, "waterfall terms")?;
    for _ in 0..count {
        let offset = decoder.position();
        let term = decode_uncertainty_term(decoder)?;
        if terms
            .last()
            .is_some_and(|previous| previous.kind() >= term.kind())
        {
            return Err(VvCodecError::at(
                offset,
                "waterfall uncertainty kinds are duplicated or out of canonical order",
            ));
        }
        terms.push(term);
    }
    let encoded_total = decoder.f64()?;
    let waterfall = decode_model(
        decoder,
        "uncertainty waterfall",
        UncertaintyWaterfall::try_new(qoi, unit, mode, terms),
    )?;
    if waterfall.total().to_bits() != encoded_total.to_bits() {
        return Err(VvCodecError::at(
            decoder.position(),
            "waterfall total is not the canonical derived value",
        ));
    }
    Ok(waterfall)
}

fn encode_validation_metric(
    encoder: &mut Encoder,
    value: &ValidationMetric,
) -> Result<(), VvCodecError> {
    encode_artifact_id(encoder, value.name())?;
    encode_qoi_id(encoder, value.qoi())?;
    encode_observation_selection(encoder, value.observations())?;
    encoder.f64(value.observed())?;
    encoder.f64(value.predicted())?;
    encoder.f64(value.experimental_uncertainty())?;
    encoder.f64(value.numerical_uncertainty())?;
    encoder.f64(value.combined_uncertainty())
}

fn decode_validation_metric(decoder: &mut Decoder<'_>) -> Result<ValidationMetric, VvCodecError> {
    let name = decode_artifact_id(decoder)?;
    let qoi = decode_qoi_id(decoder)?;
    let observations = decode_observation_selection(decoder)?;
    let observed = decoder.f64()?;
    let predicted = decoder.f64()?;
    let experimental_uncertainty = decoder.f64()?;
    let numerical_uncertainty = decoder.f64()?;
    let encoded_combined = decoder.f64()?;
    let metric = decode_model(
        decoder,
        "validation metric",
        ValidationMetric::try_new(
            name,
            qoi,
            observations,
            observed,
            predicted,
            experimental_uncertainty,
            numerical_uncertainty,
        ),
    )?;
    if metric.combined_uncertainty().to_bits() != encoded_combined.to_bits() {
        return Err(VvCodecError::at(
            decoder.position(),
            "combined validation uncertainty is not the canonical derived value",
        ));
    }
    Ok(metric)
}

fn encode_posterior_check(
    encoder: &mut Encoder,
    value: &PosteriorPredictiveCheck,
) -> Result<(), VvCodecError> {
    encode_artifact_id(encoder, value.name())?;
    encode_qoi_id(encoder, value.qoi())?;
    encode_observation_selection(encoder, value.observations())?;
    encoder.f64(value.tail_probability())?;
    encoder.f64(value.minimum_tail_probability())?;
    encoder.hash(value.artifact_hash())
}

fn decode_posterior_check(
    decoder: &mut Decoder<'_>,
) -> Result<PosteriorPredictiveCheck, VvCodecError> {
    let name = decode_artifact_id(decoder)?;
    let qoi = decode_qoi_id(decoder)?;
    let observations = decode_observation_selection(decoder)?;
    let tail_probability = decoder.f64()?;
    let minimum_tail_probability = decoder.f64()?;
    let artifact_hash = decoder.hash()?;
    decode_model(
        decoder,
        "posterior-predictive check",
        PosteriorPredictiveCheck::try_new(
            name,
            qoi,
            observations,
            tail_probability,
            minimum_tail_probability,
            artifact_hash,
        ),
    )
}

fn encode_evidence_axis(encoder: &mut Encoder, value: EvidenceAxis) -> Result<(), VvCodecError> {
    encoder.u8(match value {
        EvidenceAxis::CodeVerification => 0,
        EvidenceAxis::SolutionVerification => 1,
        EvidenceAxis::NumericalUncertainty => 2,
        EvidenceAxis::ParameterDataUncertainty => 3,
        EvidenceAxis::ModelFormValidation => 4,
        EvidenceAxis::PredictionDomainRelevance => 5,
        EvidenceAxis::ComparisonToExperiment => 6,
    })
}

fn decode_evidence_axis(decoder: &mut Decoder<'_>) -> Result<EvidenceAxis, VvCodecError> {
    let offset = decoder.position();
    match decoder.u8()? {
        0 => Ok(EvidenceAxis::CodeVerification),
        1 => Ok(EvidenceAxis::SolutionVerification),
        2 => Ok(EvidenceAxis::NumericalUncertainty),
        3 => Ok(EvidenceAxis::ParameterDataUncertainty),
        4 => Ok(EvidenceAxis::ModelFormValidation),
        5 => Ok(EvidenceAxis::PredictionDomainRelevance),
        6 => Ok(EvidenceAxis::ComparisonToExperiment),
        tag => Err(VvCodecError::at(
            offset,
            format!("unknown evidence-axis tag {tag}"),
        )),
    }
}

fn encode_evidence_axis_status(
    encoder: &mut Encoder,
    value: &EvidenceAxisStatus,
) -> Result<(), VvCodecError> {
    match value {
        EvidenceAxisStatus::Present { artifacts } => {
            encoder.u8(0)?;
            encoder.count(artifacts.len())?;
            for artifact in artifacts {
                encoder.hash(*artifact)?;
            }
            Ok(())
        }
        EvidenceAxisStatus::Missing { reason } => {
            encoder.u8(1)?;
            encoder.string(reason)
        }
        EvidenceAxisStatus::Refused { rule, reason } => {
            encoder.u8(2)?;
            encode_rule(encoder, *rule)?;
            encoder.string(reason)
        }
    }
}

fn decode_evidence_axis_status(
    decoder: &mut Decoder<'_>,
) -> Result<EvidenceAxisStatus, VvCodecError> {
    let offset = decoder.position();
    match decoder.u8()? {
        0 => {
            let count = decoder.count()?;
            let mut artifacts = bounded_vec(decoder, count, "evidence-axis artifacts")?;
            for _ in 0..count {
                let offset = decoder.position();
                let artifact = decoder.hash()?;
                ensure_strictly_increasing(
                    artifacts.last(),
                    &artifact,
                    offset,
                    "evidence-axis artifact",
                )?;
                artifacts.push(artifact);
            }
            Ok(EvidenceAxisStatus::Present { artifacts })
        }
        1 => Ok(EvidenceAxisStatus::Missing {
            reason: decoder.string()?,
        }),
        2 => Ok(EvidenceAxisStatus::Refused {
            rule: decode_rule(decoder)?,
            reason: decoder.string()?,
        }),
        tag => Err(VvCodecError::at(
            offset,
            format!("unknown evidence-axis-status tag {tag}"),
        )),
    }
}

fn encode_evidence_axes(encoder: &mut Encoder, value: &EvidenceAxes) -> Result<(), VvCodecError> {
    encoder.count(value.axes().len())?;
    for (axis, status) in value.axes() {
        encode_evidence_axis(encoder, *axis)?;
        encode_evidence_axis_status(encoder, status)?;
    }
    Ok(())
}

fn decode_evidence_axes(decoder: &mut Decoder<'_>) -> Result<EvidenceAxes, VvCodecError> {
    let count = decoder.count()?;
    let mut rows = bounded_vec(decoder, count, "evidence axes")?;
    for _ in 0..count {
        let offset = decoder.position();
        let axis = decode_evidence_axis(decoder)?;
        if rows.last().is_some_and(|(previous, _)| *previous >= axis) {
            return Err(VvCodecError::at(
                offset,
                "evidence axes are duplicated or out of canonical order",
            ));
        }
        rows.push((axis, decode_evidence_axis_status(decoder)?));
    }
    decode_model(decoder, "evidence axes", EvidenceAxes::try_new(rows))
}

fn encode_prediction_assessment(
    encoder: &mut Encoder,
    value: &PredictionAssessment,
) -> Result<(), VvCodecError> {
    encode_header(encoder, value.header())?;
    encode_artifact_ref(encoder, value.context())?;
    encode_artifact_ref(encoder, value.validation_plan())?;
    encode_qoi_id(encoder, value.qoi())?;
    encoder.count(value.dependencies().len())?;
    for dependency in value.dependencies() {
        encode_evidence_dependency(encoder, dependency)?;
    }
    encode_waterfall(encoder, value.waterfall())?;
    encoder.count(value.validation_metrics().len())?;
    for metric in value.validation_metrics() {
        encode_validation_metric(encoder, metric)?;
    }
    encoder.count(value.posterior_checks().len())?;
    for check in value.posterior_checks() {
        encode_posterior_check(encoder, check)?;
    }
    encode_applicability_point(encoder, value.applicability_point())?;
    encode_applicability_decision(encoder, value.applicability())?;
    encode_evidence_axes(encoder, value.evidence_axes())?;
    encoder.count(value.assumption_checks().len())?;
    for (id, passed) in value.assumption_checks() {
        encode_assumption_id(encoder, id)?;
        encoder.bool(*passed)?;
    }
    Ok(())
}

fn decode_prediction_assessment(
    decoder: &mut Decoder<'_>,
) -> Result<PredictionAssessment, VvCodecError> {
    let header = decode_header(decoder)?;
    let context = decode_artifact_ref(decoder)?;
    let validation_plan = decode_artifact_ref(decoder)?;
    let qoi = decode_qoi_id(decoder)?;

    let count = decoder.count()?;
    let mut dependencies = bounded_vec(decoder, count, "prediction dependencies")?;
    for _ in 0..count {
        let offset = decoder.position();
        let dependency = decode_evidence_dependency(decoder)?;
        ensure_strictly_increasing(
            dependencies.last(),
            &dependency,
            offset,
            "prediction dependency",
        )?;
        dependencies.push(dependency);
    }
    let waterfall = decode_waterfall(decoder)?;

    let count = decoder.count()?;
    let mut validation_metrics = bounded_vec(decoder, count, "validation metrics")?;
    for _ in 0..count {
        let offset = decoder.position();
        let metric = decode_validation_metric(decoder)?;
        ensure_strictly_increasing(
            validation_metrics.last().map(ValidationMetric::name),
            metric.name(),
            offset,
            "validation metric id",
        )?;
        validation_metrics.push(metric);
    }

    let count = decoder.count()?;
    let mut posterior_checks = bounded_vec(decoder, count, "posterior-predictive checks")?;
    for _ in 0..count {
        let offset = decoder.position();
        let check = decode_posterior_check(decoder)?;
        ensure_strictly_increasing(
            posterior_checks.last().map(PosteriorPredictiveCheck::name),
            check.name(),
            offset,
            "posterior-predictive check id",
        )?;
        posterior_checks.push(check);
    }

    let applicability_point = decode_applicability_point(decoder)?;
    let applicability = decode_applicability_decision(decoder)?;
    let evidence_axes = decode_evidence_axes(decoder)?;

    let count = decoder.count()?;
    let mut assumption_checks = bounded_vec(decoder, count, "prediction assumption checks")?;
    for _ in 0..count {
        let offset = decoder.position();
        let id = decode_assumption_id(decoder)?;
        ensure_strictly_increasing(
            assumption_checks.last().map(|(id, _)| id),
            &id,
            offset,
            "prediction assumption-check id",
        )?;
        assumption_checks.push((id, decoder.bool()?));
    }

    decode_model(
        decoder,
        "prediction assessment",
        PredictionAssessment::try_new(
            header,
            context,
            validation_plan,
            qoi,
            dependencies,
            waterfall,
            validation_metrics,
            posterior_checks,
            applicability_point,
            applicability,
            evidence_axes,
            assumption_checks,
        ),
    )
}

fn encode_assumption_evidence(
    encoder: &mut Encoder,
    value: &AssumptionEvidence,
) -> Result<(), VvCodecError> {
    encoder.string(value.requirement())?;
    encoder.bool(value.artifact().is_some())?;
    if let Some(artifact) = value.artifact() {
        encode_evidence_target(encoder, artifact)?;
    }
    Ok(())
}

fn decode_assumption_evidence(
    decoder: &mut Decoder<'_>,
) -> Result<AssumptionEvidence, VvCodecError> {
    let requirement = decoder.string()?;
    let artifact = if decoder.bool()? {
        Some(decode_evidence_target(decoder)?)
    } else {
        None
    };
    decode_model(
        decoder,
        "assumption evidence",
        AssumptionEvidence::try_new(requirement, artifact),
    )
}

fn encode_runtime_monitor(
    encoder: &mut Encoder,
    value: &RuntimeMonitorSpec,
) -> Result<(), VvCodecError> {
    encoder.string(value.signal())?;
    encoder.bool(value.evidence_hash().is_some())?;
    if let Some(hash) = value.evidence_hash() {
        encoder.hash(hash)?;
    }
    Ok(())
}

fn decode_runtime_monitor(decoder: &mut Decoder<'_>) -> Result<RuntimeMonitorSpec, VvCodecError> {
    let signal = decoder.string()?;
    let evidence_hash = if decoder.bool()? {
        Some(decoder.hash()?)
    } else {
        None
    };
    decode_model(
        decoder,
        "runtime monitor",
        RuntimeMonitorSpec::try_new(signal, evidence_hash),
    )
}

fn encode_violation_effect(
    encoder: &mut Encoder,
    value: &ViolationEffect,
) -> Result<(), VvCodecError> {
    match value {
        ViolationEffect::Demote { reason } => {
            encoder.u8(0)?;
            encoder.string(reason)
        }
        ViolationEffect::EscalateOrRefuse { target_lane } => {
            encoder.u8(1)?;
            encode_artifact_id(encoder, target_lane)
        }
        ViolationEffect::Refuse { reason } => {
            encoder.u8(2)?;
            encoder.string(reason)
        }
    }
}

fn decode_violation_effect(decoder: &mut Decoder<'_>) -> Result<ViolationEffect, VvCodecError> {
    let offset = decoder.position();
    match decoder.u8()? {
        0 => Ok(ViolationEffect::Demote {
            reason: decoder.string()?,
        }),
        1 => Ok(ViolationEffect::EscalateOrRefuse {
            target_lane: decode_artifact_id(decoder)?,
        }),
        2 => Ok(ViolationEffect::Refuse {
            reason: decoder.string()?,
        }),
        tag => Err(VvCodecError::at(
            offset,
            format!("unknown assumption-violation-effect tag {tag}"),
        )),
    }
}

fn encode_review_gate(encoder: &mut Encoder, value: &ReviewGate) -> Result<(), VvCodecError> {
    match value {
        ReviewGate::Phase { gate } => {
            encoder.u8(0)?;
            encode_artifact_id(encoder, gate)
        }
        ReviewGate::EverySolve => encoder.u8(1),
        ReviewGate::EveryQuery => encoder.u8(2),
        ReviewGate::EveryUpdate => encoder.u8(3),
    }
}

fn decode_review_gate(decoder: &mut Decoder<'_>) -> Result<ReviewGate, VvCodecError> {
    let offset = decoder.position();
    match decoder.u8()? {
        0 => Ok(ReviewGate::Phase {
            gate: decode_artifact_id(decoder)?,
        }),
        1 => Ok(ReviewGate::EverySolve),
        2 => Ok(ReviewGate::EveryQuery),
        3 => Ok(ReviewGate::EveryUpdate),
        tag => Err(VvCodecError::at(
            offset,
            format!("unknown review-gate tag {tag}"),
        )),
    }
}

fn encode_assumption_row(encoder: &mut Encoder, value: &AssumptionRow) -> Result<(), VvCodecError> {
    encode_assumption_id(encoder, value.id())?;
    encoder.string(value.predicate())?;
    encoder.string(value.scope())?;
    encode_assumption_evidence(encoder, value.evidence())?;
    encode_runtime_monitor(encoder, value.monitor())?;
    encode_violation_effect(encoder, value.violation_effect())?;
    encode_artifact_id(encoder, value.owner())?;
    encode_review_gate(encoder, value.review_gate())
}

fn decode_assumption_row(decoder: &mut Decoder<'_>) -> Result<AssumptionRow, VvCodecError> {
    let id = decode_assumption_id(decoder)?;
    let predicate = decoder.string()?;
    let scope = decoder.string()?;
    let evidence = decode_assumption_evidence(decoder)?;
    let monitor = decode_runtime_monitor(decoder)?;
    let violation_effect = decode_violation_effect(decoder)?;
    let owner = decode_artifact_id(decoder)?;
    let review_gate = decode_review_gate(decoder)?;
    decode_model(
        decoder,
        "assumption row",
        AssumptionRow::try_new(
            id,
            predicate,
            scope,
            evidence,
            monitor,
            violation_effect,
            owner,
            review_gate,
        ),
    )
}

fn encode_assumptions_ledger(
    encoder: &mut Encoder,
    value: &AssumptionsLedger,
) -> Result<(), VvCodecError> {
    encode_header(encoder, value.header())?;
    encoder.count(value.rows().len())?;
    for row in value.rows().values() {
        encode_assumption_row(encoder, row)?;
    }
    Ok(())
}

fn decode_assumptions_ledger(decoder: &mut Decoder<'_>) -> Result<AssumptionsLedger, VvCodecError> {
    let header = decode_header(decoder)?;
    let count = decoder.count()?;
    let mut rows = bounded_vec(decoder, count, "assumption rows")?;
    for _ in 0..count {
        let offset = decoder.position();
        let row = decode_assumption_row(decoder)?;
        ensure_strictly_increasing(
            rows.last().map(AssumptionRow::id),
            row.id(),
            offset,
            "assumption-row id",
        )?;
        rows.push(row);
    }
    decode_model(
        decoder,
        "assumptions ledger",
        AssumptionsLedger::try_new(header, rows),
    )
}

fn expect_root(decoder: &mut Decoder<'_>, expected: u8, context: &str) -> Result<(), VvCodecError> {
    let offset = decoder.position();
    let actual = decoder.u8()?;
    if actual == expected {
        Ok(())
    } else {
        Err(VvCodecError::at(
            offset,
            format!("expected {context} root tag {expected}, found {actual}"),
        ))
    }
}

fn encode_artifact_payload(encoder: &mut Encoder, value: &VvArtifact) -> Result<(), VvCodecError> {
    encode_artifact_kind(encoder, value.kind())?;
    match value {
        VvArtifact::ContextOfUse(value) => encode_context(encoder, value),
        VvArtifact::ValidationPlan(value) => encode_validation_plan(encoder, value),
        VvArtifact::ExperimentArtifact(value) => encode_experiment(encoder, value),
        VvArtifact::CalibrationSplit(value) => encode_calibration_split(encoder, value),
        VvArtifact::SolutionVerificationReceipt(value) => {
            encode_solution_verification(encoder, value)
        }
        VvArtifact::PredictionAssessment(value) => encode_prediction_assessment(encoder, value),
        VvArtifact::AssumptionsLedger(value) => encode_assumptions_ledger(encoder, value),
    }
}

fn decode_artifact_payload(decoder: &mut Decoder<'_>) -> Result<VvArtifact, VvCodecError> {
    match decode_artifact_kind(decoder)? {
        ArtifactKind::ContextOfUse => Ok(VvArtifact::ContextOfUse(decode_context(decoder)?)),
        ArtifactKind::ValidationPlan => {
            Ok(VvArtifact::ValidationPlan(decode_validation_plan(decoder)?))
        }
        ArtifactKind::ExperimentArtifact => {
            Ok(VvArtifact::ExperimentArtifact(decode_experiment(decoder)?))
        }
        ArtifactKind::CalibrationSplit => Ok(VvArtifact::CalibrationSplit(
            decode_calibration_split(decoder)?,
        )),
        ArtifactKind::SolutionVerificationReceipt => Ok(VvArtifact::SolutionVerificationReceipt(
            decode_solution_verification(decoder)?,
        )),
        ArtifactKind::PredictionAssessment => Ok(VvArtifact::PredictionAssessment(
            decode_prediction_assessment(decoder)?,
        )),
        ArtifactKind::AssumptionsLedger => Ok(VvArtifact::AssumptionsLedger(
            decode_assumptions_ledger(decoder)?,
        )),
    }
}

fn canonical_artifact_bytes(
    kind: ArtifactKind,
    encode: impl FnOnce(&mut Encoder) -> Result<(), VvCodecError>,
) -> Result<Vec<u8>, VvCodecError> {
    let mut encoder = Encoder::new()?;
    encoder.u8(ROOT_ARTIFACT)?;
    encode_artifact_kind(&mut encoder, kind)?;
    encode(&mut encoder)?;
    Ok(encoder.finish())
}

fn content_hash_for(bytes: &[u8]) -> ContentHash {
    fs_blake3::hash_domain(VV_ARTIFACT_FAMILY, bytes)
}

fn case_content_hash_for(bytes: &[u8]) -> ContentHash {
    fs_blake3::hash_domain(VV_CASE_FAMILY, bytes)
}

impl VvArtifact {
    /// Encode one top-level V&V artifact into the exact bounded transport.
    pub fn canonical_bytes(&self) -> Result<Vec<u8>, VvCodecError> {
        let mut encoder = Encoder::new()?;
        encoder.u8(ROOT_ARTIFACT)?;
        encode_artifact_payload(&mut encoder, self)?;
        Ok(encoder.finish())
    }

    /// Decode one exact current top-level artifact transport.
    pub fn from_canonical_bytes(bytes: &[u8]) -> Result<Self, VvCodecError> {
        let mut decoder = Decoder::new(bytes)?;
        expect_root(&mut decoder, ROOT_ARTIFACT, "artifact")?;
        let artifact = decode_artifact_payload(&mut decoder)?;
        decoder.finish()?;
        if artifact.canonical_bytes()?.as_slice() != bytes {
            return Err(VvCodecError::at(
                0,
                "transport is structurally valid but not a canonical fixed point",
            ));
        }
        Ok(artifact)
    }

    /// Domain-separated identity of this exact canonical artifact.
    pub fn content_hash(&self) -> Result<ContentHash, VvCodecError> {
        Ok(content_hash_for(&self.canonical_bytes()?))
    }
}

macro_rules! concrete_artifact_codec {
    ($type:ty, $kind:expr, $variant:ident, $encode:ident) => {
        impl $type {
            /// Encode this artifact into the exact current bounded transport.
            pub fn canonical_bytes(&self) -> Result<Vec<u8>, VvCodecError> {
                canonical_artifact_bytes($kind, |encoder| $encode(encoder, self))
            }

            /// Decode this exact concrete artifact kind from canonical transport.
            pub fn from_canonical_bytes(bytes: &[u8]) -> Result<Self, VvCodecError> {
                match VvArtifact::from_canonical_bytes(bytes)? {
                    VvArtifact::$variant(value) => Ok(value),
                    other => Err(VvCodecError::at(
                        13,
                        format!(
                            "expected {} artifact, found {}",
                            $kind.slug(),
                            other.kind().slug()
                        ),
                    )),
                }
            }

            /// Domain-separated identity of this exact canonical artifact.
            pub fn content_hash(&self) -> Result<ContentHash, VvCodecError> {
                Ok(content_hash_for(&self.canonical_bytes()?))
            }
        }
    };
}

concrete_artifact_codec!(
    ContextOfUse,
    ArtifactKind::ContextOfUse,
    ContextOfUse,
    encode_context
);
concrete_artifact_codec!(
    ValidationPlan,
    ArtifactKind::ValidationPlan,
    ValidationPlan,
    encode_validation_plan
);
concrete_artifact_codec!(
    ExperimentArtifact,
    ArtifactKind::ExperimentArtifact,
    ExperimentArtifact,
    encode_experiment
);
concrete_artifact_codec!(
    CalibrationSplit,
    ArtifactKind::CalibrationSplit,
    CalibrationSplit,
    encode_calibration_split
);
concrete_artifact_codec!(
    SolutionVerificationReceipt,
    ArtifactKind::SolutionVerificationReceipt,
    SolutionVerificationReceipt,
    encode_solution_verification
);
concrete_artifact_codec!(
    PredictionAssessment,
    ArtifactKind::PredictionAssessment,
    PredictionAssessment,
    encode_prediction_assessment
);
concrete_artifact_codec!(
    AssumptionsLedger,
    ArtifactKind::AssumptionsLedger,
    AssumptionsLedger,
    encode_assumptions_ledger
);

fn encode_case_body(encoder: &mut Encoder, value: &VvCase) -> Result<(), VvCodecError> {
    encode_context(encoder, value.context())?;
    encode_validation_plan(encoder, value.validation_plan())?;

    encoder.count(value.experiments().len())?;
    for experiment in value.experiments().values() {
        encode_experiment(encoder, experiment)?;
    }

    encoder.count(value.splits().len())?;
    for split in value.splits().values() {
        encode_calibration_split(encoder, split)?;
    }

    encoder.count(value.solution_verification().len())?;
    for receipt in value.solution_verification().values() {
        encode_solution_verification(encoder, receipt)?;
    }

    encoder.count(value.predictions().len())?;
    for prediction in value.predictions().values() {
        encode_prediction_assessment(encoder, prediction)?;
    }

    encode_assumptions_ledger(encoder, value.assumptions())
}

fn decode_case_body(decoder: &mut Decoder<'_>) -> Result<VvCase, VvCodecError> {
    let context = decode_context(decoder)?;
    let validation_plan = decode_validation_plan(decoder)?;

    let count = decoder.count()?;
    let mut experiments = bounded_vec(decoder, count, "case experiments")?;
    for _ in 0..count {
        let offset = decoder.position();
        let artifact = decode_experiment(decoder)?;
        ensure_strictly_increasing(
            experiments.last().map(ExperimentArtifact::id),
            artifact.id(),
            offset,
            "case experiment id",
        )?;
        experiments.push(artifact);
    }

    let count = decoder.count()?;
    let mut splits = bounded_vec(decoder, count, "case splits")?;
    for _ in 0..count {
        let offset = decoder.position();
        let artifact = decode_calibration_split(decoder)?;
        ensure_strictly_increasing(
            splits.last().map(CalibrationSplit::id),
            artifact.id(),
            offset,
            "case split id",
        )?;
        splits.push(artifact);
    }

    let count = decoder.count()?;
    let mut solution_verification =
        bounded_vec(decoder, count, "case solution-verification receipts")?;
    for _ in 0..count {
        let offset = decoder.position();
        let artifact = decode_solution_verification(decoder)?;
        ensure_strictly_increasing(
            solution_verification
                .last()
                .map(SolutionVerificationReceipt::id),
            artifact.id(),
            offset,
            "case solution-verification id",
        )?;
        solution_verification.push(artifact);
    }

    let count = decoder.count()?;
    let mut predictions = bounded_vec(decoder, count, "case predictions")?;
    for _ in 0..count {
        let offset = decoder.position();
        let artifact = decode_prediction_assessment(decoder)?;
        ensure_strictly_increasing(
            predictions.last().map(PredictionAssessment::id),
            artifact.id(),
            offset,
            "case prediction id",
        )?;
        predictions.push(artifact);
    }

    let assumptions = decode_assumptions_ledger(decoder)?;
    decode_model(
        decoder,
        "V&V case",
        VvCase::try_new(
            context,
            validation_plan,
            experiments,
            splits,
            solution_verification,
            predictions,
            assumptions,
        ),
    )
}

impl VvCase {
    /// Encode this complete case into the exact current bounded transport.
    pub fn canonical_bytes(&self) -> Result<Vec<u8>, VvCodecError> {
        let mut encoder = Encoder::new()?;
        encoder.u8(ROOT_CASE)?;
        encode_case_body(&mut encoder, self)?;
        Ok(encoder.finish())
    }

    /// Decode only an exact current bounded transport.
    pub fn from_canonical_bytes(bytes: &[u8]) -> Result<Self, VvCodecError> {
        let mut decoder = Decoder::new(bytes)?;
        expect_root(&mut decoder, ROOT_CASE, "case")?;
        let case = decode_case_body(&mut decoder)?;
        decoder.finish()?;
        if case.canonical_bytes()?.as_slice() != bytes {
            return Err(VvCodecError::at(
                0,
                "transport is structurally valid but not a canonical fixed point",
            ));
        }
        case.validate()
            .map_err(|error| model_error(0, "V&V case", error))?;
        Ok(case)
    }

    /// Domain-separated content identity of the exact canonical transport.
    pub fn content_hash(&self) -> Result<ContentHash, VvCodecError> {
        Ok(case_content_hash_for(&self.canonical_bytes()?))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_hash(label: &str) -> ContentHash {
        fs_blake3::hash_domain(
            "org.frankensim.fs-evidence.vv-codec-test.v1",
            label.as_bytes(),
        )
    }

    fn split_ref(id: &str, content: &str) -> ArtifactRef {
        ArtifactRef::new(
            ArtifactKind::CalibrationSplit,
            ArtifactId::try_new(id).expect("test split id"),
            test_hash(content),
        )
    }

    fn encoded_blind_selection(
        selection_split: &ArtifactRef,
        release: &BlindReleaseReceipt,
    ) -> Vec<u8> {
        encoded_blind_selection_parts(
            selection_split,
            release.split(),
            release.blind_commitment(),
            release.authority_receipt_hash(),
        )
    }

    fn encoded_blind_selection_parts(
        selection_split: &ArtifactRef,
        release_split: &ArtifactRef,
        blind_commitment: ContentHash,
        authority_receipt_hash: ContentHash,
    ) -> Vec<u8> {
        // Write the release fields directly so the negative cases can express
        // transports that the safe model constructors deliberately refuse.
        let mut encoder = Encoder::new().expect("selection encoder");
        encode_artifact_ref(&mut encoder, selection_split).expect("selection split");
        encoder.count(1).expect("one selected row");
        encode_observation_id(
            &mut encoder,
            &ObservationId::try_new("blind-1").expect("test observation id"),
        )
        .expect("selected row");
        encoder.u8(1).expect("blind partition tag");
        encode_artifact_ref(&mut encoder, release_split).expect("release split");
        encoder.hash(blind_commitment).expect("blind commitment");
        encoder
            .hash(authority_receipt_hash)
            .expect("release authority receipt");
        encoder.finish()
    }

    fn encoded_validation_selection(selection_split: &ArtifactRef) -> Vec<u8> {
        let mut encoder = Encoder::new().expect("selection encoder");
        encode_artifact_ref(&mut encoder, selection_split).expect("selection split");
        encoder.count(1).expect("one selected row");
        encode_observation_id(
            &mut encoder,
            &ObservationId::try_new("validation-1").expect("test observation id"),
        )
        .expect("selected row");
        encoder.u8(0).expect("validation partition tag");
        encoder.finish()
    }

    #[test]
    fn blind_selection_decoder_rejects_crosswired_release_and_binds_release_identity() {
        let outer_split = split_ref("split-1", "split-content-a");
        for (mismatch, release_split) in [
            (
                "same id but different content hash",
                split_ref("split-1", "split-content-b"),
            ),
            (
                "different id but same content hash",
                split_ref("split-2", "split-content-a"),
            ),
        ] {
            let forged = encoded_blind_selection_parts(
                &outer_split,
                &release_split,
                test_hash("blind-commitment"),
                test_hash("release-authority"),
            );
            let mut decoder = Decoder::new(&forged).expect("forged selection prefix");
            let error = decode_observation_selection(&mut decoder)
                .expect_err("decoder must refuse a crosswired release split");
            assert_eq!(error.rule_name(), CANONICAL_RULE, "{mismatch}");
            assert!(
                error.detail().contains("selection.blind_release"),
                "{mismatch}: {error}",
            );
            assert!(
                error.detail().contains("exact kind, id, and content hash"),
                "{mismatch}: {error}",
            );
        }

        let wrong_kind = ArtifactRef::new(
            ArtifactKind::ContextOfUse,
            outer_split.id().clone(),
            outer_split.hash(),
        );
        let forged = encoded_blind_selection_parts(
            &outer_split,
            &wrong_kind,
            test_hash("blind-commitment"),
            test_hash("release-authority"),
        );
        let mut decoder = Decoder::new(&forged).expect("wrong-kind selection prefix");
        let error = decode_observation_selection(&mut decoder)
            .expect_err("blind-release decoder must refuse a non-split release authority");
        assert_eq!(error.rule_name(), CANONICAL_RULE);
        assert!(error.detail().contains("blind-release receipt"));

        let release = BlindReleaseReceipt::new(
            outer_split.clone(),
            test_hash("blind-commitment"),
            test_hash("release-authority-a"),
        )
        .expect("exact release");
        let exact = encoded_blind_selection(&outer_split, &release);
        let mut decoder = Decoder::new(&exact).expect("exact selection prefix");
        let decoded = decode_observation_selection(&mut decoder)
            .expect("exact release and selection split must decode");
        decoder.finish().expect("exact selection transport");
        assert_eq!(decoded.split(), &outer_split);
        let expected_partition = EvidencePartition::BlindHoldout {
            release: release.clone(),
        };
        assert_eq!(decoded.partition(), &expected_partition);

        let other_authority = BlindReleaseReceipt::new(
            outer_split.clone(),
            release.blind_commitment(),
            test_hash("release-authority-b"),
        )
        .expect("other release authority");
        assert_ne!(
            encoded_blind_selection(&outer_split, &other_authority),
            exact,
            "authority receipt hash must move canonical selection bytes",
        );

        let other_commitment = BlindReleaseReceipt::new(
            outer_split.clone(),
            test_hash("other-blind-commitment"),
            release.authority_receipt_hash(),
        )
        .expect("other blind commitment");
        assert_ne!(
            encoded_blind_selection(&outer_split, &other_commitment),
            exact,
            "blind commitment must move canonical selection bytes",
        );
    }

    #[test]
    fn blind_release_refuses_zero_authority_fields_in_constructors_and_raw_transport() {
        let zero = ContentHash::from_slice(&[0_u8; 32]).expect("zero digest has the fixed width");
        let exact_split = split_ref("split-1", "split-content");
        let zero_split = ArtifactRef::new(
            ArtifactKind::CalibrationSplit,
            exact_split.id().clone(),
            zero,
        );
        let commitment = test_hash("blind-commitment");
        let authority = test_hash("release-authority");

        for (field, split, candidate_commitment, candidate_authority) in [
            (
                "split content hash",
                zero_split.clone(),
                commitment,
                authority,
            ),
            ("blind commitment", exact_split.clone(), zero, authority),
            ("authority receipt", exact_split.clone(), commitment, zero),
        ] {
            let error =
                BlindReleaseReceipt::new(split.clone(), candidate_commitment, candidate_authority)
                    .expect_err("zero authority material must not construct a blind release");
            assert!(
                error.to_string().contains("non-zero exact split identity"),
                "{field}: {error}",
            );

            let forged = encoded_blind_selection_parts(
                &split,
                &split,
                candidate_commitment,
                candidate_authority,
            );
            let mut decoder = Decoder::new(&forged).expect("forged selection prefix");
            let error = decode_observation_selection(&mut decoder)
                .expect_err("zero-but-internally-equal release transport must refuse");
            assert_eq!(error.rule_name(), CANONICAL_RULE, "{field}");
            assert!(
                error.detail().contains("blind-release receipt"),
                "{field}: {error}",
            );
        }
    }

    #[test]
    fn validation_selection_refuses_zero_split_capabilities_in_raw_transport() {
        let zero = ContentHash::from_slice(&[0_u8; 32]).expect("zero digest has the fixed width");
        let zero_split = ArtifactRef::new(
            ArtifactKind::CalibrationSplit,
            ArtifactId::try_new("split-1").expect("test split id"),
            zero,
        );
        let forged = encoded_validation_selection(&zero_split);
        let mut decoder = Decoder::new(&forged).expect("forged selection prefix");
        let error = decode_observation_selection(&mut decoder)
            .expect_err("a zero split digest cannot become a typed validation capability");
        assert_eq!(error.rule_name(), CANONICAL_RULE);
        assert!(
            error.detail().contains("non-zero exact split identity"),
            "{error}",
        );
    }

    #[test]
    fn canonical_transport_byte_limit_is_exact_for_encoder_and_decoder() {
        const PREFIX_BYTES: usize =
            MAGIC.len() + core::mem::size_of::<u32>() + core::mem::size_of::<u32>();
        let exact_payload_len = MAX_VV_CANONICAL_BYTES - PREFIX_BYTES;

        let mut exact_encoder = Encoder::new().expect("canonical encoder prefix fits");
        exact_encoder
            .raw(&vec![0; exact_payload_len])
            .expect("an exact-limit transport must encode");
        assert_eq!(exact_encoder.finish().len(), MAX_VV_CANONICAL_BYTES);

        let mut oversized_encoder = Encoder::new().expect("canonical encoder prefix fits");
        let encode_error = oversized_encoder
            .raw(&vec![0; exact_payload_len + 1])
            .expect_err("one byte beyond the canonical envelope must refuse");
        assert_eq!(encode_error.offset(), PREFIX_BYTES);
        assert_eq!(encode_error.rule_name(), CANONICAL_RULE);
        assert_eq!(
            encode_error.detail(),
            format!(
                "canonical transport would require {} bytes; limit is {MAX_VV_CANONICAL_BYTES}",
                MAX_VV_CANONICAL_BYTES + 1
            )
        );

        let mut exact_transport = vec![0; MAX_VV_CANONICAL_BYTES];
        exact_transport[..MAGIC.len()].copy_from_slice(MAGIC);
        exact_transport[4..8].copy_from_slice(&VV_SCHEMA_VERSION.to_le_bytes());
        exact_transport[8..12].copy_from_slice(&VV_RULESET_VERSION.to_le_bytes());
        let exact_decoder =
            Decoder::new(&exact_transport).expect("the exact byte limit passes the size gate");
        assert_eq!(exact_decoder.position(), PREFIX_BYTES);

        let mut oversized_transport = exact_transport;
        oversized_transport.push(0);
        let decode_error = Decoder::new(&oversized_transport)
            .err()
            .expect("one byte beyond the canonical envelope must refuse before decoding");
        assert_eq!(decode_error.offset(), 0);
        assert_eq!(decode_error.rule_name(), CANONICAL_RULE);
        assert_eq!(
            decode_error.detail(),
            format!(
                "canonical transport has {} bytes; limit is {MAX_VV_CANONICAL_BYTES}",
                MAX_VV_CANONICAL_BYTES + 1
            )
        );
    }

    #[test]
    fn collection_and_nested_aggregate_limits_refuse_before_allocation() {
        assert_eq!(
            MAX_VV_TOTAL_COLLECTION_ITEMS % MAX_VV_COLLECTION_ITEMS,
            0,
            "the exact-boundary construction assumes whole maximum-size collections",
        );
        let full_collections = MAX_VV_TOTAL_COLLECTION_ITEMS / MAX_VV_COLLECTION_ITEMS;
        assert_eq!(full_collections, 16);

        let mut per_collection = Encoder::new().expect("canonical encoder prefix");
        let per_collection_offset = per_collection.bytes.len();
        let error = per_collection
            .count(MAX_VV_COLLECTION_ITEMS + 1)
            .expect_err("one entry beyond the per-collection limit must refuse");
        assert_eq!(error.offset(), per_collection_offset);
        assert_eq!(
            error.detail(),
            format!(
                "collection has {} entries; limit is {MAX_VV_COLLECTION_ITEMS}",
                MAX_VV_COLLECTION_ITEMS + 1
            )
        );
        assert_eq!(
            per_collection.bytes.len(),
            per_collection_offset,
            "the encoder must refuse before writing an inadmissible count",
        );

        let mut exact_aggregate = Encoder::new().expect("canonical encoder prefix");
        for _ in 0..full_collections {
            exact_aggregate
                .count(MAX_VV_COLLECTION_ITEMS)
                .expect("the exact aggregate item limit must encode");
        }
        assert_eq!(
            exact_aggregate.collection_items,
            MAX_VV_TOTAL_COLLECTION_ITEMS
        );
        let aggregate_offset = exact_aggregate.bytes.len();
        assert_eq!(aggregate_offset, 12 + full_collections * 8);
        let error = exact_aggregate
            .count(1)
            .expect_err("one nested item beyond the aggregate limit must refuse");
        assert_eq!(error.offset(), aggregate_offset);
        assert_eq!(
            error.detail(),
            format!("aggregate collection entries exceed {MAX_VV_TOTAL_COLLECTION_ITEMS}")
        );
        assert_eq!(
            exact_aggregate.bytes.len(),
            aggregate_offset,
            "the encoder must refuse before writing the over-budget nested count",
        );

        let mut oversized_count_transport = Encoder::new().expect("canonical encoder prefix");
        oversized_count_transport
            .u64((MAX_VV_COLLECTION_ITEMS + 1) as u64)
            .expect("raw adversarial count word fits the byte envelope");
        let oversized_count_bytes = oversized_count_transport.finish();
        let mut decoder = Decoder::new(&oversized_count_bytes)
            .expect("adversarial count transport has a valid prefix");
        let error = decoder
            .count()
            .expect_err("decoder must gate a count before bounded_vec can allocate");
        assert_eq!(error.offset(), 12);
        assert_eq!(
            error.detail(),
            format!(
                "collection length {} exceeds {MAX_VV_COLLECTION_ITEMS}",
                MAX_VV_COLLECTION_ITEMS + 1
            )
        );

        let mut nested_transport = Encoder::new().expect("canonical encoder prefix");
        for _ in 0..full_collections {
            nested_transport
                .u64(MAX_VV_COLLECTION_ITEMS as u64)
                .expect("raw exact-limit count word fits");
        }
        nested_transport.u64(1).expect("raw N+1 count word fits");
        let nested_bytes = nested_transport.finish();
        let mut decoder = Decoder::new(&nested_bytes).expect("nested-count transport prefix");
        for _ in 0..full_collections {
            assert_eq!(
                decoder.count().expect("exact aggregate count must decode"),
                MAX_VV_COLLECTION_ITEMS
            );
        }
        assert_eq!(decoder.collection_items, MAX_VV_TOTAL_COLLECTION_ITEMS);
        assert_eq!(decoder.position(), aggregate_offset);
        let error = decoder
            .count()
            .expect_err("decoder must refuse the aggregate N+1 count before allocation");
        assert_eq!(error.offset(), aggregate_offset);
        assert_eq!(
            error.detail(),
            format!("aggregate collection entries exceed {MAX_VV_TOTAL_COLLECTION_ITEMS}")
        );
    }
}
