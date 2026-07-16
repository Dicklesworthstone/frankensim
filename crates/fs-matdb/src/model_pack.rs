//! Canonical runtime transport for admitted constitutive-model cards.
//!
//! [`NormalizedModelPack`] is deliberately separate from the scalar/curve
//! [`crate::NormalizedPack`] wire format. A model card is not a property claim:
//! it carries a law identity, a dimensioned parameter block, state semantics,
//! validity, and source provenance that downstream law adapters validate again.
//! Raw NASA tables, kinetics files, licenses, and unit expressions remain an
//! offline-compiler concern.

use std::collections::{BTreeMap, BTreeSet};

use fs_blake3::{ContentHash, hash_domain};
use fs_evidence::ValidityDomain;
use fs_qty::Dims;

use crate::{
    ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, PackError, Provenance,
    ValidityBoundSide,
};

/// Current normalized model-card pack wire schema.
pub const MODEL_PACK_SCHEMA_VERSION: u32 = 1;
/// Canonical coherent-SI target basis and explicit dimension order.
pub const MODEL_PACK_TARGET_BASIS: &str = "SI-six-base[m,kg,s,K,A,mol]";

const MAGIC: &[u8; 8] = b"FSMODPK\0";
const MODEL_PACK_HASH_DOMAIN: &str = "org.frankensim.fs-matdb.normalized-model-pack.v1";
const MAX_PACK_BYTES: usize = 256 * 1024 * 1024;
const MAX_STRING_BYTES: usize = 1_048_576;
const MAX_MODELS: usize = 100_000;
const MAX_PARAMETERS: usize = 100_000;
const MAX_PARAMETERS_PER_MODEL: usize = 4_096;
const MAX_VALIDITY_AXES: usize = 100_000;
const MAX_VALIDITY_AXES_PER_MODEL: usize = 4_096;
const MAX_SOURCES: usize = 1_000_000;
const MAX_SOURCES_PER_MODEL: usize = 4_096;
const MAX_NORMALIZATIONS: usize = 300_000;
const MIN_MODEL_BYTES: usize = 49;
const MIN_PARAMETER_BYTES: usize = 18;
const MIN_VALIDITY_BYTES: usize = 20;
const MIN_NORMALIZATION_BYTES: usize = 96;

/// Exact normalized numeric field inside one constitutive-model card.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum ModelNormalizationTarget {
    /// One named law parameter.
    Parameter {
        /// Canonical identity of the reconstructed model card.
        model: ContentHash,
        /// Parameter name inside the card's ordered parameter block.
        parameter: String,
    },
    /// One endpoint of a named model-validity interval.
    ValidityBound {
        /// Canonical identity of the reconstructed model card.
        model: ContentHash,
        /// Existing validity-axis name.
        axis: String,
        /// Lower or upper endpoint.
        side: ValidityBoundSide,
    },
}

/// Auditable unit/basis transform applied to one model-card numeric field.
#[derive(Debug, Clone, PartialEq)]
pub struct ModelNormalizationReceipt {
    target: ModelNormalizationTarget,
    source_literal: ContentHash,
    dims: Dims,
    scale: f64,
    offset: f64,
    source_basis: String,
    target_basis: String,
    source_frame: Option<String>,
    target_frame: Option<String>,
}

impl ModelNormalizationReceipt {
    /// Construct an immutable model-field normalization receipt.
    #[allow(clippy::too_many_arguments)]
    #[must_use]
    pub fn new(
        target: ModelNormalizationTarget,
        source_literal: ContentHash,
        dims: Dims,
        scale: f64,
        offset: f64,
        source_basis: impl Into<String>,
        target_basis: impl Into<String>,
        source_frame: Option<String>,
        target_frame: Option<String>,
    ) -> Self {
        Self {
            target,
            source_literal,
            dims,
            scale,
            offset,
            source_basis: source_basis.into(),
            target_basis: target_basis.into(),
            source_frame,
            target_frame,
        }
    }

    /// Structurally linked model field.
    #[must_use]
    pub fn target(&self) -> &ModelNormalizationTarget {
        &self.target
    }

    /// Hash of the exact source literal and source-unit token.
    #[must_use]
    pub fn source_literal(&self) -> ContentHash {
        self.source_literal
    }

    /// Six-base dimensions of the normalized field.
    #[must_use]
    pub fn dims(&self) -> Dims {
        self.dims
    }

    /// Multiplicative term in `si = source * scale + offset`.
    #[must_use]
    pub fn scale(&self) -> f64 {
        self.scale
    }

    /// Additive term in `si = source * scale + offset`.
    #[must_use]
    pub fn offset(&self) -> f64 {
        self.offset
    }

    /// Explicit source unit/basis expression.
    #[must_use]
    pub fn source_basis(&self) -> &str {
        &self.source_basis
    }

    /// Canonical normalized basis.
    #[must_use]
    pub fn target_basis(&self) -> &str {
        &self.target_basis
    }

    /// Optional source frame identifier.
    #[must_use]
    pub fn source_frame(&self) -> Option<&str> {
        self.source_frame.as_deref()
    }

    /// Optional target frame identifier.
    #[must_use]
    pub fn target_frame(&self) -> Option<&str> {
        self.target_frame.as_deref()
    }
}

/// Runtime-loadable result of an admitted offline model-card compilation.
#[derive(Debug, Clone, PartialEq)]
pub struct NormalizedModelPack {
    pack_id: String,
    compiler: String,
    source_artifact: ContentHash,
    redistribution_terms: String,
    models: Vec<ConstitutiveModelCard>,
    normalizations: Vec<ModelNormalizationReceipt>,
}

impl NormalizedModelPack {
    /// Admit cards and complete normalization evidence into the portable profile.
    ///
    /// Cards are sorted by their full semantic content identity. Every parameter
    /// and both endpoints of every validity interval require exactly one receipt;
    /// missing, duplicate, dangling, or dimension-inconsistent receipts refuse.
    pub fn new(
        pack_id: impl Into<String>,
        compiler: impl Into<String>,
        source_artifact: ContentHash,
        redistribution_terms: impl Into<String>,
        mut models: Vec<ConstitutiveModelCard>,
        mut normalizations: Vec<ModelNormalizationReceipt>,
    ) -> Result<Self, PackError> {
        let pack_id = pack_id.into();
        let compiler = compiler.into();
        let redistribution_terms = redistribution_terms.into();
        require_text("pack_id", &pack_id)?;
        require_text("compiler", &compiler)?;
        require_text("redistribution_terms", &redistribution_terms)?;
        if models.is_empty() {
            return Err(invalid(
                "models",
                "a normalized model pack requires at least one card",
            ));
        }
        if models.len() > MAX_MODELS {
            return Err(limit("models", MAX_MODELS, models.len()));
        }
        if normalizations.len() > MAX_NORMALIZATIONS {
            return Err(limit(
                "model_normalizations",
                MAX_NORMALIZATIONS,
                normalizations.len(),
            ));
        }
        let mut total_parameters = 0usize;
        let mut total_validity_axes = 0usize;
        let mut total_sources = 0usize;
        for model in &models {
            checked_total(
                &mut total_parameters,
                model.parameters.len(),
                MAX_PARAMETERS,
                "model_parameters",
            )?;
            checked_total(
                &mut total_validity_axes,
                model.validity.bounds().len(),
                MAX_VALIDITY_AXES,
                "model_validity_axes",
            )?;
            checked_total(
                &mut total_sources,
                model.sources.len(),
                MAX_SOURCES,
                "model_sources",
            )?;
            validate_model(model)?;
        }
        estimate_pack_bytes(
            &pack_id,
            &compiler,
            &redistribution_terms,
            &models,
            &normalizations,
        )?;
        models.sort_by_cached_key(ConstitutiveModelCard::content_hash);
        for pair in models.windows(2) {
            if pair[0].content_hash() == pair[1].content_hash() {
                return Err(invalid(
                    "models",
                    format!("duplicate model card {}", pair[0].content_hash()),
                ));
            }
        }
        let models_by_hash: BTreeMap<_, _> = models
            .iter()
            .map(|model| (model.content_hash(), model))
            .collect();

        normalizations.sort_by(|left, right| left.target.cmp(&right.target));
        for pair in normalizations.windows(2) {
            if pair[0].target == pair[1].target {
                return Err(invalid(
                    "model_normalizations",
                    format!("duplicate target receipt {:?}", pair[0].target),
                ));
            }
        }
        for receipt in &normalizations {
            validate_receipt(&models_by_hash, receipt)?;
        }
        validate_complete_receipt_coverage(&models, &normalizations)?;
        validate_validity_receipt_coherence(&normalizations)?;

        Ok(Self {
            pack_id,
            compiler,
            source_artifact,
            redistribution_terms,
            models,
            normalizations,
        })
    }

    /// Stable pack name supplied by the source manifest.
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
    pub fn source_artifact(&self) -> ContentHash {
        self.source_artifact
    }

    /// Retained redistribution decision/terms.
    #[must_use]
    pub fn redistribution_terms(&self) -> &str {
        &self.redistribution_terms
    }

    /// Canonically ordered runtime model cards.
    #[must_use]
    pub fn models(&self) -> &[ConstitutiveModelCard] {
        &self.models
    }

    /// Canonically ordered, complete field-normalization receipts.
    #[must_use]
    pub fn normalizations(&self) -> &[ModelNormalizationReceipt] {
        &self.normalizations
    }

    /// Canonical binary representation consumed by L1.
    #[must_use]
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut writer = Writer::default();
        writer.bytes.extend_from_slice(MAGIC);
        writer.u32(MODEL_PACK_SCHEMA_VERSION);
        writer.string(&self.pack_id);
        writer.string(&self.compiler);
        writer.hash(self.source_artifact);
        writer.string(&self.redistribution_terms);
        writer.count(self.models.len());
        for model in &self.models {
            writer.hash(model.content_hash());
            encode_model(&mut writer, model);
        }
        writer.count(self.normalizations.len());
        for receipt in &self.normalizations {
            encode_target(&mut writer, &receipt.target);
            writer.hash(receipt.source_literal);
            writer.dims(receipt.dims);
            writer.f64(receipt.scale);
            writer.f64(receipt.offset);
            writer.string(&receipt.source_basis);
            writer.string(&receipt.target_basis);
            writer.optional_string(receipt.source_frame.as_deref());
            writer.optional_string(receipt.target_frame.as_deref());
        }
        writer.bytes
    }

    /// Domain-separated identity of the canonical pack bytes.
    #[must_use]
    pub fn content_hash(&self) -> ContentHash {
        hash_domain(MODEL_PACK_HASH_DOMAIN, &self.to_bytes())
    }

    /// Verify an externally pinned whole-pack identity before decoding.
    pub fn from_bytes_verified(expected: ContentHash, bytes: &[u8]) -> Result<Self, PackError> {
        if bytes.len() > MAX_PACK_BYTES {
            return Err(limit("model_pack_bytes", MAX_PACK_BYTES, bytes.len()));
        }
        let actual = hash_domain(MODEL_PACK_HASH_DOMAIN, bytes);
        if actual != expected {
            return Err(PackError::IdentityMismatch {
                kind: "model pack",
                expected,
                actual,
            });
        }
        Self::from_bytes(bytes)
    }

    /// Decode and semantically re-admit one canonical model pack.
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, PackError> {
        if bytes.len() > MAX_PACK_BYTES {
            return Err(limit("model_pack_bytes", MAX_PACK_BYTES, bytes.len()));
        }
        let mut reader = Reader::new(bytes);
        reader.expect(MAGIC, "normalized model-pack magic")?;
        let version = reader.u32()?;
        if version != MODEL_PACK_SCHEMA_VERSION {
            return Err(reader.malformed(format!(
                "unsupported schema version {version}; expected {MODEL_PACK_SCHEMA_VERSION}"
            )));
        }
        let pack_id = reader.string()?;
        let compiler = reader.string()?;
        let source_artifact = reader.hash()?;
        let redistribution_terms = reader.string()?;

        let model_count = reader.count("models", MAX_MODELS)?;
        reader.require_items(model_count, MIN_MODEL_BYTES, "model cards")?;
        let mut total_parameters = 0usize;
        let mut total_validity_axes = 0usize;
        let mut total_sources = 0usize;
        let mut models = Vec::with_capacity(model_count);
        for _ in 0..model_count {
            let expected = reader.hash()?;
            let model = decode_model(
                &mut reader,
                &mut total_parameters,
                &mut total_validity_axes,
                &mut total_sources,
            )?;
            let actual = model.content_hash();
            if actual != expected {
                return Err(PackError::IdentityMismatch {
                    kind: "model card",
                    expected,
                    actual,
                });
            }
            models.push(model);
        }

        let normalization_count = reader.count("model_normalizations", MAX_NORMALIZATIONS)?;
        reader.require_items(
            normalization_count,
            MIN_NORMALIZATION_BYTES,
            "model normalization receipts",
        )?;
        let mut normalizations = Vec::with_capacity(normalization_count);
        for _ in 0..normalization_count {
            normalizations.push(ModelNormalizationReceipt::new(
                decode_target(&mut reader)?,
                reader.hash()?,
                reader.dims()?,
                reader.f64()?,
                reader.f64()?,
                reader.string()?,
                reader.string()?,
                reader.optional_string()?,
                reader.optional_string()?,
            ));
        }
        reader.finish()?;
        let pack = Self::new(
            pack_id,
            compiler,
            source_artifact,
            redistribution_terms,
            models,
            normalizations,
        )?;
        if pack.to_bytes() != bytes {
            return Err(PackError::Malformed {
                at: bytes.len(),
                detail: "decoded fields do not reproduce the canonical model-pack byte stream"
                    .to_string(),
            });
        }
        Ok(pack)
    }
}

fn validate_model(model: &ConstitutiveModelCard) -> Result<(), PackError> {
    require_text("model.law", &model.law.0)?;
    if model.law_version == 0 {
        return Err(invalid(
            "model.law_version",
            "portable law versions start at one",
        ));
    }
    if model.parameters.len() > MAX_PARAMETERS_PER_MODEL {
        return Err(limit(
            "model_parameters",
            MAX_PARAMETERS_PER_MODEL,
            model.parameters.len(),
        ));
    }
    for (name, parameter) in &model.parameters {
        require_text("model.parameter.name", name)?;
        require_portable_f64("model.parameter.value", parameter.value)?;
    }
    if model.validity.bounds().len() > MAX_VALIDITY_AXES_PER_MODEL {
        return Err(limit(
            "model_validity_axes",
            MAX_VALIDITY_AXES_PER_MODEL,
            model.validity.bounds().len(),
        ));
    }
    for (axis, &(lower, upper)) in model.validity.bounds() {
        require_text("model.validity.axis", axis)?;
        require_portable_f64("model.validity.lower", lower)?;
        require_portable_f64("model.validity.upper", upper)?;
        if lower > upper {
            return Err(invalid(
                "model.validity",
                format!("axis {axis:?} has reversed normalized bounds"),
            ));
        }
    }
    if model.sources.is_empty() {
        return Err(invalid(
            "model.sources",
            "portable model cards require at least one source artifact",
        ));
    }
    if model.sources.len() > MAX_SOURCES_PER_MODEL {
        return Err(limit(
            "model_sources",
            MAX_SOURCES_PER_MODEL,
            model.sources.len(),
        ));
    }
    if !model.sources.windows(2).all(|pair| pair[0] < pair[1]) {
        return Err(invalid(
            "model.sources",
            "source hashes must be strictly increasing and deduplicated",
        ));
    }
    require_text("model.provenance.source", &model.provenance.source)?;
    require_text("model.provenance.license", &model.provenance.license)?;
    if model.provenance.artifact.is_none() {
        return Err(invalid(
            "model.provenance.artifact",
            "portable model cards require an exact source artifact hash",
        ));
    }
    model.validate()?;
    Ok(())
}

fn validate_receipt(
    models: &BTreeMap<ContentHash, &ConstitutiveModelCard>,
    receipt: &ModelNormalizationReceipt,
) -> Result<(), PackError> {
    require_portable_f64("model_normalization.scale", receipt.scale)?;
    if receipt.scale <= 0.0 {
        return Err(invalid(
            "model_normalization.scale",
            "normalization scale must be positive",
        ));
    }
    require_portable_f64("model_normalization.offset", receipt.offset)?;
    require_text("model_normalization.source_basis", &receipt.source_basis)?;
    if receipt.target_basis != MODEL_PACK_TARGET_BASIS {
        return Err(invalid(
            "model_normalization.target_basis",
            format!("target basis must be {MODEL_PACK_TARGET_BASIS:?}"),
        ));
    }
    match (&receipt.source_frame, &receipt.target_frame) {
        (None, None) => {}
        (Some(source), Some(target)) => {
            require_text("model_normalization.source_frame", source)?;
            require_text("model_normalization.target_frame", target)?;
        }
        _ => {
            return Err(invalid(
                "model_normalization.frame",
                "source and target frames must be both absent or both present",
            ));
        }
    }

    match &receipt.target {
        ModelNormalizationTarget::Parameter { model, parameter } => {
            require_text("model_normalization.parameter", parameter)?;
            let card = find_model(models, *model)?;
            let expected = card.parameters.get(parameter).ok_or_else(|| {
                invalid(
                    "model_normalization.target",
                    format!("model {model} has no parameter {parameter:?}"),
                )
            })?;
            if expected.dims != receipt.dims {
                return Err(invalid(
                    "model_normalization.dims",
                    format!(
                        "parameter {parameter:?} has dimensions {:?}, receipt declares {:?}",
                        expected.dims, receipt.dims
                    ),
                ));
            }
        }
        ModelNormalizationTarget::ValidityBound { model, axis, .. } => {
            require_text("model_normalization.validity_axis", axis)?;
            let card = find_model(models, *model)?;
            if card.validity.bound(axis).is_none() {
                return Err(invalid(
                    "model_normalization.target",
                    format!("model {model} has no validity axis {axis:?}"),
                ));
            }
        }
    }
    Ok(())
}

fn find_model<'a>(
    models: &BTreeMap<ContentHash, &'a ConstitutiveModelCard>,
    identity: ContentHash,
) -> Result<&'a ConstitutiveModelCard, PackError> {
    models.get(&identity).copied().ok_or_else(|| {
        invalid(
            "model_normalization.target",
            format!("unknown model card {identity}"),
        )
    })
}

fn validate_complete_receipt_coverage(
    models: &[ConstitutiveModelCard],
    receipts: &[ModelNormalizationReceipt],
) -> Result<(), PackError> {
    let mut expected = BTreeSet::new();
    for model in models {
        let identity = model.content_hash();
        for parameter in model.parameters.keys() {
            expected.insert(ModelNormalizationTarget::Parameter {
                model: identity,
                parameter: parameter.clone(),
            });
        }
        for axis in model.validity.bounds().keys() {
            for side in [ValidityBoundSide::Lower, ValidityBoundSide::Upper] {
                expected.insert(ModelNormalizationTarget::ValidityBound {
                    model: identity,
                    axis: axis.clone(),
                    side,
                });
            }
        }
    }
    if receipts.len() != expected.len() {
        return Err(invalid(
            "model_normalizations",
            format!(
                "receipts must cover every normalized field exactly once; expected {}, found {}",
                expected.len(),
                receipts.len()
            ),
        ));
    }
    if let Some((actual, expected)) = receipts
        .iter()
        .map(|receipt| &receipt.target)
        .zip(expected.iter())
        .find(|(actual, expected)| actual != expected)
    {
        return Err(invalid(
            "model_normalizations",
            format!("receipt target mismatch; expected {expected:?}, found {actual:?}"),
        ));
    }
    Ok(())
}

fn validate_validity_receipt_coherence(
    receipts: &[ModelNormalizationReceipt],
) -> Result<(), PackError> {
    let mut lower_by_axis: BTreeMap<(ContentHash, String), &ModelNormalizationReceipt> =
        BTreeMap::new();
    for receipt in receipts {
        let ModelNormalizationTarget::ValidityBound { model, axis, side } = &receipt.target else {
            continue;
        };
        match side {
            ValidityBoundSide::Lower => {
                lower_by_axis.insert((*model, axis.clone()), receipt);
            }
            ValidityBoundSide::Upper => {
                let lower = lower_by_axis.get(&(*model, axis.clone())).ok_or_else(|| {
                    invalid(
                        "model_normalizations",
                        format!("validity axis {axis:?} is missing its lower receipt"),
                    )
                })?;
                if lower.dims != receipt.dims
                    || lower.scale.to_bits() != receipt.scale.to_bits()
                    || lower.offset.to_bits() != receipt.offset.to_bits()
                    || lower.source_basis != receipt.source_basis
                    || lower.target_basis != receipt.target_basis
                    || lower.source_frame != receipt.source_frame
                    || lower.target_frame != receipt.target_frame
                {
                    return Err(invalid(
                        "model_normalizations",
                        format!("validity axis {axis:?} endpoints must share one exact transform"),
                    ));
                }
            }
        }
    }
    Ok(())
}

fn estimate_pack_bytes(
    pack_id: &str,
    compiler: &str,
    redistribution_terms: &str,
    models: &[ConstitutiveModelCard],
    receipts: &[ModelNormalizationReceipt],
) -> Result<usize, PackError> {
    let mut bytes = MAGIC.len() + 4;
    add_string_size(&mut bytes, pack_id)?;
    add_string_size(&mut bytes, compiler)?;
    add_size(&mut bytes, 32)?;
    add_string_size(&mut bytes, redistribution_terms)?;
    add_size(&mut bytes, 4)?;
    for model in models {
        add_size(&mut bytes, 32)?;
        add_string_size(&mut bytes, &model.law.0)?;
        add_size(&mut bytes, 4 + 4)?;
        for name in model.parameters.keys() {
            add_string_size(&mut bytes, name)?;
            add_size(&mut bytes, 8 + 6)?;
        }
        add_size(&mut bytes, 4 + 1 + 4)?;
        for axis in model.validity.bounds().keys() {
            add_string_size(&mut bytes, axis)?;
            add_size(&mut bytes, 16)?;
        }
        add_size(&mut bytes, 4)?;
        add_size(&mut bytes, model.sources.len().saturating_mul(32))?;
        add_string_size(&mut bytes, &model.provenance.source)?;
        add_string_size(&mut bytes, &model.provenance.license)?;
        add_size(
            &mut bytes,
            if model.provenance.artifact.is_some() {
                1 + 32
            } else {
                1
            },
        )?;
    }
    add_size(&mut bytes, 4)?;
    for receipt in receipts {
        match &receipt.target {
            ModelNormalizationTarget::Parameter { parameter, .. } => {
                add_size(&mut bytes, 1 + 32)?;
                add_string_size(&mut bytes, parameter)?;
            }
            ModelNormalizationTarget::ValidityBound { axis, .. } => {
                add_size(&mut bytes, 1 + 32)?;
                add_string_size(&mut bytes, axis)?;
                add_size(&mut bytes, 1)?;
            }
        }
        add_size(&mut bytes, 32 + 6 + 8 + 8)?;
        add_string_size(&mut bytes, &receipt.source_basis)?;
        add_string_size(&mut bytes, &receipt.target_basis)?;
        add_optional_string_size(&mut bytes, receipt.source_frame.as_deref())?;
        add_optional_string_size(&mut bytes, receipt.target_frame.as_deref())?;
    }
    Ok(bytes)
}

fn add_optional_string_size(bytes: &mut usize, value: Option<&str>) -> Result<(), PackError> {
    add_size(bytes, 1)?;
    if let Some(value) = value {
        add_string_size(bytes, value)?;
    }
    Ok(())
}

fn add_string_size(bytes: &mut usize, value: &str) -> Result<(), PackError> {
    add_size(bytes, 4)?;
    add_size(bytes, value.len())
}

fn add_size(bytes: &mut usize, amount: usize) -> Result<(), PackError> {
    *bytes = bytes
        .checked_add(amount)
        .ok_or_else(|| limit("model_pack_bytes", MAX_PACK_BYTES, usize::MAX))?;
    if *bytes > MAX_PACK_BYTES {
        Err(limit("model_pack_bytes", MAX_PACK_BYTES, *bytes))
    } else {
        Ok(())
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

fn require_portable_f64(field: &'static str, value: f64) -> Result<(), PackError> {
    if !value.is_finite() {
        return Err(invalid(field, "value must be finite"));
    }
    if value.to_bits() == (-0.0f64).to_bits() {
        return Err(invalid(field, "negative zero is not canonical"));
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

fn encode_model(writer: &mut Writer, model: &ConstitutiveModelCard) {
    writer.string(&model.law.0);
    writer.u32(model.law_version);
    writer.count(model.parameters.len());
    for (name, parameter) in &model.parameters {
        writer.string(name);
        writer.f64(parameter.value);
        writer.dims(parameter.dims);
    }
    writer.u32(model.state_schema_version);
    writer.u8(match model.initial_state {
        InitialStatePolicy::ZeroInternalState => 0,
        InitialStatePolicy::RequiresDeclaredState => 1,
    });
    writer.count(model.validity.bounds().len());
    for (axis, &(lower, upper)) in model.validity.bounds() {
        writer.string(axis);
        writer.f64(lower);
        writer.f64(upper);
    }
    writer.count(model.sources.len());
    for source in &model.sources {
        writer.hash(*source);
    }
    encode_provenance(writer, &model.provenance);
}

fn decode_model(
    reader: &mut Reader<'_>,
    total_parameters: &mut usize,
    total_validity_axes: &mut usize,
    total_sources: &mut usize,
) -> Result<ConstitutiveModelCard, PackError> {
    let law = LawId(reader.string()?);
    let law_version = reader.u32()?;
    let parameter_count = reader.count("model_parameters", MAX_PARAMETERS_PER_MODEL)?;
    checked_total(
        total_parameters,
        parameter_count,
        MAX_PARAMETERS,
        "model_parameters",
    )?;
    reader.require_items(parameter_count, MIN_PARAMETER_BYTES, "model parameters")?;
    let mut parameters = BTreeMap::new();
    let mut previous_parameter: Option<String> = None;
    for _ in 0..parameter_count {
        let name = reader.string()?;
        if previous_parameter
            .as_ref()
            .is_some_and(|previous| previous >= &name)
        {
            return Err(
                reader.malformed("model parameters are not strictly increasing and deduplicated")
            );
        }
        let parameter = LawParameter {
            value: reader.f64()?,
            dims: reader.dims()?,
        };
        previous_parameter = Some(name.clone());
        parameters.insert(name, parameter);
    }
    let state_schema_version = reader.u32()?;
    let initial_state = match reader.u8()? {
        0 => InitialStatePolicy::ZeroInternalState,
        1 => InitialStatePolicy::RequiresDeclaredState,
        tag => return Err(reader.malformed(format!("unknown initial-state-policy tag {tag}"))),
    };
    let validity_count = reader.count("model_validity_axes", MAX_VALIDITY_AXES_PER_MODEL)?;
    checked_total(
        total_validity_axes,
        validity_count,
        MAX_VALIDITY_AXES,
        "model_validity_axes",
    )?;
    reader.require_items(validity_count, MIN_VALIDITY_BYTES, "model validity axes")?;
    let mut validity = ValidityDomain::unconstrained();
    let mut previous_axis: Option<String> = None;
    for _ in 0..validity_count {
        let axis = reader.string()?;
        if previous_axis
            .as_ref()
            .is_some_and(|previous| previous >= &axis)
        {
            return Err(reader
                .malformed("model validity axes are not strictly increasing and deduplicated"));
        }
        let lower = reader.f64()?;
        let upper = reader.f64()?;
        previous_axis = Some(axis.clone());
        validity = validity.with(axis, lower, upper);
    }
    let source_count = reader.count("model_sources", MAX_SOURCES_PER_MODEL)?;
    checked_total(total_sources, source_count, MAX_SOURCES, "model_sources")?;
    reader.require_items(source_count, 32, "model source hashes")?;
    let mut sources = Vec::with_capacity(source_count);
    for _ in 0..source_count {
        sources.push(reader.hash()?);
    }
    Ok(ConstitutiveModelCard {
        law,
        law_version,
        parameters,
        state_schema_version,
        initial_state,
        validity,
        sources,
        provenance: decode_provenance(reader)?,
    })
}

fn checked_total(
    total: &mut usize,
    add: usize,
    maximum: usize,
    resource: &'static str,
) -> Result<(), PackError> {
    *total = total
        .checked_add(add)
        .ok_or_else(|| limit(resource, maximum, usize::MAX))?;
    if *total > maximum {
        Err(limit(resource, maximum, *total))
    } else {
        Ok(())
    }
}

fn encode_target(writer: &mut Writer, target: &ModelNormalizationTarget) {
    match target {
        ModelNormalizationTarget::Parameter { model, parameter } => {
            writer.u8(0);
            writer.hash(*model);
            writer.string(parameter);
        }
        ModelNormalizationTarget::ValidityBound { model, axis, side } => {
            writer.u8(1);
            writer.hash(*model);
            writer.string(axis);
            writer.u8(match side {
                ValidityBoundSide::Lower => 0,
                ValidityBoundSide::Upper => 1,
            });
        }
    }
}

fn decode_target(reader: &mut Reader<'_>) -> Result<ModelNormalizationTarget, PackError> {
    match reader.u8()? {
        0 => Ok(ModelNormalizationTarget::Parameter {
            model: reader.hash()?,
            parameter: reader.string()?,
        }),
        1 => {
            let model = reader.hash()?;
            let axis = reader.string()?;
            let side = match reader.u8()? {
                0 => ValidityBoundSide::Lower,
                1 => ValidityBoundSide::Upper,
                tag => {
                    return Err(
                        reader.malformed(format!("unknown model-validity-bound-side tag {tag}"))
                    );
                }
            };
            Ok(ModelNormalizationTarget::ValidityBound { model, axis, side })
        }
        tag => Err(reader.malformed(format!("unknown model-normalization-target tag {tag}"))),
    }
}

fn encode_provenance(writer: &mut Writer, provenance: &Provenance) {
    writer.string(&provenance.source);
    writer.string(&provenance.license);
    match provenance.artifact {
        None => writer.u8(0),
        Some(hash) => {
            writer.u8(1);
            writer.hash(hash);
        }
    }
}

fn decode_provenance(reader: &mut Reader<'_>) -> Result<Provenance, PackError> {
    let source = reader.string()?;
    let license = reader.string()?;
    let artifact = match reader.u8()? {
        0 => None,
        1 => Some(reader.hash()?),
        tag => return Err(reader.malformed(format!("unknown provenance-artifact tag {tag}"))),
    };
    Ok(Provenance {
        source,
        license,
        artifact,
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
        self.bytes.extend_from_slice(value.as_bytes());
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
        let minimum = count
            .checked_mul(minimum_width)
            .ok_or_else(|| self.malformed(format!("{field} byte length overflow")))?;
        if self.bytes.len().saturating_sub(self.cursor) < minimum {
            Err(self.malformed(format!("truncated {field} needs at least {minimum} bytes")))
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
                    "{} trailing bytes after canonical model pack",
                    self.bytes.len() - self.cursor
                ),
            })
        }
    }
}
