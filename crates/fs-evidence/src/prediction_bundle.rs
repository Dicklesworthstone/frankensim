//! Generic canonical prediction-input bundles and atomic sealers
//! (bead frankensim-jmh21.1, input half).
//!
//! A blind prediction is only as trustworthy as the boundary between what
//! the predictor was GIVEN and what it later CLAIMS. This module owns the
//! given side: [`PredictionExecutionInput`] binds source/build/runtime
//! identity, typed scenario and card references, parameter-distribution
//! references, the random-stream design, the admitted model rungs and
//! applicability policy, QoI identities, the evidence role, the blind
//! partition, and an access policy — and structurally CANNOT carry target
//! labels, outcomes, or computed predictions: no such field exists, the
//! canonical transport rejects unknown fields, and the hostile battery
//! proves both.
//!
//! It builds on the V&V vocabulary by REFERENCE ([`ArtifactRef`] /
//! [`ContentHash`]) rather than duplicating [`ContextOfUse`],
//! [`ValidationPlan`], [`CalibrationSplit`], or the solution-verification
//! receipts; the referenced artifacts remain their own authorities.
//!
//! Sealing is atomic: canonical bytes are written to a temporary path and
//! atomically renamed, so an interrupted publication leaves no readable
//! bundle at the sealed path. A seal establishes bytes, ordering, and
//! schema semantics ONLY; it does not prove execution correctness,
//! coverage, or physical validity (this module's no-claim).

use std::collections::BTreeMap;
use std::io::Write;
use std::path::{Path, PathBuf};

use fs_blake3::ContentHash;

use crate::vv::{ApplicabilityPolicy, ArtifactKind, ArtifactRef};

/// Versioned domain for the input bundle's semantic identity.
pub const PREDICTION_INPUT_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-evidence.prediction-execution-input.v1";
/// Canonical transport magic (input bundle).
const INPUT_MAGIC: &[u8; 4] = b"FSPI";
/// Canonical transport version byte.
const TRANSPORT_VERSION: u8 = 1;
/// Maximum accepted canonical transport size.
pub const MAX_PREDICTION_BUNDLE_BYTES: usize = 1024 * 1024;
/// Maximum entries for one bounded collection in a bundle.
pub const MAX_BUNDLE_ITEMS: usize = 1_024;
/// Maximum UTF-8 bytes for one bounded string field.
pub const MAX_BUNDLE_TEXT_BYTES: usize = 4 * 1024;

/// Typed refusal for construction, transport, and sealing.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PredictionBundleError {
    /// Stable machine slug (`prediction-input-<field>` style).
    pub rule: &'static str,
    /// Field path that refused.
    pub field: String,
    /// Human diagnosis; never part of semantic identity.
    pub detail: String,
}

impl core::fmt::Display for PredictionBundleError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "{}: {}: {}", self.rule, self.field, self.detail)
    }
}

impl std::error::Error for PredictionBundleError {}

fn refuse(
    rule: &'static str,
    field: impl Into<String>,
    detail: impl Into<String>,
) -> PredictionBundleError {
    PredictionBundleError {
        rule,
        field: field.into(),
        detail: detail.into(),
    }
}

fn checked_text(field: &str, value: &str) -> Result<(), PredictionBundleError> {
    if value.is_empty() || value.len() > MAX_BUNDLE_TEXT_BYTES {
        return Err(refuse(
            "prediction-input-text-bounds",
            field,
            format!("text must be 1..={MAX_BUNDLE_TEXT_BYTES} bytes"),
        ));
    }
    if value.chars().any(char::is_control) {
        return Err(refuse(
            "prediction-input-text-control",
            field,
            "control characters are not admissible",
        ));
    }
    Ok(())
}

/// One named deterministic random stream the executor must honour.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RandomStreamDesign {
    /// Stable stream name (e.g. `"sample-draw"`).
    pub name: String,
    /// Domain-separation string mixed into the stream's seed derivation.
    pub seed_domain: String,
    /// Root seed for this stream.
    pub seed: u64,
    /// Declared number of independent substreams.
    pub substreams: u32,
}

/// Which model rungs the executor MAY use, and what out-of-domain requests do.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ModelRungPolicy {
    /// Admitted rung identities, canonically sorted and deduplicated.
    pub allowed_rungs: Vec<String>,
    /// Reused V&V applicability policy: refuse or degrade out-of-domain.
    pub applicability: ApplicabilityPolicy,
}

/// Who may read the sealed input before the blind release.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AccessPolicy {
    /// Executor processes only; scorers see it after release.
    ExecutorOnly,
    /// Anyone may read the input; only outputs are blind-partitioned.
    Open,
}

impl AccessPolicy {
    const fn tag(self) -> u8 {
        match self {
            Self::ExecutorOnly => 0,
            Self::Open => 1,
        }
    }

    fn from_tag(tag: u8) -> Result<Self, PredictionBundleError> {
        match tag {
            0 => Ok(Self::ExecutorOnly),
            1 => Ok(Self::Open),
            other => Err(refuse(
                "prediction-input-access-policy",
                "access_policy",
                format!("unknown access-policy tag {other}"),
            )),
        }
    }
}

/// The complete, target-free execution input for one blind prediction run.
///
/// Constructed only through [`PredictionExecutionInput::try_new`], so the
/// canonical ordering and bounds invariants hold by construction. There is
/// deliberately NO field for targets, outcomes, or computed predictions,
/// and the canonical decoder refuses any bytes beyond the declared schema.
#[derive(Debug, Clone, PartialEq)]
pub struct PredictionExecutionInput {
    source_identity: BTreeMap<String, String>,
    context_of_use: ArtifactRef,
    validation_plan: ArtifactRef,
    calibration_split: ArtifactRef,
    scenarios: Vec<ArtifactRef>,
    parameter_distributions: Vec<ArtifactRef>,
    random_streams: Vec<RandomStreamDesign>,
    model_rungs: ModelRungPolicy,
    qoi_identities: Vec<String>,
    evidence_role: String,
    blind_partition: ArtifactRef,
    access_policy: AccessPolicy,
}

impl PredictionExecutionInput {
    /// Validate and canonicalize one execution input.
    ///
    /// # Errors
    /// Typed refusals for empty/oversized collections, malformed text,
    /// wrong reference kinds, duplicate identities, and empty stream or
    /// rung declarations. Refusal happens before any state is retained.
    #[expect(
        clippy::too_many_arguments,
        reason = "one transaction admits the complete input or none of it; a \
                  builder would let partially-validated state escape"
    )]
    #[expect(
        clippy::too_many_lines,
        reason = "one transaction keeps mutually dependent admission checks \
                  from admitting partial state (the ArtifactHeader::try_new \
                  precedent)"
    )]
    pub fn try_new(
        source_identity: Vec<(String, String)>,
        context_of_use: ArtifactRef,
        validation_plan: ArtifactRef,
        calibration_split: ArtifactRef,
        scenarios: Vec<ArtifactRef>,
        parameter_distributions: Vec<ArtifactRef>,
        random_streams: Vec<RandomStreamDesign>,
        model_rungs: ModelRungPolicy,
        qoi_identities: Vec<String>,
        evidence_role: String,
        blind_partition: ArtifactRef,
        access_policy: AccessPolicy,
    ) -> Result<Self, PredictionBundleError> {
        if source_identity.is_empty() || source_identity.len() > MAX_BUNDLE_ITEMS {
            return Err(refuse(
                "prediction-input-source-identity",
                "source_identity",
                format!("must declare 1..={MAX_BUNDLE_ITEMS} identity entries"),
            ));
        }
        let mut identity_map = BTreeMap::new();
        for (key, value) in source_identity {
            checked_text("source_identity.key", &key)?;
            checked_text("source_identity.value", &value)?;
            if identity_map.insert(key.clone(), value).is_some() {
                return Err(refuse(
                    "prediction-input-source-identity",
                    "source_identity",
                    format!("duplicate identity key {key:?}"),
                ));
            }
        }
        expect_kind(
            "context_of_use",
            &context_of_use,
            ArtifactKind::ContextOfUse,
        )?;
        expect_kind(
            "validation_plan",
            &validation_plan,
            ArtifactKind::ValidationPlan,
        )?;
        expect_kind(
            "calibration_split",
            &calibration_split,
            ArtifactKind::CalibrationSplit,
        )?;
        expect_kind(
            "blind_partition",
            &blind_partition,
            ArtifactKind::ExperimentArtifact,
        )?;
        bounded_refs("scenarios", &scenarios, 1)?;
        bounded_refs("parameter_distributions", &parameter_distributions, 0)?;
        if random_streams.is_empty() || random_streams.len() > MAX_BUNDLE_ITEMS {
            return Err(refuse(
                "prediction-input-random-streams",
                "random_streams",
                format!("must declare 1..={MAX_BUNDLE_ITEMS} streams"),
            ));
        }
        let mut stream_names = Vec::new();
        for stream in &random_streams {
            checked_text("random_streams.name", &stream.name)?;
            checked_text("random_streams.seed_domain", &stream.seed_domain)?;
            if stream.substreams == 0 {
                return Err(refuse(
                    "prediction-input-random-streams",
                    "random_streams.substreams",
                    "a stream must declare at least one substream",
                ));
            }
            stream_names.push(stream.name.clone());
        }
        stream_names.sort_unstable();
        if stream_names.windows(2).any(|pair| pair[0] == pair[1]) {
            return Err(refuse(
                "prediction-input-random-streams",
                "random_streams",
                "duplicate stream names are not admissible",
            ));
        }
        let mut random_streams = random_streams;
        random_streams.sort_by(|a, b| a.name.cmp(&b.name));
        if model_rungs.allowed_rungs.is_empty()
            || model_rungs.allowed_rungs.len() > MAX_BUNDLE_ITEMS
        {
            return Err(refuse(
                "prediction-input-model-rungs",
                "model_rungs.allowed_rungs",
                format!("must admit 1..={MAX_BUNDLE_ITEMS} rungs"),
            ));
        }
        for rung in &model_rungs.allowed_rungs {
            checked_text("model_rungs.allowed_rungs", rung)?;
        }
        let mut model_rungs = model_rungs;
        model_rungs.allowed_rungs.sort_unstable();
        model_rungs.allowed_rungs.dedup();
        if qoi_identities.is_empty() || qoi_identities.len() > MAX_BUNDLE_ITEMS {
            return Err(refuse(
                "prediction-input-qoi",
                "qoi_identities",
                format!("must declare 1..={MAX_BUNDLE_ITEMS} QoI identities"),
            ));
        }
        for qoi in &qoi_identities {
            checked_text("qoi_identities", qoi)?;
        }
        let mut qoi_identities = qoi_identities;
        qoi_identities.sort_unstable();
        if qoi_identities.windows(2).any(|pair| pair[0] == pair[1]) {
            return Err(refuse(
                "prediction-input-qoi",
                "qoi_identities",
                "duplicate QoI identities are not admissible",
            ));
        }
        checked_text("evidence_role", &evidence_role)?;
        Ok(Self {
            source_identity: identity_map,
            context_of_use,
            validation_plan,
            calibration_split,
            scenarios,
            parameter_distributions,
            random_streams,
            model_rungs,
            qoi_identities,
            evidence_role,
            blind_partition,
            access_policy,
        })
    }

    /// Canonical transport bytes: versioned, bounded, deterministic.
    ///
    /// # Errors
    /// Refuses when the encoding would exceed [`MAX_PREDICTION_BUNDLE_BYTES`].
    pub fn canonical_bytes(&self) -> Result<Vec<u8>, PredictionBundleError> {
        let mut bytes = Vec::new();
        bytes.extend_from_slice(INPUT_MAGIC);
        bytes.push(TRANSPORT_VERSION);
        push_map(&mut bytes, &self.source_identity);
        push_ref(&mut bytes, &self.context_of_use);
        push_ref(&mut bytes, &self.validation_plan);
        push_ref(&mut bytes, &self.calibration_split);
        push_ref_list(&mut bytes, &self.scenarios);
        push_ref_list(&mut bytes, &self.parameter_distributions);
        push_u32(
            &mut bytes,
            cast_len("random_streams", self.random_streams.len())?,
        );
        for stream in &self.random_streams {
            push_text(&mut bytes, &stream.name);
            push_text(&mut bytes, &stream.seed_domain);
            bytes.extend_from_slice(&stream.seed.to_le_bytes());
            bytes.extend_from_slice(&stream.substreams.to_le_bytes());
        }
        push_u32(
            &mut bytes,
            cast_len("model_rungs", self.model_rungs.allowed_rungs.len())?,
        );
        for rung in &self.model_rungs.allowed_rungs {
            push_text(&mut bytes, rung);
        }
        bytes.push(applicability_tag(self.model_rungs.applicability));
        push_u32(
            &mut bytes,
            cast_len("qoi_identities", self.qoi_identities.len())?,
        );
        for qoi in &self.qoi_identities {
            push_text(&mut bytes, qoi);
        }
        push_text(&mut bytes, &self.evidence_role);
        push_ref(&mut bytes, &self.blind_partition);
        bytes.push(self.access_policy.tag());
        if bytes.len() > MAX_PREDICTION_BUNDLE_BYTES {
            return Err(refuse(
                "prediction-input-transport-bounds",
                "canonical_bytes",
                format!("encoding exceeds {MAX_PREDICTION_BUNDLE_BYTES} bytes"),
            ));
        }
        Ok(bytes)
    }

    /// Decode canonical transport bytes, refusing trailing or foreign data.
    ///
    /// # Errors
    /// Typed refusals for magic/version mismatch, truncation, bound
    /// violations, unknown tags, and — critically for the no-target
    /// guarantee — ANY bytes beyond the declared schema.
    pub fn from_canonical_bytes(bytes: &[u8]) -> Result<Self, PredictionBundleError> {
        if bytes.len() > MAX_PREDICTION_BUNDLE_BYTES {
            return Err(refuse(
                "prediction-input-transport-bounds",
                "from_canonical_bytes",
                "transport exceeds the declared maximum",
            ));
        }
        let mut reader = Reader { bytes, offset: 0 };
        let magic = reader.take(4)?;
        if magic != INPUT_MAGIC {
            return Err(refuse(
                "prediction-input-transport-magic",
                "magic",
                "not a prediction-execution-input transport",
            ));
        }
        if reader.take(1)?[0] != TRANSPORT_VERSION {
            return Err(refuse(
                "prediction-input-transport-version",
                "version",
                "unknown transport version",
            ));
        }
        let source_identity = reader.map_entries()?;
        let context_of_use = reader.artifact_ref()?;
        let validation_plan = reader.artifact_ref()?;
        let calibration_split = reader.artifact_ref()?;
        let scenarios = reader.ref_list()?;
        let parameter_distributions = reader.ref_list()?;
        let stream_count = reader.u32_bounded("random_streams")?;
        let mut random_streams = Vec::new();
        for _ in 0..stream_count {
            let name = reader.text()?;
            let seed_domain = reader.text()?;
            let seed = u64::from_le_bytes(reader.take(8)?.try_into().expect("8 bytes"));
            let substreams = u32::from_le_bytes(reader.take(4)?.try_into().expect("4 bytes"));
            random_streams.push(RandomStreamDesign {
                name,
                seed_domain,
                seed,
                substreams,
            });
        }
        let rung_count = reader.u32_bounded("model_rungs")?;
        let mut allowed_rungs = Vec::new();
        for _ in 0..rung_count {
            allowed_rungs.push(reader.text()?);
        }
        let applicability = applicability_from_tag(reader.take(1)?[0])?;
        let qoi_count = reader.u32_bounded("qoi_identities")?;
        let mut qoi_identities = Vec::new();
        for _ in 0..qoi_count {
            qoi_identities.push(reader.text()?);
        }
        let evidence_role = reader.text()?;
        let blind_partition = reader.artifact_ref()?;
        let access_policy = AccessPolicy::from_tag(reader.take(1)?[0])?;
        if reader.offset != bytes.len() {
            return Err(refuse(
                "prediction-input-transport-trailing",
                "from_canonical_bytes",
                format!(
                    "{} undeclared trailing byte(s); a target or outcome smuggled \
                     past the schema would land here, so trailing data refuses",
                    bytes.len() - reader.offset
                ),
            ));
        }
        Self::try_new(
            source_identity,
            context_of_use,
            validation_plan,
            calibration_split,
            scenarios,
            parameter_distributions,
            random_streams,
            ModelRungPolicy {
                allowed_rungs,
                applicability,
            },
            qoi_identities,
            evidence_role,
            blind_partition,
            access_policy,
        )
    }

    /// Semantic identity: the canonical bytes hashed in the versioned
    /// domain. Every semantic edit mints a new identity.
    ///
    /// # Errors
    /// Propagates canonical-encoding refusals.
    pub fn identity(&self) -> Result<ContentHash, PredictionBundleError> {
        Ok(fs_blake3::hash_domain(
            PREDICTION_INPUT_IDENTITY_DOMAIN,
            &self.canonical_bytes()?,
        ))
    }

    /// Referenced blind partition (read access for downstream joins).
    #[must_use]
    pub const fn blind_partition(&self) -> &ArtifactRef {
        &self.blind_partition
    }

    /// Declared access policy.
    #[must_use]
    pub const fn access_policy(&self) -> AccessPolicy {
        self.access_policy
    }
}

fn expect_kind(
    field: &str,
    reference: &ArtifactRef,
    expected: ArtifactKind,
) -> Result<(), PredictionBundleError> {
    if reference.kind() == expected {
        Ok(())
    } else {
        Err(refuse(
            "prediction-input-reference-kind",
            field,
            format!(
                "expected {expected:?} reference, found {:?}",
                reference.kind()
            ),
        ))
    }
}

fn bounded_refs(
    field: &str,
    references: &[ArtifactRef],
    minimum: usize,
) -> Result<(), PredictionBundleError> {
    if references.len() < minimum || references.len() > MAX_BUNDLE_ITEMS {
        return Err(refuse(
            "prediction-input-reference-bounds",
            field,
            format!("must declare {minimum}..={MAX_BUNDLE_ITEMS} references"),
        ));
    }
    let mut hashes: Vec<[u8; 32]> = references
        .iter()
        .map(|reference| *reference.hash().as_bytes())
        .collect();
    hashes.sort_unstable();
    if hashes.windows(2).any(|pair| pair[0] == pair[1]) {
        return Err(refuse(
            "prediction-input-reference-bounds",
            field,
            "duplicate referenced hashes are not admissible",
        ));
    }
    Ok(())
}

fn applicability_tag(policy: ApplicabilityPolicy) -> u8 {
    match policy {
        ApplicabilityPolicy::Demote => 1,
        ApplicabilityPolicy::Refuse => 0,
    }
}

fn applicability_from_tag(tag: u8) -> Result<ApplicabilityPolicy, PredictionBundleError> {
    match tag {
        0 => Ok(ApplicabilityPolicy::Refuse),
        1 => Ok(ApplicabilityPolicy::Demote),
        other => Err(refuse(
            "prediction-input-applicability",
            "model_rungs.applicability",
            format!("unknown applicability tag {other}"),
        )),
    }
}

fn cast_len(field: &str, len: usize) -> Result<u32, PredictionBundleError> {
    u32::try_from(len).map_err(|_| {
        refuse(
            "prediction-input-transport-bounds",
            field,
            "collection length exceeds transport width",
        )
    })
}

fn push_u32(bytes: &mut Vec<u8>, value: u32) {
    bytes.extend_from_slice(&value.to_le_bytes());
}

fn push_text(bytes: &mut Vec<u8>, value: &str) {
    push_u32(
        bytes,
        u32::try_from(value.len()).expect("bounded by admission"),
    );
    bytes.extend_from_slice(value.as_bytes());
}

fn push_map(bytes: &mut Vec<u8>, map: &BTreeMap<String, String>) {
    push_u32(
        bytes,
        u32::try_from(map.len()).expect("bounded by admission"),
    );
    for (key, value) in map {
        push_text(bytes, key);
        push_text(bytes, value);
    }
}

fn push_ref(bytes: &mut Vec<u8>, reference: &ArtifactRef) {
    bytes.push(kind_tag(reference.kind()));
    push_text(bytes, reference.id().as_str());
    bytes.extend_from_slice(reference.hash().as_bytes());
}

fn push_ref_list(bytes: &mut Vec<u8>, references: &[ArtifactRef]) {
    push_u32(
        bytes,
        u32::try_from(references.len()).expect("bounded by admission"),
    );
    for reference in references {
        push_ref(bytes, reference);
    }
}

const fn kind_tag(kind: ArtifactKind) -> u8 {
    match kind {
        ArtifactKind::ContextOfUse => 0,
        ArtifactKind::ValidationPlan => 1,
        ArtifactKind::ExperimentArtifact => 2,
        ArtifactKind::CalibrationSplit => 3,
        ArtifactKind::SolutionVerificationReceipt => 4,
        ArtifactKind::PredictionAssessment => 5,
        ArtifactKind::AssumptionsLedger => 6,
    }
}

fn kind_from_tag(tag: u8) -> Result<ArtifactKind, PredictionBundleError> {
    Ok(match tag {
        0 => ArtifactKind::ContextOfUse,
        1 => ArtifactKind::ValidationPlan,
        2 => ArtifactKind::ExperimentArtifact,
        3 => ArtifactKind::CalibrationSplit,
        4 => ArtifactKind::SolutionVerificationReceipt,
        5 => ArtifactKind::PredictionAssessment,
        6 => ArtifactKind::AssumptionsLedger,
        other => {
            return Err(refuse(
                "prediction-input-reference-kind",
                "artifact_ref.kind",
                format!("unknown artifact-kind tag {other}"),
            ));
        }
    })
}

struct Reader<'a> {
    bytes: &'a [u8],
    offset: usize,
}

impl<'a> Reader<'a> {
    fn take(&mut self, count: usize) -> Result<&'a [u8], PredictionBundleError> {
        let end = self
            .offset
            .checked_add(count)
            .filter(|end| *end <= self.bytes.len());
        let Some(end) = end else {
            return Err(refuse(
                "prediction-input-transport-truncated",
                "transport",
                format!("truncated at offset {}", self.offset),
            ));
        };
        let slice = &self.bytes[self.offset..end];
        self.offset = end;
        Ok(slice)
    }

    fn u32_bounded(&mut self, field: &str) -> Result<u32, PredictionBundleError> {
        let value = u32::from_le_bytes(self.take(4)?.try_into().expect("4 bytes"));
        if value as usize > MAX_BUNDLE_ITEMS {
            return Err(refuse(
                "prediction-input-transport-bounds",
                field,
                format!("declared {value} entries, cap is {MAX_BUNDLE_ITEMS}"),
            ));
        }
        Ok(value)
    }

    fn text(&mut self) -> Result<String, PredictionBundleError> {
        let len = u32::from_le_bytes(self.take(4)?.try_into().expect("4 bytes")) as usize;
        if len > MAX_BUNDLE_TEXT_BYTES {
            return Err(refuse(
                "prediction-input-text-bounds",
                "transport",
                "declared text exceeds the bound",
            ));
        }
        String::from_utf8(self.take(len)?.to_vec()).map_err(|_| {
            refuse(
                "prediction-input-text-bounds",
                "transport",
                "text is not valid UTF-8",
            )
        })
    }

    fn map_entries(&mut self) -> Result<Vec<(String, String)>, PredictionBundleError> {
        let count = self.u32_bounded("source_identity")?;
        let mut entries = Vec::new();
        for _ in 0..count {
            let key = self.text()?;
            let value = self.text()?;
            entries.push((key, value));
        }
        Ok(entries)
    }

    fn artifact_ref(&mut self) -> Result<ArtifactRef, PredictionBundleError> {
        let kind = kind_from_tag(self.take(1)?[0])?;
        let id_text = self.text()?;
        let id = crate::vv::ArtifactId::try_new(&id_text).map_err(|error| {
            refuse(
                "prediction-input-reference-kind",
                "artifact_ref.id",
                format!("{error}"),
            )
        })?;
        let hash_bytes: [u8; 32] = self.take(32)?.try_into().expect("32 bytes");
        Ok(ArtifactRef::new(
            kind,
            id,
            ContentHash::from_slice(&hash_bytes).expect("32 bytes"),
        ))
    }

    fn ref_list(&mut self) -> Result<Vec<ArtifactRef>, PredictionBundleError> {
        let count = self.u32_bounded("references")?;
        let mut references = Vec::new();
        for _ in 0..count {
            references.push(self.artifact_ref()?);
        }
        Ok(references)
    }
}

/// Atomically seal canonical bytes at `path`.
///
/// Bytes land in `<path>.partial.<pid>` first and are atomically renamed,
/// so a crash mid-write leaves NO readable bundle at the sealed path — an
/// interrupted publication is unscoreable by construction. Returns the
/// sealed identity.
///
/// # Errors
/// Encoding refusals, an already-sealed path (seals are immutable; a new
/// semantic edit mints a new identity and belongs at a new path), and I/O
/// failures with the partial path cleaned up.
pub fn seal_prediction_input(
    input: &PredictionExecutionInput,
    path: &Path,
) -> Result<ContentHash, PredictionBundleError> {
    let bytes = input.canonical_bytes()?;
    let identity = input.identity()?;
    seal_bytes(&bytes, path)?;
    Ok(identity)
}

/// Shared atomic publication: partial write, fsync, rename; an existing
/// sealed path refuses (seals are immutable).
fn seal_bytes(bytes: &[u8], path: &Path) -> Result<(), PredictionBundleError> {
    if path.exists() {
        return Err(refuse(
            "prediction-input-seal-immutable",
            path.to_string_lossy(),
            "a sealed bundle is immutable; a semantic edit mints a new identity",
        ));
    }
    let partial: PathBuf = path.with_extension(format!("partial.{}", std::process::id()));
    let write = (|| -> std::io::Result<()> {
        let mut file = std::fs::File::create(&partial)?;
        file.write_all(bytes)?;
        file.sync_all()?;
        std::fs::rename(&partial, path)
    })();
    if let Err(error) = write {
        let _ = std::fs::remove_file(&partial);
        return Err(refuse(
            "prediction-input-seal-io",
            path.to_string_lossy(),
            format!("sealing failed: {error}"),
        ));
    }
    Ok(())
}

/// Load and independently verify a sealed input from artifact bytes alone.
///
/// # Errors
/// Transport refusals, and an identity mismatch when `expected` is given —
/// stale, truncated, swapped, or bit-flipped bundles all land here.
pub fn load_sealed_input(
    path: &Path,
    expected: Option<ContentHash>,
) -> Result<(PredictionExecutionInput, ContentHash), PredictionBundleError> {
    let bytes = std::fs::read(path).map_err(|error| {
        refuse(
            "prediction-input-seal-io",
            path.to_string_lossy(),
            format!("cannot read sealed bundle: {error}"),
        )
    })?;
    let input = PredictionExecutionInput::from_canonical_bytes(&bytes)?;
    let identity = fs_blake3::hash_domain(PREDICTION_INPUT_IDENTITY_DOMAIN, &bytes);
    if let Some(expected) = expected
        && expected != identity
    {
        return Err(refuse(
            "prediction-input-seal-identity",
            path.to_string_lossy(),
            "sealed bytes do not match the expected identity (stale, swapped, or tampered)",
        ));
    }
    Ok((input, identity))
}

/// Versioned domain for the output bundle's semantic identity.
pub const PREDICTION_OUTPUT_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-evidence.prediction-output-bundle.v1";
/// Canonical transport magic (output bundle).
const OUTPUT_MAGIC: &[u8; 4] = b"FSPO";

/// Family of one produced output artifact.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputFamily {
    /// Time-series trajectory artifact.
    Trajectory,
    /// Detected-event artifact.
    Event,
    /// Energy/consistency-track artifact.
    Energy,
    /// Registered aggregate distribution or effect.
    Aggregate,
}

impl OutputFamily {
    const fn tag(self) -> u8 {
        match self {
            Self::Trajectory => 0,
            Self::Event => 1,
            Self::Energy => 2,
            Self::Aggregate => 3,
        }
    }

    fn from_tag(tag: u8) -> Result<Self, PredictionBundleError> {
        Ok(match tag {
            0 => Self::Trajectory,
            1 => Self::Event,
            2 => Self::Energy,
            3 => Self::Aggregate,
            other => {
                return Err(refuse(
                    "prediction-output-artifact-family",
                    "artifacts.family",
                    format!("unknown output-family tag {other}"),
                ));
            }
        })
    }
}

/// Content-addressed reference to one produced output artifact.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OutputArtifactRef {
    /// Artifact family.
    pub family: OutputFamily,
    /// Stable producer-scoped identity.
    pub id: String,
    /// Digest of the artifact's canonical bytes.
    pub hash: ContentHash,
}

/// Exact sample accounting: the denominators a scorer divides by.
///
/// The partition is total by construction: `succeeded + refused + failed`
/// must equal `requested`, so a dropped sample cannot vanish silently.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SampleAccounting {
    /// Samples the input design requested.
    pub requested: u64,
    /// Samples that produced admitted artifacts.
    pub succeeded: u64,
    /// Samples refused by declared policy (e.g. applicability).
    pub refused: u64,
    /// Samples that failed outside declared policy.
    pub failed: u64,
}

/// The complete claimed output of one blind prediction run.
///
/// Binds EXACTLY ONE sealed input root; a scorer joins output to input by
/// identity and [`PredictionOutputBundle::verify_against_input`] refuses a
/// root mismatch, so input-before-output ordering is checkable, not
/// conventional.
#[derive(Debug, Clone, PartialEq)]
pub struct PredictionOutputBundle {
    input_root: ContentHash,
    accounting: SampleAccounting,
    artifacts: Vec<OutputArtifactRef>,
    registered_aggregates: Vec<OutputArtifactRef>,
    work_units: u128,
    budget_exceeded: bool,
    discrepancy: Option<OutputArtifactRef>,
    numerical_evidence: ArtifactRef,
    checker_instructions: Vec<String>,
}

impl PredictionOutputBundle {
    /// Validate and canonicalize one output bundle.
    ///
    /// # Errors
    /// Typed refusals for a non-total sample partition, zero requested
    /// samples, unbounded collections, duplicate artifact hashes, a
    /// non-verification numerical-evidence reference, and malformed text.
    #[expect(
        clippy::too_many_arguments,
        reason = "one transaction admits the complete output or none of it                   (the input-side precedent)"
    )]
    pub fn try_new(
        input_root: ContentHash,
        accounting: SampleAccounting,
        artifacts: Vec<OutputArtifactRef>,
        registered_aggregates: Vec<OutputArtifactRef>,
        work_units: u128,
        budget_exceeded: bool,
        discrepancy: Option<OutputArtifactRef>,
        numerical_evidence: ArtifactRef,
        checker_instructions: Vec<String>,
    ) -> Result<Self, PredictionBundleError> {
        if accounting.requested == 0 {
            return Err(refuse(
                "prediction-output-accounting",
                "accounting.requested",
                "a run that requested zero samples has nothing to claim",
            ));
        }
        let partition = accounting
            .succeeded
            .checked_add(accounting.refused)
            .and_then(|sum| sum.checked_add(accounting.failed));
        if partition != Some(accounting.requested) {
            return Err(refuse(
                "prediction-output-accounting",
                "accounting",
                format!(
                    "succeeded {} + refused {} + failed {} must equal requested {}                      exactly; a sample cannot vanish from the denominator",
                    accounting.succeeded,
                    accounting.refused,
                    accounting.failed,
                    accounting.requested
                ),
            ));
        }
        if accounting.succeeded > 0 && artifacts.is_empty() {
            return Err(refuse(
                "prediction-output-artifacts",
                "artifacts",
                "succeeded samples without any produced artifact are unscoreable",
            ));
        }
        bounded_output_refs("artifacts", &artifacts)?;
        bounded_output_refs("registered_aggregates", &registered_aggregates)?;
        if let Some(reference) = &discrepancy {
            checked_text("discrepancy.id", &reference.id)?;
        }
        expect_kind(
            "numerical_evidence",
            &numerical_evidence,
            ArtifactKind::SolutionVerificationReceipt,
        )?;
        if checker_instructions.is_empty() || checker_instructions.len() > MAX_BUNDLE_ITEMS {
            return Err(refuse(
                "prediction-output-checker",
                "checker_instructions",
                format!("must declare 1..={MAX_BUNDLE_ITEMS} instructions"),
            ));
        }
        for instruction in &checker_instructions {
            checked_text("checker_instructions", instruction)?;
        }
        Ok(Self {
            input_root,
            accounting,
            artifacts,
            registered_aggregates,
            work_units,
            budget_exceeded,
            discrepancy,
            numerical_evidence,
            checker_instructions,
        })
    }

    /// Canonical transport bytes (versioned, bounded, deterministic).
    ///
    /// # Errors
    /// Refuses when the encoding would exceed [`MAX_PREDICTION_BUNDLE_BYTES`].
    pub fn canonical_bytes(&self) -> Result<Vec<u8>, PredictionBundleError> {
        let mut bytes = Vec::new();
        bytes.extend_from_slice(OUTPUT_MAGIC);
        bytes.push(TRANSPORT_VERSION);
        bytes.extend_from_slice(self.input_root.as_bytes());
        bytes.extend_from_slice(&self.accounting.requested.to_le_bytes());
        bytes.extend_from_slice(&self.accounting.succeeded.to_le_bytes());
        bytes.extend_from_slice(&self.accounting.refused.to_le_bytes());
        bytes.extend_from_slice(&self.accounting.failed.to_le_bytes());
        push_output_refs(&mut bytes, &self.artifacts)?;
        push_output_refs(&mut bytes, &self.registered_aggregates)?;
        bytes.extend_from_slice(&self.work_units.to_le_bytes());
        bytes.push(u8::from(self.budget_exceeded));
        match &self.discrepancy {
            None => bytes.push(0),
            Some(reference) => {
                bytes.push(1);
                push_output_ref(&mut bytes, reference);
            }
        }
        push_ref(&mut bytes, &self.numerical_evidence);
        push_u32(
            &mut bytes,
            cast_len("checker_instructions", self.checker_instructions.len())?,
        );
        for instruction in &self.checker_instructions {
            push_text(&mut bytes, instruction);
        }
        if bytes.len() > MAX_PREDICTION_BUNDLE_BYTES {
            return Err(refuse(
                "prediction-output-transport-bounds",
                "canonical_bytes",
                format!("encoding exceeds {MAX_PREDICTION_BUNDLE_BYTES} bytes"),
            ));
        }
        Ok(bytes)
    }

    /// Decode canonical transport bytes, refusing trailing or foreign data.
    ///
    /// # Errors
    /// Typed refusals mirroring the input decoder, with the same
    /// trailing-byte rule closing the outcome-smuggling route.
    pub fn from_canonical_bytes(bytes: &[u8]) -> Result<Self, PredictionBundleError> {
        if bytes.len() > MAX_PREDICTION_BUNDLE_BYTES {
            return Err(refuse(
                "prediction-output-transport-bounds",
                "from_canonical_bytes",
                "transport exceeds the declared maximum",
            ));
        }
        let mut reader = Reader { bytes, offset: 0 };
        if reader.take(4)? != OUTPUT_MAGIC {
            return Err(refuse(
                "prediction-output-transport-magic",
                "magic",
                "not a prediction-output-bundle transport",
            ));
        }
        if reader.take(1)?[0] != TRANSPORT_VERSION {
            return Err(refuse(
                "prediction-output-transport-version",
                "version",
                "unknown transport version",
            ));
        }
        let root_bytes: [u8; 32] = reader.take(32)?.try_into().expect("32 bytes");
        let input_root = ContentHash::from_slice(&root_bytes).expect("32 bytes");
        let requested = u64::from_le_bytes(reader.take(8)?.try_into().expect("8 bytes"));
        let succeeded = u64::from_le_bytes(reader.take(8)?.try_into().expect("8 bytes"));
        let refused = u64::from_le_bytes(reader.take(8)?.try_into().expect("8 bytes"));
        let failed = u64::from_le_bytes(reader.take(8)?.try_into().expect("8 bytes"));
        let artifacts = reader.output_refs()?;
        let registered_aggregates = reader.output_refs()?;
        let work_units = u128::from_le_bytes(reader.take(16)?.try_into().expect("16 bytes"));
        let budget_exceeded = match reader.take(1)?[0] {
            0 => false,
            1 => true,
            other => {
                return Err(refuse(
                    "prediction-output-transport-bounds",
                    "budget_exceeded",
                    format!("boolean tag must be 0 or 1, found {other}"),
                ));
            }
        };
        let discrepancy = match reader.take(1)?[0] {
            0 => None,
            1 => Some(reader.output_ref()?),
            other => {
                return Err(refuse(
                    "prediction-output-transport-bounds",
                    "discrepancy",
                    format!("option tag must be 0 or 1, found {other}"),
                ));
            }
        };
        let numerical_evidence = reader.artifact_ref()?;
        let instruction_count = reader.u32_bounded("checker_instructions")?;
        let mut checker_instructions = Vec::new();
        for _ in 0..instruction_count {
            checker_instructions.push(reader.text()?);
        }
        if reader.offset != bytes.len() {
            return Err(refuse(
                "prediction-output-transport-trailing",
                "from_canonical_bytes",
                format!(
                    "{} undeclared trailing byte(s)",
                    bytes.len() - reader.offset
                ),
            ));
        }
        Self::try_new(
            input_root,
            SampleAccounting {
                requested,
                succeeded,
                refused,
                failed,
            },
            artifacts,
            registered_aggregates,
            work_units,
            budget_exceeded,
            discrepancy,
            numerical_evidence,
            checker_instructions,
        )
    }

    /// Semantic identity in the versioned output domain.
    ///
    /// # Errors
    /// Propagates canonical-encoding refusals.
    pub fn identity(&self) -> Result<ContentHash, PredictionBundleError> {
        Ok(fs_blake3::hash_domain(
            PREDICTION_OUTPUT_IDENTITY_DOMAIN,
            &self.canonical_bytes()?,
        ))
    }

    /// The single input root this output claims to answer.
    #[must_use]
    pub const fn input_root(&self) -> ContentHash {
        self.input_root
    }

    /// Exact sample accounting.
    #[must_use]
    pub const fn accounting(&self) -> SampleAccounting {
        self.accounting
    }

    /// Enforce input-before-output ordering: this output is scoreable only
    /// against the exact sealed input it was produced from.
    ///
    /// # Errors
    /// Refuses when the bound root differs from `sealed_input_identity`.
    pub fn verify_against_input(
        &self,
        sealed_input_identity: ContentHash,
    ) -> Result<(), PredictionBundleError> {
        if self.input_root == sealed_input_identity {
            Ok(())
        } else {
            Err(refuse(
                "prediction-output-input-root",
                "input_root",
                "output binds a different input root; scoring it against this                  input would break the blind-prediction join",
            ))
        }
    }
}

fn bounded_output_refs(
    field: &str,
    references: &[OutputArtifactRef],
) -> Result<(), PredictionBundleError> {
    if references.len() > MAX_BUNDLE_ITEMS {
        return Err(refuse(
            "prediction-output-artifacts",
            field,
            format!("must declare at most {MAX_BUNDLE_ITEMS} references"),
        ));
    }
    let mut hashes: Vec<[u8; 32]> = references
        .iter()
        .map(|reference| *reference.hash.as_bytes())
        .collect();
    hashes.sort_unstable();
    if hashes.windows(2).any(|pair| pair[0] == pair[1]) {
        return Err(refuse(
            "prediction-output-artifacts",
            field,
            "duplicate referenced hashes are not admissible",
        ));
    }
    for reference in references {
        checked_text(field, &reference.id)?;
    }
    Ok(())
}

fn push_output_ref(bytes: &mut Vec<u8>, reference: &OutputArtifactRef) {
    bytes.push(reference.family.tag());
    push_text(bytes, &reference.id);
    bytes.extend_from_slice(reference.hash.as_bytes());
}

fn push_output_refs(
    bytes: &mut Vec<u8>,
    references: &[OutputArtifactRef],
) -> Result<(), PredictionBundleError> {
    push_u32(bytes, cast_len("output_refs", references.len())?);
    for reference in references {
        push_output_ref(bytes, reference);
    }
    Ok(())
}

impl Reader<'_> {
    fn output_ref(&mut self) -> Result<OutputArtifactRef, PredictionBundleError> {
        let family = OutputFamily::from_tag(self.take(1)?[0])?;
        let id = self.text()?;
        let hash_bytes: [u8; 32] = self.take(32)?.try_into().expect("32 bytes");
        Ok(OutputArtifactRef {
            family,
            id,
            hash: ContentHash::from_slice(&hash_bytes).expect("32 bytes"),
        })
    }

    fn output_refs(&mut self) -> Result<Vec<OutputArtifactRef>, PredictionBundleError> {
        let count = self.u32_bounded("output_refs")?;
        let mut references = Vec::new();
        for _ in 0..count {
            references.push(self.output_ref()?);
        }
        Ok(references)
    }
}

/// Atomically seal an output bundle (same partial-then-rename contract as
/// [`seal_prediction_input`]; an interrupted publication is unscoreable).
///
/// # Errors
/// Encoding refusals, an already-sealed path, and I/O failures.
pub fn seal_prediction_output(
    output: &PredictionOutputBundle,
    path: &Path,
) -> Result<ContentHash, PredictionBundleError> {
    let bytes = output.canonical_bytes()?;
    let identity = output.identity()?;
    seal_bytes(&bytes, path)?;
    Ok(identity)
}

/// Load and independently verify a sealed output from artifact bytes alone.
///
/// # Errors
/// Transport refusals and identity mismatch when `expected` is given.
pub fn load_sealed_output(
    path: &Path,
    expected: Option<ContentHash>,
) -> Result<(PredictionOutputBundle, ContentHash), PredictionBundleError> {
    let bytes = std::fs::read(path).map_err(|error| {
        refuse(
            "prediction-output-seal-io",
            path.to_string_lossy(),
            format!("cannot read sealed bundle: {error}"),
        )
    })?;
    let output = PredictionOutputBundle::from_canonical_bytes(&bytes)?;
    let identity = fs_blake3::hash_domain(PREDICTION_OUTPUT_IDENTITY_DOMAIN, &bytes);
    if let Some(expected) = expected
        && expected != identity
    {
        return Err(refuse(
            "prediction-output-seal-identity",
            path.to_string_lossy(),
            "sealed bytes do not match the expected identity (stale, swapped, or tampered)",
        ));
    }
    Ok((output, identity))
}
