//! Canonical identity envelopes for resumable long-running jobs (bead h61n).
//!
//! This module binds the request metadata needed to distinguish a hybrid
//! machine run from a theorem check and to refuse resume-schema guessing. It
//! deliberately stops before governor admission: constructing a request does
//! not authenticate a session grant, reserve a budget, execute work, decode a
//! checkpoint, or establish theorem or scientific correctness.

use core::fmt;

use fs_blake3::{ContentHash, hash_domain};
use fs_package::{
    ReceiptSchemaCatalog, ReceiptSchemaCatalogError, ReceiptSchemaDescriptor,
    ReceiptTransportProfile,
};

/// Semantic identity version of a canonical long-job request.
pub const LONG_JOB_REQUEST_IDENTITY_VERSION: u32 = 1;
/// Domain-separated identity of a canonical long-job request.
pub const LONG_JOB_REQUEST_IDENTITY_DOMAIN: &str = "org.frankensim.fs-session.long-job-request.v1";
/// Maximum UTF-8 bytes in an exact long-job operator.
pub const MAX_LONG_JOB_OPERATOR_BYTES: usize = 128;
/// Maximum UTF-8 bytes in a canonical model-family name.
pub const MAX_LONG_JOB_MODEL_FAMILY_BYTES: usize = 128;
const LONG_JOB_REQUEST_DIGEST_BYTES: usize = 32;

const FIELD_JOB_KIND: u8 = 1;
const FIELD_OPERATOR: u8 = 2;
const FIELD_CORE_NANOSECONDS: u8 = 3;
const FIELD_MEMORY_BYTES: u8 = 4;
const FIELD_WALL_NANOSECONDS: u8 = 5;
const FIELD_MAX_PARALLEL_CORES: u8 = 6;
const FIELD_CANONICAL_PROGRAM_HASH: u8 = 7;
const FIELD_MODEL_FAMILY: u8 = 8;
const FIELD_MODEL_VERSION: u8 = 9;
const FIELD_STATE_SCHEMA_VERSION: u8 = 10;
const FIELD_MODEL_INSTANCE_HASH: u8 = 11;
const FIELD_CONTRACT_HASH: u8 = 12;
const FIELD_CODE_HASH: u8 = 13;
const FIELD_RECEIPT_CATALOG_HASH: u8 = 14;
const FIELD_RESUME_FAMILY: u8 = 15;
const FIELD_RESUME_WIRE_SCHEMA_VERSION: u8 = 16;
const FIELD_RESUME_DESCRIPTOR_HASH: u8 = 17;
const LONG_JOB_KIND_HYBRID_MACHINE_TAG: u8 = 1;
const LONG_JOB_KIND_THEOREM_CHECK_TAG: u8 = 2;

/// Closed class of long-running jobs declared by this identity schema.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LongJobKind {
    /// A versioned hybrid-machine execution.
    HybridMachine,
    /// A bounded theorem-checking or theorem-search execution.
    TheoremCheck,
}

impl LongJobKind {
    const fn tag(self) -> u8 {
        match self {
            Self::HybridMachine => LONG_JOB_KIND_HYBRID_MACHINE_TAG,
            Self::TheoremCheck => LONG_JOB_KIND_THEOREM_CHECK_TAG,
        }
    }
}

/// Exact integer resource request for one long-running job.
///
/// These values are request metadata, not a reservation or grant. A later
/// governor integration must compare them with authenticated live authority
/// and reserve the resources atomically.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LongJobBudget {
    core_nanoseconds: u64,
    memory_bytes: u64,
    wall_nanoseconds: u64,
    max_parallel_cores: u64,
}

impl LongJobBudget {
    /// Construct a nonzero exact budget request.
    pub fn try_new(
        core_nanoseconds: u64,
        memory_bytes: u64,
        wall_nanoseconds: u64,
        max_parallel_cores: u64,
    ) -> Result<Self, LongJobRequestError> {
        require_nonzero_integer("core-nanoseconds", core_nanoseconds)?;
        require_nonzero_integer("memory-bytes", memory_bytes)?;
        require_nonzero_integer("wall-nanoseconds", wall_nanoseconds)?;
        require_nonzero_integer("max-parallel-cores", max_parallel_cores)?;
        Ok(Self {
            core_nanoseconds,
            memory_bytes,
            wall_nanoseconds,
            max_parallel_cores,
        })
    }

    /// Requested aggregate CPU time in nanoseconds.
    #[must_use]
    pub const fn core_nanoseconds(&self) -> u64 {
        self.core_nanoseconds
    }

    /// Requested resident-memory ceiling in bytes.
    #[must_use]
    pub const fn memory_bytes(&self) -> u64 {
        self.memory_bytes
    }

    /// Requested wall-clock ceiling in nanoseconds.
    #[must_use]
    pub const fn wall_nanoseconds(&self) -> u64 {
        self.wall_nanoseconds
    }

    /// Requested maximum simultaneous cores.
    #[must_use]
    pub const fn max_parallel_cores(&self) -> u64 {
        self.max_parallel_cores
    }
}

/// Caller-declared exact receipt-family row for resumable state bytes.
///
/// Catalog membership proves that this row exists, not that its owner semantics
/// are compatible with a particular job kind or model. That compatibility is
/// deliberately reserved for an owner-specific admission adapter.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DeclaredResumeSchema {
    family_id: String,
    wire_schema_version: u32,
    descriptor_hash: ContentHash,
}

impl DeclaredResumeSchema {
    /// Capture an already validated owner descriptor without weakening it to
    /// a family/version guess.
    #[must_use]
    pub fn from_descriptor(descriptor: &ReceiptSchemaDescriptor) -> Self {
        Self {
            family_id: descriptor.family_id().to_owned(),
            wire_schema_version: descriptor.wire_schema_version(),
            descriptor_hash: descriptor.content_hash(),
        }
    }

    /// Globally qualified receipt-family id.
    #[must_use]
    pub fn family_id(&self) -> &str {
        &self.family_id
    }

    /// Exact owner wire schema.
    #[must_use]
    pub const fn wire_schema_version(&self) -> u32 {
        self.wire_schema_version
    }

    /// Exact descriptor identity.
    #[must_use]
    pub const fn descriptor_hash(&self) -> ContentHash {
        self.descriptor_hash
    }
}

/// Model and resume semantics that must survive cancellation and restart.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ResumableModelIdentity {
    family: String,
    model_version: u32,
    state_schema_version: u32,
    model_instance_hash: ContentHash,
    contract_hash: ContentHash,
    code_hash: ContentHash,
    resume_schema: DeclaredResumeSchema,
}

impl ResumableModelIdentity {
    /// Construct one bounded model identity from an exact resume descriptor.
    #[allow(clippy::too_many_arguments)] // Every identity dimension stays explicit at the boundary.
    pub fn try_new(
        family: impl AsRef<str>,
        model_version: u32,
        state_schema_version: u32,
        model_instance_hash: ContentHash,
        contract_hash: ContentHash,
        code_hash: ContentHash,
        resume_descriptor: &ReceiptSchemaDescriptor,
    ) -> Result<Self, LongJobRequestError> {
        let family = family.as_ref();
        validate_model_family(family)?;
        require_nonzero_integer("model-version", u64::from(model_version))?;
        require_nonzero_integer("state-schema-version", u64::from(state_schema_version))?;
        require_nonzero_hash("model-instance-hash", model_instance_hash)?;
        require_nonzero_hash("contract-hash", contract_hash)?;
        require_nonzero_hash("code-hash", code_hash)?;
        Ok(Self {
            family: family.to_owned(),
            model_version,
            state_schema_version,
            model_instance_hash,
            contract_hash,
            code_hash,
            resume_schema: DeclaredResumeSchema::from_descriptor(resume_descriptor),
        })
    }

    /// Canonical model family.
    #[must_use]
    pub fn family(&self) -> &str {
        &self.family
    }

    /// Exact model semantic version.
    #[must_use]
    pub const fn model_version(&self) -> u32 {
        self.model_version
    }

    /// Exact resumable-state schema version.
    #[must_use]
    pub const fn state_schema_version(&self) -> u32 {
        self.state_schema_version
    }

    /// Identity of the exact model instance and canonical parameters.
    #[must_use]
    pub const fn model_instance_hash(&self) -> ContentHash {
        self.model_instance_hash
    }

    /// Identity of the governing model contract.
    #[must_use]
    pub const fn contract_hash(&self) -> ContentHash {
        self.contract_hash
    }

    /// Identity of the executing code semantics.
    #[must_use]
    pub const fn code_hash(&self) -> ContentHash {
        self.code_hash
    }

    /// Exact descriptor expected to govern resumable state bytes.
    #[must_use]
    pub const fn resume_schema(&self) -> &DeclaredResumeSchema {
        &self.resume_schema
    }
}

/// Canonical identity envelope for one resumable long-running request.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LongJobRequest {
    kind: LongJobKind,
    operator: String,
    budget: LongJobBudget,
    canonical_program_hash: ContentHash,
    model: ResumableModelIdentity,
    receipt_catalog_hash: ContentHash,
}

impl LongJobRequest {
    /// Construct an exact request by verifying canonical catalog bytes against
    /// a caller-supplied expected pin. Lookup is exact; no schema fallback is
    /// attempted.
    ///
    /// This validates metadata only. It does not authenticate capability,
    /// reserve budget, inspect a live cancellation gate, or execute work.
    pub fn try_new(
        kind: LongJobKind,
        operator: impl AsRef<str>,
        budget: LongJobBudget,
        canonical_program_hash: ContentHash,
        model: ResumableModelIdentity,
        receipt_catalog_bytes: &[u8],
        expected_catalog_hash: ContentHash,
    ) -> Result<Self, LongJobRequestError> {
        let operator = operator.as_ref();
        validate_exact_operator(operator)?;
        require_nonzero_hash("canonical-program-hash", canonical_program_hash)?;
        require_nonzero_hash("receipt-catalog-hash", expected_catalog_hash)?;

        let resume = model.resume_schema();
        let receipt_catalog =
            ReceiptSchemaCatalog::from_bytes_verified(receipt_catalog_bytes, expected_catalog_hash)
                .map_err(|source| catalog_error(&model, expected_catalog_hash, source))?;
        let descriptor = receipt_catalog
            .require_exact(
                resume.family_id(),
                resume.wire_schema_version(),
                resume.descriptor_hash(),
            )
            .map_err(|source| catalog_error(&model, expected_catalog_hash, source))?;
        match descriptor.transport() {
            ReceiptTransportProfile::DigestOnly => {
                return Err(LongJobRequestError::ResumeTransportUnavailable {
                    catalog_hash: expected_catalog_hash,
                    model_family: model.family.clone(),
                    state_schema_version: model.state_schema_version,
                    family: resume.family_id().to_owned(),
                    wire_schema_version: resume.wire_schema_version(),
                    descriptor_hash: resume.descriptor_hash(),
                });
            }
            ReceiptTransportProfile::CanonicalBytes { maximum_bytes }
                if maximum_bytes > budget.memory_bytes =>
            {
                return Err(LongJobRequestError::ResumeTransportExceedsMemoryBudget {
                    catalog_hash: expected_catalog_hash,
                    model_family: model.family.clone(),
                    state_schema_version: model.state_schema_version,
                    family: resume.family_id().to_owned(),
                    wire_schema_version: resume.wire_schema_version(),
                    descriptor_hash: resume.descriptor_hash(),
                    maximum_bytes,
                    requested_memory_bytes: budget.memory_bytes,
                });
            }
            ReceiptTransportProfile::CanonicalBytes { .. } => {}
        }

        Ok(Self {
            kind,
            operator: operator.to_owned(),
            budget,
            canonical_program_hash,
            model,
            receipt_catalog_hash: expected_catalog_hash,
        })
    }

    /// Closed long-job class.
    #[must_use]
    pub const fn kind(&self) -> LongJobKind {
        self.kind
    }

    /// Exact operator requested for later capability admission.
    #[must_use]
    pub fn operator(&self) -> &str {
        &self.operator
    }

    /// Exact integer resource request.
    #[must_use]
    pub const fn budget(&self) -> LongJobBudget {
        self.budget
    }

    /// Caller-claimed hash of the canonical program bytes.
    #[must_use]
    pub const fn canonical_program_hash(&self) -> ContentHash {
        self.canonical_program_hash
    }

    /// Model and resumable-state semantics.
    #[must_use]
    pub const fn model(&self) -> &ResumableModelIdentity {
        &self.model
    }

    /// Caller-supplied expected receipt-catalog pin verified at construction.
    #[must_use]
    pub const fn receipt_catalog_hash(&self) -> ContentHash {
        self.receipt_catalog_hash
    }

    /// Domain-separated identity over every request field.
    #[must_use]
    pub fn content_hash(&self) -> ContentHash {
        long_job_request_hash_with_schema(
            self,
            LONG_JOB_REQUEST_IDENTITY_VERSION,
            LONG_JOB_REQUEST_IDENTITY_DOMAIN,
        )
    }

    fn identity_preimage(&self) -> Vec<u8> {
        let mut encoder = Encoder::new();
        encoder.u8(FIELD_JOB_KIND);
        encoder.u8(self.kind.tag());
        encoder.u8(FIELD_OPERATOR);
        encoder.string(&self.operator);
        encoder.u8(FIELD_CORE_NANOSECONDS);
        encoder.u64(self.budget.core_nanoseconds);
        encoder.u8(FIELD_MEMORY_BYTES);
        encoder.u64(self.budget.memory_bytes);
        encoder.u8(FIELD_WALL_NANOSECONDS);
        encoder.u64(self.budget.wall_nanoseconds);
        encoder.u8(FIELD_MAX_PARALLEL_CORES);
        encoder.u64(self.budget.max_parallel_cores);
        encoder.u8(FIELD_CANONICAL_PROGRAM_HASH);
        encoder.hash(self.canonical_program_hash);
        encoder.u8(FIELD_MODEL_FAMILY);
        encoder.string(&self.model.family);
        encoder.u8(FIELD_MODEL_VERSION);
        encoder.u32(self.model.model_version);
        encoder.u8(FIELD_STATE_SCHEMA_VERSION);
        encoder.u32(self.model.state_schema_version);
        encoder.u8(FIELD_MODEL_INSTANCE_HASH);
        encoder.hash(self.model.model_instance_hash);
        encoder.u8(FIELD_CONTRACT_HASH);
        encoder.hash(self.model.contract_hash);
        encoder.u8(FIELD_CODE_HASH);
        encoder.hash(self.model.code_hash);
        encoder.u8(FIELD_RECEIPT_CATALOG_HASH);
        encoder.hash(self.receipt_catalog_hash);
        encoder.u8(FIELD_RESUME_FAMILY);
        encoder.string(&self.model.resume_schema.family_id);
        encoder.u8(FIELD_RESUME_WIRE_SCHEMA_VERSION);
        encoder.u32(self.model.resume_schema.wire_schema_version);
        encoder.u8(FIELD_RESUME_DESCRIPTOR_HASH);
        encoder.hash(self.model.resume_schema.descriptor_hash);
        encoder.finish()
    }
}

fn long_job_request_hash_with_schema(
    request: &LongJobRequest,
    identity_version: u32,
    domain: &str,
) -> ContentHash {
    let fields = request.identity_preimage();
    let mut preimage = Vec::with_capacity(4 + fields.len());
    preimage.extend_from_slice(&identity_version.to_le_bytes());
    preimage.extend_from_slice(&fields);
    hash_domain(domain, &preimage)
}

fn long_job_request_identity_transport_is_current(identity_version: u32, digest: &[u8]) -> bool {
    identity_version == LONG_JOB_REQUEST_IDENTITY_VERSION
        && digest.len() == LONG_JOB_REQUEST_DIGEST_BYTES
}

/// Typed refusal at the long-job identity boundary.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LongJobRequestError {
    /// A string or semantic value is not canonical.
    InvalidField {
        /// Stable field name.
        field: &'static str,
        /// Stable requirement, excluding caller-controlled input bytes.
        requirement: &'static str,
    },
    /// One exact byte/count bound was exceeded.
    ResourceLimit {
        /// Stable resource name.
        resource: &'static str,
        /// Configured maximum.
        limit: u64,
        /// Exact observed amount.
        observed: u64,
    },
    /// Catalog transport verification or exact descriptor lookup failed.
    ReceiptCatalog {
        /// Caller-supplied expected catalog identity.
        catalog_hash: ContentHash,
        /// Exact model family requesting the row.
        model_family: String,
        /// Exact model state schema requesting the row.
        state_schema_version: u32,
        /// Declared receipt family.
        resume_family: String,
        /// Declared receipt wire schema.
        resume_wire_schema_version: u32,
        /// Declared descriptor identity.
        resume_descriptor_hash: ContentHash,
        /// Bounded catalog-layer refusal.
        source: ReceiptSchemaCatalogError,
    },
    /// The exact descriptor retains only a digest and cannot govern resumable
    /// bytes.
    ResumeTransportUnavailable {
        /// Caller-supplied expected catalog identity.
        catalog_hash: ContentHash,
        /// Exact model family requesting the row.
        model_family: String,
        /// Exact model state schema requesting the row.
        state_schema_version: u32,
        /// Exact receipt family.
        family: String,
        /// Exact wire schema.
        wire_schema_version: u32,
        /// Exact descriptor identity.
        descriptor_hash: ContentHash,
    },
    /// The descriptor's maximum canonical bytes exceed the request's entire
    /// memory ceiling.
    ResumeTransportExceedsMemoryBudget {
        /// Caller-supplied expected catalog identity.
        catalog_hash: ContentHash,
        /// Exact model family requesting the row.
        model_family: String,
        /// Exact model state schema requesting the row.
        state_schema_version: u32,
        /// Exact receipt family.
        family: String,
        /// Exact wire schema.
        wire_schema_version: u32,
        /// Exact descriptor identity.
        descriptor_hash: ContentHash,
        /// Maximum canonical bytes declared by the descriptor.
        maximum_bytes: u64,
        /// Whole-job memory requested by this envelope.
        requested_memory_bytes: u64,
    },
}

impl fmt::Display for LongJobRequestError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidField { field, requirement } => {
                write!(f, "invalid long-job {field}: {requirement}")
            }
            Self::ResourceLimit {
                resource,
                limit,
                observed,
            } => write!(
                f,
                "long-job {resource} limit {limit} exceeded by exact observation {observed}"
            ),
            Self::ReceiptCatalog {
                catalog_hash,
                model_family,
                state_schema_version,
                resume_family,
                resume_wire_schema_version,
                resume_descriptor_hash,
                source,
            } => write!(
                f,
                "receipt catalog {catalog_hash} refused model {model_family} state schema v{state_schema_version} resume row {resume_family} wire v{resume_wire_schema_version} descriptor {resume_descriptor_hash}: {source}"
            ),
            Self::ResumeTransportUnavailable {
                catalog_hash,
                model_family,
                state_schema_version,
                family,
                wire_schema_version,
                descriptor_hash,
            } => write!(
                f,
                "receipt catalog {catalog_hash} row {family} wire v{wire_schema_version} descriptor {descriptor_hash} for model {model_family} state schema v{state_schema_version} is digest-only; resumable jobs require bounded canonical bytes"
            ),
            Self::ResumeTransportExceedsMemoryBudget {
                catalog_hash,
                model_family,
                state_schema_version,
                family,
                wire_schema_version,
                descriptor_hash,
                maximum_bytes,
                requested_memory_bytes,
            } => write!(
                f,
                "receipt catalog {catalog_hash} row {family} wire v{wire_schema_version} descriptor {descriptor_hash} for model {model_family} state schema v{state_schema_version} admits {maximum_bytes} canonical bytes, exceeding the requested {requested_memory_bytes}-byte job memory ceiling"
            ),
        }
    }
}

impl std::error::Error for LongJobRequestError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::ReceiptCatalog { source, .. } => Some(source),
            _ => None,
        }
    }
}

fn catalog_error(
    model: &ResumableModelIdentity,
    catalog_hash: ContentHash,
    source: ReceiptSchemaCatalogError,
) -> LongJobRequestError {
    LongJobRequestError::ReceiptCatalog {
        catalog_hash,
        model_family: model.family.clone(),
        state_schema_version: model.state_schema_version,
        resume_family: model.resume_schema.family_id.clone(),
        resume_wire_schema_version: model.resume_schema.wire_schema_version,
        resume_descriptor_hash: model.resume_schema.descriptor_hash,
        source,
    }
}

fn validate_exact_operator(value: &str) -> Result<(), LongJobRequestError> {
    validate_bounded_name("operator", value, MAX_LONG_JOB_OPERATOR_BYTES)?;
    if value.contains('*') || !fs_ir::admission::valid_operator_pattern(value) {
        return Err(LongJobRequestError::InvalidField {
            field: "operator",
            requirement: "must be an exact canonical operator name with no wildcard",
        });
    }
    Ok(())
}

fn validate_model_family(value: &str) -> Result<(), LongJobRequestError> {
    validate_bounded_name("model-family", value, MAX_LONG_JOB_MODEL_FAMILY_BYTES)?;
    if !value.contains('.')
        || !value.is_ascii()
        || !value.bytes().all(|byte| {
            byte.is_ascii_lowercase() || byte.is_ascii_digit() || matches!(byte, b'.' | b'_' | b'-')
        })
        || !value
            .as_bytes()
            .first()
            .is_some_and(|byte| byte.is_ascii_alphanumeric())
        || !value
            .as_bytes()
            .last()
            .is_some_and(|byte| byte.is_ascii_alphanumeric())
        || value.split(['.', '_', '-']).any(str::is_empty)
        || value.split(['.', '_', '-']).any(is_placeholder_component)
    {
        return Err(LongJobRequestError::InvalidField {
            field: "model-family",
            requirement: "must be a qualified lowercase ASCII model family with no adjacent separators or placeholder components",
        });
    }
    Ok(())
}

fn is_placeholder_component(component: &str) -> bool {
    matches!(
        component,
        "todo" | "tbd" | "placeholder" | "pending" | "unknown" | "none" | "na"
    )
}

fn validate_bounded_name(
    field: &'static str,
    value: &str,
    maximum_bytes: usize,
) -> Result<(), LongJobRequestError> {
    if value.len() > maximum_bytes {
        return Err(LongJobRequestError::ResourceLimit {
            resource: field,
            limit: maximum_bytes as u64,
            observed: value.len() as u64,
        });
    }
    if value.is_empty() {
        return Err(LongJobRequestError::InvalidField {
            field,
            requirement: "must not be empty",
        });
    }
    Ok(())
}

fn require_nonzero_integer(field: &'static str, value: u64) -> Result<(), LongJobRequestError> {
    if value == 0 {
        Err(LongJobRequestError::InvalidField {
            field,
            requirement: "must be nonzero",
        })
    } else {
        Ok(())
    }
}

fn require_nonzero_hash(
    field: &'static str,
    value: ContentHash,
) -> Result<(), LongJobRequestError> {
    if value.as_bytes().iter().all(|byte| *byte == 0) {
        Err(LongJobRequestError::InvalidField {
            field,
            requirement: "must be a nonzero 32-byte content identity",
        })
    } else {
        Ok(())
    }
}

struct Encoder {
    bytes: Vec<u8>,
}

impl Encoder {
    fn new() -> Self {
        Self {
            bytes: Vec::with_capacity(512),
        }
    }

    fn u8(&mut self, value: u8) {
        self.bytes.push(value);
    }

    fn u32(&mut self, value: u32) {
        self.bytes.extend_from_slice(&value.to_le_bytes());
    }

    fn u64(&mut self, value: u64) {
        self.bytes.extend_from_slice(&value.to_le_bytes());
    }

    fn string(&mut self, value: &str) {
        self.u64(value.len() as u64);
        self.bytes.extend_from_slice(value.as_bytes());
    }

    fn hash(&mut self, value: ContentHash) {
        self.bytes.extend_from_slice(value.as_bytes());
    }

    fn finish(self) -> Vec<u8> {
        self.bytes
    }
}

#[allow(dead_code)]
fn classify_long_job_request_identity_fields(request: &LongJobRequest) {
    let LongJobRequest {
        kind,
        operator,
        budget,
        canonical_program_hash,
        model,
        receipt_catalog_hash,
    } = request;
    let kind_variant = match kind {
        LongJobKind::HybridMachine => LONG_JOB_KIND_HYBRID_MACHINE_TAG,
        LongJobKind::TheoremCheck => LONG_JOB_KIND_THEOREM_CHECK_TAG,
    };
    let LongJobBudget {
        core_nanoseconds,
        memory_bytes,
        wall_nanoseconds,
        max_parallel_cores,
    } = budget;
    let ResumableModelIdentity {
        family,
        model_version,
        state_schema_version,
        model_instance_hash,
        contract_hash,
        code_hash,
        resume_schema,
    } = model;
    let DeclaredResumeSchema {
        family_id,
        wire_schema_version,
        descriptor_hash,
    } = resume_schema;
    let _ = (
        kind_variant,
        operator,
        core_nanoseconds,
        memory_bytes,
        wall_nanoseconds,
        max_parallel_cores,
        canonical_program_hash,
        family,
        model_version,
        state_schema_version,
        model_instance_hash,
        contract_hash,
        code_hash,
        receipt_catalog_hash,
        family_id,
        wire_schema_version,
        descriptor_hash,
    );
}

/// Owner-local declaration consumed by `xtask check-identities`.
#[allow(dead_code)]
pub const LONG_JOB_REQUEST_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-session:long-job-request",
    "version_const=LONG_JOB_REQUEST_IDENTITY_VERSION",
    "version=1",
    "domain=org.frankensim.fs-session.long-job-request.v1",
    "domain_const=LONG_JOB_REQUEST_IDENTITY_DOMAIN",
    "encoder=LongJobRequest::content_hash",
    "encoder_helpers=long_job_request_hash_with_schema,LongJobRequest::identity_preimage,LongJobKind::tag,Encoder::new,Encoder::u8,Encoder::u32,Encoder::u64,Encoder::string,Encoder::hash,Encoder::finish",
    "schema_constants=LONG_JOB_REQUEST_IDENTITY_VERSION,LONG_JOB_REQUEST_IDENTITY_DOMAIN,MAX_LONG_JOB_OPERATOR_BYTES,MAX_LONG_JOB_MODEL_FAMILY_BYTES,LONG_JOB_REQUEST_DIGEST_BYTES,FIELD_JOB_KIND,FIELD_OPERATOR,FIELD_CORE_NANOSECONDS,FIELD_MEMORY_BYTES,FIELD_WALL_NANOSECONDS,FIELD_MAX_PARALLEL_CORES,FIELD_CANONICAL_PROGRAM_HASH,FIELD_MODEL_FAMILY,FIELD_MODEL_VERSION,FIELD_STATE_SCHEMA_VERSION,FIELD_MODEL_INSTANCE_HASH,FIELD_CONTRACT_HASH,FIELD_CODE_HASH,FIELD_RECEIPT_CATALOG_HASH,FIELD_RESUME_FAMILY,FIELD_RESUME_WIRE_SCHEMA_VERSION,FIELD_RESUME_DESCRIPTOR_HASH,LONG_JOB_KIND_HYBRID_MACHINE_TAG,LONG_JOB_KIND_THEOREM_CHECK_TAG",
    "schema_functions=LongJobRequest::try_new,LongJobBudget::try_new,ResumableModelIdentity::try_new,DeclaredResumeSchema::from_descriptor,catalog_error,validate_exact_operator,validate_model_family,is_placeholder_component,validate_bounded_name,require_nonzero_integer,require_nonzero_hash,crates/fs-package/src/receipt_catalog.rs#ReceiptSchemaCatalog::from_bytes_verified,crates/fs-package/src/receipt_catalog.rs#ReceiptSchemaCatalog::require_exact,crates/fs-blake3/src/lib.rs#hash_domain,crates/fs-ir/src/admission.rs#valid_operator_pattern",
    "schema_types=crates/fs-package/src/receipt_catalog.rs#ReceiptSchemaDescriptor",
    "schema_dependencies=fs-package:receipt-schema-catalog,fs-package:receipt-schema-descriptor",
    "digest=blake3-256-domain-separated",
    "encoding=typed-binary",
    "sources=LongJobRequest,LongJobKind,LongJobBudget,ResumableModelIdentity,DeclaredResumeSchema",
    "source_fields=LongJobRequest.kind:derived:expanded-into-LongJobKind,LongJobRequest.operator:semantic,LongJobRequest.budget:derived:expanded-into-LongJobBudget,LongJobRequest.canonical_program_hash:semantic,LongJobRequest.model:derived:expanded-into-ResumableModelIdentity,LongJobRequest.receipt_catalog_hash:semantic,LongJobKind.variant:semantic,LongJobBudget.core_nanoseconds:semantic,LongJobBudget.memory_bytes:semantic,LongJobBudget.wall_nanoseconds:semantic,LongJobBudget.max_parallel_cores:semantic,ResumableModelIdentity.family:semantic,ResumableModelIdentity.model_version:semantic,ResumableModelIdentity.state_schema_version:semantic,ResumableModelIdentity.model_instance_hash:semantic,ResumableModelIdentity.contract_hash:semantic,ResumableModelIdentity.code_hash:semantic,ResumableModelIdentity.resume_schema:derived:expanded-into-DeclaredResumeSchema,DeclaredResumeSchema.family_id:semantic,DeclaredResumeSchema.wire_schema_version:semantic,DeclaredResumeSchema.descriptor_hash:semantic",
    "source_bindings=LongJobRequest.operator>operator-byte-count+operator-utf8,LongJobRequest.canonical_program_hash>canonical-program-hash,LongJobRequest.receipt_catalog_hash>receipt-catalog-hash,LongJobKind.variant>job-kind,LongJobBudget.core_nanoseconds>core-nanoseconds,LongJobBudget.memory_bytes>memory-bytes,LongJobBudget.wall_nanoseconds>wall-nanoseconds,LongJobBudget.max_parallel_cores>max-parallel-cores,ResumableModelIdentity.family>model-family-byte-count+model-family-utf8,ResumableModelIdentity.model_version>model-version,ResumableModelIdentity.state_schema_version>state-schema-version,ResumableModelIdentity.model_instance_hash>model-instance-hash,ResumableModelIdentity.contract_hash>contract-hash,ResumableModelIdentity.code_hash>code-hash,DeclaredResumeSchema.family_id>resume-family-byte-count+resume-family-utf8,DeclaredResumeSchema.wire_schema_version>resume-wire-schema-version,DeclaredResumeSchema.descriptor_hash>resume-descriptor-hash",
    "external_semantic_fields=identity-domain,identity-version,canonical-field-order,field-tag-u8,length-count-u64-le,fixed-numeric-little-endian",
    "semantic_fields=identity-domain,identity-version,canonical-field-order,field-tag-u8,length-count-u64-le,fixed-numeric-little-endian,job-kind,operator-byte-count,operator-utf8,core-nanoseconds,memory-bytes,wall-nanoseconds,max-parallel-cores,canonical-program-hash,model-family-byte-count,model-family-utf8,model-version,state-schema-version,model-instance-hash,contract-hash,code-hash,receipt-catalog-hash,resume-family-byte-count,resume-family-utf8,resume-wire-schema-version,resume-descriptor-hash",
    "excluded_fields=none",
    "consumers=LongJobRequest::content_hash,future-Governor::admit_long_job,fs-ledger-long-job-adapter",
    "mutations=identity-domain:crates/fs-session/src/long_job.rs#long_job_request_version_and_domain_move_identity,identity-version:crates/fs-session/src/long_job.rs#long_job_request_version_and_domain_move_identity,canonical-field-order:crates/fs-session/tests/long_job.rs#long_job_request_identity_matches_independent_preimage,field-tag-u8:crates/fs-session/tests/long_job.rs#long_job_request_identity_matches_independent_preimage,length-count-u64-le:crates/fs-session/tests/long_job.rs#long_job_request_identity_matches_independent_preimage,fixed-numeric-little-endian:crates/fs-session/tests/long_job.rs#long_job_request_identity_matches_independent_preimage,job-kind:crates/fs-session/src/long_job.rs#long_job_kind_tags_are_stable_and_move_identity,operator-byte-count:crates/fs-session/tests/long_job.rs#long_job_kind_and_every_field_move_identity,operator-utf8:crates/fs-session/tests/long_job.rs#long_job_kind_and_every_field_move_identity,core-nanoseconds:crates/fs-session/tests/long_job.rs#long_job_kind_and_every_field_move_identity,memory-bytes:crates/fs-session/tests/long_job.rs#long_job_kind_and_every_field_move_identity,wall-nanoseconds:crates/fs-session/tests/long_job.rs#long_job_kind_and_every_field_move_identity,max-parallel-cores:crates/fs-session/tests/long_job.rs#long_job_kind_and_every_field_move_identity,canonical-program-hash:crates/fs-session/tests/long_job.rs#long_job_kind_and_every_field_move_identity,model-family-byte-count:crates/fs-session/tests/long_job.rs#long_job_kind_and_every_field_move_identity,model-family-utf8:crates/fs-session/tests/long_job.rs#long_job_kind_and_every_field_move_identity,model-version:crates/fs-session/tests/long_job.rs#long_job_kind_and_every_field_move_identity,state-schema-version:crates/fs-session/tests/long_job.rs#long_job_kind_and_every_field_move_identity,model-instance-hash:crates/fs-session/tests/long_job.rs#long_job_kind_and_every_field_move_identity,contract-hash:crates/fs-session/tests/long_job.rs#long_job_kind_and_every_field_move_identity,code-hash:crates/fs-session/tests/long_job.rs#long_job_kind_and_every_field_move_identity,receipt-catalog-hash:crates/fs-session/src/long_job.rs#private_resume_and_catalog_fields_move_independently,resume-family-byte-count:crates/fs-session/src/long_job.rs#private_resume_and_catalog_fields_move_independently,resume-family-utf8:crates/fs-session/src/long_job.rs#private_resume_and_catalog_fields_move_independently,resume-wire-schema-version:crates/fs-session/src/long_job.rs#private_resume_and_catalog_fields_move_independently,resume-descriptor-hash:crates/fs-session/src/long_job.rs#private_resume_and_catalog_fields_move_independently",
    "nonsemantic_mutations=none",
    "field_guard=classify_long_job_request_identity_fields",
    "transport_guard=long_job_request_identity_transport_is_current",
    "version_guard=crates/fs-session/src/long_job.rs#long_job_request_version_and_domain_move_identity",
    "coupling_surface=fs-session:long-job-request",
];

#[cfg(test)]
mod tests {
    use super::*;

    fn hash(byte: u8) -> ContentHash {
        ContentHash([byte; 32])
    }

    fn request() -> LongJobRequest {
        let descriptor = ReceiptSchemaDescriptor::try_new(
            "fs-ir:hybrid-machine-checkpoint",
            1,
            1,
            "org.frankensim.fs-ir.hybrid-machine-checkpoint.v1",
            ReceiptTransportProfile::CanonicalBytes {
                maximum_bytes: 1_048_576,
            },
            hash(1),
        )
        .expect("descriptor");
        let catalog = ReceiptSchemaCatalog::try_new(vec![descriptor.clone()]).expect("catalog");
        let model = ResumableModelIdentity::try_new(
            "machine.hybrid",
            1,
            1,
            hash(2),
            hash(3),
            hash(4),
            &descriptor,
        )
        .expect("model");
        LongJobRequest::try_new(
            LongJobKind::HybridMachine,
            "machine.run",
            LongJobBudget::try_new(10, 2_000_000, 30, 2).expect("budget"),
            hash(5),
            model,
            &catalog.to_bytes(),
            catalog.content_hash(),
        )
        .expect("request")
    }

    #[test]
    fn long_job_request_version_and_domain_move_identity() {
        let request = request();
        let current = request.content_hash();
        assert_ne!(
            current,
            long_job_request_hash_with_schema(
                &request,
                LONG_JOB_REQUEST_IDENTITY_VERSION + 1,
                LONG_JOB_REQUEST_IDENTITY_DOMAIN,
            )
        );
        assert_ne!(
            current,
            long_job_request_hash_with_schema(
                &request,
                LONG_JOB_REQUEST_IDENTITY_VERSION,
                "org.frankensim.fs-session.long-job-request.alternate.v1",
            )
        );
        assert!(long_job_request_identity_transport_is_current(
            LONG_JOB_REQUEST_IDENTITY_VERSION,
            current.as_bytes(),
        ));
        assert!(!long_job_request_identity_transport_is_current(
            LONG_JOB_REQUEST_IDENTITY_VERSION + 1,
            current.as_bytes(),
        ));
        assert!(!long_job_request_identity_transport_is_current(
            LONG_JOB_REQUEST_IDENTITY_VERSION,
            &current.as_bytes()[..LONG_JOB_REQUEST_DIGEST_BYTES - 1],
        ));
    }

    #[test]
    fn long_job_kind_tags_are_stable_and_move_identity() {
        assert_eq!(LongJobKind::HybridMachine.tag(), 1);
        assert_eq!(LongJobKind::TheoremCheck.tag(), 2);

        let baseline = request();
        let mut moved = baseline.clone();
        moved.kind = LongJobKind::TheoremCheck;
        assert_ne!(baseline.content_hash(), moved.content_hash());
    }

    #[test]
    fn private_resume_and_catalog_fields_move_independently() {
        let baseline = request();
        let baseline_hash = baseline.content_hash();

        let mut moved = baseline.clone();
        moved.receipt_catalog_hash = hash(10);
        assert_ne!(baseline_hash, moved.content_hash(), "catalog pin");

        let mut moved = baseline.clone();
        moved.model.resume_schema.family_id = "fs-ir:hybrid-machine-checkpoinu".to_owned();
        assert_eq!(
            moved.model.resume_schema.family_id.len(),
            baseline.model.resume_schema.family_id.len(),
        );
        assert_ne!(baseline_hash, moved.content_hash(), "resume family bytes");

        let mut moved = baseline.clone();
        moved.model.resume_schema.family_id.push_str("-next");
        assert_ne!(
            moved.model.resume_schema.family_id.len(),
            baseline.model.resume_schema.family_id.len(),
        );
        assert_ne!(baseline_hash, moved.content_hash(), "resume family length");

        let mut moved = baseline.clone();
        moved.model.resume_schema.wire_schema_version += 1;
        assert_ne!(baseline_hash, moved.content_hash(), "resume wire schema");

        let mut moved = baseline.clone();
        moved.model.resume_schema.descriptor_hash = hash(11);
        assert_ne!(baseline_hash, moved.content_hash(), "resume descriptor");
    }
}
