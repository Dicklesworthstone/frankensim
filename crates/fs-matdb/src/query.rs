//! The query path (bead 5hmy, PR-4 of 5): every material answer is
//! `Evidence<PropertySample>` PLUS a [`PropertyUsageReceipt`] — never a
//! bare number.
//!
//! Discipline, in order:
//! - the query point is validated (finite, named axes);
//! - only claims whose [`fs_evidence::ValidityDomain`] CONTAINS the
//!   point are candidates — evaluation outside validity is a typed
//!   refusal, never a silent extrapolation;
//! - selection among candidates is an EXPLICIT policy; conflicting
//!   claims are never averaged into an invented canonical value —
//!   ambiguity refuses and names the candidates;
//! - the evidence slices map honestly: the datum's band is the STATED
//!   uncertainty (statistical slice); `Unstated` uncertainty maps to an
//!   explicit numerical no-claim, not a manufactured certificate; the
//!   claim's validity and in-domain fact live in the model slice; and
//!   the receipt records what was considered, selected, and decided.

use std::{
    collections::{BTreeMap, BTreeSet},
    fmt,
};

use fs_blake3::{ContentHash, hash_domain};
use fs_evidence::{
    Evidence, ModelEvidence, NumericalCertificate, ProvenanceHash, SensitivitySummary,
    StatisticalCertificate,
};
use fs_qty::Dims;

use crate::{
    ClaimId, ClaimSet, InterpolationPolicy, MatDbError, PropertyClaim, PropertyValue,
    UncertaintyModel,
};

/// Semantic version of the portable property-usage receipt identity.
///
/// Version 2 replaces the historical unframed v1 preimage. V1 did not bind
/// collection counts, so moving one claim id between adjacent collections
/// could preserve the same byte stream. No v1 decoder is provided because
/// those bytes cannot recover the missing boundaries.
pub const PROPERTY_USAGE_RECEIPT_IDENTITY_VERSION: u32 = 2;
/// Exact BLAKE3 domain for the portable property-usage receipt identity.
pub const PROPERTY_USAGE_RECEIPT_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-matdb.property-usage-receipt.v2";
/// Closed wire/schema version retained inside every portable receipt.
pub const PROPERTY_USAGE_RECEIPT_SCHEMA_VERSION: u32 = 2;

const PROPERTY_USAGE_RECEIPT_MAGIC: &[u8; 8] = b"FSMATUR\0";
const FIELD_PROPERTY: u8 = 1;
const FIELD_QUERY_POINT: u8 = 2;
const FIELD_CONSIDERED: u8 = 3;
const FIELD_IN_DOMAIN: u8 = 4;
const FIELD_SELECTED: u8 = 5;
const FIELD_POLICY: u8 = 6;
const FIELD_DECISION: u8 = 7;
const FIELD_OBSERVATION_BACKED: u8 = 8;
const FIELD_EVALUATOR_VERSION: u8 = 9;
const FIELD_SOURCE_HASHES: u8 = 10;

const DECISION_CONSTANT_WITHIN_VALIDITY: u8 = 1;
const DECISION_EXACT_SCALAR: u8 = 2;
const DECISION_EXACT_TABULATED: u8 = 3;
const DECISION_LINEAR_INSIDE: u8 = 4;

/// Maximum canonical bytes accepted for one portable property-usage receipt.
pub const MAX_PROPERTY_USAGE_RECEIPT_BYTES: usize = 1024 * 1024;
/// Maximum UTF-8 bytes in the property name.
pub const MAX_PROPERTY_USAGE_PROPERTY_BYTES: usize = 256;
/// Maximum named coordinates in one query point.
pub const MAX_PROPERTY_USAGE_QUERY_AXES: usize = 256;
/// Maximum UTF-8 bytes in one query-axis name.
pub const MAX_PROPERTY_USAGE_AXIS_BYTES: usize = 128;
/// Maximum claim ids in either the considered or in-domain collection.
pub const MAX_PROPERTY_USAGE_CLAIM_IDS: usize = 8_192;
/// Maximum UTF-8 bytes in the closed selection-policy tag.
pub const MAX_PROPERTY_USAGE_POLICY_BYTES: usize = 64;
/// Maximum selected-claim/source artifact hashes in one receipt.
pub const MAX_PROPERTY_USAGE_SOURCE_HASHES: usize = 8_192;

/// The query evaluator's semantic version (recorded in every receipt;
/// bumped when selection or evaluation semantics change).
pub const MATDB_EVALUATOR_VERSION: u32 = 1;

/// A fail-closed refusal at the portable property-usage receipt boundary.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PropertyUsageReceiptError {
    /// The retained schema is not the only schema this build can interpret.
    UnsupportedSchemaVersion {
        /// Version found in the receipt.
        found: u32,
        /// Version understood by this build.
        supported: u32,
    },
    /// A semantic field cannot be represented by the canonical v2 profile.
    InvalidField {
        /// Stable field/rule name.
        field: &'static str,
        /// Teaching detail.
        detail: String,
    },
    /// A collection or byte field exceeded its public processing budget.
    ResourceLimit {
        /// Stable resource name.
        resource: &'static str,
        /// Configured maximum.
        limit: usize,
        /// Exact observed count or byte length.
        observed: u64,
    },
    /// The closed selection-policy tag is not owned by this evaluator.
    UnknownPolicyTag {
        /// Foreign tag retained by the bytes.
        tag: String,
    },
    /// A decision discriminant is not part of schema v2.
    UnknownDecisionTag {
        /// Foreign discriminant.
        tag: u8,
        /// Byte offset immediately after the discriminant.
        at: usize,
    },
    /// The binary envelope is truncated, malformed, or non-canonical.
    Malformed {
        /// Byte offset at which decoding refused.
        at: usize,
        /// Stable diagnostic detail.
        detail: String,
    },
    /// The retained in-band identity did not reproduce from decoded fields.
    IdentityMismatch {
        /// Identity retained by the transport.
        expected: ContentHash,
        /// Identity recomputed from the decoded receipt.
        actual: ContentHash,
    },
    /// An external caller-pinned identity did not match the admitted receipt.
    ExternalIdentityMismatch {
        /// Identity required by the caller.
        expected: ContentHash,
        /// Identity reproduced by the receipt.
        actual: ContentHash,
    },
}

impl fmt::Display for PropertyUsageReceiptError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnsupportedSchemaVersion { found, supported } => write!(
                f,
                "property-usage receipt schema {found} is unsupported; this build requires {supported}"
            ),
            Self::InvalidField { field, detail } => {
                write!(
                    f,
                    "property-usage receipt field '{field}' refused: {detail}"
                )
            }
            Self::ResourceLimit {
                resource,
                limit,
                observed,
            } => write!(
                f,
                "property-usage receipt resource '{resource}' exceeds {limit} (observed {observed})"
            ),
            Self::UnknownPolicyTag { tag } => {
                write!(f, "property-usage receipt policy tag '{tag}' is unknown")
            }
            Self::UnknownDecisionTag { tag, at } => write!(
                f,
                "property-usage receipt decision tag {tag} is unknown at byte {at}"
            ),
            Self::Malformed { at, detail } => {
                write!(f, "malformed property-usage receipt at byte {at}: {detail}")
            }
            Self::IdentityMismatch { expected, actual } => write!(
                f,
                "property-usage receipt identity mismatch: encoded {}, reconstructed {}",
                expected.to_hex(),
                actual.to_hex()
            ),
            Self::ExternalIdentityMismatch { expected, actual } => write!(
                f,
                "property-usage receipt external identity mismatch: required {}, reconstructed {}",
                expected.to_hex(),
                actual.to_hex()
            ),
        }
    }
}

impl std::error::Error for PropertyUsageReceiptError {}

/// A named-axis query point ("T" → 293.15, "normal_pressure" → 2.0e5).
/// Axis names match [`fs_evidence::ValidityDomain`] axis names.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct QueryPoint {
    axes: BTreeMap<String, f64>,
}

impl QueryPoint {
    /// An empty point (only unconstrained claims can answer it).
    #[must_use]
    pub fn new() -> QueryPoint {
        QueryPoint::default()
    }

    /// Set one named axis.
    ///
    /// # Errors
    /// [`MatDbError::NonFiniteQueryPoint`] for a non-finite coordinate.
    pub fn with(mut self, axis: impl Into<String>, value: f64) -> Result<QueryPoint, MatDbError> {
        let axis = axis.into();
        if !value.is_finite() {
            return Err(MatDbError::NonFiniteQueryPoint {
                axis,
                bits: value.to_bits(),
            });
        }
        self.axes.insert(axis, value);
        Ok(self)
    }

    /// The named coordinates.
    #[must_use]
    pub fn axes(&self) -> &BTreeMap<String, f64> {
        &self.axes
    }
}

/// How the answer chooses among in-domain candidate claims. Fusion is
/// explicit: no policy invents a canonical value from disagreeing
/// claims — ambiguity is a typed refusal naming the candidates.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SelectionPolicy {
    /// Exactly one in-domain claim may exist; two or more refuse.
    SingleClaimOnly,
    /// Observation-backed claims outrank citation-only claims; the
    /// surviving set must still be a singleton.
    PreferObservationBacked,
}

impl SelectionPolicy {
    /// Stable receipt tag.
    #[must_use]
    pub fn tag(self) -> &'static str {
        match self {
            SelectionPolicy::SingleClaimOnly => "single-claim-only",
            SelectionPolicy::PreferObservationBacked => "prefer-observation-backed",
        }
    }

    /// The policy a receipt tag names. The pinned-claim tag is NOT a
    /// `SelectionPolicy` (its pin travels in the receipt's `selected`
    /// field), so this refuses it; receipt admission uses
    /// [`admitted_policy_tag`].
    ///
    /// # Errors
    /// [`MatDbError::UnknownPolicyTag`] for a tag no policy owns.
    pub fn from_tag(tag: &str) -> Result<SelectionPolicy, MatDbError> {
        match tag {
            "single-claim-only" => Ok(SelectionPolicy::SingleClaimOnly),
            "prefer-observation-backed" => Ok(SelectionPolicy::PreferObservationBacked),
            other => Err(MatDbError::UnknownPolicyTag {
                tag: other.to_string(),
            }),
        }
    }
}

/// Receipt tag of the pinned-claim selection: the caller names the exact
/// claim id, so a card with coexisting conflicting claims resolves only by
/// an EXPLICIT recorded choice — never an auto-pick. The pin itself is the
/// receipt's `selected` field; replay authority additionally requires the
/// verifier to cross-check `selected` against the externally recorded pin,
/// exactly as it must trust-check the query point.
pub const PINNED_CLAIM_POLICY_TAG: &str = "pinned-claim";

/// The closed set of policy tags a portable receipt may carry, mapped to
/// their static spellings. `None` is a foreign tag and refuses upstream.
fn admitted_policy_tag(tag: &str) -> Option<&'static str> {
    match tag {
        "single-claim-only" => Some("single-claim-only"),
        "prefer-observation-backed" => Some("prefer-observation-backed"),
        PINNED_CLAIM_POLICY_TAG => Some(PINNED_CLAIM_POLICY_TAG),
        _ => None,
    }
}

/// How one query selects among in-domain candidates: a named policy or an
/// explicit caller-supplied claim pin.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ClaimSelection {
    /// Selection under a named [`SelectionPolicy`].
    Policy(SelectionPolicy),
    /// Selection pinned to one exact claim id.
    Pinned(ClaimId),
}

impl ClaimSelection {
    fn tag(self) -> &'static str {
        match self {
            ClaimSelection::Policy(policy) => policy.tag(),
            ClaimSelection::Pinned(_) => PINNED_CLAIM_POLICY_TAG,
        }
    }
}

/// How the selected claim was evaluated at the point (successful paths
/// only — extrapolation never succeeds, it refuses).
#[derive(Debug, Clone, PartialEq)]
pub enum EvaluationDecision {
    /// A scalar claim valid across its whole validity box.
    ConstantWithinValidity,
    /// A single tabulated scalar (no abscissa involved).
    ExactScalar,
    /// A tabulated value hit exactly (bit-equal abscissa).
    ExactTabulated {
        /// The matched abscissa.
        at: f64,
    },
    /// Piecewise-linear interpolation strictly inside the knot span.
    LinearInside {
        /// Left bracketing knot abscissa.
        x_lo: f64,
        /// Right bracketing knot abscissa.
        x_hi: f64,
    },
}

impl EvaluationDecision {
    fn encode(&self, encoder: &mut ReceiptEncoder) {
        match self {
            EvaluationDecision::ConstantWithinValidity => {
                encoder.u8(DECISION_CONSTANT_WITHIN_VALIDITY);
            }
            EvaluationDecision::ExactScalar => {
                encoder.u8(DECISION_EXACT_SCALAR);
            }
            EvaluationDecision::ExactTabulated { at } => {
                encoder.u8(DECISION_EXACT_TABULATED);
                encoder.f64(*at);
            }
            EvaluationDecision::LinearInside { x_lo, x_hi } => {
                encoder.u8(DECISION_LINEAR_INSIDE);
                encoder.f64(*x_lo);
                encoder.f64(*x_hi);
            }
        }
    }

    fn decode(reader: &mut ReceiptReader<'_>) -> Result<Self, PropertyUsageReceiptError> {
        let tag = reader.u8()?;
        match tag {
            DECISION_CONSTANT_WITHIN_VALIDITY => Ok(Self::ConstantWithinValidity),
            DECISION_EXACT_SCALAR => Ok(Self::ExactScalar),
            DECISION_EXACT_TABULATED => Ok(Self::ExactTabulated { at: reader.f64()? }),
            DECISION_LINEAR_INSIDE => Ok(Self::LinearInside {
                x_lo: reader.f64()?,
                x_hi: reader.f64()?,
            }),
            tag => Err(PropertyUsageReceiptError::UnknownDecisionTag {
                tag,
                at: reader.position(),
            }),
        }
    }
}

/// The evaluated sample a query returns (inside `Evidence<_>`).
#[derive(Debug, Clone, PartialEq)]
pub struct PropertySample {
    /// The evaluated SI value.
    pub value: f64,
    /// The value's dimensions.
    pub dims: Dims,
    /// The stated uncertainty model it inherits from its claim.
    pub uncertainty: UncertaintyModel,
}

/// The receipt every answer carries: what was asked, what was
/// considered, what was selected under which policy, how it was
/// evaluated, and which sources are load-bearing.
///
/// V2 deliberately does not claim query-time propagation of joint-correlation
/// blocks, retention of a complete evaluated sample, or tensor/frame-transform
/// semantics. Those fields require a later schema/domain rotation rather than
/// being smuggled into an opaque extension payload.
#[derive(Debug, Clone, PartialEq)]
pub struct PropertyUsageReceipt {
    /// Closed portable receipt schema. Only
    /// [`PROPERTY_USAGE_RECEIPT_SCHEMA_VERSION`] is admitted.
    pub schema_version: u32,
    /// The property name queried.
    pub property: String,
    /// The query point (named axes, canonical order).
    pub query_point: Vec<(String, f64)>,
    /// Every claim whose key matched the name, in insertion order —
    /// including out-of-domain claims (the receipt shows what was NOT
    /// eligible, so a narrow answer cannot masquerade as consensus).
    pub considered: Vec<ClaimId>,
    /// The in-domain candidates after validity filtering.
    pub in_domain: Vec<ClaimId>,
    /// The selected claim.
    pub selected: ClaimId,
    /// The selection policy's stable tag.
    pub policy: &'static str,
    /// How the value was produced.
    pub decision: EvaluationDecision,
    /// Whether the selected claim is observation-backed (specimen and
    /// process context exist). Citation-only answers can never be
    /// Validated-class downstream.
    pub observation_backed: bool,
    /// The evaluator's semantic version.
    pub evaluator_version: u32,
    /// Content hashes of the selected claim and its observations.
    pub source_hashes: Vec<ContentHash>,
}

/// Exhaustive owner-type classifier for the portable receipt identity. Adding
/// a receipt field must make identity governance fail until it is deliberately
/// classified and bound.
#[allow(dead_code)]
fn classify_property_usage_receipt_identity_fields(source: &PropertyUsageReceipt) {
    let PropertyUsageReceipt {
        schema_version,
        property,
        query_point,
        considered,
        in_domain,
        selected,
        policy,
        decision,
        observation_backed,
        evaluator_version,
        source_hashes,
    } = source;
    let _ = (
        schema_version,
        property,
        query_point,
        considered,
        in_domain,
        selected,
        policy,
        decision,
        observation_backed,
        evaluator_version,
        source_hashes,
    );
}

/// Owner-local declaration consumed by `xtask check-identities`.
#[allow(dead_code)]
pub const PROPERTY_USAGE_RECEIPT_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-matdb:property-usage-receipt",
    "version_const=PROPERTY_USAGE_RECEIPT_IDENTITY_VERSION",
    "version=2",
    "domain=org.frankensim.fs-matdb.property-usage-receipt.v2",
    "domain_const=PROPERTY_USAGE_RECEIPT_IDENTITY_DOMAIN",
    "encoder=PropertyUsageReceipt::try_content_hash",
    "encoder_helpers=PropertyUsageReceipt::content_hash,PropertyUsageReceipt::identity_preimage,property_usage_receipt_hash,EvaluationDecision::encode,ReceiptEncoder::new,ReceiptEncoder::bytes,ReceiptEncoder::u8,ReceiptEncoder::u32,ReceiptEncoder::u64,ReceiptEncoder::count,ReceiptEncoder::string,ReceiptEncoder::hash,ReceiptEncoder::f64,ReceiptEncoder::boolean,ReceiptEncoder::finish",
    "schema_constants=PROPERTY_USAGE_RECEIPT_IDENTITY_VERSION,PROPERTY_USAGE_RECEIPT_IDENTITY_DOMAIN,PROPERTY_USAGE_RECEIPT_SCHEMA_VERSION,PROPERTY_USAGE_RECEIPT_MAGIC,FIELD_PROPERTY,FIELD_QUERY_POINT,FIELD_CONSIDERED,FIELD_IN_DOMAIN,FIELD_SELECTED,FIELD_POLICY,FIELD_DECISION,FIELD_OBSERVATION_BACKED,FIELD_EVALUATOR_VERSION,FIELD_SOURCE_HASHES,DECISION_CONSTANT_WITHIN_VALIDITY,DECISION_EXACT_SCALAR,DECISION_EXACT_TABULATED,DECISION_LINEAR_INSIDE,MATDB_EVALUATOR_VERSION,MAX_PROPERTY_USAGE_RECEIPT_BYTES,MAX_PROPERTY_USAGE_PROPERTY_BYTES,MAX_PROPERTY_USAGE_QUERY_AXES,MAX_PROPERTY_USAGE_AXIS_BYTES,MAX_PROPERTY_USAGE_CLAIM_IDS,MAX_PROPERTY_USAGE_POLICY_BYTES,MAX_PROPERTY_USAGE_SOURCE_HASHES",
    "schema_functions=PropertyUsageReceipt::validate_portable,PropertyUsageReceipt::to_bytes,PropertyUsageReceipt::from_bytes,PropertyUsageReceipt::from_bytes_verified,EvaluationDecision::decode,SelectionPolicy::tag,SelectionPolicy::from_tag,admitted_policy_tag,ClaimSelection::tag,query_point_exact_eq,evaluation_decision_exact_eq,ReceiptReader::new,ReceiptReader::take,ReceiptReader::require_remaining_items,ReceiptReader::expect,ReceiptReader::expect_tag,ReceiptReader::u8,ReceiptReader::u32,ReceiptReader::u64,ReceiptReader::count,ReceiptReader::string,ReceiptReader::hash,ReceiptReader::f64,ReceiptReader::boolean,ReceiptReader::finish,invalid_field,resource_limit,check_bytes,check_count,crates/fs-blake3/src/lib.rs#hash_domain",
    "schema_types=crates/fs-matdb/src/lib.rs#ClaimSet,crates/fs-matdb/src/lib.rs#InterpolationPolicy,crates/fs-matdb/src/lib.rs#MatDbError,crates/fs-matdb/src/lib.rs#PropertyClaim,crates/fs-matdb/src/lib.rs#PropertyValue,crates/fs-matdb/src/lib.rs#UncertaintyModel,crates/fs-matdb/src/query.rs#ClaimSelection",
    "schema_dependencies=none",
    "digest=blake3-256-domain-separated",
    "encoding=canonical-transport-exact-bits",
    "sources=PropertyUsageReceipt",
    "source_fields=PropertyUsageReceipt.schema_version:semantic,PropertyUsageReceipt.property:semantic,PropertyUsageReceipt.query_point:semantic,PropertyUsageReceipt.considered:semantic,PropertyUsageReceipt.in_domain:semantic,PropertyUsageReceipt.selected:semantic,PropertyUsageReceipt.policy:semantic,PropertyUsageReceipt.decision:semantic,PropertyUsageReceipt.observation_backed:semantic,PropertyUsageReceipt.evaluator_version:semantic,PropertyUsageReceipt.source_hashes:semantic",
    "source_bindings=PropertyUsageReceipt.schema_version>wire-schema-version,PropertyUsageReceipt.property>property-byte-count+property-utf8,PropertyUsageReceipt.query_point>query-point-count+query-point-order+query-axis-byte-count+query-axis-utf8+query-coordinate-f64-exact-bits,PropertyUsageReceipt.considered>considered-count+considered-order+considered-claim-id,PropertyUsageReceipt.in_domain>in-domain-count+in-domain-order+in-domain-claim-id,PropertyUsageReceipt.selected>selected-claim-id,PropertyUsageReceipt.policy>policy-byte-count+policy-utf8,PropertyUsageReceipt.decision>decision-tag+decision-at-f64-exact-bits+decision-x-lo-f64-exact-bits+decision-x-hi-f64-exact-bits,PropertyUsageReceipt.observation_backed>observation-backed-flag,PropertyUsageReceipt.evaluator_version>evaluator-version,PropertyUsageReceipt.source_hashes>source-hash-count+source-hash-order+source-hash",
    "external_semantic_fields=identity-domain,identity-version,wire-magic,canonical-field-order,field-tag-u8,length-count-u64-le,fixed-numeric-little-endian,in-band-identity",
    "semantic_fields=identity-domain,identity-version,wire-magic,canonical-field-order,field-tag-u8,length-count-u64-le,fixed-numeric-little-endian,in-band-identity,wire-schema-version,property-byte-count,property-utf8,query-point-count,query-point-order,query-axis-byte-count,query-axis-utf8,query-coordinate-f64-exact-bits,considered-count,considered-order,considered-claim-id,in-domain-count,in-domain-order,in-domain-claim-id,selected-claim-id,policy-byte-count,policy-utf8,decision-tag,decision-at-f64-exact-bits,decision-x-lo-f64-exact-bits,decision-x-hi-f64-exact-bits,observation-backed-flag,evaluator-version,source-hash-count,source-hash-order,source-hash",
    "excluded_fields=none",
    "consumers=PropertyUsageReceipt::try_content_hash,PropertyUsageReceipt::content_hash,PropertyUsageReceipt::to_bytes,PropertyUsageReceipt::from_bytes,PropertyUsageReceipt::from_bytes_verified,ClaimSet::query,ClaimSet::query_pinned,ClaimSet::verify_receipt,MaterialAnswer",
    "mutations=identity-domain:crates/fs-matdb/tests/query.rs#property_usage_receipt_v2_round_trips_and_replays_exactly,identity-version:crates/fs-matdb/tests/query.rs#property_usage_receipt_v2_round_trips_and_replays_exactly,wire-magic:crates/fs-matdb/tests/query.rs#property_usage_receipt_v2_decoder_fails_closed,canonical-field-order:crates/fs-matdb/tests/query.rs#property_usage_receipt_v2_round_trips_and_replays_exactly,field-tag-u8:crates/fs-matdb/tests/query.rs#property_usage_receipt_v2_decoder_fails_closed,length-count-u64-le:crates/fs-matdb/tests/query.rs#property_usage_receipt_v2_binds_collection_boundaries,fixed-numeric-little-endian:crates/fs-matdb/tests/query.rs#property_usage_receipt_v2_round_trips_and_replays_exactly,in-band-identity:crates/fs-matdb/tests/query.rs#property_usage_receipt_v2_decoder_fails_closed,wire-schema-version:crates/fs-matdb/tests/query.rs#property_usage_receipt_identity_fields_move_independently,property-byte-count:crates/fs-matdb/tests/query.rs#property_usage_receipt_identity_fields_move_independently,property-utf8:crates/fs-matdb/tests/query.rs#property_usage_receipt_identity_fields_move_independently,query-point-count:crates/fs-matdb/tests/query.rs#property_usage_receipt_identity_fields_move_independently,query-point-order:crates/fs-matdb/tests/query.rs#property_usage_receipt_identity_fields_move_independently,query-axis-byte-count:crates/fs-matdb/tests/query.rs#property_usage_receipt_identity_fields_move_independently,query-axis-utf8:crates/fs-matdb/tests/query.rs#property_usage_receipt_identity_fields_move_independently,query-coordinate-f64-exact-bits:crates/fs-matdb/tests/query.rs#property_usage_receipt_identity_fields_move_independently,considered-count:crates/fs-matdb/tests/query.rs#property_usage_receipt_v2_binds_collection_boundaries,considered-order:crates/fs-matdb/tests/query.rs#property_usage_receipt_identity_fields_move_independently,considered-claim-id:crates/fs-matdb/tests/query.rs#property_usage_receipt_identity_fields_move_independently,in-domain-count:crates/fs-matdb/tests/query.rs#property_usage_receipt_v2_binds_collection_boundaries,in-domain-order:crates/fs-matdb/tests/query.rs#property_usage_receipt_identity_fields_move_independently,in-domain-claim-id:crates/fs-matdb/tests/query.rs#property_usage_receipt_identity_fields_move_independently,selected-claim-id:crates/fs-matdb/tests/query.rs#property_usage_receipt_identity_fields_move_independently,policy-byte-count:crates/fs-matdb/tests/query.rs#property_usage_receipt_identity_fields_move_independently,policy-utf8:crates/fs-matdb/tests/query.rs#property_usage_receipt_identity_fields_move_independently,decision-tag:crates/fs-matdb/tests/query.rs#property_usage_receipt_identity_fields_move_independently,decision-at-f64-exact-bits:crates/fs-matdb/tests/query.rs#property_usage_receipt_identity_fields_move_independently,decision-x-lo-f64-exact-bits:crates/fs-matdb/tests/query.rs#property_usage_receipt_identity_fields_move_independently,decision-x-hi-f64-exact-bits:crates/fs-matdb/tests/query.rs#property_usage_receipt_identity_fields_move_independently,observation-backed-flag:crates/fs-matdb/tests/query.rs#property_usage_receipt_identity_fields_move_independently,evaluator-version:crates/fs-matdb/tests/query.rs#property_usage_receipt_identity_fields_move_independently,source-hash-count:crates/fs-matdb/tests/query.rs#property_usage_receipt_identity_fields_move_independently,source-hash-order:crates/fs-matdb/tests/query.rs#property_usage_receipt_identity_fields_move_independently,source-hash:crates/fs-matdb/tests/query.rs#property_usage_receipt_identity_fields_move_independently",
    "nonsemantic_mutations=none",
    "field_guard=classify_property_usage_receipt_identity_fields",
    "transport_guard=PropertyUsageReceipt::from_bytes_verified",
    "version_guard=crates/fs-matdb/tests/query.rs#property_usage_receipt_v2_round_trips_and_replays_exactly",
    "coupling_surface=fs-matdb:property-usage-receipt",
];

impl PropertyUsageReceipt {
    /// Canonical v2 receipt identity over every field and every collection
    /// boundary. Identity alone does not prove that the receipt replays against
    /// a particular [`ClaimSet`]; call [`ClaimSet::verify_receipt`] for that
    /// semantic check.
    #[must_use]
    pub fn content_hash(&self) -> ContentHash {
        property_usage_receipt_hash(
            PROPERTY_USAGE_RECEIPT_IDENTITY_VERSION,
            PROPERTY_USAGE_RECEIPT_IDENTITY_DOMAIN,
            &self.identity_preimage(),
        )
    }

    /// Admit this public-field value to the portable v2 profile before
    /// minting its authoritative identity. [`Self::content_hash`] remains an
    /// infallible structural digest for in-memory mutation diagnostics; only
    /// this method and [`Self::to_bytes`] establish portable admission.
    pub fn try_content_hash(&self) -> Result<ContentHash, PropertyUsageReceiptError> {
        self.validate_portable()?;
        Ok(self.content_hash())
    }

    /// Encode the only supported portable receipt schema.
    ///
    /// The fixed-width in-band identity is appended after the canonical field
    /// stream. Shape/resource admission happens before a second copy of any
    /// caller-owned variable field is allocated.
    pub fn to_bytes(&self) -> Result<Vec<u8>, PropertyUsageReceiptError> {
        let identity = self.try_content_hash()?;
        let mut bytes = self.identity_preimage();
        bytes.extend_from_slice(identity.as_bytes());
        if bytes.len() > MAX_PROPERTY_USAGE_RECEIPT_BYTES {
            return Err(resource_limit(
                "receipt-bytes",
                MAX_PROPERTY_USAGE_RECEIPT_BYTES,
                bytes.len(),
            ));
        }
        Ok(bytes)
    }

    /// Decode, structurally admit, hash-verify, and byte-reproduce a portable
    /// v2 receipt. This does not replay it against a material claim set.
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, PropertyUsageReceiptError> {
        if bytes.len() > MAX_PROPERTY_USAGE_RECEIPT_BYTES {
            return Err(resource_limit(
                "receipt-bytes",
                MAX_PROPERTY_USAGE_RECEIPT_BYTES,
                bytes.len(),
            ));
        }
        let mut reader = ReceiptReader::new(bytes);
        reader.expect(PROPERTY_USAGE_RECEIPT_MAGIC, "wire magic")?;
        let schema_version = reader.u32()?;
        if schema_version != PROPERTY_USAGE_RECEIPT_SCHEMA_VERSION {
            return Err(PropertyUsageReceiptError::UnsupportedSchemaVersion {
                found: schema_version,
                supported: PROPERTY_USAGE_RECEIPT_SCHEMA_VERSION,
            });
        }

        reader.expect_tag(FIELD_PROPERTY, "property")?;
        let property = reader.string("property-bytes", MAX_PROPERTY_USAGE_PROPERTY_BYTES)?;

        reader.expect_tag(FIELD_QUERY_POINT, "query-point")?;
        let query_count = reader.count("query-axes", MAX_PROPERTY_USAGE_QUERY_AXES)?;
        reader.require_remaining_items(query_count, 16, "query axes")?;
        let mut query_point = Vec::with_capacity(query_count);
        for _ in 0..query_count {
            let axis = reader.string("axis-bytes", MAX_PROPERTY_USAGE_AXIS_BYTES)?;
            let value = reader.f64()?;
            query_point.push((axis, value));
        }

        reader.expect_tag(FIELD_CONSIDERED, "considered")?;
        let considered_count =
            reader.count("considered-claim-ids", MAX_PROPERTY_USAGE_CLAIM_IDS)?;
        reader.require_remaining_items(considered_count, 32, "considered claim ids")?;
        let mut considered = Vec::with_capacity(considered_count);
        for _ in 0..considered_count {
            considered.push(ClaimId(reader.hash()?));
        }

        reader.expect_tag(FIELD_IN_DOMAIN, "in-domain")?;
        let in_domain_count = reader.count("in-domain-claim-ids", MAX_PROPERTY_USAGE_CLAIM_IDS)?;
        reader.require_remaining_items(in_domain_count, 32, "in-domain claim ids")?;
        let mut in_domain = Vec::with_capacity(in_domain_count);
        for _ in 0..in_domain_count {
            in_domain.push(ClaimId(reader.hash()?));
        }

        reader.expect_tag(FIELD_SELECTED, "selected")?;
        let selected = ClaimId(reader.hash()?);

        reader.expect_tag(FIELD_POLICY, "policy")?;
        let policy_text = reader.string("policy-bytes", MAX_PROPERTY_USAGE_POLICY_BYTES)?;
        let Some(policy) = admitted_policy_tag(&policy_text) else {
            return Err(PropertyUsageReceiptError::UnknownPolicyTag { tag: policy_text });
        };

        reader.expect_tag(FIELD_DECISION, "decision")?;
        let decision = EvaluationDecision::decode(&mut reader)?;

        reader.expect_tag(FIELD_OBSERVATION_BACKED, "observation-backed")?;
        let observation_backed = reader.boolean("observation-backed")?;

        reader.expect_tag(FIELD_EVALUATOR_VERSION, "evaluator-version")?;
        let evaluator_version = reader.u32()?;

        reader.expect_tag(FIELD_SOURCE_HASHES, "source-hashes")?;
        let source_count = reader.count("source-hashes", MAX_PROPERTY_USAGE_SOURCE_HASHES)?;
        reader.require_remaining_items(source_count, 32, "source hashes")?;
        let mut source_hashes = Vec::with_capacity(source_count);
        for _ in 0..source_count {
            source_hashes.push(reader.hash()?);
        }

        let encoded_identity = reader.hash()?;
        reader.finish()?;

        let receipt = Self {
            schema_version,
            property,
            query_point,
            considered,
            in_domain,
            selected,
            policy,
            decision,
            observation_backed,
            evaluator_version,
            source_hashes,
        };
        receipt.validate_portable()?;
        let actual = receipt.content_hash();
        if encoded_identity != actual {
            return Err(PropertyUsageReceiptError::IdentityMismatch {
                expected: encoded_identity,
                actual,
            });
        }
        let reproduced = receipt.to_bytes()?;
        if reproduced != bytes {
            return Err(PropertyUsageReceiptError::Malformed {
                at: bytes.len(),
                detail: "decoded fields do not reproduce the canonical byte stream".to_string(),
            });
        }
        Ok(receipt)
    }

    /// Decode with an external content identity pinned by a ledger, package,
    /// or other retention authority.
    pub fn from_bytes_verified(
        bytes: &[u8],
        expected: ContentHash,
    ) -> Result<Self, PropertyUsageReceiptError> {
        let receipt = Self::from_bytes(bytes)?;
        let actual = receipt.content_hash();
        if actual != expected {
            return Err(PropertyUsageReceiptError::ExternalIdentityMismatch { expected, actual });
        }
        Ok(receipt)
    }

    fn identity_preimage(&self) -> Vec<u8> {
        let mut encoder = ReceiptEncoder::new();
        encoder.bytes(PROPERTY_USAGE_RECEIPT_MAGIC);
        encoder.u32(self.schema_version);

        encoder.u8(FIELD_PROPERTY);
        encoder.string(&self.property);

        encoder.u8(FIELD_QUERY_POINT);
        encoder.count(self.query_point.len());
        for (axis, value) in &self.query_point {
            encoder.string(axis);
            encoder.f64(*value);
        }

        encoder.u8(FIELD_CONSIDERED);
        encoder.count(self.considered.len());
        for id in &self.considered {
            encoder.hash(id.0);
        }

        encoder.u8(FIELD_IN_DOMAIN);
        encoder.count(self.in_domain.len());
        for id in &self.in_domain {
            encoder.hash(id.0);
        }

        encoder.u8(FIELD_SELECTED);
        encoder.hash(self.selected.0);

        encoder.u8(FIELD_POLICY);
        encoder.string(self.policy);

        encoder.u8(FIELD_DECISION);
        self.decision.encode(&mut encoder);

        encoder.u8(FIELD_OBSERVATION_BACKED);
        encoder.boolean(self.observation_backed);

        encoder.u8(FIELD_EVALUATOR_VERSION);
        encoder.u32(self.evaluator_version);

        encoder.u8(FIELD_SOURCE_HASHES);
        encoder.count(self.source_hashes.len());
        for hash in &self.source_hashes {
            encoder.hash(*hash);
        }
        encoder.finish()
    }

    #[allow(clippy::too_many_lines)] // the closed portable profile IS this one straight-line check sequence
    fn validate_portable(&self) -> Result<(), PropertyUsageReceiptError> {
        if self.schema_version != PROPERTY_USAGE_RECEIPT_SCHEMA_VERSION {
            return Err(PropertyUsageReceiptError::UnsupportedSchemaVersion {
                found: self.schema_version,
                supported: PROPERTY_USAGE_RECEIPT_SCHEMA_VERSION,
            });
        }
        if self.property.trim().is_empty() {
            return Err(invalid_field("property", "must not be blank"));
        }
        check_bytes(
            "property-bytes",
            self.property.len(),
            MAX_PROPERTY_USAGE_PROPERTY_BYTES,
        )?;
        check_count(
            "query-axes",
            self.query_point.len(),
            MAX_PROPERTY_USAGE_QUERY_AXES,
        )?;
        let mut prior_axis: Option<&str> = None;
        for (axis, value) in &self.query_point {
            if axis.trim().is_empty() {
                return Err(invalid_field("query-axis", "must not be blank"));
            }
            check_bytes("axis-bytes", axis.len(), MAX_PROPERTY_USAGE_AXIS_BYTES)?;
            if prior_axis.is_some_and(|prior| prior >= axis.as_str()) {
                return Err(invalid_field(
                    "query-point",
                    "axis names must be strictly increasing and unique",
                ));
            }
            if !value.is_finite() {
                return Err(invalid_field(
                    "query-coordinate",
                    format!("must be finite (bits {:#018x})", value.to_bits()),
                ));
            }
            prior_axis = Some(axis);
        }
        check_count(
            "considered-claim-ids",
            self.considered.len(),
            MAX_PROPERTY_USAGE_CLAIM_IDS,
        )?;
        if self.considered.is_empty() {
            return Err(invalid_field(
                "considered",
                "must retain at least the selected claim",
            ));
        }
        let considered_unique: BTreeSet<_> = self.considered.iter().copied().collect();
        if considered_unique.len() != self.considered.len() {
            return Err(invalid_field("considered", "claim ids must be unique"));
        }
        check_count(
            "in-domain-claim-ids",
            self.in_domain.len(),
            MAX_PROPERTY_USAGE_CLAIM_IDS,
        )?;
        if self.in_domain.is_empty() {
            return Err(invalid_field(
                "in-domain",
                "must retain at least the selected claim",
            ));
        }
        let mut considered_cursor = 0;
        for candidate in &self.in_domain {
            let Some(relative) = self.considered[considered_cursor..]
                .iter()
                .position(|considered| considered == candidate)
            else {
                return Err(invalid_field(
                    "in-domain",
                    "must be an order-preserving subsequence of considered claims",
                ));
            };
            considered_cursor += relative + 1;
        }
        if !self.in_domain.contains(&self.selected) {
            return Err(invalid_field(
                "selected",
                "must be present in the in-domain collection",
            ));
        }
        check_bytes(
            "policy-bytes",
            self.policy.len(),
            MAX_PROPERTY_USAGE_POLICY_BYTES,
        )?;
        if admitted_policy_tag(self.policy).is_none() {
            return Err(PropertyUsageReceiptError::UnknownPolicyTag {
                tag: self.policy.to_string(),
            });
        }
        match &self.decision {
            EvaluationDecision::ExactTabulated { at } if !at.is_finite() => {
                return Err(invalid_field(
                    "decision.at",
                    format!("must be finite (bits {:#018x})", at.to_bits()),
                ));
            }
            EvaluationDecision::LinearInside { x_lo, x_hi }
                if !x_lo.is_finite() || !x_hi.is_finite() || x_lo >= x_hi =>
            {
                return Err(invalid_field(
                    "decision.linear-inside",
                    "bounds must be finite and strictly increasing",
                ));
            }
            _ => {}
        }
        check_count(
            "source-hashes",
            self.source_hashes.len(),
            MAX_PROPERTY_USAGE_SOURCE_HASHES,
        )?;
        if self.source_hashes.first().copied() != Some(self.selected.0) {
            return Err(invalid_field(
                "source-hashes",
                "first source hash must be the selected claim identity",
            ));
        }
        if self.observation_backed != (self.source_hashes.len() > 1) {
            return Err(invalid_field(
                "observation-backed",
                "must agree with the presence of observation source hashes",
            ));
        }
        Ok(())
    }
}

fn property_usage_receipt_hash(version: u32, domain: &str, wire_preimage: &[u8]) -> ContentHash {
    let mut identity_preimage = Vec::with_capacity(4 + wire_preimage.len());
    identity_preimage.extend_from_slice(&version.to_le_bytes());
    identity_preimage.extend_from_slice(wire_preimage);
    hash_domain(domain, &identity_preimage)
}

fn invalid_field(field: &'static str, detail: impl Into<String>) -> PropertyUsageReceiptError {
    PropertyUsageReceiptError::InvalidField {
        field,
        detail: detail.into(),
    }
}

fn resource_limit(
    resource: &'static str,
    limit: usize,
    observed: usize,
) -> PropertyUsageReceiptError {
    PropertyUsageReceiptError::ResourceLimit {
        resource,
        limit,
        observed: u64::try_from(observed).unwrap_or(u64::MAX),
    }
}

fn check_bytes(
    resource: &'static str,
    observed: usize,
    limit: usize,
) -> Result<(), PropertyUsageReceiptError> {
    if observed > limit {
        Err(resource_limit(resource, limit, observed))
    } else {
        Ok(())
    }
}

fn check_count(
    resource: &'static str,
    observed: usize,
    limit: usize,
) -> Result<(), PropertyUsageReceiptError> {
    check_bytes(resource, observed, limit)
}

struct ReceiptEncoder {
    bytes: Vec<u8>,
}

impl ReceiptEncoder {
    fn new() -> Self {
        Self {
            bytes: Vec::with_capacity(512),
        }
    }

    fn bytes(&mut self, bytes: &[u8]) {
        self.bytes.extend_from_slice(bytes);
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

    fn count(&mut self, value: usize) {
        self.u64(u64::try_from(value).unwrap_or(u64::MAX));
    }

    fn string(&mut self, value: &str) {
        self.count(value.len());
        self.bytes(value.as_bytes());
    }

    fn hash(&mut self, value: ContentHash) {
        self.bytes(value.as_bytes());
    }

    fn f64(&mut self, value: f64) {
        self.u64(value.to_bits());
    }

    fn boolean(&mut self, value: bool) {
        self.u8(u8::from(value));
    }

    fn finish(self) -> Vec<u8> {
        self.bytes
    }
}

struct ReceiptReader<'a> {
    bytes: &'a [u8],
    cursor: usize,
}

impl<'a> ReceiptReader<'a> {
    fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, cursor: 0 }
    }

    fn position(&self) -> usize {
        self.cursor
    }

    fn malformed(&self, detail: impl Into<String>) -> PropertyUsageReceiptError {
        PropertyUsageReceiptError::Malformed {
            at: self.cursor,
            detail: detail.into(),
        }
    }

    fn take(&mut self, length: usize) -> Result<&'a [u8], PropertyUsageReceiptError> {
        let end = self
            .cursor
            .checked_add(length)
            .ok_or_else(|| self.malformed("byte offset overflow"))?;
        let slice = self
            .bytes
            .get(self.cursor..end)
            .ok_or_else(|| self.malformed(format!("truncated field requires {length} bytes")))?;
        self.cursor = end;
        Ok(slice)
    }

    fn require_remaining_items(
        &self,
        count: usize,
        minimum_width: usize,
        field: &str,
    ) -> Result<(), PropertyUsageReceiptError> {
        let minimum = count
            .checked_mul(minimum_width)
            .ok_or_else(|| self.malformed(format!("{field} byte length overflow")))?;
        if self.bytes.len().saturating_sub(self.cursor) < minimum {
            Err(self.malformed(format!(
                "truncated {field} requires at least {minimum} bytes"
            )))
        } else {
            Ok(())
        }
    }

    fn expect(&mut self, expected: &[u8], field: &str) -> Result<(), PropertyUsageReceiptError> {
        let actual = self.take(expected.len())?;
        if actual == expected {
            Ok(())
        } else {
            Err(self.malformed(format!("invalid {field}")))
        }
    }

    fn expect_tag(&mut self, expected: u8, field: &str) -> Result<(), PropertyUsageReceiptError> {
        let actual = self.u8()?;
        if actual == expected {
            Ok(())
        } else {
            Err(self.malformed(format!(
                "invalid {field} field tag {actual}; expected {expected}"
            )))
        }
    }

    fn u8(&mut self) -> Result<u8, PropertyUsageReceiptError> {
        Ok(self.take(1)?[0])
    }

    fn u32(&mut self) -> Result<u32, PropertyUsageReceiptError> {
        let bytes: [u8; 4] = self
            .take(4)?
            .try_into()
            .map_err(|_| self.malformed("u32 width"))?;
        Ok(u32::from_le_bytes(bytes))
    }

    fn u64(&mut self) -> Result<u64, PropertyUsageReceiptError> {
        let bytes: [u8; 8] = self
            .take(8)?
            .try_into()
            .map_err(|_| self.malformed("u64 width"))?;
        Ok(u64::from_le_bytes(bytes))
    }

    fn count(
        &mut self,
        resource: &'static str,
        limit: usize,
    ) -> Result<usize, PropertyUsageReceiptError> {
        let observed = self.u64()?;
        if observed > limit as u64 {
            return Err(PropertyUsageReceiptError::ResourceLimit {
                resource,
                limit,
                observed,
            });
        }
        usize::try_from(observed).map_err(|_| PropertyUsageReceiptError::ResourceLimit {
            resource,
            limit,
            observed,
        })
    }

    fn string(
        &mut self,
        resource: &'static str,
        limit: usize,
    ) -> Result<String, PropertyUsageReceiptError> {
        let length = self.count(resource, limit)?;
        let start = self.cursor;
        let bytes = self.take(length)?;
        let text =
            std::str::from_utf8(bytes).map_err(|error| PropertyUsageReceiptError::Malformed {
                at: start + error.valid_up_to(),
                detail: format!("{resource} is not UTF-8"),
            })?;
        Ok(text.to_string())
    }

    fn hash(&mut self) -> Result<ContentHash, PropertyUsageReceiptError> {
        let mut hash = [0_u8; 32];
        hash.copy_from_slice(self.take(32)?);
        Ok(ContentHash(hash))
    }

    fn f64(&mut self) -> Result<f64, PropertyUsageReceiptError> {
        Ok(f64::from_bits(self.u64()?))
    }

    fn boolean(&mut self, field: &str) -> Result<bool, PropertyUsageReceiptError> {
        match self.u8()? {
            0 => Ok(false),
            1 => Ok(true),
            tag => Err(self.malformed(format!("invalid {field} boolean tag {tag}"))),
        }
    }

    fn finish(&self) -> Result<(), PropertyUsageReceiptError> {
        if self.cursor == self.bytes.len() {
            Ok(())
        } else {
            Err(self.malformed(format!("{} trailing bytes", self.bytes.len() - self.cursor)))
        }
    }
}

/// A complete material answer: the evidence-carried sample plus its
/// usage receipt.
#[derive(Debug, Clone, PartialEq)]
pub struct MaterialAnswer {
    /// The sample with its honest evidence slices.
    pub evidence: Evidence<PropertySample>,
    /// The usage receipt.
    pub receipt: PropertyUsageReceipt,
}

/// Narrow the in-domain candidates under the caller's selection. Policy
/// arms may keep several candidates (the caller refuses ambiguity); the
/// pinned arm refuses here with the pin-specific diagnosis.
fn select_candidates<'p, 'c>(
    property: &str,
    selection: ClaimSelection,
    considered: &[ClaimId],
    in_domain_pairs: &[&'p (ClaimId, &'c PropertyClaim)],
) -> Result<Vec<&'p (ClaimId, &'c PropertyClaim)>, MatDbError> {
    match selection {
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly) => Ok(in_domain_pairs.to_vec()),
        ClaimSelection::Policy(SelectionPolicy::PreferObservationBacked) => {
            let backed: Vec<_> = in_domain_pairs
                .iter()
                .filter(|(_, claim)| !claim.observations.is_empty())
                .copied()
                .collect();
            if backed.is_empty() {
                Ok(in_domain_pairs.to_vec())
            } else {
                Ok(backed)
            }
        }
        ClaimSelection::Pinned(pinned) => {
            let matched: Vec<_> = in_domain_pairs
                .iter()
                .filter(|(id, _)| *id == pinned)
                .copied()
                .collect();
            if matched.is_empty() {
                return Err(if considered.contains(&pinned) {
                    MatDbError::PinnedClaimOutOfDomain {
                        property: property.to_string(),
                        pinned,
                    }
                } else {
                    MatDbError::PinnedClaimUnknown {
                        property: property.to_string(),
                        pinned,
                    }
                });
            }
            Ok(matched)
        }
    }
}

/// Honest slice mapping. Numerical: the stated band as an ESTIMATE around
/// the value (never Exact/Enclosure — a datum is not interval-certified
/// numerics), or an explicit no-claim for Unstated uncertainty.
/// Statistical: the stated half-width.
fn honest_certificates(
    uncertainty: &UncertaintyModel,
    value: f64,
) -> (NumericalCertificate, StatisticalCertificate) {
    match *uncertainty {
        UncertaintyModel::Unstated => (
            NumericalCertificate::no_claim(),
            StatisticalCertificate::None,
        ),
        UncertaintyModel::HalfWidth {
            half_width,
            confidence,
        } => (
            NumericalCertificate::estimate(value - half_width, value + half_width),
            StatisticalCertificate::HalfWidth {
                half_width,
                confidence,
            },
        ),
        UncertaintyModel::RelativeHalfWidth {
            fraction,
            confidence,
        } => {
            let half_width = fraction * value.abs();
            (
                NumericalCertificate::estimate(value - half_width, value + half_width),
                StatisticalCertificate::HalfWidth {
                    half_width,
                    confidence,
                },
            )
        }
    }
}

impl ClaimSet {
    /// Answer a property query at a point under an explicit selection
    /// policy.
    ///
    /// # Errors
    /// [`MatDbError::UnknownProperty`] when no claim carries the name;
    /// [`MatDbError::NoClaimInDomain`] when the point is outside every
    /// claim's validity (THE extrapolation refusal);
    /// [`MatDbError::AmbiguousSelection`] when the policy cannot narrow
    /// the candidates to one claim (fusion must be explicit);
    /// [`MatDbError::MissingQueryAxis`] / [`MatDbError::OutsideKnotSpan`]
    /// / [`MatDbError::UnsupportedEvaluation`] from curve evaluation.
    pub fn query(
        &self,
        property: &str,
        point: &QueryPoint,
        policy: SelectionPolicy,
    ) -> Result<MaterialAnswer, MatDbError> {
        self.query_selected(property, point, ClaimSelection::Policy(policy))
    }

    /// Answer a property query with selection pinned to one exact claim
    /// id: the EXPLICIT resolution for a card carrying coexisting
    /// conflicting claims. The receipt records the
    /// [`PINNED_CLAIM_POLICY_TAG`] policy with the pin as its `selected`
    /// claim; a verifier must cross-check `selected` against the
    /// externally recorded pin, since the receipt alone cannot prove who
    /// chose it.
    ///
    /// # Errors
    /// [`MatDbError::PinnedClaimUnknown`] when the pin names no claim
    /// under the property; [`MatDbError::PinnedClaimOutOfDomain`] when
    /// the pinned claim does not cover the point (a pin never bypasses
    /// the extrapolation refusal); plus every refusal [`Self::query`]
    /// can return except [`MatDbError::AmbiguousSelection`].
    pub fn query_pinned(
        &self,
        property: &str,
        point: &QueryPoint,
        pinned: ClaimId,
    ) -> Result<MaterialAnswer, MatDbError> {
        self.query_selected(property, point, ClaimSelection::Pinned(pinned))
    }

    fn query_selected(
        &self,
        property: &str,
        point: &QueryPoint,
        selection: ClaimSelection,
    ) -> Result<MaterialAnswer, MatDbError> {
        let considered_pairs = self.claims_for(property);
        if considered_pairs.is_empty() {
            return Err(MatDbError::UnknownProperty {
                property: property.to_string(),
            });
        }
        let considered: Vec<ClaimId> = considered_pairs.iter().map(|(id, _)| *id).collect();
        let in_domain_pairs: Vec<_> = considered_pairs
            .iter()
            .filter(|(_, claim)| claim.validity.contains(point.axes()))
            .collect();
        if in_domain_pairs.is_empty() {
            return Err(MatDbError::NoClaimInDomain {
                property: property.to_string(),
                considered: considered.len(),
            });
        }
        let in_domain: Vec<ClaimId> = in_domain_pairs.iter().map(|(id, _)| *id).collect();
        let selected_pairs = select_candidates(property, selection, &considered, &in_domain_pairs)?;
        if selected_pairs.len() != 1 {
            return Err(MatDbError::AmbiguousSelection {
                property: property.to_string(),
                candidates: selected_pairs.iter().map(|(id, _)| *id).collect(),
            });
        }
        let (selected_id, claim) = selected_pairs[0];
        let (value, decision) = evaluate(&claim.value, claim.interpolation, point)?;

        let mut source_hashes = vec![selected_id.0];
        for observation in &claim.observations {
            source_hashes.push(observation.0);
        }
        let receipt = PropertyUsageReceipt {
            schema_version: PROPERTY_USAGE_RECEIPT_SCHEMA_VERSION,
            property: property.to_string(),
            query_point: point
                .axes()
                .iter()
                .map(|(axis, &v)| (axis.clone(), v))
                .collect(),
            considered,
            in_domain,
            selected: *selected_id,
            policy: selection.tag(),
            decision,
            observation_backed: !claim.observations.is_empty(),
            evaluator_version: MATDB_EVALUATOR_VERSION,
            source_hashes,
        };

        let (numerical, statistical) = honest_certificates(&claim.uncertainty, value);
        let model = ModelEvidence {
            cards: vec![format!("fs-matdb:{property}")],
            assumptions: vec![format!(
                "claim provenance: {} ({})",
                claim.provenance.source, claim.provenance.license
            )],
            validity: claim.validity.clone(),
            discrepancy_rel: 0.0,
            in_domain: true,
        };
        let receipt_identity = receipt
            .try_content_hash()
            .map_err(|error| MatDbError::PropertyUsageReceiptNotPortable { error })?;
        let mut provenance_prefix = [0_u8; 8];
        provenance_prefix.copy_from_slice(&receipt_identity.0[..8]);
        let provenance = ProvenanceHash(u64::from_le_bytes(provenance_prefix));
        let evidence = Evidence {
            value: PropertySample {
                value,
                dims: claim.value.dims(),
                uncertainty: claim.uncertainty.clone(),
            },
            qoi: value,
            numerical,
            statistical,
            model,
            sensitivity: SensitivitySummary::default(),
            provenance,
            adjoint_ref: None,
        };
        Ok(MaterialAnswer { evidence, receipt })
    }
}

impl ClaimSet {
    /// Verify a receipt against this claim set: the receipt-completeness
    /// battery's door (bead 5hmy, PR-5). The query is REPLAYED from the
    /// receipt's own fields and every field must reproduce — a receipt
    /// with any deleted, substituted, or stale field fails with a typed
    /// refusal naming the first divergent field. A receipt that
    /// verifies is exactly as trustworthy as the claim set it is
    /// checked against, no more.
    ///
    /// # Errors
    /// [`MatDbError::EvaluatorVersionDrift`] before any replay (a
    /// receipt from another evaluator version cannot be adjudicated
    /// here); [`MatDbError::PropertyUsageReceiptNotPortable`] when the
    /// public receipt fields violate the closed canonical profile;
    /// [`MatDbError::UnknownPolicyTag`] for a foreign policy tag; the
    /// replay's own refusals (e.g. a tampered query point may refuse as
    /// out-of-domain); and
    /// [`MatDbError::ReceiptMismatch`] naming the field that failed to
    /// reproduce.
    pub fn verify_receipt(&self, receipt: &PropertyUsageReceipt) -> Result<(), MatDbError> {
        if receipt.schema_version != PROPERTY_USAGE_RECEIPT_SCHEMA_VERSION {
            return Err(MatDbError::ReceiptSchemaVersionDrift {
                receipt: receipt.schema_version,
                current: PROPERTY_USAGE_RECEIPT_SCHEMA_VERSION,
            });
        }
        if receipt.evaluator_version != MATDB_EVALUATOR_VERSION {
            return Err(MatDbError::EvaluatorVersionDrift {
                receipt: receipt.evaluator_version,
                current: MATDB_EVALUATOR_VERSION,
            });
        }
        receipt.validate_portable().map_err(|error| match error {
            PropertyUsageReceiptError::UnknownPolicyTag { tag } => {
                MatDbError::UnknownPolicyTag { tag }
            }
            error => MatDbError::PropertyUsageReceiptNotPortable { error },
        })?;
        let selection = if receipt.policy == PINNED_CLAIM_POLICY_TAG {
            // The pin travels as the receipt's own `selected` field; the
            // replay proves that pin resolves and evaluates identically,
            // while WHO pinned it is the verifier's external cross-check.
            ClaimSelection::Pinned(receipt.selected)
        } else {
            ClaimSelection::Policy(SelectionPolicy::from_tag(receipt.policy)?)
        };
        let mut point = QueryPoint::new();
        for (axis, value) in &receipt.query_point {
            point = point.with(axis.clone(), *value)?;
        }
        let replayed = self.query_selected(&receipt.property, &point, selection)?;
        let fresh = &replayed.receipt;
        for (field, matches) in [
            (
                "query-point",
                query_point_exact_eq(&fresh.query_point, &receipt.query_point),
            ),
            ("considered", fresh.considered == receipt.considered),
            ("in_domain", fresh.in_domain == receipt.in_domain),
            ("selected", fresh.selected == receipt.selected),
            (
                "decision",
                evaluation_decision_exact_eq(&fresh.decision, &receipt.decision),
            ),
            (
                "observation_backed",
                fresh.observation_backed == receipt.observation_backed,
            ),
            (
                "source_hashes",
                fresh.source_hashes == receipt.source_hashes,
            ),
        ] {
            if !matches {
                return Err(MatDbError::ReceiptMismatch { field });
            }
        }
        Ok(())
    }
}

fn query_point_exact_eq(left: &[(String, f64)], right: &[(String, f64)]) -> bool {
    left.len() == right.len()
        && left
            .iter()
            .zip(right)
            .all(|((left_axis, left_value), (right_axis, right_value))| {
                left_axis == right_axis && left_value.to_bits() == right_value.to_bits()
            })
}

fn evaluation_decision_exact_eq(left: &EvaluationDecision, right: &EvaluationDecision) -> bool {
    match (left, right) {
        (
            EvaluationDecision::ConstantWithinValidity,
            EvaluationDecision::ConstantWithinValidity,
        )
        | (EvaluationDecision::ExactScalar, EvaluationDecision::ExactScalar) => true,
        (
            EvaluationDecision::ExactTabulated { at: left },
            EvaluationDecision::ExactTabulated { at: right },
        ) => left.to_bits() == right.to_bits(),
        (
            EvaluationDecision::LinearInside {
                x_lo: left_lo,
                x_hi: left_hi,
            },
            EvaluationDecision::LinearInside {
                x_lo: right_lo,
                x_hi: right_hi,
            },
        ) => left_lo.to_bits() == right_lo.to_bits() && left_hi.to_bits() == right_hi.to_bits(),
        _ => false,
    }
}

/// Evaluate a claim payload at a point under its interpolation policy.
fn evaluate(
    value: &PropertyValue,
    policy: InterpolationPolicy,
    point: &QueryPoint,
) -> Result<(f64, EvaluationDecision), MatDbError> {
    match (value, policy) {
        (PropertyValue::Scalar { value, .. }, InterpolationPolicy::ConstantWithinValidity) => {
            Ok((*value, EvaluationDecision::ConstantWithinValidity))
        }
        (PropertyValue::Scalar { value, .. }, InterpolationPolicy::TabulatedOnly) => {
            Ok((*value, EvaluationDecision::ExactScalar))
        }
        (PropertyValue::Scalar { .. }, InterpolationPolicy::LinearInside) => {
            Err(MatDbError::UnsupportedEvaluation {
                reason: "a scalar claim has no knot span to interpolate inside",
            })
        }
        (
            PropertyValue::Curve {
                abscissa, knots, ..
            },
            InterpolationPolicy::LinearInside,
        ) => {
            let x = *point
                .axes()
                .get(abscissa)
                .ok_or_else(|| MatDbError::MissingQueryAxis {
                    axis: abscissa.clone(),
                })?;
            let first = knots[0].0;
            let last = knots[knots.len() - 1].0;
            if x < first || x > last {
                return Err(MatDbError::OutsideKnotSpan {
                    axis: abscissa.clone(),
                    requested: x,
                    lo: first,
                    hi: last,
                });
            }
            if let Some(&(kx, ky)) = knots.iter().find(|&&(kx, _)| kx.to_bits() == x.to_bits()) {
                return Ok((ky, EvaluationDecision::ExactTabulated { at: kx }));
            }
            let window = knots
                .windows(2)
                .find(|w| w[0].0 <= x && x <= w[1].0)
                .expect("span containment guarantees a bracketing window");
            let (x0, y0) = window[0];
            let (x1, y1) = window[1];
            let t = (x - x0) / (x1 - x0);
            let y = y0 + t * (y1 - y0);
            Ok((y, EvaluationDecision::LinearInside { x_lo: x0, x_hi: x1 }))
        }
        (
            PropertyValue::Curve {
                abscissa, knots, ..
            },
            InterpolationPolicy::TabulatedOnly,
        ) => {
            let x = *point
                .axes()
                .get(abscissa)
                .ok_or_else(|| MatDbError::MissingQueryAxis {
                    axis: abscissa.clone(),
                })?;
            knots
                .iter()
                .find(|&&(kx, _)| kx.to_bits() == x.to_bits())
                .map(|&(kx, ky)| (ky, EvaluationDecision::ExactTabulated { at: kx }))
                .ok_or(MatDbError::UnsupportedEvaluation {
                    reason: "tabulated-only claim has no knot at the requested abscissa",
                })
        }
        (PropertyValue::Curve { .. }, InterpolationPolicy::ConstantWithinValidity) => {
            Err(MatDbError::UnsupportedEvaluation {
                reason: "a curve claim cannot be evaluated as a validity-wide constant",
            })
        }
    }
}
