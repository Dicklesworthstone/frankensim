//! WRITE-TIME enforcement of the three-color schema (Proposal 3,
//! bead qmao.1): the [`ColorGraph`] accepts only writes whose claimed
//! color exactly matches what its evidence derives: Estimated leaves
//! enter through a validated weak-source gate; positive leaves are minted
//! from typed certificate or anchoring origins authenticated by an injected
//! source authority (or a separately scoped authenticated source
//! waiver); derived colors are recomputed from their parents. An
//! estimated result CANNOT be written as verified (the laundering
//! refusal). Validated claims are re-checked against the CURRENT
//! execution state and every regime exit AUTO-DEMOTES. Certificate artifacts,
//! admitting policy identities, all demotions, and authenticated
//! operation-bound waivers participate in the node provenance hash and cannot
//! be quietly dropped later.
//!
//! The color enum and pairwise algebra live in fs-evidence (usable by
//! every layer); this module is the HELM-side gatekeeper over
//! already-colored values. Rows are canonical JSON lines ready for the
//! event stream; a dedicated schema table is a CONTRACT no-claim.

use crate::hash::{ContentHash, hash_bytes};
use fs_evidence::{
    COLOR_ALGEBRA_VERSION, Color, ColorPayloadError, ColorRank, Demotion, IntervalOp,
    MAX_COLOR_IDENTITY_BYTES, NumericalCertificate, ValidityDomain,
    color_identity_reason as identity_reason, color_leaf_identity_reason, compose,
    demotion_estimator_identity, regime_demotion, validate_color_payload, verified_from,
};
use std::collections::BTreeMap;

/// Maximum number of direct parents accepted by one color derivation.
/// This also bounds lineage vectors presented to waiver verifiers.
pub const MAX_COLOR_PARENTS: usize = 1_024;

/// Maximum number of distinct historical waiver authorities copied into one
/// node. [`MAX_WAIVER_CLOSURE_BYTES`] independently bounds their aggregate
/// retained payload before cloning and row serialization.
pub const MAX_WAIVER_DEPENDENCIES: usize = 1_024;

/// Maximum aggregate retained bytes for the complete waiver-authority closure
/// copied into one node. The count limit alone is insufficient because each
/// dependency carries a signed color payload and artifact lineage.
pub const MAX_WAIVER_CLOSURE_BYTES: usize = 8 * 1024 * 1024;

/// Maximum UTF-8 byte length of a canonical color-graph node identity.
/// Node names use the shared fs-evidence identity grammar.
pub const MAX_COLOR_NODE_NAME_BYTES: usize = MAX_COLOR_IDENTITY_BYTES;
const _: () = assert!(MAX_COLOR_NODE_NAME_BYTES == MAX_COLOR_IDENTITY_BYTES);

const MAX_WAIVER_REASON_BYTES: usize = 4_096;
const MAX_WAIVER_SIGNATURE_BYTES: usize = 4_096;
const MAX_CLAIMED_COLOR_BYTES: usize = 1_048_576;
/// Maximum number of distinct axes in any admitted Validated color, including
/// the aggregate intersection schema of a multi-parent derivation.
pub const MAX_VALIDITY_AXES: usize = 1_024;
/// Current schema of `color-write` event rows.
pub const COLOR_WRITE_ROW_SCHEMA_VERSION: u32 = 7;
/// Current schema of exact regime-demotion event rows.
pub const COLOR_DEMOTION_ROW_SCHEMA_VERSION: u32 = 1;
const COLOR_NODE_HASH_ENCODING_VERSION: u8 = 9;

fn is_placeholder_token(value: &str) -> bool {
    [
        "-",
        "?",
        "n/a",
        "na",
        "none",
        "not run",
        "pending",
        "placeholder",
        "tbd",
        "todo",
        "unknown",
    ]
    .iter()
    .any(|placeholder| value.eq_ignore_ascii_case(placeholder))
}

fn validate_node_name(name: &str) -> Result<(), ColorWriteError> {
    if let Some(reason) = identity_reason(name) {
        Err(ColorWriteError::InvalidNodeName { reason })
    } else {
        Ok(())
    }
}

/// A human ANNOTATION (ticket, memo, name, rationale). It travels in
/// provenance but AUTHORIZES NOTHING (bead qmao.1.1): presence of
/// caller-created strings is not proof. The only path past a
/// laundering refusal is an authenticated [`WaiverGrant`] through
/// [`ColorGraph::derive_waived`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Waiver {
    /// Waiver identifier (ticket, memo).
    pub id: String,
    /// The human who accepts responsibility.
    pub signer: String,
    /// Why.
    pub reason: String,
}

/// The canonical scope string a color-claim grant must carry.
pub const WAIVER_SCOPE_COLOR_UPGRADE: &str = "color-upgrade";

/// The canonical scope string a SOURCE-color grant must carry (bead
/// gp3.16). Distinct from [`WAIVER_SCOPE_COLOR_UPGRADE`] so a grant
/// authorizing a derived upgrade can never mint a positive leaf, and
/// vice versa.
pub const WAIVER_SCOPE_SOURCE_COLOR: &str = "source-color";

/// TYPED origin evidence for a positive-colored SOURCE leaf (bead
/// gp3.16). Mirrors the package schema-v7 claim-origin vocabulary
/// (fs-package `ClaimOrigin::SourceCertificate` / `AnchoredSource`)
/// without coupling this layer upward: the semantics agree, the types
/// live here. The origin is an INPUT that re-derives the color, not a
/// memo riding alongside it — a Verified leaf is minted through
/// [`fs_evidence::verified_from`] on the carried certificate, and a
/// Validated leaf must name its anchoring dataset exactly. Public data
/// is not authority: [`ColorGraph::source_with_origin`] additionally
/// requires an injected [`SourceOriginVerifier`], whose atomic policy decision
/// is retained in the node provenance.
#[derive(Debug, Clone, PartialEq)]
pub enum SourceOrigin {
    /// A Verified leaf's minting certificate plus the producer identity
    /// (e.g. "fs-solver/ivp-cert"). The color is RE-DERIVED via
    /// [`fs_evidence::verified_from`]; anything weaker than an
    /// exact/enclosure certificate refuses, and the certificate's
    /// interval must match the claimed color bit-exactly.
    Certificate {
        /// Non-blank producer identity.
        producer: String,
        /// Content address of the retained certificate artifact. Two proof
        /// objects yielding the same interval remain distinct, subpoenable
        /// pieces of evidence.
        certificate_hash: ContentHash,
        /// The interval certificate that mints the color.
        certificate: NumericalCertificate,
    },
    /// A Validated leaf's anchoring dataset by identity + content hash.
    /// The id must equal the color's named dataset exactly.
    Anchoring {
        /// The anchoring dataset identity.
        dataset_id: String,
        /// Content hash of the dataset artifact.
        content_hash: ContentHash,
        /// The exact regime attested by that dataset. Carrying it in the
        /// origin lets the gate rederive the complete Validated color
        /// instead of accepting a caller-asserted validity box.
        regime: ValidityDomain,
    },
}

/// Atomic result of an injected admission-policy decision. The decision and
/// the policy identity come from one callback, so a mutable verifier cannot
/// accept under one trust configuration and report another configuration in a
/// second call. Fields are private: capability implementations construct
/// decisions through [`Self::accept`] and [`Self::reject`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PolicyDecision {
    accepted: bool,
    policy_fingerprint: ContentHash,
}

impl PolicyDecision {
    /// Accept under the exact policy identified by `policy_fingerprint`.
    #[must_use]
    pub const fn accept(policy_fingerprint: ContentHash) -> Self {
        Self {
            accepted: true,
            policy_fingerprint,
        }
    }

    /// Reject under the exact policy identified by `policy_fingerprint`.
    #[must_use]
    pub const fn reject(policy_fingerprint: ContentHash) -> Self {
        Self {
            accepted: false,
            policy_fingerprint,
        }
    }

    /// Whether the request was accepted.
    #[must_use]
    pub const fn accepted(self) -> bool {
        self.accepted
    }

    /// Stable identity of the trust roots and decision semantics used.
    #[must_use]
    pub const fn policy_fingerprint(self) -> ContentHash {
        self.policy_fingerprint
    }
}

/// The exact source-admission question presented to a trusted origin
/// verifier. Fields are private/read-only so a verifier cannot observe a
/// partial request: node identity, claimed color, and the complete origin
/// always travel together.
#[derive(Debug, Clone, Copy)]
pub struct SourceOriginRequest<'a> {
    node_name: &'a str,
    claimed_color: &'a Color,
    origin: &'a SourceOrigin,
}

impl<'a> SourceOriginRequest<'a> {
    /// Build the exact request a source gate will present.
    #[must_use]
    pub fn new(node_name: &'a str, claimed_color: &'a Color, origin: &'a SourceOrigin) -> Self {
        Self {
            node_name,
            claimed_color,
            origin,
        }
    }

    /// Node identity covered by this admission.
    #[must_use]
    pub fn node_name(&self) -> &str {
        self.node_name
    }

    /// Exact claimed color covered by this admission.
    #[must_use]
    pub fn claimed_color(&self) -> &Color {
        self.claimed_color
    }

    /// Complete certificate or anchor presented for admission.
    #[must_use]
    pub fn origin(&self) -> &SourceOrigin {
        self.origin
    }

    /// Domain-separated, versioned, length-prefixed identity bytes for
    /// capability implementations that authenticate a request by MAC or
    /// signature. Every floating-point field is bit-exact.
    #[must_use]
    pub fn canonical_bytes(&self) -> Vec<u8> {
        let mut out = vec![1u8];
        push_field(&mut out, SOURCE_ORIGIN_REQUEST_DOMAIN);
        push_field(&mut out, self.node_name.as_bytes());
        push_field(&mut out, &self.claimed_color.canonical_bytes());
        push_source_origin(&mut out, self.origin);
        out
    }
}

/// Capability that authenticates a typed source origin. Merely constructing
/// public certificate fields or writing down a dataset hash is not authority;
/// the injected verifier must resolve and accept the whole request.
pub trait SourceOriginVerifier {
    /// Atomic acceptance and policy identity for this exact request.
    fn verify(&self, request: &SourceOriginRequest<'_>) -> PolicyDecision;
}

/// Fail-closed default when no source-certificate or dataset authority is
/// wired into the admission path.
#[derive(Debug, Default, Clone, Copy)]
pub struct NoSourceOriginVerifier;

impl SourceOriginVerifier for NoSourceOriginVerifier {
    fn verify(&self, _request: &SourceOriginRequest<'_>) -> PolicyDecision {
        PolicyDecision::reject(hash_bytes(
            b"frankensim/fs-ledger/source-origin-policy/v1/deny-all",
        ))
    }
}

/// Why a typed source origin failed to mint the claimed color
/// (structured, teaching — the forged-source refusals).
#[derive(Debug, Clone, PartialEq)]
pub enum SourceOriginRejection {
    /// The origin kind does not fit the color (a certificate cannot
    /// anchor a Validated claim; a dataset cannot certify an interval).
    OriginKindMismatch {
        /// The claimed color's stable name.
        color: &'static str,
    },
    /// [`fs_evidence::verified_from`] refused the certificate
    /// (estimate/no-claim kind, NaN or inverted bounds).
    CertificateRefused {
        /// The evidence-layer refusal, verbatim.
        why: String,
    },
    /// The certificate re-derives a DIFFERENT Verified color than
    /// claimed (bit-exact comparison).
    CertificateMismatch,
    /// The origin names a different dataset than the Validated color.
    DatasetMismatch {
        /// The dataset the origin names.
        origin: String,
        /// The dataset the color names.
        color: String,
    },
    /// The anchoring origin carries a different regime than the claimed
    /// Validated color.
    RegimeMismatch,
    /// Estimated leaves state their own dispersion; they carry no
    /// origin and no waiver (use [`ColorGraph::source`]).
    EstimatedNeedsNoOrigin,
    /// The producer identity is blank, placeholder text, or padded.
    BlankProducer,
    /// The anchoring dataset identity is blank, placeholder text, or padded.
    BlankDataset,
    /// The anchoring regime is empty or contains an unusable axis.
    InvalidRegime {
        /// Empty for an undeclared regime; otherwise the malformed axis.
        axis: String,
    },
    /// The injected source-origin capability did not authenticate the
    /// complete node/color/origin request under the named policy.
    VerifierRefused {
        /// Stable identity of the rejecting policy.
        policy_fingerprint: ContentHash,
    },
    /// The injected capability panicked. External trust code cannot unwind
    /// through the ledger write gate or leave a partially admitted node.
    VerifierPanicked,
}

impl core::fmt::Display for SourceOriginRejection {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::OriginKindMismatch { color } => {
                write!(f, "origin kind cannot mint a {color} source")
            }
            Self::CertificateRefused { why } => write!(f, "certificate refused: {why}"),
            Self::CertificateMismatch => {
                f.write_str("certificate interval differs from the claimed Verified color")
            }
            Self::DatasetMismatch { origin, color } => write!(
                f,
                "anchoring dataset `{origin}` differs from claimed dataset `{color}`"
            ),
            Self::RegimeMismatch => {
                f.write_str("anchoring regime differs from the claimed Validated regime")
            }
            Self::EstimatedNeedsNoOrigin => {
                f.write_str("Estimated sources must use `source` without an origin")
            }
            Self::BlankProducer => {
                f.write_str("certificate producer identity is blank, placeholder text, or padded")
            }
            Self::BlankDataset => {
                f.write_str("anchoring dataset identity is blank, placeholder text, or padded")
            }
            Self::InvalidRegime { axis } if axis.is_empty() => {
                f.write_str("anchoring regime declares no bounded axes")
            }
            Self::InvalidRegime { axis } => {
                write!(f, "anchoring regime axis `{axis}` has invalid bounds")
            }
            Self::VerifierRefused { policy_fingerprint } => write!(
                f,
                "source-origin policy {} refused the complete admission request",
                policy_fingerprint.to_hex()
            ),
            Self::VerifierPanicked => {
                f.write_str("source-origin verifier panicked; admission failed closed")
            }
        }
    }
}

impl std::error::Error for SourceOriginRejection {}

impl SourceOrigin {
    fn derive_color(&self) -> Result<Color, SourceOriginRejection> {
        match self {
            SourceOrigin::Certificate {
                producer,
                certificate,
                ..
            } => {
                if identity_reason(producer).is_some() {
                    return Err(SourceOriginRejection::BlankProducer);
                }
                verified_from(certificate).map_err(|error| {
                    SourceOriginRejection::CertificateRefused {
                        why: error.to_string(),
                    }
                })
            }
            SourceOrigin::Anchoring {
                dataset_id, regime, ..
            } => {
                if color_leaf_identity_reason(dataset_id).is_some() {
                    return Err(SourceOriginRejection::BlankDataset);
                }
                if regime.bounds().is_empty() {
                    return Err(SourceOriginRejection::InvalidRegime {
                        axis: String::new(),
                    });
                }
                if let Some((axis, _)) = regime.bounds().iter().find(|(axis, (lo, hi))| {
                    identity_reason(axis).is_some() || !lo.is_finite() || !hi.is_finite() || lo > hi
                }) {
                    return Err(SourceOriginRejection::InvalidRegime { axis: axis.clone() });
                }
                Ok(Color::Validated {
                    regime: regime.clone(),
                    dataset: dataset_id.clone(),
                })
            }
        }
    }
}

/// Why a color payload is structurally unusable even if a policy authority
/// signs the exact bytes. Waivers may authorize a claim-strength exception;
/// they never authorize malformed epistemic data.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ColorStructureRejection {
    /// A machine identity is blank, reserved placeholder text, or padded.
    InvalidIdentity {
        /// Stable payload field (`dataset`, `axis`, or `estimator`).
        field: &'static str,
        /// The offending identity when it is useful for localization.
        value: String,
        /// Stable reason (`blank`, `placeholder`, or `surrounding whitespace`).
        reason: &'static str,
    },
    /// A Verified interval contains NaN or is inverted. Ordered infinite
    /// endpoints remain sound, possibly vacuous enclosures.
    InvalidVerifiedInterval {
        /// Stable field-level reason.
        reason: &'static str,
    },
    /// A Validated color has no regime axes or one malformed axis.
    InvalidValidatedRegime {
        /// Empty for a wholly undeclared regime; otherwise the offending axis.
        axis: String,
        /// Stable field-level reason.
        reason: &'static str,
    },
    /// An Estimated dispersion is NaN or negative.
    InvalidEstimatedDispersion {
        /// Stable field-level reason.
        reason: &'static str,
    },
}

impl core::fmt::Display for ColorStructureRejection {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidIdentity {
                field,
                value,
                reason,
            } => write!(f, "{field} identity {value:?} is invalid: {reason}"),
            Self::InvalidVerifiedInterval { reason } => {
                write!(f, "Verified interval is invalid: {reason}")
            }
            Self::InvalidValidatedRegime { axis, reason } if axis.is_empty() => {
                write!(f, "Validated regime is invalid: {reason}")
            }
            Self::InvalidValidatedRegime { axis, reason } => {
                write!(f, "Validated regime axis {axis:?} is invalid: {reason}")
            }
            Self::InvalidEstimatedDispersion { reason } => {
                write!(f, "Estimated dispersion is invalid: {reason}")
            }
        }
    }
}

impl std::error::Error for ColorStructureRejection {}

const WAIVER_PAYLOAD_DOMAIN: &[u8] = b"frankensim/fs-ledger/color-waiver";
const COLOR_NODE_HASH_DOMAIN: &[u8] = b"frankensim/fs-ledger/color-node/v2";
const SOURCE_ORIGIN_REQUEST_DOMAIN: &[u8] = b"frankensim/fs-ledger/source-origin-request";

fn interval_op_tag(op: IntervalOp) -> u8 {
    match op {
        IntervalOp::Add => 1,
        IntervalOp::Mul => 2,
        IntervalOp::Hull => 3,
    }
}

fn interval_op_name(op: IntervalOp) -> &'static str {
    match op {
        IntervalOp::Add => "add",
        IntervalOp::Mul => "mul",
        IntervalOp::Hull => "hull",
    }
}

fn numerical_kind_tag(kind: fs_evidence::NumericalKind) -> u8 {
    match kind {
        fs_evidence::NumericalKind::Exact => 1,
        fs_evidence::NumericalKind::Enclosure => 2,
        fs_evidence::NumericalKind::Estimate => 3,
        fs_evidence::NumericalKind::NoClaim => 4,
    }
}

fn numerical_kind_name(kind: fs_evidence::NumericalKind) -> &'static str {
    match kind {
        fs_evidence::NumericalKind::Exact => "exact",
        fs_evidence::NumericalKind::Enclosure => "enclosure",
        fs_evidence::NumericalKind::Estimate => "estimate",
        fs_evidence::NumericalKind::NoClaim => "no-claim",
    }
}

fn human_text_reason(value: &str) -> Option<&'static str> {
    if value.len() > MAX_WAIVER_REASON_BYTES {
        return Some("too-long");
    }
    let trimmed = value.trim();
    if trimmed.is_empty() {
        Some("blank")
    } else if trimmed != value {
        Some("surrounding-whitespace")
    } else if value.chars().any(|ch| {
        ch.is_control()
            || matches!(
                ch,
                '\u{200e}'
                    | '\u{200f}'
                    | '\u{202a}'..='\u{202e}'
                    | '\u{2066}'..='\u{2069}'
            )
    }) {
        Some("control-character")
    } else if is_placeholder_token(value) {
        Some("placeholder")
    } else {
        None
    }
}

fn waiver_annotation_reason(waiver: &Waiver) -> Option<(&'static str, &'static str)> {
    for (field, value) in [
        ("waiver_id", waiver.id.as_str()),
        ("signer", waiver.signer.as_str()),
    ] {
        if let Some(reason) = identity_reason(value) {
            return Some((field, reason));
        }
    }
    human_text_reason(&waiver.reason).map(|reason| ("reason", reason))
}

fn validate_color_structure(color: &Color) -> Result<(), ColorStructureRejection> {
    if let Color::Validated { regime, .. } = color
        && regime.bounds().len() > MAX_VALIDITY_AXES
    {
        return Err(ColorStructureRejection::InvalidValidatedRegime {
            axis: String::new(),
            reason: "validity regime exceeds the axis limit",
        });
    }
    validate_color_payload(color).map_err(|error| match error {
        ColorPayloadError::InvalidIdentity {
            field,
            value,
            reason,
        } => ColorStructureRejection::InvalidIdentity {
            field,
            value,
            reason,
        },
        ColorPayloadError::InvalidVerifiedInterval { reason } => {
            ColorStructureRejection::InvalidVerifiedInterval { reason }
        }
        ColorPayloadError::InvalidValidatedRegime { axis, reason } => {
            ColorStructureRejection::InvalidValidatedRegime { axis, reason }
        }
        ColorPayloadError::InvalidEstimatedDispersion { reason } => {
            ColorStructureRejection::InvalidEstimatedDispersion { reason }
        }
    })
}

fn validate_source_origin_resource_limits(origin: &SourceOrigin) -> Result<(), ColorWriteError> {
    if let SourceOrigin::Anchoring { regime, .. } = origin
        && regime.bounds().len() > MAX_VALIDITY_AXES
    {
        return Err(ColorWriteError::InvalidColor {
            rejection: ColorStructureRejection::InvalidValidatedRegime {
                axis: String::new(),
                reason: "validity regime exceeds the axis limit",
            },
        });
    }
    Ok(())
}

fn estimated_payload_error(color: &Color) -> Option<(&'static str, &'static str)> {
    match validate_color_structure(color) {
        Err(ColorStructureRejection::InvalidIdentity {
            field: "estimator",
            reason,
            ..
        }) => Some(("estimator", reason)),
        Err(ColorStructureRejection::InvalidEstimatedDispersion { reason }) => {
            Some(("dispersion", reason))
        }
        Ok(()) | Err(_) => None,
    }
}

fn estimated_source_payload_error(color: &Color) -> Option<(&'static str, &'static str)> {
    if let Color::Estimated { estimator, .. } = color
        && let Some(why) = color_leaf_identity_reason(estimator)
    {
        return Some(("estimator", why));
    }
    estimated_payload_error(color)
}

fn push_source_origin(out: &mut Vec<u8>, origin: &SourceOrigin) {
    match origin {
        SourceOrigin::Certificate {
            producer,
            certificate_hash,
            certificate,
        } => {
            out.push(1);
            push_field(out, producer.as_bytes());
            out.extend_from_slice(certificate_hash.as_bytes());
            out.push(numerical_kind_tag(certificate.kind));
            out.extend_from_slice(&certificate.lo.to_bits().to_le_bytes());
            out.extend_from_slice(&certificate.hi.to_bits().to_le_bytes());
        }
        SourceOrigin::Anchoring {
            dataset_id,
            content_hash,
            regime,
        } => {
            out.push(2);
            push_field(out, dataset_id.as_bytes());
            out.extend_from_slice(content_hash.as_bytes());
            let color = Color::Validated {
                regime: regime.clone(),
                dataset: dataset_id.clone(),
            };
            push_field(out, &color.canonical_bytes());
        }
    }
}

fn source_origin_canonical_bytes(origin: &SourceOrigin) -> Vec<u8> {
    let mut bytes = Vec::new();
    push_source_origin(&mut bytes, origin);
    bytes
}

/// An AUTHENTICATED waiver: a versioned, length-prefixed payload bound
/// to the exact node identity, evidence lineage, claimed color, scope,
/// signer key, and expiry — plus signature bytes over that payload.
/// Verification happens through a caller-supplied [`WaiverVerifier`]
/// capability; the grant travels whole in the provenance hash.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct WaiverGrant {
    /// The human annotation riding along (never itself authorizing).
    pub annotation: Waiver,
    /// Issuer key identity the verifier resolves.
    pub key_id: String,
    /// Must equal [`WAIVER_SCOPE_COLOR_UPGRADE`] for color upgrades.
    pub scope: String,
    /// The node name this grant is bound to.
    pub node_name: String,
    /// The exact versioned [`Color::canonical_bytes`] being authorized.
    pub claimed_color: Vec<u8>,
    /// The exact parent provenance hashes, in write order — binds the
    /// grant to one evidence lineage (replay to another node fails).
    pub parent_hashes: Vec<ContentHash>,
    /// Last day the grant is valid (days since 2026-01-01).
    pub expires_day: u32,
    /// Signature bytes over [`WaiverGrant::signing_payload`].
    pub signature: Vec<u8>,
}

impl WaiverGrant {
    /// Canonical signing payload, DOMAIN-SEPARATED, VERSIONED, and
    /// LENGTH-PREFIXED (no delimiters, so adversarial text cannot collide
    /// structurally): version byte 3, domain string, operation tag, then each
    /// field as u64-LE length + bytes, parent count + raw 32-byte hashes, and
    /// expiry as u32 LE. Version 3 binds the operation as well as the full
    /// bit-exact color payload, so an Add grant cannot authorize Mul. The
    /// signature is NOT part of its own payload.
    #[must_use]
    pub fn signing_payload(&self, op: IntervalOp) -> Vec<u8> {
        let mut out = vec![3u8];
        push_field(&mut out, WAIVER_PAYLOAD_DOMAIN);
        out.push(interval_op_tag(op));
        for field in [
            self.key_id.as_str(),
            self.scope.as_str(),
            self.node_name.as_str(),
        ] {
            push_field(&mut out, field.as_bytes());
        }
        push_field(&mut out, &self.claimed_color);
        for field in [
            self.annotation.id.as_str(),
            self.annotation.signer.as_str(),
            self.annotation.reason.as_str(),
        ] {
            push_field(&mut out, field.as_bytes());
        }
        push_len(&mut out, self.parent_hashes.len());
        for h in &self.parent_hashes {
            out.extend_from_slice(h.as_bytes());
        }
        out.extend_from_slice(&self.expires_day.to_le_bytes());
        out
    }

    /// Canonical signing payload for a SOURCE-color grant (bead gp3.16):
    /// version byte 4, operation tag 0 (a leaf has no composition
    /// operation), otherwise field-for-field identical to
    /// [`WaiverGrant::signing_payload`]. A v3 derive payload can never
    /// collide with a v4 source payload (distinct version bytes), so a
    /// signature over one cannot authorize the other.
    #[must_use]
    pub fn signing_payload_source(&self) -> Vec<u8> {
        let mut out = vec![4u8];
        push_field(&mut out, WAIVER_PAYLOAD_DOMAIN);
        out.push(0); // no operation: source leaf
        for field in [
            self.key_id.as_str(),
            self.scope.as_str(),
            self.node_name.as_str(),
        ] {
            push_field(&mut out, field.as_bytes());
        }
        push_field(&mut out, &self.claimed_color);
        for field in [
            self.annotation.id.as_str(),
            self.annotation.signer.as_str(),
            self.annotation.reason.as_str(),
        ] {
            push_field(&mut out, field.as_bytes());
        }
        push_len(&mut out, self.parent_hashes.len());
        for h in &self.parent_hashes {
            out.extend_from_slice(h.as_bytes());
        }
        out.extend_from_slice(&self.expires_day.to_le_bytes());
        out
    }

    fn signing_payload_for(&self, operation: Option<IntervalOp>) -> Vec<u8> {
        operation.map_or_else(
            || self.signing_payload_source(),
            |op| self.signing_payload(op),
        )
    }

    fn payload_version(operation: Option<IntervalOp>) -> u8 {
        if operation.is_some() { 3 } else { 4 }
    }
}

fn push_len(out: &mut Vec<u8>, len: usize) {
    let len = u64::try_from(len).expect("a Rust allocation length always fits in u64");
    out.extend_from_slice(&len.to_le_bytes());
}

fn push_field(out: &mut Vec<u8>, bytes: &[u8]) {
    push_len(out, bytes.len());
    out.extend_from_slice(bytes);
}

/// The signature-verification CAPABILITY (injected; this crate ships
/// no cryptography). Implementations resolve `key_id` and check
/// `signature` over `payload`.
pub trait WaiverVerifier {
    /// Atomic authentication decision and policy identity for the exact
    /// `key_id`, `payload`, and `signature` tuple.
    fn verify(&self, key_id: &str, payload: &[u8], signature: &[u8]) -> PolicyDecision;
}

/// The in-tree default: NO verifier exists, so NOTHING authenticates
/// (the no-crypto no-claim — fail closed until a Franken-compliant
/// signature capability is wired in).
#[derive(Debug, Default, Clone, Copy)]
pub struct NoWaiverVerifier;

impl WaiverVerifier for NoWaiverVerifier {
    fn verify(&self, _key_id: &str, _payload: &[u8], _signature: &[u8]) -> PolicyDecision {
        PolicyDecision::reject(hash_bytes(
            b"frankensim/fs-ledger/waiver-policy/v1/deny-all",
        ))
    }
}

/// Why a grant failed to authorize (structured, teaching).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum WaiverRejection {
    /// Authenticated metadata is malformed. A verifier cannot legitimize a
    /// blank, placeholder, padded, hostile, or oversized authority field.
    InvalidField {
        /// Stable field name.
        field: &'static str,
        /// Stable structural reason.
        reason: &'static str,
    },
    /// A bounded authority collection exceeded its declared limit.
    ResourceLimitExceeded {
        /// Stable resource name.
        resource: &'static str,
        /// Maximum accepted entries or bytes.
        limit: usize,
        /// Offered entries or bytes.
        actual: usize,
    },
    /// Scope is not [`WAIVER_SCOPE_COLOR_UPGRADE`].
    ScopeMismatch,
    /// The grant names a different node.
    NodeMismatch,
    /// The grant authorizes a different color than claimed.
    ColorMismatch,
    /// The grant's parent hashes differ from the actual lineage
    /// (replay to another node / tampered evidence).
    LineageMismatch,
    /// Expired as of the supplied day.
    Expired,
    /// The verifier refused the signature (wrong key, tampered payload,
    /// rotated-out key, or no verifier capability at all) under the named
    /// policy.
    VerifierRefused {
        /// Stable identity of the rejecting signature policy.
        policy_fingerprint: ContentHash,
    },
    /// The injected signature capability panicked; admission failed closed.
    VerifierPanicked,
}

impl core::fmt::Display for WaiverRejection {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidField { field, reason } => {
                write!(f, "grant field `{field}` is invalid ({reason})")
            }
            Self::ResourceLimitExceeded {
                resource,
                limit,
                actual,
            } => write!(
                f,
                "grant resource `{resource}` exceeds limit {limit} (offered {actual})"
            ),
            Self::ScopeMismatch => f.write_str("scope does not authorize this write kind"),
            Self::NodeMismatch => f.write_str("grant names a different node"),
            Self::ColorMismatch => f.write_str("grant authorizes a different color"),
            Self::LineageMismatch => {
                f.write_str("grant parent hashes differ from the actual lineage")
            }
            Self::Expired => f.write_str("grant was expired at the admission date"),
            Self::VerifierRefused { policy_fingerprint } => write!(
                f,
                "signature policy {} refused the grant",
                policy_fingerprint.to_hex()
            ),
            Self::VerifierPanicked => {
                f.write_str("signature verifier panicked; admission failed closed")
            }
        }
    }
}

impl std::error::Error for WaiverRejection {}

fn validate_waiver_grant(grant: &WaiverGrant) -> Result<(), WaiverRejection> {
    if let Some((field, reason)) = waiver_annotation_reason(&grant.annotation) {
        return Err(WaiverRejection::InvalidField { field, reason });
    }
    for (field, value) in [
        ("key_id", grant.key_id.as_str()),
        ("scope", grant.scope.as_str()),
        ("node_name", grant.node_name.as_str()),
    ] {
        if let Some(reason) = identity_reason(value) {
            return Err(WaiverRejection::InvalidField { field, reason });
        }
    }
    if grant.claimed_color.is_empty() {
        return Err(WaiverRejection::InvalidField {
            field: "claimed_color",
            reason: "blank",
        });
    }
    if grant.claimed_color.len() > MAX_CLAIMED_COLOR_BYTES {
        return Err(WaiverRejection::ResourceLimitExceeded {
            resource: "claimed_color_bytes",
            limit: MAX_CLAIMED_COLOR_BYTES,
            actual: grant.claimed_color.len(),
        });
    }
    if grant.parent_hashes.len() > MAX_COLOR_PARENTS {
        return Err(WaiverRejection::ResourceLimitExceeded {
            resource: "parent_hashes",
            limit: MAX_COLOR_PARENTS,
            actual: grant.parent_hashes.len(),
        });
    }
    if grant.signature.is_empty() {
        return Err(WaiverRejection::InvalidField {
            field: "signature",
            reason: "blank",
        });
    }
    if grant.signature.len() > MAX_WAIVER_SIGNATURE_BYTES {
        return Err(WaiverRejection::ResourceLimitExceeded {
            resource: "signature_bytes",
            limit: MAX_WAIVER_SIGNATURE_BYTES,
            actual: grant.signature.len(),
        });
    }
    Ok(())
}

fn json_f64(value: f64) -> String {
    if value.is_finite() {
        value.to_string()
    } else {
        format!("\"non-finite:{value}\"")
    }
}

fn json_string(value: &str) -> String {
    use core::fmt::Write as _;
    let mut out = String::with_capacity(value.len() + 2);
    out.push('"');
    for ch in value.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if u32::from(c) < 0x20 => {
                let _ = write!(out, "\\u{:04x}", u32::from(c));
            }
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

fn hex_bytes(bytes: &[u8]) -> String {
    use core::fmt::Write as _;
    let mut out = String::with_capacity(bytes.len() * 2);
    for byte in bytes {
        let _ = write!(out, "{byte:02x}");
    }
    out
}

fn optional_hex_bytes_json(bytes: Option<&[u8]>) -> String {
    bytes.map_or("null".to_string(), |bytes| {
        format!("\"{}\"", hex_bytes(bytes))
    })
}

fn parent_hashes_json(parent_hashes: &[ContentHash]) -> String {
    parent_hashes
        .iter()
        .map(|hash| format!("\"{}\"", hash.to_hex()))
        .collect::<Vec<_>>()
        .join(",")
}

fn waiver_json(waiver: Option<&Waiver>) -> String {
    waiver.map_or("null".to_string(), |waiver| {
        format!(
            "{{\"id\":{},\"signer\":{},\"reason\":{}}}",
            json_string(&waiver.id),
            json_string(&waiver.signer),
            json_string(&waiver.reason)
        )
    })
}

fn grant_json(grant: Option<&WaiverGrant>, operation: Option<IntervalOp>) -> String {
    grant.map_or("null".to_string(), |grant| {
        let signing_payload = grant.signing_payload_for(operation);
        let payload_version = WaiverGrant::payload_version(operation);
        format!(
            "{{\"payload_version\":{payload_version},\"key_id\":{},\"scope\":{},\"node_name\":{},\
             \"claimed_color_hex\":\"{}\",\"parent_hashes\":[{}],\"expires_day\":{},\
             \"signing_payload_hex\":\"{}\",\"signature_hex\":\"{}\",\
             \"authorized\":true}}",
            json_string(&grant.key_id),
            json_string(&grant.scope),
            json_string(&grant.node_name),
            hex_bytes(&grant.claimed_color),
            parent_hashes_json(&grant.parent_hashes),
            grant.expires_day,
            hex_bytes(&signing_payload),
            hex_bytes(&grant.signature),
        )
    })
}

fn waiver_dependencies_json(dependencies: &[WaiverDependency]) -> String {
    dependencies
        .iter()
        .map(|dependency| {
            let operation = dependency
                .operation
                .map_or("null".to_string(), |op| json_string(interval_op_name(op)));
            format!(
                "{{\"authorizing_node\":{},\"operation\":{},\"policy_fingerprint\":\"{}\",\
                 \"admission_day\":{},\"waiver\":{},\"grant\":{}}}",
                dependency.authorizing_node,
                operation,
                dependency.policy_fingerprint.to_hex(),
                dependency.admission_day,
                waiver_json(Some(&dependency.grant.annotation)),
                grant_json(Some(&dependency.grant), dependency.operation),
            )
        })
        .collect::<Vec<_>>()
        .join(",")
}

fn origin_json(origin: Option<&SourceOrigin>) -> String {
    origin.map_or("null".to_string(), |origin| match origin {
        SourceOrigin::Certificate {
            producer,
            certificate_hash,
            certificate,
        } => format!(
            "{{\"kind\":\"certificate\",\"producer\":{},\"certificate_kind\":\
             {},\"certificate_hash\":\"{}\",\"lo\":{},\"hi\":{}}}",
            json_string(producer),
            json_string(numerical_kind_name(certificate.kind)),
            certificate_hash.to_hex(),
            json_f64(certificate.lo),
            json_f64(certificate.hi)
        ),
        SourceOrigin::Anchoring {
            dataset_id,
            content_hash,
            regime,
        } => format!(
            "{{\"kind\":\"anchoring\",\"dataset\":{},\"content_hash\":\"{}\",\
             \"regime\":{}}}",
            json_string(dataset_id),
            content_hash.to_hex(),
            Color::Validated {
                regime: regime.clone(),
                dataset: dataset_id.clone(),
            }
            .payload_json()
        ),
    })
}

fn optional_hash_json(hash: Option<ContentHash>) -> String {
    hash.map_or("null".to_string(), |hash| format!("\"{}\"", hash.to_hex()))
}

fn optional_u32_json(value: Option<u32>) -> String {
    value.map_or("null".to_string(), |value| value.to_string())
}

/// One regime-exit demotion observed while folding a derived node.
/// The parent POSITION is part of the record because a legal parent
/// list may contain the same node more than once; an id alone would
/// make replay ambiguous. Entries are stored in ascending position.
#[derive(Debug, Clone, PartialEq)]
pub struct ColorDemotion {
    parent_index: usize,
    parent_id: u64,
    reason: Demotion,
}

impl ColorDemotion {
    /// Position in the derived node's parent list.
    #[must_use]
    pub fn parent_index(&self) -> usize {
        self.parent_index
    }

    /// Id found at [`Self::parent_index`].
    #[must_use]
    pub fn parent_id(&self) -> u64 {
        self.parent_id
    }

    /// The regime-exit diagnosis.
    #[must_use]
    pub fn reason(&self) -> &Demotion {
        &self.reason
    }
}

/// One authenticated waiver on which this node depends transitively. The
/// authorizing node id and its original operation make the full signed grant
/// independently resolvable during replay; entries are canonicalized by
/// ascending authorizing-node id and never duplicated.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct WaiverDependency {
    authorizing_node: u64,
    operation: Option<IntervalOp>,
    grant: WaiverGrant,
    policy_fingerprint: ContentHash,
    admission_day: u32,
}

impl WaiverDependency {
    /// Node at which the grant was originally authenticated.
    #[must_use]
    pub fn authorizing_node(&self) -> u64 {
        self.authorizing_node
    }

    /// Operation covered by the original grant (`None` for a waived source).
    #[must_use]
    pub fn operation(&self) -> Option<IntervalOp> {
        self.operation
    }

    /// Complete signed grant retained from the authorizing node.
    #[must_use]
    pub fn grant(&self) -> &WaiverGrant {
        &self.grant
    }

    /// Policy that authenticated the grant at its authorizing node.
    #[must_use]
    pub fn policy_fingerprint(&self) -> ContentHash {
        self.policy_fingerprint
    }

    /// Historical day on which the grant was admitted.
    #[must_use]
    pub fn admission_day(&self) -> u32 {
        self.admission_day
    }
}

/// One colored ledger node. Fields are PRIVATE and read-only (bead
/// gp3.16): a written node cannot be edited after the gate accepted
/// it — the only mutation surface on the graph is the gated write
/// methods, so provenance hashes always describe what they cover.
#[derive(Debug, Clone)]
pub struct ColorNode {
    id: u64,
    name: String,
    color: Color,
    parents: Vec<u64>,
    operation: Option<IntervalOp>,
    /// EVERY regime demotion that fired while folding the parents, as
    /// canonical order (ascending parent position in the write's
    /// parent list).
    demotions: Vec<ColorDemotion>,
    origin: Option<SourceOrigin>,
    origin_policy_fingerprint: Option<ContentHash>,
    waiver: Option<Waiver>,
    grant: Option<WaiverGrant>,
    waiver_policy_fingerprint: Option<ContentHash>,
    waiver_admission_day: Option<u32>,
    waiver_dependencies: Vec<WaiverDependency>,
    hash: ContentHash,
}

impl ColorNode {
    /// Node id (write order).
    #[must_use]
    pub fn id(&self) -> u64 {
        self.id
    }

    /// Human name.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// The color declaration as written, without applying waiver taint.
    /// Scientific admission code must use [`Self::scientific_color`]; this
    /// accessor is intentionally named as unverified so a transitive waiver
    /// cannot disappear behind an innocuous `color()` call.
    #[must_use]
    pub fn declared_color_unverified(&self) -> &Color {
        &self.color
    }

    /// The color available for ordinary scientific admission. Any direct or
    /// inherited waiver makes the declaration unavailable; callers that elect
    /// to consume waived evidence must do so through explicit waiver policy.
    #[must_use]
    pub fn scientific_color(&self) -> Option<&Color> {
        (!self.depends_on_waiver()).then_some(&self.color)
    }

    /// Parent node ids.
    #[must_use]
    pub fn parents(&self) -> &[u64] {
        &self.parents
    }

    /// Composition operation (`None` only for source nodes).
    #[must_use]
    pub fn operation(&self) -> Option<IntervalOp> {
        self.operation
    }

    /// Every demotion that fired at this write, in canonical parent-list
    /// order.
    #[must_use]
    pub fn demotions(&self) -> &[ColorDemotion] {
        &self.demotions
    }

    /// The typed source origin, when this is a positive-colored leaf.
    #[must_use]
    pub fn origin(&self) -> Option<&SourceOrigin> {
        self.origin.as_ref()
    }

    /// Policy that authenticated this node's typed source origin.
    #[must_use]
    pub fn origin_policy_fingerprint(&self) -> Option<ContentHash> {
        self.origin_policy_fingerprint
    }

    /// The human annotation, when one was recorded (never authorizing).
    #[must_use]
    pub fn waiver(&self) -> Option<&Waiver> {
        self.waiver.as_ref()
    }

    /// The authenticated grant, when one authorized this write.
    #[must_use]
    pub fn grant(&self) -> Option<&WaiverGrant> {
        self.grant.as_ref()
    }

    /// Policy that authenticated this node's direct waiver grant.
    #[must_use]
    pub fn waiver_policy_fingerprint(&self) -> Option<ContentHash> {
        self.waiver_policy_fingerprint
    }

    /// Historical day on which this node's direct waiver was admitted.
    #[must_use]
    pub fn waiver_admission_day(&self) -> Option<u32> {
        self.waiver_admission_day
    }

    /// Every earlier waiver on which this node transitively depends, sorted by
    /// authorizing node id. The grant that authorized THIS node, if any,
    /// remains available separately through [`Self::grant`].
    #[must_use]
    pub fn waiver_dependencies(&self) -> &[WaiverDependency] {
        &self.waiver_dependencies
    }

    /// Whether this node's scientific claim depends on any authenticated
    /// waiver, either at this write or transitively through its parents.
    #[must_use]
    pub fn depends_on_waiver(&self) -> bool {
        self.grant.is_some() || !self.waiver_dependencies.is_empty()
    }

    /// Provenance hash (name, color bytes, parent hashes, origin and its
    /// policy, waiver and its admission context, and transitive waiver
    /// dependencies).
    #[must_use]
    pub fn hash(&self) -> ContentHash {
        self.hash
    }
}

/// Teaching errors at the write gate.
#[derive(Debug, Clone, PartialEq)]
pub enum ColorWriteError {
    /// The durable node identity is blank, padded, a placeholder, contains
    /// disallowed characters, or exceeds [`MAX_COLOR_NODE_NAME_BYTES`].
    InvalidNodeName {
        /// Stable shared-grammar refusal reason.
        reason: &'static str,
    },
    /// A bounded graph resource exceeded its declared limit before append.
    ResourceLimitExceeded {
        /// Stable resource name.
        resource: &'static str,
        /// Maximum accepted units (entries or bytes, as named by `resource`).
        limit: usize,
        /// Offered units.
        actual: usize,
    },
    /// The claimed color outranks what the parents support.
    LaunderingRefused {
        /// The claimed rank.
        claimed: ColorRank,
        /// The rank the composition algebra derived.
        derived: ColorRank,
        /// The parents that cap the rank.
        offending_parents: Vec<u64>,
    },
    /// A non-waived claim differs from the exact color algebra result.
    ClaimMismatch {
        /// The color the caller attempted to write.
        claimed: Color,
        /// The exact color derived from the parents and operation.
        derived: Color,
    },
    /// A referenced parent does not exist.
    UnknownParent {
        /// The offending id.
        id: u64,
    },
    /// A sealed parent carries a grant without the policy/day context needed
    /// to preserve its authority transitively. This cannot arise through the
    /// public write API, but is kept fallible for future persisted imports.
    InvalidParentAuthority {
        /// Parent whose direct grant context is incomplete.
        id: u64,
    },
    /// Derivations need at least one parent.
    NoParents,
    /// A non-authorizing human waiver annotation is malformed or exceeds its
    /// audit-metadata bounds.
    InvalidWaiverAnnotation {
        /// Stable annotation field.
        field: &'static str,
        /// Shared structural refusal reason.
        reason: &'static str,
    },
    /// A waiver grant failed authentication or binding checks; the
    /// promotion is refused (fail closed).
    WaiverRefused {
        /// The structured reason.
        rejection: WaiverRejection,
    },
    /// A positive-colored LEAF (Validated or Verified) was written
    /// without typed origin evidence or an authenticated grant — the
    /// source-laundering refusal (bead gp3.16).
    SourceOriginRequired {
        /// The rank the leaf claimed.
        rank: ColorRank,
    },
    /// The typed origin evidence failed to mint the claimed color.
    SourceOriginRefused {
        /// The structured reason.
        rejection: SourceOriginRejection,
    },
    /// A direct Estimated leaf has a malformed identity or dispersion.
    InvalidEstimatedSource {
        /// `"estimator"` or `"dispersion"`.
        field: &'static str,
        /// Stable field-level refusal.
        why: &'static str,
    },
    /// A color payload is structurally malformed. Authentication cannot waive
    /// finite/order, regime, identity, or dispersion invariants.
    InvalidColor {
        /// The exact structural refusal.
        rejection: ColorStructureRejection,
    },
}

impl core::fmt::Display for ColorWriteError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            ColorWriteError::InvalidNodeName { reason } => write!(
                f,
                "color graph node name is invalid ({reason}); use a canonical non-placeholder \
                 machine identity of at most {MAX_COLOR_NODE_NAME_BYTES} bytes"
            ),
            ColorWriteError::ResourceLimitExceeded {
                resource,
                limit,
                actual,
            } => write!(
                f,
                "color graph resource `{resource}` exceeds limit {limit} (offered {actual}); \
                 split the derivation or reduce the named retained resource before admission"
            ),
            ColorWriteError::LaunderingRefused {
                claimed,
                derived,
                offending_parents,
            } => write!(
                f,
                "laundering refused: the write claims {claimed:?} but the parents \
                 support at most {derived:?} (capped by nodes {offending_parents:?}); \
                 estimates cannot become certificates by assertion — an authenticated \
                 WaiverGrant via derive_waived is the only path past this refusal, and \
                 it travels whole in provenance"
            ),
            ColorWriteError::ClaimMismatch { claimed, derived } => write!(
                f,
                "color claim mismatch: the write claims {} with payload {} but the exact \
                 parent composition derives {} with payload {}; rank alone is insufficient \
                 because narrowing an interval, widening a regime, or shrinking dispersion \
                 can strengthen a claim — omit the claim to write the derived color, or use \
                 an authenticated WaiverGrant",
                claimed.name(),
                claimed.payload_json(),
                derived.name(),
                derived.payload_json(),
            ),
            ColorWriteError::UnknownParent { id } => {
                write!(f, "parent node {id} does not exist in this color graph")
            }
            ColorWriteError::InvalidParentAuthority { id } => write!(
                f,
                "parent node {id} has a waiver grant without its admitting policy and day"
            ),
            ColorWriteError::NoParents => {
                write!(f, "derived nodes need parents; use `source` for leaves")
            }
            ColorWriteError::InvalidWaiverAnnotation { field, reason } => write!(
                f,
                "waiver annotation field `{field}` is invalid ({reason}); annotations authorize \
                 nothing but still must be bounded, canonical, and safe to display in an audit"
            ),
            ColorWriteError::WaiverRefused { rejection } => write!(
                f,
                "waiver refused ({rejection}): promotion requires an authenticated \
                 grant bound to this node, lineage, color, and scope, unexpired, with \
                 a signature the verifier capability accepts — fail closed otherwise"
            ),
            ColorWriteError::SourceOriginRequired { rank } => write!(
                f,
                "source origin required: a {rank:?} leaf cannot state its color by \
                 assertion — carry typed origin evidence authenticated by a source \
                 verifier (a minting certificate for verified, the anchoring dataset \
                 for validated) via \
                 `source_with_origin`, or an authenticated source-color grant via \
                 `source_waived`; estimates need neither"
            ),
            ColorWriteError::SourceOriginRefused { rejection } => write!(
                f,
                "source origin refused ({rejection}): the typed evidence must \
                 actually mint the claimed color and its verifier must authenticate \
                 the complete request — a certificate re-derives the \
                 verified interval bit-exactly, an anchoring names the validated \
                 dataset exactly — forged or mismatched origins fail closed"
            ),
            ColorWriteError::InvalidEstimatedSource { field, why } => write!(
                f,
                "Estimated source field `{field}` refused: {why}; use a meaningful estimator \
                 identity and a nonnegative dispersion (positive infinity is the explicit \
                 no-spread-claim sentinel)"
            ),
            ColorWriteError::InvalidColor { rejection } => write!(
                f,
                "color payload refused ({rejection}): a waiver may authorize claim policy, but \
                 cannot authorize structurally invalid epistemic data"
            ),
        }
    }
}

impl std::error::Error for ColorWriteError {}

fn waiver_grant_retained_bytes(grant: &WaiverGrant) -> Option<usize> {
    // Fixed framing, hashes, integer fields, and JSON punctuation. Variable
    // fields are counted once here; row rendering is a bounded constant-factor
    // expansion (hex plus the canonical signing payload).
    let mut bytes = 256usize;
    for len in [
        grant.annotation.id.len(),
        grant.annotation.signer.len(),
        grant.annotation.reason.len(),
        grant.key_id.len(),
        grant.scope.len(),
        grant.node_name.len(),
        grant.claimed_color.len(),
        grant.signature.len(),
    ] {
        bytes = bytes.checked_add(len)?;
    }
    bytes.checked_add(grant.parent_hashes.len().checked_mul(32)?)
}

fn add_waiver_closure_bytes(total: &mut usize, grant: &WaiverGrant) -> Result<(), ColorWriteError> {
    let actual = waiver_grant_retained_bytes(grant)
        .and_then(|bytes| total.checked_add(bytes))
        .unwrap_or(usize::MAX);
    if actual > MAX_WAIVER_CLOSURE_BYTES {
        return Err(ColorWriteError::ResourceLimitExceeded {
            resource: "waiver_closure_bytes",
            limit: MAX_WAIVER_CLOSURE_BYTES,
            actual,
        });
    }
    *total = actual;
    Ok(())
}

fn validate_waiver_closure_bytes(
    dependencies: &[WaiverDependency],
    direct_grant: Option<&WaiverGrant>,
) -> Result<(), ColorWriteError> {
    let mut total = 0usize;
    for dependency in dependencies {
        add_waiver_closure_bytes(&mut total, &dependency.grant)?;
    }
    if let Some(grant) = direct_grant {
        add_waiver_closure_bytes(&mut total, grant)?;
    }
    Ok(())
}

fn merge_fold_validity_axes(
    aggregate: &mut BTreeMap<String, (f64, f64)>,
    regime: &ValidityDomain,
) -> Result<bool, ColorWriteError> {
    // Detect an empty intersection before enforcing the union-width budget.
    // Otherwise lexical key order could falsely reject on a new axis that
    // happens to precede a shared, disjoint axis; a disjoint Validated pair
    // honestly becomes Estimated and retains no aggregate regime.
    for (axis, &(lo, hi)) in regime.bounds() {
        if let Some(&(aggregate_lo, aggregate_hi)) = aggregate.get(axis) {
            let intersection_lo = aggregate_lo.max(lo);
            let intersection_hi = aggregate_hi.min(hi);
            if !intersection_lo.is_finite()
                || !intersection_hi.is_finite()
                || intersection_lo > intersection_hi
            {
                return Ok(false);
            }
        }
    }
    let missing = regime
        .bounds()
        .keys()
        .filter(|axis| !aggregate.contains_key(axis.as_str()))
        .count();
    let actual = aggregate.len().saturating_add(missing);
    if actual > MAX_VALIDITY_AXES {
        return Err(ColorWriteError::ResourceLimitExceeded {
            resource: "derived_validity_axes",
            limit: MAX_VALIDITY_AXES,
            actual,
        });
    }
    for (axis, &(lo, hi)) in regime.bounds() {
        aggregate
            .entry(axis.clone())
            .and_modify(|(aggregate_lo, aggregate_hi)| {
                *aggregate_lo = (*aggregate_lo).max(lo);
                *aggregate_hi = (*aggregate_hi).min(hi);
            })
            .or_insert((lo, hi));
    }
    Ok(true)
}

fn preflight_fold_color(
    aggregate: &mut BTreeMap<String, (f64, f64)>,
    estimated_absorbed: &mut bool,
    color: &Color,
) -> Result<(), ColorWriteError> {
    if *estimated_absorbed {
        return Ok(());
    }
    match color {
        Color::Verified { .. } => {}
        Color::Validated { regime, .. } => {
            if !merge_fold_validity_axes(aggregate, regime)? {
                *estimated_absorbed = true;
                aggregate.clear();
            }
        }
        Color::Estimated { .. } => {
            *estimated_absorbed = true;
            aggregate.clear();
        }
    }
    Ok(())
}

struct NodeWriteMetadata {
    operation: Option<IntervalOp>,
    demotions: Vec<ColorDemotion>,
    origin: Option<SourceOrigin>,
    origin_policy_fingerprint: Option<ContentHash>,
    waiver: Option<Waiver>,
    grant: Option<WaiverGrant>,
    waiver_policy_fingerprint: Option<ContentHash>,
    waiver_admission_day: Option<u32>,
    waiver_dependencies: Vec<WaiverDependency>,
}

struct NodeHashMetadata<'a> {
    operation: Option<IntervalOp>,
    demotions: &'a [ColorDemotion],
    origin: Option<&'a SourceOrigin>,
    origin_policy_fingerprint: Option<ContentHash>,
    waiver: Option<&'a Waiver>,
    grant: Option<&'a WaiverGrant>,
    waiver_policy_fingerprint: Option<ContentHash>,
    waiver_admission_day: Option<u32>,
    waiver_dependencies: &'a [WaiverDependency],
}

impl<'a> From<&'a NodeWriteMetadata> for NodeHashMetadata<'a> {
    fn from(metadata: &'a NodeWriteMetadata) -> Self {
        Self {
            operation: metadata.operation,
            demotions: &metadata.demotions,
            origin: metadata.origin.as_ref(),
            origin_policy_fingerprint: metadata.origin_policy_fingerprint,
            waiver: metadata.waiver.as_ref(),
            grant: metadata.grant.as_ref(),
            waiver_policy_fingerprint: metadata.waiver_policy_fingerprint,
            waiver_admission_day: metadata.waiver_admission_day,
            waiver_dependencies: &metadata.waiver_dependencies,
        }
    }
}

/// The write-time color gatekeeper (append-only, deterministic).
#[derive(Debug, Default)]
pub struct ColorGraph {
    nodes: Vec<ColorNode>,
    rows: Vec<String>,
}

impl ColorGraph {
    /// Empty graph.
    #[must_use]
    pub fn new() -> Self {
        ColorGraph::default()
    }

    /// The nodes written so far.
    #[must_use]
    pub fn nodes(&self) -> &[ColorNode] {
        &self.nodes
    }

    /// The canonical JSON rows (one per write, plus demotion events).
    #[must_use]
    pub fn rows(&self) -> &[String] {
        &self.rows
    }

    fn inherited_waiver_dependencies(
        &self,
        parents: &[u64],
        direct_grant: Option<&WaiverGrant>,
    ) -> Result<Vec<WaiverDependency>, ColorWriteError> {
        if parents.len() > MAX_COLOR_PARENTS {
            return Err(ColorWriteError::ResourceLimitExceeded {
                resource: "parents",
                limit: MAX_COLOR_PARENTS,
                actual: parents.len(),
            });
        }
        let mut dependencies = BTreeMap::<u64, WaiverDependency>::new();
        let mut retained_bytes = 0usize;
        if let Some(grant) = direct_grant {
            add_waiver_closure_bytes(&mut retained_bytes, grant)?;
        }
        for parent_id in parents {
            let parent = self
                .node(*parent_id)
                .ok_or(ColorWriteError::UnknownParent { id: *parent_id })?;
            for dependency in &parent.waiver_dependencies {
                if !dependencies.contains_key(&dependency.authorizing_node) {
                    if dependencies.len() == MAX_WAIVER_DEPENDENCIES {
                        return Err(ColorWriteError::ResourceLimitExceeded {
                            resource: "waiver_dependencies",
                            limit: MAX_WAIVER_DEPENDENCIES,
                            actual: dependencies.len() + 1,
                        });
                    }
                    add_waiver_closure_bytes(&mut retained_bytes, &dependency.grant)?;
                    dependencies.insert(dependency.authorizing_node, dependency.clone());
                }
            }
            if let Some(grant) = &parent.grant {
                let (Some(policy_fingerprint), Some(admission_day)) = (
                    parent.waiver_policy_fingerprint,
                    parent.waiver_admission_day,
                ) else {
                    return Err(ColorWriteError::InvalidParentAuthority { id: parent.id });
                };
                if !dependencies.contains_key(&parent.id) {
                    if dependencies.len() == MAX_WAIVER_DEPENDENCIES {
                        return Err(ColorWriteError::ResourceLimitExceeded {
                            resource: "waiver_dependencies",
                            limit: MAX_WAIVER_DEPENDENCIES,
                            actual: dependencies.len() + 1,
                        });
                    }
                    add_waiver_closure_bytes(&mut retained_bytes, grant)?;
                    dependencies.insert(
                        parent.id,
                        WaiverDependency {
                            authorizing_node: parent.id,
                            operation: parent.operation,
                            grant: grant.clone(),
                            policy_fingerprint,
                            admission_day,
                        },
                    );
                }
            }
        }
        Ok(dependencies.into_values().collect())
    }

    /// Provenance hash over DOMAIN-SEPARATED, VERSIONED v9,
    /// LENGTH-PREFIXED encoding. V9 binds color-algebra v2 in both the hash
    /// domain and [`Color::canonical_bytes`]. V8 binds certificate artifact identity,
    /// direct source/waiver policy fingerprints, waiver admission days, and
    /// those fields in the canonical transitive waiver dependency closure.
    /// V7 first bound that dependency closure. V6 bound every regime demotion and the
    /// correct source/derived waiver payload. V5 bound the typed SOURCE ORIGIN (bead
    /// gp3.16) so a forged or substituted origin changes the node
    /// identity and every downstream hash. V4 bound source/derived
    /// status and the exact [`IntervalOp`]; v3 added
    /// [`Color::canonical_bytes`]; the former v2 representation used
    /// rounded display JSON. Length-prefixing prevents adversarial text from
    /// colliding structurally. Color-write row schema v7 persists the exact
    /// color and origin bytes consumed here, so the hash input is reconstructible
    /// without treating display JSON as canonical.
    fn node_hash(
        &self,
        name: &str,
        color: &Color,
        parents: &[u64],
        metadata: &NodeHashMetadata<'_>,
    ) -> ContentHash {
        let color_bytes = color.canonical_bytes();
        let origin_bytes = metadata.origin.map(source_origin_canonical_bytes);
        self.node_hash_from_canonical_payloads(
            name,
            &color_bytes,
            parents,
            metadata,
            origin_bytes.as_deref(),
        )
    }

    fn node_hash_from_canonical_payloads(
        &self,
        name: &str,
        color_bytes: &[u8],
        parents: &[u64],
        metadata: &NodeHashMetadata<'_>,
        origin_bytes: Option<&[u8]>,
    ) -> ContentHash {
        let mut buf = vec![COLOR_NODE_HASH_ENCODING_VERSION];
        push_field(&mut buf, COLOR_NODE_HASH_DOMAIN);
        match metadata.operation {
            Some(op) => {
                buf.push(1);
                buf.push(interval_op_tag(op));
            }
            None => buf.push(0),
        }
        push_field(&mut buf, name.as_bytes());
        push_field(&mut buf, color_bytes);
        push_len(&mut buf, parents.len());
        for &p in parents {
            let parent = self
                .node(p)
                .expect("node_hash parents are validated before append");
            push_field(&mut buf, parent.hash.as_bytes());
        }
        push_len(&mut buf, metadata.demotions.len());
        for demotion in metadata.demotions {
            push_len(&mut buf, demotion.parent_index);
            buf.extend_from_slice(&demotion.parent_id.to_le_bytes());
            push_field(&mut buf, demotion.reason.dataset.as_bytes());
            push_field(&mut buf, demotion.reason.axis.as_bytes());
            buf.extend_from_slice(&demotion.reason.value.to_bits().to_le_bytes());
        }
        match origin_bytes {
            Some(origin_bytes) => {
                buf.push(1);
                buf.extend_from_slice(origin_bytes);
            }
            None => buf.push(0),
        }
        match metadata.origin_policy_fingerprint {
            Some(policy) => {
                buf.push(1);
                buf.extend_from_slice(policy.as_bytes());
            }
            None => buf.push(0),
        }
        push_len(&mut buf, metadata.waiver_dependencies.len());
        for dependency in metadata.waiver_dependencies {
            buf.extend_from_slice(&dependency.authorizing_node.to_le_bytes());
            match dependency.operation {
                Some(op) => {
                    buf.push(1);
                    buf.push(interval_op_tag(op));
                }
                None => buf.push(0),
            }
            push_field(
                &mut buf,
                &dependency.grant.signing_payload_for(dependency.operation),
            );
            push_field(&mut buf, &dependency.grant.signature);
            buf.extend_from_slice(dependency.policy_fingerprint.as_bytes());
            buf.extend_from_slice(&dependency.admission_day.to_le_bytes());
        }
        match metadata.waiver {
            Some(w) => {
                buf.push(1);
                push_field(&mut buf, w.id.as_bytes());
                push_field(&mut buf, w.signer.as_bytes());
                push_field(&mut buf, w.reason.as_bytes());
            }
            None => buf.push(0),
        }
        match metadata.grant {
            Some(g) => {
                buf.push(1);
                let payload = g.signing_payload_for(metadata.operation);
                push_field(&mut buf, &payload);
                push_field(&mut buf, &g.signature);
            }
            None => buf.push(0),
        }
        match metadata.waiver_policy_fingerprint {
            Some(policy) => {
                buf.push(1);
                buf.extend_from_slice(policy.as_bytes());
            }
            None => buf.push(0),
        }
        match metadata.waiver_admission_day {
            Some(day) => {
                buf.push(1);
                buf.extend_from_slice(&day.to_le_bytes());
            }
            None => buf.push(0),
        }
        hash_bytes(&buf)
    }

    fn push_node(
        &mut self,
        name: &str,
        color: Color,
        parents: Vec<u64>,
        metadata: NodeWriteMetadata,
    ) -> u64 {
        let id = self.nodes.len() as u64;
        let hash = self.node_hash(name, &color, &parents, &NodeHashMetadata::from(&metadata));
        let NodeWriteMetadata {
            operation,
            demotions,
            origin,
            origin_policy_fingerprint,
            waiver,
            grant,
            waiver_policy_fingerprint,
            waiver_admission_day,
            waiver_dependencies,
        } = metadata;
        // EVERY demotion is an event row, in canonical (parent write
        // order) sequence, each naming the demoted parent — losing all
        // but the first demotion loses decision-relevant diagnostics
        // (bead gp3.16).
        for demotion in &demotions {
            let d = &demotion.reason;
            self.rows.push(format!(
                "{{\"event\":\"demotion\",\"schema_version\":{COLOR_DEMOTION_ROW_SCHEMA_VERSION},\"node\":{id},\"parent_index\":{},\
                 \"parent\":{},\
                 \"dataset\":{},\"axis\":{},\"value\":{},\"value_bits\":\"{:016x}\"}}",
                demotion.parent_index,
                demotion.parent_id,
                json_string(&d.dataset),
                json_string(&d.axis),
                json_f64(d.value),
                d.value.to_bits(),
            ));
        }
        let operation_json =
            operation.map_or("null".to_string(), |op| json_string(interval_op_name(op)));
        let color_canonical_hex = hex_bytes(&color.canonical_bytes());
        let origin_canonical_bytes = origin.as_ref().map(source_origin_canonical_bytes);
        let origin_canonical_hex = optional_hex_bytes_json(origin_canonical_bytes.as_deref());
        self.rows.push(format!(
            "{{\"event\":\"color-write\",\"schema_version\":{COLOR_WRITE_ROW_SCHEMA_VERSION},\
             \"node_hash_version\":{COLOR_NODE_HASH_ENCODING_VERSION},\
             \"color_algebra_version\":{COLOR_ALGEBRA_VERSION},\"node\":{id},\
             \"name\":{},\"operation\":{},\"color\":\"{}\",\"payload\":{},\
             \"color_canonical_hex\":\"{}\",\"parents\":{:?},\"origin\":{},\
             \"origin_canonical_hex\":{},\"origin_policy_fingerprint\":{},\
             \"waiver_dependencies\":[{}],\"waiver\":{},\"grant\":{},\
             \"waiver_policy_fingerprint\":{},\"waiver_admission_day\":{},\
             \"hash\":\"{}\"}}",
            json_string(name),
            operation_json,
            color.name(),
            color.payload_json(),
            color_canonical_hex,
            parents,
            origin_json(origin.as_ref()),
            origin_canonical_hex,
            optional_hash_json(origin_policy_fingerprint),
            waiver_dependencies_json(&waiver_dependencies),
            waiver_json(waiver.as_ref()),
            grant_json(grant.as_ref(), operation),
            optional_hash_json(waiver_policy_fingerprint),
            optional_u32_json(waiver_admission_day),
            hash.to_hex()
        ));
        self.nodes.push(ColorNode {
            id,
            name: name.to_string(),
            color,
            parents,
            operation,
            demotions,
            origin,
            origin_policy_fingerprint,
            waiver,
            grant,
            waiver_policy_fingerprint,
            waiver_admission_day,
            waiver_dependencies,
            hash,
        });
        id
    }

    /// Write an ESTIMATED leaf (a surrogate, a heuristic, an estimator
    /// output). Estimates state their own dispersion and need no
    /// origin, but the estimator identity must be meaningful and the
    /// dispersion must be nonnegative/non-NaN (positive infinity is the
    /// explicit no-spread-claim sentinel). POSITIVE colors (Validated, Verified) are REFUSED here
    /// (bead gp3.16): a leaf cannot assert a certificate into
    /// existence — carry the minting evidence via
    /// [`ColorGraph::source_with_origin`] or an authenticated grant via
    /// [`ColorGraph::source_waived`].
    ///
    /// # Errors
    /// [`ColorWriteError::SourceOriginRequired`] for positive colors;
    /// [`ColorWriteError::InvalidEstimatedSource`] for malformed estimates.
    pub fn source(&mut self, name: &str, color: Color) -> Result<u64, ColorWriteError> {
        validate_node_name(name)?;
        if !matches!(color, Color::Estimated { .. }) {
            return Err(ColorWriteError::SourceOriginRequired { rank: color.rank() });
        }
        if let Some((field, why)) = estimated_source_payload_error(&color) {
            return Err(ColorWriteError::InvalidEstimatedSource { field, why });
        }
        Ok(self.push_node(
            name,
            color,
            Vec::new(),
            NodeWriteMetadata {
                operation: None,
                demotions: Vec::new(),
                origin: None,
                origin_policy_fingerprint: None,
                waiver: None,
                grant: None,
                waiver_policy_fingerprint: None,
                waiver_admission_day: None,
                waiver_dependencies: Vec::new(),
            },
        ))
    }

    /// Write a POSITIVE-colored leaf from TYPED origin evidence (bead
    /// gp3.16). The origin is the minting INPUT, not a memo: a Verified
    /// claim is re-derived from the carried certificate through
    /// [`fs_evidence::verified_from`] and must match bit-exactly; a
    /// Validated claim is reconstructed from the origin's anchoring
    /// dataset and exact regime. Because all evidence fields are public
    /// data, the injected [`SourceOriginVerifier`] must also authenticate
    /// the complete node/color/origin request. The origin participates in the provenance hash — substituting it
    /// later changes the node identity and every downstream hash.
    ///
    /// # Errors
    /// [`ColorWriteError::SourceOriginRefused`] with the structured
    /// forged-source reason, or [`ColorWriteError::InvalidColor`] when the
    /// rederived claim exceeds structural/resource limits.
    pub fn source_with_origin(
        &mut self,
        name: &str,
        color: &Color,
        origin: SourceOrigin,
        verifier: &dyn SourceOriginVerifier,
    ) -> Result<u64, ColorWriteError> {
        validate_node_name(name)?;
        let refuse = |rejection| Err(ColorWriteError::SourceOriginRefused { rejection });
        if matches!(color, Color::Estimated { .. }) {
            return refuse(SourceOriginRejection::EstimatedNeedsNoOrigin);
        }
        validate_color_structure(color)
            .map_err(|rejection| ColorWriteError::InvalidColor { rejection })?;
        validate_source_origin_resource_limits(&origin)?;
        let derived = origin
            .derive_color()
            .map_err(|rejection| ColorWriteError::SourceOriginRefused { rejection })?;
        if derived.canonical_bytes() != color.canonical_bytes() {
            let rejection = match (&derived, color) {
                (Color::Verified { .. }, Color::Verified { .. }) => {
                    SourceOriginRejection::CertificateMismatch
                }
                (
                    Color::Validated {
                        dataset: origin_dataset,
                        ..
                    },
                    Color::Validated {
                        dataset: color_dataset,
                        ..
                    },
                ) if origin_dataset != color_dataset => SourceOriginRejection::DatasetMismatch {
                    origin: origin_dataset.clone(),
                    color: color_dataset.clone(),
                },
                (Color::Validated { .. }, Color::Validated { .. }) => {
                    SourceOriginRejection::RegimeMismatch
                }
                _ => SourceOriginRejection::OriginKindMismatch {
                    color: color.name(),
                },
            };
            return refuse(rejection);
        }
        let request = SourceOriginRequest::new(name, color, &origin);
        let decision =
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| verifier.verify(&request)))
                .map_err(|_| ColorWriteError::SourceOriginRefused {
                    rejection: SourceOriginRejection::VerifierPanicked,
                })?;
        if !decision.accepted() {
            return refuse(SourceOriginRejection::VerifierRefused {
                policy_fingerprint: decision.policy_fingerprint(),
            });
        }
        Ok(self.push_node(
            name,
            derived,
            Vec::new(),
            NodeWriteMetadata {
                operation: None,
                demotions: Vec::new(),
                origin: Some(origin),
                origin_policy_fingerprint: Some(decision.policy_fingerprint()),
                waiver: None,
                grant: None,
                waiver_policy_fingerprint: None,
                waiver_admission_day: None,
                waiver_dependencies: Vec::new(),
            },
        ))
    }

    /// Write a POSITIVE-colored leaf authorized by an AUTHENTICATED
    /// [`WaiverGrant`] carrying the SOURCE-COLOR scope (bead gp3.16) —
    /// the human-responsibility door when typed origin evidence does
    /// not exist. The grant must name THIS node, authorize exactly the
    /// claimed color bytes, carry an EMPTY lineage (a leaf has no
    /// parents — a grant minted for a derive cannot be replayed here),
    /// be unexpired, and verify over the v4 source signing payload.
    /// Fail closed on any mismatch.
    ///
    /// # Errors
    /// [`ColorWriteError::InvalidColor`] if the claimed payload is malformed;
    /// [`ColorWriteError::WaiverRefused`] with the structured
    /// rejection; [`ColorWriteError::SourceOriginRequired`] doctrine
    /// does not apply here (this IS the waiver path), but Estimated
    /// leaves are refused via
    /// [`SourceOriginRejection::EstimatedNeedsNoOrigin`].
    pub fn source_waived(
        &mut self,
        name: &str,
        color: Color,
        grant: WaiverGrant,
        verifier: &dyn WaiverVerifier,
        today_day: u32,
    ) -> Result<u64, ColorWriteError> {
        validate_node_name(name)?;
        validate_color_structure(&color)
            .map_err(|rejection| ColorWriteError::InvalidColor { rejection })?;
        if color.rank() < ColorRank::Validated {
            return Err(ColorWriteError::SourceOriginRefused {
                rejection: SourceOriginRejection::EstimatedNeedsNoOrigin,
            });
        }
        let refuse = |rejection| Err(ColorWriteError::WaiverRefused { rejection });
        validate_waiver_grant(&grant)
            .map_err(|rejection| ColorWriteError::WaiverRefused { rejection })?;
        validate_waiver_closure_bytes(&[], Some(&grant))?;
        if grant.scope != WAIVER_SCOPE_SOURCE_COLOR {
            return refuse(WaiverRejection::ScopeMismatch);
        }
        if grant.node_name != name {
            return refuse(WaiverRejection::NodeMismatch);
        }
        if grant.claimed_color != color.canonical_bytes() {
            return refuse(WaiverRejection::ColorMismatch);
        }
        if !grant.parent_hashes.is_empty() {
            return refuse(WaiverRejection::LineageMismatch);
        }
        if today_day > grant.expires_day {
            return refuse(WaiverRejection::Expired);
        }
        let decision = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            verifier.verify(
                &grant.key_id,
                &grant.signing_payload_source(),
                &grant.signature,
            )
        }))
        .map_err(|_| ColorWriteError::WaiverRefused {
            rejection: WaiverRejection::VerifierPanicked,
        })?;
        if !decision.accepted() {
            return refuse(WaiverRejection::VerifierRefused {
                policy_fingerprint: decision.policy_fingerprint(),
            });
        }
        Ok(self.push_node(
            name,
            color,
            Vec::new(),
            NodeWriteMetadata {
                operation: None,
                demotions: Vec::new(),
                origin: None,
                origin_policy_fingerprint: None,
                waiver: Some(grant.annotation.clone()),
                grant: Some(grant),
                waiver_policy_fingerprint: Some(decision.policy_fingerprint()),
                waiver_admission_day: Some(today_day),
                waiver_dependencies: Vec::new(),
            },
        ))
    }

    /// Regime re-checks + composition fold shared by the derive paths.
    /// EVERY demotion is collected (bead gp3.16), with both parent id
    /// and position in canonical ascending-position order. Retaining only the first demotion loses
    /// decision-relevant diagnostics when several parents exit their
    /// regimes at once. Effective Validated axes are preflighted into one
    /// bounded map before any parent color is cloned or composed.
    fn fold_parents(
        &self,
        parents: &[u64],
        op: IntervalOp,
        state: &BTreeMap<String, f64>,
    ) -> Result<(Color, Vec<ColorDemotion>), ColorWriteError> {
        if parents.is_empty() {
            return Err(ColorWriteError::NoParents);
        }
        if parents.len() > MAX_COLOR_PARENTS {
            return Err(ColorWriteError::ResourceLimitExceeded {
                resource: "parents",
                limit: MAX_COLOR_PARENTS,
                actual: parents.len(),
            });
        }
        let parent_nodes = parents
            .iter()
            .map(|&id| self.node(id).ok_or(ColorWriteError::UnknownParent { id }))
            .collect::<Result<Vec<_>, _>>()?;
        let mut demotions = Vec::new();
        let mut aggregate_validity = BTreeMap::new();
        let mut estimated_absorbed = false;
        for (parent_index, (parent_id, parent)) in parents
            .iter()
            .copied()
            .zip(parent_nodes.iter().copied())
            .enumerate()
        {
            if let Some(reason) = regime_demotion(&parent.color, state) {
                demotions.push(ColorDemotion {
                    parent_index,
                    parent_id,
                    reason,
                });
                estimated_absorbed = true;
                aggregate_validity.clear();
            } else {
                preflight_fold_color(
                    &mut aggregate_validity,
                    &mut estimated_absorbed,
                    &parent.color,
                )?;
            }
        }
        let mut next_demotion = 0usize;
        let mut derived = None;
        for (parent_index, parent) in parent_nodes.into_iter().enumerate() {
            let effective = if demotions
                .get(next_demotion)
                .is_some_and(|demotion| demotion.parent_index == parent_index)
            {
                let reason = &demotions[next_demotion].reason;
                next_demotion += 1;
                Color::Estimated {
                    estimator: demotion_estimator_identity(&reason.dataset, &reason.axis),
                    dispersion: f64::INFINITY,
                }
            } else {
                parent.color.clone()
            };
            derived = Some(match derived {
                None => effective,
                Some(current) => compose(&current, &effective, op),
            });
        }
        derived
            .map(|color| (color, demotions))
            .ok_or(ColorWriteError::NoParents)
    }

    fn laundering_error(
        &self,
        parents: &[u64],
        state: &BTreeMap<String, f64>,
        claimed: ColorRank,
        cap: ColorRank,
    ) -> ColorWriteError {
        let offending: Vec<u64> = parents
            .iter()
            .copied()
            .filter(|&p| {
                let parent = self
                    .node(p)
                    .expect("laundering parents were validated by fold_parents");
                let effective_rank = if regime_demotion(&parent.color, state).is_some() {
                    ColorRank::Estimated
                } else {
                    parent.color.rank()
                };
                effective_rank <= cap
            })
            .collect();
        ColorWriteError::LaunderingRefused {
            claimed,
            derived: cap,
            offending_parents: offending,
        }
    }

    /// Write a DERIVED node: the composition algebra folds the parent
    /// colors (with regime re-checks against `state`, auto-demoting on
    /// exit), and any explicit claimed color must equal that exact result.
    /// Rank-only weakening is not accepted because the payload may still
    /// narrow an interval, widen a regime, or shrink dispersion.
    /// The `waiver` argument is a HUMAN ANNOTATION only (bead
    /// qmao.1.1): it is recorded and hashed but authorizes NOTHING —
    /// an upgrade claim is refused here regardless. The authorized
    /// path is [`ColorGraph::derive_waived`].
    ///
    /// # Errors
    /// [`ColorWriteError`] teaching errors; the laundering refusal
    /// names the capping parents.
    pub fn derive(
        &mut self,
        name: &str,
        parents: &[u64],
        op: IntervalOp,
        claimed: Option<Color>,
        state: &BTreeMap<String, f64>,
        waiver: Option<Waiver>,
    ) -> Result<u64, ColorWriteError> {
        validate_node_name(name)?;
        if let Some(claimed) = &claimed {
            validate_color_structure(claimed)
                .map_err(|rejection| ColorWriteError::InvalidColor { rejection })?;
        }
        if let Some(waiver) = &waiver
            && let Some((field, reason)) = waiver_annotation_reason(waiver)
        {
            return Err(ColorWriteError::InvalidWaiverAnnotation { field, reason });
        }
        let (derived, demotions) = self.fold_parents(parents, op, state)?;
        let waiver_dependencies = self.inherited_waiver_dependencies(parents, None)?;
        let written = match claimed {
            None => derived,
            Some(c) if c.canonical_bytes() == derived.canonical_bytes() => c,
            Some(c) if c.rank() > derived.rank() => {
                return Err(self.laundering_error(parents, state, c.rank(), derived.rank()));
            }
            Some(c) => {
                return Err(ColorWriteError::ClaimMismatch {
                    claimed: c,
                    derived,
                });
            }
        };
        validate_color_structure(&written)
            .map_err(|rejection| ColorWriteError::InvalidColor { rejection })?;
        Ok(self.push_node(
            name,
            written,
            parents.to_vec(),
            NodeWriteMetadata {
                operation: Some(op),
                demotions,
                origin: None,
                origin_policy_fingerprint: None,
                waiver,
                grant: None,
                waiver_policy_fingerprint: None,
                waiver_admission_day: None,
                waiver_dependencies,
            },
        ))
    }

    /// Write a DERIVED node whose claim is authorized by an AUTHENTICATED
    /// [`WaiverGrant`] (bead qmao.1.1):
    /// the grant must carry the color-upgrade scope, name THIS node,
    /// authorize exactly the claimed color, bind the exact parent
    /// provenance hashes and exact operation (replay to another node,
    /// lineage, or operation fails), be unexpired
    /// as of `today_day`, and carry a signature the `verifier`
    /// capability accepts over the canonical length-prefixed payload.
    /// Any failure refuses the write (fail closed) — with the in-tree
    /// [`NoWaiverVerifier`] every promotion is refused (the no-crypto
    /// no-claim).
    ///
    /// # Errors
    /// [`ColorWriteError::InvalidColor`] if the claimed payload is malformed;
    /// [`ColorWriteError::WaiverRefused`] with the structured
    /// rejection, plus the ordinary derive errors.
    #[allow(clippy::too_many_arguments)] // the authorization surface is the point
    pub fn derive_waived(
        &mut self,
        name: &str,
        parents: &[u64],
        op: IntervalOp,
        claimed: Color,
        state: &BTreeMap<String, f64>,
        grant: WaiverGrant,
        verifier: &dyn WaiverVerifier,
        today_day: u32,
    ) -> Result<u64, ColorWriteError> {
        validate_node_name(name)?;
        validate_color_structure(&claimed)
            .map_err(|rejection| ColorWriteError::InvalidColor { rejection })?;
        validate_waiver_grant(&grant)
            .map_err(|rejection| ColorWriteError::WaiverRefused { rejection })?;
        let (_derived, demotions) = self.fold_parents(parents, op, state)?;
        let waiver_dependencies = self.inherited_waiver_dependencies(parents, Some(&grant))?;
        let refuse = |rejection| Err(ColorWriteError::WaiverRefused { rejection });
        if grant.scope != WAIVER_SCOPE_COLOR_UPGRADE {
            return refuse(WaiverRejection::ScopeMismatch);
        }
        if grant.node_name != name {
            return refuse(WaiverRejection::NodeMismatch);
        }
        if grant.claimed_color != claimed.canonical_bytes() {
            return refuse(WaiverRejection::ColorMismatch);
        }
        let lineage: Vec<ContentHash> = parents
            .iter()
            .map(|&p| {
                self.node(p)
                    .expect("waived parents were validated by fold_parents")
                    .hash
            })
            .collect();
        if grant.parent_hashes != lineage {
            return refuse(WaiverRejection::LineageMismatch);
        }
        if today_day > grant.expires_day {
            return refuse(WaiverRejection::Expired);
        }
        let decision = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            verifier.verify(&grant.key_id, &grant.signing_payload(op), &grant.signature)
        }))
        .map_err(|_| ColorWriteError::WaiverRefused {
            rejection: WaiverRejection::VerifierPanicked,
        })?;
        if !decision.accepted() {
            return refuse(WaiverRejection::VerifierRefused {
                policy_fingerprint: decision.policy_fingerprint(),
            });
        }
        Ok(self.push_node(
            name,
            claimed,
            parents.to_vec(),
            NodeWriteMetadata {
                operation: Some(op),
                demotions,
                origin: None,
                origin_policy_fingerprint: None,
                waiver: Some(grant.annotation.clone()),
                grant: Some(grant),
                waiver_policy_fingerprint: Some(decision.policy_fingerprint()),
                waiver_admission_day: Some(today_day),
                waiver_dependencies,
            },
        ))
    }

    /// The node by id — CHECKED (bead gp3.16): an invalid public id is
    /// a caller error to surface, not a panic to detonate.
    #[must_use]
    pub fn node(&self, id: u64) -> Option<&ColorNode> {
        self.nodes.get(usize::try_from(id).ok()?)
    }

    fn replay_error(node: &ColorNode, why: impl Into<String>) -> ColorReplayError {
        ColorReplayError {
            node: node.id,
            why: why.into(),
        }
    }

    fn validate_replay_waiver_dependencies(
        &self,
        node: &ColorNode,
    ) -> Result<(), ColorReplayError> {
        if node.waiver_dependencies.len() > MAX_WAIVER_DEPENDENCIES {
            return Err(Self::replay_error(
                node,
                "waiver dependency closure exceeds the replay limit",
            ));
        }
        validate_waiver_closure_bytes(&node.waiver_dependencies, node.grant.as_ref()).map_err(
            |error| Self::replay_error(node, format!("waiver closure is oversized: {error}")),
        )?;
        let mut previous = None;
        for dependency in &node.waiver_dependencies {
            if dependency.authorizing_node >= node.id {
                return Err(Self::replay_error(
                    node,
                    "waiver dependency authorizing node is self/forward, not strictly earlier",
                ));
            }
            if previous.is_some_and(|previous| previous >= dependency.authorizing_node) {
                return Err(Self::replay_error(
                    node,
                    "waiver dependencies are duplicated or not in canonical ascending order",
                ));
            }
            previous = Some(dependency.authorizing_node);
            let Some(authorizing_node) = self.node(dependency.authorizing_node) else {
                return Err(Self::replay_error(
                    node,
                    "waiver dependency authorizing node is missing",
                ));
            };
            if authorizing_node.operation != dependency.operation
                || authorizing_node.grant.as_ref() != Some(&dependency.grant)
                || authorizing_node.waiver_policy_fingerprint != Some(dependency.policy_fingerprint)
                || authorizing_node.waiver_admission_day != Some(dependency.admission_day)
            {
                return Err(Self::replay_error(
                    node,
                    "waiver dependency differs from its historical authorizing node",
                ));
            }
            validate_waiver_grant(&dependency.grant).map_err(|rejection| {
                Self::replay_error(
                    node,
                    format!("waiver dependency has invalid grant metadata: {rejection}"),
                )
            })?;
        }

        let expected = if node.parents.is_empty() {
            Ok(Vec::new())
        } else {
            self.inherited_waiver_dependencies(&node.parents, node.grant.as_ref())
        }
        .map_err(|error| Self::replay_error(node, format!("invalid parent authority: {error}")))?;
        if node.waiver_dependencies != expected {
            return Err(Self::replay_error(
                node,
                "waiver dependency closure differs from the exact parent-derived union",
            ));
        }
        Ok(())
    }

    fn validate_replay_demotions(&self, node: &ColorNode) -> Result<(), ColorReplayError> {
        let mut previous_index = None;
        for demotion in &node.demotions {
            if previous_index.is_some_and(|previous| previous >= demotion.parent_index) {
                return Err(Self::replay_error(
                    node,
                    "demotions are not in unique ascending parent-position order",
                ));
            }
            previous_index = Some(demotion.parent_index);
            if node.parents.get(demotion.parent_index) != Some(&demotion.parent_id) {
                return Err(Self::replay_error(
                    node,
                    "demotion parent position and id disagree",
                ));
            }
            let Some(Color::Validated { regime, dataset }) = self
                .node(demotion.parent_id)
                .map(ColorNode::declared_color_unverified)
            else {
                return Err(Self::replay_error(
                    node,
                    "demotion names a parent that is not Validated",
                ));
            };
            if dataset != &demotion.reason.dataset {
                return Err(Self::replay_error(
                    node,
                    "demotion dataset differs from its parent anchor",
                ));
            }
            let value = demotion.reason.value;
            if regime.bounds().is_empty() {
                if demotion.reason.axis != "<undeclared-regime>" || value.is_finite() {
                    return Err(Self::replay_error(
                        node,
                        "empty-regime demotion is not the canonical sentinel",
                    ));
                }
            } else if let Some((lo, hi)) = regime.bound(&demotion.reason.axis) {
                if lo.is_finite()
                    && hi.is_finite()
                    && lo <= hi
                    && value.is_finite()
                    && value >= lo
                    && value <= hi
                {
                    return Err(Self::replay_error(
                        node,
                        "demotion value remains inside its parent regime",
                    ));
                }
            } else {
                return Err(Self::replay_error(
                    node,
                    "demotion axis is absent from its parent regime",
                ));
            }
        }
        Ok(())
    }

    fn validate_replay_source(node: &ColorNode) -> Result<(), ColorReplayError> {
        if node.operation.is_some() || !node.demotions.is_empty() {
            return Err(Self::replay_error(
                node,
                "source leaf carries an operation or demotion",
            ));
        }
        match (&node.color, &node.origin, &node.grant) {
            (Color::Estimated { .. }, None, None) => {
                if node.origin_policy_fingerprint.is_some()
                    || node.waiver.is_some()
                    || node.waiver_policy_fingerprint.is_some()
                    || node.waiver_admission_day.is_some()
                {
                    return Err(Self::replay_error(
                        node,
                        "Estimated source carries orphan authority or human-waiver metadata",
                    ));
                }
                if let Some((field, why)) = estimated_source_payload_error(&node.color) {
                    Err(Self::replay_error(
                        node,
                        format!("Estimated source field `{field}` is invalid: {why}"),
                    ))
                } else {
                    Ok(())
                }
            }
            (Color::Estimated { .. }, _, _) => Err(Self::replay_error(
                node,
                "Estimated leaf must not carry source authority",
            )),
            (_, Some(origin), None) => {
                if node.origin_policy_fingerprint.is_none()
                    || node.waiver_policy_fingerprint.is_some()
                    || node.waiver_admission_day.is_some()
                {
                    return Err(Self::replay_error(
                        node,
                        "typed-origin source lacks its policy or carries waiver context",
                    ));
                }
                let derived = origin.derive_color().map_err(|rejection| {
                    Self::replay_error(
                        node,
                        format!("typed source origin no longer mints: {rejection}"),
                    )
                })?;
                if derived.canonical_bytes() != node.color.canonical_bytes() {
                    return Err(Self::replay_error(
                        node,
                        "typed source origin does not rederive the stored color",
                    ));
                }
                if node.waiver.is_some() {
                    return Err(Self::replay_error(
                        node,
                        "typed-origin source also carries an unrelated waiver",
                    ));
                }
                Ok(())
            }
            (_, None, Some(grant)) => {
                validate_waiver_grant(grant).map_err(|rejection| {
                    Self::replay_error(
                        node,
                        format!("source grant metadata is invalid: {rejection}"),
                    )
                })?;
                if node.origin_policy_fingerprint.is_some()
                    || node.waiver_policy_fingerprint.is_none()
                    || node.waiver_admission_day.is_none()
                    || node
                        .waiver_admission_day
                        .is_some_and(|day| day > grant.expires_day)
                    || grant.scope != WAIVER_SCOPE_SOURCE_COLOR
                    || grant.node_name != node.name
                    || grant.claimed_color != node.color.canonical_bytes()
                    || !grant.parent_hashes.is_empty()
                    || node.waiver.as_ref() != Some(&grant.annotation)
                {
                    return Err(Self::replay_error(
                        node,
                        "source grant fields do not bind the stored leaf",
                    ));
                }
                Ok(())
            }
            (_, Some(_), Some(_)) => Err(Self::replay_error(
                node,
                "source leaf carries both typed origin and waiver authority",
            )),
            (_, None, None) => Err(Self::replay_error(
                node,
                "positive-colored leaf carries neither typed origin nor grant",
            )),
        }
    }

    fn validate_replay_derived(&self, node: &ColorNode) -> Result<(), ColorReplayError> {
        let Some(op) = node.operation else {
            return Err(Self::replay_error(
                node,
                "derived node lacks a composition operation",
            ));
        };
        if node.origin.is_some() {
            return Err(Self::replay_error(
                node,
                "derived node carries a source-only origin",
            ));
        }
        if node.origin_policy_fingerprint.is_some() {
            return Err(Self::replay_error(
                node,
                "derived node carries source admission policy context",
            ));
        }
        let mut aggregate_validity = BTreeMap::new();
        let mut estimated_absorbed = false;
        let mut derived = None;
        let mut next_demotion = 0usize;
        for (index, parent) in node.parents.iter().enumerate() {
            let Some(parent_node) = self.node(*parent) else {
                return Err(Self::replay_error(node, "derived parent is missing"));
            };
            let effective = if node
                .demotions
                .get(next_demotion)
                .is_some_and(|demotion| demotion.parent_index == index)
            {
                let reason = &node.demotions[next_demotion].reason;
                next_demotion += 1;
                Color::Estimated {
                    estimator: demotion_estimator_identity(&reason.dataset, &reason.axis),
                    dispersion: f64::INFINITY,
                }
            } else {
                parent_node.color.clone()
            };
            preflight_fold_color(&mut aggregate_validity, &mut estimated_absorbed, &effective)
                .map_err(|error| {
                    Self::replay_error(node, format!("derived validity preflight failed: {error}"))
                })?;
            derived = Some(match derived {
                None => effective,
                Some(current) => compose(&current, &effective, op),
            });
        }
        let Some(derived) = derived else {
            return Err(Self::replay_error(node, "derived node has no parents"));
        };
        if let Some(grant) = &node.grant {
            validate_waiver_grant(grant).map_err(|rejection| {
                Self::replay_error(
                    node,
                    format!("derived grant metadata is invalid: {rejection}"),
                )
            })?;
            let mut lineage = Vec::with_capacity(node.parents.len());
            for parent in &node.parents {
                let Some(parent_node) = self.node(*parent) else {
                    return Err(Self::replay_error(node, "derived parent is missing"));
                };
                lineage.push(parent_node.hash);
            }
            if node.waiver_policy_fingerprint.is_none()
                || node.waiver_admission_day.is_none()
                || node
                    .waiver_admission_day
                    .is_some_and(|day| day > grant.expires_day)
                || grant.scope != WAIVER_SCOPE_COLOR_UPGRADE
                || grant.node_name != node.name
                || grant.claimed_color != node.color.canonical_bytes()
                || grant.parent_hashes != lineage
                || node.waiver.as_ref() != Some(&grant.annotation)
            {
                return Err(Self::replay_error(
                    node,
                    "derived grant fields do not bind the stored node",
                ));
            }
        } else {
            if node.waiver_policy_fingerprint.is_some() || node.waiver_admission_day.is_some() {
                return Err(Self::replay_error(
                    node,
                    "ordinary derived node carries orphan waiver admission context",
                ));
            }
            if derived.canonical_bytes() != node.color.canonical_bytes() {
                return Err(Self::replay_error(
                    node,
                    "written color does not rederive from parents and demotions",
                ));
            }
        }
        Ok(())
    }

    fn verify_replay_node(
        &self,
        position: usize,
        node: &ColorNode,
    ) -> Result<(), ColorReplayError> {
        if usize::try_from(node.id).ok() != Some(position) {
            return Err(Self::replay_error(
                node,
                "stored id differs from append position",
            ));
        }
        if let Some(reason) = identity_reason(&node.name) {
            return Err(Self::replay_error(
                node,
                format!("stored node name is invalid: {reason}"),
            ));
        }
        if let Err(rejection) = validate_color_structure(&node.color) {
            return Err(Self::replay_error(
                node,
                format!("stored color is structurally invalid: {rejection}"),
            ));
        }
        if let Some(waiver) = &node.waiver
            && let Some((field, reason)) = waiver_annotation_reason(waiver)
        {
            return Err(Self::replay_error(
                node,
                format!("stored waiver annotation field `{field}` is invalid ({reason})"),
            ));
        }
        if node.parents.iter().any(|parent| {
            usize::try_from(*parent)
                .ok()
                .is_none_or(|parent| parent >= position)
        }) {
            return Err(Self::replay_error(
                node,
                "parent id is missing or does not precede the derived node",
            ));
        }
        if node.parents.len() > MAX_COLOR_PARENTS {
            return Err(Self::replay_error(
                node,
                "parent list exceeds the replay limit",
            ));
        }
        self.validate_replay_waiver_dependencies(node)?;
        self.validate_replay_demotions(node)?;
        let metadata = NodeHashMetadata {
            operation: node.operation,
            demotions: &node.demotions,
            origin: node.origin.as_ref(),
            origin_policy_fingerprint: node.origin_policy_fingerprint,
            waiver: node.waiver.as_ref(),
            grant: node.grant.as_ref(),
            waiver_policy_fingerprint: node.waiver_policy_fingerprint,
            waiver_admission_day: node.waiver_admission_day,
            waiver_dependencies: &node.waiver_dependencies,
        };
        if self.node_hash(&node.name, &node.color, &node.parents, &metadata) != node.hash {
            return Err(Self::replay_error(
                node,
                "provenance hash does not rederive from the stored fields",
            ));
        }
        if node.parents.is_empty() {
            Self::validate_replay_source(node)
        } else {
            self.validate_replay_derived(node)
        }
    }

    /// IN-MEMORY STRUCTURAL REPLAY AUDIT (bead gp3.16): rederive every node from its stored
    /// inputs and refuse on any divergence. For each derived node the
    /// recorded demotions reconstruct the effective parent colors
    /// (a demotion determines the bounded, length-framed estimator identity
    /// exactly), the composition
    /// algebra re-folds them, and — for unwaived writes — the written
    /// color must match bit-exactly. Every node's provenance hash is
    /// recomputed and compared, so the graph's whole hash chain is
    /// re-earned, never trusted. Positive-colored leaves must carry
    /// their typed origin or a structurally bound historical grant (the
    /// sealed-source invariant, re-checked). This method does not parse
    /// persisted rows, resolve policy fingerprints, re-run external authority
    /// capabilities, or apply a new current-day expiry decision.
    ///
    /// # Errors
    /// [`ColorReplayError`] naming the first diverging node.
    pub fn verify_replay(&self) -> Result<(), ColorReplayError> {
        for (position, node) in self.nodes.iter().enumerate() {
            self.verify_replay_node(position, node)?;
        }
        Ok(())
    }
}

/// A replay-audit divergence: the first node whose stored state does
/// not rederive from its inputs.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ColorReplayError {
    /// The diverging node id.
    pub node: u64,
    /// What failed to rederive.
    pub why: String,
}

impl core::fmt::Display for ColorReplayError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "color replay audit failed at node {}: {}",
            self.node, self.why
        )
    }
}

impl std::error::Error for ColorReplayError {}

/// Evidence that a persisted color-row stream was independently reconstructed
/// and rehashed without consulting a [`ColorGraph`] or any in-memory [`Color`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ColorRowVerification {
    node_count: usize,
    demotion_count: usize,
    terminal_hash: Option<ContentHash>,
}

impl ColorRowVerification {
    /// Number of schema-v7 `color-write` rows whose hashes were re-earned.
    #[must_use]
    pub fn node_count(self) -> usize {
        self.node_count
    }

    /// Number of schema-v1 demotion rows incorporated into those preimages.
    #[must_use]
    pub fn demotion_count(self) -> usize {
        self.demotion_count
    }

    /// Hash of the final accepted node, or `None` for an empty stream.
    #[must_use]
    pub fn terminal_hash(self) -> Option<ContentHash> {
        self.terminal_hash
    }
}

/// A fail-closed refusal from [`verify_color_row_stream`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ColorRowVerificationError {
    /// A row is not strict JSON or violates the canonical row shape/order.
    Malformed {
        /// Zero-based position in the supplied stream. A value equal to the
        /// stream length denotes unfinished demotion rows at end of input.
        row: usize,
        /// Node named by the row, when it was decoded before the refusal.
        node: Option<u64>,
        /// Exact structural reason.
        why: String,
    },
    /// Prior or future color/demotion row schemas are never interpreted as the
    /// current preimage layout.
    UnsupportedSchemaVersion {
        /// Zero-based row position.
        row: usize,
        /// Node named by the row.
        node: u64,
        /// `"color-write"` or `"demotion"`.
        event: &'static str,
        /// Version found in storage.
        found: u64,
        /// Only version accepted by this build.
        expected: u64,
    },
    /// The row names a node-hash encoding this build cannot reconstruct.
    UnsupportedNodeHashVersion {
        /// Zero-based row position.
        row: usize,
        /// Node named by the row.
        node: u64,
        /// Version found in storage.
        found: u64,
        /// Only version accepted by this build.
        expected: u64,
    },
    /// The row names a color algebra this build cannot interpret.
    UnsupportedColorAlgebraVersion {
        /// Zero-based row position.
        row: usize,
        /// Node named by the row.
        node: u64,
        /// Version found in storage.
        found: u64,
        /// Only version accepted by this build.
        expected: u64,
    },
    /// Every preimage field decoded, but the independently computed hash did
    /// not match the row's claimed hash.
    HashMismatch {
        /// Zero-based row position.
        row: usize,
        /// Node named by the row.
        node: u64,
        /// Hash claimed by the persisted row.
        stored: ContentHash,
        /// Hash independently reconstructed from persisted fields.
        reconstructed: ContentHash,
    },
}

impl core::fmt::Display for ColorRowVerificationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Malformed { row, node, why } => {
                if let Some(node) = node {
                    write!(f, "color row {row} for node {node} is malformed: {why}")
                } else {
                    write!(f, "color row {row} is malformed: {why}")
                }
            }
            Self::UnsupportedSchemaVersion {
                row,
                node,
                event,
                found,
                expected,
            } => write!(
                f,
                "color row {row} for node {node} has unsupported {event} schema {found}; expected {expected}"
            ),
            Self::UnsupportedNodeHashVersion {
                row,
                node,
                found,
                expected,
            } => write!(
                f,
                "color row {row} for node {node} has unsupported node-hash version {found}; expected {expected}"
            ),
            Self::UnsupportedColorAlgebraVersion {
                row,
                node,
                found,
                expected,
            } => write!(
                f,
                "color row {row} for node {node} has unsupported color-algebra version {found}; expected {expected}"
            ),
            Self::HashMismatch {
                row,
                node,
                stored,
                reconstructed,
            } => write!(
                f,
                "color row {row} for node {node} claims hash {stored}, but persisted fields reconstruct {reconstructed}"
            ),
        }
    }
}

impl std::error::Error for ColorRowVerificationError {}

#[derive(Debug, Clone, PartialEq, Eq)]
enum RowJson {
    Null,
    Bool(bool),
    Number(String),
    String(String),
    Array(Vec<RowJson>),
    Object(BTreeMap<String, RowJson>),
}

struct RowJsonParser<'a> {
    bytes: &'a [u8],
    cursor: usize,
}

impl<'a> RowJsonParser<'a> {
    fn new(row: &'a str) -> Self {
        Self {
            bytes: row.as_bytes(),
            cursor: 0,
        }
    }

    fn parse(mut self) -> Result<RowJson, String> {
        self.skip_whitespace();
        let value = self.parse_value(0)?;
        self.skip_whitespace();
        if self.cursor != self.bytes.len() {
            return Err("trailing bytes after JSON value".to_string());
        }
        Ok(value)
    }

    fn parse_value(&mut self, depth: usize) -> Result<RowJson, String> {
        if depth > 64 {
            return Err("JSON nesting exceeds 64 levels".to_string());
        }
        self.skip_whitespace();
        match self.bytes.get(self.cursor).copied() {
            Some(b'n') => {
                self.consume_literal(b"null")?;
                Ok(RowJson::Null)
            }
            Some(b't') => {
                self.consume_literal(b"true")?;
                Ok(RowJson::Bool(true))
            }
            Some(b'f') => {
                self.consume_literal(b"false")?;
                Ok(RowJson::Bool(false))
            }
            Some(b'"') => self.parse_string().map(RowJson::String),
            Some(b'[') => self.parse_array(depth + 1),
            Some(b'{') => self.parse_object(depth + 1),
            Some(b'-' | b'0'..=b'9') => self.parse_number().map(RowJson::Number),
            Some(_) => Err("unexpected JSON token".to_string()),
            None => Err("unexpected end of JSON".to_string()),
        }
    }

    fn consume_literal(&mut self, literal: &[u8]) -> Result<(), String> {
        if self.bytes.get(self.cursor..self.cursor + literal.len()) == Some(literal) {
            self.cursor += literal.len();
            Ok(())
        } else {
            Err("malformed JSON literal".to_string())
        }
    }

    fn parse_array(&mut self, depth: usize) -> Result<RowJson, String> {
        self.cursor += 1;
        self.skip_whitespace();
        let mut values = Vec::new();
        if self.bytes.get(self.cursor) == Some(&b']') {
            self.cursor += 1;
            return Ok(RowJson::Array(values));
        }
        loop {
            values.push(self.parse_value(depth)?);
            self.skip_whitespace();
            match self.bytes.get(self.cursor) {
                Some(b',') => {
                    self.cursor += 1;
                    self.skip_whitespace();
                }
                Some(b']') => {
                    self.cursor += 1;
                    return Ok(RowJson::Array(values));
                }
                _ => return Err("JSON array lacks comma or closing bracket".to_string()),
            }
        }
    }

    fn parse_object(&mut self, depth: usize) -> Result<RowJson, String> {
        self.cursor += 1;
        self.skip_whitespace();
        let mut fields = BTreeMap::new();
        if self.bytes.get(self.cursor) == Some(&b'}') {
            self.cursor += 1;
            return Ok(RowJson::Object(fields));
        }
        loop {
            if self.bytes.get(self.cursor) != Some(&b'"') {
                return Err("JSON object key is not a string".to_string());
            }
            let key = self.parse_string()?;
            self.skip_whitespace();
            if self.bytes.get(self.cursor) != Some(&b':') {
                return Err("JSON object key lacks a colon".to_string());
            }
            self.cursor += 1;
            let value = self.parse_value(depth)?;
            if fields.insert(key, value).is_some() {
                return Err("JSON object contains a duplicate key".to_string());
            }
            self.skip_whitespace();
            match self.bytes.get(self.cursor) {
                Some(b',') => {
                    self.cursor += 1;
                    self.skip_whitespace();
                }
                Some(b'}') => {
                    self.cursor += 1;
                    return Ok(RowJson::Object(fields));
                }
                _ => return Err("JSON object lacks comma or closing brace".to_string()),
            }
        }
    }

    fn parse_string(&mut self) -> Result<String, String> {
        self.cursor += 1;
        let mut value = String::new();
        loop {
            let start = self.cursor;
            while self
                .bytes
                .get(self.cursor)
                .is_some_and(|byte| *byte != b'"' && *byte != b'\\' && *byte >= 0x20)
            {
                self.cursor += 1;
            }
            if self.cursor > start {
                let segment = core::str::from_utf8(&self.bytes[start..self.cursor])
                    .map_err(|_| "JSON string is not UTF-8".to_string())?;
                value.push_str(segment);
            }
            match self.bytes.get(self.cursor).copied() {
                Some(b'"') => {
                    self.cursor += 1;
                    return Ok(value);
                }
                Some(b'\\') => {
                    self.cursor += 1;
                    match self.bytes.get(self.cursor).copied() {
                        Some(b'"') => value.push('"'),
                        Some(b'\\') => value.push('\\'),
                        Some(b'/') => value.push('/'),
                        Some(b'b') => value.push('\u{0008}'),
                        Some(b'f') => value.push('\u{000c}'),
                        Some(b'n') => value.push('\n'),
                        Some(b'r') => value.push('\r'),
                        Some(b't') => value.push('\t'),
                        Some(b'u') => {
                            self.cursor += 1;
                            let first = self.parse_hex_quad()?;
                            let scalar = if (0xd800..=0xdbff).contains(&first) {
                                if self.bytes.get(self.cursor..self.cursor + 2) != Some(b"\\u") {
                                    return Err("high surrogate lacks a low surrogate".to_string());
                                }
                                self.cursor += 2;
                                let second = self.parse_hex_quad()?;
                                if !(0xdc00..=0xdfff).contains(&second) {
                                    return Err("invalid low surrogate".to_string());
                                }
                                0x1_0000
                                    + ((u32::from(first) - 0xd800) << 10)
                                    + (u32::from(second) - 0xdc00)
                            } else {
                                if (0xdc00..=0xdfff).contains(&first) {
                                    return Err("unpaired low surrogate".to_string());
                                }
                                u32::from(first)
                            };
                            value.push(
                                char::from_u32(scalar)
                                    .ok_or_else(|| "invalid Unicode scalar".to_string())?,
                            );
                            continue;
                        }
                        _ => return Err("invalid JSON string escape".to_string()),
                    }
                    self.cursor += 1;
                }
                Some(_) => return Err("JSON string contains an unescaped control".to_string()),
                None => return Err("unterminated JSON string".to_string()),
            }
        }
    }

    fn parse_hex_quad(&mut self) -> Result<u16, String> {
        let digits = self
            .bytes
            .get(self.cursor..self.cursor + 4)
            .ok_or_else(|| "short Unicode escape".to_string())?;
        let digits =
            core::str::from_utf8(digits).map_err(|_| "Unicode escape is not ASCII".to_string())?;
        let value = u16::from_str_radix(digits, 16)
            .map_err(|_| "Unicode escape is not hexadecimal".to_string())?;
        self.cursor += 4;
        Ok(value)
    }

    fn parse_number(&mut self) -> Result<String, String> {
        let start = self.cursor;
        if self.bytes.get(self.cursor) == Some(&b'-') {
            self.cursor += 1;
        }
        match self.bytes.get(self.cursor).copied() {
            Some(b'0') => self.cursor += 1,
            Some(b'1'..=b'9') => {
                self.cursor += 1;
                while self.bytes.get(self.cursor).is_some_and(u8::is_ascii_digit) {
                    self.cursor += 1;
                }
            }
            _ => return Err("JSON number has an invalid integer part".to_string()),
        }
        if self.bytes.get(self.cursor) == Some(&b'.') {
            self.cursor += 1;
            let fraction_start = self.cursor;
            while self.bytes.get(self.cursor).is_some_and(u8::is_ascii_digit) {
                self.cursor += 1;
            }
            if self.cursor == fraction_start {
                return Err("JSON number has an empty fraction".to_string());
            }
        }
        if matches!(self.bytes.get(self.cursor), Some(b'e' | b'E')) {
            self.cursor += 1;
            if matches!(self.bytes.get(self.cursor), Some(b'+' | b'-')) {
                self.cursor += 1;
            }
            let exponent_start = self.cursor;
            while self.bytes.get(self.cursor).is_some_and(u8::is_ascii_digit) {
                self.cursor += 1;
            }
            if self.cursor == exponent_start {
                return Err("JSON number has an empty exponent".to_string());
            }
        }
        core::str::from_utf8(&self.bytes[start..self.cursor])
            .map(str::to_string)
            .map_err(|_| "JSON number is not ASCII".to_string())
    }

    fn skip_whitespace(&mut self) {
        while matches!(
            self.bytes.get(self.cursor),
            Some(b' ' | b'\n' | b'\r' | b'\t')
        ) {
            self.cursor += 1;
        }
    }
}

fn row_object(value: &RowJson) -> Result<&BTreeMap<String, RowJson>, String> {
    match value {
        RowJson::Object(fields) => Ok(fields),
        _ => Err("expected a JSON object".to_string()),
    }
}

fn row_array(value: &RowJson) -> Result<&[RowJson], String> {
    match value {
        RowJson::Array(values) => Ok(values),
        _ => Err("expected a JSON array".to_string()),
    }
}

fn row_field<'a>(
    fields: &'a BTreeMap<String, RowJson>,
    field: &str,
) -> Result<&'a RowJson, String> {
    fields
        .get(field)
        .ok_or_else(|| format!("missing `{field}` field"))
}

fn row_string<'a>(value: &'a RowJson, field: &str) -> Result<&'a str, String> {
    match value {
        RowJson::String(value) => Ok(value),
        _ => Err(format!("`{field}` is not a string")),
    }
}

fn row_u64(value: &RowJson, field: &str) -> Result<u64, String> {
    match value {
        RowJson::Number(value) => value
            .parse()
            .map_err(|_| format!("`{field}` is not a canonical unsigned integer")),
        _ => Err(format!("`{field}` is not a number")),
    }
}

fn row_u32(value: &RowJson, field: &str) -> Result<u32, String> {
    u32::try_from(row_u64(value, field)?).map_err(|_| format!("`{field}` does not fit in u32"))
}

fn row_bool(value: &RowJson, field: &str) -> Result<bool, String> {
    match value {
        RowJson::Bool(value) => Ok(*value),
        _ => Err(format!("`{field}` is not a Boolean")),
    }
}

fn row_lower_hex(value: &str, field: &str, max_bytes: usize) -> Result<Vec<u8>, String> {
    if value.len() > max_bytes.saturating_mul(2) {
        return Err(format!(
            "`{field}` exceeds the {max_bytes}-byte decoded limit"
        ));
    }
    if !value.len().is_multiple_of(2)
        || !value
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
    {
        return Err(format!(
            "`{field}` is not even-length lowercase hexadecimal"
        ));
    }
    let (pairs, remainder) = value.as_bytes().as_chunks::<2>();
    debug_assert!(remainder.is_empty());
    pairs
        .iter()
        .map(|pair| {
            let digits = core::str::from_utf8(pair)
                .map_err(|_| format!("`{field}` is not ASCII hexadecimal"))?;
            u8::from_str_radix(digits, 16).map_err(|_| format!("`{field}` is not hexadecimal"))
        })
        .collect()
}

fn row_fixed_lower_hex<const N: usize>(value: &str, field: &str) -> Result<[u8; N], String> {
    if value.len() != N.saturating_mul(2) {
        return Err(format!("`{field}` must encode exactly {N} bytes"));
    }
    let bytes = row_lower_hex(value, field, N)?;
    bytes
        .try_into()
        .map_err(|_| format!("`{field}` must encode exactly {N} bytes"))
}

fn row_hash(value: &RowJson, field: &str) -> Result<ContentHash, String> {
    let encoded = row_string(value, field)?;
    Ok(ContentHash(row_fixed_lower_hex::<32>(encoded, field)?))
}

fn row_optional_hash(value: &RowJson, field: &str) -> Result<Option<ContentHash>, String> {
    match value {
        RowJson::Null => Ok(None),
        _ => row_hash(value, field).map(Some),
    }
}

fn row_operation(value: &RowJson, field: &str) -> Result<Option<u8>, String> {
    match value {
        RowJson::Null => Ok(None),
        RowJson::String(operation) => match operation.as_str() {
            "add" => Ok(Some(interval_op_tag(IntervalOp::Add))),
            "mul" => Ok(Some(interval_op_tag(IntervalOp::Mul))),
            "hull" => Ok(Some(interval_op_tag(IntervalOp::Hull))),
            _ => Err(format!("`{field}` names an unknown interval operation")),
        },
        _ => Err(format!("`{field}` is neither null nor a string")),
    }
}

fn row_waiver(value: &RowJson, field: &str) -> Result<Option<Waiver>, String> {
    if matches!(value, RowJson::Null) {
        return Ok(None);
    }
    let object = row_object(value).map_err(|_| format!("`{field}` is not an object or null"))?;
    let waiver = Waiver {
        id: row_string(row_field(object, "id")?, "id")?.to_string(),
        signer: row_string(row_field(object, "signer")?, "signer")?.to_string(),
        reason: row_string(row_field(object, "reason")?, "reason")?.to_string(),
    };
    if let Some((invalid_field, reason)) = waiver_annotation_reason(&waiver) {
        return Err(format!(
            "`{field}.{invalid_field}` violates the waiver grammar ({reason})"
        ));
    }
    Ok(Some(waiver))
}

struct PersistedGrant {
    signing_payload: Vec<u8>,
    signature: Vec<u8>,
    retained_bytes: usize,
}

#[allow(clippy::too_many_lines)] // Complete persisted grant reconstruction is one invariant.
fn row_grant(
    value: &RowJson,
    field: &str,
    operation: Option<u8>,
    waiver: Option<&Waiver>,
) -> Result<Option<PersistedGrant>, String> {
    if matches!(value, RowJson::Null) {
        return Ok(None);
    }
    let waiver = waiver.ok_or_else(|| format!("`{field}` lacks its waiver annotation"))?;
    let object = row_object(value).map_err(|_| format!("`{field}` is not an object or null"))?;
    let expected_version = if operation.is_some() { 3 } else { 4 };
    let payload_version = row_u64(row_field(object, "payload_version")?, "payload_version")?;
    if payload_version != expected_version {
        return Err(format!(
            "`{field}.payload_version` is {payload_version}, expected {expected_version}"
        ));
    }
    if !row_bool(row_field(object, "authorized")?, "authorized")? {
        return Err(format!("`{field}.authorized` is not true"));
    }
    let key_id = row_string(row_field(object, "key_id")?, "key_id")?;
    let scope = row_string(row_field(object, "scope")?, "scope")?;
    let node_name = row_string(row_field(object, "node_name")?, "node_name")?;
    let claimed_color = row_lower_hex(
        row_string(row_field(object, "claimed_color_hex")?, "claimed_color_hex")?,
        "claimed_color_hex",
        MAX_CLAIMED_COLOR_BYTES,
    )?;
    if claimed_color.is_empty() || claimed_color.len() > MAX_CLAIMED_COLOR_BYTES {
        return Err(format!(
            "`{field}.claimed_color_hex` has an invalid decoded length"
        ));
    }
    let parent_values = row_array(row_field(object, "parent_hashes")?)?;
    if parent_values.len() > MAX_COLOR_PARENTS {
        return Err(format!("`{field}.parent_hashes` exceeds the parent limit"));
    }
    let parent_hashes = parent_values
        .iter()
        .map(|value| row_hash(value, "parent_hashes[]"))
        .collect::<Result<Vec<_>, _>>()?;
    let expires_day = row_u32(row_field(object, "expires_day")?, "expires_day")?;
    let signing_payload = row_lower_hex(
        row_string(
            row_field(object, "signing_payload_hex")?,
            "signing_payload_hex",
        )?,
        "signing_payload_hex",
        MAX_WAIVER_CLOSURE_BYTES,
    )?;
    let signature = row_lower_hex(
        row_string(row_field(object, "signature_hex")?, "signature_hex")?,
        "signature_hex",
        MAX_WAIVER_SIGNATURE_BYTES,
    )?;
    if signature.is_empty() || signature.len() > MAX_WAIVER_SIGNATURE_BYTES {
        return Err(format!(
            "`{field}.signature_hex` has an invalid decoded length"
        ));
    }
    let retained_bytes = [
        waiver.id.len(),
        waiver.signer.len(),
        waiver.reason.len(),
        key_id.len(),
        scope.len(),
        node_name.len(),
        claimed_color.len(),
        signature.len(),
    ]
    .into_iter()
    .try_fold(256usize, usize::checked_add)
    .and_then(|bytes| {
        parent_hashes
            .len()
            .checked_mul(32)
            .and_then(|hash_bytes| bytes.checked_add(hash_bytes))
    })
    .ok_or_else(|| format!("`{field}` retained-byte count overflowed"))?;

    let mut reconstructed = vec![u8::try_from(expected_version).expect("version is 3 or 4")];
    push_field(&mut reconstructed, WAIVER_PAYLOAD_DOMAIN);
    reconstructed.push(operation.unwrap_or(0));
    for text in [key_id, scope, node_name] {
        push_field(&mut reconstructed, text.as_bytes());
    }
    push_field(&mut reconstructed, &claimed_color);
    for text in [
        waiver.id.as_str(),
        waiver.signer.as_str(),
        waiver.reason.as_str(),
    ] {
        push_field(&mut reconstructed, text.as_bytes());
    }
    push_len(&mut reconstructed, parent_hashes.len());
    for hash in parent_hashes {
        reconstructed.extend_from_slice(hash.as_bytes());
    }
    reconstructed.extend_from_slice(&expires_day.to_le_bytes());
    if signing_payload != reconstructed {
        return Err(format!(
            "`{field}.signing_payload_hex` does not reconstruct from persisted grant fields"
        ));
    }
    Ok(Some(PersistedGrant {
        signing_payload,
        signature,
        retained_bytes,
    }))
}

#[derive(Debug)]
struct PersistedDemotion {
    parent_index: usize,
    parent_id: u64,
    dataset: String,
    axis: String,
    value_bits: u64,
}

#[derive(Clone, Copy)]
struct RowErrorContext {
    row: usize,
    node: Option<u64>,
}

impl RowErrorContext {
    fn malformed(self, why: impl Into<String>) -> ColorRowVerificationError {
        ColorRowVerificationError::Malformed {
            row: self.row,
            node: self.node,
            why: why.into(),
        }
    }

    fn decode<T>(self, result: Result<T, String>) -> Result<T, ColorRowVerificationError> {
        result.map_err(|why| self.malformed(why))
    }
}

/// Independently parse, reconstruct, and rehash a canonical persisted color
/// row stream. This path consumes only JSON rows: it never consults a
/// [`ColorGraph`], [`ColorNode`], or in-memory [`Color`]. Parent hashes come
/// from earlier accepted `color-write` rows, while immediately preceding
/// demotion rows are folded into the named node's exact v9 preimage.
///
/// Schema-v7 color rows, schema-v1 demotion rows, node-hash v9, and the current
/// color algebra are accepted literally. Every other version is refused; no
/// legacy display JSON is guessed or silently migrated. Exact lowercase
/// `color_canonical_hex`, `origin_canonical_hex`, IEEE-754 demotion bits,
/// grant payloads/signatures, policy fingerprints, and authority closure are
/// the persisted preimage inputs.
///
/// # Errors
/// [`ColorRowVerificationError`] names the first malformed, unsupported, or
/// hash-divergent row.
#[allow(clippy::too_many_lines)] // The field order mirrors the security-critical v9 preimage.
pub fn verify_color_row_stream(
    rows: &[String],
) -> Result<ColorRowVerification, ColorRowVerificationError> {
    let mut hashes = Vec::<ContentHash>::new();
    let mut pending_demotions = Vec::<PersistedDemotion>::new();
    let mut demotion_count = 0usize;

    for (row_index, row) in rows.iter().enumerate() {
        let initial = RowErrorContext {
            row: row_index,
            node: None,
        };
        let parsed = initial.decode(RowJsonParser::new(row).parse())?;
        let object = initial.decode(row_object(&parsed))?;
        let event = initial.decode(
            row_field(object, "event")
                .and_then(|value| row_string(value, "event"))
                .map(str::to_string),
        )?;
        let node =
            initial.decode(row_field(object, "node").and_then(|value| row_u64(value, "node")))?;
        let context = RowErrorContext {
            row: row_index,
            node: Some(node),
        };
        let expected_node = u64::try_from(hashes.len())
            .map_err(|_| context.malformed("verified node count no longer fits in u64"))?;

        match event.as_str() {
            "demotion" => {
                let schema = context.decode(
                    row_field(object, "schema_version")
                        .and_then(|value| row_u64(value, "schema_version")),
                )?;
                if schema != u64::from(COLOR_DEMOTION_ROW_SCHEMA_VERSION) {
                    return Err(ColorRowVerificationError::UnsupportedSchemaVersion {
                        row: row_index,
                        node,
                        event: "demotion",
                        found: schema,
                        expected: u64::from(COLOR_DEMOTION_ROW_SCHEMA_VERSION),
                    });
                }
                if node != expected_node {
                    return Err(context.malformed(format!(
                        "demotion must precede the next node {expected_node}, not node {node}"
                    )));
                }
                let parent_index_u64 = context.decode(
                    row_field(object, "parent_index")
                        .and_then(|value| row_u64(value, "parent_index")),
                )?;
                let parent_index = usize::try_from(parent_index_u64)
                    .map_err(|_| context.malformed("parent_index does not fit in usize"))?;
                if pending_demotions
                    .last()
                    .is_some_and(|previous| previous.parent_index >= parent_index)
                {
                    return Err(
                        context.malformed("demotion parent indexes are not strictly increasing")
                    );
                }
                let parent_id = context.decode(
                    row_field(object, "parent").and_then(|value| row_u64(value, "parent")),
                )?;
                let dataset = context.decode(
                    row_field(object, "dataset")
                        .and_then(|value| row_string(value, "dataset"))
                        .map(str::to_string),
                )?;
                let axis = context.decode(
                    row_field(object, "axis")
                        .and_then(|value| row_string(value, "axis"))
                        .map(str::to_string),
                )?;
                let value_bits = context.decode(
                    row_field(object, "value_bits")
                        .and_then(|value| row_string(value, "value_bits"))
                        .and_then(|value| row_fixed_lower_hex::<8>(value, "value_bits")),
                )?;
                match context.decode(row_field(object, "value"))? {
                    RowJson::Number(_) | RowJson::String(_) => {}
                    _ => {
                        return Err(context.malformed(
                            "demotion display value is neither a number nor a tagged string",
                        ));
                    }
                }
                pending_demotions.push(PersistedDemotion {
                    parent_index,
                    parent_id,
                    dataset,
                    axis,
                    value_bits: u64::from_be_bytes(value_bits),
                });
                demotion_count += 1;
            }
            "color-write" => {
                let schema = context.decode(
                    row_field(object, "schema_version")
                        .and_then(|value| row_u64(value, "schema_version")),
                )?;
                if schema != u64::from(COLOR_WRITE_ROW_SCHEMA_VERSION) {
                    return Err(ColorRowVerificationError::UnsupportedSchemaVersion {
                        row: row_index,
                        node,
                        event: "color-write",
                        found: schema,
                        expected: u64::from(COLOR_WRITE_ROW_SCHEMA_VERSION),
                    });
                }
                if node != expected_node {
                    return Err(context.malformed(format!(
                        "node ids must be dense in stream order; expected {expected_node}"
                    )));
                }
                let hash_version = context.decode(
                    row_field(object, "node_hash_version")
                        .and_then(|value| row_u64(value, "node_hash_version")),
                )?;
                if hash_version != u64::from(COLOR_NODE_HASH_ENCODING_VERSION) {
                    return Err(ColorRowVerificationError::UnsupportedNodeHashVersion {
                        row: row_index,
                        node,
                        found: hash_version,
                        expected: u64::from(COLOR_NODE_HASH_ENCODING_VERSION),
                    });
                }
                let algebra_version = context.decode(
                    row_field(object, "color_algebra_version")
                        .and_then(|value| row_u64(value, "color_algebra_version")),
                )?;
                if algebra_version != u64::from(COLOR_ALGEBRA_VERSION) {
                    return Err(ColorRowVerificationError::UnsupportedColorAlgebraVersion {
                        row: row_index,
                        node,
                        found: algebra_version,
                        expected: u64::from(COLOR_ALGEBRA_VERSION),
                    });
                }

                let operation = context.decode(
                    row_field(object, "operation")
                        .and_then(|value| row_operation(value, "operation")),
                )?;
                let name = context.decode(
                    row_field(object, "name")
                        .and_then(|value| row_string(value, "name"))
                        .map(str::to_string),
                )?;
                if let Some(reason) = identity_reason(&name) {
                    return Err(context.malformed(format!(
                        "node name violates the canonical identity grammar ({reason})"
                    )));
                }
                let color_name = context.decode(
                    row_field(object, "color").and_then(|value| row_string(value, "color")),
                )?;
                if !matches!(color_name, "verified" | "validated" | "estimated") {
                    return Err(context.malformed("color names an unknown evidence rank"));
                }
                context.decode(row_field(object, "payload").and_then(row_object))?;
                let color_bytes = context.decode(
                    row_field(object, "color_canonical_hex")
                        .and_then(|value| row_string(value, "color_canonical_hex"))
                        .and_then(|value| {
                            row_lower_hex(value, "color_canonical_hex", MAX_CLAIMED_COLOR_BYTES)
                        }),
                )?;
                if color_bytes.is_empty() || color_bytes.len() > MAX_CLAIMED_COLOR_BYTES {
                    return Err(
                        context.malformed("color_canonical_hex has an invalid decoded length")
                    );
                }
                let rank_tag = match color_name {
                    "verified" => 0,
                    "validated" => 1,
                    "estimated" => 2,
                    _ => unreachable!("color rank was validated above"),
                };
                if color_bytes.first() != Some(&(COLOR_ALGEBRA_VERSION as u8))
                    || color_bytes.get(1) != Some(&rank_tag)
                {
                    return Err(context
                        .malformed("canonical color header disagrees with the row algebra/rank"));
                }
                let parent_values =
                    context.decode(row_field(object, "parents").and_then(row_array))?;
                if parent_values.len() > MAX_COLOR_PARENTS {
                    return Err(context.malformed("parent list exceeds the verifier limit"));
                }
                let parents = parent_values
                    .iter()
                    .map(|value| row_u64(value, "parents[]"))
                    .collect::<Result<Vec<_>, _>>()
                    .map_err(|why| context.malformed(why))?;
                if parents.is_empty() != operation.is_none() {
                    return Err(context.malformed(
                        "source/derived parent shape disagrees with operation presence",
                    ));
                }
                for &parent in &parents {
                    let parent_index = usize::try_from(parent)
                        .map_err(|_| context.malformed("parent id does not fit in usize"))?;
                    if parent_index >= hashes.len() {
                        return Err(context.malformed(format!(
                            "parent {parent} does not name an earlier accepted row"
                        )));
                    }
                }
                for demotion in &pending_demotions {
                    if parents.get(demotion.parent_index) != Some(&demotion.parent_id) {
                        return Err(context.malformed(
                            "demotion parent position/id disagrees with the color-write row",
                        ));
                    }
                }

                let origin_bytes =
                    match context.decode(row_field(object, "origin_canonical_hex"))? {
                        RowJson::Null => None,
                        value => {
                            let bytes = context.decode(
                                row_string(value, "origin_canonical_hex").and_then(|value| {
                                    row_lower_hex(
                                        value,
                                        "origin_canonical_hex",
                                        MAX_CLAIMED_COLOR_BYTES,
                                    )
                                }),
                            )?;
                            if bytes.is_empty() {
                                return Err(context
                                    .malformed("origin_canonical_hex is empty instead of null"));
                            }
                            Some(bytes)
                        }
                    };
                let origin_present = match context.decode(row_field(object, "origin"))? {
                    RowJson::Null => false,
                    value => {
                        context.decode(row_object(value))?;
                        true
                    }
                };
                if origin_present != origin_bytes.is_some() {
                    return Err(context.malformed(
                        "origin display object disagrees with exact origin-byte presence",
                    ));
                }
                let origin_policy = context.decode(
                    row_field(object, "origin_policy_fingerprint")
                        .and_then(|value| row_optional_hash(value, "origin_policy_fingerprint")),
                )?;
                if origin_bytes.is_some() != origin_policy.is_some() {
                    return Err(context.malformed(
                        "exact source origin and admitting policy must appear together",
                    ));
                }
                if origin_bytes.is_some() && !parents.is_empty() {
                    return Err(
                        context.malformed("derived row carries source-only origin metadata")
                    );
                }

                let dependency_values =
                    context.decode(row_field(object, "waiver_dependencies").and_then(row_array))?;
                if dependency_values.len() > MAX_WAIVER_DEPENDENCIES {
                    return Err(
                        context.malformed("waiver dependency list exceeds the verifier limit")
                    );
                }
                let mut dependencies = Vec::with_capacity(dependency_values.len());
                let mut previous_authority = None;
                let mut closure_bytes = 0usize;
                for dependency_value in dependency_values {
                    let dependency = context.decode(row_object(dependency_value))?;
                    let authorizing_node = context.decode(
                        row_field(dependency, "authorizing_node")
                            .and_then(|value| row_u64(value, "authorizing_node")),
                    )?;
                    if authorizing_node >= node
                        || previous_authority.is_some_and(|previous| previous >= authorizing_node)
                    {
                        return Err(context.malformed(
                            "waiver dependencies are not unique earlier nodes in ascending order",
                        ));
                    }
                    previous_authority = Some(authorizing_node);
                    let dependency_operation = context.decode(
                        row_field(dependency, "operation")
                            .and_then(|value| row_operation(value, "dependency.operation")),
                    )?;
                    let policy = context.decode(
                        row_field(dependency, "policy_fingerprint")
                            .and_then(|value| row_hash(value, "policy_fingerprint")),
                    )?;
                    let admission_day = context.decode(
                        row_field(dependency, "admission_day")
                            .and_then(|value| row_u32(value, "admission_day")),
                    )?;
                    let dependency_waiver = context.decode(
                        row_field(dependency, "waiver")
                            .and_then(|value| row_waiver(value, "dependency.waiver")),
                    )?;
                    let grant = context
                        .decode(row_field(dependency, "grant").and_then(|value| {
                            row_grant(
                                value,
                                "dependency.grant",
                                dependency_operation,
                                dependency_waiver.as_ref(),
                            )
                        }))?
                        .ok_or_else(|| {
                            context.malformed("waiver dependency carries a null grant")
                        })?;
                    closure_bytes = closure_bytes
                        .checked_add(grant.retained_bytes)
                        .ok_or_else(|| context.malformed("waiver closure byte count overflowed"))?;
                    if closure_bytes > MAX_WAIVER_CLOSURE_BYTES {
                        return Err(
                            context.malformed("waiver closure exceeds the verifier byte limit")
                        );
                    }
                    dependencies.push((
                        authorizing_node,
                        dependency_operation,
                        grant,
                        policy,
                        admission_day,
                    ));
                }

                let waiver = context.decode(
                    row_field(object, "waiver").and_then(|value| row_waiver(value, "waiver")),
                )?;
                let grant = context.decode(
                    row_field(object, "grant")
                        .and_then(|value| row_grant(value, "grant", operation, waiver.as_ref())),
                )?;
                if let Some(grant) = &grant {
                    closure_bytes = closure_bytes
                        .checked_add(grant.retained_bytes)
                        .ok_or_else(|| context.malformed("waiver closure byte count overflowed"))?;
                    if closure_bytes > MAX_WAIVER_CLOSURE_BYTES {
                        return Err(
                            context.malformed("waiver closure exceeds the verifier byte limit")
                        );
                    }
                }
                let waiver_policy = context.decode(
                    row_field(object, "waiver_policy_fingerprint")
                        .and_then(|value| row_optional_hash(value, "waiver_policy_fingerprint")),
                )?;
                let waiver_day = match context.decode(row_field(object, "waiver_admission_day"))? {
                    RowJson::Null => None,
                    value => Some(context.decode(row_u32(value, "waiver_admission_day"))?),
                };
                if grant.is_some() != waiver_policy.is_some()
                    || grant.is_some() != waiver_day.is_some()
                {
                    return Err(context.malformed(
                        "direct grant, policy fingerprint, and admission day must appear together",
                    ));
                }
                if parents.is_empty() && !dependencies.is_empty() {
                    return Err(context
                        .malformed("source row carries an impossible transitive waiver closure"));
                }

                let mut preimage = vec![COLOR_NODE_HASH_ENCODING_VERSION];
                push_field(&mut preimage, COLOR_NODE_HASH_DOMAIN);
                if let Some(operation) = operation {
                    preimage.push(1);
                    preimage.push(operation);
                } else {
                    preimage.push(0);
                }
                push_field(&mut preimage, name.as_bytes());
                push_field(&mut preimage, &color_bytes);
                push_len(&mut preimage, parents.len());
                for parent in &parents {
                    let index = usize::try_from(*parent)
                        .map_err(|_| context.malformed("parent id does not fit in usize"))?;
                    push_field(&mut preimage, hashes[index].as_bytes());
                }
                push_len(&mut preimage, pending_demotions.len());
                for demotion in &pending_demotions {
                    push_len(&mut preimage, demotion.parent_index);
                    preimage.extend_from_slice(&demotion.parent_id.to_le_bytes());
                    push_field(&mut preimage, demotion.dataset.as_bytes());
                    push_field(&mut preimage, demotion.axis.as_bytes());
                    preimage.extend_from_slice(&demotion.value_bits.to_le_bytes());
                }
                if let Some(origin_bytes) = &origin_bytes {
                    preimage.push(1);
                    preimage.extend_from_slice(origin_bytes);
                } else {
                    preimage.push(0);
                }
                if let Some(policy) = origin_policy {
                    preimage.push(1);
                    preimage.extend_from_slice(policy.as_bytes());
                } else {
                    preimage.push(0);
                }
                push_len(&mut preimage, dependencies.len());
                for (authorizing_node, operation, grant, policy, admission_day) in dependencies {
                    preimage.extend_from_slice(&authorizing_node.to_le_bytes());
                    if let Some(operation) = operation {
                        preimage.push(1);
                        preimage.push(operation);
                    } else {
                        preimage.push(0);
                    }
                    push_field(&mut preimage, &grant.signing_payload);
                    push_field(&mut preimage, &grant.signature);
                    preimage.extend_from_slice(policy.as_bytes());
                    preimage.extend_from_slice(&admission_day.to_le_bytes());
                }
                if let Some(waiver) = &waiver {
                    preimage.push(1);
                    push_field(&mut preimage, waiver.id.as_bytes());
                    push_field(&mut preimage, waiver.signer.as_bytes());
                    push_field(&mut preimage, waiver.reason.as_bytes());
                } else {
                    preimage.push(0);
                }
                if let Some(grant) = grant {
                    preimage.push(1);
                    push_field(&mut preimage, &grant.signing_payload);
                    push_field(&mut preimage, &grant.signature);
                } else {
                    preimage.push(0);
                }
                if let Some(policy) = waiver_policy {
                    preimage.push(1);
                    preimage.extend_from_slice(policy.as_bytes());
                } else {
                    preimage.push(0);
                }
                if let Some(day) = waiver_day {
                    preimage.push(1);
                    preimage.extend_from_slice(&day.to_le_bytes());
                } else {
                    preimage.push(0);
                }

                let reconstructed = hash_bytes(&preimage);
                let stored = context
                    .decode(row_field(object, "hash").and_then(|value| row_hash(value, "hash")))?;
                if reconstructed != stored {
                    return Err(ColorRowVerificationError::HashMismatch {
                        row: row_index,
                        node,
                        stored,
                        reconstructed,
                    });
                }
                hashes.push(stored);
                pending_demotions.clear();
            }
            _ => {
                return Err(context.malformed(format!("unsupported color-row event {event:?}")));
            }
        }
    }

    if !pending_demotions.is_empty() {
        return Err(ColorRowVerificationError::Malformed {
            row: rows.len(),
            node: u64::try_from(hashes.len()).ok(),
            why: "stream ended before pending demotions' color-write row".to_string(),
        });
    }
    Ok(ColorRowVerification {
        node_count: hashes.len(),
        demotion_count,
        terminal_hash: hashes.last().copied(),
    })
}

/// Stable identity of the ledger color-admission policy (bead 6pf9): binds
/// the row schema and algebra versions the authority mints receipts under.
#[must_use]
pub fn color_admission_policy_fingerprint() -> ContentHash {
    hash_bytes(
        format!(
            "fs-ledger/color-admission-policy/v1/row-schema={COLOR_WRITE_ROW_SCHEMA_VERSION}/algebra={COLOR_ALGEBRA_VERSION}"
        )
        .as_bytes(),
    )
}

/// Why the admission authority refused to mint a receipt for a node.
#[derive(Debug, Clone, PartialEq)]
pub enum ColorAdmissionRefusal {
    /// No node with this id exists in the graph (unresolved evidence).
    UnknownNode {
        /// The requested node id.
        node: u64,
    },
    /// The node's scientific claim depends on an authenticated waiver,
    /// directly or transitively. Waived evidence never converts to admitted
    /// scientific evidence.
    WaiverTainted {
        /// The refused node id.
        node: u64,
    },
    /// Only positive ranks (Verified/Validated) receive admission receipts.
    NotPositive {
        /// The refused node id.
        node: u64,
        /// The node's declared rank.
        rank: ColorRank,
    },
    /// The node failed the structural replay audit: its stored state does
    /// not rederive from its inputs.
    ReplayDivergence(ColorReplayError),
    /// The node's Validated regime excludes the current execution state.
    /// Out-of-regime evidence must be re-derived (and demoted) before any
    /// admission decision.
    OutOfRegime {
        /// The refused node id.
        node: u64,
        /// The demotion the regime check derived.
        demotion: Demotion,
    },
}

impl core::fmt::Display for ColorAdmissionRefusal {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::UnknownNode { node } => write!(f, "admission refused: node {node} not found"),
            Self::WaiverTainted { node } => write!(
                f,
                "admission refused: node {node} depends on an authenticated waiver"
            ),
            Self::NotPositive { node, rank } => write!(
                f,
                "admission refused: node {node} declares {rank:?}, not positive evidence"
            ),
            Self::ReplayDivergence(error) => write!(f, "admission refused: {error}"),
            Self::OutOfRegime { node, .. } => write!(
                f,
                "admission refused: node {node} regime excludes the current state"
            ),
        }
    }
}

impl std::error::Error for ColorAdmissionRefusal {}

impl ColorGraph {
    /// Mint an admission receipt for one node (bead 6pf9). Receipts convert
    /// declared colors into [`fs_evidence::AdmittedColor`] through
    /// [`LedgerColorAdmissionVerifier`]; minting refuses waiver-tainted,
    /// non-positive, and replay-divergent nodes, so waived, unresolved, and
    /// tampered evidence never acquires scientific admission. Regime
    /// awareness lives in [`Self::admission_receipt_in_regime`].
    ///
    /// # Errors
    /// [`ColorAdmissionRefusal`] naming the refusing gate.
    pub fn admission_receipt(
        &self,
        id: u64,
    ) -> Result<fs_evidence::AdmissionReceipt, ColorAdmissionRefusal> {
        let node = self
            .node(id)
            .ok_or(ColorAdmissionRefusal::UnknownNode { node: id })?;
        let color = node
            .scientific_color()
            .ok_or(ColorAdmissionRefusal::WaiverTainted { node: id })?;
        let rank = color.rank();
        if rank == ColorRank::Estimated {
            return Err(ColorAdmissionRefusal::NotPositive { node: id, rank });
        }
        let position =
            usize::try_from(id).map_err(|_| ColorAdmissionRefusal::UnknownNode { node: id })?;
        self.verify_replay_node(position, node)
            .map_err(ColorAdmissionRefusal::ReplayDivergence)?;
        Ok(fs_evidence::AdmissionReceipt::from_parts(
            node.hash(),
            COLOR_WRITE_ROW_SCHEMA_VERSION,
            COLOR_ALGEBRA_VERSION,
            color_admission_policy_fingerprint(),
        ))
    }

    /// [`Self::admission_receipt`] with a regime gate: a Validated node whose
    /// regime excludes the CURRENT execution state is refused with the exact
    /// demotion the regime check derived. Regime exit demotes structurally —
    /// it never converts to admitted evidence at the stale rank.
    ///
    /// # Errors
    /// [`ColorAdmissionRefusal`] naming the refusing gate.
    pub fn admission_receipt_in_regime(
        &self,
        id: u64,
        state: &BTreeMap<String, f64>,
    ) -> Result<fs_evidence::AdmissionReceipt, ColorAdmissionRefusal> {
        let node = self
            .node(id)
            .ok_or(ColorAdmissionRefusal::UnknownNode { node: id })?;
        if let Some(color) = node.scientific_color()
            && let Some(demotion) = regime_demotion(color, state)
        {
            return Err(ColorAdmissionRefusal::OutOfRegime { node: id, demotion });
        }
        self.admission_receipt(id)
    }
}

/// The ledger-side admission oracle (bead 6pf9): authenticates a
/// (candidate, receipt) pair by re-deriving it from the graph's replay-
/// audited node state. Acceptance requires the receipt's node hash to name
/// a live node whose provenance hash re-derives, whose scientific color is
/// bit-exactly the candidate (canonical bytes, not display JSON), whose
/// receipt versions match this build, and whose policy fingerprint is this
/// authority's own.
#[derive(Debug, Clone, Copy)]
pub struct LedgerColorAdmissionVerifier<'g> {
    graph: &'g ColorGraph,
}

impl<'g> LedgerColorAdmissionVerifier<'g> {
    /// Wrap the graph that will re-derive admissions.
    #[must_use]
    pub fn new(graph: &'g ColorGraph) -> Self {
        LedgerColorAdmissionVerifier { graph }
    }
}

impl fs_evidence::AdmissionVerifier for LedgerColorAdmissionVerifier<'_> {
    fn verify(
        &self,
        candidate: &Color,
        receipt: &fs_evidence::AdmissionReceipt,
    ) -> fs_evidence::AdmissionDecision {
        let policy = color_admission_policy_fingerprint();
        if receipt.row_schema_version() != COLOR_WRITE_ROW_SCHEMA_VERSION
            || receipt.color_algebra_version() != COLOR_ALGEBRA_VERSION
            || receipt.policy_fingerprint() != policy
        {
            return fs_evidence::AdmissionDecision::reject(policy);
        }
        let Some((position, node)) = self
            .graph
            .nodes()
            .iter()
            .enumerate()
            .find(|(_, node)| node.hash() == receipt.node_hash())
        else {
            return fs_evidence::AdmissionDecision::reject(policy);
        };
        let Some(scientific) = node.scientific_color() else {
            return fs_evidence::AdmissionDecision::reject(policy);
        };
        if scientific.canonical_bytes() != candidate.canonical_bytes() {
            return fs_evidence::AdmissionDecision::reject(policy);
        }
        if self.graph.verify_replay_node(position, node).is_err() {
            return fs_evidence::AdmissionDecision::reject(policy);
        }
        fs_evidence::AdmissionDecision::accept(policy)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    struct AllowFixtureSource;

    impl SourceOriginVerifier for AllowFixtureSource {
        fn verify(&self, request: &SourceOriginRequest<'_>) -> PolicyDecision {
            let policy = hash_bytes(b"fs-ledger/internal-fixture/source-policy/v1");
            let accepted = matches!(
                (
                    request.node_name(),
                    request.claimed_color(),
                    request.origin()
                ),
                (
                    "anchored",
                    Color::Validated { .. },
                    SourceOrigin::Anchoring { .. }
                ) | (
                    "certified",
                    Color::Verified { .. },
                    SourceOrigin::Certificate { .. }
                )
            );
            if accepted {
                PolicyDecision::accept(policy)
            } else {
                PolicyDecision::reject(policy)
            }
        }
    }

    fn strict_row_string<'a>(row: &'a str, key: &str) -> Option<&'a str> {
        let marker = format!("\"{key}\":\"");
        let value = row.get(row.find(&marker)? + marker.len()..)?;
        value.get(..value.find('"')?)
    }

    fn strict_lower_hex(value: &str) -> Option<Vec<u8>> {
        if !value.len().is_multiple_of(2)
            || !value
                .bytes()
                .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
        {
            return None;
        }
        let (pairs, remainder) = value.as_bytes().as_chunks::<2>();
        debug_assert!(remainder.is_empty());
        pairs
            .iter()
            .map(|pair| {
                core::str::from_utf8(pair)
                    .ok()
                    .and_then(|digits| u8::from_str_radix(digits, 16).ok())
            })
            .collect()
    }

    #[test]
    #[allow(clippy::too_many_lines)] // One canonical schema-v7 reconstruction and tamper story.
    fn schema_v7_rows_reconstruct_exact_node_hash_payloads() {
        let mut graph = ColorGraph::new();
        let estimated_base = graph
            .source(
                "estimated-base",
                Color::Estimated {
                    estimator: "rom-v1".to_string(),
                    dispersion: 0.1,
                },
            )
            .expect("base estimated source");
        let estimated = graph
            .source(
                "estimated",
                Color::Estimated {
                    estimator: "rom-v1".to_string(),
                    dispersion: 0.1_f64.next_up(),
                },
            )
            .expect("estimated source");
        let certified_color = Color::Verified {
            lo: 1.0_f64.next_up(),
            hi: 2.0_f64.next_down(),
        };
        graph
            .source_with_origin(
                "certified",
                &certified_color,
                SourceOrigin::Certificate {
                    producer: "fixture-certifier".to_string(),
                    certificate_hash: hash_bytes(b"certificate artifact"),
                    certificate: NumericalCertificate::enclosure(
                        1.0_f64.next_up(),
                        2.0_f64.next_down(),
                    ),
                },
                &AllowFixtureSource,
            )
            .expect("certified source");
        let regime = ValidityDomain::unconstrained().with(
            "re",
            1_000.0_f64.next_up(),
            100_000.0_f64.next_down(),
        );
        let anchored_color = Color::Validated {
            regime: regime.clone(),
            dataset: "campaign-a".to_string(),
        };
        graph
            .source_with_origin(
                "anchored",
                &anchored_color,
                SourceOrigin::Anchoring {
                    dataset_id: "campaign-a".to_string(),
                    content_hash: hash_bytes(b"anchoring artifact"),
                    regime,
                },
                &AllowFixtureSource,
            )
            .expect("anchored source");
        graph
            .derive(
                "derived",
                &[estimated],
                IntervalOp::Hull,
                None,
                &BTreeMap::new(),
                None,
            )
            .expect("ordinary derivation");

        let base_node = graph.node(estimated_base).expect("base node");
        let adjacent_node = graph.node(estimated).expect("adjacent node");
        assert_eq!(
            base_node.color.payload_json(),
            adjacent_node.color.payload_json(),
            "display JSON deliberately rounds adjacent floats"
        );
        assert_ne!(
            base_node.color.canonical_bytes(),
            adjacent_node.color.canonical_bytes(),
        );
        assert_ne!(base_node.hash, adjacent_node.hash);
        let base_row = graph
            .rows()
            .iter()
            .find(|row| row.contains("\"name\":\"estimated-base\""))
            .expect("base color-write row");
        let adjacent_row = graph
            .rows()
            .iter()
            .find(|row| row.contains("\"name\":\"estimated\""))
            .expect("adjacent color-write row");
        assert_ne!(
            strict_row_string(base_row, "color_canonical_hex"),
            strict_row_string(adjacent_row, "color_canonical_hex"),
        );

        assert!(strict_lower_hex("AA").is_none());
        assert!(strict_lower_hex("0").is_none());
        for node in graph.nodes() {
            let row = graph
                .rows()
                .iter()
                .find(|row| {
                    row.contains("\"event\":\"color-write\"")
                        && row.contains(&format!("\"name\":\"{}\"", node.name))
                })
                .expect("one color-write row per node");
            assert!(row.contains("\"schema_version\":7"));
            assert!(row.contains("\"node_hash_version\":9"));
            let color_bytes = strict_lower_hex(
                strict_row_string(row, "color_canonical_hex").expect("exact color bytes"),
            )
            .expect("strict lowercase color hex");
            assert_eq!(color_bytes, node.color.canonical_bytes());
            let origin_bytes = strict_row_string(row, "origin_canonical_hex")
                .map(|encoded| strict_lower_hex(encoded).expect("strict lowercase origin hex"));
            assert_eq!(
                origin_bytes,
                node.origin.as_ref().map(source_origin_canonical_bytes)
            );
            let metadata = NodeHashMetadata {
                operation: node.operation,
                demotions: &node.demotions,
                origin: node.origin.as_ref(),
                origin_policy_fingerprint: node.origin_policy_fingerprint,
                waiver: node.waiver.as_ref(),
                grant: node.grant.as_ref(),
                waiver_policy_fingerprint: node.waiver_policy_fingerprint,
                waiver_admission_day: node.waiver_admission_day,
                waiver_dependencies: &node.waiver_dependencies,
            };
            assert_eq!(
                graph.node_hash_from_canonical_payloads(
                    &node.name,
                    &color_bytes,
                    &node.parents,
                    &metadata,
                    origin_bytes.as_deref(),
                ),
                node.hash,
            );
            let mut tampered = color_bytes.clone();
            tampered[0] ^= 1;
            assert_ne!(
                graph.node_hash_from_canonical_payloads(
                    &node.name,
                    &tampered,
                    &node.parents,
                    &metadata,
                    origin_bytes.as_deref(),
                ),
                node.hash,
            );
            if let Some(mut tampered_origin) = origin_bytes.clone() {
                tampered_origin[0] ^= 1;
                assert_ne!(
                    graph.node_hash_from_canonical_payloads(
                        &node.name,
                        &color_bytes,
                        &node.parents,
                        &metadata,
                        Some(&tampered_origin),
                    ),
                    node.hash,
                );
            }
        }
    }

    #[test]
    fn estimated_source_cannot_reroot_a_reserved_derived_identity() {
        let mut graph = ColorGraph::new();
        let error = graph
            .source(
                "rerooted",
                Color::Estimated {
                    estimator: "derived:composed:deadbeef".to_string(),
                    dispersion: f64::INFINITY,
                },
            )
            .expect_err("derived identities require retained parent lineage");
        assert!(matches!(
            error,
            ColorWriteError::InvalidEstimatedSource {
                field: "estimator",
                why: "derived-identity-requires-lineage",
            }
        ));
    }

    #[test]
    fn replay_rejects_hash_bound_source_origin_tamper() {
        let regime = ValidityDomain::unconstrained().with("re", 1e3, 1e5);
        let color = Color::Validated {
            regime: regime.clone(),
            dataset: "campaign-a".to_string(),
        };
        let mut graph = ColorGraph::new();
        let id = graph
            .source_with_origin(
                "anchored",
                &color,
                SourceOrigin::Anchoring {
                    dataset_id: "campaign-a".to_string(),
                    content_hash: hash_bytes(b"original artifact"),
                    regime,
                },
                &AllowFixtureSource,
            )
            .expect("valid anchor");
        graph.verify_replay().expect("untampered graph");

        let SourceOrigin::Anchoring { content_hash, .. } = graph.nodes
            [usize::try_from(id).expect("small id")]
        .origin
        .as_mut()
        .expect("origin") else {
            panic!("expected anchoring origin");
        };
        *content_hash = hash_bytes(b"substituted artifact");
        let error = graph.verify_replay().expect_err("tamper must diverge");
        assert_eq!(error.node, id);
        assert!(error.why.contains("provenance hash"));
    }

    #[test]
    fn replay_rejects_hash_consistent_estimated_leaf_tampering_only_at_sources() {
        let clean_source = || {
            let mut graph = ColorGraph::new();
            graph
                .source(
                    "estimate",
                    Color::Estimated {
                        estimator: "rom-v1".to_string(),
                        dispersion: 0.1,
                    },
                )
                .expect("valid Estimated source");
            graph
        };

        let mut rerooted = clean_source();
        rerooted.nodes[0].color = Color::Estimated {
            estimator: "derived:v2:composed:6:rom-v1".to_string(),
            dispersion: 0.1,
        };
        rehash_node(&mut rerooted, 0);
        let error = rerooted
            .verify_replay()
            .expect_err("a hash-consistent derived identity still needs its lineage");
        assert!(error.why.contains("derived-identity-requires-lineage"));

        let mut annotated_source = clean_source();
        annotated_source.nodes[0].waiver = Some(Waiver {
            id: "human-note".to_string(),
            signer: "reviewer".to_string(),
            reason: "an annotation is not source authority".to_string(),
        });
        rehash_node(&mut annotated_source, 0);
        let error = annotated_source
            .verify_replay()
            .expect_err("a hash-consistent orphan source annotation must refuse");
        assert!(
            error
                .why
                .contains("orphan authority or human-waiver metadata")
        );

        let mut derived = clean_source();
        derived
            .derive(
                "annotated-derived",
                &[0],
                IntervalOp::Hull,
                None,
                &BTreeMap::new(),
                Some(Waiver {
                    id: "review-note".to_string(),
                    signer: "reviewer".to_string(),
                    reason: "human context on a real operation".to_string(),
                }),
            )
            .expect("ordinary derived annotations remain legal");
        derived
            .verify_replay()
            .expect("source-only metadata rules must not reject derived operation nodes");

        derived.nodes[1]
            .waiver
            .as_mut()
            .expect("ordinary annotation")
            .reason
            .push('\u{202e}');
        rehash_node(&mut derived, 1);
        let error = derived
            .verify_replay()
            .expect_err("hash-consistent hostile annotations must fail replay");
        assert!(error.why.contains("waiver annotation"));
        assert!(error.why.contains("control-character"));
    }

    #[test]
    fn waiver_closure_byte_budget_refuses_before_clone_amplification() {
        let grant = WaiverGrant {
            annotation: Waiver {
                id: "large-authority".to_string(),
                signer: "fixture-authority".to_string(),
                reason: "aggregate closure accounting fixture".to_string(),
            },
            key_id: "fixture-key".to_string(),
            scope: WAIVER_SCOPE_COLOR_UPGRADE.to_string(),
            node_name: "large-node".to_string(),
            claimed_color: vec![0; MAX_CLAIMED_COLOR_BYTES],
            parent_hashes: Vec::new(),
            expires_day: u32::MAX,
            signature: vec![1],
        };
        let mut retained = 0usize;
        let mut admitted = 0usize;
        loop {
            match add_waiver_closure_bytes(&mut retained, &grant) {
                Ok(()) => admitted += 1,
                Err(ColorWriteError::ResourceLimitExceeded {
                    resource: "waiver_closure_bytes",
                    limit,
                    actual,
                }) => {
                    assert_eq!(limit, MAX_WAIVER_CLOSURE_BYTES);
                    assert!(actual > limit);
                    break;
                }
                other => panic!("unexpected closure accounting result: {other:?}"),
            }
        }
        assert!(admitted > 0 && admitted < MAX_WAIVER_DEPENDENCIES);
        assert!(retained <= MAX_WAIVER_CLOSURE_BYTES);
    }

    #[test]
    fn disjoint_axis_preflight_precedes_union_cap_in_any_key_order() {
        let mut aggregate = BTreeMap::new();
        aggregate.insert("z-shared".to_string(), (0.0, 1.0));
        for index in 0..(MAX_VALIDITY_AXES - 1) {
            aggregate.insert(format!("bounded-axis-{index:04}"), (0.0, 1.0));
        }
        assert_eq!(aggregate.len(), MAX_VALIDITY_AXES);
        let regime = ValidityDomain::unconstrained()
            .with("a-new-axis", 0.0, 1.0)
            .with("z-shared", 2.0, 3.0);
        assert!(matches!(
            merge_fold_validity_axes(&mut aggregate, &regime),
            Ok(false)
        ));
        assert_eq!(aggregate.len(), MAX_VALIDITY_AXES);
        assert!(!aggregate.contains_key("a-new-axis"));
    }

    #[test]
    fn replay_rederives_the_owner_defined_bounded_demotion_identity() {
        let regime = ValidityDomain::unconstrained().with("re", 1e3, 1e5);
        let color = Color::Validated {
            regime: regime.clone(),
            dataset: "campaign-a".to_string(),
        };
        let mut graph = ColorGraph::new();
        let source = graph
            .source_with_origin(
                "anchored",
                &color,
                SourceOrigin::Anchoring {
                    dataset_id: "campaign-a".to_string(),
                    content_hash: hash_bytes(b"anchoring dataset"),
                    regime,
                },
                &AllowFixtureSource,
            )
            .expect("valid anchor");
        let state = BTreeMap::from([("re".to_string(), 5e2)]);
        let derived = graph
            .derive(
                "outside-regime",
                &[source],
                IntervalOp::Hull,
                None,
                &state,
                None,
            )
            .expect("derive with automatic demotion");
        assert!(matches!(
            graph.node(derived).map(|node| &node.color),
            Some(Color::Estimated { estimator, dispersion })
                if estimator == &demotion_estimator_identity("campaign-a", "re")
                    && dispersion.is_infinite()
        ));
        graph
            .verify_replay()
            .expect("replay must share fs-evidence's demotion identity grammar");
    }

    #[test]
    fn replay_rejects_hash_consistent_malformed_waived_color() {
        let color = Color::Verified {
            lo: f64::NAN,
            hi: 1.0,
        };
        let annotation = Waiver {
            id: "historical-waiver".to_string(),
            signer: "fixture-authority".to_string(),
            reason: "replay must not trust historical admission".to_string(),
        };
        let grant = WaiverGrant {
            annotation: annotation.clone(),
            key_id: "fixture-key".to_string(),
            scope: WAIVER_SCOPE_SOURCE_COLOR.to_string(),
            node_name: "historical-malformed".to_string(),
            claimed_color: color.canonical_bytes(),
            parent_hashes: Vec::new(),
            expires_day: u32::MAX,
            signature: vec![1],
        };
        let mut graph = ColorGraph::new();
        let id = graph.push_node(
            "historical-malformed",
            color,
            Vec::new(),
            NodeWriteMetadata {
                operation: None,
                demotions: Vec::new(),
                origin: None,
                origin_policy_fingerprint: None,
                waiver: Some(annotation),
                grant: Some(grant),
                waiver_policy_fingerprint: Some(hash_bytes(b"historical waiver policy")),
                waiver_admission_day: Some(1),
                waiver_dependencies: Vec::new(),
            },
        );

        let error = graph
            .verify_replay()
            .expect_err("replay must reject a hash-consistent malformed color");
        assert_eq!(error.node, id);
        assert!(error.why.contains("structurally invalid"));
        assert!(error.why.contains("bounds contain NaN"));
    }

    fn historical_waiver_dependency_graph() -> (ColorGraph, u64) {
        let color = Color::Verified { lo: 1.0, hi: 2.0 };
        let annotation = Waiver {
            id: "historical-waiver".to_string(),
            signer: "fixture-authority".to_string(),
            reason: "historical policy exception".to_string(),
        };
        let grant = WaiverGrant {
            annotation: annotation.clone(),
            key_id: "fixture-key".to_string(),
            scope: WAIVER_SCOPE_SOURCE_COLOR.to_string(),
            node_name: "waived-source".to_string(),
            claimed_color: color.canonical_bytes(),
            parent_hashes: Vec::new(),
            expires_day: u32::MAX,
            signature: vec![1],
        };
        let mut graph = ColorGraph::new();
        let source = graph.push_node(
            "waived-source",
            color.clone(),
            Vec::new(),
            NodeWriteMetadata {
                operation: None,
                demotions: Vec::new(),
                origin: None,
                origin_policy_fingerprint: None,
                waiver: Some(annotation),
                grant: Some(grant),
                waiver_policy_fingerprint: Some(hash_bytes(b"historical waiver policy")),
                waiver_admission_day: Some(1),
                waiver_dependencies: Vec::new(),
            },
        );
        let dependencies = graph
            .inherited_waiver_dependencies(&[source], None)
            .expect("complete historical authority");
        let child = graph.push_node(
            "ordinary-child",
            color,
            vec![source],
            NodeWriteMetadata {
                operation: Some(IntervalOp::Hull),
                demotions: Vec::new(),
                origin: None,
                origin_policy_fingerprint: None,
                waiver: None,
                grant: None,
                waiver_policy_fingerprint: None,
                waiver_admission_day: None,
                waiver_dependencies: dependencies,
            },
        );
        graph.verify_replay().expect("fixture graph replays");
        (graph, child)
    }

    fn rehash_node(graph: &mut ColorGraph, id: u64) -> ContentHash {
        let index = usize::try_from(id).expect("small fixture id");
        let old_hash = graph.nodes[index].hash;
        let new_hash = {
            let node = &graph.nodes[index];
            let metadata = NodeHashMetadata {
                operation: node.operation,
                demotions: &node.demotions,
                origin: node.origin.as_ref(),
                origin_policy_fingerprint: node.origin_policy_fingerprint,
                waiver: node.waiver.as_ref(),
                grant: node.grant.as_ref(),
                waiver_policy_fingerprint: node.waiver_policy_fingerprint,
                waiver_admission_day: node.waiver_admission_day,
                waiver_dependencies: &node.waiver_dependencies,
            };
            graph.node_hash(&node.name, &node.color, &node.parents, &metadata)
        };
        graph.nodes[index].hash = new_hash;
        old_hash
    }

    #[test]
    fn replay_rejects_hash_consistent_missing_admission_context() {
        let regime = ValidityDomain::unconstrained().with("re", 1e3, 1e5);
        let color = Color::Validated {
            regime: regime.clone(),
            dataset: "campaign-a".to_string(),
        };
        let mut source_graph = ColorGraph::new();
        let source = source_graph
            .source_with_origin(
                "anchored",
                &color,
                SourceOrigin::Anchoring {
                    dataset_id: "campaign-a".to_string(),
                    content_hash: hash_bytes(b"anchor artifact"),
                    regime,
                },
                &AllowFixtureSource,
            )
            .expect("valid source");
        let source_index = usize::try_from(source).expect("small source id");
        source_graph.nodes[source_index].origin_policy_fingerprint = None;
        rehash_node(&mut source_graph, source);
        let error = source_graph
            .verify_replay()
            .expect_err("typed source without policy must refuse");
        assert!(error.why.contains("lacks its policy"));

        let (mut waiver_graph, _) = historical_waiver_dependency_graph();
        waiver_graph.nodes[0].waiver_admission_day = None;
        rehash_node(&mut waiver_graph, 0);
        let error = waiver_graph
            .verify_replay()
            .expect_err("direct waiver without admission day must refuse");
        assert!(error.why.contains("source grant fields"));
    }

    #[test]
    fn replay_rejects_hash_consistent_waiver_dependency_tamper() {
        let (mut omitted, child) = historical_waiver_dependency_graph();
        omitted.nodes[child as usize].waiver_dependencies.clear();
        let old_hash = rehash_node(&mut omitted, child);
        assert_ne!(old_hash, omitted.nodes[child as usize].hash);
        let error = omitted
            .verify_replay()
            .expect_err("omitted dependency must refuse");
        assert!(error.why.contains("parent-derived union"));

        let (mut mutated, child) = historical_waiver_dependency_graph();
        mutated.nodes[child as usize].waiver_dependencies[0]
            .grant
            .annotation
            .reason
            .push_str(" edited");
        let old_hash = rehash_node(&mut mutated, child);
        assert_ne!(old_hash, mutated.nodes[child as usize].hash);
        let error = mutated
            .verify_replay()
            .expect_err("mutated dependency must refuse");
        assert!(error.why.contains("historical authorizing node"));

        let (mut policy_mutated, child) = historical_waiver_dependency_graph();
        let index = usize::try_from(child).expect("small child id");
        policy_mutated.nodes[index].waiver_dependencies[0].policy_fingerprint =
            hash_bytes(b"substituted waiver policy");
        rehash_node(&mut policy_mutated, child);
        let error = policy_mutated
            .verify_replay()
            .expect_err("dependency policy substitution must refuse");
        assert!(error.why.contains("historical authorizing node"));

        let (mut duplicated, child) = historical_waiver_dependency_graph();
        let duplicate = duplicated.nodes[child as usize].waiver_dependencies[0].clone();
        duplicated.nodes[child as usize]
            .waiver_dependencies
            .push(duplicate);
        let old_hash = rehash_node(&mut duplicated, child);
        assert_ne!(old_hash, duplicated.nodes[child as usize].hash);
        let error = duplicated
            .verify_replay()
            .expect_err("duplicate dependency must refuse");
        assert!(error.why.contains("duplicated"));

        let (mut self_referential, child) = historical_waiver_dependency_graph();
        self_referential.nodes[child as usize].waiver_dependencies[0].authorizing_node = child;
        let old_hash = rehash_node(&mut self_referential, child);
        assert_ne!(old_hash, self_referential.nodes[child as usize].hash);
        let error = self_referential
            .verify_replay()
            .expect_err("self/forward dependency must refuse");
        assert!(error.why.contains("self/forward"));
    }

    // ── Admission battery (bead 6pf9, stage S1) ────────────────────────────

    use fs_evidence::{
        AdmissionRejection, AdmissionVerifier as _, AdmittedColor, NoAdmissionVerifier,
        no_admission_policy,
    };

    /// One clean graph with a certified Verified source, an anchored
    /// Validated source, and an Estimated source.
    fn admission_fixture() -> (ColorGraph, u64, u64, u64) {
        let mut graph = ColorGraph::new();
        let certified_color = Color::Verified { lo: 1.0, hi: 2.0 };
        let certified = graph
            .source_with_origin(
                "certified",
                &certified_color,
                SourceOrigin::Certificate {
                    producer: "fixture-certifier".to_string(),
                    certificate_hash: hash_bytes(b"admission certificate"),
                    certificate: NumericalCertificate::enclosure(1.0, 2.0),
                },
                &AllowFixtureSource,
            )
            .expect("certified source");
        let regime = ValidityDomain::unconstrained().with("re", 1e3, 1e5);
        let anchored = graph
            .source_with_origin(
                "anchored",
                &Color::Validated {
                    regime: regime.clone(),
                    dataset: "campaign-a".to_string(),
                },
                SourceOrigin::Anchoring {
                    dataset_id: "campaign-a".to_string(),
                    content_hash: hash_bytes(b"admission anchor"),
                    regime,
                },
                &AllowFixtureSource,
            )
            .expect("anchored source");
        let estimated = graph
            .source(
                "estimated",
                Color::Estimated {
                    estimator: "rom-v1".to_string(),
                    dispersion: 0.25,
                },
            )
            .expect("estimated source");
        (graph, certified, anchored, estimated)
    }

    #[test]
    fn admission_mints_for_clean_positive_nodes_and_the_oracle_accepts() {
        let (graph, certified, anchored, _) = admission_fixture();
        let verifier = LedgerColorAdmissionVerifier::new(&graph);
        for id in [certified, anchored] {
            let receipt = graph.admission_receipt(id).expect("mint receipt");
            assert_eq!(
                receipt.policy_fingerprint(),
                color_admission_policy_fingerprint()
            );
            let declared = graph
                .node(id)
                .expect("node")
                .scientific_color()
                .expect("unwaived")
                .clone();
            let admitted = AdmittedColor::from_receipt(declared.clone(), receipt, &verifier)
                .expect("oracle admits its own receipt");
            assert_eq!(admitted.admitted_color(), &declared);
            assert!(admitted.rank() >= ColorRank::Validated);
        }
    }

    #[test]
    fn deny_all_default_refuses_even_a_genuine_receipt() {
        let (graph, certified, _, _) = admission_fixture();
        let receipt = graph.admission_receipt(certified).expect("mint receipt");
        let declared = graph
            .node(certified)
            .expect("node")
            .scientific_color()
            .expect("unwaived")
            .clone();
        let error = AdmittedColor::from_receipt(declared, receipt, &NoAdmissionVerifier)
            .expect_err("deny-all refuses genuine evidence");
        assert_eq!(
            error,
            AdmissionRejection::Refused {
                policy: no_admission_policy()
            }
        );
    }

    #[test]
    fn waived_estimated_and_unknown_nodes_never_mint_receipts() {
        let (waived_graph, child) = historical_waiver_dependency_graph();
        // Both the directly waived source and its transitively tainted child.
        assert_eq!(
            waived_graph.admission_receipt(0),
            Err(ColorAdmissionRefusal::WaiverTainted { node: 0 })
        );
        assert_eq!(
            waived_graph.admission_receipt(child),
            Err(ColorAdmissionRefusal::WaiverTainted { node: child })
        );

        let (graph, _, _, estimated) = admission_fixture();
        assert_eq!(
            graph.admission_receipt(estimated),
            Err(ColorAdmissionRefusal::NotPositive {
                node: estimated,
                rank: ColorRank::Estimated
            })
        );
        assert_eq!(
            graph.admission_receipt(999),
            Err(ColorAdmissionRefusal::UnknownNode { node: 999 })
        );
    }

    #[test]
    fn forged_receipts_and_substituted_candidates_are_refused() {
        let (graph, certified, anchored, _) = admission_fixture();
        let verifier = LedgerColorAdmissionVerifier::new(&graph);
        let policy = color_admission_policy_fingerprint();
        let declared = graph
            .node(certified)
            .expect("node")
            .scientific_color()
            .expect("unwaived")
            .clone();

        // Receipt naming a node hash the graph never produced.
        let forged = fs_evidence::AdmissionReceipt::from_parts(
            hash_bytes(b"no such node"),
            COLOR_WRITE_ROW_SCHEMA_VERSION,
            COLOR_ALGEBRA_VERSION,
            policy,
        );
        assert!(!verifier.verify(&declared, &forged).accepted());

        // Genuine receipt, substituted candidate payload.
        let receipt = graph.admission_receipt(certified).expect("mint receipt");
        let substituted = Color::Verified { lo: 0.0, hi: 9.0 };
        assert!(!verifier.verify(&substituted, &receipt).accepted());

        // Genuine candidate, receipt for a DIFFERENT node.
        let anchored_receipt = graph.admission_receipt(anchored).expect("mint receipt");
        assert!(!verifier.verify(&declared, &anchored_receipt).accepted());

        // Wrong policy fingerprint and wrong schema version.
        let wrong_policy = fs_evidence::AdmissionReceipt::from_parts(
            receipt.node_hash(),
            COLOR_WRITE_ROW_SCHEMA_VERSION,
            COLOR_ALGEBRA_VERSION,
            hash_bytes(b"someone else's policy"),
        );
        assert!(!verifier.verify(&declared, &wrong_policy).accepted());
        let wrong_schema = fs_evidence::AdmissionReceipt::from_parts(
            receipt.node_hash(),
            COLOR_WRITE_ROW_SCHEMA_VERSION + 1,
            COLOR_ALGEBRA_VERSION,
            policy,
        );
        assert!(!verifier.verify(&declared, &wrong_schema).accepted());
    }

    #[test]
    fn tampered_node_state_refuses_at_mint_and_at_verify() {
        let (mut graph, certified, _, _) = admission_fixture();
        let declared = graph
            .node(certified)
            .expect("node")
            .scientific_color()
            .expect("unwaived")
            .clone();
        // Mint against clean state, then tamper the stored color without
        // recomputing the provenance hash.
        let receipt = graph.admission_receipt(certified).expect("mint receipt");
        graph.nodes[usize::try_from(certified).expect("small id")].color =
            Color::Verified { lo: 0.5, hi: 9.5 };

        let refusal = graph
            .admission_receipt(certified)
            .expect_err("tampered node must not mint");
        assert!(matches!(
            refusal,
            ColorAdmissionRefusal::ReplayDivergence(_)
        ));

        let verifier = LedgerColorAdmissionVerifier::new(&graph);
        assert!(
            !verifier.verify(&declared, &receipt).accepted(),
            "a pre-tamper receipt must die with the tampered graph"
        );
    }

    #[test]
    fn out_of_regime_state_refuses_admission_structurally() {
        let (graph, _, anchored, _) = admission_fixture();
        let outside = BTreeMap::from([("re".to_string(), 5e2)]);
        let refusal = graph
            .admission_receipt_in_regime(anchored, &outside)
            .expect_err("out-of-regime must refuse");
        assert!(matches!(
            refusal,
            ColorAdmissionRefusal::OutOfRegime { node, .. } if node == anchored
        ));

        let inside = BTreeMap::from([("re".to_string(), 5e4)]);
        graph
            .admission_receipt_in_regime(anchored, &inside)
            .expect("in-regime state mints");
    }

    #[test]
    fn stale_algebra_receipts_refuse_before_the_capability_decides() {
        let (graph, certified, _, _) = admission_fixture();
        let verifier = LedgerColorAdmissionVerifier::new(&graph);
        let genuine = graph.admission_receipt(certified).expect("mint receipt");
        let declared = graph
            .node(certified)
            .expect("node")
            .scientific_color()
            .expect("unwaived")
            .clone();
        let stale = fs_evidence::AdmissionReceipt::from_parts(
            genuine.node_hash(),
            genuine.row_schema_version(),
            COLOR_ALGEBRA_VERSION + 1,
            genuine.policy_fingerprint(),
        );
        let error = AdmittedColor::from_receipt(declared, stale, &verifier)
            .expect_err("stale algebra must refuse");
        assert_eq!(
            error,
            AdmissionRejection::StaleAlgebra {
                receipt: COLOR_ALGEBRA_VERSION + 1,
                current: COLOR_ALGEBRA_VERSION
            }
        );
    }
}
