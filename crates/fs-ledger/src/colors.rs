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
    Color, ColorRank, Demotion, IntervalOp, NumericalCertificate, ValidityDomain, check_regime,
    compose, verified_from,
};
use std::collections::BTreeMap;

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
/// gp3.16). Mirrors the schema-v6 claim-origin vocabulary
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
                if identity_reason(dataset_id).is_some() {
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
    /// A Verified interval is non-finite or inverted.
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
const COLOR_NODE_HASH_DOMAIN: &[u8] = b"frankensim/fs-ledger/color-node";
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

fn identity_reason(value: &str) -> Option<&'static str> {
    let trimmed = value.trim();
    if trimmed.is_empty() {
        Some("blank")
    } else if trimmed != value {
        Some("surrounding-whitespace")
    } else if is_placeholder_token(value) {
        Some("placeholder")
    } else {
        None
    }
}

fn validate_color_structure(color: &Color) -> Result<(), ColorStructureRejection> {
    match color {
        Color::Verified { lo, hi } => {
            if !lo.is_finite() || !hi.is_finite() {
                Err(ColorStructureRejection::InvalidVerifiedInterval {
                    reason: "bounds must be finite",
                })
            } else if lo > hi {
                Err(ColorStructureRejection::InvalidVerifiedInterval {
                    reason: "lower bound exceeds upper bound",
                })
            } else {
                Ok(())
            }
        }
        Color::Validated { regime, dataset } => {
            if let Some(reason) = identity_reason(dataset) {
                return Err(ColorStructureRejection::InvalidIdentity {
                    field: "dataset",
                    value: dataset.clone(),
                    reason,
                });
            }
            if regime.bounds().is_empty() {
                return Err(ColorStructureRejection::InvalidValidatedRegime {
                    axis: String::new(),
                    reason: "at least one bounded axis is required",
                });
            }
            for (axis, (lo, hi)) in regime.bounds() {
                if let Some(reason) = identity_reason(axis) {
                    return Err(ColorStructureRejection::InvalidIdentity {
                        field: "axis",
                        value: axis.clone(),
                        reason,
                    });
                }
                if !lo.is_finite() || !hi.is_finite() {
                    return Err(ColorStructureRejection::InvalidValidatedRegime {
                        axis: axis.clone(),
                        reason: "bounds must be finite",
                    });
                }
                if lo > hi {
                    return Err(ColorStructureRejection::InvalidValidatedRegime {
                        axis: axis.clone(),
                        reason: "lower bound exceeds upper bound",
                    });
                }
            }
            Ok(())
        }
        Color::Estimated {
            estimator,
            dispersion,
        } => {
            if let Some(reason) = identity_reason(estimator) {
                return Err(ColorStructureRejection::InvalidIdentity {
                    field: "estimator",
                    value: estimator.clone(),
                    reason,
                });
            }
            if dispersion.is_nan() || *dispersion < 0.0 {
                Err(ColorStructureRejection::InvalidEstimatedDispersion {
                    reason: "value is NaN or negative",
                })
            } else {
                Ok(())
            }
        }
    }
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

    /// The color as WRITTEN (post demotion, post waiver).
    #[must_use]
    pub fn color(&self) -> &Color {
        &self.color
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
    ) -> Result<Vec<WaiverDependency>, ColorWriteError> {
        let mut dependencies = BTreeMap::<u64, WaiverDependency>::new();
        for parent_id in parents {
            let parent = self
                .node(*parent_id)
                .ok_or(ColorWriteError::UnknownParent { id: *parent_id })?;
            for dependency in &parent.waiver_dependencies {
                dependencies
                    .entry(dependency.authorizing_node)
                    .or_insert_with(|| dependency.clone());
            }
            if let Some(grant) = &parent.grant {
                let (Some(policy_fingerprint), Some(admission_day)) = (
                    parent.waiver_policy_fingerprint,
                    parent.waiver_admission_day,
                ) else {
                    return Err(ColorWriteError::InvalidParentAuthority { id: parent.id });
                };
                dependencies
                    .entry(parent.id)
                    .or_insert_with(|| WaiverDependency {
                        authorizing_node: parent.id,
                        operation: parent.operation,
                        grant: grant.clone(),
                        policy_fingerprint,
                        admission_day,
                    });
            }
        }
        Ok(dependencies.into_values().collect())
    }

    /// Provenance hash over DOMAIN-SEPARATED, VERSIONED v8,
    /// LENGTH-PREFIXED encoding. V8 binds certificate artifact identity,
    /// direct source/waiver policy fingerprints, waiver admission days, and
    /// those fields in the canonical transitive waiver dependency closure.
    /// V7 first bound that dependency closure. V6 bound every regime demotion and the
    /// correct source/derived waiver payload. V5 bound the typed SOURCE ORIGIN (bead
    /// gp3.16) so a forged or substituted origin changes the node
    /// identity and every downstream hash. V4 bound source/derived
    /// status and the exact [`IntervalOp`]; v3 added
    /// [`Color::canonical_bytes`]; the former v2 representation used
    /// rounded display JSON. Length-prefixing prevents adversarial text
    /// from colliding structurally.
    fn node_hash(
        &self,
        name: &str,
        color: &Color,
        parents: &[u64],
        metadata: &NodeHashMetadata<'_>,
    ) -> ContentHash {
        let mut buf = vec![8u8]; // encoding version
        push_field(&mut buf, COLOR_NODE_HASH_DOMAIN);
        match metadata.operation {
            Some(op) => {
                buf.push(1);
                buf.push(interval_op_tag(op));
            }
            None => buf.push(0),
        }
        push_field(&mut buf, name.as_bytes());
        push_field(&mut buf, &color.canonical_bytes());
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
        match metadata.origin {
            Some(origin) => {
                buf.push(1);
                push_source_origin(&mut buf, origin);
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
                "{{\"event\":\"demotion\",\"node\":{id},\"parent_index\":{},\
                 \"parent\":{},\
                 \"dataset\":{},\"axis\":{},\"value\":{}}}",
                demotion.parent_index,
                demotion.parent_id,
                json_string(&d.dataset),
                json_string(&d.axis),
                json_f64(d.value)
            ));
        }
        let operation_json =
            operation.map_or("null".to_string(), |op| json_string(interval_op_name(op)));
        self.rows.push(format!(
            "{{\"event\":\"color-write\",\"schema_version\":5,\"node\":{id},\
             \"name\":{},\"operation\":{},\"color\":\"{}\",\"payload\":{},\
             \"parents\":{:?},\"origin\":{},\"origin_policy_fingerprint\":{},\
             \"waiver_dependencies\":[{}],\"waiver\":{},\"grant\":{},\
             \"waiver_policy_fingerprint\":{},\"waiver_admission_day\":{},\
             \"hash\":\"{}\"}}",
            json_string(name),
            operation_json,
            color.name(),
            color.payload_json(),
            parents,
            origin_json(origin.as_ref()),
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
        if !matches!(color, Color::Estimated { .. }) {
            return Err(ColorWriteError::SourceOriginRequired { rank: color.rank() });
        }
        if let Some((field, why)) = estimated_payload_error(&color) {
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
    /// forged-source reason.
    pub fn source_with_origin(
        &mut self,
        name: &str,
        color: &Color,
        origin: SourceOrigin,
        verifier: &dyn SourceOriginVerifier,
    ) -> Result<u64, ColorWriteError> {
        let refuse = |rejection| Err(ColorWriteError::SourceOriginRefused { rejection });
        if matches!(color, Color::Estimated { .. }) {
            return refuse(SourceOriginRejection::EstimatedNeedsNoOrigin);
        }
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
        validate_color_structure(&color)
            .map_err(|rejection| ColorWriteError::InvalidColor { rejection })?;
        if color.rank() < ColorRank::Validated {
            return Err(ColorWriteError::SourceOriginRefused {
                rejection: SourceOriginRejection::EstimatedNeedsNoOrigin,
            });
        }
        let refuse = |rejection| Err(ColorWriteError::WaiverRefused { rejection });
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
    /// regimes at once.
    fn fold_parents(
        &self,
        parents: &[u64],
        op: IntervalOp,
        state: &BTreeMap<String, f64>,
    ) -> Result<(Color, Vec<ColorDemotion>), ColorWriteError> {
        if parents.is_empty() {
            return Err(ColorWriteError::NoParents);
        }
        let mut demotions = Vec::new();
        let mut effective: Vec<Color> = Vec::with_capacity(parents.len());
        for (parent_index, &p) in parents.iter().enumerate() {
            let parent = self
                .node(p)
                .ok_or(ColorWriteError::UnknownParent { id: p })?;
            let (c, d) = check_regime(&parent.color, state);
            if let Some(reason) = d {
                demotions.push(ColorDemotion {
                    parent_index,
                    parent_id: p,
                    reason,
                });
            }
            effective.push(c);
        }
        let mut derived = effective[0].clone();
        for c in &effective[1..] {
            derived = compose(&derived, c, op);
        }
        Ok((derived, demotions))
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
                let (eff, _) = check_regime(&parent.color, state);
                eff.rank() <= cap
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
        let (derived, demotions) = self.fold_parents(parents, op, state)?;
        let waiver_dependencies = self.inherited_waiver_dependencies(parents)?;
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
        validate_color_structure(&claimed)
            .map_err(|rejection| ColorWriteError::InvalidColor { rejection })?;
        let (_derived, demotions) = self.fold_parents(parents, op, state)?;
        let waiver_dependencies = self.inherited_waiver_dependencies(parents)?;
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
        }

        let expected = if node.parents.is_empty() {
            Ok(Vec::new())
        } else {
            self.inherited_waiver_dependencies(&node.parents)
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
            let Some(Color::Validated { regime, dataset }) =
                self.node(demotion.parent_id).map(ColorNode::color)
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
                    || node.waiver_policy_fingerprint.is_some()
                    || node.waiver_admission_day.is_some()
                {
                    return Err(Self::replay_error(
                        node,
                        "Estimated source carries orphan admission context",
                    ));
                }
                if let Some((field, why)) = estimated_payload_error(&node.color) {
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
        let mut effective = Vec::with_capacity(node.parents.len());
        for (index, parent) in node.parents.iter().enumerate() {
            let Some(parent_node) = self.node(*parent) else {
                return Err(Self::replay_error(node, "derived parent is missing"));
            };
            effective.push(
                node.demotions
                    .iter()
                    .find(|demotion| demotion.parent_index == index)
                    .map_or_else(
                        || parent_node.color.clone(),
                        |demotion| Color::Estimated {
                            estimator: format!(
                                "regime-exit:{}@{}",
                                demotion.reason.dataset, demotion.reason.axis
                            ),
                            dispersion: f64::INFINITY,
                        },
                    ),
            );
        }
        let Some((first, remaining)) = effective.split_first() else {
            return Err(Self::replay_error(node, "derived node has no parents"));
        };
        let mut derived = first.clone();
        for color in remaining {
            derived = compose(&derived, color, op);
        }
        if let Some(grant) = &node.grant {
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
        if let Err(rejection) = validate_color_structure(&node.color) {
            return Err(Self::replay_error(
                node,
                format!("stored color is structurally invalid: {rejection}"),
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

    /// REPLAY AUDIT (bead gp3.16): rederive every node from its stored
    /// inputs and refuse on any divergence. For each derived node the
    /// recorded demotions reconstruct the effective parent colors
    /// (a demotion determines the demoted form exactly:
    /// `estimated{regime-exit:dataset@axis, ∞}`), the composition
    /// algebra re-folds them, and — for unwaived writes — the written
    /// color must match bit-exactly. Every node's provenance hash is
    /// recomputed and compared, so the graph's whole hash chain is
    /// re-earned, never trusted. Positive-colored leaves must carry
    /// their typed origin or an authenticated grant (the sealed-source
    /// invariant, re-checked).
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

#[cfg(test)]
mod tests {
    use super::*;

    struct AllowFixtureSource;

    impl SourceOriginVerifier for AllowFixtureSource {
        fn verify(&self, request: &SourceOriginRequest<'_>) -> PolicyDecision {
            let policy = hash_bytes(b"fs-ledger/internal-fixture/source-policy/v1");
            if request.node_name() == "anchored"
                && matches!(request.claimed_color(), Color::Validated { .. })
                && matches!(request.origin(), SourceOrigin::Anchoring { .. })
            {
                PolicyDecision::accept(policy)
            } else {
                PolicyDecision::reject(policy)
            }
        }
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
        assert!(error.why.contains("bounds must be finite"));
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
            .inherited_waiver_dependencies(&[source])
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
}
