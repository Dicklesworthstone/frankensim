//! fs-package — machine-checkable evidence packages (plan addendum,
//! Proposal 12). Layer: L6.
//!
//! When FrankenSim asserts "this design meets spec", the assertion travels as a
//! transport-complete, CONTENT-ADDRESSED bundle: the color-typed claims, the raw
//! certificate data behind each (carried in the [`fs_evidence::Color`]),
//! provenance (code version + the constellation lockfile), and a Merkle root
//! over the package identity so any tamper is detectable. A standalone,
//! open-source CHECKER re-verifies the package WITHOUT re-running every solver.
//! It is not self-authenticating: source, anchor, falsifier, derivation, waiver,
//! and signature artifacts require explicit [`VerificationCapabilities`] through
//! [`EvidencePackage::verify_with`].
//!
//! Completeness is enforced, not assumed: a validated-color claim that is
//! missing its regime tag OR its anchoring dataset FAILS verification (an
//! unfalsifiable "validated" claim is worse than none). An all-estimated
//! package is still valid and round-trips — honesty about low confidence is
//! not a defect.
//!
//! The Merkle tree uses the in-house BLAKE3 content hash from [`fs_blake3`]
//! (pure safe Rust, zero deps — Franken-compliant), with every leaf and node
//! DOMAIN-SEPARATED under `fs-package:v7:…` tags, yielding a
//! 32-byte [`ContentHash`] root; signature bytes are DETACHED and OPTIONAL, and
//! become authenticated only through an injected purpose-bound policy.
//! Everything is deterministic: the same package yields the same root and
//! JSON.

use fs_blake3::hash_domain;
use fs_evidence::{
    Color, ColorPayloadError, ColorRank, IntervalOp, compose, validate_color_payload,
};
use origin::{identity_reason, is_placeholder_token, validate_origin_shape};

pub use fs_blake3::ContentHash;

pub mod coverage;
pub mod origin;
pub use coverage::{
    ConceptPresence, CoverageStatus, PackageCoverageReport, PackagePresenceReport,
    package_coverage, package_coverage_with, package_presence, package_presence_with,
    verified_package_coverage, verified_package_presence,
};
pub use origin::{
    AnchoredSourceRequest, AnchoredSourceVerifier, ClaimOrigin, DerivationRequest,
    DerivationVerifier, FalsifierRequest, FalsifierVerifier, NoAnchoredSourceVerifier,
    NoDerivationVerifier, NoFalsifierVerifier, NoSignatureVerifier, NoSourceCertificateVerifier,
    NoWaiverVerifier, OriginError, PolicyFingerprint, SignatureIntent, SignaturePurpose,
    SignatureRequest, SignatureVerification, SignatureVerifier, SourceCertificateRequest,
    SourceCertificateVerifier, VerificationCapabilities, VerificationDecision, WaiverGrant,
    WaiverVerification, WaiverVerifier, signature_subject_hash,
};

/// A COMPOSITION RECEIPT (schema v3, bead xfxq): this claim's color was
/// derived from earlier claims in the package, and the standalone
/// checker re-runs the derivation — `compose` folded over the parents'
/// colors in order must EQUAL the claimed color exactly. Parents are
/// indices into the package's claim list and must precede this claim
/// (a DAG by construction).
#[derive(Debug, Clone, PartialEq)]
pub struct CompositionReceipt {
    /// Exact fs-evidence color algebra used to derive this claim.
    pub color_algebra_version: u32,
    /// Parent claim indices, in fold order (each < this claim's index).
    pub parents: Vec<usize>,
    /// The ledger operation the derivation used.
    pub op: IntervalOp,
    /// Canonical 64-hex address of the derivation proof artifact.
    pub artifact_hash: String,
}

/// One falsifier's adversarial record against a claim (schema v3):
/// negative results travel WITH the claim; a refuted claim fails
/// verification outright.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FalsifierRecord {
    /// Stable registered identity of the falsifier that ran (meaningful,
    /// non-placeholder text).
    pub name: String,
    /// Adversarial attempts executed (strictly positive).
    pub attempts: u64,
    /// Did it refute the claim?
    pub refuted: bool,
    /// Meaningful, non-placeholder outcome summary.
    pub detail: String,
    /// Canonical 64-hex content address of the executable falsifier artifact
    /// and retained results represented by this record.
    pub artifact_hash: String,
}

/// An anchoring-dataset identity (schema v3): the reference data behind
/// a validated claim, by stable id and content hash — not just a name.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AnchorRecord {
    /// Stable, non-blank dataset identity.
    pub dataset_id: String,
    /// Canonical 64-character lowercase hex hash of the dataset artifact.
    pub content_hash: String,
}

/// One claim in an evidence package: a statement plus its epistemic color
/// (which carries the certificate data — an interval, a regime+dataset, or an
/// estimator+dispersion), and optionally its composition receipt,
/// falsifier records, and dataset anchors (schema v3).
#[derive(Debug, Clone, PartialEq)]
pub struct Claim {
    /// SEALED (schema v6): fields are crate-private so a claim
    /// can only exist through the origin-typed constructors below — public
    /// `Color::Verified` alone can no longer mint a checker-passing claim.
    pub(crate) id: String,
    pub(crate) statement: String,
    pub(crate) color: Color,
    pub(crate) receipt: Option<CompositionReceipt>,
    pub(crate) falsifiers: Vec<FalsifierRecord>,
    pub(crate) anchors: Vec<AnchorRecord>,
    pub(crate) origin: ClaimOrigin,
}

impl Claim {
    fn sealed(
        id: impl Into<String>,
        statement: impl Into<String>,
        color: Color,
        origin: ClaimOrigin,
    ) -> Claim {
        Claim {
            id: id.into(),
            statement: statement.into(),
            color,
            receipt: None,
            falsifiers: Vec::new(),
            anchors: Vec::new(),
            origin,
        }
    }

    /// A VERIFIED claim from a named producer's certificate artifact.
    #[must_use]
    pub fn from_certificate(
        id: impl Into<String>,
        statement: impl Into<String>,
        lo: f64,
        hi: f64,
        producer: impl Into<String>,
        certificate_hash: impl Into<String>,
    ) -> Claim {
        Claim::sealed(
            id,
            statement,
            Color::Verified { lo, hi },
            ClaimOrigin::SourceCertificate {
                producer: producer.into(),
                certificate_hash: certificate_hash.into(),
            },
        )
    }

    /// A VALIDATED claim anchored to its reference dataset: the origin
    /// names the color's dataset and a matching content-hash anchor
    /// record is attached automatically.
    #[must_use]
    pub fn anchored(
        id: impl Into<String>,
        statement: impl Into<String>,
        regime: fs_evidence::ValidityDomain,
        dataset: impl Into<String>,
        content_hash: impl Into<String>,
    ) -> Claim {
        let dataset = dataset.into();
        let content_hash = content_hash.into();
        let mut claim = Claim::sealed(
            id,
            statement,
            Color::Validated {
                regime,
                dataset: dataset.clone(),
            },
            ClaimOrigin::AnchoredSource {
                dataset_id: dataset.clone(),
                content_hash: content_hash.clone(),
            },
        );
        claim.anchors.push(AnchorRecord {
            dataset_id: dataset,
            content_hash,
        });
        claim
    }

    /// An ESTIMATED claim from a named estimator.
    #[must_use]
    pub fn estimated(
        id: impl Into<String>,
        statement: impl Into<String>,
        estimator: impl Into<String>,
        dispersion: f64,
    ) -> Claim {
        let estimator = estimator.into();
        Claim::sealed(
            id,
            statement,
            Color::Estimated {
                estimator: estimator.clone(),
                dispersion,
            },
            ClaimOrigin::EstimatedSource { estimator },
        )
    }

    /// A DERIVED claim: its color must re-derive bit-exactly from the
    /// named parents under `op` (the checker re-runs the fold).
    #[must_use]
    pub fn derived(
        id: impl Into<String>,
        statement: impl Into<String>,
        color: Color,
        parents: Vec<usize>,
        op: IntervalOp,
        artifact_hash: impl Into<String>,
    ) -> Claim {
        let mut claim = Claim::sealed(id, statement, color, ClaimOrigin::Derived);
        claim.receipt = Some(CompositionReceipt {
            color_algebra_version: fs_evidence::COLOR_ALGEBRA_VERSION,
            parents,
            op,
            artifact_hash: artifact_hash.into(),
        });
        claim
    }

    /// A WAIVED claim: any color, authorized only by an explicit,
    /// expiring, MAC'd grant that an INJECTED verifier must accept.
    #[must_use]
    pub fn waived(
        id: impl Into<String>,
        statement: impl Into<String>,
        color: Color,
        grant: WaiverGrant,
    ) -> Claim {
        Claim::sealed(
            id,
            statement,
            color,
            ClaimOrigin::AuthenticatedWaiver(grant),
        )
    }

    /// Read-only accessors (the sealed fields' public view).
    #[must_use]
    pub fn id(&self) -> &str {
        &self.id
    }
    /// The human-readable claim text.
    #[must_use]
    pub fn statement(&self) -> &str {
        &self.statement
    }
    /// The epistemic color + certificate payload.
    #[must_use]
    pub fn declared_color_unverified(&self) -> &Color {
        &self.color
    }
    /// Domain-separated hash of the complete raw claim declaration. This is an
    /// artifact address, not an admission decision.
    #[must_use]
    pub fn declared_content_hash_unverified(&self) -> ContentHash {
        hash_domain("fs-package:v7:claim", self.canonical().as_bytes())
    }

    /// Domain-separated subject hash for external falsifier/derivation
    /// artifacts. External artifact addresses and waiver MAC bytes are omitted
    /// so an artifact may embed this digest without a fixed-point cycle; the
    /// package root separately binds the complete declaration.
    #[must_use]
    pub fn declared_verification_subject_hash_unverified(&self) -> ContentHash {
        let mut subject = self.clone();
        if let Some(receipt) = &mut subject.receipt {
            receipt.artifact_hash.clear();
        }
        for falsifier in &mut subject.falsifiers {
            falsifier.artifact_hash.clear();
        }
        hash_domain(
            "fs-package:v7:claim-verification-subject",
            subject.authorization_canonical().as_bytes(),
        )
    }
    /// The composition receipt, when derived.
    #[must_use]
    pub fn declared_receipt_unverified(&self) -> Option<&CompositionReceipt> {
        self.receipt.as_ref()
    }
    /// Attached falsifier records.
    #[must_use]
    pub fn declared_falsifiers_unverified(&self) -> &[FalsifierRecord] {
        &self.falsifiers
    }
    /// Attached anchor records.
    #[must_use]
    pub fn declared_anchors_unverified(&self) -> &[AnchorRecord] {
        &self.anchors
    }
    /// Where this claim's certificate came from.
    #[must_use]
    pub fn declared_origin_unverified(&self) -> &ClaimOrigin {
        &self.origin
    }

    /// Attach a falsifier record (builder style).
    #[must_use]
    pub fn with_falsifier(mut self, rec: FalsifierRecord) -> Claim {
        self.falsifiers.push(rec);
        self
    }

    /// Attach a dataset anchor (builder style).
    #[must_use]
    pub fn with_anchor(
        mut self,
        dataset_id: impl Into<String>,
        content_hash: impl Into<String>,
    ) -> Claim {
        self.anchors.push(AnchorRecord {
            dataset_id: dataset_id.into(),
            content_hash: content_hash.into(),
        });
        self
    }

    /// Whether this validated claim carries a canonical content-hash anchor
    /// for the exact dataset named by its color. Other color classes return
    /// `false` because they have no validated dataset to anchor.
    #[must_use]
    pub fn has_declared_matching_validated_anchor_unverified(&self) -> bool {
        let Color::Validated { dataset, .. } = &self.color else {
            return false;
        };
        let required_origin_hash = match &self.origin {
            ClaimOrigin::AnchoredSource {
                dataset_id,
                content_hash,
            } if dataset_id == dataset => Some(content_hash.as_str()),
            ClaimOrigin::AnchoredSource { .. } => return false,
            ClaimOrigin::Derived | ClaimOrigin::AuthenticatedWaiver(_) => None,
            ClaimOrigin::SourceCertificate { .. } | ClaimOrigin::EstimatedSource { .. } => {
                return false;
            }
        };
        self.anchors.iter().any(|anchor| {
            anchor.dataset_id == *dataset
                && is_canonical_content_hash(&anchor.content_hash)
                && required_origin_hash.is_none_or(|hash| anchor.content_hash == hash)
        })
    }

    /// Whether this claim is a certificate-class result subject to the
    /// no-falsifier-no-ship release rule. Estimated claims remain explicitly
    /// low-assurance rather than being promoted into this class.
    #[must_use]
    pub fn declared_requires_release_falsifier_unverified(&self) -> bool {
        matches!(
            &self.color,
            Color::Verified { .. } | Color::Validated { .. }
        )
    }

    /// Whether this raw declaration can satisfy the release gate's minimum
    /// substantive-evidence requirement after its admission receipt classifies
    /// it as scientific. This predicate does not authenticate the declaration.
    ///
    /// A finite `Verified` interval is informative even when wide. An ordered
    /// interval with an infinite endpoint remains a sound, explicit vacuous
    /// enclosure, but it cannot by itself justify release. `Validated` evidence
    /// remains eligible because package verification separately authenticates
    /// its bounded regime and exact anchoring dataset.
    #[must_use]
    pub fn declared_is_release_scientific_evidence_unverified(&self) -> bool {
        matches!(
            &self.color,
            Color::Verified { lo, hi } if lo.is_finite() && hi.is_finite()
        ) || matches!(&self.color, Color::Validated { .. })
    }

    /// Whether release admission must find a matching content-hash dataset
    /// anchor for this claim.
    #[must_use]
    pub fn declared_requires_validated_anchor_unverified(&self) -> bool {
        matches!(&self.color, Color::Validated { .. })
    }

    /// The schema-v7 canonical body (id, statement, color, receipt,
    /// falsifiers, anchors), excluding the claim origin.
    fn canonical_body(&self) -> String {
        use core::fmt::Write as _;
        let mut out = String::from("claim|");
        push_atom(&mut out, &self.id);
        push_atom(&mut out, &self.statement);
        match &self.color {
            Color::Verified { lo, hi } => {
                out.push_str("verified|");
                let _ = write!(out, "{}|{}|", lo.to_bits(), hi.to_bits());
            }
            Color::Validated { regime, dataset } => {
                out.push_str("validated|");
                for (k, (lo, hi)) in regime.bounds() {
                    push_atom(&mut out, k);
                    let _ = write!(out, "{}|{}|", lo.to_bits(), hi.to_bits());
                }
                push_atom(&mut out, dataset);
            }
            Color::Estimated {
                estimator,
                dispersion,
            } => {
                out.push_str("estimated|");
                push_atom(&mut out, estimator);
                let _ = write!(out, "{}|", dispersion.to_bits());
            }
        }
        // Schema-v3 fields bind into the content address too.
        match &self.receipt {
            Some(r) => {
                let _ = write!(
                    out,
                    "receipt:color-algebra-v{}:{}|",
                    r.color_algebra_version,
                    op_name(r.op)
                );
                for &p in &r.parents {
                    let _ = write!(out, "{p}|");
                }
                push_atom(&mut out, &r.artifact_hash);
            }
            None => out.push_str("no-receipt|"),
        }
        for fr in &self.falsifiers {
            out.push_str("falsifier|");
            push_atom(&mut out, &fr.name);
            let _ = write!(out, "{}|{}|", fr.attempts, fr.refuted);
            push_atom(&mut out, &fr.detail);
            push_atom(&mut out, &fr.artifact_hash);
        }
        for a in &self.anchors {
            out.push_str("anchor|");
            push_atom(&mut out, &a.dataset_id);
            push_atom(&mut out, &a.content_hash);
        }
        out
    }

    /// Full canonical string (schema v7): the claim body plus the origin.
    fn canonical(&self) -> String {
        let mut out = self.canonical_body();
        out.push_str("origin|");
        for part in self.origin.canonical_parts() {
            push_atom(&mut out, &part);
        }
        out
    }

    /// Canonical authorization context. It differs from the content-address
    /// form only by omitting waiver MAC bytes, which makes it possible to
    /// compute a stable message before installing the final authenticator.
    fn authorization_canonical(&self) -> String {
        use core::fmt::Write as _;

        let mut out = self.canonical_body();
        out.push_str("origin|");
        match &self.origin {
            ClaimOrigin::AuthenticatedWaiver(grant) => {
                push_atom(&mut out, self.origin.kind());
                push_atom(&mut out, &grant.waiver_id);
                let _ = write!(out, "{}|", grant.expiry_day);
            }
            _ => {
                for part in self.origin.canonical_parts() {
                    push_atom(&mut out, &part);
                }
            }
        }
        out
    }
}

/// Stable op name for hashing/JSON.
fn op_name(op: IntervalOp) -> &'static str {
    match op {
        IntervalOp::Add => "add",
        IntervalOp::Mul => "mul",
        IntervalOp::Hull => "hull",
    }
}

fn op_parse(name: &str) -> Option<IntervalOp> {
    match name {
        "add" => Some(IntervalOp::Add),
        "mul" => Some(IntervalOp::Mul),
        "hull" => Some(IntervalOp::Hull),
        _ => None,
    }
}

/// Where a package came from — enough to reproduce it.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Provenance {
    /// The code version / commit that produced the claims.
    pub code_version: String,
    /// The pinned dependency constellation (lockfile digest).
    pub constellation_lock: String,
}

impl Provenance {
    /// Provenance.
    #[must_use]
    pub fn new(
        code_version: impl Into<String>,
        constellation_lock: impl Into<String>,
    ) -> Provenance {
        Provenance {
            code_version: code_version.into(),
            constellation_lock: constellation_lock.into(),
        }
    }
}

/// A transport-complete, content-addressed evidence bundle whose external
/// artifacts require explicit verification capabilities.
#[derive(Debug, Clone, PartialEq)]
pub struct EvidencePackage {
    /// The format version (stability promise for external checkers).
    pub format_version: u32,
    /// The claims, in order.
    claims: Vec<Claim>,
    /// Provenance.
    pub provenance: Provenance,
    /// An OPTIONAL detached signature over a canonical typed signature subject.
    pub signature: Option<String>,
}

/// The by-color budget pie over a package's ADMITTED claims.
///
/// The first three buckets contain scientific claims only. A directly waived
/// claim and every derived descendant of one are counted exclusively in
/// [`ColorBreakdown::waived`], irrespective of their underlying color.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct ColorBreakdown {
    /// Verified-color claims.
    pub verified: usize,
    /// Validated-color claims.
    pub validated: usize,
    /// Estimated-color claims.
    pub estimated: usize,
    /// Waiver-dependent claims, including the complete derived descendant cone.
    pub waived: usize,
}

/// Whether a verified claim is admitted as scientific evidence or remains
/// dependent on one or more administrative waivers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdmissionClass {
    /// No waiver occurs in this claim's transitive derivation ancestry.
    Scientific,
    /// This claim is directly waived or derived from a waiver-dependent claim.
    WaiverDependent,
}

/// Stable origin class captured in a verification receipt.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdmissionOriginKind {
    /// Externally re-verified certificate artifact.
    SourceCertificate,
    /// Externally re-verified anchoring dataset.
    AnchoredSource,
    /// Self-described estimate identity.
    EstimatedSource,
    /// Re-derived composition receipt.
    Derived,
    /// Authenticated administrative waiver.
    AuthenticatedWaiver,
}

impl AdmissionOriginKind {
    const fn tag(self) -> &'static str {
        match self {
            AdmissionOriginKind::SourceCertificate => "source-certificate",
            AdmissionOriginKind::AnchoredSource => "anchored-source",
            AdmissionOriginKind::EstimatedSource => "estimated-source",
            AdmissionOriginKind::Derived => "derived",
            AdmissionOriginKind::AuthenticatedWaiver => "authenticated-waiver",
        }
    }
}

/// Auditable admission decision for one claim in topological package order.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ClaimAdmission {
    /// Stable claim position.
    claim_index: usize,
    /// Stable claim identity.
    claim_id: String,
    /// Origin mechanism whose semantics were admitted for this claim.
    origin_kind: AdmissionOriginKind,
    /// Scientific or waiver-dependent admission class.
    class: AdmissionClass,
    /// Registry entry when this claim is directly waived.
    direct_waiver: Option<usize>,
    /// Immediate receipt parents that are waiver-dependent. Traversing these
    /// edges reaches every transitive waiver without copying waiver-id strings
    /// into every descendant.
    waiver_parents: Vec<usize>,
}

impl ClaimAdmission {
    /// Stable claim position in package order.
    #[must_use]
    pub const fn claim_index(&self) -> usize {
        self.claim_index
    }
    /// Stable claim identity.
    #[must_use]
    pub fn claim_id(&self) -> &str {
        &self.claim_id
    }
    /// Origin mechanism admitted for this claim.
    #[must_use]
    pub const fn origin_kind(&self) -> AdmissionOriginKind {
        self.origin_kind
    }
    /// Scientific or waiver-dependent admission class.
    #[must_use]
    pub const fn class(&self) -> AdmissionClass {
        self.class
    }
    /// Interned direct-waiver index, when directly waived.
    #[must_use]
    pub const fn direct_waiver(&self) -> Option<usize> {
        self.direct_waiver
    }
    /// Immediate waiver-dependent parent claim indices.
    #[must_use]
    pub fn waiver_parents(&self) -> &[usize] {
        &self.waiver_parents
    }
}

/// One interned waiver identity in a verification receipt.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReceiptWaiver {
    /// Stable registry position referenced by `ClaimAdmission::direct_waiver`.
    registry_index: usize,
    /// Claim that directly carries this waiver.
    claim_index: usize,
    /// Authenticated waiver identity, stored exactly once in the receipt.
    waiver_id: String,
}

impl ReceiptWaiver {
    /// Stable position in the receipt waiver registry.
    #[must_use]
    pub const fn registry_index(&self) -> usize {
        self.registry_index
    }
    /// Claim index that directly carries this waiver.
    #[must_use]
    pub const fn claim_index(&self) -> usize {
        self.claim_index
    }
    /// Authenticated waiver identity.
    #[must_use]
    pub fn waiver_id(&self) -> &str {
        &self.waiver_id
    }
}

/// Stable identities of all external policies made available to one decision.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct VerificationPolicyFingerprints {
    /// Source-certificate policy, when invoked.
    source_certificates: Option<PolicyFingerprint>,
    /// Anchoring-dataset policy, when invoked.
    anchored_sources: Option<PolicyFingerprint>,
    /// Falsifier-artifact policy, when invoked.
    falsifiers: Option<PolicyFingerprint>,
    /// Derivation-artifact policy, when invoked.
    derivations: Option<PolicyFingerprint>,
    /// Waiver-authorization policy, when invoked.
    waivers: Option<PolicyFingerprint>,
    /// Detached-signature policy, when invoked.
    signatures: Option<PolicyFingerprint>,
}

impl VerificationPolicyFingerprints {
    /// Source-certificate policy identity, when invoked.
    #[must_use]
    pub const fn source_certificates(&self) -> Option<PolicyFingerprint> {
        self.source_certificates
    }
    /// Anchoring-source policy identity, when invoked.
    #[must_use]
    pub const fn anchored_sources(&self) -> Option<PolicyFingerprint> {
        self.anchored_sources
    }
    /// Falsifier policy identity, when invoked.
    #[must_use]
    pub const fn falsifiers(&self) -> Option<PolicyFingerprint> {
        self.falsifiers
    }
    /// Derivation policy identity, when invoked.
    #[must_use]
    pub const fn derivations(&self) -> Option<PolicyFingerprint> {
        self.derivations
    }
    /// Waiver policy identity, when invoked.
    #[must_use]
    pub const fn waivers(&self) -> Option<PolicyFingerprint> {
        self.waivers
    }
    /// Signature policy identity, when invoked.
    #[must_use]
    pub const fn signatures(&self) -> Option<PolicyFingerprint> {
        self.signatures
    }
}

/// Read-only authenticated signature payload. Private fields prevent safe
/// downstream code from substituting a different signature or purpose into a
/// genuine authentication decision.
///
/// ```compile_fail
/// use fs_package::{AuthenticatedSignature, SignaturePurpose};
///
/// // Authentication payloads can only be created by package verification.
/// let forged = AuthenticatedSignature {
///     signature: "forged".to_string(),
///     purpose: SignaturePurpose::PackageRootAttestation,
/// };
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AuthenticatedSignature {
    signature: String,
    purpose: SignaturePurpose,
}

impl AuthenticatedSignature {
    /// Detached signature bytes authenticated by the recorded policy.
    #[must_use]
    pub fn signature(&self) -> &str {
        &self.signature
    }

    /// Exact domain/gate purpose authenticated by the recorded policy.
    #[must_use]
    pub const fn purpose(&self) -> SignaturePurpose {
        self.purpose
    }
}

/// Detached-signature decision made during package verification. Positive
/// authority is meaningful only inside a sealed [`VerificationReceipt`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SignatureStatus {
    /// No detached signature was present.
    Unsigned,
    /// Signature bytes may have been present, but the enclosing transport was
    /// refused before bounded authentication. Raw rejected bytes are omitted.
    Refused {
        /// Stable bounded refusal reason.
        reason: &'static str,
    },
    /// Signature bytes were present but no verifier was supplied.
    Unverified(String),
    /// The injected verifier authenticated the canonical typed signature subject.
    Authenticated(AuthenticatedSignature),
}

/// Replayable record of one package admission decision.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VerificationReceipt {
    /// Exact package root to which this decision applies.
    package_root: ContentHash,
    /// External policy identities used for the decision.
    policy_fingerprints: VerificationPolicyFingerprints,
    /// Explicit waiver clock context, when a waiver policy was installed.
    waiver_day: Option<u64>,
    /// Detached-signature decision.
    signature: SignatureStatus,
    /// Per-claim scientific/waiver-dependent decisions in package order.
    admissions: Vec<ClaimAdmission>,
    /// Interned direct waiver identities referenced by admissions.
    waiver_registry: Vec<ReceiptWaiver>,
    /// Domain-separated digest over every field above.
    receipt_hash: ContentHash,
}

impl VerificationReceipt {
    /// Exact package root admitted by this decision.
    #[must_use]
    pub const fn package_root(&self) -> ContentHash {
        self.package_root
    }
    /// External policy identities actually invoked by this decision.
    #[must_use]
    pub const fn policy_fingerprints(&self) -> &VerificationPolicyFingerprints {
        &self.policy_fingerprints
    }
    /// Waiver clock day, only when a waiver policy was invoked.
    #[must_use]
    pub const fn waiver_day(&self) -> Option<u64> {
        self.waiver_day
    }
    /// Detached-signature decision.
    #[must_use]
    pub fn signature(&self) -> &SignatureStatus {
        &self.signature
    }
    /// Per-claim admission decisions in package order.
    #[must_use]
    pub fn admissions(&self) -> &[ClaimAdmission] {
        &self.admissions
    }
    /// Interned direct-waiver identities.
    #[must_use]
    pub fn waiver_registry(&self) -> &[ReceiptWaiver] {
        &self.waiver_registry
    }
    /// Stored domain-separated receipt digest.
    #[must_use]
    pub const fn receipt_hash(&self) -> ContentHash {
        self.receipt_hash
    }

    /// Canonical pre-signature context for release approval. Producers obtain
    /// this from an unsigned verification pass, construct a concrete
    /// [`SignaturePurpose::ReleaseApproval`], sign [`signature_subject_hash`],
    /// attach the bytes, and then run the final release gate.
    #[must_use]
    pub fn release_admission_context(&self) -> ContentHash {
        release_admission_context_hash(
            self.package_root,
            self.policy_fingerprints,
            self.waiver_day,
            &self.admissions,
            &self.waiver_registry,
        )
    }
    /// Recompute the integrity digest over every receipt field.
    ///
    /// This is an unkeyed integrity check, not independent authenticity; trust
    /// still comes from replaying the named policies against the package.
    #[must_use]
    pub fn recomputed_hash(&self) -> ContentHash {
        verification_receipt_hash(
            self.package_root,
            self.policy_fingerprints,
            self.waiver_day,
            &self.signature,
            &self.admissions,
            &self.waiver_registry,
        )
    }

    /// Whether the stored receipt digest matches all receipt fields.
    #[must_use]
    pub fn validate_hash(&self) -> bool {
        self.receipt_hash == self.recomputed_hash()
    }
}

/// The result of verifying a package.
#[derive(Debug, Clone, PartialEq)]
pub struct PackageReport {
    /// The recomputed content address (domain-separated BLAKE3 Merkle root).
    merkle_root: ContentHash,
    /// The by-color budget pie.
    breakdown: ColorBreakdown,
    /// The number of claims.
    claims: usize,
    /// Scientific magnitude accounting with waiver-dependent values excluded.
    magnitude_budget: MagnitudeBudget,
    /// Policy-bound, replayable per-claim admission receipt.
    receipt: VerificationReceipt,
}

impl PackageReport {
    /// Recomputed package content root.
    #[must_use]
    pub const fn merkle_root(&self) -> ContentHash {
        self.merkle_root
    }
    /// Admitted scientific/waiver-dependent color counts.
    #[must_use]
    pub const fn breakdown(&self) -> &ColorBreakdown {
        &self.breakdown
    }
    /// Number of admitted claims.
    #[must_use]
    pub const fn claims(&self) -> usize {
        self.claims
    }
    /// Admitted magnitude budget with waiver values excluded.
    #[must_use]
    pub const fn magnitude_budget(&self) -> MagnitudeBudget {
        self.magnitude_budget
    }
    /// Policy-bound verification receipt.
    #[must_use]
    pub const fn receipt(&self) -> &VerificationReceipt {
        &self.receipt
    }
}

/// A checked package paired inseparably with the report and policy receipt that
/// admitted its gated evidence. This does not imply authorship: the nested
/// `SignatureStatus` may still be `Unsigned` or `Unverified`.
#[derive(Debug, Clone, PartialEq)]
pub struct VerifiedPackage {
    package: EvidencePackage,
    report: PackageReport,
}

impl VerifiedPackage {
    /// Structurally parsed package. Its raw declarations remain explicitly
    /// unverified when accessed outside the admission view below.
    #[must_use]
    pub fn package(&self) -> &EvidencePackage {
        &self.package
    }

    /// Successful verification report and policy-bound receipt.
    #[must_use]
    pub fn report(&self) -> &PackageReport {
        &self.report
    }

    /// Recheck structural binding between the retained package and report.
    #[must_use]
    pub fn validate_binding(&self) -> bool {
        let (admissions, waiver_registry) = self.package.admission_decisions();
        let breakdown = self.package.admitted_color_breakdown(&admissions);
        let magnitude = self.package.magnitude_budget_from(&admissions);
        self.report.receipt.validate_hash()
            && self.package.try_merkle_root() == Ok(self.report.merkle_root)
            && self.report.receipt.package_root() == self.report.merkle_root
            && self.report.claims == self.package.claims.len()
            && self.report.receipt.admissions == admissions
            && self.report.receipt.waiver_registry == waiver_registry
            && self.report.breakdown == breakdown
            && magnitude_budgets_bitwise_equal(self.report.magnitude_budget, magnitude)
    }

    /// Claims paired with their admission decisions in topological order.
    #[must_use]
    pub fn admitted_claims(&self) -> impl ExactSizeIterator<Item = AdmittedClaim<'_>> {
        self.package
            .claims
            .iter()
            .zip(&self.report.receipt.admissions)
            .map(|(claim, admission)| AdmittedClaim { claim, admission })
    }
}

fn magnitude_budgets_bitwise_equal(left: MagnitudeBudget, right: MagnitudeBudget) -> bool {
    left.verified_width.to_bits() == right.verified_width.to_bits()
        && left.estimated_dispersion.to_bits() == right.estimated_dispersion.to_bits()
        && left.validated_unquantified == right.validated_unquantified
        && left.waived_unquantified == right.waived_unquantified
        && left.quantified_total.to_bits() == right.quantified_total.to_bits()
}

/// Read-only admitted view of one claim.
#[derive(Debug, Clone, Copy)]
pub struct AdmittedClaim<'a> {
    claim: &'a Claim,
    admission: &'a ClaimAdmission,
}

impl AdmittedClaim<'_> {
    /// Claim identity.
    #[must_use]
    pub fn id(&self) -> &str {
        &self.claim.id
    }

    /// Human-readable assertion bound to the admitted origin/artifacts.
    #[must_use]
    pub fn statement(&self) -> &str {
        &self.claim.statement
    }

    /// Auditable admission decision.
    #[must_use]
    pub fn admission(&self) -> &ClaimAdmission {
        self.admission
    }

    /// Scientific color only when the claim has no waiver in its transitive
    /// derivation ancestry.
    #[must_use]
    pub fn scientific_color(&self) -> Option<&Color> {
        (self.admission.class == AdmissionClass::Scientific).then_some(&self.claim.color)
    }

    /// Raw declared color, explicitly not a scientific-admission assertion.
    #[must_use]
    pub fn declared_color_unverified(&self) -> &Color {
        &self.claim.color
    }

    /// Origin mechanism admitted by the policy recorded in the package receipt.
    #[must_use]
    pub fn origin(&self) -> &ClaimOrigin {
        &self.claim.origin
    }

    /// Composition receipt admitted by the derivation policy, when derived.
    #[must_use]
    pub fn derivation(&self) -> Option<&CompositionReceipt> {
        self.claim.receipt.as_ref()
    }

    /// Attached falsifier records, all authenticated by the recorded falsifier
    /// policy when the slice is non-empty.
    #[must_use]
    pub fn falsifiers(&self) -> &[FalsifierRecord] {
        &self.claim.falsifiers
    }

    /// Attached anchor declarations bound into the admitted claim subject.
    /// Dataset quality remains limited by the origin/derivation policy.
    #[must_use]
    pub fn anchors(&self) -> &[AnchorRecord] {
        &self.claim.anchors
    }
}

/// A structured verification failure.
#[derive(Debug, Clone, PartialEq)]
pub enum PackageError {
    /// Required reproducibility provenance is blank.
    IncompleteProvenance {
        /// What is missing (`"code_version"` or `"constellation_lock"`).
        missing: &'static str,
    },
    /// A machine identity is padded or uses a reserved placeholder token.
    InvalidIdentity {
        /// Claim identity when the field belongs to a claim.
        claim: Option<String>,
        /// Stable field path.
        field: &'static str,
        /// `"placeholder"` or `"surrounding-whitespace"`.
        reason: &'static str,
    },
    /// A claim id is blank or duplicates an earlier claim id.
    InvalidClaimId {
        /// The claim's position in the package.
        index: usize,
        /// The invalid id (blank for the blank-id case).
        id: String,
        /// Why it is invalid (`"blank"` or `"duplicate"`).
        reason: &'static str,
    },
    /// A claim has no meaningful human-readable assertion.
    InvalidClaimStatement {
        /// The claim id.
        claim: String,
        /// Why it is invalid (`"blank"` or `"placeholder"`).
        reason: &'static str,
    },
    /// A validated claim is missing part of its evidence.
    IncompleteValidatedClaim {
        /// The claim id.
        claim: String,
        /// What is missing (`"regime"` or `"dataset"`).
        missing: &'static str,
    },
    /// A verified claim's certificate interval is not a finite `[lo <= hi]`.
    IncompleteVerifiedClaim {
        /// The claim id.
        claim: String,
    },
    /// A validated claim has a malformed validity-domain axis.
    InvalidValidatedRegime {
        /// The claim id.
        claim: String,
        /// The malformed axis name (blank for a blank name).
        axis: String,
    },
    /// An estimated claim is missing its estimator identity.
    IncompleteEstimatedClaim {
        /// The claim id.
        claim: String,
        /// What is missing (`"estimator"`).
        missing: &'static str,
    },
    /// An estimated claim's dispersion is NaN or negative. Positive infinity
    /// is the lower-layer algebra's explicit no-quantitative-claim sentinel.
    InvalidEstimatedDispersion {
        /// The claim id.
        claim: String,
    },
    /// Finite claim magnitudes overflowed while deriving the package budget.
    MagnitudeOverflow {
        /// Claim at which the finite subtotal became non-finite.
        claim: String,
        /// Budget component (`"verified_width"` or `"estimated_dispersion"`).
        component: &'static str,
    },
    /// The in-memory package cannot fit the standalone checker's bounded
    /// transport envelope.
    TransportLimit {
        /// Field or container that exceeded its limit.
        what: String,
        /// Configured upper bound.
        limit: usize,
    },
    /// The declared format version is unsupported.
    UnsupportedFormat {
        /// The version found.
        found: u32,
    },
    /// A derived receipt names a different color-algebra semantics.
    UnsupportedColorAlgebra {
        /// Derived claim carrying the receipt.
        claim: String,
        /// Algebra version found.
        found: u32,
        /// Algebra version this checker executes.
        supported: u32,
    },
    /// A composition receipt does not re-derive the claimed color: the
    /// checker re-ran `compose` over the parents and got a different
    /// result — a forged or stale derivation (schema v3).
    ReceiptMismatch {
        /// The claim id.
        claim: String,
    },
    /// A receipt references a parent at or after the claim itself (the
    /// derivation DAG must point strictly backwards), or out of range.
    BadReceiptParent {
        /// The claim id.
        claim: String,
        /// The offending parent index.
        parent: usize,
    },
    /// A derivation receipt lacks a canonical proof-artifact address.
    InvalidDerivationArtifact {
        /// Derived claim.
        claim: String,
    },
    /// A derivation proof artifact could not be authenticated.
    DerivationRefused {
        /// Derived claim.
        claim: String,
        /// Why verification refused.
        why: &'static str,
        /// Atomic rejecting policy identity, absent for missing capability or
        /// callback panic.
        policy_fingerprint: Option<PolicyFingerprint>,
    },
    /// Schema v6: an origin whose fields fail shape validation.
    InvalidOrigin {
        /// The claim.
        claim: String,
        /// The field-level refusal.
        why: String,
    },
    /// Schema v6: an origin inconsistent with its claim's color class
    /// (raw colors, unrelated anchors, estimator mismatches, Derived
    /// without a receipt or a receipt without Derived).
    OriginMismatch {
        /// The claim.
        claim: String,
        /// The origin kind tag.
        origin: &'static str,
    },
    /// Schema v6: a source-certificate artifact could not be authenticated.
    SourceCertificateRefused {
        /// The claim.
        claim: String,
        /// Declared certificate producer.
        producer: String,
        /// Why verification refused.
        why: &'static str,
        /// Atomic rejecting policy identity, absent for missing capability or
        /// callback panic.
        policy_fingerprint: Option<PolicyFingerprint>,
    },
    /// Schema v6: an anchoring dataset could not be authenticated.
    AnchoredSourceRefused {
        /// The claim.
        claim: String,
        /// Declared dataset identity.
        dataset: String,
        /// Why verification refused.
        why: &'static str,
        /// Atomic rejecting policy identity, absent for missing capability or
        /// callback panic.
        policy_fingerprint: Option<PolicyFingerprint>,
    },
    /// Schema v6: a falsifier artifact could not be authenticated.
    FalsifierRefused {
        /// Target claim.
        claim: String,
        /// Falsifier record identity.
        falsifier: String,
        /// Why verification refused.
        why: &'static str,
        /// Atomic rejecting policy identity, absent for missing capability or
        /// callback panic.
        policy_fingerprint: Option<PolicyFingerprint>,
    },
    /// Schema v6: a waiver grant that is expired or that the injected
    /// verifier rejected (or no capability was injected at all).
    WaiverRefused {
        /// The claim.
        claim: String,
        /// The waiver id.
        waiver: String,
        /// Why.
        why: &'static str,
        /// Atomic rejecting policy identity, absent for missing capability,
        /// expiry, unavailable context, or callback panic.
        policy_fingerprint: Option<PolicyFingerprint>,
    },
    /// An installed verifier could not provide its policy identity.
    PolicyFingerprintRefused {
        /// Capability kind.
        capability: &'static str,
        /// Why the policy callback failed.
        why: &'static str,
        /// First policy identity observed for this capability kind.
        previous: PolicyFingerprint,
        /// Conflicting identity returned by a later decision.
        observed: PolicyFingerprint,
    },
    /// A supplied detached-signature verifier rejected or could not evaluate
    /// the signature.
    SignatureRefused {
        /// Why signature authentication refused.
        why: &'static str,
        /// Atomic rejecting policy identity, absent for callback panic.
        policy_fingerprint: Option<PolicyFingerprint>,
    },
    /// Detached signature bytes are blank, padded, placeholder, or contain
    /// control characters before policy evaluation.
    InvalidSignature {
        /// Stable shape refusal.
        why: &'static str,
    },
    /// Two claims reuse one waiver authorization identity.
    DuplicateWaiverId {
        /// Duplicated waiver id.
        waiver: String,
        /// Claim that first used it.
        first_claim: String,
        /// Later claim that reused it.
        duplicate_claim: String,
    },
    /// A waiver-MAC builder targeted a non-waiver claim or missing index.
    InvalidWaiverTarget {
        /// Requested claim index.
        index: usize,
    },
    /// A falsifier REFUTED this claim; a refuted claim cannot verify.
    RefutedClaim {
        /// The claim id.
        claim: String,
        /// The refuting falsifier.
        falsifier: String,
    },
    /// A falsifier record is not meaningful evidence: identities and outcome
    /// details must be non-blank and non-placeholder, and at least one
    /// adversarial attempt must have run.
    InvalidFalsifierRecord {
        /// The claim id.
        claim: String,
        /// Position of the malformed record within the claim.
        falsifier: usize,
        /// The invalid field (`"name"`, `"attempts"`, or `"detail"`).
        field: &'static str,
    },
    /// An anchoring-dataset record lacks a stable identity or a canonical
    /// content hash.
    InvalidAnchorRecord {
        /// The claim id.
        claim: String,
        /// Position of the malformed record within the claim.
        anchor: usize,
        /// The invalid field (`"dataset_id"` or `"content_hash"`).
        field: &'static str,
    },
}

impl core::fmt::Display for PackageError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(formatter, "package verification refused: {self:?}")
    }
}

impl core::error::Error for PackageError {}

/// The one format version this build understands. v2 (bead qmao.6.1)
/// added complete color payloads + the strict parser + root
/// recomputation; v3 (bead xfxq) added composition receipts (checker
/// re-runs the derivation), falsifier records (refuted claims fail),
/// and dataset anchors; v4 (bead 7uq9) replaces the 64-bit FNV-1a
/// content address with a domain-separated 32-byte BLAKE3 root
/// ([`ContentHash`]); v5 sealed claims with typed origins; v6 adds
/// capability-gated anchoring datasets and signatures, policy-bound
/// verification receipts, transitive waiver admission, and waiver-separated
/// magnitude accounting; v7 binds the color-algebra version into every
/// composition receipt. Earlier transports are refused by version.
pub const FORMAT_VERSION: u32 = 7;
const _: () = assert!(FORMAT_VERSION == fs_crosswalk::SUPPORTED_PACKAGE_FORMAT);

/// Domain-separated integrity digest for the independently distributed
/// checker's gate-context receipt. This is integrity-only, not a signature.
#[must_use]
pub fn hash_checker_decision(payload: &[u8]) -> ContentHash {
    hash_domain("fs-package:v7:checker-decision", payload)
}

fn verify_attached_records(claim: &Claim) -> Result<(), PackageError> {
    for (falsifier, record) in claim.falsifiers.iter().enumerate() {
        let field = if identity_reason(&record.name).is_some() {
            Some("name")
        } else if record.attempts == 0 {
            Some("attempts")
        } else if is_blank_or_placeholder(&record.detail) {
            Some("detail")
        } else if !is_canonical_content_hash(&record.artifact_hash) {
            Some("artifact_hash")
        } else {
            None
        };
        if let Some(field) = field {
            return Err(PackageError::InvalidFalsifierRecord {
                claim: claim.id.clone(),
                falsifier,
                field,
            });
        }
    }
    for (anchor, record) in claim.anchors.iter().enumerate() {
        let field = if identity_reason(&record.dataset_id).is_some() {
            Some("dataset_id")
        } else if !is_canonical_content_hash(&record.content_hash) {
            Some("content_hash")
        } else {
            None
        };
        if let Some(field) = field {
            return Err(PackageError::InvalidAnchorRecord {
                claim: claim.id.clone(),
                anchor,
                field,
            });
        }
    }
    Ok(())
}

fn verify_color_payload(claim: &Claim) -> Result<(), PackageError> {
    match &claim.color {
        Color::Verified { .. } => {}
        Color::Validated { regime, dataset } => {
            if regime.bounds().is_empty() {
                return Err(PackageError::IncompleteValidatedClaim {
                    claim: claim.id.clone(),
                    missing: "regime",
                });
            }
            if dataset.trim().is_empty() {
                return Err(PackageError::IncompleteValidatedClaim {
                    claim: claim.id.clone(),
                    missing: "dataset",
                });
            }
        }
        Color::Estimated { estimator, .. } => {
            if estimator.trim().is_empty() {
                return Err(PackageError::IncompleteEstimatedClaim {
                    claim: claim.id.clone(),
                    missing: "estimator",
                });
            }
        }
    }
    validate_color_payload(&claim.color).map_err(|error| match error {
        ColorPayloadError::InvalidIdentity { field, reason, .. } => PackageError::InvalidIdentity {
            claim: Some(claim.id.clone()),
            field: match field {
                "dataset" => "color.dataset",
                "estimator" => "color.estimator",
                _ => "color.regime.axis",
            },
            reason,
        },
        ColorPayloadError::InvalidVerifiedInterval { .. } => {
            PackageError::IncompleteVerifiedClaim {
                claim: claim.id.clone(),
            }
        }
        ColorPayloadError::InvalidValidatedRegime { axis, .. } => {
            PackageError::InvalidValidatedRegime {
                claim: claim.id.clone(),
                axis,
            }
        }
        ColorPayloadError::InvalidEstimatedDispersion { .. } => {
            PackageError::InvalidEstimatedDispersion {
                claim: claim.id.clone(),
            }
        }
    })
}

/// Equality for receipt re-derivation is an identity check, not a numerical
/// comparison. IEEE-754 signed zero is semantically distinct in canonical
/// evidence bytes and therefore must remain distinct here as well.
fn colors_bitwise_equal(left: &Color, right: &Color) -> bool {
    match (left, right) {
        (
            Color::Verified {
                lo: left_lo,
                hi: left_hi,
            },
            Color::Verified { lo, hi },
        ) => left_lo.to_bits() == lo.to_bits() && left_hi.to_bits() == hi.to_bits(),
        (
            Color::Validated {
                regime: left_regime,
                dataset: left_dataset,
            },
            Color::Validated { regime, dataset },
        ) => {
            left_dataset == dataset
                && left_regime.bounds().len() == regime.bounds().len()
                && left_regime.bounds().iter().zip(regime.bounds()).all(
                    |((left_axis, (left_lo, left_hi)), (axis, (lo, hi)))| {
                        left_axis == axis
                            && left_lo.to_bits() == lo.to_bits()
                            && left_hi.to_bits() == hi.to_bits()
                    },
                )
        }
        (
            Color::Estimated {
                estimator: left_estimator,
                dispersion: left_dispersion,
            },
            Color::Estimated {
                estimator,
                dispersion,
            },
        ) => left_estimator == estimator && left_dispersion.to_bits() == dispersion.to_bits(),
        _ => false,
    }
}

fn verify_origin_binding(claim: &Claim) -> Result<(), PackageError> {
    validate_origin_shape(&claim.id, &claim.origin, &is_canonical_content_hash).map_err(
        |error| PackageError::InvalidOrigin {
            claim: error.claim,
            why: error.why,
        },
    )?;
    let consistent = match (&claim.origin, &claim.color) {
        (ClaimOrigin::SourceCertificate { .. }, Color::Verified { .. })
        | (ClaimOrigin::Derived | ClaimOrigin::AuthenticatedWaiver(_), _) => true,
        (
            ClaimOrigin::AnchoredSource {
                dataset_id,
                content_hash,
            },
            Color::Validated { dataset, .. },
        ) => {
            dataset_id == dataset
                && claim.anchors.iter().any(|anchor| {
                    anchor.dataset_id == *dataset_id && anchor.content_hash == *content_hash
                })
        }
        (ClaimOrigin::EstimatedSource { estimator: from }, Color::Estimated { estimator, .. }) => {
            from == estimator
        }
        _ => false,
    };
    if !consistent || matches!(claim.origin, ClaimOrigin::Derived) != claim.receipt.is_some() {
        return Err(PackageError::OriginMismatch {
            claim: claim.id.clone(),
            origin: claim.origin.kind(),
        });
    }
    Ok(())
}

fn add_color_transport(
    index: usize,
    color: &Color,
    bytes: &mut usize,
    nodes: &mut usize,
) -> Result<(), PackageError> {
    match color {
        Color::Verified { .. } => {}
        Color::Validated { regime, dataset } => {
            check_transport_count("validated regime axes", regime.bounds().len())?;
            add_transport_nodes(nodes, regime.bounds().len().saturating_mul(3))?;
            add_transport_text(bytes, &format!("claims[{index}].dataset"), dataset)?;
            for axis in regime.bounds().keys() {
                add_transport_text(bytes, &format!("claims[{index}].regime axis"), axis)?;
                add_transport_bytes(bytes, 64)?;
            }
        }
        Color::Estimated { estimator, .. } => {
            add_transport_text(bytes, &format!("claims[{index}].estimator"), estimator)?;
        }
    }
    Ok(())
}

fn add_record_transport(
    claim: &Claim,
    bytes: &mut usize,
    nodes: &mut usize,
) -> Result<(), PackageError> {
    if let Some(receipt) = &claim.receipt {
        check_transport_count("receipt parents", receipt.parents.len())?;
        add_transport_nodes(nodes, receipt.parents.len().saturating_add(4))?;
        add_transport_bytes(
            bytes,
            32usize
                .saturating_mul(receipt.parents.len())
                .saturating_add(80),
        )?;
        add_transport_text(bytes, "receipt.artifact_hash", &receipt.artifact_hash)?;
    }
    check_transport_count("falsifiers", claim.falsifiers.len())?;
    check_transport_count("anchors", claim.anchors.len())?;
    add_transport_nodes(
        nodes,
        claim
            .falsifiers
            .len()
            .saturating_mul(6)
            .saturating_add(claim.anchors.len().saturating_mul(4)),
    )?;
    add_transport_bytes(
        bytes,
        claim
            .falsifiers
            .len()
            .saturating_mul(160)
            .saturating_add(claim.anchors.len().saturating_mul(128)),
    )?;
    for falsifier in &claim.falsifiers {
        add_transport_text(bytes, "falsifier.name", &falsifier.name)?;
        add_transport_text(bytes, "falsifier.detail", &falsifier.detail)?;
        add_transport_text(bytes, "falsifier.artifact_hash", &falsifier.artifact_hash)?;
    }
    for anchor in &claim.anchors {
        add_transport_text(bytes, "anchor.dataset_id", &anchor.dataset_id)?;
        add_transport_text(bytes, "anchor.content_hash", &anchor.content_hash)?;
    }
    Ok(())
}

fn add_origin_transport(
    origin: &ClaimOrigin,
    bytes: &mut usize,
    nodes: &mut usize,
) -> Result<(), PackageError> {
    match origin {
        ClaimOrigin::SourceCertificate {
            producer,
            certificate_hash,
        } => {
            add_transport_nodes(nodes, 3)?;
            add_transport_text(bytes, "origin.producer", producer)?;
            add_transport_text(bytes, "origin.certificate_hash", certificate_hash)?;
        }
        ClaimOrigin::AnchoredSource {
            dataset_id,
            content_hash,
        } => {
            add_transport_nodes(nodes, 3)?;
            add_transport_text(bytes, "origin.dataset_id", dataset_id)?;
            add_transport_text(bytes, "origin.content_hash", content_hash)?;
        }
        ClaimOrigin::EstimatedSource { estimator } => {
            add_transport_nodes(nodes, 2)?;
            add_transport_text(bytes, "origin.estimator", estimator)?;
        }
        ClaimOrigin::Derived => add_transport_nodes(nodes, 1)?,
        ClaimOrigin::AuthenticatedWaiver(grant) => {
            add_transport_nodes(nodes, 4)?;
            add_transport_bytes(bytes, 32)?;
            add_transport_text(bytes, "origin.waiver_id", &grant.waiver_id)?;
            add_transport_text(bytes, "origin.mac", &grant.mac)?;
        }
    }
    Ok(())
}

fn add_claim_transport(
    index: usize,
    claim: &Claim,
    bytes: &mut usize,
    nodes: &mut usize,
) -> Result<(), PackageError> {
    add_transport_bytes(bytes, 256)?;
    add_transport_nodes(nodes, 11)?;
    add_transport_text(bytes, &format!("claims[{index}].id"), &claim.id)?;
    add_transport_text(
        bytes,
        &format!("claims[{index}].statement"),
        &claim.statement,
    )?;
    add_color_transport(index, &claim.color, bytes, nodes)?;
    add_record_transport(claim, bytes, nodes)?;
    add_origin_transport(&claim.origin, bytes, nodes)
}

impl EvidencePackage {
    /// An empty package at the current format version.
    #[must_use]
    pub fn new(provenance: Provenance) -> EvidencePackage {
        EvidencePackage {
            format_version: FORMAT_VERSION,
            claims: Vec::new(),
            provenance,
            signature: None,
        }
    }

    /// Add a claim (builder style).
    #[must_use]
    pub fn with_claim(mut self, claim: Claim) -> EvidencePackage {
        self.claims.push(claim);
        self
    }

    /// Raw claim declarations before external origin/falsifier/signature
    /// admission. Callers needing scientific values must use
    /// [`VerifiedPackage::admitted_claims`].
    #[must_use]
    pub fn declared_claims_unverified(&self) -> &[Claim] {
        &self.claims
    }

    /// Whether raw declarations are safe to scan for bounded diagnostics.
    ///
    /// This checks the complete transport envelope and structural semantics but
    /// performs no external authentication and grants no scientific authority.
    /// It exists so standalone preflight tools can inventory declaration-level
    /// blockers after a fail-closed capability refusal without rescanning an
    /// oversized or malformed builder.
    #[must_use]
    pub fn is_structurally_inspectable_unverified(&self) -> bool {
        self.verify_structural().is_ok()
    }

    /// Attach a detached signature (builder style).
    #[must_use]
    pub fn signed(mut self, signature: impl Into<String>) -> EvidencePackage {
        self.signature = Some(signature.into());
        self
    }

    /// Compute the BLAKE3 Merkle content address after enforcing the standalone
    /// transport byte, node, string, and container limits. Every leaf and
    /// internal node is domain-separated under `fs-package:v7:...`; detached
    /// signatures are excluded so signing does not change the address.
    ///
    /// # Errors
    /// [`PackageError::TransportLimit`] when an in-memory builder exceeds the
    /// bounded package envelope.
    pub fn try_merkle_root(&self) -> Result<ContentHash, PackageError> {
        self.verify_transport_limits()?;
        Ok(self.merkle_root_unchecked())
    }

    fn merkle_root_unchecked(&self) -> ContentHash {
        let mut level: Vec<ContentHash> = Vec::with_capacity(self.claims.len() + 1);
        level.push(hash_domain(
            "fs-package:v7:header",
            self.package_header().as_bytes(),
        ));
        level.extend(
            self.claims
                .iter()
                .map(|c| hash_domain("fs-package:v7:claim", c.canonical().as_bytes())),
        );
        while level.len() > 1 {
            let mut next = Vec::with_capacity(level.len().div_ceil(2));
            for pair in level.chunks(2) {
                match pair {
                    [a, b] => next.push(combine(a, b)),
                    [a] => next.push(*a), // odd node carries up
                    _ => {}
                }
            }
            level = next;
        }
        match level.as_slice() {
            [root] => *root,
            [] => hash_domain("fs-package:v7:empty-internal-level", b""),
            _ => hash_domain("fs-package:v7:invalid-internal-level", b""),
        }
    }

    fn package_header(&self) -> String {
        use core::fmt::Write as _;
        let mut out = String::from("package|");
        let _ = write!(
            out,
            "format:{}|claims:{}|",
            self.format_version,
            self.claims.len()
        );
        push_atom(&mut out, &self.provenance.code_version);
        push_atom(&mut out, &self.provenance.constellation_lock);
        out
    }

    fn authorization_context(&self) -> ContentHash {
        let mut canonical = self.package_header();
        for claim in &self.claims {
            push_atom(&mut canonical, &claim.authorization_canonical());
        }
        hash_domain("fs-package:v7:authorization-context", canonical.as_bytes())
    }

    /// Stable, domain-separated bytes authenticated by a waiver MAC at
    /// `claim_index`. The context binds package provenance, ordered claims,
    /// target index, waiver id, and expiry. Detached signatures and every
    /// waiver MAC are intentionally excluded.
    ///
    /// # Errors
    /// [`PackageError::TransportLimit`] when the package exceeds the bounded
    /// authorization envelope, or [`PackageError::InvalidWaiverTarget`] when
    /// the index does not name a waiver-origin claim.
    pub fn waiver_message(&self, claim_index: usize) -> Result<Vec<u8>, PackageError> {
        self.verify_transport_limits()?;
        self.waiver_message_with_context(claim_index, self.authorization_context())
            .ok_or(PackageError::InvalidWaiverTarget { index: claim_index })
    }

    fn waiver_message_with_context(
        &self,
        claim_index: usize,
        authorization_context: ContentHash,
    ) -> Option<Vec<u8>> {
        use core::fmt::Write as _;

        let claim = self.claims.get(claim_index)?;
        let ClaimOrigin::AuthenticatedWaiver(grant) = &claim.origin else {
            return None;
        };
        let mut message = String::from("fs-package:v7:waiver-authorization|");
        push_atom(&mut message, &authorization_context.to_hex());
        let _ = write!(message, "claim-index:{claim_index}|");
        push_atom(&mut message, &claim.canonical_body());
        push_atom(&mut message, &grant.waiver_id);
        let _ = write!(message, "expiry-day:{}|", grant.expiry_day);
        Some(message.into_bytes())
    }

    /// Install the final authenticator for one waiver-origin claim. The
    /// corresponding [`EvidencePackage::waiver_message`] is stable before and
    /// after this operation because MAC bytes are excluded from its context.
    ///
    /// # Errors
    /// [`PackageError::InvalidWaiverTarget`] when `claim_index` is absent or
    /// names a claim with another origin kind.
    pub fn with_waiver_mac(
        mut self,
        claim_index: usize,
        mac: impl Into<String>,
    ) -> Result<EvidencePackage, PackageError> {
        let Some(claim) = self.claims.get_mut(claim_index) else {
            return Err(PackageError::InvalidWaiverTarget { index: claim_index });
        };
        let ClaimOrigin::AuthenticatedWaiver(grant) = &mut claim.origin else {
            return Err(PackageError::InvalidWaiverTarget { index: claim_index });
        };
        grant.mac = mac.into();
        Ok(self)
    }

    fn admission_decisions(&self) -> (Vec<ClaimAdmission>, Vec<ReceiptWaiver>) {
        let mut waiver_dependent = Vec::with_capacity(self.claims.len());
        let mut decisions = Vec::with_capacity(self.claims.len());
        let mut waiver_registry = Vec::new();
        for (claim_index, claim) in self.claims.iter().enumerate() {
            let mut direct_waiver = None;
            let mut waiver_parents = Vec::new();
            let class = match &claim.origin {
                ClaimOrigin::AuthenticatedWaiver(grant) => {
                    let registry_index = waiver_registry.len();
                    waiver_registry.push(ReceiptWaiver {
                        registry_index,
                        claim_index,
                        waiver_id: grant.waiver_id.clone(),
                    });
                    direct_waiver = Some(registry_index);
                    AdmissionClass::WaiverDependent
                }
                ClaimOrigin::Derived => {
                    if let Some(receipt) = &claim.receipt {
                        for &parent in &receipt.parents {
                            if waiver_dependent.get(parent).copied().unwrap_or(true) {
                                waiver_parents.push(parent);
                            }
                        }
                    }
                    if waiver_parents.is_empty() {
                        AdmissionClass::Scientific
                    } else {
                        AdmissionClass::WaiverDependent
                    }
                }
                ClaimOrigin::SourceCertificate { .. }
                | ClaimOrigin::AnchoredSource { .. }
                | ClaimOrigin::EstimatedSource { .. } => AdmissionClass::Scientific,
            };
            waiver_dependent.push(class == AdmissionClass::WaiverDependent);
            decisions.push(ClaimAdmission {
                claim_index,
                claim_id: claim.id.clone(),
                origin_kind: match &claim.origin {
                    ClaimOrigin::SourceCertificate { .. } => AdmissionOriginKind::SourceCertificate,
                    ClaimOrigin::AnchoredSource { .. } => AdmissionOriginKind::AnchoredSource,
                    ClaimOrigin::EstimatedSource { .. } => AdmissionOriginKind::EstimatedSource,
                    ClaimOrigin::Derived => AdmissionOriginKind::Derived,
                    ClaimOrigin::AuthenticatedWaiver(_) => AdmissionOriginKind::AuthenticatedWaiver,
                },
                class,
                direct_waiver,
                waiver_parents,
            });
        }
        (decisions, waiver_registry)
    }

    fn admitted_color_breakdown(&self, admissions: &[ClaimAdmission]) -> ColorBreakdown {
        let mut b = ColorBreakdown::default();
        for (c, admission) in self.claims.iter().zip(admissions) {
            if admission.class == AdmissionClass::WaiverDependent {
                b.waived += 1;
                continue;
            }
            match c.color.rank() {
                ColorRank::Verified => b.verified += 1,
                ColorRank::Validated => b.validated += 1,
                ColorRank::Estimated => b.estimated += 1,
            }
        }
        b
    }

    /// The by-color budget pie, available only after fail-closed verification
    /// with no external capabilities.
    ///
    /// # Errors
    /// Any refusal from [`EvidencePackage::verify`].
    pub fn color_breakdown(&self) -> Result<ColorBreakdown, PackageError> {
        self.verify().map(|report| report.breakdown)
    }

    /// The by-color budget pie after verification with explicit capabilities.
    ///
    /// # Errors
    /// Any refusal from [`EvidencePackage::verify_with`].
    pub fn color_breakdown_with(
        &self,
        capabilities: &VerificationCapabilities<'_>,
    ) -> Result<ColorBreakdown, PackageError> {
        self.verify_with(capabilities)
            .map(|report| report.breakdown)
    }

    fn verify_claim(&self, index: usize, claim: &Claim) -> Result<(), PackageError> {
        // Schema-v3 semantic re-verification (solver-free): refuted falsifiers
        // fail and composition receipts are independently re-derived.
        verify_attached_records(claim)?;
        if let Some(fr) = claim.falsifiers.iter().find(|f| f.refuted) {
            return Err(PackageError::RefutedClaim {
                claim: claim.id.clone(),
                falsifier: fr.name.clone(),
            });
        }
        if let Some(receipt) = &claim.receipt {
            if receipt.color_algebra_version != fs_evidence::COLOR_ALGEBRA_VERSION {
                return Err(PackageError::UnsupportedColorAlgebra {
                    claim: claim.id.clone(),
                    found: receipt.color_algebra_version,
                    supported: fs_evidence::COLOR_ALGEBRA_VERSION,
                });
            }
            if !is_canonical_content_hash(&receipt.artifact_hash) {
                return Err(PackageError::InvalidDerivationArtifact {
                    claim: claim.id.clone(),
                });
            }
            let mut derived: Option<Color> = None;
            for &parent in &receipt.parents {
                if parent >= index {
                    return Err(PackageError::BadReceiptParent {
                        claim: claim.id.clone(),
                        parent,
                    });
                }
                let parent_color = &self.claims[parent].color;
                derived = Some(match derived {
                    None => parent_color.clone(),
                    Some(current) => compose(&current, parent_color, receipt.op),
                });
            }
            if !matches!(derived, Some(ref color) if colors_bitwise_equal(color, &claim.color)) {
                return Err(PackageError::ReceiptMismatch {
                    claim: claim.id.clone(),
                });
            }
        }
        verify_color_payload(claim)?;
        verify_origin_binding(claim)
    }

    /// Re-verify structural semantics and every capability-gated origin.
    /// Source-certificate hashes are artifact addresses, not proof by
    /// themselves; waiver origins likewise require an authenticator plus an
    /// explicit date. Missing capabilities always fail closed.
    ///
    /// # Errors
    /// Any structural [`PackageError`],
    /// [`PackageError::SourceCertificateRefused`], or
    /// [`PackageError::WaiverRefused`].
    #[allow(clippy::too_many_lines)] // one ordered, fail-closed transcript over every claim authority
    pub fn verify_with(
        &self,
        capabilities: &VerificationCapabilities<'_>,
    ) -> Result<PackageReport, PackageError> {
        self.verify_structural()?;
        let merkle_root = self.merkle_root_unchecked();
        let mut policy_fingerprints = VerificationPolicyFingerprints::default();
        // The authorization context serializes and hashes the whole package.
        // Compute it once so W waiver claims remain O(package size + W), not
        // O(W * package size).
        let waiver_context = (self.waiver_claims() > 0).then(|| self.authorization_context());
        let claim_hashes: Vec<ContentHash> = self
            .claims
            .iter()
            .map(Claim::declared_content_hash_unverified)
            .collect();
        let claim_subject_hashes: Vec<ContentHash> = self
            .claims
            .iter()
            .map(Claim::declared_verification_subject_hash_unverified)
            .collect();
        for (claim_index, claim) in self.claims.iter().enumerate() {
            if let (ClaimOrigin::Derived, Some(receipt)) = (&claim.origin, &claim.receipt) {
                let Some(verifier) = capabilities.derivations else {
                    return Err(PackageError::DerivationRefused {
                        claim: claim.id.clone(),
                        why: "derivation capability missing",
                        policy_fingerprint: None,
                    });
                };
                let Some(artifact_hash) = ContentHash::from_hex(&receipt.artifact_hash) else {
                    return Err(PackageError::InvalidDerivationArtifact {
                        claim: claim.id.clone(),
                    });
                };
                let parent_claim_hashes: Vec<ContentHash> = receipt
                    .parents
                    .iter()
                    .map(|&parent| claim_hashes[parent])
                    .collect();
                let request = DerivationRequest {
                    package_provenance: &self.provenance,
                    package_root: merkle_root,
                    claim_index,
                    claim_id: &claim.id,
                    statement: &claim.statement,
                    color: &claim.color,
                    child_subject_hash: claim_subject_hashes[claim_index],
                    anchors: &claim.anchors,
                    op: receipt.op,
                    parent_indices: &receipt.parents,
                    parent_claim_hashes: &parent_claim_hashes,
                    artifact_hash,
                };
                let decision = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                    verifier.verify(&request)
                }))
                .map_err(|_| PackageError::DerivationRefused {
                    claim: claim.id.clone(),
                    why: "verifier callback panicked",
                    policy_fingerprint: None,
                })?;
                bind_policy_fingerprint(
                    &mut policy_fingerprints.derivations,
                    decision.policy_fingerprint(),
                    "derivations",
                )?;
                if !decision.accepted() {
                    return Err(PackageError::DerivationRefused {
                        claim: claim.id.clone(),
                        why: "rejected by the injected verifier",
                        policy_fingerprint: Some(decision.policy_fingerprint()),
                    });
                }
            }
            self.verify_validated_anchors(
                claim_index,
                claim,
                capabilities,
                &mut policy_fingerprints,
            )?;
            match (&claim.origin, &claim.color) {
                (
                    ClaimOrigin::SourceCertificate {
                        producer,
                        certificate_hash,
                    },
                    Color::Verified { lo, hi },
                ) => {
                    let Some(verifier) = capabilities.source_certificates else {
                        return Err(PackageError::SourceCertificateRefused {
                            claim: claim.id.clone(),
                            producer: producer.clone(),
                            why: "source-certificate capability missing",
                            policy_fingerprint: None,
                        });
                    };
                    let Some(certificate_hash) = ContentHash::from_hex(certificate_hash) else {
                        return Err(PackageError::InvalidOrigin {
                            claim: claim.id.clone(),
                            why: "source-certificate hash is not canonical".to_string(),
                        });
                    };
                    let request = SourceCertificateRequest {
                        package_provenance: &self.provenance,
                        claim_index,
                        claim_id: &claim.id,
                        statement: &claim.statement,
                        lo: *lo,
                        hi: *hi,
                        producer,
                        certificate_hash,
                    };
                    let decision = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                        verifier.verify(&request)
                    }))
                    .map_err(|_| PackageError::SourceCertificateRefused {
                        claim: claim.id.clone(),
                        producer: producer.clone(),
                        why: "verifier callback panicked",
                        policy_fingerprint: None,
                    })?;
                    bind_policy_fingerprint(
                        &mut policy_fingerprints.source_certificates,
                        decision.policy_fingerprint(),
                        "source-certificates",
                    )?;
                    if !decision.accepted() {
                        return Err(PackageError::SourceCertificateRefused {
                            claim: claim.id.clone(),
                            producer: producer.clone(),
                            why: "rejected by the injected verifier",
                            policy_fingerprint: Some(decision.policy_fingerprint()),
                        });
                    }
                }
                (ClaimOrigin::AuthenticatedWaiver(grant), _) => {
                    let Some(waivers) = capabilities.waivers else {
                        return Err(PackageError::WaiverRefused {
                            claim: claim.id.clone(),
                            waiver: grant.waiver_id.clone(),
                            why: "waiver capability missing",
                            policy_fingerprint: None,
                        });
                    };
                    if grant.expiry_day < waivers.today_day {
                        return Err(PackageError::WaiverRefused {
                            claim: claim.id.clone(),
                            waiver: grant.waiver_id.clone(),
                            why: "expired",
                            policy_fingerprint: None,
                        });
                    }
                    let Some(message) = waiver_context
                        .and_then(|context| self.waiver_message_with_context(claim_index, context))
                    else {
                        return Err(PackageError::WaiverRefused {
                            claim: claim.id.clone(),
                            waiver: grant.waiver_id.clone(),
                            why: "authorization message unavailable",
                            policy_fingerprint: None,
                        });
                    };
                    let decision = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                        waivers.verifier.verify(&grant.mac, &message)
                    }))
                    .map_err(|_| PackageError::WaiverRefused {
                        claim: claim.id.clone(),
                        waiver: grant.waiver_id.clone(),
                        why: "verifier callback panicked",
                        policy_fingerprint: None,
                    })?;
                    bind_policy_fingerprint(
                        &mut policy_fingerprints.waivers,
                        decision.policy_fingerprint(),
                        "waivers",
                    )?;
                    if !decision.accepted() {
                        return Err(PackageError::WaiverRefused {
                            claim: claim.id.clone(),
                            waiver: grant.waiver_id.clone(),
                            why: "rejected by the injected verifier",
                            policy_fingerprint: Some(decision.policy_fingerprint()),
                        });
                    }
                }
                _ => {}
            }
            self.verify_falsifiers(
                claim_index,
                claim,
                claim_subject_hashes[claim_index],
                merkle_root,
                capabilities,
                &mut policy_fingerprints,
            )?;
        }
        let (admissions, waiver_registry) = self.admission_decisions();
        let breakdown = self.admitted_color_breakdown(&admissions);
        let magnitude_budget = self.magnitude_budget_from(&admissions);
        let waiver_day = policy_fingerprints
            .waivers
            .and_then(|_| capabilities.waivers.map(|waiver| waiver.today_day));
        let admission_context = release_admission_context_hash(
            merkle_root,
            policy_fingerprints,
            waiver_day,
            &admissions,
            &waiver_registry,
        );
        let (signature, signature_policy) =
            self.verify_signature(merkle_root, admission_context, capabilities)?;
        policy_fingerprints.signatures = signature_policy;
        let receipt_hash = verification_receipt_hash(
            merkle_root,
            policy_fingerprints,
            waiver_day,
            &signature,
            &admissions,
            &waiver_registry,
        );
        Ok(PackageReport {
            merkle_root,
            breakdown,
            claims: self.claims.len(),
            magnitude_budget,
            receipt: VerificationReceipt {
                package_root: merkle_root,
                policy_fingerprints,
                waiver_day,
                signature,
                admissions,
                waiver_registry,
                receipt_hash,
            },
        })
    }

    fn verify_validated_anchors(
        &self,
        claim_index: usize,
        claim: &Claim,
        capabilities: &VerificationCapabilities<'_>,
        policies: &mut VerificationPolicyFingerprints,
    ) -> Result<(), PackageError> {
        let Color::Validated { regime, dataset } = &claim.color else {
            return Ok(());
        };
        match &claim.origin {
            ClaimOrigin::AnchoredSource {
                dataset_id,
                content_hash,
            } => self.verify_validated_anchor(
                claim_index,
                claim,
                regime,
                dataset_id,
                content_hash,
                capabilities,
                policies,
            ),
            ClaimOrigin::Derived => {
                let mut matching = claim
                    .anchors
                    .iter()
                    .filter(|anchor| anchor.dataset_id == *dataset)
                    .peekable();
                if matching.peek().is_none() {
                    return Err(PackageError::AnchoredSourceRefused {
                        claim: claim.id.clone(),
                        dataset: dataset.clone(),
                        why: "matching validated anchor missing",
                        policy_fingerprint: None,
                    });
                }
                for anchor in matching {
                    self.verify_validated_anchor(
                        claim_index,
                        claim,
                        regime,
                        &anchor.dataset_id,
                        &anchor.content_hash,
                        capabilities,
                        policies,
                    )?;
                }
                Ok(())
            }
            // An administrative waiver never turns its declared color into
            // scientific validation authority. The other alternatives are
            // rejected as color/origin mismatches by structural verification.
            ClaimOrigin::AuthenticatedWaiver(_)
            | ClaimOrigin::SourceCertificate { .. }
            | ClaimOrigin::EstimatedSource { .. } => Ok(()),
        }
    }

    #[allow(clippy::too_many_arguments)] // exact typed anchor subject; no ambient authority inputs
    fn verify_validated_anchor(
        &self,
        claim_index: usize,
        claim: &Claim,
        regime: &fs_evidence::ValidityDomain,
        dataset_id: &str,
        content_hash: &str,
        capabilities: &VerificationCapabilities<'_>,
        policies: &mut VerificationPolicyFingerprints,
    ) -> Result<(), PackageError> {
        let Some(verifier) = capabilities.anchored_sources else {
            return Err(PackageError::AnchoredSourceRefused {
                claim: claim.id.clone(),
                dataset: dataset_id.to_string(),
                why: "anchored-source capability missing",
                policy_fingerprint: None,
            });
        };
        let Some(content_hash) = ContentHash::from_hex(content_hash) else {
            return Err(PackageError::InvalidOrigin {
                claim: claim.id.clone(),
                why: "validated anchor hash is not canonical".to_string(),
            });
        };
        let request = AnchoredSourceRequest {
            package_provenance: &self.provenance,
            claim_index,
            claim_id: &claim.id,
            statement: &claim.statement,
            regime,
            dataset_id,
            content_hash,
        };
        let decision =
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| verifier.verify(&request)))
                .map_err(|_| PackageError::AnchoredSourceRefused {
                    claim: claim.id.clone(),
                    dataset: dataset_id.to_string(),
                    why: "verifier callback panicked",
                    policy_fingerprint: None,
                })?;
        bind_policy_fingerprint(
            &mut policies.anchored_sources,
            decision.policy_fingerprint(),
            "anchored-sources",
        )?;
        if !decision.accepted() {
            return Err(PackageError::AnchoredSourceRefused {
                claim: claim.id.clone(),
                dataset: dataset_id.to_string(),
                why: "rejected by the injected verifier",
                policy_fingerprint: Some(decision.policy_fingerprint()),
            });
        }
        Ok(())
    }

    fn verify_falsifiers(
        &self,
        claim_index: usize,
        claim: &Claim,
        claim_subject_hash: ContentHash,
        merkle_root: ContentHash,
        capabilities: &VerificationCapabilities<'_>,
        policies: &mut VerificationPolicyFingerprints,
    ) -> Result<(), PackageError> {
        for (falsifier_index, falsifier) in claim.falsifiers.iter().enumerate() {
            let Some(verifier) = capabilities.falsifiers else {
                return Err(PackageError::FalsifierRefused {
                    claim: claim.id.clone(),
                    falsifier: falsifier.name.clone(),
                    why: "falsifier capability missing",
                    policy_fingerprint: None,
                });
            };
            let Some(artifact_hash) = ContentHash::from_hex(&falsifier.artifact_hash) else {
                return Err(PackageError::InvalidFalsifierRecord {
                    claim: claim.id.clone(),
                    falsifier: falsifier_index,
                    field: "artifact_hash",
                });
            };
            let request = FalsifierRequest {
                package_provenance: &self.provenance,
                package_root: merkle_root,
                claim_index,
                claim_id: &claim.id,
                statement: &claim.statement,
                color: &claim.color,
                origin: &claim.origin,
                claim_subject_hash,
                falsifier_index,
                name: &falsifier.name,
                attempts: falsifier.attempts,
                refuted: falsifier.refuted,
                detail: &falsifier.detail,
                artifact_hash,
            };
            let decision = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                verifier.verify(&request)
            }))
            .map_err(|_| PackageError::FalsifierRefused {
                claim: claim.id.clone(),
                falsifier: falsifier.name.clone(),
                why: "verifier callback panicked",
                policy_fingerprint: None,
            })?;
            bind_policy_fingerprint(
                &mut policies.falsifiers,
                decision.policy_fingerprint(),
                "falsifiers",
            )?;
            if !decision.accepted() {
                return Err(PackageError::FalsifierRefused {
                    claim: claim.id.clone(),
                    falsifier: falsifier.name.clone(),
                    why: "rejected by the injected verifier",
                    policy_fingerprint: Some(decision.policy_fingerprint()),
                });
            }
        }
        Ok(())
    }

    fn verify_signature(
        &self,
        merkle_root: ContentHash,
        admission_context: ContentHash,
        capabilities: &VerificationCapabilities<'_>,
    ) -> Result<(SignatureStatus, Option<PolicyFingerprint>), PackageError> {
        let Some(signature) = &self.signature else {
            return Ok((SignatureStatus::Unsigned, None));
        };
        let Some(verification) = capabilities.signatures else {
            return Ok((SignatureStatus::Unverified(signature.clone()), None));
        };
        let purpose = match verification.intent {
            SignatureIntent::PackageRootAttestation => SignaturePurpose::PackageRootAttestation,
            SignatureIntent::ReleaseApproval {
                checker_protocol,
                expected_root,
            } => {
                if expected_root != merkle_root {
                    return Err(PackageError::SignatureRefused {
                        why: "release-approval purpose names a different package root",
                        policy_fingerprint: None,
                    });
                }
                SignaturePurpose::ReleaseApproval {
                    checker_protocol,
                    expected_root,
                    admission_context,
                }
            }
        };
        let request = SignatureRequest {
            package_root: merkle_root,
            signature,
            purpose,
        };
        let decision = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            verification.verifier.verify(&request)
        }))
        .map_err(|_| PackageError::SignatureRefused {
            why: "verifier callback panicked",
            policy_fingerprint: None,
        })?;
        if !decision.accepted() {
            return Err(PackageError::SignatureRefused {
                why: "rejected by the injected verifier",
                policy_fingerprint: Some(decision.policy_fingerprint()),
            });
        }
        Ok((
            SignatureStatus::Authenticated(AuthenticatedSignature {
                signature: signature.clone(),
                purpose,
            }),
            Some(decision.policy_fingerprint()),
        ))
    }

    /// The number of waiver-origin claims (a checker without an injected
    /// waiver capability must fail closed when this is non-zero).
    #[must_use]
    pub fn waiver_claims(&self) -> usize {
        self.claims
            .iter()
            .filter(|c| matches!(c.origin, ClaimOrigin::AuthenticatedWaiver(_)))
            .count()
    }

    /// Re-verify with NO external trust capabilities. Only an empty package or
    /// ungated Estimated-source claims without falsifier records can pass.
    /// Source certificates, anchored sources, derivations, waivers, and attached
    /// artifact records all fail closed without their exact capabilities.
    ///
    /// # Errors
    /// [`PackageError`] on an unsupported format or an incomplete claim.
    pub fn verify(&self) -> Result<PackageReport, PackageError> {
        self.verify_with(&VerificationCapabilities::deny_all())
    }

    /// Consume an in-memory package and retain it with the exact report/receipt
    /// that admitted it, avoiding a raw package/report split.
    ///
    /// # Errors
    /// Any refusal from [`EvidencePackage::verify_with`].
    pub fn into_verified_with(
        self,
        capabilities: &VerificationCapabilities<'_>,
    ) -> Result<VerifiedPackage, PackageError> {
        let report = self.verify_with(capabilities)?;
        Ok(VerifiedPackage {
            package: self,
            report,
        })
    }

    /// Consume a deny-all-admissible package into its receipt-bearing view.
    ///
    /// # Errors
    /// Any refusal from [`EvidencePackage::verify`].
    pub fn into_verified(self) -> Result<VerifiedPackage, PackageError> {
        self.into_verified_with(&VerificationCapabilities::deny_all())
    }

    fn verify_structural(&self) -> Result<(), PackageError> {
        if self.format_version != FORMAT_VERSION {
            return Err(PackageError::UnsupportedFormat {
                found: self.format_version,
            });
        }
        self.verify_transport_limits()?;
        if let Some(reason) = identity_reason(&self.provenance.code_version) {
            if reason == "blank" {
                return Err(PackageError::IncompleteProvenance {
                    missing: "code_version",
                });
            }
            return Err(PackageError::InvalidIdentity {
                claim: None,
                field: "provenance.code_version",
                reason,
            });
        }
        if let Some(reason) = identity_reason(&self.provenance.constellation_lock) {
            if reason == "blank" {
                return Err(PackageError::IncompleteProvenance {
                    missing: "constellation_lock",
                });
            }
            return Err(PackageError::InvalidIdentity {
                claim: None,
                field: "provenance.constellation_lock",
                reason,
            });
        }
        if let Some(signature) = &self.signature {
            let why = if signature.trim().is_empty() {
                Some("blank")
            } else if signature.trim() != signature {
                Some("surrounding-whitespace")
            } else if signature.chars().any(char::is_control) {
                Some("control-character")
            } else if is_placeholder(signature) {
                Some("placeholder")
            } else {
                None
            };
            if let Some(why) = why {
                return Err(PackageError::InvalidSignature { why });
            }
        }
        let mut claim_ids = std::collections::BTreeSet::new();
        let mut waiver_ids = std::collections::BTreeMap::new();
        for (index, c) in self.claims.iter().enumerate() {
            if let Some(reason) = identity_reason(&c.id) {
                return Err(PackageError::InvalidClaimId {
                    index,
                    id: c.id.clone(),
                    reason,
                });
            }
            if !claim_ids.insert(c.id.as_str()) {
                return Err(PackageError::InvalidClaimId {
                    index,
                    id: c.id.clone(),
                    reason: "duplicate",
                });
            }
            let statement = c.statement.trim();
            if statement.is_empty() {
                return Err(PackageError::InvalidClaimStatement {
                    claim: c.id.clone(),
                    reason: "blank",
                });
            }
            if is_placeholder(statement) {
                return Err(PackageError::InvalidClaimStatement {
                    claim: c.id.clone(),
                    reason: "placeholder",
                });
            }
            self.verify_claim(index, c)?;
            if let ClaimOrigin::AuthenticatedWaiver(grant) = &c.origin
                && let Some(first_claim) =
                    waiver_ids.insert(grant.waiver_id.as_str(), c.id.as_str())
            {
                return Err(PackageError::DuplicateWaiverId {
                    waiver: grant.waiver_id.clone(),
                    first_claim: first_claim.to_string(),
                    duplicate_claim: c.id.clone(),
                });
            }
        }
        self.verify_finite_magnitude_sums()?;
        Ok(())
    }

    fn verify_transport_limits(&self) -> Result<(), PackageError> {
        self.transport_usage().map(|_| ())
    }

    fn transport_usage(&self) -> Result<(usize, usize), PackageError> {
        check_transport_count("claims", self.claims.len())?;
        let mut bytes = 512usize;
        let mut nodes = 13usize;
        add_transport_text(
            &mut bytes,
            "provenance.code_version",
            &self.provenance.code_version,
        )?;
        add_transport_text(
            &mut bytes,
            "provenance.constellation_lock",
            &self.provenance.constellation_lock,
        )?;
        if let Some(signature) = &self.signature {
            add_transport_text(&mut bytes, "signature", signature)?;
        }
        for (index, claim) in self.claims.iter().enumerate() {
            add_claim_transport(index, claim, &mut bytes, &mut nodes)?;
            if bytes > MAX_PACKAGE_BYTES {
                return Err(PackageError::TransportLimit {
                    what: "serialized package size".to_string(),
                    limit: MAX_PACKAGE_BYTES,
                });
            }
            if nodes > MAX_JSON_NODES {
                return Err(PackageError::TransportLimit {
                    what: "serialized JSON nodes".to_string(),
                    limit: MAX_JSON_NODES,
                });
            }
        }
        Ok((bytes, nodes))
    }

    fn verify_finite_magnitude_sums(&self) -> Result<(), PackageError> {
        let mut verified_width = 0.0f64;
        let mut estimated_finite = 0.0f64;
        let (admissions, _) = self.admission_decisions();
        for (claim, admission) in self.claims.iter().zip(&admissions) {
            if admission.class == AdmissionClass::WaiverDependent {
                continue;
            }
            match &claim.color {
                Color::Verified { lo, hi } => {
                    if lo.is_infinite() || hi.is_infinite() {
                        // An explicit infinite endpoint is a sound, vacuous
                        // enclosure and remains visible as +inf in the actual
                        // magnitude budget. This finite-overflow guard only
                        // rejects finite inputs whose arithmetic overflows.
                        continue;
                    }
                    let width = hi - lo;
                    let next = verified_width + width;
                    if !width.is_finite() || !next.is_finite() {
                        return Err(PackageError::MagnitudeOverflow {
                            claim: claim.id.clone(),
                            component: "verified_width",
                        });
                    }
                    verified_width = next;
                }
                Color::Estimated { dispersion, .. } if dispersion.is_finite() => {
                    let next = estimated_finite + dispersion;
                    if !next.is_finite() {
                        return Err(PackageError::MagnitudeOverflow {
                            claim: claim.id.clone(),
                            component: "estimated_dispersion",
                        });
                    }
                    estimated_finite = next;
                }
                Color::Estimated { .. } | Color::Validated { .. } => {}
            }
        }
        if !(verified_width + estimated_finite).is_finite() {
            return Err(PackageError::MagnitudeOverflow {
                claim: "<aggregate>".to_string(),
                component: "quantified_total",
            });
        }
        Ok(())
    }

    /// The per-claim uncertainty MAGNITUDE attribution (bead qmao.6.1):
    /// the budget pie over error magnitudes, not claim counts. Verified
    /// claims contribute their interval width, estimated claims their
    /// dispersion; validated claims carry regional trust with no
    /// numeric bound and are reported as an unquantified COUNT rather
    /// than laundered into a number.
    pub fn magnitude_budget(&self) -> Result<MagnitudeBudget, PackageError> {
        self.verify().map(|report| report.magnitude_budget)
    }

    /// Scientific magnitude attribution after verification with explicit
    /// external capabilities.
    ///
    /// # Errors
    /// Any package, origin, falsifier, waiver, or signature refusal.
    pub fn magnitude_budget_with(
        &self,
        capabilities: &VerificationCapabilities<'_>,
    ) -> Result<MagnitudeBudget, PackageError> {
        self.verify_with(capabilities)
            .map(|report| report.magnitude_budget)
    }

    fn declared_magnitude_budget_unverified(&self) -> MagnitudeBudget {
        let (admissions, _) = self.admission_decisions();
        self.magnitude_budget_from(&admissions)
    }

    fn magnitude_budget_from(&self, admissions: &[ClaimAdmission]) -> MagnitudeBudget {
        let mut b = MagnitudeBudget::default();
        for (c, admission) in self.claims.iter().zip(admissions) {
            if admission.class == AdmissionClass::WaiverDependent {
                b.waived_unquantified += 1;
                continue;
            }
            match &c.color {
                Color::Verified { lo, hi } => {
                    if lo.is_infinite() || hi.is_infinite() {
                        b.verified_width = f64::INFINITY;
                    } else {
                        b.verified_width += hi - lo;
                    }
                }
                Color::Validated { .. } => b.validated_unquantified += 1,
                Color::Estimated { dispersion, .. } => b.estimated_dispersion += dispersion,
            }
        }
        b.quantified_total = b.verified_width + b.estimated_dispersion;
        b
    }

    /// Emit the package as deterministic, self-describing JSON —
    /// schema v7: complete color payloads, algebra-versioned receipts, and typed origins
    /// (floats as bit-exact hex),
    /// provenance, signature, the 64-hex BLAKE3 content root, and the
    /// magnitude budget. [`EvidencePackage::from_json`] round-trips this
    /// semantically and refuses anything else. Serialization is fallible so an
    /// untrusted in-memory builder cannot bypass the standalone transport
    /// envelope.
    ///
    /// # Errors
    /// [`PackageError::TransportLimit`] when the builder cannot be serialized
    /// within the package limits.
    pub fn to_json(&self) -> Result<String, PackageError> {
        self.verify_transport_limits()?;
        Ok(self.to_json_unchecked())
    }

    fn to_json_unchecked(&self) -> String {
        use core::fmt::Write as _;
        let mut out = String::new();
        let _ = write!(
            out,
            "{{\"format_version\":{},\"merkle_root\":\"{}\",\"provenance\":{{\"code_version\":\"{}\",\"constellation_lock\":\"{}\"}},\"signature\":",
            self.format_version,
            self.merkle_root_unchecked(),
            json_escape(&self.provenance.code_version),
            json_escape(&self.provenance.constellation_lock),
        );
        match &self.signature {
            Some(s) => {
                let _ = write!(out, "\"{}\"", json_escape(s));
            }
            None => out.push_str("null"),
        }
        let mb = self.declared_magnitude_budget_unverified();
        let _ = write!(
            out,
            ",\"magnitude_budget\":{{\"verified_width_bits\":\"{:016x}\",\"estimated_dispersion_bits\":\"{:016x}\",\"validated_unquantified\":{},\"waived_unquantified\":{}}}",
            mb.verified_width.to_bits(),
            mb.estimated_dispersion.to_bits(),
            mb.validated_unquantified,
            mb.waived_unquantified
        );
        out.push_str(",\"claims\":[");
        for (i, c) in self.claims.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            push_claim_json(&mut out, c);
        }
        out.push_str("]}");
        out
    }
}

fn is_blank_or_placeholder(text: &str) -> bool {
    let text = text.trim();
    text.is_empty() || is_placeholder(text)
}

fn is_placeholder(text: &str) -> bool {
    is_placeholder_token(text)
}

fn check_transport_count(what: &str, count: usize) -> Result<(), PackageError> {
    if count > MAX_JSON_CONTAINER_ITEMS {
        return Err(PackageError::TransportLimit {
            what: what.to_string(),
            limit: MAX_JSON_CONTAINER_ITEMS,
        });
    }
    Ok(())
}

fn escaped_json_len(value: &str) -> usize {
    value
        .chars()
        .map(|ch| match ch {
            '"' | '\\' | '\n' | '\r' | '\t' => 2,
            c if c.is_control() => 6,
            c => c.len_utf8(),
        })
        .sum()
}

fn add_transport_text(total: &mut usize, what: &str, value: &str) -> Result<(), PackageError> {
    if value.len() > MAX_JSON_STRING_BYTES {
        return Err(PackageError::TransportLimit {
            what: what.to_string(),
            limit: MAX_JSON_STRING_BYTES,
        });
    }
    add_transport_bytes(total, escaped_json_len(value).saturating_add(2))
}

fn add_transport_bytes(total: &mut usize, amount: usize) -> Result<(), PackageError> {
    *total = total
        .checked_add(amount)
        .ok_or_else(|| PackageError::TransportLimit {
            what: "serialized package size".to_string(),
            limit: MAX_PACKAGE_BYTES,
        })?;
    if *total > MAX_PACKAGE_BYTES {
        return Err(PackageError::TransportLimit {
            what: "serialized package size".to_string(),
            limit: MAX_PACKAGE_BYTES,
        });
    }
    Ok(())
}

fn add_transport_nodes(total: &mut usize, amount: usize) -> Result<(), PackageError> {
    *total = total
        .checked_add(amount)
        .ok_or_else(|| PackageError::TransportLimit {
            what: "serialized JSON nodes".to_string(),
            limit: MAX_JSON_NODES,
        })?;
    if *total > MAX_JSON_NODES {
        return Err(PackageError::TransportLimit {
            what: "serialized JSON nodes".to_string(),
            limit: MAX_JSON_NODES,
        });
    }
    Ok(())
}

pub(crate) fn is_canonical_content_hash(hash: &str) -> bool {
    hash.len() == 64
        && hash
            .as_bytes()
            .iter()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(byte))
}

fn push_origin_json(out: &mut String, origin: &ClaimOrigin) {
    use core::fmt::Write as _;

    match origin {
        ClaimOrigin::SourceCertificate {
            producer,
            certificate_hash,
        } => {
            let _ = write!(
                out,
                "{{\"kind\":\"source-certificate\",\"producer\":\"{}\",\"certificate_hash\":\"{}\"}}",
                json_escape(producer),
                json_escape(certificate_hash)
            );
        }
        ClaimOrigin::AnchoredSource {
            dataset_id,
            content_hash,
        } => {
            let _ = write!(
                out,
                "{{\"kind\":\"anchored-source\",\"dataset_id\":\"{}\",\"content_hash\":\"{}\"}}",
                json_escape(dataset_id),
                json_escape(content_hash)
            );
        }
        ClaimOrigin::EstimatedSource { estimator } => {
            let _ = write!(
                out,
                "{{\"kind\":\"estimated-source\",\"estimator\":\"{}\"}}",
                json_escape(estimator)
            );
        }
        ClaimOrigin::Derived => out.push_str("{\"kind\":\"derived\"}"),
        ClaimOrigin::AuthenticatedWaiver(grant) => {
            let _ = write!(
                out,
                "{{\"kind\":\"authenticated-waiver\",\"waiver_id\":\"{}\",\"expiry_day\":{},\"mac\":\"{}\"}}",
                json_escape(&grant.waiver_id),
                grant.expiry_day,
                json_escape(&grant.mac)
            );
        }
    }
}

fn push_claim_json(out: &mut String, claim: &Claim) {
    use core::fmt::Write as _;
    let _ = write!(
        out,
        "{{\"id\":\"{}\",\"statement\":\"{}\",\"color\":",
        json_escape(&claim.id),
        json_escape(&claim.statement),
    );
    match &claim.color {
        Color::Verified { lo, hi } => {
            let _ = write!(
                out,
                "{{\"kind\":\"verified\",\"lo_bits\":\"{:016x}\",\"hi_bits\":\"{:016x}\"}}",
                lo.to_bits(),
                hi.to_bits()
            );
        }
        Color::Validated { regime, dataset } => {
            let _ = write!(out, "{{\"kind\":\"validated\",\"regime\":{{");
            for (index, (axis, (lo, hi))) in regime.bounds().iter().enumerate() {
                if index > 0 {
                    out.push(',');
                }
                let _ = write!(
                    out,
                    "\"{}\":[\"{:016x}\",\"{:016x}\"]",
                    json_escape(axis),
                    lo.to_bits(),
                    hi.to_bits()
                );
            }
            let _ = write!(out, "}},\"dataset\":\"{}\"}}", json_escape(dataset));
        }
        Color::Estimated {
            estimator,
            dispersion,
        } => {
            let _ = write!(
                out,
                "{{\"kind\":\"estimated\",\"estimator\":\"{}\",\"dispersion_bits\":\"{:016x}\"}}",
                json_escape(estimator),
                dispersion.to_bits()
            );
        }
    }
    match &claim.receipt {
        Some(receipt) => {
            let _ = write!(
                out,
                ",\"receipt\":{{\"color_algebra_version\":{},\"op\":\"{}\",\"parents\":{:?},\"artifact_hash\":\"{}\"}}",
                receipt.color_algebra_version,
                op_name(receipt.op),
                receipt.parents,
                json_escape(&receipt.artifact_hash)
            );
        }
        None => out.push_str(",\"receipt\":null"),
    }
    out.push_str(",\"falsifiers\":[");
    for (index, falsifier) in claim.falsifiers.iter().enumerate() {
        if index > 0 {
            out.push(',');
        }
        let _ = write!(
            out,
            "{{\"name\":\"{}\",\"attempts\":{},\"refuted\":{},\"detail\":\"{}\",\"artifact_hash\":\"{}\"}}",
            json_escape(&falsifier.name),
            falsifier.attempts,
            falsifier.refuted,
            json_escape(&falsifier.detail),
            json_escape(&falsifier.artifact_hash)
        );
    }
    out.push_str("],\"anchors\":[");
    for (index, anchor) in claim.anchors.iter().enumerate() {
        if index > 0 {
            out.push(',');
        }
        let _ = write!(
            out,
            "{{\"dataset_id\":\"{}\",\"content_hash\":\"{}\"}}",
            json_escape(&anchor.dataset_id),
            json_escape(&anchor.content_hash)
        );
    }
    out.push_str("],\"origin\":");
    push_origin_json(out, &claim.origin);
    out.push('}');
}

/// The magnitude budget (see [`EvidencePackage::magnitude_budget`]).
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct MagnitudeBudget {
    /// Σ (hi − lo) over verified claims.
    pub verified_width: f64,
    /// Σ dispersion over estimated claims.
    pub estimated_dispersion: f64,
    /// Validated claims (regional trust, no numeric bound — counted,
    /// never converted into a fake magnitude).
    pub validated_unquantified: usize,
    /// Directly waived claims and their complete derived descendant cone.
    /// Their underlying colors contribute to no scientific magnitude subtotal.
    pub waived_unquantified: usize,
    /// verified_width + estimated_dispersion (reconciles with the
    /// parts by construction; the parser re-derives and refuses drift).
    pub quantified_total: f64,
}

fn bind_policy_fingerprint(
    slot: &mut Option<PolicyFingerprint>,
    fingerprint: PolicyFingerprint,
    capability: &'static str,
) -> Result<(), PackageError> {
    match slot {
        Some(previous) if *previous != fingerprint => Err(PackageError::PolicyFingerprintRefused {
            capability,
            why: "policy fingerprint changed during package verification",
            previous: *previous,
            observed: fingerprint,
        }),
        Some(_) => Ok(()),
        None => {
            *slot = Some(fingerprint);
            Ok(())
        }
    }
}

fn release_admission_context_hash(
    package_root: ContentHash,
    mut policies: VerificationPolicyFingerprints,
    waiver_day: Option<u64>,
    admissions: &[ClaimAdmission],
    waiver_registry: &[ReceiptWaiver],
) -> ContentHash {
    // Signature verification happens after this context is formed. Excluding
    // its own policy avoids a circular subject while binding every scientific
    // and administrative decision the release approval is meant to authorize.
    policies.signatures = None;
    let provisional_receipt = verification_receipt_hash(
        package_root,
        policies,
        waiver_day,
        &SignatureStatus::Unsigned,
        admissions,
        waiver_registry,
    );
    hash_domain(
        "fs-package:v7:release-admission-context",
        provisional_receipt.as_bytes(),
    )
}

fn verification_receipt_hash(
    package_root: ContentHash,
    policies: VerificationPolicyFingerprints,
    waiver_day: Option<u64>,
    signature: &SignatureStatus,
    admissions: &[ClaimAdmission],
    waiver_registry: &[ReceiptWaiver],
) -> ContentHash {
    fn stable_usize(value: usize) -> [u8; 8] {
        u64::try_from(value)
            .expect("schema-v7 transport bounds fit u64")
            .to_le_bytes()
    }

    fn push_atom_hash(hasher: &mut fs_blake3::Blake3, value: &[u8]) {
        hasher.update(&stable_usize(value.len()));
        hasher.update(value);
    }

    fn push_policy(hasher: &mut fs_blake3::Blake3, label: &[u8], fingerprint: Option<ContentHash>) {
        push_atom_hash(hasher, label);
        match fingerprint {
            Some(fingerprint) => push_atom_hash(hasher, fingerprint.as_bytes()),
            None => push_atom_hash(hasher, b"none"),
        }
    }

    let mut hasher = fs_blake3::Blake3::new();
    hasher.update(b"fs-package:v7:verification-receipt\0");
    push_atom_hash(&mut hasher, package_root.as_bytes());
    push_policy(
        &mut hasher,
        b"source-certificates",
        policies.source_certificates,
    );
    push_policy(&mut hasher, b"anchored-sources", policies.anchored_sources);
    push_policy(&mut hasher, b"falsifiers", policies.falsifiers);
    push_policy(&mut hasher, b"derivations", policies.derivations);
    push_policy(&mut hasher, b"waivers", policies.waivers);
    push_policy(&mut hasher, b"signatures", policies.signatures);
    match waiver_day {
        Some(day) => push_atom_hash(&mut hasher, &day.to_le_bytes()),
        None => push_atom_hash(&mut hasher, b"none"),
    }
    match signature {
        SignatureStatus::Unsigned => push_atom_hash(&mut hasher, b"signature:unsigned"),
        SignatureStatus::Refused { reason } => {
            push_atom_hash(&mut hasher, b"signature:refused");
            push_atom_hash(&mut hasher, reason.as_bytes());
        }
        SignatureStatus::Unverified(value) => {
            push_atom_hash(&mut hasher, b"signature:unverified");
            push_atom_hash(&mut hasher, value.as_bytes());
        }
        SignatureStatus::Authenticated(authenticated) => {
            push_atom_hash(&mut hasher, b"signature:authenticated");
            push_atom_hash(&mut hasher, authenticated.signature().as_bytes());
            match authenticated.purpose() {
                SignaturePurpose::PackageRootAttestation => {
                    push_atom_hash(&mut hasher, b"package-root-attestation");
                }
                SignaturePurpose::ReleaseApproval {
                    checker_protocol,
                    expected_root,
                    admission_context,
                } => {
                    push_atom_hash(&mut hasher, b"release-approval");
                    push_atom_hash(&mut hasher, &checker_protocol.to_le_bytes());
                    push_atom_hash(&mut hasher, expected_root.as_bytes());
                    push_atom_hash(&mut hasher, admission_context.as_bytes());
                }
            }
        }
    }
    push_atom_hash(&mut hasher, b"waiver-registry");
    push_atom_hash(&mut hasher, &stable_usize(waiver_registry.len()));
    for waiver in waiver_registry {
        push_atom_hash(&mut hasher, b"waiver-entry");
        push_atom_hash(&mut hasher, &stable_usize(waiver.registry_index));
        push_atom_hash(&mut hasher, &stable_usize(waiver.claim_index));
        push_atom_hash(&mut hasher, waiver.waiver_id.as_bytes());
    }
    push_atom_hash(&mut hasher, b"admissions");
    push_atom_hash(&mut hasher, &stable_usize(admissions.len()));
    for admission in admissions {
        push_atom_hash(&mut hasher, b"admission-entry");
        push_atom_hash(&mut hasher, &stable_usize(admission.claim_index));
        push_atom_hash(&mut hasher, admission.claim_id.as_bytes());
        push_atom_hash(&mut hasher, admission.origin_kind.tag().as_bytes());
        push_atom_hash(
            &mut hasher,
            match admission.class {
                AdmissionClass::Scientific => b"scientific",
                AdmissionClass::WaiverDependent => b"waiver-dependent",
            },
        );
        match admission.direct_waiver {
            Some(index) => push_atom_hash(&mut hasher, &stable_usize(index)),
            None => push_atom_hash(&mut hasher, b"none"),
        }
        push_atom_hash(&mut hasher, b"waiver-parents");
        push_atom_hash(&mut hasher, &stable_usize(admission.waiver_parents.len()));
        for parent in &admission.waiver_parents {
            push_atom_hash(&mut hasher, &stable_usize(*parent));
        }
    }
    let inner = hasher.finalize();
    hash_domain("fs-package:v7:verification-receipt", inner.as_bytes())
}

/// Combine two child hashes into a domain-separated parent node hash.
fn combine(a: &ContentHash, b: &ContentHash) -> ContentHash {
    let mut buf = [0u8; 64];
    buf[..32].copy_from_slice(a.as_bytes());
    buf[32..].copy_from_slice(b.as_bytes());
    hash_domain("fs-package:v7:node", &buf)
}

fn push_atom(out: &mut String, value: &str) {
    use core::fmt::Write as _;
    let _ = write!(out, "{}:", value.len());
    out.push_str(value);
    out.push('|');
}

/// Minimal JSON string escaping.
fn json_escape(s: &str) -> String {
    use core::fmt::Write as _;
    let mut out = String::with_capacity(s.len());
    for ch in s.chars() {
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
    out
}

// ---------------------------------------------------------------------------
// Strict schema-v7 parser: the package is a PROOF
// ARTIFACT, so parsing fails closed — unknown fields, missing fields,
// wrong types, bad hex, NaN/inverted certificates, a magnitude budget
// that does not re-derive, or an embedded root that does not recompute
// from the parsed fields are each a structured refusal.
// ---------------------------------------------------------------------------

/// A structured parse failure.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParseError {
    /// What was being parsed.
    pub what: String,
    /// Why it refused.
    pub why: String,
}

impl core::fmt::Display for ParseError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "package parse refused at {}: {}", self.what, self.why)
    }
}

impl core::error::Error for ParseError {}

/// Maximum serialized package size accepted by the standalone checker.
pub const MAX_PACKAGE_BYTES: usize = 64 * 1024 * 1024;
/// Maximum JSON nesting depth before schema mapping.
pub const MAX_JSON_DEPTH: usize = 64;
/// Maximum total JSON values in one package.
pub const MAX_JSON_NODES: usize = 1_000_000;
/// Maximum decoded bytes in a JSON string or object key.
pub const MAX_JSON_STRING_BYTES: usize = 1024 * 1024;
/// Maximum members in any one object or array.
pub const MAX_JSON_CONTAINER_ITEMS: usize = 100_000;
/// Numeric fields in this schema are bounded integers; longer tokens are
/// hostile or malformed even before exact conversion.
pub const MAX_JSON_NUMBER_BYTES: usize = 128;

/// Minimal JSON value for the strict mapper.
#[derive(Debug, Clone, PartialEq)]
enum Jv {
    Null,
    Bool(bool),
    Str(String),
    /// Raw decimal spelling. Integer-valued schema fields must never pass
    /// through `f64`, which cannot represent every `u64` exactly.
    Num(String),
    Arr(Vec<Jv>),
    Obj(Vec<(String, Jv)>),
}

impl Jv {
    fn kind(&self) -> &'static str {
        match self {
            Self::Null => "null",
            Self::Bool(_) => "boolean",
            Self::Str(_) => "string",
            Self::Num(_) => "number",
            Self::Arr(_) => "array",
            Self::Obj(_) => "object",
        }
    }
}

struct Jp<'a> {
    b: &'a [u8],
    at: usize,
    nodes: usize,
}

impl Jp<'_> {
    fn err(&self, what: &str, why: impl Into<String>) -> ParseError {
        ParseError {
            what: format!("{what} (byte {})", self.at),
            why: why.into(),
        }
    }

    fn ws(&mut self) {
        while self
            .b
            .get(self.at)
            .is_some_and(|c| matches!(c, b' ' | b'\t' | b'\n' | b'\r'))
        {
            self.at += 1;
        }
    }

    fn eat(&mut self, c: u8, what: &str) -> Result<(), ParseError> {
        self.ws();
        if self.b.get(self.at) == Some(&c) {
            self.at += 1;
            Ok(())
        } else {
            Err(self.err(what, format!("expected {:?}", char::from(c))))
        }
    }

    fn value(&mut self) -> Result<Jv, ParseError> {
        self.value_at(0)
    }

    fn object_at(&mut self, depth: usize) -> Result<Jv, ParseError> {
        self.at += 1;
        let mut fields = Vec::new();
        let mut keys = std::collections::BTreeSet::new();
        self.ws();
        if self.b.get(self.at) == Some(&b'}') {
            self.at += 1;
            return Ok(Jv::Obj(fields));
        }
        loop {
            let key = self.string()?;
            self.eat(b':', "object")?;
            let value = self.value_at(depth + 1)?;
            if !keys.insert(key.clone()) {
                return Err(self.err("object", format!("duplicate key {key:?}")));
            }
            fields.push((key, value));
            if fields.len() > MAX_JSON_CONTAINER_ITEMS {
                return Err(self.err(
                    "object",
                    format!("object member count exceeds limit {MAX_JSON_CONTAINER_ITEMS}"),
                ));
            }
            self.ws();
            match self.b.get(self.at) {
                Some(b',') => {
                    self.at += 1;
                    self.ws();
                }
                Some(b'}') => {
                    self.at += 1;
                    return Ok(Jv::Obj(fields));
                }
                _ => return Err(self.err("object", "expected ',' or '}'")),
            }
        }
    }

    fn array_at(&mut self, depth: usize) -> Result<Jv, ParseError> {
        self.at += 1;
        let mut items = Vec::new();
        self.ws();
        if self.b.get(self.at) == Some(&b']') {
            self.at += 1;
            return Ok(Jv::Arr(items));
        }
        loop {
            items.push(self.value_at(depth + 1)?);
            if items.len() > MAX_JSON_CONTAINER_ITEMS {
                return Err(self.err(
                    "array",
                    format!("array element count exceeds limit {MAX_JSON_CONTAINER_ITEMS}"),
                ));
            }
            self.ws();
            match self.b.get(self.at) {
                Some(b',') => self.at += 1,
                Some(b']') => {
                    self.at += 1;
                    return Ok(Jv::Arr(items));
                }
                _ => return Err(self.err("array", "expected ',' or ']'")),
            }
        }
    }

    fn value_at(&mut self, depth: usize) -> Result<Jv, ParseError> {
        if depth > MAX_JSON_DEPTH {
            return Err(self.err(
                "value",
                format!("nesting depth exceeds limit {MAX_JSON_DEPTH}"),
            ));
        }
        self.nodes = self.nodes.checked_add(1).ok_or_else(|| {
            self.err(
                "value",
                "JSON node counter overflowed before schema mapping",
            )
        })?;
        if self.nodes > MAX_JSON_NODES {
            return Err(self.err(
                "value",
                format!("JSON node count exceeds limit {MAX_JSON_NODES}"),
            ));
        }
        self.ws();
        match self.b.get(self.at) {
            Some(b'"') => Ok(Jv::Str(self.string()?)),
            Some(b'{') => self.object_at(depth),
            Some(b'[') => self.array_at(depth),
            Some(b'n') => {
                if self.b[self.at..].starts_with(b"null") {
                    self.at += 4;
                    Ok(Jv::Null)
                } else {
                    Err(self.err("literal", "unknown literal"))
                }
            }
            Some(b't') => {
                if self.b[self.at..].starts_with(b"true") {
                    self.at += 4;
                    Ok(Jv::Bool(true))
                } else {
                    Err(self.err("literal", "unknown literal"))
                }
            }
            Some(b'f') => {
                if self.b[self.at..].starts_with(b"false") {
                    self.at += 5;
                    Ok(Jv::Bool(false))
                } else {
                    Err(self.err("literal", "unknown literal"))
                }
            }
            Some(c) if c.is_ascii_digit() || *c == b'-' => {
                let start = self.at;
                while self.b.get(self.at).is_some_and(|c| {
                    c.is_ascii_digit() || matches!(c, b'-' | b'+' | b'.' | b'e' | b'E')
                }) {
                    self.at += 1;
                }
                if self.at - start > MAX_JSON_NUMBER_BYTES {
                    return Err(self.err(
                        "number",
                        format!("number token exceeds {MAX_JSON_NUMBER_BYTES} bytes"),
                    ));
                }
                let text = core::str::from_utf8(&self.b[start..self.at]).unwrap_or("");
                text.parse::<f64>()
                    .map(|_| Jv::Num(text.to_string()))
                    .map_err(|_| self.err("number", format!("bad number {text:?}")))
            }
            _ => Err(self.err("value", "unexpected byte or end of input")),
        }
    }

    fn string(&mut self) -> Result<String, ParseError> {
        self.ws();
        if self.b.get(self.at) != Some(&b'"') {
            return Err(self.err("string", "expected '\"'"));
        }
        self.at += 1;
        let mut out = String::new();
        loop {
            match self.b.get(self.at) {
                None => return Err(self.err("string", "unterminated")),
                Some(b'"') => {
                    self.at += 1;
                    return Ok(out);
                }
                Some(b'\\') => {
                    self.at += 1;
                    match self.b.get(self.at) {
                        Some(b'"') => out.push('"'),
                        Some(b'\\') => out.push('\\'),
                        Some(b'n') => out.push('\n'),
                        Some(b'r') => out.push('\r'),
                        Some(b't') => out.push('\t'),
                        Some(b'u') => {
                            let hex = self
                                .b
                                .get(self.at + 1..self.at + 5)
                                .and_then(|h| core::str::from_utf8(h).ok())
                                .and_then(|h| u32::from_str_radix(h, 16).ok())
                                .and_then(char::from_u32)
                                .ok_or_else(|| self.err("string", "bad \\u escape"))?;
                            out.push(hex);
                            self.at += 4;
                        }
                        _ => return Err(self.err("string", "bad escape")),
                    }
                    if out.len() > MAX_JSON_STRING_BYTES {
                        return Err(self.err(
                            "string",
                            format!("decoded string exceeds {MAX_JSON_STRING_BYTES} bytes"),
                        ));
                    }
                    self.at += 1;
                }
                Some(&c) if c < 0x20 => {
                    return Err(self.err("string", "unescaped control character"));
                }
                Some(&c) => {
                    // Multi-byte UTF-8 passes through byte-wise.
                    let len = if c < 0x80 {
                        1
                    } else if c >> 5 == 0b110 {
                        2
                    } else if c >> 4 == 0b1110 {
                        3
                    } else {
                        4
                    };
                    let chunk = self
                        .b
                        .get(self.at..self.at + len)
                        .and_then(|ch| core::str::from_utf8(ch).ok())
                        .ok_or_else(|| self.err("string", "invalid UTF-8"))?;
                    out.push_str(chunk);
                    if out.len() > MAX_JSON_STRING_BYTES {
                        return Err(self.err(
                            "string",
                            format!("decoded string exceeds {MAX_JSON_STRING_BYTES} bytes"),
                        ));
                    }
                    self.at += len;
                }
            }
        }
    }
}

fn obj_fields(v: Jv, what: &str) -> Result<Vec<(String, Jv)>, ParseError> {
    match v {
        Jv::Obj(f) => Ok(f),
        other => Err(ParseError {
            what: what.to_string(),
            why: format!("expected an object, got {}", other.kind()),
        }),
    }
}

/// Take field `key` from `fields`; strict mappers call this for every
/// expected key and then refuse leftovers.
fn take_field(fields: &mut Vec<(String, Jv)>, key: &str, what: &str) -> Result<Jv, ParseError> {
    let idx = fields
        .iter()
        .position(|(k, _)| k == key)
        .ok_or(ParseError {
            what: what.to_string(),
            why: format!("missing required field {key:?}"),
        })?;
    Ok(fields.remove(idx).1)
}

fn no_leftovers(fields: &[(String, Jv)], what: &str) -> Result<(), ParseError> {
    if let Some((k, _)) = fields.first() {
        return Err(ParseError {
            what: what.to_string(),
            why: format!("unknown field {k:?} (schema v7 is closed — fail closed)"),
        });
    }
    Ok(())
}

fn as_str(v: Jv, what: &str) -> Result<String, ParseError> {
    match v {
        Jv::Str(s) => Ok(s),
        other => Err(ParseError {
            what: what.to_string(),
            why: format!("expected a string, got {}", other.kind()),
        }),
    }
}

fn hex_u64(v: Jv, what: &str) -> Result<u64, ParseError> {
    let hex = as_str(v, what)?;
    if hex.len() != 16
        || !hex
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
    {
        return Err(ParseError {
            what: what.to_string(),
            why: format!("expected 16 hex digits, got {hex:?}"),
        });
    }
    Ok(u64::from_str_radix(&hex, 16).expect("validated hexadecimal u64"))
}

fn bits_f64(v: Jv, what: &str, must_be_finite: bool) -> Result<f64, ParseError> {
    let value = f64::from_bits(hex_u64(v, what)?);
    if must_be_finite && !value.is_finite() {
        return Err(ParseError {
            what: what.to_string(),
            why: format!("non-finite value {value} where a finite certificate is required"),
        });
    }
    Ok(value)
}

fn decimal_u64(v: Jv, what: &str) -> Result<u64, ParseError> {
    match v {
        Jv::Num(text)
            if !text.is_empty()
                && text.bytes().all(|b| b.is_ascii_digit())
                && (text == "0" || !text.starts_with('0')) =>
        {
            text.parse::<u64>().map_err(|_| ParseError {
                what: what.to_string(),
                why: format!("unsigned integer out of range: {text:?}"),
            })
        }
        other => Err(ParseError {
            what: what.to_string(),
            why: format!("expected an unsigned decimal integer, got {}", other.kind()),
        }),
    }
}

fn decimal_usize(v: Jv, what: &str) -> Result<usize, ParseError> {
    let value = decimal_u64(v, what)?;
    usize::try_from(value).map_err(|_| ParseError {
        what: what.to_string(),
        why: format!("unsigned integer {value} does not fit usize"),
    })
}

fn parse_package_fields(text: &str) -> Result<Vec<(String, Jv)>, ParseError> {
    if text.len() > MAX_PACKAGE_BYTES {
        return Err(ParseError {
            what: "package".to_string(),
            why: format!("input exceeds the {MAX_PACKAGE_BYTES}-byte package limit"),
        });
    }
    let mut parser = Jp {
        b: text.as_bytes(),
        at: 0,
        nodes: 0,
    };
    let root = parser.value()?;
    parser.ws();
    if parser.at != parser.b.len() {
        return Err(ParseError {
            what: "package".to_string(),
            why: "trailing bytes after the package object".to_string(),
        });
    }
    obj_fields(root, "package")
}

fn parse_format_version(fields: &mut Vec<(String, Jv)>) -> Result<u32, ParseError> {
    let raw = decimal_u64(
        take_field(fields, "format_version", "package")?,
        "format_version",
    )?;
    let version = u32::try_from(raw).map_err(|_| ParseError {
        what: "format_version".to_string(),
        why: format!("version {raw} does not fit u32"),
    })?;
    if version != FORMAT_VERSION {
        return Err(ParseError {
            what: "format_version".to_string(),
            why: format!("unsupported version {version} (this build reads {FORMAT_VERSION})"),
        });
    }
    Ok(version)
}

fn parse_provenance(fields: &mut Vec<(String, Jv)>) -> Result<Provenance, ParseError> {
    let mut provenance = obj_fields(take_field(fields, "provenance", "package")?, "provenance")?;
    let parsed = Provenance {
        code_version: as_str(
            take_field(&mut provenance, "code_version", "provenance")?,
            "code_version",
        )?,
        constellation_lock: as_str(
            take_field(&mut provenance, "constellation_lock", "provenance")?,
            "constellation_lock",
        )?,
    };
    no_leftovers(&provenance, "provenance")?;
    Ok(parsed)
}

fn parse_signature(fields: &mut Vec<(String, Jv)>) -> Result<Option<String>, ParseError> {
    match take_field(fields, "signature", "package")? {
        Jv::Null => Ok(None),
        Jv::Str(signature) => Ok(Some(signature)),
        other => Err(ParseError {
            what: "signature".to_string(),
            why: format!("expected a string or null, got {}", other.kind()),
        }),
    }
}

fn parse_magnitude_budget(fields: &mut Vec<(String, Jv)>) -> Result<MagnitudeBudget, ParseError> {
    let mut budget = obj_fields(
        take_field(fields, "magnitude_budget", "package")?,
        "magnitude_budget",
    )?;
    let verified_width = bits_f64(
        take_field(&mut budget, "verified_width_bits", "magnitude_budget")?,
        "verified_width_bits",
        false,
    )?;
    let estimated_dispersion = bits_f64(
        take_field(&mut budget, "estimated_dispersion_bits", "magnitude_budget")?,
        "estimated_dispersion_bits",
        false,
    )?;
    let validated_unquantified = decimal_usize(
        take_field(&mut budget, "validated_unquantified", "magnitude_budget")?,
        "validated_unquantified",
    )?;
    let waived_unquantified = decimal_usize(
        take_field(&mut budget, "waived_unquantified", "magnitude_budget")?,
        "waived_unquantified",
    )?;
    no_leftovers(&budget, "magnitude_budget")?;
    Ok(MagnitudeBudget {
        verified_width,
        estimated_dispersion,
        validated_unquantified,
        waived_unquantified,
        quantified_total: verified_width + estimated_dispersion,
    })
}

fn parse_claims(fields: &mut Vec<(String, Jv)>) -> Result<Vec<Claim>, ParseError> {
    let values = match take_field(fields, "claims", "package")? {
        Jv::Arr(items) => items,
        other => {
            return Err(ParseError {
                what: "claims".to_string(),
                why: format!("expected an array, got {}", other.kind()),
            });
        }
    };
    values
        .into_iter()
        .enumerate()
        .map(|(index, value)| parse_claim(value, index))
        .collect()
}

fn verify_declarations(
    package: &EvidencePackage,
    declared_root: ContentHash,
    declared_budget: MagnitudeBudget,
) -> Result<(), ParseError> {
    let recomputed_budget = package.declared_magnitude_budget_unverified();
    if recomputed_budget.verified_width.to_bits() != declared_budget.verified_width.to_bits()
        || recomputed_budget.estimated_dispersion.to_bits()
            != declared_budget.estimated_dispersion.to_bits()
        || recomputed_budget.validated_unquantified != declared_budget.validated_unquantified
        || recomputed_budget.waived_unquantified != declared_budget.waived_unquantified
    {
        return Err(ParseError {
            what: "magnitude_budget".to_string(),
            why: "declared budget does not re-derive from the claims (tamper or drift)".to_string(),
        });
    }
    let recomputed_root = package.try_merkle_root().map_err(|error| ParseError {
        what: "merkle_root".to_string(),
        why: format!("parsed package exceeds transport limits: {error:?}"),
    })?;
    if recomputed_root != declared_root {
        return Err(ParseError {
            what: "merkle_root".to_string(),
            why: format!(
                "embedded root {declared_root} does not recompute from the parsed fields \
                 (got {recomputed_root}) — tampered or forged content"
            ),
        });
    }
    Ok(())
}

/// Parse the embedded content root: exactly 64 hex chars (schema v7).
/// A 16-hex value is the legacy v3 FNV root and is named in the refusal.
fn parse_declared_root(fields: &mut Vec<(String, Jv)>) -> Result<ContentHash, ParseError> {
    let raw = match take_field(fields, "merkle_root", "package")? {
        Jv::Str(s) => s,
        other => {
            return Err(ParseError {
                what: "merkle_root".to_string(),
                why: format!("expected a hex string, got {}", other.kind()),
            });
        }
    };
    let canonical = raw.len() == 64
        && raw
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte));
    if !canonical {
        return Err(ParseError {
            what: "merkle_root".to_string(),
            why: match raw.len() {
                16 => "a 16-hex root is the legacy v3 FNV content address; schema v7 requires \
                       the 64-hex lowercase BLAKE3 root"
                    .to_string(),
                64 => "the 64-char root must use lowercase hexadecimal".to_string(),
                n => format!("expected exactly 64 lowercase hex chars, got {n} chars"),
            },
        });
    }
    ContentHash::from_hex(&raw).ok_or_else(|| ParseError {
        what: "merkle_root".to_string(),
        why: match raw.len() {
            16 => "a 16-hex root is the legacy v3 FNV content address; schema v7 requires \
                   the 64-hex lowercase BLAKE3 root"
                .to_string(),
            64 => "the 64-char root contains non-lowercase-hex characters".to_string(),
            n => format!("expected exactly 64 lowercase hex chars, got {n} chars"),
        },
    })
}

impl EvidencePackage {
    /// Parse the deterministic schema-v7 FrankenSim JSON profile structurally:
    /// every field
    /// mapped, unknown fields refused, floats reconstructed bit-exactly,
    /// the magnitude budget re-derived and compared, and the embedded
    /// content root recomputed from the parsed fields — a package whose
    /// root does not recompute is tampered or forged, and never loads.
    /// Capability-gated sources, anchors, falsifiers, derivations, waivers, and
    /// signatures are retained but not authenticated; call
    /// [`EvidencePackage::from_json_with`] or
    /// [`EvidencePackage::verify_with`] before using them as evidence.
    /// Earlier schema versions (v3's 16-hex FNV root) are refused by
    /// version before any field is interpreted.
    ///
    /// # Errors
    /// [`ParseError`] naming the field and the refusal.
    pub fn from_json(text: &str) -> Result<EvidencePackage, ParseError> {
        let mut fields = parse_package_fields(text)?;
        let format_version = parse_format_version(&mut fields)?;
        let declared_root = parse_declared_root(&mut fields)?;
        let provenance = parse_provenance(&mut fields)?;
        let signature = parse_signature(&mut fields)?;
        let declared_budget = parse_magnitude_budget(&mut fields)?;
        let claims = parse_claims(&mut fields)?;
        no_leftovers(&fields, "package")?;
        let pkg = EvidencePackage {
            format_version,
            claims,
            provenance,
            signature,
        };
        verify_declarations(&pkg, declared_root, declared_budget)?;
        pkg.verify_structural().map_err(|error| ParseError {
            what: "package semantics".to_string(),
            why: error.to_string(),
        })?;
        Ok(pkg)
    }

    /// Parse a package and authenticate every capability-gated origin before
    /// returning it.
    ///
    /// # Errors
    /// [`ParseError`] for syntax, transport, integrity, semantic, artifact,
    /// signature, or waiver refusal.
    pub fn from_json_with(
        text: &str,
        capabilities: &VerificationCapabilities<'_>,
    ) -> Result<VerifiedPackage, ParseError> {
        let package = Self::from_json(text)?;
        let report = package
            .verify_with(capabilities)
            .map_err(|error| ParseError {
                what: "package verification capabilities".to_string(),
                why: error.to_string(),
            })?;
        Ok(VerifiedPackage { package, report })
    }
}

fn parse_claim(v: Jv, index: usize) -> Result<Claim, ParseError> {
    let what = format!("claims[{index}]");
    let mut f = obj_fields(v, &what)?;
    let id = as_str(take_field(&mut f, "id", &what)?, &what)?;
    let statement = as_str(take_field(&mut f, "statement", &what)?, &what)?;
    let color = parse_color(take_field(&mut f, "color", &what)?, &what)?;
    let receipt_v = take_field(&mut f, "receipt", &what)?;
    let falsifiers_v = take_field(&mut f, "falsifiers", &what)?;
    let anchors_v = take_field(&mut f, "anchors", &what)?;
    let origin_v = take_field(&mut f, "origin", &what)?;
    no_leftovers(&f, &what)?;
    Ok(Claim {
        id,
        statement,
        color,
        receipt: parse_receipt(receipt_v, &what)?,
        falsifiers: parse_falsifiers(falsifiers_v, &what)?,
        anchors: parse_anchors(anchors_v, &what)?,
        origin: parse_origin(origin_v, &what)?,
    })
}

fn parse_origin(value: Jv, what: &str) -> Result<ClaimOrigin, ParseError> {
    let mut fields = obj_fields(value, what)?;
    let kind = as_str(take_field(&mut fields, "kind", what)?, what)?;
    let origin = match kind.as_str() {
        "source-certificate" => ClaimOrigin::SourceCertificate {
            producer: as_str(take_field(&mut fields, "producer", what)?, what)?,
            certificate_hash: as_str(take_field(&mut fields, "certificate_hash", what)?, what)?,
        },
        "anchored-source" => ClaimOrigin::AnchoredSource {
            dataset_id: as_str(take_field(&mut fields, "dataset_id", what)?, what)?,
            content_hash: as_str(take_field(&mut fields, "content_hash", what)?, what)?,
        },
        "estimated-source" => ClaimOrigin::EstimatedSource {
            estimator: as_str(take_field(&mut fields, "estimator", what)?, what)?,
        },
        "derived" => ClaimOrigin::Derived,
        "authenticated-waiver" => ClaimOrigin::AuthenticatedWaiver(WaiverGrant {
            waiver_id: as_str(take_field(&mut fields, "waiver_id", what)?, what)?,
            expiry_day: decimal_u64(take_field(&mut fields, "expiry_day", what)?, "expiry_day")?,
            mac: as_str(take_field(&mut fields, "mac", what)?, what)?,
        }),
        other => {
            return Err(ParseError {
                what: what.to_string(),
                why: format!("unknown origin kind {other:?} — fail closed"),
            });
        }
    };
    no_leftovers(&fields, "claim origin")?;
    Ok(origin)
}

fn parse_color(value: Jv, what: &str) -> Result<Color, ParseError> {
    let mut fields = obj_fields(value, what)?;
    let kind = as_str(take_field(&mut fields, "kind", what)?, what)?;
    let color = match kind.as_str() {
        "verified" => {
            let lo = bits_f64(take_field(&mut fields, "lo_bits", what)?, what, false)?;
            let hi = bits_f64(take_field(&mut fields, "hi_bits", what)?, what, false)?;
            if lo > hi {
                return Err(ParseError {
                    what: what.to_string(),
                    why: format!("verified interval inverted: {lo} > {hi}"),
                });
            }
            Color::Verified { lo, hi }
        }
        "validated" => {
            let regime_fields = obj_fields(take_field(&mut fields, "regime", what)?, what)?;
            let mut domain = fs_evidence::ValidityDomain::unconstrained();
            for (param, bounds) in regime_fields {
                let Jv::Arr(pair) = bounds else {
                    return Err(ParseError {
                        what: what.to_string(),
                        why: format!("regime {param:?} must be a [lo_bits, hi_bits] pair"),
                    });
                };
                let [lo_v, hi_v]: [Jv; 2] = pair.try_into().map_err(|_| ParseError {
                    what: what.to_string(),
                    why: format!("regime {param:?} must have exactly two bounds"),
                })?;
                let lo = bits_f64(lo_v, what, true)?;
                let hi = bits_f64(hi_v, what, true)?;
                if param.trim().is_empty() {
                    return Err(ParseError {
                        what: what.to_string(),
                        why: "regime axis name must be non-blank".to_string(),
                    });
                }
                if lo > hi {
                    return Err(ParseError {
                        what: what.to_string(),
                        why: format!("regime axis {param:?} has inverted bounds: {lo} > {hi}"),
                    });
                }
                domain = domain.with(param, lo, hi);
            }
            let dataset = as_str(take_field(&mut fields, "dataset", what)?, what)?;
            Color::Validated {
                regime: domain,
                dataset,
            }
        }
        "estimated" => {
            let estimator = as_str(take_field(&mut fields, "estimator", what)?, what)?;
            let dispersion = bits_f64(
                take_field(&mut fields, "dispersion_bits", what)?,
                what,
                false,
            )?;
            if dispersion.is_nan() || dispersion < 0.0 {
                return Err(ParseError {
                    what: what.to_string(),
                    why: format!("NaN or negative dispersion {dispersion}"),
                });
            }
            Color::Estimated {
                estimator,
                dispersion,
            }
        }
        other => {
            return Err(ParseError {
                what: what.to_string(),
                why: format!("unknown color kind {other:?} — fail closed"),
            });
        }
    };
    no_leftovers(&fields, "claim color")?;
    validate_color_payload(&color).map_err(|error| ParseError {
        what: what.to_string(),
        why: error.to_string(),
    })?;
    Ok(color)
}

fn parse_receipt(value: Jv, what: &str) -> Result<Option<CompositionReceipt>, ParseError> {
    let Jv::Obj(mut fields) = value else {
        return match value {
            Jv::Null => Ok(None),
            other => Err(ParseError {
                what: what.to_string(),
                why: format!("receipt must be an object or null, got {}", other.kind()),
            }),
        };
    };
    let color_algebra_version = u32::try_from(decimal_u64(
        take_field(&mut fields, "color_algebra_version", what)?,
        "color_algebra_version",
    )?)
    .map_err(|_| ParseError {
        what: what.to_string(),
        why: "color_algebra_version exceeds u32".to_string(),
    })?;
    if color_algebra_version != fs_evidence::COLOR_ALGEBRA_VERSION {
        return Err(ParseError {
            what: what.to_string(),
            why: format!(
                "unsupported color algebra {color_algebra_version} (this build reads {})",
                fs_evidence::COLOR_ALGEBRA_VERSION
            ),
        });
    }
    let op_name = as_str(take_field(&mut fields, "op", what)?, what)?;
    let op = op_parse(&op_name).ok_or_else(|| ParseError {
        what: what.to_string(),
        why: format!("unknown receipt op {op_name:?} — fail closed"),
    })?;
    let parents = match take_field(&mut fields, "parents", what)? {
        Jv::Arr(items) => items
            .into_iter()
            .map(|value| decimal_usize(value, what))
            .collect::<Result<Vec<usize>, ParseError>>()?,
        other => {
            return Err(ParseError {
                what: what.to_string(),
                why: format!("receipt parents must be an array, got {}", other.kind()),
            });
        }
    };
    let artifact_hash = as_str(take_field(&mut fields, "artifact_hash", what)?, what)?;
    no_leftovers(&fields, "claim receipt")?;
    Ok(Some(CompositionReceipt {
        color_algebra_version,
        parents,
        op,
        artifact_hash,
    }))
}

fn parse_falsifiers(value: Jv, what: &str) -> Result<Vec<FalsifierRecord>, ParseError> {
    let Jv::Arr(items) = value else {
        return Err(ParseError {
            what: what.to_string(),
            why: "falsifiers must be an array".to_string(),
        });
    };
    items
        .into_iter()
        .map(|value| {
            let mut fields = obj_fields(value, what)?;
            let name = as_str(take_field(&mut fields, "name", what)?, what)?;
            let attempts = decimal_u64(take_field(&mut fields, "attempts", what)?, what)?;
            let refuted = match take_field(&mut fields, "refuted", what)? {
                Jv::Bool(value) => value,
                other => {
                    return Err(ParseError {
                        what: what.to_string(),
                        why: format!("falsifier refuted must be a bool, got {}", other.kind()),
                    });
                }
            };
            let detail = as_str(take_field(&mut fields, "detail", what)?, what)?;
            let artifact_hash = as_str(take_field(&mut fields, "artifact_hash", what)?, what)?;
            no_leftovers(&fields, "falsifier record")?;
            Ok(FalsifierRecord {
                name,
                attempts,
                refuted,
                detail,
                artifact_hash,
            })
        })
        .collect()
}

fn parse_anchors(value: Jv, what: &str) -> Result<Vec<AnchorRecord>, ParseError> {
    let Jv::Arr(items) = value else {
        return Err(ParseError {
            what: what.to_string(),
            why: "anchors must be an array".to_string(),
        });
    };
    items
        .into_iter()
        .map(|value| {
            let mut fields = obj_fields(value, what)?;
            let record = AnchorRecord {
                dataset_id: as_str(take_field(&mut fields, "dataset_id", what)?, what)?,
                content_hash: as_str(take_field(&mut fields, "content_hash", what)?, what)?,
            };
            no_leftovers(&fields, "anchor record")?;
            Ok(record)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn magnitude_budget_never_subtracts_same_signed_infinite_endpoints() {
        for value in [f64::NEG_INFINITY, f64::INFINITY] {
            let claim = Claim {
                id: "vacuous".to_string(),
                statement: "sound but vacuous enclosure".to_string(),
                color: Color::Verified {
                    lo: value,
                    hi: value,
                },
                receipt: None,
                falsifiers: Vec::new(),
                anchors: Vec::new(),
                origin: ClaimOrigin::Derived,
            };
            let package = EvidencePackage::new(Provenance::new("commit", "lock")).with_claim(claim);
            let admission = ClaimAdmission {
                claim_index: 0,
                claim_id: "vacuous".to_string(),
                origin_kind: AdmissionOriginKind::Derived,
                class: AdmissionClass::Scientific,
                direct_waiver: None,
                waiver_parents: Vec::new(),
            };
            let budget = package.magnitude_budget_from(&[admission]);
            assert!(budget.verified_width.is_infinite());
            assert!(!budget.verified_width.is_nan());
            assert!(budget.quantified_total.is_infinite());
        }
    }

    #[test]
    fn anchored_origin_requires_the_exact_attached_content_hash() {
        let origin_hash = "11".repeat(32);
        let unrelated_hash = "22".repeat(32);
        let claim = Claim {
            id: "validated".to_string(),
            statement: "matches reference data".to_string(),
            color: Color::Validated {
                regime: fs_evidence::ValidityDomain::unconstrained().with("Re", 1.0, 2.0),
                dataset: "wind-tunnel".to_string(),
            },
            receipt: None,
            falsifiers: Vec::new(),
            anchors: vec![AnchorRecord {
                dataset_id: "wind-tunnel".to_string(),
                content_hash: unrelated_hash,
            }],
            origin: ClaimOrigin::AnchoredSource {
                dataset_id: "wind-tunnel".to_string(),
                content_hash: origin_hash,
            },
        };
        assert!(!claim.has_declared_matching_validated_anchor_unverified());
        let package = EvidencePackage::new(Provenance::new("commit", "lock")).with_claim(claim);
        assert!(matches!(
            package.verify(),
            Err(PackageError::OriginMismatch {
                origin: "anchored-source",
                ..
            })
        ));
    }

    #[test]
    fn transport_node_accounting_matches_emitted_json_shapes() {
        let empty = EvidencePackage::new(Provenance::new("commit", "lock"));
        assert_eq!(
            empty.transport_usage().expect("bounded empty package").1,
            13
        );

        let estimated = EvidencePackage::new(Provenance::new("commit", "lock"))
            .with_claim(Claim::estimated("estimate", "bounded", "estimator", 1.0));
        assert_eq!(
            estimated
                .transport_usage()
                .expect("bounded estimated package")
                .1,
            26
        );

        let verified = EvidencePackage::new(Provenance::new("commit", "lock")).with_claim(
            Claim::from_certificate(
                "certificate",
                "bounded",
                0.0,
                1.0,
                "producer",
                "11".repeat(32),
            ),
        );
        assert_eq!(
            verified
                .transport_usage()
                .expect("bounded certificate package")
                .1,
            27
        );
    }
}
