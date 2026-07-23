//! CLAIM ORIGINS (schema v6): where a claim's certificate
//! CAME FROM — the missing half of "machine-checkable".
//!
//! Schema v4 made the content address collision-resistant, but content
//! consistency is not evidence origin: `Color::Verified { lo, hi }` is
//! public algebra (fs-evidence composition needs it), so any producer
//! could mint a finite interval and the standalone checker would pass
//! it. v5 sealed claims AT THE PACKAGE BOUNDARY; v6 additionally makes
//! every external trust decision an explicit, fingerprinted capability.
//! Every claim must carry
//! a [`ClaimOrigin`] consistent with its color, bound into the content
//! address, and re-derivable by the checker — while the Color algebra
//! itself stays public and untouched.
//!
//! The five origins and their re-derivation obligations:
//! - [`ClaimOrigin::SourceCertificate`] — a named producer plus the
//!   64-hex content hash of its certificate artifact (solver
//!   certificate, proof object). Shape is checked locally; a positive
//!   verdict requires an injected [`SourceCertificateVerifier`] to
//!   establish the exact typed claim request. The artifact hash makes
//!   the certificate subpoenable without shipping it.
//! - [`ClaimOrigin::AnchoredSource`] — a validated claim's reference
//!   dataset by id + content hash; must MATCH the color's named
//!   dataset exactly and be accepted by an injected
//!   [`AnchoredSourceVerifier`] over the exact claim and regime.
//! - [`ClaimOrigin::EstimatedSource`] — the estimator identity; must
//!   match the color's estimator string exactly.
//! - [`ClaimOrigin::Derived`] — a composition receipt; the checker
//!   re-runs `compose` over the parents and the result must equal the
//!   claimed color bit-exactly (the v3 receipt machinery, now the
//!   origin itself). A derived `Validated` claim must also carry a matching
//!   dataset anchor, and every matching anchor is independently admitted by
//!   [`AnchoredSourceVerifier`]; the derivation verifier cannot authorize it.
//! - [`ClaimOrigin::AuthenticatedWaiver`] — an explicit, expiring,
//!   MAC'd grant. NEVER self-authorizing: verification requires an
//!   INJECTED [`WaiverVerifier`] capability plus a date context; the
//!   in-tree default refuses everything (the fs-ledger fail-closed
//!   pattern). The MAC binds the claim's canonical bytes, so a waiver
//!   replayed onto a different claim fails.

use core::fmt;
use fs_evidence::color_leaf_identity_reason;

use crate::{ContentHash, Provenance, SemanticWitness};

/// Semantic version of the typed detached-signature subject identity.
pub const SIGNATURE_SUBJECT_IDENTITY_VERSION: u32 = 9;
/// Exact BLAKE3 domain for the typed detached-signature subject.
pub const SIGNATURE_SUBJECT_IDENTITY_DOMAIN: &str = "fs-package:v9:signature-subject";
const _: () = assert!(SIGNATURE_SUBJECT_IDENTITY_VERSION == crate::FORMAT_VERSION);

/// Owner-local declaration consumed by `xtask check-identities`.
#[allow(dead_code)]
pub const SIGNATURE_SUBJECT_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-package:signature-subject",
    "version_const=SIGNATURE_SUBJECT_IDENTITY_VERSION",
    "version=9",
    "domain=fs-package:v9:signature-subject",
    "domain_const=SIGNATURE_SUBJECT_IDENTITY_DOMAIN",
    "encoder=signature_subject_hash",
    "encoder_helpers=signature_subject_hash_with_domain,SignatureSubjectFields::from_inputs",
    "schema_constants=SIGNATURE_SUBJECT_IDENTITY_VERSION,SIGNATURE_SUBJECT_IDENTITY_DOMAIN,crates/fs-package/src/lib.rs#FORMAT_VERSION",
    "schema_functions=crates/fs-package/src/lib.rs#admit_retained_content_hash",
    "schema_dependencies=fs-package:package-root,fs-package:release-admission-context",
    "digest=blake3-derive-key",
    "encoding=typed-binary",
    "sources=SignatureSubjectFields",
    "source_fields=SignatureSubjectFields.package_root:semantic,SignatureSubjectFields.purpose_tag:semantic,SignatureSubjectFields.checker_protocol:semantic,SignatureSubjectFields.expected_root:semantic,SignatureSubjectFields.admission_context:semantic,SignatureSubjectFields.semantic_context:semantic",
    "source_bindings=SignatureSubjectFields.package_root>package-root,SignatureSubjectFields.purpose_tag>purpose-tag,SignatureSubjectFields.checker_protocol>checker-protocol-presence-and-value,SignatureSubjectFields.expected_root>expected-root-presence-and-value,SignatureSubjectFields.admission_context>admission-context-presence-and-value,SignatureSubjectFields.semantic_context>semantic-context-presence-and-value",
    "external_semantic_fields=identity-version,digest-domain",
    "semantic_fields=identity-version,digest-domain,package-root,purpose-tag,checker-protocol-presence-and-value,expected-root-presence-and-value,admission-context-presence-and-value,semantic-context-presence-and-value",
    "excluded_fields=detached-signature-bytes:signature-output-not-subject-input",
    "consumers=signature_subject_hash,SignatureRequest::subject_hash,SignatureVerifier::verify,release-signature-producers",
    "mutations=identity-version:crates/fs-package/tests/package.rs#package_identity_versions_and_transports_fail_closed,digest-domain:crates/fs-package/src/origin.rs#signature_subject_identity_fields_move_independently,package-root:crates/fs-package/src/origin.rs#signature_subject_identity_fields_move_independently,purpose-tag:crates/fs-package/src/origin.rs#signature_subject_identity_fields_move_independently,checker-protocol-presence-and-value:crates/fs-package/src/origin.rs#signature_subject_identity_fields_move_independently,expected-root-presence-and-value:crates/fs-package/src/origin.rs#signature_subject_identity_fields_move_independently,admission-context-presence-and-value:crates/fs-package/src/origin.rs#signature_subject_identity_fields_move_independently,semantic-context-presence-and-value:crates/fs-package/src/origin.rs#signature_subject_identity_fields_move_independently",
    "nonsemantic_mutations=detached-signature-bytes:crates/fs-package/src/origin.rs#signature_subject_excludes_detached_signature_bytes",
    "field_guard=classify_signature_subject_identity_fields",
    "transport_guard=admit_retained_signature_subject_hash",
    "version_guard=crates/fs-package/tests/package.rs#package_identity_versions_and_transports_fail_closed",
    "coupling_surface=fs-package:signature-subject",
];

/// A stable identity for one external verification policy.
///
/// Implementations must change this fingerprint whenever their trust roots,
/// accepted artifact set, validation algorithm, or other decision semantics
/// change. Verification receipts bind these values so replay cannot silently
/// substitute a different policy.
pub type PolicyFingerprint = ContentHash;

/// Atomic result of one external policy decision. Acceptance and the exact
/// policy identity are returned by the same callback, preventing a mutable
/// verifier from making a decision under one policy while separately reporting
/// another fingerprint.
///
/// ```compile_fail
/// use fs_package::{ContentHash, VerificationDecision};
///
/// // Decisions must be constructed atomically through `accept` or `reject`.
/// let forged = VerificationDecision {
///     accepted: true,
///     policy_fingerprint: ContentHash([0; 32]),
/// };
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct VerificationDecision {
    /// Whether this exact request was accepted.
    accepted: bool,
    /// Stable identity of the policy that made this decision.
    policy_fingerprint: PolicyFingerprint,
}

impl VerificationDecision {
    /// Accepted under `policy_fingerprint`.
    #[must_use]
    pub const fn accept(policy_fingerprint: PolicyFingerprint) -> Self {
        Self {
            accepted: true,
            policy_fingerprint,
        }
    }

    /// Rejected under `policy_fingerprint`.
    #[must_use]
    pub const fn reject(policy_fingerprint: PolicyFingerprint) -> Self {
        Self {
            accepted: false,
            policy_fingerprint,
        }
    }

    /// Whether the policy accepted this exact request.
    #[must_use]
    pub const fn accepted(self) -> bool {
        self.accepted
    }

    /// Stable identity of the policy that made this decision.
    #[must_use]
    pub const fn policy_fingerprint(self) -> PolicyFingerprint {
        self.policy_fingerprint
    }
}

/// An explicit waiver grant that travels WITH its claim (schema v6).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct WaiverGrant {
    /// Stable, non-blank waiver identity (audit key).
    pub waiver_id: String,
    /// Last day (days since the Unix epoch) this waiver is valid.
    pub expiry_day: u64,
    /// Authenticator over the waiver id, expiry, and the CLAIM'S
    /// canonical bytes (replay onto another claim changes the message).
    /// Opaque here: only an injected [`WaiverVerifier`] can accept it.
    pub mac: String,
}

/// The waiver-verification CAPABILITY (injected; fs-package ships no
/// cryptography — the same fail-closed pattern as the checker's
/// [`SignatureVerifier`] and fs-ledger's waivers).
pub trait WaiverVerifier {
    /// Return a decision whose `accepted` bit is true iff `mac` authenticates
    /// the package-owned, domain-separated `message`. The message already binds
    /// the waiver id and expiry; passing them separately would let an
    /// implementation accidentally authenticate only a subset of the
    /// authorization context.
    fn verify(&self, mac: &str, message: &[u8]) -> VerificationDecision;
}

/// The in-tree default: nothing authenticates. A package whose claims
/// carry waiver origins can NEVER verify without an explicitly
/// injected capability and date context.
#[derive(Debug, Default, Clone, Copy)]
pub struct NoWaiverVerifier;

impl WaiverVerifier for NoWaiverVerifier {
    fn verify(&self, _mac: &str, _message: &[u8]) -> VerificationDecision {
        VerificationDecision::reject(fs_blake3::hash_domain(
            "fs-package:v8:policy",
            b"deny-all-waivers",
        ))
    }
}

/// Typed input to an injected source-certificate verifier.
///
/// The certificate hash is only an artifact address. Acceptance requires a
/// capability that obtains or otherwise recognizes that artifact and checks
/// that it establishes THIS exact claim under THIS package provenance.
#[derive(Debug, Clone, Copy)]
pub struct SourceCertificateRequest<'a> {
    /// Package provenance under which the certificate is being admitted.
    pub package_provenance: &'a Provenance,
    /// Recomputed root binding the certificate decision to this package.
    pub package_root: ContentHash,
    /// Stable position of the claim in the package.
    pub claim_index: usize,
    /// Claim identity.
    pub claim_id: &'a str,
    /// Human-readable assertion bound to the certificate.
    pub statement: &'a str,
    /// Source-specific claim subject excluding every external artifact address,
    /// including this source certificate's address.
    pub claim_subject_hash: ContentHash,
    /// Certified interval lower bound.
    pub lo: f64,
    /// Certified interval upper bound.
    pub hi: f64,
    /// Declared certificate producer.
    pub producer: &'a str,
    /// Parsed content address of the certificate artifact.
    pub certificate_hash: ContentHash,
    /// Inline portable witness, when the source certificate carries one.
    pub semantic_witness: Option<&'a SemanticWitness>,
}

/// Capability that re-verifies a source certificate artifact against the
/// exact typed claim request. `fs-package` deliberately has no permissive
/// built-in implementation.
pub trait SourceCertificateVerifier {
    /// True only when the addressed artifact establishes the supplied claim.
    fn verify(&self, request: &SourceCertificateRequest<'_>) -> VerificationDecision;
}

/// The in-tree source-certificate default: no artifact is trusted.
#[derive(Debug, Default, Clone, Copy)]
pub struct NoSourceCertificateVerifier;

impl SourceCertificateVerifier for NoSourceCertificateVerifier {
    fn verify(&self, _request: &SourceCertificateRequest<'_>) -> VerificationDecision {
        VerificationDecision::reject(fs_blake3::hash_domain(
            "fs-package:v8:policy",
            b"deny-all-source-certificates",
        ))
    }
}

/// Typed input to an injected anchoring-dataset verifier.
///
/// A content hash is an address, not evidence that the referenced dataset is
/// appropriate for this exact validation claim. The verifier receives every
/// semantic field whose substitution could change that decision.
#[derive(Debug, Clone, Copy)]
pub struct AnchoredSourceRequest<'a> {
    /// Package provenance under which the dataset is admitted.
    pub package_provenance: &'a Provenance,
    /// Stable position of the claim in the package.
    pub claim_index: usize,
    /// Claim identity.
    pub claim_id: &'a str,
    /// Human-readable assertion bound to the dataset decision.
    pub statement: &'a str,
    /// Exact validity regime claimed for the dataset comparison.
    pub regime: &'a fs_evidence::ValidityDomain,
    /// Dataset identity named by both the color and origin.
    pub dataset_id: &'a str,
    /// Parsed content address of the anchoring dataset.
    pub content_hash: ContentHash,
}

/// Capability that re-verifies an anchoring dataset against the exact claim.
pub trait AnchoredSourceVerifier {
    /// True only when the addressed dataset supports the supplied request.
    fn verify(&self, request: &AnchoredSourceRequest<'_>) -> VerificationDecision;
}

/// The in-tree anchoring default: no dataset address is trusted by itself.
#[derive(Debug, Default, Clone, Copy)]
pub struct NoAnchoredSourceVerifier;

impl AnchoredSourceVerifier for NoAnchoredSourceVerifier {
    fn verify(&self, _request: &AnchoredSourceRequest<'_>) -> VerificationDecision {
        VerificationDecision::reject(fs_blake3::hash_domain(
            "fs-package:v8:policy",
            b"deny-all-anchored-sources",
        ))
    }
}

/// Typed input to an injected falsifier-artifact verifier.
#[derive(Debug, Clone, Copy)]
pub struct FalsifierRequest<'a> {
    /// Package provenance under which the falsifier is admitted.
    pub package_provenance: &'a Provenance,
    /// Recomputed package root, binding sibling context and package identity.
    pub package_root: ContentHash,
    /// Stable claim position.
    pub claim_index: usize,
    /// Claim identity and human assertion targeted by the falsifier.
    pub claim_id: &'a str,
    /// Human-readable claim assertion.
    pub statement: &'a str,
    /// Exact declared color of the target claim.
    pub color: &'a fs_evidence::Color,
    /// Exact declared origin of the target claim.
    pub origin: &'a ClaimOrigin,
    /// Domain-separated claim subject hash excluding external artifact
    /// addresses and waiver MAC bytes, avoiding content-address fixed points.
    pub claim_subject_hash: ContentHash,
    /// Stable record position within the claim.
    pub falsifier_index: usize,
    /// Registered falsifier identity.
    pub name: &'a str,
    /// Number of represented adversarial attempts.
    pub attempts: u64,
    /// Whether the artifact refuted the claim.
    pub refuted: bool,
    /// Human-readable outcome summary.
    pub detail: &'a str,
    /// Parsed content address of the falsifier artifact.
    pub artifact_hash: ContentHash,
}

/// Capability that authenticates a falsifier artifact against its exact claim.
pub trait FalsifierVerifier {
    /// Atomic acceptance and policy identity for this exact request.
    fn verify(&self, request: &FalsifierRequest<'_>) -> VerificationDecision;
}

/// Typed input to an injected derivation-artifact verifier.
#[derive(Debug, Clone, Copy)]
pub struct DerivationRequest<'a> {
    /// Package provenance and root binding the complete sibling context.
    pub package_provenance: &'a Provenance,
    /// Recomputed package root.
    pub package_root: ContentHash,
    /// Exact child identity, assertion, and declared color.
    pub claim_index: usize,
    /// Derived child identity.
    pub claim_id: &'a str,
    /// Derived child assertion.
    pub statement: &'a str,
    /// Derived child declared color.
    pub color: &'a fs_evidence::Color,
    /// Child subject hash excluding external artifact addresses and waiver MAC
    /// bytes, avoiding content-address fixed points.
    pub child_subject_hash: ContentHash,
    /// Exact attached anchor declarations, if any.
    pub anchors: &'a [crate::AnchorRecord],
    /// Exact fold operation and ordered parent identities/content hashes.
    pub op: fs_evidence::IntervalOp,
    /// Ordered parent positions in the package.
    pub parent_indices: &'a [usize],
    /// Full content hashes of the ordered parent claims.
    pub parent_claim_hashes: &'a [ContentHash],
    /// Parsed content address of the derivation proof artifact.
    pub artifact_hash: ContentHash,
}

/// Capability that authenticates a derivation artifact against the exact child
/// and ordered parents. Re-folding the color remains a separate package check.
pub trait DerivationVerifier {
    /// Atomic acceptance and policy identity for this exact request.
    fn verify(&self, request: &DerivationRequest<'_>) -> VerificationDecision;
}

/// The in-tree derivation default: a receipt and hash never authenticate themselves.
#[derive(Debug, Default, Clone, Copy)]
pub struct NoDerivationVerifier;

impl DerivationVerifier for NoDerivationVerifier {
    fn verify(&self, _request: &DerivationRequest<'_>) -> VerificationDecision {
        VerificationDecision::reject(fs_blake3::hash_domain(
            "fs-package:v8:policy",
            b"deny-all-derivations",
        ))
    }
}

/// The in-tree falsifier default: an artifact address never authenticates itself.
#[derive(Debug, Default, Clone, Copy)]
pub struct NoFalsifierVerifier;

impl FalsifierVerifier for NoFalsifierVerifier {
    fn verify(&self, _request: &FalsifierRequest<'_>) -> VerificationDecision {
        VerificationDecision::reject(fs_blake3::hash_domain(
            "fs-package:v8:policy",
            b"deny-all-falsifiers",
        ))
    }
}

/// Domain-separated purpose of one detached-signature decision.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SignaturePurpose {
    /// Policy-authenticated attestation over the package root only.
    PackageRootAttestation,
    /// Explicit release-gate approval under one checker protocol, expected
    /// root, and exact scientific admission context. This cannot be substituted
    /// by a generic root attestation or replayed under different trust policies.
    ReleaseApproval {
        /// Independently distributed checker protocol.
        checker_protocol: u32,
        /// Root expected by the release gate.
        expected_root: ContentHash,
        /// Domain-separated digest over non-signature policy fingerprints,
        /// waiver clock, admissions, and the compact waiver graph.
        admission_context: ContentHash,
        /// Caller-supplied semantic-checker/plugin-set context.
        semantic_context: ContentHash,
    },
}

/// Requested signature decision before `fs-package` computes the concrete
/// scientific admission context.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SignatureIntent {
    /// Authenticate the package root only.
    PackageRootAttestation,
    /// Authenticate a release decision under an explicit checker protocol and
    /// expected root and caller-supplied semantic-checker context.
    /// `fs-package` adds the computed admission-context digest.
    ReleaseApproval {
        /// Independently distributed checker protocol.
        checker_protocol: u32,
        /// Root expected by the release gate.
        expected_root: ContentHash,
        /// Exact semantic-checker/plugin-set context approved by the caller.
        semantic_context: ContentHash,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct SignatureSubjectFields {
    package_root: ContentHash,
    purpose_tag: u8,
    checker_protocol: Option<u32>,
    expected_root: Option<ContentHash>,
    admission_context: Option<ContentHash>,
    semantic_context: Option<ContentHash>,
}

impl SignatureSubjectFields {
    fn from_inputs(package_root: ContentHash, purpose: SignaturePurpose) -> Self {
        match purpose {
            SignaturePurpose::PackageRootAttestation => Self {
                package_root,
                purpose_tag: 0,
                checker_protocol: None,
                expected_root: None,
                admission_context: None,
                semantic_context: None,
            },
            SignaturePurpose::ReleaseApproval {
                checker_protocol,
                expected_root,
                admission_context,
                semantic_context,
            } => Self {
                package_root,
                purpose_tag: 1,
                checker_protocol: Some(checker_protocol),
                expected_root: Some(expected_root),
                admission_context: Some(admission_context),
                semantic_context: Some(semantic_context),
            },
        }
    }
}

#[allow(dead_code)]
fn classify_signature_subject_identity_fields(fields: &SignatureSubjectFields) {
    let SignatureSubjectFields {
        package_root,
        purpose_tag,
        checker_protocol,
        expected_root,
        admission_context,
        semantic_context,
    } = fields;
    let _ = (
        package_root,
        purpose_tag,
        checker_protocol,
        expected_root,
        admission_context,
        semantic_context,
    );
}

/// Typed signature subject.
#[derive(Debug, Clone, Copy)]
pub struct SignatureRequest<'a> {
    /// Recomputed package root.
    pub package_root: ContentHash,
    /// Detached signature bytes.
    pub signature: &'a str,
    /// Gate/domain purpose for this authentication decision.
    pub purpose: SignaturePurpose,
}

impl SignatureRequest<'_> {
    /// Canonical domain-separated digest that signature bytes must authenticate.
    #[must_use]
    pub fn subject_hash(&self) -> ContentHash {
        signature_subject_hash(self.package_root, self.purpose)
    }
}

/// Canonical digest for producing or verifying detached signature bytes.
///
/// Release subjects include separate scientific admission and semantic-checker
/// context digests, preventing replay under different verifier policies,
/// waiver time, admission decisions, or plugin sets.
#[must_use]
pub fn signature_subject_hash(package_root: ContentHash, purpose: SignaturePurpose) -> ContentHash {
    signature_subject_hash_with_domain(
        &SignatureSubjectFields::from_inputs(package_root, purpose),
        SIGNATURE_SUBJECT_IDENTITY_DOMAIN,
    )
}

fn signature_subject_hash_with_domain(
    fields: &SignatureSubjectFields,
    domain: &str,
) -> ContentHash {
    let mut subject = Vec::with_capacity(140);
    subject.extend_from_slice(fields.package_root.as_bytes());
    match (
        fields.purpose_tag,
        fields.checker_protocol,
        fields.expected_root,
        fields.admission_context,
        fields.semantic_context,
    ) {
        (0, None, None, None, None) => {
            subject.extend_from_slice(b"package-root-attestation");
        }
        (
            1,
            Some(checker_protocol),
            Some(expected_root),
            Some(admission_context),
            Some(semantic_context),
        ) => {
            subject.extend_from_slice(b"release-approval");
            subject.extend_from_slice(&checker_protocol.to_le_bytes());
            subject.extend_from_slice(expected_root.as_bytes());
            subject.extend_from_slice(admission_context.as_bytes());
            subject.extend_from_slice(semantic_context.as_bytes());
        }
        invalid => {
            subject.extend_from_slice(b"invalid-purpose-shape");
            subject.push(invalid.0);
            if let Some(value) = invalid.1 {
                subject.extend_from_slice(&value.to_le_bytes());
            }
            for value in [invalid.2, invalid.3, invalid.4] {
                match value {
                    Some(value) => {
                        subject.push(1);
                        subject.extend_from_slice(value.as_bytes());
                    }
                    None => subject.push(0),
                }
            }
        }
    }
    fs_blake3::hash_domain(domain, &subject)
}

/// Admit a retained signature-subject digest only under the exact schema
/// version and fixed-width binary transport.
#[must_use]
pub fn admit_retained_signature_subject_hash(version: u32, bytes: &[u8]) -> Option<ContentHash> {
    crate::admit_retained_content_hash(version, SIGNATURE_SUBJECT_IDENTITY_VERSION, bytes)
}

/// The signature-verification capability. `fs-package` deliberately ships no
/// cryptography; callers inject the policy used to authenticate the exact
/// typed signature subject.
pub trait SignatureVerifier {
    /// Atomic acceptance and policy identity for this exact request.
    fn verify(&self, request: &SignatureRequest<'_>) -> VerificationDecision;
}

/// The in-tree signature default: no signature authenticates.
#[derive(Debug, Default, Clone, Copy)]
pub struct NoSignatureVerifier;

impl SignatureVerifier for NoSignatureVerifier {
    fn verify(&self, _request: &SignatureRequest<'_>) -> VerificationDecision {
        VerificationDecision::reject(fs_blake3::hash_domain(
            "fs-package:v8:policy",
            b"deny-all-signatures",
        ))
    }
}

/// Signature verifier plus its explicit purpose context.
#[derive(Clone, Copy)]
pub struct SignatureVerification<'a> {
    /// Authenticator implementation.
    pub verifier: &'a dyn SignatureVerifier,
    /// Domain/gate intent; package verification materializes the concrete
    /// purpose after scientific admission.
    pub intent: SignatureIntent,
}

/// Waiver authentication capability plus its explicit clock context.
#[derive(Clone, Copy)]
pub struct WaiverVerification<'a> {
    /// Authenticator implementation.
    pub verifier: &'a dyn WaiverVerifier,
    /// Current day, as days since the Unix epoch.
    pub today_day: u64,
}

/// External capabilities available for one package-verification decision.
/// Missing capabilities fail closed only for origin kinds that require them.
#[derive(Clone, Copy)]
pub struct VerificationCapabilities<'a> {
    /// Source-certificate artifact verifier.
    pub source_certificates: Option<&'a dyn SourceCertificateVerifier>,
    /// Anchoring-dataset artifact verifier.
    pub anchored_sources: Option<&'a dyn AnchoredSourceVerifier>,
    /// Falsifier-artifact verifier.
    pub falsifiers: Option<&'a dyn FalsifierVerifier>,
    /// Derivation-artifact verifier.
    pub derivations: Option<&'a dyn DerivationVerifier>,
    /// Waiver authenticator and clock context.
    pub waivers: Option<WaiverVerification<'a>>,
    /// Detached-signature verifier.
    pub signatures: Option<SignatureVerification<'a>>,
}

impl<'a> VerificationCapabilities<'a> {
    /// Deny every external source, anchor, falsifier, derivation, waiver, and
    /// signature capability.
    #[must_use]
    pub const fn deny_all() -> Self {
        Self {
            source_certificates: None,
            anchored_sources: None,
            falsifiers: None,
            derivations: None,
            waivers: None,
            signatures: None,
        }
    }

    /// Install an anchoring-dataset verification capability.
    #[must_use]
    pub const fn with_anchored_sources(mut self, verifier: &'a dyn AnchoredSourceVerifier) -> Self {
        self.anchored_sources = Some(verifier);
        self
    }

    /// Install a falsifier-artifact verification capability.
    #[must_use]
    pub const fn with_falsifiers(mut self, verifier: &'a dyn FalsifierVerifier) -> Self {
        self.falsifiers = Some(verifier);
        self
    }

    /// Install a derivation-artifact verification capability.
    #[must_use]
    pub const fn with_derivations(mut self, verifier: &'a dyn DerivationVerifier) -> Self {
        self.derivations = Some(verifier);
        self
    }

    /// Install a source-certificate verification capability.
    #[must_use]
    pub const fn with_source_certificates(
        mut self,
        verifier: &'a dyn SourceCertificateVerifier,
    ) -> Self {
        self.source_certificates = Some(verifier);
        self
    }

    /// Install a waiver verifier together with the decision's clock context.
    #[must_use]
    pub const fn with_waivers(mut self, verifier: &'a dyn WaiverVerifier, today_day: u64) -> Self {
        self.waivers = Some(WaiverVerification {
            verifier,
            today_day,
        });
        self
    }

    /// Install a detached-signature verification capability.
    #[must_use]
    pub const fn with_signatures(mut self, verifier: &'a dyn SignatureVerifier) -> Self {
        self.signatures = Some(SignatureVerification {
            verifier,
            intent: SignatureIntent::PackageRootAttestation,
        });
        self
    }

    /// Install signature verification for explicit release approval.
    #[must_use]
    pub const fn with_release_signatures(
        mut self,
        verifier: &'a dyn SignatureVerifier,
        checker_protocol: u32,
        expected_root: ContentHash,
        semantic_context: ContentHash,
    ) -> Self {
        self.signatures = Some(SignatureVerification {
            verifier,
            intent: SignatureIntent::ReleaseApproval {
                checker_protocol,
                expected_root,
                semantic_context,
            },
        });
        self
    }
}

impl Default for VerificationCapabilities<'_> {
    fn default() -> Self {
        Self::deny_all()
    }
}

pub(crate) fn is_placeholder_token(text: &str) -> bool {
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
    .any(|placeholder| text.eq_ignore_ascii_case(placeholder))
}

/// Reject identities whose canonical spelling would be ambiguous or
/// meaningless. Human-readable descriptions use a separate, looser policy.
pub(crate) fn identity_reason(text: &str) -> Option<&'static str> {
    let trimmed = text.trim();
    if trimmed.is_empty() {
        Some("blank")
    } else if trimmed != text {
        Some("surrounding-whitespace")
    } else if text.chars().any(char::is_control) {
        Some("control-character")
    } else if !text.bytes().all(|byte| {
        byte.is_ascii_alphanumeric()
            || matches!(byte, b'-' | b'_' | b'.' | b'/' | b':' | b'@' | b'+' | b'=')
    }) {
        Some("invalid-character")
    } else if is_placeholder_token(text) {
        Some("placeholder")
    } else {
        None
    }
}

/// Where a claim's certificate came from (schema v6). Bound into the
/// content address and re-derived by the standalone checker.
#[derive(Debug, Clone, PartialEq)]
pub enum ClaimOrigin {
    /// A named producer's certificate artifact (64-hex content hash).
    SourceCertificate {
        /// Non-blank producer identity (e.g. "fs-solver/ivp-cert").
        producer: String,
        /// Canonical 64-hex lowercase content hash of the certificate.
        certificate_hash: String,
    },
    /// The validated color's reference dataset, by id + content hash.
    AnchoredSource {
        /// Must equal the color's named dataset exactly.
        dataset_id: String,
        /// Canonical 64-hex lowercase content hash of the dataset.
        content_hash: String,
    },
    /// The estimated color's estimator identity.
    EstimatedSource {
        /// Must equal the color's estimator string exactly.
        estimator: String,
    },
    /// Derived from earlier claims: the composition receipt IS the
    /// origin (parents by index, fold op) — re-run by the checker.
    Derived,
    /// An explicit, expiring, MAC'd waiver (see [`WaiverGrant`]).
    AuthenticatedWaiver(WaiverGrant),
}

impl ClaimOrigin {
    /// Stable kind tag for hashing, JSON, and refusal messages.
    #[must_use]
    pub const fn kind(&self) -> &'static str {
        match self {
            ClaimOrigin::SourceCertificate { .. } => "source-certificate",
            ClaimOrigin::AnchoredSource { .. } => "anchored-source",
            ClaimOrigin::EstimatedSource { .. } => "estimated-source",
            ClaimOrigin::Derived => "derived",
            ClaimOrigin::AuthenticatedWaiver(_) => "authenticated-waiver",
        }
    }

    /// The canonical atom sequence bound into the claim's content
    /// hash (length-prefixed strings via the caller's `push_atom`
    /// discipline; this returns the ordered raw parts).
    #[must_use]
    pub fn canonical_parts(&self) -> Vec<String> {
        match self {
            ClaimOrigin::SourceCertificate {
                producer,
                certificate_hash,
            } => vec![
                self.kind().to_string(),
                producer.clone(),
                certificate_hash.clone(),
            ],
            ClaimOrigin::AnchoredSource {
                dataset_id,
                content_hash,
            } => vec![
                self.kind().to_string(),
                dataset_id.clone(),
                content_hash.clone(),
            ],
            ClaimOrigin::EstimatedSource { estimator } => {
                vec![self.kind().to_string(), estimator.clone()]
            }
            ClaimOrigin::Derived => vec![self.kind().to_string()],
            ClaimOrigin::AuthenticatedWaiver(grant) => vec![
                self.kind().to_string(),
                grant.waiver_id.clone(),
                grant.expiry_day.to_string(),
                grant.mac.clone(),
            ],
        }
    }
}

/// A structured origin-validation refusal (field-level, teaching).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OriginError {
    /// Which claim.
    pub claim: String,
    /// The refusal.
    pub why: String,
}

impl fmt::Display for OriginError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "claim '{}': {}", self.claim, self.why)
    }
}

impl core::error::Error for OriginError {}

/// Shape-level validation shared by construction and parsing: non-blank
/// identities, canonical 64-hex hashes where required. Color-class
/// consistency and re-derivation live with the package verifier (they
/// need the claim's color, siblings, and the injected capabilities).
///
/// # Errors
/// [`OriginError`] naming the field.
pub fn validate_origin_shape(
    claim_id: &str,
    origin: &ClaimOrigin,
    is_canonical_hash: &dyn Fn(&str) -> bool,
) -> Result<(), OriginError> {
    let refuse = |why: String| {
        Err(OriginError {
            claim: claim_id.to_string(),
            why,
        })
    };
    match origin {
        ClaimOrigin::SourceCertificate {
            producer,
            certificate_hash,
        } => {
            if let Some(reason) = identity_reason(producer) {
                return refuse(format!(
                    "source-certificate origin has an invalid producer ({reason})"
                ));
            }
            if !is_canonical_hash(certificate_hash) {
                return refuse(
                    "source-certificate origin needs a canonical 64-hex certificate hash"
                        .to_string(),
                );
            }
            Ok(())
        }
        ClaimOrigin::AnchoredSource {
            dataset_id,
            content_hash,
        } => {
            if let Some(reason) = color_leaf_identity_reason(dataset_id) {
                return refuse(format!(
                    "anchored-source origin has an invalid dataset id ({reason})"
                ));
            }
            if !is_canonical_hash(content_hash) {
                return refuse(
                    "anchored-source origin needs a canonical 64-hex dataset hash".to_string(),
                );
            }
            Ok(())
        }
        ClaimOrigin::EstimatedSource { estimator } => {
            if let Some(reason) = color_leaf_identity_reason(estimator) {
                return refuse(format!(
                    "estimated-source origin has an invalid estimator ({reason})"
                ));
            }
            Ok(())
        }
        ClaimOrigin::Derived => Ok(()),
        ClaimOrigin::AuthenticatedWaiver(grant) => {
            if let Some(reason) = identity_reason(&grant.waiver_id) {
                return refuse(format!("waiver origin has an invalid waiver id ({reason})"));
            }
            if grant.mac.trim().is_empty() {
                return refuse("waiver origin has a blank authenticator".to_string());
            }
            Ok(())
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hex64() -> String {
        "0123456789abcdef".repeat(4)
    }

    fn canonical(h: &str) -> bool {
        h.len() == 64
            && h.bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
    }

    #[test]
    fn shape_validation_fails_closed_per_field() {
        let ok = |o: &ClaimOrigin| validate_origin_shape("c", o, &canonical).is_ok();
        assert!(ok(&ClaimOrigin::SourceCertificate {
            producer: "fs-solver/ivp".to_string(),
            certificate_hash: hex64(),
        }));
        assert!(!ok(&ClaimOrigin::SourceCertificate {
            producer: "  ".to_string(),
            certificate_hash: hex64(),
        }));
        assert!(!ok(&ClaimOrigin::SourceCertificate {
            producer: "p".to_string(),
            certificate_hash: "deadbeef".to_string(),
        }));
        assert!(!ok(&ClaimOrigin::AnchoredSource {
            dataset_id: String::new(),
            content_hash: hex64(),
        }));
        assert!(!ok(&ClaimOrigin::EstimatedSource {
            estimator: " ".to_string(),
        }));
        assert!(ok(&ClaimOrigin::Derived));
        assert!(!ok(&ClaimOrigin::AuthenticatedWaiver(WaiverGrant {
            waiver_id: "w1".to_string(),
            expiry_day: 20_000,
            mac: "  ".to_string(),
        })));
    }

    #[test]
    fn canonical_parts_are_kind_prefixed_and_distinct() {
        let a = ClaimOrigin::EstimatedSource {
            estimator: "surrogate-v2".to_string(),
        };
        let b = ClaimOrigin::SourceCertificate {
            producer: "surrogate-v2".to_string(),
            certificate_hash: hex64(),
        };
        assert_eq!(a.canonical_parts()[0], "estimated-source");
        assert_ne!(a.canonical_parts(), b.canonical_parts());
        // The waiver's expiry and mac are bound (tamper moves the parts).
        let w1 = ClaimOrigin::AuthenticatedWaiver(WaiverGrant {
            waiver_id: "w".to_string(),
            expiry_day: 1,
            mac: "m".to_string(),
        });
        let w2 = ClaimOrigin::AuthenticatedWaiver(WaiverGrant {
            waiver_id: "w".to_string(),
            expiry_day: 2,
            mac: "m".to_string(),
        });
        assert_ne!(w1.canonical_parts(), w2.canonical_parts());
    }

    #[test]
    fn the_default_waiver_verifier_refuses_everything() {
        let grant = WaiverGrant {
            waiver_id: "w1".to_string(),
            expiry_day: u64::MAX,
            mac: "anything".to_string(),
        };
        assert!(!NoWaiverVerifier.verify(&grant.mac, b"message").accepted());
    }

    #[test]
    fn the_default_source_verifier_refuses_everything() {
        let provenance = Provenance::new("v", "lock");
        let request = SourceCertificateRequest {
            package_provenance: &provenance,
            package_root: ContentHash([0; 32]),
            claim_index: 0,
            claim_id: "c",
            statement: "bounded",
            claim_subject_hash: ContentHash([0; 32]),
            lo: 0.0,
            hi: 1.0,
            producer: "solver/cert",
            certificate_hash: ContentHash([0; 32]),
            semantic_witness: None,
        };
        assert!(!NoSourceCertificateVerifier.verify(&request).accepted());
    }

    #[test]
    fn the_default_anchor_and_signature_verifiers_refuse_everything() {
        let provenance = Provenance::new("v", "lock");
        let regime = fs_evidence::ValidityDomain::unconstrained().with("Re", 1.0, 2.0);
        let request = AnchoredSourceRequest {
            package_provenance: &provenance,
            claim_index: 0,
            claim_id: "c",
            statement: "matches",
            regime: &regime,
            dataset_id: "dataset",
            content_hash: ContentHash([0; 32]),
        };
        assert!(!NoAnchoredSourceVerifier.verify(&request).accepted());
        let signature = SignatureRequest {
            package_root: ContentHash([0; 32]),
            signature: "signature",
            purpose: SignaturePurpose::PackageRootAttestation,
        };
        assert!(!NoSignatureVerifier.verify(&signature).accepted());
    }

    #[test]
    fn signature_subject_identity_fields_move_independently() {
        let fields = SignatureSubjectFields {
            package_root: ContentHash([1; 32]),
            purpose_tag: 1,
            checker_protocol: Some(4),
            expected_root: Some(ContentHash([2; 32])),
            admission_context: Some(ContentHash([3; 32])),
            semantic_context: Some(ContentHash([4; 32])),
        };
        let baseline =
            signature_subject_hash_with_domain(&fields, SIGNATURE_SUBJECT_IDENTITY_DOMAIN);
        assert_ne!(
            baseline,
            signature_subject_hash_with_domain(
                &fields,
                "fs-package:v9:alternate-signature-subject"
            ),
            "digest-domain",
        );

        let mut changed = fields;
        changed.package_root = ContentHash([9; 32]);
        assert_ne!(
            baseline,
            signature_subject_hash_with_domain(&changed, SIGNATURE_SUBJECT_IDENTITY_DOMAIN),
            "package-root",
        );

        let mut changed = fields;
        changed.purpose_tag = 0;
        assert_ne!(
            baseline,
            signature_subject_hash_with_domain(&changed, SIGNATURE_SUBJECT_IDENTITY_DOMAIN),
            "purpose-tag",
        );

        let mut changed = fields;
        changed.checker_protocol = None;
        assert_ne!(
            baseline,
            signature_subject_hash_with_domain(&changed, SIGNATURE_SUBJECT_IDENTITY_DOMAIN),
            "checker-protocol-presence-and-value",
        );

        let mut changed = fields;
        changed.expected_root = Some(ContentHash([9; 32]));
        assert_ne!(
            baseline,
            signature_subject_hash_with_domain(&changed, SIGNATURE_SUBJECT_IDENTITY_DOMAIN),
            "expected-root-presence-and-value",
        );

        let mut changed = fields;
        changed.admission_context = None;
        assert_ne!(
            baseline,
            signature_subject_hash_with_domain(&changed, SIGNATURE_SUBJECT_IDENTITY_DOMAIN),
            "admission-context-presence-and-value",
        );

        let mut changed = fields;
        changed.semantic_context = Some(ContentHash([9; 32]));
        assert_ne!(
            baseline,
            signature_subject_hash_with_domain(&changed, SIGNATURE_SUBJECT_IDENTITY_DOMAIN),
            "semantic-context-presence-and-value",
        );
    }

    #[test]
    fn signature_subject_excludes_detached_signature_bytes() {
        let purpose = SignaturePurpose::ReleaseApproval {
            checker_protocol: 4,
            expected_root: ContentHash([2; 32]),
            admission_context: ContentHash([3; 32]),
            semantic_context: ContentHash([4; 32]),
        };
        let left = SignatureRequest {
            package_root: ContentHash([1; 32]),
            signature: "detached-a",
            purpose,
        };
        let right = SignatureRequest {
            package_root: left.package_root,
            signature: "detached-b",
            purpose,
        };
        assert_eq!(left.subject_hash(), right.subject_hash());
    }

    #[test]
    fn machine_identities_reject_placeholders_and_padding() {
        assert_eq!(identity_reason("todo"), Some("placeholder"));
        assert_eq!(identity_reason(" producer"), Some("surrounding-whitespace"));
        assert_eq!(identity_reason("producer"), None);
        assert_eq!(
            identity_reason("prod\u{202e}recudor"),
            Some("invalid-character")
        );
        assert_eq!(identity_reason("prod\0ucer"), Some("control-character"));
    }
}
