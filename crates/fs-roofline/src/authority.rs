//! Promotion-authority verification (bead fz2.7): WHO authorized a
//! baseline, checkable — not just a free-text `promoted_by` annotation.
//!
//! A [`crate::BaselineAxes`] record is tamper-EVIDENT (its canonical
//! content hash moves when any field moves) but was previously only
//! operator-TRUSTED: nothing proved the named operator authorized the
//! promotion, and a locally editable store was an explicit trust root.
//! This module closes that gap with the workspace's capability pattern:
//!
//! - A [`PromotionAttestation`] travels with the baseline: a KEY
//!   IDENTITY plus an opaque signature over the record's
//!   [`crate::BaselineAxes::content_hash`] — which already binds the
//!   canonical schema/domain hash, the sorted source receipt
//!   identities, the band/drift policy, the machine identity, and the
//!   promotion time. Signing the hash signs all of them.
//! - A [`PromotionAuthorityVerifier`] CAPABILITY interprets signatures
//!   and answers with a typed [`KeyVerdict`] (authorized / wrong
//!   signature / unknown key / revoked key). The in-tree default,
//!   [`NoPromotionAuthority`], refuses everything — verification only
//!   ever happens against an explicitly injected authority.
//! - `StaticKeyRegistry` is the unit-test-only deterministic registry: it
//!   tracks which key identities are authorized or REVOKED and checks
//!   domain-separated keyed-hash tags. It is an operator-governed
//!   registry with tamper-evident tags, NOT unforgeable cryptography —
//!   anyone with this code can mint a tag, so unforgeability requires
//!   an external verifier implementation (the no-crypto no-claim,
//!   exactly like fs-checker signatures and fs-package waivers).
//!   Rotation = authorize a new key id; revocation = mark the old one
//!   revoked. An old attestation then re-verifies as
//!   [`KeyVerdict::RevokedKey`], while an explicit attestation of the SAME
//!   immutable record under the new key is a valid re-endorsement. No record
//!   mutation or fictitious re-promotion is required.

use fs_blake3::{ContentHash, hash_domain};
use std::collections::{BTreeMap, BTreeSet};

/// Domain tag for registry keyed-hash attestation tags.
pub const PROMOTION_AUTHORITY_DOMAIN: &str = "frankensim.fs-roofline.promotion-authority.v1";

/// Domain for the exact authority policy bound into every decision.
pub const PROMOTION_AUTHORITY_POLICY_DOMAIN: &str =
    "frankensim.fs-roofline.promotion-authority-policy.v1";

/// Semantic version of the configured promotion-authority policy identity.
pub const PROMOTION_AUTHORITY_POLICY_IDENTITY_VERSION: u32 = 1;

/// Owner-local authority-policy declaration consumed by `xtask check-identities`.
#[allow(dead_code)]
pub const PROMOTION_AUTHORITY_POLICY_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-roofline:promotion-authority-policy",
    "version_const=PROMOTION_AUTHORITY_POLICY_IDENTITY_VERSION",
    "version=1",
    "domain=frankensim.fs-roofline.promotion-authority-policy.v1",
    "domain_const=PROMOTION_AUTHORITY_POLICY_DOMAIN",
    "encoder=promotion_authority_policy_receipt",
    "encoder_helpers=promotion_authority_policy_receipt_with_domain",
    "schema_functions=ConfiguredPromotionAuthority::from_text,ConfiguredPromotionAuthority::policy_receipt,ConfiguredPromotionAuthority::verify,PromotionAuthorityDecision::new,PromotionAuthorityDecision::policy_receipt,crates/fs-blake3/src/lib.rs#hash_domain",
    "schema_constants=PROMOTION_AUTHORITY_POLICY_IDENTITY_VERSION,PROMOTION_AUTHORITY_POLICY_DOMAIN,MAX_PROMOTION_AUTHORITY_POLICY_BYTES,MAX_PROMOTION_AUTHORITY_POLICY_ENTRIES,MAX_PROMOTION_AUTHORITY_POLICY_LINE_BYTES,MAX_PROMOTION_AUTHORITY_FIELD_BYTES",
    "schema_dependencies=none",
    "digest=fs-blake3",
    "encoding=canonical-transport-exact-bits",
    "sources=PromotionAuthorityPolicyIdentityInput",
    "source_fields=PromotionAuthorityPolicyIdentityInput.canonical_bytes:semantic",
    "source_bindings=PromotionAuthorityPolicyIdentityInput.canonical_bytes>canonical-policy-bytes",
    "external_semantic_fields=digest-domain,identity-version",
    "semantic_fields=digest-domain,identity-version,canonical-policy-bytes",
    "excluded_fields=none",
    "consumers=ConfiguredPromotionAuthority::from_text,ConfiguredPromotionAuthority::policy_receipt,ConfiguredPromotionAuthority::verify,PromotionAuthorityDecision::policy_receipt,AttestedBaselineStore::policy_for_run,AttestedAxisBaselinePolicy,AxisAdmissionSnapshot::receipt_json",
    "mutations=digest-domain:crates/fs-roofline/src/authority.rs#configured_authority_policy_identity_fields_move_independently,identity-version:crates/fs-roofline/src/authority.rs#configured_authority_policy_identity_versions_fail_closed,canonical-policy-bytes:crates/fs-roofline/src/authority.rs#configured_authority_policy_identity_fields_move_independently",
    "nonsemantic_mutations=none",
    "field_guard=classify_promotion_authority_policy_identity_fields",
    "transport_guard=ConfiguredPromotionAuthority::from_text",
    "version_guard=crates/fs-roofline/src/authority.rs#configured_authority_policy_identity_versions_fail_closed",
    "coupling_surface=fs-roofline:promotion-authority-policy",
];

/// Maximum accepted size of a configured authority policy.
pub const MAX_PROMOTION_AUTHORITY_POLICY_BYTES: usize = 1024 * 1024;
const MAX_PROMOTION_AUTHORITY_POLICY_LINE_BYTES: usize = 16 * 1024;
const MAX_PROMOTION_AUTHORITY_FIELD_BYTES: usize = 4096;
const MAX_PROMOTION_AUTHORITY_POLICY_ENTRIES: usize = 4096;

/// A signed promotion authorization that travels WITH its baseline.
/// The signature covers the record's content hash and nothing else —
/// editing ANY signed field (operator, axes, receipts, policy, time)
/// moves the hash and invalidates the attestation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PromotionAttestation {
    key_id: String,
    signature: String,
}

impl PromotionAttestation {
    /// Construct an attestation value. Store/policy boundaries validate both
    /// fields with [`Self::well_formed`].
    #[must_use]
    pub fn new(key_id: impl Into<String>, signature: impl Into<String>) -> PromotionAttestation {
        PromotionAttestation {
            key_id: key_id.into(),
            signature: signature.into(),
        }
    }

    /// The signing key's stable identity (bound into ledger receipts).
    #[must_use]
    pub fn key_id(&self) -> &str {
        &self.key_id
    }

    /// The opaque signature (only an injected verifier interprets it).
    #[must_use]
    pub fn signature(&self) -> &str {
        &self.signature
    }

    /// Structurally usable: both parts are bounded, non-blank, have no
    /// surrounding whitespace, and contain no control characters.
    #[must_use]
    pub fn well_formed(&self) -> bool {
        validate_authority_field("key id", &self.key_id).is_ok()
            && validate_authority_field("signature", &self.signature).is_ok()
    }
}

/// The typed answer of a promotion-authority verifier. Every variant
/// except [`KeyVerdict::Authorized`] fails closed with its own
/// teaching name.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KeyVerdict {
    /// The key is currently authorized and the signature verifies.
    Authorized,
    /// The key exists and is authorized, but the signature does not
    /// verify over this record (forged or edited record).
    WrongSignature,
    /// No such key in the authority's registry.
    UnknownKey,
    /// The key was rotated out: its old attestations need explicit
    /// re-endorsement under a currently authorized key.
    RevokedKey,
}

impl KeyVerdict {
    /// Stable name for receipts and refusal messages.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            KeyVerdict::Authorized => "authorized",
            KeyVerdict::WrongSignature => "wrong-signature",
            KeyVerdict::UnknownKey => "unknown-key",
            KeyVerdict::RevokedKey => "revoked-key",
        }
    }
}

/// One atomic answer from a promotion authority. The verdict and the exact
/// policy receipt are returned together so callers cannot accidentally bind a
/// verdict observed under one policy to the identity of another policy.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PromotionAuthorityDecision {
    verdict: KeyVerdict,
    policy_receipt: ContentHash,
}

impl PromotionAuthorityDecision {
    /// Bind `verdict` to the exact policy that produced it.
    #[must_use]
    pub const fn new(verdict: KeyVerdict, policy_receipt: ContentHash) -> Self {
        Self {
            verdict,
            policy_receipt,
        }
    }

    /// Typed authority verdict.
    #[must_use]
    pub const fn verdict(self) -> KeyVerdict {
        self.verdict
    }

    /// Content identity of the exact authority policy used for this decision.
    #[must_use]
    pub const fn policy_receipt(self) -> ContentHash {
        self.policy_receipt
    }
}

/// The promotion-authority CAPABILITY (injected; fs-roofline ships no
/// cryptography).
pub trait PromotionAuthorityVerifier {
    /// Judge `signature` by `key_id` over `message` (the baseline's
    /// content-hash bytes).
    #[must_use]
    fn verify(&self, key_id: &str, signature: &str, message: &[u8]) -> PromotionAuthorityDecision;
}

/// The fail-closed default: every key is unknown; nothing verifies.
#[derive(Debug, Default, Clone, Copy)]
pub struct NoPromotionAuthority;

impl PromotionAuthorityVerifier for NoPromotionAuthority {
    fn verify(
        &self,
        _key_id: &str,
        _signature: &str,
        _message: &[u8],
    ) -> PromotionAuthorityDecision {
        PromotionAuthorityDecision::new(
            KeyVerdict::UnknownKey,
            hash_domain(
                PROMOTION_AUTHORITY_POLICY_DOMAIN,
                b"no-promotion-authority-v1",
            ),
        )
    }
}

/// Deterministic in-tree authority registry: authorized/revoked key
/// identities with domain-separated keyed-hash tags.
///
/// Trust class (documented, honest): operator-governed and
/// tamper-EVIDENT, not unforgeable — the tag function is public, so
/// this registry proves "the record matches what a listed, unrevoked
/// key id attested" only against accidental or unprivileged edits. An
/// adversary with repo access can mint tags; unforgeable signatures
/// require an EXTERNAL [`PromotionAuthorityVerifier`] implementation.
#[derive(Debug, Default, Clone)]
#[cfg(test)]
pub struct StaticKeyRegistry {
    /// key id → revoked?
    keys: BTreeMap<String, bool>,
}

#[cfg(test)]
impl StaticKeyRegistry {
    /// An empty registry (verifies nothing).
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Authorize a key id (idempotent; un-revokes a revoked id only
    /// via this explicit call — that IS a rotation decision).
    pub fn authorize(&mut self, key_id: impl Into<String>) {
        self.keys.insert(key_id.into(), false);
    }

    /// Revoke a key id (unknown ids are recorded as revoked too — a
    /// revocation list survives registry rebuilds).
    pub fn revoke(&mut self, key_id: impl Into<String>) {
        self.keys.insert(key_id.into(), true);
    }

    /// The deterministic tag a listed key id produces over `message`
    /// (public by design — see the trust-class note above).
    #[must_use]
    pub fn tag(key_id: &str, message: &[u8]) -> String {
        let mut preimage = Vec::with_capacity(key_id.len() + 1 + message.len());
        preimage.extend_from_slice(key_id.as_bytes());
        preimage.push(0);
        preimage.extend_from_slice(message);
        hash_domain(PROMOTION_AUTHORITY_DOMAIN, &preimage).to_hex()
    }

    fn policy_receipt(&self) -> ContentHash {
        let mut preimage = b"static-key-registry-v1".to_vec();
        for (key_id, revoked) in &self.keys {
            preimage.extend_from_slice(&(key_id.len() as u64).to_le_bytes());
            preimage.extend_from_slice(key_id.as_bytes());
            preimage.push(u8::from(*revoked));
        }
        hash_domain(PROMOTION_AUTHORITY_POLICY_DOMAIN, &preimage)
    }
}

#[cfg(test)]
impl PromotionAuthorityVerifier for StaticKeyRegistry {
    fn verify(&self, key_id: &str, signature: &str, message: &[u8]) -> PromotionAuthorityDecision {
        let verdict = match self.keys.get(key_id) {
            None => KeyVerdict::UnknownKey,
            Some(true) => KeyVerdict::RevokedKey,
            Some(false) => {
                if signature == StaticKeyRegistry::tag(key_id, message) {
                    KeyVerdict::Authorized
                } else {
                    KeyVerdict::WrongSignature
                }
            }
        };
        PromotionAuthorityDecision::new(verdict, self.policy_receipt())
    }
}

/// A malformed configured promotion-authority policy.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PromotionAuthorityConfigError {
    detail: String,
}

impl PromotionAuthorityConfigError {
    /// Stable teaching diagnostic.
    #[must_use]
    pub fn detail(&self) -> &str {
        &self.detail
    }
}

impl core::fmt::Display for PromotionAuthorityConfigError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "promotion authority policy refused: {}", self.detail)
    }
}

impl core::error::Error for PromotionAuthorityConfigError {}

type ExactAuthorizations = BTreeMap<String, BTreeMap<ContentHash, String>>;

struct PromotionAuthorityPolicyIdentityInput<'a> {
    canonical_bytes: &'a [u8],
}

fn promotion_authority_policy_receipt(canonical_bytes: &[u8]) -> ContentHash {
    promotion_authority_policy_receipt_with_domain(
        PROMOTION_AUTHORITY_POLICY_DOMAIN,
        &PromotionAuthorityPolicyIdentityInput { canonical_bytes },
    )
}

fn promotion_authority_policy_receipt_with_domain(
    domain: &str,
    input: &PromotionAuthorityPolicyIdentityInput<'_>,
) -> ContentHash {
    hash_domain(domain, input.canonical_bytes)
}

#[allow(dead_code)]
fn classify_promotion_authority_policy_identity_fields(
    input: &PromotionAuthorityPolicyIdentityInput<'_>,
) {
    let PromotionAuthorityPolicyIdentityInput { canonical_bytes: _ } = input;
}

/// A bounded, immutable exact-message authority for shipped entrypoints.
///
/// Its canonical text format is sorted UTF-8 TSV with a mandatory final
/// newline for every non-empty policy:
///
/// ```text
/// authorize\t<key-id>\t<64-lowercase-hex-message>\t<opaque-signature>
/// revoke\t<key-id>
/// ```
///
/// An authorization names one exact `(key, message, signature)` tuple. A
/// revoked key cannot also have authorizations. The exact canonical bytes are
/// hashed into every [`PromotionAuthorityDecision`]. This is a protected-file
/// trust root, not a cryptographic verifier: the policy must be supplied by a
/// separately protected deployment channel.
#[derive(Debug, Clone)]
pub struct ConfiguredPromotionAuthority {
    authorized: ExactAuthorizations,
    revoked: BTreeSet<String>,
    policy_receipt: ContentHash,
}

impl ConfiguredPromotionAuthority {
    /// Parse an immutable exact-message policy.
    ///
    /// # Errors
    /// Refuses oversized, non-canonical, ambiguous, duplicated, or malformed
    /// input. Empty text is a valid deny-all configured policy.
    #[allow(clippy::too_many_lines)] // Keeping the strict grammar's refusal order local is auditable.
    pub fn from_text(text: &str) -> Result<Self, PromotionAuthorityConfigError> {
        if text.len() > MAX_PROMOTION_AUTHORITY_POLICY_BYTES {
            return Err(policy_config_error(format!(
                "policy exceeds the {MAX_PROMOTION_AUTHORITY_POLICY_BYTES}-byte bound"
            )));
        }
        if text.contains('\r') {
            return Err(policy_config_error(
                "policy must use canonical LF line endings",
            ));
        }
        if !text.is_empty() && !text.ends_with('\n') {
            return Err(policy_config_error(
                "non-empty policy must end with a newline",
            ));
        }

        let mut authorized = ExactAuthorizations::new();
        let mut revoked = BTreeSet::new();
        let mut previous: Option<&str> = None;
        let lines = text.strip_suffix('\n').unwrap_or(text);
        if !text.is_empty() && lines.is_empty() {
            return Err(policy_config_error(
                "a deny-all policy must be the canonical empty byte string",
            ));
        }
        if !lines.is_empty() {
            for (index, line) in lines.split('\n').enumerate() {
                let line_number = index + 1;
                if line.is_empty() {
                    return Err(policy_config_error(format!("line {line_number} is blank")));
                }
                if line.len() > MAX_PROMOTION_AUTHORITY_POLICY_LINE_BYTES {
                    return Err(policy_config_error(format!(
                        "line {line_number} exceeds the {MAX_PROMOTION_AUTHORITY_POLICY_LINE_BYTES}-byte bound"
                    )));
                }
                if index >= MAX_PROMOTION_AUTHORITY_POLICY_ENTRIES {
                    return Err(policy_config_error(format!(
                        "policy exceeds the {MAX_PROMOTION_AUTHORITY_POLICY_ENTRIES}-entry bound"
                    )));
                }
                if previous.is_some_and(|prior| prior >= line) {
                    return Err(policy_config_error(format!(
                        "line {line_number} is not in strict canonical byte order"
                    )));
                }
                previous = Some(line);

                let mut fields = line.split('\t');
                let record = (
                    fields.next(),
                    fields.next(),
                    fields.next(),
                    fields.next(),
                    fields.next(),
                );
                match record {
                    (Some("authorize"), Some(key_id), Some(message_hex), Some(signature), None) => {
                        validate_authority_field("key id", key_id).map_err(|detail| {
                            policy_config_error(format!("line {line_number}: {detail}"))
                        })?;
                        validate_authority_field("signature", signature).map_err(|detail| {
                            policy_config_error(format!("line {line_number}: {detail}"))
                        })?;
                        if message_hex.len() != 64
                            || !message_hex
                                .bytes()
                                .all(|byte| matches!(byte, b'0'..=b'9' | b'a'..=b'f'))
                        {
                            return Err(policy_config_error(format!(
                                "line {line_number}: message must be 64 lowercase hex characters"
                            )));
                        }
                        let message = ContentHash::from_hex(message_hex).ok_or_else(|| {
                            policy_config_error(format!(
                                "line {line_number}: malformed message hash"
                            ))
                        })?;
                        if revoked.contains(key_id) {
                            return Err(policy_config_error(format!(
                                "line {line_number}: key {key_id:?} is both authorized and revoked"
                            )));
                        }
                        let messages = authorized.entry(key_id.to_string()).or_default();
                        if messages.insert(message, signature.to_string()).is_some() {
                            return Err(policy_config_error(format!(
                                "line {line_number}: duplicate key/message authorization"
                            )));
                        }
                    }
                    (Some("revoke"), Some(key_id), None, None, None) => {
                        validate_authority_field("key id", key_id).map_err(|detail| {
                            policy_config_error(format!("line {line_number}: {detail}"))
                        })?;
                        if authorized.contains_key(key_id) {
                            return Err(policy_config_error(format!(
                                "line {line_number}: key {key_id:?} is both authorized and revoked"
                            )));
                        }
                        if !revoked.insert(key_id.to_string()) {
                            return Err(policy_config_error(format!(
                                "line {line_number}: duplicate revocation"
                            )));
                        }
                    }
                    _ => {
                        return Err(policy_config_error(format!(
                            "line {line_number} must be authorize/key/message/signature or revoke/key TSV"
                        )));
                    }
                }
            }
        }

        Ok(Self {
            authorized,
            revoked,
            policy_receipt: promotion_authority_policy_receipt(text.as_bytes()),
        })
    }

    /// Identity of the exact canonical policy bytes.
    #[must_use]
    pub const fn policy_receipt(&self) -> ContentHash {
        self.policy_receipt
    }
}

impl PromotionAuthorityVerifier for ConfiguredPromotionAuthority {
    fn verify(&self, key_id: &str, signature: &str, message: &[u8]) -> PromotionAuthorityDecision {
        let verdict = if self.revoked.contains(key_id) {
            KeyVerdict::RevokedKey
        } else if let Some(messages) = self.authorized.get(key_id) {
            let exact = ContentHash::from_slice(message)
                .and_then(|message| messages.get(&message))
                .is_some_and(|expected| expected == signature);
            if exact {
                KeyVerdict::Authorized
            } else {
                KeyVerdict::WrongSignature
            }
        } else {
            KeyVerdict::UnknownKey
        };
        PromotionAuthorityDecision::new(verdict, self.policy_receipt)
    }
}

fn policy_config_error(detail: impl Into<String>) -> PromotionAuthorityConfigError {
    PromotionAuthorityConfigError {
        detail: detail.into(),
    }
}

fn validate_authority_field(field: &str, value: &str) -> Result<(), String> {
    if value.is_empty() || value.trim() != value {
        return Err(format!(
            "{field} must be non-blank without surrounding whitespace"
        ));
    }
    if value.len() > MAX_PROMOTION_AUTHORITY_FIELD_BYTES {
        return Err(format!(
            "{field} exceeds the {MAX_PROMOTION_AUTHORITY_FIELD_BYTES}-byte bound"
        ));
    }
    if value.chars().any(char::is_control) {
        return Err(format!("{field} must not contain control characters"));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_default_authority_knows_no_keys() {
        let first = NoPromotionAuthority.verify("any", "sig", b"m");
        let second = NoPromotionAuthority.verify("other", "different", b"message");
        assert_eq!(first.verdict(), KeyVerdict::UnknownKey);
        assert_eq!(second.verdict(), KeyVerdict::UnknownKey);
        assert_eq!(first.policy_receipt(), second.policy_receipt());
    }

    #[test]
    fn registry_verdicts_cover_the_four_outcomes() {
        let mut registry = StaticKeyRegistry::new();
        registry.authorize("ops/2026-q3");
        let message = b"baseline-content-hash";
        let good = StaticKeyRegistry::tag("ops/2026-q3", message);
        assert_eq!(
            registry.verify("ops/2026-q3", &good, message).verdict(),
            KeyVerdict::Authorized
        );
        assert_eq!(
            registry.verify("ops/2026-q3", "forged", message).verdict(),
            KeyVerdict::WrongSignature
        );
        assert_eq!(
            registry.verify("ops/2026-q2", &good, message).verdict(),
            KeyVerdict::UnknownKey
        );
        registry.revoke("ops/2026-q3");
        assert_eq!(
            registry.verify("ops/2026-q3", &good, message).verdict(),
            KeyVerdict::RevokedKey
        );
        // A tag is message-bound: any other message fails.
        registry.authorize("ops/2026-q3");
        assert_eq!(
            registry
                .verify("ops/2026-q3", &good, b"other message")
                .verdict(),
            KeyVerdict::WrongSignature
        );
        // And key-bound: another key's tag never transfers.
        registry.authorize("ops/2026-q4");
        assert_eq!(
            registry.verify("ops/2026-q4", &good, message).verdict(),
            KeyVerdict::WrongSignature
        );
    }

    #[test]
    fn static_policy_receipt_moves_with_registry_policy() {
        let mut registry = StaticKeyRegistry::new();
        let empty = registry.verify("none", "sig", b"message").policy_receipt();
        registry.authorize("ops/a");
        let authorized = registry.verify("none", "sig", b"message").policy_receipt();
        registry.revoke("ops/a");
        let revoked = registry.verify("none", "sig", b"message").policy_receipt();
        assert_ne!(empty, authorized);
        assert_ne!(authorized, revoked);
    }

    #[test]
    fn configured_authority_binds_exact_messages_and_policy_bytes() {
        let message = ContentHash([0x2a; 32]);
        let text = format!(
            "authorize\tops/a\t{}\tsignature-a\nrevoke\tops/old\n",
            message.to_hex()
        );
        let authority = ConfiguredPromotionAuthority::from_text(&text).expect("canonical policy");
        let accepted = authority.verify("ops/a", "signature-a", message.as_bytes());
        assert_eq!(accepted.verdict(), KeyVerdict::Authorized);
        assert_eq!(accepted.policy_receipt(), authority.policy_receipt());
        assert_eq!(
            authority
                .verify("ops/a", "wrong", message.as_bytes())
                .verdict(),
            KeyVerdict::WrongSignature
        );
        assert_eq!(
            authority
                .verify("ops/a", "signature-a", &[0x2b; 32])
                .verdict(),
            KeyVerdict::WrongSignature
        );
        assert_eq!(
            authority
                .verify("ops/old", "anything", message.as_bytes())
                .verdict(),
            KeyVerdict::RevokedKey
        );
        assert_eq!(
            authority
                .verify("ops/unknown", "anything", message.as_bytes())
                .verdict(),
            KeyVerdict::UnknownKey
        );
        let empty = ConfiguredPromotionAuthority::from_text("").expect("empty deny policy");
        assert_ne!(empty.policy_receipt(), authority.policy_receipt());
    }

    #[test]
    fn configured_authority_policy_identity_fields_move_independently() {
        let input = PromotionAuthorityPolicyIdentityInput {
            canonical_bytes: b"revoke\tops/a\n",
        };
        let current = promotion_authority_policy_receipt(input.canonical_bytes);
        let changed_domain = promotion_authority_policy_receipt_with_domain(
            "frankensim.fs-roofline.promotion-authority-policy-shadow.v1",
            &input,
        );
        let changed_bytes = promotion_authority_policy_receipt(b"revoke\tops/b\n");
        assert_ne!(current, changed_domain, "the digest domain is semantic");
        assert_ne!(
            current, changed_bytes,
            "the exact policy bytes are semantic"
        );
    }

    #[test]
    fn configured_authority_policy_identity_versions_fail_closed() {
        assert_eq!(PROMOTION_AUTHORITY_POLICY_IDENTITY_VERSION, 1);
        assert!(PROMOTION_AUTHORITY_POLICY_DOMAIN.ends_with(".v1"));
        let input = PromotionAuthorityPolicyIdentityInput {
            canonical_bytes: b"revoke\tops/a\n",
        };
        let current = promotion_authority_policy_receipt(input.canonical_bytes);
        let future = promotion_authority_policy_receipt_with_domain(
            "frankensim.fs-roofline.promotion-authority-policy.v2",
            &input,
        );
        let configured = ConfiguredPromotionAuthority::from_text("revoke\tops/a\n")
            .expect("canonical configured authority");
        assert_eq!(configured.policy_receipt(), current);
        assert_ne!(configured.policy_receipt(), future);
    }

    #[test]
    fn configured_authority_parser_refuses_ambiguous_encodings() {
        let message = ContentHash([0x2a; 32]).to_hex();
        assert!(ConfiguredPromotionAuthority::from_text("\n").is_err());
        assert!(ConfiguredPromotionAuthority::from_text("revoke\tops/a").is_err());
        assert!(
            ConfiguredPromotionAuthority::from_text(&format!(
                "authorize\tops/a\t{}\tsig\n",
                message.to_uppercase()
            ))
            .is_err()
        );
        assert!(
            ConfiguredPromotionAuthority::from_text(&format!(
                "revoke\tops/z\nauthorize\tops/a\t{message}\tsig\n"
            ))
            .is_err(),
            "lines must be in canonical byte order"
        );
        assert!(
            ConfiguredPromotionAuthority::from_text(&format!(
                "authorize\tops/a\t{message}\tsig\nrevoke\tops/a\n"
            ))
            .is_err(),
            "a key cannot be both active and revoked"
        );
        let duplicate = format!("authorize\tops/a\t{message}\tsig\n");
        assert!(
            ConfiguredPromotionAuthority::from_text(&format!("{duplicate}{duplicate}")).is_err(),
            "duplicate configuration entries are ambiguous"
        );
    }
}
