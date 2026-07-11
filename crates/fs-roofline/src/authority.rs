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
//! - [`StaticKeyRegistry`] is the in-tree deterministic registry: it
//!   tracks which key identities are authorized or REVOKED and checks
//!   domain-separated keyed-hash tags. It is an operator-governed
//!   registry with tamper-evident tags, NOT unforgeable cryptography —
//!   anyone with this code can mint a tag, so unforgeability requires
//!   an external verifier implementation (the no-crypto no-claim,
//!   exactly like fs-checker signatures and fs-package waivers).
//!   Rotation = authorize a new key id; revocation = mark the old one
//!   revoked; previously attested records then re-verify as
//!   [`KeyVerdict::RevokedKey`] and demand re-promotion.

use fs_blake3::hash_domain;
use std::collections::BTreeMap;

/// Domain tag for registry keyed-hash attestation tags.
pub const PROMOTION_AUTHORITY_DOMAIN: &str = "frankensim.fs-roofline.promotion-authority.v1";

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
    /// A new attestation; both parts must be non-blank.
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

    /// Structurally usable: both parts non-blank.
    #[must_use]
    pub fn well_formed(&self) -> bool {
        !self.key_id.trim().is_empty() && !self.signature.trim().is_empty()
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
    /// The key was rotated out: records it signed need re-promotion
    /// under a currently authorized key.
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

/// The promotion-authority CAPABILITY (injected; fs-roofline ships no
/// cryptography).
pub trait PromotionAuthorityVerifier {
    /// Judge `signature` by `key_id` over `message` (the baseline's
    /// content-hash bytes).
    fn verify(&self, key_id: &str, signature: &str, message: &[u8]) -> KeyVerdict;
}

/// The fail-closed default: every key is unknown; nothing verifies.
#[derive(Debug, Default, Clone, Copy)]
pub struct NoPromotionAuthority;

impl PromotionAuthorityVerifier for NoPromotionAuthority {
    fn verify(&self, _key_id: &str, _signature: &str, _message: &[u8]) -> KeyVerdict {
        KeyVerdict::UnknownKey
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
pub struct StaticKeyRegistry {
    /// key id → revoked?
    keys: BTreeMap<String, bool>,
}

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
}

impl PromotionAuthorityVerifier for StaticKeyRegistry {
    fn verify(&self, key_id: &str, signature: &str, message: &[u8]) -> KeyVerdict {
        match self.keys.get(key_id) {
            None => KeyVerdict::UnknownKey,
            Some(true) => KeyVerdict::RevokedKey,
            Some(false) => {
                if signature == StaticKeyRegistry::tag(key_id, message) {
                    KeyVerdict::Authorized
                } else {
                    KeyVerdict::WrongSignature
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_default_authority_knows_no_keys() {
        assert_eq!(
            NoPromotionAuthority.verify("any", "sig", b"m"),
            KeyVerdict::UnknownKey
        );
    }

    #[test]
    fn registry_verdicts_cover_the_four_outcomes() {
        let mut registry = StaticKeyRegistry::new();
        registry.authorize("ops/2026-q3");
        let message = b"baseline-content-hash";
        let good = StaticKeyRegistry::tag("ops/2026-q3", message);
        assert_eq!(
            registry.verify("ops/2026-q3", &good, message),
            KeyVerdict::Authorized
        );
        assert_eq!(
            registry.verify("ops/2026-q3", "forged", message),
            KeyVerdict::WrongSignature
        );
        assert_eq!(
            registry.verify("ops/2026-q2", &good, message),
            KeyVerdict::UnknownKey
        );
        registry.revoke("ops/2026-q3");
        assert_eq!(
            registry.verify("ops/2026-q3", &good, message),
            KeyVerdict::RevokedKey
        );
        // A tag is message-bound: any other message fails.
        registry.authorize("ops/2026-q3");
        assert_eq!(
            registry.verify("ops/2026-q3", &good, b"other message"),
            KeyVerdict::WrongSignature
        );
        // And key-bound: another key's tag never transfers.
        registry.authorize("ops/2026-q4");
        assert_eq!(
            registry.verify("ops/2026-q4", &good, message),
            KeyVerdict::WrongSignature
        );
    }
}
