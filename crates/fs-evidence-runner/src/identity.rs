//! Presented digest values and nominal Runner V2 root references.
//!
//! The wrappers here validate syntax and nominal role/domain separation only.
//! They do not establish existence, byte possession, content equivalence,
//! lifecycle completion, durability, verification, admission, or authority.

use crate::catalog::DigestRoleV2;

#[allow(
    dead_code,
    reason = "reserved for the sealed family-domain registry owned by the next schema leaf"
)]
const DOMAIN_PREFIX: &str = "org.frankensim.fs-evidence-runner.";
#[allow(
    dead_code,
    reason = "reserved for the sealed family-domain registry owned by the next schema leaf"
)]
const DOMAIN_SUFFIX: &str = ".v1";
#[allow(
    dead_code,
    reason = "reserved for the sealed family-domain registry owned by the next schema leaf"
)]
const DOMAIN_MAX_BYTES: usize = 128;
const DIGEST_BYTES: usize = 32;
const LOWER_HEX_BYTES: usize = DIGEST_BYTES * 2;

/// Deterministic construction failures for presented identity data.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum IdentityError {
    /// A domain did not follow the one frozen lowercase `.v1` generation rule.
    InvalidDomain,
    /// A digest byte slice did not contain exactly 32 bytes.
    WrongDigestLength {
        /// Observed byte length.
        observed: usize,
        /// Required byte length.
        expected: usize,
    },
    /// A textual digest did not contain exactly 64 ASCII hexadecimal bytes.
    WrongLowerHexLength {
        /// Observed byte length.
        observed: usize,
        /// Required byte length.
        expected: usize,
    },
    /// A textual digest used an uppercase or non-hexadecimal byte.
    NonCanonicalLowerHex {
        /// Zero-based byte offset.
        index: usize,
        /// Rejected byte.
        byte: u8,
    },
    /// A digest carried a role other than the nominal wrapper's exact role.
    WrongRole {
        /// Required role.
        expected: DigestRoleV2,
        /// Presented role.
        observed: DigestRoleV2,
    },
    /// A digest carried a domain other than the nominal wrapper's exact
    /// registered domain.
    WrongDomain {
        /// Required registered domain.
        expected: &'static str,
        /// Presented domain.
        observed: String,
    },
}

/// Sealed witness that a canonical domain was admitted by a crate-owned
/// registry or frozen descriptor.
///
/// Fields and registration construction are crate-private so a public caller
/// cannot turn a bare string into a registered digest domain.
///
/// ```
/// use fs_evidence_runner::identity::SourceIdentityRootV2;
///
/// let witness = SourceIdentityRootV2::DESCRIPTOR.domain_witness();
/// assert_eq!(
///     witness.as_str(),
///     "org.frankensim.fs-evidence-runner.source-identity.v1"
/// );
/// ```
///
/// ```compile_fail
/// use fs_evidence_runner::identity::DigestDomainV1;
///
/// let forged = DigestDomainV1("org.frankensim.fs-evidence-runner.forged.v1");
/// ```
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct DigestDomainV1(&'static str);

impl DigestDomainV1 {
    const fn from_frozen_descriptor(domain: &'static str) -> Self {
        Self(domain)
    }

    /// Checked registration entry point for crate-owned sealed registries.
    #[allow(
        dead_code,
        reason = "reserved for the sealed family-domain registry owned by the next schema leaf"
    )]
    pub(crate) fn from_registered(domain: &'static str) -> Result<Self, IdentityError> {
        validate_domain(domain)?;
        Ok(Self(domain))
    }

    /// Returns the exact canonical registered domain.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        self.0
    }
}

/// A role-, domain-, and width-bound presented digest value.
///
/// The all-zero byte string is valid and remains distinct from typed absence.
/// A private-field domain witness prevents a bare string from claiming sealed
/// registry membership.
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub struct DigestValueV2 {
    role: DigestRoleV2,
    domain: DigestDomainV1,
    bytes: [u8; DIGEST_BYTES],
}

impl DigestValueV2 {
    /// Validates the exact 32-byte shape against a sealed domain witness.
    pub fn new(
        role: DigestRoleV2,
        domain: DigestDomainV1,
        bytes: &[u8],
    ) -> Result<Self, IdentityError> {
        if bytes.len() != DIGEST_BYTES {
            return Err(IdentityError::WrongDigestLength {
                observed: bytes.len(),
                expected: DIGEST_BYTES,
            });
        }
        let mut exact = [0_u8; DIGEST_BYTES];
        exact.copy_from_slice(bytes);
        Ok(Self {
            role,
            domain,
            bytes: exact,
        })
    }

    /// Preserves an exact 32-byte array under a sealed domain witness.
    pub const fn from_array(
        role: DigestRoleV2,
        domain: DigestDomainV1,
        bytes: [u8; DIGEST_BYTES],
    ) -> Self {
        Self {
            role,
            domain,
            bytes,
        }
    }

    /// Parses an exact 64-byte lowercase hexadecimal presentation.
    pub fn parse_lower_hex(
        role: DigestRoleV2,
        domain: DigestDomainV1,
        lower_hex: &str,
    ) -> Result<Self, IdentityError> {
        let bytes = decode_lower_hex(lower_hex)?;
        Ok(Self::from_array(role, domain, bytes))
    }

    /// Returns the closed digest role.
    #[must_use]
    pub const fn role(&self) -> DigestRoleV2 {
        self.role
    }

    /// Returns the exact canonical domain.
    #[must_use]
    pub fn domain(&self) -> &str {
        self.domain.as_str()
    }

    /// Returns the sealed registered-domain witness.
    #[must_use]
    pub const fn domain_witness(&self) -> DigestDomainV1 {
        self.domain
    }

    /// Returns the exact 32 presented bytes.
    #[must_use]
    pub const fn bytes(&self) -> &[u8; DIGEST_BYTES] {
        &self.bytes
    }

    /// Produces the canonical lowercase 64-byte hexadecimal presentation.
    #[must_use]
    pub fn to_lower_hex(&self) -> String {
        encode_lower_hex(&self.bytes)
    }
}

#[allow(
    dead_code,
    reason = "reserved for the sealed family-domain registry owned by the next schema leaf"
)]
fn validate_domain(domain: &str) -> Result<(), IdentityError> {
    if domain.len() > DOMAIN_MAX_BYTES {
        return Err(IdentityError::InvalidDomain);
    }
    let Some(without_prefix) = domain.strip_prefix(DOMAIN_PREFIX) else {
        return Err(IdentityError::InvalidDomain);
    };
    let Some(schema) = without_prefix.strip_suffix(DOMAIN_SUFFIX) else {
        return Err(IdentityError::InvalidDomain);
    };
    if schema.is_empty() {
        return Err(IdentityError::InvalidDomain);
    }

    let mut previous_was_separator = true;
    for byte in schema.bytes() {
        if byte.is_ascii_lowercase() || byte.is_ascii_digit() {
            previous_was_separator = false;
        } else if byte == b'-' && !previous_was_separator {
            previous_was_separator = true;
        } else {
            return Err(IdentityError::InvalidDomain);
        }
    }
    if previous_was_separator {
        return Err(IdentityError::InvalidDomain);
    }
    Ok(())
}

fn decode_lower_hex(text: &str) -> Result<[u8; DIGEST_BYTES], IdentityError> {
    if text.len() != LOWER_HEX_BYTES {
        return Err(IdentityError::WrongLowerHexLength {
            observed: text.len(),
            expected: LOWER_HEX_BYTES,
        });
    }

    let input = text.as_bytes();
    let mut output = [0_u8; DIGEST_BYTES];
    for index in 0..DIGEST_BYTES {
        let high_index = index * 2;
        let low_index = high_index + 1;
        let high =
            lower_hex_nibble(input[high_index]).ok_or(IdentityError::NonCanonicalLowerHex {
                index: high_index,
                byte: input[high_index],
            })?;
        let low =
            lower_hex_nibble(input[low_index]).ok_or(IdentityError::NonCanonicalLowerHex {
                index: low_index,
                byte: input[low_index],
            })?;
        output[index] = (high << 4) | low;
    }
    Ok(output)
}

const fn lower_hex_nibble(byte: u8) -> Option<u8> {
    match byte {
        b'0'..=b'9' => Some(byte - b'0'),
        b'a'..=b'f' => Some(byte - b'a' + 10),
        _ => None,
    }
}

fn encode_lower_hex(bytes: &[u8; DIGEST_BYTES]) -> String {
    const HEX: &[u8; 16] = b"0123456789abcdef";
    let mut output = String::with_capacity(LOWER_HEX_BYTES);
    for byte in bytes {
        output.push(char::from(HEX[usize::from(byte >> 4)]));
        output.push(char::from(HEX[usize::from(byte & 0x0f)]));
    }
    output
}

/// Frozen nominal role/domain descriptor for one presented root wrapper.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct PresentedIdentityDescriptorV1 {
    schema_name: &'static str,
    domain: DigestDomainV1,
    role: DigestRoleV2,
}

impl PresentedIdentityDescriptorV1 {
    const fn new(schema_name: &'static str, domain: &'static str, role: DigestRoleV2) -> Self {
        Self {
            schema_name,
            domain: DigestDomainV1::from_frozen_descriptor(domain),
            role,
        }
    }

    /// Returns the exact kebab-case schema name.
    #[must_use]
    pub const fn schema_name(self) -> &'static str {
        self.schema_name
    }

    /// Returns the exact registered `.v1` domain.
    #[must_use]
    pub const fn domain(self) -> &'static str {
        self.domain.as_str()
    }

    /// Returns the sealed witness for this exact registered domain.
    #[must_use]
    pub const fn domain_witness(self) -> DigestDomainV1 {
        self.domain
    }

    /// Returns the exact expected digest role.
    #[must_use]
    pub const fn role(self) -> DigestRoleV2 {
        self.role
    }
}

macro_rules! define_presented_root {
    (
        $(#[$meta:meta])*
        $name:ident,
        $schema_name:literal,
        $domain:literal,
        $role:path
    ) => {
        $(#[$meta])*
        #[derive(Clone, Debug, Eq, Hash, PartialEq)]
        pub struct $name {
            digest: DigestValueV2,
        }

        impl $name {
            /// The frozen nominal schema, domain, and role descriptor.
            pub const DESCRIPTOR: PresentedIdentityDescriptorV1 =
                PresentedIdentityDescriptorV1::new($schema_name, $domain, $role);

            /// Parses a public, explicitly non-authoritative textual
            /// presentation after checking its exact nominal role and domain.
            pub fn parse_presented(
                role: DigestRoleV2,
                domain: &str,
                lower_hex: &str,
            ) -> Result<Self, IdentityError> {
                if role != Self::DESCRIPTOR.role() {
                    return Err(IdentityError::WrongRole {
                        expected: Self::DESCRIPTOR.role(),
                        observed: role,
                    });
                }
                if domain != Self::DESCRIPTOR.domain() {
                    return Err(IdentityError::WrongDomain {
                        expected: Self::DESCRIPTOR.domain(),
                        observed: domain.to_owned(),
                    });
                }
                Self::from_digest(DigestValueV2::parse_lower_hex(
                    role,
                    Self::DESCRIPTOR.domain_witness(),
                    lower_hex,
                )?)
            }

            /// Admits an exact presented digest for crate-internal schema
            /// assembly after enforcing this wrapper's nominal role/domain.
            ///
            /// This is intentionally crate-private: semantic constructors
            /// remain with their designated owning modules and phases.
            pub(crate) fn from_digest(digest: DigestValueV2) -> Result<Self, IdentityError> {
                if digest.role() != Self::DESCRIPTOR.role() {
                    return Err(IdentityError::WrongRole {
                        expected: Self::DESCRIPTOR.role(),
                        observed: digest.role(),
                    });
                }
                if digest.domain() != Self::DESCRIPTOR.domain() {
                    return Err(IdentityError::WrongDomain {
                        expected: Self::DESCRIPTOR.domain(),
                        observed: digest.domain().to_owned(),
                    });
                }
                Ok(Self { digest })
            }

            /// Returns the checked presented digest without changing role.
            #[must_use]
            pub const fn digest(&self) -> &DigestValueV2 {
                &self.digest
            }

            /// Returns the exact nominal digest role.
            #[must_use]
            pub const fn role(&self) -> DigestRoleV2 {
                self.digest.role()
            }

            /// Returns the exact 32 presented bytes.
            #[must_use]
            pub const fn bytes(&self) -> &[u8; DIGEST_BYTES] {
                self.digest.bytes()
            }

            /// Returns the exact registered nominal domain.
            #[must_use]
            pub fn domain(&self) -> &str {
                self.digest.domain()
            }
        }
    };
}

define_presented_root!(
    /// Presented source identity reference.
    SourceIdentityRootV2,
    "source-identity",
    "org.frankensim.fs-evidence-runner.source-identity.v1",
    DigestRoleV2::Source
);
define_presented_root!(
    /// Presented build identity reference.
    BuildIdentityRootV2,
    "build-identity",
    "org.frankensim.fs-evidence-runner.build-identity.v1",
    DigestRoleV2::Build
);
define_presented_root!(
    /// Presented toolchain identity reference.
    ToolchainIdentityRootV2,
    "toolchain-identity",
    "org.frankensim.fs-evidence-runner.toolchain-identity.v1",
    DigestRoleV2::Toolchain
);
define_presented_root!(
    /// Presented case-manifest reference.
    CaseManifestRootV2,
    "case-manifest",
    "org.frankensim.fs-evidence-runner.case-manifest.v1",
    DigestRoleV2::CaseManifest
);
define_presented_root!(
    /// Presented encoded-artifact reference.
    ArtifactEncodedRootV2,
    "artifact-encoded",
    "org.frankensim.fs-evidence-runner.artifact-encoded.v1",
    DigestRoleV2::ArtifactEncoded
);
define_presented_root!(
    /// Presented decoded-content artifact reference.
    ArtifactContentRootV2,
    "artifact-content",
    "org.frankensim.fs-evidence-runner.artifact-content.v1",
    DigestRoleV2::ArtifactContent
);
define_presented_root!(
    /// Presented complete stored-object reference.
    StoredObjectRootV2,
    "stored-object",
    "org.frankensim.fs-evidence-runner.stored-object.v1",
    DigestRoleV2::StoredObject
);
define_presented_root!(
    /// Presented artifact-inventory reference.
    ArtifactInventoryRootV2,
    "artifact-inventory",
    "org.frankensim.fs-evidence-runner.artifact-inventory.v1",
    DigestRoleV2::ArtifactInventory
);
define_presented_root!(
    /// Presented lifecycle-log reference.
    LifecycleLogRootV2,
    "lifecycle-log",
    "org.frankensim.fs-evidence-runner.lifecycle-log.v1",
    DigestRoleV2::LifecycleLog
);
define_presented_root!(
    /// Presented run-summary reference.
    RunSummaryRootV2,
    "run-summary",
    "org.frankensim.fs-evidence-runner.run-summary.v1",
    DigestRoleV2::RunSummary
);
define_presented_root!(
    /// Presented run-terminal-record reference.
    RunTerminalRecordRootV2,
    "run-terminal-record",
    "org.frankensim.fs-evidence-runner.run-terminal-record.v1",
    DigestRoleV2::RunTerminal
);
define_presented_root!(
    /// Presented bundle-manifest reference.
    BundleManifestRootV2,
    "bundle-manifest",
    "org.frankensim.fs-evidence-runner.bundle-manifest.v1",
    DigestRoleV2::BundleManifest
);
define_presented_root!(
    /// Presented publication-commit reference with no durability claim.
    PresentedPublicationCommitRefV2,
    "presented-publication-commit-ref",
    "org.frankensim.fs-evidence-runner.presented-publication-commit-ref.v1",
    DigestRoleV2::DurablePublication
);
define_presented_root!(
    /// Presented durable-publication identity reference.
    DurablePublicationIdentityV2,
    "durable-publication-identity",
    "org.frankensim.fs-evidence-runner.durable-publication-identity.v1",
    DigestRoleV2::DurablePublication
);
define_presented_root!(
    /// Presented seal reference with no seal-construction authority.
    SealRootV2,
    "seal",
    "org.frankensim.fs-evidence-runner.seal.v1",
    DigestRoleV2::Seal
);
define_presented_root!(
    /// Presented published-bundle-receipt reference.
    PublishedBundleReceiptRootV2,
    "published-bundle-receipt",
    "org.frankensim.fs-evidence-runner.published-bundle-receipt.v1",
    DigestRoleV2::PublishedBundleReceipt
);
define_presented_root!(
    /// Presented authority-scope reference with no authority claim.
    AuthorityScopeRootV2,
    "authority-scope",
    "org.frankensim.fs-evidence-runner.authority-scope.v1",
    DigestRoleV2::ClaimScope
);
define_presented_root!(
    /// Presented external-mutation-set reference.
    ExternalMutationSetRootV2,
    "external-mutation-set",
    "org.frankensim.fs-evidence-runner.external-mutation-set.v1",
    DigestRoleV2::Policy
);
define_presented_root!(
    /// Presented artifact-set reference.
    ArtifactSetRootV2,
    "artifact-set",
    "org.frankensim.fs-evidence-runner.artifact-set.v1",
    DigestRoleV2::ArtifactInventory
);
define_presented_root!(
    /// Presented resource-identity reference.
    ResourceIdentityRootV2,
    "resource-identity",
    "org.frankensim.fs-evidence-runner.resource-identity.v1",
    DigestRoleV2::Policy
);
define_presented_root!(
    /// Presented Runner-limits schema reference.
    RunnerLimitsSchemaRootV2,
    "runner-limits-schema",
    "org.frankensim.fs-evidence-runner.runner-limits-schema.v1",
    DigestRoleV2::Policy
);
define_presented_root!(
    /// Presented Runner-limits value reference.
    RunnerLimitsRootV2,
    "runner-limits",
    "org.frankensim.fs-evidence-runner.runner-limits.v1",
    DigestRoleV2::Policy
);
define_presented_root!(
    /// Presented Runner-budgets schema reference.
    RunnerBudgetsSchemaRootV2,
    "runner-budgets-schema",
    "org.frankensim.fs-evidence-runner.runner-budgets-schema.v1",
    DigestRoleV2::Policy
);
define_presented_root!(
    /// Presented Runner-budgets value reference.
    RunnerBudgetsRootV2,
    "runner-budgets",
    "org.frankensim.fs-evidence-runner.runner-budgets.v1",
    DigestRoleV2::Policy
);
define_presented_root!(
    /// Presented root-capability-policy reference.
    RootCapabilityPolicyRootV2,
    "root-capability-policy",
    "org.frankensim.fs-evidence-runner.root-capability-policy.v1",
    DigestRoleV2::Policy
);
define_presented_root!(
    /// Presented no-claim-scope reference.
    NoClaimScopeRootV1,
    "no-claim-scope",
    "org.frankensim.fs-evidence-runner.no-claim-scope.v1",
    DigestRoleV2::ClaimScope
);
define_presented_root!(
    /// Presented cancellation stop reference with no stop/drain claim.
    CancelledStopRootV2,
    "cancelled-stop",
    "org.frankensim.fs-evidence-runner.cancelled-stop.v1",
    DigestRoleV2::RunTerminal
);
define_presented_root!(
    /// Presented timeout stop reference with no stop/drain claim.
    TimedOutStopRootV2,
    "timed-out-stop",
    "org.frankensim.fs-evidence-runner.timed-out-stop.v1",
    DigestRoleV2::RunTerminal
);
define_presented_root!(
    /// Presented controlled-internal-error drain reference with no drain claim.
    DrainedInternalErrorRootV2,
    "drained-internal-error",
    "org.frankensim.fs-evidence-runner.drained-internal-error.v1",
    DigestRoleV2::RunTerminal
);

/// Exact, frozen presented-root descriptor inventory in declaration order.
pub const ALL_PRESENTED_IDENTITY_DESCRIPTORS_V1: [PresentedIdentityDescriptorV1; 29] = [
    SourceIdentityRootV2::DESCRIPTOR,
    BuildIdentityRootV2::DESCRIPTOR,
    ToolchainIdentityRootV2::DESCRIPTOR,
    CaseManifestRootV2::DESCRIPTOR,
    ArtifactEncodedRootV2::DESCRIPTOR,
    ArtifactContentRootV2::DESCRIPTOR,
    StoredObjectRootV2::DESCRIPTOR,
    ArtifactInventoryRootV2::DESCRIPTOR,
    LifecycleLogRootV2::DESCRIPTOR,
    RunSummaryRootV2::DESCRIPTOR,
    RunTerminalRecordRootV2::DESCRIPTOR,
    BundleManifestRootV2::DESCRIPTOR,
    PresentedPublicationCommitRefV2::DESCRIPTOR,
    DurablePublicationIdentityV2::DESCRIPTOR,
    SealRootV2::DESCRIPTOR,
    PublishedBundleReceiptRootV2::DESCRIPTOR,
    AuthorityScopeRootV2::DESCRIPTOR,
    ExternalMutationSetRootV2::DESCRIPTOR,
    ArtifactSetRootV2::DESCRIPTOR,
    ResourceIdentityRootV2::DESCRIPTOR,
    RunnerLimitsSchemaRootV2::DESCRIPTOR,
    RunnerLimitsRootV2::DESCRIPTOR,
    RunnerBudgetsSchemaRootV2::DESCRIPTOR,
    RunnerBudgetsRootV2::DESCRIPTOR,
    RootCapabilityPolicyRootV2::DESCRIPTOR,
    NoClaimScopeRootV1::DESCRIPTOR,
    CancelledStopRootV2::DESCRIPTOR,
    TimedOutStopRootV2::DESCRIPTOR,
    DrainedInternalErrorRootV2::DESCRIPTOR,
];

#[cfg(test)]
mod tests {
    use std::collections::BTreeSet;

    use super::*;

    #[test]
    fn digest_width_domain_and_all_zero_presence_are_checked() {
        let domain = SourceIdentityRootV2::DESCRIPTOR.domain_witness();
        assert_eq!(
            DigestValueV2::new(DigestRoleV2::Source, domain, &[0; 31]),
            Err(IdentityError::WrongDigestLength {
                observed: 31,
                expected: 32
            })
        );
        let zero = DigestValueV2::new(DigestRoleV2::Source, domain, &[0; 32])
            .expect("all-zero presented digest is syntactically valid");
        assert_eq!(zero.bytes(), &[0; 32]);
        assert_eq!(zero.domain(), domain.as_str());
        assert_eq!(zero.role(), DigestRoleV2::Source);

        for invalid in [
            "org.frankensim.fs-evidence-runner.Source-identity.v1",
            "org.frankensim.fs-evidence-runner.source--identity.v1",
            "org.frankensim.fs-evidence-runner.-source.v1",
            "org.frankensim.fs-evidence-runner.source-.v1",
            "org.frankensim.fs-evidence-runner.source-identity.v2",
            "other.source-identity.v1",
        ] {
            assert_eq!(
                DigestDomainV1::from_registered(invalid),
                Err(IdentityError::InvalidDomain),
                "{invalid}"
            );
        }
    }

    #[test]
    fn lowercase_hex_form_is_exact_and_round_trips_every_byte() {
        let mut bytes = [0_u8; 32];
        for (index, byte) in bytes.iter_mut().enumerate() {
            *byte = u8::try_from(index * 7).expect("bounded fixture");
        }
        let domain = BuildIdentityRootV2::DESCRIPTOR.domain_witness();
        let value = DigestValueV2::from_array(DigestRoleV2::Build, domain, bytes);
        let text = value.to_lower_hex();
        assert_eq!(text.len(), 64);
        assert!(
            text.bytes()
                .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
        );
        assert_eq!(
            DigestValueV2::parse_lower_hex(DigestRoleV2::Build, domain, &text).expect("round trip"),
            value
        );
        assert_eq!(
            DigestValueV2::parse_lower_hex(DigestRoleV2::Build, domain, &text[..63]),
            Err(IdentityError::WrongLowerHexLength {
                observed: 63,
                expected: 64
            })
        );
        let uppercase = "A0".repeat(32);
        assert_eq!(
            DigestValueV2::parse_lower_hex(DigestRoleV2::Build, domain, &uppercase),
            Err(IdentityError::NonCanonicalLowerHex {
                index: 0,
                byte: b'A'
            })
        );
    }

    #[test]
    fn descriptor_inventory_is_complete_unique_and_generation_conformant() {
        assert_eq!(ALL_PRESENTED_IDENTITY_DESCRIPTORS_V1.len(), 29);
        let expected = [
            ("source-identity", DigestRoleV2::Source),
            ("build-identity", DigestRoleV2::Build),
            ("toolchain-identity", DigestRoleV2::Toolchain),
            ("case-manifest", DigestRoleV2::CaseManifest),
            ("artifact-encoded", DigestRoleV2::ArtifactEncoded),
            ("artifact-content", DigestRoleV2::ArtifactContent),
            ("stored-object", DigestRoleV2::StoredObject),
            ("artifact-inventory", DigestRoleV2::ArtifactInventory),
            ("lifecycle-log", DigestRoleV2::LifecycleLog),
            ("run-summary", DigestRoleV2::RunSummary),
            ("run-terminal-record", DigestRoleV2::RunTerminal),
            ("bundle-manifest", DigestRoleV2::BundleManifest),
            (
                "presented-publication-commit-ref",
                DigestRoleV2::DurablePublication,
            ),
            (
                "durable-publication-identity",
                DigestRoleV2::DurablePublication,
            ),
            ("seal", DigestRoleV2::Seal),
            (
                "published-bundle-receipt",
                DigestRoleV2::PublishedBundleReceipt,
            ),
            ("authority-scope", DigestRoleV2::ClaimScope),
            ("external-mutation-set", DigestRoleV2::Policy),
            ("artifact-set", DigestRoleV2::ArtifactInventory),
            ("resource-identity", DigestRoleV2::Policy),
            ("runner-limits-schema", DigestRoleV2::Policy),
            ("runner-limits", DigestRoleV2::Policy),
            ("runner-budgets-schema", DigestRoleV2::Policy),
            ("runner-budgets", DigestRoleV2::Policy),
            ("root-capability-policy", DigestRoleV2::Policy),
            ("no-claim-scope", DigestRoleV2::ClaimScope),
            ("cancelled-stop", DigestRoleV2::RunTerminal),
            ("timed-out-stop", DigestRoleV2::RunTerminal),
            ("drained-internal-error", DigestRoleV2::RunTerminal),
        ];
        let mut names = BTreeSet::new();
        let mut domains = BTreeSet::new();
        for (descriptor, (expected_name, expected_role)) in ALL_PRESENTED_IDENTITY_DESCRIPTORS_V1
            .into_iter()
            .zip(expected)
        {
            assert_eq!(descriptor.schema_name(), expected_name);
            assert_eq!(descriptor.role(), expected_role);
            assert!(names.insert(descriptor.schema_name()));
            assert!(domains.insert(descriptor.domain()));
            assert_eq!(
                descriptor.domain(),
                format!(
                    "org.frankensim.fs-evidence-runner.{}.v1",
                    descriptor.schema_name()
                )
            );
            validate_domain(descriptor.domain()).expect("frozen descriptor domain is valid");
        }
    }

    #[test]
    fn wrappers_reject_role_and_domain_substitution() {
        let source = DigestValueV2::new(
            DigestRoleV2::Source,
            SourceIdentityRootV2::DESCRIPTOR.domain_witness(),
            &[7; 32],
        )
        .expect("valid");
        assert!(SourceIdentityRootV2::from_digest(source.clone()).is_ok());
        assert!(matches!(
            BuildIdentityRootV2::from_digest(source),
            Err(IdentityError::WrongRole {
                expected: DigestRoleV2::Build,
                observed: DigestRoleV2::Source
            })
        ));

        let wrong_domain = DigestValueV2::new(
            DigestRoleV2::Source,
            BuildIdentityRootV2::DESCRIPTOR.domain_witness(),
            &[7; 32],
        )
        .expect("syntactically valid presentation");
        assert!(matches!(
            SourceIdentityRootV2::from_digest(wrong_domain),
            Err(IdentityError::WrongDomain { .. })
        ));
    }

    macro_rules! assert_presented_wrapper {
        ($type:ty, $byte:expr) => {{
            let descriptor = <$type>::DESCRIPTOR;
            let text = format!("{:02x}", $byte).repeat(32);
            let root = <$type>::parse_presented(descriptor.role(), descriptor.domain(), &text)
                .expect("exact presented wrapper parses");
            assert_eq!(root.bytes(), &[$byte; 32]);
            assert_eq!(root.domain(), descriptor.domain());
            assert_eq!(root.digest().role(), descriptor.role());
        }};
    }

    #[test]
    fn every_nominal_wrapper_has_a_checked_presented_parser() {
        assert_presented_wrapper!(SourceIdentityRootV2, 0);
        assert_presented_wrapper!(BuildIdentityRootV2, 1);
        assert_presented_wrapper!(ToolchainIdentityRootV2, 2);
        assert_presented_wrapper!(CaseManifestRootV2, 3);
        assert_presented_wrapper!(ArtifactEncodedRootV2, 4);
        assert_presented_wrapper!(ArtifactContentRootV2, 5);
        assert_presented_wrapper!(StoredObjectRootV2, 6);
        assert_presented_wrapper!(ArtifactInventoryRootV2, 7);
        assert_presented_wrapper!(LifecycleLogRootV2, 8);
        assert_presented_wrapper!(RunSummaryRootV2, 9);
        assert_presented_wrapper!(RunTerminalRecordRootV2, 10);
        assert_presented_wrapper!(BundleManifestRootV2, 11);
        assert_presented_wrapper!(PresentedPublicationCommitRefV2, 12);
        assert_presented_wrapper!(DurablePublicationIdentityV2, 13);
        assert_presented_wrapper!(SealRootV2, 14);
        assert_presented_wrapper!(PublishedBundleReceiptRootV2, 15);
        assert_presented_wrapper!(AuthorityScopeRootV2, 16);
        assert_presented_wrapper!(ExternalMutationSetRootV2, 17);
        assert_presented_wrapper!(ArtifactSetRootV2, 18);
        assert_presented_wrapper!(ResourceIdentityRootV2, 19);
        assert_presented_wrapper!(RunnerLimitsSchemaRootV2, 20);
        assert_presented_wrapper!(RunnerLimitsRootV2, 21);
        assert_presented_wrapper!(RunnerBudgetsSchemaRootV2, 22);
        assert_presented_wrapper!(RunnerBudgetsRootV2, 23);
        assert_presented_wrapper!(RootCapabilityPolicyRootV2, 24);
        assert_presented_wrapper!(NoClaimScopeRootV1, 25);
        assert_presented_wrapper!(CancelledStopRootV2, 26);
        assert_presented_wrapper!(TimedOutStopRootV2, 27);
        assert_presented_wrapper!(DrainedInternalErrorRootV2, 28);
    }

    #[test]
    fn wrapper_parser_checks_nominal_metadata_before_text() {
        let descriptor = SourceIdentityRootV2::DESCRIPTOR;
        assert!(matches!(
            SourceIdentityRootV2::parse_presented(
                DigestRoleV2::Build,
                descriptor.domain(),
                "not-even-hex"
            ),
            Err(IdentityError::WrongRole { .. })
        ));
        assert!(matches!(
            SourceIdentityRootV2::parse_presented(
                descriptor.role(),
                BuildIdentityRootV2::DESCRIPTOR.domain(),
                "not-even-hex"
            ),
            Err(IdentityError::WrongDomain { .. })
        ));
    }

    #[test]
    fn each_digest_input_field_moves_presented_identity() {
        let source_domain = SourceIdentityRootV2::DESCRIPTOR.domain_witness();
        let base =
            DigestValueV2::new(DigestRoleV2::Source, source_domain, &[1; 32]).expect("valid");
        let role_changed =
            DigestValueV2::new(DigestRoleV2::Build, source_domain, &[1; 32]).expect("valid");
        let domain_changed = DigestValueV2::new(
            DigestRoleV2::Source,
            BuildIdentityRootV2::DESCRIPTOR.domain_witness(),
            &[1; 32],
        )
        .expect("valid");
        let bytes_changed =
            DigestValueV2::new(DigestRoleV2::Source, source_domain, &[2; 32]).expect("valid");
        assert_ne!(base, role_changed);
        assert_ne!(base, domain_changed);
        assert_ne!(base, bytes_changed);
    }
}
