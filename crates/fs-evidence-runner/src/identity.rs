//! Presented digest values and nominal Runner V2 root references.
//!
//! The wrappers here validate syntax and nominal role/domain separation only.
//! They do not establish existence, byte possession, content equivalence,
//! lifecycle completion, durability, verification, admission, or authority.

use crate::canonical::CanonicalFrameV1;
use crate::catalog::DigestRoleV2;
use crate::construction::{
    ConstructionClosedSemanticV2, ConstructionErrorKindV2, ConstructionErrorV2,
    ConstructionFixedObservationV2, ConstructionObservedDataClassV2, ConstructionObservedV2,
};
use fs_blake3::ContentHash;

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
/// ```compile_fail,E0423
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
///
/// The generic digest constructor preserves its exact role, registered domain,
/// and bytes:
///
/// ```
/// use fs_evidence_runner::catalog::DigestRoleV2;
/// use fs_evidence_runner::identity::{DigestValueV2, SourceIdentityRootV2};
///
/// let digest = DigestValueV2::from_array(
///     DigestRoleV2::Source,
///     SourceIdentityRootV2::DESCRIPTOR.domain_witness(),
///     [0_u8; 32],
/// );
///
/// assert_eq!(digest.role(), DigestRoleV2::Source);
/// assert_eq!(digest.domain(), SourceIdentityRootV2::DESCRIPTOR.domain());
/// assert_eq!(digest.bytes(), &[0_u8; 32]);
/// ```
///
/// Even a role- and domain-matching generic digest cannot coerce into a
/// nominal wrapper:
///
/// ```compile_fail,E0308
/// use fs_evidence_runner::catalog::DigestRoleV2;
/// use fs_evidence_runner::identity::{DigestValueV2, SourceIdentityRootV2};
///
/// let digest = DigestValueV2::from_array(
///     DigestRoleV2::Source,
///     SourceIdentityRootV2::DESCRIPTOR.domain_witness(),
///     [0_u8; 32],
/// );
/// let _source: SourceIdentityRootV2 = digest;
/// ```
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
    #[must_use]
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
    for (index, output_byte) in output.iter_mut().enumerate() {
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
        *output_byte = (high << 4) | low;
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

impl ConstructionClosedSemanticV2 for PresentedIdentityDescriptorV1 {
    fn construction_stable_name(&self) -> &'static str {
        self.schema_name()
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
    ///
    /// Public construction is explicitly presented and non-authoritative:
    ///
    /// ```
    /// use fs_evidence_runner::catalog::DigestRoleV2;
    /// use fs_evidence_runner::identity::SourceIdentityRootV2;
    ///
    /// let source = SourceIdentityRootV2::parse_presented(
    ///     DigestRoleV2::Source,
    ///     SourceIdentityRootV2::DESCRIPTOR.domain(),
    ///     &"11".repeat(32),
    /// )
    /// .unwrap();
    /// assert_eq!(source.role(), DigestRoleV2::Source);
    /// ```
    ///
    /// The generic-digest semantic constructor is not public:
    ///
    /// ```compile_fail,E0624
    /// use fs_evidence_runner::catalog::DigestRoleV2;
    /// use fs_evidence_runner::identity::{DigestValueV2, SourceIdentityRootV2};
    ///
    /// let digest = DigestValueV2::from_array(
    ///     DigestRoleV2::Source,
    ///     SourceIdentityRootV2::DESCRIPTOR.domain_witness(),
    ///     [0_u8; 32],
    /// );
    /// let _source = SourceIdentityRootV2::from_digest(digest);
    /// ```
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
    ///
    /// Presented syntax can be checked without claiming lifecycle completion:
    ///
    /// ```
    /// use fs_evidence_runner::catalog::DigestRoleV2;
    /// use fs_evidence_runner::identity::LifecycleLogRootV2;
    ///
    /// let root = LifecycleLogRootV2::parse_presented(
    ///     DigestRoleV2::LifecycleLog,
    ///     LifecycleLogRootV2::DESCRIPTOR.domain(),
    ///     &"22".repeat(32),
    /// )
    /// .unwrap();
    /// assert_eq!(root.bytes(), &[0x22_u8; 32]);
    /// ```
    ///
    /// Lifecycle construction remains with the lifecycle owner:
    ///
    /// ```compile_fail,E0624
    /// use fs_evidence_runner::catalog::DigestRoleV2;
    /// use fs_evidence_runner::identity::{DigestValueV2, LifecycleLogRootV2};
    ///
    /// let digest = DigestValueV2::from_array(
    ///     DigestRoleV2::LifecycleLog,
    ///     LifecycleLogRootV2::DESCRIPTOR.domain_witness(),
    ///     [0x22_u8; 32],
    /// );
    /// let _root = LifecycleLogRootV2::from_digest(digest);
    /// ```
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
    ///
    /// The public parser checks only the presented nominal value:
    ///
    /// ```
    /// use fs_evidence_runner::catalog::DigestRoleV2;
    /// use fs_evidence_runner::identity::DurablePublicationIdentityV2;
    ///
    /// let root = DurablePublicationIdentityV2::parse_presented(
    ///     DigestRoleV2::DurablePublication,
    ///     DurablePublicationIdentityV2::DESCRIPTOR.domain(),
    ///     &"33".repeat(32),
    /// )
    /// .unwrap();
    /// assert_eq!(root.role(), DigestRoleV2::DurablePublication);
    /// ```
    ///
    /// Presented bytes cannot invoke the private durability constructor:
    ///
    /// ```compile_fail,E0624
    /// use fs_evidence_runner::catalog::DigestRoleV2;
    /// use fs_evidence_runner::identity::{
    ///     DigestValueV2, DurablePublicationIdentityV2,
    /// };
    ///
    /// let digest = DigestValueV2::from_array(
    ///     DigestRoleV2::DurablePublication,
    ///     DurablePublicationIdentityV2::DESCRIPTOR.domain_witness(),
    ///     [0x33_u8; 32],
    /// );
    /// let _root = DurablePublicationIdentityV2::from_digest(digest);
    /// ```
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
    ///
    /// Parsing exposes the nominal bytes but grants no authority:
    ///
    /// ```
    /// use fs_evidence_runner::catalog::DigestRoleV2;
    /// use fs_evidence_runner::identity::AuthorityScopeRootV2;
    ///
    /// let root = AuthorityScopeRootV2::parse_presented(
    ///     DigestRoleV2::ClaimScope,
    ///     AuthorityScopeRootV2::DESCRIPTOR.domain(),
    ///     &"44".repeat(32),
    /// )
    /// .unwrap();
    /// assert_eq!(root.domain(), AuthorityScopeRootV2::DESCRIPTOR.domain());
    /// ```
    ///
    /// Authority construction remains private to its owning phase:
    ///
    /// ```compile_fail,E0624
    /// use fs_evidence_runner::catalog::DigestRoleV2;
    /// use fs_evidence_runner::identity::{AuthorityScopeRootV2, DigestValueV2};
    ///
    /// let digest = DigestValueV2::from_array(
    ///     DigestRoleV2::ClaimScope,
    ///     AuthorityScopeRootV2::DESCRIPTOR.domain_witness(),
    ///     [0x44_u8; 32],
    /// );
    /// let _root = AuthorityScopeRootV2::from_digest(digest);
    /// ```
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
define_presented_root!(
    /// Presented comparison-expression reference with no evaluator claim.
    ComparisonExprRootV2,
    "comparison-expression",
    "org.frankensim.fs-evidence-runner.comparison-expression.v1",
    DigestRoleV2::Spec
);
define_presented_root!(
    /// Presented effect-expectation reference with no evaluator claim.
    EffectExpectationRootV2,
    "effect-expectation",
    "org.frankensim.fs-evidence-runner.effect-expectation.v1",
    DigestRoleV2::Spec
);
define_presented_root!(
    /// Presented observation-batch reference with no execution claim.
    PresentedObservationBatchRootV2,
    "presented-observation-batch",
    "org.frankensim.fs-evidence-runner.presented-observation-batch.v1",
    DigestRoleV2::Run
);
define_presented_root!(
    /// Presented observation-set reference with no execution claim.
    ObservationSetRootV2,
    "observation-set",
    "org.frankensim.fs-evidence-runner.observation-set.v1",
    DigestRoleV2::Run
);
define_presented_root!(
    /// Presented effect-snapshot input with no execution claim.
    PresentedEffectSnapshotRootV2,
    "presented-effect-snapshot",
    "org.frankensim.fs-evidence-runner.presented-effect-snapshot.v1",
    DigestRoleV2::Run
);
define_presented_root!(
    /// Presented effect-snapshot requirements with no evaluator claim.
    EffectSnapshotRequirementsRootV2,
    "effect-snapshot-requirements",
    "org.frankensim.fs-evidence-runner.effect-snapshot-requirements.v1",
    DigestRoleV2::Spec
);
define_presented_root!(
    /// Presented effect-snapshot result with no execution claim.
    EffectSnapshotRootV2,
    "effect-snapshot",
    "org.frankensim.fs-evidence-runner.effect-snapshot.v1",
    DigestRoleV2::Run
);
define_presented_root!(
    /// Presented evaluator registry view with no policy-admission claim.
    EvaluationRegistryViewRootV2,
    "evaluation-registry-view",
    "org.frankensim.fs-evidence-runner.evaluation-registry-view.v1",
    DigestRoleV2::Policy
);
define_presented_root!(
    /// Presented evaluation context with no execution claim.
    EvaluationContextRootV2,
    "evaluation-context",
    "org.frankensim.fs-evidence-runner.evaluation-context.v1",
    DigestRoleV2::Run
);
define_presented_root!(
    /// Presented evaluation error with no conformance claim.
    EvaluationErrorRootV2,
    "evaluation-error",
    "org.frankensim.fs-evidence-runner.evaluation-error.v1",
    DigestRoleV2::Run
);
define_presented_root!(
    /// Presented mismatch detail with no conformance claim.
    MismatchDetailRootV2,
    "mismatch-detail",
    "org.frankensim.fs-evidence-runner.mismatch-detail.v1",
    DigestRoleV2::Run
);
define_presented_root!(
    /// Presented case-conformance specification with no execution claim.
    CaseConformanceSpecRootV2,
    "case-conformance-spec",
    "org.frankensim.fs-evidence-runner.case-conformance-spec.v1",
    DigestRoleV2::Spec
);
define_presented_root!(
    /// Presented raw case outcome with no conformance claim.
    CaseRawOutcomeRootV2,
    "case-raw-outcome",
    "org.frankensim.fs-evidence-runner.case-raw-outcome.v1",
    DigestRoleV2::Run
);
define_presented_root!(
    /// Presented case-conformance verdict with no authority claim.
    CaseConformanceVerdictRootV2,
    "case-conformance-verdict",
    "org.frankensim.fs-evidence-runner.case-conformance-verdict.v1",
    DigestRoleV2::Run
);

/// Exact, frozen presented-root descriptor inventory in declaration order.
pub const ALL_PRESENTED_IDENTITY_DESCRIPTORS_V1: [PresentedIdentityDescriptorV1; 43] = [
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
    ComparisonExprRootV2::DESCRIPTOR,
    EffectExpectationRootV2::DESCRIPTOR,
    PresentedObservationBatchRootV2::DESCRIPTOR,
    ObservationSetRootV2::DESCRIPTOR,
    PresentedEffectSnapshotRootV2::DESCRIPTOR,
    EffectSnapshotRequirementsRootV2::DESCRIPTOR,
    EffectSnapshotRootV2::DESCRIPTOR,
    EvaluationRegistryViewRootV2::DESCRIPTOR,
    EvaluationContextRootV2::DESCRIPTOR,
    EvaluationErrorRootV2::DESCRIPTOR,
    MismatchDetailRootV2::DESCRIPTOR,
    CaseConformanceSpecRootV2::DESCRIPTOR,
    CaseRawOutcomeRootV2::DESCRIPTOR,
    CaseConformanceVerdictRootV2::DESCRIPTOR,
];

/// Canonical identity domain for the non-wire constructor-owner handoff table.
pub const CONSTRUCTOR_OWNER_HANDOFF_PROJECTION_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.constructor-owner-handoff-projection.v1";

/// Sole semantic-constructor owner for a nominal presented-root schema.
///
/// These values are local planning and handoff data, not wire discriminants.
/// They neither construct a semantic identity nor prove that a downstream
/// owner has completed its continuation.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub enum ConstructorOwnerV1 {
    /// This base-schema leaf owns the complete semantic constructor.
    BaseSchema,
    /// Invocation/run identity leaf owns the complete semantic constructor.
    InvocationAndRunIdentity,
    /// RunnerSpec/family registry owns the complete semantic constructor.
    RunnerSpecAndFamilyRegistry,
    /// Artifact schema leaf owns the complete semantic constructor.
    ArtifactSchema,
    /// Lifecycle phase owns the complete semantic constructor.
    Lifecycle,
    /// Physical capability phase owns the complete semantic constructor.
    CapabilityAcquisition,
    /// Durable-publication phase owns the complete semantic constructor.
    DurablePublication,
    /// Authority-coherence phase owns the complete semantic constructor.
    AuthorityCoherence,
    /// Closed comparison/effect/case-conformance phase owns the constructor.
    EvaluatorAndCaseConformance,
    /// Canonical byte/codec phase owns byte-derived artifact constructors.
    CanonicalBytesAndCodec,
}

impl ConstructorOwnerV1 {
    /// Every owner represented by the frozen handoff table.
    pub const ALL: [Self; 10] = [
        Self::BaseSchema,
        Self::InvocationAndRunIdentity,
        Self::RunnerSpecAndFamilyRegistry,
        Self::ArtifactSchema,
        Self::Lifecycle,
        Self::CapabilityAcquisition,
        Self::DurablePublication,
        Self::AuthorityCoherence,
        Self::EvaluatorAndCaseConformance,
        Self::CanonicalBytesAndCodec,
    ];

    /// Stable non-wire owner code.
    #[must_use]
    pub const fn code(self) -> u16 {
        match self {
            Self::BaseSchema => 1,
            Self::InvocationAndRunIdentity => 2,
            Self::RunnerSpecAndFamilyRegistry => 3,
            Self::ArtifactSchema => 4,
            Self::Lifecycle => 5,
            Self::CapabilityAcquisition => 6,
            Self::DurablePublication => 7,
            Self::AuthorityCoherence => 8,
            Self::EvaluatorAndCaseConformance => 9,
            Self::CanonicalBytesAndCodec => 10,
        }
    }

    /// Resolve one exact constructor-owner code.
    ///
    /// # Errors
    ///
    /// Unknown and reserved codes refuse rather than selecting a default
    /// owner.
    pub fn from_code(code: u16) -> Result<Self, ConstructionErrorV2> {
        Self::ALL
            .into_iter()
            .find(|owner| owner.code() == code)
            .ok_or_else(|| {
                ConstructionErrorV2::new(
                    ConstructionErrorKindV2::UnknownCode,
                    "constructor_owner_handoff.presented.owner_code",
                    "one exact registered constructor-owner code",
                    code,
                )
            })
    }

    /// Stable lowercase owner name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::BaseSchema => "base-schema",
            Self::InvocationAndRunIdentity => "invocation-and-run-identity",
            Self::RunnerSpecAndFamilyRegistry => "runner-spec-and-family-registry",
            Self::ArtifactSchema => "artifact-schema",
            Self::Lifecycle => "lifecycle",
            Self::CapabilityAcquisition => "capability-acquisition",
            Self::DurablePublication => "durable-publication",
            Self::AuthorityCoherence => "authority-coherence",
            Self::EvaluatorAndCaseConformance => "evaluator-and-case-conformance",
            Self::CanonicalBytesAndCodec => "canonical-bytes-and-codec",
        }
    }

    /// Exact owning Bead or phase identifier.
    #[must_use]
    pub const fn owner_id(self) -> &'static str {
        match self {
            Self::BaseSchema => "frankensim-epic-foundations-huq.24.1.1.1",
            Self::InvocationAndRunIdentity => "frankensim-epic-foundations-huq.24.1.1.2",
            Self::RunnerSpecAndFamilyRegistry => "frankensim-epic-foundations-huq.24.1.1.3",
            Self::ArtifactSchema => "frankensim-epic-foundations-huq.24.1.1.4",
            Self::Lifecycle => "frankensim-epic-foundations-huq.24.1.2",
            Self::CapabilityAcquisition => "frankensim-epic-foundations-huq.24.2.2.1",
            Self::DurablePublication => "frankensim-epic-foundations-huq.24.2.2",
            Self::AuthorityCoherence => "frankensim-epic-foundations-huq.24.5",
            Self::EvaluatorAndCaseConformance => "frankensim-epic-foundations-huq.24.1.1.3.1",
            Self::CanonicalBytesAndCodec => "frankensim-epic-foundations-huq.24.2.1.2",
        }
    }
}

impl ConstructionClosedSemanticV2 for ConstructorOwnerV1 {
    fn construction_stable_name(&self) -> &'static str {
        self.name()
    }
}

/// One source-defined nominal-wrapper to semantic-constructor-owner handoff.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct ConstructorOwnerHandoffEntryV1 {
    descriptor: PresentedIdentityDescriptorV1,
    owner: ConstructorOwnerV1,
}

impl ConstructorOwnerHandoffEntryV1 {
    const fn new(descriptor: PresentedIdentityDescriptorV1, owner: ConstructorOwnerV1) -> Self {
        Self { descriptor, owner }
    }

    /// Exact nominal wrapper/domain/role descriptor.
    #[must_use]
    pub const fn descriptor(self) -> PresentedIdentityDescriptorV1 {
        self.descriptor
    }

    /// Current or sole downstream semantic-constructor owner.
    #[must_use]
    pub const fn owner(self) -> ConstructorOwnerV1 {
        self.owner
    }
}

/// One untrusted, non-wire constructor-owner handoff row.
///
/// Every nominal descriptor and owner witness is retained independently so
/// exact reconstruction can reject stale names and Bead identifiers even when
/// a numeric code still resolves. Construction of this presentation does not
/// admit it or establish any downstream ownership.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PresentedConstructorOwnerHandoffEntryV1 {
    schema_name: String,
    domain: String,
    role_code: u16,
    owner_code: u16,
    owner_name: String,
    owner_id: String,
}

impl PresentedConstructorOwnerHandoffEntryV1 {
    /// Construct one explicitly presented row without admitting it.
    pub fn new(
        schema_name: impl Into<String>,
        domain: impl Into<String>,
        role_code: u16,
        owner_code: u16,
        owner_name: impl Into<String>,
        owner_id: impl Into<String>,
    ) -> Self {
        Self {
            schema_name: schema_name.into(),
            domain: domain.into(),
            role_code,
            owner_code,
            owner_name: owner_name.into(),
            owner_id: owner_id.into(),
        }
    }

    /// Exact presented nominal schema name.
    #[must_use]
    pub fn schema_name(&self) -> &str {
        &self.schema_name
    }

    /// Exact presented nominal domain.
    #[must_use]
    pub fn domain(&self) -> &str {
        &self.domain
    }

    /// Raw presented digest-role code.
    #[must_use]
    pub const fn role_code(&self) -> u16 {
        self.role_code
    }

    /// Raw presented constructor-owner code.
    #[must_use]
    pub const fn owner_code(&self) -> u16 {
        self.owner_code
    }

    /// Presented stable owner name.
    #[must_use]
    pub fn owner_name(&self) -> &str {
        &self.owner_name
    }

    /// Presented owning Bead or phase identifier.
    #[must_use]
    pub fn owner_id(&self) -> &str {
        &self.owner_id
    }
}

impl From<ConstructorOwnerHandoffEntryV1> for PresentedConstructorOwnerHandoffEntryV1 {
    fn from(entry: ConstructorOwnerHandoffEntryV1) -> Self {
        Self::new(
            entry.descriptor().schema_name(),
            entry.descriptor().domain(),
            entry.descriptor().role().code(),
            entry.owner().code(),
            entry.owner().name(),
            entry.owner().owner_id(),
        )
    }
}

/// Exact source-defined constructor-owner handoff inventory.
pub const FROZEN_CONSTRUCTOR_OWNER_HANDOFF_ENTRIES_V1: [ConstructorOwnerHandoffEntryV1; 43] = [
    ConstructorOwnerHandoffEntryV1::new(
        SourceIdentityRootV2::DESCRIPTOR,
        ConstructorOwnerV1::InvocationAndRunIdentity,
    ),
    ConstructorOwnerHandoffEntryV1::new(
        BuildIdentityRootV2::DESCRIPTOR,
        ConstructorOwnerV1::InvocationAndRunIdentity,
    ),
    ConstructorOwnerHandoffEntryV1::new(
        ToolchainIdentityRootV2::DESCRIPTOR,
        ConstructorOwnerV1::InvocationAndRunIdentity,
    ),
    ConstructorOwnerHandoffEntryV1::new(
        CaseManifestRootV2::DESCRIPTOR,
        ConstructorOwnerV1::RunnerSpecAndFamilyRegistry,
    ),
    ConstructorOwnerHandoffEntryV1::new(
        ArtifactEncodedRootV2::DESCRIPTOR,
        ConstructorOwnerV1::CanonicalBytesAndCodec,
    ),
    ConstructorOwnerHandoffEntryV1::new(
        ArtifactContentRootV2::DESCRIPTOR,
        ConstructorOwnerV1::CanonicalBytesAndCodec,
    ),
    ConstructorOwnerHandoffEntryV1::new(
        StoredObjectRootV2::DESCRIPTOR,
        ConstructorOwnerV1::CanonicalBytesAndCodec,
    ),
    ConstructorOwnerHandoffEntryV1::new(
        ArtifactInventoryRootV2::DESCRIPTOR,
        ConstructorOwnerV1::ArtifactSchema,
    ),
    ConstructorOwnerHandoffEntryV1::new(
        LifecycleLogRootV2::DESCRIPTOR,
        ConstructorOwnerV1::Lifecycle,
    ),
    ConstructorOwnerHandoffEntryV1::new(
        RunSummaryRootV2::DESCRIPTOR,
        ConstructorOwnerV1::Lifecycle,
    ),
    ConstructorOwnerHandoffEntryV1::new(
        RunTerminalRecordRootV2::DESCRIPTOR,
        ConstructorOwnerV1::Lifecycle,
    ),
    ConstructorOwnerHandoffEntryV1::new(
        BundleManifestRootV2::DESCRIPTOR,
        ConstructorOwnerV1::ArtifactSchema,
    ),
    ConstructorOwnerHandoffEntryV1::new(
        PresentedPublicationCommitRefV2::DESCRIPTOR,
        ConstructorOwnerV1::DurablePublication,
    ),
    ConstructorOwnerHandoffEntryV1::new(
        DurablePublicationIdentityV2::DESCRIPTOR,
        ConstructorOwnerV1::DurablePublication,
    ),
    ConstructorOwnerHandoffEntryV1::new(
        SealRootV2::DESCRIPTOR,
        ConstructorOwnerV1::DurablePublication,
    ),
    ConstructorOwnerHandoffEntryV1::new(
        PublishedBundleReceiptRootV2::DESCRIPTOR,
        ConstructorOwnerV1::DurablePublication,
    ),
    ConstructorOwnerHandoffEntryV1::new(
        AuthorityScopeRootV2::DESCRIPTOR,
        ConstructorOwnerV1::AuthorityCoherence,
    ),
    ConstructorOwnerHandoffEntryV1::new(
        ExternalMutationSetRootV2::DESCRIPTOR,
        ConstructorOwnerV1::EvaluatorAndCaseConformance,
    ),
    ConstructorOwnerHandoffEntryV1::new(
        ArtifactSetRootV2::DESCRIPTOR,
        ConstructorOwnerV1::EvaluatorAndCaseConformance,
    ),
    ConstructorOwnerHandoffEntryV1::new(
        ResourceIdentityRootV2::DESCRIPTOR,
        ConstructorOwnerV1::CapabilityAcquisition,
    ),
    ConstructorOwnerHandoffEntryV1::new(
        RunnerLimitsSchemaRootV2::DESCRIPTOR,
        ConstructorOwnerV1::BaseSchema,
    ),
    ConstructorOwnerHandoffEntryV1::new(
        RunnerLimitsRootV2::DESCRIPTOR,
        ConstructorOwnerV1::BaseSchema,
    ),
    ConstructorOwnerHandoffEntryV1::new(
        RunnerBudgetsSchemaRootV2::DESCRIPTOR,
        ConstructorOwnerV1::BaseSchema,
    ),
    ConstructorOwnerHandoffEntryV1::new(
        RunnerBudgetsRootV2::DESCRIPTOR,
        ConstructorOwnerV1::BaseSchema,
    ),
    ConstructorOwnerHandoffEntryV1::new(
        RootCapabilityPolicyRootV2::DESCRIPTOR,
        ConstructorOwnerV1::BaseSchema,
    ),
    ConstructorOwnerHandoffEntryV1::new(
        NoClaimScopeRootV1::DESCRIPTOR,
        ConstructorOwnerV1::BaseSchema,
    ),
    ConstructorOwnerHandoffEntryV1::new(
        CancelledStopRootV2::DESCRIPTOR,
        ConstructorOwnerV1::Lifecycle,
    ),
    ConstructorOwnerHandoffEntryV1::new(
        TimedOutStopRootV2::DESCRIPTOR,
        ConstructorOwnerV1::Lifecycle,
    ),
    ConstructorOwnerHandoffEntryV1::new(
        DrainedInternalErrorRootV2::DESCRIPTOR,
        ConstructorOwnerV1::Lifecycle,
    ),
    ConstructorOwnerHandoffEntryV1::new(
        ComparisonExprRootV2::DESCRIPTOR,
        ConstructorOwnerV1::EvaluatorAndCaseConformance,
    ),
    ConstructorOwnerHandoffEntryV1::new(
        EffectExpectationRootV2::DESCRIPTOR,
        ConstructorOwnerV1::EvaluatorAndCaseConformance,
    ),
    ConstructorOwnerHandoffEntryV1::new(
        PresentedObservationBatchRootV2::DESCRIPTOR,
        ConstructorOwnerV1::EvaluatorAndCaseConformance,
    ),
    ConstructorOwnerHandoffEntryV1::new(
        ObservationSetRootV2::DESCRIPTOR,
        ConstructorOwnerV1::EvaluatorAndCaseConformance,
    ),
    ConstructorOwnerHandoffEntryV1::new(
        PresentedEffectSnapshotRootV2::DESCRIPTOR,
        ConstructorOwnerV1::EvaluatorAndCaseConformance,
    ),
    ConstructorOwnerHandoffEntryV1::new(
        EffectSnapshotRequirementsRootV2::DESCRIPTOR,
        ConstructorOwnerV1::EvaluatorAndCaseConformance,
    ),
    ConstructorOwnerHandoffEntryV1::new(
        EffectSnapshotRootV2::DESCRIPTOR,
        ConstructorOwnerV1::EvaluatorAndCaseConformance,
    ),
    ConstructorOwnerHandoffEntryV1::new(
        EvaluationRegistryViewRootV2::DESCRIPTOR,
        ConstructorOwnerV1::EvaluatorAndCaseConformance,
    ),
    ConstructorOwnerHandoffEntryV1::new(
        EvaluationContextRootV2::DESCRIPTOR,
        ConstructorOwnerV1::EvaluatorAndCaseConformance,
    ),
    ConstructorOwnerHandoffEntryV1::new(
        EvaluationErrorRootV2::DESCRIPTOR,
        ConstructorOwnerV1::EvaluatorAndCaseConformance,
    ),
    ConstructorOwnerHandoffEntryV1::new(
        MismatchDetailRootV2::DESCRIPTOR,
        ConstructorOwnerV1::EvaluatorAndCaseConformance,
    ),
    ConstructorOwnerHandoffEntryV1::new(
        CaseConformanceSpecRootV2::DESCRIPTOR,
        ConstructorOwnerV1::EvaluatorAndCaseConformance,
    ),
    ConstructorOwnerHandoffEntryV1::new(
        CaseRawOutcomeRootV2::DESCRIPTOR,
        ConstructorOwnerV1::EvaluatorAndCaseConformance,
    ),
    ConstructorOwnerHandoffEntryV1::new(
        CaseConformanceVerdictRootV2::DESCRIPTOR,
        ConstructorOwnerV1::EvaluatorAndCaseConformance,
    ),
];

/// Bounded, source-defined, non-wire nominal-constructor handoff projection.
///
/// This projection exact-set binds every nominal wrapper descriptor to its
/// current or sole downstream semantic-constructor owner. It does not expose
/// any semantic, lifecycle, durability, stop/drain, seal, or authority
/// constructor.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ConstructorOwnerHandoffProjectionV1 {
    entries: Box<[ConstructorOwnerHandoffEntryV1]>,
    root: ContentHash,
}

impl ConstructorOwnerHandoffProjectionV1 {
    /// Reconstruct the frozen mapping as an exact set.
    ///
    /// Reordering is accepted and canonicalized. Missing, extra, duplicate,
    /// descriptor-mutated, or owner-mutated rows refuse.
    pub fn reconstruct_exact_set(
        entries: &[ConstructorOwnerHandoffEntryV1],
    ) -> Result<Self, ConstructionErrorV2> {
        if entries.len() != FROZEN_CONSTRUCTOR_OWNER_HANDOFF_ENTRIES_V1.len() {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::OutOfRange,
                "constructor_owner_handoff.entries",
                "exactly one row for every nominal identity descriptor",
                entries.len(),
            ));
        }

        let mut schema_names = std::collections::BTreeSet::new();
        for entry in entries {
            if !schema_names.insert(entry.descriptor.schema_name()) {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Duplicate,
                    "constructor_owner_handoff.descriptor",
                    "one row per nominal schema name",
                    ConstructionObservedV2::closed(&entry.descriptor),
                ));
            }
        }

        let mut canonical = Vec::with_capacity(entries.len());
        for expected in FROZEN_CONSTRUCTOR_OWNER_HANDOFF_ENTRIES_V1 {
            let Some(observed) = entries
                .iter()
                .find(|entry| entry.descriptor == expected.descriptor)
            else {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Missing,
                    "constructor_owner_handoff.descriptor",
                    "the exact frozen schema/domain/role descriptor",
                    ConstructionObservedV2::closed(&expected.descriptor),
                ));
            };
            if observed.owner != expected.owner {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Incompatible,
                    "constructor_owner_handoff.owner",
                    "the exact current or sole downstream constructor owner",
                    ConstructionObservedV2::closed(&observed.owner),
                ));
            }
            canonical.push(expected);
        }
        let root = constructor_owner_handoff_root(&canonical)?;
        Ok(Self {
            entries: canonical.into_boxed_slice(),
            root,
        })
    }

    /// Reconstruct the exact frozen closeout sequence.
    ///
    /// This first preserves every exact-set validation and canonicalization
    /// rule, then additionally refuses a caller-presented permutation. The
    /// sequence form is used by source-closure, E2E, coverage, and log proof;
    /// callers that semantically own an unordered set retain
    /// [`Self::reconstruct_exact_set`].
    pub fn reconstruct_exact_sequence(
        entries: &[ConstructorOwnerHandoffEntryV1],
    ) -> Result<Self, ConstructionErrorV2> {
        let projection = Self::reconstruct_exact_set(entries)?;
        if entries != FROZEN_CONSTRUCTOR_OWNER_HANDOFF_ENTRIES_V1 {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::OutOfOrder,
                "constructor_owner_handoff.entries",
                "the exact frozen 43-row closeout sequence",
                ConstructionObservedV2::fixed(
                    ConstructionFixedObservationV2::ExactSetDifferentOrder,
                ),
            ));
        }
        Ok(projection)
    }

    /// Materialize the exact frozen handoff as independently presented rows.
    ///
    /// The result is suitable for serialization-facing or E2E fixtures, but
    /// remains non-authoritative until passed back through
    /// [`Self::reconstruct_presented_exact_sequence`].
    #[must_use]
    pub fn frozen_presented_sequence() -> Box<[PresentedConstructorOwnerHandoffEntryV1]> {
        FROZEN_CONSTRUCTOR_OWNER_HANDOFF_ENTRIES_V1
            .into_iter()
            .map(PresentedConstructorOwnerHandoffEntryV1::from)
            .collect()
    }

    /// Reconstruct the exact ordered handoff from raw presented fields.
    ///
    /// Unknown codes, stale owner names or Bead identifiers, descriptor
    /// mutation, duplication, omission, insertion, and reordering all refuse
    /// with a field-specific diagnostic. No numeric owner code is sufficient
    /// without its matching stable name and owner identifier.
    pub fn reconstruct_presented_exact_sequence(
        entries: &[PresentedConstructorOwnerHandoffEntryV1],
    ) -> Result<Self, ConstructionErrorV2> {
        if entries.len() != FROZEN_CONSTRUCTOR_OWNER_HANDOFF_ENTRIES_V1.len() {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::OutOfRange,
                "constructor_owner_handoff.presented.entries",
                "the exact frozen 43-row presented sequence",
                entries.len(),
            ));
        }

        let mut schema_names = std::collections::BTreeSet::new();
        for entry in entries {
            if !schema_names.insert(entry.schema_name()) {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Duplicate,
                    "constructor_owner_handoff.presented.schema_name",
                    "one row per exact nominal schema name",
                    ConstructionObservedV2::fixed(
                        ConstructionFixedObservationV2::RepeatedPresentedSchemaName,
                    ),
                ));
            }
        }

        for (observed, expected) in entries
            .iter()
            .zip(FROZEN_CONSTRUCTOR_OWNER_HANDOFF_ENTRIES_V1)
        {
            let expected_descriptor = expected.descriptor();
            if observed.schema_name() != expected_descriptor.schema_name() {
                let is_known_out_of_order = FROZEN_CONSTRUCTOR_OWNER_HANDOFF_ENTRIES_V1
                    .iter()
                    .any(|entry| entry.descriptor().schema_name() == observed.schema_name());
                return Err(ConstructionErrorV2::new(
                    if is_known_out_of_order {
                        ConstructionErrorKindV2::OutOfOrder
                    } else {
                        ConstructionErrorKindV2::UnknownCode
                    },
                    "constructor_owner_handoff.presented.schema_name",
                    expected_descriptor.schema_name(),
                    ConstructionObservedV2::fixed(
                        ConstructionFixedObservationV2::DifferentOrOutOfOrderSchemaName,
                    ),
                ));
            }
            if observed.domain() != expected_descriptor.domain() {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Incompatible,
                    "constructor_owner_handoff.presented.domain",
                    expected_descriptor.domain(),
                    ConstructionObservedV2::fixed(
                        ConstructionFixedObservationV2::DifferentNominalDomain,
                    ),
                ));
            }
            let observed_role = DigestRoleV2::from_code(observed.role_code()).map_err(|_| {
                ConstructionErrorV2::new(
                    ConstructionErrorKindV2::UnknownCode,
                    "constructor_owner_handoff.presented.role_code",
                    "one exact closed digest-role code",
                    observed.role_code(),
                )
            })?;
            if observed_role != expected_descriptor.role() {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Incompatible,
                    "constructor_owner_handoff.presented.role_code",
                    "the exact digest-role code for this nominal schema",
                    observed.role_code(),
                ));
            }
            let observed_owner = ConstructorOwnerV1::from_code(observed.owner_code())?;
            if observed_owner != expected.owner() {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Incompatible,
                    "constructor_owner_handoff.presented.owner_code",
                    "the exact constructor-owner code for this nominal schema",
                    observed.owner_code(),
                ));
            }
            if observed.owner_name() != expected.owner().name() {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Incompatible,
                    "constructor_owner_handoff.presented.owner_name",
                    expected.owner().name(),
                    ConstructionObservedV2::fixed(
                        ConstructionFixedObservationV2::StaleOrSubstitutedOwnerName,
                    ),
                ));
            }
            if observed.owner_id() != expected.owner().owner_id() {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Incompatible,
                    "constructor_owner_handoff.presented.owner_id",
                    expected.owner().owner_id(),
                    ConstructionObservedV2::fixed(
                        ConstructionFixedObservationV2::StaleOrSubstitutedOwnerIdentifier,
                    ),
                ));
            }
        }

        Self::reconstruct_exact_sequence(&FROZEN_CONSTRUCTOR_OWNER_HANDOFF_ENTRIES_V1)
    }

    /// Construct the exact crate-owned frozen mapping.
    #[must_use]
    pub fn frozen() -> Self {
        Self::reconstruct_exact_set(&FROZEN_CONSTRUCTOR_OWNER_HANDOFF_ENTRIES_V1)
            .expect("the source-defined constructor-owner mapping is internally valid")
    }

    /// Canonically ordered exact mapping.
    #[must_use]
    pub fn entries(&self) -> &[ConstructorOwnerHandoffEntryV1] {
        &self.entries
    }

    /// Look up the sole constructor owner by exact nominal schema name.
    ///
    /// # Errors
    ///
    /// Unknown or mutated schema names refuse rather than selecting a default
    /// owner.
    pub fn owner_for_schema_name(
        &self,
        schema_name: &str,
    ) -> Result<ConstructorOwnerV1, ConstructionErrorV2> {
        self.entries
            .iter()
            .find(|entry| entry.descriptor.schema_name() == schema_name)
            .map(|entry| entry.owner)
            .ok_or_else(|| {
                ConstructionErrorV2::new_redacted(
                    ConstructionErrorKindV2::UnknownCode,
                    "constructor_owner_handoff.schema_name",
                    "one exact frozen nominal schema name",
                    ConstructionObservedDataClassV2::CallerControlledText,
                )
            })
    }

    /// Domain-separated canonical identity of every schema, domain, role, and
    /// owner row.
    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }
}

fn constructor_owner_handoff_root(
    entries: &[ConstructorOwnerHandoffEntryV1],
) -> Result<ContentHash, ConstructionErrorV2> {
    let mut frame = CanonicalFrameV1::new(b"FSCOHPR\x01", 32 * 1024)?;
    frame.push_u32(
        "constructor_owner_handoff.count",
        u32::try_from(entries.len()).map_err(|_| {
            ConstructionErrorV2::new(
                ConstructionErrorKindV2::TooLarge,
                "constructor_owner_handoff.count",
                "a count representable as u32",
                entries.len(),
            )
        })?,
    )?;
    for entry in entries {
        frame.push_str(
            "constructor_owner_handoff.schema_name",
            entry.descriptor.schema_name(),
        )?;
        frame.push_str(
            "constructor_owner_handoff.domain",
            entry.descriptor.domain(),
        )?;
        frame.push_u16(
            "constructor_owner_handoff.role",
            entry.descriptor.role().code(),
        )?;
        frame.push_u16("constructor_owner_handoff.owner", entry.owner.code())?;
        frame.push_str("constructor_owner_handoff.owner_id", entry.owner.owner_id())?;
    }
    Ok(frame.root(CONSTRUCTOR_OWNER_HANDOFF_PROJECTION_DOMAIN_V1))
}

/// Domain for the exact root-free evaluator-member guard projection.
pub const ROOT_FREE_EVALUATOR_MEMBER_GUARD_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.root-free-evaluator-member-guard.v1";

/// The four evaluator-owned nested member schemas that deliberately have no
/// standalone nominal root.
///
/// Attempting to import the forbidden planning wrappers fails:
///
/// ```compile_fail,E0432
/// use fs_evidence_runner::identity::{
///     CaseEvaluationKeyRootV2,
///     EvaluationPrecontextRootV2,
///     ExpectedEvaluationErrorFrameRootV2,
///     ExpectedMismatchFrameRootV2,
/// };
/// ```
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
#[repr(u16)]
pub enum RootFreeEvaluatorMemberV1 {
    /// Root-free case-evaluation key nested under evaluator/conformance data.
    CaseEvaluationKey = 1,
    /// Root-free evaluation precontext nested under evaluator/conformance data.
    EvaluationPrecontext = 2,
    /// Root-free expected mismatch frame.
    ExpectedMismatchFrame = 3,
    /// Root-free expected evaluation-error frame.
    ExpectedEvaluationErrorFrame = 4,
}

impl RootFreeEvaluatorMemberV1 {
    /// Exact member inventory and order.
    pub const ALL: [Self; 4] = [
        Self::CaseEvaluationKey,
        Self::EvaluationPrecontext,
        Self::ExpectedMismatchFrame,
        Self::ExpectedEvaluationErrorFrame,
    ];

    /// Exact non-wire guard code.
    #[must_use]
    pub const fn code(self) -> u16 {
        self as u16
    }

    /// Exact downstream Rust schema name.
    #[must_use]
    pub const fn schema_name(self) -> &'static str {
        match self {
            Self::CaseEvaluationKey => "CaseEvaluationKeyV2",
            Self::EvaluationPrecontext => "EvaluationPrecontextV2",
            Self::ExpectedMismatchFrame => "ExpectedMismatchFrameV2",
            Self::ExpectedEvaluationErrorFrame => "ExpectedEvaluationErrorFrameV2",
        }
    }
}

impl ConstructionClosedSemanticV2 for RootFreeEvaluatorMemberV1 {
    fn construction_stable_name(&self) -> &'static str {
        self.schema_name()
    }
}

/// One exact source-owned guard descriptor for a downstream root-free member.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct RootFreeEvaluatorMemberGuardDescriptorV1 {
    member: RootFreeEvaluatorMemberV1,
    forbidden_root_schema: &'static str,
    dependency_rank: u16,
    predecessors: &'static [RootFreeEvaluatorMemberV1],
    allowed_nominal_inputs: &'static [&'static str],
}

impl RootFreeEvaluatorMemberGuardDescriptorV1 {
    const fn new(
        member: RootFreeEvaluatorMemberV1,
        forbidden_root_schema: &'static str,
        dependency_rank: u16,
        predecessors: &'static [RootFreeEvaluatorMemberV1],
        allowed_nominal_inputs: &'static [&'static str],
    ) -> Self {
        Self {
            member,
            forbidden_root_schema,
            dependency_rank,
            predecessors,
            allowed_nominal_inputs,
        }
    }

    /// Exact root-free member schema.
    #[must_use]
    pub const fn member(self) -> RootFreeEvaluatorMemberV1 {
        self.member
    }

    /// Sole downstream semantic owner.
    #[must_use]
    pub const fn owner(self) -> ConstructorOwnerV1 {
        let _ = self;
        ConstructorOwnerV1::EvaluatorAndCaseConformance
    }

    /// Exact forbidden standalone wrapper spelling.
    #[must_use]
    pub const fn forbidden_root_schema(self) -> &'static str {
        self.forbidden_root_schema
    }

    /// Topological dependency rank.
    #[must_use]
    pub const fn dependency_rank(self) -> u16 {
        self.dependency_rank
    }

    /// Exact root-free-member predecessors.
    #[must_use]
    pub const fn predecessors(self) -> &'static [RootFreeEvaluatorMemberV1] {
        self.predecessors
    }

    /// Exact nominal input roots allowed to flow into this nested member.
    #[must_use]
    pub const fn allowed_nominal_inputs(self) -> &'static [&'static str] {
        self.allowed_nominal_inputs
    }
}

const PRECONTEXT_PREDECESSORS_V1: [RootFreeEvaluatorMemberV1; 1] =
    [RootFreeEvaluatorMemberV1::CaseEvaluationKey];
const EXPECTED_FRAME_PREDECESSORS_V1: [RootFreeEvaluatorMemberV1; 1] =
    [RootFreeEvaluatorMemberV1::EvaluationPrecontext];
const CASE_KEY_INPUTS_V1: [&str; 1] = ["evaluation-registry-view"];
const PRECONTEXT_INPUTS_V1: [&str; 6] = [
    "comparison-expression",
    "effect-snapshot-requirements",
    "presented-observation-batch",
    "observation-set",
    "presented-effect-snapshot",
    "effect-snapshot",
];

/// Actual-result and authority roots that can never be folded into a
/// root-free expectation member.
pub const FORBIDDEN_ROOT_FREE_ACTUAL_INPUTS_V1: [&str; 11] = [
    "evaluation-context",
    "evaluation-error",
    "mismatch-detail",
    "case-manifest",
    "case-raw-outcome",
    "case-conformance-verdict",
    "lifecycle-log",
    "run-summary",
    "run-terminal-record",
    "authority-scope",
    "resource-identity",
];

/// Exact source-defined four-row root-free member guard inventory.
pub const FROZEN_ROOT_FREE_EVALUATOR_MEMBER_GUARDS_V1: [RootFreeEvaluatorMemberGuardDescriptorV1;
    4] = [
    RootFreeEvaluatorMemberGuardDescriptorV1::new(
        RootFreeEvaluatorMemberV1::CaseEvaluationKey,
        "CaseEvaluationKeyRootV2",
        0,
        &[],
        &CASE_KEY_INPUTS_V1,
    ),
    RootFreeEvaluatorMemberGuardDescriptorV1::new(
        RootFreeEvaluatorMemberV1::EvaluationPrecontext,
        "EvaluationPrecontextRootV2",
        1,
        &PRECONTEXT_PREDECESSORS_V1,
        &PRECONTEXT_INPUTS_V1,
    ),
    RootFreeEvaluatorMemberGuardDescriptorV1::new(
        RootFreeEvaluatorMemberV1::ExpectedMismatchFrame,
        "ExpectedMismatchFrameRootV2",
        2,
        &EXPECTED_FRAME_PREDECESSORS_V1,
        &[],
    ),
    RootFreeEvaluatorMemberGuardDescriptorV1::new(
        RootFreeEvaluatorMemberV1::ExpectedEvaluationErrorFrame,
        "ExpectedEvaluationErrorFrameRootV2",
        2,
        &EXPECTED_FRAME_PREDECESSORS_V1,
        &[],
    ),
];

/// Untrusted, non-wire row presented for exact root-free guard
/// reconstruction.
///
/// Optional root/domain/role/parser fields exist only so a checker can reject
/// fabricated widening attempts. A successfully reconstructed guard always
/// has all four absent and both wildcard booleans false.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PresentedRootFreeEvaluatorMemberGuardV1 {
    member: RootFreeEvaluatorMemberV1,
    owner: ConstructorOwnerV1,
    forbidden_root_schema: String,
    dependency_rank: u16,
    predecessors: Box<[RootFreeEvaluatorMemberV1]>,
    allowed_nominal_inputs: Box<[String]>,
    standalone_root_schema: Option<String>,
    digest_role: Option<DigestRoleV2>,
    standalone_domain: Option<String>,
    presented_parser: Option<String>,
    generic_digest_allowed: bool,
    wildcard_members_allowed: bool,
}

impl PresentedRootFreeEvaluatorMemberGuardV1 {
    /// Construct an explicitly presented guard row. This does not admit it.
    #[allow(
        clippy::too_many_arguments,
        reason = "every widening surface is explicit so exact reconstruction can reject it independently"
    )]
    pub fn new(
        member: RootFreeEvaluatorMemberV1,
        owner: ConstructorOwnerV1,
        forbidden_root_schema: impl Into<String>,
        dependency_rank: u16,
        predecessors: Vec<RootFreeEvaluatorMemberV1>,
        allowed_nominal_inputs: Vec<String>,
        standalone_root_schema: Option<String>,
        digest_role: Option<DigestRoleV2>,
        standalone_domain: Option<String>,
        presented_parser: Option<String>,
        generic_digest_allowed: bool,
        wildcard_members_allowed: bool,
    ) -> Self {
        Self {
            member,
            owner,
            forbidden_root_schema: forbidden_root_schema.into(),
            dependency_rank,
            predecessors: predecessors.into_boxed_slice(),
            allowed_nominal_inputs: allowed_nominal_inputs.into_boxed_slice(),
            standalone_root_schema,
            digest_role,
            standalone_domain,
            presented_parser,
            generic_digest_allowed,
            wildcard_members_allowed,
        }
    }

    fn from_frozen(descriptor: RootFreeEvaluatorMemberGuardDescriptorV1) -> Self {
        Self::new(
            descriptor.member(),
            descriptor.owner(),
            descriptor.forbidden_root_schema(),
            descriptor.dependency_rank(),
            descriptor.predecessors().to_vec(),
            descriptor
                .allowed_nominal_inputs()
                .iter()
                .map(|value| (*value).to_owned())
                .collect(),
            None,
            None,
            None,
            None,
            false,
            false,
        )
    }

    /// Presented member identity.
    #[must_use]
    pub const fn member(&self) -> RootFreeEvaluatorMemberV1 {
        self.member
    }
}

/// Exact, deterministic guard projection for four downstream root-free
/// evaluator members.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RootFreeEvaluatorMemberGuardProjectionV1 {
    rows: Box<[RootFreeEvaluatorMemberGuardDescriptorV1]>,
    root: ContentHash,
}

impl RootFreeEvaluatorMemberGuardProjectionV1 {
    /// Construct the exact frozen guard.
    #[must_use]
    pub fn frozen() -> Self {
        let rows = FROZEN_ROOT_FREE_EVALUATOR_MEMBER_GUARDS_V1
            .iter()
            .copied()
            .map(PresentedRootFreeEvaluatorMemberGuardV1::from_frozen)
            .collect::<Vec<_>>();
        Self::reconstruct_exact(&rows)
            .expect("the source-defined root-free evaluator guard is internally valid")
    }

    /// Reconstruct the exact ordered four-row guard.
    ///
    /// Missing, extra, duplicate, reordered, owner-mutated, reverse-edge,
    /// fabricated-root, generic-digest, parser, wildcard, and actual-root
    /// presentations refuse.
    pub fn reconstruct_exact(
        presented: &[PresentedRootFreeEvaluatorMemberGuardV1],
    ) -> Result<Self, ConstructionErrorV2> {
        if presented.len() != FROZEN_ROOT_FREE_EVALUATOR_MEMBER_GUARDS_V1.len() {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::OutOfRange,
                "root_free_guard.count",
                "exactly four ordered root-free evaluator member rows",
                presented.len(),
            ));
        }

        let mut seen = std::collections::BTreeSet::new();
        for row in presented {
            if !seen.insert(row.member) {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Duplicate,
                    "root_free_guard.member",
                    "one row per exact root-free member",
                    ConstructionObservedV2::closed(&row.member),
                ));
            }
        }

        for (index, (row, expected)) in presented
            .iter()
            .zip(FROZEN_ROOT_FREE_EVALUATOR_MEMBER_GUARDS_V1)
            .enumerate()
        {
            if row.member != expected.member {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::OutOfOrder,
                    "root_free_guard.member",
                    "the exact four-row controlling member order",
                    index,
                ));
            }
            if row.owner != ConstructorOwnerV1::EvaluatorAndCaseConformance {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Incompatible,
                    "root_free_guard.owner",
                    "EvaluatorAndCaseConformance",
                    ConstructionObservedV2::closed(&row.owner),
                ));
            }
            if row.forbidden_root_schema != expected.forbidden_root_schema {
                return Err(ConstructionErrorV2::new_redacted(
                    ConstructionErrorKindV2::Incompatible,
                    "root_free_guard.forbidden_root_schema",
                    expected.forbidden_root_schema,
                    ConstructionObservedDataClassV2::CallerControlledText,
                ));
            }
            if row.dependency_rank != expected.dependency_rank {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Incompatible,
                    "root_free_guard.dependency_rank",
                    "the exact acyclic dependency rank",
                    row.dependency_rank,
                ));
            }
            for predecessor in &row.predecessors {
                if predecessor.code() >= row.member.code() {
                    return Err(ConstructionErrorV2::new(
                        ConstructionErrorKindV2::Incompatible,
                        "root_free_guard.predecessors",
                        "strictly earlier root-free predecessors only",
                        ConstructionObservedV2::closed(predecessor),
                    ));
                }
            }
            if row.predecessors.as_ref() != expected.predecessors {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Incompatible,
                    "root_free_guard.predecessors",
                    "the exact dependency-DAG predecessor sequence",
                    row.predecessors.len(),
                ));
            }
            let expected_inputs = expected.allowed_nominal_inputs;
            if row.allowed_nominal_inputs.len() != expected_inputs.len()
                || row
                    .allowed_nominal_inputs
                    .iter()
                    .map(String::as_str)
                    .ne(expected_inputs.iter().copied())
            {
                return Err(ConstructionErrorV2::new_redacted(
                    ConstructionErrorKindV2::Incompatible,
                    "root_free_guard.allowed_nominal_inputs",
                    "the exact forward-only nominal input-root sequence",
                    ConstructionObservedDataClassV2::BulkPayload,
                ));
            }
            if row.standalone_root_schema.is_some()
                || row.digest_role.is_some()
                || row.standalone_domain.is_some()
                || row.presented_parser.is_some()
            {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Unexpected,
                    "root_free_guard.standalone_identity_surface",
                    "no root, role, domain, or presented parser",
                    ConstructionObservedV2::closed(&row.member),
                ));
            }
            if row.generic_digest_allowed || row.wildcard_members_allowed {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Unexpected,
                    "root_free_guard.widening_surface",
                    "no generic digest conversion or wildcard member",
                    ConstructionObservedV2::closed(&row.member),
                ));
            }
        }

        let rows = FROZEN_ROOT_FREE_EVALUATOR_MEMBER_GUARDS_V1;
        let root = root_free_evaluator_member_guard_root(&rows)?;
        Ok(Self {
            rows: rows.into(),
            root,
        })
    }

    /// Exact four-row descriptor inventory.
    #[must_use]
    pub fn rows(&self) -> &[RootFreeEvaluatorMemberGuardDescriptorV1] {
        &self.rows
    }

    /// Domain-separated guard root.
    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }
}

fn root_free_evaluator_member_guard_root(
    rows: &[RootFreeEvaluatorMemberGuardDescriptorV1],
) -> Result<ContentHash, ConstructionErrorV2> {
    let mut frame = CanonicalFrameV1::new(b"FSRFMGRD\x01", 32 * 1024)?;
    frame.push_u32(
        "root_free_guard.count",
        u32::try_from(rows.len()).map_err(|_| {
            ConstructionErrorV2::new(
                ConstructionErrorKindV2::TooLarge,
                "root_free_guard.count",
                "a count representable as u32",
                rows.len(),
            )
        })?,
    )?;
    for row in rows {
        frame.push_u16("root_free_guard.member", row.member.code())?;
        frame.push_str("root_free_guard.schema_name", row.member.schema_name())?;
        frame.push_u16("root_free_guard.owner", row.owner().code())?;
        frame.push_str("root_free_guard.owner_id", row.owner().owner_id())?;
        frame.push_str(
            "root_free_guard.forbidden_root_schema",
            row.forbidden_root_schema,
        )?;
        frame.push_u16("root_free_guard.dependency_rank", row.dependency_rank)?;
        frame.push_u32(
            "root_free_guard.predecessor_count",
            u32::try_from(row.predecessors.len()).expect("bounded frozen predecessor count"),
        )?;
        for predecessor in row.predecessors {
            frame.push_u16("root_free_guard.predecessor", predecessor.code())?;
        }
        frame.push_u32(
            "root_free_guard.allowed_input_count",
            u32::try_from(row.allowed_nominal_inputs.len())
                .expect("bounded frozen allowed-input count"),
        )?;
        for input in row.allowed_nominal_inputs {
            frame.push_str("root_free_guard.allowed_input", input)?;
        }
        frame.push_u32(
            "root_free_guard.forbidden_actual_input_count",
            u32::try_from(FORBIDDEN_ROOT_FREE_ACTUAL_INPUTS_V1.len())
                .expect("bounded frozen forbidden-input count"),
        )?;
        for forbidden in FORBIDDEN_ROOT_FREE_ACTUAL_INPUTS_V1 {
            frame.push_str("root_free_guard.forbidden_actual_input", forbidden)?;
        }
        frame.push_u16("root_free_guard.standalone_identity_surface", 0)?;
        frame.push_u16("root_free_guard.generic_digest_allowed", 0)?;
        frame.push_u16("root_free_guard.wildcard_members_allowed", 0)?;
    }
    Ok(frame.root(ROOT_FREE_EVALUATOR_MEMBER_GUARD_DOMAIN_V1))
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeSet;

    use super::*;

    fn presented_root_free_guard_rows() -> Vec<PresentedRootFreeEvaluatorMemberGuardV1> {
        FROZEN_ROOT_FREE_EVALUATOR_MEMBER_GUARDS_V1
            .iter()
            .copied()
            .map(PresentedRootFreeEvaluatorMemberGuardV1::from_frozen)
            .collect()
    }

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
        assert_eq!(ALL_PRESENTED_IDENTITY_DESCRIPTORS_V1.len(), 43);
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
            ("comparison-expression", DigestRoleV2::Spec),
            ("effect-expectation", DigestRoleV2::Spec),
            ("presented-observation-batch", DigestRoleV2::Run),
            ("observation-set", DigestRoleV2::Run),
            ("presented-effect-snapshot", DigestRoleV2::Run),
            ("effect-snapshot-requirements", DigestRoleV2::Spec),
            ("effect-snapshot", DigestRoleV2::Run),
            ("evaluation-registry-view", DigestRoleV2::Policy),
            ("evaluation-context", DigestRoleV2::Run),
            ("evaluation-error", DigestRoleV2::Run),
            ("mismatch-detail", DigestRoleV2::Run),
            ("case-conformance-spec", DigestRoleV2::Spec),
            ("case-raw-outcome", DigestRoleV2::Run),
            ("case-conformance-verdict", DigestRoleV2::Run),
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
        assert_presented_wrapper!(ComparisonExprRootV2, 29);
        assert_presented_wrapper!(EffectExpectationRootV2, 30);
        assert_presented_wrapper!(PresentedObservationBatchRootV2, 31);
        assert_presented_wrapper!(ObservationSetRootV2, 32);
        assert_presented_wrapper!(PresentedEffectSnapshotRootV2, 33);
        assert_presented_wrapper!(EffectSnapshotRequirementsRootV2, 34);
        assert_presented_wrapper!(EffectSnapshotRootV2, 35);
        assert_presented_wrapper!(EvaluationRegistryViewRootV2, 36);
        assert_presented_wrapper!(EvaluationContextRootV2, 37);
        assert_presented_wrapper!(EvaluationErrorRootV2, 38);
        assert_presented_wrapper!(MismatchDetailRootV2, 39);
        assert_presented_wrapper!(CaseConformanceSpecRootV2, 40);
        assert_presented_wrapper!(CaseRawOutcomeRootV2, 41);
        assert_presented_wrapper!(CaseConformanceVerdictRootV2, 42);
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

    #[test]
    fn constructor_owner_handoff_exactly_covers_every_nominal_descriptor() {
        let projection = ConstructorOwnerHandoffProjectionV1::frozen();
        assert_eq!(
            projection.entries().len(),
            ALL_PRESENTED_IDENTITY_DESCRIPTORS_V1.len()
        );
        for (entry, descriptor) in projection
            .entries()
            .iter()
            .zip(ALL_PRESENTED_IDENTITY_DESCRIPTORS_V1)
        {
            assert_eq!(entry.descriptor(), descriptor);
            assert_eq!(
                projection
                    .owner_for_schema_name(descriptor.schema_name())
                    .expect("known handoff schema"),
                entry.owner()
            );
            assert!(
                !entry.owner().name().is_empty() && !entry.owner().owner_id().is_empty(),
                "{} must name its sole semantic-constructor owner",
                descriptor.schema_name()
            );
        }
        let represented = projection
            .entries()
            .iter()
            .map(|entry| entry.owner())
            .collect::<BTreeSet<_>>();
        assert_eq!(
            represented,
            ConstructorOwnerV1::ALL.into_iter().collect::<BTreeSet<_>>()
        );
        assert_eq!(
            projection
                .owner_for_schema_name("unknown-schema")
                .expect_err("no default owner")
                .kind(),
            ConstructionErrorKindV2::UnknownCode
        );
        for schema in ["artifact-encoded", "artifact-content", "stored-object"] {
            assert_eq!(
                projection
                    .owner_for_schema_name(schema)
                    .expect("byte-root handoff"),
                ConstructorOwnerV1::CanonicalBytesAndCodec
            );
        }
        for schema in [
            "external-mutation-set",
            "artifact-set",
            "comparison-expression",
            "effect-expectation",
            "presented-observation-batch",
            "observation-set",
            "presented-effect-snapshot",
            "effect-snapshot-requirements",
            "effect-snapshot",
            "evaluation-registry-view",
            "evaluation-context",
            "evaluation-error",
            "mismatch-detail",
            "case-conformance-spec",
            "case-raw-outcome",
            "case-conformance-verdict",
        ] {
            assert_eq!(
                projection
                    .owner_for_schema_name(schema)
                    .expect("evaluator/conformance handoff"),
                ConstructorOwnerV1::EvaluatorAndCaseConformance
            );
        }
        for superseded in [
            "observation-input",
            "effect-snapshot-input",
            "case-evaluation-registry",
        ] {
            assert_eq!(
                projection
                    .owner_for_schema_name(superseded)
                    .expect_err("superseded planning name must not survive")
                    .kind(),
                ConstructionErrorKindV2::UnknownCode
            );
        }
    }

    #[test]
    fn constructor_owner_handoff_reconstruction_is_set_exact_and_order_stable() {
        let frozen = ConstructorOwnerHandoffProjectionV1::frozen();
        let mut reordered = frozen.entries().to_vec();
        reordered.reverse();
        let reconstructed = ConstructorOwnerHandoffProjectionV1::reconstruct_exact_set(&reordered)
            .expect("set reconstruction canonicalizes order");
        assert_eq!(reconstructed.entries(), frozen.entries());
        assert_eq!(reconstructed.root(), frozen.root());

        let missing = &frozen.entries()[1..];
        assert_eq!(
            ConstructorOwnerHandoffProjectionV1::reconstruct_exact_set(missing)
                .expect_err("missing descriptor")
                .kind(),
            ConstructionErrorKindV2::OutOfRange
        );

        let mut duplicate = frozen.entries().to_vec();
        duplicate[1] = duplicate[0];
        assert_eq!(
            ConstructorOwnerHandoffProjectionV1::reconstruct_exact_set(&duplicate)
                .expect_err("duplicate descriptor")
                .kind(),
            ConstructionErrorKindV2::Duplicate
        );

        let mut owner_mutated = frozen.entries().to_vec();
        owner_mutated[0] = ConstructorOwnerHandoffEntryV1::new(
            owner_mutated[0].descriptor(),
            ConstructorOwnerV1::BaseSchema,
        );
        assert_eq!(
            ConstructorOwnerHandoffProjectionV1::reconstruct_exact_set(&owner_mutated)
                .expect_err("owner mutation")
                .kind(),
            ConstructionErrorKindV2::Incompatible
        );
    }

    #[test]
    fn constructor_owner_handoff_ordered_closeout_rejects_reordering_without_weakening_exact_set() {
        let frozen = ConstructorOwnerHandoffProjectionV1::frozen();
        assert_eq!(
            ConstructorOwnerHandoffProjectionV1::reconstruct_exact_sequence(frozen.entries())
                .expect("the frozen closeout sequence"),
            frozen
        );

        let mut reordered = frozen.entries().to_vec();
        reordered.reverse();
        assert_eq!(
            ConstructorOwnerHandoffProjectionV1::reconstruct_exact_set(&reordered)
                .expect("the unordered exact-set API remains permutation-invariant"),
            frozen
        );
        let error = ConstructorOwnerHandoffProjectionV1::reconstruct_exact_sequence(&reordered)
            .expect_err("the ordered closeout API rejects the same permutation");
        assert_eq!(error.kind(), ConstructionErrorKindV2::OutOfOrder);
        assert_eq!(error.field(), "constructor_owner_handoff.entries");
    }

    #[test]
    fn presented_constructor_owner_handoff_refuses_every_stale_or_unknown_identity_field() {
        let frozen = ConstructorOwnerHandoffProjectionV1::frozen();
        let presented = ConstructorOwnerHandoffProjectionV1::frozen_presented_sequence();
        assert_eq!(presented.len(), 43);
        assert_eq!(
            ConstructorOwnerHandoffProjectionV1::reconstruct_presented_exact_sequence(&presented)
                .expect("exact raw presentation"),
            frozen
        );

        let first = &presented[0];
        assert_eq!(
            (
                first.schema_name(),
                first.domain(),
                first.role_code(),
                first.owner_code(),
                first.owner_name(),
                first.owner_id(),
            ),
            (
                SourceIdentityRootV2::DESCRIPTOR.schema_name(),
                SourceIdentityRootV2::DESCRIPTOR.domain(),
                DigestRoleV2::Source.code(),
                ConstructorOwnerV1::InvocationAndRunIdentity.code(),
                ConstructorOwnerV1::InvocationAndRunIdentity.name(),
                ConstructorOwnerV1::InvocationAndRunIdentity.owner_id(),
            )
        );

        let mutations = [
            (
                PresentedConstructorOwnerHandoffEntryV1::new(
                    "unknown-source-identity",
                    first.domain(),
                    first.role_code(),
                    first.owner_code(),
                    first.owner_name(),
                    first.owner_id(),
                ),
                ConstructionErrorKindV2::UnknownCode,
                "constructor_owner_handoff.presented.schema_name",
            ),
            (
                PresentedConstructorOwnerHandoffEntryV1::new(
                    first.schema_name(),
                    "org.frankensim.fs-evidence-runner.source-identity.stale.v1",
                    first.role_code(),
                    first.owner_code(),
                    first.owner_name(),
                    first.owner_id(),
                ),
                ConstructionErrorKindV2::Incompatible,
                "constructor_owner_handoff.presented.domain",
            ),
            (
                PresentedConstructorOwnerHandoffEntryV1::new(
                    first.schema_name(),
                    first.domain(),
                    u16::MAX,
                    first.owner_code(),
                    first.owner_name(),
                    first.owner_id(),
                ),
                ConstructionErrorKindV2::UnknownCode,
                "constructor_owner_handoff.presented.role_code",
            ),
            (
                PresentedConstructorOwnerHandoffEntryV1::new(
                    first.schema_name(),
                    first.domain(),
                    DigestRoleV2::Build.code(),
                    first.owner_code(),
                    first.owner_name(),
                    first.owner_id(),
                ),
                ConstructionErrorKindV2::Incompatible,
                "constructor_owner_handoff.presented.role_code",
            ),
            (
                PresentedConstructorOwnerHandoffEntryV1::new(
                    first.schema_name(),
                    first.domain(),
                    first.role_code(),
                    u16::MAX,
                    first.owner_name(),
                    first.owner_id(),
                ),
                ConstructionErrorKindV2::UnknownCode,
                "constructor_owner_handoff.presented.owner_code",
            ),
            (
                PresentedConstructorOwnerHandoffEntryV1::new(
                    first.schema_name(),
                    first.domain(),
                    first.role_code(),
                    ConstructorOwnerV1::BaseSchema.code(),
                    ConstructorOwnerV1::BaseSchema.name(),
                    ConstructorOwnerV1::BaseSchema.owner_id(),
                ),
                ConstructionErrorKindV2::Incompatible,
                "constructor_owner_handoff.presented.owner_code",
            ),
            (
                PresentedConstructorOwnerHandoffEntryV1::new(
                    first.schema_name(),
                    first.domain(),
                    first.role_code(),
                    first.owner_code(),
                    "stale-owner-name",
                    first.owner_id(),
                ),
                ConstructionErrorKindV2::Incompatible,
                "constructor_owner_handoff.presented.owner_name",
            ),
            (
                PresentedConstructorOwnerHandoffEntryV1::new(
                    first.schema_name(),
                    first.domain(),
                    first.role_code(),
                    first.owner_code(),
                    first.owner_name(),
                    "frankensim-stale-owner",
                ),
                ConstructionErrorKindV2::Incompatible,
                "constructor_owner_handoff.presented.owner_id",
            ),
        ];
        for (mutation, expected_kind, expected_field) in mutations {
            let mut rows = presented.to_vec();
            rows[0] = mutation;
            let error =
                ConstructorOwnerHandoffProjectionV1::reconstruct_presented_exact_sequence(&rows)
                    .expect_err("raw-field mutation must refuse");
            assert_eq!(error.kind(), expected_kind);
            assert_eq!(error.field(), expected_field);
        }

        let mut reordered = presented.to_vec();
        reordered.swap(0, 1);
        let reorder_error =
            ConstructorOwnerHandoffProjectionV1::reconstruct_presented_exact_sequence(&reordered)
                .expect_err("known schema at the wrong sequence position");
        assert_eq!(reorder_error.kind(), ConstructionErrorKindV2::OutOfOrder);
        assert_eq!(
            reorder_error.field(),
            "constructor_owner_handoff.presented.schema_name"
        );

        let mut duplicate = presented.to_vec();
        duplicate[1] = duplicate[0].clone();
        let duplicate_error =
            ConstructorOwnerHandoffProjectionV1::reconstruct_presented_exact_sequence(&duplicate)
                .expect_err("duplicate schema name");
        assert_eq!(duplicate_error.kind(), ConstructionErrorKindV2::Duplicate);
        assert_eq!(
            duplicate_error.field(),
            "constructor_owner_handoff.presented.schema_name"
        );

        let missing_error =
            ConstructorOwnerHandoffProjectionV1::reconstruct_presented_exact_sequence(
                &presented[1..],
            )
            .expect_err("missing row");
        assert_eq!(missing_error.kind(), ConstructionErrorKindV2::OutOfRange);
        assert_eq!(
            missing_error.field(),
            "constructor_owner_handoff.presented.entries"
        );
        let mut extra = presented.to_vec();
        extra.push(presented[0].clone());
        let extra_error =
            ConstructorOwnerHandoffProjectionV1::reconstruct_presented_exact_sequence(&extra)
                .expect_err("extra row");
        assert_eq!(extra_error.kind(), ConstructionErrorKindV2::OutOfRange);
        assert_eq!(
            extra_error.field(),
            "constructor_owner_handoff.presented.entries"
        );
    }

    #[test]
    fn owner_schema_domain_and_role_mutations_move_the_handoff_root_and_refuse() {
        let frozen = ConstructorOwnerHandoffProjectionV1::frozen();
        let original = frozen.entries()[0];
        let descriptor = original.descriptor();
        let mutations = [
            ConstructorOwnerHandoffEntryV1::new(descriptor, ConstructorOwnerV1::BaseSchema),
            ConstructorOwnerHandoffEntryV1::new(
                PresentedIdentityDescriptorV1::new(
                    "source-identity-mutated",
                    descriptor.domain(),
                    descriptor.role(),
                ),
                original.owner(),
            ),
            ConstructorOwnerHandoffEntryV1::new(
                PresentedIdentityDescriptorV1::new(
                    descriptor.schema_name(),
                    "org.frankensim.fs-evidence-runner.source-identity-mutated.v1",
                    descriptor.role(),
                ),
                original.owner(),
            ),
            ConstructorOwnerHandoffEntryV1::new(
                PresentedIdentityDescriptorV1::new(
                    descriptor.schema_name(),
                    descriptor.domain(),
                    DigestRoleV2::Build,
                ),
                original.owner(),
            ),
        ];

        for mutation in mutations {
            let mut rows = frozen.entries().to_vec();
            rows[0] = mutation;
            assert_ne!(
                constructor_owner_handoff_root(&rows).expect("bounded mutation root"),
                frozen.root()
            );
            assert!(
                ConstructorOwnerHandoffProjectionV1::reconstruct_exact_set(&rows).is_err(),
                "owner, schema, domain, and role mutations must refuse exact reconstruction"
            );
        }
    }

    #[test]
    fn root_free_guard_inventory_is_exact_ordered_and_root_stable() {
        let guard = RootFreeEvaluatorMemberGuardProjectionV1::frozen();
        assert_eq!(guard.rows(), &FROZEN_ROOT_FREE_EVALUATOR_MEMBER_GUARDS_V1);
        assert_eq!(
            guard
                .rows()
                .iter()
                .map(|row| row.member())
                .collect::<Vec<_>>(),
            RootFreeEvaluatorMemberV1::ALL
        );
        assert_eq!(
            guard
                .rows()
                .iter()
                .map(|row| row.member().schema_name())
                .collect::<Vec<_>>(),
            [
                "CaseEvaluationKeyV2",
                "EvaluationPrecontextV2",
                "ExpectedMismatchFrameV2",
                "ExpectedEvaluationErrorFrameV2",
            ]
        );
        assert_eq!(
            guard
                .rows()
                .iter()
                .map(|row| row.forbidden_root_schema())
                .collect::<Vec<_>>(),
            [
                "CaseEvaluationKeyRootV2",
                "EvaluationPrecontextRootV2",
                "ExpectedMismatchFrameRootV2",
                "ExpectedEvaluationErrorFrameRootV2",
            ]
        );
        assert_eq!(
            guard
                .rows()
                .iter()
                .map(|row| row.dependency_rank())
                .collect::<Vec<_>>(),
            [0, 1, 2, 2]
        );
        assert_eq!(guard.rows()[0].predecessors(), &[]);
        assert_eq!(
            guard.rows()[1].predecessors(),
            &[RootFreeEvaluatorMemberV1::CaseEvaluationKey]
        );
        assert_eq!(
            guard.rows()[2].predecessors(),
            &[RootFreeEvaluatorMemberV1::EvaluationPrecontext]
        );
        assert_eq!(
            guard.rows()[3].predecessors(),
            &[RootFreeEvaluatorMemberV1::EvaluationPrecontext]
        );
        assert_eq!(
            guard.rows()[0].allowed_nominal_inputs(),
            &["evaluation-registry-view"]
        );
        assert_eq!(
            guard.rows()[1].allowed_nominal_inputs(),
            &[
                "comparison-expression",
                "effect-snapshot-requirements",
                "presented-observation-batch",
                "observation-set",
                "presented-effect-snapshot",
                "effect-snapshot",
            ]
        );
        assert!(guard.rows()[2].allowed_nominal_inputs().is_empty());
        assert!(guard.rows()[3].allowed_nominal_inputs().is_empty());
        assert_eq!(FORBIDDEN_ROOT_FREE_ACTUAL_INPUTS_V1.len(), 11);
        assert_eq!(ALL_PRESENTED_IDENTITY_DESCRIPTORS_V1.len(), 43);
        let handoff = ConstructorOwnerHandoffProjectionV1::frozen();
        assert_eq!(
            ALL_PRESENTED_IDENTITY_DESCRIPTORS_V1[29..]
                .iter()
                .filter(|descriptor| {
                    handoff.owner_for_schema_name(descriptor.schema_name())
                        == Ok(ConstructorOwnerV1::EvaluatorAndCaseConformance)
                })
                .count(),
            14
        );

        let reconstructed = RootFreeEvaluatorMemberGuardProjectionV1::reconstruct_exact(
            &presented_root_free_guard_rows(),
        )
        .expect("the exact presented inventory reconstructs");
        assert_eq!(reconstructed, guard);
        assert_eq!(
            root_free_evaluator_member_guard_root(guard.rows()).expect("bounded frozen guard root"),
            guard.root()
        );
    }

    #[test]
    fn root_free_guard_rejects_missing_extra_duplicate_and_reordered_rows() {
        let rows = presented_root_free_guard_rows();

        assert_eq!(
            RootFreeEvaluatorMemberGuardProjectionV1::reconstruct_exact(&rows[1..])
                .expect_err("missing row must refuse")
                .kind(),
            ConstructionErrorKindV2::OutOfRange
        );

        let mut extra = rows.clone();
        extra.push(rows[0].clone());
        assert_eq!(
            RootFreeEvaluatorMemberGuardProjectionV1::reconstruct_exact(&extra)
                .expect_err("extra row must refuse")
                .kind(),
            ConstructionErrorKindV2::OutOfRange
        );

        let mut duplicate = rows.clone();
        duplicate[1] = duplicate[0].clone();
        assert_eq!(
            RootFreeEvaluatorMemberGuardProjectionV1::reconstruct_exact(&duplicate)
                .expect_err("duplicate member must refuse")
                .kind(),
            ConstructionErrorKindV2::Duplicate
        );

        let mut reordered = rows;
        reordered.swap(2, 3);
        assert_eq!(
            RootFreeEvaluatorMemberGuardProjectionV1::reconstruct_exact(&reordered)
                .expect_err("reordered members must refuse")
                .kind(),
            ConstructionErrorKindV2::OutOfOrder
        );
    }

    #[test]
    fn root_free_guard_rejects_owner_edge_and_context_collapsing_mutants() {
        let rows = presented_root_free_guard_rows();

        let mut owner_mutated = rows.clone();
        owner_mutated[0].owner = ConstructorOwnerV1::BaseSchema;
        let owner_error =
            RootFreeEvaluatorMemberGuardProjectionV1::reconstruct_exact(&owner_mutated)
                .expect_err("owner mutation must refuse");
        assert_eq!(owner_error.kind(), ConstructionErrorKindV2::Incompatible);
        assert_eq!(owner_error.field(), "root_free_guard.owner");

        let mut reverse_edge = rows.clone();
        reverse_edge[1].predecessors =
            vec![RootFreeEvaluatorMemberV1::ExpectedMismatchFrame].into_boxed_slice();
        let reverse_error =
            RootFreeEvaluatorMemberGuardProjectionV1::reconstruct_exact(&reverse_edge)
                .expect_err("reverse dependency must refuse");
        assert_eq!(reverse_error.kind(), ConstructionErrorKindV2::Incompatible);
        assert_eq!(reverse_error.field(), "root_free_guard.predecessors");

        let mut cycle = rows.clone();
        cycle[0].predecessors =
            vec![RootFreeEvaluatorMemberV1::EvaluationPrecontext].into_boxed_slice();
        let cycle_error = RootFreeEvaluatorMemberGuardProjectionV1::reconstruct_exact(&cycle)
            .expect_err("cycle-forming dependency must refuse");
        assert_eq!(cycle_error.kind(), ConstructionErrorKindV2::Incompatible);
        assert_eq!(cycle_error.field(), "root_free_guard.predecessors");

        for forbidden in FORBIDDEN_ROOT_FREE_ACTUAL_INPUTS_V1 {
            let mut context_collapsed = rows.clone();
            context_collapsed[1].allowed_nominal_inputs =
                vec![forbidden.to_owned()].into_boxed_slice();
            let error =
                RootFreeEvaluatorMemberGuardProjectionV1::reconstruct_exact(&context_collapsed)
                    .expect_err("actual-result or authority input must refuse");
            assert_eq!(error.kind(), ConstructionErrorKindV2::Incompatible);
            assert_eq!(error.field(), "root_free_guard.allowed_nominal_inputs");
        }
    }

    #[test]
    fn root_free_guard_rejects_every_fabricated_identity_and_widening_surface() {
        let rows = presented_root_free_guard_rows();

        let mutations = [
            {
                let mut mutation = rows.clone();
                mutation[0].standalone_root_schema = Some("CaseEvaluationKeyRootV2".to_owned());
                mutation
            },
            {
                let mut mutation = rows.clone();
                mutation[0].digest_role = Some(DigestRoleV2::Spec);
                mutation
            },
            {
                let mut mutation = rows.clone();
                mutation[0].standalone_domain =
                    Some("org.frankensim.fs-evidence-runner.case-evaluation-key.v1".to_owned());
                mutation
            },
            {
                let mut mutation = rows.clone();
                mutation[0].presented_parser = Some("parse_presented".to_owned());
                mutation
            },
        ];
        for mutation in mutations {
            let error = RootFreeEvaluatorMemberGuardProjectionV1::reconstruct_exact(&mutation)
                .expect_err("fabricated standalone identity surface must refuse");
            assert_eq!(error.kind(), ConstructionErrorKindV2::Unexpected);
            assert_eq!(error.field(), "root_free_guard.standalone_identity_surface");
        }

        let mut generic_digest = rows.clone();
        generic_digest[0].generic_digest_allowed = true;
        let generic_error =
            RootFreeEvaluatorMemberGuardProjectionV1::reconstruct_exact(&generic_digest)
                .expect_err("generic digest widening must refuse");
        assert_eq!(generic_error.kind(), ConstructionErrorKindV2::Unexpected);
        assert_eq!(generic_error.field(), "root_free_guard.widening_surface");

        let mut wildcard = rows;
        wildcard[0].wildcard_members_allowed = true;
        let wildcard_error = RootFreeEvaluatorMemberGuardProjectionV1::reconstruct_exact(&wildcard)
            .expect_err("wildcard member widening must refuse");
        assert_eq!(wildcard_error.kind(), ConstructionErrorKindV2::Unexpected);
        assert_eq!(wildcard_error.field(), "root_free_guard.widening_surface");
    }
}
