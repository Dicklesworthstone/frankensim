//! Exact standards-edition and source-lineage manifest.
//!
//! The source-reference and rule-provenance APIs deliberately retain only
//! identifiers, locators, hashes, access policy, and derived-rule provenance;
//! no API parameter accepts a protected standards body. A downstream standards
//! compiler must obtain licensed bytes out of band, hash them, and present that
//! observed hash at the rule-binding boundary.

use crate::ContentHash;
use core::fmt;
use fs_blake3::hash_domain;

/// Canonical wire schema for [`StandardManifest`].
pub const STANDARD_MANIFEST_SCHEMA_VERSION: u32 = 1;
/// Domain for the canonical manifest identity.
pub const STANDARD_MANIFEST_IDENTITY_DOMAIN: &str = "org.frankensim.fs-vvreg.standard-manifest.v1";
/// Domain for one canonical source-row identity.
pub const STANDARD_SOURCE_IDENTITY_DOMAIN: &str = "org.frankensim.fs-vvreg.standard-source.v1";
/// Domain for a sealed derived-rule provenance identity.
pub const STANDARD_RULE_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-vvreg.standard-rule-provenance.v1";

const WIRE_MAGIC: [u8; 4] = *b"FSMF";

/// Hard maximum number of standard editions in one manifest.
pub const MAX_STANDARD_RECORDS: usize = 4_096;
/// Hard maximum amendment/corrigendum rows on one edition.
pub const MAX_EDITION_CHANGES: usize = 64;
/// Hard maximum UTF-8 bytes in one manifest string.
pub const MAX_STANDARD_STRING_BYTES: usize = 4_096;
/// Hard maximum canonical manifest bytes.
pub const MAX_STANDARD_MANIFEST_BYTES: usize = 16 * 1_024 * 1_024;

/// Explicit decode/build envelope for a standards manifest.
///
/// Callers may tighten the compile-time ceilings but cannot relax them.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ManifestLimits {
    /// Maximum edition rows.
    pub records: usize,
    /// Maximum amendment/corrigendum rows per edition.
    pub changes_per_record: usize,
    /// Maximum bytes in any one UTF-8 string.
    pub string_bytes: usize,
    /// Maximum complete canonical payload bytes.
    pub manifest_bytes: usize,
}

impl ManifestLimits {
    /// Conservative default envelope.
    pub const DEFAULT: Self = Self {
        records: MAX_STANDARD_RECORDS,
        changes_per_record: MAX_EDITION_CHANGES,
        string_bytes: MAX_STANDARD_STRING_BYTES,
        manifest_bytes: MAX_STANDARD_MANIFEST_BYTES,
    };
}

impl Default for ManifestLimits {
    fn default() -> Self {
        Self::DEFAULT
    }
}

/// Stable lowercase identifier for a standards family.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct StandardId(String);

impl StandardId {
    /// Validate and retain a lowercase ASCII slug such as `iso-12100`.
    ///
    /// # Errors
    ///
    /// Refuses blank, malformed, or overlong identifiers.
    pub fn try_new(value: &str) -> Result<Self, ManifestError> {
        validate_identifier(value, "standard_id", MAX_STANDARD_STRING_BYTES)?;
        Ok(Self(value.to_string()))
    }

    /// Borrow the canonical slug.
    #[must_use]
    pub fn as_str(&self) -> &str {
        &self.0
    }
}

/// Exact `(standard, part, edition)` lookup key.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct StandardEditionKey {
    standard_id: StandardId,
    part: Option<String>,
    edition: String,
}

impl StandardEditionKey {
    /// Construct one exact edition key. `part=None` means a non-partitioned
    /// standard, not an unknown part.
    ///
    /// # Errors
    ///
    /// Refuses malformed identifiers, blank editions, or overlong strings.
    pub fn try_new(
        standard_id: &str,
        part: Option<&str>,
        edition: &str,
    ) -> Result<Self, ManifestError> {
        let standard_id = StandardId::try_new(standard_id)?;
        let part = match part {
            Some(part) => {
                validate_compact_text(part, "part", MAX_STANDARD_STRING_BYTES)?;
                Some(part.to_string())
            }
            None => None,
        };
        validate_compact_text(edition, "edition", MAX_STANDARD_STRING_BYTES)?;
        Ok(Self {
            standard_id,
            part,
            edition: edition.to_string(),
        })
    }

    /// Standards-family identifier.
    #[must_use]
    pub const fn standard_id(&self) -> &StandardId {
        &self.standard_id
    }

    /// Exact part identifier, when the standard is partitioned.
    #[must_use]
    pub fn part(&self) -> Option<&str> {
        self.part.as_deref()
    }

    /// Exact edition/revision label.
    #[must_use]
    pub fn edition(&self) -> &str {
        &self.edition
    }
}

/// Kind of official change incorporated after an edition was published.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EditionChangeKind {
    /// Amendment.
    Amendment,
    /// Technical corrigendum or erratum.
    Corrigendum,
}

impl EditionChangeKind {
    const fn tag(self) -> u8 {
        match self {
            Self::Amendment => 1,
            Self::Corrigendum => 2,
        }
    }
}

/// One ordered amendment or corrigendum incorporated into an exact edition.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EditionChange {
    kind: EditionChangeKind,
    designation: String,
}

impl EditionChange {
    /// Construct an exact change designation such as `Amd 1:2024`.
    ///
    /// # Errors
    ///
    /// Refuses blank, control-bearing, or overlong designations.
    pub fn try_new(kind: EditionChangeKind, designation: &str) -> Result<Self, ManifestError> {
        validate_compact_text(designation, "change.designation", MAX_STANDARD_STRING_BYTES)?;
        Ok(Self {
            kind,
            designation: designation.to_string(),
        })
    }

    /// Amendment or corrigendum.
    #[must_use]
    pub const fn kind(&self) -> EditionChangeKind {
        self.kind
    }

    /// Exact publisher designation.
    #[must_use]
    pub fn designation(&self) -> &str {
        &self.designation
    }
}

/// Jurisdiction and application profile attached to one source row.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct JurisdictionProfile {
    jurisdiction: String,
    profile: String,
}

impl JurisdictionProfile {
    /// Construct an explicit jurisdiction/profile pair.
    ///
    /// # Errors
    ///
    /// Refuses blank, control-bearing, or overlong values.
    pub fn try_new(jurisdiction: &str, profile: &str) -> Result<Self, ManifestError> {
        validate_compact_text(jurisdiction, "jurisdiction", MAX_STANDARD_STRING_BYTES)?;
        validate_compact_text(profile, "profile", MAX_STANDARD_STRING_BYTES)?;
        Ok(Self {
            jurisdiction: jurisdiction.to_string(),
            profile: profile.to_string(),
        })
    }

    /// Jurisdiction label.
    #[must_use]
    pub fn jurisdiction(&self) -> &str {
        &self.jurisdiction
    }

    /// Application/profile label.
    #[must_use]
    pub fn profile(&self) -> &str {
        &self.profile
    }
}

/// External reference to protected or open source bytes.
///
/// The type intentionally has no `text` or `body` field. `catalog` identifies
/// the publisher/repository and `locator` tells an authorized caller where to
/// obtain bytes out of band.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProtectedTextReference {
    catalog: String,
    locator: String,
}

impl ProtectedTextReference {
    /// Construct an external source reference.
    ///
    /// # Errors
    ///
    /// Refuses blank, control-bearing, or overlong fields.
    pub fn try_new(catalog: &str, locator: &str) -> Result<Self, ManifestError> {
        validate_compact_text(catalog, "source.catalog", MAX_STANDARD_STRING_BYTES)?;
        validate_compact_text(locator, "source.locator", MAX_STANDARD_STRING_BYTES)?;
        Ok(Self {
            catalog: catalog.to_string(),
            locator: locator.to_string(),
        })
    }

    /// Publisher or repository catalog.
    #[must_use]
    pub fn catalog(&self) -> &str {
        &self.catalog
    }

    /// External clause/document locator; never source body text.
    #[must_use]
    pub fn locator(&self) -> &str {
        &self.locator
    }
}

/// Content-addressed source state.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SourcePin {
    /// Exact hash of the complete authorized source bytes.
    Pinned(ContentHash),
    /// Known target whose source bytes have not been authenticated.
    Unpinned,
}

/// License or access-rights policy for source bytes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum StandardLicense {
    /// Openly redistributable bytes under an SPDX expression.
    Spdx {
        /// SPDX expression.
        id: String,
    },
    /// Protected bytes under named publisher terms.
    Restricted {
        /// Human-readable terms identifier, not the protected terms body.
        terms: String,
    },
    /// Bytes supplied locally by an authorized user.
    UserSupplied {
        /// Access-policy identifier acknowledged by the supplying user.
        policy: String,
    },
}

/// Current accessibility of the referenced source.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SourceAccess {
    /// The configured source is currently accessible under its policy.
    Available,
    /// Access was revoked; the historical row remains inspectable.
    Revoked {
        /// Bounded reason code, not secret or protected prose.
        reason_code: String,
    },
}

/// Lifecycle state of an exact edition.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum StandardStatus {
    /// Current for its declared jurisdiction/profile.
    Current,
    /// Withdrawn without a declared successor.
    Withdrawn {
        /// Bounded reason code.
        reason_code: String,
    },
    /// Superseded by another exact key in the same manifest.
    Superseded {
        /// Exact successor; support never transfers by family name alone.
        by: StandardEditionKey,
    },
}

impl StandardStatus {
    const fn tag(&self) -> u8 {
        match self {
            Self::Current => 1,
            Self::Withdrawn { .. } => 2,
            Self::Superseded { .. } => 3,
        }
    }
}

/// Untrusted authoring input for one source-manifest row.
///
/// Admission through [`StandardManifest::try_new`] validates every field,
/// canonicalizes row order, checks collisions, and proves the complete
/// supersession graph before a [`StandardSourceRecord`] can exist.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StandardSourceDraft {
    /// Exact standard/part/edition key.
    pub key: StandardEditionKey,
    /// Human title; no protected normative text.
    pub title: String,
    /// Ordered incorporated amendments/corrigenda.
    pub changes: Vec<EditionChange>,
    /// Current, withdrawn, or exact supersession state.
    pub status: StandardStatus,
    /// Jurisdiction and application profile.
    pub jurisdiction: JurisdictionProfile,
    /// External reference to source bytes.
    pub source: ProtectedTextReference,
    /// Exact source hash or explicit unpinned state.
    pub source_pin: SourcePin,
    /// License/access-rights policy.
    pub license: StandardLicense,
    /// Current access state.
    pub access: SourceAccess,
    /// Explicit no-claim boundary for the row.
    pub no_claim: String,
}

/// Admitted, canonical source-manifest row.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StandardSourceRecord {
    draft: StandardSourceDraft,
    identity: ContentHash,
}

impl StandardSourceRecord {
    /// Exact edition key.
    #[must_use]
    pub const fn key(&self) -> &StandardEditionKey {
        &self.draft.key
    }

    /// Human title.
    #[must_use]
    pub fn title(&self) -> &str {
        &self.draft.title
    }

    /// Ordered amendment/corrigendum chain.
    #[must_use]
    pub fn changes(&self) -> &[EditionChange] {
        &self.draft.changes
    }

    /// Edition lifecycle state.
    #[must_use]
    pub const fn status(&self) -> &StandardStatus {
        &self.draft.status
    }

    /// Jurisdiction and profile.
    #[must_use]
    pub const fn jurisdiction(&self) -> &JurisdictionProfile {
        &self.draft.jurisdiction
    }

    /// Protected-text reference.
    #[must_use]
    pub const fn source(&self) -> &ProtectedTextReference {
        &self.draft.source
    }

    /// Exact source pin state.
    #[must_use]
    pub const fn source_pin(&self) -> SourcePin {
        self.draft.source_pin
    }

    /// License policy.
    #[must_use]
    pub const fn license(&self) -> &StandardLicense {
        &self.draft.license
    }

    /// Current access state.
    #[must_use]
    pub const fn access(&self) -> &SourceAccess {
        &self.draft.access
    }

    /// Explicit no-claim boundary.
    #[must_use]
    pub fn no_claim(&self) -> &str {
        &self.draft.no_claim
    }

    /// Canonical identity of every semantic row field.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
}

/// Why manifest construction or canonical decoding refused.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ManifestError {
    /// A required string is blank.
    EmptyField {
        /// Stable field path.
        field: &'static str,
    },
    /// An identifier is not a lowercase ASCII slug.
    InvalidIdentifier {
        /// Stable field path.
        field: &'static str,
    },
    /// A compact field contains a control character.
    ControlCharacter {
        /// Stable field path.
        field: &'static str,
    },
    /// One string exceeds the explicit byte cap.
    StringLimit {
        /// Stable field path.
        field: &'static str,
        /// Offered byte length.
        len: usize,
        /// Admitted cap.
        limit: usize,
    },
    /// A caller tried to relax a compile-time resource ceiling.
    LimitExceedsHardMaximum {
        /// Limit field.
        field: &'static str,
        /// Offered limit.
        offered: usize,
        /// Compile-time maximum.
        maximum: usize,
    },
    /// Edition row count exceeds the envelope.
    RecordLimit {
        /// Offered row count.
        count: usize,
        /// Admitted cap.
        limit: usize,
    },
    /// Amendment/corrigendum count exceeds the envelope.
    ChangeLimit {
        /// Owning exact edition.
        key: StandardEditionKey,
        /// Offered change count.
        count: usize,
        /// Admitted cap.
        limit: usize,
    },
    /// The canonical byte payload exceeds its envelope.
    ManifestByteLimit {
        /// Required bytes.
        required: usize,
        /// Admitted cap.
        limit: usize,
    },
    /// Checked canonical-size arithmetic overflowed.
    SizeOverflow,
    /// A fallible allocation refused.
    AllocationFailed {
        /// Allocation purpose.
        what: &'static str,
        /// Required element/byte count.
        required: usize,
    },
    /// Two rows declare the same exact edition key.
    EditionCollision {
        /// Colliding key.
        key: StandardEditionKey,
    },
    /// One edition repeats the same typed change designation.
    DuplicateChange {
        /// Owning key.
        key: StandardEditionKey,
        /// Repeated designation.
        designation: String,
    },
    /// A row supersedes itself.
    SelfSupersession {
        /// Invalid key.
        key: StandardEditionKey,
    },
    /// A supersession target is absent.
    UnknownSupersessionTarget {
        /// Historical source key.
        key: StandardEditionKey,
        /// Missing target.
        target: StandardEditionKey,
    },
    /// Supersession edges form a cycle.
    SupersessionCycle {
        /// First source used to expose the cycle.
        key: StandardEditionKey,
    },
    /// A pinned hash is the all-zero sentinel.
    ZeroHash {
        /// Hash role.
        field: &'static str,
    },
    /// Canonical wire magic is wrong or truncated.
    InvalidMagic,
    /// Wire schema is unknown.
    UnknownSchema {
        /// Decoded version.
        version: u32,
    },
    /// Canonical bytes end before a named field completes.
    UnexpectedEof {
        /// Field being decoded.
        field: &'static str,
    },
    /// A wire enum tag is unknown.
    UnknownTag {
        /// Tagged field.
        field: &'static str,
        /// Rejected byte.
        tag: u8,
    },
    /// A framed string is not UTF-8.
    InvalidUtf8 {
        /// Field being decoded.
        field: &'static str,
    },
    /// Bytes remain after the declared manifest.
    TrailingBytes {
        /// Remaining byte count.
        remaining: usize,
    },
    /// The input decodes semantically but is not in canonical sorted encoding.
    NonCanonicalEncoding,
}

impl fmt::Display for ManifestError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyField { field } => {
                write!(formatter, "standard manifest field {field} is blank")
            }
            Self::InvalidIdentifier { field } => write!(
                formatter,
                "standard manifest field {field} must be a lowercase ASCII slug"
            ),
            Self::ControlCharacter { field } => {
                write!(
                    formatter,
                    "standard manifest field {field} contains a control character"
                )
            }
            Self::StringLimit { field, len, limit } => write!(
                formatter,
                "standard manifest field {field} is {len} bytes, exceeding {limit}"
            ),
            Self::LimitExceedsHardMaximum {
                field,
                offered,
                maximum,
            } => write!(
                formatter,
                "standard manifest limit {field}={offered} exceeds hard maximum {maximum}"
            ),
            Self::RecordLimit { count, limit } => write!(
                formatter,
                "standard manifest has {count} rows, exceeding {limit}"
            ),
            Self::ChangeLimit { key, count, limit } => write!(
                formatter,
                "standard edition {}:{} has {count} changes, exceeding {limit}",
                key.standard_id().as_str(),
                key.edition()
            ),
            Self::ManifestByteLimit { required, limit } => write!(
                formatter,
                "standard manifest requires {required} canonical bytes, exceeding {limit}"
            ),
            Self::SizeOverflow => {
                formatter.write_str("standard manifest size arithmetic overflowed")
            }
            Self::AllocationFailed { what, required } => write!(
                formatter,
                "standard manifest could not reserve {required} units for {what}"
            ),
            Self::EditionCollision { key } => write!(
                formatter,
                "standard edition collision for {} part {:?} edition {}",
                key.standard_id().as_str(),
                key.part(),
                key.edition()
            ),
            Self::DuplicateChange { key, designation } => write!(
                formatter,
                "standard edition {}:{} repeats change {designation}",
                key.standard_id().as_str(),
                key.edition()
            ),
            Self::SelfSupersession { key } => write!(
                formatter,
                "standard edition {}:{} supersedes itself",
                key.standard_id().as_str(),
                key.edition()
            ),
            Self::UnknownSupersessionTarget { key, target } => write!(
                formatter,
                "standard edition {}:{} names absent successor {}:{}",
                key.standard_id().as_str(),
                key.edition(),
                target.standard_id().as_str(),
                target.edition()
            ),
            Self::SupersessionCycle { key } => write!(
                formatter,
                "standard supersession cycle is reachable from {}:{}",
                key.standard_id().as_str(),
                key.edition()
            ),
            Self::ZeroHash { field } => {
                write!(
                    formatter,
                    "standard manifest {field} cannot use the zero hash sentinel"
                )
            }
            Self::InvalidMagic => formatter.write_str("standard manifest wire magic is invalid"),
            Self::UnknownSchema { version } => {
                write!(
                    formatter,
                    "standard manifest wire schema {version} is unsupported"
                )
            }
            Self::UnexpectedEof { field } => {
                write!(formatter, "standard manifest ended while decoding {field}")
            }
            Self::UnknownTag { field, tag } => {
                write!(formatter, "standard manifest {field} tag {tag} is unknown")
            }
            Self::InvalidUtf8 { field } => {
                write!(formatter, "standard manifest field {field} is not UTF-8")
            }
            Self::TrailingBytes { remaining } => {
                write!(
                    formatter,
                    "standard manifest has {remaining} trailing bytes"
                )
            }
            Self::NonCanonicalEncoding => formatter
                .write_str("standard manifest bytes are semantically valid but non-canonical"),
        }
    }
}

impl std::error::Error for ManifestError {}

/// Versioned, sorted, content-addressed standards source manifest.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StandardManifest {
    records: Vec<StandardSourceRecord>,
    canonical: Vec<u8>,
    identity: ContentHash,
    limits: ManifestLimits,
}

impl StandardManifest {
    /// Validate drafts, canonicalize row order, verify supersession closure,
    /// and seal the canonical wire identity.
    ///
    /// # Errors
    ///
    /// Returns a deterministic [`ManifestError`] before publishing a partial
    /// manifest.
    pub fn try_new(
        mut drafts: Vec<StandardSourceDraft>,
        limits: ManifestLimits,
    ) -> Result<Self, ManifestError> {
        validate_limits(limits)?;
        if drafts.len() > limits.records {
            return Err(ManifestError::RecordLimit {
                count: drafts.len(),
                limit: limits.records,
            });
        }
        for draft in &drafts {
            validate_draft(draft, limits)?;
        }
        drafts.sort_by(|left, right| left.key.cmp(&right.key));
        for pair in drafts.windows(2) {
            if pair[0].key == pair[1].key {
                return Err(ManifestError::EditionCollision {
                    key: pair[0].key.clone(),
                });
            }
        }
        validate_supersession_graph(&drafts)?;

        let canonical = encode_manifest(&drafts, limits)?;
        let identity = hash_domain(STANDARD_MANIFEST_IDENTITY_DOMAIN, &canonical);
        let mut records = Vec::new();
        records
            .try_reserve_exact(drafts.len())
            .map_err(|_| ManifestError::AllocationFailed {
                what: "source records",
                required: drafts.len(),
            })?;
        for draft in drafts {
            let bytes = encode_record(&draft, limits)?;
            records.push(StandardSourceRecord {
                identity: hash_domain(STANDARD_SOURCE_IDENTITY_DOMAIN, &bytes),
                draft,
            });
        }
        Ok(Self {
            records,
            canonical,
            identity,
            limits,
        })
    }

    /// Decode and re-admit canonical wire bytes under `limits`.
    ///
    /// # Errors
    ///
    /// Refuses malformed, oversized, unknown-version, non-canonical, or
    /// semantically invalid bytes.
    pub fn decode_canonical(bytes: &[u8], limits: ManifestLimits) -> Result<Self, ManifestError> {
        validate_limits(limits)?;
        if bytes.len() > limits.manifest_bytes {
            return Err(ManifestError::ManifestByteLimit {
                required: bytes.len(),
                limit: limits.manifest_bytes,
            });
        }
        let mut reader = WireReader::new(bytes, limits);
        let magic = reader.read_array::<4>("magic")?;
        if magic != WIRE_MAGIC {
            return Err(ManifestError::InvalidMagic);
        }
        let version = reader.read_u32("schema")?;
        if version != STANDARD_MANIFEST_SCHEMA_VERSION {
            return Err(ManifestError::UnknownSchema { version });
        }
        let record_count = usize::try_from(reader.read_u32("record_count")?)
            .map_err(|_| ManifestError::SizeOverflow)?;
        if record_count > limits.records {
            return Err(ManifestError::RecordLimit {
                count: record_count,
                limit: limits.records,
            });
        }
        let mut drafts = Vec::new();
        drafts
            .try_reserve_exact(record_count)
            .map_err(|_| ManifestError::AllocationFailed {
                what: "decoded drafts",
                required: record_count,
            })?;
        for _ in 0..record_count {
            drafts.push(reader.read_draft()?);
        }
        if reader.remaining() != 0 {
            return Err(ManifestError::TrailingBytes {
                remaining: reader.remaining(),
            });
        }
        let manifest = Self::try_new(drafts, limits)?;
        if manifest.canonical.as_slice() != bytes {
            return Err(ManifestError::NonCanonicalEncoding);
        }
        Ok(manifest)
    }

    /// Canonically sorted source rows.
    #[must_use]
    pub fn records(&self) -> &[StandardSourceRecord] {
        &self.records
    }

    /// Exact-key lookup; an unknown or wrong edition never falls back.
    #[must_use]
    pub fn record(&self, key: &StandardEditionKey) -> Option<&StandardSourceRecord> {
        self.records
            .binary_search_by(|record| record.key().cmp(key))
            .ok()
            .map(|index| &self.records[index])
    }

    /// Canonical schema bytes retained at admission.
    #[must_use]
    pub fn canonical_bytes(&self) -> &[u8] {
        &self.canonical
    }

    /// Domain-separated identity of the complete manifest.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Envelope used to admit this manifest.
    #[must_use]
    pub const fn limits(&self) -> ManifestLimits {
        self.limits
    }

    /// Emit deterministic metadata-only rows. The output retains exact keys,
    /// locators, access state, and identities but never source body text or
    /// restricted license prose.
    #[must_use]
    pub fn redacted_rows(&self) -> Vec<String> {
        self.records.iter().map(redacted_row).collect()
    }

    /// Bind one derived rule to an exact admitted source row.
    ///
    /// # Errors
    ///
    /// Refuses unknown/wrong editions, unpinned or revoked sources, hash
    /// mismatch, stale editions under current-only use, unread sources, or an
    /// invalid/zero derived identity.
    pub fn bind_rule(
        &self,
        request: RuleBindingRequest,
    ) -> Result<RuleProvenance, RuleBindingError> {
        validate_identifier(&request.rule_id, "rule_id", self.limits.string_bytes)
            .map_err(RuleBindingError::InvalidRequest)?;
        validate_compact_text(request.clause.as_str(), "clause", self.limits.string_bytes)
            .map_err(RuleBindingError::InvalidRequest)?;
        if request.derived_rule_hash.as_bytes() == &[0; 32] {
            return Err(RuleBindingError::InvalidRequest(ManifestError::ZeroHash {
                field: "derived_rule_hash",
            }));
        }
        if request.reference_state == ReferenceState::Unread {
            return Err(RuleBindingError::UnreadSource);
        }
        let record =
            self.record(&request.source)
                .ok_or_else(|| RuleBindingError::UnknownEdition {
                    key: request.source.clone(),
                })?;
        let SourcePin::Pinned(expected_source_hash) = record.source_pin() else {
            return Err(RuleBindingError::UnpinnedSource {
                key: request.source,
            });
        };
        if matches!(record.access(), SourceAccess::Revoked { .. }) {
            return Err(RuleBindingError::AccessRevoked {
                key: request.source,
            });
        }
        if expected_source_hash != request.observed_source_hash {
            return Err(RuleBindingError::SourceHashMismatch {
                key: request.source,
                expected: expected_source_hash,
                observed: request.observed_source_hash,
            });
        }
        let historical = !matches!(record.status(), StandardStatus::Current);
        if historical && request.use_policy == RuleUsePolicy::CurrentOnly {
            return Err(RuleBindingError::HistoricalEdition {
                key: request.source,
            });
        }

        let source_identity = record.identity();
        let identity = rule_identity(
            self.identity,
            source_identity,
            expected_source_hash,
            &request,
            historical,
        )
        .map_err(RuleBindingError::InvalidRequest)?;
        Ok(RuleProvenance {
            source: request.source,
            clause: request.clause,
            rule_id: request.rule_id,
            source_hash: expected_source_hash,
            derived_rule_hash: request.derived_rule_hash,
            reference_state: request.reference_state,
            use_policy: request.use_policy,
            historical,
            source_identity,
            manifest_identity: self.identity,
            identity,
        })
    }
}

/// Exact clause/subclause/table/figure locator without source body text.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ClauseLocator(String);

impl ClauseLocator {
    /// Construct a compact exact locator.
    ///
    /// # Errors
    ///
    /// Refuses blank, control-bearing, or overlong locators.
    pub fn try_new(value: &str) -> Result<Self, ManifestError> {
        validate_compact_text(value, "clause", MAX_STANDARD_STRING_BYTES)?;
        Ok(Self(value.to_string()))
    }

    /// Borrow the locator.
    #[must_use]
    pub fn as_str(&self) -> &str {
        &self.0
    }
}

/// Explicit source-engagement state for one derived rule.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReferenceState {
    /// Source has not been read; rule binding refuses.
    Unread,
    /// Source was read at the named locator.
    Read,
    /// Rule was derived from the source.
    Derived,
    /// Derivation was reproduced independently.
    Reproduced,
}

impl ReferenceState {
    const fn tag(self) -> u8 {
        match self {
            Self::Unread => 0,
            Self::Read => 1,
            Self::Derived => 2,
            Self::Reproduced => 3,
        }
    }
}

/// Whether a rule may bind historical editions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RuleUsePolicy {
    /// Only a row explicitly marked current may bind.
    CurrentOnly,
    /// Withdrawn/superseded rows may bind for historical replay only.
    HistoricalReplay,
}

impl RuleUsePolicy {
    const fn tag(self) -> u8 {
        match self {
            Self::CurrentOnly => 1,
            Self::HistoricalReplay => 2,
        }
    }
}

/// Untrusted request to bind one derived rule.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RuleBindingRequest {
    /// Stable lowercase rule identifier.
    pub rule_id: String,
    /// Exact source edition.
    pub source: StandardEditionKey,
    /// Exact clause locator.
    pub clause: ClauseLocator,
    /// Hash observed from the authorized source bytes.
    pub observed_source_hash: ContentHash,
    /// Hash of canonical derived-rule bytes, never protected source bytes.
    pub derived_rule_hash: ContentHash,
    /// Read/derived/reproduced state.
    pub reference_state: ReferenceState,
    /// Current-use versus historical-replay policy.
    pub use_policy: RuleUsePolicy,
}

/// Sealed exact provenance for one derived rule.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RuleProvenance {
    source: StandardEditionKey,
    clause: ClauseLocator,
    rule_id: String,
    source_hash: ContentHash,
    derived_rule_hash: ContentHash,
    reference_state: ReferenceState,
    use_policy: RuleUsePolicy,
    historical: bool,
    source_identity: ContentHash,
    manifest_identity: ContentHash,
    identity: ContentHash,
}

impl RuleProvenance {
    /// Exact source edition.
    #[must_use]
    pub const fn source(&self) -> &StandardEditionKey {
        &self.source
    }

    /// Exact source locator.
    #[must_use]
    pub const fn clause(&self) -> &ClauseLocator {
        &self.clause
    }

    /// Stable derived-rule id.
    #[must_use]
    pub fn rule_id(&self) -> &str {
        &self.rule_id
    }

    /// Exact authenticated source hash.
    #[must_use]
    pub const fn source_hash(&self) -> ContentHash {
        self.source_hash
    }

    /// Canonical derived-rule hash.
    #[must_use]
    pub const fn derived_rule_hash(&self) -> ContentHash {
        self.derived_rule_hash
    }

    /// Explicit engagement state.
    #[must_use]
    pub const fn reference_state(&self) -> ReferenceState {
        self.reference_state
    }

    /// Admission policy attached to this rule use.
    #[must_use]
    pub const fn use_policy(&self) -> RuleUsePolicy {
        self.use_policy
    }

    /// Whether this provenance is historical-only.
    #[must_use]
    pub const fn historical(&self) -> bool {
        self.historical
    }

    /// Canonical source-row identity.
    #[must_use]
    pub const fn source_identity(&self) -> ContentHash {
        self.source_identity
    }

    /// Exact manifest identity used for admission.
    #[must_use]
    pub const fn manifest_identity(&self) -> ContentHash {
        self.manifest_identity
    }

    /// Domain-separated identity of the complete provenance record.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
}

/// Why derived-rule binding refused.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RuleBindingError {
    /// Request validation failed.
    InvalidRequest(ManifestError),
    /// Exact edition is absent; no family/part fallback is attempted.
    UnknownEdition {
        /// Missing key.
        key: StandardEditionKey,
    },
    /// Exact source bytes are not pinned.
    UnpinnedSource {
        /// Refused key.
        key: StandardEditionKey,
    },
    /// Source access was revoked.
    AccessRevoked {
        /// Refused key.
        key: StandardEditionKey,
    },
    /// Observed source bytes do not match the manifest pin.
    SourceHashMismatch {
        /// Refused key.
        key: StandardEditionKey,
        /// Manifest pin.
        expected: ContentHash,
        /// Observed source hash.
        observed: ContentHash,
    },
    /// Historical edition was offered to a current-only rule.
    HistoricalEdition {
        /// Refused key.
        key: StandardEditionKey,
    },
    /// An unread source cannot support derived-rule provenance.
    UnreadSource,
}

impl fmt::Display for RuleBindingError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidRequest(error) => {
                write!(formatter, "rule binding request invalid: {error}")
            }
            Self::UnknownEdition { key } => write!(
                formatter,
                "unknown exact standard edition {} part {:?} edition {}",
                key.standard_id().as_str(),
                key.part(),
                key.edition()
            ),
            Self::UnpinnedSource { key } => write!(
                formatter,
                "standard edition {}:{} has no authenticated source hash",
                key.standard_id().as_str(),
                key.edition()
            ),
            Self::AccessRevoked { key } => write!(
                formatter,
                "access to standard edition {}:{} is revoked",
                key.standard_id().as_str(),
                key.edition()
            ),
            Self::SourceHashMismatch { key, .. } => write!(
                formatter,
                "observed bytes do not match standard edition {}:{}",
                key.standard_id().as_str(),
                key.edition()
            ),
            Self::HistoricalEdition { key } => write!(
                formatter,
                "standard edition {}:{} is historical and cannot support current rules",
                key.standard_id().as_str(),
                key.edition()
            ),
            Self::UnreadSource => formatter.write_str(
                "a derived rule cannot bind an unread source; declare read, derived, or reproduced",
            ),
        }
    }
}

impl std::error::Error for RuleBindingError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::InvalidRequest(error) => Some(error),
            _ => None,
        }
    }
}

fn validate_identifier(
    value: &str,
    field: &'static str,
    max_bytes: usize,
) -> Result<(), ManifestError> {
    validate_string_limit(value, field, max_bytes)?;
    let bytes = value.as_bytes();
    let Some((&first, rest)) = bytes.split_first() else {
        return Err(ManifestError::EmptyField { field });
    };
    if !first.is_ascii_lowercase()
        || !bytes
            .last()
            .is_some_and(|byte| byte.is_ascii_alphanumeric())
        || !rest
            .iter()
            .all(|byte| byte.is_ascii_lowercase() || byte.is_ascii_digit() || *byte == b'-')
        || bytes.windows(2).any(|pair| pair == b"--")
    {
        return Err(ManifestError::InvalidIdentifier { field });
    }
    Ok(())
}

fn validate_limits(limits: ManifestLimits) -> Result<(), ManifestError> {
    for (field, offered, maximum) in [
        ("records", limits.records, MAX_STANDARD_RECORDS),
        (
            "changes_per_record",
            limits.changes_per_record,
            MAX_EDITION_CHANGES,
        ),
        (
            "string_bytes",
            limits.string_bytes,
            MAX_STANDARD_STRING_BYTES,
        ),
        (
            "manifest_bytes",
            limits.manifest_bytes,
            MAX_STANDARD_MANIFEST_BYTES,
        ),
    ] {
        if offered > maximum {
            return Err(ManifestError::LimitExceedsHardMaximum {
                field,
                offered,
                maximum,
            });
        }
    }
    Ok(())
}

fn validate_compact_text(
    value: &str,
    field: &'static str,
    max_bytes: usize,
) -> Result<(), ManifestError> {
    validate_string_limit(value, field, max_bytes)?;
    if value.trim().is_empty() {
        return Err(ManifestError::EmptyField { field });
    }
    if value.chars().any(char::is_control) {
        return Err(ManifestError::ControlCharacter { field });
    }
    Ok(())
}

fn validate_string_limit(
    value: &str,
    field: &'static str,
    max_bytes: usize,
) -> Result<(), ManifestError> {
    if value.len() > max_bytes {
        return Err(ManifestError::StringLimit {
            field,
            len: value.len(),
            limit: max_bytes,
        });
    }
    Ok(())
}

fn validate_key(key: &StandardEditionKey, limits: ManifestLimits) -> Result<(), ManifestError> {
    validate_identifier(key.standard_id.as_str(), "standard_id", limits.string_bytes)?;
    if let Some(part) = key.part() {
        validate_compact_text(part, "part", limits.string_bytes)?;
    }
    validate_compact_text(key.edition(), "edition", limits.string_bytes)
}

fn validate_draft(
    draft: &StandardSourceDraft,
    limits: ManifestLimits,
) -> Result<(), ManifestError> {
    validate_key(&draft.key, limits)?;
    validate_compact_text(&draft.title, "title", limits.string_bytes)?;
    if draft.changes.len() > limits.changes_per_record {
        return Err(ManifestError::ChangeLimit {
            key: draft.key.clone(),
            count: draft.changes.len(),
            limit: limits.changes_per_record,
        });
    }
    for (index, change) in draft.changes.iter().enumerate() {
        validate_compact_text(
            change.designation(),
            "change.designation",
            limits.string_bytes,
        )?;
        if draft.changes[..index].iter().any(|prior| {
            prior.kind() == change.kind() && prior.designation() == change.designation()
        }) {
            return Err(ManifestError::DuplicateChange {
                key: draft.key.clone(),
                designation: change.designation().to_string(),
            });
        }
    }
    match &draft.status {
        StandardStatus::Current => {}
        StandardStatus::Withdrawn { reason_code } => {
            validate_identifier(reason_code, "status.reason_code", limits.string_bytes)?;
        }
        StandardStatus::Superseded { by } => validate_key(by, limits)?,
    }
    validate_compact_text(
        draft.jurisdiction.jurisdiction(),
        "jurisdiction",
        limits.string_bytes,
    )?;
    validate_compact_text(draft.jurisdiction.profile(), "profile", limits.string_bytes)?;
    validate_compact_text(
        draft.source.catalog(),
        "source.catalog",
        limits.string_bytes,
    )?;
    validate_compact_text(
        draft.source.locator(),
        "source.locator",
        limits.string_bytes,
    )?;
    if let SourcePin::Pinned(hash) = draft.source_pin
        && hash.as_bytes() == &[0; 32]
    {
        return Err(ManifestError::ZeroHash {
            field: "source_hash",
        });
    }
    match &draft.license {
        StandardLicense::Spdx { id } => {
            validate_compact_text(id, "license.spdx", limits.string_bytes)?;
        }
        StandardLicense::Restricted { terms } => {
            validate_compact_text(terms, "license.terms", limits.string_bytes)?;
        }
        StandardLicense::UserSupplied { policy } => {
            validate_compact_text(policy, "license.policy", limits.string_bytes)?;
        }
    }
    if let SourceAccess::Revoked { reason_code } = &draft.access {
        validate_identifier(reason_code, "access.reason_code", limits.string_bytes)?;
    }
    validate_compact_text(&draft.no_claim, "no_claim", limits.string_bytes)
}

fn validate_supersession_graph(drafts: &[StandardSourceDraft]) -> Result<(), ManifestError> {
    for draft in drafts {
        let StandardStatus::Superseded { by } = &draft.status else {
            continue;
        };
        if by == &draft.key {
            return Err(ManifestError::SelfSupersession {
                key: draft.key.clone(),
            });
        }
        if drafts.binary_search_by(|row| row.key.cmp(by)).is_err() {
            return Err(ManifestError::UnknownSupersessionTarget {
                key: draft.key.clone(),
                target: by.clone(),
            });
        }
    }
    // White/gray/black traversal over the functional supersession graph.
    // Exact targets were checked above, so each edge lookup is deterministic.
    let mut state = Vec::new();
    state
        .try_reserve_exact(drafts.len())
        .map_err(|_| ManifestError::AllocationFailed {
            what: "supersession traversal state",
            required: drafts.len(),
        })?;
    state.resize(drafts.len(), 0_u8);
    for start in 0..drafts.len() {
        if state[start] != 0 {
            continue;
        }
        let mut cursor = start;
        while state[cursor] == 0 {
            state[cursor] = 1;
            let StandardStatus::Superseded { by } = &drafts[cursor].status else {
                break;
            };
            cursor = drafts
                .binary_search_by(|row| row.key.cmp(by))
                .map_err(|_| ManifestError::UnknownSupersessionTarget {
                    key: drafts[start].key.clone(),
                    target: by.clone(),
                })?;
        }
        if state[cursor] == 1 && matches!(&drafts[cursor].status, StandardStatus::Superseded { .. })
        {
            return Err(ManifestError::SupersessionCycle {
                key: drafts[start].key.clone(),
            });
        }

        cursor = start;
        while state[cursor] == 1 {
            state[cursor] = 2;
            let StandardStatus::Superseded { by } = &drafts[cursor].status else {
                break;
            };
            cursor = drafts
                .binary_search_by(|row| row.key.cmp(by))
                .map_err(|_| ManifestError::UnknownSupersessionTarget {
                    key: drafts[start].key.clone(),
                    target: by.clone(),
                })?;
        }
    }
    Ok(())
}

fn checked_add(left: usize, right: usize) -> Result<usize, ManifestError> {
    left.checked_add(right).ok_or(ManifestError::SizeOverflow)
}

fn string_wire_len(value: &str) -> Result<usize, ManifestError> {
    checked_add(4, value.len())
}

fn key_wire_len(key: &StandardEditionKey) -> Result<usize, ManifestError> {
    let mut len = string_wire_len(key.standard_id().as_str())?;
    len = checked_add(len, 1)?;
    if let Some(part) = key.part() {
        len = checked_add(len, string_wire_len(part)?)?;
    }
    checked_add(len, string_wire_len(key.edition())?)
}

#[allow(clippy::too_many_lines)] // Mirrors every canonical field exactly once.
fn record_wire_len(draft: &StandardSourceDraft) -> Result<usize, ManifestError> {
    let mut len = key_wire_len(&draft.key)?;
    len = checked_add(len, string_wire_len(&draft.title)?)?;
    len = checked_add(len, 4)?;
    for change in &draft.changes {
        len = checked_add(len, 1)?;
        len = checked_add(len, string_wire_len(change.designation())?)?;
    }
    len = checked_add(len, 1)?;
    match &draft.status {
        StandardStatus::Current => {}
        StandardStatus::Withdrawn { reason_code } => {
            len = checked_add(len, string_wire_len(reason_code)?)?;
        }
        StandardStatus::Superseded { by } => {
            len = checked_add(len, key_wire_len(by)?)?;
        }
    }
    len = checked_add(len, string_wire_len(draft.jurisdiction.jurisdiction())?)?;
    len = checked_add(len, string_wire_len(draft.jurisdiction.profile())?)?;
    len = checked_add(len, string_wire_len(draft.source.catalog())?)?;
    len = checked_add(len, string_wire_len(draft.source.locator())?)?;
    len = checked_add(len, 1)?;
    if matches!(draft.source_pin, SourcePin::Pinned(_)) {
        len = checked_add(len, 32)?;
    }
    len = checked_add(len, 1)?;
    match &draft.license {
        StandardLicense::Spdx { id } => len = checked_add(len, string_wire_len(id)?)?,
        StandardLicense::Restricted { terms } => {
            len = checked_add(len, string_wire_len(terms)?)?;
        }
        StandardLicense::UserSupplied { policy } => {
            len = checked_add(len, string_wire_len(policy)?)?;
        }
    }
    len = checked_add(len, 1)?;
    if let SourceAccess::Revoked { reason_code } = &draft.access {
        len = checked_add(len, string_wire_len(reason_code)?)?;
    }
    checked_add(len, string_wire_len(&draft.no_claim)?)
}

fn manifest_wire_len(drafts: &[StandardSourceDraft]) -> Result<usize, ManifestError> {
    let mut len = 12usize;
    for draft in drafts {
        len = checked_add(len, record_wire_len(draft)?)?;
    }
    Ok(len)
}

fn push_u32(out: &mut Vec<u8>, value: usize) -> Result<(), ManifestError> {
    let value = u32::try_from(value).map_err(|_| ManifestError::SizeOverflow)?;
    out.extend_from_slice(&value.to_le_bytes());
    Ok(())
}

fn push_string(out: &mut Vec<u8>, value: &str) -> Result<(), ManifestError> {
    push_u32(out, value.len())?;
    out.extend_from_slice(value.as_bytes());
    Ok(())
}

fn push_key(out: &mut Vec<u8>, key: &StandardEditionKey) -> Result<(), ManifestError> {
    push_string(out, key.standard_id().as_str())?;
    match key.part() {
        Some(part) => {
            out.push(1);
            push_string(out, part)?;
        }
        None => out.push(0),
    }
    push_string(out, key.edition())
}

#[allow(clippy::too_many_lines)] // Fixed field order is the wire contract.
fn push_record(out: &mut Vec<u8>, draft: &StandardSourceDraft) -> Result<(), ManifestError> {
    push_key(out, &draft.key)?;
    push_string(out, &draft.title)?;
    push_u32(out, draft.changes.len())?;
    for change in &draft.changes {
        out.push(change.kind().tag());
        push_string(out, change.designation())?;
    }
    out.push(draft.status.tag());
    match &draft.status {
        StandardStatus::Current => {}
        StandardStatus::Withdrawn { reason_code } => push_string(out, reason_code)?,
        StandardStatus::Superseded { by } => push_key(out, by)?,
    }
    push_string(out, draft.jurisdiction.jurisdiction())?;
    push_string(out, draft.jurisdiction.profile())?;
    push_string(out, draft.source.catalog())?;
    push_string(out, draft.source.locator())?;
    match draft.source_pin {
        SourcePin::Unpinned => out.push(0),
        SourcePin::Pinned(hash) => {
            out.push(1);
            out.extend_from_slice(hash.as_bytes());
        }
    }
    match &draft.license {
        StandardLicense::Spdx { id } => {
            out.push(1);
            push_string(out, id)?;
        }
        StandardLicense::Restricted { terms } => {
            out.push(2);
            push_string(out, terms)?;
        }
        StandardLicense::UserSupplied { policy } => {
            out.push(3);
            push_string(out, policy)?;
        }
    }
    match &draft.access {
        SourceAccess::Available => out.push(1),
        SourceAccess::Revoked { reason_code } => {
            out.push(2);
            push_string(out, reason_code)?;
        }
    }
    push_string(out, &draft.no_claim)
}

fn encode_manifest(
    drafts: &[StandardSourceDraft],
    limits: ManifestLimits,
) -> Result<Vec<u8>, ManifestError> {
    let required = manifest_wire_len(drafts)?;
    if required > limits.manifest_bytes {
        return Err(ManifestError::ManifestByteLimit {
            required,
            limit: limits.manifest_bytes,
        });
    }
    let mut out = Vec::new();
    out.try_reserve_exact(required)
        .map_err(|_| ManifestError::AllocationFailed {
            what: "canonical manifest bytes",
            required,
        })?;
    out.extend_from_slice(&WIRE_MAGIC);
    out.extend_from_slice(&STANDARD_MANIFEST_SCHEMA_VERSION.to_le_bytes());
    push_u32(&mut out, drafts.len())?;
    for draft in drafts {
        push_record(&mut out, draft)?;
    }
    debug_assert_eq!(out.len(), required);
    Ok(out)
}

fn encode_record(
    draft: &StandardSourceDraft,
    limits: ManifestLimits,
) -> Result<Vec<u8>, ManifestError> {
    let required = record_wire_len(draft)?;
    if required > limits.manifest_bytes {
        return Err(ManifestError::ManifestByteLimit {
            required,
            limit: limits.manifest_bytes,
        });
    }
    let mut out = Vec::new();
    out.try_reserve_exact(required)
        .map_err(|_| ManifestError::AllocationFailed {
            what: "canonical source row",
            required,
        })?;
    push_record(&mut out, draft)?;
    debug_assert_eq!(out.len(), required);
    Ok(out)
}

struct WireReader<'a> {
    bytes: &'a [u8],
    offset: usize,
    limits: ManifestLimits,
}

impl<'a> WireReader<'a> {
    const fn new(bytes: &'a [u8], limits: ManifestLimits) -> Self {
        Self {
            bytes,
            offset: 0,
            limits,
        }
    }

    fn remaining(&self) -> usize {
        self.bytes.len() - self.offset
    }

    fn read_array<const N: usize>(
        &mut self,
        field: &'static str,
    ) -> Result<[u8; N], ManifestError> {
        let end = self
            .offset
            .checked_add(N)
            .ok_or(ManifestError::SizeOverflow)?;
        let bytes = self
            .bytes
            .get(self.offset..end)
            .ok_or(ManifestError::UnexpectedEof { field })?;
        self.offset = end;
        bytes
            .try_into()
            .map_err(|_| ManifestError::UnexpectedEof { field })
    }

    fn read_u8(&mut self, field: &'static str) -> Result<u8, ManifestError> {
        Ok(self.read_array::<1>(field)?[0])
    }

    fn read_u32(&mut self, field: &'static str) -> Result<u32, ManifestError> {
        Ok(u32::from_le_bytes(self.read_array::<4>(field)?))
    }

    fn read_hash(&mut self, field: &'static str) -> Result<ContentHash, ManifestError> {
        Ok(ContentHash(self.read_array::<32>(field)?))
    }

    fn read_string(&mut self, field: &'static str) -> Result<String, ManifestError> {
        let len =
            usize::try_from(self.read_u32(field)?).map_err(|_| ManifestError::SizeOverflow)?;
        if len > self.limits.string_bytes {
            return Err(ManifestError::StringLimit {
                field,
                len,
                limit: self.limits.string_bytes,
            });
        }
        let end = self
            .offset
            .checked_add(len)
            .ok_or(ManifestError::SizeOverflow)?;
        let bytes = self
            .bytes
            .get(self.offset..end)
            .ok_or(ManifestError::UnexpectedEof { field })?;
        let value =
            core::str::from_utf8(bytes).map_err(|_| ManifestError::InvalidUtf8 { field })?;
        let mut owned = String::new();
        owned
            .try_reserve_exact(len)
            .map_err(|_| ManifestError::AllocationFailed {
                what: field,
                required: len,
            })?;
        owned.push_str(value);
        self.offset = end;
        Ok(owned)
    }

    fn read_key(&mut self) -> Result<StandardEditionKey, ManifestError> {
        let standard_id = self.read_string("standard_id")?;
        let part = match self.read_u8("part.present")? {
            0 => None,
            1 => Some(self.read_string("part")?),
            tag => {
                return Err(ManifestError::UnknownTag {
                    field: "part.present",
                    tag,
                });
            }
        };
        let edition = self.read_string("edition")?;
        StandardEditionKey::try_new(&standard_id, part.as_deref(), &edition)
    }

    #[allow(clippy::too_many_lines)] // Mirrors the canonical record grammar.
    fn read_draft(&mut self) -> Result<StandardSourceDraft, ManifestError> {
        let key = self.read_key()?;
        let title = self.read_string("title")?;
        let change_count = usize::try_from(self.read_u32("change_count")?)
            .map_err(|_| ManifestError::SizeOverflow)?;
        if change_count > self.limits.changes_per_record {
            return Err(ManifestError::ChangeLimit {
                key,
                count: change_count,
                limit: self.limits.changes_per_record,
            });
        }
        let mut changes = Vec::new();
        changes
            .try_reserve_exact(change_count)
            .map_err(|_| ManifestError::AllocationFailed {
                what: "decoded edition changes",
                required: change_count,
            })?;
        for _ in 0..change_count {
            let kind = match self.read_u8("change.kind")? {
                1 => EditionChangeKind::Amendment,
                2 => EditionChangeKind::Corrigendum,
                tag => {
                    return Err(ManifestError::UnknownTag {
                        field: "change.kind",
                        tag,
                    });
                }
            };
            let designation = self.read_string("change.designation")?;
            changes.push(EditionChange::try_new(kind, &designation)?);
        }
        let status = match self.read_u8("status")? {
            1 => StandardStatus::Current,
            2 => StandardStatus::Withdrawn {
                reason_code: self.read_string("status.reason_code")?,
            },
            3 => StandardStatus::Superseded {
                by: self.read_key()?,
            },
            tag => {
                return Err(ManifestError::UnknownTag {
                    field: "status",
                    tag,
                });
            }
        };
        let jurisdiction = self.read_string("jurisdiction")?;
        let profile = self.read_string("profile")?;
        let catalog = self.read_string("source.catalog")?;
        let locator = self.read_string("source.locator")?;
        let source_pin = match self.read_u8("source_pin")? {
            0 => SourcePin::Unpinned,
            1 => SourcePin::Pinned(self.read_hash("source_hash")?),
            tag => {
                return Err(ManifestError::UnknownTag {
                    field: "source_pin",
                    tag,
                });
            }
        };
        let license = match self.read_u8("license")? {
            1 => StandardLicense::Spdx {
                id: self.read_string("license.spdx")?,
            },
            2 => StandardLicense::Restricted {
                terms: self.read_string("license.terms")?,
            },
            3 => StandardLicense::UserSupplied {
                policy: self.read_string("license.policy")?,
            },
            tag => {
                return Err(ManifestError::UnknownTag {
                    field: "license",
                    tag,
                });
            }
        };
        let access = match self.read_u8("access")? {
            1 => SourceAccess::Available,
            2 => SourceAccess::Revoked {
                reason_code: self.read_string("access.reason_code")?,
            },
            tag => {
                return Err(ManifestError::UnknownTag {
                    field: "access",
                    tag,
                });
            }
        };
        let no_claim = self.read_string("no_claim")?;
        Ok(StandardSourceDraft {
            key,
            title,
            changes,
            status,
            jurisdiction: JurisdictionProfile::try_new(&jurisdiction, &profile)?,
            source: ProtectedTextReference::try_new(&catalog, &locator)?,
            source_pin,
            license,
            access,
            no_claim,
        })
    }
}

fn rule_identity(
    manifest_identity: ContentHash,
    source_identity: ContentHash,
    source_hash: ContentHash,
    request: &RuleBindingRequest,
    historical: bool,
) -> Result<ContentHash, ManifestError> {
    let mut required = checked_add(4, 32 * 3)?;
    required = checked_add(required, key_wire_len(&request.source)?)?;
    required = checked_add(required, string_wire_len(request.clause.as_str())?)?;
    required = checked_add(required, string_wire_len(&request.rule_id)?)?;
    required = checked_add(required, 32 + 3)?;
    let mut bytes = Vec::new();
    bytes
        .try_reserve_exact(required)
        .map_err(|_| ManifestError::AllocationFailed {
            what: "derived-rule identity bytes",
            required,
        })?;
    bytes.extend_from_slice(&STANDARD_MANIFEST_SCHEMA_VERSION.to_le_bytes());
    bytes.extend_from_slice(manifest_identity.as_bytes());
    bytes.extend_from_slice(source_identity.as_bytes());
    bytes.extend_from_slice(source_hash.as_bytes());
    push_key(&mut bytes, &request.source)?;
    push_string(&mut bytes, request.clause.as_str())?;
    push_string(&mut bytes, &request.rule_id)?;
    bytes.extend_from_slice(request.derived_rule_hash.as_bytes());
    bytes.push(request.reference_state.tag());
    bytes.push(request.use_policy.tag());
    bytes.push(u8::from(historical));
    debug_assert_eq!(bytes.len(), required);
    Ok(hash_domain(STANDARD_RULE_IDENTITY_DOMAIN, &bytes))
}

fn redacted_row(record: &StandardSourceRecord) -> String {
    use fmt::Write as _;

    let mut row = String::new();
    row.push_str("{\"standard_id\":");
    crate::push_json_str(&mut row, record.key().standard_id().as_str());
    row.push_str(",\"part\":");
    match record.key().part() {
        Some(part) => crate::push_json_str(&mut row, part),
        None => row.push_str("null"),
    }
    row.push_str(",\"edition\":");
    crate::push_json_str(&mut row, record.key().edition());
    row.push_str(",\"status\":");
    let _ = write!(row, "{}", record.status().tag());
    row.push_str(",\"source_record_identity\":");
    crate::push_json_str(&mut row, &record.identity().to_hex());
    row.push_str(",\"source_hash\":");
    match record.source_pin() {
        SourcePin::Pinned(hash) => crate::push_json_str(&mut row, &hash.to_hex()),
        SourcePin::Unpinned => row.push_str("null"),
    }
    row.push_str(",\"catalog\":");
    crate::push_json_str(&mut row, record.source().catalog());
    row.push_str(",\"locator\":");
    crate::push_json_str(&mut row, record.source().locator());
    row.push_str(",\"access\":");
    crate::push_json_str(
        &mut row,
        match record.access() {
            SourceAccess::Available => "available",
            SourceAccess::Revoked { .. } => "revoked",
        },
    );
    row.push('}');
    row
}
