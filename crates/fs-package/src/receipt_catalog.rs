//! Content-addressed metadata for versioned receipt families (bead h61n).
//!
//! This catalog answers one narrow question: which exact receipt codec did a
//! package or ledger row declare? It binds a globally qualified family id,
//! wire schema, owner identity version/domain, transport profile, and an
//! owner-supplied codec fingerprint. Lookup is exact and never falls forward
//! or backward across versions.
//!
//! Catalog membership is not decoder authority. This module neither imports
//! receipt owners nor interprets payload bytes. A higher L6 adapter that
//! depends on both this crate and the owner must still dispatch the exact
//! decoder, verify the owner receipt identity, and replay the receipt against
//! known semantics before making any scientific or recovery claim.

use std::collections::BTreeMap;

use fs_blake3::{ContentHash, hash_domain};

/// Wire schema understood by this catalog implementation.
pub const RECEIPT_SCHEMA_CATALOG_VERSION: u32 = 1;
/// Semantic identity version of one receipt-schema descriptor.
pub const RECEIPT_SCHEMA_DESCRIPTOR_IDENTITY_VERSION: u32 = 1;
/// Domain-separated identity of one receipt-schema descriptor.
pub const RECEIPT_SCHEMA_DESCRIPTOR_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-package.receipt-schema-descriptor.v1";
/// Semantic identity version of a complete canonical catalog.
pub const RECEIPT_SCHEMA_CATALOG_IDENTITY_VERSION: u32 = 1;
/// Domain-separated identity of a complete canonical catalog.
pub const RECEIPT_SCHEMA_CATALOG_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-package.receipt-schema-catalog.v1";

/// Maximum UTF-8 bytes in a globally qualified receipt-family id.
pub const MAX_RECEIPT_FAMILY_ID_BYTES: usize = 128;
/// Maximum UTF-8 bytes in an owner receipt identity domain.
pub const MAX_RECEIPT_IDENTITY_DOMAIN_BYTES: usize = 256;
/// Maximum canonical receipt bytes a replayable descriptor may declare.
pub const MAX_RECEIPT_TRANSPORT_BYTES: u64 = 64 * 1024 * 1024;
/// Maximum receipt-family versions in one catalog.
pub const MAX_RECEIPT_SCHEMA_ENTRIES: usize = 4_096;
/// Maximum canonical catalog transport bytes, including its in-band identity.
pub const MAX_RECEIPT_SCHEMA_CATALOG_BYTES: usize = 4 * 1024 * 1024;

const CATALOG_MAGIC: &[u8; 8] = b"FSPRCAT\0";
const FIELD_FAMILY: u8 = 1;
const FIELD_WIRE_SCHEMA_VERSION: u8 = 2;
const FIELD_OWNER_IDENTITY_VERSION: u8 = 3;
const FIELD_OWNER_IDENTITY_DOMAIN: u8 = 4;
const FIELD_TRANSPORT: u8 = 5;
const FIELD_OWNER_SCHEMA_FINGERPRINT: u8 = 6;
const FIELD_DESCRIPTOR_HASH: u8 = 7;
const TRANSPORT_DIGEST_ONLY: u8 = 1;
const TRANSPORT_CANONICAL_BYTES: u8 = 2;
// Tags, fixed widths, and the shortest accepted `a:b` / `a.b` identities.
const MIN_ENCODED_RECEIPT_SCHEMA_DESCRIPTOR_BYTES: usize =
    (1 + 8 + 3) + (1 + 4) + (1 + 4) + (1 + 8 + 3) + (1 + 1 + 8) + (1 + 32) + (1 + 32);

/// How a retained receipt may be transported under one exact owner schema.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReceiptTransportProfile {
    /// Only the fixed-width owner receipt identity is retained. No replayable
    /// owner bytes are claimed by this row.
    DigestOnly,
    /// Canonical owner bytes are retained under an explicit processing cap.
    CanonicalBytes {
        /// Maximum owner bytes admitted by this schema.
        maximum_bytes: u64,
    },
}

impl ReceiptTransportProfile {
    /// Maximum canonical byte payload, absent for digest-only retention.
    #[must_use]
    pub const fn maximum_bytes(self) -> Option<u64> {
        match self {
            Self::DigestOnly => None,
            Self::CanonicalBytes { maximum_bytes } => Some(maximum_bytes),
        }
    }

    /// Whether the row declares that canonical owner bytes can be retained.
    #[must_use]
    pub const fn declares_canonical_bytes(self) -> bool {
        matches!(self, Self::CanonicalBytes { .. })
    }

    const fn tag(self) -> u8 {
        match self {
            Self::DigestOnly => TRANSPORT_DIGEST_ONLY,
            Self::CanonicalBytes { .. } => TRANSPORT_CANONICAL_BYTES,
        }
    }

    const fn encoded_limit(self) -> u64 {
        match self {
            Self::DigestOnly => 0,
            Self::CanonicalBytes { maximum_bytes } => maximum_bytes,
        }
    }
}

/// One exact, owner-fingerprinted receipt codec declaration.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReceiptSchemaDescriptor {
    family_id: String,
    wire_schema_version: u32,
    owner_identity_version: u32,
    owner_identity_domain: String,
    transport: ReceiptTransportProfile,
    owner_schema_fingerprint: ContentHash,
}

#[allow(dead_code)]
fn classify_receipt_schema_descriptor_identity_fields(descriptor: &ReceiptSchemaDescriptor) {
    let ReceiptSchemaDescriptor {
        family_id,
        wire_schema_version,
        owner_identity_version,
        owner_identity_domain,
        transport,
        owner_schema_fingerprint,
    } = descriptor;
    match transport {
        ReceiptTransportProfile::DigestOnly => {}
        ReceiptTransportProfile::CanonicalBytes { maximum_bytes } => {
            let _ = maximum_bytes;
        }
    }
    let _ = (
        family_id,
        wire_schema_version,
        owner_identity_version,
        owner_identity_domain,
        owner_schema_fingerprint,
    );
}

impl ReceiptSchemaDescriptor {
    /// Construct one bounded, canonical descriptor.
    ///
    /// # Errors
    /// [`ReceiptSchemaCatalogError::InvalidField`] or
    /// [`ReceiptSchemaCatalogError::ResourceLimit`] names the first refused
    /// field. Versions and owner fingerprints must be nonzero.
    pub fn try_new(
        family_id: impl AsRef<str>,
        wire_schema_version: u32,
        owner_identity_version: u32,
        owner_identity_domain: impl AsRef<str>,
        transport: ReceiptTransportProfile,
        owner_schema_fingerprint: ContentHash,
    ) -> Result<Self, ReceiptSchemaCatalogError> {
        let family_id = family_id.as_ref();
        let owner_identity_domain = owner_identity_domain.as_ref();
        validate_descriptor_fields(
            family_id,
            wire_schema_version,
            owner_identity_version,
            owner_identity_domain,
            transport,
            owner_schema_fingerprint,
        )?;
        Ok(Self {
            family_id: family_id.to_owned(),
            wire_schema_version,
            owner_identity_version,
            owner_identity_domain: owner_identity_domain.to_owned(),
            transport,
            owner_schema_fingerprint,
        })
    }

    /// Globally qualified family id.
    #[must_use]
    pub fn family_id(&self) -> &str {
        &self.family_id
    }

    /// Exact owner wire schema. Catalog lookup never performs fallback.
    #[must_use]
    pub const fn wire_schema_version(&self) -> u32 {
        self.wire_schema_version
    }

    /// Exact owner receipt-identity version.
    #[must_use]
    pub const fn owner_identity_version(&self) -> u32 {
        self.owner_identity_version
    }

    /// Exact owner receipt-identity domain.
    #[must_use]
    pub fn owner_identity_domain(&self) -> &str {
        &self.owner_identity_domain
    }

    /// Declared transport profile.
    #[must_use]
    pub const fn transport(&self) -> ReceiptTransportProfile {
        self.transport
    }

    /// Owner-supplied fingerprint of the exact codec/schema semantics.
    #[must_use]
    pub const fn owner_schema_fingerprint(&self) -> ContentHash {
        self.owner_schema_fingerprint
    }

    /// Domain-separated descriptor identity over every field.
    #[must_use]
    pub fn content_hash(&self) -> ContentHash {
        descriptor_hash_with_schema(
            self,
            RECEIPT_SCHEMA_DESCRIPTOR_IDENTITY_VERSION,
            RECEIPT_SCHEMA_DESCRIPTOR_IDENTITY_DOMAIN,
        )
    }

    /// Admit a retained descriptor identity only under the exact schema and
    /// fixed-width digest transport.
    #[must_use]
    pub fn admit_retained_content_hash(version: u32, bytes: &[u8]) -> Option<ContentHash> {
        admit_retained_hash(version, RECEIPT_SCHEMA_DESCRIPTOR_IDENTITY_VERSION, bytes)
    }

    fn validate(&self) -> Result<(), ReceiptSchemaCatalogError> {
        validate_descriptor_fields(
            &self.family_id,
            self.wire_schema_version,
            self.owner_identity_version,
            &self.owner_identity_domain,
            self.transport,
            self.owner_schema_fingerprint,
        )
    }

    fn encode_wire(&self, encoder: &mut Encoder) {
        encoder.u8(FIELD_FAMILY);
        encoder.string(&self.family_id);
        encoder.u8(FIELD_WIRE_SCHEMA_VERSION);
        encoder.u32(self.wire_schema_version);
        encoder.u8(FIELD_OWNER_IDENTITY_VERSION);
        encoder.u32(self.owner_identity_version);
        encoder.u8(FIELD_OWNER_IDENTITY_DOMAIN);
        encoder.string(&self.owner_identity_domain);
        encoder.u8(FIELD_TRANSPORT);
        encoder.u8(self.transport.tag());
        encoder.u64(self.transport.encoded_limit());
        encoder.u8(FIELD_OWNER_SCHEMA_FINGERPRINT);
        encoder.hash(self.owner_schema_fingerprint);
        encoder.u8(FIELD_DESCRIPTOR_HASH);
        encoder.hash(self.content_hash());
    }

    fn identity_preimage(&self) -> Vec<u8> {
        let mut encoder = Encoder::new();
        encoder.u8(FIELD_FAMILY);
        encoder.string(&self.family_id);
        encoder.u8(FIELD_WIRE_SCHEMA_VERSION);
        encoder.u32(self.wire_schema_version);
        encoder.u8(FIELD_OWNER_IDENTITY_VERSION);
        encoder.u32(self.owner_identity_version);
        encoder.u8(FIELD_OWNER_IDENTITY_DOMAIN);
        encoder.string(&self.owner_identity_domain);
        encoder.u8(FIELD_TRANSPORT);
        encoder.u8(self.transport.tag());
        encoder.u64(self.transport.encoded_limit());
        encoder.u8(FIELD_OWNER_SCHEMA_FINGERPRINT);
        encoder.hash(self.owner_schema_fingerprint);
        encoder.finish()
    }
}

fn descriptor_hash_with_schema(
    descriptor: &ReceiptSchemaDescriptor,
    identity_version: u32,
    domain: &str,
) -> ContentHash {
    let fields = descriptor.identity_preimage();
    let mut preimage = Vec::with_capacity(4 + fields.len());
    preimage.extend_from_slice(&identity_version.to_le_bytes());
    preimage.extend_from_slice(&fields);
    hash_domain(domain, &preimage)
}

/// A canonical set of exact receipt-family schema rows.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReceiptSchemaCatalog {
    catalog_version: u32,
    entries: Vec<ReceiptSchemaDescriptor>,
}

#[allow(dead_code)]
fn classify_receipt_schema_catalog_identity_fields(catalog: &ReceiptSchemaCatalog) {
    let ReceiptSchemaCatalog {
        catalog_version,
        entries,
    } = catalog;
    let _ = (catalog_version, entries);
}

impl ReceiptSchemaCatalog {
    /// Build a canonical set. Caller order is nonsemantic; rows are sorted by
    /// `(family UTF-8 bytes, wire schema version)`.
    ///
    /// # Errors
    /// Duplicate keys, reused owner domains/fingerprints, invalid descriptors,
    /// or catalog limits refuse explicitly.
    pub fn try_new(
        mut entries: Vec<ReceiptSchemaDescriptor>,
    ) -> Result<Self, ReceiptSchemaCatalogError> {
        if entries.len() > MAX_RECEIPT_SCHEMA_ENTRIES {
            return Err(resource_limit(
                "schema-entries",
                MAX_RECEIPT_SCHEMA_ENTRIES,
                entries.len(),
            ));
        }
        for entry in &entries {
            entry.validate()?;
        }
        entries.sort_by(|left, right| {
            left.family_id
                .as_bytes()
                .cmp(right.family_id.as_bytes())
                .then(left.wire_schema_version.cmp(&right.wire_schema_version))
        });

        let mut domains: BTreeMap<&str, (&str, u32)> = BTreeMap::new();
        let mut fingerprints: BTreeMap<ContentHash, (&str, u32)> = BTreeMap::new();
        let mut prior: Option<(&str, u32)> = None;
        for entry in &entries {
            let key = (entry.family_id.as_str(), entry.wire_schema_version);
            if prior == Some(key) {
                return Err(ReceiptSchemaCatalogError::DuplicateSchema {
                    family: entry.family_id.clone(),
                    wire_schema_version: entry.wire_schema_version,
                });
            }
            prior = Some(key);
            if let Some(&(family, version)) = domains.get(entry.owner_identity_domain.as_str()) {
                if (family, version) != key {
                    return Err(ReceiptSchemaCatalogError::ReusedOwnerIdentityDomain {
                        domain: entry.owner_identity_domain.clone(),
                        first_family: family.to_string(),
                        first_wire_schema_version: version,
                        duplicate_family: entry.family_id.clone(),
                        duplicate_wire_schema_version: entry.wire_schema_version,
                    });
                }
            } else {
                domains.insert(&entry.owner_identity_domain, key);
            }
            if let Some(&(family, version)) = fingerprints.get(&entry.owner_schema_fingerprint) {
                if (family, version) != key {
                    return Err(ReceiptSchemaCatalogError::ReusedOwnerSchemaFingerprint {
                        fingerprint: entry.owner_schema_fingerprint,
                        first_family: family.to_string(),
                        first_wire_schema_version: version,
                        duplicate_family: entry.family_id.clone(),
                        duplicate_wire_schema_version: entry.wire_schema_version,
                    });
                }
            } else {
                fingerprints.insert(entry.owner_schema_fingerprint, key);
            }
        }

        let catalog = Self {
            catalog_version: RECEIPT_SCHEMA_CATALOG_VERSION,
            entries,
        };
        let encoded_len = catalog.wire_preimage().len().saturating_add(32);
        if encoded_len > MAX_RECEIPT_SCHEMA_CATALOG_BYTES {
            return Err(resource_limit(
                "catalog-bytes",
                MAX_RECEIPT_SCHEMA_CATALOG_BYTES,
                encoded_len,
            ));
        }
        Ok(catalog)
    }

    /// Catalog wire schema version.
    #[must_use]
    pub const fn catalog_version(&self) -> u32 {
        self.catalog_version
    }

    /// Canonically ordered exact schema rows.
    #[must_use]
    pub fn entries(&self) -> &[ReceiptSchemaDescriptor] {
        &self.entries
    }

    /// Stable identity of this exact canonical set.
    #[must_use]
    pub fn content_hash(&self) -> ContentHash {
        catalog_hash_with_schema(
            self,
            RECEIPT_SCHEMA_CATALOG_IDENTITY_VERSION,
            RECEIPT_SCHEMA_CATALOG_IDENTITY_DOMAIN,
        )
    }

    /// Canonical binary catalog plus its fixed-width in-band identity.
    #[must_use]
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut bytes = self.wire_preimage();
        bytes.extend_from_slice(self.content_hash().as_bytes());
        bytes
    }

    /// Decode, validate, hash-check, and byte-reproduce the only supported
    /// catalog schema.
    ///
    /// # Errors
    /// Counts and strings are capped before allocation. Unknown tags,
    /// noncanonical order, and stale or self-inconsistent identities fail
    /// closed. Whole-catalog substitution is prevented only when the caller
    /// supplies an independently trusted pin to [`Self::from_bytes_verified`].
    #[allow(clippy::too_many_lines)] // one ordered fail-closed decode transcript
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, ReceiptSchemaCatalogError> {
        if bytes.len() > MAX_RECEIPT_SCHEMA_CATALOG_BYTES {
            return Err(resource_limit(
                "catalog-bytes",
                MAX_RECEIPT_SCHEMA_CATALOG_BYTES,
                bytes.len(),
            ));
        }
        let mut reader = Reader::new(bytes);
        reader.expect(CATALOG_MAGIC, "catalog magic")?;
        let catalog_version = reader.u32()?;
        if catalog_version != RECEIPT_SCHEMA_CATALOG_VERSION {
            return Err(ReceiptSchemaCatalogError::UnsupportedCatalogVersion {
                found: catalog_version,
                supported: RECEIPT_SCHEMA_CATALOG_VERSION,
            });
        }
        let entry_count = reader.count("schema-entries", MAX_RECEIPT_SCHEMA_ENTRIES)?;
        reader.require_remaining_items(
            entry_count,
            MIN_ENCODED_RECEIPT_SCHEMA_DESCRIPTOR_BYTES,
            "schema entries",
        )?;
        let mut entries = Vec::with_capacity(entry_count);
        let mut prior: Option<(String, u32)> = None;
        for _ in 0..entry_count {
            reader.expect_tag(FIELD_FAMILY, "family")?;
            let family_id = reader.string("family-id", MAX_RECEIPT_FAMILY_ID_BYTES)?;
            reader.expect_tag(FIELD_WIRE_SCHEMA_VERSION, "wire-schema-version")?;
            let wire_schema_version = reader.u32()?;
            reader.expect_tag(FIELD_OWNER_IDENTITY_VERSION, "owner-identity-version")?;
            let owner_identity_version = reader.u32()?;
            reader.expect_tag(FIELD_OWNER_IDENTITY_DOMAIN, "owner-identity-domain")?;
            let owner_identity_domain =
                reader.string("owner-identity-domain", MAX_RECEIPT_IDENTITY_DOMAIN_BYTES)?;
            reader.expect_tag(FIELD_TRANSPORT, "transport")?;
            let transport_tag = reader.u8()?;
            let transport_limit = reader.u64()?;
            let transport = match (transport_tag, transport_limit) {
                (TRANSPORT_DIGEST_ONLY, 0) => ReceiptTransportProfile::DigestOnly,
                (TRANSPORT_DIGEST_ONLY, _) => {
                    return Err(reader.malformed(
                        "digest-only transport must encode a zero canonical-byte limit",
                    ));
                }
                (TRANSPORT_CANONICAL_BYTES, maximum_bytes) => {
                    ReceiptTransportProfile::CanonicalBytes { maximum_bytes }
                }
                (tag, _) => {
                    return Err(ReceiptSchemaCatalogError::UnknownTransportTag {
                        tag,
                        at: reader.position().saturating_sub(9),
                    });
                }
            };
            reader.expect_tag(FIELD_OWNER_SCHEMA_FINGERPRINT, "owner-schema-fingerprint")?;
            let owner_schema_fingerprint = reader.hash()?;
            reader.expect_tag(FIELD_DESCRIPTOR_HASH, "descriptor-hash")?;
            let expected_descriptor_hash = reader.hash()?;

            if let Some((prior_family, prior_version)) = &prior {
                if prior_family.as_bytes() > family_id.as_bytes()
                    || (prior_family == &family_id && *prior_version >= wire_schema_version)
                {
                    return Err(ReceiptSchemaCatalogError::NonCanonicalOrder {
                        previous_family: prior_family.clone(),
                        previous_wire_schema_version: *prior_version,
                        found_family: family_id,
                        found_wire_schema_version: wire_schema_version,
                    });
                }
            }

            let descriptor = ReceiptSchemaDescriptor::try_new(
                family_id,
                wire_schema_version,
                owner_identity_version,
                owner_identity_domain,
                transport,
                owner_schema_fingerprint,
            )?;
            let actual_descriptor_hash = descriptor.content_hash();
            if actual_descriptor_hash != expected_descriptor_hash {
                return Err(ReceiptSchemaCatalogError::IdentityMismatch {
                    scope: "descriptor",
                    expected: expected_descriptor_hash,
                    actual: actual_descriptor_hash,
                });
            }
            prior = Some((descriptor.family_id.clone(), descriptor.wire_schema_version));
            entries.push(descriptor);
        }
        let expected_catalog_hash = reader.hash()?;
        reader.finish()?;
        let catalog = Self::try_new(entries)?;
        let actual_catalog_hash = catalog.content_hash();
        if actual_catalog_hash != expected_catalog_hash {
            return Err(ReceiptSchemaCatalogError::IdentityMismatch {
                scope: "catalog",
                expected: expected_catalog_hash,
                actual: actual_catalog_hash,
            });
        }
        if catalog.to_bytes() != bytes {
            return Err(ReceiptSchemaCatalogError::Malformed {
                at: bytes.len(),
                detail: "decoded catalog does not reproduce the canonical bytes".to_string(),
            });
        }
        Ok(catalog)
    }

    /// Decode under an independently trusted package/ledger catalog pin.
    pub fn from_bytes_verified(
        bytes: &[u8],
        expected: ContentHash,
    ) -> Result<Self, ReceiptSchemaCatalogError> {
        let catalog = Self::from_bytes(bytes)?;
        let actual = catalog.content_hash();
        if actual != expected {
            return Err(ReceiptSchemaCatalogError::ExternalIdentityMismatch { expected, actual });
        }
        Ok(catalog)
    }

    /// Resolve one exact family/wire-schema tuple and require the caller's
    /// descriptor identity. No nearest-version or compatibility fallback is
    /// performed.
    pub fn require_exact(
        &self,
        family: &str,
        wire_schema_version: u32,
        expected_descriptor: ContentHash,
    ) -> Result<&ReceiptSchemaDescriptor, ReceiptSchemaCatalogError> {
        validate_machine_identity("family-id", family, MAX_RECEIPT_FAMILY_ID_BYTES, true)?;
        let mut family_exists = false;
        for descriptor in &self.entries {
            if descriptor.family_id == family {
                family_exists = true;
                if descriptor.wire_schema_version == wire_schema_version {
                    let actual = descriptor.content_hash();
                    if actual != expected_descriptor {
                        return Err(ReceiptSchemaCatalogError::DescriptorMismatch {
                            family: family.to_string(),
                            wire_schema_version,
                            expected: expected_descriptor,
                            actual,
                        });
                    }
                    return Ok(descriptor);
                }
            }
        }
        if family_exists {
            Err(ReceiptSchemaCatalogError::UnsupportedWireSchema {
                family: family.to_string(),
                found: wire_schema_version,
            })
        } else {
            Err(ReceiptSchemaCatalogError::UnknownFamily {
                family: family.to_string(),
            })
        }
    }

    /// Admit a retained catalog identity only under the exact schema and
    /// fixed-width digest transport.
    #[must_use]
    pub fn admit_retained_content_hash(version: u32, bytes: &[u8]) -> Option<ContentHash> {
        admit_retained_hash(version, RECEIPT_SCHEMA_CATALOG_IDENTITY_VERSION, bytes)
    }

    fn wire_preimage(&self) -> Vec<u8> {
        let mut encoder = Encoder::new();
        encoder.bytes(CATALOG_MAGIC);
        encoder.u32(self.catalog_version);
        encoder.count(self.entries.len());
        for descriptor in &self.entries {
            descriptor.encode_wire(&mut encoder);
        }
        encoder.finish()
    }
}

fn catalog_hash_with_schema(
    catalog: &ReceiptSchemaCatalog,
    identity_version: u32,
    domain: &str,
) -> ContentHash {
    let wire = catalog.wire_preimage();
    let mut preimage = Vec::with_capacity(4 + wire.len());
    preimage.extend_from_slice(&identity_version.to_le_bytes());
    preimage.extend_from_slice(&wire);
    hash_domain(domain, &preimage)
}

/// Typed refusal at the receipt-schema metadata boundary.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ReceiptSchemaCatalogError {
    /// The retained catalog wire schema is not supported by this build.
    UnsupportedCatalogVersion {
        /// Version found.
        found: u32,
        /// Version supported.
        supported: u32,
    },
    /// One descriptor field is not canonical or meaningful.
    InvalidField {
        /// Stable field/rule name.
        field: &'static str,
        /// Teaching diagnostic.
        detail: String,
    },
    /// A byte/count/declared transport budget was exceeded.
    ResourceLimit {
        /// Stable resource name.
        resource: &'static str,
        /// Configured maximum.
        limit: u64,
        /// Exact observed amount.
        observed: u64,
    },
    /// Two rows own the same family and wire schema.
    DuplicateSchema {
        /// Duplicated family.
        family: String,
        /// Duplicated wire schema.
        wire_schema_version: u32,
    },
    /// One owner identity domain was silently aliased across two schema keys.
    ReusedOwnerIdentityDomain {
        /// Reused domain.
        domain: String,
        /// First family.
        first_family: String,
        /// First wire schema.
        first_wire_schema_version: u32,
        /// Later family.
        duplicate_family: String,
        /// Later wire schema.
        duplicate_wire_schema_version: u32,
    },
    /// One codec fingerprint was silently aliased across two schema keys.
    ReusedOwnerSchemaFingerprint {
        /// Reused fingerprint.
        fingerprint: ContentHash,
        /// First family.
        first_family: String,
        /// First wire schema.
        first_wire_schema_version: u32,
        /// Later family.
        duplicate_family: String,
        /// Later wire schema.
        duplicate_wire_schema_version: u32,
    },
    /// Encoded rows are not in canonical key order.
    NonCanonicalOrder {
        /// Previous family.
        previous_family: String,
        /// Previous wire schema.
        previous_wire_schema_version: u32,
        /// Out-of-order family.
        found_family: String,
        /// Out-of-order wire schema.
        found_wire_schema_version: u32,
    },
    /// No row owns the requested family.
    UnknownFamily {
        /// Requested family.
        family: String,
    },
    /// The family exists, but not at this exact wire schema.
    UnsupportedWireSchema {
        /// Requested family.
        family: String,
        /// Unsupported requested schema.
        found: u32,
    },
    /// The family/version exists but the caller pinned another descriptor.
    DescriptorMismatch {
        /// Requested family.
        family: String,
        /// Requested wire schema.
        wire_schema_version: u32,
        /// Caller-pinned descriptor.
        expected: ContentHash,
        /// Catalog descriptor.
        actual: ContentHash,
    },
    /// A closed transport discriminant was unknown.
    UnknownTransportTag {
        /// Unknown tag.
        tag: u8,
        /// Byte offset of the tag.
        at: usize,
    },
    /// The binary envelope is truncated, malformed, or noncanonical.
    Malformed {
        /// Byte offset at refusal.
        at: usize,
        /// Teaching diagnostic.
        detail: String,
    },
    /// An in-band descriptor or catalog identity did not reproduce.
    IdentityMismatch {
        /// `"descriptor"` or `"catalog"`.
        scope: &'static str,
        /// Retained identity.
        expected: ContentHash,
        /// Reconstructed identity.
        actual: ContentHash,
    },
    /// A caller-pinned catalog identity did not match the admitted bytes.
    ExternalIdentityMismatch {
        /// Caller requirement.
        expected: ContentHash,
        /// Reconstructed identity.
        actual: ContentHash,
    },
}

impl core::fmt::Display for ReceiptSchemaCatalogError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::UnsupportedCatalogVersion { found, supported } => write!(
                f,
                "receipt-schema catalog version {found} is unsupported; expected {supported}"
            ),
            Self::InvalidField { field, detail } => {
                write!(f, "receipt-schema field '{field}' refused: {detail}")
            }
            Self::ResourceLimit {
                resource,
                limit,
                observed,
            } => write!(
                f,
                "receipt-schema resource '{resource}' exceeds {limit} (observed {observed})"
            ),
            Self::DuplicateSchema {
                family,
                wire_schema_version,
            } => write!(
                f,
                "receipt family '{family}' wire schema {wire_schema_version} is duplicated"
            ),
            Self::ReusedOwnerIdentityDomain {
                domain,
                first_family,
                first_wire_schema_version,
                duplicate_family,
                duplicate_wire_schema_version,
            } => write!(
                f,
                "owner identity domain '{domain}' aliases {first_family}@{first_wire_schema_version} and {duplicate_family}@{duplicate_wire_schema_version}"
            ),
            Self::ReusedOwnerSchemaFingerprint {
                fingerprint,
                first_family,
                first_wire_schema_version,
                duplicate_family,
                duplicate_wire_schema_version,
            } => write!(
                f,
                "owner schema fingerprint {} aliases {first_family}@{first_wire_schema_version} and {duplicate_family}@{duplicate_wire_schema_version}",
                fingerprint.to_hex()
            ),
            Self::NonCanonicalOrder {
                previous_family,
                previous_wire_schema_version,
                found_family,
                found_wire_schema_version,
            } => write!(
                f,
                "receipt schemas are not canonical: {found_family}@{found_wire_schema_version} follows {previous_family}@{previous_wire_schema_version}"
            ),
            Self::UnknownFamily { family } => {
                write!(f, "receipt family '{family}' is not catalogued")
            }
            Self::UnsupportedWireSchema { family, found } => write!(
                f,
                "receipt family '{family}' has no exact wire schema {found}"
            ),
            Self::DescriptorMismatch {
                family,
                wire_schema_version,
                expected,
                actual,
            } => write!(
                f,
                "receipt schema {family}@{wire_schema_version} descriptor mismatch: required {}, catalogued {}",
                expected.to_hex(),
                actual.to_hex()
            ),
            Self::UnknownTransportTag { tag, at } => {
                write!(f, "unknown receipt transport tag {tag} at byte {at}")
            }
            Self::Malformed { at, detail } => {
                write!(f, "malformed receipt-schema catalog at byte {at}: {detail}")
            }
            Self::IdentityMismatch {
                scope,
                expected,
                actual,
            } => write!(
                f,
                "receipt-schema {scope} identity mismatch: retained {}, reconstructed {}",
                expected.to_hex(),
                actual.to_hex()
            ),
            Self::ExternalIdentityMismatch { expected, actual } => write!(
                f,
                "receipt-schema catalog external identity mismatch: required {}, reconstructed {}",
                expected.to_hex(),
                actual.to_hex()
            ),
        }
    }
}

impl std::error::Error for ReceiptSchemaCatalogError {}

fn validate_descriptor_fields(
    family_id: &str,
    wire_schema_version: u32,
    owner_identity_version: u32,
    owner_identity_domain: &str,
    transport: ReceiptTransportProfile,
    owner_schema_fingerprint: ContentHash,
) -> Result<(), ReceiptSchemaCatalogError> {
    validate_machine_identity("family-id", family_id, MAX_RECEIPT_FAMILY_ID_BYTES, true)?;
    if wire_schema_version == 0 {
        return Err(invalid_field("wire-schema-version", "must be nonzero"));
    }
    if owner_identity_version == 0 {
        return Err(invalid_field("owner-identity-version", "must be nonzero"));
    }
    validate_machine_identity(
        "owner-identity-domain",
        owner_identity_domain,
        MAX_RECEIPT_IDENTITY_DOMAIN_BYTES,
        false,
    )?;
    if owner_schema_fingerprint == ContentHash([0; 32]) {
        return Err(invalid_field(
            "owner-schema-fingerprint",
            "zero is not an owner codec identity",
        ));
    }
    if let ReceiptTransportProfile::CanonicalBytes { maximum_bytes } = transport {
        if maximum_bytes == 0 {
            return Err(invalid_field(
                "maximum-transport-bytes",
                "must be nonzero for canonical bytes",
            ));
        }
        if maximum_bytes > MAX_RECEIPT_TRANSPORT_BYTES {
            return Err(ReceiptSchemaCatalogError::ResourceLimit {
                resource: "maximum-transport-bytes",
                limit: MAX_RECEIPT_TRANSPORT_BYTES,
                observed: maximum_bytes,
            });
        }
    }
    Ok(())
}

fn validate_machine_identity(
    field: &'static str,
    value: &str,
    maximum_bytes: usize,
    require_colon: bool,
) -> Result<(), ReceiptSchemaCatalogError> {
    if value.is_empty() {
        return Err(invalid_field(field, "must not be empty"));
    }
    if value.trim() != value {
        return Err(invalid_field(
            field,
            "must not contain surrounding whitespace",
        ));
    }
    if value.len() > maximum_bytes {
        return Err(resource_limit(field, maximum_bytes, value.len()));
    }
    if !value.is_ascii() {
        return Err(invalid_field(
            field,
            "must be ASCII to avoid confusable machine identities",
        ));
    }
    if !value.bytes().all(|byte| {
        byte.is_ascii_lowercase()
            || byte.is_ascii_digit()
            || matches!(byte, b'.' | b'_' | b'-' | b'/' | b':')
    }) {
        return Err(invalid_field(
            field,
            "must use lowercase ASCII letters, digits, '.', '_', '-', '/', or ':'",
        ));
    }
    if !value
        .as_bytes()
        .first()
        .is_some_and(|byte| byte.is_ascii_alphanumeric())
        || !value
            .as_bytes()
            .last()
            .is_some_and(|byte| byte.is_ascii_alphanumeric())
    {
        return Err(invalid_field(
            field,
            "must start and end with an ASCII letter or digit",
        ));
    }
    if require_colon && !value.contains(':') {
        return Err(invalid_field(
            field,
            "must be globally qualified with an owner ':' separator",
        ));
    }
    if !require_colon && !value.contains('.') && !value.contains(':') {
        return Err(invalid_field(field, "must be a qualified identity domain"));
    }
    if value.split(['.', ':', '/', '_', '-']).any(str::is_empty) {
        return Err(invalid_field(
            field,
            "must not contain adjacent identity separators",
        ));
    }
    if value
        .split(['.', ':', '/', '_', '-'])
        .any(is_placeholder_component)
    {
        return Err(invalid_field(
            field,
            "contains a reserved placeholder component",
        ));
    }
    Ok(())
}

fn is_placeholder_component(component: &str) -> bool {
    matches!(
        component,
        "todo" | "tbd" | "placeholder" | "pending" | "unknown" | "none" | "na"
    )
}

fn invalid_field(field: &'static str, detail: impl Into<String>) -> ReceiptSchemaCatalogError {
    ReceiptSchemaCatalogError::InvalidField {
        field,
        detail: detail.into(),
    }
}

fn resource_limit(
    resource: &'static str,
    limit: usize,
    observed: usize,
) -> ReceiptSchemaCatalogError {
    ReceiptSchemaCatalogError::ResourceLimit {
        resource,
        limit: u64::try_from(limit).unwrap_or(u64::MAX),
        observed: u64::try_from(observed).unwrap_or(u64::MAX),
    }
}

fn admit_retained_hash(version: u32, supported: u32, bytes: &[u8]) -> Option<ContentHash> {
    if version != supported || bytes.len() != 32 {
        return None;
    }
    let mut hash = [0_u8; 32];
    hash.copy_from_slice(bytes);
    Some(ContentHash(hash))
}

struct Encoder {
    bytes: Vec<u8>,
}

impl Encoder {
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

    fn finish(self) -> Vec<u8> {
        self.bytes
    }
}

struct Reader<'a> {
    bytes: &'a [u8],
    cursor: usize,
}

impl<'a> Reader<'a> {
    fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, cursor: 0 }
    }

    fn position(&self) -> usize {
        self.cursor
    }

    fn malformed(&self, detail: impl Into<String>) -> ReceiptSchemaCatalogError {
        ReceiptSchemaCatalogError::Malformed {
            at: self.cursor,
            detail: detail.into(),
        }
    }

    fn take(&mut self, length: usize) -> Result<&'a [u8], ReceiptSchemaCatalogError> {
        let end = self
            .cursor
            .checked_add(length)
            .ok_or_else(|| self.malformed("byte offset overflow"))?;
        let bytes = self
            .bytes
            .get(self.cursor..end)
            .ok_or_else(|| self.malformed(format!("truncated field requires {length} bytes")))?;
        self.cursor = end;
        Ok(bytes)
    }

    fn expect(&mut self, expected: &[u8], field: &str) -> Result<(), ReceiptSchemaCatalogError> {
        let actual = self.take(expected.len())?;
        if actual == expected {
            Ok(())
        } else {
            Err(self.malformed(format!("invalid {field}")))
        }
    }

    fn expect_tag(&mut self, expected: u8, field: &str) -> Result<(), ReceiptSchemaCatalogError> {
        let actual = self.u8()?;
        if actual == expected {
            Ok(())
        } else {
            Err(self.malformed(format!(
                "invalid {field} field tag {actual}; expected {expected}"
            )))
        }
    }

    fn require_remaining_items(
        &self,
        count: usize,
        minimum_width: usize,
        field: &str,
    ) -> Result<(), ReceiptSchemaCatalogError> {
        let required = count
            .checked_mul(minimum_width)
            .ok_or_else(|| self.malformed(format!("{field} byte length overflow")))?;
        if self.bytes.len().saturating_sub(self.cursor) < required {
            Err(self.malformed(format!(
                "truncated {field} requires at least {required} bytes"
            )))
        } else {
            Ok(())
        }
    }

    fn u8(&mut self) -> Result<u8, ReceiptSchemaCatalogError> {
        Ok(self.take(1)?[0])
    }

    fn u32(&mut self) -> Result<u32, ReceiptSchemaCatalogError> {
        let bytes: [u8; 4] = self
            .take(4)?
            .try_into()
            .map_err(|_| self.malformed("u32 width"))?;
        Ok(u32::from_le_bytes(bytes))
    }

    fn u64(&mut self) -> Result<u64, ReceiptSchemaCatalogError> {
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
    ) -> Result<usize, ReceiptSchemaCatalogError> {
        let observed = self.u64()?;
        let limit = u64::try_from(limit).unwrap_or(u64::MAX);
        if observed > limit {
            return Err(ReceiptSchemaCatalogError::ResourceLimit {
                resource,
                limit,
                observed,
            });
        }
        usize::try_from(observed).map_err(|_| ReceiptSchemaCatalogError::ResourceLimit {
            resource,
            limit,
            observed,
        })
    }

    fn string(
        &mut self,
        resource: &'static str,
        limit: usize,
    ) -> Result<String, ReceiptSchemaCatalogError> {
        let length = self.count(resource, limit)?;
        let start = self.cursor;
        let bytes = self.take(length)?;
        let text =
            std::str::from_utf8(bytes).map_err(|error| ReceiptSchemaCatalogError::Malformed {
                at: start + error.valid_up_to(),
                detail: format!("{resource} is not UTF-8"),
            })?;
        Ok(text.to_string())
    }

    fn hash(&mut self) -> Result<ContentHash, ReceiptSchemaCatalogError> {
        let mut hash = [0_u8; 32];
        hash.copy_from_slice(self.take(32)?);
        Ok(ContentHash(hash))
    }

    fn finish(&self) -> Result<(), ReceiptSchemaCatalogError> {
        if self.cursor == self.bytes.len() {
            Ok(())
        } else {
            Err(self.malformed(format!("{} trailing bytes", self.bytes.len() - self.cursor)))
        }
    }
}

/// Owner-local descriptor identity declaration consumed by
/// `xtask check-identities`.
#[allow(dead_code)]
pub const RECEIPT_SCHEMA_DESCRIPTOR_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-package:receipt-schema-descriptor",
    "version_const=RECEIPT_SCHEMA_DESCRIPTOR_IDENTITY_VERSION",
    "version=1",
    "domain=org.frankensim.fs-package.receipt-schema-descriptor.v1",
    "domain_const=RECEIPT_SCHEMA_DESCRIPTOR_IDENTITY_DOMAIN",
    "encoder=ReceiptSchemaDescriptor::content_hash",
    "encoder_helpers=descriptor_hash_with_schema,ReceiptSchemaDescriptor::identity_preimage,ReceiptTransportProfile::tag,ReceiptTransportProfile::encoded_limit,Encoder::new,Encoder::u8,Encoder::u32,Encoder::u64,Encoder::count,Encoder::string,Encoder::hash,Encoder::finish",
    "schema_constants=RECEIPT_SCHEMA_DESCRIPTOR_IDENTITY_VERSION,RECEIPT_SCHEMA_DESCRIPTOR_IDENTITY_DOMAIN,FIELD_FAMILY,FIELD_WIRE_SCHEMA_VERSION,FIELD_OWNER_IDENTITY_VERSION,FIELD_OWNER_IDENTITY_DOMAIN,FIELD_TRANSPORT,FIELD_OWNER_SCHEMA_FINGERPRINT,TRANSPORT_DIGEST_ONLY,TRANSPORT_CANONICAL_BYTES,MAX_RECEIPT_FAMILY_ID_BYTES,MAX_RECEIPT_IDENTITY_DOMAIN_BYTES,MAX_RECEIPT_TRANSPORT_BYTES",
    "schema_functions=ReceiptSchemaDescriptor::try_new,ReceiptSchemaDescriptor::validate,ReceiptSchemaDescriptor::admit_retained_content_hash,validate_descriptor_fields,validate_machine_identity,is_placeholder_component,invalid_field,resource_limit,admit_retained_hash,crates/fs-blake3/src/lib.rs#hash_domain",
    "schema_types=std#u64",
    "schema_dependencies=none",
    "digest=blake3-256-domain-separated",
    "encoding=typed-binary",
    "sources=ReceiptSchemaDescriptor",
    "source_fields=ReceiptSchemaDescriptor.family_id:semantic,ReceiptSchemaDescriptor.wire_schema_version:semantic,ReceiptSchemaDescriptor.owner_identity_version:semantic,ReceiptSchemaDescriptor.owner_identity_domain:semantic,ReceiptSchemaDescriptor.transport:semantic,ReceiptSchemaDescriptor.owner_schema_fingerprint:semantic",
    "source_bindings=ReceiptSchemaDescriptor.family_id>family-byte-count+family-utf8,ReceiptSchemaDescriptor.wire_schema_version>wire-schema-version,ReceiptSchemaDescriptor.owner_identity_version>owner-identity-version,ReceiptSchemaDescriptor.owner_identity_domain>identity-domain-byte-count+identity-domain-utf8,ReceiptSchemaDescriptor.transport>transport-tag+maximum-transport-bytes,ReceiptSchemaDescriptor.owner_schema_fingerprint>owner-schema-fingerprint",
    "external_semantic_fields=identity-version,digest-domain,canonical-field-order,field-tag-u8,length-count-u64-le",
    "semantic_fields=identity-version,digest-domain,canonical-field-order,field-tag-u8,length-count-u64-le,family-byte-count,family-utf8,wire-schema-version,owner-identity-version,identity-domain-byte-count,identity-domain-utf8,transport-tag,maximum-transport-bytes,owner-schema-fingerprint",
    "excluded_fields=none",
    "consumers=ReceiptSchemaDescriptor::content_hash,ReceiptSchemaDescriptor::admit_retained_content_hash,ReceiptSchemaCatalog::content_hash,ReceiptSchemaCatalog::require_exact,receipt-schema-ledger-and-package-adapters",
    "mutations=identity-version:crates/fs-package/src/receipt_catalog.rs#descriptor_identity_schema_moves_version_and_domain,digest-domain:crates/fs-package/src/receipt_catalog.rs#descriptor_identity_schema_moves_version_and_domain,canonical-field-order:crates/fs-package/tests/receipt_catalog.rs#receipt_schema_catalog_round_trips_and_locks_the_independent_preimage,field-tag-u8:crates/fs-package/tests/receipt_catalog.rs#receipt_schema_catalog_decoder_refuses_hostile_and_noncanonical_bytes,length-count-u64-le:crates/fs-package/tests/receipt_catalog.rs#receipt_schema_catalog_round_trips_and_locks_the_independent_preimage,family-byte-count:crates/fs-package/tests/receipt_catalog.rs#receipt_schema_descriptor_identity_binds_every_field,family-utf8:crates/fs-package/tests/receipt_catalog.rs#receipt_schema_descriptor_identity_binds_every_field,wire-schema-version:crates/fs-package/tests/receipt_catalog.rs#receipt_schema_descriptor_identity_binds_every_field,owner-identity-version:crates/fs-package/tests/receipt_catalog.rs#receipt_schema_descriptor_identity_binds_every_field,identity-domain-byte-count:crates/fs-package/tests/receipt_catalog.rs#receipt_schema_descriptor_identity_binds_every_field,identity-domain-utf8:crates/fs-package/tests/receipt_catalog.rs#receipt_schema_descriptor_identity_binds_every_field,transport-tag:crates/fs-package/tests/receipt_catalog.rs#receipt_schema_descriptor_identity_binds_every_field,maximum-transport-bytes:crates/fs-package/tests/receipt_catalog.rs#receipt_schema_descriptor_identity_binds_every_field,owner-schema-fingerprint:crates/fs-package/tests/receipt_catalog.rs#receipt_schema_descriptor_identity_binds_every_field",
    "nonsemantic_mutations=none",
    "field_guard=classify_receipt_schema_descriptor_identity_fields",
    "transport_guard=ReceiptSchemaDescriptor::admit_retained_content_hash",
    "version_guard=crates/fs-package/src/receipt_catalog.rs#descriptor_identity_schema_moves_version_and_domain",
    "coupling_surface=fs-package:receipt-schema-descriptor",
];

/// Owner-local catalog identity declaration consumed by
/// `xtask check-identities`.
#[allow(dead_code)]
pub const RECEIPT_SCHEMA_CATALOG_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-package:receipt-schema-catalog",
    "version_const=RECEIPT_SCHEMA_CATALOG_IDENTITY_VERSION",
    "version=1",
    "domain=org.frankensim.fs-package.receipt-schema-catalog.v1",
    "domain_const=RECEIPT_SCHEMA_CATALOG_IDENTITY_DOMAIN",
    "encoder=ReceiptSchemaCatalog::content_hash",
    "encoder_helpers=catalog_hash_with_schema,ReceiptSchemaCatalog::wire_preimage,ReceiptSchemaDescriptor::encode_wire,ReceiptSchemaDescriptor::content_hash,ReceiptTransportProfile::tag,ReceiptTransportProfile::encoded_limit,Encoder::new,Encoder::bytes,Encoder::u8,Encoder::u32,Encoder::u64,Encoder::count,Encoder::string,Encoder::hash,Encoder::finish",
    "schema_constants=RECEIPT_SCHEMA_CATALOG_IDENTITY_VERSION,RECEIPT_SCHEMA_CATALOG_IDENTITY_DOMAIN,RECEIPT_SCHEMA_CATALOG_VERSION,CATALOG_MAGIC,FIELD_FAMILY,FIELD_WIRE_SCHEMA_VERSION,FIELD_OWNER_IDENTITY_VERSION,FIELD_OWNER_IDENTITY_DOMAIN,FIELD_TRANSPORT,FIELD_OWNER_SCHEMA_FINGERPRINT,FIELD_DESCRIPTOR_HASH,MIN_ENCODED_RECEIPT_SCHEMA_DESCRIPTOR_BYTES,MAX_RECEIPT_SCHEMA_ENTRIES,MAX_RECEIPT_SCHEMA_CATALOG_BYTES",
    "schema_functions=ReceiptSchemaCatalog::try_new,ReceiptSchemaCatalog::to_bytes,ReceiptSchemaCatalog::from_bytes,ReceiptSchemaCatalog::from_bytes_verified,ReceiptSchemaCatalog::require_exact,ReceiptSchemaCatalog::admit_retained_content_hash,Reader::new,Reader::take,Reader::expect,Reader::expect_tag,Reader::require_remaining_items,Reader::u8,Reader::u32,Reader::u64,Reader::count,Reader::string,Reader::hash,Reader::finish,admit_retained_hash,crates/fs-blake3/src/lib.rs#hash_domain",
    "schema_dependencies=fs-package:receipt-schema-descriptor",
    "digest=blake3-256-domain-separated",
    "encoding=canonical-transport-exact-bits",
    "sources=ReceiptSchemaCatalog",
    "source_fields=ReceiptSchemaCatalog.catalog_version:semantic,ReceiptSchemaCatalog.entries:semantic",
    "source_bindings=ReceiptSchemaCatalog.catalog_version>catalog-wire-version,ReceiptSchemaCatalog.entries>entry-count+entry-order+descriptor-fields+descriptor-hashes",
    "external_semantic_fields=identity-version,digest-domain,wire-magic,length-count-u64-le,in-band-identity",
    "semantic_fields=identity-version,digest-domain,wire-magic,length-count-u64-le,in-band-identity,catalog-wire-version,entry-count,entry-order,descriptor-fields,descriptor-hashes",
    "excluded_fields=input-entry-order:constructor-canonicalizes-entries-before-identity",
    "consumers=ReceiptSchemaCatalog::content_hash,ReceiptSchemaCatalog::to_bytes,ReceiptSchemaCatalog::from_bytes,ReceiptSchemaCatalog::from_bytes_verified,ReceiptSchemaCatalog::require_exact,ReceiptSchemaCatalog::admit_retained_content_hash,receipt-schema-ledger-and-package-adapters",
    "mutations=identity-version:crates/fs-package/src/receipt_catalog.rs#catalog_identity_schema_moves_version_and_domain,digest-domain:crates/fs-package/src/receipt_catalog.rs#catalog_identity_schema_moves_version_and_domain,wire-magic:crates/fs-package/tests/receipt_catalog.rs#receipt_schema_catalog_decoder_refuses_hostile_and_noncanonical_bytes,length-count-u64-le:crates/fs-package/tests/receipt_catalog.rs#receipt_schema_catalog_round_trips_and_locks_the_independent_preimage,in-band-identity:crates/fs-package/tests/receipt_catalog.rs#receipt_schema_catalog_decoder_refuses_hostile_and_noncanonical_bytes,catalog-wire-version:crates/fs-package/tests/receipt_catalog.rs#receipt_schema_catalog_decoder_refuses_hostile_and_noncanonical_bytes,entry-count:crates/fs-package/tests/receipt_catalog.rs#receipt_schema_catalog_round_trips_and_locks_the_independent_preimage,entry-order:crates/fs-package/tests/receipt_catalog.rs#receipt_schema_catalog_decoder_refuses_hostile_and_noncanonical_bytes,descriptor-fields:crates/fs-package/tests/receipt_catalog.rs#receipt_schema_descriptor_identity_binds_every_field,descriptor-hashes:crates/fs-package/tests/receipt_catalog.rs#receipt_schema_catalog_decoder_refuses_hostile_and_noncanonical_bytes",
    "nonsemantic_mutations=input-entry-order:crates/fs-package/tests/receipt_catalog.rs#receipt_schema_catalog_input_order_is_nonsemantic_and_lookup_is_exact",
    "field_guard=classify_receipt_schema_catalog_identity_fields",
    "transport_guard=ReceiptSchemaCatalog::from_bytes_verified",
    "version_guard=crates/fs-package/src/receipt_catalog.rs#catalog_identity_schema_moves_version_and_domain",
    "coupling_surface=fs-package:receipt-schema-catalog",
];

#[cfg(test)]
mod tests {
    use super::*;

    fn descriptor() -> ReceiptSchemaDescriptor {
        ReceiptSchemaDescriptor::try_new(
            "fs-matdb:property-usage-receipt",
            2,
            2,
            "org.frankensim.fs-matdb.property-usage-receipt.v2",
            ReceiptTransportProfile::CanonicalBytes {
                maximum_bytes: 1024 * 1024,
            },
            ContentHash([0x51; 32]),
        )
        .expect("descriptor")
    }

    #[test]
    fn descriptor_identity_schema_moves_version_and_domain() {
        let descriptor = descriptor();
        let canonical = descriptor.content_hash();
        assert_ne!(
            canonical,
            descriptor_hash_with_schema(
                &descriptor,
                RECEIPT_SCHEMA_DESCRIPTOR_IDENTITY_VERSION + 1,
                RECEIPT_SCHEMA_DESCRIPTOR_IDENTITY_DOMAIN,
            )
        );
        assert_ne!(
            canonical,
            descriptor_hash_with_schema(
                &descriptor,
                RECEIPT_SCHEMA_DESCRIPTOR_IDENTITY_VERSION,
                "org.frankensim.fs-package.receipt-schema-descriptor.foreign",
            )
        );
    }

    #[test]
    fn catalog_identity_schema_moves_version_and_domain() {
        let catalog = ReceiptSchemaCatalog::try_new(vec![descriptor()]).expect("catalog");
        let canonical = catalog.content_hash();
        assert_ne!(
            canonical,
            catalog_hash_with_schema(
                &catalog,
                RECEIPT_SCHEMA_CATALOG_IDENTITY_VERSION + 1,
                RECEIPT_SCHEMA_CATALOG_IDENTITY_DOMAIN,
            )
        );
        assert_ne!(
            canonical,
            catalog_hash_with_schema(
                &catalog,
                RECEIPT_SCHEMA_CATALOG_IDENTITY_VERSION,
                "org.frankensim.fs-package.receipt-schema-catalog.foreign",
            )
        );
    }
}
