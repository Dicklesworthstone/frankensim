//! ident — CANONICAL REPLAY IDENTITY (bead gp3.14): one versioned,
//! typed, length-prefixed encoding for every replay-bearing artifact
//! (plans, certificates, snapshots, evidence packages, golden metric
//! streams), replacing ad hoc delimiter concatenation.
//!
//! WHY: fresh-eyes failures showed identities that under-bind — fields
//! joined with `|` or bare concatenation let DISTINCT inputs share an
//! identity (`("ab","c")` vs `("a","bc")`), and unversioned encodings
//! cannot evolve without silently re-keying history. The fs-cheb /
//! flagship golden misses are the concrete evidence: arithmetic
//! schedule changes propagated to semantic consumers with no
//! dependency-aware identity naming the affected artifacts.
//!
//! THE ENCODING (schema v1). A field is
//! `tag(1B) | key_len(u64 LE) | key | val_len(u64 LE) | val_bytes`;
//! the stream is framed by [`REPLAY_IDENTITY_DOMAIN`]'s exact bytes, the
//! schema version, and the length-prefixed artifact kind. Properties the
//! battery certifies:
//! - length prefixes kill delimiter/split collisions;
//! - the TYPE TAG is hashed, so `str "1"`, `u64 1`, and `bytes b"1"`
//!   have different roots (type confusion changes the identity);
//! - field ORDER is semantic (ordered operation parameters);
//! - floats travel as bit patterns (`-0.0` and `0.0` differ; NaN is
//!   representable and stable) — never as formatted text;
//! - the schema version is part of the root, and consumers verify
//!   declared versions FAIL-CLOSED via [`check_version`].
//!
//! FIELD DISCIPLINE (the producer inventory's three classes):
//! - SEMANTIC fields go through the typed `push_*` methods and bind
//!   the root: algorithm + schema versions, deterministic mode, full
//!   logical RNG identity, machine/ISA class where the claim needs
//!   it, budgets and units, representation choices, ordered operation
//!   parameters, certificate regimes, parent artifact roots
//!   ([`IdentityBuilder::child`]), and dependency implementation
//!   identities.
//! - PROVENANCE-ONLY / DELIBERATELY-EXCLUDED fields (wall-clock,
//!   hostnames, transient handles) are declared via
//!   [`IdentityBuilder::exclude`]: never hashed, but RECORDED so the
//!   exclusion is documented in code and testable (mutation coverage
//!   asserts they do not move the root).
//!
//! The 64-bit root is [`crate::fnv1a64`] over the canonical bytes —
//! the house digest until the BLAKE3-class ledger hash supersedes it;
//! [`ReplayIdentity::canonical_bytes`] exposes the exact stream so
//! stronger digests can bind the SAME encoding.

use core::fmt;

/// The current identity schema version. Bump ONLY with a migration
/// note in the producing crate's CONTRACT (changing the encoding
/// re-keys every root — a justified-golden-bump event by definition).
pub const IDENT_SCHEMA_VERSION: u32 = 1;

/// Domain prefix framed into every canonical replay identity.
pub const REPLAY_IDENTITY_DOMAIN: &str = "fsid";

/// Owner-local declaration consumed by `xtask check-identities`.
pub const REPLAY_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-obs:replay-identity-frame",
    "version_const=IDENT_SCHEMA_VERSION",
    "version=1",
    "domain=fsid",
    "domain_const=REPLAY_IDENTITY_DOMAIN",
    "encoder=IdentityBuilder::finish",
    "encoder_helpers=IdentityBuilder::new",
    "schema_constants=IDENT_SCHEMA_VERSION,REPLAY_IDENTITY_DOMAIN,TAG_STR,TAG_U64,TAG_I64,TAG_F64_BITS,TAG_BYTES,TAG_CHILD,TAG_FLAG",
    "schema_functions=IdentityBuilder::str,IdentityBuilder::u64,IdentityBuilder::i64,IdentityBuilder::f64_bits,IdentityBuilder::bytes,IdentityBuilder::flag,IdentityBuilder::child,IdentityBuilder::child_root64,BoundedIdentityBuilder::new,BoundedIdentityBuilder::field,BoundedIdentityBuilder::str,BoundedIdentityBuilder::u64,BoundedIdentityBuilder::i64,BoundedIdentityBuilder::f64_bits,BoundedIdentityBuilder::bytes,BoundedIdentityBuilder::flag,BoundedIdentityBuilder::child,BoundedIdentityBuilder::child_root64,BoundedIdentityBuilder::finish,framed_length,encoded_field_len,enforce_canonical_budget,reserve_canonical_bytes,append_len_bytes,push_len_bytes,crates/fs-obs/src/lib.rs#fnv1a64",
    "schema_dependencies=none",
    "digest=fnv1a64",
    "encoding=typed-binary",
    "sources=ReplayIdentity",
    "source_fields=ReplayIdentity.version:semantic,ReplayIdentity.kind:semantic,ReplayIdentity.root:derived:fnv-root-of-canonical-bytes,ReplayIdentity.bytes:semantic,ReplayIdentity.excluded:nonsemantic:audit-only",
    "source_bindings=ReplayIdentity.version>schema-version,ReplayIdentity.kind>kind,ReplayIdentity.bytes>field-order+type-tag+field-key+exact-value-bytes+child-version-root",
    "external_semantic_fields=domain-prefix",
    "semantic_fields=domain-prefix,schema-version,kind,field-order,type-tag,field-key,exact-value-bytes,child-version-root",
    "excluded_fields=max-canonical-bytes:admission-budget-only",
    "consumers=IdentityBuilder,BoundedIdentityBuilder,fs-substrate,fs-evidence,fs-flagship-e2e,fs-report,fs-geom,fs-verify",
    "mutations=domain-prefix:crates/fs-obs/tests/ident.rs#ident_006_unknown_schema_versions_fail_closed,schema-version:crates/fs-obs/tests/ident.rs#ident_006_unknown_schema_versions_fail_closed,kind:crates/fs-obs/tests/ident.rs#ident_001_every_semantic_field_moves_the_root,field-order:crates/fs-obs/tests/ident.rs#ident_004_type_confusion_and_order_are_semantic,type-tag:crates/fs-obs/tests/ident.rs#ident_004_type_confusion_and_order_are_semantic,field-key:crates/fs-obs/tests/ident.rs#ident_003_delimiter_and_split_collisions_refused,exact-value-bytes:crates/fs-obs/tests/ident.rs#ident_004_type_confusion_and_order_are_semantic,child-version-root:crates/fs-obs/tests/ident.rs#ident_005_dependency_children_propagate",
    "nonsemantic_mutations=ReplayIdentity.excluded:crates/fs-obs/tests/ident.rs#ident_002_documented_exclusions_do_not_move_the_root,max-canonical-bytes:crates/fs-obs/tests/ident.rs#ident_007_bounded_builder_is_byte_exact_and_refuses_at_the_cap",
    "field_guard=classify_replay_identity_fields",
    "transport_guard=IdentityBuilder::field",
    "version_guard=crates/fs-obs/tests/ident.rs#ident_006_unknown_schema_versions_fail_closed",
    "coupling_surface=fs-obs:replay-identity-frame",
];

/// Typed field tags — hashed, so type confusion changes the root.
const TAG_STR: u8 = 0x01;
const TAG_U64: u8 = 0x02;
const TAG_I64: u8 = 0x03;
const TAG_F64_BITS: u8 = 0x04;
const TAG_BYTES: u8 = 0x05;
const TAG_CHILD: u8 = 0x06;
const TAG_FLAG: u8 = 0x07;

/// A finished canonical identity: the versioned root plus the exact
/// byte stream it binds (for stronger digests and forensic replay).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReplayIdentity {
    version: u32,
    kind: String,
    root: u64,
    bytes: Vec<u8>,
    excluded: Vec<(&'static str, &'static str)>,
}

#[allow(dead_code)]
fn classify_replay_identity_fields(identity: &ReplayIdentity) {
    let ReplayIdentity {
        version: _,
        kind: _,
        root: _,
        bytes: _,
        excluded: _,
    } = identity;
}

impl ReplayIdentity {
    /// The schema version this identity was produced under.
    #[must_use]
    pub fn version(&self) -> u32 {
        self.version
    }

    /// The artifact kind the identity names.
    #[must_use]
    pub fn kind(&self) -> &str {
        &self.kind
    }

    /// The 64-bit root (FNV-1a over the canonical bytes).
    #[must_use]
    pub fn root(&self) -> u64 {
        self.root
    }

    /// The exact canonical byte stream (feed to a stronger digest to
    /// bind the same encoding).
    #[must_use]
    pub fn canonical_bytes(&self) -> &[u8] {
        &self.bytes
    }

    /// The documented exclusions `(field, why)` — the audit trail for
    /// "deliberately not part of the identity".
    #[must_use]
    pub fn exclusions(&self) -> &[(&'static str, &'static str)] {
        &self.excluded
    }

    /// Canonical display form:
    /// `<REPLAY_IDENTITY_DOMAIN>-v<version>:<kind>:<root hex>`.
    #[must_use]
    pub fn hex(&self) -> String {
        format!(
            "{}-v{}:{}:{:016x}",
            REPLAY_IDENTITY_DOMAIN, self.version, self.kind, self.root
        )
    }
}

/// Why an identity could not be accepted (Decalogue P10).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IdentError {
    /// The declared schema version is not the supported one — a
    /// verifier that hashes by the wrong rules would mint or miss
    /// identities silently, so this FAILS CLOSED.
    UnknownSchemaVersion {
        /// The version the artifact declared.
        declared: u32,
        /// The version this build supports.
        supported: u32,
    },
}

/// Why a bounded canonical identity could not be constructed.
///
/// This is separate from [`IdentError`]: `IdentError` validates identities
/// that already exist, while this type reports producer-side resource and
/// allocation refusal before any partial identity can escape.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IdentityBuildError {
    /// The resulting canonical stream would exceed the producer's explicit
    /// byte budget.
    CanonicalBytesExceeded {
        /// Canonical bytes required after the attempted append.
        requested: usize,
        /// Maximum canonical bytes admitted by the producer.
        limit: usize,
    },
    /// A key, value, or kind length cannot be represented by schema v1's
    /// unsigned 64-bit framing.
    FramedLengthNotRepresentable {
        /// Native length that could not be framed.
        length: usize,
    },
    /// Checked arithmetic could not represent the resulting canonical length.
    CanonicalLengthOverflow,
    /// Capacity reservation failed before the attempted append mutated the
    /// canonical stream.
    AllocationFailed {
        /// Canonical bytes retained or targeted when reservation failed.
        requested: usize,
    },
}

impl fmt::Display for IdentError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            IdentError::UnknownSchemaVersion {
                declared,
                supported,
            } => write!(
                f,
                "replay-identity schema v{declared} is not supported (this build knows \
                 v{supported}); refusing to hash by the wrong rules — upgrade the verifier \
                 or re-produce the artifact"
            ),
        }
    }
}

impl fmt::Display for IdentityBuildError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::CanonicalBytesExceeded { requested, limit } => write!(
                f,
                "canonical replay identity requires {requested} bytes; limit is {limit}"
            ),
            Self::FramedLengthNotRepresentable { length } => write!(
                f,
                "canonical replay identity length {length} does not fit schema-v1 u64 framing"
            ),
            Self::CanonicalLengthOverflow => {
                f.write_str("canonical replay identity length arithmetic overflowed")
            }
            Self::AllocationFailed { requested } => write!(
                f,
                "canonical replay identity allocation failed while retaining or targeting \
                 {requested} canonical bytes"
            ),
        }
    }
}

impl core::error::Error for IdentError {}

impl core::error::Error for IdentityBuildError {}

/// Verify a DECLARED schema version before trusting any root computed
/// under it.
///
/// # Errors
/// [`IdentError::UnknownSchemaVersion`] for anything other than the
/// exact supported version — unknown versions fail closed.
pub fn check_version(declared: u32) -> Result<(), IdentError> {
    if declared == IDENT_SCHEMA_VERSION {
        Ok(())
    } else {
        Err(IdentError::UnknownSchemaVersion {
            declared,
            supported: IDENT_SCHEMA_VERSION,
        })
    }
}

/// Builder for one artifact's canonical identity. Field ORDER is
/// semantic; every `push` method is a typed, length-prefixed append.
#[derive(Debug, Clone)]
pub struct IdentityBuilder {
    kind: String,
    buf: Vec<u8>,
    excluded: Vec<(&'static str, &'static str)>,
}

/// Allocation-fallible canonical identity builder with an explicit byte cap.
///
/// The schema-v1 byte stream is exactly the one emitted by
/// [`IdentityBuilder`]. Each consuming append validates framing and the final
/// length, reserves the complete field, and only then mutates the buffer. An
/// errored builder is consumed, so safe code cannot finish a partial identity.
#[derive(Debug)]
pub struct BoundedIdentityBuilder {
    kind: String,
    buf: Vec<u8>,
    excluded: Vec<(&'static str, &'static str)>,
    max_canonical_bytes: usize,
}

impl IdentityBuilder {
    /// Start an identity for one artifact kind (e.g. `"solver-snapshot"`,
    /// `"evidence-package"`). The kind and schema version are framed
    /// into the hashed stream.
    #[must_use]
    pub fn new(kind: &str) -> Self {
        let mut buf = Vec::with_capacity(128);
        buf.extend_from_slice(REPLAY_IDENTITY_DOMAIN.as_bytes());
        buf.extend_from_slice(&IDENT_SCHEMA_VERSION.to_le_bytes());
        push_len_bytes(&mut buf, kind.as_bytes());
        IdentityBuilder {
            kind: kind.to_string(),
            buf,
            excluded: Vec::new(),
        }
    }

    fn field(&mut self, tag: u8, key: &str, value: &[u8]) {
        self.buf.push(tag);
        push_len_bytes(&mut self.buf, key.as_bytes());
        push_len_bytes(&mut self.buf, value);
    }

    /// A semantic string field (algorithm names, modes, units).
    #[must_use]
    pub fn str(mut self, key: &str, value: &str) -> Self {
        self.field(TAG_STR, key, value.as_bytes());
        self
    }

    /// A semantic unsigned integer (sizes, seeds, counters, versions).
    #[must_use]
    pub fn u64(mut self, key: &str, value: u64) -> Self {
        self.field(TAG_U64, key, &value.to_le_bytes());
        self
    }

    /// A semantic signed integer.
    #[must_use]
    pub fn i64(mut self, key: &str, value: i64) -> Self {
        self.field(TAG_I64, key, &value.to_le_bytes());
        self
    }

    /// A semantic float, bound by BIT PATTERN (`-0.0 != 0.0`, NaN
    /// stable) — never formatted text.
    #[must_use]
    pub fn f64_bits(mut self, key: &str, value: f64) -> Self {
        self.field(TAG_F64_BITS, key, &value.to_bits().to_le_bytes());
        self
    }

    /// A semantic raw-bytes field (payload digests, packed tables).
    #[must_use]
    pub fn bytes(mut self, key: &str, value: &[u8]) -> Self {
        self.field(TAG_BYTES, key, value);
        self
    }

    /// A semantic boolean (deterministic mode, feature switches).
    #[must_use]
    pub fn flag(mut self, key: &str, value: bool) -> Self {
        self.field(TAG_FLAG, key, &[u8::from(value)]);
        self
    }

    /// A parent artifact root or dependency implementation identity —
    /// the dependency-aware edge: when an upstream identity changes,
    /// every identity that bound it as a child changes with it, naming
    /// the downstream goldens that need re-verification.
    #[must_use]
    pub fn child(mut self, key: &str, root: &ReplayIdentity) -> Self {
        let mut val = Vec::with_capacity(12);
        val.extend_from_slice(&root.version.to_le_bytes());
        val.extend_from_slice(&root.root.to_le_bytes());
        self.field(TAG_CHILD, key, &val);
        self
    }

    /// A raw 64-bit child root (for pre-existing hashes — golden
    /// constants, fnv content hashes — that are not yet
    /// [`ReplayIdentity`] values).
    #[must_use]
    pub fn child_root64(mut self, key: &str, root: u64) -> Self {
        self.field(TAG_CHILD, key, &root.to_le_bytes());
        self
    }

    /// DOCUMENT a field as deliberately excluded from the identity
    /// (provenance-only: wall-clock, hostnames, transient handles).
    /// Never hashed; recorded so the exclusion is auditable and the
    /// mutation battery can assert it does not move the root.
    #[must_use]
    pub fn exclude(mut self, key: &'static str, why: &'static str) -> Self {
        self.excluded.push((key, why));
        self
    }

    /// Finish: bind the canonical bytes into the versioned root.
    #[must_use]
    pub fn finish(self) -> ReplayIdentity {
        let root = crate::fnv1a64(&self.buf);
        ReplayIdentity {
            version: IDENT_SCHEMA_VERSION,
            kind: self.kind,
            root,
            bytes: self.buf,
            excluded: self.excluded,
        }
    }
}

impl BoundedIdentityBuilder {
    /// Start a schema-v1 identity whose complete canonical stream may retain at
    /// most `max_canonical_bytes` bytes.
    ///
    /// # Errors
    /// [`IdentityBuildError`] when the kind cannot be framed, its header exceeds
    /// the byte cap, or capacity cannot be reserved.
    pub fn new(kind: &str, max_canonical_bytes: usize) -> Result<Self, IdentityBuildError> {
        let kind_len = framed_length(kind.len())?;
        let header_len = REPLAY_IDENTITY_DOMAIN
            .len()
            .checked_add(core::mem::size_of::<u32>())
            .and_then(|length| length.checked_add(core::mem::size_of::<u64>()))
            .and_then(|length| length.checked_add(kind.len()))
            .ok_or(IdentityBuildError::CanonicalLengthOverflow)?;
        enforce_canonical_budget(header_len, max_canonical_bytes)?;

        let mut owned_kind = String::new();
        owned_kind.try_reserve_exact(kind.len()).map_err(|_| {
            IdentityBuildError::AllocationFailed {
                requested: header_len,
            }
        })?;
        owned_kind.push_str(kind);

        let mut buf = Vec::new();
        reserve_canonical_bytes(&mut buf, header_len, header_len)?;
        buf.extend_from_slice(REPLAY_IDENTITY_DOMAIN.as_bytes());
        buf.extend_from_slice(&IDENT_SCHEMA_VERSION.to_le_bytes());
        append_len_bytes(&mut buf, kind_len, kind.as_bytes());
        Ok(Self {
            kind: owned_kind,
            buf,
            excluded: Vec::new(),
            max_canonical_bytes,
        })
    }

    fn field(mut self, tag: u8, key: &str, value: &[u8]) -> Result<Self, IdentityBuildError> {
        let key_len = framed_length(key.len())?;
        let value_len = framed_length(value.len())?;
        let added = encoded_field_len(key.len(), value.len())?;
        let requested = self
            .buf
            .len()
            .checked_add(added)
            .ok_or(IdentityBuildError::CanonicalLengthOverflow)?;
        enforce_canonical_budget(requested, self.max_canonical_bytes)?;
        reserve_canonical_bytes(&mut self.buf, added, requested)?;
        self.buf.push(tag);
        append_len_bytes(&mut self.buf, key_len, key.as_bytes());
        append_len_bytes(&mut self.buf, value_len, value);
        Ok(self)
    }

    /// Append a semantic string field.
    ///
    /// # Errors
    /// [`IdentityBuildError`] on framing, byte-budget, or allocation refusal.
    pub fn str(self, key: &str, value: &str) -> Result<Self, IdentityBuildError> {
        self.field(TAG_STR, key, value.as_bytes())
    }

    /// Append a semantic unsigned integer.
    ///
    /// # Errors
    /// [`IdentityBuildError`] on framing, byte-budget, or allocation refusal.
    pub fn u64(self, key: &str, value: u64) -> Result<Self, IdentityBuildError> {
        self.field(TAG_U64, key, &value.to_le_bytes())
    }

    /// Append a semantic signed integer.
    ///
    /// # Errors
    /// [`IdentityBuildError`] on framing, byte-budget, or allocation refusal.
    pub fn i64(self, key: &str, value: i64) -> Result<Self, IdentityBuildError> {
        self.field(TAG_I64, key, &value.to_le_bytes())
    }

    /// Append a semantic float by exact bit pattern.
    ///
    /// # Errors
    /// [`IdentityBuildError`] on framing, byte-budget, or allocation refusal.
    pub fn f64_bits(self, key: &str, value: f64) -> Result<Self, IdentityBuildError> {
        self.field(TAG_F64_BITS, key, &value.to_bits().to_le_bytes())
    }

    /// Append a semantic raw-bytes field.
    ///
    /// # Errors
    /// [`IdentityBuildError`] on framing, byte-budget, or allocation refusal.
    pub fn bytes(self, key: &str, value: &[u8]) -> Result<Self, IdentityBuildError> {
        self.field(TAG_BYTES, key, value)
    }

    /// Append a semantic boolean.
    ///
    /// # Errors
    /// [`IdentityBuildError`] on framing, byte-budget, or allocation refusal.
    pub fn flag(self, key: &str, value: bool) -> Result<Self, IdentityBuildError> {
        self.field(TAG_FLAG, key, &[u8::from(value)])
    }

    /// Append a schema-versioned child identity.
    ///
    /// # Errors
    /// [`IdentityBuildError`] on framing, byte-budget, or allocation refusal.
    pub fn child(self, key: &str, root: &ReplayIdentity) -> Result<Self, IdentityBuildError> {
        let mut value = [0_u8; 12];
        value[..4].copy_from_slice(&root.version.to_le_bytes());
        value[4..].copy_from_slice(&root.root.to_le_bytes());
        self.field(TAG_CHILD, key, &value)
    }

    /// Append a raw 64-bit child root.
    ///
    /// # Errors
    /// [`IdentityBuildError`] on framing, byte-budget, or allocation refusal.
    pub fn child_root64(self, key: &str, root: u64) -> Result<Self, IdentityBuildError> {
        self.field(TAG_CHILD, key, &root.to_le_bytes())
    }

    /// Record a deliberately excluded provenance field without hashing it.
    ///
    /// # Errors
    /// [`IdentityBuildError::AllocationFailed`] when the audit record cannot be
    /// retained. Canonical-byte accounting is unchanged because exclusions are
    /// not hashed.
    pub fn exclude(
        mut self,
        key: &'static str,
        why: &'static str,
    ) -> Result<Self, IdentityBuildError> {
        self.excluded
            .try_reserve(1)
            .map_err(|_| IdentityBuildError::AllocationFailed {
                requested: self.buf.len(),
            })?;
        self.excluded.push((key, why));
        Ok(self)
    }

    /// Finish the already-reserved canonical stream and compute its v1 root.
    #[must_use]
    pub fn finish(self) -> ReplayIdentity {
        let root = crate::fnv1a64(&self.buf);
        ReplayIdentity {
            version: IDENT_SCHEMA_VERSION,
            kind: self.kind,
            root,
            bytes: self.buf,
            excluded: self.excluded,
        }
    }
}

fn framed_length(length: usize) -> Result<u64, IdentityBuildError> {
    u64::try_from(length).map_err(|_| IdentityBuildError::FramedLengthNotRepresentable { length })
}

fn encoded_field_len(key_len: usize, value_len: usize) -> Result<usize, IdentityBuildError> {
    1_usize
        .checked_add(core::mem::size_of::<u64>())
        .and_then(|length| length.checked_add(key_len))
        .and_then(|length| length.checked_add(core::mem::size_of::<u64>()))
        .and_then(|length| length.checked_add(value_len))
        .ok_or(IdentityBuildError::CanonicalLengthOverflow)
}

fn enforce_canonical_budget(requested: usize, limit: usize) -> Result<(), IdentityBuildError> {
    if requested > limit {
        Err(IdentityBuildError::CanonicalBytesExceeded { requested, limit })
    } else {
        Ok(())
    }
}

fn reserve_canonical_bytes(
    buf: &mut Vec<u8>,
    additional: usize,
    requested: usize,
) -> Result<(), IdentityBuildError> {
    buf.try_reserve_exact(additional)
        .map_err(|_| IdentityBuildError::AllocationFailed { requested })
}

fn append_len_bytes(buf: &mut Vec<u8>, framed_len: u64, bytes: &[u8]) {
    buf.extend_from_slice(&framed_len.to_le_bytes());
    buf.extend_from_slice(bytes);
}

fn push_len_bytes(buf: &mut Vec<u8>, bytes: &[u8]) {
    buf.extend_from_slice(&(bytes.len() as u64).to_le_bytes());
    buf.extend_from_slice(bytes);
}

#[cfg(test)]
mod bounded_tests {
    use super::*;

    #[test]
    fn encoded_field_length_overflow_is_structured() {
        assert_eq!(
            encoded_field_len(usize::MAX, 1),
            Err(IdentityBuildError::CanonicalLengthOverflow)
        );
    }

    #[test]
    fn capacity_overflow_is_structured_before_mutation() {
        let mut buf = Vec::new();
        assert_eq!(
            reserve_canonical_bytes(&mut buf, usize::MAX, usize::MAX),
            Err(IdentityBuildError::AllocationFailed {
                requested: usize::MAX,
            })
        );
        assert!(buf.is_empty());
    }

    #[test]
    fn native_lengths_are_converted_without_truncation() {
        if usize::BITS <= u64::BITS {
            let expected = u64::try_from(usize::MAX).expect("guarded native length fits u64");
            assert_eq!(framed_length(usize::MAX), Ok(expected));
        } else {
            assert_eq!(
                framed_length(usize::MAX),
                Err(IdentityBuildError::FramedLengthNotRepresentable { length: usize::MAX })
            );
        }
    }
}
