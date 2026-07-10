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
//! the stream is framed by magic `fsid`, the schema version, and the
//! length-prefixed artifact kind. Properties the battery certifies:
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

/// Frame magic: identifies a canonical identity stream.
const MAGIC: &[u8; 4] = b"fsid";

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

    /// Canonical display form: `fsid-v<version>:<kind>:<root hex>`.
    #[must_use]
    pub fn hex(&self) -> String {
        format!("fsid-v{}:{}:{:016x}", self.version, self.kind, self.root)
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

impl core::error::Error for IdentError {}

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

impl IdentityBuilder {
    /// Start an identity for one artifact kind (e.g. `"solver-snapshot"`,
    /// `"evidence-package"`). The kind and schema version are framed
    /// into the hashed stream.
    #[must_use]
    pub fn new(kind: &str) -> Self {
        let mut buf = Vec::with_capacity(128);
        buf.extend_from_slice(MAGIC);
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

fn push_len_bytes(buf: &mut Vec<u8>, bytes: &[u8]) {
    buf.extend_from_slice(&(bytes.len() as u64).to_le_bytes());
    buf.extend_from_slice(bytes);
}
