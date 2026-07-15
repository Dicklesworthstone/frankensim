//! Strong, schema-typed identities over the workspace BLAKE3 implementation.
//!
//! This module deliberately does not replace [`crate::ContentHash`]. That
//! compatibility digest is still used by persisted ledgers and by identity
//! dialects whose migration is tracked separately. New authority-bearing
//! code should instead use the nominal types here:
//!
//! - [`ContentId`] names exact raw bytes and says nothing about meaning or
//!   authenticity;
//! - [`SemanticId`], [`WireContentId`], [`EvidenceNodeId`], [`EntityId`],
//!   [`SourceByteId`], [`SourceId`], [`ModelId`], [`CheckerId`],
//!   [`VerifierId`], and [`KeyPolicyId`] cannot be interchanged, even when
//!   their 32 digest bytes happen to match;
//! - the schema marker type is part of each Rust type, so two domains or two
//!   schema versions cannot be mixed accidentally;
//! - [`legacy::LegacyProvenanceV1`] retains an old FNV value without any
//!   widening or conversion into a strong identity;
//! - [`AuthorityRef`] separates presented, verified, and admitted authority
//!   from content/semantic consistency.
//!
//! # Compile-time separation
//!
//! Raw content cannot be passed as a semantic identity:
//!
//! ```compile_fail
//! use fs_blake3::identity::{ContentId, SemanticId};
//!
//! enum Demo {}
//! impl fs_blake3::identity::CanonicalSchema for Demo {
//!     const DOMAIN: &'static str = "org.example.demo.v1";
//!     const NAME: &'static str = "demo";
//!     const VERSION: u32 = 1;
//!     const CONTEXT: &'static str = "example";
//!     const FIELDS: &'static [fs_blake3::identity::FieldSpec] = &[];
//! }
//!
//! fn needs_semantic(_: SemanticId<Demo>) {}
//! needs_semantic(ContentId::of_bytes(b"demo"));
//! ```
//!
//! Semantic domains are nominal, not runtime strings:
//!
//! ```compile_fail
//! use fs_blake3::identity::{CanonicalSchema, FieldSpec, SemanticId};
//!
//! enum A {}
//! enum B {}
//! impl CanonicalSchema for A {
//!     const DOMAIN: &'static str = "org.example.a.v1";
//!     const NAME: &'static str = "a";
//!     const VERSION: u32 = 1;
//!     const CONTEXT: &'static str = "example";
//!     const FIELDS: &'static [FieldSpec] = &[];
//! }
//! impl CanonicalSchema for B {
//!     const DOMAIN: &'static str = "org.example.b.v1";
//!     const NAME: &'static str = "b";
//!     const VERSION: u32 = 1;
//!     const CONTEXT: &'static str = "example";
//!     const FIELDS: &'static [FieldSpec] = &[];
//! }
//!
//! fn needs_b(_: SemanticId<B>) {}
//! fn misuse(value: SemanticId<A>) { needs_b(value); }
//! ```
//!
//! Presented authority is not admitted authority:
//!
//! ```compile_fail
//! use fs_blake3::identity::{
//!     Admitted, AuthorityRef, CanonicalSchema, Presented, StrongIdentity,
//! };
//!
//! fn needs_admitted<I, V, P>(_: AuthorityRef<I, V, P, Admitted>)
//! where
//!     I: StrongIdentity,
//!     V: CanonicalSchema,
//!     P: CanonicalSchema,
//! {}
//! fn misuse<I, V, P>(value: AuthorityRef<I, V, P, Presented>)
//! where
//!     I: StrongIdentity,
//!     V: CanonicalSchema,
//!     P: CanonicalSchema,
//! {
//!     needs_admitted(value);
//! }
//! ```

use core::fmt;
use core::hash::{Hash, Hasher};
use core::marker::PhantomData;

use crate::{Blake3, ContentHash, derive_key_hasher, hash_bytes};

/// Version of the canonical binary frame defined by this module.
pub const CANONICAL_FRAME_VERSION: u32 = 1;

/// BLAKE3 derive-key context for complete canonical identity frames.
pub const CANONICAL_IDENTITY_HASH_DOMAIN: &str =
    "org.frankensim.fs-blake3.canonical-identity-frame.v1";

/// BLAKE3 derive-key context for non-recursive schema descriptors.
pub const SCHEMA_ID_HASH_DOMAIN: &str = "org.frankensim.fs-blake3.schema-id.v1";

const CANONICAL_MAGIC: &[u8; 8] = b"FSID\0\0\0\x01";
// v2 (bead sj31i.52.10): field descriptors bind the expected child
// role/schema recursively, so every v1 schema id is a DIFFERENT value
// with an explicit no-authority crosswalk boundary (no silent
// reinterpretation of v1 identities as child-bound).
const SCHEMA_MAGIC: &[u8; 8] = b"FSSCHEM\x02";
const FIELD_MARKER: u8 = 0xf0;
const END_MARKER: u8 = 0xff;
const FLOAT_POLICY_FINITE_EXACT_BITS: u8 = 1;

/// Owner-local declaration consumed by `xtask check-identities`.
pub const SCHEMA_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-blake3:schema-id",
    "version_const=CANONICAL_FRAME_VERSION",
    "version=1",
    "domain=org.frankensim.fs-blake3.schema-id.v1",
    "domain_const=SCHEMA_ID_HASH_DOMAIN",
    "encoder=SchemaId::for_schema",
    "encoder_helpers=write_schema_descriptor,hash_len_bytes",
    "schema_constants=CANONICAL_FRAME_VERSION,SCHEMA_ID_HASH_DOMAIN,SCHEMA_MAGIC",
    "schema_functions=SchemaId::for_schema,FieldSpec::name,FieldSpec::wire_type,FieldSpec::presence,Presence::tag,WireType::tag,crates/fs-blake3/src/lib.rs#Blake3::finalize,crates/fs-blake3/src/lib.rs#Blake3::update,crates/fs-blake3/src/lib.rs#derive_key_hasher",
    "schema_dependencies=none",
    "digest=fs-blake3",
    "encoding=typed-binary",
    "sources=SchemaDescriptorSource",
    "source_fields=SchemaDescriptorSource.domain:semantic,SchemaDescriptorSource.name:semantic,SchemaDescriptorSource.version:semantic,SchemaDescriptorSource.context:semantic,SchemaDescriptorSource.fields:semantic",
    "source_bindings=SchemaDescriptorSource.domain>domain,SchemaDescriptorSource.name>schema-name,SchemaDescriptorSource.version>schema-version,SchemaDescriptorSource.context>context,SchemaDescriptorSource.fields>declared-field-count+field-order+ordered-field-name+wire-type+presence",
    "external_semantic_fields=schema-descriptor-magic,schema-descriptor-version",
    "semantic_fields=schema-descriptor-magic,schema-descriptor-version,domain,schema-name,schema-version,context,declared-field-count,field-order,ordered-field-name,wire-type,presence",
    "excluded_fields=none",
    "consumers=CanonicalEncoder,IdentityReceipt,StrongIdentity,SchemaId",
    "mutations=schema-descriptor-magic:crates/fs-blake3/tests/identity.rs#schema_descriptor_and_every_header_field_move_identity,schema-descriptor-version:crates/fs-blake3/tests/identity.rs#schema_descriptor_and_every_header_field_move_identity,domain:crates/fs-blake3/tests/identity.rs#schema_descriptor_and_every_header_field_move_identity,schema-name:crates/fs-blake3/tests/identity.rs#schema_descriptor_and_every_header_field_move_identity,schema-version:crates/fs-blake3/tests/identity.rs#schema_descriptor_and_every_header_field_move_identity,context:crates/fs-blake3/tests/identity.rs#schema_descriptor_and_every_header_field_move_identity,declared-field-count:crates/fs-blake3/tests/identity.rs#schema_descriptor_and_every_header_field_move_identity,field-order:crates/fs-blake3/tests/identity.rs#schema_descriptor_and_every_header_field_move_identity,ordered-field-name:crates/fs-blake3/tests/identity.rs#schema_descriptor_and_every_header_field_move_identity,wire-type:crates/fs-blake3/tests/identity.rs#schema_descriptor_and_every_header_field_move_identity,presence:crates/fs-blake3/tests/identity.rs#schema_descriptor_and_every_header_field_move_identity",
    "nonsemantic_mutations=none",
    "field_guard=classify_schema_descriptor_fields",
    "transport_guard=SchemaId::for_schema",
    "version_guard=crates/fs-blake3/tests/identity.rs#schema_versions_are_nominal_and_digest_distinct",
    "coupling_surface=fs-blake3:schema-id",
];

/// Owner-local declaration consumed by `xtask check-identities`.
pub const CANONICAL_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-blake3:canonical-identity-frame",
    "version_const=CANONICAL_FRAME_VERSION",
    "version=1",
    "domain=org.frankensim.fs-blake3.canonical-identity-frame.v1",
    "domain_const=CANONICAL_IDENTITY_HASH_DOMAIN",
    "encoder=CanonicalEncoder::finish",
    "encoder_helpers=CanonicalEncoder::new_internal,CanonicalEncoder::begin_field,CanonicalEncoder::append,CanonicalEncoder::bytes_stream,CanonicalEncoder::ordered_bytes,CanonicalEncoder::canonical_set,CanonicalEncoder::child,CanonicalEncoder::ordered_children",
    "schema_constants=CANONICAL_FRAME_VERSION,CANONICAL_IDENTITY_HASH_DOMAIN,CANONICAL_MAGIC,FIELD_MARKER,END_MARKER,FLOAT_POLICY_FINITE_EXACT_BITS",
    "schema_functions=CanonicalEncoder::finish,CanonicalEncoder::new_internal,CanonicalEncoder::begin_field,CanonicalEncoder::append,CanonicalEncoder::utf8,CanonicalEncoder::bytes,CanonicalEncoder::u64,CanonicalEncoder::i64,CanonicalEncoder::flag,CanonicalEncoder::finite_f64,CanonicalEncoder::optional_bytes,CanonicalEncoder::variant,CanonicalEncoder::bytes_stream,CanonicalEncoder::ordered_bytes,CanonicalEncoder::canonical_set,CanonicalEncoder::child,CanonicalEncoder::ordered_children,SchemaId::for_schema,Presence::tag,WireType::tag,IdentityRole::tag,crates/fs-blake3/src/lib.rs#Blake3::new,crates/fs-blake3/src/lib.rs#Blake3::finalize,crates/fs-blake3/src/lib.rs#Blake3::update,crates/fs-blake3/src/lib.rs#ContentHash::as_bytes,crates/fs-blake3/src/lib.rs#derive_key_hasher",
    "schema_dependencies=fs-blake3:schema-id",
    "digest=fs-blake3",
    "encoding=typed-binary",
    "sources=CanonicalIdentityHeaderSource",
    "source_fields=CanonicalIdentityHeaderSource.role:semantic,CanonicalIdentityHeaderSource.domain:semantic,CanonicalIdentityHeaderSource.schema_name:semantic,CanonicalIdentityHeaderSource.schema_id:semantic,CanonicalIdentityHeaderSource.version:semantic,CanonicalIdentityHeaderSource.context:semantic,CanonicalIdentityHeaderSource.fields:semantic",
    "source_bindings=CanonicalIdentityHeaderSource.role>role-tag,CanonicalIdentityHeaderSource.domain>domain,CanonicalIdentityHeaderSource.schema_name>schema-name,CanonicalIdentityHeaderSource.schema_id>schema-id,CanonicalIdentityHeaderSource.version>semantic-version,CanonicalIdentityHeaderSource.context>context,CanonicalIdentityHeaderSource.fields>declared-field-count+ordered-field-schema",
    "external_semantic_fields=canonical-magic,canonical-frame-version,float-policy,canonical-field-stream",
    "semantic_fields=canonical-magic,canonical-frame-version,role-tag,domain,schema-name,schema-id,semantic-version,context,float-policy,declared-field-count,ordered-field-schema,canonical-field-stream",
    "excluded_fields=display-json-debug-text:display-transport-only,admission-budgets:admission-budget-only,cancellation-schedule:execution-schedule-only",
    "consumers=CanonicalEncoder,IdentityReceipt,StrongIdentity,AuthorityRef,IdentityAuditRecord",
    "mutations=canonical-magic:crates/fs-blake3/tests/identity.rs#manual_frame_parity_and_header_mutation_sensitivity,canonical-frame-version:crates/fs-blake3/tests/identity.rs#manual_frame_parity_and_header_mutation_sensitivity,role-tag:crates/fs-blake3/tests/identity.rs#roles_domains_versions_and_raw_content_are_separate,domain:crates/fs-blake3/tests/identity.rs#manual_frame_parity_and_header_mutation_sensitivity,schema-name:crates/fs-blake3/tests/identity.rs#manual_frame_parity_and_header_mutation_sensitivity,schema-id:crates/fs-blake3/tests/identity.rs#manual_frame_parity_and_header_mutation_sensitivity,semantic-version:crates/fs-blake3/tests/identity.rs#manual_frame_parity_and_header_mutation_sensitivity,context:crates/fs-blake3/tests/identity.rs#manual_frame_parity_and_header_mutation_sensitivity,float-policy:crates/fs-blake3/tests/identity.rs#manual_frame_parity_and_header_mutation_sensitivity,declared-field-count:crates/fs-blake3/tests/identity.rs#manual_frame_parity_and_header_mutation_sensitivity,ordered-field-schema:crates/fs-blake3/tests/identity.rs#manual_frame_parity_and_header_mutation_sensitivity,canonical-field-stream:crates/fs-blake3/tests/identity.rs#every_semantic_field_is_mutation_sensitive",
    "nonsemantic_mutations=display-json-debug-text:crates/fs-blake3/tests/identity.rs#display_and_debug_are_not_hash_inputs,admission-budgets:crates/fs-blake3/tests/identity.rs#budgets_do_not_move_an_admitted_identity,cancellation-schedule:crates/fs-blake3/tests/identity.rs#stream_partition_and_non_cancelling_probes_are_invariant",
    "field_guard=classify_canonical_identity_header_fields",
    "transport_guard=CanonicalEncoder::new_internal",
    "version_guard=crates/fs-blake3/tests/identity.rs#schema_versions_are_nominal_and_digest_distinct",
    "coupling_surface=fs-blake3:canonical-identity-frame",
];

/// Type-level description of one registered canonical identity schema.
///
/// Implementations should be zero-sized marker types with hardcoded,
/// globally unique, versioned constants. Runtime domain strings are not
/// accepted by [`CanonicalEncoder`].
pub trait CanonicalSchema: 'static {
    /// Globally unique, versioned semantic domain.
    const DOMAIN: &'static str;
    /// Stable human-readable schema name used in receipts and the frame.
    const NAME: &'static str;
    /// Semantic schema version. Unknown versions are different marker types.
    const VERSION: u32;
    /// Stable purpose/context string; never host, clock, or display text.
    const CONTEXT: &'static str;
    /// Complete top-level field schema in exact canonical order.
    const FIELDS: &'static [FieldSpec];
}

/// Whether a field is required or explicitly optional.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum Presence {
    /// Exactly one value must be encoded.
    Required = 1,
    /// A presence tag is encoded before the optional value.
    Optional = 2,
}

impl Presence {
    /// Stable v1 binary tag. Changing a tag requires a frame-version bump.
    #[must_use]
    pub const fn tag(self) -> u8 {
        match self {
            Self::Required => 1,
            Self::Optional => 2,
        }
    }
}

/// Canonical wire grammar for a top-level field.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum WireType {
    /// Exact UTF-8 bytes with no implicit normalization.
    Utf8 = 1,
    /// Arbitrary bytes.
    Bytes = 2,
    /// Unsigned 64-bit little-endian integer.
    U64 = 3,
    /// Signed 64-bit little-endian integer.
    I64 = 4,
    /// Boolean encoded as exactly zero or one.
    Bool = 5,
    /// Finite IEEE-754 bits; signed zero is preserved.
    FiniteF64 = 6,
    /// Numeric variant tag followed by a length-framed byte payload.
    Variant = 7,
    /// Ordered length-framed byte sequence.
    OrderedBytes = 8,
    /// Strictly increasing, duplicate-free byte set.
    CanonicalSet = 9,
    /// One full typed child identity.
    Child = 10,
    /// Ordered sequence of one typed child role/schema.
    OrderedChildren = 11,
}

impl WireType {
    /// Stable v1 binary tag. Changing a tag requires a frame-version bump.
    #[must_use]
    pub const fn tag(self) -> u8 {
        match self {
            Self::Utf8 => 1,
            Self::Bytes => 2,
            Self::U64 => 3,
            Self::I64 => 4,
            Self::Bool => 5,
            Self::FiniteF64 => 6,
            Self::Variant => 7,
            Self::OrderedBytes => 8,
            Self::CanonicalSet => 9,
            Self::Child => 10,
            Self::OrderedChildren => 11,
        }
    }
}

/// The parent-declared binding for a child field (bead sj31i.52.10):
/// the EXACT expected child role and complete schema identity —
/// domain, name, version, context, and the child's full recursive
/// field schema. A parent schema admits only this child type; wrong
/// role/domain/name/version/context/nested schema refuses at encode
/// time, and the binding is part of the parent schema-id preimage, so
/// changing the expected child type changes the parent [`SchemaId`].
///
/// Field-schema comparison uses `&'static` slice POINTER identity: the
/// binding built by [`ChildSpec::for_identity`] captures the child
/// schema's own `FIELDS` static, which is exactly what the encoder
/// sees for a matching type — while a structurally identical but
/// DISTINCT schema type stays non-confusable.
#[derive(Debug, Clone, Copy)]
pub struct ChildSpec {
    role: IdentityRole,
    domain: &'static str,
    name: &'static str,
    version: u32,
    context: &'static str,
    fields: &'static [FieldSpec],
}

impl ChildSpec {
    /// The binding for exactly the identity type `J`.
    #[must_use]
    pub const fn for_identity<J: StrongIdentity>() -> Self {
        Self {
            role: J::ROLE,
            domain: <J::Schema as CanonicalSchema>::DOMAIN,
            name: <J::Schema as CanonicalSchema>::NAME,
            version: <J::Schema as CanonicalSchema>::VERSION,
            context: <J::Schema as CanonicalSchema>::CONTEXT,
            fields: <J::Schema as CanonicalSchema>::FIELDS,
        }
    }

    /// Check an encoder-supplied identity type against this binding,
    /// returning the first mismatched dimension.
    fn matches<J: StrongIdentity>(&self) -> Result<(), &'static str> {
        if self.role.tag() != J::ROLE.tag() {
            return Err("child role");
        }
        if self.domain != <J::Schema as CanonicalSchema>::DOMAIN {
            return Err("child schema domain");
        }
        if self.name != <J::Schema as CanonicalSchema>::NAME {
            return Err("child schema name");
        }
        if self.version != <J::Schema as CanonicalSchema>::VERSION {
            return Err("child schema version");
        }
        if self.context != <J::Schema as CanonicalSchema>::CONTEXT {
            return Err("child schema context");
        }
        if !core::ptr::eq(self.fields, <J::Schema as CanonicalSchema>::FIELDS) {
            return Err("child field schema");
        }
        Ok(())
    }
}

impl PartialEq for ChildSpec {
    fn eq(&self, other: &Self) -> bool {
        // Pointer identity on the recursive tail keeps equality total
        // even for (pathological) cyclic `&'static` schema graphs.
        self.role.tag() == other.role.tag()
            && self.domain == other.domain
            && self.name == other.name
            && self.version == other.version
            && self.context == other.context
            && core::ptr::eq(self.fields, other.fields)
    }
}

impl Eq for ChildSpec {}

impl Hash for ChildSpec {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.role.tag().hash(state);
        self.domain.hash(state);
        self.name.hash(state);
        self.version.hash(state);
        self.context.hash(state);
        (self.fields.as_ptr() as usize).hash(state);
        self.fields.len().hash(state);
    }
}

/// One field in a [`CanonicalSchema`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct FieldSpec {
    name: &'static str,
    wire_type: WireType,
    presence: Presence,
    child: Option<&'static ChildSpec>,
}

impl FieldSpec {
    /// Declare one required field. Child wire types MUST declare their
    /// expected child binding through [`FieldSpec::child_of`] /
    /// [`FieldSpec::ordered_children_of`]; an unbound child field is a
    /// COMPILE-TIME error in const schema declarations:
    ///
    /// ```compile_fail
    /// use fs_blake3::identity::{FieldSpec, WireType};
    ///
    /// const F: FieldSpec = FieldSpec::required("lineage", WireType::Child);
    /// ```
    #[must_use]
    pub const fn required(name: &'static str, wire_type: WireType) -> Self {
        match wire_type {
            WireType::Child | WireType::OrderedChildren => {
                panic!(
                    "child fields must declare their expected child schema via child_of/ordered_children_of"
                )
            }
            _ => {}
        }
        Self {
            name,
            wire_type,
            presence: Presence::Required,
            child: None,
        }
    }

    /// Declare one required child field bound to exactly `spec`.
    #[must_use]
    pub const fn child_of(name: &'static str, spec: &'static ChildSpec) -> Self {
        Self {
            name,
            wire_type: WireType::Child,
            presence: Presence::Required,
            child: Some(spec),
        }
    }

    /// Declare one required ordered-children field bound to exactly
    /// `spec` (empty collections still validate against the binding).
    #[must_use]
    pub const fn ordered_children_of(name: &'static str, spec: &'static ChildSpec) -> Self {
        Self {
            name,
            wire_type: WireType::OrderedChildren,
            presence: Presence::Required,
            child: Some(spec),
        }
    }

    /// Declare one explicitly optional byte field.
    ///
    /// Canonical-frame v1 deliberately exposes no generic optional constructor:
    /// optional presence is representable only for [`WireType::Bytes`], the
    /// wire grammar implemented by [`CanonicalEncoder::optional_bytes`].
    ///
    /// ```compile_fail
    /// use fs_blake3::identity::{FieldSpec, WireType};
    ///
    /// let _ = FieldSpec::optional("value", WireType::U64);
    /// ```
    #[must_use]
    pub const fn optional_bytes(name: &'static str) -> Self {
        Self {
            name,
            wire_type: WireType::Bytes,
            presence: Presence::Optional,
            child: None,
        }
    }

    /// The declared child binding, when this is a child field.
    #[must_use]
    pub const fn child_spec(self) -> Option<&'static ChildSpec> {
        self.child
    }

    /// Stable field name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        self.name
    }

    /// Declared wire type.
    #[must_use]
    pub const fn wire_type(self) -> WireType {
        self.wire_type
    }

    /// Required/optional policy.
    #[must_use]
    pub const fn presence(self) -> Presence {
        self.presence
    }
}

/// Caller key for the next exact schema field.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Field {
    ordinal: u32,
    name: &'static str,
}

impl Field {
    /// Construct a field selector. The encoder checks both values against the
    /// static schema before hashing any field bytes.
    #[must_use]
    pub const fn new(ordinal: u32, name: &'static str) -> Self {
        Self { ordinal, name }
    }
}

/// Semantic role encoded into every typed identity frame.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum IdentityRole {
    /// General normalized semantic identity.
    Semantic = 1,
    /// Exact versioned canonical transport bytes.
    WireContent = 2,
    /// Ordered evidence node identity.
    EvidenceNode = 3,
    /// Durable entity/lineage identity.
    Entity = 4,
    /// Exact source-byte identity under a source schema.
    SourceBytes = 5,
    /// Source record identity.
    Source = 6,
    /// Model identity.
    Model = 7,
    /// Checker identity.
    Checker = 8,
    /// Schema descriptor identity.
    Schema = 9,
    /// Verifier implementation/policy identity.
    Verifier = 10,
    /// Key-policy identity.
    KeyPolicy = 11,
    /// Normalized problem meaning.
    ProblemSemantic = 12,
}

impl IdentityRole {
    /// Stable v1 binary tag. Changing a tag requires a frame-version bump.
    #[must_use]
    pub const fn tag(self) -> u8 {
        match self {
            Self::Semantic => 1,
            Self::WireContent => 2,
            Self::EvidenceNode => 3,
            Self::Entity => 4,
            Self::SourceBytes => 5,
            Self::Source => 6,
            Self::Model => 7,
            Self::Checker => 8,
            Self::Schema => 9,
            Self::Verifier => 10,
            Self::KeyPolicy => 11,
            Self::ProblemSemantic => 12,
        }
    }
}

/// Exact raw bytes under plain BLAKE3 mode.
///
/// Equality proves only that the digests match under the BLAKE3
/// collision-resistance assumption. It does not prove origin, authority, or
/// semantic equivalence.
#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct ContentId(ContentHash);

impl ContentId {
    /// Hash exact bytes in plain BLAKE3 mode.
    #[must_use]
    pub fn of_bytes(bytes: &[u8]) -> Self {
        Self(hash_bytes(bytes))
    }

    /// Parse a retained raw content ID. Parsing is not verification.
    #[must_use]
    pub fn parse_slice(bytes: &[u8]) -> Option<Self> {
        ContentHash::from_slice(bytes).map(Self)
    }

    /// Parse 64 hexadecimal digits. Parsing is not verification.
    #[must_use]
    pub fn parse_hex(value: &str) -> Option<Self> {
        ContentHash::from_hex(value).map(Self)
    }

    /// Exact digest bytes.
    #[must_use]
    pub fn as_bytes(&self) -> &[u8; 32] {
        self.0.as_bytes()
    }

    /// Lowercase hexadecimal rendering.
    #[must_use]
    pub fn to_hex(self) -> String {
        self.0.to_hex()
    }
}

impl fmt::Display for ContentId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(&self.0, f)
    }
}

impl fmt::Debug for ContentId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "ContentId({self})")
    }
}

mod strong_identity_sealed {
    pub trait Sealed {}
}

/// Public behavior shared by canonical strong identities.
///
/// This trait intentionally has no constructor and no conversion to
/// [`ContentHash`] or [`ContentId`]. Its private supertrait closes the role
/// universe to the nominal wrappers defined in this module.
#[allow(private_bounds)] // Intentional sealed-trait boundary.
pub trait StrongIdentity:
    strong_identity_sealed::Sealed + Copy + Eq + Ord + Hash + fmt::Debug + fmt::Display + 'static
{
    /// Static schema marker carried by this Rust type.
    type Schema: CanonicalSchema;
    /// Non-interchangeable semantic role.
    const ROLE: IdentityRole;
    /// Exact 32 digest bytes.
    fn as_bytes(&self) -> &[u8; 32];
    /// Strict typed parsing of retained digest bytes. This does not add trust.
    fn parse_slice(bytes: &[u8]) -> Option<Self>;
    /// Lowercase hexadecimal rendering.
    fn to_hex(self) -> String;
}

macro_rules! strong_identity {
    ($(#[$meta:meta])* $name:ident, $role:expr) => {
        $(#[$meta])*
        pub struct $name<D: CanonicalSchema> {
            digest: ContentHash,
            marker: PhantomData<fn() -> D>,
        }

        impl<D: CanonicalSchema> $name<D> {
            fn from_digest(digest: ContentHash) -> Self {
                Self { digest, marker: PhantomData }
            }

            /// Parse 64 hexadecimal digits under this exact role/schema type.
            /// Parsing is not verification or authority admission.
            #[must_use]
            pub fn parse_hex(value: &str) -> Option<Self> {
                ContentHash::from_hex(value).map(Self::from_digest)
            }
        }

        impl<D: CanonicalSchema> Copy for $name<D> {}
        impl<D: CanonicalSchema> Clone for $name<D> {
            fn clone(&self) -> Self { *self }
        }
        impl<D: CanonicalSchema> PartialEq for $name<D> {
            fn eq(&self, other: &Self) -> bool { self.digest == other.digest }
        }
        impl<D: CanonicalSchema> Eq for $name<D> {}
        impl<D: CanonicalSchema> PartialOrd for $name<D> {
            fn partial_cmp(&self, other: &Self) -> Option<core::cmp::Ordering> {
                Some(self.cmp(other))
            }
        }
        impl<D: CanonicalSchema> Ord for $name<D> {
            fn cmp(&self, other: &Self) -> core::cmp::Ordering {
                self.digest.cmp(&other.digest)
            }
        }
        impl<D: CanonicalSchema> Hash for $name<D> {
            fn hash<H: Hasher>(&self, state: &mut H) { self.digest.hash(state); }
        }
        impl<D: CanonicalSchema> fmt::Display for $name<D> {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                fmt::Display::fmt(&self.digest, f)
            }
        }
        impl<D: CanonicalSchema> fmt::Debug for $name<D> {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                write!(f, "{}<{}>({})", stringify!($name), D::DOMAIN, self.digest)
            }
        }
        impl<D: CanonicalSchema> StrongIdentity for $name<D> {
            type Schema = D;
            const ROLE: IdentityRole = $role;
            fn as_bytes(&self) -> &[u8; 32] { self.digest.as_bytes() }
            fn parse_slice(bytes: &[u8]) -> Option<Self> {
                ContentHash::from_slice(bytes).map(Self::from_digest)
            }
            fn to_hex(self) -> String { self.digest.to_hex() }
        }
        impl<D: CanonicalSchema> strong_identity_sealed::Sealed for $name<D> {}
    };
}

strong_identity!(
    /// Normalized semantic identity under schema `D`.
    SemanticId,
    IdentityRole::Semantic
);
strong_identity!(
    /// Exact versioned canonical transport identity under schema `D`.
    WireContentId,
    IdentityRole::WireContent
);
strong_identity!(
    /// Ordered evidence-node identity under schema `D`.
    EvidenceNodeId,
    IdentityRole::EvidenceNode
);
strong_identity!(
    /// Durable entity/lineage identity under schema `D`.
    EntityId,
    IdentityRole::Entity
);
strong_identity!(
    /// Exact source-byte identity under schema `D`.
    SourceByteId,
    IdentityRole::SourceBytes
);
strong_identity!(
    /// Source record identity under schema `D`.
    SourceId,
    IdentityRole::Source
);
strong_identity!(
    /// Model identity under schema `D`.
    ModelId,
    IdentityRole::Model
);
strong_identity!(
    /// Checker identity under schema `D`.
    CheckerId,
    IdentityRole::Checker
);
strong_identity!(
    /// Verifier implementation/policy identity under schema `D`.
    VerifierId,
    IdentityRole::Verifier
);
strong_identity!(
    /// Key-policy identity under schema `D`.
    KeyPolicyId,
    IdentityRole::KeyPolicy
);
strong_identity!(
    /// Normalized problem-meaning identity under schema `D`.
    ProblemSemanticId,
    IdentityRole::ProblemSemantic
);

/// Non-recursive identity of the complete static descriptor for schema `D`.
///
/// The descriptor is hashed directly under [`SCHEMA_ID_HASH_DOMAIN`], so a
/// canonical identity frame can safely include this value without defining a
/// schema in terms of a frame that already requires itself.
pub struct SchemaId<D: CanonicalSchema> {
    digest: ContentHash,
    marker: PhantomData<fn() -> D>,
}

impl<D: CanonicalSchema> SchemaId<D> {
    fn from_digest(digest: ContentHash) -> Self {
        Self {
            digest,
            marker: PhantomData,
        }
    }

    /// Compute the schema descriptor identity without allocation.
    ///
    /// This names the descriptor exactly as declared; it does not admit the
    /// descriptor for canonical construction. [`CanonicalEncoder`] separately
    /// validates descriptor structure, resource limits, and cancellation.
    #[must_use]
    pub fn for_schema() -> Self {
        let source = SchemaDescriptorSource {
            domain: D::DOMAIN,
            name: D::NAME,
            version: D::VERSION,
            context: D::CONTEXT,
            fields: D::FIELDS,
        };
        let mut hasher = derive_key_hasher(SCHEMA_ID_HASH_DOMAIN);
        match write_schema_descriptor(&source, |bytes| {
            hasher.update(bytes);
            Ok::<(), core::convert::Infallible>(())
        }) {
            Ok(()) => {}
            Err(never) => match never {},
        }
        Self::from_digest(hasher.finalize())
    }

    /// Parse 64 hexadecimal digits under this exact schema marker.
    /// Parsing does not prove that the value equals [`Self::for_schema`].
    #[must_use]
    pub fn parse_hex(value: &str) -> Option<Self> {
        ContentHash::from_hex(value).map(Self::from_digest)
    }
}

impl<D: CanonicalSchema> Copy for SchemaId<D> {}
impl<D: CanonicalSchema> Clone for SchemaId<D> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<D: CanonicalSchema> PartialEq for SchemaId<D> {
    fn eq(&self, other: &Self) -> bool {
        self.digest == other.digest
    }
}
impl<D: CanonicalSchema> Eq for SchemaId<D> {}
impl<D: CanonicalSchema> PartialOrd for SchemaId<D> {
    fn partial_cmp(&self, other: &Self) -> Option<core::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl<D: CanonicalSchema> Ord for SchemaId<D> {
    fn cmp(&self, other: &Self) -> core::cmp::Ordering {
        self.digest.cmp(&other.digest)
    }
}
impl<D: CanonicalSchema> Hash for SchemaId<D> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.digest.hash(state);
    }
}
impl<D: CanonicalSchema> fmt::Display for SchemaId<D> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(&self.digest, f)
    }
}
impl<D: CanonicalSchema> fmt::Debug for SchemaId<D> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "SchemaId<{}>({})", D::DOMAIN, self.digest)
    }
}
impl<D: CanonicalSchema> StrongIdentity for SchemaId<D> {
    type Schema = D;
    const ROLE: IdentityRole = IdentityRole::Schema;

    fn as_bytes(&self) -> &[u8; 32] {
        self.digest.as_bytes()
    }

    fn parse_slice(bytes: &[u8]) -> Option<Self> {
        ContentHash::from_slice(bytes).map(Self::from_digest)
    }

    fn to_hex(self) -> String {
        self.digest.to_hex()
    }
}
impl<D: CanonicalSchema> strong_identity_sealed::Sealed for SchemaId<D> {}

struct SchemaDescriptorSource<'a> {
    domain: &'a str,
    name: &'a str,
    version: u32,
    context: &'a str,
    fields: &'a [FieldSpec],
}

#[allow(dead_code)]
fn classify_schema_descriptor_fields(source: &SchemaDescriptorSource<'_>) {
    let SchemaDescriptorSource {
        domain: _,
        name: _,
        version: _,
        context: _,
        fields: _,
    } = source;
}

/// Nested child-binding descriptors deeper than this are POISON-tagged
/// in the schema-id preimage (still deterministic and well-defined);
/// the encoder separately refuses to construct under such bindings.
const MAX_SCHEMA_CHILD_DEPTH: u32 = 16;

fn write_schema_descriptor<E>(
    source: &SchemaDescriptorSource<'_>,
    mut update: impl FnMut(&[u8]) -> Result<(), E>,
) -> Result<(), E> {
    write_schema_descriptor_at(source, &mut update, 0)
}

fn write_schema_descriptor_at<E>(
    source: &SchemaDescriptorSource<'_>,
    update: &mut impl FnMut(&[u8]) -> Result<(), E>,
    depth: u32,
) -> Result<(), E> {
    update(SCHEMA_MAGIC)?;
    update(&CANONICAL_FRAME_VERSION.to_le_bytes())?;
    hash_len_bytes(update, source.domain.as_bytes())?;
    hash_len_bytes(update, source.name.as_bytes())?;
    update(&source.version.to_le_bytes())?;
    hash_len_bytes(update, source.context.as_bytes())?;
    update(&(source.fields.len() as u64).to_le_bytes())?;
    for field in source.fields {
        hash_len_bytes(update, field.name.as_bytes())?;
        update(&[field.wire_type.tag(), field.presence.tag()])?;
        // bead sj31i.52.10: the expected child binding is part of the
        // parent schema identity, recursively — changing the expected
        // child type changes the parent SchemaId.
        match field.child {
            None => update(&[0u8])?,
            Some(_) if depth >= MAX_SCHEMA_CHILD_DEPTH => update(&[2u8])?,
            Some(child) => {
                update(&[1u8, child.role.tag()])?;
                write_schema_descriptor_at(
                    &SchemaDescriptorSource {
                        domain: child.domain,
                        name: child.name,
                        version: child.version,
                        context: child.context,
                        fields: child.fields,
                    },
                    update,
                    depth + 1,
                )?;
            }
        }
    }
    Ok(())
}

fn hash_len_bytes<E>(
    update: &mut impl FnMut(&[u8]) -> Result<(), E>,
    value: &[u8],
) -> Result<(), E> {
    update(&(value.len() as u64).to_le_bytes())?;
    update(value)
}

/// Explicit resource envelope for one canonical identity operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct CanonicalLimits {
    max_canonical_bytes: u64,
    max_field_bytes: u64,
    max_fields: u32,
    max_collection_items: u64,
    cancellation_poll_bytes: u32,
}

impl CanonicalLimits {
    /// Construct an explicit resource envelope.
    #[must_use]
    pub const fn new(
        max_canonical_bytes: u64,
        max_field_bytes: u64,
        max_fields: u32,
        max_collection_items: u64,
        cancellation_poll_bytes: u32,
    ) -> Self {
        Self {
            max_canonical_bytes,
            max_field_bytes,
            max_fields,
            max_collection_items,
            cancellation_poll_bytes,
        }
    }

    /// Maximum complete canonical-frame bytes.
    #[must_use]
    pub const fn max_canonical_bytes(self) -> u64 {
        self.max_canonical_bytes
    }

    /// Maximum payload bytes for one field or collection item.
    #[must_use]
    pub const fn max_field_bytes(self) -> u64 {
        self.max_field_bytes
    }

    /// Maximum top-level fields.
    #[must_use]
    pub const fn max_fields(self) -> u32 {
        self.max_fields
    }

    /// Maximum items in one collection and chunks in one streamed byte field.
    #[must_use]
    pub const fn max_collection_items(self) -> u64 {
        self.max_collection_items
    }

    /// Maximum payload bytes between cancellation polls.
    #[must_use]
    pub const fn cancellation_poll_bytes(self) -> u32 {
        self.cancellation_poll_bytes
    }
}

impl Default for CanonicalLimits {
    fn default() -> Self {
        Self::new(1 << 20, 1 << 18, 256, 16_384, 4096)
    }
}

/// Resource dimension that refused canonical construction.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum LimitKind {
    /// Complete canonical frame.
    CanonicalBytes,
    /// One field or collection item.
    FieldBytes,
    /// Top-level field count.
    Fields,
    /// Collection item count.
    CollectionItems,
    /// Non-semantic chunk count in one streamed byte field.
    StreamChunks,
}

/// Fail-closed canonical construction error.
///
/// Every fallible encoder operation consumes the encoder. An error therefore
/// leaves no value on which [`CanonicalEncoder::finish`] could be called.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CanonicalError {
    /// A resource envelope is internally invalid.
    InvalidLimits(&'static str),
    /// A schema descriptor is empty or otherwise invalid.
    InvalidSchemaDescriptor(&'static str),
    /// Checked length arithmetic overflowed.
    LengthOverflow,
    /// An explicit resource budget was exceeded.
    LimitExceeded {
        /// Refused resource dimension.
        kind: LimitKind,
        /// Requested total.
        requested: u64,
        /// Configured limit.
        limit: u64,
    },
    /// Caller selected a field other than the exact next schema field.
    FieldOrder {
        /// Expected ordinal.
        expected: u32,
        /// Supplied ordinal.
        actual: u32,
    },
    /// Caller field name differed from the static schema.
    FieldName,
    /// Field method did not match the declared wire grammar.
    WireType,
    /// Required/optional encoding did not match the static schema.
    Presence,
    /// Finish was attempted before every declared field was encoded.
    MissingFields {
        /// Declared count.
        expected: u32,
        /// Encoded count.
        actual: u32,
    },
    /// A streamed length or collection count did not match its declaration.
    DeclaredLengthMismatch {
        /// Declared value.
        declared: u64,
        /// Observed value.
        observed: u64,
    },
    /// A generic semantic float was NaN or infinite.
    NonFiniteFloat {
        /// Exact refused IEEE-754 bits.
        bits: u64,
    },
    /// A set item duplicated the preceding item.
    DuplicateSetItem {
        /// Zero-based item index.
        index: u64,
    },
    /// A set item was smaller than its predecessor.
    NonCanonicalSetOrder {
        /// Zero-based item index.
        index: u64,
    },
    /// Caller-supplied cancellation was observed; no receipt was published.
    Cancelled {
        /// Canonical bytes absorbed before cancellation was observed.
        absorbed_bytes: u64,
    },
    /// A child field was declared without its expected-child binding.
    ChildBindingMissing {
        /// The unbound field.
        field: &'static str,
    },
    /// The supplied child identity does not match the parent-declared
    /// binding.
    ChildBindingMismatch {
        /// The bound field.
        field: &'static str,
        /// First mismatched dimension.
        what: &'static str,
    },
}

impl fmt::Display for CanonicalError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidLimits(reason) => {
                write!(f, "invalid canonical identity limits: {reason}")
            }
            Self::InvalidSchemaDescriptor(reason) => {
                write!(f, "invalid canonical identity schema: {reason}")
            }
            Self::LengthOverflow => f.write_str("canonical identity length arithmetic overflowed"),
            Self::LimitExceeded {
                kind,
                requested,
                limit,
            } => write!(
                f,
                "canonical identity {kind:?} request {requested} exceeds limit {limit}"
            ),
            Self::FieldOrder { expected, actual } => write!(
                f,
                "canonical identity expected field {expected}, received field {actual}"
            ),
            Self::FieldName => f.write_str("canonical identity field name does not match schema"),
            Self::WireType => f.write_str("canonical identity wire type does not match schema"),
            Self::Presence => f.write_str("canonical identity presence does not match schema"),
            Self::MissingFields { expected, actual } => write!(
                f,
                "canonical identity has {actual} fields; schema requires {expected}"
            ),
            Self::DeclaredLengthMismatch { declared, observed } => write!(
                f,
                "canonical identity declared {declared} bytes/items but observed {observed}"
            ),
            Self::NonFiniteFloat { bits } => {
                write!(
                    f,
                    "canonical identity refuses non-finite f64 bits 0x{bits:016x}"
                )
            }
            Self::DuplicateSetItem { index } => {
                write!(f, "canonical identity set item {index} is a duplicate")
            }
            Self::NonCanonicalSetOrder { index } => write!(
                f,
                "canonical identity set item {index} is not in canonical order"
            ),
            Self::Cancelled { absorbed_bytes } => write!(
                f,
                "canonical identity cancelled after absorbing {absorbed_bytes} bytes"
            ),
            Self::ChildBindingMissing { field } => write!(
                f,
                "child field `{field}` declares no expected child schema; bind it \
                 with FieldSpec::child_of or ordered_children_of"
            ),
            Self::ChildBindingMismatch { field, what } => write!(
                f,
                "child field `{field}` refuses this identity type: {what} does not \
                 match the parent-declared binding"
            ),
        }
    }
}

impl core::error::Error for CanonicalError {}

/// Caller-supplied cancellation checkpoint.
///
/// This leaf crate cannot depend on `fs-exec` because `fs-exec` already
/// depends on it. Downstream code adapts its `Cx` to this one-method trait.
pub trait CancellationProbe {
    /// Return true when construction must stop without publishing an ID.
    fn is_cancelled(&mut self) -> bool;
}

impl<F> CancellationProbe for F
where
    F: FnMut() -> bool,
{
    fn is_cancelled(&mut self) -> bool {
        self()
    }
}

/// Explicit probe for synchronous, non-cancellable call sites.
#[derive(Debug, Default, Clone, Copy)]
pub struct NeverCancel;

impl CancellationProbe for NeverCancel {
    fn is_cancelled(&mut self) -> bool {
        false
    }
}

struct CanonicalIdentityHeaderSource<'a> {
    role: IdentityRole,
    domain: &'a str,
    schema_name: &'a str,
    schema_id: [u8; 32],
    version: u32,
    context: &'a str,
    fields: &'a [FieldSpec],
}

#[allow(dead_code)]
fn classify_canonical_identity_header_fields(source: &CanonicalIdentityHeaderSource<'_>) {
    let CanonicalIdentityHeaderSource {
        role: _,
        domain: _,
        schema_name: _,
        schema_id: _,
        version: _,
        context: _,
        fields: _,
    } = source;
}

/// Transactional, bounded, streaming canonical identity encoder.
///
/// The encoder retains only two BLAKE3 states, fixed metadata, and counters:
/// it never buffers the canonical preimage. One hasher produces the typed
/// derive-key root and the other produces a plain [`ContentId`] for collision
/// adjudication. Every fallible operation consumes `self`; only [`finish`](Self::finish)
/// publishes either root.
pub struct CanonicalEncoder<I, C> {
    semantic_hasher: Blake3,
    preimage_hasher: Blake3,
    make_identity: fn(ContentHash) -> I,
    role: IdentityRole,
    schema_id: [u8; 32],
    limits: CanonicalLimits,
    cancellation: C,
    canonical_bytes: u64,
    next_field: u32,
    collection_items: u64,
}

impl<I, C> fmt::Debug for CanonicalEncoder<I, C> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("CanonicalEncoder")
            .field("role", &self.role)
            .field("schema_id", &"<typed-schema-id>")
            .field("limits", &self.limits)
            .field("canonical_bytes", &self.canonical_bytes)
            .field("next_field", &self.next_field)
            .field("collection_items", &self.collection_items)
            .finish_non_exhaustive()
    }
}

macro_rules! encoder_constructor {
    ($name:ident, $role:expr) => {
        impl<D, C> CanonicalEncoder<$name<D>, C>
        where
            D: CanonicalSchema,
            C: CancellationProbe,
        {
            /// Start an encoder for this exact role and static schema.
            ///
            /// Header/schema limits and cancellation are checked before a
            /// usable encoder is returned.
            pub fn new(limits: CanonicalLimits, cancellation: C) -> Result<Self, CanonicalError> {
                Self::new_internal::<D>($role, $name::<D>::from_digest, limits, cancellation)
            }
        }
    };
}

encoder_constructor!(SemanticId, IdentityRole::Semantic);
encoder_constructor!(WireContentId, IdentityRole::WireContent);
encoder_constructor!(EvidenceNodeId, IdentityRole::EvidenceNode);
encoder_constructor!(EntityId, IdentityRole::Entity);
encoder_constructor!(SourceByteId, IdentityRole::SourceBytes);
encoder_constructor!(SourceId, IdentityRole::Source);
encoder_constructor!(ModelId, IdentityRole::Model);
encoder_constructor!(CheckerId, IdentityRole::Checker);
encoder_constructor!(VerifierId, IdentityRole::Verifier);
encoder_constructor!(KeyPolicyId, IdentityRole::KeyPolicy);
encoder_constructor!(ProblemSemanticId, IdentityRole::ProblemSemantic);

impl<I, C> CanonicalEncoder<I, C>
where
    C: CancellationProbe,
{
    fn new_internal<D: CanonicalSchema>(
        role: IdentityRole,
        make_identity: fn(ContentHash) -> I,
        limits: CanonicalLimits,
        cancellation: C,
    ) -> Result<Self, CanonicalError> {
        validate_limits(limits)?;
        let mut encoder = Self {
            semantic_hasher: derive_key_hasher(CANONICAL_IDENTITY_HASH_DOMAIN),
            preimage_hasher: Blake3::new(),
            make_identity,
            role,
            schema_id: [0; 32],
            limits,
            cancellation,
            canonical_bytes: 0,
            next_field: 0,
            collection_items: 0,
        };
        encoder.checkpoint()?;
        encoder.validate_schema::<D>()?;
        let provisional_source = CanonicalIdentityHeaderSource {
            role,
            domain: D::DOMAIN,
            schema_name: D::NAME,
            schema_id: [0; 32],
            version: D::VERSION,
            context: D::CONTEXT,
            fields: D::FIELDS,
        };
        let header_bytes = encoder.canonical_header_len(&provisional_source)?;
        enforce_limit(
            LimitKind::CanonicalBytes,
            header_bytes,
            limits.max_canonical_bytes,
        )?;
        let schema_id = encoder.compute_schema_id::<D>()?;
        encoder.schema_id = schema_id;
        let source = CanonicalIdentityHeaderSource {
            role,
            domain: D::DOMAIN,
            schema_name: D::NAME,
            schema_id,
            version: D::VERSION,
            context: D::CONTEXT,
            fields: D::FIELDS,
        };
        encoder.write_header(&source)?;
        debug_assert_eq!(encoder.canonical_bytes, header_bytes);
        Ok(encoder)
    }

    fn write_header(
        &mut self,
        source: &CanonicalIdentityHeaderSource<'_>,
    ) -> Result<(), CanonicalError> {
        self.append(CANONICAL_MAGIC)?;
        self.append(&CANONICAL_FRAME_VERSION.to_le_bytes())?;
        self.append(&[source.role.tag(), FLOAT_POLICY_FINITE_EXACT_BITS])?;
        self.append_len_bytes(source.domain.as_bytes())?;
        self.append_len_bytes(source.schema_name.as_bytes())?;
        self.append(&source.schema_id)?;
        self.append(&source.version.to_le_bytes())?;
        self.append_len_bytes(source.context.as_bytes())?;
        let field_count =
            u32::try_from(source.fields.len()).map_err(|_| CanonicalError::LengthOverflow)?;
        self.append(&field_count.to_le_bytes())?;
        for (ordinal, field) in source.fields.iter().copied().enumerate() {
            let ordinal = u32::try_from(ordinal).map_err(|_| CanonicalError::LengthOverflow)?;
            self.append(&ordinal.to_le_bytes())?;
            self.append_len_bytes(field.name.as_bytes())?;
            self.append(&[field.wire_type.tag(), field.presence.tag()])?;
        }
        Ok(())
    }

    fn checkpoint(&mut self) -> Result<(), CanonicalError> {
        if self.cancellation.is_cancelled() {
            Err(CanonicalError::Cancelled {
                absorbed_bytes: self.canonical_bytes,
            })
        } else {
            Ok(())
        }
    }

    fn validate_schema<D: CanonicalSchema>(&mut self) -> Result<(), CanonicalError> {
        self.checkpoint()?;
        if D::DOMAIN.is_empty() || D::NAME.is_empty() || D::CONTEXT.is_empty() {
            return Err(CanonicalError::InvalidSchemaDescriptor(
                "domain, schema name, and context must be non-empty",
            ));
        }
        if D::VERSION == 0 {
            return Err(CanonicalError::InvalidSchemaDescriptor(
                "semantic version zero is reserved",
            ));
        }
        let field_count =
            u64::try_from(D::FIELDS.len()).map_err(|_| CanonicalError::LengthOverflow)?;
        enforce_limit(
            LimitKind::Fields,
            field_count,
            u64::from(self.limits.max_fields),
        )?;
        for descriptor in [D::DOMAIN, D::NAME, D::CONTEXT] {
            self.checkpoint()?;
            enforce_limit(
                LimitKind::FieldBytes,
                as_u64(descriptor.len())?,
                self.limits.max_field_bytes,
            )?;
        }
        for (index, field) in D::FIELDS.iter().copied().enumerate() {
            self.checkpoint()?;
            if field.name.is_empty() {
                return Err(CanonicalError::InvalidSchemaDescriptor(
                    "field names must be non-empty",
                ));
            }
            enforce_limit(
                LimitKind::FieldBytes,
                as_u64(field.name.len())?,
                self.limits.max_field_bytes,
            )?;
            for previous in &D::FIELDS[..index] {
                if self.compare_canonical_slices(previous.name.as_bytes(), field.name.as_bytes())?
                    == core::cmp::Ordering::Equal
                {
                    return Err(CanonicalError::InvalidSchemaDescriptor(
                        "field names must be unique",
                    ));
                }
            }
        }
        Ok(())
    }

    fn canonical_header_len(
        &mut self,
        source: &CanonicalIdentityHeaderSource<'_>,
    ) -> Result<u64, CanonicalError> {
        self.checkpoint()?;
        let mut total = checked_sum(&[
            as_u64(CANONICAL_MAGIC.len())?,
            u64::from(u32::BITS / 8),
            2,
            u64::from(u64::BITS / 8),
            as_u64(source.domain.len())?,
            u64::from(u64::BITS / 8),
            as_u64(source.schema_name.len())?,
            32,
            u64::from(u32::BITS / 8),
            u64::from(u64::BITS / 8),
            as_u64(source.context.len())?,
            u64::from(u32::BITS / 8),
        ])?;
        for field in source.fields {
            self.checkpoint()?;
            total = checked_add(
                total,
                checked_sum(&[
                    u64::from(u32::BITS / 8),
                    u64::from(u64::BITS / 8),
                    as_u64(field.name.len())?,
                    2,
                ])?,
            )?;
        }
        Ok(total)
    }

    fn auxiliary_update(
        &mut self,
        hasher: &mut Blake3,
        mut bytes: &[u8],
    ) -> Result<(), CanonicalError> {
        let stride = usize::try_from(self.limits.cancellation_poll_bytes)
            .map_err(|_| CanonicalError::LengthOverflow)?;
        while !bytes.is_empty() {
            self.checkpoint()?;
            let take = stride.min(bytes.len());
            let (chunk, remainder) = bytes.split_at(take);
            hasher.update(chunk);
            bytes = remainder;
        }
        Ok(())
    }

    fn compute_schema_id<D: CanonicalSchema>(&mut self) -> Result<[u8; 32], CanonicalError> {
        let source = SchemaDescriptorSource {
            domain: D::DOMAIN,
            name: D::NAME,
            version: D::VERSION,
            context: D::CONTEXT,
            fields: D::FIELDS,
        };
        let mut hasher = derive_key_hasher(SCHEMA_ID_HASH_DOMAIN);
        write_schema_descriptor(&source, |bytes| self.auxiliary_update(&mut hasher, bytes))?;
        Ok(*hasher.finalize().as_bytes())
    }

    fn compare_canonical_slices(
        &mut self,
        left: &[u8],
        right: &[u8],
    ) -> Result<core::cmp::Ordering, CanonicalError> {
        let stride = usize::try_from(self.limits.cancellation_poll_bytes)
            .map_err(|_| CanonicalError::LengthOverflow)?;
        let common_len = left.len().min(right.len());
        let mut offset = 0usize;
        self.checkpoint()?;
        while offset < common_len {
            self.checkpoint()?;
            let end = offset.saturating_add(stride).min(common_len);
            match left[offset..end].cmp(&right[offset..end]) {
                core::cmp::Ordering::Equal => offset = end,
                ordering => return Ok(ordering),
            }
        }
        Ok(left.len().cmp(&right.len()))
    }

    fn append(&mut self, mut bytes: &[u8]) -> Result<(), CanonicalError> {
        let length = as_u64(bytes.len())?;
        let requested = checked_add(self.canonical_bytes, length)?;
        enforce_limit(
            LimitKind::CanonicalBytes,
            requested,
            self.limits.max_canonical_bytes,
        )?;

        let stride = usize::try_from(self.limits.cancellation_poll_bytes)
            .map_err(|_| CanonicalError::LengthOverflow)?;
        while !bytes.is_empty() {
            self.checkpoint()?;
            let take = stride.min(bytes.len());
            let (chunk, remainder) = bytes.split_at(take);
            self.semantic_hasher.update(chunk);
            self.preimage_hasher.update(chunk);
            self.canonical_bytes = checked_add(self.canonical_bytes, as_u64(take)?)?;
            bytes = remainder;
        }
        Ok(())
    }

    fn append_len_bytes(&mut self, bytes: &[u8]) -> Result<(), CanonicalError> {
        self.append(&as_u64(bytes.len())?.to_le_bytes())?;
        self.append(bytes)
    }

    /// bead sj31i.52.10: a child field admits ONLY the parent-declared
    /// child role and complete schema identity.
    fn validate_child_binding<J: StrongIdentity>(spec: FieldSpec) -> Result<(), CanonicalError> {
        let Some(expected) = spec.child else {
            return Err(CanonicalError::ChildBindingMissing { field: spec.name });
        };
        expected
            .matches::<J>()
            .map_err(|what| CanonicalError::ChildBindingMismatch {
                field: spec.name,
                what,
            })
    }

    fn validate_field<D: CanonicalSchema>(
        &self,
        field: Field,
        wire_type: WireType,
        presence: Presence,
    ) -> Result<FieldSpec, CanonicalError> {
        if field.ordinal != self.next_field {
            return Err(CanonicalError::FieldOrder {
                expected: self.next_field,
                actual: field.ordinal,
            });
        }
        let index = usize::try_from(field.ordinal).map_err(|_| CanonicalError::LengthOverflow)?;
        let Some(expected) = D::FIELDS.get(index).copied() else {
            return Err(CanonicalError::FieldOrder {
                expected: self.next_field,
                actual: field.ordinal,
            });
        };
        if expected.name != field.name {
            return Err(CanonicalError::FieldName);
        }
        if expected.wire_type != wire_type {
            return Err(CanonicalError::WireType);
        }
        if expected.presence != presence {
            return Err(CanonicalError::Presence);
        }
        Ok(expected)
    }

    fn field_prefix_len(spec: FieldSpec) -> Result<u64, CanonicalError> {
        checked_sum(&[
            1,
            u64::from(u32::BITS / 8),
            u64::from(u64::BITS / 8),
            as_u64(spec.name.len())?,
            2,
        ])
    }

    fn ensure_additional(&self, additional: u64) -> Result<(), CanonicalError> {
        let requested = checked_add(self.canonical_bytes, additional)?;
        enforce_limit(
            LimitKind::CanonicalBytes,
            requested,
            self.limits.max_canonical_bytes,
        )
    }

    fn ensure_field_bytes(&self, requested: u64) -> Result<(), CanonicalError> {
        enforce_limit(
            LimitKind::FieldBytes,
            requested,
            self.limits.max_field_bytes,
        )
    }

    fn begin_field<D: CanonicalSchema>(
        &mut self,
        field: Field,
        wire_type: WireType,
        presence: Presence,
    ) -> Result<(), CanonicalError> {
        let spec = self.validate_field::<D>(field, wire_type, presence)?;
        self.append(&[FIELD_MARKER])?;
        self.append(&field.ordinal.to_le_bytes())?;
        self.append_len_bytes(spec.name.as_bytes())?;
        self.append(&[wire_type.tag(), presence.tag()])
    }

    fn complete_field(&mut self) -> Result<(), CanonicalError> {
        self.next_field = self
            .next_field
            .checked_add(1)
            .ok_or(CanonicalError::LengthOverflow)?;
        Ok(())
    }

    fn add_collection_items(&mut self, count: u64) -> Result<(), CanonicalError> {
        self.collection_items = checked_add(self.collection_items, count)?;
        Ok(())
    }
}

fn validate_limits(limits: CanonicalLimits) -> Result<(), CanonicalError> {
    if limits.cancellation_poll_bytes == 0 {
        return Err(CanonicalError::InvalidLimits(
            "cancellation_poll_bytes must be positive",
        ));
    }
    Ok(())
}

fn as_u64(value: usize) -> Result<u64, CanonicalError> {
    u64::try_from(value).map_err(|_| CanonicalError::LengthOverflow)
}

fn checked_add(left: u64, right: u64) -> Result<u64, CanonicalError> {
    left.checked_add(right)
        .ok_or(CanonicalError::LengthOverflow)
}

fn checked_sum(values: &[u64]) -> Result<u64, CanonicalError> {
    values
        .iter()
        .try_fold(0u64, |sum, value| checked_add(sum, *value))
}

fn enforce_limit(kind: LimitKind, requested: u64, limit: u64) -> Result<(), CanonicalError> {
    if requested > limit {
        Err(CanonicalError::LimitExceeded {
            kind,
            requested,
            limit,
        })
    } else {
        Ok(())
    }
}

impl<I, C> CanonicalEncoder<I, C>
where
    I: StrongIdentity,
    C: CancellationProbe,
{
    /// Encode one required exact UTF-8 field. No Unicode normalization,
    /// locale transform, JSON rendering, or display formatting is applied.
    pub fn utf8(mut self, field: Field, value: &str) -> Result<Self, CanonicalError> {
        let spec = self.validate_field::<I::Schema>(field, WireType::Utf8, Presence::Required)?;
        let length = as_u64(value.len())?;
        self.ensure_field_bytes(length)?;
        self.ensure_additional(checked_add(
            Self::field_prefix_len(spec)?,
            checked_add(u64::from(u64::BITS / 8), length)?,
        )?)?;
        self.begin_field::<I::Schema>(field, WireType::Utf8, Presence::Required)?;
        self.append_len_bytes(value.as_bytes())?;
        self.complete_field()?;
        Ok(self)
    }

    /// Encode one required byte slice.
    pub fn bytes(self, field: Field, value: &[u8]) -> Result<Self, CanonicalError> {
        if value.is_empty() {
            return self.bytes_stream(field, 0, core::iter::empty());
        }
        self.bytes_stream(field, as_u64(value.len())?, core::iter::once(value))
    }

    /// Encode a required byte field from chunks without retaining the field.
    ///
    /// `declared_len` is admitted against field and complete-frame budgets
    /// before `chunks` is read. Too few or too many bytes consume and refuse
    /// the encoder, so no partial identity can be finished.
    pub fn bytes_stream<'a, T>(
        mut self,
        field: Field,
        declared_len: u64,
        chunks: T,
    ) -> Result<Self, CanonicalError>
    where
        T: IntoIterator<Item = &'a [u8]>,
    {
        let spec = self.validate_field::<I::Schema>(field, WireType::Bytes, Presence::Required)?;
        self.ensure_field_bytes(declared_len)?;
        self.ensure_additional(checked_add(
            Self::field_prefix_len(spec)?,
            checked_add(u64::from(u64::BITS / 8), declared_len)?,
        )?)?;
        self.begin_field::<I::Schema>(field, WireType::Bytes, Presence::Required)?;
        self.append(&declared_len.to_le_bytes())?;
        let mut observed = 0u64;
        let mut chunk_count = 0u64;
        for chunk in chunks {
            self.checkpoint()?;
            chunk_count = checked_add(chunk_count, 1)?;
            enforce_limit(
                LimitKind::StreamChunks,
                chunk_count,
                self.limits.max_collection_items,
            )?;
            observed = checked_add(observed, as_u64(chunk.len())?)?;
            if observed > declared_len {
                return Err(CanonicalError::DeclaredLengthMismatch {
                    declared: declared_len,
                    observed,
                });
            }
            self.append(chunk)?;
        }
        if observed != declared_len {
            return Err(CanonicalError::DeclaredLengthMismatch {
                declared: declared_len,
                observed,
            });
        }
        self.complete_field()?;
        Ok(self)
    }

    /// Encode one required little-endian `u64`.
    pub fn u64(mut self, field: Field, value: u64) -> Result<Self, CanonicalError> {
        self.fixed_field(field, WireType::U64, &value.to_le_bytes())?;
        Ok(self)
    }

    /// Encode one required little-endian `i64`.
    pub fn i64(mut self, field: Field, value: i64) -> Result<Self, CanonicalError> {
        self.fixed_field(field, WireType::I64, &value.to_le_bytes())?;
        Ok(self)
    }

    /// Encode one required boolean as exactly zero or one.
    pub fn flag(mut self, field: Field, value: bool) -> Result<Self, CanonicalError> {
        self.fixed_field(field, WireType::Bool, &[u8::from(value)])?;
        Ok(self)
    }

    /// Encode one finite `f64` by its exact IEEE-754 little-endian bits.
    ///
    /// `+0.0` and `-0.0` are intentionally distinct. A schema that normalizes
    /// signed zero must do so before this call and use its own schema version.
    /// Every NaN payload and both infinities refuse before field bytes mutate
    /// the hash state.
    pub fn finite_f64(mut self, field: Field, value: f64) -> Result<Self, CanonicalError> {
        if !value.is_finite() {
            return Err(CanonicalError::NonFiniteFloat {
                bits: value.to_bits(),
            });
        }
        self.fixed_field(field, WireType::FiniteF64, &value.to_bits().to_le_bytes())?;
        Ok(self)
    }

    /// Encode an explicitly optional byte field.
    ///
    /// `None`, `Some(&[])`, and an absent schema field are three different
    /// states; the last is refused at [`finish`](Self::finish).
    pub fn optional_bytes(
        mut self,
        field: Field,
        value: Option<&[u8]>,
    ) -> Result<Self, CanonicalError> {
        let spec = self.validate_field::<I::Schema>(field, WireType::Bytes, Presence::Optional)?;
        let value_len = value.map_or(0, <[u8]>::len);
        let value_len = as_u64(value_len)?;
        self.ensure_field_bytes(value_len)?;
        let payload_len = if value.is_some() {
            checked_sum(&[1, u64::from(u64::BITS / 8), value_len])?
        } else {
            1
        };
        self.ensure_additional(checked_add(Self::field_prefix_len(spec)?, payload_len)?)?;
        self.begin_field::<I::Schema>(field, WireType::Bytes, Presence::Optional)?;
        match value {
            None => self.append(&[0])?,
            Some(bytes) => {
                self.append(&[1])?;
                self.append_len_bytes(bytes)?;
            }
        }
        self.complete_field()?;
        Ok(self)
    }

    /// Encode a numeric variant tag and exact byte payload.
    pub fn variant(
        mut self,
        field: Field,
        variant: u32,
        payload: &[u8],
    ) -> Result<Self, CanonicalError> {
        let spec =
            self.validate_field::<I::Schema>(field, WireType::Variant, Presence::Required)?;
        let payload_len = as_u64(payload.len())?;
        self.ensure_field_bytes(payload_len)?;
        let encoded_payload = checked_sum(&[
            u64::from(u32::BITS / 8),
            u64::from(u64::BITS / 8),
            payload_len,
        ])?;
        self.ensure_additional(checked_add(Self::field_prefix_len(spec)?, encoded_payload)?)?;
        self.begin_field::<I::Schema>(field, WireType::Variant, Presence::Required)?;
        self.append(&variant.to_le_bytes())?;
        self.append_len_bytes(payload)?;
        self.complete_field()?;
        Ok(self)
    }

    /// Encode a caller-ordered byte sequence.
    ///
    /// Sequence order is semantic. The encoder never sorts or allocates.
    pub fn ordered_bytes<'a, T>(
        mut self,
        field: Field,
        declared_count: u64,
        values: T,
    ) -> Result<Self, CanonicalError>
    where
        T: IntoIterator<Item = &'a [u8]>,
    {
        let spec =
            self.validate_field::<I::Schema>(field, WireType::OrderedBytes, Presence::Required)?;
        enforce_limit(
            LimitKind::CollectionItems,
            declared_count,
            self.limits.max_collection_items,
        )?;
        self.ensure_additional(checked_add(
            Self::field_prefix_len(spec)?,
            u64::from(u64::BITS / 8),
        )?)?;
        self.begin_field::<I::Schema>(field, WireType::OrderedBytes, Presence::Required)?;
        self.append(&declared_count.to_le_bytes())?;
        let mut observed = 0u64;
        let mut field_payload = u64::from(u64::BITS / 8);
        for value in values {
            observed = checked_add(observed, 1)?;
            if observed > declared_count {
                return Err(CanonicalError::DeclaredLengthMismatch {
                    declared: declared_count,
                    observed,
                });
            }
            let value_len = as_u64(value.len())?;
            self.ensure_field_bytes(value_len)?;
            field_payload = checked_add(
                field_payload,
                checked_add(u64::from(u64::BITS / 8), value_len)?,
            )?;
            self.ensure_field_bytes(field_payload)?;
            self.ensure_additional(checked_add(u64::from(u64::BITS / 8), value_len)?)?;
            self.append_len_bytes(value)?;
        }
        if observed != declared_count {
            return Err(CanonicalError::DeclaredLengthMismatch {
                declared: declared_count,
                observed,
            });
        }
        self.add_collection_items(observed)?;
        self.complete_field()?;
        Ok(self)
    }

    /// Encode a strictly lexicographically increasing, duplicate-free set.
    ///
    /// The core refuses unsorted input instead of secretly allocating and
    /// guessing a domain's collation rules.
    pub fn canonical_set<'a, T>(
        mut self,
        field: Field,
        declared_count: u64,
        values: T,
    ) -> Result<Self, CanonicalError>
    where
        T: IntoIterator<Item = &'a [u8]>,
    {
        let spec =
            self.validate_field::<I::Schema>(field, WireType::CanonicalSet, Presence::Required)?;
        enforce_limit(
            LimitKind::CollectionItems,
            declared_count,
            self.limits.max_collection_items,
        )?;
        self.ensure_additional(checked_add(
            Self::field_prefix_len(spec)?,
            u64::from(u64::BITS / 8),
        )?)?;
        self.begin_field::<I::Schema>(field, WireType::CanonicalSet, Presence::Required)?;
        self.append(&declared_count.to_le_bytes())?;
        let mut observed = 0u64;
        let mut field_payload = u64::from(u64::BITS / 8);
        let mut previous: Option<&'a [u8]> = None;
        for value in values {
            observed = checked_add(observed, 1)?;
            if observed > declared_count {
                return Err(CanonicalError::DeclaredLengthMismatch {
                    declared: declared_count,
                    observed,
                });
            }
            let value_len = as_u64(value.len())?;
            self.ensure_field_bytes(value_len)?;
            let next_field_payload = checked_add(
                field_payload,
                checked_add(u64::from(u64::BITS / 8), value_len)?,
            )?;
            self.ensure_field_bytes(next_field_payload)?;
            self.ensure_additional(checked_add(u64::from(u64::BITS / 8), value_len)?)?;
            // Admit the item before scanning a hostile equal prefix for order.
            if let Some(before) = previous {
                match self.compare_canonical_slices(before, value)? {
                    core::cmp::Ordering::Equal => {
                        return Err(CanonicalError::DuplicateSetItem {
                            index: observed - 1,
                        });
                    }
                    core::cmp::Ordering::Greater => {
                        return Err(CanonicalError::NonCanonicalSetOrder {
                            index: observed - 1,
                        });
                    }
                    core::cmp::Ordering::Less => {}
                }
            }
            self.append_len_bytes(value)?;
            field_payload = next_field_payload;
            previous = Some(value);
        }
        if observed != declared_count {
            return Err(CanonicalError::DeclaredLengthMismatch {
                declared: declared_count,
                observed,
            });
        }
        self.add_collection_items(observed)?;
        self.complete_field()?;
        Ok(self)
    }

    /// Encode one full typed child identity, including role and schema.
    pub fn child<J>(mut self, field: Field, child: J) -> Result<Self, CanonicalError>
    where
        J: StrongIdentity,
    {
        let spec = self.validate_field::<I::Schema>(field, WireType::Child, Presence::Required)?;
        Self::validate_child_binding::<J>(spec)?;
        self.validate_schema::<J::Schema>()?;
        let child_schema_id = self.compute_schema_id::<J::Schema>()?;
        let child_len = typed_child_len::<J>()?;
        self.ensure_field_bytes(child_len)?;
        self.ensure_additional(checked_add(Self::field_prefix_len(spec)?, child_len)?)?;
        self.begin_field::<I::Schema>(field, WireType::Child, Presence::Required)?;
        self.append_typed_child(child, child_schema_id)?;
        self.complete_field()?;
        Ok(self)
    }

    /// Encode an ordered sequence of children of one exact role/schema type.
    pub fn ordered_children<J, T>(
        mut self,
        field: Field,
        declared_count: u64,
        children: T,
    ) -> Result<Self, CanonicalError>
    where
        J: StrongIdentity,
        T: IntoIterator<Item = J>,
    {
        let spec =
            self.validate_field::<I::Schema>(field, WireType::OrderedChildren, Presence::Required)?;
        Self::validate_child_binding::<J>(spec)?;
        enforce_limit(
            LimitKind::CollectionItems,
            declared_count,
            self.limits.max_collection_items,
        )?;
        self.validate_schema::<J::Schema>()?;
        let child_schema_id = self.compute_schema_id::<J::Schema>()?;
        let descriptor_len = typed_child_descriptor_len::<J>()?;
        let payload_len = checked_sum(&[
            u64::from(u64::BITS / 8),
            descriptor_len,
            declared_count
                .checked_mul(32)
                .ok_or(CanonicalError::LengthOverflow)?,
        ])?;
        self.ensure_field_bytes(payload_len)?;
        self.ensure_additional(checked_add(Self::field_prefix_len(spec)?, payload_len)?)?;
        self.begin_field::<I::Schema>(field, WireType::OrderedChildren, Presence::Required)?;
        self.append(&declared_count.to_le_bytes())?;
        self.append_typed_child_descriptor::<J>(child_schema_id)?;
        let mut observed = 0u64;
        for child in children {
            observed = checked_add(observed, 1)?;
            if observed > declared_count {
                return Err(CanonicalError::DeclaredLengthMismatch {
                    declared: declared_count,
                    observed,
                });
            }
            self.append(child.as_bytes())?;
        }
        if observed != declared_count {
            return Err(CanonicalError::DeclaredLengthMismatch {
                declared: declared_count,
                observed,
            });
        }
        self.add_collection_items(observed)?;
        self.complete_field()?;
        Ok(self)
    }

    fn fixed_field(
        &mut self,
        field: Field,
        wire_type: WireType,
        bytes: &[u8],
    ) -> Result<(), CanonicalError> {
        let spec = self.validate_field::<I::Schema>(field, wire_type, Presence::Required)?;
        let length = as_u64(bytes.len())?;
        self.ensure_field_bytes(length)?;
        self.ensure_additional(checked_add(Self::field_prefix_len(spec)?, length)?)?;
        self.begin_field::<I::Schema>(field, wire_type, Presence::Required)?;
        self.append(bytes)?;
        self.complete_field()
    }

    fn append_typed_child<J: StrongIdentity>(
        &mut self,
        child: J,
        child_schema_id: [u8; 32],
    ) -> Result<(), CanonicalError> {
        self.append_typed_child_descriptor::<J>(child_schema_id)?;
        self.append(child.as_bytes())
    }

    fn append_typed_child_descriptor<J: StrongIdentity>(
        &mut self,
        child_schema_id: [u8; 32],
    ) -> Result<(), CanonicalError> {
        self.append(&[J::ROLE.tag()])?;
        self.append_len_bytes(J::Schema::DOMAIN.as_bytes())?;
        self.append_len_bytes(J::Schema::NAME.as_bytes())?;
        self.append(&child_schema_id)?;
        self.append(&J::Schema::VERSION.to_le_bytes())?;
        self.append_len_bytes(J::Schema::CONTEXT.as_bytes())
    }

    /// Finish the exact declared field set and publish both roots.
    ///
    /// The final cancellation checkpoint is the publication linearization
    /// point. If it refuses, neither the typed root nor the preimage root is
    /// returned.
    pub fn finish(mut self) -> Result<IdentityReceipt<I>, CanonicalError> {
        let expected =
            u32::try_from(I::Schema::FIELDS.len()).map_err(|_| CanonicalError::LengthOverflow)?;
        if self.next_field != expected {
            return Err(CanonicalError::MissingFields {
                expected,
                actual: self.next_field,
            });
        }
        self.ensure_additional(checked_sum(&[1, u64::from(u32::BITS / 8)])?)?;
        self.append(&[END_MARKER])?;
        self.append(&self.next_field.to_le_bytes())?;
        self.checkpoint()?;
        let id = (self.make_identity)(self.semantic_hasher.finalize());
        let canonical_preimage = ContentId(self.preimage_hasher.finalize());
        Ok(IdentityReceipt {
            id,
            canonical_preimage,
            schema_id: self.schema_id,
            canonical_bytes: self.canonical_bytes,
            field_count: self.next_field,
            collection_items: self.collection_items,
            limits: self.limits,
        })
    }
}

fn typed_child_descriptor_len<J: StrongIdentity>() -> Result<u64, CanonicalError> {
    checked_sum(&[
        1,
        u64::from(u64::BITS / 8),
        as_u64(J::Schema::DOMAIN.len())?,
        u64::from(u64::BITS / 8),
        as_u64(J::Schema::NAME.len())?,
        32,
        u64::from(u32::BITS / 8),
        u64::from(u64::BITS / 8),
        as_u64(J::Schema::CONTEXT.len())?,
    ])
}

fn typed_child_len<J: StrongIdentity>() -> Result<u64, CanonicalError> {
    checked_add(typed_child_descriptor_len::<J>()?, 32)
}

/// Successfully published typed identity plus its exact canonical-frame root.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct IdentityReceipt<I: StrongIdentity> {
    id: I,
    canonical_preimage: ContentId,
    schema_id: [u8; 32],
    canonical_bytes: u64,
    field_count: u32,
    collection_items: u64,
    limits: CanonicalLimits,
}

impl<I: StrongIdentity> IdentityReceipt<I> {
    /// Typed role/schema-specific identity.
    #[must_use]
    pub const fn id(self) -> I {
        self.id
    }

    /// Plain BLAKE3 root of the complete canonical frame.
    #[must_use]
    pub const fn canonical_preimage(self) -> ContentId {
        self.canonical_preimage
    }

    /// Exact number of frame bytes absorbed.
    #[must_use]
    pub const fn canonical_bytes(self) -> u64 {
        self.canonical_bytes
    }

    /// Exact number of encoded top-level fields.
    #[must_use]
    pub const fn field_count(self) -> u32 {
        self.field_count
    }

    /// Total successfully encoded collection items.
    #[must_use]
    pub const fn collection_items(self) -> u64 {
        self.collection_items
    }

    /// Admission budgets used by the producer. They are evidence metadata,
    /// not hash inputs.
    #[must_use]
    pub const fn limits(self) -> CanonicalLimits {
        self.limits
    }

    /// Identity of the complete static schema descriptor.
    #[must_use]
    pub fn schema_id(self) -> SchemaId<I::Schema> {
        SchemaId::from_digest(ContentHash(self.schema_id))
    }

    /// Fixed-size, payload-free audit record for an unanchored identity.
    #[must_use]
    pub fn audit_record(self) -> IdentityAuditRecord {
        IdentityAuditRecord::from_receipt(self)
    }
}

/// Trust state retained by bounded identity audit records.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum TrustState {
    /// Digest/semantic consistency only; no external anchor was presented.
    Unanchored,
    /// External data was presented but has not been verified.
    Presented,
    /// A verifier capability accepted the presentation; policy admission is
    /// still separate.
    Verified,
    /// A separate admission capability accepted the verified authority.
    Admitted,
}

/// Explicit boundary on what a receipt or authority record does not prove.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum NoClaimState {
    /// Digest and canonical semantics only; external trust is still required.
    ExternalTrustRequired,
    /// Authority state does not prove scientific/model correctness.
    ScientificCorrectnessNotProven,
}

/// Fixed-size, payload-free identity record suitable for bounded logging.
///
/// It never retains source payloads, canonical bytes, signatures, JSON, debug
/// text, hostnames, or clocks.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct IdentityAuditRecord {
    id: [u8; 32],
    canonical_preimage: ContentId,
    role: IdentityRole,
    domain: &'static str,
    schema_name: &'static str,
    schema_id: [u8; 32],
    version: u32,
    context: &'static str,
    canonical_bytes: u64,
    field_count: u32,
    collection_items: u64,
    limits: CanonicalLimits,
    trust: TrustState,
    anchor: Option<ContentId>,
    verifier: Option<[u8; 32]>,
    key_policy: Option<[u8; 32]>,
    no_claim: NoClaimState,
}

impl IdentityAuditRecord {
    fn from_receipt<I: StrongIdentity>(receipt: IdentityReceipt<I>) -> Self {
        Self {
            id: *receipt.id.as_bytes(),
            canonical_preimage: receipt.canonical_preimage,
            role: I::ROLE,
            domain: I::Schema::DOMAIN,
            schema_name: I::Schema::NAME,
            schema_id: receipt.schema_id,
            version: I::Schema::VERSION,
            context: I::Schema::CONTEXT,
            canonical_bytes: receipt.canonical_bytes,
            field_count: receipt.field_count,
            collection_items: receipt.collection_items,
            limits: receipt.limits,
            trust: TrustState::Unanchored,
            anchor: None,
            verifier: None,
            key_policy: None,
            no_claim: NoClaimState::ExternalTrustRequired,
        }
    }

    /// Typed digest bytes; `role`, domain, and schema must travel with them.
    #[must_use]
    pub const fn id(self) -> [u8; 32] {
        self.id
    }

    /// Plain root of the complete canonical frame.
    #[must_use]
    pub const fn canonical_preimage(self) -> ContentId {
        self.canonical_preimage
    }

    /// Non-interchangeable identity role.
    #[must_use]
    pub const fn role(self) -> IdentityRole {
        self.role
    }

    /// Registered static domain.
    #[must_use]
    pub const fn domain(self) -> &'static str {
        self.domain
    }

    /// Registered static schema name.
    #[must_use]
    pub const fn schema_name(self) -> &'static str {
        self.schema_name
    }

    /// Schema descriptor digest bytes.
    #[must_use]
    pub const fn schema_id(self) -> [u8; 32] {
        self.schema_id
    }

    /// Semantic schema version.
    #[must_use]
    pub const fn version(self) -> u32 {
        self.version
    }

    /// Static purpose/context.
    #[must_use]
    pub const fn context(self) -> &'static str {
        self.context
    }

    /// Complete canonical frame size.
    #[must_use]
    pub const fn canonical_bytes(self) -> u64 {
        self.canonical_bytes
    }

    /// Encoded field count.
    #[must_use]
    pub const fn field_count(self) -> u32 {
        self.field_count
    }

    /// Encoded collection item count.
    #[must_use]
    pub const fn collection_items(self) -> u64 {
        self.collection_items
    }

    /// Producer admission budgets.
    #[must_use]
    pub const fn limits(self) -> CanonicalLimits {
        self.limits
    }

    /// Trust state; presence alone is never admitted trust.
    #[must_use]
    pub const fn trust(self) -> TrustState {
        self.trust
    }

    /// External anchor bytes, present only when authority data was supplied.
    #[must_use]
    pub const fn anchor(self) -> Option<ContentId> {
        self.anchor
    }

    /// Verifier ID bytes, present only after an authority reference exists.
    #[must_use]
    pub const fn verifier(self) -> Option<[u8; 32]> {
        self.verifier
    }

    /// Key-policy ID bytes, present only after an authority reference exists.
    #[must_use]
    pub const fn key_policy(self) -> Option<[u8; 32]> {
        self.key_policy
    }

    /// Explicit no-claim boundary.
    #[must_use]
    pub const fn no_claim(self) -> NoClaimState {
        self.no_claim
    }
}

/// One retained observation used to adjudicate a claimed typed identity.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct ByteObservation {
    content_id: ContentId,
    length: u64,
}

impl ByteObservation {
    /// Construct an observation from an independently retained byte root and
    /// exact length. This data is untrusted until adjudicated.
    #[must_use]
    pub const fn new(content_id: ContentId, length: u64) -> Self {
        Self { content_id, length }
    }

    /// Retained byte root.
    #[must_use]
    pub const fn content_id(self) -> ContentId {
        self.content_id
    }

    /// Retained exact byte length.
    #[must_use]
    pub const fn length(self) -> u64 {
        self.length
    }
}

/// A typed identity presented with its independent canonical-byte observation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ObservedIdentity<I: StrongIdentity> {
    id: I,
    bytes: ByteObservation,
}

impl<I: StrongIdentity> ObservedIdentity<I> {
    /// Use a producer receipt as an observation.
    #[must_use]
    pub const fn from_receipt(receipt: IdentityReceipt<I>) -> Self {
        Self {
            id: receipt.id,
            bytes: ByteObservation::new(receipt.canonical_preimage, receipt.canonical_bytes),
        }
    }

    /// Present parsed/untrusted retained data for adjudication.
    #[must_use]
    pub const fn presented(id: I, bytes: ByteObservation) -> Self {
        Self { id, bytes }
    }

    /// Claimed typed identity.
    #[must_use]
    pub const fn id(self) -> I {
        self.id
    }

    /// Independent canonical-byte observation.
    #[must_use]
    pub const fn bytes(self) -> ByteObservation {
        self.bytes
    }
}

/// Typed refusal for one claimed ID backed by different byte observations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SameIdDifferentBytes<I: StrongIdentity> {
    id: I,
    first: ByteObservation,
    second: ByteObservation,
}

impl<I: StrongIdentity> SameIdDifferentBytes<I> {
    /// Refused typed ID.
    #[must_use]
    pub const fn id(self) -> I {
        self.id
    }

    /// First observation; it is not privileged over the second.
    #[must_use]
    pub const fn first(self) -> ByteObservation {
        self.first
    }

    /// Second observation; it is not privileged over the first.
    #[must_use]
    pub const fn second(self) -> ByteObservation {
        self.second
    }
}

/// Result of comparing two observations in one exact typed namespace.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IdentityAdjudication<I: StrongIdentity> {
    /// Typed IDs differ.
    DistinctIds,
    /// Typed ID, byte root, and byte length all match.
    SameObservation,
    /// One typed ID was presented for different retained byte observations.
    Refused(SameIdDifferentBytes<I>),
}

/// Compare two observations without first-wins/last-wins behavior.
///
/// The refusal path relies on independently retaining byte roots and lengths.
/// It cannot detect a collision after all distinguishing observations were
/// discarded, nor can finite testing prove BLAKE3 collision resistance.
#[must_use]
pub fn adjudicate<I: StrongIdentity>(
    first: ObservedIdentity<I>,
    second: ObservedIdentity<I>,
) -> IdentityAdjudication<I> {
    if first.id != second.id {
        IdentityAdjudication::DistinctIds
    } else if first.bytes == second.bytes {
        IdentityAdjudication::SameObservation
    } else {
        IdentityAdjudication::Refused(SameIdDifferentBytes {
            id: first.id,
            first: first.bytes,
            second: second.bytes,
        })
    }
}

/// Presented external anchor data. Its presence is not verification or trust.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct ExternalAnchorRef(ContentId);

impl ExternalAnchorRef {
    /// Mark exact external anchor bytes as presented, without trusting them.
    #[must_use]
    pub const fn presented(content_id: ContentId) -> Self {
        Self(content_id)
    }

    /// Presented anchor content ID.
    #[must_use]
    pub const fn content_id(self) -> ContentId {
        self.0
    }
}

/// Authority typestate: data has merely been presented.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Presented;

/// Authority typestate: an injected verifier accepted the exact presentation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Verified;

/// Authority typestate: a separate policy capability admitted the verifier's
/// decision for the exact subject/policy/context.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Admitted;

/// State marker for [`AuthorityRef`].
pub trait AuthorityState: Copy + 'static {
    /// Runtime log state corresponding to this typestate.
    const TRUST_STATE: TrustState;
}

impl AuthorityState for Presented {
    const TRUST_STATE: TrustState = TrustState::Presented;
}
impl AuthorityState for Verified {
    const TRUST_STATE: TrustState = TrustState::Verified;
}
impl AuthorityState for Admitted {
    const TRUST_STATE: TrustState = TrustState::Admitted;
}

/// Explicit authority data for one typed subject.
///
/// No `Deref` or conversion to the subject is implemented. The state, exact
/// verifier ID, and exact key-policy ID must remain visible at decision sites.
pub struct AuthorityRef<I, V, P, S>
where
    I: StrongIdentity,
    V: CanonicalSchema,
    P: CanonicalSchema,
    S: AuthorityState,
{
    receipt: IdentityReceipt<I>,
    anchor: ExternalAnchorRef,
    verifier: VerifierId<V>,
    key_policy: KeyPolicyId<P>,
    state: PhantomData<fn() -> S>,
}

impl<I, V, P, S> fmt::Debug for AuthorityRef<I, V, P, S>
where
    I: StrongIdentity,
    V: CanonicalSchema,
    P: CanonicalSchema,
    S: AuthorityState,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("AuthorityRef")
            .field("receipt", &self.receipt)
            .field("anchor", &self.anchor)
            .field("verifier", &self.verifier)
            .field("key_policy", &self.key_policy)
            .field("trust", &S::TRUST_STATE)
            .finish()
    }
}

impl<I, V, P, S> PartialEq for AuthorityRef<I, V, P, S>
where
    I: StrongIdentity,
    V: CanonicalSchema,
    P: CanonicalSchema,
    S: AuthorityState,
{
    fn eq(&self, other: &Self) -> bool {
        self.receipt == other.receipt
            && self.anchor == other.anchor
            && self.verifier == other.verifier
            && self.key_policy == other.key_policy
    }
}

impl<I, V, P, S> Eq for AuthorityRef<I, V, P, S>
where
    I: StrongIdentity,
    V: CanonicalSchema,
    P: CanonicalSchema,
    S: AuthorityState,
{
}

impl<I, V, P, S> AuthorityRef<I, V, P, S>
where
    I: StrongIdentity,
    V: CanonicalSchema,
    P: CanonicalSchema,
    S: AuthorityState,
{
    /// Exact subject receipt.
    #[must_use]
    pub const fn receipt(&self) -> IdentityReceipt<I> {
        self.receipt
    }

    /// Presented external anchor.
    #[must_use]
    pub const fn anchor(&self) -> ExternalAnchorRef {
        self.anchor
    }

    /// Exact verifier identity.
    #[must_use]
    pub const fn verifier(&self) -> VerifierId<V> {
        self.verifier
    }

    /// Exact key-policy identity.
    #[must_use]
    pub const fn key_policy(&self) -> KeyPolicyId<P> {
        self.key_policy
    }

    /// Runtime state corresponding to typestate `S`.
    #[must_use]
    pub const fn trust_state(&self) -> TrustState {
        S::TRUST_STATE
    }

    /// Fixed-size, payload-free audit record retaining trust and verifier data.
    #[must_use]
    pub fn audit_record(&self) -> IdentityAuditRecord {
        let mut record = self.receipt.audit_record();
        record.trust = S::TRUST_STATE;
        record.anchor = Some(self.anchor.content_id());
        record.verifier = Some(*self.verifier.as_bytes());
        record.key_policy = Some(*self.key_policy.as_bytes());
        record.no_claim = if S::TRUST_STATE == TrustState::Admitted {
            NoClaimState::ScientificCorrectnessNotProven
        } else {
            NoClaimState::ExternalTrustRequired
        };
        record
    }
}

impl<I, V, P> AuthorityRef<I, V, P, Presented>
where
    I: StrongIdentity,
    V: CanonicalSchema,
    P: CanonicalSchema,
{
    /// Present external anchor/verifier/policy data. This always yields the
    /// untrusted [`Presented`] state.
    #[must_use]
    pub const fn present(
        receipt: IdentityReceipt<I>,
        anchor: ExternalAnchorRef,
        verifier: VerifierId<V>,
        key_policy: KeyPolicyId<P>,
    ) -> Self {
        Self {
            receipt,
            anchor,
            verifier,
            key_policy,
            state: PhantomData,
        }
    }

    /// Ask an injected verifier capability to accept the exact presentation.
    /// The presentation is consumed on both success and refusal.
    pub fn verify<A>(self, capability: &A) -> Result<AuthorityRef<I, V, P, Verified>, A::Error>
    where
        A: AuthorityVerifier<I, V, P>,
    {
        capability.verify(&self)?;
        Ok(AuthorityRef {
            receipt: self.receipt,
            anchor: self.anchor,
            verifier: self.verifier,
            key_policy: self.key_policy,
            state: PhantomData,
        })
    }
}

impl<I, V, P> AuthorityRef<I, V, P, Verified>
where
    I: StrongIdentity,
    V: CanonicalSchema,
    P: CanonicalSchema,
{
    /// Ask a separate admission capability to trust this verifier decision for
    /// the exact subject, anchor, key policy, and context.
    pub fn admit<A>(self, capability: &A) -> Result<AuthorityRef<I, V, P, Admitted>, A::Error>
    where
        A: AuthorityAdmitter<I, V, P>,
    {
        capability.admit(&self)?;
        Ok(AuthorityRef {
            receipt: self.receipt,
            anchor: self.anchor,
            verifier: self.verifier,
            key_policy: self.key_policy,
            state: PhantomData,
        })
    }
}

/// Injected capability that validates presented external evidence.
pub trait AuthorityVerifier<I, V, P>
where
    I: StrongIdentity,
    V: CanonicalSchema,
    P: CanonicalSchema,
{
    /// Structured verifier refusal.
    type Error;
    /// Verify the exact subject, canonical preimage, anchor, verifier ID, and
    /// key-policy ID. A successful return does not itself admit the verifier.
    fn verify(&self, presented: &AuthorityRef<I, V, P, Presented>) -> Result<(), Self::Error>;
}

/// Separate policy capability that admits a verified authority decision.
pub trait AuthorityAdmitter<I, V, P>
where
    I: StrongIdentity,
    V: CanonicalSchema,
    P: CanonicalSchema,
{
    /// Structured policy refusal.
    type Error;
    /// Admit the exact verified subject/verifier/key-policy/context binding.
    fn admit(&self, verified: &AuthorityRef<I, V, P, Verified>) -> Result<(), Self::Error>;
}

/// Quarantined legacy identity types. They deliberately have no conversion,
/// widening, equality bridge, or child-identity implementation.
pub mod legacy {
    /// Exact historical FNV-1a `u64` provenance value.
    #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
    pub struct LegacyProvenanceV1(u64);

    impl LegacyProvenanceV1 {
        /// Retain an exact historical value without claiming strong identity.
        #[must_use]
        pub const fn new(value: u64) -> Self {
            Self(value)
        }

        /// Exact legacy value for replay/crosswalk lookup only.
        #[must_use]
        pub const fn value(self) -> u64 {
            self.0
        }
    }
}
