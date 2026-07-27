//! Conformance tests for schema-typed canonical identities, bounded encoding,
//! collision refusal, and explicit authority-state transitions.

use std::cell::Cell;
use std::rc::Rc;

use fs_blake3::identity::{
    AuthorityAdmitter, AuthorityRef, AuthorityVerifier, ByteObservation, CANONICAL_FRAME_VERSION,
    CANONICAL_IDENTITY_HASH_DOMAIN, CanonicalEncoder, CanonicalError, CanonicalLimits,
    CanonicalSchema, CheckerId, ChildSpec, ContentId, EntityId, EvidenceNodeId, ExternalAnchorRef,
    Field, FieldSpec, IdentityAdjudication, IdentityRole, KeyPolicyId, LimitKind, ModelId,
    NeverCancel, NoClaimState, ObservedIdentity, OrderedBytesStreamDisposition,
    OrderedBytesStreamError, OrderedBytesStreamPhase, Presence, ProblemSemanticId,
    SCHEMA_ID_HASH_DOMAIN, SchemaDescriptorKind, SchemaId, SemanticId, SourceByteId, SourceId,
    StrongIdentity, TrustState, VerifierId, WireContentId, WireType, adjudicate,
    legacy::LegacyProvenanceV1,
};
use fs_blake3::{ContentHash, hash_bytes, hash_domain};

const LIMITS: CanonicalLimits = CanonicalLimits::new(64 * 1024, 16 * 1024, 64, 1024, 7);

enum LeafV1 {}
impl CanonicalSchema for LeafV1 {
    const DOMAIN: &'static str = "org.frankensim.test.identity.leaf.v1";
    const NAME: &'static str = "identity-test-leaf";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "G0 leaf identity fixture";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required("value", WireType::U64)];
}

enum LeafV2 {}
impl CanonicalSchema for LeafV2 {
    const DOMAIN: &'static str = "org.frankensim.test.identity.leaf.v1";
    const NAME: &'static str = "identity-test-leaf";
    const VERSION: u32 = 2;
    const CONTEXT: &'static str = "G0 leaf identity fixture";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required("value", WireType::U64)];
}

enum LeafV1StructuralTwin {}
impl CanonicalSchema for LeafV1StructuralTwin {
    const DOMAIN: &'static str = LeafV1::DOMAIN;
    const NAME: &'static str = LeafV1::NAME;
    const VERSION: u32 = LeafV1::VERSION;
    const CONTEXT: &'static str = LeafV1::CONTEXT;
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required("value", WireType::U64)];
}

enum OtherDomain {}
impl CanonicalSchema for OtherDomain {
    const DOMAIN: &'static str = "org.frankensim.test.identity.other-domain.v1";
    const NAME: &'static str = "identity-test-leaf";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "G0 leaf identity fixture";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required("value", WireType::U64)];
}

enum OtherName {}
impl CanonicalSchema for OtherName {
    const DOMAIN: &'static str = LeafV1::DOMAIN;
    const NAME: &'static str = "identity-test-leaf-renamed";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = LeafV1::CONTEXT;
    const FIELDS: &'static [FieldSpec] = LeafV1::FIELDS;
}

enum OtherContext {}
impl CanonicalSchema for OtherContext {
    const DOMAIN: &'static str = LeafV1::DOMAIN;
    const NAME: &'static str = LeafV1::NAME;
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "different purpose";
    const FIELDS: &'static [FieldSpec] = LeafV1::FIELDS;
}

enum OtherFieldName {}
impl CanonicalSchema for OtherFieldName {
    const DOMAIN: &'static str = LeafV1::DOMAIN;
    const NAME: &'static str = LeafV1::NAME;
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = LeafV1::CONTEXT;
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required("other", WireType::U64)];
}

enum OtherWireType {}
impl CanonicalSchema for OtherWireType {
    const DOMAIN: &'static str = LeafV1::DOMAIN;
    const NAME: &'static str = LeafV1::NAME;
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = LeafV1::CONTEXT;
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required("value", WireType::I64)];
}

enum OtherPresence {}
impl CanonicalSchema for OtherPresence {
    const DOMAIN: &'static str = LeafV1::DOMAIN;
    const NAME: &'static str = LeafV1::NAME;
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = LeafV1::CONTEXT;
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::optional_bytes("value")];
}

enum RequiredBytesPresence {}
impl CanonicalSchema for RequiredBytesPresence {
    const DOMAIN: &'static str = LeafV1::DOMAIN;
    const NAME: &'static str = LeafV1::NAME;
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = LeafV1::CONTEXT;
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required("value", WireType::Bytes)];
}

enum HeaderOneField {}
impl CanonicalSchema for HeaderOneField {
    const DOMAIN: &'static str = "org.frankensim.test.identity.header-shape.v1";
    const NAME: &'static str = "identity-test-header-shape";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "G0 header field schema fixture";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required("a", WireType::U64)];
}

enum HeaderTwoFields {}
impl CanonicalSchema for HeaderTwoFields {
    const DOMAIN: &'static str = HeaderOneField::DOMAIN;
    const NAME: &'static str = HeaderOneField::NAME;
    const VERSION: u32 = HeaderOneField::VERSION;
    const CONTEXT: &'static str = HeaderOneField::CONTEXT;
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("a", WireType::U64),
        FieldSpec::required("b", WireType::I64),
    ];
}

enum HeaderTwoFieldsReordered {}
impl CanonicalSchema for HeaderTwoFieldsReordered {
    const DOMAIN: &'static str = HeaderOneField::DOMAIN;
    const NAME: &'static str = HeaderOneField::NAME;
    const VERSION: u32 = HeaderOneField::VERSION;
    const CONTEXT: &'static str = HeaderOneField::CONTEXT;
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("b", WireType::I64),
        FieldSpec::required("a", WireType::U64),
    ];
}

enum OptionalLeaf {}
impl CanonicalSchema for OptionalLeaf {
    const DOMAIN: &'static str = "org.frankensim.test.identity.optional.v1";
    const NAME: &'static str = "identity-test-optional";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "G0 optional fixture";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::optional_bytes("payload")];
}

enum FloatLeaf {}
impl CanonicalSchema for FloatLeaf {
    const DOMAIN: &'static str = "org.frankensim.test.identity.float.v1";
    const NAME: &'static str = "identity-test-float";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "G0 finite exact-bit float fixture";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required("value", WireType::FiniteF64)];
}

enum BytesLeaf {}
impl CanonicalSchema for BytesLeaf {
    const DOMAIN: &'static str = "org.frankensim.test.identity.bytes.v1";
    const NAME: &'static str = "identity-test-bytes";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "G0 streaming bytes fixture";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required("payload", WireType::Bytes)];
}

enum DescriptorDomainFloor {}
impl CanonicalSchema for DescriptorDomainFloor {
    const DOMAIN: &'static str =
        "org.frankensim.test.identity.schema-descriptor-domain-is-the-binding-floor.v1";
    const NAME: &'static str = "domain-floor";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "G0 domain floor";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required("value", WireType::U64)];
}

enum DescriptorNameFloor {}
impl CanonicalSchema for DescriptorNameFloor {
    const DOMAIN: &'static str = "org.frankensim.test.floor.name.v1";
    const NAME: &'static str = "schema-name-is-the-unique-longest-descriptor-floor";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "G0 name floor";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required("value", WireType::U64)];
}

enum DescriptorContextFloor {}
impl CanonicalSchema for DescriptorContextFloor {
    const DOMAIN: &'static str = "org.frankensim.test.floor.context.v1";
    const NAME: &'static str = "context-floor";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "G0 context is the unique longest schema descriptor floor";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required("value", WireType::U64)];
}

enum DescriptorFieldNameFloor {}
impl CanonicalSchema for DescriptorFieldNameFloor {
    const DOMAIN: &'static str = "org.frankensim.test.floor.field-name.v1";
    const NAME: &'static str = "field-name-floor";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "G0 field-name floor";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required(
        "field_name_is_the_unique_longest_schema_descriptor_floor",
        WireType::U64,
    )];
}

enum DescriptorNestedChildFloor {}
impl CanonicalSchema for DescriptorNestedChildFloor {
    const DOMAIN: &'static str = "org.frankensim.test.floor.nested-child.v1";
    const NAME: &'static str = "nested-child-floor";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str =
        "G0 nested child context is the unique longest recursive schema descriptor floor";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required("value", WireType::U64)];
}

static DESCRIPTOR_NESTED_CHILD_SPEC: ChildSpec =
    ChildSpec::for_identity::<SemanticId<DescriptorNestedChildFloor>>();

enum DescriptorNestedParent {}
impl CanonicalSchema for DescriptorNestedParent {
    const DOMAIN: &'static str = "org.frankensim.test.floor.nested-parent.v1";
    const NAME: &'static str = "nested-parent-floor";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "G0 nested parent";
    const FIELDS: &'static [FieldSpec] =
        &[FieldSpec::child_of("child", &DESCRIPTOR_NESTED_CHILD_SPEC)];
}

enum OtherSourceBytes {}
impl CanonicalSchema for OtherSourceBytes {
    const DOMAIN: &'static str = "org.frankensim.test.identity.other-source-bytes.v1";
    const NAME: &'static str = BytesLeaf::NAME;
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = BytesLeaf::CONTEXT;
    const FIELDS: &'static [FieldSpec] = BytesLeaf::FIELDS;
}

enum VariantLeaf {}
impl CanonicalSchema for VariantLeaf {
    const DOMAIN: &'static str = "org.frankensim.test.identity.variant.v1";
    const NAME: &'static str = "identity-test-variant";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "G0 variant fixture";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required("choice", WireType::Variant)];
}

enum SequenceLeaf {}
impl CanonicalSchema for SequenceLeaf {
    const DOMAIN: &'static str = "org.frankensim.test.identity.sequence.v1";
    const NAME: &'static str = "identity-test-sequence";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "G0 ordered sequence fixture";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required("items", WireType::OrderedBytes)];
}

enum TinySequenceLeaf {}
impl CanonicalSchema for TinySequenceLeaf {
    const DOMAIN: &'static str = "d";
    const NAME: &'static str = "n";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "c";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required("i", WireType::OrderedBytes)];
}

enum TinySetLeaf {}
impl CanonicalSchema for TinySetLeaf {
    const DOMAIN: &'static str = "d";
    const NAME: &'static str = "n";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "c";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required("i", WireType::CanonicalSet)];
}

enum SetLeaf {}
impl CanonicalSchema for SetLeaf {
    const DOMAIN: &'static str = "org.frankensim.test.identity.set.v1";
    const NAME: &'static str = "identity-test-set";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "G0 canonical set fixture";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required("items", WireType::CanonicalSet)];
}

enum ChildParent {}
impl CanonicalSchema for ChildParent {
    const DOMAIN: &'static str = "org.frankensim.test.identity.child-parent.v1";
    const NAME: &'static str = "identity-test-child-parent";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "G0 typed child fixture";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::child_of("child", &LEAF_SEMANTIC_CHILD)];
}

/// bead sj31i.52.10: the fixtures bind exactly a `SemanticId<LeafV1>`
/// child; every other role/schema now refuses instead of merely
/// hashing differently.
static LEAF_SEMANTIC_CHILD: ChildSpec = ChildSpec::for_identity::<SemanticId<LeafV1>>();

enum ChildrenParent {}
impl CanonicalSchema for ChildrenParent {
    const DOMAIN: &'static str = "org.frankensim.test.identity.children-parent.v1";
    const NAME: &'static str = "identity-test-children-parent";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "G0 ordered children fixture";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::ordered_children_of(
        "children",
        &LEAF_SEMANTIC_CHILD,
    )];
}

enum AllFields {}
impl CanonicalSchema for AllFields {
    const DOMAIN: &'static str = "org.frankensim.test.identity.all-fields.v1";
    const NAME: &'static str = "identity-test-all-fields";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "G0 complete field mutation fixture";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("text", WireType::Utf8),
        FieldSpec::required("bytes", WireType::Bytes),
        FieldSpec::required("unsigned", WireType::U64),
        FieldSpec::required("signed", WireType::I64),
        FieldSpec::required("flag", WireType::Bool),
        FieldSpec::required("float", WireType::FiniteF64),
        FieldSpec::optional_bytes("optional"),
        FieldSpec::required("variant", WireType::Variant),
        FieldSpec::required("sequence", WireType::OrderedBytes),
        FieldSpec::required("set", WireType::CanonicalSet),
        FieldSpec::child_of("child", &LEAF_SEMANTIC_CHILD),
        FieldSpec::ordered_children_of("children", &LEAF_SEMANTIC_CHILD),
    ];
}

enum VerifierSchema {}
impl CanonicalSchema for VerifierSchema {
    const DOMAIN: &'static str = "org.frankensim.test.identity.verifier.v1";
    const NAME: &'static str = "identity-test-verifier";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "G0 authority verifier fixture";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required("version", WireType::U64)];
}

enum PolicySchema {}
impl CanonicalSchema for PolicySchema {
    const DOMAIN: &'static str = "org.frankensim.test.identity.key-policy.v1";
    const NAME: &'static str = "identity-test-key-policy";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "G0 authority key-policy fixture";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required("version", WireType::U64)];
}

enum EmptyDomainSchema {}
impl CanonicalSchema for EmptyDomainSchema {
    const DOMAIN: &'static str = "";
    const NAME: &'static str = "identity-test-invalid-empty-domain";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "G0 invalid schema fixture";
    const FIELDS: &'static [FieldSpec] = &[];
}

enum ZeroVersionSchema {}
impl CanonicalSchema for ZeroVersionSchema {
    const DOMAIN: &'static str = "org.frankensim.test.identity.invalid-zero-version.v1";
    const NAME: &'static str = "identity-test-invalid-zero-version";
    const VERSION: u32 = 0;
    const CONTEXT: &'static str = "G0 invalid schema fixture";
    const FIELDS: &'static [FieldSpec] = &[];
}

enum DuplicateFieldSchema {}
impl CanonicalSchema for DuplicateFieldSchema {
    const DOMAIN: &'static str = "org.frankensim.test.identity.invalid-duplicate-field.v1";
    const NAME: &'static str = "identity-test-invalid-duplicate-field";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "G0 invalid schema fixture";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("same", WireType::U64),
        FieldSpec::required("same", WireType::I64),
    ];
}

enum NestedChildSchema {}
impl CanonicalSchema for NestedChildSchema {
    const DOMAIN: &'static str = "org.frankensim.test.identity.nested-child.v1";
    const NAME: &'static str = "identity-test-nested-child";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "G0 recursively valid child-schema fixture";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::child_of("leaf", &LEAF_SEMANTIC_CHILD),
        FieldSpec::required("tag", WireType::U64),
    ];
}

enum NestedChildSchemaDrift {}
impl CanonicalSchema for NestedChildSchemaDrift {
    const DOMAIN: &'static str = NestedChildSchema::DOMAIN;
    const NAME: &'static str = NestedChildSchema::NAME;
    const VERSION: u32 = NestedChildSchema::VERSION;
    const CONTEXT: &'static str = NestedChildSchema::CONTEXT;
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::child_of("leaf", &LEAF_SEMANTIC_CHILD),
        FieldSpec::required("tag", WireType::I64),
    ];
}

static NESTED_CHILD_SPEC: ChildSpec = ChildSpec::for_identity::<SemanticId<NestedChildSchema>>();

enum NestedParentSchema {}
impl CanonicalSchema for NestedParentSchema {
    const DOMAIN: &'static str = "org.frankensim.test.identity.nested-parent.v1";
    const NAME: &'static str = "identity-test-nested-parent";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "G0 recursive schema-admission fixture";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::child_of("child", &NESTED_CHILD_SPEC)];
}

static EMPTY_DOMAIN_CHILD_SPEC: ChildSpec =
    ChildSpec::for_identity::<SemanticId<EmptyDomainSchema>>();

enum InvalidNestedChildSchema {}
impl CanonicalSchema for InvalidNestedChildSchema {
    const DOMAIN: &'static str = "org.frankensim.test.identity.invalid-nested-child.v1";
    const NAME: &'static str = "identity-test-invalid-nested-child";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "G0 invalid recursive schema fixture";
    const FIELDS: &'static [FieldSpec] =
        &[FieldSpec::child_of("invalid", &EMPTY_DOMAIN_CHILD_SPEC)];
}

static INVALID_NESTED_CHILD_SPEC: ChildSpec =
    ChildSpec::for_identity::<SemanticId<InvalidNestedChildSchema>>();

enum InvalidNestedParentSchema {}
impl CanonicalSchema for InvalidNestedParentSchema {
    const DOMAIN: &'static str = "org.frankensim.test.identity.invalid-nested-parent.v1";
    const NAME: &'static str = "identity-test-invalid-nested-parent";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "G0 invalid recursive schema parent fixture";
    const FIELDS: &'static [FieldSpec] =
        &[FieldSpec::child_of("child", &INVALID_NESTED_CHILD_SPEC)];
}

enum Depth17LeafSchema {}
impl CanonicalSchema for Depth17LeafSchema {
    const DOMAIN: &'static str = "org.frankensim.test.identity.depth-17-leaf.v1";
    const NAME: &'static str = "identity-test-depth-17-leaf";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "G0 finite over-depth schema leaf";
    const FIELDS: &'static [FieldSpec] = &[];
}

static DEPTH_17_SPEC: ChildSpec = ChildSpec::for_identity::<SemanticId<Depth17LeafSchema>>();

macro_rules! depth_schema {
    ($schema:ident, $binding:ident, $child:ident) => {
        enum $schema {}
        impl CanonicalSchema for $schema {
            const DOMAIN: &'static str =
                concat!("org.frankensim.test.identity.", stringify!($schema), ".v1");
            const NAME: &'static str = stringify!($schema);
            const VERSION: u32 = 1;
            const CONTEXT: &'static str = "G0 finite over-depth schema chain";
            const FIELDS: &'static [FieldSpec] = &[FieldSpec::child_of("child", &$child)];
        }
        static $binding: ChildSpec = ChildSpec::for_identity::<SemanticId<$schema>>();
    };
}

depth_schema!(Depth16Schema, DEPTH_16_SPEC, DEPTH_17_SPEC);
depth_schema!(Depth15Schema, DEPTH_15_SPEC, DEPTH_16_SPEC);
depth_schema!(Depth14Schema, DEPTH_14_SPEC, DEPTH_15_SPEC);
depth_schema!(Depth13Schema, DEPTH_13_SPEC, DEPTH_14_SPEC);
depth_schema!(Depth12Schema, DEPTH_12_SPEC, DEPTH_13_SPEC);
depth_schema!(Depth11Schema, DEPTH_11_SPEC, DEPTH_12_SPEC);
depth_schema!(Depth10Schema, DEPTH_10_SPEC, DEPTH_11_SPEC);
depth_schema!(Depth09Schema, DEPTH_09_SPEC, DEPTH_10_SPEC);
depth_schema!(Depth08Schema, DEPTH_08_SPEC, DEPTH_09_SPEC);
depth_schema!(Depth07Schema, DEPTH_07_SPEC, DEPTH_08_SPEC);
depth_schema!(Depth06Schema, DEPTH_06_SPEC, DEPTH_07_SPEC);
depth_schema!(Depth05Schema, DEPTH_05_SPEC, DEPTH_06_SPEC);
depth_schema!(Depth04Schema, DEPTH_04_SPEC, DEPTH_05_SPEC);
depth_schema!(Depth03Schema, DEPTH_03_SPEC, DEPTH_04_SPEC);
depth_schema!(Depth02Schema, DEPTH_02_SPEC, DEPTH_03_SPEC);
depth_schema!(Depth01Schema, DEPTH_01_SPEC, DEPTH_02_SPEC);

enum OverDepthParentSchema {}
impl CanonicalSchema for OverDepthParentSchema {
    const DOMAIN: &'static str = "org.frankensim.test.identity.over-depth-parent.v1";
    const NAME: &'static str = "identity-test-over-depth-parent";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "G0 finite over-depth schema root";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::child_of("child", &DEPTH_01_SPEC)];
}

enum BranchLeafSchema {}
impl CanonicalSchema for BranchLeafSchema {
    const DOMAIN: &'static str = "org.frankensim.test.identity.branch-leaf.v1";
    const NAME: &'static str = "identity-test-branch-leaf";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "G0 recursive expansion-cap leaf";
    const FIELDS: &'static [FieldSpec] = &[];
}

static BRANCH_LEAF_SPEC: ChildSpec = ChildSpec::for_identity::<SemanticId<BranchLeafSchema>>();

macro_rules! branch_schema {
    ($schema:ident, $binding:ident, $child:ident) => {
        enum $schema {}
        impl CanonicalSchema for $schema {
            const DOMAIN: &'static str =
                concat!("org.frankensim.test.identity.", stringify!($schema), ".v1");
            const NAME: &'static str = stringify!($schema);
            const VERSION: u32 = 1;
            const CONTEXT: &'static str = "G0 recursive expansion-cap branch";
            const FIELDS: &'static [FieldSpec] = &[
                FieldSpec::child_of("left", &$child),
                FieldSpec::child_of("right", &$child),
            ];
        }
        static $binding: ChildSpec = ChildSpec::for_identity::<SemanticId<$schema>>();
    };
}

branch_schema!(Branch01Schema, BRANCH_01_SPEC, BRANCH_LEAF_SPEC);
branch_schema!(Branch02Schema, BRANCH_02_SPEC, BRANCH_01_SPEC);
branch_schema!(Branch03Schema, BRANCH_03_SPEC, BRANCH_02_SPEC);
branch_schema!(Branch04Schema, BRANCH_04_SPEC, BRANCH_03_SPEC);

enum Branch05Schema {}
impl CanonicalSchema for Branch05Schema {
    const DOMAIN: &'static str = "org.frankensim.test.identity.branch-05.v1";
    const NAME: &'static str = "identity-test-branch-05";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "G0 recursive expansion-cap root";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::child_of("left", &BRANCH_04_SPEC),
        FieldSpec::child_of("right", &BRANCH_04_SPEC),
    ];
}

fn limits_with_max_field_bytes(max_field_bytes: u64) -> CanonicalLimits {
    CanonicalLimits::new(
        LIMITS.max_canonical_bytes(),
        max_field_bytes,
        LIMITS.max_fields(),
        LIMITS.max_collection_items(),
        LIMITS.cancellation_poll_bytes(),
    )
}

fn top_level_descriptor_floor<D: CanonicalSchema>() -> (SchemaDescriptorKind, u64) {
    let longest_field_name = D::FIELDS
        .iter()
        .copied()
        .map(|field| field.name().len())
        .max()
        .unwrap_or(0);
    let candidates = [
        (SchemaDescriptorKind::Domain, D::DOMAIN.len()),
        (SchemaDescriptorKind::Name, D::NAME.len()),
        (SchemaDescriptorKind::Context, D::CONTEXT.len()),
        (SchemaDescriptorKind::FieldName, longest_field_name),
    ];
    let mut binding = SchemaDescriptorKind::Domain;
    let mut floor = 0usize;
    let mut binding_count = 0usize;
    for (candidate, bytes) in candidates {
        if bytes > floor {
            binding = candidate;
            floor = bytes;
            binding_count = 1;
        } else if bytes == floor {
            binding_count += 1;
        }
    }
    assert_eq!(
        binding_count, 1,
        "the fixture must have one unique binding descriptor term"
    );
    (
        binding,
        u64::try_from(floor).expect("static descriptor length fits u64"),
    )
}

fn descriptor_floor_receipt<D: CanonicalSchema>(
    max_field_bytes: u64,
) -> Result<fs_blake3::identity::IdentityReceipt<SemanticId<D>>, CanonicalError> {
    CanonicalEncoder::<SemanticId<D>, _>::new(
        limits_with_max_field_bytes(max_field_bytes),
        NeverCancel,
    )?
    .u64(Field::new(0, D::FIELDS[0].name()), 17)?
    .finish()
}

fn assert_descriptor_floor_boundaries<D: CanonicalSchema>(expected_binding: SchemaDescriptorKind) {
    let (binding, floor) = top_level_descriptor_floor::<D>();
    assert_eq!(binding, expected_binding);
    let supplied = floor.checked_sub(1).expect("descriptor floor is positive");
    let error = CanonicalEncoder::<SemanticId<D>, _>::new(
        limits_with_max_field_bytes(supplied),
        NeverCancel,
    )
    .expect_err("one-under the schema descriptor floor must refuse construction");
    assert_eq!(
        error,
        CanonicalError::MaxFieldBytesBelowSchemaFloor {
            binding,
            floor,
            supplied,
        }
    );
    assert!(
        !matches!(error, CanonicalError::LimitExceeded { .. }),
        "schema admission and payload-limit refusals must remain distinct"
    );
    assert_eq!(
        error.to_string(),
        format!(
            "canonical identity max_field_bytes {supplied} is below schema-descriptor \
             floor {floor}, bound by {binding:?}"
        )
    );

    let exact = descriptor_floor_receipt::<D>(floor).expect("exact descriptor floor admits");
    let one_over =
        descriptor_floor_receipt::<D>(floor + 1).expect("one over descriptor floor admits");
    let baseline = descriptor_floor_receipt::<D>(LIMITS.max_field_bytes())
        .expect("baseline descriptor limit admits");
    for admitted in [one_over, baseline] {
        assert_eq!(exact.id(), admitted.id());
        assert_eq!(exact.canonical_preimage(), admitted.canonical_preimage());
        assert_eq!(exact.canonical_bytes(), admitted.canonical_bytes());
    }
}

fn leaf(value: u64) -> fs_blake3::identity::IdentityReceipt<SemanticId<LeafV1>> {
    CanonicalEncoder::<SemanticId<LeafV1>, _>::new(LIMITS, NeverCancel)
        .expect("valid static leaf schema")
        .u64(Field::new(0, "value"), value)
        .expect("leaf field")
        .finish()
        .expect("leaf receipt")
}

fn leaf_v2(value: u64) -> fs_blake3::identity::IdentityReceipt<SemanticId<LeafV2>> {
    CanonicalEncoder::<SemanticId<LeafV2>, _>::new(LIMITS, NeverCancel)
        .expect("valid static leaf schema")
        .u64(Field::new(0, "value"), value)
        .expect("leaf field")
        .finish()
        .expect("leaf receipt")
}

fn other_domain(value: u64) -> fs_blake3::identity::IdentityReceipt<SemanticId<OtherDomain>> {
    CanonicalEncoder::<SemanticId<OtherDomain>, _>::new(LIMITS, NeverCancel)
        .expect("valid static leaf schema")
        .u64(Field::new(0, "value"), value)
        .expect("leaf field")
        .finish()
        .expect("leaf receipt")
}

#[allow(clippy::too_many_arguments)] // One explicit value per canonical wire form.
fn build_all(
    text: &str,
    bytes: &[u8],
    unsigned: u64,
    signed: i64,
    flag: bool,
    float: f64,
    optional: Option<&[u8]>,
    variant_tag: u32,
    variant_payload: &[u8],
    sequence: &[&[u8]],
    set: &[&[u8]],
    child_value: u64,
    children: [u64; 2],
) -> fs_blake3::identity::IdentityReceipt<SemanticId<AllFields>> {
    let child = leaf(child_value);
    let child_a = leaf(children[0]);
    let child_b = leaf(children[1]);
    CanonicalEncoder::<SemanticId<AllFields>, _>::new(LIMITS, NeverCancel)
        .expect("valid all-fields schema")
        .utf8(Field::new(0, "text"), text)
        .expect("text")
        .bytes(Field::new(1, "bytes"), bytes)
        .expect("bytes")
        .u64(Field::new(2, "unsigned"), unsigned)
        .expect("unsigned")
        .i64(Field::new(3, "signed"), signed)
        .expect("signed")
        .flag(Field::new(4, "flag"), flag)
        .expect("flag")
        .finite_f64(Field::new(5, "float"), float)
        .expect("finite float")
        .optional_bytes(Field::new(6, "optional"), optional)
        .expect("optional")
        .variant(Field::new(7, "variant"), variant_tag, variant_payload)
        .expect("variant")
        .ordered_bytes(
            Field::new(8, "sequence"),
            sequence.len() as u64,
            sequence.iter().copied(),
        )
        .expect("sequence")
        .canonical_set(Field::new(9, "set"), set.len() as u64, set.iter().copied())
        .expect("set")
        .child(Field::new(10, "child"), child.id())
        .expect("child")
        .ordered_children(Field::new(11, "children"), 2, [child_a.id(), child_b.id()])
        .expect("children")
        .finish()
        .expect("all-fields receipt")
}

fn all_baseline() -> fs_blake3::identity::IdentityReceipt<SemanticId<AllFields>> {
    build_all(
        "alpha",
        b"a|b\0c",
        7,
        -9,
        true,
        -0.0,
        Some(b"present"),
        3,
        b"variant",
        &[b"a", b"bc"],
        &[b"a", b"b"],
        11,
        [13, 17],
    )
}

#[test]
fn official_and_independent_blake3_vectors_remain_foundational() {
    assert_eq!(
        hash_bytes(b"").to_hex(),
        "af1349b9f5f9a1a6a0404dea36dcc9499bcb25c9adc112b7cc9a93cae41f3262"
    );
    assert_eq!(
        hash_bytes(b"abc").to_hex(),
        "6437b3ac38465133ffb63b75273a8db548c558465d79db03fd359c6cd5bd9d85"
    );
    let data: Vec<u8> = (0..1025).map(|index| (index % 251) as u8).collect();
    assert_eq!(
        hash_domain("BLAKE3 2019-12-27 16:29:52 test vectors context", &data).to_hex(),
        "effaa245f065fbf82ac186839a249707c3bddf6d3fdda22d1b95a3c970379bcb"
    );
}

#[test]
fn canonical_v1_tags_are_explicit_and_exhaustive() {
    assert_eq!(Presence::Required.tag(), 1);
    assert_eq!(Presence::Optional.tag(), 2);
    for (wire, tag) in [
        (WireType::Utf8, 1),
        (WireType::Bytes, 2),
        (WireType::U64, 3),
        (WireType::I64, 4),
        (WireType::Bool, 5),
        (WireType::FiniteF64, 6),
        (WireType::Variant, 7),
        (WireType::OrderedBytes, 8),
        (WireType::CanonicalSet, 9),
        (WireType::Child, 10),
        (WireType::OrderedChildren, 11),
    ] {
        assert_eq!(wire.tag(), tag);
    }
    for (role, tag) in [
        (IdentityRole::Semantic, 1),
        (IdentityRole::WireContent, 2),
        (IdentityRole::EvidenceNode, 3),
        (IdentityRole::Entity, 4),
        (IdentityRole::SourceBytes, 5),
        (IdentityRole::Source, 6),
        (IdentityRole::Model, 7),
        (IdentityRole::Checker, 8),
        (IdentityRole::Schema, 9),
        (IdentityRole::Verifier, 10),
        (IdentityRole::KeyPolicy, 11),
        (IdentityRole::ProblemSemantic, 12),
    ] {
        assert_eq!(role.tag(), tag);
    }
}

fn push_len_bytes(out: &mut Vec<u8>, bytes: &[u8]) {
    out.extend_from_slice(&(bytes.len() as u64).to_le_bytes());
    out.extend_from_slice(bytes);
}

const MANUAL_MAX_SCHEMA_CHILD_DEPTH: u32 = 16;

#[derive(Clone, Copy)]
struct ManualDescriptorSource<'a> {
    domain: &'a str,
    name: &'a str,
    version: u32,
    context: &'a str,
    fields: &'a [FieldSpec],
}

fn push_manual_descriptor_at(
    out: &mut Vec<u8>,
    source: ManualDescriptorSource<'_>,
    depth: u32,
    poison_offsets: &mut Vec<usize>,
) {
    out.extend_from_slice(b"FSSCHEM\x02");
    out.extend_from_slice(&CANONICAL_FRAME_VERSION.to_le_bytes());
    push_len_bytes(out, source.domain.as_bytes());
    push_len_bytes(out, source.name.as_bytes());
    out.extend_from_slice(&source.version.to_le_bytes());
    push_len_bytes(out, source.context.as_bytes());
    out.extend_from_slice(&(source.fields.len() as u64).to_le_bytes());
    for field in source.fields {
        push_len_bytes(out, field.name().as_bytes());
        out.extend_from_slice(&[field.wire_type().tag(), field.presence().tag()]);
        // bead sj31i.52.10: the v2 descriptor binds child fields
        // recursively (0 = unbound scalar, 1 + role + descriptor,
        // 2 = a child edge beyond the supported recursive depth).
        match field.child_spec() {
            None => out.push(0),
            Some(_) if depth >= MANUAL_MAX_SCHEMA_CHILD_DEPTH => {
                poison_offsets.push(out.len());
                out.push(2);
            }
            Some(child) => {
                out.extend_from_slice(&[1, child.role().tag()]);
                push_manual_descriptor_at(
                    out,
                    ManualDescriptorSource {
                        domain: child.domain(),
                        name: child.name(),
                        version: child.version(),
                        context: child.context(),
                        fields: child.fields(),
                    },
                    depth + 1,
                    poison_offsets,
                );
            }
        }
    }
}

fn manual_schema_descriptor_with_poison_offsets<D: CanonicalSchema>() -> (Vec<u8>, Vec<usize>) {
    let mut out = Vec::new();
    let mut poison_offsets = Vec::new();
    push_manual_descriptor_at(
        &mut out,
        ManualDescriptorSource {
            domain: D::DOMAIN,
            name: D::NAME,
            version: D::VERSION,
            context: D::CONTEXT,
            fields: D::FIELDS,
        },
        0,
        &mut poison_offsets,
    );
    (out, poison_offsets)
}

fn manual_schema_descriptor<D: CanonicalSchema>() -> Vec<u8> {
    manual_schema_descriptor_with_poison_offsets::<D>().0
}

fn manual_leaf_frame(role: IdentityRole, value: u64) -> Vec<u8> {
    let mut out = Vec::new();
    out.extend_from_slice(b"FSID\0\0\0\x01");
    out.extend_from_slice(&CANONICAL_FRAME_VERSION.to_le_bytes());
    out.extend_from_slice(&[role.tag(), 1]);
    push_len_bytes(&mut out, LeafV1::DOMAIN.as_bytes());
    push_len_bytes(&mut out, LeafV1::NAME.as_bytes());
    out.extend_from_slice(SchemaId::<LeafV1>::for_schema().as_bytes());
    out.extend_from_slice(&LeafV1::VERSION.to_le_bytes());
    push_len_bytes(&mut out, LeafV1::CONTEXT.as_bytes());
    out.extend_from_slice(&1u32.to_le_bytes());
    out.extend_from_slice(&0u32.to_le_bytes());
    push_len_bytes(&mut out, b"value");
    out.extend_from_slice(&[WireType::U64.tag(), Presence::Required.tag()]);
    out.push(0xf0);
    out.extend_from_slice(&0u32.to_le_bytes());
    push_len_bytes(&mut out, b"value");
    out.extend_from_slice(&[WireType::U64.tag(), Presence::Required.tag()]);
    out.extend_from_slice(&value.to_le_bytes());
    out.push(0xff);
    out.extend_from_slice(&1u32.to_le_bytes());
    out
}

fn manual_sequence_frame(rows: &[&[u8]]) -> Vec<u8> {
    let mut out = Vec::new();
    out.extend_from_slice(b"FSID\0\0\0\x01");
    out.extend_from_slice(&CANONICAL_FRAME_VERSION.to_le_bytes());
    out.extend_from_slice(&[IdentityRole::Semantic.tag(), 1]);
    push_len_bytes(&mut out, SequenceLeaf::DOMAIN.as_bytes());
    push_len_bytes(&mut out, SequenceLeaf::NAME.as_bytes());
    out.extend_from_slice(SchemaId::<SequenceLeaf>::for_schema().as_bytes());
    out.extend_from_slice(&SequenceLeaf::VERSION.to_le_bytes());
    push_len_bytes(&mut out, SequenceLeaf::CONTEXT.as_bytes());
    out.extend_from_slice(&1u32.to_le_bytes());
    out.extend_from_slice(&0u32.to_le_bytes());
    push_len_bytes(&mut out, b"items");
    out.extend_from_slice(&[WireType::OrderedBytes.tag(), Presence::Required.tag()]);
    out.push(0xf0);
    out.extend_from_slice(&0u32.to_le_bytes());
    push_len_bytes(&mut out, b"items");
    out.extend_from_slice(&[WireType::OrderedBytes.tag(), Presence::Required.tag()]);
    out.extend_from_slice(&(rows.len() as u64).to_le_bytes());
    for row in rows {
        push_len_bytes(&mut out, row);
    }
    out.push(0xff);
    out.extend_from_slice(&1u32.to_le_bytes());
    out
}

#[test]
fn manual_frame_parity_and_header_mutation_sensitivity() {
    let descriptor = manual_schema_descriptor::<LeafV1>();
    assert_eq!(
        *SchemaId::<LeafV1>::for_schema().as_bytes(),
        *hash_domain(SCHEMA_ID_HASH_DOMAIN, &descriptor).as_bytes()
    );

    let receipt = leaf(42);
    let frame = manual_leaf_frame(IdentityRole::Semantic, 42);
    assert_eq!(receipt.canonical_bytes(), frame.len() as u64);
    assert_eq!(receipt.canonical_preimage(), ContentId::of_bytes(&frame));
    assert_eq!(
        *receipt.id().as_bytes(),
        *hash_domain(CANONICAL_IDENTITY_HASH_DOMAIN, &frame).as_bytes()
    );

    let domain_len = 14usize;
    let domain = domain_len + 8;
    let schema_name_len = domain + LeafV1::DOMAIN.len();
    let schema_name = schema_name_len + 8;
    let schema_id = schema_name + LeafV1::NAME.len();
    let semantic_version = schema_id + 32;
    let context_len = semantic_version + 4;
    let context = context_len + 8;
    let declared_field_count = context + LeafV1::CONTEXT.len();
    let field_ordinal = declared_field_count + 4;
    let field_name_len = field_ordinal + 4;
    let field_name = field_name_len + 8;
    let field_wire = field_name + "value".len();
    let field_presence = field_wire + 1;
    let field_stream = field_presence + 1;
    let baseline = hash_domain(CANONICAL_IDENTITY_HASH_DOMAIN, &frame);
    for (index, field) in [
        (0usize, "canonical magic"),
        (8, "frame version"),
        (12, "role"),
        (13, "float policy"),
        (domain_len, "domain length"),
        (domain, "domain"),
        (schema_name_len, "schema-name length"),
        (schema_name, "schema name"),
        (schema_id, "schema ID"),
        (semantic_version, "semantic version"),
        (context_len, "context length"),
        (context, "context"),
        (declared_field_count, "declared field count"),
        (field_ordinal, "field ordinal"),
        (field_name_len, "field-name length"),
        (field_name, "field name"),
        (field_wire, "field wire type"),
        (field_presence, "field presence"),
        (field_stream, "canonical field stream"),
        (frame.len() - 5, "end marker"),
    ] {
        let mut moved = frame.clone();
        moved[index] ^= 1;
        assert_ne!(
            hash_domain(CANONICAL_IDENTITY_HASH_DOMAIN, &moved),
            baseline,
            "{field} byte {index} must move the root under the collision assumption"
        );
    }
}

#[test]
fn schema_descriptor_and_every_header_field_move_identity() {
    let baseline = SchemaId::<LeafV1>::for_schema().to_hex();
    let all_fields_descriptor = manual_schema_descriptor::<AllFields>();
    let all_fields_root = hash_domain(SCHEMA_ID_HASH_DOMAIN, &all_fields_descriptor);
    assert_eq!(
        *SchemaId::<AllFields>::for_schema().as_bytes(),
        *all_fields_root.as_bytes()
    );
    for index in 0..all_fields_descriptor.len() {
        let mut moved = all_fields_descriptor.clone();
        moved[index] ^= 1;
        assert_ne!(
            all_fields_root,
            hash_domain(SCHEMA_ID_HASH_DOMAIN, &moved),
            "schema descriptor byte {index} must move its root under the collision assumption"
        );
    }
    for moved in [
        SchemaId::<OtherDomain>::for_schema().to_hex(),
        SchemaId::<OtherName>::for_schema().to_hex(),
        SchemaId::<LeafV2>::for_schema().to_hex(),
        SchemaId::<OtherContext>::for_schema().to_hex(),
        SchemaId::<OtherFieldName>::for_schema().to_hex(),
        SchemaId::<OtherWireType>::for_schema().to_hex(),
        SchemaId::<OtherPresence>::for_schema().to_hex(),
        SchemaId::<HeaderTwoFields>::for_schema().to_hex(),
        SchemaId::<HeaderTwoFieldsReordered>::for_schema().to_hex(),
    ] {
        assert_ne!(baseline, moved);
    }
    assert_ne!(
        SchemaId::<HeaderOneField>::for_schema().to_hex(),
        SchemaId::<HeaderTwoFields>::for_schema().to_hex(),
        "declared field count is semantic"
    );
    assert_ne!(
        SchemaId::<HeaderTwoFields>::for_schema().to_hex(),
        SchemaId::<HeaderTwoFieldsReordered>::for_schema().to_hex(),
        "field order is semantic"
    );
}

#[test]
fn over_depth_schema_poison_tag_is_identity_bearing() {
    let (descriptor, poison_offsets) =
        manual_schema_descriptor_with_poison_offsets::<OverDepthParentSchema>();
    assert_eq!(
        poison_offsets.len(),
        1,
        "the finite seventeen-edge fixture has exactly one poisoned child edge"
    );
    let poison_offset = poison_offsets[0];
    assert_eq!(descriptor[poison_offset], 2);

    let baseline = hash_domain(SCHEMA_ID_HASH_DOMAIN, &descriptor);
    assert_eq!(
        SchemaId::<OverDepthParentSchema>::for_schema().as_bytes(),
        baseline.as_bytes(),
        "the independent depth-aware descriptor must reproduce SchemaId"
    );

    let mut moved = descriptor;
    moved[poison_offset] ^= 1;
    assert_ne!(
        baseline,
        hash_domain(SCHEMA_ID_HASH_DOMAIN, &moved),
        "the child-depth poison tag is schema-identity-bearing"
    );
}

#[test]
fn required_and_optional_byte_presence_have_distinct_schema_identities() {
    assert_ne!(
        SchemaId::<RequiredBytesPresence>::for_schema().to_hex(),
        SchemaId::<OtherPresence>::for_schema().to_hex()
    );
}

#[test]
fn schema_versions_are_nominal_and_digest_distinct() {
    let first = leaf(5);
    let second = leaf_v2(5);
    assert_ne!(first.id().to_hex(), second.id().to_hex());
    assert_ne!(
        first.canonical_preimage().to_hex(),
        second.canonical_preimage().to_hex()
    );
    assert_ne!(
        SchemaId::<LeafV1>::for_schema().to_hex(),
        SchemaId::<LeafV2>::for_schema().to_hex()
    );
}

#[test]
#[allow(clippy::too_many_lines)] // One exhaustive role-to-nominal-type matrix.
fn roles_domains_versions_and_raw_content_are_separate() {
    assert_eq!(SemanticId::<LeafV1>::ROLE, IdentityRole::Semantic);
    assert_eq!(WireContentId::<LeafV1>::ROLE, IdentityRole::WireContent);
    assert_eq!(EvidenceNodeId::<LeafV1>::ROLE, IdentityRole::EvidenceNode);
    assert_eq!(EntityId::<LeafV1>::ROLE, IdentityRole::Entity);
    assert_eq!(SourceByteId::<LeafV1>::ROLE, IdentityRole::SourceBytes);
    assert_eq!(SourceId::<LeafV1>::ROLE, IdentityRole::Source);
    assert_eq!(ModelId::<LeafV1>::ROLE, IdentityRole::Model);
    assert_eq!(CheckerId::<LeafV1>::ROLE, IdentityRole::Checker);
    assert_eq!(SchemaId::<LeafV1>::ROLE, IdentityRole::Schema);
    assert_eq!(VerifierId::<LeafV1>::ROLE, IdentityRole::Verifier);
    assert_eq!(KeyPolicyId::<LeafV1>::ROLE, IdentityRole::KeyPolicy);
    assert_eq!(
        ProblemSemanticId::<LeafV1>::ROLE,
        IdentityRole::ProblemSemantic
    );

    let semantic = leaf(8);
    let source = CanonicalEncoder::<SourceId<LeafV1>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .u64(Field::new(0, "value"), 8)
        .unwrap()
        .finish()
        .unwrap();
    let model = CanonicalEncoder::<ModelId<LeafV1>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .u64(Field::new(0, "value"), 8)
        .unwrap()
        .finish()
        .unwrap();
    let checker = CanonicalEncoder::<CheckerId<LeafV1>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .u64(Field::new(0, "value"), 8)
        .unwrap()
        .finish()
        .unwrap();
    let evidence = CanonicalEncoder::<EvidenceNodeId<LeafV1>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .u64(Field::new(0, "value"), 8)
        .unwrap()
        .finish()
        .unwrap();
    let entity = CanonicalEncoder::<EntityId<LeafV1>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .u64(Field::new(0, "value"), 8)
        .unwrap()
        .finish()
        .unwrap();
    let wire = CanonicalEncoder::<WireContentId<LeafV1>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .u64(Field::new(0, "value"), 8)
        .unwrap()
        .finish()
        .unwrap();
    let problem = CanonicalEncoder::<ProblemSemanticId<LeafV1>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .u64(Field::new(0, "value"), 8)
        .unwrap()
        .finish()
        .unwrap();
    let source_bytes = CanonicalEncoder::<SourceByteId<LeafV1>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .u64(Field::new(0, "value"), 8)
        .unwrap()
        .finish()
        .unwrap();
    let verifier = CanonicalEncoder::<VerifierId<LeafV1>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .u64(Field::new(0, "value"), 8)
        .unwrap()
        .finish()
        .unwrap();
    let key_policy = CanonicalEncoder::<KeyPolicyId<LeafV1>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .u64(Field::new(0, "value"), 8)
        .unwrap()
        .finish()
        .unwrap();
    let mut roots = vec![
        semantic.id().to_hex(),
        source.id().to_hex(),
        model.id().to_hex(),
        checker.id().to_hex(),
        evidence.id().to_hex(),
        entity.id().to_hex(),
        wire.id().to_hex(),
        problem.id().to_hex(),
        source_bytes.id().to_hex(),
        verifier.id().to_hex(),
        key_policy.id().to_hex(),
        SchemaId::<LeafV1>::for_schema().to_hex(),
    ];
    roots.sort();
    roots.dedup();
    assert_eq!(roots.len(), 12, "role/domain is part of every typed root");

    assert_ne!(semantic.id().to_hex(), other_domain(8).id().to_hex());
    let raw_frame = ContentId::of_bytes(&manual_leaf_frame(IdentityRole::Semantic, 8));
    assert_eq!(raw_frame, semantic.canonical_preimage());
    assert_ne!(raw_frame.to_hex(), semantic.id().to_hex());
}

#[test]
fn same_source_bytes_across_domains_keep_content_but_move_semantic_id() {
    let source_bytes = b"same exact input";
    let raw = ContentId::of_bytes(source_bytes);
    let first = CanonicalEncoder::<SourceByteId<BytesLeaf>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .bytes(Field::new(0, "payload"), source_bytes)
        .unwrap()
        .finish()
        .unwrap();

    let second = CanonicalEncoder::<SourceByteId<OtherSourceBytes>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .bytes(Field::new(0, "payload"), source_bytes)
        .unwrap()
        .finish()
        .unwrap();
    assert_eq!(raw, ContentId::of_bytes(source_bytes));
    assert_ne!(first.id().to_hex(), second.id().to_hex());
    assert_ne!(
        first.canonical_preimage().to_hex(),
        second.canonical_preimage().to_hex()
    );
}

#[test]
#[allow(clippy::too_many_lines)] // One ordered schema/header refusal matrix.
fn field_schema_order_presence_and_variants_are_non_confusable() {
    let none = CanonicalEncoder::<SemanticId<OptionalLeaf>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .optional_bytes(Field::new(0, "payload"), None)
        .unwrap()
        .finish()
        .unwrap();
    let empty = CanonicalEncoder::<SemanticId<OptionalLeaf>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .optional_bytes(Field::new(0, "payload"), Some(b""))
        .unwrap()
        .finish()
        .unwrap();
    assert_ne!(none.id(), empty.id());

    let first = CanonicalEncoder::<SemanticId<VariantLeaf>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .variant(Field::new(0, "choice"), 1, b"same")
        .unwrap()
        .finish()
        .unwrap();
    let second = CanonicalEncoder::<SemanticId<VariantLeaf>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .variant(Field::new(0, "choice"), 2, b"same")
        .unwrap()
        .finish()
        .unwrap();
    assert_ne!(first.id(), second.id());

    let wrong_ordinal = CanonicalEncoder::<SemanticId<LeafV1>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .u64(Field::new(1, "value"), 1);
    assert!(matches!(
        wrong_ordinal,
        Err(CanonicalError::FieldOrder {
            expected: 0,
            actual: 1
        })
    ));
    let wrong_name = CanonicalEncoder::<SemanticId<LeafV1>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .u64(Field::new(0, "other"), 1);
    assert!(matches!(wrong_name, Err(CanonicalError::FieldName)));
    let wrong_type = CanonicalEncoder::<SemanticId<LeafV1>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .i64(Field::new(0, "value"), 1);
    assert!(matches!(wrong_type, Err(CanonicalError::WireType)));
    let wrong_presence = CanonicalEncoder::<SemanticId<OptionalLeaf>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .bytes(Field::new(0, "payload"), b"present");
    assert!(matches!(wrong_presence, Err(CanonicalError::Presence)));
    let missing = CanonicalEncoder::<SemanticId<LeafV1>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .finish();
    assert!(matches!(
        missing,
        Err(CanonicalError::MissingFields {
            expected: 1,
            actual: 0
        })
    ));

    let one = CanonicalEncoder::<SemanticId<HeaderOneField>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .u64(Field::new(0, "a"), 7)
        .unwrap()
        .finish()
        .unwrap();
    let two = CanonicalEncoder::<SemanticId<HeaderTwoFields>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .u64(Field::new(0, "a"), 7)
        .unwrap()
        .i64(Field::new(1, "b"), -9)
        .unwrap()
        .finish()
        .unwrap();
    let reordered =
        CanonicalEncoder::<SemanticId<HeaderTwoFieldsReordered>, _>::new(LIMITS, NeverCancel)
            .unwrap()
            .i64(Field::new(0, "b"), -9)
            .unwrap()
            .u64(Field::new(1, "a"), 7)
            .unwrap()
            .finish()
            .unwrap();
    assert_ne!(one.id().to_hex(), two.id().to_hex());
    assert_ne!(two.id().to_hex(), reordered.id().to_hex());

    let frame = manual_leaf_frame(IdentityRole::Semantic, 1);
    let declared_count_offset = 8
        + 4
        + 2
        + 8
        + LeafV1::DOMAIN.len()
        + 8
        + LeafV1::NAME.len()
        + 32
        + 4
        + 8
        + LeafV1::CONTEXT.len();
    let field_ordinal_offset = declared_count_offset + 4;
    let field_wire_offset = field_ordinal_offset + 4 + 8 + "value".len();
    for index in [
        declared_count_offset,
        field_ordinal_offset,
        field_wire_offset,
    ] {
        let mut moved = frame.clone();
        moved[index] ^= 1;
        assert_ne!(
            hash_domain(CANONICAL_IDENTITY_HASH_DOMAIN, &frame),
            hash_domain(CANONICAL_IDENTITY_HASH_DOMAIN, &moved)
        );
    }
}

#[test]
fn finite_float_preserves_signed_zero_and_refuses_nonfinite() {
    let encode = |value| {
        CanonicalEncoder::<SemanticId<FloatLeaf>, _>::new(LIMITS, NeverCancel)
            .unwrap()
            .finite_f64(Field::new(0, "value"), value)
            .and_then(CanonicalEncoder::finish)
    };
    let positive = encode(0.0).unwrap();
    let negative = encode(-0.0).unwrap();
    assert_ne!(positive.id(), negative.id());
    for finite in [f64::MIN, f64::MAX, f64::MIN_POSITIVE, f64::from_bits(1)] {
        assert!(encode(finite).is_ok());
    }
    for bits in [
        f64::INFINITY.to_bits(),
        f64::NEG_INFINITY.to_bits(),
        f64::NAN.to_bits(),
        0x7ff0_0000_0000_0001,
        0xfff8_0000_0000_0001,
    ] {
        assert_eq!(
            encode(f64::from_bits(bits)),
            Err(CanonicalError::NonFiniteFloat { bits })
        );
    }
}

#[test]
#[allow(clippy::too_many_lines)] // Every independent semantic field moves exactly once.
fn every_semantic_field_is_mutation_sensitive() {
    let baseline = all_baseline().id().to_hex();
    let mutations = [
        build_all(
            "beta",
            b"a|b\0c",
            7,
            -9,
            true,
            -0.0,
            Some(b"present"),
            3,
            b"variant",
            &[b"a", b"bc"],
            &[b"a", b"b"],
            11,
            [13, 17],
        )
        .id()
        .to_hex(),
        build_all(
            "alpha",
            b"different",
            7,
            -9,
            true,
            -0.0,
            Some(b"present"),
            3,
            b"variant",
            &[b"a", b"bc"],
            &[b"a", b"b"],
            11,
            [13, 17],
        )
        .id()
        .to_hex(),
        build_all(
            "alpha",
            b"a|b\0c",
            8,
            -9,
            true,
            -0.0,
            Some(b"present"),
            3,
            b"variant",
            &[b"a", b"bc"],
            &[b"a", b"b"],
            11,
            [13, 17],
        )
        .id()
        .to_hex(),
        build_all(
            "alpha",
            b"a|b\0c",
            7,
            -8,
            true,
            -0.0,
            Some(b"present"),
            3,
            b"variant",
            &[b"a", b"bc"],
            &[b"a", b"b"],
            11,
            [13, 17],
        )
        .id()
        .to_hex(),
        build_all(
            "alpha",
            b"a|b\0c",
            7,
            -9,
            false,
            -0.0,
            Some(b"present"),
            3,
            b"variant",
            &[b"a", b"bc"],
            &[b"a", b"b"],
            11,
            [13, 17],
        )
        .id()
        .to_hex(),
        build_all(
            "alpha",
            b"a|b\0c",
            7,
            -9,
            true,
            0.0,
            Some(b"present"),
            3,
            b"variant",
            &[b"a", b"bc"],
            &[b"a", b"b"],
            11,
            [13, 17],
        )
        .id()
        .to_hex(),
        build_all(
            "alpha",
            b"a|b\0c",
            7,
            -9,
            true,
            -0.0,
            None,
            3,
            b"variant",
            &[b"a", b"bc"],
            &[b"a", b"b"],
            11,
            [13, 17],
        )
        .id()
        .to_hex(),
        build_all(
            "alpha",
            b"a|b\0c",
            7,
            -9,
            true,
            -0.0,
            Some(b"present"),
            4,
            b"variant",
            &[b"a", b"bc"],
            &[b"a", b"b"],
            11,
            [13, 17],
        )
        .id()
        .to_hex(),
        build_all(
            "alpha",
            b"a|b\0c",
            7,
            -9,
            true,
            -0.0,
            Some(b"present"),
            3,
            b"moved",
            &[b"a", b"bc"],
            &[b"a", b"b"],
            11,
            [13, 17],
        )
        .id()
        .to_hex(),
        build_all(
            "alpha",
            b"a|b\0c",
            7,
            -9,
            true,
            -0.0,
            Some(b"present"),
            3,
            b"variant",
            &[b"bc", b"a"],
            &[b"a", b"b"],
            11,
            [13, 17],
        )
        .id()
        .to_hex(),
        build_all(
            "alpha",
            b"a|b\0c",
            7,
            -9,
            true,
            -0.0,
            Some(b"present"),
            3,
            b"variant",
            &[b"a", b"bc"],
            &[b"a", b"c"],
            11,
            [13, 17],
        )
        .id()
        .to_hex(),
        build_all(
            "alpha",
            b"a|b\0c",
            7,
            -9,
            true,
            -0.0,
            Some(b"present"),
            3,
            b"variant",
            &[b"a", b"bc"],
            &[b"a", b"b"],
            12,
            [13, 17],
        )
        .id()
        .to_hex(),
        build_all(
            "alpha",
            b"a|b\0c",
            7,
            -9,
            true,
            -0.0,
            Some(b"present"),
            3,
            b"variant",
            &[b"a", b"bc"],
            &[b"a", b"b"],
            11,
            [17, 13],
        )
        .id()
        .to_hex(),
    ];
    assert!(mutations.iter().all(|moved| moved != &baseline));
}

#[test]
fn canonical_collections_bind_order_and_refuse_ambiguity() {
    let ab = CanonicalEncoder::<SemanticId<SequenceLeaf>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .ordered_bytes(Field::new(0, "items"), 2, [b"a".as_slice(), b"bc"])
        .unwrap()
        .finish()
        .unwrap();
    let ba = CanonicalEncoder::<SemanticId<SequenceLeaf>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .ordered_bytes(Field::new(0, "items"), 2, [b"bc".as_slice(), b"a"])
        .unwrap()
        .finish()
        .unwrap();
    assert_ne!(ab.id(), ba.id());

    let duplicate = CanonicalEncoder::<SemanticId<SetLeaf>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .canonical_set(Field::new(0, "items"), 2, [b"a".as_slice(), b"a"]);
    assert!(matches!(
        duplicate,
        Err(CanonicalError::DuplicateSetItem { index: 1 })
    ));
    let unsorted = CanonicalEncoder::<SemanticId<SetLeaf>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .canonical_set(Field::new(0, "items"), 2, [b"b".as_slice(), b"a"]);
    assert!(matches!(
        unsorted,
        Err(CanonicalError::NonCanonicalSetOrder { index: 1 })
    ));
}

#[test]
fn empty_canonical_set_charges_count_prefix_to_field_budget_before_iteration() {
    let below_count_prefix = CanonicalLimits::new(1024, 7, 4, 4, 4);
    let refused =
        CanonicalEncoder::<SemanticId<TinySetLeaf>, _>::new(below_count_prefix, NeverCancel)
            .unwrap()
            .canonical_set(Field::new(0, "i"), 0, PanicChunks);
    assert!(matches!(
        refused,
        Err(CanonicalError::LimitExceeded {
            kind: LimitKind::FieldBytes,
            requested: 8,
            limit: 7,
        })
    ));

    let exact_count_prefix = CanonicalLimits::new(1024, 8, 4, 4, 4);
    let admitted =
        CanonicalEncoder::<SemanticId<TinySetLeaf>, _>::new(exact_count_prefix, NeverCancel)
            .unwrap()
            .canonical_set(Field::new(0, "i"), 0, std::iter::empty::<&'static [u8]>())
            .unwrap()
            .finish()
            .unwrap();
    assert_eq!(admitted.collection_items(), 0);
    assert_eq!(admitted.limits(), exact_count_prefix);
}

#[test]
fn typed_children_bind_role_schema_order_and_full_digest() {
    let a = leaf(1);
    let b = leaf(2);
    let ab = CanonicalEncoder::<EvidenceNodeId<ChildrenParent>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .ordered_children(Field::new(0, "children"), 2, [a.id(), b.id()])
        .unwrap()
        .finish()
        .unwrap();
    let ba = CanonicalEncoder::<EvidenceNodeId<ChildrenParent>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .ordered_children(Field::new(0, "children"), 2, [b.id(), a.id()])
        .unwrap()
        .finish()
        .unwrap();
    assert_ne!(ab.id(), ba.id());

    let shared_digest = *a.id().as_bytes();
    let semantic_id = SemanticId::<LeafV1>::parse_slice(&shared_digest).unwrap();
    let source_id = SourceId::<LeafV1>::parse_slice(&shared_digest).unwrap();
    let other_schema_id = SemanticId::<OtherDomain>::parse_slice(&shared_digest).unwrap();
    // The declared binding admits exactly SemanticId<LeafV1>...
    CanonicalEncoder::<SemanticId<ChildParent>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .child(Field::new(0, "child"), semantic_id)
        .unwrap()
        .finish()
        .unwrap();
    // ...and REFUSES every other role/schema (bead sj31i.52.10): a
    // wrong child no longer merely hashes differently — it cannot be
    // constructed under this parent at all.
    let wrong_role = CanonicalEncoder::<SemanticId<ChildParent>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .child(Field::new(0, "child"), source_id)
        .expect_err("a Source child under a Semantic binding refuses");
    assert!(matches!(
        wrong_role,
        CanonicalError::ChildBindingMismatch {
            field: "child",
            what: "child role",
        }
    ));
    let wrong_schema = CanonicalEncoder::<SemanticId<ChildParent>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .child(Field::new(0, "child"), other_schema_id)
        .expect_err("an OtherDomain child under a LeafV1 binding refuses");
    assert!(matches!(
        wrong_schema,
        CanonicalError::ChildBindingMismatch { field: "child", .. }
    ));
}

enum ChildParentAltVersion {}
impl CanonicalSchema for ChildParentAltVersion {
    const DOMAIN: &'static str = "org.frankensim.test.identity.child-parent.v1";
    const NAME: &'static str = "identity-test-child-parent";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "G0 typed child fixture";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::child_of("child", &LEAF_V2_SEMANTIC_CHILD)];
}

enum LeafV2Marker {}
impl CanonicalSchema for LeafV2Marker {
    const DOMAIN: &'static str = LeafV1::DOMAIN;
    const NAME: &'static str = LeafV1::NAME;
    const VERSION: u32 = 2;
    const CONTEXT: &'static str = LeafV1::CONTEXT;
    const FIELDS: &'static [FieldSpec] = LeafV1::FIELDS;
}

static LEAF_V2_SEMANTIC_CHILD: ChildSpec = ChildSpec::for_identity::<SemanticId<LeafV2Marker>>();

/// bead sj31i.52.10: the expected child binding is part of the parent
/// schema identity — two parents identical except for the expected
/// child VERSION have different schema ids — and empty ordered-children
/// collections still validate their binding.
#[test]
fn child_bindings_are_schema_id_bearing_and_checked_when_empty() {
    let bound_v1 = *SchemaId::<ChildParent>::for_schema().as_bytes();
    let bound_v2 = *SchemaId::<ChildParentAltVersion>::for_schema().as_bytes();
    assert_ne!(
        bound_v1, bound_v2,
        "changing the expected child descriptor version must change the parent schema id"
    );

    // An EMPTY ordered-children collection still validates the binding:
    // the right type admits with zero items, the wrong type refuses
    // before any item is read.
    CanonicalEncoder::<EvidenceNodeId<ChildrenParent>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .ordered_children::<SemanticId<LeafV1>, _>(Field::new(0, "children"), 0, [])
        .unwrap()
        .finish()
        .unwrap();
    let refusal = CanonicalEncoder::<EvidenceNodeId<ChildrenParent>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .ordered_children::<SemanticId<OtherDomain>, _>(Field::new(0, "children"), 0, [])
        .expect_err("an empty collection of the wrong child schema still refuses");
    assert!(matches!(
        refusal,
        CanonicalError::ChildBindingMismatch {
            field: "children",
            ..
        }
    ));
}

#[test]
fn structural_child_bindings_admit_equivalent_markers_and_refuse_nested_drift() {
    assert_eq!(
        SchemaId::<LeafV1>::for_schema().as_bytes(),
        SchemaId::<LeafV1StructuralTwin>::for_schema().as_bytes(),
        "distinct markers with the same complete descriptor name the same schema"
    );
    let equivalent = SemanticId::<LeafV1StructuralTwin>::parse_slice(&[0x31; 32])
        .expect("typed equivalent-marker digest");
    CanonicalEncoder::<SemanticId<ChildParent>, _>::new(LIMITS, NeverCancel)
        .expect("valid parent schema")
        .child(Field::new(0, "child"), equivalent)
        .expect("role plus structural descriptor is the admission relation")
        .finish()
        .expect("equivalent-marker parent receipt");

    let nested_drift = SemanticId::<NestedChildSchemaDrift>::parse_slice(&[0x32; 32])
        .expect("typed nested-drift digest");
    let refusal = CanonicalEncoder::<SemanticId<NestedParentSchema>, _>::new(LIMITS, NeverCancel)
        .expect("valid recursively bound parent schema")
        .child(Field::new(0, "child"), nested_drift)
        .expect_err("a nested wire-type change must not satisfy the child binding");
    assert!(matches!(
        refusal,
        CanonicalError::ChildBindingMismatch {
            field: "child",
            what: "child field schema",
        }
    ));
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum RowFixtureError {
    LengthSource(&'static str),
    Producer(&'static str),
}

fn streamed_sequence(
    rows: &[&[u8]],
    limits: CanonicalLimits,
    schedule: u8,
) -> fs_blake3::identity::IdentityReceipt<SemanticId<SequenceLeaf>> {
    let lengths = rows
        .iter()
        .map(|row| Ok::<u64, RowFixtureError>(row.len() as u64));
    CanonicalEncoder::<SemanticId<SequenceLeaf>, _>::new(limits, NeverCancel)
        .unwrap()
        .ordered_bytes_stream(
            Field::new(0, "items"),
            rows.len() as u64,
            lengths,
            |row_index, mut sink| {
                let row = rows[row_index as usize];
                match schedule {
                    0 => sink
                        .write(row)
                        .expect("whole-row fixture must fit its declaration"),
                    1 => {
                        for byte in row {
                            sink.write(core::slice::from_ref(byte))
                                .expect("byte chunks must fit their declaration");
                        }
                    }
                    2 => {
                        sink.write(b"")
                            .expect("empty chunks are legal and non-semantic");
                        for chunk in row.chunks(3) {
                            sink.write(chunk)
                                .expect("uneven chunks must fit their declaration");
                        }
                        sink.write(b"")
                            .expect("trailing empty chunks are non-semantic");
                    }
                    _ => panic!("unknown fixture schedule"),
                }
                Ok(())
            },
        )
        .unwrap()
        .finish()
        .unwrap()
}

#[test]
fn stream_partition_and_non_cancelling_probes_are_invariant() {
    let data: Vec<u8> = (0..4097).map(|index| (index % 251) as u8).collect();
    let one = CanonicalEncoder::<SemanticId<BytesLeaf>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .bytes(Field::new(0, "payload"), &data)
        .unwrap()
        .finish()
        .unwrap();
    let chunks: Vec<&[u8]> = data.chunks(13).collect();
    let split = CanonicalEncoder::<SemanticId<BytesLeaf>, _>::new(LIMITS, || false)
        .unwrap()
        .bytes_stream(Field::new(0, "payload"), data.len() as u64, chunks)
        .unwrap()
        .finish()
        .unwrap();
    assert_eq!(one.id(), split.id());
    assert_eq!(one.canonical_preimage(), split.canonical_preimage());
    assert_eq!(one.canonical_bytes(), split.canonical_bytes());
}

#[test]
fn ordered_row_stream_matches_eager_and_independent_frame() {
    let rows: [&[u8]; 5] = [b"", b"a", b"b\0c", &[0xff, 0x80], b"longer-row"];
    let eager = CanonicalEncoder::<SemanticId<SequenceLeaf>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .ordered_bytes(Field::new(0, "items"), rows.len() as u64, rows)
        .unwrap()
        .finish()
        .unwrap();
    let streamed = streamed_sequence(&rows, LIMITS, 2);
    assert_eq!(streamed, eager);
    assert_eq!(streamed.collection_items(), rows.len() as u64);

    let frame = manual_sequence_frame(&rows);
    assert_eq!(streamed.canonical_bytes(), frame.len() as u64);
    assert_eq!(streamed.canonical_preimage(), ContentId::of_bytes(&frame));
    assert_eq!(
        streamed.id().as_bytes(),
        hash_domain(CANONICAL_IDENTITY_HASH_DOMAIN, &frame).as_bytes()
    );

    let no_rows: [&[u8]; 0] = [];
    let empty = streamed_sequence(&no_rows, LIMITS, 0);
    let empty_frame = manual_sequence_frame(&no_rows);
    assert_eq!(
        empty.canonical_preimage(),
        ContentId::of_bytes(&empty_frame)
    );
    assert_eq!(empty.collection_items(), 0);

    let one_empty: [&[u8]; 1] = [b""];
    let singleton = streamed_sequence(&one_empty, LIMITS, 2);
    assert_ne!(empty.id(), singleton.id());
    assert_eq!(singleton.collection_items(), 1);

    let cancel_flag = Cell::new(false);
    let borrowed_probe =
        CanonicalEncoder::<SemanticId<SequenceLeaf>, _>::new(LIMITS, || cancel_flag.get())
            .unwrap()
            .ordered_bytes_stream(
                Field::new(0, "items"),
                rows.len() as u64,
                rows.iter()
                    .map(|row| Ok::<u64, RowFixtureError>(row.len() as u64)),
                |row_index, mut sink| {
                    sink.write(rows[row_index as usize])
                        .map_err(|_| RowFixtureError::Producer("mapped-canonical-source"))?;
                    Ok(())
                },
            )
            .unwrap()
            .finish()
            .unwrap();
    assert_eq!(borrowed_probe, eager);
}

#[test]
fn ordered_row_stream_chunk_partition_and_schedule_are_nonsemantic() {
    let rows: [&[u8]; 3] = [b"abcdef", b"", b"ghijklmnop"];
    let whole = streamed_sequence(&rows, LIMITS, 0);
    let bytes = streamed_sequence(&rows, LIMITS, 1);
    let uneven = streamed_sequence(&rows, LIMITS, 2);
    assert_eq!(whole.id(), bytes.id());
    assert_eq!(whole.id(), uneven.id());
    assert_eq!(whole.canonical_preimage(), bytes.canonical_preimage());
    assert_eq!(whole.canonical_preimage(), uneven.canonical_preimage());
    assert_eq!(whole.canonical_bytes(), uneven.canonical_bytes());
    assert_eq!(whole.collection_items(), 3);
    assert_eq!(uneven.collection_items(), 3);

    for stride in [1, 7, 4096] {
        let limits = CanonicalLimits::new(64 * 1024, 16 * 1024, 64, 1024, stride);
        let scheduled = streamed_sequence(&rows, limits, 2);
        assert_eq!(scheduled.id(), whole.id());
        assert_eq!(scheduled.canonical_preimage(), whole.canonical_preimage());
        assert_eq!(scheduled.canonical_bytes(), whole.canonical_bytes());
    }

    let moved_boundaries: [&[u8]; 2] = [b"abc", b"def"];
    let one_row: [&[u8]; 1] = [b"abcdef"];
    let reordered: [&[u8]; 3] = [b"ghijklmnop", b"", b"abcdef"];
    assert_ne!(
        streamed_sequence(&moved_boundaries, LIMITS, 0).id(),
        streamed_sequence(&one_row, LIMITS, 0).id()
    );
    assert_ne!(whole.id(), streamed_sequence(&reordered, LIMITS, 0).id());
}

struct PanicChunks;
impl IntoIterator for PanicChunks {
    type Item = &'static [u8];
    type IntoIter = std::iter::Empty<&'static [u8]>;

    fn into_iter(self) -> Self::IntoIter {
        panic!("hostile declared length/count must refuse before reading the iterator")
    }
}

struct PanicLengths;
impl IntoIterator for PanicLengths {
    type Item = Result<u64, RowFixtureError>;
    type IntoIter = std::iter::Empty<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        panic!("hostile declared count must refuse before reading row lengths")
    }
}

struct ArmBeforeSecond<'a> {
    first: Option<&'a [u8]>,
    second: Option<&'a [u8]>,
    armed: &'a Cell<bool>,
}

impl<'a> Iterator for ArmBeforeSecond<'a> {
    type Item = &'a [u8];

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(first) = self.first.take() {
            return Some(first);
        }
        let second = self.second.take()?;
        self.armed.set(true);
        Some(second)
    }
}

#[test]
#[allow(clippy::too_many_lines)] // One hostile-length/count transactional refusal matrix.
fn hostile_lengths_counts_and_stream_mismatch_refuse_before_publication() {
    let tiny = CanonicalLimits::new(1024, 64, 4, 2, 4);
    let huge = CanonicalEncoder::<SemanticId<BytesLeaf>, _>::new(tiny, NeverCancel)
        .unwrap()
        .bytes_stream(Field::new(0, "payload"), u64::MAX, PanicChunks);
    assert!(matches!(
        huge,
        Err(CanonicalError::LimitExceeded {
            kind: LimitKind::FieldBytes,
            requested: u64::MAX,
            limit: 64
        })
    ));

    let huge_count = CanonicalEncoder::<SemanticId<SequenceLeaf>, _>::new(tiny, NeverCancel)
        .unwrap()
        .ordered_bytes(Field::new(0, "items"), u64::MAX, PanicChunks);
    assert!(matches!(
        huge_count,
        Err(CanonicalError::LimitExceeded {
            kind: LimitKind::CollectionItems,
            requested: u64::MAX,
            limit: 2
        })
    ));

    let armed = Cell::new(false);
    let oversized_set_item = [b'a'; 65];
    let set_values = ArmBeforeSecond {
        first: Some(b"a"),
        second: Some(&oversized_set_item),
        armed: &armed,
    };
    let oversized_set = CanonicalEncoder::<SemanticId<SetLeaf>, _>::new(tiny, || {
        assert!(
            !armed.get(),
            "oversized set item must refuse before ordering comparison"
        );
        false
    })
    .unwrap()
    .canonical_set(Field::new(0, "items"), 2, set_values);
    assert!(matches!(
        oversized_set,
        Err(CanonicalError::LimitExceeded {
            kind: LimitKind::FieldBytes,
            requested: 65,
            limit: 64
        })
    ));

    let short = CanonicalEncoder::<SemanticId<BytesLeaf>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .bytes_stream(Field::new(0, "payload"), 4, [b"123".as_slice()]);
    assert!(matches!(
        short,
        Err(CanonicalError::DeclaredLengthMismatch {
            declared: 4,
            observed: 3
        })
    ));
    let long = CanonicalEncoder::<SemanticId<BytesLeaf>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .bytes_stream(Field::new(0, "payload"), 3, [b"1234".as_slice()]);
    assert!(matches!(
        long,
        Err(CanonicalError::DeclaredLengthMismatch {
            declared: 3,
            observed: 4
        })
    ));

    let empty = CanonicalEncoder::<SemanticId<BytesLeaf>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .bytes_stream(Field::new(0, "payload"), 0, std::iter::empty::<&[u8]>())
        .unwrap()
        .finish();
    assert!(empty.is_ok());
    let no_stream_chunks = CanonicalLimits::new(
        LIMITS.max_canonical_bytes(),
        LIMITS.max_field_bytes(),
        LIMITS.max_fields(),
        0,
        LIMITS.cancellation_poll_bytes(),
    );
    let direct_empty =
        CanonicalEncoder::<SemanticId<BytesLeaf>, _>::new(no_stream_chunks, NeverCancel)
            .unwrap()
            .bytes(Field::new(0, "payload"), b"")
            .and_then(CanonicalEncoder::finish);
    assert!(direct_empty.is_ok());
    let too_many_empty_chunks =
        CanonicalEncoder::<SemanticId<BytesLeaf>, _>::new(tiny, NeverCancel)
            .unwrap()
            .bytes_stream(
                Field::new(0, "payload"),
                0,
                std::iter::repeat_n(b"".as_slice(), 3),
            );
    assert!(matches!(
        too_many_empty_chunks,
        Err(CanonicalError::LimitExceeded {
            kind: LimitKind::StreamChunks,
            requested: 3,
            limit: 2
        })
    ));

    let short_sequence = CanonicalEncoder::<SemanticId<SequenceLeaf>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .ordered_bytes(Field::new(0, "items"), 2, [b"a".as_slice()]);
    assert!(matches!(
        short_sequence,
        Err(CanonicalError::DeclaredLengthMismatch {
            declared: 2,
            observed: 1
        })
    ));
    let short_set = CanonicalEncoder::<SemanticId<SetLeaf>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .canonical_set(Field::new(0, "items"), 2, [b"a".as_slice()]);
    assert!(matches!(
        short_set,
        Err(CanonicalError::DeclaredLengthMismatch {
            declared: 2,
            observed: 1
        })
    ));
    let short_children =
        CanonicalEncoder::<EvidenceNodeId<ChildrenParent>, _>::new(LIMITS, NeverCancel)
            .unwrap()
            .ordered_children(Field::new(0, "children"), 2, [leaf(1).id()]);
    assert!(matches!(
        short_children,
        Err(CanonicalError::DeclaredLengthMismatch {
            declared: 2,
            observed: 1
        })
    ));
}

#[test]
#[allow(
    clippy::too_many_lines,
    reason = "one exact underwrite, overwrite, and row-count refusal matrix"
)]
fn ordered_row_stream_refuses_incomplete_excess_and_wrong_counts() {
    let underwrite = CanonicalEncoder::<SemanticId<SequenceLeaf>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .ordered_bytes_stream(
            Field::new(0, "items"),
            1,
            [Ok::<u64, RowFixtureError>(3)],
            |_, mut sink| {
                sink.write(b"ab").unwrap();
                Ok(())
            },
        );
    let Err(OrderedBytesStreamError::Canonical { source, diagnostic }) = underwrite else {
        panic!("short row must consume and refuse the encoder")
    };
    assert_eq!(
        source,
        CanonicalError::DeclaredLengthMismatch {
            declared: 3,
            observed: 2,
        }
    );
    assert_eq!(diagnostic.schema_domain(), SequenceLeaf::DOMAIN);
    assert_eq!(diagnostic.schema_name(), SequenceLeaf::NAME);
    assert_eq!(diagnostic.field_ordinal(), 0);
    assert_eq!(diagnostic.field_name(), "items");
    assert_eq!(diagnostic.phase(), OrderedBytesStreamPhase::RowCompletion);
    assert_eq!(diagnostic.row_index(), Some(0));
    assert_eq!(diagnostic.declared_rows(), 1);
    assert_eq!(diagnostic.completed_rows(), 0);
    assert_eq!(diagnostic.declared_row_bytes(), Some(3));
    assert_eq!(diagnostic.written_row_bytes(), 2);
    assert_eq!(diagnostic.prior_collection_items(), 0);
    assert_eq!(diagnostic.stream_chunks(), 1);
    assert_eq!(
        diagnostic.disposition(),
        OrderedBytesStreamDisposition::EncoderConsumedNoPublication
    );

    let overwrite = CanonicalEncoder::<SemanticId<SequenceLeaf>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .ordered_bytes_stream(
            Field::new(0, "items"),
            1,
            [Ok::<u64, RowFixtureError>(3)],
            |_, mut sink| {
                sink.write(b"ab").unwrap();
                let first = sink.write(b"cd").unwrap_err();
                assert_eq!(
                    first,
                    CanonicalError::DeclaredLengthMismatch {
                        declared: 3,
                        observed: 4,
                    }
                );
                assert_eq!(sink.write(b"ignored"), Err(first));
                Err(RowFixtureError::Producer("must-not-mask-poison"))
            },
        );
    let Err(OrderedBytesStreamError::Canonical { source, diagnostic }) = overwrite else {
        panic!("ignored overwrite error must remain sticky")
    };
    assert_eq!(
        source,
        CanonicalError::DeclaredLengthMismatch {
            declared: 3,
            observed: 4,
        }
    );
    assert_eq!(diagnostic.phase(), OrderedBytesStreamPhase::RowChunk);
    assert_eq!(diagnostic.written_row_bytes(), 2);
    assert_eq!(diagnostic.stream_chunks(), 2);

    let calls = Cell::new(0u64);
    let too_few = CanonicalEncoder::<SemanticId<SequenceLeaf>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .ordered_bytes_stream(
            Field::new(0, "items"),
            2,
            [Ok::<u64, RowFixtureError>(1)],
            |_, mut sink| {
                calls.set(calls.get() + 1);
                sink.write(b"a").unwrap();
                Ok(())
            },
        );
    assert!(matches!(
        too_few,
        Err(OrderedBytesStreamError::Canonical {
            source: CanonicalError::DeclaredLengthMismatch {
                declared: 2,
                observed: 1,
            },
            ..
        })
    ));
    assert_eq!(calls.get(), 1);

    calls.set(0);
    let too_many = CanonicalEncoder::<SemanticId<SequenceLeaf>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .ordered_bytes_stream(
            Field::new(0, "items"),
            1,
            [
                Ok::<u64, RowFixtureError>(1),
                Ok::<u64, RowFixtureError>(999),
            ],
            |_, mut sink| {
                calls.set(calls.get() + 1);
                sink.write(b"a").unwrap();
                Ok(())
            },
        );
    assert!(matches!(
        too_many,
        Err(OrderedBytesStreamError::Canonical {
            source: CanonicalError::DeclaredLengthMismatch {
                declared: 1,
                observed: 2,
            },
            ..
        })
    ));
    assert_eq!(calls.get(), 1, "an excess declaration has no row callback");
}

#[test]
fn ordered_row_stream_preserves_length_and_payload_producer_errors() {
    let calls = Cell::new(0u64);
    let zero_row_source_failure =
        CanonicalEncoder::<SemanticId<SequenceLeaf>, _>::new(LIMITS, NeverCancel)
            .unwrap()
            .ordered_bytes_stream(
                Field::new(0, "items"),
                0,
                [Err(RowFixtureError::LengthSource("zero-row-source"))],
                |_, _| -> Result<(), RowFixtureError> {
                    panic!("zero declared rows have no payload producer")
                },
            );
    assert!(matches!(
        zero_row_source_failure,
        Err(OrderedBytesStreamError::Producer {
            source: RowFixtureError::LengthSource("zero-row-source"),
            ..
        })
    ));

    let length_failure = CanonicalEncoder::<SemanticId<SequenceLeaf>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .ordered_bytes_stream(
            Field::new(0, "items"),
            2,
            [
                Ok::<u64, RowFixtureError>(1),
                Err(RowFixtureError::LengthSource("sentinel-length")),
            ],
            |_, mut sink| {
                calls.set(calls.get() + 1);
                sink.write(b"a").unwrap();
                Ok(())
            },
        );
    let Err(OrderedBytesStreamError::Producer { source, diagnostic }) = length_failure else {
        panic!("fallible row-length source must be preserved")
    };
    assert_eq!(source, RowFixtureError::LengthSource("sentinel-length"));
    assert_eq!(diagnostic.phase(), OrderedBytesStreamPhase::RowDeclaration);
    assert_eq!(diagnostic.row_index(), Some(1));
    assert_eq!(diagnostic.completed_rows(), 1);
    assert_eq!(calls.get(), 1);

    let before_payload = CanonicalEncoder::<SemanticId<SequenceLeaf>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .ordered_bytes_stream(
            Field::new(0, "items"),
            1,
            [Ok::<u64, RowFixtureError>(3)],
            |_, _| Err(RowFixtureError::Producer("before-payload")),
        );
    let Err(OrderedBytesStreamError::Producer { source, diagnostic }) = before_payload else {
        panic!("a producer may fail before its first payload chunk")
    };
    assert_eq!(source, RowFixtureError::Producer("before-payload"));
    assert_eq!(diagnostic.phase(), OrderedBytesStreamPhase::RowProducer);
    assert_eq!(diagnostic.written_row_bytes(), 0);

    let partial_failure = CanonicalEncoder::<SemanticId<SequenceLeaf>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .ordered_bytes_stream(
            Field::new(0, "items"),
            1,
            [Ok::<u64, RowFixtureError>(3)],
            |_, mut sink| {
                sink.write(b"a").unwrap();
                Err(RowFixtureError::Producer("sentinel-payload"))
            },
        );
    let Err(OrderedBytesStreamError::Producer { source, diagnostic }) = partial_failure else {
        panic!("producer error must outrank a derived underwrite")
    };
    assert_eq!(source, RowFixtureError::Producer("sentinel-payload"));
    assert_eq!(diagnostic.phase(), OrderedBytesStreamPhase::RowProducer);
    assert_eq!(diagnostic.declared_row_bytes(), Some(3));
    assert_eq!(diagnostic.written_row_bytes(), 1);
    assert!(diagnostic.canonical_bytes() > 1);
    assert_eq!(
        diagnostic.disposition(),
        OrderedBytesStreamDisposition::EncoderConsumedNoPublication
    );
}

#[test]
#[allow(clippy::too_many_lines)] // One transactional limit/refusal matrix.
fn ordered_row_stream_limits_refuse_before_hostile_producers() {
    let below_count_prefix = CanonicalLimits::new(1024, 7, 4, 4, 4);
    let eager_empty =
        CanonicalEncoder::<SemanticId<TinySequenceLeaf>, _>::new(below_count_prefix, NeverCancel)
            .unwrap()
            .ordered_bytes(Field::new(0, "i"), 0, std::iter::empty::<&'static [u8]>());
    assert!(matches!(
        eager_empty,
        Err(CanonicalError::LimitExceeded {
            kind: LimitKind::FieldBytes,
            requested: 8,
            limit: 7,
        })
    ));
    let streamed_empty =
        CanonicalEncoder::<SemanticId<TinySequenceLeaf>, _>::new(below_count_prefix, NeverCancel)
            .unwrap()
            .ordered_bytes_stream(
                Field::new(0, "i"),
                0,
                PanicLengths,
                |_, _| -> Result<(), RowFixtureError> {
                    panic!("count-prefix refusal must precede row production")
                },
            );
    assert!(matches!(
        streamed_empty,
        Err(OrderedBytesStreamError::Canonical {
            source: CanonicalError::LimitExceeded {
                kind: LimitKind::FieldBytes,
                requested: 8,
                limit: 7,
            },
            ..
        })
    ));

    let tiny = CanonicalLimits::new(64 * 1024, 64, 4, 2, 4);
    let count_refusal = CanonicalEncoder::<SemanticId<SequenceLeaf>, _>::new(tiny, NeverCancel)
        .unwrap()
        .ordered_bytes_stream(
            Field::new(0, "items"),
            3,
            PanicLengths,
            |_, _| -> Result<(), RowFixtureError> {
                panic!("count refusal must precede row production")
            },
        );
    assert!(matches!(
        count_refusal,
        Err(OrderedBytesStreamError::Canonical {
            source: CanonicalError::LimitExceeded {
                kind: LimitKind::CollectionItems,
                requested: 3,
                limit: 2,
            },
            ..
        })
    ));

    let at_exact_row_and_chunk_cap =
        CanonicalEncoder::<SemanticId<SequenceLeaf>, _>::new(tiny, NeverCancel)
            .unwrap()
            .ordered_bytes_stream(
                Field::new(0, "items"),
                2,
                [Ok::<u64, RowFixtureError>(0), Ok::<u64, RowFixtureError>(0)],
                |_, mut sink| {
                    sink.write(b"").unwrap();
                    Ok(())
                },
            )
            .unwrap()
            .finish()
            .unwrap();
    assert_eq!(at_exact_row_and_chunk_cap.collection_items(), 2);

    let producer_called = Cell::new(false);
    let item_refusal = CanonicalEncoder::<SemanticId<SequenceLeaf>, _>::new(tiny, NeverCancel)
        .unwrap()
        .ordered_bytes_stream(
            Field::new(0, "items"),
            1,
            [Ok::<u64, RowFixtureError>(65)],
            |_, _| {
                producer_called.set(true);
                Ok(())
            },
        );
    assert!(matches!(
        item_refusal,
        Err(OrderedBytesStreamError::Canonical {
            source: CanonicalError::LimitExceeded {
                kind: LimitKind::FieldBytes,
                requested: 65,
                limit: 64,
            },
            ..
        })
    ));
    assert!(!producer_called.get());

    let aggregate_refusal = CanonicalEncoder::<SemanticId<SequenceLeaf>, _>::new(tiny, NeverCancel)
        .unwrap()
        .ordered_bytes_stream(
            Field::new(0, "items"),
            1,
            [Ok::<u64, RowFixtureError>(50)],
            |_, _| {
                producer_called.set(true);
                Ok(())
            },
        );
    assert!(matches!(
        aggregate_refusal,
        Err(OrderedBytesStreamError::Canonical {
            source: CanonicalError::LimitExceeded {
                kind: LimitKind::FieldBytes,
                requested: 66,
                limit: 64,
            },
            ..
        })
    ));
    assert!(!producer_called.get());

    let arithmetic_limits = CanonicalLimits::new(u64::MAX, u64::MAX, 4, 2, 4);
    let arithmetic_refusal =
        CanonicalEncoder::<SemanticId<SequenceLeaf>, _>::new(arithmetic_limits, NeverCancel)
            .unwrap()
            .ordered_bytes_stream(
                Field::new(0, "items"),
                1,
                [Ok::<u64, RowFixtureError>(u64::MAX)],
                |_, _| {
                    producer_called.set(true);
                    Ok(())
                },
            );
    assert!(matches!(
        arithmetic_refusal,
        Err(OrderedBytesStreamError::Canonical {
            source: CanonicalError::LengthOverflow,
            ..
        })
    ));
    assert!(!producer_called.get());

    let chunk_refusal = CanonicalEncoder::<SemanticId<SequenceLeaf>, _>::new(tiny, NeverCancel)
        .unwrap()
        .ordered_bytes_stream(
            Field::new(0, "items"),
            1,
            [Ok::<u64, RowFixtureError>(0)],
            |_, mut sink| {
                sink.write(b"").unwrap();
                sink.write(b"").unwrap();
                let source = sink.write(b"").unwrap_err();
                assert!(matches!(
                    source,
                    CanonicalError::LimitExceeded {
                        kind: LimitKind::StreamChunks,
                        requested: 3,
                        limit: 2,
                    }
                ));
                Ok(())
            },
        );
    let Err(OrderedBytesStreamError::Canonical { source, diagnostic }) = chunk_refusal else {
        panic!("ignored chunk-limit error must poison the encoder")
    };
    assert!(matches!(
        source,
        CanonicalError::LimitExceeded {
            kind: LimitKind::StreamChunks,
            requested: 3,
            limit: 2,
        }
    ));
    assert_eq!(diagnostic.stream_chunks(), 3);
    assert_eq!(diagnostic.written_row_bytes(), 0);

    let payload = [7u8; 40];
    let baseline_rows: [&[u8]; 1] = [&payload];
    let baseline = streamed_sequence(&baseline_rows, LIMITS, 0);
    let prefinish_bytes = baseline.canonical_bytes() - 5;
    let exact = CanonicalLimits::new(prefinish_bytes, 128, 4, 2, 4);
    let exact_field = CanonicalEncoder::<SemanticId<SequenceLeaf>, _>::new(exact, NeverCancel)
        .unwrap()
        .ordered_bytes_stream(
            Field::new(0, "items"),
            1,
            [Ok::<u64, RowFixtureError>(payload.len() as u64)],
            |_, mut sink| {
                sink.write(&payload).unwrap();
                Ok(())
            },
        );
    assert!(exact_field.is_ok());

    let below = CanonicalLimits::new(prefinish_bytes - 1, 128, 4, 2, 4);
    let frame_refusal = CanonicalEncoder::<SemanticId<SequenceLeaf>, _>::new(below, NeverCancel)
        .unwrap()
        .ordered_bytes_stream(
            Field::new(0, "items"),
            1,
            [Ok::<u64, RowFixtureError>(payload.len() as u64)],
            |_, _| {
                producer_called.set(true);
                Ok(())
            },
        );
    assert!(matches!(
        frame_refusal,
        Err(OrderedBytesStreamError::Canonical {
            source: CanonicalError::LimitExceeded {
                kind: LimitKind::CanonicalBytes,
                ..
            },
            ..
        })
    ));
    assert!(!producer_called.get());
}

#[test]
fn schema_descriptor_floor_is_typed_recursive_and_identity_neutral() {
    assert_descriptor_floor_boundaries::<DescriptorDomainFloor>(SchemaDescriptorKind::Domain);
    assert_descriptor_floor_boundaries::<DescriptorNameFloor>(SchemaDescriptorKind::Name);
    assert_descriptor_floor_boundaries::<DescriptorContextFloor>(SchemaDescriptorKind::Context);
    assert_descriptor_floor_boundaries::<DescriptorFieldNameFloor>(SchemaDescriptorKind::FieldName);

    let recursive_floor = u64::try_from(DescriptorNestedChildFloor::CONTEXT.len())
        .expect("static child descriptor length fits u64");
    let supplied = recursive_floor
        .checked_sub(1)
        .expect("recursive descriptor floor is positive");
    let recursive_error = CanonicalEncoder::<SemanticId<DescriptorNestedParent>, _>::new(
        limits_with_max_field_bytes(supplied),
        NeverCancel,
    )
    .expect_err("one-under a nested descriptor floor must refuse construction");
    assert_eq!(
        recursive_error,
        CanonicalError::MaxFieldBytesBelowSchemaFloor {
            binding: SchemaDescriptorKind::Context,
            floor: recursive_floor,
            supplied,
        }
    );
    CanonicalEncoder::<SemanticId<DescriptorNestedParent>, _>::new(
        limits_with_max_field_bytes(recursive_floor),
        NeverCancel,
    )
    .expect("the exact recursive descriptor floor admits construction");
}

#[test]
fn field_bytes_limit_exceeded_is_reserved_for_payloads_after_schema_admission() {
    let (_, descriptor_floor) = top_level_descriptor_floor::<BytesLeaf>();
    let payload_bytes = descriptor_floor
        .checked_add(1)
        .expect("fixture descriptor floor plus one fits u64");
    let payload_len = usize::try_from(payload_bytes).expect("fixture payload length fits usize");
    let payload = vec![0x5au8; payload_len];
    let refusal = CanonicalEncoder::<SemanticId<BytesLeaf>, _>::new(
        limits_with_max_field_bytes(descriptor_floor),
        NeverCancel,
    )
    .expect("the exact schema descriptor floor admits construction")
    .bytes(Field::new(0, "payload"), &payload);
    assert!(matches!(
        refusal,
        Err(CanonicalError::LimitExceeded {
            kind: LimitKind::FieldBytes,
            requested,
            limit,
        }) if requested == payload_bytes && limit == descriptor_floor
    ));
}

#[test]
fn invalid_parent_and_child_schema_descriptors_refuse_before_publication() {
    let invalid_limits = CanonicalLimits::new(1024, 256, 4, 4, 0);
    assert!(matches!(
        CanonicalEncoder::<SemanticId<LeafV1>, _>::new(invalid_limits, NeverCancel),
        Err(CanonicalError::InvalidLimits(_))
    ));
    let empty_domain =
        CanonicalEncoder::<SemanticId<EmptyDomainSchema>, _>::new(LIMITS, NeverCancel);
    assert!(matches!(
        empty_domain,
        Err(CanonicalError::InvalidSchemaDescriptor(_))
    ));
    let zero_version =
        CanonicalEncoder::<SemanticId<ZeroVersionSchema>, _>::new(LIMITS, NeverCancel);
    assert!(matches!(
        zero_version,
        Err(CanonicalError::InvalidSchemaDescriptor(_))
    ));
    let duplicate_fields =
        CanonicalEncoder::<SemanticId<DuplicateFieldSchema>, _>::new(LIMITS, NeverCancel);
    assert!(matches!(
        duplicate_fields,
        Err(CanonicalError::InvalidSchemaDescriptor(_))
    ));

    let invalid_child =
        SemanticId::<EmptyDomainSchema>::parse_slice(&[0x5au8; 32]).expect("typed parse");
    let nested = CanonicalEncoder::<SemanticId<ChildParent>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .child(Field::new(0, "child"), invalid_child);
    assert!(matches!(
        nested,
        Err(CanonicalError::InvalidSchemaDescriptor(_))
    ));

    CanonicalEncoder::<SemanticId<NestedParentSchema>, _>::new(LIMITS, NeverCancel)
        .expect("valid recursive child descriptors admit");
    let one_field_per_descriptor = CanonicalLimits::new(
        LIMITS.max_canonical_bytes(),
        LIMITS.max_field_bytes(),
        1,
        LIMITS.max_collection_items(),
        LIMITS.cancellation_poll_bytes(),
    );
    let nested_field_cap = CanonicalEncoder::<SemanticId<NestedParentSchema>, _>::new(
        one_field_per_descriptor,
        NeverCancel,
    );
    assert!(matches!(
        nested_field_cap,
        Err(CanonicalError::LimitExceeded {
            kind: LimitKind::Fields,
            requested: 2,
            limit: 1,
        })
    ));

    let invalid_grandchild =
        CanonicalEncoder::<SemanticId<InvalidNestedParentSchema>, _>::new(LIMITS, NeverCancel);
    assert!(matches!(
        invalid_grandchild,
        Err(CanonicalError::InvalidSchemaDescriptor(
            "domain, schema name, and context must be non-empty"
        ))
    ));

    CanonicalEncoder::<SemanticId<Depth01Schema>, _>::new(LIMITS, NeverCancel)
        .expect("the exact sixteen-child-edge descriptor boundary admits");
    let over_depth =
        CanonicalEncoder::<SemanticId<OverDepthParentSchema>, _>::new(LIMITS, NeverCancel);
    assert!(matches!(
        over_depth,
        Err(CanonicalError::InvalidSchemaDescriptor(
            "child schema nesting exceeds the supported depth"
        ))
    ));

    let bounded_expansion = CanonicalLimits::new(
        LIMITS.max_canonical_bytes(),
        LIMITS.max_field_bytes(),
        2,
        LIMITS.max_collection_items(),
        LIMITS.cancellation_poll_bytes(),
    );
    let branch_expansion =
        CanonicalEncoder::<SemanticId<Branch05Schema>, _>::new(bounded_expansion, NeverCancel);
    assert!(matches!(
        branch_expansion,
        Err(CanonicalError::LimitExceeded {
            kind: LimitKind::Fields,
            requested: 36,
            limit: 34,
        })
    ));
}

#[test]
fn budgets_do_not_move_an_admitted_identity() {
    let baseline = leaf(99);
    for limits in [
        CanonicalLimits::new(
            baseline.canonical_bytes(),
            LIMITS.max_field_bytes(),
            LIMITS.max_fields(),
            LIMITS.max_collection_items(),
            LIMITS.cancellation_poll_bytes(),
        ),
        CanonicalLimits::new(
            LIMITS.max_canonical_bytes(),
            64,
            LIMITS.max_fields(),
            LIMITS.max_collection_items(),
            LIMITS.cancellation_poll_bytes(),
        ),
        CanonicalLimits::new(
            LIMITS.max_canonical_bytes(),
            LIMITS.max_field_bytes(),
            1,
            LIMITS.max_collection_items(),
            LIMITS.cancellation_poll_bytes(),
        ),
        CanonicalLimits::new(
            LIMITS.max_canonical_bytes(),
            LIMITS.max_field_bytes(),
            LIMITS.max_fields(),
            0,
            LIMITS.cancellation_poll_bytes(),
        ),
        CanonicalLimits::new(
            LIMITS.max_canonical_bytes(),
            LIMITS.max_field_bytes(),
            LIMITS.max_fields(),
            LIMITS.max_collection_items(),
            1,
        ),
    ] {
        let admitted = CanonicalEncoder::<SemanticId<LeafV1>, _>::new(limits, NeverCancel)
            .unwrap()
            .u64(Field::new(0, "value"), 99)
            .unwrap()
            .finish()
            .unwrap();
        assert_eq!(baseline.id(), admitted.id());
        assert_eq!(baseline.canonical_preimage(), admitted.canonical_preimage());
        assert_eq!(baseline.canonical_bytes(), admitted.canonical_bytes());
    }
    let below = CanonicalLimits::new(
        baseline.canonical_bytes() - 1,
        LIMITS.max_field_bytes(),
        LIMITS.max_fields(),
        LIMITS.max_collection_items(),
        LIMITS.cancellation_poll_bytes(),
    );
    let refused = CanonicalEncoder::<SemanticId<LeafV1>, _>::new(below, NeverCancel)
        .and_then(|encoder| encoder.u64(Field::new(0, "value"), 99))
        .and_then(CanonicalEncoder::finish);
    assert!(matches!(
        refused,
        Err(CanonicalError::LimitExceeded {
            kind: LimitKind::CanonicalBytes,
            ..
        })
    ));
}

#[derive(Clone)]
struct CountProbe {
    calls: Rc<Cell<usize>>,
    cancel_at: Option<usize>,
}

impl fs_blake3::identity::CancellationProbe for CountProbe {
    fn is_cancelled(&mut self) -> bool {
        let next = self.calls.get() + 1;
        self.calls.set(next);
        self.cancel_at == Some(next)
    }
}

fn cancellable_leaf(
    probe: CountProbe,
) -> Result<fs_blake3::identity::IdentityReceipt<SemanticId<LeafV1>>, CanonicalError> {
    CanonicalEncoder::<SemanticId<LeafV1>, _>::new(LIMITS, probe)?
        .u64(Field::new(0, "value"), 123)?
        .finish()
}

fn cancellable_all_schema(probe: CountProbe) -> Result<(), CanonicalError> {
    CanonicalEncoder::<SemanticId<AllFields>, _>::new(LIMITS, probe)?;
    Ok(())
}

fn cancellable_nested_schema(probe: CountProbe) -> Result<(), CanonicalError> {
    CanonicalEncoder::<SemanticId<NestedParentSchema>, _>::new(LIMITS, probe)?;
    Ok(())
}

fn cancellable_equal_prefix_set(
    probe: CountProbe,
) -> Result<fs_blake3::identity::IdentityReceipt<SemanticId<SetLeaf>>, CanonicalError> {
    let first = [b'a'; 256];
    let mut second = first;
    second[255] = b'b';
    CanonicalEncoder::<SemanticId<SetLeaf>, _>::new(LIMITS, probe)?
        .canonical_set(
            Field::new(0, "items"),
            2,
            [first.as_slice(), second.as_slice()],
        )?
        .finish()
}

fn cancellable_ordered_rows(
    probe: CountProbe,
) -> Result<fs_blake3::identity::IdentityReceipt<SemanticId<SequenceLeaf>>, CanonicalError> {
    let rows: [&[u8]; 3] = [b"abcdefghijk", b"", b"lmnopqrstuvwxyz"];
    let lengths = rows
        .iter()
        .map(|row| Ok::<u64, CanonicalError>(row.len() as u64));
    let encoder = CanonicalEncoder::<SemanticId<SequenceLeaf>, _>::new(LIMITS, probe)?;
    let encoder = encoder
        .ordered_bytes_stream(
            Field::new(0, "items"),
            rows.len() as u64,
            lengths,
            |row_index, mut sink| -> Result<(), CanonicalError> {
                sink.write(b"")?;
                for chunk in rows[row_index as usize].chunks(3) {
                    sink.write(chunk)?;
                }
                Ok(())
            },
        )
        .map_err(|error| match error {
            OrderedBytesStreamError::Canonical { source, .. }
            | OrderedBytesStreamError::Producer { source, .. } => source,
        })?;
    encoder.finish()
}

#[test]
fn cancellation_at_every_checkpoint_publishes_no_partial_identity() {
    let calls = Rc::new(Cell::new(0));
    let baseline = cancellable_leaf(CountProbe {
        calls: Rc::clone(&calls),
        cancel_at: None,
    })
    .unwrap();
    let total_calls = calls.get();
    assert!(total_calls > 4);
    for cancel_at in 1..=total_calls {
        let calls = Rc::new(Cell::new(0));
        let result = cancellable_leaf(CountProbe {
            calls,
            cancel_at: Some(cancel_at),
        });
        assert!(matches!(result, Err(CanonicalError::Cancelled { .. })));
    }
    let retry = cancellable_leaf(CountProbe {
        calls: Rc::new(Cell::new(0)),
        cancel_at: None,
    })
    .unwrap();
    assert_eq!(baseline, retry);
}

#[test]
fn cancellation_covers_schema_validation_and_long_set_comparisons() {
    let schema_calls = Rc::new(Cell::new(0));
    cancellable_all_schema(CountProbe {
        calls: Rc::clone(&schema_calls),
        cancel_at: None,
    })
    .unwrap();
    for cancel_at in 1..=schema_calls.get() {
        let result = cancellable_all_schema(CountProbe {
            calls: Rc::new(Cell::new(0)),
            cancel_at: Some(cancel_at),
        });
        assert!(matches!(result, Err(CanonicalError::Cancelled { .. })));
    }

    let nested_calls = Rc::new(Cell::new(0));
    cancellable_nested_schema(CountProbe {
        calls: Rc::clone(&nested_calls),
        cancel_at: None,
    })
    .unwrap();
    for cancel_at in 1..=nested_calls.get() {
        let result = cancellable_nested_schema(CountProbe {
            calls: Rc::new(Cell::new(0)),
            cancel_at: Some(cancel_at),
        });
        assert!(matches!(result, Err(CanonicalError::Cancelled { .. })));
    }

    let set_calls = Rc::new(Cell::new(0));
    cancellable_equal_prefix_set(CountProbe {
        calls: Rc::clone(&set_calls),
        cancel_at: None,
    })
    .unwrap();
    for cancel_at in 1..=set_calls.get() {
        let result = cancellable_equal_prefix_set(CountProbe {
            calls: Rc::new(Cell::new(0)),
            cancel_at: Some(cancel_at),
        });
        assert!(matches!(result, Err(CanonicalError::Cancelled { .. })));
    }
}

#[test]
fn cancellation_at_every_ordered_row_checkpoint_publishes_nothing() {
    let calls = Rc::new(Cell::new(0));
    let baseline = cancellable_ordered_rows(CountProbe {
        calls: Rc::clone(&calls),
        cancel_at: None,
    })
    .unwrap();
    let total_calls = calls.get();
    assert!(total_calls > 10);

    for cancel_at in 1..=total_calls {
        let result = cancellable_ordered_rows(CountProbe {
            calls: Rc::new(Cell::new(0)),
            cancel_at: Some(cancel_at),
        });
        assert!(matches!(result, Err(CanonicalError::Cancelled { .. })));
    }

    let retry = cancellable_ordered_rows(CountProbe {
        calls: Rc::new(Cell::new(0)),
        cancel_at: None,
    })
    .unwrap();
    assert_eq!(retry, baseline);
}

#[test]
fn synthetic_collision_adjudication_preserves_both_observations() {
    let synthetic = SemanticId::<LeafV1>::parse_hex(&"11".repeat(32)).unwrap();
    let first_bytes = ByteObservation::new(ContentId::of_bytes(b"first"), 5);
    let same_root_other_length = ByteObservation::new(first_bytes.content_id(), 6);
    let other_root_same_length = ByteObservation::new(ContentId::of_bytes(b"second"), 5);
    let first = ObservedIdentity::presented(synthetic, first_bytes);
    let second = ObservedIdentity::presented(synthetic, same_root_other_length);
    let IdentityAdjudication::Refused(length_refusal) = adjudicate(first, second) else {
        panic!("same typed ID with different bytes must refuse")
    };
    assert_eq!(length_refusal.id(), synthetic);
    assert_eq!(length_refusal.first(), first_bytes);
    assert_eq!(length_refusal.second(), same_root_other_length);

    let third = ObservedIdentity::presented(synthetic, other_root_same_length);
    let IdentityAdjudication::Refused(root_refusal) = adjudicate(first, third) else {
        panic!("same typed ID with a different byte root must refuse")
    };
    assert_eq!(root_refusal.id(), synthetic);
    assert_eq!(root_refusal.first(), first_bytes);
    assert_eq!(root_refusal.second(), other_root_same_length);
    assert_eq!(
        adjudicate(first, first),
        IdentityAdjudication::SameObservation
    );
    let distinct = ObservedIdentity::presented(leaf(7).id(), first_bytes);
    assert_eq!(
        adjudicate(first, distinct),
        IdentityAdjudication::DistinctIds
    );
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum AuthorityError {
    Mismatch,
}

#[derive(Clone, Copy)]
struct ExactVerifier {
    subject: SemanticId<LeafV1>,
    preimage: ContentId,
    anchor: ExternalAnchorRef,
    verifier: VerifierId<VerifierSchema>,
    policy: KeyPolicyId<PolicySchema>,
}

impl AuthorityVerifier<SemanticId<LeafV1>, VerifierSchema, PolicySchema> for ExactVerifier {
    type Error = AuthorityError;

    fn verify(
        &self,
        presented: &AuthorityRef<
            SemanticId<LeafV1>,
            VerifierSchema,
            PolicySchema,
            fs_blake3::identity::Presented,
        >,
    ) -> Result<(), Self::Error> {
        let receipt = presented.receipt();
        if receipt.id() == self.subject
            && receipt.canonical_preimage() == self.preimage
            && presented.anchor() == self.anchor
            && presented.verifier() == self.verifier
            && presented.key_policy() == self.policy
        {
            Ok(())
        } else {
            Err(AuthorityError::Mismatch)
        }
    }
}

#[derive(Clone, Copy)]
struct ExactAdmitter {
    subject: SemanticId<LeafV1>,
    preimage: ContentId,
    anchor: ExternalAnchorRef,
    verifier: VerifierId<VerifierSchema>,
    policy: KeyPolicyId<PolicySchema>,
}

impl AuthorityAdmitter<SemanticId<LeafV1>, VerifierSchema, PolicySchema> for ExactAdmitter {
    type Error = AuthorityError;

    fn admit(
        &self,
        verified: &AuthorityRef<
            SemanticId<LeafV1>,
            VerifierSchema,
            PolicySchema,
            fs_blake3::identity::Verified,
        >,
    ) -> Result<(), Self::Error> {
        let receipt = verified.receipt();
        if receipt.id() == self.subject
            && receipt.canonical_preimage() == self.preimage
            && verified.anchor() == self.anchor
            && verified.verifier() == self.verifier
            && verified.key_policy() == self.policy
        {
            Ok(())
        } else {
            Err(AuthorityError::Mismatch)
        }
    }
}

fn verifier_id() -> fs_blake3::identity::IdentityReceipt<VerifierId<VerifierSchema>> {
    CanonicalEncoder::<VerifierId<VerifierSchema>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .u64(Field::new(0, "version"), 1)
        .unwrap()
        .finish()
        .unwrap()
}

fn policy_id() -> fs_blake3::identity::IdentityReceipt<KeyPolicyId<PolicySchema>> {
    CanonicalEncoder::<KeyPolicyId<PolicySchema>, _>::new(LIMITS, NeverCancel)
        .unwrap()
        .u64(Field::new(0, "version"), 1)
        .unwrap()
        .finish()
        .unwrap()
}

#[test]
fn authority_presence_verification_and_admission_are_distinct() {
    let subject = leaf(33);
    let verifier = verifier_id().id();
    let policy = policy_id().id();
    let anchor = ExternalAnchorRef::presented(ContentId::of_bytes(b"external anchor"));
    let presented = AuthorityRef::present(subject, anchor, verifier, policy);
    assert_eq!(presented.trust_state(), TrustState::Presented);
    assert_eq!(
        presented.audit_record().no_claim(),
        NoClaimState::ExternalTrustRequired
    );

    let capability = ExactVerifier {
        subject: subject.id(),
        preimage: subject.canonical_preimage(),
        anchor,
        verifier,
        policy,
    };
    let verified = presented.verify(&capability).unwrap();
    assert_eq!(verified.trust_state(), TrustState::Verified);
    let admitter = ExactAdmitter {
        subject: subject.id(),
        preimage: subject.canonical_preimage(),
        anchor,
        verifier,
        policy,
    };
    let admitted = verified.admit(&admitter).unwrap();
    assert_eq!(admitted.trust_state(), TrustState::Admitted);
    let log = admitted.audit_record();
    assert_eq!(log.id(), *subject.id().as_bytes());
    assert_eq!(log.domain(), LeafV1::DOMAIN);
    assert_eq!(log.schema_name(), LeafV1::NAME);
    assert_eq!(log.version(), LeafV1::VERSION);
    assert_eq!(log.context(), LeafV1::CONTEXT);
    assert_eq!(log.anchor(), Some(anchor.content_id()));
    assert_eq!(log.verifier(), Some(*verifier.as_bytes()));
    assert_eq!(log.key_policy(), Some(*policy.as_bytes()));
    assert_eq!(log.trust(), TrustState::Admitted);
    assert_eq!(log.no_claim(), NoClaimState::ScientificCorrectnessNotProven);

    let presented = || AuthorityRef::present(subject, anchor, verifier, policy);
    let wrong_subject = ExactVerifier {
        subject: leaf(34).id(),
        ..capability
    };
    assert!(matches!(
        presented().verify(&wrong_subject),
        Err(AuthorityError::Mismatch)
    ));
    let wrong_preimage = ExactVerifier {
        preimage: ContentId::of_bytes(b"wrong preimage"),
        ..capability
    };
    assert!(matches!(
        presented().verify(&wrong_preimage),
        Err(AuthorityError::Mismatch)
    ));
    let wrong_anchor = ExactVerifier {
        anchor: ExternalAnchorRef::presented(ContentId::of_bytes(b"wrong anchor")),
        ..capability
    };
    assert!(matches!(
        presented().verify(&wrong_anchor),
        Err(AuthorityError::Mismatch)
    ));
    let wrong_verifier = ExactVerifier {
        verifier: VerifierId::<VerifierSchema>::parse_slice(&[0x22; 32]).unwrap(),
        ..capability
    };
    assert!(matches!(
        presented().verify(&wrong_verifier),
        Err(AuthorityError::Mismatch)
    ));
    let wrong_policy = ExactVerifier {
        policy: KeyPolicyId::<PolicySchema>::parse_slice(&[0x33; 32]).unwrap(),
        ..capability
    };
    assert!(matches!(
        presented().verify(&wrong_policy),
        Err(AuthorityError::Mismatch)
    ));

    let verified = presented().verify(&capability).unwrap();
    let refusing_admitter = ExactAdmitter {
        policy: KeyPolicyId::<PolicySchema>::parse_slice(&[0x44; 32]).unwrap(),
        ..admitter
    };
    assert!(matches!(
        verified.admit(&refusing_admitter),
        Err(AuthorityError::Mismatch)
    ));
}

#[test]
fn bounded_audit_records_retain_metadata_without_payloads() {
    let receipt = all_baseline();
    let record = receipt.audit_record();
    assert_eq!(record.id(), *receipt.id().as_bytes());
    assert_eq!(record.canonical_preimage(), receipt.canonical_preimage());
    assert_eq!(record.role(), IdentityRole::Semantic);
    assert_eq!(record.schema_id(), *receipt.schema_id().as_bytes());
    assert_eq!(record.canonical_bytes(), receipt.canonical_bytes());
    assert_eq!(record.field_count(), AllFields::FIELDS.len() as u32);
    assert_eq!(record.collection_items(), 6);
    assert_eq!(record.limits(), LIMITS);
    assert_eq!(record.trust(), TrustState::Unanchored);
    assert_eq!(record.anchor(), None);
    assert_eq!(record.verifier(), None);
    assert_eq!(record.key_policy(), None);
    let rendered = format!("{record:?}");
    assert!(!rendered.contains("a|b"));
    assert!(!rendered.contains("variant"));
    assert!(!rendered.contains("present"));
}

#[test]
fn display_and_debug_are_not_hash_inputs() {
    let receipt = leaf(77);
    let before = receipt.id();
    let _display = receipt.id().to_string();
    let _debug = format!("{:?}", receipt.id());
    let _audit = format!("{:?}", receipt.audit_record());
    assert_eq!(before, receipt.id());
}

#[test]
fn legacy_values_remain_exact_and_quarantined() {
    let old = LegacyProvenanceV1::new(0x0123_4567_89ab_cdef);
    assert_eq!(old.value(), 0x0123_4567_89ab_cdef);
    assert_ne!(
        format!("{old:?}"),
        leaf(old.value()).id().to_string(),
        "legacy replay value is not reinterpreted as a strong ID"
    );
}

#[test]
fn strict_typed_parsing_preserves_role_without_adding_trust() {
    let receipt = leaf(55);
    let parsed = SemanticId::<LeafV1>::parse_slice(receipt.id().as_bytes()).unwrap();
    assert_eq!(parsed, receipt.id());
    assert!(SemanticId::<LeafV1>::parse_slice(&[0u8; 31]).is_none());
    assert_eq!(
        SemanticId::<LeafV1>::parse_hex(&receipt.id().to_hex()),
        Some(receipt.id())
    );
    let raw = ContentId::parse_hex(&receipt.canonical_preimage().to_hex()).unwrap();
    assert_eq!(raw, receipt.canonical_preimage());
    assert_eq!(receipt.audit_record().trust(), TrustState::Unanchored);
}

#[test]
fn compatibility_content_hash_remains_source_compatible_and_untyped() {
    let compatibility = hash_bytes(b"compatibility");
    let direct = ContentHash(*compatibility.as_bytes());
    assert_eq!(compatibility, direct);
    assert_eq!(
        *ContentId::of_bytes(b"compatibility").as_bytes(),
        *compatibility.as_bytes()
    );
}
