//! Promotion-authority battery (bead `sj31i.52.9.1`). Generic foreign
//! capabilities may still drive policy-relative
//! `Presented → Verified → Admitted`, but public IDs/observations and a
//! configuration charter never authenticate which code ran. V3 promotion
//! requires a non-copyable root to own and execute its verifier/admitter, burn
//! replay state before calls, and bind raw evidence back to the exact live root
//! (or execute an explicit predecessor-bound crosswalk).

use fs_blake3::identity::{
    AuthorityAdmitter, AuthorityRef, AuthorityVerifier, ByteObservation, CanonicalEncoder,
    CanonicalLimits, CanonicalSchema, ChildSpec, ContentId, ExternalAnchorRef, Field, FieldSpec,
    IdentityReceipt, IdentityRole, KeyPolicyId, LEGACY_PROMOTION_ROOT_CHARTER_V1_DOMAIN,
    LEGACY_PROMOTION_ROOT_CHARTER_V2_DOMAIN, MAX_PROMOTION_CONTEXT_BYTES, NeverCancel,
    ObservedIdentity, OwnerPromotionAdmitter, OwnerPromotionCapabilities, OwnerPromotionVerifier,
    PROMOTION_ADMISSION_DECISION_IDENTITY_DOMAIN, PROMOTION_ADMISSION_DECISION_IDENTITY_VERSION,
    PROMOTION_DECISION_IDENTITY_DOMAIN, PROMOTION_DECISION_IDENTITY_VERSION,
    PROMOTION_REQUEST_IDENTITY_DOMAIN, PROMOTION_REQUEST_IDENTITY_VERSION,
    PROMOTION_ROOT_CHARTER_DOMAIN, PROMOTION_ROOT_CHARTER_IDENTITY_VERSION,
    PROMOTION_VERIFICATION_DECISION_IDENTITY_DOMAIN,
    PROMOTION_VERIFICATION_DECISION_IDENTITY_VERSION, Presented, PromotionAdmissionRequest,
    PromotionCapabilityDescriptor, PromotionCapabilityStage, PromotionCapabilityVerdict,
    PromotionDecisionDisposition, PromotionDecisionRequest, PromotionDecisionScope,
    PromotionRefusal, PromotionTrustRoot, PromotionWitness, SchemaId, SemanticId, StrongIdentity,
    Verified, VerifierId, WireType, legacy,
};
use fs_blake3::{ContentHash, hash_bytes, hash_domain};

const LIMITS: CanonicalLimits = CanonicalLimits::new(64 * 1024, 16 * 1024, 64, 1024, 7);
const EXPECTED_CURRENT_CHARTER_VERSION: u32 = 3;
const EXPECTED_CURRENT_CHARTER_DOMAIN: &str = "org.frankensim.fs-blake3.promotion-root-charter.v3";
const EXPECTED_LEGACY_CHARTER_V2_DOMAIN: &str =
    "org.frankensim.fs-blake3.promotion-root-charter.v2";
const EXPECTED_LEGACY_CHARTER_V1_DOMAIN: &str =
    "org.frankensim.fs-blake3.promotion-root-charter.v1";

struct SubjectV1;
impl CanonicalSchema for SubjectV1 {
    const DOMAIN: &'static str = "org.frankensim.test.promotion.subject.v1";
    const NAME: &'static str = "promotion-test-subject";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "G0 promotion subject fixture";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required("value", WireType::U64)];
}

struct VerifierSchemaV1;
impl CanonicalSchema for VerifierSchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.test.promotion.verifier.v1";
    const NAME: &'static str = "promotion-test-verifier";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "G0 promotion verifier fixture";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required("key", WireType::U64)];
}

struct PolicySchemaV1;
impl CanonicalSchema for PolicySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.test.promotion.policy.v1";
    const NAME: &'static str = "promotion-test-policy";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "G0 promotion policy fixture";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required("rule", WireType::U64)];
}

const CHARTER_SHARED_DOMAIN: &str = "org.frankensim.test.promotion.charter-shared.v1";

struct CharterLeafV1;
impl CanonicalSchema for CharterLeafV1 {
    const DOMAIN: &'static str = "org.frankensim.test.promotion.charter-leaf.v1";
    const NAME: &'static str = "promotion-charter-leaf";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "G0 promotion charter leaf";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required("leaf", WireType::U64)];
}

struct CharterLeafBindingVariant;
impl CanonicalSchema for CharterLeafBindingVariant {
    const DOMAIN: &'static str = CharterLeafV1::DOMAIN;
    const NAME: &'static str = CharterLeafV1::NAME;
    const VERSION: u32 = CharterLeafV1::VERSION;
    const CONTEXT: &'static str = CharterLeafV1::CONTEXT;
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required("other-leaf", WireType::U64)];
}

static CHARTER_LEAF_V1: ChildSpec = ChildSpec::for_identity::<SemanticId<CharterLeafV1>>();
static CHARTER_LEAF_BINDING_VARIANT: ChildSpec =
    ChildSpec::for_identity::<SemanticId<CharterLeafBindingVariant>>();

const CHARTER_FIELDS_V1: &[FieldSpec] = &[
    FieldSpec::required("value", WireType::U64),
    FieldSpec::child_of("child", &CHARTER_LEAF_V1),
];
const CHARTER_FIELDS_VARIANT: &[FieldSpec] = &[
    FieldSpec::required("changed-value", WireType::U64),
    FieldSpec::child_of("child", &CHARTER_LEAF_V1),
];
const CHARTER_CHILD_BINDING_VARIANT: &[FieldSpec] = &[
    FieldSpec::required("value", WireType::U64),
    FieldSpec::child_of("child", &CHARTER_LEAF_BINDING_VARIANT),
];

struct CharterSchemaV1;
impl CanonicalSchema for CharterSchemaV1 {
    const DOMAIN: &'static str = CHARTER_SHARED_DOMAIN;
    const NAME: &'static str = "promotion-charter-schema";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "G0 promotion charter schema";
    const FIELDS: &'static [FieldSpec] = CHARTER_FIELDS_V1;
}

struct CharterDomainVariant;
impl CanonicalSchema for CharterDomainVariant {
    const DOMAIN: &'static str = "org.frankensim.test.promotion.charter-other-domain.v1";
    const NAME: &'static str = CharterSchemaV1::NAME;
    const VERSION: u32 = CharterSchemaV1::VERSION;
    const CONTEXT: &'static str = CharterSchemaV1::CONTEXT;
    const FIELDS: &'static [FieldSpec] = CHARTER_FIELDS_V1;
}

struct CharterNameVariant;
impl CanonicalSchema for CharterNameVariant {
    const DOMAIN: &'static str = CHARTER_SHARED_DOMAIN;
    const NAME: &'static str = "promotion-charter-schema-other-name";
    const VERSION: u32 = CharterSchemaV1::VERSION;
    const CONTEXT: &'static str = CharterSchemaV1::CONTEXT;
    const FIELDS: &'static [FieldSpec] = CHARTER_FIELDS_V1;
}

struct CharterVersionVariant;
impl CanonicalSchema for CharterVersionVariant {
    const DOMAIN: &'static str = CHARTER_SHARED_DOMAIN;
    const NAME: &'static str = CharterSchemaV1::NAME;
    const VERSION: u32 = 2;
    const CONTEXT: &'static str = CharterSchemaV1::CONTEXT;
    const FIELDS: &'static [FieldSpec] = CHARTER_FIELDS_V1;
}

struct CharterContextVariant;
impl CanonicalSchema for CharterContextVariant {
    const DOMAIN: &'static str = CHARTER_SHARED_DOMAIN;
    const NAME: &'static str = CharterSchemaV1::NAME;
    const VERSION: u32 = CharterSchemaV1::VERSION;
    const CONTEXT: &'static str = "G0 promotion charter schema other context";
    const FIELDS: &'static [FieldSpec] = CHARTER_FIELDS_V1;
}

struct CharterFieldsVariant;
impl CanonicalSchema for CharterFieldsVariant {
    const DOMAIN: &'static str = CHARTER_SHARED_DOMAIN;
    const NAME: &'static str = CharterSchemaV1::NAME;
    const VERSION: u32 = CharterSchemaV1::VERSION;
    const CONTEXT: &'static str = CharterSchemaV1::CONTEXT;
    const FIELDS: &'static [FieldSpec] = CHARTER_FIELDS_VARIANT;
}

struct CharterChildBindingVariant;
impl CanonicalSchema for CharterChildBindingVariant {
    const DOMAIN: &'static str = CHARTER_SHARED_DOMAIN;
    const NAME: &'static str = CharterSchemaV1::NAME;
    const VERSION: u32 = CharterSchemaV1::VERSION;
    const CONTEXT: &'static str = CharterSchemaV1::CONTEXT;
    const FIELDS: &'static [FieldSpec] = CHARTER_CHILD_BINDING_VARIANT;
}

struct CharterOverDepthTailA;
impl CanonicalSchema for CharterOverDepthTailA {
    const DOMAIN: &'static str = "org.frankensim.test.promotion.charter-over-depth-tail-a.v1";
    const NAME: &'static str = "promotion-charter-over-depth-tail-a";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "G0 promotion charter divergent tail A";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required("tail-a", WireType::U64)];
}

struct CharterOverDepthTailB;
impl CanonicalSchema for CharterOverDepthTailB {
    const DOMAIN: &'static str = "org.frankensim.test.promotion.charter-over-depth-tail-b.v1";
    const NAME: &'static str = "promotion-charter-over-depth-tail-b";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "G0 promotion charter divergent tail B";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required("tail-b", WireType::Bytes)];
}

static CHARTER_OVER_DEPTH_TAIL_A: ChildSpec =
    ChildSpec::for_identity::<SemanticId<CharterOverDepthTailA>>();
static CHARTER_OVER_DEPTH_TAIL_B: ChildSpec =
    ChildSpec::for_identity::<SemanticId<CharterOverDepthTailB>>();

macro_rules! charter_over_depth_schema {
    ($schema:ident, $binding:ident, $child:ident, $level:literal) => {
        struct $schema;
        impl CanonicalSchema for $schema {
            const DOMAIN: &'static str = concat!(
                "org.frankensim.test.promotion.charter-over-depth-",
                $level,
                ".v1"
            );
            const NAME: &'static str = concat!("promotion-charter-over-depth-", $level);
            const VERSION: u32 = 1;
            const CONTEXT: &'static str = "G0 promotion charter finite over-depth chain";
            const FIELDS: &'static [FieldSpec] = &[FieldSpec::child_of("child", &$child)];
        }
        static $binding: ChildSpec = ChildSpec::for_identity::<SemanticId<$schema>>();
    };
}

charter_over_depth_schema!(
    CharterOverDepth16A,
    CHARTER_OVER_DEPTH_16_A,
    CHARTER_OVER_DEPTH_TAIL_A,
    "16"
);
charter_over_depth_schema!(
    CharterOverDepth16B,
    CHARTER_OVER_DEPTH_16_B,
    CHARTER_OVER_DEPTH_TAIL_B,
    "16"
);
charter_over_depth_schema!(
    CharterOverDepth15A,
    CHARTER_OVER_DEPTH_15_A,
    CHARTER_OVER_DEPTH_16_A,
    "15"
);
charter_over_depth_schema!(
    CharterOverDepth15B,
    CHARTER_OVER_DEPTH_15_B,
    CHARTER_OVER_DEPTH_16_B,
    "15"
);
charter_over_depth_schema!(
    CharterOverDepth14A,
    CHARTER_OVER_DEPTH_14_A,
    CHARTER_OVER_DEPTH_15_A,
    "14"
);
charter_over_depth_schema!(
    CharterOverDepth14B,
    CHARTER_OVER_DEPTH_14_B,
    CHARTER_OVER_DEPTH_15_B,
    "14"
);
charter_over_depth_schema!(
    CharterOverDepth13A,
    CHARTER_OVER_DEPTH_13_A,
    CHARTER_OVER_DEPTH_14_A,
    "13"
);
charter_over_depth_schema!(
    CharterOverDepth13B,
    CHARTER_OVER_DEPTH_13_B,
    CHARTER_OVER_DEPTH_14_B,
    "13"
);
charter_over_depth_schema!(
    CharterOverDepth12A,
    CHARTER_OVER_DEPTH_12_A,
    CHARTER_OVER_DEPTH_13_A,
    "12"
);
charter_over_depth_schema!(
    CharterOverDepth12B,
    CHARTER_OVER_DEPTH_12_B,
    CHARTER_OVER_DEPTH_13_B,
    "12"
);
charter_over_depth_schema!(
    CharterOverDepth11A,
    CHARTER_OVER_DEPTH_11_A,
    CHARTER_OVER_DEPTH_12_A,
    "11"
);
charter_over_depth_schema!(
    CharterOverDepth11B,
    CHARTER_OVER_DEPTH_11_B,
    CHARTER_OVER_DEPTH_12_B,
    "11"
);
charter_over_depth_schema!(
    CharterOverDepth10A,
    CHARTER_OVER_DEPTH_10_A,
    CHARTER_OVER_DEPTH_11_A,
    "10"
);
charter_over_depth_schema!(
    CharterOverDepth10B,
    CHARTER_OVER_DEPTH_10_B,
    CHARTER_OVER_DEPTH_11_B,
    "10"
);
charter_over_depth_schema!(
    CharterOverDepth09A,
    CHARTER_OVER_DEPTH_09_A,
    CHARTER_OVER_DEPTH_10_A,
    "09"
);
charter_over_depth_schema!(
    CharterOverDepth09B,
    CHARTER_OVER_DEPTH_09_B,
    CHARTER_OVER_DEPTH_10_B,
    "09"
);
charter_over_depth_schema!(
    CharterOverDepth08A,
    CHARTER_OVER_DEPTH_08_A,
    CHARTER_OVER_DEPTH_09_A,
    "08"
);
charter_over_depth_schema!(
    CharterOverDepth08B,
    CHARTER_OVER_DEPTH_08_B,
    CHARTER_OVER_DEPTH_09_B,
    "08"
);
charter_over_depth_schema!(
    CharterOverDepth07A,
    CHARTER_OVER_DEPTH_07_A,
    CHARTER_OVER_DEPTH_08_A,
    "07"
);
charter_over_depth_schema!(
    CharterOverDepth07B,
    CHARTER_OVER_DEPTH_07_B,
    CHARTER_OVER_DEPTH_08_B,
    "07"
);
charter_over_depth_schema!(
    CharterOverDepth06A,
    CHARTER_OVER_DEPTH_06_A,
    CHARTER_OVER_DEPTH_07_A,
    "06"
);
charter_over_depth_schema!(
    CharterOverDepth06B,
    CHARTER_OVER_DEPTH_06_B,
    CHARTER_OVER_DEPTH_07_B,
    "06"
);
charter_over_depth_schema!(
    CharterOverDepth05A,
    CHARTER_OVER_DEPTH_05_A,
    CHARTER_OVER_DEPTH_06_A,
    "05"
);
charter_over_depth_schema!(
    CharterOverDepth05B,
    CHARTER_OVER_DEPTH_05_B,
    CHARTER_OVER_DEPTH_06_B,
    "05"
);
charter_over_depth_schema!(
    CharterOverDepth04A,
    CHARTER_OVER_DEPTH_04_A,
    CHARTER_OVER_DEPTH_05_A,
    "04"
);
charter_over_depth_schema!(
    CharterOverDepth04B,
    CHARTER_OVER_DEPTH_04_B,
    CHARTER_OVER_DEPTH_05_B,
    "04"
);
charter_over_depth_schema!(
    CharterOverDepth03A,
    CHARTER_OVER_DEPTH_03_A,
    CHARTER_OVER_DEPTH_04_A,
    "03"
);
charter_over_depth_schema!(
    CharterOverDepth03B,
    CHARTER_OVER_DEPTH_03_B,
    CHARTER_OVER_DEPTH_04_B,
    "03"
);
charter_over_depth_schema!(
    CharterOverDepth02A,
    CHARTER_OVER_DEPTH_02_A,
    CHARTER_OVER_DEPTH_03_A,
    "02"
);
charter_over_depth_schema!(
    CharterOverDepth02B,
    CHARTER_OVER_DEPTH_02_B,
    CHARTER_OVER_DEPTH_03_B,
    "02"
);
charter_over_depth_schema!(
    CharterOverDepth01A,
    CHARTER_OVER_DEPTH_01_A,
    CHARTER_OVER_DEPTH_02_A,
    "01"
);
charter_over_depth_schema!(
    CharterOverDepth01B,
    CHARTER_OVER_DEPTH_01_B,
    CHARTER_OVER_DEPTH_02_B,
    "01"
);

struct CharterOverDepthRootA;
impl CanonicalSchema for CharterOverDepthRootA {
    const DOMAIN: &'static str = "org.frankensim.test.promotion.charter-over-depth-root.v1";
    const NAME: &'static str = "promotion-charter-over-depth-root";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "G0 promotion charter over-depth root";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::child_of("child", &CHARTER_OVER_DEPTH_01_A)];
}

struct CharterOverDepthRootB;
impl CanonicalSchema for CharterOverDepthRootB {
    const DOMAIN: &'static str = CharterOverDepthRootA::DOMAIN;
    const NAME: &'static str = CharterOverDepthRootA::NAME;
    const VERSION: u32 = CharterOverDepthRootA::VERSION;
    const CONTEXT: &'static str = CharterOverDepthRootA::CONTEXT;
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::child_of("child", &CHARTER_OVER_DEPTH_01_B)];
}

type Subject = SemanticId<SubjectV1>;
type PresentedAuthority = AuthorityRef<Subject, VerifierSchemaV1, PolicySchemaV1, Presented>;
type VerifiedAuthority = AuthorityRef<Subject, VerifierSchemaV1, PolicySchemaV1, Verified>;
type Witness = PromotionWitness<Subject, VerifierSchemaV1, PolicySchemaV1>;
type Root = PromotionTrustRoot<VerifierSchemaV1, PolicySchemaV1>;

const FIXED_CHARTER_CONTEXT: &str = "promotion-charter-fixed-context";

fn push_charter_reference_field(preimage: &mut Vec<u8>, bytes: &[u8]) {
    let length = u64::try_from(bytes.len()).expect("test fields fit u64 framing");
    preimage.extend_from_slice(&length.to_le_bytes());
    preimage.extend_from_slice(bytes);
}

fn push_capability_descriptor_reference(
    preimage: &mut Vec<u8>,
    descriptor: Option<PromotionCapabilityDescriptor>,
) {
    match descriptor {
        None => push_charter_reference_field(preimage, &[0]),
        Some(descriptor) => {
            push_charter_reference_field(preimage, &[1]);
            push_charter_reference_field(
                preimage,
                descriptor.implementation().content_id().as_bytes(),
            );
            push_charter_reference_field(
                preimage,
                &descriptor.implementation().length().to_le_bytes(),
            );
            push_charter_reference_field(
                preimage,
                descriptor.configuration().content_id().as_bytes(),
            );
            push_charter_reference_field(
                preimage,
                &descriptor.configuration().length().to_le_bytes(),
            );
            push_charter_reference_field(preimage, &descriptor.protocol_version().to_le_bytes());
        }
    }
}

fn push_identity_reference_field(preimage: &mut Vec<u8>, bytes: &[u8], framed: bool) {
    if framed {
        push_charter_reference_field(preimage, bytes);
    } else {
        preimage.extend_from_slice(bytes);
    }
}

fn changed_digest_bytes(mut bytes: [u8; 32]) -> [u8; 32] {
    bytes[0] ^= 0xA5;
    bytes
}

#[derive(Clone, Copy)]
struct RequestIdentityReference<'a> {
    domain: &'a str,
    version: u32,
    subject_role: u8,
    subject_schema: [u8; 32],
    subject_id: [u8; 32],
    subject_preimage: ContentId,
    subject_bytes: u64,
    anchor: ContentId,
    verifier_role: u8,
    verifier_schema: [u8; 32],
    verifier_id: [u8; 32],
    verifier_observation: ByteObservation,
    key_policy_role: u8,
    key_policy_schema: [u8; 32],
    key_policy_id: [u8; 32],
    key_policy_observation: ByteObservation,
    root_charter: [u8; 32],
    root_context: &'a str,
    attempt: ContentId,
    decision_context: ContentId,
    epoch: u64,
    sequence: u64,
    predecessor: Option<[u8; 32]>,
}

fn request_identity_reference(input: &RequestIdentityReference<'_>, framed: bool) -> ContentHash {
    let mut preimage = Vec::new();
    push_identity_reference_field(&mut preimage, &input.version.to_le_bytes(), framed);
    push_identity_reference_field(&mut preimage, &[input.subject_role], framed);
    push_identity_reference_field(&mut preimage, &input.subject_schema, framed);
    push_identity_reference_field(&mut preimage, &input.subject_id, framed);
    push_identity_reference_field(&mut preimage, input.subject_preimage.as_bytes(), framed);
    push_identity_reference_field(&mut preimage, &input.subject_bytes.to_le_bytes(), framed);
    push_identity_reference_field(&mut preimage, input.anchor.as_bytes(), framed);
    push_identity_reference_field(&mut preimage, &[input.verifier_role], framed);
    push_identity_reference_field(&mut preimage, &input.verifier_schema, framed);
    push_identity_reference_field(&mut preimage, &input.verifier_id, framed);
    push_identity_reference_field(
        &mut preimage,
        input.verifier_observation.content_id().as_bytes(),
        framed,
    );
    push_identity_reference_field(
        &mut preimage,
        &input.verifier_observation.length().to_le_bytes(),
        framed,
    );
    push_identity_reference_field(&mut preimage, &[input.key_policy_role], framed);
    push_identity_reference_field(&mut preimage, &input.key_policy_schema, framed);
    push_identity_reference_field(&mut preimage, &input.key_policy_id, framed);
    push_identity_reference_field(
        &mut preimage,
        input.key_policy_observation.content_id().as_bytes(),
        framed,
    );
    push_identity_reference_field(
        &mut preimage,
        &input.key_policy_observation.length().to_le_bytes(),
        framed,
    );
    push_identity_reference_field(&mut preimage, &input.root_charter, framed);
    push_identity_reference_field(&mut preimage, input.root_context.as_bytes(), framed);
    push_identity_reference_field(&mut preimage, input.attempt.as_bytes(), framed);
    push_identity_reference_field(&mut preimage, input.decision_context.as_bytes(), framed);
    push_identity_reference_field(&mut preimage, &input.epoch.to_le_bytes(), framed);
    push_identity_reference_field(&mut preimage, &input.sequence.to_le_bytes(), framed);
    match input.predecessor {
        None => push_identity_reference_field(&mut preimage, &[0], framed),
        Some(predecessor) => {
            push_identity_reference_field(&mut preimage, &[1], framed);
            push_identity_reference_field(&mut preimage, &predecessor, framed);
        }
    }
    hash_domain(input.domain, &preimage)
}

fn push_decision_descriptor_reference(
    preimage: &mut Vec<u8>,
    descriptor: PromotionCapabilityDescriptor,
    framed: bool,
) {
    push_identity_reference_field(
        preimage,
        descriptor.implementation().content_id().as_bytes(),
        framed,
    );
    push_identity_reference_field(
        preimage,
        &descriptor.implementation().length().to_le_bytes(),
        framed,
    );
    push_identity_reference_field(
        preimage,
        descriptor.configuration().content_id().as_bytes(),
        framed,
    );
    push_identity_reference_field(
        preimage,
        &descriptor.configuration().length().to_le_bytes(),
        framed,
    );
    push_identity_reference_field(
        preimage,
        &descriptor.protocol_version().to_le_bytes(),
        framed,
    );
}

fn descriptor_identity_variants(
    descriptor: PromotionCapabilityDescriptor,
) -> [PromotionCapabilityDescriptor; 5] {
    [
        PromotionCapabilityDescriptor::new(
            ByteObservation::new(
                ContentId::of_bytes(b"other capability implementation root"),
                descriptor.implementation().length(),
            ),
            descriptor.configuration(),
            descriptor.protocol_version(),
        ),
        PromotionCapabilityDescriptor::new(
            ByteObservation::new(
                descriptor.implementation().content_id(),
                descriptor.implementation().length() + 1,
            ),
            descriptor.configuration(),
            descriptor.protocol_version(),
        ),
        PromotionCapabilityDescriptor::new(
            descriptor.implementation(),
            ByteObservation::new(
                ContentId::of_bytes(b"other capability configuration root"),
                descriptor.configuration().length(),
            ),
            descriptor.protocol_version(),
        ),
        PromotionCapabilityDescriptor::new(
            descriptor.implementation(),
            ByteObservation::new(
                descriptor.configuration().content_id(),
                descriptor.configuration().length() + 1,
            ),
            descriptor.protocol_version(),
        ),
        PromotionCapabilityDescriptor::new(
            descriptor.implementation(),
            descriptor.configuration(),
            descriptor.protocol_version() + 1,
        ),
    ]
}

fn verification_decision_reference(
    domain: &str,
    version: u32,
    stage_tag: u8,
    request: [u8; 32],
    descriptor: PromotionCapabilityDescriptor,
    statement: ContentId,
    framed: bool,
) -> ContentHash {
    let mut preimage = Vec::new();
    push_identity_reference_field(&mut preimage, &version.to_le_bytes(), framed);
    push_identity_reference_field(&mut preimage, &[stage_tag], framed);
    push_identity_reference_field(&mut preimage, &request, framed);
    push_decision_descriptor_reference(&mut preimage, descriptor, framed);
    push_identity_reference_field(&mut preimage, statement.as_bytes(), framed);
    hash_domain(domain, &preimage)
}

#[allow(clippy::too_many_arguments)]
fn admission_decision_reference(
    domain: &str,
    version: u32,
    stage_tag: u8,
    request: [u8; 32],
    verification: [u8; 32],
    descriptor: PromotionCapabilityDescriptor,
    statement: ContentId,
    framed: bool,
) -> ContentHash {
    let mut preimage = Vec::new();
    push_identity_reference_field(&mut preimage, &version.to_le_bytes(), framed);
    push_identity_reference_field(&mut preimage, &[stage_tag], framed);
    push_identity_reference_field(&mut preimage, &request, framed);
    push_identity_reference_field(&mut preimage, &verification, framed);
    push_decision_descriptor_reference(&mut preimage, descriptor, framed);
    push_identity_reference_field(&mut preimage, statement.as_bytes(), framed);
    hash_domain(domain, &preimage)
}

fn current_charter_reference<V, P>(
    verifier: ObservedIdentity<VerifierId<V>>,
    key_policy: ObservedIdentity<KeyPolicyId<P>>,
    context: &str,
) -> ContentHash
where
    V: CanonicalSchema,
    P: CanonicalSchema,
{
    current_charter_reference_with_roles(
        verifier,
        key_policy,
        context,
        <VerifierId<V> as StrongIdentity>::ROLE,
        <KeyPolicyId<P> as StrongIdentity>::ROLE,
    )
}

fn current_charter_reference_with_roles<V, P>(
    verifier: ObservedIdentity<VerifierId<V>>,
    key_policy: ObservedIdentity<KeyPolicyId<P>>,
    context: &str,
    verifier_identity_role: IdentityRole,
    key_policy_identity_role: IdentityRole,
) -> ContentHash
where
    V: CanonicalSchema,
    P: CanonicalSchema,
{
    charter_reference_with_roles_and_capabilities(
        verifier,
        key_policy,
        context,
        verifier_identity_role,
        key_policy_identity_role,
        0,
        0,
        None,
        None,
    )
}

#[allow(clippy::too_many_arguments)]
fn charter_reference_with_roles_and_capabilities<V, P>(
    verifier: ObservedIdentity<VerifierId<V>>,
    key_policy: ObservedIdentity<KeyPolicyId<P>>,
    context: &str,
    verifier_identity_role: IdentityRole,
    key_policy_identity_role: IdentityRole,
    capability_mode: u8,
    decision_epoch: u64,
    verifier_capability: Option<PromotionCapabilityDescriptor>,
    admission_capability: Option<PromotionCapabilityDescriptor>,
) -> ContentHash
where
    V: CanonicalSchema,
    P: CanonicalSchema,
{
    let verifier_schema = SchemaId::<V>::for_schema();
    let key_policy_schema = SchemaId::<P>::for_schema();
    let verifier_id = verifier.id();
    let verifier_observation = verifier.bytes();
    let verifier_observation_root = verifier_observation.content_id();
    let key_policy_id = key_policy.id();
    let key_policy_observation = key_policy.bytes();
    let key_policy_observation_root = key_policy_observation.content_id();
    let identity_version = EXPECTED_CURRENT_CHARTER_VERSION.to_le_bytes();
    let verifier_role = [verifier_identity_role.tag()];
    let verifier_length = verifier_observation.length().to_le_bytes();
    let key_policy_role = [key_policy_identity_role.tag()];
    let key_policy_length = key_policy_observation.length().to_le_bytes();
    let capability_mode = [capability_mode];
    let decision_epoch = decision_epoch.to_le_bytes();
    let mut preimage = Vec::new();
    for field in [
        identity_version.as_slice(),
        verifier_role.as_slice(),
        V::DOMAIN.as_bytes(),
        verifier_schema.as_bytes(),
        verifier_id.as_bytes(),
        verifier_observation_root.as_bytes(),
        verifier_length.as_slice(),
        key_policy_role.as_slice(),
        P::DOMAIN.as_bytes(),
        key_policy_schema.as_bytes(),
        key_policy_id.as_bytes(),
        key_policy_observation_root.as_bytes(),
        key_policy_length.as_slice(),
        capability_mode.as_slice(),
        decision_epoch.as_slice(),
    ] {
        push_charter_reference_field(&mut preimage, field);
    }
    push_capability_descriptor_reference(&mut preimage, verifier_capability);
    push_capability_descriptor_reference(&mut preimage, admission_capability);
    push_charter_reference_field(&mut preimage, context.as_bytes());
    hash_domain(EXPECTED_CURRENT_CHARTER_DOMAIN, &preimage)
}

fn owner_charter_reference<V, P>(
    verifier: ObservedIdentity<VerifierId<V>>,
    key_policy: ObservedIdentity<KeyPolicyId<P>>,
    context: &str,
    decision_epoch: u64,
    verifier_capability: PromotionCapabilityDescriptor,
    admission_capability: PromotionCapabilityDescriptor,
) -> ContentHash
where
    V: CanonicalSchema,
    P: CanonicalSchema,
{
    charter_reference_with_roles_and_capabilities(
        verifier,
        key_policy,
        context,
        <VerifierId<V> as StrongIdentity>::ROLE,
        <KeyPolicyId<P> as StrongIdentity>::ROLE,
        1,
        decision_epoch,
        Some(verifier_capability),
        Some(admission_capability),
    )
}

fn legacy_charter_reference<V, P>(
    verifier: ObservedIdentity<VerifierId<V>>,
    key_policy: ObservedIdentity<KeyPolicyId<P>>,
    context: &str,
) -> ContentHash
where
    V: CanonicalSchema,
    P: CanonicalSchema,
{
    let verifier_id = verifier.id();
    let verifier_observation = verifier.bytes();
    let verifier_observation_root = verifier_observation.content_id();
    let key_policy_id = key_policy.id();
    let key_policy_observation = key_policy.bytes();
    let key_policy_observation_root = key_policy_observation.content_id();
    let verifier_length = verifier_observation.length().to_le_bytes();
    let key_policy_length = key_policy_observation.length().to_le_bytes();
    let mut preimage = Vec::new();
    for field in [
        V::DOMAIN.as_bytes(),
        P::DOMAIN.as_bytes(),
        verifier_id.as_bytes(),
        verifier_observation_root.as_bytes(),
        verifier_length.as_slice(),
        key_policy_id.as_bytes(),
        key_policy_observation_root.as_bytes(),
        key_policy_length.as_slice(),
        context.as_bytes(),
    ] {
        push_charter_reference_field(&mut preimage, field);
    }
    hash_domain(EXPECTED_LEGACY_CHARTER_V1_DOMAIN, &preimage)
}

fn legacy_v2_charter_reference<V, P>(
    verifier: ObservedIdentity<VerifierId<V>>,
    key_policy: ObservedIdentity<KeyPolicyId<P>>,
    context: &str,
) -> ContentHash
where
    V: CanonicalSchema,
    P: CanonicalSchema,
{
    let verifier_schema = SchemaId::<V>::for_schema();
    let key_policy_schema = SchemaId::<P>::for_schema();
    let verifier_id = verifier.id();
    let verifier_observation = verifier.bytes();
    let key_policy_id = key_policy.id();
    let key_policy_observation = key_policy.bytes();
    let identity_version = 2u32.to_le_bytes();
    let verifier_role = [<VerifierId<V> as StrongIdentity>::ROLE.tag()];
    let key_policy_role = [<KeyPolicyId<P> as StrongIdentity>::ROLE.tag()];
    let verifier_length = verifier_observation.length().to_le_bytes();
    let key_policy_length = key_policy_observation.length().to_le_bytes();
    let mut preimage = Vec::new();
    for field in [
        identity_version.as_slice(),
        verifier_role.as_slice(),
        V::DOMAIN.as_bytes(),
        verifier_schema.as_bytes(),
        verifier_id.as_bytes(),
        verifier_observation.content_id().as_bytes(),
        verifier_length.as_slice(),
        key_policy_role.as_slice(),
        P::DOMAIN.as_bytes(),
        key_policy_schema.as_bytes(),
        key_policy_id.as_bytes(),
        key_policy_observation.content_id().as_bytes(),
        key_policy_length.as_slice(),
        context.as_bytes(),
    ] {
        push_charter_reference_field(&mut preimage, field);
    }
    hash_domain(EXPECTED_LEGACY_CHARTER_V2_DOMAIN, &preimage)
}

fn fixed_charter_bindings<V, P>() -> (
    ObservedIdentity<VerifierId<V>>,
    ObservedIdentity<KeyPolicyId<P>>,
)
where
    V: CanonicalSchema,
    P: CanonicalSchema,
{
    let verifier_digest = hash_bytes(b"fixed promotion charter verifier id");
    let policy_digest = hash_bytes(b"fixed promotion charter key-policy id");
    let verifier = VerifierId::<V>::parse_slice(verifier_digest.as_bytes())
        .expect("fixed verifier digest parses under every exact schema type");
    let key_policy = KeyPolicyId::<P>::parse_slice(policy_digest.as_bytes())
        .expect("fixed key-policy digest parses under every exact schema type");
    (
        ObservedIdentity::presented(
            verifier,
            ByteObservation::new(
                ContentId::of_bytes(b"fixed promotion charter verifier bytes"),
                41,
            ),
        ),
        ObservedIdentity::presented(
            key_policy,
            ByteObservation::new(
                ContentId::of_bytes(b"fixed promotion charter key-policy bytes"),
                43,
            ),
        ),
    )
}

fn fixed_charter_root<V, P>() -> PromotionTrustRoot<V, P>
where
    V: CanonicalSchema,
    P: CanonicalSchema,
{
    let (verifier, key_policy) = fixed_charter_bindings::<V, P>();
    PromotionTrustRoot::configure(verifier, key_policy, FIXED_CHARTER_CONTEXT)
        .expect("fixed charter root configures")
}

fn subject_receipt(value: u64) -> IdentityReceipt<Subject> {
    CanonicalEncoder::<Subject, _>::new(LIMITS, NeverCancel)
        .expect("valid subject schema")
        .u64(Field::new(0, "value"), value)
        .expect("subject field")
        .finish()
        .expect("subject receipt")
}

fn verifier_receipt(key: u64) -> IdentityReceipt<VerifierId<VerifierSchemaV1>> {
    CanonicalEncoder::<VerifierId<VerifierSchemaV1>, _>::new(LIMITS, NeverCancel)
        .expect("valid verifier schema")
        .u64(Field::new(0, "key"), key)
        .expect("verifier field")
        .finish()
        .expect("verifier receipt")
}

fn policy_receipt(rule: u64) -> IdentityReceipt<KeyPolicyId<PolicySchemaV1>> {
    CanonicalEncoder::<KeyPolicyId<PolicySchemaV1>, _>::new(LIMITS, NeverCancel)
        .expect("valid policy schema")
        .u64(Field::new(0, "rule"), rule)
        .expect("policy field")
        .finish()
        .expect("policy receipt")
}

fn anchor() -> ExternalAnchorRef {
    ExternalAnchorRef::presented(ContentId::of_bytes(b"promotion-test-anchor"))
}

/// The adversary: accepts everything it is shown.
struct PermitAll;

impl AuthorityVerifier<Subject, VerifierSchemaV1, PolicySchemaV1> for PermitAll {
    type Error = core::convert::Infallible;
    fn verify(&self, _presented: &PresentedAuthority) -> Result<(), Self::Error> {
        Ok(())
    }
}

impl AuthorityAdmitter<Subject, VerifierSchemaV1, PolicySchemaV1> for PermitAll {
    type Error = core::convert::Infallible;
    fn admit(&self, _verified: &VerifiedAuthority) -> Result<(), Self::Error> {
        Ok(())
    }
}

fn permit_all_admitted(
    subject: IdentityReceipt<Subject>,
    verifier: IdentityReceipt<VerifierId<VerifierSchemaV1>>,
    policy: IdentityReceipt<KeyPolicyId<PolicySchemaV1>>,
) -> AuthorityRef<Subject, VerifierSchemaV1, PolicySchemaV1, fs_blake3::identity::Admitted> {
    permit_all_admitted_with_anchor(subject, verifier, policy, anchor())
}

fn permit_all_admitted_with_anchor(
    subject: IdentityReceipt<Subject>,
    verifier: IdentityReceipt<VerifierId<VerifierSchemaV1>>,
    policy: IdentityReceipt<KeyPolicyId<PolicySchemaV1>>,
    presented_anchor: ExternalAnchorRef,
) -> AuthorityRef<Subject, VerifierSchemaV1, PolicySchemaV1, fs_blake3::identity::Admitted> {
    AuthorityRef::present(subject, presented_anchor, verifier.id(), policy.id())
        .verify(&PermitAll)
        .expect("permit-all verifies anything")
        .admit(&PermitAll)
        .expect("permit-all admits anything")
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum OwnerBehavior {
    Exact,
    Refuse,
    Cancel,
    Panic,
}

#[derive(Debug, Clone, Copy)]
struct ExactOwnerVerifier {
    descriptor: PromotionCapabilityDescriptor,
    expected_subject: ContentId,
    expected_anchor: ContentId,
    behavior: OwnerBehavior,
}

impl OwnerPromotionVerifier for ExactOwnerVerifier {
    fn descriptor(&self) -> PromotionCapabilityDescriptor {
        self.descriptor
    }

    fn verify(&self, request: &PromotionDecisionRequest) -> PromotionCapabilityVerdict {
        match self.behavior {
            OwnerBehavior::Panic => panic!("injected owner verifier fault"),
            OwnerBehavior::Cancel => PromotionCapabilityVerdict::Cancelled {
                reason: ContentId::of_bytes(b"owner verifier cancelled"),
            },
            OwnerBehavior::Refuse => PromotionCapabilityVerdict::Refuse {
                reason: ContentId::of_bytes(b"owner verifier forced refusal"),
            },
            OwnerBehavior::Exact => {
                if request.subject_role() != IdentityRole::Semantic
                    || request.subject_schema() != *SchemaId::<SubjectV1>::for_schema().as_bytes()
                    || request.subject_preimage() != self.expected_subject
                    || request.anchor() != self.expected_anchor
                    || request.verifier_role() != IdentityRole::Verifier
                    || request.verifier_schema()
                        != *SchemaId::<VerifierSchemaV1>::for_schema().as_bytes()
                    || request.key_policy_role() != IdentityRole::KeyPolicy
                    || request.key_policy_schema()
                        != *SchemaId::<PolicySchemaV1>::for_schema().as_bytes()
                    || request.root_context() != "promotion-test"
                    || request.scope().decision_context()
                        != ContentId::of_bytes(b"promotion-test/decision-context/v1")
                {
                    PromotionCapabilityVerdict::Refuse {
                        reason: ContentId::of_bytes(b"owner verifier exact binding mismatch"),
                    }
                } else {
                    PromotionCapabilityVerdict::Approve {
                        statement: ContentId::of_bytes(b"owner verifier approved exact request"),
                    }
                }
            }
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct ExactOwnerAdmitter {
    descriptor: PromotionCapabilityDescriptor,
    behavior: OwnerBehavior,
}

impl OwnerPromotionAdmitter for ExactOwnerAdmitter {
    fn descriptor(&self) -> PromotionCapabilityDescriptor {
        self.descriptor
    }

    fn admit(&self, request: &PromotionAdmissionRequest) -> PromotionCapabilityVerdict {
        let verification_bound_statement =
            ContentId::of_bytes(request.verification_decision().as_bytes());
        match self.behavior {
            OwnerBehavior::Panic => panic!("injected owner admission fault"),
            OwnerBehavior::Cancel => PromotionCapabilityVerdict::Cancelled {
                reason: ContentId::of_bytes(b"owner admission cancelled"),
            },
            OwnerBehavior::Refuse => PromotionCapabilityVerdict::Refuse {
                reason: ContentId::of_bytes(b"owner admission forced refusal"),
            },
            OwnerBehavior::Exact => {
                if request.verification_capability() != owner_verifier_descriptor()
                    || request.verification_statement()
                        != ContentId::of_bytes(b"owner verifier approved exact request")
                {
                    PromotionCapabilityVerdict::Refuse {
                        reason: ContentId::of_bytes(
                            b"owner admission verification evidence mismatch",
                        ),
                    }
                } else {
                    PromotionCapabilityVerdict::Approve {
                        statement: verification_bound_statement,
                    }
                }
            }
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct ForeignPermitAllVerifier {
    descriptor: PromotionCapabilityDescriptor,
}

impl OwnerPromotionVerifier for ForeignPermitAllVerifier {
    fn descriptor(&self) -> PromotionCapabilityDescriptor {
        self.descriptor
    }

    fn verify(&self, _request: &PromotionDecisionRequest) -> PromotionCapabilityVerdict {
        PromotionCapabilityVerdict::Approve {
            statement: ContentId::of_bytes(b"foreign verifier permit-all statement"),
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct ForeignPermitAllAdmitter {
    descriptor: PromotionCapabilityDescriptor,
}

impl OwnerPromotionAdmitter for ForeignPermitAllAdmitter {
    fn descriptor(&self) -> PromotionCapabilityDescriptor {
        self.descriptor
    }

    fn admit(&self, _request: &PromotionAdmissionRequest) -> PromotionCapabilityVerdict {
        PromotionCapabilityVerdict::Approve {
            statement: ContentId::of_bytes(b"foreign admission permit-all statement"),
        }
    }
}

type OwnerRoot = PromotionTrustRoot<
    VerifierSchemaV1,
    PolicySchemaV1,
    OwnerPromotionCapabilities<ExactOwnerVerifier, ExactOwnerAdmitter>,
>;

fn capability_descriptor(
    implementation: &[u8],
    configuration: &[u8],
    version: u32,
) -> PromotionCapabilityDescriptor {
    PromotionCapabilityDescriptor::new(
        ByteObservation::new(
            ContentId::of_bytes(implementation),
            implementation.len() as u64,
        ),
        ByteObservation::new(
            ContentId::of_bytes(configuration),
            configuration.len() as u64,
        ),
        version,
    )
}

fn owner_verifier_descriptor() -> PromotionCapabilityDescriptor {
    capability_descriptor(
        b"owner promotion verifier executable v1",
        b"owner promotion verifier immutable config v1",
        1,
    )
}

fn owner_admitter_descriptor() -> PromotionCapabilityDescriptor {
    capability_descriptor(
        b"owner promotion admission executable v1",
        b"owner promotion admission immutable config v1",
        1,
    )
}

fn decision_scope(epoch: u64, sequence: u64) -> PromotionDecisionScope {
    PromotionDecisionScope::fresh(
        ContentId::of_bytes(&sequence.to_le_bytes()),
        ContentId::of_bytes(b"promotion-test/decision-context/v1"),
        epoch,
        sequence,
    )
}

fn owner_root(
    subject: IdentityReceipt<Subject>,
    verifier: IdentityReceipt<VerifierId<VerifierSchemaV1>>,
    policy: IdentityReceipt<KeyPolicyId<PolicySchemaV1>>,
    epoch: u64,
    verifier_behavior: OwnerBehavior,
    admission_behavior: OwnerBehavior,
) -> OwnerRoot {
    Root::configure_owner_executed(
        ObservedIdentity::from_receipt(verifier),
        ObservedIdentity::from_receipt(policy),
        "promotion-test",
        epoch,
        ExactOwnerVerifier {
            descriptor: owner_verifier_descriptor(),
            expected_subject: subject.canonical_preimage(),
            expected_anchor: anchor().content_id(),
            behavior: verifier_behavior,
        },
        ExactOwnerAdmitter {
            descriptor: owner_admitter_descriptor(),
            behavior: admission_behavior,
        },
    )
    .expect("owner root configures")
}

#[test]
fn permit_all_reaches_policy_relative_admission_but_never_promotion() {
    // The adversary presents ITS OWN verifier/policy identities and
    // sails through the generic ladder — that is the documented
    // policy-relative lane.
    let rogue_verifier = verifier_receipt(0xBAD);
    let rogue_policy = policy_receipt(0xBAD);
    let admitted = permit_all_admitted(subject_receipt(7), rogue_verifier, rogue_policy);

    // The domain owner's root was configured independently for the REAL
    // verifier and policy; the rogue admission cannot cross it.
    let root = Root::configure(
        ObservedIdentity::from_receipt(verifier_receipt(0x600D)),
        ObservedIdentity::from_receipt(policy_receipt(0x600D)),
        "promotion-test",
    )
    .expect("root configures");
    let refusal = root
        .admit_for_promotion(
            &admitted,
            ObservedIdentity::from_receipt(rogue_verifier).bytes(),
            ObservedIdentity::from_receipt(rogue_policy).bytes(),
        )
        .expect_err("a foreign verifier must never mint promotion");
    assert_eq!(refusal, PromotionRefusal::ForeignVerifier);

    // The more dangerous bypass presents the OWNER's exact IDs and byte
    // observations. A configuration-only compatibility root still has no
    // executable owner capability and therefore cannot mint.
    let owner_verifier = verifier_receipt(0x600D);
    let owner_policy = policy_receipt(0x600D);
    let exact_id_spoof = permit_all_admitted(subject_receipt(8), owner_verifier, owner_policy);
    assert_eq!(
        root.admit_for_promotion(
            &exact_id_spoof,
            ObservedIdentity::from_receipt(owner_verifier).bytes(),
            ObservedIdentity::from_receipt(owner_policy).bytes(),
        )
        .expect_err("public owner IDs and observations are not executable authority"),
        PromotionRefusal::OwnerCapabilitiesUnavailable
    );
}

#[test]
fn the_owner_executed_root_mints_replay_evidence_and_binds_live_authority() {
    let subject = subject_receipt(11);
    let verifier = verifier_receipt(0x600D);
    let policy = policy_receipt(0x600D);
    let admitted = permit_all_admitted(subject, verifier, policy);
    let mut root = owner_root(
        subject,
        verifier,
        policy,
        7,
        OwnerBehavior::Exact,
        OwnerBehavior::Exact,
    );
    let witness: Witness = root
        .decide_for_promotion(
            &admitted,
            ObservedIdentity::from_receipt(verifier).bytes(),
            ObservedIdentity::from_receipt(policy).bytes(),
            decision_scope(7, 1),
        )
        .expect("both stored owner capabilities approve");
    let bound = root
        .bind_witness(&witness)
        .expect("raw evidence binds only to its exact live owner root");
    assert_eq!(bound.subject(), subject);
    assert_eq!(bound.decision_id(), witness.decision_id());
    // Exact subject receipt/preimage, anchor, verifier and policy
    // observations, and context remain bound.
    assert_eq!(witness.subject(), subject);
    assert_eq!(witness.anchor(), anchor());
    assert_eq!(
        witness.verifier().bytes(),
        ObservedIdentity::from_receipt(verifier).bytes()
    );
    assert_eq!(
        witness.key_policy().bytes(),
        ObservedIdentity::from_receipt(policy).bytes()
    );
    assert_eq!(witness.context(), "promotion-test");
    assert_eq!(witness.scope(), decision_scope(7, 1));
    assert_eq!(witness.root_charter(), root.charter());
    // Bounded audit retains canonical decision identities, not payloads.
    let audit = witness.audit();
    assert_eq!(bound.audit(), audit);
    assert_eq!(audit.subject_role, IdentityRole::Semantic);
    assert_eq!(
        audit.subject_schema,
        *SchemaId::<SubjectV1>::for_schema().as_bytes()
    );
    assert_eq!(audit.subject_id, *subject.id().as_bytes());
    assert_eq!(audit.subject_preimage, subject.canonical_preimage());
    assert_eq!(audit.subject_bytes, subject.canonical_bytes());
    assert_eq!(audit.anchor, anchor().content_id());
    assert_eq!(audit.verifier_domain, VerifierSchemaV1::DOMAIN);
    assert_eq!(audit.verifier_role, IdentityRole::Verifier);
    assert_eq!(audit.verifier_id, *verifier.id().as_bytes());
    assert_eq!(audit.key_policy_domain, PolicySchemaV1::DOMAIN);
    assert_eq!(audit.key_policy_role, IdentityRole::KeyPolicy);
    assert_eq!(audit.key_policy_id, *policy.id().as_bytes());
    assert_eq!(
        audit.verifier_observation,
        ObservedIdentity::from_receipt(verifier).bytes()
    );
    assert_eq!(audit.context, "promotion-test");
    assert_eq!(
        audit.verification_statement,
        witness.verification_statement()
    );
    assert_eq!(audit.admission_statement, witness.admission_statement());
    assert_eq!(audit.disposition, PromotionDecisionDisposition::Approved);
    assert_eq!(
        witness.disposition(),
        PromotionDecisionDisposition::Approved
    );
    assert_eq!(audit.request_id, witness.request_id());
    assert_eq!(audit.verification_decision, witness.verification_decision());
    assert_eq!(audit.admission_decision, witness.admission_decision());
    assert_eq!(audit.decision_id, witness.decision_id());
}

#[test]
fn exact_owner_ids_and_observations_cannot_impersonate_owner_execution() {
    let authorized_subject = subject_receipt(41);
    let spoofed_subject = subject_receipt(42);
    let verifier = verifier_receipt(0x600D);
    let policy = policy_receipt(0x600D);
    let spoofed = permit_all_admitted(spoofed_subject, verifier, policy);
    let mut owner = owner_root(
        authorized_subject,
        verifier,
        policy,
        9,
        OwnerBehavior::Exact,
        OwnerBehavior::Exact,
    );

    let refusal = owner
        .decide_for_promotion(
            &spoofed,
            ObservedIdentity::from_receipt(verifier).bytes(),
            ObservedIdentity::from_receipt(policy).bytes(),
            decision_scope(9, 1),
        )
        .expect_err("owner verifier must execute and reject the wrong subject");
    assert!(matches!(
        refusal,
        PromotionRefusal::CapabilityRefused {
            stage: PromotionCapabilityStage::Verification,
            ..
        }
    ));

    // Burn-before-call: fixing the subject does not allow the refused scope to
    // be replayed.
    let authorized = permit_all_admitted(authorized_subject, verifier, policy);
    assert_eq!(
        owner
            .decide_for_promotion(
                &authorized,
                ObservedIdentity::from_receipt(verifier).bytes(),
                ObservedIdentity::from_receipt(policy).bytes(),
                decision_scope(9, 1),
            )
            .expect_err("refused sequence was already burned"),
        PromotionRefusal::StaleOrReplayedDecision {
            last_attempted: 1,
            presented: 1,
        }
    );
}

#[test]
fn same_id_different_bytes_refuses_with_both_observations_retained() {
    let verifier = verifier_receipt(0x600D);
    let policy = policy_receipt(0x600D);
    let admitted = permit_all_admitted(subject_receipt(3), verifier, policy);
    let root = Root::configure(
        ObservedIdentity::from_receipt(verifier),
        ObservedIdentity::from_receipt(policy),
        "promotion-test",
    )
    .expect("root configures");
    // The adversary claims the trusted verifier ID over DIFFERENT
    // canonical bytes (a would-be second preimage). The root refuses
    // and retains both observations, neither privileged.
    let forged = ByteObservation::new(ContentId::of_bytes(b"not-the-verifier-bytes"), 999);
    let refusal = root
        .admit_for_promotion(
            &admitted,
            forged,
            ObservedIdentity::from_receipt(policy).bytes(),
        )
        .expect_err("same ID over different bytes must refuse");
    let PromotionRefusal::VerifierObservationMismatch {
        configured,
        presented,
    } = refusal
    else {
        panic!("expected a verifier observation mismatch, got {refusal:?}");
    };
    assert_eq!(configured, ObservedIdentity::from_receipt(verifier).bytes());
    assert_eq!(presented, forged);

    // Same discipline on the key-policy axis.
    let refusal = root
        .admit_for_promotion(
            &admitted,
            ObservedIdentity::from_receipt(verifier).bytes(),
            forged,
        )
        .expect_err("same policy ID over different bytes must refuse");
    assert!(matches!(
        refusal,
        PromotionRefusal::KeyPolicyObservationMismatch { .. }
    ));
}

#[test]
fn foreign_key_policy_refuses_even_with_the_trusted_verifier() {
    let verifier = verifier_receipt(0x600D);
    let trusted_policy = policy_receipt(0x600D);
    let rogue_policy = policy_receipt(0xBAD);
    let admitted = permit_all_admitted(subject_receipt(5), verifier, rogue_policy);
    let root = Root::configure(
        ObservedIdentity::from_receipt(verifier),
        ObservedIdentity::from_receipt(trusted_policy),
        "promotion-test",
    )
    .expect("root configures");
    assert_eq!(
        root.admit_for_promotion(
            &admitted,
            ObservedIdentity::from_receipt(verifier).bytes(),
            ObservedIdentity::from_receipt(rogue_policy).bytes(),
        )
        .expect_err("a foreign key policy must never mint promotion"),
        PromotionRefusal::ForeignKeyPolicy
    );
}

#[test]
fn an_empty_context_never_configures_a_root() {
    assert_eq!(
        Root::configure(
            ObservedIdentity::from_receipt(verifier_receipt(1)),
            ObservedIdentity::from_receipt(policy_receipt(1)),
            "",
        )
        .expect_err("empty context refuses"),
        PromotionRefusal::EmptyContext
    );

    let oversized_context: &'static str =
        Box::leak("x".repeat(MAX_PROMOTION_CONTEXT_BYTES + 1).into_boxed_str());
    assert_eq!(
        Root::configure(
            ObservedIdentity::from_receipt(verifier_receipt(1)),
            ObservedIdentity::from_receipt(policy_receipt(1)),
            oversized_context,
        )
        .expect_err("oversized context refuses before any transcript is derived"),
        PromotionRefusal::ContextTooLong {
            maximum_bytes: MAX_PROMOTION_CONTEXT_BYTES,
            presented_bytes: MAX_PROMOTION_CONTEXT_BYTES + 1,
        }
    );

    let legacy_verifier = ObservedIdentity::from_receipt(verifier_receipt(1));
    let legacy_policy = ObservedIdentity::from_receipt(policy_receipt(1));
    let legacy_v2 = legacy::promotion_root_charter_v2_for_replay(
        legacy_verifier,
        legacy_policy,
        oversized_context,
    )
    .expect("v3 bounds must not retroactively change exact historical v2 replay");
    assert_eq!(
        legacy_v2.as_bytes(),
        legacy_v2_charter_reference(legacy_verifier, legacy_policy, oversized_context).as_bytes()
    );
}

#[test]
fn owner_execution_requires_nonzero_epoch_and_protocol_versions() {
    let subject = subject_receipt(12);
    let verifier = verifier_receipt(0x600D);
    let policy = policy_receipt(0x600D);
    let configure = |epoch, verifier_descriptor, admission_descriptor| {
        Root::configure_owner_executed(
            ObservedIdentity::from_receipt(verifier),
            ObservedIdentity::from_receipt(policy),
            "promotion-test",
            epoch,
            ExactOwnerVerifier {
                descriptor: verifier_descriptor,
                expected_subject: subject.canonical_preimage(),
                expected_anchor: anchor().content_id(),
                behavior: OwnerBehavior::Exact,
            },
            ExactOwnerAdmitter {
                descriptor: admission_descriptor,
                behavior: OwnerBehavior::Exact,
            },
        )
    };
    assert_eq!(
        configure(0, owner_verifier_descriptor(), owner_admitter_descriptor())
            .expect_err("zero cannot identify an owner policy epoch"),
        PromotionRefusal::InvalidDecisionEpoch
    );
    assert_eq!(
        configure(
            7,
            capability_descriptor(b"owner verifier", b"owner verifier config", 0),
            owner_admitter_descriptor(),
        )
        .expect_err("zero cannot identify a verifier decision protocol"),
        PromotionRefusal::InvalidCapabilityProtocolVersion {
            stage: PromotionCapabilityStage::Verification,
        }
    );
    assert_eq!(
        configure(
            7,
            owner_verifier_descriptor(),
            capability_descriptor(b"owner admitter", b"owner admitter config", 0),
        )
        .expect_err("zero cannot identify an admission decision protocol"),
        PromotionRefusal::InvalidCapabilityProtocolVersion {
            stage: PromotionCapabilityStage::Admission,
        }
    );
}

// ── Root-charter provenance (beads sj31i.52.9 + sj31i.52.11) ──────────────

#[test]
#[allow(clippy::too_many_lines)]
fn charter_v3_matches_independent_reference_and_configuration_axes_move() {
    assert_eq!(
        PROMOTION_ROOT_CHARTER_IDENTITY_VERSION,
        EXPECTED_CURRENT_CHARTER_VERSION
    );
    assert_eq!(
        PROMOTION_ROOT_CHARTER_DOMAIN,
        EXPECTED_CURRENT_CHARTER_DOMAIN
    );
    let verifier = verifier_receipt(0x600D);
    let policy = policy_receipt(0x600D);
    let verifier_observed = ObservedIdentity::from_receipt(verifier);
    let policy_observed = ObservedIdentity::from_receipt(policy);
    let make = |v, p, ctx| {
        Root::configure(v, p, ctx)
            .expect("root configures")
            .charter()
    };
    let baseline = make(verifier_observed, policy_observed, "promotion-test");
    assert_eq!(
        baseline.as_bytes(),
        current_charter_reference(verifier_observed, policy_observed, "promotion-test").as_bytes(),
        "the streamed implementation must match an independent buffered v3 grammar"
    );
    assert_ne!(
        baseline.as_bytes(),
        current_charter_reference_with_roles(
            verifier_observed,
            policy_observed,
            "promotion-test",
            IdentityRole::KeyPolicy,
            IdentityRole::Verifier,
        )
        .as_bytes(),
        "the explicit verifier and key-policy role tags are identity-bearing"
    );

    // Byte-identical configuration => identical charter (the fingerprint is
    // configuration-relative, including a Copy of a configuration-only root).
    let rebuilt = make(verifier_observed, policy_observed, "promotion-test");
    assert_eq!(baseline, rebuilt);
    let root = Root::configure(verifier_observed, policy_observed, "promotion-test")
        .expect("root configures");
    let copied = root;
    assert_eq!(copied.charter(), root.charter());

    // Every retained non-schema axis moves independently. In particular,
    // observation root and exact length are separate inputs rather than one
    // opaque debug record.
    let verifier_bytes = verifier_observed.bytes();
    let policy_bytes = policy_observed.bytes();
    let other_verifier_id = make(
        ObservedIdentity::presented(verifier_receipt(0xBAD).id(), verifier_bytes),
        policy_observed,
        "promotion-test",
    );
    let other_verifier_root = make(
        ObservedIdentity::presented(
            verifier.id(),
            ByteObservation::new(
                ContentId::of_bytes(b"other-verifier-bytes"),
                verifier_bytes.length(),
            ),
        ),
        policy_observed,
        "promotion-test",
    );
    let other_verifier_length = make(
        ObservedIdentity::presented(
            verifier.id(),
            ByteObservation::new(verifier_bytes.content_id(), verifier_bytes.length() + 1),
        ),
        policy_observed,
        "promotion-test",
    );
    let other_policy_id = make(
        verifier_observed,
        ObservedIdentity::presented(policy_receipt(0xBAD).id(), policy_bytes),
        "promotion-test",
    );
    let other_policy_root = make(
        verifier_observed,
        ObservedIdentity::presented(
            policy.id(),
            ByteObservation::new(
                ContentId::of_bytes(b"other-policy-bytes"),
                policy_bytes.length(),
            ),
        ),
        "promotion-test",
    );
    let other_policy_length = make(
        verifier_observed,
        ObservedIdentity::presented(
            policy.id(),
            ByteObservation::new(policy_bytes.content_id(), policy_bytes.length() + 1),
        ),
        "promotion-test",
    );
    let other_context = make(verifier_observed, policy_observed, "promotion-test-other");
    let charters = [
        baseline,
        other_verifier_id,
        other_verifier_root,
        other_verifier_length,
        other_policy_id,
        other_policy_root,
        other_policy_length,
        other_context,
    ];
    for (i, a) in charters.iter().enumerate() {
        for (j, b) in charters.iter().enumerate() {
            if i != j {
                assert_ne!(a, b, "axes {i} and {j} must yield distinct charters");
            }
        }
    }
}

#[test]
#[allow(clippy::too_many_lines)]
fn owner_capability_descriptors_epoch_and_mode_move_v3_charter() {
    let subject = subject_receipt(17);
    let verifier = verifier_receipt(0x600D);
    let policy = policy_receipt(0x600D);
    let profile = Root::configure(
        ObservedIdentity::from_receipt(verifier),
        ObservedIdentity::from_receipt(policy),
        "promotion-test",
    )
    .expect("configuration profile")
    .charter();
    let make = |epoch, verifier_descriptor, admission_descriptor| {
        Root::configure_owner_executed(
            ObservedIdentity::from_receipt(verifier),
            ObservedIdentity::from_receipt(policy),
            "promotion-test",
            epoch,
            ExactOwnerVerifier {
                descriptor: verifier_descriptor,
                expected_subject: subject.canonical_preimage(),
                expected_anchor: anchor().content_id(),
                behavior: OwnerBehavior::Exact,
            },
            ExactOwnerAdmitter {
                descriptor: admission_descriptor,
                behavior: OwnerBehavior::Exact,
            },
        )
        .expect("owner root")
        .charter()
    };
    let verifier_descriptor = owner_verifier_descriptor();
    let admission_descriptor = owner_admitter_descriptor();
    let with_implementation_length = |descriptor: PromotionCapabilityDescriptor| {
        PromotionCapabilityDescriptor::new(
            ByteObservation::new(
                descriptor.implementation().content_id(),
                descriptor.implementation().length() + 1,
            ),
            descriptor.configuration(),
            descriptor.protocol_version(),
        )
    };
    let with_configuration_length = |descriptor: PromotionCapabilityDescriptor| {
        PromotionCapabilityDescriptor::new(
            descriptor.implementation(),
            ByteObservation::new(
                descriptor.configuration().content_id(),
                descriptor.configuration().length() + 1,
            ),
            descriptor.protocol_version(),
        )
    };
    let baseline = make(7, verifier_descriptor, admission_descriptor);
    assert_eq!(
        baseline.as_bytes(),
        owner_charter_reference(
            ObservedIdentity::from_receipt(verifier),
            ObservedIdentity::from_receipt(policy),
            "promotion-test",
            7,
            verifier_descriptor,
            admission_descriptor,
        )
        .as_bytes(),
        "owner-executed v3 charter must match the independent buffered grammar"
    );
    let variants = [
        make(8, verifier_descriptor, admission_descriptor),
        make(
            7,
            capability_descriptor(
                b"other verifier executable",
                b"owner promotion verifier immutable config v1",
                1,
            ),
            admission_descriptor,
        ),
        make(
            7,
            with_implementation_length(verifier_descriptor),
            admission_descriptor,
        ),
        make(
            7,
            capability_descriptor(
                b"owner promotion verifier executable v1",
                b"other verifier configuration",
                1,
            ),
            admission_descriptor,
        ),
        make(
            7,
            with_configuration_length(verifier_descriptor),
            admission_descriptor,
        ),
        make(
            7,
            capability_descriptor(
                b"owner promotion verifier executable v1",
                b"owner promotion verifier immutable config v1",
                2,
            ),
            admission_descriptor,
        ),
        make(
            7,
            verifier_descriptor,
            capability_descriptor(
                b"other admission executable",
                b"owner promotion admission immutable config v1",
                1,
            ),
        ),
        make(
            7,
            verifier_descriptor,
            with_implementation_length(admission_descriptor),
        ),
        make(
            7,
            verifier_descriptor,
            capability_descriptor(
                b"owner promotion admission executable v1",
                b"other admission configuration",
                1,
            ),
        ),
        make(
            7,
            verifier_descriptor,
            with_configuration_length(admission_descriptor),
        ),
        make(
            7,
            verifier_descriptor,
            capability_descriptor(
                b"owner promotion admission executable v1",
                b"owner promotion admission immutable config v1",
                2,
            ),
        ),
    ];
    assert_ne!(baseline, profile, "capability mode is identity-bearing");
    for variant in variants {
        assert_ne!(baseline, variant);
    }
}

#[test]
fn schema_descriptor_axes_move_current_charter_under_reused_domains() {
    let baseline = fixed_charter_root::<CharterSchemaV1, CharterSchemaV1>().charter();
    let variants = [
        (
            "verifier domain",
            fixed_charter_root::<CharterDomainVariant, CharterSchemaV1>().charter(),
        ),
        (
            "verifier name",
            fixed_charter_root::<CharterNameVariant, CharterSchemaV1>().charter(),
        ),
        (
            "verifier version",
            fixed_charter_root::<CharterVersionVariant, CharterSchemaV1>().charter(),
        ),
        (
            "verifier context",
            fixed_charter_root::<CharterContextVariant, CharterSchemaV1>().charter(),
        ),
        (
            "verifier fields",
            fixed_charter_root::<CharterFieldsVariant, CharterSchemaV1>().charter(),
        ),
        (
            "verifier recursive child binding",
            fixed_charter_root::<CharterChildBindingVariant, CharterSchemaV1>().charter(),
        ),
        (
            "key-policy domain",
            fixed_charter_root::<CharterSchemaV1, CharterDomainVariant>().charter(),
        ),
        (
            "key-policy name",
            fixed_charter_root::<CharterSchemaV1, CharterNameVariant>().charter(),
        ),
        (
            "key-policy version",
            fixed_charter_root::<CharterSchemaV1, CharterVersionVariant>().charter(),
        ),
        (
            "key-policy context",
            fixed_charter_root::<CharterSchemaV1, CharterContextVariant>().charter(),
        ),
        (
            "key-policy fields",
            fixed_charter_root::<CharterSchemaV1, CharterFieldsVariant>().charter(),
        ),
        (
            "key-policy recursive child binding",
            fixed_charter_root::<CharterSchemaV1, CharterChildBindingVariant>().charter(),
        ),
    ];

    for (axis, charter) in variants {
        assert_ne!(
            charter, baseline,
            "changing the {axis} must move the current charter"
        );
    }
}

#[test]
fn over_depth_schema_collisions_refuse_before_charter_authority() {
    assert_eq!(
        SchemaId::<CharterOverDepthRootA>::for_schema().as_bytes(),
        SchemaId::<CharterOverDepthRootB>::for_schema().as_bytes(),
        "the depth poison deliberately collapses divergent tails"
    );

    let (verifier_a, policy) = fixed_charter_bindings::<CharterOverDepthRootA, CharterSchemaV1>();
    let legacy_a =
        legacy::promotion_root_charter_v1_for_replay(verifier_a, policy, FIXED_CHARTER_CONTEXT)
            .expect("historical v1 replay bypasses only the current schema-depth guard");
    assert_eq!(
        PromotionTrustRoot::<CharterOverDepthRootA, CharterSchemaV1>::configure(
            verifier_a,
            policy,
            FIXED_CHARTER_CONTEXT,
        )
        .expect_err("a poison-tagged verifier schema cannot mint a current charter"),
        PromotionRefusal::SchemaNestingExceedsCharter {
            role: IdentityRole::Verifier,
            maximum_depth: 16,
        }
    );

    let (verifier_b, policy) = fixed_charter_bindings::<CharterOverDepthRootB, CharterSchemaV1>();
    let legacy_b =
        legacy::promotion_root_charter_v1_for_replay(verifier_b, policy, FIXED_CHARTER_CONTEXT)
            .expect("the divergent historical v1 tail remains replayable");
    assert_eq!(
        PromotionTrustRoot::<CharterOverDepthRootB, CharterSchemaV1>::configure(
            verifier_b,
            policy,
            FIXED_CHARTER_CONTEXT,
        )
        .expect_err("the divergent colliding verifier tail also refuses"),
        PromotionRefusal::SchemaNestingExceedsCharter {
            role: IdentityRole::Verifier,
            maximum_depth: 16,
        }
    );
    assert_eq!(
        legacy_a, legacy_b,
        "faithful v1 replay preserves the historical same-domain collapse"
    );

    let (verifier, over_depth_policy) =
        fixed_charter_bindings::<CharterSchemaV1, CharterOverDepthRootA>();
    assert_eq!(
        PromotionTrustRoot::<CharterSchemaV1, CharterOverDepthRootA>::configure(
            verifier,
            over_depth_policy,
            FIXED_CHARTER_CONTEXT,
        )
        .expect_err("a poison-tagged key-policy schema cannot mint a current charter"),
        PromotionRefusal::SchemaNestingExceedsCharter {
            role: IdentityRole::KeyPolicy,
            maximum_depth: 16,
        }
    );
}

#[test]
#[allow(clippy::too_many_lines)]
fn legacy_v1_and_v2_replay_are_exact_and_nominally_quarantined() {
    assert_eq!(
        LEGACY_PROMOTION_ROOT_CHARTER_V1_DOMAIN,
        EXPECTED_LEGACY_CHARTER_V1_DOMAIN
    );
    assert_eq!(
        LEGACY_PROMOTION_ROOT_CHARTER_V2_DOMAIN,
        EXPECTED_LEGACY_CHARTER_V2_DOMAIN
    );
    let (verifier, key_policy) = fixed_charter_bindings::<CharterSchemaV1, CharterSchemaV1>();
    let baseline_root = PromotionTrustRoot::<CharterSchemaV1, CharterSchemaV1>::configure(
        verifier,
        key_policy,
        FIXED_CHARTER_CONTEXT,
    )
    .expect("baseline root configures");
    let baseline_legacy = baseline_root.legacy_v1_charter_for_replay();
    let baseline_v2 = baseline_root.legacy_v2_charter_for_replay();
    assert_eq!(
        baseline_v2.as_bytes(),
        legacy_v2_charter_reference(verifier, key_policy, FIXED_CHARTER_CONTEXT).as_bytes(),
        "v2 remains exact replay evidence without current authority"
    );
    let reference = legacy_charter_reference(verifier, key_policy, FIXED_CHARTER_CONTEXT);
    assert_eq!(
        baseline_legacy.as_bytes(),
        reference.as_bytes(),
        "the replay wrapper must preserve the exact historical v1 grammar"
    );
    assert_eq!(baseline_legacy.to_string(), baseline_legacy.to_hex());
    assert_eq!(baseline_v2.to_string(), baseline_v2.to_hex());

    // Same-domain schema changes collapsed in v1 because the incomplete
    // grammar did not bind SchemaId. They remain reproducible only through the
    // nominal legacy wrapper; current v3 charters distinguish every case.
    let legacy_collisions = [
        fixed_charter_root::<CharterNameVariant, CharterSchemaV1>().legacy_v1_charter_for_replay(),
        fixed_charter_root::<CharterVersionVariant, CharterSchemaV1>()
            .legacy_v1_charter_for_replay(),
        fixed_charter_root::<CharterContextVariant, CharterSchemaV1>()
            .legacy_v1_charter_for_replay(),
        fixed_charter_root::<CharterFieldsVariant, CharterSchemaV1>()
            .legacy_v1_charter_for_replay(),
        fixed_charter_root::<CharterChildBindingVariant, CharterSchemaV1>()
            .legacy_v1_charter_for_replay(),
        fixed_charter_root::<CharterSchemaV1, CharterNameVariant>().legacy_v1_charter_for_replay(),
        fixed_charter_root::<CharterSchemaV1, CharterVersionVariant>()
            .legacy_v1_charter_for_replay(),
        fixed_charter_root::<CharterSchemaV1, CharterContextVariant>()
            .legacy_v1_charter_for_replay(),
        fixed_charter_root::<CharterSchemaV1, CharterFieldsVariant>()
            .legacy_v1_charter_for_replay(),
        fixed_charter_root::<CharterSchemaV1, CharterChildBindingVariant>()
            .legacy_v1_charter_for_replay(),
    ];
    for replay in legacy_collisions {
        assert_eq!(replay, baseline_legacy);
    }

    let v2_schema_distinctions = [
        fixed_charter_root::<CharterNameVariant, CharterSchemaV1>().legacy_v2_charter_for_replay(),
        fixed_charter_root::<CharterVersionVariant, CharterSchemaV1>()
            .legacy_v2_charter_for_replay(),
        fixed_charter_root::<CharterContextVariant, CharterSchemaV1>()
            .legacy_v2_charter_for_replay(),
        fixed_charter_root::<CharterFieldsVariant, CharterSchemaV1>()
            .legacy_v2_charter_for_replay(),
        fixed_charter_root::<CharterChildBindingVariant, CharterSchemaV1>()
            .legacy_v2_charter_for_replay(),
        fixed_charter_root::<CharterSchemaV1, CharterNameVariant>().legacy_v2_charter_for_replay(),
        fixed_charter_root::<CharterSchemaV1, CharterVersionVariant>()
            .legacy_v2_charter_for_replay(),
        fixed_charter_root::<CharterSchemaV1, CharterContextVariant>()
            .legacy_v2_charter_for_replay(),
        fixed_charter_root::<CharterSchemaV1, CharterFieldsVariant>()
            .legacy_v2_charter_for_replay(),
        fixed_charter_root::<CharterSchemaV1, CharterChildBindingVariant>()
            .legacy_v2_charter_for_replay(),
    ];
    for replay in v2_schema_distinctions {
        assert_ne!(
            replay, baseline_v2,
            "faithful v2 replay retains its complete-schema distinction"
        );
    }

    assert_ne!(
        fixed_charter_root::<CharterNameVariant, CharterSchemaV1>().charter(),
        baseline_root.charter()
    );
    assert_ne!(
        fixed_charter_root::<CharterSchemaV1, CharterNameVariant>().charter(),
        baseline_root.charter()
    );
}

#[test]
fn unscoped_policy_relative_admission_is_a_permanent_downgrade_boundary() {
    let subject = subject_receipt(11);
    let verifier = verifier_receipt(0x600D);
    let policy = policy_receipt(0x600D);
    let admitted = permit_all_admitted(subject, verifier, policy);
    let root = owner_root(
        subject,
        verifier,
        policy,
        7,
        OwnerBehavior::Exact,
        OwnerBehavior::Exact,
    );
    assert_eq!(
        root.admit_for_promotion(
            &admitted,
            ObservedIdentity::from_receipt(verifier).bytes(),
            ObservedIdentity::from_receipt(policy).bytes(),
        )
        .expect_err("the v2-shaped unscoped API never invokes owner code"),
        PromotionRefusal::UnscopedPromotionForbidden
    );
}

#[test]
fn same_public_charter_with_foreign_code_cannot_bind_to_the_owner_root() {
    let subject = subject_receipt(71);
    let verifier = verifier_receipt(0x600D);
    let policy = policy_receipt(0x600D);
    let admitted = permit_all_admitted(subject, verifier, policy);
    let mut owner = owner_root(
        subject,
        verifier,
        policy,
        7,
        OwnerBehavior::Exact,
        OwnerBehavior::Exact,
    );
    let owner_charter = owner.charter();
    let owner_witness = owner
        .decide_for_promotion(
            &admitted,
            ObservedIdentity::from_receipt(verifier).bytes(),
            ObservedIdentity::from_receipt(policy).bytes(),
            decision_scope(7, 1),
        )
        .expect("owner capabilities execute before authority binding");

    // The adversary repeats EVERY public input and even lies with the owner's
    // descriptors, while executing different permit-all capability types.
    let mut foreign = Root::configure_owner_executed(
        ObservedIdentity::from_receipt(verifier),
        ObservedIdentity::from_receipt(policy),
        "promotion-test",
        7,
        ForeignPermitAllVerifier {
            descriptor: owner_verifier_descriptor(),
        },
        ForeignPermitAllAdmitter {
            descriptor: owner_admitter_descriptor(),
        },
    )
    .expect("foreign root can describe itself but gains no owner instance");
    assert_eq!(
        foreign.charter(),
        owner_charter,
        "self-asserted descriptors can collide intentionally; the charter is not authority"
    );
    let foreign_witness = foreign
        .decide_for_promotion(
            &admitted,
            ObservedIdentity::from_receipt(verifier).bytes(),
            ObservedIdentity::from_receipt(policy).bytes(),
            decision_scope(7, 1),
        )
        .expect("foreign root executes only its own permit-all capability");
    assert_eq!(foreign_witness.root_charter(), owner_charter);
    assert_eq!(
        foreign_witness.request_id(),
        owner_witness.request_id(),
        "the roots received the same exact canonical request"
    );
    assert_ne!(
        foreign_witness.verification_decision(),
        owner_witness.verification_decision(),
        "different executed verifier statements must move the stage transcript"
    );
    assert_ne!(
        foreign_witness.admission_decision(),
        owner_witness.admission_decision(),
        "admission binds the exact preceding verifier transcript and its own statement"
    );
    assert_ne!(
        foreign_witness.decision_id(),
        owner_witness.decision_id(),
        "the committed identity distinguishes the actually returned transcripts"
    );
    assert_eq!(
        owner
            .bind_witness(&foreign_witness)
            .expect_err("matching public charter cannot spoof the private live owner root"),
        PromotionRefusal::ForeignRootInstance
    );
}

#[test]
#[allow(clippy::too_many_lines)]
fn decision_identity_binds_every_request_scope_axis() {
    let subject = subject_receipt(78);
    let verifier = verifier_receipt(0x600D);
    let policy = policy_receipt(0x600D);
    let decide = |subject, presented_anchor, scope| {
        let admitted = permit_all_admitted_with_anchor(subject, verifier, policy, presented_anchor);
        let mut root = Root::configure_owner_executed(
            ObservedIdentity::from_receipt(verifier),
            ObservedIdentity::from_receipt(policy),
            "promotion-test",
            7,
            ForeignPermitAllVerifier {
                descriptor: owner_verifier_descriptor(),
            },
            ForeignPermitAllAdmitter {
                descriptor: owner_admitter_descriptor(),
            },
        )
        .expect("permissive transcript fixture configures");
        root.decide_for_promotion(
            &admitted,
            ObservedIdentity::from_receipt(verifier).bytes(),
            ObservedIdentity::from_receipt(policy).bytes(),
            scope,
        )
        .expect("permissive owner fixture exposes transcript mutation axes")
    };
    let attempt = ContentId::of_bytes(b"request-axis-attempt");
    let decision_context = ContentId::of_bytes(b"promotion-test/decision-context/v1");
    let baseline_scope = PromotionDecisionScope::fresh(attempt, decision_context, 7, 1);
    let baseline = decide(subject, anchor(), baseline_scope);
    let variants = [
        decide(subject_receipt(79), anchor(), baseline_scope),
        decide(
            subject,
            ExternalAnchorRef::presented(ContentId::of_bytes(b"request-axis-other-anchor")),
            baseline_scope,
        ),
        decide(
            subject,
            anchor(),
            PromotionDecisionScope::fresh(
                ContentId::of_bytes(b"request-axis-other-attempt"),
                decision_context,
                7,
                1,
            ),
        ),
        decide(
            subject,
            anchor(),
            PromotionDecisionScope::fresh(
                attempt,
                ContentId::of_bytes(b"request-axis-other-context"),
                7,
                1,
            ),
        ),
        decide(
            subject,
            anchor(),
            PromotionDecisionScope::fresh(attempt, decision_context, 7, 2),
        ),
    ];
    for variant in variants {
        assert_eq!(
            variant.root_charter(),
            baseline.root_charter(),
            "request-axis fixtures deliberately share one public root charter"
        );
        assert_ne!(variant.request_id(), baseline.request_id());
        assert_ne!(
            variant.verification_decision(),
            baseline.verification_decision()
        );
        assert_ne!(variant.admission_decision(), baseline.admission_decision());
        assert_ne!(variant.decision_id(), baseline.decision_id());
    }

    // Independently replay the complete request grammar, then perturb every
    // declared request field. This catches an implementation that happens to
    // move on high-level fixture changes while silently omitting a nested
    // role/schema/observation/length axis.
    let audit = baseline.audit();
    let reference = RequestIdentityReference {
        domain: PROMOTION_REQUEST_IDENTITY_DOMAIN,
        version: PROMOTION_REQUEST_IDENTITY_VERSION,
        subject_role: audit.subject_role.tag(),
        subject_schema: audit.subject_schema,
        subject_id: audit.subject_id,
        subject_preimage: audit.subject_preimage,
        subject_bytes: audit.subject_bytes,
        anchor: audit.anchor,
        verifier_role: audit.verifier_role.tag(),
        verifier_schema: audit.verifier_schema,
        verifier_id: audit.verifier_id,
        verifier_observation: audit.verifier_observation,
        key_policy_role: audit.key_policy_role.tag(),
        key_policy_schema: audit.key_policy_schema,
        key_policy_id: audit.key_policy_id,
        key_policy_observation: audit.key_policy_observation,
        root_charter: *audit.root_charter.as_bytes(),
        root_context: audit.context,
        attempt: audit.scope.attempt(),
        decision_context: audit.scope.decision_context(),
        epoch: audit.scope.epoch(),
        sequence: audit.scope.sequence(),
        predecessor: None,
    };
    let reference_id = request_identity_reference(&reference, true);
    assert_eq!(baseline.request_id().as_bytes(), reference_id.as_bytes());
    assert_ne!(
        reference_id,
        request_identity_reference(
            &RequestIdentityReference {
                domain: "org.frankensim.fs-blake3.promotion-request.other.v1",
                ..reference
            },
            true,
        ),
        "request digest domain is identity-bearing"
    );
    assert_ne!(
        reference_id,
        request_identity_reference(&reference, false),
        "request length framing is identity-bearing"
    );

    let field_variants = [
        RequestIdentityReference {
            version: reference.version + 1,
            ..reference
        },
        RequestIdentityReference {
            subject_role: IdentityRole::Entity.tag(),
            ..reference
        },
        RequestIdentityReference {
            subject_schema: changed_digest_bytes(reference.subject_schema),
            ..reference
        },
        RequestIdentityReference {
            subject_id: changed_digest_bytes(reference.subject_id),
            ..reference
        },
        RequestIdentityReference {
            subject_preimage: ContentId::of_bytes(b"other request subject preimage"),
            ..reference
        },
        RequestIdentityReference {
            subject_bytes: reference.subject_bytes + 1,
            ..reference
        },
        RequestIdentityReference {
            anchor: ContentId::of_bytes(b"other request anchor"),
            ..reference
        },
        RequestIdentityReference {
            verifier_role: IdentityRole::Checker.tag(),
            ..reference
        },
        RequestIdentityReference {
            verifier_schema: changed_digest_bytes(reference.verifier_schema),
            ..reference
        },
        RequestIdentityReference {
            verifier_id: changed_digest_bytes(reference.verifier_id),
            ..reference
        },
        RequestIdentityReference {
            verifier_observation: ByteObservation::new(
                ContentId::of_bytes(b"other verifier observation root"),
                reference.verifier_observation.length(),
            ),
            ..reference
        },
        RequestIdentityReference {
            verifier_observation: ByteObservation::new(
                reference.verifier_observation.content_id(),
                reference.verifier_observation.length() + 1,
            ),
            ..reference
        },
        RequestIdentityReference {
            key_policy_role: IdentityRole::Model.tag(),
            ..reference
        },
        RequestIdentityReference {
            key_policy_schema: changed_digest_bytes(reference.key_policy_schema),
            ..reference
        },
        RequestIdentityReference {
            key_policy_id: changed_digest_bytes(reference.key_policy_id),
            ..reference
        },
        RequestIdentityReference {
            key_policy_observation: ByteObservation::new(
                ContentId::of_bytes(b"other policy observation root"),
                reference.key_policy_observation.length(),
            ),
            ..reference
        },
        RequestIdentityReference {
            key_policy_observation: ByteObservation::new(
                reference.key_policy_observation.content_id(),
                reference.key_policy_observation.length() + 1,
            ),
            ..reference
        },
        RequestIdentityReference {
            root_charter: changed_digest_bytes(reference.root_charter),
            ..reference
        },
        RequestIdentityReference {
            root_context: "other request root context",
            ..reference
        },
        RequestIdentityReference {
            attempt: ContentId::of_bytes(b"other request attempt"),
            ..reference
        },
        RequestIdentityReference {
            decision_context: ContentId::of_bytes(b"other request decision context"),
            ..reference
        },
        RequestIdentityReference {
            epoch: reference.epoch + 1,
            ..reference
        },
        RequestIdentityReference {
            sequence: reference.sequence + 1,
            ..reference
        },
        RequestIdentityReference {
            predecessor: Some(*baseline.decision_id().as_bytes()),
            ..reference
        },
    ];
    for variant in field_variants {
        assert_ne!(reference_id, request_identity_reference(&variant, true));
    }
    let predecessor_a = RequestIdentityReference {
        predecessor: Some(*baseline.decision_id().as_bytes()),
        ..reference
    };
    let predecessor_b = RequestIdentityReference {
        predecessor: Some(changed_digest_bytes(*baseline.decision_id().as_bytes())),
        ..reference
    };
    assert_ne!(
        request_identity_reference(&predecessor_a, true),
        request_identity_reference(&predecessor_b, true),
        "a present predecessor's exact decision ID is identity-bearing"
    );
}

#[test]
#[allow(clippy::too_many_lines)]
fn stage_decision_identities_match_independent_complete_grammars() {
    let subject = subject_receipt(81);
    let verifier = verifier_receipt(0x600D);
    let policy = policy_receipt(0x600D);
    let admitted = permit_all_admitted(subject, verifier, policy);
    let mut root = owner_root(
        subject,
        verifier,
        policy,
        7,
        OwnerBehavior::Exact,
        OwnerBehavior::Exact,
    );
    let witness = root
        .decide_for_promotion(
            &admitted,
            ObservedIdentity::from_receipt(verifier).bytes(),
            ObservedIdentity::from_receipt(policy).bytes(),
            decision_scope(7, 1),
        )
        .expect("owner stage-decision fixture");
    let request = *witness.request_id().as_bytes();

    let verification = verification_decision_reference(
        PROMOTION_VERIFICATION_DECISION_IDENTITY_DOMAIN,
        PROMOTION_VERIFICATION_DECISION_IDENTITY_VERSION,
        1,
        request,
        witness.verifier_capability(),
        witness.verification_statement(),
        true,
    );
    assert_eq!(
        witness.verification_decision().as_bytes(),
        verification.as_bytes()
    );
    for variant in descriptor_identity_variants(witness.verifier_capability()) {
        assert_ne!(
            verification,
            verification_decision_reference(
                PROMOTION_VERIFICATION_DECISION_IDENTITY_DOMAIN,
                PROMOTION_VERIFICATION_DECISION_IDENTITY_VERSION,
                1,
                request,
                variant,
                witness.verification_statement(),
                true,
            )
        );
    }
    for variant in [
        verification_decision_reference(
            "org.frankensim.fs-blake3.promotion-verification-decision.other.v1",
            PROMOTION_VERIFICATION_DECISION_IDENTITY_VERSION,
            1,
            request,
            witness.verifier_capability(),
            witness.verification_statement(),
            true,
        ),
        verification_decision_reference(
            PROMOTION_VERIFICATION_DECISION_IDENTITY_DOMAIN,
            PROMOTION_VERIFICATION_DECISION_IDENTITY_VERSION + 1,
            1,
            request,
            witness.verifier_capability(),
            witness.verification_statement(),
            true,
        ),
        verification_decision_reference(
            PROMOTION_VERIFICATION_DECISION_IDENTITY_DOMAIN,
            PROMOTION_VERIFICATION_DECISION_IDENTITY_VERSION,
            2,
            request,
            witness.verifier_capability(),
            witness.verification_statement(),
            true,
        ),
        verification_decision_reference(
            PROMOTION_VERIFICATION_DECISION_IDENTITY_DOMAIN,
            PROMOTION_VERIFICATION_DECISION_IDENTITY_VERSION,
            1,
            changed_digest_bytes(request),
            witness.verifier_capability(),
            witness.verification_statement(),
            true,
        ),
        verification_decision_reference(
            PROMOTION_VERIFICATION_DECISION_IDENTITY_DOMAIN,
            PROMOTION_VERIFICATION_DECISION_IDENTITY_VERSION,
            1,
            request,
            witness.verifier_capability(),
            ContentId::of_bytes(b"other verification statement"),
            true,
        ),
        verification_decision_reference(
            PROMOTION_VERIFICATION_DECISION_IDENTITY_DOMAIN,
            PROMOTION_VERIFICATION_DECISION_IDENTITY_VERSION,
            1,
            request,
            witness.verifier_capability(),
            witness.verification_statement(),
            false,
        ),
    ] {
        assert_ne!(verification, variant);
    }

    let verification_id = *witness.verification_decision().as_bytes();
    let admission = admission_decision_reference(
        PROMOTION_ADMISSION_DECISION_IDENTITY_DOMAIN,
        PROMOTION_ADMISSION_DECISION_IDENTITY_VERSION,
        2,
        request,
        verification_id,
        witness.admission_capability(),
        witness.admission_statement(),
        true,
    );
    assert_eq!(
        witness.admission_decision().as_bytes(),
        admission.as_bytes()
    );
    for variant in descriptor_identity_variants(witness.admission_capability()) {
        assert_ne!(
            admission,
            admission_decision_reference(
                PROMOTION_ADMISSION_DECISION_IDENTITY_DOMAIN,
                PROMOTION_ADMISSION_DECISION_IDENTITY_VERSION,
                2,
                request,
                verification_id,
                variant,
                witness.admission_statement(),
                true,
            )
        );
    }
    for variant in [
        admission_decision_reference(
            "org.frankensim.fs-blake3.promotion-admission-decision.other.v1",
            PROMOTION_ADMISSION_DECISION_IDENTITY_VERSION,
            2,
            request,
            verification_id,
            witness.admission_capability(),
            witness.admission_statement(),
            true,
        ),
        admission_decision_reference(
            PROMOTION_ADMISSION_DECISION_IDENTITY_DOMAIN,
            PROMOTION_ADMISSION_DECISION_IDENTITY_VERSION + 1,
            2,
            request,
            verification_id,
            witness.admission_capability(),
            witness.admission_statement(),
            true,
        ),
        admission_decision_reference(
            PROMOTION_ADMISSION_DECISION_IDENTITY_DOMAIN,
            PROMOTION_ADMISSION_DECISION_IDENTITY_VERSION,
            1,
            request,
            verification_id,
            witness.admission_capability(),
            witness.admission_statement(),
            true,
        ),
        admission_decision_reference(
            PROMOTION_ADMISSION_DECISION_IDENTITY_DOMAIN,
            PROMOTION_ADMISSION_DECISION_IDENTITY_VERSION,
            2,
            changed_digest_bytes(request),
            verification_id,
            witness.admission_capability(),
            witness.admission_statement(),
            true,
        ),
        admission_decision_reference(
            PROMOTION_ADMISSION_DECISION_IDENTITY_DOMAIN,
            PROMOTION_ADMISSION_DECISION_IDENTITY_VERSION,
            2,
            request,
            changed_digest_bytes(verification_id),
            witness.admission_capability(),
            witness.admission_statement(),
            true,
        ),
        admission_decision_reference(
            PROMOTION_ADMISSION_DECISION_IDENTITY_DOMAIN,
            PROMOTION_ADMISSION_DECISION_IDENTITY_VERSION,
            2,
            request,
            verification_id,
            witness.admission_capability(),
            ContentId::of_bytes(b"other admission statement"),
            true,
        ),
        admission_decision_reference(
            PROMOTION_ADMISSION_DECISION_IDENTITY_DOMAIN,
            PROMOTION_ADMISSION_DECISION_IDENTITY_VERSION,
            2,
            request,
            verification_id,
            witness.admission_capability(),
            witness.admission_statement(),
            false,
        ),
    ] {
        assert_ne!(admission, variant);
    }
}

#[test]
#[allow(clippy::too_many_lines)] // the sealed-registry version guard: one grammar, one straight-line replay
fn final_decision_identity_matches_independent_disposition_grammar() {
    let subject = subject_receipt(80);
    let verifier = verifier_receipt(0x600D);
    let policy = policy_receipt(0x600D);
    let admitted = permit_all_admitted(subject, verifier, policy);
    let mut root = owner_root(
        subject,
        verifier,
        policy,
        7,
        OwnerBehavior::Exact,
        OwnerBehavior::Exact,
    );
    let witness = root
        .decide_for_promotion(
            &admitted,
            ObservedIdentity::from_receipt(verifier).bytes(),
            ObservedIdentity::from_receipt(policy).bytes(),
            decision_scope(7, 1),
        )
        .expect("owner decision fixture");
    let reference = |domain: &str,
                     version: u32,
                     request: [u8; 32],
                     verification: [u8; 32],
                     admission: [u8; 32],
                     disposition_tag: u8,
                     framed: bool| {
        let mut preimage = Vec::new();
        let version = version.to_le_bytes();
        let disposition = [disposition_tag];
        for field in [
            version.as_slice(),
            request.as_slice(),
            verification.as_slice(),
            admission.as_slice(),
            disposition.as_slice(),
        ] {
            if framed {
                push_charter_reference_field(&mut preimage, field);
            } else {
                preimage.extend_from_slice(field);
            }
        }
        hash_domain(domain, &preimage)
    };
    let request = *witness.request_id().as_bytes();
    let verification = *witness.verification_decision().as_bytes();
    let admission = *witness.admission_decision().as_bytes();
    let baseline = reference(
        PROMOTION_DECISION_IDENTITY_DOMAIN,
        PROMOTION_DECISION_IDENTITY_VERSION,
        request,
        verification,
        admission,
        1,
        true,
    );
    assert_eq!(
        witness.decision_id().as_bytes(),
        baseline.as_bytes(),
        "committed decision must match the independent version/request/stages/disposition grammar"
    );
    assert_ne!(
        baseline,
        reference(
            PROMOTION_DECISION_IDENTITY_DOMAIN,
            PROMOTION_DECISION_IDENTITY_VERSION + 1,
            request,
            verification,
            admission,
            1,
            true,
        )
    );
    assert_ne!(
        baseline,
        reference(
            "org.frankensim.fs-blake3.promotion-decision.other.v1",
            PROMOTION_DECISION_IDENTITY_VERSION,
            request,
            verification,
            admission,
            1,
            true,
        )
    );
    assert_ne!(
        baseline,
        reference(
            PROMOTION_DECISION_IDENTITY_DOMAIN,
            PROMOTION_DECISION_IDENTITY_VERSION,
            request,
            verification,
            admission,
            2,
            true,
        ),
        "the approved disposition tag is identity-bearing"
    );
    assert_ne!(
        baseline,
        reference(
            PROMOTION_DECISION_IDENTITY_DOMAIN,
            PROMOTION_DECISION_IDENTITY_VERSION,
            request,
            verification,
            admission,
            1,
            false,
        ),
        "length framing is identity-bearing"
    );
    for variant in [
        reference(
            PROMOTION_DECISION_IDENTITY_DOMAIN,
            PROMOTION_DECISION_IDENTITY_VERSION,
            changed_digest_bytes(request),
            verification,
            admission,
            1,
            true,
        ),
        reference(
            PROMOTION_DECISION_IDENTITY_DOMAIN,
            PROMOTION_DECISION_IDENTITY_VERSION,
            request,
            changed_digest_bytes(verification),
            admission,
            1,
            true,
        ),
        reference(
            PROMOTION_DECISION_IDENTITY_DOMAIN,
            PROMOTION_DECISION_IDENTITY_VERSION,
            request,
            verification,
            changed_digest_bytes(admission),
            1,
            true,
        ),
    ] {
        assert_ne!(baseline, variant);
    }
}

#[test]
fn swapped_capability_descriptors_and_wrong_contexts_refuse() {
    let subject = subject_receipt(72);
    let verifier = verifier_receipt(0x600D);
    let policy = policy_receipt(0x600D);
    let admitted = permit_all_admitted(subject, verifier, policy);
    let owner = owner_root(
        subject,
        verifier,
        policy,
        7,
        OwnerBehavior::Exact,
        OwnerBehavior::Exact,
    );
    let mut swapped = Root::configure_owner_executed(
        ObservedIdentity::from_receipt(verifier),
        ObservedIdentity::from_receipt(policy),
        "promotion-test",
        7,
        ExactOwnerVerifier {
            descriptor: owner_admitter_descriptor(),
            expected_subject: subject.canonical_preimage(),
            expected_anchor: anchor().content_id(),
            behavior: OwnerBehavior::Exact,
        },
        ForeignPermitAllAdmitter {
            descriptor: owner_verifier_descriptor(),
        },
    )
    .expect("typed stages configure, but swapped descriptors remain visible");
    assert_ne!(swapped.charter(), owner.charter());
    let swapped_witness = swapped
        .decide_for_promotion(
            &admitted,
            ObservedIdentity::from_receipt(verifier).bytes(),
            ObservedIdentity::from_receipt(policy).bytes(),
            decision_scope(7, 1),
        )
        .expect("swapped test capability is internally permissive");
    assert_eq!(
        owner
            .bind_witness(&swapped_witness)
            .expect_err("swapped executable stages cannot bind as owner authority"),
        PromotionRefusal::ForeignRootInstance
    );

    let mut wrong_context = Root::configure_owner_executed(
        ObservedIdentity::from_receipt(verifier),
        ObservedIdentity::from_receipt(policy),
        "promotion-test-other",
        7,
        ExactOwnerVerifier {
            descriptor: owner_verifier_descriptor(),
            expected_subject: subject.canonical_preimage(),
            expected_anchor: anchor().content_id(),
            behavior: OwnerBehavior::Exact,
        },
        ExactOwnerAdmitter {
            descriptor: owner_admitter_descriptor(),
            behavior: OwnerBehavior::Exact,
        },
    )
    .expect("wrong-context root is descriptive but owner verifier still executes");
    assert!(matches!(
        wrong_context
            .decide_for_promotion(
                &admitted,
                ObservedIdentity::from_receipt(verifier).bytes(),
                ObservedIdentity::from_receipt(policy).bytes(),
                decision_scope(7, 1),
            )
            .expect_err("stored verifier rejects the wrong root context"),
        PromotionRefusal::CapabilityRefused {
            stage: PromotionCapabilityStage::Verification,
            ..
        }
    ));
}

#[test]
fn wrong_anchor_epoch_and_decision_context_refuse_without_publication() {
    let subject = subject_receipt(73);
    let verifier = verifier_receipt(0x600D);
    let policy = policy_receipt(0x600D);
    let mut root = owner_root(
        subject,
        verifier,
        policy,
        7,
        OwnerBehavior::Exact,
        OwnerBehavior::Exact,
    );
    let wrong_anchor = ExternalAnchorRef::presented(ContentId::of_bytes(b"wrong-anchor"));
    let admitted = permit_all_admitted_with_anchor(subject, verifier, policy, wrong_anchor);
    assert!(matches!(
        root.decide_for_promotion(
            &admitted,
            ObservedIdentity::from_receipt(verifier).bytes(),
            ObservedIdentity::from_receipt(policy).bytes(),
            decision_scope(7, 1),
        )
        .expect_err("owner verifier rejects wrong anchor"),
        PromotionRefusal::CapabilityRefused {
            stage: PromotionCapabilityStage::Verification,
            ..
        }
    ));

    let admitted = permit_all_admitted(subject, verifier, policy);
    assert_eq!(
        root.decide_for_promotion(
            &admitted,
            ObservedIdentity::from_receipt(verifier).bytes(),
            ObservedIdentity::from_receipt(policy).bytes(),
            decision_scope(8, 2),
        )
        .expect_err("wrong epoch refuses before owner execution"),
        PromotionRefusal::WrongDecisionEpoch {
            configured: 7,
            presented: 8,
        }
    );
    let wrong_context_scope = PromotionDecisionScope::fresh(
        ContentId::of_bytes(b"wrong-context-attempt"),
        ContentId::of_bytes(b"wrong-decision-context"),
        7,
        2,
    );
    assert!(matches!(
        root.decide_for_promotion(
            &admitted,
            ObservedIdentity::from_receipt(verifier).bytes(),
            ObservedIdentity::from_receipt(policy).bytes(),
            wrong_context_scope,
        )
        .expect_err("owner verifier rejects wrong decision context"),
        PromotionRefusal::CapabilityRefused {
            stage: PromotionCapabilityStage::Verification,
            ..
        }
    ));
}

#[test]
#[allow(clippy::too_many_lines)]
fn stale_replay_and_cancelled_attempts_are_burned_before_capability_calls() {
    let subject = subject_receipt(74);
    let verifier = verifier_receipt(0x600D);
    let policy = policy_receipt(0x600D);
    let admitted = permit_all_admitted(subject, verifier, policy);
    let mut root = owner_root(
        subject,
        verifier,
        policy,
        7,
        OwnerBehavior::Cancel,
        OwnerBehavior::Exact,
    );
    let stable_charter = root.charter();
    assert!(matches!(
        root.decide_for_promotion(
            &admitted,
            ObservedIdentity::from_receipt(verifier).bytes(),
            ObservedIdentity::from_receipt(policy).bytes(),
            decision_scope(7, 10),
        )
        .expect_err("cancelled verifier publishes no witness"),
        PromotionRefusal::CapabilityCancelled {
            stage: PromotionCapabilityStage::Verification,
            ..
        }
    ));
    assert_eq!(
        root.decide_for_promotion(
            &admitted,
            ObservedIdentity::from_receipt(verifier).bytes(),
            ObservedIdentity::from_receipt(policy).bytes(),
            decision_scope(7, 10),
        )
        .expect_err("cancelled sequence was burned"),
        PromotionRefusal::StaleOrReplayedDecision {
            last_attempted: 10,
            presented: 10,
        }
    );
    assert_eq!(
        root.charter(),
        stable_charter,
        "burned replay state is deliberately not a portable charter input"
    );
    assert_eq!(
        root.decide_for_promotion(
            &admitted,
            ObservedIdentity::from_receipt(verifier).bytes(),
            ObservedIdentity::from_receipt(policy).bytes(),
            decision_scope(7, 9),
        )
        .expect_err("older sequence is stale"),
        PromotionRefusal::StaleOrReplayedDecision {
            last_attempted: 10,
            presented: 9,
        }
    );

    let mut admission_refusal = owner_root(
        subject,
        verifier,
        policy,
        7,
        OwnerBehavior::Exact,
        OwnerBehavior::Refuse,
    );
    assert!(matches!(
        admission_refusal
            .decide_for_promotion(
                &admitted,
                ObservedIdentity::from_receipt(verifier).bytes(),
                ObservedIdentity::from_receipt(policy).bytes(),
                decision_scope(7, 1),
            )
            .expect_err("admission refusal occurs after owner verification"),
        PromotionRefusal::CapabilityRefused {
            stage: PromotionCapabilityStage::Admission,
            ..
        }
    ));
    assert!(matches!(
        admission_refusal
            .decide_for_promotion(
                &admitted,
                ObservedIdentity::from_receipt(verifier).bytes(),
                ObservedIdentity::from_receipt(policy).bytes(),
                decision_scope(7, 1),
            )
            .expect_err("admission-refused sequence was burned"),
        PromotionRefusal::StaleOrReplayedDecision { .. }
    ));

    let mut admission_cancellation = owner_root(
        subject,
        verifier,
        policy,
        7,
        OwnerBehavior::Exact,
        OwnerBehavior::Cancel,
    );
    assert!(matches!(
        admission_cancellation
            .decide_for_promotion(
                &admitted,
                ObservedIdentity::from_receipt(verifier).bytes(),
                ObservedIdentity::from_receipt(policy).bytes(),
                decision_scope(7, 4),
            )
            .expect_err("admission cancellation occurs after owner verification"),
        PromotionRefusal::CapabilityCancelled {
            stage: PromotionCapabilityStage::Admission,
            ..
        }
    ));
    assert_eq!(
        admission_cancellation
            .decide_for_promotion(
                &admitted,
                ObservedIdentity::from_receipt(verifier).bytes(),
                ObservedIdentity::from_receipt(policy).bytes(),
                decision_scope(7, 4),
            )
            .expect_err("admission-cancelled sequence was burned"),
        PromotionRefusal::StaleOrReplayedDecision {
            last_attempted: 4,
            presented: 4,
        }
    );
}

#[test]
fn verifier_and_admission_panics_poison_without_partial_witnesses() {
    let subject = subject_receipt(75);
    let verifier = verifier_receipt(0x600D);
    let policy = policy_receipt(0x600D);
    let admitted = permit_all_admitted(subject, verifier, policy);
    for (verifier_behavior, admission_behavior, stage) in [
        (
            OwnerBehavior::Panic,
            OwnerBehavior::Exact,
            PromotionCapabilityStage::Verification,
        ),
        (
            OwnerBehavior::Exact,
            OwnerBehavior::Panic,
            PromotionCapabilityStage::Admission,
        ),
    ] {
        let mut root = owner_root(
            subject,
            verifier,
            policy,
            7,
            verifier_behavior,
            admission_behavior,
        );
        let stable_charter = root.charter();
        assert_eq!(
            root.decide_for_promotion(
                &admitted,
                ObservedIdentity::from_receipt(verifier).bytes(),
                ObservedIdentity::from_receipt(policy).bytes(),
                decision_scope(7, 1),
            )
            .expect_err("panic is caught before publication"),
            PromotionRefusal::CapabilityPanicked { stage }
        );
        assert_eq!(
            root.decide_for_promotion(
                &admitted,
                ObservedIdentity::from_receipt(verifier).bytes(),
                ObservedIdentity::from_receipt(policy).bytes(),
                decision_scope(7, 2),
            )
            .expect_err("faulted capability state is never reused"),
            PromotionRefusal::RootPoisoned
        );
        assert_eq!(
            root.charter(),
            stable_charter,
            "runtime poison state is deliberately not a portable charter input"
        );
    }
}

#[test]
#[allow(clippy::too_many_lines)]
fn equivalent_reconstruction_requires_an_owner_executed_crosswalk() {
    let subject = subject_receipt(76);
    let verifier = verifier_receipt(0x600D);
    let policy = policy_receipt(0x600D);
    let admitted = permit_all_admitted(subject, verifier, policy);
    let mut source_root = owner_root(
        subject,
        verifier,
        policy,
        7,
        OwnerBehavior::Exact,
        OwnerBehavior::Exact,
    );
    let source_witness = source_root
        .decide_for_promotion(
            &admitted,
            ObservedIdentity::from_receipt(verifier).bytes(),
            ObservedIdentity::from_receipt(policy).bytes(),
            decision_scope(7, 1),
        )
        .expect("source owner decision");
    let source_bound = source_root
        .bind_witness(&source_witness)
        .expect("source live binding");

    let mut rebuilt = owner_root(
        subject,
        verifier,
        policy,
        7,
        OwnerBehavior::Exact,
        OwnerBehavior::Exact,
    );
    assert_eq!(rebuilt.charter(), source_root.charter());
    assert_eq!(
        rebuilt
            .bind_witness(&source_witness)
            .expect_err("same public reconstruction is a different live root"),
        PromotionRefusal::ForeignRootInstance
    );

    let forged_crosswalk_scope = PromotionDecisionScope::crosswalk(
        ContentId::of_bytes(b"forged-direct-crosswalk-attempt"),
        ContentId::of_bytes(b"promotion-test/decision-context/v1"),
        7,
        1,
        source_bound.decision_id(),
    );
    assert_eq!(
        rebuilt
            .decide_for_promotion(
                &admitted,
                ObservedIdentity::from_receipt(verifier).bytes(),
                ObservedIdentity::from_receipt(policy).bytes(),
                forged_crosswalk_scope,
            )
            .expect_err("a predecessor ID without a live-bound source is not a crosswalk"),
        PromotionRefusal::CrosswalkPredecessorMismatch
    );

    let bad_scope = PromotionDecisionScope::crosswalk(
        ContentId::of_bytes(b"bad-crosswalk-attempt"),
        ContentId::of_bytes(b"promotion-test/decision-context/v1"),
        7,
        1,
        // A different real decision-shaped value, not caller text.
        {
            let other_subject = subject_receipt(77);
            let other_admitted = permit_all_admitted(other_subject, verifier, policy);
            let mut other_root = owner_root(
                other_subject,
                verifier,
                policy,
                7,
                OwnerBehavior::Exact,
                OwnerBehavior::Exact,
            );
            other_root
                .decide_for_promotion(
                    &other_admitted,
                    ObservedIdentity::from_receipt(verifier).bytes(),
                    ObservedIdentity::from_receipt(policy).bytes(),
                    decision_scope(7, 1),
                )
                .expect("other decision")
                .decision_id()
        },
    );
    assert_eq!(
        rebuilt
            .crosswalk_witness(&source_bound, bad_scope)
            .expect_err("wrong predecessor cannot start owner code"),
        PromotionRefusal::CrosswalkPredecessorMismatch
    );

    let crosswalk_scope = PromotionDecisionScope::crosswalk(
        ContentId::of_bytes(b"valid-crosswalk-attempt"),
        ContentId::of_bytes(b"promotion-test/decision-context/v1"),
        7,
        1,
        source_bound.decision_id(),
    );
    let rebuilt_witness = rebuilt
        .crosswalk_witness(&source_bound, crosswalk_scope)
        .expect("target owner re-executes both capabilities for exact predecessor");
    let rebuilt_bound = rebuilt
        .bind_witness(&rebuilt_witness)
        .expect("crosswalk result binds to the target live root");
    assert_eq!(rebuilt_bound.subject(), subject);
    assert_ne!(rebuilt_bound.decision_id(), source_bound.decision_id());
}
