//! Promotion-authority battery (bead sj31i.52.9): a foreign
//! permit-everything capability can drive the generic
//! `Presented → Verified → Admitted` ladder — that admission is
//! explicitly policy-relative — but it can never mint the opaque
//! [`PromotionWitness`]; promotion-capable admission exists only
//! through an independently configured [`PromotionTrustRoot`] that
//! re-adjudicates exact verifier/key-policy identities AND their
//! canonical-byte observations. Compile-fail proofs (private-field
//! construction, sealed typestates) live as `compile_fail` doctests on
//! [`PromotionWitness`].

use fs_blake3::identity::{
    AuthorityAdmitter, AuthorityRef, AuthorityVerifier, ByteObservation, CanonicalEncoder,
    CanonicalLimits, CanonicalSchema, ChildSpec, ContentId, ExternalAnchorRef, Field, FieldSpec,
    IdentityReceipt, IdentityRole, KeyPolicyId, LEGACY_PROMOTION_ROOT_CHARTER_V1_DOMAIN,
    NeverCancel, ObservedIdentity, PROMOTION_ROOT_CHARTER_DOMAIN,
    PROMOTION_ROOT_CHARTER_IDENTITY_VERSION, Presented, PromotionRefusal, PromotionTrustRoot,
    PromotionWitness, SchemaId, SemanticId, StrongIdentity, Verified, VerifierId, WireType, legacy,
};
use fs_blake3::{ContentHash, hash_bytes, hash_domain};

const LIMITS: CanonicalLimits = CanonicalLimits::new(64 * 1024, 16 * 1024, 64, 1024, 7);
const EXPECTED_CURRENT_CHARTER_VERSION: u32 = 2;
const EXPECTED_CURRENT_CHARTER_DOMAIN: &str = "org.frankensim.fs-blake3.promotion-root-charter.v2";
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
        context.as_bytes(),
    ] {
        push_charter_reference_field(&mut preimage, field);
    }
    hash_domain(EXPECTED_CURRENT_CHARTER_DOMAIN, &preimage)
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
    AuthorityRef::present(subject, anchor(), verifier.id(), policy.id())
        .verify(&PermitAll)
        .expect("permit-all verifies anything")
        .admit(&PermitAll)
        .expect("permit-all admits anything")
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
}

#[test]
fn the_configured_root_mints_a_fully_bound_witness() {
    let subject = subject_receipt(11);
    let verifier = verifier_receipt(0x600D);
    let policy = policy_receipt(0x600D);
    let admitted = permit_all_admitted(subject, verifier, policy);
    let root = Root::configure(
        ObservedIdentity::from_receipt(verifier),
        ObservedIdentity::from_receipt(policy),
        "promotion-test",
    )
    .expect("root configures");
    let witness: Witness = root
        .admit_for_promotion(
            &admitted,
            ObservedIdentity::from_receipt(verifier).bytes(),
            ObservedIdentity::from_receipt(policy).bytes(),
        )
        .expect("the exact configured binding promotes");
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
    // Bounded audit: namespaces plus observation roots/lengths only.
    let audit = witness.audit();
    assert_eq!(audit.verifier_domain, VerifierSchemaV1::DOMAIN);
    assert_eq!(audit.key_policy_domain, PolicySchemaV1::DOMAIN);
    assert_eq!(
        audit.verifier_observation,
        ObservedIdentity::from_receipt(verifier).bytes()
    );
    assert_eq!(audit.context, "promotion-test");
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
}

// ── Root-charter provenance (beads sj31i.52.9 + sj31i.52.11) ──────────────

#[test]
#[allow(clippy::too_many_lines)]
fn charter_v2_matches_independent_reference_and_existing_axes_move() {
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
        "the streamed implementation must match an independent buffered v2 grammar"
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

    // Byte-identical configuration => identical charter (provenance is
    // configuration-relative: two roots making identical decisions ARE the
    // same policy — including a Copy of the root).
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
fn legacy_v1_replay_is_exact_and_nominally_quarantined() {
    assert_eq!(
        LEGACY_PROMOTION_ROOT_CHARTER_V1_DOMAIN,
        EXPECTED_LEGACY_CHARTER_V1_DOMAIN
    );
    let (verifier, key_policy) = fixed_charter_bindings::<CharterSchemaV1, CharterSchemaV1>();
    let baseline_root = PromotionTrustRoot::<CharterSchemaV1, CharterSchemaV1>::configure(
        verifier,
        key_policy,
        FIXED_CHARTER_CONTEXT,
    )
    .expect("baseline root configures");
    let baseline_legacy = baseline_root.legacy_v1_charter_for_replay();
    let reference = legacy_charter_reference(verifier, key_policy, FIXED_CHARTER_CONTEXT);
    assert_eq!(
        baseline_legacy.as_bytes(),
        reference.as_bytes(),
        "the replay wrapper must preserve the exact historical v1 grammar"
    );
    assert_eq!(baseline_legacy.to_string(), baseline_legacy.to_hex());

    // Same-domain schema changes collapsed in v1 because the incomplete
    // grammar did not bind SchemaId. They remain reproducible only through the
    // nominal legacy wrapper; current v2 charters distinguish every case.
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
fn witnesses_carry_the_minting_roots_exact_charter() {
    let subject = subject_receipt(11);
    let verifier = verifier_receipt(0x600D);
    let policy = policy_receipt(0x600D);
    let admitted = permit_all_admitted(subject, verifier, policy);
    let root = Root::configure(
        ObservedIdentity::from_receipt(verifier),
        ObservedIdentity::from_receipt(policy),
        "promotion-test",
    )
    .expect("root configures");
    let witness: Witness = root
        .admit_for_promotion(
            &admitted,
            ObservedIdentity::from_receipt(verifier).bytes(),
            ObservedIdentity::from_receipt(policy).bytes(),
        )
        .expect("the exact configured binding promotes");
    assert_eq!(witness.root_charter(), root.charter());
}

#[test]
fn a_self_configured_root_mints_witnesses_with_a_visibly_foreign_charter() {
    // The domain owner pins the charter of ITS root (in a golden, a
    // CONTRACT constant, or a composition-root check).
    let owner_verifier = verifier_receipt(0x600D);
    let owner_policy = policy_receipt(0x600D);
    let owner_root = Root::configure(
        ObservedIdentity::from_receipt(owner_verifier),
        ObservedIdentity::from_receipt(owner_policy),
        "promotion-test",
    )
    .expect("owner root configures");
    let pinned = owner_root.charter();

    // The adversary can freely configure its OWN root around its rogue
    // permit-everything verifier and mint a structurally genuine witness —
    // configure() is public and that is by design.
    let rogue_verifier = verifier_receipt(0xBAD);
    let rogue_policy = policy_receipt(0xBAD);
    let rogue_admitted = permit_all_admitted(subject_receipt(7), rogue_verifier, rogue_policy);
    let rogue_root = Root::configure(
        ObservedIdentity::from_receipt(rogue_verifier),
        ObservedIdentity::from_receipt(rogue_policy),
        "promotion-test",
    )
    .expect("rogue root configures");
    let rogue_witness = rogue_root
        .admit_for_promotion(
            &rogue_admitted,
            ObservedIdentity::from_receipt(rogue_verifier).bytes(),
            ObservedIdentity::from_receipt(rogue_policy).bytes(),
        )
        .expect("a self-configured root promotes its own binding");

    // ...but the witness carries the rogue configuration's charter, so the
    // pinning consumption boundary refuses it. THIS is the root-provenance
    // closure: self-configured and cloned-then-reconfigured roots are
    // distinguishable exactly where witnesses are consumed.
    assert_ne!(rogue_witness.root_charter(), pinned);

    // An honestly-minted witness from the owner's configuration passes the
    // same pin.
    let honest_admitted = permit_all_admitted(subject_receipt(8), owner_verifier, owner_policy);
    let honest_witness = owner_root
        .admit_for_promotion(
            &honest_admitted,
            ObservedIdentity::from_receipt(owner_verifier).bytes(),
            ObservedIdentity::from_receipt(owner_policy).bytes(),
        )
        .expect("owner binding promotes");
    assert_eq!(honest_witness.root_charter(), pinned);
}
