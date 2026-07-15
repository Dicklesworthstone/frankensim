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
    CanonicalLimits, CanonicalSchema, ContentId, ExternalAnchorRef, Field, FieldSpec,
    IdentityReceipt, KeyPolicyId, NeverCancel, ObservedIdentity, Presented, PromotionRefusal,
    PromotionTrustRoot, PromotionWitness, SemanticId, Verified, VerifierId, WireType,
};

const LIMITS: CanonicalLimits = CanonicalLimits::new(64 * 1024, 16 * 1024, 64, 1024, 7);

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

type Subject = SemanticId<SubjectV1>;
type PresentedAuthority = AuthorityRef<Subject, VerifierSchemaV1, PolicySchemaV1, Presented>;
type VerifiedAuthority = AuthorityRef<Subject, VerifierSchemaV1, PolicySchemaV1, Verified>;
type Witness = PromotionWitness<Subject, VerifierSchemaV1, PolicySchemaV1>;
type Root = PromotionTrustRoot<VerifierSchemaV1, PolicySchemaV1>;

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
