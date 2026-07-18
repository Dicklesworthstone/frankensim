//! G0/G3/G5 conformance battery for proposition-bound spectral admission and truth.

#![allow(clippy::wildcard_imports)]
#![allow(
    clippy::too_many_lines,
    reason = "each long G0/G3/G5 fixture keeps one end-to-end authority and refusal narrative auditable"
)]

use fs_blake3::identity::{
    AuthorityAdmitter, AuthorityRef, AuthorityVerifier, ByteObservation, CanonicalSchema,
    ContentId, ExternalAnchorRef, IdentityAdjudication, IdentityReceipt, NoClaimState,
    ObservedIdentity, PromotionRefusal, PromotionRootCharter, StrongIdentity as _, TrustState,
    adjudicate,
};
use fs_qty::{Angle, Dims, QtyAny, Time};
use fs_spectral::admission::*;
use fs_spectral::truth::*;

fn subject(byte: u8) -> SpectralSubjectId {
    SpectralSubjectId::from_bytes([byte; 32])
}

fn form_id(byte: u8) -> SpectralFormId {
    SpectralFormId::from_bytes([byte; 32])
}

fn norm_id(byte: u8) -> SpectralNormId {
    SpectralNormId::from_bytes([byte; 32])
}

fn scaling_id(byte: u8) -> SpectralScalingId {
    SpectralScalingId::from_bytes([byte; 32])
}

fn scaling_map(byte: u8) -> SpectralScalingMapId {
    SpectralScalingMapId::from_bytes([byte; 32])
}

fn function_id(byte: u8) -> SpectralFunctionId {
    SpectralFunctionId::from_bytes([byte; 32])
}

fn continuation_id(byte: u8) -> SpectralContinuationId {
    SpectralContinuationId::from_bytes([byte; 32])
}

fn gauge_artifact_id(byte: u8) -> SpectralGaugeArtifactId {
    SpectralGaugeArtifactId::from_bytes([byte; 32])
}

fn quotient_map_id(byte: u8) -> SpectralQuotientMapId {
    SpectralQuotientMapId::from_bytes([byte; 32])
}

fn cluster_id(byte: u8) -> SpectralClusterIdV1 {
    SpectralClusterIdV1::from_bytes([byte; 32])
}

fn indexed_cluster_id(index: usize) -> SpectralClusterIdV1 {
    let mut bytes = [0_u8; 32];
    bytes[..8].copy_from_slice(&u64::try_from(index).unwrap().to_le_bytes());
    SpectralClusterIdV1::from_bytes(bytes)
}

fn region_id(byte: u8) -> SpectralRegionId {
    SpectralRegionId::from_bytes([byte; 32])
}

fn scaling(dims: Dims, scale: f64, identity_byte: u8) -> SpectralScalingContextV1 {
    SpectralScalingContextV1::new(
        scaling_id(identity_byte),
        dims,
        scale,
        scaling_map(identity_byte.wrapping_add(1)),
        scaling_map(identity_byte.wrapping_add(2)),
        scaling_map(identity_byte.wrapping_add(3)),
        scaling_map(identity_byte.wrapping_add(4)),
    )
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum AuthorityError {
    Mismatch,
}

#[derive(Clone, Copy)]
struct ExactAuthority {
    proposition: SpectralPropositionId,
    preimage: ContentId,
    anchor: ExternalAnchorRef,
    verifier: IdentityReceipt<SpectralAuthorityVerifierIdV1>,
    policy: IdentityReceipt<SpectralAuthorityPolicyIdV1>,
}

impl
    AuthorityVerifier<
        SpectralPropositionId,
        SpectralVerifierIdentitySchemaV1,
        SpectralAuthorityPolicySchemaV1,
    > for ExactAuthority
{
    type Error = AuthorityError;

    fn verify(&self, presented: &PresentedSpectralAuthorityV1) -> Result<(), Self::Error> {
        let receipt = presented.receipt();
        if receipt.id() == self.proposition
            && receipt.canonical_preimage() == self.preimage
            && presented.anchor() == self.anchor
            && presented.verifier() == self.verifier.id()
            && presented.key_policy() == self.policy.id()
        {
            Ok(())
        } else {
            Err(AuthorityError::Mismatch)
        }
    }
}

impl
    AuthorityAdmitter<
        SpectralPropositionId,
        SpectralVerifierIdentitySchemaV1,
        SpectralAuthorityPolicySchemaV1,
    > for ExactAuthority
{
    type Error = AuthorityError;

    fn admit(&self, verified: &VerifiedSpectralAuthorityV1) -> Result<(), Self::Error> {
        let receipt = verified.receipt();
        if receipt.id() == self.proposition
            && receipt.canonical_preimage() == self.preimage
            && verified.anchor() == self.anchor
            && verified.verifier() == self.verifier.id()
            && verified.key_policy() == self.policy.id()
        {
            Ok(())
        } else {
            Err(AuthorityError::Mismatch)
        }
    }
}

fn exact_authority(receipt: IdentityReceipt<SpectralPropositionId>, seed: u8) -> ExactAuthority {
    ExactAuthority {
        proposition: receipt.id(),
        preimage: receipt.canonical_preimage(),
        anchor: ExternalAnchorRef::presented(ContentId::of_bytes(&[b's', b'p', b'e', b'c', seed])),
        verifier: spectral_verifier_receipt(b"fs-spectral-test-exact-verifier-v1").unwrap(),
        policy: spectral_authority_policy_receipt(b"fs-spectral-test-admission-policy-v1").unwrap(),
    }
}

fn present(
    receipt: IdentityReceipt<SpectralPropositionId>,
    authority: ExactAuthority,
) -> PresentedSpectralAuthorityV1 {
    AuthorityRef::present(
        receipt,
        authority.anchor,
        authority.verifier.id(),
        authority.policy.id(),
    )
}

fn policy_relative_and_promotion(
    receipt: IdentityReceipt<SpectralPropositionId>,
    authority: ExactAuthority,
) -> (
    AdmittedSpectralAuthorityV1,
    SpectralPromotionWitnessV1,
    PromotionRootCharter,
) {
    let admitted = present(receipt, authority)
        .verify(&authority)
        .unwrap()
        .admit(&authority)
        .unwrap();
    let root = spectral_promotion_trust_root(authority.verifier, authority.policy).unwrap();
    let promotion = root
        .admit_for_promotion(
            &admitted,
            ObservedIdentity::from_receipt(authority.verifier).bytes(),
            ObservedIdentity::from_receipt(authority.policy).bytes(),
        )
        .unwrap();
    (admitted, promotion, root.charter())
}

fn admit_receipt(
    receipt: IdentityReceipt<SpectralPropositionId>,
    seed: u8,
) -> AdmittedSpectralWitnessV1 {
    let (admitted, promotion, charter) =
        policy_relative_and_promotion(receipt, exact_authority(receipt, seed));
    AdmittedSpectralWitnessV1::from_authority(&admitted, promotion, charter).unwrap()
}

/// An untrusted caller-controlled verifier/admitter pair that accepts every
/// presented binding. It can reach only policy-relative admission unless a
/// separately configured promotion root also accepts the binding.
struct PermitAll;

impl
    AuthorityVerifier<
        SpectralPropositionId,
        SpectralVerifierIdentitySchemaV1,
        SpectralAuthorityPolicySchemaV1,
    > for PermitAll
{
    type Error = core::convert::Infallible;

    fn verify(&self, _presented: &PresentedSpectralAuthorityV1) -> Result<(), Self::Error> {
        Ok(())
    }
}

impl
    AuthorityAdmitter<
        SpectralPropositionId,
        SpectralVerifierIdentitySchemaV1,
        SpectralAuthorityPolicySchemaV1,
    > for PermitAll
{
    type Error = core::convert::Infallible;

    fn admit(&self, _verified: &VerifiedSpectralAuthorityV1) -> Result<(), Self::Error> {
        Ok(())
    }
}

#[derive(Clone, Copy)]
struct Axes {
    subject: SpectralSubjectId,
    scalar: SpectralScalarFieldV1,
    class: SpectralProblemClassV1,
    scaling: SpectralScalingContextV1,
    domain: SpectralMetricV1,
    codomain: SpectralMetricV1,
}

fn axes(
    seed: u8,
    class: SpectralProblemClassV1,
    scalar: SpectralScalarFieldV1,
    dims: Dims,
    dimension: u32,
) -> Axes {
    let metric = SpectralMetricV1::euclidean(dimension);
    Axes {
        subject: subject(seed),
        scalar,
        class,
        scaling: scaling(dims, 1.0, seed.wrapping_add(20)),
        domain: metric,
        codomain: metric,
    }
}

fn standard_axes(seed: u8, dimension: u32) -> Axes {
    axes(
        seed,
        SpectralProblemClassV1::new(
            SpectralRepresentationV1::StandardLinear,
            DescriptorRoleV1::Ordinary,
            SpectralOperatorOriginV1::Direct,
        ),
        SpectralScalarFieldV1::Real,
        Dims::NONE,
        dimension,
    )
}

#[allow(clippy::too_many_arguments)]
fn problem_spec(
    axes: Axes,
    structures: Vec<StructureClaimV1>,
    regularity: Vec<RegularityClaimV1>,
    spaces: SpectralSpaceContextV1,
    ordering: SpectralOrderingV1,
    scope: CompletenessScopeV1,
) -> SpectralProblemSpecV1 {
    SpectralProblemSpecV1::new(
        axes.subject,
        axes.scalar,
        axes.class,
        StructureProfileV1::new(structures),
        axes.scaling,
        spaces,
        RegularityProfileV1::new(regularity),
        ordering,
        scope,
    )
}

fn default_spec(
    axes: Axes,
    structures: Vec<StructureClaimV1>,
    regularity: Vec<RegularityClaimV1>,
    ordering: SpectralOrderingV1,
    scope: CompletenessScopeV1,
) -> SpectralProblemSpecV1 {
    problem_spec(
        axes,
        structures,
        regularity,
        SpectralSpaceContextV1::new(
            axes.domain,
            axes.codomain,
            GaugeConventionV1::Unknown,
            ZeroPaddingConventionV1::Unknown,
        ),
        ordering,
        scope,
    )
}

fn structure_claim_with_seed(
    axes: Axes,
    property: StructurePropertyV1,
    support: StructureSupportV1,
    disposition: WitnessDispositionV1,
    tolerance: f64,
    norm: SpectralNormId,
    seed: u8,
) -> StructureClaimV1 {
    let receipt = structure_proposition_receipt(
        axes.subject,
        axes.scalar,
        axes.class,
        axes.scaling,
        axes.domain,
        axes.codomain,
        property,
        support,
        disposition,
        tolerance,
        norm,
    )
    .unwrap();
    StructureClaimV1::new(
        property,
        support,
        disposition,
        tolerance,
        norm,
        admit_receipt(receipt, seed),
    )
}

fn structure_claim(
    axes: Axes,
    property: StructurePropertyV1,
    support: StructureSupportV1,
    disposition: WitnessDispositionV1,
    tolerance: f64,
    seed: u8,
) -> StructureClaimV1 {
    structure_claim_with_seed(
        axes,
        property,
        support,
        disposition,
        tolerance,
        norm_id(200),
        seed,
    )
}

fn regularity_claim(
    axes: Axes,
    class: RegularityClassV1,
    disposition: WitnessDispositionV1,
    seed: u8,
) -> RegularityClaimV1 {
    let receipt = regularity_proposition_receipt(
        axes.subject,
        axes.scalar,
        axes.class,
        axes.scaling,
        axes.domain,
        axes.codomain,
        class,
        disposition,
    )
    .unwrap();
    RegularityClaimV1::new(class, disposition, admit_receipt(receipt, seed))
}

fn truth_witness(
    problem: SpectralProblemId,
    proposition: SpectralTruthPropositionV1,
    seed: u8,
) -> AdmittedSpectralWitnessV1 {
    admit_receipt(
        truth_proposition_receipt(problem, &proposition).unwrap(),
        seed,
    )
}

fn exact_cluster(
    problem: SpectralProblemId,
    id: SpectralClusterIdV1,
    enclosure: SpectralEnclosureV1,
    multiplicity: u32,
    seed: u8,
) -> SpectralClusterV1 {
    let witness = truth_witness(
        problem,
        SpectralTruthPropositionV1::Multiplicity {
            cluster: id,
            enclosure,
            kind: MultiplicityKindV1::Algebraic,
            assertion: MultiplicityAssertionV1::Exact,
            lower: multiplicity,
            upper: Some(multiplicity),
        },
        seed,
    );
    SpectralClusterV1::new(
        id,
        SpectralLocalizationV1::candidate(enclosure),
        MultiplicityClaimV1::Exact {
            value: multiplicity,
            witness,
        },
        MultiplicityClaimV1::Unknown,
        if multiplicity == 1 {
            InternalClusterStateV1::Simple
        } else {
            InternalClusterStateV1::Unknown {
                reason: UnknownSeparationReasonV1::MissingEvidence,
            }
        },
    )
    .unwrap()
}

fn candidate_cluster(id: SpectralClusterIdV1) -> SpectralClusterV1 {
    SpectralClusterV1::new(
        id,
        SpectralLocalizationV1::candidate(SpectralEnclosureV1::Real(
            FiniteIntervalV1::new(0.0, 1.0).unwrap(),
        )),
        MultiplicityClaimV1::Unknown,
        MultiplicityClaimV1::Unknown,
        InternalClusterStateV1::NoClaim,
    )
    .unwrap()
}

fn has_admission_issue(
    report: &SpectralAdmissionReportV1,
    predicate: impl Fn(&SpectralAdmissionIssueV1) -> bool,
) -> bool {
    report.issues().iter().any(predicate)
}

fn has_truth_issue(
    report: &SpectralTruthReportV1,
    predicate: impl Fn(&SpectralTruthErrorV1) -> bool,
) -> bool {
    report.issues().iter().any(predicate)
}

#[test]
fn legacy_and_unknown_schema_versions_fail_closed_before_identity_admission() {
    let axes = standard_axes(0, 2);
    assert_eq!(SPECTRAL_PROBLEM_SCHEMA_VERSION, 2);
    assert_eq!(
        SpectralProblemIdentitySchemaV2::DOMAIN,
        "org.frankensim.fs-spectral.problem-semantic.v2"
    );
    for found in [1, SPECTRAL_PROBLEM_SCHEMA_VERSION + 1] {
        let report = validate_problem(SpectralProblemSpecV1::with_schema_version(
            found,
            axes.subject,
            axes.scalar,
            axes.class,
            StructureProfileV1::new(Vec::new()),
            axes.scaling,
            SpectralSpaceContextV1::new(
                axes.domain,
                axes.codomain,
                GaugeConventionV1::Unknown,
                ZeroPaddingConventionV1::Unknown,
            ),
            RegularityProfileV1::new(Vec::new()),
            SpectralOrderingV1::SetValued,
            CompletenessScopeV1::CandidateOnly,
        ))
        .unwrap_err();
        assert_eq!(
            report.issues(),
            &[SpectralAdmissionIssueV1::UnsupportedSchemaVersion {
                found,
                supported: SPECTRAL_PROBLEM_SCHEMA_VERSION,
            }]
        );
    }
}

#[test]
fn authority_typestate_checks_subject_preimage_anchor_verifier_and_policy() {
    let axes = standard_axes(1, 4);
    let receipt = regularity_proposition_receipt(
        axes.subject,
        axes.scalar,
        axes.class,
        axes.scaling,
        axes.domain,
        axes.codomain,
        RegularityClassV1::FiniteDimensional,
        WitnessDispositionV1::Witnessed,
    )
    .unwrap();
    let authority = exact_authority(receipt, 1);

    let wrong_preimage = ExactAuthority {
        preimage: ContentId::of_bytes(b"wrong-preimage"),
        ..authority
    };
    assert_eq!(
        present(receipt, authority).verify(&wrong_preimage),
        Err(AuthorityError::Mismatch)
    );
    let wrong_anchor = ExactAuthority {
        anchor: ExternalAnchorRef::presented(ContentId::of_bytes(b"wrong-anchor")),
        ..authority
    };
    assert_eq!(
        present(receipt, authority).verify(&wrong_anchor),
        Err(AuthorityError::Mismatch)
    );
    let wrong_verifier = ExactAuthority {
        verifier: spectral_verifier_receipt(b"wrong-verifier").unwrap(),
        ..authority
    };
    assert_eq!(
        present(receipt, authority).verify(&wrong_verifier),
        Err(AuthorityError::Mismatch)
    );
    let wrong_policy = ExactAuthority {
        policy: spectral_authority_policy_receipt(b"wrong-policy").unwrap(),
        ..authority
    };
    assert_eq!(
        present(receipt, authority).verify(&wrong_policy),
        Err(AuthorityError::Mismatch)
    );

    let verified = present(receipt, authority).verify(&authority).unwrap();
    assert_eq!(verified.admit(&wrong_policy), Err(AuthorityError::Mismatch));
    let witness = admit_receipt(receipt, 1);
    assert!(witness.matches_receipt(receipt));
    assert_eq!(witness.audit().trust(), TrustState::Admitted);
    assert_eq!(
        witness.audit().no_claim(),
        NoClaimState::ScientificCorrectnessNotProven
    );
    let promotion = witness.promotion_audit();
    assert_eq!(
        promotion.verifier_domain,
        SpectralVerifierIdentitySchemaV1::DOMAIN
    );
    assert_eq!(
        promotion.key_policy_domain,
        SpectralAuthorityPolicySchemaV1::DOMAIN
    );
    assert_eq!(
        promotion.verifier_observation,
        ObservedIdentity::from_receipt(authority.verifier).bytes()
    );
    assert_eq!(
        promotion.key_policy_observation,
        ObservedIdentity::from_receipt(authority.policy).bytes()
    );
    assert_eq!(promotion.context, SPECTRAL_PROMOTION_CONTEXT_V1);
}

#[test]
fn configured_root_refuses_permit_all_bindings_outside_its_configuration() {
    let axes = standard_axes(41, 4);
    let receipt = regularity_proposition_receipt(
        axes.subject,
        axes.scalar,
        axes.class,
        axes.scaling,
        axes.domain,
        axes.codomain,
        RegularityClassV1::FiniteDimensional,
        WitnessDispositionV1::Witnessed,
    )
    .unwrap();
    let trusted_verifier =
        spectral_verifier_receipt(b"fs-spectral-test-exact-verifier-v1").unwrap();
    let trusted_policy =
        spectral_authority_policy_receipt(b"fs-spectral-test-admission-policy-v1").unwrap();
    let rogue_verifier = spectral_verifier_receipt(b"foreign-permit-all-verifier").unwrap();
    let rogue_policy = spectral_authority_policy_receipt(b"foreign-permit-all-policy").unwrap();
    let anchor = ExternalAnchorRef::presented(ContentId::of_bytes(b"foreign-permit-all"));
    let root = spectral_promotion_trust_root(trusted_verifier, trusted_policy).unwrap();

    let rogue = AuthorityRef::present(receipt, anchor, rogue_verifier.id(), rogue_policy.id())
        .verify(&PermitAll)
        .unwrap()
        .admit(&PermitAll)
        .unwrap();
    assert_eq!(
        root.admit_for_promotion(
            &rogue,
            ObservedIdentity::from_receipt(rogue_verifier).bytes(),
            ObservedIdentity::from_receipt(rogue_policy).bytes(),
        )
        .unwrap_err(),
        PromotionRefusal::ForeignVerifier
    );

    let trusted_verifier_rogue_policy =
        AuthorityRef::present(receipt, anchor, trusted_verifier.id(), rogue_policy.id())
            .verify(&PermitAll)
            .unwrap()
            .admit(&PermitAll)
            .unwrap();
    assert_eq!(
        root.admit_for_promotion(
            &trusted_verifier_rogue_policy,
            ObservedIdentity::from_receipt(trusted_verifier).bytes(),
            ObservedIdentity::from_receipt(rogue_policy).bytes(),
        )
        .unwrap_err(),
        PromotionRefusal::ForeignKeyPolicy
    );

    let trusted_ids =
        AuthorityRef::present(receipt, anchor, trusted_verifier.id(), trusted_policy.id())
            .verify(&PermitAll)
            .unwrap()
            .admit(&PermitAll)
            .unwrap();
    let forged = ByteObservation::new(ContentId::of_bytes(b"different-canonical-bytes"), 999);
    let verifier_refusal = root
        .admit_for_promotion(
            &trusted_ids,
            forged,
            ObservedIdentity::from_receipt(trusted_policy).bytes(),
        )
        .unwrap_err();
    let PromotionRefusal::VerifierObservationMismatch {
        configured,
        presented,
    } = verifier_refusal
    else {
        panic!("expected verifier observation mismatch, got {verifier_refusal:?}");
    };
    assert_eq!(
        configured,
        ObservedIdentity::from_receipt(trusted_verifier).bytes()
    );
    assert_eq!(presented, forged);

    let policy_refusal = root
        .admit_for_promotion(
            &trusted_ids,
            ObservedIdentity::from_receipt(trusted_verifier).bytes(),
            forged,
        )
        .unwrap_err();
    let PromotionRefusal::KeyPolicyObservationMismatch {
        configured,
        presented,
    } = policy_refusal
    else {
        panic!("expected key-policy observation mismatch, got {policy_refusal:?}");
    };
    assert_eq!(
        configured,
        ObservedIdentity::from_receipt(trusted_policy).bytes()
    );
    assert_eq!(presented, forged);
}

#[test]
fn favorable_witness_pairing_refuses_every_mismatched_promotion_axis() {
    let axes = standard_axes(42, 4);
    let receipt = regularity_proposition_receipt(
        axes.subject,
        axes.scalar,
        axes.class,
        axes.scaling,
        axes.domain,
        axes.codomain,
        RegularityClassV1::FiniteDimensional,
        WitnessDispositionV1::Witnessed,
    )
    .unwrap();
    let authority = exact_authority(receipt, 42);
    let (admitted, promotion, pinned) = policy_relative_and_promotion(receipt, authority);
    AdmittedSpectralWitnessV1::from_authority(&admitted, promotion, pinned).unwrap();

    let other_receipt = regularity_proposition_receipt(
        axes.subject,
        axes.scalar,
        axes.class,
        axes.scaling,
        axes.domain,
        axes.codomain,
        RegularityClassV1::RegularPencil,
        WitnessDispositionV1::Witnessed,
    )
    .unwrap();
    let (_, other_subject, _) =
        policy_relative_and_promotion(other_receipt, exact_authority(other_receipt, 42));
    assert_eq!(
        AdmittedSpectralWitnessV1::from_authority(&admitted, other_subject, pinned),
        Err(SpectralPromotionBindingErrorV1::Subject)
    );

    let (_, other_anchor, _) = policy_relative_and_promotion(receipt, exact_authority(receipt, 43));
    assert_eq!(
        AdmittedSpectralWitnessV1::from_authority(&admitted, other_anchor, pinned),
        Err(SpectralPromotionBindingErrorV1::Anchor)
    );

    let foreign_verifier = ExactAuthority {
        verifier: spectral_verifier_receipt(b"other-exact-verifier").unwrap(),
        ..authority
    };
    let (_, other_verifier, _) = policy_relative_and_promotion(receipt, foreign_verifier);
    assert_eq!(
        AdmittedSpectralWitnessV1::from_authority(&admitted, other_verifier, pinned),
        Err(SpectralPromotionBindingErrorV1::Verifier)
    );

    let foreign_policy = ExactAuthority {
        policy: spectral_authority_policy_receipt(b"other-exact-policy").unwrap(),
        ..authority
    };
    let (_, other_policy, _) = policy_relative_and_promotion(receipt, foreign_policy);
    assert_eq!(
        AdmittedSpectralWitnessV1::from_authority(&admitted, other_policy, pinned),
        Err(SpectralPromotionBindingErrorV1::KeyPolicy)
    );

    let wrong_context_root = SpectralPromotionTrustRootV1::configure(
        ObservedIdentity::from_receipt(authority.verifier),
        ObservedIdentity::from_receipt(authority.policy),
        "org.frankensim.fs-spectral.wrong-promotion-context",
    )
    .unwrap();
    let wrong_context = wrong_context_root
        .admit_for_promotion(
            &admitted,
            ObservedIdentity::from_receipt(authority.verifier).bytes(),
            ObservedIdentity::from_receipt(authority.policy).bytes(),
        )
        .unwrap();
    assert_eq!(
        AdmittedSpectralWitnessV1::from_authority(&admitted, wrong_context, pinned),
        Err(SpectralPromotionBindingErrorV1::Context)
    );
}

#[test]
fn proposition_relabeling_and_problem_replay_fail_closed() {
    let axes = standard_axes(2, 4);
    let self_adjoint = structure_claim(
        axes,
        StructurePropertyV1::SelfAdjoint,
        StructureSupportV1::InnerProduct(axes.domain.id()),
        WitnessDispositionV1::Witnessed,
        0.0,
        2,
    );
    let relabeled = StructureClaimV1::new(
        StructurePropertyV1::Normal,
        self_adjoint.support(),
        self_adjoint.disposition(),
        self_adjoint.tolerance(),
        self_adjoint.norm(),
        *self_adjoint.witness(),
    );
    let report = validate_problem(default_spec(
        axes,
        vec![relabeled],
        Vec::new(),
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap_err();
    assert!(has_admission_issue(&report, |issue| matches!(
        issue,
        SpectralAdmissionIssueV1::WitnessPropositionMismatch { .. }
    )));

    let changed_subject = Axes {
        subject: subject(99),
        ..axes
    };
    let report = validate_problem(default_spec(
        changed_subject,
        vec![self_adjoint],
        Vec::new(),
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap_err();
    assert!(has_admission_issue(&report, |issue| matches!(
        issue,
        SpectralAdmissionIssueV1::WitnessPropositionMismatch { .. }
    )));

    let rebound_scaling = SpectralScalingContextV1::new(
        axes.scaling.id(),
        axes.scaling.spectral_dims(),
        2.0,
        axes.scaling.left_map(),
        axes.scaling.right_map(),
        axes.scaling.operator_map(),
        axes.scaling.inverse_map(),
    );
    let report = validate_problem(default_spec(
        Axes {
            scaling: rebound_scaling,
            ..axes
        },
        vec![self_adjoint],
        Vec::new(),
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap_err();
    assert!(
        has_admission_issue(&report, |issue| matches!(
            issue,
            SpectralAdmissionIssueV1::WitnessPropositionMismatch { .. }
        )),
        "same scaling ID rebound to changed scale escaped: {:#?}",
        report.issues()
    );
}

#[test]
fn specialized_methods_require_exact_structure_and_typed_support() {
    let axes = standard_axes(3, 4);
    let finite = regularity_claim(
        axes,
        RegularityClassV1::FiniteDimensional,
        WitnessDispositionV1::Witnessed,
        3,
    );
    let approximate = structure_claim(
        axes,
        StructurePropertyV1::SelfAdjoint,
        StructureSupportV1::InnerProduct(axes.domain.id()),
        WitnessDispositionV1::Witnessed,
        1.0e-8,
        4,
    );
    let problem = validate_problem(default_spec(
        axes,
        vec![approximate],
        vec![finite],
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap();
    let report =
        assess_method_class(&problem, SpectralMethodClassV1::SelfAdjointLanczos).unwrap_err();
    assert_eq!(
        report.issues(),
        &[SpectralAdmissionIssueV1::ExactStructureWitnessRequired {
            method: SpectralMethodClassV1::SelfAdjointLanczos,
            property: StructurePropertyV1::SelfAdjoint,
        }]
    );

    let exact = structure_claim(
        axes,
        StructurePropertyV1::SelfAdjoint,
        StructureSupportV1::InnerProduct(axes.domain.id()),
        WitnessDispositionV1::Witnessed,
        0.0,
        5,
    );
    let problem = validate_problem(default_spec(
        axes,
        vec![exact],
        vec![finite],
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap();
    assert!(assess_method_class(&problem, SpectralMethodClassV1::SelfAdjointLanczos).is_ok());

    let contradicted = structure_claim(
        axes,
        StructurePropertyV1::SelfAdjoint,
        StructureSupportV1::InnerProduct(axes.domain.id()),
        WitnessDispositionV1::Contradicted,
        0.0,
        6,
    );
    let problem = validate_problem(default_spec(
        axes,
        vec![contradicted],
        vec![finite],
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap();
    let report =
        assess_method_class(&problem, SpectralMethodClassV1::SelfAdjointLanczos).unwrap_err();
    assert!(has_admission_issue(&report, |issue| matches!(
        issue,
        SpectralAdmissionIssueV1::ContradictedMethodObligation {
            property: StructurePropertyV1::SelfAdjoint,
            ..
        }
    )));

    let hamiltonian = structure_claim(
        axes,
        StructurePropertyV1::Hamiltonian,
        StructureSupportV1::SymplecticForm(form_id(7)),
        WitnessDispositionV1::Witnessed,
        0.0,
        7,
    );
    let contradicted_other_form = structure_claim(
        axes,
        StructurePropertyV1::Hamiltonian,
        StructureSupportV1::SymplecticForm(form_id(8)),
        WitnessDispositionV1::Contradicted,
        0.0,
        8,
    );
    let problem = validate_problem(default_spec(
        axes,
        vec![hamiltonian, contradicted_other_form],
        vec![finite],
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap();
    let admitted = assess_method_class(
        &problem,
        SpectralMethodClassV1::HamiltonianStructurePreserving,
    );
    assert_eq!(
        admitted.unwrap().selected_support(),
        Some(StructureSupportV1::SymplecticForm(form_id(7)))
    );
    assert!(assess_method_class(&problem, SpectralMethodClassV1::SelfAdjointLanczos).is_err());
}

#[test]
fn generalized_symplectic_and_krein_method_obligations_are_live() {
    let generalized_class = SpectralProblemClassV1::new(
        SpectralRepresentationV1::GeneralizedPencil,
        DescriptorRoleV1::Ordinary,
        SpectralOperatorOriginV1::Direct,
    );
    let generalized = axes(
        70,
        generalized_class,
        SpectralScalarFieldV1::Complex,
        Dims::NONE,
        4,
    );
    let hdp = structure_claim(
        generalized,
        StructurePropertyV1::HermitianDefinitePencil,
        StructureSupportV1::InnerProduct(generalized.domain.id()),
        WitnessDispositionV1::Witnessed,
        0.0,
        70,
    );
    let problem = validate_problem(default_spec(
        generalized,
        vec![hdp],
        Vec::new(),
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap();
    assert!(
        assess_method_class(
            &problem,
            SpectralMethodClassV1::GeneralizedSelfAdjointLanczos,
        )
        .is_ok()
    );

    let symplectic_axes = standard_axes(72, 4);
    let symplectic_finite = regularity_claim(
        symplectic_axes,
        RegularityClassV1::FiniteDimensional,
        WitnessDispositionV1::Witnessed,
        72,
    );
    let symplectic = structure_claim(
        symplectic_axes,
        StructurePropertyV1::Symplectic,
        StructureSupportV1::SymplecticForm(form_id(72)),
        WitnessDispositionV1::Witnessed,
        0.0,
        73,
    );
    let problem = validate_problem(default_spec(
        symplectic_axes,
        vec![symplectic],
        vec![symplectic_finite],
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap();
    assert!(
        assess_method_class(
            &problem,
            SpectralMethodClassV1::SymplecticStructurePreserving,
        )
        .is_ok()
    );
    let j_self_adjoint = structure_claim(
        symplectic_axes,
        StructurePropertyV1::JSelfAdjoint,
        StructureSupportV1::KreinForm(form_id(73)),
        WitnessDispositionV1::Witnessed,
        0.0,
        74,
    );
    let positive_problem = validate_problem(default_spec(
        symplectic_axes,
        vec![j_self_adjoint],
        vec![symplectic_finite],
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap();
    let report = assess_method_class(&positive_problem, SpectralMethodClassV1::KreinJOrthogonal)
        .unwrap_err();
    assert!(has_admission_issue(&report, |issue| matches!(
        issue,
        SpectralAdmissionIssueV1::IndefiniteMetricRequired {
            method: SpectralMethodClassV1::KreinJOrthogonal,
        }
    )));

    let metric_id = SpectralMetricId::from_bytes([74; 32]);
    let metric_witness = admit_receipt(
        metric_proposition_receipt(
            metric_id,
            4,
            MetricDefinitenessPropositionV1::Indefinite {
                positive: 2,
                negative: 2,
            },
        )
        .unwrap(),
        75,
    );
    let indefinite_metric = SpectralMetricV1::new(
        metric_id,
        4,
        MetricDefinitenessV1::Indefinite {
            positive: 2,
            negative: 2,
            witness: metric_witness,
        },
    );
    let indefinite_axes = Axes {
        domain: indefinite_metric,
        codomain: indefinite_metric,
        ..symplectic_axes
    };
    let indefinite_finite = regularity_claim(
        indefinite_axes,
        RegularityClassV1::FiniteDimensional,
        WitnessDispositionV1::Witnessed,
        76,
    );
    let j_self_adjoint = structure_claim(
        indefinite_axes,
        StructurePropertyV1::JSelfAdjoint,
        StructureSupportV1::KreinForm(form_id(74)),
        WitnessDispositionV1::Witnessed,
        0.0,
        77,
    );
    let metric_self_adjoint = structure_claim(
        indefinite_axes,
        StructurePropertyV1::SelfAdjoint,
        StructureSupportV1::InnerProduct(indefinite_metric.id()),
        WitnessDispositionV1::Witnessed,
        0.0,
        78,
    );
    let problem = validate_problem(default_spec(
        indefinite_axes,
        vec![j_self_adjoint, metric_self_adjoint],
        vec![indefinite_finite],
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap();
    assert!(assess_method_class(&problem, SpectralMethodClassV1::KreinJOrthogonal).is_ok());
    let report =
        assess_method_class(&problem, SpectralMethodClassV1::SelfAdjointLanczos).unwrap_err();
    assert!(has_admission_issue(&report, |issue| matches!(
        issue,
        SpectralAdmissionIssueV1::PositiveDefiniteMetricRequired {
            method: SpectralMethodClassV1::SelfAdjointLanczos,
        }
    )));

    let nonnormal = structure_claim(
        symplectic_axes,
        StructurePropertyV1::Nonnormal,
        StructureSupportV1::InnerProduct(symplectic_axes.domain.id()),
        WitnessDispositionV1::Witnessed,
        0.0,
        79,
    );
    let problem = validate_problem(default_spec(
        symplectic_axes,
        vec![nonnormal],
        vec![symplectic_finite],
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap();
    assert!(assess_method_class(&problem, SpectralMethodClassV1::SelfAdjointLanczos).is_err());
}

#[test]
fn real_scalar_and_conjugate_pairs_do_not_forge_real_spectrum_ordering() {
    let axes = standard_axes(8, 4);
    let bare = validate_problem(default_spec(
        axes,
        Vec::new(),
        Vec::new(),
        SpectralOrderingV1::RealAscending,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap_err();
    assert!(
        bare.issues()
            .contains(&SpectralAdmissionIssueV1::OrderingUnavailable)
    );

    let conjugate = structure_claim(
        axes,
        StructurePropertyV1::RealConjugatePairs,
        StructureSupportV1::Conjugation(form_id(8)),
        WitnessDispositionV1::Witnessed,
        0.0,
        8,
    );
    assert!(
        validate_problem(default_spec(
            axes,
            vec![conjugate],
            Vec::new(),
            SpectralOrderingV1::RealAscending,
            CompletenessScopeV1::CandidateOnly,
        ))
        .is_err()
    );

    let approximate = structure_claim(
        axes,
        StructurePropertyV1::RealSpectrum,
        StructureSupportV1::FormFree,
        WitnessDispositionV1::Witnessed,
        1.0e-12,
        9,
    );
    assert!(
        validate_problem(default_spec(
            axes,
            vec![approximate],
            Vec::new(),
            SpectralOrderingV1::RealAscending,
            CompletenessScopeV1::CandidateOnly,
        ))
        .is_err()
    );

    let exact = structure_claim(
        axes,
        StructurePropertyV1::RealSpectrum,
        StructureSupportV1::FormFree,
        WitnessDispositionV1::Witnessed,
        0.0,
        10,
    );
    assert!(
        validate_problem(default_spec(
            axes,
            vec![exact],
            Vec::new(),
            SpectralOrderingV1::RealAscending,
            CompletenessScopeV1::CandidateOnly,
        ))
        .is_ok()
    );

    let self_adjoint = structure_claim(
        axes,
        StructurePropertyV1::SelfAdjoint,
        StructureSupportV1::InnerProduct(axes.domain.id()),
        WitnessDispositionV1::Witnessed,
        0.0,
        11,
    );
    let theorem_closed = validate_problem(default_spec(
        axes,
        vec![self_adjoint],
        Vec::new(),
        SpectralOrderingV1::RealAscending,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap();
    assert!(theorem_closed.requires_real_spectrum_truth());
}

#[test]
fn descriptor_polynomial_and_floquet_obligations_cross_route_exactly() {
    let descriptor_class = SpectralProblemClassV1::new(
        SpectralRepresentationV1::MatrixPolynomial { grade: 2 },
        DescriptorRoleV1::Descriptor {
            infinity_policy: InfiniteEigenvaluePolicyV1::ExcludeWithCount,
        },
        SpectralOperatorOriginV1::Direct,
    );
    let descriptor_axes = axes(
        11,
        descriptor_class,
        SpectralScalarFieldV1::Complex,
        Dims::NONE,
        3,
    );
    let descriptor_problem = validate_problem(default_spec(
        descriptor_axes,
        Vec::new(),
        vec![
            regularity_claim(
                descriptor_axes,
                RegularityClassV1::RegularDescriptor,
                WitnessDispositionV1::Witnessed,
                11,
            ),
            regularity_claim(
                descriptor_axes,
                RegularityClassV1::RegularPolynomial { grade: 2 },
                WitnessDispositionV1::Witnessed,
                12,
            ),
        ],
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::FullFinite {
            algebraic_cardinality: 6,
            infinity_policy: InfiniteEigenvaluePolicyV1::ExcludeWithCount,
        },
    ))
    .unwrap();
    assert!(
        assess_method_class(&descriptor_problem, SpectralMethodClassV1::DescriptorPencil).is_ok()
    );
    assert!(
        assess_method_class(&descriptor_problem, SpectralMethodClassV1::PolynomialKrylov).is_ok()
    );
    assert!(
        assess_method_class(&descriptor_problem, SpectralMethodClassV1::GeneralArnoldi).is_err()
    );

    let floquet_class = SpectralProblemClassV1::new(
        SpectralRepresentationV1::StandardLinear,
        DescriptorRoleV1::Ordinary,
        SpectralOperatorOriginV1::MonodromyFloquet {
            period: Time::new(0.02),
            parameter: FloquetParameterV1::Multiplier,
            branch: FloquetBranchConventionV1::MultipliersOnly,
        },
    );
    let floquet_axes = axes(
        13,
        floquet_class,
        SpectralScalarFieldV1::Complex,
        Dims::NONE,
        4,
    );
    let floquet_problem = validate_problem(default_spec(
        floquet_axes,
        Vec::new(),
        vec![regularity_claim(
            floquet_axes,
            RegularityClassV1::WellPosedMonodromy,
            WitnessDispositionV1::Witnessed,
            13,
        )],
        SpectralOrderingV1::MagnitudeAscending {
            tie_break: ComplexTieBreakV1::RealThenImagThenLineage,
        },
        CompletenessScopeV1::Partial { requested: 2 },
    ))
    .unwrap();
    assert!(assess_method_class(&floquet_problem, SpectralMethodClassV1::MonodromyArnoldi).is_ok());
    assert!(assess_method_class(&floquet_problem, SpectralMethodClassV1::GeneralArnoldi).is_err());

    let invalid_class = SpectralProblemClassV1::new(
        SpectralRepresentationV1::StandardLinear,
        DescriptorRoleV1::Ordinary,
        SpectralOperatorOriginV1::MonodromyFloquet {
            period: Time::new(0.02),
            parameter: FloquetParameterV1::Multiplier,
            branch: FloquetBranchConventionV1::ContinuousFrom {
                continuation: continuation_id(13),
                anchor_phase: Angle::new(0.0),
            },
        },
    );
    let invalid_axes = axes(
        14,
        invalid_class,
        SpectralScalarFieldV1::Complex,
        Dims::NONE,
        4,
    );
    let report = validate_problem(default_spec(
        invalid_axes,
        Vec::new(),
        Vec::new(),
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap_err();
    assert!(
        report
            .issues()
            .contains(&SpectralAdmissionIssueV1::FloquetSemanticMismatch)
    );
}

#[test]
fn operator_function_no_claim_is_classifiable_but_not_executable() {
    let class = SpectralProblemClassV1::new(
        SpectralRepresentationV1::StandardLinear,
        DescriptorRoleV1::Ordinary,
        SpectralOperatorOriginV1::AnalyticOperatorFunction {
            function: function_id(15),
            branch_policy: OperatorFunctionBranchPolicyV1::NoClaim,
        },
    );
    let axes = axes(15, class, SpectralScalarFieldV1::Complex, Dims::NONE, 4);
    let problem = validate_problem(default_spec(
        axes,
        Vec::new(),
        vec![regularity_claim(
            axes,
            RegularityClassV1::AnalyticOperatorFunction,
            WitnessDispositionV1::Witnessed,
            15,
        )],
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap();
    let report =
        assess_method_class(&problem, SpectralMethodClassV1::OperatorFunctionKrylov).unwrap_err();
    assert!(
        report
            .issues()
            .contains(&SpectralAdmissionIssueV1::MethodOriginMismatch {
                method: SpectralMethodClassV1::OperatorFunctionKrylov,
            })
    );
}

#[test]
fn metric_gauge_and_zero_padding_evidence_cannot_be_rebound() {
    let axes = standard_axes(16, 4);
    let metric_receipt = metric_proposition_receipt(
        axes.domain.id(),
        4,
        MetricDefinitenessPropositionV1::PositiveDefinite {
            lower: 0.5,
            upper: 2.0,
        },
    )
    .unwrap();
    let positive_metric = SpectralMetricV1::new(
        axes.domain.id(),
        4,
        MetricDefinitenessV1::PositiveDefinite {
            lower: 0.5,
            upper: 2.0,
            witness: admit_receipt(metric_receipt, 16),
        },
    );
    let conflicting_metric =
        SpectralMetricV1::new(axes.domain.id(), 4, MetricDefinitenessV1::Euclidean);
    let report = validate_problem(problem_spec(
        axes,
        Vec::new(),
        Vec::new(),
        SpectralSpaceContextV1::new(
            positive_metric,
            conflicting_metric,
            GaugeConventionV1::Unknown,
            ZeroPaddingConventionV1::Unknown,
        ),
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap_err();
    assert!(has_admission_issue(&report, |issue| matches!(
        issue,
        SpectralAdmissionIssueV1::MetricIdentityConflict { .. }
    )));

    let fabricated_euclidean = SpectralMetricV1::new(
        SpectralMetricId::from_bytes([0xEF; 32]),
        4,
        MetricDefinitenessV1::Euclidean,
    );
    let report = validate_problem(problem_spec(
        axes,
        Vec::new(),
        Vec::new(),
        SpectralSpaceContextV1::new(
            fabricated_euclidean,
            fabricated_euclidean,
            GaugeConventionV1::Unknown,
            ZeroPaddingConventionV1::Unknown,
        ),
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap_err();
    assert!(
        has_admission_issue(&report, |issue| matches!(
            issue,
            SpectralAdmissionIssueV1::InvalidMetric { .. }
        )),
        "caller-fabricated Euclidean metric was not refused: {:#?}",
        report.issues()
    );

    let fixed_gauge = gauge_artifact_id(17);
    let gauge_receipt = gauge_proposition_receipt(
        axes.subject,
        axes.scalar,
        axes.class,
        axes.scaling,
        axes.domain,
        axes.codomain,
        GaugePropositionV1::Fixed {
            nullity: 1,
            gauge: fixed_gauge,
        },
    )
    .unwrap();
    let gauge_witness = admit_receipt(gauge_receipt, 17);
    let valid_spaces = SpectralSpaceContextV1::new(
        axes.domain,
        axes.codomain,
        GaugeConventionV1::Fixed {
            nullity: 1,
            gauge: fixed_gauge,
            witness: gauge_witness,
        },
        ZeroPaddingConventionV1::Unknown,
    );
    let fixed_problem = validate_problem(problem_spec(
        axes,
        Vec::new(),
        Vec::new(),
        valid_spaces,
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap();
    assert_eq!(fixed_problem.known_algebraic_cardinality(), Some(4));

    let report = validate_problem(problem_spec(
        axes,
        Vec::new(),
        Vec::new(),
        SpectralSpaceContextV1::new(
            axes.domain,
            axes.codomain,
            GaugeConventionV1::Fixed {
                nullity: 1,
                gauge: gauge_artifact_id(18),
                witness: gauge_witness,
            },
            ZeroPaddingConventionV1::Unknown,
        ),
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap_err();
    assert!(has_admission_issue(&report, |issue| matches!(
        issue,
        SpectralAdmissionIssueV1::WitnessPropositionMismatch { .. }
    )));

    let quotient_a = quotient_map_id(18);
    let rebound_spaces = SpectralSpaceContextV1::new(
        axes.domain,
        axes.codomain,
        GaugeConventionV1::Quotiented {
            nullity: 1,
            quotient: quotient_a,
            witness: gauge_witness,
        },
        ZeroPaddingConventionV1::Unknown,
    );
    let report = validate_problem(problem_spec(
        axes,
        Vec::new(),
        Vec::new(),
        rebound_spaces,
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap_err();
    assert!(has_admission_issue(&report, |issue| matches!(
        issue,
        SpectralAdmissionIssueV1::WitnessPropositionMismatch { .. }
    )));

    let quotient_witness_a = admit_receipt(
        gauge_proposition_receipt(
            axes.subject,
            axes.scalar,
            axes.class,
            axes.scaling,
            axes.domain,
            axes.codomain,
            GaugePropositionV1::Quotiented {
                nullity: 1,
                quotient: quotient_a,
            },
        )
        .unwrap(),
        18,
    );
    let quotient_problem = validate_problem(problem_spec(
        axes,
        Vec::new(),
        Vec::new(),
        SpectralSpaceContextV1::new(
            axes.domain,
            axes.codomain,
            GaugeConventionV1::Quotiented {
                nullity: 1,
                quotient: quotient_a,
                witness: quotient_witness_a,
            },
            ZeroPaddingConventionV1::Unknown,
        ),
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::FullFinite {
            algebraic_cardinality: 4,
            infinity_policy: InfiniteEigenvaluePolicyV1::NoClaim,
        },
    ))
    .unwrap();
    assert_eq!(quotient_problem.known_algebraic_cardinality(), Some(4));

    let quotient_b = quotient_map_id(19);
    let report = validate_problem(problem_spec(
        axes,
        Vec::new(),
        Vec::new(),
        SpectralSpaceContextV1::new(
            axes.domain,
            axes.codomain,
            GaugeConventionV1::Quotiented {
                nullity: 1,
                quotient: quotient_b,
                witness: quotient_witness_a,
            },
            ZeroPaddingConventionV1::Unknown,
        ),
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::FullFinite {
            algebraic_cardinality: 4,
            infinity_policy: InfiniteEigenvaluePolicyV1::NoClaim,
        },
    ))
    .unwrap_err();
    assert!(has_admission_issue(&report, |issue| matches!(
        issue,
        SpectralAdmissionIssueV1::WitnessPropositionMismatch { .. }
    )));

    let quotient_witness_b = admit_receipt(
        gauge_proposition_receipt(
            axes.subject,
            axes.scalar,
            axes.class,
            axes.scaling,
            axes.domain,
            axes.codomain,
            GaugePropositionV1::Quotiented {
                nullity: 1,
                quotient: quotient_b,
            },
        )
        .unwrap(),
        19,
    );
    let induced_structure = structure_claim(
        axes,
        StructurePropertyV1::RealSpectrum,
        StructureSupportV1::FormFree,
        WitnessDispositionV1::Witnessed,
        0.0,
        20,
    );
    let induced_regularity = regularity_claim(
        axes,
        RegularityClassV1::FiniteDimensional,
        WitnessDispositionV1::Witnessed,
        20,
    );
    let quotient_with_claims = |quotient, witness| {
        validate_problem(problem_spec(
            axes,
            vec![induced_structure],
            vec![induced_regularity],
            SpectralSpaceContextV1::new(
                axes.domain,
                axes.codomain,
                GaugeConventionV1::Quotiented {
                    nullity: 1,
                    quotient,
                    witness,
                },
                ZeroPaddingConventionV1::Unknown,
            ),
            SpectralOrderingV1::SetValued,
            CompletenessScopeV1::CandidateOnly,
        ))
        .unwrap()
    };
    let induced_a = quotient_with_claims(quotient_a, quotient_witness_a);
    let induced_b = quotient_with_claims(quotient_b, quotient_witness_b);
    assert_ne!(induced_a.problem_id(), induced_b.problem_id());

    let quotient_a_padding = admit_receipt(
        zero_padding_proposition_receipt(
            axes.subject,
            axes.scalar,
            axes.class,
            axes.scaling,
            axes.domain,
            axes.codomain,
            GaugeContextV1::Quotiented {
                nullity: 1,
                quotient: quotient_a,
            },
            ZeroPaddingPropositionV1::Omitted { count: 1 },
        )
        .unwrap(),
        20,
    );
    let report = validate_problem(problem_spec(
        axes,
        Vec::new(),
        Vec::new(),
        SpectralSpaceContextV1::new(
            axes.domain,
            axes.codomain,
            GaugeConventionV1::Quotiented {
                nullity: 1,
                quotient: quotient_b,
                witness: quotient_witness_b,
            },
            ZeroPaddingConventionV1::Omitted {
                count: 1,
                witness: quotient_a_padding,
            },
        ),
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap_err();
    assert!(has_admission_issue(&report, |issue| matches!(
        issue,
        SpectralAdmissionIssueV1::WitnessPropositionMismatch { .. }
    )));

    let quotient_problem_b = validate_problem(problem_spec(
        axes,
        Vec::new(),
        Vec::new(),
        SpectralSpaceContextV1::new(
            axes.domain,
            axes.codomain,
            GaugeConventionV1::Quotiented {
                nullity: 1,
                quotient: quotient_b,
                witness: quotient_witness_b,
            },
            ZeroPaddingConventionV1::Unknown,
        ),
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::FullFinite {
            algebraic_cardinality: 4,
            infinity_policy: InfiniteEigenvaluePolicyV1::NoClaim,
        },
    ))
    .unwrap();
    assert_ne!(
        quotient_problem.problem_id(),
        quotient_problem_b.problem_id(),
        "the exact quotient map must participate in semantic identity"
    );

    let large_nullity = 7;
    let large_quotient = quotient_map_id(20);
    let large_quotient_witness = admit_receipt(
        gauge_proposition_receipt(
            axes.subject,
            axes.scalar,
            axes.class,
            axes.scaling,
            axes.domain,
            axes.codomain,
            GaugePropositionV1::Quotiented {
                nullity: large_nullity,
                quotient: large_quotient,
            },
        )
        .unwrap(),
        20,
    );
    let large_padding_witness = admit_receipt(
        zero_padding_proposition_receipt(
            axes.subject,
            axes.scalar,
            axes.class,
            axes.scaling,
            axes.domain,
            axes.codomain,
            GaugeContextV1::Quotiented {
                nullity: large_nullity,
                quotient: large_quotient,
            },
            ZeroPaddingPropositionV1::ExplicitlyPadded {
                count: large_nullity,
            },
        )
        .unwrap(),
        21,
    );
    let large_pre_reduction_nullity = validate_problem(problem_spec(
        axes,
        Vec::new(),
        Vec::new(),
        SpectralSpaceContextV1::new(
            axes.domain,
            axes.codomain,
            GaugeConventionV1::Quotiented {
                nullity: large_nullity,
                quotient: large_quotient,
                witness: large_quotient_witness,
            },
            ZeroPaddingConventionV1::ExplicitlyPadded {
                count: large_nullity,
                witness: large_padding_witness,
            },
        ),
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap();
    assert_eq!(
        large_pre_reduction_nullity.known_algebraic_cardinality(),
        Some(4)
    );
    assert!(assess_gap_semantics(&large_pre_reduction_nullity).is_ok());
}

#[test]
fn inner_product_structure_requires_one_shared_admitted_operator_space() {
    let standard = standard_axes(17, 4);
    let codomain_id = SpectralMetricId::from_bytes([0xA1; 32]);
    let codomain_receipt = metric_proposition_receipt(
        codomain_id,
        4,
        MetricDefinitenessPropositionV1::PositiveDefinite {
            lower: 0.5,
            upper: 2.0,
        },
    )
    .unwrap();
    let cross_space = Axes {
        codomain: SpectralMetricV1::new(
            codomain_id,
            4,
            MetricDefinitenessV1::PositiveDefinite {
                lower: 0.5,
                upper: 2.0,
                witness: admit_receipt(codomain_receipt, 18),
            },
        ),
        ..standard
    };
    let self_adjoint = structure_claim(
        cross_space,
        StructurePropertyV1::SelfAdjoint,
        StructureSupportV1::InnerProduct(cross_space.domain.id()),
        WitnessDispositionV1::Witnessed,
        0.0,
        17,
    );
    let report = validate_problem(default_spec(
        cross_space,
        vec![self_adjoint],
        Vec::new(),
        SpectralOrderingV1::RealAscending,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap_err();
    assert!(
        has_admission_issue(&report, |issue| matches!(
            issue,
            SpectralAdmissionIssueV1::InvalidStructureSupport {
                property: StructurePropertyV1::SelfAdjoint,
                support: StructureSupportV1::InnerProduct(found),
            } if *found == cross_space.domain.id()
        )),
        "a one-sided metric witness must not admit cross-space self-adjointness: {report:?}"
    );
    assert!(
        report
            .issues()
            .contains(&SpectralAdmissionIssueV1::OrderingUnavailable),
        "cross-space self-adjointness must not manufacture real-spectrum ordering authority: {report:?}"
    );

    for (property, seed) in [
        (StructurePropertyV1::Normal, 19),
        (StructurePropertyV1::Nonnormal, 20),
    ] {
        let claim = structure_claim(
            cross_space,
            property,
            StructureSupportV1::InnerProduct(cross_space.domain.id()),
            WitnessDispositionV1::Witnessed,
            0.0,
            seed,
        );
        let report = validate_problem(default_spec(
            cross_space,
            vec![claim],
            Vec::new(),
            SpectralOrderingV1::SetValued,
            CompletenessScopeV1::CandidateOnly,
        ))
        .unwrap_err();
        assert!(
            has_admission_issue(&report, |issue| matches!(
                issue,
                SpectralAdmissionIssueV1::InvalidStructureSupport {
                    property: found,
                    ..
                } if *found == property
            )),
            "{property:?} must be bound to a common domain/codomain inner product: {report:?}"
        );
    }

    let generalized_class = SpectralProblemClassV1::new(
        SpectralRepresentationV1::GeneralizedPencil,
        DescriptorRoleV1::Ordinary,
        SpectralOperatorOriginV1::Direct,
    );
    let generalized = axes(
        21,
        generalized_class,
        SpectralScalarFieldV1::Real,
        Dims::NONE,
        4,
    );
    let pencil_codomain_id = SpectralMetricId::from_bytes([0xA2; 32]);
    let pencil_codomain_receipt = metric_proposition_receipt(
        pencil_codomain_id,
        4,
        MetricDefinitenessPropositionV1::PositiveDefinite {
            lower: 0.25,
            upper: 4.0,
        },
    )
    .unwrap();
    let cross_pencil = Axes {
        codomain: SpectralMetricV1::new(
            pencil_codomain_id,
            4,
            MetricDefinitenessV1::PositiveDefinite {
                lower: 0.25,
                upper: 4.0,
                witness: admit_receipt(pencil_codomain_receipt, 21),
            },
        ),
        ..generalized
    };
    let hermitian_definite = structure_claim(
        cross_pencil,
        StructurePropertyV1::HermitianDefinitePencil,
        StructureSupportV1::InnerProduct(cross_pencil.domain.id()),
        WitnessDispositionV1::Witnessed,
        0.0,
        22,
    );
    let report = validate_problem(default_spec(
        cross_pencil,
        vec![hermitian_definite],
        Vec::new(),
        SpectralOrderingV1::RealAscending,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap_err();
    assert!(
        has_admission_issue(&report, |issue| matches!(
            issue,
            SpectralAdmissionIssueV1::InvalidStructureSupport {
                property: StructurePropertyV1::HermitianDefinitePencil,
                ..
            }
        )),
        "one positive endpoint must not admit a Hermitian-definite pencil: {report:?}"
    );
    assert!(
        report
            .issues()
            .contains(&SpectralAdmissionIssueV1::OrderingUnavailable),
        "the invalid pencil witness must not manufacture real-spectrum ordering authority: {report:?}"
    );
    assert!(
        report.issues().contains(
            &SpectralAdmissionIssueV1::OrdinaryFiniteSpectrumWitnessRequired {
                required: RegularityClassV1::InvertiblePencilWeight,
            },
        ),
        "the invalid pencil witness must not manufacture invertibility or regularity: {report:?}"
    );
}

#[test]
fn unresolved_metric_definiteness_cannot_admit_adjoint_structure() {
    let base = standard_axes(20, 3);
    let unresolved = SpectralMetricV1::new(
        SpectralMetricId::from_bytes([0xA3; 32]),
        3,
        MetricDefinitenessV1::Unknown,
    );
    let axes = Axes {
        domain: unresolved,
        codomain: unresolved,
        ..base
    };
    let claim = structure_claim(
        axes,
        StructurePropertyV1::Normal,
        StructureSupportV1::InnerProduct(unresolved.id()),
        WitnessDispositionV1::Witnessed,
        0.0,
        20,
    );
    let report = validate_problem(default_spec(
        axes,
        vec![claim],
        Vec::new(),
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap_err();
    assert!(
        has_admission_issue(&report, |issue| matches!(
            issue,
            SpectralAdmissionIssueV1::InvalidStructureSupport {
                property: StructurePropertyV1::Normal,
                ..
            }
        )),
        "unknown definiteness cannot establish the nondegenerate inner product needed for an adjoint: {report:?}"
    );
}

#[test]
fn exact_normality_complements_cannot_both_be_contradicted() {
    let axes = standard_axes(23, 3);
    let support = StructureSupportV1::InnerProduct(axes.domain.id());
    let normal_refuted = structure_claim(
        axes,
        StructurePropertyV1::Normal,
        support,
        WitnessDispositionV1::Contradicted,
        0.0,
        23,
    );
    let nonnormal_refuted = structure_claim(
        axes,
        StructurePropertyV1::Nonnormal,
        support,
        WitnessDispositionV1::Contradicted,
        0.0,
        24,
    );
    let report = validate_problem(default_spec(
        axes,
        vec![normal_refuted, nonnormal_refuted],
        Vec::new(),
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap_err();
    assert!(
        report
            .issues()
            .contains(&SpectralAdmissionIssueV1::ComplementaryStructureConflict { support }),
        "logical complements Normal and Nonnormal cannot both be exactly false: {report:?}"
    );

    let nonnormal_witnessed = structure_claim(
        axes,
        StructurePropertyV1::Nonnormal,
        support,
        WitnessDispositionV1::Witnessed,
        0.0,
        25,
    );
    assert!(
        validate_problem(default_spec(
            axes,
            vec![normal_refuted, nonnormal_witnessed],
            Vec::new(),
            SpectralOrderingV1::SetValued,
            CompletenessScopeV1::CandidateOnly,
        ))
        .is_ok(),
        "a witnessed Nonnormal proposition is consistent with refuted Normal"
    );
}

#[test]
fn gap_interpretation_requires_explicit_gauge_and_zero_serialization_semantics() {
    let axes = standard_axes(17, 4);
    let unknown = validate_problem(default_spec(
        axes,
        Vec::new(),
        Vec::new(),
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap();
    let report = assess_gap_semantics(&unknown).unwrap_err();
    assert_eq!(
        report.issues(),
        &[
            SpectralAdmissionIssueV1::GapGaugeConventionRequired,
            SpectralAdmissionIssueV1::GapZeroPaddingConventionRequired,
        ],
        "unknown structural-zero semantics did not fail closed"
    );

    let gauge = admit_receipt(
        gauge_proposition_receipt(
            axes.subject,
            axes.scalar,
            axes.class,
            axes.scaling,
            axes.domain,
            axes.codomain,
            GaugePropositionV1::None,
        )
        .unwrap(),
        17,
    );
    let zero_padding = admit_receipt(
        zero_padding_proposition_receipt(
            axes.subject,
            axes.scalar,
            axes.class,
            axes.scaling,
            axes.domain,
            axes.codomain,
            GaugeContextV1::CertifiedNone,
            ZeroPaddingPropositionV1::NonePresent,
        )
        .unwrap(),
        18,
    );
    let explicit = validate_problem(problem_spec(
        axes,
        Vec::new(),
        Vec::new(),
        SpectralSpaceContextV1::new(
            axes.domain,
            axes.codomain,
            GaugeConventionV1::CertifiedNone { witness: gauge },
            ZeroPaddingConventionV1::CertifiedNonePresent {
                witness: zero_padding,
            },
        ),
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap();
    let token = assess_gap_semantics(&explicit).unwrap();
    assert_eq!(token.problem_id(), explicit.problem_id());
    assert!(matches!(
        token.gauge(),
        GaugeConventionV1::CertifiedNone { .. }
    ));
    assert!(matches!(
        token.zero_padding(),
        ZeroPaddingConventionV1::CertifiedNonePresent { .. }
    ));

    let fixed_gauge = gauge_artifact_id(19);
    let fixed = admit_receipt(
        gauge_proposition_receipt(
            axes.subject,
            axes.scalar,
            axes.class,
            axes.scaling,
            axes.domain,
            axes.codomain,
            GaugePropositionV1::Fixed {
                nullity: 1,
                gauge: fixed_gauge,
            },
        )
        .unwrap(),
        19,
    );
    let mismatched_padding = admit_receipt(
        zero_padding_proposition_receipt(
            axes.subject,
            axes.scalar,
            axes.class,
            axes.scaling,
            axes.domain,
            axes.codomain,
            GaugeContextV1::Fixed {
                nullity: 1,
                gauge: fixed_gauge,
            },
            ZeroPaddingPropositionV1::ExplicitlyPadded { count: 2 },
        )
        .unwrap(),
        20,
    );
    let inconsistent = validate_problem(problem_spec(
        axes,
        Vec::new(),
        Vec::new(),
        SpectralSpaceContextV1::new(
            axes.domain,
            axes.codomain,
            GaugeConventionV1::Fixed {
                nullity: 1,
                gauge: fixed_gauge,
                witness: fixed,
            },
            ZeroPaddingConventionV1::ExplicitlyPadded {
                count: 2,
                witness: mismatched_padding,
            },
        ),
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap();
    assert_eq!(
        assess_gap_semantics(&inconsistent).unwrap_err().issues(),
        &[SpectralAdmissionIssueV1::GapStructuralZeroCountMismatch {
            gauge_nullity: 1,
            declared_zero_count: 2,
        }]
    );
}

#[test]
fn normalization_is_dimension_checked_reversible_and_overflow_safe() {
    let scale = scaling(Dims([0, 0, -1, 0, 0, 0]), 2.0, 18);
    assert_eq!(
        scale
            .normalize(QtyAny::new(6.0, Dims([0, 0, -1, 0, 0, 0])))
            .unwrap(),
        3.0
    );
    let denormalized = scale.denormalize(3.0).unwrap();
    assert_eq!(denormalized.value, 6.0);
    assert_eq!(denormalized.dims, Dims([0, 0, -1, 0, 0, 0]));
    assert!(matches!(
        scale.normalize(QtyAny::new(6.0, Dims::NONE)),
        Err(SpectralAdmissionIssueV1::UnitMismatch { .. })
    ));
    assert!(matches!(
        scale.normalize(QtyAny::new(f64::INFINITY, Dims([0, 0, -1, 0, 0, 0]))),
        Err(SpectralAdmissionIssueV1::NonFinite { .. })
    ));
    let tiny = scaling(Dims::NONE, f64::MIN_POSITIVE, 19);
    assert!(matches!(
        tiny.normalize(QtyAny::new(f64::MAX, Dims::NONE)),
        Err(SpectralAdmissionIssueV1::NonFinite {
            field: AdmissionFieldV1::NormalizedSpectralValue
        })
    ));
    let huge = scaling(Dims::NONE, f64::MAX, 20);
    assert_eq!(
        huge.normalize(QtyAny::new(f64::MIN_POSITIVE, Dims::NONE)),
        Err(SpectralAdmissionIssueV1::Underflow {
            field: AdmissionFieldV1::NormalizedSpectralValue,
        })
    );
    assert_eq!(
        tiny.denormalize(f64::MIN_POSITIVE),
        Err(SpectralAdmissionIssueV1::Underflow {
            field: AdmissionFieldV1::SpectralValue,
        })
    );
}

#[test]
fn admission_resource_caps_fire_before_quadratic_claim_analysis() {
    let axes = standard_axes(20, 4);
    // Promotion-bearing schema-v2 witnesses make this canonical set larger
    // than 64 KiB. The problem-specific field envelope must still admit the
    // promised 256-claim boundary without broadening verifier/policy inputs.
    let at_limit_claims = (0..MAX_STRUCTURE_CLAIMS_V1)
        .map(|index| {
            let byte = u8::try_from(index).unwrap();
            structure_claim_with_seed(
                axes,
                StructurePropertyV1::SelfAdjoint,
                StructureSupportV1::InnerProduct(axes.domain.id()),
                WitnessDispositionV1::Witnessed,
                0.0,
                norm_id(byte),
                byte,
            )
        })
        .collect();
    let at_limit = validate_problem(default_spec(
        axes,
        at_limit_claims,
        Vec::new(),
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap();
    assert_eq!(at_limit.structure_claims().len(), MAX_STRUCTURE_CLAIMS_V1);

    let claim = structure_claim(
        axes,
        StructurePropertyV1::SelfAdjoint,
        StructureSupportV1::InnerProduct(axes.domain.id()),
        WitnessDispositionV1::Witnessed,
        0.0,
        20,
    );

    let over_limit = validate_problem(default_spec(
        axes,
        vec![claim; MAX_STRUCTURE_CLAIMS_V1 + 1],
        Vec::new(),
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap_err();
    assert_eq!(
        over_limit.issues(),
        &[SpectralAdmissionIssueV1::TooManyClaims {
            profile: ClaimProfileV1::Structure,
            found: MAX_STRUCTURE_CLAIMS_V1 + 1,
            limit: MAX_STRUCTURE_CLAIMS_V1,
        }]
    );

    let descriptor_class = SpectralProblemClassV1::new(
        SpectralRepresentationV1::GeneralizedPencil,
        DescriptorRoleV1::Descriptor {
            infinity_policy: InfiniteEigenvaluePolicyV1::IncludeProjective,
        },
        SpectralOperatorOriginV1::MonodromyFloquet {
            period: Time::new(0.02),
            parameter: FloquetParameterV1::Multiplier,
            branch: FloquetBranchConventionV1::MultipliersOnly,
        },
    );
    let descriptor_axes = crate::axes(
        21,
        descriptor_class,
        SpectralScalarFieldV1::Complex,
        Dims::NONE,
        4,
    );
    let regularity_at_limit = vec![
        regularity_claim(
            descriptor_axes,
            RegularityClassV1::FiniteDimensional,
            WitnessDispositionV1::Witnessed,
            21,
        ),
        regularity_claim(
            descriptor_axes,
            RegularityClassV1::RegularPencil,
            WitnessDispositionV1::Witnessed,
            22,
        ),
        regularity_claim(
            descriptor_axes,
            RegularityClassV1::InvertiblePencilWeight,
            WitnessDispositionV1::Witnessed,
            23,
        ),
        regularity_claim(
            descriptor_axes,
            RegularityClassV1::RegularDescriptor,
            WitnessDispositionV1::Witnessed,
            24,
        ),
        regularity_claim(
            descriptor_axes,
            RegularityClassV1::WellPosedMonodromy,
            WitnessDispositionV1::Witnessed,
            25,
        ),
    ];
    assert_eq!(regularity_at_limit.len(), MAX_REGULARITY_CLAIMS_V1);
    let at_limit = validate_problem(default_spec(
        descriptor_axes,
        Vec::new(),
        regularity_at_limit,
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap();
    assert_eq!(at_limit.regularity_claims().len(), MAX_REGULARITY_CLAIMS_V1);

    let duplicate = regularity_claim(
        descriptor_axes,
        RegularityClassV1::FiniteDimensional,
        WitnessDispositionV1::Witnessed,
        26,
    );
    let over_limit = validate_problem(default_spec(
        descriptor_axes,
        Vec::new(),
        vec![duplicate; MAX_REGULARITY_CLAIMS_V1 + 1],
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap_err();
    assert_eq!(
        over_limit.issues(),
        &[SpectralAdmissionIssueV1::TooManyClaims {
            profile: ClaimProfileV1::Regularity,
            found: MAX_REGULARITY_CLAIMS_V1 + 1,
            limit: MAX_REGULARITY_CLAIMS_V1,
        }]
    );
}

#[test]
fn problem_identity_is_permutation_stable_and_semantic_axis_sensitive() {
    let class = SpectralProblemClassV1::new(
        SpectralRepresentationV1::MatrixPolynomial { grade: 2 },
        DescriptorRoleV1::Descriptor {
            infinity_policy: InfiniteEigenvaluePolicyV1::ExcludeWithCount,
        },
        SpectralOperatorOriginV1::Direct,
    );
    let axes = axes(21, class, SpectralScalarFieldV1::Complex, Dims::NONE, 2);
    let gyroscopic = structure_claim(
        axes,
        StructurePropertyV1::Gyroscopic,
        StructureSupportV1::FormFree,
        WitnessDispositionV1::Witnessed,
        0.0,
        21,
    );
    let palindromic = structure_claim(
        axes,
        StructurePropertyV1::Palindromic {
            parity: PalindromicParityV1::Palindromic,
            involution: PolynomialInvolutionV1::ConjugateTranspose,
        },
        StructureSupportV1::FormFree,
        WitnessDispositionV1::Witnessed,
        0.0,
        22,
    );
    let descriptor = regularity_claim(
        axes,
        RegularityClassV1::RegularDescriptor,
        WitnessDispositionV1::Witnessed,
        23,
    );
    let polynomial = regularity_claim(
        axes,
        RegularityClassV1::RegularPolynomial { grade: 2 },
        WitnessDispositionV1::Witnessed,
        24,
    );
    let scope = CompletenessScopeV1::FullFinite {
        algebraic_cardinality: 4,
        infinity_policy: InfiniteEigenvaluePolicyV1::ExcludeWithCount,
    };
    let first = validate_problem(default_spec(
        axes,
        vec![gyroscopic, palindromic],
        vec![descriptor, polynomial],
        SpectralOrderingV1::SetValued,
        scope,
    ))
    .unwrap();
    let permuted = validate_problem(default_spec(
        axes,
        vec![palindromic, gyroscopic],
        vec![polynomial, descriptor],
        SpectralOrderingV1::SetValued,
        scope,
    ))
    .unwrap();
    assert_eq!(first.problem_id(), permuted.problem_id());
    assert_eq!(first.identity_receipt(), permuted.identity_receipt());
    assert_eq!(
        first.problem_id().as_bytes(),
        first.identity_receipt().id().as_bytes(),
        "the narrow problem ID and retained producer receipt must name the same canonical observation"
    );
    assert_eq!(first.identity_receipt().field_count(), 9);
    assert_eq!(first.identity_receipt().collection_items(), 4);
    assert!(first.identity_receipt().canonical_bytes() > 0);
    let retained = ObservedIdentity::from_receipt(first.identity_receipt());
    let synthetic_collision = ObservedIdentity::presented(
        first.identity_receipt().id(),
        ByteObservation::new(
            first.identity_receipt().canonical_preimage(),
            first.identity_receipt().canonical_bytes() + 1,
        ),
    );
    assert!(
        matches!(
            adjudicate(retained, synthetic_collision),
            IdentityAdjudication::Refused(_)
        ),
        "retaining the producer observation must permit same-digest/different-bytes refusal"
    );
    assert_eq!(first.spec(), permuted.spec());
    assert_eq!(first.problem_id().to_hex().len(), 64);

    let retargeted = validate_problem(default_spec(
        axes,
        vec![gyroscopic, palindromic],
        vec![descriptor, polynomial],
        SpectralOrderingV1::MagnitudeAscending {
            tie_break: ComplexTieBreakV1::RealThenImagThenLineage,
        },
        scope,
    ))
    .unwrap();
    assert_ne!(first.problem_id(), retargeted.problem_id());
    assert_ne!(first.identity_receipt(), retargeted.identity_receipt());

    let reanchored_gyroscopic = structure_claim_with_seed(
        axes,
        StructurePropertyV1::Gyroscopic,
        StructureSupportV1::FormFree,
        WitnessDispositionV1::Witnessed,
        0.0,
        norm_id(200),
        99,
    );
    let reanchored = validate_problem(default_spec(
        axes,
        vec![reanchored_gyroscopic, palindromic],
        vec![descriptor, polynomial],
        SpectralOrderingV1::SetValued,
        scope,
    ))
    .unwrap();
    assert_ne!(first.problem_id(), reanchored.problem_id());
    assert_ne!(first.identity_receipt(), reanchored.identity_receipt());

    let descriptor_receipt = regularity_proposition_receipt(
        axes.subject,
        axes.scalar,
        axes.class,
        axes.scaling,
        axes.domain,
        axes.codomain,
        RegularityClassV1::RegularDescriptor,
        WitnessDispositionV1::Witnessed,
    )
    .unwrap();
    let authority = exact_authority(descriptor_receipt, 23);
    let admitted = present(descriptor_receipt, authority)
        .verify(&authority)
        .unwrap()
        .admit(&authority)
        .unwrap();
    let changed_observation =
        ByteObservation::new(ContentId::of_bytes(b"independent-verifier-observation"), 37);
    let observation_root = SpectralPromotionTrustRootV1::configure(
        ObservedIdentity::presented(authority.verifier.id(), changed_observation),
        ObservedIdentity::from_receipt(authority.policy),
        SPECTRAL_PROMOTION_CONTEXT_V1,
    )
    .unwrap();
    let reobserved_promotion = observation_root
        .admit_for_promotion(
            &admitted,
            changed_observation,
            ObservedIdentity::from_receipt(authority.policy).bytes(),
        )
        .unwrap();
    let reobserved_descriptor = RegularityClaimV1::new(
        RegularityClassV1::RegularDescriptor,
        WitnessDispositionV1::Witnessed,
        AdmittedSpectralWitnessV1::from_authority(
            &admitted,
            reobserved_promotion,
            observation_root.charter(),
        )
        .unwrap(),
    );
    let reobserved = validate_problem(default_spec(
        axes,
        vec![gyroscopic, palindromic],
        vec![reobserved_descriptor, polynomial],
        SpectralOrderingV1::SetValued,
        scope,
    ))
    .unwrap();
    assert_ne!(first.problem_id(), reobserved.problem_id());

    let changed_length = ByteObservation::new(changed_observation.content_id(), 38);
    let length_root = SpectralPromotionTrustRootV1::configure(
        ObservedIdentity::presented(authority.verifier.id(), changed_length),
        ObservedIdentity::from_receipt(authority.policy),
        SPECTRAL_PROMOTION_CONTEXT_V1,
    )
    .unwrap();
    let length_promotion = length_root
        .admit_for_promotion(
            &admitted,
            changed_length,
            ObservedIdentity::from_receipt(authority.policy).bytes(),
        )
        .unwrap();
    let changed_length_descriptor = RegularityClaimV1::new(
        RegularityClassV1::RegularDescriptor,
        WitnessDispositionV1::Witnessed,
        AdmittedSpectralWitnessV1::from_authority(
            &admitted,
            length_promotion,
            length_root.charter(),
        )
        .unwrap(),
    );
    let changed_length_problem = validate_problem(default_spec(
        axes,
        vec![gyroscopic, palindromic],
        vec![changed_length_descriptor, polynomial],
        SpectralOrderingV1::SetValued,
        scope,
    ))
    .unwrap();
    assert_ne!(reobserved.problem_id(), changed_length_problem.problem_id());

    let changed_policy_observation =
        ByteObservation::new(ContentId::of_bytes(b"independent-policy-observation"), 29);
    let policy_root = SpectralPromotionTrustRootV1::configure(
        ObservedIdentity::from_receipt(authority.verifier),
        ObservedIdentity::presented(authority.policy.id(), changed_policy_observation),
        SPECTRAL_PROMOTION_CONTEXT_V1,
    )
    .unwrap();
    let policy_promotion = policy_root
        .admit_for_promotion(
            &admitted,
            ObservedIdentity::from_receipt(authority.verifier).bytes(),
            changed_policy_observation,
        )
        .unwrap();
    let changed_policy_descriptor = RegularityClaimV1::new(
        RegularityClassV1::RegularDescriptor,
        WitnessDispositionV1::Witnessed,
        AdmittedSpectralWitnessV1::from_authority(
            &admitted,
            policy_promotion,
            policy_root.charter(),
        )
        .unwrap(),
    );
    let changed_policy_problem = validate_problem(default_spec(
        axes,
        vec![gyroscopic, palindromic],
        vec![changed_policy_descriptor, polynomial],
        SpectralOrderingV1::SetValued,
        scope,
    ))
    .unwrap();
    assert_ne!(first.problem_id(), changed_policy_problem.problem_id());
}

#[test]
fn truth_evidence_is_bound_to_problem_result_set_and_proposition_family() {
    let axes = standard_axes(25, 2);
    let problem = validate_problem(default_spec(
        axes,
        Vec::new(),
        Vec::new(),
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap();
    let clusters = vec![candidate_cluster(cluster_id(25))];
    let result_set = spectral_result_set_receipt(&clusters).unwrap().id();
    let estimate = truth_witness(
        problem.problem_id(),
        SpectralTruthPropositionV1::ResultEstimate { result_set },
        25,
    );
    let draft = SpectralTruthDraftV1::new(
        SpectralResultAuthorityV1::Estimated { witness: estimate },
        SpectralCoverageV1::Candidates,
        clusters.clone(),
        ScopeBoundaryStateV1::NoClaim,
        SpectralTerminationV1::Completed,
    );
    assert!(SpectralTruthV1::new(&problem, draft.clone()).is_ok());

    let changed_axes = Axes {
        subject: subject(26),
        ..axes
    };
    let changed_problem = validate_problem(default_spec(
        changed_axes,
        Vec::new(),
        Vec::new(),
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap();
    let report = SpectralTruthV1::new(&changed_problem, draft).unwrap_err();
    assert!(has_truth_issue(&report, |issue| matches!(
        issue,
        SpectralTruthErrorV1::WitnessPropositionMismatch { .. }
    )));

    let estimate_receipt = truth_proposition_receipt(
        problem.problem_id(),
        &SpectralTruthPropositionV1::ResultEstimate { result_set },
    )
    .unwrap();
    let certified_receipt = truth_proposition_receipt(
        problem.problem_id(),
        &SpectralTruthPropositionV1::ResultCertifiedEnclosure { result_set },
    )
    .unwrap();
    assert_ne!(estimate_receipt.id(), certified_receipt.id());
}

#[test]
fn certified_empty_region_is_exact_and_replay_resistant() {
    let region = region_id(27);
    let axes = standard_axes(27, 2);
    let problem = validate_problem(default_spec(
        axes,
        Vec::new(),
        Vec::new(),
        SpectralOrderingV1::NamedRegion { region },
        CompletenessScopeV1::Region {
            region,
            boundary: RegionBoundaryPolicyV1::Closed,
        },
    ))
    .unwrap();
    let result_set = spectral_result_set_receipt(&[]).unwrap().id();
    let boundary = RegionBoundaryStateV1::IntersectionsResolved {
        included: Vec::new(),
        excluded_algebraic: 0,
        witness: truth_witness(
            problem.problem_id(),
            SpectralTruthPropositionV1::RegionBoundaryIntersections {
                result_set,
                included: Vec::new(),
                excluded_algebraic: 0,
            },
            27,
        ),
    };
    let coverage = SpectralCoverageV1::RegionComplete {
        algebraic_cardinality: 0,
        witness: truth_witness(
            problem.problem_id(),
            SpectralTruthPropositionV1::RegionCompleteness {
                result_set,
                algebraic_cardinality: 0,
            },
            28,
        ),
    };
    let draft = SpectralTruthDraftV1::new(
        SpectralResultAuthorityV1::NoClaim,
        coverage,
        Vec::new(),
        ScopeBoundaryStateV1::Region(boundary),
        SpectralTerminationV1::Completed,
    );
    let truth = SpectralTruthV1::new(&problem, draft.clone()).unwrap();
    assert!(truth.clusters().is_empty());

    let other_region = region_id(28);
    let other_problem = validate_problem(default_spec(
        axes,
        Vec::new(),
        Vec::new(),
        SpectralOrderingV1::NamedRegion {
            region: other_region,
        },
        CompletenessScopeV1::Region {
            region: other_region,
            boundary: RegionBoundaryPolicyV1::Closed,
        },
    ))
    .unwrap();
    assert!(SpectralTruthV1::new(&other_problem, draft).is_err());

    let descriptor_class = SpectralProblemClassV1::new(
        SpectralRepresentationV1::GeneralizedPencil,
        DescriptorRoleV1::Descriptor {
            infinity_policy: InfiniteEigenvaluePolicyV1::ExcludeWithCount,
        },
        SpectralOperatorOriginV1::Direct,
    );
    let descriptor_axes = crate::axes(
        28,
        descriptor_class,
        SpectralScalarFieldV1::Complex,
        Dims::NONE,
        2,
    );
    let empty_region_draft = |problem: &ValidatedSpectralProblemV1, seed: u8| {
        let result_set = spectral_result_set_receipt(&[]).unwrap().id();
        SpectralTruthDraftV1::new(
            SpectralResultAuthorityV1::NoClaim,
            SpectralCoverageV1::RegionComplete {
                algebraic_cardinality: 0,
                witness: truth_witness(
                    problem.problem_id(),
                    SpectralTruthPropositionV1::RegionCompleteness {
                        result_set,
                        algebraic_cardinality: 0,
                    },
                    seed,
                ),
            },
            Vec::new(),
            ScopeBoundaryStateV1::Region(RegionBoundaryStateV1::IntersectionsResolved {
                included: Vec::new(),
                excluded_algebraic: 0,
                witness: truth_witness(
                    problem.problem_id(),
                    SpectralTruthPropositionV1::RegionBoundaryIntersections {
                        result_set,
                        included: Vec::new(),
                        excluded_algebraic: 0,
                    },
                    seed.wrapping_add(1),
                ),
            }),
            SpectralTerminationV1::Completed,
        )
    };
    let unestablished = validate_problem(default_spec(
        descriptor_axes,
        Vec::new(),
        Vec::new(),
        SpectralOrderingV1::NamedRegion { region },
        CompletenessScopeV1::Region {
            region,
            boundary: RegionBoundaryPolicyV1::Closed,
        },
    ))
    .unwrap();
    let report =
        SpectralTruthV1::new(&unestablished, empty_region_draft(&unestablished, 29)).unwrap_err();
    assert!(
        report
            .issues()
            .contains(&SpectralTruthErrorV1::DiscreteSpectrumRegularityNotEstablished)
    );

    let established = validate_problem(default_spec(
        descriptor_axes,
        Vec::new(),
        vec![
            regularity_claim(
                descriptor_axes,
                RegularityClassV1::RegularPencil,
                WitnessDispositionV1::Witnessed,
                30,
            ),
            regularity_claim(
                descriptor_axes,
                RegularityClassV1::RegularDescriptor,
                WitnessDispositionV1::Witnessed,
                31,
            ),
        ],
        SpectralOrderingV1::NamedRegion { region },
        CompletenessScopeV1::Region {
            region,
            boundary: RegionBoundaryPolicyV1::Closed,
        },
    ))
    .unwrap();
    assert!(SpectralTruthV1::new(&established, empty_region_draft(&established, 32)).is_ok());
}

#[test]
fn partial_cluster_closure_overrun_never_splits_repeated_boundary_cluster() {
    let axes = standard_axes(29, 4);
    let problem = validate_problem(default_spec(
        axes,
        Vec::new(),
        Vec::new(),
        SpectralOrderingV1::MagnitudeAscending {
            tie_break: ComplexTieBreakV1::RealThenImagThenLineage,
        },
        CompletenessScopeV1::Partial { requested: 2 },
    ))
    .unwrap();
    let first = exact_cluster(
        problem.problem_id(),
        cluster_id(29),
        SpectralEnclosureV1::Real(FiniteIntervalV1::new(1.0, 1.0).unwrap()),
        1,
        29,
    );
    let boundary_cluster = exact_cluster(
        problem.problem_id(),
        cluster_id(30),
        SpectralEnclosureV1::Real(FiniteIntervalV1::new(2.0, 2.0).unwrap()),
        2,
        30,
    );
    let clusters = vec![first, boundary_cluster];
    let result_set = spectral_result_set_receipt(&clusters).unwrap().id();
    let status = PartialCoverageStatusV1::ClusterClosureOverrun {
        boundary_cluster: boundary_cluster.id(),
        preceding_algebraic: 1,
    };
    let boundary = ScopeBoundaryStateV1::Partial(PartialBoundaryStateV1::ClusterClosed {
        cluster: boundary_cluster.id(),
        witness: truth_witness(
            problem.problem_id(),
            SpectralTruthPropositionV1::PartialBoundaryClusterClosed {
                result_set,
                cluster: boundary_cluster.id(),
            },
            31,
        ),
    });
    let coverage = SpectralCoverageV1::Partial {
        returned_algebraic: 3,
        status,
        witness: truth_witness(
            problem.problem_id(),
            SpectralTruthPropositionV1::PartialCoverage {
                result_set,
                returned_algebraic: 3,
                status,
            },
            32,
        ),
    };
    let draft = SpectralTruthDraftV1::new(
        SpectralResultAuthorityV1::NoClaim,
        coverage,
        clusters.clone(),
        boundary,
        SpectralTerminationV1::Completed,
    );
    assert!(SpectralTruthV1::new(&problem, draft.clone()).is_ok());

    let no_boundary = SpectralTruthDraftV1::new(
        SpectralResultAuthorityV1::NoClaim,
        coverage,
        clusters,
        ScopeBoundaryStateV1::NoClaim,
        SpectralTerminationV1::Completed,
    );
    let report = SpectralTruthV1::new(&problem, no_boundary).unwrap_err();
    assert!(
        report
            .issues()
            .contains(&SpectralTruthErrorV1::InvalidPartialCoverage)
    );
}

#[test]
fn ordinary_full_spectrum_requires_exact_count_and_full_boundary() {
    let axes = standard_axes(33, 2);
    let problem = validate_problem(default_spec(
        axes,
        Vec::new(),
        Vec::new(),
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::FullFinite {
            algebraic_cardinality: 2,
            infinity_policy: InfiniteEigenvaluePolicyV1::NoClaim,
        },
    ))
    .unwrap();
    let clusters = vec![
        exact_cluster(
            problem.problem_id(),
            cluster_id(33),
            SpectralEnclosureV1::Real(FiniteIntervalV1::new(1.0, 1.0).unwrap()),
            1,
            33,
        ),
        exact_cluster(
            problem.problem_id(),
            cluster_id(34),
            SpectralEnclosureV1::Real(FiniteIntervalV1::new(2.0, 2.0).unwrap()),
            1,
            34,
        ),
    ];
    let result_set = spectral_result_set_receipt(&clusters).unwrap().id();
    let coverage = SpectralCoverageV1::FullFinite {
        finite_algebraic: 2,
        infinity: InfinityAccountingV1::NotApplicable,
        witness: truth_witness(
            problem.problem_id(),
            SpectralTruthPropositionV1::FullCompleteness {
                result_set,
                finite_algebraic: 2,
                infinity: InfinityAccountingStatementV1::NotApplicable,
            },
            35,
        ),
    };
    let draft = SpectralTruthDraftV1::new(
        SpectralResultAuthorityV1::NoClaim,
        coverage,
        clusters.clone(),
        ScopeBoundaryStateV1::FullSpectrum,
        SpectralTerminationV1::Completed,
    );
    assert!(SpectralTruthV1::new(&problem, draft.clone()).is_ok());

    let report = SpectralTruthV1::new(
        &problem,
        SpectralTruthDraftV1::new(
            SpectralResultAuthorityV1::NoClaim,
            coverage,
            clusters,
            ScopeBoundaryStateV1::NoClaim,
            SpectralTerminationV1::Completed,
        ),
    )
    .unwrap_err();
    assert!(
        report
            .issues()
            .contains(&SpectralTruthErrorV1::BoundaryCoverageMismatch)
    );
}

fn descriptor_full_problem(
    seed: u8,
    policy: InfiniteEigenvaluePolicyV1,
) -> ValidatedSpectralProblemV1 {
    let class = SpectralProblemClassV1::new(
        SpectralRepresentationV1::GeneralizedPencil,
        DescriptorRoleV1::Descriptor {
            infinity_policy: policy,
        },
        SpectralOperatorOriginV1::Direct,
    );
    let axes = axes(seed, class, SpectralScalarFieldV1::Complex, Dims::NONE, 2);
    let regularity = vec![
        regularity_claim(
            axes,
            RegularityClassV1::RegularPencil,
            WitnessDispositionV1::Witnessed,
            seed,
        ),
        regularity_claim(
            axes,
            RegularityClassV1::RegularDescriptor,
            WitnessDispositionV1::Witnessed,
            seed.wrapping_add(1),
        ),
    ];
    validate_problem(default_spec(
        axes,
        Vec::new(),
        regularity,
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::FullFinite {
            algebraic_cardinality: 2,
            infinity_policy: policy,
        },
    ))
    .unwrap()
}

fn excluded_full_draft(
    problem: &ValidatedSpectralProblemV1,
    finite_algebraic: u32,
    excluded_algebraic: u32,
    seed: u8,
) -> SpectralTruthDraftV1 {
    let clusters = vec![exact_cluster(
        problem.problem_id(),
        cluster_id(seed),
        SpectralEnclosureV1::Real(FiniteIntervalV1::new(1.0, 1.0).unwrap()),
        finite_algebraic,
        seed,
    )];
    let result_set = spectral_result_set_receipt(&clusters).unwrap().id();
    let infinity = InfinityAccountingV1::ExcludedWithCount {
        algebraic: excluded_algebraic,
        witness: truth_witness(
            problem.problem_id(),
            SpectralTruthPropositionV1::ExcludedInfinity {
                result_set,
                algebraic: excluded_algebraic,
            },
            seed.wrapping_add(1),
        ),
    };
    SpectralTruthDraftV1::new(
        SpectralResultAuthorityV1::NoClaim,
        SpectralCoverageV1::FullFinite {
            finite_algebraic,
            infinity,
            witness: truth_witness(
                problem.problem_id(),
                SpectralTruthPropositionV1::FullCompleteness {
                    result_set,
                    finite_algebraic,
                    infinity: InfinityAccountingStatementV1::ExcludedWithCount {
                        algebraic: excluded_algebraic,
                    },
                },
                seed.wrapping_add(2),
            ),
        },
        clusters,
        ScopeBoundaryStateV1::FullSpectrum,
        SpectralTerminationV1::Completed,
    )
}

fn satisfied_partial_draft(problem: &ValidatedSpectralProblemV1, seed: u8) -> SpectralTruthDraftV1 {
    let clusters = vec![exact_cluster(
        problem.problem_id(),
        cluster_id(seed),
        SpectralEnclosureV1::Real(FiniteIntervalV1::new(1.0, 1.0).unwrap()),
        1,
        seed,
    )];
    let result_set = spectral_result_set_receipt(&clusters).unwrap().id();
    let lower = PositiveFiniteV1::new(0.5).unwrap();
    let norm = norm_id(seed);
    let status = PartialCoverageStatusV1::Satisfied;
    SpectralTruthDraftV1::new(
        SpectralResultAuthorityV1::NoClaim,
        SpectralCoverageV1::Partial {
            returned_algebraic: 1,
            status,
            witness: truth_witness(
                problem.problem_id(),
                SpectralTruthPropositionV1::PartialCoverage {
                    result_set,
                    returned_algebraic: 1,
                    status,
                },
                seed.wrapping_add(1),
            ),
        },
        clusters,
        ScopeBoundaryStateV1::Partial(PartialBoundaryStateV1::Separated {
            lower,
            norm,
            witness: truth_witness(
                problem.problem_id(),
                SpectralTruthPropositionV1::PartialBoundarySeparated {
                    result_set,
                    lower,
                    norm,
                },
                seed.wrapping_add(2),
            ),
        }),
        SpectralTerminationV1::Completed,
    )
}

#[test]
fn incomplete_partial_truth_retains_bounded_evidence_without_claiming_completion() {
    let axes = standard_axes(35, 3);
    let problem = validate_problem(default_spec(
        axes,
        Vec::new(),
        Vec::new(),
        SpectralOrderingV1::MagnitudeAscending {
            tie_break: ComplexTieBreakV1::RealThenImagThenLineage,
        },
        CompletenessScopeV1::Partial { requested: 2 },
    ))
    .unwrap();
    let clusters = vec![exact_cluster(
        problem.problem_id(),
        cluster_id(35),
        SpectralEnclosureV1::Real(FiniteIntervalV1::new(1.0, 1.0).unwrap()),
        1,
        35,
    )];
    let result_set = spectral_result_set_receipt(&clusters).unwrap().id();
    let status = PartialCoverageStatusV1::Incomplete;
    let truth = SpectralTruthV1::new(
        &problem,
        SpectralTruthDraftV1::new(
            SpectralResultAuthorityV1::NoClaim,
            SpectralCoverageV1::Partial {
                returned_algebraic: 1,
                status,
                witness: truth_witness(
                    problem.problem_id(),
                    SpectralTruthPropositionV1::PartialCoverage {
                        result_set,
                        returned_algebraic: 1,
                        status,
                    },
                    36,
                ),
            },
            clusters,
            ScopeBoundaryStateV1::Partial(PartialBoundaryStateV1::Unknown {
                reason: UnknownSeparationReasonV1::BudgetExhausted,
            }),
            SpectralTerminationV1::BudgetExhausted,
        ),
    )
    .unwrap();
    assert!(matches!(
        truth.coverage(),
        SpectralCoverageV1::Partial {
            returned_algebraic: 1,
            status: PartialCoverageStatusV1::Incomplete,
            ..
        }
    ));
    assert_eq!(truth.termination(), SpectralTerminationV1::BudgetExhausted);
}

#[test]
fn descriptor_full_truth_accounts_for_included_and_excluded_infinity() {
    let included_problem =
        descriptor_full_problem(36, InfiniteEigenvaluePolicyV1::IncludeProjective);
    let finite = exact_cluster(
        included_problem.problem_id(),
        cluster_id(36),
        SpectralEnclosureV1::ComplexBox {
            real: FiniteIntervalV1::new(1.0, 1.0).unwrap(),
            imag: FiniteIntervalV1::new(0.0, 0.0).unwrap(),
        },
        1,
        36,
    );
    let projective = exact_cluster(
        included_problem.problem_id(),
        cluster_id(37),
        SpectralEnclosureV1::ProjectiveInfinity,
        1,
        37,
    );
    let included_clusters = vec![finite, projective];
    let included_set = spectral_result_set_receipt(&included_clusters)
        .unwrap()
        .id();
    let infinity = InfinityAccountingV1::Included {
        algebraic: 1,
        cluster: Some(projective.id()),
        witness: truth_witness(
            included_problem.problem_id(),
            SpectralTruthPropositionV1::IncludedInfinity {
                result_set: included_set,
                algebraic: 1,
                cluster: Some(projective.id()),
            },
            38,
        ),
    };
    let included_coverage = SpectralCoverageV1::FullFinite {
        finite_algebraic: 1,
        infinity,
        witness: truth_witness(
            included_problem.problem_id(),
            SpectralTruthPropositionV1::FullCompleteness {
                result_set: included_set,
                finite_algebraic: 1,
                infinity: InfinityAccountingStatementV1::Included {
                    algebraic: 1,
                    cluster: Some(projective.id()),
                },
            },
            39,
        ),
    };
    assert!(
        SpectralTruthV1::new(
            &included_problem,
            SpectralTruthDraftV1::new(
                SpectralResultAuthorityV1::NoClaim,
                included_coverage,
                included_clusters,
                ScopeBoundaryStateV1::FullSpectrum,
                SpectralTerminationV1::Completed,
            ),
        )
        .is_ok()
    );

    let excluded_problem =
        descriptor_full_problem(40, InfiniteEigenvaluePolicyV1::ExcludeWithCount);
    let excluded_clusters = vec![exact_cluster(
        excluded_problem.problem_id(),
        cluster_id(40),
        SpectralEnclosureV1::Real(FiniteIntervalV1::new(1.0, 1.0).unwrap()),
        1,
        40,
    )];
    let excluded_set = spectral_result_set_receipt(&excluded_clusters)
        .unwrap()
        .id();
    let infinity = InfinityAccountingV1::ExcludedWithCount {
        algebraic: 1,
        witness: truth_witness(
            excluded_problem.problem_id(),
            SpectralTruthPropositionV1::ExcludedInfinity {
                result_set: excluded_set,
                algebraic: 1,
            },
            41,
        ),
    };
    let excluded_coverage = SpectralCoverageV1::FullFinite {
        finite_algebraic: 1,
        infinity,
        witness: truth_witness(
            excluded_problem.problem_id(),
            SpectralTruthPropositionV1::FullCompleteness {
                result_set: excluded_set,
                finite_algebraic: 1,
                infinity: InfinityAccountingStatementV1::ExcludedWithCount { algebraic: 1 },
            },
            42,
        ),
    };
    assert!(
        SpectralTruthV1::new(
            &excluded_problem,
            SpectralTruthDraftV1::new(
                SpectralResultAuthorityV1::NoClaim,
                excluded_coverage,
                excluded_clusters,
                ScopeBoundaryStateV1::FullSpectrum,
                SpectralTerminationV1::Completed,
            ),
        )
        .is_ok()
    );
}

#[test]
fn full_completeness_requires_admitted_regularity_and_consumes_theorem_closure() {
    let finite_space = standard_axes(42, 2);
    let report = validate_problem(default_spec(
        finite_space,
        Vec::new(),
        vec![regularity_claim(
            finite_space,
            RegularityClassV1::FiniteDimensional,
            WitnessDispositionV1::Contradicted,
            42,
        )],
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap_err();
    assert!(has_admission_issue(&report, |issue| matches!(
        issue,
        SpectralAdmissionIssueV1::RegularityMismatch
    )));

    let generalized_class = SpectralProblemClassV1::new(
        SpectralRepresentationV1::GeneralizedPencil,
        DescriptorRoleV1::Descriptor {
            infinity_policy: InfiniteEigenvaluePolicyV1::ExcludeWithCount,
        },
        SpectralOperatorOriginV1::Direct,
    );
    let generalized = axes(
        43,
        generalized_class,
        SpectralScalarFieldV1::Complex,
        Dims::NONE,
        2,
    );
    let scope = CompletenessScopeV1::FullFinite {
        algebraic_cardinality: 2,
        infinity_policy: InfiniteEigenvaluePolicyV1::ExcludeWithCount,
    };

    let report = validate_problem(default_spec(
        generalized,
        Vec::new(),
        vec![regularity_claim(
            generalized,
            RegularityClassV1::FiniteDimensional,
            WitnessDispositionV1::Contradicted,
            43,
        )],
        SpectralOrderingV1::SetValued,
        scope,
    ))
    .unwrap_err();
    assert!(has_admission_issue(&report, |issue| matches!(
        issue,
        SpectralAdmissionIssueV1::RegularityMismatch
    )));

    let contradictory_equation = validate_problem(default_spec(
        generalized,
        Vec::new(),
        vec![
            regularity_claim(
                generalized,
                RegularityClassV1::RegularPencil,
                WitnessDispositionV1::Contradicted,
                43,
            ),
            regularity_claim(
                generalized,
                RegularityClassV1::RegularDescriptor,
                WitnessDispositionV1::Witnessed,
                44,
            ),
        ],
        SpectralOrderingV1::SetValued,
        scope,
    ))
    .unwrap();
    assert_eq!(contradictory_equation.known_algebraic_cardinality(), None);
    let report = SpectralTruthV1::new(
        &contradictory_equation,
        excluded_full_draft(&contradictory_equation, 1, 1, 45),
    )
    .unwrap_err();
    assert!(report.issues().contains(
        &SpectralTruthErrorV1::FullCompletenessRegularityNotEstablished {
            requested: 2,
            established: None,
        }
    ));

    let contradictory_descriptor = validate_problem(default_spec(
        generalized,
        Vec::new(),
        vec![
            regularity_claim(
                generalized,
                RegularityClassV1::RegularPencil,
                WitnessDispositionV1::Witnessed,
                46,
            ),
            regularity_claim(
                generalized,
                RegularityClassV1::RegularDescriptor,
                WitnessDispositionV1::Contradicted,
                47,
            ),
        ],
        SpectralOrderingV1::SetValued,
        scope,
    ))
    .unwrap();
    assert_eq!(contradictory_descriptor.known_algebraic_cardinality(), None);
    let report = SpectralTruthV1::new(
        &contradictory_descriptor,
        excluded_full_draft(&contradictory_descriptor, 1, 1, 48),
    )
    .unwrap_err();
    assert!(report.issues().contains(
        &SpectralTruthErrorV1::FullCompletenessRegularityNotEstablished {
            requested: 2,
            established: None,
        }
    ));

    let unestablished = validate_problem(default_spec(
        generalized,
        Vec::new(),
        Vec::new(),
        SpectralOrderingV1::SetValued,
        scope,
    ))
    .unwrap();
    assert_eq!(unestablished.known_algebraic_cardinality(), None);
    let report = SpectralTruthV1::new(
        &unestablished,
        excluded_full_draft(&unestablished, 1, 1, 49),
    )
    .unwrap_err();
    assert!(report.issues().contains(
        &SpectralTruthErrorV1::FullCompletenessRegularityNotEstablished {
            requested: 2,
            established: None,
        }
    ));

    let hdp = structure_claim(
        generalized,
        StructurePropertyV1::HermitianDefinitePencil,
        StructureSupportV1::InnerProduct(generalized.domain.id()),
        WitnessDispositionV1::Witnessed,
        0.0,
        50,
    );
    let theorem_backed = validate_problem(default_spec(
        generalized,
        vec![hdp],
        vec![regularity_claim(
            generalized,
            RegularityClassV1::RegularDescriptor,
            WitnessDispositionV1::Witnessed,
            51,
        )],
        SpectralOrderingV1::SetValued,
        scope,
    ))
    .unwrap();
    assert_eq!(theorem_backed.known_algebraic_cardinality(), Some(2));
    assert!(
        SpectralTruthV1::new(
            &theorem_backed,
            excluded_full_draft(&theorem_backed, 2, 0, 52),
        )
        .is_ok()
    );

    let weight_backed = validate_problem(default_spec(
        generalized,
        Vec::new(),
        vec![
            regularity_claim(
                generalized,
                RegularityClassV1::InvertiblePencilWeight,
                WitnessDispositionV1::Witnessed,
                59,
            ),
            regularity_claim(
                generalized,
                RegularityClassV1::RegularDescriptor,
                WitnessDispositionV1::Witnessed,
                60,
            ),
        ],
        SpectralOrderingV1::SetValued,
        scope,
    ))
    .unwrap();
    assert_eq!(weight_backed.known_algebraic_cardinality(), Some(2));
    let report = SpectralTruthV1::new(
        &weight_backed,
        excluded_full_draft(&weight_backed, 1, 1, 61),
    )
    .unwrap_err();
    assert!(
        report
            .issues()
            .contains(&SpectralTruthErrorV1::InfinityAccountingMismatch)
    );

    let polynomial_class = SpectralProblemClassV1::new(
        SpectralRepresentationV1::MatrixPolynomial { grade: 2 },
        DescriptorRoleV1::Descriptor {
            infinity_policy: InfiniteEigenvaluePolicyV1::ExcludeWithCount,
        },
        SpectralOperatorOriginV1::Direct,
    );
    let polynomial = axes(
        53,
        polynomial_class,
        SpectralScalarFieldV1::Complex,
        Dims::NONE,
        2,
    );
    let polynomial_scope = CompletenessScopeV1::FullFinite {
        algebraic_cardinality: 4,
        infinity_policy: InfiniteEigenvaluePolicyV1::ExcludeWithCount,
    };
    let contradicted_polynomial = validate_problem(default_spec(
        polynomial,
        Vec::new(),
        vec![
            regularity_claim(
                polynomial,
                RegularityClassV1::RegularPolynomial { grade: 2 },
                WitnessDispositionV1::Contradicted,
                53,
            ),
            regularity_claim(
                polynomial,
                RegularityClassV1::RegularDescriptor,
                WitnessDispositionV1::Witnessed,
                54,
            ),
        ],
        SpectralOrderingV1::SetValued,
        polynomial_scope,
    ))
    .unwrap();
    assert_eq!(contradicted_polynomial.known_algebraic_cardinality(), None);
    let report = SpectralTruthV1::new(
        &contradicted_polynomial,
        excluded_full_draft(&contradicted_polynomial, 3, 1, 55),
    )
    .unwrap_err();
    assert!(report.issues().contains(
        &SpectralTruthErrorV1::FullCompletenessRegularityNotEstablished {
            requested: 4,
            established: None,
        }
    ));

    let leading_backed = validate_problem(default_spec(
        polynomial,
        Vec::new(),
        vec![
            regularity_claim(
                polynomial,
                RegularityClassV1::InvertiblePolynomialLeadingCoefficient { grade: 2 },
                WitnessDispositionV1::Witnessed,
                56,
            ),
            regularity_claim(
                polynomial,
                RegularityClassV1::RegularDescriptor,
                WitnessDispositionV1::Witnessed,
                57,
            ),
        ],
        SpectralOrderingV1::SetValued,
        polynomial_scope,
    ))
    .unwrap();
    assert_eq!(leading_backed.known_algebraic_cardinality(), Some(4));
    let report = SpectralTruthV1::new(
        &leading_backed,
        excluded_full_draft(&leading_backed, 3, 1, 58),
    )
    .unwrap_err();
    assert!(
        report
            .issues()
            .contains(&SpectralTruthErrorV1::InfinityAccountingMismatch)
    );
    assert!(
        SpectralTruthV1::new(
            &leading_backed,
            excluded_full_draft(&leading_backed, 4, 0, 62),
        )
        .is_ok()
    );
}

#[test]
fn invertibility_excludes_favorable_projective_truth_in_every_coverage_mode() {
    let generalized_class = SpectralProblemClassV1::new(
        SpectralRepresentationV1::GeneralizedPencil,
        DescriptorRoleV1::Descriptor {
            infinity_policy: InfiniteEigenvaluePolicyV1::IncludeProjective,
        },
        SpectralOperatorOriginV1::Direct,
    );
    let generalized = axes(
        90,
        generalized_class,
        SpectralScalarFieldV1::Complex,
        Dims::NONE,
        2,
    );
    let generalized_regularity = vec![
        regularity_claim(
            generalized,
            RegularityClassV1::InvertiblePencilWeight,
            WitnessDispositionV1::Witnessed,
            90,
        ),
        regularity_claim(
            generalized,
            RegularityClassV1::RegularDescriptor,
            WitnessDispositionV1::Witnessed,
            91,
        ),
    ];
    let generalized_candidates = validate_problem(default_spec(
        generalized,
        Vec::new(),
        generalized_regularity.clone(),
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap();

    let raw_id = cluster_id(90);
    let raw_projective = SpectralClusterV1::new(
        raw_id,
        SpectralLocalizationV1::candidate(SpectralEnclosureV1::ProjectiveInfinity),
        MultiplicityClaimV1::Unknown,
        MultiplicityClaimV1::Unknown,
        InternalClusterStateV1::NoClaim,
    )
    .unwrap();
    assert!(
        SpectralTruthV1::new(
            &generalized_candidates,
            SpectralTruthDraftV1::new(
                SpectralResultAuthorityV1::Candidate,
                SpectralCoverageV1::Candidates,
                vec![raw_projective],
                ScopeBoundaryStateV1::NoClaim,
                SpectralTerminationV1::Completed,
            ),
        )
        .is_ok(),
        "an explicitly non-authoritative projective diagnostic must remain representable"
    );

    let enclosed_id = cluster_id(91);
    let enclosed_projective = SpectralClusterV1::new(
        enclosed_id,
        SpectralLocalizationV1::enclosed(
            SpectralEnclosureV1::ProjectiveInfinity,
            truth_witness(
                generalized_candidates.problem_id(),
                SpectralTruthPropositionV1::ClusterLocalization {
                    cluster: enclosed_id,
                    authority: LocalizationAuthorityV1::Enclosed,
                    enclosure: SpectralEnclosureV1::ProjectiveInfinity,
                },
                92,
            ),
        ),
        MultiplicityClaimV1::Unknown,
        MultiplicityClaimV1::Unknown,
        InternalClusterStateV1::NoClaim,
    )
    .unwrap();
    let report = SpectralTruthV1::new(
        &generalized_candidates,
        SpectralTruthDraftV1::new(
            SpectralResultAuthorityV1::NoClaim,
            SpectralCoverageV1::Candidates,
            vec![enclosed_projective],
            ScopeBoundaryStateV1::NoClaim,
            SpectralTerminationV1::Completed,
        ),
    )
    .unwrap_err();
    assert!(
        report
            .issues()
            .contains(&SpectralTruthErrorV1::ProjectiveInfinityExcludedByRegularity)
    );

    let generalized_partial = validate_problem(default_spec(
        generalized,
        Vec::new(),
        generalized_regularity,
        SpectralOrderingV1::Projective {
            chart: SpectralProjectiveChartId::from_bytes([90; 32]),
            infinity: ProjectiveInfinityPlacementV1::Last,
            tie_break: ComplexTieBreakV1::RealThenImagThenLineage,
        },
        CompletenessScopeV1::Partial { requested: 1 },
    ))
    .unwrap();
    let projective = exact_cluster(
        generalized_partial.problem_id(),
        cluster_id(92),
        SpectralEnclosureV1::ProjectiveInfinity,
        1,
        93,
    );
    let result_set = spectral_result_set_receipt(&[projective]).unwrap().id();
    let status = PartialCoverageStatusV1::Satisfied;
    let lower = PositiveFiniteV1::new(0.25).unwrap();
    let norm = norm_id(92);
    let report = SpectralTruthV1::new(
        &generalized_partial,
        SpectralTruthDraftV1::new(
            SpectralResultAuthorityV1::NoClaim,
            SpectralCoverageV1::Partial {
                returned_algebraic: 1,
                status,
                witness: truth_witness(
                    generalized_partial.problem_id(),
                    SpectralTruthPropositionV1::PartialCoverage {
                        result_set,
                        returned_algebraic: 1,
                        status,
                    },
                    94,
                ),
            },
            vec![projective],
            ScopeBoundaryStateV1::Partial(PartialBoundaryStateV1::Separated {
                lower,
                norm,
                witness: truth_witness(
                    generalized_partial.problem_id(),
                    SpectralTruthPropositionV1::PartialBoundarySeparated {
                        result_set,
                        lower,
                        norm,
                    },
                    95,
                ),
            }),
            SpectralTerminationV1::Completed,
        ),
    )
    .unwrap_err();
    assert!(
        report
            .issues()
            .contains(&SpectralTruthErrorV1::ProjectiveInfinityExcludedByRegularity)
    );

    let polynomial_class = SpectralProblemClassV1::new(
        SpectralRepresentationV1::MatrixPolynomial { grade: 2 },
        DescriptorRoleV1::Descriptor {
            infinity_policy: InfiniteEigenvaluePolicyV1::IncludeProjective,
        },
        SpectralOperatorOriginV1::Direct,
    );
    let polynomial = axes(
        96,
        polynomial_class,
        SpectralScalarFieldV1::Complex,
        Dims::NONE,
        2,
    );
    let polynomial_candidates = validate_problem(default_spec(
        polynomial,
        Vec::new(),
        vec![
            regularity_claim(
                polynomial,
                RegularityClassV1::InvertiblePolynomialLeadingCoefficient { grade: 2 },
                WitnessDispositionV1::Witnessed,
                96,
            ),
            regularity_claim(
                polynomial,
                RegularityClassV1::RegularDescriptor,
                WitnessDispositionV1::Witnessed,
                97,
            ),
        ],
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap();
    let projective = SpectralClusterV1::new(
        cluster_id(96),
        SpectralLocalizationV1::candidate(SpectralEnclosureV1::ProjectiveInfinity),
        MultiplicityClaimV1::Unknown,
        MultiplicityClaimV1::Unknown,
        InternalClusterStateV1::NoClaim,
    )
    .unwrap();
    let result_set = spectral_result_set_receipt(&[projective]).unwrap().id();
    let report = SpectralTruthV1::new(
        &polynomial_candidates,
        SpectralTruthDraftV1::new(
            SpectralResultAuthorityV1::Estimated {
                witness: truth_witness(
                    polynomial_candidates.problem_id(),
                    SpectralTruthPropositionV1::ResultEstimate { result_set },
                    98,
                ),
            },
            SpectralCoverageV1::Candidates,
            vec![projective],
            ScopeBoundaryStateV1::NoClaim,
            SpectralTerminationV1::Completed,
        ),
    )
    .unwrap_err();
    assert!(
        report
            .issues()
            .contains(&SpectralTruthErrorV1::ProjectiveInfinityExcludedByRegularity)
    );
}

#[test]
fn projective_clusters_and_no_result_claims_fail_closed() {
    let axes = standard_axes(43, 2);
    let problem = validate_problem(default_spec(
        axes,
        Vec::new(),
        Vec::new(),
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap();
    let projective = candidate_cluster(cluster_id(43));
    let projective = SpectralClusterV1::new(
        projective.id(),
        SpectralLocalizationV1::candidate(SpectralEnclosureV1::ProjectiveInfinity),
        MultiplicityClaimV1::Unknown,
        MultiplicityClaimV1::Unknown,
        InternalClusterStateV1::NoClaim,
    )
    .unwrap();
    let report = SpectralTruthV1::new(
        &problem,
        SpectralTruthDraftV1::new(
            SpectralResultAuthorityV1::Candidate,
            SpectralCoverageV1::Candidates,
            vec![projective],
            ScopeBoundaryStateV1::NoClaim,
            SpectralTerminationV1::Completed,
        ),
    )
    .unwrap_err();
    assert!(
        report
            .issues()
            .contains(&SpectralTruthErrorV1::ProjectiveClusterNotAdmitted)
    );

    let report = SpectralTruthV1::new(
        &problem,
        SpectralTruthDraftV1::new(
            SpectralResultAuthorityV1::Candidate,
            SpectralCoverageV1::NoResult,
            Vec::new(),
            ScopeBoundaryStateV1::NoClaim,
            SpectralTerminationV1::NumericalFailure,
        ),
    )
    .unwrap_err();
    assert!(
        report
            .issues()
            .contains(&SpectralTruthErrorV1::NoResultHasClaims)
    );
}

#[test]
fn cluster_truth_rejects_impossible_multiplicity_and_internal_states() {
    let id = cluster_id(44);
    let interval = SpectralEnclosureV1::Real(FiniteIntervalV1::new(0.0, 1.0).unwrap());
    assert_eq!(
        SpectralClusterV1::new(
            id,
            SpectralLocalizationV1::candidate(interval),
            MultiplicityClaimV1::Unknown,
            MultiplicityClaimV1::Unknown,
            InternalClusterStateV1::Simple,
        ),
        Err(SpectralTruthErrorV1::InvalidInternalClusterState)
    );

    let axes = standard_axes(44, 2);
    let problem = validate_problem(default_spec(
        axes,
        Vec::new(),
        Vec::new(),
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap();
    let exact_one = exact_cluster(problem.problem_id(), id, interval, 1, 44);
    let truth = SpectralTruthV1::new(
        &problem,
        SpectralTruthDraftV1::new(
            SpectralResultAuthorityV1::NoClaim,
            SpectralCoverageV1::Candidates,
            vec![exact_one],
            ScopeBoundaryStateV1::NoClaim,
            SpectralTerminationV1::Completed,
        ),
    )
    .unwrap();
    let exact_one = &truth.clusters()[0];
    assert_eq!(exact_one.defectivity(), DefectivityStateV1::Unknown);
    assert!(matches!(
        exact_one.internal(),
        InternalClusterStateV1::Simple
    ));

    let repeated_id = cluster_id(45);
    let algebraic = MultiplicityClaimV1::Exact {
        value: 2,
        witness: truth_witness(
            problem.problem_id(),
            SpectralTruthPropositionV1::Multiplicity {
                cluster: repeated_id,
                enclosure: interval,
                kind: MultiplicityKindV1::Algebraic,
                assertion: MultiplicityAssertionV1::Exact,
                lower: 2,
                upper: Some(2),
            },
            45,
        ),
    };
    let geometric = MultiplicityClaimV1::Exact {
        value: 1,
        witness: truth_witness(
            problem.problem_id(),
            SpectralTruthPropositionV1::Multiplicity {
                cluster: repeated_id,
                enclosure: interval,
                kind: MultiplicityKindV1::Geometric,
                assertion: MultiplicityAssertionV1::Exact,
                lower: 1,
                upper: Some(1),
            },
            46,
        ),
    };
    let repeated = SpectralClusterV1::new(
        repeated_id,
        SpectralLocalizationV1::candidate(interval),
        algebraic,
        geometric,
        InternalClusterStateV1::Unknown {
            reason: UnknownSeparationReasonV1::MissingEvidence,
        },
    )
    .unwrap();
    let truth = SpectralTruthV1::new(
        &problem,
        SpectralTruthDraftV1::new(
            SpectralResultAuthorityV1::NoClaim,
            SpectralCoverageV1::Candidates,
            vec![repeated],
            ScopeBoundaryStateV1::NoClaim,
            SpectralTerminationV1::Completed,
        ),
    )
    .unwrap();
    assert_eq!(
        truth.clusters()[0].defectivity(),
        DefectivityStateV1::ProvenDefective
    );
}

#[test]
fn truth_resource_caps_precede_sorting_hashing_and_reference_scans() {
    let axes = standard_axes(45, 2);
    let problem = validate_problem(default_spec(
        axes,
        Vec::new(),
        Vec::new(),
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap();
    let clusters: Vec<_> = (0..=MAX_SPECTRAL_CLUSTERS_V1)
        .map(|index| candidate_cluster(indexed_cluster_id(index)))
        .collect();
    let report = SpectralTruthV1::new(
        &problem,
        SpectralTruthDraftV1::new(
            SpectralResultAuthorityV1::NoClaim,
            SpectralCoverageV1::Candidates,
            clusters,
            ScopeBoundaryStateV1::NoClaim,
            SpectralTerminationV1::Completed,
        ),
    )
    .unwrap_err();
    assert_eq!(
        report.issues(),
        &[SpectralTruthErrorV1::TooManyClusters {
            found: MAX_SPECTRAL_CLUSTERS_V1 + 1,
            limit: MAX_SPECTRAL_CLUSTERS_V1,
        }]
    );

    let arbitrary = admit_receipt(
        regularity_proposition_receipt(
            axes.subject,
            axes.scalar,
            axes.class,
            axes.scaling,
            axes.domain,
            axes.codomain,
            RegularityClassV1::FiniteDimensional,
            WitnessDispositionV1::Witnessed,
        )
        .unwrap(),
        45,
    );
    let references = vec![cluster_id(45); MAX_REGION_BOUNDARY_REFERENCES_V1 + 1];
    let report = SpectralTruthV1::new(
        &problem,
        SpectralTruthDraftV1::new(
            SpectralResultAuthorityV1::NoClaim,
            SpectralCoverageV1::Candidates,
            Vec::new(),
            ScopeBoundaryStateV1::Region(RegionBoundaryStateV1::IntersectionsResolved {
                included: references,
                excluded_algebraic: 0,
                witness: arbitrary,
            }),
            SpectralTerminationV1::Completed,
        ),
    )
    .unwrap_err();
    assert_eq!(
        report.issues(),
        &[SpectralTruthErrorV1::TooManyBoundaryReferences {
            found: MAX_REGION_BOUNDARY_REFERENCES_V1 + 1,
            limit: MAX_REGION_BOUNDARY_REFERENCES_V1,
        }]
    );
}

#[test]
fn deterministic_reports_canonicalize_cluster_and_reference_order() {
    let region = region_id(46);
    let axes = standard_axes(46, 2);
    let problem = validate_problem(default_spec(
        axes,
        Vec::new(),
        Vec::new(),
        SpectralOrderingV1::NamedRegion { region },
        CompletenessScopeV1::Region {
            region,
            boundary: RegionBoundaryPolicyV1::Closed,
        },
    ))
    .unwrap();
    let first = candidate_cluster(cluster_id(46));
    let second = candidate_cluster(cluster_id(47));
    let clusters = vec![first, second];
    let result_set = spectral_result_set_receipt(&clusters).unwrap().id();
    let included = vec![second.id(), first.id(), first.id()];
    let witness = truth_witness(
        problem.problem_id(),
        SpectralTruthPropositionV1::RegionBoundaryIntersections {
            result_set,
            included: included.clone(),
            excluded_algebraic: 0,
        },
        46,
    );
    let make = |clusters: Vec<SpectralClusterV1>, included: Vec<SpectralClusterIdV1>| {
        SpectralTruthDraftV1::new(
            SpectralResultAuthorityV1::NoClaim,
            SpectralCoverageV1::Candidates,
            clusters,
            ScopeBoundaryStateV1::Region(RegionBoundaryStateV1::IntersectionsResolved {
                included,
                excluded_algebraic: 0,
                witness,
            }),
            SpectralTerminationV1::Completed,
        )
    };
    let first_report =
        SpectralTruthV1::new(&problem, make(vec![first, second], included.clone())).unwrap_err();
    let second_report = SpectralTruthV1::new(
        &problem,
        make(
            vec![second, first],
            vec![first.id(), second.id(), first.id()],
        ),
    )
    .unwrap_err();
    assert_eq!(first_report, second_report);
    assert!(
        first_report
            .issues()
            .contains(&SpectralTruthErrorV1::DuplicateSeparationReference)
    );

    let duplicate_id = cluster_id(48);
    let duplicate_a = candidate_cluster(duplicate_id);
    let duplicate_b = SpectralClusterV1::new(
        duplicate_id,
        SpectralLocalizationV1::candidate(SpectralEnclosureV1::Real(
            FiniteIntervalV1::new(2.0, 3.0).unwrap(),
        )),
        MultiplicityClaimV1::Unknown,
        MultiplicityClaimV1::Unknown,
        InternalClusterStateV1::NoClaim,
    )
    .unwrap();
    let duplicate_draft = |clusters| {
        SpectralTruthDraftV1::new(
            SpectralResultAuthorityV1::NoClaim,
            SpectralCoverageV1::Candidates,
            clusters,
            ScopeBoundaryStateV1::NoClaim,
            SpectralTerminationV1::Completed,
        )
    };
    let first = SpectralTruthV1::new(&problem, duplicate_draft(vec![duplicate_a, duplicate_b]))
        .unwrap_err();
    let second = SpectralTruthV1::new(&problem, duplicate_draft(vec![duplicate_b, duplicate_a]))
        .unwrap_err();
    assert_eq!(first, second);
    assert!(
        first
            .issues()
            .contains(&SpectralTruthErrorV1::DuplicateCluster)
    );

    let partial_problem = validate_problem(default_spec(
        standard_axes(49, 3),
        Vec::new(),
        Vec::new(),
        SpectralOrderingV1::MagnitudeAscending {
            tie_break: ComplexTieBreakV1::RealThenImagThenLineage,
        },
        CompletenessScopeV1::Partial { requested: 2 },
    ))
    .unwrap();
    let duplicate_id = cluster_id(49);
    let enclosure = SpectralEnclosureV1::Real(FiniteIntervalV1::new(0.0, 1.0).unwrap());
    let multiplicity_one =
        exact_cluster(partial_problem.problem_id(), duplicate_id, enclosure, 1, 49);
    let multiplicity_two =
        exact_cluster(partial_problem.problem_id(), duplicate_id, enclosure, 2, 50);
    let result_set = spectral_result_set_receipt(&[multiplicity_one, multiplicity_two])
        .unwrap()
        .id();
    let status = PartialCoverageStatusV1::ClusterClosureOverrun {
        boundary_cluster: duplicate_id,
        preceding_algebraic: 1,
    };
    let coverage_witness = truth_witness(
        partial_problem.problem_id(),
        SpectralTruthPropositionV1::PartialCoverage {
            result_set,
            returned_algebraic: 3,
            status,
        },
        51,
    );
    let boundary_witness = truth_witness(
        partial_problem.problem_id(),
        SpectralTruthPropositionV1::PartialBoundaryClusterClosed {
            result_set,
            cluster: duplicate_id,
        },
        52,
    );
    let duplicate_partial_draft = |clusters| {
        SpectralTruthDraftV1::new(
            SpectralResultAuthorityV1::NoClaim,
            SpectralCoverageV1::Partial {
                returned_algebraic: 3,
                status,
                witness: coverage_witness,
            },
            clusters,
            ScopeBoundaryStateV1::Partial(PartialBoundaryStateV1::ClusterClosed {
                cluster: duplicate_id,
                witness: boundary_witness,
            }),
            SpectralTerminationV1::Completed,
        )
    };
    let one_first = SpectralTruthV1::new(
        &partial_problem,
        duplicate_partial_draft(vec![multiplicity_one, multiplicity_two]),
    )
    .unwrap_err();
    let two_first = SpectralTruthV1::new(
        &partial_problem,
        duplicate_partial_draft(vec![multiplicity_two, multiplicity_one]),
    )
    .unwrap_err();
    assert_eq!(one_first, two_first);
    assert!(
        one_first
            .issues()
            .contains(&SpectralTruthErrorV1::InvalidPartialCoverage)
    );
}

#[test]
fn exact_theorem_closure_is_support_and_definiteness_aware() {
    let axes = standard_axes(47, 4);
    let self_adjoint = structure_claim(
        axes,
        StructurePropertyV1::SelfAdjoint,
        StructureSupportV1::InnerProduct(axes.domain.id()),
        WitnessDispositionV1::Witnessed,
        0.0,
        47,
    );
    let nonnormal = structure_claim(
        axes,
        StructurePropertyV1::Nonnormal,
        StructureSupportV1::InnerProduct(axes.domain.id()),
        WitnessDispositionV1::Witnessed,
        0.0,
        48,
    );
    let report = validate_problem(default_spec(
        axes,
        vec![self_adjoint, nonnormal],
        Vec::new(),
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap_err();
    assert!(has_admission_issue(&report, |issue| matches!(
        issue,
        SpectralAdmissionIssueV1::StructureTheoremConflict {
            premise: StructurePropertyV1::SelfAdjoint,
            consequence: StructurePropertyV1::Normal,
            ..
        }
    )));

    let generalized_class = SpectralProblemClassV1::new(
        SpectralRepresentationV1::GeneralizedPencil,
        DescriptorRoleV1::Ordinary,
        SpectralOperatorOriginV1::Direct,
    );
    let generalized = crate::axes(
        48,
        generalized_class,
        SpectralScalarFieldV1::Complex,
        Dims::NONE,
        4,
    );
    let support = StructureSupportV1::InnerProduct(generalized.domain.id());
    let hermitian_definite = structure_claim(
        generalized,
        StructurePropertyV1::HermitianDefinitePencil,
        support,
        WitnessDispositionV1::Witnessed,
        0.0,
        48,
    );
    let contradicted_self_adjoint = structure_claim(
        generalized,
        StructurePropertyV1::SelfAdjoint,
        support,
        WitnessDispositionV1::Contradicted,
        0.0,
        49,
    );
    let report = validate_problem(default_spec(
        generalized,
        vec![hermitian_definite, contradicted_self_adjoint],
        Vec::new(),
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap_err();
    assert!(has_admission_issue(&report, |issue| matches!(
        issue,
        SpectralAdmissionIssueV1::StructureTheoremConflict {
            premise: StructurePropertyV1::HermitianDefinitePencil,
            consequence: StructurePropertyV1::SelfAdjoint,
            support: found,
        } if *found == support
    )));

    let raw = standard_axes(49, 4);
    let unknown_metric = SpectralMetricV1::new(
        SpectralMetricId::from_bytes([49; 32]),
        4,
        MetricDefinitenessV1::Unknown,
    );
    let indefinite_axes = Axes {
        domain: unknown_metric,
        codomain: unknown_metric,
        ..raw
    };
    let self_adjoint = structure_claim(
        indefinite_axes,
        StructurePropertyV1::SelfAdjoint,
        StructureSupportV1::InnerProduct(unknown_metric.id()),
        WitnessDispositionV1::Witnessed,
        0.0,
        49,
    );
    let nonreal = structure_claim(
        indefinite_axes,
        StructurePropertyV1::RealSpectrum,
        StructureSupportV1::FormFree,
        WitnessDispositionV1::Contradicted,
        0.0,
        50,
    );
    assert!(
        validate_problem(default_spec(
            indefinite_axes,
            vec![self_adjoint, nonreal],
            Vec::new(),
            SpectralOrderingV1::SetValued,
            CompletenessScopeV1::CandidateOnly,
        ))
        .is_ok(),
        "indefinite/unknown-metric self-adjointness must not forge a real-spectrum theorem"
    );

    let singular_id = SpectralMetricId::from_bytes([50; 32]);
    let singular_witness = admit_receipt(
        metric_proposition_receipt(
            singular_id,
            4,
            MetricDefinitenessPropositionV1::Singular { rank: 3 },
        )
        .unwrap(),
        51,
    );
    let singular_metric = SpectralMetricV1::new(
        singular_id,
        4,
        MetricDefinitenessV1::Singular {
            rank: 3,
            witness: singular_witness,
        },
    );
    let singular_axes = Axes {
        domain: singular_metric,
        codomain: singular_metric,
        ..raw
    };
    let singular_self_adjoint = structure_claim(
        singular_axes,
        StructurePropertyV1::SelfAdjoint,
        StructureSupportV1::InnerProduct(singular_id),
        WitnessDispositionV1::Witnessed,
        0.0,
        52,
    );
    let report = validate_problem(default_spec(
        singular_axes,
        vec![singular_self_adjoint],
        Vec::new(),
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap_err();
    assert!(has_admission_issue(&report, |issue| matches!(
        issue,
        SpectralAdmissionIssueV1::InvalidStructureSupport {
            property: StructurePropertyV1::SelfAdjoint,
            support: StructureSupportV1::InnerProduct(found),
        } if *found == singular_id
    )));
}

#[test]
fn structure_tolerance_claims_respect_nested_defect_sets() {
    let axes = standard_axes(50, 3);
    let support = StructureSupportV1::InnerProduct(axes.domain.id());
    let loose_witness = structure_claim_with_seed(
        axes,
        StructurePropertyV1::Normal,
        support,
        WitnessDispositionV1::Witnessed,
        1.0e-3,
        norm_id(50),
        50,
    );
    let tight_contradiction = structure_claim_with_seed(
        axes,
        StructurePropertyV1::Normal,
        support,
        WitnessDispositionV1::Contradicted,
        1.0e-4,
        norm_id(50),
        51,
    );
    assert!(
        validate_problem(default_spec(
            axes,
            vec![loose_witness, tight_contradiction],
            Vec::new(),
            SpectralOrderingV1::SetValued,
            CompletenessScopeV1::CandidateOnly,
        ))
        .is_ok(),
        "a property may hold at a loose tolerance while failing at a tighter one"
    );

    let tight_witness = structure_claim_with_seed(
        axes,
        StructurePropertyV1::Normal,
        support,
        WitnessDispositionV1::Witnessed,
        1.0e-4,
        norm_id(50),
        52,
    );
    let loose_contradiction = structure_claim_with_seed(
        axes,
        StructurePropertyV1::Normal,
        support,
        WitnessDispositionV1::Contradicted,
        1.0e-3,
        norm_id(50),
        53,
    );
    let report = validate_problem(default_spec(
        axes,
        vec![tight_witness, loose_contradiction],
        Vec::new(),
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap_err();
    assert!(has_admission_issue(&report, |issue| matches!(
        issue,
        SpectralAdmissionIssueV1::ContradictoryStructure {
            property: StructurePropertyV1::Normal,
            ..
        }
    )));

    let generalized_class = SpectralProblemClassV1::new(
        SpectralRepresentationV1::GeneralizedPencil,
        DescriptorRoleV1::Ordinary,
        SpectralOperatorOriginV1::Direct,
    );
    let generalized = crate::axes(
        54,
        generalized_class,
        SpectralScalarFieldV1::Complex,
        Dims::NONE,
        3,
    );
    let support = StructureSupportV1::InnerProduct(generalized.domain.id());
    let exact_hdp = structure_claim_with_seed(
        generalized,
        StructurePropertyV1::HermitianDefinitePencil,
        support,
        WitnessDispositionV1::Witnessed,
        0.0,
        norm_id(54),
        54,
    );
    let cross_norm_contradiction = structure_claim_with_seed(
        generalized,
        StructurePropertyV1::HermitianDefinitePencil,
        support,
        WitnessDispositionV1::Contradicted,
        1.0e-3,
        norm_id(55),
        55,
    );
    let report = validate_problem(default_spec(
        generalized,
        vec![exact_hdp, cross_norm_contradiction],
        Vec::new(),
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap_err();
    assert!(has_admission_issue(&report, |issue| matches!(
        issue,
        SpectralAdmissionIssueV1::ContradictoryStructure {
            property: StructurePropertyV1::HermitianDefinitePencil,
            ..
        }
    )));
}

#[test]
fn equation_specific_structure_and_ordinary_finite_semantics_fail_closed() {
    let standard = standard_axes(51, 2);
    let palindromic = structure_claim(
        standard,
        StructurePropertyV1::Palindromic {
            parity: PalindromicParityV1::Palindromic,
            involution: PolynomialInvolutionV1::ConjugateTranspose,
        },
        StructureSupportV1::FormFree,
        WitnessDispositionV1::Witnessed,
        0.0,
        51,
    );
    let report = validate_problem(default_spec(
        standard,
        vec![palindromic],
        Vec::new(),
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap_err();
    assert!(has_admission_issue(&report, |issue| matches!(
        issue,
        SpectralAdmissionIssueV1::StructureRepresentationMismatch {
            property: StructurePropertyV1::Palindromic { .. },
            representation: SpectralRepresentationV1::StandardLinear,
        }
    )));

    let class = SpectralProblemClassV1::new(
        SpectralRepresentationV1::MatrixPolynomial { grade: 2 },
        DescriptorRoleV1::Ordinary,
        SpectralOperatorOriginV1::Direct,
    );
    let polynomial = axes(52, class, SpectralScalarFieldV1::Complex, Dims::NONE, 2);
    let regular = regularity_claim(
        polynomial,
        RegularityClassV1::RegularPolynomial { grade: 2 },
        WitnessDispositionV1::Witnessed,
        52,
    );
    let report = validate_problem(default_spec(
        polynomial,
        Vec::new(),
        vec![regular],
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap_err();
    assert!(has_admission_issue(&report, |issue| matches!(
        issue,
        SpectralAdmissionIssueV1::OrdinaryFiniteSpectrumWitnessRequired {
            required: RegularityClassV1::InvertiblePolynomialLeadingCoefficient { grade: 2 }
        }
    )));

    let leading = regularity_claim(
        polynomial,
        RegularityClassV1::InvertiblePolynomialLeadingCoefficient { grade: 2 },
        WitnessDispositionV1::Witnessed,
        53,
    );
    let leading_problem = validate_problem(default_spec(
        polynomial,
        Vec::new(),
        vec![leading],
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap();
    assert_eq!(leading_problem.known_algebraic_cardinality(), Some(4));
    assert!(assess_method_class(&leading_problem, SpectralMethodClassV1::PolynomialKrylov).is_ok());
    let contradicted_regular = regularity_claim(
        polynomial,
        RegularityClassV1::RegularPolynomial { grade: 2 },
        WitnessDispositionV1::Contradicted,
        54,
    );
    let report = validate_problem(default_spec(
        polynomial,
        Vec::new(),
        vec![leading, contradicted_regular],
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap_err();
    assert!(has_admission_issue(&report, |issue| matches!(
        issue,
        SpectralAdmissionIssueV1::RegularityTheoremConflict {
            premise: RegularityClassV1::InvertiblePolynomialLeadingCoefficient { grade: 2 },
            consequence: RegularityClassV1::RegularPolynomial { grade: 2 },
        }
    )));

    let pencil_class = SpectralProblemClassV1::new(
        SpectralRepresentationV1::GeneralizedPencil,
        DescriptorRoleV1::Ordinary,
        SpectralOperatorOriginV1::Direct,
    );
    let raw_pencil = axes(
        53,
        pencil_class,
        SpectralScalarFieldV1::Complex,
        Dims::NONE,
        2,
    );
    let hdp = structure_claim(
        raw_pencil,
        StructurePropertyV1::HermitianDefinitePencil,
        StructureSupportV1::InnerProduct(raw_pencil.domain.id()),
        WitnessDispositionV1::Witnessed,
        0.0,
        55,
    );
    let contradicted_regular = regularity_claim(
        raw_pencil,
        RegularityClassV1::RegularPencil,
        WitnessDispositionV1::Contradicted,
        56,
    );
    let report = validate_problem(default_spec(
        raw_pencil,
        vec![hdp],
        vec![contradicted_regular],
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap_err();
    assert!(has_admission_issue(&report, |issue| matches!(
        issue,
        SpectralAdmissionIssueV1::StructureRegularityTheoremConflict {
            premise: StructurePropertyV1::HermitianDefinitePencil,
            consequence: RegularityClassV1::RegularPencil,
            ..
        }
    )));

    let unknown_metric = SpectralMetricV1::new(
        SpectralMetricId::from_bytes([53; 32]),
        2,
        MetricDefinitenessV1::Unknown,
    );
    let pencil = Axes {
        domain: unknown_metric,
        codomain: unknown_metric,
        ..raw_pencil
    };
    let hdp = structure_claim(
        pencil,
        StructurePropertyV1::HermitianDefinitePencil,
        StructureSupportV1::InnerProduct(unknown_metric.id()),
        WitnessDispositionV1::Witnessed,
        0.0,
        54,
    );
    let report = validate_problem(default_spec(
        pencil,
        vec![hdp],
        Vec::new(),
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap_err();
    assert!(has_admission_issue(&report, |issue| matches!(
        issue,
        SpectralAdmissionIssueV1::InvalidStructureSupport {
            property: StructurePropertyV1::HermitianDefinitePencil,
            ..
        }
    )));
}

#[test]
fn real_spectrum_truth_and_total_cardinality_do_not_depend_on_ordering() {
    let standard = standard_axes(54, 2);
    let report = validate_problem(default_spec(
        standard,
        Vec::new(),
        Vec::new(),
        SpectralOrderingV1::MagnitudeAscending {
            tie_break: ComplexTieBreakV1::RealThenImagThenLineage,
        },
        CompletenessScopeV1::Partial { requested: 3 },
    ))
    .unwrap_err();
    assert!(has_admission_issue(&report, |issue| matches!(
        issue,
        SpectralAdmissionIssueV1::DimensionMismatch { left: 3, right: 2 }
    )));

    let real_spectrum = structure_claim(
        standard,
        StructurePropertyV1::RealSpectrum,
        StructureSupportV1::FormFree,
        WitnessDispositionV1::Witnessed,
        0.0,
        54,
    );
    let problem = validate_problem(default_spec(
        standard,
        vec![real_spectrum],
        Vec::new(),
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap();
    let off_axis = SpectralEnclosureV1::ComplexBox {
        real: FiniteIntervalV1::new(0.0, 1.0).unwrap(),
        imag: FiniteIntervalV1::new(1.0, 2.0).unwrap(),
    };
    let id = cluster_id(54);
    let localization = SpectralLocalizationV1::enclosed(
        off_axis,
        truth_witness(
            problem.problem_id(),
            SpectralTruthPropositionV1::ClusterLocalization {
                cluster: id,
                authority: LocalizationAuthorityV1::Enclosed,
                enclosure: off_axis,
            },
            55,
        ),
    );
    let cluster = SpectralClusterV1::new(
        id,
        localization,
        MultiplicityClaimV1::Unknown,
        MultiplicityClaimV1::Unknown,
        InternalClusterStateV1::NoClaim,
    )
    .unwrap();
    let report = SpectralTruthV1::new(
        &problem,
        SpectralTruthDraftV1::new(
            SpectralResultAuthorityV1::NoClaim,
            SpectralCoverageV1::Candidates,
            vec![cluster],
            ScopeBoundaryStateV1::NoClaim,
            SpectralTerminationV1::Completed,
        ),
    )
    .unwrap_err();
    assert!(has_truth_issue(&report, |issue| matches!(
        issue,
        SpectralTruthErrorV1::RealSpectrumEnclosureConflict
    )));

    let descriptor_class = SpectralProblemClassV1::new(
        SpectralRepresentationV1::GeneralizedPencil,
        DescriptorRoleV1::Descriptor {
            infinity_policy: InfiniteEigenvaluePolicyV1::IncludeProjective,
        },
        SpectralOperatorOriginV1::Direct,
    );
    let descriptor = axes(
        57,
        descriptor_class,
        SpectralScalarFieldV1::Complex,
        Dims::NONE,
        2,
    );
    let projective_problem = validate_problem(default_spec(
        descriptor,
        vec![structure_claim(
            descriptor,
            StructurePropertyV1::RealSpectrum,
            StructureSupportV1::FormFree,
            WitnessDispositionV1::Witnessed,
            0.0,
            57,
        )],
        vec![
            regularity_claim(
                descriptor,
                RegularityClassV1::RegularPencil,
                WitnessDispositionV1::Witnessed,
                58,
            ),
            regularity_claim(
                descriptor,
                RegularityClassV1::RegularDescriptor,
                WitnessDispositionV1::Witnessed,
                59,
            ),
        ],
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap();
    assert!(!projective_problem.projective_infinity_is_excluded());
    let projective_id = cluster_id(57);
    let projective = SpectralClusterV1::new(
        projective_id,
        SpectralLocalizationV1::enclosed(
            SpectralEnclosureV1::ProjectiveInfinity,
            truth_witness(
                projective_problem.problem_id(),
                SpectralTruthPropositionV1::ClusterLocalization {
                    cluster: projective_id,
                    authority: LocalizationAuthorityV1::Enclosed,
                    enclosure: SpectralEnclosureV1::ProjectiveInfinity,
                },
                60,
            ),
        ),
        MultiplicityClaimV1::Unknown,
        MultiplicityClaimV1::Unknown,
        InternalClusterStateV1::NoClaim,
    )
    .unwrap();
    assert!(
        SpectralTruthV1::new(
            &projective_problem,
            SpectralTruthDraftV1::new(
                SpectralResultAuthorityV1::NoClaim,
                SpectralCoverageV1::Candidates,
                vec![projective],
                ScopeBoundaryStateV1::NoClaim,
                SpectralTerminationV1::Completed,
            ),
        )
        .is_ok(),
        "projective infinity lies on the extended real line unless invertibility excludes it"
    );

    let no_structure = validate_problem(default_spec(
        standard,
        Vec::new(),
        Vec::new(),
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap();
    let clusters = (0..3)
        .map(|index| {
            exact_cluster(
                no_structure.problem_id(),
                indexed_cluster_id(index),
                SpectralEnclosureV1::Real(
                    FiniteIntervalV1::new(index as f64, index as f64).unwrap(),
                ),
                1,
                56_u8.wrapping_add(index as u8),
            )
        })
        .collect();
    let report = SpectralTruthV1::new(
        &no_structure,
        SpectralTruthDraftV1::new(
            SpectralResultAuthorityV1::NoClaim,
            SpectralCoverageV1::Candidates,
            clusters,
            ScopeBoundaryStateV1::NoClaim,
            SpectralTerminationV1::Completed,
        ),
    )
    .unwrap_err();
    assert!(has_truth_issue(&report, |issue| matches!(
        issue,
        SpectralTruthErrorV1::CoverageCardinalityMismatch
    )));

    let id = cluster_id(59);
    let enclosure = SpectralEnclosureV1::Real(FiniteIntervalV1::new(0.0, 1.0).unwrap());
    let lower_bound = MultiplicityClaimV1::LowerBound {
        value: 3,
        witness: truth_witness(
            no_structure.problem_id(),
            SpectralTruthPropositionV1::Multiplicity {
                cluster: id,
                enclosure,
                kind: MultiplicityKindV1::Algebraic,
                assertion: MultiplicityAssertionV1::LowerBound,
                lower: 3,
                upper: None,
            },
            59,
        ),
    };
    let lower_bound_cluster = SpectralClusterV1::new(
        id,
        SpectralLocalizationV1::candidate(enclosure),
        lower_bound,
        MultiplicityClaimV1::Unknown,
        InternalClusterStateV1::NoClaim,
    )
    .unwrap();
    let report = SpectralTruthV1::new(
        &no_structure,
        SpectralTruthDraftV1::new(
            SpectralResultAuthorityV1::NoClaim,
            SpectralCoverageV1::Candidates,
            vec![lower_bound_cluster],
            ScopeBoundaryStateV1::NoClaim,
            SpectralTerminationV1::Completed,
        ),
    )
    .unwrap_err();
    assert!(has_truth_issue(&report, |issue| matches!(
        issue,
        SpectralTruthErrorV1::CoverageCardinalityMismatch
    )));

    let id = cluster_id(60);
    let geometric = MultiplicityClaimV1::LowerBound {
        value: 3,
        witness: truth_witness(
            no_structure.problem_id(),
            SpectralTruthPropositionV1::Multiplicity {
                cluster: id,
                enclosure,
                kind: MultiplicityKindV1::Geometric,
                assertion: MultiplicityAssertionV1::LowerBound,
                lower: 3,
                upper: None,
            },
            60,
        ),
    };
    let geometric_cluster = SpectralClusterV1::new(
        id,
        SpectralLocalizationV1::candidate(enclosure),
        MultiplicityClaimV1::Unknown,
        geometric,
        InternalClusterStateV1::NoClaim,
    )
    .unwrap();
    let report = SpectralTruthV1::new(
        &no_structure,
        SpectralTruthDraftV1::new(
            SpectralResultAuthorityV1::NoClaim,
            SpectralCoverageV1::Candidates,
            vec![geometric_cluster],
            ScopeBoundaryStateV1::NoClaim,
            SpectralTerminationV1::Completed,
        ),
    )
    .unwrap_err();
    assert!(has_truth_issue(&report, |issue| matches!(
        issue,
        SpectralTruthErrorV1::GeometricCapacityExceeded {
            minimum: 3,
            dimension: 2,
            ..
        }
    )));

    let geometric_clusters = (0..2)
        .map(|index| {
            let id = indexed_cluster_id(100 + index);
            let geometric = MultiplicityClaimV1::LowerBound {
                value: 2,
                witness: truth_witness(
                    no_structure.problem_id(),
                    SpectralTruthPropositionV1::Multiplicity {
                        cluster: id,
                        enclosure,
                        kind: MultiplicityKindV1::Geometric,
                        assertion: MultiplicityAssertionV1::LowerBound,
                        lower: 2,
                        upper: None,
                    },
                    61_u8.wrapping_add(index as u8),
                ),
            };
            SpectralClusterV1::new(
                id,
                SpectralLocalizationV1::candidate(enclosure),
                MultiplicityClaimV1::Unknown,
                geometric,
                InternalClusterStateV1::NoClaim,
            )
            .unwrap()
        })
        .collect();
    let report = SpectralTruthV1::new(
        &no_structure,
        SpectralTruthDraftV1::new(
            SpectralResultAuthorityV1::NoClaim,
            SpectralCoverageV1::Candidates,
            geometric_clusters,
            ScopeBoundaryStateV1::NoClaim,
            SpectralTerminationV1::Completed,
        ),
    )
    .unwrap_err();
    assert!(has_truth_issue(&report, |issue| matches!(
        issue,
        SpectralTruthErrorV1::CoverageCardinalityMismatch
    )));

    let polynomial_class = SpectralProblemClassV1::new(
        SpectralRepresentationV1::MatrixPolynomial { grade: 2 },
        DescriptorRoleV1::Ordinary,
        SpectralOperatorOriginV1::Direct,
    );
    let polynomial_axes = axes(
        62,
        polynomial_class,
        SpectralScalarFieldV1::Complex,
        Dims::NONE,
        4,
    );
    let quotient = quotient_map_id(62);
    let quotient_witness = admit_receipt(
        gauge_proposition_receipt(
            polynomial_axes.subject,
            polynomial_axes.scalar,
            polynomial_axes.class,
            polynomial_axes.scaling,
            polynomial_axes.domain,
            polynomial_axes.codomain,
            GaugePropositionV1::Quotiented {
                nullity: 3,
                quotient,
            },
        )
        .unwrap(),
        62,
    );
    let leading = regularity_claim(
        polynomial_axes,
        RegularityClassV1::InvertiblePolynomialLeadingCoefficient { grade: 2 },
        WitnessDispositionV1::Witnessed,
        63,
    );
    let quotient_problem = validate_problem(problem_spec(
        polynomial_axes,
        Vec::new(),
        vec![leading],
        SpectralSpaceContextV1::new(
            polynomial_axes.domain,
            polynomial_axes.codomain,
            GaugeConventionV1::Quotiented {
                nullity: 3,
                quotient,
                witness: quotient_witness,
            },
            ZeroPaddingConventionV1::Unknown,
        ),
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap();
    assert_eq!(quotient_problem.known_algebraic_cardinality(), Some(8));
    let quotient_cluster_id = cluster_id(62);
    let quotient_enclosure = SpectralEnclosureV1::Real(FiniteIntervalV1::new(0.0, 1.0).unwrap());
    let geometric = MultiplicityClaimV1::LowerBound {
        value: 5,
        witness: truth_witness(
            quotient_problem.problem_id(),
            SpectralTruthPropositionV1::Multiplicity {
                cluster: quotient_cluster_id,
                enclosure: quotient_enclosure,
                kind: MultiplicityKindV1::Geometric,
                assertion: MultiplicityAssertionV1::LowerBound,
                lower: 5,
                upper: None,
            },
            64,
        ),
    };
    let quotient_cluster = SpectralClusterV1::new(
        quotient_cluster_id,
        SpectralLocalizationV1::candidate(quotient_enclosure),
        MultiplicityClaimV1::Unknown,
        geometric,
        InternalClusterStateV1::NoClaim,
    )
    .unwrap();
    let report = SpectralTruthV1::new(
        &quotient_problem,
        SpectralTruthDraftV1::new(
            SpectralResultAuthorityV1::NoClaim,
            SpectralCoverageV1::Candidates,
            vec![quotient_cluster],
            ScopeBoundaryStateV1::NoClaim,
            SpectralTerminationV1::Completed,
        ),
    )
    .unwrap_err();
    assert!(has_truth_issue(&report, |issue| matches!(
        issue,
        SpectralTruthErrorV1::GeometricCapacityExceeded {
            minimum: 5,
            dimension: 4,
            ..
        }
    )));
}

#[test]
fn public_truth_receipt_caps_are_reachable_at_the_advertised_limit() {
    let problem = validate_problem(default_spec(
        standard_axes(60, 3 * u32::try_from(MAX_SPECTRAL_CLUSTERS_V1).unwrap()),
        Vec::new(),
        Vec::new(),
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap();
    let enclosure = SpectralEnclosureV1::ComplexBox {
        real: FiniteIntervalV1::new(-1.0, 1.0).unwrap(),
        imag: FiniteIntervalV1::new(-1.0, 1.0).unwrap(),
    };
    let separation = PositiveFiniteV1::new(1.0).unwrap();
    let norm = norm_id(60);
    let clusters: Vec<_> = (0..MAX_SPECTRAL_CLUSTERS_V1)
        .map(|index| {
            let id = indexed_cluster_id(index);
            let algebraic_statement = MultiplicityStatementV1::Bounds { lower: 2, upper: 3 };
            let geometric_statement = MultiplicityStatementV1::Bounds { lower: 1, upper: 2 };
            let algebraic = MultiplicityClaimV1::Bounds {
                lower: 2,
                upper: 3,
                witness: truth_witness(
                    problem.problem_id(),
                    SpectralTruthPropositionV1::Multiplicity {
                        cluster: id,
                        enclosure,
                        kind: MultiplicityKindV1::Algebraic,
                        assertion: MultiplicityAssertionV1::Bounds,
                        lower: 2,
                        upper: Some(3),
                    },
                    60,
                ),
            };
            let geometric = MultiplicityClaimV1::Bounds {
                lower: 1,
                upper: 2,
                witness: truth_witness(
                    problem.problem_id(),
                    SpectralTruthPropositionV1::Multiplicity {
                        cluster: id,
                        enclosure,
                        kind: MultiplicityKindV1::Geometric,
                        assertion: MultiplicityAssertionV1::Bounds,
                        lower: 1,
                        upper: Some(2),
                    },
                    60,
                ),
            };
            let internal = InternalClusterStateV1::Resolved {
                lower: separation,
                norm,
                witness: truth_witness(
                    problem.problem_id(),
                    SpectralTruthPropositionV1::InternalResolution {
                        cluster: id,
                        enclosure,
                        algebraic: algebraic_statement,
                        geometric: geometric_statement,
                        lower: separation,
                        norm,
                    },
                    60,
                ),
            };
            SpectralClusterV1::new(
                id,
                SpectralLocalizationV1::candidate(enclosure),
                algebraic,
                geometric,
                internal,
            )
            .unwrap()
        })
        .collect();
    spectral_result_set_receipt(&clusters).unwrap();

    let result_set = spectral_result_set_receipt(&[]).unwrap().id();
    let references: Vec<_> = (0..MAX_REGION_BOUNDARY_REFERENCES_V1)
        .map(indexed_cluster_id)
        .collect();
    truth_proposition_receipt(
        problem.problem_id(),
        &SpectralTruthPropositionV1::RegionBoundaryIntersections {
            result_set,
            included: references,
            excluded_algebraic: 0,
        },
    )
    .unwrap();
}

#[test]
fn multiplicity_and_projective_truth_cannot_be_replayed_or_duplicated() {
    let problem = validate_problem(default_spec(
        standard_axes(61, 2),
        Vec::new(),
        Vec::new(),
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap();
    let original_id = cluster_id(61);
    let original_enclosure = SpectralEnclosureV1::Real(FiniteIntervalV1::new(0.0, 1.0).unwrap());
    let stale_multiplicity_witness = truth_witness(
        problem.problem_id(),
        SpectralTruthPropositionV1::Multiplicity {
            cluster: original_id,
            enclosure: original_enclosure,
            kind: MultiplicityKindV1::Algebraic,
            assertion: MultiplicityAssertionV1::Exact,
            lower: 1,
            upper: Some(1),
        },
        61,
    );
    let rebound = SpectralClusterV1::new(
        original_id,
        SpectralLocalizationV1::candidate(SpectralEnclosureV1::Real(
            FiniteIntervalV1::new(2.0, 3.0).unwrap(),
        )),
        MultiplicityClaimV1::Exact {
            value: 1,
            witness: stale_multiplicity_witness,
        },
        MultiplicityClaimV1::Unknown,
        InternalClusterStateV1::Simple,
    )
    .unwrap();
    let report = SpectralTruthV1::new(
        &problem,
        SpectralTruthDraftV1::new(
            SpectralResultAuthorityV1::NoClaim,
            SpectralCoverageV1::Candidates,
            vec![rebound],
            ScopeBoundaryStateV1::NoClaim,
            SpectralTerminationV1::Completed,
        ),
    )
    .unwrap_err();
    assert!(has_truth_issue(&report, |issue| matches!(
        issue,
        SpectralTruthErrorV1::WitnessPropositionMismatch { .. }
    )));

    let foreign_localization = SpectralClusterV1::new(
        original_id,
        SpectralLocalizationV1::enclosed(original_enclosure, stale_multiplicity_witness),
        MultiplicityClaimV1::Unknown,
        MultiplicityClaimV1::Unknown,
        InternalClusterStateV1::NoClaim,
    )
    .unwrap();
    let report = SpectralTruthV1::new(
        &problem,
        SpectralTruthDraftV1::new(
            SpectralResultAuthorityV1::NoClaim,
            SpectralCoverageV1::Candidates,
            vec![foreign_localization],
            ScopeBoundaryStateV1::NoClaim,
            SpectralTerminationV1::Completed,
        ),
    )
    .unwrap_err();
    assert!(has_truth_issue(&report, |issue| matches!(
        issue,
        SpectralTruthErrorV1::WitnessPropositionMismatch { .. }
    )));

    let descriptor = descriptor_full_problem(62, InfiniteEigenvaluePolicyV1::IncludeProjective);
    let make_projective = |id| {
        SpectralClusterV1::new(
            id,
            SpectralLocalizationV1::candidate(SpectralEnclosureV1::ProjectiveInfinity),
            MultiplicityClaimV1::Unknown,
            MultiplicityClaimV1::Unknown,
            InternalClusterStateV1::NoClaim,
        )
        .unwrap()
    };
    let report = SpectralTruthV1::new(
        &descriptor,
        SpectralTruthDraftV1::new(
            SpectralResultAuthorityV1::NoClaim,
            SpectralCoverageV1::Candidates,
            vec![
                make_projective(cluster_id(62)),
                make_projective(cluster_id(63)),
            ],
            ScopeBoundaryStateV1::NoClaim,
            SpectralTerminationV1::Completed,
        ),
    )
    .unwrap_err();
    assert!(has_truth_issue(&report, |issue| matches!(
        issue,
        SpectralTruthErrorV1::MultipleProjectiveClusters
    )));
}

#[test]
fn region_boundary_policy_is_enforced_even_for_candidate_coverage() {
    let region = region_id(64);
    let axes = standard_axes(64, 2);
    let problem = validate_problem(default_spec(
        axes,
        Vec::new(),
        Vec::new(),
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::Region {
            region,
            boundary: RegionBoundaryPolicyV1::Closed,
        },
    ))
    .unwrap();
    let clusters = Vec::new();
    let result_set = spectral_result_set_receipt(&clusters).unwrap().id();
    let witness = truth_witness(
        problem.problem_id(),
        SpectralTruthPropositionV1::RegionBoundaryIntersections {
            result_set,
            included: Vec::new(),
            excluded_algebraic: 1,
        },
        64,
    );
    let report = SpectralTruthV1::new(
        &problem,
        SpectralTruthDraftV1::new(
            SpectralResultAuthorityV1::NoClaim,
            SpectralCoverageV1::Candidates,
            clusters,
            ScopeBoundaryStateV1::Region(RegionBoundaryStateV1::IntersectionsResolved {
                included: Vec::new(),
                excluded_algebraic: 1,
                witness,
            }),
            SpectralTerminationV1::Completed,
        ),
    )
    .unwrap_err();
    assert!(has_truth_issue(&report, |issue| matches!(
        issue,
        SpectralTruthErrorV1::BoundaryScopeMismatch
    )));

    let open_problem = validate_problem(default_spec(
        axes,
        Vec::new(),
        Vec::new(),
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::Region {
            region,
            boundary: RegionBoundaryPolicyV1::Open,
        },
    ))
    .unwrap();
    let cluster = candidate_cluster(cluster_id(64));
    let open_clusters = vec![cluster];
    let open_set = spectral_result_set_receipt(&open_clusters).unwrap().id();
    let open_witness = truth_witness(
        open_problem.problem_id(),
        SpectralTruthPropositionV1::RegionBoundaryIntersections {
            result_set: open_set,
            included: vec![cluster.id()],
            excluded_algebraic: 0,
        },
        65,
    );
    let report = SpectralTruthV1::new(
        &open_problem,
        SpectralTruthDraftV1::new(
            SpectralResultAuthorityV1::NoClaim,
            SpectralCoverageV1::Candidates,
            open_clusters,
            ScopeBoundaryStateV1::Region(RegionBoundaryStateV1::IntersectionsResolved {
                included: vec![cluster.id()],
                excluded_algebraic: 0,
                witness: open_witness,
            }),
            SpectralTerminationV1::Completed,
        ),
    )
    .unwrap_err();
    assert!(has_truth_issue(&report, |issue| matches!(
        issue,
        SpectralTruthErrorV1::BoundaryScopeMismatch
    )));

    let empty_set = spectral_result_set_receipt(&[]).unwrap().id();
    let excessive_exclusion = truth_witness(
        open_problem.problem_id(),
        SpectralTruthPropositionV1::RegionBoundaryIntersections {
            result_set: empty_set,
            included: Vec::new(),
            excluded_algebraic: 3,
        },
        66,
    );
    let report = SpectralTruthV1::new(
        &open_problem,
        SpectralTruthDraftV1::new(
            SpectralResultAuthorityV1::NoClaim,
            SpectralCoverageV1::Candidates,
            Vec::new(),
            ScopeBoundaryStateV1::Region(RegionBoundaryStateV1::IntersectionsResolved {
                included: Vec::new(),
                excluded_algebraic: 3,
                witness: excessive_exclusion,
            }),
            SpectralTerminationV1::Completed,
        ),
    )
    .unwrap_err();
    assert!(has_truth_issue(&report, |issue| matches!(
        issue,
        SpectralTruthErrorV1::CoverageCardinalityMismatch
    )));

    let complete_cluster = exact_cluster(
        open_problem.problem_id(),
        cluster_id(66),
        SpectralEnclosureV1::Real(FiniteIntervalV1::new(0.0, 1.0).unwrap()),
        2,
        67,
    );
    let complete_set = spectral_result_set_receipt(&[complete_cluster])
        .unwrap()
        .id();
    let complete_boundary = truth_witness(
        open_problem.problem_id(),
        SpectralTruthPropositionV1::RegionBoundaryIntersections {
            result_set: complete_set,
            included: Vec::new(),
            excluded_algebraic: 1,
        },
        68,
    );
    let complete_coverage = truth_witness(
        open_problem.problem_id(),
        SpectralTruthPropositionV1::RegionCompleteness {
            result_set: complete_set,
            algebraic_cardinality: 2,
        },
        69,
    );
    let report = SpectralTruthV1::new(
        &open_problem,
        SpectralTruthDraftV1::new(
            SpectralResultAuthorityV1::NoClaim,
            SpectralCoverageV1::RegionComplete {
                algebraic_cardinality: 2,
                witness: complete_coverage,
            },
            vec![complete_cluster],
            ScopeBoundaryStateV1::Region(RegionBoundaryStateV1::IntersectionsResolved {
                included: Vec::new(),
                excluded_algebraic: 1,
                witness: complete_boundary,
            }),
            SpectralTerminationV1::Completed,
        ),
    )
    .unwrap_err();
    assert!(has_truth_issue(&report, |issue| matches!(
        issue,
        SpectralTruthErrorV1::CoverageCardinalityMismatch
    )));

    let refuse_problem = validate_problem(default_spec(
        axes,
        Vec::new(),
        Vec::new(),
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::Region {
            region,
            boundary: RegionBoundaryPolicyV1::RefuseIntersection,
        },
    ))
    .unwrap();
    let empty_set = spectral_result_set_receipt(&[]).unwrap().id();
    let refuse_witness = truth_witness(
        refuse_problem.problem_id(),
        SpectralTruthPropositionV1::RegionBoundaryIntersections {
            result_set: empty_set,
            included: Vec::new(),
            excluded_algebraic: 0,
        },
        66,
    );
    let report = SpectralTruthV1::new(
        &refuse_problem,
        SpectralTruthDraftV1::new(
            SpectralResultAuthorityV1::NoClaim,
            SpectralCoverageV1::Candidates,
            Vec::new(),
            ScopeBoundaryStateV1::Region(RegionBoundaryStateV1::IntersectionsResolved {
                included: Vec::new(),
                excluded_algebraic: 0,
                witness: refuse_witness,
            }),
            SpectralTerminationV1::Completed,
        ),
    )
    .unwrap_err();
    assert!(has_truth_issue(&report, |issue| matches!(
        issue,
        SpectralTruthErrorV1::BoundaryScopeMismatch
    )));
}

#[test]
fn signed_zero_is_canonical_in_sealed_problem_and_truth_objects() {
    let axes = standard_axes(65, 2);
    let make_problem = |zero: f64| {
        validate_problem(default_spec(
            axes,
            Vec::new(),
            Vec::new(),
            SpectralOrderingV1::NearestShift {
                real: zero,
                imag: -zero,
                tie_break: ComplexTieBreakV1::ImagThenRealThenLineage,
            },
            CompletenessScopeV1::Partial { requested: 1 },
        ))
        .unwrap()
    };
    let positive = make_problem(0.0);
    let negative = make_problem(-0.0);
    assert_eq!(positive.problem_id(), negative.problem_id());
    assert_eq!(positive.spec(), negative.spec());
    let SpectralOrderingV1::NearestShift { real, imag, .. } = negative.spec().ordering() else {
        panic!("fixture must retain nearest-shift ordering");
    };
    assert_eq!(real.to_bits(), 0.0f64.to_bits());
    assert_eq!(imag.to_bits(), 0.0f64.to_bits());

    let truth_problem = validate_problem(default_spec(
        axes,
        Vec::new(),
        Vec::new(),
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap();
    let result_set = spectral_result_set_receipt(&[]).unwrap().id();
    let norm = norm_id(65);
    let witness = truth_witness(
        truth_problem.problem_id(),
        SpectralTruthPropositionV1::ResultResidualBound {
            result_set,
            upper: 0.0,
            norm,
        },
        65,
    );
    let truth = SpectralTruthV1::new(
        &truth_problem,
        SpectralTruthDraftV1::new(
            SpectralResultAuthorityV1::ResidualBounded {
                upper: -0.0,
                norm,
                witness,
            },
            SpectralCoverageV1::Candidates,
            Vec::new(),
            ScopeBoundaryStateV1::NoClaim,
            SpectralTerminationV1::Completed,
        ),
    )
    .unwrap();
    let SpectralResultAuthorityV1::ResidualBounded { upper, .. } = truth.authority() else {
        panic!("fixture must retain residual authority");
    };
    assert_eq!(upper.to_bits(), 0.0f64.to_bits());
}

#[test]
fn projective_partial_prefix_requires_chart_and_infinity_placement() {
    let class = SpectralProblemClassV1::new(
        SpectralRepresentationV1::GeneralizedPencil,
        DescriptorRoleV1::Descriptor {
            infinity_policy: InfiniteEigenvaluePolicyV1::IncludeProjective,
        },
        SpectralOperatorOriginV1::Direct,
    );
    let projective_axes = axes(79, class, SpectralScalarFieldV1::Complex, Dims::NONE, 3);
    let report = validate_problem(default_spec(
        projective_axes,
        Vec::new(),
        Vec::new(),
        SpectralOrderingV1::MagnitudeAscending {
            tie_break: ComplexTieBreakV1::RealThenImagThenLineage,
        },
        CompletenessScopeV1::Partial { requested: 1 },
    ))
    .unwrap_err();
    assert!(has_admission_issue(&report, |issue| matches!(
        issue,
        SpectralAdmissionIssueV1::ProjectivePrefixOrderingRequired
    )));

    let literal_real = structure_claim(
        projective_axes,
        StructurePropertyV1::RealSpectrum,
        StructureSupportV1::FormFree,
        WitnessDispositionV1::Witnessed,
        0.0,
        79,
    );
    let report = validate_problem(default_spec(
        projective_axes,
        vec![literal_real],
        Vec::new(),
        SpectralOrderingV1::RealAscending,
        CompletenessScopeV1::Partial { requested: 1 },
    ))
    .unwrap_err();
    assert!(
        report
            .issues()
            .contains(&SpectralAdmissionIssueV1::OrderingUnavailable)
    );
    assert!(has_admission_issue(&report, |issue| matches!(
        issue,
        SpectralAdmissionIssueV1::ProjectivePrefixOrderingRequired
    )));

    for (offset, policy) in [
        InfiniteEigenvaluePolicyV1::NoClaim,
        InfiniteEigenvaluePolicyV1::ExcludeWithCount,
    ]
    .into_iter()
    .enumerate()
    {
        let unaccounted_class = SpectralProblemClassV1::new(
            SpectralRepresentationV1::GeneralizedPencil,
            DescriptorRoleV1::Descriptor {
                infinity_policy: policy,
            },
            SpectralOperatorOriginV1::Direct,
        );
        let unaccounted_axes = axes(
            90_u8.wrapping_add(offset as u8),
            unaccounted_class,
            SpectralScalarFieldV1::Complex,
            Dims::NONE,
            3,
        );
        let report = validate_problem(default_spec(
            unaccounted_axes,
            Vec::new(),
            Vec::new(),
            SpectralOrderingV1::MagnitudeAscending {
                tie_break: ComplexTieBreakV1::RealThenImagThenLineage,
            },
            CompletenessScopeV1::Partial { requested: 1 },
        ))
        .unwrap_err();
        assert!(has_admission_issue(&report, |issue| matches!(
            issue,
            SpectralAdmissionIssueV1::InfinityPolicyMismatch
        )));
    }

    let unestablished = validate_problem(default_spec(
        projective_axes,
        Vec::new(),
        Vec::new(),
        SpectralOrderingV1::Projective {
            chart: SpectralProjectiveChartId::from_bytes([79; 32]),
            infinity: ProjectiveInfinityPlacementV1::Last,
            tie_break: ComplexTieBreakV1::ImagThenRealThenLineage,
        },
        CompletenessScopeV1::Partial { requested: 1 },
    ))
    .unwrap();
    let report = SpectralTruthV1::new(&unestablished, satisfied_partial_draft(&unestablished, 79))
        .unwrap_err();
    assert!(
        report
            .issues()
            .contains(&SpectralTruthErrorV1::DiscreteSpectrumRegularityNotEstablished)
    );

    let hdp = structure_claim(
        projective_axes,
        StructurePropertyV1::HermitianDefinitePencil,
        StructureSupportV1::InnerProduct(projective_axes.domain.id()),
        WitnessDispositionV1::Witnessed,
        0.0,
        80,
    );
    let regular_descriptor = regularity_claim(
        projective_axes,
        RegularityClassV1::RegularDescriptor,
        WitnessDispositionV1::Witnessed,
        81,
    );
    let report = validate_problem(default_spec(
        projective_axes,
        vec![hdp],
        vec![regular_descriptor],
        SpectralOrderingV1::Projective {
            chart: SpectralProjectiveChartId::from_bytes([80; 32]),
            infinity: ProjectiveInfinityPlacementV1::Last,
            tie_break: ComplexTieBreakV1::RealThenImagThenLineage,
        },
        CompletenessScopeV1::Partial { requested: 4 },
    ))
    .unwrap_err();
    assert!(has_admission_issue(&report, |issue| matches!(
        issue,
        SpectralAdmissionIssueV1::DimensionMismatch { left: 4, right: 3 }
    )));

    let established = validate_problem(default_spec(
        projective_axes,
        vec![hdp],
        vec![regular_descriptor],
        SpectralOrderingV1::Projective {
            chart: SpectralProjectiveChartId::from_bytes([81; 32]),
            infinity: ProjectiveInfinityPlacementV1::Last,
            tie_break: ComplexTieBreakV1::RealThenImagThenLineage,
        },
        CompletenessScopeV1::Partial { requested: 1 },
    ))
    .unwrap();
    assert!(SpectralTruthV1::new(&established, satisfied_partial_draft(&established, 81)).is_ok());

    let finite_real_prefix = validate_problem(default_spec(
        projective_axes,
        vec![hdp],
        vec![regular_descriptor],
        SpectralOrderingV1::RealAscending,
        CompletenessScopeV1::Partial { requested: 1 },
    ))
    .unwrap();
    assert!(finite_real_prefix.requires_real_spectrum_truth());
    assert!(finite_real_prefix.projective_infinity_is_excluded());
    assert!(
        SpectralTruthV1::new(
            &finite_real_prefix,
            satisfied_partial_draft(&finite_real_prefix, 82),
        )
        .is_ok()
    );

    let polynomial_class = SpectralProblemClassV1::new(
        SpectralRepresentationV1::MatrixPolynomial { grade: 2 },
        DescriptorRoleV1::Descriptor {
            infinity_policy: InfiniteEigenvaluePolicyV1::IncludeProjective,
        },
        SpectralOperatorOriginV1::Direct,
    );
    let polynomial = axes(
        82,
        polynomial_class,
        SpectralScalarFieldV1::Complex,
        Dims::NONE,
        3,
    );
    let leading = regularity_claim(
        polynomial,
        RegularityClassV1::InvertiblePolynomialLeadingCoefficient { grade: 2 },
        WitnessDispositionV1::Witnessed,
        82,
    );
    let descriptor = regularity_claim(
        polynomial,
        RegularityClassV1::RegularDescriptor,
        WitnessDispositionV1::Witnessed,
        83,
    );
    let report = validate_problem(default_spec(
        polynomial,
        Vec::new(),
        vec![leading, descriptor],
        SpectralOrderingV1::Projective {
            chart: SpectralProjectiveChartId::from_bytes([82; 32]),
            infinity: ProjectiveInfinityPlacementV1::Last,
            tie_break: ComplexTieBreakV1::RealThenImagThenLineage,
        },
        CompletenessScopeV1::Partial { requested: 7 },
    ))
    .unwrap_err();
    assert!(has_admission_issue(&report, |issue| matches!(
        issue,
        SpectralAdmissionIssueV1::DimensionMismatch { left: 7, right: 6 }
    )));

    let finite_magnitude_prefix = validate_problem(default_spec(
        polynomial,
        Vec::new(),
        vec![leading, descriptor],
        SpectralOrderingV1::MagnitudeAscending {
            tie_break: ComplexTieBreakV1::RealThenImagThenLineage,
        },
        CompletenessScopeV1::Partial { requested: 1 },
    ))
    .unwrap();
    assert!(finite_magnitude_prefix.projective_infinity_is_excluded());
    assert!(
        SpectralTruthV1::new(
            &finite_magnitude_prefix,
            satisfied_partial_draft(&finite_magnitude_prefix, 83),
        )
        .is_ok()
    );

    let contradicted_polynomial = validate_problem(default_spec(
        polynomial,
        Vec::new(),
        vec![
            regularity_claim(
                polynomial,
                RegularityClassV1::RegularPolynomial { grade: 2 },
                WitnessDispositionV1::Contradicted,
                84,
            ),
            descriptor,
        ],
        SpectralOrderingV1::Projective {
            chart: SpectralProjectiveChartId::from_bytes([84; 32]),
            infinity: ProjectiveInfinityPlacementV1::Last,
            tie_break: ComplexTieBreakV1::RealThenImagThenLineage,
        },
        CompletenessScopeV1::Partial { requested: 1 },
    ))
    .unwrap();
    assert_eq!(contradicted_polynomial.known_algebraic_cardinality(), None);
    let report = SpectralTruthV1::new(
        &contradicted_polynomial,
        satisfied_partial_draft(&contradicted_polynomial, 84),
    )
    .unwrap_err();
    assert!(
        report
            .issues()
            .contains(&SpectralTruthErrorV1::DiscreteSpectrumRegularityNotEstablished)
    );
}

#[test]
fn internal_separation_receipts_bind_membership_and_both_multiplicity_axes() {
    assert_eq!(SpectralResultSetIdentitySchemaV2::VERSION, 2);
    assert_eq!(
        SpectralResultSetIdentitySchemaV2::DOMAIN,
        "org.frankensim.fs-spectral.result-set.v2"
    );
    let problem = validate_problem(default_spec(
        standard_axes(80, 2),
        Vec::new(),
        Vec::new(),
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap();
    let id = cluster_id(80);
    let enclosure = SpectralEnclosureV1::Real(FiniteIntervalV1::new(1.0, 1.0).unwrap());
    let exact_multiplicity = MultiplicityClaimV1::Exact {
        value: 2,
        witness: truth_witness(
            problem.problem_id(),
            SpectralTruthPropositionV1::Multiplicity {
                cluster: id,
                enclosure,
                kind: MultiplicityKindV1::Algebraic,
                assertion: MultiplicityAssertionV1::Exact,
                lower: 2,
                upper: Some(2),
            },
            80,
        ),
    };
    let degeneracy_witness = truth_witness(
        problem.problem_id(),
        SpectralTruthPropositionV1::InternalDegeneracy {
            cluster: id,
            enclosure,
            algebraic: MultiplicityStatementV1::Exact { value: 2 },
            geometric: MultiplicityStatementV1::Unknown,
        },
        81,
    );
    let cluster = SpectralClusterV1::new(
        id,
        SpectralLocalizationV1::candidate(enclosure),
        exact_multiplicity,
        MultiplicityClaimV1::Unknown,
        InternalClusterStateV1::ProvenDegenerate {
            witness: degeneracy_witness,
        },
    )
    .unwrap();
    assert!(
        SpectralTruthV1::new(
            &problem,
            SpectralTruthDraftV1::new(
                SpectralResultAuthorityV1::NoClaim,
                SpectralCoverageV1::Candidates,
                vec![cluster],
                ScopeBoundaryStateV1::NoClaim,
                SpectralTerminationV1::Completed,
            ),
        )
        .is_ok()
    );

    let bounded_multiplicity = MultiplicityClaimV1::Bounds {
        lower: 2,
        upper: 2,
        witness: truth_witness(
            problem.problem_id(),
            SpectralTruthPropositionV1::Multiplicity {
                cluster: id,
                enclosure,
                kind: MultiplicityKindV1::Algebraic,
                assertion: MultiplicityAssertionV1::Bounds,
                lower: 2,
                upper: Some(2),
            },
            82,
        ),
    };
    let rebound = SpectralClusterV1::new(
        id,
        SpectralLocalizationV1::candidate(enclosure),
        bounded_multiplicity,
        MultiplicityClaimV1::Unknown,
        InternalClusterStateV1::ProvenDegenerate {
            witness: degeneracy_witness,
        },
    )
    .unwrap();
    let report = SpectralTruthV1::new(
        &problem,
        SpectralTruthDraftV1::new(
            SpectralResultAuthorityV1::NoClaim,
            SpectralCoverageV1::Candidates,
            vec![rebound],
            ScopeBoundaryStateV1::NoClaim,
            SpectralTerminationV1::Completed,
        ),
    )
    .unwrap_err();
    assert!(has_truth_issue(&report, |issue| matches!(
        issue,
        SpectralTruthErrorV1::WitnessPropositionMismatch { .. }
    )));

    let rebound_id = cluster_id(82);
    let rebound_membership = SpectralClusterV1::new(
        rebound_id,
        SpectralLocalizationV1::candidate(enclosure),
        MultiplicityClaimV1::Exact {
            value: 2,
            witness: truth_witness(
                problem.problem_id(),
                SpectralTruthPropositionV1::Multiplicity {
                    cluster: rebound_id,
                    enclosure,
                    kind: MultiplicityKindV1::Algebraic,
                    assertion: MultiplicityAssertionV1::Exact,
                    lower: 2,
                    upper: Some(2),
                },
                83,
            ),
        },
        MultiplicityClaimV1::Unknown,
        InternalClusterStateV1::ProvenDegenerate {
            witness: degeneracy_witness,
        },
    )
    .unwrap();
    let report = SpectralTruthV1::new(
        &problem,
        SpectralTruthDraftV1::new(
            SpectralResultAuthorityV1::NoClaim,
            SpectralCoverageV1::Candidates,
            vec![rebound_membership],
            ScopeBoundaryStateV1::NoClaim,
            SpectralTerminationV1::Completed,
        ),
    )
    .unwrap_err();
    assert!(has_truth_issue(&report, |issue| matches!(
        issue,
        SpectralTruthErrorV1::WitnessPropositionMismatch { .. }
    )));

    let resolved_id = cluster_id(83);
    let resolved_algebraic = MultiplicityClaimV1::Exact {
        value: 2,
        witness: truth_witness(
            problem.problem_id(),
            SpectralTruthPropositionV1::Multiplicity {
                cluster: resolved_id,
                enclosure,
                kind: MultiplicityKindV1::Algebraic,
                assertion: MultiplicityAssertionV1::Exact,
                lower: 2,
                upper: Some(2),
            },
            84,
        ),
    };
    let lower = PositiveFiniteV1::new(0.125).unwrap();
    let norm = norm_id(83);
    let resolved = SpectralClusterV1::new(
        resolved_id,
        SpectralLocalizationV1::candidate(enclosure),
        resolved_algebraic,
        MultiplicityClaimV1::Unknown,
        InternalClusterStateV1::Resolved {
            lower,
            norm,
            witness: truth_witness(
                problem.problem_id(),
                SpectralTruthPropositionV1::InternalResolution {
                    cluster: resolved_id,
                    enclosure,
                    algebraic: MultiplicityStatementV1::Exact { value: 2 },
                    geometric: MultiplicityStatementV1::Unknown,
                    lower,
                    norm,
                },
                85,
            ),
        },
    )
    .unwrap();
    let truth = SpectralTruthV1::new(
        &problem,
        SpectralTruthDraftV1::new(
            SpectralResultAuthorityV1::NoClaim,
            SpectralCoverageV1::Candidates,
            vec![resolved],
            ScopeBoundaryStateV1::NoClaim,
            SpectralTerminationV1::Completed,
        ),
    )
    .unwrap();
    assert!(matches!(
        truth.clusters()[0].internal(),
        InternalClusterStateV1::Resolved {
            lower: found_lower,
            norm: found_norm,
            ..
        } if found_lower == lower && found_norm == norm
    ));

    assert_eq!(
        SpectralClusterV1::new(
            cluster_id(81),
            SpectralLocalizationV1::candidate(enclosure),
            resolved_algebraic,
            MultiplicityClaimV1::Unknown,
            InternalClusterStateV1::UndefinedSeparation {
                reason: UndefinedSeparationReasonV1::ProjectiveInfinityInAffineCoordinates,
            },
        ),
        Err(SpectralTruthErrorV1::InvalidInternalClusterState),
        "a finite enclosure has a mathematically defined affine separation proposition"
    );

    let projective_problem =
        descriptor_full_problem(86, InfiniteEigenvaluePolicyV1::IncludeProjective);
    let projective_id = cluster_id(86);
    let projective_enclosure = SpectralEnclosureV1::ProjectiveInfinity;
    assert_eq!(
        SpectralClusterV1::new(
            projective_id,
            SpectralLocalizationV1::candidate(projective_enclosure),
            MultiplicityClaimV1::Unknown,
            MultiplicityClaimV1::Unknown,
            InternalClusterStateV1::UndefinedSeparation {
                reason: UndefinedSeparationReasonV1::ProjectiveInfinityInAffineCoordinates,
            },
        ),
        Err(SpectralTruthErrorV1::InvalidInternalClusterState),
        "undefined internal separation requires evidence that the projective cluster is repeated"
    );
    let singleton_id = cluster_id(85);
    let singleton_multiplicity = MultiplicityClaimV1::Exact {
        value: 1,
        witness: truth_witness(
            projective_problem.problem_id(),
            SpectralTruthPropositionV1::Multiplicity {
                cluster: singleton_id,
                enclosure: projective_enclosure,
                kind: MultiplicityKindV1::Algebraic,
                assertion: MultiplicityAssertionV1::Exact,
                lower: 1,
                upper: Some(1),
            },
            85,
        ),
    };
    assert_eq!(
        SpectralClusterV1::new(
            singleton_id,
            SpectralLocalizationV1::candidate(projective_enclosure),
            singleton_multiplicity,
            MultiplicityClaimV1::Unknown,
            InternalClusterStateV1::UndefinedSeparation {
                reason: UndefinedSeparationReasonV1::ProjectiveInfinityInAffineCoordinates,
            },
        ),
        Err(SpectralTruthErrorV1::InvalidInternalClusterState),
        "a singleton has vacuous internal separation and must use Simple rather than UndefinedSeparation"
    );
    let replayed_multiplicity = MultiplicityClaimV1::Exact {
        value: 2,
        witness: truth_witness(
            problem.problem_id(),
            SpectralTruthPropositionV1::Multiplicity {
                cluster: resolved_id,
                enclosure,
                kind: MultiplicityKindV1::Algebraic,
                assertion: MultiplicityAssertionV1::Exact,
                lower: 2,
                upper: Some(2),
            },
            86,
        ),
    };
    let replayed_undefined = SpectralClusterV1::new(
        projective_id,
        SpectralLocalizationV1::candidate(projective_enclosure),
        replayed_multiplicity,
        MultiplicityClaimV1::Unknown,
        InternalClusterStateV1::UndefinedSeparation {
            reason: UndefinedSeparationReasonV1::ProjectiveInfinityInAffineCoordinates,
        },
    )
    .unwrap();
    let report = SpectralTruthV1::new(
        &projective_problem,
        SpectralTruthDraftV1::new(
            SpectralResultAuthorityV1::NoClaim,
            SpectralCoverageV1::Candidates,
            vec![replayed_undefined],
            ScopeBoundaryStateV1::NoClaim,
            SpectralTerminationV1::Completed,
        ),
    )
    .unwrap_err();
    assert!(
        has_truth_issue(&report, |issue| matches!(
            issue,
            SpectralTruthErrorV1::WitnessPropositionMismatch { .. }
        )),
        "undefined separation must not launder replayed cluster/multiplicity evidence: {report:?}"
    );

    let projective_multiplicity = MultiplicityClaimV1::Exact {
        value: 2,
        witness: truth_witness(
            projective_problem.problem_id(),
            SpectralTruthPropositionV1::Multiplicity {
                cluster: projective_id,
                enclosure: projective_enclosure,
                kind: MultiplicityKindV1::Algebraic,
                assertion: MultiplicityAssertionV1::Exact,
                lower: 2,
                upper: Some(2),
            },
            87,
        ),
    };
    let no_claim = SpectralClusterV1::new(
        projective_id,
        SpectralLocalizationV1::candidate(projective_enclosure),
        projective_multiplicity,
        MultiplicityClaimV1::Unknown,
        InternalClusterStateV1::NoClaim,
    )
    .unwrap();
    let unknown = SpectralClusterV1::new(
        projective_id,
        SpectralLocalizationV1::candidate(projective_enclosure),
        projective_multiplicity,
        MultiplicityClaimV1::Unknown,
        InternalClusterStateV1::Unknown {
            reason: UnknownSeparationReasonV1::MissingEvidence,
        },
    )
    .unwrap();
    let undefined = SpectralClusterV1::new(
        projective_id,
        SpectralLocalizationV1::candidate(projective_enclosure),
        projective_multiplicity,
        MultiplicityClaimV1::Unknown,
        InternalClusterStateV1::UndefinedSeparation {
            reason: UndefinedSeparationReasonV1::ProjectiveInfinityInAffineCoordinates,
        },
    )
    .unwrap();
    let no_claim_receipt = spectral_result_set_receipt(&[no_claim]).unwrap();
    let unknown_receipt = spectral_result_set_receipt(&[unknown]).unwrap();
    let undefined_receipt = spectral_result_set_receipt(&[undefined]).unwrap();
    assert_ne!(no_claim_receipt, unknown_receipt);
    assert_ne!(no_claim_receipt, undefined_receipt);
    assert_ne!(unknown_receipt, undefined_receipt);

    let truth = SpectralTruthV1::new(
        &projective_problem,
        SpectralTruthDraftV1::new(
            SpectralResultAuthorityV1::NoClaim,
            SpectralCoverageV1::Candidates,
            vec![undefined],
            ScopeBoundaryStateV1::NoClaim,
            SpectralTerminationV1::Completed,
        ),
    )
    .unwrap();
    assert_eq!(truth.result_set_identity_receipt(), undefined_receipt);
    assert_eq!(truth.result_set_id(), undefined_receipt.id());
    let retained = ObservedIdentity::from_receipt(truth.result_set_identity_receipt());
    let synthetic_collision = ObservedIdentity::presented(
        truth.result_set_id(),
        ByteObservation::new(
            undefined_receipt.canonical_preimage(),
            undefined_receipt.canonical_bytes() + 1,
        ),
    );
    assert!(
        matches!(
            adjudicate(retained, synthetic_collision),
            IdentityAdjudication::Refused(_)
        ),
        "the validated truth must retain enough result-set observation data for collision adjudication"
    );
    assert!(matches!(
        truth.clusters()[0].internal(),
        InternalClusterStateV1::UndefinedSeparation {
            reason: UndefinedSeparationReasonV1::ProjectiveInfinityInAffineCoordinates,
        }
    ));
}

/// The rogue-pair hole the charter pin closes (bead sj31i.52.9): a
/// permit-everything admission PAIRED with a witness from a SELF-CONFIGURED
/// root carries the same rogue identities on both sides, so every
/// subject/anchor/verifier/policy/context check passes — only the pinned
/// root charter exposes the foreign configuration.
#[test]
fn self_configured_root_witnesses_fail_the_pinned_charter() {
    let axes = standard_axes(43, 4);
    let receipt = regularity_proposition_receipt(
        axes.subject,
        axes.scalar,
        axes.class,
        axes.scaling,
        axes.domain,
        axes.codomain,
        RegularityClassV1::FiniteDimensional,
        WitnessDispositionV1::Witnessed,
    )
    .unwrap();
    let trusted_verifier =
        spectral_verifier_receipt(b"fs-spectral-test-exact-verifier-v1").unwrap();
    let trusted_policy =
        spectral_authority_policy_receipt(b"fs-spectral-test-admission-policy-v1").unwrap();
    let pinned = spectral_promotion_trust_root(trusted_verifier, trusted_policy)
        .unwrap()
        .charter();

    // The adversary self-configures a root around its OWN rogue identities
    // and mints a matching admission + witness pair.
    let rogue_verifier = spectral_verifier_receipt(b"foreign-permit-all-verifier").unwrap();
    let rogue_policy = spectral_authority_policy_receipt(b"foreign-permit-all-policy").unwrap();
    let anchor = ExternalAnchorRef::presented(ContentId::of_bytes(b"self-configured-pair"));
    let rogue_admitted =
        AuthorityRef::present(receipt, anchor, rogue_verifier.id(), rogue_policy.id())
            .verify(&PermitAll)
            .unwrap()
            .admit(&PermitAll)
            .unwrap();
    let rogue_root = spectral_promotion_trust_root(rogue_verifier, rogue_policy).unwrap();
    let rogue_witness = rogue_root
        .admit_for_promotion(
            &rogue_admitted,
            ObservedIdentity::from_receipt(rogue_verifier).bytes(),
            ObservedIdentity::from_receipt(rogue_policy).bytes(),
        )
        .expect("a self-configured root promotes its own binding by design");

    // Every identity axis matches between the pair; only the charter differs.
    assert_eq!(
        AdmittedSpectralWitnessV1::from_authority(&rogue_admitted, rogue_witness, pinned),
        Err(SpectralPromotionBindingErrorV1::RootCharter)
    );
    // The pair passes when the consumer pins the ROGUE charter — proving the
    // refusal above is exactly the provenance discrimination, not an
    // incidental binding mismatch.
    AdmittedSpectralWitnessV1::from_authority(&rogue_admitted, rogue_witness, rogue_root.charter())
        .expect("identity axes all match within the rogue pair");
}
