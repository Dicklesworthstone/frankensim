//! G0/G3/G5 conformance tests for RB.8a physical-domain adapter admission.

#![allow(clippy::wildcard_imports)]

use fs_blake3::identity::{
    AuthorityAdmitter, AuthorityRef, AuthorityVerifier, ContentId, ExternalAnchorRef,
    IdentityReceipt,
};
use fs_qty::{Dims, Time};
use fs_spectral::adapter::*;
use fs_spectral::admission::*;

fn bytes(seed: u8) -> [u8; 32] {
    [seed; 32]
}

fn map(seed: u8) -> SpectralAdapterMapIdV1 {
    SpectralAdapterMapIdV1::from_bytes(bytes(seed))
}

fn no_claim(seed: u8) -> SpectralAdapterNoClaimIdV1 {
    SpectralAdapterNoClaimIdV1::from_bytes(bytes(seed))
}

fn scaling(seed: u8, dims: Dims) -> SpectralScalingContextV1 {
    SpectralScalingContextV1::new(
        SpectralScalingId::from_bytes(bytes(seed)),
        dims,
        1.0,
        SpectralScalingMapId::from_bytes(bytes(seed.wrapping_add(1))),
        SpectralScalingMapId::from_bytes(bytes(seed.wrapping_add(2))),
        SpectralScalingMapId::from_bytes(bytes(seed.wrapping_add(3))),
        SpectralScalingMapId::from_bytes(bytes(seed.wrapping_add(4))),
    )
}

fn problem(
    seed: u8,
    class: SpectralProblemClassV1,
    scalar: SpectralScalarFieldV1,
    dims: Dims,
) -> ValidatedSpectralProblemV1 {
    let metric = SpectralMetricV1::euclidean(4);
    validate_problem(SpectralProblemSpecV1::new(
        SpectralSubjectId::from_bytes(bytes(seed)),
        scalar,
        class,
        StructureProfileV1::new(Vec::new()),
        scaling(seed.wrapping_add(20), dims),
        SpectralSpaceContextV1::new(
            metric,
            metric,
            GaugeConventionV1::Unknown,
            ZeroPaddingConventionV1::Unknown,
        ),
        RegularityProfileV1::new(Vec::new()),
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .expect("test problem must admit")
}

fn standard_problem(seed: u8) -> ValidatedSpectralProblemV1 {
    problem(
        seed,
        SpectralProblemClassV1::new(
            SpectralRepresentationV1::StandardLinear,
            DescriptorRoleV1::Ordinary,
            SpectralOperatorOriginV1::Direct,
        ),
        SpectralScalarFieldV1::Real,
        Dims::NONE,
    )
}

fn descriptor_problem(seed: u8) -> ValidatedSpectralProblemV1 {
    problem(
        seed,
        SpectralProblemClassV1::new(
            SpectralRepresentationV1::GeneralizedPencil,
            DescriptorRoleV1::Descriptor {
                infinity_policy: InfiniteEigenvaluePolicyV1::IncludeProjective,
            },
            SpectralOperatorOriginV1::Direct,
        ),
        SpectralScalarFieldV1::Complex,
        Dims::NONE,
    )
}

fn periodic_problem(seed: u8) -> ValidatedSpectralProblemV1 {
    problem(
        seed,
        SpectralProblemClassV1::new(
            SpectralRepresentationV1::StandardLinear,
            DescriptorRoleV1::Ordinary,
            SpectralOperatorOriginV1::MonodromyFloquet {
                period: Time::new(0.25),
                parameter: FloquetParameterV1::Multiplier,
                branch: FloquetBranchConventionV1::MultipliersOnly,
            },
        ),
        SpectralScalarFieldV1::Complex,
        Dims::NONE,
    )
}

#[derive(Clone, Copy)]
struct SourceOptions {
    domain: PhysicalSourceDomainV1,
    model_version: PhysicalModelVersionIdV1,
    operator_class: PhysicalOperatorClassV1,
    state_dimension: u32,
    dual_dimension: u32,
    target_scaling: SpectralScalingId,
    target_domain_metric: SpectralMetricId,
    target_codomain_metric: SpectralMetricId,
    frame: FrameCrosswalkV1,
    constraints: AdapterBindingV1<PhysicalConstraintSetIdV1>,
    nullspace: AdapterBindingV1<PhysicalNullspaceIdV1>,
    parameters: AdapterBindingV1<PhysicalParameterSchemaIdV1>,
    boundaries: AdapterBindingV1<PhysicalBoundarySchemaIdV1>,
    linearization: LinearizationContextV1,
    phase_section: PhaseSectionContextV1,
    structure: AdapterBindingV1<PhysicalStructureWitnessIdV1>,
}

impl SourceOptions {
    fn standard(target: &ValidatedSpectralProblemV1) -> Self {
        let model_version = PhysicalModelVersionIdV1::from_bytes(bytes(2));
        Self {
            domain: PhysicalSourceDomainV1::Control,
            model_version,
            operator_class: PhysicalOperatorClassV1::StandardLinear,
            state_dimension: target.spec().spaces().domain().dimension(),
            dual_dimension: target.spec().spaces().codomain().dimension(),
            target_scaling: target.spec().scaling().id(),
            target_domain_metric: target.spec().spaces().domain().id(),
            target_codomain_metric: target.spec().spaces().codomain().id(),
            frame: FrameCrosswalkV1 {
                source_frame: PhysicalFrameIdV1::from_bytes(bytes(7)),
                target_frame: PhysicalFrameIdV1::from_bytes(bytes(8)),
                map: map(9),
                kind: FrameMapKindV1::ExactTransform,
            },
            constraints: AdapterBindingV1::NotApplicable {
                justification: no_claim(10),
            },
            nullspace: AdapterBindingV1::NotApplicable {
                justification: no_claim(11),
            },
            parameters: AdapterBindingV1::Retained {
                artifact: PhysicalParameterSchemaIdV1::from_bytes(bytes(12)),
                map: map(13),
            },
            boundaries: AdapterBindingV1::Retained {
                artifact: PhysicalBoundarySchemaIdV1::from_bytes(bytes(14)),
                map: map(15),
            },
            linearization: LinearizationContextV1::Frozen {
                point: PhysicalLinearizationPointIdV1::from_bytes(bytes(16)),
                model_version,
            },
            phase_section: PhaseSectionContextV1::NotPeriodic {
                justification: no_claim(17),
            },
            structure: AdapterBindingV1::NotApplicable {
                justification: no_claim(18),
            },
        }
    }
}

fn source(options: SourceOptions) -> PhysicalSourceContextV1 {
    PhysicalSourceContextV1::new(
        options.domain,
        PhysicalSourceArtifactIdV1::from_bytes(bytes(1)),
        options.model_version,
        options.operator_class,
        PhysicalStateSpaceIdV1::from_bytes(bytes(3)),
        PhysicalDualSpaceIdV1::from_bytes(bytes(4)),
        options.state_dimension,
        options.dual_dimension,
        UnitCrosswalkV1 {
            source_units: PhysicalUnitSystemIdV1::from_bytes(bytes(5)),
            target_scaling: options.target_scaling,
            map: map(6),
        },
        options.frame,
        MetricNormCrosswalkV1 {
            source_metric: PhysicalMetricIdV1::from_bytes(bytes(19)),
            source_norm: PhysicalNormIdV1::from_bytes(bytes(20)),
            target_domain_metric: options.target_domain_metric,
            target_codomain_metric: options.target_codomain_metric,
            metric_map: map(21),
            norm_map: map(22),
            target_norm: SpectralNormId::from_bytes(bytes(23)),
        },
        options.constraints,
        options.nullspace,
        options.parameters,
        options.boundaries,
        options.linearization,
        options.phase_section,
        options.structure,
    )
}

fn qoi(seed: u8) -> SpectralQoiCrosswalkV1 {
    SpectralQoiCrosswalkV1::new(
        PhysicalQoiIdV1::from_bytes(bytes(seed)),
        SpectralQoiIdV1::from_bytes(bytes(seed.wrapping_add(1))),
        map(seed.wrapping_add(2)),
        ReverseInterpretationV1::Partial {
            map: map(seed.wrapping_add(3)),
            no_claim: no_claim(seed.wrapping_add(4)),
        },
    )
}

fn adapter_spec(
    target: &ValidatedSpectralProblemV1,
    options: SourceOptions,
    qois: Vec<SpectralQoiCrosswalkV1>,
) -> SpectralAdapterSpecV1 {
    SpectralAdapterSpecV1::new(
        source(options),
        target.problem_id(),
        AdapterFidelityV1::ExactOneWay {
            forward: map(24),
            no_inverse: no_claim(25),
        },
        ReverseInterpretationV1::Unavailable {
            no_claim: no_claim(26),
        },
        qois,
    )
}

fn assert_issue(
    result: Result<ValidatedSpectralAdapterV1, SpectralAdapterReportV1>,
    expected: SpectralAdapterIssueV1,
) {
    let report = result.expect_err("adapter must refuse adversarial input");
    assert!(
        report.issues().contains(&expected),
        "missing {expected:?} in {:?}",
        report.issues()
    );
}

#[test]
fn canonical_qoi_order_and_replay_are_identity_stable() {
    let target = standard_problem(40);
    let options = SourceOptions::standard(&target);
    let first = validate_adapter_v1(
        adapter_spec(&target, options, vec![qoi(40), qoi(30)]),
        &target,
    )
    .unwrap();
    let reordered = validate_adapter_v1(
        adapter_spec(&target, options, vec![qoi(30), qoi(40)]),
        &target,
    )
    .unwrap();

    assert_eq!(first.adapter_id(), reordered.adapter_id());
    assert_eq!(
        first.identity_receipt().canonical_preimage(),
        reordered.identity_receipt().canonical_preimage()
    );
    assert_eq!(first.qois()[0].source(), qoi(30).source());

    let replay = validate_adapter_v1(first.spec().clone(), &target).unwrap();
    assert_eq!(replay.identity_receipt(), first.identity_receipt());

    let mut changed_version = options;
    changed_version.model_version = PhysicalModelVersionIdV1::from_bytes(bytes(88));
    changed_version.linearization = LinearizationContextV1::Frozen {
        point: PhysicalLinearizationPointIdV1::from_bytes(bytes(16)),
        model_version: changed_version.model_version,
    };
    let changed_version = validate_adapter_v1(
        adapter_spec(&target, changed_version, vec![qoi(30), qoi(40)]),
        &target,
    )
    .unwrap();
    assert_ne!(first.adapter_id(), changed_version.adapter_id());

    let mut changed_frame = options;
    changed_frame.frame.source_frame = PhysicalFrameIdV1::from_bytes(bytes(89));
    let changed_frame = validate_adapter_v1(
        adapter_spec(&target, changed_frame, vec![qoi(30), qoi(40)]),
        &target,
    )
    .unwrap();
    assert_ne!(first.adapter_id(), changed_frame.adapter_id());
}

#[test]
fn missing_or_wrong_units_metrics_frames_and_context_refuse() {
    let target = standard_problem(41);
    let options = SourceOptions::standard(&target);

    let mut wrong_metric = options;
    wrong_metric.target_domain_metric = SpectralMetricId::from_bytes(bytes(90));
    assert_issue(
        validate_adapter_v1(adapter_spec(&target, wrong_metric, Vec::new()), &target),
        SpectralAdapterIssueV1::MetricMismatch,
    );

    let mut wrong_units = options;
    wrong_units.target_scaling = SpectralScalingId::from_bytes(bytes(91));
    assert_issue(
        validate_adapter_v1(adapter_spec(&target, wrong_units, Vec::new()), &target),
        SpectralAdapterIssueV1::ScalingMismatch,
    );

    let mut wrong_frame = options;
    wrong_frame.frame.kind = FrameMapKindV1::Identity;
    assert_issue(
        validate_adapter_v1(adapter_spec(&target, wrong_frame, Vec::new()), &target),
        SpectralAdapterIssueV1::FrameMismatch,
    );

    let mut lossy_frame = options;
    lossy_frame.frame.kind = FrameMapKindV1::Lossy;
    assert_issue(
        validate_adapter_v1(adapter_spec(&target, lossy_frame, Vec::new()), &target),
        SpectralAdapterIssueV1::InadmissibleFrameMap,
    );

    let mut stale = options;
    stale.linearization = LinearizationContextV1::Frozen {
        point: PhysicalLinearizationPointIdV1::from_bytes(bytes(16)),
        model_version: PhysicalModelVersionIdV1::from_bytes(bytes(92)),
    };
    assert_issue(
        validate_adapter_v1(adapter_spec(&target, stale, Vec::new()), &target),
        SpectralAdapterIssueV1::StaleLinearization,
    );

    let mut missing_boundary = options;
    missing_boundary.boundaries = AdapterBindingV1::Unknown;
    assert_issue(
        validate_adapter_v1(adapter_spec(&target, missing_boundary, Vec::new()), &target),
        SpectralAdapterIssueV1::UnknownBinding {
            field: AdapterFieldV1::Boundaries,
        },
    );
}

#[test]
fn lossy_ambiguous_duplicate_and_oversized_crosswalks_refuse() {
    let target = standard_problem(42);
    let options = SourceOptions::standard(&target);
    let base_source = source(options);

    let lossy = SpectralAdapterSpecV1::new(
        base_source.clone(),
        target.problem_id(),
        AdapterFidelityV1::Lossy {
            reason: no_claim(93),
        },
        ReverseInterpretationV1::Unavailable {
            no_claim: no_claim(26),
        },
        Vec::new(),
    );
    assert_issue(
        validate_adapter_v1(lossy, &target),
        SpectralAdapterIssueV1::InadmissibleFidelity,
    );

    let reverse_overclaim = SpectralAdapterSpecV1::new(
        base_source.clone(),
        target.problem_id(),
        AdapterFidelityV1::ExactOneWay {
            forward: map(24),
            no_inverse: no_claim(25),
        },
        ReverseInterpretationV1::Exact { map: map(101) },
        Vec::new(),
    );
    assert_issue(
        validate_adapter_v1(reverse_overclaim, &target),
        SpectralAdapterIssueV1::ReverseInterpretationOverclaim,
    );

    let quotient_kernel = PhysicalNullspaceIdV1::from_bytes(bytes(102));
    let quotient_mismatch = SpectralAdapterSpecV1::new(
        base_source.clone(),
        target.problem_id(),
        AdapterFidelityV1::ExactQuotient {
            forward: map(103),
            kernel: quotient_kernel,
            no_inverse: no_claim(104),
        },
        ReverseInterpretationV1::Unavailable {
            no_claim: no_claim(105),
        },
        Vec::new(),
    );
    assert_issue(
        validate_adapter_v1(quotient_mismatch, &target),
        SpectralAdapterIssueV1::QuotientKernelMismatch,
    );

    let mut quotient_options = options;
    quotient_options.nullspace = AdapterBindingV1::Retained {
        artifact: quotient_kernel,
        map: map(106),
    };
    let quotient = SpectralAdapterSpecV1::new(
        source(quotient_options),
        target.problem_id(),
        AdapterFidelityV1::ExactQuotient {
            forward: map(103),
            kernel: quotient_kernel,
            no_inverse: no_claim(104),
        },
        ReverseInterpretationV1::Unavailable {
            no_claim: no_claim(105),
        },
        Vec::new(),
    );
    validate_adapter_v1(quotient, &target).unwrap();

    let ambiguous = SpectralAdapterSpecV1::new(
        base_source,
        target.problem_id(),
        AdapterFidelityV1::Ambiguous,
        ReverseInterpretationV1::Unavailable {
            no_claim: no_claim(26),
        },
        Vec::new(),
    );
    assert_issue(
        validate_adapter_v1(ambiguous, &target),
        SpectralAdapterIssueV1::InadmissibleFidelity,
    );

    assert_issue(
        validate_adapter_v1(
            adapter_spec(&target, options, vec![qoi(50), qoi(50)]),
            &target,
        ),
        SpectralAdapterIssueV1::DuplicateQoi,
    );

    assert_issue(
        validate_adapter_v1(
            adapter_spec(
                &target,
                options,
                vec![qoi(60); MAX_SPECTRAL_ADAPTER_QOIS_V1 + 1],
            ),
            &target,
        ),
        SpectralAdapterIssueV1::TooManyQois {
            found: MAX_SPECTRAL_ADAPTER_QOIS_V1 + 1,
            limit: MAX_SPECTRAL_ADAPTER_QOIS_V1,
        },
    );
}

#[test]
fn descriptor_constraints_and_periodic_phase_are_cross_checked() {
    let descriptor = descriptor_problem(43);
    let mut descriptor_options = SourceOptions::standard(&descriptor);
    descriptor_options.operator_class =
        PhysicalOperatorClassV1::GeneralizedPencil { descriptor: true };
    assert_issue(
        validate_adapter_v1(
            adapter_spec(&descriptor, descriptor_options, Vec::new()),
            &descriptor,
        ),
        SpectralAdapterIssueV1::DescriptorConstraintsMissing,
    );
    descriptor_options.constraints = AdapterBindingV1::Retained {
        artifact: PhysicalConstraintSetIdV1::from_bytes(bytes(94)),
        map: map(95),
    };
    validate_adapter_v1(
        adapter_spec(&descriptor, descriptor_options, Vec::new()),
        &descriptor,
    )
    .unwrap();

    let periodic = periodic_problem(44);
    let mut periodic_options = SourceOptions::standard(&periodic);
    periodic_options.domain = PhysicalSourceDomainV1::PeriodicDynamics;
    assert_issue(
        validate_adapter_v1(
            adapter_spec(&periodic, periodic_options, Vec::new()),
            &periodic,
        ),
        SpectralAdapterIssueV1::PhaseOriginMismatch,
    );
    periodic_options.phase_section = PhaseSectionContextV1::Periodic {
        artifact: PhysicalPhaseSectionIdV1::from_bytes(bytes(96)),
        map: map(97),
    };
    validate_adapter_v1(
        adapter_spec(&periodic, periodic_options, Vec::new()),
        &periodic,
    )
    .unwrap();
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
    verifier: SpectralAuthorityVerifierIdV1,
    policy: SpectralAuthorityPolicyIdV1,
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
            && presented.verifier() == self.verifier
            && presented.key_policy() == self.policy
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
            && verified.verifier() == self.verifier
            && verified.key_policy() == self.policy
        {
            Ok(())
        } else {
            Err(AuthorityError::Mismatch)
        }
    }
}

fn admitted_witness(receipt: IdentityReceipt<SpectralPropositionId>) -> AdmittedSpectralWitnessV1 {
    let authority = ExactAuthority {
        proposition: receipt.id(),
        preimage: receipt.canonical_preimage(),
        anchor: ExternalAnchorRef::presented(ContentId::of_bytes(b"adapter-structure-test")),
        verifier: spectral_verifier_receipt(b"adapter-exact-verifier")
            .unwrap()
            .id(),
        policy: spectral_authority_policy_receipt(b"adapter-test-policy")
            .unwrap()
            .id(),
    };
    let presented = AuthorityRef::present(
        receipt,
        authority.anchor,
        authority.verifier,
        authority.policy,
    );
    let admitted = presented
        .verify(&authority)
        .unwrap()
        .admit(&authority)
        .unwrap();
    AdmittedSpectralWitnessV1::from_authority(admitted)
}

fn structured_problem(seed: u8) -> ValidatedSpectralProblemV1 {
    let metric = SpectralMetricV1::euclidean(4);
    let class = SpectralProblemClassV1::new(
        SpectralRepresentationV1::StandardLinear,
        DescriptorRoleV1::Ordinary,
        SpectralOperatorOriginV1::Direct,
    );
    let scaling = scaling(seed.wrapping_add(20), Dims::NONE);
    let subject = SpectralSubjectId::from_bytes(bytes(seed));
    let support = StructureSupportV1::InnerProduct(metric.id());
    let norm = SpectralNormId::from_bytes(bytes(98));
    let receipt = structure_proposition_receipt(
        subject,
        SpectralScalarFieldV1::Real,
        class,
        scaling,
        metric,
        metric,
        StructurePropertyV1::SelfAdjoint,
        support,
        WitnessDispositionV1::Witnessed,
        0.0,
        norm,
    )
    .unwrap();
    let claim = StructureClaimV1::new(
        StructurePropertyV1::SelfAdjoint,
        support,
        WitnessDispositionV1::Witnessed,
        0.0,
        norm,
        admitted_witness(receipt),
    );
    validate_problem(SpectralProblemSpecV1::new(
        subject,
        SpectralScalarFieldV1::Real,
        class,
        StructureProfileV1::new(vec![claim]),
        scaling,
        SpectralSpaceContextV1::new(
            metric,
            metric,
            GaugeConventionV1::Unknown,
            ZeroPaddingConventionV1::Unknown,
        ),
        RegularityProfileV1::new(Vec::new()),
        SpectralOrderingV1::SetValued,
        CompletenessScopeV1::CandidateOnly,
    ))
    .unwrap()
}

#[test]
fn target_structure_claims_require_a_retained_source_witness() {
    let target = structured_problem(45);
    let mut options = SourceOptions::standard(&target);
    assert_issue(
        validate_adapter_v1(adapter_spec(&target, options, Vec::new()), &target),
        SpectralAdapterIssueV1::StructureNotRetained,
    );
    options.structure = AdapterBindingV1::Retained {
        artifact: PhysicalStructureWitnessIdV1::from_bytes(bytes(99)),
        map: map(100),
    };
    validate_adapter_v1(adapter_spec(&target, options, Vec::new()), &target).unwrap();
}
