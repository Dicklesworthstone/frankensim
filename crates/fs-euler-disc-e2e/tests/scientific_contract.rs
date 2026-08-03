//! G0/G3 conformance and hostile tests for the Euler-disc scientific contract.

#![allow(clippy::too_many_lines)] // Each named battery is a hand-maintained semantic oracle.

use std::{
    collections::{BTreeMap, BTreeSet},
    mem::size_of,
};

use fs_blake3::ContentHash;
use fs_euler_disc_e2e::contract::{
    MAX_EULER_CONTRACT_BYTES, MAX_EULER_GRAPH_BYTES, MAX_EULER_TEXT_BYTES,
};
use fs_euler_disc_e2e::protocol::{MAX_CONTRACT_CHECK_RECEIPT_BYTES, MAX_EVIDENCE_PACKET_BYTES};
use fs_euler_disc_e2e::{
    AssessmentDisposition, AuthorityCeiling, CLAIM_POLICY_ASSESSMENT_LOG_DOMAIN,
    ClaimEvidencePacket, ClaimPolicyAssessment, ClaimPolicyAssessmentLog, ContractCheckReceipt,
    ContractError, ContractIdentity, DeclaredEvidenceAccessClass, EULER_ASSESSMENT_IDENTITY_DOMAIN,
    EULER_CLAIM_GRAPH_IDENTITY_DOMAIN, EULER_CLAIM_POLICY_SCHEMA_VERSION,
    EULER_CONTRACT_IDENTITY_DOMAIN, EULER_CONTRACT_SCHEMA_VERSION,
    EULER_EVIDENCE_PACKET_IDENTITY_DOMAIN, EULER_OWNER_MATRIX_IDENTITY_DOMAIN,
    EULER_OWNER_MATRIX_SCHEMA_VERSION, EULER_PREREQUISITE_RECEIPT_IDENTITY_DOMAIN,
    EULER_PROTOCOL_SCHEMA_VERSION, EulerClaimGraph, EulerClaimKind, EulerClaimSpec,
    EulerContextExtension, EulerScientificContract, EvidenceAuthorityClass,
    EvidenceAuthorityDeclaration, EvidenceRecord, EvidenceRequirement, FROZEN_CLAIM_GRAPH_HASH_HEX,
    FROZEN_CONTEXT_HASH_HEX, FROZEN_CONTRACT_IDENTITY_HEX, HYPOTHESIS_SOURCE_DECLARATION_DOMAIN,
    HypothesisSource, MAX_EULER_NO_CLAIMS, MAX_OWNER_MATRIX_BYTES, MAX_PROTOCOL_ID_BYTES,
    MAX_VALIDITY_DOMAIN_AXES, MAX_VALIDITY_DOMAIN_CANONICAL_BYTES, OwnerMatrix, OwnerRole,
    OwnerRow, ProtocolBudget, ProtocolSeed, ReportedScientificDisposition, ScientificRisk,
    admit_frozen_contract, build_frozen_contract, check_frozen_contract,
};
use fs_evidence::vv::{
    AcceptanceCriterion, ApplicabilityDomain, ApplicabilityPoint, CategoricalDomainAxis,
    ContextOfUse, NumericDomainAxis, QoiSpec, UnitId,
};
use fs_evidence::{Color, ValidityDomain};
use fs_govern::evidence_contract::NoClaimBoundary;
use fs_ir::campaign::{
    CampaignClaim, CampaignClaimId, ClaimDependency, EvidenceGapId, EvidenceUse,
};

fn hash(label: &str) -> ContentHash {
    fs_blake3::hash_domain(
        "org.frankensim.fs-euler-disc-e2e.test-artifact.v1",
        label.as_bytes(),
    )
}

fn point(contract: &EulerScientificContract) -> ApplicabilityPoint {
    point_at_fraction(contract, 0.5)
}

fn point_at_fraction(contract: &EulerScientificContract, fraction: f64) -> ApplicabilityPoint {
    let numeric = contract
        .context()
        .applicability()
        .numeric()
        .iter()
        .map(|(axis, domain)| {
            let (lo, hi) = domain.bounds();
            (axis.clone(), lo + (hi - lo) * fraction)
        })
        .collect();
    let categorical = contract
        .context()
        .applicability()
        .categorical()
        .iter()
        .map(|(axis, domain)| {
            (
                axis.clone(),
                domain
                    .allowed()
                    .iter()
                    .next()
                    .expect("nonempty categorical domain")
                    .clone(),
            )
        })
        .collect();
    ApplicabilityPoint::try_new(numeric, categorical).expect("in-domain point")
}

fn covering_regime(contract: &EulerScientificContract) -> ValidityDomain {
    contract.context().applicability().numeric().iter().fold(
        ValidityDomain::unconstrained(),
        |regime, (axis, domain)| {
            let (lo, hi) = domain.bounds();
            regime.with(axis.as_str(), lo, hi)
        },
    )
}

fn validity_domain_with_exact_canonical_bytes(
    axis_count: usize,
    target_bytes: usize,
) -> ValidityDomain {
    // Shared Color v2 encodes the domain as one u64 count and, per row,
    // u64-length + axis bytes + two (u64-length + f64-bits) fields.
    let fixed_bytes = size_of::<u64>() + axis_count * 40;
    let prefixes = (0..axis_count)
        .map(|index| format!("axis-{index:02}-"))
        .collect::<Vec<_>>();
    let minimum_name_bytes = prefixes.iter().map(String::len).sum::<usize>();
    let target_name_bytes = target_bytes
        .checked_sub(fixed_bytes)
        .expect("target must cover fixed validity-domain framing");
    assert!(target_name_bytes >= minimum_name_bytes);

    let maximum_name_bytes = axis_count * fs_evidence::MAX_COLOR_IDENTITY_BYTES;
    assert!(target_name_bytes <= maximum_name_bytes);
    let mut remaining_extra = target_name_bytes - minimum_name_bytes;
    let mut regime = ValidityDomain::unconstrained();
    for prefix in prefixes {
        let extra_capacity = fs_evidence::MAX_COLOR_IDENTITY_BYTES - prefix.len();
        let extra = remaining_extra.min(extra_capacity);
        remaining_extra -= extra;
        regime = regime.with(format!("{prefix}{}", "x".repeat(extra)), 0.0, 1.0);
    }
    assert_eq!(
        remaining_extra, 0,
        "requested exact domain size is feasible"
    );

    let observed_bytes = size_of::<u64>()
        + regime
            .bounds()
            .keys()
            .map(|axis| 40 + axis.len())
            .sum::<usize>();
    assert_eq!(observed_bytes, target_bytes);
    regime
}

fn declared_access_class(requirement: EvidenceRequirement) -> DeclaredEvidenceAccessClass {
    match requirement {
        EvidenceRequirement::CalibrationPartition => DeclaredEvidenceAccessClass::Calibration,
        EvidenceRequirement::PhysicalValidation
        | EvidenceRequirement::RivalMechanismDiscrimination => {
            DeclaredEvidenceAccessClass::Validation
        }
        EvidenceRequirement::BlindHoldout => DeclaredEvidenceAccessClass::BlindHoldout,
        _ => DeclaredEvidenceAccessClass::NotApplicable,
    }
}

fn authority(
    contract: &EulerScientificContract,
    claim: EulerClaimKind,
    requirement: EvidenceRequirement,
    weak: bool,
) -> EvidenceAuthorityDeclaration {
    match (requirement.authority_class(), weak) {
        (EvidenceAuthorityClass::StructuralProcess, false) => {
            EvidenceAuthorityDeclaration::StructuralProcess {
                receipt_hash: hash(&format!("process:{}:{}", claim.id(), requirement.code())),
            }
        }
        (
            EvidenceAuthorityClass::StructuralProcess | EvidenceAuthorityClass::VerifiedNumerics,
            true,
        ) => EvidenceAuthorityDeclaration::VerifiedNumerics {
            color: Color::Estimated {
                estimator: format!("weak-{}-{}", claim.id(), requirement.code()),
                dispersion: 1.0,
            },
        },
        (EvidenceAuthorityClass::VerifiedNumerics, false) => {
            EvidenceAuthorityDeclaration::VerifiedNumerics {
                color: Color::Verified { lo: 0.0, hi: 0.0 },
            }
        }
        (EvidenceAuthorityClass::ValidatedPhysical, false) => {
            EvidenceAuthorityDeclaration::ValidatedPhysical {
                color: Color::Validated {
                    regime: covering_regime(contract),
                    dataset: format!("dataset-{}-{}", claim.id(), requirement.code()),
                },
            }
        }
        (EvidenceAuthorityClass::ValidatedPhysical, true) => {
            EvidenceAuthorityDeclaration::ValidatedPhysical {
                color: Color::Estimated {
                    estimator: format!("weak-{}-{}", claim.id(), requirement.code()),
                    dispersion: 1.0,
                },
            }
        }
    }
}

fn record_with(
    contract: &EulerScientificContract,
    claim: EulerClaimKind,
    requirement: EvidenceRequirement,
    weak: bool,
    artifact_hash: ContentHash,
    access_class: DeclaredEvidenceAccessClass,
) -> EvidenceRecord {
    let qois = contract
        .claim_graph()
        .claim(claim)
        .expect("frozen claim")
        .campaign()
        .qois
        .clone();
    EvidenceRecord::try_new(
        contract.identity(),
        claim,
        requirement,
        qois,
        authority(contract, claim, requirement, weak),
        artifact_hash,
        format!("source-{}-{}", claim.id(), requirement.code()),
        requirement.source_schema(),
        requirement.source_kind(),
        hash(&format!(
            "schema-admission:{}:{}",
            claim.id(),
            requirement.code()
        )),
        access_class,
        true,
    )
    .expect("evidence record")
}

fn records(contract: &EulerScientificContract, claim: EulerClaimKind) -> Vec<EvidenceRecord> {
    contract
        .claim_graph()
        .claim(claim)
        .expect("frozen claim")
        .requirements()
        .iter()
        .map(|requirement| {
            record_with(
                contract,
                claim,
                *requirement,
                false,
                hash(&format!("{}:{}", claim.id(), requirement.code())),
                declared_access_class(*requirement),
            )
        })
        .collect()
}

fn validated_record_with_regime(
    contract: &EulerScientificContract,
    regime: ValidityDomain,
    dataset: &str,
    label: &str,
) -> Result<EvidenceRecord, ContractError> {
    let claim = EulerClaimKind::BlindTrajectoryPrediction;
    let requirement = EvidenceRequirement::PhysicalValidation;
    let qois = contract
        .claim_graph()
        .claim(claim)
        .expect("frozen blind-prediction claim")
        .campaign()
        .qois
        .clone();
    EvidenceRecord::try_new(
        contract.identity(),
        claim,
        requirement,
        qois,
        EvidenceAuthorityDeclaration::ValidatedPhysical {
            color: Color::Validated {
                regime,
                dataset: dataset.to_owned(),
            },
        },
        hash(&format!("validity-artifact-{label}")),
        format!("validity-source-{label}"),
        requirement.source_schema(),
        requirement.source_kind(),
        hash(&format!("validity-schema-receipt-{label}")),
        DeclaredEvidenceAccessClass::Validation,
        true,
    )
}

fn packet(
    contract: &EulerScientificContract,
    claim: EulerClaimKind,
    evidence: Vec<EvidenceRecord>,
    target_fitted: bool,
    scientific: ReportedScientificDisposition,
    expected: AssessmentDisposition,
) -> ClaimEvidencePacket {
    packet_at_point(
        contract,
        claim,
        point(contract),
        evidence,
        target_fitted,
        scientific,
        expected,
    )
}

#[allow(clippy::too_many_arguments)]
fn packet_at_point(
    contract: &EulerScientificContract,
    claim: EulerClaimKind,
    applicability_point: ApplicabilityPoint,
    evidence: Vec<EvidenceRecord>,
    target_fitted: bool,
    scientific: ReportedScientificDisposition,
    expected: AssessmentDisposition,
) -> ClaimEvidencePacket {
    let units = claim_units(contract, claim);
    packet_with_policy_declarations(
        contract,
        claim,
        applicability_point,
        evidence,
        true,
        target_fitted,
        scientific,
        expected,
        units,
    )
}

fn claim_units(contract: &EulerScientificContract, claim: EulerClaimKind) -> Vec<String> {
    contract
        .claim_graph()
        .claim(claim)
        .expect("frozen claim")
        .campaign()
        .qois
        .iter()
        .map(|id| contract.context().qois()[id].unit().as_str().to_owned())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect()
}

#[allow(clippy::too_many_arguments)]
fn packet_with_policy_declarations(
    contract: &EulerScientificContract,
    claim: EulerClaimKind,
    applicability_point: ApplicabilityPoint,
    evidence: Vec<EvidenceRecord>,
    no_claims_accepted: bool,
    target_fitted: bool,
    scientific: ReportedScientificDisposition,
    expected: AssessmentDisposition,
    units: Vec<String>,
) -> ClaimEvidencePacket {
    packet_from_parts(
        contract.identity(),
        format!("case-{}-{}", claim.id(), scientific.code()),
        claim,
        applicability_point,
        evidence,
        no_claims_accepted,
        target_fitted,
        scientific,
        expected,
        units,
        ProtocolSeed::not_applicable("deterministic-fixture").expect("seed"),
        ProtocolBudget::try_new(60_000, 64 * 1024 * 1024, 0.05).expect("budget"),
    )
}

#[allow(clippy::too_many_arguments)]
fn packet_from_parts(
    contract_identity: ContractIdentity,
    case_id: impl Into<String>,
    claim: EulerClaimKind,
    applicability_point: ApplicabilityPoint,
    evidence: Vec<EvidenceRecord>,
    no_claims_accepted: bool,
    target_fitted: bool,
    scientific: ReportedScientificDisposition,
    expected: AssessmentDisposition,
    units: Vec<String>,
    seed: ProtocolSeed,
    budget: ProtocolBudget,
) -> ClaimEvidencePacket {
    ClaimEvidencePacket::try_new(
        contract_identity,
        case_id,
        hash("shared-design-set"),
        hash("shared-aggregate-qoi-derivation-receipt"),
        claim,
        applicability_point,
        evidence,
        no_claims_accepted,
        target_fitted,
        scientific,
        expected,
        units,
        seed,
        budget,
    )
    .expect("packet")
}

const TOPOLOGICAL_CLAIM_ORDER: [EulerClaimKind; 9] = [
    EulerClaimKind::NumericalTrajectoryVerification,
    EulerClaimKind::CalibratedReproduction,
    EulerClaimKind::BlindTrajectoryPrediction,
    EulerClaimKind::EventOrCrossoverPrediction,
    EulerClaimKind::QualitativeEffectDirection,
    EulerClaimKind::Ranking,
    EulerClaimKind::NonlinearOptimumInterval,
    EulerClaimKind::EnergyChannelAttribution,
    EulerClaimKind::MechanismAttribution,
];

fn assess_with_frozen_prerequisites(
    contract: &EulerScientificContract,
    target: ClaimEvidencePacket,
) -> ClaimPolicyAssessment {
    let admitted = admit_frozen_contract(contract.clone()).expect("exact frozen admission");
    let target_claim = target.claim();
    let shared_point = target.point().clone();
    let mut target = Some(target);
    let mut prior = BTreeMap::<EulerClaimKind, ClaimPolicyAssessment>::new();
    for kind in TOPOLOGICAL_CLAIM_ORDER {
        let prerequisites = contract
            .claim_graph()
            .dependencies()
            .iter()
            .filter(|dependency| dependency.dependent.as_str() == kind.id())
            .map(|dependency| {
                let prerequisite = EulerClaimKind::ALL
                    .into_iter()
                    .find(|candidate| candidate.id() == dependency.prerequisite.as_str())
                    .expect("known prerequisite");
                prior[&prerequisite]
                    .as_prerequisite_for(kind, dependency.use_kind)
                    .expect("positive prerequisite receipt")
            })
            .collect::<Vec<_>>();
        let current = if kind == target_claim {
            target.take().expect("target is consumed exactly once")
        } else {
            packet_at_point(
                contract,
                kind,
                shared_point.clone(),
                records(contract, kind),
                false,
                ReportedScientificDisposition::Positive,
                AssessmentDisposition::ReferenceCompleteCandidate,
            )
        };
        let assessment = current
            .assess(&admitted, &prerequisites)
            .expect("structural assessment");
        if kind == target_claim {
            return assessment;
        }
        assert_eq!(
            assessment.disposition(),
            AssessmentDisposition::ReferenceCompleteCandidate,
            "upstream {}: {:?}",
            kind.id(),
            assessment.reasons()
        );
        prior.insert(kind, assessment);
    }
    panic!("target claim missing from topological order")
}

fn synthetic_reported_positive_assessment_map(
    contract: &EulerScientificContract,
) -> BTreeMap<EulerClaimKind, ClaimPolicyAssessment> {
    let admitted = admit_frozen_contract(contract.clone()).expect("exact frozen admission");
    let mut prior = BTreeMap::<EulerClaimKind, ClaimPolicyAssessment>::new();
    for kind in TOPOLOGICAL_CLAIM_ORDER {
        let prerequisites = contract
            .claim_graph()
            .dependencies()
            .iter()
            .filter(|dependency| dependency.dependent.as_str() == kind.id())
            .map(|dependency| {
                let prerequisite = EulerClaimKind::ALL
                    .into_iter()
                    .find(|candidate| candidate.id() == dependency.prerequisite.as_str())
                    .expect("known prerequisite");
                prior[&prerequisite]
                    .as_prerequisite_for(kind, dependency.use_kind)
                    .expect("positive prerequisite receipt")
            })
            .collect::<Vec<_>>();
        let assessment = packet(
            contract,
            kind,
            records(contract, kind),
            false,
            ReportedScientificDisposition::Positive,
            AssessmentDisposition::ReferenceCompleteCandidate,
        )
        .assess(&admitted, &prerequisites)
        .expect("positive structural assessment");
        assert_eq!(
            assessment.disposition(),
            AssessmentDisposition::ReferenceCompleteCandidate,
            "{}: {:?}",
            kind.id(),
            assessment.reasons()
        );
        prior.insert(kind, assessment);
    }
    prior
}

fn golden_requirements(kind: EulerClaimKind) -> &'static [EvidenceRequirement] {
    use EvidenceRequirement as E;
    match kind {
        EulerClaimKind::NumericalTrajectoryVerification => &[
            E::CodeVerification,
            E::SolutionVerification,
            E::IndependentReconstruction,
        ],
        EulerClaimKind::CalibratedReproduction => &[
            E::CodeVerification,
            E::SolutionVerification,
            E::CalibrationPartition,
            E::ApplicabilityCheck,
            E::UncertaintyClosure,
            E::IndependentReconstruction,
        ],
        EulerClaimKind::BlindTrajectoryPrediction
        | EulerClaimKind::EventOrCrossoverPrediction
        | EulerClaimKind::QualitativeEffectDirection
        | EulerClaimKind::Ranking
        | EulerClaimKind::NonlinearOptimumInterval => &[
            E::CodeVerification,
            E::SolutionVerification,
            E::PhysicalValidation,
            E::BlindHoldout,
            E::PreregisteredAnalysis,
            E::ApplicabilityCheck,
            E::UncertaintyClosure,
            E::MultiplicityControl,
            E::IndependentReconstruction,
        ],
        EulerClaimKind::EnergyChannelAttribution => &[
            E::CodeVerification,
            E::SolutionVerification,
            E::PhysicalValidation,
            E::BlindHoldout,
            E::PreregisteredAnalysis,
            E::ApplicabilityCheck,
            E::UncertaintyClosure,
            E::MultiplicityControl,
            E::EnergyBalanceClosure,
            E::IndependentReconstruction,
        ],
        EulerClaimKind::MechanismAttribution => &[
            E::CodeVerification,
            E::SolutionVerification,
            E::PhysicalValidation,
            E::BlindHoldout,
            E::PreregisteredAnalysis,
            E::ApplicabilityCheck,
            E::UncertaintyClosure,
            E::MultiplicityControl,
            E::EnergyBalanceClosure,
            E::IndependentReconstruction,
            E::RivalMechanismDiscrimination,
        ],
    }
}

fn golden_qois(kind: EulerClaimKind) -> &'static [&'static str] {
    match kind {
        EulerClaimKind::NumericalTrajectoryVerification => &["numerical-trajectory-error"],
        EulerClaimKind::CalibratedReproduction | EulerClaimKind::BlindTrajectoryPrediction => {
            &["normalized-trajectory-discrepancy"]
        }
        EulerClaimKind::EventOrCrossoverPrediction => {
            &["event-class-disposition", "event-time-error"]
        }
        EulerClaimKind::QualitativeEffectDirection => &["qualitative-effect-disposition"],
        EulerClaimKind::Ranking => &["configuration-ranking-disposition"],
        EulerClaimKind::NonlinearOptimumInterval => {
            &["optimum-containment-disposition", "optimum-interval-width"]
        }
        EulerClaimKind::EnergyChannelAttribution => {
            &["energy-balance-residual", "energy-channel-fraction-error"]
        }
        EulerClaimKind::MechanismAttribution => &[
            "energy-channel-fraction-error",
            "rival-mechanism-disposition",
        ],
    }
}

#[test]
fn g0_frozen_policy_matches_an_independent_nine_claim_oracle() {
    let contract = build_frozen_contract().expect("frozen contract");
    assert_eq!(contract.claim_graph().claims().len(), 9);
    for kind in EulerClaimKind::ALL {
        let claim = contract.claim_graph().claim(kind).expect("claim");
        assert_eq!(claim.campaign().id.as_str(), kind.id());
        assert_eq!(claim.requirements(), golden_requirements(kind));
        assert_eq!(
            claim
                .campaign()
                .qois
                .iter()
                .map(fs_evidence::vv::QoiId::as_str)
                .collect::<Vec<_>>(),
            golden_qois(kind)
        );
    }

    let expected_dependencies = [
        (
            "blind-trajectory-prediction",
            "energy-channel-attribution",
            EvidenceUse::ValidationInput,
        ),
        (
            "blind-trajectory-prediction",
            "event-or-crossover-prediction",
            EvidenceUse::ValidationInput,
        ),
        (
            "blind-trajectory-prediction",
            "mechanism-attribution",
            EvidenceUse::ValidationInput,
        ),
        (
            "blind-trajectory-prediction",
            "nonlinear-optimum-interval",
            EvidenceUse::ValidationInput,
        ),
        (
            "blind-trajectory-prediction",
            "qualitative-effect-direction",
            EvidenceUse::ValidationInput,
        ),
        (
            "blind-trajectory-prediction",
            "ranking",
            EvidenceUse::ValidationInput,
        ),
        (
            "calibrated-reproduction",
            "blind-trajectory-prediction",
            EvidenceUse::CalibrationInput,
        ),
        (
            "energy-channel-attribution",
            "mechanism-attribution",
            EvidenceUse::ValidationInput,
        ),
        (
            "numerical-trajectory-verification",
            "blind-trajectory-prediction",
            EvidenceUse::ValidationInput,
        ),
        (
            "numerical-trajectory-verification",
            "calibrated-reproduction",
            EvidenceUse::ValidationInput,
        ),
    ];
    let actual_dependencies = contract
        .claim_graph()
        .dependencies()
        .iter()
        .map(|edge| {
            (
                edge.prerequisite.as_str(),
                edge.dependent.as_str(),
                edge.use_kind,
            )
        })
        .collect::<Vec<_>>();
    assert_eq!(actual_dependencies, expected_dependencies);

    let exact_no_claims = [
        "A successful blind case is local to its exact declared Context of Use and applicability domain.",
        "Agreement in an exponent, event time, or stop time does not identify an energy-loss mechanism.",
        "Deterministic software verification does not establish physical validation.",
        "Fitting or selecting against protected target outcomes is calibrated reproduction, not emergent prediction.",
        "Geometric similarity does not establish dynamic similarity across scale, material, support, or environment.",
        "Negative and inconclusive results are retained terminal outcomes and are never erased or promoted.",
        "Transcript and publication sources generate hypotheses only; they are not validation evidence.",
    ];
    assert_eq!(
        contract
            .no_claims()
            .entries()
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>(),
        exact_no_claims
    );

    let owner_oracle = [
        (
            OwnerRole::ContextOfUse,
            "fs-evidence",
            "org.frankensim.fs-evidence.vv-artifact.v3",
            AuthorityCeiling::StructuralContextDeclaration,
        ),
        (
            OwnerRole::VvEvidenceArtifact,
            "fs-evidence",
            "org.frankensim.fs-evidence.vv-artifact.v3",
            AuthorityCeiling::StructuralEvidenceReferenceOnly,
        ),
        (
            OwnerRole::VvSchemaAdmissionReceipt,
            "fs-evidence",
            "org.frankensim.fs-evidence.vv-schema-admission-receipt.v2",
            AuthorityCeiling::StructuralSchemaAdmissionOnly,
        ),
        (
            OwnerRole::CampaignClaimVocabulary,
            "fs-ir",
            "fs-ir:experiment-campaign-schema-v1",
            AuthorityCeiling::CampaignVocabularyOnly,
        ),
        (
            OwnerRole::NoClaimBoundary,
            "fs-govern",
            "fs-govern:authority-algebra-v2",
            AuthorityCeiling::CanonicalNoClaimBoundary,
        ),
        (
            OwnerRole::HypothesisSourceDeclaration,
            "fs-euler-disc-e2e",
            "org.frankensim.fs-euler-disc-e2e.hypothesis-source-declaration.v1",
            AuthorityCeiling::HypothesisOnly,
        ),
        (
            OwnerRole::EulerClaimGraph,
            "fs-euler-disc-e2e",
            "org.frankensim.fs-euler-disc-e2e.claim-graph.v1",
            AuthorityCeiling::StructuralClaimPolicyOnly,
        ),
        (
            OwnerRole::EulerScientificContract,
            "fs-euler-disc-e2e",
            "org.frankensim.fs-euler-disc-e2e.scientific-contract.v1",
            AuthorityCeiling::CandidateEligibilityOnly,
        ),
        (
            OwnerRole::ClaimEvidencePacket,
            "fs-euler-disc-e2e",
            "org.frankensim.fs-euler-disc-e2e.claim-evidence-packet.v1",
            AuthorityCeiling::CandidateEligibilityOnly,
        ),
        (
            OwnerRole::PrerequisiteAssessmentReceipt,
            "fs-euler-disc-e2e",
            "org.frankensim.fs-euler-disc-e2e.prerequisite-assessment-receipt.v1",
            AuthorityCeiling::StructuralDependencyOnly,
        ),
        (
            OwnerRole::ClaimPolicyAssessment,
            "fs-euler-disc-e2e",
            "org.frankensim.fs-euler-disc-e2e.claim-policy-assessment.v1",
            AuthorityCeiling::CandidateEligibilityOnly,
        ),
        (
            OwnerRole::ContractCheckReceipt,
            "fs-euler-disc-e2e",
            "org.frankensim.fs-euler-disc-e2e.contract-check-receipt.v1",
            AuthorityCeiling::StructuralCheckOnly,
        ),
        (
            OwnerRole::ClaimPolicyAssessmentLog,
            "fs-euler-disc-e2e",
            "org.frankensim.fs-euler-disc-e2e.claim-policy-assessment-log.v1",
            AuthorityCeiling::DiagnosticRetentionOnly,
        ),
        (
            OwnerRole::OwnerMatrixRegistry,
            "fs-euler-disc-e2e",
            "org.frankensim.fs-euler-disc-e2e.owner-matrix.v1",
            AuthorityCeiling::StructuralRoutingRegistryOnly,
        ),
        (
            OwnerRole::AggregateQoiDerivationReceipt,
            "fs-euler-disc-e2e",
            "org.frankensim.fs-euler-disc-e2e.aggregate-qoi-derivation-receipt.v1",
            AuthorityCeiling::StructuralEvidenceReferenceOnly,
        ),
    ];
    assert_eq!(contract.owner_matrix().rows().len(), owner_oracle.len());
    for (role, owner, schema, ceiling) in owner_oracle {
        let row = &contract.owner_matrix().rows()[&role];
        assert_eq!(row.owner_crate(), owner);
        assert_eq!(row.source_schema(), schema);
        assert_eq!(row.authority_ceiling(), ceiling);
    }
}

#[test]
fn owner_matrix_transport_is_exact_versioned_and_fail_closed() {
    let contract = build_frozen_contract().expect("frozen contract");
    let matrix = contract.owner_matrix();
    assert_eq!(matrix.schema_version(), EULER_OWNER_MATRIX_SCHEMA_VERSION);
    let bytes = matrix.canonical_bytes().expect("owner-matrix bytes");
    assert!(bytes.len() <= MAX_OWNER_MATRIX_BYTES);
    assert_eq!(
        matrix.identity().as_hash(),
        fs_blake3::hash_domain(EULER_OWNER_MATRIX_IDENTITY_DOMAIN, &bytes)
    );
    assert_eq!(
        OwnerMatrix::from_canonical_bytes(&bytes).expect("owner-matrix fixed point"),
        *matrix
    );
    matrix.verify_identity().expect("owner-matrix identity");

    let contract_bytes = contract.canonical_bytes().expect("contract bytes");
    assert!(
        contract_bytes
            .windows(matrix.identity().as_hash().as_bytes().len())
            .any(|window| window == matrix.identity().as_hash().as_bytes()),
        "the composite contract must embed the owner-matrix identity"
    );
    assert!(
        contract_bytes
            .windows(bytes.len())
            .any(|window| window == bytes.as_slice()),
        "the composite contract must embed the exact owner-matrix transport"
    );
    let identity_offset = contract_bytes
        .windows(matrix.identity().as_hash().as_bytes().len())
        .position(|window| window == matrix.identity().as_hash().as_bytes())
        .expect("embedded owner-matrix identity offset");
    let mut mismatched_contract = contract_bytes.clone();
    mismatched_contract[identity_offset] ^= 1;
    let error = EulerScientificContract::from_canonical_bytes(&mismatched_contract)
        .expect_err("contract must refuse owner-matrix identity/bytes mismatch");
    assert_eq!(error.code(), "EulerContractOwnerMatrixIdentityMismatch");

    let mut wrong_version = bytes.clone();
    wrong_version[8..12].copy_from_slice(&(EULER_OWNER_MATRIX_SCHEMA_VERSION + 1).to_le_bytes());
    let error = OwnerMatrix::from_canonical_bytes(&wrong_version)
        .expect_err("unknown owner-matrix version must refuse");
    assert_eq!(error.code(), "EulerOwnerMatrixUnsupportedVersion");
    for version in [0, 2, u32::MAX] {
        assert!(OwnerMatrix::migration_policy(version).is_err());
    }

    let alternate_ceiling = |expected| {
        if expected == AuthorityCeiling::HypothesisOnly {
            AuthorityCeiling::StructuralCheckOnly
        } else {
            AuthorityCeiling::HypothesisOnly
        }
    };
    for role in OwnerRole::ALL {
        let wrong_owner = OwnerRow::try_new(
            role,
            format!("{}-counterfeit", role.expected_owner_crate()),
            role.expected_source_schema(),
            role.expected_authority_ceiling(),
        )
        .expect_err("every owner substitution must refuse");
        assert_eq!(wrong_owner.code(), "EulerContractGenericSchemaFork");

        let wrong_schema = OwnerRow::try_new(
            role,
            role.expected_owner_crate(),
            format!("{}-counterfeit", role.expected_source_schema()),
            role.expected_authority_ceiling(),
        )
        .expect_err("every opaque routing-address substitution must refuse");
        assert_eq!(wrong_schema.code(), "EulerContractGenericSchemaFork");

        let wrong_ceiling = OwnerRow::try_new(
            role,
            role.expected_owner_crate(),
            role.expected_source_schema(),
            alternate_ceiling(role.expected_authority_ceiling()),
        )
        .expect_err("every authority-ceiling substitution must refuse");
        assert_eq!(wrong_ceiling.code(), "EulerContractGenericSchemaFork");
    }
}

#[test]
fn g0_public_hostile_input_byte_caps_have_exact_and_maximum_plus_one_boundaries() {
    let exact_text = "x".repeat(MAX_EULER_TEXT_BYTES);
    HypothesisSource::try_new(exact_text, "bounded-locator")
        .expect("the exact public text-byte maximum must be admitted");
    let error = HypothesisSource::try_new("x".repeat(MAX_EULER_TEXT_BYTES + 1), "bounded-locator")
        .expect_err("maximum-plus-one public text bytes must refuse");
    assert_eq!(error.code(), "EulerContractInvalidText");

    let mut owner_matrix_bytes = vec![0_u8; MAX_OWNER_MATRIX_BYTES];
    let error = OwnerMatrix::from_canonical_bytes(&owner_matrix_bytes)
        .expect_err("exact-limit hostile owner-matrix bytes must reach structural decoding");
    assert_eq!(error.code(), "EulerOwnerMatrixMalformedCanonical");
    owner_matrix_bytes.push(0);
    let error = OwnerMatrix::from_canonical_bytes(&owner_matrix_bytes)
        .expect_err("maximum-plus-one owner-matrix bytes must refuse at preflight");
    assert_eq!(error.code(), "EulerOwnerMatrixTooLarge");

    let mut contract_bytes = vec![0_u8; MAX_EULER_CONTRACT_BYTES];
    let error = EulerScientificContract::from_canonical_bytes(&contract_bytes)
        .expect_err("exact-limit hostile contract bytes must reach structural decoding");
    assert_eq!(error.code(), "EulerContractMalformedCanonical");
    contract_bytes.push(0);
    let error = EulerScientificContract::from_canonical_bytes(&contract_bytes)
        .expect_err("maximum-plus-one contract bytes must refuse at preflight");
    assert_eq!(error.code(), "EulerContractTooLarge");

    let mut receipt_bytes = vec![0_u8; MAX_CONTRACT_CHECK_RECEIPT_BYTES];
    let error = ContractCheckReceipt::from_canonical_bytes(&receipt_bytes)
        .expect_err("exact-limit hostile receipt bytes must reach structural decoding");
    assert_eq!(error.code(), "EulerContractCheckReceiptMagic");
    receipt_bytes.push(0);
    let error = ContractCheckReceipt::from_canonical_bytes(&receipt_bytes)
        .expect_err("maximum-plus-one receipt bytes must refuse at preflight");
    assert_eq!(error.code(), "EulerContractCheckReceiptTooLarge");
}

#[test]
fn g0_public_packet_byte_cap_preflight_is_exact_and_bounds_oversized_refusals() {
    const LENGTH_FRAME_BYTES: usize = size_of::<u32>();
    const MAX_POINT_ID_BYTES: usize = fs_evidence::vv::MAX_VV_ID_BYTES;
    const MAX_CATEGORICAL_ROW_BYTES: usize = 2 * LENGTH_FRAME_BYTES + 2 * MAX_POINT_ID_BYTES;
    const MIN_FIXED_AXIS_ROW_BYTES: usize = 2 * LENGTH_FRAME_BYTES + MAX_POINT_ID_BYTES + 1;

    let contract = build_frozen_contract().expect("frozen contract");
    let claim = EulerClaimKind::BlindTrajectoryPrediction;
    // Keep one Validated Color in the packet so the exact-size preflight must
    // include the nested Color-v2 dataset, regime, field frames, and local
    // authority tag rather than accounting only for point and text rows.
    let physical_record = records(&contract, claim)
        .into_iter()
        .find(|record| record.requirement() == EvidenceRequirement::PhysicalValidation)
        .expect("validated physical record");
    let empty_point = ApplicabilityPoint::try_new(Vec::new(), Vec::new()).expect("empty point");
    let construct = |point| {
        ClaimEvidencePacket::try_new(
            contract.identity(),
            "packet-byte-cap-boundary",
            hash("packet-boundary-design-set"),
            hash("packet-boundary-aggregate-qoi-derivation-receipt"),
            claim,
            point,
            vec![physical_record.clone()],
            true,
            false,
            ReportedScientificDisposition::Positive,
            AssessmentDisposition::ReferenceCompleteCandidate,
            claim_units(&contract, claim),
            ProtocolSeed::Fixed { value: 7 },
            ProtocolBudget::try_new(1, 1, 0.0).expect("budget"),
        )
    };
    let base = construct(empty_point).expect("bounded base packet");
    let base_len = base.canonical_bytes().expect("base packet bytes").len();
    let exact_point_row_bytes = MAX_EVIDENCE_PACKET_BYTES
        .checked_sub(base_len)
        .expect("the base packet is below the packet byte cap");

    let point_with_exact_categorical_row_bytes = |target_row_bytes: usize| {
        let row_count = target_row_bytes.div_ceil(MAX_CATEGORICAL_ROW_BYTES);
        assert!(row_count <= fs_evidence::vv::MAX_VV_ITEMS);
        assert!(
            target_row_bytes >= row_count * MIN_FIXED_AXIS_ROW_BYTES,
            "the requested packet-boundary payload is representable"
        );
        let mut remaining_extra = target_row_bytes - row_count * MIN_FIXED_AXIS_ROW_BYTES;
        let mut categorical = Vec::with_capacity(row_count);
        for index in 0..row_count {
            let prefix = format!("packet-cap-axis-{index:04}-");
            assert!(prefix.len() <= MAX_POINT_ID_BYTES);
            let axis = format!("{prefix}{}", "a".repeat(MAX_POINT_ID_BYTES - prefix.len()));
            let extra = remaining_extra.min(MAX_POINT_ID_BYTES - 1);
            remaining_extra -= extra;
            categorical.push((
                fs_evidence::vv::AxisId::try_new(axis).expect("bounded unique packet axis"),
                "v".repeat(1 + extra),
            ));
        }
        assert_eq!(remaining_extra, 0, "the exact row budget was consumed");
        ApplicabilityPoint::try_new(Vec::new(), categorical).expect("bounded categorical point")
    };

    let exact = construct(point_with_exact_categorical_row_bytes(
        exact_point_row_bytes,
    ))
    .expect("the exact packet byte maximum must be admitted");
    assert_eq!(
        exact.canonical_bytes().expect("exact packet bytes").len(),
        MAX_EVIDENCE_PACKET_BYTES
    );
    exact
        .verify_identity()
        .expect("the exact-boundary packet identity must remain a fixed point");

    let maximum_plus_one = construct(point_with_exact_categorical_row_bytes(
        exact_point_row_bytes + 1,
    ))
    .expect_err("maximum-plus-one packet bytes must refuse during size preflight");
    let largest_generic_point = construct(point_with_exact_categorical_row_bytes(
        fs_evidence::vv::MAX_VV_ITEMS * MAX_CATEGORICAL_ROW_BYTES,
    ))
    .expect_err("the largest generic point must refuse with bounded packet diagnostics");
    for error in [&maximum_plus_one, &largest_generic_point] {
        assert_eq!(error.code(), "EulerProtocolPacketTooLarge");
        assert_eq!(
            error.detail(),
            "canonical evidence packet exceeds its byte budget"
        );
        assert!(error.detail().len() < MAX_PROTOCOL_ID_BYTES);
    }
}

#[test]
fn g0_public_evidence_gap_text_refusal_is_bounded_without_clone_amplification() {
    let contract = build_frozen_contract().expect("frozen contract");
    let kind = EulerClaimKind::NumericalTrajectoryVerification;
    let template = contract
        .claim_graph()
        .claim(kind)
        .expect("frozen claim with an evidence gap");
    assert!(!template.campaign().evidence_gaps.is_empty());

    let attempt = |expected_field: bool, value: String| {
        let mut campaign = template.campaign().clone();
        let gap = campaign
            .evidence_gaps
            .first_mut()
            .expect("frozen claim has an evidence gap");
        if expected_field {
            gap.expected_evidence = value;
        } else {
            gap.description = value;
        }
        EulerClaimSpec::try_new(kind, campaign, template.requirements().to_vec())
    };

    for expected_field in [true, false] {
        attempt(expected_field, "x".repeat(MAX_EULER_TEXT_BYTES))
            .expect("the exact evidence-gap text limit must remain admissible");
        let plus_one = attempt(expected_field, "x".repeat(MAX_EULER_TEXT_BYTES + 1))
            .expect_err("maximum-plus-one evidence-gap text must refuse");
        let very_large = attempt(expected_field, "x".repeat(1024 * 1024))
            .expect_err("very large evidence-gap text must refuse without an extra full clone");
        assert_eq!(plus_one.code(), "EulerContractInvalidText");
        assert_eq!(very_large.code(), "EulerContractInvalidText");
        assert_eq!(very_large.detail(), plus_one.detail());
        assert!(plus_one.detail().len() < MAX_EULER_TEXT_BYTES);
    }
}

#[test]
fn g0_negative_result_erasure_routes_only_to_terminal_retention() {
    let contract = build_frozen_contract().expect("frozen contract");
    let terminal = "retain-as-terminal-non-promotion";
    assert!(
        contract
            .extension()
            .decision_alternatives()
            .iter()
            .any(|alternative| alternative == terminal)
    );
    let risk = contract
        .extension()
        .risks()
        .iter()
        .find(|risk| risk.code() == "negative-result-erasure")
        .expect("negative-result-erasure risk");
    assert_eq!(risk.decision_alternative(), terminal);
    assert_ne!(
        risk.decision_alternative(),
        "advance-to-separate-candidate-review"
    );
    assert_eq!(risk.affected_claims(), EulerClaimKind::ALL);
}

#[test]
fn g0_every_claim_has_synthetic_structural_positive_missing_and_weakest_cases() {
    let contract = build_frozen_contract().expect("frozen contract");
    for kind in EulerClaimKind::ALL {
        let complete = records(&contract, kind);
        let positive = assess_with_frozen_prerequisites(
            &contract,
            packet(
                &contract,
                kind,
                complete.clone(),
                false,
                ReportedScientificDisposition::Positive,
                AssessmentDisposition::ReferenceCompleteCandidate,
            ),
        );
        assert_eq!(
            positive.disposition(),
            AssessmentDisposition::ReferenceCompleteCandidate,
            "{}: {:?}",
            kind.id(),
            positive.reasons()
        );

        for missing in golden_requirements(kind) {
            let without = complete
                .iter()
                .filter(|record| record.requirement() != *missing)
                .cloned()
                .collect();
            let assessment = assess_with_frozen_prerequisites(
                &contract,
                packet(
                    &contract,
                    kind,
                    without,
                    false,
                    ReportedScientificDisposition::Positive,
                    AssessmentDisposition::Refused,
                ),
            );
            assert_eq!(assessment.disposition(), AssessmentDisposition::Refused);
            assert!(
                assessment
                    .reasons()
                    .contains(&format!("missing-evidence:{}", missing.code()))
            );
        }

        for weakest in golden_requirements(kind) {
            let weakened = golden_requirements(kind)
                .iter()
                .map(|requirement| {
                    record_with(
                        &contract,
                        kind,
                        *requirement,
                        requirement == weakest,
                        hash(&format!("weak:{}:{}", kind.id(), requirement.code())),
                        declared_access_class(*requirement),
                    )
                })
                .collect();
            let assessment = assess_with_frozen_prerequisites(
                &contract,
                packet(
                    &contract,
                    kind,
                    weakened,
                    false,
                    ReportedScientificDisposition::Positive,
                    AssessmentDisposition::DemotedCandidate,
                ),
            );
            assert_eq!(
                assessment.disposition(),
                AssessmentDisposition::DemotedCandidate,
                "{} / {}: {:?}",
                kind.id(),
                weakest.code(),
                assessment.reasons()
            );
        }
    }
}

#[test]
fn g0_unexpected_evidence_binds_its_observed_authority_slot_shape() {
    let contract = build_frozen_contract().expect("frozen contract");
    let claim = EulerClaimKind::NumericalTrajectoryVerification;
    let requirement = EvidenceRequirement::PhysicalValidation;
    let claim_spec = contract.claim_graph().claim(claim).expect("frozen claim");
    assert!(!claim_spec.requirements().contains(&requirement));
    let qois = claim_spec.campaign().qois.clone();
    let role_receipt = hash("unexpected-structural-role-receipt");
    let unexpected = EvidenceRecord::try_new(
        contract.identity(),
        claim,
        requirement,
        qois,
        EvidenceAuthorityDeclaration::StructuralProcess {
            receipt_hash: role_receipt,
        },
        hash("unexpected-structural-artifact"),
        "unexpected-structural-source",
        requirement.source_schema(),
        requirement.source_kind(),
        hash("unexpected-structural-schema-receipt"),
        declared_access_class(requirement),
        true,
    )
    .expect("unexpected evidence row with mismatched authority class");
    let mut evidence = records(&contract, claim);
    evidence.push(unexpected);
    let assessment = assess_with_frozen_prerequisites(
        &contract,
        packet(
            &contract,
            claim,
            evidence,
            false,
            ReportedScientificDisposition::Positive,
            AssessmentDisposition::Refused,
        ),
    );
    assert_eq!(assessment.disposition(), AssessmentDisposition::Refused);
    assert!(
        assessment
            .reasons()
            .contains(&"unexpected-evidence:physical-validation".to_owned())
    );
    assert!(assessment.reasons().contains(
        &"weak-authority:physical-validation:requires-validated-physical:observed-structural-process"
            .to_owned()
    ));
    assert!(assessment.log().json_line().contains(&format!(
        "evidence:physical-validation:role-receipt:{}",
        role_receipt.to_hex()
    )));
    fs_euler_disc_e2e::ClaimPolicyAssessmentLog::from_json_line(assessment.log().json_line())
        .expect("unexpected mismatched-authority evidence log must round-trip");
}

#[test]
fn g0_unexpected_evidence_hypothesis_collisions_refuse_and_round_trip() {
    let contract = build_frozen_contract().expect("frozen contract");
    let claim = EulerClaimKind::NumericalTrajectoryVerification;
    let requirement = EvidenceRequirement::PhysicalValidation;
    let claim_spec = contract.claim_graph().claim(claim).expect("frozen claim");
    assert!(!claim_spec.requirements().contains(&requirement));
    let qois = claim_spec.campaign().qois.clone();
    let hypothesis_hash = contract
        .extension()
        .hypothesis_sources()
        .first()
        .expect("frozen hypothesis source")
        .declaration_hash();

    for slot in ["artifact", "schema-admission-receipt", "role-receipt"] {
        let selected_hash = |candidate: &str, ordinary_label: &str| {
            if slot == candidate {
                hypothesis_hash
            } else {
                hash(ordinary_label)
            }
        };
        let unexpected = EvidenceRecord::try_new(
            contract.identity(),
            claim,
            requirement,
            qois.clone(),
            EvidenceAuthorityDeclaration::StructuralProcess {
                receipt_hash: selected_hash(
                    "role-receipt",
                    "unexpected-collision-ordinary-role-receipt",
                ),
            },
            selected_hash("artifact", "unexpected-collision-ordinary-artifact"),
            format!("unexpected-collision-source-{slot}"),
            requirement.source_schema(),
            requirement.source_kind(),
            selected_hash(
                "schema-admission-receipt",
                "unexpected-collision-ordinary-schema-receipt",
            ),
            declared_access_class(requirement),
            true,
        )
        .expect("unexpected evidence row with one hypothesis collision");
        let mut evidence = records(&contract, claim);
        evidence.push(unexpected);
        let assessment = assess_with_frozen_prerequisites(
            &contract,
            packet(
                &contract,
                claim,
                evidence,
                false,
                ReportedScientificDisposition::Positive,
                AssessmentDisposition::Refused,
            ),
        );
        let collision_reason =
            format!("hypothesis-source-cannot-satisfy-evidence:physical-validation:{slot}");
        assert_eq!(assessment.disposition(), AssessmentDisposition::Refused);
        assert!(assessment.reasons().contains(&collision_reason));
        fs_euler_disc_e2e::ClaimPolicyAssessmentLog::from_json_line(assessment.log().json_line())
            .unwrap_or_else(|error| panic!("unexpected collision {slot} must round-trip: {error}"));
    }
}

#[test]
fn g0_direct_claim_dag_receipts_fail_closed_on_every_binding() {
    let contract = build_frozen_contract().expect("frozen contract");
    let admitted = admit_frozen_contract(contract.clone()).expect("exact frozen admission");
    let positive = synthetic_reported_positive_assessment_map(&contract);

    // Every one of the ten frozen direct edges is necessary. This loop is an
    // edge-by-edge negative oracle, not merely a test that a dependent has at
    // least one generic ancestor.
    for dependent in TOPOLOGICAL_CLAIM_ORDER {
        let edges = contract
            .claim_graph()
            .dependencies()
            .iter()
            .filter(|edge| edge.dependent.as_str() == dependent.id())
            .collect::<Vec<_>>();
        if edges.is_empty() {
            continue;
        }
        let receipts = edges
            .iter()
            .map(|edge| {
                let prerequisite = EulerClaimKind::ALL
                    .into_iter()
                    .find(|kind| kind.id() == edge.prerequisite.as_str())
                    .expect("known prerequisite");
                positive[&prerequisite]
                    .as_prerequisite_for(dependent, edge.use_kind)
                    .expect("exact direct receipt")
            })
            .collect::<Vec<_>>();
        let target = packet(
            &contract,
            dependent,
            records(&contract, dependent),
            false,
            ReportedScientificDisposition::Positive,
            AssessmentDisposition::ReferenceCompleteCandidate,
        );
        let baseline = target
            .assess(&admitted, &receipts)
            .expect("baseline assessment");
        assert_eq!(
            baseline.disposition(),
            AssessmentDisposition::ReferenceCompleteCandidate
        );
        for missing_index in 0..receipts.len() {
            let without = receipts
                .iter()
                .enumerate()
                .filter(|(index, _)| *index != missing_index)
                .map(|(_, receipt)| receipt.clone())
                .collect::<Vec<_>>();
            let refused = target
                .assess(&admitted, &without)
                .expect("missing-edge assessment");
            assert_eq!(refused.disposition(), AssessmentDisposition::Refused);
            assert!(refused.reasons().iter().any(|reason| {
                reason
                    == &format!(
                        "missing-prerequisite-receipt:{}:{}",
                        receipts[missing_index].prerequisite().id(),
                        match receipts[missing_index].use_kind() {
                            EvidenceUse::CalibrationInput => "calibration-input",
                            EvidenceUse::ValidationInput => "validation-input",
                        }
                    )
            }));
        }
    }

    let numerical = &positive[&EulerClaimKind::NumericalTrajectoryVerification];
    let calibrated = &positive[&EulerClaimKind::CalibratedReproduction];
    let blind = &positive[&EulerClaimKind::BlindTrajectoryPrediction];

    let calibrated_packet = packet(
        &contract,
        EulerClaimKind::CalibratedReproduction,
        records(&contract, EulerClaimKind::CalibratedReproduction),
        false,
        ReportedScientificDisposition::Positive,
        AssessmentDisposition::ReferenceCompleteCandidate,
    );
    let exact_numerical = numerical
        .as_prerequisite_for(
            EulerClaimKind::CalibratedReproduction,
            EvidenceUse::ValidationInput,
        )
        .expect("exact numerical edge");

    let wrong_use = numerical
        .as_prerequisite_for(
            EulerClaimKind::CalibratedReproduction,
            EvidenceUse::CalibrationInput,
        )
        .expect("wrong-use receipt remains a declaration");
    let refused = calibrated_packet
        .assess(&admitted, std::slice::from_ref(&wrong_use))
        .expect("wrong-use assessment");
    assert_eq!(refused.disposition(), AssessmentDisposition::Refused);
    assert!(
        refused.reasons().contains(
            &"unexpected-prerequisite-receipt:numerical-trajectory-verification:calibration-input"
                .to_owned()
        )
    );
    assert!(
        refused.reasons().contains(
            &"missing-prerequisite-receipt:numerical-trajectory-verification:validation-input"
                .to_owned()
        )
    );

    let wrong_dependent = numerical
        .as_prerequisite_for(EulerClaimKind::Ranking, EvidenceUse::ValidationInput)
        .expect("wrong-dependent receipt remains a declaration");
    let refused = calibrated_packet
        .assess(&admitted, std::slice::from_ref(&wrong_dependent))
        .expect("wrong-dependent assessment");
    assert_eq!(refused.disposition(), AssessmentDisposition::Refused);
    assert!(
        refused.reasons().contains(
            &"unexpected-prerequisite-receipt:numerical-trajectory-verification:validation-input"
                .to_owned()
        )
    );
    assert!(
        refused.reasons().contains(
            &"missing-prerequisite-receipt:numerical-trajectory-verification:validation-input"
                .to_owned()
        )
    );

    let invalid_forward = calibrated_packet
        .assess(&admitted, &[wrong_use.clone(), wrong_dependent.clone()])
        .expect("forward invalid receipt order");
    let invalid_reverse = calibrated_packet
        .assess(&admitted, &[wrong_dependent, wrong_use])
        .expect("reverse invalid receipt order");
    assert_eq!(
        invalid_forward, invalid_reverse,
        "invalid receipt permutations must retain the same first divergence, log, and assessment identity"
    );
    assert_eq!(
        invalid_forward.log().identity(),
        invalid_reverse.log().identity()
    );
    assert_eq!(invalid_forward.identity(), invalid_reverse.identity());

    let duplicate = calibrated_packet
        .assess(
            &admitted,
            &[exact_numerical.clone(), exact_numerical.clone()],
        )
        .expect("duplicate receipt assessment");
    assert_eq!(duplicate.disposition(), AssessmentDisposition::Refused);
    assert!(
        duplicate
            .reasons()
            .iter()
            .any(|reason| reason.starts_with("duplicate-prerequisite-receipt:"))
    );

    let maximum_receipts =
        vec![exact_numerical.clone(); fs_euler_disc_e2e::protocol::MAX_PREREQUISITE_RECEIPTS];
    assert_eq!(
        calibrated_packet
            .assess(&admitted, &maximum_receipts)
            .expect("the exact prerequisite cardinality boundary is admitted")
            .disposition(),
        AssessmentDisposition::Refused,
        "bounded duplicate receipts remain an ordinary refusal"
    );
    let maximum_plus_one =
        vec![exact_numerical.clone(); fs_euler_disc_e2e::protocol::MAX_PREREQUISITE_RECEIPTS + 1];
    let error = calibrated_packet
        .assess(&admitted, &maximum_plus_one)
        .expect_err("maximum-plus-one prerequisite receipts must refuse before accumulation");
    assert_eq!(error.code(), "EulerProtocolPrerequisiteCardinality");

    let event_packet = packet(
        &contract,
        EulerClaimKind::EventOrCrossoverPrediction,
        records(&contract, EulerClaimKind::EventOrCrossoverPrediction),
        false,
        ReportedScientificDisposition::Positive,
        AssessmentDisposition::ReferenceCompleteCandidate,
    );
    let indirect = numerical
        .as_prerequisite_for(
            EulerClaimKind::EventOrCrossoverPrediction,
            EvidenceUse::ValidationInput,
        )
        .expect("indirect receipt declaration");
    let refused = event_packet
        .assess(&admitted, &[indirect])
        .expect("indirect receipt assessment");
    assert_eq!(refused.disposition(), AssessmentDisposition::Refused);
    assert!(refused.reasons().iter().any(|reason| {
        reason == "missing-prerequisite-receipt:blind-trajectory-prediction:validation-input"
    }));

    let numerical_alt_point = packet_at_point(
        &contract,
        EulerClaimKind::NumericalTrajectoryVerification,
        point_at_fraction(&contract, 0.25),
        records(&contract, EulerClaimKind::NumericalTrajectoryVerification),
        false,
        ReportedScientificDisposition::Positive,
        AssessmentDisposition::ReferenceCompleteCandidate,
    )
    .assess(&admitted, &[])
    .expect("alternate-point numerical assessment");
    let point_mismatch = numerical_alt_point
        .as_prerequisite_for(
            EulerClaimKind::CalibratedReproduction,
            EvidenceUse::ValidationInput,
        )
        .expect("alternate-point receipt");
    let refused = calibrated_packet
        .assess(&admitted, &[point_mismatch])
        .expect("point-mismatch assessment");
    assert_eq!(refused.disposition(), AssessmentDisposition::Refused);
    assert!(
        refused
            .reasons()
            .iter()
            .any(|reason| { reason.starts_with("prerequisite-applicability-point-mismatch:") })
    );

    let blind_packet = packet(
        &contract,
        EulerClaimKind::BlindTrajectoryPrediction,
        records(&contract, EulerClaimKind::BlindTrajectoryPrediction),
        false,
        ReportedScientificDisposition::Positive,
        AssessmentDisposition::ReferenceCompleteCandidate,
    );
    let blind_receipts = [
        calibrated
            .as_prerequisite_for(
                EulerClaimKind::BlindTrajectoryPrediction,
                EvidenceUse::CalibrationInput,
            )
            .expect("calibration receipt"),
        numerical
            .as_prerequisite_for(
                EulerClaimKind::BlindTrajectoryPrediction,
                EvidenceUse::ValidationInput,
            )
            .expect("verification receipt"),
    ];
    let forward = blind_packet
        .assess(&admitted, &blind_receipts)
        .expect("forward receipt order");
    let reverse = blind_packet
        .assess(
            &admitted,
            &[blind_receipts[1].clone(), blind_receipts[0].clone()],
        )
        .expect("reverse receipt order");
    assert_eq!(forward, reverse, "receipt order must be canonicalized");

    for outcome in [
        ReportedScientificDisposition::Negative,
        ReportedScientificDisposition::Inconclusive,
    ] {
        let terminal = packet(
            &contract,
            EulerClaimKind::NumericalTrajectoryVerification,
            records(&contract, EulerClaimKind::NumericalTrajectoryVerification),
            false,
            outcome,
            AssessmentDisposition::RetainedTerminal,
        )
        .assess(&admitted, &[])
        .expect("reported terminal assessment");
        assert!(
            terminal
                .as_prerequisite_for(
                    EulerClaimKind::CalibratedReproduction,
                    EvidenceUse::ValidationInput,
                )
                .is_err(),
            "reported {outcome:?} cannot mint a prerequisite receipt"
        );
    }

    let requirement = golden_requirements(EulerClaimKind::NumericalTrajectoryVerification)[0];
    let demoted_records = golden_requirements(EulerClaimKind::NumericalTrajectoryVerification)
        .iter()
        .map(|candidate| {
            record_with(
                &contract,
                EulerClaimKind::NumericalTrajectoryVerification,
                *candidate,
                *candidate == requirement,
                hash(&format!("dag-demoted:{}", candidate.code())),
                declared_access_class(*candidate),
            )
        })
        .collect();
    let demoted = packet(
        &contract,
        EulerClaimKind::NumericalTrajectoryVerification,
        demoted_records,
        false,
        ReportedScientificDisposition::Positive,
        AssessmentDisposition::DemotedCandidate,
    )
    .assess(&admitted, &[])
    .expect("demoted assessment");
    assert_eq!(
        demoted.disposition(),
        AssessmentDisposition::DemotedCandidate
    );
    assert!(
        demoted
            .as_prerequisite_for(
                EulerClaimKind::CalibratedReproduction,
                EvidenceUse::ValidationInput,
            )
            .is_err()
    );

    let refused_source = packet(
        &contract,
        EulerClaimKind::NumericalTrajectoryVerification,
        Vec::new(),
        false,
        ReportedScientificDisposition::Positive,
        AssessmentDisposition::Refused,
    )
    .assess(&admitted, &[])
    .expect("refused source assessment");
    assert!(
        refused_source
            .as_prerequisite_for(
                EulerClaimKind::CalibratedReproduction,
                EvidenceUse::ValidationInput,
            )
            .is_err()
    );

    // The direct blind receipt itself verifies and binds the source assessment
    // identity; this ensures the happy edge is not only inferred by absence of
    // diagnostics.
    blind
        .as_prerequisite_for(
            EulerClaimKind::EventOrCrossoverPrediction,
            EvidenceUse::ValidationInput,
        )
        .expect("direct blind receipt")
        .verify()
        .expect("direct receipt identity");
}

#[test]
fn g0_declared_access_class_leakage_target_fitting_and_role_aliasing_refuse() {
    let contract = build_frozen_contract().expect("frozen contract");
    for claim in EulerClaimKind::ALL {
        for &requirement in golden_requirements(claim) {
            let expected = declared_access_class(requirement);
            let wrong_access = match expected {
                DeclaredEvidenceAccessClass::NotApplicable
                | DeclaredEvidenceAccessClass::BlindHoldout => {
                    DeclaredEvidenceAccessClass::Calibration
                }
                DeclaredEvidenceAccessClass::Calibration => DeclaredEvidenceAccessClass::Validation,
                DeclaredEvidenceAccessClass::Validation => {
                    DeclaredEvidenceAccessClass::BlindHoldout
                }
            };
            let wrong = records(&contract, claim)
                .into_iter()
                .map(|record| {
                    if record.requirement() == requirement {
                        record_with(
                            &contract,
                            claim,
                            requirement,
                            false,
                            record.artifact_hash(),
                            wrong_access,
                        )
                    } else {
                        record
                    }
                })
                .collect();
            let assessment = assess_with_frozen_prerequisites(
                &contract,
                packet(
                    &contract,
                    claim,
                    wrong,
                    false,
                    ReportedScientificDisposition::Positive,
                    AssessmentDisposition::Refused,
                ),
            );
            let expected_reason = format!(
                "access-class-mismatch:{}:expected-{}:observed-{}",
                requirement.code(),
                expected.code(),
                wrong_access.code()
            );
            assert_eq!(assessment.disposition(), AssessmentDisposition::Refused);
            assert!(
                assessment.reasons().contains(&expected_reason),
                "{claim:?} / {requirement:?} did not retain {expected_reason}: {:?}",
                assessment.reasons()
            );
        }
    }

    for kind in EulerClaimKind::ALL
        .into_iter()
        .filter(|kind| kind.forbids_target_fitting())
    {
        let assessment = assess_with_frozen_prerequisites(
            &contract,
            packet(
                &contract,
                kind,
                records(&contract, kind),
                true,
                ReportedScientificDisposition::Positive,
                AssessmentDisposition::Refused,
            ),
        );
        assert_eq!(assessment.disposition(), AssessmentDisposition::Refused);
        assert!(
            assessment
                .reasons()
                .contains(&"protected-target-fitting-invalidates-emergent-claim".to_owned())
        );
    }
    let calibrated = EulerClaimKind::CalibratedReproduction;
    assert_eq!(
        assess_with_frozen_prerequisites(
            &contract,
            packet(
                &contract,
                calibrated,
                records(&contract, calibrated),
                true,
                ReportedScientificDisposition::Positive,
                AssessmentDisposition::ReferenceCompleteCandidate,
            ),
        )
        .disposition(),
        AssessmentDisposition::ReferenceCompleteCandidate
    );

    let numerical = EulerClaimKind::NumericalTrajectoryVerification;
    let shared = hash("shared-cross-role-artifact");
    let aliased = vec![
        record_with(
            &contract,
            numerical,
            EvidenceRequirement::CodeVerification,
            false,
            shared,
            DeclaredEvidenceAccessClass::NotApplicable,
        ),
        record_with(
            &contract,
            numerical,
            EvidenceRequirement::SolutionVerification,
            false,
            shared,
            DeclaredEvidenceAccessClass::NotApplicable,
        ),
    ];
    let error = ClaimEvidencePacket::try_new(
        contract.identity(),
        "aliased-case",
        hash("aliased-case-design-set"),
        hash("aliased-case-aggregate-qoi-derivation-receipt"),
        numerical,
        point(&contract),
        aliased,
        true,
        false,
        ReportedScientificDisposition::Positive,
        AssessmentDisposition::Refused,
        vec!["1".to_owned()],
        ProtocolSeed::Fixed { value: 7 },
        ProtocolBudget::try_new(1, 1, 0.0).expect("budget"),
    )
    .expect_err("cross-role alias must refuse");
    assert_eq!(error.code(), "EulerProtocolCrossRoleEvidenceAlias");

    let complete = records(&contract, numerical);
    let code = complete
        .iter()
        .find(|record| record.requirement() == EvidenceRequirement::CodeVerification)
        .expect("code record");
    let solution = complete
        .iter()
        .find(|record| record.requirement() == EvidenceRequirement::SolutionVerification)
        .expect("solution record");
    for (label, case_id, design_set_identity, derivation_receipt_identity) in [
        (
            "design set aliasing an evidence artifact",
            "design-set-evidence-alias-case",
            code.artifact_hash(),
            hash("design-alias-distinct-aggregate-qoi-derivation-receipt"),
        ),
        (
            "aggregate-QoI derivation receipt aliasing an evidence schema receipt",
            "derivation-evidence-alias-case",
            hash("derivation-alias-distinct-design-set"),
            code.schema_admission_receipt_hash(),
        ),
    ] {
        let error = ClaimEvidencePacket::try_new(
            contract.identity(),
            case_id,
            design_set_identity,
            derivation_receipt_identity,
            numerical,
            point(&contract),
            complete.clone(),
            true,
            false,
            ReportedScientificDisposition::Positive,
            AssessmentDisposition::Refused,
            vec!["1".to_owned()],
            ProtocolSeed::Fixed { value: 7 },
            ProtocolBudget::try_new(1, 1, 0.0).expect("budget"),
        )
        .expect_err(label);
        assert_eq!(error.code(), "EulerProtocolCrossRoleEvidenceAlias");
    }
    let nested_alias = EvidenceRecord::try_new(
        code.contract_identity(),
        code.claim(),
        code.requirement(),
        code.qois().to_vec(),
        EvidenceAuthorityDeclaration::StructuralProcess {
            receipt_hash: solution.artifact_hash(),
        },
        code.artifact_hash(),
        code.source_id(),
        code.source_schema(),
        code.source_kind(),
        code.schema_admission_receipt_hash(),
        code.access_class(),
        code.independent(),
    )
    .expect("record-local nested hash remains nonzero");
    let error = ClaimEvidencePacket::try_new(
        contract.identity(),
        "nested-aliased-case",
        hash("nested-aliased-case-design-set"),
        hash("nested-aliased-case-aggregate-qoi-derivation-receipt"),
        numerical,
        point(&contract),
        vec![nested_alias, solution.clone()],
        true,
        false,
        ReportedScientificDisposition::Positive,
        AssessmentDisposition::Refused,
        vec!["1".to_owned()],
        ProtocolSeed::Fixed { value: 7 },
        ProtocolBudget::try_new(1, 1, 0.0).expect("budget"),
    )
    .expect_err("a nested receipt hash cannot alias another role's artifact");
    assert_eq!(error.code(), "EulerProtocolCrossRoleEvidenceAlias");
}

#[test]
fn g0_negative_and_inconclusive_are_retained_terminal_non_promotions() {
    let contract = build_frozen_contract().expect("frozen contract");
    for outcome in [
        ReportedScientificDisposition::Negative,
        ReportedScientificDisposition::Inconclusive,
    ] {
        for kind in EulerClaimKind::ALL {
            let assessment = assess_with_frozen_prerequisites(
                &contract,
                packet(
                    &contract,
                    kind,
                    records(&contract, kind),
                    false,
                    outcome,
                    AssessmentDisposition::RetainedTerminal,
                ),
            );
            assert_eq!(
                assessment.disposition(),
                AssessmentDisposition::RetainedTerminal
            );
            assert_eq!(assessment.reported_scientific_disposition(), outcome);
        }

        let replace_record = |claim: EulerClaimKind,
                              selected: EvidenceRequirement,
                              replacement_authority: Option<EvidenceAuthorityDeclaration>,
                              independent: bool| {
            records(&contract, claim)
                .into_iter()
                .map(|record| {
                    if record.requirement() != selected {
                        return record;
                    }
                    EvidenceRecord::try_new(
                        record.contract_identity(),
                        record.claim(),
                        record.requirement(),
                        record.qois().to_vec(),
                        replacement_authority
                            .clone()
                            .unwrap_or_else(|| record.authority().clone()),
                        record.artifact_hash(),
                        record.source_id(),
                        record.source_schema(),
                        record.source_kind(),
                        record.schema_admission_receipt_hash(),
                        record.access_class(),
                        independent,
                    )
                    .expect("locally well-formed weakness fixture")
                })
                .collect::<Vec<_>>()
        };
        let numerical = EulerClaimKind::NumericalTrajectoryVerification;
        let blind = EulerClaimKind::BlindTrajectoryPrediction;
        let estimated = |label: &str| Color::Estimated {
            estimator: format!("terminal-{label}-{}", outcome.code()),
            dispersion: 1.0,
        };
        let (excluded_axis, excluded_bound) = contract
            .context()
            .applicability()
            .numeric()
            .iter()
            .find_map(|(axis, domain)| {
                let (lo, hi) = domain.bounds();
                (lo < hi).then_some((axis.as_str(), lo))
            })
            .expect("frozen context has a nondegenerate numeric applicability axis");
        let noncovering_regime =
            covering_regime(&contract).with(excluded_axis, excluded_bound, excluded_bound);
        let weakness_fixtures = vec![
            (
                "authority-class",
                numerical,
                replace_record(
                    numerical,
                    EvidenceRequirement::CodeVerification,
                    Some(EvidenceAuthorityDeclaration::VerifiedNumerics {
                        color: estimated("authority-class"),
                    }),
                    true,
                ),
                "weak-authority:code-verification:requires-structural-process:observed-verified-numerics"
                    .to_owned(),
            ),
            (
                "verified-color-kind",
                numerical,
                replace_record(
                    numerical,
                    EvidenceRequirement::SolutionVerification,
                    Some(EvidenceAuthorityDeclaration::VerifiedNumerics {
                        color: estimated("verified-color-kind"),
                    }),
                    true,
                ),
                "weak-authority:solution-verification:requires-finite-verified-color".to_owned(),
            ),
            (
                "vacuous-verified-enclosure",
                numerical,
                replace_record(
                    numerical,
                    EvidenceRequirement::SolutionVerification,
                    Some(EvidenceAuthorityDeclaration::VerifiedNumerics {
                        color: Color::Verified {
                            lo: f64::NEG_INFINITY,
                            hi: f64::INFINITY,
                        },
                    }),
                    true,
                ),
                "weak-authority:solution-verification:verified-enclosure-is-vacuous".to_owned(),
            ),
            (
                "validated-color-kind",
                blind,
                replace_record(
                    blind,
                    EvidenceRequirement::PhysicalValidation,
                    Some(EvidenceAuthorityDeclaration::ValidatedPhysical {
                        color: estimated("validated-color-kind"),
                    }),
                    true,
                ),
                "weak-authority:physical-validation:requires-validated-color".to_owned(),
            ),
            (
                "validity-domain",
                blind,
                replace_record(
                    blind,
                    EvidenceRequirement::PhysicalValidation,
                    Some(EvidenceAuthorityDeclaration::ValidatedPhysical {
                        color: Color::Validated {
                            regime: noncovering_regime,
                            dataset: format!("terminal-validity-domain-{}", outcome.code()),
                        },
                    }),
                    true,
                ),
                "weak-validity-domain:physical-validation:does-not-cover-case".to_owned(),
            ),
            (
                "independence",
                numerical,
                replace_record(
                    numerical,
                    EvidenceRequirement::IndependentReconstruction,
                    None,
                    false,
                ),
                "weak-independence:independent-reconstruction:independent-evidence-required"
                    .to_owned(),
            ),
        ];

        for (weakness, claim, evidence, expected_reason) in weakness_fixtures {
            let assessment = assess_with_frozen_prerequisites(
                &contract,
                packet(
                    &contract,
                    claim,
                    evidence,
                    false,
                    outcome,
                    AssessmentDisposition::RetainedTerminal,
                ),
            );
            assert_eq!(
                assessment.disposition(),
                AssessmentDisposition::RetainedTerminal,
                "{outcome:?} with {weakness} weakness must remain terminal: {:?}",
                assessment.reasons()
            );
            assert_eq!(assessment.reasons(), &[expected_reason]);
            fs_euler_disc_e2e::ClaimPolicyAssessmentLog::from_json_line(
                assessment.log().json_line(),
            )
            .unwrap_or_else(|error| {
                panic!("{outcome:?} with {weakness} weakness must round-trip: {error}")
            });

            let forged_demotion = assessment
                .log()
                .json_line()
                .replacen(
                    "\"expected_disposition\":\"retained-terminal-non-promotion\",\"observed_disposition\":\"retained-terminal-non-promotion\"",
                    "\"expected_disposition\":\"demoted-candidate\",\"observed_disposition\":\"demoted-candidate\"",
                    1,
                )
                .replacen(
                    "\"authority_state\":\"terminal-non-promotion\"",
                    "\"authority_state\":\"demoted-below-requested-claim\"",
                    1,
                );
            let error = fs_euler_disc_e2e::ClaimPolicyAssessmentLog::from_json_line(
                forged_demotion,
            )
            .expect_err(
                "a nonpositive terminal outcome cannot be relabelled as a demoted candidate",
            );
            assert_eq!(error.code(), "EulerProtocolMalformedAssessmentLog");
        }
    }
}

#[test]
fn g0_assessment_log_binds_packet_contract_identity_bidirectionally() {
    let contract = build_frozen_contract().expect("frozen contract");
    let admitted = admit_frozen_contract(contract.clone()).expect("structural admission");
    let claim = EulerClaimKind::NumericalTrajectoryVerification;
    let foreign_contract = ContractIdentity::from_hash(hash("foreign-packet-contract"));
    let foreign_records = records(&contract, claim)
        .into_iter()
        .map(|record| {
            EvidenceRecord::try_new(
                foreign_contract,
                record.claim(),
                record.requirement(),
                record.qois().to_vec(),
                record.authority().clone(),
                record.artifact_hash(),
                record.source_id(),
                record.source_schema(),
                record.source_kind(),
                record.schema_admission_receipt_hash(),
                record.access_class(),
                record.independent(),
            )
            .expect("foreign packet evidence remains internally contract-bound")
        })
        .collect();
    let foreign_packet = packet_from_parts(
        foreign_contract,
        "foreign-contract-packet",
        claim,
        point(&contract),
        foreign_records,
        true,
        false,
        ReportedScientificDisposition::Positive,
        AssessmentDisposition::Refused,
        claim_units(&contract, claim),
        ProtocolSeed::not_applicable("deterministic-fixture").expect("seed"),
        ProtocolBudget::try_new(60_000, 64 * 1024 * 1024, 0.05).expect("budget"),
    );
    let assessment = foreign_packet
        .assess(&admitted, &[])
        .expect("a stale packet contract must remain a retained refusal");
    assert_eq!(assessment.disposition(), AssessmentDisposition::Refused);
    assert!(
        assessment
            .reasons()
            .binary_search(&"contract-identity-mismatch".to_owned())
            .is_ok()
    );
    let line = assessment.log().json_line();
    let frozen_hex = contract.identity().as_hash().to_hex();
    let foreign_hex = foreign_contract.as_hash().to_hex();
    assert!(line.contains(&format!("\"contract_identity\":\"{frozen_hex}\"")));
    assert!(line.contains(&format!("\"packet_contract_identity\":\"{foreign_hex}\"")));
    ClaimPolicyAssessmentLog::from_json_line(line)
        .expect("the truthful packet-contract mismatch must round-trip");

    let alternate_reason = assessment
        .reasons()
        .iter()
        .find(|reason| reason.as_str() != "contract-identity-mismatch")
        .expect("foreign evidence rows also retain a stale binding reason");
    let moved_first_divergence = line.replacen(
        "\"first_divergence\":\"contract-identity-mismatch\"",
        &format!("\"first_divergence\":\"{alternate_reason}\""),
        1,
    );
    assert_ne!(
        moved_first_divergence, line,
        "the hostile fixture must move the exact first-divergence field"
    );
    let reasons_prefix = "\"reasons\":[\"contract-identity-mismatch\",";
    assert!(
        moved_first_divergence.contains(reasons_prefix),
        "the hostile fixture requires the exact canonically first reason row"
    );
    let missing_reason = moved_first_divergence.replacen(reasons_prefix, "\"reasons\":[", 1);
    ClaimPolicyAssessmentLog::from_json_line(missing_reason)
        .expect_err("a mismatched packet contract cannot omit its exact reason");

    let false_reason = line.replacen(
        &format!("\"packet_contract_identity\":\"{foreign_hex}\""),
        &format!("\"packet_contract_identity\":\"{frozen_hex}\""),
        1,
    );
    ClaimPolicyAssessmentLog::from_json_line(false_reason)
        .expect_err("a matching packet contract cannot retain a false mismatch reason");

    let matching_refusal = assess_with_frozen_prerequisites(
        &contract,
        packet(
            &contract,
            EulerClaimKind::BlindTrajectoryPrediction,
            records(&contract, EulerClaimKind::BlindTrajectoryPrediction),
            true,
            ReportedScientificDisposition::Positive,
            AssessmentDisposition::Refused,
        ),
    );
    assert!(
        matching_refusal
            .reasons()
            .iter()
            .any(|reason| { reason == "protected-target-fitting-invalidates-emergent-claim" })
    );
    let unreported_mismatch = matching_refusal.log().json_line().replacen(
        &format!("\"packet_contract_identity\":\"{frozen_hex}\""),
        &format!("\"packet_contract_identity\":\"{foreign_hex}\""),
        1,
    );
    ClaimPolicyAssessmentLog::from_json_line(unreported_mismatch)
        .expect_err("a forged packet-contract mismatch requires its exact retained reason");
}

#[test]
fn g0_exact_hypothesis_hash_and_cross_role_substitutions_refuse() {
    let contract = build_frozen_contract().expect("frozen contract");
    let numerical = EulerClaimKind::NumericalTrajectoryVerification;
    let source_hash = contract.extension().hypothesis_sources()[0].declaration_hash();
    let substituted = records(&contract, numerical)
        .into_iter()
        .map(|record| {
            if record.requirement() == EvidenceRequirement::CodeVerification {
                record_with(
                    &contract,
                    numerical,
                    EvidenceRequirement::CodeVerification,
                    false,
                    source_hash,
                    DeclaredEvidenceAccessClass::NotApplicable,
                )
            } else {
                record
            }
        })
        .collect();
    let assessment = assess_with_frozen_prerequisites(
        &contract,
        packet(
            &contract,
            numerical,
            substituted,
            false,
            ReportedScientificDisposition::Positive,
            AssessmentDisposition::Refused,
        ),
    );
    assert_eq!(assessment.disposition(), AssessmentDisposition::Refused);
    assert!(
        assessment
            .reasons()
            .iter()
            .any(|reason| reason.starts_with("hypothesis-source-cannot-satisfy-evidence:"))
    );

    for nested_slot in ["role-receipt", "schema-admission-receipt"] {
        let substituted = records(&contract, numerical)
            .into_iter()
            .map(|record| {
                if record.requirement() != EvidenceRequirement::CodeVerification {
                    return record;
                }
                let authority = if nested_slot == "role-receipt" {
                    EvidenceAuthorityDeclaration::StructuralProcess {
                        receipt_hash: source_hash,
                    }
                } else {
                    record.authority().clone()
                };
                let schema_receipt = if nested_slot == "schema-admission-receipt" {
                    source_hash
                } else {
                    record.schema_admission_receipt_hash()
                };
                EvidenceRecord::try_new(
                    record.contract_identity(),
                    record.claim(),
                    record.requirement(),
                    record.qois().to_vec(),
                    authority,
                    record.artifact_hash(),
                    record.source_id(),
                    record.source_schema(),
                    record.source_kind(),
                    schema_receipt,
                    record.access_class(),
                    record.independent(),
                )
                .expect("nested hypothesis hash is locally nonzero")
            })
            .collect();
        let assessment = assess_with_frozen_prerequisites(
            &contract,
            packet(
                &contract,
                numerical,
                substituted,
                false,
                ReportedScientificDisposition::Positive,
                AssessmentDisposition::Refused,
            ),
        );
        assert_eq!(assessment.disposition(), AssessmentDisposition::Refused);
        assert!(assessment.reasons().iter().any(|reason| {
            reason
                == &format!(
                    "hypothesis-source-cannot-satisfy-evidence:code-verification:{nested_slot}"
                )
        }));
    }

    // Packet construction can enforce local hash-role uniqueness, but only a
    // contract-bound assessment knows the frozen hypothesis registry. Exercise
    // both opaque top-level roles at that semantic boundary and require the
    // strict retained-log reader to reproduce the exact refusal.
    for (role, design_set_identity, aggregate_derivation_identity) in [
        (
            "design-set",
            source_hash,
            hash("hypothesis-collision-distinct-aggregate-derivation"),
        ),
        (
            "aggregate-qoi-derivation-receipt",
            hash("hypothesis-collision-distinct-design-set"),
            source_hash,
        ),
    ] {
        let collision = ClaimEvidencePacket::try_new(
            contract.identity(),
            format!("hypothesis-collision-{role}"),
            design_set_identity,
            aggregate_derivation_identity,
            numerical,
            point(&contract),
            records(&contract, numerical),
            true,
            false,
            ReportedScientificDisposition::Positive,
            AssessmentDisposition::Refused,
            claim_units(&contract, numerical),
            ProtocolSeed::not_applicable("deterministic-fixture").expect("seed"),
            ProtocolBudget::try_new(60_000, 64 * 1024 * 1024, 0.05).expect("budget"),
        )
        .expect("opaque packet roles remain locally distinct and well formed");
        let assessment = assess_with_frozen_prerequisites(&contract, collision);
        let expected_reason = format!("hypothesis-source-cannot-satisfy-packet-role:{role}");
        assert_eq!(assessment.disposition(), AssessmentDisposition::Refused);
        assert_eq!(assessment.reasons(), &[expected_reason]);
        fs_euler_disc_e2e::ClaimPolicyAssessmentLog::from_json_line(assessment.log().json_line())
            .unwrap_or_else(|error| panic!("top-level collision {role} must round-trip: {error}"));
    }

    let blind = EulerClaimKind::BlindTrajectoryPrediction;
    let software_as_physical = records(&contract, blind)
        .into_iter()
        .map(|record| {
            if record.requirement() == EvidenceRequirement::PhysicalValidation {
                EvidenceRecord::try_new(
                    contract.identity(),
                    blind,
                    EvidenceRequirement::PhysicalValidation,
                    record.qois().to_vec(),
                    EvidenceAuthorityDeclaration::VerifiedNumerics {
                        color: Color::Verified { lo: 0.0, hi: 0.0 },
                    },
                    record.artifact_hash(),
                    record.source_id(),
                    EvidenceRequirement::PhysicalValidation.source_schema(),
                    EvidenceRequirement::PhysicalValidation.source_kind(),
                    record.schema_admission_receipt_hash(),
                    DeclaredEvidenceAccessClass::Validation,
                    true,
                )
                .expect("structurally valid but scientifically weak row")
            } else {
                record
            }
        })
        .collect();
    assert_eq!(
        assess_with_frozen_prerequisites(
            &contract,
            packet(
                &contract,
                blind,
                software_as_physical,
                false,
                ReportedScientificDisposition::Positive,
                AssessmentDisposition::DemotedCandidate,
            ),
        )
        .disposition(),
        AssessmentDisposition::DemotedCandidate
    );

    let mechanism = EulerClaimKind::MechanismAttribution;
    let without_rival = records(&contract, mechanism)
        .into_iter()
        .filter(|record| record.requirement() != EvidenceRequirement::RivalMechanismDiscrimination)
        .collect();
    assert_eq!(
        assess_with_frozen_prerequisites(
            &contract,
            packet(
                &contract,
                mechanism,
                without_rival,
                false,
                ReportedScientificDisposition::Positive,
                AssessmentDisposition::Refused,
            ),
        )
        .disposition(),
        AssessmentDisposition::Refused
    );
}

fn context_with(
    base: &EulerScientificContract,
    decision: String,
    qois: Vec<QoiSpec>,
    applicability: ApplicabilityDomain,
) -> ContextOfUse {
    ContextOfUse::try_new(
        base.context().header().clone(),
        decision,
        qois,
        applicability,
        base.context().applicability_policy(),
    )
    .expect("mutated context")
}

fn contract_with_context(
    base: &EulerScientificContract,
    context: ContextOfUse,
) -> EulerScientificContract {
    EulerScientificContract::try_new(
        context,
        base.extension().clone(),
        base.claim_graph().clone(),
        base.no_claims().clone(),
        base.owner_matrix().clone(),
    )
    .expect("mutated contract")
}

fn extension_with(
    base: &EulerScientificContract,
    apparatus: String,
    environment: String,
    frame: String,
) -> EulerContextExtension {
    EulerContextExtension::try_new(
        base.extension().users().to_vec(),
        apparatus,
        environment,
        frame,
        base.extension().decision_alternatives().to_vec(),
        base.extension().risks().to_vec(),
        base.extension().hypothesis_sources().to_vec(),
    )
    .expect("extension")
}

fn contract_with_extension(
    base: &EulerScientificContract,
    extension: EulerContextExtension,
) -> Result<EulerScientificContract, fs_euler_disc_e2e::ContractError> {
    EulerScientificContract::try_new(
        base.context().clone(),
        extension,
        base.claim_graph().clone(),
        base.no_claims().clone(),
        base.owner_matrix().clone(),
    )
}

#[test]
fn g3_every_required_context_semantic_moves_identity_and_stales_receipts() {
    let base = build_frozen_contract().expect("frozen contract");
    let receipt = check_frozen_contract(&base).expect("receipt");
    assert!(receipt.passed());
    let base_qois = base.context().qois().values().cloned().collect::<Vec<_>>();
    let mut mutations = Vec::new();

    mutations.push(contract_with_context(
        &base,
        context_with(
            &base,
            format!("{} {}", base.context().decision(), "Changed decision."),
            base_qois.clone(),
            base.context().applicability().clone(),
        ),
    ));

    let renamed_qois = base_qois
        .iter()
        .map(|qoi| {
            if qoi.id().as_str() == "numerical-trajectory-error" {
                QoiSpec::try_new(
                    qoi.id().clone(),
                    "Renamed independent numerical trajectory discrepancy",
                    qoi.unit().clone(),
                    qoi.acceptance().clone(),
                )
                .expect("renamed QoI")
            } else {
                qoi.clone()
            }
        })
        .collect();
    mutations.push(contract_with_context(
        &base,
        context_with(
            &base,
            base.context().decision().to_owned(),
            renamed_qois,
            base.context().applicability().clone(),
        ),
    ));

    let threshold_qois = base_qois
        .iter()
        .map(|qoi| {
            if qoi.id().as_str() == "numerical-trajectory-error" {
                QoiSpec::try_new(
                    qoi.id().clone(),
                    qoi.name(),
                    qoi.unit().clone(),
                    AcceptanceCriterion::ClosedRange {
                        lo: 0.0,
                        hi: 2.0e-8,
                    },
                )
                .expect("changed threshold")
            } else {
                qoi.clone()
            }
        })
        .collect();
    mutations.push(contract_with_context(
        &base,
        context_with(
            &base,
            base.context().decision().to_owned(),
            threshold_qois,
            base.context().applicability().clone(),
        ),
    ));

    let numeric = base
        .context()
        .applicability()
        .numeric()
        .values()
        .map(|axis| {
            let (lo, hi) = axis.bounds();
            NumericDomainAxis::try_new(
                axis.axis().clone(),
                axis.unit().clone(),
                lo,
                if axis.axis().as_str() == "outer-radius" {
                    hi + 1.0
                } else {
                    hi
                },
            )
            .expect("changed domain")
        })
        .collect();
    let categorical = base
        .context()
        .applicability()
        .categorical()
        .values()
        .cloned()
        .collect::<Vec<CategoricalDomainAxis>>();
    mutations.push(contract_with_context(
        &base,
        context_with(
            &base,
            base.context().decision().to_owned(),
            base_qois,
            ApplicabilityDomain::try_new(numeric, categorical).expect("changed applicability"),
        ),
    ));

    mutations.push(
        contract_with_extension(
            &base,
            extension_with(
                &base,
                format!(
                    "{} {}",
                    base.extension().apparatus_population(),
                    "Changed apparatus population."
                ),
                base.extension().environment_population().to_owned(),
                base.extension().observation_frame().to_owned(),
            ),
        )
        .expect("apparatus mutation"),
    );
    mutations.push(
        contract_with_extension(
            &base,
            extension_with(
                &base,
                base.extension().apparatus_population().to_owned(),
                format!(
                    "{} {}",
                    base.extension().environment_population(),
                    "Changed environment population."
                ),
                base.extension().observation_frame().to_owned(),
            ),
        )
        .expect("environment mutation"),
    );

    let mut users = base.extension().users().to_vec();
    users.push("Independent scientific auditor".to_owned());
    mutations.push(
        contract_with_extension(
            &base,
            EulerContextExtension::try_new(
                users,
                base.extension().apparatus_population(),
                base.extension().environment_population(),
                base.extension().observation_frame(),
                base.extension().decision_alternatives().to_vec(),
                base.extension().risks().to_vec(),
                base.extension().hypothesis_sources().to_vec(),
            )
            .expect("user mutation"),
        )
        .expect("user contract mutation"),
    );

    let mut alternatives = base.extension().decision_alternatives().to_vec();
    alternatives.push("Defer mechanism attribution pending rival tests".to_owned());
    mutations.push(
        contract_with_extension(
            &base,
            EulerContextExtension::try_new(
                base.extension().users().to_vec(),
                base.extension().apparatus_population(),
                base.extension().environment_population(),
                base.extension().observation_frame(),
                alternatives,
                base.extension().risks().to_vec(),
                base.extension().hypothesis_sources().to_vec(),
            )
            .expect("decision-alternative mutation"),
        )
        .expect("decision-alternative contract mutation"),
    );

    let mut risks = base.extension().risks().to_vec();
    let first_risk = risks.first().expect("frozen risk registry").clone();
    risks[0] = ScientificRisk::try_new(
        first_risk.code(),
        format!("{} {}", first_risk.consequence(), "Changed consequence."),
        first_risk.severity(),
        first_risk.affected_claims().to_vec(),
        first_risk.decision_alternative(),
    )
    .expect("risk mutation");
    mutations.push(
        contract_with_extension(
            &base,
            EulerContextExtension::try_new(
                base.extension().users().to_vec(),
                base.extension().apparatus_population(),
                base.extension().environment_population(),
                base.extension().observation_frame(),
                base.extension().decision_alternatives().to_vec(),
                risks,
                base.extension().hypothesis_sources().to_vec(),
            )
            .expect("risk-registry mutation"),
        )
        .expect("risk-registry contract mutation"),
    );

    let mut sources = base.extension().hypothesis_sources().to_vec();
    let first_source = sources.first().expect("frozen hypothesis sources").clone();
    sources[0] = HypothesisSource::try_new(
        first_source.id(),
        format!("{}#identity-mutation", first_source.locator()),
    )
    .expect("hypothesis-source mutation");
    mutations.push(
        contract_with_extension(
            &base,
            EulerContextExtension::try_new(
                base.extension().users().to_vec(),
                base.extension().apparatus_population(),
                base.extension().environment_population(),
                base.extension().observation_frame(),
                base.extension().decision_alternatives().to_vec(),
                base.extension().risks().to_vec(),
                sources,
            )
            .expect("hypothesis-source registry mutation"),
        )
        .expect("hypothesis-source contract mutation"),
    );

    let mut no_claims = base.no_claims().entries().to_vec();
    no_claims.push("Additional conservative no-claim for mutation proof.".to_owned());
    let refs = no_claims.iter().map(String::as_str).collect::<Vec<_>>();
    mutations.push(
        EulerScientificContract::try_new(
            base.context().clone(),
            base.extension().clone(),
            base.claim_graph().clone(),
            NoClaimBoundary::new(&refs).expect("no claims"),
            base.owner_matrix().clone(),
        )
        .expect("no-claim mutation"),
    );

    for mutation in mutations {
        assert_ne!(mutation.identity(), base.identity());
        assert!(receipt.verify_subject(&mutation).is_err());
        let check = check_frozen_contract(&mutation).expect("mutated check");
        assert!(!check.passed());
        assert!(
            check
                .issues()
                .contains(&"not-the-literal-frozen-v1-contract".to_owned())
        );
        assert!(check.issues().iter().all(|issue| {
            issue == "not-the-literal-frozen-v1-contract"
                || issue == "not-the-literal-frozen-v1-context"
        }));
    }
}

#[test]
fn g0_unit_and_frame_mismatches_refuse_before_identity_publication() {
    let base = build_frozen_contract().expect("frozen contract");
    let qois = base
        .context()
        .qois()
        .values()
        .map(|qoi| {
            if qoi.id().as_str() == "event-time-error" {
                QoiSpec::try_new(
                    qoi.id().clone(),
                    qoi.name(),
                    UnitId::try_new("ms").expect("unit"),
                    qoi.acceptance().clone(),
                )
                .expect("unit-mutated qoi")
            } else {
                qoi.clone()
            }
        })
        .collect();
    let context = context_with(
        &base,
        base.context().decision().to_owned(),
        qois,
        base.context().applicability().clone(),
    );
    let error = EulerScientificContract::try_new(
        context,
        base.extension().clone(),
        base.claim_graph().clone(),
        base.no_claims().clone(),
        base.owner_matrix().clone(),
    )
    .expect_err("undeclared QoI unit");
    assert_eq!(error.code(), "EulerContractUnitCoverage");

    let numeric = base
        .context()
        .applicability()
        .numeric()
        .values()
        .map(|axis| {
            let (lo, hi) = axis.bounds();
            let unit = if axis.axis().as_str() == "outer-radius" {
                UnitId::try_new("cm").expect("undeclared but valid axis unit")
            } else {
                axis.unit().clone()
            };
            NumericDomainAxis::try_new(axis.axis().clone(), unit, lo, hi)
                .expect("unit-mutated numeric applicability axis")
        })
        .collect();
    let categorical = base
        .context()
        .applicability()
        .categorical()
        .values()
        .cloned()
        .collect();
    let applicability =
        ApplicabilityDomain::try_new(numeric, categorical).expect("unit-mutated applicability");
    let context = context_with(
        &base,
        base.context().decision().to_owned(),
        base.context().qois().values().cloned().collect(),
        applicability,
    );
    let error = EulerScientificContract::try_new(
        context,
        base.extension().clone(),
        base.claim_graph().clone(),
        base.no_claims().clone(),
        base.owner_matrix().clone(),
    )
    .expect_err("undeclared numeric applicability-axis unit");
    assert_eq!(error.code(), "EulerContractUnitCoverage");

    let wrong_frame = extension_with(
        &base,
        base.extension().apparatus_population().to_owned(),
        base.extension().environment_population().to_owned(),
        "body-fixed-frame".to_owned(),
    );
    let error = contract_with_extension(&base, wrong_frame).expect_err("frame mismatch");
    assert_eq!(error.code(), "EulerContractFrameMismatch");
}

#[test]
fn assessment_log_unit_set_binding_is_delimiter_unambiguous() {
    let contract = build_frozen_contract().expect("frozen contract");
    let claim = EulerClaimKind::EnergyChannelAttribution;
    assert_eq!(
        claim_units(&contract, claim),
        vec!["1".to_owned(), "j".to_owned()]
    );

    let wrong_units = assess_with_frozen_prerequisites(
        &contract,
        packet_with_policy_declarations(
            &contract,
            claim,
            point(&contract),
            records(&contract, claim),
            true,
            false,
            ReportedScientificDisposition::Positive,
            AssessmentDisposition::Refused,
            vec!["1+j".to_owned()],
        ),
    );
    assert_eq!(wrong_units.disposition(), AssessmentDisposition::Refused);
    assert!(
        wrong_units
            .reasons()
            .contains(&"claim-unit-set-mismatch:expected-1,j:observed-1+j".to_owned())
    );
    fs_euler_disc_e2e::ClaimPolicyAssessmentLog::from_json_line(wrong_units.log().json_line())
        .expect("the unambiguous comma-framed mismatch reason must round-trip");

    let positive = assess_with_frozen_prerequisites(
        &contract,
        packet_with_policy_declarations(
            &contract,
            claim,
            point(&contract),
            records(&contract, claim),
            true,
            false,
            ReportedScientificDisposition::Positive,
            AssessmentDisposition::ReferenceCompleteCandidate,
            claim_units(&contract, claim),
        ),
    );
    assert_eq!(
        positive.disposition(),
        AssessmentDisposition::ReferenceCompleteCandidate
    );
    let forged =
        positive
            .log()
            .json_line()
            .replacen("\"units\":[\"1\",\"j\"]", "\"units\":[\"1+j\"]", 1);
    assert_ne!(forged, positive.log().json_line());
    let error = fs_euler_disc_e2e::ClaimPolicyAssessmentLog::from_json_line(forged)
        .expect_err("a delimiter-collision unit vector must not decode as a positive log");
    assert_eq!(error.code(), "EulerProtocolMalformedAssessmentLog");
}

#[test]
fn g0_policy_declaration_guards_are_retained_and_fail_closed() {
    let contract = build_frozen_contract().expect("frozen contract");
    let admitted = admit_frozen_contract(contract.clone()).expect("structural admission");
    let claim = EulerClaimKind::NumericalTrajectoryVerification;

    let no_claim_opt_out = packet_with_policy_declarations(
        &contract,
        claim,
        point(&contract),
        records(&contract, claim),
        false,
        false,
        ReportedScientificDisposition::Positive,
        AssessmentDisposition::Refused,
        claim_units(&contract, claim),
    )
    .assess(&admitted, &[])
    .expect("no-claim opt-out assessment");
    assert_eq!(
        no_claim_opt_out.disposition(),
        AssessmentDisposition::Refused
    );
    assert!(
        no_claim_opt_out
            .reasons()
            .contains(&"binding-no-claims-not-accepted".to_owned())
    );

    let outside = packet_with_policy_declarations(
        &contract,
        claim,
        point_at_fraction(&contract, 1.25),
        records(&contract, claim),
        true,
        false,
        ReportedScientificDisposition::Positive,
        AssessmentDisposition::Refused,
        claim_units(&contract, claim),
    )
    .assess(&admitted, &[])
    .expect("out-of-domain assessment");
    assert_eq!(outside.disposition(), AssessmentDisposition::Refused);
    assert!(
        outside
            .reasons()
            .iter()
            .any(|reason| reason.starts_with("out-of-domain-numeric:"))
    );

    let wrong_units = packet_with_policy_declarations(
        &contract,
        claim,
        point(&contract),
        records(&contract, claim),
        true,
        false,
        ReportedScientificDisposition::Positive,
        AssessmentDisposition::Refused,
        vec!["deliberately-wrong-unit".to_owned()],
    )
    .assess(&admitted, &[])
    .expect("wrong-unit assessment");
    assert_eq!(wrong_units.disposition(), AssessmentDisposition::Refused);
    assert!(
        wrong_units
            .reasons()
            .iter()
            .any(|reason| reason.starts_with("claim-unit-set-mismatch:"))
    );

    let maximal_units = (0..18)
        .map(|index| {
            let prefix = format!("u{index:02}");
            format!("{prefix}{}", "x".repeat(256 - prefix.len()))
        })
        .collect::<Vec<_>>();
    assert!(
        maximal_units
            .iter()
            .all(|unit| unit.len() == 256 && unit.len() <= MAX_EULER_TEXT_BYTES)
    );
    let maximal_unit_mismatch = packet_with_policy_declarations(
        &contract,
        claim,
        point(&contract),
        records(&contract, claim),
        true,
        false,
        ReportedScientificDisposition::Positive,
        AssessmentDisposition::Refused,
        maximal_units,
    )
    .assess(&admitted, &[])
    .expect("maximum canonical unit-set mismatch remains a policy refusal");
    assert_eq!(
        maximal_unit_mismatch.disposition(),
        AssessmentDisposition::Refused
    );
    assert_eq!(maximal_unit_mismatch.reasons().len(), 1);
    let long_reason = &maximal_unit_mismatch.reasons()[0];
    assert!(long_reason.starts_with("claim-unit-set-mismatch:"));
    assert!(
        long_reason.len() > MAX_EULER_TEXT_BYTES,
        "the regression must cross the generic single-text bound"
    );
    assert!(
        maximal_unit_mismatch
            .log()
            .json_line()
            .contains(&format!("\"first_divergence\":\"{long_reason}\""))
    );
    let decoded = fs_euler_disc_e2e::ClaimPolicyAssessmentLog::from_json_line(
        maximal_unit_mismatch.log().json_line(),
    )
    .expect("strict reader admits its maximum-unit generated refusal");
    assert_eq!(decoded, maximal_unit_mismatch.log().clone());

    let expectation_mismatch = packet_with_policy_declarations(
        &contract,
        claim,
        point(&contract),
        records(&contract, claim),
        true,
        false,
        ReportedScientificDisposition::Positive,
        AssessmentDisposition::Refused,
        claim_units(&contract, claim),
    )
    .assess(&admitted, &[])
    .expect("expectation-mismatch assessment");
    assert_eq!(
        expectation_mismatch.disposition(),
        AssessmentDisposition::ReferenceCompleteCandidate
    );
    assert!(expectation_mismatch.reasons().contains(
        &"expected-disposition-mismatch:expected-refused:observed-reference-complete-candidate-unreadmitted"
            .to_owned()
    ));

    let weak_independence = records(&contract, claim)
        .into_iter()
        .map(|record| {
            if record.requirement() != EvidenceRequirement::IndependentReconstruction {
                return record;
            }
            EvidenceRecord::try_new(
                record.contract_identity(),
                record.claim(),
                record.requirement(),
                record.qois().to_vec(),
                record.authority().clone(),
                record.artifact_hash(),
                record.source_id(),
                record.source_schema(),
                record.source_kind(),
                record.schema_admission_receipt_hash(),
                record.access_class(),
                false,
            )
            .expect("structurally formed non-independent declaration")
        })
        .collect();
    let weak_independence = packet(
        &contract,
        claim,
        weak_independence,
        false,
        ReportedScientificDisposition::Positive,
        AssessmentDisposition::DemotedCandidate,
    )
    .assess(&admitted, &[])
    .expect("weak-independence assessment");
    assert_eq!(
        weak_independence.disposition(),
        AssessmentDisposition::DemotedCandidate
    );
    assert!(weak_independence.reasons().contains(
        &"weak-independence:independent-reconstruction:independent-evidence-required".to_owned()
    ));
}

#[test]
fn g0_maximum_source_id_assessment_log_round_trips_and_plus_one_refuses() {
    let contract = build_frozen_contract().expect("frozen contract");
    let admitted = admit_frozen_contract(contract.clone()).expect("structural admission");
    let claim = EulerClaimKind::NumericalTrajectoryVerification;
    let selected = EvidenceRequirement::CodeVerification;
    let maximal_source_id = "s".repeat(MAX_PROTOCOL_ID_BYTES);
    let evidence = records(&contract, claim)
        .into_iter()
        .map(|record| {
            if record.requirement() != selected {
                return record;
            }
            EvidenceRecord::try_new(
                record.contract_identity(),
                record.claim(),
                record.requirement(),
                record.qois().to_vec(),
                record.authority().clone(),
                record.artifact_hash(),
                maximal_source_id.clone(),
                record.source_schema(),
                record.source_kind(),
                record.schema_admission_receipt_hash(),
                record.access_class(),
                record.independent(),
            )
            .expect("the exact source-id byte maximum must be admitted")
        })
        .collect::<Vec<_>>();
    let assessment = packet(
        &contract,
        claim,
        evidence,
        false,
        ReportedScientificDisposition::Positive,
        AssessmentDisposition::ReferenceCompleteCandidate,
    )
    .assess(&admitted, &[])
    .expect("a maximum-source-id packet must remain assessable");
    assert_eq!(
        assessment.disposition(),
        AssessmentDisposition::ReferenceCompleteCandidate
    );
    let decoded =
        fs_euler_disc_e2e::ClaimPolicyAssessmentLog::from_json_line(assessment.log().json_line())
            .expect("the strict log reader must admit the writer's maximum source-id row");
    assert_eq!(decoded, assessment.log().clone());

    let record = records(&contract, claim)
        .into_iter()
        .find(|record| record.requirement() == selected)
        .expect("selected evidence record");
    let error = EvidenceRecord::try_new(
        record.contract_identity(),
        record.claim(),
        record.requirement(),
        record.qois().to_vec(),
        record.authority().clone(),
        record.artifact_hash(),
        "s".repeat(MAX_PROTOCOL_ID_BYTES + 1),
        record.source_schema(),
        record.source_kind(),
        record.schema_admission_receipt_hash(),
        record.access_class(),
        record.independent(),
    )
    .expect_err("source-id maximum plus one must refuse before log production");
    assert_eq!(error.code(), "EulerProtocolInvalidIdentity");
}

#[test]
fn g0_canonical_permutations_preserve_the_exact_contract_identity() {
    let base = build_frozen_contract().expect("frozen contract");
    let mut users = base.extension().users().to_vec();
    users.reverse();
    let mut alternatives = base.extension().decision_alternatives().to_vec();
    alternatives.reverse();
    let mut risks = base.extension().risks().to_vec();
    risks.reverse();
    let mut sources = base.extension().hypothesis_sources().to_vec();
    sources.reverse();
    let extension = EulerContextExtension::try_new(
        users,
        base.extension().apparatus_population(),
        base.extension().environment_population(),
        base.extension().observation_frame(),
        alternatives,
        risks,
        sources,
    )
    .expect("permuted extension");

    let mut claims = base
        .claim_graph()
        .claims()
        .values()
        .map(|claim| {
            let mut campaign = claim.campaign().clone();
            campaign.qois.reverse();
            campaign.evidence_gaps.reverse();
            let mut requirements = claim.requirements().to_vec();
            requirements.reverse();
            EulerClaimSpec::try_new(claim.kind(), campaign, requirements).expect("permuted claim")
        })
        .collect::<Vec<_>>();
    claims.reverse();
    let mut dependencies = base.claim_graph().dependencies().to_vec();
    dependencies.reverse();
    let graph = EulerClaimGraph::try_new(claims, dependencies).expect("permuted graph");

    let mut owner_rows = base
        .owner_matrix()
        .rows()
        .values()
        .cloned()
        .collect::<Vec<_>>();
    owner_rows.reverse();
    let owner_matrix = OwnerMatrix::try_new(owner_rows).expect("permuted owners");
    let mut no_claims = base.no_claims().entries().to_vec();
    no_claims.reverse();
    let no_claim_refs = no_claims.iter().map(String::as_str).collect::<Vec<_>>();
    let mut qois = base.context().qois().values().cloned().collect::<Vec<_>>();
    qois.reverse();
    let context = context_with(
        &base,
        base.context().decision().to_owned(),
        qois,
        base.context().applicability().clone(),
    );
    let permuted = EulerScientificContract::try_new(
        context,
        extension,
        graph,
        NoClaimBoundary::new(&no_claim_refs).expect("permuted no claims"),
        owner_matrix,
    )
    .expect("permuted contract");
    assert_eq!(permuted.identity(), base.identity());
    assert_eq!(
        permuted.canonical_bytes().expect("permuted bytes"),
        base.canonical_bytes().expect("base bytes")
    );
}

#[test]
fn g0_graph_rejects_policy_weakening_duplicate_endpoint_roles_and_cycles() {
    let base = build_frozen_contract().expect("frozen contract");
    let numerical = base
        .claim_graph()
        .claim(EulerClaimKind::NumericalTrajectoryVerification)
        .expect("claim");
    let error = EulerClaimSpec::try_new(
        numerical.kind(),
        numerical.campaign().clone(),
        vec![EvidenceRequirement::CodeVerification],
    )
    .expect_err("policy weakening");
    assert_eq!(error.code(), "EulerContractClaimEvidencePolicyMismatch");

    let claims = base.claim_graph().claims().values().cloned().collect();
    let prerequisite =
        CampaignClaimId::try_new(EulerClaimKind::NumericalTrajectoryVerification.id())
            .expect("claim id");
    let dependent =
        CampaignClaimId::try_new(EulerClaimKind::CalibratedReproduction.id()).expect("claim id");
    let collision = vec![
        ClaimDependency {
            prerequisite: prerequisite.clone(),
            dependent: dependent.clone(),
            use_kind: EvidenceUse::CalibrationInput,
        },
        ClaimDependency {
            prerequisite: prerequisite.clone(),
            dependent: dependent.clone(),
            use_kind: EvidenceUse::ValidationInput,
        },
    ];
    let error = EulerClaimGraph::try_new(claims, collision).expect_err("role collision");
    assert_eq!(error.code(), "EulerContractDependencyRoleCollision");

    let claims = base.claim_graph().claims().values().cloned().collect();
    let cycle = vec![
        ClaimDependency {
            prerequisite: prerequisite.clone(),
            dependent: dependent.clone(),
            use_kind: EvidenceUse::ValidationInput,
        },
        ClaimDependency {
            prerequisite: dependent,
            dependent: prerequisite,
            use_kind: EvidenceUse::ValidationInput,
        },
    ];
    let error = EulerClaimGraph::try_new(claims, cycle).expect_err("cycle");
    assert_eq!(error.code(), "EulerContractClaimCycle");
}

#[test]
fn g0_owner_counterfeits_and_cross_schema_evidence_refuse() {
    let error = OwnerRow::try_new(
        OwnerRole::ContextOfUse,
        "fs-evidence",
        "counterfeit-context-schema-v99",
        AuthorityCeiling::StructuralContextDeclaration,
    )
    .expect_err("counterfeit owner schema");
    assert_eq!(error.code(), "EulerContractGenericSchemaFork");

    let contract = build_frozen_contract().expect("frozen contract");
    let claim = EulerClaimKind::BlindTrajectoryPrediction;
    let qois = contract
        .claim_graph()
        .claim(claim)
        .expect("claim")
        .campaign()
        .qois
        .clone();
    let error = EvidenceRecord::try_new(
        contract.identity(),
        claim,
        EvidenceRequirement::BlindHoldout,
        qois,
        EvidenceAuthorityDeclaration::StructuralProcess {
            receipt_hash: hash("blind-process-receipt"),
        },
        hash("wrong-schema"),
        "wrong-schema-source",
        "counterfeit-vv-artifact-family-v99",
        EvidenceRequirement::BlindHoldout.source_kind(),
        hash("schema-admission-receipt"),
        DeclaredEvidenceAccessClass::BlindHoldout,
        true,
    )
    .expect_err("calibration schema cannot satisfy blind role");
    assert_eq!(error.code(), "EulerProtocolSourceSchemaMismatch");
}

#[test]
fn g0_transport_decoders_refuse_malformed_truncated_and_oversized_inputs() {
    let contract = build_frozen_contract().expect("frozen contract");
    let graph_bytes = contract
        .claim_graph()
        .canonical_bytes()
        .expect("graph bytes");
    for end in 0..graph_bytes.len() {
        assert!(
            EulerClaimGraph::from_canonical_bytes(&graph_bytes[..end]).is_err(),
            "truncation at {end} unexpectedly decoded"
        );
    }
    let mut trailing = graph_bytes.clone();
    trailing.push(0);
    assert!(EulerClaimGraph::from_canonical_bytes(&trailing).is_err());
    let mut wrong_magic = graph_bytes.clone();
    wrong_magic[0] ^= 0xff;
    assert!(EulerClaimGraph::from_canonical_bytes(&wrong_magic).is_err());
    let mut unknown_kind = graph_bytes.clone();
    unknown_kind[16] = 0xff;
    assert_eq!(
        EulerClaimGraph::from_canonical_bytes(&unknown_kind)
            .expect_err("unknown claim kind")
            .code(),
        "EulerContractUnknownClaimKind"
    );
    let mut graph_at_byte_cap = vec![0; MAX_EULER_GRAPH_BYTES];
    let error = EulerClaimGraph::from_canonical_bytes(&graph_at_byte_cap)
        .expect_err("exact-cap hostile graph bytes must reach structural decoding");
    assert_eq!(error.code(), "EulerContractMalformedCanonical");
    graph_at_byte_cap.push(0);
    let error = EulerClaimGraph::from_canonical_bytes(&graph_at_byte_cap)
        .expect_err("maximum-plus-one graph bytes must refuse at preflight");
    assert_eq!(error.code(), "EulerContractGraphTooLarge");

    let bytes = contract.canonical_bytes().expect("contract bytes");
    for end in [0, 1, 7, 8, 11, bytes.len() - 1] {
        assert!(EulerScientificContract::from_canonical_bytes(&bytes[..end]).is_err());
    }
    let mut wrong_version = bytes.clone();
    wrong_version[8..12].copy_from_slice(&2_u32.to_le_bytes());
    assert_eq!(
        EulerScientificContract::from_canonical_bytes(&wrong_version)
            .expect_err("future contract")
            .code(),
        "EulerContractUnsupportedVersion"
    );
    let mut trailing = bytes;
    trailing.push(0);
    assert!(EulerScientificContract::from_canonical_bytes(&trailing).is_err());

    for length in 0..256_usize {
        let hostile = (0..length)
            .map(|index| ((index * 73 + length * 19) & 0xff) as u8)
            .collect::<Vec<_>>();
        let result = std::panic::catch_unwind(|| EulerClaimGraph::from_canonical_bytes(&hostile));
        assert!(result.is_ok(), "decoder panicked for length {length}");
        assert!(result.expect("no panic").is_err());
    }
}

#[test]
fn g0_check_receipts_logs_and_domain_separation_are_exact_and_bounded() {
    let contract = build_frozen_contract().expect("frozen contract");
    let admitted = admit_frozen_contract(contract.clone()).expect("structural admission");
    admitted
        .receipt()
        .verify_subject(admitted.contract())
        .expect("exact receipt");
    assert_ne!(
        contract.context_hash(),
        contract.claim_graph().content_hash().expect("graph")
    );
    assert_ne!(contract.context_hash(), contract.identity().as_hash());
    assert_ne!(
        contract.claim_graph().content_hash().expect("graph"),
        contract.identity().as_hash()
    );

    let claim = EulerClaimKind::BlindTrajectoryPrediction;
    let packet = packet(
        &contract,
        claim,
        records(&contract, claim),
        false,
        ReportedScientificDisposition::Positive,
        AssessmentDisposition::ReferenceCompleteCandidate,
    );
    let packet_for_log = packet.clone();
    let first = assess_with_frozen_prerequisites(&contract, packet.clone());
    let second = assess_with_frozen_prerequisites(&contract, packet);
    assert_eq!(first, second);
    assert_eq!(first.log().json_line().as_bytes().last(), Some(&b'\n'));
    assert!(first.log().json_line().len() <= 32 * 1024);
    for field in [
        "\"schema_version\"",
        "\"protocol_id\"",
        "\"contract_identity\"",
        "\"packet_contract_identity\"",
        "\"case_id\"",
        "\"packet_source_id\"",
        "\"packet_source_schema\"",
        "\"evidence_sources\"",
        "\"units\"",
        "\"seed\"",
        "\"budgets\"",
        "\"expected_disposition\"",
        "\"observed_disposition\"",
        "\"reported_scientific_disposition\"",
        "\"first_divergence\"",
        "\"reasons\"",
        "\"authority_state\"",
        "\"no_claim_state\"",
        "\"relative_artifacts\"",
        "\"reproduction_command\"",
        "\"redaction\"",
    ] {
        assert!(first.log().json_line().contains(field), "missing {field}");
    }
    assert!(!first.log().json_line().contains("youtube.com"));
    assert!(!first.log().json_line().contains("raw_observation"));
    assert!(!first.log().json_line().contains("source_prose"));
    assert!(first.log().json_line().contains(
        "\"redaction\":\"bounded-structured-protocol-metadata-no-raw-payload-or-artifact-bytes\""
    ));
    let mut expected_evidence_entries = Vec::new();
    for record in packet_for_log.records().values() {
        assert!(first.log().json_line().contains(&format!(
            "{}:{}:{}:{}",
            record.requirement().code(),
            record.source_kind().slug(),
            record.source_schema(),
            record.source_id()
        )));
        for (slot, value) in [
            ("artifact", record.artifact_hash()),
            (
                "schema-admission-receipt",
                record.schema_admission_receipt_hash(),
            ),
        ] {
            expected_evidence_entries.push(format!(
                "evidence:{}:{slot}:{}",
                record.requirement().code(),
                value.to_hex()
            ));
            assert!(
                first.log().json_line().contains(&format!(
                    "evidence:{}:{slot}:{}",
                    record.requirement().code(),
                    value.to_hex()
                )),
                "retained log must name every {slot} identity for {}",
                record.requirement().code()
            );
        }
        let role_prefix = format!("evidence:{}:role-receipt:", record.requirement().code());
        match record.authority() {
            EvidenceAuthorityDeclaration::StructuralProcess { receipt_hash } => {
                expected_evidence_entries.push(format!("{role_prefix}{}", receipt_hash.to_hex()));
                assert!(
                    first
                        .log()
                        .json_line()
                        .contains(&format!("{role_prefix}{}", receipt_hash.to_hex())),
                    "structural evidence must retain its role-receipt identity"
                );
            }
            EvidenceAuthorityDeclaration::VerifiedNumerics { .. }
            | EvidenceAuthorityDeclaration::ValidatedPhysical { .. } => assert!(
                !first.log().json_line().contains(&role_prefix),
                "a role-receipt entry is absent only when the authority declaration has no such slot"
            ),
        }
    }
    expected_evidence_entries.sort();
    let mut prior_offset = 0;
    for entry in expected_evidence_entries {
        let relative = first.log().json_line()[prior_offset..]
            .find(&entry)
            .expect("every labeled evidence identity must be retained in canonical order");
        prior_offset += relative + entry.len();
    }
    let prerequisite_route =
        contract.owner_matrix().rows()[&OwnerRole::PrerequisiteAssessmentReceipt].source_schema();
    assert!(first.log().json_line().contains(prerequisite_route));
    assert_eq!(
        prerequisite_route,
        EULER_PREREQUISITE_RECEIPT_IDENTITY_DOMAIN
    );
    assert_eq!(
        contract.owner_matrix().rows()[&OwnerRole::ClaimPolicyAssessmentLog].source_schema(),
        CLAIM_POLICY_ASSESSMENT_LOG_DOMAIN
    );
    let log_declaration =
        fs_euler_disc_e2e::protocol::CLAIM_POLICY_ASSESSMENT_LOG_IDENTITY_SCHEMA_DECLARATION;
    assert!(log_declaration.iter().any(|line| {
        *line
            == "schema_dependencies=fs-euler-disc-e2e:claim-evidence-packet,fs-euler-disc-e2e:scientific-contract,fs-euler-disc-e2e:owner-matrix"
    }));
    assert!(!log_declaration.iter().any(|line| {
        line.starts_with("schema_dependencies=")
            && line.contains("fs-euler-disc-e2e:prerequisite-assessment-receipt")
    }));
    assert!(first.log().json_line().contains(&format!(
        "\"reproduction_command\":\"{}\"",
        fs_euler_disc_e2e::protocol::ASSESSMENT_LOG_REPRODUCTION_COMMAND
    )));
    assert!(
        include_str!("../../../scripts/ci/euler_disc_contract_e2e.sh")
            .contains(fs_euler_disc_e2e::protocol::ASSESSMENT_LOG_REPRODUCTION_COMMAND)
    );
    assert_eq!(
        first.log().identity(),
        fs_blake3::hash_domain(
            fs_euler_disc_e2e::protocol::CLAIM_POLICY_ASSESSMENT_LOG_DOMAIN,
            first.log().json_line().as_bytes()
        )
    );
}

#[test]
fn g0_frozen_contract_digests_are_literal_anchors() {
    let contract = build_frozen_contract().expect("frozen contract");
    let actual = (
        contract.context_hash().to_hex(),
        contract
            .claim_graph()
            .content_hash()
            .expect("claim-graph identity")
            .to_hex(),
        contract.identity().as_hash().to_hex(),
    );
    let reviewed = (
        FROZEN_CONTEXT_HASH_HEX.to_owned(),
        FROZEN_CLAIM_GRAPH_HASH_HEX.to_owned(),
        FROZEN_CONTRACT_IDENTITY_HEX.to_owned(),
    );
    assert_eq!(
        actual, reviewed,
        "the Context, claim graph, and complete contract anchors move as one reviewed set"
    );
}

#[test]
fn g0_derived_numeric_qois_apply_ranges_to_the_qoi_not_error_of_error() {
    let contract = build_frozen_contract().expect("frozen contract");
    let expected = [
        ("energy-balance-residual", -0.001, 0.001),
        ("energy-channel-fraction-error", 0.0, 0.05),
        ("event-time-error", 0.0, 0.5),
        ("normalized-trajectory-discrepancy", 0.0, 0.05),
        ("numerical-trajectory-error", 0.0, 1.0e-8),
        ("optimum-interval-width", 0.0, 0.25),
    ];
    for (id, expected_lo, expected_hi) in expected {
        let qoi = &contract.context().qois()
            [&fs_evidence::vv::QoiId::try_new(id).expect("reviewed QoI id")];
        let AcceptanceCriterion::ClosedRange { lo, hi } = qoi.acceptance() else {
            panic!("derived numeric QoI {id} must use a direct closed range");
        };
        assert_eq!((*lo, *hi), (expected_lo, expected_hi));
    }
}

#[test]
fn g0_all_categorical_claim_dispositions_are_explicit_qois() {
    let contract = build_frozen_contract().expect("frozen contract");
    let expected = [
        (
            "event-class-disposition",
            "matches-preregistered-event-class",
        ),
        (
            "qualitative-effect-disposition",
            "matches-preregistered-direction",
        ),
        (
            "configuration-ranking-disposition",
            "matches-preregistered-order-or-tie-rule",
        ),
        (
            "optimum-containment-disposition",
            "contains-preregistered-optimum-under-exact-score",
        ),
        (
            "rival-mechanism-disposition",
            "discriminates-preregistered-rival",
        ),
    ];
    for (id, expected_category) in expected {
        let qoi = &contract.context().qois()
            [&fs_evidence::vv::QoiId::try_new(id).expect("reviewed QoI id")];
        assert_eq!(qoi.unit().as_str(), "1");
        let AcceptanceCriterion::CategoryEquals { expected } = qoi.acceptance() else {
            panic!("categorical disposition QoI {id} must use CategoryEquals");
        };
        assert_eq!(expected, expected_category);
    }
}

#[test]
fn hypothesis_source_identity_mutation_battery() {
    let base = HypothesisSource::try_new("source-ab", "locator-c").expect("source");
    let changed_id = HypothesisSource::try_new("source-ac", "locator-c").expect("source");
    let changed_locator = HypothesisSource::try_new("source-ab", "locator-d").expect("source");
    let framed_left = HypothesisSource::try_new("ab", "c").expect("source");
    let framed_right = HypothesisSource::try_new("a", "bc").expect("source");
    for source in [
        &base,
        &changed_id,
        &changed_locator,
        &framed_left,
        &framed_right,
    ] {
        source
            .verify_identity()
            .expect("exact declaration identity");
    }
    assert_ne!(base.declaration_hash(), changed_id.declaration_hash());
    assert_ne!(base.declaration_hash(), changed_locator.declaration_hash());
    assert_ne!(
        framed_left.declaration_hash(),
        framed_right.declaration_hash()
    );
}

#[test]
fn claim_graph_identity_mutation_battery() {
    let contract = build_frozen_contract().expect("frozen contract");
    let graph = contract.claim_graph();
    let bytes = graph.canonical_bytes().expect("graph bytes");
    let identity = graph.content_hash().expect("graph identity");
    assert_eq!(
        identity,
        fs_blake3::hash_domain(EULER_CLAIM_GRAPH_IDENTITY_DOMAIN, &bytes)
    );
    let decoded = EulerClaimGraph::from_canonical_bytes(&bytes).expect("graph fixed point");
    assert_eq!(decoded.content_hash().expect("decoded identity"), identity);

    for index in [0, 8, 12, bytes.len() / 2, bytes.len() - 1] {
        let mut mutated = bytes.clone();
        mutated[index] ^= 1;
        let raw_identity = fs_blake3::hash_domain(EULER_CLAIM_GRAPH_IDENTITY_DOMAIN, &mutated);
        assert_ne!(raw_identity, identity);
        if let Ok(candidate) = EulerClaimGraph::from_canonical_bytes(&mutated) {
            assert_ne!(
                candidate.content_hash().expect("candidate identity"),
                identity
            );
        }
    }
    g0_graph_rejects_policy_weakening_duplicate_endpoint_roles_and_cycles();
    g0_canonical_permutations_preserve_the_exact_contract_identity();
}

#[test]
fn scientific_contract_identity_mutation_battery() {
    let contract = build_frozen_contract().expect("frozen contract");
    let bytes = contract.canonical_bytes().expect("contract bytes");
    assert_eq!(
        contract.identity().as_hash(),
        fs_blake3::hash_domain(EULER_CONTRACT_IDENTITY_DOMAIN, &bytes)
    );
    assert_eq!(
        EulerScientificContract::from_canonical_bytes(&bytes)
            .expect("contract fixed point")
            .identity(),
        contract.identity()
    );
    g3_every_required_context_semantic_moves_identity_and_stales_receipts();
}

#[test]
fn claim_evidence_packet_identity_mutation_battery() {
    let contract = build_frozen_contract().expect("frozen contract");
    let claim = EulerClaimKind::NumericalTrajectoryVerification;
    let base_point = point(&contract);
    let base_budget = ProtocolBudget::try_new(60_000, 64 * 1024 * 1024, 0.05).expect("base budget");
    let base = packet_from_parts(
        contract.identity(),
        "identity-case-a",
        claim,
        base_point.clone(),
        Vec::new(),
        true,
        false,
        ReportedScientificDisposition::Positive,
        AssessmentDisposition::Refused,
        vec!["identity-unit-a".to_owned()],
        ProtocolSeed::Fixed { value: 7 },
        base_budget,
    );

    let alternate_context = context_with(
        &contract,
        format!("{} {}", contract.context().decision(), "Identity mutation."),
        contract.context().qois().values().cloned().collect(),
        contract.context().applicability().clone(),
    );
    let alternate_contract = contract_with_context(&contract, alternate_context);
    let packet_with_bindings = |design_set_identity, aggregate_qoi_derivation_receipt_identity| {
        ClaimEvidencePacket::try_new(
            contract.identity(),
            "identity-case-a",
            design_set_identity,
            aggregate_qoi_derivation_receipt_identity,
            claim,
            base_point.clone(),
            Vec::new(),
            true,
            false,
            ReportedScientificDisposition::Positive,
            AssessmentDisposition::Refused,
            vec!["identity-unit-a".to_owned()],
            ProtocolSeed::Fixed { value: 7 },
            base_budget,
        )
        .expect("identity binding variant")
    };
    let variants = vec![
        (
            "contract-identity",
            packet_from_parts(
                alternate_contract.identity(),
                "identity-case-a",
                claim,
                base_point.clone(),
                Vec::new(),
                true,
                false,
                ReportedScientificDisposition::Positive,
                AssessmentDisposition::Refused,
                vec!["identity-unit-a".to_owned()],
                ProtocolSeed::Fixed { value: 7 },
                base_budget,
            ),
        ),
        (
            "case-id",
            packet_from_parts(
                contract.identity(),
                "identity-case-b",
                claim,
                base_point.clone(),
                Vec::new(),
                true,
                false,
                ReportedScientificDisposition::Positive,
                AssessmentDisposition::Refused,
                vec!["identity-unit-a".to_owned()],
                ProtocolSeed::Fixed { value: 7 },
                base_budget,
            ),
        ),
        (
            "design-set-identity",
            packet_with_bindings(
                hash("alternate-design-set"),
                hash("shared-aggregate-qoi-derivation-receipt"),
            ),
        ),
        (
            "aggregate-qoi-derivation-receipt-identity",
            packet_with_bindings(
                hash("shared-design-set"),
                hash("alternate-aggregate-qoi-derivation-receipt"),
            ),
        ),
        (
            "claim-kind",
            packet_from_parts(
                contract.identity(),
                "identity-case-a",
                EulerClaimKind::CalibratedReproduction,
                base_point.clone(),
                Vec::new(),
                true,
                false,
                ReportedScientificDisposition::Positive,
                AssessmentDisposition::Refused,
                vec!["identity-unit-a".to_owned()],
                ProtocolSeed::Fixed { value: 7 },
                base_budget,
            ),
        ),
        (
            "applicability-point-anchor",
            packet_from_parts(
                contract.identity(),
                "identity-case-a",
                claim,
                point_at_fraction(&contract, 0.25),
                Vec::new(),
                true,
                false,
                ReportedScientificDisposition::Positive,
                AssessmentDisposition::Refused,
                vec!["identity-unit-a".to_owned()],
                ProtocolSeed::Fixed { value: 7 },
                base_budget,
            ),
        ),
        (
            "evidence-registry",
            packet_from_parts(
                contract.identity(),
                "identity-case-a",
                claim,
                base_point.clone(),
                records(&contract, claim),
                true,
                false,
                ReportedScientificDisposition::Positive,
                AssessmentDisposition::Refused,
                vec!["identity-unit-a".to_owned()],
                ProtocolSeed::Fixed { value: 7 },
                base_budget,
            ),
        ),
        (
            "no-claim-acceptance",
            packet_from_parts(
                contract.identity(),
                "identity-case-a",
                claim,
                base_point.clone(),
                Vec::new(),
                false,
                false,
                ReportedScientificDisposition::Positive,
                AssessmentDisposition::Refused,
                vec!["identity-unit-a".to_owned()],
                ProtocolSeed::Fixed { value: 7 },
                base_budget,
            ),
        ),
        (
            "target-fitting-state",
            packet_from_parts(
                contract.identity(),
                "identity-case-a",
                claim,
                base_point.clone(),
                Vec::new(),
                true,
                true,
                ReportedScientificDisposition::Positive,
                AssessmentDisposition::Refused,
                vec!["identity-unit-a".to_owned()],
                ProtocolSeed::Fixed { value: 7 },
                base_budget,
            ),
        ),
        (
            "reported-scientific-disposition",
            packet_from_parts(
                contract.identity(),
                "identity-case-a",
                claim,
                base_point.clone(),
                Vec::new(),
                true,
                false,
                ReportedScientificDisposition::Negative,
                AssessmentDisposition::Refused,
                vec!["identity-unit-a".to_owned()],
                ProtocolSeed::Fixed { value: 7 },
                base_budget,
            ),
        ),
        (
            "expected-disposition",
            packet_from_parts(
                contract.identity(),
                "identity-case-a",
                claim,
                base_point.clone(),
                Vec::new(),
                true,
                false,
                ReportedScientificDisposition::Positive,
                AssessmentDisposition::DemotedCandidate,
                vec!["identity-unit-a".to_owned()],
                ProtocolSeed::Fixed { value: 7 },
                base_budget,
            ),
        ),
        (
            "unit-set",
            packet_from_parts(
                contract.identity(),
                "identity-case-a",
                claim,
                base_point.clone(),
                Vec::new(),
                true,
                false,
                ReportedScientificDisposition::Positive,
                AssessmentDisposition::Refused,
                vec!["identity-unit-b".to_owned()],
                ProtocolSeed::Fixed { value: 7 },
                base_budget,
            ),
        ),
        (
            "seed-declaration",
            packet_from_parts(
                contract.identity(),
                "identity-case-a",
                claim,
                base_point.clone(),
                Vec::new(),
                true,
                false,
                ReportedScientificDisposition::Positive,
                AssessmentDisposition::Refused,
                vec!["identity-unit-a".to_owned()],
                ProtocolSeed::Fixed { value: 8 },
                base_budget,
            ),
        ),
        (
            "protocol-budget",
            packet_from_parts(
                contract.identity(),
                "identity-case-a",
                claim,
                base_point.clone(),
                Vec::new(),
                true,
                false,
                ReportedScientificDisposition::Positive,
                AssessmentDisposition::Refused,
                vec!["identity-unit-a".to_owned()],
                ProtocolSeed::Fixed { value: 7 },
                ProtocolBudget::try_new(60_001, 64 * 1024 * 1024, 0.05).expect("alternate budget"),
            ),
        ),
    ];
    base.verify_identity().expect("base packet identity");
    for (field, variant) in variants {
        variant.verify_identity().expect("variant packet identity");
        assert_ne!(
            variant.identity(),
            base.identity(),
            "semantic packet field {field} did not move identity"
        );
    }

    let bytes = base.canonical_bytes().expect("packet bytes");
    assert_eq!(
        base.identity(),
        fs_blake3::hash_domain(EULER_EVIDENCE_PACKET_IDENTITY_DOMAIN, &bytes)
    );
    assert_ne!(
        base.identity(),
        fs_blake3::hash_domain(
            "org.frankensim.fs-euler-disc-e2e.claim-evidence-packet.v2",
            &bytes
        ),
        "identity domain/version must move the packet identity"
    );
    for (field, mutated) in [
        ("transport-magic", {
            let mut mutated = bytes.clone();
            mutated[0] ^= 1;
            mutated
        }),
        ("protocol-schema-version", {
            let mut mutated = bytes.clone();
            mutated[8..12].copy_from_slice(&2_u32.to_le_bytes());
            mutated
        }),
        ("fixed-numeric-little-endian", {
            let mut mutated = bytes.clone();
            mutated[8..12].copy_from_slice(&1_u32.to_be_bytes());
            mutated
        }),
        ("canonical-field-order", {
            let mut mutated = Vec::with_capacity(bytes.len());
            mutated.extend_from_slice(&bytes[..8]);
            mutated.extend_from_slice(&bytes[12..44]);
            mutated.extend_from_slice(&bytes[8..12]);
            mutated.extend_from_slice(&bytes[44..]);
            mutated
        }),
        ("length-framing", {
            let mut mutated = Vec::with_capacity(bytes.len() - 4);
            mutated.extend_from_slice(&bytes[..44]);
            mutated.extend_from_slice(&bytes[48..]);
            mutated
        }),
    ] {
        assert_ne!(
            base.identity(),
            fs_blake3::hash_domain(EULER_EVIDENCE_PACKET_IDENTITY_DOMAIN, &mutated),
            "encoding semantic {field} did not move packet identity"
        );
    }
}

#[test]
fn prerequisite_receipt_identity_mutation_battery() {
    let contract = build_frozen_contract().expect("frozen contract");
    let positive = synthetic_reported_positive_assessment_map(&contract);
    let numerical = &positive[&EulerClaimKind::NumericalTrajectoryVerification];
    let receipts = [
        numerical
            .as_prerequisite_for(
                EulerClaimKind::CalibratedReproduction,
                EvidenceUse::ValidationInput,
            )
            .expect("direct receipt"),
        numerical
            .as_prerequisite_for(
                EulerClaimKind::BlindTrajectoryPrediction,
                EvidenceUse::ValidationInput,
            )
            .expect("different dependent"),
        numerical
            .as_prerequisite_for(
                EulerClaimKind::CalibratedReproduction,
                EvidenceUse::CalibrationInput,
            )
            .expect("different use"),
    ];
    for receipt in &receipts {
        receipt.verify().expect("receipt identity");
    }
    assert_ne!(receipts[0].identity(), receipts[1].identity());
    assert_ne!(receipts[0].identity(), receipts[2].identity());
    assert_eq!(
        receipts[0].identity(),
        fs_blake3::hash_domain(
            EULER_PREREQUISITE_RECEIPT_IDENTITY_DOMAIN,
            &receipts[0].canonical_bytes().expect("receipt bytes")
        )
    );
}

#[test]
fn claim_policy_assessment_identity_mutation_battery() {
    let contract = build_frozen_contract().expect("frozen contract");
    let claim = EulerClaimKind::NumericalTrajectoryVerification;
    let positive = assess_with_frozen_prerequisites(
        &contract,
        packet(
            &contract,
            claim,
            records(&contract, claim),
            false,
            ReportedScientificDisposition::Positive,
            AssessmentDisposition::ReferenceCompleteCandidate,
        ),
    );
    let negative = assess_with_frozen_prerequisites(
        &contract,
        packet(
            &contract,
            claim,
            records(&contract, claim),
            false,
            ReportedScientificDisposition::Negative,
            AssessmentDisposition::RetainedTerminal,
        ),
    );
    positive.verify_identity().expect("positive identity");
    negative.verify_identity().expect("negative identity");
    assert_ne!(positive.identity(), negative.identity());
    assert_ne!(positive.log().identity(), negative.log().identity());
    assert_ne!(
        fs_blake3::hash_domain(EULER_ASSESSMENT_IDENTITY_DOMAIN, b"assessment-a"),
        fs_blake3::hash_domain(EULER_ASSESSMENT_IDENTITY_DOMAIN, b"assessment-b")
    );
}

#[test]
fn contract_check_receipt_identity_mutation_battery() {
    let contract = build_frozen_contract().expect("frozen contract");
    let receipt = check_frozen_contract(&contract).expect("check receipt");
    let bytes = receipt.canonical_bytes().expect("receipt bytes");
    let decoded = ContractCheckReceipt::from_canonical_bytes(&bytes).expect("receipt fixed point");
    assert_eq!(decoded, receipt);
    decoded.verify_identity().expect("decoded identity");
    decoded.verify_subject(&contract).expect("decoded subject");

    // A self-consistent transport is not an independent check. Forge the
    // three subject hashes in a genuine passing receipt so that it binds a
    // valid but non-frozen contract, then prove that verification re-runs the
    // literal-anchor checker instead of trusting the transported pass bit.
    let mutated = contract_with_context(
        &contract,
        context_with(
            &contract,
            format!("{} {}", contract.context().decision(), "Mutated decision."),
            contract.context().qois().values().cloned().collect(),
            contract.context().applicability().clone(),
        ),
    );
    let mut forged_bytes = bytes.clone();
    let checker_length = u32::from_le_bytes(
        forged_bytes[12..16]
            .try_into()
            .expect("checker length framing"),
    ) as usize;
    let subject_start = 16 + checker_length;
    let context_start = subject_start + 32;
    let graph_start = context_start + 32;
    forged_bytes[subject_start..context_start]
        .copy_from_slice(mutated.identity().as_hash().as_bytes());
    forged_bytes[context_start..graph_start].copy_from_slice(mutated.context_hash().as_bytes());
    forged_bytes[graph_start..graph_start + 32].copy_from_slice(
        mutated
            .claim_graph()
            .content_hash()
            .expect("mutated graph hash")
            .as_bytes(),
    );
    let forged = ContractCheckReceipt::from_canonical_bytes(&forged_bytes)
        .expect("self-consistent forged transport decodes");
    forged.verify_identity().expect("forged transport identity");
    assert!(forged.passed(), "transport still carries a forged pass bit");
    let error = forged
        .verify_subject(&mutated)
        .expect_err("literal-anchor recheck must reject forged passing receipt");
    assert_eq!(error.code(), "EulerContractStaleCheckReceipt");
    assert!(
        !check_frozen_contract(&mutated)
            .expect("mutated structural check")
            .passed()
    );

    for index in [0, 8, 12, bytes.len() / 2, bytes.len() - 1] {
        let mut mutated = bytes.clone();
        mutated[index] ^= 1;
        match ContractCheckReceipt::from_canonical_bytes(&mutated) {
            Ok(candidate) => {
                candidate.verify_identity().expect("candidate identity");
                assert_ne!(candidate.identity(), receipt.identity());
                assert!(candidate.verify_subject(&contract).is_err());
            }
            Err(error) => assert!(!error.code().is_empty()),
        }
    }
    let mut trailing = bytes.clone();
    trailing.push(0);
    assert!(ContractCheckReceipt::from_canonical_bytes(&trailing).is_err());
    for end in 0..bytes.len() {
        assert!(ContractCheckReceipt::from_canonical_bytes(&bytes[..end]).is_err());
    }
}

#[test]
fn claim_policy_assessment_log_identity_mutation_battery() {
    let contract = build_frozen_contract().expect("frozen contract");
    let claim = EulerClaimKind::NumericalTrajectoryVerification;
    let positive = assess_with_frozen_prerequisites(
        &contract,
        packet(
            &contract,
            claim,
            records(&contract, claim),
            false,
            ReportedScientificDisposition::Positive,
            AssessmentDisposition::ReferenceCompleteCandidate,
        ),
    );
    let refused = assess_with_frozen_prerequisites(
        &contract,
        packet(
            &contract,
            EulerClaimKind::BlindTrajectoryPrediction,
            records(&contract, EulerClaimKind::BlindTrajectoryPrediction),
            true,
            ReportedScientificDisposition::Positive,
            AssessmentDisposition::Refused,
        ),
    );
    let blind_positive = assess_with_frozen_prerequisites(
        &contract,
        packet(
            &contract,
            EulerClaimKind::BlindTrajectoryPrediction,
            records(&contract, EulerClaimKind::BlindTrajectoryPrediction),
            false,
            ReportedScientificDisposition::Positive,
            AssessmentDisposition::ReferenceCompleteCandidate,
        ),
    );
    let multi_refused = assess_with_frozen_prerequisites(
        &contract,
        packet(
            &contract,
            EulerClaimKind::QualitativeEffectDirection,
            Vec::new(),
            true,
            ReportedScientificDisposition::Positive,
            AssessmentDisposition::Refused,
        ),
    );
    let weak_requirement = golden_requirements(claim)[0];
    let demoted = assess_with_frozen_prerequisites(
        &contract,
        packet(
            &contract,
            claim,
            golden_requirements(claim)
                .iter()
                .map(|requirement| {
                    record_with(
                        &contract,
                        claim,
                        *requirement,
                        *requirement == weak_requirement,
                        hash(&format!("log-demotion:{}", requirement.code())),
                        declared_access_class(*requirement),
                    )
                })
                .collect(),
            false,
            ReportedScientificDisposition::Positive,
            AssessmentDisposition::DemotedCandidate,
        ),
    );
    let admitted = admit_frozen_contract(contract.clone()).expect("structural admission");
    let single_refused = packet_with_policy_declarations(
        &contract,
        claim,
        point(&contract),
        records(&contract, claim),
        true,
        false,
        ReportedScientificDisposition::Positive,
        AssessmentDisposition::Refused,
        vec!["wrong-unit".to_owned()],
    )
    .assess(&admitted, &[])
    .expect("single-reason refusal");
    let calibrated_refused = assess_with_frozen_prerequisites(
        &contract,
        packet_with_policy_declarations(
            &contract,
            EulerClaimKind::CalibratedReproduction,
            point(&contract),
            records(&contract, EulerClaimKind::CalibratedReproduction),
            true,
            false,
            ReportedScientificDisposition::Positive,
            AssessmentDisposition::Refused,
            vec!["wrong-unit".to_owned()],
        ),
    );
    let numeric_axis = contract
        .context()
        .applicability()
        .numeric()
        .keys()
        .next()
        .expect("frozen context has a numeric axis")
        .clone();
    let missing_numeric_point = ApplicabilityPoint::try_new(
        point(&contract)
            .numeric()
            .iter()
            .filter(|(axis, _)| *axis != &numeric_axis)
            .map(|(axis, value)| (axis.clone(), *value))
            .collect(),
        point(&contract)
            .categorical()
            .iter()
            .map(|(axis, value)| (axis.clone(), value.clone()))
            .collect(),
    )
    .expect("point with one omitted numeric context axis");
    let missing_numeric_refused = assess_with_frozen_prerequisites(
        &contract,
        packet_at_point(
            &contract,
            claim,
            missing_numeric_point,
            records(&contract, claim),
            false,
            ReportedScientificDisposition::Positive,
            AssessmentDisposition::Refused,
        ),
    );
    let categorical_axis = contract
        .context()
        .applicability()
        .categorical()
        .keys()
        .next()
        .expect("frozen context has a categorical axis")
        .clone();
    let missing_categorical_point = ApplicabilityPoint::try_new(
        point(&contract)
            .numeric()
            .iter()
            .map(|(axis, value)| (axis.clone(), *value))
            .collect(),
        point(&contract)
            .categorical()
            .iter()
            .filter(|(axis, _)| *axis != &categorical_axis)
            .map(|(axis, value)| (axis.clone(), value.clone()))
            .collect(),
    )
    .expect("point with one omitted categorical context axis");
    let missing_categorical_refused = assess_with_frozen_prerequisites(
        &contract,
        packet_at_point(
            &contract,
            claim,
            missing_categorical_point,
            records(&contract, claim),
            false,
            ReportedScientificDisposition::Positive,
            AssessmentDisposition::Refused,
        ),
    );
    let access_requirement = EvidenceRequirement::CodeVerification;
    let access_refused = assess_with_frozen_prerequisites(
        &contract,
        packet(
            &contract,
            claim,
            golden_requirements(claim)
                .iter()
                .map(|requirement| {
                    record_with(
                        &contract,
                        claim,
                        *requirement,
                        false,
                        hash(&format!("log-access-refusal:{}", requirement.code())),
                        if *requirement == access_requirement {
                            DeclaredEvidenceAccessClass::Calibration
                        } else {
                            declared_access_class(*requirement)
                        },
                    )
                })
                .collect(),
            false,
            ReportedScientificDisposition::Positive,
            AssessmentDisposition::Refused,
        ),
    );
    let weak_refused = packet_with_policy_declarations(
        &contract,
        claim,
        point(&contract),
        golden_requirements(claim)
            .iter()
            .map(|requirement| {
                record_with(
                    &contract,
                    claim,
                    *requirement,
                    *requirement == weak_requirement,
                    hash(&format!("log-weak-refusal:{}", requirement.code())),
                    declared_access_class(*requirement),
                )
            })
            .collect(),
        true,
        false,
        ReportedScientificDisposition::Positive,
        AssessmentDisposition::Refused,
        vec!["wrong-unit".to_owned()],
    )
    .assess(&admitted, &[])
    .expect("hard refusal with mismatched authority");
    assert_eq!(
        demoted.disposition(),
        AssessmentDisposition::DemotedCandidate
    );
    assert_eq!(single_refused.disposition(), AssessmentDisposition::Refused);
    assert_eq!(single_refused.reasons().len(), 1);
    assert_eq!(calibrated_refused.reasons().len(), 1);
    assert_eq!(missing_numeric_refused.reasons().len(), 1);
    assert_eq!(missing_categorical_refused.reasons().len(), 1);
    assert_eq!(access_refused.reasons().len(), 1);
    positive.log().verify_identity().expect("positive log");
    refused.log().verify_identity().expect("refused log");
    demoted
        .log()
        .verify_identity()
        .expect("mismatched-authority demotion log");
    weak_refused
        .log()
        .verify_identity()
        .expect("mismatched-authority refusal log");
    let decoded =
        fs_euler_disc_e2e::ClaimPolicyAssessmentLog::from_json_line(positive.log().json_line())
            .expect("exact canonical assessment log");
    assert_eq!(decoded, positive.log().clone());
    fs_euler_disc_e2e::ClaimPolicyAssessmentLog::from_json_line(demoted.log().json_line())
        .expect("mismatched-authority demotion must round-trip");
    fs_euler_disc_e2e::ClaimPolicyAssessmentLog::from_json_line(weak_refused.log().json_line())
        .expect("mismatched-authority refusal must round-trip");
    let assert_malformed = |label: &str, candidate: String| {
        let Err(error) = fs_euler_disc_e2e::ClaimPolicyAssessmentLog::from_json_line(candidate)
        else {
            panic!("hostile {label} mutation was admitted");
        };
        assert_eq!(
            error.code(),
            "EulerProtocolMalformedAssessmentLog",
            "unexpected refusal for {label}: {error}"
        );
    };
    let wrong_version =
        positive
            .log()
            .json_line()
            .replacen("\"schema_version\":1", "\"schema_version\":2", 1);
    assert_malformed("wrong version", wrong_version);
    let wrong_domain = positive.log().json_line().replacen(
        CLAIM_POLICY_ASSESSMENT_LOG_DOMAIN,
        "org.frankensim.fs-euler-disc-e2e.claim-policy-assessment-log.v2",
        1,
    );
    assert_malformed("wrong domain", wrong_domain);
    let mut trailing_line = positive.log().json_line().to_owned();
    trailing_line.push_str("{}\n");
    assert_malformed("second JSON line", trailing_line);

    let original = positive.log().json_line();
    assert!(original.contains("\"normalized_accuracy_limit_bits\":"));
    assert!(!original.contains("\"accuracy_limit_bits\":"));
    let replace_field_value =
        |line: &str, field: &str, following_field: &str, replacement: &str| {
            let marker = format!("\"{field}\":");
            let start = line.find(&marker).expect("field marker") + marker.len();
            let following = format!(",\"{following_field}\":");
            let end = start + line[start..].find(&following).expect("following field");
            let mut mutated = line.to_owned();
            mutated.replace_range(start..end, replacement);
            mutated
        };
    let remove_array_entry = |line: &str, entry: &str| {
        let leading = format!("\"{entry}\",");
        if line.contains(&leading) {
            return line.replacen(&leading, "", 1);
        }
        let trailing = format!(",\"{entry}\"");
        assert!(line.contains(&trailing), "array entry must be present");
        line.replacen(&trailing, "", 1)
    };
    let add_refusal_reasons = |line: &str, existing_reason: &str, added_reasons: &[&str]| {
        let existing = format!("\"reasons\":[\"{existing_reason}\"]");
        assert!(
            line.contains(&existing),
            "refusal must have one baseline reason"
        );
        let mut reasons = Vec::with_capacity(added_reasons.len() + 1);
        reasons.push(existing_reason);
        reasons.extend_from_slice(added_reasons);
        reasons.sort_unstable();
        let replacement = reasons
            .iter()
            .map(|reason| format!("\"{reason}\""))
            .collect::<Vec<_>>()
            .join(",");
        line.replacen(&existing, &format!("\"reasons\":[{replacement}]"), 1)
    };
    let add_refusal_reason = |line: &str, existing_reason: &str, added_reason: &str| {
        add_refusal_reasons(line, existing_reason, &[added_reason])
    };

    assert_malformed(
        "substituted frozen contract identity",
        original.replacen(
            &format!("\"contract_identity\":\"{}\"", contract.identity()),
            &format!("\"contract_identity\":\"{}\"", "1".repeat(64)),
            1,
        ),
    );

    let packet_artifact = format!("packet:{}", positive.packet_identity().to_hex());
    let design_set_artifact = format!("design-set:{}", positive.design_set_identity().to_hex());
    let aggregate_qoi_derivation_artifact = format!(
        "aggregate-qoi-derivation:{}:{}",
        fs_euler_disc_e2e::contract::EULER_AGGREGATE_QOI_DERIVATION_RECEIPT_SCHEMA,
        positive
            .aggregate_qoi_derivation_receipt_identity()
            .to_hex()
    );
    assert!(original.contains(&design_set_artifact));
    assert!(original.contains(&aggregate_qoi_derivation_artifact));
    assert_malformed(
        "candidate missing its exact design-set artifact",
        remove_array_entry(original, &design_set_artifact),
    );
    assert_malformed(
        "candidate missing its aggregate-QoI derivation receipt",
        remove_array_entry(original, &aggregate_qoi_derivation_artifact),
    );
    assert_malformed(
        "aggregate-QoI derivation artifact with wrong owner route",
        original.replacen(
            fs_euler_disc_e2e::contract::EULER_AGGREGATE_QOI_DERIVATION_RECEIPT_SCHEMA,
            "org.frankensim.fs-euler-disc-e2e.wrong-qoi-derivation-route.v1",
            1,
        ),
    );
    assert_malformed(
        "non-refused zero-edge claim with an injected prerequisite",
        original.replacen(
            &packet_artifact,
            &format!(
                "{packet_artifact}\",\"prerequisite:{}:{}:{}",
                EulerClaimKind::BlindTrajectoryPrediction.id(),
                EULER_PREREQUISITE_RECEIPT_IDENTITY_DOMAIN,
                hash("injected-prerequisite").to_hex()
            ),
            1,
        ),
    );
    let refused_packet_artifact = format!("packet:{}", single_refused.packet_identity().to_hex());
    assert_malformed(
        "refused zero-edge claim with an unreasoned injected prerequisite",
        single_refused.log().json_line().replacen(
            &refused_packet_artifact,
            &format!(
                "{refused_packet_artifact}\",\"prerequisite:{}:{}:{}",
                EulerClaimKind::BlindTrajectoryPrediction.id(),
                EULER_PREREQUISITE_RECEIPT_IDENTITY_DOMAIN,
                hash("injected-refused-prerequisite").to_hex()
            ),
            1,
        ),
    );
    let injected_unexpected_claim = EulerClaimKind::BlindTrajectoryPrediction;
    let injected_with_uncorrelated_malformed = single_refused.log().json_line().replacen(
        &refused_packet_artifact,
        &format!(
            "{refused_packet_artifact}\",\"prerequisite:{}:{}:{}",
            injected_unexpected_claim.id(),
            EULER_PREREQUISITE_RECEIPT_IDENTITY_DOMAIN,
            hash("injected-malformed-and-stale-prerequisite").to_hex()
        ),
        1,
    );
    assert_malformed(
        "claim-specific stale unexpected prerequisite hidden behind an uncorrelated malformed reason",
        add_refusal_reasons(
            &injected_with_uncorrelated_malformed,
            single_refused
                .reasons()
                .first()
                .expect("unit refusal has one reason"),
            &[
                "malformed-prerequisite-receipt:EulerProtocolPrerequisiteReceiptIdentityMismatch",
                &format!(
                    "stale-prerequisite-receipt:{}",
                    injected_unexpected_claim.id()
                ),
            ],
        ),
    );
    let calibrated_refusal_reason = calibrated_refused
        .reasons()
        .first()
        .expect("calibrated unit refusal has one reason");
    let reachable_malformed_prerequisite = add_refusal_reasons(
        calibrated_refused.log().json_line(),
        calibrated_refusal_reason,
        &["malformed-prerequisite-receipt:EulerProtocolPrerequisiteReceiptIdentityMismatch"],
    );
    fs_euler_disc_e2e::ClaimPolicyAssessmentLog::from_json_line(&reachable_malformed_prerequisite)
        .expect("reachable malformed-prerequisite code with a retained row must remain admissible");
    assert_malformed(
        "producer-unreachable prerequisite failure code",
        add_refusal_reasons(
            calibrated_refused.log().json_line(),
            calibrated_refusal_reason,
            &["malformed-prerequisite-receipt:EulerContractCheckReceiptTooLarge"],
        ),
    );
    let blind_prerequisite = blind_positive
        .log()
        .json_line()
        .split('"')
        .find(|value| value.starts_with("prerequisite:"))
        .expect("blind candidate retains direct prerequisites");
    assert_malformed(
        "non-refused blind claim with a removed prerequisite",
        remove_array_entry(blind_positive.log().json_line(), blind_prerequisite),
    );

    let packet_only_artifacts = format!("[\"packet:{}\"]", positive.packet_identity().to_hex());
    let without_evidence_sources = replace_field_value(original, "evidence_sources", "units", "[]");
    assert_malformed(
        "candidate with amputated evidence source registry",
        without_evidence_sources.clone(),
    );
    assert_malformed(
        "candidate with packet-only relative artifact registry",
        replace_field_value(
            original,
            "relative_artifacts",
            "reproduction_command",
            &packet_only_artifacts,
        ),
    );
    assert_malformed(
        "candidate with both evidence registries amputated",
        replace_field_value(
            &without_evidence_sources,
            "relative_artifacts",
            "reproduction_command",
            &packet_only_artifacts,
        ),
    );

    let source_record = records(&contract, claim)
        .into_iter()
        .next()
        .expect("positive claim has evidence");
    let canonical_source = format!(
        "{}:{}:{}:{}",
        source_record.requirement().code(),
        source_record.source_kind().slug(),
        source_record.source_schema(),
        source_record.source_id()
    );
    assert!(original.contains(&canonical_source));
    assert_malformed(
        "evidence source with unknown requirement",
        original.replacen(
            &canonical_source,
            &canonical_source.replacen(
                source_record.requirement().code(),
                "code-verification-unknown",
                1,
            ),
            1,
        ),
    );
    assert_malformed(
        "evidence source with wrong artifact kind",
        original.replacen(
            &canonical_source,
            &canonical_source.replacen(source_record.source_kind().slug(), "context-of-use", 1),
            1,
        ),
    );
    assert_malformed(
        "evidence source with wrong schema route",
        original.replacen(
            &canonical_source,
            &canonical_source.replacen(
                source_record.source_schema(),
                "org.frankensim.wrong-evidence-schema.v1",
                1,
            ),
            1,
        ),
    );
    let duplicate_source = format!("{canonical_source}-duplicate");
    assert!(canonical_source.as_str() < duplicate_source.as_str());
    assert_malformed(
        "duplicate evidence requirement source row",
        original.replacen(
            &canonical_source,
            &format!("{canonical_source}\",\"{duplicate_source}"),
            1,
        ),
    );

    let artifact_entry = format!(
        "evidence:{}:artifact:{}",
        source_record.requirement().code(),
        source_record.artifact_hash().to_hex()
    );
    let schema_entry = format!(
        "evidence:{}:schema-admission-receipt:{}",
        source_record.requirement().code(),
        source_record.schema_admission_receipt_hash().to_hex()
    );
    let role_entry = match source_record.authority() {
        EvidenceAuthorityDeclaration::StructuralProcess { receipt_hash } => format!(
            "evidence:{}:role-receipt:{}",
            source_record.requirement().code(),
            receipt_hash.to_hex()
        ),
        EvidenceAuthorityDeclaration::VerifiedNumerics { .. }
        | EvidenceAuthorityDeclaration::ValidatedPhysical { .. } => {
            panic!("first frozen numerical-verification record must be structural")
        }
    };
    assert!(original.contains(&artifact_entry));
    assert!(original.contains(&schema_entry));
    assert!(original.contains(&role_entry));
    let evidence_artifact_identity = source_record.artifact_hash().to_hex();
    let design_set_identity = positive.design_set_identity().to_hex();
    let aggregate_derivation_identity = positive
        .aggregate_qoi_derivation_receipt_identity()
        .to_hex();
    assert_malformed(
        "design-set identity aliased into an evidence slot",
        original.replace(
            design_set_identity.as_str(),
            evidence_artifact_identity.as_str(),
        ),
    );
    assert_malformed(
        "aggregate-QoI derivation identity aliased into an evidence slot",
        original.replace(
            aggregate_derivation_identity.as_str(),
            evidence_artifact_identity.as_str(),
        ),
    );
    let hypothesis_identity = contract.extension().hypothesis_sources()[0]
        .declaration_hash()
        .to_hex();
    for (role, ordinary_identity) in [
        ("design-set", design_set_identity.as_str()),
        (
            "aggregate-QoI derivation-receipt",
            aggregate_derivation_identity.as_str(),
        ),
    ] {
        assert_malformed(
            &format!("unreported hypothesis-source collision in the {role} packet role"),
            original.replace(ordinary_identity, hypothesis_identity.as_str()),
        );
    }
    let hypothesis_collision_entry = format!(
        "evidence:{}:artifact:{}",
        source_record.requirement().code(),
        hypothesis_identity
    );
    assert_malformed(
        "unreported hypothesis-source collision in a positive log",
        original.replacen(&artifact_entry, &hypothesis_collision_entry, 1),
    );

    let refused_reason = refused
        .reasons()
        .first()
        .expect("target-fitted refusal has one reason");
    assert_eq!(refused.reasons().len(), 1);
    let false_collision_reason = format!(
        "hypothesis-source-cannot-satisfy-evidence:{}:artifact",
        source_record.requirement().code()
    );
    assert_malformed(
        "false hypothesis-source collision reason",
        add_refusal_reason(
            refused.log().json_line(),
            refused_reason,
            &false_collision_reason,
        ),
    );
    for role in ["design-set", "aggregate-qoi-derivation-receipt"] {
        assert_malformed(
            &format!("false hypothesis-source collision reason for the {role} packet role"),
            add_refusal_reason(
                refused.log().json_line(),
                refused_reason,
                &format!("hypothesis-source-cannot-satisfy-packet-role:{role}"),
            ),
        );
    }
    for false_reason in [
        "claim-not-present-in-contract".to_owned(),
        format!(
            "source-schema-mismatch:{}",
            source_record.requirement().code()
        ),
        format!(
            "source-kind-mismatch:{}",
            source_record.requirement().code()
        ),
    ] {
        assert_malformed(
            &format!("line-provably false reason {false_reason}"),
            add_refusal_reason(refused.log().json_line(), refused_reason, &false_reason),
        );
    }

    let single_reason = single_refused
        .reasons()
        .first()
        .expect("unit refusal has one reason");
    for false_reason in [
        "missing-prerequisite-receipt:blind-trajectory-prediction:validation-input",
        "stale-prerequisite-receipt:blind-trajectory-prediction",
        "prerequisite-design-set-mismatch:blind-trajectory-prediction",
    ] {
        assert_malformed(
            &format!("prerequisite reason without a compatible frozen-DAG artifact {false_reason}"),
            add_refusal_reason(
                single_refused.log().json_line(),
                single_reason,
                false_reason,
            ),
        );
    }
    let calibrated_reason = calibrated_refused
        .reasons()
        .first()
        .expect("calibrated unit refusal has one reason");
    assert_malformed(
        "stale sole prerequisite without its evaluator-required missing-edge reason",
        add_refusal_reason(
            calibrated_refused.log().json_line(),
            calibrated_reason,
            &format!(
                "stale-prerequisite-receipt:{}",
                EulerClaimKind::NumericalTrajectoryVerification.id()
            ),
        ),
    );
    let missing_numeric_reason = missing_numeric_refused
        .reasons()
        .first()
        .expect("missing numeric axis refusal has one reason");
    assert_eq!(
        missing_numeric_reason,
        &format!("missing-numeric-axis:{}", numeric_axis.as_str())
    );
    assert_malformed(
        "one numeric axis reported both missing and out of domain",
        add_refusal_reason(
            missing_numeric_refused.log().json_line(),
            missing_numeric_reason,
            &format!("out-of-domain-numeric:{}", numeric_axis.as_str()),
        ),
    );
    let missing_categorical_reason = missing_categorical_refused
        .reasons()
        .first()
        .expect("missing categorical axis refusal has one reason");
    assert_eq!(
        missing_categorical_reason,
        &format!("missing-categorical-axis:{}", categorical_axis.as_str())
    );
    assert_malformed(
        "one categorical axis reported both missing and out of domain",
        add_refusal_reason(
            missing_categorical_refused.log().json_line(),
            missing_categorical_reason,
            &format!("out-of-domain-category:{}", categorical_axis.as_str()),
        ),
    );
    let access_reason = access_refused
        .reasons()
        .first()
        .expect("access-class refusal has one reason");
    assert_eq!(
        access_reason,
        "access-class-mismatch:code-verification:expected-not-applicable:observed-calibration"
    );
    assert_malformed(
        "one evidence role reported two observed access classes",
        add_refusal_reason(
            access_refused.log().json_line(),
            access_reason,
            "access-class-mismatch:code-verification:expected-not-applicable:observed-validation",
        ),
    );
    assert_malformed(
        "evidence artifact with zero identity",
        original.replacen(
            &artifact_entry,
            &format!(
                "evidence:{}:artifact:{}",
                source_record.requirement().code(),
                "0".repeat(64)
            ),
            1,
        ),
    );
    assert_malformed(
        "evidence artifact with wrong-width identity",
        original.replacen(
            &artifact_entry,
            &format!(
                "evidence:{}:artifact:{}",
                source_record.requirement().code(),
                "1".repeat(62)
            ),
            1,
        ),
    );
    assert_malformed(
        "evidence row missing schema-admission receipt slot",
        remove_array_entry(original, &schema_entry),
    );
    assert_malformed(
        "complete structural evidence row missing role-receipt slot",
        remove_array_entry(original, &role_entry),
    );
    assert_malformed(
        "refused structural evidence row missing role-receipt slot",
        remove_array_entry(single_refused.log().json_line(), &role_entry),
    );
    let weak_schema_prefix = format!(
        "evidence:{}:schema-admission-receipt:",
        weak_requirement.code()
    );
    for (label, log) in [
        ("demoted", demoted.log().json_line()),
        ("refused", weak_refused.log().json_line()),
    ] {
        assert!(log.contains(&weak_schema_prefix));
        assert!(!log.contains(&format!(
            "evidence:{}:role-receipt:",
            weak_requirement.code()
        )));
        assert_malformed(
            &format!("{label} mismatched numerical authority with injected role receipt"),
            log.replacen(
                &weak_schema_prefix,
                &format!(
                    "evidence:{}:role-receipt:{}\",\"{weak_schema_prefix}",
                    weak_requirement.code(),
                    hash(&format!("injected-{label}-weak-role-receipt")).to_hex()
                ),
                1,
            ),
        );
    }
    let numerical_requirement = golden_requirements(claim)
        .iter()
        .copied()
        .find(|requirement| {
            requirement.authority_class() == EvidenceAuthorityClass::VerifiedNumerics
        })
        .expect("numerical-verification claim has a verified-numerics role");
    let numerical_schema_prefix = format!(
        "evidence:{}:schema-admission-receipt:",
        numerical_requirement.code()
    );
    assert!(demoted.log().json_line().contains(&numerical_schema_prefix));
    assert_malformed(
        "demoted numerical evidence row with injected role-receipt slot",
        demoted.log().json_line().replacen(
            &numerical_schema_prefix,
            &format!(
                "evidence:{}:role-receipt:{}\",\"{numerical_schema_prefix}",
                numerical_requirement.code(),
                hash("injected-numerical-role-receipt").to_hex()
            ),
            1,
        ),
    );
    assert_malformed(
        "one logical hash reused across evidence slots",
        original.replacen(
            &schema_entry,
            &format!(
                "evidence:{}:schema-admission-receipt:{}",
                source_record.requirement().code(),
                source_record.artifact_hash().to_hex()
            ),
            1,
        ),
    );

    assert!(
        refused
            .log()
            .json_line()
            .contains(EULER_PREREQUISITE_RECEIPT_IDENTITY_DOMAIN)
    );
    assert_malformed(
        "prerequisite artifact with wrong owner route",
        refused.log().json_line().replacen(
            EULER_PREREQUISITE_RECEIPT_IDENTITY_DOMAIN,
            "org.frankensim.fs-euler-disc-e2e.wrong-prerequisite-route.v1",
            1,
        ),
    );
    let prerequisite_entries = refused
        .log()
        .json_line()
        .split('"')
        .filter(|value| value.starts_with("prerequisite:"))
        .collect::<Vec<_>>();
    assert!(
        prerequisite_entries.len() >= 2,
        "blind refusal must retain both direct prerequisites"
    );
    let shared_hash = prerequisite_entries[0]
        .rsplit_once(':')
        .expect("prerequisite entry has an identity")
        .1;
    let second_prefix = prerequisite_entries[1]
        .rsplit_once(':')
        .expect("prerequisite entry has an identity")
        .0;
    assert_malformed(
        "one prerequisite identity reused under two claim labels",
        refused.log().json_line().replacen(
            prerequisite_entries[1],
            &format!("{second_prefix}:{shared_hash}"),
            1,
        ),
    );

    let refusal_reason = single_refused
        .reasons()
        .first()
        .expect("single refusal reason");
    assert_malformed(
        "refusal carrying invented hard reason",
        single_refused
            .log()
            .json_line()
            .replace(refusal_reason, "invented-hard-reason"),
    );
    assert_malformed(
        "refusal carrying unknown typed reason operand",
        single_refused
            .log()
            .json_line()
            .replace(refusal_reason, "missing-evidence:not-a-v1-role"),
    );

    assert_malformed(
        "extra field",
        original.replacen("}\n", ",\"unexpected_field\":false}\n", 1),
    );
    assert_malformed(
        "duplicate field",
        original.replacen(
            ",\"packet_identity\":",
            &format!(
                ",\"contract_identity\":\"{}\",\"packet_identity\":",
                contract.identity()
            ),
            1,
        ),
    );
    assert_malformed(
        "missing field",
        original.replacen(
            ",\"artifact_resolution_state\":\"logical-content-identities-only-not-persisted-by-this-crate\"",
            "",
            1,
        ),
    );
    assert_malformed(
        "reordered fields",
        original.replacen(
            ",\"target_fitted\":false,\"applicability_state\":\"campaign-anchor-point-plus-content-bound-design-set\"",
            ",\"applicability_state\":\"campaign-anchor-point-plus-content-bound-design-set\",\"target_fitted\":false",
            1,
        ),
    );
    assert_malformed(
        "wrong JSON type",
        original.replacen(
            &format!("\"packet_source_schema\":\"{EULER_EVIDENCE_PACKET_IDENTITY_DOMAIN}\""),
            "\"packet_source_schema\":1",
            1,
        ),
    );
    assert_malformed(
        "noncanonical integer",
        original.replacen("\"schema_version\":1", "\"schema_version\":01", 1),
    );
    assert_malformed(
        "legacy accuracy-limit wire key",
        original.replacen(
            "\"normalized_accuracy_limit_bits\":",
            "\"accuracy_limit_bits\":",
            1,
        ),
    );
    assert_malformed(
        "noncanonical string escape",
        original.replacen("\"protocol_id\":\"o", "\"protocol_id\":\"\\u006f", 1),
    );
    assert_malformed(
        "short packet identity",
        original.replacen(
            &format!(
                "\"packet_identity\":\"{}\"",
                positive.packet_identity().to_hex()
            ),
            &format!("\"packet_identity\":\"{}\"", "0".repeat(62)),
            1,
        ),
    );
    assert_malformed(
        "zero packet identity",
        original.replacen(
            &format!(
                "\"packet_identity\":\"{}\"",
                positive.packet_identity().to_hex()
            ),
            &format!("\"packet_identity\":\"{}\"", "0".repeat(64)),
            1,
        ),
    );
    assert_malformed(
        "zero design-set identity",
        original.replacen(
            &format!(
                "\"design_set_identity\":\"{}\"",
                positive.design_set_identity().to_hex()
            ),
            &format!("\"design_set_identity\":\"{}\"", "0".repeat(64)),
            1,
        ),
    );
    assert_malformed(
        "zero aggregate-QoI derivation-receipt identity",
        original.replacen(
            &format!(
                "\"aggregate_qoi_derivation_receipt_identity\":\"{}\"",
                positive
                    .aggregate_qoi_derivation_receipt_identity()
                    .to_hex()
            ),
            &format!(
                "\"aggregate_qoi_derivation_receipt_identity\":\"{}\"",
                "0".repeat(64)
            ),
            1,
        ),
    );
    assert_malformed(
        "design-set and aggregate-QoI derivation-receipt identity alias",
        original.replace(
            positive
                .aggregate_qoi_derivation_receipt_identity()
                .to_hex()
                .as_str(),
            positive.design_set_identity().to_hex().as_str(),
        ),
    );
    assert_malformed(
        "unknown claim",
        original.replacen(
            &format!("\"claim\":\"{}\"", positive.claim().id()),
            "\"claim\":\"unknown-claim\"",
            1,
        ),
    );
    assert_malformed(
        "unknown disposition",
        original.replacen(
            "\"expected_disposition\":\"reference-complete-candidate-unreadmitted\"",
            "\"expected_disposition\":\"unknown\"",
            1,
        ),
    );
    assert_malformed(
        "authority mismatch",
        original.replacen(
            "\"authority_state\":\"unreadmitted-reference-candidate-only\"",
            "\"authority_state\":\"no-candidate-authority\"",
            1,
        ),
    );
    assert_malformed(
        "missing expected-disposition mismatch reason",
        original.replacen(
            "\"expected_disposition\":\"reference-complete-candidate-unreadmitted\"",
            "\"expected_disposition\":\"refused\"",
            1,
        ),
    );
    assert_malformed(
        "candidate with nonpositive reported outcome",
        original.replacen(
            "\"reported_scientific_disposition\":\"positive\"",
            "\"reported_scientific_disposition\":\"negative\"",
            1,
        ),
    );
    assert_malformed(
        "unbound no-claim state",
        original.replacen(
            "\"no_claim_state\":\"accepted\"",
            "\"no_claim_state\":\"not-accepted\"",
            1,
        ),
    );
    assert_malformed(
        "unbound target-fitting state",
        multi_refused.log().json_line().replacen(
            "\"target_fitted\":true",
            "\"target_fitted\":false",
            1,
        ),
    );
    assert_malformed(
        "candidate carrying a policy reason",
        original.replacen(
            "\"first_divergence\":null,\"reasons\":[]",
            "\"first_divergence\":\"invented-policy-reason\",\"reasons\":[\"invented-policy-reason\"]",
            1,
        ),
    );
    assert_malformed(
        "noncanonical expected-disposition reason",
        original.replacen(
            "\"first_divergence\":null,\"reasons\":[]",
            "\"first_divergence\":\"expected-disposition-mismatch:bogus\",\"reasons\":[\"expected-disposition-mismatch:bogus\"]",
            1,
        ),
    );
    let demotion_reason = demoted.reasons().first().expect("demotion reason");
    assert!(demotion_reason.starts_with("weak-"));
    assert_malformed(
        "demotion carrying a hard reason",
        demoted
            .log()
            .json_line()
            .replace(demotion_reason, &format!("hard-{demotion_reason}")),
    );
    let refusal_reason = single_refused
        .reasons()
        .first()
        .expect("single refusal reason");
    assert!(!refusal_reason.starts_with("weak-"));
    assert_malformed(
        "refusal carrying only a weak reason",
        single_refused
            .log()
            .json_line()
            .replace(refusal_reason, &format!("weak-{refusal_reason}")),
    );
    let first_divergence = multi_refused
        .log()
        .json_line()
        .split("\"first_divergence\":\"")
        .nth(1)
        .and_then(|tail| tail.split('"').next())
        .expect("refused log carries first divergence");
    assert_malformed(
        "unbound first divergence",
        multi_refused.log().json_line().replacen(
            &format!("\"first_divergence\":\"{first_divergence}\""),
            "\"first_divergence\":\"not-a-retained-reason\"",
            1,
        ),
    );
    let swap_first_two = |line: &str, field: &str| {
        let marker = format!("\"{field}\":[\"");
        let first_start = line.find(&marker).expect("array field") + marker.len();
        let first_end = first_start
            + line[first_start..]
                .find("\",\"")
                .expect("at least two array values");
        let second_start = first_end + 3;
        let second_end = second_start
            + line[second_start..]
                .find('"')
                .expect("second array value terminator");
        let first = &line[first_start..first_end];
        let second = &line[second_start..second_end];
        let mut swapped = line.to_owned();
        swapped.replace_range(first_start..second_end, &format!("{second}\",\"{first}"));
        swapped
    };
    assert_malformed(
        "unsorted reasons",
        swap_first_two(multi_refused.log().json_line(), "reasons"),
    );
    assert_malformed(
        "unsorted relative artifacts",
        swap_first_two(original, "relative_artifacts"),
    );
    let mut without_lf = original.to_owned();
    assert_eq!(without_lf.pop(), Some('\n'));
    assert_malformed("missing terminal LF", without_lf);
    assert_malformed(
        "embedded carriage return",
        original.replacen(",\"claim\":", "\r,\"claim\":", 1),
    );
    assert_malformed(
        "CRLF terminator",
        format!(
            "{}\r\n",
            original.strip_suffix('\n').expect("canonical terminal LF")
        ),
    );
    assert_ne!(positive.log().identity(), refused.log().identity());
    assert_eq!(
        positive.log().identity(),
        fs_blake3::hash_domain(
            CLAIM_POLICY_ASSESSMENT_LOG_DOMAIN,
            positive.log().json_line().as_bytes()
        )
    );
}

#[test]
fn assessment_log_preserves_detection_order_separately_from_canonical_reason_order() {
    let contract = build_frozen_contract().expect("frozen contract");
    let claim = EulerClaimKind::QualitativeEffectDirection;
    let assessment = assess_with_frozen_prerequisites(
        &contract,
        packet(
            &contract,
            claim,
            Vec::new(),
            true,
            ReportedScientificDisposition::Positive,
            AssessmentDisposition::Refused,
        ),
    );
    assert_eq!(assessment.disposition(), AssessmentDisposition::Refused);
    assert!(assessment.reasons()[0].starts_with("missing-evidence:"));
    assert!(
        assessment.log().json_line().contains(
            "\"first_divergence\":\"protected-target-fitting-invalidates-emergent-claim\""
        )
    );
    let log = assessment.log().json_line();
    let mut reason_cursor = log
        .find("\"reasons\":[")
        .expect("assessment log carries canonical reasons");
    for reason in assessment.reasons() {
        let relative = log[reason_cursor..]
            .find(reason)
            .unwrap_or_else(|| panic!("assessment log omitted canonical reason {reason}"));
        reason_cursor += relative + reason.len();
    }
}

#[test]
fn euler_identity_versions_and_domains_fail_closed() {
    assert_eq!(EULER_CONTRACT_SCHEMA_VERSION, 1);
    assert_eq!(EULER_CLAIM_POLICY_SCHEMA_VERSION, 1);
    assert_eq!(EULER_OWNER_MATRIX_SCHEMA_VERSION, 1);
    assert_eq!(EULER_PROTOCOL_SCHEMA_VERSION, 1);
    fs_euler_disc_e2e::protocol::protocol_migration_policy(1).expect("current protocol version");
    for unsupported in [0, 2, u32::MAX] {
        assert!(
            fs_euler_disc_e2e::protocol::protocol_migration_policy(unsupported).is_err(),
            "protocol version {unsupported} must refuse"
        );
    }
    let identity_domains = [
        HYPOTHESIS_SOURCE_DECLARATION_DOMAIN,
        EULER_CLAIM_GRAPH_IDENTITY_DOMAIN,
        EULER_OWNER_MATRIX_IDENTITY_DOMAIN,
        EULER_CONTRACT_IDENTITY_DOMAIN,
        EULER_EVIDENCE_PACKET_IDENTITY_DOMAIN,
        EULER_PREREQUISITE_RECEIPT_IDENTITY_DOMAIN,
        EULER_ASSESSMENT_IDENTITY_DOMAIN,
        fs_euler_disc_e2e::CONTRACT_CHECK_RECEIPT_DOMAIN,
        CLAIM_POLICY_ASSESSMENT_LOG_DOMAIN,
    ];
    assert_eq!(
        identity_domains.into_iter().collect::<BTreeSet<_>>().len(),
        identity_domains.len()
    );
    assert!(
        identity_domains
            .iter()
            .all(|domain| domain.strip_suffix(".v1").is_some())
    );
    let payload = b"same semantic bytes";
    assert_eq!(
        identity_domains
            .iter()
            .map(|domain| fs_blake3::hash_domain(domain, payload))
            .collect::<BTreeSet<_>>()
            .len(),
        identity_domains.len()
    );
    let routed_schemas =
        [fs_euler_disc_e2e::contract::EULER_AGGREGATE_QOI_DERIVATION_RECEIPT_SCHEMA];
    assert!(
        routed_schemas
            .iter()
            .all(|schema| schema.strip_suffix(".v1").is_some())
    );
    assert!(
        routed_schemas
            .iter()
            .all(|schema| !identity_domains.contains(schema))
    );
    for version in [0, 2, u32::MAX] {
        assert!(EulerScientificContract::migration_policy(version).is_err());
        assert!(OwnerMatrix::migration_policy(version).is_err());
    }
}

#[test]
fn g0_schema_migration_and_manifest_direction_are_explicit() {
    assert!(EulerScientificContract::migration_policy(1).is_ok());
    assert!(OwnerMatrix::migration_policy(1).is_ok());
    for version in [0, 2, u32::MAX] {
        assert!(EulerScientificContract::migration_policy(version).is_err());
    }

    let root_manifest = include_str!("../../../Cargo.toml");
    let crate_manifest = include_str!("../Cargo.toml");
    assert!(root_manifest.contains("\"crates/fs-euler-disc-e2e\""));
    for dependency in ["fs-blake3", "fs-evidence", "fs-govern", "fs-ir"] {
        assert!(
            crate_manifest.contains(&format!("{dependency} = {{ path = \"../{dependency}\" }}"))
        );
    }
    for forbidden in [
        "fs-contact",
        "fs-solid",
        "fs-lbm",
        "fs-opt",
        "fs-render",
        "fs-session",
    ] {
        assert!(!crate_manifest.contains(&format!("{forbidden} =")));
    }
}

#[test]
fn g0_validity_domain_axis_cardinality_has_exact_maximum_boundaries() {
    let contract = build_frozen_contract().expect("frozen contract");
    let exact = (0..MAX_VALIDITY_DOMAIN_AXES)
        .fold(ValidityDomain::unconstrained(), |regime, index| {
            regime.with(format!("axis-{index:02}"), 0.0, 1.0)
        });
    validated_record_with_regime(&contract, exact, "dataset-axis-max", "axis-max")
        .expect("the exact validity-axis maximum must be admitted");

    let maximum_plus_one = (0..=MAX_VALIDITY_DOMAIN_AXES)
        .fold(ValidityDomain::unconstrained(), |regime, index| {
            regime.with(format!("axis-{index:02}"), 0.0, 1.0)
        });
    let error = validated_record_with_regime(
        &contract,
        maximum_plus_one,
        "dataset-axis-max-plus-one",
        "axis-max-plus-one",
    )
    .expect_err("maximum-plus-one validity axes must refuse before serialization");
    assert_eq!(error.code(), "EulerProtocolValidityDomainCardinality");
}

#[test]
fn g0_validity_domain_canonical_bytes_are_preflighted_at_exact_boundary() {
    let contract = build_frozen_contract().expect("frozen contract");
    let dataset = "dataset-validity-byte-boundary";
    let exact = validity_domain_with_exact_canonical_bytes(
        MAX_VALIDITY_DOMAIN_AXES,
        MAX_VALIDITY_DOMAIN_CANONICAL_BYTES,
    );
    let exact_color = Color::Validated {
        regime: exact.clone(),
        dataset: dataset.to_owned(),
    };
    assert_eq!(
        exact_color.canonical_bytes().len(),
        2 + size_of::<u64>() + dataset.len() + MAX_VALIDITY_DOMAIN_CANONICAL_BYTES,
        "the local preflight must match the shared Color v2 wire exactly"
    );
    validated_record_with_regime(&contract, exact, dataset, "byte-max")
        .expect("the exact canonical validity-domain byte maximum must be admitted");

    let maximum_plus_one = validity_domain_with_exact_canonical_bytes(
        MAX_VALIDITY_DOMAIN_AXES,
        MAX_VALIDITY_DOMAIN_CANONICAL_BYTES + 1,
    );
    let error =
        validated_record_with_regime(&contract, maximum_plus_one, dataset, "byte-max-plus-one")
            .expect_err(
                "maximum-plus-one canonical validity bytes must refuse before serialization",
            );
    assert_eq!(error.code(), "EulerProtocolValidityDomainTooLarge");
}

#[test]
fn g0_validated_regime_requires_every_context_and_regime_axis_at_inclusive_bounds() {
    let contract = build_frozen_contract().expect("frozen contract");
    let claim = EulerClaimKind::BlindTrajectoryPrediction;
    let extra_axis =
        fs_evidence::vv::AxisId::try_new("surface-film-thickness").expect("extra validity axis");
    let regime = covering_regime(&contract).with(extra_axis.as_str(), 0.25, 0.75);
    let records_with_regime = |label: &str| {
        records(&contract, claim)
            .into_iter()
            .map(|record| {
                if record.requirement() == EvidenceRequirement::PhysicalValidation {
                    validated_record_with_regime(
                        &contract,
                        regime.clone(),
                        &format!("dataset-{label}"),
                        label,
                    )
                    .expect("validated record")
                } else {
                    record
                }
            })
            .collect::<Vec<_>>()
    };
    let point_with_extra = |value: f64| {
        let mut numeric = contract
            .context()
            .applicability()
            .numeric()
            .iter()
            .map(|(axis, domain)| {
                let (lo, hi) = domain.bounds();
                (axis.clone(), lo + (hi - lo) * 0.5)
            })
            .collect::<Vec<_>>();
        numeric.push((extra_axis.clone(), value));
        let categorical = contract
            .context()
            .applicability()
            .categorical()
            .iter()
            .map(|(axis, domain)| {
                (
                    axis.clone(),
                    domain
                        .allowed()
                        .iter()
                        .next()
                        .expect("allowed value")
                        .clone(),
                )
            })
            .collect();
        ApplicabilityPoint::try_new(numeric, categorical).expect("finite point")
    };

    let missing = assess_with_frozen_prerequisites(
        &contract,
        packet_at_point(
            &contract,
            claim,
            point(&contract),
            records_with_regime("missing-extra-axis"),
            false,
            ReportedScientificDisposition::Positive,
            AssessmentDisposition::DemotedCandidate,
        ),
    );
    assert_eq!(
        missing.disposition(),
        AssessmentDisposition::DemotedCandidate
    );
    assert!(
        missing
            .reasons()
            .contains(&"weak-validity-domain:physical-validation:does-not-cover-case".to_owned())
    );

    let outside = assess_with_frozen_prerequisites(
        &contract,
        packet_at_point(
            &contract,
            claim,
            point_with_extra(0.750_000_1),
            records_with_regime("outside-extra-axis"),
            false,
            ReportedScientificDisposition::Positive,
            AssessmentDisposition::DemotedCandidate,
        ),
    );
    assert_eq!(
        outside.disposition(),
        AssessmentDisposition::DemotedCandidate
    );
    assert!(
        outside
            .reasons()
            .contains(&"weak-validity-domain:physical-validation:does-not-cover-case".to_owned())
    );

    for (label, boundary) in [("lower-boundary", 0.25), ("upper-boundary", 0.75)] {
        let assessment = assess_with_frozen_prerequisites(
            &contract,
            packet_at_point(
                &contract,
                claim,
                point_with_extra(boundary),
                records_with_regime(label),
                false,
                ReportedScientificDisposition::Positive,
                AssessmentDisposition::ReferenceCompleteCandidate,
            ),
        );
        assert_eq!(
            assessment.disposition(),
            AssessmentDisposition::ReferenceCompleteCandidate,
            "inclusive {label} must be covered: {:?}",
            assessment.reasons()
        );
    }
}

#[test]
fn g0_evidence_constructors_and_packet_bounds_refuse_hostile_inputs() {
    use fs_evidence::vv::ArtifactKind;

    let contract = build_frozen_contract().expect("frozen contract");
    let admitted = admit_frozen_contract(contract.clone()).expect("structural admission");
    let claim = EulerClaimKind::NumericalTrajectoryVerification;
    let complete = records(&contract, claim);
    let structural = complete
        .iter()
        .find(|record| record.requirement() == EvidenceRequirement::CodeVerification)
        .expect("structural record");
    let numerical = complete
        .iter()
        .find(|record| record.requirement() == EvidenceRequirement::SolutionVerification)
        .expect("numerical record");
    let zero = ContentHash([0; 32]);

    let rebuild = |record: &EvidenceRecord,
                   authority: EvidenceAuthorityDeclaration,
                   artifact_hash: ContentHash,
                   source_kind: ArtifactKind,
                   schema_receipt: ContentHash,
                   qois: Vec<fs_evidence::vv::QoiId>| {
        EvidenceRecord::try_new(
            record.contract_identity(),
            record.claim(),
            record.requirement(),
            qois,
            authority,
            artifact_hash,
            record.source_id(),
            record.source_schema(),
            source_kind,
            schema_receipt,
            record.access_class(),
            record.independent(),
        )
    };

    let error = rebuild(
        structural,
        structural.authority().clone(),
        zero,
        structural.source_kind(),
        structural.schema_admission_receipt_hash(),
        structural.qois().to_vec(),
    )
    .expect_err("zero artifact identity");
    assert_eq!(error.code(), "EulerProtocolZeroIdentity");

    let error = rebuild(
        structural,
        structural.authority().clone(),
        structural.artifact_hash(),
        structural.source_kind(),
        zero,
        structural.qois().to_vec(),
    )
    .expect_err("zero schema-admission receipt identity");
    assert_eq!(error.code(), "EulerProtocolZeroIdentity");

    let error = rebuild(
        structural,
        EvidenceAuthorityDeclaration::StructuralProcess { receipt_hash: zero },
        structural.artifact_hash(),
        structural.source_kind(),
        structural.schema_admission_receipt_hash(),
        structural.qois().to_vec(),
    )
    .expect_err("zero process receipt identity");
    assert_eq!(error.code(), "EulerProtocolZeroIdentity");

    let error = rebuild(
        structural,
        structural.authority().clone(),
        structural.artifact_hash(),
        ArtifactKind::ContextOfUse,
        structural.schema_admission_receipt_hash(),
        structural.qois().to_vec(),
    )
    .expect_err("wrong generic container kind");
    assert_eq!(error.code(), "EulerProtocolSourceKindMismatch");

    let error = rebuild(
        numerical,
        EvidenceAuthorityDeclaration::VerifiedNumerics {
            color: Color::Verified {
                lo: f64::NAN,
                hi: 0.0,
            },
        },
        numerical.artifact_hash(),
        numerical.source_kind(),
        numerical.schema_admission_receipt_hash(),
        numerical.qois().to_vec(),
    )
    .expect_err("malformed numerical color");
    assert_eq!(error.code(), "EulerProtocolMalformedColor");

    for (wall_ms, memory_bytes, accuracy) in [
        (0, 1, 0.0),
        (1, 0, 0.0),
        (1, 1, f64::NAN),
        (1, 1, f64::INFINITY),
        (1, 1, -1.0),
    ] {
        let error =
            ProtocolBudget::try_new(wall_ms, memory_bytes, accuracy).expect_err("invalid budget");
        assert_eq!(error.code(), "EulerProtocolInvalidBudget");
    }
    assert_eq!(
        ProtocolBudget::try_new(1, 1, -0.0)
            .expect("signed zero normalizes")
            .normalized_accuracy_limit()
            .to_bits(),
        0.0_f64.to_bits()
    );

    let error = ClaimEvidencePacket::try_new(
        contract.identity(),
        "malformed-direct-seed-variant",
        hash("malformed-seed-design-set"),
        hash("malformed-seed-aggregate-qoi-derivation-receipt"),
        claim,
        point(&contract),
        complete.clone(),
        true,
        false,
        ReportedScientificDisposition::Positive,
        AssessmentDisposition::Refused,
        claim_units(&contract, claim),
        ProtocolSeed::NotApplicable {
            reason: " not-a-canonical-machine-reason ".to_owned(),
        },
        ProtocolBudget::try_new(1, 1, 0.0).expect("budget"),
    )
    .expect_err("public enum construction cannot bypass seed validation");
    assert_eq!(error.code(), "EulerProtocolInvalidIdentity");

    let too_many = vec![structural.clone(); fs_euler_disc_e2e::protocol::MAX_EVIDENCE_RECORDS + 1];
    let error = ClaimEvidencePacket::try_new(
        contract.identity(),
        "too-many-evidence-rows",
        hash("too-many-rows-design-set"),
        hash("too-many-rows-aggregate-qoi-derivation-receipt"),
        claim,
        point(&contract),
        too_many,
        true,
        false,
        ReportedScientificDisposition::Positive,
        AssessmentDisposition::Refused,
        claim_units(&contract, claim),
        ProtocolSeed::not_applicable("bounded-cardinality-test").expect("seed"),
        ProtocolBudget::try_new(1, 1, 0.0).expect("budget"),
    )
    .expect_err("maximum plus one evidence rows");
    assert_eq!(error.code(), "EulerProtocolEvidenceCardinality");

    let stale_qoi = fs_evidence::vv::QoiId::try_new("deliberately-wrong-qoi").expect("qoi");
    let qoi_mismatch = complete
        .into_iter()
        .map(|record| {
            if record.requirement() != EvidenceRequirement::CodeVerification {
                return record;
            }
            rebuild(
                &record,
                record.authority().clone(),
                record.artifact_hash(),
                record.source_kind(),
                record.schema_admission_receipt_hash(),
                vec![stale_qoi.clone()],
            )
            .expect("locally well-formed stale QoI declaration")
        })
        .collect();
    let assessment = packet(
        &contract,
        claim,
        qoi_mismatch,
        false,
        ReportedScientificDisposition::Positive,
        AssessmentDisposition::Refused,
    )
    .assess(&admitted, &[])
    .expect("QoI mismatch assessment");
    assert_eq!(assessment.disposition(), AssessmentDisposition::Refused);
    assert!(
        assessment
            .reasons()
            .contains(&"qoi-binding-mismatch:code-verification".to_owned())
    );
}

#[test]
fn g0_public_evidence_record_bounds_color_identity_refusals() {
    const VERY_LARGE_IDENTITY_BYTES: usize = fs_evidence::MAX_COLOR_IDENTITY_BYTES * 4_096;

    let contract = build_frozen_contract().expect("frozen contract");
    let physical_records = records(&contract, EulerClaimKind::BlindTrajectoryPrediction);
    let physical = physical_records
        .iter()
        .find(|record| record.requirement() == EvidenceRequirement::PhysicalValidation)
        .expect("physical-validation record");
    let numerical_records = records(&contract, EulerClaimKind::NumericalTrajectoryVerification);
    let numerical = numerical_records
        .iter()
        .find(|record| record.requirement() == EvidenceRequirement::SolutionVerification)
        .expect("solution-verification record");

    // Exercise the public record constructor, not only the private authority
    // validator: removing the delegation in `EvidenceRecord::try_new` must make
    // these hostile-boundary cases fail.
    let rebuild = |record: &EvidenceRecord, authority: EvidenceAuthorityDeclaration| {
        EvidenceRecord::try_new(
            record.contract_identity(),
            record.claim(),
            record.requirement(),
            record.qois().to_vec(),
            authority,
            record.artifact_hash(),
            record.source_id(),
            record.source_schema(),
            record.source_kind(),
            record.schema_admission_receipt_hash(),
            record.access_class(),
            record.independent(),
        )
    };
    let validated =
        |dataset: String, regime: ValidityDomain| EvidenceAuthorityDeclaration::ValidatedPhysical {
            color: Color::Validated { regime, dataset },
        };
    let estimated = |estimator: String| EvidenceAuthorityDeclaration::VerifiedNumerics {
        color: Color::Estimated {
            estimator,
            dispersion: 0.125,
        },
    };
    let assert_bounded_refusal =
        |plus_one: ContractError, very_large: ContractError, expected_detail: &str| {
            assert_eq!(plus_one.code(), "EulerProtocolMalformedColor");
            assert_eq!(very_large.code(), "EulerProtocolMalformedColor");
            assert_eq!(plus_one.detail(), expected_detail);
            assert_eq!(very_large.detail(), plus_one.detail());
            assert!(plus_one.detail().len() < fs_evidence::MAX_COLOR_IDENTITY_BYTES);
        };

    rebuild(
        physical,
        validated(
            "d".repeat(fs_evidence::MAX_COLOR_IDENTITY_BYTES),
            covering_regime(&contract),
        ),
    )
    .expect("a dataset identity exactly at the shared byte limit must remain admissible");
    let dataset_plus_one = rebuild(
        physical,
        validated(
            "d".repeat(fs_evidence::MAX_COLOR_IDENTITY_BYTES + 1),
            covering_regime(&contract),
        ),
    )
    .expect_err("a dataset identity one byte over the limit must refuse publicly");
    let dataset_very_large = rebuild(
        physical,
        validated(
            "d".repeat(VERY_LARGE_IDENTITY_BYTES),
            covering_regime(&contract),
        ),
    )
    .expect_err("a very large dataset identity must refuse publicly with bounded detail");
    assert_bounded_refusal(
        dataset_plus_one,
        dataset_very_large,
        "validated-color dataset identity exceeds the v1 byte limit of 256",
    );

    rebuild(
        numerical,
        estimated("e".repeat(fs_evidence::MAX_COLOR_IDENTITY_BYTES)),
    )
    .expect("an estimator identity exactly at the shared byte limit must remain admissible");
    let estimator_plus_one = rebuild(
        numerical,
        estimated("e".repeat(fs_evidence::MAX_COLOR_IDENTITY_BYTES + 1)),
    )
    .expect_err("an estimator identity one byte over the limit must refuse publicly");
    let estimator_very_large = rebuild(numerical, estimated("e".repeat(VERY_LARGE_IDENTITY_BYTES)))
        .expect_err("a very large estimator identity must refuse publicly with bounded detail");
    assert_bounded_refusal(
        estimator_plus_one,
        estimator_very_large,
        "estimated-color estimator identity exceeds the v1 byte limit of 256",
    );

    let regime =
        |axis_bytes: usize| ValidityDomain::unconstrained().with("a".repeat(axis_bytes), 0.0, 1.0);
    rebuild(
        physical,
        validated(
            "dataset-for-exact-axis".to_owned(),
            regime(fs_evidence::MAX_COLOR_IDENTITY_BYTES),
        ),
    )
    .expect("a regime axis identity exactly at the shared byte limit must remain admissible");
    let axis_plus_one = rebuild(
        physical,
        validated(
            "dataset-for-axis-plus-one".to_owned(),
            regime(fs_evidence::MAX_COLOR_IDENTITY_BYTES + 1),
        ),
    )
    .expect_err("a regime axis identity one byte over the limit must refuse publicly");
    let axis_very_large = rebuild(
        physical,
        validated(
            "dataset-for-large-axis".to_owned(),
            regime(VERY_LARGE_IDENTITY_BYTES),
        ),
    )
    .expect_err("a very large regime axis must refuse publicly with bounded detail");
    assert_bounded_refusal(
        axis_plus_one,
        axis_very_large,
        "validated-color regime axis identity exceeds the v1 byte limit of 256",
    );
}

#[test]
fn g0_public_claim_packet_bounds_direct_seed_reason_refusals() {
    const VERY_LARGE_REASON_BYTES: usize = MAX_PROTOCOL_ID_BYTES * 4_096;

    let contract = build_frozen_contract().expect("frozen contract");
    let claim = EulerClaimKind::NumericalTrajectoryVerification;
    // Construct the public enum variant directly so this test proves that the
    // packet constructor itself retains the validation boundary.
    let construct = |case_id: &str, reason: String| {
        ClaimEvidencePacket::try_new(
            contract.identity(),
            case_id,
            hash("direct-seed-design-set"),
            hash("direct-seed-aggregate-qoi-derivation-receipt"),
            claim,
            point(&contract),
            records(&contract, claim),
            true,
            false,
            ReportedScientificDisposition::Positive,
            AssessmentDisposition::ReferenceCompleteCandidate,
            claim_units(&contract, claim),
            ProtocolSeed::NotApplicable { reason },
            ProtocolBudget::try_new(1, 1, 0.0).expect("budget"),
        )
    };

    construct(
        "direct-seed-reason-exact-limit",
        "s".repeat(MAX_PROTOCOL_ID_BYTES),
    )
    .expect("a direct seed reason exactly at the byte limit must remain admissible");
    let plus_one = construct(
        "direct-seed-reason-plus-one",
        "s".repeat(MAX_PROTOCOL_ID_BYTES + 1),
    )
    .expect_err("a direct seed reason one byte over the limit must refuse publicly");
    let very_large = construct(
        "direct-seed-reason-very-large",
        "s".repeat(VERY_LARGE_REASON_BYTES),
    )
    .expect_err("a very large direct seed reason must refuse publicly with bounded detail");

    assert_eq!(plus_one.code(), "EulerProtocolInvalidIdentity");
    assert_eq!(very_large.code(), "EulerProtocolInvalidIdentity");
    assert_eq!(
        plus_one.detail(),
        "packet.seed.reason must be a bounded canonical machine identity"
    );
    assert_eq!(very_large.detail(), plus_one.detail());
    assert!(plus_one.detail().len() < MAX_PROTOCOL_ID_BYTES);
}

#[test]
fn g0_no_claim_constructor_and_decoder_share_the_exact_72_row_boundary() {
    let base = build_frozen_contract().expect("frozen contract");
    let make_boundary = |count: usize| {
        let mut rows = base.no_claims().entries().to_vec();
        rows.extend(
            (rows.len()..count)
                .map(|index| format!("Additional narrowing no-claim boundary {index:03}.")),
        );
        let refs = rows.iter().map(String::as_str).collect::<Vec<_>>();
        NoClaimBoundary::new(&refs).expect("bounded no-claim set")
    };

    let exact = EulerScientificContract::try_new(
        base.context().clone(),
        base.extension().clone(),
        base.claim_graph().clone(),
        make_boundary(MAX_EULER_NO_CLAIMS),
        base.owner_matrix().clone(),
    )
    .expect("the exact no-claim maximum must publish an identity");
    let exact_bytes = exact.canonical_bytes().expect("exact-boundary bytes");
    assert_eq!(
        EulerScientificContract::from_canonical_bytes(&exact_bytes)
            .expect("exact-boundary fixed point"),
        exact
    );

    let error = EulerScientificContract::try_new(
        base.context().clone(),
        base.extension().clone(),
        base.claim_graph().clone(),
        make_boundary(MAX_EULER_NO_CLAIMS + 1),
        base.owner_matrix().clone(),
    )
    .expect_err("maximum-plus-one no-claim rows must refuse before identity publication");
    assert_eq!(error.code(), "EulerContractNoClaimCardinality");
}

#[test]
fn g0_empty_boundary_and_malformed_local_rows_refuse() {
    let base = build_frozen_contract().expect("frozen contract");
    let error = EulerContextExtension::try_new(
        vec![],
        base.extension().apparatus_population(),
        base.extension().environment_population(),
        base.extension().observation_frame(),
        base.extension().decision_alternatives().to_vec(),
        base.extension().risks().to_vec(),
        base.extension().hypothesis_sources().to_vec(),
    )
    .expect_err("empty users");
    assert_eq!(error.code(), "EulerContractCardinality");

    let error = ScientificRisk::try_new(
        "bad-severity",
        "Bad severity must refuse.",
        0,
        vec![EulerClaimKind::Ranking],
        base.extension().decision_alternatives()[0].clone(),
    )
    .expect_err("severity zero");
    assert_eq!(error.code(), "EulerContractInvalidRisk");

    let error = OwnerMatrix::try_new(vec![]).expect_err("empty owner matrix");
    assert_eq!(error.code(), "EulerContractOwnerMatrixIncomplete");

    let mut campaign: CampaignClaim = base
        .claim_graph()
        .claim(EulerClaimKind::MechanismAttribution)
        .expect("mechanism")
        .campaign()
        .clone();
    campaign.qois =
        vec![fs_evidence::vv::QoiId::try_new("event-time-error").expect("wrong but valid qoi")];
    let error = EulerClaimSpec::try_new(
        EulerClaimKind::MechanismAttribution,
        campaign,
        golden_requirements(EulerClaimKind::MechanismAttribution).to_vec(),
    )
    .expect_err("stop time alone cannot identify mechanism");
    assert_eq!(error.code(), "EulerContractClaimQoiPolicyMismatch");

    let mut duplicate_gap_claims = base
        .claim_graph()
        .claims()
        .values()
        .cloned()
        .collect::<Vec<_>>();
    let shared_gap = EvidenceGapId::try_new("globally-duplicated-gap").expect("gap id");
    for claim in duplicate_gap_claims.iter_mut().take(2) {
        let mut campaign = claim.campaign().clone();
        campaign.evidence_gaps[0].id = shared_gap.clone();
        *claim = EulerClaimSpec::try_new(claim.kind(), campaign, claim.requirements().to_vec())
            .expect("claim-local gap remains valid");
    }
    let error =
        EulerClaimGraph::try_new(duplicate_gap_claims, vec![]).expect_err("global duplicate gap");
    assert_eq!(error.code(), "EulerContractDuplicate");
}
