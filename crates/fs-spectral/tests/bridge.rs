//! G0/G3/G5 statement-schema tests for the Maslov--Krein--Evans bridge.

#![allow(clippy::too_many_lines, clippy::wildcard_imports)]

use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_spectral::bridge::*;

macro_rules! id {
    ($ty:ident, $byte:expr) => {
        $ty::from_bytes([$byte; 32])
    };
}

fn with_cx<R>(cancelled: bool, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    if cancelled {
        gate.request();
    }
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x4D_4B_45,
                kernel_id: 11,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn convention(orientation: CountOrientationV1) -> SignedCountConventionV1 {
    SignedCountConventionV1 {
        orientation,
        endpoints: EndpointConventionV1 {
            left: EndpointWeightV1::Excluded,
            right: EndpointWeightV1::Excluded,
        },
    }
}

fn transform(byte: u8) -> SignedCountTransformV1 {
    SignedCountTransformV1::new(id!(BridgeCorrespondenceMapIdV1, byte), 1, 0)
        .expect("positive sign is valid")
}

const ALL_HYPOTHESES: [BridgeHypothesisKindV1; 17] = [
    BridgeHypothesisKindV1::SymplecticPath,
    BridgeHypothesisKindV1::LagrangianFredholmPair,
    BridgeHypothesisKindV1::CrossingFormsControlled,
    BridgeHypothesisKindV1::EndpointConventionResolved,
    BridgeHypothesisKindV1::PontryaginPathContinuous,
    BridgeHypothesisKindV1::KreinFormNondegenerate,
    BridgeHypothesisKindV1::NeutralSignaturePolicyClosed,
    BridgeHypothesisKindV1::AnalyticFredholmFamily,
    BridgeHypothesisKindV1::EssentialSpectrumExcluded,
    BridgeHypothesisKindV1::ContourAdmissible,
    BridgeHypothesisKindV1::EvansNormalizationFixed,
    BridgeHypothesisKindV1::PeriodicMonodromyCorrespondence,
    BridgeHypothesisKindV1::SpatialDichotomy,
    BridgeHypothesisKindV1::ParameterDirectionAligned,
    BridgeHypothesisKindV1::MultiplicityPreserved,
    BridgeHypothesisKindV1::CorrespondenceMapsExact,
    BridgeHypothesisKindV1::MachineInterpretationExact,
];

const ALL_FALSIFIERS: [BridgeFalsifierKindV1; 10] = [
    BridgeFalsifierKindV1::DegenerateCrossing,
    BridgeFalsifierKindV1::TangentialCrossing,
    BridgeFalsifierKindV1::EndpointCrossing,
    BridgeFalsifierKindV1::NeutralKreinSignature,
    BridgeFalsifierKindV1::NonFredholmFamily,
    BridgeFalsifierKindV1::ContourContact,
    BridgeFalsifierKindV1::EssentialSpectrumContact,
    BridgeFalsifierKindV1::OrientationMutation,
    BridgeFalsifierKindV1::MultiplicityMutation,
    BridgeFalsifierKindV1::StatementMutation,
];

fn hypotheses() -> Vec<BridgeHypothesisV1> {
    ALL_HYPOTHESES
        .iter()
        .enumerate()
        .map(|(index, &kind)| BridgeHypothesisV1 {
            kind,
            state: BridgeHypothesisStateV1::WitnessReferenced {
                witness: BridgeHypothesisWitnessIdV1::from_bytes([(80 + index) as u8; 32]),
            },
        })
        .collect()
}

fn falsifiers() -> Vec<BridgeFalsifierV1> {
    ALL_FALSIFIERS
        .iter()
        .enumerate()
        .map(|(index, &kind)| BridgeFalsifierV1 {
            kind,
            artifact: BridgeFalsifierArtifactIdV1::from_bytes([(120 + index) as u8; 32]),
            response: BridgeFalsifierResponseV1::RefuseNode,
        })
        .collect()
}

fn node(
    byte: u8,
    scope: BridgeTheoremScopeV1,
    conclusion: BridgeConclusionV1,
    proof: BridgeProofStateV1,
) -> BridgeTheoremNodeV1 {
    BridgeTheoremNodeV1 {
        id: id!(BridgeTheoremNodeIdV1, byte),
        scope,
        domain: BridgeDomainIdV1::from_bytes([byte.wrapping_add(20); 32]),
        operator_family: BridgeOperatorFamilyIdV1::from_bytes([byte.wrapping_add(30); 32]),
        parameterization: BridgeParameterizationIdV1::from_bytes([byte.wrapping_add(40); 32]),
        conclusion,
        hypotheses: hypotheses(),
        falsifiers: falsifiers(),
        proof,
    }
}

fn statement_only(byte: u8) -> BridgeProofStateV1 {
    BridgeProofStateV1::StatementOnly {
        no_claim: id!(BridgeNoClaimIdV1, byte),
    }
}

fn proof_reference(byte: u8) -> BridgeProofStateV1 {
    BridgeProofStateV1::ProofArtifactReferenced {
        artifact: id!(BridgeProofArtifactIdV1, byte),
        verifier: id!(BridgeVerifierIdV1, 21),
        formal_system: id!(BridgeFormalSystemVersionIdV1, 23),
        no_claim: id!(BridgeNoClaimIdV1, 222),
    }
}

fn fixture() -> BridgeLatticeSpecV1 {
    let positive = convention(CountOrientationV1::Positive);
    let multiplicity = id!(BridgeMultiplicityRuleIdV1, 18);
    let maslov = MaslovIndexObjectV1 {
        path: id!(LagrangianPathIdV1, 1),
        reference: id!(LagrangianReferenceIdV1, 2),
        symplectic_form: id!(BridgeSymplecticFormIdV1, 3),
        crossing_forms: id!(CrossingFormFamilyIdV1, 4),
        convention: positive,
        multiplicity,
    };
    let krein = KreinSpectralFlowObjectV1 {
        path: id!(PontryaginPathIdV1, 5),
        krein_form: id!(BridgeKreinFormIdV1, 6),
        neutral_policy: NeutralKreinPolicyV1::Excluded {
            witness: id!(BridgeHypothesisWitnessIdV1, 7),
        },
        convention: positive,
        multiplicity,
    };
    let evans = EvansWindingObjectV1 {
        family: id!(AnalyticFredholmFamilyIdV1, 8),
        analytic_domain: id!(EvansAnalyticDomainIdV1, 9),
        contour: id!(EvansContourIdV1, 10),
        essential_spectrum_exclusion: id!(EssentialSpectrumExclusionIdV1, 11),
        normalization: id!(EvansNormalizationIdV1, 12),
        convention: positive,
        multiplicity,
    };
    let machine = Some(MachineInstabilityCountObjectV1 {
        model: id!(MachineInstabilityModelIdV1, 13),
        parameterization: id!(BridgeParameterizationIdV1, 14),
        convention: positive,
        multiplicity,
    });

    let finite = node(
        31,
        BridgeTheoremScopeV1::ClassicalFiniteHamiltonian { phase_dimension: 4 },
        BridgeConclusionV1::PairwiseEquality {
            left: BridgeCountKindV1::Maslov,
            right: BridgeCountKindV1::Krein,
            left_to_common: transform(41),
            right_to_common: transform(42),
        },
        proof_reference(151),
    );
    let periodic = node(
        32,
        BridgeTheoremScopeV1::PeriodicMonodromy {
            state_dimension: 4,
            monodromy_map: id!(BridgeCorrespondenceMapIdV1, 43),
        },
        BridgeConclusionV1::TripleEquality {
            maslov_to_common: transform(44),
            krein_to_common: transform(45),
            evans_to_common: transform(46),
        },
        statement_only(152),
    );
    let spatial = node(
        33,
        BridgeTheoremScopeV1::SpatialDynamicsEvans {
            stable_dimension: 2,
            unstable_dimension: 2,
            spatial_map: id!(BridgeCorrespondenceMapIdV1, 47),
        },
        BridgeConclusionV1::TripleEquality {
            maslov_to_common: transform(48),
            krein_to_common: transform(49),
            evans_to_common: transform(50),
        },
        proof_reference(153),
    );
    let maximal = node(
        34,
        BridgeTheoremScopeV1::MaximalMaslovKreinEvans {
            finite_to_periodic: id!(BridgeCorrespondenceMapIdV1, 51),
            periodic_to_evans: id!(BridgeCorrespondenceMapIdV1, 52),
        },
        BridgeConclusionV1::TripleEquality {
            maslov_to_common: transform(53),
            krein_to_common: transform(54),
            evans_to_common: transform(55),
        },
        statement_only(154),
    );
    let mut machine_node = node(
        35,
        BridgeTheoremScopeV1::MachineInstabilityCorollary {
            interpretation: id!(BridgeCorrespondenceMapIdV1, 56),
        },
        BridgeConclusionV1::SpectralToMachine {
            spectral: BridgeCountKindV1::Evans,
            spectral_to_machine: transform(57),
            interpretation: id!(BridgeCorrespondenceMapIdV1, 56),
        },
        statement_only(155),
    );
    machine_node.parameterization = id!(BridgeParameterizationIdV1, 14);

    let edge = |stronger: u8, weaker: u8, map: u8| BridgeImplicationV1 {
        stronger: id!(BridgeTheoremNodeIdV1, stronger),
        weaker: id!(BridgeTheoremNodeIdV1, weaker),
        projection: id!(BridgeCorrespondenceMapIdV1, map),
        state: BridgeImplicationStateV1::WitnessReferenced {
            witness: BridgeHypothesisWitnessIdV1::from_bytes([map.wrapping_add(90); 32]),
        },
    };
    BridgeLatticeSpecV1::new(
        id!(BridgeStatementVersionIdV1, 19),
        maslov,
        krein,
        evans,
        machine,
        vec![finite, periodic, spatial, maximal, machine_node],
        vec![
            edge(34, 32, 61),
            edge(34, 33, 62),
            edge(32, 31, 63),
            edge(33, 31, 64),
            edge(34, 35, 65),
        ],
        BridgeTrustedComputingBaseV1 {
            tcb: id!(BridgeTcbIdV1, 20),
            verifier: id!(BridgeVerifierIdV1, 21),
            policy: id!(BridgeProofPolicyIdV1, 22),
            formal_system: id!(BridgeFormalSystemVersionIdV1, 23),
            no_claim: id!(BridgeNoClaimIdV1, 24),
        },
        BridgeValidationBudgetV1 {
            max_nodes: 16,
            max_implications: 32,
            max_hypotheses_per_node: 32,
            max_falsifiers_per_node: 16,
        },
    )
}

fn rebuild(
    original: &BridgeLatticeSpecV1,
    schema_version: u32,
    nodes: Vec<BridgeTheoremNodeV1>,
    implications: Vec<BridgeImplicationV1>,
    machine: Option<MachineInstabilityCountObjectV1>,
) -> BridgeLatticeSpecV1 {
    BridgeLatticeSpecV1::with_schema_version(
        schema_version,
        original.statement_version(),
        original.maslov(),
        original.krein(),
        original.evans(),
        machine,
        nodes,
        implications,
        original.tcb(),
        original.budget(),
    )
}

#[test]
fn bridge_replay_is_permutation_invariant_and_unauthoritative() {
    with_cx(false, |cx| {
        let base = fixture();
        let validated = validate_bridge_lattice_v1(base.clone(), cx).expect("fixture is valid");
        assert_eq!(
            validated.scientific_authority(),
            BridgeScientificAuthorityV1::ScientificCorrectnessNotProven
        );
        assert!(validated.dispositions().iter().any(|record| {
            record.node == id!(BridgeTheoremNodeIdV1, 31)
                && record.disposition == BridgeNodeDispositionV1::ReferencedNotVerified
        }));
        assert!(validated.dispositions().iter().any(|record| {
            record.node == id!(BridgeTheoremNodeIdV1, 34)
                && record.disposition == BridgeNodeDispositionV1::ConjectureOnly
        }));

        let mut nodes = base.nodes().to_vec();
        nodes.reverse();
        for node in &mut nodes {
            node.hypotheses.reverse();
            node.falsifiers.reverse();
        }
        let mut implications = base.implications().to_vec();
        implications.reverse();
        let permuted = rebuild(
            &base,
            BRIDGE_LATTICE_SCHEMA_VERSION_V1,
            nodes,
            implications,
            base.machine(),
        );
        let replay = validate_bridge_lattice_v1(permuted, cx).expect("permutation is canonical");
        assert_eq!(validated.lattice_id(), replay.lattice_id());
        assert_eq!(
            validated.identity_receipt().canonical_preimage(),
            replay.identity_receipt().canonical_preimage()
        );
    });
}

#[test]
fn executable_orientation_and_endpoint_transforms_are_exact() {
    let source = SignedCountConventionV1 {
        orientation: CountOrientationV1::Positive,
        endpoints: EndpointConventionV1 {
            left: EndpointWeightV1::Excluded,
            right: EndpointWeightV1::Excluded,
        },
    };
    let target = SignedCountConventionV1 {
        orientation: CountOrientationV1::Negative,
        endpoints: EndpointConventionV1 {
            left: EndpointWeightV1::Full,
            right: EndpointWeightV1::Excluded,
        },
    };
    let derived = derive_convention_transform_v1(
        id!(BridgeCorrespondenceMapIdV1, 170),
        source,
        target,
        EndpointSignatureTraceV1 {
            left: None,
            right: Some(3),
        },
    )
    .expect("canonical-right signature resolves the reversed initial endpoint");
    assert_eq!(derived.sign(), -1);
    assert_eq!(derived.doubled_offset(), -6);
    assert_eq!(
        derived.apply(DoubledSignedCountV1(8)),
        Some(DoubledSignedCountV1(-14))
    );
    assert_eq!(
        derive_convention_transform_v1(
            id!(BridgeCorrespondenceMapIdV1, 170),
            source,
            target,
            EndpointSignatureTraceV1 {
                left: None,
                right: None,
            },
        ),
        Err(ConventionTransformErrorV1::RightEndpointUnresolved)
    );

    let asymmetric_source = SignedCountConventionV1 {
        orientation: CountOrientationV1::Positive,
        endpoints: EndpointConventionV1 {
            left: EndpointWeightV1::Full,
            right: EndpointWeightV1::Excluded,
        },
    };
    let asymmetric_target = SignedCountConventionV1 {
        orientation: CountOrientationV1::Negative,
        endpoints: asymmetric_source.endpoints,
    };
    let swapped = derive_convention_transform_v1(
        id!(BridgeCorrespondenceMapIdV1, 172),
        asymmetric_source,
        asymmetric_target,
        EndpointSignatureTraceV1 {
            left: Some(3),
            right: Some(5),
        },
    )
    .expect("both canonical endpoint signatures are resolved");
    assert_eq!(swapped.sign(), -1);
    assert_eq!(swapped.doubled_offset(), -4);
    let overflowing =
        SignedCountTransformV1::new(id!(BridgeCorrespondenceMapIdV1, 171), 1, i64::MAX)
            .expect("transform shape is valid");
    assert_eq!(overflowing.apply(DoubledSignedCountV1(1)), None);
}

#[test]
fn duplicate_node_refusals_have_canonical_issue_order() {
    with_cx(false, |cx| {
        let base = fixture();
        let mut malformed = base
            .nodes()
            .iter()
            .find(|node| node.id == id!(BridgeTheoremNodeIdV1, 31))
            .expect("finite fixture node exists")
            .clone();
        malformed.scope = BridgeTheoremScopeV1::ClassicalFiniteHamiltonian { phase_dimension: 3 };

        let mut first_nodes = base.nodes().to_vec();
        first_nodes.insert(0, malformed.clone());
        let first = rebuild(
            &base,
            BRIDGE_LATTICE_SCHEMA_VERSION_V1,
            first_nodes,
            base.implications().to_vec(),
            base.machine(),
        );
        let first = validate_bridge_lattice_v1(first, cx).expect_err("duplicate node refuses");

        let mut second_nodes = base.nodes().to_vec();
        second_nodes.push(malformed.clone());
        second_nodes.reverse();
        let second = rebuild(
            &base,
            BRIDGE_LATTICE_SCHEMA_VERSION_V1,
            second_nodes,
            base.implications().to_vec(),
            base.machine(),
        );
        let second = validate_bridge_lattice_v1(second, cx).expect_err("permutation refuses");

        assert_eq!(first.issues(), second.issues());
        assert!(
            first
                .issues()
                .contains(&BridgeValidationIssueV1::DuplicateNode {
                    node: id!(BridgeTheoremNodeIdV1, 31),
                })
        );
        assert!(
            first
                .issues()
                .contains(&BridgeValidationIssueV1::InvalidScopeDimension {
                    node: id!(BridgeTheoremNodeIdV1, 31),
                },)
        );

        let hypothesis = malformed.hypotheses[0];
        let mut shorter = malformed.clone();
        shorter.hypotheses = vec![hypothesis; 33];
        let mut longer = malformed;
        longer.hypotheses = vec![hypothesis; 34];
        let first = rebuild(
            &base,
            BRIDGE_LATTICE_SCHEMA_VERSION_V1,
            vec![longer.clone(), shorter.clone()],
            Vec::new(),
            base.machine(),
        );
        let second = rebuild(
            &base,
            BRIDGE_LATTICE_SCHEMA_VERSION_V1,
            vec![shorter, longer],
            Vec::new(),
            base.machine(),
        );
        let first = validate_bridge_lattice_v1(first, cx).expect_err("oversized node refuses");
        let second = validate_bridge_lattice_v1(second, cx).expect_err("permutation refuses");
        assert_eq!(first.issues(), second.issues());
        assert_eq!(
            first.issues(),
            &[BridgeValidationIssueV1::TooManyHypotheses {
                node: id!(BridgeTheoremNodeIdV1, 31),
                found: 33,
                limit: 32,
            }]
        );
    });
}

#[test]
fn periodic_maslov_scope_requires_even_state_dimension() {
    with_cx(false, |cx| {
        let base = fixture();
        let mut nodes = base.nodes().to_vec();
        let periodic = nodes
            .iter_mut()
            .find(|node| node.id == id!(BridgeTheoremNodeIdV1, 32))
            .expect("periodic fixture node exists");
        periodic.scope = match periodic.scope {
            BridgeTheoremScopeV1::PeriodicMonodromy { monodromy_map, .. } => {
                BridgeTheoremScopeV1::PeriodicMonodromy {
                    state_dimension: 3,
                    monodromy_map,
                }
            }
            _ => unreachable!("fixture node has periodic scope"),
        };
        let malformed = rebuild(
            &base,
            BRIDGE_LATTICE_SCHEMA_VERSION_V1,
            nodes,
            base.implications().to_vec(),
            base.machine(),
        );
        let report =
            validate_bridge_lattice_v1(malformed, cx).expect_err("odd symplectic scope refuses");
        assert!(
            report
                .issues()
                .contains(&BridgeValidationIssueV1::InvalidScopeDimension {
                    node: id!(BridgeTheoremNodeIdV1, 32),
                },)
        );

        let mut nodes = base.nodes().to_vec();
        let periodic = nodes
            .iter_mut()
            .find(|node| node.id == id!(BridgeTheoremNodeIdV1, 32))
            .expect("periodic fixture node exists");
        periodic.scope = match periodic.scope {
            BridgeTheoremScopeV1::PeriodicMonodromy { monodromy_map, .. } => {
                BridgeTheoremScopeV1::PeriodicMonodromy {
                    state_dimension: 3,
                    monodromy_map,
                }
            }
            _ => unreachable!("fixture node has periodic scope"),
        };
        periodic.conclusion = BridgeConclusionV1::PairwiseEquality {
            left: BridgeCountKindV1::Krein,
            right: BridgeCountKindV1::Evans,
            left_to_common: transform(173),
            right_to_common: transform(174),
        };
        let general_periodic = rebuild(
            &base,
            BRIDGE_LATTICE_SCHEMA_VERSION_V1,
            nodes,
            base.implications().to_vec(),
            base.machine(),
        );
        validate_bridge_lattice_v1(general_periodic, cx)
            .expect("odd general periodic scope without Maslov remains representable");
    });
}

#[test]
fn required_hypotheses_and_falsifiers_fail_closed() {
    with_cx(false, |cx| {
        let base = fixture();
        let mut nodes = base.nodes().to_vec();
        let spatial = nodes
            .iter_mut()
            .find(|node| node.id == id!(BridgeTheoremNodeIdV1, 33))
            .expect("spatial node exists");
        spatial
            .hypotheses
            .retain(|item| item.kind != BridgeHypothesisKindV1::EssentialSpectrumExcluded);
        spatial
            .falsifiers
            .retain(|item| item.kind != BridgeFalsifierKindV1::ContourContact);
        let malformed = rebuild(
            &base,
            BRIDGE_LATTICE_SCHEMA_VERSION_V1,
            nodes,
            base.implications().to_vec(),
            base.machine(),
        );
        let report = validate_bridge_lattice_v1(malformed, cx).expect_err("missing gates refuse");
        assert!(
            report
                .issues()
                .contains(&BridgeValidationIssueV1::MissingHypothesis {
                    node: id!(BridgeTheoremNodeIdV1, 33),
                    kind: BridgeHypothesisKindV1::EssentialSpectrumExcluded,
                })
        );
        assert!(
            report
                .issues()
                .contains(&BridgeValidationIssueV1::MissingFalsifier {
                    node: id!(BridgeTheoremNodeIdV1, 33),
                    kind: BridgeFalsifierKindV1::ContourContact,
                })
        );
    });
}

#[test]
fn neutral_tangential_and_essential_spectrum_falsifiers_are_mandatory() {
    with_cx(false, |cx| {
        for missing in [
            BridgeFalsifierKindV1::NeutralKreinSignature,
            BridgeFalsifierKindV1::TangentialCrossing,
            BridgeFalsifierKindV1::EssentialSpectrumContact,
        ] {
            let base = fixture();
            let mut nodes = base.nodes().to_vec();
            let maximal = nodes
                .iter_mut()
                .find(|node| node.id == id!(BridgeTheoremNodeIdV1, 34))
                .expect("maximal node exists");
            maximal.falsifiers.retain(|item| item.kind != missing);
            let malformed = rebuild(
                &base,
                BRIDGE_LATTICE_SCHEMA_VERSION_V1,
                nodes,
                base.implications().to_vec(),
                base.machine(),
            );
            let report =
                validate_bridge_lattice_v1(malformed, cx).expect_err("falsifier gap refuses");
            assert!(
                report
                    .issues()
                    .contains(&BridgeValidationIssueV1::MissingFalsifier {
                        node: id!(BridgeTheoremNodeIdV1, 34),
                        kind: missing,
                    })
            );
        }
    });
}

#[test]
fn unresolved_neutral_krein_type_cannot_masquerade_as_a_witness() {
    with_cx(false, |cx| {
        let base = fixture();
        let mut krein = base.krein();
        krein.neutral_policy = NeutralKreinPolicyV1::Unresolved {
            no_claim: id!(BridgeNoClaimIdV1, 175),
        };
        let overclaim = BridgeLatticeSpecV1::new(
            base.statement_version(),
            base.maslov(),
            krein,
            base.evans(),
            base.machine(),
            base.nodes().to_vec(),
            base.implications().to_vec(),
            base.tcb(),
            base.budget(),
        );
        let report =
            validate_bridge_lattice_v1(overclaim, cx).expect_err("neutral overclaim refuses");
        assert!(
            report
                .issues()
                .contains(&BridgeValidationIssueV1::NeutralPolicyOverclaim {
                    node: id!(BridgeTheoremNodeIdV1, 31),
                })
        );

        let mut nodes = base.nodes().to_vec();
        for node in &mut nodes {
            if node
                .conclusion
                .count_kinds()
                .contains(&BridgeCountKindV1::Krein)
            {
                let neutral = node
                    .hypotheses
                    .iter_mut()
                    .find(|item| item.kind == BridgeHypothesisKindV1::NeutralSignaturePolicyClosed)
                    .expect("fixture carries the neutral hypothesis");
                neutral.state = BridgeHypothesisStateV1::Unresolved {
                    no_claim: id!(BridgeNoClaimIdV1, 176),
                };
            }
        }
        let honest = BridgeLatticeSpecV1::new(
            base.statement_version(),
            base.maslov(),
            krein,
            base.evans(),
            base.machine(),
            nodes,
            base.implications().to_vec(),
            base.tcb(),
            base.budget(),
        );
        let honest = validate_bridge_lattice_v1(honest, cx)
            .expect("unresolved neutral type remains an honest conjecture");
        assert!(honest.dispositions().iter().any(|record| {
            record.node == id!(BridgeTheoremNodeIdV1, 31)
                && record.disposition == BridgeNodeDispositionV1::ConjectureOnly
        }));
    });
}

#[test]
fn machine_count_requires_an_explicit_distinct_object_and_map() {
    with_cx(false, |cx| {
        let base = fixture();
        let malformed = rebuild(
            &base,
            BRIDGE_LATTICE_SCHEMA_VERSION_V1,
            base.nodes().to_vec(),
            base.implications().to_vec(),
            None,
        );
        let report = validate_bridge_lattice_v1(malformed, cx).expect_err("machine object missing");
        assert!(
            report
                .issues()
                .contains(&BridgeValidationIssueV1::MachineObjectMissing {
                    node: id!(BridgeTheoremNodeIdV1, 35),
                })
        );

        let mut nodes = base.nodes().to_vec();
        let machine = nodes
            .iter_mut()
            .find(|node| node.id == id!(BridgeTheoremNodeIdV1, 35))
            .expect("machine node exists");
        machine.conclusion = BridgeConclusionV1::SpectralToMachine {
            spectral: BridgeCountKindV1::MachineInstability,
            spectral_to_machine: transform(57),
            interpretation: id!(BridgeCorrespondenceMapIdV1, 56),
        };
        let malformed = rebuild(
            &base,
            BRIDGE_LATTICE_SCHEMA_VERSION_V1,
            nodes,
            base.implications().to_vec(),
            base.machine(),
        );
        let report = validate_bridge_lattice_v1(malformed, cx).expect_err("self-label refuses");
        assert!(
            report
                .issues()
                .contains(&BridgeValidationIssueV1::InvalidMachineSource {
                    node: id!(BridgeTheoremNodeIdV1, 35),
                })
        );

        let mut nodes = base.nodes().to_vec();
        let finite = nodes
            .iter_mut()
            .find(|node| node.id == id!(BridgeTheoremNodeIdV1, 31))
            .expect("finite node exists");
        finite.conclusion = BridgeConclusionV1::PairwiseEquality {
            left: BridgeCountKindV1::Maslov,
            right: BridgeCountKindV1::MachineInstability,
            left_to_common: transform(58),
            right_to_common: transform(59),
        };
        let bypass = rebuild(
            &base,
            BRIDGE_LATTICE_SCHEMA_VERSION_V1,
            nodes,
            base.implications().to_vec(),
            base.machine(),
        );
        let report =
            validate_bridge_lattice_v1(bypass, cx).expect_err("generic equality cannot bypass map");
        assert!(report.issues().contains(
            &BridgeValidationIssueV1::MachineCountNeedsExplicitRelation {
                node: id!(BridgeTheoremNodeIdV1, 31),
            }
        ));
    });
}

#[test]
fn implication_cycles_and_missing_weaker_scopes_refuse() {
    with_cx(false, |cx| {
        let base = fixture();
        let mut implications = base.implications().to_vec();
        implications.push(BridgeImplicationV1 {
            stronger: id!(BridgeTheoremNodeIdV1, 31),
            weaker: id!(BridgeTheoremNodeIdV1, 34),
            projection: id!(BridgeCorrespondenceMapIdV1, 180),
            state: BridgeImplicationStateV1::Unresolved {
                no_claim: id!(BridgeNoClaimIdV1, 181),
            },
        });
        let cyclic = rebuild(
            &base,
            BRIDGE_LATTICE_SCHEMA_VERSION_V1,
            base.nodes().to_vec(),
            implications,
            base.machine(),
        );
        let report = validate_bridge_lattice_v1(cyclic, cx).expect_err("cycle refuses");
        assert!(
            report
                .issues()
                .contains(&BridgeValidationIssueV1::ImplicationCycle)
        );

        let nodes: Vec<_> = base
            .nodes()
            .iter()
            .filter(|node| node.id != id!(BridgeTheoremNodeIdV1, 33))
            .cloned()
            .collect();
        let implications: Vec<_> = base
            .implications()
            .iter()
            .filter(|edge| {
                edge.stronger != id!(BridgeTheoremNodeIdV1, 33)
                    && edge.weaker != id!(BridgeTheoremNodeIdV1, 33)
            })
            .copied()
            .collect();
        let incomplete = rebuild(
            &base,
            BRIDGE_LATTICE_SCHEMA_VERSION_V1,
            nodes,
            implications,
            base.machine(),
        );
        let report = validate_bridge_lattice_v1(incomplete, cx).expect_err("scope gap refuses");
        assert!(
            report
                .issues()
                .contains(&BridgeValidationIssueV1::MissingScope {
                    scope: "spatial-dynamics-Evans",
                })
        );
    });
}

#[test]
fn implication_scope_order_and_projection_branches_fail_closed() {
    with_cx(false, |cx| {
        let base = fixture();
        let mut reversed = base.implications().to_vec();
        let finite_projection = reversed
            .iter_mut()
            .find(|edge| {
                edge.stronger == id!(BridgeTheoremNodeIdV1, 32)
                    && edge.weaker == id!(BridgeTheoremNodeIdV1, 31)
            })
            .expect("periodic-to-finite fixture edge exists");
        core::mem::swap(
            &mut finite_projection.stronger,
            &mut finite_projection.weaker,
        );
        let reversed = rebuild(
            &base,
            BRIDGE_LATTICE_SCHEMA_VERSION_V1,
            base.nodes().to_vec(),
            reversed,
            base.machine(),
        );
        let report = validate_bridge_lattice_v1(reversed, cx)
            .expect_err("lower-to-higher extension edge refuses");
        assert!(report.issues().contains(
            &BridgeValidationIssueV1::ImplicationScopeOrderMismatch {
                stronger: id!(BridgeTheoremNodeIdV1, 31),
                weaker: id!(BridgeTheoremNodeIdV1, 32),
            },
        ));

        let machine_only: Vec<_> = base
            .implications()
            .iter()
            .filter(|edge| edge.weaker == id!(BridgeTheoremNodeIdV1, 35))
            .copied()
            .collect();
        let disconnected = rebuild(
            &base,
            BRIDGE_LATTICE_SCHEMA_VERSION_V1,
            base.nodes().to_vec(),
            machine_only,
            base.machine(),
        );
        let report = validate_bridge_lattice_v1(disconnected, cx)
            .expect_err("machine edge cannot replace spectral projection branches");
        for scope in ["periodic-monodromy", "spatial-dynamics-Evans"] {
            assert!(report.issues().contains(
                &BridgeValidationIssueV1::MaximalProjectionCoverageMissing {
                    node: id!(BridgeTheoremNodeIdV1, 34),
                    scope,
                },
            ));
        }

        let mut broken_periodic_descent = base.implications().to_vec();
        let periodic_to_finite = broken_periodic_descent
            .iter_mut()
            .find(|edge| {
                edge.stronger == id!(BridgeTheoremNodeIdV1, 32)
                    && edge.weaker == id!(BridgeTheoremNodeIdV1, 31)
            })
            .expect("periodic-to-finite fixture edge exists");
        periodic_to_finite.state = BridgeImplicationStateV1::Refuted {
            counterexample: id!(BridgeCounterexampleIdV1, 221),
        };
        let broken_periodic = rebuild(
            &base,
            BRIDGE_LATTICE_SCHEMA_VERSION_V1,
            base.nodes().to_vec(),
            broken_periodic_descent,
            base.machine(),
        );
        let report = validate_bridge_lattice_v1(broken_periodic, cx)
            .expect_err("a visible branch must still reach the classical scope");
        assert!(report.issues().contains(
            &BridgeValidationIssueV1::MaximalProjectionCoverageMissing {
                node: id!(BridgeTheoremNodeIdV1, 34),
                scope: "periodic-monodromy",
            },
        ));
        assert!(!report.issues().contains(
            &BridgeValidationIssueV1::MaximalProjectionCoverageMissing {
                node: id!(BridgeTheoremNodeIdV1, 34),
                scope: "spatial-dynamics-Evans",
            },
        ));

        let mut serial_branches = base.implications().to_vec();
        let maximal_to_spatial = serial_branches
            .iter_mut()
            .find(|edge| {
                edge.stronger == id!(BridgeTheoremNodeIdV1, 34)
                    && edge.weaker == id!(BridgeTheoremNodeIdV1, 33)
            })
            .expect("maximal-to-spatial fixture edge exists");
        maximal_to_spatial.stronger = id!(BridgeTheoremNodeIdV1, 32);
        let serial = rebuild(
            &base,
            BRIDGE_LATTICE_SCHEMA_VERSION_V1,
            base.nodes().to_vec(),
            serial_branches,
            base.machine(),
        );
        let report = validate_bridge_lattice_v1(serial, cx)
            .expect_err("a serial periodic-to-spatial path is not two direct branches");
        assert!(report.issues().contains(
            &BridgeValidationIssueV1::MaximalProjectionCoverageMissing {
                node: id!(BridgeTheoremNodeIdV1, 34),
                scope: "spatial-dynamics-Evans",
            },
        ));
        assert!(!report.issues().contains(
            &BridgeValidationIssueV1::MaximalProjectionCoverageMissing {
                node: id!(BridgeTheoremNodeIdV1, 34),
                scope: "periodic-monodromy",
            },
        ));

        let mut refuted_nodes = base.nodes().to_vec();
        let spatial = refuted_nodes
            .iter_mut()
            .find(|node| node.id == id!(BridgeTheoremNodeIdV1, 33))
            .expect("spatial fixture node exists");
        spatial.proof = BridgeProofStateV1::Refuted {
            counterexample: id!(BridgeCounterexampleIdV1, 219),
        };
        let refuted_node = rebuild(
            &base,
            BRIDGE_LATTICE_SCHEMA_VERSION_V1,
            refuted_nodes,
            base.implications().to_vec(),
            base.machine(),
        );
        let report = validate_bridge_lattice_v1(refuted_node, cx)
            .expect_err("a refuted theorem node provides no positive branch coverage");
        assert!(report.issues().contains(
            &BridgeValidationIssueV1::MaximalProjectionCoverageMissing {
                node: id!(BridgeTheoremNodeIdV1, 34),
                scope: "spatial-dynamics-Evans",
            },
        ));

        let mut refuted_branches = base.implications().to_vec();
        for edge in &mut refuted_branches {
            if edge.stronger == id!(BridgeTheoremNodeIdV1, 34)
                && matches!(
                    edge.weaker,
                    node if node == id!(BridgeTheoremNodeIdV1, 32)
                        || node == id!(BridgeTheoremNodeIdV1, 33)
                )
            {
                edge.state = BridgeImplicationStateV1::Refuted {
                    counterexample: id!(BridgeCounterexampleIdV1, 220),
                };
            }
        }
        let refuted = rebuild(
            &base,
            BRIDGE_LATTICE_SCHEMA_VERSION_V1,
            base.nodes().to_vec(),
            refuted_branches,
            base.machine(),
        );
        let report = validate_bridge_lattice_v1(refuted, cx)
            .expect_err("refuted branch edges provide no positive coverage");
        for scope in ["periodic-monodromy", "spatial-dynamics-Evans"] {
            assert!(report.issues().contains(
                &BridgeValidationIssueV1::MaximalProjectionCoverageMissing {
                    node: id!(BridgeTheoremNodeIdV1, 34),
                    scope,
                },
            ));
        }
    });
}

#[test]
fn refuted_hypothesis_stays_visible_as_refuted_not_truth() {
    with_cx(false, |cx| {
        let base = fixture();
        let mut nodes = base.nodes().to_vec();
        let spatial = nodes
            .iter_mut()
            .find(|node| node.id == id!(BridgeTheoremNodeIdV1, 33))
            .expect("spatial node exists");
        let mut replacement = spatial.clone();
        replacement.id = id!(BridgeTheoremNodeIdV1, 36);
        let hypothesis = spatial
            .hypotheses
            .iter_mut()
            .find(|item| item.kind == BridgeHypothesisKindV1::SpatialDichotomy)
            .expect("spatial dichotomy exists");
        hypothesis.state = BridgeHypothesisStateV1::Refuted {
            counterexample: id!(BridgeCounterexampleIdV1, 190),
        };
        nodes.push(replacement);
        let mut implications = base.implications().to_vec();
        let mut maximal_to_replacement = *implications
            .iter()
            .find(|edge| {
                edge.stronger == id!(BridgeTheoremNodeIdV1, 34)
                    && edge.weaker == id!(BridgeTheoremNodeIdV1, 33)
            })
            .expect("maximal-to-spatial fixture edge exists");
        maximal_to_replacement.weaker = id!(BridgeTheoremNodeIdV1, 36);
        maximal_to_replacement.projection = id!(BridgeCorrespondenceMapIdV1, 191);
        let mut replacement_to_finite = *implications
            .iter()
            .find(|edge| {
                edge.stronger == id!(BridgeTheoremNodeIdV1, 33)
                    && edge.weaker == id!(BridgeTheoremNodeIdV1, 31)
            })
            .expect("spatial-to-finite fixture edge exists");
        replacement_to_finite.stronger = id!(BridgeTheoremNodeIdV1, 36);
        replacement_to_finite.projection = id!(BridgeCorrespondenceMapIdV1, 192);
        implications.extend([maximal_to_replacement, replacement_to_finite]);
        let refuted = rebuild(
            &base,
            BRIDGE_LATTICE_SCHEMA_VERSION_V1,
            nodes,
            implications,
            base.machine(),
        );
        let validated = validate_bridge_lattice_v1(refuted, cx).expect("refutation is valid data");
        assert!(validated.dispositions().iter().any(|record| {
            record.node == id!(BridgeTheoremNodeIdV1, 33)
                && record.disposition == BridgeNodeDispositionV1::Refuted
        }));
        assert_eq!(
            validated.scientific_authority(),
            BridgeScientificAuthorityV1::ScientificCorrectnessNotProven
        );
    });
}

#[test]
fn proof_references_must_match_the_declared_tcb() {
    with_cx(false, |cx| {
        let base = fixture();
        let mut nodes = base.nodes().to_vec();
        let finite = nodes
            .iter_mut()
            .find(|node| node.id == id!(BridgeTheoremNodeIdV1, 31))
            .expect("finite node exists");
        finite.proof = BridgeProofStateV1::ProofArtifactReferenced {
            artifact: id!(BridgeProofArtifactIdV1, 195),
            verifier: id!(BridgeVerifierIdV1, 196),
            formal_system: base.tcb().formal_system,
            no_claim: id!(BridgeNoClaimIdV1, 197),
        };
        let malformed = rebuild(
            &base,
            BRIDGE_LATTICE_SCHEMA_VERSION_V1,
            nodes,
            base.implications().to_vec(),
            base.machine(),
        );
        let report =
            validate_bridge_lattice_v1(malformed, cx).expect_err("foreign verifier refuses");
        assert!(
            report
                .issues()
                .contains(&BridgeValidationIssueV1::ProofTcbMismatch {
                    node: id!(BridgeTheoremNodeIdV1, 31),
                })
        );
    });
}

#[test]
fn schema_statement_and_convention_mutations_are_byte_visible() {
    with_cx(false, |cx| {
        let base = fixture();
        let valid = validate_bridge_lattice_v1(base.clone(), cx).expect("fixture is valid");

        let wrong_schema = rebuild(
            &base,
            99,
            base.nodes().to_vec(),
            base.implications().to_vec(),
            base.machine(),
        );
        let report =
            validate_bridge_lattice_v1(wrong_schema, cx).expect_err("schema mutation refuses");
        assert!(
            report
                .issues()
                .contains(&BridgeValidationIssueV1::UnsupportedSchemaVersion {
                    found: 99,
                    supported: BRIDGE_LATTICE_SCHEMA_VERSION_V1,
                })
        );

        let mut maslov = base.maslov();
        maslov.convention.orientation = CountOrientationV1::Negative;
        let mutated = BridgeLatticeSpecV1::new(
            base.statement_version(),
            maslov,
            base.krein(),
            base.evans(),
            base.machine(),
            base.nodes().to_vec(),
            base.implications().to_vec(),
            base.tcb(),
            base.budget(),
        );
        let mutated = validate_bridge_lattice_v1(mutated, cx).expect("mutation remains explicit");
        assert_ne!(valid.lattice_id(), mutated.lattice_id());

        let mut evans = base.evans();
        evans.contour = id!(EvansContourIdV1, 251);
        let contour_mutation = BridgeLatticeSpecV1::new(
            base.statement_version(),
            base.maslov(),
            base.krein(),
            evans,
            base.machine(),
            base.nodes().to_vec(),
            base.implications().to_vec(),
            base.tcb(),
            base.budget(),
        );
        let contour_mutation = validate_bridge_lattice_v1(contour_mutation, cx)
            .expect("contour deformation remains explicit statement data");
        assert_ne!(valid.lattice_id(), contour_mutation.lattice_id());

        let mut nodes = base.nodes().to_vec();
        nodes[0].proof = statement_only(250);
        let reviewer_mutation = rebuild(
            &base,
            BRIDGE_LATTICE_SCHEMA_VERSION_V1,
            nodes,
            base.implications().to_vec(),
            base.machine(),
        );
        let reviewer_mutation =
            validate_bridge_lattice_v1(reviewer_mutation, cx).expect("review mutation is explicit");
        assert_ne!(valid.lattice_id(), reviewer_mutation.lattice_id());
    });
}

#[test]
fn budgets_and_cancellation_refuse_before_identity_work() {
    with_cx(false, |cx| {
        let base = fixture();
        let invalid_budget = BridgeLatticeSpecV1::new(
            base.statement_version(),
            base.maslov(),
            base.krein(),
            base.evans(),
            base.machine(),
            base.nodes().to_vec(),
            base.implications().to_vec(),
            base.tcb(),
            BridgeValidationBudgetV1 {
                max_nodes: 0,
                ..base.budget()
            },
        );
        let report =
            validate_bridge_lattice_v1(invalid_budget, cx).expect_err("zero budget refuses");
        assert_eq!(report.issues(), &[BridgeValidationIssueV1::InvalidBudget]);

        let mut nodes = base.nodes().to_vec();
        let mut extra = nodes[0].clone();
        extra.id = id!(BridgeTheoremNodeIdV1, 200);
        nodes.push(extra);
        let constrained = BridgeLatticeSpecV1::new(
            base.statement_version(),
            base.maslov(),
            base.krein(),
            base.evans(),
            base.machine(),
            nodes,
            base.implications().to_vec(),
            base.tcb(),
            BridgeValidationBudgetV1 {
                max_nodes: 5,
                ..base.budget()
            },
        );
        let report =
            validate_bridge_lattice_v1(constrained, cx).expect_err("caller cap refuses early");
        assert_eq!(
            report.issues(),
            &[BridgeValidationIssueV1::TooManyNodes { found: 6, limit: 5 }]
        );
    });

    with_cx(true, |cx| {
        let report =
            validate_bridge_lattice_v1(fixture(), cx).expect_err("pre-cancelled context refuses");
        assert_eq!(report.issues(), &[BridgeValidationIssueV1::Cancelled]);
    });
}
