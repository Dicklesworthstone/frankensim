//! G0/G3 examples and bounded-cancellation regression for RD.X1 statements.

#![allow(clippy::wildcard_imports)]
#![allow(
    clippy::too_many_lines,
    reason = "each RD.X1 fixture keeps one topology or hypothesis-deletion narrative explicit"
)]

use asupersync::types::Budget;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::derived::*;
use fs_geom::exit_path::*;
use std::sync::OnceLock;

fn bytes(seed: u8) -> [u8; 32] {
    [seed; 32]
}

fn derived_units(system: DerivedUnitSystemIdV1, seed: u8) -> UnitBindingV1 {
    UnitBindingV1 {
        system,
        quantity: DerivedQuantityKindIdV1::from_bytes(bytes(seed)),
        scale_to_canonical: 1.0,
    }
}

fn derived_complex(
    seed: u8,
    chart: ConfigurationChartIdV1,
    role: DerivedComplexRoleV1,
) -> FiniteDerivedComplexV1 {
    FiniteDerivedComplexV1 {
        id: DerivedComplexIdV1::from_bytes(bytes(seed)),
        chart,
        role,
        spaces: vec![
            GradedSpaceV1 {
                degree: 0,
                dimension: 2,
                quantity: DerivedQuantityKindIdV1::from_bytes(bytes(seed.wrapping_add(1))),
            },
            GradedSpaceV1 {
                degree: 1,
                dimension: 1,
                quantity: DerivedQuantityKindIdV1::from_bytes(bytes(seed.wrapping_add(2))),
            },
        ],
        differentials: vec![ComplexDifferentialV1 {
            from_degree: 0,
            to_degree: 1,
            map: DerivedLinearMapIdV1::from_bytes(bytes(seed.wrapping_add(3))),
            square_zero_witness: DerivedWitnessIdV1::from_bytes(bytes(seed.wrapping_add(4))),
        }],
        resolution: FiniteResolutionV1 {
            id: DerivedResolutionIdV1::from_bytes(bytes(seed.wrapping_add(5))),
            min_degree: 0,
            max_degree: 1,
            max_basis_dimension: 2,
            truncation_order: 0,
            remainder: None,
        },
        computability: FiniteComputabilityV1::ExactFinite {
            kernel: DerivedWitnessIdV1::from_bytes(bytes(seed.wrapping_add(6))),
        },
    }
}

fn admitted_subject_ir() -> DerivedGeometryIrV1 {
    let chart = ConfigurationChartIdV1::from_bytes(bytes(6));
    let frame = DerivedFrameIdV1::from_bytes(bytes(4));
    let units = DerivedUnitSystemIdV1::from_bytes(bytes(5));
    let equality = EqualityConstraintIdV1::from_bytes(bytes(20));
    let local_model = DerivedLocalModelIdV1::from_bytes(bytes(40));
    let stratum = StratumIdV1::from_bytes(bytes(50));
    DerivedGeometryIrV1 {
        schema_version: DERIVED_GEOMETRY_SCHEMA_VERSION_V1,
        subject: DerivedSubjectIdV1::from_bytes(bytes(1)),
        model_version: DerivedModelVersionIdV1::from_bytes(bytes(2)),
        category: GeometricCategoryV1::Semialgebraic,
        coefficients: CoefficientSystemV1::RationalReal,
        frame,
        unit_system: units,
        locality: LocalityScopeV1::GermAt {
            chart,
            point: DerivedWitnessIdV1::from_bytes(bytes(7)),
        },
        compactness: CompactnessV1::RelativelyCompact {
            witness: DerivedWitnessIdV1::from_bytes(bytes(8)),
        },
        charts: vec![ConfigurationChartV1 {
            id: chart,
            class: ConfigurationChartClassV1::Semialgebraic,
            coordinate_dimension: 2,
            ambient_dimension: 2,
            frame,
            coordinates: derived_units(units, 10),
            locality: LocalityScopeV1::GermAt {
                chart,
                point: DerivedWitnessIdV1::from_bytes(bytes(7)),
            },
            compactness: CompactnessV1::RelativelyCompact {
                witness: DerivedWitnessIdV1::from_bytes(bytes(8)),
            },
            regularity: RegularityClassV1::Polynomial,
            computability: FiniteComputabilityV1::ExactFinite {
                kernel: DerivedWitnessIdV1::from_bytes(bytes(9)),
            },
        }],
        equalities: vec![EqualityConstraintGermV1 {
            id: equality,
            chart,
            codomain_dimension: 1,
            equation: LocalFunctionEncodingV1::Polynomial {
                polynomial: PolynomialIdV1::from_bytes(bytes(21)),
                variables: 2,
                degree: 2,
            },
            regularity: RegularityClassV1::Polynomial,
            units: derived_units(units, 11),
            computability: FiniteComputabilityV1::ExactFinite {
                kernel: DerivedWitnessIdV1::from_bytes(bytes(22)),
            },
        }],
        inequalities: Vec::new(),
        boundaries: Vec::new(),
        contacts: Vec::new(),
        constitutive_data: Vec::new(),
        complexes: vec![
            derived_complex(30, chart, DerivedComplexRoleV1::Tangent),
            derived_complex(31, chart, DerivedComplexRoleV1::Cotangent),
            derived_complex(32, chart, DerivedComplexRoleV1::DeformationObstruction),
        ],
        local_models: vec![DerivedLocalModelV1 {
            id: local_model,
            chart,
            class: DerivedLocalModelClassV1::RegularCompleteIntersection,
            equalities: vec![equality],
            active_inequalities: Vec::new(),
            active_contacts: Vec::new(),
            constitutive_data: Vec::new(),
            tangent_complex: DerivedComplexIdV1::from_bytes(bytes(30)),
            cotangent_complex: DerivedComplexIdV1::from_bytes(bytes(31)),
            deformation_complex: DerivedComplexIdV1::from_bytes(bytes(32)),
            virtual_dimension: 1,
            locality: LocalityScopeV1::GermAt {
                chart,
                point: DerivedWitnessIdV1::from_bytes(bytes(7)),
            },
            presentation: PresentationScopeV1::Literal {
                no_claim: DerivedNoClaimIdV1::from_bytes(bytes(41)),
            },
        }],
        stratification: StratificationV1 {
            id: StratificationIdV1::from_bytes(bytes(3)),
            class: StratificationClassV1::FiniteIncidence,
            strata: vec![StratumSpecV1 {
                id: stratum,
                chart,
                local_model,
                dimension: 1,
                active_inequalities: Vec::new(),
                active_contacts: Vec::new(),
                relative_boundary: None,
                regularity: RegularityClassV1::Polynomial,
                compactness: CompactnessV1::RelativelyCompact {
                    witness: DerivedWitnessIdV1::from_bytes(bytes(51)),
                },
            }],
            incidences: Vec::new(),
            local_links: Vec::new(),
        },
        proof_state: DerivedProofStateV1::StructuralNoClaim {
            no_claim: DerivedNoClaimIdV1::from_bytes(bytes(60)),
        },
    }
}

static ADMITTED_SUBJECT: OnceLock<AdmittedDerivedGeometryV1> = OnceLock::new();

fn admitted_subject() -> &'static AdmittedDerivedGeometryV1 {
    ADMITTED_SUBJECT.get_or_init(|| {
        with_cx(false, |cx| {
            admit_derived_geometry_v1(
                admitted_subject_ir(),
                DerivedAdmissionBudgetV1::STANDARD,
                cx,
            )
            .expect("RD.X1 subject fixture admits through RD.1a")
        })
    })
}

fn with_cx<R>(cancelled: bool, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    if cancelled {
        gate.request();
    }
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x4558_4954,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn witness(seed: u8) -> ExitPathWitnessIdV1 {
    ExitPathWitnessIdV1::from_bytes(bytes(seed))
}

fn falsifier(seed: u8, kind: ExitPathFalsifierKindV1) -> ExitPathFalsifierV1 {
    ExitPathFalsifierV1 {
        id: ExitPathFalsifierIdV1::from_bytes(bytes(seed)),
        kind,
        left: ExitPathCountermodelIdV1::from_bytes(bytes(seed.wrapping_add(1))),
        right: ExitPathCountermodelIdV1::from_bytes(bytes(seed.wrapping_add(2))),
        witness: witness(seed.wrapping_add(3)),
    }
}

fn required_falsifiers() -> Vec<ExitPathFalsifierV1> {
    vec![
        falsifier(40, ExitPathFalsifierKindV1::DirectionReversal),
        falsifier(20, ExitPathFalsifierKindV1::SameIncidenceDifferentLink),
        falsifier(50, ExitPathFalsifierKindV1::HypothesisDeletion),
        falsifier(30, ExitPathFalsifierKindV1::SameIncidenceDifferentMonodromy),
    ]
}

fn regular_cell_ir() -> ExitPathFamilyIrV1 {
    let subject = admitted_subject();
    let subject_ir = subject.ir();
    let card = ExitPathTheoremCardIdV1::from_bytes(bytes(60));
    ExitPathFamilyIrV1::new(
        subject.id(),
        subject_ir.model_version,
        subject_ir.stratification.id,
        subject_ir.frame,
        subject_ir.unit_system,
        ExitStratifiedSpaceClassV1::FiniteRegularCell,
        ExitPathDirectionV1::Exit,
        ConstructibleVarianceV1::SheafContravariant,
        ExitPathConventionIdV1::from_bytes(bytes(6)),
        StratifiedPathEquivalenceV1::EndpointFixed {
            relation: StratifiedPathEquivalenceIdV1::from_bytes(bytes(7)),
        },
        CoefficientSystemV1::RationalReal,
        ExitPathHypothesesV1 {
            links: ConicalLinkHypothesisV1::Contractible {
                links: ExitLinkCatalogIdV1::from_bytes(bytes(8)),
                witness: witness(9),
            },
            monodromy: MonodromyHypothesisV1::Trivial {
                witness: witness(10),
            },
            constructibility: ConstructibilityHypothesisV1::LocallyConstantOnStrata {
                witness: witness(11),
            },
            properness: ExitPropernessHypothesisV1::Compact {
                witness: witness(12),
            },
            refinement: RefinementHypothesisV1::Identity {
                refinement: StratifiedRefinementIdV1::from_bytes(bytes(13)),
            },
            homotopy: HomotopyFidelityV1::IncidenceOnly,
        },
        required_falsifiers(),
        ExitPathTheoremStateV1::Preregistered { card },
        ExitPathTcbV1 {
            tcb: ExitPathTcbIdV1::from_bytes(bytes(61)),
            checker: ExitPathCheckerIdV1::from_bytes(bytes(62)),
            theorem_card: card,
        },
        ExitPathBudgetV1 {
            max_truncation: 3,
            max_referenced_artifact_slots: 10_000,
            max_implication_checks: 6,
            declared_wall_seconds: 5.0,
        },
    )
}

fn validate(ir: ExitPathFamilyIrV1) -> ValidatedExitPathFamilyV1 {
    with_cx(false, |cx| {
        validate_exit_path_family_v1(ir, admitted_subject(), cx).expect("fixture admits")
    })
}

fn assert_issue(ir: ExitPathFamilyIrV1, expected: ExitPathSemanticIssueV1) {
    let report = with_cx(false, |cx| {
        validate_exit_path_family_v1(ir, admitted_subject(), cx).expect_err("fixture must refuse")
    });
    assert!(
        report.issues().contains(&expected),
        "missing {expected:?} in {:?}",
        report.issues()
    );
}

fn state(
    admitted: &ValidatedExitPathFamilyV1,
    approximation: ExitPathApproximationV1,
) -> ExitPathNodeStateV1 {
    admitted
        .theorem_lattice()
        .iter()
        .find(|node| node.approximation == approximation)
        .expect("lattice node")
        .state
}

#[test]
fn regular_cells_admit_poset_but_do_not_erase_higher_unknowns() {
    let admitted = validate(regular_cell_ir());
    assert_eq!(
        state(&admitted, ExitPathApproximationV1::IncidencePoset),
        ExitPathNodeStateV1::SufficientStatement
    );
    assert_eq!(
        state(
            &admitted,
            ExitPathApproximationV1::StratumGroupoidEnrichedExitCategory,
        ),
        ExitPathNodeStateV1::Unknown {
            reason: ExitPathUnknownReasonV1::HomotopyDataInsufficient,
        }
    );
    assert_eq!(
        state(&admitted, ExitPathApproximationV1::FullHigherCategory),
        ExitPathNodeStateV1::Unknown {
            reason: ExitPathUnknownReasonV1::FullHigherDataInsufficient,
        }
    );
    assert_eq!(
        admitted.scientific_authority(),
        ExitPathScientificAuthorityV1::ScientificCorrectnessNotProven
    );
}

#[test]
fn cone_or_cusp_link_data_selects_groupoid_enriched_category_not_poset() {
    let mut cone = regular_cell_ir();
    cone.space_class = ExitStratifiedSpaceClassV1::ConicalSemialgebraic;
    cone.path_equivalence = StratifiedPathEquivalenceV1::HigherThrough {
        degree: 1,
        relation: StratifiedPathEquivalenceIdV1::from_bytes(bytes(70)),
    };
    cone.hypotheses.links = ConicalLinkHypothesisV1::RetainedThrough {
        links: ExitLinkCatalogIdV1::from_bytes(bytes(71)),
        degree: 1,
        witness: witness(72),
    };
    cone.hypotheses.monodromy = MonodromyHypothesisV1::Groupoids {
        groupoids: StratumGroupoidCatalogIdV1::from_bytes(bytes(73)),
        witness: witness(74),
    };
    cone.hypotheses.homotopy = HomotopyFidelityV1::RetainedThrough {
        degree: 1,
        coherence: HigherCoherenceArtifactIdV1::from_bytes(bytes(75)),
        witness: witness(76),
    };
    let admitted = validate(cone);
    assert_eq!(
        state(&admitted, ExitPathApproximationV1::IncidencePoset),
        ExitPathNodeStateV1::Unknown {
            reason: ExitPathUnknownReasonV1::LinkDataInsufficient,
        }
    );
    assert_eq!(
        state(
            &admitted,
            ExitPathApproximationV1::StratumGroupoidEnrichedExitCategory,
        ),
        ExitPathNodeStateV1::SufficientStatement
    );
}

#[test]
fn circular_strata_and_local_systems_require_richer_fidelity() {
    let mut circle = regular_cell_ir();
    circle.path_equivalence = StratifiedPathEquivalenceV1::HigherThrough {
        degree: 2,
        relation: StratifiedPathEquivalenceIdV1::from_bytes(bytes(80)),
    };
    circle.hypotheses.monodromy = MonodromyHypothesisV1::LocalSystemsThrough {
        groupoids: StratumGroupoidCatalogIdV1::from_bytes(bytes(81)),
        local_systems: LocalSystemCatalogIdV1::from_bytes(bytes(82)),
        degree: 2,
        witness: witness(83),
    };
    circle.hypotheses.homotopy = HomotopyFidelityV1::RetainedThrough {
        degree: 2,
        coherence: HigherCoherenceArtifactIdV1::from_bytes(bytes(84)),
        witness: witness(85),
    };
    let admitted = validate(circle);
    assert_eq!(
        state(&admitted, ExitPathApproximationV1::IncidencePoset),
        ExitPathNodeStateV1::Unknown {
            reason: ExitPathUnknownReasonV1::MonodromyDataInsufficient,
        }
    );
    assert_eq!(
        state(
            &admitted,
            ExitPathApproximationV1::SimplicialCategory {
                max_simplex_dimension: 2,
            }
        ),
        ExitPathNodeStateV1::SufficientStatement
    );
}

#[test]
fn finite_higher_and_full_higher_nodes_are_separate() {
    let mut finite = regular_cell_ir();
    finite.path_equivalence = StratifiedPathEquivalenceV1::HigherThrough {
        degree: 3,
        relation: StratifiedPathEquivalenceIdV1::from_bytes(bytes(90)),
    };
    finite.hypotheses.links = ConicalLinkHypothesisV1::RetainedThrough {
        links: ExitLinkCatalogIdV1::from_bytes(bytes(91)),
        degree: 3,
        witness: witness(92),
    };
    finite.hypotheses.monodromy = MonodromyHypothesisV1::LocalSystemsThrough {
        groupoids: StratumGroupoidCatalogIdV1::from_bytes(bytes(93)),
        local_systems: LocalSystemCatalogIdV1::from_bytes(bytes(94)),
        degree: 3,
        witness: witness(95),
    };
    finite.hypotheses.homotopy = HomotopyFidelityV1::RetainedThrough {
        degree: 3,
        coherence: HigherCoherenceArtifactIdV1::from_bytes(bytes(96)),
        witness: witness(97),
    };
    let finite_admitted = validate(finite.clone());
    assert_eq!(
        state(
            &finite_admitted,
            ExitPathApproximationV1::HigherTruncation { degree: 3 }
        ),
        ExitPathNodeStateV1::SufficientStatement
    );
    assert!(matches!(
        state(
            &finite_admitted,
            ExitPathApproximationV1::FullHigherCategory
        ),
        ExitPathNodeStateV1::Unknown { .. }
    ));

    let mut weak_equivalence = finite.clone();
    weak_equivalence.path_equivalence = StratifiedPathEquivalenceV1::EndpointFixed {
        relation: StratifiedPathEquivalenceIdV1::from_bytes(bytes(109)),
    };
    assert_eq!(
        state(
            &validate(weak_equivalence),
            ExitPathApproximationV1::HigherTruncation { degree: 3 }
        ),
        ExitPathNodeStateV1::Unknown {
            reason: ExitPathUnknownReasonV1::PathEquivalenceDataInsufficient,
        }
    );

    let mut weak_monodromy = finite.clone();
    weak_monodromy.hypotheses.monodromy = MonodromyHypothesisV1::Trivial {
        witness: witness(119),
    };
    assert_eq!(
        state(
            &validate(weak_monodromy),
            ExitPathApproximationV1::HigherTruncation { degree: 3 }
        ),
        ExitPathNodeStateV1::Unknown {
            reason: ExitPathUnknownReasonV1::MonodromyDataInsufficient,
        }
    );

    finite.path_equivalence = StratifiedPathEquivalenceV1::FullHigher {
        relation: StratifiedPathEquivalenceIdV1::from_bytes(bytes(98)),
        coherence: HigherCoherenceArtifactIdV1::from_bytes(bytes(99)),
    };
    finite.hypotheses.links = ConicalLinkHypothesisV1::FullHigher {
        links: ExitLinkCatalogIdV1::from_bytes(bytes(100)),
        coherence: HigherCoherenceArtifactIdV1::from_bytes(bytes(101)),
        witness: witness(102),
    };
    finite.hypotheses.monodromy = MonodromyHypothesisV1::FullHigher {
        groupoids: StratumGroupoidCatalogIdV1::from_bytes(bytes(103)),
        local_systems: LocalSystemCatalogIdV1::from_bytes(bytes(104)),
        coherence: HigherCoherenceArtifactIdV1::from_bytes(bytes(105)),
        witness: witness(106),
    };
    finite.hypotheses.homotopy = HomotopyFidelityV1::FullHigher {
        coherence: HigherCoherenceArtifactIdV1::from_bytes(bytes(107)),
        witness: witness(108),
    };
    assert_eq!(
        state(
            &validate(finite),
            ExitPathApproximationV1::FullHigherCategory
        ),
        ExitPathNodeStateV1::SufficientStatement
    );
}

#[test]
fn hypothesis_deletion_degrades_every_node_to_unknown() {
    let mut deleted = regular_cell_ir();
    deleted.hypotheses.properness = ExitPropernessHypothesisV1::Unknown {
        no_claim: ExitPathNoClaimIdV1::from_bytes(bytes(110)),
    };
    let admitted = validate(deleted);
    assert!(admitted.theorem_lattice().iter().all(|node| {
        node.state
            == ExitPathNodeStateV1::Unknown {
                reason: ExitPathUnknownReasonV1::CommonHypothesisMissing,
            }
    }));
}

#[test]
fn every_snapshot_field_family_affects_identity() {
    let baseline = validate(regular_cell_ir());

    let mut second_subject_ir = admitted_subject_ir();
    second_subject_ir.model_version = DerivedModelVersionIdV1::from_bytes(bytes(125));
    let second_subject = with_cx(false, |cx| {
        admit_derived_geometry_v1(second_subject_ir, DerivedAdmissionBudgetV1::STANDARD, cx)
            .expect("second RD.1a subject admits")
    });
    let mut rebound = regular_cell_ir();
    rebound.geometry = second_subject.id();
    rebound.model_version = second_subject.ir().model_version;
    rebound.stratification = second_subject.ir().stratification.id;
    rebound.frame = second_subject.ir().frame;
    rebound.units = second_subject.ir().unit_system;
    let rebound = with_cx(false, |cx| {
        validate_exit_path_family_v1(rebound, &second_subject, cx).expect("rebound family admits")
    });
    assert_ne!(baseline.snapshot_id(), rebound.snapshot_id());

    let mut entrance = regular_cell_ir();
    entrance.direction = ExitPathDirectionV1::Entrance;
    let entrance = validate(entrance);
    assert_ne!(baseline.snapshot_id(), entrance.snapshot_id());

    let mut cosheaf = regular_cell_ir();
    cosheaf.variance = ConstructibleVarianceV1::CosheafCovariant;
    assert_ne!(baseline.snapshot_id(), validate(cosheaf).snapshot_id());

    let mut complex = regular_cell_ir();
    complex.constructible_coefficients = CoefficientSystemV1::RationalComplex;
    assert_ne!(baseline.snapshot_id(), validate(complex).snapshot_id());

    let mut refined = regular_cell_ir();
    refined.hypotheses.refinement = RefinementHypothesisV1::CommonRefinement {
        refinement: StratifiedRefinementIdV1::from_bytes(bytes(111)),
        forward: RefinementMapIdV1::from_bytes(bytes(112)),
        reverse: RefinementMapIdV1::from_bytes(bytes(113)),
        witness: witness(114),
    };
    assert_ne!(baseline.snapshot_id(), validate(refined).snapshot_id());

    let mut path = regular_cell_ir();
    path.path_equivalence = StratifiedPathEquivalenceV1::EndpointFixed {
        relation: StratifiedPathEquivalenceIdV1::from_bytes(bytes(115)),
    };
    assert_ne!(baseline.snapshot_id(), validate(path).snapshot_id());

    let mut link = regular_cell_ir();
    link.hypotheses.links = ConicalLinkHypothesisV1::Contractible {
        links: ExitLinkCatalogIdV1::from_bytes(bytes(8)),
        witness: witness(116),
    };
    assert_ne!(baseline.snapshot_id(), validate(link).snapshot_id());

    let mut monodromy = regular_cell_ir();
    monodromy.hypotheses.monodromy = MonodromyHypothesisV1::Trivial {
        witness: witness(117),
    };
    assert_ne!(baseline.snapshot_id(), validate(monodromy).snapshot_id());

    let mut lifecycle = regular_cell_ir();
    lifecycle.theorem_state = ExitPathTheoremStateV1::Candidate {
        card: lifecycle.tcb.theorem_card,
        witness: witness(118),
    };
    assert_ne!(baseline.snapshot_id(), validate(lifecycle).snapshot_id());

    let mut tcb = regular_cell_ir();
    tcb.tcb.checker = ExitPathCheckerIdV1::from_bytes(bytes(119));
    assert_ne!(baseline.snapshot_id(), validate(tcb).snapshot_id());

    let mut budget = regular_cell_ir();
    budget.budget.max_implication_checks += 1;
    assert_ne!(baseline.snapshot_id(), validate(budget).snapshot_id());

    let mut truncation = regular_cell_ir();
    truncation.budget.max_truncation = 2;
    truncation.budget.max_implication_checks = 5;
    assert_ne!(baseline.snapshot_id(), validate(truncation).snapshot_id());

    let mut falsifier = regular_cell_ir();
    falsifier.falsifiers[0].witness = witness(120);
    assert_ne!(baseline.snapshot_id(), validate(falsifier).snapshot_id());
}

#[test]
fn falsifier_set_is_canonical_and_mandatory() {
    let first = validate(regular_cell_ir());
    let mut reordered = regular_cell_ir();
    reordered.falsifiers.reverse();
    assert_eq!(
        first.identity_receipt(),
        validate(reordered).identity_receipt()
    );

    let mut missing = regular_cell_ir();
    missing
        .falsifiers
        .retain(|item| item.kind != ExitPathFalsifierKindV1::SameIncidenceDifferentMonodromy);
    assert_issue(
        missing,
        ExitPathSemanticIssueV1::MissingFalsifier {
            kind: ExitPathFalsifierKindV1::SameIncidenceDifferentMonodromy,
        },
    );

    let mut degenerate = regular_cell_ir();
    degenerate.falsifiers[0].right = degenerate.falsifiers[0].left;
    assert_issue(degenerate, ExitPathSemanticIssueV1::DegenerateFalsifier);
}

#[test]
fn theorem_state_tcb_and_refutation_references_fail_closed() {
    let mut card = regular_cell_ir();
    card.tcb.theorem_card = ExitPathTheoremCardIdV1::from_bytes(bytes(120));
    assert_issue(card, ExitPathSemanticIssueV1::TheoremCardMismatch);

    let mut refuted = regular_cell_ir();
    refuted.theorem_state = ExitPathTheoremStateV1::RefutationRecorded {
        approximation: ExitPathApproximationV1::IncidencePoset,
        falsifier: ExitPathFalsifierIdV1::from_bytes(bytes(121)),
    };
    assert_issue(refuted, ExitPathSemanticIssueV1::UnknownRecordedFalsifier);

    let mut absent_node = regular_cell_ir();
    absent_node.theorem_state = ExitPathTheoremStateV1::RefutationRecorded {
        approximation: ExitPathApproximationV1::HigherTruncation { degree: 9 },
        falsifier: absent_node.falsifiers[0].id,
    };
    assert_issue(
        absent_node,
        ExitPathSemanticIssueV1::UnknownRecordedApproximation,
    );

    let mut retained_refutation = regular_cell_ir();
    retained_refutation.hypotheses.homotopy = HomotopyFidelityV1::RetainedThrough {
        degree: 1,
        coherence: HigherCoherenceArtifactIdV1::from_bytes(bytes(122)),
        witness: witness(123),
    };
    let falsifier = retained_refutation.falsifiers[0].id;
    retained_refutation.theorem_state = ExitPathTheoremStateV1::RefutationRecorded {
        approximation: ExitPathApproximationV1::IncidencePoset,
        falsifier,
    };
    let admitted = validate(retained_refutation);
    assert_eq!(
        state(&admitted, ExitPathApproximationV1::IncidencePoset),
        ExitPathNodeStateV1::RefutationRecorded { falsifier }
    );
    assert_eq!(
        state(
            &admitted,
            ExitPathApproximationV1::StratumGroupoidEnrichedExitCategory,
        ),
        ExitPathNodeStateV1::SufficientStatement
    );
}

fn with_schema(ir: ExitPathFamilyIrV1, schema_version: u32) -> ExitPathFamilyIrV1 {
    ExitPathFamilyIrV1::with_schema_version(
        schema_version,
        ir.geometry,
        ir.model_version,
        ir.stratification,
        ir.frame,
        ir.units,
        ir.space_class,
        ir.direction,
        ir.variance,
        ir.convention,
        ir.path_equivalence,
        ir.constructible_coefficients,
        ir.hypotheses,
        ir.falsifiers,
        ir.theorem_state,
        ir.tcb,
        ir.budget,
    )
}

#[test]
fn schema_degrees_budgets_caps_and_cancellation_refuse() {
    assert_issue(
        with_schema(regular_cell_ir(), EXIT_PATH_SCHEMA_VERSION_V1 + 1),
        ExitPathSemanticIssueV1::UnsupportedSchemaVersion {
            found: EXIT_PATH_SCHEMA_VERSION_V1 + 1,
            supported: EXIT_PATH_SCHEMA_VERSION_V1,
        },
    );

    let mut mismatched_subject = regular_cell_ir();
    mismatched_subject.model_version = DerivedModelVersionIdV1::from_bytes(bytes(124));
    assert_issue(
        mismatched_subject,
        ExitPathSemanticIssueV1::SubjectBindingMismatch {
            field: ExitPathSubjectBindingFieldV1::ModelVersion,
        },
    );

    let mut missing_identity = regular_cell_ir();
    missing_identity.hypotheses.links = ConicalLinkHypothesisV1::Contractible {
        links: ExitLinkCatalogIdV1::from_bytes([0; 32]),
        witness: witness(9),
    };
    assert_issue(
        missing_identity,
        ExitPathSemanticIssueV1::MissingIdentity {
            field: ExitPathIdentityFieldV1::LinkData,
        },
    );

    let mut degree = regular_cell_ir();
    degree.budget.max_truncation = MAX_EXIT_PATH_TRUNCATION_V1 + 1;
    assert_issue(
        degree,
        ExitPathSemanticIssueV1::InvalidValue {
            field: ExitPathFieldV1::MaxTruncation,
        },
    );

    let mut budget = regular_cell_ir();
    budget.budget.max_implication_checks = 1;
    assert_issue(
        budget,
        ExitPathSemanticIssueV1::InvalidValue {
            field: ExitPathFieldV1::ImplicationBudget,
        },
    );

    let mut exact_artifacts = regular_cell_ir();
    let required = exact_artifacts.required_referenced_artifact_slots();
    exact_artifacts.budget.max_referenced_artifact_slots = required;
    validate(exact_artifacts);

    let mut insufficient_artifacts = regular_cell_ir();
    insufficient_artifacts.budget.max_referenced_artifact_slots = required - 1;
    assert_issue(
        insufficient_artifacts,
        ExitPathSemanticIssueV1::ReferencedArtifactBudgetExceeded {
            required,
            available: required - 1,
        },
    );

    let mut capped = regular_cell_ir();
    let sample = capped.falsifiers[0];
    capped.falsifiers = vec![sample; MAX_EXIT_PATH_FALSIFIERS_V1 + 1];
    assert_issue(
        capped,
        ExitPathSemanticIssueV1::TooManyFalsifiers {
            found: MAX_EXIT_PATH_FALSIFIERS_V1 + 1,
            limit: MAX_EXIT_PATH_FALSIFIERS_V1,
        },
    );

    let report = with_cx(true, |cx| {
        validate_exit_path_family_v1(regular_cell_ir(), admitted_subject(), cx)
            .expect_err("cancelled")
    });
    assert_eq!(report.issues(), &[ExitPathSemanticIssueV1::Cancelled]);
}
