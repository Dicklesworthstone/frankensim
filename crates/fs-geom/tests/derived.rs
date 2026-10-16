//! RD.1a G0/G3/G5 fixtures for finite derived/stratified object admission.

#![allow(clippy::too_many_lines, clippy::wildcard_imports)]

use asupersync::types::Budget;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::derived::*;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    with_gate_cx(&gate, f)
}

fn with_gate_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 0x5244_3161,
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

fn subject(seed: u8) -> DerivedSubjectIdV1 {
    DerivedSubjectIdV1::from_bytes([seed; 32])
}

fn model_version(seed: u8) -> DerivedModelVersionIdV1 {
    DerivedModelVersionIdV1::from_bytes([seed; 32])
}

fn chart_id(seed: u8) -> ConfigurationChartIdV1 {
    ConfigurationChartIdV1::from_bytes([seed; 32])
}

fn frame(seed: u8) -> DerivedFrameIdV1 {
    DerivedFrameIdV1::from_bytes([seed; 32])
}

fn unit_system(seed: u8) -> DerivedUnitSystemIdV1 {
    DerivedUnitSystemIdV1::from_bytes([seed; 32])
}

fn quantity(seed: u8) -> DerivedQuantityKindIdV1 {
    DerivedQuantityKindIdV1::from_bytes([seed; 32])
}

fn polynomial(seed: u8) -> PolynomialIdV1 {
    PolynomialIdV1::from_bytes([seed; 32])
}

fn analytic(seed: u8) -> AnalyticProgramIdV1 {
    AnalyticProgramIdV1::from_bytes([seed; 32])
}

fn witness(seed: u8) -> DerivedWitnessIdV1 {
    DerivedWitnessIdV1::from_bytes([seed; 32])
}

fn no_claim(seed: u8) -> DerivedNoClaimIdV1 {
    DerivedNoClaimIdV1::from_bytes([seed; 32])
}

fn equality_id(seed: u8) -> EqualityConstraintIdV1 {
    EqualityConstraintIdV1::from_bytes([seed; 32])
}

fn inequality_id(seed: u8) -> InequalityConstraintIdV1 {
    InequalityConstraintIdV1::from_bytes([seed; 32])
}

fn boundary_id(seed: u8) -> RelativeBoundaryIdV1 {
    RelativeBoundaryIdV1::from_bytes([seed; 32])
}

fn contact_id(seed: u8) -> ContactConstraintIdV1 {
    ContactConstraintIdV1::from_bytes([seed; 32])
}

fn constitutive_id(seed: u8) -> ConstitutiveDatumIdV1 {
    ConstitutiveDatumIdV1::from_bytes([seed; 32])
}

fn complex_id(seed: u8) -> DerivedComplexIdV1 {
    DerivedComplexIdV1::from_bytes([seed; 32])
}

fn map(seed: u8) -> DerivedLinearMapIdV1 {
    DerivedLinearMapIdV1::from_bytes([seed; 32])
}

fn resolution(seed: u8) -> DerivedResolutionIdV1 {
    DerivedResolutionIdV1::from_bytes([seed; 32])
}

fn local_model_id(seed: u8) -> DerivedLocalModelIdV1 {
    DerivedLocalModelIdV1::from_bytes([seed; 32])
}

fn stratification_id(seed: u8) -> StratificationIdV1 {
    StratificationIdV1::from_bytes([seed; 32])
}

fn stratum_id(seed: u8) -> StratumIdV1 {
    StratumIdV1::from_bytes([seed; 32])
}

fn link_id(seed: u8) -> LocalLinkIdV1 {
    LocalLinkIdV1::from_bytes([seed; 32])
}

fn units(quantity_seed: u8) -> UnitBindingV1 {
    UnitBindingV1 {
        system: unit_system(3),
        quantity: quantity(quantity_seed),
        scale_to_canonical: 1.0,
    }
}

fn exact(seed: u8) -> FiniteComputabilityV1 {
    FiniteComputabilityV1::ExactFinite {
        kernel: witness(seed),
    }
}

fn polynomial_encoding(seed: u8) -> LocalFunctionEncodingV1 {
    LocalFunctionEncodingV1::Polynomial {
        polynomial: polynomial(seed),
        variables: 2,
        degree: 2,
    }
}

fn finite_complex(seed: u8, role: DerivedComplexRoleV1) -> FiniteDerivedComplexV1 {
    FiniteDerivedComplexV1 {
        id: complex_id(seed),
        chart: chart_id(4),
        role,
        spaces: vec![
            GradedSpaceV1 {
                degree: 0,
                dimension: 2,
                quantity: quantity(seed.wrapping_add(1)),
            },
            GradedSpaceV1 {
                degree: 1,
                dimension: 1,
                quantity: quantity(seed.wrapping_add(2)),
            },
        ],
        differentials: vec![ComplexDifferentialV1 {
            from_degree: 0,
            to_degree: 1,
            map: map(seed.wrapping_add(3)),
            square_zero_witness: witness(seed.wrapping_add(4)),
        }],
        resolution: FiniteResolutionV1 {
            id: resolution(seed.wrapping_add(5)),
            min_degree: 0,
            max_degree: 1,
            max_basis_dimension: 2,
            truncation_order: 0,
            remainder: None,
        },
        computability: exact(seed.wrapping_add(6)),
    }
}

fn regular_linkage() -> DerivedGeometryIrV1 {
    let chart = chart_id(4);
    let equality = equality_id(20);
    let local_model = local_model_id(40);
    let stratum = stratum_id(50);
    DerivedGeometryIrV1 {
        schema_version: DERIVED_GEOMETRY_SCHEMA_VERSION_V1,
        subject: subject(1),
        model_version: model_version(2),
        category: GeometricCategoryV1::Semialgebraic,
        coefficients: CoefficientSystemV1::RationalReal,
        frame: frame(5),
        unit_system: unit_system(3),
        locality: LocalityScopeV1::GermAt {
            chart,
            point: witness(6),
        },
        compactness: CompactnessV1::RelativelyCompact {
            witness: witness(7),
        },
        charts: vec![ConfigurationChartV1 {
            id: chart,
            class: ConfigurationChartClassV1::Semialgebraic,
            coordinate_dimension: 2,
            ambient_dimension: 2,
            frame: frame(5),
            coordinates: units(8),
            locality: LocalityScopeV1::GermAt {
                chart,
                point: witness(6),
            },
            compactness: CompactnessV1::RelativelyCompact {
                witness: witness(7),
            },
            regularity: RegularityClassV1::Polynomial,
            computability: exact(9),
        }],
        equalities: vec![EqualityConstraintGermV1 {
            id: equality,
            chart,
            codomain_dimension: 1,
            equation: polynomial_encoding(21),
            regularity: RegularityClassV1::Polynomial,
            units: units(10),
            computability: exact(22),
        }],
        inequalities: Vec::new(),
        boundaries: Vec::new(),
        contacts: Vec::new(),
        constitutive_data: Vec::new(),
        complexes: vec![
            finite_complex(30, DerivedComplexRoleV1::Tangent),
            finite_complex(31, DerivedComplexRoleV1::Cotangent),
            finite_complex(32, DerivedComplexRoleV1::DeformationObstruction),
        ],
        local_models: vec![DerivedLocalModelV1 {
            id: local_model,
            chart,
            class: DerivedLocalModelClassV1::RegularCompleteIntersection,
            equalities: vec![equality],
            active_inequalities: Vec::new(),
            active_contacts: Vec::new(),
            constitutive_data: Vec::new(),
            tangent_complex: complex_id(30),
            cotangent_complex: complex_id(31),
            deformation_complex: complex_id(32),
            virtual_dimension: 1,
            locality: LocalityScopeV1::GermAt {
                chart,
                point: witness(6),
            },
            presentation: PresentationScopeV1::Literal {
                no_claim: no_claim(41),
            },
        }],
        stratification: StratificationV1 {
            id: stratification_id(60),
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
                    witness: witness(51),
                },
            }],
            incidences: Vec::new(),
            local_links: Vec::new(),
        },
        proof_state: DerivedProofStateV1::StructuralNoClaim {
            no_claim: no_claim(61),
        },
    }
}

fn contact_corner() -> DerivedGeometryIrV1 {
    let mut ir = regular_linkage();
    let chart = chart_id(4);
    let boundary_stratum = stratum_id(50);
    let parent_stratum = stratum_id(70);
    let parent_model = local_model_id(71);
    let inequality = inequality_id(72);
    let contact = contact_id(73);
    let constitutive = constitutive_id(74);

    ir.inequalities.push(InequalityConstraintGermV1 {
        id: inequality,
        chart,
        sense: InequalitySenseV1::NonNegative,
        function: polynomial_encoding(75),
        state: ActiveSetStateV1::Active {
            witness: witness(76),
        },
        normal_cone: NormalConeClassV1::Ray,
        units: units(77),
        computability: exact(78),
    });
    ir.boundaries.extend([
        RelativeBoundaryV1 {
            id: boundary_id(79),
            chart,
            parent: parent_stratum,
            boundary: boundary_stratum,
            orientation: BoundaryOrientationV1::Outward,
            witness: witness(80),
            units: units(81),
        },
        RelativeBoundaryV1 {
            id: boundary_id(82),
            chart,
            parent: parent_stratum,
            boundary: boundary_stratum,
            orientation: BoundaryOrientationV1::Inward,
            witness: witness(83),
            units: units(84),
        },
    ]);
    ir.contacts.push(ContactConstraintV1 {
        id: contact,
        chart,
        side_a: boundary_id(79),
        side_b: boundary_id(82),
        gap: polynomial_encoding(85),
        state: ActiveSetStateV1::Active {
            witness: witness(86),
        },
        normal_cone: NormalConeClassV1::Polyhedral { generators: 2 },
        law: ContactLawV1::Coulomb {
            friction_coefficient: 0.3,
        },
        units: units(87),
        computability: exact(88),
    });
    ir.constitutive_data.push(ConstitutiveDatumV1 {
        id: constitutive,
        chart,
        role: ConstitutiveRoleV1::Dissipation,
        state_dimension: 1,
        law: polynomial_encoding(89),
        units: units(90),
        computability: exact(91),
    });

    let boundary_model = &mut ir.local_models[0];
    boundary_model.class = DerivedLocalModelClassV1::ContactCorner;
    boundary_model.active_inequalities.push(inequality);
    boundary_model.active_contacts.push(contact);
    boundary_model.constitutive_data.push(constitutive);
    ir.local_models.push(DerivedLocalModelV1 {
        id: parent_model,
        chart,
        class: DerivedLocalModelClassV1::GeneralFiniteDerived,
        equalities: Vec::new(),
        active_inequalities: Vec::new(),
        active_contacts: Vec::new(),
        constitutive_data: Vec::new(),
        tangent_complex: complex_id(30),
        cotangent_complex: complex_id(31),
        deformation_complex: complex_id(32),
        virtual_dimension: 2,
        locality: LocalityScopeV1::CompactNeighborhood {
            chart,
            witness: witness(92),
        },
        presentation: PresentationScopeV1::Literal {
            no_claim: no_claim(93),
        },
    });

    let stratum = &mut ir.stratification.strata[0];
    stratum.active_inequalities.push(inequality);
    stratum.active_contacts.push(contact);
    stratum.relative_boundary = Some(boundary_id(79));
    ir.stratification.strata.push(StratumSpecV1 {
        id: parent_stratum,
        chart,
        local_model: parent_model,
        dimension: 2,
        active_inequalities: Vec::new(),
        active_contacts: Vec::new(),
        relative_boundary: None,
        regularity: RegularityClassV1::Polynomial,
        compactness: CompactnessV1::RelativelyCompact {
            witness: witness(94),
        },
    });
    ir.stratification.incidences.push(StratumIncidenceV1 {
        lower: boundary_stratum,
        upper: parent_stratum,
        codimension: 1,
        witness: witness(95),
    });
    ir.stratification.local_links.push(LocalLinkV1 {
        id: link_id(96),
        stratum: boundary_stratum,
        ambient_stratum: parent_stratum,
        dimension: 0,
        compactness_witness: witness(97),
        topology: LocalLinkTopologyV1::FiniteComplex {
            resolution: resolution(98),
            witness: witness(99),
        },
    });
    ir
}

fn add_second_chart(ir: &mut DerivedGeometryIrV1, seed: u8) -> ConfigurationChartIdV1 {
    let id = chart_id(seed);
    let mut chart = ir.charts[0].clone();
    chart.id = id;
    chart.locality = LocalityScopeV1::GermAt {
        chart: id,
        point: witness(seed.wrapping_add(1)),
    };
    ir.charts.push(chart);
    id
}

fn assert_issue(
    result: Result<AdmittedDerivedGeometryV1, DerivedAdmissionReportV1>,
    expected: impl Fn(&DerivedAdmissionIssueV1) -> bool,
) {
    let report = result.expect_err("adversarial object must refuse");
    assert!(
        report.issues().iter().any(expected),
        "expected issue missing from {:?}",
        report.issues()
    );
}

#[test]
fn regular_linkage_replays_canonically_and_binds_versions() {
    with_cx(|cx| {
        let source = regular_linkage();
        let admitted =
            admit_derived_geometry_v1(source.clone(), DerivedAdmissionBudgetV1::STANDARD, cx)
                .expect("regular finite linkage must admit");
        let replay = admit_derived_geometry_v1(
            admitted.ir().clone(),
            DerivedAdmissionBudgetV1::STANDARD,
            cx,
        )
        .expect("canonical schema replay must admit");
        assert_eq!(admitted.admission_receipt(), replay.admission_receipt());
        assert_eq!(admitted.id(), replay.id());

        let mut reordered = source.clone();
        reordered.complexes.reverse();
        let reordered =
            admit_derived_geometry_v1(reordered, DerivedAdmissionBudgetV1::STANDARD, cx)
                .expect("input order is not semantic");
        assert_eq!(admitted.id(), reordered.id());
        assert_eq!(
            admitted.admission_receipt().canonical_preimage(),
            reordered.admission_receipt().canonical_preimage()
        );

        let mut changed_version = source;
        changed_version.model_version = model_version(111);
        let changed_version =
            admit_derived_geometry_v1(changed_version, DerivedAdmissionBudgetV1::STANDARD, cx)
                .expect("new immutable version remains structurally admissible");
        assert_ne!(admitted.id(), changed_version.id());
        assert_eq!(
            admitted.admission_receipt().limits().max_canonical_bytes(),
            DerivedAdmissionBudgetV1::STANDARD.max_canonical_bytes
        );
    });
}

#[test]
fn redundant_cusp_and_node_presentations_are_distinct_and_deterministic() {
    with_cx(|cx| {
        let mut redundant = regular_linkage();
        redundant.equalities.push(EqualityConstraintGermV1 {
            id: equality_id(19),
            chart: chart_id(4),
            codomain_dimension: 1,
            equation: polynomial_encoding(18),
            regularity: RegularityClassV1::Polynomial,
            units: units(17),
            computability: exact(16),
        });
        redundant.local_models[0].class = DerivedLocalModelClassV1::RedundantPresentation;
        redundant.local_models[0].equalities = vec![equality_id(20), equality_id(19)];

        let first =
            admit_derived_geometry_v1(redundant.clone(), DerivedAdmissionBudgetV1::STANDARD, cx)
                .expect("redundant equations remain an explicit admitted presentation");
        redundant.equalities.reverse();
        redundant.local_models[0].equalities.reverse();
        let reordered =
            admit_derived_geometry_v1(redundant.clone(), DerivedAdmissionBudgetV1::STANDARD, cx)
                .expect("canonical singular presentation order must replay");
        assert_eq!(first.id(), reordered.id());

        let mut cusp = redundant.clone();
        cusp.local_models[0].class = DerivedLocalModelClassV1::Cusp;
        let cusp = admit_derived_geometry_v1(cusp, DerivedAdmissionBudgetV1::STANDARD, cx)
            .expect("finite polynomial cusp local model must admit structurally");
        let mut node = redundant;
        node.local_models[0].class = DerivedLocalModelClassV1::Node;
        let node = admit_derived_geometry_v1(node, DerivedAdmissionBudgetV1::STANDARD, cx)
            .expect("finite polynomial node local model must admit structurally");
        assert_ne!(first.id(), cusp.id());
        assert_ne!(cusp.id(), node.id());
    });
}

#[test]
fn boundary_contact_corner_keeps_roles_units_and_active_sets_distinct() {
    with_cx(|cx| {
        let corner = contact_corner();
        let admitted =
            admit_derived_geometry_v1(corner.clone(), DerivedAdmissionBudgetV1::STANDARD, cx)
                .expect("finite contact corner and compact link must admit");
        assert_eq!(admitted.ir().boundaries.len(), 2);
        assert_eq!(admitted.ir().contacts.len(), 1);
        assert_eq!(admitted.ir().constitutive_data.len(), 1);

        let mut mixed_units = corner.clone();
        mixed_units.inequalities[0].units.system = unit_system(120);
        assert_issue(
            admit_derived_geometry_v1(mixed_units, DerivedAdmissionBudgetV1::STANDARD, cx),
            |issue| {
                matches!(
                    issue,
                    DerivedAdmissionIssueV1::MixedUnitSystem {
                        kind: DerivedObjectKindV1::Inequality
                    }
                )
            },
        );

        let mut candidate = corner.clone();
        candidate.inequalities[0].state = ActiveSetStateV1::Candidate {
            no_claim: no_claim(121),
        };
        assert_issue(
            admit_derived_geometry_v1(candidate, DerivedAdmissionBudgetV1::STANDARD, cx),
            |issue| {
                matches!(
                    issue,
                    DerivedAdmissionIssueV1::ActiveSetMismatch {
                        kind: DerivedObjectKindV1::Inequality
                    }
                )
            },
        );

        let mut negative_zero = corner.clone();
        negative_zero.contacts[0].law = ContactLawV1::Coulomb {
            friction_coefficient: -0.0,
        };
        assert_issue(
            admit_derived_geometry_v1(negative_zero, DerivedAdmissionBudgetV1::STANDARD, cx),
            |issue| {
                matches!(
                    issue,
                    DerivedAdmissionIssueV1::InvalidDimension {
                        kind: DerivedObjectKindV1::Contact,
                        field: "friction_coefficient"
                    }
                )
            },
        );

        let mut negative_zero_unit = corner.clone();
        negative_zero_unit.charts[0].coordinates.scale_to_canonical = -0.0;
        assert_issue(
            admit_derived_geometry_v1(negative_zero_unit, DerivedAdmissionBudgetV1::STANDARD, cx),
            |issue| {
                matches!(
                    issue,
                    DerivedAdmissionIssueV1::InvalidUnitScale {
                        kind: DerivedObjectKindV1::Chart
                    }
                )
            },
        );

        let mut complex_contact = corner;
        complex_contact.category = GeometricCategoryV1::Algebraic;
        complex_contact.coefficients = CoefficientSystemV1::RationalComplex;
        complex_contact.charts[0].class = ConfigurationChartClassV1::Algebraic;
        assert_issue(
            admit_derived_geometry_v1(complex_contact, DerivedAdmissionBudgetV1::STANDARD, cx),
            |issue| matches!(issue, DerivedAdmissionIssueV1::OrderedSemanticsRequiresReal),
        );
    });
}

#[test]
fn real_complex_category_changes_identity_without_laundering_order() {
    with_cx(|cx| {
        let real =
            admit_derived_geometry_v1(regular_linkage(), DerivedAdmissionBudgetV1::STANDARD, cx)
                .unwrap();
        let mut complex = regular_linkage();
        complex.category = GeometricCategoryV1::Algebraic;
        complex.coefficients = CoefficientSystemV1::RationalComplex;
        complex.charts[0].class = ConfigurationChartClassV1::Algebraic;
        let complex = admit_derived_geometry_v1(complex, DerivedAdmissionBudgetV1::STANDARD, cx)
            .expect("finite algebraic complex object without ordered constraints admits");
        assert_ne!(real.id(), complex.id());

        let mut invalid = regular_linkage();
        invalid.coefficients = CoefficientSystemV1::RationalComplex;
        assert_issue(
            admit_derived_geometry_v1(invalid, DerivedAdmissionBudgetV1::STANDARD, cx),
            |issue| matches!(issue, DerivedAdmissionIssueV1::CategoryCoefficientMismatch),
        );
    });
}

#[test]
fn restricted_analytic_program_admits_but_opaque_encoding_refuses() {
    with_cx(|cx| {
        let mut analytic_ir = regular_linkage();
        analytic_ir.category = GeometricCategoryV1::RestrictedAnalytic;
        analytic_ir.charts[0].class = ConfigurationChartClassV1::RestrictedAnalytic;
        analytic_ir.charts[0].regularity = RegularityClassV1::Analytic;
        analytic_ir.equalities[0].equation = LocalFunctionEncodingV1::RestrictedAnalytic {
            program: analytic(130),
            primitives: witness(131),
            derivative_order: 3,
        };
        analytic_ir.equalities[0].regularity = RegularityClassV1::Analytic;
        let admitted =
            admit_derived_geometry_v1(analytic_ir.clone(), DerivedAdmissionBudgetV1::STANDARD, cx)
                .expect("restricted analytic program with finite primitives must admit");
        assert_eq!(
            admitted.ir().category,
            GeometricCategoryV1::RestrictedAnalytic
        );

        let mut opaque = analytic_ir;
        opaque.equalities[0].equation = LocalFunctionEncodingV1::OpaqueExternal {
            no_claim: no_claim(132),
        };
        opaque.equalities[0].computability = FiniteComputabilityV1::ExternalOpaque {
            no_claim: no_claim(133),
        };
        assert_issue(
            admit_derived_geometry_v1(opaque, DerivedAdmissionBudgetV1::STANDARD, cx),
            |issue| {
                matches!(
                    issue,
                    DerivedAdmissionIssueV1::UnsupportedFunctionEncoding {
                        kind: DerivedObjectKindV1::Equality
                    }
                )
            },
        );
    });
}

#[test]
fn unbounded_infinite_or_malformed_local_models_refuse() {
    with_cx(|cx| {
        let mut unbounded = regular_linkage();
        unbounded.locality = LocalityScopeV1::GlobalUnbounded;
        unbounded.compactness = CompactnessV1::Unbounded;
        unbounded.charts[0].locality = LocalityScopeV1::GlobalUnbounded;
        unbounded.charts[0].compactness = CompactnessV1::Unbounded;
        assert_issue(
            admit_derived_geometry_v1(unbounded, DerivedAdmissionBudgetV1::STANDARD, cx),
            |issue| {
                matches!(
                    issue,
                    DerivedAdmissionIssueV1::UnsupportedLocality {
                        kind: DerivedObjectKindV1::Global
                    }
                )
            },
        );

        let mut infinite = regular_linkage();
        infinite.complexes[0].computability = FiniteComputabilityV1::Infinite;
        assert_issue(
            admit_derived_geometry_v1(infinite, DerivedAdmissionBudgetV1::STANDARD, cx),
            |issue| {
                matches!(
                    issue,
                    DerivedAdmissionIssueV1::UnsupportedComputability {
                        kind: DerivedObjectKindV1::Complex
                    }
                )
            },
        );

        let mut wrong_role = regular_linkage();
        wrong_role.complexes[0].role = DerivedComplexRoleV1::Cotangent;
        assert_issue(
            admit_derived_geometry_v1(wrong_role, DerivedAdmissionBudgetV1::STANDARD, cx),
            |issue| {
                matches!(
                    issue,
                    DerivedAdmissionIssueV1::ComplexRoleMismatch {
                        field: "tangent_complex"
                    }
                )
            },
        );

        let mut sparse_degrees = regular_linkage();
        sparse_degrees.complexes[0].spaces[1].degree = 2;
        sparse_degrees.complexes[0].differentials[0].to_degree = 2;
        sparse_degrees.complexes[0].resolution.max_degree = 2;
        assert_issue(
            admit_derived_geometry_v1(sparse_degrees, DerivedAdmissionBudgetV1::STANDARD, cx),
            |issue| {
                matches!(
                    issue,
                    DerivedAdmissionIssueV1::InvalidComplex {
                        field: "noncontiguous_degree"
                    }
                )
            },
        );

        let mut oversized_resolution = regular_linkage();
        oversized_resolution.complexes[0]
            .resolution
            .max_basis_dimension = DERIVED_GEOMETRY_HARD_MAX_DIMENSION_V1 + 1;
        assert_issue(
            admit_derived_geometry_v1(oversized_resolution, DerivedAdmissionBudgetV1::STANDARD, cx),
            |issue| {
                matches!(
                    issue,
                    DerivedAdmissionIssueV1::InvalidComplex {
                        field: "resolution_bounds"
                    }
                )
            },
        );
    });
}

#[test]
fn cross_chart_boundaries_incidences_and_constitutive_refs_refuse() {
    with_cx(|cx| {
        let mut chart_locality = regular_linkage();
        let other_chart = add_second_chart(&mut chart_locality, 147);
        chart_locality.charts[0].locality = LocalityScopeV1::GermAt {
            chart: other_chart,
            point: witness(148),
        };
        assert_issue(
            admit_derived_geometry_v1(chart_locality, DerivedAdmissionBudgetV1::STANDARD, cx),
            |issue| {
                matches!(
                    issue,
                    DerivedAdmissionIssueV1::CrossChartReference {
                        kind: DerivedObjectKindV1::Chart,
                        field: "self_locality_chart"
                    }
                )
            },
        );

        let mut equality = regular_linkage();
        let other_chart = add_second_chart(&mut equality, 149);
        equality.equalities[0].chart = other_chart;
        assert_issue(
            admit_derived_geometry_v1(equality, DerivedAdmissionBudgetV1::STANDARD, cx),
            |issue| {
                matches!(
                    issue,
                    DerivedAdmissionIssueV1::CrossChartReference {
                        kind: DerivedObjectKindV1::LocalModel,
                        field: "equality"
                    }
                )
            },
        );

        let mut boundary = contact_corner();
        let other_chart = add_second_chart(&mut boundary, 150);
        boundary.boundaries[0].chart = other_chart;
        assert_issue(
            admit_derived_geometry_v1(boundary, DerivedAdmissionBudgetV1::STANDARD, cx),
            |issue| {
                matches!(
                    issue,
                    DerivedAdmissionIssueV1::CrossChartReference {
                        kind: DerivedObjectKindV1::Boundary,
                        field: "parent" | "boundary"
                    }
                )
            },
        );

        let mut incidence = contact_corner();
        let other_chart = add_second_chart(&mut incidence, 152);
        incidence.stratification.strata[1].chart = other_chart;
        assert_issue(
            admit_derived_geometry_v1(incidence, DerivedAdmissionBudgetV1::STANDARD, cx),
            |issue| {
                matches!(
                    issue,
                    DerivedAdmissionIssueV1::CrossChartReference {
                        kind: DerivedObjectKindV1::Incidence,
                        field: "lower_upper"
                    }
                )
            },
        );

        let mut constitutive = contact_corner();
        let other_chart = add_second_chart(&mut constitutive, 154);
        constitutive.constitutive_data[0].chart = other_chart;
        assert_issue(
            admit_derived_geometry_v1(constitutive, DerivedAdmissionBudgetV1::STANDARD, cx),
            |issue| {
                matches!(
                    issue,
                    DerivedAdmissionIssueV1::CrossChartReference {
                        kind: DerivedObjectKindV1::LocalModel,
                        field: "constitutive_data"
                    }
                )
            },
        );

        let mut active_inequality = contact_corner();
        let other_chart = add_second_chart(&mut active_inequality, 156);
        active_inequality.inequalities[0].chart = other_chart;
        assert_issue(
            admit_derived_geometry_v1(active_inequality, DerivedAdmissionBudgetV1::STANDARD, cx),
            |issue| {
                matches!(
                    issue,
                    DerivedAdmissionIssueV1::CrossChartReference {
                        kind: DerivedObjectKindV1::LocalModel,
                        field: "active_inequality"
                    }
                )
            },
        );

        let mut active_contact = contact_corner();
        let other_chart = add_second_chart(&mut active_contact, 157);
        active_contact.contacts[0].chart = other_chart;
        assert_issue(
            admit_derived_geometry_v1(active_contact, DerivedAdmissionBudgetV1::STANDARD, cx),
            |issue| {
                matches!(
                    issue,
                    DerivedAdmissionIssueV1::CrossChartReference {
                        kind: DerivedObjectKindV1::LocalModel,
                        field: "active_contact"
                    }
                )
            },
        );

        let mut tangent_complex = regular_linkage();
        let other_chart = add_second_chart(&mut tangent_complex, 158);
        tangent_complex.complexes[0].chart = other_chart;
        assert_issue(
            admit_derived_geometry_v1(tangent_complex, DerivedAdmissionBudgetV1::STANDARD, cx),
            |issue| {
                matches!(
                    issue,
                    DerivedAdmissionIssueV1::CrossChartReference {
                        kind: DerivedObjectKindV1::LocalModel,
                        field: "tangent_complex"
                    }
                )
            },
        );

        let mut stratum_model = regular_linkage();
        let other_chart = add_second_chart(&mut stratum_model, 159);
        stratum_model.stratification.strata[0].chart = other_chart;
        assert_issue(
            admit_derived_geometry_v1(stratum_model, DerivedAdmissionBudgetV1::STANDARD, cx),
            |issue| {
                matches!(
                    issue,
                    DerivedAdmissionIssueV1::CrossChartReference {
                        kind: DerivedObjectKindV1::Stratum,
                        field: "local_model"
                    }
                )
            },
        );
    });
}

#[test]
fn schema_budget_cancellation_and_proof_scope_fail_closed() {
    with_cx(|cx| {
        let mut wrong_schema = regular_linkage();
        wrong_schema.schema_version = 2;
        assert_issue(
            admit_derived_geometry_v1(wrong_schema, DerivedAdmissionBudgetV1::STANDARD, cx),
            |issue| {
                matches!(
                    issue,
                    DerivedAdmissionIssueV1::UnsupportedSchemaVersion {
                        found: 2,
                        supported: DERIVED_GEOMETRY_SCHEMA_VERSION_V1
                    }
                )
            },
        );

        let tiny = DerivedAdmissionBudgetV1 {
            max_objects: 1,
            ..DerivedAdmissionBudgetV1::STANDARD
        };
        assert_issue(
            admit_derived_geometry_v1(regular_linkage(), tiny, cx),
            |issue| {
                matches!(
                    issue,
                    DerivedAdmissionIssueV1::ResourceLimit {
                        kind: DerivedObjectKindV1::Complex,
                        ..
                    }
                )
            },
        );

        let mut bad_proof = regular_linkage();
        bad_proof.proof_state = DerivedProofStateV1::ExternallyChecked {
            theorem: DerivedTheoremIdV1::from_bytes([140; 32]),
            checker: DerivedCheckerIdV1::from_bytes([141; 32]),
            receipt: witness(142),
            scope: DerivedProofScopeV1::LocalModel(local_model_id(143)),
        };
        assert_issue(
            admit_derived_geometry_v1(bad_proof, DerivedAdmissionBudgetV1::STANDARD, cx),
            |issue| matches!(issue, DerivedAdmissionIssueV1::InvalidProofState),
        );

        let structural =
            admit_derived_geometry_v1(regular_linkage(), DerivedAdmissionBudgetV1::STANDARD, cx)
                .expect("structural no-claim object must admit");
        let mut externally_checked = regular_linkage();
        externally_checked.proof_state = DerivedProofStateV1::ExternallyChecked {
            theorem: DerivedTheoremIdV1::from_bytes([144; 32]),
            checker: DerivedCheckerIdV1::from_bytes([145; 32]),
            receipt: witness(146),
            scope: DerivedProofScopeV1::Object,
        };
        let externally_checked =
            admit_derived_geometry_v1(externally_checked, DerivedAdmissionBudgetV1::STANDARD, cx)
                .expect("complete external proof metadata must admit without elevating authority");
        assert_ne!(structural.id(), externally_checked.id());
    });

    let gate = CancelGate::new();
    gate.request();
    with_gate_cx(&gate, |cx| {
        assert_issue(
            admit_derived_geometry_v1(regular_linkage(), DerivedAdmissionBudgetV1::STANDARD, cx),
            |issue| {
                matches!(
                    issue,
                    DerivedAdmissionIssueV1::Cancelled {
                        stage: "preflight",
                        completed: 0
                    }
                )
            },
        );
    });
}

#[test]
fn local_link_and_incidence_mutations_refuse() {
    with_cx(|cx| {
        let mut wrong_link = contact_corner();
        wrong_link.stratification.local_links[0].dimension = 1;
        assert_issue(
            admit_derived_geometry_v1(wrong_link, DerivedAdmissionBudgetV1::STANDARD, cx),
            |issue| {
                matches!(
                    issue,
                    DerivedAdmissionIssueV1::InvalidLocalLink { field: "dimension" }
                )
            },
        );

        let mut overflowing_link = contact_corner();
        overflowing_link.stratification.local_links[0].dimension = u32::MAX;
        assert_issue(
            admit_derived_geometry_v1(overflowing_link, DerivedAdmissionBudgetV1::STANDARD, cx),
            |issue| {
                matches!(
                    issue,
                    DerivedAdmissionIssueV1::InvalidLocalLink { field: "dimension" }
                )
            },
        );

        let mut wrong_incidence = contact_corner();
        wrong_incidence.stratification.incidences[0].codimension = 2;
        assert_issue(
            admit_derived_geometry_v1(wrong_incidence, DerivedAdmissionBudgetV1::STANDARD, cx),
            |issue| {
                matches!(
                    issue,
                    DerivedAdmissionIssueV1::InvalidStratification {
                        field: "incidence_dimension"
                    }
                )
            },
        );

        let mut missing_witness = contact_corner();
        missing_witness.stratification.incidences[0].witness = witness(0);
        let report =
            admit_derived_geometry_v1(missing_witness, DerivedAdmissionBudgetV1::STANDARD, cx)
                .expect_err("zero incidence witness must refuse");
        assert_eq!(
            report.issues(),
            &[DerivedAdmissionIssueV1::MissingIdentity {
                kind: DerivedObjectKindV1::Incidence,
                field: "witness",
            }]
        );
    });
}

#[test]
fn duplicate_reference_reports_are_permutation_stable() {
    with_cx(|cx| {
        let mut original = contact_corner();
        let other_chart = add_second_chart(&mut original, 160);
        let mut duplicate = original.constitutive_data[0].clone();
        duplicate.chart = other_chart;
        duplicate.state_dimension = 0;
        original.constitutive_data.push(duplicate);

        let mut reversed = original.clone();
        reversed.constitutive_data.reverse();

        let original_report =
            admit_derived_geometry_v1(original, DerivedAdmissionBudgetV1::STANDARD, cx)
                .expect_err("ambiguous duplicate reference must refuse");
        let reversed_report =
            admit_derived_geometry_v1(reversed, DerivedAdmissionBudgetV1::STANDARD, cx)
                .expect_err("permuted ambiguous duplicate reference must refuse");

        assert_eq!(original_report, reversed_report);
        assert_eq!(
            original_report.issues(),
            &[DerivedAdmissionIssueV1::DuplicateIdentity {
                kind: DerivedObjectKindV1::Constitutive
            }]
        );
    });
}

#[test]
fn preflight_rank_refusal_is_permutation_stable() {
    with_cx(|cx| {
        let mut original = regular_linkage();
        original.complexes[0].spaces[0].dimension = 3;
        original.complexes[0].spaces[1].dimension = 3;
        original.complexes[0].resolution.max_basis_dimension = 3;
        original.complexes[1].spaces[0].dimension = 2;
        original.complexes[1].spaces[1].dimension = 2;

        let mut reversed = original.clone();
        reversed.complexes.reverse();
        let budget = DerivedAdmissionBudgetV1 {
            max_total_rank: 4,
            ..DerivedAdmissionBudgetV1::STANDARD
        };

        let original_report = admit_derived_geometry_v1(original, budget, cx)
            .expect_err("rank excess must refuse before sorting");
        let reversed_report = admit_derived_geometry_v1(reversed, budget, cx)
            .expect_err("permuted rank excess must refuse identically");

        assert_eq!(original_report, reversed_report);
        assert_eq!(
            original_report.issues(),
            &[DerivedAdmissionIssueV1::ResourceLimit {
                kind: DerivedObjectKindV1::Complex,
                requested: 5,
                limit: 4,
            }]
        );
    });
}
