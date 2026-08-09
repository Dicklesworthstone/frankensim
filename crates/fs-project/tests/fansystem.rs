//! Fan-system declaration and lowering battery (bead frn2i.1): strict
//! validation, fail-closed lowering, analytic composition checks, identity
//! mutation, and refusal classes.

use fs_project::fansystem::{
    FAN_SYSTEM_DECL_VERSION, FanBankDecl, FanSystemDecl, FanSystemTopology, RatedPointAdmission,
    RatedPointDecl, lower_fan_system,
};
use fs_project::spec::{FanCurveDecl, FanCurvePoint, FanToleranceBasis, dims};
use fs_qty::QtyAny;

fn flow_q(value: f64) -> QtyAny {
    QtyAny::new(value, dims::VOLUMETRIC_FLOW)
}

fn pascals(value: f64) -> QtyAny {
    QtyAny::new(value, dims::PRESSURE)
}

fn curve(points: &[(f64, f64)], tolerance: f64, min_flow: f64) -> FanCurveDecl {
    FanCurveDecl {
        points: points
            .iter()
            .map(|&(q, p)| FanCurvePoint {
                flow: flow_q(q),
                static_pressure: pascals(p),
            })
            .collect(),
        pressure_tolerance_rel: tolerance,
        tolerance_basis: FanToleranceBasis::Manufacturer,
        source: "Synthetic fixture curve; not manufacturer performance data".to_string(),
        source_id: "fixture-fan-curve-v1".to_string(),
        min_flow: flow_q(min_flow),
    }
}

fn bank_decl(id: &str, curve: FanCurveDecl, count: usize, speed_ratio: f64) -> FanBankDecl {
    FanBankDecl {
        bank_id: id.to_string(),
        curve,
        count,
        arrangement: fs_airflow::FanArrangement::Series,
        speed_ratio,
        speed_ratio_domain: (0.5, 2.0),
        rated_point: None,
    }
}

fn single_system() -> FanSystemDecl {
    FanSystemDecl {
        version: FAN_SYSTEM_DECL_VERSION,
        banks: vec![bank_decl(
            "bank-a",
            curve(&[(0.0, 60.0), (0.008, 50.0), (0.02, 0.0)], 0.08, 0.002),
            2,
            1.0,
        )],
        topology: FanSystemTopology::Single,
    }
}

#[test]
fn single_bank_declaration_validates_and_lowers_lossless() {
    let decl = single_system();
    decl.validate().expect("valid declaration");
    let lowered = lower_fan_system(&decl).expect("lowering");
    assert_eq!(lowered.members.len(), 1);
    assert_eq!(lowered.members[0].0, "bank-a");
    assert_eq!(lowered.members[0].1.count(), 2);
    assert_eq!(lowered.system_bank.count(), 2);
    assert!(
        lowered
            .declaration_identity
            .starts_with("fan-system-decl:v1:")
    );
}

#[test]
fn validation_refuses_every_inferred_default_class() {
    let mut decl = single_system();
    decl.banks[0].count = 0;
    assert_eq!(
        decl.validate().expect_err("zero count").code,
        "fan-bank-count"
    );

    let mut decl = single_system();
    decl.banks[0].speed_ratio = 4.0;
    assert_eq!(
        decl.validate().expect_err("out of domain").code,
        "fan-speed-out-of-domain"
    );

    let mut decl = single_system();
    decl.banks[0].speed_ratio_domain = (0.0, 2.0);
    assert_eq!(
        decl.validate().expect_err("bad domain").code,
        "fan-speed-domain"
    );

    let mut decl = single_system();
    decl.banks[0].bank_id = String::new();
    assert_eq!(
        decl.validate().expect_err("empty identity").code,
        "fan-bank-identity"
    );

    let mut decl = single_system();
    decl.version = 99;
    assert_eq!(
        decl.validate().expect_err("stale version").code,
        "fan-system-version"
    );

    let mut decl = single_system();
    decl.banks.push(decl.banks[0].clone());
    assert_eq!(
        decl.validate().expect_err("duplicate bank").code,
        "fan-bank-duplicate"
    );
}

#[test]
fn curve_validation_covers_monotonicity_units_tolerance_and_stall() {
    let bad_flow = curve(&[(0.0, 60.0), (0.0, 50.0), (0.02, 0.0)], 0.08, 0.002);
    let mut decl = single_system();
    decl.banks[0].curve = bad_flow;
    assert_eq!(
        decl.validate().expect_err("flat flow").code,
        "fan-curve-monotonicity"
    );

    let rising = curve(&[(0.0, 10.0), (0.01, 20.0), (0.02, 30.0)], 0.08, 0.002);
    let mut decl = single_system();
    decl.banks[0].curve = rising;
    assert_eq!(
        decl.validate().expect_err("pressure rise").code,
        "fan-curve-pressure-rise"
    );

    let mut decl = single_system();
    decl.banks[0].curve.pressure_tolerance_rel = 1.5;
    assert_eq!(
        decl.validate().expect_err("tolerance").code,
        "fan-curve-tolerance"
    );

    let mut decl = single_system();
    decl.banks[0].curve.min_flow = flow_q(0.5);
    assert_eq!(
        decl.validate().expect_err("stall").code,
        "fan-stall-boundary"
    );

    let mut decl = single_system();
    decl.banks[0].curve.points[0].flow = pascals(1.0);
    assert_eq!(decl.validate().expect_err("units").code, "fan-curve-units");

    let mut decl = single_system();
    decl.banks[0].curve.source = String::new();
    assert_eq!(
        decl.validate().expect_err("provenance").code,
        "fan-curve-provenance"
    );
}

#[test]
fn rated_point_checked_against_curve_or_typed_correlation_only() {
    let mut decl = single_system();
    // Curve p(Q) = 60 - 1250 Q on [0, 0.008] then to 0 at 0.02; at
    // Q = 0.008 the curve pressure is exactly 50.
    decl.banks[0].rated_point = Some(RatedPointDecl {
        flow: flow_q(0.008),
        static_pressure: pascals(49.0),
        admission: RatedPointAdmission::CheckedWithinDeclaredTolerance,
    });
    decl.validate().expect("within 8 percent tolerance");

    let mut decl = single_system();
    decl.banks[0].rated_point = Some(RatedPointDecl {
        flow: flow_q(0.008),
        static_pressure: pascals(20.0),
        admission: RatedPointAdmission::CheckedWithinDeclaredTolerance,
    });
    assert_eq!(
        decl.validate().expect_err("beyond tolerance").code,
        "fan-rated-mismatch"
    );

    // The same disagreement is legal as correlation-only evidence.
    let mut decl = single_system();
    decl.banks[0].rated_point = Some(RatedPointDecl {
        flow: flow_q(0.008),
        static_pressure: pascals(20.0),
        admission: RatedPointAdmission::CorrelationOnly,
    });
    decl.validate().expect("correlation-only admission");

    let mut decl = single_system();
    decl.banks[0].rated_point = Some(RatedPointDecl {
        flow: flow_q(0.5),
        static_pressure: pascals(1.0),
        admission: RatedPointAdmission::CheckedWithinDeclaredTolerance,
    });
    assert_eq!(
        decl.validate().expect_err("off curve").code,
        "fan-rated-off-curve"
    );
}

#[test]
fn multi_bank_topology_must_be_explicit_and_total() {
    let mut decl = single_system();
    decl.banks.push(bank_decl(
        "bank-b",
        curve(&[(0.0, 40.0), (0.01, 0.0)], 0.05, 0.001),
        1,
        1.1,
    ));
    decl.topology = FanSystemTopology::Single;
    assert_eq!(
        decl.validate().expect_err("ambiguous single").code,
        "fan-system-topology"
    );

    decl.topology = FanSystemTopology::Series(vec!["bank-a".to_string(), "ghost".to_string()]);
    assert_eq!(
        decl.validate().expect_err("orphan member").code,
        "fan-system-orphan"
    );

    decl.topology = FanSystemTopology::Series(vec!["bank-a".to_string()]);
    assert_eq!(
        decl.validate().expect_err("short member list").code,
        "fan-system-topology"
    );

    decl.topology = FanSystemTopology::Parallel(vec![
        "bank-a".to_string(),
        "bank-b".to_string(),
        "bank-a".to_string(),
    ]);
    assert_eq!(
        decl.validate().expect_err("duplicate member").code,
        "fan-system-duplicate-member"
    );

    decl.topology = FanSystemTopology::Series(vec!["bank-b".to_string(), "bank-a".to_string()]);
    decl.validate().expect("complete explicit topology");
}

#[test]
fn composite_lowering_matches_hand_composition() {
    let mut decl = single_system();
    decl.banks[0].count = 1;
    decl.banks.push(bank_decl(
        "bank-b",
        curve(&[(0.0, 30.0), (0.02, 0.0)], 0.05, 0.001),
        1,
        1.0,
    ));
    decl.topology = FanSystemTopology::Series(vec!["bank-a".to_string(), "bank-b".to_string()]);
    let lowered = lower_fan_system(&decl).expect("lowering");
    // bank-a curve: p = 60 - 1250 Q for Q in [0, 0.008] then (60-50)/(0.008-0.02)
    // slope; bank-b: p = 30 - 1500 Q on [0, 0.02]. Shared domain lo =
    // max(0.002, 0.001) = 0.002, hi = min(0.02, 0.02) = 0.02.
    let composite = &lowered.system_bank;
    let points: Vec<(f64, f64)> = composite
        .curve()
        .points()
        .iter()
        .map(|point| (point.flow.value(), point.pressure.value()))
        .collect();
    let &(q0, p0) = points.first().expect("first");
    assert!((q0 - 0.002).abs() <= 1e-15);
    // At Q = 0.002: a gives 57.5, b gives 27; composite 84.5.
    assert!((p0 - 84.5).abs() <= 1e-9, "composite at domain low: {p0}");
    assert_eq!(composite.count(), 1);
    assert_eq!(composite.speed_ratio().to_bits(), 1.0_f64.to_bits());
}

#[test]
fn declaration_identity_is_deterministic_and_mutation_sensitive() {
    let decl = single_system();
    assert_eq!(decl.identity(), single_system().identity());
    let mut mutated = single_system();
    mutated.banks[0].speed_ratio = 1.1;
    assert_ne!(decl.identity(), mutated.identity());
    let mut mutated_topology = single_system();
    mutated_topology.topology = FanSystemTopology::Series(vec!["bank-a".to_string()]);
    assert_ne!(decl.identity(), mutated_topology.identity());
}

#[test]
fn wire_units_and_nonfinite_values_refuse() {
    let mut decl = single_system();
    decl.banks[0].rated_point = Some(RatedPointDecl {
        flow: pascals(1.0),
        static_pressure: pascals(1.0),
        admission: RatedPointAdmission::CorrelationOnly,
    });
    assert_eq!(
        decl.validate().expect_err("wrong units").code,
        "fan-rated-units"
    );

    let mut decl = single_system();
    decl.banks[0].speed_ratio = f64::NAN;
    assert_eq!(
        decl.validate().expect_err("non-finite").code,
        "fan-speed-out-of-domain"
    );
}

// ---------------------------------------------------------------------------
// Schema v2 wire surface: round trips, migration, and hostile subsections.
// ---------------------------------------------------------------------------

mod wire_surface {
    use super::*;
    use fs_project::fansystem::FanSystemTopology;
    use fs_project::{
        AirflowLeakage, Budgets, ConsequenceClass, Cooling, DecisionGate, EntityDecl, Envelope,
        Fan, GeometryArtifact, GeometryAssignment, InterfaceCardBinding, InterfaceState,
        MaterialBinding, MeshSelector, Metadata, OutputRequest, PowerDissipation, ProjectSpec,
        RequirementDirection, RequirementSeverity, RequirementSource, RequirementSourceKind,
        SafetyFactorPolicy, Seeds, SolverSettings, ThermalLimit, UnitsDoctrine, Vent, Versions,
        migrate_envelope, parse_sexpr, print_sexpr,
    };
    use fs_qty::QtyAny;
    use fs_scenario::EntityDeclaration as _;

    fn kelvin(value: f64) -> QtyAny {
        QtyAny::new(value, dims::TEMPERATURE)
    }

    fn watts(value: f64) -> QtyAny {
        QtyAny::new(value, dims::POWER)
    }

    fn reference_assembly() -> Vec<EntityDecl> {
        vec![
            EntityDecl::Assembly {
                name: "enclosure-asm".to_string(),
                display: "Enclosure".to_string(),
                expect_id: None,
            },
            EntityDecl::Part {
                parent: "enclosure-asm".to_string(),
                name: "board".to_string(),
                display: "Main board".to_string(),
                expect_id: None,
            },
            EntityDecl::Region {
                parent: "board".to_string(),
                name: "cpu".to_string(),
                display: "CPU".to_string(),
                expect_id: None,
            },
            EntityDecl::Region {
                parent: "board".to_string(),
                name: "sink-base".to_string(),
                display: "Heat sink base".to_string(),
                expect_id: None,
            },
            EntityDecl::Interface {
                parent: "enclosure-asm".to_string(),
                name: "cpu-sink-tim".to_string(),
                display: "CPU to sink TIM".to_string(),
                from: "cpu".to_string(),
                to: "sink-base".to_string(),
                expect_id: None,
            },
        ]
    }

    fn spec_with_fan_system() -> ProjectSpec {
        ProjectSpec {
            metadata: Some(Metadata {
                name: "fan-system-fixture".to_string(),
                created: "2026-08-09".to_string(),
                context_of_use: "fixture context of use".to_string(),
                intended_decision: "fixture intended decision".to_string(),
                decision_gate: DecisionGate::DesignSelection,
                consequence: ConsequenceClass::Reliability,
            }),
            versions: Some(Versions {
                schema: fs_project::FSIM_VERSION,
                constellation: "0".repeat(64),
                workspace: "e5c8061f4faed986b831b8978d0c8d1812e960fb".to_string(),
            }),
            seeds: Some(Seeds { root: 0x5EED_0001 }),
            budgets: Some(Budgets {
                solve_time: QtyAny::new(3600.0, dims::TIME),
                memory_bytes: 8 * 1024 * 1024 * 1024,
                accuracy_rel: 0.02,
            }),
            capabilities: Some(vec!["thermal.conduction-solve".to_string()]),
            units: Some(UnitsDoctrine {
                storage: "si-base".to_string(),
                display: "engineering".to_string(),
            }),
            geometry: Some(vec![GeometryArtifact {
                role: "enclosure".to_string(),
                format: "stl".to_string(),
                source_hash: 0x00ab_cdef_0123_4567,
                parser_version: "0.0.1".to_string(),
            }]),
            assignments: Some(vec![GeometryAssignment {
                artifact: "enclosure".to_string(),
                target: "cpu".to_string(),
                length_unit: "m".to_string(),
                selector: MeshSelector::NamedGroup {
                    name: "CPU".to_string(),
                },
                allow_overlap: false,
            }]),
            assembly: Some(reference_assembly()),
            materials: Some(vec![MaterialBinding {
                region: "board".to_string(),
                card: "ab".repeat(32),
                claim: None,
                state: "fr4/nominal".to_string(),
                temp_lo: kelvin(233.15),
                temp_hi: kelvin(398.15),
                source: "matdb".to_string(),
            }]),
            interface_cards: Some(vec![InterfaceCardBinding {
                interface: "cpu-sink-tim".to_string(),
                card: "cd".repeat(32),
                claim: None,
                source: "matdb".to_string(),
                state: InterfaceState::Tim {
                    thickness: QtyAny::new(100e-6, dims::LENGTH),
                    thickness_half_width: QtyAny::new(10e-6, dims::LENGTH),
                },
            }]),
            perfect_contacts: None,
            power: Some(vec![PowerDissipation {
                region: "cpu".to_string(),
                watts: watts(35.0),
                duty: 1.0,
            }]),
            cooling: Some(Cooling {
                fans: vec![Fan {
                    name: "intake-1".to_string(),
                    flow: flow_q(0.012),
                    static_pressure: pascals(45.0),
                    curve: None,
                }],
                vents: vec![Vent {
                    region: "sink-base".to_string(),
                    area: QtyAny::new(0.004, dims::AREA),
                }],
                leakage: watts(2.5),
                airflow_leakage: None,
                fan_system: Some(super::single_system()),
            }),
            envelope: Some(Envelope {
                ambient_lo: kelvin(273.15),
                ambient_hi: kelvin(318.15),
                pressure: QtyAny::new(101_325.0, dims::PRESSURE),
            }),
            requirements: Some(vec![ThermalLimit {
                qoi: "t-junction-max".to_string(),
                class: "junction".to_string(),
                region: "cpu".to_string(),
                direction: RequirementDirection::AtMost,
                limit: kelvin(378.15),
                margin: kelvin(10.0),
                source: RequirementSource {
                    kind: RequirementSourceKind::Datasheet,
                    document: "cpu-thermal-specification".to_string(),
                    version: "rev-7".to_string(),
                    locator: "table-5:tj-max".to_string(),
                },
                safety_factor: SafetyFactorPolicy {
                    factor: 1.1,
                    source: RequirementSource {
                        kind: RequirementSourceKind::InternalPolicy,
                        document: "thermal-derating-policy".to_string(),
                        version: "2026.1".to_string(),
                        locator: "section-4.2".to_string(),
                    },
                },
                severity: RequirementSeverity::ReliabilityDerating,
            }]),
            solver: Some(SolverSettings {
                fidelity: "auto".to_string(),
                tolerance_rel: 1e-6,
            }),
            outputs: Some(vec![OutputRequest {
                name: "t-junction-max".to_string(),
                kind: "scalar".to_string(),
            }]),
        }
    }

    #[test]
    fn fan_system_wire_round_trip_is_canonical() {
        let spec = spec_with_fan_system();
        let rendered = print_sexpr(&spec).expect("renders");
        assert!(rendered.contains("(fan-system"));
        assert!(rendered.contains("(topology single)"));
        let decoded = parse_sexpr(&rendered).expect("reparses");
        assert_eq!(decoded.spec, spec);
        let rerendered = print_sexpr(&decoded.spec).expect("re-renders");
        assert_eq!(rendered, rerendered);
    }

    #[test]
    fn v1_envelopes_migrate_to_v2_with_a_receipted_rewrite() {
        let mut v1_spec = spec_with_fan_system();
        v1_spec.cooling.as_mut().expect("cooling").fan_system = None;
        let v1_rendered = print_sexpr(&v1_spec).expect("v1 renders");
        let v1 = v1_rendered.replacen("(fsim-project :version 2", "(fsim-project :version 1", 1);
        assert_ne!(v1, v1_rendered, "the version rewrite must bite");

        let refusal = parse_sexpr(&v1).expect_err("v1 must not parse at v2 directly");
        assert_eq!(refusal.code, "fsim-unsupported-version");

        let migrated = migrate_envelope(&v1, 1).expect("registered v1 rule migrates");
        assert_eq!(migrated.decoded.canonical, v1_rendered);
        assert!(
            migrated
                .receipt
                .verifies(v1.as_bytes(), migrated.decoded.canonical.as_bytes())
        );
        assert_eq!(migrated.receipt.source_version, 1);
        assert_eq!(migrated.receipt.target_version, 2);
        assert_eq!(migrated.receipt.rule.label(), "cooling-fan-system-v2");
    }

    #[test]
    fn hostile_fan_system_subsections_refuse_with_structured_findings() {
        let rendered = print_sexpr(&spec_with_fan_system()).expect("renders");
        // Both hostile shapes must refuse with a typed, stable refusal —
        // either recognition findings or the canonical-bytes gate; never
        // silently admitted.
        for hostile in [
            rendered.replacen("(topology single)", "(topology sideways)", 1),
            rendered.replacen("(banks (bank", "(banks (sank", 1),
        ] {
            let refusal = parse_sexpr(&hostile).expect_err("hostile subsection refuses");
            assert!(
                matches!(
                    refusal.code,
                    "project-recognition-violations" | "fsim-non-canonical"
                ),
                "unexpected refusal code {}",
                refusal.code
            );
        }
    }
}
