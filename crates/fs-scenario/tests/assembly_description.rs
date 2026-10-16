//! IRP M6.1-W1a: G0 schema/refusals, G3 rename, G5 unchanged IR pins.
//! All comparisons use exact bytes/identities or zero structural violations.
use fs_blake3::{ContentHash, hash_bytes};
use fs_qty::{Dims, QtyAny};
use fs_scenario::assembly::*;
use fs_scenario::{
    DatumFeature, EntityCatalog, EntityDeclaration, Environment, FrameId, GeometryFingerprint,
    InterfacePair, PlacementBasis, Scenario, Tolerance, ToleranceKind, ToleranceSource,
};

const B: AssemblyBudget = DEFAULT_ASSEMBLY_BUDGET;
const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);
const PRESSURE: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const LEGACY: &str = concat!(
    "(scenario \"legacy-v1\" 7 (environment ",
    "(qty 0 1 0 -2 0 0) (qty 0 1 0 -2 0 0) (qty -9.80665 1 0 -2 0 0) ",
    "(qty 293.15 0 0 0 1 0) (qty 101325 -1 1 -2 0 0)) ",
    "(frames) (bcs) (cases) (combos) (ensembles) (contacts))"
);
fn hash(s: &str) -> ContentHash {
    hash_bytes(s.as_bytes())
}
fn measured() -> ToleranceSource {
    ToleranceSource::Measured {
        ledger_row: "test-fixture:measurement".into(),
        method: "test-fixture:calipers".into(),
    }
}
fn reconstructed() -> ToleranceSource {
    ToleranceSource::Reconstructed {
        ledger_row: "test-fixture:scan".into(),
        route: "test-fixture:calibrated-reconstruction".into(),
    }
}
fn bound() -> DimensionBound {
    DimensionBound {
        property: ParameterKey::new("coordinate", LENGTH),
        lower: QtyAny::new(-1e-4, LENGTH),
        upper: QtyAny::new(1e-4, LENGTH),
        source: reconstructed(),
    }
}
fn fixture() -> (Scenario, AssemblyInput) {
    let scenario = Scenario::new("test-fixture:all-kinds", 7, Environment::earth_lab());
    let mut catalog = EntityCatalog::new();
    let root = catalog
        .declare(EntityDeclaration::assembly("test-fixture:root"))
        .unwrap();
    let mut parts = Vec::new();
    let mut sides = Vec::new();
    for (index, kind) in [
        PartKind::String,
        PartKind::Shell,
        PartKind::Solid,
        PartKind::Duct,
        PartKind::Cavity,
        PartKind::ValvePattern,
        PartKind::DrivenBody,
    ]
    .into_iter()
    .enumerate()
    {
        let fingerprint =
            GeometryFingerprint::of_bytes(format!("test-fixture:geometry-{index}").as_bytes());
        let entity = catalog
            .declare(
                EntityDeclaration::part(root, &format!("test-fixture:part-{index}"))
                    .with_fingerprint(fingerprint),
            )
            .unwrap();
        sides.push(
            catalog
                .declare(EntityDeclaration::surface(entity, "test-fixture:surface"))
                .unwrap(),
        );
        catalog
            .declare_placement(
                entity,
                FrameId(0),
                if index == 0 {
                    PlacementBasis::AsBuilt {
                        registration_ref: hash("test-fixture:registration"),
                    }
                } else {
                    PlacementBasis::Nominal
                },
            )
            .unwrap();
        let geometry = if index % 2 == 0 {
            GeometrySource::Imported {
                fingerprint,
                representation: match index {
                    0 => Representation::Mesh,
                    2 => Representation::Nurbs,
                    4 => Representation::Sdf,
                    _ => Representation::Voxel,
                },
                ledger_row: "test-fixture:source".into(),
            }
        } else {
            GeometrySource::Parametric {
                generator: "test-fixture:generator".into(),
                params: vec![(
                    ParameterKey::new("length", LENGTH),
                    QtyAny::new(0.5, LENGTH),
                )],
            }
        };
        parts.push(PartDecl {
            entity,
            kind,
            geometry,
            materials: vec![MaterialRef {
                card: hash("test-fixture:material"),
                query_point: vec![(
                    ParameterKey::new("temperature", Dims([0, 0, 0, 1, 0, 0])),
                    QtyAny::new(293.15, Dims([0, 0, 0, 1, 0, 0])),
                )],
                frame: FrameId(0),
                layer: Some(0),
                unknown: vec![
                    UnknownProperty::Unknown {
                        property: ParameterKey::new("loss_factor", Dims::NONE),
                        bound: None,
                    },
                    UnknownProperty::Unknown {
                        property: ParameterKey::new("coordinate", LENGTH),
                        bound: Some(bound()),
                    },
                ],
            }],
            uncertainty: vec![bound()],
        });
    }
    let datum = catalog
        .declare_datum(
            parts[0].entity,
            "test-fixture:datum",
            DatumFeature::Plane,
            &[],
        )
        .unwrap();
    let secondary = catalog
        .declare_datum(
            parts[0].entity,
            "test-fixture:secondary",
            DatumFeature::Axis,
            &[datum],
        )
        .unwrap();
    catalog
        .declare_tolerance(Tolerance::new(
            parts[0].entity,
            ToleranceKind::Position,
            QtyAny::new(1e-4, LENGTH),
            vec![datum, secondary],
            measured(),
        ))
        .unwrap();
    catalog
        .declare_tolerance(Tolerance::new(
            parts[1].entity,
            ToleranceKind::Flatness,
            QtyAny::new(2e-4, LENGTH),
            vec![],
            reconstructed(),
        ))
        .unwrap();
    let mut interfaces = Vec::new();
    for (index, kind) in [
        InterfaceKind::PowerPort {
            schema: hash("test-fixture:port"),
        },
        InterfaceKind::Contact {
            law: hash("test-fixture:contact"),
        },
        InterfaceKind::Glue {
            law: hash("test-fixture:glue"),
        },
    ]
    .into_iter()
    .enumerate()
    {
        let pair = if index == 0 {
            InterfacePair::ordered(sides[index], sides[index + 1])
        } else {
            InterfacePair::unordered(sides[index], sides[index + 1])
        };
        let entity = catalog
            .declare(EntityDeclaration::interface(
                root,
                &format!("test-fixture:interface-{index}"),
                pair,
            ))
            .unwrap();
        interfaces.push(InterfaceDecl { entity, kind });
    }
    let path = ParameterPath {
        entity: parts[5].entity,
        parameter: ParameterKey::new("ToneHole.open_fraction", Dims::NONE),
    };
    let configurations = vec![
        ConfigurationState {
            name: "test-fixture:declared".into(),
            assignments: vec![(path.clone(), QtyAny::new(0.25, Dims::NONE))],
            provenance: StateProvenance::Declared {
                rationale: "synthetic schema fixture, not a specimen".into(),
            },
        },
        ConfigurationState {
            name: "test-fixture:measured".into(),
            assignments: vec![(path.clone(), QtyAny::new(0.5, Dims::NONE))],
            provenance: StateProvenance::PlayingMeasured {
                ledger_row: "test-fixture:player".into(),
            },
        },
        ConfigurationState {
            name: "test-fixture:speech".into(),
            assignments: vec![(path, QtyAny::new(1.0, Dims::NONE))],
            provenance: StateProvenance::SpeechImaging {
                ledger_row: "test-fixture:speech".into(),
                subject: "test-fixture:subject".into(),
                token: "test-fixture:vowel".into(),
                playing_envelope: None,
            },
        },
    ];
    let observers = vec![ObserverDecl {
        entity: root,
        frame: FrameId(0),
        quantity: ParameterKey::new("pressure", PRESSURE),
        position: [
            QtyAny::new(0.0, LENGTH),
            QtyAny::new(1.0, LENGTH),
            QtyAny::new(0.0, LENGTH),
        ],
    }];
    (
        scenario,
        AssemblyInput {
            root,
            catalog,
            parts,
            interfaces,
            configurations,
            observers,
            reference_state: Some(ReferenceState::Assembled),
        },
    )
}
fn event(test: &str) {
    println!(
        "{{\"bead\":\"frankensim-irp-v5-instrument-realism-pexfp.176\",\"test\":\"{test}\",\"tolerance\":0,\"violations\":0,\"verdict\":\"pass\"}}"
    );
}

#[test]
fn ad_001_round_trip_byte_stable() {
    let (s, input) = fixture();
    for state in [
        ReferenceState::Unloaded,
        ReferenceState::Strung,
        ReferenceState::Assembled,
        ReferenceState::Played,
    ] {
        let mut input = input.clone();
        input.reference_state = Some(state);
        let description = AssemblyDescription::new(&s, input, B).unwrap();
        let bytes = description.to_bytes(B).unwrap();
        let decoded = AssemblyDescription::from_bytes(&bytes, &s, B).unwrap();
        assert_eq!(description, decoded);
        assert_eq!(bytes, decoded.to_bytes(B).unwrap());
        assert_eq!(description.content_hash(B).unwrap(), hash_bytes(&bytes));
        assert_eq!(
            decoded.content_hash(B).unwrap(),
            description.content_hash(B).unwrap()
        );
        assert_eq!(decoded.input().parts.len(), 7);
        assert_eq!(decoded.input().interfaces.len(), 3);
        assert_eq!(decoded.input().catalog.tolerances().len(), 2);
    }
    event("ad_001_round_trip_byte_stable");
}
#[test]
fn ad_002_scenario_ir_pins_unmoved() {
    let old = fs_scenario::ir::parse_ir(LEGACY).unwrap();
    let before = fs_scenario::ir::write_ir(old.scenario());
    let expected =
        ContentHash::from_hex("f74494d0152c0587c611cdab5dcc2f96baca05bfe50fa13a18f4a64ada15bffc")
            .unwrap();
    assert_eq!(hash_bytes(before.as_bytes()), expected);
    let (_, input) = fixture();
    let (description, receipt) = AssemblyDescription::attach_to_ir(LEGACY, input, B).unwrap();
    assert_eq!(fs_scenario::ir::SCENARIO_IR_VERSION, 2);
    assert_eq!(description.scenario_ir_hash(), expected);
    assert_eq!(fs_scenario::ir::write_ir(old.scenario()), before);
    assert!(receipt.verifies(
        LEGACY.as_bytes(),
        before.as_bytes(),
        &description.to_bytes(B).unwrap()
    ));
    assert!(
        receipt
            .rules()
            .contains(&AssemblyAttachmentRule::AppendMoleZero)
    );
    let (_, input) = fixture();
    let (attached, current) = AssemblyDescription::attach_to_ir(&before, input, B).unwrap();
    assert_eq!(
        current.rules(),
        [
            AssemblyAttachmentRule::PreserveScenarioV2,
            AssemblyAttachmentRule::AttachAssemblyV1
        ]
    );
    assert!(current.verifies(
        before.as_bytes(),
        before.as_bytes(),
        &attached.to_bytes(B).unwrap()
    ));
    event("ad_002_scenario_ir_pins_unmoved");
}
#[test]
fn ad_003_binds_scenario_hash_and_refuses_mismatch() {
    let (s, input) = fixture();
    let a = AssemblyDescription::new(&s, input, B).unwrap();
    a.verify_scenario(&s).unwrap();
    let mut changed = s.clone();
    changed.seed += 1;
    assert!(matches!(
        a.verify_scenario(&changed),
        Err(AssemblyError::ScenarioHashMismatch { .. })
    ));
    assert!(matches!(
        AssemblyDescription::from_bytes(&a.to_bytes(B).unwrap(), &changed, B),
        Err(AssemblyError::ScenarioHashMismatch { .. })
    ));
    event("ad_003_binds_scenario_hash_and_refuses_mismatch");
}
#[test]
fn ad_004_unknown_material_property_is_explicit() {
    let (s, mut input) = fixture();
    let m = &input.parts[0].materials[0];
    assert!(matches!(
        m.unknown[0],
        UnknownProperty::Unknown { bound: None, .. }
    ));
    assert_eq!(
        m.require_known("loss_factor"),
        Err(AssemblyError::UnknownPropertyDefaulted {
            property: "loss_factor".into()
        })
    );
    m.require_known("density").unwrap();
    input.parts[0].materials[0].query_point.push((
        ParameterKey::new("loss_factor", Dims::NONE),
        QtyAny::new(0.01, Dims::NONE),
    ));
    assert_eq!(
        AssemblyDescription::new(&s, input, B).unwrap_err(),
        AssemblyError::UnknownPropertyDefaulted {
            property: "loss_factor".into()
        }
    );
    event("ad_004_unknown_material_property_is_explicit");
}
#[test]
fn ad_005_reference_state_required() {
    let (s, mut input) = fixture();
    input.reference_state = None;
    assert_eq!(
        AssemblyDescription::new(&s, input, B).unwrap_err(),
        AssemblyError::MissingReferenceState
    );
    event("ad_005_reference_state_required");
}
#[test]
fn ad_006_configuration_assignment_dims_checked() {
    let (s, mut input) = fixture();
    input.configurations[0].assignments[0].1 = QtyAny::new(1000.0, PRESSURE);
    assert_eq!(
        AssemblyDescription::new(&s, input, B).unwrap_err(),
        AssemblyError::DimensionMismatch {
            parameter: "ToneHole.open_fraction".into(),
            expected: Dims::NONE,
            found: PRESSURE
        }
    );
    event("ad_006_configuration_assignment_dims_checked");
}
#[test]
fn ad_007_budget_cap_and_cap_plus_one() {
    let (s, input) = fixture();
    let a = AssemblyDescription::new(&s, input.clone(), B).unwrap();
    let bytes = a.to_bytes(B).unwrap();
    let profiles = [
        (
            "parts",
            AssemblyBudget {
                max_parts: input.parts.len(),
                ..B
            },
        ),
        (
            "interfaces",
            AssemblyBudget {
                max_interfaces: input.interfaces.len(),
                ..B
            },
        ),
        (
            "configurations",
            AssemblyBudget {
                max_configurations: input.configurations.len(),
                ..B
            },
        ),
        (
            "assignments per state",
            AssemblyBudget {
                max_assignments_per_state: 1,
                ..B
            },
        ),
        (
            "wire bytes",
            AssemblyBudget {
                max_wire_bytes: bytes.len(),
                ..B
            },
        ),
    ];
    for (resource, at_cap) in profiles {
        let accepted = AssemblyDescription::new(&s, input.clone(), at_cap).unwrap();
        assert_eq!(accepted.to_bytes(at_cap).unwrap(), bytes);
        assert!(AssemblyDescription::from_bytes(&bytes, &s, at_cap).is_ok());
        let mut over = at_cap;
        match resource {
            "parts" => over.max_parts -= 1,
            "interfaces" => over.max_interfaces -= 1,
            "configurations" => over.max_configurations -= 1,
            "assignments per state" => over.max_assignments_per_state -= 1,
            _ => over.max_wire_bytes -= 1,
        }
        let cap = match resource {
            "parts" => over.max_parts,
            "interfaces" => over.max_interfaces,
            "configurations" => over.max_configurations,
            "assignments per state" => over.max_assignments_per_state,
            _ => over.max_wire_bytes,
        };
        assert_eq!(
            AssemblyDescription::new(&s, input.clone(), over).unwrap_err(),
            AssemblyError::BudgetExceeded { resource, cap }
        );
        assert_eq!(
            AssemblyDescription::from_bytes(&bytes, &s, over).unwrap_err(),
            AssemblyError::BudgetExceeded { resource, cap }
        );
    }
    event("ad_007_budget_cap_and_cap_plus_one");
}
#[test]
fn ad_008_rename_is_receipt_not_orphan() {
    let (s, mut input) = fixture();
    let path = input.configurations[0].assignments[0].0.clone();
    let before = input.catalog.receipt_root();
    input
        .catalog
        .rename(path.entity, "test-fixture:display-renamed")
        .unwrap();
    assert_ne!(before, input.catalog.receipt_root());
    assert!(input.catalog.verify_receipts());
    assert_eq!(path, input.configurations[0].assignments[0].0);
    assert_eq!(input.catalog.get(path.entity).unwrap().id(), path.entity);
    let a = AssemblyDescription::new(&s, input, B).unwrap();
    let back = AssemblyDescription::from_bytes(&a.to_bytes(B).unwrap(), &s, B).unwrap();
    assert_eq!(back.input().catalog, a.input().catalog);
    assert_eq!(back.input().configurations[0].assignments[0].0, path);
    event("ad_008_rename_is_receipt_not_orphan");
}
#[test]
fn ad_009_decode_refuses_trailing_and_unknown_version() {
    let (s, input) = fixture();
    let a = AssemblyDescription::new(&s, input, B).unwrap();
    let bytes = a.to_bytes(B).unwrap();
    let mut trailing = bytes.clone();
    trailing.push(0);
    assert_eq!(
        AssemblyDescription::from_bytes(&trailing, &s, B).unwrap_err(),
        AssemblyError::TrailingBytes
    );
    let mut wrong = bytes.clone();
    wrong[8..12].copy_from_slice(&2_u32.to_le_bytes());
    assert_eq!(
        AssemblyDescription::from_bytes(&wrong, &s, B).unwrap_err(),
        AssemblyError::UnknownVersion { found: 2 }
    );
    for end in [0, 7, 11, bytes.len() - 1] {
        assert!(matches!(
            AssemblyDescription::from_bytes(&bytes[..end], &s, B),
            Err(AssemblyError::Wire { .. })
        ));
    }
    event("ad_009_decode_refuses_trailing_and_unknown_version");
}
#[test]
fn ad_010_tolerance_authority_and_metadata_refusals() {
    for source in [
        ToleranceSource::Measured {
            ledger_row: String::new(),
            method: "calipers".into(),
        },
        ToleranceSource::Measured {
            ledger_row: "test-fixture:source".into(),
            method: String::new(),
        },
        ToleranceSource::Reconstructed {
            ledger_row: String::new(),
            route: "scan".into(),
        },
        ToleranceSource::Reconstructed {
            ledger_row: "test-fixture:source".into(),
            route: String::new(),
        },
    ] {
        let (s, mut input) = fixture();
        input.parts[0].uncertainty[0].source = source;
        assert!(matches!(
            AssemblyDescription::new(&s, input, B),
            Err(AssemblyError::Invalid { .. })
        ));
    }
    let (s, mut input) = fixture();
    input.parts[0].materials[0].frame = FrameId(99);
    assert!(matches!(
        AssemblyDescription::new(&s, input, B),
        Err(AssemblyError::Invalid { field: "frame", .. })
    ));
    let (s, mut input) = fixture();
    if let GeometrySource::Imported { fingerprint, .. } = &mut input.parts[0].geometry {
        *fingerprint = GeometryFingerprint::of_bytes(b"test-fixture:unrelated-chart");
    } else {
        panic!("fixture part must cite an imported chart");
    }
    assert!(matches!(
        AssemblyDescription::new(&s, input, B),
        Err(AssemblyError::Invalid {
            field: "geometry fingerprint",
            ..
        })
    ));
    let (s, mut input) = fixture();
    input.parts[0].uncertainty[0].upper.value = f64::NAN;
    assert!(matches!(
        AssemblyDescription::new(&s, input, B),
        Err(AssemblyError::Invalid {
            field: "quantity",
            ..
        })
    ));
    let (s, mut input) = fixture();
    let repeated = input.configurations[0].assignments[0].clone();
    input.configurations[0].assignments.push(repeated);
    assert!(matches!(
        AssemblyDescription::new(&s, input, B),
        Err(AssemblyError::Invalid {
            field: "assignment",
            ..
        })
    ));
    event("ad_010_tolerance_authority_and_metadata_refusals");
}

#[test]
fn ad_011_import_and_retirement_receipts_survive_round_trip() {
    use fs_scenario::entity::{
        Correspondence, ImportRevision, ImportScope, ImportedEntity, MatchBasis, RebindEvent,
    };
    let (s, mut input) = fixture();
    let f = GeometryFingerprint::of_bytes(b"test-fixture:import");
    let old = input
        .catalog
        .declare(EntityDeclaration::part(input.root, "test-fixture:import-old").with_fingerprint(f))
        .unwrap();
    let declaration =
        EntityDeclaration::part(input.root, "test-fixture:import-new").with_fingerprint(f);
    let new = declaration.identity();
    let revision = ImportRevision {
        event: RebindEvent::Import,
        scope: ImportScope::Partial,
        entities: vec![ImportedEntity {
            declaration: declaration.clone(),
            correspondence: Correspondence::Auto,
        }],
    };
    input.catalog.apply_import(&revision).unwrap();
    assert!(!input.catalog.is_active(old));
    input.catalog.apply_import(&revision).unwrap();
    let changed = EntityDeclaration::part(input.root, "test-fixture:import-new")
        .with_fingerprint(GeometryFingerprint::of_bytes(b"test-fixture:changed"));
    input
        .catalog
        .apply_import(&ImportRevision {
            event: RebindEvent::Remesh,
            scope: ImportScope::Partial,
            entities: vec![ImportedEntity {
                declaration: changed.clone(),
                correspondence: Correspondence::Auto,
            }],
        })
        .unwrap();
    let asserted = EntityDeclaration::part(input.root, "test-fixture:asserted");
    let last = asserted.identity();
    input
        .catalog
        .apply_import(&ImportRevision {
            event: RebindEvent::RevisionMigration,
            scope: ImportScope::Partial,
            entities: vec![ImportedEntity {
                declaration: asserted,
                correspondence: Correspondence::Asserted(changed.identity()),
            }],
        })
        .unwrap();
    input.catalog.retire_entity(last).unwrap();
    assert!(!input.catalog.is_active(new));
    for basis in [
        MatchBasis::GeometryFingerprint,
        MatchBasis::DeclaredPath,
        MatchBasis::Asserted,
        MatchBasis::Absent,
        MatchBasis::Unchanged,
    ] {
        assert!(input.catalog.receipts().iter().any(|r| r.basis() == basis));
    }
    let a = AssemblyDescription::new(&s, input, B).unwrap();
    let bytes = a.to_bytes(B).unwrap();
    let back = AssemblyDescription::from_bytes(&bytes, &s, B).unwrap();
    assert_eq!(back.input().catalog, a.input().catalog);
    assert_eq!(back.to_bytes(B).unwrap(), bytes);
    event("ad_011_import_and_retirement_receipts_survive_round_trip");
}
