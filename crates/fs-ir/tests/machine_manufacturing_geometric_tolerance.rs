//! Machine-IR datum-backed geometric-tolerance admission (Gauntlet G0/G3/G5).

use core::num::NonZeroU64;

use fs_blake3::ContentHash;
use fs_ir::machine::manufacturing::ManufacturingArtifactRefV1;
use fs_ir::machine::manufacturing::datum_system::{
    AdmittedMachineDatumSystemV1, DatumFeatureBindingV1, DatumFeatureIdV1, DatumFeatureTargetV1,
    DatumReferenceFrameIdV1, DatumReferenceFrameV1, MachineDatumSystemDraftV1,
};
use fs_ir::machine::manufacturing::geometric_tolerance::{
    GeometricCharacteristicV1, GeometricToleranceControlIdV1, GeometricToleranceControlV1,
    GeometricToleranceDatumUseIssueV1, GeometricToleranceLengthErrorV1,
    GeometricToleranceLengthUnitV1, GeometricToleranceLengthV1,
    GeometricTolerancePresentationRefV1, GeometricToleranceSemanticSourceRefV1,
    GeometricToleranceSpecificationRefV1, MAX_MACHINE_GEOMETRIC_TOLERANCE_CONTROLS_V1,
    MachineGeometricToleranceAdmissionErrorV1, MachineGeometricToleranceDraftV1,
};
use fs_ir::machine::{
    AdmittedMachineGraph, BodyId, MachineGraphDraft, MaterialBinding, MaterialCardRef,
    MaterialTarget, ModelRef, SubsystemId, SubsystemSpec, SurfacePatchId,
};

fn nz(value: u64) -> NonZeroU64 {
    NonZeroU64::new(value).expect("fixture version is nonzero")
}

fn body(key: &str) -> BodyId {
    BodyId::new(key).expect("fixture body key is canonical")
}

fn patch(key: &str) -> SurfacePatchId {
    SurfacePatchId::new(key).expect("fixture surface-patch key is canonical")
}

fn control_id(key: &str) -> GeometricToleranceControlIdV1 {
    GeometricToleranceControlIdV1::new(key).expect("fixture control key is canonical")
}

fn datum_id(key: &str) -> DatumFeatureIdV1 {
    DatumFeatureIdV1::new(key).expect("fixture datum key is canonical")
}

fn frame_id(key: &str) -> DatumReferenceFrameIdV1 {
    DatumReferenceFrameIdV1::new(key).expect("fixture frame key is canonical")
}

fn length(value: f64, unit: GeometricToleranceLengthUnitV1) -> GeometricToleranceLengthV1 {
    GeometricToleranceLengthV1::try_new(value, unit).expect("fixture length is positive")
}

fn artifact(namespace: &str, byte: u8) -> ManufacturingArtifactRefV1 {
    ManufacturingArtifactRefV1::new(namespace, nz(1), ContentHash([byte; 32]))
        .expect("fixture artifact coordinate is canonical")
}

fn specification(byte: u8) -> GeometricToleranceSpecificationRefV1 {
    GeometricToleranceSpecificationRefV1::new(artifact("tolerance/specification", byte))
}

fn semantic_source(byte: u8) -> GeometricToleranceSemanticSourceRefV1 {
    GeometricToleranceSemanticSourceRefV1::new(artifact("tolerance/semantic-source", byte))
}

fn presentation(byte: u8) -> GeometricTolerancePresentationRefV1 {
    GeometricTolerancePresentationRefV1::new(artifact("tolerance/presentation", byte))
}

fn material(target: BodyId, key: &str, byte: u8) -> MaterialBinding {
    MaterialBinding {
        target: MaterialTarget::Body(target),
        material: MaterialCardRef::new(key, nz(1), [byte; 32])
            .expect("fixture material is canonical"),
    }
}

fn admitted_graph(model_byte: u8) -> AdmittedMachineGraph {
    let part = body("body/part");
    let alternate = body("body/part-alternate");
    let other = body("body/other");
    MachineGraphDraft {
        clocks: Vec::new(),
        subsystems: vec![
            SubsystemSpec {
                id: SubsystemId::new("subsystem/part").expect("canonical subsystem"),
                model: ModelRef::new("models/tolerance-part", nz(1), [model_byte; 32])
                    .expect("canonical model"),
                bodies: vec![part.clone(), alternate.clone()],
                surface_patches: vec![
                    patch("surface/part/datum-a"),
                    patch("surface/part/datum-b"),
                    patch("surface/part/control-a"),
                    patch("surface/part/control-b"),
                    patch("surface/part/control-c"),
                ],
                contact_features: Vec::new(),
                state_slots: Vec::new(),
            },
            SubsystemSpec {
                id: SubsystemId::new("subsystem/other").expect("canonical subsystem"),
                model: ModelRef::new("models/tolerance-other", nz(1), [0x52; 32])
                    .expect("canonical model"),
                bodies: vec![other.clone()],
                surface_patches: vec![
                    patch("surface/other/datum-a"),
                    patch("surface/other/control-a"),
                ],
                contact_features: Vec::new(),
                state_slots: Vec::new(),
            },
        ],
        terminals: Vec::new(),
        ports: Vec::new(),
        relations: Vec::new(),
        materials: vec![
            material(part, "materials/part", 1),
            material(alternate, "materials/part-alternate", 2),
            material(other, "materials/other", 3),
        ],
        interfaces: Vec::new(),
    }
    .admit()
    .expect("geometric-tolerance fixture graph admits")
}

fn surface_datum(id: &str, declared_body: &str, surface: &str) -> DatumFeatureBindingV1 {
    DatumFeatureBindingV1::new(
        datum_id(id),
        body(declared_body),
        DatumFeatureTargetV1::SurfacePatch(patch(surface)),
    )
}

fn reference_frame(id: &str, primary: &str, secondary: Option<&str>) -> DatumReferenceFrameV1 {
    DatumReferenceFrameV1::new(
        frame_id(id),
        datum_id(primary),
        secondary.map(datum_id),
        None,
    )
}

fn admitted_datum(
    graph: &AdmittedMachineGraph,
    part_primary_id: &str,
) -> AdmittedMachineDatumSystemV1 {
    MachineDatumSystemDraftV1 {
        datum_features: vec![
            surface_datum(part_primary_id, "body/part", "surface/part/datum-a"),
            surface_datum("datum/b", "body/part", "surface/part/datum-b"),
            surface_datum("datum/other", "body/other", "surface/other/datum-a"),
        ],
        reference_frames: vec![
            reference_frame("datum-frame/part", part_primary_id, Some("datum/b")),
            reference_frame(
                "datum-frame/part-alternate",
                "datum/b",
                Some(part_primary_id),
            ),
            reference_frame("datum-frame/other", "datum/other", None),
        ],
    }
    .admit_against(graph)
    .expect("fixture datum catalog admits")
}

#[allow(clippy::too_many_arguments)]
fn control(
    id: &str,
    declared_body: &str,
    controlled_patch: &str,
    characteristic: GeometricCharacteristicV1,
    zone_width: GeometricToleranceLengthV1,
    datum_frame: Option<&str>,
    specification_byte: u8,
    source_byte: u8,
    presentation_byte: Option<u8>,
) -> GeometricToleranceControlV1 {
    GeometricToleranceControlV1::new(
        control_id(id),
        body(declared_body),
        patch(controlled_patch),
        characteristic,
        zone_width,
        datum_frame.map(frame_id),
        specification(specification_byte),
        semantic_source(source_byte),
        presentation_byte.map(presentation),
    )
}

#[derive(Clone)]
struct ControlFixture {
    id: &'static str,
    declared_body: &'static str,
    controlled_patch: &'static str,
    characteristic: GeometricCharacteristicV1,
    zone_width: GeometricToleranceLengthV1,
    datum_frame: Option<&'static str>,
    specification_byte: u8,
    source_byte: u8,
    presentation_byte: Option<u8>,
}

impl ControlFixture {
    fn baseline() -> Self {
        Self {
            id: "tolerance/single",
            declared_body: "body/part",
            controlled_patch: "surface/part/control-a",
            characteristic: GeometricCharacteristicV1::Flatness,
            zone_width: length(1.0, GeometricToleranceLengthUnitV1::Millimetre),
            datum_frame: None,
            specification_byte: 0x61,
            source_byte: 0x62,
            presentation_byte: Some(0x63),
        }
    }

    fn build(self) -> GeometricToleranceControlV1 {
        control(
            self.id,
            self.declared_body,
            self.controlled_patch,
            self.characteristic,
            self.zone_width,
            self.datum_frame,
            self.specification_byte,
            self.source_byte,
            self.presentation_byte,
        )
    }
}

fn singleton(control: GeometricToleranceControlV1) -> MachineGeometricToleranceDraftV1 {
    MachineGeometricToleranceDraftV1 {
        controls: vec![control],
    }
}

fn valid_draft() -> MachineGeometricToleranceDraftV1 {
    MachineGeometricToleranceDraftV1 {
        controls: vec![
            control(
                "tolerance/perpendicularity",
                "body/part",
                "surface/part/control-c",
                GeometricCharacteristicV1::Perpendicularity,
                length(10.0, GeometricToleranceLengthUnitV1::Micrometre),
                Some("datum-frame/part"),
                0x71,
                0x72,
                Some(0x73),
            ),
            ControlFixture::baseline().build(),
            control(
                "tolerance/parallelism",
                "body/part",
                "surface/part/control-b",
                GeometricCharacteristicV1::Parallelism,
                length(0.02, GeometricToleranceLengthUnitV1::Millimetre),
                Some("datum-frame/part"),
                0x74,
                0x75,
                None,
            ),
        ],
    }
}

#[test]
#[allow(clippy::too_many_lines)]
fn mgt_001_controls_are_graph_datum_bound_order_invariant_and_identity_complete() {
    let graph = admitted_graph(0x41);
    let datum = admitted_datum(&graph, "datum/a");
    let draft = valid_draft();
    let baseline = draft
        .clone()
        .admit_against(&graph, &datum)
        .expect("flatness and orientation controls admit");
    let mut reordered = draft;
    reordered.controls.reverse();
    let reordered = reordered
        .admit_against(&graph, &datum)
        .expect("caller control order is non-semantic");

    assert_eq!(baseline.graph(), graph.identity());
    assert_eq!(baseline.datum(), datum.identity());
    assert_eq!(baseline.identity(), reordered.identity());
    assert_eq!(
        baseline.identity_receipt().canonical_preimage(),
        reordered.identity_receipt().canonical_preimage()
    );
    assert_eq!(
        baseline
            .controls()
            .iter()
            .map(|entry| entry.id().canonical_key())
            .collect::<Vec<_>>(),
        [
            "tolerance/parallelism",
            "tolerance/perpendicularity",
            "tolerance/single",
        ]
    );
    assert!(
        baseline.controls()[0]
            .characteristic()
            .requires_datum_frame()
    );
    assert_eq!(
        baseline.controls()[0]
            .datum_frame()
            .expect("parallelism retains datum")
            .canonical_key(),
        "datum-frame/part"
    );
    assert_eq!(
        baseline.controls()[2].characteristic(),
        GeometricCharacteristicV1::Flatness
    );
    assert!(baseline.controls()[2].datum_frame().is_none());
    assert_eq!(
        baseline.controls()[2]
            .specification()
            .artifact()
            .namespace(),
        "tolerance/specification"
    );
    assert_eq!(
        baseline.controls()[2]
            .semantic_source()
            .artifact()
            .namespace(),
        "tolerance/semantic-source"
    );
    assert_eq!(
        baseline.controls()[2]
            .presentation()
            .expect("presentation is retained")
            .artifact()
            .namespace(),
        "tolerance/presentation"
    );

    let one_metre = length(1.0, GeometricToleranceLengthUnitV1::Metre);
    let thousand_millimetres = length(1_000.0, GeometricToleranceLengthUnitV1::Millimetre);
    assert_eq!(
        one_metre.metres().to_bits(),
        thousand_millimetres.metres().to_bits()
    );
    assert_ne!(one_metre, thousand_millimetres);
    assert_eq!(one_metre.unit().symbol(), "m");
    assert_eq!(
        thousand_millimetres.submitted_value().to_bits(),
        1_000.0_f64.to_bits()
    );

    let base_single = singleton(ControlFixture::baseline().build())
        .admit_against(&graph, &datum)
        .expect("baseline singleton admits");
    let mutations = [
        ControlFixture {
            id: "tolerance/single-renamed",
            ..ControlFixture::baseline()
        },
        ControlFixture {
            declared_body: "body/part-alternate",
            ..ControlFixture::baseline()
        },
        ControlFixture {
            controlled_patch: "surface/part/control-b",
            ..ControlFixture::baseline()
        },
        ControlFixture {
            characteristic: GeometricCharacteristicV1::Parallelism,
            datum_frame: Some("datum-frame/part"),
            ..ControlFixture::baseline()
        },
        ControlFixture {
            zone_width: length(1_000.0, GeometricToleranceLengthUnitV1::Micrometre),
            ..ControlFixture::baseline()
        },
        ControlFixture {
            specification_byte: 0x64,
            ..ControlFixture::baseline()
        },
        ControlFixture {
            source_byte: 0x65,
            ..ControlFixture::baseline()
        },
        ControlFixture {
            presentation_byte: Some(0x66),
            ..ControlFixture::baseline()
        },
        ControlFixture {
            presentation_byte: None,
            ..ControlFixture::baseline()
        },
    ];
    for mutation in mutations {
        let changed = singleton(mutation.build())
            .admit_against(&graph, &datum)
            .expect("identity mutation remains structurally admissible");
        assert_ne!(base_single.identity(), changed.identity());
    }

    let parallel = singleton(
        ControlFixture {
            characteristic: GeometricCharacteristicV1::Parallelism,
            datum_frame: Some("datum-frame/part"),
            ..ControlFixture::baseline()
        }
        .build(),
    )
    .admit_against(&graph, &datum)
    .expect("parallelism baseline admits");
    let perpendicular = singleton(
        ControlFixture {
            characteristic: GeometricCharacteristicV1::Perpendicularity,
            datum_frame: Some("datum-frame/part"),
            ..ControlFixture::baseline()
        }
        .build(),
    )
    .admit_against(&graph, &datum)
    .expect("perpendicularity with the same frame admits");
    let alternate_frame = singleton(
        ControlFixture {
            characteristic: GeometricCharacteristicV1::Parallelism,
            datum_frame: Some("datum-frame/part-alternate"),
            ..ControlFixture::baseline()
        }
        .build(),
    )
    .admit_against(&graph, &datum)
    .expect("parallelism with the alternate same-body frame admits");
    assert_ne!(parallel.identity(), perpendicular.identity());
    assert_ne!(parallel.identity(), alternate_frame.identity());

    let other_graph = admitted_graph(0x42);
    let other_graph_datum = admitted_datum(&other_graph, "datum/a");
    let graph_changed = singleton(ControlFixture::baseline().build())
        .admit_against(&other_graph, &other_graph_datum)
        .expect("same declaration admits against changed graph");
    assert_ne!(base_single.identity(), graph_changed.identity());

    let changed_datum = admitted_datum(&graph, "datum/a-renamed");
    assert_ne!(datum.identity(), changed_datum.identity());
    let datum_changed = singleton(ControlFixture::baseline().build())
        .admit_against(&graph, &changed_datum)
        .expect("same declaration admits against changed datum identity");
    assert_ne!(base_single.identity(), datum_changed.identity());
}

#[test]
#[allow(clippy::too_many_lines)]
fn mgt_002_admission_refuses_invalid_lengths_aliases_selectors_and_datum_use() {
    let graph = admitted_graph(0x41);
    let datum = admitted_datum(&graph, "datum/a");

    assert_eq!(
        GeometricToleranceLengthV1::try_new(f64::NAN, GeometricToleranceLengthUnitV1::Metre,),
        Err(GeometricToleranceLengthErrorV1::NonFinite)
    );
    assert_eq!(
        GeometricToleranceLengthV1::try_new(0.0, GeometricToleranceLengthUnitV1::Metre),
        Err(GeometricToleranceLengthErrorV1::NonPositive)
    );
    assert_eq!(
        GeometricToleranceLengthV1::try_new(-1.0, GeometricToleranceLengthUnitV1::Metre),
        Err(GeometricToleranceLengthErrorV1::NonPositive)
    );
    assert_eq!(
        GeometricToleranceLengthV1::try_new(
            f64::from_bits(1),
            GeometricToleranceLengthUnitV1::Nanometre,
        ),
        Err(GeometricToleranceLengthErrorV1::SiUnderflow)
    );
    assert_eq!(
        GeometricToleranceLengthErrorV1::SiUnderflow.code(),
        "GeometricToleranceLengthSiUnderflow"
    );

    let other_graph = admitted_graph(0x42);
    let mismatched_datum = admitted_datum(&other_graph, "datum/a");
    assert_eq!(
        singleton(ControlFixture::baseline().build()).admit_against(&graph, &mismatched_datum),
        Err(
            MachineGeometricToleranceAdmissionErrorV1::GraphDatumMismatch {
                graph: graph.identity(),
                datum_graph: other_graph.identity(),
            }
        )
    );
    assert_eq!(
        MachineGeometricToleranceDraftV1 {
            controls: Vec::new(),
        }
        .admit_against(&graph, &datum),
        Err(MachineGeometricToleranceAdmissionErrorV1::NoControls)
    );

    let duplicate = ControlFixture::baseline().build();
    assert_eq!(
        MachineGeometricToleranceDraftV1 {
            controls: vec![
                duplicate,
                ControlFixture {
                    controlled_patch: "surface/part/control-b",
                    ..ControlFixture::baseline()
                }
                .build(),
            ],
        }
        .admit_against(&graph, &datum),
        Err(
            MachineGeometricToleranceAdmissionErrorV1::DuplicateControl {
                control: control_id("tolerance/single"),
            }
        )
    );

    assert_eq!(
        MachineGeometricToleranceDraftV1 {
            controls: vec![
                ControlFixture {
                    id: "tolerance/selector-a",
                    ..ControlFixture::baseline()
                }
                .build(),
                ControlFixture {
                    id: "tolerance/selector-b",
                    declared_body: "body/part-alternate",
                    ..ControlFixture::baseline()
                }
                .build(),
            ],
        }
        .admit_against(&graph, &datum),
        Err(
            MachineGeometricToleranceAdmissionErrorV1::DuplicateControlSelector {
                first: control_id("tolerance/selector-a"),
                duplicate: control_id("tolerance/selector-b"),
            }
        )
    );

    assert_eq!(
        singleton(
            ControlFixture {
                declared_body: "body/missing",
                ..ControlFixture::baseline()
            }
            .build()
        )
        .admit_against(&graph, &datum),
        Err(MachineGeometricToleranceAdmissionErrorV1::UnknownBody {
            control: control_id("tolerance/single"),
            body: body("body/missing"),
        })
    );
    assert_eq!(
        singleton(
            ControlFixture {
                controlled_patch: "surface/part/missing",
                ..ControlFixture::baseline()
            }
            .build()
        )
        .admit_against(&graph, &datum),
        Err(
            MachineGeometricToleranceAdmissionErrorV1::UnknownSurfacePatch {
                control: control_id("tolerance/single"),
                patch: patch("surface/part/missing"),
            }
        )
    );
    assert!(matches!(
        singleton(
            ControlFixture {
                controlled_patch: "surface/other/control-a",
                ..ControlFixture::baseline()
            }
            .build()
        )
        .admit_against(&graph, &datum),
        Err(
            MachineGeometricToleranceAdmissionErrorV1::SurfaceOwnerMismatch {
                control,
                body,
                patch,
                body_owner,
                patch_owner,
            }
        ) if control == control_id("tolerance/single")
            && body == crate::body("body/part")
            && patch == crate::patch("surface/other/control-a")
            && body_owner == SubsystemId::new("subsystem/part").unwrap()
            && patch_owner == SubsystemId::new("subsystem/other").unwrap()
    ));

    assert_eq!(
        singleton(
            ControlFixture {
                datum_frame: Some("datum-frame/part"),
                ..ControlFixture::baseline()
            }
            .build()
        )
        .admit_against(&graph, &datum),
        Err(MachineGeometricToleranceAdmissionErrorV1::InvalidDatumUse {
            control: control_id("tolerance/single"),
            issue: GeometricToleranceDatumUseIssueV1::FlatnessHasDatum,
        })
    );
    assert_eq!(
        singleton(
            ControlFixture {
                characteristic: GeometricCharacteristicV1::Parallelism,
                ..ControlFixture::baseline()
            }
            .build()
        )
        .admit_against(&graph, &datum),
        Err(MachineGeometricToleranceAdmissionErrorV1::InvalidDatumUse {
            control: control_id("tolerance/single"),
            issue: GeometricToleranceDatumUseIssueV1::OrientationMissingDatum,
        })
    );
    assert_eq!(
        singleton(
            ControlFixture {
                characteristic: GeometricCharacteristicV1::Perpendicularity,
                datum_frame: Some("datum-frame/missing"),
                ..ControlFixture::baseline()
            }
            .build()
        )
        .admit_against(&graph, &datum),
        Err(
            MachineGeometricToleranceAdmissionErrorV1::UnknownDatumFrame {
                control: control_id("tolerance/single"),
                frame: frame_id("datum-frame/missing"),
            }
        )
    );
    assert_eq!(
        singleton(
            ControlFixture {
                declared_body: "body/other",
                controlled_patch: "surface/other/control-a",
                characteristic: GeometricCharacteristicV1::Parallelism,
                datum_frame: Some("datum-frame/part"),
                ..ControlFixture::baseline()
            }
            .build()
        )
        .admit_against(&graph, &datum),
        Err(
            MachineGeometricToleranceAdmissionErrorV1::DatumBodyMismatch {
                control: control_id("tolerance/single"),
                frame: frame_id("datum-frame/part"),
                controlled_body: body("body/other"),
                datum_body: body("body/part"),
            }
        )
    );
    assert_eq!(
        MachineGeometricToleranceAdmissionErrorV1::NoControls.code(),
        "MachineGeometricToleranceNoControls"
    );
}

fn boundary_graph() -> AdmittedMachineGraph {
    let target = body("body/boundary");
    let mut patches = vec![patch("surface/boundary/datum")];
    patches.extend(
        (0..MAX_MACHINE_GEOMETRIC_TOLERANCE_CONTROLS_V1)
            .map(|index| patch(&format!("surface/boundary/control-{index:04}"))),
    );
    MachineGraphDraft {
        clocks: Vec::new(),
        subsystems: vec![SubsystemSpec {
            id: SubsystemId::new("subsystem/boundary").expect("canonical subsystem"),
            model: ModelRef::new("models/tolerance-boundary", nz(1), [0x81; 32])
                .expect("canonical model"),
            bodies: vec![target.clone()],
            surface_patches: patches,
            contact_features: Vec::new(),
            state_slots: Vec::new(),
        }],
        terminals: Vec::new(),
        ports: Vec::new(),
        relations: Vec::new(),
        materials: vec![material(target, "materials/boundary", 0x82)],
        interfaces: Vec::new(),
    }
    .admit()
    .expect("exact-cap geometric-tolerance graph admits")
}

fn boundary_datum(graph: &AdmittedMachineGraph) -> AdmittedMachineDatumSystemV1 {
    MachineDatumSystemDraftV1 {
        datum_features: vec![surface_datum(
            "datum/boundary",
            "body/boundary",
            "surface/boundary/datum",
        )],
        reference_frames: vec![reference_frame(
            "datum-frame/boundary",
            "datum/boundary",
            None,
        )],
    }
    .admit_against(graph)
    .expect("boundary datum catalog admits")
}

fn boundary_draft() -> MachineGeometricToleranceDraftV1 {
    MachineGeometricToleranceDraftV1 {
        controls: (0..MAX_MACHINE_GEOMETRIC_TOLERANCE_CONTROLS_V1)
            .map(|index| {
                control(
                    &format!("tolerance/boundary-{index:04}"),
                    "body/boundary",
                    &format!("surface/boundary/control-{index:04}"),
                    GeometricCharacteristicV1::Flatness,
                    length(10.0, GeometricToleranceLengthUnitV1::Micrometre),
                    None,
                    0x83,
                    0x84,
                    Some(0x85),
                )
            })
            .collect(),
    }
}

#[test]
fn mgt_003_exact_resource_cap_admits_and_one_over_refuses_before_deduplication() {
    let graph = boundary_graph();
    let datum = boundary_datum(&graph);
    let exact = boundary_draft();
    let admitted = exact
        .clone()
        .admit_against(&graph, &datum)
        .expect("exact geometric-tolerance control cap admits");
    assert_eq!(
        admitted.controls().len(),
        MAX_MACHINE_GEOMETRIC_TOLERANCE_CONTROLS_V1
    );

    let mut too_many = exact;
    let repeated = too_many.controls[0].clone();
    too_many.controls.push(repeated);
    assert_eq!(
        too_many.admit_against(&graph, &datum),
        Err(MachineGeometricToleranceAdmissionErrorV1::ControlLimit {
            actual: MAX_MACHINE_GEOMETRIC_TOLERANCE_CONTROLS_V1 + 1,
            max: MAX_MACHINE_GEOMETRIC_TOLERANCE_CONTROLS_V1,
        })
    );
}

#[test]
fn mgt_004_identical_input_replays_the_complete_receipt() {
    let graph = admitted_graph(0x41);
    let datum = admitted_datum(&graph, "datum/a");
    let first = valid_draft()
        .admit_against(&graph, &datum)
        .expect("first replay admits");
    let second = valid_draft()
        .admit_against(&graph, &datum)
        .expect("second replay admits");
    assert_eq!(first.identity(), second.identity());
    assert_eq!(first.identity_receipt(), second.identity_receipt());
    assert_eq!(
        first.identity_receipt().canonical_preimage(),
        second.identity_receipt().canonical_preimage()
    );
}
