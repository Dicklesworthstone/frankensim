//! Machine-IR ordered structural assembly admission (Gauntlet G0/G3/G5).

use core::num::NonZeroU64;

use fs_blake3::ContentHash;
use fs_ir::machine::manufacturing::ManufacturingArtifactRefV1;
use fs_ir::machine::manufacturing::assembly::{
    AssemblyEndpointRoleV1, AssemblyExecutionEvidenceRefV1, AssemblyFeatureSelectorV1,
    AssemblyJoinKindV1, AssemblyOperationIdV1, AssemblyOperationModeV1, AssemblyOperationV1,
    AssemblyPathRefV1, AssemblyPreloadErrorV1, AssemblyPreloadUnitV1, AssemblyPreloadUseIssueV1,
    AssemblyPreloadV1, AssemblyProcedureRefV1, MAX_MACHINE_ASSEMBLY_OPERATIONS_V1,
    MachineAssemblyAdmissionErrorV1, MachineAssemblyDraftV1,
};
use fs_ir::machine::{
    AdmittedMachineGraph, BodyId, ContactFeatureId, MachineGraphDraft, MaterialBinding,
    MaterialCardRef, MaterialTarget, ModelRef, SubsystemId, SubsystemSpec,
};

fn nz(value: u64) -> NonZeroU64 {
    NonZeroU64::new(value).expect("fixture version is nonzero")
}

fn body(key: &str) -> BodyId {
    BodyId::new(key).expect("fixture body key is canonical")
}

fn feature(key: &str) -> ContactFeatureId {
    ContactFeatureId::new(key).expect("fixture contact-feature key is canonical")
}

fn operation_id(key: &str) -> AssemblyOperationIdV1 {
    AssemblyOperationIdV1::new(key).expect("fixture operation key is canonical")
}

fn selector(body_key: &str, feature_key: &str) -> AssemblyFeatureSelectorV1 {
    AssemblyFeatureSelectorV1::new(body(body_key), feature(feature_key))
}

fn preload(value: f64, unit: AssemblyPreloadUnitV1) -> AssemblyPreloadV1 {
    AssemblyPreloadV1::try_new(value, unit).expect("fixture preload is positive and finite")
}

fn artifact(namespace: &str, byte: u8) -> ManufacturingArtifactRefV1 {
    ManufacturingArtifactRefV1::new(namespace, nz(1), ContentHash([byte; 32]))
        .expect("fixture artifact coordinate is canonical")
}

fn procedure(byte: u8) -> AssemblyProcedureRefV1 {
    AssemblyProcedureRefV1::new(artifact("assembly/procedure", byte))
}

fn path(byte: u8) -> AssemblyPathRefV1 {
    AssemblyPathRefV1::new(artifact("assembly/path", byte))
}

fn execution_evidence(byte: u8) -> AssemblyExecutionEvidenceRefV1 {
    AssemblyExecutionEvidenceRefV1::new(artifact("assembly/execution-evidence", byte))
}

fn material(target: BodyId, key: &str, byte: u8) -> MaterialBinding {
    MaterialBinding {
        target: MaterialTarget::Body(target),
        material: MaterialCardRef::new(key, nz(1), [byte; 32])
            .expect("fixture material is canonical"),
    }
}

fn admitted_graph(model_byte: u8) -> AdmittedMachineGraph {
    let base = body("body/base");
    let base_alternate = body("body/base-alternate");
    let parts = (0..6)
        .map(|index| body(&format!("body/part-{index}")))
        .collect::<Vec<_>>();
    let other = body("body/other");

    let mut assembly_bodies = vec![base.clone(), base_alternate.clone()];
    assembly_bodies.extend(parts.iter().cloned());
    let mut assembly_features = (0..8)
        .map(|index| feature(&format!("contact/base/join-{index}")))
        .collect::<Vec<_>>();
    assembly_features.push(feature("contact/base-alternate/join-0"));
    for index in 0..6 {
        assembly_features.push(feature(&format!("contact/part-{index}/join")));
    }
    assembly_features.push(feature("contact/part-0/join-alternate"));

    let mut materials = vec![
        material(base, "materials/base", 1),
        material(base_alternate, "materials/base-alternate", 2),
    ];
    materials.extend(parts.into_iter().enumerate().map(|(index, part)| {
        material(
            part,
            &format!("materials/part-{index}"),
            u8::try_from(index + 3).expect("fixture byte fits"),
        )
    }));
    materials.push(material(other.clone(), "materials/other", 0x21));

    MachineGraphDraft {
        clocks: Vec::new(),
        subsystems: vec![
            SubsystemSpec {
                id: SubsystemId::new("subsystem/assembly").expect("canonical subsystem"),
                model: ModelRef::new("models/assembly", nz(1), [model_byte; 32])
                    .expect("canonical model"),
                bodies: assembly_bodies,
                surface_patches: Vec::new(),
                contact_features: assembly_features,
                state_slots: Vec::new(),
            },
            SubsystemSpec {
                id: SubsystemId::new("subsystem/other").expect("canonical subsystem"),
                model: ModelRef::new("models/assembly-other", nz(1), [0x52; 32])
                    .expect("canonical model"),
                bodies: vec![other],
                surface_patches: Vec::new(),
                contact_features: vec![feature("contact/other/join")],
                state_slots: Vec::new(),
            },
        ],
        terminals: Vec::new(),
        ports: Vec::new(),
        relations: Vec::new(),
        materials,
        interfaces: Vec::new(),
    }
    .admit()
    .expect("assembly fixture graph admits")
}

#[allow(clippy::too_many_arguments)]
fn operation(
    id: &str,
    ordinal: u32,
    mode: AssemblyOperationModeV1,
    base_body: &str,
    base_feature: &str,
    incoming_body: &str,
    incoming_feature: &str,
    join_kind: AssemblyJoinKindV1,
    declared_preload: Option<AssemblyPreloadV1>,
    procedure_byte: u8,
    path_byte: u8,
    evidence_byte: u8,
) -> AssemblyOperationV1 {
    AssemblyOperationV1::new(
        operation_id(id),
        ordinal,
        mode,
        selector(base_body, base_feature),
        selector(incoming_body, incoming_feature),
        join_kind,
        declared_preload,
        procedure(procedure_byte),
        path(path_byte),
        execution_evidence(evidence_byte),
    )
}

#[derive(Clone)]
struct OperationFixture {
    id: &'static str,
    ordinal: u32,
    mode: AssemblyOperationModeV1,
    base_body: &'static str,
    base_feature: &'static str,
    incoming_body: &'static str,
    incoming_feature: &'static str,
    join_kind: AssemblyJoinKindV1,
    preload: Option<AssemblyPreloadV1>,
    procedure_byte: u8,
    path_byte: u8,
    evidence_byte: u8,
}

impl OperationFixture {
    fn baseline() -> Self {
        Self {
            id: "assembly/single",
            ordinal: 0,
            mode: AssemblyOperationModeV1::AttachIncoming,
            base_body: "body/base",
            base_feature: "contact/base/join-0",
            incoming_body: "body/part-0",
            incoming_feature: "contact/part-0/join",
            join_kind: AssemblyJoinKindV1::PreloadedBolt,
            preload: Some(preload(2.0, AssemblyPreloadUnitV1::Kilonewton)),
            procedure_byte: 0x61,
            path_byte: 0x62,
            evidence_byte: 0x63,
        }
    }

    fn build(self) -> AssemblyOperationV1 {
        operation(
            self.id,
            self.ordinal,
            self.mode,
            self.base_body,
            self.base_feature,
            self.incoming_body,
            self.incoming_feature,
            self.join_kind,
            self.preload,
            self.procedure_byte,
            self.path_byte,
            self.evidence_byte,
        )
    }
}

fn singleton(operation: AssemblyOperationV1) -> MachineAssemblyDraftV1 {
    MachineAssemblyDraftV1 {
        initial_body: body("body/base"),
        operations: vec![operation],
    }
}

fn valid_draft() -> MachineAssemblyDraftV1 {
    let declarations = [
        (
            AssemblyJoinKindV1::PreloadedBolt,
            Some(preload(2.0, AssemblyPreloadUnitV1::Kilonewton)),
        ),
        (AssemblyJoinKindV1::Weld, None),
        (AssemblyJoinKindV1::AdhesiveBond, None),
        (AssemblyJoinKindV1::Key, None),
        (AssemblyJoinKindV1::Spline, None),
        (AssemblyJoinKindV1::InterferenceFit, None),
    ];
    let mut operations = declarations
        .into_iter()
        .enumerate()
        .map(|(index, (join_kind, declared_preload))| {
            let ordinal = u32::try_from(index).expect("fixture ordinal fits u32");
            operation(
                &format!("assembly/op-{index}"),
                ordinal,
                AssemblyOperationModeV1::AttachIncoming,
                "body/base",
                &format!("contact/base/join-{index}"),
                &format!("body/part-{index}"),
                &format!("contact/part-{index}/join"),
                join_kind,
                declared_preload,
                0x70 + u8::try_from(index).expect("fixture byte fits"),
                0x80 + u8::try_from(index).expect("fixture byte fits"),
                0x90 + u8::try_from(index).expect("fixture byte fits"),
            )
        })
        .collect::<Vec<_>>();
    operations.reverse();
    MachineAssemblyDraftV1 {
        initial_body: body("body/base"),
        operations,
    }
}

#[test]
#[allow(clippy::too_many_lines)]
fn mas_001_order_families_units_and_semantic_fields_reach_identity() {
    let graph = admitted_graph(0x41);
    let draft = valid_draft();
    let baseline = draft
        .clone()
        .admit_against(&graph)
        .expect("all six joining families admit");
    let mut reordered = draft;
    reordered.operations.reverse();
    let reordered = reordered
        .admit_against(&graph)
        .expect("caller collection order is non-semantic");

    assert_eq!(baseline.graph(), graph.identity());
    assert_eq!(baseline.initial_body(), &body("body/base"));
    assert_eq!(baseline.identity(), reordered.identity());
    assert_eq!(
        baseline.identity_receipt().canonical_preimage(),
        reordered.identity_receipt().canonical_preimage()
    );
    assert_eq!(
        baseline
            .operations()
            .iter()
            .map(|entry| (entry.ordinal(), entry.join_kind()))
            .collect::<Vec<_>>(),
        [
            (0, AssemblyJoinKindV1::PreloadedBolt),
            (1, AssemblyJoinKindV1::Weld),
            (2, AssemblyJoinKindV1::AdhesiveBond),
            (3, AssemblyJoinKindV1::Key),
            (4, AssemblyJoinKindV1::Spline),
            (5, AssemblyJoinKindV1::InterferenceFit),
        ]
    );
    assert!(AssemblyJoinKindV1::PreloadedBolt.requires_preload());
    assert_eq!(
        AssemblyJoinKindV1::InterferenceFit.name(),
        "interference-fit"
    );
    let retained_preload = baseline.operations()[0]
        .preload()
        .expect("preloaded bolt retains target");
    assert_eq!(retained_preload.unit().symbol(), "kN");
    assert_eq!(retained_preload.newtons().to_bits(), 2_000.0_f64.to_bits());
    assert_eq!(
        baseline.operations()[0].procedure().artifact().namespace(),
        "assembly/procedure"
    );
    assert_eq!(
        baseline.operations()[0].path().artifact().namespace(),
        "assembly/path"
    );
    assert_eq!(
        baseline.operations()[0]
            .execution_evidence()
            .artifact()
            .namespace(),
        "assembly/execution-evidence"
    );

    let thousand_newtons = preload(1_000.0, AssemblyPreloadUnitV1::Newton);
    let one_kilonewton = preload(1.0, AssemblyPreloadUnitV1::Kilonewton);
    assert_eq!(
        thousand_newtons.newtons().to_bits(),
        one_kilonewton.newtons().to_bits()
    );
    assert_ne!(thousand_newtons, one_kilonewton);

    let base_single = singleton(OperationFixture::baseline().build())
        .admit_against(&graph)
        .expect("baseline singleton admits");
    let mutations = [
        OperationFixture {
            id: "assembly/single-renamed",
            ..OperationFixture::baseline()
        },
        OperationFixture {
            base_feature: "contact/base/join-6",
            ..OperationFixture::baseline()
        },
        OperationFixture {
            incoming_body: "body/part-1",
            incoming_feature: "contact/part-1/join",
            ..OperationFixture::baseline()
        },
        OperationFixture {
            incoming_feature: "contact/part-0/join-alternate",
            ..OperationFixture::baseline()
        },
        OperationFixture {
            preload: Some(preload(3.0, AssemblyPreloadUnitV1::Kilonewton)),
            ..OperationFixture::baseline()
        },
        OperationFixture {
            preload: Some(preload(2_000.0, AssemblyPreloadUnitV1::Newton)),
            ..OperationFixture::baseline()
        },
        OperationFixture {
            procedure_byte: 0x64,
            ..OperationFixture::baseline()
        },
        OperationFixture {
            path_byte: 0x65,
            ..OperationFixture::baseline()
        },
        OperationFixture {
            evidence_byte: 0x66,
            ..OperationFixture::baseline()
        },
    ];
    for mutation in mutations {
        let changed = singleton(mutation.build())
            .admit_against(&graph)
            .expect("identity mutation remains structurally admissible");
        assert_ne!(base_single.identity(), changed.identity());
    }

    let weld = singleton(
        OperationFixture {
            join_kind: AssemblyJoinKindV1::Weld,
            preload: None,
            ..OperationFixture::baseline()
        }
        .build(),
    )
    .admit_against(&graph)
    .expect("weld singleton admits");
    let adhesive = singleton(
        OperationFixture {
            join_kind: AssemblyJoinKindV1::AdhesiveBond,
            preload: None,
            ..OperationFixture::baseline()
        }
        .build(),
    )
    .admit_against(&graph)
    .expect("adhesive singleton admits");
    assert_ne!(weld.identity(), adhesive.identity());

    let ordered = MachineAssemblyDraftV1 {
        initial_body: body("body/base"),
        operations: vec![
            OperationFixture {
                id: "assembly/order-a",
                join_kind: AssemblyJoinKindV1::Weld,
                preload: None,
                ..OperationFixture::baseline()
            }
            .build(),
            OperationFixture {
                id: "assembly/order-b",
                ordinal: 1,
                base_feature: "contact/base/join-1",
                incoming_body: "body/part-1",
                incoming_feature: "contact/part-1/join",
                join_kind: AssemblyJoinKindV1::Weld,
                preload: None,
                ..OperationFixture::baseline()
            }
            .build(),
        ],
    }
    .admit_against(&graph)
    .expect("two independent incoming bodies admit");
    let swapped = MachineAssemblyDraftV1 {
        initial_body: body("body/base"),
        operations: vec![
            OperationFixture {
                id: "assembly/order-a",
                ordinal: 1,
                join_kind: AssemblyJoinKindV1::Weld,
                preload: None,
                ..OperationFixture::baseline()
            }
            .build(),
            OperationFixture {
                id: "assembly/order-b",
                ordinal: 0,
                base_feature: "contact/base/join-1",
                incoming_body: "body/part-1",
                incoming_feature: "contact/part-1/join",
                join_kind: AssemblyJoinKindV1::Weld,
                preload: None,
                ..OperationFixture::baseline()
            }
            .build(),
        ],
    }
    .admit_against(&graph)
    .expect("swapped total order remains structurally available");
    assert_ne!(ordered.identity(), swapped.identity());

    let other_graph = admitted_graph(0x42);
    let graph_changed = singleton(OperationFixture::baseline().build())
        .admit_against(&other_graph)
        .expect("same declarations admit against changed graph");
    assert_ne!(base_single.identity(), graph_changed.identity());

    let alternate_initial = MachineAssemblyDraftV1 {
        initial_body: body("body/base-alternate"),
        operations: vec![
            OperationFixture {
                base_body: "body/base-alternate",
                base_feature: "contact/base-alternate/join-0",
                ..OperationFixture::baseline()
            }
            .build(),
        ],
    }
    .admit_against(&graph)
    .expect("alternate initial body and matching endpoint admit");
    assert_ne!(base_single.identity(), alternate_initial.identity());

    let cross_subsystem = MachineAssemblyDraftV1 {
        initial_body: body("body/base"),
        operations: vec![
            OperationFixture {
                incoming_body: "body/other",
                incoming_feature: "contact/other/join",
                join_kind: AssemblyJoinKindV1::Weld,
                preload: None,
                ..OperationFixture::baseline()
            }
            .build(),
        ],
    }
    .admit_against(&graph)
    .expect("cross-subsystem endpoints remain explicit and structurally admissible");
    assert_eq!(cross_subsystem.operations().len(), 1);

    let continued = MachineAssemblyDraftV1 {
        initial_body: body("body/base"),
        operations: vec![
            OperationFixture {
                id: "assembly/continue-attach",
                join_kind: AssemblyJoinKindV1::Weld,
                preload: None,
                ..OperationFixture::baseline()
            }
            .build(),
            OperationFixture {
                id: "assembly/continue-existing",
                ordinal: 1,
                mode: AssemblyOperationModeV1::ContinueExisting,
                base_feature: "contact/base/join-1",
                incoming_feature: "contact/part-0/join-alternate",
                join_kind: AssemblyJoinKindV1::AdhesiveBond,
                preload: None,
                ..OperationFixture::baseline()
            }
            .build(),
        ],
    }
    .admit_against(&graph)
    .expect("a distinct continuation joint between attached bodies admits");
    assert_eq!(
        continued.operations()[1].mode(),
        AssemblyOperationModeV1::ContinueExisting
    );
    let reversed_continuation = MachineAssemblyDraftV1 {
        initial_body: body("body/base"),
        operations: vec![
            OperationFixture {
                id: "assembly/continue-attach",
                join_kind: AssemblyJoinKindV1::Weld,
                preload: None,
                ..OperationFixture::baseline()
            }
            .build(),
            OperationFixture {
                id: "assembly/continue-existing",
                ordinal: 1,
                mode: AssemblyOperationModeV1::ContinueExisting,
                base_body: "body/part-0",
                base_feature: "contact/part-0/join-alternate",
                incoming_body: "body/base",
                incoming_feature: "contact/base/join-1",
                join_kind: AssemblyJoinKindV1::AdhesiveBond,
                preload: None,
                ..OperationFixture::baseline()
            }
            .build(),
        ],
    }
    .admit_against(&graph)
    .expect("reversed role order remains structurally admissible");
    assert_ne!(continued.identity(), reversed_continuation.identity());
}

#[test]
#[allow(clippy::too_many_lines)]
fn mas_002_admission_refuses_invalid_order_selectors_reuse_and_preload() {
    let graph = admitted_graph(0x41);

    assert_eq!(
        AssemblyPreloadV1::try_new(f64::NAN, AssemblyPreloadUnitV1::Newton),
        Err(AssemblyPreloadErrorV1::NonFinite)
    );
    assert_eq!(
        AssemblyPreloadV1::try_new(0.0, AssemblyPreloadUnitV1::Newton),
        Err(AssemblyPreloadErrorV1::NonPositive)
    );
    assert_eq!(
        AssemblyPreloadV1::try_new(-1.0, AssemblyPreloadUnitV1::Newton),
        Err(AssemblyPreloadErrorV1::NonPositive)
    );
    assert_eq!(
        AssemblyPreloadV1::try_new(f64::MAX, AssemblyPreloadUnitV1::Kilonewton),
        Err(AssemblyPreloadErrorV1::SiNonFinite)
    );
    assert_eq!(
        AssemblyPreloadErrorV1::SiNonFinite.code(),
        "AssemblyPreloadSiNonFinite"
    );

    assert_eq!(
        MachineAssemblyDraftV1 {
            initial_body: body("body/base"),
            operations: Vec::new(),
        }
        .admit_against(&graph),
        Err(MachineAssemblyAdmissionErrorV1::NoOperations)
    );
    assert_eq!(
        MachineAssemblyDraftV1 {
            initial_body: body("body/missing"),
            operations: vec![OperationFixture::baseline().build()],
        }
        .admit_against(&graph),
        Err(MachineAssemblyAdmissionErrorV1::UnknownInitialBody {
            body: body("body/missing"),
        })
    );

    assert_eq!(
        MachineAssemblyDraftV1 {
            initial_body: body("body/base"),
            operations: vec![
                OperationFixture::baseline().build(),
                OperationFixture {
                    ordinal: 1,
                    base_feature: "contact/base/join-1",
                    incoming_body: "body/part-1",
                    incoming_feature: "contact/part-1/join",
                    ..OperationFixture::baseline()
                }
                .build(),
            ],
        }
        .admit_against(&graph),
        Err(MachineAssemblyAdmissionErrorV1::DuplicateOperation {
            operation: operation_id("assembly/single"),
        })
    );
    assert_eq!(
        MachineAssemblyDraftV1 {
            initial_body: body("body/base"),
            operations: vec![
                OperationFixture {
                    id: "assembly/ordinal-a",
                    ..OperationFixture::baseline()
                }
                .build(),
                OperationFixture {
                    id: "assembly/ordinal-b",
                    base_feature: "contact/base/join-1",
                    incoming_body: "body/part-1",
                    incoming_feature: "contact/part-1/join",
                    ..OperationFixture::baseline()
                }
                .build(),
            ],
        }
        .admit_against(&graph),
        Err(MachineAssemblyAdmissionErrorV1::DuplicateOrdinal {
            ordinal: 0,
            first: operation_id("assembly/ordinal-a"),
            duplicate: operation_id("assembly/ordinal-b"),
        })
    );
    assert_eq!(
        singleton(
            OperationFixture {
                ordinal: 1,
                ..OperationFixture::baseline()
            }
            .build()
        )
        .admit_against(&graph),
        Err(MachineAssemblyAdmissionErrorV1::OrdinalGap {
            operation: operation_id("assembly/single"),
            expected: 0,
            actual: 1,
        })
    );

    assert_eq!(
        singleton(
            OperationFixture {
                base_body: "body/missing",
                ..OperationFixture::baseline()
            }
            .build()
        )
        .admit_against(&graph),
        Err(MachineAssemblyAdmissionErrorV1::UnknownBody {
            operation: operation_id("assembly/single"),
            role: AssemblyEndpointRoleV1::Base,
            body: body("body/missing"),
        })
    );
    assert_eq!(
        singleton(
            OperationFixture {
                incoming_feature: "contact/part-0/missing",
                ..OperationFixture::baseline()
            }
            .build()
        )
        .admit_against(&graph),
        Err(MachineAssemblyAdmissionErrorV1::UnknownFeature {
            operation: operation_id("assembly/single"),
            role: AssemblyEndpointRoleV1::Incoming,
            feature: feature("contact/part-0/missing"),
        })
    );
    assert!(matches!(
        singleton(
            OperationFixture {
                base_feature: "contact/other/join",
                ..OperationFixture::baseline()
            }
            .build()
        )
        .admit_against(&graph),
        Err(MachineAssemblyAdmissionErrorV1::FeatureOwnerMismatch {
            operation,
            role: AssemblyEndpointRoleV1::Base,
            body,
            feature,
            body_owner,
            feature_owner,
        }) if operation == operation_id("assembly/single")
            && body == crate::body("body/base")
            && feature == crate::feature("contact/other/join")
            && body_owner == SubsystemId::new("subsystem/assembly").unwrap()
            && feature_owner == SubsystemId::new("subsystem/other").unwrap()
    ));
    assert_eq!(
        singleton(
            OperationFixture {
                incoming_body: "body/base",
                incoming_feature: "contact/base/join-0",
                ..OperationFixture::baseline()
            }
            .build()
        )
        .admit_against(&graph),
        Err(MachineAssemblyAdmissionErrorV1::SameFeature {
            operation: operation_id("assembly/single"),
            feature: feature("contact/base/join-0"),
        })
    );
    assert_eq!(
        singleton(
            OperationFixture {
                incoming_body: "body/base",
                incoming_feature: "contact/base/join-1",
                ..OperationFixture::baseline()
            }
            .build()
        )
        .admit_against(&graph),
        Err(MachineAssemblyAdmissionErrorV1::SameBody {
            operation: operation_id("assembly/single"),
            body: body("body/base"),
        })
    );

    assert_eq!(
        singleton(
            OperationFixture {
                preload: None,
                ..OperationFixture::baseline()
            }
            .build()
        )
        .admit_against(&graph),
        Err(MachineAssemblyAdmissionErrorV1::InvalidPreloadUse {
            operation: operation_id("assembly/single"),
            issue: AssemblyPreloadUseIssueV1::PreloadedBoltMissing,
        })
    );
    assert_eq!(
        singleton(
            OperationFixture {
                join_kind: AssemblyJoinKindV1::Weld,
                ..OperationFixture::baseline()
            }
            .build()
        )
        .admit_against(&graph),
        Err(MachineAssemblyAdmissionErrorV1::InvalidPreloadUse {
            operation: operation_id("assembly/single"),
            issue: AssemblyPreloadUseIssueV1::NonBoltHasPreload,
        })
    );

    assert_eq!(
        singleton(
            OperationFixture {
                base_body: "body/part-1",
                base_feature: "contact/part-1/join",
                ..OperationFixture::baseline()
            }
            .build()
        )
        .admit_against(&graph),
        Err(MachineAssemblyAdmissionErrorV1::BaseUnavailable {
            operation: operation_id("assembly/single"),
            body: body("body/part-1"),
        })
    );
    assert_eq!(
        singleton(
            OperationFixture {
                mode: AssemblyOperationModeV1::ContinueExisting,
                join_kind: AssemblyJoinKindV1::Weld,
                preload: None,
                ..OperationFixture::baseline()
            }
            .build()
        )
        .admit_against(&graph),
        Err(
            MachineAssemblyAdmissionErrorV1::ContinuationIncomingUnavailable {
                operation: operation_id("assembly/single"),
                body: body("body/part-0"),
            }
        )
    );

    let first = OperationFixture {
        id: "assembly/attach-first",
        join_kind: AssemblyJoinKindV1::Weld,
        preload: None,
        ..OperationFixture::baseline()
    }
    .build();
    let second_attaches_same_body = OperationFixture {
        id: "assembly/attach-again",
        ordinal: 1,
        base_feature: "contact/base/join-1",
        incoming_feature: "contact/part-0/join-alternate",
        join_kind: AssemblyJoinKindV1::Weld,
        preload: None,
        ..OperationFixture::baseline()
    }
    .build();
    assert_eq!(
        MachineAssemblyDraftV1 {
            initial_body: body("body/base"),
            operations: vec![first, second_attaches_same_body],
        }
        .admit_against(&graph),
        Err(MachineAssemblyAdmissionErrorV1::IncomingAlreadyAttached {
            operation: operation_id("assembly/attach-again"),
            body: body("body/part-0"),
        })
    );

    let first = OperationFixture {
        id: "assembly/reuse-first",
        join_kind: AssemblyJoinKindV1::Weld,
        preload: None,
        ..OperationFixture::baseline()
    }
    .build();
    let repeated_feature = OperationFixture {
        id: "assembly/reuse-second",
        ordinal: 1,
        incoming_body: "body/part-1",
        incoming_feature: "contact/part-1/join",
        join_kind: AssemblyJoinKindV1::Weld,
        preload: None,
        ..OperationFixture::baseline()
    }
    .build();
    assert_eq!(
        MachineAssemblyDraftV1 {
            initial_body: body("body/base"),
            operations: vec![first, repeated_feature],
        }
        .admit_against(&graph),
        Err(MachineAssemblyAdmissionErrorV1::FeatureReuse {
            feature: feature("contact/base/join-0"),
            first: operation_id("assembly/reuse-first"),
            duplicate: operation_id("assembly/reuse-second"),
        })
    );
    assert_eq!(
        MachineAssemblyAdmissionErrorV1::NoOperations.code(),
        "MachineAssemblyNoOperations"
    );
}

fn boundary_graph() -> AdmittedMachineGraph {
    let base = body("body/boundary-base");
    let incoming = body("body/boundary-incoming");

    let mut features = (0..MAX_MACHINE_ASSEMBLY_OPERATIONS_V1)
        .map(|index| feature(&format!("contact/boundary-base/join-{index:04}")))
        .collect::<Vec<_>>();
    features.extend(
        (0..MAX_MACHINE_ASSEMBLY_OPERATIONS_V1)
            .map(|index| feature(&format!("contact/boundary-incoming/join-{index:04}"))),
    );

    MachineGraphDraft {
        clocks: Vec::new(),
        subsystems: vec![SubsystemSpec {
            id: SubsystemId::new("subsystem/boundary").expect("canonical subsystem"),
            model: ModelRef::new("models/assembly-boundary", nz(1), [0xA2; 32])
                .expect("canonical model"),
            bodies: vec![base.clone(), incoming.clone()],
            surface_patches: Vec::new(),
            contact_features: features,
            state_slots: Vec::new(),
        }],
        terminals: Vec::new(),
        ports: Vec::new(),
        relations: Vec::new(),
        materials: vec![
            material(base, "materials/boundary-base", 0xA1),
            material(incoming, "materials/boundary-incoming", 0xA6),
        ],
        interfaces: Vec::new(),
    }
    .admit()
    .expect("exact-cap assembly graph admits")
}

fn boundary_draft() -> MachineAssemblyDraftV1 {
    MachineAssemblyDraftV1 {
        initial_body: body("body/boundary-base"),
        operations: (0..MAX_MACHINE_ASSEMBLY_OPERATIONS_V1)
            .map(|index| {
                operation(
                    &format!("assembly/boundary-{index:04}"),
                    u32::try_from(index).expect("operation cap fits u32"),
                    if index == 0 {
                        AssemblyOperationModeV1::AttachIncoming
                    } else {
                        AssemblyOperationModeV1::ContinueExisting
                    },
                    "body/boundary-base",
                    &format!("contact/boundary-base/join-{index:04}"),
                    "body/boundary-incoming",
                    &format!("contact/boundary-incoming/join-{index:04}"),
                    AssemblyJoinKindV1::Weld,
                    None,
                    0xA3,
                    0xA4,
                    0xA5,
                )
            })
            .collect(),
    }
}

#[test]
fn mas_003_exact_resource_cap_admits_and_one_over_refuses_before_deduplication() {
    let graph = boundary_graph();
    let exact = boundary_draft();
    let admitted = exact
        .clone()
        .admit_against(&graph)
        .expect("exact assembly-operation cap admits");
    assert_eq!(
        admitted.operations().len(),
        MAX_MACHINE_ASSEMBLY_OPERATIONS_V1
    );

    let mut too_many = exact;
    let repeated = too_many.operations[0].clone();
    too_many.operations.push(repeated);
    assert_eq!(
        too_many.admit_against(&graph),
        Err(MachineAssemblyAdmissionErrorV1::OperationLimit {
            actual: MAX_MACHINE_ASSEMBLY_OPERATIONS_V1 + 1,
            max: MAX_MACHINE_ASSEMBLY_OPERATIONS_V1,
        })
    );
}

#[test]
fn mas_004_identical_input_replays_the_complete_receipt() {
    let graph = admitted_graph(0x41);
    let first = valid_draft()
        .admit_against(&graph)
        .expect("first replay admits");
    let second = valid_draft()
        .admit_against(&graph)
        .expect("second replay admits");
    assert_eq!(first.identity(), second.identity());
    assert_eq!(first.identity_receipt(), second.identity_receipt());
    assert_eq!(
        first.identity_receipt().canonical_preimage(),
        second.identity_receipt().canonical_preimage()
    );
}
