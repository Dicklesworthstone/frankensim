//! Machine-IR E7 as-built manufacturing lineage seed (Gauntlet G0/G3/G5).

use core::num::NonZeroU64;

use std::collections::BTreeSet;

use fs_blake3::ContentHash;
use fs_ir::machine::manufacturing::{
    MAX_MANUFACTURING_PROCESS_STEPS_V1, MachineManufacturingDraftV1, ManufacturingAdmissionErrorV1,
    ManufacturingArtifactRefV1, ManufacturingProcessKindV1, ManufacturingProcessStepV1,
    ManufacturingReferenceErrorV1, ManufacturingStepIdV1,
};
use fs_ir::machine::{
    AdmittedMachineGraph, BodyId, DependentKind, LineageEvent, LineageRecord, LineageRefusal,
    LineageRelation, MachineGraphDraft, MaterialBinding, MaterialCardRef, MaterialTarget, ModelRef,
    SubsystemId, SubsystemSpec,
};

fn nz(value: u64) -> NonZeroU64 {
    NonZeroU64::new(value).expect("test schema version is nonzero")
}

fn artifact(namespace: &str, byte: u8) -> ManufacturingArtifactRefV1 {
    artifact_version(namespace, 1, byte)
}

fn artifact_version(namespace: &str, version: u64, byte: u8) -> ManufacturingArtifactRefV1 {
    ManufacturingArtifactRefV1::new(namespace, nz(version), ContentHash([byte; 32]))
        .expect("valid manufacturing artifact")
}

fn admitted_graph(byte: u8) -> AdmittedMachineGraph {
    let subsystem = SubsystemId::new("subsystem/manufacturing").unwrap();
    let plant = BodyId::new("body/plant").unwrap();
    let fixture = BodyId::new("body/fixture").unwrap();
    MachineGraphDraft {
        clocks: Vec::new(),
        subsystems: vec![SubsystemSpec {
            id: subsystem,
            model: ModelRef::new("models/manufacturing", nz(1), [byte; 32]).unwrap(),
            bodies: vec![plant.clone(), fixture.clone()],
            surface_patches: Vec::new(),
            contact_features: Vec::new(),
            state_slots: Vec::new(),
        }],
        terminals: Vec::new(),
        ports: Vec::new(),
        relations: Vec::new(),
        materials: vec![
            MaterialBinding {
                target: MaterialTarget::Body(plant),
                material: MaterialCardRef::new("materials/plant", nz(1), [1; 32]).unwrap(),
            },
            MaterialBinding {
                target: MaterialTarget::Body(fixture),
                material: MaterialCardRef::new("materials/fixture", nz(1), [2; 32]).unwrap(),
            },
        ],
        interfaces: Vec::new(),
    }
    .admit()
    .expect("minimal two-body manufacturing graph admits")
}

fn process_artifacts(byte: u8) -> [ManufacturingArtifactRefV1; 5] {
    [
        artifact("manufacturing/process-specification", byte),
        artifact("manufacturing/input-material-state", byte.wrapping_add(1)),
        artifact("manufacturing/microstructure-state", byte.wrapping_add(2)),
        artifact("manufacturing/residual-stress-state", byte.wrapping_add(3)),
        artifact("manufacturing/property-state", byte.wrapping_add(4)),
    ]
}

fn step_with(
    id: &str,
    body: &str,
    predecessor: Option<&str>,
    process: ManufacturingProcessKindV1,
    artifacts: [ManufacturingArtifactRefV1; 5],
) -> ManufacturingProcessStepV1 {
    let [
        process_specification,
        input_material_state,
        microstructure_state,
        residual_stress_state,
        property_state,
    ] = artifacts;
    ManufacturingProcessStepV1::new(
        ManufacturingStepIdV1::new(id).expect("valid step id"),
        BodyId::new(body).expect("valid body id"),
        predecessor.map(|value| ManufacturingStepIdV1::new(value).expect("valid predecessor")),
        process,
        process_specification,
        input_material_state,
        microstructure_state,
        residual_stress_state,
        property_state,
    )
}

fn step(
    id: &str,
    predecessor: Option<&str>,
    process: ManufacturingProcessKindV1,
    byte: u8,
) -> ManufacturingProcessStepV1 {
    step_with(
        id,
        "body/plant",
        predecessor,
        process,
        process_artifacts(byte),
    )
}

fn draft(process_steps: Vec<ManufacturingProcessStepV1>) -> MachineManufacturingDraftV1 {
    MachineManufacturingDraftV1 {
        correlation_model: artifact("manufacturing/process-correlation-model", 90),
        process_steps,
    }
}

fn valid_steps() -> Vec<ManufacturingProcessStepV1> {
    vec![
        step(
            "process/plant/cast",
            None,
            ManufacturingProcessKindV1::Casting,
            10,
        ),
        step(
            "process/plant/machine",
            Some("process/plant/cast"),
            ManufacturingProcessKindV1::Machining,
            20,
        ),
        step(
            "process/plant/heat-treat",
            Some("process/plant/machine"),
            ManufacturingProcessKindV1::HeatTreatment,
            30,
        ),
    ]
}

#[test]
fn mm_001_process_lineage_is_graph_bound_order_invariant_and_identity_complete() {
    let graph = admitted_graph(10);
    let baseline = draft(valid_steps())
        .admit_against(&graph)
        .expect("valid manufacturing state");

    let mut reversed = valid_steps();
    reversed.reverse();
    let replay = draft(reversed)
        .admit_against(&graph)
        .expect("caller order is not semantic");

    assert_eq!(baseline.identity(), replay.identity());
    assert_eq!(
        baseline.identity_receipt().canonical_preimage(),
        replay.identity_receipt().canonical_preimage()
    );
    assert_eq!(baseline.graph(), graph.identity());
    assert_eq!(
        baseline
            .process_steps()
            .iter()
            .map(|entry| entry.id().canonical_key())
            .collect::<Vec<_>>(),
        [
            "process/plant/cast",
            "process/plant/machine",
            "process/plant/heat-treat"
        ]
    );
    assert_eq!(
        baseline.process_steps()[2].process(),
        ManufacturingProcessKindV1::HeatTreatment
    );
    assert_eq!(
        baseline.process_steps()[2]
            .predecessor()
            .expect("non-root step")
            .canonical_key(),
        "process/plant/machine"
    );

    let single_identity = |correlation_model: ManufacturingArtifactRefV1,
                           step: ManufacturingProcessStepV1| {
        MachineManufacturingDraftV1 {
            correlation_model,
            process_steps: vec![step],
        }
        .admit_against(&graph)
        .expect("identity mutation remains structurally admissible")
        .identity()
    };
    let base_artifacts = process_artifacts(100);
    let base_correlation = artifact("manufacturing/process-correlation-model", 90);
    let single_baseline = single_identity(
        base_correlation.clone(),
        step_with(
            "process/plant/identity",
            "body/plant",
            None,
            ManufacturingProcessKindV1::Casting,
            base_artifacts.clone(),
        ),
    );

    for correlation in [
        artifact("manufacturing/process-correlation-model-alternate", 90),
        artifact_version("manufacturing/process-correlation-model", 2, 90),
        artifact("manufacturing/process-correlation-model", 91),
    ] {
        assert_ne!(
            single_baseline,
            single_identity(
                correlation,
                step_with(
                    "process/plant/identity",
                    "body/plant",
                    None,
                    ManufacturingProcessKindV1::Casting,
                    base_artifacts.clone(),
                ),
            )
        );
    }
    for changed_step in [
        step_with(
            "process/plant/identity-alternate",
            "body/plant",
            None,
            ManufacturingProcessKindV1::Casting,
            base_artifacts.clone(),
        ),
        step_with(
            "process/plant/identity",
            "body/fixture",
            None,
            ManufacturingProcessKindV1::Casting,
            base_artifacts.clone(),
        ),
        step_with(
            "process/plant/identity",
            "body/plant",
            None,
            ManufacturingProcessKindV1::Forging,
            base_artifacts.clone(),
        ),
    ] {
        assert_ne!(
            single_baseline,
            single_identity(base_correlation.clone(), changed_step)
        );
    }
    for role in 0..base_artifacts.len() {
        for replacement in [
            artifact("manufacturing/alternate-coordinate", 100 + role as u8),
            artifact_version(base_artifacts[role].namespace(), 2, 100 + role as u8),
            artifact(base_artifacts[role].namespace(), 0xa0 + role as u8),
        ] {
            let mut changed = base_artifacts.clone();
            changed[role] = replacement;
            assert_ne!(
                single_baseline,
                single_identity(
                    base_correlation.clone(),
                    step_with(
                        "process/plant/identity",
                        "body/plant",
                        None,
                        ManufacturingProcessKindV1::Casting,
                        changed,
                    ),
                )
            );
        }
    }

    let reordered_chain = draft(vec![
        step(
            "process/plant/cast",
            None,
            ManufacturingProcessKindV1::Casting,
            10,
        ),
        step(
            "process/plant/machine",
            Some("process/plant/heat-treat"),
            ManufacturingProcessKindV1::Machining,
            20,
        ),
        step(
            "process/plant/heat-treat",
            Some("process/plant/cast"),
            ManufacturingProcessKindV1::HeatTreatment,
            30,
        ),
    ])
    .admit_against(&graph)
    .expect("alternate predecessor chain remains valid");
    assert_ne!(baseline.identity(), reordered_chain.identity());

    let changed_graph = draft(valid_steps())
        .admit_against(&admitted_graph(11))
        .expect("same manufacturing draft binds a different admitted graph");
    assert_ne!(baseline.identity(), changed_graph.identity());

    println!(
        "{{\"suite\":\"fs-ir/machine-manufacturing\",\"case\":\"mm-001\",\
         \"verdict\":\"pass\",\"steps\":3,\"correlation_model\":\"{}\",\
         \"detail\":\"process ordering is reconstructed from predecessors and every \
         material-state/correlation coordinate is receipt-semantic\"}}",
        baseline.correlation_model().namespace()
    );
}

#[test]
fn mm_002_malformed_or_foreign_process_chains_refuse_without_identity() {
    let graph = admitted_graph(10);

    let unknown = ManufacturingProcessStepV1::new(
        ManufacturingStepIdV1::new("process/foreign/cast").unwrap(),
        BodyId::new("body/foreign").unwrap(),
        None,
        ManufacturingProcessKindV1::Casting,
        artifact("manufacturing/process-specification", 1),
        artifact("manufacturing/input-material-state", 2),
        artifact("manufacturing/microstructure-state", 3),
        artifact("manufacturing/residual-stress-state", 4),
        artifact("manufacturing/property-state", 5),
    );
    assert!(matches!(
        draft(vec![unknown]).admit_against(&graph),
        Err(ManufacturingAdmissionErrorV1::UnknownBody { .. })
    ));

    let duplicate = step(
        "process/plant/cast",
        None,
        ManufacturingProcessKindV1::Casting,
        50,
    );
    assert!(matches!(
        draft(vec![duplicate.clone(), duplicate]).admit_against(&graph),
        Err(ManufacturingAdmissionErrorV1::DuplicateStep { .. })
    ));

    let duplicate_a = step(
        "process/plant/duplicate-a",
        None,
        ManufacturingProcessKindV1::Casting,
        51,
    );
    let duplicate_b = step(
        "process/plant/duplicate-b",
        None,
        ManufacturingProcessKindV1::Forging,
        52,
    );
    let foreign = step_with(
        "process/foreign/unknown",
        "body/foreign",
        None,
        ManufacturingProcessKindV1::Casting,
        process_artifacts(53),
    );
    let malformed = vec![
        duplicate_b.clone(),
        foreign,
        duplicate_a.clone(),
        duplicate_b,
        duplicate_a,
    ];
    let mut reversed_malformed = malformed.clone();
    reversed_malformed.reverse();
    let first_refusal = draft(malformed).admit_against(&graph).unwrap_err();
    let reordered_refusal = draft(reversed_malformed).admit_against(&graph).unwrap_err();
    assert_eq!(first_refusal, reordered_refusal);
    assert!(matches!(
        first_refusal,
        ManufacturingAdmissionErrorV1::DuplicateStep { ref step }
            if step.canonical_key() == "process/plant/duplicate-a"
    ));

    assert!(matches!(
        draft(vec![step(
            "process/plant/machine",
            Some("process/plant/missing"),
            ManufacturingProcessKindV1::Machining,
            60,
        )])
        .admit_against(&graph),
        Err(ManufacturingAdmissionErrorV1::MissingPredecessor { .. })
    ));

    let fixture_root = step_with(
        "process/fixture/root",
        "body/fixture",
        None,
        ManufacturingProcessKindV1::Casting,
        process_artifacts(61),
    );
    let cross_body = step(
        "process/plant/cross-body",
        Some("process/fixture/root"),
        ManufacturingProcessKindV1::Machining,
        62,
    );
    assert!(matches!(
        draft(vec![fixture_root, cross_body]).admit_against(&graph),
        Err(ManufacturingAdmissionErrorV1::CrossBodyPredecessor { .. })
    ));

    assert!(matches!(
        draft(vec![
            step(
                "process/plant/root-a",
                None,
                ManufacturingProcessKindV1::Casting,
                63,
            ),
            step(
                "process/plant/root-b",
                None,
                ManufacturingProcessKindV1::Forging,
                64,
            ),
        ])
        .admit_against(&graph),
        Err(ManufacturingAdmissionErrorV1::RootCardinality { roots: 2, .. })
    ));

    let root = step(
        "process/plant/root",
        None,
        ManufacturingProcessKindV1::Casting,
        70,
    );
    let fork_a = step(
        "process/plant/fork-a",
        Some("process/plant/root"),
        ManufacturingProcessKindV1::Machining,
        71,
    );
    let fork_b = step(
        "process/plant/fork-b",
        Some("process/plant/root"),
        ManufacturingProcessKindV1::Coating,
        72,
    );
    assert!(matches!(
        draft(vec![root, fork_a, fork_b]).admit_against(&graph),
        Err(ManufacturingAdmissionErrorV1::ForkedProcessChain { .. })
    ));

    let cycle_a = step(
        "process/plant/cycle-a",
        Some("process/plant/cycle-b"),
        ManufacturingProcessKindV1::Forging,
        73,
    );
    let cycle_b = step(
        "process/plant/cycle-b",
        Some("process/plant/cycle-a"),
        ManufacturingProcessKindV1::Machining,
        74,
    );
    assert!(matches!(
        draft(vec![cycle_a, cycle_b]).admit_against(&graph),
        Err(ManufacturingAdmissionErrorV1::RootCardinality { roots: 0, .. })
    ));

    let disconnected = vec![
        step(
            "process/plant/lone-root",
            None,
            ManufacturingProcessKindV1::Casting,
            75,
        ),
        step(
            "process/plant/disconnected-a",
            Some("process/plant/disconnected-b"),
            ManufacturingProcessKindV1::Forging,
            76,
        ),
        step(
            "process/plant/disconnected-b",
            Some("process/plant/disconnected-a"),
            ManufacturingProcessKindV1::Machining,
            77,
        ),
    ];
    assert!(matches!(
        draft(disconnected).admit_against(&graph),
        Err(ManufacturingAdmissionErrorV1::DisconnectedProcessChain { .. })
    ));

    assert!(matches!(
        ManufacturingArtifactRefV1::new("manufacturing/zero", nz(1), ContentHash([0; 32])),
        Err(ManufacturingReferenceErrorV1::ZeroDigest)
    ));

    println!(
        "{{\"suite\":\"fs-ir/machine-manufacturing\",\"case\":\"mm-002\",\
         \"verdict\":\"pass\",\"detail\":\"foreign bodies, duplicate or dangling \
         steps, forks, cycles, and zero artifact identities refuse before publication\"}}"
    );
}

#[test]
fn mm_003_manufacturing_history_rebinds_only_across_unambiguous_body_lineage() {
    let graph = admitted_graph(10);
    let state = draft(valid_steps())
        .admit_against(&graph)
        .expect("valid manufacturing state");
    let source = BodyId::new("body/plant").unwrap();
    let source_element: fs_ir::machine::MachineElementId = source.clone().into();
    let expected_keys: BTreeSet<String> = state
        .process_steps()
        .iter()
        .map(|step| step.id().canonical_key().to_owned())
        .collect();

    let successor = BodyId::new("body/plant-worn").unwrap();
    let successor_element: fs_ir::machine::MachineElementId = successor.clone().into();
    let wear = LineageRecord::admit(
        LineageEvent::Wear,
        vec![
            LineageRelation::new(source.clone().into(), vec![successor.clone().into()])
                .expect("one-to-one body lineage"),
        ],
        state.lineage_dependents(),
    )
    .expect("unambiguous body lineage rebinds manufacturing history");
    let expected_rebindings: BTreeSet<_> = expected_keys
        .iter()
        .cloned()
        .map(|key| (key, source_element.clone(), successor_element.clone()))
        .collect();
    let actual_rebindings: BTreeSet<_> = wear
        .rebindings()
        .iter()
        .map(|binding| {
            assert_eq!(
                binding.dependent().kind(),
                DependentKind::ManufacturingState
            );
            (
                binding.dependent().canonical_key().to_owned(),
                binding.dependent().source().clone(),
                binding.target().clone(),
            )
        })
        .collect();
    assert_eq!(actual_rebindings, expected_rebindings);

    let mut reversed_dependents = state.lineage_dependents();
    reversed_dependents.reverse();
    let replay = LineageRecord::admit(
        LineageEvent::Wear,
        vec![
            LineageRelation::new(source.clone().into(), vec![successor.clone().into()])
                .expect("one-to-one body lineage"),
        ],
        reversed_dependents,
    )
    .expect("dependent caller order is not semantic");
    assert_eq!(wear.identity(), replay.identity());

    let split = LineageRecord::admit(
        LineageEvent::Split,
        vec![
            LineageRelation::new(
                source.into(),
                vec![
                    BodyId::new("body/plant-left").unwrap().into(),
                    BodyId::new("body/plant-right").unwrap().into(),
                ],
            )
            .expect("one-to-many body lineage"),
        ],
        state.lineage_dependents(),
    )
    .expect_err("ambiguous body split must invalidate manufacturing history");
    let LineageRefusal::Ambiguous(invalidation) = split else {
        panic!("expected complete ambiguity invalidation receipt");
    };
    let exact_invalidated: BTreeSet<_> = invalidation
        .invalidated_dependents()
        .iter()
        .map(|binding| {
            assert_eq!(binding.kind(), DependentKind::ManufacturingState);
            (binding.canonical_key().to_owned(), binding.source().clone())
        })
        .collect();
    let expected_invalidated: BTreeSet<_> = expected_keys
        .iter()
        .cloned()
        .map(|key| (key, source_element.clone()))
        .collect();
    assert_eq!(exact_invalidated, expected_invalidated);
    assert_eq!(
        invalidation.considered_dependents(),
        invalidation.invalidated_dependents()
    );

    println!(
        "{{\"suite\":\"fs-ir/machine-manufacturing\",\"case\":\"mm-003\",\
         \"verdict\":\"pass\",\"rebound\":{},\"invalidated\":{},\
         \"detail\":\"durable-body process history follows one-to-one wear and is \
         completely invalidated on an ambiguous split\"}}",
        wear.rebindings().len(),
        invalidation.invalidated_dependents().len()
    );
}

#[test]
fn mm_004_maximum_legal_rows_fit_the_published_identity_envelope() {
    let graph = admitted_graph(10);
    let long_namespace = format!("manufacturing/{}", "a".repeat(114));
    assert_eq!(long_namespace.len(), 128);
    let coordinate = artifact_version(&long_namespace, u64::MAX, 0x5a);
    let mut process_steps = Vec::with_capacity(MAX_MANUFACTURING_PROCESS_STEPS_V1);
    let mut predecessor = None;
    for index in 0..MAX_MANUFACTURING_PROCESS_STEPS_V1 {
        let key = format!("process/plant/{}-{index:04}", "a".repeat(109));
        assert_eq!(key.len(), 128);
        let id = ManufacturingStepIdV1::new(key).expect("maximum-length step ID");
        process_steps.push(ManufacturingProcessStepV1::new(
            id.clone(),
            BodyId::new("body/plant").unwrap(),
            predecessor,
            ManufacturingProcessKindV1::Machining,
            coordinate.clone(),
            coordinate.clone(),
            coordinate.clone(),
            coordinate.clone(),
            coordinate.clone(),
        ));
        predecessor = Some(id);
    }

    let admitted = MachineManufacturingDraftV1 {
        correlation_model: coordinate.clone(),
        process_steps,
    }
    .admit_against(&graph)
    .expect("the exact published step and key maxima fit the identity envelope");
    assert_eq!(
        admitted.process_steps().len(),
        MAX_MANUFACTURING_PROCESS_STEPS_V1
    );

    let mut over_limit = admitted.process_steps().to_vec();
    over_limit.push(over_limit.last().expect("nonempty maximum").clone());
    assert!(matches!(
        MachineManufacturingDraftV1 {
            correlation_model: coordinate,
            process_steps: over_limit,
        }
        .admit_against(&graph),
        Err(ManufacturingAdmissionErrorV1::ProcessStepLimit {
            actual,
            max: MAX_MANUFACTURING_PROCESS_STEPS_V1,
        }) if actual == MAX_MANUFACTURING_PROCESS_STEPS_V1 + 1
    ));

    println!(
        "{{\"suite\":\"fs-ir/machine-manufacturing\",\"case\":\"mm-004\",\
         \"verdict\":\"pass\",\"steps\":{},\"key_bytes\":128,\
         \"detail\":\"the exact public maximum encodes and limit plus one refuses \
         before identity work\"}}",
        MAX_MANUFACTURING_PROCESS_STEPS_V1
    );
}
