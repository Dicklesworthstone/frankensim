//! Machine-IR engineering-workflow and FMI/SSP quarantine conformance (G0/G3/G5).

use core::num::NonZeroU64;

use fs_blake3::{ContentHash, hash_domain};
use fs_evidence::ColorRank;
use fs_ir::machine::interop::*;
use fs_ir::machine::{
    AdmittedMachineGraph, BodyId, ClockId, ClockSpec, FrameBinding, MachineClock,
    MachineGraphDraft, MaterialBinding, MaterialCardRef, MaterialTarget, ModelRef,
    OrientationParity, RelationId, RelationMode, RelationSpec, StateSlotId, SubsystemId,
    SubsystemSpec, TerminalCausality, TerminalId, TerminalQuantitySpec, TerminalShape,
    TerminalSpec,
};
use fs_qty::Dims;

fn nz(value: u64) -> NonZeroU64 {
    NonZeroU64::new(value).expect("test value is nonzero")
}

fn digest(label: &str) -> ContentHash {
    hash_domain(
        "org.frankensim.fs-ir.machine-interop-test.v1",
        label.as_bytes(),
    )
}

fn artifact(namespace: &str, label: &str) -> InteropArtifactRefV1 {
    InteropArtifactRefV1::new(namespace, nz(1), digest(label)).expect("valid artifact ref")
}

fn terminal(
    key: &str,
    owner: &SubsystemId,
    causality: TerminalCausality,
    clock: &ClockId,
) -> TerminalSpec {
    TerminalSpec {
        id: TerminalId::new(key).expect("valid terminal id"),
        owner: owner.clone(),
        quantity: TerminalQuantitySpec::Dimensional(Dims([0; 6])),
        shape: TerminalShape::Scalar,
        causality,
        clock: clock.clone(),
        frame: FrameBinding::new("world/mechanical", OrientationParity::Preserving)
            .expect("valid frame"),
    }
}

fn admitted_graph(model_byte: u8) -> AdmittedMachineGraph {
    let clock = ClockId::new("clock/continuous").expect("valid clock id");
    let subsystem = SubsystemId::new("subsystem/plant").expect("valid subsystem id");
    let state = StateSlotId::new("state/position").expect("valid state id");
    let body = BodyId::new("body/plant").expect("valid body id");
    let output = terminal(
        "terminal/position-output",
        &subsystem,
        TerminalCausality::Output,
        &clock,
    );
    let input = terminal(
        "terminal/position-input",
        &subsystem,
        TerminalCausality::Input,
        &clock,
    );

    MachineGraphDraft {
        clocks: vec![ClockSpec {
            id: clock,
            clock: MachineClock::Continuous,
        }],
        subsystems: vec![SubsystemSpec {
            id: subsystem,
            model: ModelRef::new("models/interop-fixture", nz(1), [model_byte; 32])
                .expect("valid model ref"),
            bodies: vec![body.clone()],
            surface_patches: Vec::new(),
            contact_features: Vec::new(),
            state_slots: vec![state.clone()],
        }],
        terminals: vec![output.clone(), input.clone()],
        ports: Vec::new(),
        relations: vec![RelationSpec {
            id: RelationId::new("relation/position-state").expect("valid relation id"),
            source: output.id,
            target: input.id,
            mode: RelationMode::Stateful { state_slot: state },
        }],
        materials: vec![MaterialBinding {
            target: MaterialTarget::Body(body),
            material: MaterialCardRef::new("materials/interop-fixture", nz(1), [0x44; 32])
                .expect("valid material ref"),
        }],
        interfaces: Vec::new(),
    }
    .admit()
    .expect("valid interop fixture graph")
}

fn workflow_draft(graph: &AdmittedMachineGraph, suffix: &str) -> MachineWorkflowDraftV1 {
    let steps = MachineWorkflowStageV1::ORDERED
        .into_iter()
        .map(|stage| MachineWorkflowStepV1 {
            stage,
            artifact: artifact(
                &format!("workflow/{}", stage.name()),
                &format!("{}-{suffix}", stage.name()),
            ),
        })
        .collect();
    MachineWorkflowDraftV1 {
        machine_graph_id: graph.identity(),
        steps,
    }
}

fn output(name: &str, target: ForeignOutputTargetV1, label: &str) -> ForeignOutputBindingV1 {
    ForeignOutputBindingV1 {
        name: name.to_owned(),
        target,
        artifact: artifact("external/output", label),
    }
}

fn foreign_draft(
    graph: &AdmittedMachineGraph,
    workflow: &AdmittedMachineWorkflowV1,
    claimed_rank: ColorRank,
    standard: InterchangeStandardV1,
    outputs: Vec<ForeignOutputBindingV1>,
) -> ForeignExecutionDraftV1 {
    ForeignExecutionDraftV1 {
        machine_graph_id: graph.identity(),
        workflow_id: workflow.identity(),
        standard,
        model_description: artifact("external/model-description", "model-description"),
        adapter: artifact("external/isolated-adapter", "adapter"),
        isolation_receipt: artifact("external/isolation-receipt", "isolation-receipt"),
        claimed_rank,
        outputs,
    }
}

fn terminal_target() -> ForeignOutputTargetV1 {
    ForeignOutputTargetV1::Terminal(
        TerminalId::new("terminal/position-output").expect("valid target terminal"),
    )
}

fn state_target() -> ForeignOutputTargetV1 {
    ForeignOutputTargetV1::StateSlot(
        StateSlotId::new("state/position").expect("valid target state"),
    )
}

#[test]
fn workflow_admission_requires_every_stage_in_exact_order_and_graph() {
    let graph = admitted_graph(0x11);
    let workflow = workflow_draft(&graph, "base")
        .admit_against(&graph)
        .expect("complete workflow admits");
    let replay = workflow_draft(&graph, "base")
        .admit_against(&graph)
        .expect("same workflow replays");

    assert_eq!(workflow.machine_graph_id(), graph.identity());
    assert_eq!(workflow.steps().len(), MACHINE_WORKFLOW_STAGE_COUNT_V1);
    assert_eq!(workflow.identity(), replay.identity());
    assert_eq!(workflow.identity_receipt().id(), workflow.identity());
    assert_eq!(
        workflow
            .steps()
            .iter()
            .map(|step| step.stage)
            .collect::<Vec<_>>(),
        MachineWorkflowStageV1::ORDERED.to_vec()
    );
    assert_ne!(
        workflow.identity(),
        workflow_draft(&graph, "changed-artifact")
            .admit_against(&graph)
            .expect("changed workflow still structurally admits")
            .identity()
    );

    let mut missing = workflow_draft(&graph, "missing");
    missing.steps.pop();
    assert_eq!(
        missing.admit_against(&graph),
        Err(MachineWorkflowRefusalV1::StageCount {
            actual: MACHINE_WORKFLOW_STAGE_COUNT_V1 - 1,
            expected: MACHINE_WORKFLOW_STAGE_COUNT_V1,
        })
    );

    let mut reordered = workflow_draft(&graph, "reordered");
    reordered.steps.swap(3, 4);
    let refusal = reordered
        .admit_against(&graph)
        .expect_err("reordered workflow must refuse");
    assert!(matches!(
        &refusal,
        MachineWorkflowRefusalV1::StageOrder {
            index: 3,
            expected: MachineWorkflowStageV1::DeclareScenario,
            actual: MachineWorkflowStageV1::ChooseDecisionContext,
        }
    ));
    assert_eq!(refusal.code(), "MachineWorkflowStageOrder");
    assert!(refusal.fix().contains("canonical order"));

    let other_graph = admitted_graph(0x12);
    assert!(matches!(
        workflow_draft(&graph, "wrong-graph").admit_against(&other_graph),
        Err(MachineWorkflowRefusalV1::MachineGraphMismatch { .. })
    ));
}

#[test]
fn foreign_outputs_cannot_launder_validated_or_verified_authority() {
    let graph = admitted_graph(0x21);
    let workflow = workflow_draft(&graph, "authority")
        .admit_against(&graph)
        .expect("workflow admits");

    for claimed in [ColorRank::Validated, ColorRank::Verified] {
        assert_eq!(
            foreign_draft(
                &graph,
                &workflow,
                claimed,
                InterchangeStandardV1::Fmi302,
                vec![output("position", terminal_target(), "position")],
            )
            .admit_against(&workflow, &graph),
            Err(ForeignExecutionRefusalV1::EvidenceLaundering {
                claimed,
                maximum: ColorRank::Estimated,
            })
        );
    }

    let admitted = foreign_draft(
        &graph,
        &workflow,
        ColorRank::Estimated,
        InterchangeStandardV1::Fmi302,
        vec![output("position", terminal_target(), "position")],
    )
    .admit_against(&workflow, &graph)
    .expect("Estimated foreign output admits");

    assert_eq!(admitted.evidence_rank(), ColorRank::Estimated);
    assert!(admitted.estimated_dispersion().is_infinite());
    assert_eq!(
        admitted.no_authority_policy(),
        FOREIGN_EXECUTION_NO_AUTHORITY_POLICY_V1
    );
    assert!(
        admitted
            .estimator_identity()
            .starts_with("external-adapter:")
    );
    assert_eq!(admitted.standard().name(), "fmi");
    assert_eq!(admitted.standard().version(), "3.0.2");
    assert_eq!(admitted.machine_graph_id(), graph.identity());
    assert_eq!(admitted.workflow_id(), workflow.identity());
    assert_eq!(admitted.identity_receipt().id(), admitted.identity());
}

#[test]
fn foreign_receipt_canonicalizes_output_order_and_binds_output_and_standard() {
    let graph = admitted_graph(0x31);
    let workflow = workflow_draft(&graph, "canonical")
        .admit_against(&graph)
        .expect("workflow admits");
    let first = vec![
        output("state", state_target(), "state-output"),
        output("terminal", terminal_target(), "terminal-output"),
    ];
    let second = vec![first[1].clone(), first[0].clone()];
    let left = foreign_draft(
        &graph,
        &workflow,
        ColorRank::Estimated,
        InterchangeStandardV1::Ssp20,
        first,
    )
    .admit_against(&workflow, &graph)
    .expect("first order admits");
    let right = foreign_draft(
        &graph,
        &workflow,
        ColorRank::Estimated,
        InterchangeStandardV1::Ssp20,
        second,
    )
    .admit_against(&workflow, &graph)
    .expect("second order admits");

    assert_eq!(left.identity(), right.identity());
    assert_eq!(
        left.outputs()
            .iter()
            .map(|binding| binding.name.as_str())
            .collect::<Vec<_>>(),
        vec!["state", "terminal"]
    );
    assert_eq!(left.standard().version(), "2.0");

    let changed = foreign_draft(
        &graph,
        &workflow,
        ColorRank::Estimated,
        InterchangeStandardV1::Ssp20,
        vec![
            output("state", state_target(), "changed-state-output"),
            output("terminal", terminal_target(), "terminal-output"),
        ],
    )
    .admit_against(&workflow, &graph)
    .expect("changed output still admits");
    assert_ne!(left.identity(), changed.identity());

    let fmi = foreign_draft(
        &graph,
        &workflow,
        ColorRank::Estimated,
        InterchangeStandardV1::Fmi302,
        vec![
            output("state", state_target(), "state-output"),
            output("terminal", terminal_target(), "terminal-output"),
        ],
    )
    .admit_against(&workflow, &graph)
    .expect("FMI coordinate admits");
    assert_ne!(left.identity(), fmi.identity());
}

#[test]
fn foreign_output_admission_refuses_ambiguous_unbounded_and_foreign_targets() {
    let graph = admitted_graph(0x41);
    let workflow = workflow_draft(&graph, "refusals")
        .admit_against(&graph)
        .expect("workflow admits");

    assert_eq!(
        foreign_draft(
            &graph,
            &workflow,
            ColorRank::Estimated,
            InterchangeStandardV1::Fmi302,
            Vec::new(),
        )
        .admit_against(&workflow, &graph),
        Err(ForeignExecutionRefusalV1::NoOutputs)
    );

    let template = output("same", terminal_target(), "same");
    assert_eq!(
        foreign_draft(
            &graph,
            &workflow,
            ColorRank::Estimated,
            InterchangeStandardV1::Fmi302,
            vec![template; MAX_FOREIGN_OUTPUT_BINDINGS_V1 + 1],
        )
        .admit_against(&workflow, &graph),
        Err(ForeignExecutionRefusalV1::TooManyOutputs {
            actual: MAX_FOREIGN_OUTPUT_BINDINGS_V1 + 1,
            max: MAX_FOREIGN_OUTPUT_BINDINGS_V1,
        })
    );

    assert!(matches!(
        foreign_draft(
            &graph,
            &workflow,
            ColorRank::Estimated,
            InterchangeStandardV1::Fmi302,
            vec![output(&"a".repeat(129), terminal_target(), "long-name")],
        )
        .admit_against(&workflow, &graph),
        Err(ForeignExecutionRefusalV1::InvalidOutputName { index: 0, .. })
    ));

    assert!(matches!(
        foreign_draft(
            &graph,
            &workflow,
            ColorRank::Estimated,
            InterchangeStandardV1::Fmi302,
            vec![
                output("duplicate", terminal_target(), "first"),
                output("duplicate", state_target(), "second"),
            ],
        )
        .admit_against(&workflow, &graph),
        Err(ForeignExecutionRefusalV1::DuplicateOutputName {
            first_index: 0,
            duplicate_index: 1,
            ..
        })
    ));

    assert!(matches!(
        foreign_draft(
            &graph,
            &workflow,
            ColorRank::Estimated,
            InterchangeStandardV1::Fmi302,
            vec![
                output("first", terminal_target(), "first"),
                output("second", terminal_target(), "second"),
            ],
        )
        .admit_against(&workflow, &graph),
        Err(ForeignExecutionRefusalV1::DuplicateOutputTarget {
            first_index: 0,
            duplicate_index: 1,
            ..
        })
    ));

    let unknown = ForeignOutputTargetV1::Terminal(
        TerminalId::new("terminal/not-in-graph").expect("valid foreign terminal id"),
    );
    assert!(matches!(
        foreign_draft(
            &graph,
            &workflow,
            ColorRank::Estimated,
            InterchangeStandardV1::Fmi302,
            vec![output("unknown", unknown, "unknown")],
        )
        .admit_against(&workflow, &graph),
        Err(ForeignExecutionRefusalV1::UnknownOutputTarget { index: 0, .. })
    ));
}

#[test]
fn interop_references_are_bounded_versioned_and_nonzero() {
    let reference = artifact("external/adapter", "adapter-reference");
    assert_eq!(reference.namespace(), "external/adapter");
    assert_eq!(reference.schema_version(), nz(1));
    assert_eq!(reference.content_hash(), digest("adapter-reference"));

    assert!(matches!(
        InteropArtifactRefV1::new("Bad/Namespace", nz(1), digest("bad")),
        Err(InteropReferenceErrorV1::Namespace(_))
    ));
    assert_eq!(
        InteropArtifactRefV1::new("external/adapter", nz(1), ContentHash([0; 32])),
        Err(InteropReferenceErrorV1::ZeroDigest)
    );
}
