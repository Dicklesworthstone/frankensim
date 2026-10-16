//! Machine-IR engineering-workflow and FMI/SSP quarantine conformance (G0/G3/G5).

use core::fmt::Write as _;
use core::num::NonZeroU64;

use fs_blake3::{ContentHash, hash_bytes, hash_domain};
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

fn bind_model_description(
    mut draft: ForeignExecutionDraftV1,
    bytes: &[u8],
) -> ForeignExecutionDraftV1 {
    draft.model_description =
        InteropArtifactRefV1::new("external/model-description", nz(1), hash_bytes(bytes))
            .expect("raw model-description digest is nonzero");
    draft
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

const FMI_MODEL_DESCRIPTION: &[u8] = b"\xef\xbb\xbf<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n\
<!-- bounded structural fixture -->\n\
<fmiModelDescription instantiationToken=\"fixture-token\" modelName=\"Fixture.Model\" fmiVersion=\"3.0\">\n\
  <CoSimulation modelIdentifier=\"Fixture_Model\"/>\n\
</fmiModelDescription>";

const SSP_MODEL_DESCRIPTION_DEFAULT_NAMESPACE: &[u8] = br#"<!-- bounded structural fixture -->
<SystemStructureDescription name="Fixture System" version="2.0" xmlns="http://ssp-standard.org/SSP1/SystemStructureDescription">
  <System name="Root"/>
</SystemStructureDescription>"#;

const SSP_MODEL_DESCRIPTION_PREFIXED_NAMESPACE: &[u8] = br#"<?xml version='1.0' encoding='utf-8'?>
<?fixture deterministic?>
<structure:SystemStructureDescription xmlns:structure="http://ssp-standard.org/SSP1/SystemStructureDescription" version="2.0" name="Fixture System">
  <structure:System name="Root"/>
</structure:SystemStructureDescription>"#;

const SSP_MODEL_DESCRIPTION_IMPLICIT_UTF8: &[u8] = br#"<?xml version="1.0"?>
<ssd:SystemStructureDescription name="Fixture System" xmlns:ssd="http://ssp-standard.org/SSP1/SystemStructureDescription" version="2.0">
  <ssd:System name="Root"/>
</ssd:SystemStructureDescription>"#;

#[test]
fn model_description_preflight_retains_exact_bytes_and_checked_type() {
    let graph = admitted_graph(0x0f);
    let workflow = workflow_draft(&graph, "model-preflight")
        .admit_against(&graph)
        .expect("workflow admits");

    for (standard, bytes) in [
        (InterchangeStandardV1::Fmi302, FMI_MODEL_DESCRIPTION),
        (
            InterchangeStandardV1::Ssp20,
            SSP_MODEL_DESCRIPTION_DEFAULT_NAMESPACE,
        ),
        (
            InterchangeStandardV1::Ssp20,
            SSP_MODEL_DESCRIPTION_PREFIXED_NAMESPACE,
        ),
        (
            InterchangeStandardV1::Ssp20,
            SSP_MODEL_DESCRIPTION_IMPLICIT_UTF8,
        ),
    ] {
        let draft = bind_model_description(
            foreign_draft(
                &graph,
                &workflow,
                ColorRank::Estimated,
                standard,
                vec![output("position", terminal_target(), "position")],
            ),
            bytes,
        );
        let nominal = draft
            .clone()
            .admit_against(&workflow, &graph)
            .expect("nominal compatibility path remains available");
        let checked = draft
            .preflight_model_description(bytes)
            .expect("official root coordinate passes bounded preflight");
        let preflight = checked.preflight();
        assert_eq!(preflight.standard(), standard);
        assert_eq!(preflight.content_hash(), hash_bytes(bytes));
        assert_eq!(preflight.byte_len(), bytes.len());
        assert_eq!(preflight.root_attribute_count(), 3);

        let admitted = checked
            .admit_against(&workflow, &graph)
            .expect("preflighted quarantine receipt admits");
        assert_eq!(admitted.preflight(), preflight);
        assert_eq!(admitted.execution().identity(), nominal.identity());
        assert_eq!(admitted.execution().standard(), standard);
        assert_eq!(
            admitted.execution().model_description().content_hash(),
            hash_bytes(bytes)
        );
        assert_eq!(admitted.execution().evidence_rank(), ColorRank::Estimated);
    }
}

#[test]
fn model_description_preflight_enforces_standard_root_coordinates() {
    let graph = admitted_graph(0x1f);
    let workflow = workflow_draft(&graph, "model-coordinate-refusals")
        .admit_against(&graph)
        .expect("workflow admits");
    let error_for = |standard, bytes: &[u8]| {
        bind_model_description(
            foreign_draft(
                &graph,
                &workflow,
                ColorRank::Estimated,
                standard,
                vec![output("position", terminal_target(), "position")],
            ),
            bytes,
        )
        .preflight_model_description(bytes)
        .expect_err("fixture must refuse")
    };

    let fmi_patch_in_xml = br#"<?xml version="1.0" encoding="UTF-8"?>
<fmiModelDescription fmiVersion="3.0.2" modelName="Fixture" instantiationToken="token">
</fmiModelDescription>"#;
    assert_eq!(
        error_for(InterchangeStandardV1::Fmi302, fmi_patch_in_xml),
        ModelDescriptionPreflightRefusalV1::WrongAttributeValue {
            name: "fmiVersion".to_owned(),
            expected: "3.0",
            actual: "3.0.2".to_owned(),
        }
    );

    let ssp_root_declared_as_fmi = br#"<?xml version="1.0" encoding="UTF-8"?>
<SystemStructureDescription xmlns="http://ssp-standard.org/SSP1/SystemStructureDescription" version="2.0" name="Fixture">
</SystemStructureDescription>"#;
    assert_eq!(
        error_for(InterchangeStandardV1::Fmi302, ssp_root_declared_as_fmi),
        ModelDescriptionPreflightRefusalV1::WrongRoot {
            expected: "fmiModelDescription",
            actual: "SystemStructureDescription".to_owned(),
        }
    );

    let fmi_without_declaration =
        br#"<fmiModelDescription fmiVersion="3.0" modelName="Fixture" instantiationToken="token">
</fmiModelDescription>"#;
    assert_eq!(
        error_for(InterchangeStandardV1::Fmi302, fmi_without_declaration),
        ModelDescriptionPreflightRefusalV1::MissingXmlDeclaration
    );

    let ssp_wrong_namespace = br#"<ssd:SystemStructureDescription xmlns:ssd="https://example.invalid/ssp" version="2.0" name="Fixture">
</ssd:SystemStructureDescription>"#;
    assert_eq!(
        error_for(InterchangeStandardV1::Ssp20, ssp_wrong_namespace),
        ModelDescriptionPreflightRefusalV1::WrongAttributeValue {
            name: "xmlns:ssd".to_owned(),
            expected: "http://ssp-standard.org/SSP1/SystemStructureDescription",
            actual: "https://example.invalid/ssp".to_owned(),
        }
    );

    let duplicate_version = br#"<SystemStructureDescription xmlns="http://ssp-standard.org/SSP1/SystemStructureDescription" version="2.0" name="Fixture" version="2.0">
</SystemStructureDescription>"#;
    assert_eq!(
        error_for(InterchangeStandardV1::Ssp20, duplicate_version),
        ModelDescriptionPreflightRefusalV1::DuplicateAttribute {
            name: "version".to_owned(),
        }
    );
}

#[test]
fn model_description_preflight_refuses_unsafe_unbounded_or_miskeyed_bytes() {
    let graph = admitted_graph(0x2f);
    let workflow = workflow_draft(&graph, "model-byte-refusals")
        .admit_against(&graph)
        .expect("workflow admits");
    let draft = || {
        foreign_draft(
            &graph,
            &workflow,
            ColorRank::Estimated,
            InterchangeStandardV1::Ssp20,
            vec![output("position", terminal_target(), "position")],
        )
    };

    assert_eq!(
        draft()
            .preflight_model_description(SSP_MODEL_DESCRIPTION_DEFAULT_NAMESPACE)
            .expect_err("nominal digest must not impersonate inspected bytes"),
        ModelDescriptionPreflightRefusalV1::ContentHashMismatch {
            declared: draft().model_description.content_hash(),
            actual: hash_bytes(SSP_MODEL_DESCRIPTION_DEFAULT_NAMESPACE),
        }
    );

    let invalid_utf8 = [0xff];
    assert_eq!(
        bind_model_description(draft(), &invalid_utf8)
            .preflight_model_description(&invalid_utf8)
            .expect_err("invalid UTF-8 must refuse"),
        ModelDescriptionPreflightRefusalV1::InvalidUtf8 { valid_up_to: 0 }
    );

    let with_doctype = br#"<!DOCTYPE SystemStructureDescription [<!ENTITY x "expanded">]>
<SystemStructureDescription xmlns="http://ssp-standard.org/SSP1/SystemStructureDescription" version="2.0" name="Fixture">
</SystemStructureDescription>"#;
    assert_eq!(
        bind_model_description(draft(), with_doctype)
            .preflight_model_description(with_doctype)
            .expect_err("DTD declarations must refuse"),
        ModelDescriptionPreflightRefusalV1::ForbiddenDeclaration {
            kind: "DOCTYPE",
            offset: 0,
        }
    );

    let nul = b"<SystemStructureDescription\0>";
    assert_eq!(
        bind_model_description(draft(), nul)
            .preflight_model_description(nul)
            .expect_err("NUL must refuse"),
        ModelDescriptionPreflightRefusalV1::NulByte { offset: 27 }
    );

    let oversized = vec![b' '; MAX_INTEROP_MODEL_DESCRIPTION_BYTES_V1 + 1];
    assert_eq!(
        draft()
            .preflight_model_description(&oversized)
            .expect_err("oversized input must refuse before hashing or parsing"),
        ModelDescriptionPreflightRefusalV1::TooLarge {
            actual: MAX_INTEROP_MODEL_DESCRIPTION_BYTES_V1 + 1,
            max: MAX_INTEROP_MODEL_DESCRIPTION_BYTES_V1,
        }
    );

    let mut oversized_prefix = vec![b' '; MAX_INTEROP_ROOT_PREFIX_BYTES_V1];
    oversized_prefix.extend_from_slice(SSP_MODEL_DESCRIPTION_DEFAULT_NAMESPACE);
    assert_eq!(
        bind_model_description(draft(), &oversized_prefix)
            .preflight_model_description(&oversized_prefix)
            .expect_err("root prefix cap must refuse deterministically"),
        ModelDescriptionPreflightRefusalV1::RootPrefixTooLarge {
            max: MAX_INTEROP_ROOT_PREFIX_BYTES_V1,
        }
    );

    let mut excessive_attributes = String::from(
        "<SystemStructureDescription xmlns=\"http://ssp-standard.org/SSP1/SystemStructureDescription\" version=\"2.0\" name=\"Fixture\"",
    );
    for index in 0..(MAX_INTEROP_ROOT_ATTRIBUTES_V1 - 2) {
        write!(excessive_attributes, " a{index}=\"x\"").expect("writing to String cannot fail");
    }
    excessive_attributes.push('>');
    assert_eq!(
        bind_model_description(draft(), excessive_attributes.as_bytes())
            .preflight_model_description(excessive_attributes.as_bytes())
            .expect_err("attribute cap must refuse deterministically"),
        ModelDescriptionPreflightRefusalV1::TooManyAttributes {
            actual: MAX_INTEROP_ROOT_ATTRIBUTES_V1 + 1,
            max: MAX_INTEROP_ROOT_ATTRIBUTES_V1,
        }
    );
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
fn foreign_receipt_binds_nominal_coordinates_and_refuses_identity_rebinding() {
    let graph = admitted_graph(0x38);
    let workflow = workflow_draft(&graph, "identity-matrix")
        .admit_against(&graph)
        .expect("workflow admits");
    let base_draft = foreign_draft(
        &graph,
        &workflow,
        ColorRank::Estimated,
        InterchangeStandardV1::Fmi302,
        vec![output("position", terminal_target(), "position")],
    );
    let base = base_draft
        .clone()
        .admit_against(&workflow, &graph)
        .expect("base receipt admits");
    let base_id = base.identity();
    let admit_id = |draft: ForeignExecutionDraftV1| {
        draft
            .admit_against(&workflow, &graph)
            .expect("single-coordinate mutation remains structurally admissible")
            .identity()
    };

    assert_eq!(base.model_description(), &base_draft.model_description);
    assert_eq!(base.adapter(), &base_draft.adapter);
    assert_eq!(base.isolation_receipt(), &base_draft.isolation_receipt);

    let mut changed = base_draft.clone();
    changed.model_description = artifact("external/model-description", "changed-model");
    assert_ne!(base_id, admit_id(changed));

    let mut changed = base_draft.clone();
    changed.adapter = InteropArtifactRefV1::new(
        "external/alternate-adapter",
        changed.adapter.schema_version(),
        changed.adapter.content_hash(),
    )
    .expect("alternate adapter coordinate is canonical");
    assert_ne!(base_id, admit_id(changed));

    let mut changed = base_draft.clone();
    changed.isolation_receipt = InteropArtifactRefV1::new(
        changed.isolation_receipt.namespace(),
        nz(2),
        changed.isolation_receipt.content_hash(),
    )
    .expect("alternate isolation schema version is valid");
    assert_ne!(base_id, admit_id(changed));

    let mut changed = base_draft.clone();
    changed.outputs[0].name = "renamed-position".to_owned();
    assert_ne!(base_id, admit_id(changed));

    let mut changed = base_draft.clone();
    changed.outputs[0].target = state_target();
    assert_ne!(base_id, admit_id(changed));

    let other_graph = admitted_graph(0x39);
    let other_workflow = workflow_draft(&other_graph, "other-graph")
        .admit_against(&other_graph)
        .expect("other workflow admits");
    let other = foreign_draft(
        &other_graph,
        &other_workflow,
        ColorRank::Estimated,
        InterchangeStandardV1::Fmi302,
        vec![output("position", terminal_target(), "position")],
    )
    .admit_against(&other_workflow, &other_graph)
    .expect("same foreign coordinates admit against another exact graph");
    assert_ne!(base_id, other.identity());

    assert_eq!(
        foreign_draft(
            &other_graph,
            &other_workflow,
            ColorRank::Estimated,
            InterchangeStandardV1::Fmi302,
            vec![output("position", terminal_target(), "position")],
        )
        .admit_against(&other_workflow, &graph),
        Err(ForeignExecutionRefusalV1::WorkflowGraphMismatch)
    );

    let mut wrong_graph = base_draft.clone();
    wrong_graph.machine_graph_id = other_graph.identity();
    assert_eq!(
        wrong_graph.admit_against(&workflow, &graph),
        Err(ForeignExecutionRefusalV1::MachineGraphMismatch)
    );

    let alternate_workflow = workflow_draft(&graph, "alternate-workflow")
        .admit_against(&graph)
        .expect("alternate workflow admits against the same graph");
    let alternate = foreign_draft(
        &graph,
        &alternate_workflow,
        ColorRank::Estimated,
        InterchangeStandardV1::Fmi302,
        vec![output("position", terminal_target(), "position")],
    )
    .admit_against(&alternate_workflow, &graph)
    .expect("same coordinates admit against the alternate exact workflow");
    assert_ne!(base_id, alternate.identity());

    let mut wrong_workflow = base_draft;
    wrong_workflow.workflow_id = alternate_workflow.identity();
    assert_eq!(
        wrong_workflow.admit_against(&workflow, &graph),
        Err(ForeignExecutionRefusalV1::WorkflowIdentityMismatch)
    );
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
