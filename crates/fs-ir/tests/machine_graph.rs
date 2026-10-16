//! Machine-IR E0 PR-2 graph admission conformance (Gauntlet G0/G3).

use core::num::NonZeroU64;

use fs_blake3::identity::StrongIdentity;
use fs_ir::machine::{
    BodyId, ClockId, ClockSpec, ContactFeatureId, FrameBinding, InterfaceBinding, InterfaceCardRef,
    InterfaceId, InterfaceOrientation, MAX_MACHINE_GRAPH_CLOCKS, MachineClock, MachineGraphDraft,
    MachineGraphRefusal, MachineGraphRule, MaterialBinding, MaterialCardRef, MaterialTarget,
    ModelRef, OrientationParity, PortEnergyRole, PortId, PortSpec, RelationId, RelationMode,
    RelationSpec, SolvePolicyRef, StateSlotId, SubsystemId, SubsystemSpec, SurfacePatchId,
    TerminalCausality, TerminalId, TerminalQuantitySpec, TerminalShape, TerminalSpec,
};
use fs_qty::Dims;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

const PRESSURE_DIMS: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const VOLUME_FLOW_DIMS: Dims = Dims([3, 0, -1, 0, 0, 0]);

fn nz(value: u64) -> NonZeroU64 {
    NonZeroU64::new(value).expect("test value is nonzero")
}

fn digest(byte: u8) -> [u8; 32] {
    [byte; 32]
}

fn model(namespace: &str, version: u64, byte: u8) -> ModelRef {
    ModelRef::new(namespace, nz(version), digest(byte)).expect("valid model ref")
}

fn material(namespace: &str, version: u64, byte: u8) -> MaterialCardRef {
    MaterialCardRef::new(namespace, nz(version), digest(byte)).expect("valid material ref")
}

fn interface_card(namespace: &str, version: u64, byte: u8) -> InterfaceCardRef {
    InterfaceCardRef::new(namespace, nz(version), digest(byte)).expect("valid interface ref")
}

fn policy(namespace: &str, version: u64, byte: u8) -> SolvePolicyRef {
    SolvePolicyRef::new(namespace, nz(version), digest(byte)).expect("valid policy ref")
}

fn frame(orientation: OrientationParity) -> FrameBinding {
    FrameBinding::new("world/mechanical", orientation).expect("valid frame")
}

fn pressure() -> TerminalQuantitySpec {
    TerminalQuantitySpec::Semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static))
}

fn terminal(
    key: &str,
    owner: &SubsystemId,
    quantity: TerminalQuantitySpec,
    causality: TerminalCausality,
    clock: &ClockId,
    orientation: OrientationParity,
) -> TerminalSpec {
    TerminalSpec {
        id: TerminalId::new(key).expect("valid terminal id"),
        owner: owner.clone(),
        quantity,
        shape: TerminalShape::Scalar,
        causality,
        clock: clock.clone(),
        frame: frame(orientation),
    }
}

fn rules(refusal: &MachineGraphRefusal) -> Vec<MachineGraphRule> {
    refusal
        .findings()
        .iter()
        .map(|finding| finding.rule())
        .collect()
}

// Keeping the complete fixture together makes every identity-bearing field
// visible to mutation tests without hiding defaults in nested builders.
#[allow(clippy::too_many_lines)]
fn valid_graph() -> MachineGraphDraft {
    let clock = ClockId::new("clock/mechanical").expect("valid clock id");
    let auxiliary_clock = ClockId::new("clock/auxiliary").expect("valid clock id");
    let subsystem_a = SubsystemId::new("subsystem/source").expect("valid subsystem id");
    let subsystem_b = SubsystemId::new("subsystem/load").expect("valid subsystem id");
    let body_a = BodyId::new("body/source").expect("valid body id");
    let body_b = BodyId::new("body/load").expect("valid body id");
    let body_a_aux = BodyId::new("body/source-aux").expect("valid body id");
    let body_b_aux = BodyId::new("body/load-aux").expect("valid body id");

    let a_effort = terminal(
        "terminal/source-effort",
        &subsystem_a,
        pressure(),
        TerminalCausality::Output,
        &clock,
        OrientationParity::Preserving,
    );
    let a_flow = terminal(
        "terminal/source-flow",
        &subsystem_a,
        TerminalQuantitySpec::Dimensional(VOLUME_FLOW_DIMS),
        TerminalCausality::Output,
        &clock,
        OrientationParity::Preserving,
    );
    let b_effort = terminal(
        "terminal/load-effort",
        &subsystem_b,
        pressure(),
        TerminalCausality::Input,
        &clock,
        OrientationParity::Preserving,
    );
    let b_flow = terminal(
        "terminal/load-flow",
        &subsystem_b,
        TerminalQuantitySpec::Dimensional(VOLUME_FLOW_DIMS),
        TerminalCausality::Input,
        &clock,
        OrientationParity::Preserving,
    );
    let a_aux_effort = terminal(
        "terminal/source-aux-effort",
        &subsystem_a,
        pressure(),
        TerminalCausality::Output,
        &auxiliary_clock,
        OrientationParity::Preserving,
    );
    let a_aux_flow = terminal(
        "terminal/source-aux-flow",
        &subsystem_a,
        TerminalQuantitySpec::Dimensional(VOLUME_FLOW_DIMS),
        TerminalCausality::Output,
        &auxiliary_clock,
        OrientationParity::Preserving,
    );
    let b_aux_effort = terminal(
        "terminal/load-aux-effort",
        &subsystem_b,
        pressure(),
        TerminalCausality::Input,
        &auxiliary_clock,
        OrientationParity::Preserving,
    );
    let b_aux_flow = terminal(
        "terminal/load-aux-flow",
        &subsystem_b,
        TerminalQuantitySpec::Dimensional(VOLUME_FLOW_DIMS),
        TerminalCausality::Input,
        &auxiliary_clock,
        OrientationParity::Preserving,
    );

    let port_a = PortSpec {
        id: PortId::new("port/source").expect("valid port id"),
        owner: subsystem_a.clone(),
        effort: a_effort.id.clone(),
        flow: a_flow.id.clone(),
        energy_role: PortEnergyRole::OutOfSubsystem,
    };
    let port_b = PortSpec {
        id: PortId::new("port/load").expect("valid port id"),
        owner: subsystem_b.clone(),
        effort: b_effort.id.clone(),
        flow: b_flow.id.clone(),
        energy_role: PortEnergyRole::IntoSubsystem,
    };
    let port_a_aux = PortSpec {
        id: PortId::new("port/source-aux").expect("valid port id"),
        owner: subsystem_a.clone(),
        effort: a_aux_effort.id.clone(),
        flow: a_aux_flow.id.clone(),
        energy_role: PortEnergyRole::OutOfSubsystem,
    };
    let port_b_aux = PortSpec {
        id: PortId::new("port/load-aux").expect("valid port id"),
        owner: subsystem_b.clone(),
        effort: b_aux_effort.id.clone(),
        flow: b_aux_flow.id.clone(),
        energy_role: PortEnergyRole::IntoSubsystem,
    };

    MachineGraphDraft {
        clocks: vec![
            ClockSpec {
                id: clock,
                clock: MachineClock::Periodic {
                    period_ns: nz(1_000_000),
                    phase_ns: 0,
                },
            },
            ClockSpec {
                id: auxiliary_clock,
                clock: MachineClock::Periodic {
                    period_ns: nz(4_000_000),
                    phase_ns: 1_000_000,
                },
            },
        ],
        subsystems: vec![
            SubsystemSpec {
                id: subsystem_a,
                model: model("models/source", 1, 1),
                bodies: vec![body_a.clone(), body_a_aux.clone()],
                surface_patches: vec![
                    SurfacePatchId::new("surface/source-main").expect("valid surface id"),
                    SurfacePatchId::new("surface/source-aux").expect("valid surface id"),
                ],
                contact_features: vec![
                    ContactFeatureId::new("contact/source-main").expect("valid contact id"),
                    ContactFeatureId::new("contact/source-aux").expect("valid contact id"),
                ],
                state_slots: Vec::new(),
            },
            SubsystemSpec {
                id: subsystem_b,
                model: model("models/load", 1, 2),
                bodies: vec![body_b.clone(), body_b_aux.clone()],
                surface_patches: vec![
                    SurfacePatchId::new("surface/load-main").expect("valid surface id"),
                    SurfacePatchId::new("surface/load-aux").expect("valid surface id"),
                ],
                contact_features: vec![
                    ContactFeatureId::new("contact/load-main").expect("valid contact id"),
                    ContactFeatureId::new("contact/load-aux").expect("valid contact id"),
                ],
                state_slots: Vec::new(),
            },
        ],
        terminals: vec![
            a_effort.clone(),
            a_flow.clone(),
            b_effort.clone(),
            b_flow.clone(),
            a_aux_effort.clone(),
            a_aux_flow.clone(),
            b_aux_effort.clone(),
            b_aux_flow.clone(),
        ],
        ports: vec![
            port_a.clone(),
            port_b.clone(),
            port_a_aux.clone(),
            port_b_aux.clone(),
        ],
        relations: vec![
            RelationSpec {
                id: RelationId::new("relation/effort").expect("valid relation id"),
                source: a_effort.id,
                target: b_effort.id,
                mode: RelationMode::Algebraic { solve_policy: None },
            },
            RelationSpec {
                id: RelationId::new("relation/flow").expect("valid relation id"),
                source: a_flow.id,
                target: b_flow.id,
                mode: RelationMode::Algebraic { solve_policy: None },
            },
            RelationSpec {
                id: RelationId::new("relation/aux-effort").expect("valid relation id"),
                source: a_aux_effort.id.clone(),
                target: b_aux_effort.id.clone(),
                mode: RelationMode::Algebraic { solve_policy: None },
            },
            RelationSpec {
                id: RelationId::new("relation/aux-flow").expect("valid relation id"),
                source: a_aux_flow.id.clone(),
                target: b_aux_flow.id.clone(),
                mode: RelationMode::Algebraic { solve_policy: None },
            },
        ],
        materials: vec![
            MaterialBinding {
                target: MaterialTarget::Body(body_a),
                material: material("materials/source", 1, 11),
            },
            MaterialBinding {
                target: MaterialTarget::Body(body_b),
                material: material("materials/load", 1, 12),
            },
            MaterialBinding {
                target: MaterialTarget::Body(body_a_aux),
                material: material("materials/source-aux", 1, 13),
            },
            MaterialBinding {
                target: MaterialTarget::Body(body_b_aux),
                material: material("materials/load-aux", 1, 14),
            },
        ],
        interfaces: vec![
            InterfaceBinding {
                id: InterfaceId::new("interface/source-load").expect("valid interface id"),
                negative: port_a.id,
                positive: port_b.id,
                interface: interface_card("interfaces/hydraulic", 1, 21),
                orientation: InterfaceOrientation::Aligned,
            },
            InterfaceBinding {
                id: InterfaceId::new("interface/source-load-aux").expect("valid interface id"),
                negative: port_a_aux.id,
                positive: port_b_aux.id,
                interface: interface_card("interfaces/hydraulic-aux", 1, 22),
                orientation: InterfaceOrientation::Aligned,
            },
        ],
    }
}

fn cycle_graph() -> MachineGraphDraft {
    let clock = ClockId::new("clock/control").expect("valid clock id");
    let subsystem_a = SubsystemId::new("subsystem/a").expect("valid subsystem id");
    let subsystem_b = SubsystemId::new("subsystem/b").expect("valid subsystem id");
    let dims = TerminalQuantitySpec::Dimensional(Dims::NONE);
    let a_out = terminal(
        "terminal/a-out",
        &subsystem_a,
        dims,
        TerminalCausality::Output,
        &clock,
        OrientationParity::Preserving,
    );
    let a_in = terminal(
        "terminal/a-in",
        &subsystem_a,
        dims,
        TerminalCausality::Input,
        &clock,
        OrientationParity::Preserving,
    );
    let b_out = terminal(
        "terminal/b-out",
        &subsystem_b,
        dims,
        TerminalCausality::Output,
        &clock,
        OrientationParity::Preserving,
    );
    let b_in = terminal(
        "terminal/b-in",
        &subsystem_b,
        dims,
        TerminalCausality::Input,
        &clock,
        OrientationParity::Preserving,
    );

    MachineGraphDraft {
        clocks: vec![ClockSpec {
            id: clock,
            clock: MachineClock::Continuous,
        }],
        subsystems: vec![
            SubsystemSpec {
                id: subsystem_a,
                model: model("models/a", 1, 31),
                bodies: Vec::new(),
                surface_patches: Vec::new(),
                contact_features: Vec::new(),
                state_slots: Vec::new(),
            },
            SubsystemSpec {
                id: subsystem_b,
                model: model("models/b", 1, 32),
                bodies: Vec::new(),
                surface_patches: Vec::new(),
                contact_features: Vec::new(),
                state_slots: Vec::new(),
            },
        ],
        terminals: vec![a_out.clone(), a_in.clone(), b_out.clone(), b_in.clone()],
        ports: Vec::new(),
        relations: vec![
            RelationSpec {
                id: RelationId::new("relation/a-to-b").expect("valid relation id"),
                source: a_out.id,
                target: b_in.id,
                mode: RelationMode::Algebraic { solve_policy: None },
            },
            RelationSpec {
                id: RelationId::new("relation/b-to-a").expect("valid relation id"),
                source: b_out.id,
                target: a_in.id,
                mode: RelationMode::Algebraic { solve_policy: None },
            },
        ],
        materials: Vec::new(),
        interfaces: Vec::new(),
    }
}

#[test]
fn g0_valid_graph_admits_with_structured_decision() {
    let decision = valid_graph().admit_with_decision();
    assert_eq!(decision.code(), "MachineGraphAdmitted");
    assert_eq!(decision.submitted_counts().subsystems, 2);
    assert_eq!(decision.submitted_counts().terminals, 8);
    assert_eq!(decision.submitted_counts().ports, 4);
    assert_eq!(decision.submitted_counts().relations, 4);
    assert!(decision.result().is_ok());
}

#[test]
fn g3_permuting_every_collection_preserves_graph_identity() {
    let expected = valid_graph().admit().expect("baseline graph admits");
    let mut permuted = valid_graph();
    permuted.clocks.reverse();
    permuted.subsystems.reverse();
    for subsystem in &mut permuted.subsystems {
        subsystem.bodies.reverse();
        subsystem.surface_patches.reverse();
        subsystem.contact_features.reverse();
        subsystem.state_slots.reverse();
    }
    permuted.terminals.reverse();
    permuted.ports.reverse();
    permuted.relations.reverse();
    permuted.materials.reverse();
    permuted.interfaces.reverse();
    let actual = permuted.admit().expect("permuted graph admits");
    assert_eq!(actual.identity(), expected.identity());
    assert_eq!(actual.identity_receipt(), expected.identity_receipt());
}

#[test]
fn g0_equal_dimensions_do_not_erase_semantic_kind_or_affine_role() {
    let mut pressure_vs_stress = valid_graph();
    pressure_vs_stress.terminals[2].quantity =
        TerminalQuantitySpec::Semantic(SemanticType::new(QuantityKind::Stress, ValueForm::Static));
    let refusal = pressure_vs_stress
        .admit()
        .expect_err("pressure and stress must not connect");
    assert!(rules(&refusal).contains(&MachineGraphRule::RelationQuantityGap));

    let mut absolute_vs_delta = valid_graph();
    absolute_vs_delta.terminals[0].quantity = TerminalQuantitySpec::Semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    absolute_vs_delta.terminals[2].quantity = TerminalQuantitySpec::Semantic(SemanticType::new(
        QuantityKind::TemperatureDifference,
        ValueForm::Static,
    ));
    let refusal = absolute_vs_delta
        .admit()
        .expect_err("absolute and difference temperature must not connect");
    assert!(rules(&refusal).contains(&MachineGraphRule::RelationQuantityGap));
}

#[test]
fn g0_semantic_form_and_periodic_phase_refuse_before_identity() {
    let mut bad_form = valid_graph();
    bad_form.terminals[0].quantity = TerminalQuantitySpec::Semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Instantaneous,
    ));
    let decision = bad_form.admit_with_decision();
    let refusal = decision.result().expect_err("invalid form refuses");
    assert_eq!(decision.code(), "MachineGraphRefused");
    assert!(rules(refusal).contains(&MachineGraphRule::UnsupportedTerminalForm));

    let mut bad_phase = valid_graph();
    bad_phase.clocks[0].clock = MachineClock::Periodic {
        period_ns: nz(10),
        phase_ns: 10,
    };
    let refusal = bad_phase
        .admit()
        .expect_err("phase equal to period refuses");
    assert!(rules(&refusal).contains(&MachineGraphRule::InvalidClockPhase));
}

#[test]
fn g0_relation_checks_causality_clock_frame_and_orientation() {
    let mut causality = valid_graph();
    causality.terminals[2].causality = TerminalCausality::Output;
    let refusal = causality.admit().expect_err("output-to-output refuses");
    assert!(rules(&refusal).contains(&MachineGraphRule::RelationCausalityGap));

    let mut clock = valid_graph();
    let slow = ClockId::new("clock/slow").expect("valid clock id");
    clock.clocks.push(ClockSpec {
        id: slow.clone(),
        clock: MachineClock::Periodic {
            period_ns: nz(2_000_000),
            phase_ns: 0,
        },
    });
    clock.terminals[2].clock = slow.clone();
    clock.terminals[3].clock = slow;
    let refusal = clock.admit().expect_err("clock gap refuses");
    assert!(rules(&refusal).contains(&MachineGraphRule::RelationClockGap));

    let mut frame_gap = valid_graph();
    frame_gap.terminals[2].frame =
        FrameBinding::new("world/load", OrientationParity::Preserving).expect("valid frame");
    frame_gap.terminals[3].frame =
        FrameBinding::new("world/load", OrientationParity::Preserving).expect("valid frame");
    let refusal = frame_gap.admit().expect_err("frame gap refuses");
    assert!(rules(&refusal).contains(&MachineGraphRule::RelationFrameGap));

    let mut orientation = valid_graph();
    orientation.terminals[2].frame = frame(OrientationParity::Reversing);
    orientation.terminals[3].frame = frame(OrientationParity::Reversing);
    orientation.interfaces.remove(0);
    let refusal = orientation
        .admit()
        .expect_err("unmediated relation orientation gap remains explicit");
    assert!(rules(&refusal).contains(&MachineGraphRule::RelationOrientationGap));
}

#[test]
fn g0_source_closure_refuses_missing_and_multiple_producers() {
    let mut missing = valid_graph();
    missing.relations.remove(0);
    let refusal = missing.admit().expect_err("unclosed input refuses");
    assert!(rules(&refusal).contains(&MachineGraphRule::MissingSourceClosure));

    let mut multiple = valid_graph();
    let mut duplicate = multiple.relations[0].clone();
    duplicate.id = RelationId::new("relation/effort-backup").expect("valid relation id");
    multiple.relations.push(duplicate);
    let refusal = multiple.admit().expect_err("multiple sources refuse");
    assert!(rules(&refusal).contains(&MachineGraphRule::MultipleInputSources));
}

#[test]
fn g0_explicit_external_input_needs_no_internal_relation() {
    let mut graph = cycle_graph();
    let removed = graph.relations.remove(0);
    let target = graph
        .terminals
        .iter_mut()
        .find(|terminal| terminal.id == removed.target)
        .expect("removed relation target exists");
    target.causality = TerminalCausality::ExternalInput;
    graph
        .admit()
        .expect("explicit boundary input is source-closed by declaration");
}

#[test]
fn g0_ownership_clock_port_and_interface_references_fail_closed() {
    let mut ownership = valid_graph();
    let shared_body = ownership.subsystems[0].bodies[0].clone();
    ownership.subsystems[1].bodies.push(shared_body);
    let refusal = ownership
        .admit()
        .expect_err("one durable body cannot have two owners");
    assert!(rules(&refusal).contains(&MachineGraphRule::DuplicateOwnership));

    let mut clock = valid_graph();
    clock.terminals[0].clock = ClockId::new("clock/missing").expect("valid clock id");
    let refusal = clock.admit().expect_err("dangling terminal clock refuses");
    assert!(rules(&refusal).contains(&MachineGraphRule::UnknownTerminalClock));

    let mut power = valid_graph();
    power.terminals[1].quantity = TerminalQuantitySpec::Dimensional(Dims::NONE);
    let refusal = power.admit().expect_err("non-power port refuses");
    assert!(rules(&refusal).contains(&MachineGraphRule::PortPowerDimensionMismatch));

    let mut shape = valid_graph();
    shape.terminals[1].shape = TerminalShape::Vector { components: nz(2) };
    let refusal = shape.admit().expect_err("effort/flow shape gap refuses");
    assert!(rules(&refusal).contains(&MachineGraphRule::PortShapeMismatch));

    let mut interface = valid_graph();
    interface.interfaces[0].positive = PortId::new("port/missing").expect("valid port id");
    let refusal = interface
        .admit()
        .expect_err("dangling interface port refuses");
    assert!(rules(&refusal).contains(&MachineGraphRule::UnknownInterfacePort));

    let mut causality = valid_graph();
    causality.terminals[2].causality = TerminalCausality::Output;
    causality.terminals[3].causality = TerminalCausality::Output;
    let refusal = causality
        .admit()
        .expect_err("output-to-output interface refuses");
    assert!(rules(&refusal).contains(&MachineGraphRule::InterfaceCausalityGap));

    let mut relation = valid_graph();
    relation.relations.remove(0);
    let refusal = relation
        .admit()
        .expect_err("interface without exact relation refuses");
    assert!(rules(&refusal).contains(&MachineGraphRule::InterfaceRelationGap));
}

#[test]
fn g3_refusal_findings_are_permutation_invariant() {
    let mut first = valid_graph();
    first.relations.remove(0);
    first.materials.push(first.materials[0].clone());
    let mut second = first.clone();
    second.clocks.reverse();
    second.subsystems.reverse();
    second.terminals.reverse();
    second.ports.reverse();
    second.relations.reverse();
    second.materials.reverse();
    second.interfaces.reverse();

    let first = first.admit().expect_err("fixture intentionally refuses");
    let second = second.admit().expect_err("fixture intentionally refuses");
    assert_eq!(first.findings(), second.findings());

    let mut conflicting_first = valid_graph();
    let mut conflicting_terminal = conflicting_first.terminals[0].clone();
    conflicting_terminal.quantity =
        TerminalQuantitySpec::Semantic(SemanticType::new(QuantityKind::Stress, ValueForm::Static));
    conflicting_first.terminals.push(conflicting_terminal);
    let mut conflicting_second = conflicting_first.clone();
    conflicting_second.terminals.reverse();
    let conflicting_first = conflicting_first
        .admit()
        .expect_err("conflicting duplicate terminal refuses");
    let conflicting_second = conflicting_second
        .admit()
        .expect_err("permuted conflicting duplicate terminal refuses");
    assert_eq!(conflicting_first.findings(), conflicting_second.findings());
}

#[test]
fn g0_algebraic_cycle_requires_explicit_policy_boundary() {
    let refusal = cycle_graph()
        .admit()
        .expect_err("unguarded algebraic cycle refuses");
    assert!(rules(&refusal).contains(&MachineGraphRule::AlgebraicLoopWithoutSolvePolicy));

    let mut governed = cycle_graph();
    governed.relations[0].mode = RelationMode::Algebraic {
        solve_policy: Some(policy("policies/newton", 1, 41)),
    };
    governed
        .admit()
        .expect("one explicit solve boundary breaks the ungoverned cycle");
}

#[test]
fn g0_stateful_dependency_breaks_cycle_and_accounts_for_state() {
    let mut unaccounted = cycle_graph();
    let state = StateSlotId::new("state/b-delay").expect("valid state id");
    unaccounted.subsystems[1].state_slots.push(state.clone());
    let refusal = unaccounted
        .clone()
        .admit()
        .expect_err("unwritten declared state refuses");
    assert!(rules(&refusal).contains(&MachineGraphRule::UnaccountedState));

    unaccounted.relations[0].mode = RelationMode::Stateful { state_slot: state };
    unaccounted
        .admit()
        .expect("stateful edge breaks feedthrough and names one writer");
}

#[test]
fn g0_state_owner_and_writer_conflicts_are_structured() {
    let mut graph = cycle_graph();
    let state = StateSlotId::new("state/a-delay").expect("valid state id");
    graph.subsystems[0].state_slots.push(state.clone());
    graph.relations[0].mode = RelationMode::Stateful {
        state_slot: state.clone(),
    };
    let refusal = graph
        .admit()
        .expect_err("target cannot write foreign state");
    assert!(rules(&refusal).contains(&MachineGraphRule::StateOwnerMismatch));

    let mut multiple = cycle_graph();
    let state = StateSlotId::new("state/b-delay").expect("valid state id");
    multiple.subsystems[1].state_slots.push(state.clone());
    multiple.relations[0].mode = RelationMode::Stateful {
        state_slot: state.clone(),
    };
    let mut second = multiple.relations[0].clone();
    second.id = RelationId::new("relation/a-to-b-backup").expect("valid relation id");
    second.mode = RelationMode::Stateful { state_slot: state };
    multiple.relations.push(second);
    let refusal = multiple.admit().expect_err("multiple writers refuse");
    assert!(rules(&refusal).contains(&MachineGraphRule::MultipleStateWriters));
}

#[test]
fn g0_material_and_interface_closure_refuse_conflicts() {
    let mut missing_material = valid_graph();
    missing_material.materials.remove(0);
    let refusal = missing_material
        .admit()
        .expect_err("every owned body needs material binding");
    assert!(rules(&refusal).contains(&MachineGraphRule::MissingBodyMaterial));

    let mut duplicate_material = valid_graph();
    duplicate_material
        .materials
        .push(duplicate_material.materials[0].clone());
    let refusal = duplicate_material
        .admit()
        .expect_err("duplicate material target refuses");
    assert!(rules(&refusal).contains(&MachineGraphRule::DuplicateMaterialBinding));

    let mut orientation = valid_graph();
    orientation.interfaces[0].orientation = InterfaceOrientation::Opposed;
    let refusal = orientation
        .admit()
        .expect_err("declared opposed orientation must be witnessed");
    assert!(rules(&refusal).contains(&MachineGraphRule::InterfaceOrientationGap));

    let mut energy = valid_graph();
    energy.ports[1].energy_role = PortEnergyRole::OutOfSubsystem;
    let refusal = energy
        .admit()
        .expect_err("interface energy roles must complement");
    assert!(rules(&refusal).contains(&MachineGraphRule::InterfaceEnergyRoleGap));
}

#[test]
fn g0_interface_admits_conjugate_effort_flow_causality() {
    let mut graph = valid_graph();
    graph.terminals[1].causality = TerminalCausality::Input;
    graph.terminals[3].causality = TerminalCausality::Output;
    let source = graph.terminals[3].id.clone();
    let target = graph.terminals[1].id.clone();
    let flow_relation = graph
        .relations
        .iter_mut()
        .find(|relation| relation.id.canonical_key() == "relation/flow")
        .expect("flow relation exists");
    flow_relation.source = source;
    flow_relation.target = target;
    flow_relation.mode = RelationMode::Algebraic {
        solve_policy: Some(policy("policies/interface-coupling", 1, 66)),
    };
    graph
        .admit()
        .expect("each conjugate pair has one producer and one consumer");
}

#[test]
#[allow(clippy::too_many_lines)]
fn g3_semantic_inputs_move_identity_without_caller_order_leaks() {
    let baseline = valid_graph().admit().expect("baseline admits").identity();

    let mut model_version = valid_graph();
    model_version.subsystems[0].model = model("models/source", 2, 1);
    assert_ne!(
        model_version
            .admit()
            .expect("model mutation admits")
            .identity(),
        baseline
    );

    let mut clock_period = valid_graph();
    clock_period.clocks[0].clock = MachineClock::Periodic {
        period_ns: nz(2_000_000),
        phase_ns: 0,
    };
    assert_ne!(
        clock_period
            .admit()
            .expect("clock mutation admits")
            .identity(),
        baseline
    );

    let mut semantic_kind = valid_graph();
    semantic_kind.terminals[0].quantity =
        TerminalQuantitySpec::Semantic(SemanticType::new(QuantityKind::Stress, ValueForm::Static));
    semantic_kind.terminals[2].quantity = semantic_kind.terminals[0].quantity;
    assert_ne!(
        semantic_kind
            .admit()
            .expect("compatible semantic mutation admits")
            .identity(),
        baseline
    );

    let mut explicit_no_kind = valid_graph();
    explicit_no_kind.terminals[0].quantity = TerminalQuantitySpec::Dimensional(PRESSURE_DIMS);
    explicit_no_kind.terminals[2].quantity = TerminalQuantitySpec::Dimensional(PRESSURE_DIMS);
    assert_ne!(
        explicit_no_kind
            .admit()
            .expect("dimension-only mutation admits")
            .identity(),
        baseline
    );

    let mut solve_policy = valid_graph();
    solve_policy.relations[0].mode = RelationMode::Algebraic {
        solve_policy: Some(policy("policies/direct", 1, 42)),
    };
    let first_policy_identity = solve_policy
        .admit()
        .expect("policy mutation admits")
        .identity();
    assert_ne!(first_policy_identity, baseline);
    let mut other_policy = valid_graph();
    other_policy.relations[0].mode = RelationMode::Algebraic {
        solve_policy: Some(policy("policies/direct", 2, 43)),
    };
    assert_ne!(
        other_policy
            .admit()
            .expect("second policy mutation admits")
            .identity(),
        first_policy_identity
    );

    let mut material_card = valid_graph();
    material_card.materials[0].material = material("materials/source", 1, 77);
    assert_ne!(
        material_card
            .admit()
            .expect("material-card mutation admits")
            .identity(),
        baseline
    );

    let mut interface_card_graph = valid_graph();
    interface_card_graph.interfaces[0].interface = interface_card("interfaces/hydraulic", 1, 78);
    assert_ne!(
        interface_card_graph
            .admit()
            .expect("interface-card mutation admits")
            .identity(),
        baseline
    );

    let mut oriented = valid_graph();
    oriented.terminals[2].frame = frame(OrientationParity::Reversing);
    oriented.terminals[3].frame = frame(OrientationParity::Reversing);
    oriented.interfaces[0].orientation = InterfaceOrientation::Opposed;
    assert_ne!(
        oriented
            .admit()
            .expect("opposed interface relation admits")
            .identity(),
        baseline
    );
}

#[test]
fn g0_public_resource_envelope_refuses_before_graph_work() {
    let mut graph = valid_graph();
    let clock = graph.clocks[0].clone();
    graph.clocks = vec![clock; MAX_MACHINE_GRAPH_CLOCKS + 1];
    let decision = graph.admit_with_decision();
    assert_eq!(
        decision.submitted_counts().clocks,
        MAX_MACHINE_GRAPH_CLOCKS + 1
    );
    let refusal = decision.result().expect_err("clock limit refuses");
    assert_eq!(refusal.findings().len(), 1);
    assert_eq!(
        refusal.findings()[0].rule(),
        MachineGraphRule::ResourceLimit
    );
}

#[test]
fn g0_external_references_and_graph_ids_are_role_separated() {
    assert!(ModelRef::new("models/a", nz(1), [0; 32]).is_err());
    assert!(MaterialCardRef::new("Materials/A", nz(1), digest(1)).is_err());

    let subsystem = SubsystemId::new("shared/key").expect("valid subsystem id");
    let clock = ClockId::new("shared/key").expect("valid clock id");
    let relation = RelationId::new("shared/key").expect("valid relation id");
    let interface = InterfaceId::new("shared/key").expect("valid interface id");
    assert_ne!(subsystem.identity().as_bytes(), clock.identity().as_bytes());
    assert_ne!(
        subsystem.identity().as_bytes(),
        relation.identity().as_bytes()
    );
    assert_ne!(clock.identity().as_bytes(), relation.identity().as_bytes());
    assert_ne!(
        interface.identity().as_bytes(),
        subsystem.identity().as_bytes()
    );
    assert_ne!(interface.identity().as_bytes(), clock.identity().as_bytes());
    assert_ne!(
        interface.identity().as_bytes(),
        relation.identity().as_bytes()
    );
}
