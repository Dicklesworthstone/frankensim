//! Machine-IR E0 PR-3 behavior admission conformance (Gauntlet G0/G3).

use core::num::NonZeroU64;

use fs_blake3::identity::StrongIdentity;
use fs_ir::machine::semantics::{
    BodyMotion, ConditionBinding, ConditionHistoryRef, ConditionSource, ConditionTarget,
    ConditionValueRef, CorrelationModelRef, CrossingSemantics, DependenceMember, DependenceModel,
    DependenceSpec, DistributionRef, EventDependency, EventId, EventOrder, EventSpec,
    FiniteNonNegative, GuardOrientation, GuardRef, HistoryContinuity,
    MAX_MACHINE_BEHAVIOR_NESTED_REFERENCES, MAX_MACHINE_BEHAVIOR_STATE_CONTRACTS,
    MachineBehaviorDraft, MachineBehaviorRefusal, MachineBehaviorRule, MotionBinding,
    MotionPathRef, NoClaimRef, OutcomeSetRef, ParameterRef, ResetMapRef, ResetSemantics,
    SimultaneityGroupRef, StateSlotContract, ToleranceId, ToleranceLawRef, ToleranceSemantics,
    ToleranceSpec, ToleranceTarget,
};
use fs_ir::machine::{
    BodyId, ClockId, ClockSpec, FrameBinding, MachineClock, MachineGraphDraft, MaterialBinding,
    MaterialCardRef, MaterialTarget, ModelRef, OrientationParity, RelationId, RelationMode,
    RelationSpec, StateSlotId, SubsystemId, SubsystemSpec, TerminalCausality, TerminalId,
    TerminalQuantitySpec, TerminalShape, TerminalSpec,
};
use fs_qty::Dims;

fn nz(value: u64) -> NonZeroU64 {
    NonZeroU64::new(value).expect("test value is nonzero")
}

fn digest(byte: u8) -> [u8; 32] {
    [byte; 32]
}

fn frame() -> FrameBinding {
    FrameBinding::new("world/mechanical", OrientationParity::Preserving).expect("valid frame")
}

fn quantity() -> TerminalQuantitySpec {
    TerminalQuantitySpec::Dimensional(Dims::NONE)
}

fn model(byte: u8) -> ModelRef {
    ModelRef::new("models/behavior-fixture", nz(1), digest(byte)).expect("valid model ref")
}

fn material(byte: u8) -> MaterialCardRef {
    MaterialCardRef::new("materials/behavior-fixture", nz(1), digest(byte))
        .expect("valid material ref")
}

fn condition_value(byte: u8) -> ConditionValueRef {
    ConditionValueRef::new("values/condition", nz(1), digest(byte)).expect("valid value ref")
}

fn history(byte: u8) -> ConditionHistoryRef {
    ConditionHistoryRef::new("signals/history", nz(1), digest(byte)).expect("valid history ref")
}

fn distribution(byte: u8) -> DistributionRef {
    DistributionRef::new("distributions/marginal", nz(1), digest(byte))
        .expect("valid distribution ref")
}

fn guard(byte: u8) -> GuardRef {
    GuardRef::new("guards/state-threshold", nz(1), digest(byte)).expect("valid guard ref")
}

fn no_claim(byte: u8) -> NoClaimRef {
    NoClaimRef::new("claims/crossing-unknown", nz(1), digest(byte)).expect("valid no-claim ref")
}

fn reset(byte: u8) -> ResetMapRef {
    ResetMapRef::new("resets/state-map", nz(1), digest(byte)).expect("valid reset ref")
}

fn motion_path(byte: u8) -> MotionPathRef {
    MotionPathRef::new("motions/body-path", nz(1), digest(byte)).expect("valid motion ref")
}

fn outcomes(byte: u8) -> OutcomeSetRef {
    OutcomeSetRef::new("resets/outcome-set", nz(1), digest(byte)).expect("valid outcome ref")
}

fn simultaneity(byte: u8) -> SimultaneityGroupRef {
    SimultaneityGroupRef::new("events/simultaneity", nz(1), digest(byte))
        .expect("valid simultaneity ref")
}

fn parameter(name: &str, byte: u8) -> ParameterRef {
    ParameterRef::new(name, nz(1), digest(byte)).expect("valid parameter ref")
}

fn tolerance_law(byte: u8) -> ToleranceLawRef {
    ToleranceLawRef::new("tolerances/additive", nz(1), digest(byte)).expect("valid tolerance law")
}

fn correlation(byte: u8) -> CorrelationModelRef {
    CorrelationModelRef::new("correlations/joint-law", nz(1), digest(byte))
        .expect("valid correlation model")
}

fn scalar(value: f64) -> FiniteNonNegative {
    FiniteNonNegative::new(value).expect("valid nonnegative finite scalar")
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
        quantity: quantity(),
        shape: TerminalShape::Scalar,
        causality,
        clock: clock.clone(),
        frame: frame(),
    }
}

#[allow(clippy::too_many_lines)]
fn valid_graph() -> MachineGraphDraft {
    let continuous = ClockId::new("clock/continuous").expect("valid clock id");
    let event = ClockId::new("clock/events").expect("valid clock id");
    let subsystem = SubsystemId::new("subsystem/plant").expect("valid subsystem id");
    let state_a = StateSlotId::new("state/position").expect("valid state id");
    let state_b = StateSlotId::new("state/velocity").expect("valid state id");
    let body = BodyId::new("body/plant").expect("valid body id");

    let output_a = terminal(
        "terminal/state-a-source",
        &subsystem,
        TerminalCausality::Output,
        &continuous,
    );
    let input_a = terminal(
        "terminal/state-a-sink",
        &subsystem,
        TerminalCausality::Input,
        &continuous,
    );
    let output_b = terminal(
        "terminal/state-b-source",
        &subsystem,
        TerminalCausality::Output,
        &continuous,
    );
    let input_b = terminal(
        "terminal/state-b-sink",
        &subsystem,
        TerminalCausality::Input,
        &continuous,
    );
    let boundary = terminal(
        "terminal/external-command",
        &subsystem,
        TerminalCausality::ExternalInput,
        &continuous,
    );
    let guard_output = terminal(
        "terminal/guard-observation",
        &subsystem,
        TerminalCausality::Output,
        &continuous,
    );

    MachineGraphDraft {
        clocks: vec![
            ClockSpec {
                id: continuous,
                clock: MachineClock::Continuous,
            },
            ClockSpec {
                id: event,
                clock: MachineClock::EventDriven,
            },
        ],
        subsystems: vec![SubsystemSpec {
            id: subsystem,
            model: model(1),
            bodies: vec![body.clone()],
            surface_patches: Vec::new(),
            contact_features: Vec::new(),
            state_slots: vec![state_a.clone(), state_b.clone()],
        }],
        terminals: vec![
            output_a.clone(),
            input_a.clone(),
            output_b.clone(),
            input_b.clone(),
            boundary,
            guard_output,
        ],
        ports: Vec::new(),
        relations: vec![
            RelationSpec {
                id: RelationId::new("relation/state-a").expect("valid relation id"),
                source: output_a.id,
                target: input_a.id,
                mode: RelationMode::Stateful {
                    state_slot: state_a,
                },
            },
            RelationSpec {
                id: RelationId::new("relation/state-b").expect("valid relation id"),
                source: output_b.id,
                target: input_b.id,
                mode: RelationMode::Stateful {
                    state_slot: state_b,
                },
            },
        ],
        materials: vec![MaterialBinding {
            target: MaterialTarget::Body(body),
            material: material(2),
        }],
        interfaces: Vec::new(),
    }
}

#[allow(clippy::too_many_lines)]
fn valid_behavior() -> MachineBehaviorDraft {
    let subsystem = SubsystemId::new("subsystem/plant").expect("valid subsystem id");
    let continuous = ClockId::new("clock/continuous").expect("valid clock id");
    let event_clock = ClockId::new("clock/events").expect("valid clock id");
    let state_a = StateSlotId::new("state/position").expect("valid state id");
    let state_b = StateSlotId::new("state/velocity").expect("valid state id");
    let body = BodyId::new("body/plant").expect("valid body id");
    let boundary = TerminalId::new("terminal/external-command").expect("valid terminal id");
    let guard_terminal = TerminalId::new("terminal/guard-observation").expect("valid terminal id");
    let event = EventId::new("event/contact").expect("valid event id");
    let random_a = ToleranceId::new("tolerance/position").expect("valid tolerance id");
    let random_b = ToleranceId::new("tolerance/velocity").expect("valid tolerance id");

    MachineBehaviorDraft {
        state_contracts: vec![
            StateSlotContract {
                id: state_a.clone(),
                owner: subsystem.clone(),
                quantity: quantity(),
                shape: TerminalShape::Scalar,
                clock: continuous.clone(),
                frame: frame(),
            },
            StateSlotContract {
                id: state_b.clone(),
                owner: subsystem,
                quantity: quantity(),
                shape: TerminalShape::Scalar,
                clock: continuous.clone(),
                frame: frame(),
            },
        ],
        conditions: vec![
            ConditionBinding {
                target: ConditionTarget::Initial(state_a.clone()),
                quantity: quantity(),
                shape: TerminalShape::Scalar,
                clock: continuous.clone(),
                frame: frame(),
                source: ConditionSource::Fixed(condition_value(10)),
            },
            ConditionBinding {
                target: ConditionTarget::Initial(state_b.clone()),
                quantity: quantity(),
                shape: TerminalShape::Scalar,
                clock: continuous.clone(),
                frame: frame(),
                source: ConditionSource::Distribution(distribution(11)),
            },
            ConditionBinding {
                target: ConditionTarget::Boundary(boundary.clone()),
                quantity: quantity(),
                shape: TerminalShape::Scalar,
                clock: continuous.clone(),
                frame: frame(),
                source: ConditionSource::History {
                    history: history(12),
                    continuity: HistoryContinuity::ResetAtEvents {
                        events: vec![event.clone()],
                    },
                },
            },
        ],
        motions: vec![MotionBinding {
            body: body.clone(),
            clock: continuous,
            reference_frame: frame(),
            motion: BodyMotion::Static,
        }],
        events: vec![EventSpec {
            id: event,
            clock: event_clock,
            guard: guard(20),
            orientation: GuardOrientation::NegativeToPositive,
            crossing: CrossingSemantics::Unknown(no_claim(21)),
            dependencies: vec![
                EventDependency::State(state_a.clone()),
                EventDependency::Terminal(guard_terminal),
            ],
            reset: ResetSemantics::Deterministic {
                map: reset(22),
                writes: vec![state_a.clone(), state_b.clone()],
            },
            order: EventOrder::TotalPriority { microstep: 0 },
        }],
        tolerances: vec![
            ToleranceSpec {
                id: random_a.clone(),
                target: ToleranceTarget::Element(state_a.into()),
                parameter: parameter("parameters/position", 30),
                quantity: quantity(),
                shape: TerminalShape::Scalar,
                semantics: ToleranceSemantics::Random {
                    scale: scalar(0.1),
                    law: tolerance_law(31),
                    marginal: distribution(32),
                },
            },
            ToleranceSpec {
                id: random_b.clone(),
                target: ToleranceTarget::Element(state_b.clone().into()),
                parameter: parameter("parameters/velocity", 33),
                quantity: quantity(),
                shape: TerminalShape::Scalar,
                semantics: ToleranceSemantics::Random {
                    scale: scalar(0.2),
                    law: tolerance_law(34),
                    marginal: distribution(35),
                },
            },
            ToleranceSpec {
                id: ToleranceId::new("tolerance/body-clearance").expect("valid tolerance id"),
                target: ToleranceTarget::Element(body.into()),
                parameter: parameter("parameters/body-clearance", 36),
                quantity: quantity(),
                shape: TerminalShape::Scalar,
                semantics: ToleranceSemantics::Bounded {
                    minus: scalar(0.01),
                    plus: scalar(0.02),
                    law: tolerance_law(37),
                },
            },
        ],
        dependences: vec![DependenceSpec {
            members: vec![
                DependenceMember::Condition(ConditionTarget::Initial(state_b)),
                DependenceMember::Tolerance(random_a),
                DependenceMember::Tolerance(random_b),
            ],
            model: DependenceModel::Correlated(correlation(40)),
        }],
    }
}

fn permutation_fixture() -> (MachineGraphDraft, MachineBehaviorDraft) {
    let mut graph = valid_graph();
    let second_body = BodyId::new("body/plant-aux").expect("valid body id");
    graph.subsystems[0].bodies.push(second_body.clone());
    graph.materials.push(MaterialBinding {
        target: MaterialTarget::Body(second_body.clone()),
        material: material(70),
    });

    let mut behavior = valid_behavior();
    behavior.motions.push(MotionBinding {
        body: second_body,
        clock: ClockId::new("clock/continuous").expect("valid clock id"),
        reference_frame: frame(),
        motion: BodyMotion::Prescribed {
            path: motion_path(71),
        },
    });

    let mut second_event = behavior.events[0].clone();
    second_event.id = EventId::new("event/contact-aux").expect("valid event id");
    second_event.guard = guard(72);
    second_event.reset = ResetSemantics::Terminal {
        relation: reset(73),
    };
    second_event.order = EventOrder::TotalPriority { microstep: 1 };
    let second_event_id = second_event.id.clone();
    behavior.events.push(second_event);
    let boundary = behavior
        .conditions
        .iter_mut()
        .find(|condition| matches!(&condition.target, ConditionTarget::Boundary(_)))
        .expect("fixture has one boundary condition");
    let ConditionSource::History {
        continuity: HistoryContinuity::ResetAtEvents { events },
        ..
    } = &mut boundary.source
    else {
        panic!("fixture boundary is a reset-delimited history");
    };
    events.push(second_event_id);

    let boundary_terminal =
        TerminalId::new("terminal/external-command").expect("valid terminal id");
    let random_boundary =
        ToleranceId::new("tolerance/external-command").expect("valid tolerance id");
    behavior.tolerances.push(ToleranceSpec {
        id: random_boundary.clone(),
        target: ToleranceTarget::Element(boundary_terminal.into()),
        parameter: parameter("parameters/external-command", 74),
        quantity: quantity(),
        shape: TerminalShape::Scalar,
        semantics: ToleranceSemantics::Random {
            scale: scalar(0.3),
            law: tolerance_law(75),
            marginal: distribution(76),
        },
    });
    behavior.dependences[0]
        .members
        .push(DependenceMember::Tolerance(random_boundary));

    (graph, behavior)
}

fn rules(refusal: &MachineBehaviorRefusal) -> Vec<MachineBehaviorRule> {
    refusal
        .findings()
        .iter()
        .map(|finding| finding.rule())
        .collect()
}

#[test]
fn g0_fully_populated_behavior_admits_with_structured_decision() {
    let graph = valid_graph().admit().expect("base graph admits");
    let decision = valid_behavior().admit_with_decision(&graph);
    assert_eq!(decision.code(), "MachineBehaviorAdmitted");
    assert_eq!(decision.submitted_counts().state_contracts, 2);
    assert_eq!(decision.submitted_counts().conditions, 3);
    assert_eq!(decision.submitted_counts().events, 1);
    assert_eq!(decision.submitted_counts().tolerances, 3);
    let admitted = decision.result().expect("behavior admits");
    assert_eq!(admitted.base_graph(), graph.identity());
}

#[test]
fn g3_collection_and_nested_order_do_not_leak_into_identity() {
    let (graph, behavior) = permutation_fixture();
    let graph = graph.admit().expect("base graph admits");
    let expected = behavior
        .clone()
        .admit_against(&graph)
        .expect("baseline behavior admits");
    let mut permuted = behavior;
    permuted.state_contracts.reverse();
    permuted.conditions.reverse();
    for condition in &mut permuted.conditions {
        if let ConditionSource::History {
            continuity: HistoryContinuity::ResetAtEvents { events },
            ..
        } = &mut condition.source
        {
            events.reverse();
        }
    }
    permuted.motions.reverse();
    permuted.events.reverse();
    for event in &mut permuted.events {
        event.dependencies.reverse();
        if let ResetSemantics::Deterministic { writes, .. } = &mut event.reset {
            writes.reverse();
        }
    }
    permuted.tolerances.reverse();
    permuted.dependences[0].members.reverse();
    let actual = permuted
        .admit_against(&graph)
        .expect("permuted behavior admits");
    assert_eq!(actual.identity(), expected.identity());
    assert_eq!(actual.identity_receipt(), expected.identity_receipt());
}

#[test]
fn g0_state_and_boundary_closure_fail_closed() {
    let graph = valid_graph().admit().expect("base graph admits");

    let mut missing_contract = valid_behavior();
    missing_contract.state_contracts.pop();
    let refusal = missing_contract
        .admit_against(&graph)
        .expect_err("missing state contract refuses");
    assert!(rules(&refusal).contains(&MachineBehaviorRule::MissingStateContract));

    let mut missing_initial = valid_behavior();
    missing_initial.conditions.remove(1);
    let refusal = missing_initial
        .admit_against(&graph)
        .expect_err("missing initial condition refuses");
    assert!(rules(&refusal).contains(&MachineBehaviorRule::MissingInitialCondition));

    let mut missing_boundary = valid_behavior();
    missing_boundary
        .conditions
        .retain(|condition| !matches!(&condition.target, ConditionTarget::Boundary(_)));
    let refusal = missing_boundary
        .admit_against(&graph)
        .expect_err("missing external source refuses");
    assert!(rules(&refusal).contains(&MachineBehaviorRule::MissingBoundaryCondition));

    let mut duplicate = valid_behavior();
    duplicate.conditions.push(duplicate.conditions[0].clone());
    let refusal = duplicate
        .admit_against(&graph)
        .expect_err("duplicate initial source refuses");
    assert!(rules(&refusal).contains(&MachineBehaviorRule::DuplicateCondition));
}

#[test]
fn g0_conditions_check_type_clock_frame_causality_and_history_events() {
    let graph = valid_graph().admit().expect("base graph admits");

    let mut quantity_gap = valid_behavior();
    quantity_gap.conditions[0].quantity =
        TerminalQuantitySpec::Dimensional(Dims([1, 0, 0, 0, 0, 0]));
    let refusal = quantity_gap
        .admit_against(&graph)
        .expect_err("condition quantity gap refuses");
    assert!(rules(&refusal).contains(&MachineBehaviorRule::ConditionQuantityGap));

    let mut shape_gap = valid_behavior();
    shape_gap.conditions[0].shape = TerminalShape::Vector { components: nz(2) };
    let refusal = shape_gap
        .admit_against(&graph)
        .expect_err("condition shape gap refuses");
    assert!(rules(&refusal).contains(&MachineBehaviorRule::ConditionShapeGap));

    let mut clock_gap = valid_behavior();
    clock_gap.conditions[0].clock = ClockId::new("clock/events").expect("valid clock id");
    let refusal = clock_gap
        .admit_against(&graph)
        .expect_err("condition clock gap refuses");
    assert!(rules(&refusal).contains(&MachineBehaviorRule::ConditionClockGap));

    let mut frame_gap = valid_behavior();
    frame_gap.conditions[0].frame =
        FrameBinding::new("world/offset", OrientationParity::Preserving).expect("valid frame");
    let refusal = frame_gap
        .admit_against(&graph)
        .expect_err("condition frame gap refuses");
    assert!(rules(&refusal).contains(&MachineBehaviorRule::ConditionFrameGap));

    let mut initial_history = valid_behavior();
    initial_history.conditions[0].source = ConditionSource::History {
        history: history(50),
        continuity: HistoryContinuity::Continuous,
    };
    let refusal = initial_history
        .admit_against(&graph)
        .expect_err("state history is not an initial value");
    assert!(rules(&refusal).contains(&MachineBehaviorRule::InvalidInitialSource));

    let mut unknown_reset_event = valid_behavior();
    let ConditionSource::History {
        continuity: HistoryContinuity::ResetAtEvents { events },
        ..
    } = &mut unknown_reset_event.conditions[2].source
    else {
        panic!("fixture boundary is a reset-delimited history");
    };
    events.push(EventId::new("event/missing").expect("valid event id"));
    let refusal = unknown_reset_event
        .admit_against(&graph)
        .expect_err("unknown history event refuses");
    assert!(rules(&refusal).contains(&MachineBehaviorRule::UnknownHistoryEvent));

    let mut wrong_target = valid_behavior();
    wrong_target.conditions[2].target = ConditionTarget::Boundary(
        TerminalId::new("terminal/guard-observation").expect("valid terminal id"),
    );
    let refusal = wrong_target
        .admit_against(&graph)
        .expect_err("output terminal cannot receive boundary source");
    assert!(rules(&refusal).contains(&MachineBehaviorRule::BoundaryCausalityGap));
}

#[test]
fn g0_motion_event_and_reset_semantics_are_structurally_closed() {
    let graph = valid_graph().admit().expect("base graph admits");

    let mut missing_motion = valid_behavior();
    missing_motion.motions.clear();
    let refusal = missing_motion
        .admit_against(&graph)
        .expect_err("omitted static motion refuses");
    assert!(rules(&refusal).contains(&MachineBehaviorRule::MissingMotion));

    let mut wrong_clock = valid_behavior();
    wrong_clock.events[0].clock = ClockId::new("clock/continuous").expect("valid clock id");
    let refusal = wrong_clock
        .admit_against(&graph)
        .expect_err("event needs event-driven clock");
    assert!(rules(&refusal).contains(&MachineBehaviorRule::NonEventClock));

    let mut duplicate_dependency = valid_behavior();
    let repeated_dependency = duplicate_dependency.events[0].dependencies[0].clone();
    duplicate_dependency.events[0]
        .dependencies
        .push(repeated_dependency);
    let refusal = duplicate_dependency
        .admit_against(&graph)
        .expect_err("duplicate guard dependency refuses");
    assert!(rules(&refusal).contains(&MachineBehaviorRule::DuplicateEventDependency));

    let mut empty_dependencies = valid_behavior();
    empty_dependencies.events[0].dependencies.clear();
    let refusal = empty_dependencies
        .admit_against(&graph)
        .expect_err("empty guard dependency set refuses");
    assert!(rules(&refusal).contains(&MachineBehaviorRule::EmptyEventDependencies));

    let mut unknown_dependency = valid_behavior();
    unknown_dependency.events[0]
        .dependencies
        .push(EventDependency::State(
            StateSlotId::new("state/missing").expect("valid state id"),
        ));
    let refusal = unknown_dependency
        .admit_against(&graph)
        .expect_err("unknown guard dependency refuses");
    assert!(rules(&refusal).contains(&MachineBehaviorRule::UnknownEventDependency));

    let mut duplicate_write = valid_behavior();
    let ResetSemantics::Deterministic { writes, .. } = &mut duplicate_write.events[0].reset else {
        panic!("fixture has deterministic reset");
    };
    writes.push(writes[0].clone());
    let refusal = duplicate_write
        .admit_against(&graph)
        .expect_err("duplicate reset write refuses");
    assert!(rules(&refusal).contains(&MachineBehaviorRule::DuplicateResetWrite));

    let mut empty_writes = valid_behavior();
    let ResetSemantics::Deterministic { writes, .. } = &mut empty_writes.events[0].reset else {
        panic!("fixture has deterministic reset");
    };
    writes.clear();
    let refusal = empty_writes
        .admit_against(&graph)
        .expect_err("empty reset write set refuses");
    assert!(rules(&refusal).contains(&MachineBehaviorRule::EmptyResetWrites));

    let mut unknown_write = valid_behavior();
    let ResetSemantics::Deterministic { writes, .. } = &mut unknown_write.events[0].reset else {
        panic!("fixture has deterministic reset");
    };
    writes.push(StateSlotId::new("state/missing").expect("valid state id"));
    let refusal = unknown_write
        .admit_against(&graph)
        .expect_err("unknown reset write refuses");
    assert!(rules(&refusal).contains(&MachineBehaviorRule::UnknownResetState));
}

#[test]
fn g0_superdense_order_and_set_valued_resets_are_explicit() {
    let graph = valid_graph().admit().expect("base graph admits");
    let mut duplicate_priority = valid_behavior();
    let mut second = duplicate_priority.events[0].clone();
    second.id = EventId::new("event/contact-backup").expect("valid event id");
    second.guard = guard(60);
    second.reset = ResetSemantics::SetValued {
        relation: reset(61),
        outcomes: outcomes(62),
        writes: vec![StateSlotId::new("state/position").expect("valid state id")],
    };
    duplicate_priority.events.push(second);
    let refusal = duplicate_priority
        .clone()
        .admit_against(&graph)
        .expect_err("same-clock duplicate microstep refuses");
    assert!(rules(&refusal).contains(&MachineBehaviorRule::DuplicateEventPriority));

    let mut mixed = duplicate_priority.clone();
    let mixed_group = simultaneity(65);
    for event in &mut mixed.events {
        event.order = EventOrder::SetValued {
            group: mixed_group.clone(),
        };
    }
    let mut priority_event = mixed.events[0].clone();
    priority_event.id = EventId::new("event/priority-third").expect("valid event id");
    priority_event.guard = guard(66);
    priority_event.order = EventOrder::TotalPriority { microstep: 2 };
    mixed.events.push(priority_event);
    let refusal = mixed
        .admit_against(&graph)
        .expect_err("one clock cannot mix priority and set-valued collision policy");
    assert!(rules(&refusal).contains(&MachineBehaviorRule::MixedEventOrderPolicy));

    let group = simultaneity(63);
    for event in &mut duplicate_priority.events {
        event.order = EventOrder::SetValued {
            group: group.clone(),
        };
    }
    duplicate_priority
        .admit_against(&graph)
        .expect("two same-clock events may retain set-valued simultaneity");

    let mut multiple_groups = valid_behavior();
    for (key, byte) in [
        ("event/group-a-peer", 67),
        ("event/group-b-first", 68),
        ("event/group-b-peer", 69),
    ] {
        let mut event = multiple_groups.events[0].clone();
        event.id = EventId::new(key).expect("valid event id");
        event.guard = guard(byte);
        multiple_groups.events.push(event);
    }
    let group_a = simultaneity(70);
    let group_b = simultaneity(71);
    for (index, event) in multiple_groups.events.iter_mut().enumerate() {
        event.order = EventOrder::SetValued {
            group: if index < 2 {
                group_a.clone()
            } else {
                group_b.clone()
            },
        };
    }
    let refusal = multiple_groups
        .admit_against(&graph)
        .expect_err("one clock cannot name multiple collision groups");
    assert!(rules(&refusal).contains(&MachineBehaviorRule::MultipleSimultaneityGroups));

    let mut singleton = valid_behavior();
    singleton.events[0].order = EventOrder::SetValued {
        group: simultaneity(64),
    };
    let refusal = singleton
        .admit_against(&graph)
        .expect_err("a set-valued group must name multiple events");
    assert!(rules(&refusal).contains(&MachineBehaviorRule::SingletonSimultaneityGroup));
}

#[test]
#[allow(clippy::too_many_lines)] // One matrix keeps joint closure failures visible together.
fn g0_tolerances_and_dependence_never_infer_independence() {
    assert_eq!(scalar(-0.0), scalar(0.0));
    assert!(FiniteNonNegative::new(-f64::from_bits(1)).is_err());
    assert!(FiniteNonNegative::new(f64::NAN).is_err());

    let graph = valid_graph().admit().expect("base graph admits");
    let mut zero = valid_behavior();
    let ToleranceSemantics::Random { scale, .. } = &mut zero.tolerances[0].semantics else {
        panic!("fixture tolerance is random");
    };
    *scale = scalar(0.0);
    let refusal = zero
        .admit_against(&graph)
        .expect_err("zero random scale refuses");
    assert!(rules(&refusal).contains(&MachineBehaviorRule::ZeroTolerance));

    let mut unknown_target = valid_behavior();
    unknown_target.tolerances[0].target = ToleranceTarget::Element(
        StateSlotId::new("state/missing")
            .expect("valid state id")
            .into(),
    );
    let refusal = unknown_target
        .admit_against(&graph)
        .expect_err("unknown tolerance target refuses");
    assert!(rules(&refusal).contains(&MachineBehaviorRule::UnknownToleranceTarget));

    let mut tolerance_type_gap = valid_behavior();
    tolerance_type_gap.tolerances[0].quantity =
        TerminalQuantitySpec::Dimensional(Dims([1, 0, 0, 0, 0, 0]));
    tolerance_type_gap.tolerances[0].shape = TerminalShape::Vector { components: nz(2) };
    let refusal = tolerance_type_gap
        .admit_against(&graph)
        .expect_err("state tolerance contract mismatch refuses");
    let refusal_rules = rules(&refusal);
    assert!(refusal_rules.contains(&MachineBehaviorRule::ToleranceQuantityGap));
    assert!(refusal_rules.contains(&MachineBehaviorRule::ToleranceShapeGap));

    let mut duplicate_binding = valid_behavior();
    let mut second = duplicate_binding.tolerances[2].clone();
    second.id = ToleranceId::new("tolerance/body-clearance-backup").expect("valid tolerance id");
    duplicate_binding.tolerances.push(second);
    let refusal = duplicate_binding
        .admit_against(&graph)
        .expect_err("one exact parameter has one v1 tolerance law");
    assert!(rules(&refusal).contains(&MachineBehaviorRule::DuplicateToleranceBinding));

    let mut missing = valid_behavior();
    missing.dependences.clear();
    let refusal = missing
        .admit_against(&graph)
        .expect_err("random dependence cannot be inferred");
    assert!(rules(&refusal).contains(&MachineBehaviorRule::MissingDependence));

    let mut missing_condition = valid_behavior();
    missing_condition.dependences[0]
        .members
        .retain(|member| !matches!(member, DependenceMember::Condition(_)));
    let refusal = missing_condition
        .admit_against(&graph)
        .expect_err("distribution-valued condition needs joint dependence");
    assert!(rules(&refusal).contains(&MachineBehaviorRule::MissingDependence));

    let mut fixed_condition_member = valid_behavior();
    fixed_condition_member.dependences[0]
        .members
        .push(DependenceMember::Condition(ConditionTarget::Initial(
            StateSlotId::new("state/position").expect("valid state id"),
        )));
    let refusal = fixed_condition_member
        .admit_against(&graph)
        .expect_err("fixed condition is not a random axis");
    assert!(rules(&refusal).contains(&MachineBehaviorRule::NonRandomConditionInDependence));

    let mut unknown_condition_member = valid_behavior();
    unknown_condition_member.dependences[0]
        .members
        .push(DependenceMember::Condition(ConditionTarget::Initial(
            StateSlotId::new("state/missing").expect("valid state id"),
        )));
    let refusal = unknown_condition_member
        .admit_against(&graph)
        .expect_err("unknown condition target is not a random axis");
    assert!(rules(&refusal).contains(&MachineBehaviorRule::UnknownConditionMember));

    let mut bounded_member = valid_behavior();
    let bounded = bounded_member.tolerances[2].id.clone();
    bounded_member.dependences[0]
        .members
        .push(DependenceMember::Tolerance(bounded));
    let refusal = bounded_member
        .admit_against(&graph)
        .expect_err("bounded tolerance is not a random correlation member");
    assert!(rules(&refusal).contains(&MachineBehaviorRule::BoundedToleranceInDependence));

    let mut singleton = valid_behavior();
    singleton.dependences[0].members.truncate(1);
    let refusal = singleton
        .admit_against(&graph)
        .expect_err("correlation needs at least two members");
    assert!(rules(&refusal).contains(&MachineBehaviorRule::CorrelatedGroupTooSmall));

    let mut independent = valid_behavior();
    independent.dependences[0].model = DependenceModel::Independent;
    independent
        .admit_against(&graph)
        .expect("independence is accepted only when explicitly declared");

    let mut multiple = valid_behavior();
    let repeated_dependence = multiple.dependences[0].clone();
    multiple.dependences.push(repeated_dependence);
    let refusal = multiple
        .admit_against(&graph)
        .expect_err("v1 admits one global joint model only");
    assert!(rules(&refusal).contains(&MachineBehaviorRule::MultipleDependenceModels));
    assert!(rules(&refusal).contains(&MachineBehaviorRule::DuplicateDependenceCoverage));
}

#[test]
fn g3_selected_refs_scalars_motion_and_base_graph_move_identity() {
    let graph = valid_graph().admit().expect("base graph admits");
    let baseline = valid_behavior()
        .admit_against(&graph)
        .expect("baseline behavior admits")
        .identity();
    let assert_moves = |behavior: MachineBehaviorDraft, label: &str| {
        assert_ne!(
            behavior
                .admit_against(&graph)
                .unwrap_or_else(|error| panic!("{label} must remain admissible: {error}"))
                .identity(),
            baseline,
            "{label} must move behavior identity"
        );
    };

    let mut changed_guard = valid_behavior();
    changed_guard.events[0].guard = guard(99);
    assert_moves(changed_guard, "guard reference");

    let mut changed_condition = valid_behavior();
    changed_condition.conditions[0].source = ConditionSource::Fixed(condition_value(98));
    assert_moves(changed_condition, "condition value reference");

    let mut changed_crossing = valid_behavior();
    changed_crossing.events[0].crossing = CrossingSemantics::Unknown(no_claim(97));
    assert_moves(changed_crossing, "crossing no-claim reference");

    let mut changed_reset = valid_behavior();
    let ResetSemantics::Deterministic { map, writes } = &mut changed_reset.events[0].reset else {
        panic!("fixture reset is deterministic");
    };
    *map = reset(96);
    writes.pop();
    assert_moves(changed_reset, "reset map and write set");

    let mut changed_order = valid_behavior();
    changed_order.events[0].order = EventOrder::TotalPriority { microstep: 1 };
    assert_moves(changed_order, "event microstep");

    let mut prescribed = valid_behavior();
    prescribed.motions[0].motion = BodyMotion::Prescribed {
        path: motion_path(95),
    };
    assert_moves(prescribed, "motion path");

    let mut changed_tolerance = valid_behavior();
    let ToleranceSemantics::Bounded { plus, .. } = &mut changed_tolerance.tolerances[2].semantics
    else {
        panic!("fixture tolerance is bounded");
    };
    *plus = scalar(f64::from_bits(0.02f64.to_bits() + 1));
    assert_moves(changed_tolerance, "adjacent tolerance width");

    let mut changed_parameter = valid_behavior();
    changed_parameter.tolerances[0].parameter = parameter("parameters/position-alt", 94);
    let ToleranceSemantics::Random { law, marginal, .. } =
        &mut changed_parameter.tolerances[0].semantics
    else {
        panic!("fixture tolerance is random");
    };
    *law = tolerance_law(93);
    *marginal = distribution(92);
    assert_moves(changed_parameter, "parameter law and marginal references");

    let mut changed_correlation = valid_behavior();
    changed_correlation.dependences[0].model = DependenceModel::Correlated(correlation(91));
    assert_moves(changed_correlation, "correlation model reference");

    let mut changed_dependence = valid_behavior();
    changed_dependence.dependences[0].model = DependenceModel::Independent;
    assert_moves(changed_dependence, "dependence model kind");

    let mut changed_graph = valid_graph();
    changed_graph.subsystems[0].model = model(98);
    let changed_graph = changed_graph.admit().expect("changed graph admits");
    assert_ne!(
        valid_behavior()
            .admit_against(&changed_graph)
            .expect("same behavior binds changed graph")
            .identity(),
        baseline
    );
}

#[test]
fn g3_refusal_findings_are_permutation_invariant() {
    let graph = valid_graph().admit().expect("base graph admits");
    let mut first = valid_behavior();
    first.conditions.push(first.conditions[0].clone());
    let repeated_dependency = first.events[0].dependencies[0].clone();
    first.events[0].dependencies.push(repeated_dependency);
    let repeated_member = first.dependences[0].members[0].clone();
    first.dependences[0].members.push(repeated_member);
    let mut second = first.clone();
    second.state_contracts.reverse();
    second.conditions.reverse();
    second.events.reverse();
    second.events[0].dependencies.reverse();
    second.tolerances.reverse();
    second.dependences.reverse();
    second.dependences[0].members.reverse();
    let first = first
        .admit_against(&graph)
        .expect_err("fixture intentionally refuses");
    let second = second
        .admit_against(&graph)
        .expect_err("permuted fixture intentionally refuses");
    assert_eq!(first.findings(), second.findings());
}

#[test]
fn g0_resource_envelope_and_nominal_roles_fail_closed() {
    let graph = valid_graph().admit().expect("base graph admits");
    let mut behavior = valid_behavior();
    let contract = behavior.state_contracts[0].clone();
    behavior.state_contracts = vec![contract; MAX_MACHINE_BEHAVIOR_STATE_CONTRACTS + 1];
    let decision = behavior.admit_with_decision(&graph);
    assert_eq!(
        decision.submitted_counts().state_contracts,
        MAX_MACHINE_BEHAVIOR_STATE_CONTRACTS + 1
    );
    let refusal = decision.result().expect_err("resource limit refuses first");
    assert_eq!(refusal.findings().len(), 1);
    assert_eq!(
        refusal.findings()[0].rule(),
        MachineBehaviorRule::ResourceLimit
    );

    let mut nested = valid_behavior();
    let dependency = nested.events[0].dependencies[0].clone();
    nested.events[0].dependencies = vec![dependency; MAX_MACHINE_BEHAVIOR_NESTED_REFERENCES + 1];
    let decision = nested.admit_with_decision(&graph);
    assert!(decision.submitted_counts().nested_references > MAX_MACHINE_BEHAVIOR_NESTED_REFERENCES);
    let refusal = decision
        .result()
        .expect_err("aggregate nested-reference limit refuses first");
    assert_eq!(refusal.findings().len(), 1);
    assert_eq!(
        refusal.findings()[0].rule(),
        MachineBehaviorRule::ResourceLimit
    );

    assert!(GuardRef::new("Guards/Bad", nz(1), digest(1)).is_err());
    assert!(DistributionRef::new("distributions/a", nz(1), [0; 32]).is_err());
    let event = EventId::new("shared/key").expect("valid event id");
    let tolerance = ToleranceId::new("shared/key").expect("valid tolerance id");
    assert_ne!(event.identity().as_bytes(), tolerance.identity().as_bytes());
}
