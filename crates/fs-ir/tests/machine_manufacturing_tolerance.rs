//! Machine-IR/fs-toleralloc manufacturing-axis crosswalk (Gauntlet G0/G3/G5).

use core::num::NonZeroU64;

use std::collections::BTreeSet;

use fs_blake3::ContentHash;
use fs_evidence::ColorRank;
use fs_ir::machine::manufacturing::tolerance_axis::MachineToleranceAxisCrosswalkErrorV1;
use fs_ir::machine::manufacturing::{
    AdmittedMachineManufacturingStateV1, MachineManufacturingDraftV1, ManufacturingArtifactRefV1,
    ManufacturingProcessKindV1, ManufacturingProcessStepV1, ManufacturingStepIdV1,
};
use fs_ir::machine::semantics::{
    AdmittedMachineBehavior, BodyMotion, ConditionBinding, ConditionSource, ConditionTarget,
    ConditionValueRef, CorrelationModelRef, DependenceMember, DependenceModel, DependenceSpec,
    DistributionRef, FiniteNonNegative, MachineBehaviorDraft, MotionBinding, ParameterRef,
    ToleranceId, ToleranceLawRef, ToleranceSemantics, ToleranceSpec, ToleranceTarget,
};
use fs_ir::machine::{
    AdmittedMachineGraph, BodyId, ClockId, ClockSpec, DependentKind, FrameBinding, LineageEvent,
    LineageRecord, LineageRefusal, LineageRelation, MachineClock, MachineElementId,
    MachineGraphDraft, MaterialBinding, MaterialCardRef, MaterialTarget, ModelRef,
    OrientationParity, SubsystemId, SubsystemSpec, TerminalCausality, TerminalId,
    TerminalQuantitySpec, TerminalShape, TerminalSpec,
};
use fs_qty::Dims;
use fs_toleralloc::{
    AdmittedCorrelationModel, CorrelatedStackReceipt, CorrelatedStackTerm,
    CorrelationAdmissionError, propagate_correlated_stack,
};

const CORRELATION_NAMESPACE: &str = "correlations/manufacturing-stack";

fn nz(value: u64) -> NonZeroU64 {
    NonZeroU64::new(value).expect("test schema version is nonzero")
}

fn scalar(value: f64) -> FiniteNonNegative {
    FiniteNonNegative::new(value).expect("finite nonnegative scalar")
}

fn frame() -> FrameBinding {
    FrameBinding::new("frame/manufacturing", OrientationParity::Preserving).expect("valid frame")
}

fn quantity() -> TerminalQuantitySpec {
    TerminalQuantitySpec::Dimensional(Dims::NONE)
}

fn model_ref(byte: u8) -> CorrelationModelRef {
    CorrelationModelRef::new(CORRELATION_NAMESPACE, nz(1), [byte; 32])
        .expect("valid Machine correlation reference")
}

fn artifact(namespace: &str, byte: u8) -> ManufacturingArtifactRefV1 {
    ManufacturingArtifactRefV1::new(namespace, nz(1), ContentHash([byte; 32]))
        .expect("valid manufacturing artifact")
}

fn distribution(byte: u8) -> DistributionRef {
    DistributionRef::new("distributions/manufacturing-axis", nz(1), [byte; 32])
        .expect("valid distribution")
}

fn graph(model_byte: u8) -> AdmittedMachineGraph {
    let subsystem = SubsystemId::new("subsystem/gearbox").unwrap();
    let gear = BodyId::new("body/gear").unwrap();
    let fixture = BodyId::new("body/fixture").unwrap();
    let clock = ClockId::new("clock/continuous").unwrap();
    MachineGraphDraft {
        clocks: vec![ClockSpec {
            id: clock.clone(),
            clock: MachineClock::Continuous,
        }],
        subsystems: vec![SubsystemSpec {
            id: subsystem.clone(),
            model: ModelRef::new("models/gearbox", nz(1), [model_byte; 32]).unwrap(),
            bodies: vec![gear.clone(), fixture.clone()],
            surface_patches: Vec::new(),
            contact_features: Vec::new(),
            state_slots: Vec::new(),
        }],
        terminals: vec![TerminalSpec {
            id: TerminalId::new("terminal/ambient").unwrap(),
            owner: subsystem,
            quantity: quantity(),
            shape: TerminalShape::Scalar,
            causality: TerminalCausality::ExternalInput,
            clock,
            frame: frame(),
        }],
        ports: Vec::new(),
        relations: Vec::new(),
        materials: vec![
            MaterialBinding {
                target: MaterialTarget::Body(gear),
                material: MaterialCardRef::new("materials/gear", nz(1), [1; 32]).unwrap(),
            },
            MaterialBinding {
                target: MaterialTarget::Body(fixture),
                material: MaterialCardRef::new("materials/fixture", nz(1), [2; 32]).unwrap(),
            },
        ],
        interfaces: Vec::new(),
    }
    .admit()
    .expect("crosswalk fixture graph admits")
}

fn tolerance(index: usize) -> ToleranceSpec {
    ToleranceSpec {
        id: ToleranceId::new(format!("tolerance/axis-{index:03}")).unwrap(),
        target: ToleranceTarget::Element(MachineElementId::Body(BodyId::new("body/gear").unwrap())),
        parameter: ParameterRef::new(
            format!("parameters/manufacturing-axis-{index:03}"),
            nz(1),
            [30_u8.wrapping_add(index as u8); 32],
        )
        .unwrap(),
        quantity: quantity(),
        shape: TerminalShape::Scalar,
        semantics: ToleranceSemantics::Random {
            scale: scalar(0.01 + index as f64 * 0.001),
            law: ToleranceLawRef::new(
                "tolerances/manufacturing-additive",
                nz(1),
                [60_u8.wrapping_add(index as u8); 32],
            )
            .unwrap(),
            marginal: distribution(90_u8.wrapping_add(index as u8)),
        },
    }
}

fn behavior_draft(axis_count: usize, dependence: Option<DependenceModel>) -> MachineBehaviorDraft {
    let clock = ClockId::new("clock/continuous").unwrap();
    let tolerances: Vec<_> = (0..axis_count).map(tolerance).collect();
    let members = tolerances
        .iter()
        .map(|specification| DependenceMember::Tolerance(specification.id.clone()))
        .collect();
    MachineBehaviorDraft {
        state_contracts: Vec::new(),
        conditions: vec![ConditionBinding {
            target: ConditionTarget::Boundary(TerminalId::new("terminal/ambient").unwrap()),
            quantity: quantity(),
            shape: TerminalShape::Scalar,
            clock: clock.clone(),
            frame: frame(),
            source: ConditionSource::Fixed(
                ConditionValueRef::new("values/ambient", nz(1), [10; 32]).unwrap(),
            ),
        }],
        motions: vec![
            MotionBinding {
                body: BodyId::new("body/gear").unwrap(),
                clock: clock.clone(),
                reference_frame: frame(),
                motion: BodyMotion::Static,
            },
            MotionBinding {
                body: BodyId::new("body/fixture").unwrap(),
                clock,
                reference_frame: frame(),
                motion: BodyMotion::Static,
            },
        ],
        events: Vec::new(),
        tolerances,
        dependences: dependence
            .map(|model| vec![DependenceSpec { members, model }])
            .unwrap_or_default(),
    }
}

fn behavior(
    graph: &AdmittedMachineGraph,
    axis_count: usize,
    model_digest: u8,
) -> AdmittedMachineBehavior {
    behavior_draft(
        axis_count,
        Some(DependenceModel::Correlated(model_ref(model_digest))),
    )
    .admit_against(graph)
    .expect("tolerance-only behavior admits")
}

fn manufacturing(
    graph: &AdmittedMachineGraph,
    process_byte: u8,
) -> AdmittedMachineManufacturingStateV1 {
    manufacturing_from_steps(
        graph,
        vec![process_step(
            "process/gear/machine",
            "body/gear",
            process_byte,
        )],
    )
}

fn manufacturing_from_steps(
    graph: &AdmittedMachineGraph,
    process_steps: Vec<ManufacturingProcessStepV1>,
) -> AdmittedMachineManufacturingStateV1 {
    MachineManufacturingDraftV1 {
        // This content digest intentionally differs from the toleralloc
        // semantic digest. The crosswalk requires a separate link artifact.
        correlation_model: artifact("manufacturing/correlation-model-content", 200),
        process_steps,
    }
    .admit_against(graph)
    .expect("body manufacturing history admits")
}

fn process_step(id: &str, body: &str, byte: u8) -> ManufacturingProcessStepV1 {
    ManufacturingProcessStepV1::new(
        ManufacturingStepIdV1::new(id).unwrap(),
        BodyId::new(body).unwrap(),
        None,
        ManufacturingProcessKindV1::Machining,
        artifact("manufacturing/process-specification", byte),
        artifact("manufacturing/input-material-state", byte.wrapping_add(1)),
        artifact("manufacturing/microstructure-state", byte.wrapping_add(2)),
        artifact("manufacturing/residual-stress-state", byte.wrapping_add(3)),
        artifact("manufacturing/property-state", byte.wrapping_add(4)),
    )
}

fn identity_factor(dimension: usize) -> Vec<f64> {
    let mut factor = vec![0.0; dimension * dimension];
    for index in 0..dimension {
        factor[index * dimension + index] = 1.0;
    }
    factor
}

fn correlation_model(digest: u8, dimension: usize, factor: Vec<f64>) -> AdmittedCorrelationModel {
    AdmittedCorrelationModel::try_new(
        CORRELATION_NAMESPACE,
        nz(1),
        [digest; 32],
        dimension,
        factor,
    )
    .expect("valid correlation model")
}

fn terms_for_behavior(behavior: &AdmittedMachineBehavior) -> Vec<CorrelatedStackTerm> {
    behavior.dependences()[0]
        .members
        .iter()
        .enumerate()
        .map(|(index, member)| {
            let name = match member {
                DependenceMember::Tolerance(tolerance) => tolerance.canonical_key().to_owned(),
                DependenceMember::Condition(ConditionTarget::Initial(state)) => {
                    format!("condition/initial/{}", state.canonical_key())
                }
                DependenceMember::Condition(ConditionTarget::Boundary(terminal)) => {
                    format!("condition/boundary/{}", terminal.canonical_key())
                }
            };
            CorrelatedStackTerm {
                name,
                signed_sensitivity: index as f64 + 1.0,
                sensitivity_color: ColorRank::Validated,
                // Deliberately differs from the Machine random scale: the
                // crosswalk is structural and does not infer scale-to-sigma.
                standard_deviation: 0.25 + index as f64 * 0.125,
            }
        })
        .collect()
}

fn stack(
    model: &AdmittedCorrelationModel,
    terms: &[CorrelatedStackTerm],
) -> CorrelatedStackReceipt {
    propagate_correlated_stack(model, terms).expect("finite correlated stack propagates")
}

fn coordinate_link(byte: u8) -> ManufacturingArtifactRefV1 {
    artifact("manufacturing/correlation-coordinate-link", byte)
}

#[test]
#[allow(clippy::too_many_lines)] // One identity and lineage mutation matrix.
fn mmt_001_crosswalk_binds_exact_receipts_identity_and_body_lineage() {
    let graph = graph(1);
    let behavior = behavior(&graph, 2, 90);
    let manufacturing = manufacturing(&graph, 20);
    let model = correlation_model(90, 2, identity_factor(2));
    let terms = terms_for_behavior(&behavior);
    let stack = stack(&model, &terms);
    let link = coordinate_link(210);
    let admitted = manufacturing
        .bind_correlated_tolerance_axes(&behavior, &stack, link.clone())
        .expect("closed tolerance-axis crosswalk admits");

    assert_eq!(admitted.graph(), graph.identity());
    assert_eq!(admitted.behavior(), behavior.identity());
    assert_eq!(admitted.manufacturing_state(), manufacturing.identity());
    assert_eq!(admitted.correlation_coordinate_link(), &link);
    assert_eq!(admitted.stack_receipt(), &stack);
    assert_eq!(admitted.tolerance_axes().len(), 2);
    for axis in admitted.tolerance_axes() {
        assert_eq!(
            axis.tolerance().canonical_key(),
            terms[axis.position()].name
        );
        assert_eq!(axis.body().canonical_key(), "body/gear");
        let specification = behavior
            .tolerances()
            .iter()
            .find(|candidate| candidate.id == *axis.tolerance())
            .expect("bound tolerance remains present");
        let ToleranceSemantics::Random { scale, .. } = &specification.semantics else {
            panic!("dependence tolerance remains random");
        };
        assert_ne!(
            scale.get().to_bits(),
            terms[axis.position()].standard_deviation.to_bits(),
            "structural crosswalk must not require Machine scale == stack sigma"
        );
    }

    let replay = manufacturing
        .bind_correlated_tolerance_axes(&behavior, &stack, link.clone())
        .expect("deterministic replay admits");
    assert_eq!(admitted.identity(), replay.identity());
    assert_eq!(
        admitted.identity_receipt().canonical_preimage(),
        replay.identity_receipt().canonical_preimage()
    );

    let correlated_model = correlation_model(90, 2, vec![1.0, 0.0, 0.5, 0.75_f64.sqrt()]);
    let correlated_stack = crate::stack(&correlated_model, &terms);
    let factor_changed = manufacturing
        .bind_correlated_tolerance_axes(&behavior, &correlated_stack, link.clone())
        .expect("same digest with different exact factor remains structurally bindable");
    assert_ne!(admitted.identity(), factor_changed.identity());

    for mutation in 0..3 {
        let mut changed_terms = terms.clone();
        match mutation {
            0 => changed_terms[0].signed_sensitivity = 1.5,
            1 => changed_terms[0].sensitivity_color = ColorRank::Verified,
            2 => changed_terms[0].standard_deviation = 0.5,
            _ => unreachable!(),
        }
        let changed_stack = crate::stack(&model, &changed_terms);
        let changed = manufacturing
            .bind_correlated_tolerance_axes(&behavior, &changed_stack, link.clone())
            .expect("term mutation remains structurally admissible");
        assert_ne!(admitted.identity(), changed.identity());
    }

    let changed_link = manufacturing
        .bind_correlated_tolerance_axes(&behavior, &stack, coordinate_link(211))
        .expect("alternate explicit coordinate link admits");
    assert_ne!(admitted.identity(), changed_link.identity());

    let mut changed_behavior_draft =
        behavior_draft(2, Some(DependenceModel::Correlated(model_ref(90))));
    let ToleranceSemantics::Random { marginal, .. } =
        &mut changed_behavior_draft.tolerances[0].semantics
    else {
        unreachable!()
    };
    *marginal = distribution(222);
    let changed_behavior = changed_behavior_draft
        .admit_against(&graph)
        .expect("changed marginal admits");
    let changed_behavior_terms = terms_for_behavior(&changed_behavior);
    let changed_behavior_stack = crate::stack(&model, &changed_behavior_terms);
    let behavior_changed = manufacturing
        .bind_correlated_tolerance_axes(&changed_behavior, &changed_behavior_stack, link.clone())
        .expect("changed behavior remains crosswalkable");
    assert_ne!(admitted.identity(), behavior_changed.identity());

    let changed_manufacturing = crate::manufacturing(&graph, 21);
    let manufacturing_changed = changed_manufacturing
        .bind_correlated_tolerance_axes(&behavior, &stack, link)
        .expect("changed process receipt remains crosswalkable");
    assert_ne!(admitted.identity(), manufacturing_changed.identity());

    let source = BodyId::new("body/gear").unwrap();
    let successor = BodyId::new("body/gear-worn").unwrap();
    let exact_keys: BTreeSet<_> = admitted
        .tolerance_axes()
        .iter()
        .map(|axis| axis.tolerance().canonical_key().to_owned())
        .collect();
    let wear = LineageRecord::admit(
        LineageEvent::Wear,
        vec![
            LineageRelation::new(
                MachineElementId::Body(source.clone()),
                vec![MachineElementId::Body(successor.clone())],
            )
            .unwrap(),
        ],
        admitted.lineage_dependents(),
    )
    .expect("one-to-one body wear rebinds every tolerance axis");
    let rebound: BTreeSet<_> = wear
        .rebindings()
        .iter()
        .map(|binding| {
            assert_eq!(
                binding.dependent().kind(),
                DependentKind::ManufacturingToleranceAxis
            );
            (
                binding.dependent().canonical_key().to_owned(),
                binding.dependent().source().clone(),
                binding.target().clone(),
            )
        })
        .collect();
    let expected_rebound: BTreeSet<_> = exact_keys
        .iter()
        .cloned()
        .map(|key| {
            (
                key,
                MachineElementId::Body(source.clone()),
                MachineElementId::Body(successor.clone()),
            )
        })
        .collect();
    assert_eq!(rebound, expected_rebound);

    let split = LineageRecord::admit(
        LineageEvent::Split,
        vec![
            LineageRelation::new(
                MachineElementId::Body(source.clone()),
                vec![
                    MachineElementId::Body(BodyId::new("body/gear-left").unwrap()),
                    MachineElementId::Body(BodyId::new("body/gear-right").unwrap()),
                ],
            )
            .unwrap(),
        ],
        admitted.lineage_dependents(),
    )
    .expect_err("ambiguous body split invalidates every tolerance axis");
    let LineageRefusal::Ambiguous(invalidation) = split else {
        panic!("expected ambiguity invalidation")
    };
    let invalidated: BTreeSet<_> = invalidation
        .invalidated_dependents()
        .iter()
        .map(|binding| {
            assert_eq!(binding.kind(), DependentKind::ManufacturingToleranceAxis);
            (binding.canonical_key().to_owned(), binding.source().clone())
        })
        .collect();
    let expected_invalidated: BTreeSet<_> = exact_keys
        .into_iter()
        .map(|key| (key, MachineElementId::Body(source.clone())))
        .collect();
    assert_eq!(invalidated, expected_invalidated);

    println!(
        "{{\"suite\":\"fs-ir/machine-manufacturing-tolerance\",\"case\":\"mmt-001\",\
         \"verdict\":\"pass\",\"axes\":2,\"detail\":\"exact factor/terms/receipts and \
         body-axis attachments are structurally bound without scale-to-sigma authority\"}}"
    );
}

#[test]
#[allow(clippy::too_many_lines)] // One complete structured-refusal matrix.
fn mmt_002_semantic_gaps_refuse_before_identity_publication() {
    let graph = graph(1);
    let behavior = behavior(&graph, 2, 90);
    let manufacturing = manufacturing(&graph, 20);
    let model = correlation_model(90, 2, identity_factor(2));
    let terms = terms_for_behavior(&behavior);
    let stack = stack(&model, &terms);

    let other_manufacturing = crate::manufacturing(&crate::graph(2), 20);
    assert!(matches!(
        other_manufacturing.bind_correlated_tolerance_axes(&behavior, &stack, coordinate_link(210)),
        Err(MachineToleranceAxisCrosswalkErrorV1::GraphMismatch { .. })
    ));

    let no_dependence = behavior_draft(0, None)
        .admit_against(&graph)
        .expect("behavior without random axes needs no dependence");
    assert!(matches!(
        manufacturing.bind_correlated_tolerance_axes(&no_dependence, &stack, coordinate_link(210)),
        Err(MachineToleranceAxisCrosswalkErrorV1::MissingDependence)
    ));

    let independent = behavior_draft(2, Some(DependenceModel::Independent))
        .admit_against(&graph)
        .expect("explicit independent behavior admits");
    assert!(matches!(
        manufacturing.bind_correlated_tolerance_axes(&independent, &stack, coordinate_link(210)),
        Err(MachineToleranceAxisCrosswalkErrorV1::IndependentDependence)
    ));

    let model3 = correlation_model(90, 3, identity_factor(3));
    let terms3 = vec![
        terms[0].clone(),
        terms[1].clone(),
        CorrelatedStackTerm {
            name: "tolerance/extra-axis".to_owned(),
            signed_sensitivity: 3.0,
            sensitivity_color: ColorRank::Validated,
            standard_deviation: 0.5,
        },
    ];
    let stack3 = crate::stack(&model3, &terms3);
    assert!(matches!(
        manufacturing.bind_correlated_tolerance_axes(&behavior, &stack3, coordinate_link(210)),
        Err(
            MachineToleranceAxisCrosswalkErrorV1::AxisDimensionMismatch {
                behavior: 2,
                stack: 3
            }
        )
    ));

    let mut mixed_draft = behavior_draft(2, Some(DependenceModel::Correlated(model_ref(90))));
    mixed_draft.conditions[0].source = ConditionSource::Distribution(distribution(150));
    mixed_draft.dependences[0]
        .members
        .push(DependenceMember::Condition(ConditionTarget::Boundary(
            TerminalId::new("terminal/ambient").unwrap(),
        )));
    let mixed = mixed_draft
        .admit_against(&graph)
        .expect("mixed random condition/tolerance behavior admits structurally");
    let mixed_terms = terms_for_behavior(&mixed);
    let mixed_stack = crate::stack(&model3, &mixed_terms);
    assert!(matches!(
        manufacturing.bind_correlated_tolerance_axes(&mixed, &mixed_stack, coordinate_link(210)),
        Err(MachineToleranceAxisCrosswalkErrorV1::ConditionAxisUnsupported { .. })
    ));

    let wrong_model = correlation_model(91, 2, identity_factor(2));
    let wrong_model_stack = crate::stack(&wrong_model, &terms);
    assert!(matches!(
        manufacturing.bind_correlated_tolerance_axes(
            &behavior,
            &wrong_model_stack,
            coordinate_link(210)
        ),
        Err(MachineToleranceAxisCrosswalkErrorV1::BehaviorCorrelationModelMismatch)
    ));

    let mut non_scalar_draft = behavior_draft(2, Some(DependenceModel::Correlated(model_ref(90))));
    non_scalar_draft.tolerances[0].shape = TerminalShape::Vector { components: nz(2) };
    let non_scalar = non_scalar_draft
        .admit_against(&graph)
        .expect("body vector tolerance is structurally admitted by behavior v1");
    let non_scalar_terms = terms_for_behavior(&non_scalar);
    let non_scalar_stack = crate::stack(&model, &non_scalar_terms);
    assert!(matches!(
        manufacturing.bind_correlated_tolerance_axes(
            &non_scalar,
            &non_scalar_stack,
            coordinate_link(210)
        ),
        Err(MachineToleranceAxisCrosswalkErrorV1::NonScalarTolerance { .. })
    ));

    let mut subsystem_draft = behavior_draft(2, Some(DependenceModel::Correlated(model_ref(90))));
    subsystem_draft.tolerances[0].target =
        ToleranceTarget::Subsystem(SubsystemId::new("subsystem/gearbox").unwrap());
    let subsystem = subsystem_draft
        .admit_against(&graph)
        .expect("subsystem tolerance is valid behavior but not body manufacturing v1");
    let subsystem_terms = terms_for_behavior(&subsystem);
    let subsystem_stack = crate::stack(&model, &subsystem_terms);
    assert!(matches!(
        manufacturing.bind_correlated_tolerance_axes(
            &subsystem,
            &subsystem_stack,
            coordinate_link(210)
        ),
        Err(MachineToleranceAxisCrosswalkErrorV1::UnsupportedToleranceTarget { .. })
    ));

    let mut unmanufactured_draft =
        behavior_draft(2, Some(DependenceModel::Correlated(model_ref(90))));
    unmanufactured_draft.tolerances[0].target =
        ToleranceTarget::Element(MachineElementId::Body(BodyId::new("body/fixture").unwrap()));
    let unmanufactured = unmanufactured_draft
        .admit_against(&graph)
        .expect("fixture tolerance is valid behavior");
    let unmanufactured_terms = terms_for_behavior(&unmanufactured);
    let unmanufactured_stack = crate::stack(&model, &unmanufactured_terms);
    assert!(matches!(
        manufacturing.bind_correlated_tolerance_axes(
            &unmanufactured,
            &unmanufactured_stack,
            coordinate_link(210)
        ),
        Err(MachineToleranceAxisCrosswalkErrorV1::MissingManufacturingHistory { .. })
    ));

    let mut wrong_names = terms.clone();
    wrong_names[0].name = "tolerance/wrong-position".to_owned();
    let wrong_name_stack = crate::stack(&model, &wrong_names);
    assert!(matches!(
        manufacturing.bind_correlated_tolerance_axes(
            &behavior,
            &wrong_name_stack,
            coordinate_link(210)
        ),
        Err(MachineToleranceAxisCrosswalkErrorV1::AxisNameMismatch { .. })
    ));

    println!(
        "{{\"suite\":\"fs-ir/machine-manufacturing-tolerance\",\"case\":\"mmt-002\",\
         \"verdict\":\"pass\",\"detail\":\"graph/dependence/model/dimension/axis/lineage \
         gaps refuse without a crosswalk identity\"}}"
    );
}

#[test]
#[allow(clippy::too_many_lines)] // Both body contexts stay visible in one event.
fn mmt_003_partial_ambiguity_is_attachment_local_and_requires_readmission() {
    let graph = graph(1);
    let mut draft = behavior_draft(2, Some(DependenceModel::Correlated(model_ref(90))));
    draft.tolerances[0].target =
        ToleranceTarget::Element(MachineElementId::Body(BodyId::new("body/fixture").unwrap()));
    let behavior = draft
        .admit_against(&graph)
        .expect("two-body tolerance behavior admits");
    let manufacturing = manufacturing_from_steps(
        &graph,
        vec![
            process_step("process/gear/machine", "body/gear", 20),
            process_step("process/fixture/machine", "body/fixture", 30),
        ],
    );
    let model = correlation_model(90, 2, identity_factor(2));
    let terms = terms_for_behavior(&behavior);
    let stack = stack(&model, &terms);
    let admitted = manufacturing
        .bind_correlated_tolerance_axes(&behavior, &stack, coordinate_link(210))
        .expect("both body attachments have manufacturing history");
    assert_eq!(admitted.tolerance_axes().len(), 2);
    assert_eq!(
        admitted
            .tolerance_axes()
            .iter()
            .map(|axis| axis.body().canonical_key())
            .collect::<BTreeSet<_>>(),
        BTreeSet::from(["body/fixture", "body/gear"])
    );

    let gear_axis = admitted
        .tolerance_axes()
        .iter()
        .find(|axis| axis.body().canonical_key() == "body/gear")
        .expect("gear axis exists");
    let refusal = LineageRecord::admit(
        LineageEvent::Remesh,
        vec![
            LineageRelation::new(
                MachineElementId::Body(BodyId::new("body/gear").unwrap()),
                vec![
                    MachineElementId::Body(BodyId::new("body/gear-left").unwrap()),
                    MachineElementId::Body(BodyId::new("body/gear-right").unwrap()),
                ],
            )
            .unwrap(),
            LineageRelation::new(
                MachineElementId::Body(BodyId::new("body/fixture").unwrap()),
                vec![MachineElementId::Body(
                    BodyId::new("body/fixture-refurbished").unwrap(),
                )],
            )
            .unwrap(),
        ],
        admitted.lineage_dependents(),
    )
    .expect_err("one ambiguous body refuses the aggregate lineage event");
    let invalidation = refusal
        .invalidation()
        .expect("typed ambiguity invalidation");
    assert_eq!(invalidation.considered_dependents().len(), 2);
    assert_eq!(invalidation.ambiguous_relations().len(), 1);
    assert_eq!(invalidation.invalidated_dependents().len(), 1);
    assert_eq!(
        invalidation.invalidated_dependents()[0].canonical_key(),
        gear_axis.tolerance().canonical_key()
    );
    assert_eq!(
        admitted.behavior(),
        behavior.identity(),
        "the old receipt remains immutable and must not masquerade as a successor"
    );

    println!(
        "{{\"suite\":\"fs-ir/machine-manufacturing-tolerance\",\"case\":\"mmt-003\",\
         \"verdict\":\"pass\",\"considered\":2,\"invalidated\":1,\"detail\":\"partial \
         ambiguity is attachment-local and successor crosswalks require explicit readmission\"}}"
    );
}

#[test]
fn mmt_004_exact_128_axis_owner_boundary_is_preserved() {
    let graph = graph(1);
    let behavior = behavior(&graph, 128, 90);
    let manufacturing = manufacturing(&graph, 20);
    let model = correlation_model(90, 128, identity_factor(128));
    let terms = terms_for_behavior(&behavior);
    let stack = stack(&model, &terms);
    let admitted = manufacturing
        .bind_correlated_tolerance_axes(&behavior, &stack, coordinate_link(210))
        .expect("owner's exact maximum axis dimension remains crosswalkable");
    assert_eq!(admitted.tolerance_axes().len(), 128);
    assert_eq!(admitted.tolerance_axes()[0].position(), 0);
    assert_eq!(admitted.tolerance_axes()[127].position(), 127);

    assert!(matches!(
        AdmittedCorrelationModel::try_new(
            CORRELATION_NAMESPACE,
            nz(1),
            [90; 32],
            129,
            identity_factor(129),
        ),
        Err(CorrelationAdmissionError::InvalidDimension {
            dimension: 129,
            max: 128
        })
    ));

    println!(
        "{{\"suite\":\"fs-ir/machine-manufacturing-tolerance\",\"case\":\"mmt-004\",\
         \"verdict\":\"pass\",\"axes\":128,\"detail\":\"crosswalk inherits the exact \
         fs-toleralloc owner cap and cannot widen it\"}}"
    );
}
