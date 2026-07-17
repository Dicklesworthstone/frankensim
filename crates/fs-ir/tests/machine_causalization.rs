//! Equation-variable hypergraph and causalization receipt conformance.
//!
//! The names record the Gauntlet tier exercised by each test. Assertions keep
//! structural semantics, provenance, hybrid theorem commitments, cancellation,
//! and no-claim boundaries independently observable in failure logs.

use core::num::NonZeroU64;

use fs_blake3::identity::{
    CanonicalError, CanonicalSchema, LimitKind, Presence, StrongIdentity, WireType,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ir::machine::causalization::*;
use fs_ir::machine::semantics::{
    AdmittedMachineBehavior, ConditionBinding, ConditionSource, ConditionTarget, ConditionValueRef,
    MachineBehaviorDraft, StateSlotContract,
};
use fs_ir::machine::{
    AdmittedMachineGraph, ClockId, ClockSpec, FrameBinding, MAX_MACHINE_ENTITY_KEY_BYTES,
    MachineClock, MachineElementId, MachineGraphDraft, ModelRef, OrientationParity, RelationId,
    RelationMode, RelationSpec, StateSlotId, SubsystemId, SubsystemSpec, TerminalCausality,
    TerminalId, TerminalQuantitySpec, TerminalShape, TerminalSpec,
};
use fs_qty::Dims;

fn nz(value: u64) -> NonZeroU64 {
    NonZeroU64::new(value).expect("test fixture value is nonzero")
}

fn schema_field_count<S: CanonicalSchema>() -> u32 {
    u32::try_from(S::FIELDS.len()).expect("test schema field count fits u32")
}

fn schema_layout<S: CanonicalSchema>() -> Vec<(&'static str, WireType, Presence)> {
    S::FIELDS
        .iter()
        .map(|field| (field.name(), field.wire_type(), field.presence()))
        .collect()
}

macro_rules! cref {
    ($ty:ident, $namespace:literal, $byte:expr) => {
        $ty::new($namespace, nz(1), [$byte; 32]).expect("valid causal reference")
    };
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x5eed,
                kernel_id: 2,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn with_cancelled_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    gate.request();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x5eed,
                kernel_id: 2,
                tile: 0,
                iteration: 1,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

#[allow(clippy::too_many_arguments)]
fn incidence_spec(
    equation: EquationId,
    variable: VariableId,
    derivative_order: u16,
    solve_participation: SolveParticipation,
    coefficient_dimensions: Dims,
    term: SignalContract,
    operator: Option<IncidenceOperatorRef>,
    clock_relation: IncidenceClockRelation,
    activation: ActivationDomain,
) -> Result<IncidenceSpec, CanonicalError> {
    with_cx(|cx| {
        IncidenceSpec::new(
            equation,
            variable,
            derivative_order,
            solve_participation,
            coefficient_dimensions,
            term,
            operator,
            clock_relation,
            activation,
            cx,
        )
    })
}

fn maximum_matching_binding(
    graph: &AdmittedCausalGraph,
    domain: CausalReceiptDomain,
    matching: &[CausalMatchingPair],
    certificate: MaximumMatchingCertificateRef,
    checker: CausalCheckerRef,
) -> Result<MaximumMatchingBinding, MaximumMatchingBindingError> {
    with_cx(|cx| MaximumMatchingBinding::new(graph, domain, matching, certificate, checker, cx))
}

fn conditional_coverage_binding(
    graph: &AdmittedCausalGraph,
    outcomes: &[ConditionalCausalOutcome],
    certificate: ConditionalCoverageRef,
    checker: CausalCheckerRef,
) -> Result<ConditionalCoverageBinding, ConditionalCoverageBindingError> {
    with_cx(|cx| {
        ConditionalCoverageBinding::for_mode_cells(graph, outcomes, certificate, checker, cx)
    })
}

fn conditional_outcome(
    child: &AdmittedCausalizationReceipt,
) -> Result<ConditionalCausalOutcome, ConditionalOutcomeError> {
    with_cx(|cx| ConditionalCausalOutcome::from_mode_cell(child, cx))
}

fn minimal_machine() -> (AdmittedMachineGraph, SubsystemId, ClockId) {
    let owner = SubsystemId::new("subsystem/plant").expect("valid subsystem ID");
    let clock = ClockId::new("clock/continuous").expect("valid clock ID");
    let graph = MachineGraphDraft {
        clocks: vec![ClockSpec {
            id: clock.clone(),
            clock: MachineClock::Continuous,
        }],
        subsystems: vec![SubsystemSpec {
            id: owner.clone(),
            model: ModelRef::new("models/causal-test", nz(1), [1; 32]).expect("valid model ref"),
            bodies: Vec::new(),
            surface_patches: Vec::new(),
            contact_features: Vec::new(),
            state_slots: Vec::new(),
        }],
        terminals: Vec::new(),
        ports: Vec::new(),
        relations: Vec::new(),
        materials: Vec::new(),
        interfaces: Vec::new(),
    }
    .admit()
    .expect("minimal Machine graph admits");
    (graph, owner, clock)
}

fn state_machine(model_byte: u8) -> (AdmittedMachineGraph, SubsystemId, ClockId, StateSlotId) {
    let owner = SubsystemId::new("subsystem/state-plant").expect("valid subsystem ID");
    let clock = ClockId::new("clock/state-continuous").expect("valid clock ID");
    let state = StateSlotId::new("state/position").expect("valid state ID");
    let source = TerminalSpec {
        id: TerminalId::new("terminal/state-source").expect("valid state source terminal"),
        owner: owner.clone(),
        quantity: TerminalQuantitySpec::Dimensional(Dims::NONE),
        shape: TerminalShape::Scalar,
        causality: TerminalCausality::Output,
        clock: clock.clone(),
        frame: frame(),
    };
    let sink = TerminalSpec {
        id: TerminalId::new("terminal/state-sink").expect("valid state sink terminal"),
        owner: owner.clone(),
        quantity: TerminalQuantitySpec::Dimensional(Dims::NONE),
        shape: TerminalShape::Scalar,
        causality: TerminalCausality::Input,
        clock: clock.clone(),
        frame: frame(),
    };
    let graph = MachineGraphDraft {
        clocks: vec![ClockSpec {
            id: clock.clone(),
            clock: MachineClock::Continuous,
        }],
        subsystems: vec![SubsystemSpec {
            id: owner.clone(),
            model: ModelRef::new("models/causal-state-test", nz(1), [model_byte; 32])
                .expect("valid state model ref"),
            bodies: Vec::new(),
            surface_patches: Vec::new(),
            contact_features: Vec::new(),
            state_slots: vec![state.clone()],
        }],
        terminals: vec![source.clone(), sink.clone()],
        ports: Vec::new(),
        relations: vec![RelationSpec {
            id: RelationId::new("relation/state-position").expect("valid state relation"),
            source: source.id,
            target: sink.id,
            mode: RelationMode::Stateful {
                state_slot: state.clone(),
            },
        }],
        materials: Vec::new(),
        interfaces: Vec::new(),
    }
    .admit()
    .expect("minimal state Machine graph admits");
    (graph, owner, clock, state)
}

fn state_behavior(
    machine: &AdmittedMachineGraph,
    owner: &SubsystemId,
    clock: &ClockId,
    state: &StateSlotId,
    initial_value_byte: u8,
    shape: TerminalShape,
) -> AdmittedMachineBehavior {
    MachineBehaviorDraft {
        state_contracts: vec![StateSlotContract {
            id: state.clone(),
            owner: owner.clone(),
            quantity: TerminalQuantitySpec::Dimensional(Dims::NONE),
            shape,
            clock: clock.clone(),
            frame: frame(),
        }],
        conditions: vec![ConditionBinding {
            target: ConditionTarget::Initial(state.clone()),
            quantity: TerminalQuantitySpec::Dimensional(Dims::NONE),
            shape,
            clock: clock.clone(),
            frame: frame(),
            source: ConditionSource::Fixed(
                ConditionValueRef::new("test/state-initial-value", nz(1), [initial_value_byte; 32])
                    .expect("valid state initial-value ref"),
            ),
        }],
        motions: Vec::new(),
        events: Vec::new(),
        tolerances: Vec::new(),
        dependences: Vec::new(),
    }
    .admit_against(machine)
    .expect("minimal state behavior admits")
}

fn frame() -> FrameBinding {
    FrameBinding::new("world/mechanical", OrientationParity::Preserving)
        .expect("valid frame binding")
}

fn signal(clock: &ClockId) -> SignalContract {
    SignalContract {
        quantity: TerminalQuantitySpec::Dimensional(Dims::NONE),
        shape: TerminalShape::Scalar,
        clock: clock.clone(),
        frame: frame(),
    }
}

fn lineage(owner: &SubsystemId, namespace: &'static str, byte: u8) -> NodeLineage {
    NodeLineage::new(
        NodeOrigin::Machine(MachineNodeOrigin::Subsystem(owner.clone())),
        CausalOwner::Subsystem(owner.clone()),
        NormalizedNodeSemanticRef::new(namespace, nz(1), [byte; 32])
            .expect("valid normalized node meaning"),
    )
}

fn extraction(byte: u8) -> CausalExtractionContext {
    CausalExtractionContext {
        extractor: CausalExtractorRef::new("test/causal-extractor", nz(1), [byte; 32])
            .expect("extractor ref"),
        coverage: cref!(CausalExtractionCoverageRef, "test/coverage", 11),
        evidence: CausalExtractionEvidence::Unverified,
        budget: cref!(CausalBudgetRef, "test/extraction-budget", 12),
        capabilities: cref!(CausalCapabilityRef, "test/extraction-capabilities", 13),
        seed_policy: CausalSeedPolicy::NoRandomness,
        determinism: CausalDeterminism::Deterministic,
    }
}

fn analysis() -> CausalAnalysisContext {
    CausalAnalysisContext {
        analyzer: cref!(CausalAnalyzerRef, "test/causal-analyzer", 20),
        budget: cref!(CausalBudgetRef, "test/analysis-budget", 21),
        capabilities: cref!(CausalCapabilityRef, "test/analysis-capabilities", 22),
        seed_policy: CausalSeedPolicy::NoRandomness,
        determinism: CausalDeterminism::Deterministic,
    }
}

fn equation(
    owner: &SubsystemId,
    clock: &ClockId,
    namespace: &'static str,
    byte: u8,
) -> EquationSpec {
    let lineage = lineage(owner, namespace, byte);
    EquationSpec {
        id: EquationId::derive(&lineage).expect("equation identity"),
        diagnostic_label: format!("equation-{byte}").into_boxed_str(),
        lineage,
        owner: CausalOwner::Subsystem(owner.clone()),
        supports: vec![CausalSupport::Lumped],
        residual: signal(clock),
        role: EquationRole::Constraint,
        solve_participation: EquationParticipation::Matching,
        activation: ActivationDomain::Always,
    }
}

fn variable(
    owner: &SubsystemId,
    clock: &ClockId,
    namespace: &'static str,
    byte: u8,
) -> VariableSpec {
    let lineage = lineage(owner, namespace, byte);
    VariableSpec {
        id: VariableId::derive(&lineage).expect("variable identity"),
        diagnostic_label: format!("variable-{byte}").into_boxed_str(),
        lineage,
        owner: CausalOwner::Subsystem(owner.clone()),
        supports: vec![CausalSupport::Lumped],
        value: signal(clock),
        role: VariableRole::Algebraic,
        solve_participation: SolveParticipation::Unknown,
        port_schema_crosswalk: None,
        activation: ActivationDomain::Always,
    }
}

fn state_variable(owner: &SubsystemId, clock: &ClockId, state: &StateSlotId) -> VariableSpec {
    let lineage = NodeLineage::new(
        NodeOrigin::Machine(MachineNodeOrigin::Element(MachineElementId::StateSlot(
            state.clone(),
        ))),
        CausalOwner::Subsystem(owner.clone()),
        NormalizedNodeSemanticRef::new("test/variable/state-position", nz(1), [29; 32])
            .expect("valid state-variable meaning"),
    );
    VariableSpec {
        id: VariableId::derive(&lineage).expect("state variable identity"),
        diagnostic_label: "state-position".into(),
        lineage,
        owner: CausalOwner::Subsystem(owner.clone()),
        supports: vec![CausalSupport::MachineElement(MachineElementId::StateSlot(
            state.clone(),
        ))],
        value: signal(clock),
        role: VariableRole::State,
        solve_participation: SolveParticipation::Unknown,
        port_schema_crosswalk: None,
        activation: ActivationDomain::Always,
    }
}

fn state_causal_draft(
    owner: &SubsystemId,
    clock: &ClockId,
    state: &StateSlotId,
) -> CausalGraphDraft {
    let mut equation = equation(owner, clock, "test/equation/state-position", 28);
    equation.role = EquationRole::StateUpdate;
    let variable = state_variable(owner, clock, state);
    let incidence = incidence(&equation, &variable);
    CausalGraphDraft {
        units: CausalUnitConvention::SiBaseDimensions,
        scope: CausalGraphScope::CompleteMachineModel,
        extraction: extraction(10),
        equations: vec![equation],
        variables: vec![variable],
        conditions: Vec::new(),
        incidences: vec![incidence],
    }
}

fn incidence(equation: &EquationSpec, variable: &VariableSpec) -> IncidenceSpec {
    incidence_spec(
        equation.id.clone(),
        variable.id.clone(),
        0,
        SolveParticipation::Unknown,
        Dims::NONE,
        equation.residual.clone(),
        None,
        IncidenceClockRelation::SameClock,
        ActivationDomain::Always,
    )
    .expect("incidence identity")
}

fn minimal_causal_draft(owner: &SubsystemId, clock: &ClockId) -> CausalGraphDraft {
    let equation = equation(owner, clock, "test/equation/constraint", 30);
    let variable = variable(owner, clock, "test/variable/x", 31);
    let incidence = incidence(&equation, &variable);
    CausalGraphDraft {
        units: CausalUnitConvention::SiBaseDimensions,
        scope: CausalGraphScope::CompleteMachineModel,
        extraction: extraction(10),
        equations: vec![equation],
        variables: vec![variable],
        conditions: Vec::new(),
        incidences: vec![incidence],
    }
}

fn admit_minimal() -> (AdmittedMachineGraph, AdmittedCausalGraph) {
    let (machine, owner, clock) = minimal_machine();
    let graph = with_cx(|cx| minimal_causal_draft(&owner, &clock).admit_against(&machine, cx))
        .expect("minimal causal graph admits");
    (machine, graph)
}

fn complete_receipt(graph: &AdmittedCausalGraph) -> CausalizationReceiptDraft {
    let equation = &graph.equations()[0];
    let variable = &graph.variables()[0];
    let incidence = &graph.incidences()[0];
    CausalizationReceiptDraft {
        structure: graph.structure_identity_receipt(),
        artifact: graph.artifact_identity_receipt(),
        analysis: analysis(),
        domain: CausalReceiptDomain::UnconditionalGraph,
        determination: DeterminationClass::WellDetermined,
        structural_rank: StructuralRankState::FullRelativeToMinSide,
        conditionality: Conditionality::Unconditional,
        matching: vec![CausalMatchingPair {
            incidence: incidence.id.clone(),
            equation: equation.id.clone(),
            variable: DerivativeVariableKey {
                variable: variable.id.clone(),
                derivative_order: 0,
            },
        }],
        unmatched_equations: Vec::new(),
        unmatched_variables: Vec::new(),
        conditional_outcomes: Vec::new(),
        maximum_matching_certificate: None,
        conditional_coverage: None,
        unknown_axes: Vec::new(),
        evidence: CausalReceiptEvidence::Unverified,
    }
}

fn graph_rules(refusal: &CausalGraphRefusal) -> Vec<CausalGraphRule> {
    refusal
        .findings()
        .iter()
        .map(CausalGraphFinding::rule)
        .collect()
}

fn assert_graph_rule(refusal: &CausalGraphRefusal, expected: CausalGraphRule) {
    let rules = graph_rules(refusal);
    assert!(
        rules.contains(&expected),
        "expected rule={} ({expected:?}); actual_rules={rules:?}; findings={:#?}",
        expected.code(),
        refusal.findings()
    );
}

fn assert_graph_rules_exact(refusal: &CausalGraphRefusal, expected: &[CausalGraphRule]) {
    let rules = graph_rules(refusal);
    assert_eq!(
        rules.as_slice(),
        expected,
        "unexpected graph-rule set; actual_rules={rules:?}; findings={:#?}",
        refusal.findings()
    );
}

fn receipt_rules(refusal: &CausalReceiptRefusal) -> Vec<CausalReceiptRule> {
    refusal
        .findings()
        .iter()
        .map(CausalReceiptFinding::rule)
        .collect()
}

fn assert_receipt_rule(refusal: &CausalReceiptRefusal, expected: CausalReceiptRule) {
    let rules = receipt_rules(refusal);
    assert!(
        rules.contains(&expected),
        "expected rule={} ({expected:?}); actual_rules={rules:?}; findings={:#?}",
        expected.code(),
        refusal.findings()
    );
}

fn assert_receipt_finding(
    refusal: &CausalReceiptRefusal,
    rule: CausalReceiptRule,
    subject: &CausalReceiptSubject,
) {
    assert!(
        refusal
            .findings()
            .iter()
            .any(|finding| finding.rule() == rule && finding.subject() == subject),
        "expected finding rule={} ({rule:?}), subject={subject:?}; findings={:#?}",
        rule.code(),
        refusal.findings()
    );
}

fn assert_receipt_rules_exact(refusal: &CausalReceiptRefusal, expected: &[CausalReceiptRule]) {
    let rules = receipt_rules(refusal);
    assert_eq!(
        rules.as_slice(),
        expected,
        "unexpected receipt-rule set; actual_rules={rules:?}; findings={:#?}",
        refusal.findings()
    );
}

#[test]
fn g0_minimal_graph_and_receipt_admit_with_complete_identity_receipts() {
    let (_machine, graph) = admit_minimal();
    assert_eq!(
        CausalStructureIdentitySchemaV1::FIELDS[1].wire_type(),
        WireType::Child
    );
    assert_eq!(
        CausalStructureIdentitySchemaV1::FIELDS[2].wire_type(),
        WireType::Bytes
    );
    assert_eq!(
        CausalStructureIdentitySchemaV1::FIELDS[8].wire_type(),
        WireType::OrderedChildren
    );
    assert_eq!(
        CausalStructureIdentitySchemaV1::FIELDS[9].wire_type(),
        WireType::OrderedBytes
    );
    assert_eq!(
        graph.structure_identity_receipt().field_count(),
        schema_field_count::<CausalStructureIdentitySchemaV1>()
    );
    assert_eq!(
        graph.artifact_identity_receipt().field_count(),
        schema_field_count::<CausalGraphArtifactIdentitySchemaV1>()
    );
    assert_eq!(graph.equations()[0].id.identity_receipt().field_count(), 1);
    assert_eq!(graph.variables()[0].id.identity_receipt().field_count(), 1);
    assert_eq!(graph.incidences()[0].id.identity_receipt().field_count(), 1);

    let decision = with_receipt_decision(complete_receipt(&graph), &graph);
    assert_eq!(decision.code(), "CausalReceiptAdmitted");
    assert!(decision.submitted_counts().complete);
    assert_eq!(decision.submitted_counts().matching, 1);
    let receipt = decision.into_result().expect("closed receipt admits");
    assert_eq!(
        CausalOutcomeIdentitySchemaV1::FIELDS[1].wire_type(),
        WireType::Child
    );
    assert_eq!(
        CausalOutcomeIdentitySchemaV1::FIELDS[2].wire_type(),
        WireType::Bytes
    );
    assert_eq!(
        CausalizationReceiptIdentitySchemaV1::FIELDS[1].wire_type(),
        WireType::Child
    );
    assert_eq!(
        CausalizationReceiptIdentitySchemaV1::FIELDS[2].wire_type(),
        WireType::Bytes
    );
    assert_eq!(
        CausalizationReceiptIdentitySchemaV1::FIELDS[3].wire_type(),
        WireType::Child
    );
    assert_eq!(
        CausalizationReceiptIdentitySchemaV1::FIELDS[4].wire_type(),
        WireType::Bytes
    );
    assert_eq!(
        CausalizationReceiptIdentitySchemaV1::FIELDS[16].wire_type(),
        WireType::Child
    );
    assert_eq!(
        CausalizationReceiptIdentitySchemaV1::FIELDS[17].wire_type(),
        WireType::Bytes
    );
    assert_eq!(
        receipt.outcome_identity_receipt().field_count(),
        schema_field_count::<CausalOutcomeIdentitySchemaV1>()
    );
    assert_eq!(
        receipt.identity_receipt().field_count(),
        schema_field_count::<CausalizationReceiptIdentitySchemaV1>()
    );
    assert_eq!(receipt.domain(), &CausalReceiptDomain::UnconditionalGraph);
    assert_eq!(receipt.matching().len(), 1);
}

#[test]
#[allow(clippy::too_many_lines)] // Exact candidate top-level layouts remain reviewable together.
fn g0_causal_identity_family_candidate_top_level_layout_is_explicit() {
    assert_eq!(CAUSAL_STRUCTURE_IDENTITY_SCHEMA_VERSION_V1, 1);
    assert_eq!(CAUSAL_GRAPH_ARTIFACT_IDENTITY_SCHEMA_VERSION_V1, 1);
    assert_eq!(CAUSAL_OUTCOME_IDENTITY_SCHEMA_VERSION_V1, 1);
    assert_eq!(CAUSALIZATION_RECEIPT_IDENTITY_SCHEMA_VERSION_V1, 1);
    assert_eq!(
        CausalStructureIdentitySchemaV1::VERSION,
        CAUSAL_STRUCTURE_IDENTITY_SCHEMA_VERSION_V1
    );
    assert_eq!(
        CausalGraphArtifactIdentitySchemaV1::VERSION,
        CAUSAL_GRAPH_ARTIFACT_IDENTITY_SCHEMA_VERSION_V1
    );
    assert_eq!(
        CausalOutcomeIdentitySchemaV1::VERSION,
        CAUSAL_OUTCOME_IDENTITY_SCHEMA_VERSION_V1
    );
    assert_eq!(
        CausalizationReceiptIdentitySchemaV1::VERSION,
        CAUSALIZATION_RECEIPT_IDENTITY_SCHEMA_VERSION_V1
    );
    assert_eq!(
        CausalStructureIdentitySchemaV1::DOMAIN,
        CAUSAL_STRUCTURE_IDENTITY_DOMAIN_V1
    );
    assert_eq!(
        CausalGraphArtifactIdentitySchemaV1::DOMAIN,
        CAUSAL_GRAPH_ARTIFACT_IDENTITY_DOMAIN_V1
    );
    assert_eq!(
        CausalOutcomeIdentitySchemaV1::DOMAIN,
        CAUSAL_OUTCOME_IDENTITY_DOMAIN_V1
    );
    assert_eq!(
        CausalizationReceiptIdentitySchemaV1::DOMAIN,
        CAUSALIZATION_RECEIPT_IDENTITY_DOMAIN_V1
    );
    assert_eq!(
        schema_layout::<CausalStructureIdentitySchemaV1>(),
        vec![
            (
                "causal-structure-schema-version",
                WireType::U64,
                Presence::Required
            ),
            ("machine-graph-id", WireType::Child, Presence::Required),
            (
                "machine-graph-receipt-adjudication",
                WireType::Bytes,
                Presence::Required
            ),
            ("unit-convention", WireType::Variant, Presence::Required),
            ("extraction-scope", WireType::Bytes, Presence::Required),
            ("equations", WireType::OrderedBytes, Presence::Required),
            ("variables", WireType::OrderedBytes, Presence::Required),
            (
                "activation-conditions",
                WireType::OrderedBytes,
                Presence::Required
            ),
            ("incidences", WireType::OrderedChildren, Presence::Required),
            (
                "incidence-receipt-adjudications",
                WireType::OrderedBytes,
                Presence::Required
            ),
        ]
    );
    assert_eq!(
        schema_layout::<CausalGraphArtifactIdentitySchemaV1>(),
        vec![
            (
                "causal-graph-artifact-schema-version",
                WireType::U64,
                Presence::Required
            ),
            ("causal-structure-id", WireType::Child, Presence::Required),
            (
                "causal-structure-receipt-adjudication",
                WireType::Bytes,
                Presence::Required
            ),
            (
                "machine-behavior-receipt-adjudication",
                WireType::Bytes,
                Presence::Optional
            ),
            ("extraction-context", WireType::Bytes, Presence::Required),
            (
                "equation-lineage",
                WireType::OrderedBytes,
                Presence::Required
            ),
            (
                "variable-lineage",
                WireType::OrderedBytes,
                Presence::Required
            ),
            (
                "activation-condition-provenance",
                WireType::OrderedBytes,
                Presence::Required
            ),
            (
                "incidence-provenance",
                WireType::OrderedBytes,
                Presence::Required
            ),
        ]
    );
    assert_eq!(
        schema_layout::<CausalOutcomeIdentitySchemaV1>(),
        vec![
            (
                "causal-outcome-schema-version",
                WireType::U64,
                Presence::Required
            ),
            ("causal-structure-id", WireType::Child, Presence::Required),
            (
                "causal-structure-receipt-adjudication",
                WireType::Bytes,
                Presence::Required
            ),
            ("analysis-domain", WireType::Bytes, Presence::Required),
            ("outcome-axes", WireType::Bytes, Presence::Required),
            ("matching", WireType::OrderedBytes, Presence::Required),
            (
                "unmatched-equations",
                WireType::OrderedBytes,
                Presence::Required
            ),
            (
                "unmatched-variables",
                WireType::OrderedBytes,
                Presence::Required
            ),
            (
                "conditional-outcome-semantics",
                WireType::OrderedBytes,
                Presence::Required
            ),
        ]
    );
    assert_eq!(
        schema_layout::<CausalizationReceiptIdentitySchemaV1>(),
        vec![
            (
                "causalization-receipt-schema-version",
                WireType::U64,
                Presence::Required
            ),
            ("causal-structure-id", WireType::Child, Presence::Required),
            (
                "causal-structure-receipt-adjudication",
                WireType::Bytes,
                Presence::Required
            ),
            (
                "causal-graph-artifact-id",
                WireType::Child,
                Presence::Required
            ),
            (
                "causal-graph-artifact-receipt-adjudication",
                WireType::Bytes,
                Presence::Required
            ),
            ("analysis-context", WireType::Bytes, Presence::Required),
            ("analysis-domain", WireType::Bytes, Presence::Required),
            ("outcome-axes", WireType::Bytes, Presence::Required),
            ("matching", WireType::OrderedBytes, Presence::Required),
            (
                "unmatched-equations",
                WireType::OrderedBytes,
                Presence::Required
            ),
            (
                "unmatched-variables",
                WireType::OrderedBytes,
                Presence::Required
            ),
            (
                "conditional-outcomes",
                WireType::OrderedBytes,
                Presence::Required
            ),
            (
                "maximum-matching-certificate",
                WireType::Bytes,
                Presence::Optional
            ),
            ("conditional-coverage", WireType::Bytes, Presence::Optional),
            ("unknown-axes", WireType::OrderedBytes, Presence::Required),
            ("evidence-state", WireType::Variant, Presence::Required),
            (
                "normalized-causal-outcome-id",
                WireType::Child,
                Presence::Required
            ),
            (
                "normalized-causal-outcome-receipt-adjudication",
                WireType::Bytes,
                Presence::Required
            ),
        ]
    );
}

#[test]
#[allow(clippy::too_many_lines)] // One representative cross-family quotient mutation matrix.
fn g0_causal_identity_family_semantic_mutations_move_receipts() {
    let (machine, owner, clock) = minimal_machine();
    let baseline_graph =
        with_cx(|cx| minimal_causal_draft(&owner, &clock).admit_against(&machine, cx))
            .expect("baseline causal graph admits");
    let baseline_receipt =
        with_cx(|cx| complete_receipt(&baseline_graph).admit_against(&baseline_graph, cx))
            .expect("baseline causal receipt admits");

    let mut scoped_draft = minimal_causal_draft(&owner, &clock);
    scoped_draft.scope = CausalGraphScope::Partial {
        boundary: cref!(CausalPartialScopeRef, "test/identity-scope", 134),
    };
    let scoped_graph =
        with_cx(|cx| scoped_draft.admit_against(&machine, cx)).expect("partial-scope graph admits");
    assert_ne!(
        baseline_graph.structure_identity(),
        scoped_graph.structure_identity()
    );
    assert_ne!(
        baseline_graph
            .structure_identity_receipt()
            .canonical_preimage(),
        scoped_graph
            .structure_identity_receipt()
            .canonical_preimage()
    );
    let scoped_receipt =
        with_cx(|cx| complete_receipt(&scoped_graph).admit_against(&scoped_graph, cx))
            .expect("partial-scope causal receipt admits");
    assert_ne!(
        baseline_receipt.outcome_identity(),
        scoped_receipt.outcome_identity(),
        "normalized structure is part of normalized outcome semantics"
    );
    assert_ne!(
        baseline_receipt
            .outcome_identity_receipt()
            .canonical_preimage(),
        scoped_receipt
            .outcome_identity_receipt()
            .canonical_preimage()
    );
    assert_ne!(baseline_receipt.identity(), scoped_receipt.identity());
    assert_ne!(
        baseline_receipt.identity_receipt().canonical_preimage(),
        scoped_receipt.identity_receipt().canonical_preimage()
    );

    let mut provenance_draft = minimal_causal_draft(&owner, &clock);
    provenance_draft.extraction = extraction(135);
    let provenance_graph = with_cx(|cx| provenance_draft.admit_against(&machine, cx))
        .expect("provenance-substituted graph admits");
    assert_eq!(
        baseline_graph.structure_identity(),
        provenance_graph.structure_identity(),
        "producer context must remain outside normalized structure"
    );
    assert_eq!(
        baseline_graph
            .structure_identity_receipt()
            .canonical_preimage(),
        provenance_graph
            .structure_identity_receipt()
            .canonical_preimage(),
        "producer context must preserve normalized structure bytes"
    );
    assert_ne!(
        baseline_graph.artifact_identity(),
        provenance_graph.artifact_identity(),
        "producer context must move the graph artifact"
    );
    assert_ne!(
        baseline_graph
            .artifact_identity_receipt()
            .canonical_preimage(),
        provenance_graph
            .artifact_identity_receipt()
            .canonical_preimage(),
        "producer context must move graph-artifact bytes"
    );
    let provenance_receipt =
        with_cx(|cx| complete_receipt(&provenance_graph).admit_against(&provenance_graph, cx))
            .expect("provenance-substituted causal receipt admits");
    assert_eq!(
        baseline_receipt.outcome_identity(),
        provenance_receipt.outcome_identity(),
        "graph-artifact provenance must remain outside normalized outcome semantics"
    );
    assert_eq!(
        baseline_receipt
            .outcome_identity_receipt()
            .canonical_preimage(),
        provenance_receipt
            .outcome_identity_receipt()
            .canonical_preimage(),
        "equal normalized outcomes must have byte-identical canonical preimages"
    );
    assert_ne!(
        baseline_receipt.identity(),
        provenance_receipt.identity(),
        "graph-artifact provenance must move the complete evidence receipt"
    );
    assert_ne!(
        baseline_receipt.identity_receipt().canonical_preimage(),
        provenance_receipt.identity_receipt().canonical_preimage()
    );

    let mut analysis_variant = complete_receipt(&baseline_graph);
    analysis_variant.analysis.analyzer =
        cref!(CausalAnalyzerRef, "test/identity-analysis-variant", 136);
    let analysis_receipt = with_cx(|cx| analysis_variant.admit_against(&baseline_graph, cx))
        .expect("analysis-context variant admits");
    assert_eq!(
        baseline_receipt.outcome_identity(),
        analysis_receipt.outcome_identity(),
        "producer analysis context must remain outside normalized outcome"
    );
    assert_eq!(
        baseline_receipt
            .outcome_identity_receipt()
            .canonical_preimage(),
        analysis_receipt
            .outcome_identity_receipt()
            .canonical_preimage()
    );
    assert_ne!(
        baseline_receipt.identity(),
        analysis_receipt.identity(),
        "producer analysis context must move the full receipt"
    );
    assert_ne!(
        baseline_receipt.identity_receipt().canonical_preimage(),
        analysis_receipt.identity_receipt().canonical_preimage()
    );

    let mut evidence_variant = complete_receipt(&baseline_graph);
    evidence_variant.evidence = CausalReceiptEvidence::CheckerReferenced(cref!(
        CausalCheckerRef,
        "test/identity-evidence-variant",
        137
    ));
    let evidence_receipt = with_cx(|cx| evidence_variant.admit_against(&baseline_graph, cx))
        .expect("evidence-state variant admits");
    assert_eq!(
        baseline_receipt.outcome_identity(),
        evidence_receipt.outcome_identity(),
        "checker provenance must remain outside normalized outcome semantics"
    );
    assert_eq!(
        baseline_receipt
            .outcome_identity_receipt()
            .canonical_preimage(),
        evidence_receipt
            .outcome_identity_receipt()
            .canonical_preimage()
    );
    assert_ne!(
        baseline_receipt.identity(),
        evidence_receipt.identity(),
        "checker provenance must move the complete evidence receipt"
    );
    assert_ne!(
        baseline_receipt.identity_receipt().canonical_preimage(),
        evidence_receipt.identity_receipt().canonical_preimage()
    );
}

#[test]
fn g0_state_aware_admission_binds_exact_behavior_contract_and_provenance() {
    let (machine, owner, clock, state) = state_machine(1);
    let behavior_a = state_behavior(&machine, &owner, &clock, &state, 2, TerminalShape::Scalar);
    let behavior_b = state_behavior(&machine, &owner, &clock, &state, 3, TerminalShape::Scalar);
    let draft = state_causal_draft(&owner, &clock, &state);

    let graph_a = with_cx(|cx| {
        draft
            .clone()
            .admit_against_behavior(&machine, &behavior_a, cx)
    })
    .expect("state graph with its exact behavior contract admits");
    let graph_b = with_cx(|cx| {
        draft
            .clone()
            .admit_against_behavior(&machine, &behavior_b, cx)
    })
    .expect("a second conformant behavior provenance admits");
    assert_eq!(graph_a.machine_behavior(), Some(behavior_a.identity()));
    assert_eq!(graph_b.machine_behavior(), Some(behavior_b.identity()));
    assert_eq!(
        graph_a.structure_identity(),
        graph_b.structure_identity(),
        "initial-value provenance must not change normalized causal structure"
    );
    assert_ne!(
        graph_a.artifact_identity(),
        graph_b.artifact_identity(),
        "distinct behavior provenance must move the causal graph artifact"
    );

    let refusal = with_cx(|cx| draft.clone().admit_against(&machine, cx))
        .expect_err("a state graph cannot omit its behavior overlay");
    assert_graph_rules_exact(&refusal, &[CausalGraphRule::StateBehaviorMismatch]);

    let incompatible_behavior = state_behavior(
        &machine,
        &owner,
        &clock,
        &state,
        2,
        TerminalShape::Vector { components: nz(2) },
    );
    let refusal = with_cx(|cx| {
        draft
            .clone()
            .admit_against_behavior(&machine, &incompatible_behavior, cx)
    })
    .expect_err("a state graph cannot substitute a shape-incompatible state contract");
    assert_graph_rules_exact(&refusal, &[CausalGraphRule::StateBehaviorMismatch]);

    let (foreign_machine, foreign_owner, foreign_clock, foreign_state) = state_machine(4);
    let foreign_behavior = state_behavior(
        &foreign_machine,
        &foreign_owner,
        &foreign_clock,
        &foreign_state,
        2,
        TerminalShape::Scalar,
    );
    let refusal = with_cx(|cx| draft.admit_against_behavior(&machine, &foreign_behavior, cx))
        .expect_err("a behavior overlay bound to another Machine graph must refuse");
    assert_graph_rules_exact(&refusal, &[CausalGraphRule::StateBehaviorMismatch]);
}

// A helper cannot return `Cx` because it borrows its arena. Tests that need a
// decision use this closure-shaped equivalent instead.
fn with_receipt_decision(
    draft: CausalizationReceiptDraft,
    graph: &AdmittedCausalGraph,
) -> CausalReceiptAdmissionDecision {
    with_cx(|cx| draft.admit_with_decision(graph, cx))
}

#[test]
fn g5_collection_permutations_and_labels_do_not_move_identity() {
    let (machine, owner, clock) = minimal_machine();
    let mut draft = minimal_causal_draft(&owner, &clock);
    let second_equation = equation(&owner, &clock, "test/equation/second", 32);
    let second_variable = variable(&owner, &clock, "test/variable/y", 33);
    let second_incidence = incidence(&second_equation, &second_variable);
    draft.equations.push(second_equation);
    draft.variables.push(second_variable);
    draft.incidences.push(second_incidence);
    let replay = draft.clone();
    let first = with_cx(|cx| draft.admit_against(&machine, cx)).expect("first graph");

    let mut permuted = replay;
    permuted.equations.reverse();
    permuted.variables.reverse();
    permuted.incidences.reverse();
    for equation in &mut permuted.equations {
        equation.diagnostic_label = "presentation-only rename".into();
    }
    for variable in &mut permuted.variables {
        variable.diagnostic_label = "another presentation-only rename".into();
    }
    let second = with_cx(|cx| permuted.admit_against(&machine, cx)).expect("replay graph");
    assert_eq!(first.structure_identity(), second.structure_identity());
    assert_eq!(first.artifact_identity(), second.artifact_identity());
}

#[test]
fn g3_normalized_structure_is_separate_from_producer_provenance() {
    let (machine, owner, clock) = minimal_machine();
    let first = with_cx(|cx| minimal_causal_draft(&owner, &clock).admit_against(&machine, cx))
        .expect("first graph");
    let mut changed = minimal_causal_draft(&owner, &clock);
    changed.extraction = extraction(99);
    let second = with_cx(|cx| changed.admit_against(&machine, cx)).expect("second graph");
    assert_eq!(first.structure_identity(), second.structure_identity());
    assert_ne!(first.artifact_identity(), second.artifact_identity());
}

#[test]
fn g3_audited_predicate_provenance_is_separate_from_normalized_structure() {
    let (machine, owner, clock) = minimal_machine();
    let first = with_cx(|cx| hybrid_draft(&owner, &clock).admit_against(&machine, cx))
        .expect("first audited-predicate graph");

    let mut changed_source = hybrid_draft(&owner, &clock);
    let alternate_source = cref!(SourceArtifactRef, "test/mode-predicate/alternate", 58);
    let ActivationConditionSource::AuditedPredicate(source_escape) =
        &mut changed_source.conditions[0].source
    else {
        panic!("hybrid fixture must use an audited predicate");
    };
    source_escape.source = alternate_source.clone();
    source_escape.audited_source = alternate_source;
    let source_variant = with_cx(|cx| changed_source.admit_against(&machine, cx))
        .expect("source-substituted audited-predicate graph");

    let mut changed_audit = hybrid_draft(&owner, &clock);
    let ActivationConditionSource::AuditedPredicate(audit_escape) =
        &mut changed_audit.conditions[0].source
    else {
        panic!("hybrid fixture must use an audited predicate");
    };
    audit_escape.audit = cref!(
        EscapeHatchAuditRef,
        "test/mode-predicate-audit/alternate",
        59
    );
    let audit_variant = with_cx(|cx| changed_audit.admit_against(&machine, cx))
        .expect("audit-substituted audited-predicate graph");

    assert_eq!(
        first.structure_identity(),
        source_variant.structure_identity(),
        "producer-specific opaque predicate artifacts must not move normalized structure"
    );
    assert_eq!(
        first.structure_identity_receipt().canonical_preimage(),
        source_variant
            .structure_identity_receipt()
            .canonical_preimage(),
        "source substitution must preserve the complete normalized preimage, not merely its digest"
    );
    assert_eq!(
        first.structure_identity(),
        audit_variant.structure_identity(),
        "predicate audit provenance must not move normalized structure"
    );
    assert_eq!(
        first.structure_identity_receipt().canonical_preimage(),
        audit_variant
            .structure_identity_receipt()
            .canonical_preimage(),
        "audit substitution must preserve the complete normalized preimage, not merely its digest"
    );
    assert_ne!(
        first.artifact_identity(),
        source_variant.artifact_identity(),
        "the provenance artifact must bind the exact predicate source"
    );
    assert_ne!(
        first.artifact_identity_receipt().canonical_preimage(),
        source_variant
            .artifact_identity_receipt()
            .canonical_preimage(),
        "source substitution must move the complete provenance preimage"
    );
    assert_ne!(
        first.artifact_identity(),
        audit_variant.artifact_identity(),
        "the provenance artifact must bind the exact predicate audit"
    );
    assert_ne!(
        first.artifact_identity_receipt().canonical_preimage(),
        audit_variant
            .artifact_identity_receipt()
            .canonical_preimage(),
        "audit substitution must move the complete provenance preimage"
    );

    let mut mismatched = hybrid_draft(&owner, &clock);
    let ActivationConditionSource::AuditedPredicate(mismatched_escape) =
        &mut mismatched.conditions[0].source
    else {
        panic!("hybrid fixture must use an audited predicate");
    };
    mismatched_escape.audited_source = cref!(
        SourceArtifactRef,
        "test/mode-predicate/mismatched-audit-target",
        60
    );
    let refusal = with_cx(|cx| mismatched.admit_against(&machine, cx))
        .expect_err("a declared audit target that differs from its source must refuse");
    assert_graph_rules_exact(&refusal, &[CausalGraphRule::EscapeAuditMismatch]);
}

#[test]
fn g0_matching_binds_the_exact_incidence_and_derivative_vertex() {
    let (_machine, graph) = admit_minimal();
    let mut wrong_order = complete_receipt(&graph);
    wrong_order.matching[0].variable.derivative_order = 1;
    wrong_order.unmatched_equations = vec![graph.equations()[0].id.clone()];
    wrong_order.unmatched_variables = vec![DerivativeVariableKey {
        variable: graph.variables()[0].id.clone(),
        derivative_order: 0,
    }];
    wrong_order.determination = DeterminationClass::Mixed;
    wrong_order.structural_rank = StructuralRankState::Unknown;
    wrong_order.unknown_axes = vec![CausalUnknownAxisState {
        axis: CausalOutcomeAxis::StructuralRank,
        reason: CausalUnknownReason::IncompleteMetadata,
        resume_checkpoint: None,
    }];
    let refusal = with_cx(|cx| wrong_order.clone().admit_against(&graph, cx))
        .expect_err("wrong derivative endpoint refuses");
    assert_receipt_rule(&refusal, CausalReceiptRule::UnknownMatchingEndpoint);
    assert_receipt_rule(&refusal, CausalReceiptRule::NonIncidenceMatch);

    wrong_order.unmatched_equations.clear();
    wrong_order.unmatched_variables.clear();
    wrong_order.determination = DeterminationClass::Unknown;
    wrong_order.unknown_axes.push(CausalUnknownAxisState {
        axis: CausalOutcomeAxis::Determination,
        reason: CausalUnknownReason::IncompleteMetadata,
        resume_checkpoint: None,
    });
    let refusal = with_cx(|cx| wrong_order.admit_against(&graph, cx))
        .expect_err("invalid matching cannot suppress either complement diagnostic");
    assert_receipt_finding(
        &refusal,
        CausalReceiptRule::UnmatchedSetMismatch,
        &CausalReceiptSubject::UnmatchedEquations,
    );
    assert_receipt_finding(
        &refusal,
        CausalReceiptRule::UnmatchedSetMismatch,
        &CausalReceiptSubject::UnmatchedVariables,
    );
}

#[test]
fn g0_duplicate_variable_match_does_not_consume_the_other_equation() {
    let (machine, owner, clock) = minimal_machine();
    let mut graph_draft = minimal_causal_draft(&owner, &clock);
    let second_equation = equation(&owner, &clock, "test/equation/duplicate-variable", 137);
    let second_incidence = incidence(&second_equation, &graph_draft.variables[0]);
    graph_draft.equations.push(second_equation);
    graph_draft.incidences.push(second_incidence);
    let graph = with_cx(|cx| graph_draft.admit_against(&machine, cx))
        .expect("two equations sharing one unknown variable admit structurally");

    let mut matching = graph
        .incidences()
        .iter()
        .map(|incidence| CausalMatchingPair {
            incidence: incidence.id.clone(),
            equation: incidence.equation.clone(),
            variable: DerivativeVariableKey {
                variable: incidence.variable.clone(),
                derivative_order: incidence.derivative_order,
            },
        })
        .collect::<Vec<_>>();
    matching.sort();
    let unmatched_equations = vec![matching[1].equation.clone()];
    let refusal = with_cx(|cx| {
        CausalizationReceiptDraft {
            structure: graph.structure_identity_receipt(),
            artifact: graph.artifact_identity_receipt(),
            analysis: analysis(),
            domain: CausalReceiptDomain::UnconditionalGraph,
            determination: DeterminationClass::OverDetermined,
            structural_rank: StructuralRankState::FullRelativeToMinSide,
            conditionality: Conditionality::Unconditional,
            matching,
            unmatched_equations,
            unmatched_variables: Vec::new(),
            conditional_outcomes: Vec::new(),
            maximum_matching_certificate: None,
            conditional_coverage: None,
            unknown_axes: Vec::new(),
            evidence: CausalReceiptEvidence::Unverified,
        }
        .admit_against(&graph, cx)
    })
    .expect_err("one unknown variable cannot be matched to two equations");
    assert_receipt_rule(&refusal, CausalReceiptRule::DuplicateMatchingEndpoint);
    assert!(
        !refusal
            .findings()
            .iter()
            .any(|finding| finding.rule() == CausalReceiptRule::UnmatchedSetMismatch),
        "the rejected second edge must leave its equation in the exact unmatched complement: {:#?}",
        refusal.findings()
    );
}

#[test]
fn g0_duplicate_equation_match_does_not_consume_the_other_variable() {
    let (machine, owner, clock) = minimal_machine();
    let mut graph_draft = minimal_causal_draft(&owner, &clock);
    let second_variable = variable(&owner, &clock, "test/variable/duplicate-equation", 142);
    let second_incidence = incidence(&graph_draft.equations[0], &second_variable);
    graph_draft.variables.push(second_variable);
    graph_draft.incidences.push(second_incidence);
    let graph = with_cx(|cx| graph_draft.admit_against(&machine, cx))
        .expect("one equation incident to two unknown variables admits structurally");

    let mut matching = graph
        .incidences()
        .iter()
        .map(|incidence| CausalMatchingPair {
            incidence: incidence.id.clone(),
            equation: incidence.equation.clone(),
            variable: DerivativeVariableKey {
                variable: incidence.variable.clone(),
                derivative_order: incidence.derivative_order,
            },
        })
        .collect::<Vec<_>>();
    matching.sort();
    let unmatched_variables = vec![matching[1].variable.clone()];
    let refusal = with_cx(|cx| {
        CausalizationReceiptDraft {
            structure: graph.structure_identity_receipt(),
            artifact: graph.artifact_identity_receipt(),
            analysis: analysis(),
            domain: CausalReceiptDomain::UnconditionalGraph,
            determination: DeterminationClass::UnderDetermined,
            structural_rank: StructuralRankState::FullRelativeToMinSide,
            conditionality: Conditionality::Unconditional,
            matching,
            unmatched_equations: Vec::new(),
            unmatched_variables,
            conditional_outcomes: Vec::new(),
            maximum_matching_certificate: None,
            conditional_coverage: None,
            unknown_axes: Vec::new(),
            evidence: CausalReceiptEvidence::Unverified,
        }
        .admit_against(&graph, cx)
    })
    .expect_err("one equation cannot be matched to two unknown variables");
    assert_receipt_rule(&refusal, CausalReceiptRule::DuplicateMatchingEndpoint);
    assert!(
        !refusal
            .findings()
            .iter()
            .any(|finding| finding.rule() == CausalReceiptRule::UnmatchedSetMismatch),
        "the rejected second edge must leave its variable in the exact unmatched complement: {:#?}",
        refusal.findings()
    );
}

#[test]
fn g0_empty_graph_refuses_vacuous_well_or_full_claims() {
    let (machine, _owner, _clock) = minimal_machine();
    let graph = with_cx(|cx| {
        CausalGraphDraft {
            units: CausalUnitConvention::SiBaseDimensions,
            scope: CausalGraphScope::CompleteMachineModel,
            extraction: extraction(10),
            equations: Vec::new(),
            variables: Vec::new(),
            conditions: Vec::new(),
            incidences: Vec::new(),
        }
        .admit_against(&machine, cx)
    })
    .expect("empty structural graph is representable");
    let false_claim = CausalizationReceiptDraft {
        structure: graph.structure_identity_receipt(),
        artifact: graph.artifact_identity_receipt(),
        analysis: analysis(),
        domain: CausalReceiptDomain::UnconditionalGraph,
        determination: DeterminationClass::WellDetermined,
        structural_rank: StructuralRankState::FullRelativeToMinSide,
        conditionality: Conditionality::Unconditional,
        matching: Vec::new(),
        unmatched_equations: Vec::new(),
        unmatched_variables: Vec::new(),
        conditional_outcomes: Vec::new(),
        maximum_matching_certificate: None,
        conditional_coverage: None,
        unknown_axes: Vec::new(),
        evidence: CausalReceiptEvidence::Unverified,
    };
    let refusal = with_cx(|cx| false_claim.clone().admit_against(&graph, cx))
        .expect_err("empty graph cannot mint vacuous Well/Full authority");
    assert_receipt_rule(&refusal, CausalReceiptRule::OutcomeAxisMismatch);

    let mut concrete_empty = false_claim;
    concrete_empty.determination = DeterminationClass::EmptyProjection;
    concrete_empty.structural_rank = StructuralRankState::NotApplicable;
    let admitted = with_cx(|cx| concrete_empty.admit_against(&graph, cx))
        .expect("a structurally empty projection is a concrete non-rank result");
    assert_eq!(
        admitted.determination(),
        DeterminationClass::EmptyProjection
    );
    assert_eq!(
        admitted.structural_rank(),
        StructuralRankState::NotApplicable
    );
}

#[test]
fn g0_one_sided_bipartitions_use_directed_determination_not_empty_projection() {
    let (machine, owner, clock) = minimal_machine();
    let isolated_equation = equation(&owner, &clock, "test/equation/isolated", 140);
    let equation_graph = with_cx(|cx| {
        CausalGraphDraft {
            units: CausalUnitConvention::SiBaseDimensions,
            scope: CausalGraphScope::CompleteMachineModel,
            extraction: extraction(10),
            equations: vec![isolated_equation.clone()],
            variables: Vec::new(),
            conditions: Vec::new(),
            incidences: Vec::new(),
        }
        .admit_against(&machine, cx)
    })
    .expect("an equation-only structural projection is representable");
    let equation_receipt = CausalizationReceiptDraft {
        structure: equation_graph.structure_identity_receipt(),
        artifact: equation_graph.artifact_identity_receipt(),
        analysis: analysis(),
        domain: CausalReceiptDomain::UnconditionalGraph,
        determination: DeterminationClass::OverDetermined,
        structural_rank: StructuralRankState::NotApplicable,
        conditionality: Conditionality::Unconditional,
        matching: Vec::new(),
        unmatched_equations: vec![isolated_equation.id],
        unmatched_variables: Vec::new(),
        conditional_outcomes: Vec::new(),
        maximum_matching_certificate: None,
        conditional_coverage: None,
        unknown_axes: Vec::new(),
        evidence: CausalReceiptEvidence::Unverified,
    };
    with_cx(|cx| equation_receipt.clone().admit_against(&equation_graph, cx))
        .expect("equation-only projection is concretely overdetermined");
    let mut false_equation_empty = equation_receipt;
    false_equation_empty.determination = DeterminationClass::EmptyProjection;
    let refusal = with_cx(|cx| false_equation_empty.admit_against(&equation_graph, cx))
        .expect_err("an equation-only projection is not empty");
    assert_receipt_rule(&refusal, CausalReceiptRule::OutcomeAxisMismatch);

    let isolated_variable = variable(&owner, &clock, "test/variable/isolated", 141);
    let variable_graph = with_cx(|cx| {
        CausalGraphDraft {
            units: CausalUnitConvention::SiBaseDimensions,
            scope: CausalGraphScope::CompleteMachineModel,
            extraction: extraction(10),
            equations: Vec::new(),
            variables: vec![isolated_variable.clone()],
            conditions: Vec::new(),
            incidences: Vec::new(),
        }
        .admit_against(&machine, cx)
    })
    .expect("an isolated order-zero unknown is a structural vertex");
    let variable_receipt = CausalizationReceiptDraft {
        structure: variable_graph.structure_identity_receipt(),
        artifact: variable_graph.artifact_identity_receipt(),
        analysis: analysis(),
        domain: CausalReceiptDomain::UnconditionalGraph,
        determination: DeterminationClass::UnderDetermined,
        structural_rank: StructuralRankState::NotApplicable,
        conditionality: Conditionality::Unconditional,
        matching: Vec::new(),
        unmatched_equations: Vec::new(),
        unmatched_variables: vec![DerivativeVariableKey {
            variable: isolated_variable.id,
            derivative_order: 0,
        }],
        conditional_outcomes: Vec::new(),
        maximum_matching_certificate: None,
        conditional_coverage: None,
        unknown_axes: Vec::new(),
        evidence: CausalReceiptEvidence::Unverified,
    };
    with_cx(|cx| variable_receipt.clone().admit_against(&variable_graph, cx))
        .expect("variable-only projection is concretely underdetermined");
    let mut false_variable_empty = variable_receipt;
    false_variable_empty.determination = DeterminationClass::EmptyProjection;
    let refusal = with_cx(|cx| false_variable_empty.admit_against(&variable_graph, cx))
        .expect_err("a variable-only projection is not empty");
    assert_receipt_rule(&refusal, CausalReceiptRule::OutcomeAxisMismatch);
}

#[test]
fn g0_non_saturating_maximum_claim_is_bound_to_graph_domain_and_witness() {
    let (machine, owner, clock) = minimal_machine();
    let e0 = equation(&owner, &clock, "test/equation/a", 40);
    let e1 = equation(&owner, &clock, "test/equation/b", 41);
    let v0 = variable(&owner, &clock, "test/variable/a", 42);
    let v1 = variable(&owner, &clock, "test/variable/b", 43);
    let i0 = incidence(&e0, &v0);
    let graph = with_cx(|cx| {
        CausalGraphDraft {
            units: CausalUnitConvention::SiBaseDimensions,
            scope: CausalGraphScope::CompleteMachineModel,
            extraction: extraction(10),
            equations: vec![e0, e1],
            variables: vec![v0, v1],
            conditions: Vec::new(),
            incidences: vec![i0],
        }
        .admit_against(&machine, cx)
    })
    .expect("sparse graph admits");
    let pair = CausalMatchingPair {
        incidence: graph.incidences()[0].id.clone(),
        equation: graph.incidences()[0].equation.clone(),
        variable: DerivativeVariableKey {
            variable: graph.incidences()[0].variable.clone(),
            derivative_order: 0,
        },
    };
    let domain = CausalReceiptDomain::UnconditionalGraph;
    let maximum = maximum_matching_binding(
        &graph,
        domain.clone(),
        core::slice::from_ref(&pair),
        cref!(MaximumMatchingCertificateRef, "test/maximum-matching", 44),
        cref!(CausalCheckerRef, "test/matching-checker", 45),
    )
    .expect("bound maximum theorem");
    let unmatched_equation = graph
        .equations()
        .iter()
        .find(|equation| equation.id != pair.equation)
        .expect("second equation")
        .id
        .clone();
    let unmatched_variable = graph
        .variables()
        .iter()
        .find(|variable| variable.id != pair.variable.variable)
        .expect("second variable")
        .id
        .clone();
    let draft = CausalizationReceiptDraft {
        structure: graph.structure_identity_receipt(),
        artifact: graph.artifact_identity_receipt(),
        analysis: analysis(),
        domain,
        determination: DeterminationClass::Mixed,
        structural_rank: StructuralRankState::Deficient,
        conditionality: Conditionality::Unconditional,
        matching: vec![pair],
        unmatched_equations: vec![unmatched_equation],
        unmatched_variables: vec![DerivativeVariableKey {
            variable: unmatched_variable,
            derivative_order: 0,
        }],
        conditional_outcomes: Vec::new(),
        maximum_matching_certificate: Some(maximum),
        conditional_coverage: None,
        unknown_axes: Vec::new(),
        evidence: CausalReceiptEvidence::CheckerReferenced(cref!(
            CausalCheckerRef,
            "test/matching-checker",
            45
        )),
    };
    let mut checker_substitution = draft.clone();
    checker_substitution.evidence = CausalReceiptEvidence::CheckerReferenced(cref!(
        CausalCheckerRef,
        "test/unrelated-checker",
        46
    ));
    let refusal = with_cx(|cx| checker_substitution.admit_against(&graph, cx))
        .expect_err("unrelated checker cannot validate the bound theorem");
    assert_receipt_rule(&refusal, CausalReceiptRule::OutcomeAxisMismatch);
    with_cx(|cx| draft.admit_against(&graph, cx)).expect("bound maximum claim admits");
}

fn branch(condition: &ActivationConditionRef, branch: &ActivationBranchRef) -> ActivationDomain {
    ActivationDomain::Conditional {
        cubes: vec![ActivationCube {
            selections: vec![ConditionBranchSelection {
                condition: condition.clone(),
                branch: branch.clone(),
            }],
        }],
    }
}

fn hybrid_draft(owner: &SubsystemId, clock: &ClockId) -> CausalGraphDraft {
    let equation = equation(owner, clock, "test/equation/hybrid", 50);
    let mut unknown = variable(owner, clock, "test/variable/hybrid", 51);
    unknown.solve_participation = SolveParticipation::ModeDependent;
    let mut parameter = variable(owner, clock, "test/variable/mode", 52);
    parameter.role = VariableRole::Parameter;
    parameter.solve_participation = SolveParticipation::ConditionOnly;
    let condition = cref!(ActivationConditionRef, "test/condition/mode", 53);
    let branch_a = cref!(ActivationBranchRef, "test/branch/a", 54);
    let branch_b = cref!(ActivationBranchRef, "test/branch/b", 55);
    let source = cref!(SourceArtifactRef, "test/mode-predicate", 56);
    let condition_spec = ActivationConditionSpec {
        condition: condition.clone(),
        source: ActivationConditionSource::AuditedPredicate(AuditedEscapeHatch {
            source: source.clone(),
            audit: cref!(EscapeHatchAuditRef, "test/mode-predicate-audit", 57),
            audited_source: source,
        }),
        branches: vec![branch_b.clone(), branch_a.clone()],
        dependencies: vec![parameter.id.clone()],
    };
    let unknown_incidence = incidence_spec(
        equation.id.clone(),
        unknown.id.clone(),
        0,
        SolveParticipation::Unknown,
        Dims::NONE,
        equation.residual.clone(),
        None,
        IncidenceClockRelation::SameClock,
        branch(&condition, &branch_a),
    )
    .expect("unknown branch incidence");
    let known_incidence = incidence_spec(
        equation.id.clone(),
        unknown.id.clone(),
        0,
        SolveParticipation::KnownInput,
        Dims::NONE,
        equation.residual.clone(),
        None,
        IncidenceClockRelation::SameClock,
        branch(&condition, &branch_b),
    )
    .expect("known branch incidence");
    CausalGraphDraft {
        units: CausalUnitConvention::SiBaseDimensions,
        scope: CausalGraphScope::CompleteMachineModel,
        extraction: extraction(10),
        equations: vec![equation],
        variables: vec![unknown, parameter],
        conditions: vec![condition_spec],
        incidences: vec![known_incidence, unknown_incidence],
    }
}

fn hybrid_graph(
    machine: &AdmittedMachineGraph,
    owner: &SubsystemId,
    clock: &ClockId,
) -> AdmittedCausalGraph {
    with_cx(|cx| hybrid_draft(owner, clock).admit_against(machine, cx))
        .expect("hybrid graph admits")
}

fn cartesian_hybrid_graph(
    machine: &AdmittedMachineGraph,
    owner: &SubsystemId,
    clock: &ClockId,
    secondary_branches: usize,
) -> AdmittedCausalGraph {
    let mut draft = hybrid_draft(owner, clock);
    let condition = cref!(
        ActivationConditionRef,
        "test/condition/z-secondary-cartesian",
        104
    );
    let source = cref!(SourceArtifactRef, "test/secondary-cartesian-source", 105);
    let branches = (0..secondary_branches)
        .map(|index| {
            ActivationBranchRef::new(
                format!("test/branch/z-secondary-{index:05}"),
                nz(1),
                [106; 32],
            )
            .expect("secondary branch reference")
        })
        .collect::<Vec<_>>();
    let mut marker = variable(owner, clock, "test/variable/secondary-mode-marker", 112);
    marker.role = VariableRole::Parameter;
    marker.solve_participation = SolveParticipation::KnownInput;
    marker.activation = branch(
        &condition,
        branches
            .first()
            .expect("Cartesian fixture has a secondary branch"),
    );
    draft.variables.push(marker);
    draft.conditions.push(ActivationConditionSpec {
        condition,
        source: ActivationConditionSource::AuditedPredicate(AuditedEscapeHatch {
            source: source.clone(),
            audit: cref!(EscapeHatchAuditRef, "test/secondary-cartesian-audit", 107),
            audited_source: source,
        }),
        branches,
        dependencies: vec![draft.variables[1].id.clone()],
    });
    with_cx(|cx| draft.admit_against(machine, cx)).expect("Cartesian hybrid graph admits")
}

fn cartesian_mode_assignment(
    graph: &AdmittedCausalGraph,
    ordinal: usize,
) -> Vec<ConditionBranchSelection> {
    let total = graph
        .conditions()
        .iter()
        .fold(1usize, |product, condition| {
            product
                .checked_mul(condition.branches.len())
                .expect("test Cartesian domain fits usize")
        });
    assert!(ordinal < total, "test Cartesian ordinal is in range");
    graph
        .conditions()
        .iter()
        .enumerate()
        .map(|(condition_index, condition)| {
            let stride =
                graph.conditions()[condition_index + 1..]
                    .iter()
                    .fold(1usize, |product, later| {
                        product
                            .checked_mul(later.branches.len())
                            .expect("test Cartesian stride fits usize")
                    });
            ConditionBranchSelection {
                condition: condition.condition.clone(),
                branch: condition.branches[(ordinal / stride) % condition.branches.len()].clone(),
            }
        })
        .collect()
}

#[test]
fn g5_nested_dnf_and_condition_table_permutations_do_not_move_identity() {
    let (machine, owner, clock) = minimal_machine();
    let mut draft = hybrid_draft(&owner, &clock);
    let parameter = draft.variables[1].id.clone();
    let condition = cref!(ActivationConditionRef, "test/condition/secondary", 66);
    let branch_c = cref!(ActivationBranchRef, "test/branch/c", 67);
    let branch_d = cref!(ActivationBranchRef, "test/branch/d", 68);
    let source = cref!(SourceArtifactRef, "test/secondary-predicate", 69);
    draft.conditions.push(ActivationConditionSpec {
        condition: condition.clone(),
        source: ActivationConditionSource::AuditedPredicate(AuditedEscapeHatch {
            source: source.clone(),
            audit: cref!(EscapeHatchAuditRef, "test/secondary-predicate-audit", 70),
            audited_source: source,
        }),
        branches: vec![branch_d.clone(), branch_c.clone()],
        dependencies: vec![parameter],
    });

    for incidence in &mut draft.incidences {
        let primary = match &incidence.activation {
            ActivationDomain::Conditional { cubes } => cubes[0].selections[0].clone(),
            ActivationDomain::Always => panic!("hybrid fixture incidence is conditional"),
        };
        let activation = ActivationDomain::Conditional {
            cubes: vec![
                ActivationCube {
                    selections: vec![
                        ConditionBranchSelection {
                            condition: condition.clone(),
                            branch: branch_d.clone(),
                        },
                        primary.clone(),
                    ],
                },
                ActivationCube {
                    selections: vec![
                        primary,
                        ConditionBranchSelection {
                            condition: condition.clone(),
                            branch: branch_c.clone(),
                        },
                    ],
                },
            ],
        };
        *incidence = incidence_spec(
            incidence.equation.clone(),
            incidence.variable.clone(),
            incidence.derivative_order,
            incidence.solve_participation,
            incidence.coefficient_dimensions,
            incidence.term.clone(),
            incidence.operator.clone(),
            incidence.clock_relation.clone(),
            activation,
        )
        .expect("two-condition incidence identity");
    }

    let replay = draft.clone();
    let baseline = with_cx(|cx| draft.admit_against(&machine, cx))
        .expect("canonical two-condition graph admits");
    let mut permuted = replay;
    permuted.conditions.reverse();
    for condition in &mut permuted.conditions {
        condition.branches.reverse();
        condition.dependencies.reverse();
    }
    permuted.incidences.reverse();
    for incidence in &mut permuted.incidences {
        let ActivationDomain::Conditional { cubes } = &mut incidence.activation else {
            panic!("two-condition fixture incidence is conditional");
        };
        cubes.reverse();
        for cube in cubes {
            cube.selections.reverse();
        }
    }
    let reordered = with_cx(|cx| permuted.admit_against(&machine, cx))
        .expect("nested-permuted two-condition graph admits");
    assert_eq!(
        baseline.structure_identity(),
        reordered.structure_identity(),
        "nested DNF or condition-table order moved normalized structure identity"
    );
    assert_eq!(
        baseline.artifact_identity(),
        reordered.artifact_identity(),
        "nested DNF or condition-table order moved provenance artifact identity"
    );
}

fn mode_assignment(
    graph: &AdmittedCausalGraph,
    branch_index: usize,
) -> Vec<ConditionBranchSelection> {
    vec![ConditionBranchSelection {
        condition: graph.conditions()[0].condition.clone(),
        branch: graph.conditions()[0].branches[branch_index].clone(),
    }]
}

#[test]
fn g0_mode_dependent_participation_requires_total_finite_domain_coverage() {
    let (machine, owner, clock) = minimal_machine();
    let mut draft = hybrid_draft(&owner, &clock);
    draft
        .incidences
        .retain(|incidence| incidence.solve_participation == SolveParticipation::Unknown);
    let refusal = with_cx(|cx| draft.admit_against(&machine, cx))
        .expect_err("one uncovered mode must refuse");
    assert_graph_rule(&refusal, CausalGraphRule::DerivativeParticipationMismatch);
}

#[test]
fn g0_activation_condition_dependencies_must_be_globally_available() {
    let (machine, owner, clock) = minimal_machine();
    let mut draft = hybrid_draft(&owner, &clock);
    let condition = draft.conditions[0].condition.clone();
    let branch = draft.conditions[0].branches[0].clone();
    draft.variables[1].activation = crate::branch(&condition, &branch);
    let refusal = with_cx(|cx| draft.admit_against(&machine, cx))
        .expect_err("a global condition cannot depend on a conditionally unavailable value");
    assert_graph_rule(&refusal, CausalGraphRule::InvalidActivationCondition);
}

#[test]
fn g0_guard_backed_conditions_require_always_available_guard_and_incidences() {
    let (machine, owner, clock) = minimal_machine();
    let mut draft = hybrid_draft(&owner, &clock);
    let condition = draft.conditions[0].condition.clone();
    let dependency = draft.variables[1].clone();
    let mut guard = equation(&owner, &clock, "test/equation/mode-guard", 119);
    guard.role = EquationRole::Guard;
    guard.solve_participation = EquationParticipation::ConditionOnly;
    let guard_incidence = incidence_spec(
        guard.id.clone(),
        dependency.id,
        0,
        SolveParticipation::ConditionOnly,
        Dims::NONE,
        guard.residual.clone(),
        None,
        IncidenceClockRelation::SameClock,
        ActivationDomain::Always,
    )
    .expect("guard dependency incidence");
    let guard_id = guard.id.clone();
    draft.conditions[0].source = ActivationConditionSource::GuardEquation {
        equation: guard.id.clone(),
        obligation: cref!(
            GuardSolveObligationRef,
            "test/guard-root-solve-obligation",
            120
        ),
    };
    draft.equations.push(guard);
    draft.incidences.push(guard_incidence);

    with_cx(|cx| draft.clone().admit_against(&machine, cx))
        .expect("always-available guard-backed condition admits structurally");
    let conditional_branch = draft.conditions[0].branches[0].clone();

    let mut conditional_incidence_draft = draft.clone();
    let guard_incidence = conditional_incidence_draft
        .incidences
        .iter_mut()
        .find(|incidence| incidence.equation == guard_id)
        .expect("guard fixture incidence");
    *guard_incidence = incidence_spec(
        guard_incidence.equation.clone(),
        guard_incidence.variable.clone(),
        guard_incidence.derivative_order,
        guard_incidence.solve_participation,
        guard_incidence.coefficient_dimensions,
        guard_incidence.term.clone(),
        guard_incidence.operator.clone(),
        guard_incidence.clock_relation.clone(),
        branch(&condition, &conditional_branch),
    )
    .expect("conditional guard incidence identity");
    let refusal = with_cx(|cx| conditional_incidence_draft.admit_against(&machine, cx))
        .expect_err("a global condition cannot depend on a conditionally available guard edge");
    assert_graph_rule(&refusal, CausalGraphRule::InvalidActivationCondition);

    let guard = draft
        .equations
        .iter_mut()
        .find(|equation| equation.role == EquationRole::Guard)
        .expect("guard fixture equation");
    guard.activation = branch(&condition, &conditional_branch);
    let refusal = with_cx(|cx| draft.admit_against(&machine, cx))
        .expect_err("a global condition cannot depend on a conditionally available guard");
    assert_graph_rule(&refusal, CausalGraphRule::InvalidActivationCondition);
}

fn mode_cell_receipt(
    graph: &AdmittedCausalGraph,
    assignment: Vec<ConditionBranchSelection>,
) -> AdmittedCausalizationReceipt {
    let active_unknown = graph
        .incidences()
        .iter()
        .find(|incidence| incidence.solve_participation == SolveParticipation::Unknown)
        .expect("unknown branch incidence");
    let unknown_branch = match &active_unknown.activation {
        ActivationDomain::Conditional { cubes } => cubes.iter().any(|cube| {
            cube.selections
                .iter()
                .all(|selection| assignment.binary_search(selection).is_ok())
        }),
        ActivationDomain::Always => unreachable!("fixture incidence is conditional"),
    };
    let (determination, rank, matching, unmatched_equations) = if unknown_branch {
        (
            DeterminationClass::WellDetermined,
            StructuralRankState::FullRelativeToMinSide,
            vec![CausalMatchingPair {
                incidence: active_unknown.id.clone(),
                equation: active_unknown.equation.clone(),
                variable: DerivativeVariableKey {
                    variable: active_unknown.variable.clone(),
                    derivative_order: 0,
                },
            }],
            Vec::new(),
        )
    } else {
        (
            DeterminationClass::OverDetermined,
            StructuralRankState::NotApplicable,
            Vec::new(),
            vec![graph.equations()[0].id.clone()],
        )
    };
    with_cx(|cx| {
        CausalizationReceiptDraft {
            structure: graph.structure_identity_receipt(),
            artifact: graph.artifact_identity_receipt(),
            analysis: analysis(),
            domain: CausalReceiptDomain::ModeCell { assignment },
            determination,
            structural_rank: rank,
            conditionality: Conditionality::Unconditional,
            matching,
            unmatched_equations,
            unmatched_variables: Vec::new(),
            conditional_outcomes: Vec::new(),
            maximum_matching_certificate: None,
            conditional_coverage: None,
            unknown_axes: Vec::new(),
            evidence: CausalReceiptEvidence::Unverified,
        }
        .admit_against(graph, cx)
    })
    .expect("mode-cell receipt admits")
}

fn incomplete_mode_cell_receipt(
    graph: &AdmittedCausalGraph,
    assignment: Vec<ConditionBranchSelection>,
) -> AdmittedCausalizationReceipt {
    let active_unknown = graph
        .incidences()
        .iter()
        .find(|incidence| incidence.solve_participation == SolveParticipation::Unknown)
        .expect("hybrid fixture has an unknown structural incidence");
    let unknown_active = match &active_unknown.activation {
        ActivationDomain::Always => true,
        ActivationDomain::Conditional { cubes } => cubes.iter().any(|cube| {
            cube.selections
                .iter()
                .all(|selection| assignment.binary_search(selection).is_ok())
        }),
    };
    with_cx(|cx| {
        CausalizationReceiptDraft {
            structure: graph.structure_identity_receipt(),
            artifact: graph.artifact_identity_receipt(),
            analysis: analysis(),
            domain: CausalReceiptDomain::ModeCell { assignment },
            determination: DeterminationClass::Unknown,
            structural_rank: StructuralRankState::Unknown,
            conditionality: Conditionality::Unconditional,
            matching: Vec::new(),
            unmatched_equations: vec![active_unknown.equation.clone()],
            unmatched_variables: if unknown_active {
                vec![DerivativeVariableKey {
                    variable: active_unknown.variable.clone(),
                    derivative_order: active_unknown.derivative_order,
                }]
            } else {
                Vec::new()
            },
            conditional_outcomes: Vec::new(),
            maximum_matching_certificate: None,
            conditional_coverage: None,
            unknown_axes: vec![
                CausalUnknownAxisState {
                    axis: CausalOutcomeAxis::Determination,
                    reason: CausalUnknownReason::IncompleteMetadata,
                    resume_checkpoint: None,
                },
                CausalUnknownAxisState {
                    axis: CausalOutcomeAxis::StructuralRank,
                    reason: CausalUnknownReason::IncompleteMetadata,
                    resume_checkpoint: None,
                },
            ],
            evidence: CausalReceiptEvidence::Unverified,
        }
        .admit_against(graph, cx)
    })
    .expect("incomplete hybrid mode-cell receipt admits")
}

#[test]
fn g3_hybrid_summary_uses_typed_children_and_bound_coverage() {
    let (machine, owner, clock) = minimal_machine();
    let graph = hybrid_graph(&machine, &owner, &clock);
    let child_a = mode_cell_receipt(&graph, mode_assignment(&graph, 0));
    let child_b = mode_cell_receipt(&graph, mode_assignment(&graph, 1));
    let outcomes = vec![
        conditional_outcome(&child_b).expect("mode child"),
        conditional_outcome(&child_a).expect("mode child"),
    ];
    let coverage = conditional_coverage_binding(
        &graph,
        &outcomes,
        cref!(ConditionalCoverageRef, "test/hybrid-coverage", 58),
        cref!(CausalCheckerRef, "test/hybrid-checker", 59),
    )
    .expect("coverage binding");
    let summary = CausalizationReceiptDraft {
        structure: graph.structure_identity_receipt(),
        artifact: graph.artifact_identity_receipt(),
        analysis: analysis(),
        domain: CausalReceiptDomain::HybridSummary,
        determination: DeterminationClass::Unknown,
        structural_rank: StructuralRankState::Unknown,
        conditionality: Conditionality::Conditional,
        matching: Vec::new(),
        unmatched_equations: Vec::new(),
        unmatched_variables: Vec::new(),
        conditional_outcomes: outcomes,
        maximum_matching_certificate: None,
        conditional_coverage: Some(coverage),
        unknown_axes: vec![
            CausalUnknownAxisState {
                axis: CausalOutcomeAxis::StructuralRank,
                reason: CausalUnknownReason::NonUniformAcrossModes,
                resume_checkpoint: None,
            },
            CausalUnknownAxisState {
                axis: CausalOutcomeAxis::Determination,
                reason: CausalUnknownReason::NonUniformAcrossModes,
                resume_checkpoint: None,
            },
        ],
        evidence: CausalReceiptEvidence::CheckerReferenced(cref!(
            CausalCheckerRef,
            "test/hybrid-checker",
            59
        )),
    };
    let mut checker_substitution = summary.clone();
    checker_substitution.evidence = CausalReceiptEvidence::CheckerReferenced(cref!(
        CausalCheckerRef,
        "test/unrelated-hybrid-checker",
        63
    ));
    let refusal = with_cx(|cx| checker_substitution.admit_against(&graph, cx))
        .expect_err("unrelated checker cannot validate coverage");
    assert_receipt_rule(&refusal, CausalReceiptRule::ConditionalCoverageMismatch);
    let receipt = with_cx(|cx| summary.admit_against(&graph, cx))
        .expect("heterogeneous hybrid summary admits honestly");
    assert_eq!(receipt.conditional_outcomes().len(), 2);
    assert_eq!(receipt.unknown_axes().len(), 2);
}

#[test]
#[allow(clippy::too_many_lines)] // Both concrete on/off cells and their summary stay visible.
fn g3_hybrid_summary_retains_a_concrete_empty_off_mode() {
    let (machine, owner, clock) = minimal_machine();
    let mut draft = hybrid_draft(&owner, &clock);
    let active_domain = draft
        .incidences
        .iter()
        .find(|incidence| incidence.solve_participation == SolveParticipation::Unknown)
        .expect("hybrid fixture has an active unknown incidence")
        .activation
        .clone();
    draft.equations[0].activation = active_domain.clone();
    draft.variables[0].activation = active_domain;
    draft
        .incidences
        .retain(|incidence| incidence.solve_participation == SolveParticipation::Unknown);
    let graph = with_cx(|cx| draft.admit_against(&machine, cx))
        .expect("hybrid graph with a structurally empty off mode admits");
    let incidence = graph
        .incidences()
        .first()
        .expect("on mode retains its structural incidence");
    let active_branch = match &incidence.activation {
        ActivationDomain::Conditional { cubes } => cubes[0].selections[0].branch.clone(),
        ActivationDomain::Always => panic!("off-mode fixture incidence is conditional"),
    };

    let mut children = Vec::new();
    for ordinal in 0..2 {
        let assignment = mode_assignment(&graph, ordinal);
        let active = assignment[0].branch == active_branch;
        let child = with_cx(|cx| {
            CausalizationReceiptDraft {
                structure: graph.structure_identity_receipt(),
                artifact: graph.artifact_identity_receipt(),
                analysis: analysis(),
                domain: CausalReceiptDomain::ModeCell { assignment },
                determination: if active {
                    DeterminationClass::WellDetermined
                } else {
                    DeterminationClass::EmptyProjection
                },
                structural_rank: if active {
                    StructuralRankState::FullRelativeToMinSide
                } else {
                    StructuralRankState::NotApplicable
                },
                conditionality: Conditionality::Unconditional,
                matching: active
                    .then(|| CausalMatchingPair {
                        incidence: incidence.id.clone(),
                        equation: incidence.equation.clone(),
                        variable: DerivativeVariableKey {
                            variable: incidence.variable.clone(),
                            derivative_order: incidence.derivative_order,
                        },
                    })
                    .into_iter()
                    .collect(),
                unmatched_equations: Vec::new(),
                unmatched_variables: Vec::new(),
                conditional_outcomes: Vec::new(),
                maximum_matching_certificate: None,
                conditional_coverage: None,
                unknown_axes: Vec::new(),
                evidence: CausalReceiptEvidence::Unverified,
            }
            .admit_against(&graph, cx)
        })
        .expect("both on and empty off mode cells admit concretely");
        children.push(child);
    }
    assert!(children.iter().any(|child| {
        child.determination() == DeterminationClass::EmptyProjection
            && child.structural_rank() == StructuralRankState::NotApplicable
    }));
    let outcomes = children
        .iter()
        .map(|child| conditional_outcome(child).expect("concrete mode child"))
        .collect::<Vec<_>>();
    let checker = cref!(CausalCheckerRef, "test/off-mode-summary-checker", 120);
    let coverage = conditional_coverage_binding(
        &graph,
        &outcomes,
        cref!(
            ConditionalCoverageRef,
            "test/off-mode-summary-coverage",
            121
        ),
        checker.clone(),
    )
    .expect("complete coverage accepts the concrete empty mode");
    let summary = CausalizationReceiptDraft {
        structure: graph.structure_identity_receipt(),
        artifact: graph.artifact_identity_receipt(),
        analysis: analysis(),
        domain: CausalReceiptDomain::HybridSummary,
        determination: DeterminationClass::Unknown,
        structural_rank: StructuralRankState::Unknown,
        conditionality: Conditionality::Conditional,
        matching: Vec::new(),
        unmatched_equations: Vec::new(),
        unmatched_variables: Vec::new(),
        conditional_outcomes: outcomes,
        maximum_matching_certificate: None,
        conditional_coverage: Some(coverage),
        unknown_axes: vec![
            CausalUnknownAxisState {
                axis: CausalOutcomeAxis::Determination,
                reason: CausalUnknownReason::NonUniformAcrossModes,
                resume_checkpoint: None,
            },
            CausalUnknownAxisState {
                axis: CausalOutcomeAxis::StructuralRank,
                reason: CausalUnknownReason::NonUniformAcrossModes,
                resume_checkpoint: None,
            },
        ],
        evidence: CausalReceiptEvidence::CheckerReferenced(checker),
    };
    let admitted = with_cx(|cx| summary.admit_against(&graph, cx))
        .expect("hybrid summary retains both the on and empty off mode");
    assert_eq!(admitted.conditional_outcomes().len(), 2);
    assert!(
        admitted.conditional_outcomes().iter().any(|outcome| {
            outcome.determination() == DeterminationClass::EmptyProjection
                && outcome.structural_rank() == StructuralRankState::NotApplicable
                && outcome.unknown_axes().is_empty()
        }),
        "the admitted summary must retain the concrete empty off-mode child"
    );
}

#[test]
#[allow(clippy::too_many_lines)] // Progress mutation, coverage refusal, and parent identity are one law.
fn g3_unknown_hybrid_summary_retains_child_resume_state_without_claiming_coverage() {
    let (machine, owner, clock) = minimal_machine();
    let graph = hybrid_graph(&machine, &owner, &clock);
    let incidence = graph
        .incidences()
        .iter()
        .find(|incidence| incidence.solve_participation == SolveParticipation::Unknown)
        .expect("hybrid fixture has an unknown structural incidence");
    let active_branch = match &incidence.activation {
        ActivationDomain::Conditional { cubes } => cubes[0].selections[0].branch.clone(),
        ActivationDomain::Always => panic!("hybrid fixture incidence is conditional"),
    };
    let assignment = (0..2)
        .map(|ordinal| mode_assignment(&graph, ordinal))
        .find(|assignment| assignment[0].branch == active_branch)
        .expect("one mode activates the unknown incidence");
    let child_draft = CausalizationReceiptDraft {
        structure: graph.structure_identity_receipt(),
        artifact: graph.artifact_identity_receipt(),
        analysis: analysis(),
        domain: CausalReceiptDomain::ModeCell { assignment },
        determination: DeterminationClass::Unknown,
        structural_rank: StructuralRankState::Unknown,
        conditionality: Conditionality::Unconditional,
        matching: Vec::new(),
        unmatched_equations: vec![incidence.equation.clone()],
        unmatched_variables: vec![DerivativeVariableKey {
            variable: incidence.variable.clone(),
            derivative_order: incidence.derivative_order,
        }],
        conditional_outcomes: Vec::new(),
        maximum_matching_certificate: None,
        conditional_coverage: None,
        unknown_axes: vec![
            CausalUnknownAxisState {
                axis: CausalOutcomeAxis::Determination,
                reason: CausalUnknownReason::BudgetExhausted,
                resume_checkpoint: Some(cref!(
                    CausalCheckpointRef,
                    "test/mode-determination-resume",
                    126
                )),
            },
            CausalUnknownAxisState {
                axis: CausalOutcomeAxis::StructuralRank,
                reason: CausalUnknownReason::BudgetExhausted,
                resume_checkpoint: Some(cref!(CausalCheckpointRef, "test/mode-rank-resume", 127)),
            },
        ],
        evidence: CausalReceiptEvidence::Unverified,
    };
    let mut changed_progress = child_draft.clone();
    changed_progress.unknown_axes[0].reason = CausalUnknownReason::Cancelled;
    changed_progress.unknown_axes[0].resume_checkpoint = Some(cref!(
        CausalCheckpointRef,
        "test/mode-determination-cancelled",
        128
    ));
    let child_a = with_cx(|cx| child_draft.admit_against(&graph, cx))
        .expect("incomplete mode-cell receipt admits with exact resume state");
    let child_b = with_cx(|cx| changed_progress.admit_against(&graph, cx))
        .expect("alternate honest progress state admits");
    assert_eq!(
        child_a.outcome_identity(),
        child_b.outcome_identity(),
        "progress metadata is excluded from producer-independent outcome semantics"
    );
    assert_eq!(
        child_a.outcome_identity_receipt().canonical_preimage(),
        child_b.outcome_identity_receipt().canonical_preimage(),
        "progress-only mutations must preserve the normalized outcome bytes"
    );
    assert_ne!(
        child_a.identity(),
        child_b.identity(),
        "progress metadata remains bound by the full evidence receipt"
    );
    assert_ne!(
        child_a.identity_receipt().canonical_preimage(),
        child_b.identity_receipt().canonical_preimage(),
        "progress-only mutations must move the complete receipt bytes"
    );
    let outcome_a = conditional_outcome(&child_a).expect("incomplete mode child projects");
    let outcome_b = conditional_outcome(&child_b).expect("alternate mode child projects");
    assert_eq!(outcome_a.unknown_axes(), child_a.unknown_axes());
    assert_eq!(outcome_b.unknown_axes(), child_b.unknown_axes());
    assert_eq!(
        conditional_coverage_binding(
            &graph,
            core::slice::from_ref(&outcome_a),
            cref!(ConditionalCoverageRef, "test/incomplete-mode-coverage", 129),
            cref!(CausalCheckerRef, "test/incomplete-mode-checker", 130),
        ),
        Err(ConditionalCoverageBindingError::NonConcreteChild { outcome_index: 0 }),
        "an incomplete child must never enter a complete-coverage commitment"
    );

    let mut foreign_draft = hybrid_draft(&owner, &clock);
    foreign_draft.extraction = extraction(131);
    let foreign_graph = with_cx(|cx| foreign_draft.admit_against(&machine, cx))
        .expect("provenance-distinct foreign hybrid graph admits");
    let foreign_child = mode_cell_receipt(&foreign_graph, mode_assignment(&foreign_graph, 0));
    let foreign_outcome = conditional_outcome(&foreign_child).expect("foreign mode child");
    let first_order = [outcome_a.clone(), foreign_outcome.clone()];
    let reverse_order = [foreign_outcome, outcome_a.clone()];
    let invalid_a = conditional_coverage_binding(
        &graph,
        &first_order,
        cref!(ConditionalCoverageRef, "test/invalid-order-coverage", 132),
        cref!(CausalCheckerRef, "test/invalid-order-checker", 133),
    )
    .expect_err("foreign and nonconcrete children refuse complete coverage");
    let invalid_b = conditional_coverage_binding(
        &graph,
        &reverse_order,
        cref!(ConditionalCoverageRef, "test/invalid-order-coverage", 132),
        cref!(CausalCheckerRef, "test/invalid-order-checker", 133),
    )
    .expect_err("reversed invalid child multiset also refuses");
    assert_eq!(
        invalid_a, invalid_b,
        "invalid coverage diagnostics must use fixed precedence and canonical indices"
    );
    assert!(matches!(
        invalid_a,
        ConditionalCoverageBindingError::ForeignGraph { .. }
    ));

    let summary_for = |outcome: ConditionalCausalOutcome| CausalizationReceiptDraft {
        structure: graph.structure_identity_receipt(),
        artifact: graph.artifact_identity_receipt(),
        analysis: analysis(),
        domain: CausalReceiptDomain::HybridSummary,
        determination: DeterminationClass::Unknown,
        structural_rank: StructuralRankState::Unknown,
        conditionality: Conditionality::Unknown,
        matching: Vec::new(),
        unmatched_equations: Vec::new(),
        unmatched_variables: Vec::new(),
        conditional_outcomes: vec![outcome],
        maximum_matching_certificate: None,
        conditional_coverage: None,
        unknown_axes: vec![
            CausalUnknownAxisState {
                axis: CausalOutcomeAxis::Determination,
                reason: CausalUnknownReason::IncompleteMetadata,
                resume_checkpoint: None,
            },
            CausalUnknownAxisState {
                axis: CausalOutcomeAxis::StructuralRank,
                reason: CausalUnknownReason::IncompleteMetadata,
                resume_checkpoint: None,
            },
            CausalUnknownAxisState {
                axis: CausalOutcomeAxis::Conditionality,
                reason: CausalUnknownReason::IncompleteMetadata,
                resume_checkpoint: None,
            },
        ],
        evidence: CausalReceiptEvidence::Unverified,
    };
    let summary_a = with_cx(|cx| summary_for(outcome_a).admit_against(&graph, cx))
        .expect("unknown summary retains incomplete child without coverage authority");
    let summary_b = with_cx(|cx| summary_for(outcome_b).admit_against(&graph, cx))
        .expect("unknown summary retains alternate progress state");
    assert_eq!(summary_a.outcome_identity(), summary_b.outcome_identity());
    assert_eq!(
        summary_a.outcome_identity_receipt().canonical_preimage(),
        summary_b.outcome_identity_receipt().canonical_preimage()
    );
    assert_ne!(summary_a.identity(), summary_b.identity());
    assert_ne!(
        summary_a.identity_receipt().canonical_preimage(),
        summary_b.identity_receipt().canonical_preimage()
    );
    assert_eq!(summary_a.conditional_outcomes()[0].unknown_axes().len(), 2);
}

#[test]
#[allow(clippy::too_many_lines)] // Concrete and incomplete child semantics are one aggregation law.
fn g3_partial_hybrid_summary_does_not_confuse_missing_evidence_with_nonuniformity() {
    let (machine, owner, clock) = minimal_machine();
    let graph = hybrid_graph(&machine, &owner, &clock);
    let active_unknown = graph
        .incidences()
        .iter()
        .find(|incidence| incidence.solve_participation == SolveParticipation::Unknown)
        .expect("hybrid fixture has an unknown structural incidence");
    let active_branch = match &active_unknown.activation {
        ActivationDomain::Conditional { cubes } => cubes[0].selections[0].branch.clone(),
        ActivationDomain::Always => panic!("hybrid fixture incidence is conditional"),
    };
    let (active_assignment, inactive_assignment) = (0..2)
        .map(|ordinal| mode_assignment(&graph, ordinal))
        .partition::<Vec<_>, _>(|assignment| assignment[0].branch == active_branch);
    let active_assignment = active_assignment
        .into_iter()
        .next()
        .expect("one mode activates the unknown incidence");
    let inactive_assignment = inactive_assignment
        .into_iter()
        .next()
        .expect("one mode excludes the unknown incidence");

    let concrete_child = mode_cell_receipt(&graph, active_assignment);
    let incomplete_child = with_cx(|cx| {
        CausalizationReceiptDraft {
            structure: graph.structure_identity_receipt(),
            artifact: graph.artifact_identity_receipt(),
            analysis: analysis(),
            domain: CausalReceiptDomain::ModeCell {
                assignment: inactive_assignment,
            },
            determination: DeterminationClass::Unknown,
            structural_rank: StructuralRankState::Unknown,
            conditionality: Conditionality::Unconditional,
            matching: Vec::new(),
            unmatched_equations: vec![graph.equations()[0].id.clone()],
            unmatched_variables: Vec::new(),
            conditional_outcomes: Vec::new(),
            maximum_matching_certificate: None,
            conditional_coverage: None,
            unknown_axes: vec![
                CausalUnknownAxisState {
                    axis: CausalOutcomeAxis::Determination,
                    reason: CausalUnknownReason::IncompleteMetadata,
                    resume_checkpoint: None,
                },
                CausalUnknownAxisState {
                    axis: CausalOutcomeAxis::StructuralRank,
                    reason: CausalUnknownReason::IncompleteMetadata,
                    resume_checkpoint: None,
                },
            ],
            evidence: CausalReceiptEvidence::Unverified,
        }
        .admit_against(&graph, cx)
    })
    .expect("honestly incomplete inactive mode admits");
    let CausalReceiptDomain::ModeCell {
        assignment: concrete_assignment,
    } = concrete_child.domain()
    else {
        panic!("fixture child is a mode cell");
    };
    let expected_concrete_assignment = concrete_assignment.clone();
    let CausalReceiptDomain::ModeCell {
        assignment: incomplete_assignment,
    } = incomplete_child.domain()
    else {
        panic!("fixture child is a mode cell");
    };
    let expected_incomplete_assignment = incomplete_assignment.clone();
    let outcomes = vec![
        conditional_outcome(&concrete_child).expect("concrete mode child projects"),
        conditional_outcome(&incomplete_child).expect("incomplete mode child projects"),
    ];
    assert!(
        outcomes
            .iter()
            .any(|outcome| outcome.determination() != DeterminationClass::Unknown)
    );
    assert!(
        outcomes
            .iter()
            .any(|outcome| outcome.determination() == DeterminationClass::Unknown)
    );
    assert!(matches!(
        conditional_coverage_binding(
            &graph,
            &outcomes,
            cref!(ConditionalCoverageRef, "test/partial-mode-coverage", 138),
            cref!(CausalCheckerRef, "test/partial-mode-checker", 139),
        ),
        Err(ConditionalCoverageBindingError::NonConcreteChild { .. })
    ));

    let summary = CausalizationReceiptDraft {
        structure: graph.structure_identity_receipt(),
        artifact: graph.artifact_identity_receipt(),
        analysis: analysis(),
        domain: CausalReceiptDomain::HybridSummary,
        determination: DeterminationClass::Unknown,
        structural_rank: StructuralRankState::Unknown,
        conditionality: Conditionality::Unknown,
        matching: Vec::new(),
        unmatched_equations: Vec::new(),
        unmatched_variables: Vec::new(),
        conditional_outcomes: outcomes,
        maximum_matching_certificate: None,
        conditional_coverage: None,
        unknown_axes: vec![
            CausalUnknownAxisState {
                axis: CausalOutcomeAxis::Determination,
                reason: CausalUnknownReason::IncompleteMetadata,
                resume_checkpoint: None,
            },
            CausalUnknownAxisState {
                axis: CausalOutcomeAxis::StructuralRank,
                reason: CausalUnknownReason::IncompleteMetadata,
                resume_checkpoint: None,
            },
            CausalUnknownAxisState {
                axis: CausalOutcomeAxis::Conditionality,
                reason: CausalUnknownReason::IncompleteMetadata,
                resume_checkpoint: None,
            },
        ],
        evidence: CausalReceiptEvidence::Unverified,
    };
    let admitted = with_cx(|cx| summary.admit_against(&graph, cx)).expect(
        "partial summary retains concrete and incomplete children without inventing nonuniformity",
    );
    assert_eq!(admitted.conditional_outcomes().len(), 2);
    let retained_concrete = admitted
        .conditional_outcomes()
        .iter()
        .find(|outcome| outcome.assignment() == expected_concrete_assignment.as_slice())
        .expect("admitted summary retains the concrete assignment");
    assert_ne!(
        retained_concrete.determination(),
        DeterminationClass::Unknown
    );
    assert!(retained_concrete.unknown_axes().is_empty());
    let retained_incomplete = admitted
        .conditional_outcomes()
        .iter()
        .find(|outcome| outcome.assignment() == expected_incomplete_assignment.as_slice())
        .expect("admitted summary retains the incomplete assignment");
    assert_eq!(
        retained_incomplete.determination(),
        DeterminationClass::Unknown
    );
    assert_eq!(retained_incomplete.unknown_axes().len(), 2);
    assert!(
        retained_incomplete
            .unknown_axes()
            .iter()
            .all(|state| state.reason == CausalUnknownReason::IncompleteMetadata)
    );
    assert!(
        admitted
            .unknown_axes()
            .iter()
            .all(|state| state.reason != CausalUnknownReason::NonUniformAcrossModes),
        "one concrete value plus missing evidence is not a disagreement theorem"
    );
}

#[test]
#[allow(clippy::too_many_lines)] // Agreeing and disagreeing concrete folds share one fixture.
fn g3_partial_hybrid_aggregation_uses_only_concrete_disagreements() {
    let (machine, owner, clock) = minimal_machine();
    let graph = cartesian_hybrid_graph(&machine, &owner, &clock, 3);
    let concrete_outcomes = (0..6)
        .map(|ordinal| {
            let child = mode_cell_receipt(&graph, cartesian_mode_assignment(&graph, ordinal));
            conditional_outcome(&child).expect("concrete Cartesian mode child projects")
        })
        .collect::<Vec<_>>();
    let anchor = concrete_outcomes[0].clone();
    let agreeing = concrete_outcomes
        .iter()
        .filter(|outcome| {
            outcome.determination() == anchor.determination()
                && outcome.structural_rank() == anchor.structural_rank()
        })
        .take(3)
        .cloned()
        .collect::<Vec<_>>();
    assert_eq!(
        agreeing.len(),
        3,
        "the 2x3 fixture supplies three concrete cells with equal axes"
    );
    let disagreeing = concrete_outcomes
        .iter()
        .find(|outcome| {
            outcome.determination() != anchor.determination()
                && outcome.structural_rank() != anchor.structural_rank()
        })
        .cloned()
        .expect("the other primary mode supplies a concrete axis disagreement");
    let incomplete_child = incomplete_mode_cell_receipt(&graph, agreeing[2].assignment().to_vec());
    let incomplete =
        conditional_outcome(&incomplete_child).expect("incomplete Cartesian child projects");

    let summary_for = |conditional_outcomes: Vec<ConditionalCausalOutcome>, nonuniform: bool| {
        let aggregate_reason = if nonuniform {
            CausalUnknownReason::NonUniformAcrossModes
        } else {
            CausalUnknownReason::IncompleteMetadata
        };
        CausalizationReceiptDraft {
            structure: graph.structure_identity_receipt(),
            artifact: graph.artifact_identity_receipt(),
            analysis: analysis(),
            domain: CausalReceiptDomain::HybridSummary,
            determination: DeterminationClass::Unknown,
            structural_rank: StructuralRankState::Unknown,
            conditionality: Conditionality::Unknown,
            matching: Vec::new(),
            unmatched_equations: Vec::new(),
            unmatched_variables: Vec::new(),
            conditional_outcomes,
            maximum_matching_certificate: None,
            conditional_coverage: None,
            unknown_axes: vec![
                CausalUnknownAxisState {
                    axis: CausalOutcomeAxis::Determination,
                    reason: aggregate_reason,
                    resume_checkpoint: None,
                },
                CausalUnknownAxisState {
                    axis: CausalOutcomeAxis::StructuralRank,
                    reason: aggregate_reason,
                    resume_checkpoint: None,
                },
                CausalUnknownAxisState {
                    axis: CausalOutcomeAxis::Conditionality,
                    reason: CausalUnknownReason::IncompleteMetadata,
                    resume_checkpoint: None,
                },
            ],
            evidence: CausalReceiptEvidence::Unverified,
        }
    };

    let agreeing_summary = with_cx(|cx| {
        summary_for(
            vec![agreeing[0].clone(), agreeing[1].clone(), incomplete.clone()],
            false,
        )
        .admit_against(&graph, cx)
    })
    .expect("two agreeing concrete cells plus one unknown cell are not nonuniform");
    assert!(
        agreeing_summary
            .unknown_axes()
            .iter()
            .all(|state| state.reason != CausalUnknownReason::NonUniformAcrossModes)
    );

    let disagreeing_summary = with_cx(|cx| {
        summary_for(vec![agreeing[0].clone(), disagreeing, incomplete], true)
            .admit_against(&graph, cx)
    })
    .expect("two disagreeing concrete cells require nonuniformity despite an unknown peer");
    assert_eq!(
        disagreeing_summary
            .unknown_axes()
            .iter()
            .filter(|state| state.reason == CausalUnknownReason::NonUniformAcrossModes)
            .count(),
        2,
        "only the two disagreeing concrete axes are classified as nonuniform"
    );
}

#[test]
#[allow(clippy::too_many_lines)]
fn g0_cartesian_coverage_requires_every_mode_cell_and_concrete_child() {
    let (machine, owner, clock) = minimal_machine();
    let graph = cartesian_hybrid_graph(&machine, &owner, &clock, 3);
    assert_eq!(
        graph
            .conditions()
            .iter()
            .map(|condition| condition.branches.len())
            .collect::<Vec<_>>(),
        vec![2, 3],
        "fixture must exercise asymmetric 2x3 Cartesian strides"
    );
    let mut outcomes: Vec<_> = (0..6)
        .map(|ordinal| {
            let child = mode_cell_receipt(&graph, cartesian_mode_assignment(&graph, ordinal));
            conditional_outcome(&child).expect("Cartesian mode child")
        })
        .collect();
    outcomes.swap(0, 5);
    outcomes.swap(1, 3);
    let certificate = cref!(
        ConditionalCoverageRef,
        "test/cartesian-complete-coverage",
        108
    );
    let checker = cref!(CausalCheckerRef, "test/cartesian-complete-checker", 109);
    conditional_coverage_binding(&graph, &outcomes, certificate.clone(), checker.clone())
        .expect("shuffled complete 2x3 Cartesian cover binds");

    let mut missing = outcomes.clone();
    assert!(missing.pop().is_some());
    assert_eq!(
        conditional_coverage_binding(&graph, &missing, certificate.clone(), checker.clone(),),
        Err(ConditionalCoverageBindingError::IncompleteCartesianCover {
            submitted: 5,
            expected: 6,
        }),
        "a self-consistent proper subset is not a coverage theorem"
    );

    let mut duplicated = outcomes.clone();
    duplicated[5] = duplicated[0].clone();
    assert_eq!(
        conditional_coverage_binding(&graph, &duplicated, certificate.clone(), checker.clone(),),
        Err(ConditionalCoverageBindingError::DuplicateAssignment { index: 5 }),
        "a duplicate cell cannot replace an omitted Cartesian cell"
    );
    assert_eq!(
        conditional_coverage_binding(&graph, &[], certificate.clone(), checker.clone()),
        Err(ConditionalCoverageBindingError::EmptyOutcomeSet)
    );

    let foreign_graph = hybrid_graph(&machine, &owner, &clock);
    let foreign_child = mode_cell_receipt(&foreign_graph, mode_assignment(&foreign_graph, 0));
    let foreign_outcome = conditional_outcome(&foreign_child).expect("foreign mode child");
    assert_eq!(
        conditional_coverage_binding(
            &graph,
            core::slice::from_ref(&foreign_outcome),
            certificate.clone(),
            checker.clone(),
        ),
        Err(ConditionalCoverageBindingError::ForeignGraph { outcome_index: 0 })
    );

    let assignment = cartesian_mode_assignment(&graph, 0);
    let concrete = mode_cell_receipt(&graph, assignment.clone());
    let unknown = with_cx(|cx| {
        CausalizationReceiptDraft {
            structure: graph.structure_identity_receipt(),
            artifact: graph.artifact_identity_receipt(),
            analysis: analysis(),
            domain: CausalReceiptDomain::ModeCell { assignment },
            determination: DeterminationClass::Unknown,
            structural_rank: StructuralRankState::Unknown,
            conditionality: Conditionality::Unconditional,
            matching: concrete.matching().to_vec(),
            unmatched_equations: concrete.unmatched_equations().to_vec(),
            unmatched_variables: concrete.unmatched_variables().to_vec(),
            conditional_outcomes: Vec::new(),
            maximum_matching_certificate: None,
            conditional_coverage: None,
            unknown_axes: vec![
                CausalUnknownAxisState {
                    axis: CausalOutcomeAxis::Determination,
                    reason: CausalUnknownReason::NotAnalyzed,
                    resume_checkpoint: None,
                },
                CausalUnknownAxisState {
                    axis: CausalOutcomeAxis::StructuralRank,
                    reason: CausalUnknownReason::NotAnalyzed,
                    resume_checkpoint: None,
                },
            ],
            evidence: CausalReceiptEvidence::Unverified,
        }
        .admit_against(&graph, cx)
    })
    .expect("unknown-axis mode receipt admits honestly");
    let unknown_outcome = conditional_outcome(&unknown).expect("unknown mode child");
    assert_eq!(
        conditional_coverage_binding(
            &graph,
            core::slice::from_ref(&unknown_outcome),
            certificate,
            checker,
        ),
        Err(ConditionalCoverageBindingError::NonConcreteChild { outcome_index: 0 })
    );
}

#[test]
fn g0_explicit_cartesian_coverage_refuses_domains_beyond_public_envelope() {
    let (machine, owner, clock) = minimal_machine();
    let graph = cartesian_hybrid_graph(&machine, &owner, &clock, 2_049);
    let child = mode_cell_receipt(&graph, cartesian_mode_assignment(&graph, 0));
    let outcome = conditional_outcome(&child).expect("large-domain mode child");
    assert_eq!(
        conditional_coverage_binding(
            &graph,
            core::slice::from_ref(&outcome),
            cref!(ConditionalCoverageRef, "test/large-domain-coverage", 110),
            cref!(CausalCheckerRef, "test/large-domain-checker", 111),
        ),
        Err(ConditionalCoverageBindingError::ExplicitDomainTooLarge {
            required_outcomes: 4_098,
            max_outcomes: MAX_CAUSAL_CONDITIONAL_OUTCOMES,
            required_selections: 8_196,
            max_selections: MAX_CAUSAL_CONDITIONAL_SELECTIONS,
        })
    );
}

#[test]
fn g0_hybrid_summary_rejects_concrete_axes_contradicting_children() {
    let (machine, owner, clock) = minimal_machine();
    let graph = hybrid_graph(&machine, &owner, &clock);
    let child_a = mode_cell_receipt(&graph, mode_assignment(&graph, 0));
    let child_b = mode_cell_receipt(&graph, mode_assignment(&graph, 1));
    let outcomes = vec![
        conditional_outcome(&child_a).expect("mode child"),
        conditional_outcome(&child_b).expect("mode child"),
    ];
    let coverage = conditional_coverage_binding(
        &graph,
        &outcomes,
        cref!(ConditionalCoverageRef, "test/hybrid-coverage", 60),
        cref!(CausalCheckerRef, "test/hybrid-checker", 61),
    )
    .expect("coverage binding");
    let false_summary = CausalizationReceiptDraft {
        structure: graph.structure_identity_receipt(),
        artifact: graph.artifact_identity_receipt(),
        analysis: analysis(),
        domain: CausalReceiptDomain::HybridSummary,
        determination: DeterminationClass::WellDetermined,
        structural_rank: StructuralRankState::FullRelativeToMinSide,
        conditionality: Conditionality::Conditional,
        matching: Vec::new(),
        unmatched_equations: Vec::new(),
        unmatched_variables: Vec::new(),
        conditional_outcomes: outcomes,
        maximum_matching_certificate: None,
        conditional_coverage: Some(coverage),
        unknown_axes: Vec::new(),
        evidence: CausalReceiptEvidence::CheckerReferenced(cref!(
            CausalCheckerRef,
            "test/hybrid-checker",
            61
        )),
    };
    let refusal = with_cx(|cx| false_summary.admit_against(&graph, cx))
        .expect_err("contradictory concrete summary axes refuse");
    assert_receipt_rule(&refusal, CausalReceiptRule::OutcomeAxisMismatch);
}

#[test]
#[allow(clippy::too_many_lines)]
fn g0_uniform_theorem_axes_obey_exact_bipartition_semantics() {
    let (machine, owner, clock) = minimal_machine();
    let graph = hybrid_graph(&machine, &owner, &clock);
    let certificate = cref!(ConditionalCoverageRef, "test/uniform-axis-coverage", 64);
    let checker = cref!(CausalCheckerRef, "test/uniform-axis-checker", 65);

    for (determination, rank) in [
        (
            DeterminationClass::WellDetermined,
            StructuralRankState::Deficient,
        ),
        (
            DeterminationClass::Mixed,
            StructuralRankState::FullRelativeToMinSide,
        ),
        (
            DeterminationClass::WellDetermined,
            StructuralRankState::NotApplicable,
        ),
        (
            DeterminationClass::UnderDetermined,
            StructuralRankState::Deficient,
        ),
        (
            DeterminationClass::OverDetermined,
            StructuralRankState::Deficient,
        ),
        (
            DeterminationClass::Mixed,
            StructuralRankState::NotApplicable,
        ),
        (
            DeterminationClass::EmptyProjection,
            StructuralRankState::FullRelativeToMinSide,
        ),
        (
            DeterminationClass::EmptyProjection,
            StructuralRankState::Deficient,
        ),
    ] {
        let refusal = ConditionalCoverageBinding::for_uniform_theorem(
            &graph,
            determination,
            rank,
            certificate.clone(),
            checker.clone(),
        );
        assert_eq!(
            refusal,
            Err(ConditionalCoverageBindingError::IncompatibleUniformClaim),
            "incompatible uniform axes unexpectedly admitted: determination={determination:?}, rank={rank:?}"
        );
    }

    for (determination, rank) in [
        (
            DeterminationClass::Unknown,
            StructuralRankState::FullRelativeToMinSide,
        ),
        (
            DeterminationClass::WellDetermined,
            StructuralRankState::Unknown,
        ),
        (DeterminationClass::Unknown, StructuralRankState::Unknown),
    ] {
        let refusal = ConditionalCoverageBinding::for_uniform_theorem(
            &graph,
            determination,
            rank,
            certificate.clone(),
            checker.clone(),
        );
        assert_eq!(
            refusal,
            Err(ConditionalCoverageBindingError::NonConcreteUniformClaim),
            "non-concrete uniform axes returned the wrong refusal: determination={determination:?}, rank={rank:?}"
        );
    }

    for (determination, rank) in [
        (
            DeterminationClass::EmptyProjection,
            StructuralRankState::NotApplicable,
        ),
        (
            DeterminationClass::WellDetermined,
            StructuralRankState::FullRelativeToMinSide,
        ),
        (
            DeterminationClass::UnderDetermined,
            StructuralRankState::FullRelativeToMinSide,
        ),
        (
            DeterminationClass::OverDetermined,
            StructuralRankState::FullRelativeToMinSide,
        ),
        (DeterminationClass::Mixed, StructuralRankState::Deficient),
        (
            DeterminationClass::UnderDetermined,
            StructuralRankState::NotApplicable,
        ),
        (
            DeterminationClass::OverDetermined,
            StructuralRankState::NotApplicable,
        ),
    ] {
        ConditionalCoverageBinding::for_uniform_theorem(
            &graph,
            determination,
            rank,
            certificate.clone(),
            checker.clone(),
        )
        .unwrap_or_else(|error| {
            panic!(
                "compatible uniform axes refused: determination={determination:?}, rank={rank:?}, error={error}"
            )
        });
    }

    let mut uniform_draft = hybrid_draft(&owner, &clock);
    for incidence in &mut uniform_draft.incidences {
        *incidence = incidence_spec(
            incidence.equation.clone(),
            incidence.variable.clone(),
            incidence.derivative_order,
            SolveParticipation::Unknown,
            incidence.coefficient_dimensions,
            incidence.term.clone(),
            incidence.operator.clone(),
            incidence.clock_relation.clone(),
            incidence.activation.clone(),
        )
        .expect("uniform unknown occurrence identity");
    }
    let uniform_graph = with_cx(|cx| uniform_draft.admit_against(&machine, cx))
        .expect("two-branch uniformly unknown graph admits");
    let coverage = ConditionalCoverageBinding::for_uniform_theorem(
        &uniform_graph,
        DeterminationClass::WellDetermined,
        StructuralRankState::FullRelativeToMinSide,
        certificate,
        checker.clone(),
    )
    .expect("compatible uniform theorem binding");
    let receipt = with_cx(|cx| {
        CausalizationReceiptDraft {
            structure: uniform_graph.structure_identity_receipt(),
            artifact: uniform_graph.artifact_identity_receipt(),
            analysis: analysis(),
            domain: CausalReceiptDomain::HybridSummary,
            determination: DeterminationClass::WellDetermined,
            structural_rank: StructuralRankState::FullRelativeToMinSide,
            conditionality: Conditionality::Unconditional,
            matching: Vec::new(),
            unmatched_equations: Vec::new(),
            unmatched_variables: Vec::new(),
            conditional_outcomes: Vec::new(),
            maximum_matching_certificate: None,
            conditional_coverage: Some(coverage),
            unknown_axes: Vec::new(),
            evidence: CausalReceiptEvidence::CheckerReferenced(checker),
        }
        .admit_against(&uniform_graph, cx)
    })
    .expect("schema-compatible uniform theorem commitment admits");
    assert_eq!(receipt.determination(), DeterminationClass::WellDetermined);
    assert_eq!(
        receipt.structural_rank(),
        StructuralRankState::FullRelativeToMinSide
    );
}

#[test]
fn g4_oversized_conditional_children_refuse_before_nested_scans_or_clones() {
    let (machine, owner, clock) = minimal_machine();
    let graph = hybrid_graph(&machine, &owner, &clock);
    let selection = mode_assignment(&graph, 0)
        .pop()
        .expect("mode fixture has one selection");
    assert_eq!(
        maximum_matching_binding(
            &graph,
            CausalReceiptDomain::ModeCell {
                assignment: vec![selection; MAX_CAUSAL_CONDITIONS + 1],
            },
            &[],
            cref!(
                MaximumMatchingCertificateRef,
                "test/oversized-mode-maximum",
                68
            ),
            cref!(CausalCheckerRef, "test/oversized-mode-checker", 69),
        ),
        Err(MaximumMatchingBindingError::InvalidDomain(
            MaximumMatchingDomainError::AssignmentCardinality {
                submitted: MAX_CAUSAL_CONDITIONS + 1,
                expected: graph.conditions().len(),
                max: MAX_CAUSAL_CONDITIONS,
            }
        )),
        "maximum binding must cap a mode assignment before sorting it"
    );
    let child = mode_cell_receipt(&graph, mode_assignment(&graph, 0));
    let outcome = conditional_outcome(&child).expect("mode child");
    let oversized = vec![outcome; MAX_CAUSAL_CONDITIONAL_OUTCOMES + 1];
    assert_eq!(
        conditional_coverage_binding(
            &graph,
            &oversized,
            cref!(ConditionalCoverageRef, "test/oversized-coverage", 66),
            cref!(CausalCheckerRef, "test/oversized-checker", 67),
        ),
        Err(ConditionalCoverageBindingError::OutcomeSetLimit {
            submitted: MAX_CAUSAL_CONDITIONAL_OUTCOMES + 1,
            max: MAX_CAUSAL_CONDITIONAL_OUTCOMES,
        }),
        "coverage binding must apply its outer cap before cloning assignments"
    );
    let decision = with_cx(|cx| {
        CausalizationReceiptDraft {
            structure: graph.structure_identity_receipt(),
            artifact: graph.artifact_identity_receipt(),
            analysis: analysis(),
            domain: CausalReceiptDomain::HybridSummary,
            determination: DeterminationClass::Unknown,
            structural_rank: StructuralRankState::Unknown,
            conditionality: Conditionality::Unknown,
            matching: Vec::new(),
            unmatched_equations: Vec::new(),
            unmatched_variables: Vec::new(),
            conditional_outcomes: oversized,
            maximum_matching_certificate: None,
            conditional_coverage: None,
            unknown_axes: Vec::new(),
            evidence: CausalReceiptEvidence::Unverified,
        }
        .admit_with_decision(&graph, cx)
    });
    assert!(!decision.submitted_counts().complete);
    assert_eq!(
        decision.submitted_counts().conditional_outcomes,
        MAX_CAUSAL_CONDITIONAL_OUTCOMES + 1
    );
    let refusal = decision
        .into_result()
        .expect_err("oversized child set refuses before receipt canonicalization");
    assert_receipt_rules_exact(&refusal, &[CausalReceiptRule::ResourceLimit]);
}

#[test]
fn g0_duplicate_nominal_ids_refuse_before_caller_order_can_select_a_payload() {
    let (machine, owner, clock) = minimal_machine();
    let mut draft = minimal_causal_draft(&owner, &clock);
    let mut duplicate = draft.equations[0].clone();
    duplicate.role = EquationRole::Balance;
    draft.equations.push(duplicate);
    let refusal = with_cx(|cx| draft.admit_against(&machine, cx))
        .expect_err("duplicate nominal equation refuses");
    assert_graph_rules_exact(&refusal, &[CausalGraphRule::DuplicateEquation]);
}

#[test]
fn g0_guard_equations_are_condition_only_and_must_define_a_condition() {
    let (machine, owner, clock) = minimal_machine();
    let mut draft = minimal_causal_draft(&owner, &clock);
    draft.equations[0].role = EquationRole::Guard;
    let refusal =
        with_cx(|cx| draft.admit_against(&machine, cx)).expect_err("orphan matching guard refuses");
    assert_graph_rule(&refusal, CausalGraphRule::InvalidActivationCondition);
}

#[test]
fn g0_foreign_ownership_derivative_units_and_clock_contracts_refuse_exactly() {
    let (machine, owner, clock) = minimal_machine();

    let mut foreign_owner = minimal_causal_draft(&owner, &clock);
    foreign_owner.equations[0].owner = CausalOwner::Subsystem(
        SubsystemId::new("subsystem/not-in-machine").expect("valid foreign subsystem ID"),
    );
    let refusal = with_cx(|cx| foreign_owner.admit_against(&machine, cx))
        .expect_err("foreign node owner must refuse");
    assert_graph_rule(&refusal, CausalGraphRule::UnknownOwner);

    let mut derivative = minimal_causal_draft(&owner, &clock);
    let equation = derivative.equations[0].clone();
    let variable = derivative.variables[0].clone();
    derivative.incidences[0] = incidence_spec(
        equation.id.clone(),
        variable.id.clone(),
        MAX_CAUSAL_DERIVATIVE_ORDER + 1,
        SolveParticipation::Unknown,
        Dims::NONE,
        equation.residual.clone(),
        None,
        IncidenceClockRelation::SameClock,
        ActivationDomain::Always,
    )
    .expect("out-of-policy derivative remains canonically representable");
    let refusal = with_cx(|cx| derivative.admit_against(&machine, cx))
        .expect_err("v1 derivative-order overflow must refuse");
    assert_graph_rule(&refusal, CausalGraphRule::DerivativeOrderLimit);

    let mut wrong_units = minimal_causal_draft(&owner, &clock);
    let equation = wrong_units.equations[0].clone();
    let variable = wrong_units.variables[0].clone();
    wrong_units.incidences[0] = incidence_spec(
        equation.id.clone(),
        variable.id.clone(),
        0,
        SolveParticipation::Unknown,
        Dims([1, 0, 0, 0, 0, 0]),
        equation.residual.clone(),
        Some(cref!(
            IncidenceOperatorRef,
            "test/dimension-changing-operator",
            70
        )),
        IncidenceClockRelation::SameClock,
        ActivationDomain::Always,
    )
    .expect("dimensionally inconsistent row remains canonically representable");
    let refusal = with_cx(|cx| wrong_units.admit_against(&machine, cx))
        .expect_err("dimensionally open incidence must refuse");
    assert_graph_rule(&refusal, CausalGraphRule::IncidenceUnitMismatch);

    let mut wrong_clock = minimal_causal_draft(&owner, &clock);
    let equation = wrong_clock.equations[0].clone();
    let variable = wrong_clock.variables[0].clone();
    let bridge = cref!(ClockBridgeRef, "test/identity-clock-bridge", 71);
    wrong_clock.incidences[0] = incidence_spec(
        equation.id.clone(),
        variable.id.clone(),
        0,
        SolveParticipation::Unknown,
        Dims::NONE,
        equation.residual.clone(),
        None,
        IncidenceClockRelation::AuditedBridge {
            source: clock.clone(),
            target: clock.clone(),
            bridge: bridge.clone(),
            audit: cref!(ClockBridgeAuditRef, "test/identity-clock-audit", 72),
            audited_bridge: bridge,
        },
        ActivationDomain::Always,
    )
    .expect("invalid same-clock bridge remains canonically representable");
    let refusal = with_cx(|cx| wrong_clock.admit_against(&machine, cx))
        .expect_err("bridge declaration cannot disguise a same-clock edge");
    assert_graph_rule(&refusal, CausalGraphRule::IncidenceClockMismatch);
}

fn indexed_test_digest(domain: u8, index: usize) -> [u8; 32] {
    let mut digest = [0_u8; 32];
    digest[..8].copy_from_slice(
        &u64::try_from(index)
            .expect("test index fits u64")
            .to_le_bytes(),
    );
    digest[31] = domain;
    digest
}

fn maximal_test_namespace(role: &str, tag: char, index: usize) -> String {
    let stem = format!("test/{role}/{tag}{index:08x}/");
    let namespace = format!(
        "{stem}{}",
        tag.to_string()
            .repeat(MAX_MACHINE_ENTITY_KEY_BYTES - stem.len())
    );
    assert_eq!(
        namespace.len(),
        MAX_MACHINE_ENTITY_KEY_BYTES,
        "role={role}; index={index}"
    );
    namespace
}

fn maximum_legal_activation_envelope() -> ActivationDomain {
    assert_eq!(
        MAX_CAUSAL_ACTIVATION_SELECTIONS % MAX_CAUSAL_CUBES_PER_ACTIVATION,
        0,
        "aggregate selection cap must decompose across the maximum cube count"
    );
    let conditions_per_cube = MAX_CAUSAL_ACTIVATION_SELECTIONS / MAX_CAUSAL_CUBES_PER_ACTIVATION;
    let conditions: Vec<_> = (0..conditions_per_cube)
        .map(|index| {
            ActivationConditionRef::new(
                maximal_test_namespace("condition", 'c', index),
                nz(1),
                indexed_test_digest(73, index),
            )
            .expect("maximum-length unique condition reference")
        })
        .collect();
    let branches: Vec<_> = (0..MAX_CAUSAL_CUBES_PER_ACTIVATION)
        .map(|index| {
            ActivationBranchRef::new(
                maximal_test_namespace("branch", 'b', index),
                nz(1),
                indexed_test_digest(74, index),
            )
            .expect("maximum-length unique branch reference")
        })
        .collect();
    let fixed_branch = ActivationBranchRef::new(
        maximal_test_namespace("branch", 'b', MAX_CAUSAL_CUBES_PER_ACTIVATION),
        nz(1),
        indexed_test_digest(74, MAX_CAUSAL_CUBES_PER_ACTIVATION),
    )
    .expect("maximum-length fixed branch reference");
    ActivationDomain::Conditional {
        cubes: branches
            .into_iter()
            .map(|branch| ActivationCube {
                selections: conditions
                    .iter()
                    .cloned()
                    .enumerate()
                    .map(|(condition_index, condition)| ConditionBranchSelection {
                        condition,
                        branch: if condition_index == 0 {
                            branch.clone()
                        } else {
                            fixed_branch.clone()
                        },
                    })
                    .collect(),
            })
            .collect(),
    }
}

#[test]
#[allow(clippy::too_many_lines)]
fn g4_oversized_single_activation_refuses_before_canonical_sorting() {
    let (machine, owner, clock) = minimal_machine();
    let mut draft = minimal_causal_draft(&owner, &clock);

    let at_cap = incidence_spec(
        draft.equations[0].id.clone(),
        draft.variables[0].id.clone(),
        0,
        SolveParticipation::Unknown,
        Dims::NONE,
        draft.equations[0].residual.clone(),
        None,
        IncidenceClockRelation::SameClock,
        maximum_legal_activation_envelope(),
    )
    .expect("standalone incidence identity must realize the maximum activation sub-envelope");
    assert_eq!(at_cap.id.identity_receipt().field_count(), 1);
    let canonical_ref_bytes = 8 + MAX_MACHINE_ENTITY_KEY_BYTES + 8 + 32;
    let activation_preimage_bytes = 1
        + 8
        + MAX_CAUSAL_CUBES_PER_ACTIVATION * 8
        + MAX_CAUSAL_ACTIVATION_SELECTIONS * 2 * canonical_ref_bytes;
    assert!(
        at_cap.id.identity_receipt().canonical_bytes()
            > u64::try_from(activation_preimage_bytes).expect("activation envelope fits u64"),
        "canonical receipt must retain all max-length activation references plus fixed incidence/framing bytes"
    );
    drop(at_cap);

    let ActivationDomain::Conditional {
        cubes: mut aggregate_oversized_cubes,
    } = maximum_legal_activation_envelope()
    else {
        unreachable!("maximum envelope is conditional")
    };
    aggregate_oversized_cubes[0]
        .selections
        .push(ConditionBranchSelection {
            condition: ActivationConditionRef::new(
                maximal_test_namespace("condition", 'c', 64),
                nz(1),
                indexed_test_digest(73, 64),
            )
            .expect("one-over-cap condition reference"),
            branch: ActivationBranchRef::new(
                maximal_test_namespace("branch", 'b', MAX_CAUSAL_CUBES_PER_ACTIVATION),
                nz(1),
                indexed_test_digest(74, MAX_CAUSAL_CUBES_PER_ACTIVATION),
            )
            .expect("one-over-cap branch reference"),
        });
    let aggregate_error = incidence_spec(
        draft.equations[0].id.clone(),
        draft.variables[0].id.clone(),
        0,
        SolveParticipation::Unknown,
        Dims::NONE,
        draft.equations[0].residual.clone(),
        None,
        IncidenceClockRelation::SameClock,
        ActivationDomain::Conditional {
            cubes: aggregate_oversized_cubes,
        },
    )
    .expect_err("aggregate selection cap must refuse before sorting or identity encoding");
    assert_eq!(
        aggregate_error,
        CanonicalError::LimitExceeded {
            kind: LimitKind::CollectionItems,
            requested: u64::try_from(MAX_CAUSAL_ACTIVATION_SELECTIONS + 1)
                .expect("test cap fits u64"),
            limit: u64::try_from(MAX_CAUSAL_ACTIVATION_SELECTIONS).expect("test cap fits u64"),
        }
    );

    let condition = cref!(ActivationConditionRef, "test/resource-condition", 73);
    let branch = cref!(ActivationBranchRef, "test/resource-branch", 74);
    let cube = ActivationCube {
        selections: vec![ConditionBranchSelection { condition, branch }],
    };
    let oversized_activation = ActivationDomain::Conditional {
        cubes: vec![cube; MAX_CAUSAL_CUBES_PER_ACTIVATION + 1],
    };
    draft.variables[0].activation = oversized_activation.clone();
    let constructor_error = incidence_spec(
        draft.equations[0].id.clone(),
        draft.variables[0].id.clone(),
        0,
        SolveParticipation::Unknown,
        Dims::NONE,
        draft.equations[0].residual.clone(),
        None,
        IncidenceClockRelation::SameClock,
        oversized_activation,
    )
    .expect_err("incidence constructor must refuse before sorting an oversized DNF");
    assert_eq!(
        constructor_error,
        CanonicalError::LimitExceeded {
            kind: LimitKind::CollectionItems,
            requested: u64::try_from(MAX_CAUSAL_CUBES_PER_ACTIVATION + 1)
                .expect("test cap fits u64"),
            limit: u64::try_from(MAX_CAUSAL_CUBES_PER_ACTIVATION).expect("test cap fits u64"),
        }
    );
    let decision = with_cx(|cx| draft.admit_with_decision(&machine, cx));
    assert!(
        !decision.submitted_counts().complete,
        "first exceeded nested cap must short-circuit remaining telemetry"
    );
    assert_eq!(
        decision.submitted_counts().max_activation_cubes_per_row,
        MAX_CAUSAL_CUBES_PER_ACTIVATION + 1
    );
    let refusal = decision
        .into_result()
        .expect_err("one oversized DNF row must fail before sorting");
    assert_graph_rules_exact(&refusal, &[CausalGraphRule::ResourceLimit]);
}

#[test]
fn g4_single_cube_selection_boundary_is_exact_and_pre_sort() {
    let (machine, owner, clock) = minimal_machine();
    let draft = minimal_causal_draft(&owner, &clock);
    let branch = ActivationBranchRef::new(
        maximal_test_namespace("branch", 'q', 0),
        nz(1),
        indexed_test_digest(90, 0),
    )
    .expect("maximum-length branch reference");
    let mut selections: Vec<_> = (0..=MAX_CAUSAL_SELECTIONS_PER_CUBE)
        .map(|index| ConditionBranchSelection {
            condition: ActivationConditionRef::new(
                maximal_test_namespace("condition", 'p', index),
                nz(1),
                indexed_test_digest(91, index),
            )
            .expect("maximum-length unique condition reference"),
            branch: branch.clone(),
        })
        .collect();
    let one_over = selections
        .pop()
        .expect("one-over fixture contains a final selection");
    let at_cap_activation = ActivationDomain::Conditional {
        cubes: vec![ActivationCube {
            selections: selections.clone(),
        }],
    };
    incidence_spec(
        draft.equations[0].id.clone(),
        draft.variables[0].id.clone(),
        0,
        SolveParticipation::Unknown,
        Dims::NONE,
        draft.equations[0].residual.clone(),
        None,
        IncidenceClockRelation::SameClock,
        at_cap_activation,
    )
    .expect("one cube at the exact public selection cap must be representable");

    selections.push(one_over);
    let oversized_activation = ActivationDomain::Conditional {
        cubes: vec![ActivationCube { selections }],
    };
    let constructor_error = incidence_spec(
        draft.equations[0].id.clone(),
        draft.variables[0].id.clone(),
        0,
        SolveParticipation::Unknown,
        Dims::NONE,
        draft.equations[0].residual.clone(),
        None,
        IncidenceClockRelation::SameClock,
        oversized_activation.clone(),
    )
    .expect_err("one cube above its selection cap must refuse before sorting");
    assert_eq!(
        constructor_error,
        CanonicalError::LimitExceeded {
            kind: LimitKind::CollectionItems,
            requested: u64::try_from(MAX_CAUSAL_SELECTIONS_PER_CUBE + 1)
                .expect("test cap fits u64"),
            limit: u64::try_from(MAX_CAUSAL_SELECTIONS_PER_CUBE).expect("test cap fits u64"),
        }
    );

    let mut hostile_draft = draft;
    hostile_draft.variables[0].activation = oversized_activation;
    let decision = with_cx(|cx| hostile_draft.admit_with_decision(&machine, cx));
    assert!(
        !decision.submitted_counts().complete,
        "per-cube overflow must short-circuit telemetry before graph canonicalization"
    );
    assert_eq!(
        decision
            .submitted_counts()
            .max_activation_selections_per_cube,
        MAX_CAUSAL_SELECTIONS_PER_CUBE + 1
    );
    let refusal = decision
        .into_result()
        .expect_err("per-cube overflow must refuse the graph");
    assert_graph_rules_exact(&refusal, &[CausalGraphRule::ResourceLimit]);
}

#[test]
fn g4_oversized_outer_graph_collection_refuses_before_nested_telemetry() {
    let (machine, owner, clock) = minimal_machine();
    let mut template_draft = hybrid_draft(&owner, &clock);
    let condition = template_draft
        .conditions
        .pop()
        .expect("hybrid fixture has one condition");
    let mut draft = minimal_causal_draft(&owner, &clock);
    draft.conditions = vec![condition; MAX_CAUSAL_CONDITIONS + 1];
    let decision = with_cx(|cx| draft.admit_with_decision(&machine, cx));
    assert!(
        !decision.submitted_counts().complete,
        "outer cap must refuse before nested telemetry traversal"
    );
    assert_eq!(
        decision.submitted_counts().conditions,
        MAX_CAUSAL_CONDITIONS + 1
    );
    let refusal = decision
        .into_result()
        .expect_err("oversized condition collection must refuse immediately");
    assert_graph_rules_exact(&refusal, &[CausalGraphRule::ResourceLimit]);
}

#[test]
fn g0_mode_dependent_totality_requires_unknown_or_known_resolution() {
    let (machine, owner, clock) = minimal_machine();
    let mut draft = hybrid_draft(&owner, &clock);
    for incidence in &mut draft.incidences {
        *incidence = incidence_spec(
            incidence.equation.clone(),
            incidence.variable.clone(),
            incidence.derivative_order,
            SolveParticipation::ConditionOnly,
            incidence.coefficient_dimensions,
            incidence.term.clone(),
            incidence.operator.clone(),
            incidence.clock_relation.clone(),
            incidence.activation.clone(),
        )
        .expect("condition-only occurrence identity");
    }
    let refusal = with_cx(|cx| draft.admit_against(&machine, cx))
        .expect_err("auxiliary reads cannot discharge solve-participation totality");
    assert_graph_rule(&refusal, CausalGraphRule::DerivativeParticipationMismatch);
}

#[test]
fn g0_finite_domain_dnf_implication_proves_consensus_and_finds_counterexample() {
    let (machine, owner, clock) = minimal_machine();
    let mut draft = minimal_causal_draft(&owner, &clock);
    let mut parameter = variable(&owner, &clock, "test/variable/dnf-parameter", 75);
    parameter.role = VariableRole::Parameter;
    parameter.solve_participation = SolveParticipation::ConditionOnly;

    let condition_a = cref!(ActivationConditionRef, "test/dnf/condition-a", 76);
    let a_on = cref!(ActivationBranchRef, "test/dnf/a-on", 77);
    let a_off = cref!(ActivationBranchRef, "test/dnf/a-off", 78);
    let condition_b = cref!(ActivationConditionRef, "test/dnf/condition-b", 79);
    let b_zero = cref!(ActivationBranchRef, "test/dnf/b-zero", 80);
    let b_one = cref!(ActivationBranchRef, "test/dnf/b-one", 81);
    let source_a = cref!(SourceArtifactRef, "test/dnf/source-a", 82);
    let source_b = cref!(SourceArtifactRef, "test/dnf/source-b", 83);
    draft.conditions = vec![
        ActivationConditionSpec {
            condition: condition_a.clone(),
            source: ActivationConditionSource::AuditedPredicate(AuditedEscapeHatch {
                source: source_a.clone(),
                audit: cref!(EscapeHatchAuditRef, "test/dnf/audit-a", 84),
                audited_source: source_a,
            }),
            branches: vec![a_on.clone(), a_off],
            dependencies: vec![parameter.id.clone()],
        },
        ActivationConditionSpec {
            condition: condition_b.clone(),
            source: ActivationConditionSource::AuditedPredicate(AuditedEscapeHatch {
                source: source_b.clone(),
                audit: cref!(EscapeHatchAuditRef, "test/dnf/audit-b", 85),
                audited_source: source_b,
            }),
            branches: vec![b_zero.clone(), b_one.clone()],
            dependencies: vec![parameter.id.clone()],
        },
    ];
    draft.variables.push(parameter);

    let consensus_cube = |branch_b: &ActivationBranchRef| ActivationCube {
        selections: vec![
            ConditionBranchSelection {
                condition: condition_a.clone(),
                branch: a_on.clone(),
            },
            ConditionBranchSelection {
                condition: condition_b.clone(),
                branch: branch_b.clone(),
            },
        ],
    };
    draft.variables[0].activation = ActivationDomain::Conditional {
        cubes: vec![consensus_cube(&b_zero), consensus_cube(&b_one)],
    };
    let equation = draft.equations[0].clone();
    let causal_variable = draft.variables[0].clone();
    draft.incidences[0] = incidence_spec(
        equation.id.clone(),
        causal_variable.id.clone(),
        0,
        SolveParticipation::Unknown,
        Dims::NONE,
        equation.residual.clone(),
        None,
        IncidenceClockRelation::SameClock,
        branch(&condition_a, &a_on),
    )
    .expect("consensus antecedent incidence");

    let mut missing_branch = draft.clone();
    let ActivationDomain::Conditional { cubes } = &mut missing_branch.variables[0].activation
    else {
        unreachable!("fixture activation is conditional");
    };
    assert!(cubes.pop().is_some(), "fixture removes one B-domain branch");
    let refusal = with_cx(|cx| missing_branch.admit_against(&machine, cx))
        .expect_err("A does not imply A-and-B0 when B1 remains in the finite domain");
    assert_graph_rule(&refusal, CausalGraphRule::ActivationMismatch);

    with_cx(|cx| draft.admit_against(&machine, cx))
        .expect("(A and B0) or (A and B1) exactly covers A over the declared B domain");
}

#[test]
fn g0_schema_migration_binds_complete_same_family_native_receipt() {
    let (_machine, graph) = admit_minimal();
    let predecessor = HistoricalCausalIdentityReceipt::new(
        CausalMigrationArtifactKind::Structure,
        0,
        [1; 32],
        [2; 32],
        [3; 32],
        0,
        0,
        0,
    )
    .expect("legacy predecessor receipt");
    let migration = with_cx(|cx| {
        CausalSchemaMigrationDraft::for_structure(
            predecessor,
            &graph,
            cref!(CausalMigrationRef, "test/causal-migration", 62),
        )
        .admit(cx)
    })
    .expect("same-family migration admits");
    let target = graph.structure_identity_receipt();
    assert_eq!(migration.target_identity(), *target.id().as_bytes());
    assert_eq!(
        migration.target_canonical_preimage(),
        *target.canonical_preimage().as_bytes()
    );
    assert_eq!(
        migration.target_schema_identity(),
        *target.schema_id().as_bytes()
    );
    assert_eq!(migration.target_canonical_bytes(), target.canonical_bytes());
    assert_eq!(migration.target_field_count(), target.field_count());
    assert_eq!(
        migration.target_collection_items(),
        target.collection_items()
    );
}

#[test]
fn g0_schema_migration_refuses_incomplete_history_and_family_or_version_substitution() {
    assert_eq!(
        HistoricalCausalIdentityReceipt::new(
            CausalMigrationArtifactKind::Structure,
            0,
            [0; 32],
            [2; 32],
            [3; 32],
            0,
            0,
            0,
        ),
        Err(HistoricalReceiptError::ZeroDigest)
    );
    assert_eq!(
        HistoricalCausalIdentityReceipt::new(
            CausalMigrationArtifactKind::Structure,
            1,
            [1; 32],
            [2; 32],
            [3; 32],
            0,
            0,
            0,
        ),
        Err(HistoricalReceiptError::IncompleteCanonicalMetadata)
    );

    let (_machine, graph) = admit_minimal();
    let wrong_family = HistoricalCausalIdentityReceipt::new(
        CausalMigrationArtifactKind::GraphArtifact,
        0,
        [4; 32],
        [5; 32],
        [6; 32],
        0,
        0,
        0,
    )
    .expect("legacy wrong-family receipt");
    assert_eq!(
        with_cx(|cx| {
            CausalSchemaMigrationDraft::for_structure(
                wrong_family,
                &graph,
                cref!(CausalMigrationRef, "test/wrong-family-migration", 115),
            )
            .admit(cx)
        }),
        Err(CausalMigrationError::ArtifactKindMismatch)
    );
    let same_version = HistoricalCausalIdentityReceipt::new(
        CausalMigrationArtifactKind::Structure,
        CAUSAL_STRUCTURE_IDENTITY_SCHEMA_VERSION_V1,
        [7; 32],
        [8; 32],
        [9; 32],
        1,
        1,
        0,
    )
    .expect("canonical same-version predecessor");
    assert_eq!(
        with_cx(|cx| {
            CausalSchemaMigrationDraft::for_structure(
                same_version,
                &graph,
                cref!(CausalMigrationRef, "test/not-older-migration", 116),
            )
            .admit(cx)
        }),
        Err(CausalMigrationError::PredecessorNotOlder {
            predecessor: CAUSAL_STRUCTURE_IDENTITY_SCHEMA_VERSION_V1,
            target: CAUSAL_STRUCTURE_IDENTITY_SCHEMA_VERSION_V1,
        })
    );
}

#[test]
fn g0_schema_migration_binds_graph_artifact_and_causalization_receipt_targets() {
    let (_machine, graph) = admit_minimal();
    let graph_predecessor = HistoricalCausalIdentityReceipt::new(
        CausalMigrationArtifactKind::GraphArtifact,
        0,
        [10; 32],
        [11; 32],
        [12; 32],
        0,
        0,
        0,
    )
    .expect("legacy graph-artifact predecessor");
    let graph_migration = with_cx(|cx| {
        CausalSchemaMigrationDraft::for_graph_artifact(
            graph_predecessor,
            &graph,
            cref!(CausalMigrationRef, "test/graph-artifact-migration", 117),
        )
        .admit(cx)
    })
    .expect("graph-artifact migration admits");
    let graph_target = graph.artifact_identity_receipt();
    assert_eq!(
        graph_migration.kind(),
        CausalMigrationArtifactKind::GraphArtifact
    );
    assert_eq!(
        graph_migration.target_identity(),
        *graph_target.id().as_bytes()
    );
    assert_eq!(
        graph_migration.target_canonical_preimage(),
        *graph_target.canonical_preimage().as_bytes()
    );
    assert_eq!(
        graph_migration.target_schema_identity(),
        *graph_target.schema_id().as_bytes()
    );
    assert_eq!(
        graph_migration.target_canonical_bytes(),
        graph_target.canonical_bytes()
    );
    assert_eq!(
        graph_migration.target_field_count(),
        graph_target.field_count()
    );
    assert_eq!(
        graph_migration.target_collection_items(),
        graph_target.collection_items()
    );

    let causalization = with_cx(|cx| complete_receipt(&graph).admit_against(&graph, cx))
        .expect("causalization target admits");
    let receipt_predecessor = HistoricalCausalIdentityReceipt::new(
        CausalMigrationArtifactKind::CausalizationReceipt,
        0,
        [13; 32],
        [14; 32],
        [15; 32],
        0,
        0,
        0,
    )
    .expect("legacy causalization predecessor");
    let receipt_migration = with_cx(|cx| {
        CausalSchemaMigrationDraft::for_causalization_receipt(
            receipt_predecessor,
            &causalization,
            cref!(CausalMigrationRef, "test/causalization-migration", 118),
        )
        .admit(cx)
    })
    .expect("causalization migration admits");
    let receipt_target = causalization.identity_receipt();
    assert_eq!(
        receipt_migration.kind(),
        CausalMigrationArtifactKind::CausalizationReceipt
    );
    assert_eq!(
        receipt_migration.target_identity(),
        *receipt_target.id().as_bytes()
    );
    assert_eq!(
        receipt_migration.target_canonical_preimage(),
        *receipt_target.canonical_preimage().as_bytes()
    );
    assert_eq!(
        receipt_migration.target_schema_identity(),
        *receipt_target.schema_id().as_bytes()
    );
    assert_eq!(
        receipt_migration.target_canonical_bytes(),
        receipt_target.canonical_bytes()
    );
    assert_eq!(
        receipt_migration.target_field_count(),
        receipt_target.field_count()
    );
    assert_eq!(
        receipt_migration.target_collection_items(),
        receipt_target.collection_items()
    );
}

#[test]
#[allow(clippy::too_many_lines)]
fn g0_theorem_binding_constructors_reject_uninhabitable_domains() {
    let (_, unconditional_graph) = admit_minimal();
    let unconditional_receipt = with_cx(|cx| {
        complete_receipt(&unconditional_graph).admit_against(&unconditional_graph, cx)
    })
    .expect("unconditional receipt admits");
    assert_eq!(
        conditional_outcome(&unconditional_receipt),
        Err(ConditionalOutcomeError::NotModeCell)
    );
    assert_eq!(
        maximum_matching_binding(
            &unconditional_graph,
            CausalReceiptDomain::ModeCell {
                assignment: Vec::new(),
            },
            &[],
            cref!(MaximumMatchingCertificateRef, "test/empty-mode-maximum", 80),
            cref!(CausalCheckerRef, "test/empty-mode-checker", 81),
        ),
        Err(MaximumMatchingBindingError::InvalidDomain(
            MaximumMatchingDomainError::ConditionFreeGraph
        ))
    );
    assert_eq!(
        ConditionalCoverageBinding::for_uniform_theorem(
            &unconditional_graph,
            DeterminationClass::WellDetermined,
            StructuralRankState::FullRelativeToMinSide,
            cref!(ConditionalCoverageRef, "test/empty-domain-coverage", 82),
            cref!(CausalCheckerRef, "test/empty-domain-checker", 83),
        ),
        Err(ConditionalCoverageBindingError::InvalidGraphDomain)
    );

    let (machine, owner, clock) = minimal_machine();
    let conditional_graph = hybrid_graph(&machine, &owner, &clock);
    assert_eq!(
        maximum_matching_binding(
            &conditional_graph,
            CausalReceiptDomain::HybridSummary,
            &[],
            cref!(MaximumMatchingCertificateRef, "test/hybrid-maximum", 113),
            cref!(CausalCheckerRef, "test/hybrid-maximum-checker", 114),
        ),
        Err(MaximumMatchingBindingError::HybridSummaryDomain)
    );
    assert_eq!(
        maximum_matching_binding(
            &conditional_graph,
            CausalReceiptDomain::UnconditionalGraph,
            &[],
            cref!(
                MaximumMatchingCertificateRef,
                "test/conditional-maximum",
                84
            ),
            cref!(CausalCheckerRef, "test/conditional-checker", 85),
        ),
        Err(MaximumMatchingBindingError::InvalidDomain(
            MaximumMatchingDomainError::ConditionalGraph
        ))
    );
    let mut wrong_condition = mode_assignment(&conditional_graph, 0);
    wrong_condition[0].condition = cref!(ActivationConditionRef, "test/foreign-condition", 90);
    assert_eq!(
        maximum_matching_binding(
            &conditional_graph,
            CausalReceiptDomain::ModeCell {
                assignment: wrong_condition,
            },
            &[],
            cref!(
                MaximumMatchingCertificateRef,
                "test/wrong-condition-maximum",
                91
            ),
            cref!(CausalCheckerRef, "test/wrong-condition-checker", 92),
        ),
        Err(MaximumMatchingBindingError::InvalidDomain(
            MaximumMatchingDomainError::InvalidSelection { index: 0 }
        ))
    );
    let mut wrong_branch = mode_assignment(&conditional_graph, 0);
    wrong_branch[0].branch = cref!(ActivationBranchRef, "test/foreign-branch", 93);
    assert_eq!(
        maximum_matching_binding(
            &conditional_graph,
            CausalReceiptDomain::ModeCell {
                assignment: wrong_branch,
            },
            &[],
            cref!(
                MaximumMatchingCertificateRef,
                "test/wrong-branch-maximum",
                94
            ),
            cref!(CausalCheckerRef, "test/wrong-branch-checker", 95),
        ),
        Err(MaximumMatchingBindingError::InvalidDomain(
            MaximumMatchingDomainError::InvalidSelection { index: 0 }
        ))
    );

    let mode_domain = CausalReceiptDomain::ModeCell {
        assignment: mode_assignment(&conditional_graph, 0),
    };
    let mode_receipt =
        mode_cell_receipt(&conditional_graph, mode_assignment(&conditional_graph, 0));
    let valid_pair = mode_receipt.matching()[0].clone();
    maximum_matching_binding(
        &conditional_graph,
        mode_domain.clone(),
        core::slice::from_ref(&valid_pair),
        cref!(MaximumMatchingCertificateRef, "test/inhabited-maximum", 96),
        cref!(CausalCheckerRef, "test/inhabited-checker", 97),
    )
    .expect("exact active graph matching witness binds");

    let mut wrong_order = valid_pair.clone();
    wrong_order.variable.derivative_order = wrong_order.variable.derivative_order.saturating_add(1);
    assert_eq!(
        maximum_matching_binding(
            &conditional_graph,
            mode_domain.clone(),
            core::slice::from_ref(&wrong_order),
            cref!(
                MaximumMatchingCertificateRef,
                "test/wrong-order-maximum",
                98
            ),
            cref!(CausalCheckerRef, "test/wrong-order-checker", 99),
        ),
        Err(MaximumMatchingBindingError::InvalidMatchingSet(
            MaximumMatchingWitnessError::EndpointMismatch { index: 0 }
        ))
    );
    let foreign_pair = complete_receipt(&unconditional_graph).matching[0].clone();
    assert_eq!(
        maximum_matching_binding(
            &conditional_graph,
            mode_domain.clone(),
            core::slice::from_ref(&foreign_pair),
            cref!(
                MaximumMatchingCertificateRef,
                "test/foreign-pair-maximum",
                100
            ),
            cref!(CausalCheckerRef, "test/foreign-pair-checker", 101),
        ),
        Err(MaximumMatchingBindingError::InvalidMatchingSet(
            MaximumMatchingWitnessError::ForeignIncidence { index: 0 }
        ))
    );
    let inactive_or_known = conditional_graph
        .incidences()
        .iter()
        .find(|incidence| incidence.id != valid_pair.incidence)
        .expect("hybrid graph has a second structural incidence");
    let inactive_pair = CausalMatchingPair {
        incidence: inactive_or_known.id.clone(),
        equation: inactive_or_known.equation.clone(),
        variable: DerivativeVariableKey {
            variable: inactive_or_known.variable.clone(),
            derivative_order: inactive_or_known.derivative_order,
        },
    };
    assert_eq!(
        maximum_matching_binding(
            &conditional_graph,
            mode_domain,
            core::slice::from_ref(&inactive_pair),
            cref!(
                MaximumMatchingCertificateRef,
                "test/inactive-pair-maximum",
                102
            ),
            cref!(CausalCheckerRef, "test/inactive-pair-checker", 103),
        ),
        Err(MaximumMatchingBindingError::InvalidMatchingSet(
            MaximumMatchingWitnessError::NonUnknownIncidence { index: 0 }
        ))
    );
}

#[test]
fn g4_pre_cancelled_public_identity_constructors_publish_nothing() {
    let (machine, owner, clock) = minimal_machine();
    let draft = minimal_causal_draft(&owner, &clock);
    let equation = draft.equations[0].clone();
    let variable = draft.variables[0].clone();
    let incidence_error = with_cancelled_cx(|cx| {
        IncidenceSpec::new(
            equation.id,
            variable.id,
            0,
            SolveParticipation::Unknown,
            Dims::NONE,
            equation.residual,
            None,
            IncidenceClockRelation::SameClock,
            ActivationDomain::Always,
            cx,
        )
    })
    .expect_err("pre-cancelled incidence construction must publish no identity");
    assert_eq!(
        incidence_error,
        CanonicalError::Cancelled { absorbed_bytes: 0 }
    );

    let graph = with_cx(|cx| draft.admit_against(&machine, cx))
        .expect("graph admits before theorem-constructor cancellation");
    let receipt = complete_receipt(&graph);
    let maximum_error = with_cancelled_cx(|cx| {
        MaximumMatchingBinding::new(
            &graph,
            CausalReceiptDomain::UnconditionalGraph,
            &receipt.matching,
            cref!(MaximumMatchingCertificateRef, "test/cancelled-maximum", 86),
            cref!(CausalCheckerRef, "test/cancelled-maximum-checker", 87),
            cx,
        )
    })
    .expect_err("pre-cancelled maximum binding must publish nothing");
    assert_eq!(
        maximum_error,
        MaximumMatchingBindingError::Identity(CanonicalError::Cancelled { absorbed_bytes: 0 })
    );

    let predecessor = HistoricalCausalIdentityReceipt::new(
        CausalMigrationArtifactKind::Structure,
        0,
        [122; 32],
        [123; 32],
        [124; 32],
        0,
        0,
        0,
    )
    .expect("legacy predecessor receipt");
    let migration = CausalSchemaMigrationDraft::for_structure(
        predecessor,
        &graph,
        cref!(CausalMigrationRef, "test/cancelled-migration", 125),
    );
    let migration_error = with_cancelled_cx(|cx| migration.admit(cx))
        .expect_err("pre-cancelled migration must publish no identity");
    assert_eq!(
        migration_error,
        CausalMigrationError::Identity(CanonicalError::Cancelled { absorbed_bytes: 0 })
    );

    let conditional_graph = hybrid_graph(&machine, &owner, &clock);
    let child = mode_cell_receipt(&conditional_graph, mode_assignment(&conditional_graph, 0));
    assert_eq!(
        with_cancelled_cx(|cx| ConditionalCausalOutcome::from_mode_cell(&child, cx)),
        Err(ConditionalOutcomeError::Cancelled),
        "pre-cancelled child construction must copy no assignment"
    );
    let outcome = conditional_outcome(&child).expect("mode child");
    let coverage_error = with_cancelled_cx(|cx| {
        ConditionalCoverageBinding::for_mode_cells(
            &conditional_graph,
            core::slice::from_ref(&outcome),
            cref!(ConditionalCoverageRef, "test/cancelled-coverage", 88),
            cref!(CausalCheckerRef, "test/cancelled-coverage-checker", 89),
            cx,
        )
    })
    .expect_err("pre-cancelled coverage binding must publish nothing");
    assert_eq!(
        coverage_error,
        ConditionalCoverageBindingError::Identity(CanonicalError::Cancelled { absorbed_bytes: 0 })
    );
}

#[test]
fn g4_pre_cancelled_graph_and_receipt_publish_no_identity() {
    let (machine, owner, clock) = minimal_machine();
    let graph_decision = with_cancelled_cx(|cx| {
        minimal_causal_draft(&owner, &clock).admit_with_decision(&machine, cx)
    });
    assert!(
        !graph_decision.submitted_counts().complete,
        "pre-cancellation must interrupt nested telemetry before graph work"
    );
    let graph_refusal = graph_decision
        .into_result()
        .expect_err("pre-cancelled graph admission refuses");
    assert_graph_rules_exact(&graph_refusal, &[CausalGraphRule::Cancelled]);

    let graph = with_cx(|cx| minimal_causal_draft(&owner, &clock).admit_against(&machine, cx))
        .expect("graph admits before receipt cancellation");
    let receipt_decision =
        with_cancelled_cx(|cx| complete_receipt(&graph).admit_with_decision(&graph, cx));
    assert!(
        !receipt_decision.submitted_counts().complete,
        "pre-cancellation must interrupt receipt telemetry before admission"
    );
    let receipt_refusal = receipt_decision
        .into_result()
        .expect_err("pre-cancelled receipt admission refuses");
    assert_receipt_rules_exact(&receipt_refusal, &[CausalReceiptRule::Cancelled]);
}
