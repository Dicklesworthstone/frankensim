//! Machine-IR E0 PR-4 assurance admission conformance (Gauntlet G0/G3).

use core::num::NonZeroU64;

use fs_blake3::ContentHash;
use fs_evidence::vv::*;
use fs_ir::machine::assurance::*;
use fs_ir::machine::assurance_codec::{
    MAX_MACHINE_ASSURANCE_VV_CASE_RECEIPTS, MachineAssuranceAstAdmissionError,
    MachineAssuranceCodecRule, admit_machine_assurance_ast_v1, parse_machine_assurance_program_v1,
    parse_machine_assurance_v1, write_machine_assurance_program_v1, write_machine_assurance_v1,
};
use fs_ir::machine::semantics::{
    AdmittedMachineBehavior, BodyMotion, ConditionBinding, ConditionSource, ConditionTarget,
    ConditionValueRef, MachineBehaviorDraft, MotionBinding, StateSlotContract,
};
use fs_ir::machine::{
    AdmittedMachineGraph, BodyId, ClockId, ClockSpec, FrameBinding, MachineClock,
    MachineGraphDraft, MaterialBinding, MaterialCardRef, MaterialTarget, ModelRef,
    OrientationParity, RelationId, RelationMode, RelationSpec, StateSlotId, SubsystemId,
    SubsystemSpec, TerminalCausality, TerminalId, TerminalQuantitySpec, TerminalShape,
    TerminalSpec,
};
use fs_ir::{Node, NodeKind, VersionedProgram};
use fs_qty::Dims;

fn nz(value: u64) -> NonZeroU64 {
    NonZeroU64::new(value).expect("test value is nonzero")
}

fn digest(byte: u8) -> [u8; 32] {
    [byte; 32]
}

fn hash(label: &str) -> ContentHash {
    fs_blake3::hash_domain(
        "org.frankensim.fs-ir.machine-assurance-test.v1",
        label.as_bytes(),
    )
}

fn artifact_id(value: &str) -> ArtifactId {
    ArtifactId::try_new(value).expect("valid fixture artifact id")
}

fn qoi_id(value: &str) -> QoiId {
    QoiId::try_new(value).expect("valid fixture QoI id")
}

fn observation_id(value: &str) -> ObservationId {
    ObservationId::try_new(value).expect("valid fixture observation id")
}

fn axis_id(value: &str) -> AxisId {
    AxisId::try_new(value).expect("valid fixture axis id")
}

fn unit_id(value: &str) -> UnitId {
    UnitId::try_new(value).expect("valid fixture unit id")
}

fn header(id: &str, units: &[&str]) -> ArtifactHeader {
    ArtifactHeader::try_new(
        artifact_id(id),
        units.iter().copied().map(unit_id).collect(),
        SeedDeclaration::Fixed(0x5eed),
        DeclaredBudget::Limit(1.0e-6),
        DeclaredBudget::Limit(10_000),
        DeclaredBudget::Limit(1 << 20),
        vec![("fs-evidence".to_owned(), "1.0.0".to_owned())],
        vec!["vv-artifacts".to_owned()],
    )
    .expect("valid fixture header")
}

fn external_target(label: &str) -> EvidenceTarget {
    EvidenceTarget::External {
        family: artifact_id("fixture-family"),
        id: artifact_id(label),
        hash: hash(label),
    }
}

fn artifact_reference<T>(artifact: &T) -> ArtifactRef
where
    T: Clone + Into<VvArtifact>,
{
    let artifact: VvArtifact = artifact.clone().into();
    ArtifactRef::new(
        artifact.kind(),
        artifact.id().clone(),
        artifact.content_hash().expect("canonical artifact hash"),
    )
}

fn diagnostic_record(rule: VvRule) -> DiagnosticRecord {
    DiagnosticRecord::try_new(true, hash(rule.slug()), format!("{} passed", rule.slug()))
        .expect("diagnostic fixture")
}

fn diagnostic_plan() -> DiagnosticPlan {
    DiagnosticPlan::new(
        diagnostic_record(VvRule::DiagnosticObservability),
        diagnostic_record(VvRule::DiagnosticIdentifiability),
        diagnostic_record(VvRule::DiagnosticConfounding),
        diagnostic_record(VvRule::DiagnosticInverseCrime),
    )
}

fn categorical_axes() -> EvidenceAxes {
    EvidenceAxes::try_new(
        EvidenceAxis::ALL
            .into_iter()
            .map(|axis| {
                (
                    axis,
                    EvidenceAxisStatus::Missing {
                        reason: "fixture makes no positive evidence-color claim".to_owned(),
                    },
                )
            })
            .collect(),
    )
    .expect("complete evidence axes")
}

fn uncertainty_term(
    kind: PredictionUncertaintyKind,
    magnitude: f64,
    source: EvidenceTarget,
) -> UncertaintyTerm {
    UncertaintyTerm::try_new(kind, magnitude, source).expect("valid uncertainty term")
}

#[allow(clippy::too_many_lines)]
fn admitted_vv_case(acceptance_hi: f64) -> AdmittedVvCase {
    admitted_vv_case_with_extra_observations(acceptance_hi, 0)
}

#[allow(clippy::too_many_lines)]
fn admitted_vv_case_with_extra_observations(
    acceptance_hi: f64,
    extra_observations: usize,
) -> AdmittedVvCase {
    assert!(extra_observations <= 4_093);
    let qoi = qoi_id("mass");
    let unit = unit_id("kg");
    let load_axis = axis_id("load");
    let regime_axis = axis_id("regime");
    let context = ContextOfUse::try_new(
        header("context-1", &["kg", "unitless"]),
        "Decide whether retained mass satisfies the release criterion.",
        vec![
            QoiSpec::try_new(
                qoi.clone(),
                "retained mass",
                unit.clone(),
                AcceptanceCriterion::ClosedRange {
                    lo: 9.0,
                    hi: acceptance_hi,
                },
            )
            .expect("QoI fixture"),
        ],
        ApplicabilityDomain::try_new(
            vec![
                NumericDomainAxis::try_new(load_axis.clone(), unit_id("unitless"), 1.0, 10.0)
                    .expect("numeric applicability fixture"),
            ],
            vec![
                CategoricalDomainAxis::try_new(
                    regime_axis.clone(),
                    vec!["nominal".to_owned(), "hot".to_owned()],
                )
                .expect("categorical applicability fixture"),
            ],
        )
        .expect("applicability fixture"),
        ApplicabilityPolicy::Demote,
    )
    .expect("context fixture");
    let context_ref = artifact_reference(&context);

    let source_bytes = hash("source-bytes");
    let instrument = artifact_id("instrument-1");
    let clock = artifact_id("clock-1");
    let typed_manifest_entry = |id: ObservationId, locator_hash: ContentHash| {
        let acquisition_channel = artifact_id(&format!("acquisition-channel-{}", id.as_str()));
        let source = ObservationSourceRef::try_new(
            source_bytes,
            "org.frankensim.fs-ir.machine-assurance-test.row-locator.v1",
            1,
            locator_hash,
            hash(&format!("extraction-receipt-{}", id.as_str())),
        )
        .expect("typed observation source");
        let row = ObservationManifestRow::try_new(
            source,
            qoi.clone(),
            instrument.clone(),
            acquisition_channel,
            clock.clone(),
        )
        .expect("typed observation manifest row");
        (id, row)
    };
    let mut manifest_rows = vec![
        typed_manifest_entry(observation_id("cal-1"), hash("row-cal-1")),
        typed_manifest_entry(observation_id("val-1"), hash("row-val-1")),
        typed_manifest_entry(observation_id("blind-1"), hash("row-blind-1")),
    ];
    let mut calibration_ids = vec![observation_id("cal-1")];
    let mut validation_ids = vec![observation_id("val-1")];
    let mut blind_rows = vec![(observation_id("blind-1"), hash("row-blind-1"))];
    for index in 0..extra_observations {
        let id = observation_id(&format!("extra-{index:04}"));
        let source = hash(&format!("row-extra-{index:04}"));
        manifest_rows.push(typed_manifest_entry(id.clone(), source));
        match index % 3 {
            0 => calibration_ids.push(id),
            1 => validation_ids.push(id),
            _ => blind_rows.push((id, source)),
        }
    }
    let replicates = u32::try_from(manifest_rows.len()).expect("bounded observation fixture");

    let experiment = ExperimentArtifact::try_new(
        header("experiment-1", &["kg"]),
        artifact_id("experiment-1-dataset"),
        ExperimentOrigin::Physical {
            apparatus_id: artifact_id("apparatus-1"),
            facility_id: artifact_id("facility-1"),
        },
        vec![qoi.clone()],
        ObservationManifest::try_new(manifest_rows).expect("injective observation manifest"),
        vec![InstrumentCalibration::new(
            instrument,
            hash("instrument-calibration"),
            true,
        )],
        ClockSynchronization::SingleClock { clock_id: clock },
        RepeatabilitySummary::try_new(
            replicates,
            CovarianceMatrix::try_new(1, vec![0.25]).expect("scalar covariance"),
        )
        .expect("repeatability fixture"),
        DataAuthenticity::new(source_bytes, hash("custody"), true),
    )
    .expect("experiment fixture");
    let experiment_ref = artifact_reference(&experiment);

    let split = CalibrationSplit::try_new(
        header("split-1", &["unitless"]),
        experiment_ref.clone(),
        hash("preregistered-analysis"),
        calibration_ids,
        validation_ids,
        blind_rows,
    )
    .expect("calibration split fixture");
    let split_ref = artifact_reference(&split);

    let validation_plan = ValidationPlan::try_new(
        header("validation-plan-1", &["kg"]),
        context_ref.clone(),
        vec![
            QoiValidationPlan::try_new(
                qoi.clone(),
                vec![experiment_ref.clone()],
                split_ref.clone(),
                vec![
                    ValidationMetricSpec::IntervalAgreement,
                    ValidationMetricSpec::PosteriorPredictive {
                        minimum_tail_probability: 0.05,
                    },
                ],
                diagnostic_plan(),
            )
            .expect("QoI validation plan"),
        ],
    )
    .expect("validation plan fixture");
    let validation_plan_ref = artifact_reference(&validation_plan);

    let numerical =
        |label| NumericalUncertainty::try_new(0.01, hash(label)).expect("numerical uncertainty");
    let solution = SolutionVerificationReceipt::try_new(
        header("solution-1", &["kg"]),
        artifact_id("solve-1"),
        qoi.clone(),
        unit.clone(),
        numerical("mesh-bound"),
        numerical("time-bound"),
        numerical("nonlinear-bound"),
        numerical("iterative-bound"),
    )
    .expect("solution verification fixture");
    let numerical_floor = solution.combined_half_width();
    let solution_ref = artifact_reference(&solution);
    let validation_selection = split
        .validation_selection(split_ref, vec![observation_id("val-1")])
        .expect("validation selection");

    let solution_source = EvidenceTarget::VvArtifact(solution_ref);
    let model_source = external_target("model-discrepancy");
    let parameter_source = external_target("parameter-data");
    let data_source = external_target("measurement-data");
    let aleatory_source = external_target("aleatory-model");
    let epistemic_source = external_target("epistemic-model");
    let mut dependencies = vec![
        EvidenceDependency::physical_validation(
            qoi.clone(),
            experiment_ref,
            validation_selection.clone(),
        ),
        EvidenceDependency::new(
            qoi.clone(),
            DependencyRole::SolutionVerification,
            solution_source.clone(),
        ),
        EvidenceDependency::new(
            qoi.clone(),
            DependencyRole::ModelDiscrepancy,
            model_source.clone(),
        ),
        EvidenceDependency::new(
            qoi.clone(),
            DependencyRole::ParameterData,
            parameter_source.clone(),
        ),
        EvidenceDependency::new(
            qoi.clone(),
            DependencyRole::ParameterData,
            data_source.clone(),
        ),
        EvidenceDependency::new(
            qoi.clone(),
            DependencyRole::ParameterData,
            aleatory_source.clone(),
        ),
        EvidenceDependency::new(
            qoi.clone(),
            DependencyRole::ModelDiscrepancy,
            epistemic_source.clone(),
        ),
    ];
    dependencies.sort_by_key(|dependency| dependency.target().hash());

    let waterfall = UncertaintyWaterfall::try_new(
        qoi.clone(),
        unit,
        WaterfallMode::GuaranteedBound,
        vec![
            uncertainty_term(PredictionUncertaintyKind::ModelForm, 0.1, model_source),
            uncertainty_term(PredictionUncertaintyKind::Parameter, 0.1, parameter_source),
            uncertainty_term(
                PredictionUncertaintyKind::Numerical,
                numerical_floor,
                solution_source,
            ),
            uncertainty_term(PredictionUncertaintyKind::Data, 0.1, data_source),
            uncertainty_term(PredictionUncertaintyKind::Aleatory, 0.1, aleatory_source),
            uncertainty_term(PredictionUncertaintyKind::Epistemic, 0.1, epistemic_source),
        ],
    )
    .expect("uncertainty waterfall");
    let metric = ValidationMetric::try_new(
        artifact_id("interval-agreement"),
        qoi.clone(),
        validation_selection.clone(),
        9.9,
        9.95,
        0.05,
        numerical_floor,
    )
    .expect("validation metric");
    let posterior = PosteriorPredictiveCheck::try_new(
        artifact_id("posterior-check"),
        qoi.clone(),
        validation_selection,
        0.5,
        0.05,
        hash("posterior-check-artifact"),
    )
    .expect("posterior check");

    let mut assumptions = AssumptionsLedger::try_program_seed(header("assumptions", &["unitless"]))
        .expect("program assumptions");
    for row in assumptions.rows().values().cloned().collect::<Vec<_>>() {
        let label = format!("{}-evidence", row.id().as_str());
        assumptions
            .replace_row(
                row.with_evidence(external_target(&label))
                    .with_monitor_evidence(hash(&format!("{label}-monitor"))),
            )
            .expect("attach assumption evidence");
    }
    let assumption_checks = assumptions
        .rows()
        .keys()
        .cloned()
        .map(|id| (id, true))
        .collect();
    let prediction = PredictionAssessment::try_new(
        header("prediction-1", &["kg"]),
        context_ref,
        validation_plan_ref,
        qoi,
        dependencies,
        waterfall,
        vec![metric],
        vec![posterior],
        ApplicabilityPoint::try_new(
            vec![(load_axis, 5.0)],
            vec![(regime_axis, "nominal".to_owned())],
        )
        .expect("applicability point"),
        ApplicabilityDecision::InDomain,
        categorical_axes(),
        assumption_checks,
    )
    .expect("prediction fixture");

    VvCase::try_new(
        context,
        validation_plan,
        vec![experiment],
        vec![split],
        vec![solution],
        vec![prediction],
        assumptions,
    )
    .expect("closed V&V case")
    .admit()
    .expect("V&V case admits")
}

fn frame() -> FrameBinding {
    FrameBinding::new("world/mechanical", OrientationParity::Preserving).expect("valid frame")
}

fn quantity() -> TerminalQuantitySpec {
    TerminalQuantitySpec::Dimensional(Dims([0, 1, 0, 0, 0, 0]))
}

fn model(byte: u8) -> ModelRef {
    ModelRef::new("models/mass-plant", nz(1), digest(byte)).expect("valid model ref")
}

fn material(byte: u8) -> MaterialCardRef {
    MaterialCardRef::new("materials/mass-plant", nz(1), digest(byte)).expect("valid material ref")
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

fn valid_graph() -> MachineGraphDraft {
    let continuous = ClockId::new("clock/continuous").expect("valid clock id");
    let sampled = ClockId::new("clock/sampled").expect("valid clock id");
    let subsystem = SubsystemId::new("subsystem/plant").expect("valid subsystem id");
    let state = StateSlotId::new("state/mass").expect("valid state id");
    let body = BodyId::new("body/plant").expect("valid body id");
    let source = terminal(
        "terminal/mass-source",
        &subsystem,
        TerminalCausality::Output,
        &continuous,
    );
    let sink = terminal(
        "terminal/mass-sink",
        &subsystem,
        TerminalCausality::Input,
        &continuous,
    );
    let sensor_output = terminal(
        "terminal/mass-sensor",
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
                id: sampled,
                clock: MachineClock::Periodic {
                    period_ns: nz(1_000_000),
                    phase_ns: 0,
                },
            },
        ],
        subsystems: vec![SubsystemSpec {
            id: subsystem,
            model: model(1),
            bodies: vec![body.clone()],
            surface_patches: Vec::new(),
            contact_features: Vec::new(),
            state_slots: vec![state.clone()],
        }],
        terminals: vec![source.clone(), sink.clone(), sensor_output],
        ports: Vec::new(),
        relations: vec![RelationSpec {
            id: RelationId::new("relation/mass-state").expect("valid relation id"),
            source: source.id,
            target: sink.id,
            mode: RelationMode::Stateful { state_slot: state },
        }],
        materials: vec![MaterialBinding {
            target: MaterialTarget::Body(body),
            material: material(2),
        }],
        interfaces: Vec::new(),
    }
}

fn valid_behavior() -> MachineBehaviorDraft {
    let subsystem = SubsystemId::new("subsystem/plant").expect("valid subsystem id");
    let state = StateSlotId::new("state/mass").expect("valid state id");
    let continuous = ClockId::new("clock/continuous").expect("valid clock id");
    MachineBehaviorDraft {
        state_contracts: vec![StateSlotContract {
            id: state.clone(),
            owner: subsystem,
            quantity: quantity(),
            shape: TerminalShape::Scalar,
            clock: continuous.clone(),
            frame: frame(),
        }],
        conditions: vec![ConditionBinding {
            target: ConditionTarget::Initial(state),
            quantity: quantity(),
            shape: TerminalShape::Scalar,
            clock: continuous.clone(),
            frame: frame(),
            source: ConditionSource::Fixed(
                ConditionValueRef::new("values/initial-mass", nz(1), digest(3))
                    .expect("valid condition value"),
            ),
        }],
        motions: vec![MotionBinding {
            body: BodyId::new("body/plant").expect("valid body id"),
            clock: continuous,
            reference_frame: frame(),
            motion: BodyMotion::Static,
        }],
        events: Vec::new(),
        tolerances: Vec::new(),
        dependences: Vec::new(),
    }
}

fn artifact_ref(admitted: &AdmittedVvCase, kind: ArtifactKind, id: &ArtifactId) -> ArtifactRef {
    let hash = admitted
        .receipt()
        .artifact_hashes()
        .get(&(kind, id.clone()))
        .copied()
        .expect("fixture artifact has admitted hash");
    ArtifactRef::new(kind, id.clone(), hash)
}

fn refs(admitted: &AdmittedVvCase) -> (ArtifactRef, ArtifactRef, ArtifactRef) {
    let case = admitted.case();
    (
        artifact_ref(admitted, ArtifactKind::ContextOfUse, case.context().id()),
        artifact_ref(
            admitted,
            ArtifactKind::ValidationPlan,
            case.validation_plan().id(),
        ),
        artifact_ref(
            admitted,
            ArtifactKind::ExperimentArtifact,
            case.experiments()
                .keys()
                .next()
                .expect("fixture experiment"),
        ),
    )
}

macro_rules! aref {
    ($name:ident, $namespace:literal, $byte:expr) => {
        $name::new($namespace, nz(1), digest($byte)).expect("valid assurance reference")
    };
}

fn valid_assurance(admitted: &AdmittedVvCase) -> MachineAssuranceDraft {
    let (context, validation_plan, experiment) = refs(admitted);
    let qoi = qoi_id("mass");
    let qoi_key = ContextQoiKey {
        context: context.id().clone(),
        qoi: qoi.clone(),
    };
    let sensor = SensorId::new("sensor/mass").expect("valid sensor id");
    let subsystem = SubsystemId::new("subsystem/plant").expect("valid subsystem id");
    let baseline = FidelityRungId::new("fidelity/plant-baseline").expect("valid rung id");
    let hazard = HazardId::new("hazard/mass-loss").expect("valid hazard id");
    MachineAssuranceDraft {
        sensors: vec![SensorSpec {
            id: sensor.clone(),
            owner: subsystem.clone(),
            target: ObservationTarget::State(
                StateSlotId::new("state/mass").expect("valid state id"),
            ),
            quantity: quantity(),
            shape: TerminalShape::Scalar,
            clock: ClockId::new("clock/continuous").expect("valid clock id"),
            frame: frame(),
            timing: ObservationTiming::Direct,
            model: aref!(SensorModelRef, "sensors/mass-model", 10),
            calibration: aref!(CalibrationRef, "calibrations/mass", 11),
            exposure: SensorExposure::PlantSignal {
                output: TerminalId::new("terminal/mass-sensor").expect("valid terminal id"),
            },
        }],
        experiments: vec![ExperimentSpec {
            id: ExperimentId::new("experiment/mass-release").expect("valid experiment id"),
            artifact: experiment,
            context: context.clone(),
            instruments: vec![SensorInstrumentBinding {
                sensor: sensor.clone(),
                instrument: artifact_id("instrument-1"),
            }],
            qois: vec![qoi.clone()],
        }],
        contexts: vec![ContextBinding {
            context: context.clone(),
            validation_plan,
            qois: vec![QoiBinding {
                id: qoi,
                inputs: vec![QoiInput {
                    target: QoiTarget::Sensor(sensor),
                    quantity: quantity(),
                    shape: TerminalShape::Scalar,
                }],
                unit: unit_id("kg"),
                definition: aref!(QoiDefinitionRef, "qois/retained-mass", 12),
                unit_bridge: aref!(UnitQuantityBridgeRef, "units/kg-to-mass", 13),
            }],
            budget: aref!(DecisionBudgetRef, "budgets/release", 14),
        }],
        hazards: vec![HazardSpec {
            id: hazard.clone(),
            context: context.clone(),
            scope: vec![MachineScope::WholeMachine],
            requirement: aref!(SafetyRequirementRef, "requirements/retain-mass", 15),
            operating_envelope: aref!(OperatingEnvelopeRef, "envelopes/release", 16),
            safety_case: aref!(SafetyCaseRef, "safety-cases/mass-loss", 17),
            assumptions: vec![AssumptionId::try_new("A-001").expect("seed assumption")],
            fault_coverage: FaultCoverage::Modeled,
        }],
        faults: vec![FaultSpec {
            id: FaultId::new("fault/leak").expect("valid fault id"),
            affected: vec![MachineScope::WholeMachine],
            hazards: vec![hazard],
            model: aref!(FaultModelRef, "faults/leak-model", 18),
            containment: aref!(FaultContainmentRef, "containment/leak", 19),
            injection: aref!(FaultInjectionRef, "injections/leak", 20),
        }],
        accounting_windows: vec![AccountingWindow {
            id: AccountingWindowId::new("accounting/mass-window").expect("valid accounting id"),
            context,
            clock: ClockId::new("clock/continuous").expect("valid clock id"),
            balance: BalanceKind::Mass,
            quantity: quantity(),
            boundary: aref!(AccountingBoundaryRef, "boundaries/plant", 21),
            interval: aref!(AccountingIntervalRef, "intervals/release", 22),
            entries: vec![AccountingEntry {
                target: AccountingTarget::State(
                    StateSlotId::new("state/mass").expect("valid state id"),
                ),
                role: AccountingRole::Storage,
                orientation: AccountingOrientation::StoredIncreasePositive,
                policy: aref!(AccountingPolicyRef, "accounting/storage", 23),
                loss_ownership: None,
            }],
            audit_policy: aref!(AccountingPolicyRef, "accounting/window-audit", 24),
        }],
        fidelity: FidelityPolicy {
            baselines: vec![baseline.clone()],
            rungs: vec![FidelityRung {
                id: baseline.clone(),
                subsystem,
                model: model(1),
                model_crosswalk: aref!(ModelCrosswalkRef, "crosswalks/baseline", 25),
                validity_domain: aref!(ValidityDomainRef, "validity/baseline", 26),
                cost_error_model: aref!(CostErrorModelRef, "cost-error/baseline", 27),
                falsifiers: vec![aref!(FalsifierRef, "falsifiers/mass", 28)],
                qois: vec![qoi_key],
            }],
            escalations: vec![EscalationSpec {
                from: baseline,
                trigger: aref!(EscalationTriggerRef, "triggers/baseline-exit", 29),
                action: EscalationAction::Refuse(aref!(
                    NoClaimRef,
                    "no-claims/no-higher-fidelity",
                    30
                )),
            }],
            fixed_replay: aref!(FixedReplayRef, "replay/fixed-baseline", 31),
        },
    }
}

fn admitted_stack(
    acceptance_hi: f64,
) -> (
    AdmittedMachineGraph,
    AdmittedMachineBehavior,
    AdmittedVvCase,
) {
    let graph = valid_graph().admit().expect("machine graph admits");
    let behavior = valid_behavior()
        .admit_against(&graph)
        .expect("machine behavior admits");
    (graph, behavior, admitted_vv_case(acceptance_hi))
}

fn rules(refusal: &MachineAssuranceRefusal) -> Vec<MachineAssuranceRule> {
    refusal
        .findings()
        .iter()
        .map(|finding| finding.rule())
        .collect()
}

#[test]
fn g0_fully_populated_assurance_admits_with_exact_receipts() {
    let (graph, behavior, case) = admitted_stack(11.0);
    let decision = valid_assurance(&case).admit_with_decision(&graph, &behavior, &[case.clone()]);
    assert_eq!(decision.code(), "MachineAssuranceAdmitted");
    assert_eq!(decision.submitted_counts().sensors, 1);
    assert_eq!(decision.submitted_counts().experiments, 1);
    assert_eq!(decision.submitted_counts().contexts, 1);
    assert_eq!(decision.submitted_counts().hazards, 1);
    assert_eq!(decision.submitted_counts().faults, 1);
    let admitted = decision.result().expect("assurance admits");
    assert_eq!(admitted.base_graph(), graph.identity());
    assert_eq!(admitted.base_behavior(), behavior.identity());
    assert_eq!(admitted.vv_cases().len(), 1);
    assert_eq!(
        admitted.vv_cases()[0].receipt_hash,
        case.receipt().receipt_hash()
    );
    assert_eq!(admitted.vv_cases()[0].case_hash, case.receipt().case_hash());
}

#[test]
fn g3_vv_case_acceptance_and_reference_drift_move_or_refuse() {
    let (graph, behavior, case) = admitted_stack(11.0);
    let baseline = valid_assurance(&case)
        .admit_against(&graph, &behavior, &[case.clone()])
        .expect("baseline admits");

    let changed_case = admitted_vv_case(12.0);
    let changed = valid_assurance(&changed_case)
        .admit_against(&graph, &behavior, &[changed_case.clone()])
        .expect("changed accepted criterion remains structurally admissible");
    assert_ne!(changed.identity(), baseline.identity());
    assert_ne!(
        changed.vv_cases()[0].receipt_hash,
        baseline.vv_cases()[0].receipt_hash
    );

    let mut stale = valid_assurance(&case);
    stale.contexts[0].context = ArtifactRef::new(
        ArtifactKind::ContextOfUse,
        stale.contexts[0].context.id().clone(),
        hash("stale-context"),
    );
    let refusal = stale
        .admit_against(&graph, &behavior, &[case.clone()])
        .expect_err("stale context hash refuses");
    assert!(rules(&refusal).contains(&MachineAssuranceRule::InvalidArtifactReference));

    let refusal = valid_assurance(&case)
        .admit_against(&graph, &behavior, &[])
        .expect_err("missing admitted V&V case refuses");
    assert!(rules(&refusal).contains(&MachineAssuranceRule::MissingVvCase));

    let mut wrong_unit = valid_assurance(&case);
    wrong_unit.contexts[0].qois[0].unit = unit_id("g");
    let refusal = wrong_unit
        .admit_against(&graph, &behavior, &[case.clone()])
        .expect_err("machine QoI unit must equal the admitted Context-of-Use unit");
    assert!(rules(&refusal).contains(&MachineAssuranceRule::QoiUnitMismatch));

    let mut missing_experiment = valid_assurance(&case);
    missing_experiment.experiments.clear();
    let refusal = missing_experiment
        .admit_against(&graph, &behavior, &[case.clone()])
        .expect_err("every admitted case experiment needs one machine binding");
    assert!(rules(&refusal).contains(&MachineAssuranceRule::MissingExperimentBinding));
}

#[test]
fn g0_sensor_clock_and_experiment_instrument_closure_fail_closed() {
    let (graph, behavior, case) = admitted_stack(11.0);

    let mut direct_mismatch = valid_assurance(&case);
    direct_mismatch.sensors[0].clock = ClockId::new("clock/sampled").expect("valid sampled clock");
    direct_mismatch.sensors[0].exposure = SensorExposure::ExperimentOnly;
    let refusal = direct_mismatch
        .admit_against(&graph, &behavior, &[case.clone()])
        .expect_err("unbridged clock mismatch refuses");
    assert!(rules(&refusal).contains(&MachineAssuranceRule::SensorClockGap));

    let mut bridged = valid_assurance(&case);
    bridged.sensors[0].clock = ClockId::new("clock/sampled").expect("valid sampled clock");
    bridged.sensors[0].timing = ObservationTiming::ModeledResampling {
        bridge: aref!(SamplingBridgeRef, "sampling/continuous-to-periodic", 90),
    };
    bridged.sensors[0].exposure = SensorExposure::ExperimentOnly;
    bridged
        .admit_against(&graph, &behavior, &[case.clone()])
        .expect("explicit resampling bridge admits");

    let mut unknown_instrument = valid_assurance(&case);
    unknown_instrument.experiments[0].instruments[0].instrument = artifact_id("instrument-missing");
    let refusal = unknown_instrument
        .admit_against(&graph, &behavior, &[case.clone()])
        .expect_err("instrument mismatch refuses");
    let observed = rules(&refusal);
    assert!(observed.contains(&MachineAssuranceRule::UnknownExperimentInstrument));
    assert!(observed.contains(&MachineAssuranceRule::ExperimentInstrumentSetMismatch));

    let mut input_as_sensor_output = valid_assurance(&case);
    input_as_sensor_output.sensors[0].exposure = SensorExposure::PlantSignal {
        output: TerminalId::new("terminal/mass-sink").expect("valid terminal id"),
    };
    let refusal = input_as_sensor_output
        .admit_against(&graph, &behavior, &[case.clone()])
        .expect_err("an input terminal cannot expose a plant sensor signal");
    assert!(rules(&refusal).contains(&MachineAssuranceRule::SensorOutputCausalityGap));

    let mut duplicate_sensor = valid_assurance(&case);
    let mut duplicate = duplicate_sensor.experiments[0].instruments[0].clone();
    duplicate.instrument = artifact_id("instrument-missing");
    duplicate_sensor.experiments[0].instruments.push(duplicate);
    let refusal = duplicate_sensor
        .admit_against(&graph, &behavior, &[case.clone()])
        .expect_err("one sensor cannot impersonate two instruments");
    assert!(rules(&refusal).contains(&MachineAssuranceRule::DuplicateExperimentSensor));

    let mut duplicate_artifact = valid_assurance(&case);
    let mut second_binding = duplicate_artifact.experiments[0].clone();
    second_binding.id =
        ExperimentId::new("experiment/mass-release-alias").expect("valid experiment id");
    duplicate_artifact.experiments.push(second_binding);
    let refusal = duplicate_artifact
        .admit_against(&graph, &behavior, &[case.clone()])
        .expect_err("one evidence artifact cannot have two local experiment identities");
    assert!(rules(&refusal).contains(&MachineAssuranceRule::DuplicateExperimentArtifactBinding));
}

#[test]
fn g0_hazard_assumptions_and_fault_coverage_are_honest() {
    let (graph, behavior, case) = admitted_stack(11.0);

    let mut uncovered = valid_assurance(&case);
    uncovered.faults.clear();
    let refusal = uncovered
        .admit_against(&graph, &behavior, &[case.clone()])
        .expect_err("modeled hazard without a fault refuses");
    assert!(rules(&refusal).contains(&MachineAssuranceRule::UncoveredModeledHazard));

    let mut honestly_unmodeled = valid_assurance(&case);
    honestly_unmodeled.hazards[0].fault_coverage =
        FaultCoverage::Unmodeled(aref!(NoClaimRef, "no-claims/fault-model-absent", 91));
    honestly_unmodeled.faults.clear();
    honestly_unmodeled
        .admit_against(&graph, &behavior, &[case.clone()])
        .expect("explicitly unmodeled hazard permits an empty fault set");

    let mut contradictory = valid_assurance(&case);
    contradictory.hazards[0].fault_coverage =
        FaultCoverage::Unmodeled(aref!(NoClaimRef, "no-claims/fault-model-absent", 92));
    let refusal = contradictory
        .admit_against(&graph, &behavior, &[case.clone()])
        .expect_err("unmodeled hazard cannot also receive a fault edge");
    assert!(rules(&refusal).contains(&MachineAssuranceRule::ContradictoryFaultCoverage));

    let mut unknown_assumption = valid_assurance(&case);
    unknown_assumption.hazards[0].assumptions =
        vec![AssumptionId::try_new("hazard-assumption-missing").expect("valid assumption id")];
    let refusal = unknown_assumption
        .admit_against(&graph, &behavior, &[case.clone()])
        .expect_err("hazard assumption must come from the exact V&V case");
    assert!(rules(&refusal).contains(&MachineAssuranceRule::UnknownHazardAssumption));
}

#[test]
fn g0_accounting_checks_balance_sign_target_and_loss_ownership() {
    let (graph, behavior, case) = admitted_stack(11.0);

    let mut wrong_balance = valid_assurance(&case);
    wrong_balance.accounting_windows[0].balance = BalanceKind::Energy;
    let refusal = wrong_balance
        .admit_against(&graph, &behavior, &[case.clone()])
        .expect_err("mass quantity cannot claim an energy balance");
    assert!(rules(&refusal).contains(&MachineAssuranceRule::AccountingBalanceQuantityGap));

    let mut wrong_sign = valid_assurance(&case);
    wrong_sign.accounting_windows[0].entries[0].orientation =
        AccountingOrientation::IntoBoundaryPositive;
    let refusal = wrong_sign
        .admit_against(&graph, &behavior, &[case.clone()])
        .expect_err("storage requires stored-increase sign semantics");
    assert!(rules(&refusal).contains(&MachineAssuranceRule::InvalidAccountingRole));

    let mut missing_loss_owner = valid_assurance(&case);
    missing_loss_owner.accounting_windows[0].entries[0] = AccountingEntry {
        target: AccountingTarget::Relation(
            RelationId::new("relation/mass-state").expect("valid relation id"),
        ),
        role: AccountingRole::Dissipation,
        orientation: AccountingOrientation::NonnegativeLoss,
        policy: aref!(AccountingPolicyRef, "accounting/dissipation", 93),
        loss_ownership: None,
    };
    let refusal = missing_loss_owner
        .admit_against(&graph, &behavior, &[case.clone()])
        .expect_err("dissipation needs exact loss ownership");
    assert!(rules(&refusal).contains(&MachineAssuranceRule::MissingLossOwnership));

    let mut duplicate_target = valid_assurance(&case);
    let duplicate_entry = duplicate_target.accounting_windows[0].entries[0].clone();
    duplicate_target.accounting_windows[0]
        .entries
        .push(duplicate_entry);
    let refusal = duplicate_target
        .admit_against(&graph, &behavior, &[case.clone()])
        .expect_err("one window cannot count one target twice");
    assert!(rules(&refusal).contains(&MachineAssuranceRule::DuplicateAccountingEntry));

    let mut source_and_state = valid_assurance(&case);
    source_and_state.accounting_windows[0]
        .entries
        .push(AccountingEntry {
            target: AccountingTarget::Relation(
                RelationId::new("relation/mass-state").expect("valid relation id"),
            ),
            role: AccountingRole::IncludedSource,
            orientation: AccountingOrientation::IntoBoundaryPositive,
            policy: aref!(AccountingPolicyRef, "accounting/included-source", 99),
            loss_ownership: None,
        });
    source_and_state
        .admit_against(&graph, &behavior, &[case.clone()])
        .expect("a relation source and its state storage are distinct contributions");

    let mut aliased_storage = valid_assurance(&case);
    aliased_storage.accounting_windows[0]
        .entries
        .push(AccountingEntry {
            target: AccountingTarget::Relation(
                RelationId::new("relation/mass-state").expect("valid relation id"),
            ),
            role: AccountingRole::Storage,
            orientation: AccountingOrientation::StoredIncreasePositive,
            policy: aref!(AccountingPolicyRef, "accounting/aliased-storage", 100),
            loss_ownership: None,
        });
    let refusal = aliased_storage
        .admit_against(&graph, &behavior, &[case.clone()])
        .expect_err("stateful relation storage aliases its exact state storage");
    assert!(rules(&refusal).contains(&MachineAssuranceRule::OverlappingAccountingTarget));
}

#[test]
fn g0_fidelity_baselines_and_termination_fail_closed() {
    let (graph, behavior, case) = admitted_stack(11.0);

    let mut wrong_model = valid_assurance(&case);
    wrong_model.fidelity.rungs[0].model = model(99);
    let refusal = wrong_model
        .admit_against(&graph, &behavior, &[case.clone()])
        .expect_err("baseline must equal the graph model");
    assert!(rules(&refusal).contains(&MachineAssuranceRule::BaselineModelMismatch));

    let mut missing_transition = valid_assurance(&case);
    missing_transition.fidelity.escalations.clear();
    let refusal = missing_transition
        .admit_against(&graph, &behavior, &[case.clone()])
        .expect_err("each rung needs an explicit outgoing action");
    assert!(rules(&refusal).contains(&MachineAssuranceRule::MissingEscalation));

    let mut cyclic = valid_assurance(&case);
    let baseline = cyclic.fidelity.rungs[0].id.clone();
    cyclic.fidelity.escalations[0].action = EscalationAction::Escalate {
        target: baseline,
        transfer: aref!(StateTransferRef, "transfers/self", 94),
        crosswalk: aref!(ModelCrosswalkRef, "crosswalks/self", 95),
    };
    let refusal = cyclic
        .admit_against(&graph, &behavior, &[case.clone()])
        .expect_err("self cycle cannot masquerade as escalation");
    let observed = rules(&refusal);
    assert!(observed.contains(&MachineAssuranceRule::SelfEscalation));
    assert!(observed.contains(&MachineAssuranceRule::FidelityEscalationCycle));

    let mut unreachable = valid_assurance(&case);
    let mut extra = unreachable.fidelity.rungs[0].clone();
    extra.id = FidelityRungId::new("fidelity/plant-unreachable").expect("valid rung id");
    extra.model = model(96);
    unreachable.fidelity.rungs.push(extra.clone());
    unreachable.fidelity.escalations.push(EscalationSpec {
        from: extra.id,
        trigger: aref!(EscalationTriggerRef, "triggers/unreachable", 97),
        action: EscalationAction::Refuse(aref!(NoClaimRef, "no-claims/unreachable", 98)),
    });
    let refusal = unreachable
        .admit_against(&graph, &behavior, &[case.clone()])
        .expect_err("every rung must be reachable from the subsystem baseline");
    assert!(rules(&refusal).contains(&MachineAssuranceRule::UnreachableFidelityRung));

    let mut drops_qoi = permutation_assurance(&case);
    let high = FidelityRungId::new("fidelity/plant-high").expect("valid rung id");
    drops_qoi
        .fidelity
        .rungs
        .iter_mut()
        .find(|rung| rung.id == high)
        .expect("high rung")
        .qois
        .clear();
    let refusal = drops_qoi
        .admit_against(&graph, &behavior, &[case.clone()])
        .expect_err("escalation cannot drop a decision QoI");
    assert!(rules(&refusal).contains(&MachineAssuranceRule::FidelityQoiDrop));
}

#[test]
fn g0_resource_refusal_precedes_deep_graph_or_evidence_work() {
    let (graph, behavior, case) = admitted_stack(11.0);
    let mut oversized = valid_assurance(&case);
    oversized.sensors = vec![oversized.sensors[0].clone(); MAX_MACHINE_ASSURANCE_SENSORS + 1];
    let decision = oversized.admit_with_decision(&graph, &behavior, &[]);
    assert_eq!(
        decision.submitted_counts().sensors,
        MAX_MACHINE_ASSURANCE_SENSORS + 1
    );
    let refusal = decision.result().expect_err("public sensor cap refuses");
    assert_eq!(rules(refusal), vec![MachineAssuranceRule::ResourceLimit]);
}

#[test]
fn g0_aggregate_vv_observation_manifests_hit_the_preflight_cap() {
    let graph = valid_graph().admit().expect("machine graph admits");
    let behavior = valid_behavior()
        .admit_against(&graph)
        .expect("machine behavior admits");
    let case = admitted_vv_case_with_extra_observations(11.0, 4_093);
    let draft = valid_assurance(&case);
    let vv_cases = vec![case; 17];
    let decision = draft.admit_with_decision(&graph, &behavior, &vv_cases);
    assert_eq!(decision.submitted_counts().vv_cases, 17);
    assert!(
        decision.submitted_counts().vv_nested_references > MAX_MACHINE_ASSURANCE_NESTED_REFERENCES
    );
    let refusal = decision
        .result()
        .expect_err("aggregate nested V&V content must refuse before receipt verification");
    assert_eq!(rules(refusal), vec![MachineAssuranceRule::ResourceLimit]);
}

fn permutation_assurance(case: &AdmittedVvCase) -> MachineAssuranceDraft {
    let mut draft = valid_assurance(case);
    draft.sensors.push(SensorSpec {
        id: SensorId::new("sensor/mass-backup").expect("valid sensor id"),
        owner: SubsystemId::new("subsystem/plant").expect("valid subsystem id"),
        target: ObservationTarget::State(StateSlotId::new("state/mass").expect("valid state id")),
        quantity: quantity(),
        shape: TerminalShape::Scalar,
        clock: ClockId::new("clock/continuous").expect("valid clock id"),
        frame: frame(),
        timing: ObservationTiming::Direct,
        model: aref!(SensorModelRef, "sensors/mass-backup", 110),
        calibration: aref!(CalibrationRef, "calibrations/mass-backup", 111),
        exposure: SensorExposure::ExperimentOnly,
    });
    draft.contexts[0].qois[0].inputs.push(QoiInput {
        target: QoiTarget::State(StateSlotId::new("state/mass").expect("valid state id")),
        quantity: quantity(),
        shape: TerminalShape::Scalar,
    });
    let context = draft.contexts[0].context.clone();
    draft.hazards[0].scope.push(MachineScope::Subsystem(
        SubsystemId::new("subsystem/plant").expect("valid subsystem id"),
    ));
    draft.hazards[0]
        .assumptions
        .push(AssumptionId::try_new("A-002").expect("seed assumption"));
    draft.hazards.push(HazardSpec {
        id: HazardId::new("hazard/uncatalogued").expect("valid hazard id"),
        context: context.clone(),
        scope: vec![MachineScope::WholeMachine],
        requirement: aref!(SafetyRequirementRef, "requirements/catalogue", 112),
        operating_envelope: aref!(OperatingEnvelopeRef, "envelopes/catalogue", 113),
        safety_case: aref!(SafetyCaseRef, "safety-cases/catalogue", 114),
        assumptions: vec![AssumptionId::try_new("A-003").expect("seed assumption")],
        fault_coverage: FaultCoverage::Unmodeled(aref!(
            NoClaimRef,
            "no-claims/uncatalogued-faults",
            115
        )),
    });
    draft.faults[0].affected.push(MachineScope::Subsystem(
        SubsystemId::new("subsystem/plant").expect("valid subsystem id"),
    ));
    let mut second_window = draft.accounting_windows[0].clone();
    second_window.id =
        AccountingWindowId::new("accounting/mass-window-backup").expect("valid accounting id");
    second_window.boundary = aref!(AccountingBoundaryRef, "boundaries/plant-backup", 117);
    second_window.interval = aref!(AccountingIntervalRef, "intervals/release-backup", 118);
    draft.accounting_windows.push(second_window);

    let high = FidelityRungId::new("fidelity/plant-high").expect("valid rung id");
    let qois = draft.fidelity.rungs[0].qois.clone();
    let baseline = draft.fidelity.rungs[0].id.clone();
    draft.fidelity.rungs[0]
        .falsifiers
        .push(aref!(FalsifierRef, "falsifiers/mass-secondary", 119));
    draft.fidelity.rungs.push(FidelityRung {
        id: high.clone(),
        subsystem: SubsystemId::new("subsystem/plant").expect("valid subsystem id"),
        model: model(120),
        model_crosswalk: aref!(ModelCrosswalkRef, "crosswalks/high", 121),
        validity_domain: aref!(ValidityDomainRef, "validity/high", 122),
        cost_error_model: aref!(CostErrorModelRef, "cost-error/high", 123),
        falsifiers: vec![
            aref!(FalsifierRef, "falsifiers/high-a", 124),
            aref!(FalsifierRef, "falsifiers/high-b", 125),
        ],
        qois,
    });
    draft.fidelity.escalations[0] = EscalationSpec {
        from: baseline,
        trigger: aref!(EscalationTriggerRef, "triggers/escalate-high", 126),
        action: EscalationAction::Escalate {
            target: high.clone(),
            transfer: aref!(StateTransferRef, "transfers/baseline-to-high", 127),
            crosswalk: aref!(ModelCrosswalkRef, "crosswalks/baseline-to-high", 128),
        },
    };
    draft.fidelity.escalations.push(EscalationSpec {
        from: high,
        trigger: aref!(EscalationTriggerRef, "triggers/high-exit", 129),
        action: EscalationAction::Refuse(aref!(NoClaimRef, "no-claims/high-exit", 130)),
    });
    draft
}

#[test]
fn g3_outer_and_nested_permutations_preserve_assurance_identity() {
    let (graph, behavior, case) = admitted_stack(11.0);
    let draft = permutation_assurance(&case);
    let expected = draft
        .clone()
        .admit_against(&graph, &behavior, &[case.clone()])
        .expect("rich assurance fixture admits");
    let mut permuted = draft;
    permuted.sensors.reverse();
    permuted.experiments.reverse();
    for experiment in &mut permuted.experiments {
        experiment.instruments.reverse();
        experiment.qois.reverse();
    }
    permuted.contexts.reverse();
    for context in &mut permuted.contexts {
        context.qois.reverse();
        for qoi in &mut context.qois {
            qoi.inputs.reverse();
        }
    }
    permuted.hazards.reverse();
    for hazard in &mut permuted.hazards {
        hazard.scope.reverse();
        hazard.assumptions.reverse();
    }
    permuted.faults.reverse();
    for fault in &mut permuted.faults {
        fault.affected.reverse();
        fault.hazards.reverse();
    }
    permuted.accounting_windows.reverse();
    for window in &mut permuted.accounting_windows {
        window.entries.reverse();
    }
    permuted.fidelity.baselines.reverse();
    permuted.fidelity.rungs.reverse();
    for rung in &mut permuted.fidelity.rungs {
        rung.falsifiers.reverse();
        rung.qois.reverse();
    }
    permuted.fidelity.escalations.reverse();
    let actual = permuted
        .admit_against(&graph, &behavior, &[case.clone()])
        .expect("permuted assurance admits");
    assert_eq!(actual.identity(), expected.identity());
    assert_eq!(actual.identity_receipt(), expected.identity_receipt());
}

#[test]
fn g3_behavior_and_policy_artifacts_move_identity_independently() {
    let (graph, behavior, case) = admitted_stack(11.0);
    let baseline = valid_assurance(&case)
        .admit_against(&graph, &behavior, &[case.clone()])
        .expect("baseline admits");

    let mut replay_changed = valid_assurance(&case);
    replay_changed.fidelity.fixed_replay = aref!(FixedReplayRef, "replay/alternate", 140);
    let replay_changed = replay_changed
        .admit_against(&graph, &behavior, &[case.clone()])
        .expect("alternate replay admits");
    assert_ne!(replay_changed.identity(), baseline.identity());

    let mut trigger_changed = valid_assurance(&case);
    trigger_changed.fidelity.escalations[0].trigger =
        aref!(EscalationTriggerRef, "triggers/alternate", 141);
    let trigger_changed = trigger_changed
        .admit_against(&graph, &behavior, &[case.clone()])
        .expect("alternate trigger admits");
    assert_ne!(trigger_changed.identity(), baseline.identity());

    let mut changed_graph_draft = valid_graph();
    changed_graph_draft.subsystems[0].model = model(142);
    let changed_graph = changed_graph_draft.admit().expect("changed graph admits");
    let changed_behavior = valid_behavior()
        .admit_against(&changed_graph)
        .expect("changed behavior admits");
    let refusal = valid_assurance(&case)
        .admit_against(&graph, &changed_behavior, &[case.clone()])
        .expect_err("behavior from a different graph refuses");
    assert!(rules(&refusal).contains(&MachineAssuranceRule::BehaviorGraphMismatch));
}

fn list_items_mut(node: &mut Node) -> &mut Vec<Node> {
    let NodeKind::List(items) = &mut node.kind else {
        panic!("fixture node must be a list")
    };
    items
}

#[test]
fn g0_assurance_codec_roundtrips_versioned_sexpr_and_json_with_live_authority() {
    let (graph, behavior, case) = admitted_stack(11.0);
    let admitted = valid_assurance(&case)
        .admit_against(&graph, &behavior, &[case.clone()])
        .expect("assurance fixture admits");
    let program = write_machine_assurance_program_v1(&admitted).expect("assurance encodes");
    let sexpr = program
        .print_sexpr_checked()
        .expect("canonical assurance s-expression");
    let json = program
        .print_json_checked()
        .expect("canonical assurance JSON");
    let from_sexpr = VersionedProgram::parse_sexpr(&sexpr).expect("s-expression reparses");
    let from_json = VersionedProgram::parse_json(&json).expect("JSON reparses");

    for artifact in [&from_sexpr, &from_json] {
        let decoded =
            parse_machine_assurance_program_v1(artifact).expect("assurance program decodes");
        assert_eq!(decoded.base_graph(), graph.identity());
        assert_eq!(decoded.base_behavior(), behavior.identity());
        let readmitted = decoded
            .admit_against(&graph, &behavior, &[case.clone()])
            .expect("decoded assurance admits with the live V&V case");
        assert_eq!(readmitted.identity(), admitted.identity());
        assert!(
            write_machine_assurance_v1(&readmitted)
                .expect("readmitted assurance re-encodes")
                .same_shape(artifact.program())
        );
    }
}

#[test]
fn g0_serialized_vv_commitment_cannot_mint_assurance_authority() {
    let (graph, behavior, case) = admitted_stack(11.0);
    let admitted = valid_assurance(&case)
        .admit_against(&graph, &behavior, &[case.clone()])
        .expect("assurance fixture admits");
    let node = write_machine_assurance_v1(&admitted).expect("assurance encodes");
    let error = parse_machine_assurance_v1(&node)
        .expect("assurance syntax decodes")
        .admit_against(&graph, &behavior, &[])
        .expect_err("serialized commitments are not admitted V&V authority");
    let MachineAssuranceAstAdmissionError::Assurance(refusal) = error else {
        panic!("missing live authority must reach semantic assurance admission")
    };
    assert!(rules(&refusal).contains(&MachineAssuranceRule::MissingVvCase));
}

#[test]
fn g0_changed_missing_and_duplicate_vv_commitments_refuse_after_real_admission() {
    let (graph, behavior, case) = admitted_stack(11.0);
    let admitted = valid_assurance(&case)
        .admit_against(&graph, &behavior, &[case.clone()])
        .expect("assurance fixture admits");
    let canonical = write_machine_assurance_v1(&admitted).expect("assurance encodes");

    let mut changed = canonical.clone();
    let root = list_items_mut(&mut changed);
    let commitments = list_items_mut(&mut root[3]);
    let binding = list_items_mut(&mut commitments[1]);
    binding[2] = Node::synthetic(NodeKind::Str("99".to_owned()));
    let changed_error =
        admit_machine_assurance_ast_v1(&changed, &graph, &behavior, &[case.clone()])
            .expect_err("changed schema commitment refuses");
    assert!(matches!(
        changed_error,
        MachineAssuranceAstAdmissionError::VvCaseBindingMismatch { .. }
    ));

    let mut missing = canonical.clone();
    list_items_mut(&mut list_items_mut(&mut missing)[3]).truncate(1);
    let missing_error =
        admit_machine_assurance_ast_v1(&missing, &graph, &behavior, &[case.clone()])
            .expect_err("missing commitment refuses");
    assert!(matches!(
        missing_error,
        MachineAssuranceAstAdmissionError::VvCaseBindingMismatch { .. }
    ));

    let mut duplicate = canonical;
    let commitments = list_items_mut(&mut list_items_mut(&mut duplicate)[3]);
    commitments.push(commitments[1].clone());
    let duplicate_error =
        admit_machine_assurance_ast_v1(&duplicate, &graph, &behavior, &[case.clone()])
            .expect_err("duplicate commitment refuses");
    assert!(matches!(
        duplicate_error,
        MachineAssuranceAstAdmissionError::VvCaseBindingMismatch { .. }
    ));
}

#[test]
fn g0_assurance_base_binding_refuses_before_invalid_draft_semantics() {
    let (graph, behavior, case) = admitted_stack(11.0);
    let admitted = valid_assurance(&case)
        .admit_against(&graph, &behavior, &[case.clone()])
        .expect("assurance fixture admits");
    let mut node = write_machine_assurance_v1(&admitted).expect("assurance encodes");
    let root = list_items_mut(&mut node);
    list_items_mut(&mut root[1])[1] =
        Node::synthetic(NodeKind::Str(format!("{:02x}", 0x55).repeat(32)));
    list_items_mut(&mut root[4]).truncate(1);

    let error = admit_machine_assurance_ast_v1(&node, &graph, &behavior, &[case])
        .expect_err("wrong base graph refuses before empty sensors");
    assert!(matches!(
        error,
        MachineAssuranceAstAdmissionError::BaseGraphMismatch { .. }
    ));
}

#[test]
fn g3_rich_assurance_codec_preserves_every_nested_permutation_identity() {
    let (graph, behavior, case) = admitted_stack(11.0);
    let admitted = permutation_assurance(&case)
        .admit_against(&graph, &behavior, &[case.clone()])
        .expect("rich assurance admits");
    let node = write_machine_assurance_v1(&admitted).expect("rich assurance encodes");
    let readmitted = parse_machine_assurance_v1(&node)
        .expect("rich assurance decodes")
        .admit_against(&graph, &behavior, &[case])
        .expect("rich assurance readmits");
    assert_eq!(readmitted.identity(), admitted.identity());
    assert_eq!(readmitted.identity_receipt(), admitted.identity_receipt());
}

#[test]
fn g0_assurance_codec_reports_exact_digest_path_and_resource_first_precedence() {
    let (graph, behavior, case) = admitted_stack(11.0);
    let admitted = valid_assurance(&case)
        .admit_against(&graph, &behavior, &[case])
        .expect("assurance fixture admits");
    let canonical = write_machine_assurance_v1(&admitted).expect("assurance encodes");

    let mut uppercase = canonical.clone();
    let root = list_items_mut(&mut uppercase);
    let commitments = list_items_mut(&mut root[3]);
    let binding = list_items_mut(&mut commitments[1]);
    binding[4] = Node::synthetic(NodeKind::Str("AA".repeat(32)));
    let reference_error = parse_machine_assurance_v1(&uppercase)
        .expect_err("uppercase content hash is not canonical assurance syntax");
    assert_eq!(
        reference_error.rule(),
        MachineAssuranceCodecRule::InvalidReference
    );
    assert_eq!(reference_error.path(), "$[3][1][4]");

    let mut oversized = canonical;
    let commitments = list_items_mut(&mut list_items_mut(&mut oversized)[3]);
    let repeated = commitments[1].clone();
    commitments.extend((1..MAX_MACHINE_ASSURANCE_VV_CASE_RECEIPTS).map(|_| repeated.clone()));
    commitments.push(Node::synthetic(NodeKind::Float(f64::NAN)));
    let resource_error = parse_machine_assurance_v1(&oversized)
        .expect_err("public receipt cap must precede recursive AST validation");
    assert_eq!(
        resource_error.rule(),
        MachineAssuranceCodecRule::ResourceLimit
    );
    assert_eq!(resource_error.path(), "$[3]");
}
