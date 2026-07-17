//! Shared admitted Machine graph/behavior/assurance fixture for PR-5 tests.

use core::num::NonZeroU64;

use fs_blake3::ContentHash;
use fs_evidence::vv::*;
use fs_ir::machine::assurance::*;
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

pub(crate) fn admitted_machine_stack() -> (
    AdmittedMachineGraph,
    AdmittedMachineBehavior,
    AdmittedMachineAssurance,
) {
    let graph = valid_graph().admit().expect("machine graph admits");
    let behavior = valid_behavior()
        .admit_against(&graph)
        .expect("machine behavior admits");
    let case = admitted_vv_case(11.0);
    let assurance = valid_assurance(&case)
        .admit_against(&graph, &behavior, &[case])
        .expect("machine assurance admits");
    (graph, behavior, assurance)
}
