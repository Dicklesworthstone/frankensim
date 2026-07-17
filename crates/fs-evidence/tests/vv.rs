use fs_blake3::ContentHash;
use fs_evidence::vv::*;

fn hash(label: &str) -> ContentHash {
    fs_blake3::hash_domain("org.frankensim.fs-evidence.vv-test.v1", label.as_bytes())
}

fn push_identity_string(bytes: &mut Vec<u8>, value: &str) {
    bytes.extend_from_slice(&(value.len() as u64).to_le_bytes());
    bytes.extend_from_slice(value.as_bytes());
}

fn blind_holdout_identity_for(
    domain: &str,
    preregistration_hash: ContentHash,
    rows_in_wire_order: &[(ObservationId, ContentHash)],
) -> ContentHash {
    let mut bytes = Vec::new();
    bytes.extend_from_slice(preregistration_hash.as_bytes());
    bytes.extend_from_slice(&(rows_in_wire_order.len() as u64).to_le_bytes());
    for (id, source) in rows_in_wire_order {
        push_identity_string(&mut bytes, id.as_str());
        bytes.extend_from_slice(source.as_bytes());
    }
    fs_blake3::hash_domain(domain, &bytes)
}

fn schema_admission_receipt_preimage_in_wire_order(
    schema_version: u32,
    ruleset_version: u32,
    case_hash: ContentHash,
    context_id: &ArtifactId,
    qois_in_wire_order: &[QoiId],
    artifact_hashes_in_wire_order: &[((ArtifactKind, ArtifactId), ContentHash)],
) -> Vec<u8> {
    let mut bytes = Vec::new();
    bytes.extend_from_slice(&schema_version.to_le_bytes());
    bytes.extend_from_slice(&ruleset_version.to_le_bytes());
    bytes.extend_from_slice(case_hash.as_bytes());
    push_identity_string(&mut bytes, context_id.as_str());
    bytes.extend_from_slice(&(qois_in_wire_order.len() as u64).to_le_bytes());
    for qoi in qois_in_wire_order {
        push_identity_string(&mut bytes, qoi.as_str());
    }
    bytes.extend_from_slice(&(artifact_hashes_in_wire_order.len() as u64).to_le_bytes());
    for ((kind, id), hash) in artifact_hashes_in_wire_order {
        bytes.push(kind.canonical_wire_tag());
        push_identity_string(&mut bytes, kind.slug());
        push_identity_string(&mut bytes, id.as_str());
        bytes.extend_from_slice(hash.as_bytes());
    }
    bytes
}

fn schema_admission_receipt_preimage(
    schema_version: u32,
    ruleset_version: u32,
    case_hash: ContentHash,
    context_id: &ArtifactId,
    qois_in_wire_order: &[QoiId],
    artifact_hashes: &[((ArtifactKind, ArtifactId), ContentHash)],
) -> Vec<u8> {
    let mut artifact_hashes_in_wire_order = artifact_hashes.to_vec();
    artifact_hashes_in_wire_order.sort_by(
        |((left_kind, left_id), _), ((right_kind, right_id), _)| {
            left_kind
                .canonical_wire_tag()
                .cmp(&right_kind.canonical_wire_tag())
                .then_with(|| left_kind.slug().cmp(right_kind.slug()))
                .then_with(|| left_id.cmp(right_id))
        },
    );
    schema_admission_receipt_preimage_in_wire_order(
        schema_version,
        ruleset_version,
        case_hash,
        context_id,
        qois_in_wire_order,
        &artifact_hashes_in_wire_order,
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

fn observation_source_ref(source: &str) -> ObservationSourceRef {
    observation_source_ref_with(
        "source-bytes",
        "org.frankensim.fs-evidence.test-row-locator.v1",
        1,
        source,
        &format!("extraction-receipt-{source}"),
    )
}

fn observation_source_ref_with(
    dataset_source_bytes: &str,
    locator_domain: &str,
    locator_contract_version: u32,
    locator: &str,
    extraction_receipt: &str,
) -> ObservationSourceRef {
    ObservationSourceRef::try_new(
        hash(dataset_source_bytes),
        locator_domain,
        locator_contract_version,
        hash(locator),
        hash(extraction_receipt),
    )
    .expect("valid typed observation-source reference")
}

fn manifest_row(source: &str) -> ObservationManifestRow {
    manifest_row_with(source, "length", "instrument-1", "channel-1", "clock-1")
}

fn manifest_row_with(
    source: &str,
    qoi: &str,
    instrument: &str,
    acquisition_channel: &str,
    clock: &str,
) -> ObservationManifestRow {
    manifest_row_with_source(
        observation_source_ref(source),
        qoi,
        instrument,
        acquisition_channel,
        clock,
    )
}

fn manifest_row_with_source(
    source: ObservationSourceRef,
    qoi: &str,
    instrument: &str,
    acquisition_channel: &str,
    clock: &str,
) -> ObservationManifestRow {
    ObservationManifestRow::try_new(
        source,
        qoi_id(qoi),
        artifact_id(instrument),
        artifact_id(acquisition_channel),
        artifact_id(clock),
    )
    .expect("valid typed observation-manifest row")
}

fn axis_id(value: &str) -> AxisId {
    AxisId::try_new(value).expect("valid fixture axis id")
}

fn unit_id(value: &str) -> UnitId {
    UnitId::try_new(value).expect("valid fixture unit id")
}

fn header(id: &str, unit: &str) -> ArtifactHeader {
    header_with_units(id, &[unit])
}

fn header_with_units(id: &str, units: &[&str]) -> ArtifactHeader {
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

fn reference(kind: ArtifactKind, id: &str) -> ArtifactRef {
    ArtifactRef::new(kind, artifact_id(id), hash(id))
}

fn external_target(label: &str) -> EvidenceTarget {
    EvidenceTarget::External {
        family: artifact_id("fixture-family"),
        id: artifact_id(label),
        hash: hash(label),
    }
}

fn assert_rule<T>(result: Result<T, VvErrors>, expected: VvRule) {
    let error = match result {
        Ok(_) => panic!("expected {} refusal", expected.slug()),
        Err(error) => error,
    };
    assert!(
        error
            .violations()
            .iter()
            .any(|violation| violation.rule() == expected),
        "expected {}, got {error}",
        expected.slug(),
    );
}

fn assert_only_rule_field<T>(result: Result<T, VvErrors>, rule: VvRule, field: &str) {
    let error = match result {
        Ok(_) => panic!("expected {} refusal at {field}", rule.slug()),
        Err(error) => error,
    };
    assert_eq!(
        error.violations().len(),
        1,
        "expected one isolated {} refusal at {field}, got {error}",
        rule.slug(),
    );
    let violation = &error.violations()[0];
    assert_eq!(violation.rule(), rule, "unexpected refusal rule: {error}");
    assert_eq!(
        violation.field(),
        field,
        "unexpected refusal field: {error}"
    );
}

fn covariance_1() -> CovarianceMatrix {
    CovarianceMatrix::try_new(1, vec![0.25]).expect("positive scalar covariance")
}

fn experiment(
    id: &str,
    origin: ExperimentOrigin,
    calibration_current: bool,
    authenticated: bool,
) -> Result<ExperimentArtifact, VvErrors> {
    ExperimentArtifact::try_new(
        header(id, "m"),
        artifact_id(&format!("{id}-dataset")),
        origin,
        vec![qoi_id("length")],
        ObservationManifest::try_new(vec![
            (observation_id("cal-1"), manifest_row("row-cal-1")),
            (observation_id("val-1"), manifest_row("row-val-1")),
            (observation_id("blind-1"), manifest_row("row-blind-1")),
        ])
        .expect("injective manifest fixture"),
        vec![InstrumentCalibration::new(
            artifact_id("instrument-1"),
            hash("instrument-calibration"),
            calibration_current,
        )],
        ClockSynchronization::SingleClock {
            clock_id: artifact_id("clock-1"),
        },
        RepeatabilitySummary::try_new(3, covariance_1()).expect("repeatability fixture"),
        DataAuthenticity::new(hash("source-bytes"), hash("custody"), authenticated),
    )
}

fn experiment_with_manifest(
    manifest: ObservationManifest,
    qois: Vec<QoiId>,
    instruments: Vec<InstrumentCalibration>,
    clocks: ClockSynchronization,
) -> Result<ExperimentArtifact, VvErrors> {
    ExperimentArtifact::try_new(
        header("manifest-experiment", "m"),
        artifact_id("manifest-dataset"),
        ExperimentOrigin::Physical {
            apparatus_id: artifact_id("apparatus-1"),
            facility_id: artifact_id("facility-1"),
        },
        qois,
        manifest,
        instruments,
        clocks,
        RepeatabilitySummary::try_new(3, covariance_1()).expect("repeatability fixture"),
        DataAuthenticity::new(hash("source-bytes"), hash("custody"), true),
    )
}

fn experiment_with_authenticity(
    authenticity: DataAuthenticity,
) -> Result<ExperimentArtifact, VvErrors> {
    ExperimentArtifact::try_new(
        header("authenticity-experiment", "m"),
        artifact_id("authenticity-dataset"),
        ExperimentOrigin::Physical {
            apparatus_id: artifact_id("apparatus-1"),
            facility_id: artifact_id("facility-1"),
        },
        vec![qoi_id("length")],
        ObservationManifest::try_new(vec![(observation_id("row-1"), manifest_row("row-1"))])
            .expect("authenticity manifest fixture"),
        vec![InstrumentCalibration::new(
            artifact_id("instrument-1"),
            hash("instrument-calibration"),
            true,
        )],
        ClockSynchronization::SingleClock {
            clock_id: artifact_id("clock-1"),
        },
        RepeatabilitySummary::try_new(3, covariance_1()).expect("repeatability fixture"),
        authenticity,
    )
}

fn split_with(
    calibration: Vec<&str>,
    validation: Vec<&str>,
    blind: Vec<&str>,
) -> Result<CalibrationSplit, VvErrors> {
    CalibrationSplit::try_new(
        header("split-1", "unitless"),
        reference(ArtifactKind::ExperimentArtifact, "experiment-1"),
        hash("preregistration"),
        calibration.into_iter().map(observation_id).collect(),
        validation.into_iter().map(observation_id).collect(),
        blind
            .into_iter()
            .map(|id| (observation_id(id), hash(&format!("row-{id}"))))
            .collect(),
    )
}

fn uncertainty_terms(magnitudes: [f64; 6]) -> Vec<UncertaintyTerm> {
    PredictionUncertaintyKind::ALL
        .into_iter()
        .zip(magnitudes)
        .map(|(kind, magnitude)| {
            UncertaintyTerm::try_new(kind, magnitude, external_target(&format!("{kind:?}")))
                .expect("valid uncertainty term")
        })
        .collect()
}

fn categorical_axes() -> EvidenceAxes {
    EvidenceAxes::try_new(
        EvidenceAxis::ALL
            .into_iter()
            .map(|axis| {
                (
                    axis,
                    EvidenceAxisStatus::Missing {
                        reason: "fixture explicitly makes no positive claim".to_owned(),
                    },
                )
            })
            .collect(),
    )
    .expect("complete categorical axes")
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

fn rebuild_prediction(
    base: &PredictionAssessment,
    context: ArtifactRef,
    validation_plan: ArtifactRef,
    dependencies: Vec<EvidenceDependency>,
) -> PredictionAssessment {
    PredictionAssessment::try_new(
        base.header().clone(),
        context,
        validation_plan,
        base.qoi().clone(),
        dependencies,
        base.waterfall().clone(),
        base.validation_metrics().to_vec(),
        base.posterior_checks().to_vec(),
        base.applicability_point().clone(),
        base.applicability().clone(),
        base.evidence_axes().clone(),
        base.assumption_checks()
            .iter()
            .map(|(id, passed)| (id.clone(), *passed))
            .collect(),
    )
    .expect("rebuilt prediction fixture")
}

fn rebuild_prediction_with_selection(
    base: &PredictionAssessment,
    observations: ObservationSelection,
) -> PredictionAssessment {
    let dependencies = base
        .dependencies()
        .iter()
        .map(|dependency| {
            if dependency.role() != DependencyRole::PhysicalValidation {
                return dependency.clone();
            }
            let EvidenceTarget::VvArtifact(experiment) = dependency.target() else {
                panic!("physical validation fixture must target an experiment artifact");
            };
            EvidenceDependency::physical_validation(
                dependency.qoi().clone(),
                experiment.clone(),
                observations.clone(),
            )
        })
        .collect();
    let validation_metrics = base
        .validation_metrics()
        .iter()
        .map(|metric| {
            ValidationMetric::try_new(
                metric.name().clone(),
                metric.qoi().clone(),
                observations.clone(),
                metric.observed(),
                metric.predicted(),
                metric.experimental_uncertainty(),
                metric.numerical_uncertainty(),
            )
            .expect("selection-rebound validation metric")
        })
        .collect();
    let posterior_checks = base
        .posterior_checks()
        .iter()
        .map(|check| {
            PosteriorPredictiveCheck::try_new(
                check.name().clone(),
                check.qoi().clone(),
                observations.clone(),
                check.tail_probability(),
                check.minimum_tail_probability(),
                check.artifact_hash(),
            )
            .expect("selection-rebound posterior check")
        })
        .collect();
    PredictionAssessment::try_new(
        base.header().clone(),
        base.context().clone(),
        base.validation_plan().clone(),
        base.qoi().clone(),
        dependencies,
        base.waterfall().clone(),
        validation_metrics,
        posterior_checks,
        base.applicability_point().clone(),
        base.applicability().clone(),
        base.evidence_axes().clone(),
        base.assumption_checks()
            .iter()
            .map(|(id, passed)| (id.clone(), *passed))
            .collect(),
    )
    .expect("prediction with rebound observation selection")
}

fn case_replacing_prediction(base: &VvCase, prediction: PredictionAssessment) -> VvCase {
    VvCase::try_new(
        base.context().clone(),
        base.validation_plan().clone(),
        base.experiments().values().cloned().collect(),
        base.splits().values().cloned().collect(),
        base.solution_verification().values().cloned().collect(),
        vec![prediction],
        base.assumptions().clone(),
    )
    .expect("case with one replaced prediction fixture")
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum OriginKind {
    Physical,
    SyntheticHighFidelity,
    SecondImplementation,
}

#[derive(Clone, Copy, Debug)]
struct CaseKnobs {
    origin: OriginKind,
    failing_diagnostic: Option<VvRule>,
    outside_domain: bool,
    context_domain_lo: f64,
    reported_domain_lo: Option<f64>,
    failing_assumption: Option<&'static str>,
    force_in_domain_claim: bool,
    process_as_uncertainty: bool,
    tamper_seed: Option<&'static str>,
    predicted: f64,
    posterior_minimum: f64,
    posterior_tail: f64,
    blind_source_label: &'static str,
}

impl Default for CaseKnobs {
    fn default() -> Self {
        Self {
            origin: OriginKind::Physical,
            failing_diagnostic: None,
            outside_domain: false,
            context_domain_lo: 1.0,
            reported_domain_lo: None,
            failing_assumption: None,
            force_in_domain_claim: false,
            process_as_uncertainty: false,
            tamper_seed: None,
            // Genuinely agrees: |9.9 - 9.95| = 0.05 within the combined
            // uncertainty (~0.09); bead gt1k3 now DERIVES this outcome,
            // so a laundering fixture would refuse.
            predicted: 9.95,
            posterior_minimum: 0.05,
            posterior_tail: 0.5,
            blind_source_label: "row-blind-1",
        }
    }
}

fn diagnostic_record(rule: VvRule, failing: Option<VvRule>) -> DiagnosticRecord {
    let passed = failing != Some(rule);
    DiagnosticRecord::try_new(
        passed,
        hash(rule.slug()),
        if passed {
            format!("{} passed independently", rule.slug())
        } else {
            format!("{} was not established", rule.slug())
        },
    )
    .expect("diagnostic fixture")
}

fn diagnostic_plan(failing: Option<VvRule>) -> DiagnosticPlan {
    DiagnosticPlan::new(
        diagnostic_record(VvRule::DiagnosticObservability, failing),
        diagnostic_record(VvRule::DiagnosticIdentifiability, failing),
        diagnostic_record(VvRule::DiagnosticConfounding, failing),
        diagnostic_record(VvRule::DiagnosticInverseCrime, failing),
    )
}

fn uncertainty_term(
    kind: PredictionUncertaintyKind,
    magnitude: f64,
    source: EvidenceTarget,
) -> UncertaintyTerm {
    UncertaintyTerm::try_new(kind, magnitude, source).expect("valid closed-case uncertainty term")
}

fn closed_case(knobs: CaseKnobs) -> VvCase {
    let qoi = qoi_id("length");
    let unit = unit_id("m");
    let re_axis = axis_id("Re");
    let regime_axis = axis_id("regime");
    let applicability = ApplicabilityDomain::try_new(
        vec![
            NumericDomainAxis::try_new(
                re_axis.clone(),
                unit_id("unitless"),
                knobs.context_domain_lo,
                10.0,
            )
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
    .expect("mixed applicability fixture");
    let context = ContextOfUse::try_new(
        header_with_units("context-1", &["m", "unitless"]),
        "Decide whether the measured length satisfies the release criterion.",
        vec![
            QoiSpec::try_new(
                qoi.clone(),
                "released length",
                unit.clone(),
                AcceptanceCriterion::ClosedRange { lo: 9.0, hi: 11.0 },
            )
            .expect("QoI fixture"),
        ],
        applicability,
        ApplicabilityPolicy::Demote,
    )
    .expect("context fixture");
    let context_ref = artifact_reference(&context);

    let origin = match knobs.origin {
        OriginKind::Physical => ExperimentOrigin::Physical {
            apparatus_id: artifact_id("apparatus-1"),
            facility_id: artifact_id("facility-1"),
        },
        OriginKind::SyntheticHighFidelity => ExperimentOrigin::SyntheticHighFidelity {
            producer: artifact_id("high-fidelity-code"),
        },
        OriginKind::SecondImplementation => ExperimentOrigin::SecondImplementation {
            producer: artifact_id("second-implementation"),
        },
    };
    let experiment = experiment("experiment-1", origin, true, true).expect("experiment fixture");
    let experiment_ref = artifact_reference(&experiment);

    let split = CalibrationSplit::try_new(
        header("split-1", "unitless"),
        experiment_ref.clone(),
        hash("preregistered-analysis"),
        vec![observation_id("cal-1")],
        vec![observation_id("val-1")],
        vec![(observation_id("blind-1"), hash(knobs.blind_source_label))],
    )
    .expect("closed-case split");
    let split_ref = artifact_reference(&split);

    let plan_row = QoiValidationPlan::try_new(
        qoi.clone(),
        vec![experiment_ref.clone()],
        split_ref.clone(),
        vec![
            ValidationMetricSpec::IntervalAgreement,
            ValidationMetricSpec::PosteriorPredictive {
                minimum_tail_probability: 0.05,
            },
        ],
        diagnostic_plan(knobs.failing_diagnostic),
    )
    .expect("QoI validation-plan row");
    let validation_plan = ValidationPlan::try_new(
        header("validation-plan-1", "m"),
        context_ref.clone(),
        vec![plan_row],
    )
    .expect("validation plan fixture");
    let validation_plan_ref = artifact_reference(&validation_plan);

    let numerical = |label| {
        NumericalUncertainty::try_new(0.01, hash(label)).expect("solution uncertainty fixture")
    };
    let solution = SolutionVerificationReceipt::try_new(
        header("solution-1", "m"),
        artifact_id("solve-1"),
        qoi.clone(),
        unit.clone(),
        numerical("mesh-bound"),
        numerical("time-bound"),
        numerical("nonlinear-bound"),
        numerical("iterative-bound"),
    )
    .expect("solution-verification fixture");
    let numerical_floor = solution.combined_half_width();
    let solution_ref = artifact_reference(&solution);

    let validation_selection = split
        .validation_selection(split_ref.clone(), vec![observation_id("val-1")])
        .expect("validation selection fixture");
    let physical_dependency = EvidenceDependency::physical_validation(
        qoi.clone(),
        experiment_ref.clone(),
        validation_selection.clone(),
    );
    let solution_source = EvidenceTarget::VvArtifact(solution_ref.clone());
    let model_source = external_target(if knobs.process_as_uncertainty {
        "process-conformance"
    } else {
        "model-discrepancy"
    });
    let parameter_source = external_target("parameter-data");
    let data_source = external_target("measurement-data");
    let aleatory_source = external_target("aleatory-model");
    let epistemic_source = external_target("epistemic-model");

    let mut dependencies = vec![
        physical_dependency,
        EvidenceDependency::new(
            qoi.clone(),
            DependencyRole::SolutionVerification,
            solution_source.clone(),
        ),
        EvidenceDependency::new(
            qoi.clone(),
            if knobs.process_as_uncertainty {
                DependencyRole::ProcessConformance
            } else {
                DependencyRole::ModelDiscrepancy
            },
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
        unit.clone(),
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
    .expect("closed-case uncertainty waterfall");
    let metric = ValidationMetric::try_new(
        artifact_id("interval-agreement"),
        qoi.clone(),
        validation_selection.clone(),
        9.9,
        knobs.predicted,
        0.05,
        numerical_floor,
    )
    .expect("validation metric fixture");
    let posterior = PosteriorPredictiveCheck::try_new(
        artifact_id("posterior-check"),
        qoi.clone(),
        validation_selection,
        knobs.posterior_tail,
        knobs.posterior_minimum,
        hash("posterior-check-artifact"),
    )
    .expect("posterior-predictive fixture");

    let mut assumptions = AssumptionsLedger::try_program_seed(header("assumptions", "unitless"))
        .expect("program assumption fixture");
    let retained_rows = assumptions.rows().values().cloned().collect::<Vec<_>>();
    for row in retained_rows {
        let label = format!("{}-evidence", row.id().as_str());
        assumptions
            .replace_row(
                row.with_evidence(external_target(&label))
                    .with_monitor_evidence(hash(&format!("{}-monitor", label))),
            )
            .expect("attach retained assumption evidence");
    }
    if let Some(id) = knobs.tamper_seed {
        let id = AssumptionId::try_new(id).expect("tampered seed id");
        let original = assumptions
            .rows()
            .get(&id)
            .expect("seed row to tamper")
            .clone();
        let tampered = AssumptionRow::try_new(
            original.id().clone(),
            format!("Tampered: {}", original.predicate()),
            original.scope(),
            original.evidence().clone(),
            original.monitor().clone(),
            original.violation_effect().clone(),
            original.owner().clone(),
            original.review_gate().clone(),
        )
        .expect("structurally valid tampered seed row");
        assumptions
            .replace_row(tampered)
            .expect("replace existing seed row");
    }

    let applicability_point = ApplicabilityPoint::try_new(
        vec![(
            re_axis.clone(),
            if knobs.outside_domain { 20.0 } else { 5.0 },
        )],
        vec![(regime_axis, "nominal".to_owned())],
    )
    .expect("applicability point fixture");
    let assumption_checks = assumptions
        .rows()
        .keys()
        .cloned()
        .map(|id| {
            let passed = knobs.failing_assumption != Some(id.as_str());
            (id, passed)
        })
        .collect::<Vec<_>>();
    let mut domain_violations = Vec::new();
    if knobs.outside_domain {
        domain_violations.push(DomainViolation::Numeric {
            axis: re_axis,
            value: 20.0,
            lo: knobs.reported_domain_lo.unwrap_or(knobs.context_domain_lo),
            hi: 10.0,
        });
    }
    let assumption_forces_refusal = knobs.failing_assumption.is_some_and(|id| {
        let id = AssumptionId::try_new(id).expect("failing assumption id");
        domain_violations.push(DomainViolation::Assumption { id: id.clone() });
        matches!(
            assumptions
                .rows()
                .get(&id)
                .expect("failing assumption row")
                .violation_effect(),
            ViolationEffect::EscalateOrRefuse { .. } | ViolationEffect::Refuse { .. }
        )
    });
    let honest_applicability = if domain_violations.is_empty() {
        ApplicabilityDecision::InDomain
    } else if assumption_forces_refusal {
        ApplicabilityDecision::Refused {
            violations: domain_violations,
        }
    } else {
        ApplicabilityDecision::Demoted {
            violations: domain_violations,
        }
    };
    let applicability_decision = if knobs.force_in_domain_claim {
        ApplicabilityDecision::InDomain
    } else {
        honest_applicability
    };

    let prediction = PredictionAssessment::try_new(
        header("prediction-1", "m"),
        context_ref,
        validation_plan_ref,
        qoi,
        dependencies,
        waterfall,
        vec![metric],
        vec![posterior],
        applicability_point,
        applicability_decision,
        categorical_axes(),
        assumption_checks,
    )
    .expect("prediction assessment fixture");

    VvCase::try_new(
        context,
        validation_plan,
        vec![experiment],
        vec![split],
        vec![solution],
        vec![prediction],
        assumptions,
    )
    .expect("closed V&V case fixture")
}

#[test]
fn vv_rule_slugs_are_stable_and_unique() {
    let expected = [
        (
            VvRule::SplitPartitionsDisjoint,
            "vv-split-partitions-disjoint",
        ),
        (
            VvRule::SplitBlindHoldoutSealed,
            "vv-split-blind-holdout-sealed",
        ),
        (VvRule::ColorCategoricalOnly, "vv-color-categorical-only"),
        (
            VvRule::ValidationRequiresPhysicalReferent,
            "vv-validation-requires-physical-referent",
        ),
        (VvRule::QoiDependencyClosed, "vv-qoi-dependency-closed"),
        (VvRule::QoiDependencyIsolated, "vv-qoi-dependency-isolated"),
        (VvRule::WaterfallArithmetic, "vv-waterfall-arithmetic"),
        (
            VvRule::WaterfallDependenceDeclared,
            "vv-waterfall-dependence-declared",
        ),
        (
            VvRule::ProcessConformanceSeparate,
            "vv-process-conformance-separate",
        ),
        (VvRule::AssumptionA001, "vv-assumption-a001"),
        (VvRule::AssumptionA008, "vv-assumption-a008"),
        (VvRule::ReceiptBinding, "vv-receipt-binding"),
    ];
    let mut slugs = std::collections::BTreeSet::new();
    for (rule, slug) in expected {
        assert_eq!(rule.slug(), slug);
        assert!(slugs.insert(slug), "duplicate rule slug {slug}");
    }
}

#[test]
fn calibration_validation_and_blind_partitions_are_pairwise_disjoint() {
    for (calibration, validation, blind) in [
        (vec!["shared"], vec!["shared"], vec!["blind"]),
        (vec!["shared"], vec!["validation"], vec!["shared"]),
        (vec!["calibration"], vec!["shared"], vec!["shared"]),
    ] {
        assert_rule(
            split_with(calibration, validation, blind),
            VvRule::SplitPartitionsDisjoint,
        );
    }
}

#[test]
fn blind_holdout_requires_exact_release_and_calibration_never_validates() {
    let split =
        split_with(vec!["cal-1"], vec!["val-1"], vec!["blind-1"]).expect("valid disjoint split");
    let split_ref = artifact_reference(&split);
    let forged_split_ref = ArtifactRef::new(
        ArtifactKind::CalibrationSplit,
        split.id().clone(),
        hash("forged-same-id-split-content"),
    );

    let forged_validation_error = split
        .validation_selection(forged_split_ref.clone(), vec![observation_id("val-1")])
        .expect_err("a same-id forged split hash cannot mint a validation capability");
    assert!(
        forged_validation_error
            .violations()
            .iter()
            .any(|violation| {
                violation.rule() == VvRule::SplitPartitionsDisjoint
                    && violation.field() == "selection.split"
            })
    );
    let forged_release = BlindReleaseReceipt::new(
        forged_split_ref.clone(),
        split.blind_commitment(),
        hash("release-authority"),
    )
    .expect("a release receipt remains structurally constructible before split admission");
    let forged_blind_error = split
        .blind_selection(
            forged_split_ref,
            vec![observation_id("blind-1")],
            forged_release,
        )
        .expect_err("a same-id forged split hash cannot mint a blind-release capability");
    assert!(forged_blind_error.violations().iter().any(|violation| {
        violation.rule() == VvRule::SplitPartitionsDisjoint
            && violation.field() == "selection.split"
    }));

    assert_rule(
        split.validation_selection(split_ref.clone(), vec![observation_id("cal-1")]),
        VvRule::ValidationRequiresPhysicalReferent,
    );

    let wrong_release = BlindReleaseReceipt::new(
        split_ref.clone(),
        hash("wrong-commitment"),
        hash("release-authority"),
    )
    .expect("structurally well-formed but wrongly bound release");
    assert_rule(
        split.blind_selection(
            split_ref.clone(),
            vec![observation_id("blind-1")],
            wrong_release,
        ),
        VvRule::SplitBlindHoldoutSealed,
    );

    let release = BlindReleaseReceipt::new(
        split_ref.clone(),
        split.blind_commitment(),
        hash("release-authority"),
    )
    .expect("valid release receipt");
    let selection = split
        .blind_selection(split_ref, vec![observation_id("blind-1")], release)
        .expect("correctly released blind row");
    assert!(matches!(
        selection.partition(),
        EvidencePartition::BlindHoldout { .. }
    ));
}

#[test]
fn released_blind_selection_roundtrips_and_release_fields_move_artifact_identity() {
    let case = closed_case(CaseKnobs::default());
    let split = case.splits().values().next().expect("closed-case split");
    let split_ref = artifact_reference(split);
    let prediction = case
        .predictions()
        .values()
        .next()
        .expect("closed-case prediction");

    let blind_prediction_for_authority = |authority: &str| {
        let release =
            BlindReleaseReceipt::new(split_ref.clone(), split.blind_commitment(), hash(authority))
                .expect("exact blind release fixture");
        let selection = split
            .blind_selection(split_ref.clone(), vec![observation_id("blind-1")], release)
            .expect("released blind selection fixture");
        rebuild_prediction_with_selection(prediction, selection)
    };

    let baseline_prediction = blind_prediction_for_authority("release-authority-a");
    case_replacing_prediction(&case, baseline_prediction.clone())
        .validate()
        .expect("an exactly released blind selection is valid in the complete case");
    let baseline: VvArtifact = baseline_prediction.into();
    let baseline_bytes = baseline.canonical_bytes().expect("blind prediction bytes");
    let decoded = VvArtifact::from_canonical_bytes(&baseline_bytes)
        .expect("safe released-blind artifact roundtrip");
    assert_eq!(decoded, baseline);
    assert_eq!(
        decoded.content_hash().expect("decoded blind identity"),
        baseline.content_hash().expect("baseline blind identity"),
    );

    let other_authority: VvArtifact = blind_prediction_for_authority("release-authority-b").into();
    assert_ne!(
        other_authority
            .canonical_bytes()
            .expect("other authority bytes"),
        baseline_bytes,
        "release authority is retained in canonical transport",
    );
    assert_ne!(
        other_authority
            .content_hash()
            .expect("other authority identity"),
        baseline.content_hash().expect("baseline blind identity"),
        "release authority is part of prediction-artifact identity",
    );

    let mut forged_commitment_bytes = baseline_bytes;
    let old_commitment = split.blind_commitment();
    let new_commitment = hash("forged-other-blind-commitment");
    let commitment_offsets = forged_commitment_bytes
        .windows(old_commitment.as_bytes().len())
        .enumerate()
        .filter_map(|(offset, window)| (window == old_commitment.as_bytes()).then_some(offset))
        .collect::<Vec<_>>();
    assert_eq!(
        commitment_offsets.len(),
        3,
        "fixture carries the release commitment once in its dependency, metric, and posterior check",
    );
    for offset in commitment_offsets {
        forged_commitment_bytes[offset..offset + new_commitment.as_bytes().len()]
            .copy_from_slice(new_commitment.as_bytes());
    }
    let forged_commitment_artifact = VvArtifact::from_canonical_bytes(&forged_commitment_bytes)
        .expect("standalone decoding binds but cannot authenticate the enclosing split commitment");
    assert_ne!(
        forged_commitment_artifact
            .content_hash()
            .expect("forged commitment artifact identity"),
        baseline.content_hash().expect("baseline blind identity"),
        "blind commitment is part of prediction-artifact identity",
    );
    let VvArtifact::PredictionAssessment(forged_prediction) = forged_commitment_artifact else {
        panic!("forged fixture remains a prediction assessment");
    };
    let error = case_replacing_prediction(&case, forged_prediction)
        .validate()
        .expect_err("whole-case validation must authenticate the release commitment");
    let sealed_fields = error
        .violations()
        .iter()
        .filter(|violation| violation.rule() == VvRule::SplitBlindHoldoutSealed)
        .map(VvViolation::field)
        .collect::<std::collections::BTreeSet<_>>();
    assert_eq!(
        sealed_fields,
        std::collections::BTreeSet::from([
            "prediction.physical_validation.observations",
            "prediction.posterior_checks.observations",
            "prediction.validation_metrics.observations",
        ]),
        "every consumer of the forged release must fail closed at its exact field",
    );
}

#[test]
fn stale_blind_split_hash_roundtrips_standalone_but_all_case_consumers_refuse() {
    let case = closed_case(CaseKnobs::default());
    let split = case.splits().values().next().expect("closed-case split");
    let split_ref = artifact_reference(split);
    let prediction = case
        .predictions()
        .values()
        .next()
        .expect("closed-case prediction");
    let release = BlindReleaseReceipt::new(
        split_ref.clone(),
        split.blind_commitment(),
        hash("release-authority"),
    )
    .expect("exact blind release fixture");
    let selection = split
        .blind_selection(split_ref.clone(), vec![observation_id("blind-1")], release)
        .expect("released blind selection fixture");
    let blind_prediction = rebuild_prediction_with_selection(prediction, selection);
    case_replacing_prediction(&case, blind_prediction.clone())
        .validate()
        .expect("baseline blind prediction is valid in the complete case");

    let mut stale_bytes = blind_prediction
        .canonical_bytes()
        .expect("blind prediction bytes");
    let live_split_hash = split_ref.hash();
    let stale_split_hash = hash("stale-same-id-split-content");
    assert_ne!(stale_split_hash, live_split_hash);
    let split_hash_offsets = stale_bytes
        .windows(live_split_hash.as_bytes().len())
        .enumerate()
        .filter_map(|(offset, window)| (window == live_split_hash.as_bytes()).then_some(offset))
        .collect::<Vec<_>>();
    assert_eq!(
        split_hash_offsets.len(),
        6,
        "each dependency, validation metric, and posterior check carries both the outer selection and nested release split reference",
    );
    for offset in split_hash_offsets {
        stale_bytes[offset..offset + stale_split_hash.as_bytes().len()]
            .copy_from_slice(stale_split_hash.as_bytes());
    }

    let stale_prediction = PredictionAssessment::from_canonical_bytes(&stale_bytes)
        .expect("standalone decode accepts mutually consistent stale outer and release split refs");
    let error = case_replacing_prediction(&case, stale_prediction)
        .validate()
        .expect_err("whole-case validation must authenticate every selection split reference");
    let stale_fields = error
        .violations()
        .iter()
        .filter(|violation| violation.rule() == VvRule::SplitPartitionsDisjoint)
        .map(VvViolation::field)
        .collect::<std::collections::BTreeSet<_>>();
    assert_eq!(
        stale_fields,
        std::collections::BTreeSet::from([
            "prediction.physical_validation.observations",
            "prediction.posterior_checks.observations",
            "prediction.validation_metrics.observations",
        ]),
        "dependency, metric, and posterior consumers must all refuse their stale split capability at the exact field",
    );
}

#[test]
fn lineage_debug_surfaces_are_bounded_and_redact_authority_material() {
    fn push_hash_fragments(fragments: &mut Vec<String>, hash: ContentHash) {
        fragments.push(hash.to_string());
        fragments.push(format!("{hash:?}"));
    }

    let case = closed_case(CaseKnobs::default());
    let experiment = case
        .experiments()
        .values()
        .next()
        .expect("closed-case experiment");
    let split = case.splits().values().next().expect("closed-case split");
    let blind_id = split
        .blind_sources()
        .keys()
        .next()
        .expect("closed-case blind row")
        .clone();
    let split_ref = artifact_reference(split);
    let release_authority = hash("debug-release-authority");
    let release = BlindReleaseReceipt::new(
        split_ref.clone(),
        split.blind_commitment(),
        release_authority,
    )
    .expect("exact blind release");
    let selection = split
        .blind_selection(split_ref, vec![blind_id.clone()], release.clone())
        .expect("exact blind selection");
    let admitted = case.clone().admit().expect("closed case admits");
    let receipt = admitted.receipt();

    let mut sensitive_fragments = vec![
        case.context().id().as_str().to_string(),
        experiment.id().as_str().to_string(),
        experiment.dataset_id().as_str().to_string(),
        split.id().as_str().to_string(),
        split.experiment().id().as_str().to_string(),
        selection.split().id().as_str().to_string(),
        release.split().id().as_str().to_string(),
        receipt.context_id().as_str().to_string(),
    ];
    let mut debug_surfaces = vec![
        ("manifest", format!("{:?}", experiment.manifest())),
        ("authenticity", format!("{:?}", experiment.authenticity())),
        ("experiment", format!("{experiment:?}")),
        ("split", format!("{split:?}")),
        ("release", format!("{release:?}")),
        ("selection", format!("{selection:?}")),
        ("partition", format!("{:?}", selection.partition())),
        ("case", format!("{case:?}")),
        ("receipt", format!("{receipt:?}")),
        ("admitted case", format!("{admitted:?}")),
    ];

    for (row_id, row) in experiment.manifest().rows() {
        let source = row.source_ref();
        let locator = source.locator_identity();
        sensitive_fragments.extend([
            row_id.as_str().to_string(),
            row.qoi().as_str().to_string(),
            row.instrument().as_str().to_string(),
            row.acquisition_channel().as_str().to_string(),
            row.clock().as_str().to_string(),
            source.locator_domain().to_string(),
        ]);
        for hash in [
            source.dataset_source_bytes_hash(),
            source.locator_hash(),
            source.extraction_receipt_hash(),
        ] {
            push_hash_fragments(&mut sensitive_fragments, hash);
        }
        debug_surfaces.extend([
            ("source", format!("{source:?}")),
            ("locator", format!("{locator:?}")),
            ("manifest row", format!("{row:?}")),
        ]);
    }

    match experiment.origin() {
        ExperimentOrigin::Physical {
            apparatus_id,
            facility_id,
        } => sensitive_fragments.extend([
            apparatus_id.as_str().to_string(),
            facility_id.as_str().to_string(),
        ]),
        ExperimentOrigin::SyntheticHighFidelity { producer }
        | ExperimentOrigin::SecondImplementation { producer } => {
            sensitive_fragments.push(producer.as_str().to_string());
        }
    }
    for instrument in experiment.instruments() {
        sensitive_fragments.push(instrument.instrument_id().as_str().to_string());
        push_hash_fragments(&mut sensitive_fragments, instrument.certificate_hash());
    }
    match experiment.clocks() {
        ClockSynchronization::SingleClock { clock_id } => {
            sensitive_fragments.push(clock_id.as_str().to_string());
        }
        ClockSynchronization::Synchronized {
            clock_ids,
            method,
            evidence_hash,
            ..
        } => {
            sensitive_fragments.extend(
                clock_ids
                    .iter()
                    .map(|clock_id| clock_id.as_str().to_string()),
            );
            sensitive_fragments.push(method.clone());
            push_hash_fragments(&mut sensitive_fragments, *evidence_hash);
        }
    }
    for hash in [
        experiment.authenticity().source_bytes_hash(),
        experiment.authenticity().custody_receipt_hash(),
        experiment.observations_hash(),
        split.experiment().hash(),
        split.preregistration_hash(),
        split.blind_commitment(),
        selection.split().hash(),
        release.split().hash(),
        release.blind_commitment(),
        release.authority_receipt_hash(),
        release_authority,
        receipt.case_hash(),
        receipt.receipt_hash(),
    ] {
        push_hash_fragments(&mut sensitive_fragments, hash);
    }
    sensitive_fragments.extend(
        split
            .calibration_ids()
            .iter()
            .chain(split.validation_ids())
            .map(|id| id.as_str().to_string()),
    );
    for (id, source_hash) in split.blind_sources() {
        sensitive_fragments.push(id.as_str().to_string());
        push_hash_fragments(&mut sensitive_fragments, *source_hash);
    }
    sensitive_fragments.extend(selection.ids().iter().map(|id| id.as_str().to_string()));
    sensitive_fragments.extend(receipt.qois().iter().map(|qoi| qoi.as_str().to_string()));
    for ((_, id), hash) in receipt.artifact_hashes() {
        sensitive_fragments.push(id.as_str().to_string());
        push_hash_fragments(&mut sensitive_fragments, *hash);
    }
    sensitive_fragments.sort();
    sensitive_fragments.dedup();

    for artifact in [
        VvArtifact::from(case.context().clone()),
        VvArtifact::from(case.validation_plan().clone()),
        VvArtifact::from(experiment.clone()),
        VvArtifact::from(split.clone()),
        VvArtifact::from(
            case.solution_verification()
                .values()
                .next()
                .expect("closed-case solution verification")
                .clone(),
        ),
        VvArtifact::from(
            case.predictions()
                .values()
                .next()
                .expect("closed-case prediction")
                .clone(),
        ),
        VvArtifact::from(case.assumptions().clone()),
    ] {
        debug_surfaces.push((artifact.kind().slug(), format!("{artifact:?}")));
    }
    for (surface, debug) in debug_surfaces {
        assert!(
            debug.len() < 1_024,
            "{surface} Debug must remain bounded: {} bytes",
            debug.len(),
        );
        for sensitive in &sensitive_fragments {
            assert!(
                !debug.contains(sensitive.as_str()),
                "{surface} Debug leaked `{sensitive}`: {debug}",
            );
        }
    }
    let representative_source = experiment
        .manifest()
        .rows()
        .values()
        .next()
        .expect("closed-case manifest source")
        .source_ref();
    assert!(format!("{representative_source:?}").contains("locator_contract_version"));
    assert!(format!("{:?}", experiment.manifest()).contains("row_count"));
    assert!(format!("{experiment:?}").contains("observation_count"));
    assert!(format!("{split:?}").contains("blind_holdout_count"));
    assert!(format!("{selection:?}").contains("observation_count"));
    assert!(format!("{case:?}").contains("experiment_count"));
    assert!(format!("{receipt:?}").contains("binding_present"));
    assert!(format!("{admitted:?}").contains("receipt_binding_present"));
}

#[test]
fn experiment_origin_and_traceability_remain_distinct() {
    let physical = experiment(
        "physical-experiment",
        ExperimentOrigin::Physical {
            apparatus_id: artifact_id("apparatus-1"),
            facility_id: artifact_id("facility-1"),
        },
        true,
        true,
    )
    .expect("physical experiment");
    assert!(physical.origin().is_physical());

    let synthetic = experiment(
        "synthetic-experiment",
        ExperimentOrigin::SyntheticHighFidelity {
            producer: artifact_id("reference-code"),
        },
        true,
        true,
    )
    .expect("synthetic discrepancy evidence remains representable");
    assert!(!synthetic.origin().is_physical());

    assert_rule(
        experiment(
            "stale-calibration",
            ExperimentOrigin::Physical {
                apparatus_id: artifact_id("apparatus-1"),
                facility_id: artifact_id("facility-1"),
            },
            false,
            true,
        ),
        VvRule::ExperimentInstrumentCalibration,
    );
    assert_rule(
        experiment(
            "unauthenticated-data",
            ExperimentOrigin::Physical {
                apparatus_id: artifact_id("apparatus-1"),
                facility_id: artifact_id("facility-1"),
            },
            true,
            false,
        ),
        VvRule::ExperimentDataAuthenticity,
    );
}

#[test]
fn experiment_authenticity_and_row_dataset_closure_refuse_sentinels_and_cross_wiring() {
    for (field, authenticity) in [
        (
            "experiment.authenticity.source_bytes_hash",
            DataAuthenticity::new(ContentHash([0; 32]), hash("custody"), true),
        ),
        (
            "experiment.authenticity.custody_receipt_hash",
            DataAuthenticity::new(hash("source-bytes"), ContentHash([0; 32]), true),
        ),
    ] {
        let error = experiment_with_authenticity(authenticity)
            .expect_err("all-zero provenance hashes must refuse direct construction");
        assert!(
            error.violations().iter().any(|violation| {
                violation.rule() == VvRule::ExperimentDataAuthenticity && violation.field() == field
            }),
            "expected exact authenticity field {field}, got {error}",
        );
    }

    let foreign_dataset_manifest = ObservationManifest::try_new(vec![(
        observation_id("row-1"),
        manifest_row_with_source(
            observation_source_ref_with(
                "foreign-source-bytes",
                "org.frankensim.fs-evidence.test-row-locator.v1",
                1,
                "row-1",
                "extraction-receipt-row-1",
            ),
            "length",
            "instrument-1",
            "channel-1",
            "clock-1",
        ),
    )])
    .expect("foreign dataset source remains structurally typed");
    let error = experiment_with_manifest(
        foreign_dataset_manifest,
        vec![qoi_id("length")],
        vec![InstrumentCalibration::new(
            artifact_id("instrument-1"),
            hash("instrument-calibration"),
            true,
        )],
        ClockSynchronization::SingleClock {
            clock_id: artifact_id("clock-1"),
        },
    )
    .expect_err("a row from different source bytes must refuse experiment closure");
    assert!(error.violations().iter().any(|violation| {
        violation.rule() == VvRule::ExperimentDataAuthenticity
            && violation.field() == "experiment.manifest.dataset_source_bytes_hash"
    }));
}

#[test]
fn experiment_codec_refuses_forged_authenticity_and_typed_row_sources() {
    let artifact = experiment(
        "codec-authenticity",
        ExperimentOrigin::Physical {
            apparatus_id: artifact_id("apparatus-1"),
            facility_id: artifact_id("facility-1"),
        },
        true,
        true,
    )
    .expect("valid experiment codec fixture");
    let bytes = artifact.canonical_bytes().expect("canonical experiment");

    for (field, offset) in [
        (
            "experiment.authenticity.source_bytes_hash",
            bytes.len() - 65,
        ),
        (
            "experiment.authenticity.custody_receipt_hash",
            bytes.len() - 33,
        ),
    ] {
        let mut forged = bytes.clone();
        forged[offset..offset + 32].fill(0);
        let error = ExperimentArtifact::from_canonical_bytes(&forged)
            .expect_err("zeroed authenticity hash must refuse canonical decode");
        assert!(
            error.detail().contains(field),
            "expected {field} decode refusal, got {error}",
        );
    }

    let mut foreign_dataset = bytes.clone();
    let dataset_needle = hash("source-bytes");
    let dataset_offset = foreign_dataset
        .windows(32)
        .position(|window| window == dataset_needle.as_bytes())
        .expect("first row source carries the dataset bytes hash");
    foreign_dataset[dataset_offset..dataset_offset + 32]
        .copy_from_slice(hash("foreign-source-bytes").as_bytes());
    let error = ExperimentArtifact::from_canonical_bytes(&foreign_dataset)
        .expect_err("cross-dataset row source must refuse canonical decode");
    assert!(
        error
            .detail()
            .contains("experiment.manifest.dataset_source_bytes_hash"),
        "expected row/dataset closure refusal, got {error}",
    );

    let mut missing_extraction_receipt = bytes;
    let receipt_needle = hash("extraction-receipt-row-cal-1");
    let receipt_offset = missing_extraction_receipt
        .windows(32)
        .position(|window| window == receipt_needle.as_bytes())
        .expect("row source carries the extraction receipt");
    missing_extraction_receipt[receipt_offset..receipt_offset + 32].fill(0);
    let error = ExperimentArtifact::from_canonical_bytes(&missing_extraction_receipt)
        .expect_err("zero extraction receipt must refuse canonical decode");
    assert!(
        error
            .detail()
            .contains("experiment.manifest.extraction_receipt_hash"),
        "expected extraction-receipt refusal, got {error}",
    );
}

#[test]
fn experiment_manifest_row_metadata_is_codec_and_identity_bearing() {
    let artifact = experiment(
        "typed-row-experiment",
        ExperimentOrigin::Physical {
            apparatus_id: artifact_id("apparatus-1"),
            facility_id: artifact_id("facility-1"),
        },
        true,
        true,
    )
    .expect("fully bound experiment");
    let row = artifact
        .manifest()
        .row(&observation_id("val-1"))
        .expect("validation row binding");
    assert_eq!(
        row.source_ref().dataset_source_bytes_hash(),
        hash("source-bytes")
    );
    assert_eq!(
        row.source_ref().locator_domain(),
        "org.frankensim.fs-evidence.test-row-locator.v1"
    );
    assert_eq!(row.source_ref().locator_contract_version(), 1);
    assert_eq!(row.source_ref().locator_hash(), hash("row-val-1"));
    assert_eq!(
        row.source_ref().extraction_receipt_hash(),
        hash("extraction-receipt-row-val-1")
    );
    assert_eq!(row.locator_hash(), hash("row-val-1"));
    assert_eq!(row.qoi(), &qoi_id("length"));
    assert_eq!(row.instrument(), &artifact_id("instrument-1"));
    assert_eq!(row.acquisition_channel(), &artifact_id("channel-1"));
    assert_eq!(row.clock(), &artifact_id("clock-1"));
    assert!(artifact.clocks().contains_clock(row.clock()));
    assert!(
        !artifact
            .clocks()
            .contains_clock(&artifact_id("undeclared-clock"))
    );
    let unsorted_direct_variant = ClockSynchronization::Synchronized {
        clock_ids: vec![artifact_id("clock-z"), artifact_id("clock-a")],
        method: "direct public variant fixture".to_string(),
        max_skew_seconds: 0.0,
        evidence_hash: hash("direct-clock-topology-evidence"),
    };
    assert!(
        unsorted_direct_variant.contains_clock(&artifact_id("clock-a")),
        "public membership must remain correct before constructor canonicalization",
    );

    let bytes = artifact
        .canonical_bytes()
        .expect("canonical experiment bytes");
    let decoded = ExperimentArtifact::from_canonical_bytes(&bytes)
        .expect("typed manifest row survives canonical decode");
    assert_eq!(decoded, artifact);
    assert_eq!(
        decoded
            .canonical_bytes()
            .expect("canonical experiment re-encode"),
        bytes
    );

    let original = ObservationManifest::try_new(vec![(
        observation_id("row-1"),
        manifest_row_with(
            "same-source",
            "length",
            "instrument-1",
            "channel-1",
            "clock-1",
        ),
    )])
    .expect("original manifest");
    let rebound = ObservationManifest::try_new(vec![(
        observation_id("row-1"),
        manifest_row_with(
            "same-source",
            "length",
            "instrument-1",
            "channel-2",
            "clock-1",
        ),
    )])
    .expect("channel-rebound manifest");
    assert_ne!(
        original.canonical_hash(),
        rebound.canonical_hash(),
        "acquisition-channel authority must be identity-bearing"
    );

    let source_manifest_hash = |source| {
        ObservationManifest::try_new(vec![(
            observation_id("row-1"),
            manifest_row_with_source(source, "length", "instrument-1", "channel-1", "clock-1"),
        )])
        .expect("source mutation remains a structurally valid manifest")
        .canonical_hash()
    };
    let baseline_source_hash = source_manifest_hash(observation_source_ref("same-source"));
    for (field, source) in [
        (
            "dataset source bytes",
            observation_source_ref_with(
                "different-source-bytes",
                "org.frankensim.fs-evidence.test-row-locator.v1",
                1,
                "same-source",
                "extraction-receipt-same-source",
            ),
        ),
        (
            "locator domain",
            observation_source_ref_with(
                "source-bytes",
                "org.frankensim.fs-evidence.alternate-row-locator.v1",
                1,
                "same-source",
                "extraction-receipt-same-source",
            ),
        ),
        (
            "locator contract version",
            observation_source_ref_with(
                "source-bytes",
                "org.frankensim.fs-evidence.test-row-locator.v1",
                2,
                "same-source",
                "extraction-receipt-same-source",
            ),
        ),
        (
            "locator hash",
            observation_source_ref_with(
                "source-bytes",
                "org.frankensim.fs-evidence.test-row-locator.v1",
                1,
                "different-source",
                "extraction-receipt-same-source",
            ),
        ),
        (
            "extraction receipt",
            observation_source_ref_with(
                "source-bytes",
                "org.frankensim.fs-evidence.test-row-locator.v1",
                1,
                "same-source",
                "different-extraction-receipt",
            ),
        ),
    ] {
        assert_ne!(
            baseline_source_hash,
            source_manifest_hash(source),
            "typed source field {field} must be manifest-identity-bearing",
        );
    }
}

#[test]
fn observation_manifest_identity_version_and_domain_are_exact() {
    assert_eq!(VV_OBSERVATION_MANIFEST_IDENTITY_VERSION, 3);
    assert_eq!(VV_OBSERVATION_MANIFEST_IDENTITY_VERSION, VV_SCHEMA_VERSION);
    assert_eq!(
        VV_OBSERVATION_MANIFEST_IDENTITY_DOMAIN,
        "org.frankensim.fs-evidence.vv-observation-manifest.v3"
    );
    assert_ne!(
        fs_blake3::hash_domain(VV_OBSERVATION_MANIFEST_IDENTITY_DOMAIN, b"version-guard"),
        fs_blake3::hash_domain(
            "org.frankensim.fs-evidence.vv-observation-manifest.v4",
            b"version-guard",
        ),
        "rotating the identity era must rotate the digest domain",
    );
}

#[test]
fn observation_manifest_identity_preimage_is_exact_and_independently_reproducible() {
    let manifest = ObservationManifest::try_new(vec![
        (
            observation_id("row-a"),
            manifest_row_with_source(
                observation_source_ref_with(
                    "source-bytes",
                    "org.frankensim.fs-evidence.preimage-locator.v1",
                    0x0102_0304,
                    "locator-a",
                    "receipt-a",
                ),
                "length",
                "instrument-a",
                "channel-a",
                "clock-a",
            ),
        ),
        (
            observation_id("row-b"),
            manifest_row_with_source(
                observation_source_ref_with(
                    "source-bytes",
                    "org.frankensim.fs-evidence.preimage-locator.v1",
                    0x0a0b_0c0d,
                    "locator-b",
                    "receipt-b",
                ),
                "temperature",
                "instrument-b",
                "channel-b",
                "clock-b",
            ),
        ),
    ])
    .expect("independent manifest-preimage fixture");

    let mut preimage = Vec::new();
    preimage.extend_from_slice(&(manifest.rows().len() as u64).to_le_bytes());
    for (id, row) in manifest.rows() {
        push_identity_string(&mut preimage, id.as_str());
        let source = row.source_ref();
        preimage.extend_from_slice(source.dataset_source_bytes_hash().as_bytes());
        push_identity_string(&mut preimage, source.locator_domain());
        preimage.extend_from_slice(&source.locator_contract_version().to_le_bytes());
        preimage.extend_from_slice(source.locator_hash().as_bytes());
        preimage.extend_from_slice(source.extraction_receipt_hash().as_bytes());
        for identity in [
            row.qoi().as_str(),
            row.instrument().as_str(),
            row.acquisition_channel().as_str(),
            row.clock().as_str(),
        ] {
            push_identity_string(&mut preimage, identity);
        }
    }

    assert_eq!(
        manifest.canonical_hash(),
        fs_blake3::hash_domain(VV_OBSERVATION_MANIFEST_IDENTITY_DOMAIN, &preimage),
        "manifest identity must equal the independently reconstructed u64/u32 little-endian, length-framed preimage",
    );
}

#[test]
fn observation_manifest_identity_fields_move_independently() {
    let one_row_hash = |id: &str,
                        source: ObservationSourceRef,
                        qoi: &str,
                        instrument: &str,
                        channel: &str,
                        clock: &str| {
        ObservationManifest::try_new(vec![(
            observation_id(id),
            manifest_row_with_source(source, qoi, instrument, channel, clock),
        )])
        .expect("identity mutation remains a valid one-row manifest")
        .canonical_hash()
    };
    let baseline = one_row_hash(
        "row-1",
        observation_source_ref("locator-1"),
        "length",
        "instrument-1",
        "channel-1",
        "clock-1",
    );
    let source_hash = |source| {
        one_row_hash(
            "row-1",
            source,
            "length",
            "instrument-1",
            "channel-1",
            "clock-1",
        )
    };
    for (field, mutated) in [
        (
            "dataset source bytes",
            source_hash(observation_source_ref_with(
                "alternate-source-bytes",
                "org.frankensim.fs-evidence.test-row-locator.v1",
                1,
                "locator-1",
                "extraction-receipt-locator-1",
            )),
        ),
        (
            "locator domain",
            source_hash(observation_source_ref_with(
                "source-bytes",
                "org.frankensim.fs-evidence.alternate-row-locator.v1",
                1,
                "locator-1",
                "extraction-receipt-locator-1",
            )),
        ),
        (
            "locator contract version",
            source_hash(observation_source_ref_with(
                "source-bytes",
                "org.frankensim.fs-evidence.test-row-locator.v1",
                2,
                "locator-1",
                "extraction-receipt-locator-1",
            )),
        ),
        (
            "locator hash",
            source_hash(observation_source_ref_with(
                "source-bytes",
                "org.frankensim.fs-evidence.test-row-locator.v1",
                1,
                "locator-2",
                "extraction-receipt-locator-1",
            )),
        ),
        (
            "extraction receipt",
            source_hash(observation_source_ref_with(
                "source-bytes",
                "org.frankensim.fs-evidence.test-row-locator.v1",
                1,
                "locator-1",
                "alternate-extraction-receipt",
            )),
        ),
        (
            "observation id",
            one_row_hash(
                "row-2",
                observation_source_ref("locator-1"),
                "length",
                "instrument-1",
                "channel-1",
                "clock-1",
            ),
        ),
        (
            "QoI",
            one_row_hash(
                "row-1",
                observation_source_ref("locator-1"),
                "temperature",
                "instrument-1",
                "channel-1",
                "clock-1",
            ),
        ),
        (
            "instrument",
            one_row_hash(
                "row-1",
                observation_source_ref("locator-1"),
                "length",
                "instrument-2",
                "channel-1",
                "clock-1",
            ),
        ),
        (
            "acquisition channel",
            one_row_hash(
                "row-1",
                observation_source_ref("locator-1"),
                "length",
                "instrument-1",
                "channel-2",
                "clock-1",
            ),
        ),
        (
            "clock",
            one_row_hash(
                "row-1",
                observation_source_ref("locator-1"),
                "length",
                "instrument-1",
                "channel-1",
                "clock-2",
            ),
        ),
    ] {
        assert_ne!(baseline, mutated, "{field} must move manifest identity");
    }

    let first = (
        observation_id("row-1"),
        manifest_row_with(
            "locator-1",
            "length",
            "instrument-1",
            "channel-1",
            "clock-1",
        ),
    );
    let second = (
        observation_id("row-2"),
        manifest_row_with(
            "locator-2",
            "length",
            "instrument-1",
            "channel-1",
            "clock-1",
        ),
    );
    let forward = ObservationManifest::try_new(vec![first.clone(), second.clone()])
        .expect("forward input order")
        .canonical_hash();
    let reverse = ObservationManifest::try_new(vec![second, first])
        .expect("reverse input order")
        .canonical_hash();
    assert_eq!(
        forward, reverse,
        "caller order is canonicalized by observation identity"
    );
    assert_ne!(baseline, forward, "row count must move manifest identity");

    let naive_left = format!("{}{}", "length", "sensor");
    let naive_right = format!("{}{}", "lengths", "ensor");
    assert_eq!(
        naive_left, naive_right,
        "fixture must collide without framing"
    );
    let framed_left = one_row_hash(
        "row-1",
        observation_source_ref("locator-1"),
        "length",
        "sensor",
        "channel-1",
        "clock-1",
    );
    let framed_right = one_row_hash(
        "row-1",
        observation_source_ref("locator-1"),
        "lengths",
        "ensor",
        "channel-1",
        "clock-1",
    );
    assert_ne!(
        framed_left, framed_right,
        "length framing must separate concatenation aliases"
    );
}

#[test]
fn vv_artifact_identity_version_and_domain_are_exact() {
    assert_eq!(VV_ARTIFACT_IDENTITY_VERSION, 3);
    assert_eq!(VV_ARTIFACT_IDENTITY_VERSION, VV_SCHEMA_VERSION);
    assert_eq!(
        VV_ARTIFACT_FAMILY,
        "org.frankensim.fs-evidence.vv-artifact.v3"
    );

    let artifact = VvArtifact::from(
        experiment(
            "artifact-identity-version-guard",
            ExperimentOrigin::Physical {
                apparatus_id: artifact_id("apparatus-1"),
                facility_id: artifact_id("facility-1"),
            },
            true,
            true,
        )
        .expect("valid experiment artifact"),
    );
    let bytes = artifact.canonical_bytes().expect("canonical artifact");
    assert_eq!(
        artifact.content_hash().expect("artifact content identity"),
        fs_blake3::hash_domain(VV_ARTIFACT_FAMILY, &bytes),
    );
    assert_ne!(
        fs_blake3::hash_domain(VV_ARTIFACT_FAMILY, &bytes),
        fs_blake3::hash_domain("org.frankensim.fs-evidence.vv-artifact.v4", &bytes),
        "a new artifact identity era must rotate the digest domain",
    );
    for stale_schema in [1_u32, 2, 4] {
        let mut stale = bytes.clone();
        stale[4..8].copy_from_slice(&stale_schema.to_le_bytes());
        let error = VvArtifact::from_canonical_bytes(&stale)
            .expect_err("non-v3 canonical transport must refuse");
        assert_eq!(error.offset(), 4);
        assert!(error.detail().contains("unsupported V&V schema version"));
    }
}

#[test]
fn artifact_header_accuracy_normalizes_signed_zero_before_wire_identity() {
    let base = closed_case(CaseKnobs::default()).context().clone();
    let context_with_accuracy = |accuracy: f64| {
        let base_header = base.header();
        let header = ArtifactHeader::try_new(
            base_header.id().clone(),
            base_header.units().to_vec(),
            base_header.seed().clone(),
            DeclaredBudget::Limit(accuracy),
            base_header.time_ms().clone(),
            base_header.memory_bytes().clone(),
            base_header
                .versions()
                .iter()
                .map(|(component, version)| (component.clone(), version.clone()))
                .collect(),
            base_header.capabilities().iter().cloned().collect(),
        )
        .expect("finite non-negative header accuracy");
        ContextOfUse::try_new(
            header,
            base.decision(),
            base.qois().values().cloned().collect(),
            base.applicability().clone(),
            base.applicability_policy(),
        )
        .expect("context with signed-zero accuracy fixture")
    };

    let positive_zero = context_with_accuracy(0.0);
    let negative_zero = context_with_accuracy(-0.0);
    for context in [&positive_zero, &negative_zero] {
        let DeclaredBudget::Limit(accuracy) = context.header().accuracy() else {
            panic!("fixture must retain an explicit accuracy limit");
        };
        assert_eq!(
            accuracy.to_bits(),
            0.0_f64.to_bits(),
            "accepted header accuracy zero is stored canonically",
        );
    }
    let canonical_bytes = positive_zero
        .canonical_bytes()
        .expect("positive-zero context bytes");
    assert_eq!(
        canonical_bytes,
        negative_zero
            .canonical_bytes()
            .expect("negative-zero context bytes"),
        "mathematically equal header accuracy zeros cannot mint distinct artifact bytes",
    );
    assert_eq!(
        positive_zero
            .content_hash()
            .expect("positive-zero context identity"),
        negative_zero
            .content_hash()
            .expect("negative-zero context identity"),
        "mathematically equal header accuracy zeros cannot mint distinct artifact identities",
    );

    let mut seed_and_budget_prefix = Vec::new();
    seed_and_budget_prefix.push(0);
    seed_and_budget_prefix.extend_from_slice(&0x5eed_u64.to_le_bytes());
    seed_and_budget_prefix.push(0);
    seed_and_budget_prefix.extend_from_slice(&0.0_f64.to_bits().to_le_bytes());
    seed_and_budget_prefix.push(0);
    seed_and_budget_prefix.extend_from_slice(&10_000_u64.to_le_bytes());
    let prefix_offsets = canonical_bytes
        .windows(seed_and_budget_prefix.len())
        .enumerate()
        .filter_map(|(offset, window)| {
            (window == seed_and_budget_prefix.as_slice()).then_some(offset)
        })
        .collect::<Vec<_>>();
    assert_eq!(
        prefix_offsets.len(),
        1,
        "fixture exposes one fixed-seed/accuracy/time header sequence",
    );
    let prefix_offset = prefix_offsets[0];
    let accuracy_offset = prefix_offset + 1 + core::mem::size_of::<u64>() + 1;
    let mut noncanonical_negative_zero = canonical_bytes;
    noncanonical_negative_zero[accuracy_offset..accuracy_offset + core::mem::size_of::<f64>()]
        .copy_from_slice(&(-0.0_f64).to_bits().to_le_bytes());
    let error = ContextOfUse::from_canonical_bytes(&noncanonical_negative_zero)
        .expect_err("decoder must refuse a noncanonical negative-zero header alias");
    assert_eq!(error.rule_name(), "vv-canonical-identity");
    assert!(error.detail().contains("not a canonical fixed point"));
}

#[test]
fn vv_artifact_identity_fields_move_independently() {
    let artifact = VvArtifact::from(
        experiment(
            "artifact-identity-fields",
            ExperimentOrigin::Physical {
                apparatus_id: artifact_id("apparatus-1"),
                facility_id: artifact_id("facility-1"),
            },
            true,
            true,
        )
        .expect("valid experiment artifact"),
    );
    let bytes = artifact.canonical_bytes().expect("canonical artifact");
    assert!(bytes.len() > 14, "artifact transport includes a payload");
    let baseline = artifact.content_hash().expect("artifact content identity");
    let moved = |field: &str, mutated: &[u8]| {
        assert_ne!(
            baseline,
            fs_blake3::hash_domain(VV_ARTIFACT_FAMILY, mutated),
            "{field} must move exact artifact identity",
        );
    };

    let mut magic = bytes.clone();
    magic[0] ^= 1;
    moved("transport magic", &magic);
    let mut schema = bytes.clone();
    schema[4] ^= 1;
    moved("wire schema version", &schema);
    let mut ruleset = bytes.clone();
    ruleset[8] ^= 1;
    moved("ruleset version", &ruleset);
    let mut root = bytes.clone();
    root[12] ^= 1;
    moved("root tag", &root);
    let mut kind = bytes.clone();
    kind[13] ^= 1;
    moved("artifact kind", &kind);
    let mut payload = bytes.clone();
    let last = payload.len() - 1;
    payload[last] ^= 1;
    moved("artifact payload", &payload);

    let mut big_endian_schema = bytes.clone();
    big_endian_schema[4..8].copy_from_slice(&VV_SCHEMA_VERSION.to_be_bytes());
    moved("fixed numeric little endian", &big_endian_schema);

    let dataset_id = b"artifact-identity-fields-dataset";
    let dataset_offset = bytes
        .windows(dataset_id.len())
        .position(|window| window == dataset_id)
        .expect("experiment payload carries its dataset identity");
    let length_offset = dataset_offset
        .checked_sub(8)
        .expect("dataset identity has a u64 length prefix");
    let mut length_framing = bytes.clone();
    length_framing[length_offset] ^= 1;
    moved("length framing", &length_framing);

    let mut field_order = bytes.clone();
    field_order.swap(dataset_offset, dataset_offset + 1);
    moved("canonical field order", &field_order);

    let decoded = VvArtifact::from_canonical_bytes(&bytes).expect("canonical artifact decodes");
    assert_eq!(
        decoded
            .canonical_bytes()
            .expect("canonical artifact re-encodes"),
        bytes,
        "the admitted transport is an exact canonical fixed point",
    );
}

#[test]
fn vv_artifact_valid_semantic_mutations_cover_all_seven_variants_and_concrete_wrappers() {
    let baseline = closed_case(CaseKnobs::default());
    let context_variant = closed_case(CaseKnobs {
        context_domain_lo: 2.0,
        ..CaseKnobs::default()
    });
    let plan_variant = closed_case(CaseKnobs {
        failing_diagnostic: Some(VvRule::DiagnosticConfounding),
        ..CaseKnobs::default()
    });
    let experiment_variant = closed_case(CaseKnobs {
        origin: OriginKind::SyntheticHighFidelity,
        ..CaseKnobs::default()
    });
    let split_variant = closed_case(CaseKnobs {
        blind_source_label: "row-blind-alternate",
        ..CaseKnobs::default()
    });
    let prediction_variant = closed_case(CaseKnobs {
        predicted: 9.96,
        ..CaseKnobs::default()
    });
    let assumptions_variant = closed_case(CaseKnobs {
        tamper_seed: Some("A-001"),
        ..CaseKnobs::default()
    });
    let numerical = |label| {
        NumericalUncertainty::try_new(0.02, hash(label))
            .expect("artifact identity solution variant")
    };
    let solution_variant = SolutionVerificationReceipt::try_new(
        header("solution-variant", "m"),
        artifact_id("solve-variant"),
        qoi_id("length"),
        unit_id("m"),
        numerical("mesh-variant"),
        numerical("time-variant"),
        numerical("nonlinear-variant"),
        numerical("iterative-variant"),
    )
    .expect("valid solution-verification variant");

    let mut seen = Vec::new();
    macro_rules! assert_artifact_pair {
        ($label:literal, $baseline:expr, $variant:expr) => {{
            let concrete_baseline = (*$baseline).clone();
            let concrete_variant = (*$variant).clone();
            let wrapped_baseline = VvArtifact::from(concrete_baseline.clone());
            let wrapped_variant = VvArtifact::from(concrete_variant.clone());
            let baseline_bytes = wrapped_baseline
                .canonical_bytes()
                .expect("baseline wrapper canonical bytes");
            let variant_bytes = wrapped_variant
                .canonical_bytes()
                .expect("variant wrapper canonical bytes");

            assert_eq!(
                wrapped_baseline.kind(),
                wrapped_variant.kind(),
                "{} mutation must stay within one artifact family",
                $label,
            );
            assert_eq!(
                concrete_baseline
                    .canonical_bytes()
                    .expect("baseline concrete canonical bytes"),
                baseline_bytes,
                "{} concrete and enum wrappers must share one exact transport",
                $label,
            );
            assert_eq!(
                concrete_variant
                    .canonical_bytes()
                    .expect("variant concrete canonical bytes"),
                variant_bytes,
                "{} variant concrete and enum wrappers must share one exact transport",
                $label,
            );
            assert_eq!(
                concrete_baseline
                    .content_hash()
                    .expect("baseline concrete identity"),
                wrapped_baseline
                    .content_hash()
                    .expect("baseline wrapper identity"),
                "{} concrete and wrapper identities must agree",
                $label,
            );
            assert_eq!(
                concrete_variant
                    .content_hash()
                    .expect("variant concrete identity"),
                wrapped_variant
                    .content_hash()
                    .expect("variant wrapper identity"),
                "{} variant concrete and wrapper identities must agree",
                $label,
            );
            assert_ne!(
                baseline_bytes, variant_bytes,
                "{} valid semantic mutation must move canonical transport",
                $label,
            );
            assert_ne!(
                wrapped_baseline
                    .content_hash()
                    .expect("baseline wrapper identity"),
                wrapped_variant
                    .content_hash()
                    .expect("variant wrapper identity"),
                "{} valid semantic mutation must move content identity",
                $label,
            );
            assert_eq!(
                VvArtifact::from_canonical_bytes(&variant_bytes)
                    .expect("variant wrapper fixed point"),
                wrapped_variant,
                "{} wrapper transport must round-trip exactly",
                $label,
            );
            seen.push(wrapped_baseline.kind());
        }};
    }

    assert_artifact_pair!("context", baseline.context(), context_variant.context());
    assert_artifact_pair!(
        "validation plan",
        baseline.validation_plan(),
        plan_variant.validation_plan()
    );
    assert_artifact_pair!(
        "experiment",
        baseline.experiments().values().next().expect("experiment"),
        experiment_variant
            .experiments()
            .values()
            .next()
            .expect("experiment variant")
    );
    assert_artifact_pair!(
        "calibration split",
        baseline.splits().values().next().expect("split"),
        split_variant
            .splits()
            .values()
            .next()
            .expect("split variant")
    );
    assert_artifact_pair!(
        "solution verification",
        baseline
            .solution_verification()
            .values()
            .next()
            .expect("solution verification"),
        &solution_variant
    );
    assert_artifact_pair!(
        "prediction assessment",
        baseline.predictions().values().next().expect("prediction"),
        prediction_variant
            .predictions()
            .values()
            .next()
            .expect("prediction variant")
    );
    assert_artifact_pair!(
        "assumptions ledger",
        baseline.assumptions(),
        assumptions_variant.assumptions()
    );

    seen.sort();
    seen.dedup();
    assert_eq!(
        seen,
        vec![
            ArtifactKind::ContextOfUse,
            ArtifactKind::ValidationPlan,
            ArtifactKind::ExperimentArtifact,
            ArtifactKind::CalibrationSplit,
            ArtifactKind::SolutionVerificationReceipt,
            ArtifactKind::PredictionAssessment,
            ArtifactKind::AssumptionsLedger,
        ],
        "the semantic mutation matrix must cover every top-level artifact family",
    );
}

#[test]
fn vv_case_identity_version_domain_and_stage_separation_are_exact() {
    assert_eq!(VV_CASE_IDENTITY_VERSION, 3);
    assert_eq!(VV_CASE_IDENTITY_VERSION, VV_SCHEMA_VERSION);
    assert_eq!(VV_CASE_FAMILY, "org.frankensim.fs-evidence.vv-case.v3");

    let case = closed_case(CaseKnobs::default());
    let bytes = case.canonical_bytes().expect("canonical case");
    let case_hash = case.content_hash().expect("complete-case identity");
    assert_eq!(case_hash, fs_blake3::hash_domain(VV_CASE_FAMILY, &bytes));
    assert_ne!(
        case_hash,
        fs_blake3::hash_domain(VV_ARTIFACT_FAMILY, &bytes),
        "an admission-authoritative complete case must not reuse the individual-artifact domain",
    );
    assert_ne!(
        fs_blake3::hash_domain(VV_CASE_FAMILY, b"version-guard"),
        fs_blake3::hash_domain("org.frankensim.fs-evidence.vv-case.v4", b"version-guard",),
        "a new complete-case identity era must rotate the digest domain",
    );
}

#[test]
fn vv_case_identity_preimage_fields_move_independently() {
    enum CaseFieldMutation {
        Context(ContextOfUse),
        ValidationPlan(ValidationPlan),
        Experiments(Vec<ExperimentArtifact>),
        Splits(Vec<CalibrationSplit>),
        SolutionVerification(Vec<SolutionVerificationReceipt>),
        Predictions(Vec<PredictionAssessment>),
        Assumptions(AssumptionsLedger),
    }

    let baseline = closed_case(CaseKnobs::default());
    let baseline_bytes = baseline.canonical_bytes().expect("canonical case");
    let baseline_hash = baseline.content_hash().expect("complete-case identity");
    assert_eq!(&baseline_bytes[..4], b"FSVV");
    assert_eq!(
        &baseline_bytes[4..8],
        &VV_SCHEMA_VERSION.to_le_bytes(),
        "wire schema version is fixed little-endian state",
    );
    assert_eq!(
        &baseline_bytes[8..12],
        &VV_RULESET_VERSION.to_le_bytes(),
        "ruleset version is fixed little-endian state",
    );
    assert_eq!(
        baseline_bytes[12], 1,
        "complete cases use the case root tag"
    );

    let position = |needle: &[u8]| {
        baseline_bytes
            .windows(needle.len())
            .position(|window| window == needle)
            .expect("fixture identity occurs in the canonical case")
    };
    let preregistration = hash("preregistered-analysis");
    let ordered_artifact_ids = [
        b"Decide whether the measured length satisfies the release criterion.".as_slice(),
        b"validation-plan-1".as_slice(),
        b"experiment-1-dataset".as_slice(),
        preregistration.as_bytes(),
        b"solve-1".as_slice(),
        b"prediction-1".as_slice(),
        b"assumptions".as_slice(),
    ]
    .map(|needle| position(needle));
    assert!(
        ordered_artifact_ids
            .windows(2)
            .all(|pair| pair[0] < pair[1]),
        "the complete-case encoder has one exact top-level field order",
    );

    let context_variant = closed_case(CaseKnobs {
        context_domain_lo: 2.0,
        ..CaseKnobs::default()
    });
    let plan_variant = closed_case(CaseKnobs {
        failing_diagnostic: Some(VvRule::DiagnosticConfounding),
        ..CaseKnobs::default()
    });
    let experiment_variant = closed_case(CaseKnobs {
        origin: OriginKind::SyntheticHighFidelity,
        ..CaseKnobs::default()
    });
    let split_variant = closed_case(CaseKnobs {
        blind_source_label: "row-blind-alternate",
        ..CaseKnobs::default()
    });
    let prediction_variant = closed_case(CaseKnobs {
        predicted: 9.96,
        ..CaseKnobs::default()
    });
    let assumptions_variant = closed_case(CaseKnobs {
        tamper_seed: Some("A-001"),
        ..CaseKnobs::default()
    });
    let numerical = |label| {
        NumericalUncertainty::try_new(0.02, hash(label))
            .expect("complete-case identity solution variant")
    };
    let solution_variant = SolutionVerificationReceipt::try_new(
        header("solution-variant", "m"),
        artifact_id("solve-variant"),
        qoi_id("length"),
        unit_id("m"),
        numerical("mesh-variant"),
        numerical("time-variant"),
        numerical("nonlinear-variant"),
        numerical("iterative-variant"),
    )
    .expect("typed solution-verification variant");

    for (field, mutation) in [
        (
            "context artifact",
            CaseFieldMutation::Context(context_variant.context().clone()),
        ),
        (
            "validation-plan artifact",
            CaseFieldMutation::ValidationPlan(plan_variant.validation_plan().clone()),
        ),
        (
            "experiment artifact registry",
            CaseFieldMutation::Experiments(
                experiment_variant.experiments().values().cloned().collect(),
            ),
        ),
        (
            "calibration-split artifact registry",
            CaseFieldMutation::Splits(split_variant.splits().values().cloned().collect()),
        ),
        (
            "solution-verification artifact registry",
            CaseFieldMutation::SolutionVerification(vec![solution_variant.clone()]),
        ),
        (
            "prediction-assessment artifact registry",
            CaseFieldMutation::Predictions(
                prediction_variant.predictions().values().cloned().collect(),
            ),
        ),
        (
            "assumptions-ledger artifact",
            CaseFieldMutation::Assumptions(assumptions_variant.assumptions().clone()),
        ),
    ] {
        let mut context = baseline.context().clone();
        let mut validation_plan = baseline.validation_plan().clone();
        let mut experiments = baseline.experiments().values().cloned().collect::<Vec<_>>();
        let mut splits = baseline.splits().values().cloned().collect::<Vec<_>>();
        let mut solution_verification = baseline
            .solution_verification()
            .values()
            .cloned()
            .collect::<Vec<_>>();
        let mut predictions = baseline.predictions().values().cloned().collect::<Vec<_>>();
        let mut assumptions = baseline.assumptions().clone();
        match mutation {
            CaseFieldMutation::Context(value) => context = value,
            CaseFieldMutation::ValidationPlan(value) => validation_plan = value,
            CaseFieldMutation::Experiments(value) => experiments = value,
            CaseFieldMutation::Splits(value) => splits = value,
            CaseFieldMutation::SolutionVerification(value) => solution_verification = value,
            CaseFieldMutation::Predictions(value) => predictions = value,
            CaseFieldMutation::Assumptions(value) => assumptions = value,
        }
        let variant = VvCase::try_new(
            context,
            validation_plan,
            experiments,
            splits,
            solution_verification,
            predictions,
            assumptions,
        )
        .expect("one-field complete-case identity variant");
        assert_ne!(
            baseline_hash,
            variant
                .content_hash()
                .expect("variant complete-case identity"),
            "{field} must move complete-case identity",
        );
    }

    let moved_preimage = |field: &str, mutated: &[u8]| {
        assert_ne!(
            baseline_hash,
            fs_blake3::hash_domain(VV_CASE_FAMILY, mutated),
            "{field} must move the exact complete-case preimage",
        );
    };
    let mut magic = baseline_bytes.clone();
    magic[0] ^= 1;
    moved_preimage("transport magic", &magic);
    let mut schema = baseline_bytes.clone();
    schema[4] ^= 1;
    moved_preimage("wire schema version", &schema);
    let mut ruleset = baseline_bytes.clone();
    ruleset[8] ^= 1;
    moved_preimage("ruleset version", &ruleset);
    let mut root = baseline_bytes.clone();
    root[12] ^= 1;
    moved_preimage("root tag", &root);
    let context_id = b"context-1";
    let context_id_at = position(context_id);
    let context_length_at = context_id_at
        .checked_sub(8)
        .expect("context id has a u64 length prefix");
    let mut framing = baseline_bytes.clone();
    framing[context_length_at] ^= 1;
    moved_preimage("length framing", &framing);

    let decoded =
        VvCase::from_canonical_bytes(&baseline_bytes).expect("canonical complete case decodes");
    assert_eq!(decoded, baseline);
    assert_eq!(
        decoded.canonical_bytes().expect("complete-case re-encode"),
        baseline_bytes,
        "complete-case identity transport is an exact fixed point",
    );
    assert_eq!(
        decoded
            .content_hash()
            .expect("decoded complete-case identity"),
        baseline_hash,
    );
}

#[test]
fn experiment_manifest_cross_wiring_refuses_at_admission() {
    let calibration = || {
        vec![InstrumentCalibration::new(
            artifact_id("instrument-1"),
            hash("instrument-calibration"),
            true,
        )]
    };
    let single_clock = || ClockSynchronization::SingleClock {
        clock_id: artifact_id("clock-1"),
    };

    let wrong_qoi = ObservationManifest::try_new(vec![(
        observation_id("row-1"),
        manifest_row_with(
            "row-1",
            "temperature",
            "instrument-1",
            "channel-1",
            "clock-1",
        ),
    )])
    .expect("structurally typed but semantically cross-wired QoI");
    assert_only_rule_field(
        experiment_with_manifest(
            wrong_qoi,
            vec![qoi_id("length")],
            calibration(),
            single_clock(),
        ),
        VvRule::QoiDependencyClosed,
        "experiment.manifest.qoi",
    );

    let wrong_instrument = ObservationManifest::try_new(vec![(
        observation_id("row-1"),
        manifest_row_with("row-1", "length", "instrument-2", "channel-1", "clock-1"),
    )])
    .expect("structurally typed but uncalibrated instrument");
    assert_only_rule_field(
        experiment_with_manifest(
            wrong_instrument,
            vec![qoi_id("length")],
            calibration(),
            single_clock(),
        ),
        VvRule::ExperimentInstrumentCalibration,
        "experiment.manifest.instrument",
    );

    let zero_certificate =
        ObservationManifest::try_new(vec![(observation_id("row-1"), manifest_row("row-1"))])
            .expect("row bound to the declared instrument");
    assert_only_rule_field(
        experiment_with_manifest(
            zero_certificate,
            vec![qoi_id("length")],
            vec![InstrumentCalibration::new(
                artifact_id("instrument-1"),
                ContentHash([0; 32]),
                true,
            )],
            single_clock(),
        ),
        VvRule::ExperimentInstrumentCalibration,
        "experiment.instruments",
    );

    let wrong_clock = ObservationManifest::try_new(vec![(
        observation_id("row-1"),
        manifest_row_with("row-1", "length", "instrument-1", "channel-1", "clock-2"),
    )])
    .expect("structurally typed but unsynchronized clock");
    assert_only_rule_field(
        experiment_with_manifest(
            wrong_clock,
            vec![qoi_id("length")],
            calibration(),
            single_clock(),
        ),
        VvRule::ExperimentClockSynchronization,
        "experiment.manifest.clock",
    );
}

#[test]
fn prediction_context_and_plan_stale_reference_matrix_refuses_exact_fields() {
    let baseline = closed_case(CaseKnobs::default());
    let prediction = baseline
        .predictions()
        .values()
        .next()
        .expect("baseline prediction");
    let stale_context = ArtifactRef::new(
        prediction.context().kind(),
        prediction.context().id().clone(),
        hash("stale-context-content"),
    );
    let stale_plan = ArtifactRef::new(
        prediction.validation_plan().kind(),
        prediction.validation_plan().id().clone(),
        hash("stale-validation-plan-content"),
    );
    assert_ne!(stale_context.hash(), prediction.context().hash());
    assert_ne!(stale_plan.hash(), prediction.validation_plan().hash());

    for (field, context, plan) in [
        (
            "prediction.context",
            stale_context,
            prediction.validation_plan().clone(),
        ),
        (
            "prediction.validation_plan",
            prediction.context().clone(),
            stale_plan,
        ),
    ] {
        let rebuilt = rebuild_prediction(
            prediction,
            context,
            plan,
            prediction.dependencies().to_vec(),
        );
        assert_only_rule_field(
            case_replacing_prediction(&baseline, rebuilt).validate(),
            VvRule::QoiDependencyClosed,
            field,
        );
    }
}

#[test]
fn validation_plan_stale_context_hash_refuses_before_prediction_consumption() {
    let baseline = closed_case(CaseKnobs::default());
    let stale_context = ArtifactRef::new(
        ArtifactKind::ContextOfUse,
        baseline.context().id().clone(),
        hash("stale-plan-context-content"),
    );
    let stale_plan = ValidationPlan::try_new(
        baseline.validation_plan().header().clone(),
        stale_context,
        baseline
            .validation_plan()
            .by_qoi()
            .values()
            .cloned()
            .collect(),
    )
    .expect("same-kind stale context remains structurally typed");
    let stale_plan_ref = artifact_reference(&stale_plan);
    let prediction = baseline
        .predictions()
        .values()
        .next()
        .expect("baseline prediction");
    let rebuilt_prediction = rebuild_prediction(
        prediction,
        prediction.context().clone(),
        stale_plan_ref,
        prediction.dependencies().to_vec(),
    );
    let stale_case = VvCase::try_new(
        baseline.context().clone(),
        stale_plan,
        baseline.experiments().values().cloned().collect(),
        baseline.splits().values().cloned().collect(),
        baseline.solution_verification().values().cloned().collect(),
        vec![rebuilt_prediction],
        baseline.assumptions().clone(),
    )
    .expect("case remains structurally closed around the stale plan");

    assert_only_rule_field(
        stale_case.validate(),
        VvRule::QoiDependencyClosed,
        "validation_plan.context",
    );
}

#[test]
fn validation_plan_stale_split_hash_refuses_at_the_plan_edge() {
    let baseline = closed_case(CaseKnobs::default());
    let plan_row = baseline
        .validation_plan()
        .by_qoi()
        .values()
        .next()
        .expect("baseline validation-plan row");
    let stale_split = ArtifactRef::new(
        ArtifactKind::CalibrationSplit,
        plan_row.split().id().clone(),
        hash("stale-plan-split-content"),
    );
    let rebuilt_row = QoiValidationPlan::try_new(
        plan_row.qoi().clone(),
        plan_row.experiments().to_vec(),
        stale_split,
        plan_row.metrics().to_vec(),
        plan_row.diagnostics().clone(),
    )
    .expect("same-kind stale split remains structurally typed");
    let rebuilt_plan = ValidationPlan::try_new(
        baseline.validation_plan().header().clone(),
        baseline.validation_plan().context().clone(),
        vec![rebuilt_row],
    )
    .expect("stale split remains an unresolved plan reference");
    let prediction = baseline
        .predictions()
        .values()
        .next()
        .expect("baseline prediction");
    let rebuilt_prediction = rebuild_prediction(
        prediction,
        prediction.context().clone(),
        artifact_reference(&rebuilt_plan),
        prediction.dependencies().to_vec(),
    );
    let stale_case = VvCase::try_new(
        baseline.context().clone(),
        rebuilt_plan,
        baseline.experiments().values().cloned().collect(),
        baseline.splits().values().cloned().collect(),
        baseline.solution_verification().values().cloned().collect(),
        vec![rebuilt_prediction],
        baseline.assumptions().clone(),
    )
    .expect("case retains the unresolved same-id plan edge");

    let error = stale_case
        .validate()
        .expect_err("stale validation-plan split content identity must refuse");
    let plan_edge_violations = error
        .violations()
        .iter()
        .filter(|violation| {
            violation.rule() == VvRule::SplitPartitionsDisjoint
                && violation.field() == "validation_plan.split"
        })
        .count();
    assert_eq!(
        plan_edge_violations, 1,
        "the stale split must be diagnosed exactly once at validation_plan.split: {error}",
    );
    let split_fields = error
        .violations()
        .iter()
        .filter(|violation| violation.rule() == VvRule::SplitPartitionsDisjoint)
        .map(VvViolation::field)
        .collect::<std::collections::BTreeSet<_>>();
    assert_eq!(
        split_fields,
        std::collections::BTreeSet::from([
            "prediction.physical_validation.observations",
            "prediction.posterior_checks.observations",
            "prediction.validation_metrics.observations",
            "validation_plan.split",
        ]),
        "plan and all three prediction consumers must diagnose the same stale split authority: {error}",
    );
}

#[test]
fn calibration_split_stale_experiment_hash_refuses_at_the_split_edge() {
    let baseline = closed_case(CaseKnobs::default());
    let split = baseline.splits().values().next().expect("baseline split");
    let stale_experiment = ArtifactRef::new(
        ArtifactKind::ExperimentArtifact,
        split.experiment().id().clone(),
        hash("stale-split-experiment-content"),
    );
    let rebuilt_split = CalibrationSplit::try_new(
        split.header().clone(),
        stale_experiment,
        split.preregistration_hash(),
        split.calibration_ids().iter().cloned().collect(),
        split.validation_ids().iter().cloned().collect(),
        split
            .blind_sources()
            .iter()
            .map(|(id, source)| (id.clone(), *source))
            .collect(),
    )
    .expect("same-kind stale experiment remains structurally typed");
    let rebuilt_split_ref = artifact_reference(&rebuilt_split);
    let plan_row = baseline
        .validation_plan()
        .by_qoi()
        .values()
        .next()
        .expect("baseline validation-plan row");
    let rebuilt_row = QoiValidationPlan::try_new(
        plan_row.qoi().clone(),
        plan_row.experiments().to_vec(),
        rebuilt_split_ref,
        plan_row.metrics().to_vec(),
        plan_row.diagnostics().clone(),
    )
    .expect("plan follows the rebuilt split artifact");
    let rebuilt_plan = ValidationPlan::try_new(
        baseline.validation_plan().header().clone(),
        baseline.validation_plan().context().clone(),
        vec![rebuilt_row],
    )
    .expect("rebuilt plan");
    let prediction = baseline
        .predictions()
        .values()
        .next()
        .expect("baseline prediction");
    let rebuilt_prediction = rebuild_prediction(
        prediction,
        prediction.context().clone(),
        artifact_reference(&rebuilt_plan),
        prediction.dependencies().to_vec(),
    );
    let stale_case = VvCase::try_new(
        baseline.context().clone(),
        rebuilt_plan,
        baseline.experiments().values().cloned().collect(),
        vec![rebuilt_split],
        baseline.solution_verification().values().cloned().collect(),
        vec![rebuilt_prediction],
        baseline.assumptions().clone(),
    )
    .expect("case retains the unresolved same-id split edge");
    let error = stale_case
        .validate()
        .expect_err("stale split-to-experiment content identity must refuse");
    let split_edge_violations = error
        .violations()
        .iter()
        .filter(|violation| {
            violation.rule() == VvRule::SplitPartitionsDisjoint
                && violation.field() == "split.experiment"
        })
        .count();
    assert_eq!(
        split_edge_violations, 1,
        "the stale experiment must be diagnosed exactly once at split.experiment: {error}",
    );
}

#[test]
fn physical_dependency_cross_wired_to_another_planned_experiment_refuses_exactly() {
    let baseline = closed_case(CaseKnobs::default());
    let first_experiment = baseline
        .experiments()
        .values()
        .next()
        .expect("first physical experiment")
        .clone();
    let second_experiment = experiment(
        "experiment-2",
        ExperimentOrigin::Physical {
            apparatus_id: artifact_id("apparatus-2"),
            facility_id: artifact_id("facility-2"),
        },
        true,
        true,
    )
    .expect("second physical experiment");
    let first_experiment_ref = artifact_reference(&first_experiment);
    let second_experiment_ref = artifact_reference(&second_experiment);
    let baseline_plan_row = baseline
        .validation_plan()
        .by_qoi()
        .get(&qoi_id("length"))
        .expect("baseline plan row");
    let expanded_plan_row = QoiValidationPlan::try_new(
        baseline_plan_row.qoi().clone(),
        vec![first_experiment_ref, second_experiment_ref.clone()],
        baseline_plan_row.split().clone(),
        baseline_plan_row.metrics().to_vec(),
        baseline_plan_row.diagnostics().clone(),
    )
    .expect("both physical experiments are declared by the plan");
    let expanded_plan = ValidationPlan::try_new(
        baseline.validation_plan().header().clone(),
        baseline.validation_plan().context().clone(),
        vec![expanded_plan_row],
    )
    .expect("expanded validation plan");
    let expanded_plan_ref = artifact_reference(&expanded_plan);

    let baseline_prediction = baseline
        .predictions()
        .values()
        .next()
        .expect("baseline prediction");
    let dependencies = baseline_prediction
        .dependencies()
        .iter()
        .map(|dependency| {
            if dependency.role() == DependencyRole::PhysicalValidation {
                EvidenceDependency::physical_validation(
                    dependency.qoi().clone(),
                    second_experiment_ref.clone(),
                    dependency
                        .observations()
                        .expect("physical dependency selection")
                        .clone(),
                )
            } else {
                dependency.clone()
            }
        })
        .collect();
    let cross_wired_prediction = rebuild_prediction(
        baseline_prediction,
        baseline_prediction.context().clone(),
        expanded_plan_ref,
        dependencies,
    );
    let cross_wired_case = VvCase::try_new(
        baseline.context().clone(),
        expanded_plan,
        vec![first_experiment, second_experiment],
        baseline.splits().values().cloned().collect(),
        baseline.solution_verification().values().cloned().collect(),
        vec![cross_wired_prediction],
        baseline.assumptions().clone(),
    )
    .expect("structurally closed cross-wired case");

    assert_only_rule_field(
        cross_wired_case.validate(),
        VvRule::ValidationRequiresPhysicalReferent,
        "prediction.physical_validation.observations",
    );
}

#[test]
fn clocks_repeatability_and_covariance_fail_closed() {
    assert_rule(
        ClockSynchronization::synchronized(
            vec![artifact_id("clock-1"), artifact_id("clock-1")],
            "PTP",
            1.0e-6,
            hash("clock-sync"),
        ),
        VvRule::ExperimentClockSynchronization,
    );
    assert_rule(
        CovarianceMatrix::try_new(2, vec![1.0, 2.0, 1.0]),
        VvRule::ExperimentRepeatabilityCovariance,
    );
    assert_rule(
        RepeatabilitySummary::try_new(1, covariance_1()),
        VvRule::ExperimentRepeatabilityCovariance,
    );
}

#[test]
fn repeatability_covariance_binds_explicit_qoi_axis_order_on_wire_and_identity() {
    let covariance = |lower_triangle| {
        CovarianceMatrix::try_new(2, lower_triangle).expect("positive-definite two-QoI covariance")
    };
    let manifest = || {
        ObservationManifest::try_new(vec![
            (
                observation_id("length-row"),
                manifest_row_with(
                    "length-locator",
                    "length",
                    "instrument-1",
                    "channel-1",
                    "clock-1",
                ),
            ),
            (
                observation_id("temperature-row"),
                manifest_row_with(
                    "temperature-locator",
                    "temperature",
                    "instrument-1",
                    "channel-1",
                    "clock-1",
                ),
            ),
        ])
        .expect("two-QoI observation manifest")
    };
    let experiment_for = |qoi_order: Vec<QoiId>, lower_triangle: Vec<f64>| {
        ExperimentArtifact::try_new(
            header_with_units("covariance-axis-experiment", &["m", "K"]),
            artifact_id("covariance-axis-dataset"),
            ExperimentOrigin::Physical {
                apparatus_id: artifact_id("apparatus-1"),
                facility_id: artifact_id("facility-1"),
            },
            qoi_order,
            manifest(),
            vec![InstrumentCalibration::new(
                artifact_id("instrument-1"),
                hash("instrument-calibration"),
                true,
            )],
            ClockSynchronization::SingleClock {
                clock_id: artifact_id("clock-1"),
            },
            RepeatabilitySummary::try_new(3, covariance(lower_triangle))
                .expect("unbound repeatability"),
            DataAuthenticity::new(hash("source-bytes"), hash("custody"), true),
        )
        .expect("axis-bound experiment")
    };

    let length_first = experiment_for(
        vec![qoi_id("length"), qoi_id("temperature")],
        vec![1.0, 0.25, 4.0],
    );
    let temperature_first_equivalent = experiment_for(
        vec![qoi_id("temperature"), qoi_id("length")],
        vec![4.0, 0.25, 1.0],
    );
    let temperature_first_relabelled = experiment_for(
        vec![qoi_id("temperature"), qoi_id("length")],
        vec![1.0, 0.25, 4.0],
    );
    assert_eq!(
        length_first.repeatability().qoi_order(),
        &[qoi_id("length"), qoi_id("temperature")],
    );
    assert_eq!(
        temperature_first_equivalent.repeatability().qoi_order(),
        &[qoi_id("length"), qoi_id("temperature")],
        "paired tensor permutations normalize to sorted QoI axes",
    );
    assert_eq!(
        length_first.qois(),
        temperature_first_equivalent.qois(),
        "the experiment QoI set is intentionally unchanged",
    );
    assert_eq!(
        length_first
            .canonical_bytes()
            .expect("length-first transport"),
        temperature_first_equivalent
            .canonical_bytes()
            .expect("equivalent temperature-first transport"),
        "simultaneous axis and tensor permutation is representation-only",
    );
    assert_eq!(
        length_first.content_hash().expect("length-first identity"),
        temperature_first_equivalent
            .content_hash()
            .expect("equivalent temperature-first identity"),
        "equivalent labeled tensors must share canonical identity",
    );
    assert_ne!(
        length_first
            .canonical_bytes()
            .expect("length-first transport"),
        temperature_first_relabelled
            .canonical_bytes()
            .expect("relabelled temperature-first transport"),
        "relabeling axes without permuting the matrix changes tensor semantics",
    );
    assert_ne!(
        length_first.content_hash().expect("length-first identity"),
        temperature_first_relabelled
            .content_hash()
            .expect("relabelled temperature-first identity"),
        "a semantically different labeled tensor must move identity",
    );
    let bytes = temperature_first_relabelled
        .canonical_bytes()
        .expect("temperature-first transport");
    let decoded = ExperimentArtifact::from_canonical_bytes(&bytes)
        .expect("explicit covariance order decodes");
    assert_eq!(decoded, temperature_first_relabelled);
    assert_eq!(
        decoded.repeatability().qoi_order(),
        &[qoi_id("length"), qoi_id("temperature")],
        "decode retains the canonical labeled-tensor order",
    );

    assert_only_rule_field(
        RepeatabilitySummary::try_new_for_qois(
            3,
            vec![qoi_id("length"), qoi_id("length")],
            covariance(vec![1.0, 0.25, 4.0]),
        ),
        VvRule::ExperimentRepeatabilityCovariance,
        "experiment.repeatability.qoi_order",
    );
    assert_only_rule_field(
        RepeatabilitySummary::try_new_for_qois(
            3,
            vec![qoi_id("length")],
            covariance(vec![1.0, 0.25, 4.0]),
        ),
        VvRule::ExperimentRepeatabilityCovariance,
        "experiment.repeatability.qoi_order",
    );
}

#[test]
fn covariance_signed_zero_has_one_canonical_representation() {
    let positive_zero = CovarianceMatrix::try_new(2, vec![1.0, 0.0, 4.0]).expect("PSD covariance");
    let negative_zero = CovarianceMatrix::try_new(2, vec![1.0, -0.0, 4.0]).expect("PSD covariance");

    assert_eq!(positive_zero, negative_zero);
    assert_eq!(
        positive_zero
            .lower_triangle()
            .iter()
            .map(|value| value.to_bits())
            .collect::<Vec<_>>(),
        negative_zero
            .lower_triangle()
            .iter()
            .map(|value| value.to_bits())
            .collect::<Vec<_>>(),
        "scientifically identical signed zeros must not mint distinct artifact bytes",
    );
    assert_eq!(
        negative_zero.lower_triangle()[1].to_bits(),
        0.0f64.to_bits()
    );

    let experiment_with_zero = |zero: f64| {
        ExperimentArtifact::try_new(
            header_with_units("covariance-zero-experiment", &["m", "K"]),
            artifact_id("covariance-zero-dataset"),
            ExperimentOrigin::Physical {
                apparatus_id: artifact_id("apparatus-1"),
                facility_id: artifact_id("facility-1"),
            },
            vec![qoi_id("length"), qoi_id("temperature")],
            ObservationManifest::try_new(vec![
                (
                    observation_id("length-row"),
                    manifest_row_with(
                        "length-locator",
                        "length",
                        "instrument-1",
                        "channel-1",
                        "clock-1",
                    ),
                ),
                (
                    observation_id("temperature-row"),
                    manifest_row_with(
                        "temperature-locator",
                        "temperature",
                        "instrument-1",
                        "channel-1",
                        "clock-1",
                    ),
                ),
            ])
            .expect("two-QoI observation manifest"),
            vec![InstrumentCalibration::new(
                artifact_id("instrument-1"),
                hash("instrument-calibration"),
                true,
            )],
            ClockSynchronization::SingleClock {
                clock_id: artifact_id("clock-1"),
            },
            RepeatabilitySummary::try_new(
                3,
                CovarianceMatrix::try_new(2, vec![1.0, zero, 4.0])
                    .expect("PSD covariance with a mathematical zero"),
            )
            .expect("repeatability summary"),
            DataAuthenticity::new(hash("source-bytes"), hash("custody"), true),
        )
        .expect("axis-bound experiment")
    };
    let positive_experiment = experiment_with_zero(0.0);
    let negative_experiment = experiment_with_zero(-0.0);
    let canonical = positive_experiment
        .canonical_bytes()
        .expect("positive-zero experiment transport");
    assert_eq!(
        canonical,
        negative_experiment
            .canonical_bytes()
            .expect("negative-zero experiment transport"),
        "the complete artifact transport must normalize covariance signed zero",
    );
    assert_eq!(
        positive_experiment
            .content_hash()
            .expect("positive-zero experiment identity"),
        negative_experiment
            .content_hash()
            .expect("negative-zero experiment identity"),
        "the complete artifact identity must normalize covariance signed zero",
    );

    let mut covariance_needle = Vec::new();
    covariance_needle.extend_from_slice(&2_u64.to_le_bytes());
    covariance_needle.extend_from_slice(&3_u64.to_le_bytes());
    covariance_needle.extend_from_slice(&1.0_f64.to_bits().to_le_bytes());
    covariance_needle.extend_from_slice(&0.0_f64.to_bits().to_le_bytes());
    covariance_needle.extend_from_slice(&4.0_f64.to_bits().to_le_bytes());
    let offsets = canonical
        .windows(covariance_needle.len())
        .enumerate()
        .filter_map(|(offset, window)| (window == covariance_needle.as_slice()).then_some(offset))
        .collect::<Vec<_>>();
    assert_eq!(offsets.len(), 1, "fixture has one exact covariance payload");
    let mut forged = canonical;
    let covariance_zero_offset = offsets[0] + 3 * core::mem::size_of::<u64>();
    forged[covariance_zero_offset..covariance_zero_offset + core::mem::size_of::<u64>()]
        .copy_from_slice(&(-0.0_f64).to_bits().to_le_bytes());
    let error = ExperimentArtifact::from_canonical_bytes(&forged)
        .expect_err("a forged negative-zero wire alias must not be a canonical fixed point");
    assert_eq!(error.offset(), 0);
    assert!(error.detail().contains("not a canonical fixed point"));
}

#[test]
fn covariance_canonical_axis_transport_refuses_an_order_sensitive_ldlt_fixed_point() {
    // This nearly singular PSD tensor passes the deliberately fail-closed
    // floating LDL^T predicate in declared [a, c, b] order but not after the
    // exact canonical [a, b, c] permutation. Admission must refuse it instead
    // of emitting canonical bytes that the decoder would reject.
    let declared_order_covariance = CovarianceMatrix::try_new(
        3,
        vec![
            2_310_095_392.549_797,
            -297_069.917_183_430_46,
            38.290_817_715_431_54,
            595_465.321_408_436,
            -76.887_196_473_647_75,
            154.592_072_504_607_76,
        ],
    )
    .expect("the witness is accepted in its declared axis order");

    assert_only_rule_field(
        RepeatabilitySummary::try_new_for_qois(
            3,
            vec![qoi_id("a"), qoi_id("c"), qoi_id("b")],
            declared_order_covariance,
        ),
        VvRule::ExperimentRepeatabilityCovariance,
        "experiment.covariance",
    );
}

#[test]
fn diagnostic_failures_round_trip_as_evidence_instead_of_disappearing() {
    let passed =
        DiagnosticRecord::try_new(true, hash("diagnostic-pass"), "independent check passed")
            .expect("positive diagnostic");
    let failed = DiagnosticRecord::try_new(
        false,
        hash("diagnostic-fail"),
        "inverse-crime independence was not established",
    )
    .expect("adverse diagnostics remain representable");
    let plan = DiagnosticPlan::new(passed.clone(), passed.clone(), passed, failed);
    assert!(plan.observability().passed());
    assert!(!plan.inverse_crime().passed());
    assert_eq!(
        plan.inverse_crime().detail(),
        "inverse-crime independence was not established"
    );
}

#[test]
fn validation_metrics_include_experimental_and_numerical_uncertainty() {
    let split = split_with(vec!["cal-1"], vec!["val-1"], vec!["blind-1"]).expect("valid split");
    let split_ref = artifact_reference(&split);
    let selection = split
        .validation_selection(split_ref, vec![observation_id("val-1")])
        .expect("held-out validation selection");
    let metric = ValidationMetric::try_new(
        artifact_id("normalized-discrepancy"),
        qoi_id("length"),
        selection,
        9.8,
        10.0,
        0.3,
        0.4,
    )
    .expect("metric with both uncertainty sources");
    assert_eq!(metric.experimental_uncertainty(), 0.3);
    assert_eq!(metric.numerical_uncertainty(), 0.4);
    assert!(metric.combined_uncertainty() >= 0.7);
    assert!(metric.combined_uncertainty() - 0.7 < 1.0e-12);
}

#[test]
fn solution_verification_receipt_accounts_for_all_four_numerical_sources() {
    let component =
        |name| NumericalUncertainty::try_new(0.1, hash(name)).expect("numerical component");
    let receipt = SolutionVerificationReceipt::try_new(
        header("solution-verification", "m"),
        artifact_id("solve-1"),
        qoi_id("length"),
        unit_id("m"),
        component("mesh"),
        component("time"),
        component("nonlinear"),
        component("iterative"),
    )
    .expect("complete solution-verification receipt");
    assert!(receipt.combined_half_width() >= 0.4);
    assert!(receipt.combined_half_width() - 0.4 < 1.0e-12);
}

#[test]
fn categorical_evidence_axes_are_complete_categories_not_scores() {
    assert_rule(
        EvidenceAxes::try_new(vec![(
            EvidenceAxis::CodeVerification,
            EvidenceAxisStatus::Missing {
                reason: "not run".to_owned(),
            },
        )]),
        VvRule::ColorCategoricalOnly,
    );

    let axes = categorical_axes();
    assert_eq!(axes.axes().len(), EvidenceAxis::ALL.len());
    assert!(
        axes.axes()
            .values()
            .all(|status| matches!(status, EvidenceAxisStatus::Missing { .. }))
    );
}

#[test]
fn bound_and_probabilistic_waterfalls_recompute_declared_arithmetic() {
    let bound = UncertaintyWaterfall::try_new(
        qoi_id("length"),
        unit_id("m"),
        WaterfallMode::GuaranteedBound,
        uncertainty_terms([1.0, 2.0, 3.0, 4.0, 5.0, 6.0]),
    )
    .expect("complete bound waterfall");
    assert!(bound.total() >= 21.0);
    assert!(bound.total() - 21.0 < 1.0e-12);

    let mut identity = vec![0.0; 36];
    for index in 0..6 {
        identity[index * 6 + index] = 1.0;
    }
    let independent = UncertaintyWaterfall::try_new(
        qoi_id("length"),
        unit_id("m"),
        WaterfallMode::Probabilistic {
            confidence: 0.95,
            dependence: CorrelationMatrix::try_new(6, identity).expect("identity dependence"),
        },
        uncertainty_terms([3.0, 4.0, 0.0, 0.0, 0.0, 0.0]),
    )
    .expect("independent probabilistic waterfall");
    assert!(independent.total() >= 5.0);
    assert!(independent.total() - 5.0 < 1.0e-12);

    let mut correlated = vec![0.0; 36];
    for index in 0..6 {
        correlated[index * 6 + index] = 1.0;
    }
    correlated[1] = 1.0;
    correlated[6] = 1.0;
    let correlated = UncertaintyWaterfall::try_new(
        qoi_id("length"),
        unit_id("m"),
        WaterfallMode::Probabilistic {
            confidence: 0.95,
            dependence: CorrelationMatrix::try_new(6, correlated)
                .expect("positive-semidefinite correlated block"),
        },
        uncertainty_terms([3.0, 4.0, 0.0, 0.0, 0.0, 0.0]),
    )
    .expect("correlated probabilistic waterfall");
    assert!(correlated.total() >= 7.0);
    assert!(correlated.total() - 7.0 < 1.0e-12);
}

#[test]
fn waterfall_categories_and_dependence_are_mandatory() {
    let mut incomplete = uncertainty_terms([1.0; 6]);
    incomplete.pop();
    assert_rule(
        UncertaintyWaterfall::try_new(
            qoi_id("length"),
            unit_id("m"),
            WaterfallMode::GuaranteedBound,
            incomplete,
        ),
        VvRule::WaterfallModeDeclared,
    );

    assert_rule(
        CorrelationMatrix::try_new(2, vec![1.0, 0.9, 0.0, 1.0]),
        VvRule::WaterfallDependenceDeclared,
    );
    assert_rule(
        UncertaintyTerm::try_new(
            PredictionUncertaintyKind::Data,
            f64::NAN,
            external_target("bad-data"),
        ),
        VvRule::WaterfallArithmetic,
    );
}

#[test]
fn prediction_constructor_refuses_cross_qoi_evidence() {
    let qoi = qoi_id("length");
    let other = qoi_id("temperature");
    let dependency = EvidenceDependency::new(
        other,
        DependencyRole::CodeVerification,
        external_target("code-check"),
    );
    let waterfall = UncertaintyWaterfall::try_new(
        qoi.clone(),
        unit_id("m"),
        WaterfallMode::GuaranteedBound,
        uncertainty_terms([0.0; 6]),
    )
    .expect("waterfall fixture");

    assert_rule(
        PredictionAssessment::try_new(
            header("prediction-1", "m"),
            reference(ArtifactKind::ContextOfUse, "context-1"),
            reference(ArtifactKind::ValidationPlan, "validation-plan-1"),
            qoi,
            vec![dependency],
            waterfall,
            Vec::new(),
            Vec::new(),
            ApplicabilityPoint::try_new(Vec::new(), Vec::new()).expect("empty applicability point"),
            ApplicabilityDecision::InDomain,
            categorical_axes(),
            Vec::new(),
        ),
        VvRule::QoiDependencyIsolated,
    );
}

#[test]
fn applicability_inputs_refuse_reversed_and_nonfinite_values() {
    assert_rule(
        NumericDomainAxis::try_new(axis_id("Re"), unit_id("unitless"), 10.0, 1.0),
        VvRule::AssumptionDomainEnforced,
    );
    assert_rule(
        ApplicabilityPoint::try_new(vec![(axis_id("Re"), f64::NAN)], Vec::new()),
        VvRule::ApplicabilityDecision,
    );
}

#[test]
fn program_assumptions_seed_is_exact_and_operational() {
    let ledger = AssumptionsLedger::try_program_seed(header("assumptions", "unitless"))
        .expect("program assumptions seed");
    let expected = [
        "A-001", "A-002", "A-003", "A-004", "A-005", "A-006", "A-007", "A-008",
    ];
    assert_eq!(ledger.rows().len(), expected.len());
    for id in expected {
        let id = AssumptionId::try_new(id).expect("seed id");
        let row = ledger.rows().get(&id).expect("required seed row");
        assert!(!row.predicate().is_empty());
        assert!(!row.scope().is_empty());
        assert!(!row.evidence().requirement().is_empty());
        assert!(!row.monitor().signal().is_empty());
        assert!(!row.owner().as_str().is_empty());
    }
}

#[test]
fn closed_case_is_a_canonical_fixed_point_with_a_bound_receipt() {
    let case = closed_case(CaseKnobs::default());
    case.validate().expect("closed physical case validates");

    let bytes = case.canonical_bytes().expect("canonical case bytes");
    let decoded = VvCase::from_canonical_bytes(&bytes).expect("decode canonical case");
    assert_eq!(decoded, case);
    assert_eq!(
        decoded.canonical_bytes().expect("re-encoded case"),
        bytes,
        "decode/encode must be a byte-for-byte fixed point"
    );
    assert_eq!(
        decoded.content_hash().expect("decoded case hash"),
        case.content_hash().expect("source case hash")
    );

    let artifacts = case.artifacts();
    for artifact in &artifacts {
        let artifact_bytes = artifact
            .canonical_bytes()
            .expect("canonical artifact bytes");
        let artifact_round_trip =
            VvArtifact::from_canonical_bytes(&artifact_bytes).expect("decode canonical artifact");
        assert_eq!(&artifact_round_trip, artifact);
        assert_eq!(
            artifact_round_trip
                .canonical_bytes()
                .expect("re-encoded artifact"),
            artifact_bytes
        );
        assert_eq!(
            artifact_round_trip
                .content_hash()
                .expect("round-trip artifact hash"),
            artifact.content_hash().expect("source artifact hash")
        );
    }

    let admitted = case.clone().admit().expect("closed case admits");
    let receipt = admitted.receipt();
    assert_eq!(receipt.schema_version(), VV_SCHEMA_VERSION);
    assert_eq!(receipt.ruleset_version(), VV_RULESET_VERSION);
    assert_eq!(
        receipt.case_hash(),
        admitted.case().content_hash().expect("admitted case hash")
    );
    assert_eq!(receipt.context_id(), admitted.case().context().id());
    assert_eq!(receipt.qois().len(), 1);
    assert_eq!(receipt.artifact_hashes().len(), artifacts.len());
    assert!(receipt.has_valid_binding());
    receipt
        .verify_case(admitted.case())
        .expect("receipt re-verifies its exact case");
}

#[test]
fn schema_admission_receipt_identity_version_domain_fields_and_verification_are_exact() {
    assert_eq!(VV_SCHEMA_ADMISSION_RECEIPT_IDENTITY_VERSION, 2);
    assert_eq!(
        VV_SCHEMA_ADMISSION_RECEIPT_IDENTITY_DOMAIN,
        "org.frankensim.fs-evidence.vv-schema-admission-receipt.v2"
    );
    assert_eq!(
        [
            ArtifactKind::ContextOfUse,
            ArtifactKind::ValidationPlan,
            ArtifactKind::ExperimentArtifact,
            ArtifactKind::CalibrationSplit,
            ArtifactKind::SolutionVerificationReceipt,
            ArtifactKind::PredictionAssessment,
            ArtifactKind::AssumptionsLedger,
        ]
        .map(ArtifactKind::canonical_wire_tag),
        [0, 1, 2, 3, 4, 5, 6],
        "receipt row order is governed by explicit stable wire tags, not enum declaration order",
    );

    let admitted = closed_case(CaseKnobs::default())
        .admit()
        .expect("baseline case admits");
    let receipt = admitted.receipt();
    let canonical_qois = receipt.qois().iter().cloned().collect::<Vec<_>>();
    let canonical_artifacts = receipt
        .artifact_hashes()
        .iter()
        .map(|(key, value)| (key.clone(), *value))
        .collect::<Vec<_>>();
    let baseline_preimage = schema_admission_receipt_preimage(
        receipt.schema_version(),
        receipt.ruleset_version(),
        receipt.case_hash(),
        receipt.context_id(),
        &canonical_qois,
        &canonical_artifacts,
    );
    assert_eq!(
        receipt.receipt_hash(),
        fs_blake3::hash_domain(
            VV_SCHEMA_ADMISSION_RECEIPT_IDENTITY_DOMAIN,
            &baseline_preimage,
        ),
        "the published v2 domain and exact canonical fields reproduce the receipt identity",
    );
    for other_era in [
        "org.frankensim.fs-evidence.vv-schema-admission-receipt.v1",
        "org.frankensim.fs-evidence.vv-schema-admission-receipt.v3",
    ] {
        assert_ne!(
            receipt.receipt_hash(),
            fs_blake3::hash_domain(other_era, &baseline_preimage),
            "the corrected v2 receipt preimage must stay separated from {other_era}",
        );
    }

    let moved = |field: &str, preimage: Vec<u8>| {
        assert_ne!(
            receipt.receipt_hash(),
            fs_blake3::hash_domain(VV_SCHEMA_ADMISSION_RECEIPT_IDENTITY_DOMAIN, &preimage),
            "{field} must move schema-admission receipt identity",
        );
    };
    let first_artifact_kind_tag_offset = 4
        + 4
        + 32
        + 8
        + receipt.context_id().as_str().len()
        + 8
        + canonical_qois
            .iter()
            .map(|qoi| 8 + qoi.as_str().len())
            .sum::<usize>()
        + 8;
    assert_eq!(
        baseline_preimage[first_artifact_kind_tag_offset],
        canonical_artifacts
            .first()
            .expect("closed case has artifact rows")
            .0
            .0
            .canonical_wire_tag(),
        "the stable artifact-kind tag is an explicit receipt-preimage byte",
    );
    let mut changed_kind_tag = baseline_preimage.clone();
    changed_kind_tag[first_artifact_kind_tag_offset] ^= 0x80;
    moved("artifact kind order tag", changed_kind_tag);
    moved(
        "wire schema version",
        schema_admission_receipt_preimage(
            receipt.schema_version() + 1,
            receipt.ruleset_version(),
            receipt.case_hash(),
            receipt.context_id(),
            &canonical_qois,
            &canonical_artifacts,
        ),
    );
    moved(
        "ruleset version",
        schema_admission_receipt_preimage(
            receipt.schema_version(),
            receipt.ruleset_version() + 1,
            receipt.case_hash(),
            receipt.context_id(),
            &canonical_qois,
            &canonical_artifacts,
        ),
    );
    moved(
        "case hash",
        schema_admission_receipt_preimage(
            receipt.schema_version(),
            receipt.ruleset_version(),
            hash("other-case"),
            receipt.context_id(),
            &canonical_qois,
            &canonical_artifacts,
        ),
    );
    moved(
        "context id",
        schema_admission_receipt_preimage(
            receipt.schema_version(),
            receipt.ruleset_version(),
            receipt.case_hash(),
            &artifact_id("other-context"),
            &canonical_qois,
            &canonical_artifacts,
        ),
    );

    let mut changed_qois = canonical_qois.clone();
    changed_qois.push(qoi_id("temperature"));
    changed_qois.sort();
    moved(
        "QoI count and identity",
        schema_admission_receipt_preimage(
            receipt.schema_version(),
            receipt.ruleset_version(),
            receipt.case_hash(),
            receipt.context_id(),
            &changed_qois,
            &canonical_artifacts,
        ),
    );
    let mut reversed_qois = changed_qois.clone();
    reversed_qois.reverse();
    assert_ne!(
        schema_admission_receipt_preimage(
            receipt.schema_version(),
            receipt.ruleset_version(),
            receipt.case_hash(),
            receipt.context_id(),
            &changed_qois,
            &canonical_artifacts,
        ),
        schema_admission_receipt_preimage(
            receipt.schema_version(),
            receipt.ruleset_version(),
            receipt.case_hash(),
            receipt.context_id(),
            &reversed_qois,
            &canonical_artifacts,
        ),
        "QoI canonical order is part of the receipt preimage",
    );

    let mut changed_artifact_hash = canonical_artifacts.clone();
    changed_artifact_hash[0].1 = hash("other-artifact-content");
    moved(
        "artifact hash",
        schema_admission_receipt_preimage(
            receipt.schema_version(),
            receipt.ruleset_version(),
            receipt.case_hash(),
            receipt.context_id(),
            &canonical_qois,
            &changed_artifact_hash,
        ),
    );
    let mut changed_artifact_id = canonical_artifacts.clone();
    changed_artifact_id[0].0.1 = artifact_id("other-artifact-id");
    moved(
        "artifact id",
        schema_admission_receipt_preimage(
            receipt.schema_version(),
            receipt.ruleset_version(),
            receipt.case_hash(),
            receipt.context_id(),
            &canonical_qois,
            &changed_artifact_id,
        ),
    );
    let mut changed_artifact_kind = canonical_artifacts.clone();
    let original_kind = changed_artifact_kind[0].0.0;
    changed_artifact_kind[0].0.0 = if original_kind == ArtifactKind::AssumptionsLedger {
        ArtifactKind::ContextOfUse
    } else {
        ArtifactKind::AssumptionsLedger
    };
    moved(
        "artifact kind",
        schema_admission_receipt_preimage(
            receipt.schema_version(),
            receipt.ruleset_version(),
            receipt.case_hash(),
            receipt.context_id(),
            &canonical_qois,
            &changed_artifact_kind,
        ),
    );
    moved(
        "artifact count",
        schema_admission_receipt_preimage(
            receipt.schema_version(),
            receipt.ruleset_version(),
            receipt.case_hash(),
            receipt.context_id(),
            &canonical_qois,
            &canonical_artifacts[..canonical_artifacts.len() - 1],
        ),
    );
    let mut reversed_artifacts = canonical_artifacts.clone();
    reversed_artifacts.reverse();
    assert_eq!(
        baseline_preimage,
        schema_admission_receipt_preimage(
            receipt.schema_version(),
            receipt.ruleset_version(),
            receipt.case_hash(),
            receipt.context_id(),
            &canonical_qois,
            &reversed_artifacts,
        ),
        "caller or map insertion order is representation-only because receipt rows are re-sorted by stable kind tag and artifact id",
    );
    assert_ne!(
        baseline_preimage,
        schema_admission_receipt_preimage_in_wire_order(
            receipt.schema_version(),
            receipt.ruleset_version(),
            receipt.case_hash(),
            receipt.context_id(),
            &canonical_qois,
            &reversed_artifacts,
        ),
        "the explicit canonical row-order rule itself is identity-bearing",
    );

    let mut changed_length_framing = baseline_preimage.clone();
    changed_length_framing[40] ^= 1;
    moved("length framing", changed_length_framing);
    let mut big_endian_schema = baseline_preimage;
    big_endian_schema[..4].copy_from_slice(&receipt.schema_version().to_be_bytes());
    moved("fixed numeric little endian", big_endian_schema);

    assert!(receipt.has_valid_binding());
    receipt
        .verify_case(admitted.case())
        .expect("receipt verifies the exact admitted case");
    let changed_case = closed_case(CaseKnobs {
        context_domain_lo: 2.0,
        ..CaseKnobs::default()
    });
    assert_only_rule_field(
        receipt.verify_case(&changed_case),
        VvRule::ReceiptBinding,
        "receipt",
    );
}

#[test]
fn canonical_transport_refuses_wrong_domain_and_trailing_bytes() {
    let bytes = closed_case(CaseKnobs::default())
        .canonical_bytes()
        .expect("canonical fixture");

    let mut wrong_domain = bytes.clone();
    wrong_domain[0] ^= 0xff;
    let error = VvCase::from_canonical_bytes(&wrong_domain)
        .expect_err("wrong transport domain must refuse");
    assert_eq!(error.rule_name(), "vv-canonical-identity");

    for stale_schema in [1_u32, 2] {
        let mut stale = bytes.clone();
        stale[4..8].copy_from_slice(&stale_schema.to_le_bytes());
        let error = VvCase::from_canonical_bytes(&stale)
            .expect_err("pre-v3 canonical transport must refuse");
        assert_eq!(error.rule_name(), "vv-canonical-identity");
        assert_eq!(error.offset(), 4);
        assert!(error.detail().contains("unsupported V&V schema version"));
    }

    let mut trailing = bytes;
    trailing.push(0);
    let error = VvCase::from_canonical_bytes(&trailing).expect_err("trailing bytes must refuse");
    assert_eq!(error.rule_name(), "vv-canonical-identity");
}

#[test]
fn canonical_artifact_transport_refuses_unknown_stable_kind_tag_at_exact_byte() {
    let artifact = VvArtifact::from(closed_case(CaseKnobs::default()).context().clone());
    let mut bytes = artifact
        .canonical_bytes()
        .expect("canonical context artifact");
    assert_eq!(&bytes[..4], b"FSVV");
    assert_eq!(bytes[12], 0, "artifact transport root tag");
    assert_eq!(
        bytes[13],
        ArtifactKind::ContextOfUse.canonical_wire_tag(),
        "stable artifact-kind tag immediately follows the root",
    );

    bytes[13] = u8::MAX;
    let error = VvArtifact::from_canonical_bytes(&bytes)
        .expect_err("an unregistered artifact-kind wire tag must refuse");
    assert_eq!(error.offset(), 13);
    assert_eq!(error.detail(), "unknown artifact-kind tag 255");
}

#[test]
fn canonical_case_decode_refuses_forged_derived_applicability() {
    let forged = closed_case(CaseKnobs {
        outside_domain: true,
        force_in_domain_claim: true,
        ..CaseKnobs::default()
    });
    let bytes = forged
        .canonical_bytes()
        .expect("structurally encodable forged case");

    let error = VvCase::from_canonical_bytes(&bytes)
        .expect_err("decode must run closed-case semantic validation");
    assert_eq!(error.rule_name(), "vv-canonical-identity");
    assert!(
        error.to_string().contains("vv-applicability-decision"),
        "semantic refusal must preserve the governing rule: {error}"
    );
}

#[test]
fn derived_applicability_refuses_signed_zero_aliases() {
    let honest = closed_case(CaseKnobs {
        outside_domain: true,
        context_domain_lo: 0.0,
        ..CaseKnobs::default()
    });
    honest
        .validate()
        .expect("bit-exact derived applicability remains valid");

    let forged = closed_case(CaseKnobs {
        outside_domain: true,
        context_domain_lo: 0.0,
        reported_domain_lo: Some(-0.0),
        ..CaseKnobs::default()
    });
    assert_rule(forged.validate(), VvRule::ApplicabilityDecision);

    let bytes = forged
        .canonical_bytes()
        .expect("signed-zero alias remains structurally encodable");
    let error = VvCase::from_canonical_bytes(&bytes)
        .expect_err("canonical decode must reject a derived signed-zero alias");
    assert_eq!(error.rule_name(), "vv-canonical-identity");
    assert!(
        error.to_string().contains("vv-applicability-decision"),
        "semantic refusal must preserve the governing rule: {error}"
    );
}

#[test]
fn validation_metric_specs_refuse_signed_zero_duplicates_and_noncanonical_order() {
    let experiment_ref = reference(ArtifactKind::ExperimentArtifact, "experiment-1");
    let split_ref = reference(ArtifactKind::CalibrationSplit, "split-1");
    assert_rule(
        QoiValidationPlan::try_new(
            qoi_id("length"),
            vec![experiment_ref.clone()],
            split_ref.clone(),
            vec![
                ValidationMetricSpec::NormalizedDiscrepancy { maximum: 0.0 },
                ValidationMetricSpec::NormalizedDiscrepancy { maximum: 1.0 },
                ValidationMetricSpec::NormalizedDiscrepancy { maximum: -0.0 },
            ],
            diagnostic_plan(None),
        ),
        VvRule::ValidationMetricUncertainty,
    );

    let negative_zero = QoiValidationPlan::try_new(
        qoi_id("length"),
        vec![experiment_ref.clone()],
        split_ref.clone(),
        vec![ValidationMetricSpec::NormalizedDiscrepancy { maximum: -0.0 }],
        diagnostic_plan(None),
    )
    .expect("one exact negative-zero threshold remains representable");
    let ValidationMetricSpec::NormalizedDiscrepancy { maximum } = &negative_zero.metrics()[0]
    else {
        panic!("fixture must retain its normalized-discrepancy metric");
    };
    assert_eq!(maximum.to_bits(), (-0.0f64).to_bits());

    let row = QoiValidationPlan::try_new(
        qoi_id("length"),
        vec![experiment_ref],
        split_ref,
        vec![
            ValidationMetricSpec::NormalizedDiscrepancy { maximum: 0.0 },
            ValidationMetricSpec::NormalizedDiscrepancy { maximum: 1.0 },
            ValidationMetricSpec::NormalizedDiscrepancy { maximum: 2.0 },
        ],
        diagnostic_plan(None),
    )
    .expect("canonical metric fixture");
    let plan = ValidationPlan::try_new(
        header("validation-plan-zero", "unitless"),
        reference(ArtifactKind::ContextOfUse, "context-1"),
        vec![row],
    )
    .expect("validation plan fixture");
    let mut bytes = plan.canonical_bytes().expect("canonical validation plan");
    let mut encoded_metric_tail = Vec::new();
    encoded_metric_tail.push(1);
    encoded_metric_tail.extend_from_slice(&1.0f64.to_bits().to_le_bytes());
    encoded_metric_tail.push(1);
    encoded_metric_tail.extend_from_slice(&2.0f64.to_bits().to_le_bytes());
    let offsets = bytes
        .windows(encoded_metric_tail.len())
        .enumerate()
        .filter_map(|(offset, window)| (window == encoded_metric_tail.as_slice()).then_some(offset))
        .collect::<Vec<_>>();
    assert_eq!(
        offsets.len(),
        1,
        "fixture must expose one threshold to mutate"
    );
    let encoded_two_offset = offsets[0] + 1 + std::mem::size_of::<f64>() + 1;
    bytes[encoded_two_offset..encoded_two_offset + std::mem::size_of::<f64>()]
        .copy_from_slice(&(-0.0f64).to_bits().to_le_bytes());

    let error = ValidationPlan::from_canonical_bytes(&bytes)
        .expect_err("decoder must reject a signed-zero alias out of semantic order");
    assert_eq!(error.rule_name(), "vv-canonical-identity");
}

#[test]
fn admission_receipt_refuses_a_semantically_different_case() {
    let case = closed_case(CaseKnobs::default());
    let (_, receipt) = case
        .clone()
        .admit()
        .expect("source case admits")
        .into_parts();
    let changed = closed_case(CaseKnobs {
        predicted: 9.94,
        ..CaseKnobs::default()
    });
    changed.validate().expect("changed case remains admissible");
    assert_rule(receipt.verify_case(&changed), VvRule::ReceiptBinding);
}

#[test]
fn synthetic_or_second_implementation_data_cannot_validate_physics() {
    for origin in [
        OriginKind::SyntheticHighFidelity,
        OriginKind::SecondImplementation,
    ] {
        let case = closed_case(CaseKnobs {
            origin,
            ..CaseKnobs::default()
        });
        assert_rule(case.validate(), VvRule::ValidationRequiresPhysicalReferent);
    }
}

#[test]
fn every_adverse_diagnostic_has_its_exact_refusal_rule() {
    for rule in [
        VvRule::DiagnosticObservability,
        VvRule::DiagnosticIdentifiability,
        VvRule::DiagnosticConfounding,
        VvRule::DiagnosticInverseCrime,
    ] {
        let case = closed_case(CaseKnobs {
            failing_diagnostic: Some(rule),
            ..CaseKnobs::default()
        });
        assert_rule(case.validate(), rule);
    }
}

#[test]
fn process_conformance_cannot_substitute_for_uncertainty_or_validation() {
    let case = closed_case(CaseKnobs {
        process_as_uncertainty: true,
        ..CaseKnobs::default()
    });
    assert_rule(case.validate(), VvRule::ProcessConformanceSeparate);
}

#[test]
fn context_domain_exit_must_demote_instead_of_silently_extrapolating() {
    let honest = closed_case(CaseKnobs {
        outside_domain: true,
        ..CaseKnobs::default()
    });
    honest
        .validate()
        .expect("domain exit with the derived Demoted result remains honest evidence");
    assert!(matches!(
        honest
            .predictions()
            .values()
            .next()
            .expect("prediction")
            .applicability(),
        ApplicabilityDecision::Demoted { .. }
    ));

    let laundered = closed_case(CaseKnobs {
        outside_domain: true,
        force_in_domain_claim: true,
        ..CaseKnobs::default()
    });
    assert_rule(laundered.validate(), VvRule::ApplicabilityDecision);
}

#[test]
fn every_seed_assumption_domain_exit_demotes_or_refuses_by_policy() {
    for id in [
        "A-001", "A-002", "A-003", "A-004", "A-005", "A-006", "A-007", "A-008",
    ] {
        let honest = closed_case(CaseKnobs {
            failing_assumption: Some(id),
            ..CaseKnobs::default()
        });
        honest
            .validate()
            .unwrap_or_else(|error| panic!("honest {id} disposition refused: {error}"));
        let assessment = honest.predictions().values().next().expect("prediction");
        assert!(
            matches!(
                assessment.applicability(),
                ApplicabilityDecision::Demoted { .. } | ApplicabilityDecision::Refused { .. }
            ),
            "{id} failure retained an in-domain claim"
        );

        let laundered = closed_case(CaseKnobs {
            failing_assumption: Some(id),
            force_in_domain_claim: true,
            ..CaseKnobs::default()
        });
        assert_rule(laundered.validate(), VvRule::ApplicabilityDecision);
    }
}

#[test]
fn all_eight_seed_rows_have_individual_semantic_tripwires() {
    let rows = [
        ("A-001", VvRule::AssumptionA001),
        ("A-002", VvRule::AssumptionA002),
        ("A-003", VvRule::AssumptionA003),
        ("A-004", VvRule::AssumptionA004),
        ("A-005", VvRule::AssumptionA005),
        ("A-006", VvRule::AssumptionA006),
        ("A-007", VvRule::AssumptionA007),
        ("A-008", VvRule::AssumptionA008),
    ];
    for (id, rule) in rows {
        let case = closed_case(CaseKnobs {
            tamper_seed: Some(id),
            ..CaseKnobs::default()
        });
        assert_rule(case.validate(), rule);
    }
}

#[test]
fn context_qoi_without_prediction_is_not_a_closed_claim_graph() {
    let case = closed_case(CaseKnobs::default());
    let missing_prediction = VvCase::try_new(
        case.context().clone(),
        case.validation_plan().clone(),
        case.experiments().values().cloned().collect(),
        case.splits().values().cloned().collect(),
        case.solution_verification().values().cloned().collect(),
        Vec::new(),
        case.assumptions().clone(),
    )
    .expect("structurally representable incomplete case");
    assert_rule(missing_prediction.validate(), VvRule::QoiDependencyClosed);
}

#[test]
fn preregistered_metric_outcomes_are_derived_not_asserted() {
    // bead gt1k3: interval agreement is EVALUATED from the artifact
    // numbers — a prediction whose discrepancy exceeds the combined
    // uncertainty refuses even though its metric list is non-empty.
    assert_rule(
        closed_case(CaseKnobs {
            predicted: 10.25,
            ..CaseKnobs::default()
        })
        .validate(),
        VvRule::ValidationMetricUncertainty,
    );
}

#[test]
fn weakened_posterior_threshold_refuses() {
    // bead gt1k3: the check's own threshold must be EXACTLY the
    // preregistered plan value; a weakened copy cannot substitute.
    assert_rule(
        closed_case(CaseKnobs {
            posterior_minimum: 0.001,
            ..CaseKnobs::default()
        })
        .validate(),
        VvRule::ValidationMetricUncertainty,
    );
}

#[test]
fn failed_posterior_check_refuses() {
    // bead gt1k3: PosteriorPredictiveCheck::passed() is REQUIRED — a
    // tail probability below the preregistered minimum refuses.
    assert_rule(
        closed_case(CaseKnobs {
            posterior_tail: 0.01,
            ..CaseKnobs::default()
        })
        .validate(),
        VvRule::ValidationMetricUncertainty,
    );
}

#[test]
fn observation_manifest_refuses_aliased_source_rows() {
    // bead xl3yi/i94v.3.3.1: two ids pointing at one immutable raw locator are
    // unrepresentable even if other provenance is relabelled.
    let aliased = ObservationManifest::try_new(vec![
        (observation_id("cal-1"), manifest_row("row-shared")),
        (observation_id("val-1"), manifest_row("row-shared")),
    ]);
    assert!(aliased.is_err(), "aliasing one source row must refuse");
    for (field, source) in [
        (
            "dataset source bytes",
            ObservationSourceRef::try_new(
                ContentHash([0; 32]),
                "org.frankensim.fs-evidence.test-row-locator.v1",
                1,
                hash("locator"),
                hash("extraction-receipt"),
            ),
        ),
        (
            "locator contract version",
            ObservationSourceRef::try_new(
                hash("source-bytes"),
                "org.frankensim.fs-evidence.test-row-locator.v1",
                0,
                hash("locator"),
                hash("extraction-receipt"),
            ),
        ),
        (
            "locator hash",
            ObservationSourceRef::try_new(
                hash("source-bytes"),
                "org.frankensim.fs-evidence.test-row-locator.v1",
                1,
                ContentHash([0; 32]),
                hash("extraction-receipt"),
            ),
        ),
        (
            "extraction receipt",
            ObservationSourceRef::try_new(
                hash("source-bytes"),
                "org.frankensim.fs-evidence.test-row-locator.v1",
                1,
                hash("locator"),
                ContentHash([0; 32]),
            ),
        ),
        (
            "locator domain",
            ObservationSourceRef::try_new(
                hash("source-bytes"),
                "unknown",
                1,
                hash("locator"),
                hash("extraction-receipt"),
            ),
        ),
    ] {
        assert!(source.is_err(), "invalid {field} must refuse");
    }
    // Genuinely distinct replicate rows with equal VALUES stay valid:
    // identity is the full typed source, never the row value.
    ObservationManifest::try_new(vec![
        (observation_id("rep-1"), manifest_row("locator-run-1")),
        (observation_id("rep-2"), manifest_row("locator-run-2")),
    ])
    .expect("distinct locators admit regardless of value equality");

    // A new extraction receipt remains artifact-identity-bearing, but cannot
    // relabel one immutable raw locator into a second observation or holdout.
    let receipt_relabel = ObservationManifest::try_new(vec![
        (
            observation_id("extract-1"),
            manifest_row_with_source(
                observation_source_ref_with(
                    "source-bytes",
                    "org.frankensim.fs-evidence.test-row-locator.v1",
                    1,
                    "shared-locator",
                    "extraction-receipt-1",
                ),
                "length",
                "instrument-1",
                "channel-1",
                "clock-1",
            ),
        ),
        (
            observation_id("extract-2"),
            manifest_row_with_source(
                observation_source_ref_with(
                    "source-bytes",
                    "org.frankensim.fs-evidence.test-row-locator.v1",
                    1,
                    "shared-locator",
                    "extraction-receipt-2",
                ),
                "length",
                "instrument-1",
                "channel-1",
                "clock-1",
            ),
        ),
    ]);
    assert_rule(receipt_relabel, VvRule::SchemaIdentity);
}

#[test]
fn observation_manifest_locator_injectivity_is_metadata_and_receipt_independent() {
    let shared_source = observation_source_ref_with(
        "source-bytes",
        "org.frankensim.fs-evidence.test-row-locator.v1",
        1,
        "shared-locator",
        "extraction-receipt-1",
    );
    let baseline_row = manifest_row_with_source(
        shared_source.clone(),
        "length",
        "instrument-1",
        "channel-1",
        "clock-1",
    );
    let receipt_relabel = observation_source_ref_with(
        "source-bytes",
        "org.frankensim.fs-evidence.test-row-locator.v1",
        1,
        "shared-locator",
        "extraction-receipt-2",
    );
    assert_eq!(
        shared_source.locator_identity(),
        receipt_relabel.locator_identity(),
        "extraction evidence must not mint a second raw locator identity",
    );

    for (label, relabelled_row) in [
        (
            "QoI",
            manifest_row_with_source(
                shared_source.clone(),
                "temperature",
                "instrument-1",
                "channel-1",
                "clock-1",
            ),
        ),
        (
            "instrument",
            manifest_row_with_source(
                shared_source.clone(),
                "length",
                "instrument-2",
                "channel-1",
                "clock-1",
            ),
        ),
        (
            "acquisition channel",
            manifest_row_with_source(
                shared_source.clone(),
                "length",
                "instrument-1",
                "channel-2",
                "clock-1",
            ),
        ),
        (
            "clock",
            manifest_row_with_source(
                shared_source.clone(),
                "length",
                "instrument-1",
                "channel-1",
                "clock-2",
            ),
        ),
        (
            "extraction receipt",
            manifest_row_with_source(
                receipt_relabel,
                "length",
                "instrument-1",
                "channel-1",
                "clock-1",
            ),
        ),
    ] {
        let manifest = ObservationManifest::try_new(vec![
            (observation_id("row-1"), baseline_row.clone()),
            (observation_id("row-2"), relabelled_row),
        ]);
        let error = manifest.expect_err("one raw locator cannot be relabelled");
        assert_eq!(
            error.violations().len(),
            1,
            "{label} relabelling must have one isolated refusal: {error}",
        );
        assert_eq!(error.violations()[0].rule(), VvRule::SchemaIdentity);
        assert_eq!(error.violations()[0].field(), "experiment.manifest");
    }

    assert_only_rule_field(
        ObservationManifest::try_new(vec![
            (observation_id("duplicate-id"), baseline_row),
            (
                observation_id("duplicate-id"),
                manifest_row_with_source(
                    observation_source_ref("distinct-locator"),
                    "length",
                    "instrument-1",
                    "channel-1",
                    "clock-1",
                ),
            ),
        ]),
        VvRule::SchemaIdentity,
        "experiment.manifest",
    );
}

#[test]
fn observation_manifest_bare_locator_hash_is_scoped_by_dataset_domain_and_version() {
    let baseline = observation_source_ref_with(
        "source-bytes",
        "org.frankensim.fs-evidence.test-row-locator.v1",
        1,
        "shared-bare-locator-hash",
        "extraction-receipt-1",
    );
    for (label, scoped_source) in [
        (
            "dataset bytes",
            observation_source_ref_with(
                "other-source-bytes",
                "org.frankensim.fs-evidence.test-row-locator.v1",
                1,
                "shared-bare-locator-hash",
                "extraction-receipt-2",
            ),
        ),
        (
            "locator domain",
            observation_source_ref_with(
                "source-bytes",
                "org.frankensim.fs-evidence.other-row-locator.v1",
                1,
                "shared-bare-locator-hash",
                "extraction-receipt-2",
            ),
        ),
        (
            "locator contract version",
            observation_source_ref_with(
                "source-bytes",
                "org.frankensim.fs-evidence.test-row-locator.v1",
                2,
                "shared-bare-locator-hash",
                "extraction-receipt-2",
            ),
        ),
    ] {
        assert_eq!(
            baseline.locator_hash(),
            scoped_source.locator_hash(),
            "the fixture isolates {label} from the bare locator digest",
        );
        assert_ne!(
            baseline.locator_identity(),
            scoped_source.locator_identity(),
            "{label} is part of the scoped raw-locator identity",
        );
        ObservationManifest::try_new(vec![
            (
                observation_id("row-1"),
                manifest_row_with_source(
                    baseline.clone(),
                    "length",
                    "instrument-1",
                    "channel-1",
                    "clock-1",
                ),
            ),
            (
                observation_id("row-2"),
                manifest_row_with_source(
                    scoped_source,
                    "length",
                    "instrument-1",
                    "channel-1",
                    "clock-1",
                ),
            ),
        ])
        .unwrap_or_else(|error| panic!("{label} must scope otherwise equal locators: {error}"));
    }
}

#[test]
fn canonical_experiment_decode_preserves_duplicate_raw_locator_refusal() {
    let artifact = experiment(
        "duplicate-raw-locator-transport",
        ExperimentOrigin::Physical {
            apparatus_id: artifact_id("apparatus-1"),
            facility_id: artifact_id("facility-1"),
        },
        true,
        true,
    )
    .expect("baseline experiment has distinct raw locators");
    let mut bytes = artifact
        .canonical_bytes()
        .expect("baseline experiment canonical bytes");
    let validation_locator = hash("row-val-1");
    let calibration_locator = hash("row-cal-1");
    let locator_offset = bytes
        .windows(validation_locator.as_bytes().len())
        .position(|window| window == validation_locator.as_bytes())
        .expect("canonical manifest carries the validation locator hash");
    bytes[locator_offset..locator_offset + validation_locator.as_bytes().len()]
        .copy_from_slice(calibration_locator.as_bytes());

    let error = ExperimentArtifact::from_canonical_bytes(&bytes)
        .expect_err("a second id cannot smuggle an existing raw locator through transport");
    assert_eq!(error.rule_name(), "vv-canonical-identity");
    assert!(
        error
            .detail()
            .contains("experiment manifest refused by the model")
            && error.detail().contains(VvRule::SchemaIdentity.slug())
            && error
                .detail()
                .contains("distinct observation ids cannot alias one immutable raw locator"),
        "canonical refusal must retain the exact model rule and cause: {error}",
    );
}

#[test]
fn observation_manifest_refuses_oversized_canonical_hash_preimage() {
    let padding = "x".repeat(220);
    let locator_domain = format!("locator-{padding}");
    let qoi = format!("qoi-{padding}");
    let instrument = format!("instrument-{padding}");
    let channel = format!("channel-{padding}");
    let clock = format!("clock-{padding}");
    let rows = (0..MAX_VV_ITEMS)
        .map(|index| {
            let row_id = format!("row-{index:04}-{padding}");
            let source = ObservationSourceRef::try_new(
                hash("oversized-source-bytes"),
                locator_domain.clone(),
                1,
                hash(&format!("locator-{index}")),
                hash(&format!("receipt-{index}")),
            )
            .expect("bounded row source");
            (
                observation_id(&row_id),
                manifest_row_with_source(source, &qoi, &instrument, &channel, &clock),
            )
        })
        .collect();
    assert_rule(
        ObservationManifest::try_new(rows),
        VvRule::SchemaCardinality,
    );
}

#[test]
fn blind_row_repointed_at_seen_source_refuses() {
    // bead xl3yi: the id sets still cover exactly, but blind-1 binds
    // cal-1's SOURCE row — the case-level manifest cross-check refuses.
    assert_rule(
        closed_case(CaseKnobs {
            blind_source_label: "row-cal-1",
            ..CaseKnobs::default()
        })
        .validate(),
        VvRule::SplitBlindHoldoutSealed,
    );
}

#[test]
fn blind_commitment_binds_source_identities() {
    // bead xl3yi: identical id partitions over different source rows
    // seal DIFFERENT commitments.
    let split_a =
        split_with(vec!["cal-1"], vec!["val-1"], vec!["blind-1"]).expect("baseline split");
    let split_b = CalibrationSplit::try_new(
        header("split-1", "unitless"),
        reference(ArtifactKind::ExperimentArtifact, "experiment-1"),
        hash("preregistration"),
        vec![observation_id("cal-1")],
        vec![observation_id("val-1")],
        vec![(observation_id("blind-1"), hash("row-blind-ALTERNATE"))],
    )
    .expect("repointed split");
    assert_ne!(
        split_a.blind_commitment(),
        split_b.blind_commitment(),
        "source identity is commitment-bearing"
    );
}

#[test]
fn blind_holdout_identity_version_domain_and_fields_are_exact() {
    assert_eq!(VV_BLIND_HOLDOUT_IDENTITY_VERSION, 2);
    assert_eq!(
        VV_BLIND_HOLDOUT_IDENTITY_DOMAIN,
        "org.frankensim.fs-evidence.vv-blind-holdout.v2"
    );

    let baseline = CalibrationSplit::try_new(
        header("blind-identity", "unitless"),
        reference(ArtifactKind::ExperimentArtifact, "experiment-1"),
        hash("blind-preregistration"),
        vec![observation_id("cal-1")],
        vec![observation_id("val-1")],
        vec![
            (observation_id("blind-a"), hash("source-a")),
            (observation_id("blind-b"), hash("source-b")),
        ],
    )
    .expect("baseline blind-holdout commitment");
    let canonical_rows = baseline
        .blind_sources()
        .iter()
        .map(|(id, source)| (id.clone(), *source))
        .collect::<Vec<_>>();
    assert_eq!(
        baseline.blind_commitment(),
        blind_holdout_identity_for(
            VV_BLIND_HOLDOUT_IDENTITY_DOMAIN,
            baseline.preregistration_hash(),
            &canonical_rows,
        ),
        "the published v2 domain and exact canonical preimage reproduce the commitment",
    );
    assert_ne!(
        baseline.blind_commitment(),
        blind_holdout_identity_for(
            "org.frankensim.fs-evidence.vv-blind-holdout.v3",
            baseline.preregistration_hash(),
            &canonical_rows,
        ),
        "a new blind-holdout identity era must rotate the domain",
    );

    let preregistration_variant = CalibrationSplit::try_new(
        header("blind-identity", "unitless"),
        reference(ArtifactKind::ExperimentArtifact, "experiment-1"),
        hash("other-preregistration"),
        vec![observation_id("cal-1")],
        vec![observation_id("val-1")],
        canonical_rows.clone(),
    )
    .expect("preregistration identity variant");
    let observation_id_variant = CalibrationSplit::try_new(
        header("blind-identity", "unitless"),
        reference(ArtifactKind::ExperimentArtifact, "experiment-1"),
        baseline.preregistration_hash(),
        vec![observation_id("cal-1")],
        vec![observation_id("val-1")],
        vec![
            (observation_id("blind-a"), hash("source-a")),
            (observation_id("blind-c"), hash("source-b")),
        ],
    )
    .expect("observation identity variant");
    let source_variant = CalibrationSplit::try_new(
        header("blind-identity", "unitless"),
        reference(ArtifactKind::ExperimentArtifact, "experiment-1"),
        baseline.preregistration_hash(),
        vec![observation_id("cal-1")],
        vec![observation_id("val-1")],
        vec![
            (observation_id("blind-a"), hash("source-a")),
            (observation_id("blind-b"), hash("source-c")),
        ],
    )
    .expect("source locator variant");
    let row_count_variant = CalibrationSplit::try_new(
        header("blind-identity", "unitless"),
        reference(ArtifactKind::ExperimentArtifact, "experiment-1"),
        baseline.preregistration_hash(),
        vec![observation_id("cal-1")],
        vec![observation_id("val-1")],
        vec![(observation_id("blind-a"), hash("source-a"))],
    )
    .expect("blind row-count variant");
    for (field, commitment) in [
        (
            "preregistration hash",
            preregistration_variant.blind_commitment(),
        ),
        (
            "observation id and length framing",
            observation_id_variant.blind_commitment(),
        ),
        ("source locator hash", source_variant.blind_commitment()),
        ("blind row count", row_count_variant.blind_commitment()),
    ] {
        assert_ne!(
            baseline.blind_commitment(),
            commitment,
            "{field} must move the blind-holdout commitment",
        );
    }

    let mut reversed_rows = canonical_rows.clone();
    reversed_rows.reverse();
    assert_ne!(
        blind_holdout_identity_for(
            VV_BLIND_HOLDOUT_IDENTITY_DOMAIN,
            baseline.preregistration_hash(),
            &canonical_rows,
        ),
        blind_holdout_identity_for(
            VV_BLIND_HOLDOUT_IDENTITY_DOMAIN,
            baseline.preregistration_hash(),
            &reversed_rows,
        ),
        "canonical row order is commitment-bearing",
    );
    let input_permuted = CalibrationSplit::try_new(
        header("blind-identity", "unitless"),
        reference(ArtifactKind::ExperimentArtifact, "experiment-1"),
        baseline.preregistration_hash(),
        vec![observation_id("cal-1")],
        vec![observation_id("val-1")],
        reversed_rows,
    )
    .expect("permuted constructor input");
    assert_eq!(
        baseline.blind_commitment(),
        input_permuted.blind_commitment(),
        "constructor input order must normalize to the canonical observation-id order",
    );
}

#[test]
fn blind_holdout_identity_ignores_noncommitment_split_fields() {
    let preregistration = hash("blind-noncommitment-preregistration");
    let blind_sources = vec![
        (observation_id("blind-a"), hash("source-a")),
        (observation_id("blind-b"), hash("source-b")),
    ];
    let commitment = |split: Result<CalibrationSplit, VvErrors>| {
        split
            .expect("noncommitment variant remains a valid calibration split")
            .blind_commitment()
    };
    let baseline = commitment(CalibrationSplit::try_new(
        header("blind-noncommitment-a", "unitless"),
        reference(ArtifactKind::ExperimentArtifact, "experiment-a"),
        preregistration,
        vec![observation_id("cal-a")],
        vec![observation_id("val-a")],
        blind_sources.clone(),
    ));

    for (field, variant) in [
        (
            "CalibrationSplit.header",
            commitment(CalibrationSplit::try_new(
                header("blind-noncommitment-b", "m"),
                reference(ArtifactKind::ExperimentArtifact, "experiment-a"),
                preregistration,
                vec![observation_id("cal-a")],
                vec![observation_id("val-a")],
                blind_sources.clone(),
            )),
        ),
        (
            "CalibrationSplit.experiment",
            commitment(CalibrationSplit::try_new(
                header("blind-noncommitment-a", "unitless"),
                reference(ArtifactKind::ExperimentArtifact, "experiment-b"),
                preregistration,
                vec![observation_id("cal-a")],
                vec![observation_id("val-a")],
                blind_sources.clone(),
            )),
        ),
        (
            "CalibrationSplit.calibration",
            commitment(CalibrationSplit::try_new(
                header("blind-noncommitment-a", "unitless"),
                reference(ArtifactKind::ExperimentArtifact, "experiment-a"),
                preregistration,
                vec![observation_id("cal-b")],
                vec![observation_id("val-a")],
                blind_sources.clone(),
            )),
        ),
        (
            "CalibrationSplit.validation",
            commitment(CalibrationSplit::try_new(
                header("blind-noncommitment-a", "unitless"),
                reference(ArtifactKind::ExperimentArtifact, "experiment-a"),
                preregistration,
                vec![observation_id("cal-a")],
                vec![observation_id("val-b")],
                blind_sources,
            )),
        ),
    ] {
        assert_eq!(
            baseline, variant,
            "{field} is part of complete split identity but not the narrower blind-release commitment",
        );
    }
}

#[test]
fn calibration_split_refuses_missing_preregistration_in_model_and_transport() {
    assert_rule(
        CalibrationSplit::try_new(
            header("split-zero-preregistration", "unitless"),
            reference(ArtifactKind::ExperimentArtifact, "experiment-1"),
            ContentHash([0; 32]),
            vec![observation_id("cal-1")],
            vec![observation_id("val-1")],
            vec![(observation_id("blind-1"), hash("row-blind-1"))],
        ),
        VvRule::SplitBlindHoldoutSealed,
    );

    let split = split_with(vec!["cal-1"], vec!["val-1"], vec!["blind-1"])
        .expect("baseline preregistered split");
    let mut bytes = split.canonical_bytes().expect("split canonical bytes");
    let preregistration = hash("preregistration");
    let offset = bytes
        .windows(32)
        .position(|window| window == preregistration.as_bytes())
        .expect("canonical split carries preregistration identity");
    bytes[offset..offset + 32].fill(0);
    let error = CalibrationSplit::from_canonical_bytes(&bytes)
        .expect_err("transport cannot smuggle a zero preregistration identity");
    assert!(
        error
            .detail()
            .contains("calibration split refused by the model")
            && error.detail().contains("split.preregistration_hash"),
        "expected model-level split refusal after decode, got {error}",
    );
}
