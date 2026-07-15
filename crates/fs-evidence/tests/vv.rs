use fs_blake3::ContentHash;
use fs_evidence::vv::*;

fn hash(label: &str) -> ContentHash {
    fs_blake3::hash_domain("org.frankensim.fs-evidence.vv-test.v1", label.as_bytes())
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
            (observation_id("cal-1"), hash("row-cal-1")),
            (observation_id("val-1"), hash("row-val-1")),
            (observation_id("blind-1"), hash("row-blind-1")),
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
    let split_ref = ArtifactRef::new(
        ArtifactKind::CalibrationSplit,
        split.id().clone(),
        hash("split-1"),
    );

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
    let split_ref = ArtifactRef::new(
        ArtifactKind::CalibrationSplit,
        split.id().clone(),
        hash("split-1"),
    );
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
fn canonical_transport_refuses_wrong_domain_and_trailing_bytes() {
    let bytes = closed_case(CaseKnobs::default())
        .canonical_bytes()
        .expect("canonical fixture");

    let mut wrong_domain = bytes.clone();
    wrong_domain[0] ^= 0xff;
    let error = VvCase::from_canonical_bytes(&wrong_domain)
        .expect_err("wrong transport domain must refuse");
    assert_eq!(error.rule_name(), "vv-canonical-identity");

    let mut trailing = bytes;
    trailing.push(0);
    let error = VvCase::from_canonical_bytes(&trailing).expect_err("trailing bytes must refuse");
    assert_eq!(error.rule_name(), "vv-canonical-identity");
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
    // bead xl3yi: two ids pointing at ONE source row are
    // unrepresentable at the experiment.
    let aliased = ObservationManifest::try_new(vec![
        (observation_id("cal-1"), hash("row-shared")),
        (observation_id("val-1"), hash("row-shared")),
    ]);
    assert!(aliased.is_err(), "aliasing one source row must refuse");
    let zero = ObservationManifest::try_new(vec![(observation_id("cal-1"), ContentHash([0; 32]))]);
    assert!(zero.is_err(), "the all-zero source sentinel must refuse");
    // Genuinely distinct replicate rows with equal VALUES stay valid:
    // identity is the source locator, never the row value.
    ObservationManifest::try_new(vec![
        (observation_id("rep-1"), hash("locator-run-1")),
        (observation_id("rep-2"), hash("locator-run-2")),
    ])
    .expect("distinct locators admit regardless of value equality");
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
