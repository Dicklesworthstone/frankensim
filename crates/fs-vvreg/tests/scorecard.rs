//! G0/G5 battery for the deterministic public V&V scorecard.
//!
//! The synthetic-corpus tests pin known-error cell arithmetic, loud NO-DATA
//! semantics, reference-uncertainty rendering, and fail-closed refusals. The
//! real-corpus test is the e2e lane: it logs per-cell inputs plus the gap
//! list and pins the honesty header of the committed artifact. Determinism
//! (G5) is byte-identity of both renders across independent builds.

use fs_qty::{Dims, QtyAny};
use fs_vvreg::ContentHash;
use fs_vvreg::adversarial::{AdversarialOutcome, DominantUncertainty, adversarial_registry};
use fs_vvreg::corpus::{
    AcceptanceRecord, AcquisitionProvenance, Availability, CalibrationRecord, ContextRange,
    CorpusArtifact, CorpusEnvelope, CorpusLicense, CorpusRegistry, DatasetDraft, DatasetPartition,
    EnvironmentCondition, EvidenceLevel, GeometryRecord, LEVEL_C_COOLING_QOIS,
    MeasurementUncertainty, PayloadRetention, PreprocessingLineage, PreprocessingStep,
    RedistributionPolicy, RetentionClass, RetentionPolicy, SensorPlacement, SensorRecord, corpus,
};
use fs_vvreg::scorecard::{
    FalseAcceptanceCell, MAX_SCORECARD_RUN_RECORDS, ReferenceUncertainty, ScorecardError,
    ScorecardRunRecord, build_scorecard,
};

const TEMPERATURE: Dims = Dims([0, 0, 0, 1, 0, 0]);
const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);

fn hash(byte: u8) -> ContentHash {
    ContentHash([byte; 32])
}

fn artifact(byte: u8, locator: &str) -> CorpusArtifact {
    CorpusArtifact {
        digest: hash(byte),
        byte_len: 64,
        media_type: "text/csv".to_string(),
        locator: locator.to_string(),
    }
}

fn complete_draft(id: &str) -> DatasetDraft {
    let raw = hash(1);
    DatasetDraft {
        id: Some(id.to_string()),
        title: Some("Complete published-experiment probe".to_string()),
        raw_payload: Some(PayloadRetention::OriginalRaw(artifact(
            1,
            "data/probe/raw.csv",
        ))),
        sensors: Some(vec![SensorRecord {
            id: "temperature-sensor".to_string(),
            instrument_id: Availability::Available("instrument-serial-42".to_string()),
            raw_channel: "temperature".to_string(),
            quantity_dims: TEMPERATURE,
            calibration: Availability::Available(CalibrationRecord {
                certificate_id: "calibration-2026".to_string(),
                certificate_hash: hash(2),
                issued_on: "2026-01-02".to_string(),
                valid_through: Some("2027-01-02".to_string()),
            }),
            placement: Availability::Available(SensorPlacement {
                frame: "probe-frame".to_string(),
                coordinates: [
                    QtyAny::new(0.0, LENGTH),
                    QtyAny::new(0.1, LENGTH),
                    QtyAny::new(0.2, LENGTH),
                ],
                uncertainty: [
                    QtyAny::new(1e-4, LENGTH),
                    QtyAny::new(1e-4, LENGTH),
                    QtyAny::new(1e-4, LENGTH),
                ],
            }),
            uncertainty: MeasurementUncertainty::Bounded {
                half_width: QtyAny::new(0.2, TEMPERATURE),
            },
        }]),
        geometry: Some(Availability::Available(GeometryRecord {
            nominal: artifact(3, "data/probe/nominal.txt"),
            as_built: Some(artifact(4, "data/probe/as-built.txt")),
            frame: "probe-frame".to_string(),
        })),
        environment: Some(Availability::Available(vec![EnvironmentCondition {
            name: "ambient_temperature".to_string(),
            value: QtyAny::new(298.15, TEMPERATURE),
            uncertainty: QtyAny::new(0.1, TEMPERATURE),
        }])),
        partition: Some(DatasetPartition::Validation),
        preprocessing: Some(PreprocessingLineage::Complete(vec![PreprocessingStep {
            ordinal: 0,
            operation: "identity-import".to_string(),
            version: "1".to_string(),
            input: raw,
            output: raw,
        }])),
        final_artifact: Some(raw),
        context_of_use: Some(vec![ContextRange {
            name: "ambient_temperature".to_string(),
            lo: QtyAny::new(290.0, TEMPERATURE),
            hi: QtyAny::new(310.0, TEMPERATURE),
        }]),
        license: Some(Availability::Available(CorpusLicense {
            identifier: "CC-BY-4.0".to_string(),
            terms: "Attribution required".to_string(),
            redistribution: RedistributionPolicy::Allowed,
        })),
        provenance: Some(AcquisitionProvenance {
            measured_by: "Metrology Team".to_string(),
            organization: "Probe Laboratory".to_string(),
            measured_on: Availability::Available("2026-02-03".to_string()),
            source_record: "lab-book-42".to_string(),
        }),
        retention: Some(RetentionPolicy {
            class: RetentionClass::Years(20),
            preserve_raw: true,
            preserve_calibration: true,
            policy_id: "lab-retention-v1".to_string(),
        }),
        acceptance_envelopes: Some(vec![AcceptanceRecord {
            metric: "surface_temperature".to_string(),
            dims: TEMPERATURE,
            envelope: CorpusEnvelope::Tolerance {
                atol: 0.5,
                rtol: 0.01,
            },
            regime: vec![ContextRange {
                name: "ambient_temperature".to_string(),
                lo: QtyAny::new(295.0, TEMPERATURE),
                hi: QtyAny::new(305.0, TEMPERATURE),
            }],
        }]),
        evidence_level: Some(EvidenceLevel::PublishedExperiment),
    }
}

fn synthetic_registry() -> CorpusRegistry {
    let experiment = complete_draft("synthetic-experiment");
    let mut analytic = complete_draft("synthetic-analytic");
    analytic.evidence_level = Some(EvidenceLevel::Analytic);
    CorpusRegistry::build(vec![experiment, analytic]).expect("synthetic drafts admit cleanly")
}

fn run(
    dataset_id: &str,
    predicted: f64,
    uncertainty: ReferenceUncertainty,
    identity_byte: u8,
) -> ScorecardRunRecord {
    ScorecardRunRecord::try_new(
        dataset_id,
        "surface_temperature",
        predicted,
        300.0,
        uncertainty,
        hash(identity_byte),
    )
    .expect("test run record is valid")
}

#[test]
fn known_errors_populate_the_declared_cell() {
    let registry = synthetic_registry();
    let runs = [
        run(
            "synthetic-experiment",
            300.5,
            ReferenceUncertainty::Bounded { half_width: 0.25 },
            7,
        ),
        run(
            "synthetic-experiment",
            310.0,
            ReferenceUncertainty::Unstated,
            8,
        ),
    ];
    let scorecard = build_scorecard(&registry, adversarial_registry(), &runs, &[])
        .expect("synthetic scorecard builds");

    assert!(!scorecard.corpus_seeded());
    let cells = scorecard.cells();
    assert_eq!(cells.len(), 1);
    let cell = &cells[0];
    assert_eq!(cell.qoi(), "surface_temperature");
    assert_eq!(cell.regime(), "ambient_temperature in [295, 305] K");
    assert_eq!(cell.dataset_ids().len(), 2);
    assert_eq!(cell.external_datasets(), 1);
    assert_eq!(cell.runs().len(), 2);
    assert_eq!(cell.envelope_pass(), 1);
    assert_eq!(cell.envelope_fail(), 1);
    assert_eq!(cell.envelope_unpinned(), 0);
    assert_eq!(cell.runs()[0].signed_error(), 0.5);
    assert_eq!(cell.runs()[1].signed_error(), 10.0);

    let markdown = scorecard.render_markdown();
    assert!(markdown.contains("n=2 min=0.5 max=10 max_abs=10"));
    assert!(markdown.contains("pass=1 fail=1 unpinned=0"));
}

#[test]
fn empty_cell_renders_no_data_not_zero() {
    let registry = synthetic_registry();
    let scorecard = build_scorecard(&registry, adversarial_registry(), &[], &[])
        .expect("empty scorecard builds");

    let cell = &scorecard.cells()[0];
    assert!(cell.runs().is_empty());
    assert_eq!(cell.false_acceptance(), FalseAcceptanceCell::NoData);

    let markdown = scorecard.render_markdown();
    assert!(markdown.contains("| NO-DATA | NO-DATA | NO-DATA | NO-DATA |"));
    assert!(markdown.contains("false_acceptance_total: NO-DATA (0 executed challenges)"));
    // The prediction-error column must be NO-DATA, never a zero-count
    // summary ("| n=0 ...") and never a zero false-acceptance claim.
    assert!(!markdown.contains("| n=0"));
    assert!(!markdown.contains("0 of 0 executed"));
    assert!(!markdown.contains("false_acceptance_total: 0"));

    let json = scorecard.render_json();
    assert!(json.contains("\"envelope\":{\"status\":\"no-data\"}"));
    assert!(json.contains("\"false_acceptance\":{\"status\":\"no-data\"}"));
    assert!(json.contains("\"false_acceptance_total\":{\"status\":\"no-data\"}"));
    assert!(json.contains("\"interval_coverage\":{\"status\":\"no-data\""));
}

#[test]
fn reference_uncertainty_is_rendered_with_the_error() {
    let registry = synthetic_registry();
    let runs = [
        run(
            "synthetic-experiment",
            300.5,
            ReferenceUncertainty::Bounded { half_width: 0.25 },
            7,
        ),
        run(
            "synthetic-analytic",
            310.0,
            ReferenceUncertainty::Unstated,
            8,
        ),
    ];
    let scorecard = build_scorecard(&registry, adversarial_registry(), &runs, &[])
        .expect("synthetic scorecard builds");

    let markdown = scorecard.render_markdown();
    assert!(markdown.contains("### Run detail"));
    assert!(markdown.contains("+/-0.25"));
    assert!(markdown.contains("uncertainty unstated"));

    let json = scorecard.render_json();
    assert!(json.contains(
        "\"reference_uncertainty\":{\"status\":\"bounded\",\"half_width\":{\"display\":\"0.25\""
    ));
    assert!(json.contains("\"reference_uncertainty\":{\"status\":\"unstated\"}"));
}

#[test]
fn scorecard_regenerates_byte_identically() {
    let build = || {
        let registry = synthetic_registry();
        let runs = [
            run(
                "synthetic-experiment",
                300.5,
                ReferenceUncertainty::Bounded { half_width: 0.25 },
                7,
            ),
            run(
                "synthetic-analytic",
                310.0,
                ReferenceUncertainty::Unstated,
                8,
            ),
        ];
        let assessment = adversarial_registry()
            .assess(
                "biot-extremes-lumped-breakdown",
                AdversarialOutcome::Refused {
                    dominant: DominantUncertainty::SpatialTemperature,
                },
            )
            .expect("assessment admits");
        let scorecard = build_scorecard(&registry, adversarial_registry(), &runs, &[assessment])
            .expect("synthetic scorecard builds");
        (
            scorecard.render_markdown(),
            scorecard.render_json(),
            scorecard.identity(),
        )
    };
    let (markdown_a, json_a, identity_a) = build();
    let (markdown_b, json_b, identity_b) = build();
    assert_eq!(markdown_a, markdown_b);
    assert_eq!(json_a, json_b);
    assert_eq!(identity_a, identity_b);
}

#[test]
fn real_corpus_e2e_logs_cells_and_gaps() {
    let scorecard = build_scorecard(corpus(), adversarial_registry(), &[], &[])
        .expect("real-corpus scorecard builds");

    assert!(scorecard.corpus_seeded());
    assert!(!scorecard.cells().is_empty());
    for cell in scorecard.cells() {
        println!(
            "cell qoi={} regime={} refs={} external={} runs={}",
            cell.qoi(),
            cell.regime(),
            cell.dataset_ids().len(),
            cell.external_datasets(),
            cell.runs().len()
        );
    }
    for gap in scorecard.known_gaps() {
        println!("gap {gap}");
    }

    // Every scoped cooling QoI is visible somewhere: as a populated cell or
    // as a loud gap row. Gap cells are impossible to hide.
    for qoi in LEVEL_C_COOLING_QOIS {
        let in_cells = scorecard.cells().iter().any(|cell| cell.qoi() == *qoi);
        let in_gaps = scorecard
            .known_gaps()
            .iter()
            .any(|gap| gap.starts_with(&format!("qoi={qoi} ")));
        assert!(
            in_cells || in_gaps,
            "cooling QoI {qoi} must appear as a cell or a gap"
        );
    }
    // The current corpus has analytic-only cells, so the gap list must be
    // loud rather than empty.
    assert!(!scorecard.known_gaps().is_empty());

    let markdown = scorecard.render_markdown();
    assert!(markdown.contains("corpus_authority: seeded"));
    assert!(markdown.contains("interval_coverage: NO-DATA"));
    assert!(markdown.contains("## Known gaps"));
    assert!(markdown.contains("## Adversarial regime limitations"));
    assert!(markdown.contains("false_acceptance_total: NO-DATA (0 executed challenges)"));
}

#[test]
fn seeded_false_acceptance_surfaces_in_the_right_cell() {
    let assessment = adversarial_registry()
        .assess(
            "contact-dominated-two-layer-stack",
            AdversarialOutcome::Prediction {
                absolute_error: 5.0,
                allowed_error: 1.0,
                dominant: DominantUncertainty::ContactResistance,
            },
        )
        .expect("assessment admits");
    assert!(assessment.is_false_acceptance());

    let scorecard = build_scorecard(corpus(), adversarial_registry(), &[], &[assessment])
        .expect("real-corpus scorecard builds");

    assert_eq!(scorecard.executed_assessments(), 1);
    assert_eq!(scorecard.false_acceptance_total(), 1);

    let bound_cells: Vec<_> = scorecard
        .cells()
        .iter()
        .filter(|cell| {
            cell.dataset_ids()
                .iter()
                .any(|id| id == "thermal-a-contact-series")
        })
        .collect();
    assert!(
        !bound_cells.is_empty(),
        "the retained challenge dataset must own at least one cell"
    );
    for cell in &bound_cells {
        assert_eq!(
            cell.false_acceptance(),
            FalseAcceptanceCell::Counted {
                executed: 1,
                false_acceptances: 1,
            },
            "false acceptance must surface in cell qoi={} regime={}",
            cell.qoi(),
            cell.regime()
        );
    }
    let unbound = scorecard
        .cells()
        .iter()
        .find(|cell| {
            !cell
                .dataset_ids()
                .iter()
                .any(|id| id == "thermal-a-contact-series")
        })
        .expect("the corpus has cells beyond the challenge dataset");
    assert_eq!(unbound.false_acceptance(), FalseAcceptanceCell::NoData);

    let markdown = scorecard.render_markdown();
    assert!(markdown.contains("false_acceptance_total: 1 of 1 executed"));
    assert!(markdown.contains("false_acceptance_count: 1"));
    assert!(markdown.contains("| 1 of 1 executed |"));
}

#[test]
fn run_records_fail_closed() {
    let registry = synthetic_registry();

    assert_eq!(
        ScorecardRunRecord::try_new(
            "synthetic-experiment",
            "surface_temperature",
            f64::NAN,
            300.0,
            ReferenceUncertainty::Unstated,
            hash(9),
        ),
        Err(ScorecardError::InvalidRunField {
            field: "predicted",
            reason: "must be finite",
        })
    );
    assert_eq!(
        ScorecardRunRecord::try_new(
            "synthetic-experiment",
            "surface_temperature",
            300.0,
            300.0,
            ReferenceUncertainty::Bounded { half_width: -1.0 },
            hash(9),
        ),
        Err(ScorecardError::InvalidRunField {
            field: "reference_uncertainty",
            reason: "half-width must be finite and non-negative",
        })
    );
    assert_eq!(
        ScorecardRunRecord::try_new(
            "",
            "surface_temperature",
            300.0,
            300.0,
            ReferenceUncertainty::Unstated,
            hash(9),
        ),
        Err(ScorecardError::InvalidRunField {
            field: "dataset_id",
            reason: "is blank",
        })
    );

    let unknown_dataset = ScorecardRunRecord::try_new(
        "no-such-dataset",
        "surface_temperature",
        300.0,
        300.0,
        ReferenceUncertainty::Unstated,
        hash(9),
    )
    .expect("locally valid record");
    assert_eq!(
        build_scorecard(&registry, adversarial_registry(), &[unknown_dataset], &[]),
        Err(ScorecardError::UnknownDataset {
            dataset_id: "no-such-dataset".to_string(),
        })
    );

    let unknown_metric = ScorecardRunRecord::try_new(
        "synthetic-experiment",
        "no-such-metric",
        300.0,
        300.0,
        ReferenceUncertainty::Unstated,
        hash(9),
    )
    .expect("locally valid record");
    assert_eq!(
        build_scorecard(&registry, adversarial_registry(), &[unknown_metric], &[]),
        Err(ScorecardError::UnknownMetric {
            dataset_id: "synthetic-experiment".to_string(),
            metric: "no-such-metric".to_string(),
        })
    );

    let template = run(
        "synthetic-experiment",
        300.5,
        ReferenceUncertainty::Unstated,
        7,
    );
    let too_many = vec![template; MAX_SCORECARD_RUN_RECORDS + 1];
    assert_eq!(
        build_scorecard(&registry, adversarial_registry(), &too_many, &[]),
        Err(ScorecardError::ResourceLimit {
            limit: MAX_SCORECARD_RUN_RECORDS,
            observed: MAX_SCORECARD_RUN_RECORDS + 1,
        })
    );
}

#[test]
fn foreign_and_duplicate_assessments_are_refused() {
    let registry = synthetic_registry();
    let assessment = adversarial_registry()
        .assess(
            "biot-extremes-lumped-breakdown",
            AdversarialOutcome::Refused {
                dominant: DominantUncertainty::SpatialTemperature,
            },
        )
        .expect("assessment admits");

    assert_eq!(
        build_scorecard(
            &registry,
            adversarial_registry(),
            &[],
            &[assessment.clone(), assessment.clone()]
        ),
        Err(ScorecardError::DuplicateAssessment {
            case_id: "biot-extremes-lumped-breakdown".to_string(),
        })
    );

    let small_registry =
        fs_vvreg::adversarial::AdversarialRegistry::build(vec![adversarial_registry().cases()[0]])
            .expect("single-case registry builds");
    assert_eq!(
        build_scorecard(&registry, &small_registry, &[], &[assessment]),
        Err(ScorecardError::ForeignAssessment {
            case_id: "biot-extremes-lumped-breakdown".to_string(),
        })
    );
}
