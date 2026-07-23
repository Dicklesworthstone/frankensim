//! G0/G3 battery for the evidence-bearing V&V corpus schema.

use fs_evidence::{ColorRank, NumericalKind};
use fs_qty::{Dims, QtyAny};
use fs_vvreg::ContentHash;
use fs_vvreg::corpus::{
    AcceptanceRecord, AcquisitionProvenance, Availability, CalibrationRecord, ContextRange,
    ContextValue, CorpusArtifact, CorpusDataset, CorpusEnvelope, CorpusError, CorpusLicense,
    CorpusQueryRefusal, CorpusRegistry, DatasetDraft, DatasetField, DatasetPartition,
    EnvironmentCondition, EvidenceLevel, GeometryRecord, LEVEL_C_COOLING_QOIS,
    MeasurementUncertainty, PayloadRetention, PreprocessingLineage, PreprocessingStep,
    RedistributionPolicy, RetentionClass, RetentionPolicy, SensorPlacement, SensorRecord,
    admit_dataset, corpus,
};
use fs_vvreg::partition::{DatasetPurpose, PartitionLedger, PartitionRefusal};
use fs_vvreg::thermal_level_a::thermal_level_a_cases;

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

fn available_calibration_mut(sensor: &mut SensorRecord) -> &mut CalibrationRecord {
    let Availability::Available(calibration) = &mut sensor.calibration else {
        panic!("test draft must carry available calibration")
    };
    calibration
}

fn available_placement_mut(sensor: &mut SensorRecord) -> &mut SensorPlacement {
    let Availability::Available(placement) = &mut sensor.placement else {
        panic!("test draft must carry available placement")
    };
    placement
}

fn complete_preprocessing_mut(draft: &mut DatasetDraft) -> &mut Vec<PreprocessingStep> {
    let Some(PreprocessingLineage::Complete(steps)) = &mut draft.preprocessing else {
        panic!("test draft must carry complete preprocessing")
    };
    steps
}

#[test]
fn every_top_level_mandatory_field_has_a_typed_refusal() {
    macro_rules! missing {
        ($field:ident, $expected:expr) => {{
            let mut draft = complete_draft("missing-probe");
            draft.$field = None;
            assert_eq!(
                admit_dataset(draft),
                Err(CorpusError::MissingField { field: $expected })
            );
        }};
    }

    missing!(id, DatasetField::Id);
    missing!(title, DatasetField::Title);
    missing!(raw_payload, DatasetField::RawPayload);
    missing!(sensors, DatasetField::Sensors);
    missing!(geometry, DatasetField::Geometry);
    missing!(environment, DatasetField::Environment);
    missing!(partition, DatasetField::Partition);
    missing!(preprocessing, DatasetField::Preprocessing);
    missing!(final_artifact, DatasetField::FinalArtifact);
    missing!(context_of_use, DatasetField::ContextOfUse);
    missing!(license, DatasetField::License);
    missing!(provenance, DatasetField::Provenance);
    missing!(retention, DatasetField::Retention);
    missing!(acceptance_envelopes, DatasetField::AcceptanceEnvelopes);
    missing!(evidence_level, DatasetField::EvidenceLevel);
}

#[test]
fn nested_sensor_calibration_placement_and_uncertainty_fail_closed() {
    let mut missing_calibration = complete_draft("bad-calibration");
    available_calibration_mut(&mut missing_calibration.sensors.as_mut().unwrap()[0])
        .certificate_id = " ".to_string();
    assert!(matches!(
        admit_dataset(missing_calibration),
        Err(CorpusError::InvalidField {
            field: DatasetField::Sensors,
            ..
        })
    ));

    let mut bad_placement = complete_draft("bad-placement");
    available_placement_mut(&mut bad_placement.sensors.as_mut().unwrap()[0]).uncertainty[0] =
        QtyAny::new(-1.0, LENGTH);
    assert!(matches!(
        admit_dataset(bad_placement),
        Err(CorpusError::InvalidField {
            field: DatasetField::Sensors,
            ..
        })
    ));

    let mut bad_uncertainty = complete_draft("bad-uncertainty");
    bad_uncertainty.sensors.as_mut().unwrap()[0].uncertainty = MeasurementUncertainty::Bounded {
        half_width: QtyAny::new(0.2, LENGTH),
    };
    assert!(matches!(
        admit_dataset(bad_uncertainty),
        Err(CorpusError::InvalidField {
            field: DatasetField::Sensors,
            ..
        })
    ));
}

#[test]
fn covariance_uses_squared_quantity_dimensions() {
    let mut valid = complete_draft("covariance-valid");
    valid.sensors.as_mut().unwrap()[0].uncertainty = MeasurementUncertainty::CovarianceDiagonal {
        variance: QtyAny::new(0.04, Dims([0, 0, 0, 2, 0, 0])),
    };
    assert!(admit_dataset(valid).is_ok());

    let mut invalid = complete_draft("covariance-invalid");
    invalid.sensors.as_mut().unwrap()[0].uncertainty = MeasurementUncertainty::CovarianceDiagonal {
        variance: QtyAny::new(0.04, TEMPERATURE),
    };
    assert!(matches!(
        admit_dataset(invalid),
        Err(CorpusError::InvalidField {
            field: DatasetField::Sensors,
            ..
        })
    ));
}

#[test]
fn raw_to_final_preprocessing_lineage_is_exact_and_gap_free() {
    let mut bad_input = complete_draft("bad-lineage-input");
    complete_preprocessing_mut(&mut bad_input)[0].input = hash(9);
    assert_eq!(
        admit_dataset(bad_input),
        Err(CorpusError::BrokenLineage {
            step: 0,
            reason: "input hash does not equal the preceding retained artifact"
        })
    );

    let mut bad_final = complete_draft("bad-lineage-final");
    bad_final.final_artifact = Some(hash(9));
    assert_eq!(
        admit_dataset(bad_final),
        Err(CorpusError::BrokenLineage {
            step: 1,
            reason: "final artifact does not equal the last transform output"
        })
    );
}

#[test]
fn retention_cannot_drop_raw_or_calibration_evidence() {
    for drop_raw in [true, false] {
        let mut draft = complete_draft(if drop_raw {
            "drop-raw"
        } else {
            "drop-calibration"
        });
        let retention = draft.retention.as_mut().unwrap();
        if drop_raw {
            retention.preserve_raw = false;
        } else {
            retention.preserve_calibration = false;
        }
        assert!(matches!(
            admit_dataset(draft),
            Err(CorpusError::InvalidField {
                field: DatasetField::Retention,
                ..
            })
        ));
    }
}

#[test]
fn explicit_unstated_uncertainty_demotes_experimental_data() {
    let admitted = admit_dataset(complete_draft("stated-uncertainty")).unwrap();
    assert_eq!(admitted.physical_claim_cap(), ColorRank::Validated);

    let mut unstated = complete_draft("unstated-uncertainty");
    unstated.sensors.as_mut().unwrap()[0].uncertainty = MeasurementUncertainty::Unstated;
    let admitted = admit_dataset(unstated).unwrap();
    assert_eq!(admitted.physical_claim_cap(), ColorRank::Estimated);

    let mut cross_code = complete_draft("cross-code");
    cross_code.evidence_level = Some(EvidenceLevel::CrossCode);
    let admitted = admit_dataset(cross_code).unwrap();
    assert_eq!(admitted.physical_claim_cap(), ColorRank::Estimated);
}

#[test]
fn every_explicit_authority_gap_demotes_experimental_data() {
    let assert_demoted = |draft: DatasetDraft| {
        assert_eq!(
            admit_dataset(draft).unwrap().physical_claim_cap(),
            ColorRank::Estimated
        );
    };

    let mut derived = complete_draft("derived-only");
    derived.raw_payload = Some(PayloadRetention::DerivedOnly {
        retained: artifact(1, "data/probe/raw.csv"),
        reason: "original acquisition unavailable".to_string(),
    });
    assert_demoted(derived);

    let mut geometry = complete_draft("geometry-unavailable");
    geometry.geometry = Some(Availability::Unavailable {
        reason: "geometry records unavailable".to_string(),
    });
    assert_demoted(geometry);

    let mut environment = complete_draft("environment-unavailable");
    environment.environment = Some(Availability::Unavailable {
        reason: "environment records unavailable".to_string(),
    });
    assert_demoted(environment);

    let mut lineage = complete_draft("lineage-unreplayable");
    lineage.preprocessing = Some(PreprocessingLineage::Unreplayable {
        retained_input: hash(1),
        retained_output: hash(1),
        reason: "transform tool unavailable".to_string(),
    });
    assert_demoted(lineage);

    let mut license = complete_draft("license-unavailable");
    license.license = Some(Availability::Unavailable {
        reason: "redistribution terms unresolved".to_string(),
    });
    assert_demoted(license);

    let mut date = complete_draft("date-unavailable");
    date.provenance.as_mut().unwrap().measured_on = Availability::Unavailable {
        reason: "acquisition date unavailable".to_string(),
    };
    assert_demoted(date);

    let mut instrument = complete_draft("instrument-unavailable");
    instrument.sensors.as_mut().unwrap()[0].instrument_id = Availability::Unavailable {
        reason: "instrument identity unavailable".to_string(),
    };
    assert_demoted(instrument);

    let mut calibration = complete_draft("calibration-unavailable");
    calibration.sensors.as_mut().unwrap()[0].calibration = Availability::Unavailable {
        reason: "calibration unavailable".to_string(),
    };
    assert_demoted(calibration);

    let mut placement = complete_draft("placement-unavailable");
    placement.sensors.as_mut().unwrap()[0].placement = Availability::Unavailable {
        reason: "placement unavailable".to_string(),
    };
    assert_demoted(placement);

    let mut envelope = complete_draft("envelope-unpinned");
    envelope.acceptance_envelopes.as_mut().unwrap()[0].envelope = CorpusEnvelope::Unpinned {
        basis: "no scalar acceptance rule is defensible".to_string(),
    };
    assert_demoted(envelope);
}

#[test]
fn canonical_dataset_round_trip_is_bit_exact_and_tamper_evident() {
    let dataset = admit_dataset(complete_draft("round-trip")).unwrap();
    let bytes = dataset.encode();
    assert_eq!(CorpusDataset::decode(&bytes).unwrap(), dataset);
    assert_eq!(CorpusDataset::decode(&bytes).unwrap().encode(), bytes);

    let mut bad_magic = bytes.clone();
    bad_magic[0] ^= 1;
    assert_eq!(
        CorpusDataset::decode(&bad_magic),
        Err(CorpusError::BadMagic)
    );

    let mut trailing = bytes;
    trailing.push(0);
    assert_eq!(
        CorpusDataset::decode(&trailing),
        Err(CorpusError::TrailingBytes { count: 1 })
    );
}

#[test]
fn caller_registry_is_deterministic_but_has_no_query_authority() {
    let a = complete_draft("dataset-a");
    let b = complete_draft("dataset-b");
    let forward = CorpusRegistry::build(vec![a.clone(), b.clone()]).unwrap();
    let reverse = CorpusRegistry::build(vec![b, a]).unwrap();
    assert_eq!(forward.digest(), reverse.digest());
    let partitions = PartitionLedger::capture(&forward);
    assert!(matches!(
        forward.query(
            &partitions,
            "dataset-a",
            DatasetPurpose::Validation,
            &[ContextValue {
                name: "ambient_temperature".to_string(),
                value: QtyAny::new(300.0, TEMPERATURE),
            }],
        ),
        Err(PartitionRefusal::Corpus(
            CorpusQueryRefusal::UnauthoritativeRegistry
        ))
    ));
}

#[test]
fn seeded_query_requires_exact_partition_and_complete_in_domain_context() {
    let context = [ContextValue {
        name: "reference_cost_work_units".to_string(),
        value: QtyAny::dimensionless(250.0),
    }];
    let partitions = PartitionLedger::capture(corpus());
    let evidence = corpus()
        .query(
            &partitions,
            "fs-benchmark-cht-query-v1",
            DatasetPurpose::Validation,
            &context,
        )
        .unwrap();
    assert_eq!(evidence.value.id(), "fs-benchmark-cht-query-v1");
    assert_eq!(evidence.value.evidence_level(), EvidenceLevel::CrossCode);
    assert_eq!(evidence.value.physical_claim_cap(), ColorRank::Estimated);
    assert_eq!(evidence.numerical.kind, NumericalKind::NoClaim);
    assert!(evidence.statistical.rel_width(1.0).is_infinite());
    assert_eq!(
        evidence.model.validity.bound("reference_cost_work_units"),
        Some((250.0, 250.0))
    );
    assert!(evidence.model.discrepancy_rel.is_infinite());

    assert!(matches!(
        corpus().query(
            &partitions,
            "fs-benchmark-cht-query-v1",
            DatasetPurpose::Calibration,
            &context,
        ),
        Err(PartitionRefusal::PurposeMismatch {
            dataset_id,
            declared: DatasetPartition::Validation,
            attempted: DatasetPurpose::Calibration,
        }) if dataset_id == "fs-benchmark-cht-query-v1"
    ));
    assert!(matches!(
        corpus().query(
            &partitions,
            "fs-benchmark-cht-query-v1",
            DatasetPurpose::Validation,
            &[],
        ),
        Err(PartitionRefusal::Corpus(CorpusQueryRefusal::MissingContext {
            name,
        })) if name == "reference_cost_work_units"
    ));
}

#[test]
fn context_refusals_distinguish_range_dimension_unknown_and_duplicate_errors() {
    let partitions = PartitionLedger::capture(corpus());
    let out_of_range = [ContextValue {
        name: "reference_cost_work_units".to_string(),
        value: QtyAny::dimensionless(500.0),
    }];
    assert!(matches!(
        corpus().query(
            &partitions,
            "fs-benchmark-cht-query-v1",
            DatasetPurpose::Validation,
            &out_of_range,
        ),
        Err(PartitionRefusal::Corpus(
            CorpusQueryRefusal::OutOfContext { .. }
        ))
    ));

    let wrong_dims = [ContextValue {
        name: "reference_cost_work_units".to_string(),
        value: QtyAny::new(1.0, LENGTH),
    }];
    assert!(matches!(
        corpus().query(
            &partitions,
            "fs-benchmark-cht-query-v1",
            DatasetPurpose::Validation,
            &wrong_dims,
        ),
        Err(PartitionRefusal::Corpus(
            CorpusQueryRefusal::ContextDimensionMismatch { .. }
        ))
    ));

    let unknown = [ContextValue {
        name: "wind_speed".to_string(),
        value: QtyAny::new(0.0, Dims([1, 0, -1, 0, 0, 0])),
    }];
    assert!(matches!(
        corpus().query(
            &partitions,
            "fs-benchmark-cht-query-v1",
            DatasetPurpose::Validation,
            &unknown,
        ),
        Err(PartitionRefusal::Corpus(CorpusQueryRefusal::UnknownContext {
            name,
        })) if name == "wind_speed"
    ));

    let duplicate = [
        ContextValue {
            name: "reference_cost_work_units".to_string(),
            value: QtyAny::dimensionless(250.0),
        },
        ContextValue {
            name: "reference_cost_work_units".to_string(),
            value: QtyAny::dimensionless(250.0),
        },
    ];
    assert!(matches!(
        corpus().query(
            &partitions,
            "fs-benchmark-cht-query-v1",
            DatasetPurpose::Validation,
            &duplicate,
        ),
        Err(PartitionRefusal::Corpus(CorpusQueryRefusal::DuplicateContext {
            name,
        })) if name == "reference_cost_work_units"
    ));
}

#[test]
fn seeded_raw_payload_is_bound_to_the_retained_file() {
    let bytes = include_bytes!("../../../data/vv-corpus/fs-benchmark-cht-query-v1/raw-sensors.csv");
    let dataset = corpus().dataset("fs-benchmark-cht-query-v1").unwrap();
    assert_eq!(
        dataset.raw_payload().artifact().digest,
        fs_blake3::hash_bytes(bytes)
    );
    assert_eq!(
        dataset.raw_payload().artifact().byte_len,
        bytes.len() as u64
    );
}

#[test]
fn martin_moyce_is_retained_without_inventing_missing_authority() {
    let bytes = include_bytes!("../../../data/reference/martin-moyce-1952.jsonl");
    let dataset = corpus().dataset("martin-moyce-1952-square-column").unwrap();
    assert_eq!(
        dataset.raw_payload().artifact().digest,
        fs_blake3::hash_bytes(bytes)
    );
    assert!(matches!(
        dataset.raw_payload(),
        PayloadRetention::DerivedOnly { .. }
    ));
    assert_eq!(dataset.evidence_level(), EvidenceLevel::PublishedExperiment);
    assert_eq!(dataset.physical_claim_cap(), ColorRank::Estimated);

    let evidence = corpus()
        .query(
            &PartitionLedger::capture(corpus()),
            "martin-moyce-1952-square-column",
            DatasetPurpose::Validation,
            &[ContextValue {
                name: "t_star".to_string(),
                value: QtyAny::dimensionless(1.0),
            }],
        )
        .unwrap();
    assert_eq!(evidence.numerical.kind, NumericalKind::NoClaim);
    assert!(evidence.statistical.rel_width(1.0).is_infinite());
    assert_eq!(evidence.model.validity.bound("t_star"), Some((0.41, 2.95)));
    assert!(evidence.model.discrepancy_rel.is_infinite());
}

#[test]
fn published_electronics_thermal_level_c_rows_are_retained_and_fail_closed() {
    let pires_source =
        include_bytes!("../../../data/vv-corpus/level-c/pires-fonseca-2024/source.pdf");
    let pires_final =
        include_bytes!("../../../data/vv-corpus/level-c/pires-fonseca-2024/digitized.tsv");
    let nunes_source = include_bytes!("../../../data/vv-corpus/level-c/nunes-2023/source.pdf");
    let nunes_final = include_bytes!("../../../data/vv-corpus/level-c/nunes-2023/digitized.tsv");
    let markal_source =
        include_bytes!("../../../data/vv-corpus/level-c/markal-kul-2026/source.pdf");
    let markal_final =
        include_bytes!("../../../data/vv-corpus/level-c/markal-kul-2026/supplementary.zip");

    for (id, source, final_artifact) in [
        (
            "pires-fonseca-2024-flat-strip-fins",
            pires_source.as_slice(),
            pires_final.as_slice(),
        ),
        (
            "nunes-2023-micro-pin-fin",
            nunes_source.as_slice(),
            nunes_final.as_slice(),
        ),
        (
            "markal-kul-2026-fin-distribution",
            markal_final.as_slice(),
            markal_final.as_slice(),
        ),
    ] {
        let dataset = corpus().dataset(id).unwrap();
        assert_eq!(dataset.evidence_level(), EvidenceLevel::PublishedExperiment);
        assert_eq!(dataset.partition(), DatasetPartition::Validation);
        assert_eq!(dataset.physical_claim_cap(), ColorRank::Estimated);
        assert!(matches!(
            dataset.raw_payload(),
            PayloadRetention::DerivedOnly { .. }
        ));
        assert!(matches!(dataset.geometry(), Availability::Available(_)));
        assert!(matches!(dataset.environment(), Availability::Available(_)));
        assert!(matches!(dataset.license(), Availability::Available(_)));
        assert!(matches!(
            dataset.preprocessing(),
            PreprocessingLineage::Complete(_)
        ));
        assert_eq!(
            dataset.raw_payload().artifact().digest,
            fs_blake3::hash_bytes(source)
        );
        assert_eq!(
            dataset.raw_payload().artifact().byte_len,
            source.len() as u64
        );
        assert_eq!(
            dataset.final_artifact(),
            fs_blake3::hash_bytes(final_artifact)
        );
        assert!(!dataset.acceptance_envelopes().is_empty());
    }

    let markal = corpus()
        .dataset("markal-kul-2026-fin-distribution")
        .unwrap();
    let Availability::Available(markal_geometry) = markal.geometry() else {
        panic!("Markal-Kul nominal geometry must be retained")
    };
    assert_eq!(
        markal_geometry.nominal.digest,
        fs_blake3::hash_bytes(markal_source)
    );

    let acquisition_log = include_str!("../../../data/vv-corpus/level-c/acquisition-log-v1.tsv");
    assert_eq!(
        acquisition_log
            .lines()
            .filter(|line| line.starts_with("admit\t"))
            .count(),
        3
    );
    assert_eq!(
        acquisition_log
            .lines()
            .filter(|line| line.starts_with("reject\t"))
            .count(),
        5
    );
}

fn tsv_row<'a>(tsv: &'a str, series: &str, x: &str) -> Vec<&'a str> {
    tsv.lines()
        .skip(1)
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .find(|columns| columns[0] == series && columns[1] == x)
        .unwrap_or_else(|| panic!("missing digitized row {series}/{x}"))
}

#[test]
fn digitization_subsamples_stay_within_declared_half_widths() {
    let pires = include_str!("../../../data/vv-corpus/level-c/pires-fonseca-2024/digitized.tsv");
    let pires_row = tsv_row(pires, "flat-plate-fins", "2000");
    let stored_re = pires_row[1].parse::<f64>().unwrap();
    let stored_nu = pires_row[2].parse::<f64>().unwrap();
    let re_half_width = pires_row[3].parse::<f64>().unwrap();
    let nu_half_width = pires_row[4].parse::<f64>().unwrap();
    // Independent cursor placement on the retained Figure 7 source.
    let redigitized_re = 2_025.0;
    let redigitized_nu = 5.85;
    assert!((stored_re - redigitized_re).abs() <= re_half_width);
    assert!((stored_nu - redigitized_nu).abs() <= nu_half_width);
    let published_flat_plate_correlation = 0.112 * stored_re.powf(0.52);
    assert!((stored_nu - published_flat_plate_correlation).abs() <= nu_half_width);

    let nunes = include_str!("../../../data/vv-corpus/level-c/nunes-2023/digitized.tsv");
    let nunes_row = tsv_row(nunes, "S1-G1000-subcooling20", "3.2");
    let stored_superheat = nunes_row[1].parse::<f64>().unwrap();
    let stored_heat_flux = nunes_row[2].parse::<f64>().unwrap();
    let superheat_half_width = nunes_row[3].parse::<f64>().unwrap();
    let heat_flux_half_width = nunes_row[4].parse::<f64>().unwrap();
    // Independent cursor placement on the retained Figure 6a source.
    let redigitized_superheat = 3.0;
    let redigitized_heat_flux = 50.0;
    assert!((stored_superheat - redigitized_superheat).abs() <= superheat_half_width);
    assert!((stored_heat_flux - redigitized_heat_flux).abs() <= heat_flux_half_width);
}

#[test]
fn audit_is_deterministic_and_warns_for_seed_gaps() {
    let audit = corpus().audit();
    assert!(audit.is_clean());
    assert_eq!(audit.rows().len(), 5 + thermal_level_a_cases().len());
    for row in audit.rows() {
        assert_eq!(row.mandatory_present(), 15);
        assert_eq!(row.mandatory_total(), 15);
        assert_eq!(row.status(), "WARN");
    }
    assert!(audit.errors().is_empty());
    let rendered = audit.render_table();
    assert!(
        rendered.contains(
            "fs-benchmark-cht-query-v1 | 15/15 | 0/2 | validation | B | estimated | WARN"
        )
    );
    assert!(rendered.contains(
        "martin-moyce-1952-square-column | 15/15 | 0/2 | validation | C | estimated | WARN"
    ));
    assert!(rendered.contains(
        "pires-fonseca-2024-flat-strip-fins | 15/15 | 0/2 | validation | C | estimated | WARN"
    ));
    assert!(
        rendered
            .contains("nunes-2023-micro-pin-fin | 15/15 | 0/2 | validation | C | estimated | WARN")
    );
    assert!(rendered.contains(
        "markal-kul-2026-fin-distribution | 15/15 | 0/2 | validation | C | estimated | WARN"
    ));
    assert!(
        rendered.contains("dataset=martin-moyce-1952-square-column claim_gap=raw_payload.original")
    );
    assert!(rendered.contains(
        "dataset=martin-moyce-1952-square-column claim_gap=acceptance.surge-front-position-z"
    ));
    assert_eq!(audit.qoi_coverage().len(), LEVEL_C_COOLING_QOIS.len());
    for (qoi, expected_datasets) in [
        ("average-nusselt-number", 2),
        ("component-peak-temperature", 0),
        ("convective-thermal-resistance", 1),
        ("effective-heat-flux", 1),
        ("friction-factor", 1),
        ("pressure-drop", 2),
        ("temperature-nonuniformity", 0),
        ("thermal-interface-resistance", 0),
    ] {
        let row = audit
            .qoi_coverage()
            .iter()
            .find(|row| row.qoi() == qoi)
            .unwrap();
        assert_eq!(row.level_c_datasets(), expected_datasets);
        assert_eq!(row.is_covered(), expected_datasets != 0);
        let status = if expected_datasets == 0 {
            "GAP"
        } else {
            "COVERED"
        };
        assert!(rendered.contains(&format!("{qoi} | {expected_datasets} | {status}")));
    }
    for qoi in [
        "component-peak-temperature",
        "temperature-nonuniformity",
        "thermal-interface-resistance",
    ] {
        let row = audit
            .qoi_coverage()
            .iter()
            .find(|row| row.qoi() == qoi)
            .unwrap();
        assert!(!row.is_covered());
        assert!(rendered.contains(&format!("{qoi} | 0 | GAP")));
        assert!(rendered.contains(&format!("qoi_gap={qoi} evidence_level=C datasets=0")));
    }
    assert_eq!(corpus().audit().render_table(), rendered);
}
