//! Level-E rig ingest: channel matching, calibration, the energy-balance
//! gate, and blind sealing (bead `frankensim-extreal-program-f85xj.4.5.1`;
//! software child of `frankensim-extreal-program-f85xj.4.5`).
//!
//! Every fixture here is SYNTHETIC. That is the point: the bead's caveat is
//! that the software half must not wait on hardware, so the pipeline is
//! exercised against runs whose physics we control exactly and whose failure
//! modes we can inject deliberately.

use fs_qty::{Dims, QtyAny};
use fs_vvreg::corpus::{
    Availability, CalibrationRecord, DatasetPartition, MeasurementUncertainty, SensorPlacement,
    SensorRecord,
};
use fs_vvreg::rig::{
    ChannelRole, CoolantProperties, EnergyBalance, InstrumentSpec, MAX_INSTRUMENTS,
    MAX_RETAINED_RIG_BYTES, Reading, RigError, RigRun, RigSpec, ingest,
};

const TEMPERATURE: Dims = Dims([0, 0, 0, 1, 0, 0]);
const POWER: Dims = Dims([2, 1, -3, 0, 0, 0]);
const VOLUME_FLOW: Dims = Dims([3, 0, -1, 0, 0, 0]);

// Air at room conditions, and a run that closes exactly:
//   1.2 * 1005 * 0.01 * 10 = 120.6 W
const RHO: f64 = 1.2;
const CP: f64 = 1005.0;
const FLOW: f64 = 0.01;
const INLET: f64 = 295.0;
const OUTLET: f64 = 305.0;
const EXACT_POWER: f64 = RHO * CP * FLOW * (OUTLET - INLET);

fn certificate(id: &str) -> CalibrationRecord {
    CalibrationRecord {
        certificate_id: id.to_string(),
        certificate_hash: fs_blake3::hash_domain("test.calibration", id.as_bytes()),
        issued_on: "2026-01-15".to_string(),
        valid_through: Some("2027-01-15".to_string()),
    }
}

fn sensor(id: &str, dims: Dims, half_width: Option<f64>) -> SensorRecord {
    SensorRecord {
        id: id.to_string(),
        instrument_id: Availability::Available(format!("instrument-{id}")),
        raw_channel: format!("ch_{id}"),
        quantity_dims: dims,
        calibration: Availability::Available(certificate(id)),
        placement: Availability::Unavailable {
            reason: "synthetic rig fixture declares no placement".to_string(),
        },
        uncertainty: match half_width {
            Some(value) => MeasurementUncertainty::Bounded {
                half_width: QtyAny::new(value, dims),
            },
            None => MeasurementUncertainty::Unstated,
        },
    }
}

fn instrument(id: &str, role: ChannelRole, dims: Dims, half_width: Option<f64>) -> InstrumentSpec {
    InstrumentSpec::new(sensor(id, dims, half_width), role).expect("instrument admits")
}

/// Tight bands, so an injected error is unambiguously a violation.
fn tight_instruments() -> Vec<InstrumentSpec> {
    vec![
        instrument("power", ChannelRole::HeaterPower, POWER, Some(0.1)),
        instrument(
            "t-in",
            ChannelRole::InletTemperature,
            TEMPERATURE,
            Some(0.05),
        ),
        instrument(
            "t-out",
            ChannelRole::OutletTemperature,
            TEMPERATURE,
            Some(0.05),
        ),
        instrument("flow", ChannelRole::VolumeFlow, VOLUME_FLOW, Some(1.0e-5)),
    ]
}

fn coolant() -> CoolantProperties {
    CoolantProperties::new(RHO, CP, 0.001).expect("coolant admits")
}

fn spec_with(instruments: Vec<InstrumentSpec>) -> RigSpec {
    RigSpec::new("bench-rig-a", coolant(), instruments).expect("rig admits")
}

fn run_with(power: f64, partition: DatasetPartition) -> RigRun {
    RigRun {
        run_id: "run-0001".to_string(),
        acquired_on: "2026-06-15".to_string(),
        partition,
        readings: vec![
            Reading {
                sensor_id: "power".to_string(),
                value_si: power,
            },
            Reading {
                sensor_id: "t-in".to_string(),
                value_si: INLET,
            },
            Reading {
                sensor_id: "t-out".to_string(),
                value_si: OUTLET,
            },
            Reading {
                sensor_id: "flow".to_string(),
                value_si: FLOW,
            },
        ],
    }
}

// ---------------------------------------------------------------------------
// The happy path and the physical gate.
// ---------------------------------------------------------------------------

#[test]
fn a_closing_run_ingests_and_retains_its_balance() {
    let spec = spec_with(tight_instruments());
    let run = run_with(EXACT_POWER, DatasetPartition::Validation);
    let ingested = ingest(&spec, &run).expect("a closing run ingests");

    assert_eq!(ingested.run_id(), "run-0001");
    assert_eq!(ingested.rig_id(), "bench-rig-a");
    assert_eq!(ingested.partition(), DatasetPartition::Validation);
    assert!(!ingested.sealed());
    assert_eq!(ingested.calibrated_channels(), 4);
    assert_eq!(ingested.acquired_on(), "2026-06-15");
    assert_eq!(ingested.readings().expect("unsealed readings").len(), 4);
    let restored = fs_vvreg::rig::IngestedRun::from_retained_bytes(ingested.retained_bytes())
        .expect("canonical retained bytes round-trip");
    assert_eq!(restored, ingested);

    match ingested.balance() {
        EnergyBalance::Balanced {
            power_in_w,
            power_interval_w,
            heat_out_w,
            heat_interval_w,
            residual_w,
            residual_interval_w,
        } => {
            assert!((power_in_w - EXACT_POWER).abs() < 1e-12);
            assert!((heat_out_w - EXACT_POWER).abs() < 1e-9);
            assert!(residual_w.abs() < 1e-9, "an exact fixture must close");
            assert!(power_interval_w.contains(EXACT_POWER));
            assert!(heat_interval_w.contains(EXACT_POWER));
            assert!(residual_interval_w.contains(0.0));
        }
        other => panic!("expected a balanced verdict, got {other:?}"),
    }
}

#[test]
fn a_run_that_does_not_close_refuses_and_names_the_channels() {
    // A plausible, mundane failure: the declared heater power is right but a
    // second unlogged heat path means the coolant carries far less. This is
    // exactly the class the gate exists to catch, because nothing about the
    // individual channels looks wrong.
    let spec = spec_with(tight_instruments());
    let run = run_with(EXACT_POWER * 2.0, DatasetPartition::Validation);
    match ingest(&spec, &run).expect_err("a non-closing run must refuse") {
        RigError::EnergyBalanceViolated { detail } => {
            assert!(detail.contains("power in"), "detail: {detail}");
            assert!(detail.contains("heat out"), "detail: {detail}");
            // Every channel and its band, so the reader can see which one to
            // distrust rather than only that something is wrong.
            for channel in ["power", "t-in", "t-out", "flow"] {
                assert!(detail.contains(channel), "missing {channel} in:\n{detail}");
            }
        }
        other => panic!("expected an energy-balance refusal, got {other:?}"),
    }
}

#[test]
fn a_small_error_inside_the_stated_band_is_admitted() {
    // The gate judges against the instruments' OWN declarations, so a
    // discrepancy the metrology already allows for is not an error. A gate
    // that refused this would be unusable on real data.
    let spec = spec_with(tight_instruments());
    let run = run_with(EXACT_POWER + 0.5, DatasetPartition::Validation);
    assert!(
        ingest(&spec, &run).is_ok(),
        "a residual inside the combined band must be admitted"
    );
}

#[test]
fn an_unstated_channel_band_cannot_pass_the_gate() {
    // Level E exists to carry metrology. A channel with no stated uncertainty
    // leaves nothing to adjudicate the closure against, so the run is refused
    // rather than admitted as if it had passed.
    let mut instruments = tight_instruments();
    instruments[0] = instrument("power", ChannelRole::HeaterPower, POWER, None);
    let spec = spec_with(instruments);
    let run = run_with(EXACT_POWER, DatasetPartition::Validation);
    match ingest(&spec, &run).expect_err("an unstated band cannot be adjudicated") {
        RigError::EnergyBalanceUnadjudicated { reason } => {
            assert!(reason.contains("Unstated"), "reason: {reason}");
        }
        other => panic!("expected an unadjudicated refusal, got {other:?}"),
    }
}

#[test]
fn a_zero_temperature_rise_cannot_be_adjudicated() {
    let spec = spec_with(tight_instruments());
    let mut run = run_with(EXACT_POWER, DatasetPartition::Validation);
    run.readings[2].value_si = INLET; // outlet == inlet
    match ingest(&spec, &run).expect_err("no resolvable enthalpy flow") {
        RigError::EnergyBalanceUnadjudicated { reason } => {
            assert!(reason.contains("coolant-rise"), "reason: {reason}");
        }
        other => panic!("expected an unadjudicated refusal, got {other:?}"),
    }
}

#[test]
fn multiplicative_cross_terms_are_enclosed_instead_of_linearized_away() {
    // Nominal H = 1. With 10% independent bounded bands on K, Q and dT,
    // the exact positive product hull is [0.9^3, 1.1^3] = [0.729, 1.331].
    // First-order propagation would use only +/-0.3 and wrongly reject
    // P=1.32. The interval gate must admit it because the supports overlap.
    let instruments = vec![
        instrument("power", ChannelRole::HeaterPower, POWER, Some(0.0)),
        instrument(
            "t-in",
            ChannelRole::InletTemperature,
            TEMPERATURE,
            Some(0.05),
        ),
        instrument(
            "t-out",
            ChannelRole::OutletTemperature,
            TEMPERATURE,
            Some(0.05),
        ),
        instrument("flow", ChannelRole::VolumeFlow, VOLUME_FLOW, Some(0.1)),
    ];
    let spec = RigSpec::new(
        "cross-term-rig",
        CoolantProperties::new(1.0, 1.0, 0.1).expect("coolant"),
        instruments,
    )
    .expect("spec");
    let run = RigRun {
        run_id: "cross-term-run".to_string(),
        acquired_on: "2026-06-15".to_string(),
        partition: DatasetPartition::Validation,
        readings: vec![
            Reading {
                sensor_id: "power".to_string(),
                value_si: 1.32,
            },
            Reading {
                sensor_id: "t-in".to_string(),
                value_si: 10.0,
            },
            Reading {
                sensor_id: "t-out".to_string(),
                value_si: 11.0,
            },
            Reading {
                sensor_id: "flow".to_string(),
                value_si: 1.0,
            },
        ],
    };
    let ingested = ingest(&spec, &run).expect("the exact interval overlap admits");
    match ingested.balance() {
        EnergyBalance::Balanced {
            heat_interval_w,
            residual_interval_w,
            ..
        } => {
            assert!(
                heat_interval_w.lo() <= 0.729,
                "lower endpoint must widen outward: {heat_interval_w:?}"
            );
            assert!(
                heat_interval_w.hi() >= 1.331,
                "upper endpoint must include cross terms: {heat_interval_w:?}"
            );
            assert!(residual_interval_w.contains(0.0));
        }
        other => panic!("expected balanced, got {other:?}"),
    }

    // The refusal direction of the same headline: the exact hull's lower
    // endpoint is 0.9^3 = 0.729 W, while a symmetric first-order band around
    // 1.0 W reaches down to 0.7 W. A power reading in [0.7, 0.729) is
    // provably inconsistent with the declared supports yet survives the
    // linearization, so pin the refusal to keep the low side from silently
    // regressing to a widened lower endpoint.
    let low_run = RigRun {
        run_id: "cross-term-low-run".to_string(),
        acquired_on: "2026-06-15".to_string(),
        partition: DatasetPartition::Validation,
        readings: vec![
            Reading {
                sensor_id: "power".to_string(),
                value_si: 0.71,
            },
            Reading {
                sensor_id: "t-in".to_string(),
                value_si: 10.0,
            },
            Reading {
                sensor_id: "t-out".to_string(),
                value_si: 11.0,
            },
            Reading {
                sensor_id: "flow".to_string(),
                value_si: 1.0,
            },
        ],
    };
    match ingest(&spec, &low_run).expect_err("a sub-hull power reading must refuse") {
        RigError::EnergyBalanceViolated { detail } => {
            assert!(detail.contains("power in"), "detail: {detail}");
        }
        other => panic!("expected an energy-balance refusal, got {other:?}"),
    }
}

#[test]
fn covariance_without_a_coverage_policy_cannot_impersonate_a_support_band() {
    let mut flow = sensor("flow", VOLUME_FLOW, Some(1.0e-5));
    let squared = VOLUME_FLOW
        .checked_plus(VOLUME_FLOW)
        .expect("squared flow dimensions");
    flow.uncertainty = MeasurementUncertainty::CovarianceDiagonal {
        variance: QtyAny::new(1.0e-10, squared),
    };
    let mut instruments = tight_instruments();
    instruments[3] =
        InstrumentSpec::new(flow, ChannelRole::VolumeFlow).expect("valid covariance row");
    let spec = spec_with(instruments);
    match ingest(&spec, &run_with(EXACT_POWER, DatasetPartition::Validation))
        .expect_err("one sigma is not a finite support interval")
    {
        RigError::EnergyBalanceUnadjudicated { reason } => {
            assert!(reason.contains("CovarianceDiagonal"), "reason: {reason}");
            assert!(reason.contains("flow"), "reason: {reason}");
        }
        other => panic!("expected an unadjudicated refusal, got {other:?}"),
    }
}

#[test]
fn physical_direction_must_be_resolved_for_power_flow_and_temperature_rise() {
    let mut power_spans_zero = tight_instruments();
    power_spans_zero[0] = instrument("power", ChannelRole::HeaterPower, POWER, Some(EXACT_POWER));
    let error = ingest(
        &spec_with(power_spans_zero),
        &run_with(EXACT_POWER, DatasetPartition::Validation),
    )
    .expect_err("a power support touching zero is unresolved");
    assert!(matches!(error, RigError::EnergyBalanceUnadjudicated { .. }));

    let mut flow_touches_zero = tight_instruments();
    flow_touches_zero[3] = instrument("flow", ChannelRole::VolumeFlow, VOLUME_FLOW, Some(FLOW));
    let error = ingest(
        &spec_with(flow_touches_zero),
        &run_with(EXACT_POWER, DatasetPartition::Validation),
    )
    .expect_err("a flow support touching zero is unresolved");
    assert!(matches!(error, RigError::EnergyBalanceUnadjudicated { .. }));

    // Two wrong signs must never cancel into a seemingly positive Q*dT.
    let spec = spec_with(tight_instruments());
    let mut double_negative = run_with(EXACT_POWER, DatasetPartition::Validation);
    double_negative.readings[1].value_si = OUTLET;
    double_negative.readings[2].value_si = INLET;
    double_negative.readings[3].value_si = -FLOW;
    match ingest(&spec, &double_negative).expect_err("negative flow is not the declared direction")
    {
        RigError::PhysicalDirectionViolated { field } => {
            assert_eq!(field, "volume-flow interval");
        }
        other => panic!("expected resolved wrong-direction refusal, got {other:?}"),
    }
}

// ---------------------------------------------------------------------------
// Channel matching and calibration.
// ---------------------------------------------------------------------------

#[test]
fn channel_mismatches_refuse_by_sensor_name() {
    let spec = spec_with(tight_instruments());

    let mut unknown = run_with(EXACT_POWER, DatasetPartition::Validation);
    unknown.readings.push(Reading {
        sensor_id: "stray".to_string(),
        value_si: 1.0,
    });
    assert_eq!(
        ingest(&spec, &unknown).expect_err("undeclared channel"),
        RigError::UnknownReading {
            sensor_id: "stray".to_string()
        }
    );

    let mut missing = run_with(EXACT_POWER, DatasetPartition::Validation);
    missing.readings.retain(|r| r.sensor_id != "flow");
    assert_eq!(
        ingest(&spec, &missing).expect_err("missing channel"),
        RigError::MissingReading {
            sensor_id: "flow".to_string()
        }
    );

    let mut duplicated = run_with(EXACT_POWER, DatasetPartition::Validation);
    duplicated.readings.push(Reading {
        sensor_id: "flow".to_string(),
        value_si: FLOW,
    });
    assert_eq!(
        ingest(&spec, &duplicated).expect_err("duplicate channel"),
        RigError::DuplicateReading {
            sensor_id: "flow".to_string()
        }
    );

    let mut non_finite = run_with(EXACT_POWER, DatasetPartition::Validation);
    non_finite.readings[0].value_si = f64::NAN;
    assert_eq!(
        ingest(&spec, &non_finite).expect_err("non-finite reading"),
        RigError::InvalidReading {
            sensor_id: "power".to_string(),
            field: "value_si",
            requirement: "finite",
        }
    );
}

#[test]
fn an_uncalibrated_sensor_refuses_at_specification() {
    let mut instruments = tight_instruments();
    let mut uncalibrated = sensor("power", POWER, Some(0.1));
    uncalibrated.calibration = Availability::Unavailable {
        reason: "certificate expired 2025-11".to_string(),
    };
    instruments[0] =
        InstrumentSpec::new(uncalibrated, ChannelRole::HeaterPower).expect("instrument admits");

    match RigSpec::new("bench-rig-a", coolant(), instruments)
        .expect_err("an uncalibrated channel must refuse")
    {
        RigError::Uncalibrated { sensor_id, reason } => {
            assert_eq!(sensor_id, "power");
            assert!(reason.contains("expired"), "the reason must survive");
        }
        other => panic!("expected an uncalibrated refusal, got {other:?}"),
    }
}

#[test]
fn instrument_admission_reuses_the_canonical_sensor_metrology_boundary() {
    let mut zero_certificate = sensor("power", POWER, Some(0.1));
    if let Availability::Available(calibration) = &mut zero_certificate.calibration {
        calibration.certificate_hash = fs_blake3::ContentHash([0; 32]);
    }
    match InstrumentSpec::new(zero_certificate, ChannelRole::HeaterPower)
        .expect_err("a zero certificate hash must refuse")
    {
        RigError::InvalidSensor { sensor_id, detail } => {
            assert_eq!(sensor_id, "power");
            assert!(detail.contains("nonzero"), "detail: {detail}");
        }
        other => panic!("expected canonical sensor refusal, got {other:?}"),
    }

    let mut bad_date = sensor("power", POWER, Some(0.1));
    if let Availability::Available(calibration) = &mut bad_date.calibration {
        calibration.issued_on = "2026-02-30".to_string();
    }
    assert!(matches!(
        InstrumentSpec::new(bad_date, ChannelRole::HeaterPower),
        Err(RigError::InvalidSensor { .. })
    ));

    let mut negative_band = sensor("power", POWER, Some(0.1));
    negative_band.uncertainty = MeasurementUncertainty::Bounded {
        half_width: QtyAny::new(-0.1, POWER),
    };
    assert!(matches!(
        InstrumentSpec::new(negative_band, ChannelRole::HeaterPower),
        Err(RigError::InvalidSensor { .. })
    ));

    let length = Dims([1, 0, 0, 0, 0, 0]);
    let mut bad_placement = sensor("power", POWER, Some(0.1));
    bad_placement.placement = Availability::Available(SensorPlacement {
        frame: "rig-frame".to_string(),
        coordinates: [QtyAny::new(0.0, length); 3],
        uncertainty: [
            QtyAny::new(0.0, length),
            QtyAny::new(-1.0e-3, length),
            QtyAny::new(0.0, length),
        ],
    });
    assert!(matches!(
        InstrumentSpec::new(bad_placement, ChannelRole::HeaterPower),
        Err(RigError::InvalidSensor { .. })
    ));

    let invalid_id = sensor("Power Sensor", POWER, Some(0.1));
    assert!(matches!(
        InstrumentSpec::new(invalid_id, ChannelRole::HeaterPower),
        Err(RigError::InvalidSensor { .. })
    ));
}

#[test]
fn acquisition_date_is_canonical_and_covered_by_every_calibration() {
    let spec = spec_with(tight_instruments());
    let mut malformed = run_with(EXACT_POWER, DatasetPartition::Validation);
    malformed.acquired_on = "2026-02-30".to_string();
    assert!(matches!(
        ingest(&spec, &malformed),
        Err(RigError::InvalidDate {
            field: "run acquisition date",
            ..
        })
    ));

    let mut expired = run_with(EXACT_POWER, DatasetPartition::Validation);
    expired.acquired_on = "2028-01-01".to_string();
    match ingest(&spec, &expired).expect_err("every certificate expired before this run") {
        RigError::CalibrationOutOfWindow {
            sensor_id,
            acquired_on,
            ..
        } => {
            assert_eq!(sensor_id, "flow", "sensor-id order is deterministic");
            assert_eq!(acquired_on, "2028-01-01");
        }
        other => panic!("expected calibration-window refusal, got {other:?}"),
    }
}

#[test]
fn a_role_with_the_wrong_dimensions_refuses() {
    // A thermocouple declared as the flow channel is a wiring error, and it
    // would otherwise produce a confidently wrong enthalpy.
    let error = InstrumentSpec::new(
        sensor("flow", TEMPERATURE, Some(0.05)),
        ChannelRole::VolumeFlow,
    )
    .expect_err("dimensions must match the role");
    assert!(matches!(error, RigError::RoleDimensions { .. }));
}

#[test]
fn a_rig_missing_or_repeating_a_balance_role_refuses() {
    let mut short = tight_instruments();
    short.retain(|i| i.role() != ChannelRole::VolumeFlow);
    assert_eq!(
        RigSpec::new("rig", coolant(), short).expect_err("no flow channel"),
        RigError::MissingRole {
            role: ChannelRole::VolumeFlow
        }
    );

    let mut doubled = tight_instruments();
    doubled.push(instrument(
        "flow-2",
        ChannelRole::VolumeFlow,
        VOLUME_FLOW,
        Some(1.0e-5),
    ));
    assert_eq!(
        RigSpec::new("rig", coolant(), doubled).expect_err("two flow channels"),
        RigError::DuplicateRole {
            role: ChannelRole::VolumeFlow
        }
    );

    let mut repeated = tight_instruments();
    repeated.push(instrument(
        "flow",
        ChannelRole::Auxiliary,
        VOLUME_FLOW,
        Some(1.0),
    ));
    assert_eq!(
        RigSpec::new("rig", coolant(), repeated).expect_err("duplicate sensor id"),
        RigError::DuplicateSensor {
            sensor_id: "flow".to_string()
        }
    );
}

#[test]
fn auxiliary_channels_are_carried_but_do_not_join_the_closure() {
    let mut instruments = tight_instruments();
    instruments.push(instrument(
        "ambient",
        ChannelRole::Ambient,
        TEMPERATURE,
        Some(0.5),
    ));
    let spec = spec_with(instruments);
    let mut run = run_with(EXACT_POWER, DatasetPartition::Validation);
    run.readings.push(Reading {
        sensor_id: "ambient".to_string(),
        value_si: 293.0,
    });

    let ingested = ingest(&spec, &run).expect("auxiliary channels do not disturb the closure");
    assert_eq!(ingested.calibrated_channels(), 5);
    assert_eq!(
        ingested.readings().expect("validation readings retained")[0],
        Reading {
            sensor_id: "ambient".to_string(),
            value_si: 293.0,
        },
        "retained readings are canonicalized by sensor id"
    );
    match ingested.balance() {
        EnergyBalance::Balanced { residual_w, .. } => assert!(residual_w.abs() < 1e-9),
        other => panic!("expected balanced, got {other:?}"),
    }
}

// ---------------------------------------------------------------------------
// Blind sealing.
// ---------------------------------------------------------------------------

#[test]
fn a_blind_holdout_run_is_gated_then_sealed() {
    // Both properties matter and the ORDER is what delivers them: the gate
    // runs, so a blind run is known physically sane, and the values are then
    // withheld so they cannot leak into ordinary validation.
    let spec = spec_with(tight_instruments());
    let run = run_with(EXACT_POWER, DatasetPartition::BlindHoldout);
    let ingested = ingest(&spec, &run).expect("a sane blind run ingests");

    assert!(ingested.sealed());
    assert_eq!(ingested.partition(), DatasetPartition::BlindHoldout);
    assert_eq!(
        ingested.balance(),
        &EnergyBalance::Sealed,
        "a sealed run must not surface its measured values"
    );
    assert!(
        ingested.readings().is_none(),
        "ordinary access must not reveal blind raw values"
    );
    assert!(
        !format!("{ingested:?}").contains("120.6"),
        "Debug must not leak blind readings"
    );
    let restored = fs_vvreg::rig::IngestedRun::from_retained_bytes(ingested.retained_bytes())
        .expect("explicit retained evidence round-trips");
    assert_eq!(restored, ingested);
    assert!(restored.readings().is_none());

    // And sealing does not weaken the gate: a blind run that does NOT close
    // is still refused rather than hidden.
    let bad = run_with(EXACT_POWER * 3.0, DatasetPartition::BlindHoldout);
    match ingest(&spec, &bad).expect_err("an unbalanced blind run must refuse") {
        RigError::EnergyBalanceViolated { detail } => {
            assert!(detail.contains("REDACTED"), "detail: {detail}");
            for value in [
                format!("{}", EXACT_POWER * 3.0),
                format!("{EXACT_POWER}"),
                format!("{FLOW}"),
                format!("{}", OUTLET - INLET),
            ] {
                assert!(
                    !detail.contains(&value),
                    "blind failure detail leaked measured value {value:?}: {detail}"
                );
            }
        }
        other => panic!("expected redacted energy-balance refusal, got {other:?}"),
    }
}

// ---------------------------------------------------------------------------
// Determinism and identity.
// ---------------------------------------------------------------------------

#[test]
fn ingest_is_deterministic_and_declaration_order_independent() {
    let spec = spec_with(tight_instruments());
    let mut reordered = tight_instruments();
    reordered.reverse();
    let spec_reordered = spec_with(reordered);

    let run = run_with(EXACT_POWER, DatasetPartition::Validation);
    let mut shuffled = run.clone();
    shuffled.readings.reverse();

    let a = ingest(&spec, &run).expect("a");
    let b = ingest(&spec_reordered, &shuffled).expect("b");
    assert_eq!(
        a.identity(),
        b.identity(),
        "order must not move the identity"
    );
    assert_eq!(a, b);
}

#[test]
fn v2_retained_frame_and_identities_match_the_independent_frozen_oracle() {
    // These values were generated from the documented v2 byte grammar with
    // an independent encoder plus the official BLAKE3 derive-key API. They
    // deliberately do not round-trip through fs-vvreg's own decoder, so an
    // encoder+decoder drift cannot silently bless itself.
    let ingested = ingest(
        &spec_with(tight_instruments()),
        &run_with(EXACT_POWER, DatasetPartition::Validation),
    )
    .expect("frozen fixture admits");
    assert_eq!(ingested.retained_bytes().len(), 899);
    assert_eq!(
        ingested.identity(),
        fs_blake3::ContentHash::from_hex(
            "cb4bbbba8140229e784e9d29d9389ab4ac78cf5c2a53bffc1a9e00ca2d8e1898"
        )
        .expect("valid frozen v2 identity")
    );
    assert_eq!(
        ingested.raw_readings_identity(),
        fs_blake3::ContentHash::from_hex(
            "0182eab209aea8fd9cdd37f62f5b3461dd9274806cd0feebef196769b22cd244"
        )
        .expect("valid frozen raw-readings identity")
    );
    assert_eq!(
        &ingested.retained_bytes()[..23],
        &[
            b'F', b'S', b'V', b'V', b'R', b'I', b'G', 0, 2, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 11, 0, 0,
            0,
        ],
        "magic, schema, policy tags, nested sensor schema, and rig-id length are frozen"
    );
}

#[test]
fn the_identity_moves_with_the_measured_content() {
    let spec = spec_with(tight_instruments());
    let baseline =
        ingest(&spec, &run_with(EXACT_POWER, DatasetPartition::Validation)).expect("baseline");
    let nudged = ingest(
        &spec,
        &run_with(EXACT_POWER + 0.01, DatasetPartition::Validation),
    )
    .expect("nudged");
    assert_ne!(
        baseline.identity(),
        nudged.identity(),
        "a 10 mW change must move the identity"
    );

    // Partition is part of what the record claims, so it must move it too.
    let blind = ingest(
        &spec,
        &run_with(EXACT_POWER, DatasetPartition::BlindHoldout),
    )
    .expect("blind");
    assert_ne!(baseline.identity(), blind.identity());
}

#[test]
fn identity_binds_auxiliary_raw_bits_metrology_coolant_and_acquisition_date() {
    let mut instruments = tight_instruments();
    instruments.push(instrument(
        "ambient",
        ChannelRole::Ambient,
        TEMPERATURE,
        Some(0.5),
    ));
    let spec = spec_with(instruments);
    let mut run = run_with(EXACT_POWER, DatasetPartition::Validation);
    run.readings.push(Reading {
        sensor_id: "ambient".to_string(),
        value_si: 293.0,
    });
    let baseline = ingest(&spec, &run).expect("baseline");

    let mut changed_aux = run.clone();
    changed_aux
        .readings
        .iter_mut()
        .find(|reading| reading.sensor_id == "ambient")
        .expect("ambient")
        .value_si = 294.0;
    let changed_aux = ingest(&spec, &changed_aux).expect("aux change");
    assert_ne!(baseline.identity(), changed_aux.identity());
    assert_ne!(
        baseline.raw_readings_identity(),
        changed_aux.raw_readings_identity()
    );

    let mut signed_zero = run.clone();
    signed_zero
        .readings
        .iter_mut()
        .find(|reading| reading.sensor_id == "ambient")
        .expect("ambient")
        .value_si = 0.0;
    let plus_zero = ingest(&spec, &signed_zero).expect("+0");
    signed_zero
        .readings
        .iter_mut()
        .find(|reading| reading.sensor_id == "ambient")
        .expect("ambient")
        .value_si = -0.0;
    let minus_zero = ingest(&spec, &signed_zero).expect("-0");
    assert_ne!(
        plus_zero.identity(),
        minus_zero.identity(),
        "v2 preserves exact finite reading bits"
    );

    let mut later = run.clone();
    later.acquired_on = "2026-06-16".to_string();
    let later = ingest(&spec, &later).expect("later acquisition within every certificate");
    assert_ne!(baseline.identity(), later.identity());
    assert_eq!(
        baseline.raw_readings_identity(),
        later.raw_readings_identity(),
        "the raw-reading artifact is date-independent while the run record is not"
    );

    // Equal rho*cp is not semantic equality: both declarations remain bound.
    let product_power = 2.0 * FLOW * (OUTLET - INLET);
    let spec_a = RigSpec::new(
        "product-rig",
        CoolantProperties::new(1.0, 2.0, 0.001).expect("a"),
        tight_instruments(),
    )
    .expect("spec a");
    let spec_b = RigSpec::new(
        "product-rig",
        CoolantProperties::new(2.0, 1.0, 0.001).expect("b"),
        tight_instruments(),
    )
    .expect("spec b");
    let product_run = run_with(product_power, DatasetPartition::Validation);
    let a = ingest(&spec_a, &product_run).expect("a");
    let b = ingest(&spec_b, &product_run).expect("b");
    assert_ne!(
        a.identity(),
        b.identity(),
        "rho and cp are separately identity-bearing"
    );
}

#[test]
fn every_nested_metrology_family_moves_the_v2_identity() {
    let baseline_spec = spec_with(tight_instruments());
    let run = run_with(EXACT_POWER, DatasetPartition::Validation);
    let baseline = ingest(&baseline_spec, &run).expect("baseline");

    let changed_identity = |mutated: SensorRecord| {
        let mut instruments = tight_instruments();
        instruments[0] =
            InstrumentSpec::new(mutated, ChannelRole::HeaterPower).expect("mutated sensor admits");
        ingest(&spec_with(instruments), &run)
            .expect("mutated run remains physically admissible")
            .identity()
    };

    let mut instrument_id = sensor("power", POWER, Some(0.1));
    instrument_id.instrument_id = Availability::Available("instrument-power-v2".to_string());
    assert_ne!(baseline.identity(), changed_identity(instrument_id));

    let mut raw_channel = sensor("power", POWER, Some(0.1));
    raw_channel.raw_channel = "ch_power_rewired".to_string();
    assert_ne!(baseline.identity(), changed_identity(raw_channel));

    let mut calibration = sensor("power", POWER, Some(0.1));
    if let Availability::Available(record) = &mut calibration.calibration {
        record.certificate_id = "power-certificate-v2".to_string();
    }
    assert_ne!(baseline.identity(), changed_identity(calibration));

    let mut placement = sensor("power", POWER, Some(0.1));
    placement.placement = Availability::Unavailable {
        reason: "synthetic alternate placement declaration".to_string(),
    };
    assert_ne!(baseline.identity(), changed_identity(placement));

    let uncertainty = sensor("power", POWER, Some(0.11));
    assert_ne!(baseline.identity(), changed_identity(uncertainty));
}

#[test]
fn retained_codec_refuses_truncation_policy_tampering_and_trailing_bytes() {
    let ingested = ingest(
        &spec_with(tight_instruments()),
        &run_with(EXACT_POWER, DatasetPartition::Validation),
    )
    .expect("baseline");

    let mut truncated = ingested.retained_bytes().to_vec();
    truncated.pop();
    assert!(matches!(
        fs_vvreg::rig::IngestedRun::from_retained_bytes(&truncated),
        Err(RigError::MalformedRetainedRecord { .. })
    ));

    let mut wrong_policy = ingested.retained_bytes().to_vec();
    wrong_policy[12] ^= 0x7f;
    assert!(matches!(
        fs_vvreg::rig::IngestedRun::from_retained_bytes(&wrong_policy),
        Err(RigError::MalformedRetainedRecord { .. })
    ));

    let mut trailing = ingested.retained_bytes().to_vec();
    trailing.push(0);
    assert!(matches!(
        fs_vvreg::rig::IngestedRun::from_retained_bytes(&trailing),
        Err(RigError::MalformedRetainedRecord { .. })
    ));
}

#[test]
fn admission_refuses_degenerate_declarations() {
    assert!(
        CoolantProperties::new(0.0, CP, 0.0).is_err(),
        "zero density"
    );
    assert!(
        CoolantProperties::new(RHO, -1.0, 0.0).is_err(),
        "negative specific heat"
    );
    assert!(
        CoolantProperties::new(RHO, CP, -0.1).is_err(),
        "negative relative band"
    );
    assert!(
        CoolantProperties::new(RHO, CP, 1.0).is_err(),
        "a 100% relative band loses a strictly positive property lower bound"
    );
    assert!(
        CoolantProperties::new(f64::MAX, 2.0, 0.0).is_err(),
        "overflowing property product"
    );
    assert!(
        CoolantProperties::new(f64::MIN_POSITIVE, f64::MIN_POSITIVE, 0.0).is_err(),
        "underflowing property product"
    );
    assert!(
        RigSpec::new("   ", coolant(), tight_instruments()).is_err(),
        "blank rig id"
    );
    assert!(
        RigSpec::new("rig", coolant(), Vec::new()).is_err(),
        "no instruments"
    );

    let spec = spec_with(tight_instruments());
    let mut blank = run_with(EXACT_POWER, DatasetPartition::Validation);
    blank.run_id = "  ".to_string();
    assert_eq!(
        ingest(&spec, &blank).expect_err("blank run id"),
        RigError::BlankIdentifier { field: "run id" }
    );

    let mut invalid_id = run_with(EXACT_POWER, DatasetPartition::Validation);
    invalid_id.run_id = "Run 1".to_string();
    assert!(matches!(
        ingest(&spec, &invalid_id),
        Err(RigError::InvalidIdentifier {
            field: "run id",
            ..
        })
    ));

    let mut oversized = run_with(EXACT_POWER, DatasetPartition::Validation);
    oversized.readings = (0..=MAX_INSTRUMENTS)
        .map(|index| Reading {
            sensor_id: format!("sensor-{index}"),
            value_si: 1.0,
        })
        .collect();
    assert_eq!(
        ingest(&spec, &oversized).expect_err("count is checked before sorting"),
        RigError::TooManyReadings {
            have: MAX_INSTRUMENTS + 1,
            max: MAX_INSTRUMENTS,
        }
    );

    let oversized_bytes = vec![0_u8; MAX_RETAINED_RIG_BYTES + 1];
    assert_eq!(
        fs_vvreg::rig::IngestedRun::from_retained_bytes(&oversized_bytes)
            .expect_err("retained cap"),
        RigError::RetainedRecordTooLarge {
            have: MAX_RETAINED_RIG_BYTES + 1,
            max: MAX_RETAINED_RIG_BYTES,
        }
    );
}

#[test]
fn unrepresentable_heat_interval_refuses_instead_of_wrapping_or_passing() {
    let instruments = vec![
        instrument("power", ChannelRole::HeaterPower, POWER, Some(0.0)),
        instrument(
            "t-in",
            ChannelRole::InletTemperature,
            TEMPERATURE,
            Some(0.0),
        ),
        instrument(
            "t-out",
            ChannelRole::OutletTemperature,
            TEMPERATURE,
            Some(0.0),
        ),
        instrument("flow", ChannelRole::VolumeFlow, VOLUME_FLOW, Some(0.0)),
    ];
    let spec = RigSpec::new(
        "overflow-rig",
        CoolantProperties::new(1.0e200, 1.0, 0.0).expect("finite property product"),
        instruments,
    )
    .expect("spec");
    let run = RigRun {
        run_id: "overflow-run".to_string(),
        acquired_on: "2026-06-15".to_string(),
        partition: DatasetPartition::Validation,
        readings: vec![
            Reading {
                sensor_id: "power".to_string(),
                value_si: 1.0,
            },
            Reading {
                sensor_id: "t-in".to_string(),
                value_si: 1.0,
            },
            Reading {
                sensor_id: "t-out".to_string(),
                value_si: 2.0,
            },
            Reading {
                sensor_id: "flow".to_string(),
                value_si: 1.0e200,
            },
        ],
    };
    assert!(matches!(
        ingest(&spec, &run),
        Err(RigError::InvalidScalar { .. })
    ));
}
