//! Level-E rig ingest: channel matching, calibration, the energy-balance
//! gate, and blind sealing (bead `frankensim-extreal-program-f85xj.4.5`).
//!
//! Every fixture here is SYNTHETIC. That is the point: the bead's caveat is
//! that the software half must not wait on hardware, so the pipeline is
//! exercised against runs whose physics we control exactly and whose failure
//! modes we can inject deliberately.

use fs_qty::{Dims, QtyAny};
use fs_vvreg::corpus::{
    Availability, CalibrationRecord, DatasetPartition, MeasurementUncertainty, SensorRecord,
};
use fs_vvreg::rig::{
    ChannelRole, CoolantProperties, EnergyBalance, InstrumentSpec, Reading, RigError, RigRun,
    RigSpec, ingest,
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

    match ingested.balance() {
        EnergyBalance::Balanced {
            power_in_w,
            heat_out_w,
            residual_w,
            allowed_w,
        } => {
            assert!((power_in_w - EXACT_POWER).abs() < 1e-12);
            assert!((heat_out_w - EXACT_POWER).abs() < 1e-9);
            assert!(residual_w.abs() < 1e-9, "an exact fixture must close");
            assert!(*allowed_w > 0.0, "a stated band must be positive");
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
    let mut run = run_with(0.0, DatasetPartition::Validation);
    run.readings[2].value_si = INLET; // outlet == inlet
    match ingest(&spec, &run).expect_err("no resolvable enthalpy flow") {
        RigError::EnergyBalanceUnadjudicated { reason } => {
            assert!(reason.contains("zero"), "reason: {reason}");
        }
        other => panic!("expected an unadjudicated refusal, got {other:?}"),
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
    assert!(matches!(
        ingest(&spec, &non_finite).expect_err("non-finite reading"),
        RigError::InvalidScalar { .. }
    ));
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

    // And sealing does not weaken the gate: a blind run that does NOT close
    // is still refused rather than hidden.
    let bad = run_with(EXACT_POWER * 3.0, DatasetPartition::BlindHoldout);
    assert!(
        matches!(
            ingest(&spec, &bad),
            Err(RigError::EnergyBalanceViolated { .. })
        ),
        "sealing must not become a way to smuggle in unbalanced data"
    );
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
}
