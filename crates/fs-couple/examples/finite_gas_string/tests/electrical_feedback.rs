//! Endpoint source work reaches the actual material/gas consumer and its CSV.
use super::*;
use fs_couple::acoustic_realize::finite_gas::electrical::feedback::FiniteGasElectricalFeedbackSolve;
use fs_material::conductor::OhmicDrive;

const FEED: &str = "thermochemical 300 95000 mole 1 light 1";
const ROW: &str = "electrical_feedback 0.000005 0.000000000001 64";

fn scenario(thermoelastic: bool, flowing: bool) -> input::Spec {
    let mut spec = super::flow::uniform_scenario(thermoelastic);
    spec.samples = 16;
    spec.release_sample = 2;
    spec.electrical_feedback = Some(FiniteGasElectricalFeedbackSolve {
        maximum_joule_work_j: 5e-6,
        source_tolerance_j: 1e-12,
        maximum_iterations: 64,
    });
    if !flowing {
        spec.nozzle = None;
    }
    spec
}

#[test]
fn source_feedback_input_is_explicit_bounded_and_not_ignored_by_other_paths() {
    let text = include_str!("../feedback-flow.string");
    let spec = input::parse(text).unwrap();
    let controls = spec.electrical_feedback.unwrap();
    assert_eq!(controls.maximum_joule_work_j, 5e-6);
    assert_eq!(controls.source_tolerance_j, 1e-12);
    assert_eq!(controls.maximum_iterations, 64);
    assert!(spec.nozzle.is_some());
    assert!(
        input::parse(include_str!("../heated.string"))
            .unwrap()
            .electrical_feedback
            .is_none()
    );
    assert!(input::parse(&text.replace("drive current 0.1", "drive voltage -0.25")).is_ok());
    assert!(input::parse(&text.replace("drive current 0.1", "drive current 0")).is_ok());
    for row in [
        "electrical_feedback 0 0.000000000001 64",
        "electrical_feedback -1 0.000000000001 64",
        "electrical_feedback NaN 0.000000000001 64",
        "electrical_feedback inf 0.000000000001 64",
        "electrical_feedback 0.000005 0 64",
        "electrical_feedback 0.000005 0.000006 64",
        "electrical_feedback 0.000005 0.000000000001 0",
        "electrical_feedback 0.000005 0.000000000001 129",
        "electrical_feedback 0.000005 0.000000000001",
    ] {
        assert!(
            input::parse(&text.replace(ROW, row)).is_err(),
            "accepted {row}"
        );
    }
    assert!(input::parse(&text.replace("drive current 0.1", "drive none")).is_err());
    let bonded = include_str!("../bonded.string").replace("\nheat ", &format!("\n{ROW}\nheat "));
    assert!(input::parse(&bonded).is_err());
    assert!(input::parse(&text.replace(ROW, &format!("{ROW}\n{ROW}"))).is_err());
}

#[test]
fn endpoint_current_and_voltage_work_close_the_actual_consumer_energy() {
    for (thermoelastic, flowing, drive) in [
        (false, false, OhmicDrive::Current(0.1)),
        (false, true, OhmicDrive::Voltage(-0.25)),
        (true, false, OhmicDrive::Voltage(0.25)),
        (true, true, OhmicDrive::Current(-0.1)),
    ] {
        scoped(|cx, _| {
            let mut spec = scenario(thermoelastic, flowing);
            spec.drive = Some(drive);
            let card = source_law(1.0, thermoelastic.then_some(2e-5));
            let (gas, reservoir) = super::flow::inventories(&spec, FEED).unwrap();
            let a = simulate(cx, &spec, &card, gas.clone(), reservoir.clone()).unwrap();
            let b = simulate(cx, &spec, &card, gas.clone(), reservoir.clone()).unwrap();
            assert_eq!(a.records, b.records);
            assert_eq!(a.pressure_pa, b.pressure_pa);
            assert_eq!(a.electrical_feedback, b.electrical_feedback);
            let feedback = a.electrical_feedback.as_ref().unwrap();
            assert_eq!(feedback.len(), spec.samples);
            let area = std::f64::consts::PI * spec.radius_m.powi(2);
            let (mut previous_work, mut previous_heat) = (0.0, 0.0);
            let mut endpoint_differs = false;
            for (row, source) in a.records.iter().zip(feedback) {
                // Independent source-card R(T), followed by the prescribed
                // ideal current or voltage law at the accepted temperature.
                let resistivity = 8e-7 + (2e-6 - 8e-7) * (row[4] - 280.0) / 120.0;
                let resistance = resistivity * spec.length_m / area;
                let (current, power) = match drive {
                    OhmicDrive::Current(i) => (i, i * i * resistance),
                    OhmicDrive::Voltage(v) => (v / resistance, v * v / resistance),
                };
                let work = power / f64::from(spec.sample_rate_hz);
                assert!((row[20] - resistance).abs() < 1e-12);
                assert!((row[18] - current).abs() < 1e-13);
                assert!((source[0] - previous_work - work).abs() < 1e-14);
                assert!((row[10] - source[0] - source[1]).abs() < 1e-14);
                let deposited = row[10] - previous_heat;
                assert!((deposited - work).abs() <= 1.01e-12);
                let frozen = row[17].powi(2) * row[19] / f64::from(spec.sample_rate_hz);
                endpoint_differs |= (work - frozen).abs() > 1e-12;
                assert!(source[2] >= 2.0 && source[2] <= 66.0 && source[2].fract() == 0.0);
                // These fixtures have no authored damping. The total balance
                // must debit actual source work, not the nearby deposited heat.
                let defect = row[14] - row[13] - source[0] - row[11] - row[12];
                assert!((row[15] - defect).abs() < 1e-14);
                assert!(row[15].abs() <= row[16]);
                previous_work = source[0];
                previous_heat = row[10];
            }
            assert!(endpoint_differs);
            assert!(a.records.last().unwrap()[9].abs() > 1e-12);
            if flowing {
                let flow = a.gas_flow.as_ref().unwrap();
                assert_eq!(flow.records, b.gas_flow.as_ref().unwrap().records);
                assert!(flow.records.iter().all(|row| row[7] > 0.0));
                assert!(flow.records.last().unwrap()[14] > flow.records[0][14]);
            }
            let mut explicit = spec.clone();
            explicit.electrical_feedback = None;
            let old = simulate(cx, &explicit, &card, gas, reservoir).unwrap();
            assert!(old.electrical_feedback.is_none());
            assert_ne!(old.records.last().unwrap()[4], a.records.last().unwrap()[4]);
        });
    }
}

#[test]
fn zero_drive_retains_the_existing_thermal_and_flow_trajectory() {
    scoped(|cx, _| {
        let mut implicit = scenario(false, true);
        implicit.drive = Some(OhmicDrive::Voltage(0.0));
        let card = source(1.0);
        let (gas, reservoir) = super::flow::inventories(&implicit, FEED).unwrap();
        let a = simulate(cx, &implicit, &card, gas.clone(), reservoir.clone()).unwrap();
        let mut explicit = implicit.clone();
        explicit.electrical_feedback = None;
        let b = simulate(cx, &explicit, &card, gas, reservoir).unwrap();
        assert_eq!(a.pressure_pa, b.pressure_pa);
        assert_eq!(
            a.gas_flow.as_ref().unwrap().records,
            b.gas_flow.as_ref().unwrap().records
        );
        for (a, b) in a.records.iter().zip(&b.records) {
            for index in 0..27 {
                if index != 16 {
                    assert_eq!(a[index], b[index]);
                }
            }
        }
        assert!(
            a.electrical_feedback
                .as_ref()
                .unwrap()
                .iter()
                .all(|row| *row == [0.0, 0.0, 1.0])
        );
    });
}

#[test]
fn feedback_bracket_source_exit_and_cancellation_refuse_the_complete_run() {
    scoped(|cx, gate| {
        let spec = scenario(false, true);
        let card = source(1.0);
        let (gas, reservoir) = super::flow::inventories(&spec, FEED).unwrap();
        let mut short = spec.clone();
        short
            .electrical_feedback
            .as_mut()
            .unwrap()
            .maximum_joule_work_j = 1e-9;
        assert!(simulate(cx, &short, &card, gas.clone(), reservoir.clone()).is_err());
        let mut outside_source = spec.clone();
        outside_source
            .electrical_feedback
            .as_mut()
            .unwrap()
            .maximum_joule_work_j = 1.0;
        assert!(simulate(cx, &outside_source, &card, gas.clone(), reservoir.clone()).is_err());
        gate.request();
        assert!(simulate(cx, &spec, &card, gas, reservoir).is_err());
    });
}
