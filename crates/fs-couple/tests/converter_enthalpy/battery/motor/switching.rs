//! PWM electrical work heats the original finite material owners exactly once.
use super::*;
use fs_phs::circuit::converter::switching::battery_motor::{
    BatteryMotorSwitchCommand, PreparedBatteryMotorSwitchCycle,
};

fn cycle_controls() -> SwitchCycleControls {
    SwitchCycleControls {
        absolute_energy_tolerance_j: 1e-9,
        relative_energy_tolerance: 1e-11,
        ..SwitchCycleControls::default()
    }
}
fn thermal_controls() -> CircuitEnthalpySolve {
    let mut controls = drive_controls();
    controls.maximum_circuit_trials = 20_000;
    controls
}
fn command(torques: [f64; 2], load: Option<RegulatedPowerLoad>) -> BatteryMotorSwitchCommand {
    BatteryMotorSwitchCommand {
        shaft_torque_nm: torques,
        auxiliary_load: load,
    }
}

#[test]
fn switched_ten_body_calorimetry_retains_original_native_work_and_external_shaft_loads() {
    let fixture = DriveFixture::new([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0], 0.0);
    let context = point(300.0);
    let cells = fixture.cells(&context);
    let key = key();
    let channels = drive_paths(&key);
    let network = BatteryMotorThermalNetwork {
        cells: &cells,
        links: &[],
        channels: &channels,
        bearing_cell: DRIVE_IDS[9],
    };
    let body = fixture.body();
    let (native, before, torques, load, _) = steady_drive(1.0, 0.1);
    let drive = network.bind(&body, &native, before, || false).unwrap();
    let intervals = centered_pwm(0.02, native.electrical().modulation()).unwrap();
    let commands = vec![command(torques, load); intervals.len()];
    let mut original = PreparedBatteryMotorSwitchCycle::new_controlled(
        native,
        &intervals,
        4,
        electrical_controls(),
        cycle_controls(),
        || false,
    )
    .unwrap();
    let mut expected = [0.0; 9];
    let reference = original
        .step_interval_commands_controlled(&before, &commands, &mut expected, || false)
        .unwrap();
    let next = network
        .step_switch_candidate(
            &body,
            0.02,
            &drive,
            &intervals,
            &commands,
            4,
            electrical_controls(),
            cycle_controls(),
            thermal_controls(),
            &[0.0; 10],
            || false,
        )
        .unwrap();
    assert_eq!(*next.drive().state(), expected);
    let step = &next.report().electrical;
    assert_eq!(step.heat_joules(), reference.heat_joules());
    for i in 0..10 {
        near(
            fixture.cp[i] * (next.thermal().temperature_k(DRIVE_IDS[i]).unwrap() - 300.0),
            reference.heat_joules()[i],
        );
    }
    near(step.load_work_j, 0.02 * load.unwrap().requested_power_w());
    near(
        next.report().deposited_joule_j,
        reference.heat_joules().iter().sum(),
    );
    near(next.report().non_electrical_heat_j, 0.0);
    near(next.report().total_energy_residual_j, 0.0);
    near(
        next.thermal().energy_j() - body.energy_j(),
        next.report().deposited_joule_j,
    );
    near(
        step.native.supplied,
        step.shaft_work_j.iter().sum::<f64>() - step.load_work_j,
    );
    let cost = intervals.len() * 4;
    assert_eq!(step.integration_steps, cost);
    assert_eq!(next.report().circuit_trials % cost, 0);
    assert!(next.report().circuit_trials >= 2 * cost);
}

#[test]
fn switched_endpoint_resistance_changes_actual_motor_trajectory_and_replays_at_that_material_state()
{
    let (_, _, _, _, heat_w) = steady_drive(1.0, 0.1);
    let cp = heat_w.map(|power| 0.02 * power);
    let mut currents = Vec::new();
    for alpha in [0.0, 0.0039] {
        let fixture = DriveFixture::new(cp, alpha);
        let context = point(300.0);
        let cells = fixture.cells(&context);
        let key = key();
        let channels = drive_paths(&key);
        let network = BatteryMotorThermalNetwork {
            cells: &cells,
            links: &[],
            channels: &channels,
            bearing_cell: DRIVE_IDS[9],
        };
        let body = fixture.body();
        let (native, before, torques, load, _) = steady_drive(1.0, 0.1);
        let drive = network.bind(&body, &native, before, || false).unwrap();
        let intervals = centered_pwm(0.02, native.electrical().modulation()).unwrap();
        let commands = vec![command(torques, load); intervals.len()];
        let result = network
            .step_switch_candidate(
                &body,
                0.02,
                &drive,
                &intervals,
                &commands,
                4,
                electrical_controls(),
                cycle_controls(),
                thermal_controls(),
                &[0.0; 10],
                || false,
            )
            .unwrap();
        for i in 0..10 {
            let temperature = result.thermal().temperature_k(DRIVE_IDS[i]).unwrap();
            near(
                cp[i] * (temperature - 300.0),
                result.report().electrical.heat_joules()[i],
            );
            if i < 9 {
                near(
                    result.drive().channels()[i].resistance_ohm(),
                    DRIVE_R[i] * (1.0 + alpha * (temperature - 300.0)),
                );
                assert_eq!(
                    drive.channels()[i].material().unwrap().card_identity(),
                    result.drive().channels()[i]
                        .material()
                        .unwrap()
                        .card_identity()
                );
            }
        }
        // Re-run the ORIGINAL native cycle from the original state using the
        // material coefficients actually accepted at the caloric endpoint.
        let accepted = result.drive();
        let original =
            AveragedBuckBoost::new(accepted.converter_parameters(), accepted.modulation())
                .unwrap()
                .with_battery_source(
                    BatterySource::new(
                        accepted.battery_curve().clone(),
                        accepted.battery_parameters(),
                    )
                    .unwrap(),
                    accepted.lead_resistance_ohm(),
                )
                .unwrap()
                .with_dc_motor(DcMotor::new(accepted.motor_parameters()).unwrap())
                .unwrap();
        let mut replay = PreparedBatteryMotorSwitchCycle::new_controlled(
            original,
            &intervals,
            4,
            electrical_controls(),
            cycle_controls(),
            || false,
        )
        .unwrap();
        let mut replayed = [0.0; 9];
        let step = replay
            .step_interval_commands_controlled(&before, &commands, &mut replayed, || false)
            .unwrap();
        assert_eq!(replayed, *accepted.state());
        assert_eq!(step.heat_joules(), result.report().electrical.heat_joules());
        near(result.report().total_energy_residual_j, 0.0);
        currents.push(accepted.state()[6]);
    }
    assert!(
        (currents[1] - currents[0]).abs() > 1e-9,
        "R(T) must reach actual motor flux"
    );
}

#[test]
fn switched_clock_budget_and_late_cancellation_preserve_the_accepted_history() {
    let fixture = DriveFixture::new([1.0; 10], 0.0);
    let context = point(300.0);
    let cells = fixture.cells(&context);
    let key = key();
    let channels = drive_paths(&key);
    let network = BatteryMotorThermalNetwork {
        cells: &cells,
        links: &[],
        channels: &channels,
        bearing_cell: DRIVE_IDS[9],
    };
    let body = fixture
        .body()
        .advance_to(
            0.1,
            StepConfig {
                max_dt_s: 0.1,
                energy_tolerance_j: 1e-8,
            },
            || Ok(()),
        )
        .unwrap()
        .0;
    let (native, before, torques, load, _) = steady_drive(1.0, 0.1);
    let drive = network.bind(&body, &native, before, || false).unwrap();
    let intervals = centered_pwm(0.02, native.electrical().modulation()).unwrap();
    let commands = vec![command(torques, load); intervals.len()];
    let target = body.time_s() + 0.02;
    assert_ne!((target - body.time_s()).to_bits(), 0.02f64.to_bits());
    let mut polls = 0;
    let accepted = network
        .step_switch_candidate(
            &body,
            target,
            &drive,
            &intervals,
            &commands,
            4,
            electrical_controls(),
            cycle_controls(),
            thermal_controls(),
            &[0.0; 10],
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    near(accepted.thermal().time_s(), target);
    assert!(
        network
            .step_switch_candidate(
                &body,
                target.next_up(),
                &drive,
                &intervals,
                &commands,
                4,
                electrical_controls(),
                cycle_controls(),
                thermal_controls(),
                &[0.0; 10],
                || false,
            )
            .is_err(),
        "a different represented endpoint cannot borrow a clock allowance"
    );
    let mut short = thermal_controls();
    short.maximum_circuit_trials = intervals.len() * 4 - 1;
    assert!(matches!(
        network.step_switch_candidate(
            &body,
            target,
            &drive,
            &intervals,
            &commands,
            4,
            electrical_controls(),
            cycle_controls(),
            short,
            &[0.0; 10],
            || false,
        ),
        Err(CircuitEnthalpyError::CircuitBudget)
    ));
    let mut count = 0;
    assert!(matches!(
        network.step_switch_candidate(
            &body,
            target,
            &drive,
            &intervals,
            &commands,
            4,
            electrical_controls(),
            cycle_controls(),
            thermal_controls(),
            &[0.0; 10],
            || {
                count += 1;
                count >= polls
            },
        ),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    assert_eq!(drive.state(), &before);
    assert_eq!(drive.time_s(), 0.1);
    assert_eq!(body.time_s(), 0.1);
    for id in DRIVE_IDS {
        near(body.temperature_k(id).unwrap(), 300.0);
    }
}
