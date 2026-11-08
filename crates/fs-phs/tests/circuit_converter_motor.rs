//! Finite chemical energy, reciprocal armature power and both physical shafts.
use fs_phs::circuit::battery::{BatteryOcvCurve, BatteryParameters, BatterySource};
use fs_phs::circuit::converter::battery_fed::motor::{BatteryMotorDrive, BatteryMotorStep};
use fs_phs::circuit::converter::power_load::{
    PowerLoadControls, PowerLoadMode, RegulatedPowerLoad,
};
use fs_phs::circuit::converter::{AveragedBuckBoost, BridgeModulation, BuckBoostParameters};
use fs_phs::circuit::motor::{DcMotor, DcMotorParameters};
use fs_phs::{PreparedStepError, SolveBudget};

fn motor_parameters() -> DcMotorParameters {
    DcMotorParameters {
        inductance_h: 0.04,
        winding_resistance_ohm: 0.3,
        rotor_inertia_kg_m2: 0.2,
        stator_inertia_kg_m2: 0.3,
        coupling_nm_per_a: 0.7,
        bearing_drag_nm_s_per_rad: 0.005,
    }
}

fn drive(capacity_c: f64, current_a: f64, speeds: [f64; 2]) -> (BatteryMotorDrive, [f64; 9]) {
    let battery = BatterySource::new(
        BatteryOcvCurve::new(capacity_c, &[(0.0, 10.0), (0.5, 12.0), (1.0, 13.0)]).unwrap(),
        BatteryParameters {
            series_resistance_ohm: 0.08,
            polarization_resistance_ohm: 0.1,
            polarization_capacitance_f: 20.0,
            bus_capacitance_f: 0.3,
        },
    )
    .unwrap();
    let electrical = AveragedBuckBoost::new(
        BuckBoostParameters {
            inductance_h: 0.03,
            bus_capacitance_f: [0.2, 0.25],
            winding_resistance_ohm: 0.03,
            input_on_resistance_ohm: [0.01, 0.02],
            output_on_resistance_ohm: [0.01, 0.02],
        },
        BridgeModulation {
            input_high_fraction: 0.5,
            output_high_fraction: 1.0,
        },
    )
    .unwrap()
    .with_battery_source(battery, 0.04)
    .unwrap();
    let electrical_state = electrical
        .state(0.75, 0.0, 12.5, 0.0, [12.5, 6.25])
        .unwrap();
    let motor = DcMotor::new(motor_parameters()).unwrap();
    let motor_state = motor.state(current_a, speeds).unwrap();
    let model = electrical.with_dc_motor(motor).unwrap();
    let initial = model.state(&electrical_state, &motor_state).unwrap();
    (model, initial)
}

fn controls() -> PowerLoadControls {
    PowerLoadControls {
        solve: SolveBudget {
            relative: 1e-13,
            absolute: 1e-13,
            max_iterations: 40,
        },
        backtracking_halvings: 8,
        absolute_energy_tolerance_j: 2e-10,
        relative_energy_tolerance: 2e-11,
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}

// Reconstruct the winding/shaft work and physical loss channels from the
// ORIGINAL coordinates. Neither native dissipated energy nor a balance
// residual is used as an input to these independently calculated quantities.
fn check_interval(
    initial: &[f64; 9],
    next: &[f64; 9],
    torques: [f64; 2],
    dt: f64,
    step: &BatteryMotorStep,
) {
    let p = motor_parameters();
    let current = f64::midpoint(initial[6], next[6]) / p.inductance_h;
    let speeds = [
        f64::midpoint(initial[7], next[7]) / p.rotor_inertia_kg_m2,
        f64::midpoint(initial[8], next[8]) / p.stator_inertia_kg_m2,
    ];
    let relative = speeds[0] - speeds[1];
    let output_v = f64::midpoint(initial[5], next[5]) / 0.25;
    let magnetic = (next[6] - initial[6]) * current;
    let kinetic = [
        (next[7] - initial[7]) * speeds[0],
        (next[8] - initial[8]) * speeds[1],
    ];
    let copper = dt * p.winding_resistance_ohm * current * current;
    let bearing = dt * p.bearing_drag_nm_s_per_rad * relative * relative;
    let electrical_work = dt * output_v * current;
    let shaft_work = [dt * torques[0] * speeds[0], dt * torques[1] * speeds[1]];
    close(step.magnetic_energy_change_j, magnetic, 2e-12);
    for i in 0..2 {
        close(step.kinetic_energy_change_j[i], kinetic[i], 2e-12);
        close(step.shaft_work_j[i], shaft_work[i], 2e-12);
        close(step.discrete_shaft_speed_rad_s[i], speeds[i], 2e-12);
    }
    close(step.motor_electrical_work_j, electrical_work, 2e-12);
    close(step.electrical.external_work_j, -electrical_work, 2e-12);
    close(step.copper_heat_j, copper, 2e-12);
    close(step.bearing_heat_j, bearing, 2e-12);
    close(step.discrete_motor_current_a, current, 2e-12);
    close(
        step.discrete_back_emf_v,
        p.coupling_nm_per_a * relative,
        2e-12,
    );
    close(
        step.discrete_electromagnetic_torque_nm,
        p.coupling_nm_per_a * current,
        2e-12,
    );
    close(
        step.discrete_bearing_torque_nm,
        -p.bearing_drag_nm_s_per_rad * relative,
        2e-12,
    );
    close(
        magnetic + kinetic.iter().sum::<f64>() + copper + bearing,
        electrical_work + shaft_work.iter().sum::<f64>(),
        3e-10,
    );
    close(
        (next[7] - initial[7]) + (next[8] - initial[8]),
        dt * (torques[0] + torques[1]),
        3e-11,
    );
    close(step.angular_momentum_residual_nms, 0.0, 3e-11);
    let heat = step.heat_joules();
    assert!(heat.into_iter().all(|value| value >= 0.0));
    close(heat[8], copper, 2e-12);
    close(heat[9], bearing, 2e-12);
    close(step.native.dissipated, heat.iter().sum(), 3e-11);
    close(
        step.native.supplied,
        shaft_work.iter().sum::<f64>() - step.electrical.load_work_j,
        3e-11,
    );
    close(step.electrical.balance_residual_j(), 0.0, 3e-10);
    close(step.balance_residual_j(), 0.0, 3e-10);
}

#[test]
fn finite_battery_accelerates_rotor_and_finite_stator_through_the_converter() {
    let (model, initial) = drive(100.0, 0.0, [0.0, 0.0]);
    assert_eq!(model.system().state_dim(), 9);
    assert_eq!(model.system().port_dim(), 4);
    let initial_energy = model.system().hamiltonian(&initial);
    let dt = 0.002;
    let mut prepared = model.prepare_controlled(dt, controls(), || false).unwrap();
    let mut state = initial;
    let mut heat = [0.0; 10];
    let mut armature_work = 0.0;
    for _ in 0..500 {
        let mut next = [0.0; 9];
        let step = prepared
            .step_controlled(&state, &[0.0; 2], None, &mut next, || false)
            .unwrap();
        check_interval(&state, &next, [0.0; 2], dt, &step);
        assert_eq!(step.native.supplied, 0.0);
        for (total, increment) in heat.iter_mut().zip(step.heat_joules()) {
            *total += increment;
        }
        armature_work += step.motor_electrical_work_j;
        state = next;
    }
    assert!(
        state[0] < initial[0] - 0.1,
        "chemical charge must actually supply the drive"
    );
    assert!(
        state[7] > 0.1 && state[8] < -0.1,
        "both finite inertias must respond"
    );
    close(state[7] + state[8], 0.0, 3e-11);
    assert!(armature_work > 1.0 && heat[8] > 0.0 && heat[9] > 0.0);
    assert!(heat[5] > 0.0 && heat[6] > 0.0 && heat[7] > 0.0);
    close(
        prepared.system().hamiltonian(&state) - initial_energy + heat.iter().sum::<f64>(),
        0.0,
        3e-9,
    );
}

#[test]
fn finite_shaft_kinetic_energy_regenerates_into_chemical_storage_on_the_same_path() {
    let (model, initial) = drive(100.0, 0.0, [14.4, -9.6]);
    let initial_energy = model.system().hamiltonian(&initial);
    let dt = 0.002;
    let mut prepared = model.prepare_controlled(dt, controls(), || false).unwrap();
    let mut state = initial;
    let mut heat = 0.0;
    let mut motor_work = 0.0;
    let mut chemical_change = 0.0;
    let mut charging_steps = 0;
    for _ in 0..500 {
        let mut next = [0.0; 9];
        let step = prepared
            .step_controlled(&state, &[0.0; 2], None, &mut next, || false)
            .unwrap();
        check_interval(&state, &next, [0.0; 2], dt, &step);
        if step.electrical.battery.discrete_cell_current_a < 0.0
            && step.electrical.lead_current_a < 0.0
            && step.discrete_motor_current_a < 0.0
        {
            assert!(step.motor_electrical_work_j < 0.0);
            assert!(step.electrical.battery.terminal_work_j > 0.0);
            charging_steps += 1;
        }
        assert_eq!(step.native.supplied, 0.0);
        heat += step.heat_joules().iter().sum::<f64>();
        motor_work += step.motor_electrical_work_j;
        chemical_change += step.electrical.battery.chemical_energy_change_j;
        state = next;
    }
    assert!(charging_steps > 100);
    assert!(state[0] > initial[0] + 0.05 && chemical_change > 0.5);
    assert!(state[7] < initial[7] && state[8] > initial[8]);
    assert!(motor_work < -1.0 && heat > 0.0);
    close(state[7] + state[8], initial[7] + initial[8], 3e-11);
    close(
        prepared.system().hamiltonian(&state) - initial_energy + heat,
        0.0,
        3e-9,
    );
}

#[test]
fn torque_and_auxiliary_demand_change_without_resetting_any_drive_state() {
    let (model, initial) = drive(100.0, 1.0, [4.0, -8.0 / 3.0]);
    let initial_energy = model.system().hamiltonian(&initial);
    let dt = 0.002;
    let mut prepared = model.prepare_controlled(dt, controls(), || false).unwrap();
    let mut state = initial;
    let mut heat = 0.0;
    let mut shaft_work = 0.0;
    let mut load_work = 0.0;
    let phases = [
        (
            [0.0, 0.0],
            RegulatedPowerLoad::try_new(0.2, 2.0, 0.2, 100.0).unwrap(),
            PowerLoadMode::ConstantPower,
        ),
        (
            [-0.1, 0.1],
            RegulatedPowerLoad::try_new(20.0, 0.3, 0.2, 100.0).unwrap(),
            PowerLoadMode::CurrentLimited,
        ),
        (
            [-0.2, 0.05],
            RegulatedPowerLoad::try_new(1.0, 1.0, 20.0, 100.0).unwrap(),
            PowerLoadMode::VoltageFoldback,
        ),
    ];
    for (torques, load, mode) in phases {
        for _ in 0..100 {
            let mut next = [0.0; 9];
            let step = prepared
                .step_controlled(&state, &torques, Some(load), &mut next, || false)
                .unwrap();
            check_interval(&state, &next, torques, dt, &step);
            assert_eq!(step.electrical.load.mode, mode);
            assert!(step.electrical.load_work_j > 0.0);
            if mode != PowerLoadMode::ConstantPower {
                assert!(step.electrical.curtailed_load_energy_j > 0.0);
            }
            shaft_work += step.shaft_work_j.iter().sum::<f64>();
            load_work += step.electrical.load_work_j;
            heat += step.heat_joules().iter().sum::<f64>();
            state = next;
        }
    }
    assert!(shaft_work < 0.0 && load_work > 0.1 && heat > 0.0);
    assert!(state[0] < initial[0]);
    close(state[7] + state[8], initial[7] + initial[8] - 0.03, 3e-11);
    close(
        prepared.system().hamiltonian(&state) - initial_energy + heat + load_work,
        shaft_work,
        3e-9,
    );
}

#[test]
fn unsupported_load_empty_source_and_cancellation_publish_no_partial_drive_state() {
    let (model, initial) = drive(100.0, 0.0, [0.0, 0.0]);
    let mut prepared = model
        .prepare_controlled(0.002, controls(), || false)
        .unwrap();
    let mut next = [91.0; 9];
    let unsupported = RegulatedPowerLoad::try_new(1.0, 1.0, 0.5, 1.0).unwrap();
    assert!(
        prepared
            .step_controlled(&initial, &[0.0; 2], Some(unsupported), &mut next, || false)
            .is_err()
    );
    assert_eq!(next, [91.0; 9]);
    let mut polls = 0;
    assert!(matches!(
        prepared.step_controlled(&initial, &[0.0; 2], None, &mut next, || {
            polls += 1;
            polls >= 3
        }),
        Err(PreparedStepError::Cancelled)
    ));
    assert_eq!(next, [91.0; 9]);
    let first = prepared
        .step_controlled(&initial, &[0.0; 2], None, &mut next, || false)
        .unwrap();
    let accepted = next;
    let retry = prepared
        .step_controlled(&initial, &[0.0; 2], None, &mut next, || false)
        .unwrap();
    assert_eq!(next, accepted);
    assert_eq!(first.native, retry.native);
    assert_eq!(first.electrical, retry.electrical);
    assert_eq!(first.heat_joules(), retry.heat_joules());

    // A supported depleted chemical state with a loaded, sagging bus must not
    // be silently clamped when its discharge would leave the authored OCV law.
    let (small, _) = drive(0.02, 2.0, [0.0, 0.0]);
    let electrical = small
        .electrical()
        .state(0.0, 0.1, 9.7, 2.0, [9.66, 4.72])
        .unwrap();
    let motor = DcMotor::new(motor_parameters()).unwrap();
    let depleted = small
        .state(&electrical, &motor.state(2.0, [0.0; 2]).unwrap())
        .unwrap();
    let mut prepared = small
        .prepare_controlled(0.002, controls(), || false)
        .unwrap();
    next = [91.0; 9];
    assert!(
        prepared
            .step_controlled(&depleted, &[0.0; 2], None, &mut next, || false)
            .is_err()
    );
    assert_eq!(next, [91.0; 9]);
}
