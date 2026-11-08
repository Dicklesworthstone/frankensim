//! Real outer speed, machine-owned torque allocation and inner current control
//! around finite battery/inverter/magnetic/two-shaft storage. Numerical fixture
//! only: no universal stability, overcurrent or experimental motor claim.
use fs_control::speed::{
    ShaftPowerLimits, SpeedCommand, SpeedController, SpeedControllerParameters,
    SpeedControllerState,
};
use fs_control::synchronous::{
    CurrentControllerParameters, CurrentControllerState, SynchronousCurrentController,
    SynchronousMeasurement,
};
use fs_phs::circuit::battery::{BatteryOcvCurve, BatteryParameters, BatterySource};
use fs_phs::circuit::synchronous::inverter::battery::{
    BatteryInverterControls, PreparedBatterySynchronousInverter,
};
use fs_phs::circuit::synchronous::inverter::{AveragedSynchronousInverter, InverterParameters};
use fs_phs::circuit::synchronous::{SynchronousMotor, SynchronousMotorParameters};
use fs_phs::{PreparedStepError, SolveBudget};

const TS: f64 = 2e-4;
const CURRENT_LIMIT: f64 = 8.0;

fn machine() -> SynchronousMotorParameters {
    SynchronousMotorParameters {
        d_inductance_h: 0.01,
        q_inductance_h: 0.014,
        winding_resistance_ohm: 0.2,
        pm_flux_linkage_wb: 0.08,
        pole_pairs: 3,
        rotor_inertia_kg_m2: 0.025,
        stator_inertia_kg_m2: 0.07,
        bearing_drag_nm_s_per_rad: 0.002,
    }
}

fn capabilities() -> ShaftPowerLimits {
    ShaftPowerLimits {
        positive_torque_limit_nm: 3.0,
        negative_torque_limit_nm: 3.0,
        motoring_power_w: 40.0,
        regeneration_power_w: 15.0,
    }
}

fn fixture() -> (
    PreparedBatterySynchronousInverter,
    SpeedController,
    SynchronousCurrentController,
    [f64; 9],
) {
    let p = machine();
    let motor = SynchronousMotor::new(p).unwrap();
    let motor_state = motor.state([0.0; 2], 0.1, [0.0; 2]).unwrap();
    let battery = BatterySource::new(
        BatteryOcvCurve::new(20.0, &[(0.0, 20.0), (0.5, 24.0), (1.0, 26.0)]).unwrap(),
        BatteryParameters {
            series_resistance_ohm: 0.05,
            polarization_resistance_ohm: 0.08,
            polarization_capacitance_f: 5.0,
            bus_capacitance_f: 0.1,
        },
    )
    .unwrap();
    let inverter = AveragedSynchronousInverter::new(
        motor,
        InverterParameters {
            dc_capacitance_f: 0.15,
            high_on_resistance_ohm: [0.01, 0.012, 0.014],
            low_on_resistance_ohm: [0.014, 0.01, 0.012],
        },
    )
    .unwrap()
    .with_battery_source(battery, 0.03)
    .unwrap();
    let state = inverter.state(0.75, 0.0, 25.0, 25.0, &motor_state).unwrap();
    let plant = inverter
        .prepare_controlled(
            BatteryInverterControls {
                solve: SolveBudget {
                    relative: 2e-11,
                    absolute: 2e-12,
                    max_iterations: 40,
                },
                backtracking_halvings: 8,
                absolute_energy_tolerance_j: 2e-9,
                relative_energy_tolerance: 2e-10,
            },
            || false,
        )
        .unwrap();
    let speed = SpeedController::new(SpeedControllerParameters {
        sample_period_s: TS,
        proportional_gain_nm_s_per_rad: 0.4,
        integral_gain_nm_per_rad: 4.0,
        anti_windup_time_s: 0.05,
    })
    .unwrap();
    let current = SynchronousCurrentController::new(CurrentControllerParameters {
        sample_period_s: TS,
        pole_pairs: p.pole_pairs,
        proportional_gain_ohm: [2.0; 2],
        integral_gain_ohm_per_s: [200.0; 2],
        anti_windup_time_s: [0.002; 2],
        winding_resistance_ohm: p.winding_resistance_ohm,
        current_limit_a: CURRENT_LIMIT,
    })
    .unwrap();
    (plant, speed, current, state)
}

fn measurement(plant: &PreparedBatterySynchronousInverter, x: &[f64; 9]) -> SynchronousMeasurement {
    let point = plant.machine_operating_point(x).unwrap();
    SynchronousMeasurement {
        current_alpha_beta_a: point.current_alpha_beta_a,
        flux_alpha_beta_wb: [x[4], x[5]],
        relative_mechanical_angle_rad: x[6],
        shaft_speed_rad_s: point.shaft_speed_rad_s,
        dc_bus_voltage_v: x[3] / 0.15,
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}

// Independent authored SI energy: integrate the two OCV segments and use the
// original capacitor, flux and momentum coordinates, with no native H callback.
fn physical_energy(x: &[f64; 9]) -> f64 {
    let p = machine();
    let chemical = if x[0] <= 10.0 {
        20.0 * x[0] + 0.2 * x[0] * x[0]
    } else {
        let q = x[0] - 10.0;
        220.0 + 24.0 * q + 0.1 * q * q
    };
    let (s, c) = (f64::from(p.pole_pairs) * x[6]).sin_cos();
    let flux_d = c * x[4] + s * x[5] - p.pm_flux_linkage_wb;
    let flux_q = -s * x[4] + c * x[5];
    chemical
        + 0.5
            * (x[1] * x[1] / 5.0
                + x[2] * x[2] / 0.1
                + x[3] * x[3] / 0.15
                + flux_d * flux_d / p.d_inductance_h
                + flux_q * flux_q / p.q_inductance_h
                + x[7] * x[7] / p.rotor_inertia_kg_m2
                + x[8] * x[8] / p.stator_inertia_kg_m2)
}

#[test]
fn finite_battery_speed_loop_rejects_load_and_regeneratively_brakes_with_both_shafts() {
    let (mut plant, speed, current, mut state) = fixture();
    let initial = state;
    let mut accepted_speed = SpeedControllerState::default();
    let mut accepted_current = CurrentControllerState::default();
    let mut total_heat = 0.0;
    let mut total_shaft_work = 0.0;
    let mut current_limited = 0;
    let mut regeneration_limited = 0;
    let mut phase_charge_change = [0.0; 3];
    let mut phase_bridge_work = [0.0; 3];
    let mut tail_error = [0.0_f64; 3];
    for (phase, (target, torques, samples)) in [
        (15.0, [0.0, 0.0], 4000),
        (15.0, [-0.2, 0.2], 2000),
        (5.0, [0.0, 0.0], 3000),
    ]
    .into_iter()
    .enumerate()
    {
        let initial_charge = state[0];
        for sample in 0..samples {
            let measured = measurement(&plant, &state);
            let demand = speed
                .demand_controlled(
                    &accepted_speed,
                    measured.shaft_speed_rad_s,
                    SpeedCommand {
                        relative_speed_rad_s: target,
                        feedforward_torque_nm: 0.0,
                    },
                    capabilities(),
                    || false,
                )
                .unwrap();
            let allocation = plant
                .current_reference_for_torque_controlled(
                    demand.requested_torque_nm,
                    0.0,
                    CURRENT_LIMIT,
                    || false,
                )
                .unwrap();
            current_limited += usize::from(allocation.torque_limited);
            regeneration_limited +=
                usize::from(demand.power_limited && demand.requested_shaft_power_w < 0.0);
            let outer = demand
                .finish_controlled(allocation.achievable_torque_nm, || false)
                .unwrap();
            let inner = current
                .propose_controlled(
                    &accepted_current,
                    &measured,
                    allocation.current_dq_a,
                    || false,
                )
                .unwrap();
            assert!(!inner.current_limited);
            assert!(outer.achievable_shaft_power_w <= capabilities().motoring_power_w + 1e-12);
            assert!(outer.achievable_shaft_power_w >= -capabilities().regeneration_power_w - 1e-12);
            let mut next = [0.0; 9];
            let step = plant
                .step_controlled(
                    &state,
                    inner.duty_high_fraction,
                    torques,
                    TS,
                    &mut next,
                    || false,
                )
                .unwrap();
            let heat = step.heat_joules().iter().sum::<f64>();
            let shaft_work = step.inverter.shaft_work_j.iter().sum::<f64>();
            close(
                physical_energy(&next) - physical_energy(&state) + heat,
                shaft_work,
                3e-9,
            );
            close(step.native.dissipated, heat, 2e-10);
            close(step.native.supplied, shaft_work, 2e-11);
            close(step.balance_residual_j(), 0.0, 3e-9);
            close(
                step.battery.terminal_work_j + step.inverter.source_work_j + step.lead_heat_j,
                0.0,
                2e-11,
            );
            let shaft_impulse = step.inverter.electromagnetic_angular_impulse_nm_s
                + step.inverter.bearing_angular_impulse_nm_s;
            close(next[7] - state[7], shaft_impulse + TS * torques[0], 5e-11);
            close(next[8] - state[8], -shaft_impulse + TS * torques[1], 5e-11);
            assert!(step.heat_joules().iter().all(|heat| *heat >= 0.0));
            total_heat += heat;
            total_shaft_work += shaft_work;
            phase_bridge_work[phase] += step.inverter.bridge_work_j;
            // Exactly one publication point for both PI states and all nine
            // physical states, after the original native solve has accepted.
            accepted_speed = outer.next_state;
            accepted_current = inner.next_state;
            state = next;
            if sample >= samples - 100 {
                let actual = state[7] / machine().rotor_inertia_kg_m2
                    - state[8] / machine().stator_inertia_kg_m2;
                tail_error[phase] = tail_error[phase].max((actual - target).abs());
            }
        }
        phase_charge_change[phase] = state[0] - initial_charge;
    }
    assert!(
        tail_error.into_iter().all(|error| error < 0.5),
        "speed errors [rad/s]: {tail_error:?}"
    );
    assert!(current_limited > 0 && regeneration_limited > 0);
    assert!(
        phase_charge_change[0] < -0.01 && phase_charge_change[1] < 0.0,
        "finite chemical charge changes: {phase_charge_change:?}"
    );
    assert!(
        phase_charge_change[2] > 0.001 && phase_bridge_work[2] < -0.05,
        "braking must return energy: charge {phase_charge_change:?}, bridge work {phase_bridge_work:?}"
    );
    assert!(state[7] > 0.0 && state[8] < 0.0);
    close(state[7] + state[8], initial[7] + initial[8], 3e-9);
    assert!(total_heat > 0.1 && total_shaft_work < 0.0);
    close(
        physical_energy(&state) - physical_energy(&initial) + total_heat,
        total_shaft_work,
        3e-7,
    );
}

#[test]
fn rejected_and_cancelled_native_intervals_do_not_accept_either_control_state() {
    let (mut plant, speed, current, state) = fixture();
    let outer_state = SpeedControllerState::default();
    let inner_state = CurrentControllerState::default();
    let measured = measurement(&plant, &state);
    let demand = speed
        .demand_controlled(
            &outer_state,
            measured.shaft_speed_rad_s,
            SpeedCommand {
                relative_speed_rad_s: 15.0,
                feedforward_torque_nm: 0.0,
            },
            capabilities(),
            || false,
        )
        .unwrap();
    let allocation = plant
        .current_reference_for_torque_controlled(
            demand.requested_torque_nm,
            0.0,
            CURRENT_LIMIT,
            || false,
        )
        .unwrap();
    let outer = demand
        .finish_controlled(allocation.achievable_torque_nm, || false)
        .unwrap();
    let inner = current
        .propose_controlled(&inner_state, &measured, allocation.current_dq_a, || false)
        .unwrap();
    assert_ne!(outer.next_state, outer_state);
    assert_ne!(inner.next_state, inner_state);
    let held = plant.high_fraction();
    let mut next = [71.0; 9];
    let mut unsupported = state;
    unsupported[0] = -0.1;
    assert!(
        plant
            .step_controlled(
                &unsupported,
                inner.duty_high_fraction,
                [0.0; 2],
                TS,
                &mut next,
                || false
            )
            .is_err()
    );
    assert_eq!(next, [71.0; 9]);
    let mut polls = 0;
    assert!(matches!(
        plant.step_controlled(
            &state,
            inner.duty_high_fraction,
            [0.0; 2],
            TS,
            &mut next,
            || {
                polls += 1;
                polls >= 3
            }
        ),
        Err(PreparedStepError::Cancelled)
    ));
    assert_eq!(next, [71.0; 9]);
    assert_eq!(plant.high_fraction(), held);
    assert_eq!(outer_state, SpeedControllerState::default());
    assert_eq!(inner_state, CurrentControllerState::default());
    let first = plant
        .step_controlled(
            &state,
            inner.duty_high_fraction,
            [0.0; 2],
            TS,
            &mut next,
            || false,
        )
        .unwrap();
    let accepted = next;
    let (mut fresh, _, _, fresh_state) = fixture();
    let replay = fresh
        .step_controlled(
            &fresh_state,
            inner.duty_high_fraction,
            [0.0; 2],
            TS,
            &mut next,
            || false,
        )
        .unwrap();
    assert_eq!(next, accepted);
    assert_eq!(first.native, replay.native);
    assert_eq!(first.battery, replay.battery);
    assert_eq!(first.inverter, replay.inverter);
    assert_eq!(first.heat_joules(), replay.heat_joules());
}
