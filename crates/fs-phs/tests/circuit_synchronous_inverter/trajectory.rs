//! Independent finite-bus dq-current trajectory, including physical abc loss.
use fs_phs::SolveBudget;
use fs_phs::circuit::synchronous::inverter::{AveragedSynchronousInverter, InverterParameters};
use fs_phs::circuit::synchronous::{SynchronousMotor, SynchronousMotorParameters};

const INTERVAL: f64 = 0.02;
const INITIAL: [f64; 6] = [4.0, 0.15, -0.8, 0.1, 2.0, -0.5];
const DUTIES: [[f64; 3]; 4] = [
    [0.5, 0.75, 0.25],
    [0.5, 1.0, 0.0],
    [0.35, 1.0, 0.05],
    [0.2, 0.85, 0.45],
];
const SOURCES: [[f64; 3]; 4] = [
    [-0.15, 0.12, -0.02],
    [0.25, 0.08, -0.03],
    [0.1, -0.15, 0.04],
    [-0.05, -0.12, 0.02],
];

fn machine_parameters() -> SynchronousMotorParameters {
    SynchronousMotorParameters {
        d_inductance_h: 0.02,
        q_inductance_h: 0.035,
        winding_resistance_ohm: 0.12,
        pm_flux_linkage_wb: 0.3,
        pole_pairs: 2,
        rotor_inertia_kg_m2: 0.15,
        stator_inertia_kg_m2: 0.3,
        bearing_drag_nm_s_per_rad: 0.001,
    }
}

fn bridge_parameters() -> InverterParameters {
    InverterParameters {
        dc_capacitance_f: 0.15,
        high_on_resistance_ohm: [0.04, 0.06, 0.03],
        low_on_resistance_ohm: [0.025, 0.035, 0.05],
    }
}

fn torque(id: f64, iq: f64) -> f64 {
    let p = machine_parameters();
    f64::from(p.pole_pairs)
        * (p.pm_flux_linkage_wb * iq + (p.d_inductance_h - p.q_inductance_h) * id * iq)
}

// Independent coordinates [Vdc,id,iq,theta,omega_r,omega_s, Wsource,Wr,Ws,
// Qcopper,Qbearing,Qchannel, electromagnetic impulse, ideal-bridge DC work].
// No native storage, effort, discrete-gradient or phase-transform callback is
// used to construct this continuous physical reference.
fn rhs(y: [f64; 14], duty: [f64; 3], sources: [f64; 3]) -> [f64; 14] {
    let p = machine_parameters();
    let bridge = bridge_parameters();
    let (s, c) = (f64::from(p.pole_pairs) * y[3]).sin_cos();
    let alpha = c * y[1] - s * y[2];
    let beta = s * y[1] + c * y[2];
    let abc = [
        (2.0_f64 / 3.0).sqrt() * alpha,
        -alpha / 6.0_f64.sqrt() + beta / 2.0_f64.sqrt(),
        -alpha / 6.0_f64.sqrt() - beta / 2.0_f64.sqrt(),
    ];
    let channel: [f64; 3] = core::array::from_fn(|k| {
        duty[k] * bridge.high_on_resistance_ohm[k]
            + (1.0 - duty[k]) * bridge.low_on_resistance_ohm[k]
    });
    let phase_voltage: [f64; 3] = core::array::from_fn(|k| duty[k] * y[0] - channel[k] * abc[k]);
    // Floating star point: the balanced Clarke rows eliminate common mode;
    // unequal channel drops remain in the physical three-wire winding voltage.
    let va =
        (2.0_f64 / 3.0).sqrt() * (phase_voltage[0] - 0.5 * (phase_voltage[1] + phase_voltage[2]));
    let vb = (phase_voltage[1] - phase_voltage[2]) / 2.0_f64.sqrt();
    let vd = c * va + s * vb;
    let vq = -s * va + c * vb;
    let relative = y[4] - y[5];
    let electrical_speed = f64::from(p.pole_pairs) * relative;
    let electromagnetic = torque(y[1], y[2]);
    let bearing = p.bearing_drag_nm_s_per_rad * relative;
    let dc_current: f64 = duty.iter().zip(abc).map(|(d, i)| d * i).sum();
    let channel_loss: f64 = channel.iter().zip(abc).map(|(r, i)| r * i * i).sum();
    [
        (sources[0] - dc_current) / bridge.dc_capacitance_f,
        (vd - p.winding_resistance_ohm * y[1] + electrical_speed * p.q_inductance_h * y[2])
            / p.d_inductance_h,
        (vq - p.winding_resistance_ohm * y[2]
            - electrical_speed * (p.d_inductance_h * y[1] + p.pm_flux_linkage_wb))
            / p.q_inductance_h,
        relative,
        (sources[1] + electromagnetic - bearing) / p.rotor_inertia_kg_m2,
        (sources[2] - electromagnetic + bearing) / p.stator_inertia_kg_m2,
        sources[0] * y[0],
        sources[1] * y[4],
        sources[2] * y[5],
        p.winding_resistance_ohm * (y[1] * y[1] + y[2] * y[2]),
        bearing * relative,
        channel_loss,
        electromagnetic,
        y[0] * dc_current,
    ]
}

fn reference() -> [[f64; 14]; 4] {
    let mut y = [0.0; 14];
    y[..6].copy_from_slice(&INITIAL);
    let mut boundaries = [[0.0; 14]; 4];
    let dt = INTERVAL / 1024.0;
    for segment in 0..4 {
        // Every discontinuity is a common boundary of reference and native
        // grids. This is a declared held-duty/load trajectory, not a controller.
        for _ in 0..1024 {
            let k1 = rhs(y, DUTIES[segment], SOURCES[segment]);
            let k2 = rhs(
                core::array::from_fn(|i| y[i] + 0.5 * dt * k1[i]),
                DUTIES[segment],
                SOURCES[segment],
            );
            let k3 = rhs(
                core::array::from_fn(|i| y[i] + 0.5 * dt * k2[i]),
                DUTIES[segment],
                SOURCES[segment],
            );
            let k4 = rhs(
                core::array::from_fn(|i| y[i] + dt * k3[i]),
                DUTIES[segment],
                SOURCES[segment],
            );
            y = core::array::from_fn(|i| {
                y[i] + dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0
            });
        }
        boundaries[segment] = y;
    }
    boundaries
}

fn physical_state(x: [f64; 6]) -> [f64; 6] {
    let p = machine_parameters();
    let (s, c) = (f64::from(p.pole_pairs) * x[3]).sin_cos();
    [
        x[0] / bridge_parameters().dc_capacitance_f,
        (c * x[1] + s * x[2] - p.pm_flux_linkage_wb) / p.d_inductance_h,
        (-s * x[1] + c * x[2]) / p.q_inductance_h,
        x[3],
        x[4] / p.rotor_inertia_kg_m2,
        x[5] / p.stator_inertia_kg_m2,
    ]
}

fn energy(x: &[f64]) -> f64 {
    let p = machine_parameters();
    0.5 * (bridge_parameters().dc_capacitance_f * x[0] * x[0]
        + p.d_inductance_h * x[1] * x[1]
        + p.q_inductance_h * x[2] * x[2]
        + p.rotor_inertia_kg_m2 * x[4] * x[4]
        + p.stator_inertia_kg_m2 * x[5] * x[5])
}

#[test]
fn changing_duties_and_loads_refine_in_current_torque_heat_and_signed_work() {
    let reference = reference();
    let final_reference = reference[3];
    let supplied: f64 = final_reference[6..9].iter().sum();
    let dissipated: f64 = final_reference[9..12].iter().sum();
    assert!((energy(&final_reference) - energy(&INITIAL) + dissipated - supplied).abs() < 2e-12);
    let mut previous: Option<[f64; 4]> = None;
    for steps_per_segment in [4, 8, 16] {
        let p = machine_parameters();
        let machine = SynchronousMotor::new(p).unwrap();
        let motor_state = machine
            .state(
                [INITIAL[1], INITIAL[2]],
                INITIAL[3],
                [INITIAL[4], INITIAL[5]],
            )
            .unwrap();
        let inverter = AveragedSynchronousInverter::new(machine, bridge_parameters()).unwrap();
        let mut state = inverter.state(INITIAL[0], &motor_state).unwrap();
        let mut prepared = inverter
            .prepare_controlled(
                SolveBudget {
                    relative: 1e-12,
                    absolute: 1e-13,
                    max_iterations: 30,
                },
                8,
                || false,
            )
            .unwrap();
        let mut totals = [0.0; 8];
        let mut errors = [0.0; 4]; // state, signed work, losses, torque
        let dt = INTERVAL / f64::from(steps_per_segment);
        for segment in 0..4 {
            let before = totals;
            for _ in 0..steps_per_segment {
                let (mut next, mut ports) = ([0.0; 6], [0.0; 3]);
                let receipt = prepared
                    .step_controlled(
                        &state,
                        DUTIES[segment],
                        &SOURCES[segment],
                        dt,
                        &mut next,
                        &mut ports,
                        || false,
                    )
                    .unwrap();
                totals[0] += receipt.source_work_j;
                totals[1] += receipt.shaft_work_j[0];
                totals[2] += receipt.shaft_work_j[1];
                totals[3] += receipt.copper_heat_j;
                totals[4] += receipt.bearing_heat_j;
                totals[5] += receipt.channel_heat_j.iter().flatten().sum::<f64>();
                totals[6] += receipt.electromagnetic_angular_impulse_nm_s;
                totals[7] += receipt.bridge_work_j;
                state = next;
            }
            // The first two held vectors regenerate into the finite capacitor;
            // the last two consume DC work. Source and mechanical work also
            // change sign, so a passive-only or prescribed-voltage substitute
            // cannot satisfy this trajectory.
            assert_eq!(totals[7] < before[7], segment < 2);
            assert_eq!(totals[0] > before[0], segment == 1 || segment == 2);
            assert_eq!(totals[1] > before[1], segment < 2);
            let actual = physical_state(state);
            let expected = reference[segment];
            assert!(actual[0] > 0.0);
            errors[0] += actual
                .iter()
                .zip(expected)
                .map(|(a, b)| (a - b).abs())
                .sum::<f64>();
            errors[1] += (0..3)
                .map(|i| (totals[i] - expected[6 + i]).abs())
                .sum::<f64>();
            errors[1] += (totals[7] - expected[13]).abs();
            errors[2] += (3..6)
                .map(|i| (totals[i] - expected[6 + i]).abs())
                .sum::<f64>();
            let endpoint_torque = prepared
                .machine_operating_point(&state)
                .unwrap()
                .electromagnetic_torque_nm;
            errors[3] += (totals[6] - expected[12]).abs()
                + INTERVAL * (endpoint_torque - torque(expected[1], expected[2])).abs();
        }
        for heat in &totals[3..6] {
            assert!(*heat > 0.0);
        }
        let final_state = physical_state(state);
        let balance = energy(&final_state) - energy(&INITIAL) + totals[3..6].iter().sum::<f64>()
            - totals[..3].iter().sum::<f64>();
        assert!(
            balance.abs() < 2e-10,
            "original physical work balance {balance:e}"
        );
        let external_impulse = INTERVAL * SOURCES.iter().map(|u| u[1] + u[2]).sum::<f64>();
        let initial_momentum =
            p.rotor_inertia_kg_m2 * INITIAL[4] + p.stator_inertia_kg_m2 * INITIAL[5];
        assert!((state[4] + state[5] - initial_momentum - external_impulse).abs() < 2e-11);
        assert!((final_state[4] - INITIAL[4]).abs() > 0.01);
        assert!((final_state[5] - INITIAL[5]).abs() > 0.01);
        if let Some(old) = previous {
            for (quantity, (coarse, fine)) in old.iter().zip(errors).enumerate() {
                let ratio = coarse / fine;
                assert!(
                    ratio > 3.6 && ratio < 4.4,
                    "quantity {quantity}, errors {coarse:e}/{fine:e}, refinement {ratio}"
                );
            }
        }
        previous = Some(errors);
    }
}
