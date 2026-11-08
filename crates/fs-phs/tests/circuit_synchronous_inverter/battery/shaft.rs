//! Independent torsional oscillator and complete dq/circuit/load trajectory.
use super::{LEAD, battery_parameters, bridge_parameters, close, controls, machine, plateau};
use fs_phs::circuit::synchronous::{
    SynchronousMotor, SynchronousMotorParameters,
    inverter::{
        AveragedSynchronousInverter,
        battery::shaft::{ElasticShaftParameters, PreparedElasticShaftDrive},
    },
};

fn shaft() -> ElasticShaftParameters {
    ElasticShaftParameters {
        load_inertia_kg_m2: 0.045,
        stiffness_nm_per_rad: 12.0,
        damping_nm_s_per_rad: 0.03,
    }
}

// Physical coordinates: chemical charge, polarization/battery/DC voltages,
// dq currents, motor relative angle, rotor/stator speeds, twist, load speed.
const INITIAL: [f64; 11] = [60.0, 0.1, 23.5, 23.2, -0.3, 1.5, 0.07, 5.0, -0.5, 0.03, 1.0];

fn prepare(
    p: SynchronousMotorParameters,
    s: ElasticShaftParameters,
    y: [f64; 11],
) -> (PreparedElasticShaftDrive, [f64; 11]) {
    let motor = SynchronousMotor::new(p).unwrap();
    let mechanical = motor.state([y[4], y[5]], y[6], [y[7], y[8]]).unwrap();
    let drive = AveragedSynchronousInverter::new(motor, bridge_parameters())
        .unwrap()
        .with_battery_source(plateau(), LEAD)
        .unwrap();
    let initial = drive
        .state(y[0] / 100.0, y[1], y[2], y[3], &mechanical)
        .unwrap();
    let model = drive.with_elastic_shaft(s).unwrap();
    let initial = model.state(&initial, y[9], y[10]).unwrap();
    (
        model.prepare_controlled(controls(), || false).unwrap(),
        initial,
    )
}

// Independent closed-form relative motion; absolute speeds follow the actual
// conserved rotor+load angular momentum. No pHS matrix or storage is inspected.
fn oscillator(
    t: f64,
    p: SynchronousMotorParameters,
    s: ElasticShaftParameters,
    z0: f64,
    v0: f64,
    momentum: f64,
) -> [f64; 3] {
    let mobility = 1.0 / p.rotor_inertia_kg_m2 + 1.0 / s.load_inertia_kg_m2;
    let alpha = 0.5 * s.damping_nm_s_per_rad * mobility;
    let beta = (s.stiffness_nm_per_rad * mobility - alpha * alpha).sqrt();
    let b = (v0 + alpha * z0) / beta;
    let (sn, cs) = (beta * t).sin_cos();
    let envelope = (-alpha * t).exp();
    let z = envelope * (z0 * cs + b * sn);
    let velocity = envelope * (beta * (-z0 * sn + b * cs) - alpha * (z0 * cs + b * sn));
    let total_inertia = p.rotor_inertia_kg_m2 + s.load_inertia_kg_m2;
    [
        z,
        (momentum + s.load_inertia_kg_m2 * velocity) / total_inertia,
        (momentum - p.rotor_inertia_kg_m2 * velocity) / total_inertia,
    ]
}

#[test]
fn finite_shaft_matches_two_inertia_cayley_motion_and_refines_to_damped_oscillator() {
    let p = SynchronousMotorParameters {
        pm_flux_linkage_wb: 0.0,
        winding_resistance_ohm: 0.0,
        bearing_drag_nm_s_per_rad: 0.0,
        ..machine()
    };
    let initial = [60.0, 0.0, 24.0, 24.0, 0.0, 0.0, 0.2, 3.0, 0.0, 0.06, -1.0];
    for damping in [0.0, 0.06] {
        let s = ElasticShaftParameters {
            stiffness_nm_per_rad: 20.0,
            damping_nm_s_per_rad: damping,
            ..shaft()
        };
        let momentum = p.rotor_inertia_kg_m2 * initial[7] + s.load_inertia_kg_m2 * initial[10];
        let mobility = 1.0 / p.rotor_inertia_kg_m2 + 1.0 / s.load_inertia_kg_m2;
        let exact = oscillator(0.2, p, s, initial[9], initial[7] - initial[10], momentum);
        let mut errors = Vec::new();
        for count in [16, 32, 64] {
            let (mut model, mut state) = prepare(p, s, initial);
            let initial_native = state;
            let dt = 0.2 / f64::from(count);
            let (mut z, mut velocity, mut heat) = (initial[9], initial[7] - initial[10], 0.0);
            for _ in 0..count {
                let spring = 0.25 * dt * dt * s.stiffness_nm_per_rad * mobility;
                let drag = 0.5 * dt * damping * mobility;
                let new_velocity = ((1.0 - drag - spring) * velocity
                    - dt * s.stiffness_nm_per_rad * mobility * z)
                    / (1.0 + drag + spring);
                let new_z = z + 0.5 * dt * (velocity + new_velocity);
                let mean_velocity = 0.5 * (velocity + new_velocity);
                let mut next = [0.0; 11];
                let step = model
                    .step_controlled(&state, [0.0; 3], [0.0; 3], dt, &mut next, || false)
                    .unwrap();
                close(next[9], new_z, 3e-10);
                close(
                    next[7] / p.rotor_inertia_kg_m2 - next[10] / s.load_inertia_kg_m2,
                    new_velocity,
                    3e-9,
                );
                close(
                    step.shaft_torque_nm,
                    s.stiffness_nm_per_rad * 0.5 * (z + new_z) + damping * mean_velocity,
                    3e-9,
                );
                close(
                    step.shaft_heat_j,
                    dt * damping * mean_velocity * mean_velocity,
                    3e-10,
                );
                close(next[7] + next[10], momentum, 3e-11);
                close(next[8], 0.0, 3e-11);
                close(step.shaft_balance_residual_j(), 0.0, 3e-9);
                close(
                    step.load_kinetic_energy_change_j,
                    step.load_shaft_work_j,
                    3e-9,
                );
                close(
                    step.inverter.shaft_work_j[0],
                    -step.rotor_shaft_work_j,
                    3e-11,
                );
                assert_eq!(step.external_shaft_work_j, [0.0; 3]);
                assert!(step.heat_joules()[..11].iter().all(|h| h.abs() < 1e-20));
                for i in 0..6 {
                    close(next[i], initial_native[i], 3e-10);
                }
                heat += step.shaft_heat_j;
                (z, velocity, state) = (new_z, new_velocity, next);
            }
            let actual = [
                state[9],
                state[7] / p.rotor_inertia_kg_m2,
                state[10] / s.load_inertia_kg_m2,
            ];
            errors.push(
                (0..3)
                    .map(|i| (actual[i] - exact[i]).abs())
                    .fold(0.0, f64::max),
            );
            let mechanical_energy = |x: &[f64; 11]| {
                0.5 * (s.stiffness_nm_per_rad * x[9] * x[9]
                    + x[7] * x[7] / p.rotor_inertia_kg_m2
                    + x[10] * x[10] / s.load_inertia_kg_m2)
            };
            close(
                mechanical_energy(&state) - mechanical_energy(&initial_native) + heat,
                0.0,
                3e-9,
            );
            if damping > 0.0 {
                assert!(heat > 0.0);
            } else {
                assert_eq!(heat, 0.0);
            }
        }
        assert!(errors[0] / errors[1] > 3.4, "{errors:?}");
        assert!(errors[1] / errors[2] > 3.7, "{errors:?}");
    }
}

const SEGMENT: f64 = 0.012;
const DUTIES: [[f64; 3]; 4] = [
    [0.5, 0.8, 0.2],
    [0.3, 0.75, 0.4],
    [0.6, 0.25, 0.8],
    [0.5, 0.65, 0.35],
];
const TORQUES: [[f64; 3]; 4] = [
    [0.04, -0.02, -0.1],
    [0.0, 0.015, -0.35],
    [-0.08, 0.0, 0.15],
    [0.03, -0.01, -0.05],
];

// Independent continuous Kirchhoff + dq machine + elastic shaft equations.
// Eleven PHYSICAL coordinates, three external work integrals, seven distinct
// heat integrals. Bridge heat sums six channels; no native effort/storage is used.
fn rhs(y: [f64; 21], duty: [f64; 3], torque: [f64; 3]) -> [f64; 21] {
    let (p, s, b, bridge) = (
        machine(),
        shaft(),
        battery_parameters(),
        bridge_parameters(),
    );
    let (sn, cs) = (f64::from(p.pole_pairs) * y[6]).sin_cos();
    let (ia, ib) = (cs * y[4] - sn * y[5], sn * y[4] + cs * y[5]);
    let abc = [
        (2.0_f64 / 3.0).sqrt() * ia,
        -ia / 6.0_f64.sqrt() + ib / 2.0_f64.sqrt(),
        -ia / 6.0_f64.sqrt() - ib / 2.0_f64.sqrt(),
    ];
    let resistance: [f64; 3] = std::array::from_fn(|j| {
        duty[j] * bridge.high_on_resistance_ohm[j]
            + (1.0 - duty[j]) * bridge.low_on_resistance_ohm[j]
    });
    let pole: [f64; 3] = std::array::from_fn(|j| duty[j] * y[3] - resistance[j] * abc[j]);
    let va = (2.0_f64 / 3.0).sqrt() * (pole[0] - 0.5 * (pole[1] + pole[2]));
    let vb = (pole[1] - pole[2]) / 2.0_f64.sqrt();
    let (vd, vq) = (cs * va + sn * vb, -sn * va + cs * vb);
    let relative = y[7] - y[8];
    let electric_speed = f64::from(p.pole_pairs) * relative;
    let electromagnetic = f64::from(p.pole_pairs)
        * (p.pm_flux_linkage_wb * y[5] + (p.d_inductance_h - p.q_inductance_h) * y[4] * y[5]);
    let bearing = p.bearing_drag_nm_s_per_rad * relative;
    let strain_rate = y[7] - y[10];
    let transmitted = s.stiffness_nm_per_rad * y[9] + s.damping_nm_s_per_rad * strain_rate;
    let cell = (24.0 - y[1] - y[2]) / b.series_resistance_ohm;
    let lead = (y[2] - y[3]) / LEAD;
    let bridge_current: f64 = duty.iter().zip(abc).map(|(d, i)| d * i).sum();
    [
        -cell,
        (cell - y[1] / b.polarization_resistance_ohm) / b.polarization_capacitance_f,
        (cell - lead) / b.bus_capacitance_f,
        (lead - bridge_current) / bridge.dc_capacitance_f,
        (vd - p.winding_resistance_ohm * y[4] + electric_speed * p.q_inductance_h * y[5])
            / p.d_inductance_h,
        (vq - p.winding_resistance_ohm * y[5]
            - electric_speed * (p.d_inductance_h * y[4] + p.pm_flux_linkage_wb))
            / p.q_inductance_h,
        relative,
        (torque[0] + electromagnetic - bearing - transmitted) / p.rotor_inertia_kg_m2,
        (torque[1] - electromagnetic + bearing) / p.stator_inertia_kg_m2,
        strain_rate,
        (torque[2] + transmitted) / s.load_inertia_kg_m2,
        torque[0] * y[7],
        torque[1] * y[8],
        torque[2] * y[10],
        b.series_resistance_ohm * cell * cell,
        y[1] * y[1] / b.polarization_resistance_ohm,
        LEAD * lead * lead,
        p.winding_resistance_ohm * (y[4] * y[4] + y[5] * y[5]),
        bearing * relative,
        resistance.iter().zip(abc).map(|(r, i)| r * i * i).sum(),
        s.damping_nm_s_per_rad * strain_rate * strain_rate,
    ]
}

fn reference() -> [[f64; 21]; 4] {
    let mut y = [0.0; 21];
    y[..11].copy_from_slice(&INITIAL);
    let mut boundaries = [[0.0; 21]; 4];
    let dt = SEGMENT / 1024.0;
    for segment in 0..4 {
        for _ in 0..1024 {
            let f = |z| rhs(z, DUTIES[segment], TORQUES[segment]);
            let a = f(y);
            let b = f(std::array::from_fn(|i| y[i] + 0.5 * dt * a[i]));
            let c = f(std::array::from_fn(|i| y[i] + 0.5 * dt * b[i]));
            let d = f(std::array::from_fn(|i| y[i] + dt * c[i]));
            y = std::array::from_fn(|i| y[i] + dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0);
        }
        boundaries[segment] = y;
    }
    boundaries
}

fn physical(x: [f64; 11]) -> [f64; 11] {
    let (p, b) = (machine(), battery_parameters());
    let (sn, cs) = (f64::from(p.pole_pairs) * x[6]).sin_cos();
    [
        x[0],
        x[1] / b.polarization_capacitance_f,
        x[2] / b.bus_capacitance_f,
        x[3] / bridge_parameters().dc_capacitance_f,
        (cs * x[4] + sn * x[5] - p.pm_flux_linkage_wb) / p.d_inductance_h,
        (-sn * x[4] + cs * x[5]) / p.q_inductance_h,
        x[6],
        x[7] / p.rotor_inertia_kg_m2,
        x[8] / p.stator_inertia_kg_m2,
        x[9],
        x[10] / shaft().load_inertia_kg_m2,
    ]
}

#[test]
fn finite_battery_motor_and_elastic_load_refine_together_in_all_states_work_and_heat() {
    let expected = reference();
    let mut errors = Vec::new();
    for subdivisions in [8, 16, 32] {
        let (mut model, mut state) = prepare(machine(), shaft(), INITIAL);
        let dt = SEGMENT / f64::from(subdivisions);
        let mut sums = [0.0; 10];
        let mut error = 0.0_f64;
        for segment in 0..4 {
            for _ in 0..subdivisions {
                let old_h = model.system().hamiltonian(&state);
                let mut next = [0.0; 11];
                let step = model
                    .step_controlled(
                        &state,
                        DUTIES[segment],
                        TORQUES[segment],
                        dt,
                        &mut next,
                        || false,
                    )
                    .unwrap();
                assert_eq!(model.system().hamiltonian(&state), old_h);
                close(
                    next[7] + next[8] + next[10] - state[7] - state[8] - state[10],
                    dt * TORQUES[segment].iter().sum::<f64>(),
                    3e-10,
                );
                close(
                    next[9] - state[9],
                    dt * (step.discrete_shaft_speed_rad_s[0] - step.discrete_shaft_speed_rad_s[2]),
                    3e-10,
                );
                close(
                    step.native.supplied,
                    step.external_shaft_work_j.iter().sum(),
                    3e-10,
                );
                close(
                    step.inverter.shaft_work_j[0] + step.rotor_shaft_work_j,
                    step.external_shaft_work_j[0],
                    3e-10,
                );
                close(step.shaft_balance_residual_j(), 0.0, 3e-9);
                close(
                    step.load_kinetic_energy_change_j
                        - step.load_shaft_work_j
                        - step.external_shaft_work_j[2],
                    0.0,
                    3e-9,
                );
                close(step.balance_residual_j(), 0.0, 3e-9);
                assert!(step.heat_joules().iter().all(|q| *q >= 0.0));
                let amounts = [
                    step.external_shaft_work_j[0],
                    step.external_shaft_work_j[1],
                    step.external_shaft_work_j[2],
                    step.battery.series_heat_j,
                    step.battery.polarization_heat_j,
                    step.lead_heat_j,
                    step.inverter.copper_heat_j,
                    step.inverter.bearing_heat_j,
                    step.inverter.channel_heat_j.iter().flatten().sum(),
                    step.shaft_heat_j,
                ];
                for i in 0..10 {
                    sums[i] += amounts[i];
                }
                state = next;
            }
            let mut actual = [0.0; 21];
            actual[..11].copy_from_slice(&physical(state));
            actual[11..].copy_from_slice(&sums);
            for i in 0..21 {
                error = error.max(
                    (actual[i] - expected[segment][i]).abs() / (1.0 + expected[segment][i].abs()),
                );
            }
        }
        errors.push(error);
    }
    assert!(errors[0] / errors[1] > 3.2, "{errors:?}");
    assert!(errors[1] / errors[2] > 3.5, "{errors:?}");
    assert!(errors[2] < 0.01, "{errors:?}");
    // The authored load actually moves, stretches the shaft, and changes the
    // finite source; a report-only mechanical attachment cannot pass this path.
    assert!((expected[3][10] - INITIAL[10]).abs() > 0.1);
    assert!((expected[3][9] - INITIAL[9]).abs() > 0.01);
    assert!((expected[3][0] - INITIAL[0]).abs() > 0.01);
}

fn stamp(model: &PreparedElasticShaftDrive) -> ([f64; 3], Vec<f64>, Vec<f64>, Vec<f64>) {
    let (j, r, _) = model.system().structure();
    (
        model.high_fraction(),
        j.to_vec(),
        r.to_vec(),
        model.system().resistance_witness().unwrap().rows().to_vec(),
    )
}

#[test]
fn late_cancellation_preserves_electrical_and_elastic_structure_then_retries_identically() {
    let (mut model, initial) = prepare(machine(), shaft(), INITIAL);
    let (mut reference, _) = prepare(machine(), shaft(), INITIAL);
    let (mut state, mut reference_state) = ([0.0; 11], [0.0; 11]);
    model
        .step_controlled(&initial, DUTIES[0], TORQUES[0], 1e-4, &mut state, || false)
        .unwrap();
    reference
        .step_controlled(
            &initial,
            DUTIES[0],
            TORQUES[0],
            1e-4,
            &mut reference_state,
            || false,
        )
        .unwrap();
    let before = stamp(&model);
    let mut expected = [0.0; 11];
    let mut polls = 0;
    let expected_step = reference
        .step_controlled(
            &reference_state,
            DUTIES[1],
            TORQUES[1],
            1e-4,
            &mut expected,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut next = [123.0; 11];
    let mut calls = 0;
    assert!(
        model
            .step_controlled(&state, DUTIES[1], TORQUES[1], 1e-4, &mut next, || {
                calls += 1;
                calls == polls
            })
            .is_err()
    );
    assert_eq!(calls, polls);
    assert_eq!(next, [123.0; 11]);
    assert_eq!(stamp(&model), before);
    let step = model
        .step_controlled(&state, DUTIES[1], TORQUES[1], 1e-4, &mut next, || false)
        .unwrap();
    assert_eq!(next.map(f64::to_bits), expected.map(f64::to_bits));
    assert_eq!(step.battery, expected_step.battery);
    assert_eq!(step.inverter, expected_step.inverter);
    assert_eq!(step.heat_joules(), expected_step.heat_joules());
    assert_eq!(stamp(&model), stamp(&reference));
}
