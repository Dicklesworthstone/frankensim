//! Independent finite-source circuit oracle, motoring/regeneration and storage.
use super::{budget, close, machine};
use fs_phs::circuit::{
    battery::{BatteryOcvCurve, BatteryParameters, BatterySource},
    inductance::InductorFluxCurve,
    synchronous::{
        SynchronousMotor,
        inverter::{
            AveragedSynchronousInverter, InverterParameters,
            battery::{BatteryInverterControls, PreparedBatterySynchronousInverter},
        },
    },
};

#[path = "battery/shaft.rs"]
mod shaft;

const LEAD: f64 = 0.2;

fn battery_parameters() -> BatteryParameters {
    BatteryParameters {
        series_resistance_ohm: 0.3,
        polarization_resistance_ohm: 2.0,
        polarization_capacitance_f: 0.5,
        bus_capacitance_f: 0.04,
    }
}

fn bridge_parameters() -> InverterParameters {
    InverterParameters {
        dc_capacitance_f: 0.03,
        high_on_resistance_ohm: [0.025; 3],
        low_on_resistance_ohm: [0.015; 3],
    }
}

fn controls() -> BatteryInverterControls {
    BatteryInverterControls {
        solve: budget(),
        backtracking_halvings: 10,
        absolute_energy_tolerance_j: 2e-9,
        relative_energy_tolerance: 2e-9,
    }
}

fn plateau() -> BatterySource {
    BatterySource::new(
        BatteryOcvCurve::new(100.0, &[(0.0, 24.0), (1.0, 24.0)]).unwrap(),
        battery_parameters(),
    )
    .unwrap()
}

fn drive(generator: bool) -> (PreparedBatterySynchronousInverter, [f64; 9]) {
    let motor = SynchronousMotor::new(machine()).unwrap();
    let machine_state = motor
        .state(
            [0.0, if generator { -2.0 } else { 2.0 }],
            0.0,
            [if generator { 100.0 } else { 5.0 }, 0.0],
        )
        .unwrap();
    let model = AveragedSynchronousInverter::new(motor, bridge_parameters())
        .unwrap()
        .with_battery_source(plateau(), LEAD)
        .unwrap();
    let state = if generator {
        model.state(0.6, -0.1, 24.4, 25.0, &machine_state)
    } else {
        model.state(0.6, 0.1, 23.5, 23.2, &machine_state)
    }
    .unwrap();
    (
        model.prepare_controlled(controls(), || false).unwrap(),
        state,
    )
}

// Independent Kirchhoff equations in VOLTAGES: polarization, battery bus, DC
// bus. On a constant OCV plateau, their midpoint update is a 3x3 linear solve.
// No native storage, constitutive gradient, port map or resistance is inspected.
fn circuit_midpoint(v: [f64; 3], dt: f64) -> [f64; 3] {
    let p = battery_parameters();
    let (gs, gp, gl) = (
        1.0 / p.series_resistance_ohm,
        1.0 / p.polarization_resistance_ohm,
        1.0 / LEAD,
    );
    let (cp, cb, cd) = (
        p.polarization_capacitance_f,
        p.bus_capacitance_f,
        bridge_parameters().dc_capacitance_f,
    );
    let a = [
        [-(gs + gp) / cp, -gs / cp, 0.0],
        [-gs / cb, -(gs + gl) / cb, gl / cb],
        [0.0, gl / cd, -gl / cd],
    ];
    let source = [24.0 * gs / cp, 24.0 * gs / cb, 0.0];
    let mut matrix = [[0.0; 4]; 3];
    for i in 0..3 {
        for j in 0..3 {
            matrix[i][j] = if i == j { 1.0 } else { 0.0 } - 0.5 * dt * a[i][j];
        }
        matrix[i][3] =
            v[i] + dt * source[i] + 0.5 * dt * (0..3).map(|j| a[i][j] * v[j]).sum::<f64>();
    }
    for k in 0..3 {
        let pivot = matrix[k][k];
        for j in k..4 {
            matrix[k][j] /= pivot;
        }
        let normalized = matrix[k];
        for i in 0..3 {
            if i != k {
                let multiplier = matrix[i][k];
                for j in k..4 {
                    matrix[i][j] -= multiplier * normalized[j];
                }
            }
        }
    }
    [matrix[0][3], matrix[1][3], matrix[2][3]]
}

#[test]
fn finite_source_inrush_matches_independent_polarization_and_two_bus_circuit() {
    let motor = SynchronousMotor::new(machine()).unwrap();
    let machine_state = motor.state([0.0; 2], 0.0, [0.0; 2]).unwrap();
    let model = AveragedSynchronousInverter::new(motor, bridge_parameters())
        .unwrap()
        .with_battery_source(plateau(), LEAD)
        .unwrap();
    let mut state = model.state(0.6, 0.4, 22.0, 17.0, &machine_state).unwrap();
    let initial = state;
    let mut model = model.prepare_controlled(controls(), || false).unwrap();
    let p = battery_parameters();
    let (mut voltage, mut chemical) = ([0.4, 22.0, 17.0], state[0]);
    let dt = 0.002;
    for k in 0..20 {
        let expected = circuit_midpoint(voltage, dt);
        let mean = std::array::from_fn::<_, 3, _>(|i| 0.5 * (voltage[i] + expected[i]));
        let cell = (24.0 - mean[0] - mean[1]) / p.series_resistance_ohm;
        let lead = (mean[1] - mean[2]) / LEAD;
        chemical -= dt * cell;
        let mut next = [0.0; 9];
        let step = model
            .step_controlled(
                &state,
                [if k % 2 == 0 { 0.0 } else { 1.0 }; 3],
                [0.0; 2],
                dt,
                &mut next,
                || false,
            )
            .unwrap();
        close(next[0], chemical, 3e-11);
        close(next[1] / p.polarization_capacitance_f, expected[0], 3e-10);
        close(next[2] / p.bus_capacitance_f, expected[1], 3e-10);
        close(
            next[3] / bridge_parameters().dc_capacitance_f,
            expected[2],
            3e-10,
        );
        close(step.lead_current_a, lead, 3e-10);
        close(step.battery.discrete_cell_current_a, cell, 3e-10);
        close(
            step.battery.chemical_energy_change_j,
            -dt * 24.0 * cell,
            3e-10,
        );
        close(
            step.battery.series_heat_j,
            dt * p.series_resistance_ohm * cell * cell,
            3e-10,
        );
        close(
            step.battery.polarization_heat_j,
            dt * mean[0] * mean[0] / p.polarization_resistance_ohm,
            3e-11,
        );
        close(step.lead_heat_j, dt * LEAD * lead * lead, 3e-10);
        close(step.inverter.source_work_j, dt * mean[2] * lead, 3e-10);
        close(
            next[0] + next[2] + next[3],
            initial[0] + initial[2] + initial[3],
            3e-11,
        );
        for i in 4..9 {
            close(next[i], initial[i], 3e-11);
        }
        assert!(step.inverter.heat_joules().iter().all(|h| h.abs() < 1e-20));
        close(step.balance_residual_j(), 0.0, 3e-9);
        state = next;
        voltage = expected;
    }
    assert!(state[0] < initial[0]);
    assert!(state[3] > initial[3]);
    assert!(state[1] > initial[1]);
}

#[test]
fn changing_duty_motors_and_regenerates_through_finite_cell_and_both_buses() {
    for generator in [false, true] {
        let (mut model, mut state) = drive(generator);
        let initial = state;
        let h0 = model.system().hamiltonian(&initial);
        let (mut heat, mut bridge_work, mut machine_work) = (0.0, 0.0, 0.0);
        for k in 0..100 {
            // Authored changing q-axis command, expressed as three physical leg
            // fractions. The original machine supplies its angle, not its work.
            let angle = f64::from(machine().pole_pairs) * state[6];
            let (s, c) = angle.sin_cos();
            let magnitude = [0.15, 0.25, 0.4, 0.2][k / 25];
            let (alpha, beta) = (-s * magnitude, c * magnitude);
            let duty = [
                0.5 + (2.0_f64 / 3.0).sqrt() * alpha,
                0.5 - alpha / 6.0_f64.sqrt() + beta / 2.0_f64.sqrt(),
                0.5 - alpha / 6.0_f64.sqrt() - beta / 2.0_f64.sqrt(),
            ];
            let h_before = model.system().hamiltonian(&state);
            let mut next = [0.0; 9];
            let step = model
                .step_controlled(&state, duty, [0.0; 2], 1e-4, &mut next, || false)
                .unwrap();
            // A new modulation never changes the original Hamiltonian at an
            // unchanged state, nor borrows energy from an ideal DC supply.
            assert_eq!(h_before, model.system().hamiltonian(&state));
            assert_eq!(model.high_fraction(), duty);
            assert_eq!(step.native.supplied, 0.0);
            assert_eq!(step.inverter.shaft_work_j, [0.0; 2]);
            assert_eq!(step.lead_current_a < 0.0, generator);
            assert_eq!(step.battery.discrete_cell_current_a < 0.0, generator);
            let [vb, vd] = step.discrete_bus_voltage_v;
            close(step.lead_current_a, (vb - vd) / LEAD, 2e-12);
            close(
                step.battery.terminal_work_j + step.inverter.source_work_j + step.lead_heat_j,
                0.0,
                3e-12,
            );
            close(
                step.inverter.bridge_work_j - step.inverter.machine_electrical_work_j,
                step.inverter.channel_heat_j.iter().flatten().sum(),
                3e-12,
            );
            let impulse = step.inverter.electromagnetic_angular_impulse_nm_s
                + step.inverter.bearing_angular_impulse_nm_s;
            close(next[7] - state[7], impulse, 3e-11);
            close(next[8] - state[8], -impulse, 3e-11);
            close(
                next[6] - state[6],
                step.inverter.shaft_angle_increment_rad[0]
                    - step.inverter.shaft_angle_increment_rad[1],
                3e-11,
            );
            close(next[7] + next[8], initial[7] + initial[8], 3e-10);
            assert!(step.heat_joules().iter().all(|q| *q >= 0.0));
            let interval_heat: f64 = step.heat_joules().iter().sum();
            close(interval_heat, step.native.dissipated, 3e-11);
            close(step.battery.balance_residual_j(), 0.0, 3e-9);
            close(step.balance_residual_j(), 0.0, 3e-9);
            heat += interval_heat;
            bridge_work += step.inverter.bridge_work_j;
            machine_work += step.inverter.machine_electrical_work_j;
            state = next;
        }
        assert_eq!(state[0] > initial[0], generator);
        assert_eq!(state[7] < initial[7], generator);
        assert_eq!(bridge_work < 0.0, generator);
        assert_eq!(machine_work < 0.0, generator);
        close(model.system().hamiltonian(&state) - h0 + heat, 0.0, 3e-8);
    }
}

#[test]
fn original_saturated_flux_storage_and_crossed_ocv_knot_survive_composition() {
    let d = InductorFluxCurve::try_new(vec![
        (-0.4, -60.0),
        (-0.03, -1.5),
        (0.0, 0.0),
        (0.03, 1.5),
        (0.4, 60.0),
    ])
    .unwrap();
    let q = InductorFluxCurve::try_new(vec![
        (-0.5, -70.0),
        (-0.04, -2.0),
        (0.0, 0.0),
        (0.04, 2.0),
        (0.5, 70.0),
    ])
    .unwrap();
    let motor = SynchronousMotor::with_flux_curves(machine(), d.clone(), q.clone()).unwrap();
    let machine_state = motor.state([1.2, 2.5], 0.27, [10.0, -1.0]).unwrap();
    let battery = BatterySource::new(
        BatteryOcvCurve::new(1.0, &[(0.0, 20.0), (0.5, 24.0), (1.0, 25.0)]).unwrap(),
        battery_parameters(),
    )
    .unwrap();
    let model = AveragedSynchronousInverter::new(motor, bridge_parameters())
        .unwrap()
        .with_battery_source(battery, LEAD)
        .unwrap();
    let initial = model
        .state(0.50005, 0.1, 23.0, 23.0, &machine_state)
        .unwrap();
    let mut model = model.prepare_controlled(controls(), || false).unwrap();
    let point = model.machine_operating_point(&initial).unwrap();
    close(point.current_dq_a[0], 1.2, 3e-13);
    close(point.current_dq_a[1], 2.5, 3e-13);
    assert!((point.differential_inductance_dq_h[1] - machine().q_inductance_h).abs() > 0.01);
    let mut next = [0.0; 9];
    let step = model
        .step_controlled(
            &initial,
            [0.3, 0.8, 0.1],
            [0.04, -0.01],
            1e-4,
            &mut next,
            || false,
        )
        .unwrap();
    assert!(next[0] < 0.5 && initial[0] > 0.5);
    // Independent trapezoids on BOTH original OCV segments, not midpoint OCV.
    let chemical = -(initial[0] - 0.5) * (24.0 + (23.0 + 2.0 * initial[0])) * 0.5
        - (0.5 - next[0]) * ((20.0 + 8.0 * next[0]) + 24.0) * 0.5;
    close(step.battery.chemical_energy_change_j, chemical, 3e-13);
    let physical_energy = |x: &[f64; 9]| {
        let (s, c) = (f64::from(machine().pole_pairs) * x[6]).sin_cos();
        d.stored_energy_j(c * x[4] + s * x[5] - machine().pm_flux_linkage_wb)
            .unwrap()
            + q.stored_energy_j(-s * x[4] + c * x[5]).unwrap()
            + 0.5 * x[3] * x[3] / bridge_parameters().dc_capacitance_f
            + 0.5 * x[7] * x[7] / machine().rotor_inertia_kg_m2
            + 0.5 * x[8] * x[8] / machine().stator_inertia_kg_m2
    };
    close(
        step.inverter_storage_change_j,
        physical_energy(&next) - physical_energy(&initial),
        3e-11,
    );
    close(
        next[7] + next[8] - initial[7] - initial[8],
        1e-4 * (0.04 - 0.01),
        3e-11,
    );
    close(step.balance_residual_j(), 0.0, 3e-9);
}

fn stamp(model: &PreparedBatterySynchronousInverter) -> ([f64; 3], Vec<f64>, Vec<f64>, Vec<f64>) {
    let (j, r, _) = model.system().structure();
    (
        model.high_fraction(),
        j.to_vec(),
        r.to_vec(),
        model.system().resistance_witness().unwrap().rows().to_vec(),
    )
}

#[test]
fn refused_duty_soc_bus_and_late_cancellation_restore_the_full_native_stamp() {
    let (mut model, initial) = drive(false);
    let (mut reference, _) = drive(false);
    let duty = [0.2, 0.7, 0.4];
    let (mut warm, mut warm_reference) = ([0.0; 9], [0.0; 9]);
    model
        .step_controlled(&initial, duty, [0.0; 2], 1e-4, &mut warm, || false)
        .unwrap();
    reference
        .step_controlled(&initial, duty, [0.0; 2], 1e-4, &mut warm_reference, || {
            false
        })
        .unwrap();
    let held = stamp(&model);
    let mut unsupported = warm;
    unsupported[0] = -1e-8;
    let mut negative_bus = warm;
    negative_bus[3] = -1e-8;
    for (state, command) in [
        (&warm, [0.2, 1.01, 0.3]),
        (&unsupported, duty),
        (&negative_bus, duty),
    ] {
        let mut next = [123.0; 9];
        assert!(
            model
                .step_controlled(state, command, [0.0; 2], 1e-4, &mut next, || false)
                .is_err()
        );
        assert_eq!(next, [123.0; 9]);
        assert_eq!(stamp(&model), held);
    }
    let new_duty = [0.6, 0.8, 0.15];
    let mut expected = [0.0; 9];
    let mut polls = 0;
    let expected_step = reference
        .step_controlled(
            &warm_reference,
            new_duty,
            [0.0; 2],
            1e-4,
            &mut expected,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut next = [123.0; 9];
    let mut calls = 0;
    assert!(
        model
            .step_controlled(&warm, new_duty, [0.0; 2], 1e-4, &mut next, || {
                calls += 1;
                calls == polls
            })
            .is_err()
    );
    assert_eq!(calls, polls);
    assert_eq!(next, [123.0; 9]);
    assert_eq!(stamp(&model), held);
    let retry = model
        .step_controlled(&warm, new_duty, [0.0; 2], 1e-4, &mut next, || false)
        .unwrap();
    assert_eq!(next.map(f64::to_bits), expected.map(f64::to_bits));
    assert_eq!(retry.inverter, expected_step.inverter);
    assert_eq!(retry.battery, expected_step.battery);
    assert_eq!(stamp(&model), stamp(&reference));
}
