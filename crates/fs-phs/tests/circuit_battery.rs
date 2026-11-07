use fs_phs::circuit::battery::{BatteryOcvCurve, BatteryParameters, BatterySource};
use fs_phs::circuit::motor::{DcMotor, DcMotorParameters};
use fs_phs::{
    InterconnectionBudget, PreparedStepError, SolveBudget, StepWorkspace, interconnect_mapped,
};

fn parameters() -> BatteryParameters {
    BatteryParameters {
        series_resistance_ohm: 0.2,
        polarization_resistance_ohm: 0.3,
        polarization_capacitance_f: 2.0,
        bus_capacitance_f: 0.5,
    }
}
fn source(capacity: f64, knots: &[(f64, f64)]) -> BatterySource {
    BatterySource::new(BatteryOcvCurve::new(capacity, knots).unwrap(), parameters()).unwrap()
}
fn workspace(source: &BatterySource) -> StepWorkspace {
    let mut work = StepWorkspace::new(source.system()).unwrap();
    work.set_solve_budget(SolveBudget {
        relative: 2e-12,
        absolute: 2e-13,
        max_iterations: 30,
    })
    .unwrap();
    work
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a} != {b}, error={}", a - b);
}

#[test]
fn finite_ocv_integral_includes_voltage_plateaus_and_crossed_knots() {
    let curve =
        BatteryOcvCurve::new(100.0, &[(0.0, 3.0), (0.5, 3.0), (0.75, 3.6), (1.0, 4.2)]).unwrap();
    close(curve.energy_j(100.0).unwrap(), 330.0, 1e-12);
    close(curve.energy_change_j(40.0, 90.0).unwrap(), 169.2, 1e-12);
    close(curve.energy_change_j(90.0, 40.0).unwrap(), -169.2, 1e-12);
    close(curve.voltage_v(0.0).unwrap(), 3.0, 0.0);
    close(curve.voltage_v(100.0).unwrap(), 4.2, 1e-15);
    assert!(curve.voltage_v(-1e-15).is_err());
    assert!(curve.energy_j(100.0 + 1e-10).is_err());
    assert!(BatteryOcvCurve::new(100.0, &[(0.0, 4.0), (1.0, 3.0)]).is_err());
    assert!(BatteryOcvCurve::new(0.0, &[(0.0, 3.0), (1.0, 4.0)]).is_err());
}

#[test]
fn chemical_discrete_work_and_analytic_tangent_follow_the_original_source() {
    let source = source(100.0, &[(0.0, 3.0), (0.5, 3.0), (0.75, 3.6), (1.0, 4.2)]);
    let rule = source.system().storage().discrete_gradient_rule().unwrap();
    for (a, b) in [(10.0, 40.0), (20.0, 90.0), (90.0, 20.0), (70.0, 70.0)] {
        let before = [a, 0.2, 1.0];
        let after = [b, 0.6, 2.0];
        let mid = core::array::from_fn::<_, 3, _>(|i| f64::midpoint(before[i], after[i]));
        let mut effort = [0.0; 3];
        rule.evaluate(&before, &after, &mid, &mut effort);
        close(
            effort[0] * (b - a),
            source.curve().energy_change_j(a, b).unwrap(),
            2e-12,
        );
        let mut tangent = [0.0; 3];
        assert!(rule.end_jacobian_action(&before, &after, &[1.0, 0.0, 0.0], &mut tangent));
        let h = 1e-4;
        let mut plus = after;
        let mut minus = after;
        plus[0] += h;
        minus[0] -= h;
        let mut e_plus = [0.0; 3];
        let mut e_minus = [0.0; 3];
        rule.evaluate(
            &before,
            &plus,
            &core::array::from_fn::<_, 3, _>(|i| f64::midpoint(before[i], plus[i])),
            &mut e_plus,
        );
        rule.evaluate(
            &before,
            &minus,
            &core::array::from_fn::<_, 3, _>(|i| f64::midpoint(before[i], minus[i])),
            &mut e_minus,
        );
        close(tangent[0], (e_plus[0] - e_minus[0]) / (2.0 * h), 2e-10);
    }
    let continuous = source.system().storage().continuous_tangent().unwrap();
    let direction = [2.0, 3.0, 4.0];
    let mut action = [0.0; 3];
    // Plateau, smooth ramp and an exactly collinear interior knot.
    for (charge, slope) in [(25.0, 0.0), (60.0, 0.024), (75.0, 0.024)] {
        assert!(continuous.apply(&[charge, 0.2, 1.0], &direction, &mut action));
        close(action[0], slope * direction[0], 1e-15);
        close(action[1], 1.5, 0.0);
        close(action[2], 8.0, 0.0);
    }
    for (charge, slope) in [(0.0, 0.0), (50.0, 0.024), (100.0, 0.024)] {
        let state = [charge, 0.2, 1.0];
        // A zero direction cannot manufacture a missing two-sided derivative.
        assert!(!continuous.apply(&state, &direction, &mut action));
        assert!(!continuous.apply(&state, &[0.0; 3], &mut action));
        assert!(rule.end_jacobian_action(&state, &state, &direction, &mut action));
        close(action[0], 0.5 * slope * direction[0], 1e-15);
    }
}

#[test]
fn discharge_sags_and_regeneration_recharges_without_creating_energy() {
    let source = source(10.0, &[(0.0, 3.0), (0.5, 3.4), (1.0, 4.0)]);
    let mut work = workspace(&source);
    let mut state = source.state(0.8, 0.0, 3.76).unwrap();
    let initial_charge = state[0];
    let dt = 0.01;
    let mut total_heat = 0.0;
    for current in [-0.6, 0.6] {
        let previous_charge = state[0];
        for _ in 0..300 {
            let mut next = [0.0; 3];
            let mut voltage = [0.0];
            let record = work
                .step_into(
                    source.system(),
                    &state,
                    &[current],
                    dt,
                    &mut next,
                    &mut voltage,
                )
                .unwrap();
            let midpoint = core::array::from_fn::<_, 3, _>(|i| f64::midpoint(state[i], next[i]));
            let mut efforts = [0.0; 3];
            source
                .system()
                .storage()
                .discrete_gradient_rule()
                .unwrap()
                .evaluate(&state, &next, &midpoint, &mut efforts);
            let i_cell =
                (efforts[0] - efforts[1] - efforts[2]) / parameters().series_resistance_ohm;
            let heat = dt
                * (parameters().series_resistance_ohm * i_cell * i_cell
                    + efforts[1] * efforts[1] / parameters().polarization_resistance_ohm);
            let chemical = source.curve().energy_change_j(state[0], next[0]).unwrap();
            let polarization =
                (next[1] - state[1]) * midpoint[1] / parameters().polarization_capacitance_f;
            let bus = (next[2] - state[2]) * midpoint[2] / parameters().bus_capacitance_f;
            close(
                chemical + polarization + bus + heat,
                dt * current * voltage[0],
                3e-11,
            );
            close(
                (next[0] - state[0]) + (next[2] - state[2]),
                dt * current,
                3e-12,
            );
            close(record.dissipated, heat, 3e-12);
            assert!(heat >= 0.0);
            total_heat += heat;
            state = next;
        }
        let observation = source.observe(&state).unwrap();
        if current < 0.0 {
            assert!(state[0] < previous_charge);
            assert!(observation.terminal_voltage_v < observation.open_circuit_voltage_v);
            assert!(observation.polarization_voltage_v > 0.0);
        } else {
            assert!(state[0] > previous_charge);
            assert!(observation.terminal_voltage_v > observation.open_circuit_voltage_v);
            assert!(observation.cell_current_a < 0.0);
        }
    }
    assert!(total_heat > 0.1);
    // The bus also ends charged above its original voltage: chemical charge
    // therefore remains below its initial value after equal port coulombs.
    assert!(state[0] < initial_charge);

    let quiet = BatterySource::new(
        BatteryOcvCurve::new(1.0, &[(0.0, 4.0), (1.0, 4.0)]).unwrap(),
        parameters(),
    )
    .unwrap();
    let observed = quiet
        .observe(&quiet.state(0.5, 1e-18, 4.0).unwrap())
        .unwrap();
    let expected_current = -1e-18 / parameters().series_resistance_ohm;
    let expected_heat = parameters().series_resistance_ohm * expected_current * expected_current;
    assert!(observed.cell_current_a < 0.0 && observed.series_heat_w > 0.0);
    close(observed.cell_current_a / expected_current, 1.0, 4e-15);
    close(observed.series_heat_w / expected_heat, 1.0, 8e-15);

    let large = BatterySource::new(
        quiet.curve().clone(),
        BatteryParameters {
            series_resistance_ohm: 1e160,
            polarization_resistance_ohm: 1e160,
            polarization_capacitance_f: 1e-160,
            ..parameters()
        },
    )
    .unwrap();
    let observed = large
        .observe(&large.state(0.5, 1e160, 4.0).unwrap())
        .unwrap();
    // V_pol² overflows, but its resistor's physical heat is representable.
    close(observed.polarization_heat_w / 1e160, 1.0, 4e-15);
}

#[test]
fn constant_ocv_transient_converges_to_independent_two_rc_solution() {
    let source = source(100.0, &[(0.0, 4.0), (1.0, 4.0)]);
    let p = parameters();
    let input = -1.0;
    let duration = 0.8;
    let a = -(1.0 / p.series_resistance_ohm + 1.0 / p.polarization_resistance_ohm)
        / p.polarization_capacitance_f;
    let b = -1.0 / (p.series_resistance_ohm * p.polarization_capacitance_f);
    let c = -1.0 / (p.series_resistance_ohm * p.bus_capacitance_f);
    let d = c;
    let trace = 0.5 * (a + d);
    let eigen_split = ((0.5 * (a - d)).powi(2) + b * c).sqrt();
    let scale = (trace * duration).exp();
    let identity = (eigen_split * duration).cosh();
    let linear = (eigen_split * duration).sinh() / eigen_split;
    let steady_vp = -input * p.polarization_resistance_ohm;
    let steady_vb = 4.0 + input * (p.series_resistance_ohm + p.polarization_resistance_ohm);
    let z = [-steady_vp, 4.0 - steady_vb];
    let exact = [
        steady_vp + scale * (identity * z[0] + linear * ((a - trace) * z[0] + b * z[1])),
        steady_vb + scale * (identity * z[1] + linear * (c * z[0] + (d - trace) * z[1])),
    ];
    let run = |steps: usize| {
        let mut work = workspace(&source);
        let mut state = source.state(0.5, 0.0, 4.0).unwrap();
        for _ in 0..steps {
            let mut next = [0.0; 3];
            let mut output = [0.0];
            work.step_into(
                source.system(),
                &state,
                &[input],
                duration / steps as f64,
                &mut next,
                &mut output,
            )
            .unwrap();
            state = next;
        }
        let observed = source.observe(&state).unwrap();
        (observed.polarization_voltage_v - exact[0]).abs()
            + (observed.terminal_voltage_v - exact[1]).abs()
    };
    let coarse = run(20);
    let fine = run(40);
    assert!(
        fine < coarse * 0.27 && fine > coarse * 0.22,
        "coarse={coarse}, fine={fine}"
    );
    assert!(fine < 2e-5);
}

#[test]
fn bounded_charge_and_cancellation_refuse_transactionally() {
    let source = source(1.0, &[(0.0, 3.0), (1.0, 4.0)]);
    let mut work = workspace(&source);
    let empty = source.state(0.0, 0.0, 0.0).unwrap();
    let mut next = [77.0; 3];
    let mut voltage = [88.0];
    assert!(
        work.step_into(
            source.system(),
            &empty,
            &[0.0],
            0.01,
            &mut next,
            &mut voltage
        )
        .is_err()
    );
    assert_eq!(next, [77.0; 3]);
    assert_eq!(voltage, [88.0]);
    let charged = source.state(0.5, 0.0, 3.5).unwrap();
    assert!(matches!(
        work.step_into_controlled(
            source.system(),
            &charged,
            &[0.0],
            0.01,
            &mut next,
            &mut voltage,
            || true
        ),
        Err(PreparedStepError::Cancelled)
    ));
    assert_eq!(next, [77.0; 3]);
    assert_eq!(voltage, [88.0]);
    work.step_into(
        source.system(),
        &charged,
        &[0.0],
        0.01,
        &mut next,
        &mut voltage,
    )
    .unwrap();
    assert_eq!(next, charged);
    close(voltage[0], 3.5, 1e-14);
}

#[test]
fn native_battery_motor_connection_drives_and_regenerates_with_full_power_balance() {
    let motor_parameters = DcMotorParameters {
        inductance_h: 0.05,
        winding_resistance_ohm: 0.1,
        rotor_inertia_kg_m2: 0.2,
        stator_inertia_kg_m2: 0.3,
        coupling_nm_per_a: 0.8,
        bearing_drag_nm_s_per_rad: 0.001,
    };
    let battery_parameters = parameters();
    let dt = 0.002;
    for generating in [false, true] {
        let battery = source(100.0, &[(0.0, 3.0), (0.5, 3.6), (1.0, 4.2)]);
        let curve = battery.curve().clone();
        let battery_state = battery.state(0.6, 0.0, 3.72).unwrap();
        let motor = DcMotor::new(motor_parameters).unwrap();
        let rotor_speed = if generating { 10.0 } else { 0.0 };
        let motor_state = motor
            .state(
                0.0,
                [
                    rotor_speed,
                    -motor_parameters.rotor_inertia_kg_m2 / motor_parameters.stator_inertia_kg_m2
                        * rotor_speed,
                ],
            )
            .unwrap();
        let initial = [
            battery_state[0],
            battery_state[1],
            battery_state[2],
            motor_state[0],
            motor_state[1],
            motor_state[2],
        ];
        // Battery bus voltage drives the winding; winding current leaves the
        // battery. Only the two physical shaft-torque ports remain external.
        let system = interconnect_mapped(
            battery.into_system(),
            motor.into_system(),
            &[0],
            &[0],
            &[1.0],
            InterconnectionBudget {
                max_states: 6,
                max_total_ports: 4,
                max_dense_bytes: 1 << 20,
                max_coupling_multiply_terms: 1_000_000,
            },
        )
        .unwrap();
        assert_eq!(system.state_dim(), 6);
        assert_eq!(system.port_dim(), 2);
        let mut work = StepWorkspace::new(&system).unwrap();
        work.set_solve_budget(SolveBudget {
            relative: 2e-12,
            absolute: 2e-13,
            max_iterations: 30,
        })
        .unwrap();
        let shaft_torque = if generating {
            [0.0, 0.0]
        } else {
            [-0.02, 0.02]
        };
        let mut state = initial;
        let mut electrical_work = 0.0;
        let mut shaft_work = 0.0;
        let mut total_heat = 0.0;
        for _ in 0..100 {
            let (mut next, mut shaft_speed) = ([0.0; 6], [0.0; 2]);
            let record = work
                .step_into(
                    &system,
                    &state,
                    &shaft_torque,
                    dt,
                    &mut next,
                    &mut shaft_speed,
                )
                .unwrap();
            let midpoint = core::array::from_fn::<_, 6, _>(|i| f64::midpoint(state[i], next[i]));
            // Independent original-component accounts, including the internal
            // electrical work that disappears from the composite port ledger.
            let chemical = curve.energy_change_j(state[0], next[0]).unwrap();
            let ocv = if next[0] == state[0] {
                curve.voltage_v(state[0]).unwrap()
            } else {
                chemical / (next[0] - state[0])
            };
            let polarization_v = midpoint[1] / battery_parameters.polarization_capacitance_f;
            let bus_v = midpoint[2] / battery_parameters.bus_capacitance_f;
            let cell_current =
                (ocv - polarization_v - bus_v) / battery_parameters.series_resistance_ohm;
            let motor_current = midpoint[3] / motor_parameters.inductance_h;
            let rotor_speed = midpoint[4] / motor_parameters.rotor_inertia_kg_m2;
            let stator_speed = midpoint[5] / motor_parameters.stator_inertia_kg_m2;
            let relative_speed = rotor_speed - stator_speed;
            let battery_heat = dt
                * (battery_parameters.series_resistance_ohm * cell_current.powi(2)
                    + polarization_v.powi(2) / battery_parameters.polarization_resistance_ohm);
            let motor_heat = dt
                * (motor_parameters.winding_resistance_ohm * motor_current.powi(2)
                    + motor_parameters.bearing_drag_nm_s_per_rad * relative_speed.powi(2));
            let battery_delta =
                chemical + (next[1] - state[1]) * polarization_v + (next[2] - state[2]) * bus_v;
            let motor_delta = (next[3] - state[3]) * motor_current
                + (next[4] - state[4]) * rotor_speed
                + (next[5] - state[5]) * stator_speed;
            let transferred = dt * bus_v * motor_current;
            let external_work =
                dt * (shaft_torque[0] * rotor_speed + shaft_torque[1] * stator_speed);
            close(battery_delta + battery_heat + transferred, 0.0, 3e-10);
            close(motor_delta + motor_heat - transferred, external_work, 3e-10);
            close(
                battery_delta + motor_delta + battery_heat + motor_heat,
                external_work,
                3e-10,
            );
            close(record.delta_h, battery_delta + motor_delta, 3e-10);
            close(record.dissipated, battery_heat + motor_heat, 3e-10);
            close(record.supplied, external_work, 3e-12);
            close(
                (next[0] - state[0]) + (next[2] - state[2]),
                -dt * motor_current,
                3e-12,
            );
            close(next[4] + next[5], initial[4] + initial[5], 3e-12);
            electrical_work += transferred;
            shaft_work += external_work;
            total_heat += battery_heat + motor_heat;
            state = next;
        }
        assert!(total_heat > 0.0);
        if generating {
            assert!(state[0] > initial[0]);
            assert!(state[3] < 0.0);
            assert!(electrical_work < 0.0);
            assert!(state[4] < initial[4]);
            assert_eq!(shaft_work, 0.0);
        } else {
            assert!(state[0] < initial[0]);
            assert!(state[3] > 0.0);
            assert!(electrical_work > 0.0);
            assert!(state[4] > 0.0 && state[5] < 0.0);
            assert!(shaft_work < 0.0);
        }
    }
}
