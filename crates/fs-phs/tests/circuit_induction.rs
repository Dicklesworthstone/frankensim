//! Original-flux induction dynamics, slip losses and reciprocal finite-bus work.
use fs_phs::circuit::induction::{
    InductionMotor, InductionMotorParameters, phase_voltage_to_alpha_beta,
};
use fs_phs::{
    InterconnectionBudget, PortHamiltonian, PreparedStepError, QuadraticStorage, SolveBudget,
    StepWorkspace, interconnect_mapped,
};

fn parameters() -> InductionMotorParameters {
    InductionMotorParameters {
        stator_leakage_inductance_h: 0.004,
        rotor_leakage_inductance_h: 0.006,
        magnetizing_inductance_h: 0.05,
        stator_resistance_ohm: 0.8,
        rotor_resistance_ohm: 0.5,
        pole_pairs: 2,
        rotor_inertia_kg_m2: 0.1,
        stator_inertia_kg_m2: 0.2,
        bearing_drag_nm_s_per_rad: 0.0,
    }
}
fn budget() -> SolveBudget {
    SolveBudget {
        relative: 2e-12,
        absolute: 2e-13,
        max_iterations: 30,
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}; tolerance={tolerance:e}"
    );
}
fn rotate(angle: f64, x: [f64; 2]) -> [f64; 2] {
    let (s, c) = angle.sin_cos();
    [c * x[0] - s * x[1], s * x[0] + c * x[1]]
}

// Independent direct inverse of the original coupled winding matrix.
fn energy(p: InductionMotorParameters, x: &[f64; 7]) -> f64 {
    let ls = p.stator_leakage_inductance_h + p.magnetizing_inductance_h;
    let lr = p.rotor_leakage_inductance_h + p.magnetizing_inductance_h;
    let lm = p.magnetizing_inductance_h;
    let determinant = ls * lr - lm * lm;
    let rotor = rotate(f64::from(p.pole_pairs) * x[4], [x[2], x[3]]);
    0.5 * ((lr * (x[0] * x[0] + x[1] * x[1]) + ls * (x[2] * x[2] + x[3] * x[3])
        - 2.0 * lm * (x[0] * rotor[0] + x[1] * rotor[1]))
        / determinant
        + x[5] * x[5] / p.rotor_inertia_kg_m2
        + x[6] * x[6] / p.stator_inertia_kg_m2)
}

#[test]
fn original_winding_energy_gives_phase_power_and_reciprocal_torque() {
    for pole_pairs in [1, 3] {
        let p = InductionMotorParameters {
            pole_pairs,
            ..parameters()
        };
        let model = InductionMotor::new(p).unwrap();
        for theta in [-0.4, 0.0, 1.1] {
            let is = [2.3, -1.4];
            let ir = [-0.9, 1.7];
            let state = model.state(is, ir, theta, [9.0, -2.0]).unwrap();
            let point = model.operating_point(&state).unwrap();
            let ir_ab = rotate(f64::from(pole_pairs) * theta, ir);
            let torque = f64::from(pole_pairs)
                * p.magnetizing_inductance_h
                * (ir_ab[0] * is[1] - ir_ab[1] * is[0]);
            for i in 0..2 {
                close(point.stator_current_alpha_beta_a[i], is[i], 5e-14);
                close(point.rotor_current_dq_a[i], ir[i], 5e-14);
                close(point.rotor_current_alpha_beta_a[i], ir_ab[i], 5e-14);
            }
            close(point.electromagnetic_torque_nm, torque, 1e-14);
            close(model.system().hamiltonian(&state), energy(p, &state), 3e-14);
            let h = 1e-5;
            let mut plus = state;
            let mut minus = state;
            plus[4] += h;
            minus[4] -= h;
            close(
                -(energy(p, &plus) - energy(p, &minus)) / (2.0 * h),
                torque,
                3e-9,
            );
            let phase_voltage = [7.0, -4.0, 2.0];
            let voltage = phase_voltage_to_alpha_beta(phase_voltage).unwrap();
            close(point.stator_phase_current_a.iter().sum(), 0.0, 1e-14);
            close(
                phase_voltage
                    .iter()
                    .zip(point.stator_phase_current_a)
                    .map(|(v, i)| v * i)
                    .sum(),
                voltage[0] * is[0] + voltage[1] * is[1],
                2e-13,
            );
            close(
                point.stator_copper_power_w,
                p.stator_resistance_ohm * (is[0] * is[0] + is[1] * is[1]),
                2e-13,
            );
            close(
                point.rotor_copper_power_w,
                p.rotor_resistance_ohm * (ir[0] * ir[0] + ir[1] * ir[1]),
                2e-13,
            );
        }
    }
}

#[test]
fn discrete_storage_work_and_actual_endpoint_tangent_match_original_energy() {
    let p = parameters();
    let model = InductionMotor::new(p).unwrap();
    let storage = model.system().storage();
    let rule = storage.discrete_gradient_rule().unwrap();
    let a = model
        .state([1.2, -0.4], [-0.7, 1.3], 0.2, [2.0, -0.5])
        .unwrap();
    let direction = [0.13, -0.05, 0.06, 0.09, 0.17, 0.03, -0.02];
    for change in [0.0, 1e-12, -1e-8, 0.3, -0.7, std::f64::consts::PI] {
        let b = model
            .state([-0.7, 1.6], [1.5, -0.8], 0.2 + change, [3.0, -0.7])
            .unwrap();
        let mid: [f64; 7] = std::array::from_fn(|i| f64::midpoint(a[i], b[i]));
        let (mut dg, mut reversed) = ([0.0; 7], [0.0; 7]);
        rule.evaluate(&a, &b, &mid, &mut dg);
        rule.evaluate(&b, &a, &mid, &mut reversed);
        for i in 0..7 {
            close(dg[i], reversed[i], 3e-14);
        }
        let work: f64 = (0..7).map(|i| dg[i] * (b[i] - a[i])).sum();
        close(work, energy(p, &b) - energy(p, &a), 8e-14);
        let mut tangent = [0.0; 7];
        assert!(rule.end_jacobian_action(&a, &b, &direction, &mut tangent));
        let h = 1e-6;
        let plus: [f64; 7] = std::array::from_fn(|i| b[i] + h * direction[i]);
        let minus: [f64; 7] = std::array::from_fn(|i| b[i] - h * direction[i]);
        let (mut up, mut down) = ([0.0; 7], [0.0; 7]);
        rule.evaluate(
            &a,
            &plus,
            &std::array::from_fn::<_, 7, _>(|i| f64::midpoint(a[i], plus[i])),
            &mut up,
        );
        rule.evaluate(
            &a,
            &minus,
            &std::array::from_fn::<_, 7, _>(|i| f64::midpoint(a[i], minus[i])),
            &mut down,
        );
        for i in 0..7 {
            close(
                tangent[i],
                (up[i] - down[i]) / (2.0 * h),
                2e-8 * (1.0 + tangent[i].abs()),
            );
        }
        let mut midpoint_gradient = [0.0; 7];
        storage.gradient(&mid, &mut midpoint_gradient);
        let original_change = storage.midpoint_energy_remainder(&a, &b, &mid).unwrap()
            + (0..7)
                .map(|i| midpoint_gradient[i] * (b[i] - a[i]))
                .sum::<f64>();
        close(original_change, energy(p, &b) - energy(p, &a), 8e-14);
    }
    let mut tangent = [0.0; 7];
    assert!(
        storage
            .continuous_tangent()
            .unwrap()
            .apply(&a, &direction, &mut tangent)
    );
    let h = 1e-6;
    let (mut up, mut down) = ([0.0; 7], [0.0; 7]);
    storage.gradient(
        &std::array::from_fn::<_, 7, _>(|i| a[i] + h * direction[i]),
        &mut up,
    );
    storage.gradient(
        &std::array::from_fn::<_, 7, _>(|i| a[i] - h * direction[i]),
        &mut down,
    );
    for i in 0..7 {
        close(
            tangent[i],
            (up[i] - down[i]) / (2.0 * h),
            2e-8 * (1.0 + tangent[i].abs()),
        );
    }
}

// Independent steady-state cage phasor: 0 = (Rr+j slip Lr) Ir + j slip Lm Is.
// Is is real. Its air-gap power/rotor copper ratio equals the signed slip.
fn phasor(p: InductionMotorParameters, frequency: f64, speed: f64) -> ([f64; 2], [f64; 2], f64) {
    let current = 2.0;
    let slip = frequency - f64::from(p.pole_pairs) * speed;
    let lr = p.rotor_leakage_inductance_h + p.magnetizing_inductance_h;
    let ls = p.stator_leakage_inductance_h + p.magnetizing_inductance_h;
    let lm = p.magnetizing_inductance_h;
    let denominator = p.rotor_resistance_ohm.powi(2) + (slip * lr).powi(2);
    let rotor = [
        -slip * slip * lm * lr * current / denominator,
        -slip * lm * p.rotor_resistance_ohm * current / denominator,
    ];
    let voltage = [
        p.stator_resistance_ohm * current - frequency * lm * rotor[1],
        frequency * (ls * current + lm * rotor[0]),
    ];
    let torque = -f64::from(p.pole_pairs) * lm * rotor[1] * current;
    (rotor, voltage, torque)
}

#[test]
fn cage_slip_controls_standstill_motoring_synchronous_and_generating_losses() {
    let p = parameters();
    let model = InductionMotor::new(p).unwrap();
    let frequency = 40.0;
    for speed in [0.0, 15.0, 20.0, 25.0] {
        let (rotor, voltage, torque) = phasor(p, frequency, speed);
        let state = model.state([2.0, 0.0], rotor, 0.0, [speed, 0.0]).unwrap();
        let point = model.operating_point(&state).unwrap();
        close(point.electromagnetic_torque_nm, torque, 2e-14);
        close(
            point.rotor_copper_power_w,
            (frequency / f64::from(p.pole_pairs) - speed) * torque,
            2e-13,
        );
        close(
            2.0 * voltage[0],
            point.stator_copper_power_w + point.rotor_copper_power_w + torque * speed,
            2e-13,
        );
        if speed < 20.0 {
            assert!(torque > 0.0 && point.rotor_copper_power_w > 0.0);
        }
        if speed > 20.0 {
            assert!(torque < 0.0 && point.rotor_copper_power_w > 0.0);
        }
        if speed == 20.0 {
            close(torque, 0.0, 1e-14);
            close(point.rotor_copper_power_w, 0.0, 1e-25);
        }
        // The native flux equations reproduce the phasor time derivative in
        // each winding's OWN frame, including both speed and slip frequency.
        close(
            voltage[0] - p.stator_resistance_ohm * point.stator_current_alpha_beta_a[0],
            -frequency * state[1],
            2e-13,
        );
        close(
            voltage[1] - p.stator_resistance_ohm * point.stator_current_alpha_beta_a[1],
            frequency * state[0],
            2e-13,
        );
        let slip = frequency - f64::from(p.pole_pairs) * speed;
        close(
            -p.rotor_resistance_ohm * point.rotor_current_dq_a[0],
            -slip * state[3],
            2e-13,
        );
        close(
            -p.rotor_resistance_ohm * point.rotor_current_dq_a[1],
            slip * state[2],
            2e-13,
        );
    }
}

fn integrate_steady(p: InductionMotorParameters, speed: f64, count: usize) -> f64 {
    let (frequency, duration) = (40.0, 0.05);
    let (rotor, voltage, torque) = phasor(p, frequency, speed);
    let model = InductionMotor::new(p).unwrap();
    let mut state = model.state([2.0, 0.0], rotor, 0.0, [speed, 0.0]).unwrap();
    let initial = state;
    let mut prepared = model.prepare_controlled(budget(), 8, || false).unwrap();
    let dt = duration / count as f64;
    let mut heat = [0.0; 3];
    let mut work = 0.0;
    for k in 0..count {
        let v = rotate(frequency * (k as f64 + 0.5) * dt, voltage);
        let u = [v[0], v[1], -torque, torque];
        let (mut next, mut outputs) = ([0.0; 7], [0.0; 4]);
        let step = prepared
            .step_controlled(&state, &u, dt, &mut next, &mut outputs, || false)
            .unwrap();
        let internal =
            step.electromagnetic_angular_impulse_nm_s + step.bearing_angular_impulse_nm_s;
        close(next[5] - state[5], internal + dt * u[2], 3e-12);
        close(next[6] - state[6], -internal + dt * u[3], 3e-12);
        close(
            next[4] - state[4],
            step.shaft_angle_increment_rad[0] - step.shaft_angle_increment_rad[1],
            3e-12,
        );
        let loss = step.stator_copper_heat_j + step.rotor_copper_heat_j + step.bearing_heat_j;
        close(step.native.dissipated, loss, 3e-12);
        close(
            step.native.delta_h + loss,
            step.electrical_work_j.iter().sum::<f64>() + step.shaft_work_j.iter().sum::<f64>(),
            3e-12,
        );
        heat[0] += step.stator_copper_heat_j;
        heat[1] += step.rotor_copper_heat_j;
        heat[2] += step.bearing_heat_j;
        work += step.electrical_work_j.iter().sum::<f64>() + step.shaft_work_j.iter().sum::<f64>();
        state = next;
    }
    close(
        energy(p, &state) - energy(p, &initial) + heat.iter().sum::<f64>(),
        work,
        2e-10,
    );
    let expected_s = rotate(frequency * duration, [initial[0], initial[1]]);
    let expected_r = rotate(
        (frequency - f64::from(p.pole_pairs) * speed) * duration,
        [initial[2], initial[3]],
    );
    let expected = [
        expected_s[0],
        expected_s[1],
        expected_r[0],
        expected_r[1],
        speed * duration,
        initial[5],
        initial[6],
    ];
    let mut error = (0..7)
        .map(|i| (state[i] - expected[i]).abs() / (1.0 + expected[i].abs()))
        .fold(0.0_f64, f64::max);
    error = error.max((heat[0] - 4.0 * p.stator_resistance_ohm * duration).abs());
    error = error.max(
        (heat[1] - p.rotor_resistance_ohm * (rotor[0] * rotor[0] + rotor[1] * rotor[1]) * duration)
            .abs(),
    );
    assert_eq!(heat[2], 0.0);
    error
}

#[test]
fn native_flux_transients_converge_to_independent_rotating_cage_phasors() {
    for speed in [0.0, 15.0, 25.0] {
        let coarse = integrate_steady(parameters(), speed, 64);
        let fine = integrate_steady(parameters(), speed, 128);
        eprintln!(
            "cage speed={speed} rad/s: coarse error={coarse:.8e}, fine error={fine:.8e}, ratio={:.5}",
            coarse / fine
        );
        assert!(
            fine < 0.3 * coarse,
            "speed={speed}, coarse={coarse}, fine={fine}"
        );
        assert!(fine < 3e-5, "speed={speed}, fine={fine}");
    }
}

#[test]
fn lossless_induction_preserves_four_fluxes_total_energy_and_total_momentum() {
    let p = InductionMotorParameters {
        stator_resistance_ohm: 0.0,
        rotor_resistance_ohm: 0.0,
        ..parameters()
    };
    let model = InductionMotor::new(p).unwrap();
    let mut state = model.state([2.0, 0.0], [0.0, -3.0], 0.0, [0.0; 2]).unwrap();
    let initial = state;
    let mut prepared = model.prepare_controlled(budget(), 8, || false).unwrap();
    for _ in 0..50 {
        let (mut next, mut outputs) = ([0.0; 7], [0.0; 4]);
        let step = prepared
            .step_controlled(&state, &[0.0; 4], 0.001, &mut next, &mut outputs, || false)
            .unwrap();
        assert_eq!(next[..4], initial[..4]);
        assert_eq!(step.stator_copper_heat_j, 0.0);
        assert_eq!(step.rotor_copper_heat_j, 0.0);
        assert_eq!(step.bearing_heat_j, 0.0);
        close(next[5] + next[6], 0.0, 2e-13);
        state = next;
    }
    assert!(state[5] > 0.0 && state[6] < 0.0 && state[4] > 0.0);
    close(energy(p, &state), energy(p, &initial), 2e-11);
}

#[test]
fn zero_mutual_limit_is_four_independent_rl_windings() {
    let p = InductionMotorParameters {
        magnetizing_inductance_h: 0.0,
        ..parameters()
    };
    let model = InductionMotor::new(p).unwrap();
    let state = model
        .state([1.0, -2.0], [0.5, -0.7], 0.2, [3.0, 3.0])
        .unwrap();
    let mut prepared = model.prepare_controlled(budget(), 8, || false).unwrap();
    let (dt, u) = (0.003, [2.0, 1.0, 0.0, 0.0]);
    let (mut next, mut outputs) = ([0.0; 7], [0.0; 4]);
    let step = prepared
        .step_controlled(&state, &u, dt, &mut next, &mut outputs, || false)
        .unwrap();
    for i in 0..4 {
        let (l, r, v) = if i < 2 {
            (p.stator_leakage_inductance_h, p.stator_resistance_ohm, u[i])
        } else {
            (p.rotor_leakage_inductance_h, p.rotor_resistance_ohm, 0.0)
        };
        let rate = 0.5 * dt * r / l;
        close(
            next[i],
            ((1.0 - rate) * state[i] + dt * v) / (1.0 + rate),
            1e-14,
        );
    }
    close(next[4], state[4], 1e-14);
    close(next[5], state[5], 1e-14);
    close(next[6], state[6], 1e-14);
    assert_eq!(step.electromagnetic_angular_impulse_nm_s, 0.0);
    assert!(step.rotor_copper_heat_j > 0.0);
}

#[test]
fn finite_dc_bus_receives_reciprocal_induction_load_and_keeps_rotor_loss() {
    let capacitance = 0.2;
    let bus = PortHamiltonian::new(
        1,
        1,
        vec![0.0],
        vec![0.0],
        vec![1.0],
        Box::new(QuadraticStorage::new(vec![1.0 / capacitance], 1).unwrap()),
    )
    .unwrap();
    let p = parameters();
    let model = InductionMotor::new(p).unwrap();
    let machine_state = model
        .state([1.0, 0.0], [-0.5, 0.2], 0.0, [2.0, -1.0])
        .unwrap();
    let modulation = [0.6, -0.2];
    let system = interconnect_mapped(
        bus,
        model.into_system(),
        &[0],
        &[0, 1],
        &modulation,
        InterconnectionBudget {
            max_states: 8,
            max_total_ports: 5,
            max_dense_bytes: 1 << 20,
            max_coupling_multiply_terms: 1_000_000,
        },
    )
    .unwrap();
    assert_eq!((system.state_dim(), system.port_dim()), (8, 2));
    let mut state = [0.0; 8];
    state[0] = capacitance * 5.0;
    state[1..].copy_from_slice(&machine_state);
    let mut workspace = StepWorkspace::new(&system).unwrap();
    workspace.set_solve_budget(budget()).unwrap();
    workspace.set_backtracking_limit(8).unwrap();
    let (mut next, mut outputs) = ([0.0; 8], [0.0; 2]);
    let dt = 0.001;
    let step = workspace
        .step_into_controlled(
            &system,
            &state,
            &[0.0; 2],
            dt,
            &mut next,
            &mut outputs,
            || false,
        )
        .unwrap();
    let e = workspace.last_step_discrete_effort().unwrap();
    close(
        next[0] - state[0],
        -dt * (modulation[0] * e[1] + modulation[1] * e[2]),
        3e-13,
    );
    close(
        next[1] - state[1],
        dt * (modulation[0] * e[0] - p.stator_resistance_ohm * e[1]),
        3e-13,
    );
    close(
        next[2] - state[2],
        dt * (modulation[1] * e[0] - p.stator_resistance_ohm * e[2]),
        3e-13,
    );
    let copper = dt
        * (p.stator_resistance_ohm * (e[1] * e[1] + e[2] * e[2])
            + p.rotor_resistance_ohm * (e[3] * e[3] + e[4] * e[4]));
    close(step.dissipated, copper, 3e-13);
    close(step.delta_h + step.dissipated, 0.0, 3e-12);
    close(next[6] + next[7], state[6] + state[7], 3e-13);
}

#[test]
fn invalid_sources_and_late_cancellation_publish_nothing_and_retry() {
    let p = parameters();
    for bad in [
        InductionMotorParameters {
            stator_leakage_inductance_h: 0.0,
            ..p
        },
        InductionMotorParameters {
            rotor_leakage_inductance_h: -1.0,
            ..p
        },
        InductionMotorParameters {
            magnetizing_inductance_h: -1.0,
            ..p
        },
        InductionMotorParameters {
            rotor_resistance_ohm: -1.0,
            ..p
        },
        InductionMotorParameters { pole_pairs: 0, ..p },
        InductionMotorParameters {
            rotor_inertia_kg_m2: f64::NAN,
            ..p
        },
        InductionMotorParameters {
            stator_inertia_kg_m2: f64::INFINITY,
            ..p
        },
    ] {
        assert!(InductionMotor::new(bad).is_err());
    }
    let model = InductionMotor::new(p).unwrap();
    assert!(
        model
            .state([f64::NAN, 0.0], [0.0; 2], 0.0, [0.0; 2])
            .is_err()
    );
    let initial = model
        .state([1.0, -0.3], [-0.2, 0.6], 0.3, [3.0, -1.0])
        .unwrap();
    let mut prepared = model.prepare_controlled(budget(), 8, || false).unwrap();
    let (mut next, mut outputs) = ([42.0; 7], [43.0; 4]);
    let mut calls = 0;
    assert!(
        prepared
            .step_controlled(
                &initial,
                &[2.0, -1.0, 0.0, 0.0],
                0.002,
                &mut next,
                &mut outputs,
                || {
                    calls += 1;
                    false
                }
            )
            .is_ok()
    );
    let expected = (next, outputs);
    next = [42.0; 7];
    outputs = [43.0; 4];
    let mut trial_calls = 0;
    let result = prepared.step_controlled(
        &initial,
        &[2.0, -1.0, 0.0, 0.0],
        0.002,
        &mut next,
        &mut outputs,
        || {
            trial_calls += 1;
            trial_calls == calls
        },
    );
    assert_eq!(result, Err(PreparedStepError::Cancelled));
    assert_eq!((next, outputs), ([42.0; 7], [43.0; 4]));
    prepared
        .step_controlled(
            &initial,
            &[2.0, -1.0, 0.0, 0.0],
            0.002,
            &mut next,
            &mut outputs,
            || false,
        )
        .unwrap();
    assert_eq!((next, outputs), expected);
}
