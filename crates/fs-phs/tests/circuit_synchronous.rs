//! Native PM/saliency torque, reciprocal voltage and independent d/q dynamics.
use fs_phs::circuit::synchronous::{
    SynchronousMotor, SynchronousMotorParameters, phase_voltage_to_alpha_beta,
};
use fs_phs::{PreparedStepError, SolveBudget};

#[path = "circuit_synchronous/finite_bus.rs"]
mod finite_bus;

fn parameters() -> SynchronousMotorParameters {
    SynchronousMotorParameters {
        d_inductance_h: 0.012,
        q_inductance_h: 0.019,
        winding_resistance_ohm: 0.35,
        pm_flux_linkage_wb: 0.08,
        pole_pairs: 3,
        rotor_inertia_kg_m2: 0.025,
        stator_inertia_kg_m2: 0.07,
        bearing_drag_nm_s_per_rad: 0.002,
    }
}

fn budget() -> SolveBudget {
    SolveBudget {
        relative: 2e-11,
        absolute: 2e-12,
        max_iterations: 30,
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}; tolerance={tolerance:e}"
    );
}

// Independent original-coordinate storage; no production field/effort callback.
fn energy(p: SynchronousMotorParameters, x: &[f64; 5]) -> f64 {
    let (s, c) = (f64::from(p.pole_pairs) * x[2]).sin_cos();
    let d = c * x[0] + s * x[1] - p.pm_flux_linkage_wb;
    let q = -s * x[0] + c * x[1];
    0.5 * (d * d / p.d_inductance_h
        + q * q / p.q_inductance_h
        + x[3] * x[3] / p.rotor_inertia_kg_m2
        + x[4] * x[4] / p.stator_inertia_kg_m2)
}

#[test]
fn physical_phase_power_pm_and_reluctance_torque_share_one_normalization() {
    for (pm, ld, lq) in [
        (0.08, 0.012, 0.019),
        (0.0, 0.019, 0.012),
        (0.08, 0.012, 0.012),
    ] {
        for pole_pairs in [1, 3] {
            let p = SynchronousMotorParameters {
                pm_flux_linkage_wb: pm,
                d_inductance_h: ld,
                q_inductance_h: lq,
                pole_pairs,
                ..parameters()
            };
            let model = SynchronousMotor::new(p).unwrap();
            for current in [[-1.2, -2.3], [1.2, 2.3], [-1.2, 2.3], [1.2, -2.3]] {
                let theta = 0.37;
                let state = model.state(current, theta, [12.0, -4.0]).unwrap();
                let point = model.operating_point(&state).unwrap();
                let torque =
                    f64::from(pole_pairs) * (pm * current[1] + (ld - lq) * current[0] * current[1]);
                close(point.current_dq_a[0], current[0], 2e-13);
                close(point.current_dq_a[1], current[1], 2e-13);
                close(point.electromagnetic_torque_nm, torque, 2e-13);
                close(point.phase_current_a.iter().sum(), 0.0, 2e-14);
                close(
                    point
                        .speed_voltage_dq_v
                        .iter()
                        .zip(current)
                        .map(|(v, i)| v * i)
                        .sum(),
                    torque * 16.0,
                    2e-12,
                );
                let phase_voltage = [7.0, -2.0, 4.0];
                let ab = phase_voltage_to_alpha_beta(phase_voltage).unwrap();
                close(
                    ab.iter()
                        .zip(point.current_alpha_beta_a)
                        .map(|(v, i)| v * i)
                        .sum(),
                    phase_voltage
                        .iter()
                        .zip(point.phase_current_a)
                        .map(|(v, i)| v * i)
                        .sum(),
                    2e-13,
                );
                let dq_voltage = [3.0, -5.0];
                let voltage = model.voltage_from_dq(theta, dq_voltage).unwrap();
                close(
                    voltage
                        .iter()
                        .zip(point.current_alpha_beta_a)
                        .map(|(v, i)| v * i)
                        .sum(),
                    dq_voltage.iter().zip(current).map(|(v, i)| v * i).sum(),
                    2e-13,
                );
                // Only relative speed generates voltage, even with a spinning housing.
                let co_rotating = model.state(current, theta, [12.0, 12.0]).unwrap();
                for value in model
                    .operating_point(&co_rotating)
                    .unwrap()
                    .speed_voltage_dq_v
                {
                    close(value, 0.0, 2e-13);
                }
                // Negative energy derivative at fixed ORIGINAL flux is the torque.
                let h = 2e-6;
                let mut before = state;
                let mut after = state;
                before[2] -= h;
                after[2] += h;
                close(
                    -(energy(p, &after) - energy(p, &before)) / (2.0 * h),
                    torque,
                    2e-9,
                );
            }
        }
    }
    let model = SynchronousMotor::new(parameters()).unwrap();
    let state = model.state([0.0; 2], 0.0, [5.0, -2.0]).unwrap();
    let point = model.operating_point(&state).unwrap();
    close(point.speed_voltage_dq_v[1], 3.0 * 7.0 * 0.08, 1e-14);
    assert_eq!(point.electromagnetic_torque_nm, 0.0);
    assert!(state[0] > 0.0, "zero current retains PM flux linkage");
}

#[test]
fn source_discrete_work_is_symmetric_exact_and_has_its_actual_analytic_tangent() {
    let p = parameters();
    let model = SynchronousMotor::new(p).unwrap();
    let storage = model.system().storage();
    let rule = storage.discrete_gradient_rule().unwrap();
    let a = model.state([-0.7, 1.3], 0.2, [2.0, -0.5]).unwrap();
    for change in [0.0, 1e-12, -1e-8, 0.3, -0.7, std::f64::consts::TAU / 3.0] {
        let b = model.state([1.5, -0.8], 0.2 + change, [3.0, -0.7]).unwrap();
        let mid: [f64; 5] = std::array::from_fn(|i| f64::midpoint(a[i], b[i]));
        let (mut dg, mut reverse) = ([0.0; 5], [0.0; 5]);
        rule.evaluate(&a, &b, &mid, &mut dg);
        rule.evaluate(&b, &a, &mid, &mut reverse);
        for i in 0..5 {
            close(dg[i], reverse[i], 1e-13);
        }
        let work: f64 = (0..5).map(|i| dg[i] * (b[i] - a[i])).sum();
        close(work, energy(p, &b) - energy(p, &a), 3e-13);
        close(dg[3], mid[3] / p.rotor_inertia_kg_m2, 1e-14);
        close(dg[4], mid[4] / p.stator_inertia_kg_m2, 1e-14);
        let direction = [0.13, -0.05, 0.17, 0.3, -0.2];
        let mut analytic = [0.0; 5];
        assert!(rule.end_jacobian_action(&a, &b, &direction, &mut analytic));
        let h = 1e-6;
        let plus: [f64; 5] = std::array::from_fn(|i| b[i] + h * direction[i]);
        let minus: [f64; 5] = std::array::from_fn(|i| b[i] - h * direction[i]);
        let (mut up, mut down) = ([0.0; 5], [0.0; 5]);
        rule.evaluate(
            &a,
            &plus,
            &std::array::from_fn::<_, 5, _>(|i| f64::midpoint(a[i], plus[i])),
            &mut up,
        );
        rule.evaluate(
            &a,
            &minus,
            &std::array::from_fn::<_, 5, _>(|i| f64::midpoint(a[i], minus[i])),
            &mut down,
        );
        for i in 0..5 {
            close(
                analytic[i],
                (up[i] - down[i]) / (2.0 * h),
                2e-8 * (1.0 + analytic[i].abs()),
            );
        }
    }
    let direction = [0.03, 0.02, -0.11, 0.1, -0.2];
    let mut tangent = [0.0; 5];
    assert!(
        storage
            .continuous_tangent()
            .unwrap()
            .apply(&a, &direction, &mut tangent)
    );
    let h = 1e-6;
    let (mut up, mut down) = ([0.0; 5], [0.0; 5]);
    storage.gradient(
        &std::array::from_fn::<_, 5, _>(|i| a[i] + h * direction[i]),
        &mut up,
    );
    storage.gradient(
        &std::array::from_fn::<_, 5, _>(|i| a[i] - h * direction[i]),
        &mut down,
    );
    for i in 0..5 {
        close(
            tangent[i],
            (up[i] - down[i]) / (2.0 * h),
            2e-8 * (1.0 + tangent[i].abs()),
        );
    }
}

#[test]
fn zero_magnet_nonsalient_limit_matches_two_rl_windings_without_housing_drag() {
    let p = SynchronousMotorParameters {
        d_inductance_h: 0.04,
        q_inductance_h: 0.04,
        pm_flux_linkage_wb: 0.0,
        winding_resistance_ohm: 0.8,
        rotor_inertia_kg_m2: 1.0,
        stator_inertia_kg_m2: 2.0,
        ..parameters()
    };
    let model = SynchronousMotor::new(p).unwrap();
    let state = model.state([1.0, -2.0], 0.2, [2.0, 2.0]).unwrap();
    let current = model.operating_point(&state).unwrap().current_alpha_beta_a;
    let mut prepared = model.prepare_controlled(budget(), 8, || false).unwrap();
    let (mut next, mut output) = ([0.0; 5], [0.0; 4]);
    let (dt, voltage) = (0.015, [1.0, -2.0]);
    let step = prepared
        .step_controlled(
            &state,
            &[voltage[0], voltage[1], 0.0, 0.0],
            dt,
            &mut next,
            &mut output,
            || false,
        )
        .unwrap();
    let rate = p.winding_resistance_ohm * dt / (2.0 * p.d_inductance_h);
    let mut heat = 0.0;
    for i in 0..2 {
        let expected =
            ((1.0 - rate) * current[i] + dt * voltage[i] / p.d_inductance_h) / (1.0 + rate);
        close(next[i] / p.d_inductance_h, expected, 2e-11);
        close(output[i], 0.5 * (current[i] + expected), 2e-11);
        heat += dt * p.winding_resistance_ohm * output[i] * output[i];
    }
    for i in 2..5 {
        close(next[i], state[i], 2e-13);
    }
    close(step.copper_heat_j, heat, 2e-13);
    assert_eq!(step.bearing_heat_j, 0.0);
    assert_eq!(step.electromagnetic_angular_impulse_nm_s, 0.0);
    close(
        step.native.delta_h + step.copper_heat_j,
        step.electrical_work_j.iter().sum(),
        2e-12,
    );
}

// Independent physical d/q-current equations, not production flux gradients.
// Last five coordinates integrate copper/bearing heat and electrical/shaft work.
fn rhs(p: SynchronousMotorParameters, u: [f64; 4], z: &[f64; 10]) -> [f64; 10] {
    let (s, c) = (f64::from(p.pole_pairs) * z[2]).sin_cos();
    let vd = c * u[0] + s * u[1];
    let vq = -s * u[0] + c * u[1];
    let slip = z[3] - z[4];
    let electrical_speed = f64::from(p.pole_pairs) * slip;
    let ld = p.d_inductance_h * z[0] + p.pm_flux_linkage_wb;
    let lq = p.q_inductance_h * z[1];
    let torque = f64::from(p.pole_pairs)
        * (p.pm_flux_linkage_wb * z[1] + (p.d_inductance_h - p.q_inductance_h) * z[0] * z[1]);
    let drag = p.bearing_drag_nm_s_per_rad * slip;
    [
        (vd - p.winding_resistance_ohm * z[0] + electrical_speed * lq) / p.d_inductance_h,
        (vq - p.winding_resistance_ohm * z[1] - electrical_speed * ld) / p.q_inductance_h,
        slip,
        (torque - drag + u[2]) / p.rotor_inertia_kg_m2,
        (-torque + drag + u[3]) / p.stator_inertia_kg_m2,
        p.winding_resistance_ohm * (z[0] * z[0] + z[1] * z[1]),
        p.bearing_drag_nm_s_per_rad * slip * slip,
        vd * z[0] + vq * z[1],
        u[2] * z[3],
        u[3] * z[4],
    ]
}

fn reference(
    p: SynchronousMotorParameters,
    u: [f64; 4],
    initial: [f64; 10],
    duration: f64,
) -> [f64; 10] {
    let mut z = initial;
    let dt = duration / 4096.0;
    for _ in 0..4096 {
        let a = rhs(p, u, &z);
        let b = rhs(p, u, &std::array::from_fn(|i| z[i] + 0.5 * dt * a[i]));
        let c = rhs(p, u, &std::array::from_fn(|i| z[i] + 0.5 * dt * b[i]));
        let d = rhs(p, u, &std::array::from_fn(|i| z[i] + dt * c[i]));
        for i in 0..10 {
            z[i] += dt / 6.0 * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]);
        }
    }
    z
}

fn trajectory(
    p: SynchronousMotorParameters,
    u: [f64; 4],
    z: [f64; 10],
    duration: f64,
    count: usize,
) -> [f64; 10] {
    let model = SynchronousMotor::new(p).unwrap();
    let mut state = model.state([z[0], z[1]], z[2], [z[3], z[4]]).unwrap();
    let initial_energy = energy(p, &state);
    let mut prepared = model.prepare_controlled(budget(), 8, || false).unwrap();
    let dt = duration / count as f64;
    let mut sum = [0.0; 5];
    for _ in 0..count {
        let (mut next, mut output) = ([0.0; 5], [0.0; 4]);
        let step = prepared
            .step_controlled(&state, &u, dt, &mut next, &mut output, || false)
            .unwrap();
        let internal =
            step.electromagnetic_angular_impulse_nm_s + step.bearing_angular_impulse_nm_s;
        close(next[3] - state[3], internal + dt * u[2], 2e-11);
        close(next[4] - state[4], -internal + dt * u[3], 2e-11);
        close(
            next[2] - state[2],
            step.shaft_angle_increment_rad[0] - step.shaft_angle_increment_rad[1],
            2e-11,
        );
        close(
            step.native.dissipated,
            step.copper_heat_j + step.bearing_heat_j,
            2e-11,
        );
        sum[0] += step.copper_heat_j;
        sum[1] += step.bearing_heat_j;
        sum[2] += step.electrical_work_j.iter().sum::<f64>();
        sum[3] += step.shaft_work_j[0];
        sum[4] += step.shaft_work_j[1];
        state = next;
    }
    close(
        energy(p, &state) - initial_energy + sum[0] + sum[1] - sum[2] - sum[3] - sum[4],
        0.0,
        3e-9,
    );
    close(
        state[3] + state[4],
        p.rotor_inertia_kg_m2 * z[3] + p.stator_inertia_kg_m2 * z[4] + duration * (u[2] + u[3]),
        2e-10,
    );
    let point = prepared.operating_point(&state).unwrap();
    [
        point.current_dq_a[0],
        point.current_dq_a[1],
        state[2],
        point.shaft_speed_rad_s[0],
        point.shaft_speed_rad_s[1],
        sum[0],
        sum[1],
        sum[2],
        sum[3],
        sum[4],
    ]
}

#[test]
fn pmsm_and_reluctance_dynamics_converge_to_independent_dq_ode_with_real_heat() {
    let initial = [-1.2, 2.3, 0.24, 12.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let input = [7.0, -3.0, -0.08, 0.025];
    let duration = 0.06;
    for pm in [0.08, 0.0] {
        let p = SynchronousMotorParameters {
            pm_flux_linkage_wb: pm,
            ..parameters()
        };
        let oracle = reference(p, input, initial, duration);
        let mut errors = Vec::new();
        for count in [60, 120, 240] {
            let actual = trajectory(p, input, initial, duration, count);
            errors.push(
                (0..10)
                    .map(|i| (actual[i] - oracle[i]).abs() / (1.0 + oracle[i].abs()))
                    .fold(0.0_f64, f64::max),
            );
            assert!(actual[5] > 0.0 && actual[6] > 0.0);
        }
        assert!(
            errors[1] < 0.35 * errors[0] && errors[2] < 0.35 * errors[1],
            "PM={pm}, errors={errors:?}"
        );
        assert!(errors[2] < 2e-4, "PM={pm}, fine error={}", errors[2]);
    }
}

#[test]
fn lossless_short_circuit_exchanges_field_energy_and_opposite_shaft_impulses() {
    let p = SynchronousMotorParameters {
        winding_resistance_ohm: 0.0,
        bearing_drag_nm_s_per_rad: 0.0,
        ..parameters()
    };
    let model = SynchronousMotor::new(p).unwrap();
    let mut state = model.state([0.0, 4.0], 0.0, [0.0; 2]).unwrap();
    let initial = state;
    let stored = energy(p, &state);
    let mut prepared = model.prepare_controlled(budget(), 8, || false).unwrap();
    for _ in 0..80 {
        let (mut next, mut output) = ([0.0; 5], [0.0; 4]);
        let step = prepared
            .step_controlled(&state, &[0.0; 4], 0.001, &mut next, &mut output, || false)
            .unwrap();
        assert_eq!(next[..2], initial[..2]);
        assert_eq!(step.copper_heat_j, 0.0);
        assert_eq!(step.bearing_heat_j, 0.0);
        close(next[3] + next[4], 0.0, 1e-12);
        state = next;
    }
    assert!(state[3] > 0.0 && state[4] < 0.0 && state[2] > 0.0);
    close(energy(p, &state), stored, 2e-10);
}

#[test]
fn source_refusals_and_late_cancellation_preserve_outputs_and_retry() {
    let p = parameters();
    for bad in [
        SynchronousMotorParameters {
            d_inductance_h: 0.0,
            ..p
        },
        SynchronousMotorParameters {
            q_inductance_h: -1.0,
            ..p
        },
        SynchronousMotorParameters { pole_pairs: 0, ..p },
        SynchronousMotorParameters {
            stator_inertia_kg_m2: f64::INFINITY,
            ..p
        },
        SynchronousMotorParameters {
            pm_flux_linkage_wb: f64::NAN,
            ..p
        },
        SynchronousMotorParameters {
            winding_resistance_ohm: -1.0,
            ..p
        },
        SynchronousMotorParameters {
            bearing_drag_nm_s_per_rad: -1.0,
            ..p
        },
    ] {
        assert!(SynchronousMotor::new(bad).is_err());
    }
    assert!(phase_voltage_to_alpha_beta([0.0, f64::NAN, 0.0]).is_err());
    let model = SynchronousMotor::new(p).unwrap();
    assert!(model.state([0.0; 2], f64::NAN, [0.0; 2]).is_err());
    let state = model.state([1.0, -2.0], 0.3, [2.0, -1.0]).unwrap();
    let mut prepared = model.prepare_controlled(budget(), 8, || false).unwrap();
    let (mut next, mut output) = ([0.0; 5], [0.0; 4]);
    let input = [2.0, 3.0, -0.1, 0.2];
    let mut polls = 0;
    prepared
        .step_controlled(&state, &input, 0.002, &mut next, &mut output, || {
            polls += 1;
            false
        })
        .unwrap();
    let expected = (next, output);
    let sentinel = ([17.0; 5], [23.0; 4]);
    (next, output) = sentinel;
    let mut count = 0;
    assert!(matches!(
        prepared.step_controlled(&state, &input, 0.002, &mut next, &mut output, || {
            count += 1;
            count == polls
        }),
        Err(PreparedStepError::Cancelled)
    ));
    assert_eq!((next, output), sentinel);
    assert!(
        prepared
            .step_controlled(&state, &input, -0.1, &mut next, &mut output, || false)
            .is_err()
    );
    assert_eq!((next, output), sentinel);
    prepared
        .step_controlled(&state, &input, 0.002, &mut next, &mut output, || false)
        .unwrap();
    assert_eq!((next, output), expected);
    let mut limited = SynchronousMotor::new(p)
        .unwrap()
        .prepare_controlled(
            SolveBudget {
                max_iterations: 0,
                ..budget()
            },
            0,
            || false,
        )
        .unwrap();
    (next, output) = sentinel;
    assert!(
        limited
            .step_controlled(&state, &input, 0.002, &mut next, &mut output, || false)
            .is_err()
    );
    assert_eq!((next, output), sentinel);
}
