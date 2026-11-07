//! Original bounded flux-law saturation, reciprocal torque and finite-bus work.
use fs_phs::circuit::inductance::InductorFluxCurve;
use fs_phs::circuit::synchronous::{
    SynchronousMotor, SynchronousMotorParameters,
    inverter::{AveragedSynchronousInverter, InverterParameters},
};
use fs_phs::{PreparedStepError, SolveBudget};

const D: [(f64, f64); 7] = [
    (-0.22, -50.0),
    (-0.10, -10.0),
    (-0.04, -2.0),
    (0.0, 0.0),
    (0.04, 2.0),
    (0.10, 10.0),
    (0.22, 50.0),
];
const Q: [(f64, f64); 7] = [
    (-0.30, -60.0),
    (-0.12, -12.0),
    (-0.04, -2.0),
    (0.0, 0.0),
    (0.04, 2.0),
    (0.12, 12.0),
    (0.30, 60.0),
];

fn parameters() -> SynchronousMotorParameters {
    SynchronousMotorParameters {
        // Deliberately unrelated: curve construction must replace these values.
        d_inductance_h: 2.0,
        q_inductance_h: 3.0,
        winding_resistance_ohm: 0.35,
        pm_flux_linkage_wb: 0.08,
        pole_pairs: 3,
        rotor_inertia_kg_m2: 0.025,
        stator_inertia_kg_m2: 0.07,
        bearing_drag_nm_s_per_rad: 0.002,
    }
}

fn model(p: SynchronousMotorParameters) -> SynchronousMotor {
    SynchronousMotor::with_flux_curves(
        p,
        InductorFluxCurve::try_new(D.to_vec()).unwrap(),
        InductorFluxCurve::try_new(Q.to_vec()).unwrap(),
    )
    .unwrap()
}

fn budget() -> SolveBudget {
    SolveBudget {
        relative: 2e-11,
        absolute: 2e-12,
        max_iterations: 30,
    }
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}; tolerance={tolerance:e}"
    );
}

// Independent source-segment integration and current interpolation. No storage,
// source-curve effort, energy or discrete-gradient callback defines the oracle.
fn current(knots: &[(f64, f64)], flux: f64) -> f64 {
    let pair = knots
        .windows(2)
        .find(|pair| pair[0].0 <= flux && flux <= pair[1].0)
        .unwrap();
    pair[0].1 + (flux - pair[0].0) * (pair[1].1 - pair[0].1) / (pair[1].0 - pair[0].0)
}

fn axis_energy(knots: &[(f64, f64)], flux: f64) -> f64 {
    let (lo, hi, sign) = if flux < 0.0 {
        (flux, 0.0, -1.0)
    } else {
        (0.0, flux, 1.0)
    };
    let mut work = 0.0;
    for pair in knots.windows(2) {
        let a = lo.max(pair[0].0);
        let b = hi.min(pair[1].0);
        if b > a {
            work += 0.5 * (current(knots, a) + current(knots, b)) * (b - a);
        }
    }
    sign * work
}

fn energy(p: SynchronousMotorParameters, x: &[f64; 5]) -> f64 {
    let (s, c) = (f64::from(p.pole_pairs) * x[2]).sin_cos();
    let d = c * x[0] + s * x[1] - p.pm_flux_linkage_wb;
    let q = -s * x[0] + c * x[1];
    axis_energy(&D, d)
        + axis_energy(&Q, q)
        + 0.5 * (x[3] * x[3] / p.rotor_inertia_kg_m2 + x[4] * x[4] / p.stator_inertia_kg_m2)
}

#[test]
fn original_curve_currents_and_angle_energy_own_saturated_torque() {
    for pm in [0.08, 0.0, -0.08] {
        let p = SynchronousMotorParameters {
            pm_flux_linkage_wb: pm,
            ..parameters()
        };
        let model = model(p);
        close(model.parameters().d_inductance_h, 0.02, 1e-16);
        close(model.parameters().q_inductance_h, 0.02, 1e-16);
        for currents in [[6.0, 14.0], [-6.0, 14.0], [6.0, -14.0]] {
            let state = model.state(currents, 0.23, [5.0, -2.0]).unwrap();
            let point = model.operating_point(&state).unwrap();
            for k in 0..2 {
                close(point.current_dq_a[k], currents[k], 8e-14);
            }
            close(point.differential_inductance_dq_h[0], 0.0075, 1e-16);
            close(point.differential_inductance_dq_h[1], 0.00375, 1e-16);
            let torque = f64::from(p.pole_pairs)
                * (point.flux_dq_wb[0] * currents[1] - point.flux_dq_wb[1] * currents[0]);
            close(point.electromagnetic_torque_nm, torque, 2e-13);
            close(
                point.speed_voltage_dq_v[0] * currents[0]
                    + point.speed_voltage_dq_v[1] * currents[1],
                torque * 7.0,
                2e-12,
            );
            close(model.system().hamiltonian(&state), energy(p, &state), 2e-13);
            let h = 1e-6;
            let mut plus = state;
            let mut minus = state;
            plus[2] += h;
            minus[2] -= h;
            close(
                -(energy(p, &plus) - energy(p, &minus)) / (2.0 * h),
                torque,
                2e-8,
            );
            let nominal_torque = f64::from(p.pole_pairs) * pm * currents[1];
            assert!(
                (nominal_torque - torque).abs() > 0.1,
                "nominal equal-L metadata must not erase nonlinear reluctance torque"
            );
        }
    }
}

#[test]
fn scalar_segment_work_gives_symmetric_machine_work_and_exact_endpoint_action() {
    let p = parameters();
    let model = model(p);
    let storage = model.system().storage();
    let rule = storage.discrete_gradient_rule().unwrap();
    let a = model.state([-12.0, 8.0], -0.17, [2.0, -0.5]).unwrap();
    for angle in [0.0, 1e-12, -1e-7, 0.07, 0.4] {
        let b = model
            .state([15.0, -8.0], a[2] + angle, [3.0, -0.7])
            .unwrap();
        let mid: [f64; 5] = std::array::from_fn(|i| f64::midpoint(a[i], b[i]));
        let (mut effort, mut reverse) = ([0.0; 5], [0.0; 5]);
        rule.evaluate(&a, &b, &mid, &mut effort);
        rule.evaluate(&b, &a, &mid, &mut reverse);
        for i in 0..5 {
            close(effort[i], reverse[i], 2e-12);
        }
        let work: f64 = (0..5).map(|i| effort[i] * (b[i] - a[i])).sum();
        close(work, energy(p, &b) - energy(p, &a), 3e-12);
        close(effort[3], mid[3] / p.rotor_inertia_kg_m2, 1e-14);
        close(effort[4], mid[4] / p.stator_inertia_kg_m2, 1e-14);
        let mut gradient = [0.0; 5];
        storage.gradient(&mid, &mut gradient);
        let mid_work: f64 = (0..5).map(|i| gradient[i] * (b[i] - a[i])).sum();
        close(
            storage.midpoint_energy_remainder(&a, &b, &mid).unwrap() + mid_work,
            energy(p, &b) - energy(p, &a),
            3e-12,
        );
        let direction = [0.13, -0.05, 0.17, 0.3, -0.2];
        let mut action = [0.0; 5];
        assert!(rule.end_jacobian_action(&a, &b, &direction, &mut action));
        let h = 1e-6;
        let plus: [f64; 5] = std::array::from_fn(|i| b[i] + h * direction[i]);
        let minus: [f64; 5] = std::array::from_fn(|i| b[i] - h * direction[i]);
        let (mut up, mut down) = ([0.0; 5], [0.0; 5]);
        rule.evaluate(&a, &plus, &mid, &mut up);
        rule.evaluate(&a, &minus, &mid, &mut down);
        for i in 0..5 {
            close(
                action[i],
                (up[i] - down[i]) / (2.0 * h),
                2e-7 * (1.0 + action[i].abs()),
            );
        }
    }
}

#[test]
fn continuous_tangent_is_original_hessian_and_refuses_actual_source_kinks() {
    let p = SynchronousMotorParameters {
        pm_flux_linkage_wb: 0.0,
        ..parameters()
    };
    let model = model(p);
    let state = model.state([6.0, 7.0], 0.0, [2.0, -0.5]).unwrap();
    let storage = model.system().storage();
    let tangent = storage.continuous_tangent().unwrap();
    let direction = [0.03, 0.02, -0.11, 0.1, -0.2];
    let mut actual = [0.0; 5];
    assert!(tangent.apply(&state, &direction, &mut actual));
    let h = 1e-6;
    let (mut up, mut down) = ([0.0; 5], [0.0; 5]);
    storage.gradient(
        &std::array::from_fn::<_, 5, _>(|i| state[i] + h * direction[i]),
        &mut up,
    );
    storage.gradient(
        &std::array::from_fn::<_, 5, _>(|i| state[i] - h * direction[i]),
        &mut down,
    );
    for i in 0..5 {
        close(
            actual[i],
            (up[i] - down[i]) / (2.0 * h),
            2e-7 * (1.0 + actual[i].abs()),
        );
    }
    for d_flux in [0.04, 0.22] {
        let mut edge = state;
        edge[0] = d_flux;
        assert!(!tangent.apply(&edge, &direction, &mut actual));
        assert!(
            storage
                .discrete_gradient_rule()
                .unwrap()
                .end_jacobian_action(&edge, &edge, &direction, &mut actual),
            "bounded Newton law keeps the original scalar one-sided convention"
        );
    }
}

#[test]
fn linear_source_curves_reproduce_the_optimized_linear_owner() {
    let p = SynchronousMotorParameters {
        d_inductance_h: 0.012,
        q_inductance_h: 0.019,
        ..parameters()
    };
    let linear = SynchronousMotor::new(p).unwrap();
    let curves = SynchronousMotor::with_flux_curves(
        p,
        InductorFluxCurve::try_new(vec![
            (-1.0, -1.0 / p.d_inductance_h),
            (0.0, 0.0),
            (1.0, 1.0 / p.d_inductance_h),
        ])
        .unwrap(),
        InductorFluxCurve::try_new(vec![
            (-1.0, -1.0 / p.q_inductance_h),
            (0.0, 0.0),
            (1.0, 1.0 / p.q_inductance_h),
        ])
        .unwrap(),
    )
    .unwrap();
    let initial = linear.state([-3.0, 4.0], 0.21, [3.0, -1.0]).unwrap();
    let after = linear.state([5.0, -2.0], 0.48, [2.0, -0.8]).unwrap();
    let middle: [f64; 5] = std::array::from_fn(|i| f64::midpoint(initial[i], after[i]));
    let (mut a, mut b) = ([0.0; 5], [0.0; 5]);
    linear
        .system()
        .storage()
        .discrete_gradient_rule()
        .unwrap()
        .evaluate(&initial, &after, &middle, &mut a);
    curves
        .system()
        .storage()
        .discrete_gradient_rule()
        .unwrap()
        .evaluate(&initial, &after, &middle, &mut b);
    for i in 0..5 {
        close(a[i], b[i], 2e-13);
    }
    let mut linear = linear.prepare_controlled(budget(), 8, || false).unwrap();
    let mut curves = curves.prepare_controlled(budget(), 8, || false).unwrap();
    let (mut xa, mut xb, mut ya, mut yb) = ([0.0; 5], [0.0; 5], [0.0; 4], [0.0; 4]);
    let inputs = [7.0, -3.0, 0.4, -0.1];
    linear
        .step_controlled(&initial, &inputs, 0.001, &mut xa, &mut ya, || false)
        .unwrap();
    curves
        .step_controlled(&initial, &inputs, 0.001, &mut xb, &mut yb, || false)
        .unwrap();
    for i in 0..5 {
        close(xa[i], xb[i], 5e-13);
    }
    for i in 0..4 {
        close(ya[i], yb[i], 5e-12);
    }
}

#[test]
fn original_source_support_covers_split_corners_and_refused_newton_is_retryable() {
    let p = SynchronousMotorParameters {
        pm_flux_linkage_wb: 0.0,
        ..parameters()
    };
    let model = model(p);
    assert!(model.state([51.0, 0.0], 0.0, [0.0; 2]).is_err());
    let a = [0.0, 0.28, 0.0, 0.0, 0.0];
    let b = [-0.28, 0.0, std::f64::consts::FRAC_PI_2 / 3.0, 0.0, 0.0];
    assert!(model.operating_point(&a).is_ok());
    assert!(model.operating_point(&b).is_ok());
    let middle: [f64; 5] = std::array::from_fn(|i| f64::midpoint(a[i], b[i]));
    let mut effort = [0.0; 5];
    let rule = model.system().storage().discrete_gradient_rule().unwrap();
    rule.evaluate(&a, &b, &middle, &mut effort);
    assert!(
        effort.iter().all(|v| v.is_nan()),
        "valid endpoints cannot license an unsupported cross corner"
    );
    assert!(!rule.end_jacobian_action(&a, &b, &[1.0; 5], &mut effort));
    let initial = model.state([45.0, 0.0], 0.0, [0.0; 2]).unwrap();
    let mut prepared = model.prepare_controlled(budget(), 8, || false).unwrap();
    let (mut next, mut outputs) = ([123.0; 5], [456.0; 4]);
    assert!(
        prepared
            .step_controlled(
                &initial,
                &[2000.0, 0.0, 0.0, 0.0],
                0.02,
                &mut next,
                &mut outputs,
                || false
            )
            .is_err()
    );
    assert_eq!(next, [123.0; 5]);
    assert_eq!(outputs, [456.0; 4]);
    assert!(matches!(
        prepared.step_controlled(&initial, &[0.0; 4], 1e-4, &mut next, &mut outputs, || true),
        Err(PreparedStepError::Cancelled)
    ));
    assert_eq!(next, [123.0; 5]);
    assert_eq!(outputs, [456.0; 4]);
    prepared
        .step_controlled(&initial, &[0.0; 4], 1e-4, &mut next, &mut outputs, || false)
        .unwrap();
    let mut fresh = self::model(p)
        .prepare_controlled(budget(), 8, || false)
        .unwrap();
    let (mut expected, mut expected_outputs) = ([0.0; 5], [0.0; 4]);
    fresh
        .step_controlled(
            &initial,
            &[0.0; 4],
            1e-4,
            &mut expected,
            &mut expected_outputs,
            || false,
        )
        .unwrap();
    assert_eq!(next, expected);
    assert_eq!(outputs, expected_outputs);
}

fn bridge() -> InverterParameters {
    InverterParameters {
        dc_capacitance_f: 0.03,
        high_on_resistance_ohm: [0.03, 0.04, 0.02],
        low_on_resistance_ohm: [0.02, 0.03, 0.04],
    }
}

// Independent continuous coordinates [Vdc, lambda_d, lambda_q, theta, wr, ws,
// source/shaft work, copper+bearing+channel heat, electromagnetic impulse].
fn rhs(y: [f64; 9], duty: [f64; 3], inputs: [f64; 3]) -> [f64; 9] {
    let p = parameters();
    let bridge = bridge();
    let id = current(&D, y[1] - p.pm_flux_linkage_wb);
    let iq = current(&Q, y[2]);
    let (s, c) = (f64::from(p.pole_pairs) * y[3]).sin_cos();
    let ia = c * id - s * iq;
    let ib = s * id + c * iq;
    let abc = [
        (2.0_f64 / 3.0).sqrt() * ia,
        -ia / 6.0_f64.sqrt() + ib / 2.0_f64.sqrt(),
        -ia / 6.0_f64.sqrt() - ib / 2.0_f64.sqrt(),
    ];
    let channel: [f64; 3] = std::array::from_fn(|k| {
        duty[k] * bridge.high_on_resistance_ohm[k]
            + (1.0 - duty[k]) * bridge.low_on_resistance_ohm[k]
    });
    let phase: [f64; 3] = std::array::from_fn(|k| duty[k] * y[0] - channel[k] * abc[k]);
    let va = (2.0_f64 / 3.0).sqrt() * (phase[0] - 0.5 * (phase[1] + phase[2]));
    let vb = (phase[1] - phase[2]) / 2.0_f64.sqrt();
    let vd = c * va + s * vb;
    let vq = -s * va + c * vb;
    let relative = y[4] - y[5];
    let electrical = f64::from(p.pole_pairs) * relative;
    let torque = f64::from(p.pole_pairs) * (y[1] * iq - y[2] * id);
    let bearing = p.bearing_drag_nm_s_per_rad * relative;
    let dc_current: f64 = duty.iter().zip(abc).map(|(d, i)| d * i).sum();
    let channel_heat: f64 = channel.iter().zip(abc).map(|(r, i)| r * i * i).sum();
    [
        (inputs[0] - dc_current) / bridge.dc_capacitance_f,
        vd - p.winding_resistance_ohm * id + electrical * y[2],
        vq - p.winding_resistance_ohm * iq - electrical * y[1],
        relative,
        (inputs[1] + torque - bearing) / p.rotor_inertia_kg_m2,
        (inputs[2] - torque + bearing) / p.stator_inertia_kg_m2,
        inputs[0] * y[0] + inputs[1] * y[4] + inputs[2] * y[5],
        p.winding_resistance_ohm * (id * id + iq * iq) + bearing * relative + channel_heat,
        torque,
    ]
}

#[test]
fn retained_saturation_in_finite_bus_inverter_refines_to_independent_flux_dynamics() {
    let p = parameters();
    let duty = [[0.8, 0.35, 0.15], [0.2, 0.65, 0.85]];
    let inputs = [[0.2, 0.05, -0.02], [-0.1, -0.03, 0.01]];
    let duration = 0.003;
    let initial_model = model(p);
    let machine_state = initial_model.state([6.0, 8.0], 0.13, [2.0, -0.5]).unwrap();
    let initial_point = initial_model.operating_point(&machine_state).unwrap();
    let requested_torque = initial_point.electromagnetic_torque_nm + 0.2;
    let current_reference = initial_model
        .current_reference_for_torque_controlled(requested_torque, 6.0, 20.0, || false)
        .unwrap();
    assert!(!current_reference.torque_limited);
    let mut reference = [
        12.0,
        initial_point.flux_dq_wb[0],
        initial_point.flux_dq_wb[1],
        0.13,
        2.0,
        -0.5,
        0.0,
        0.0,
        0.0,
    ];
    for stage in 0..2 {
        let h = duration / 1024.0;
        for _ in 0..1024 {
            let k1 = rhs(reference, duty[stage], inputs[stage]);
            let k2 = rhs(
                std::array::from_fn(|i| reference[i] + 0.5 * h * k1[i]),
                duty[stage],
                inputs[stage],
            );
            let k3 = rhs(
                std::array::from_fn(|i| reference[i] + 0.5 * h * k2[i]),
                duty[stage],
                inputs[stage],
            );
            let k4 = rhs(
                std::array::from_fn(|i| reference[i] + h * k3[i]),
                duty[stage],
                inputs[stage],
            );
            reference = std::array::from_fn(|i| {
                reference[i] + h * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0
            });
        }
    }
    let mut errors = [0.0; 3];
    for (refinement, count) in [4, 8, 16].into_iter().enumerate() {
        let inverter = AveragedSynchronousInverter::new(model(p), bridge()).unwrap();
        let mut state = inverter.state(12.0, &machine_state).unwrap();
        let initial = state;
        let mut prepared = inverter.prepare_controlled(budget(), 8, || false).unwrap();
        assert_eq!(
            prepared
                .current_reference_for_torque_controlled(requested_torque, 6.0, 20.0, || false)
                .unwrap(),
            current_reference,
            "prepared inverter must retain the exact source-owned torque allocator"
        );
        let point = prepared.machine_operating_point(&state).unwrap();
        for k in 0..2 {
            close(point.current_dq_a[k], [6.0, 8.0][k], 8e-14);
        }
        close(
            point.electromagnetic_torque_nm,
            initial_point.electromagnetic_torque_nm,
            2e-13,
        );
        assert!(
            point.differential_inductance_dq_h[0] < prepared.machine_parameters().d_inductance_h
        );
        let (mut work, mut heat, mut impulse) = (0.0, 0.0, 0.0);
        let mut regenerated = false;
        let dt = duration / f64::from(count);
        for stage in 0..2 {
            for _ in 0..count {
                let (mut next, mut outputs) = ([0.0; 6], [0.0; 3]);
                let step = prepared
                    .step_controlled(
                        &state,
                        duty[stage],
                        &inputs[stage],
                        dt,
                        &mut next,
                        &mut outputs,
                        || false,
                    )
                    .unwrap();
                let channels: f64 = step.channel_heat_j.iter().flatten().sum();
                work += step.source_work_j + step.shaft_work_j.iter().sum::<f64>();
                heat += step.copper_heat_j + step.bearing_heat_j + channels;
                impulse += step.electromagnetic_angular_impulse_nm_s;
                regenerated |= step.bridge_work_j < 0.0;
                close(
                    next[4] + next[5] - state[4] - state[5],
                    dt * (inputs[stage][1] + inputs[stage][2]),
                    3e-12,
                );
                state = next;
            }
        }
        assert!(
            regenerated,
            "held reverse vector must actually return winding energy to this finite bus"
        );
        let point = prepared.machine_operating_point(&state).unwrap();
        let actual = [
            state[0] / bridge().dc_capacitance_f,
            point.flux_dq_wb[0],
            point.flux_dq_wb[1],
            state[3],
            point.shaft_speed_rad_s[0],
            point.shaft_speed_rad_s[1],
            work,
            heat,
            impulse,
        ];
        let final_machine = [state[1], state[2], state[3], state[4], state[5]];
        let initial_energy =
            0.5 * initial[0] * initial[0] / bridge().dc_capacitance_f + energy(p, &machine_state);
        let final_energy =
            0.5 * state[0] * state[0] / bridge().dc_capacitance_f + energy(p, &final_machine);
        close(final_energy - initial_energy, work - heat, 2e-10);
        let scales = [12.0, 0.1, 0.1, 1.0, 2.0, 1.0, 0.01, 0.1, 0.01];
        errors[refinement] = actual
            .iter()
            .zip(reference)
            .zip(scales)
            .map(|((a, b), scale)| ((a - b) / scale).powi(2))
            .sum::<f64>()
            .sqrt();
    }
    assert!(
        errors[1] < 0.31 * errors[0] && errors[2] < 0.31 * errors[1],
        "second-order nonlinear trajectory/heat/work refinement: {errors:?}"
    );
    assert!(
        errors[2] < 2e-4,
        "refined physical trajectory error: {errors:?}"
    );
}
