//! G0/G1: original magnetizing energy, reciprocal torque and nonlinear dynamics.
use fs_phs::circuit::{
    inductance::InductorFluxCurve,
    induction::{InductionMotor, InductionMotorParameters},
};
use fs_phs::{PreparedStepError, SolveBudget};

const POSITIVE: [(f64, f64); 5] = [
    (0.0, 0.0),
    (0.05, 1.0),
    (0.10, 3.0),
    (0.15, 8.0),
    (0.25, 25.0),
];

fn curve() -> InductorFluxCurve {
    InductorFluxCurve::try_new(
        POSITIVE[1..]
            .iter()
            .rev()
            .map(|&(f, i)| (-f, -i))
            .chain(POSITIVE)
            .collect(),
    )
    .unwrap()
}
fn parameters() -> InductionMotorParameters {
    InductionMotorParameters {
        stator_leakage_inductance_h: 0.004,
        rotor_leakage_inductance_h: 0.006,
        // Deliberately unrelated: source owns the actual constitutive law.
        magnetizing_inductance_h: 0.9,
        stator_resistance_ohm: 0.8,
        rotor_resistance_ohm: 0.5,
        pole_pairs: 2,
        rotor_inertia_kg_m2: 0.1,
        stator_inertia_kg_m2: 0.2,
        bearing_drag_nm_s_per_rad: 0.01,
    }
}
fn model(p: InductionMotorParameters) -> InductionMotor {
    InductionMotor::with_magnetizing_curve(p, curve()).unwrap()
}
fn budget() -> SolveBudget {
    SolveBudget {
        relative: 2e-12,
        absolute: 2e-13,
        max_iterations: 30,
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}; tolerance={tolerance:e}"
    );
}
fn rotate(angle: f64, x: [f64; 2]) -> [f64; 2] {
    let (s, c) = angle.sin_cos();
    [c * x[0] - s * x[1], s * x[0] + c * x[1]]
}
fn cross(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}
fn square(x: [f64; 2]) -> f64 {
    x[0] * x[0] + x[1] * x[1]
}

// Independent interpolation/integration of the ORIGINAL psi_m,I_m source.
// Recovery bisects the original constitutive equation instead of using the
// production transformed knots, gradient, energy or discrete work callbacks.
fn source_current(flux: f64) -> f64 {
    let pair = POSITIVE
        .windows(2)
        .find(|p| flux >= p[0].0 && flux <= p[1].0)
        .unwrap();
    pair[0].1 + (flux - pair[0].0) * (pair[1].1 - pair[0].1) / (pair[1].0 - pair[0].0)
}
fn source_energy(flux: f64) -> f64 {
    POSITIVE
        .windows(2)
        .filter(|p| flux > p[0].0)
        .map(|p| {
            let end = flux.min(p[1].0);
            0.5 * (p[0].1 + source_current(end)) * (end - p[0].0)
        })
        .sum()
}
struct Oracle {
    is: [f64; 2],
    ir: [f64; 2],
    torque: f64,
    energy: f64,
}
fn oracle(p: InductionMotorParameters, x: &[f64; 7]) -> Oracle {
    let ls = p.stator_leakage_inductance_h;
    let lr = p.rotor_leakage_inductance_h;
    let angle = f64::from(p.pole_pairs) * x[4];
    let r = rotate(angle, [x[2], x[3]]);
    let y = [
        (lr * x[0] + ls * r[0]) / (ls + lr),
        (lr * x[1] + ls * r[1]) / (ls + lr),
    ];
    let magnitude = y[0].hypot(y[1]);
    let ell = ls * lr / (ls + lr);
    assert!(magnitude <= 0.25 + ell * 25.0 + 1e-14);
    let (mut lo, mut hi) = (0.0, 0.25);
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        if mid + ell * source_current(mid) < magnitude {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let flux = if magnitude == 0.0 {
        0.0
    } else {
        0.5 * (lo + hi)
    };
    let psi = if magnitude == 0.0 {
        [0.0; 2]
    } else {
        y.map(|v| flux * v / magnitude)
    };
    let is = [(x[0] - psi[0]) / ls, (x[1] - psi[1]) / ls];
    let ir_ab = [(r[0] - psi[0]) / lr, (r[1] - psi[1]) / lr];
    let ir = rotate(-angle, ir_ab);
    Oracle {
        is,
        ir,
        torque: f64::from(p.pole_pairs) * cross(ir_ab, psi),
        energy: 0.5 * ls * square(is)
            + 0.5 * lr * square(ir)
            + source_energy(flux)
            + 0.5 * (x[5] * x[5] / p.rotor_inertia_kg_m2 + x[6] * x[6] / p.stator_inertia_kg_m2),
    }
}

#[test]
fn original_radial_law_owns_current_energy_torque_and_rotation_equivariance() {
    let p = parameters();
    let model = model(p);
    assert_eq!(model.magnetizing_curve(), Some(&curve()));
    close(model.parameters().magnetizing_inductance_h, 0.05, 1e-16);
    for (is, ir) in [
        ([0.0, 0.0], [0.0, 0.0]),
        ([1.2, -0.4], [-0.5, 0.8]),
        ([6.0, 2.0], [-1.0, 3.0]),
        ([-10.0, 8.0], [2.0, -4.0]),
    ] {
        let x = model.state(is, ir, 0.21, [5.0, -2.0]).unwrap();
        let point = model.operating_point(&x).unwrap();
        let expected = oracle(p, &x);
        for k in 0..2 {
            close(point.stator_current_alpha_beta_a[k], is[k], 2e-12);
            close(point.rotor_current_dq_a[k], ir[k], 2e-12);
            close(expected.is[k], is[k], 2e-12);
            close(expected.ir[k], ir[k], 2e-12);
        }
        close(point.electromagnetic_torque_nm, expected.torque, 3e-12);
        close(model.system().hamiltonian(&x), expected.energy, 2e-12);
        let h = 1e-6;
        let (mut plus, mut minus) = (x, x);
        plus[4] += h;
        minus[4] -= h;
        close(
            -(oracle(p, &plus).energy - oracle(p, &minus).energy) / (2.0 * h),
            expected.torque,
            3e-8,
        );
        let rotated = model
            .state(rotate(0.73, is), rotate(0.73, ir), 0.21, [5.0, -2.0])
            .unwrap();
        close(model.system().hamiltonian(&rotated), expected.energy, 2e-12);
        close(
            model
                .operating_point(&rotated)
                .unwrap()
                .electromagnetic_torque_nm,
            expected.torque,
            3e-12,
        );
    }
}

#[test]
fn symmetric_radial_work_and_analytic_actions_cover_knots_zero_and_equal_radii() {
    let p = parameters();
    let model = model(p);
    let storage = model.system().storage();
    let rule = storage.discrete_gradient_rule().unwrap();
    let direction = [0.13, -0.05, 0.06, 0.09, 0.17, 0.03, -0.02];
    let a = model
        .state([1.2, -0.4], [-0.7, 1.3], 0.2, [2.0, -0.5])
        .unwrap();
    let zero = model.state([0.0; 2], [0.0; 2], 0.2, [2.0, -0.5]).unwrap();
    let opposite = [-a[0], -a[1], -a[2], -a[3], a[4], a[5], a[6]];
    let mut pairs = vec![(zero, zero), (zero, a), (a, zero), (a, opposite), (a, a)];
    for change in [0.0, 1e-12, -1e-8, 0.3, -0.7, std::f64::consts::PI] {
        let b = model
            .state([-4.0, 6.0], [2.0, -1.0], 0.2 + change, [3.0, -0.7])
            .unwrap();
        pairs.push((a, b));
    }
    for (a, b) in pairs {
        let mid: [f64; 7] = std::array::from_fn(|i| f64::midpoint(a[i], b[i]));
        let (mut effort, mut reversed, mut gradient) = ([0.0; 7], [0.0; 7], [0.0; 7]);
        rule.evaluate(&a, &b, &mid, &mut effort);
        rule.evaluate(&b, &a, &mid, &mut reversed);
        for k in 0..7 {
            close(effort[k], reversed[k], 2e-12);
        }
        let work: f64 = (0..7).map(|k| effort[k] * (b[k] - a[k])).sum();
        let change = oracle(p, &b).energy - oracle(p, &a).energy;
        close(work, change, 3e-12);
        storage.gradient(&mid, &mut gradient);
        close(
            storage.midpoint_energy_remainder(&a, &b, &mid).unwrap()
                + (0..7).map(|k| gradient[k] * (b[k] - a[k])).sum::<f64>(),
            change,
            3e-12,
        );
        let mut action = [0.0; 7];
        assert!(rule.end_jacobian_action(&a, &b, &direction, &mut action));
        let h = 1e-6;
        let plus: [f64; 7] = std::array::from_fn(|k| b[k] + h * direction[k]);
        let minus: [f64; 7] = std::array::from_fn(|k| b[k] - h * direction[k]);
        let (mut up, mut down) = ([0.0; 7], [0.0; 7]);
        rule.evaluate(&a, &plus, &mid, &mut up);
        rule.evaluate(&a, &minus, &mid, &mut down);
        for k in 0..7 {
            close(
                action[k],
                (up[k] - down[k]) / (2.0 * h),
                2e-7 * (1.0 + action[k].abs()),
            );
        }
    }

    // Tiny simultaneous flux/angle increments must retain work beside the
    // finite magnetic bias; subtracting two stored energies would lose it.
    for h in [1e-6, 1e-10, 1e-12] {
        let b: [f64; 7] = std::array::from_fn(|k| a[k] + h * direction[k]);
        let mid: [f64; 7] = std::array::from_fn(|k| f64::midpoint(a[k], b[k]));
        let (mut effort, mut gradient) = ([0.0; 7], [0.0; 7]);
        rule.evaluate(&a, &b, &mid, &mut effort);
        storage.gradient(&mid, &mut gradient);
        let work: f64 = (0..7).map(|k| effort[k] * (b[k] - a[k])).sum();
        let scale: f64 = (0..7).map(|k| (effort[k] * (b[k] - a[k])).abs()).sum();
        let original = storage.midpoint_energy_remainder(&a, &b, &mid).unwrap()
            + (0..7).map(|k| gradient[k] * (b[k] - a[k])).sum::<f64>();
        close(original, work, 256.0 * f64::EPSILON * scale);
    }
}

#[test]
fn continuous_tangent_is_the_original_hessian_and_refuses_source_kinks() {
    let p = parameters();
    let model = model(p);
    let x = model
        .state([6.0, 2.0], [-1.0, 3.0], 0.21, [5.0, -2.0])
        .unwrap();
    let storage = model.system().storage();
    let direction = [0.13, -0.05, 0.06, 0.09, 0.17, 0.03, -0.02];
    let mut action = [0.0; 7];
    assert!(
        storage
            .continuous_tangent()
            .unwrap()
            .apply(&x, &direction, &mut action)
    );
    let h = 1e-6;
    let plus: [f64; 7] = std::array::from_fn(|k| x[k] + h * direction[k]);
    let minus: [f64; 7] = std::array::from_fn(|k| x[k] - h * direction[k]);
    let (mut up, mut down) = ([0.0; 7], [0.0; 7]);
    storage.gradient(&plus, &mut up);
    storage.gradient(&minus, &mut down);
    for k in 0..7 {
        close(
            action[k],
            (up[k] - down[k]) / (2.0 * h),
            2e-7 * (1.0 + action[k].abs()),
        );
    }

    let exact = InductorFluxCurve::try_new(vec![
        (-0.5, -9.0),
        (-0.25, -3.0),
        (-0.125, -1.0),
        (0.0, 0.0),
        (0.125, 1.0),
        (0.25, 3.0),
        (0.5, 9.0),
    ])
    .unwrap();
    let exact = InductionMotor::with_magnetizing_curve(
        InductionMotorParameters {
            stator_leakage_inductance_h: 0.125,
            rotor_leakage_inductance_h: 0.125,
            ..p
        },
        exact,
    )
    .unwrap();
    let tangent = exact.system().storage().continuous_tangent().unwrap();
    // y=psi+ell*I=.25+.0625*3=.4375, exactly represented.
    let mut kink = [0.4375, 0.0, 0.4375, 0.0, 0.0, 0.0, 0.0];
    assert!(!tangent.apply(&kink, &direction, &mut action));
    kink[0] = 1.0625;
    kink[2] = 1.0625;
    assert!(!tangent.apply(&kink, &direction, &mut action));
    assert!(tangent.apply(&[0.0; 7], &direction, &mut action));
}

#[test]
fn linear_source_matches_the_original_flux_machine_and_prepared_step() {
    let p = InductionMotorParameters {
        magnetizing_inductance_h: 0.05,
        ..parameters()
    };
    let linear = InductionMotor::new(p).unwrap();
    let saturated = InductionMotor::with_magnetizing_curve(
        p,
        InductorFluxCurve::try_new(vec![(-2.0, -40.0), (0.0, 0.0), (2.0, 40.0)]).unwrap(),
    )
    .unwrap();
    let a = linear
        .state([3.0, -1.0], [-0.7, 1.2], 0.2, [3.0, -1.0])
        .unwrap();
    let b = linear
        .state([-2.0, 4.0], [2.0, -1.5], 0.5, [2.0, -0.7])
        .unwrap();
    close(
        linear.system().hamiltonian(&a),
        saturated.system().hamiltonian(&a),
        2e-12,
    );
    let mid: [f64; 7] = std::array::from_fn(|k| f64::midpoint(a[k], b[k]));
    let (mut u, mut v) = ([0.0; 7], [0.0; 7]);
    linear
        .system()
        .storage()
        .discrete_gradient_rule()
        .unwrap()
        .evaluate(&a, &b, &mid, &mut u);
    saturated
        .system()
        .storage()
        .discrete_gradient_rule()
        .unwrap()
        .evaluate(&a, &b, &mid, &mut v);
    for k in 0..7 {
        close(u[k], v[k], 2e-12);
    }
    let mut linear = linear.prepare_controlled(budget(), 8, || false).unwrap();
    let mut saturated = saturated.prepare_controlled(budget(), 8, || false).unwrap();
    let (mut xa, mut xb, mut ya, mut yb) = ([0.0; 7], [0.0; 7], [0.0; 4], [0.0; 4]);
    let inputs = [7.0, -3.0, 0.4, -0.1];
    let first = linear
        .step_controlled(&a, &inputs, 0.002, &mut xa, &mut ya, || false)
        .unwrap();
    let second = saturated
        .step_controlled(&a, &inputs, 0.002, &mut xb, &mut yb, || false)
        .unwrap();
    for k in 0..7 {
        close(xa[k], xb[k], 2e-12);
    }
    for k in 0..4 {
        close(ya[k], yb[k], 2e-11);
    }
    close(first.native.delta_h, second.native.delta_h, 3e-12);
    close(first.rotor_copper_heat_j, second.rotor_copper_heat_j, 3e-12);
}

#[test]
fn unsupported_sources_corners_and_newton_trials_refuse_without_publishing() {
    let p = parameters();
    let asymmetric =
        InductorFluxCurve::try_new(vec![(-0.25, -20.0), (0.0, 0.0), (0.25, 25.0)]).unwrap();
    assert!(InductionMotor::with_magnetizing_curve(p, asymmetric).is_err());
    let model = model(p);
    assert!(model.state([26.0, 0.0], [0.0; 2], 0.0, [0.0; 2]).is_err());
    // Each actual endpoint is demagnetized; the crossed angle/flux corner is
    // outside the original radial support and cannot borrow endpoint admission.
    let balanced = self::model(InductionMotorParameters {
        rotor_leakage_inductance_h: 0.004,
        ..p
    });
    let a = [0.4, 0.0, -0.4, 0.0, 0.0, 0.0, 0.0];
    let b = [-0.4, 0.0, -0.4, 0.0, std::f64::consts::PI / 2.0, 0.0, 0.0];
    assert!(balanced.operating_point(&a).is_ok());
    assert!(balanced.operating_point(&b).is_ok());
    let mut effort = [0.0; 7];
    balanced
        .system()
        .storage()
        .discrete_gradient_rule()
        .unwrap()
        .evaluate(&a, &b, &a, &mut effort);
    assert!(effort.iter().all(|v| v.is_nan()));
    let initial = model.state([20.0, 0.0], [0.0; 2], 0.0, [0.0; 2]).unwrap();
    let mut prepared = model.prepare_controlled(budget(), 8, || false).unwrap();
    let (mut next, mut outputs) = ([123.0; 7], [456.0; 4]);
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
    assert_eq!(next, [123.0; 7]);
    assert_eq!(outputs, [456.0; 4]);
    assert!(matches!(
        prepared.step_controlled(&initial, &[0.0; 4], 1e-4, &mut next, &mut outputs, || true),
        Err(PreparedStepError::Cancelled)
    ));
    assert_eq!(next, [123.0; 7]);
    assert_eq!(outputs, [456.0; 4]);
    prepared
        .step_controlled(&initial, &[0.0; 4], 1e-4, &mut next, &mut outputs, || false)
        .unwrap();
    let mut fresh = self::model(p)
        .prepare_controlled(budget(), 8, || false)
        .unwrap();
    let (mut expected, mut expected_outputs) = ([0.0; 7], [0.0; 4]);
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

// Independent continuous original winding/shaft equations and physical accounts.
fn rhs(x: [f64; 12], inputs: [f64; 4]) -> [f64; 12] {
    let p = parameters();
    let state: [f64; 7] = x[..7].try_into().unwrap();
    let o = oracle(p, &state);
    let wr = x[5] / p.rotor_inertia_kg_m2;
    let ws = x[6] / p.stator_inertia_kg_m2;
    let bearing = p.bearing_drag_nm_s_per_rad * (wr - ws);
    [
        inputs[0] - p.stator_resistance_ohm * o.is[0],
        inputs[1] - p.stator_resistance_ohm * o.is[1],
        -p.rotor_resistance_ohm * o.ir[0],
        -p.rotor_resistance_ohm * o.ir[1],
        wr - ws,
        o.torque - bearing + inputs[2],
        -o.torque + bearing + inputs[3],
        inputs[0] * o.is[0] + inputs[1] * o.is[1] + inputs[2] * wr + inputs[3] * ws,
        p.stator_resistance_ohm * square(o.is),
        p.rotor_resistance_ohm * square(o.ir),
        bearing * (wr - ws),
        o.torque,
    ]
}

#[test]
fn saturated_transients_refine_to_independent_original_law_with_signed_work_and_heat() {
    let p = parameters();
    let model = model(p);
    let initial = model
        .state([6.0, 2.0], [-2.0, 1.0], 0.2, [4.0, -0.7])
        .unwrap();
    let inputs = [[8.0, -2.0, 0.2, -0.1], [-6.0, 3.0, -0.1, 0.05]];
    let duration = 0.003;
    let mut reference = [0.0; 12];
    reference[..7].copy_from_slice(&initial);
    for u in inputs {
        let h = duration / 2048.0;
        for _ in 0..2048 {
            let k1 = rhs(reference, u);
            let k2 = rhs(std::array::from_fn(|i| reference[i] + 0.5 * h * k1[i]), u);
            let k3 = rhs(std::array::from_fn(|i| reference[i] + 0.5 * h * k2[i]), u);
            let k4 = rhs(std::array::from_fn(|i| reference[i] + h * k3[i]), u);
            reference = std::array::from_fn(|i| {
                reference[i] + h * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0
            });
        }
    }
    let mut errors = [0.0; 3];
    for (refinement, count) in [4, 8, 16].into_iter().enumerate() {
        let mut prepared = self::model(p)
            .prepare_controlled(budget(), 8, || false)
            .unwrap();
        let mut actual = [0.0; 12];
        actual[..7].copy_from_slice(&initial);
        let mut regenerated = false;
        for u in inputs {
            let dt = duration / f64::from(count);
            for _ in 0..count {
                let state: [f64; 7] = actual[..7].try_into().unwrap();
                let (mut next, mut outputs) = ([0.0; 7], [0.0; 4]);
                let step = prepared
                    .step_controlled(&state, &u, dt, &mut next, &mut outputs, || false)
                    .unwrap();
                let electrical = step.electrical_work_j.iter().sum::<f64>();
                actual[7] += electrical + step.shaft_work_j.iter().sum::<f64>();
                actual[8] += step.stator_copper_heat_j;
                actual[9] += step.rotor_copper_heat_j;
                actual[10] += step.bearing_heat_j;
                actual[11] += step.electromagnetic_angular_impulse_nm_s;
                regenerated |= electrical < 0.0;
                close(
                    next[5] + next[6] - state[5] - state[6],
                    dt * (u[2] + u[3]),
                    2e-12,
                );
                close(step.native.balance_residual(), 0.0, 3e-11);
                actual[..7].copy_from_slice(&next);
            }
        }
        assert!(
            regenerated,
            "reverse stator voltage must recover original stored energy"
        );
        let final_state = actual[..7].try_into().unwrap();
        close(
            oracle(p, &final_state).energy - oracle(p, &initial).energy,
            actual[7] - actual[8] - actual[9] - actual[10],
            3e-11,
        );
        let scales = [
            0.1, 0.1, 0.1, 0.1, 1.0, 0.4, 0.14, 0.1, 0.1, 0.1, 0.01, 0.01,
        ];
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
        "second-order state/work/heat refinement: {errors:?}"
    );
    assert!(errors[2] < 2e-4, "refined physical trajectory: {errors:?}");
}

#[test]
fn resistance_feedback_preserves_the_original_saturation_source_and_energy() {
    let p = parameters();
    let original = model(p);
    let state = original
        .state([6.0, 2.0], [-2.0, 1.0], 0.2, [4.0, -0.7])
        .unwrap();
    let energy = original.system().hamiltonian(&state);
    let point = original.operating_point(&state).unwrap();
    let mut warm = original.prepare_controlled(budget(), 8, || false).unwrap();
    for resistance in [[0.0, 0.0], [1.6, 0.7], [0.2, 0.0]] {
        warm.set_winding_resistances(resistance).unwrap();
        assert_eq!(warm.magnetizing_curve(), Some(&curve()));
        assert_eq!(warm.system().hamiltonian(&state), energy);
        let rebound = warm.operating_point(&state).unwrap();
        assert_eq!(
            rebound.stator_current_alpha_beta_a,
            point.stator_current_alpha_beta_a
        );
        assert_eq!(rebound.rotor_current_dq_a, point.rotor_current_dq_a);
        assert_eq!(
            rebound.electromagnetic_torque_nm,
            point.electromagnetic_torque_nm
        );
        let mut fresh = model(InductionMotorParameters {
            stator_resistance_ohm: resistance[0],
            rotor_resistance_ohm: resistance[1],
            ..p
        })
        .prepare_controlled(budget(), 8, || false)
        .unwrap();
        let (mut actual, mut expected, mut y, mut z) = ([0.0; 7], [0.0; 7], [0.0; 4], [0.0; 4]);
        let inputs = [5.0, -2.0, 0.1, -0.1];
        let a = warm
            .step_controlled(&state, &inputs, 0.002, &mut actual, &mut y, || false)
            .unwrap();
        let b = fresh
            .step_controlled(&state, &inputs, 0.002, &mut expected, &mut z, || false)
            .unwrap();
        assert_eq!((actual, y, a), (expected, z, b));
    }
}
