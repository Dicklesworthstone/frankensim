//! Separate-carrier gyrostat regressions: observable small motion beside stored spin.
use fs_ga::{Quat, So3, So3Tangent, Vec3};
use fs_time::lie::{
    gyrostat::{HeldTorque, InertiaTensor3, carrier_momentum_step, driven_step, momentum_step},
    rigid_body_momentum_step,
};

fn v(x: f64, y: f64, z: f64) -> Vec3 {
    Vec3::new(x, y, z)
}
fn axis(i: usize, value: f64) -> Vec3 {
    let mut x = [0.0; 3];
    x[i] = value;
    v(x[0], x[1], x[2])
}
fn exp(value: Vec3) -> So3 {
    So3::exp(So3Tangent::new(value)).unwrap()
}
fn vector_close(actual: Vec3, expected: Vec3, tolerance: f64) {
    assert!(
        (actual - expected).norm() < tolerance,
        "{actual:?} != {expected:?}"
    );
}
fn rotation_close(actual: So3, expected: So3, tolerance: f64) {
    let error = actual.body_minus(expected).unwrap().angular.norm();
    assert!(error < tolerance, "rotation error {error:e}");
}

#[test]
fn carrier_axial_motion_survives_an_unresolvable_total_sum() {
    let r0 = exp(v(0.2, -0.3, 0.4));
    for inertia in [v(1.0, 1.5, 2.0), v(2.0, 2.0, 2.0), v(2.0, 2.0, 5.0)] {
        let moments = [inertia.x, inertia.y, inertia.z];
        for i in 0..3 {
            for sign in [-1.0, 1.0] {
                let stored = axis(i, sign * 1e30);
                let k = axis(i, 0.75 * moments[i]);
                assert_eq!(stored + k, stored, "the old representation loses k");
                for dt in [-0.7, 0.4] {
                    let expected_r = r0.compose(exp(axis(i, 0.75 * dt))).unwrap();
                    let drift = carrier_momentum_step(r0, k, inertia, stored, dt).unwrap();
                    assert_eq!(drift.1, k);
                    rotation_close(drift.0, expected_r, 2e-15);
                    let driven = driven_step(
                        r0,
                        axis(i, 0.75),
                        inertia,
                        stored,
                        stored,
                        HeldTorque::Body(Vec3::default()),
                        dt,
                    )
                    .unwrap();
                    assert_eq!(driven.1, axis(i, 0.75));
                    rotation_close(driven.0, expected_r, 2e-15);
                }
            }
        }
    }
}

#[test]
fn carrier_large_axial_spin_retains_external_acceleration() {
    let inertia = v(1.0, 1.5, 2.0);
    let r0 = So3::identity();
    let h = v(1e30, 0.0, 0.0);
    let omega0 = v(0.75, 0.0, 0.0);
    let dt = 0.4;
    let torque = v(0.5, 0.0, 0.0);
    for load in [
        HeldTorque::Body(torque),
        HeldTorque::Space(r0.rotate(torque).unwrap()),
    ] {
        let (r1, omega1) = driven_step(r0, omega0, inertia, h, h, load, dt).unwrap();
        vector_close(omega1, v(0.95, 0.0, 0.0), 3e-15);
        rotation_close(r1, r0.compose(exp(v(0.85 * dt, 0.0, 0.0))).unwrap(), 3e-15);
    }
}

#[test]
fn carrier_opposing_external_and_internal_impulses_do_not_erase_motion() {
    // The two large half-impulses cancel before reaching k. Sequentially
    // adding either to k would round away the carrier's initial unit momentum.
    let r0 = exp(v(0.2, -0.3, 0.4));
    let inertia = v(1.0, 1.5, 2.0);
    let h0 = v(1e18, 0.0, 0.0);
    let h1 = v(3e18, 0.0, 0.0);
    let (r1, w1) = driven_step(
        r0,
        v(1.0, 0.0, 0.0),
        inertia,
        h0,
        h1,
        HeldTorque::Body(v(2e18, 0.0, 0.0)),
        1.0,
    )
    .unwrap();
    assert_eq!(w1, v(1.0, 0.0, 0.0));
    rotation_close(r1, r0.compose(exp(v(1.0, 0.0, 0.0))).unwrap(), 2e-15);

    // A finite, representable axial spin-up still produces exactly the
    // opposite carrier impulse when there is no external torque.
    let before = v(1e12, 0.0, 0.0);
    let after = before + v(4.0, 0.0, 0.0);
    let (r2, w2) = driven_step(
        r0,
        v(1.0, 0.0, 0.0),
        inertia,
        before,
        after,
        HeldTorque::Space(Vec3::default()),
        0.2,
    )
    .unwrap();
    assert_eq!(w2, v(-3.0, 0.0, 0.0));
    rotation_close(r2, r0.compose(exp(v(-0.2, 0.0, 0.0))).unwrap(), 2e-15);
}

#[test]
fn carrier_balanced_dyadic_transfer_retains_the_small_endpoint_residual() {
    // All inputs are exact dyadics and tau*dt == h1-h0 mathematically.
    // Rounding the rotor half-transfer before the compensated kick loses
    // 2^-55 per endpoint and incorrectly reduces this carrier to rest.
    let tiny = 1.0 / (1_u64 << 54) as f64;
    let delta = 1.0 / (1_u64 << 27) as f64;
    let dt = 1.0 + delta;
    let torque = v(1.0 - delta, 0.0, 0.0);
    let h0 = v(tiny, 0.0, 0.0);
    let h1 = v(1.0, 0.0, 0.0);
    let inertia = v(1.0, 1.5, 2.0);
    let initial_omega = v(tiny, 0.0, 0.0);
    let native = driven_step(
        So3::identity(),
        initial_omega,
        inertia,
        h0,
        h1,
        HeldTorque::Body(torque),
        dt,
    )
    .unwrap();
    assert_eq!(native.1, initial_omega);
    rotation_close(native.0, exp(initial_omega.scale(dt)), 1e-30);
    let permutation = So3::try_from_quat(Quat {
        w: 0.5,
        x: 0.5,
        y: 0.5,
        z: 0.5,
    })
    .unwrap();
    let tensor = InertiaTensor3::from_principal(inertia, permutation).unwrap();
    let omega_body = permutation.rotate(initial_omega).unwrap();
    let body = tensor
        .driven_step(
            So3::identity(),
            omega_body,
            permutation.rotate(h0).unwrap(),
            permutation.rotate(h1).unwrap(),
            HeldTorque::Body(permutation.rotate(torque).unwrap()),
            dt,
        )
        .unwrap();
    assert_eq!(body.1, omega_body);
    rotation_close(body.0, exp(omega_body.scale(dt)), 1e-30);
}

#[test]
fn carrier_symmetric_flow_matches_independent_euler_solution() {
    let r0 = exp(v(-0.3, 0.2, 0.4));
    let k0 = v(1.2, -0.4, 0.8);
    for i in 0..3 {
        for axial in [0.7, 5.0] {
            for spin in [-0.6, 0.6] {
                let mut values = [2.0; 3];
                values[i] = axial;
                let inertia = v(values[0], values[1], values[2]);
                let h = axis(i, spin);
                let kp = [k0.x, k0.y, k0.z][i];
                let residual_rate = spin / 2.0 + kp * (1.0 / 2.0 - 1.0 / axial);
                for dt in [-1.3, 0.7, 2.0] {
                    let residual = exp(axis(i, -dt * residual_rate));
                    let exact_r = r0
                        .compose(exp((k0 + h).scale(dt / 2.0)))
                        .unwrap()
                        .compose(residual)
                        .unwrap();
                    let exact_k = residual.inverse().rotate(k0).unwrap();
                    let actual = carrier_momentum_step(r0, k0, inertia, h, dt).unwrap();
                    rotation_close(actual.0, exact_r, 4e-14);
                    vector_close(actual.1, exact_k, 4e-14);
                    vector_close(
                        actual.0.rotate(actual.1 + h).unwrap(),
                        r0.rotate(k0 + h).unwrap(),
                        4e-14,
                    );
                    let reversed =
                        carrier_momentum_step(actual.0, actual.1, inertia, h, -dt).unwrap();
                    rotation_close(reversed.0, r0, 4e-14);
                    vector_close(reversed.1, k0, 4e-14);
                }
            }
        }
    }
    let inertia = v(2.0, 2.0, 2.0);
    let h = v(-0.3, 0.7, 0.6);
    for k in [
        k0,
        h.scale(-1.0),
        h.scale(-2.0),
        h.scale(-2.0) + v(0.2, -0.1, 0.3),
    ] {
        for dt in [-1.3, 0.7] {
            let body = exp(h.scale(-dt / 2.0));
            let expected = r0
                .compose(exp((k + h).scale(dt / 2.0)))
                .unwrap()
                .compose(body)
                .unwrap();
            let actual = carrier_momentum_step(r0, k, inertia, h, dt).unwrap();
            rotation_close(actual.0, expected, 4e-14);
            vector_close(actual.1, body.inverse().rotate(k).unwrap(), 4e-14);
        }
    }
}

#[test]
fn carrier_nearly_axial_high_spin_keeps_its_small_attitude_phase() {
    let r0 = exp(v(-0.3, 0.2, 0.4));
    let k = v(1.25, -0.75, 2.0);
    let inertia = v(2.0, 2.0, 5.0);
    for spin in [-1e30, 1e30] {
        let h = v(0.0, 0.0, spin);
        let dt = 0.7;
        let actual = carrier_momentum_step(r0, k, inertia, h, dt).unwrap();
        // Transverse attitude correction is bounded by O(|k_perp|/|h|).
        // The order-one axial phase must remain even though |h+k| rounds to |h|.
        let axial = r0.compose(exp(v(0.0, 0.0, dt * k.z / inertia.z))).unwrap();
        rotation_close(actual.0, axial, 4e-15);
        // Euler's transverse carrier precession has a large rotor phase and a
        // separately observable small asymmetric-inertia phase.
        let precession = exp(v(0.0, 0.0, dt * (spin / inertia.x)))
            .compose(exp(v(
                0.0,
                0.0,
                dt * k.z * (1.0 / inertia.x - 1.0 / inertia.z),
            )))
            .unwrap();
        vector_close(actual.1, precession.rotate(k).unwrap(), 4e-15);
        assert!((actual.1.dot(actual.1) - k.dot(k)).abs() < 5e-15);
    }
}

#[test]
fn carrier_general_drift_matches_total_map_and_is_reversible() {
    let r0 = exp(v(0.2, -0.3, 0.4));
    let inertia = v(1.0, 1.5, 2.0);
    let h = v(0.7, -0.5, 0.9);
    let k = v(0.8, 0.7, 1.2);
    for dt in [-0.1, 0.03, 0.5] {
        let separate = carrier_momentum_step(r0, k, inertia, h, dt).unwrap();
        let total = momentum_step(r0, k + h, inertia, h, dt).unwrap();
        rotation_close(separate.0, total.0, 3e-15);
        vector_close(separate.1 + h, total.1, 3e-15);
        let backward = carrier_momentum_step(separate.0, separate.1, inertia, h, -dt).unwrap();
        rotation_close(backward.0, r0, 3e-15);
        vector_close(backward.1, k, 3e-15);
    }
    let initial_world = r0.rotate(k + h).unwrap();
    let mut state = (r0, k);
    for _ in 0..1000 {
        state = carrier_momentum_step(state.0, state.1, inertia, h, 0.001).unwrap();
    }
    vector_close(state.0.rotate(state.1 + h).unwrap(), initial_world, 2e-12);
}

#[test]
fn carrier_tensor_coordinates_preserve_motion_and_covariance() {
    let r0 = exp(v(0.2, -0.3, 0.4));
    let principal = v(1.0, 1.5, 2.0);
    // Quaternion (1,1,1,1)/2 permutes axes exactly, including at extreme h/k.
    let permutation = So3::try_from_quat(Quat {
        w: 0.5,
        x: 0.5,
        y: 0.5,
        z: 0.5,
    })
    .unwrap();
    let tensor = InertiaTensor3::from_principal(principal, permutation).unwrap();
    let kp = v(0.75, 0.0, 0.0);
    let hp = v(1e30, 0.0, 0.0);
    let kb = permutation.rotate(kp).unwrap();
    let hb = permutation.rotate(hp).unwrap();
    let result = tensor.carrier_momentum_step(r0, kb, hb, 0.4).unwrap();
    assert_eq!(result.1, kb);
    rotation_close(
        result.0,
        r0.compose(exp(permutation.rotate(v(0.3, 0.0, 0.0)).unwrap()))
            .unwrap(),
        3e-15,
    );
    let driven = tensor
        .driven_step(r0, kb, hb, hb, HeldTorque::Body(Vec3::default()), 0.4)
        .unwrap();
    assert_eq!(driven.1, kb);
    rotation_close(driven.0, result.0, 3e-15);

    let basis = exp(v(-0.3, 0.7, 0.2));
    let tensor = InertiaTensor3::from_principal(principal, basis).unwrap();
    for spin_scale in [1.0, 1e4] {
        let kp = v(0.6, -0.4, 0.8);
        let hp = v(0.3, -0.2, 0.7).scale(spin_scale);
        let kb = basis.rotate(kp).unwrap();
        let hb = basis.rotate(hp).unwrap();
        let dt = 0.02 / spin_scale;
        let body = tensor.carrier_momentum_step(r0, kb, hb, dt).unwrap();
        let native =
            carrier_momentum_step(r0.compose(basis).unwrap(), kp, principal, hp, dt).unwrap();
        rotation_close(body.0, native.0.compose(basis.inverse()).unwrap(), 3e-14);
        vector_close(body.1, basis.rotate(native.1).unwrap(), 3e-14);
    }
    let diagonal = InertiaTensor3::from_principal(principal, So3::identity()).unwrap();
    assert_eq!(
        diagonal.carrier_momentum_step(r0, kp, hp, 0.4).unwrap(),
        carrier_momentum_step(r0, kp, principal, hp, 0.4).unwrap()
    );
}

#[test]
fn carrier_zero_spin_keeps_legacy_rigid_arithmetic_and_invalid_states_refuse() {
    let r = exp(v(0.2, -0.3, 0.4));
    let k = v(0.8, -0.4, 0.7);
    let inertia = v(1.0, 1.5, 2.0);
    let zero = Vec3::default();
    assert_eq!(
        carrier_momentum_step(r, k, inertia, zero, 0.01).unwrap(),
        rigid_body_momentum_step(r, k, inertia, 0.01).unwrap()
    );
    let tensor = InertiaTensor3::from_principal(inertia, exp(v(0.2, 0.1, -0.4))).unwrap();
    assert_eq!(
        tensor.carrier_momentum_step(r, k, zero, 0.01).unwrap(),
        tensor.momentum_step(r, k, zero, 0.01).unwrap()
    );
    for dt in [0.0, f64::INFINITY, f64::NAN] {
        assert!(carrier_momentum_step(r, k, inertia, v(1.0, 0.0, 0.0), dt).is_err());
    }
    for bad in [f64::INFINITY, f64::NEG_INFINITY, f64::NAN] {
        assert!(
            carrier_momentum_step(r, v(bad, 0.0, 0.0), inertia, v(1.0, 0.0, 0.0), 0.01).is_err()
        );
        assert!(carrier_momentum_step(r, k, inertia, v(bad, 0.0, 0.0), 0.01).is_err());
    }
}
