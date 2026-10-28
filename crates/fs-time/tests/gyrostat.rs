//! G0/G1/G3 rotor momentum, analytic gyrostat and driven-reaction checks.
use fs_ga::{So3, So3Tangent, Vec3};
use fs_time::lie::{
    gyrostat::{HeldTorque, driven_step, momentum_step},
    rigid_body_momentum_step,
};

fn v(x: f64, y: f64, z: f64) -> Vec3 {
    Vec3::new(x, y, z)
}
fn mul(a: Vec3, b: Vec3) -> Vec3 {
    v(a.x * b.x, a.y * b.y, a.z * b.z)
}
fn omega(l: Vec3, stored: Vec3, i: Vec3) -> Vec3 {
    v(
        (l.x - stored.x) / i.x,
        (l.y - stored.y) / i.y,
        (l.z - stored.z) / i.z,
    )
}
fn energy(l: Vec3, stored: Vec3, i: Vec3) -> f64 {
    let w = omega(l, stored, i);
    0.5 * w.dot(mul(i, w))
}

#[test]
fn g0_total_momentum_and_bounded_energy_during_free_precession() {
    let inertia = v(1.0, 1.5, 2.0);
    let stored = v(0.7, -0.5, 0.9);
    let mut l = v(1.5, 0.2, 2.1);
    let mut r = So3::exp(So3Tangent::new(v(0.2, 0.3, -0.4))).unwrap();
    let initial = r.rotate(l).unwrap();
    let e0 = energy(l, stored, inertia);
    let mut max_energy_error = 0.0_f64;
    let mut max_momentum_error = 0.0_f64;
    for _ in 0..10_000 {
        (r, l) = momentum_step(r, l, inertia, stored, 0.001).unwrap();
        max_energy_error = max_energy_error.max((energy(l, stored, inertia) - e0).abs());
        max_momentum_error = max_momentum_error.max((r.rotate(l).unwrap() - initial).norm());
    }
    println!("gyrostat: energy {max_energy_error:e}; spatial momentum {max_momentum_error:e}");
    assert!(max_energy_error < 3e-6);
    assert!(max_momentum_error < 3e-10);
}

#[test]
fn g1_spherical_gyrostat_converges_to_exact_noncommuting_attitude() {
    let inertia = v(2.0, 2.0, 2.0);
    let l0 = v(1.2, -0.4, 0.8);
    let stored = v(-0.3, 0.7, 0.6);
    let r0 = So3::exp(So3Tangent::new(v(-0.1, 0.5, 0.2))).unwrap();
    let duration = 2.0;
    let spatial = So3::exp(So3Tangent::new(l0.scale(duration / 2.0))).unwrap();
    let rotor = So3::exp(So3Tangent::new(stored.scale(-duration / 2.0))).unwrap();
    let exact_r = r0.compose(spatial).unwrap().compose(rotor).unwrap();
    let exact_l = rotor.inverse().rotate(l0).unwrap();
    let run = |n: usize| {
        let (mut r, mut l) = (r0, l0);
        for _ in 0..n {
            (r, l) = momentum_step(r, l, inertia, stored, duration / n as f64).unwrap();
        }
        r.body_minus(exact_r).unwrap().angular.norm() + (l - exact_l).norm()
    };
    let errors = [run(40), run(80), run(160)];
    println!("gyrostat analytic errors {errors:?}");
    assert!(errors[2] < 5e-5);
    assert!((3.8..4.2).contains(&(errors[0] / errors[1])));
    assert!((3.8..4.2).contains(&(errors[1] / errors[2])));
}

#[test]
fn g1_motor_spinup_has_exact_opposite_carrier_reaction() {
    let inertia = v(1.0, 2.0, 3.0);
    let r0 = So3::exp(So3Tangent::new(v(0.2, -0.3, 0.4))).unwrap();
    let zero = v(0.0, 0.0, 0.0);
    for axis in 0..3 {
        let mut h = [0.0; 3];
        h[axis] = 0.7;
        let stored = v(h[0], h[1], h[2]);
        let dt = 0.4;
        let (r, w) =
            driven_step(r0, zero, inertia, zero, stored, HeldTorque::Space(zero), dt).unwrap();
        let expected_w = omega(zero, stored, inertia);
        let expected_r = r0
            .body_plus(So3Tangent::new(expected_w.scale(dt / 2.0)))
            .unwrap();
        assert!((w - expected_w).norm() < 2e-14);
        assert!(r.body_minus(expected_r).unwrap().angular.norm() < 2e-14);
        assert!(r.rotate(mul(inertia, w) + stored).unwrap().norm() < 2e-14);
    }
}

#[test]
fn g0_driven_rotors_obey_external_impulse_in_each_named_frame() {
    let inertia = v(1.0, 1.5, 2.0);
    let torque = v(0.03, -0.02, 0.01);
    let rate = v(-0.04, 0.05, 0.02);
    let dt = 0.002;
    for space in [false, true] {
        let mut r = So3::exp(So3Tangent::new(v(0.2, -0.3, 0.4))).unwrap();
        let mut w = v(0.4, 0.7, 1.1);
        let mut stored = v(0.2, -0.1, 0.3);
        let initial = r.rotate(mul(inertia, w) + stored).unwrap();
        let mut impulse = v(0.0, 0.0, 0.0);
        for _ in 0..1_000 {
            let before_r = r;
            let next_stored = stored + rate.scale(dt);
            let load = if space {
                HeldTorque::Space(torque)
            } else {
                HeldTorque::Body(torque)
            };
            (r, w) = driven_step(r, w, inertia, stored, next_stored, load, dt).unwrap();
            stored = next_stored;
            impulse = impulse
                + if space {
                    torque.scale(dt)
                } else {
                    (before_r.rotate(torque).unwrap() + r.rotate(torque).unwrap()).scale(dt / 2.0)
                };
        }
        let error = (r.rotate(mul(inertia, w) + stored).unwrap() - initial - impulse).norm();
        println!("gyrostat driven space={space}: impulse defect {error:e}");
        assert!(error < 2e-10);
    }
}

#[test]
fn g3_reversible_covariant_and_bit_replayable() {
    let inertia = v(1.0, 1.5, 2.0);
    let stored = v(0.7, -0.5, 0.9);
    let l0 = v(1.5, 0.2, 2.1);
    let r0 = So3::exp(So3Tangent::new(v(0.2, 0.3, -0.4))).unwrap();
    let frame = So3::exp(So3Tangent::new(v(-0.4, 0.1, 0.8))).unwrap();
    let run = |start| {
        let (mut r, mut l) = (start, l0);
        for _ in 0..256 {
            (r, l) = momentum_step(r, l, inertia, stored, 0.01).unwrap();
        }
        (r, l)
    };
    let (mut r, mut l) = run(r0);
    assert_eq!((r, l), run(r0));
    let transformed = run(frame.compose(r0).unwrap());
    assert_eq!(l, transformed.1);
    assert!(
        frame
            .compose(r)
            .unwrap()
            .body_minus(transformed.0)
            .unwrap()
            .angular
            .norm()
            < 2e-12
    );
    for _ in 0..256 {
        (r, l) = momentum_step(r, l, inertia, stored, -0.01).unwrap();
    }
    assert!((l - l0).norm() < 3e-12);
    assert!(r.body_minus(r0).unwrap().angular.norm() < 3e-12);
}

#[test]
fn g0_equilibrium_rigid_reduction_and_invalid_states() {
    let zero = v(0.0, 0.0, 0.0);
    let inertia = v(1.0, 1.5, 2.0);
    let stored = v(20.0, -50.0, 100.0);
    let rotation = So3::exp(So3Tangent::new(v(0.2, 0.3, -0.4))).unwrap();
    assert_eq!(
        momentum_step(rotation, stored, inertia, stored, 0.1).unwrap(),
        (rotation, stored)
    );
    assert_eq!(
        momentum_step(rotation, stored, inertia, zero, 0.01).unwrap(),
        rigid_body_momentum_step(rotation, stored, inertia, 0.01).unwrap()
    );
    for bad in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(momentum_step(rotation, stored, v(bad, 1.0, 2.0), stored, 0.01).is_err());
    }
    for bad in [0.0, f64::INFINITY, f64::NAN] {
        assert!(momentum_step(rotation, stored, inertia, stored, bad).is_err());
    }
    for bad in [f64::INFINITY, f64::NEG_INFINITY, f64::NAN] {
        let invalid = v(bad, 0.1, 0.2);
        assert!(momentum_step(rotation, invalid, inertia, stored, 0.01).is_err());
        assert!(momentum_step(rotation, stored, inertia, invalid, 0.01).is_err());
        assert!(
            driven_step(
                rotation,
                zero,
                inertia,
                stored,
                invalid,
                HeldTorque::Body(zero),
                0.01
            )
            .is_err()
        );
        assert!(
            driven_step(
                rotation,
                zero,
                inertia,
                stored,
                stored,
                HeldTorque::Body(invalid),
                0.01
            )
            .is_err()
        );
    }
}

#[test]
fn g3_driven_map_is_second_order_and_reverses_internal_transfer() {
    let inertia = v(1.0, 1.5, 2.0);
    let initial_r = So3::exp(So3Tangent::new(v(0.2, -0.3, 0.4))).unwrap();
    let initial_w = v(0.4, 0.7, 1.1);
    let initial_h = v(0.7, -0.5, 0.9);
    let rate = v(-0.1, 0.3, -0.2);
    let torque = HeldTorque::Space(v(0.2, 0.1, -0.1));
    let run = |n: usize| {
        let (mut r, mut w) = (initial_r, initial_w);
        let dt = 1.0 / n as f64;
        for k in 0..n {
            let h0 = initial_h + rate.scale(k as f64 * dt);
            let h1 = initial_h + rate.scale((k + 1) as f64 * dt);
            (r, w) = driven_step(r, w, inertia, h0, h1, torque, dt).unwrap();
        }
        (r, w)
    };
    let distance = |a: (So3, Vec3), b: (So3, Vec3)| {
        a.0.body_minus(b.0).unwrap().angular.norm() + (a.1 - b.1).norm()
    };
    let fine = run(160);
    let ratio = distance(run(40), run(80)) / distance(run(80), fine);
    assert!((3.8..4.2).contains(&ratio), "driven order ratio {ratio}");
    let (mut r, mut w) = fine;
    for k in (0..160).rev() {
        let h0 = initial_h + rate.scale((k + 1) as f64 / 160.0);
        let h1 = initial_h + rate.scale(k as f64 / 160.0);
        (r, w) = driven_step(r, w, inertia, h0, h1, torque, -1.0 / 160.0).unwrap();
    }
    assert!((w - initial_w).norm() < 2e-12);
    assert!(r.body_minus(initial_r).unwrap().angular.norm() < 2e-12);
}
