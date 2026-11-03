use fs_mbd::rigid_rattle::distance::{DistanceStepControls, DistanceStepError};
use fs_mbd::rigid_rattle::finite_plane::FinitePointPlaneConstraint;
use fs_mbd::{
    Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState, UnitQuaternion, Vec3,
    Wrench,
};

fn near(a: f64, b: f64, t: f64) {
    assert!((a - b).abs() <= t, "{a:e} != {b:e}, tolerance {t:e}");
}
fn vector(a: Vec3, b: Vec3, t: f64) {
    near(a.x, b.x, t);
    near(a.y, b.y, t);
    near(a.z, b.z, t);
}
fn controls() -> DistanceStepControls {
    DistanceStepControls {
        maximum_impulse_ns: 100.0,
        maximum_rotation_increment_rad: 0.25,
        maximum_position_iterations: 50,
        maximum_backtracks: 20,
    }
}
fn free() -> RigidBodyIntegrator {
    RigidBodyIntegrator::new(Gravity::new(Vec3::ZERO).unwrap())
}
fn inertia(m: f64, i: f64) -> MassProperties {
    MassProperties::new(m, Vec3::ZERO, Vec3::new(i, i, i)).unwrap()
}
fn state(p: MassProperties, x: Vec3, q: UnitQuaternion, v: Vec3, w: Vec3) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(x, q).unwrap(),
        v.scale(p.mass()),
        p.angular_momentum_body_checked(w).unwrap(),
    )
    .unwrap()
}
fn joint(p: [MassProperties; 2], a: Vec3, b: Vec3, n: Vec3) -> FinitePointPlaneConstraint {
    FinitePointPlaneConstraint::new(p, a, b, n, 1e-12, 1e-12).unwrap()
}
fn kinetic(s: RigidBodyState, p: MassProperties) -> f64 {
    0.5 * s.linear_momentum_world().norm_squared() / p.mass()
        + 0.5
            * s.angular_momentum_body()
                .dot(p.angular_velocity_body(s.angular_momentum_body()))
}
fn inventories(s: [RigidBodyState; 2], p: [MassProperties; 2]) -> (Vec3, Vec3, f64) {
    let mut linear = Vec3::ZERO;
    let mut angular = Vec3::ZERO;
    let mut energy = 0.0;
    for i in 0..2 {
        linear = linear.add(s[i].linear_momentum_world());
        angular = angular
            .add(
                s[i].pose()
                    .position_world()
                    .cross(s[i].linear_momentum_world()),
            )
            .add(
                s[i].pose()
                    .orientation()
                    .rotate_body_to_world(s[i].angular_momentum_body()),
            );
        energy += kinetic(s[i], p[i]);
    }
    (linear, angular, energy)
}

#[test]
fn g1_finite_carrier_shares_normal_load_by_actual_unequal_mass() {
    for masses in [[1.0, 1.0], [2.0, 5.0], [0.1, 100.0]] {
        let p = [inertia(masses[0], 0.7), inertia(masses[1], 1.3)];
        for n in [Vec3::new(0.0, 1.0, 0.0), Vec3::new(0.0, 0.6, 0.8)] {
            let c = joint(p, Vec3::ZERO, Vec3::ZERO, n);
            let x = Vec3::new(0.2, -0.3, 0.4);
            let v = Vec3::new(0.4, 0.0, 0.0).add(n.scale(-0.2));
            let initial = [
                state(p[0], x, UnitQuaternion::IDENTITY, v, Vec3::ZERO),
                state(p[1], x, UnitQuaternion::IDENTITY, v, Vec3::ZERO),
            ];
            let forces = [2.0, -0.7];
            let w = forces.map(|f| Wrench {
                force_world: n.scale(f),
                torque_body: Vec3::ZERO,
            });
            let dt = 0.03;
            let out = c
                .advance(initial, free(), w, dt, controls(), || false)
                .unwrap();
            let acceleration = (forces[0] + forces[1]) / (masses[0] + masses[1]);
            let impulse = -0.5 * dt * (forces[0] / masses[0] - forces[1] / masses[1])
                / (1.0 / masses[0] + 1.0 / masses[1]);
            for i in 0..2 {
                vector(
                    out.states_after[i].pose().position_world(),
                    x.add(v.scale(dt))
                        .add(n.scale(0.5 * acceleration * dt * dt)),
                    2e-12,
                );
                vector(
                    out.states_after[i].linear_momentum_world(),
                    v.add(n.scale(acceleration * dt)).scale(masses[i]),
                    2e-11,
                );
                vector(
                    out.states_after[i].angular_momentum_body(),
                    Vec3::ZERO,
                    3e-13,
                );
            }
            vector(
                out.position_impulse.first.impulse_world,
                n.scale(impulse),
                2e-11,
            );
            vector(
                out.velocity_impulse.first.impulse_world,
                n.scale(impulse),
                2e-11,
            );
            near(out.constraint_work_j, 0.0, 2e-12);
            near(out.final_gap_m, 0.0, 1e-12);
            near(out.final_normal_speed_mps, 0.0, 1e-12);
        }
    }
}

// Independent reduced equations for a point mass sliding on a freely rotating
// finite carrier. Relative reduced mass is mu; conserved L=(I+mu*r²)*theta_dot.
fn sliding_reference(mut y: [f64; 3], l: f64, mu: f64, i: f64, time: f64) -> [f64; 3] {
    let f = |s: [f64; 3]| {
        let w = l / (i + mu * s[0] * s[0]);
        [s[1], s[0] * w * w, w]
    };
    let h = time / 4096.0;
    for _ in 0..4096 {
        let a = f(y);
        let b = f(core::array::from_fn(|k| y[k] + 0.5 * h * a[k]));
        let c = f(core::array::from_fn(|k| y[k] + 0.5 * h * b[k]));
        let d = f(core::array::from_fn(|k| y[k] + h * c[k]));
        y = core::array::from_fn(|k| y[k] + h * (a[k] + 2.0 * b[k] + 2.0 * c[k] + d[k]) / 6.0);
    }
    y
}

#[test]
fn g1_sliding_rotating_finite_carrier_matches_independent_continuous_dynamics() {
    let p = [inertia(1.2, 0.5), inertia(2.3, 0.8)];
    let total = p[0].mass() + p[1].mass();
    let mu = p[0].mass() * p[1].mass() / total;
    let r = 0.7;
    let radial = 0.2;
    let omega = 1.1;
    let horizon = 0.6;
    let relative_x = Vec3::new(r, 0.0, 0.0);
    let relative_v = Vec3::new(radial, r * omega, 0.0);
    let initial = [
        state(
            p[0],
            relative_x.scale(p[1].mass() / total),
            UnitQuaternion::IDENTITY,
            relative_v.scale(p[1].mass() / total),
            Vec3::ZERO,
        ),
        state(
            p[1],
            relative_x.scale(-p[0].mass() / total),
            UnitQuaternion::IDENTITY,
            relative_v.scale(-p[0].mass() / total),
            Vec3::new(0.0, 0.0, omega),
        ),
    ];
    let c = joint(p, Vec3::ZERO, Vec3::ZERO, Vec3::new(0.0, 1.0, 0.0));
    let l = (0.8 + mu * r * r) * omega;
    let expected = sliding_reference([r, radial, 0.0], l, mu, 0.8, horizon);
    let expected_x = Vec3::new(expected[2].cos(), expected[2].sin(), 0.0).scale(expected[0]);
    let mut errors = Vec::new();
    for steps in [24, 48, 96] {
        let h = horizon / f64::from(steps);
        let mut s = initial;
        let mut work = 0.0;
        let before = inventories(s, p);
        for _ in 0..steps {
            let out = c
                .advance(s, free(), [Wrench::ZERO; 2], h, controls(), || false)
                .unwrap();
            work += out.constraint_work_j;
            s = out.states_after;
        }
        let observed = s[0]
            .pose()
            .position_world()
            .sub(s[1].pose().position_world());
        errors.push(
            observed.sub(expected_x).norm_squared().sqrt()
                + (s[1].angular_momentum_body().z / 0.8
                    - l / (0.8 + mu * expected[0] * expected[0]))
                    .abs(),
        );
        let after = inventories(s, p);
        vector(after.0, before.0, 2e-12);
        vector(after.1, before.1, 2e-11);
        near(after.2 - before.2, work, 2e-12);
        assert!(s[1].angular_momentum_body().z < initial[1].angular_momentum_body().z);
        s = s.map(|s| {
            RigidBodyState::new(
                s.pose(),
                s.linear_momentum_world().scale(-1.0),
                s.angular_momentum_body().scale(-1.0),
            )
            .unwrap()
        });
        for _ in 0..steps {
            s = c
                .advance(s, free(), [Wrench::ZERO; 2], h, controls(), || false)
                .unwrap()
                .states_after;
        }
        for i in 0..2 {
            vector(
                s[i].pose().position_world(),
                initial[i].pose().position_world(),
                3e-9,
            );
            vector(
                s[i].linear_momentum_world(),
                initial[i].linear_momentum_world().scale(-1.0),
                4e-9,
            );
        }
    }
    for k in 0..2 {
        assert!(
            errors[k] / errors[k + 1] > 3.8 && errors[k] / errors[k + 1] < 4.2,
            "{errors:?}"
        );
    }
}

fn full_fixture() -> ([MassProperties; 2], [RigidBodyState; 2], Vec3, Vec3, Vec3) {
    let p = [
        MassProperties::from_principal_frame(
            1.7,
            Vec3::ZERO,
            Vec3::new(0.7, 1.0, 1.2),
            UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, -1.0), 0.6).unwrap(),
        )
        .unwrap(),
        MassProperties::from_principal_frame(
            2.2,
            Vec3::ZERO,
            Vec3::new(0.9, 1.1, 1.5),
            UnitQuaternion::from_axis_angle(Vec3::new(-2.0, 1.0, 3.0), 0.4).unwrap(),
        )
        .unwrap(),
    ];
    let qa = UnitQuaternion::from_axis_angle(Vec3::new(0.4, -0.7, 1.0), 0.5).unwrap();
    let qb = UnitQuaternion::from_axis_angle(Vec3::new(-0.3, 0.2, 0.8), 0.4).unwrap();
    let arm = Vec3::new(0.5, -0.3, 0.2);
    let anchor = Vec3::new(0.2, 0.3, -0.1);
    let normal = Vec3::new(0.0, 1.0, 0.0);
    let wb = Vec3::new(0.2, -0.1, 0.4);
    let wa = Vec3::new(-0.3, 0.2, 0.1);
    let xb = Vec3::new(0.1, -0.4, 0.2);
    let r = qb.rotate_body_to_world(anchor.add(Vec3::new(0.6, 0.0, -0.2)));
    let xa = xb.add(r).sub(qa.rotate_body_to_world(arm));
    let vb = Vec3::new(0.1, 0.3, -0.2);
    let va = vb
        .add(qb.rotate_body_to_world(wb).cross(r))
        .sub(
            qa.rotate_body_to_world(wa)
                .cross(qa.rotate_body_to_world(arm)),
        )
        .add(qb.rotate_body_to_world(Vec3::new(0.3, 0.0, 0.1)));
    (
        p,
        [state(p[0], xa, qa, va, wa), state(p[1], xb, qb, vb, wb)],
        arm,
        anchor,
        normal,
    )
}

#[test]
fn g0_full_tensor_sliding_transfers_internal_momentum_and_measured_work() {
    let (p, initial, arm, anchor, n) = full_fixture();
    let c = joint(p, arm, anchor, n);
    let mut s = initial;
    let first = inventories(s, p);
    let mut work = 0.0;
    let mut free_defect = 0.0;
    for _ in 0..200 {
        let out = c
            .advance(s, free(), [Wrench::ZERO; 2], 0.002, controls(), || false)
            .unwrap();
        work += out.constraint_work_j;
        for i in 0..2 {
            free_defect += kinetic(out.free_steps[i].state_after, p[i])
                - kinetic(out.free_steps[i].state_before, p[i]);
        }
        s = out.states_after;
        let now = inventories(s, p);
        vector(now.0, first.0, 3e-12);
        vector(now.1, first.1, 1e-10);
        near(out.final_gap_m, 0.0, 1e-12);
        near(out.final_normal_speed_mps, 0.0, 1e-12);
    }
    near(inventories(s, p).2 - first.2, work + free_defect, 3e-12);
    assert!(
        s[1].angular_momentum_body()
            .sub(initial[1].angular_momentum_body())
            .norm_squared()
            > 1e-7
    );
}

fn compose(a: UnitQuaternion, b: UnitQuaternion) -> UnitQuaternion {
    let [w, x, y, z] = a.components();
    let [s, u, v, t] = b.components();
    UnitQuaternion::new(
        w * s - x * u - y * v - z * t,
        w * u + x * s + y * t - z * v,
        w * v - x * t + y * s + z * u,
        w * t + x * v - y * u + z * s,
    )
    .unwrap()
}

#[test]
fn g3_world_frames_and_common_boosts_preserve_native_pair_response() {
    let (p, s, arm, anchor, n) = full_fixture();
    let c = joint(p, arm, anchor, n);
    let h = 0.003;
    let base = c
        .advance(s, free(), [Wrench::ZERO; 2], h, controls(), || false)
        .unwrap();
    let q = UnitQuaternion::from_axis_angle(Vec3::new(1.0, -2.0, 3.0), 0.7).unwrap();
    let shift = Vec3::new(2.0, -3.0, 1.0);
    let boost = Vec3::new(9.0, -5.0, 3.0);
    let transformed = core::array::from_fn(|i| {
        RigidBodyState::new(
            Pose::new(
                q.rotate_body_to_world(s[i].pose().position_world())
                    .add(shift),
                compose(q, s[i].pose().orientation()),
            )
            .unwrap(),
            q.rotate_body_to_world(s[i].linear_momentum_world())
                .add(boost.scale(p[i].mass())),
            s[i].angular_momentum_body(),
        )
        .unwrap()
    });
    let out = c
        .advance(
            transformed,
            free(),
            [Wrench::ZERO; 2],
            h,
            controls(),
            || false,
        )
        .unwrap();
    for i in 0..2 {
        vector(
            out.states_after[i].pose().position_world(),
            q.rotate_body_to_world(base.states_after[i].pose().position_world())
                .add(shift)
                .add(boost.scale(h)),
            2e-12,
        );
        vector(
            out.states_after[i].linear_momentum_world(),
            q.rotate_body_to_world(base.states_after[i].linear_momentum_world())
                .add(boost.scale(p[i].mass())),
            2e-11,
        );
        vector(
            out.states_after[i].angular_momentum_body(),
            base.states_after[i].angular_momentum_body(),
            2e-11,
        );
    }
    near(out.constraint_work_j, base.constraint_work_j, 2e-13);
}

#[test]
fn g3_rotated_infinite_plane_does_not_depend_on_its_tangential_anchor() {
    let p = [inertia(1.0, 0.7), inertia(2.0, 1.1)];
    let rotation = UnitQuaternion::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), 0.7).unwrap();
    let normal = Vec3::new(0.0, 1.0, 0.0);
    let s = [
        state(
            p[0],
            Vec3::ZERO,
            UnitQuaternion::IDENTITY,
            Vec3::ZERO,
            Vec3::ZERO,
        ),
        state(p[1], Vec3::ZERO, rotation, Vec3::ZERO, Vec3::ZERO),
    ];
    let forces = [
        Wrench {
            force_world: rotation.rotate_body_to_world(normal),
            torque_body: Vec3::ZERO,
        },
        Wrench::ZERO,
    ];
    let reference = joint(p, Vec3::ZERO, Vec3::ZERO, normal)
        .advance(s, free(), forces, 0.01, controls(), || false)
        .unwrap();
    for offset in [1e16, 1e100, 1e300] {
        let actual = joint(p, Vec3::ZERO, Vec3::new(offset, 0.0, 0.0), normal)
            .advance(s, free(), forces, 0.01, controls(), || false)
            .unwrap();
        assert_eq!(actual, reference);
    }
}

#[test]
fn g4_hidden_rotating_carrier_speed_is_not_erased_by_world_origin_or_boost() {
    let p = [inertia(1.0, 1.0), inertia(1.0, 1.0)];
    let arm = Vec3::new(1.0, 0.0, 0.0);
    let c = joint(p, arm, Vec3::ZERO, Vec3::new(0.0, 1.0, 0.0));
    for x in [0.0, 2.0_f64.powi(54)] {
        let pos = Vec3::new(x, 0.0, 0.0);
        let s = [
            state(p[0], pos, UnitQuaternion::IDENTITY, Vec3::ZERO, Vec3::ZERO),
            state(
                p[1],
                pos,
                UnitQuaternion::IDENTITY,
                Vec3::ZERO,
                Vec3::new(0.0, 0.0, 0.2),
            ),
        ];
        assert_eq!(
            c.advance(s, free(), [Wrench::ZERO; 2], 0.01, controls(), || false),
            Err(DistanceStepError::InitialConstraint)
        );
    }
    let boost = Vec3::new(0.0, 2.0_f64.powi(54), 0.0);
    let s = [
        state(
            p[0],
            Vec3::ZERO,
            UnitQuaternion::IDENTITY,
            boost,
            Vec3::new(0.0, 0.0, 0.1),
        ),
        state(
            p[1],
            Vec3::ZERO,
            UnitQuaternion::IDENTITY,
            boost,
            Vec3::new(0.0, 0.0, 0.2),
        ),
    ];
    assert_eq!(
        c.advance(s, free(), [Wrench::ZERO; 2], 0.001, controls(), || false),
        Err(DistanceStepError::InitialConstraint)
    );
}

#[test]
fn g4_local_limits_and_every_cancellation_boundary_preserve_both_inputs() {
    let (p, s, arm, anchor, n) = full_fixture();
    let c = joint(p, arm, anchor, n);
    let mut calls = 0;
    let good = c
        .advance(s, free(), [Wrench::ZERO; 2], 0.002, controls(), || {
            calls += 1;
            false
        })
        .unwrap();
    for at in 1..=calls {
        let mut seen = 0;
        assert_eq!(
            c.advance(s, free(), [Wrench::ZERO; 2], 0.002, controls(), || {
                seen += 1;
                seen == at
            }),
            Err(DistanceStepError::Cancelled)
        );
    }
    assert_eq!(
        c.advance(s, free(), [Wrench::ZERO; 2], 0.002, controls(), || false)
            .unwrap(),
        good
    );
    assert!(FinitePointPlaneConstraint::new(p, arm, anchor, Vec3::ZERO, 1e-12, 1e-12).is_err());
    assert!(
        c.advance(s, free(), [Wrench::ZERO; 2], f64::NAN, controls(), || false)
            .is_err()
    );
    let mut cap = controls();
    cap.maximum_impulse_ns = 1e-14;
    assert!(
        c.advance(s, free(), [Wrench::ZERO; 2], 0.01, cap, || false)
            .is_err()
    );
    cap = controls();
    cap.maximum_rotation_increment_rad = 1e-8;
    assert_eq!(
        c.advance(s, free(), [Wrench::ZERO; 2], 0.01, cap, || false),
        Err(DistanceStepError::LocalBranch)
    );
}
