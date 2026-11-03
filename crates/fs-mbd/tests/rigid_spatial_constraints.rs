//! Physical source moments, finite recoil, and independent continuous references.

use fs_mbd::rigid_rattle::distance::{
    DistanceStepControls, DistanceStepError, RigidDistanceConstraint,
};
use fs_mbd::rigid_rattle::finite_plane::FinitePointPlaneConstraint;
use fs_mbd::spatial_loads::{PreparedSpatialLoads, WorldPointForce};
use fs_mbd::{
    Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState, UnitQuaternion, Vec3,
    Wrench,
};

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}; tol={tolerance:e}"
    );
}
fn vector(a: Vec3, b: Vec3, tolerance: f64) {
    near(a.x, b.x, tolerance);
    near(a.y, b.y, tolerance);
    near(a.z, b.z, tolerance);
}
fn controls() -> DistanceStepControls {
    DistanceStepControls {
        maximum_impulse_ns: 100.,
        maximum_rotation_increment_rad: 0.25,
        maximum_position_iterations: 60,
        maximum_backtracks: 20,
    }
}
fn engine() -> RigidBodyIntegrator {
    RigidBodyIntegrator::new(Gravity::ZERO)
}
fn inertia(mass: f64, moment: f64) -> MassProperties {
    MassProperties::new(mass, Vec3::ZERO, Vec3::new(moment, moment, moment)).unwrap()
}
fn state(p: MassProperties, x: Vec3, q: UnitQuaternion, v: Vec3, w: Vec3) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(x, q).unwrap(),
        v.scale(p.mass()),
        p.angular_momentum_body_checked(w).unwrap(),
    )
    .unwrap()
}
fn angle(s: RigidBodyState) -> f64 {
    let x = s
        .pose()
        .orientation()
        .rotate_body_to_world(Vec3::new(1., 0., 0.));
    x.y.atan2(x.x)
}
fn rotation(theta: f64) -> UnitQuaternion {
    UnitQuaternion::from_axis_angle(Vec3::new(0., 0., 1.), theta).unwrap()
}
fn force_couple(arm: f64, force: f64) -> PreparedSpatialLoads {
    PreparedSpatialLoads::new(
        Wrench::ZERO,
        Vec3::ZERO,
        &[
            WorldPointForce {
                arm_body: Vec3::new(arm, 0., 0.),
                force_world: Vec3::new(0., force, 0.),
            },
            WorldPointForce {
                arm_body: Vec3::new(-arm, 0., 0.),
                force_world: Vec3::new(0., -force, 0.),
            },
        ],
        || false,
    )
    .unwrap()
}
fn inventories(
    s: [RigidBodyState; 2],
    p: [MassProperties; 2],
    e: RigidBodyIntegrator,
) -> (Vec3, Vec3, f64) {
    let a = e.diagnostics(s[0], p[0]).unwrap();
    let b = e.diagnostics(s[1], p[1]).unwrap();
    (
        a.linear_momentum_world.add(b.linear_momentum_world),
        a.angular_momentum_world.add(b.angular_momentum_world),
        a.mechanical_energy + b.mechanical_energy,
    )
}
fn reverse(s: [RigidBodyState; 2]) -> [RigidBodyState; 2] {
    s.map(|s| {
        RigidBodyState::new(
            s.pose(),
            s.linear_momentum_world().scale(-1.),
            s.angular_momentum_body().scale(-1.),
        )
        .unwrap()
    })
}

// Independent scalar/vector continuous equations, never the discrete free map.
fn reference<const N: usize>(
    mut y: [f64; N],
    time: f64,
    rhs: impl Fn([f64; N]) -> [f64; N],
) -> [f64; N] {
    let h = time / 8192.;
    for _ in 0..8192 {
        let a = rhs(y);
        let b = rhs(core::array::from_fn(|i| y[i] + 0.5 * h * a[i]));
        let c = rhs(core::array::from_fn(|i| y[i] + 0.5 * h * b[i]));
        let d = rhs(core::array::from_fn(|i| y[i] + h * c[i]));
        y = core::array::from_fn(|i| y[i] + h * (a[i] + 2. * b[i] + 2. * c[i] + d[i]) / 6.);
    }
    y
}

#[test]
fn g1_rotating_point_loads_on_a_constrained_pair_match_independent_rotor() {
    let p = [inertia(1.5, 0.7), inertia(2.2, 1.1)];
    let total = p[0].mass() + p[1].mass();
    let r = Vec3::new(0.8, 0., 0.);
    let v = Vec3::new(0., 0.56, 0.);
    let initial = [
        state(
            p[0],
            r.scale(-p[1].mass() / total),
            rotation(0.4),
            v.scale(-p[1].mass() / total),
            Vec3::new(0., 0., 0.6),
        ),
        state(
            p[1],
            r.scale(p[0].mass() / total),
            UnitQuaternion::IDENTITY,
            v.scale(p[0].mass() / total),
            Vec3::ZERO,
        ),
    ];
    let joint = RigidDistanceConstraint::new(p, [Vec3::ZERO; 2], 0.8, 1e-12, 1e-12).unwrap();
    let loads = [force_couple(0.4, 0.8), PreparedSpatialLoads::ZERO];
    let time = 0.8;
    let expected = reference([0.4, 0.6], time, |y| [y[1], (0.64 / 0.7) * y[0].cos()]);
    let before = inventories(initial, p, engine());
    let mut errors = Vec::new();
    let mut finest = initial;
    for count in [16, 32, 64] {
        let h = time / f64::from(count);
        let mut s = initial;
        let (mut torque, mut work, mut residual) = (Vec3::ZERO, 0., 0.);
        for _ in 0..count {
            let out = joint
                .advance_spatial(s, engine(), loads, h, controls(), || false)
                .unwrap();
            work += out.constraint_work_j;
            for (i, step) in out.free_steps.into_iter().enumerate() {
                torque = torque.add(step.angular_impulse_world);
                work += step.force_work_j + step.couple_work_j;
                residual += step.work_energy_residual_j;
                assert_eq!(
                    step,
                    engine()
                        .step_spatial(step.step.state_before, p[i], &loads[i], h)
                        .unwrap()
                );
            }
            s = out.states_after;
        }
        let after = inventories(s, p, engine());
        vector(after.0, before.0, 2e-12);
        vector(after.1.sub(before.1), torque, 2e-11);
        near(after.2 - before.2, work + residual, 2e-12);
        errors.push(
            (angle(s[0]) - expected[0]).abs()
                + (s[0].angular_momentum_body().z / 0.7 - expected[1]).abs(),
        );
        finest = s;
        let mut reversed = reverse(s);
        for _ in 0..count {
            reversed = joint
                .advance_spatial(reversed, engine(), loads, h, controls(), || false)
                .unwrap()
                .states_after;
        }
        for i in 0..2 {
            vector(
                reversed[i].pose().position_world(),
                initial[i].pose().position_world(),
                2e-9,
            );
            vector(
                reversed[i].linear_momentum_world(),
                initial[i].linear_momentum_world().scale(-1.),
                2e-9,
            );
            vector(
                reversed[i].angular_momentum_body(),
                initial[i].angular_momentum_body().scale(-1.),
                2e-10,
            );
        }
    }
    assert!(
        errors[0] > 3.8 * errors[1] && errors[1] > 3.8 * errors[2],
        "{errors:?}"
    );
    let mut frozen = initial;
    for _ in 0..64 {
        let wrenches =
            core::array::from_fn(|i| loads[i].wrench_at(frozen[i].pose().orientation()).unwrap());
        frozen = joint
            .advance(frozen, engine(), wrenches, time / 64., controls(), || false)
            .unwrap()
            .states_after;
    }
    assert!((angle(frozen[0]) - expected[0]).abs() > 10. * errors[2]);
    assert!(angle(finest[0]) > 0.4);
}

#[test]
fn g1_loaded_sliding_carrier_matches_independent_coupled_rotor() {
    let p = [inertia(1.2, 0.5), inertia(2.3, 0.8)];
    let total = p[0].mass() + p[1].mass();
    let mu = p[0].mass() * p[1].mass() / total;
    let r = Vec3::new(0.7, 0., 0.);
    let v = Vec3::new(0.2, 0.7 * 1.1, 0.);
    let initial = [
        state(
            p[0],
            r.scale(p[1].mass() / total),
            UnitQuaternion::IDENTITY,
            v.scale(p[1].mass() / total),
            Vec3::ZERO,
        ),
        state(
            p[1],
            r.scale(-p[0].mass() / total),
            UnitQuaternion::IDENTITY,
            v.scale(-p[0].mass() / total),
            Vec3::new(0., 0., 1.1),
        ),
    ];
    let joint = FinitePointPlaneConstraint::new(
        p,
        Vec3::ZERO,
        Vec3::ZERO,
        Vec3::new(0., 1., 0.),
        1e-12,
        1e-12,
    )
    .unwrap();
    let loads = [PreparedSpatialLoads::ZERO, force_couple(0.3, 0.5)];
    let time = 0.6;
    let expected = reference([0.7, 0.2, 0., 1.1], time, |y| {
        [
            y[1],
            y[0] * y[3] * y[3],
            y[3],
            (0.3 * y[2].cos() - 2. * mu * y[0] * y[1] * y[3]) / (0.8 + mu * y[0] * y[0]),
        ]
    });
    let expected_position = Vec3::new(
        expected[0] * expected[2].cos(),
        expected[0] * expected[2].sin(),
        0.,
    );
    let mut errors = Vec::new();
    for count in [24, 48, 96] {
        let h = time / f64::from(count);
        let mut s = initial;
        for _ in 0..count {
            s = joint
                .advance_spatial(s, engine(), loads, h, controls(), || false)
                .unwrap()
                .states_after;
        }
        errors.push(
            s[0].pose()
                .position_world()
                .sub(s[1].pose().position_world())
                .sub(expected_position)
                .norm_squared()
                .sqrt()
                + (s[1].angular_momentum_body().z / 0.8 - expected[3]).abs(),
        );
        let mut reversed = reverse(s);
        for _ in 0..count {
            reversed = joint
                .advance_spatial(reversed, engine(), loads, h, controls(), || false)
                .unwrap()
                .states_after;
        }
        for i in 0..2 {
            vector(
                reversed[i].pose().position_world(),
                initial[i].pose().position_world(),
                3e-9,
            );
            vector(
                reversed[i].angular_momentum_body(),
                initial[i].angular_momentum_body().scale(-1.),
                3e-9,
            );
        }
    }
    assert!(
        errors[0] > 3.8 * errors[1] && errors[1] > 3.8 * errors[2],
        "{errors:?}"
    );
}

fn full_fixture() -> ([MassProperties; 2], [RigidBodyState; 2], Vec3, Vec3) {
    let p = [
        MassProperties::from_principal_frame(
            1.7,
            Vec3::ZERO,
            Vec3::new(0.7, 1., 1.2),
            UnitQuaternion::from_axis_angle(Vec3::new(1., 2., -1.), 0.6).unwrap(),
        )
        .unwrap(),
        MassProperties::from_principal_frame(
            2.2,
            Vec3::ZERO,
            Vec3::new(0.9, 1.1, 1.5),
            UnitQuaternion::from_axis_angle(Vec3::new(-2., 1., 3.), 0.4).unwrap(),
        )
        .unwrap(),
    ];
    let qa = UnitQuaternion::from_axis_angle(Vec3::new(0.4, -0.7, 1.), 0.5).unwrap();
    let qb = UnitQuaternion::from_axis_angle(Vec3::new(-0.3, 0.2, 0.8), 0.4).unwrap();
    let arm = Vec3::new(0.5, -0.3, 0.2);
    let anchor = Vec3::new(0.2, 0.3, -0.1);
    let wb = Vec3::new(0.2, -0.1, 0.4);
    let wa = Vec3::new(-0.3, 0.2, 0.1);
    let xb = Vec3::new(0.1, -0.4, 0.2);
    let r = qb.rotate_body_to_world(anchor.add(Vec3::new(0.6, 0., -0.2)));
    let xa = xb.add(r).sub(qa.rotate_body_to_world(arm));
    let vb = Vec3::new(0.1, 0.3, -0.2);
    let va = vb
        .add(qb.rotate_body_to_world(wb).cross(r))
        .sub(
            qa.rotate_body_to_world(wa)
                .cross(qa.rotate_body_to_world(arm)),
        )
        .add(qb.rotate_body_to_world(Vec3::new(0.3, 0., 0.1)));
    (
        p,
        [state(p[0], xa, qa, va, wa), state(p[1], xb, qb, vb, wb)],
        arm,
        anchor,
    )
}
fn full_loads(q: UnitQuaternion) -> [PreparedSpatialLoads; 2] {
    core::array::from_fn(|i| {
        let s = if i == 0 { 1. } else { -0.7 };
        PreparedSpatialLoads::new(
            Wrench {
                force_world: q.rotate_body_to_world(Vec3::new(0.13, -0.08, 0.21).scale(s)),
                torque_body: Vec3::new(0.03, -0.04, 0.02).scale(s),
            },
            q.rotate_body_to_world(Vec3::new(-0.05, 0.08, 0.02).scale(s)),
            &[WorldPointForce {
                arm_body: Vec3::new(0.4, -0.3, 0.2),
                force_world: q.rotate_body_to_world(Vec3::new(0.3, 0.1, -0.2).scale(s)),
            }],
            || false,
        )
        .unwrap()
    })
}

#[test]
fn g0_full_inertia_loaded_carrier_closes_external_and_reaction_accounts() {
    let (p, initial, arm, anchor) = full_fixture();
    let joint =
        FinitePointPlaneConstraint::new(p, arm, anchor, Vec3::new(0., 1., 0.), 1e-12, 1e-12)
            .unwrap();
    let loads = full_loads(UnitQuaternion::IDENTITY);
    let e = RigidBodyIntegrator::new(Gravity::new(Vec3::new(0.1, -0.3, -0.2)).unwrap());
    let before = inventories(initial, p, e);
    let mut s = initial;
    let (mut linear, mut angular, mut work, mut residual) = (Vec3::ZERO, Vec3::ZERO, 0., 0.);
    for _ in 0..200 {
        let out = joint
            .advance_spatial(s, e, loads, 0.002, controls(), || false)
            .unwrap();
        work += out.constraint_work_j;
        angular = angular
            .add(out.position_impulse.angular_impulse_world)
            .add(out.velocity_impulse.angular_impulse_world);
        for (i, step) in out.free_steps.into_iter().enumerate() {
            assert_eq!(
                step,
                e.step_spatial(step.step.state_before, p[i], &loads[i], 0.002)
                    .unwrap()
            );
            linear = linear.add(step.linear_impulse_world);
            angular = angular.add(step.angular_impulse_world);
            work += step.force_work_j + step.couple_work_j;
            residual += step.work_energy_residual_j;
        }
        s = out.states_after;
        let after = inventories(s, p, e);
        vector(after.0.sub(before.0), linear, 3e-12);
        vector(after.1.sub(before.1), angular, 5e-11);
        near(after.2 - before.2, work + residual, 5e-12);
    }
    assert!(angular.norm_squared() > 1e-5);
    assert!(
        s[1].angular_momentum_body()
            .sub(initial[1].angular_momentum_body())
            .norm_squared()
            > 1e-5
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
fn g3_spatially_loaded_constraint_is_covariant_under_frame_and_boost() {
    let (p, s, arm, anchor) = full_fixture();
    let joint =
        FinitePointPlaneConstraint::new(p, arm, anchor, Vec3::new(0., 1., 0.), 1e-12, 1e-12)
            .unwrap();
    let h = 0.004;
    let base = joint
        .advance_spatial(
            s,
            engine(),
            full_loads(UnitQuaternion::IDENTITY),
            h,
            controls(),
            || false,
        )
        .unwrap();
    let q = UnitQuaternion::from_axis_angle(Vec3::new(1., -2., 3.), 0.7).unwrap();
    let shift = Vec3::new(2., -3., 1.);
    let boost = Vec3::new(9., -5., 3.);
    let moved = core::array::from_fn(|i| {
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
    let out = joint
        .advance_spatial(moved, engine(), full_loads(q), h, controls(), || false)
        .unwrap();
    for i in 0..2 {
        vector(
            out.states_after[i].pose().position_world(),
            q.rotate_body_to_world(base.states_after[i].pose().position_world())
                .add(shift)
                .add(boost.scale(h)),
            3e-12,
        );
        vector(
            out.states_after[i].linear_momentum_world(),
            q.rotate_body_to_world(base.states_after[i].linear_momentum_world())
                .add(boost.scale(p[i].mass())),
            3e-11,
        );
        vector(
            out.states_after[i].angular_momentum_body(),
            base.states_after[i].angular_momentum_body(),
            3e-11,
        );
        vector(
            out.free_steps[i].linear_impulse_world,
            q.rotate_body_to_world(base.free_steps[i].linear_impulse_world),
            2e-13,
        );
        near(
            out.free_steps[i].force_work_j,
            base.free_steps[i].force_work_j
                + boost.dot(q.rotate_body_to_world(base.free_steps[i].linear_impulse_world)),
            2e-12,
        );
    }
    near(out.constraint_work_j, base.constraint_work_j, 3e-13);
}

#[test]
fn g4_spatial_constraint_bounds_and_all_cancellation_points_are_atomic() {
    let (p, s, arm, anchor) = full_fixture();
    let joint =
        FinitePointPlaneConstraint::new(p, arm, anchor, Vec3::new(0., 1., 0.), 1e-12, 1e-12)
            .unwrap();
    let loads = full_loads(UnitQuaternion::IDENTITY);
    let mut polls = 0;
    let baseline = joint
        .advance_spatial(s, engine(), loads, 0.004, controls(), || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(polls > 8);
    for stop in 0..polls {
        let mut calls = 0;
        assert_eq!(
            joint.advance_spatial(s, engine(), loads, 0.004, controls(), || {
                let halt = calls == stop;
                calls += 1;
                halt
            }),
            Err(DistanceStepError::Cancelled)
        );
    }
    assert_eq!(
        joint
            .advance_spatial(s, engine(), loads, 0.004, controls(), || false)
            .unwrap(),
        baseline
    );
    let mut limited = controls();
    limited.maximum_rotation_increment_rad = 1e-8;
    assert_eq!(
        joint.advance_spatial(s, engine(), loads, 0.004, limited, || false),
        Err(DistanceStepError::LocalBranch)
    );
    for h in [0., -0.01, f64::INFINITY, f64::NAN] {
        assert!(
            joint
                .advance_spatial(s, engine(), loads, h, controls(), || false)
                .is_err()
        );
    }
    let spring = PreparedSpatialLoads::new(
        Wrench::ZERO,
        Vec3::ZERO,
        &[
            WorldPointForce {
                arm_body: Vec3::new(1., 0., 0.),
                force_world: Vec3::new(1000., 0., 0.),
            },
            WorldPointForce {
                arm_body: Vec3::new(-1., 0., 0.),
                force_world: Vec3::new(-1000., 0., 0.),
            },
        ],
        || false,
    )
    .unwrap();
    assert_eq!(
        spring
            .wrench_at(UnitQuaternion::IDENTITY)
            .unwrap()
            .torque_body,
        Vec3::ZERO
    );
    assert!(spring.torque_norm_bound().unwrap() >= 2000.);
    for a in [0., 0.1, 1., 2., 3.] {
        let t = spring.wrench_at(rotation(a)).unwrap().torque_body;
        assert!(t.norm_squared().sqrt() <= spring.torque_norm_bound().unwrap() * (1. + 1e-15));
    }
    let p = [inertia(1., 1.), inertia(1., 1.)];
    let s = [state(
        p[0],
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
        Vec3::ZERO,
        Vec3::ZERO,
    ); 2];
    let c = FinitePointPlaneConstraint::new(
        p,
        Vec3::ZERO,
        Vec3::ZERO,
        Vec3::new(0., 1., 0.),
        1e-12,
        1e-12,
    )
    .unwrap();
    assert_eq!(
        c.advance_spatial(
            s,
            engine(),
            [spring, PreparedSpatialLoads::ZERO],
            0.02,
            controls(),
            || false
        ),
        Err(DistanceStepError::LocalBranch)
    );
}
