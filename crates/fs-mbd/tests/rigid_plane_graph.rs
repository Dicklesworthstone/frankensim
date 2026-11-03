//! Shared finite carriers, moving reaction arms and loaded sliding mechanisms.

use fs_mbd::rigid_rattle::distance::graph::{
    ArticulatedJoint, DistanceDrive, DistanceGraphControls, DistanceGraphError, DistanceLink,
    DrivenJoint, JointReactionImpulse, PointPlaneJoint, advance_articulated_graph,
    advance_spatial_driven_joint_graph,
};
use fs_mbd::rigid_rattle::distance::{DistanceStepControls, DistanceStepError};
use fs_mbd::rigid_rattle::finite_plane::FinitePointPlaneConstraint;
use fs_mbd::spatial_loads::{PreparedSpatialLoads, WorldPointForce};
use fs_mbd::{
    Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState, UnitQuaternion, Vec3,
    Wrench,
};

const X: Vec3 = Vec3::new(1., 0., 0.);
const Y: Vec3 = Vec3::new(0., 1., 0.);
const Z: Vec3 = Vec3::new(0., 0., 1.);
fn near(a: f64, b: f64, t: f64) {
    assert!((a - b).abs() <= t, "{a:e} != {b:e}, tol={t:e}");
}
fn norm(v: Vec3) -> f64 {
    v.x.hypot(v.y).hypot(v.z)
}
fn vector(a: Vec3, b: Vec3, t: f64) {
    assert!(norm(a.sub(b)) <= t, "{a:?} != {b:?}, tol={t:e}");
}
fn engine() -> RigidBodyIntegrator {
    RigidBodyIntegrator::new(Gravity::ZERO)
}
fn controls() -> DistanceGraphControls {
    DistanceGraphControls {
        step: DistanceStepControls {
            maximum_impulse_ns: 20.,
            maximum_rotation_increment_rad: 0.25,
            maximum_position_iterations: 50,
            maximum_backtracks: 20,
        },
        position_tolerance_m: 5e-13,
        velocity_tolerance_mps: 5e-13,
        minimum_scaled_pivot: 1e-10,
        maximum_bodies: 64,
        maximum_constraints: 32,
        maximum_trial_evaluations: 200,
        maximum_live_bytes: 2_000_000,
    }
}
fn state(p: MassProperties, x: Vec3, q: UnitQuaternion, v: Vec3, w: Vec3) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(x, q).unwrap(),
        v.scale(p.mass()),
        p.angular_momentum_body_checked(w).unwrap(),
    )
    .unwrap()
}
fn compare(a: RigidBodyState, b: RigidBodyState, t: f64) {
    vector(a.pose().position_world(), b.pose().position_world(), t);
    vector(a.linear_momentum_world(), b.linear_momentum_world(), t);
    vector(a.angular_momentum_body(), b.angular_momentum_body(), t);
    for axis in [X, Y, Z] {
        vector(
            a.pose().orientation().rotate_body_to_world(axis),
            b.pose().orientation().rotate_body_to_world(axis),
            t,
        );
    }
}
fn total(s: &[RigidBodyState], p: &[MassProperties]) -> (Vec3, Vec3, f64) {
    s.iter()
        .zip(p)
        .fold((Vec3::ZERO, Vec3::ZERO, 0.), |(a, b, c), (&s, &p)| {
            let d = engine().diagnostics(s, p).unwrap();
            (
                a.add(d.linear_momentum_world),
                b.add(d.angular_momentum_world),
                c + d.mechanical_energy,
            )
        })
}
fn planes() -> [ArticulatedJoint; 2] {
    [0, 1].map(|i| {
        ArticulatedJoint::PointPlane(
            PointPlaneJoint::new([i, 2], Vec3::ZERO, Vec3::ZERO, Z).unwrap(),
        )
    })
}
fn properties() -> [MassProperties; 3] {
    [(2., 1.), (3., 1.5), (5., 4.)]
        .map(|(m, i)| MassProperties::new(m, Vec3::ZERO, Vec3::new(i, i, i)).unwrap())
}

#[test]
fn g1_shared_carrier_reaction_matches_independent_coupled_acceleration_and_original_pair() {
    let p = properties();
    let x = [-0.7, 1.1];
    let s = [x[0], x[1], 0.].map(|x| Pose::new(X.scale(x), UnitQuaternion::IDENTITY).unwrap());
    let s = core::array::from_fn::<_, 3, _>(|i| {
        RigidBodyState::new(s[i], Vec3::ZERO, Vec3::ZERO).unwrap()
    });
    let f = [1.3, -0.4, 0.3];
    let torque = 0.2;
    let loads = core::array::from_fn::<_, 3, _>(|i| Wrench {
        force_world: Z.scale(f[i]),
        torque_body: if i == 2 { Y.scale(torque) } else { Vec3::ZERO },
    });
    // Direct acceleration balance, independent of the native impulse solver:
    // K_ij=delta_ij/m_i+1/M+x_i*x_j/I, K*lambda=-free_gap_acceleration.
    let k00 = 1. / p[0].mass() + 0.2 + x[0] * x[0] / 4.;
    let k11 = 1. / p[1].mass() + 0.2 + x[1] * x[1] / 4.;
    let k01 = 0.2 + x[0] * x[1] / 4.;
    let rhs =
        core::array::from_fn::<_, 2, _>(|i| -(f[i] / p[i].mass() - f[2] / 5. + x[i] * torque / 4.));
    let det = k00 * k11 - k01 * k01;
    let lambda = [
        (k11 * rhs[0] - k01 * rhs[1]) / det,
        (k00 * rhs[1] - k01 * rhs[0]) / det,
    ];
    let mut errors = Vec::new();
    for h in [0.04, 0.02, 0.01] {
        let out =
            advance_articulated_graph(&s, &p, &planes(), engine(), &loads, h, controls(), || false)
                .unwrap();
        let mut error = 0.;
        for i in 0..2 {
            let impulse = [
                &out.step.position_impulses[i],
                &out.step.velocity_impulses[i],
            ]
            .into_iter()
            .map(|r| match r {
                JointReactionImpulse::Point(r) => r.first.impulse_world.z,
                _ => unreachable!(),
            })
            .sum::<f64>();
            error += (impulse / h - lambda[i]).abs();
        }
        let end = out.step.states_after[2];
        error +=
            (end.linear_momentum_world().z / (5. * h) - (f[2] - lambda[0] - lambda[1]) / 5.).abs();
        error += (end.angular_momentum_body().y / (4. * h)
            - (torque + x[0] * lambda[0] + x[1] * lambda[1]) / 4.)
            .abs();
        errors.push(error);
        let after = total(&out.step.states_after, &p);
        vector(after.0, Z.scale(h * f.iter().sum::<f64>()), 3e-14);
    }
    assert!(errors[2] < 2e-8, "{errors:?}");
    for e in errors.windows(2) {
        // Starting from rest, the contact arms first change at order h^2.
        assert!(e[0] > 3.5 * e[1], "{errors:?}");
    }
    let pair = FinitePointPlaneConstraint::new(
        [p[0], p[2]],
        Vec3::ZERO,
        Vec3::ZERO,
        Z,
        controls().position_tolerance_m,
        controls().velocity_tolerance_mps,
    )
    .unwrap()
    .advance(
        [s[0], s[2]],
        engine(),
        [loads[0], loads[2]],
        0.03,
        controls().step,
        || false,
    )
    .unwrap();
    let one = [ArticulatedJoint::PointPlane(
        PointPlaneJoint::new([0, 1], Vec3::ZERO, Vec3::ZERO, Z).unwrap(),
    )];
    let graph = advance_articulated_graph(
        &[s[0], s[2]],
        &[p[0], p[2]],
        &one,
        engine(),
        &[loads[0], loads[2]],
        0.03,
        controls(),
        || false,
    )
    .unwrap();
    for i in 0..2 {
        compare(pair.states_after[i], graph.step.states_after[i], 2e-13);
    }
    near(pair.constraint_work_j, graph.constraint_work_j, 2e-14);
}

fn radial_states(y: [f64; 6], p: &[MassProperties; 3]) -> [RigidBodyState; 3] {
    let [r0, r1, u0, u1, theta, omega] = y;
    let q = UnitQuaternion::from_axis_angle(Y, theta).unwrap();
    let t = q.rotate_body_to_world(X);
    let dt = Y.cross(t).scale(omega);
    let center = (2. * r0 + 3. * r1) / 10.;
    let speed = (2. * u0 + 3. * u1) / 10.;
    let r = [r0 - center, r1 - center, -center];
    let u = [u0 - speed, u1 - speed, -speed];
    core::array::from_fn(|i| {
        state(
            p[i],
            t.scale(r[i]),
            if i == 2 { q } else { UnitQuaternion::IDENTITY },
            t.scale(u[i]).add(dt.scale(r[i])),
            if i == 2 { Y.scale(omega) } else { Vec3::ZERO },
        )
    })
}

#[test]
fn g1_two_sliding_masses_exchange_orbital_momentum_with_one_rotating_carrier() {
    let p = properties();
    let initial_y = [-0.6, 1.1, 0.17, -0.08, 0.2, 0.4];
    let initial = radial_states(initial_y, &p);
    // Eliminate the free system COM: B_ij=m_i*delta_ij-m_i*m_j/M_total.
    // The continuous Lagrangian gives r''=r*omega^2 and
    // omega'=-2*(r_dot^T B r)*omega/(I_carrier+r^T B r).
    let rhs = |y: [f64; 6]| {
        let br = [1.6 * y[0] - 0.6 * y[1], -0.6 * y[0] + 2.1 * y[1]];
        [
            y[2],
            y[3],
            y[0] * y[5] * y[5],
            y[1] * y[5] * y[5],
            y[5],
            -2. * (y[2] * br[0] + y[3] * br[1]) * y[5] / (4. + y[0] * br[0] + y[1] * br[1]),
        ]
    };
    let time = 0.8;
    let mut y = initial_y;
    let h = time / 8192.;
    for _ in 0..8192 {
        let a = rhs(y);
        let b = rhs(core::array::from_fn(|i| y[i] + h * a[i] / 2.));
        let c = rhs(core::array::from_fn(|i| y[i] + h * b[i] / 2.));
        let d = rhs(core::array::from_fn(|i| y[i] + h * c[i]));
        y = core::array::from_fn(|i| y[i] + h * (a[i] + 2. * b[i] + 2. * c[i] + d[i]) / 6.);
    }
    let expected = radial_states(y, &p);
    let before = total(&initial, &p);
    let mut errors = Vec::new();
    for count in [16, 32, 64] {
        let h = time / f64::from(count);
        let mut s = initial.to_vec();
        let mut work = 0.;
        for _ in 0..count {
            let out = advance_articulated_graph(
                &s,
                &p,
                &planes(),
                engine(),
                &[Wrench::ZERO; 3],
                h,
                controls(),
                || false,
            )
            .unwrap();
            work += out.constraint_work_j;
            s = out.step.states_after;
        }
        let error = s
            .iter()
            .zip(expected)
            .map(|(a, b)| {
                norm(a.pose().position_world().sub(b.pose().position_world()))
                    + norm(a.linear_momentum_world().sub(b.linear_momentum_world()))
                    + norm(a.angular_momentum_body().sub(b.angular_momentum_body()))
            })
            .sum::<f64>();
        errors.push(error);
        let after = total(&s, &p);
        vector(after.0, before.0, 2e-12);
        vector(after.1, before.1, 3e-11);
        near(after.2 - before.2, work, 3e-12);
        s = s
            .into_iter()
            .map(|s| {
                RigidBodyState::new(
                    s.pose(),
                    s.linear_momentum_world().scale(-1.),
                    s.angular_momentum_body().scale(-1.),
                )
                .unwrap()
            })
            .collect();
        for _ in 0..count {
            s = advance_articulated_graph(
                &s,
                &p,
                &planes(),
                engine(),
                &[Wrench::ZERO; 3],
                h,
                controls(),
                || false,
            )
            .unwrap()
            .step
            .states_after;
        }
        for i in 0..3 {
            let back = RigidBodyState::new(
                s[i].pose(),
                s[i].linear_momentum_world().scale(-1.),
                s[i].angular_momentum_body().scale(-1.),
            )
            .unwrap();
            compare(back, initial[i], 2e-10);
        }
    }
    assert!(errors[2] < 3e-5, "{errors:?}");
    for e in errors.windows(2) {
        assert!((3.7..4.3).contains(&(e[0] / e[1])), "{errors:?}");
    }
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
fn mixed(
    rotation: UnitQuaternion,
    boost: Vec3,
) -> (
    [MassProperties; 3],
    [RigidBodyState; 3],
    [DrivenJoint; 3],
    [PreparedSpatialLoads; 3],
) {
    let p = core::array::from_fn(|i| {
        MassProperties::from_principal_frame(
            2. + i as f64,
            Vec3::ZERO,
            Vec3::new(0.8, 1.2, 1.5).scale(1. + 0.2 * i as f64),
            UnitQuaternion::from_axis_angle(Vec3::new(1., 2., -1.), 0.3 + 0.2 * i as f64).unwrap(),
        )
        .unwrap()
    });
    let q = core::array::from_fn::<_, 3, _>(|i| {
        UnitQuaternion::from_axis_angle(Vec3::new(1., -0.5, 2.), 0.2 - 0.3 * i as f64).unwrap()
    });
    let arms = [Vec3::new(0.2, -0.1, 0.15), Vec3::new(-0.1, 0.15, 0.2)];
    let anchor = Z.scale(0.15);
    let xc = Vec3::new(0.2, -0.1, 0.3);
    let contact = [Vec3::new(-0.7, 0.3, 0.15), Vec3::new(0.8, -0.2, 0.15)];
    let x = [
        xc.add(q[2].rotate_body_to_world(contact[0]))
            .sub(q[0].rotate_body_to_world(arms[0])),
        xc.add(q[2].rotate_body_to_world(contact[1]))
            .sub(q[1].rotate_body_to_world(arms[1])),
        xc,
    ];
    let omega = Vec3::new(0.2, -0.1, 0.3);
    let s = core::array::from_fn(|i| {
        state(
            p[i],
            rotation.rotate_body_to_world(x[i]),
            compose(rotation, q[i]),
            rotation.rotate_body_to_world(omega.cross(x[i])).add(boost),
            q[i].rotate_world_to_body(omega),
        )
    });
    let joints = [
        DrivenJoint::PointPlane(PointPlaneJoint::new([0, 2], arms[0], anchor, Z).unwrap()),
        DrivenJoint::PointPlane(PointPlaneJoint::new([1, 2], arms[1], anchor, Z).unwrap()),
        DrivenJoint::Distance {
            link: DistanceLink {
                bodies: [0, 1],
                points_body_m: [Vec3::ZERO; 2],
                length_m: norm(x[1].sub(x[0])),
            },
            drive: DistanceDrive {
                range_rate_mps: 0.,
                range_acceleration_mps2: 0.13,
            },
        },
    ];
    let loads = core::array::from_fn(|i| {
        let factor = 0.2 + 0.1 * i as f64;
        PreparedSpatialLoads::new(
            Wrench {
                force_world: rotation.rotate_body_to_world(Vec3::new(0.3, -0.2, 0.1).scale(factor)),
                torque_body: Vec3::new(0.02, -0.03, 0.01).scale(factor),
            },
            rotation.rotate_body_to_world(Vec3::new(-0.01, 0.02, 0.03).scale(factor)),
            &[WorldPointForce {
                arm_body: Vec3::new(0.4, -0.2, 0.3),
                force_world: rotation
                    .rotate_body_to_world(Vec3::new(-0.1, 0.3, -0.2).scale(factor)),
            }],
            || false,
        )
        .unwrap()
    });
    (p, s, joints, loads)
}

#[test]
fn g0_g3_g4_off_center_sliding_composes_with_full_tensors_spatial_loads_and_actuation() {
    let (p, initial, original, loads) = mixed(UnitQuaternion::IDENTITY, Vec3::ZERO);
    let before = total(&initial, &p);
    let mut s = initial.to_vec();
    let mut joints = original.to_vec();
    let (mut linear, mut angular, mut work) = (Vec3::ZERO, Vec3::ZERO, 0.);
    for _ in 0..64 {
        let out = advance_spatial_driven_joint_graph(
            &s,
            &p,
            &joints,
            engine(),
            &loads,
            0.004,
            controls(),
            || false,
        )
        .unwrap();
        work += out.actuator_work_j + out.constraint_work_defect_j;
        angular = angular.add(out.constraint_angular_impulse_world_nms);
        assert_eq!(out.joints[0].actuator_work_j, 0.);
        assert_eq!(out.joints[1].actuator_work_j, 0.);
        for free in &out.step.free_steps {
            linear = linear.add(free.linear_impulse_world);
            angular = angular.add(free.angular_impulse_world);
            work += free.force_work_j + free.couple_work_j + free.work_energy_residual_j;
        }
        joints = out.joints.iter().map(|j| j.joint_after).collect();
        s = out.step.states_after;
    }
    let after = total(&s, &p);
    vector(after.0.sub(before.0), linear, 3e-12);
    vector(after.1.sub(before.1), angular, 3e-11);
    near(after.2 - before.2, work, 3e-12);

    let mut polls = 0;
    let out = advance_spatial_driven_joint_graph(
        &initial,
        &p,
        &original,
        engine(),
        &loads,
        0.004,
        controls(),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let rotation = UnitQuaternion::from_axis_angle(Vec3::new(1., 2., -1.), 0.8).unwrap();
    let boost = Vec3::new(0.7, -0.2, 0.3);
    let (rp, rs, rj, rl) = mixed(rotation, boost);
    let moved =
        advance_spatial_driven_joint_graph(&rs, &rp, &rj, engine(), &rl, 0.004, controls(), || {
            false
        })
        .unwrap();
    for i in 0..3 {
        let base = out.step.states_after[i];
        let expected = RigidBodyState::new(
            Pose::new(
                rotation
                    .rotate_body_to_world(base.pose().position_world())
                    .add(boost.scale(0.004)),
                compose(rotation, base.pose().orientation()),
            )
            .unwrap(),
            rotation
                .rotate_body_to_world(base.linear_momentum_world())
                .add(boost.scale(p[i].mass())),
            base.angular_momentum_body(),
        )
        .unwrap();
        compare(moved.step.states_after[i], expected, 3e-11);
    }
    let permuted = advance_spatial_driven_joint_graph(
        &initial,
        &p,
        &[original[2], original[1], original[0]],
        engine(),
        &loads,
        0.004,
        controls(),
        || false,
    )
    .unwrap();
    for i in 0..3 {
        compare(
            out.step.states_after[i],
            permuted.step.states_after[i],
            3e-12,
        );
    }
    let mut exact = controls();
    exact.maximum_live_bytes = out.planned_live_bytes;
    assert_eq!(
        out,
        advance_spatial_driven_joint_graph(
            &initial,
            &p,
            &original,
            engine(),
            &loads,
            0.004,
            exact,
            || false
        )
        .unwrap()
    );
    exact.maximum_live_bytes -= 1;
    assert!(matches!(
        advance_spatial_driven_joint_graph(
            &initial,
            &p,
            &original,
            engine(),
            &loads,
            0.004,
            exact,
            || false
        ),
        Err(DistanceGraphError::Limit { .. })
    ));
    for stop in 1..=polls {
        let mut n = 0;
        let refused = advance_spatial_driven_joint_graph(
            &initial,
            &p,
            &original,
            engine(),
            &loads,
            0.004,
            controls(),
            || {
                n += 1;
                n == stop
            },
        );
        assert!(
            matches!(
                refused,
                Err(DistanceGraphError::Step(DistanceStepError::Cancelled))
            ),
            "{stop}/{polls}: {refused:?}"
        );
    }
    assert_eq!(
        out,
        advance_spatial_driven_joint_graph(
            &initial,
            &p,
            &original,
            engine(),
            &loads,
            0.004,
            controls(),
            || false
        )
        .unwrap()
    );
    assert!(PointPlaneJoint::new([0, 0], Vec3::ZERO, Vec3::ZERO, Z).is_err());
    assert!(PointPlaneJoint::new([0, 1], Vec3::ZERO, Vec3::ZERO, Vec3::ZERO).is_err());
}
