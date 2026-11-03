//! Loaded finite mechanisms through the native graph and spatial-force owners.

use fs_mbd::rigid_rattle::distance::graph::{
    ArticulatedJoint, DistanceDrive, DistanceGraphControls, DistanceGraphError, DistanceLink,
    DrivenJoint, DrivenJointStep, JointReactionImpulse, RevoluteJoint, SphericalJoint,
    advance_spatial_articulated_graph, advance_spatial_driven_joint_graph,
};
use fs_mbd::rigid_rattle::distance::{DistanceStepControls, DistanceStepError};
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
fn vector(a: Vec3, b: Vec3, t: f64) {
    near(a.x, b.x, t);
    near(a.y, b.y, t);
    near(a.z, b.z, t);
}
fn magnitude(v: Vec3) -> f64 {
    v.x.hypot(v.y).hypot(v.z)
}
fn free() -> RigidBodyIntegrator {
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
fn hinge(bodies: [usize; 2], points: [Vec3; 2], axes: [Vec3; 2], scale: f64) -> ArticulatedJoint {
    ArticulatedJoint::Revolute(RevoluteJoint::new(bodies, points, axes, scale).unwrap())
}
fn total(s: &[RigidBodyState], p: &[MassProperties], e: RigidBodyIntegrator) -> (Vec3, Vec3, f64) {
    s.iter()
        .zip(p)
        .fold((Vec3::ZERO, Vec3::ZERO, 0.), |(a, b, c), (&s, &p)| {
            let d = e.diagnostics(s, p).unwrap();
            (
                a.add(d.linear_momentum_world),
                b.add(d.angular_momentum_world),
                c + d.mechanical_energy,
            )
        })
}
fn reverse(s: RigidBodyState) -> RigidBodyState {
    RigidBodyState::new(
        s.pose(),
        s.linear_momentum_world().scale(-1.),
        s.angular_momentum_body().scale(-1.),
    )
    .unwrap()
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
fn compare_state(a: RigidBodyState, b: RigidBodyState, t: f64) {
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

#[test]
fn g1_loaded_hinge_chain_matches_independent_nonlinear_rotor_and_reverses() {
    let moments = [0.7, 1.2, 2.1];
    let p = core::array::from_fn::<_, 3, _>(|i| {
        MassProperties::new(
            1. + i as f64,
            Vec3::ZERO,
            Vec3::new(moments[i], moments[i], moments[i]),
        )
        .unwrap()
    });
    let initial = core::array::from_fn::<_, 3, _>(|i| {
        state(
            p[i],
            Vec3::ZERO,
            UnitQuaternion::from_axis_angle(Y, 0.3).unwrap(),
            Vec3::ZERO,
            Y.scale(0.45),
        )
    });
    let joints = [
        hinge([0, 1], [Vec3::ZERO; 2], [Z; 2], 0.8),
        hinge([1, 2], [Vec3::ZERO; 2], [Z; 2], 1.1),
    ];
    let load = PreparedSpatialLoads::new(
        Wrench::ZERO,
        Vec3::ZERO,
        &[
            WorldPointForce {
                arm_body: X.scale(0.5),
                force_world: Z.scale(-0.9),
            },
            WorldPointForce {
                arm_body: X.scale(-0.5),
                force_world: Z.scale(0.9),
            },
        ],
        || false,
    )
    .unwrap();
    let loads = [load, PreparedSpatialLoads::ZERO, PreparedSpatialLoads::ZERO];
    // Independent continuous reduction: all transverse angular accelerations
    // coincide, so sum(I)*theta''=0.9*cos(theta). The free hinge spin is zero.
    let rhs = |y: [f64; 2]| [y[1], 0.9 * y[0].cos() / 4.];
    let time = 0.8;
    let mut expected = [0.3, 0.45];
    let h = time / 8192.;
    for _ in 0..8192 {
        let a = rhs(expected);
        let b = rhs(core::array::from_fn(|i| expected[i] + h * a[i] / 2.));
        let c = rhs(core::array::from_fn(|i| expected[i] + h * b[i] / 2.));
        let d = rhs(core::array::from_fn(|i| expected[i] + h * c[i]));
        expected =
            core::array::from_fn(|i| expected[i] + h * (a[i] + 2. * b[i] + 2. * c[i] + d[i]) / 6.);
    }
    let mut errors = Vec::new();
    for count in [16, 32, 64] {
        let h = time / f64::from(count);
        let mut s = initial.to_vec();
        let mut source_torque = Vec3::ZERO;
        let mut work = 0.;
        let mut defect = 0.;
        let before = total(&s, &p, free());
        for _ in 0..count {
            let out = advance_spatial_articulated_graph(
                &s,
                &p,
                &joints,
                free(),
                &loads,
                h,
                controls(),
                || false,
            )
            .unwrap();
            work += out.constraint_work_j;
            for (i, step) in out.step.free_steps.iter().enumerate() {
                assert_eq!(
                    *step,
                    free()
                        .step_spatial(step.step.state_before, p[i], &loads[i], h)
                        .unwrap()
                );
                source_torque = source_torque.add(step.angular_impulse_world);
                work += step.force_work_j + step.couple_work_j;
                defect += step.work_energy_residual_j;
            }
            s = out.step.states_after;
        }
        let axis = s[0].pose().orientation().rotate_body_to_world(X);
        let angle = (-axis.z).atan2(axis.x);
        errors.push(
            (angle - expected[0]).abs()
                + (s[0].angular_momentum_body().y / moments[0] - expected[1]).abs(),
        );
        for i in 1..3 {
            near(
                s[i].angular_momentum_body().y / moments[i],
                s[0].angular_momentum_body().y / moments[0],
                3e-12,
            );
        }
        let after = total(&s, &p, free());
        vector(after.0, before.0, 2e-12);
        vector(after.1.sub(before.1), source_torque, 3e-11);
        near(after.2 - before.2, work + defect, 4e-12);
        s = s.into_iter().map(reverse).collect();
        for _ in 0..count {
            s = advance_spatial_articulated_graph(
                &s,
                &p,
                &joints,
                free(),
                &loads,
                h,
                controls(),
                || false,
            )
            .unwrap()
            .step
            .states_after;
        }
        for i in 0..3 {
            compare_state(reverse(s[i]), initial[i], 4e-9);
        }
    }
    for e in errors.windows(2) {
        assert!((3.7..4.3).contains(&(e[0] / e[1])), "{errors:?}");
    }
}

#[test]
fn g1_spatially_loaded_actuator_keeps_independent_recoil_rotor_and_source_work() {
    // Two internally linked groups have masses 3 and 7. The right group feels
    // an external 2 N force; the range command drives their relative motion.
    // A separate world-held material force couple drives the left hinge rotor.
    let p = [(1., 2.), (2., 3.), (3., 4.), (4., 5.)]
        .map(|(m, i)| MassProperties::new(m, Vec3::ZERO, Vec3::new(i, i, i)).unwrap());
    let fractions = [-0.7, -0.7, 0.3, 0.3];
    let boost = Vec3::new(0.2, -0.3, 0.4);
    let initial = core::array::from_fn::<_, 4, _>(|i| {
        state(
            p[i],
            X.scale(fractions[i]),
            if i < 2 {
                UnitQuaternion::from_axis_angle(Y, 0.3).unwrap()
            } else {
                UnitQuaternion::IDENTITY
            },
            boost.add(X.scale(fractions[i] * 0.25)),
            if i < 2 { Y.scale(0.45) } else { Vec3::ZERO },
        )
    });
    let original = [
        DrivenJoint::Revolute(RevoluteJoint::new([0, 1], [Vec3::ZERO; 2], [Z; 2], 0.7).unwrap()),
        DrivenJoint::Spherical(SphericalJoint {
            bodies: [2, 3],
            points_body_m: [Vec3::ZERO; 2],
        }),
        DrivenJoint::Distance {
            link: DistanceLink {
                bodies: [1, 2],
                points_body_m: [Vec3::ZERO; 2],
                length_m: 1.,
            },
            drive: DistanceDrive {
                range_rate_mps: 0.25,
                range_acceleration_mps2: 0.4,
            },
        },
    ];
    let loads = [
        PreparedSpatialLoads::new(
            Wrench::ZERO,
            Vec3::ZERO,
            &[
                WorldPointForce {
                    arm_body: X.scale(0.5),
                    force_world: Z.scale(-0.9),
                },
                WorldPointForce {
                    arm_body: X.scale(-0.5),
                    force_world: Z.scale(0.9),
                },
            ],
            || false,
        )
        .unwrap(),
        PreparedSpatialLoads::ZERO,
        PreparedSpatialLoads::ZERO,
        PreparedSpatialLoads::new(
            Wrench {
                force_world: X.scale(2.),
                torque_body: Vec3::ZERO,
            },
            Vec3::ZERO,
            &[],
            || false,
        )
        .unwrap(),
    ];
    // Independent scalar angular equation, not the native rigid integrator:
    // (I0+I1) theta'' = 0.9 cos(theta). Translation is an exact polynomial.
    let rhs = |y: [f64; 2]| [y[1], 0.9 * y[0].cos() / 5.];
    let time = 0.6;
    let mut reference = [0.3, 0.45];
    let h = time / 8192.;
    for _ in 0..8192 {
        let a = rhs(reference);
        let b = rhs(core::array::from_fn(|i| reference[i] + h * a[i] / 2.));
        let c = rhs(core::array::from_fn(|i| reference[i] + h * b[i] / 2.));
        let d = rhs(core::array::from_fn(|i| reference[i] + h * c[i]));
        reference =
            core::array::from_fn(|i| reference[i] + h * (a[i] + 2. * b[i] + 2. * c[i] + d[i]) / 6.);
    }
    let before = total(&initial, &p, free());
    let mut errors = Vec::new();
    for count in [12, 24, 48] {
        let h = time / f64::from(count);
        let mut s = initial.to_vec();
        let mut joints = original.to_vec();
        let (mut actuation, mut constraint_defect, mut source_work, mut source_defect) =
            (0., 0., 0., 0.);
        let (mut linear, mut angular) = (Vec3::ZERO, Vec3::ZERO);
        for _ in 0..count {
            let out = advance_spatial_driven_joint_graph(
                &s,
                &p,
                &joints,
                free(),
                &loads,
                h,
                controls(),
                || false,
            )
            .unwrap();
            near(
                out.constraint_work_j,
                out.actuator_work_j + out.constraint_work_defect_j,
                2e-13,
            );
            assert_eq!(out.joints[0].actuator_work_j, 0.);
            assert_eq!(out.joints[1].actuator_work_j, 0.);
            actuation += out.actuator_work_j;
            constraint_defect += out.constraint_work_defect_j;
            angular = angular.add(out.constraint_angular_impulse_world_nms);
            for (i, step) in out.step.free_steps.iter().enumerate() {
                assert_eq!(
                    *step,
                    free()
                        .step_spatial(step.step.state_before, p[i], &loads[i], h)
                        .unwrap()
                );
                linear = linear.add(step.linear_impulse_world);
                angular = angular.add(step.angular_impulse_world);
                source_work += step.force_work_j + step.couple_work_j;
                source_defect += step.work_energy_residual_j;
            }
            joints = out.joints.iter().map(|j| j.joint_after).collect();
            s = out.step.states_after;
        }
        let length = 1. + 0.25 * time + 0.2 * time * time;
        let rate = 0.25 + 0.4 * time;
        for i in 0..4 {
            vector(
                s[i].pose().position_world(),
                boost
                    .scale(time)
                    .add(X.scale(0.1 * time * time + fractions[i] * length)),
                3e-11,
            );
            vector(
                s[i].linear_momentum_world().scale(1. / p[i].mass()),
                boost.add(X.scale(0.2 * time + fractions[i] * rate)),
                3e-11,
            );
        }
        // Reduced mass 21/10, external radial generalized force 0.6 N:
        // actuator force = (21/10)*0.4 - 0.6 = 0.24 N.
        near(actuation, 0.24 * (length - 1.), 4e-12);
        let axis = s[0].pose().orientation().rotate_body_to_world(X);
        errors.push(
            ((-axis.z).atan2(axis.x) - reference[0]).abs()
                + (s[0].angular_momentum_body().y / 2. - reference[1]).abs(),
        );
        near(
            s[1].angular_momentum_body().y / 3.,
            s[0].angular_momentum_body().y / 2.,
            3e-12,
        );
        let after = total(&s, &p, free());
        vector(after.0.sub(before.0), linear, 3e-12);
        vector(after.1.sub(before.1), angular, 3e-11);
        near(
            after.2 - before.2,
            actuation + constraint_defect + source_work + source_defect,
            8e-12,
        );
    }
    for e in errors.windows(2) {
        assert!((3.7..4.3).contains(&(e[0] / e[1])), "{errors:?}");
    }
}

#[test]
fn g4_spatial_actuation_keeps_native_parity_complete_payload_and_atomic_continuation() {
    let (p, s, passive) = mixed();
    let loads = loads(UnitQuaternion::IDENTITY);
    let joints = passive.map(|j| match j {
        ArticulatedJoint::Distance(link) => DrivenJoint::Distance {
            link,
            drive: DistanceDrive::default(),
        },
        ArticulatedJoint::Spherical(joint) => DrivenJoint::Spherical(joint),
        ArticulatedJoint::Revolute(joint) => DrivenJoint::Revolute(joint),
    });
    let ordinary = advance_spatial_articulated_graph(
        &s,
        &p,
        &passive,
        free(),
        &loads,
        0.003,
        controls(),
        || false,
    )
    .unwrap();
    let zero = advance_spatial_driven_joint_graph(
        &s,
        &p,
        &joints,
        free(),
        &loads,
        0.003,
        controls(),
        || false,
    )
    .unwrap();
    // Default radial commands deliberately retain relative COM/spin rates in
    // double-double arithmetic; the legacy static rod rounds its rate earlier.
    // Their physical motion agrees at roundoff, not necessarily bit-for-bit.
    for i in 0..s.len() {
        compare_state(
            zero.step.states_after[i],
            ordinary.step.states_after[i],
            2e-14,
        );
    }
    near(zero.constraint_work_j, ordinary.constraint_work_j, 2e-14);
    assert_eq!(zero.actuator_work_j, 0.);
    // With only passive Cartesian/angular rows there is no radial rate
    // substitution, so both native solves and complete receipts are identical.
    let passive_only = advance_spatial_articulated_graph(
        &s,
        &p,
        &passive[..2],
        free(),
        &loads,
        0.003,
        controls(),
        || false,
    )
    .unwrap();
    let driven_passive = advance_spatial_driven_joint_graph(
        &s,
        &p,
        &joints[..2],
        free(),
        &loads,
        0.003,
        controls(),
        || false,
    )
    .unwrap();
    assert!(
        driven_passive.step == passive_only.step,
        "passive native graph parity"
    );
    let mut joints = joints;
    if let DrivenJoint::Distance { drive, .. } = &mut joints[2] {
        drive.range_acceleration_mps2 = 0.12;
    }
    let mut polls = 0;
    let out = advance_spatial_driven_joint_graph(
        &s,
        &p,
        &joints,
        free(),
        &loads,
        0.003,
        controls(),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    assert_ne!(out.actuator_work_j, 0.);
    assert_eq!(
        out.planned_live_bytes,
        out.step.planned_live_bytes
            + p.len() * core::mem::size_of::<f64>()
            + joints.len() * core::mem::size_of::<DrivenJointStep>()
    );
    let mut exact = controls();
    exact.maximum_live_bytes = out.planned_live_bytes;
    assert_eq!(
        out,
        advance_spatial_driven_joint_graph(&s, &p, &joints, free(), &loads, 0.003, exact, || false)
            .unwrap()
    );
    exact.maximum_live_bytes -= 1;
    assert!(matches!(
        advance_spatial_driven_joint_graph(&s, &p, &joints, free(), &loads, 0.003, exact, || false),
        Err(DistanceGraphError::Limit { .. })
    ));
    for stop in 1..=polls {
        let mut n = 0;
        let refused = advance_spatial_driven_joint_graph(
            &s,
            &p,
            &joints,
            free(),
            &loads,
            0.003,
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
            &s,
            &p,
            &joints,
            free(),
            &loads,
            0.003,
            controls(),
            || false
        )
        .unwrap()
    );
    assert!(
        advance_spatial_driven_joint_graph(
            &s,
            &p,
            &joints,
            free(),
            &loads[..3],
            0.003,
            controls(),
            || false
        )
        .is_err()
    );
    let mut bad = joints;
    if let DrivenJoint::Distance { drive, .. } = &mut bad[2] {
        drive.range_acceleration_mps2 = f64::NAN;
    }
    assert!(
        advance_spatial_driven_joint_graph(&s, &p, &bad, free(), &loads, 0.003, controls(), || {
            false
        })
        .is_err()
    );
    let empty =
        advance_spatial_driven_joint_graph(&[], &[], &[], free(), &[], 0.003, controls(), || false)
            .unwrap();
    assert!(empty.step.states_after.is_empty());
    assert_eq!(empty.planned_live_bytes, 0);
}

fn mixed() -> (
    [MassProperties; 4],
    [RigidBodyState; 4],
    [ArticulatedJoint; 3],
) {
    let p = core::array::from_fn(|i| {
        MassProperties::from_principal_frame(
            1. + i as f64,
            Vec3::ZERO,
            Vec3::new(0.7, 1., 1.2).scale(1. + 0.2 * i as f64),
            UnitQuaternion::from_axis_angle(Vec3::new(1., 2., -1.), 0.3 + 0.2 * i as f64).unwrap(),
        )
        .unwrap()
    });
    let q = core::array::from_fn::<_, 4, _>(|i| {
        UnitQuaternion::from_axis_angle(Vec3::new(1., -0.5, 2.), 0.2 - 0.3 * i as f64).unwrap()
    });
    let x = [
        Vec3::new(-0.4, 0.2, -0.1),
        Vec3::new(0.3, 0.1, 0.2),
        Vec3::new(0.2, -0.5, 0.1),
        Vec3::new(0.9, -0.3, 0.4),
    ];
    let w = Vec3::new(0.2, -0.1, 0.3);
    let s = core::array::from_fn(|i| {
        state(
            p[i],
            x[i],
            q[i],
            Vec3::new(0.1, -0.2, 0.3).add(w.cross(x[i])),
            q[i].rotate_world_to_body(w),
        )
    });
    let hp = Vec3::new(0.1, 0., 0.05);
    let bp = Vec3::new(0.3, -0.2, 0.25);
    let axis = Vec3::new(0.3, 0.5, 0.8);
    let j = [
        hinge(
            [0, 1],
            core::array::from_fn(|i| q[i].rotate_world_to_body(hp.sub(x[i]))),
            [
                q[0].rotate_world_to_body(axis),
                q[1].rotate_world_to_body(axis),
            ],
            0.8,
        ),
        ArticulatedJoint::Spherical(SphericalJoint {
            bodies: [1, 2],
            points_body_m: [
                q[1].rotate_world_to_body(bp.sub(x[1])),
                q[2].rotate_world_to_body(bp.sub(x[2])),
            ],
        }),
        ArticulatedJoint::Distance(DistanceLink {
            bodies: [2, 3],
            points_body_m: [Vec3::ZERO; 2],
            length_m: magnitude(x[3].sub(x[2])),
        }),
    ];
    (p, s, j)
}
fn loads(q: UnitQuaternion) -> [PreparedSpatialLoads; 4] {
    core::array::from_fn(|i| {
        let f = 0.2 + 0.1 * i as f64;
        PreparedSpatialLoads::new(
            Wrench {
                force_world: q.rotate_body_to_world(Vec3::new(0.3, -0.2, 0.1).scale(f)),
                torque_body: Vec3::new(0.02, -0.03, 0.01).scale(f),
            },
            q.rotate_body_to_world(Vec3::new(-0.01, 0.02, 0.03).scale(f)),
            &[WorldPointForce {
                arm_body: Vec3::new(0.4, -0.2, 0.3),
                force_world: q.rotate_body_to_world(Vec3::new(-0.1, 0.3, -0.2).scale(f)),
            }],
            || false,
        )
        .unwrap()
    })
}

#[test]
fn g0_loaded_full_tensor_hinge_ball_rod_closes_source_and_constraint_balances() {
    let (p, initial, joints) = mixed();
    let loads = loads(UnitQuaternion::IDENTITY);
    let engine = RigidBodyIntegrator::new(Gravity::new(Vec3::new(0.1, -0.2, -0.3)).unwrap());
    let before = total(&initial, &p, engine);
    let mut s = initial.to_vec();
    let (mut linear, mut angular, mut work, mut defect) = (Vec3::ZERO, Vec3::ZERO, 0., 0.);
    let mut largest_couple: f64 = 0.;
    for _ in 0..160 {
        let out = advance_spatial_articulated_graph(
            &s,
            &p,
            &joints,
            engine,
            &loads,
            0.003,
            controls(),
            || false,
        )
        .unwrap();
        work += out.constraint_work_j;
        angular = angular.add(out.constraint_angular_impulse_world_nms);
        for event in out
            .step
            .position_impulses
            .iter()
            .chain(&out.step.velocity_impulses)
        {
            if let JointReactionImpulse::Couple(e) = event {
                largest_couple = largest_couple.max(magnitude(e.couple_on_second_world_nms));
            }
        }
        for (i, step) in out.step.free_steps.iter().enumerate() {
            assert_eq!(
                *step,
                engine
                    .step_spatial(step.step.state_before, p[i], &loads[i], 0.003)
                    .unwrap()
            );
            assert_eq!(
                out.step.states_after[i].pose(),
                step.step.state_after.pose()
            );
            linear = linear.add(step.linear_impulse_world);
            angular = angular.add(step.angular_impulse_world);
            work += step.force_work_j + step.couple_work_j;
            defect += step.work_energy_residual_j;
        }
        s = out.step.states_after;
        let after = total(&s, &p, engine);
        vector(after.0.sub(before.0), linear, 4e-12);
        vector(after.1.sub(before.1), angular, 5e-11);
        near(after.2 - before.2, work + defect, 7e-12);
    }
    assert!(largest_couple > 1e-6);
}

#[test]
fn g3_loaded_mechanism_preserves_world_frame_and_common_boost_covariance() {
    let (p, s, joints) = mixed();
    let h = 0.004;
    let base = advance_spatial_articulated_graph(
        &s,
        &p,
        &joints,
        free(),
        &loads(UnitQuaternion::IDENTITY),
        h,
        controls(),
        || false,
    )
    .unwrap();
    let q = UnitQuaternion::from_axis_angle(Vec3::new(1., -2., 3.), 0.7).unwrap();
    let shift = Vec3::new(2., -3., 1.);
    let boost = Vec3::new(9., -5., 3.);
    let moved = core::array::from_fn::<_, 4, _>(|i| {
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
    let out = advance_spatial_articulated_graph(
        &moved,
        &p,
        &joints,
        free(),
        &loads(q),
        h,
        controls(),
        || false,
    )
    .unwrap();
    for i in 0..4 {
        let expected = RigidBodyState::new(
            Pose::new(
                q.rotate_body_to_world(base.step.states_after[i].pose().position_world())
                    .add(shift)
                    .add(boost.scale(h)),
                compose(q, base.step.states_after[i].pose().orientation()),
            )
            .unwrap(),
            q.rotate_body_to_world(base.step.states_after[i].linear_momentum_world())
                .add(boost.scale(p[i].mass())),
            base.step.states_after[i].angular_momentum_body(),
        )
        .unwrap();
        compare_state(out.step.states_after[i], expected, 3e-11);
        near(
            out.step.free_steps[i].force_work_j,
            base.step.free_steps[i].force_work_j
                + boost.dot(q.rotate_body_to_world(base.step.free_steps[i].linear_impulse_world)),
            3e-12,
        );
    }
    near(out.constraint_work_j, base.constraint_work_j, 5e-13);
}

#[test]
fn g0_unlinked_spatial_bodies_keep_exact_native_results_and_complete_payload() {
    let (p, s, joints) = mixed();
    let loads = loads(UnitQuaternion::IDENTITY);
    let h = 0.003;
    let e = RigidBodyIntegrator::new(Gravity::new(Vec3::new(0., -9.81, 0.)).unwrap());
    let free_result =
        advance_spatial_articulated_graph(&s, &p, &[], e, &loads, h, controls(), || false).unwrap();
    for i in 0..4 {
        let native = e.step_spatial(s[i], p[i], &loads[i], h).unwrap();
        assert_eq!(free_result.step.free_steps[i], native);
        assert_eq!(free_result.step.states_after[i], native.step.state_after);
    }
    assert!(free_result.joint_rows.is_empty());
    near(free_result.constraint_work_j, 0., 0.);
    let base =
        advance_spatial_articulated_graph(&s, &p, &joints, e, &loads, h, controls(), || false)
            .unwrap();
    let extra = 4 * core::mem::size_of::<f64>()
        + 3 * core::mem::size_of::<fs_mbd::rigid_rattle::distance::graph::JointRows>();
    assert_eq!(
        base.planned_live_bytes,
        base.step.planned_live_bytes + extra
    );
    let mut c = controls();
    c.maximum_live_bytes = base.planned_live_bytes;
    assert_eq!(
        advance_spatial_articulated_graph(&s, &p, &joints, e, &loads, h, c, || false).unwrap(),
        base
    );
    c.maximum_live_bytes -= 1;
    assert!(matches!(
        advance_spatial_articulated_graph(&s, &p, &joints, e, &loads, h, c, || false),
        Err(DistanceGraphError::Limit {
            resource: "live bytes",
            ..
        })
    ));
    let empty =
        advance_spatial_articulated_graph(&[], &[], &[], e, &[], h, controls(), || false).unwrap();
    assert!(empty.step.states_after.is_empty());
    assert_eq!(empty.planned_live_bytes, 0);
}

#[test]
fn g4_all_spatial_graph_cancellation_points_and_load_refusals_are_atomic() {
    let (p, s, joints) = mixed();
    let loads = loads(UnitQuaternion::IDENTITY);
    let mut count = 0;
    let baseline = advance_spatial_articulated_graph(
        &s,
        &p,
        &joints,
        free(),
        &loads,
        0.004,
        controls(),
        || {
            count += 1;
            false
        },
    )
    .unwrap();
    assert!(count > 50);
    for stop in 0..count {
        let mut calls = 0;
        assert_eq!(
            advance_spatial_articulated_graph(
                &s,
                &p,
                &joints,
                free(),
                &loads,
                0.004,
                controls(),
                || {
                    let halt = calls == stop;
                    calls += 1;
                    halt
                }
            ),
            Err(DistanceGraphError::Step(DistanceStepError::Cancelled))
        );
    }
    assert_eq!(
        advance_spatial_articulated_graph(
            &s,
            &p,
            &joints,
            free(),
            &loads,
            0.004,
            controls(),
            || false
        )
        .unwrap(),
        baseline
    );
    assert!(
        advance_spatial_articulated_graph(
            &s,
            &p,
            &joints,
            free(),
            &loads[..3],
            0.004,
            controls(),
            || false
        )
        .is_err()
    );
    let mut c = controls();
    c.step.maximum_rotation_increment_rad = 1e-8;
    assert_eq!(
        advance_spatial_articulated_graph(&s, &p, &joints, free(), &loads, 0.004, c, || false),
        Err(DistanceGraphError::Step(DistanceStepError::LocalBranch))
    );
    for h in [0., -0.001, f64::INFINITY, f64::NAN] {
        assert!(
            advance_spatial_articulated_graph(
                &s,
                &p,
                &joints,
                free(),
                &loads,
                h,
                controls(),
                || false
            )
            .is_err()
        );
    }
    let mut bad = loads;
    bad[3] = PreparedSpatialLoads::new(
        Wrench {
            force_world: Vec3::ZERO,
            torque_body: Vec3::new(1e308, 1e308, 0.),
        },
        Vec3::new(1e308, 1e308, 0.),
        &[],
        || false,
    )
    .unwrap();
    assert!(
        advance_spatial_articulated_graph(&s, &p, &joints, free(), &bad, 0.004, controls(), || {
            false
        })
        .is_err()
    );
}
