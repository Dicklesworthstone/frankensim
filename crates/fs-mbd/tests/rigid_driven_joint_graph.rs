use fs_mbd::rigid_rattle::distance::graph::{
    ArticulatedJoint, DistanceDrive, DistanceGraphControls, DistanceGraphError, DistanceLink,
    DrivenJoint, DrivenJointGraphStep, JointReactionImpulse, RevoluteJoint, SphericalJoint,
    advance_articulated_graph, advance_driven_distance_graph, advance_driven_joint_graph,
};
use fs_mbd::rigid_rattle::distance::{DistanceStepControls, DistanceStepError};
use fs_mbd::{
    Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState, UnitQuaternion, Vec3,
    Wrench,
};

const X: Vec3 = Vec3::new(1., 0., 0.);
const Y: Vec3 = Vec3::new(0., 1., 0.);
const Z: Vec3 = Vec3::new(0., 0., 1.);

fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:e} != {b:e}, tolerance {tol:e}");
}
fn vector(a: Vec3, b: Vec3, tol: f64) {
    near(a.x, b.x, tol);
    near(a.y, b.y, tol);
    near(a.z, b.z, tol);
}
fn norm(v: Vec3) -> f64 {
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
            maximum_position_iterations: 40,
            maximum_backtracks: 16,
        },
        position_tolerance_m: 3e-13,
        velocity_tolerance_mps: 3e-13,
        minimum_scaled_pivot: 1e-10,
        maximum_bodies: 64,
        maximum_constraints: 32,
        maximum_trial_evaluations: 160,
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
fn center_link(a: usize, b: usize, length: f64) -> DistanceLink {
    DistanceLink {
        bodies: [a, b],
        points_body_m: [Vec3::ZERO; 2],
        length_m: length,
    }
}
fn run(
    s: &[RigidBodyState],
    p: &[MassProperties],
    joints: &[DrivenJoint],
    h: f64,
) -> DrivenJointGraphStep {
    advance_driven_joint_graph(
        s,
        p,
        joints,
        free(),
        &vec![Wrench::ZERO; s.len()],
        h,
        controls(),
        || false,
    )
    .unwrap()
}
fn totals(s: &[RigidBodyState], p: &[MassProperties]) -> (Vec3, Vec3, f64) {
    s.iter().zip(p).fold(
        (Vec3::ZERO, Vec3::ZERO, 0.),
        |(linear, angular, energy), (&s, &p)| {
            let d = free().diagnostics(s, p).unwrap();
            (
                linear.add(d.linear_momentum_world),
                angular.add(d.angular_momentum_world),
                energy + d.mechanical_energy,
            )
        },
    )
}
fn state_error(a: RigidBodyState, b: RigidBodyState) -> f64 {
    let mut result = norm(a.pose().position_world().sub(b.pose().position_world()))
        .max(norm(
            a.linear_momentum_world().sub(b.linear_momentum_world()),
        ))
        .max(norm(
            a.angular_momentum_body().sub(b.angular_momentum_body()),
        ));
    for axis in [X, Y, Z] {
        result = result.max(norm(
            a.pose()
                .orientation()
                .rotate_body_to_world(axis)
                .sub(b.pose().orientation().rotate_body_to_world(axis)),
        ));
    }
    result
}
fn continuation(out: &DrivenJointGraphStep) -> Vec<DrivenJoint> {
    out.joints.iter().map(|joint| joint.joint_after).collect()
}
fn reverse(s: RigidBodyState) -> RigidBodyState {
    RigidBodyState::new(
        s.pose(),
        s.linear_momentum_world().scale(-1.),
        s.angular_momentum_body().scale(-1.),
    )
    .unwrap()
}
fn reversed_joint(joint: DrivenJoint) -> DrivenJoint {
    match joint {
        DrivenJoint::Distance { link, drive } => DrivenJoint::Distance {
            link,
            drive: DistanceDrive {
                range_rate_mps: -drive.range_rate_mps,
                ..drive
            },
        },
        passive => passive,
    }
}

#[test]
fn g1_driven_rod_recoils_two_finite_joint_groups_and_transmits_hinge_couples() {
    // Hinge group mass 3, spherical group mass 7, so their independent reduced
    // mass is 21/10. A transverse torque also exercises the hinge couple rows.
    let p = [(1., 2.), (2., 3.), (3., 4.), (4., 5.)]
        .map(|(m, i)| MassProperties::new(m, Vec3::ZERO, Vec3::new(i, i, i)).unwrap());
    let h = 0.2;
    let rate = 0.25;
    let boost = Vec3::new(0.2, -0.3, 0.4);
    let fractions = [-0.7, -0.7, 0.3, 0.3];
    let initial = std::array::from_fn::<_, 4, _>(|i| {
        state(
            p[i],
            X.scale(fractions[i]),
            UnitQuaternion::IDENTITY,
            boost.add(X.scale(fractions[i] * rate)),
            Vec3::ZERO,
        )
    });
    let mut loads = [Wrench::ZERO; 4];
    loads[0].torque_body = Y.scale(5.);
    for acceleration in [0.4, -0.4] {
        let joints = [
            DrivenJoint::Revolute(
                RevoluteJoint::new([0, 1], [Vec3::ZERO; 2], [Z; 2], 0.7).unwrap(),
            ),
            DrivenJoint::Spherical(SphericalJoint {
                bodies: [2, 3],
                points_body_m: [Vec3::ZERO; 2],
            }),
            DrivenJoint::Distance {
                link: center_link(1, 2, 1.),
                drive: DistanceDrive {
                    range_rate_mps: rate,
                    range_acceleration_mps2: acceleration,
                },
            },
        ];
        let out = advance_driven_joint_graph(
            &initial,
            &p,
            &joints,
            free(),
            &loads,
            h,
            controls(),
            || false,
        )
        .unwrap();
        let end_rate = rate + acceleration * h;
        let end_length = 1. + rate * h + 0.5 * acceleration * h * h;
        for i in 0..4 {
            let motion = out.step.states_after[i]
                .point_kinematics(p[i], Vec3::ZERO)
                .unwrap();
            vector(
                out.step.states_after[i].pose().position_world(),
                boost.scale(h).add(X.scale(fractions[i] * end_length)),
                4e-13,
            );
            vector(
                motion.center_of_mass_velocity_world,
                boost.add(X.scale(fractions[i] * end_rate)),
                4e-13,
            );
            vector(
                motion.angular_velocity_world,
                if i < 2 { Y.scale(h) } else { Vec3::ZERO },
                4e-13,
            );
        }
        assert_eq!(
            out.joints
                .iter()
                .map(|j| (j.rows.first_row, j.rows.row_count))
                .collect::<Vec<_>>(),
            [(0, 5), (5, 3), (8, 1)]
        );
        for events in [&out.step.position_impulses, &out.step.velocity_impulses] {
            let couple = events[3..5].iter().fold(Vec3::ZERO, |sum, event| {
                let JointReactionImpulse::Couple(event) = event else {
                    panic!("native couple required")
                };
                vector(event.first.linear_increment_world, Vec3::ZERO, 0.);
                vector(event.second.linear_increment_world, Vec3::ZERO, 0.);
                sum.add(event.couple_on_second_world_nms)
            });
            vector(couple, Y.scale(1.5 * h), 5e-12);
            let JointReactionImpulse::Point(event) = events[8] else {
                panic!("native point required")
            };
            vector(
                event.second.impulse_world,
                X.scale(0.5 * 2.1 * acceleration * h),
                5e-12,
            );
        }
        let expected_work = 0.5 * 2.1 * (end_rate * end_rate - rate * rate);
        near(out.actuator_work_j, expected_work, 5e-12);
        assert_eq!(out.joints[0].actuator_work_j, 0.);
        assert_eq!(out.joints[1].actuator_work_j, 0.);
        near(out.joints[2].actuator_work_j, expected_work, 5e-12);
        near(out.constraint_work_defect_j, 0., 5e-12);
        vector(
            totals(&out.step.states_after, &p).0,
            totals(&initial, &p).0,
            3e-13,
        );
    }
}

// Independent planar free-floating mechanism: two hinged bodies, a finite
// spherical payload at body two's COM, and a rod setting the hinge opening.
// The actuator endpoints are measured from the common hinge, not either COM.
// Conservation gives H = A(phi)*theta_dot + B(phi)*phi_dot. The only reference
// integration below is scalar Simpson quadrature of that analytic connection.
fn mechanism_properties() -> [MassProperties; 3] {
    [
        MassProperties::new(2., Vec3::ZERO, Vec3::new(0.4, 0.5, 0.6)).unwrap(),
        MassProperties::new(3., Vec3::ZERO, Vec3::new(0.7, 0.8, 0.9)).unwrap(),
        MassProperties::new(1., Vec3::ZERO, Vec3::new(0.2, 0.3, 0.4)).unwrap(),
    ]
}
fn shape(t: f64) -> (f64, f64, f64, f64) {
    let c = 0.7 * 0.9;
    let l0 = (0.7_f64.powi(2) + 0.9_f64.powi(2) - 2. * c * 0.8_f64.cos()).sqrt();
    let l = l0 + 0.12 * t - 0.025 * t * t;
    let rate = 0.12 - 0.05 * t;
    let phi = ((0.7_f64.powi(2) + 0.9_f64.powi(2) - l * l) / (2. * c)).acos();
    (l, rate, phi, l * rate / (c * phi.sin()))
}
fn connection(phi: f64) -> (f64, f64) {
    let mu = 4. / 3.;
    let a = -0.4;
    let b = 0.6;
    (
        1.5 + mu * (a * a + b * b - 2. * a * b * phi.cos()),
        0.9 + mu * (b * b - a * b * phi.cos()),
    )
}
fn mechanism_reference(t: f64) -> [RigidBodyState; 3] {
    let (_, _, phi0, phidot0) = shape(0.);
    let (a0, b0) = connection(phi0);
    let angular_momentum = a0 * 0.25 + b0 * phidot0;
    let theta_dot = |time| {
        let (_, _, phi, phidot) = shape(time);
        let (a, b) = connection(phi);
        (angular_momentum - b * phidot) / a
    };
    let pieces = 4096;
    let dt = t / pieces as f64;
    let mut quadrature = theta_dot(0.) + theta_dot(t);
    for i in 1..pieces {
        quadrature += if i % 2 == 0 { 2. } else { 4. } * theta_dot(i as f64 * dt);
    }
    let theta = 0.15 + dt * quadrature / 3.;
    let (_, _, phi, phidot) = shape(t);
    let omega = theta_dot(t);
    let q =
        |angle: f64| UnitQuaternion::new((0.5 * angle).cos(), 0., 0., (0.5 * angle).sin()).unwrap();
    let q0 = q(theta);
    let q1 = q(theta + phi);
    let r0 = q0.rotate_body_to_world(X.scale(-0.4));
    let r1 = q1.rotate_body_to_world(X.scale(0.6));
    let separation = r1.sub(r0);
    let relative_velocity = Z
        .cross(r1)
        .scale(omega + phidot)
        .sub(Z.cross(r0).scale(omega));
    let boost = Vec3::new(0.1, -0.2, 0.3);
    let p = mechanism_properties();
    [
        state(
            p[0],
            boost.scale(t).add(separation.scale(-2. / 3.)),
            q0,
            boost.add(relative_velocity.scale(-2. / 3.)),
            Z.scale(omega),
        ),
        state(
            p[1],
            boost.scale(t).add(separation.scale(1. / 3.)),
            q1,
            boost.add(relative_velocity.scale(1. / 3.)),
            Z.scale(omega + phidot),
        ),
        state(
            p[2],
            boost.scale(t).add(separation.scale(1. / 3.)),
            UnitQuaternion::IDENTITY,
            boost.add(relative_velocity.scale(1. / 3.)),
            Vec3::ZERO,
        ),
    ]
}
fn mechanism_joints() -> [DrivenJoint; 3] {
    let (length, rate, _, _) = shape(0.);
    [
        DrivenJoint::Revolute(
            RevoluteJoint::new([0, 1], [X.scale(0.4), X.scale(-0.6)], [Z; 2], 0.8).unwrap(),
        ),
        DrivenJoint::Spherical(SphericalJoint {
            bodies: [1, 2],
            points_body_m: [Vec3::ZERO; 2],
        }),
        DrivenJoint::Distance {
            link: DistanceLink {
                bodies: [0, 1],
                points_body_m: [X.scale(1.1), X.scale(0.3)],
                length_m: length,
            },
            drive: DistanceDrive {
                range_rate_mps: rate,
                range_acceleration_mps2: -0.05,
            },
        },
    ]
}

#[test]
fn g1_actuated_hinge_mechanism_matches_continuous_finite_recoil_and_source_work() {
    let p = mechanism_properties();
    let initial = mechanism_reference(0.);
    let end = 0.3;
    let reference = mechanism_reference(end);
    let expected_work = totals(&reference, &p).2 - totals(&initial, &p).2;
    let mut errors = Vec::new();
    let mut work_errors = Vec::new();
    for count in [15, 30, 60] {
        let mut states = initial.to_vec();
        let mut joints = mechanism_joints().to_vec();
        let mut source_work = 0.;
        for _ in 0..count {
            let out = run(&states, &p, &joints, end / count as f64);
            source_work += out.actuator_work_j;
            assert_eq!(out.joints[0].actuator_work_j, 0.);
            assert_eq!(out.joints[1].actuator_work_j, 0.);
            joints = continuation(&out);
            states = out.step.states_after;
        }
        vector(totals(&states, &p).0, totals(&initial, &p).0, 2e-12);
        vector(totals(&states, &p).1, totals(&initial, &p).1, 2e-11);
        errors.push(
            states
                .iter()
                .zip(reference)
                .fold(0_f64, |error, (&s, r)| error.max(state_error(s, r))),
        );
        work_errors.push((source_work - expected_work).abs());
    }
    assert!(errors[2] < 2e-6, "state errors {errors:?}");
    assert!(
        errors[1] < 0.32 * errors[0] && errors[2] < 0.32 * errors[1],
        "state errors {errors:?}"
    );
    assert!(work_errors[2] < 2e-7, "source work errors {work_errors:?}");
    assert!(
        work_errors[1] < 0.35 * work_errors[0] && work_errors[2] < 0.35 * work_errors[1],
        "source work errors {work_errors:?}"
    );
}

fn full_tensor() -> ([MassProperties; 4], [RigidBodyState; 4], [DrivenJoint; 3]) {
    let p = std::array::from_fn(|i| {
        MassProperties::from_principal_frame(
            1. + i as f64,
            Vec3::ZERO,
            Vec3::new(0.7, 1., 1.2).scale(1. + 0.2 * i as f64),
            UnitQuaternion::from_axis_angle(Vec3::new(1., 2., -1.), 0.3 + 0.2 * i as f64).unwrap(),
        )
        .unwrap()
    });
    let q = std::array::from_fn::<_, 4, _>(|i| {
        UnitQuaternion::from_axis_angle(Vec3::new(1., -0.5, 2.), 0.2 - 0.3 * i as f64).unwrap()
    });
    let x = [
        Vec3::new(-0.4, 0.2, -0.1),
        Vec3::new(0.3, 0.1, 0.2),
        Vec3::new(0.2, -0.5, 0.1),
        Vec3::new(0.9, -0.3, 0.4),
    ];
    let omega = Vec3::new(0.2, -0.1, 0.3);
    let s = std::array::from_fn(|i| {
        state(
            p[i],
            x[i],
            q[i],
            Vec3::new(0.1, -0.2, 0.3).add(omega.cross(x[i])),
            q[i].rotate_world_to_body(omega),
        )
    });
    let hpoint = Vec3::new(0.1, 0., 0.05);
    let bpoint = Vec3::new(0.3, -0.2, 0.25);
    let axis = Vec3::new(0.3, 0.5, 0.8);
    let joints = [
        DrivenJoint::Revolute(
            RevoluteJoint::new(
                [0, 1],
                [
                    q[0].rotate_world_to_body(hpoint.sub(x[0])),
                    q[1].rotate_world_to_body(hpoint.sub(x[1])),
                ],
                [
                    q[0].rotate_world_to_body(axis),
                    q[1].rotate_world_to_body(axis),
                ],
                0.8,
            )
            .unwrap(),
        ),
        DrivenJoint::Spherical(SphericalJoint {
            bodies: [1, 2],
            points_body_m: [
                q[1].rotate_world_to_body(bpoint.sub(x[1])),
                q[2].rotate_world_to_body(bpoint.sub(x[2])),
            ],
        }),
        DrivenJoint::Distance {
            link: center_link(2, 3, norm(x[3].sub(x[2]))),
            drive: DistanceDrive {
                range_rate_mps: 0.,
                range_acceleration_mps2: 0.4,
            },
        },
    ];
    (p, s, joints)
}

#[test]
fn g0_off_center_full_tensor_actuation_retains_momenta_native_work_and_reversal() {
    let (p, initial, initial_joints) = full_tensor();
    let before = totals(&initial, &p);
    let mut states = initial.to_vec();
    let mut joints = initial_joints.to_vec();
    let mut work = 0.;
    let mut requested_angular = Vec3::ZERO;
    let mut largest_couple = 0_f64;
    let mut source = 0.;
    for _ in 0..100 {
        let out = run(&states, &p, &joints, 0.005);
        let mut event_defect = 0.;
        for event in out
            .step
            .position_impulses
            .iter()
            .chain(&out.step.velocity_impulses)
        {
            match event {
                JointReactionImpulse::Point(event) => {
                    vector(event.impulse_balance_world, Vec3::ZERO, 0.);
                    vector(
                        event.angular_momentum_change_world,
                        event.angular_impulse_world,
                        2e-15,
                    );
                    event_defect += event.work_energy_residual_j;
                }
                JointReactionImpulse::Couple(event) => {
                    vector(event.first.linear_increment_world, Vec3::ZERO, 0.);
                    vector(event.second.linear_increment_world, Vec3::ZERO, 0.);
                    vector(event.angular_momentum_change_world_nms, Vec3::ZERO, 2e-15);
                    event_defect += event.work_energy_residual_j;
                    largest_couple = largest_couple.max(norm(event.couple_on_second_world_nms));
                }
            }
        }
        let free_work: f64 = out
            .step
            .free_steps
            .iter()
            .map(|step| {
                step.diagnostics_after.mechanical_energy - step.diagnostics_before.mechanical_energy
            })
            .sum();
        near(
            out.actuator_work_j + out.constraint_work_defect_j,
            out.constraint_work_j,
            3e-15,
        );
        near(
            out.joints.iter().map(|j| j.constraint_work_defect_j).sum(),
            out.constraint_work_defect_j,
            3e-15,
        );
        assert_eq!(out.joints[0].actuator_work_j, 0.);
        assert_eq!(out.joints[1].actuator_work_j, 0.);
        work += out.actuator_work_j + out.constraint_work_defect_j + free_work + event_defect;
        source += out.actuator_work_j;
        requested_angular = requested_angular.add(out.constraint_angular_impulse_world_nms);
        for i in 0..4 {
            assert_eq!(
                out.step.states_after[i].pose(),
                out.step.free_steps[i].state_after.pose()
            );
        }
        let r = &out.step.residuals;
        assert!(r[3].gap_m.hypot(r[4].gap_m) <= controls().position_tolerance_m);
        assert!(
            r[3].range_rate_mps.hypot(r[4].range_rate_mps) <= controls().velocity_tolerance_mps
        );
        assert!(r[8].gap_m.abs() <= controls().position_tolerance_m);
        assert!(r[8].range_rate_mps.abs() <= controls().velocity_tolerance_mps);
        joints = continuation(&out);
        states = out.step.states_after;
    }
    let after = totals(&states, &p);
    vector(after.0, before.0, 3e-12);
    vector(after.1.sub(before.1), requested_angular, 3e-11);
    near(after.2 - before.2, work, 4e-12);
    assert!(
        source.abs() > 1e-3,
        "actuator must do appreciable work: {source}"
    );
    assert!(
        largest_couple > 1e-5,
        "hinge must transmit a nonzero couple"
    );

    let forward = run(&initial, &p, &initial_joints, 0.02);
    let backward = run(
        &forward
            .step
            .states_after
            .iter()
            .copied()
            .map(reverse)
            .collect::<Vec<_>>(),
        &p,
        &continuation(&forward)
            .into_iter()
            .map(reversed_joint)
            .collect::<Vec<_>>(),
        0.02,
    );
    for (&s, original) in backward.step.states_after.iter().zip(initial) {
        assert!(state_error(reverse(s), original) < 2e-11);
    }
    near(
        forward.actuator_work_j + backward.actuator_work_j,
        0.,
        2e-12,
    );
}

fn passive(joint: DrivenJoint) -> ArticulatedJoint {
    match joint {
        DrivenJoint::Distance { link, .. } => ArticulatedJoint::Distance(link),
        DrivenJoint::Spherical(joint) => ArticulatedJoint::Spherical(joint),
        DrivenJoint::Revolute(joint) => ArticulatedJoint::Revolute(joint),
        DrivenJoint::PointPlane(joint) => ArticulatedJoint::PointPlane(joint),
        DrivenJoint::Prismatic(joint) => ArticulatedJoint::Prismatic(joint),
        DrivenJoint::RevoluteMotor { .. } => {
            unreachable!("this passive fixture has no angular motor")
        }
    }
}

#[test]
fn g0_zero_drive_passive_only_and_rod_only_paths_keep_native_parity() {
    let (p, s, mut joints) = full_tensor();
    let DrivenJoint::Distance { link, .. } = joints[2] else {
        unreachable!()
    };
    joints[2] = DrivenJoint::Distance {
        link,
        drive: DistanceDrive::default(),
    };
    for count in [0, 2, 3] {
        let driven = run(&s, &p, &joints[..count], 0.01);
        let original = advance_articulated_graph(
            &s,
            &p,
            &joints[..count]
                .iter()
                .copied()
                .map(passive)
                .collect::<Vec<_>>(),
            free(),
            &[Wrench::ZERO; 4],
            0.01,
            controls(),
            || false,
        )
        .unwrap();
        if count < 3 {
            assert_eq!(driven.step, original.step);
        } else {
            for (&a, &b) in driven
                .step
                .states_after
                .iter()
                .zip(&original.step.states_after)
            {
                assert!(state_error(a, b) < 2e-13);
            }
            near(driven.constraint_work_j, original.constraint_work_j, 1e-15);
        }
        assert_eq!(driven.actuator_work_j, 0.);
        assert_eq!(
            driven.joints.iter().map(|j| j.rows).collect::<Vec<_>>(),
            original.joint_rows
        );
        if count < 3 {
            assert_eq!(
                driven.step.free_steps[3],
                free().step(s[3], p[3], Wrench::ZERO, 0.01).unwrap()
            );
        }
    }

    let p = [MassProperties::new(1., Vec3::ZERO, Vec3::new(1., 1., 1.)).unwrap(); 2];
    let s = [
        state(
            p[0],
            Vec3::ZERO,
            UnitQuaternion::IDENTITY,
            X.scale(-0.125),
            Vec3::ZERO,
        ),
        state(
            p[1],
            X,
            UnitQuaternion::IDENTITY,
            X.scale(0.125),
            Vec3::ZERO,
        ),
    ];
    let link = center_link(0, 1, 1.);
    let drive = DistanceDrive {
        range_rate_mps: 0.25,
        range_acceleration_mps2: 0.5,
    };
    let out = run(&s, &p, &[DrivenJoint::Distance { link, drive }], 0.125);
    let original = advance_driven_distance_graph(
        &s,
        &p,
        &[link],
        &[drive],
        free(),
        &[Wrench::ZERO; 2],
        0.125,
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(out.step.states_after, original.step.states_after);
    assert_eq!(out.step.free_steps, original.step.free_steps);
    assert_eq!(out.step.residuals, original.step.residuals);
    for (mixed, point) in out
        .step
        .position_impulses
        .iter()
        .chain(&out.step.velocity_impulses)
        .zip(
            original
                .step
                .position_impulses
                .iter()
                .chain(&original.step.velocity_impulses),
        )
    {
        assert_eq!(*mixed, JointReactionImpulse::Point(*point));
    }
    assert_eq!(out.actuator_work_j, original.actuator_work_j);
    assert_eq!(
        out.constraint_work_defect_j,
        original.constraint_work_defect_j
    );
    let continued = original.links[0];
    assert_eq!(
        out.joints[0].joint_after,
        DrivenJoint::Distance {
            link: continued.link_after,
            drive: continued.drive_after
        }
    );
}

fn compose(a: UnitQuaternion, b: UnitQuaternion) -> UnitQuaternion {
    let [w, x, y, z] = a.components();
    let [v, i, j, k] = b.components();
    UnitQuaternion::new(
        w * v - x * i - y * j - z * k,
        w * i + x * v + y * k - z * j,
        w * j - x * k + y * v + z * i,
        w * k + x * j - y * i + z * v,
    )
    .unwrap()
}
fn transform(
    s: RigidBodyState,
    p: MassProperties,
    q: UnitQuaternion,
    shift: Vec3,
    boost: Vec3,
) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(
            q.rotate_body_to_world(s.pose().position_world()).add(shift),
            compose(q, s.pose().orientation()),
        )
        .unwrap(),
        q.rotate_body_to_world(s.linear_momentum_world())
            .add(boost.scale(p.mass())),
        s.angular_momentum_body(),
    )
    .unwrap()
}
fn swap_ends(joint: DrivenJoint) -> DrivenJoint {
    match joint {
        DrivenJoint::Distance { link, drive } => DrivenJoint::Distance {
            link: DistanceLink {
                bodies: [link.bodies[1], link.bodies[0]],
                points_body_m: [link.points_body_m[1], link.points_body_m[0]],
                ..link
            },
            drive,
        },
        DrivenJoint::Spherical(joint) => DrivenJoint::Spherical(SphericalJoint {
            bodies: [joint.bodies[1], joint.bodies[0]],
            points_body_m: [joint.points_body_m[1], joint.points_body_m[0]],
        }),
        DrivenJoint::Revolute(joint) => DrivenJoint::Revolute(
            RevoluteJoint::new(
                [joint.bodies()[1], joint.bodies()[0]],
                [joint.points_body_m()[1], joint.points_body_m()[0]],
                [joint.axes_body()[1], joint.axes_body()[0]],
                joint.angular_length_m(),
            )
            .unwrap(),
        ),
        _ => unreachable!("this endpoint-swap fixture uses rods, balls and hinges"),
    }
}

#[test]
fn g3_mixed_actuation_is_frame_boost_endpoint_and_order_covariant() {
    let (p, s, joints) = full_tensor();
    let h = 0.02;
    let baseline = run(&s, &p, &joints, h);
    let frame = UnitQuaternion::from_axis_angle(Vec3::new(1., -2., 3.), 0.7).unwrap();
    let shift = Vec3::new(0.3, -0.7, 0.2);
    let boost = Vec3::new(8., -4., 2.);
    let transformed: Vec<_> = s
        .iter()
        .zip(p)
        .map(|(&s, p)| transform(s, p, frame, shift, boost))
        .collect();
    let moved = run(&transformed, &p, &joints, h);
    for i in 0..4 {
        let expected = transform(
            baseline.step.states_after[i],
            p[i],
            frame,
            shift.add(boost.scale(h)),
            boost,
        );
        assert!(state_error(moved.step.states_after[i], expected) < 5e-11);
    }
    near(moved.actuator_work_j, baseline.actuator_work_j, 3e-13);
    near(
        moved.constraint_work_defect_j,
        baseline.constraint_work_defect_j,
        3e-13,
    );
    let reversed: Vec<_> = joints.into_iter().rev().map(swap_ends).collect();
    let reordered = run(&s, &p, &reversed, h);
    for i in 0..4 {
        assert!(
            state_error(
                reordered.step.states_after[i],
                baseline.step.states_after[i]
            ) < 2e-11
        );
    }
    near(reordered.actuator_work_j, baseline.actuator_work_j, 3e-13);
    near(
        reordered.constraint_work_j,
        baseline.constraint_work_j,
        3e-13,
    );

    // A common boost can hide an endpoint's spin in its rounded point velocity.
    // A passive hinge precedes the rod, so command admission must address row 5.
    let p = [MassProperties::new(1., Vec3::ZERO, Vec3::new(1., 1., 1.)).unwrap(); 3];
    let boost = X.scale(2_f64.powi(54));
    let s = [
        state(
            p[0],
            Vec3::ZERO,
            UnitQuaternion::IDENTITY,
            boost,
            Vec3::ZERO,
        ),
        state(p[1], X, UnitQuaternion::IDENTITY, boost, Z),
        state(
            p[2],
            Vec3::ZERO,
            UnitQuaternion::IDENTITY,
            boost,
            Vec3::ZERO,
        ),
    ];
    let link = DistanceLink {
        bodies: [0, 1],
        points_body_m: [Vec3::ZERO, Y],
        length_m: 2_f64.sqrt(),
    };
    assert_eq!(
        s[0].point_kinematics(p[0], Vec3::ZERO)
            .unwrap()
            .point_velocity_world,
        s[1].point_kinematics(p[1], Y).unwrap().point_velocity_world
    );
    let joints = [
        DrivenJoint::Revolute(RevoluteJoint::new([0, 2], [Vec3::ZERO; 2], [Z; 2], 1.).unwrap()),
        DrivenJoint::Distance {
            link,
            drive: DistanceDrive::default(),
        },
    ];
    assert!(matches!(
        advance_driven_joint_graph(
            &s,
            &p,
            &joints,
            free(),
            &[Wrench::ZERO; 3],
            0.001,
            controls(),
            || false
        ),
        Err(DistanceGraphError::InitialConstraint { constraint: 5 })
    ));
}

#[test]
fn g4_actual_32_row_driven_mechanism_bounds_commands_rank_trials_and_payload() {
    let p = vec![MassProperties::new(1., Vec3::ZERO, Vec3::new(1., 1., 1.)).unwrap(); 9];
    let s: Vec<_> = (0_usize..9)
        .map(|i| {
            state(
                p[i],
                X.scale(i.saturating_sub(6) as f64),
                UnitQuaternion::IDENTITY,
                Vec3::ZERO,
                Vec3::ZERO,
            )
        })
        .collect();
    let mut joints: Vec<_> = (0..6)
        .map(|i| {
            DrivenJoint::Revolute(
                RevoluteJoint::new([i, i + 1], [Vec3::ZERO; 2], [Z; 2], 1.).unwrap(),
            )
        })
        .collect();
    for (i, acceleration) in [(6, 0.01), (7, -0.01)] {
        joints.push(DrivenJoint::Distance {
            link: center_link(i, i + 1, 1.),
            drive: DistanceDrive {
                range_rate_mps: 0.,
                range_acceleration_mps2: acceleration,
            },
        });
    }
    let out = run(&s, &p, &joints, 0.1);
    assert_eq!(out.step.residuals.len(), 32);
    assert_eq!(out.joints.len(), 8);
    assert!(out.actuator_work_j > 0.);
    let mut cap = controls();
    cap.maximum_live_bytes = out.planned_live_bytes;
    let evaluate = |joints: &[DrivenJoint], cap| {
        advance_driven_joint_graph(&s, &p, joints, free(), &[Wrench::ZERO; 9], 0.1, cap, || {
            false
        })
    };
    assert_eq!(evaluate(&joints, cap).unwrap(), out);
    cap.maximum_live_bytes -= 1;
    assert!(matches!(
        evaluate(&joints, cap),
        Err(DistanceGraphError::Limit {
            resource: "live bytes",
            ..
        })
    ));
    cap = controls();
    cap.maximum_constraints = 31;
    assert!(matches!(
        evaluate(&joints, cap),
        Err(DistanceGraphError::Limit {
            resource: "constraints",
            ..
        })
    ));
    cap = controls();
    cap.maximum_trial_evaluations = 0;
    assert!(matches!(
        evaluate(&joints, cap),
        Err(DistanceGraphError::Limit {
            resource: "trial evaluations",
            ..
        })
    ));
    assert!(matches!(
        evaluate(&[joints[0], joints[0]], controls()),
        Err(DistanceGraphError::DependentConstraints { .. })
    ));
    joints.push(joints[6]);
    assert!(matches!(
        evaluate(&joints, controls()),
        Err(DistanceGraphError::Limit {
            resource: "constraints",
            ..
        })
    ));

    // All commands use the existing interior/representation admission. The
    // failing actuator is after a five-row hinge, not authored joint index 1.
    for (length, rate, acceleration, duration) in [
        (1., -4., 8., 1.),
        (1., -2., 0., 1.),
        (2_f64.powi(20), 2_f64.powi(-40), 0., 1.),
        (1., f64::NAN, 0., 1.),
        (1., 1e308, 1e308, 1e308),
    ] {
        let pair = [
            joints[0],
            DrivenJoint::Distance {
                link: center_link(6, 7, length),
                drive: DistanceDrive {
                    range_rate_mps: rate,
                    range_acceleration_mps2: acceleration,
                },
            },
        ];
        assert!(matches!(
            advance_driven_joint_graph(
                &s,
                &p,
                &pair,
                free(),
                &[Wrench::ZERO; 9],
                duration,
                controls(),
                || false
            ),
            Err(DistanceGraphError::PrescribedMotion { constraint: 5, .. })
        ));
    }
}

#[test]
fn g4_every_mixed_actuator_cancellation_is_atomic_and_retry_is_identical() {
    let (p, s, joints) = full_tensor();
    let mut calls = 0;
    let reference = advance_driven_joint_graph(
        &s,
        &p,
        &joints,
        free(),
        &[Wrench::ZERO; 4],
        0.02,
        controls(),
        || {
            calls += 1;
            false
        },
    )
    .unwrap();
    assert!(calls > 10);
    for stop in 1..=calls {
        let mut current = 0;
        let result = advance_driven_joint_graph(
            &s,
            &p,
            &joints,
            free(),
            &[Wrench::ZERO; 4],
            0.02,
            controls(),
            || {
                current += 1;
                current == stop
            },
        );
        assert!(
            matches!(
                result,
                Err(DistanceGraphError::Step(DistanceStepError::Cancelled))
            ),
            "checkpoint {stop}/{calls}: {result:?}"
        );
        assert_eq!((p, s, joints), full_tensor());
    }
    assert_eq!(run(&s, &p, &joints, 0.02), reference);
}
