use fs_mbd::rigid_rattle::distance::graph::{
    ArticulatedGraphStep, ArticulatedJoint, DistanceGraphControls, DistanceGraphError,
    DistanceLink, GraphJoint, JointReactionImpulse, RevoluteJoint, SphericalJoint,
    advance_articulated_graph, advance_joint_graph,
};
use fs_mbd::rigid_rattle::distance::{DistanceStepControls, DistanceStepError};
use fs_mbd::{
    Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState, UnitQuaternion, Vec3,
    Wrench,
};

const Z: Vec3 = Vec3::new(0.0, 0.0, 1.0);

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}, tolerance {tolerance:e}"
    );
}
fn vector(a: Vec3, b: Vec3, tolerance: f64) {
    near(a.x, b.x, tolerance);
    near(a.y, b.y, tolerance);
    near(a.z, b.z, tolerance);
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
            maximum_impulse_ns: 20.0,
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
fn hinge(bodies: [usize; 2], points: [Vec3; 2], axes: [Vec3; 2], length: f64) -> ArticulatedJoint {
    ArticulatedJoint::Revolute(RevoluteJoint::new(bodies, points, axes, length).unwrap())
}
fn run(
    s: &[RigidBodyState],
    p: &[MassProperties],
    joints: &[ArticulatedJoint],
    h: f64,
) -> ArticulatedGraphStep {
    advance_articulated_graph(
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
fn total(s: &[RigidBodyState], p: &[MassProperties]) -> (Vec3, Vec3, f64) {
    s.iter().zip(p).fold(
        (Vec3::ZERO, Vec3::ZERO, 0.0),
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
fn state_near(a: RigidBodyState, b: RigidBodyState, tol: f64) {
    vector(a.pose().position_world(), b.pose().position_world(), tol);
    vector(a.linear_momentum_world(), b.linear_momentum_world(), tol);
    vector(a.angular_momentum_body(), b.angular_momentum_body(), tol);
    for axis in [Vec3::new(1., 0., 0.), Vec3::new(0., 1., 0.), Z] {
        vector(
            a.pose().orientation().rotate_body_to_world(axis),
            b.pose().orientation().rotate_body_to_world(axis),
            tol,
        );
    }
}
fn reverse(s: RigidBodyState) -> RigidBodyState {
    RigidBodyState::new(
        s.pose(),
        s.linear_momentum_world().scale(-1.),
        s.angular_momentum_body().scale(-1.),
    )
    .unwrap()
}

#[test]
fn g1_hinge_chain_shares_transverse_torque_and_force_but_leaves_axial_spin_free() {
    let p = [(1., 2.), (2., 3.), (4., 5.)]
        .map(|(m, i)| MassProperties::new(m, Vec3::ZERO, Vec3::new(i, i, i)).unwrap());
    let s = p.map(|p| {
        state(
            p,
            Vec3::ZERO,
            UnitQuaternion::IDENTITY,
            Vec3::ZERO,
            Vec3::ZERO,
        )
    });
    let joints = [
        hinge([0, 1], [Vec3::ZERO; 2], [Z; 2], 0.7),
        hinge([1, 2], [Vec3::ZERO; 2], [Z; 2], 1.3),
    ];
    let mut w = [Wrench::ZERO; 3];
    w[0] = Wrench {
        force_world: Vec3::new(7., -14., 21.),
        torque_body: Vec3::new(10., 0., 0.),
    };
    let h = 0.02;
    let out =
        advance_articulated_graph(&s, &p, &joints, free(), &w, h, controls(), || false).unwrap();
    for i in 0..3 {
        let motion = out.step.states_after[i]
            .point_kinematics(p[i], Vec3::ZERO)
            .unwrap();
        vector(
            motion.center_of_mass_velocity_world,
            Vec3::new(1., -2., 3.).scale(h),
            2e-13,
        );
        vector(motion.angular_velocity_world, Vec3::new(h, 0., 0.), 2e-13);
        assert_eq!(
            out.step.states_after[i].pose(),
            out.step.free_steps[i].state_after.pose()
        );
    }
    assert_eq!(
        out.joint_rows
            .iter()
            .map(|r| (r.first_row, r.row_count))
            .collect::<Vec<_>>(),
        [(0, 5), (5, 5)]
    );
    for events in [&out.step.position_impulses, &out.step.velocity_impulses] {
        for (span, expected) in [(3, 0.08), (8, 0.05)] {
            let couple = events[span..span + 2]
                .iter()
                .fold(Vec3::ZERO, |sum, event| {
                    let JointReactionImpulse::Couple(event) = event else {
                        panic!("native couple required")
                    };
                    vector(event.first.linear_impulse_world, Vec3::ZERO, 0.);
                    vector(event.second.linear_impulse_world, Vec3::ZERO, 0.);
                    sum.add(event.couple_on_second_world_nms)
                });
            vector(couple, Vec3::new(expected, 0., 0.), 2e-11);
        }
    }
    // Free hinge spin is a different motion from a welded attitude.
    let spin = [0.7, -0.4, 1.2];
    let spinning = std::array::from_fn::<_, 3, _>(|i| {
        state(
            p[i],
            Vec3::ZERO,
            UnitQuaternion::IDENTITY,
            Vec3::ZERO,
            Z.scale(spin[i]),
        )
    });
    let out = run(&spinning, &p, &joints, h);
    for i in 0..3 {
        assert_eq!(
            out.step.states_after[i],
            free()
                .step(spinning[i], p[i], Wrench::ZERO, h)
                .unwrap()
                .state_after
        );
    }
}

fn gyros() -> (
    [MassProperties; 2],
    [RigidBodyState; 2],
    [ArticulatedJoint; 1],
) {
    let p = [
        MassProperties::new(2., Vec3::ZERO, Vec3::new(2., 2., 3.)).unwrap(),
        MassProperties::new(5., Vec3::ZERO, Vec3::new(3., 3., 4.)).unwrap(),
    ];
    let s = std::array::from_fn(|i| {
        state(
            p[i],
            Vec3::ZERO,
            UnitQuaternion::IDENTITY,
            Vec3::ZERO,
            Vec3::new(0.7, -0.4, [1.3, -0.7][i]),
        )
    });
    (p, s, [hinge([0, 1], [Vec3::ZERO; 2], [Z; 2], 1.)])
}

#[test]
fn g1_two_gyroscopic_rotors_match_independent_continuous_hinge_solution_and_refine() {
    // With coincident COMs and axial inertias (A_i,A_i,C_i), a hinge has
    // H = (A0+A1) n cross n_dot + S n, S=C0*w0+C1*w1 constant.
    // Hence n_dot=(H cross n)/(A0+A1), and each attitude is independently
    // R_H(t) R_z((w_i-S/(A0+A1))*t). This is not another constrained solve.
    let (p, initial, joints) = gyros();
    let transverse = 5.;
    let axial = 3. * 1.3 + 4. * (-0.7);
    let angular = Vec3::new(3.5, -2., axial);
    let duration = 0.6;
    let precession =
        UnitQuaternion::from_axis_angle(angular, magnitude(angular) * duration / transverse)
            .unwrap();
    let exact_q = [1.3, -0.7].map(|spin| {
        compose(
            precession,
            UnitQuaternion::from_axis_angle(Z, (spin - axial / transverse) * duration).unwrap(),
        )
    });
    let errors = [24, 48, 96].map(|count| {
        let mut s = initial.to_vec();
        for _ in 0..count {
            s = run(&s, &p, &joints, duration / f64::from(count))
                .step
                .states_after;
        }
        let mut error = 0.;
        for i in 0..2 {
            for axis in [Vec3::new(1., 0., 0.), Z] {
                error += magnitude(
                    s[i].pose()
                        .orientation()
                        .rotate_body_to_world(axis)
                        .sub(exact_q[i].rotate_body_to_world(axis)),
                );
            }
            near(
                s[i].angular_momentum_body().z,
                initial[i].angular_momentum_body().z,
                2e-11,
            );
        }
        vector(total(&s, &p).1, angular, 2e-11);
        error
    });
    for pair in errors.windows(2) {
        assert!(
            (3.7..4.3).contains(&(pair[0] / pair[1])),
            "errors {errors:?}"
        );
    }
    let mut s = initial.to_vec();
    for _ in 0..48 {
        s = run(&s, &p, &joints, 0.0125).step.states_after;
    }
    s = s.into_iter().map(reverse).collect();
    for _ in 0..48 {
        s = run(&s, &p, &joints, 0.0125).step.states_after;
    }
    for i in 0..2 {
        state_near(reverse(s[i]), initial[i], 4e-10);
    }
}

fn mixed() -> (
    [MassProperties; 4],
    [RigidBodyState; 4],
    [ArticulatedJoint; 3],
) {
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
    let arms = std::array::from_fn(|i| q[i].rotate_world_to_body(hpoint.sub(x[i])));
    let joints = [
        hinge(
            [0, 1],
            arms,
            [
                q[0].rotate_world_to_body(axis),
                q[1].rotate_world_to_body(axis),
            ],
            0.8,
        ),
        ArticulatedJoint::Spherical(SphericalJoint {
            bodies: [1, 2],
            points_body_m: [
                q[1].rotate_world_to_body(bpoint.sub(x[1])),
                q[2].rotate_world_to_body(bpoint.sub(x[2])),
            ],
        }),
        ArticulatedJoint::Distance(DistanceLink {
            bodies: [2, 3],
            points_body_m: [Vec3::ZERO; 2],
            length_m: magnitude(x[3].sub(x[2])),
        }),
    ];
    (p, s, joints)
}

#[test]
fn g0_off_center_full_tensor_hinge_ball_rod_graph_retains_momentum_work_and_couples() {
    let (p, initial, joints) = mixed();
    let before = total(&initial, &p);
    let mut s = initial.to_vec();
    let mut work = 0.;
    let mut angular = Vec3::ZERO;
    let mut largest_couple: f64 = 0.;
    for _ in 0..160 {
        let out = run(&s, &p, &joints, 0.005);
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
                    near(event.work_energy_residual_j, 0., 2e-15);
                }
                JointReactionImpulse::Couple(event) => {
                    vector(event.first.linear_increment_world, Vec3::ZERO, 0.);
                    vector(event.second.linear_increment_world, Vec3::ZERO, 0.);
                    vector(event.angular_momentum_change_world_nms, Vec3::ZERO, 2e-15);
                    near(event.work_energy_residual_j, 0., 2e-15);
                    largest_couple =
                        largest_couple.max(magnitude(event.couple_on_second_world_nms));
                }
            }
        }
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
        let free_change: f64 = out
            .step
            .free_steps
            .iter()
            .map(|step| {
                step.diagnostics_after.mechanical_energy - step.diagnostics_before.mechanical_energy
            })
            .sum();
        work += out.constraint_work_j + free_change;
        angular = angular.add(out.constraint_angular_impulse_world_nms);
        s = out.step.states_after;
    }
    let after = total(&s, &p);
    vector(after.0, before.0, 3e-12);
    vector(after.1.sub(before.1), angular, 3e-11);
    near(after.2 - before.2, work, 4e-12);
    assert!(
        largest_couple > 1e-5,
        "hinge must actually transmit a couple"
    );
    assert_eq!(
        run(&initial, &p, &joints, 0.01)
            .joint_rows
            .iter()
            .map(|r| (r.first_row, r.row_count))
            .collect::<Vec<_>>(),
        [(0, 5), (5, 3), (8, 1)]
    );
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

#[test]
fn g3_hinge_mixed_graph_is_frame_boost_endpoint_and_joint_order_covariant() {
    let (p, s, joints) = mixed();
    let h = 0.03;
    let baseline = run(&s, &p, &joints, h);
    let rotation = UnitQuaternion::from_axis_angle(Vec3::new(1., 2., -1.), 0.8).unwrap();
    let shift = Vec3::new(3., -2., 4.);
    for boost in [Vec3::ZERO, Vec3::new(64., -32., 16.)] {
        let moved =
            std::array::from_fn::<_, 4, _>(|i| transform(s[i], p[i], rotation, shift, boost));
        let out = run(&moved, &p, &joints, h);
        for i in 0..4 {
            state_near(
                out.step.states_after[i],
                transform(
                    baseline.step.states_after[i],
                    p[i],
                    rotation,
                    shift.add(boost.scale(h)),
                    boost,
                ),
                6e-11,
            );
        }
        near(out.constraint_work_j, baseline.constraint_work_j, 2e-12);
    }
    let reversed = joints.map(|joint| match joint {
        ArticulatedJoint::Distance(mut joint) => {
            joint.bodies.reverse();
            joint.points_body_m.reverse();
            ArticulatedJoint::Distance(joint)
        }
        ArticulatedJoint::Spherical(mut joint) => {
            joint.bodies.reverse();
            joint.points_body_m.reverse();
            ArticulatedJoint::Spherical(joint)
        }
        ArticulatedJoint::Revolute(joint) => {
            let mut bodies = joint.bodies();
            bodies.reverse();
            let mut points = joint.points_body_m();
            points.reverse();
            let mut axes = joint.axes_body();
            axes.reverse();
            hinge(bodies, points, axes, joint.angular_length_m())
        }
        _ => unreachable!("this endpoint-swap fixture uses rods, balls and hinges"),
    });
    for joints in [[joints[2], joints[0], joints[1]], reversed] {
        let out = run(&s, &p, &joints, h);
        for i in 0..4 {
            state_near(
                out.step.states_after[i],
                baseline.step.states_after[i],
                3e-11,
            );
        }
        near(out.constraint_work_j, baseline.constraint_work_j, 2e-12);
    }
}

#[test]
fn g3_axis_normalization_and_explicit_angular_scales_preserve_the_same_motion() {
    let (p, s, joints) = gyros();
    let baseline = run(&s, &p, &joints, 0.04);
    for axis_scale in [1e-300, 1e300] {
        let joint = hinge([0, 1], [Vec3::ZERO; 2], [Z.scale(axis_scale); 2], 1.);
        assert_eq!(run(&s, &p, &[joint], 0.04), baseline);
    }
    for length in [0.01, 100.] {
        let joint = hinge([0, 1], [Vec3::ZERO; 2], [Z; 2], length);
        let mut c = controls();
        c.position_tolerance_m *= length;
        c.velocity_tolerance_mps *= length;
        c.step.maximum_impulse_ns /= length;
        let out = advance_articulated_graph(
            &s,
            &p,
            &[joint],
            free(),
            &[Wrench::ZERO; 2],
            0.04,
            c,
            || false,
        )
        .unwrap();
        for i in 0..2 {
            state_near(
                out.step.states_after[i],
                baseline.step.states_after[i],
                3e-11,
            );
        }
        near(out.constraint_work_j, baseline.constraint_work_j, 3e-13);
    }
}

#[test]
fn g4_hinge_group_caps_rank_scalar_boundary_and_native_payload_limits() {
    let p = [MassProperties::new(1., Vec3::ZERO, Vec3::new(1., 1., 1.)).unwrap(); 2];
    let s = p.map(|p| {
        state(
            p,
            Vec3::ZERO,
            UnitQuaternion::IDENTITY,
            Vec3::ZERO,
            Vec3::ZERO,
        )
    });
    let joint = hinge([0, 1], [Vec3::ZERO; 2], [Z; 2], 1.);
    let w = [
        Wrench {
            force_world: Vec3::ZERO,
            torque_body: Vec3::new(1., 1., 0.),
        },
        Wrench::ZERO,
    ];
    let out =
        advance_articulated_graph(&s, &p, &[joint], free(), &w, 0.1, controls(), || false).unwrap();
    let mut cap = controls();
    cap.step.maximum_impulse_ns = 0.03;
    // Each initial scalar couple is 0.025 Nms, but the combined vector exceeds 0.03.
    for event in &out.step.position_impulses[3..5] {
        let JointReactionImpulse::Couple(event) = event else {
            unreachable!()
        };
        assert!(magnitude(event.couple_on_second_world_nms) < 0.03);
    }
    assert!(matches!(
        advance_articulated_graph(&s, &p, &[joint], free(), &w, 0.1, cap, || false),
        Err(DistanceGraphError::Step(
            DistanceStepError::LocalBranch | DistanceStepError::NotConverged
        ))
    ));
    let exact = DistanceGraphControls {
        maximum_live_bytes: out.planned_live_bytes,
        maximum_trial_evaluations: out.step.trial_evaluations,
        ..controls()
    };
    assert_eq!(
        advance_articulated_graph(&s, &p, &[joint], free(), &w, 0.1, exact, || false).unwrap(),
        out
    );
    for c in [
        DistanceGraphControls {
            maximum_live_bytes: out.planned_live_bytes - 1,
            ..exact
        },
        DistanceGraphControls {
            maximum_trial_evaluations: out.step.trial_evaluations - 1,
            ..exact
        },
        DistanceGraphControls {
            maximum_constraints: 4,
            ..exact
        },
    ] {
        assert!(matches!(
            advance_articulated_graph(&s, &p, &[joint], free(), &w, 0.1, c, || false),
            Err(DistanceGraphError::Limit { .. })
        ));
    }
    assert!(matches!(
        advance_articulated_graph(
            &s,
            &p,
            &[joint, joint],
            free(),
            &[Wrench::ZERO; 2],
            0.1,
            controls(),
            || false
        ),
        Err(DistanceGraphError::DependentConstraints { .. })
    ));
    let p = vec![p[0]; 9];
    let s: Vec<_> = (0..9)
        .map(|i| {
            state(
                p[i],
                Vec3::new(if i > 6 { (i - 6) as f64 } else { 0. }, 0., 0.),
                UnitQuaternion::IDENTITY,
                Vec3::ZERO,
                Vec3::ZERO,
            )
        })
        .collect();
    let mut joints: Vec<_> = (0..6)
        .map(|i| hinge([i, i + 1], [Vec3::ZERO; 2], [Z; 2], 1.))
        .collect();
    for i in 6..8 {
        joints.push(ArticulatedJoint::Distance(DistanceLink {
            bodies: [i, i + 1],
            points_body_m: [Vec3::ZERO; 2],
            length_m: 1.,
        }));
    }
    let mut w = vec![Wrench::ZERO; 9];
    w[0].force_world = Vec3::new(1., 0., 0.);
    w[0].torque_body = Vec3::new(0.2, 0.1, 0.);
    let out = advance_articulated_graph(&s, &p, &joints, free(), &w, 0.001, controls(), || false)
        .unwrap();
    assert_eq!(out.step.residuals.len(), 32);
    vector(
        out.step.states_after[8].linear_momentum_world(),
        Vec3::new(0.001 / 9., 0., 0.),
        2e-13,
    );
    joints[6] = hinge([6, 7], [Vec3::ZERO; 2], [Z; 2], 1.);
    assert!(matches!(
        advance_articulated_graph(&s, &p, &joints, free(), &w, 0.001, controls(), || false),
        Err(DistanceGraphError::Limit {
            resource: "constraints",
            ..
        })
    ));
}

#[test]
fn g0_hinge_admission_rejects_opposed_axes_and_retains_local_spin_at_large_boost() {
    let (p, s, joints) = gyros();
    let q = UnitQuaternion::from_axis_angle(Vec3::new(1., 0., 0.), std::f64::consts::PI).unwrap();
    let opposed = [s[0], state(p[1], Vec3::ZERO, q, Vec3::ZERO, Vec3::ZERO)];
    assert!(matches!(
        advance_articulated_graph(
            &opposed,
            &p,
            &joints,
            free(),
            &[Wrench::ZERO; 2],
            0.01,
            controls(),
            || false
        ),
        Err(DistanceGraphError::Step(DistanceStepError::LocalBranch))
    ));
    let mut c = controls();
    c.position_tolerance_m = 0.1;
    let q = UnitQuaternion::from_axis_angle(Vec3::new(1., 1., 0.), 0.113).unwrap();
    let tilted = [
        state(p[0], Vec3::ZERO, q, Vec3::ZERO, Vec3::ZERO),
        state(
            p[1],
            Vec3::ZERO,
            UnitQuaternion::IDENTITY,
            Vec3::ZERO,
            Vec3::ZERO,
        ),
    ];
    assert!(matches!(
        advance_articulated_graph(
            &tilted,
            &p,
            &joints,
            free(),
            &[Wrench::ZERO; 2],
            0.01,
            c,
            || false
        ),
        Err(DistanceGraphError::InitialConstraint { .. })
    ));
    let far = (1_u64 << 54) as f64;
    let states = std::array::from_fn::<_, 2, _>(|i| {
        state(
            p[i],
            Vec3::new(far, 0., 0.),
            UnitQuaternion::IDENTITY,
            Vec3::new(far, 0., 0.),
            Vec3::new(if i == 0 { 0.2 } else { 0. }, 0., 0.),
        )
    });
    assert_eq!(
        states[0]
            .point_kinematics(p[0], Vec3::ZERO)
            .unwrap()
            .point_velocity_world,
        states[1]
            .point_kinematics(p[1], Vec3::ZERO)
            .unwrap()
            .point_velocity_world
    );
    assert!(matches!(
        advance_articulated_graph(
            &states,
            &p,
            &joints,
            free(),
            &[Wrench::ZERO; 2],
            0.01,
            controls(),
            || false
        ),
        Err(DistanceGraphError::InitialConstraint { .. })
    ));
    for axes in [[Vec3::ZERO, Z], [Vec3::new(f64::NAN, 0., 0.), Z]] {
        assert!(RevoluteJoint::new([0, 1], [Vec3::ZERO; 2], axes, 1.).is_err());
    }
    for length in [0., -1., f64::INFINITY] {
        assert!(RevoluteJoint::new([0, 1], [Vec3::ZERO; 2], [Z; 2], length).is_err());
    }
    assert!(RevoluteJoint::new([0, 0], [Vec3::ZERO; 2], [Z; 2], 1.).is_err());
    assert!(
        advance_articulated_graph(
            &s,
            &p,
            &[hinge(
                [0, 1],
                [Vec3::ZERO; 2],
                [Z; 2],
                controls().position_tolerance_m
            )],
            free(),
            &[Wrench::ZERO; 2],
            0.01,
            controls(),
            || false
        )
        .is_err()
    );
}

#[test]
fn g1_point_only_articulation_keeps_the_original_native_ball_and_rod_events() {
    let (p, s, joints) = mixed();
    let ArticulatedJoint::Spherical(ball) = joints[1] else {
        unreachable!()
    };
    let ArticulatedJoint::Distance(rod) = joints[2] else {
        unreachable!()
    };
    let original = [GraphJoint::Spherical(ball), GraphJoint::Distance(rod)];
    let old = advance_joint_graph(
        &s,
        &p,
        &original,
        free(),
        &[Wrench::ZERO; 4],
        0.02,
        controls(),
        || false,
    )
    .unwrap();
    let new = run(&s, &p, &joints[1..], 0.02);
    assert_eq!(new.step.states_after, old.step.states_after);
    assert_eq!(new.step.free_steps, old.step.free_steps);
    assert_eq!(new.step.residuals, old.step.residuals);
    assert_eq!(new.step.constraint_work_j, old.step.constraint_work_j);
    assert_eq!(new.constraint_work_j, old.constraint_work_j);
    assert_eq!(new.joint_rows, old.joint_rows);
    for (old, new) in old
        .step
        .position_impulses
        .iter()
        .chain(&old.step.velocity_impulses)
        .zip(
            new.step
                .position_impulses
                .iter()
                .chain(&new.step.velocity_impulses),
        )
    {
        assert_eq!(*new, JointReactionImpulse::Point(*old));
    }
    let empty = run(&s, &p, &[], 0.02);
    for i in 0..4 {
        assert_eq!(
            empty.step.states_after[i],
            free()
                .step(s[i], p[i], Wrench::ZERO, 0.02)
                .unwrap()
                .state_after
        );
    }
}

#[test]
fn g4_every_observed_hinge_cancellation_preserves_inputs_and_retries_identically() {
    let (p, s, joints) = mixed();
    let original = s;
    let mut polls = 0;
    let baseline = advance_articulated_graph(
        &s,
        &p,
        &joints,
        free(),
        &[Wrench::ZERO; 4],
        0.03,
        controls(),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    assert!(polls > 100);
    for stop in 1..=polls {
        let mut seen = 0;
        let result = advance_articulated_graph(
            &s,
            &p,
            &joints,
            free(),
            &[Wrench::ZERO; 4],
            0.03,
            controls(),
            || {
                seen += 1;
                seen == stop
            },
        );
        assert!(
            matches!(
                result,
                Err(DistanceGraphError::Step(DistanceStepError::Cancelled))
            ),
            "stop {stop}: {result:?}"
        );
        assert_eq!(s, original);
    }
    assert_eq!(run(&s, &p, &joints, 0.03), baseline);
}
