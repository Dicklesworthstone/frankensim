use fs_mbd::rigid_rattle::distance::graph::{
    DistanceGraphControls, DistanceGraphError, DistanceLink, GraphJoint, JointGraphStep,
    SphericalJoint, advance_distance_graph, advance_joint_graph,
};
use fs_mbd::rigid_rattle::distance::{DistanceStepControls, DistanceStepError};
use fs_mbd::{
    Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState, UnitQuaternion, Vec3,
    Wrench,
};

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}; tolerance {tolerance:e}"
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
            maximum_impulse_ns: 10.,
            maximum_rotation_increment_rad: 0.25,
            maximum_position_iterations: 40,
            maximum_backtracks: 16,
        },
        position_tolerance_m: 2e-13,
        velocity_tolerance_mps: 2e-13,
        minimum_scaled_pivot: 1e-10,
        maximum_bodies: 64,
        maximum_constraints: 32,
        maximum_trial_evaluations: 160,
        maximum_live_bytes: 2_000_000,
    }
}
fn properties(mass: f64) -> MassProperties {
    MassProperties::new(mass, Vec3::ZERO, Vec3::new(0.4, 0.5, 0.6).scale(mass)).unwrap()
}
fn state(p: MassProperties, x: Vec3, q: UnitQuaternion, v: Vec3, w: Vec3) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(x, q).unwrap(),
        v.scale(p.mass()),
        p.angular_momentum_body_checked(w).unwrap(),
    )
    .unwrap()
}
fn resting(p: MassProperties, x: Vec3) -> RigidBodyState {
    state(p, x, UnitQuaternion::IDENTITY, Vec3::ZERO, Vec3::ZERO)
}
fn ball(a: usize, b: usize, arms: [Vec3; 2]) -> GraphJoint {
    GraphJoint::Spherical(SphericalJoint {
        bodies: [a, b],
        points_body_m: arms,
    })
}
fn run(
    s: &[RigidBodyState],
    p: &[MassProperties],
    joints: &[GraphJoint],
    h: f64,
) -> JointGraphStep {
    advance_joint_graph(
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
fn state_near(a: RigidBodyState, b: RigidBodyState, tol: f64) {
    vector(a.pose().position_world(), b.pose().position_world(), tol);
    vector(a.linear_momentum_world(), b.linear_momentum_world(), tol);
    vector(a.angular_momentum_body(), b.angular_momentum_body(), tol);
    for axis in [
        Vec3::new(1., 0., 0.),
        Vec3::new(0., 1., 0.),
        Vec3::new(0., 0., 1.),
    ] {
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

#[test]
fn g1_ball_chain_has_analytic_three_axis_recoil_without_welding_attitudes() {
    let p = [properties(1.), properties(2.), properties(3.)];
    let x = Vec3::new(0.4, -0.2, 0.3);
    let v = Vec3::new(0.1, -0.3, 0.2);
    let s = std::array::from_fn::<_, 3, _>(|i| {
        state(
            p[i],
            x,
            UnitQuaternion::IDENTITY,
            v,
            Vec3::new(0.3, -0.2, 0.1).scale(i as f64 + 1.),
        )
    });
    let w = [
        Vec3::new(7., 1., -3.),
        Vec3::new(-2., 4., 2.),
        Vec3::new(1., -2., 4.),
    ]
    .map(|force_world| Wrench {
        force_world,
        torque_body: Vec3::ZERO,
    });
    let joints = [ball(0, 1, [Vec3::ZERO; 2]), ball(1, 2, [Vec3::ZERO; 2])];
    let h = 0.04;
    let out = advance_joint_graph(&s, &p, &joints, free(), &w, h, controls(), || false).unwrap();
    let acceleration = Vec3::new(1., 0.5, 0.5);
    for i in 0..3 {
        let next = out.step.states_after[i];
        vector(
            next.pose().position_world(),
            x.add(v.scale(h)).add(acceleration.scale(0.5 * h * h)),
            3e-14,
        );
        vector(
            next.linear_momentum_world(),
            v.add(acceleration.scale(h)).scale(p[i].mass()),
            5e-14,
        );
        let native = free().step(s[i], p[i], w[i], h).unwrap();
        assert_eq!(
            next.pose().orientation(),
            native.state_after.pose().orientation()
        );
        assert_eq!(
            next.angular_momentum_body(),
            native.state_after.angular_momentum_body()
        );
        assert_eq!(next.pose(), out.step.free_steps[i].state_after.pose());
    }
    let reaction = [
        w[0].force_world
            .sub(acceleration.scale(p[0].mass()))
            .scale(h / 2.),
        acceleration
            .scale(p[2].mass())
            .sub(w[2].force_world)
            .scale(h / 2.),
    ];
    for events in [&out.step.position_impulses, &out.step.velocity_impulses] {
        for (j, expected) in reaction.into_iter().enumerate() {
            let impulse = events[3 * j..3 * j + 3]
                .iter()
                .fold(Vec3::ZERO, |sum, event| sum.add(event.second.impulse_world));
            vector(impulse, expected, 3e-14);
        }
    }
    assert_eq!(
        out.joint_rows
            .iter()
            .map(|r| (r.first_row, r.row_count))
            .collect::<Vec<_>>(),
        [(0, 3), (3, 3)]
    );
    let free_change: f64 = out
        .step
        .free_steps
        .iter()
        .map(|step| {
            step.diagnostics_after.mechanical_energy - step.diagnostics_before.mechanical_energy
        })
        .sum();
    near(
        totals(&out.step.states_after, &p).2 - totals(&s, &p).2,
        free_change + out.constraint_work_j,
        8e-15,
    );
    vector(out.constraint_angular_impulse_world_nms, Vec3::ZERO, 3e-16);
}

#[test]
fn g1_spinning_ball_pair_matches_independent_continuous_rotation_and_refines() {
    // m0*a = m1*b: the two inward joint forces supply both COM centripetal
    // accelerations and lie along the body arms, so the exact spins stay constant.
    let p = [properties(2.), properties(3.)];
    let x = [Vec3::new(-0.6, 0., 0.), Vec3::new(0.4, 0., 0.)];
    let omega = Vec3::new(0., 0., 0.7);
    let initial = std::array::from_fn::<_, 2, _>(|i| {
        state(
            p[i],
            x[i],
            UnitQuaternion::IDENTITY,
            omega.cross(x[i]),
            omega,
        )
    });
    let joints = [ball(0, 1, [x[0].scale(-1.), x[1].scale(-1.)])];
    let duration = 0.8;
    let exact_q =
        UnitQuaternion::from_axis_angle(Vec3::new(0., 0., 1.), omega.z * duration).unwrap();
    let errors = [20, 40, 80].map(|count| {
        let mut s = initial.to_vec();
        for _ in 0..count {
            s = run(&s, &p, &joints, duration / f64::from(count))
                .step
                .states_after;
        }
        let exact_x = exact_q.rotate_body_to_world(x[0]);
        magnitude(s[0].pose().position_world().sub(exact_x))
            + magnitude(
                s[0].linear_momentum_world()
                    .sub(omega.cross(exact_x).scale(p[0].mass())),
            )
    });
    for pair in errors.windows(2) {
        assert!(
            (3.7..4.3).contains(&(pair[0] / pair[1])),
            "errors {errors:?}"
        );
    }
    let mut s = initial.to_vec();
    for _ in 0..40 {
        s = run(&s, &p, &joints, 0.01).step.states_after;
    }
    let before = totals(&initial, &p);
    let after = totals(&s, &p);
    vector(after.0, before.0, 2e-13);
    vector(after.1, before.1, 2e-12);
    s = s.into_iter().map(reverse).collect();
    for _ in 0..40 {
        s = run(&s, &p, &joints, 0.01).step.states_after;
    }
    for i in 0..2 {
        state_near(reverse(s[i]), initial[i], 2e-10);
    }
}

fn mixed() -> ([MassProperties; 3], [RigidBodyState; 3], [GraphJoint; 2]) {
    let p = std::array::from_fn(|i| {
        MassProperties::from_principal_frame(
            1. + i as f64,
            Vec3::ZERO,
            Vec3::new(0.7, 1., 1.2).scale(1. + 0.2 * i as f64),
            UnitQuaternion::from_axis_angle(Vec3::new(1., 2., -1.), 0.3 + 0.2 * i as f64).unwrap(),
        )
        .unwrap()
    });
    let q = std::array::from_fn::<_, 3, _>(|i| {
        UnitQuaternion::from_axis_angle(Vec3::new(1., -0.5, 2.), 0.2 - 0.3 * i as f64).unwrap()
    });
    let arms = [Vec3::new(0.4, -0.2, 0.15), Vec3::new(-0.3, 0.1, -0.1)];
    let point = Vec3::new(0.1, 0.2, -0.1);
    let x = [
        point.sub(q[0].rotate_body_to_world(arms[0])),
        point.sub(q[1].rotate_body_to_world(arms[1])),
        Vec3::new(0.4, 1.2, 0.5),
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
    let rod_arms = [Vec3::new(0.1, 0.05, -0.15), Vec3::new(-0.2, -0.1, 0.1)];
    let rod = DistanceLink {
        bodies: [1, 2],
        points_body_m: rod_arms,
        length_m: magnitude(
            s[2].pose()
                .point_world_from_body(rod_arms[1])
                .unwrap()
                .sub(s[1].pose().point_world_from_body(rod_arms[0]).unwrap()),
        ),
    };
    (p, s, [ball(0, 1, arms), GraphJoint::Distance(rod)])
}

#[test]
fn g0_mixed_full_tensor_graph_keeps_real_momentum_work_and_relative_rotation() {
    let (p, initial, joints) = mixed();
    let before = totals(&initial, &p);
    let mut s = initial.to_vec();
    let mut work = 0.;
    let mut couple = Vec3::ZERO;
    for _ in 0..200 {
        let out = run(&s, &p, &joints, 0.005);
        for event in out
            .step
            .position_impulses
            .iter()
            .chain(&out.step.velocity_impulses)
        {
            vector(event.impulse_balance_world, Vec3::ZERO, 0.);
            vector(
                event.angular_momentum_change_world,
                event.angular_impulse_world,
                8e-16,
            );
            near(event.work_energy_residual_j, 0., 5e-16);
        }
        for i in 0..3 {
            assert_eq!(
                out.step.states_after[i].pose(),
                out.step.free_steps[i].state_after.pose()
            );
        }
        let r = &out.step.residuals;
        assert!(r[0].gap_m.hypot(r[1].gap_m).hypot(r[2].gap_m) <= controls().position_tolerance_m);
        assert!(
            r[0].range_rate_mps
                .hypot(r[1].range_rate_mps)
                .hypot(r[2].range_rate_mps)
                <= controls().velocity_tolerance_mps
        );
        work += out.constraint_work_j;
        couple = couple.add(out.constraint_angular_impulse_world_nms);
        s = out.step.states_after;
    }
    let after = totals(&s, &p);
    vector(after.0, before.0, 2e-12);
    vector(after.1.sub(before.1), couple, 2e-11);
    near(after.2 - before.2, work, 3e-12);
    assert!(
        work.abs() > 1e-10,
        "finite-step work stays visible: {work:e}"
    );
    let old_omega = initial[0]
        .point_kinematics(p[0], Vec3::ZERO)
        .unwrap()
        .angular_velocity_world
        .sub(
            initial[1]
                .point_kinematics(p[1], Vec3::ZERO)
                .unwrap()
                .angular_velocity_world,
        );
    let new_omega = s[0]
        .point_kinematics(p[0], Vec3::ZERO)
        .unwrap()
        .angular_velocity_world
        .sub(
            s[1].point_kinematics(p[1], Vec3::ZERO)
                .unwrap()
                .angular_velocity_world,
        );
    assert!(
        magnitude(new_omega.sub(old_omega)) > 1e-3,
        "ball joint must leave relative rotation free"
    );
    assert_eq!(
        run(&initial, &p, &joints, 0.01)
            .joint_rows
            .iter()
            .map(|r| (r.first_row, r.row_count))
            .collect::<Vec<_>>(),
        [(0, 3), (3, 1)]
    );
}

#[test]
fn g3_mixed_motion_is_frame_boost_joint_order_and_endpoint_covariant() {
    let (p, s, joints) = mixed();
    let h = 0.03;
    let reference = run(&s, &p, &joints, h);
    let rotation = UnitQuaternion::from_axis_angle(Vec3::new(1., 2., -1.), 0.8).unwrap();
    let shift = Vec3::new(3., -2., 4.);
    for boost in [Vec3::ZERO, Vec3::new(64., -32., 16.)] {
        let moved =
            std::array::from_fn::<_, 3, _>(|i| transform(s[i], p[i], rotation, shift, boost));
        let out = run(&moved, &p, &joints, h);
        for i in 0..3 {
            state_near(
                out.step.states_after[i],
                transform(
                    reference.step.states_after[i],
                    p[i],
                    rotation,
                    shift.add(boost.scale(h)),
                    boost,
                ),
                3e-11,
            );
        }
        near(out.constraint_work_j, reference.constraint_work_j, 5e-13);
    }
    let reversed = joints.map(|joint| match joint {
        GraphJoint::Distance(mut link) => {
            link.bodies.reverse();
            link.points_body_m.reverse();
            GraphJoint::Distance(link)
        }
        GraphJoint::Spherical(mut link) => {
            link.bodies.reverse();
            link.points_body_m.reverse();
            GraphJoint::Spherical(link)
        }
    });
    for variant in [[joints[1], joints[0]], reversed] {
        let out = run(&s, &p, &variant, h);
        for i in 0..3 {
            state_near(
                out.step.states_after[i],
                reference.step.states_after[i],
                3e-11,
            );
        }
        near(out.constraint_work_j, reference.constraint_work_j, 5e-13);
    }
    let order = [2, 0, 1];
    let inverse = [1, 2, 0];
    let reordered = joints.map(|joint| match joint {
        GraphJoint::Distance(mut link) => {
            link.bodies = link.bodies.map(|i| inverse[i]);
            GraphJoint::Distance(link)
        }
        GraphJoint::Spherical(mut link) => {
            link.bodies = link.bodies.map(|i| inverse[i]);
            GraphJoint::Spherical(link)
        }
    });
    let out = run(&order.map(|i| s[i]), &order.map(|i| p[i]), &reordered, h);
    for i in 0..3 {
        state_near(
            out.step.states_after[inverse[i]],
            reference.step.states_after[i],
            3e-11,
        );
    }
}

#[test]
fn g1_distance_only_and_unlinked_bodies_retain_the_original_native_step() {
    let (p, s, joints) = mixed();
    let GraphJoint::Distance(rod) = joints[1] else {
        unreachable!()
    };
    for rods in [vec![], vec![rod]] {
        let old = advance_distance_graph(
            &s,
            &p,
            &rods,
            free(),
            &[Wrench::ZERO; 3],
            0.01,
            controls(),
            || false,
        )
        .unwrap();
        let out = run(
            &s,
            &p,
            &rods
                .into_iter()
                .map(GraphJoint::Distance)
                .collect::<Vec<_>>(),
            0.01,
        );
        assert_eq!(out.step, old);
        near(out.constraint_work_j, old.constraint_work_j, 3e-16);
    }
}

#[test]
fn g0_cartesian_admission_keeps_local_arms_and_spin_at_remote_origins_and_boosts() {
    let p = [properties(1.); 2];
    let far = 2_f64.powi(54);
    let s = [resting(p[0], Vec3::new(far, 0., 0.)); 2];
    let joints = [ball(0, 1, [Vec3::ZERO, Vec3::new(1., 0., 0.)])];
    assert_eq!(
        s[0].pose().point_world_from_body(Vec3::ZERO).unwrap(),
        s[1].pose()
            .point_world_from_body(Vec3::new(1., 0., 0.))
            .unwrap()
    );
    assert!(matches!(
        advance_joint_graph(
            &s,
            &p,
            &joints,
            free(),
            &[Wrench::ZERO; 2],
            0.001,
            controls(),
            || false
        ),
        Err(DistanceGraphError::InitialConstraint { .. })
    ));
    let arms = [Vec3::new(0., 1., 0.); 2];
    let s = [
        state(
            p[0],
            Vec3::ZERO,
            UnitQuaternion::IDENTITY,
            Vec3::new(far, 0., 0.),
            Vec3::ZERO,
        ),
        state(
            p[1],
            Vec3::ZERO,
            UnitQuaternion::IDENTITY,
            Vec3::new(far, 0., 0.),
            Vec3::new(0., 0., 1.),
        ),
    ];
    assert_eq!(
        s[0].point_kinematics(p[0], arms[0])
            .unwrap()
            .point_velocity_world,
        s[1].point_kinematics(p[1], arms[1])
            .unwrap()
            .point_velocity_world
    );
    assert!(matches!(
        advance_joint_graph(
            &s,
            &p,
            &[ball(0, 1, arms)],
            free(),
            &[Wrench::ZERO; 2],
            0.001,
            controls(),
            || false
        ),
        Err(DistanceGraphError::InitialConstraint { .. })
    ));
    let s = [resting(p[0], Vec3::new(far, 0., 0.)); 2];
    assert_eq!(run(&s, &p, &[ball(0, 1, arms)], 0.001).step.states_after, s);
}

#[test]
fn g4_ball_norm_caps_rank_and_full_scalar_row_boundary_are_admitted_before_publication() {
    let p = [properties(1.); 2];
    let s = [resting(p[0], Vec3::ZERO); 2];
    let joint = ball(0, 1, [Vec3::ZERO; 2]);
    let w = [
        Wrench {
            force_world: Vec3::new(36., 36., 36.),
            torque_body: Vec3::ZERO,
        },
        Wrench::ZERO,
    ];
    let out = advance_joint_graph(&s, &p, &[joint], free(), &w, 0.1, controls(), || false).unwrap();
    let cap = DistanceGraphControls {
        step: DistanceStepControls {
            maximum_impulse_ns: 1.,
            ..controls().step
        },
        ..controls()
    };
    // Each Cartesian reaction is 0.9 Ns but the physical vector exceeds 1 Ns.
    assert!(matches!(
        advance_joint_graph(&s, &p, &[joint], free(), &w, 0.1, cap, || false),
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
        advance_joint_graph(&s, &p, &[joint], free(), &w, 0.1, exact, || false).unwrap(),
        out
    );
    for limits in [
        DistanceGraphControls {
            maximum_live_bytes: out.planned_live_bytes - 1,
            ..exact
        },
        DistanceGraphControls {
            maximum_trial_evaluations: out.step.trial_evaluations - 1,
            ..exact
        },
        DistanceGraphControls {
            maximum_constraints: 2,
            ..exact
        },
    ] {
        assert!(matches!(
            advance_joint_graph(&s, &p, &[joint], free(), &w, 0.1, limits, || false),
            Err(DistanceGraphError::Limit { .. })
        ));
    }
    assert!(matches!(
        advance_joint_graph(
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
    let tolerance = controls().position_tolerance_m;
    let displaced = [
        s[0],
        resting(
            p[1],
            Vec3::new(0.8 * tolerance, 0.8 * tolerance, 0.8 * tolerance),
        ),
    ];
    assert!(matches!(
        advance_joint_graph(
            &displaced,
            &p,
            &[joint],
            free(),
            &[Wrench::ZERO; 2],
            0.1,
            controls(),
            || false
        ),
        Err(DistanceGraphError::InitialConstraint { .. })
    ));
    // Ten three-row joints plus two rods form a connected 32-row graph.
    let p = vec![properties(1.); 13];
    let s: Vec<_> = (0..13)
        .map(|i| {
            resting(
                p[i],
                Vec3::new(if i > 10 { (i - 10) as f64 } else { 0. }, 0., 0.),
            )
        })
        .collect();
    let mut joints: Vec<_> = (0..10).map(|i| ball(i, i + 1, [Vec3::ZERO; 2])).collect();
    for i in 10..12 {
        joints.push(GraphJoint::Distance(DistanceLink {
            bodies: [i, i + 1],
            points_body_m: [Vec3::ZERO; 2],
            length_m: 1.,
        }));
    }
    assert_eq!(run(&s, &p, &joints, 0.001).step.residuals.len(), 32);
    joints[10] = ball(10, 11, [Vec3::ZERO; 2]);
    assert!(matches!(
        advance_joint_graph(
            &s,
            &p,
            &joints,
            free(),
            &vec![Wrench::ZERO; s.len()],
            0.001,
            controls(),
            || false
        ),
        Err(DistanceGraphError::Limit {
            resource: "constraints",
            ..
        })
    ));
}

#[test]
fn g4_every_observed_mixed_transaction_cancellation_is_atomic_and_retryable() {
    let (p, s, joints) = mixed();
    let mut polls = 0;
    let baseline = advance_joint_graph(
        &s,
        &p,
        &joints,
        free(),
        &[Wrench::ZERO; 3],
        0.03,
        controls(),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    assert!(polls > 40);
    for stop in 1..=polls {
        let mut count = 0;
        let result = advance_joint_graph(
            &s,
            &p,
            &joints,
            free(),
            &[Wrench::ZERO; 3],
            0.03,
            controls(),
            || {
                count += 1;
                count == stop
            },
        );
        assert_eq!(
            result,
            Err(DistanceGraphError::Step(DistanceStepError::Cancelled)),
            "poll {stop}/{polls}"
        );
    }
    assert_eq!(run(&s, &p, &joints, 0.03), baseline);
    for bad in [
        ball(0, 0, [Vec3::ZERO; 2]),
        ball(0, 3, [Vec3::ZERO; 2]),
        ball(0, 1, [Vec3::new(f64::NAN, 0., 0.), Vec3::ZERO]),
    ] {
        assert!(
            advance_joint_graph(
                &s,
                &p,
                &[bad],
                free(),
                &[Wrench::ZERO; 3],
                0.03,
                controls(),
                || false
            )
            .is_err()
        );
    }
    assert_eq!(run(&s, &p, &joints, 0.03), baseline);
}
