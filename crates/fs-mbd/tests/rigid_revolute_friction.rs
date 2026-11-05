use fs_mbd::rigid_rattle::distance::DistanceStepControls;
use fs_mbd::rigid_rattle::distance::graph::revolute_friction::{
    RevoluteFrictionControls, RevoluteFrictionError, RevoluteFrictionInput, RevoluteFrictionRegime,
    RevoluteFrictionReport, apply_revolute_friction_impulse,
};
use fs_mbd::rigid_rattle::distance::graph::{
    ArticulatedJoint, DistanceGraphControls, DistanceGraphError, DistanceLink, PrismaticJoint,
    RevoluteJoint, advance_articulated_graph,
};
use fs_mbd::{
    Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState, UnitQuaternion, Vec3,
    Wrench,
};

const Z: Vec3 = Vec3::new(0.0, 0.0, 1.0);

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tolerance {tolerance:e}"
    );
}
fn vector(actual: Vec3, expected: Vec3, tolerance: f64) {
    near(actual.x, expected.x, tolerance);
    near(actual.y, expected.y, tolerance);
    near(actual.z, expected.z, tolerance);
}
fn properties(mass: f64, inertia: f64) -> MassProperties {
    MassProperties::new(mass, Vec3::ZERO, Vec3::new(inertia, inertia, inertia)).unwrap()
}
fn state(property: MassProperties, position: Vec3, velocity: Vec3, spin: Vec3) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(position, UnitQuaternion::IDENTITY).unwrap(),
        velocity.scale(property.mass()),
        property.angular_momentum_body_checked(spin).unwrap(),
    )
    .unwrap()
}
fn hinge(bodies: [usize; 2], points: [Vec3; 2]) -> ArticulatedJoint {
    ArticulatedJoint::Revolute(RevoluteJoint::new(bodies, points, [Z; 2], 1.0).unwrap())
}
fn controls() -> RevoluteFrictionControls {
    RevoluteFrictionControls {
        graph: DistanceGraphControls {
            step: DistanceStepControls {
                maximum_impulse_ns: 100.0,
                maximum_rotation_increment_rad: 0.25,
                maximum_position_iterations: 30,
                maximum_backtracks: 16,
            },
            position_tolerance_m: 1e-10,
            velocity_tolerance_mps: 1e-10,
            minimum_scaled_pivot: 1e-10,
            maximum_bodies: 32,
            maximum_constraints: 32,
            maximum_trial_evaluations: 160,
            maximum_live_bytes: 2_000_000,
        },
        maximum_axial_impulse_nms: 100.0,
        angular_rate_tolerance_rad_s: 1e-10,
        work_tolerance_j: 1e-9,
        linear_impulse_tolerance_ns: 1e-9,
        angular_impulse_tolerance_nms: 1e-9,
    }
}
fn run(
    states: &[RigidBodyState],
    properties: &[MassProperties],
    joints: &[ArticulatedJoint],
    index: usize,
    input: RevoluteFrictionInput,
) -> RevoluteFrictionReport {
    apply_revolute_friction_impulse(states, properties, joints, index, input, controls(), || {
        false
    })
    .unwrap()
}
fn kinetic(capacity: f64) -> RevoluteFrictionInput {
    RevoluteFrictionInput::Sliding {
        kinetic_capacity_nms: capacity,
    }
}
fn audit(before: &[RigidBodyState], properties: &[MassProperties], out: &RevoluteFrictionReport) {
    let integrator = RigidBodyIntegrator::new(Gravity::ZERO);
    let mut energy = 0.0;
    let mut linear = Vec3::ZERO;
    let mut angular = Vec3::ZERO;
    for ((&before, &after), &p) in before.iter().zip(&out.states_after).zip(properties) {
        assert_eq!(before.pose(), after.pose());
        let a = integrator.diagnostics(before, p).unwrap();
        let b = integrator.diagnostics(after, p).unwrap();
        energy += b.mechanical_energy - a.mechanical_energy;
        linear = linear.add(b.linear_momentum_world.sub(a.linear_momentum_world));
        angular = angular.add(b.angular_momentum_world.sub(a.angular_momentum_world));
    }
    near(energy, out.kinetic_energy_change_j, 2e-12);
    near(
        energy,
        out.drive_work_j + out.friction_work_j + out.constraint_work_j,
        2e-12,
    );
    assert!(out.friction_work_j <= 0.0);
    near(out.work_energy_residual_j, 0.0, 2e-12);
    near(out.constraint_work_j, 0.0, 2e-12);
    vector(linear, Vec3::ZERO, 2e-12);
    vector(angular, Vec3::ZERO, 2e-12);
}

#[test]
fn g1_coaxial_rotors_share_spin_and_arrest_only_with_kinetic_capacity() {
    let p = [properties(1.0, 2.0), properties(2.0, 3.0)];
    let s = [
        state(p[0], Vec3::ZERO, Vec3::ZERO, Vec3::ZERO),
        state(p[1], Vec3::ZERO, Vec3::ZERO, Z.scale(5.0)),
    ];
    let joints = [hinge([0, 1], [Vec3::ZERO; 2])];
    for capacity in [0.0, 1.0, 5.5, 6.0, 10.0] {
        let out = run(&s, &p, &joints, 0, kinetic(capacity));
        let applied = capacity.min(6.0);
        near(out.inverse_effective_inertia, 5.0 / 6.0, 2e-15);
        near(out.friction_impulse_nms, -applied, 2e-15);
        near(
            out.relative_rate_after_rad_s,
            5.0 - 5.0 * applied / 6.0,
            2e-14,
        );
        near(
            out.friction_work_j,
            -applied * (5.0 - 5.0 * applied / 12.0),
            2e-14,
        );
        assert_eq!(
            out.regime,
            if capacity < 6.0 {
                RevoluteFrictionRegime::Sliding
            } else {
                RevoluteFrictionRegime::Arrested
            }
        );
        audit(&s, &p, &out);
        // The existing RATTLE owner consumes these same finite body states.
        let next = advance_articulated_graph(
            &out.states_after,
            &p,
            &joints,
            RigidBodyIntegrator::new(Gravity::ZERO),
            &[Wrench::ZERO; 2],
            0.001,
            controls().graph,
            || false,
        )
        .unwrap();
        near(
            next.step.states_after[1]
                .point_kinematics(p[1], Vec3::ZERO)
                .unwrap()
                .angular_velocity_world
                .z
                - next.step.states_after[0]
                    .point_kinematics(p[0], Vec3::ZERO)
                    .unwrap()
                    .angular_velocity_world
                    .z,
            out.relative_rate_after_rad_s,
            3e-13,
        );
    }
}

#[test]
fn g1_static_drive_holds_without_work_and_breaks_away_at_the_authored_bound() {
    let p = [properties(1.0, 2.0), properties(2.0, 3.0)];
    let s = p.map(|p| state(p, Vec3::ZERO, Vec3::ZERO, Vec3::ZERO));
    let joints = [hinge([0, 1], [Vec3::ZERO; 2])];
    for drive in [-1.0, -0.3, 0.0, 0.3, 1.0] {
        let out = run(
            &s,
            &p,
            &joints,
            0,
            RevoluteFrictionInput::RestingDrive {
                drive_impulse_nms: drive,
                static_capacity_nms: 1.0,
                kinetic_capacity_nms: 0.4,
            },
        );
        assert_eq!(out.states_after, s);
        assert_eq!(out.regime, RevoluteFrictionRegime::Held);
        assert_eq!(out.friction_impulse_nms, -drive);
        assert_eq!(out.friction_work_j, 0.0);
        assert_eq!(out.drive_work_j, 0.0);
        audit(&s, &p, &out);
    }
    for drive in [-1.2, 1.2] {
        let out = run(
            &s,
            &p,
            &joints,
            0,
            RevoluteFrictionInput::RestingDrive {
                drive_impulse_nms: drive,
                static_capacity_nms: 1.0,
                kinetic_capacity_nms: 0.4,
            },
        );
        near(out.friction_impulse_nms, -drive.signum() * 0.4, 0.0);
        near(
            out.relative_rate_after_rad_s,
            drive.signum() * 2.0 / 3.0,
            2e-15,
        );
        near(out.drive_work_j, 0.4, 2e-15);
        near(out.friction_work_j, -2.0 / 15.0, 2e-15);
        assert_eq!(out.regime, RevoluteFrictionRegime::Sliding);
        audit(&s, &p, &out);
    }
    // Even far below the numerical endpoint tolerance, real initial spin is
    // not made eligible for the larger static capacity by a capture threshold.
    let moving = [s[0], state(p[1], Vec3::ZERO, Vec3::ZERO, Z.scale(1e-14))];
    assert!(matches!(
        apply_revolute_friction_impulse(
            &moving,
            &p,
            &joints,
            0,
            RevoluteFrictionInput::RestingDrive {
                drive_impulse_nms: 0.1,
                static_capacity_nms: 1.0,
                kinetic_capacity_nms: 0.4
            },
            controls(),
            || false
        ),
        Err(RevoluteFrictionError::NotResting { .. })
    ));
}

#[test]
fn g1_off_center_hinge_matches_independent_reduced_kinetic_energy() {
    let p = [properties(2.0, 0.5), properties(3.0, 0.9)];
    let a = 0.8;
    let b = 0.2;
    let spin = [-0.3, 1.7];
    let point_speed = (2.0 * a * spin[0] - 3.0 * b * spin[1]) / 5.0;
    let s = [
        state(
            p[0],
            Vec3::new(-a, 0.0, 0.0),
            Vec3::new(0.0, point_speed - a * spin[0], 0.0),
            Z.scale(spin[0]),
        ),
        state(
            p[1],
            Vec3::new(b, 0.0, 0.0),
            Vec3::new(0.0, point_speed + b * spin[1], 0.0),
            Z.scale(spin[1]),
        ),
    ];
    let joints = [hinge(
        [0, 1],
        [Vec3::new(a, 0.0, 0.0), Vec3::new(-b, 0.0, 0.0)],
    )];
    let out = run(&s, &p, &joints, 0, kinetic(0.25));
    // Eliminate the common hinge speed at zero total linear momentum:
    // K = (I1*w1² + I2*w2² + reduced_mass*(a*w1+b*w2)²)/2.
    // Inverting this independent 2x2 angular mass gives the exact response.
    let reduced_mass = 6.0 / 5.0;
    let h00 = 0.5 + reduced_mass * a * a;
    let h11 = 0.9 + reduced_mass * b * b;
    let h01 = reduced_mass * a * b;
    let determinant = h00 * h11 - h01 * h01;
    let inverse = (h00 + h11 + 2.0 * h01) / determinant;
    near(out.inverse_effective_inertia, inverse, 3e-14);
    assert!((inverse - (1.0 / 0.5 + 1.0 / 0.9)).abs() > 0.5);
    let expected = [
        spin[0] + 0.25 * (h11 + h01) / determinant,
        spin[1] - 0.25 * (h00 + h01) / determinant,
    ];
    for i in 0..2 {
        near(
            out.states_after[i]
                .point_kinematics(p[i], Vec3::ZERO)
                .unwrap()
                .angular_velocity_world
                .z,
            expected[i],
            2e-14,
        );
    }
    let points = [Vec3::new(a, 0.0, 0.0), Vec3::new(-b, 0.0, 0.0)];
    vector(
        out.states_after[0]
            .point_kinematics(p[0], points[0])
            .unwrap()
            .point_velocity_world,
        out.states_after[1]
            .point_kinematics(p[1], points[1])
            .unwrap()
            .point_velocity_world,
        2e-14,
    );
    assert!(out.constraint_impulses_ns.iter().any(|p| p.abs() > 0.01));
    audit(&s, &p, &out);
}

#[test]
fn g3_other_native_joint_transmits_friction_to_a_third_rotor() {
    let p = [
        properties(1.0, 2.0),
        properties(2.0, 3.0),
        properties(3.0, 5.0),
    ];
    let s = [
        state(p[0], Vec3::ZERO, Vec3::ZERO, Vec3::ZERO),
        state(p[1], Vec3::ZERO, Vec3::ZERO, Z.scale(4.0)),
        state(p[2], Vec3::ZERO, Vec3::ZERO, Z.scale(4.0)),
    ];
    let h = hinge([0, 1], [Vec3::ZERO; 2]);
    let rail = ArticulatedJoint::Prismatic(
        PrismaticJoint::new([1, 2], [Vec3::ZERO; 2], [UnitQuaternion::IDENTITY; 2], 0.7).unwrap(),
    );
    let out = run(&s, &p, &[h, rail], 0, kinetic(1.0));
    near(out.inverse_effective_inertia, 1.0 / 2.0 + 1.0 / 8.0, 2e-15);
    for (i, expected) in [0.5, 3.875, 3.875].into_iter().enumerate() {
        near(
            out.states_after[i]
                .point_kinematics(p[i], Vec3::ZERO)
                .unwrap()
                .angular_velocity_world
                .z,
            expected,
            2e-14,
        );
    }
    audit(&s, &p, &out);
    let reordered = run(&s, &p, &[rail, h], 1, kinetic(1.0));
    for i in 0..3 {
        vector(
            reordered.states_after[i].angular_momentum_body(),
            out.states_after[i].angular_momentum_body(),
            3e-14,
        );
    }
    let swapped = run(
        &s,
        &p,
        &[hinge([1, 0], [Vec3::ZERO; 2]), rail],
        0,
        kinetic(1.0),
    );
    near(
        swapped.friction_impulse_nms,
        -out.friction_impulse_nms,
        2e-15,
    );
    for i in 0..3 {
        vector(
            swapped.states_after[i].angular_momentum_body(),
            out.states_after[i].angular_momentum_body(),
            3e-14,
        );
    }
}

#[test]
fn g3_full_tensor_response_is_covariant_under_frame_rotation_and_common_boost() {
    let p = [0, 1].map(|i| {
        MassProperties::from_principal_frame(
            1.0 + i as f64,
            Vec3::ZERO,
            Vec3::new(2.0, 3.0, 4.0),
            UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, -1.0), 0.3 + i as f64 * 0.4)
                .unwrap(),
        )
        .unwrap()
    });
    let s = [
        state(p[0], Vec3::ZERO, Vec3::ZERO, Z.scale(-0.5)),
        state(p[1], Vec3::ZERO, Vec3::ZERO, Z.scale(2.0)),
    ];
    let joints = [hinge([0, 1], [Vec3::ZERO; 2])];
    let out = run(&s, &p, &joints, 0, kinetic(0.4));
    audit(&s, &p, &out);
    let rotation = UnitQuaternion::from_axis_angle(Vec3::new(2.0, -1.0, 3.0), 0.7).unwrap();
    let boost = Vec3::new(100.0, -200.0, 50.0);
    let transformed = core::array::from_fn::<_, 2, _>(|i| {
        RigidBodyState::new(
            Pose::new(Vec3::ZERO, rotation).unwrap(),
            boost.scale(p[i].mass()),
            s[i].angular_momentum_body(),
        )
        .unwrap()
    });
    let shifted = run(&transformed, &p, &joints, 0, kinetic(0.4));
    near(
        shifted.inverse_effective_inertia,
        out.inverse_effective_inertia,
        2e-14,
    );
    near(shifted.friction_work_j, out.friction_work_j, 2e-14);
    near(
        shifted.kinetic_energy_change_j,
        out.kinetic_energy_change_j,
        2e-14,
    );
    for i in 0..2 {
        vector(
            shifted.states_after[i].angular_momentum_body(),
            out.states_after[i].angular_momentum_body(),
            4e-14,
        );
        vector(
            shifted.states_after[i].linear_momentum_world(),
            boost.scale(p[i].mass()),
            0.0,
        );
    }
}

#[test]
fn g4_locked_axis_limits_and_every_cancellation_publish_no_partial_state() {
    let p = [properties(1.0, 2.0), properties(2.0, 3.0)];
    let s = [
        state(p[0], Vec3::ZERO, Vec3::ZERO, Vec3::ZERO),
        state(p[1], Vec3::ZERO, Vec3::ZERO, Z.scale(3.0)),
    ];
    let joints = [hinge([0, 1], [Vec3::ZERO; 2])];
    let initial = s;
    let mut checkpoints = 0;
    let baseline =
        apply_revolute_friction_impulse(&s, &p, &joints, 0, kinetic(0.3), controls(), || {
            checkpoints += 1;
            false
        })
        .unwrap();
    for cancel_at in 1..=checkpoints {
        let mut calls = 0;
        assert!(
            apply_revolute_friction_impulse(&s, &p, &joints, 0, kinetic(0.3), controls(), || {
                calls += 1;
                calls == cancel_at
            })
            .is_err()
        );
        assert_eq!(s, initial);
    }
    assert_eq!(run(&s, &p, &joints, 0, kinetic(0.3)), baseline);
    for malformed in [
        kinetic(-1.0),
        kinetic(f64::NAN),
        RevoluteFrictionInput::RestingDrive {
            drive_impulse_nms: 0.2,
            static_capacity_nms: 0.1,
            kinetic_capacity_nms: 0.3,
        },
    ] {
        assert!(
            apply_revolute_friction_impulse(&s, &p, &joints, 0, malformed, controls(), || false)
                .is_err()
        );
    }
    let mut bounded = controls();
    bounded.graph.maximum_live_bytes = baseline.planned_live_bytes;
    assert!(
        apply_revolute_friction_impulse(&s, &p, &joints, 0, kinetic(0.3), bounded, || false)
            .is_ok()
    );
    bounded.graph.maximum_live_bytes -= 1;
    assert!(matches!(
        apply_revolute_friction_impulse(&s, &p, &joints, 0, kinetic(0.3), bounded, || false),
        Err(RevoluteFrictionError::Native(
            DistanceGraphError::Limit { .. }
        ))
    ));
    let mut bounded = controls();
    bounded.maximum_axial_impulse_nms = 0.1;
    assert!(
        apply_revolute_friction_impulse(&s, &p, &joints, 0, kinetic(0.3), bounded, || false)
            .is_err()
    );
    let resting = p.map(|p| state(p, Vec3::ZERO, Vec3::ZERO, Vec3::ZERO));
    let locked = [
        joints[0],
        ArticulatedJoint::Distance(DistanceLink {
            bodies: [0, 1],
            points_body_m: [Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0)],
            length_m: 2.0_f64.sqrt(),
        }),
    ];
    assert!(matches!(
        apply_revolute_friction_impulse(
            &resting,
            &p,
            &locked,
            0,
            RevoluteFrictionInput::RestingDrive {
                drive_impulse_nms: 0.2,
                static_capacity_nms: 0.5,
                kinetic_capacity_nms: 0.3
            },
            controls(),
            || false
        ),
        Err(RevoluteFrictionError::Unresolved(
            "locked or ill-conditioned hinge axis"
        ))
    ));
    assert_eq!(s, initial);
}
