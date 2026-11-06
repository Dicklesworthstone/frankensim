use fs_contact::tangential::impact::multiple::CoupledFrictionImpulseControls;
use fs_contact::tangential::impact::multiple::relaxation::CoupledFrictionFlowControls;
use fs_mbd::rigid_rattle::distance::DistanceStepControls;
use fs_mbd::rigid_rattle::distance::graph::joint_friction::relaxation::{
    JointFrictionFlowControls, JointFrictionFlowMotion, JointFrictionFlowPort,
    JointFrictionFlowReport, advance_joint_friction,
};
use fs_mbd::rigid_rattle::distance::graph::{
    ArticulatedJoint, DistanceGraphControls, PrismaticJoint, RevoluteJoint,
};
use fs_mbd::{
    Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState, UnitQuaternion, Vec3,
};

const X: Vec3 = Vec3::new(1.0, 0.0, 0.0);
const Y: Vec3 = Vec3::new(0.0, 1.0, 0.0);
const Z: Vec3 = Vec3::new(0.0, 0.0, 1.0);

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
fn properties(mass: f64, inertia: f64) -> MassProperties {
    MassProperties::new(mass, Vec3::ZERO, Vec3::new(inertia, inertia, inertia)).unwrap()
}
fn state(p: MassProperties, position: Vec3, velocity: Vec3, spin: Vec3) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(position, UnitQuaternion::IDENTITY).unwrap(),
        velocity.scale(p.mass()),
        p.angular_momentum_body_checked(spin).unwrap(),
    )
    .unwrap()
}
fn controls() -> JointFrictionFlowControls {
    JointFrictionFlowControls {
        graph: DistanceGraphControls {
            step: DistanceStepControls {
                maximum_impulse_ns: 1000.0,
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
        flow: CoupledFrictionFlowControls {
            projection: CoupledFrictionImpulseControls {
                maximum_iterations: 256,
                impulse_tolerance_ns: 1e-10,
                velocity_tolerance_m_s: 1e-10,
                energy_tolerance_j: 1e-9,
                relative_matrix_tolerance: 1e-12,
                maximum_impulse_ns: 1000.0,
            },
            maximum_segments: 32,
            acceleration_tolerance_m_s2: 1e-10,
            maximum_force_n: 1000.0,
        },
        work_tolerance_j: 1e-9,
        linear_impulse_tolerance_ns: 1e-9,
        angular_impulse_tolerance_nms: 1e-9,
    }
}
fn hinge(bodies: [usize; 2], axis: Vec3, length: f64) -> ArticulatedJoint {
    ArticulatedJoint::Revolute(
        RevoluteJoint::new(bodies, [Vec3::ZERO; 2], [axis; 2], length).unwrap(),
    )
}
fn slider(bodies: [usize; 2], points: [Vec3; 2]) -> ArticulatedJoint {
    ArticulatedJoint::Prismatic(
        PrismaticJoint::new(bodies, points, [UnitQuaternion::IDENTITY; 2], 1.0).unwrap(),
    )
}
fn rotor_port(joint_index: usize, capacity_nm: f64, drive_nm: f64) -> JointFrictionFlowPort {
    JointFrictionFlowPort::Revolute {
        joint_index,
        capacity_nm,
        drive_nm,
    }
}
fn rail_port(joint_index: usize, capacity_n: f64, drive_n: f64) -> JointFrictionFlowPort {
    JointFrictionFlowPort::Prismatic {
        joint_index,
        capacity_n,
        drive_n,
    }
}
fn run(
    s: &[RigidBodyState],
    p: &[MassProperties],
    j: &[ArticulatedJoint],
    ports: &[JointFrictionFlowPort],
    t: f64,
) -> JointFrictionFlowReport {
    advance_joint_friction(s, p, j, ports, t, controls(), || false).unwrap()
}
fn spin(s: RigidBodyState, p: MassProperties) -> Vec3 {
    s.point_kinematics(p, Vec3::ZERO)
        .unwrap()
        .angular_velocity_world
}
fn velocity(s: RigidBodyState, p: MassProperties) -> Vec3 {
    s.center_of_mass_velocity_world(p).unwrap()
}
fn audit(before: &[RigidBodyState], properties: &[MassProperties], out: &JointFrictionFlowReport) {
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
    near(energy, out.kinetic_energy_change_j, 2e-10);
    near(
        energy,
        out.drive_work_j - out.dissipated_j + out.constraint_work_j,
        2e-10,
    );
    near(out.work_energy_residual_j, 0.0, 2e-10);
    near(out.constraint_work_j, 0.0, 2e-10);
    near(
        out.native_impulse_work_j,
        out.drive_work_j - out.dissipated_j,
        2e-10,
    );
    assert!(out.dissipated_j >= 0.0);
    vector(linear, Vec3::ZERO, 2e-10);
    vector(angular, Vec3::ZERO, 2e-10);
    vector(out.linear_momentum_change_ns, linear, 2e-10);
    vector(out.angular_momentum_change_nms, angular, 2e-10);
    for port in &out.ports {
        assert!(port.friction_work_j <= 0.0);
        near(port.scaled_velocity_residual_m_s, 0.0, 2e-10);
        near(port.flow_velocity_residual_m_s, 0.0, 2e-10);
        near(port.drive_impulse_residual_ns, 0.0, 2e-10);
    }
}

fn rotors(
    lengths: [f64; 2],
) -> (
    [RigidBodyState; 3],
    [MassProperties; 3],
    [ArticulatedJoint; 2],
) {
    let p = [
        properties(1.0, 2.0),
        properties(2.0, 3.0),
        properties(3.0, 5.0),
    ];
    let s =
        core::array::from_fn(|i| state(p[i], Vec3::ZERO, Vec3::ZERO, Z.scale([3.0, -1.0, 2.0][i])));
    (
        s,
        p,
        [hinge([0, 1], Z, lengths[0]), hinge([1, 2], Z, lengths[1])],
    )
}

#[test]
fn g1_three_rotors_cross_two_physical_arrest_events_and_hold() {
    // I=(2,3,5), rates=(-4,3), W=[[5/6,-1/3],[-1/3,8/15]].
    // First f=(1,-2), a=(3/2,-7/5): joint 1 stops at 15/7 s.
    // Then f=(1,5/8), a=(5/8,0): joint 0 stops at 17/5 s.
    // Both hold for the rest of the five-second interval.
    let ports = [rotor_port(0, 1.0, 0.0), rotor_port(1, 2.0, 0.0)];
    for lengths in [[1.0, 1.0], [0.25, 2.0], [3.0, 0.3]] {
        let (s, p, j) = rotors(lengths);
        let out = run(&s, &p, &j, &ports, 5.0);
        for i in 0..3 {
            vector(spin(out.states_after[i], p[i]), Z.scale(1.3), 3e-12);
        }
        for (i, (impulse, distance)) in [(17.0 / 5.0, 787.0 / 140.0), (-7.0 / 2.0, 45.0 / 14.0)]
            .into_iter()
            .enumerate()
        {
            let JointFrictionFlowMotion::Revolute {
                rate_after_rad_s,
                rotation_rad,
                angular_distance_rad,
                friction_impulse_nms,
                ..
            } = out.ports[i].motion
            else {
                panic!("hinge");
            };
            near(rate_after_rad_s, 0.0, 3e-12);
            near(friction_impulse_nms, impulse, 3e-12);
            near(angular_distance_rad, distance, 3e-12);
            near(
                rotation_rad,
                if i == 0 { -distance } else { distance },
                3e-12,
            );
        }
        near(out.dissipated_j, 241.0 / 20.0, 5e-12);
        assert_eq!(out, run(&s, &p, &j, &[ports[1], ports[0]], 5.0));
        audit(&s, &p, &out);
    }
}

#[test]
fn g1_split_held_load_phases_preserve_native_recoil_and_path_work() {
    let (s, p, j) = rotors([1.0, 1.0]);
    let ports = [rotor_port(0, 1.0, 0.0), rotor_port(1, 2.0, 0.0)];
    let whole = run(&s, &p, &j, &ports, 5.0);
    let first = run(&s, &p, &j, &ports, 1.0);
    let second = run(&first.states_after, &p, &j, &ports, 4.0);
    for i in 0..3 {
        vector(
            second.states_after[i].angular_momentum_body(),
            whole.states_after[i].angular_momentum_body(),
            3e-12,
        );
    }
    near(
        first.dissipated_j + second.dissipated_j,
        whole.dissipated_j,
        5e-12,
    );
    for i in 0..2 {
        near(
            first.flow.displacement_m[i] + second.flow.displacement_m[i],
            whole.flow.displacement_m[i],
            5e-12,
        );
        near(
            first.flow.distance_m[i] + second.flow.distance_m[i],
            whole.flow.distance_m[i],
            5e-12,
        );
    }
    audit(&first.states_after, &p, &second);
}

#[test]
fn g1_reversing_slider_integrates_path_loss_instead_of_endpoint_impulse_work() {
    // Relative mobility 1 kg^-1, initial speed 1, drive -2 N, capacity 1 N.
    // Acceleration -3 until t=1/3; thereafter -1. The endpoint is -2/3.
    // Friction impulse is +1/3 Ns, but physical friction work is -7/18 J.
    // Its global midpoint pairing would be +1/18 J and is not the path work.
    let p = [properties(2.0, 2.0); 2];
    let s = [
        state(p[0], Vec3::ZERO, Z.scale(0.5), Vec3::ZERO),
        state(p[1], Vec3::ZERO, Z.scale(-0.5), Vec3::ZERO),
    ];
    let j = [slider([0, 1], [Vec3::ZERO; 2])];
    let out = run(&s, &p, &j, &[rail_port(0, 1.0, -2.0)], 1.0);
    let JointFrictionFlowMotion::Prismatic {
        velocity_after_m_s,
        displacement_m,
        distance_m,
        friction_impulse_ns,
        drive_impulse_ns,
        ..
    } = out.ports[0].motion
    else {
        panic!("rail");
    };
    near(velocity_after_m_s, -2.0 / 3.0, 3e-12);
    near(displacement_m, -1.0 / 18.0, 3e-12);
    near(distance_m, 7.0 / 18.0, 3e-12);
    near(friction_impulse_ns, 1.0 / 3.0, 3e-12);
    near(drive_impulse_ns, -2.0, 3e-12);
    near(out.dissipated_j, 7.0 / 18.0, 3e-12);
    near(out.drive_work_j, 1.0 / 9.0, 3e-12);
    near(out.kinetic_energy_change_j, -5.0 / 18.0, 3e-12);
    vector(
        velocity(out.states_after[0], p[0]),
        Z.scale(-1.0 / 3.0),
        3e-12,
    );
    vector(
        velocity(out.states_after[1], p[1]),
        Z.scale(1.0 / 3.0),
        3e-12,
    );
    audit(&s, &p, &out);
}

#[test]
fn g1_changed_drives_hold_launch_and_arrest_the_same_native_slider() {
    let p = [properties(2.0, 2.0); 2];
    let s = [state(p[0], Vec3::ZERO, Vec3::ZERO, Vec3::ZERO); 2];
    let j = [slider([0, 1], [Vec3::ZERO; 2])];
    let held = run(&s, &p, &j, &[rail_port(0, 2.0, 1.0)], 3.0);
    assert_eq!(held.states_after, s);
    near(held.flow.drive_impulse_ns[0], 3.0, 0.0);
    near(held.flow.friction_impulse_ns[0], -3.0, 0.0);
    near(held.dissipated_j, 0.0, 0.0);
    near(held.drive_work_j, 0.0, 0.0);
    let launched = run(&held.states_after, &p, &j, &[rail_port(0, 2.0, 3.0)], 2.0);
    near(launched.flow.velocity_after_m_s[0], 2.0, 3e-12);
    near(launched.flow.displacement_m[0], 2.0, 3e-12);
    near(launched.dissipated_j, 4.0, 3e-12);
    near(launched.drive_work_j, 6.0, 3e-12);
    let stopped = run(
        &launched.states_after,
        &p,
        &j,
        &[rail_port(0, 2.0, 0.0)],
        2.0,
    );
    near(stopped.flow.velocity_after_m_s[0], 0.0, 3e-12);
    near(stopped.flow.displacement_m[0], 1.0, 3e-12);
    near(stopped.dissipated_j, 2.0, 3e-12);
    audit(&s, &p, &held);
    audit(&held.states_after, &p, &launched);
    audit(&launched.states_after, &p, &stopped);
}

fn mixed(
    length: f64,
) -> (
    [RigidBodyState; 3],
    [MassProperties; 3],
    [ArticulatedJoint; 2],
) {
    let p = [
        properties(1.0, 2.0),
        properties(1.0, 3.0),
        properties(2.0, 2.0),
    ];
    let s = [
        state(p[0], Vec3::ZERO, Vec3::ZERO, Y.scale(3.0)),
        state(p[1], Vec3::ZERO, Vec3::ZERO, Y.scale(-1.0)),
        state(
            p[2],
            Vec3::new(2.0, 0.0, 1.0),
            Vec3::new(-1.0, 0.0, 5.0),
            Y.scale(-1.0),
        ),
    ];
    (
        s,
        p,
        [
            hinge([0, 1], Y, length),
            slider([2, 1], [X.scale(-2.0), Vec3::ZERO]),
        ],
    )
}

#[test]
fn g1_mixed_hinge_and_rotating_rail_follow_the_independent_reduced_trajectory() {
    // Reduced physical W=[[2/3,1/3],[1/3,5/3]], v=(-4,3), C=(1,2).
    // The rail stops at one second; its holding force becomes -1/5 N.
    // Over two seconds p=(2,-11/5), travel=(-77/10,3/2), loss=107/10.
    for length in [1.0, 0.2] {
        let (s, p, j) = mixed(length);
        let out = run(
            &s,
            &p,
            &j,
            &[rotor_port(0, 1.0, 0.0), rail_port(1, 2.0, 0.0)],
            2.0,
        );
        let JointFrictionFlowMotion::Revolute {
            rate_after_rad_s,
            rotation_rad,
            friction_impulse_nms,
            ..
        } = out.ports[0].motion
        else {
            panic!("hinge");
        };
        near(rate_after_rad_s, -17.0 / 5.0, 3e-12);
        near(rotation_rad, -77.0 / 10.0, 3e-12);
        near(friction_impulse_nms, 2.0, 3e-12);
        let JointFrictionFlowMotion::Prismatic {
            velocity_after_m_s,
            displacement_m,
            friction_impulse_ns,
            ..
        } = out.ports[1].motion
        else {
            panic!("rail");
        };
        near(velocity_after_m_s, 0.0, 3e-12);
        near(displacement_m, 1.5, 3e-12);
        near(friction_impulse_ns, -11.0 / 5.0, 3e-12);
        near(out.dissipated_j, 107.0 / 10.0, 5e-12);
        for i in 0..3 {
            vector(
                spin(out.states_after[i], p[i]),
                Y.scale([2.0, -1.4, -1.4][i]),
                3e-12,
            );
            vector(
                velocity(out.states_after[i], p[i]),
                if i == 2 {
                    Vec3::new(-1.2, 0.0, 3.9)
                } else {
                    Vec3::new(0.2, 0.0, 1.1)
                },
                3e-12,
            );
        }
        audit(&s, &p, &out);
    }
}

#[test]
fn g1_full_tensor_eccentric_rail_preserves_rotational_recoil_through_arrest() {
    let p = [0.4, 0.6].map(|share| {
        MassProperties::from_inertia_tensor(
            5.0 * share,
            Vec3::ZERO,
            fs_ga::Mat3 {
                m: [5.0, 1.0, 0.0, 1.0, 4.0, 0.0, 0.0, 0.0, 6.0].map(|x| x * share),
            },
        )
        .unwrap()
    });
    let s = [
        state(p[0], X.scale(2.0), Z.scale(9.0 / 5.0), Vec3::ZERO),
        state(p[1], Vec3::ZERO, Z.scale(-6.0 / 5.0), Vec3::ZERO),
    ];
    let j = [slider([0, 1], [X.scale(-2.0), Vec3::ZERO])];
    let out = run(&s, &p, &j, &[rail_port(0, 1.0, 0.0)], 4.0);
    // W=1/mu+4*(I_sum^-1)_yy=5/6+20/19. Integrate a=-W to arrest.
    let mobility = 215.0 / 114.0;
    let impulse = -3.0 / mobility;
    let omega = Vec3::new(-2.0 / 19.0, 10.0 / 19.0, 0.0).scale(impulse);
    near(out.scaled_mobility_per_kg[0][0], mobility, 3e-12);
    near(out.flow.friction_impulse_ns[0], impulse, 3e-12);
    near(out.flow.displacement_m[0], 4.5 / mobility, 3e-12);
    near(out.dissipated_j, 4.5 / mobility, 3e-12);
    for i in 0..2 {
        vector(spin(out.states_after[i], p[i]), omega, 3e-12);
    }
    assert!(omega.x.abs() > 0.1);
    audit(&s, &p, &out);
}

#[test]
fn g3_native_flow_is_covariant_under_frame_origin_and_common_boost() {
    let (s, p, j) = mixed(0.7);
    let ports = [rotor_port(0, 1.0, 0.0), rail_port(1, 2.0, 0.0)];
    let base = run(&s, &p, &j, &ports, 2.0);
    let rotation = UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, 3.0), 0.61).unwrap();
    let offset = Vec3::new(5.0, -7.0, 11.0);
    let boost = Vec3::new(3.0, -4.0, 2.0);
    let transform = |s: RigidBodyState, p: MassProperties| {
        RigidBodyState::new(
            Pose::new(
                rotation
                    .rotate_body_to_world(s.pose().position_world())
                    .add(offset),
                rotation,
            )
            .unwrap(),
            rotation
                .rotate_body_to_world(s.linear_momentum_world())
                .add(boost.scale(p.mass())),
            s.angular_momentum_body(),
        )
        .unwrap()
    };
    let transformed = core::array::from_fn::<_, 3, _>(|i| transform(s[i], p[i]));
    let out = run(&transformed, &p, &j, &ports, 2.0);
    for i in 0..3 {
        let expected = transform(base.states_after[i], p[i]);
        assert_eq!(out.states_after[i].pose(), expected.pose());
        vector(
            out.states_after[i].linear_momentum_world(),
            expected.linear_momentum_world(),
            3e-11,
        );
        vector(
            out.states_after[i].angular_momentum_body(),
            expected.angular_momentum_body(),
            3e-11,
        );
    }
    near(out.dissipated_j, base.dissipated_j, 3e-11);
    for i in 0..2 {
        near(
            out.flow.displacement_m[i],
            base.flow.displacement_m[i],
            3e-11,
        );
        near(out.flow.distance_m[i], base.flow.distance_m[i], 3e-11);
        vector(
            out.ports[i].axis_world,
            rotation.rotate_body_to_world(base.ports[i].axis_world),
            3e-13,
        );
    }
    audit(&transformed, &p, &out);
}

#[test]
fn g4_flow_cancellation_limits_and_empty_identity_preserve_all_inputs() {
    let (s, p, j) = rotors([1.0, 1.0]);
    let saved = s;
    let ports = [rotor_port(0, 1.0, 0.0), rotor_port(1, 2.0, 0.0)];
    let mut polls = 0;
    advance_joint_friction(&s, &p, &j, &ports, 5.0, controls(), || {
        polls += 1;
        false
    })
    .unwrap();
    for stop in 1..=polls {
        let mut current = 0;
        let result = advance_joint_friction(&s, &p, &j, &ports, 5.0, controls(), || {
            current += 1;
            current == stop
        });
        assert!(result.is_err(), "cancellation poll {stop}/{polls}");
        assert_eq!(s, saved);
    }
    let mut limited = controls();
    limited.flow.maximum_segments = 1;
    assert!(advance_joint_friction(&s, &p, &j, &ports, 5.0, limited, || false).is_err());
    limited = controls();
    limited.flow.maximum_force_n = 0.5;
    assert!(advance_joint_friction(&s, &p, &j, &ports, 5.0, limited, || false).is_err());
    limited = controls();
    limited.graph.maximum_live_bytes = 1;
    assert!(advance_joint_friction(&s, &p, &j, &ports, 5.0, limited, || false).is_err());
    assert!(advance_joint_friction(&s, &p, &j, &ports, 0.0, controls(), || false).is_err());
    assert!(
        advance_joint_friction(&s, &p, &j, &[ports[0], ports[0]], 1.0, controls(), || false)
            .is_err()
    );
    let empty = run(&s, &p, &j, &[], 5.0);
    assert_eq!(empty.states_after, s);
    assert_eq!(empty.flow.count, 0);
    assert_eq!(empty.flow.factorizations, 0);
    let zero = run(
        &s,
        &p,
        &j,
        &[rotor_port(0, 0.0, 0.0), rotor_port(1, 0.0, 0.0)],
        2.0,
    );
    assert_eq!(zero.states_after, s);
    near(zero.flow.displacement_m[0], -8.0, 0.0);
    near(zero.flow.displacement_m[1], 6.0, 0.0);
    near(zero.dissipated_j, 0.0, 0.0);
    audit(&s, &p, &empty);
    audit(&s, &p, &zero);
}
