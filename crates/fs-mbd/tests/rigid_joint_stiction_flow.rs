//! Declared breakaway/kinetic capacities on original native hinge and rail ports.
use fs_contact::tangential::impact::multiple::CoupledFrictionImpulseControls;
use fs_contact::tangential::impact::multiple::relaxation::CoupledFrictionFlowControls;
use fs_mbd::rigid_rattle::distance::DistanceStepControls;
use fs_mbd::rigid_rattle::distance::graph::joint_friction::relaxation::trajectory::{
    JointFrictionTrajectoryAccuracy, JointFrictionTrajectoryControls,
    JointFrictionTrajectoryMotion, advance_joint_friction_trajectory,
};
use fs_mbd::rigid_rattle::distance::graph::joint_friction::relaxation::{
    JointFrictionFlowControls, JointFrictionFlowError, JointFrictionFlowMotion,
    JointFrictionFlowPort, JointFrictionFlowReport, MAX_JOINT_STICTION_MODES_PER_PHASE,
    advance_joint_friction,
};
use fs_mbd::rigid_rattle::distance::graph::{
    ArticulatedJoint, DistanceGraphControls, PrismaticJoint, RevoluteJoint,
};
use fs_mbd::{
    Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState, UnitQuaternion, Vec3,
    Wrench,
};
use std::cell::Cell;

const X: Vec3 = Vec3::new(1.0, 0.0, 0.0);
const Z: Vec3 = Vec3::new(0.0, 0.0, 1.0);

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}; tolerance={tolerance:e}"
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
fn hinge(bodies: [usize; 2], length: f64) -> ArticulatedJoint {
    ArticulatedJoint::Revolute(RevoluteJoint::new(bodies, [Vec3::ZERO; 2], [Z; 2], length).unwrap())
}
fn rail(bodies: [usize; 2], points: [Vec3; 2]) -> ArticulatedJoint {
    ArticulatedJoint::Prismatic(
        PrismaticJoint::new(bodies, points, [UnitQuaternion::IDENTITY; 2], 1.0).unwrap(),
    )
}
fn rotor_port(joint_index: usize, cs: f64, ck: f64, drive: f64) -> JointFrictionFlowPort {
    JointFrictionFlowPort::RevoluteStiction {
        joint_index,
        static_capacity_nm: cs,
        kinetic_capacity_nm: ck,
        drive_nm: drive,
    }
}
fn rail_port(cs: f64, ck: f64, drive: f64) -> JointFrictionFlowPort {
    JointFrictionFlowPort::PrismaticStiction {
        joint_index: 0,
        static_capacity_n: cs,
        kinetic_capacity_n: ck,
        drive_n: drive,
    }
}
fn slider(
    relative_velocity: f64,
) -> (
    [RigidBodyState; 2],
    [MassProperties; 2],
    [ArticulatedJoint; 1],
) {
    let p = [properties(2.0, 1.0); 2];
    let s = [
        state(
            p[0],
            Vec3::ZERO,
            Z.scale(relative_velocity / 2.0),
            Vec3::ZERO,
        ),
        state(
            p[1],
            Vec3::ZERO,
            Z.scale(-relative_velocity / 2.0),
            Vec3::ZERO,
        ),
    ];
    (s, p, [rail([0, 1], [Vec3::ZERO; 2])])
}
fn run(
    s: &[RigidBodyState],
    p: &[MassProperties],
    j: &[ArticulatedJoint],
    ports: &[JointFrictionFlowPort],
    duration: f64,
) -> JointFrictionFlowReport {
    advance_joint_friction(s, p, j, ports, duration, controls(), || false).unwrap()
}
fn spin(s: RigidBodyState, p: MassProperties) -> Vec3 {
    s.point_kinematics(p, Vec3::ZERO)
        .unwrap()
        .angular_velocity_world
}
fn audit(before: &[RigidBodyState], p: &[MassProperties], out: &JointFrictionFlowReport) {
    let integrator = RigidBodyIntegrator::new(Gravity::ZERO);
    let mut energy = 0.0;
    let mut linear = Vec3::ZERO;
    let mut angular = Vec3::ZERO;
    for ((&a, &b), &p) in before.iter().zip(&out.states_after).zip(p) {
        assert_eq!(a.pose(), b.pose());
        let a = integrator.diagnostics(a, p).unwrap();
        let b = integrator.diagnostics(b, p).unwrap();
        energy += b.mechanical_energy - a.mechanical_energy;
        linear = linear.add(b.linear_momentum_world.sub(a.linear_momentum_world));
        angular = angular.add(b.angular_momentum_world.sub(a.angular_momentum_world));
    }
    near(energy, out.kinetic_energy_change_j, 3e-10);
    near(
        energy,
        out.drive_work_j - out.dissipated_j + out.constraint_work_j,
        3e-10,
    );
    near(out.work_energy_residual_j, 0.0, 3e-10);
    near(
        out.native_impulse_work_j,
        out.drive_work_j - out.dissipated_j,
        3e-10,
    );
    vector(linear, Vec3::ZERO, 3e-10);
    vector(angular, Vec3::ZERO, 3e-10);
    vector(out.linear_momentum_change_ns, linear, 3e-10);
    vector(out.angular_momentum_change_nms, angular, 3e-10);
    for port in &out.ports {
        assert!(port.friction_work_j <= 0.0);
        assert!(port.minimum_scaled_speed_m_s >= 0.0);
        assert!(port.minimum_scaled_speed_m_s <= port.maximum_scaled_speed_m_s);
        near(port.flow_velocity_residual_m_s, 0.0, 3e-10);
        near(port.scaled_velocity_residual_m_s, 0.0, 3e-10);
        near(port.drive_impulse_residual_ns, 0.0, 3e-10);
    }
    if let Some(d) = out.stiction {
        assert!(d.mode_evaluations > 0);
        assert!(
            d.mode_evaluations
                <= MAX_JOINT_STICTION_MODES_PER_PHASE * controls().flow.maximum_segments
        );
    }
}

#[test]
fn g1_static_balance_and_lower_moving_friction_remain_distinct_at_arbitrarily_small_speed() {
    let source = [rail_port(2.0, 1.0, 1.5)];
    let (s, p, j) = slider(0.0);
    let held = run(&s, &p, &j, &source, 1.0);
    assert_eq!(held.states_after, s);
    near(held.flow.friction_impulse_ns[0], -1.5, 0.0);
    near(held.drive_work_j, 0.0, 0.0);
    near(held.dissipated_j, 0.0, 0.0);
    near(held.ports[0].maximum_scaled_speed_m_s, 0.0, 0.0);
    assert_eq!(held.stiction.unwrap().last_phase_holding_mask, 1);
    audit(&s, &p, &held);

    // Even 1e-12 is a real moving input, not a velocity capture band.
    for v in [1.0, 1e-12] {
        let (s, p, j) = slider(v);
        let out = run(&s, &p, &j, &source, 1.0);
        near(out.flow.flow_velocity_after_m_s[0], v + 0.5, 3e-12);
        near(out.flow.friction_impulse_ns[0], -1.0, 3e-12);
        near(out.flow.displacement_m[0], v + 0.25, 3e-12);
        near(out.dissipated_j, v + 0.25, 3e-12);
        near(out.drive_work_j, 1.5 * (v + 0.25), 3e-12);
        near(out.ports[0].minimum_scaled_speed_m_s, v, 3e-12);
        near(out.ports[0].maximum_scaled_speed_m_s, v + 0.5, 3e-12);
        audit(&s, &p, &out);
    }
}

#[test]
fn g1_hinge_exact_static_threshold_survives_nonbinary_work_conjugate_scales() {
    let p = [properties(1.0, 2.0); 2];
    let s = [state(p[0], Vec3::ZERO, Vec3::ZERO, Vec3::ZERO); 2];
    for length in [0.3, 0.7, 1.0, 2.3] {
        let j = [hinge([0, 1], length)];
        // Down-rounding only Cs/L would spuriously break away from this state.
        let held = run(&s, &p, &j, &[rotor_port(0, 3.7, 1.2, 3.7)], 0.5);
        assert_eq!(held.states_after, s);
        near(held.dissipated_j, 0.0, 0.0);
        near(held.flow.displacement_m[0], 0.0, 0.0);
        audit(&s, &p, &held);

        let released = run(&s, &p, &j, &[rotor_port(0, 3.7, 1.2, 4.0)], 0.5);
        let JointFrictionFlowMotion::Revolute {
            rate_after_rad_s,
            rotation_rad,
            friction_impulse_nms,
            drive_impulse_nms,
            ..
        } = released.ports[0].motion
        else {
            panic!("hinge");
        };
        near(rate_after_rad_s, 1.4, 4e-12);
        near(rotation_rad, 0.35, 4e-12);
        near(friction_impulse_nms, -0.6, 4e-12);
        near(drive_impulse_nms, 2.0, 4e-12);
        near(released.dissipated_j, 0.42, 4e-12);
        near(released.drive_work_j, 1.4, 4e-12);
        near(
            released.ports[0].maximum_scaled_speed_m_s,
            length * 1.4,
            4e-12,
        );
        assert_eq!(released.stiction.unwrap().breakaways, 1);
        audit(&s, &p, &released);
    }
}

#[test]
fn g1_arrest_recaptures_below_static_limit_and_reverses_above_it_with_path_loss() {
    let (s, p, j) = slider(1.0);
    // F=-1.5: arrest at .4 s, then +1.5 N static force holds the rest.
    // F=-3: arrest at .25 s, then lower +1 N kinetic force permits reversal.
    for (drive, end, travel, distance, drive_work) in [
        (-1.5, 0.0, 0.2, 0.2, -0.3),
        (-3.0, -1.5, -0.4375, 0.6875, 1.3125),
    ] {
        let out = run(&s, &p, &j, &[rail_port(2.0, 1.0, drive)], 1.0);
        let JointFrictionFlowMotion::Prismatic {
            velocity_after_m_s,
            displacement_m,
            distance_m,
            friction_impulse_ns,
            ..
        } = out.ports[0].motion
        else {
            panic!("rail");
        };
        near(velocity_after_m_s, end, 3e-12);
        near(displacement_m, travel, 3e-12);
        near(distance_m, distance, 3e-12);
        near(friction_impulse_ns, 0.5, 3e-12);
        near(out.dissipated_j, distance, 3e-12);
        near(out.drive_work_j, drive_work, 3e-12);
        near(out.ports[0].minimum_scaled_speed_m_s, 0.0, 0.0);
        near(
            out.ports[0].maximum_scaled_speed_m_s,
            1.0_f64.max(end.abs()),
            3e-12,
        );
        audit(&s, &p, &out);
    }
}

#[test]
fn g1_coupled_rotors_hold_one_joint_above_its_kinetic_bound_then_arrest_together() {
    let p = [
        properties(1.0, 2.0),
        properties(2.0, 3.0),
        properties(3.0, 5.0),
    ];
    let s = core::array::from_fn::<_, 3, _>(|i| {
        state(p[i], Vec3::ZERO, Vec3::ZERO, Z.scale([0.0, 1.0, 1.0][i]))
    });
    let ports = [rotor_port(0, 2.0, 1.0, 0.0), rotor_port(1, 1.0, 0.1, 0.0)];
    // W=[[5/6,-1/3],[-1/3,8/15]]. Joint 1 holds with -5/8 N m,
    // beyond its .1 kinetic bound; joint 0 decelerates at -5/8 rad/s2.
    for lengths in [[1.0, 1.0], [0.25, 2.0], [3.0, 0.3]] {
        let j = [hinge([0, 1], lengths[0]), hinge([1, 2], lengths[1])];
        let out = run(&s, &p, &j, &ports, 2.0);
        assert_eq!(out, run(&s, &p, &j, &[ports[1], ports[0]], 2.0));
        for (i, &p) in p.iter().enumerate() {
            vector(spin(out.states_after[i], p), Z.scale(0.8), 4e-12);
        }
        for (i, impulse) in [-1.6, -1.0].into_iter().enumerate() {
            let JointFrictionFlowMotion::Revolute {
                rotation_rad,
                angular_distance_rad,
                friction_impulse_nms,
                ..
            } = out.ports[i].motion
            else {
                panic!("hinge");
            };
            near(friction_impulse_nms, impulse, 4e-12);
            near(rotation_rad, if i == 0 { 0.8 } else { 0.0 }, 4e-12);
            near(angular_distance_rad, rotation_rad, 4e-12);
        }
        near(out.ports[1].friction_work_j, 0.0, 4e-12);
        near(out.dissipated_j, 0.8, 4e-12);
        audit(&s, &p, &out);
    }
}

#[test]
fn g3_eccentric_full_tensor_recoil_and_arrest_are_covariant_under_frame_and_boost() {
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
    let j = [rail([0, 1], [X.scale(-2.0), Vec3::ZERO])];
    let ports = [rail_port(2.0, 1.0, 0.0)];
    let base = run(&s, &p, &j, &ports, 4.0);
    let mobility = 215.0 / 114.0;
    let impulse = -3.0 / mobility;
    let omega = Vec3::new(-2.0 / 19.0, 10.0 / 19.0, 0.0).scale(impulse);
    near(base.scaled_mobility_per_kg[0][0], mobility, 3e-12);
    near(base.flow.friction_impulse_ns[0], impulse, 3e-12);
    near(base.dissipated_j, 4.5 / mobility, 3e-12);
    for (i, &p) in p.iter().enumerate() {
        vector(spin(base.states_after[i], p), omega, 3e-12);
    }
    audit(&s, &p, &base);

    let q = UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, 3.0), 0.61).unwrap();
    let offset = Vec3::new(5.0, -7.0, 11.0);
    let boost = Vec3::new(3.0, -4.0, 2.0);
    let transform = |s: RigidBodyState, p: MassProperties| {
        RigidBodyState::new(
            Pose::new(
                q.rotate_body_to_world(s.pose().position_world())
                    .add(offset),
                q,
            )
            .unwrap(),
            q.rotate_body_to_world(s.linear_momentum_world())
                .add(boost.scale(p.mass())),
            s.angular_momentum_body(),
        )
        .unwrap()
    };
    let transformed = core::array::from_fn::<_, 2, _>(|i| transform(s[i], p[i]));
    let out = run(&transformed, &p, &j, &ports, 4.0);
    for (i, &p) in p.iter().enumerate() {
        let expected = transform(base.states_after[i], p);
        vector(
            out.states_after[i].linear_momentum_world(),
            expected.linear_momentum_world(),
            4e-11,
        );
        vector(
            out.states_after[i].angular_momentum_body(),
            expected.angular_momentum_body(),
            4e-11,
        );
    }
    near(out.dissipated_j, base.dissipated_j, 4e-11);
    near(out.flow.distance_m[0], base.flow.distance_m[0], 4e-11);
    vector(out.ports[0].axis_world, q.rotate_body_to_world(Z), 3e-13);
    audit(&transformed, &p, &out);
}

#[test]
fn g3_equal_static_and_kinetic_variants_reuse_the_original_flow_exactly() {
    let (s, p, j) = slider(1.0);
    let original = JointFrictionFlowPort::Prismatic {
        joint_index: 0,
        capacity_n: 1.0,
        drive_n: -2.0,
    };
    let old = run(&s, &p, &j, &[original], 1.0);
    assert_eq!(old, run(&s, &p, &j, &[rail_port(1.0, 1.0, -2.0)], 1.0));
    assert!(old.stiction.is_none());
    let p = [properties(1.0, 2.0); 2];
    let s = [
        state(p[0], Vec3::ZERO, Vec3::ZERO, Z.scale(-0.5)),
        state(p[1], Vec3::ZERO, Vec3::ZERO, Z.scale(0.5)),
    ];
    for length in [0.3, 0.7, 2.3] {
        let j = [hinge([0, 1], length)];
        let original = JointFrictionFlowPort::Revolute {
            joint_index: 0,
            capacity_nm: 1.0,
            drive_nm: 0.0,
        };
        assert_eq!(
            run(&s, &p, &j, &[original], 2.0),
            run(&s, &p, &j, &[rotor_port(0, 1.0, 1.0, 0.0)], 2.0),
        );
    }
}

#[test]
fn g2_existing_adaptive_motion_uses_static_source_breakaway_and_kinetic_arrest() {
    let (s, p, j) = slider(0.0);
    let initial = [rail_port(2.0, 1.0, 1.5)];
    let c = JointFrictionTrajectoryControls {
        mechanical: controls(),
        maximum_step_s: 1.0,
        minimum_step_s: 1e-8,
        maximum_accepted_steps: 2048,
        maximum_split_trials: 8192,
        maximum_live_bytes: 5_000_000,
        accuracy: JointFrictionTrajectoryAccuracy {
            position_m: 1e-6,
            orientation_rad: 1e-6,
            linear_momentum_ns: 1e-6,
            angular_momentum_nms: 1e-6,
            scaled_travel_m: 1e-6,
            scaled_impulse_ns: 1e-6,
            work_j: 1e-6,
        },
        energy_residual_tolerance_j: 1e-7,
        constraint_work_tolerance_j: 1e-5,
    };
    let held_samples = Cell::new(0);
    let moving_samples = Cell::new(0);
    let out = advance_joint_friction_trajectory(
        &s,
        &p,
        &j,
        &initial,
        RigidBodyIntegrator::new(Gravity::ZERO),
        &[Wrench::ZERO; 2],
        0.0,
        4.0,
        &[1.0, 2.0],
        c,
        |time, states, ports| {
            assert!(matches!(
                ports[0],
                JointFrictionFlowPort::PrismaticStiction { .. }
            ));
            let drive = if time < 1.0 {
                held_samples.set(held_samples.get() + 1);
                assert_eq!(states, s);
                1.5
            } else if time < 2.0 {
                moving_samples.set(moving_samples.get() + 1);
                3.0
            } else {
                0.0
            };
            ports[0] = rail_port(2.0, 1.0, drive);
            Ok::<_, ()>(())
        },
        || false,
    )
    .unwrap();
    assert!(held_samples.get() >= 3 && moving_samples.get() >= 3);
    assert_eq!(initial, [rail_port(2.0, 1.0, 1.5)]);
    // Hold for one second; accelerate at 2 for one; decelerate at 1 for two.
    let JointFrictionTrajectoryMotion::Prismatic {
        displacement_m,
        distance_m,
        friction_impulse_ns,
        drive_impulse_ns,
    } = out.ports[0].motion
    else {
        panic!("rail");
    };
    near(displacement_m, 3.0, 3e-11);
    near(distance_m, 3.0, 3e-11);
    near(friction_impulse_ns, -4.5, 3e-11);
    near(drive_impulse_ns, 4.5, 3e-11);
    near(out.drive_work_j, 3.0, 3e-11);
    near(out.dissipated_j, 3.0, 3e-11);
    near(out.mechanical_energy_change_j, 0.0, 3e-11);
    near(out.work_energy_residual_j, 0.0, 3e-11);
    for (i, sign) in [1.0, -1.0].into_iter().enumerate() {
        vector(
            out.states_after[i].linear_momentum_world(),
            Vec3::ZERO,
            3e-11,
        );
        near(
            out.states_after[i].pose().position_world().z,
            sign * 1.5,
            4e-7,
        );
    }

    // The existing adaptive owner must still recognize the common segment
    // refusal through the stiction dispatch, and refine without raising it.
    let (s, p, j) = slider(1.0);
    let mut one_segment = c;
    one_segment.mechanical.flow.maximum_segments = 1;
    let refined = advance_joint_friction_trajectory(
        &s,
        &p,
        &j,
        &[rail_port(3.0, 2.0, 0.0)],
        RigidBodyIntegrator::new(Gravity::ZERO),
        &[Wrench::ZERO; 2],
        0.0,
        1.0,
        &[],
        one_segment,
        |_, _, _| Ok::<_, ()>(()),
        || false,
    )
    .unwrap();
    assert!(refined.rejected_intervals > 0);
    near(refined.dissipated_j, 0.5, 3e-11);
    for (i, sign) in [1.0, -1.0].into_iter().enumerate() {
        vector(
            refined.states_after[i].linear_momentum_world(),
            Vec3::ZERO,
            3e-11,
        );
        near(
            refined.states_after[i].pose().position_world().z,
            sign * 0.125,
            4e-7,
        );
    }
}

#[test]
fn g4_bad_static_bounds_force_segment_memory_limits_and_cancellation_are_atomic() {
    let (s, p, j) = slider(1.0);
    let original = s;
    for source in [
        rail_port(0.5, 1.0, 0.0),
        rail_port(-1.0, 0.0, 0.0),
        rail_port(2.0, f64::NAN, 0.0),
        rail_port(f64::INFINITY, 1.0, 0.0),
    ] {
        assert!(matches!(
            advance_joint_friction(&s, &p, &j, &[source], 1.0, controls(), || false),
            Err(JointFrictionFlowError::Invalid(_)),
        ));
    }
    let ports = [rail_port(2.0, 1.0, -3.0)];
    let mut too_few_segments = controls();
    too_few_segments.flow.maximum_segments = 1;
    let mut too_little_force = controls();
    too_little_force.flow.maximum_force_n = 0.5;
    let mut too_little_memory = controls();
    too_little_memory.graph.maximum_live_bytes = 1;
    for c in [too_few_segments, too_little_force, too_little_memory] {
        assert!(advance_joint_friction(&s, &p, &j, &ports, 1.0, c, || false).is_err());
    }
    let mut total = 0;
    let _ = advance_joint_friction(&s, &p, &j, &ports, 1.0, controls(), || {
        total += 1;
        false
    })
    .unwrap();
    assert!(total > 3);
    for stop in [1, total / 2, total] {
        let mut calls = 0;
        assert!(
            advance_joint_friction(&s, &p, &j, &ports, 1.0, controls(), || {
                calls += 1;
                calls >= stop
            })
            .is_err()
        );
        assert_eq!(calls, stop);
        assert_eq!(s, original);
    }
}
