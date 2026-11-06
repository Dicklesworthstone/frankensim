//! Actual native poses, finite carrier recoil and source-resolved dry friction.
use fs_contact::tangential::impact::multiple::CoupledFrictionImpulseControls;
use fs_contact::tangential::impact::multiple::relaxation::CoupledFrictionFlowControls;
use fs_mbd::rigid_rattle::distance::DistanceStepControls;
use fs_mbd::rigid_rattle::distance::graph::joint_friction::relaxation::trajectory::{
    JointFrictionTrajectoryAccuracy, JointFrictionTrajectoryControls, JointFrictionTrajectoryError,
    JointFrictionTrajectoryMotion, JointFrictionTrajectoryReport,
    advance_joint_friction_trajectory,
};
use fs_mbd::rigid_rattle::distance::graph::joint_friction::relaxation::{
    JointFrictionFlowControls, JointFrictionFlowPort,
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
const Y: Vec3 = Vec3::new(0.0, 1.0, 0.0);
const Z: Vec3 = Vec3::new(0.0, 0.0, 1.0);
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}, tolerance={tolerance:e}"
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
fn state(p: MassProperties, x: Vec3, q: UnitQuaternion, v: Vec3, omega: Vec3) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(x, q).unwrap(),
        v.scale(p.mass()),
        p.angular_momentum_body_checked(q.rotate_world_to_body(omega))
            .unwrap(),
    )
    .unwrap()
}
fn free() -> RigidBodyIntegrator {
    RigidBodyIntegrator::new(Gravity::ZERO)
}
fn rail(capacity: f64, drive: f64) -> JointFrictionFlowPort {
    JointFrictionFlowPort::Prismatic {
        joint_index: 0,
        capacity_n: capacity,
        drive_n: drive,
    }
}
fn slider(
    v: f64,
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
            UnitQuaternion::IDENTITY,
            Z.scale(v),
            Vec3::ZERO,
        ),
        state(
            p[1],
            Vec3::ZERO,
            UnitQuaternion::IDENTITY,
            Vec3::ZERO,
            Vec3::ZERO,
        ),
    ];
    let j = [ArticulatedJoint::Prismatic(
        PrismaticJoint::new([0, 1], [Vec3::ZERO; 2], [UnitQuaternion::IDENTITY; 2], 1.0).unwrap(),
    )];
    (s, p, j)
}
fn controls() -> JointFrictionTrajectoryControls {
    JointFrictionTrajectoryControls {
        mechanical: JointFrictionFlowControls {
            graph: DistanceGraphControls {
                step: DistanceStepControls {
                    maximum_impulse_ns: 1000.0,
                    maximum_rotation_increment_rad: 0.25,
                    maximum_position_iterations: 40,
                    maximum_backtracks: 20,
                },
                position_tolerance_m: 1e-10,
                velocity_tolerance_mps: 1e-10,
                minimum_scaled_pivot: 1e-10,
                maximum_bodies: 32,
                maximum_constraints: 32,
                maximum_trial_evaluations: 200,
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
        },
        maximum_step_s: 2.0,
        minimum_step_s: 1e-10,
        maximum_accepted_steps: 8192,
        maximum_split_trials: 32768,
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
    }
}
fn held(
    s: &[RigidBodyState],
    p: &[MassProperties],
    j: &[ArticulatedJoint],
    ports: &[JointFrictionFlowPort],
    duration: f64,
    c: JointFrictionTrajectoryControls,
) -> JointFrictionTrajectoryReport {
    advance_joint_friction_trajectory(
        s,
        p,
        j,
        ports,
        free(),
        &vec![Wrench::ZERO; s.len()],
        0.0,
        duration,
        &[],
        c,
        |_, _, _| Ok::<_, ()>(()),
        || false,
    )
    .unwrap()
}
fn audit(
    initial: &[RigidBodyState],
    p: &[MassProperties],
    integrator: RigidBodyIntegrator,
    out: &JointFrictionTrajectoryReport,
) {
    let change = initial
        .iter()
        .zip(&out.states_after)
        .zip(p)
        .map(|((&a, &b), &p)| {
            integrator.diagnostics(b, p).unwrap().mechanical_energy
                - integrator.diagnostics(a, p).unwrap().mechanical_energy
        })
        .sum::<f64>();
    near(out.mechanical_energy_change_j, change, 3e-12);
    near(
        change,
        out.wrench_work_j + out.drive_work_j - out.dissipated_j + out.constraint_work_j,
        1e-7,
    );
    near(out.work_energy_residual_j, 0.0, 1e-7);
    assert!(out.dissipated_j >= 0.0);
    for p in &out.ports {
        assert!(p.friction_work_j <= 0.0);
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
fn inverse(q: UnitQuaternion) -> UnitQuaternion {
    let [w, x, y, z] = q.components();
    UnitQuaternion::new(w, -x, -y, -z).unwrap()
}

#[test]
fn g1_arrest_advances_the_real_sliding_distance_and_finite_carrier_recoil() {
    // Reduced mass is 1 kg. Relative speed 1 arrests at t=1 after 1/2 m.
    // One unrefined two-second symmetric split would wrongly travel 1 m.
    let (s, p, j) = slider(1.0);
    let out = held(&s, &p, &j, &[rail(1.0, 0.0)], 2.0, controls());
    assert!(out.rejected_intervals > 0);
    near(out.states_after[0].pose().position_world().z, 1.25, 2e-7);
    near(out.states_after[1].pose().position_world().z, 0.75, 2e-7);
    for i in 0..2 {
        vector(
            out.states_after[i]
                .center_of_mass_velocity_world(p[i])
                .unwrap(),
            Z.scale(0.5),
            1e-12,
        );
    }
    let JointFrictionTrajectoryMotion::Prismatic {
        displacement_m,
        distance_m,
        friction_impulse_ns,
        ..
    } = out.ports[0].motion
    else {
        panic!()
    };
    near(displacement_m, 0.5, 2e-12);
    near(distance_m, 0.5, 2e-12);
    near(friction_impulse_ns, -1.0, 2e-12);
    near(out.dissipated_j, 0.5, 2e-12);
    audit(&s, &p, free(), &out);
}

#[test]
fn g1_reversal_keeps_path_loss_and_advances_poses_through_the_actual_stop() {
    let (s, p, j) = slider(1.0);
    let out = held(&s, &p, &j, &[rail(1.0, -2.0)], 1.0, controls());
    let JointFrictionTrajectoryMotion::Prismatic {
        displacement_m,
        distance_m,
        friction_impulse_ns,
        drive_impulse_ns,
    } = out.ports[0].motion
    else {
        panic!()
    };
    near(displacement_m, -1.0 / 18.0, 2e-11);
    near(distance_m, 7.0 / 18.0, 2e-11);
    near(friction_impulse_ns, 1.0 / 3.0, 2e-11);
    near(drive_impulse_ns, -2.0, 2e-11);
    near(out.dissipated_j, 7.0 / 18.0, 2e-11);
    near(out.drive_work_j, 1.0 / 9.0, 2e-11);
    near(
        out.states_after[0].pose().position_world().z,
        17.0 / 36.0,
        4e-7,
    );
    near(
        out.states_after[1].pose().position_world().z,
        19.0 / 36.0,
        4e-7,
    );
    near(
        out.states_after[0]
            .center_of_mass_velocity_world(p[0])
            .unwrap()
            .z,
        1.0 / 6.0,
        2e-11,
    );
    near(
        out.states_after[1]
            .center_of_mass_velocity_world(p[1])
            .unwrap()
            .z,
        5.0 / 6.0,
        2e-11,
    );
    audit(&s, &p, free(), &out);
}

#[test]
fn g1_coupled_rotors_advance_independent_analytical_attitudes_and_actions() {
    let p = [
        properties(1.0, 2.0),
        properties(2.0, 3.0),
        properties(3.0, 5.0),
    ];
    let s = core::array::from_fn::<_, 3, _>(|i| {
        state(
            p[i],
            Vec3::ZERO,
            UnitQuaternion::IDENTITY,
            Vec3::ZERO,
            Z.scale([3.0, -1.0, 2.0][i]),
        )
    });
    let j = [
        ArticulatedJoint::Revolute(
            RevoluteJoint::new([0, 1], [Vec3::ZERO; 2], [Z; 2], 0.7).unwrap(),
        ),
        ArticulatedJoint::Revolute(
            RevoluteJoint::new([1, 2], [Vec3::ZERO; 2], [Z; 2], 1.3).unwrap(),
        ),
    ];
    let ports = [
        JointFrictionFlowPort::Revolute {
            joint_index: 1,
            capacity_nm: 2.0,
            drive_nm: 0.0,
        },
        JointFrictionFlowPort::Revolute {
            joint_index: 0,
            capacity_nm: 1.0,
            drive_nm: 0.0,
        },
    ];
    let out = held(&s, &p, &j, &ports, 1.0, controls());
    for i in 0..3 {
        let omega = out.states_after[i]
            .point_kinematics(p[i], Vec3::ZERO)
            .unwrap()
            .angular_velocity_world;
        vector(omega, Z.scale([2.5, 0.0, 1.6][i]), 2e-11);
        let expected = UnitQuaternion::from_axis_angle(Z, [2.75, -0.5, 1.8][i]).unwrap();
        vector(
            out.states_after[i]
                .pose()
                .orientation()
                .rotate_body_to_world(X),
            expected.rotate_body_to_world(X),
            2e-11,
        );
    }
    near(out.dissipated_j, 7.85, 2e-10);
    assert_eq!(
        out.ports.iter().map(|p| p.joint_index).collect::<Vec<_>>(),
        vec![0, 1]
    );
    for (p, expected) in out.ports.iter().zip([-3.25, 2.3]) {
        let JointFrictionTrajectoryMotion::Revolute { rotation_rad, .. } = p.motion else {
            panic!()
        };
        near(rotation_rad, expected, 2e-11);
    }
    audit(&s, &p, free(), &out);
}

#[test]
fn g2_source_breakpoints_resolve_a_short_drive_pulse_that_midpoint_sampling_would_skip() {
    let (s, p, j) = slider(0.0);
    let start = 2.0;
    let on = start + 1.0 / 64.0;
    let off = start + 3.0 / 128.0;
    let samples = Cell::new(0);
    let out = advance_joint_friction_trajectory(
        &s,
        &p,
        &j,
        &[rail(1.0, 0.0)],
        free(),
        &[Wrench::ZERO; 2],
        start,
        0.25,
        &[on, off],
        controls(),
        |time, _, ports| {
            let pulse = time > on && time < off;
            if pulse {
                samples.set(samples.get() + 1);
            }
            ports[0] = rail(1.0, if pulse { 9.0 } else { 0.0 });
            Ok::<_, ()>(())
        },
        || false,
    )
    .unwrap();
    assert!(samples.get() >= 3);
    let distance = 9.0 / 4096.0;
    let JointFrictionTrajectoryMotion::Prismatic {
        displacement_m,
        distance_m,
        friction_impulse_ns,
        drive_impulse_ns,
    } = out.ports[0].motion
    else {
        panic!()
    };
    near(displacement_m, distance, 2e-11);
    near(distance_m, distance, 2e-11);
    near(drive_impulse_ns, 9.0 / 128.0, 2e-11);
    near(friction_impulse_ns, -9.0 / 128.0, 2e-11);
    near(out.drive_work_j, distance, 2e-11);
    near(out.dissipated_j, distance, 2e-11);
    near(
        out.states_after[0].pose().position_world().z,
        distance / 2.0,
        4e-7,
    );
    near(
        out.states_after[1].pose().position_world().z,
        -distance / 2.0,
        4e-7,
    );
    for i in 0..2 {
        vector(
            out.states_after[i].linear_momentum_world(),
            Vec3::ZERO,
            2e-11,
        );
    }
    near(out.target_time_s, 2.25, 0.0);
    audit(&s, &p, free(), &out);
}

#[test]
fn g2_position_dependent_spring_drive_uses_the_moved_native_configuration() {
    let (mut s, p, j) = slider(0.0);
    s[0] = state(
        p[0],
        Z.scale(0.5),
        UnitQuaternion::IDENTITY,
        Vec3::ZERO,
        Vec3::ZERO,
    );
    s[1] = state(
        p[1],
        Z.scale(-0.5),
        UnitQuaternion::IDENTITY,
        Vec3::ZERO,
        Vec3::ZERO,
    );
    let moved = Cell::new(false);
    let out = advance_joint_friction_trajectory(
        &s,
        &p,
        &j,
        &[rail(0.1, 0.0)],
        free(),
        &[Wrench::ZERO; 2],
        0.0,
        0.4,
        &[],
        controls(),
        |_, native, ports| {
            let q = native[0].pose().position_world().z - native[1].pose().position_world().z;
            moved.set(moved.get() || q < 0.99);
            ports[0] = rail(0.1, -q);
            Ok::<_, ()>(())
        },
        || false,
    )
    .unwrap();
    assert!(moved.get());
    // During the entire interval v<0 and q''=-q+0.1, mu=1.
    let expected_q = 0.1 + 0.9 * 0.4_f64.cos();
    let expected_v = -0.9 * 0.4_f64.sin();
    near(
        out.states_after[0].pose().position_world().z
            - out.states_after[1].pose().position_world().z,
        expected_q,
        2e-6,
    );
    near(
        out.states_after[0]
            .center_of_mass_velocity_world(p[0])
            .unwrap()
            .z
            - out.states_after[1]
                .center_of_mass_velocity_world(p[1])
                .unwrap()
                .z,
        expected_v,
        2e-6,
    );
    near(out.dissipated_j, 0.1 * (1.0 - expected_q), 3e-7);
    audit(&s, &p, free(), &out);
}

#[test]
fn g2_rotating_rail_with_full_body_inertias_matches_an_independent_centrifugal_friction_ode() {
    let axes = [
        UnitQuaternion::from_axis_angle(Y, 0.37).unwrap(),
        UnitQuaternion::from_axis_angle(X, -0.29).unwrap(),
    ];
    let p = [
        MassProperties::from_principal_frame(2.0, Vec3::ZERO, Vec3::new(0.6, 0.8, 0.7), axes[0])
            .unwrap(),
        MassProperties::from_principal_frame(3.0, Vec3::ZERO, Vec3::new(0.9, 1.3, 1.1), axes[1])
            .unwrap(),
    ];
    let rail_frame = UnitQuaternion::from_axis_angle(Y, core::f64::consts::FRAC_PI_2).unwrap();
    let joint = ArticulatedJoint::Prismatic(
        PrismaticJoint::new(
            [0, 1],
            [Vec3::ZERO; 2],
            [compose(axes[0], rail_frame), compose(axes[1], rail_frame)],
            0.9,
        )
        .unwrap(),
    );
    let from_reference = |y: [f64; 4]| {
        let u = Vec3::new(y[2].cos(), y[2].sin(), 0.0);
        let d = u.scale(y[1]).add(Z.cross(u).scale(y[0] * y[3]));
        core::array::from_fn::<_, 2, _>(|i| {
            let fraction = [0.6, -0.4][i];
            let q = compose(
                UnitQuaternion::from_axis_angle(Z, y[2]).unwrap(),
                inverse(axes[i]),
            );
            state(
                p[i],
                u.scale(fraction * y[0]),
                q,
                d.scale(fraction),
                Z.scale(y[3]),
            )
        })
    };
    let original = [1.2, 0.15, 0.2, 0.6];
    let s = from_reference(original);
    let time = 0.4;
    let capacity = 0.25;
    let mu = 1.2;
    let rhs = |y: [f64; 4]| {
        [
            y[1],
            y[0] * y[3] * y[3] - capacity / mu,
            y[3],
            -2.0 * mu * y[0] * y[1] * y[3] / (1.8 + mu * y[0] * y[0]),
        ]
    };
    let mut reference = original;
    let h = time / 2048.0;
    for _ in 0..2048 {
        let a = rhs(reference);
        let b = rhs(core::array::from_fn(|i| reference[i] + 0.5 * h * a[i]));
        let c = rhs(core::array::from_fn(|i| reference[i] + 0.5 * h * b[i]));
        let d = rhs(core::array::from_fn(|i| reference[i] + h * c[i]));
        reference = core::array::from_fn(|i| {
            reference[i] + h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0
        });
    }
    let expected = from_reference(reference);
    let out = held(&s, &p, &[joint], &[rail(capacity, 0.0)], time, controls());
    for i in 0..2 {
        vector(
            out.states_after[i].pose().position_world(),
            expected[i].pose().position_world(),
            2e-6,
        );
        vector(
            out.states_after[i].linear_momentum_world(),
            expected[i].linear_momentum_world(),
            3e-6,
        );
        vector(
            out.states_after[i].angular_momentum_body(),
            expected[i].angular_momentum_body(),
            3e-6,
        );
        vector(
            out.states_after[i]
                .pose()
                .orientation()
                .rotate_body_to_world(X),
            expected[i].pose().orientation().rotate_body_to_world(X),
            2e-6,
        );
    }
    near(
        out.dissipated_j,
        capacity * (reference[0] - original[0]),
        5e-7,
    );
    let initial_axis = s[1]
        .pose()
        .orientation()
        .rotate_body_to_world(compose(axes[1], rail_frame).rotate_body_to_world(Z));
    let final_axis = out.states_after[1]
        .pose()
        .orientation()
        .rotate_body_to_world(compose(axes[1], rail_frame).rotate_body_to_world(Z));
    assert!((initial_axis.x - final_axis.x).abs() > 0.03);
    audit(&s, &p, free(), &out);
}

#[test]
fn g3_native_gravity_potential_and_common_external_work_stay_separate_from_friction_heat() {
    let (s, p, j) = slider(1.0);
    let integrator = RigidBodyIntegrator::new(Gravity::new(Y.scale(-9.0)).unwrap());
    let wrenches = [Wrench {
        force_world: X.scale(4.0),
        torque_body: Vec3::ZERO,
    }; 2];
    let out = advance_joint_friction_trajectory(
        &s,
        &p,
        &j,
        &[rail(1.0, 0.0)],
        integrator,
        &wrenches,
        0.0,
        0.5,
        &[],
        controls(),
        |_, _, _| Ok::<_, ()>(()),
        || false,
    )
    .unwrap();
    for i in 0..2 {
        near(out.states_after[i].pose().position_world().x, 0.25, 1e-12);
        near(out.states_after[i].pose().position_world().y, -1.125, 1e-12);
    }
    near(out.wrench_work_j, 2.0, 2e-12);
    near(out.dissipated_j, 0.375, 2e-12);
    near(out.mechanical_energy_change_j, 1.625, 2e-11);
    audit(&s, &p, integrator, &out);
}

#[test]
fn g4_late_source_failure_cancellation_and_refinement_refusal_publish_no_partial_motion() {
    let (s, p, j) = slider(1.0);
    let saved = s;
    let calls = Cell::new(0);
    let failure = advance_joint_friction_trajectory(
        &s,
        &p,
        &j,
        &[rail(1.0, 0.0)],
        free(),
        &[Wrench::ZERO; 2],
        0.0,
        2.0,
        &[],
        controls(),
        |_, _, _| {
            calls.set(calls.get() + 1);
            if calls.get() == 5 {
                Err("source exhausted")
            } else {
                Ok(())
            }
        },
        || false,
    );
    assert!(matches!(
        failure,
        Err(JointFrictionTrajectoryError::Source("source exhausted"))
    ));
    assert_eq!(s, saved);
    let source_calls = Cell::new(0);
    let cancel = advance_joint_friction_trajectory(
        &s,
        &p,
        &j,
        &[rail(1.0, 0.0)],
        free(),
        &[Wrench::ZERO; 2],
        0.0,
        2.0,
        &[],
        controls(),
        |_, _, _| {
            source_calls.set(source_calls.get() + 1);
            Ok::<_, ()>(())
        },
        || source_calls.get() >= 4,
    );
    assert!(cancel.is_err());
    assert_eq!(s, saved);
    let mut bounded = controls();
    bounded.maximum_split_trials = 3;
    assert!(held_result(&s, &p, &j, 2.0, &[], bounded).is_err());
    bounded = controls();
    bounded.minimum_step_s = 0.6;
    assert!(held_result(&s, &p, &j, 2.0, &[], bounded).is_err());
    bounded = controls();
    bounded.maximum_live_bytes = 1;
    assert!(held_result(&s, &p, &j, 2.0, &[], bounded).is_err());
    assert!(held_result(&s, &p, &j, 2.0, &[0.8, 0.5], controls()).is_err());
    assert!(held_result(&s, &p, &j, 2.0, &[2.0], controls()).is_err());
    let mutated = advance_joint_friction_trajectory(
        &s,
        &p,
        &j,
        &[rail(1.0, 0.0)],
        free(),
        &[Wrench::ZERO; 2],
        0.0,
        1.0,
        &[],
        controls(),
        |_, _, ports| {
            ports[0] = JointFrictionFlowPort::Prismatic {
                joint_index: 1,
                capacity_n: 1.0,
                drive_n: 0.0,
            };
            Ok::<_, ()>(())
        },
        || false,
    );
    assert!(matches!(
        mutated,
        Err(JointFrictionTrajectoryError::Invalid(
            "source changed joint identity"
        ))
    ));
    assert_eq!(s, saved);
}
fn held_result(
    s: &[RigidBodyState],
    p: &[MassProperties],
    j: &[ArticulatedJoint],
    duration: f64,
    breakpoints: &[f64],
    c: JointFrictionTrajectoryControls,
) -> Result<JointFrictionTrajectoryReport, JointFrictionTrajectoryError<()>> {
    advance_joint_friction_trajectory(
        s,
        p,
        j,
        &[rail(1.0, 0.0)],
        free(),
        &[Wrench::ZERO; 2],
        0.0,
        duration,
        breakpoints,
        c,
        |_, _, _| Ok(()),
        || false,
    )
}

#[test]
fn g3_represented_clock_interval_also_drives_native_motion_and_work() {
    let (s, p, j) = slider(0.0);
    let start = 9_007_199_254_740_992.0;
    let mut c = controls();
    c.maximum_step_s = 16.0;
    // At this epoch start+9 rounds to start+8. A reported eight-second
    // endpoint must not conceal nine seconds of applied force or drift.
    let out = advance_joint_friction_trajectory(
        &s,
        &p,
        &j,
        &[rail(0.0, 1.0)],
        free(),
        &[Wrench::ZERO; 2],
        start,
        9.0,
        &[],
        c,
        |_, _, _| Ok::<_, ()>(()),
        || false,
    )
    .unwrap();
    near(out.target_time_s - start, 8.0, 0.0);
    near(out.duration_s, 8.0, 0.0);
    let JointFrictionTrajectoryMotion::Prismatic {
        displacement_m,
        drive_impulse_ns,
        ..
    } = out.ports[0].motion
    else {
        panic!()
    };
    near(drive_impulse_ns, 8.0, 1e-12);
    near(displacement_m, 32.0, 1e-12);
    near(out.states_after[0].pose().position_world().z, 16.0, 1e-12);
    near(out.states_after[1].pose().position_world().z, -16.0, 1e-12);
    near(out.drive_work_j, 32.0, 1e-12);
    audit(&s, &p, free(), &out);
}
