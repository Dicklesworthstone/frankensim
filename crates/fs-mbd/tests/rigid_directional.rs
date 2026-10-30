//! Native full-inertia directional friction: synthetic mechanical references.
use fs_contact::interface_binding::thermal_friction::elliptic::EllipticFrictionMode;
use fs_exec::CancelGate;
use fs_mbd::tangential_friction::RigidContactBinding;
use fs_mbd::tangential_friction::directional::*;
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_tribo::partial_slip::TangentFrame;
use fs_tribo::tangential::elliptic::EllipticCoulomb;
use fs_tribo::{TangentialLaw, TangentialSolveConfig, TangentialStep};

#[path = "rigid_directional/trajectory.rs"]
mod trajectory;

#[path = "rigid_directional/driven.rs"]
mod driven;

fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() < 3e-9 * a.abs().max(b.abs()).max(1.0),
        "{a:e} != {b:e}"
    );
}
fn vector(a: Vec3, b: Vec3) {
    close(a.x, b.x);
    close(a.y, b.y);
    close(a.z, b.z);
}
fn controls() -> RigidDirectionalControls {
    RigidDirectionalControls {
        maximum_bodies: 8,
        geometry_tolerance_m: 1e-10,
        velocity_tolerance_m_s: 1e-8,
        linear_momentum_tolerance_ns: 1e-9,
        angular_momentum_tolerance_nms: 1e-9,
        energy_tolerance_j: 1e-8,
        contact: TangentialSolveConfig {
            max_iterations: 128,
            maximum_force_n: 1e6,
            force_absolute_tolerance_n: 1e-11,
            force_relative_tolerance: 1e-13,
        },
    }
}
fn contact(arms: [Vec3; 2]) -> RigidDirectionalContact {
    RigidDirectionalContact {
        binding: RigidContactBinding {
            body_a: 0,
            body_b: 1,
            arm_a_body: arms[0],
            arm_b_body: arms[1],
        },
        frame: TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap(),
        normal_load_n: 10.0,
        contact_temperature_k: Some(300.0),
        nominal_pressure_pa: Some(1e4),
    }
}
fn full_pair(
    gauged: bool,
) -> (
    [RigidBodyState; 3],
    [MassProperties; 3],
    RigidDirectionalContact,
) {
    let moments = [Vec3::new(1.0, 1.5, 2.0), Vec3::new(0.8, 1.2, 1.6)];
    let rotations = [
        UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, 3.0), 0.7).unwrap(),
        UnitQuaternion::from_axis_angle(Vec3::new(2.0, -1.0, 1.0), -0.4).unwrap(),
    ];
    let positions = [Vec3::new(-0.2, -0.3, -0.7), Vec3::new(0.4, -0.1, 0.2)];
    let linear = [Vec3::new(2.0, 0.8, 0.0), Vec3::new(-0.4, 0.1, 0.0)];
    let angular = [Vec3::new(0.1, 0.2, -0.1), Vec3::new(-0.1, 0.2, 0.05)];
    let arms = [0, 1].map(|i| {
        if gauged {
            rotations[i].rotate_world_to_body(positions[i].scale(-1.0))
        } else {
            positions[i].scale(-1.0)
        }
    });
    let props = [0, 1].map(|i| {
        if gauged {
            MassProperties::new(2.0 + i as f64, Vec3::ZERO, moments[i]).unwrap()
        } else {
            MassProperties::from_principal_frame(
                2.0 + i as f64,
                Vec3::ZERO,
                moments[i],
                rotations[i],
            )
            .unwrap()
        }
    });
    let bodies = [0, 1].map(|i| {
        RigidBodyState::new(
            Pose::new(
                positions[i],
                if gauged {
                    rotations[i]
                } else {
                    UnitQuaternion::IDENTITY
                },
            )
            .unwrap(),
            linear[i],
            if gauged {
                rotations[i].rotate_world_to_body(angular[i])
            } else {
                angular[i]
            },
        )
        .unwrap()
    });
    (
        [
            bodies[0],
            bodies[1],
            RigidBodyState::new(Pose::identity(), Vec3::new(3.0, 2.0, 1.0), Vec3::ZERO).unwrap(),
        ],
        [
            props[0],
            props[1],
            MassProperties::new(5.0, Vec3::ZERO, Vec3::new(2.0, 3.0, 4.0)).unwrap(),
        ],
        contact(arms),
    )
}
fn relative(
    bodies: &[RigidBodyState],
    props: &[MassProperties],
    c: RigidDirectionalContact,
) -> Vec3 {
    bodies[0]
        .point_kinematics(props[0], c.binding.arm_a_body)
        .unwrap()
        .point_velocity_world
        .sub(
            bodies[1]
                .point_kinematics(props[1], c.binding.arm_b_body)
                .unwrap()
                .point_velocity_world,
        )
}

#[test]
fn actual_offcentre_full_tensor_mobility_drives_work_and_both_body_impulses() {
    let (mut bodies, props, c) = full_pair(false);
    let old = bodies;
    let law = EllipticCoulomb::new([0.9, 0.8], [0.2, 0.4], [0.6, 0.8]).unwrap();
    let mut state = law.initial_state();
    let report = apply_rigid_directional(
        &mut bodies,
        &props,
        c,
        &law,
        &mut state,
        0.002,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert!(report.mobility_per_kg[0][1].abs() > 1e-3);
    // Independent response columns use actual unit impulses, not the Gram factor.
    let before = relative(&old, &props, c);
    for (j, axis) in [Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0)]
        .into_iter()
        .enumerate()
    {
        let mut probe = old;
        probe[0] = old[0]
            .apply_impulse_at_body_point(props[0], c.binding.arm_a_body, axis)
            .unwrap()
            .state_after;
        probe[1] = old[1]
            .apply_impulse_at_body_point(props[1], c.binding.arm_b_body, axis.scale(-1.0))
            .unwrap()
            .state_after;
        let delta = relative(&probe, &props, c).sub(before);
        close(report.mobility_per_kg[0][j], delta.x);
        close(report.mobility_per_kg[1][j], delta.y);
    }
    close(
        report.kinetic_energy_change_j + report.contact.dissipated_j,
        0.0,
    );
    assert!(report.contact.dissipated_j > 0.0);
    close(report.audit.endpoint_velocity_residual_m_s, 0.0);
    vector(report.audit.linear_momentum_change_ns, Vec3::ZERO);
    vector(report.audit.angular_momentum_change_nms, Vec3::ZERO);
    for x in report.audit.body_linear_momentum_residual_ns {
        vector(x, Vec3::ZERO);
    }
    for x in report.audit.body_angular_momentum_residual_nms {
        vector(x, Vec3::ZERO);
    }
    assert_eq!(bodies[2], old[2]);
    for i in 0..2 {
        assert_eq!(bodies[i].pose(), old[i].pose());
    }
    assert_eq!(state, report.contact.state);
}

#[test]
fn changing_only_body_coordinate_frames_preserves_world_motion_and_heat() {
    let (mut a, pa, ca) = full_pair(false);
    let (mut b, pb, cb) = full_pair(true);
    let law = EllipticCoulomb::new([0.9, 0.8], [0.2, 0.4], [0.6, 0.8]).unwrap();
    let mut sa = law.initial_state();
    let mut sb = law.initial_state();
    let ra = apply_rigid_directional(
        &mut a,
        &pa,
        ca,
        &law,
        &mut sa,
        0.002,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let rb = apply_rigid_directional(
        &mut b,
        &pb,
        cb,
        &law,
        &mut sb,
        0.002,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    close(ra.contact.dissipated_j, rb.contact.dissipated_j);
    for i in 0..2 {
        close(ra.contact.traction_n[i], rb.contact.traction_n[i]);
        vector(a[i].linear_momentum_world(), b[i].linear_momentum_world());
        vector(
            a[i].pose()
                .orientation()
                .rotate_body_to_world(a[i].angular_momentum_body()),
            b[i].pose()
                .orientation()
                .rotate_body_to_world(b[i].angular_momentum_body()),
        );
        for j in 0..2 {
            close(ra.mobility_per_kg[i][j], rb.mobility_per_kg[i][j]);
        }
    }
    vector(relative(&a, &pa, ca), relative(&b, &pb, cb));
}

#[test]
fn exact_principal_stop_retains_heat_beside_large_common_energy_and_refuses_overstep() {
    let props = [
        MassProperties::new(2.0, Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).unwrap(),
        MassProperties::new(3.0, Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).unwrap(),
    ];
    let law = EllipticCoulomb::new([0.9, 0.8], [0.2, 0.4], [1.0, 0.0]).unwrap();
    let c = contact([Vec3::ZERO; 2]);
    for common in [0.0, 1e12] {
        let old = [
            RigidBodyState::new(
                Pose::identity(),
                Vec3::new(2.0, 0.0, 2.0 * common),
                Vec3::ZERO,
            )
            .unwrap(),
            RigidBodyState::new(
                Pose::identity(),
                Vec3::new(0.0, 0.0, 3.0 * common),
                Vec3::ZERO,
            )
            .unwrap(),
        ];
        let mut bodies = old;
        let mut state = law.initial_state();
        let report = apply_rigid_directional(
            &mut bodies,
            &props,
            c,
            &law,
            &mut state,
            0.6,
            controls(),
            &CancelGate::new(),
        )
        .unwrap();
        close(report.slip_velocity_after_m_s[0], 0.0);
        close(report.contact.dissipated_j, 0.6);
        close(report.kinetic_energy_change_j, -0.6);
        for i in 0..2 {
            close(
                bodies[i].center_of_mass_velocity_world(props[i]).unwrap().x,
                0.4,
            );
        }
        for dt in [0.7, 1.2] {
            let mut body = old;
            let mut history = law.initial_state();
            let saved = history.clone();
            assert!(matches!(
                apply_rigid_directional(
                    &mut body,
                    &props,
                    c,
                    &law,
                    &mut history,
                    dt,
                    controls(),
                    &CancelGate::new()
                ),
                Err(RigidDirectionalError::UnresolvedSlipEvent { .. })
            ));
            assert_eq!(body, old);
            assert_eq!(history, saved);
        }
    }
}

#[test]
fn source_callback_failure_and_cancel_publish_neither_body_nor_history() {
    let (mut bodies, props, c) = full_pair(false);
    let old = bodies;
    let mut history = 17usize;
    assert!(matches!(
        apply_rigid_directional_with(
            &mut bodies,
            &props,
            c,
            &mut history,
            0.002,
            controls(),
            |_, _, _, _| Err::<TangentialStep<usize>, _>("original source unavailable"),
            &CancelGate::new()
        ),
        Err(RigidDirectionalError::Source("original source unavailable"))
    ));
    assert_eq!(bodies, old);
    assert_eq!(history, 17);
    let law = EllipticCoulomb::new([0.9, 0.8], [0.2, 0.4], [0.6, 0.8]).unwrap();
    let mut state = law.initial_state();
    let initial = state.clone();
    let gate = CancelGate::new();
    assert!(matches!(
        apply_rigid_directional_with(
            &mut bodies,
            &props,
            c,
            &mut state,
            0.002,
            controls(),
            |s, input, compliance, mode| {
                assert_eq!(input.compliance_m_per_n, 0.0);
                assert_eq!(mode, EllipticFrictionMode::Kinetic);
                let step = law.advance_kinetic_with_compliance(s, input, compliance);
                gate.request();
                step
            },
            &gate
        ),
        Err(RigidDirectionalError::Cancelled)
    ));
    assert_eq!(bodies, old);
    assert_eq!(state, initial);
    let mut foreign = c;
    foreign.binding.arm_a_body = Vec3::ZERO;
    assert!(matches!(
        apply_rigid_directional(
            &mut bodies,
            &props,
            foreign,
            &law,
            &mut state,
            0.002,
            controls(),
            &CancelGate::new()
        ),
        Err(RigidDirectionalError::Invalid(_))
    ));
    assert_eq!(bodies, old);
    assert_eq!(state, initial);
}
