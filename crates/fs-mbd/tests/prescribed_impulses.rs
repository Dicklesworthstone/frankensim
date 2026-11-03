//! G0/G1/G3: supported articulated impulses against independent mass equations.

use fs_ga::{GaError, Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_mbd::articulated::{
    ArticulatedError, ArticulatedModel, BaseState, JointLimits, JointModel, Link,
    PrescribedBaseImpulseResponse, SpatialInertia, apply_prescribed_base_impulses,
    stationary_base_joint_momentum,
};

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-12,
        "{actual:e} != {expected:e}"
    );
}

fn vector(actual: Vec3, expected: Vec3) {
    assert!(
        (actual - expected).norm() < 2e-12,
        "{actual:?} != {expected:?}"
    );
}

fn inertia(mass: f64, center: Vec3, moment: f64) -> SpatialInertia {
    SpatialInertia::new(
        mass,
        center,
        Mat3 {
            m: [moment, 0.0, 0.0, 0.0, moment, 0.0, 0.0, 0.0, moment],
        },
    )
    .unwrap()
}

fn hinge(inertia: SpatialInertia, limits: Option<JointLimits>) -> ArticulatedModel {
    ArticulatedModel::new(vec![Link::new(
        "hinged body",
        None,
        Se3::identity(),
        JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), limits).unwrap(),
        inertia,
    )])
    .unwrap()
}

fn balanced(result: &PrescribedBaseImpulseResponse) {
    vector(result.linear_momentum_residual_world, Vec3::default());
    vector(result.angular_momentum_residual_world, Vec3::default());
    close(result.kinetic_work_residual, 0.0);
    close(
        result.kinetic_energy_change,
        result.momentum_after.kinetic_energy - result.momentum_before.kinetic_energy,
    );
    close(
        result.impulse_work,
        result.applied_impulse_work + result.support_work,
    );
}

#[test]
fn fixed_hinge_retains_full_inertia_and_all_six_support_reaction_components() {
    let model = hinge(
        SpatialInertia::new(
            2.0,
            Vec3::new(1.0, 0.0, 0.0),
            Mat3 {
                m: [1.0, 0.1, 0.2, 0.1, 1.5, 0.25, 0.2, 0.25, 2.0],
            },
        )
        .unwrap(),
        None,
    );
    // I_hinge,zz = I_com,zz + m*c_x² = 4. J_y=3 at x=2 gives
    // delta_qdot=6/4. The off-diagonal inertia creates a support couple;
    // the radial impulse is carried entirely by the support.
    let impulse = Wrench::new(Vec3::new(0.0, 0.0, 6.0), Vec3::new(0.4, 3.0, 0.0));
    let rotation = So3::exp(So3Tangent::new(Vec3::new(0.4, -0.2, 0.3))).unwrap();
    for pose in [
        Se3::identity(),
        Se3::from_parts(rotation, Vec3::new(3.0, -2.0, 4.0)).unwrap(),
    ] {
        let base = BaseState::stationary(pose);
        let result =
            apply_prescribed_base_impulses(&model, base, &[0.0], &[-2.0], &[0.0], &[impulse])
                .unwrap();
        close(result.joint_velocity_after[0], -0.5);
        assert_eq!(result.base_after, base);
        vector(result.support_impulse_body.force, Vec3::new(-0.4, 0.0, 0.0));
        vector(
            result.support_impulse_body.torque,
            Vec3::new(0.3, 0.375, 0.0),
        );
        let expected = result.support_impulse_body.transform_by(&pose);
        vector(result.support_impulse_world.force, expected.force);
        vector(result.support_impulse_world.torque, expected.torque);
        let expected_external = impulse.transform_by(&pose);
        vector(result.external_impulse_world.force, expected_external.force);
        vector(
            result.external_impulse_world.torque,
            expected_external.torque,
        );
        close(result.kinetic_energy_change, -7.5);
        close(result.support_work, 0.0);
        balanced(&result);
    }
}

#[test]
fn a_driven_base_supplies_signed_reaction_work_without_spurious_recoil() {
    let model = hinge(inertia(2.0, Vec3::new(1.0, 0.0, 0.0), 0.5), None);
    let impulse = Wrench::new(Vec3::new(0.0, 0.0, 6.0), Vec3::new(0.0, 3.0, 0.0));
    for base_speed in [-2.0, 2.0] {
        let base = BaseState::prescribed(
            Se3::identity(),
            Twist::new(Vec3::new(0.0, 0.0, 0.3), Vec3::new(0.0, base_speed, 0.0)),
            // Smooth acceleration is irrelevant to this zero-duration jump.
            Twist::new(Vec3::new(5.0, 1.0, 2.0), Vec3::new(-30.0, 40.0, 10.0)),
        );
        let result =
            apply_prescribed_base_impulses(&model, base, &[0.0], &[-2.0], &[0.0], &[impulse])
                .unwrap();
        assert_eq!(result.base_after, base);
        close(result.joint_velocity_after[0], 0.4);
        vector(result.support_impulse_body.force, Vec3::new(0.0, 1.8, 0.0));
        vector(result.support_impulse_body.torque, Vec3::default());
        close(result.support_work, 1.8 * base_speed);
        close(result.applied_impulse_work, 3.0 * (base_speed - 1.0));
        close(result.kinetic_energy_change, 4.8 * base_speed - 3.0);
        balanced(&result);
    }
}

#[test]
fn root_actuator_impulse_transmits_the_correct_support_force_and_torque() {
    let model = hinge(inertia(2.0, Vec3::new(1.0, 0.0, 0.0), 0.5), None);
    let base = BaseState::prescribed(
        Se3::identity(),
        Twist::new(Vec3::new(0.0, 0.0, 0.5), Vec3::new(0.0, -1.0, 0.0)),
        Twist::zero(),
    );
    let result =
        apply_prescribed_base_impulses(&model, base, &[0.0], &[1.0], &[-1.5], &[Wrench::default()])
            .unwrap();
    close(result.joint_velocity_after[0], 0.4);
    vector(result.external_impulse_world.force, Vec3::default());
    vector(result.external_impulse_world.torque, Vec3::default());
    vector(result.support_impulse_body.force, Vec3::new(0.0, -1.2, 0.0));
    vector(
        result.support_impulse_body.torque,
        Vec3::new(0.0, 0.0, -1.5),
    );
    close(result.applied_impulse_work, -1.05);
    close(result.support_work, 0.45);
    close(result.kinetic_energy_change, -0.6);
    balanced(&result);
}

#[test]
fn coupled_slider_impulses_match_the_analytic_joint_mass_system() {
    let joint = JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None).unwrap();
    let model = ArticulatedModel::new(vec![
        Link::new(
            "first",
            None,
            Se3::identity(),
            joint,
            inertia(2.0, Vec3::default(), 1.0),
        ),
        Link::new(
            "second",
            Some(0),
            Se3::identity(),
            joint,
            inertia(3.0, Vec3::default(), 1.0),
        ),
    ])
    .unwrap();
    let body_impulses = [
        Wrench::new(Vec3::default(), Vec3::new(4.0, 0.0, 0.0)),
        Wrench::new(Vec3::default(), Vec3::new(-1.0, 0.0, 0.0)),
    ];
    let result = apply_prescribed_base_impulses(
        &model,
        BaseState::stationary(Se3::identity()),
        &[0.0, 0.0],
        &[0.0, 0.0],
        &[2.0, 3.0],
        &body_impulses,
    )
    .unwrap();
    // M=[[5,3],[3,3]], joint + J_body impulses = [5,2].
    close(result.joint_velocity_after[0], 1.5);
    close(result.joint_velocity_after[1], -5.0 / 6.0);
    let momenta =
        stationary_base_joint_momentum(&model, &[0.0, 0.0], &result.joint_velocity_after).unwrap();
    close(momenta[0], 5.0);
    close(momenta[1], 2.0);
    vector(result.support_impulse_body.force, Vec3::new(2.0, 0.0, 0.0));
    vector(result.support_impulse_body.torque, Vec3::default());
    close(result.kinetic_energy_change, 35.0 / 12.0);
    balanced(&result);
}

#[test]
fn a_fully_fixed_tree_returns_the_whole_reaction_and_zero_net_work() {
    let local_pose = Se3::from_parts(
        So3::exp(So3Tangent::new(Vec3::new(0.1, 0.2, -0.3))).unwrap(),
        Vec3::new(1.0, -2.0, 3.0),
    )
    .unwrap();
    let model = ArticulatedModel::new(vec![Link::new(
        "fixed",
        None,
        local_pose,
        JointModel::FIXED,
        inertia(2.0, Vec3::new(0.2, 0.1, -0.3), 1.0),
    )])
    .unwrap();
    let base = BaseState::prescribed(
        Se3::identity(),
        Twist::new(Vec3::new(0.2, 0.3, -0.1), Vec3::new(1.0, -0.5, 0.6)),
        Twist::zero(),
    );
    let impulse = Wrench::new(Vec3::new(0.2, -0.1, 0.4), Vec3::new(-1.0, 2.0, 3.0));
    let result = apply_prescribed_base_impulses(&model, base, &[], &[], &[], &[impulse]).unwrap();
    let external = impulse.transform_by(&local_pose);
    vector(
        result.support_impulse_body.force,
        external.force.scale(-1.0),
    );
    vector(
        result.support_impulse_body.torque,
        external.torque.scale(-1.0),
    );
    assert_eq!(result.momentum_before, result.momentum_after);
    assert_eq!(result.base_after, base);
    assert!(result.joint_velocity_after.is_empty());
    close(result.support_work, -result.applied_impulse_work);
    balanced(&result);
}

#[test]
fn impulse_units_ignore_effort_limits_but_actual_endpoint_limits_remain_atomic() {
    let limits = JointLimits::new(-1.0, 1.0, 1.0, 1e-6).unwrap();
    let model = hinge(inertia(1.0, Vec3::default(), 1e-3), Some(limits));
    let base = BaseState::stationary(Se3::identity());
    let rate = [0.0];
    let result =
        apply_prescribed_base_impulses(&model, base, &[0.0], &rate, &[5e-4], &[Wrench::default()])
            .unwrap();
    close(result.joint_velocity_after[0], 0.5);
    balanced(&result);
    assert!(matches!(
        apply_prescribed_base_impulses(&model, base, &[0.0], &rate, &[2e-3], &[Wrench::default()]),
        Err(ArticulatedError::JointVelocityLimitViolation { .. })
    ));
    assert!(matches!(
        apply_prescribed_base_impulses(&model, base, &[0.0], &rate, &[], &[Wrench::default()]),
        Err(ArticulatedError::LengthMismatch {
            field: "joint_impulses",
            ..
        })
    ));
    assert!(
        apply_prescribed_base_impulses(
            &model,
            base,
            &[0.0],
            &rate,
            &[f64::NAN],
            &[Wrench::default()]
        )
        .is_err()
    );
    assert_eq!(rate, [0.0]);
    assert_eq!(base, BaseState::stationary(Se3::identity()));
}

#[test]
fn a_nonzero_joint_jump_lost_in_a_large_finite_rate_refuses_atomically() {
    let model = hinge(inertia(1.0, Vec3::default(), 1.0), None);
    let base = BaseState::stationary(Se3::identity());
    let rate = [1e16];
    // The unit joint inertia resolves each +/-0.5 jump exactly, but neither
    // increment can change this large f64 background rate.
    for impulse in [-0.5, 0.5] {
        assert!(matches!(
            apply_prescribed_base_impulses(
                &model,
                base,
                &[0.0],
                &rate,
                &[impulse],
                &[Wrench::default()],
            ),
            Err(ArticulatedError::Geometry(GaError::InvalidRepresentation {
                context: "prescribed_impulse.joint_velocity_jump",
                defect: 0.5,
                tolerance: 0.0,
            }))
        ));
    }
    assert_eq!(rate, [1e16]);
    // A blocked translational direction has an IDEAL zero joint jump. The
    // actual unchanged state and balancing support impulse remain legitimate.
    let held = apply_prescribed_base_impulses(
        &model,
        base,
        &[0.0],
        &rate,
        &[0.0],
        &[Wrench::new(Vec3::default(), Vec3::new(1.0, 0.0, 0.0))],
    )
    .unwrap();
    assert_eq!(held.joint_velocity_after, rate);
    vector(held.support_impulse_body.force, Vec3::new(-1.0, 0.0, 0.0));
    balanced(&held);
}

#[test]
fn a_resolved_joint_jump_erased_by_a_large_prescribed_twist_refuses_atomically() {
    let model = hinge(inertia(1.0, Vec3::default(), 1.0), None);
    let base = BaseState::prescribed(
        Se3::identity(),
        Twist::new(Vec3::new(0.0, 0.0, 1e16), Vec3::default()),
        Twist::zero(),
    );
    let rate = [0.0];
    // The generalized rate can store +/-0.5, but native link kinematics would
    // add it to 1e16 and report exactly unchanged physical velocity/momentum.
    for impulse in [-0.5, 0.5] {
        assert!(matches!(
            apply_prescribed_base_impulses(
                &model,
                base,
                &[0.0],
                &rate,
                &[impulse],
                &[Wrench::default()],
            ),
            Err(ArticulatedError::Geometry(GaError::InvalidRepresentation {
                context: "prescribed_impulse.link_velocity_jump",
                defect: 0.5,
                tolerance: 0.0,
            }))
        ));
    }
    assert_eq!(rate, [0.0]);
    let held = apply_prescribed_base_impulses(
        &model,
        base,
        &[0.0],
        &rate,
        &[0.0],
        &[Wrench::new(Vec3::default(), Vec3::new(1.0, 0.0, 0.0))],
    )
    .unwrap();
    assert_eq!(held.base_after, base);
    assert_eq!(held.joint_velocity_after, rate);
    vector(held.support_impulse_body.force, Vec3::new(-1.0, 0.0, 0.0));
    balanced(&held);
}
