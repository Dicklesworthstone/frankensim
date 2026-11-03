//! G0/G1/G3/G4: genuine supported mechanisms with no tangential mobility.

use fs_contact::frictional_multicontact::elliptic::EllipticFriction;
use fs_contact::frictional_multicontact::{
    FrictionalContactConfig, FrictionalContactDeterminacy, FrictionalContactError,
    NormalVelocityLaw,
};
use fs_ga::{Mat3, Se3, So3, Twist, Vec3};
use fs_mbd::articulated::{
    ArticulatedError, ArticulatedModel, BaseState, JointLimits, JointModel, Link,
    PrescribedBaseImpulseResponse, SpatialInertia,
};
use fs_mbd::articulated_contacts::friction::{
    ArticulatedFrictionalContact, ArticulatedFrictionalContactResult,
};
use fs_mbd::articulated_contacts::mixed::{
    ArticulatedContactBoundary as Boundary, MixedArticulatedContactSystem as System,
    MixedArticulatedImpulseResponse as Response, resolve_mixed_articulated_elliptic_contacts,
    resolve_mixed_articulated_frictional_contacts,
};
use fs_mbd::articulated_contacts::{ArticulatedContactError, ContactCounterpart, LinkContactPoint};

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-9 * expected.abs().max(1.0),
        "{actual:e} != {expected:e}"
    );
}

fn vector(actual: Vec3, expected: Vec3) {
    close(actual.x, expected.x);
    close(actual.y, expected.y);
    close(actual.z, expected.z);
}

fn inertia(center: Vec3) -> SpatialInertia {
    SpatialInertia::new(
        2.0,
        center,
        Mat3 {
            m: [1.0, 0.1, 0.2, 0.1, 1.5, 0.25, 0.2, 0.25, 2.0],
        },
    )
    .unwrap()
}

fn mechanism(joint: JointModel, center: Vec3) -> ArticulatedModel {
    ArticulatedModel::new(vec![Link::new(
        "supported mechanism",
        None,
        Se3::identity(),
        joint,
        inertia(center),
    )])
    .unwrap()
}

fn state<'a>(
    model: &'a ArticulatedModel,
    base: BaseState,
    position: &'a [f64],
    velocity: &'a [f64],
) -> System<'a> {
    System {
        model,
        boundary: Boundary::Prescribed(base),
        joint_position: position,
        joint_velocity: velocity,
    }
}

fn contact(
    pose: Se3,
    arm: Vec3,
    plane_velocity: Vec3,
    mu: f64,
    restitution: f64,
) -> ArticulatedFrictionalContact {
    ArticulatedFrictionalContact {
        id: 1,
        point: LinkContactPoint {
            system_index: 0,
            link_index: 0,
            point_body: arm,
        },
        other: ContactCounterpart::PrescribedPlane {
            point_world: pose.transform_point(arm).unwrap(),
            velocity_world: plane_velocity,
        },
        normal_world: Vec3::new(0.0, 1.0, 0.0),
        tangent_world: Vec3::new(1.0, 0.0, 0.0),
        activation_tolerance: 1e-9,
        friction_coefficient: mu,
        normal_law: NormalVelocityLaw::NewtonRestitution {
            coefficient: restitution,
        },
    }
}

fn supported(response: &Response) -> &PrescribedBaseImpulseResponse {
    let Response::Prescribed(response) = response else {
        panic!("lost native prescribed boundary")
    };
    response
}

fn audit<L>(result: &ArticulatedFrictionalContactResult<L, Response>, loss: f64) {
    close(result.energy_balance_residual, 0.0);
    close(result.endpoint_velocity_residual, 0.0);
    vector(result.linear_momentum_residual, Vec3::default());
    vector(result.angular_momentum_residual, Vec3::default());
    close(
        result.boundary_work,
        result.plane_work + result.support_work,
    );
    close(
        supported(&result.responses[0]).kinetic_energy_change,
        result.boundary_work - loss,
    );
}

#[test]
fn normal_slider_transmits_full_friction_wrench_and_signed_support_work() {
    let center = Vec3::new(0.3, -0.2, 0.4);
    let model = mechanism(
        JointModel::prismatic(Vec3::new(0.0, 1.0, 0.0), None).unwrap(),
        center,
    );
    let pose = Se3::from_parts(So3::identity(), Vec3::new(5.0, -3.0, 2.0)).unwrap();
    let arm = Vec3::new(0.7, -0.5, 0.9);
    let cases = [
        (Vec3::new(3.0, 0.4, -4.0), Vec3::new(0.0, 0.4, 0.0)),
        (Vec3::new(1.0, 0.4, 0.0), Vec3::new(4.0, 0.4, 4.0)),
    ];
    for (drive, belt) in cases {
        let base = BaseState::prescribed(pose, Twist::new(Vec3::default(), drive), Twist::zero());
        let systems = [state(&model, base, &[0.0], &[-2.0])];
        let slip = Vec3::new(drive.x - belt.x, 0.0, drive.z - belt.z);
        let expected = Vec3::new(-0.3 * slip.x, 6.0, -0.3 * slip.z);
        let momentum_change = Vec3::new(0.0, 6.0, 0.0);
        let support_force = momentum_change - expected;
        let support_couple = center.cross(momentum_change) - arm.cross(expected);
        for angle in [0.0_f64, 0.61] {
            let mut c = contact(pose, arm, belt, 0.25, 0.5);
            c.tangent_world = Vec3::new(angle.cos(), 0.0, angle.sin());
            let result = resolve_mixed_articulated_frictional_contacts(
                &systems,
                &[c],
                FrictionalContactConfig::default(),
                || false,
            )
            .unwrap();
            close(result.delassus[0], 0.5);
            assert!(result.delassus[1..].iter().all(|entry| *entry == 0.0));
            let response = supported(&result.responses[0]);
            close(response.joint_velocity_after[0], 1.0);
            assert_eq!(response.base_after, base);
            vector(response.external_impulse_world.force, expected);
            vector(response.support_impulse_body.force, support_force);
            vector(response.support_impulse_body.torque, support_couple);
            vector(response.linear_momentum_change_world, momentum_change);
            vector(
                response.angular_momentum_change_world,
                (pose.translation() + center).cross(momentum_change),
            );
            close(result.plane_work, expected.dot(belt));
            close(result.support_work, support_force.dot(drive));
            close(result.contact.dissipated_energy, 10.5);
            close(result.boundary_work, 9.9);
            assert!(matches!(
                result.contact.determinacy,
                FrictionalContactDeterminacy::UniqueVelocity(_)
            ));
            assert_eq!(
                result
                    .contact
                    .impulse_work_witness
                    .as_ref()
                    .unwrap()
                    .locked_tangent_substitutions
                    .len(),
                1
            );
            audit(&result, 10.5);
        }
    }
}

#[test]
fn full_tensor_hinge_keeps_radial_friction_and_all_support_reaction_components() {
    let model = mechanism(
        JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None).unwrap(),
        Vec3::new(1.0, 0.0, 0.0),
    );
    let arm = Vec3::new(2.0, 0.0, 0.0);
    for (drive, belt) in [(0.0, 3.0), (0.5, 3.5), (4.0, 1.0)] {
        let base = BaseState::prescribed(
            Se3::identity(),
            Twist::new(Vec3::default(), Vec3::new(drive, 0.0, 0.0)),
            Twist::zero(),
        );
        let systems = [state(&model, base, &[0.0], &[-2.0])];
        let c = contact(Se3::identity(), arm, Vec3::new(belt, 0.0, 0.0), 0.4, 0.5);
        let result = resolve_mixed_articulated_frictional_contacts(
            &systems,
            &[c],
            FrictionalContactConfig::default(),
            || false,
        )
        .unwrap();
        let radial = -2.4 * (drive - belt).signum();
        let response = supported(&result.responses[0]);
        close(result.delassus[0], 1.0); // r^2 / I_hinge = 4 / 4.
        close(result.contact.impulses[0][0], 6.0);
        close(result.contact.impulses[0][1], radial);
        close(response.joint_velocity_after[0], 1.0);
        vector(
            response.support_impulse_body.force,
            Vec3::new(-radial, 0.0, 0.0),
        );
        vector(
            response.support_impulse_body.torque,
            Vec3::new(0.6, 0.75, 0.0),
        );
        close(result.support_work, -radial * drive);
        close(result.plane_work, radial * belt);
        close(result.boundary_work, 7.2);
        close(result.contact.dissipated_energy, 13.2);
        audit(&result, 13.2);
    }
}

#[test]
fn static_gauge_and_tiny_nonzero_slider_slip_have_distinct_traction_authority() {
    let model = mechanism(
        JointModel::prismatic(Vec3::new(0.0, 1.0, 0.0), None).unwrap(),
        Vec3::default(),
    );
    for speed in [0.0, 1e-14] {
        let base = BaseState::prescribed(
            Se3::identity(),
            Twist::new(Vec3::default(), Vec3::new(speed, 0.0, 0.0)),
            Twist::zero(),
        );
        let systems = [state(&model, base, &[0.0], &[-1.0])];
        let result = resolve_mixed_articulated_frictional_contacts(
            &systems,
            &[contact(
                Se3::identity(),
                Vec3::default(),
                Vec3::default(),
                0.5,
                0.0,
            )],
            FrictionalContactConfig::default(),
            || false,
        )
        .unwrap();
        let traction = if speed == 0.0 { 0.0 } else { -1.0 };
        close(result.contact.impulses[0][0], 2.0);
        close(result.contact.impulses[0][1], traction);
        vector(
            supported(&result.responses[0]).support_impulse_body.force,
            Vec3::new(-traction, 0.0, 0.0),
        );
        assert!(matches!(
            result.contact.determinacy,
            FrictionalContactDeterminacy::UniqueVelocity(_)
        ));
        assert!(result.contact.impulse_work_witness.is_some());
        close(result.support_work, -traction * speed);
        audit(&result, result.contact.dissipated_energy);
    }
}

#[test]
fn directional_locked_slider_keeps_physical_ellipse_work_and_support_traction() {
    let model = mechanism(
        JointModel::prismatic(Vec3::new(0.0, 1.0, 0.0), None).unwrap(),
        Vec3::default(),
    );
    let drive = Vec3::new(3.0, 0.0, -4.0);
    let base = BaseState::prescribed(
        Se3::identity(),
        Twist::new(Vec3::default(), drive),
        Twist::zero(),
    );
    let systems = [state(&model, base, &[0.0], &[-2.0])];
    let geometry = contact(Se3::identity(), Vec3::default(), Vec3::default(), 0.0, 0.0);
    let c = ArticulatedFrictionalContact {
        id: geometry.id,
        point: geometry.point,
        other: geometry.other,
        normal_world: geometry.normal_world,
        tangent_world: geometry.tangent_world,
        activation_tolerance: geometry.activation_tolerance,
        normal_law: geometry.normal_law,
        friction_coefficient: EllipticFriction {
            coefficients_mu: [0.2, 0.5],
            principal_axis: [1.0, 0.0],
        },
    };
    let result = resolve_mixed_articulated_elliptic_contacts(
        &systems,
        &[c],
        FrictionalContactConfig::default(),
        || false,
    )
    .unwrap();
    let speed = 0.6_f64.hypot(2.0);
    let external = Vec3::new(-0.48 / speed, 4.0, 4.0 / speed);
    vector(
        supported(&result.responses[0]).external_impulse_world.force,
        external,
    );
    vector(
        supported(&result.responses[0]).support_impulse_body.force,
        Vec3::new(-external.x, 0.0, -external.z),
    );
    close(result.support_work, 4.0 * speed);
    close(result.contact.dissipated_energy, 4.0 + 4.0 * speed);
    audit(&result, result.contact.dissipated_energy);
}

#[test]
fn native_locked_slip_still_obeys_limits_cancellation_and_normal_mobility() {
    let limits = JointLimits::new(-1.0, 1.0, 1.0, 1e-12).unwrap();
    let model = mechanism(
        JointModel::prismatic(Vec3::new(0.0, 1.0, 0.0), Some(limits)).unwrap(),
        Vec3::default(),
    );
    let base = BaseState::prescribed(
        Se3::identity(),
        Twist::new(Vec3::default(), Vec3::new(1.0, 0.0, 0.0)),
        Twist::zero(),
    );
    let position = [0.0];
    let velocity = [-0.5];
    let systems = [state(&model, base, &position, &velocity)];
    let moving = contact(
        Se3::identity(),
        Vec3::default(),
        Vec3::new(0.0, 3.0, 0.0),
        0.2,
        0.0,
    );
    assert!(matches!(
        resolve_mixed_articulated_frictional_contacts(
            &systems,
            &[moving],
            FrictionalContactConfig::default(),
            || false,
        ),
        Err(ArticulatedContactError::Dynamics {
            error: ArticulatedError::JointVelocityLimitViolation { .. },
            ..
        })
    ));
    let stationary = contact(Se3::identity(), Vec3::default(), Vec3::default(), 0.2, 0.0);
    let mut polls = 0;
    let accepted = resolve_mixed_articulated_frictional_contacts(
        &systems,
        &[stationary],
        FrictionalContactConfig::default(),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    close(
        supported(&accepted.responses[0]).joint_velocity_after[0],
        0.0,
    );
    let mut seen = 0;
    assert!(matches!(
        resolve_mixed_articulated_frictional_contacts(
            &systems,
            &[stationary],
            FrictionalContactConfig::default(),
            || {
                seen += 1;
                seen == polls
            },
        ),
        Err(ArticulatedContactError::Cancelled)
    ));
    assert_eq!(velocity, [-0.5]);
    assert_eq!(position, [0.0]);
    assert_eq!(systems[0].boundary, Boundary::Prescribed(base));

    let fixed = mechanism(JointModel::FIXED, Vec3::default());
    let locked = [state(&fixed, base, &[], &[])];
    assert!(matches!(
        resolve_mixed_articulated_frictional_contacts(
            &locked,
            &[moving],
            FrictionalContactConfig::default(),
            || false,
        ),
        Err(ArticulatedContactError::FrictionLaw(
            FrictionalContactError::MobilityNotAdmitted(_)
        ))
    ));
}
