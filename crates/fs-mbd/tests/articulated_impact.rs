//! Single-contact impacts through the native articulated inertia owner.

use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3};
use fs_mbd::articulated::{
    ArticulatedError, ArticulatedModel, BaseState, FreeFloatingBaseState, JointLimits, JointModel,
    Link, SpatialInertia, forward_kinematics, free_floating_momentum,
};
use fs_mbd::articulated_impact::{
    ArticulatedImpactError, ArticulatedPointPlaneContact, resolve_articulated_plane_impact,
};

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 3e-12,
        "{actual:e} != {expected:e}"
    );
}

fn close_vector(actual: Vec3, expected: Vec3) {
    assert!(
        (actual - expected).norm() < 3e-12,
        "{actual:?} != {expected:?}"
    );
}

fn inertia(mass: f64, center: Vec3, moments: [f64; 3]) -> SpatialInertia {
    SpatialInertia::new(
        mass,
        center,
        Mat3 {
            m: [
                moments[0], 0.0, 0.0, 0.0, moments[1], 0.0, 0.0, 0.0, moments[2],
            ],
        },
    )
    .unwrap()
}

fn single_body(inertia: SpatialInertia) -> ArticulatedModel {
    ArticulatedModel::new(vec![Link::new(
        "free body",
        None,
        Se3::identity(),
        JointModel::FIXED,
        inertia,
    )])
    .unwrap()
}

fn plane_contact(
    point_body: Vec3,
    point_world: Vec3,
    restitution: f64,
) -> ArticulatedPointPlaneContact {
    ArticulatedPointPlaneContact {
        link_index: 0,
        point_body,
        plane_point_world: point_world,
        normal_world: Vec3::new(0.0, 1.0, 0.0),
        plane_velocity_world: Vec3::default(),
        activation_tolerance: 1e-8,
        restitution,
    }
}

/// G0/G1: a rotated body with an offset reference origin has contact response
/// 1/m + (r_com cross n)^T I_com^-1 (r_com cross n). The independent Cartesian
/// calculation checks the full velocity change, momentum, and restitution
/// energy loss for plastic, intermediate, and elastic impacts.
#[test]
fn off_center_free_body_matches_independent_effective_mass_and_energy() {
    let mass = 2.5;
    let center = Vec3::new(0.2, -0.1, 0.3);
    let moments = [1.0, 1.5, 2.0];
    let model = single_body(inertia(mass, center, moments));
    let rotation = So3::exp(So3Tangent::new(Vec3::new(0.3, -0.2, 0.7))).unwrap();
    let pose = Se3::from_parts(rotation, Vec3::new(1.0, -2.0, 0.5)).unwrap();
    let angular = Vec3::new(0.3, -0.2, 0.1);
    let linear = Vec3::new(0.1, -0.8, 0.2);
    let base = FreeFloatingBaseState::new(pose, Twist::new(angular, linear));
    let point = center + Vec3::new(0.4, 0.2, -0.3);
    let kinematics = forward_kinematics(
        &model,
        BaseState::prescribed(pose, base.twist_body, Twist::zero()),
        &[],
        &[],
    )
    .unwrap();
    let point_world = kinematics.world_from_link[0]
        .transform_point(point)
        .unwrap();
    let normal_world = rotation.rotate(Vec3::new(0.0, 1.0, 0.0)).unwrap();
    let inverse_mass = 1.0 / mass + 0.3 * 0.3 / moments[0] + 0.4 * 0.4 / moments[2];
    let normal_speed = (linear + angular.cross(point)).y;
    let com_velocity = linear + angular.cross(center);
    let initial_energy = 0.5 * mass * com_velocity.dot(com_velocity)
        + 0.5
            * (moments[0] * angular.x.powi(2)
                + moments[1] * angular.y.powi(2)
                + moments[2] * angular.z.powi(2));
    for restitution in [0.0, 0.4, 1.0] {
        let contact = ArticulatedPointPlaneContact {
            normal_world: normal_world.scale(7.0),
            ..plane_contact(point, point_world, restitution)
        };
        let result = resolve_articulated_plane_impact(&model, base, &[], &[], contact).unwrap();
        let impulse = -(1.0 + restitution) * normal_speed / inverse_mass;
        let delta_angular = Vec3::new(0.3 * impulse / moments[0], 0.0, 0.4 * impulse / moments[2]);
        let delta_linear = Vec3::new(0.0, impulse / mass, 0.0) - delta_angular.cross(center);
        let impulse_world = normal_world.scale(impulse);
        let expected_loss =
            0.5 * (1.0 - restitution * restitution) * normal_speed * normal_speed / inverse_mass;
        close(result.inverse_effective_mass, inverse_mass);
        close(result.relative_normal_speed_before, normal_speed);
        close(
            result.relative_normal_speed_after,
            -restitution * normal_speed,
        );
        close(result.normal_law.normal_impulse, impulse);
        close(result.normal_law.dissipated_energy, expected_loss);
        close(
            result.response.momentum_before.kinetic_energy,
            initial_energy,
        );
        close(
            result.response.momentum_after.kinetic_energy,
            initial_energy - expected_loss,
        );
        close(result.plane_work, 0.0);
        close(result.energy_balance_residual, 0.0);
        close_vector(
            result.response.base_after.twist_body.angular,
            angular + delta_angular,
        );
        close_vector(
            result.response.base_after.twist_body.linear,
            linear + delta_linear,
        );
        close_vector(result.response.external_linear_impulse_world, impulse_world);
        close_vector(
            result.response.external_angular_impulse_world,
            point_world.cross(impulse_world),
        );
        close_vector(
            result.response.momentum_after.linear_momentum_world
                - result.response.momentum_before.linear_momentum_world,
            impulse_world,
        );
        close_vector(
            result.response.momentum_after.angular_momentum_world
                - result.response.momentum_before.angular_momentum_world,
            point_world.cross(impulse_world),
        );
        assert_eq!(result.response.base_after.world_from_base, pose);
        assert!(result.response.joint_velocity_after.is_empty());
        assert!(result.normal_law.active && result.normal_law.impacted);
    }
}

fn two_link_model(child_inertia: f64, speed_limit: f64) -> ArticulatedModel {
    let limits = JointLimits::new(-2.0, 2.0, speed_limit, 1e-9).unwrap();
    let joint = JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), Some(limits)).unwrap();
    ArticulatedModel::new(vec![
        Link::new(
            "base",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(2.0, Vec3::default(), [1.0; 3]),
        ),
        Link::new(
            "child",
            Some(0),
            Se3::identity(),
            joint,
            inertia(3.0, Vec3::default(), [child_inertia; 3]),
        ),
    ])
    .unwrap()
}

/// G1: coincident-center bodies share linear recoil, while a free revolute
/// axis lets the struck child alone acquire contact spin. For contact arm x=1
/// and normal +y, w=1/(m_base+m_child)+1/I_child, not the isolated child mass.
#[test]
fn child_impact_recoils_the_free_base_through_the_coupled_mass() {
    let model = two_link_model(2.0, 5.0);
    let base = FreeFloatingBaseState::new(
        Se3::identity(),
        Twist::new(Vec3::default(), Vec3::new(0.0, -0.6, 0.0)),
    );
    let point = Vec3::new(1.0, 0.0, 0.0);
    let contact = ArticulatedPointPlaneContact {
        link_index: 1,
        ..plane_contact(point, point, 0.4)
    };
    let result = resolve_articulated_plane_impact(&model, base, &[0.0], &[-0.1], contact).unwrap();
    close(result.inverse_effective_mass, 0.7);
    close(result.normal_law.normal_impulse, 1.4);
    close(result.relative_normal_speed_before, -0.7);
    close(result.relative_normal_speed_after, 0.28);
    close_vector(
        result.response.base_after.twist_body.angular,
        Vec3::default(),
    );
    close_vector(
        result.response.base_after.twist_body.linear,
        Vec3::new(0.0, -0.32, 0.0),
    );
    close(result.response.joint_velocity_after[0], 0.6);
    close(result.response.momentum_before.kinetic_energy, 0.91);
    close(result.response.momentum_after.kinetic_energy, 0.616);
    close(result.normal_law.dissipated_energy, 0.294);
    close(result.energy_balance_residual, 0.0);
    close_vector(
        result.response.momentum_after.linear_momentum_world
            - result.response.momentum_before.linear_momentum_world,
        Vec3::new(0.0, 1.4, 0.0),
    );
    close_vector(
        result.response.momentum_after.angular_momentum_world
            - result.response.momentum_before.angular_momentum_world,
        Vec3::new(0.0, 0.0, 1.4),
    );
    let measured = free_floating_momentum(
        &model,
        result.response.base_after,
        &[0.0],
        &result.response.joint_velocity_after,
    )
    .unwrap();
    assert_eq!(measured, result.response.momentum_after);
}

/// G1/G3: a prescribed moving plane can supply or absorb energy. Its signed
/// impulse work is separate from nonnegative restitution loss; tangential
/// plane velocity supplies no frictional work in this normal-only law.
#[test]
fn moving_plane_work_closes_actual_kinetic_energy_change() {
    let model = single_body(inertia(2.0, Vec3::default(), [1.0; 3]));
    for (body_speed, plane_speed) in [(-0.3, 0.5), (-1.0, -0.4)] {
        let base = FreeFloatingBaseState::new(
            Se3::identity(),
            Twist::new(Vec3::default(), Vec3::new(0.0, body_speed, 0.0)),
        );
        let contact = ArticulatedPointPlaneContact {
            plane_velocity_world: Vec3::new(0.25, plane_speed, -0.1),
            ..plane_contact(Vec3::default(), Vec3::default(), 0.6)
        };
        let result = resolve_articulated_plane_impact(&model, base, &[], &[], contact).unwrap();
        let relative = body_speed - plane_speed;
        let impulse = -1.6 * relative * 2.0;
        let endpoint_speed = body_speed + impulse / 2.0;
        let dissipation = (1.0 - 0.6 * 0.6) * relative * relative;
        let kinetic_change = endpoint_speed * endpoint_speed - body_speed * body_speed;
        close(result.normal_law.normal_impulse, impulse);
        close(result.plane_work, impulse * plane_speed);
        close(result.normal_law.dissipated_energy, dissipation);
        close(result.relative_normal_speed_after, -0.6 * relative);
        close(
            result.response.base_after.twist_body.linear.y,
            endpoint_speed,
        );
        close(
            result.response.momentum_after.kinetic_energy
                - result.response.momentum_before.kinetic_energy,
            kinetic_change,
        );
        close(kinetic_change, result.plane_work - dissipation);
        close(result.energy_balance_residual, 0.0);
        if plane_speed > 0.0 {
            assert!(
                kinetic_change > 0.0,
                "moving source must supply actual body energy"
            );
        } else {
            assert!(
                result.plane_work < 0.0,
                "receding source must absorb signed work"
            );
        }
    }
}

/// G0/G4: an instantaneous normal impact neither acts across an open gap nor
/// repairs penetration. Active separating/resting contact is an exact no-op.
#[test]
fn open_and_separating_contacts_are_inert_and_penetration_refuses() {
    let model = single_body(inertia(2.0, Vec3::default(), [1.0; 3]));
    for (gap, speed, active) in [(1e-3, -1.0, false), (0.0, 0.5, true), (0.0, 0.0, true)] {
        let base = FreeFloatingBaseState::new(
            Se3::identity(),
            Twist::new(Vec3::default(), Vec3::new(0.0, speed, 0.0)),
        );
        let contact = plane_contact(Vec3::default(), Vec3::new(0.0, -gap, 0.0), 0.5);
        let result = resolve_articulated_plane_impact(&model, base, &[], &[], contact).unwrap();
        close(result.gap, gap);
        assert_eq!(result.normal_law.active, active);
        assert!(!result.normal_law.impacted);
        assert_eq!(result.response.base_after, base);
        assert_eq!(
            result.response.momentum_after,
            result.response.momentum_before
        );
        close(result.normal_law.normal_impulse, 0.0);
        close(result.normal_law.dissipated_energy, 0.0);
        close(result.plane_work, 0.0);
    }
    let base = FreeFloatingBaseState::stationary(Se3::identity());
    let penetrating = plane_contact(Vec3::default(), Vec3::new(0.0, 1e-12, 0.0), 0.5);
    assert!(matches!(
        resolve_articulated_plane_impact(&model, base, &[], &[], penetrating),
        Err(ArticulatedImpactError::ContactLaw(_))
    ));
    let invalid = ArticulatedPointPlaneContact {
        normal_world: Vec3::default(),
        ..penetrating
    };
    assert!(matches!(
        resolve_articulated_plane_impact(&model, base, &[], &[], invalid),
        Err(ArticulatedImpactError::InvalidNormal)
    ));
}

/// G4: unit impulse probes are an inertia calculation, not physical endpoint
/// states. A light link must admit a small impact even when a unit probe would
/// exceed its speed or actuator-effort bounds by many orders of magnitude.
#[test]
fn effective_mass_probe_does_not_enforce_fictitious_endpoint_limits() {
    let model = two_link_model(1e-4, 0.25);
    let base = FreeFloatingBaseState::new(
        Se3::identity(),
        Twist::new(Vec3::default(), Vec3::new(0.0, -1e-4, 0.0)),
    );
    let point = Vec3::new(1.0, 0.0, 0.0);
    let contact = ArticulatedPointPlaneContact {
        link_index: 1,
        ..plane_contact(point, point, 0.0)
    };
    let result = resolve_articulated_plane_impact(&model, base, &[0.0], &[0.0], contact).unwrap();
    close(result.inverse_effective_mass, 10_000.2);
    close(result.normal_law.normal_impulse, 1e-4 / 10_000.2);
    close(
        result.response.joint_velocity_after[0],
        (1e-4 / 10_000.2) / 1e-4,
    );
    close(result.relative_normal_speed_after, 0.0);
    assert!(result.response.joint_velocity_after[0].abs() < 0.25);
}

/// G4: actual post-impact joint limits still apply. Refusal cannot clip the
/// impulse, modify the configuration, or publish a partially updated base.
#[test]
fn inadmissible_actual_joint_speed_refuses_without_partial_state() {
    let model = two_link_model(2.0, 0.25);
    let base = FreeFloatingBaseState::new(
        Se3::identity(),
        Twist::new(Vec3::default(), Vec3::new(0.0, -0.6, 0.0)),
    );
    let q = [0.0];
    let velocity = [-0.1];
    let before = free_floating_momentum(&model, base, &q, &velocity).unwrap();
    let point = Vec3::new(1.0, 0.0, 0.0);
    let contact = ArticulatedPointPlaneContact {
        link_index: 1,
        ..plane_contact(point, point, 0.4)
    };
    assert!(matches!(
        resolve_articulated_plane_impact(&model, base, &q, &velocity, contact),
        Err(ArticulatedImpactError::Dynamics(
            ArticulatedError::JointVelocityLimitViolation { .. }
        ))
    ));
    assert_eq!(
        free_floating_momentum(&model, base, &q, &velocity).unwrap(),
        before
    );
    assert_eq!(base.world_from_base, Se3::identity());
    assert_eq!(base.twist_body.linear, Vec3::new(0.0, -0.6, 0.0));
    assert_eq!(q, [0.0]);
    assert_eq!(velocity, [-0.1]);
}
