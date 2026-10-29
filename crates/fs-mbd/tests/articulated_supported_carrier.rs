//! G1/G3/G4: automatically discovered supported motion retains native recoil.
use fs_contact::frictional_multicontact::FrictionalContactConfig;
use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::{FreeFloatingMotionLoads, FreeFloatingMotionState};
use fs_mbd::articulated::{
    ArticulatedModel, FreeFloatingBaseState, JointModel, Link, SpatialInertia,
};
use fs_mbd::articulated_contacts::step::ArticulatedContactStepConfig;
use fs_mbd::articulated_contacts::trajectory::ArticulatedTrajectorySystem;
use fs_mbd::articulated_contacts::trajectory::sphere_planes::CarrierPlane;
use fs_mbd::articulated_contacts::trajectory::sphere_planes::support::{
    CarrierSupportConfig, CarrierSupportError, CarrierSupportResult,
    advance_supported_carrier_planes,
};
use std::cell::Cell;

const GAP: f64 = 1e-9;
fn inertia(mass: f64, moment: f64) -> SpatialInertia {
    SpatialInertia::new(
        mass,
        Vec3::default(),
        Mat3 {
            m: [moment, 0.0, 0.0, 0.0, moment, 0.0, 0.0, 0.0, moment],
        },
    )
    .unwrap()
}
fn sphere() -> ArticulatedModel {
    ArticulatedModel::new(vec![Link::new(
        "sphere",
        None,
        Se3::identity(),
        JointModel::FIXED,
        inertia(2.0, 0.2),
    )])
    .unwrap()
}
fn carrier() -> ArticulatedModel {
    ArticulatedModel::new(vec![
        Link::new(
            "shell",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(2.0, 1.0),
        ),
        Link::new(
            "wheel",
            Some(0),
            Se3::identity(),
            JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None).unwrap(),
            inertia(1.0, 0.2),
        ),
    ])
    .unwrap()
}
fn state(position: Vec3, velocity: Vec3, dofs: usize) -> FreeFloatingMotionState {
    FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(
            Se3::from_parts(So3::identity(), position).unwrap(),
            Twist::new(Vec3::default(), velocity),
        ),
        joint_position: vec![0.0; dofs],
        joint_velocity: vec![0.0; dofs],
    }
}
fn plane(mu: f64) -> CarrierPlane {
    CarrierPlane {
        id: 5,
        point_world: Vec3::default(),
        velocity_world: Vec3::default(),
        normal_world: Vec3::new(0.0, 1.0, 0.0),
        tangent_world: Vec3::new(1.0, 0.0, 0.0),
        friction_coefficient: mu,
        restitution: 0.0,
    }
}
fn config(duration: f64) -> CarrierSupportConfig {
    CarrierSupportConfig {
        step: ArticulatedContactStepConfig {
            duration,
            angular_origin_world: Vec3::new(2.0, -1.0, 3.0),
            contact: FrictionalContactConfig::default(),
            drift_energy_tolerance: 1e-8,
            energy_tolerance: 1e-8,
            linear_momentum_tolerance: 1e-8,
            angular_momentum_tolerance: 1e-8,
        },
        maximum_step: 0.01,
        maximum_attempts: 4096,
        activation_tolerance: 1e-7,
        maximum_geometry_queries: 100_000,
    }
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 3e-8, "{a:.16e} != {b:.16e}");
}
fn vector(a: Vec3, b: Vec3) {
    assert!((a - b).norm() < 3e-8, "{a:?} != {b:?}");
}
fn velocity(state: &FreeFloatingMotionState) -> Vec3 {
    state
        .base
        .world_from_base
        .rotation()
        .rotate(state.base.twist_body.linear)
        .unwrap()
}
fn spin(state: &FreeFloatingMotionState) -> Vec3 {
    state
        .base
        .world_from_base
        .rotation()
        .rotate(state.base.twist_body.angular)
        .unwrap()
}
fn audit(result: &CarrierSupportResult) {
    assert!(result.minimum_gap >= 0.0);
    vector(result.linear_momentum_residual, Vec3::default());
    vector(result.angular_momentum_residual, Vec3::default());
    close(result.kinetic_work_residual, result.drift_energy_defect);
    close(
        result.momentum_after.kinetic_energy - result.momentum_before.kinetic_energy,
        result.actuator_impulse_work
            + result.gravity_impulse_work
            + result.body_torque_impulse_work
            + result.support_work
            - result.contact_loss
            + result.kinetic_work_residual,
    );
}

#[test]
fn internal_wheel_actuation_rolls_the_shell_with_analytic_recoil_and_gravity_support() {
    let model = carrier();
    let initial = state(Vec3::new(0.0, 1.0 + GAP, 0.0), Vec3::default(), 1);
    let wrenches = [Wrench::default(); 2];
    let result = advance_supported_carrier_planes(
        ArticulatedTrajectorySystem {
            model: &model,
            state: &initial,
            loads: FreeFloatingMotionLoads {
                generalized_force: &[0.4],
                external_body_wrench: &wrenches,
                gravity_world: Vec3::new(0.0, -9.81, 0.0),
            },
        },
        1.0,
        &[plane(0.3)],
        config(1.0),
        || false,
    )
    .unwrap();
    // Independent constrained mechanics: a=tau*R/(I_shell+M*R²),
    // omega_shell=-v/R and omega_wheel_absolute=tau*t/I_wheel.
    let v = 0.4 / (1.0 + 3.0);
    vector(velocity(&result.state_after), Vec3::new(v, 0.0, 0.0));
    vector(spin(&result.state_after), Vec3::new(0.0, 0.0, -v));
    close(result.state_after.joint_velocity[0], 2.0 + v);
    close(
        result.state_after.base.world_from_base.translation().x,
        0.5 * v * 1.01,
    );
    close(
        result.state_after.base.world_from_base.translation().y,
        1.0 + GAP,
    );
    close(result.state_after.joint_position[0], 0.5 * (2.0 + v) * 1.01);
    close(
        result.actuator_displacement_work,
        0.4 * result.state_after.joint_position[0],
    );
    close(result.gravity_displacement_work, 0.0);
    vector(
        result.external_linear_impulse_world,
        Vec3::new(3.0 * v, 0.0, 0.0),
    );
    assert!(result.contact_loss > 0.0); // The supported gravity kick has splitting loss.
    audit(&result);
}

#[test]
fn kinetic_sliding_transitions_to_unforced_rolling_without_velocity_projection() {
    let model = sphere();
    let initial = state(Vec3::new(0.0, 0.5 + GAP, 0.0), Vec3::new(1.0, 0.0, 0.0), 0);
    let wrenches = [Wrench::default()];
    let result = advance_supported_carrier_planes(
        ArticulatedTrajectorySystem {
            model: &model,
            state: &initial,
            loads: FreeFloatingMotionLoads {
                generalized_force: &[],
                external_body_wrench: &wrenches,
                gravity_world: Vec3::new(0.0, -10.0, 0.0),
            },
        },
        0.5,
        &[plane(0.2)],
        config(0.5),
        || false,
    )
    .unwrap();
    // Slip stops at 1/(mu*g*(1+M*R²/I))=1/7 seconds. Linear and
    // angular impulse balance give v=M*R²/(I+M*R²)=5/7 thereafter.
    close(velocity(&result.state_after).x, 5.0 / 7.0);
    close(spin(&result.state_after).z, -10.0 / 7.0);
    close(
        velocity(&result.state_after).x + 0.5 * spin(&result.state_after).z,
        0.0,
    );
    close(result.support_work, 0.0);
    assert!(result.contact_loss > 2.0 / 7.0);
    audit(&result);
}

#[test]
fn inelastic_landing_discovers_the_plane_and_continues_as_sustained_support() {
    let model = sphere();
    let initial = state(Vec3::new(0.0, 1.5, 0.0), Vec3::default(), 0);
    let saved = initial.clone();
    let wrenches = [Wrench::default()];
    let mut settings = config(2.0);
    settings.maximum_step = 0.05;
    let result = advance_supported_carrier_planes(
        ArticulatedTrajectorySystem {
            model: &model,
            state: &initial,
            loads: FreeFloatingMotionLoads {
                generalized_force: &[],
                external_body_wrench: &wrenches,
                gravity_world: Vec3::new(0.0, -2.0, 0.0),
            },
        },
        0.5,
        &[plane(0.0)],
        settings,
        || false,
    )
    .unwrap();
    let height = result.state_after.base.world_from_base.translation().y;
    assert!((0.5..=0.5 + settings.activation_tolerance).contains(&height));
    vector(velocity(&result.state_after), Vec3::default());
    close(result.gravity_displacement_work, 4.0 * (1.5 - height));
    assert!(result.attempted_steps > result.accepted_steps);
    assert!(result.accepted_steps >= 40);
    assert_eq!(initial, saved);
    audit(&result);
}

#[test]
fn inactive_wall_is_discovered_during_already_supported_sliding() {
    let model = sphere();
    let initial = state(Vec3::new(-1.5, 0.5 + GAP, 0.0), Vec3::new(2.0, 0.0, 0.0), 0);
    let wrenches = [Wrench::default()];
    let wall = CarrierPlane {
        id: 9,
        normal_world: Vec3::new(-1.0, 0.0, 0.0),
        tangent_world: Vec3::new(0.0, 0.0, 1.0),
        ..plane(0.0)
    };
    let mut settings = config(1.0);
    settings.maximum_step = 0.08;
    let result = advance_supported_carrier_planes(
        ArticulatedTrajectorySystem {
            model: &model,
            state: &initial,
            loads: FreeFloatingMotionLoads {
                generalized_force: &[],
                external_body_wrench: &wrenches,
                gravity_world: Vec3::new(0.0, -10.0, 0.0),
            },
        },
        0.5,
        &[plane(0.0), wall],
        settings,
        || false,
    )
    .unwrap();
    let position = result.state_after.base.world_from_base.translation();
    assert!((-0.5 - settings.activation_tolerance..=-0.5).contains(&position.x));
    close(position.y, 0.5 + GAP);
    vector(velocity(&result.state_after), Vec3::default());
    vector(
        result.external_linear_impulse_world,
        Vec3::new(-4.0, 0.0, 0.0),
    );
    assert!(result.attempted_steps > result.accepted_steps);
    assert!(result.contact_loss > 4.0);
    audit(&result);
}

#[test]
fn translating_support_supplies_original_point_work_once() {
    let model = sphere();
    let initial = state(Vec3::new(0.0, 0.5 + GAP, 0.0), Vec3::new(0.0, 0.2, 0.0), 0);
    let wrenches = [Wrench::default()];
    let moving = CarrierPlane {
        velocity_world: Vec3::new(1.0, 0.2, 0.0),
        ..plane(0.2)
    };
    let result = advance_supported_carrier_planes(
        ArticulatedTrajectorySystem {
            model: &model,
            state: &initial,
            loads: FreeFloatingMotionLoads {
                generalized_force: &[],
                external_body_wrench: &wrenches,
                gravity_world: Vec3::new(0.0, -10.0, 0.0),
            },
        },
        0.5,
        &[moving],
        config(0.5),
        || false,
    )
    .unwrap();
    close(velocity(&result.state_after).x, 2.0 / 7.0);
    close(spin(&result.state_after).z, 10.0 / 7.0);
    close(
        result.state_after.base.world_from_base.translation().y,
        0.6 + GAP,
    );
    close(result.support_work, 4.0 / 7.0 + 2.0);
    close(result.gravity_displacement_work, -2.0);
    audit(&result);
}

#[test]
fn simultaneous_corner_support_preserves_plane_order_and_world_coordinates() {
    let model = sphere();
    let initial = state(Vec3::new(0.5 + GAP, 0.5 + GAP, 0.0), Vec3::default(), 0);
    let wrenches = [Wrench::default()];
    let gravity = Vec3::new(-2.0, -3.0, 0.0);
    let planes = [
        plane(0.0),
        CarrierPlane {
            id: 2,
            normal_world: Vec3::new(1.0, 0.0, 0.0),
            tangent_world: Vec3::new(0.0, 0.0, 1.0),
            ..plane(0.0)
        },
    ];
    let run = |state: &FreeFloatingMotionState, planes: &[CarrierPlane], gravity| {
        advance_supported_carrier_planes(
            ArticulatedTrajectorySystem {
                model: &model,
                state,
                loads: FreeFloatingMotionLoads {
                    generalized_force: &[],
                    external_body_wrench: &wrenches,
                    gravity_world: gravity,
                },
            },
            0.5,
            planes,
            config(0.2),
            || false,
        )
        .unwrap()
    };
    let original = run(&initial, &planes, gravity);
    let swapped = run(&initial, &[planes[1], planes[0]], gravity);
    assert_eq!(original.state_after, swapped.state_after);
    assert_eq!(
        original.external_angular_impulse_about_origin_world,
        swapped.external_angular_impulse_about_origin_world
    );
    vector(velocity(&original.state_after), Vec3::default());
    let q = So3::exp(So3Tangent::new(Vec3::new(0.3, -0.2, 0.5))).unwrap();
    let shift = Vec3::new(1.0, -0.5, 0.7);
    let rotated = FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(
            Se3::from_parts(
                q,
                q.rotate(initial.base.world_from_base.translation())
                    .unwrap()
                    + shift,
            )
            .unwrap(),
            Twist::zero(),
        ),
        ..initial.clone()
    };
    let transformed = planes.map(|p| CarrierPlane {
        point_world: shift,
        normal_world: q.rotate(p.normal_world).unwrap(),
        tangent_world: q.rotate(p.tangent_world).unwrap(),
        ..p
    });
    let result = run(&rotated, &transformed, q.rotate(gravity).unwrap());
    vector(
        result.state_after.base.world_from_base.translation(),
        q.rotate(original.state_after.base.world_from_base.translation())
            .unwrap()
            + shift,
    );
    vector(velocity(&result.state_after), Vec3::default());
    close(result.contact_loss, original.contact_loss);
    audit(&result);
}

#[test]
fn unsupported_models_late_budgets_and_cancellation_publish_no_partial_carrier() {
    let model = sphere();
    let initial = state(Vec3::new(0.0, 0.5 + GAP, 0.0), Vec3::default(), 0);
    let saved = initial.clone();
    let wrenches = [Wrench::default()];
    let input = ArticulatedTrajectorySystem {
        model: &model,
        state: &initial,
        loads: FreeFloatingMotionLoads {
            generalized_force: &[],
            external_body_wrench: &wrenches,
            gravity_world: Vec3::new(0.0, -10.0, 0.0),
        },
    };
    let settings = config(0.2);
    let mut exhausted = settings;
    exhausted.maximum_geometry_queries = 7;
    assert!(
        advance_supported_carrier_planes(input, 0.5, &[plane(0.0)], exhausted, || false).is_err()
    );
    let count = Cell::new(0usize);
    assert!(
        advance_supported_carrier_planes(input, 0.5, &[plane(0.0)], settings, || {
            count.set(count.get() + 1);
            count.get() > 150
        })
        .is_err()
    );
    assert!(count.get() > 150);
    assert!(matches!(
        advance_supported_carrier_planes(
            input,
            0.5,
            &[CarrierPlane {
                restitution: 0.2,
                ..plane(0.0)
            }],
            settings,
            || false
        ),
        Err(CarrierSupportError::Restitution { plane_id: 5 })
    ));
    let invalid = [Wrench::new(Vec3::default(), Vec3::new(1.0, 0.0, 0.0))];
    assert!(
        advance_supported_carrier_planes(
            ArticulatedTrajectorySystem {
                loads: FreeFloatingMotionLoads {
                    external_body_wrench: &invalid,
                    ..input.loads
                },
                ..input
            },
            0.5,
            &[plane(0.0)],
            settings,
            || false
        )
        .is_err()
    );
    assert_eq!(initial, saved);
    let first =
        advance_supported_carrier_planes(input, 0.5, &[plane(0.0)], settings, || false).unwrap();
    let retry =
        advance_supported_carrier_planes(input, 0.5, &[plane(0.0)], settings, || false).unwrap();
    assert_eq!(first.state_after, retry.state_after);
    assert_eq!(first.contact_loss, retry.contact_loss);
}
