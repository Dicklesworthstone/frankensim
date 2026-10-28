//! Analytic sphere/plane trajectories through detection, impact and support.

use fs_mbd::sphere_plane::{
    SpherePlane, SpherePlaneAdvance, SpherePlaneConfig, SpherePlaneError, SpherePlaneEventKind,
    SpherePlaneMode, SpherePlaneState, advance_sphere_plane,
};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};

const MASS: f64 = 2.0;
const RADIUS: f64 = 0.5;
const INERTIA: f64 = 0.2;

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 5e-10,
        "{actual:e} != {expected:e}"
    );
}

fn close_vector(actual: Vec3, expected: Vec3) {
    assert!(
        actual.sub(expected).norm_squared().sqrt() < 5e-10,
        "{actual:?} != {expected:?}"
    );
}

fn scene() -> SpherePlane {
    SpherePlane {
        properties: MassProperties::new(MASS, Vec3::ZERO, Vec3::new(INERTIA, INERTIA, INERTIA))
            .unwrap(),
        radius: RADIUS,
        point_world: Vec3::ZERO,
        normal_world: Vec3::new(0.0, 1.0, 0.0),
        velocity_world: Vec3::ZERO,
        force_world: Vec3::new(0.0, -20.0, 0.0),
        restitution: 0.5,
        friction_coefficient: 0.2,
    }
}

fn state(
    position: Vec3,
    velocity: Vec3,
    spin_world: Vec3,
    orientation: UnitQuaternion,
) -> SpherePlaneState {
    let body = RigidBodyState::new(
        Pose::new(position, orientation).unwrap(),
        velocity.scale(MASS),
        orientation.rotate_world_to_body(spin_world.scale(INERTIA)),
    )
    .unwrap();
    SpherePlaneState::new(body, 0.0).unwrap()
}

fn velocity(body: RigidBodyState) -> Vec3 {
    body.linear_momentum_world().scale(1.0 / MASS)
}
fn spin(body: RigidBodyState) -> Vec3 {
    body.pose()
        .orientation()
        .rotate_body_to_world(body.angular_momentum_body())
        .scale(1.0 / INERTIA)
}
fn kinetic(body: RigidBodyState) -> f64 {
    body.linear_momentum_world().norm_squared() / (2.0 * MASS)
        + body.angular_momentum_body().norm_squared() / (2.0 * INERTIA)
}
fn angular_momentum(body: RigidBodyState) -> Vec3 {
    body.pose()
        .orientation()
        .rotate_body_to_world(body.angular_momentum_body())
        .add(
            body.pose()
                .position_world()
                .cross(body.linear_momentum_world()),
        )
}

fn balances(
    scene: &SpherePlane,
    before: SpherePlaneState,
    duration: f64,
    result: &SpherePlaneAdvance,
) {
    let after = result.state_after.body();
    close(result.state_after.time_s(), before.time_s() + duration);
    close_vector(
        after.linear_momentum_world(),
        before
            .body()
            .linear_momentum_world()
            .add(scene.force_world.scale(duration))
            .add(result.contact_impulse_world),
    );
    close(
        kinetic(after) - kinetic(before.body()),
        result.external_force_work + result.plane_work - result.dissipated_energy,
    );
    close(result.energy_balance_residual, 0.0);
    close_vector(result.linear_momentum_residual, Vec3::ZERO);
    close_vector(result.angular_momentum_residual, Vec3::ZERO);
    let plane = scene
        .point_world
        .add(scene.velocity_world.scale(result.state_after.time_s()));
    let gap = scene
        .normal_world
        .dot(after.pose().position_world().sub(plane))
        - RADIUS;
    close(result.gap_after, gap);
    assert!(
        gap >= -2e-12,
        "accepted sphere penetrates plane by {:e}",
        -gap
    );
    assert!(result.dissipated_energy >= 0.0);
}

/// G1/G3: the incoming parabola hits at .2s, and the remaining .1s must be
/// integrated after restitution. A normal frame boost supplies exactly J*U
/// through the moving plane without changing the impact loss or TOI.
#[test]
fn ballistic_impact_uses_exact_time_and_moving_plane_normal_work() {
    for boost in [0.0, 0.5] {
        let mut scene = scene();
        scene.velocity_world = Vec3::new(0.0, boost, 0.0);
        let initial = state(
            Vec3::new(0.0, 0.9, 0.0),
            Vec3::new(0.0, -1.0 + boost, 0.0),
            Vec3::ZERO,
            UnitQuaternion::IDENTITY,
        );
        let result =
            advance_sphere_plane(&scene, initial, 0.3, SpherePlaneConfig::default(), || false)
                .unwrap();
        let impact = result
            .events
            .iter()
            .find(|event| event.kind == SpherePlaneEventKind::Impact)
            .unwrap();
        close(impact.time_s, 0.2);
        close(impact.normal_impulse, 9.0);
        close_vector(
            result.state_after.body().pose().position_world(),
            Vec3::new(0.0, 0.6 + 0.3 * boost, 0.0),
        );
        close_vector(
            velocity(result.state_after.body()),
            Vec3::new(0.0, 0.5 + boost, 0.0),
        );
        close(result.gap_after, 0.1);
        close(result.plane_work, 9.0 * boost);
        close(result.external_force_work, 6.0 - 6.0 * boost);
        close(result.dissipated_energy, 6.75);
        assert_eq!(result.state_after.mode(), SpherePlaneMode::Flight);
        balances(&scene, initial, 0.3, &result);
    }
}

/// G1/G3: both unconstrained endpoints are above the plane, but the parabola
/// passes through it. Continuous detection must find the first crossing. A
/// departing sphere under inward force must also find its later return.
#[test]
fn ballistic_detection_catches_hidden_crossing_and_departing_return() {
    let mut scene = scene();
    scene.force_world = Vec3::new(0.0, 8.0, 0.0);
    let initial = state(
        Vec3::new(0.0, 1.5, 0.0),
        Vec3::new(0.0, -4.0, 0.0),
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
    );
    let result =
        advance_sphere_plane(&scene, initial, 2.0, SpherePlaneConfig::default(), || false).unwrap();
    let root_two = 2.0_f64.sqrt();
    let impact = result
        .events
        .iter()
        .find(|event| event.kind == SpherePlaneEventKind::Impact)
        .unwrap();
    close(impact.time_s, 1.0 - 1.0 / root_two);
    close(impact.normal_impulse, 6.0 * root_two);
    close(result.gap_after, 4.0 + 3.0 * root_two);
    close(velocity(result.state_after.body()).y, 4.0 + 3.0 * root_two);
    close(result.dissipated_energy, 6.0);
    balances(&scene, initial, 2.0, &result);

    scene.force_world = Vec3::new(0.0, -20.0, 0.0);
    let departing = state(
        Vec3::new(0.0, RADIUS, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
    );
    let returned = advance_sphere_plane(
        &scene,
        departing,
        0.25,
        SpherePlaneConfig::default(),
        || false,
    )
    .unwrap();
    let impact = returned
        .events
        .iter()
        .find(|event| event.kind == SpherePlaneEventKind::Impact)
        .unwrap();
    close(impact.time_s, 0.2);
    close(impact.normal_impulse, 3.0);
    close(returned.gap_after, 0.0125);
    close(velocity(returned.state_after.body()).y, 0.0);
    balances(&scene, departing, 0.25, &returned);
}

/// G0/G1: persistent support balances gravity without repeated impacts or
/// numerical creep. A translating normal support does mechanical work; an
/// outward force subsequently releases supported contact without adhesion.
#[test]
fn persistent_support_and_liftoff_obey_normal_force_and_work() {
    for normal_speed in [0.0, 0.4] {
        let mut scene = scene();
        scene.velocity_world = Vec3::new(0.0, normal_speed, 0.0);
        let initial = state(
            Vec3::new(0.0, RADIUS, 0.0),
            scene.velocity_world,
            Vec3::ZERO,
            UnitQuaternion::IDENTITY,
        );
        let supported =
            advance_sphere_plane(&scene, initial, 0.5, SpherePlaneConfig::default(), || false)
                .unwrap();
        close(
            supported.state_after.body().pose().position_world().y,
            RADIUS + 0.5 * normal_speed,
        );
        close_vector(velocity(supported.state_after.body()), scene.velocity_world);
        close_vector(supported.contact_impulse_world, Vec3::new(0.0, 10.0, 0.0));
        close(supported.plane_work, 10.0 * normal_speed);
        close(supported.external_force_work, -10.0 * normal_speed);
        close(supported.dissipated_energy, 0.0);
        assert_eq!(supported.state_after.mode(), SpherePlaneMode::Rolling);
        assert!(
            !supported
                .events
                .iter()
                .any(|event| event.kind == SpherePlaneEventKind::Impact)
        );
        balances(&scene, initial, 0.5, &supported);

        scene.force_world = Vec3::new(0.0, 4.0, 0.0);
        let released = advance_sphere_plane(
            &scene,
            supported.state_after,
            0.5,
            SpherePlaneConfig::default(),
            || false,
        )
        .unwrap();
        close(released.gap_after, 0.25);
        close(velocity(released.state_after.body()).y, normal_speed + 1.0);
        close_vector(released.contact_impulse_world, Vec3::ZERO);
        assert_eq!(released.state_after.mode(), SpherePlaneMode::Flight);
        assert!(
            released
                .events
                .iter()
                .any(|event| event.kind == SpherePlaneEventKind::LiftOff)
        );
        balances(&scene, supported.state_after, 0.5, &released);
    }
}

/// G1: m=2,R=.5,I=.2,mu=.2,g=10 gives t_roll=.2s, v_roll=1,
/// omega_z=-2 and x_roll=.24. Continuing to .5s gives x=.54 and loss=.56J.
/// Splitting the caller's interval must retain the supported slip history.
#[test]
fn sliding_sphere_reaches_rolling_at_the_analytic_time() {
    let scene = scene();
    let initial = state(
        Vec3::new(0.0, RADIUS, 0.0),
        Vec3::new(1.4, 0.0, 0.0),
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
    );
    // Collinear angular momentum and torque commute, so analytic attitude
    // must not spend the bounded substep budget, even for a tiny step cap.
    let result = advance_sphere_plane(
        &scene,
        initial,
        0.5,
        SpherePlaneConfig {
            max_attitude_step_s: 1e-12,
            max_attitude_substeps: 1,
            ..SpherePlaneConfig::default()
        },
        || false,
    )
    .unwrap();
    let rolling = result
        .events
        .iter()
        .find(|event| event.kind == SpherePlaneEventKind::Rolling)
        .unwrap();
    close(rolling.time_s, 0.2);
    close_vector(
        result.state_after.body().pose().position_world(),
        Vec3::new(0.54, RADIUS, 0.0),
    );
    close_vector(
        velocity(result.state_after.body()),
        Vec3::new(1.0, 0.0, 0.0),
    );
    close_vector(spin(result.state_after.body()), Vec3::new(0.0, 0.0, -2.0));
    // Integral of omega_z: -10 * (0.2^2) / 2 - 2 * 0.3 = -0.8 rad.
    let (sin_angle, cos_angle) = (-0.8_f64).sin_cos();
    let orientation = result.state_after.body().pose().orientation();
    close_vector(
        orientation.rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0)),
        Vec3::new(cos_angle, sin_angle, 0.0),
    );
    close_vector(
        orientation.rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0)),
        Vec3::new(0.0, 0.0, 1.0),
    );
    close(result.dissipated_energy, 0.56);
    close(kinetic(result.state_after.body()), 1.4);
    close_vector(
        angular_momentum(result.state_after.body()),
        angular_momentum(initial.body()),
    );
    assert_eq!(result.state_after.mode(), SpherePlaneMode::Rolling);
    balances(&scene, initial, 0.5, &result);

    let first =
        advance_sphere_plane(&scene, initial, 0.1, SpherePlaneConfig::default(), || false).unwrap();
    assert_eq!(first.state_after.mode(), SpherePlaneMode::Sliding);
    let second = advance_sphere_plane(
        &scene,
        first.state_after,
        0.4,
        SpherePlaneConfig::default(),
        || false,
    )
    .unwrap();
    close_vector(
        second.state_after.body().pose().position_world(),
        result.state_after.body().pose().position_world(),
    );
    close_vector(
        velocity(second.state_after.body()),
        velocity(result.state_after.body()),
    );
    close_vector(
        spin(second.state_after.body()),
        spin(result.state_after.body()),
    );
    close(
        first.dissipated_energy + second.dissipated_energy,
        result.dissipated_energy,
    );
    balances(&scene, first.state_after, 0.4, &second);
}

/// G1/G3: oblique spin creates slip despite zero tangential COM velocity.
/// Both tangent components must share one Coulomb disk, normal spin is
/// unaffected, and the impact at the world origin preserves total angular
/// momentum. Low and high friction independently exercise sliding/sticking.
#[test]
fn impact_friction_resolves_contact_spin_and_the_two_dimensional_cone() {
    for mu in [0.05, 0.2] {
        let mut scene = scene();
        scene.force_world = Vec3::ZERO;
        scene.restitution = 0.0;
        scene.friction_coefficient = mu;
        let initial = state(
            Vec3::new(0.0, RADIUS, 0.0),
            Vec3::new(0.0, -2.0, 0.0),
            Vec3::new(1.2, 3.0, 1.6),
            UnitQuaternion::IDENTITY,
        );
        let result =
            advance_sphere_plane(&scene, initial, 0.1, SpherePlaneConfig::default(), || false)
                .unwrap();
        let impulse = if mu < 0.1 { 0.2 } else { 4.0 / 7.0 };
        let tangent = Vec3::new(-0.8 * impulse, 0.0, 0.6 * impulse);
        let expected_velocity = tangent.scale(1.0 / MASS);
        let expected_spin = Vec3::new(1.2, 3.0, 1.6).add(
            Vec3::new(0.0, -RADIUS, 0.0)
                .cross(tangent)
                .scale(1.0 / INERTIA),
        );
        close_vector(
            result.contact_impulse_world,
            tangent.add(Vec3::new(0.0, 4.0, 0.0)),
        );
        close_vector(velocity(result.state_after.body()), expected_velocity);
        close_vector(spin(result.state_after.body()), expected_spin);
        close_vector(
            result.state_after.body().pose().position_world(),
            Vec3::new(0.0, RADIUS, 0.0).add(expected_velocity.scale(0.1)),
        );
        close(
            result.dissipated_energy,
            4.0 + if mu < 0.1 { 0.165 } else { 2.0 / 7.0 },
        );
        let slip = expected_velocity.add(expected_spin.cross(Vec3::new(0.0, -RADIUS, 0.0)));
        close(
            slip.norm_squared().sqrt(),
            if mu < 0.1 { 0.65 } else { 0.0 },
        );
        close_vector(
            angular_momentum(result.state_after.body()),
            angular_momentum(initial.body()),
        );
        close_vector(result.contact_angular_impulse_world, Vec3::ZERO);
        balances(&scene, initial, 0.1, &result);
    }
}

/// G1/G3: a conveyor supplies work to both translation and spin. A common
/// velocity boost and an exact cyclic world rotation preserve relative slip
/// and dissipation while transforming positions, momenta and source work.
#[test]
fn moving_conveyor_work_and_frame_covariance_match_analytic_response() {
    let mut scene = scene();
    scene.velocity_world = Vec3::new(1.4, 0.0, 0.0);
    let initial = state(
        Vec3::new(0.0, RADIUS, 0.0),
        Vec3::ZERO,
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
    );
    let result =
        advance_sphere_plane(&scene, initial, 0.5, SpherePlaneConfig::default(), || false).unwrap();
    close_vector(
        result.state_after.body().pose().position_world(),
        Vec3::new(0.16, RADIUS, 0.0),
    );
    close_vector(
        velocity(result.state_after.body()),
        Vec3::new(0.4, 0.0, 0.0),
    );
    close_vector(spin(result.state_after.body()), Vec3::new(0.0, 0.0, 2.0));
    close(result.plane_work, 1.12);
    close(result.dissipated_energy, 0.56);
    close(kinetic(result.state_after.body()), 0.56);
    balances(&scene, initial, 0.5, &result);

    let rotation = UnitQuaternion::new(0.5, 0.5, 0.5, 0.5).unwrap();
    let rotate = |v| rotation.rotate_body_to_world(v);
    let shift = Vec3::new(3.0, -4.0, 2.0);
    let boost = Vec3::new(0.3, -0.2, 0.1);
    let transformed_scene = SpherePlane {
        point_world: shift,
        normal_world: rotate(scene.normal_world),
        force_world: rotate(scene.force_world),
        velocity_world: rotate(scene.velocity_world.add(boost)),
        ..scene
    };
    let transformed_initial = state(
        rotate(initial.body().pose().position_world()).add(shift),
        rotate(boost),
        Vec3::ZERO,
        rotation,
    );
    let transformed = advance_sphere_plane(
        &transformed_scene,
        transformed_initial,
        0.5,
        SpherePlaneConfig::default(),
        || false,
    )
    .unwrap();
    close_vector(
        transformed.state_after.body().pose().position_world(),
        rotate(
            result
                .state_after
                .body()
                .pose()
                .position_world()
                .add(boost.scale(0.5)),
        )
        .add(shift),
    );
    close_vector(
        velocity(transformed.state_after.body()),
        rotate(velocity(result.state_after.body()).add(boost)),
    );
    close_vector(
        spin(transformed.state_after.body()),
        rotate(spin(result.state_after.body())),
    );
    close_vector(
        transformed.contact_impulse_world,
        rotate(result.contact_impulse_world),
    );
    close(transformed.dissipated_energy, result.dissipated_energy);
    close(
        transformed.plane_work,
        result.plane_work + boost.dot(result.contact_impulse_world),
    );
    balances(&transformed_scene, transformed_initial, 0.5, &transformed);
}

/// G2: noncommuting normal spin and friction torque require an attitude
/// integration, even though spherical momentum and center motion are exact.
/// Halving the declared Lie-midpoint interval exhibits second-order attitude
/// convergence while leaving the physical momenta and work unchanged.
#[test]
fn noncommuting_spin_and_friction_torque_refine_at_second_order() {
    let scene = scene();
    let initial = state(
        Vec3::new(0.0, RADIUS, 0.0),
        Vec3::new(1.4, 0.0, 0.0),
        Vec3::new(0.0, 3.0, 0.0),
        UnitQuaternion::IDENTITY,
    );
    let mut orientations = Vec::new();
    for h in [0.05, 0.025, 0.0125] {
        let config = SpherePlaneConfig {
            max_attitude_step_s: h,
            ..SpherePlaneConfig::default()
        };
        let result = advance_sphere_plane(&scene, initial, 0.15, config, || false).unwrap();
        close_vector(
            velocity(result.state_after.body()),
            Vec3::new(1.1, 0.0, 0.0),
        );
        close_vector(spin(result.state_after.body()), Vec3::new(0.0, 3.0, -1.5));
        close(result.state_after.body().pose().position_world().x, 0.1875);
        close(result.dissipated_energy, 0.525);
        balances(&scene, initial, 0.15, &result);
        orientations.push(result.state_after.body().pose().orientation());
    }
    let distance = |a: UnitQuaternion, b: UnitQuaternion| {
        [
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ]
        .into_iter()
        .map(|v| {
            a.rotate_body_to_world(v)
                .sub(b.rotate_body_to_world(v))
                .norm_squared()
        })
        .sum::<f64>()
        .sqrt()
    };
    let coarse = distance(orientations[0], orientations[1]);
    let fine = distance(orientations[1], orientations[2]);
    assert!(
        fine > 1e-9,
        "noncommuting torque must exercise attitude integration"
    );
    assert!(
        (3.7..4.3).contains(&(coarse / fine)),
        "midpoint refinement ratio {}",
        coarse / fine
    );

    // A common length-unit rescaling leaves omega and angular acceleration
    // unchanged. The raw L x torque underflows here, but the directions do
    // not commute and must retain the same midpoint attitude trajectory.
    let scale = 1e-100;
    let inertia_scale = scale * scale;
    let scaled_scene = SpherePlane {
        properties: MassProperties::new(
            MASS,
            Vec3::ZERO,
            Vec3::new(
                INERTIA * inertia_scale,
                INERTIA * inertia_scale,
                INERTIA * inertia_scale,
            ),
        )
        .unwrap(),
        radius: scene.radius * scale,
        point_world: scene.point_world.scale(scale),
        velocity_world: scene.velocity_world.scale(scale),
        force_world: scene.force_world.scale(scale),
        ..scene
    };
    let scaled_body = RigidBodyState::new(
        Pose::new(
            initial.body().pose().position_world().scale(scale),
            initial.body().pose().orientation(),
        )
        .unwrap(),
        initial.body().linear_momentum_world().scale(scale),
        initial.body().angular_momentum_body().scale(inertia_scale),
    )
    .unwrap();
    let scaled_initial = SpherePlaneState::new(scaled_body, initial.time_s()).unwrap();
    let scaled_result = advance_sphere_plane(
        &scaled_scene,
        scaled_initial,
        0.15,
        SpherePlaneConfig {
            max_attitude_step_s: 0.025,
            ..SpherePlaneConfig::default()
        },
        || false,
    )
    .unwrap();
    let scaled_error = distance(
        scaled_result.state_after.body().pose().orientation(),
        orientations[1],
    );
    assert!(
        scaled_error < 5e-10,
        "length-unit rescaling changed attitude by {scaled_error:e}"
    );
}

/// G4: unsupported inertia/forcing, initial penetration, cancellation and
/// exhausted event or attitude budgets return no partially advanced state.
#[test]
fn unsupported_physics_and_exhausted_event_budget_refuse_atomically() {
    let scene = scene();
    let initial = state(
        Vec3::new(0.0, 0.9, 0.0),
        Vec3::new(0.0, -1.0, 0.0),
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
    );
    let before = initial.body();
    let limited = SpherePlaneConfig {
        max_events: 1,
        ..SpherePlaneConfig::default()
    };
    assert!(matches!(
        advance_sphere_plane(&scene, initial, 0.7, limited, || false),
        Err(SpherePlaneError::EventBudget)
    ));
    assert!(matches!(
        advance_sphere_plane(&scene, initial, 0.3, SpherePlaneConfig::default(), || true),
        Err(SpherePlaneError::Cancelled)
    ));
    let noncollinear = SpherePlane {
        force_world: Vec3::new(0.0, -20.0, 2.0),
        ..scene
    };
    let sliding = state(
        Vec3::new(0.0, RADIUS, 0.0),
        Vec3::new(1.4, 0.0, 0.0),
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
    );
    assert!(matches!(
        advance_sphere_plane(
            &noncollinear,
            sliding,
            0.3,
            SpherePlaneConfig::default(),
            || false,
        ),
        Err(SpherePlaneError::NonCollinearSliding)
    ));
    let anisotropic = SpherePlane {
        properties: MassProperties::new(MASS, Vec3::ZERO, Vec3::new(0.2, 0.21, 0.22)).unwrap(),
        ..scene
    };
    assert!(
        advance_sphere_plane(
            &anisotropic,
            initial,
            0.3,
            SpherePlaneConfig::default(),
            || false
        )
        .is_err()
    );
    let penetrating = state(
        Vec3::new(0.0, RADIUS - 1e-10, 0.0),
        Vec3::ZERO,
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
    );
    assert!(matches!(
        advance_sphere_plane(
            &scene,
            penetrating,
            0.3,
            SpherePlaneConfig::default(),
            || false
        ),
        Err(SpherePlaneError::Penetration { .. })
    ));
    let attitude_budget = SpherePlaneConfig {
        max_attitude_step_s: 0.01,
        max_attitude_substeps: 1,
        ..SpherePlaneConfig::default()
    };
    let spinning = state(
        Vec3::new(0.0, RADIUS, 0.0),
        Vec3::new(1.4, 0.0, 0.0),
        Vec3::new(0.0, 3.0, 0.0),
        UnitQuaternion::IDENTITY,
    );
    assert!(matches!(
        advance_sphere_plane(&scene, spinning, 0.15, attitude_budget, || false),
        Err(SpherePlaneError::AttitudeBudget)
    ));
    assert_eq!(initial.body(), before);
    close(initial.time_s(), 0.0);
    let accepted =
        advance_sphere_plane(&scene, initial, 0.3, SpherePlaneConfig::default(), || false).unwrap();
    close(accepted.state_after.body().pose().position_world().y, 0.6);
}

#[test]
fn incline_rolling_and_moving_plane_static_traction_obey_exact_work() {
    // Resolve gravity along an incline into tangent and normal coordinates:
    // g_t=6, g_n=-8. A homogeneous sphere accelerates by (5/7)g_t.
    let incline = SpherePlane {
        force_world: Vec3::new(12.0, -16.0, 0.0),
        friction_coefficient: 0.3,
        ..scene()
    };
    let initial = state(
        Vec3::new(0.0, RADIUS, 0.0),
        Vec3::ZERO,
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
    );
    let result = advance_sphere_plane(&incline, initial, 0.7, SpherePlaneConfig::default(), || {
        false
    })
    .unwrap();
    assert_eq!(result.state_after.mode(), SpherePlaneMode::Rolling);
    close_vector(
        result.state_after.body().pose().position_world(),
        Vec3::new(1.05, RADIUS, 0.0),
    );
    close_vector(
        velocity(result.state_after.body()),
        Vec3::new(3.0, 0.0, 0.0),
    );
    close_vector(spin(result.state_after.body()), Vec3::new(0.0, 0.0, -6.0));
    close_vector(result.contact_impulse_world, Vec3::new(-2.4, 11.2, 0.0));
    close(kinetic(result.state_after.body()), 12.6);
    close(result.external_force_work, 12.6);
    close(result.plane_work, 0.0);
    close(result.dissipated_energy, 0.0);
    balances(&incline, initial, 0.7, &result);

    // Static traction can transfer energy through a moving support despite
    // zero slip and zero frictional dissipation.
    let moving = SpherePlane {
        velocity_world: Vec3::new(1.0, 0.0, 0.0),
        force_world: Vec3::new(4.0, -20.0, 0.0),
        ..scene()
    };
    let initial = state(
        Vec3::new(0.0, RADIUS, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
    );
    let result = advance_sphere_plane(&moving, initial, 0.7, SpherePlaneConfig::default(), || {
        false
    })
    .unwrap();
    assert_eq!(result.state_after.mode(), SpherePlaneMode::Rolling);
    close_vector(
        result.state_after.body().pose().position_world(),
        Vec3::new(1.05, RADIUS, 0.0),
    );
    close_vector(
        velocity(result.state_after.body()),
        Vec3::new(2.0, 0.0, 0.0),
    );
    close_vector(spin(result.state_after.body()), Vec3::new(0.0, 0.0, -2.0));
    close_vector(result.contact_impulse_world, Vec3::new(-0.8, 14.0, 0.0));
    close(result.external_force_work, 4.2);
    close(result.plane_work, -0.8);
    close(result.dissipated_energy, 0.0);
    close(kinetic(result.state_after.body()), 4.4);
    balances(&moving, initial, 0.7, &result);
}

#[test]
fn constant_force_sliding_arrests_or_reverses_at_analytic_times() {
    let initial = state(
        Vec3::new(0.0, RADIUS, 0.0),
        Vec3::new(1.4, 0.0, 0.0),
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
    );
    let drive = SpherePlane {
        force_world: Vec3::new(4.0, -20.0, 0.0),
        ..scene()
    };
    let result =
        advance_sphere_plane(&drive, initial, 0.5, SpherePlaneConfig::default(), || false).unwrap();
    // While sliding: a_COM=0, alpha_z=-10, hence slip arrests at .28 s.
    // Then a_COM=10/7 with static traction -8/7.
    let rolling = result
        .events
        .iter()
        .find(|event| event.kind == SpherePlaneEventKind::Rolling)
        .unwrap();
    close(rolling.time_s, 0.28);
    assert_eq!(result.state_after.mode(), SpherePlaneMode::Rolling);
    close(
        result.state_after.body().pose().position_world().x,
        0.7345714285714285,
    );
    close_vector(
        velocity(result.state_after.body()),
        Vec3::new(12.0 / 7.0, 0.0, 0.0),
    );
    close_vector(
        spin(result.state_after.body()),
        Vec3::new(0.0, 0.0, -24.0 / 7.0),
    );
    close(result.dissipated_energy, 0.784);
    close(kinetic(result.state_after.body()), 144.0 / 35.0);
    balances(&drive, initial, 0.5, &result);

    let opposing = SpherePlane {
        force_world: Vec3::new(-16.0, -20.0, 0.0),
        ..scene()
    };
    let result = advance_sphere_plane(
        &opposing,
        initial,
        0.2,
        SpherePlaneConfig::default(),
        || false,
    )
    .unwrap();
    // Initial slip acceleration is -15. At 7/75 s it vanishes, but the
    // static demand 32/7 exceeds the available traction 4. Slip reverses;
    // friction changes sign and its new slip acceleration is -1.
    let reversal = result
        .events
        .iter()
        .find(|event| event.kind == SpherePlaneEventKind::Sliding && event.time_s > 0.0)
        .unwrap();
    close(reversal.time_s, 7.0 / 75.0);
    assert_eq!(result.state_after.mode(), SpherePlaneMode::Sliding);
    close(
        result.state_after.body().pose().position_world().x,
        578.0 / 5625.0,
    );
    close_vector(
        velocity(result.state_after.body()),
        Vec3::new(-13.0 / 75.0, 0.0, 0.0),
    );
    close_vector(
        spin(result.state_after.body()),
        Vec3::new(0.0, 0.0, 2.0 / 15.0),
    );
    close_vector(
        result.contact_impulse_world,
        Vec3::new(4.0 / 75.0, 4.0, 0.0),
    );
    close(result.dissipated_energy, 1598.0 / 5625.0);
    close(kinetic(result.state_after.body()), 179.0 / 5625.0);
    balances(&opposing, initial, 0.2, &result);
}

#[test]
fn arbitrary_force_free_flight_and_zero_friction_support_are_admitted() {
    let free = SpherePlane {
        force_world: Vec3::new(4.0, -6.0, 2.0),
        ..scene()
    };
    let initial = state(
        Vec3::new(1.0, 3.0, -2.0),
        Vec3::new(0.5, -0.2, 0.3),
        Vec3::new(0.1, 0.2, 0.3),
        UnitQuaternion::IDENTITY,
    );
    let result =
        advance_sphere_plane(&free, initial, 0.2, SpherePlaneConfig::default(), || false).unwrap();
    assert_eq!(result.state_after.mode(), SpherePlaneMode::Flight);
    assert!(result.events.is_empty());
    close_vector(
        result.state_after.body().pose().position_world(),
        Vec3::new(1.14, 2.9, -1.92),
    );
    close_vector(
        velocity(result.state_after.body()),
        Vec3::new(0.9, -0.8, 0.5),
    );
    close_vector(spin(result.state_after.body()), Vec3::new(0.1, 0.2, 0.3));
    close_vector(result.contact_impulse_world, Vec3::ZERO);
    close(result.external_force_work, 1.32);
    close(result.dissipated_energy, 0.0);
    balances(&free, initial, 0.2, &result);

    let frictionless = SpherePlane {
        force_world: Vec3::new(0.0, -20.0, 2.0),
        friction_coefficient: 0.0,
        ..scene()
    };
    let initial = state(
        Vec3::new(0.0, RADIUS, 0.0),
        Vec3::new(1.4, 0.0, 0.0),
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
    );
    let result = advance_sphere_plane(
        &frictionless,
        initial,
        0.2,
        SpherePlaneConfig::default(),
        || false,
    )
    .unwrap();
    close_vector(
        result.state_after.body().pose().position_world(),
        Vec3::new(0.28, RADIUS, 0.02),
    );
    close_vector(
        velocity(result.state_after.body()),
        Vec3::new(1.4, 0.0, 0.2),
    );
    close_vector(spin(result.state_after.body()), Vec3::ZERO);
    close_vector(result.contact_impulse_world, Vec3::new(0.0, 4.0, 0.0));
    close(result.external_force_work, 0.04);
    close(result.dissipated_energy, 0.0);
    balances(&frictionless, initial, 0.2, &result);
}
