//! External motor torque through flight, impact, traction and rolling.

use super::*;
use fs_contact::tangential::impact::FrictionLaw;
use fs_mbd::sphere_plane::{
    advance_sphere_plane_with_torque, advance_sphere_plane_with_torque_and_support_friction,
};

fn torque_balances(
    scene: &SpherePlane,
    initial: SpherePlaneState,
    duration: f64,
    torque: Vec3,
    result: &SpherePlaneAdvance,
) {
    close(
        kinetic(result.state_after.body()) - kinetic(initial.body()),
        result.external_force_work + result.external_torque_work + result.plane_work
            - result.dissipated_energy,
    );
    close_vector(result.external_torque_impulse_world, torque.scale(duration));
    close_vector(
        result
            .state_after
            .body()
            .linear_momentum_world()
            .sub(initial.body().linear_momentum_world()),
        scene
            .force_world
            .scale(duration)
            .add(result.contact_impulse_world),
    );
    close(result.energy_balance_residual, 0.0);
    close_vector(result.linear_momentum_residual, Vec3::ZERO);
    close_vector(result.angular_momentum_residual, Vec3::ZERO);
    assert!(result.dissipated_energy >= 0.0);
}

#[test]
fn motor_driven_rolling_has_separate_conveyor_work_and_frame_covariance() {
    let moving = SpherePlane {
        velocity_world: Vec3::new(1.0, 0.0, 0.0),
        ..scene()
    };
    let initial = state(
        Vec3::new(0.0, RADIUS, 0.0),
        moving.velocity_world,
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
    );
    let torque = Vec3::new(0.0, 0.0, -0.7);
    let duration = 0.8;
    let motion = advance_sphere_plane_with_torque(
        &moving,
        initial,
        torque,
        duration,
        SpherePlaneConfig::default(),
        || false,
    )
    .unwrap();
    // I alpha = T + r f and a + r alpha = 0 give f=1, a=1/2, alpha=-1.
    assert_eq!(motion.state_after.mode(), SpherePlaneMode::Rolling);
    close_vector(
        velocity(motion.state_after.body()),
        Vec3::new(1.4, 0.0, 0.0),
    );
    close_vector(spin(motion.state_after.body()), Vec3::new(0.0, 0.0, -0.8));
    close_vector(
        motion.state_after.body().pose().position_world(),
        Vec3::new(0.96, RADIUS, 0.0),
    );
    close_vector(motion.contact_impulse_world, Vec3::new(0.8, 16.0, 0.0));
    close(motion.external_torque_work, 0.224);
    close(motion.plane_work, 0.8);
    close(motion.dissipated_energy, 0.0);
    torque_balances(&moving, initial, duration, torque, &motion);

    // Reverse the motor couple on the accepted rolling state. Static friction
    // transmits braking work back to the conveyor without adding dissipation.
    let brake = Vec3::new(0.0, 0.0, 0.7);
    let braking = advance_sphere_plane_with_torque(
        &moving,
        motion.state_after,
        brake,
        0.4,
        SpherePlaneConfig::default(),
        || false,
    )
    .unwrap();
    close_vector(
        velocity(braking.state_after.body()),
        Vec3::new(1.2, 0.0, 0.0),
    );
    close_vector(spin(braking.state_after.body()), Vec3::new(0.0, 0.0, -0.4));
    close(braking.external_torque_work, -0.168);
    close(braking.plane_work, -0.4);
    close(braking.dissipated_energy, 0.0);
    torque_balances(&moving, motion.state_after, 0.4, brake, &braking);

    let rotation = UnitQuaternion::new(0.5, 0.5, 0.5, 0.5).unwrap();
    let rotate = |v| rotation.rotate_body_to_world(v);
    let shift = Vec3::new(3.0, -4.0, 2.0);
    let transformed_scene = SpherePlane {
        point_world: shift,
        normal_world: rotate(moving.normal_world),
        force_world: rotate(moving.force_world),
        velocity_world: rotate(moving.velocity_world),
        ..moving
    };
    let transformed_initial = state(
        rotate(initial.body().pose().position_world()).add(shift),
        rotate(moving.velocity_world),
        Vec3::ZERO,
        rotation,
    );
    let transformed = advance_sphere_plane_with_torque(
        &transformed_scene,
        transformed_initial,
        rotate(torque),
        duration,
        SpherePlaneConfig::default(),
        || false,
    )
    .unwrap();
    close_vector(
        velocity(transformed.state_after.body()),
        rotate(velocity(motion.state_after.body())),
    );
    close_vector(
        spin(transformed.state_after.body()),
        rotate(spin(motion.state_after.body())),
    );
    close_vector(
        transformed.state_after.body().pose().position_world(),
        rotate(motion.state_after.body().pose().position_world()).add(shift),
    );
    close(
        transformed.external_torque_work,
        motion.external_torque_work,
    );
    close(transformed.plane_work, motion.plane_work);
    torque_balances(
        &transformed_scene,
        transformed_initial,
        duration,
        rotate(torque),
        &transformed,
    );

    // A normal-axis motor performs spin work without point slip or torsional loss.
    let stationary = scene();
    let initial = state(
        Vec3::new(0.0, RADIUS, 0.0),
        Vec3::ZERO,
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
    );
    let normal_torque = Vec3::new(0.0, 0.4, 0.0);
    let normal = advance_sphere_plane_with_torque(
        &stationary,
        initial,
        normal_torque,
        0.5,
        SpherePlaneConfig::default(),
        || false,
    )
    .unwrap();
    assert_eq!(normal.state_after.mode(), SpherePlaneMode::Rolling);
    close_vector(velocity(normal.state_after.body()), Vec3::ZERO);
    close_vector(spin(normal.state_after.body()), Vec3::new(0.0, 1.0, 0.0));
    close(normal.external_torque_work, 0.1);
    close(normal.dissipated_energy, 0.0);
    torque_balances(&stationary, initial, 0.5, normal_torque, &normal);
}

#[test]
fn torque_resolves_traction_breakaway_slip_reversal_and_static_hysteresis() {
    let scene = scene();
    let rest = state(
        Vec3::new(0.0, RADIUS, 0.0),
        Vec3::ZERO,
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
    );
    let torque = Vec3::new(0.0, 0.0, -4.0);
    let duration = 0.4;
    let driven = advance_sphere_plane_with_torque(
        &scene,
        rest,
        torque,
        duration,
        SpherePlaneConfig::default(),
        || false,
    )
    .unwrap();
    // Required static traction 40/7 exceeds 4 N. Sliding has f=4,
    // a=2, alpha=-10 and contact slip acceleration -3.
    assert_eq!(driven.state_after.mode(), SpherePlaneMode::Sliding);
    close_vector(
        velocity(driven.state_after.body()),
        Vec3::new(0.8, 0.0, 0.0),
    );
    close_vector(spin(driven.state_after.body()), Vec3::new(0.0, 0.0, -4.0));
    close(driven.external_torque_work, 3.2);
    close(driven.dissipated_energy, 0.96);
    torque_balances(&scene, rest, duration, torque, &driven);

    let initial = state(
        Vec3::new(0.0, RADIUS, 0.0),
        Vec3::new(1.4, 0.0, 0.0),
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
    );
    let duration = 0.2;
    let reversed = advance_sphere_plane_with_torque(
        &scene,
        initial,
        torque,
        duration,
        SpherePlaneConfig::default(),
        || false,
    )
    .unwrap();
    // Before arrest f=-4, a=-2, alpha=-30: slip reaches zero at 7/85.
    // Static demand is still supercritical, so f=4 and alpha=-10 afterwards.
    let arrest = 7.0 / 85.0;
    let remaining = duration - arrest;
    let reversal = reversed
        .events
        .iter()
        .find(|e| e.kind == SpherePlaneEventKind::Sliding && e.time_s > 0.0)
        .unwrap();
    close(reversal.time_s, arrest);
    let impulse_x = 4.0 * duration - 8.0 * arrest;
    close_vector(
        velocity(reversed.state_after.body()),
        Vec3::new(1.4 + impulse_x / MASS, 0.0, 0.0),
    );
    close_vector(
        spin(reversed.state_after.body()),
        Vec3::new(0.0, 0.0, (-4.0 * duration + RADIUS * impulse_x) / INERTIA),
    );
    let angle = -15.0 * arrest * arrest - 30.0 * arrest * remaining - 5.0 * remaining * remaining;
    close(reversed.external_torque_work, -4.0 * angle);
    close(
        reversed.dissipated_energy,
        2.8 * arrest + 6.0 * remaining * remaining,
    );
    torque_balances(&scene, initial, duration, torque, &reversed);

    let law = FrictionLaw::Coulomb {
        static_mu: 0.2,
        kinetic_mu: 0.1,
    };
    let torque = Vec3::new(0.0, 0.0, -2.0);
    let held = advance_sphere_plane_with_torque_and_support_friction(
        &scene,
        rest,
        torque,
        &law,
        0.2,
        SpherePlaneConfig::default(),
        || false,
    )
    .unwrap();
    assert_eq!(held.state_after.mode(), SpherePlaneMode::Rolling);
    close(held.dissipated_energy, 0.0);
    // The same drive sustains existing slip: the kinetic limit cannot arrest it.
    let slipping = state(
        Vec3::new(0.0, RADIUS, 0.0),
        Vec3::new(-0.25, 0.0, 0.0),
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
    );
    let moving = advance_sphere_plane_with_torque_and_support_friction(
        &scene,
        slipping,
        torque,
        &law,
        0.2,
        SpherePlaneConfig::default(),
        || false,
    )
    .unwrap();
    assert_eq!(moving.state_after.mode(), SpherePlaneMode::Sliding);
    close_vector(
        velocity(moving.state_after.body()),
        Vec3::new(-0.05, 0.0, 0.0),
    );
    close_vector(spin(moving.state_after.body()), Vec3::new(0.0, 0.0, -1.0));
    close(moving.dissipated_energy, 0.16);
    torque_balances(&scene, slipping, 0.2, torque, &moving);
}

#[test]
fn flight_torque_changes_impact_friction_and_survives_split_continuation() {
    let scene = scene();
    let initial = state(
        Vec3::new(0.0, RADIUS + 0.2, 0.0),
        Vec3::ZERO,
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
    );
    let torque = Vec3::new(0.0, 0.0, 0.4);
    let motion = advance_sphere_plane_with_torque(
        &scene,
        initial,
        torque,
        0.3,
        SpherePlaneConfig::default(),
        || false,
    )
    .unwrap();
    let impact = motion
        .events
        .iter()
        .find(|e| e.kind == SpherePlaneEventKind::Impact)
        .unwrap();
    close(impact.time_s, 0.2);
    close(impact.normal_impulse, 6.0);
    // Pre-impact omega_z=.4 creates slip=.2, hence Jx=-.2/(1/m+r²/I).
    close_vector(
        impact.tangent_impulse_world,
        Vec3::new(-4.0 / 35.0, 0.0, 0.0),
    );
    close_vector(
        velocity(motion.state_after.body()),
        Vec3::new(-2.0 / 35.0, 0.0, 0.0),
    );
    close_vector(
        spin(motion.state_after.body()),
        Vec3::new(0.0, 0.0, 11.0 / 35.0),
    );
    close_vector(
        motion.state_after.body().pose().position_world(),
        Vec3::new(-1.0 / 175.0, 0.55, 0.0),
    );
    close(motion.dissipated_energy, 3.0 + 2.0 / 175.0);
    close(
        motion.external_torque_work,
        0.4 * (0.04 + 2.0 / 175.0 + 0.01),
    );
    torque_balances(&scene, initial, 0.3, torque, &motion);

    let first = advance_sphere_plane_with_torque(
        &scene,
        initial,
        torque,
        0.12,
        SpherePlaneConfig::default(),
        || false,
    )
    .unwrap();
    let second = advance_sphere_plane_with_torque(
        &scene,
        first.state_after,
        torque,
        0.18,
        SpherePlaneConfig::default(),
        || false,
    )
    .unwrap();
    close_vector(
        second.state_after.body().pose().position_world(),
        motion.state_after.body().pose().position_world(),
    );
    close_vector(
        velocity(second.state_after.body()),
        velocity(motion.state_after.body()),
    );
    close_vector(
        spin(second.state_after.body()),
        spin(motion.state_after.body()),
    );
    close(
        first.external_torque_work + second.external_torque_work,
        motion.external_torque_work,
    );
    close(
        first.dissipated_energy + second.dissipated_energy,
        motion.dissipated_energy,
    );
}

/// Independent world-coordinate Newton/Euler RK4 reference for sliding before arrest.
fn curved_reference(duration: f64, torque: Vec3) -> (Vec3, Vec3) {
    let arm = Vec3::new(0.0, -RADIUS, 0.0);
    let rhs = |v: Vec3, omega: Vec3| {
        let slip = v.add(omega.cross(arm));
        let speed = slip.norm_squared().sqrt();
        assert!(speed > 0.1, "reference remains in resolved sliding");
        let force = slip.scale(-4.0 / speed);
        (
            force.scale(1.0 / MASS),
            torque.add(arm.cross(force)).scale(1.0 / INERTIA),
        )
    };
    let n = 4096;
    let h = duration / n as f64;
    let mut v = Vec3::new(1.4, 0.0, 0.0);
    let mut omega = Vec3::new(0.0, 1.7, 0.0);
    for _ in 0..n {
        let (a, b) = rhs(v, omega);
        let (c, d) = rhs(v.add(a.scale(0.5 * h)), omega.add(b.scale(0.5 * h)));
        let (e, f) = rhs(v.add(c.scale(0.5 * h)), omega.add(d.scale(0.5 * h)));
        let (g, j) = rhs(v.add(e.scale(h)), omega.add(f.scale(h)));
        v = v.add(a.add(c.scale(2.0)).add(e.scale(2.0)).add(g).scale(h / 6.0));
        omega = omega.add(b.add(d.scale(2.0)).add(f.scale(2.0)).add(j).scale(h / 6.0));
    }
    (v, omega)
}

#[test]
fn transverse_torque_sliding_converges_to_independent_newton_euler_motion() {
    let scene = scene();
    let torque = Vec3::new(0.3, 0.2, -0.4);
    let initial = state(
        Vec3::new(0.0, RADIUS, 0.0),
        Vec3::new(1.4, 0.0, 0.0),
        Vec3::new(0.0, 1.7, 0.0),
        UnitQuaternion::IDENTITY,
    );
    let duration = 0.08;
    let (expected_v, expected_omega) = curved_reference(duration, torque);
    let mut errors = Vec::new();
    for h in [0.02, 0.01, 0.005] {
        let motion = advance_sphere_plane_with_torque(
            &scene,
            initial,
            torque,
            duration,
            SpherePlaneConfig {
                max_sliding_step_s: h,
                max_attitude_step_s: 0.0005,
                ..SpherePlaneConfig::default()
            },
            || false,
        )
        .unwrap();
        assert!(motion.sliding_substeps >= 4);
        assert_eq!(motion.state_after.mode(), SpherePlaneMode::Sliding);
        torque_balances(&scene, initial, duration, torque, &motion);
        errors.push(
            velocity(motion.state_after.body())
                .sub(expected_v)
                .norm_squared()
                .sqrt()
                + RADIUS
                    * spin(motion.state_after.body())
                        .sub(expected_omega)
                        .norm_squared()
                        .sqrt(),
        );
    }
    for pair in errors.windows(2) {
        assert!((1.7..2.4).contains(&(pair[0] / pair[1])), "{errors:?}");
    }
}

#[test]
fn torque_zero_parity_and_late_refusals_preserve_original_state() {
    let scene = scene();
    let initial = state(
        Vec3::new(0.0, 2.0, 0.0),
        Vec3::new(0.2, 0.0, 0.0),
        Vec3::new(0.0, 2.0, 0.0),
        UnitQuaternion::IDENTITY,
    );
    let unchanged = initial;
    let old =
        advance_sphere_plane(&scene, initial, 0.1, SpherePlaneConfig::default(), || false).unwrap();
    let zero = advance_sphere_plane_with_torque(
        &scene,
        initial,
        Vec3::ZERO,
        0.1,
        SpherePlaneConfig::default(),
        || false,
    )
    .unwrap();
    assert_eq!(old.state_after, zero.state_after);
    assert_eq!(old.external_force_work, zero.external_force_work);
    assert_eq!(old.dissipated_energy, zero.dissipated_energy);
    assert_eq!(zero.external_torque_work, 0.0);
    assert_eq!(zero.external_torque_impulse_world, Vec3::ZERO);

    assert!(matches!(
        advance_sphere_plane_with_torque(
            &scene,
            initial,
            Vec3::new(f64::NAN, 0.0, 0.0),
            0.1,
            SpherePlaneConfig::default(),
            || false,
        ),
        Err(SpherePlaneError::Invalid(_))
    ));
    let torque = Vec3::new(0.0, 0.0, 0.4);
    let config = SpherePlaneConfig {
        max_attitude_step_s: 0.001,
        max_attitude_substeps: 1,
        ..SpherePlaneConfig::default()
    };
    assert!(matches!(
        advance_sphere_plane_with_torque(&scene, initial, torque, 0.1, config, || false,),
        Err(SpherePlaneError::AttitudeBudget)
    ));
    let mut polls = 0;
    let config = SpherePlaneConfig {
        max_attitude_step_s: 0.001,
        ..SpherePlaneConfig::default()
    };
    assert!(matches!(
        advance_sphere_plane_with_torque(&scene, initial, torque, 0.1, config, || {
            polls += 1;
            polls >= 9
        },),
        Err(SpherePlaneError::Cancelled)
    ));
    assert!(polls >= 9);
    assert_eq!(initial, unchanged);
    let first =
        advance_sphere_plane_with_torque(&scene, initial, torque, 0.1, config, || false).unwrap();
    let retry =
        advance_sphere_plane_with_torque(&scene, initial, torque, 0.1, config, || false).unwrap();
    assert_eq!(first.state_after, retry.state_after);
    assert_eq!(first.external_torque_work, retry.external_torque_work);
    torque_balances(&scene, initial, 0.1, torque, &first);
}
