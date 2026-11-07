use super::*;
use core::f64::consts::PI;
use fs_mbd::sphere_plane::hertz_impacts::{
    SpherePlaneHertzImpactAdvance, SpherePlaneHertzImpactControls, SpherePlaneHertzImpactError,
    advance_sphere_plane_hertz_impacts, advance_sphere_plane_hertz_impacts_with,
};
use fs_tribo::FrictionLaw;
use fs_tribo::sliding_patch::contact::HertzWrenchMode;
use fs_tribo::sliding_patch::contact::impulse::{
    HertzImpulseControls, HertzImpulseError, HertzImpulsePatch,
};

const SUPPORT: FrictionLaw = FrictionLaw::Coulomb {
    static_mu: 0.4,
    kinetic_mu: 0.2,
};

fn controls() -> SpherePlaneHertzImpactControls {
    SpherePlaneHertzImpactControls {
        motion: SpherePlaneConfig::default(),
        impulse: HertzImpulseControls {
            maximum_static_iterations: 96,
            maximum_iterations: 256,
            maximum_evaluations: 2048,
            impulse_tolerance_ns: 2e-9,
            spin_impulse_tolerance_nms: 5e-10,
            velocity_tolerance_m_s: 2e-13,
            energy_tolerance_j: 2e-10,
            relative_matrix_tolerance: 1e-13,
            maximum_impulse_ns: 20.0,
            maximum_spin_impulse_nms: 5.0,
        },
    }
}

fn audit(result: &SpherePlaneHertzImpactAdvance, initial: SpherePlaneState) {
    let m = &result.motion;
    close(
        kinetic(m.state_after.body()) - kinetic(initial.body()),
        m.external_force_work + m.external_torque_work + m.plane_work - m.dissipated_energy,
    );
    close(m.energy_balance_residual, 0.0);
    close_vector(m.linear_momentum_residual, Vec3::ZERO);
    close_vector(m.angular_momentum_residual, Vec3::ZERO);
    assert!(m.gap_after >= -2e-12);
    assert!(result.impact_loss_j >= 0.0);
    let loss: f64 = result
        .impacts
        .iter()
        .map(|event| {
            let phase = event.phase;
            -0.5 * phase.normal_impulse_ns
                * (phase.normal_velocity_before_m_s + phase.normal_velocity_after_m_s)
                + event.response.dissipated_j
        })
        .sum();
    close(result.impact_loss_j, loss);
}

#[test]
fn pure_spin_impacts_remove_actual_spin_and_advance_the_changed_attitude() {
    let mut s = scene();
    s.force_world = Vec3::ZERO;
    let patch = HertzImpulsePatch::new(0.1, 0.2).unwrap();
    for initial_spin in [8.0, -8.0, 0.05, 0.0] {
        let initial = state(
            Vec3::new(0.0, RADIUS, 0.0),
            Vec3::new(0.0, -2.0, 0.0),
            s.normal_world.scale(initial_spin),
            UnitQuaternion::IDENTITY,
        );
        let duration = 0.3;
        let result = advance_sphere_plane_hertz_impacts(
            &s,
            initial,
            Vec3::ZERO,
            &SUPPORT,
            &patch,
            duration,
            controls(),
            || false,
        )
        .unwrap();
        assert_eq!(result.impacts.len(), 1);
        let event = result.impacts[0];
        let capacity = 3.0 * PI / 16.0 * 0.1 * 0.2 * 6.0;
        let impulse = -initial_spin.signum() * (INERTIA * initial_spin.abs()).min(capacity);
        let final_spin = initial_spin + impulse / INERTIA;
        close(event.response.spin_impulse_nms, impulse);
        close_vector(
            spin(result.motion.state_after.body()),
            s.normal_world.scale(final_spin),
        );
        close_vector(
            result.motion.contact_impulse_world,
            Vec3::new(0.0, 6.0, 0.0),
        );
        close_vector(
            result.motion.events[0].spin_angular_impulse_world,
            s.normal_world.scale(impulse),
        );
        close_vector(
            result.motion.contact_angular_impulse_world,
            s.normal_world.scale(impulse),
        );
        close(
            result.impact_loss_j,
            3.0 + 0.5 * INERTIA * (initial_spin * initial_spin - final_spin * final_spin),
        );
        close_vector(
            result.motion.state_after.body().pose().position_world(),
            Vec3::new(0.0, RADIUS + duration, 0.0),
        );
        let expected =
            UnitQuaternion::from_axis_angle(s.normal_world, final_spin * duration).unwrap();
        let axis = Vec3::new(1.0, 0.0, 0.0);
        close_vector(
            result
                .motion
                .state_after
                .body()
                .pose()
                .orientation()
                .rotate_body_to_world(axis),
            expected.rotate_body_to_world(axis),
        );
        if initial_spin.abs() <= 0.05 {
            assert_eq!(event.response.mode, HertzWrenchMode::Sticking);
        }
        audit(&result, initial);
    }
}

#[test]
fn mixed_sliding_and_spin_share_one_budget_and_recoil_through_real_momenta() {
    let mut s = scene();
    s.force_world = Vec3::ZERO;
    let initial = state(
        Vec3::new(0.3, RADIUS, -0.2),
        Vec3::new(2.0, -2.0, 1.0),
        Vec3::new(0.0, 8.0, 0.0),
        UnitQuaternion::IDENTITY,
    );
    let patch = HertzImpulsePatch::new(0.1, 0.2).unwrap();
    let result = advance_sphere_plane_hertz_impacts(
        &s,
        initial,
        Vec3::ZERO,
        &SUPPORT,
        &patch,
        0.1,
        controls(),
        || false,
    )
    .unwrap();
    let event = result.impacts[0];
    let tangent = event.phase.tangent_axes_world[0]
        .scale(event.response.impulse_ns[0])
        .add(event.phase.tangent_axes_world[1].scale(event.response.impulse_ns[1]));
    let normal = s.normal_world.scale(event.phase.normal_impulse_ns);
    let free = s.normal_world.scale(event.response.spin_impulse_nms);
    let actual = result.motion.state_after.body();
    close_vector(
        velocity(actual),
        velocity(initial.body()).add(tangent.add(normal).scale(1.0 / MASS)),
    );
    close_vector(
        spin(actual),
        spin(initial.body()).add(
            s.normal_world
                .scale(-RADIUS)
                .cross(tangent)
                .add(free)
                .scale(1.0 / INERTIA),
        ),
    );
    close_vector(
        actual.pose().position_world(),
        initial
            .body()
            .pose()
            .position_world()
            .add(velocity(actual).scale(0.1)),
    );
    close_vector(
        result.motion.contact_angular_impulse_world,
        event
            .phase
            .contact_point_world_m
            .cross(tangent.add(normal))
            .add(free),
    );
    assert!(tangent.norm_squared().sqrt() < 0.2 * event.phase.normal_impulse_ns);
    assert!(
        free.norm_squared().sqrt() < 3.0 * PI / 16.0 * 0.1 * 0.2 * event.phase.normal_impulse_ns
    );
    assert!(event.response.spin_impulse_nms < 0.0);
    assert_eq!(event.response.mode, HertzWrenchMode::Sliding);
    audit(&result, initial);
}

#[test]
fn detected_motor_driven_spin_impact_is_covariant_and_continues_across_calls() {
    let mut s = scene();
    s.force_world = Vec3::ZERO;
    let position = Vec3::new(0.0, 1.0, 0.0);
    let v = Vec3::new(0.0, -2.0, 0.0);
    let omega = Vec3::new(0.0, 8.0, 0.0);
    let torque = Vec3::new(0.0, 0.2, 0.0);
    let initial = state(position, v, omega, UnitQuaternion::IDENTITY);
    let patch = HertzImpulsePatch::new(0.1, 0.2).unwrap();
    let duration = 0.4;
    let whole = advance_sphere_plane_hertz_impacts(
        &s,
        initial,
        torque,
        &SUPPORT,
        &patch,
        duration,
        controls(),
        || false,
    )
    .unwrap();
    close(whole.impacts[0].phase.time_s, 0.25);
    close(whole.impacts[0].phase.spin_before_rad_s, 8.25);
    let first = advance_sphere_plane_hertz_impacts(
        &s,
        initial,
        torque,
        &SUPPORT,
        &patch,
        0.2,
        controls(),
        || false,
    )
    .unwrap();
    assert!(first.impacts.is_empty());
    let second = advance_sphere_plane_hertz_impacts(
        &s,
        first.motion.state_after,
        torque,
        &SUPPORT,
        &patch,
        0.2,
        controls(),
        || false,
    )
    .unwrap();
    close_vector(
        spin(whole.motion.state_after.body()),
        spin(second.motion.state_after.body()),
    );
    close_vector(
        whole.motion.state_after.body().pose().position_world(),
        second.motion.state_after.body().pose().position_world(),
    );
    close(
        whole.motion.external_torque_work,
        first.motion.external_torque_work + second.motion.external_torque_work,
    );
    close(whole.impact_loss_j, second.impact_loss_j);
    audit(&whole, initial);
    let rotation = UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, -1.0), 0.7).unwrap();
    let rotate = |x| rotation.rotate_body_to_world(x);
    let boost = Vec3::new(0.7, -0.2, 0.3);
    let mut moved_scene = s;
    moved_scene.normal_world = rotate(s.normal_world);
    moved_scene.velocity_world = boost;
    let moved_initial = state(
        rotate(position),
        rotate(v).add(boost),
        rotate(omega),
        rotation,
    );
    let moved = advance_sphere_plane_hertz_impacts(
        &moved_scene,
        moved_initial,
        rotate(torque),
        &SUPPORT,
        &patch,
        duration,
        controls(),
        || false,
    )
    .unwrap();
    close(moved.impacts[0].phase.time_s, 0.25);
    close_vector(
        velocity(moved.motion.state_after.body()),
        rotate(velocity(whole.motion.state_after.body())).add(boost),
    );
    close_vector(
        spin(moved.motion.state_after.body()),
        rotate(spin(whole.motion.state_after.body())),
    );
    close_vector(
        moved.motion.state_after.body().pose().position_world(),
        rotate(whole.motion.state_after.body().pose().position_world()).add(boost.scale(duration)),
    );
    close_vector(
        moved.motion.events[0].spin_angular_impulse_world,
        rotate(whole.motion.events[0].spin_angular_impulse_world),
    );
    close(
        moved.motion.plane_work,
        moved.motion.contact_impulse_world.dot(boost),
    );
    close(
        moved.motion.external_torque_work,
        whole.motion.external_torque_work,
    );
    close(moved.impact_loss_j, whole.impact_loss_j);
    audit(&moved, moved_initial);
}

#[test]
fn sustained_support_keeps_its_separate_law_after_the_finite_patch_impact() {
    let mut s = scene();
    s.restitution = 0.0;
    let initial = state(
        Vec3::new(0.0, RADIUS, 0.0),
        Vec3::new(2.0, -2.0, 0.0),
        Vec3::new(0.0, 8.0, 0.0),
        UnitQuaternion::IDENTITY,
    );
    let patch = HertzImpulsePatch::new(0.1, 0.2).unwrap();
    let result = advance_sphere_plane_hertz_impacts(
        &s,
        initial,
        Vec3::ZERO,
        &SUPPORT,
        &patch,
        1.0,
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(result.impacts.len(), 1);
    assert_eq!(result.motion.state_after.mode(), SpherePlaneMode::Rolling);
    close(
        spin(result.motion.state_after.body()).y,
        result.impacts[0].response.spin_after_rad_s,
    );
    assert!(result.motion.dissipated_energy > result.impact_loss_j);
    assert!(
        result
            .motion
            .events
            .iter()
            .filter(|event| event.kind != SpherePlaneEventKind::Impact)
            .all(|event| event.spin_angular_impulse_world == Vec3::ZERO)
    );
    audit(&result, initial);
    let no_patch_friction = HertzImpulsePatch::new(0.1, 0.0).unwrap();
    let free = advance_sphere_plane_hertz_impacts(
        &s,
        initial,
        Vec3::ZERO,
        &SUPPORT,
        &no_patch_friction,
        1.0,
        controls(),
        || false,
    )
    .unwrap();
    close(spin(free.motion.state_after.body()).y, 8.0);
    assert_eq!(free.impacts[0].response.spin_impulse_nms, 0.0);
    assert_eq!(free.impacts[0].response.impulse_ns, [0.0; 2]);
    audit(&free, initial);
}

#[test]
fn event_sources_use_real_incoming_motion_and_commit_only_the_complete_trajectory() {
    let s = scene();
    let initial = state(
        Vec3::new(0.0, RADIUS, 0.0),
        Vec3::new(0.0, -2.0, 0.0),
        Vec3::new(0.0, 8.0, 0.0),
        UnitQuaternion::IDENTITY,
    );
    let mut history = Vec::<(f64, f64)>::new();
    let result = advance_sphere_plane_hertz_impacts_with(
        &s,
        initial,
        Vec3::ZERO,
        &SUPPORT,
        &mut history,
        0.25,
        controls(),
        |old, phase, _| {
            let mut next = old.clone();
            next.push((phase.time_s, phase.spin_before_rad_s));
            let radius = if old.is_empty() { 0.05 } else { 0.1 };
            Ok::<_, &'static str>((HertzImpulsePatch::new(radius, 0.2).unwrap(), next))
        },
        || false,
    )
    .unwrap();
    assert_eq!(history.len(), 2);
    close(history[0].0, 0.0);
    close(history[1].0, 0.2);
    close(history[1].1, result.impacts[0].response.spin_after_rad_s);
    close(result.impacts[0].patch.radius_m(), 0.05);
    close(result.impacts[1].patch.radius_m(), 0.1);
    audit(&result, initial);
    let mut refused = 0usize;
    let error = advance_sphere_plane_hertz_impacts_with(
        &s,
        initial,
        Vec3::ZERO,
        &SUPPORT,
        &mut refused,
        0.25,
        controls(),
        |old, _, _| {
            if *old == 1 {
                Err("second source unavailable")
            } else {
                Ok((HertzImpulsePatch::new(0.05, 0.2).unwrap(), old + 1))
            }
        },
        || false,
    )
    .unwrap_err();
    assert!(matches!(
        error,
        SpherePlaneHertzImpactError::Source("second source unavailable")
    ));
    assert_eq!(refused, 0);
}

#[test]
fn budgets_and_late_cancellation_leave_source_history_unchanged() {
    let s = scene();
    let initial = state(
        Vec3::new(0.0, RADIUS, 0.0),
        Vec3::new(0.0, -2.0, 0.0),
        Vec3::new(0.0, 8.0, 0.0),
        UnitQuaternion::IDENTITY,
    );
    for case in 0..3 {
        let mut controls = controls();
        if case == 0 {
            controls.motion.max_events = 1;
        }
        if case == 1 {
            controls.impulse.maximum_spin_impulse_nms = 1e-6;
        }
        let cancelled = std::cell::Cell::new(false);
        let mut history = 0usize;
        let error = advance_sphere_plane_hertz_impacts_with(
            &s,
            initial,
            Vec3::ZERO,
            &SUPPORT,
            &mut history,
            0.25,
            controls,
            |old, _, _| {
                if case == 2 && *old == 1 {
                    cancelled.set(true);
                }
                Ok::<_, ()>((HertzImpulsePatch::new(0.1, 0.2).unwrap(), old + 1))
            },
            || cancelled.get(),
        )
        .unwrap_err();
        match case {
            0 => assert!(matches!(
                error,
                SpherePlaneHertzImpactError::Mechanics(SpherePlaneError::EventBudget)
            )),
            1 => assert!(matches!(
                error,
                SpherePlaneHertzImpactError::Law(HertzImpulseError::ImpulseBudget { .. })
            )),
            _ => assert!(matches!(
                error,
                SpherePlaneHertzImpactError::Mechanics(SpherePlaneError::Cancelled)
            )),
        }
        assert_eq!(history, 0);
    }
}
