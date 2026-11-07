use super::*;
use core::f64::consts::PI;
use fs_mbd::sphere_plane::hertz_support::{
    SpherePlaneHertzState, SpherePlaneHertzSupportAdvance, SpherePlaneHertzSupportControls,
    SpherePlaneHertzSupportError, advance_sphere_plane_hertz_support,
    advance_sphere_plane_hertz_support_with,
};
use fs_tribo::sliding_patch::contact::{
    HertzWrenchConfig, HertzWrenchContact, HertzWrenchParameters,
};
use fs_tribo::sliding_patch::{HertzSlidingPatch, HertzSlidingPatchParameters};

const PATCH_RADIUS: f64 = 0.1;
const STATIC_MU: f64 = 0.4;
const KINETIC_MU: f64 = 0.2;

fn law(normal: f64) -> HertzWrenchContact {
    HertzWrenchContact::new(HertzWrenchParameters {
        normal_force_n: normal,
        contact_radius_m: PATCH_RADIUS,
        static_mu: STATIC_MU,
        kinetic_mu: KINETIC_MU,
    })
    .unwrap()
}

fn controls() -> SpherePlaneHertzSupportControls {
    SpherePlaneHertzSupportControls {
        motion: SpherePlaneConfig::default(),
        wrench: HertzWrenchConfig {
            maximum_static_iterations: 96,
            maximum_iterations: 256,
            maximum_evaluations: 4096,
            force_tolerance_n: 2e-9,
            moment_tolerance_nm: 2e-10,
            displacement_tolerance_m: 2e-13,
            energy_tolerance_j: 2e-10,
            relative_matrix_tolerance: 1e-13,
            maximum_force_n: 30.0,
            maximum_spin_moment_nm: 5.0,
        },
        maximum_support_phase_s: 10.0,
        maximum_support_phases: 1024,
        velocity_tolerance_m_s: 2e-10,
        spin_tolerance_rad_s: 2e-10,
        work_tolerance_j: 2e-10,
        maximum_absolute_static_work_residual_j: 2e-9,
    }
}

fn start(velocity: Vec3, omega: Vec3) -> SpherePlaneHertzState {
    SpherePlaneHertzState::new(state(
        Vec3::new(0.0, RADIUS, 0.0),
        velocity,
        omega,
        UnitQuaternion::IDENTITY,
    ))
}

fn audit(before: SpherePlaneHertzState, result: &SpherePlaneHertzSupportAdvance) {
    let m = &result.motion;
    close(
        kinetic(m.state_after.body()) - kinetic(before.motion().body()),
        m.external_force_work
            + m.external_torque_work
            + m.plane_work
            + m.static_contact_work_residual
            - m.dissipated_energy,
    );
    close(
        result.sustained_friction_loss_j + result.impact_loss_j,
        m.dissipated_energy,
    );
    close(m.energy_balance_residual, 0.0);
    close_vector(m.linear_momentum_residual, Vec3::ZERO);
    close_vector(m.angular_momentum_residual, Vec3::ZERO);
    assert!(result.sustained_friction_loss_j >= 0.0);
    assert_eq!(result.state_after.motion(), m.state_after);
}

#[test]
fn pure_spin_stops_at_its_true_time_and_retains_the_static_hysteresis_band() {
    let s = scene();
    let initial = start(Vec3::ZERO, s.normal_world.scale(4.0));
    let moment = 3.0 * PI / 16.0 * KINETIC_MU * 20.0 * PATCH_RADIUS;
    let stop = INERTIA * 4.0 / moment;
    let result = advance_sphere_plane_hertz_support(
        &s,
        initial,
        Vec3::ZERO,
        &law(20.0),
        5.0,
        controls(),
        || false,
    )
    .unwrap();
    assert!(result.state_after.patch_sticking());
    assert!(result.state_after.at_zero_motion());
    assert_eq!(result.support_segments, 2);
    close(result.transitions[0].time_s, stop);
    close(result.support_spin_angle_rad, 0.5 * 4.0 * stop);
    close(result.sustained_friction_loss_j, 0.5 * INERTIA * 16.0);
    close_vector(
        result.support_spin_impulse_world_nms,
        s.normal_world.scale(-4.0 * INERTIA),
    );
    close_vector(spin(result.motion.state_after.body()), Vec3::ZERO);
    let expected =
        UnitQuaternion::from_axis_angle(s.normal_world, result.support_spin_angle_rad).unwrap();
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
    audit(initial, &result);

    // Above kinetic capacity, below static capacity: rest holds, existing spin
    // accelerates. One coefficient cannot reproduce these two original phases.
    let torque = s.normal_world.scale(0.3);
    let mut holding_controls = controls();
    holding_controls.maximum_support_phase_s = 0.1;
    let held = advance_sphere_plane_hertz_support(
        &s,
        result.state_after,
        torque,
        &law(20.0),
        0.7,
        holding_controls,
        || false,
    )
    .unwrap();
    assert!(held.state_after.patch_sticking());
    close(held.sustained_friction_loss_j, 0.0);
    close(held.motion.external_torque_work, 0.0);
    close_vector(held.support_spin_impulse_world_nms, torque.scale(-0.7));
    assert!(held.absolute_static_work_residual_j > 0.0);
    assert!(held.absolute_static_work_residual_j >= held.motion.static_contact_work_residual.abs());
    audit(result.state_after, &held);
    holding_controls.maximum_absolute_static_work_residual_j =
        0.75 * held.absolute_static_work_residual_j;
    let mut held_history = 0usize;
    let mut calls = 0usize;
    let refused = advance_sphere_plane_hertz_support_with(
        &s,
        result.state_after,
        torque,
        &mut held_history,
        0.7,
        holding_controls,
        |old, _, _| {
            calls += 1;
            Ok::<_, ()>((law(20.0), old + 1))
        },
        || false,
    );
    assert!(matches!(
        refused,
        Err(SpherePlaneHertzSupportError::Mechanics(
            SpherePlaneError::Unresolved("cumulative static Hertz work budget")
        ))
    ));
    assert!(calls > 1);
    assert_eq!(held_history, 0);
    let moving = start(Vec3::ZERO, s.normal_world.scale(0.1));
    let accelerated =
        advance_sphere_plane_hertz_support(&s, moving, torque, &law(20.0), 0.7, controls(), || {
            false
        })
        .unwrap();
    assert!(!accelerated.state_after.patch_sticking());
    close(
        spin(accelerated.motion.state_after.body()).dot(s.normal_world),
        0.1 + 0.7 * (0.3 - moment) / INERTIA,
    );
    audit(moving, &accelerated);
}

#[test]
fn an_opposed_motor_crosses_zero_without_claiming_static_capture_or_chattering() {
    let s = scene();
    let initial = start(Vec3::ZERO, s.normal_world.scale(1.0));
    let moment = 3.0 * PI / 16.0 * KINETIC_MU * 20.0 * PATCH_RADIUS;
    let torque = -1.0;
    let stop = INERTIA / (1.0 + moment);
    let at_zero = advance_sphere_plane_hertz_support(
        &s,
        initial,
        s.normal_world.scale(torque),
        &law(20.0),
        stop,
        controls(),
        || false,
    )
    .unwrap();
    assert!(at_zero.state_after.at_zero_motion());
    assert!(!at_zero.state_after.patch_sticking());
    assert!(!at_zero.transitions.last().unwrap().patch_sticking);
    let tail = 0.23;
    let resumed = advance_sphere_plane_hertz_support(
        &s,
        at_zero.state_after,
        s.normal_world.scale(torque),
        &law(20.0),
        tail,
        controls(),
        || false,
    )
    .unwrap();
    let whole = advance_sphere_plane_hertz_support(
        &s,
        initial,
        s.normal_world.scale(torque),
        &law(20.0),
        stop + tail,
        controls(),
        || false,
    )
    .unwrap();
    let expected = tail * (torque + moment) / INERTIA;
    close(
        spin(whole.motion.state_after.body()).dot(s.normal_world),
        expected,
    );
    close_vector(
        spin(resumed.motion.state_after.body()),
        spin(whole.motion.state_after.body()),
    );
    close(
        at_zero.sustained_friction_loss_j + resumed.sustained_friction_loss_j,
        whole.sustained_friction_loss_j,
    );
    assert_eq!(whole.support_segments, 2);
    assert!(!whole.state_after.at_zero_motion());
    assert!(!whole.state_after.patch_sticking());
    audit(initial, &whole);
    audit(at_zero.state_after, &resumed);
}

#[test]
fn pure_shear_recovers_the_original_exact_rolling_trajectory_and_split_continuation() {
    let s = scene();
    let initial = start(Vec3::new(2.0, 0.0, 0.0), Vec3::ZERO);
    let ordinary =
        advance_sphere_plane(&s, initial.motion(), 1.0, controls().motion, || false).unwrap();
    let finite = advance_sphere_plane_hertz_support(
        &s,
        initial,
        Vec3::ZERO,
        &law(20.0),
        1.0,
        controls(),
        || false,
    )
    .unwrap();
    close_vector(
        velocity(finite.motion.state_after.body()),
        velocity(ordinary.state_after.body()),
    );
    close_vector(
        spin(finite.motion.state_after.body()),
        spin(ordinary.state_after.body()),
    );
    close_vector(
        finite.motion.state_after.body().pose().position_world(),
        ordinary.state_after.body().pose().position_world(),
    );
    close(finite.sustained_friction_loss_j, ordinary.dissipated_energy);
    close(finite.support_spin_angle_rad, 0.0);
    assert!(finite.state_after.patch_sticking());
    let first = advance_sphere_plane_hertz_support(
        &s,
        initial,
        Vec3::ZERO,
        &law(20.0),
        0.1,
        controls(),
        || false,
    )
    .unwrap();
    let second = advance_sphere_plane_hertz_support(
        &s,
        first.state_after,
        Vec3::ZERO,
        &law(20.0),
        0.9,
        controls(),
        || false,
    )
    .unwrap();
    close_vector(
        second.motion.state_after.body().pose().position_world(),
        finite.motion.state_after.body().pose().position_world(),
    );
    close(
        first.sustained_friction_loss_j + second.sustained_friction_loss_j,
        finite.sustained_friction_loss_j,
    );
    audit(initial, &finite);
}

#[test]
fn joint_static_support_and_breakaway_use_the_full_hertz_surface() {
    let initial = start(Vec3::ZERO, Vec3::ZERO);
    let mut supported = scene();
    supported.force_world.x = 7.0;
    let hold = advance_sphere_plane_hertz_support(
        &supported,
        initial,
        supported.normal_world.scale(0.1),
        &law(20.0),
        0.4,
        controls(),
        || false,
    )
    .unwrap();
    assert!(hold.state_after.patch_sticking());
    close_vector(
        velocity(hold.motion.state_after.body()),
        Vec3::new(1.0, 0.0, 0.0),
    );
    close(
        spin(hold.motion.state_after.body()).dot(supported.normal_world),
        0.0,
    );
    close(hold.sustained_friction_loss_j, 0.0);
    close_vector(
        hold.support_spin_impulse_world_nms,
        supported.normal_world.scale(-0.04),
    );
    audit(initial, &hold);

    let capacity = STATIC_MU * 20.0;
    let spin_capacity = 3.0 * PI / 16.0 * capacity * PATCH_RADIUS;
    let mut driven = scene();
    driven.force_world.x = 3.5 * 0.95 * capacity;
    let torque = driven.normal_world.scale(0.95 * spin_capacity);
    let broken = advance_sphere_plane_hertz_support(
        &driven,
        initial,
        torque,
        &law(20.0),
        0.2,
        controls(),
        || false,
    )
    .unwrap();
    assert!(!broken.state_after.patch_sticking());
    assert!(spin(broken.motion.state_after.body()).dot(driven.normal_world) > 0.0);
    let shear = broken.motion.contact_impulse_world.x.abs() / 0.2;
    let torsion = broken.support_spin_impulse_world_nms.y.abs() / 0.2;
    assert!(shear < KINETIC_MU * 20.0);
    assert!(torsion < 3.0 * PI / 16.0 * KINETIC_MU * 20.0 * PATCH_RADIUS);
    assert!(broken.sustained_friction_loss_j > 0.0);
    audit(initial, &broken);
}

// Independent continuous RK4 flow of physical COM velocity/spin/position and
// local dissipated power. It shares only the original constitutive surface;
// it has no endpoint projection, held-wrench work or native trajectory engine.
fn reference(
    s: SpherePlane,
    torque: Vec3,
    initial: SpherePlaneHertzState,
    duration: f64,
) -> [f64; 8] {
    let patch = HertzSlidingPatch::new(HertzSlidingPatchParameters {
        normal_force_n: 20.0,
        contact_radius_m: PATCH_RADIUS,
        friction_coefficient: KINETIC_MU,
    })
    .unwrap();
    let v = velocity(initial.motion().body());
    let w = spin(initial.motion().body());
    let x = initial.motion().body().pose().position_world();
    let mut q = [v.x, v.z, w.x, w.y, w.z, x.x, x.z, 0.0];
    let derivative = |q: [f64; 8]| {
        let u = [-q[1] + RADIUS * q[2], -q[0] - RADIUS * q[4]];
        let f = patch.evaluate(u, q[3]).unwrap();
        [
            (s.force_world.x - f.force_n[1]) / MASS,
            (s.force_world.z - f.force_n[0]) / MASS,
            (torque.x + RADIUS * f.force_n[0]) / INERTIA,
            (torque.y + f.spin_moment_nm) / INERTIA,
            (torque.z - RADIUS * f.force_n[1]) / INERTIA,
            q[0],
            q[1],
            f.dissipated_power_w,
        ]
    };
    let h = duration / 8192.0;
    for _ in 0..8192 {
        let a = derivative(q);
        let b = derivative(core::array::from_fn(|i| q[i] + 0.5 * h * a[i]));
        let c = derivative(core::array::from_fn(|i| q[i] + 0.5 * h * b[i]));
        let d = derivative(core::array::from_fn(|i| q[i] + h * c[i]));
        q = core::array::from_fn(|i| q[i] + h / 6.0 * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]));
    }
    q
}

#[test]
fn mixed_sliding_spin_converges_to_independent_continuous_motion_and_heat() {
    let mut s = scene();
    s.force_world.x = 0.7;
    s.force_world.z = -0.4;
    let torque = Vec3::new(0.02, 0.04, -0.01);
    let initial = start(Vec3::new(1.2, 0.0, 0.7), Vec3::new(0.0, 8.0, 0.0));
    let expected = reference(s, torque, initial, 0.2);
    let mut errors = Vec::new();
    for step in [0.04, 0.02, 0.01, 0.005] {
        let mut c = controls();
        c.motion.max_sliding_step_s = step;
        let result =
            advance_sphere_plane_hertz_support(&s, initial, torque, &law(20.0), 0.2, c, || false)
                .unwrap();
        let b = result.motion.state_after.body();
        let v = velocity(b);
        let w = spin(b);
        let x = b.pose().position_world();
        let actual = [
            v.x,
            v.z,
            w.x,
            w.y,
            w.z,
            x.x,
            x.z,
            result.sustained_friction_loss_j,
        ];
        let error = (0..8)
            .map(|i| (actual[i] - expected[i]).abs())
            .fold(0.0_f64, f64::max);
        errors.push(error);
        assert!(!result.state_after.patch_sticking());
        audit(initial, &result);
    }
    for pair in errors.windows(2) {
        assert!(pair[0] / pair[1] > 1.7, "{errors:?}");
    }
    assert!(errors[3] < 0.01, "{errors:?}");
}

#[test]
fn source_phases_follow_real_motion_and_failures_publish_no_history() {
    let s = scene();
    let initial = start(Vec3::ZERO, Vec3::new(0.0, 4.0, 0.0));
    let mut c = controls();
    c.maximum_support_phase_s = 0.2;
    let mut history = 0usize;
    let run = advance_sphere_plane_hertz_support_with(
        &s,
        initial,
        Vec3::ZERO,
        &mut history,
        0.4,
        c,
        |old, phase, _| {
            assert_eq!(phase.normal_force_n, 20.0);
            close(phase.time_s, *old as f64 * 0.2);
            assert!(phase.spin_rad_s > 0.0);
            let mut p = law(phase.normal_force_n).parameters();
            p.contact_radius_m *= 1.0 + *old as f64;
            Ok::<_, &'static str>((HertzWrenchContact::new(p).unwrap(), old + 1))
        },
        || false,
    )
    .unwrap();
    assert_eq!(history, 2);
    let moment = 3.0 * PI / 16.0 * KINETIC_MU * 20.0 * PATCH_RADIUS;
    close(
        spin(run.motion.state_after.body()).y,
        4.0 - 0.2 * 3.0 * moment / INERTIA,
    );
    assert_eq!(run.support_phases, 2);
    audit(initial, &run);
    history = 0;
    let failed = advance_sphere_plane_hertz_support_with(
        &s,
        initial,
        Vec3::ZERO,
        &mut history,
        0.4,
        c,
        |old, _, _| {
            if *old == 1 {
                Err("late source refusal")
            } else {
                Ok((law(20.0), 1))
            }
        },
        || false,
    );
    assert!(matches!(
        failed,
        Err(SpherePlaneHertzSupportError::Source("late source refusal"))
    ));
    assert_eq!(history, 0);
    let wrong_load = advance_sphere_plane_hertz_support_with(
        &s,
        initial,
        Vec3::ZERO,
        &mut history,
        0.4,
        c,
        |old, _, _| Ok::<_, ()>((law(19.0), old + 1)),
        || false,
    );
    assert!(wrong_load.is_err());
    assert_eq!(history, 0);
    let cancelled = advance_sphere_plane_hertz_support_with(
        &s,
        initial,
        Vec3::ZERO,
        &mut history,
        0.4,
        c,
        |old, _, _| Ok::<_, ()>((law(20.0), old + 1)),
        || true,
    );
    assert!(cancelled.is_err());
    assert_eq!(history, 0);
    c.motion.max_sliding_substeps = 1;
    assert!(
        advance_sphere_plane_hertz_support_with(
            &s,
            initial,
            Vec3::ZERO,
            &mut history,
            0.4,
            c,
            |old, _, _| Ok::<_, ()>((law(20.0), old + 1)),
            || false
        )
        .is_err()
    );
    assert_eq!(history, 0);
}

#[test]
fn a_detected_point_impact_hands_real_spin_to_the_separately_declared_support_patch() {
    let mut s = scene();
    s.restitution = 0.0;
    let initial = SpherePlaneHertzState::new(state(
        Vec3::new(0.0, RADIUS + 0.2, 0.0),
        Vec3::new(0.0, -1.0, 0.0),
        Vec3::new(0.0, 0.2, 0.0),
        UnitQuaternion::IDENTITY,
    ));
    let result = advance_sphere_plane_hertz_support(
        &s,
        initial,
        Vec3::ZERO,
        &law(20.0),
        0.8,
        controls(),
        || false,
    )
    .unwrap();
    let event = result
        .motion
        .events
        .iter()
        .find(|e| e.kind == SpherePlaneEventKind::Impact)
        .unwrap();
    close(event.time_s, ((5.0_f64).sqrt() - 1.0) / 10.0);
    close_vector(event.spin_angular_impulse_world, Vec3::ZERO);
    close(result.impact_loss_j, 5.0);
    close(result.sustained_friction_loss_j, 0.5 * INERTIA * 0.2 * 0.2);
    assert!(result.state_after.patch_sticking());
    audit(initial, &result);
    s.force_world = Vec3::new(0.0, 2.0, 0.0);
    let lifted = advance_sphere_plane_hertz_support(
        &s,
        result.state_after,
        Vec3::ZERO,
        &law(20.0),
        0.2,
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(lifted.state_after.motion().mode(), SpherePlaneMode::Flight);
    assert!(!lifted.state_after.patch_sticking());
    assert_eq!(lifted.support_phases, 0);
    audit(result.state_after, &lifted);
}

#[test]
fn moving_plane_work_rotation_and_shared_drift_preserve_the_relative_patch_physics() {
    let mut s = scene();
    s.restitution = 0.0;
    s.force_world.x = 0.7;
    s.force_world.z = -0.4;
    let torque = Vec3::new(0.02, 0.04, -0.01);
    // A resolved normal impact establishes supported continuation in both
    // arbitrary frames, without asking roundoff-sized initial normal velocity
    // of a freshly supplied grazing state to select the contact mode.
    let initial = SpherePlaneHertzState::new(state(
        Vec3::new(0.0, RADIUS + 0.1, 0.0),
        Vec3::new(1.2, -1.0, 0.7),
        Vec3::new(0.0, 8.0, 0.0),
        UnitQuaternion::IDENTITY,
    ));
    let duration = 0.2;
    let original = advance_sphere_plane_hertz_support(
        &s,
        initial,
        torque,
        &law(20.0),
        duration,
        controls(),
        || false,
    )
    .unwrap();
    let rotation = UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, -0.7), 0.6).unwrap();
    let r = |v| rotation.rotate_body_to_world(v);
    let drift = Vec3::new(2.0, 0.1, -0.4);
    let translation = Vec3::new(0.4, -0.8, 0.3);
    let mut changed = s;
    changed.point_world = translation;
    changed.normal_world = r(s.normal_world);
    changed.force_world = r(s.force_world);
    changed.velocity_world = drift;
    let body = initial.motion().body();
    let rotated = SpherePlaneHertzState::new(state(
        r(body.pose().position_world()).add(translation),
        r(velocity(body)).add(drift),
        r(spin(body)),
        rotation,
    ));
    let mut history = ();
    let result = advance_sphere_plane_hertz_support_with(
        &changed,
        rotated,
        r(torque),
        &mut history,
        duration,
        controls(),
        |_, phase, _| Ok::<_, ()>((law(phase.normal_force_n), ())),
        || false,
    )
    .unwrap();
    let final_body = original.motion.state_after.body();
    close_vector(
        velocity(result.motion.state_after.body()),
        r(velocity(final_body)).add(drift),
    );
    close_vector(spin(result.motion.state_after.body()), r(spin(final_body)));
    close_vector(
        result.motion.state_after.body().pose().position_world(),
        r(final_body.pose().position_world())
            .add(translation)
            .add(drift.scale(duration)),
    );
    close(
        result.sustained_friction_loss_j,
        original.sustained_friction_loss_j,
    );
    close(
        result.motion.plane_work,
        result.motion.contact_impulse_world.dot(drift),
    );
    close(
        result.motion.external_torque_work,
        original.motion.external_torque_work,
    );
    audit(rotated, &result);
}
