use super::*;
use core::cell::{Cell, RefCell};
use core::f64::consts::PI;
use fs_mbd::sphere_plane::hertz_contact::{
    SpherePlaneHertzContactAdvance, SpherePlaneHertzContactControls, SpherePlaneHertzContactError,
    advance_sphere_plane_hertz_contact, advance_sphere_plane_hertz_contact_with,
};
use fs_mbd::sphere_plane::hertz_impacts::SpherePlaneHertzImpactError;
use fs_mbd::sphere_plane::hertz_support::{
    SpherePlaneHertzState, SpherePlaneHertzSupportControls, SpherePlaneHertzSupportError,
};
use fs_tribo::sliding_patch::contact::impulse::{HertzImpulseControls, HertzImpulsePatch};
use fs_tribo::sliding_patch::contact::{
    HertzWrenchConfig, HertzWrenchContact, HertzWrenchParameters,
};

const IMPACT_RADIUS: f64 = 0.06;
const IMPACT_MU: f64 = 0.3;
const SUPPORT_RADIUS: f64 = 0.15;
const SUPPORT_MU: f64 = 0.25;

fn controls() -> SpherePlaneHertzContactControls {
    SpherePlaneHertzContactControls {
        support: SpherePlaneHertzSupportControls {
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
            maximum_support_phases: 128,
            velocity_tolerance_m_s: 2e-10,
            spin_tolerance_rad_s: 2e-10,
            work_tolerance_j: 2e-10,
            maximum_absolute_static_work_residual_j: 2e-9,
        },
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

fn support_law(normal: f64) -> HertzWrenchContact {
    HertzWrenchContact::new(HertzWrenchParameters {
        normal_force_n: normal,
        contact_radius_m: SUPPORT_RADIUS,
        static_mu: 0.45,
        kinetic_mu: SUPPORT_MU,
    })
    .unwrap()
}

fn impact_patch() -> HertzImpulsePatch {
    HertzImpulsePatch::new(IMPACT_RADIUS, IMPACT_MU).unwrap()
}

fn incoming() -> SpherePlaneHertzState {
    SpherePlaneHertzState::new(state(
        Vec3::new(0.0, RADIUS + 0.2, 0.0),
        Vec3::new(0.0, -1.0, 0.0),
        Vec3::new(0.0, 4.0, 0.0),
        UnitQuaternion::IDENTITY,
    ))
}

fn audit(initial: SpherePlaneHertzState, result: &SpherePlaneHertzContactAdvance) {
    let support = &result.support;
    let motion = &support.motion;
    close(
        kinetic(motion.state_after.body()) - kinetic(initial.motion().body()),
        motion.external_force_work
            + motion.external_torque_work
            + motion.plane_work
            + motion.static_contact_work_residual
            - motion.dissipated_energy,
    );
    close(
        motion.dissipated_energy,
        support.sustained_friction_loss_j + support.impact_loss_j,
    );
    close(motion.energy_balance_residual, 0.0);
    close_vector(motion.linear_momentum_residual, Vec3::ZERO);
    close_vector(motion.angular_momentum_residual, Vec3::ZERO);
    close(
        support.impact_loss_j,
        result
            .impacts
            .iter()
            .map(|event| {
                -0.5 * event.phase.normal_impulse_ns
                    * (event.phase.normal_velocity_before_m_s
                        + event.phase.normal_velocity_after_m_s)
                    + event.response.dissipated_j
            })
            .sum(),
    );
    assert_eq!(
        result.impacts.len(),
        motion
            .events
            .iter()
            .filter(|e| e.kind == SpherePlaneEventKind::Impact)
            .count()
    );
}

#[test]
fn distinct_impact_and_support_patches_drive_one_real_spin_stopping_trajectory() {
    let mut s = scene();
    s.restitution = 0.0;
    let initial = incoming();
    let duration = 3.0;
    let result = advance_sphere_plane_hertz_contact(
        &s,
        initial,
        Vec3::ZERO,
        &support_law(20.0),
        &impact_patch(),
        duration,
        controls(),
        || false,
    )
    .unwrap();
    let impact_time = ((5.0_f64).sqrt() - 1.0) / 10.0;
    let normal_impulse = MASS * (5.0_f64).sqrt();
    let angular_impulse = -(3.0 * PI / 16.0) * IMPACT_RADIUS * IMPACT_MU * normal_impulse;
    let spin_after_impact = 4.0 + angular_impulse / INERTIA;
    let support_moment = (3.0 * PI / 16.0) * SUPPORT_RADIUS * SUPPORT_MU * 20.0;
    let stop_duration = INERTIA * spin_after_impact / support_moment;
    assert_eq!(result.impacts.len(), 1);
    let event = result.impacts[0];
    close(event.phase.time_s, impact_time);
    close(event.patch.radius_m(), IMPACT_RADIUS);
    close(event.phase.normal_impulse_ns, normal_impulse);
    close(event.response.spin_impulse_nms, angular_impulse);
    close(
        result.support.transitions[0].time_s,
        impact_time + stop_duration,
    );
    assert!(result.support.state_after.patch_sticking());
    close_vector(spin(result.support.motion.state_after.body()), Vec3::ZERO);
    close_vector(
        result
            .support
            .motion
            .state_after
            .body()
            .pose()
            .position_world(),
        Vec3::new(0.0, RADIUS, 0.0),
    );
    close(
        result.support.sustained_friction_loss_j,
        0.5 * INERTIA * spin_after_impact.powi(2),
    );
    close(
        result.support.impact_loss_j,
        5.0 + 0.5 * INERTIA * (16.0 - spin_after_impact.powi(2)),
    );
    close(result.support.motion.dissipated_energy, 6.6);
    close(result.support.motion.external_force_work, 4.0);
    close_vector(
        result.support.support_spin_impulse_world_nms,
        s.normal_world.scale(-INERTIA * spin_after_impact),
    );
    close_vector(
        result.support.motion.contact_angular_impulse_world,
        s.normal_world.scale(-INERTIA * 4.0),
    );
    let angle = 4.0 * impact_time + 0.5 * spin_after_impact * stop_duration;
    let orientation = UnitQuaternion::from_axis_angle(s.normal_world, angle).unwrap();
    let axis = Vec3::new(1.0, 0.0, 0.0);
    close_vector(
        result
            .support
            .motion
            .state_after
            .body()
            .pose()
            .orientation()
            .rotate_body_to_world(axis),
        orientation.rotate_body_to_world(axis),
    );
    audit(initial, &result);

    // The inclusive impact event is published exactly once when an advance
    // ends at its actual detected time; the next call starts on the support law.
    let cut = event.phase.time_s;
    let first = advance_sphere_plane_hertz_contact(
        &s,
        initial,
        Vec3::ZERO,
        &support_law(20.0),
        &impact_patch(),
        cut,
        controls(),
        || false,
    )
    .unwrap();
    let second = advance_sphere_plane_hertz_contact(
        &s,
        first.support.state_after,
        Vec3::ZERO,
        &support_law(20.0),
        &impact_patch(),
        duration - cut,
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(first.impacts.len(), 1);
    assert!(second.impacts.is_empty());
    assert_eq!(first.support.support_phases, 0);
    close_vector(
        spin(first.support.motion.state_after.body()),
        s.normal_world.scale(spin_after_impact),
    );
    close(
        first.support.motion.dissipated_energy + second.support.motion.dissipated_energy,
        result.support.motion.dissipated_energy,
    );
    close_vector(
        second
            .support
            .motion
            .state_after
            .body()
            .pose()
            .orientation()
            .rotate_body_to_world(axis),
        orientation.rotate_body_to_world(axis),
    );
    audit(initial, &first);
    audit(first.support.state_after, &second);
}

#[test]
fn callbacks_observe_actual_postimpact_motion_and_both_histories_roll_back_together() {
    let mut s = scene();
    s.restitution = 0.0;
    let initial = incoming();
    let mut c = controls();
    c.support.maximum_support_phase_s = 0.15;
    let mut support_history = 20usize;
    let mut impact_history = 10usize;
    let order = RefCell::new(Vec::new());
    let after_impact = Cell::new(f64::NAN);
    let result = advance_sphere_plane_hertz_contact_with(
        &s,
        initial,
        Vec3::ZERO,
        &mut support_history,
        &mut impact_history,
        0.6,
        c,
        |old, phase, _| {
            order.borrow_mut().push(('s', *old, phase.time_s));
            assert!(after_impact.get().is_finite());
            if *old == 20 {
                close(phase.spin_rad_s, after_impact.get());
                close(spin(phase.body).y, after_impact.get());
            }
            Ok::<_, &'static str>((support_law(phase.normal_force_n), old + 1))
        },
        |old, phase, _| {
            order.borrow_mut().push(('i', *old, phase.time_s));
            close(phase.spin_before_rad_s, 4.0);
            let impulse = -(3.0 * PI / 16.0) * IMPACT_RADIUS * IMPACT_MU * phase.normal_impulse_ns;
            after_impact.set(4.0 + impulse / INERTIA);
            Ok::<_, u8>((impact_patch(), old + 1))
        },
        || false,
    )
    .unwrap();
    assert_eq!(impact_history, 11);
    assert_eq!(support_history, 20 + result.support.support_phases);
    assert_eq!(order.borrow()[0].0, 'i');
    assert!(order.borrow()[1..].iter().all(|e| e.0 == 's'));
    close(order.borrow()[0].2, order.borrow()[1].2);
    audit(initial, &result);

    support_history = 20;
    impact_history = 10;
    let impact_calls = Cell::new(0);
    let support_calls = Cell::new(0);
    let refusal = advance_sphere_plane_hertz_contact_with(
        &s,
        initial,
        Vec3::ZERO,
        &mut support_history,
        &mut impact_history,
        0.6,
        c,
        |old, phase, _| {
            support_calls.set(support_calls.get() + 1);
            if *old == 21 {
                Err("late support refusal")
            } else {
                Ok((support_law(phase.normal_force_n), old + 1))
            }
        },
        |old, _, _| {
            impact_calls.set(impact_calls.get() + 1);
            Ok::<_, u8>((impact_patch(), old + 1))
        },
        || false,
    );
    assert!(matches!(
        refusal,
        Err(SpherePlaneHertzContactError::Support(
            SpherePlaneHertzSupportError::Source("late support refusal")
        ))
    ));
    assert_eq!(support_calls.get(), 2);
    assert_eq!(impact_calls.get(), 1);
    assert_eq!((support_history, impact_history), (20, 10));

    let cancel = Cell::new(false);
    let cancelled = advance_sphere_plane_hertz_contact_with(
        &s,
        initial,
        Vec3::ZERO,
        &mut support_history,
        &mut impact_history,
        0.6,
        c,
        |old, phase, _| {
            cancel.set(true);
            Ok::<_, ()>((support_law(phase.normal_force_n), old + 1))
        },
        |old, _, _| Ok::<_, u8>((impact_patch(), old + 1)),
        || cancel.get(),
    );
    assert!(cancelled.is_err());
    assert_eq!((support_history, impact_history), (20, 10));
}

#[test]
fn late_impact_refusal_discards_the_private_bounce_and_retry_has_no_duplicate_event() {
    let s = scene();
    let initial = SpherePlaneHertzState::new(state(
        Vec3::new(0.0, RADIUS, 0.0),
        Vec3::new(0.0, -2.0, 0.0),
        Vec3::new(0.0, 4.0, 0.0),
        UnitQuaternion::IDENTITY,
    ));
    let mut support_history = 5usize;
    let mut impact_history = 0usize;
    let observed = RefCell::new(Vec::new());
    let failed = advance_sphere_plane_hertz_contact_with(
        &s,
        initial,
        Vec3::ZERO,
        &mut support_history,
        &mut impact_history,
        0.25,
        controls(),
        |_, _, _| -> Result<_, &'static str> {
            panic!("positive restitution must still be bouncing")
        },
        |old, phase, _| {
            observed.borrow_mut().push(phase.time_s);
            if *old == 1 {
                Err(41u8)
            } else {
                Ok((impact_patch(), old + 1))
            }
        },
        || false,
    );
    assert!(matches!(
        failed,
        Err(SpherePlaneHertzContactError::Impact(
            SpherePlaneHertzImpactError::Source(41)
        ))
    ));
    assert_eq!((support_history, impact_history), (5, 0));
    assert_eq!(observed.borrow().len(), 2);
    close(observed.borrow()[0], 0.0);
    close(observed.borrow()[1], 0.2);
    let retried = advance_sphere_plane_hertz_contact_with(
        &s,
        initial,
        Vec3::ZERO,
        &mut support_history,
        &mut impact_history,
        0.25,
        controls(),
        |_, _, _| -> Result<_, &'static str> {
            panic!("positive restitution must still be bouncing")
        },
        |old, _, _| Ok::<_, u8>((impact_patch(), old + 1)),
        || false,
    )
    .unwrap();
    assert_eq!((support_history, impact_history), (5, 2));
    assert_eq!(retried.impacts.len(), 2);
    assert_eq!(retried.support.support_phases, 0);
    for (event, time) in retried.impacts.iter().zip(observed.borrow().iter()) {
        close(event.phase.time_s, *time);
    }
    audit(initial, &retried);
}
