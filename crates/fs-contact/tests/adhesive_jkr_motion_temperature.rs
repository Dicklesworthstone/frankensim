//! Source-preserving temperature continuation of live original normal motion.
#[path = "adhesive_temperature/common.rs"]
mod common;

use core::f64::consts::PI;
use fs_contact::normal_patch::adhesive::{
    JkrBranch, SourceJkrSphereContact,
    motion::{JkrNormalMotion, JkrNormalMotionBudget, JkrNormalMotionError, MAX_JKR_NORMAL_EVENTS},
    temperature::{AdhesiveTemperatureError, AdhesiveTemperatures},
};
use fs_exec::CancelGate;

const MASS: f64 = 1.2e-5;

fn temperatures(a: f64, b: f64, interface: f64) -> AdhesiveTemperatures {
    AdhesiveTemperatures {
        surface_a_k: a,
        surface_b_k: b,
        interface_k: interface,
    }
}
fn scales(source: &SourceJkrSphereContact) -> (f64, f64, f64, f64) {
    let h = source.indentation_scale_m();
    let force =
        PI * source.normal_state().adhesion_energy_j_per_m2() * source.effective_radius_m() / 6.0;
    let speed = (force * h / MASS).sqrt();
    (h, force, speed, h / speed)
}
fn budget(source: &SourceJkrSphereContact) -> JkrNormalMotionBudget {
    let (h, force, speed, time) = scales(source);
    JkrNormalMotionBudget {
        maximum_step_s: 10.0 * time,
        maximum_substep_s: 0.01 * time,
        maximum_substeps: 4096,
        maximum_root_iterations: 48,
        maximum_events: MAX_JKR_NORMAL_EVENTS,
        minimum_indentation_m: -100.0 * h,
        maximum_accumulated_travel_m: 200.0 * h,
        maximum_normal_speed_m_s: 100.0 * speed,
        maximum_force_n: 1e4 * force,
        coordinate_absolute_tolerance_m: 1e-14 * h,
        coordinate_relative_tolerance: 1e-13,
        energy_absolute_tolerance_j: 1e-11 * force * h,
        energy_relative_tolerance: 1e-10,
        impulse_absolute_tolerance_n_s: 1e-12 * MASS * speed,
        impulse_relative_tolerance: 1e-11,
    }
}
fn motion(source: &SourceJkrSphereContact) -> JkrNormalMotion<'_> {
    let (h, force, _, time) = scales(source);
    JkrNormalMotion::new_with_damped_support(
        source,
        MASS,
        budget(source),
        0.05 * force / h,
        0.5 * h,
        0.2 * MASS / time,
    )
    .unwrap()
}
fn close(a: f64, b: f64, scale: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * scale,
        "{a:e} != {b:e}; scaled error {:e}",
        (a - b).abs() / scale
    );
}

#[test]
fn refreshed_materials_change_actual_normal_dynamics_without_restarting_motion_or_support() {
    let fixture = common::Fixture::new();
    let source = fixture.jkr();
    let law = motion(&source);
    let (_, force, speed, time) = scales(&source);
    let initial = law
        .initialize(source.evaluate_attached(0.0).unwrap().state(), 0.3 * speed)
        .unwrap();
    let accepted = law
        .advance(
            initial,
            0.2 * force,
            0.1 * time,
            &CancelGate::new_clock_free(),
        )
        .unwrap()
        .state();
    let old_support = law.support_energy_j(accepted).unwrap();
    let refresh = source
        .refresh_at_held_indentation(
            accepted.contact_state(),
            fixture.cards(),
            temperatures(310.0, 290.0, 315.0),
        )
        .unwrap();
    let candidate = law
        .continue_at_temperature(accepted, &refresh, &CancelGate::new_clock_free())
        .unwrap();
    assert_eq!(candidate.previous_state(), accepted);
    assert_eq!(
        candidate.state().velocity_m_s().to_bits(),
        accepted.velocity_m_s().to_bits()
    );
    assert_eq!(
        candidate.state().time_s().to_bits(),
        accepted.time_s().to_bits()
    );
    assert_eq!(candidate.state().steps(), accepted.steps());
    assert_eq!(
        candidate.state().contact_state().indentation_m().to_bits(),
        accepted.contact_state().indentation_m().to_bits()
    );
    assert_eq!(
        candidate.state().contact_state().branch(),
        accepted.contact_state().branch()
    );
    assert_eq!(
        candidate.motion().effective_mass_kg(),
        law.effective_mass_kg()
    );
    assert_eq!(candidate.motion().budget(), law.budget());
    assert_eq!(
        candidate.motion().support_stiffness_n_per_m(),
        law.support_stiffness_n_per_m()
    );
    assert_eq!(
        candidate.motion().support_zero_force_indentation_m(),
        law.support_zero_force_indentation_m()
    );
    assert_eq!(
        candidate.motion().support_damping_n_s_per_m(),
        law.support_damping_n_s_per_m()
    );
    assert_eq!(
        candidate
            .motion()
            .support_energy_j(candidate.state())
            .unwrap(),
        old_support
    );
    assert_ne!(
        candidate.report().before.normal_force_n,
        candidate.report().after.normal_force_n
    );
    assert_ne!(
        candidate.report().before.contact_area_m2,
        candidate.report().after.contact_area_m2
    );
    assert_eq!(candidate.report().before, *refresh.before());
    assert_eq!(candidate.report().after, *refresh.response());
    assert_eq!(
        candidate.report().parameter_work_j,
        refresh.parameter_work_j()
    );
    assert_eq!(
        candidate.report().storage_change_j,
        refresh.potential_change_j()
    );
    assert_eq!(candidate.report().energy_residual_j, 0.0);
    let (next, state, _) = candidate.into_parts();
    let evolved = next
        .advance(state, 0.0, 0.1 * time, &CancelGate::new_clock_free())
        .unwrap();
    let frozen = law
        .advance(accepted, 0.0, 0.1 * time, &CancelGate::new_clock_free())
        .unwrap();
    assert_ne!(
        evolved.state().velocity_m_s(),
        frozen.state().velocity_m_s()
    );
    assert_eq!(evolved.state().steps(), accepted.steps() + 1);
    assert_eq!(evolved.state().time_s(), accepted.time_s() + 0.1 * time);
    // Independently rebuild the driver at the original refresh's new source.
    // Only this comparison's clock begins at zero; its physical evolution agrees.
    let independent = JkrNormalMotion::new_with_damped_support(
        refresh.source(),
        MASS,
        law.budget(),
        law.support_stiffness_n_per_m(),
        law.support_zero_force_indentation_m(),
        law.support_damping_n_s_per_m(),
    )
    .unwrap();
    let start = independent
        .initialize(refresh.response().state(), accepted.velocity_m_s())
        .unwrap();
    let expected = independent
        .advance(start, 0.0, 0.1 * time, &CancelGate::new_clock_free())
        .unwrap();
    assert_eq!(
        evolved.state().contact_state(),
        expected.state().contact_state()
    );
    assert_eq!(
        evolved.state().velocity_m_s(),
        expected.state().velocity_m_s()
    );
    assert_eq!(evolved.normal_impulse_n_s, expected.normal_impulse_n_s);
    assert_eq!(
        evolved.damping_dissipation_j,
        expected.damping_dissipation_j
    );
}

#[test]
fn unchanged_temperature_keeps_the_exact_complete_token_and_next_trace() {
    let fixture = common::Fixture::new();
    let source = fixture.jkr();
    let law = motion(&source);
    let (_, _, speed, time) = scales(&source);
    let initial = law
        .initialize(source.evaluate_attached(0.0).unwrap().state(), 0.2 * speed)
        .unwrap();
    let accepted = law
        .advance(initial, 0.0, 0.2 * time, &CancelGate::new_clock_free())
        .unwrap()
        .state();
    let refresh = source
        .refresh_at_held_indentation(
            accepted.contact_state(),
            fixture.cards(),
            temperatures(300.0, 300.0, 300.0),
        )
        .unwrap();
    let next = law
        .continue_at_temperature(accepted, &refresh, &CancelGate::new_clock_free())
        .unwrap();
    assert_eq!(next.state(), accepted);
    assert_eq!(
        next.report().original_model_identity,
        next.report().refreshed_model_identity
    );
    assert_eq!(next.report().before, next.report().after);
    assert_eq!(next.report().parameter_work_j, 0.0);
    assert_eq!(next.report().storage_change_j, 0.0);
    let ordinary = law
        .advance_with_dissipation_trace(accepted, 0.0, 0.1 * time, &CancelGate::new_clock_free())
        .unwrap();
    let continued = next
        .motion()
        .advance_with_dissipation_trace(
            next.state(),
            0.0,
            0.1 * time,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(continued, ordinary);
}

#[test]
fn reversible_temperature_cycle_reports_signed_potential_work_and_no_motion() {
    let fixture = common::Fixture::new();
    let source = fixture.jkr();
    let law = motion(&source);
    let initial = law
        .initialize(source.evaluate_attached(0.0).unwrap().state(), 0.0)
        .unwrap();
    let warm = source
        .refresh_at_held_indentation(
            initial.contact_state(),
            fixture.cards(),
            temperatures(310.0, 290.0, 315.0),
        )
        .unwrap();
    let hot = law
        .continue_at_temperature(initial, &warm, &CancelGate::new_clock_free())
        .unwrap();
    assert!(hot.report().parameter_work_j > 0.0);
    let cool = warm
        .source()
        .refresh_at_held_indentation(
            hot.state().contact_state(),
            fixture.cards(),
            temperatures(300.0, 300.0, 300.0),
        )
        .unwrap();
    let cold = hot
        .motion()
        .continue_at_temperature(hot.state(), &cool, &CancelGate::new_clock_free())
        .unwrap();
    assert!(cold.report().parameter_work_j < 0.0);
    let (_, force, _, _) = scales(&source);
    close(
        hot.report().parameter_work_j + cold.report().parameter_work_j,
        0.0,
        force * source.indentation_scale_m(),
        1e-12,
    );
    assert_eq!(cold.state().velocity_m_s(), 0.0);
    assert_eq!(cold.state().time_s(), 0.0);
    assert_eq!(cold.state().steps(), 0);
    assert_eq!(cold.state().contact_state().indentation_m(), 0.0);
    assert_eq!(
        law.response(initial).unwrap().reversible_energy_j,
        cold.report().after.reversible_energy_j
    );
}

#[test]
fn full_owner_foreign_contact_and_held_load_refreshes_refuse_without_changing_history() {
    let fixture = common::Fixture::new();
    let source = fixture.jkr();
    let law = motion(&source);
    let initial = law
        .initialize(source.evaluate_attached(0.0).unwrap().state(), 0.0)
        .unwrap();
    let refresh = source
        .refresh_at_held_indentation(
            initial.contact_state(),
            fixture.cards(),
            temperatures(310.0, 290.0, 315.0),
        )
        .unwrap();
    let foreign = JkrNormalMotion::new(&source, 2.0 * MASS, budget(&source)).unwrap();
    let token = foreign.initialize(initial.contact_state(), 0.0).unwrap();
    assert!(matches!(
        law.continue_at_temperature(token, &refresh, &CancelGate::new_clock_free()),
        Err(JkrNormalMotionError::StateMismatch)
    ));
    let other = source
        .refresh_at_held_indentation(
            source
                .evaluate_attached(0.1 * source.indentation_scale_m())
                .unwrap()
                .state(),
            fixture.cards(),
            temperatures(310.0, 290.0, 315.0),
        )
        .unwrap();
    assert!(matches!(
        law.continue_at_temperature(initial, &other, &CancelGate::new_clock_free()),
        Err(JkrNormalMotionError::StateMismatch)
    ));
    let held = source
        .refresh_at_held_load(0.0, fixture.cards(), temperatures(300.0, 300.0, 300.0))
        .unwrap();
    let held_state = law.initialize(held.before().state(), 0.0).unwrap();
    assert!(matches!(
        law.continue_at_temperature(held_state, &held, &CancelGate::new_clock_free()),
        Err(JkrNormalMotionError::Invalid(_))
    ));
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        law.continue_at_temperature(initial, &refresh, &cancelled),
        Err(JkrNormalMotionError::Cancelled)
    ));
    assert_eq!(law.response(initial).unwrap(), *refresh.before());
    assert!(
        law.continue_at_temperature(initial, &refresh, &CancelGate::new_clock_free())
            .is_ok()
    );
}

#[test]
fn lost_equilibrium_and_refreshed_force_limits_cannot_silently_change_phase_or_policy() {
    let fixture = common::Fixture::new();
    let source = fixture.jkr();
    let law = motion(&source);
    let near = source
        .evaluate_attached(0.99 * source.minimum_attached_indentation_m())
        .unwrap();
    let initial = law.initialize(near.state(), 0.0).unwrap();
    assert!(matches!(
        source.refresh_at_held_indentation(
            initial.contact_state(),
            fixture.cards(),
            temperatures(300.0, 300.0, 320.0)
        ),
        Err(AdhesiveTemperatureError::Jkr(_))
    ));
    assert_eq!(law.response(initial).unwrap(), near);
    let before = source.evaluate_attached(0.0).unwrap();
    let mut limits = budget(&source);
    limits.maximum_force_n = 1.2 * before.normal_force_n.abs();
    let limited = JkrNormalMotion::new(&source, MASS, limits).unwrap();
    let accepted = limited.initialize(before.state(), 0.0).unwrap();
    let stronger = source
        .refresh_at_held_indentation(
            accepted.contact_state(),
            fixture.cards(),
            temperatures(300.0, 300.0, 280.0),
        )
        .unwrap();
    assert!(stronger.response().normal_force_n.abs() > limits.maximum_force_n);
    assert!(matches!(
        limited.continue_at_temperature(accepted, &stronger, &CancelGate::new_clock_free()),
        Err(JkrNormalMotionError::Range(_))
    ));
    assert_eq!(limited.response(accepted).unwrap(), before);
}

#[test]
fn detached_touch_preserves_phase_until_the_new_sources_original_motion_event() {
    let fixture = common::Fixture::new();
    let source = fixture.jkr();
    let law = motion(&source);
    let (_, _, speed, time) = scales(&source);
    let initial = law
        .initialize(source.detached_state(0.0).unwrap(), speed)
        .unwrap();
    let refresh = source
        .refresh_at_held_indentation(
            initial.contact_state(),
            fixture.cards(),
            temperatures(310.0, 290.0, 315.0),
        )
        .unwrap();
    let continued = law
        .continue_at_temperature(initial, &refresh, &CancelGate::new_clock_free())
        .unwrap();
    assert_eq!(
        continued.state().contact_state().branch(),
        JkrBranch::Detached
    );
    assert_eq!(continued.report().parameter_work_j, 0.0);
    let next = continued
        .motion()
        .advance(
            continued.state(),
            0.0,
            0.01 * time,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(next.event_count(), 1);
    assert_eq!(next.events().next().unwrap().offset_s, 0.0);
    close(
        next.released_energy_j,
        -refresh
            .source()
            .evaluate_attached(0.0)
            .unwrap()
            .reversible_energy_j,
        next.released_energy_j,
        1e-13,
    );
}
