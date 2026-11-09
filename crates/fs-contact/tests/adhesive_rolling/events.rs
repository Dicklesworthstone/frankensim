//! G1/G3/G4: original normal events, retained spring release and fresh junctions.
use super::*;

#[path = "events/temperature.rs"]
mod temperature;
use fs_contact::normal_patch::adhesive::rolling::events::*;
use fs_contact::normal_patch::adhesive::{JkrBranch, JkrJumpKind};

fn loaded(multiplier: f64) -> (SourceJkrRollingEvents, JkrRollingEventState) {
    let normal = normal("event-rolling-A");
    let delta = normal.evaluate_force_controlled(0.0).unwrap().indentation_m;
    let source = SourceJkrRollingEvents::new(normal, controls()).unwrap();
    let initial = source.initial_attached_state(delta).unwrap();
    let angle = initial
        .attached()
        .unwrap()
        .0
        .law()
        .critical_rolling_angle_rad();
    let state = source
        .advance_rolling(&initial, multiplier * angle)
        .unwrap()
        .into_parts()
        .0;
    (source, state)
}

#[test]
fn g1_detachment_releases_the_original_spinodal_and_retained_rolling_spring() {
    let (source, before) = loaded(2.5);
    let normal = source.normal_source();
    let fold = normal.minimum_attached_indentation_m();
    let step = source
        .advance_normal(&before, fold - normal.indentation_scale_m())
        .unwrap();
    let report = step.report();
    let w = normal.normal_state().adhesion_energy_j_per_m2();
    let e = normal.normal_state().elastic().reduced_modulus_pa();
    let r = normal.effective_radius_m();
    let a_s = (PI * w * r * r / (8.0 * e)).cbrt();
    let xi = before.attached().unwrap().1.elastic_shift_m();
    // At the displacement spinodal a_s^3/a_0^3=1/36, hence k_xi=pi W.
    let rolling_release = 0.5 * PI * w * xi * xi;
    let normal_release = 0.4 * PI * w * a_s * a_s;
    close(fold, -3.0 * a_s * a_s / r);
    close(report.rolling_released_energy_j, rolling_release);
    close(report.normal_released_energy_j, normal_release);
    close(
        report.rolling_parameter_work_j,
        rolling_release - before.rolling_storage_j(),
    );
    close(report.released_energy_j, normal_release + rolling_release);
    close(
        report.applied_work_j,
        -before.normal_response().reversible_energy_j - before.rolling_storage_j()
            + report.released_energy_j,
    );
    assert!(report.rolling_parameter_work_j < 0.0);
    assert_eq!(report.irreversible_rolling_loss_j, 0.0);
    assert!(before.accumulated_rolling_loss_j() > 0.0);
    assert_eq!(
        step.state().accumulated_rolling_loss_j(),
        before.accumulated_rolling_loss_j()
    );
    let jump = report.normal_transition.unwrap().jump.unwrap();
    assert_eq!(jump.kind, JkrJumpKind::Detach);
    assert_eq!(jump.indentation_m, fold);
    assert_eq!(step.state().normal_response().branch, JkrBranch::Detached);
    assert_eq!(step.state().normal_response().normal_force_n, 0.0);
    assert_eq!(step.state().rolling_storage_j(), 0.0);
    assert!(step.state().attached().is_none());
    assert_eq!(step.previous_state(), &before);
    let farther = source.advance_normal(step.state(), 2.0 * fold).unwrap();
    assert_eq!(farther.report().released_energy_j, 0.0);
    assert_eq!(farther.report().applied_work_j, 0.0);
    assert_eq!(
        farther.state().accumulated_released_energy_j(),
        report.released_energy_j
    );
}

#[test]
fn g1_reattachment_starts_stress_free_without_forgetting_prior_release_or_plastic_loss() {
    let (source, loaded) = loaded(-1.5);
    let fold = source.normal_source().minimum_attached_indentation_m();
    let separated = source.advance_normal(&loaded, fold).unwrap();
    let free = source.advance_rolling(separated.state(), 123.0).unwrap();
    assert_eq!(free.state(), separated.state());
    assert_eq!(free.report().applied_work_j, 0.0);
    assert_eq!(free.report().irreversible_rolling_loss_j, 0.0);
    assert!(free.report().rolling.is_none());
    let attached = source.advance_normal(free.state(), 0.0).unwrap();
    let report = attached.report();
    let jump = report.normal_transition.unwrap().jump.unwrap();
    assert_eq!(jump.kind, JkrJumpKind::Attach);
    assert_eq!(jump.indentation_m, 0.0);
    close(
        report.normal_released_energy_j,
        -source
            .normal_source()
            .evaluate_attached(0.0)
            .unwrap()
            .reversible_energy_j,
    );
    assert_eq!(report.rolling_released_energy_j, 0.0);
    assert_eq!(report.rolling_parameter_work_j, 0.0);
    assert_eq!(report.applied_work_j, 0.0);
    let (contact, memory) = attached.state().attached().unwrap();
    assert_eq!(memory, contact.initial_state());
    assert_eq!(memory.elastic_shift_m(), 0.0);
    assert_eq!(memory.accumulated_angle_rad(), 0.0);
    assert_eq!(attached.state().rolling_storage_j(), 0.0);
    close(
        attached.state().accumulated_normal_release_j(),
        separated.state().accumulated_normal_release_j() + report.normal_released_energy_j,
    );
    assert_eq!(
        attached.state().accumulated_rolling_release_j(),
        separated.state().accumulated_rolling_release_j()
    );
    assert_eq!(
        attached.state().accumulated_rolling_loss_j(),
        loaded.accumulated_rolling_loss_j()
    );
    let next = source
        .advance_rolling(
            attached.state(),
            0.25 * contact.law().critical_rolling_angle_rad(),
        )
        .unwrap();
    let expected = contact
        .advance(
            &contact.initial_state(),
            0.25 * contact.law().critical_rolling_angle_rad(),
        )
        .unwrap();
    assert_eq!(next.report().rolling.unwrap(), expected);
    assert_eq!(
        next.state().accumulated_released_energy_j(),
        attached.state().accumulated_released_energy_j()
    );
    assert_eq!(
        source.advance_rolling(next.state(), 0.0).unwrap().state(),
        next.state()
    );
}

#[test]
fn g3_normal_subdivision_and_threshold_overshoot_preserve_source_history_and_work() {
    let (source, before) = loaded(0.5);
    let fold = source.normal_source().minimum_attached_indentation_m();
    let end = 1.5 * fold;
    let whole = source.advance_normal(&before, end).unwrap();
    let mid = 0.5 * (fold + before.normal_response().indentation_m);
    let first = source.advance_normal(&before, mid).unwrap();
    let original = before
        .attached()
        .unwrap()
        .0
        .advance_normal(&before.attached().unwrap().1, mid)
        .unwrap();
    let (original_source, original_state) = original.into_parts();
    let actual = first.state().attached().unwrap();
    assert_eq!(actual.0, &original_source);
    assert_eq!(actual.1, original_state);
    let at_fold = source.advance_normal(first.state(), fold).unwrap();
    let last = source.advance_normal(at_fold.state(), end).unwrap();
    assert_eq!(last.state(), whole.state());
    close(
        first.report().normal_contact_work_j + at_fold.report().normal_contact_work_j,
        whole.report().normal_contact_work_j,
    );
    close(
        first.report().rolling_parameter_work_j + at_fold.report().rolling_parameter_work_j,
        whole.report().rolling_parameter_work_j,
    );
    assert_eq!(
        at_fold.report().normal_released_energy_j,
        whole.report().normal_released_energy_j
    );
    assert_eq!(
        at_fold.report().rolling_released_energy_j,
        whole.report().rolling_released_energy_j
    );
    assert_eq!(last.report().applied_work_j, 0.0);
    assert_eq!(last.report().released_energy_j, 0.0);
    let detached = source.initial_detached_state(0.0).unwrap();
    let touch = source.advance_normal(&detached, 0.0).unwrap();
    assert_eq!(
        touch.report().normal_transition.unwrap().jump.unwrap().kind,
        JkrJumpKind::Attach
    );
    assert_eq!(touch.state().normal_response().branch, JkrBranch::Attached);
}

#[test]
fn g3_zero_angular_steps_authenticate_even_the_marginal_branch_and_retain_weak_work() {
    let source = SourceJkrRollingEvents::new(normal("weak-event-rolling"), controls()).unwrap();
    let fold = source.normal_source().minimum_attached_indentation_m();
    let marginal = source.initial_attached_state(fold).unwrap();
    assert_eq!(
        source.advance_rolling(&marginal, 0.0).unwrap().state(),
        &marginal
    );
    assert!(
        source
            .advance_normal(&marginal, fold)
            .unwrap()
            .state()
            .attached()
            .is_none()
    );
    let before = source.initial_attached_state(1e-6).unwrap();
    let after = 1e-6_f64.next_up();
    let expected = source
        .normal_source()
        .reversible_work_between(1e-6, after)
        .unwrap();
    let step = source.advance_normal(&before, after).unwrap();
    assert_ne!(expected, 0.0);
    assert_eq!(step.report().normal_contact_work_j, expected);
    assert_eq!(step.report().applied_work_j, expected);
    assert_eq!(step.report().storage_change_j, expected);
    assert_eq!(step.report().energy_residual_j, 0.0);
}

#[test]
fn g4_event_radius_applicability_and_foreign_history_refuse_without_publication() {
    let (source, state) = loaded(0.25);
    let snapshot = state.clone();
    let foreign = SourceJkrRollingEvents::new(normal("foreign-event-rolling"), controls()).unwrap();
    assert!(matches!(
        foreign.advance_normal(&state, 0.0),
        Err(SourceJkrRollingError::StateMismatch)
    ));
    let changed = SourceJkrRollingEvents::new(
        source.normal_source().clone(),
        JkrRollingControls {
            critical_rolling_displacement_m: 2.0 * controls().critical_rolling_displacement_m,
            ..controls()
        },
    )
    .unwrap();
    assert!(matches!(
        changed.advance_rolling(&state, 0.0),
        Err(SourceJkrRollingError::StateMismatch)
    ));
    assert!(source.advance_normal(&state, f64::NAN).is_err());
    assert!(source.advance_rolling(&state, f64::INFINITY).is_err());
    assert!(source.initial_detached_state(1e-9).is_err());
    assert_eq!(state, snapshot);

    let normal = normal("event-limit-rolling");
    let fold = normal.minimum_attached_indentation_m();
    let a_s = normal.evaluate_attached(fold).unwrap().contact_radius_m;
    let limited = SourceJkrRollingEvents::new(
        normal,
        JkrRollingControls {
            critical_rolling_displacement_m: 1.5
                * controls().maximum_critical_displacement_to_contact_radius
                * a_s,
            ..controls()
        },
    )
    .unwrap();
    let accepted = limited.initial_attached_state(0.0).unwrap();
    let snapshot = accepted.clone();
    assert!(matches!(
        limited.advance_normal(&accepted, fold),
        Err(SourceJkrRollingError::Rolling(
            fs_tribo::rolling_loss::adhesive::AdhesiveRollingError::OutsideSmallShift { .. }
        ))
    ));
    assert_eq!(accepted, snapshot);
    assert_eq!(accepted.accumulated_released_energy_j(), 0.0);
}
