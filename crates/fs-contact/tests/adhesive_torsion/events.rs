//! G1/G3/G4: original mode-I events, retained twist release and fresh junctions.
use super::*;
use fs_contact::normal_patch::adhesive::torsion::events::*;
use fs_contact::normal_patch::adhesive::{JkrBranch, JkrJumpKind};

fn loaded(multiplier: f64) -> (SourceJkrTorsionEvents, JkrTorsionEventState) {
    let normal = normal("event-torsion-A");
    let delta = normal.evaluate_force_controlled(0.0).unwrap().indentation_m;
    let source = SourceJkrTorsionEvents::new(normal, controls()).unwrap();
    let initial = source.initial_attached_state(delta).unwrap();
    let state = source
        .advance_torsion(&initial, multiplier * controls().critical_twist_angle_rad)
        .unwrap()
        .into_parts()
        .0;
    (source, state)
}

#[test]
fn g1_separation_releases_the_original_fold_and_remaining_twist_energy() {
    let (source, before) = loaded(2.5);
    let normal = source.normal_source();
    let fold = normal.minimum_attached_indentation_m();
    let step = source
        .advance_normal(&before, fold - normal.indentation_scale_m())
        .unwrap();
    let report = step.report();
    let elastic = normal.normal_state().elastic();
    let g = [elastic.surface_a_state(), elastic.surface_b_state()]
        .map(|s| s.elastic_state().shear_modulus_pa());
    let reduced_g = (g[0].recip() + g[1].recip()).recip();
    let e = elastic.reduced_modulus_pa();
    let w = normal.normal_state().adhesion_energy_j_per_m2();
    let r = normal.effective_radius_m();
    let a_s = (PI * w * r * r / (8.0 * e)).cbrt();
    let phi = before.attached().unwrap().1.elastic_twist_rad();
    // Original normal fold plus two circular-punch compliances in series:
    // a_s^3=pi W R^2/(8 E*), U_phi=(8/3) G_bar a_s^3 phi_e^2.
    // At zero normal force, a_0^3=36 a_s^3, so U_phi falls by 36.
    let twist_release = (PI / 3.0) * (reduced_g / e) * w * r * r * phi * phi;
    let normal_release = 0.4 * PI * w * a_s * a_s;
    close(fold, -3.0 * a_s * a_s / r);
    close(report.torsion_released_energy_j, twist_release);
    close(
        report.torsion_released_energy_j,
        before.torsion_storage_j() / 36.0,
    );
    close(report.normal_released_energy_j, normal_release);
    close(
        report.torsion_parameter_work_j,
        twist_release - before.torsion_storage_j(),
    );
    close(report.released_energy_j, normal_release + twist_release);
    close(
        report.applied_work_j,
        -before.normal_response().reversible_energy_j - before.torsion_storage_j()
            + report.released_energy_j,
    );
    assert!(report.torsion_parameter_work_j < 0.0);
    assert_eq!(report.irreversible_torsion_loss_j, 0.0);
    assert!(before.accumulated_torsion_loss_j() > 0.0);
    assert_eq!(
        step.state().accumulated_torsion_loss_j(),
        before.accumulated_torsion_loss_j()
    );
    let jump = report.normal_transition.unwrap().jump.unwrap();
    assert_eq!(jump.kind, JkrJumpKind::Detach);
    assert_eq!(jump.indentation_m, fold);
    assert_eq!(step.state().normal_response().branch, JkrBranch::Detached);
    assert_eq!(step.state().normal_response().normal_force_n, 0.0);
    assert_eq!(step.state().torsion_storage_j(), 0.0);
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
fn g1_recontact_is_unstrained_and_retains_all_previous_energy_accounts() {
    let (source, loaded) = loaded(-1.5);
    let fold = source.normal_source().minimum_attached_indentation_m();
    let separated = source.advance_normal(&loaded, fold).unwrap();
    let free = source.advance_torsion(separated.state(), 123.0).unwrap();
    assert_eq!(free.state(), separated.state());
    assert_eq!(free.report().applied_work_j, 0.0);
    assert_eq!(free.report().irreversible_torsion_loss_j, 0.0);
    assert!(free.report().torsion.is_none());
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
    assert_eq!(report.torsion_released_energy_j, 0.0);
    assert_eq!(report.torsion_parameter_work_j, 0.0);
    assert_eq!(report.applied_work_j, 0.0);
    let (contact, memory) = attached.state().attached().unwrap();
    assert_eq!(memory, contact.initial_state());
    assert_eq!(memory.elastic_twist_rad(), 0.0);
    assert_eq!(memory.accumulated_angle_rad(), 0.0);
    assert_eq!(attached.state().torsion_storage_j(), 0.0);
    close(
        attached.state().accumulated_normal_release_j(),
        separated.state().accumulated_normal_release_j() + report.normal_released_energy_j,
    );
    assert_eq!(
        attached.state().accumulated_torsion_release_j(),
        separated.state().accumulated_torsion_release_j()
    );
    assert_eq!(
        attached.state().accumulated_torsion_loss_j(),
        loaded.accumulated_torsion_loss_j()
    );
    let angle = 0.25 * controls().critical_twist_angle_rad;
    let next = source.advance_torsion(attached.state(), angle).unwrap();
    let expected = contact.advance(&contact.initial_state(), angle).unwrap();
    assert_eq!(next.report().torsion.unwrap(), expected);
    assert_eq!(
        next.state().accumulated_released_energy_j(),
        attached.state().accumulated_released_energy_j()
    );
    assert_eq!(
        source.advance_torsion(next.state(), 0.0).unwrap().state(),
        next.state()
    );
}

#[test]
fn g3_normal_subdivision_and_threshold_overshoot_keep_original_history_and_work() {
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
        first.report().torsion_parameter_work_j + at_fold.report().torsion_parameter_work_j,
        whole.report().torsion_parameter_work_j,
    );
    assert_eq!(
        at_fold.report().normal_released_energy_j,
        whole.report().normal_released_energy_j
    );
    assert_eq!(
        at_fold.report().torsion_released_energy_j,
        whole.report().torsion_released_energy_j
    );
    assert_eq!(last.report().applied_work_j, 0.0);
    assert_eq!(last.report().released_energy_j, 0.0);
    let detached = source.initial_detached_state(0.0).unwrap();
    let touch = source.advance_normal(&detached, 0.0).unwrap();
    assert_eq!(
        touch.report().normal_transition.unwrap().jump.unwrap().kind,
        JkrJumpKind::Attach
    );
    let overshoot = source.advance_normal(&detached, 1e-6).unwrap();
    let staged = source.advance_normal(touch.state(), 1e-6).unwrap();
    assert_eq!(overshoot.state(), staged.state());
    assert_eq!(
        overshoot.report().normal_released_energy_j,
        touch.report().normal_released_energy_j
    );
    assert_eq!(overshoot.state().torsion_storage_j(), 0.0);
    assert_eq!(
        overshoot.state().attached().unwrap().1.elastic_twist_rad(),
        0.0
    );
}

#[test]
fn g3_marginal_zero_steps_and_weak_twist_history_remain_distinct_from_normal_events() {
    let source = SourceJkrTorsionEvents::new(normal("weak-event-torsion"), controls()).unwrap();
    let fold = source.normal_source().minimum_attached_indentation_m();
    let marginal = source.initial_attached_state(fold).unwrap();
    assert_eq!(
        source.advance_torsion(&marginal, 0.0).unwrap().state(),
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

    let critical = controls().critical_twist_angle_rad;
    let half = source.advance_torsion(&before, 0.5 * critical).unwrap();
    let weak = source
        .advance_torsion(half.state(), 1e-20 * critical)
        .unwrap();
    assert_ne!(weak.state(), half.state());
    assert_eq!(
        weak.state().torsion_storage_j(),
        half.state().torsion_storage_j()
    );
    let unload_normal = source.advance_normal(weak.state(), 0.5e-6).unwrap();
    let restore = source.advance_normal(unload_normal.state(), 1e-6).unwrap();
    assert_eq!(restore.state(), weak.state());
    assert_eq!(restore.state().accumulated_torsion_loss_j(), 0.0);
}

#[test]
fn g4_foreign_history_and_growing_area_applicability_refuse_without_release() {
    let (source, state) = loaded(0.25);
    let snapshot = state.clone();
    let foreign = SourceJkrTorsionEvents::new(normal("foreign-event-torsion"), controls()).unwrap();
    assert!(matches!(
        foreign.advance_normal(&state, 0.0),
        Err(SourceJkrTorsionError::StateMismatch)
    ));
    let changed = SourceJkrTorsionEvents::new(
        source.normal_source().clone(),
        JkrTorsionControls {
            critical_twist_angle_rad: 2.0 * controls().critical_twist_angle_rad,
            ..controls()
        },
    )
    .unwrap();
    assert!(matches!(
        changed.advance_torsion(&state, 0.0),
        Err(SourceJkrTorsionError::StateMismatch)
    ));
    assert!(source.advance_normal(&state, f64::NAN).is_err());
    assert!(source.advance_torsion(&state, f64::INFINITY).is_err());
    assert!(source.initial_detached_state(1e-9).is_err());
    assert_eq!(state, snapshot);

    let normal = normal("event-limit-torsion");
    let original = SourceJkrTorsionContact::new(normal.clone(), 0.0, controls()).unwrap();
    let grown = SourceJkrTorsionContact::new(normal.clone(), 1e-6, controls()).unwrap();
    let w = normal.normal_state().adhesion_energy_j_per_m2();
    let ratio0 = original.maximum_mode_iii_release_j_m2() / w;
    let ratio1 = grown.maximum_mode_iii_release_j_m2() / w;
    assert!(ratio1 > ratio0);
    let limited = SourceJkrTorsionEvents::new(
        normal,
        JkrTorsionControls {
            maximum_mode_iii_energy_release_to_adhesion: 0.5 * (ratio0 + ratio1),
            ..controls()
        },
    )
    .unwrap();
    let initial = limited.initial_attached_state(0.0).unwrap();
    let loaded = limited
        .advance_torsion(&initial, controls().critical_twist_angle_rad)
        .unwrap();
    let detached = limited.initial_detached_state(-1e-6).unwrap();
    for accepted in [loaded.state(), &detached] {
        let snapshot = accepted.clone();
        assert!(matches!(
            limited.advance_normal(accepted, 1e-6),
            Err(SourceJkrTorsionError::OutsideFrozenArea { .. })
        ));
        assert_eq!(accepted, &snapshot);
        assert_eq!(accepted.accumulated_released_energy_j(), 0.0);
        assert_eq!(
            limited.advance_torsion(accepted, 0.0).unwrap().state(),
            accepted
        );
    }
}
