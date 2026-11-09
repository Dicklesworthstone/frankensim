//! Loaded thermal continuation retains complete junction-lifecycle accounts.
use super::super::temperature::{Fixture, thermal_fixture, thermal_source, warm_temperatures};
use super::*;
use fs_contact::normal_patch::adhesive::JkrError;

fn complete_history(data: &Fixture) -> (SourceJkrTorsionEvents, JkrTorsionEventState) {
    let contact = thermal_source(data);
    let source =
        SourceJkrTorsionEvents::new(contact.normal_source().clone(), contact.controls()).unwrap();
    let angle = contact.controls().critical_twist_angle_rad;
    let initial = source.initial_attached_state(0.0).unwrap();
    let loaded = source.advance_torsion(&initial, 2.5 * angle).unwrap();
    let separated = source
        .advance_normal(
            loaded.state(),
            source.normal_source().minimum_attached_indentation_m(),
        )
        .unwrap();
    let joined = source
        .advance_normal(
            separated.state(),
            0.5 * source.normal_source().indentation_scale_m(),
        )
        .unwrap();
    let loaded = source
        .advance_torsion(joined.state(), -2.5 * angle)
        .unwrap();
    (source, loaded.into_parts().0)
}

#[test]
fn thermal_refresh_retains_yielded_torsion_history_and_later_normal_events() {
    let data = thermal_fixture();
    let (source, before) = complete_history(&data);
    assert!(before.accumulated_normal_release_j() > 0.0);
    assert!(before.accumulated_torsion_release_j() > 0.0);
    assert!(before.accumulated_torsion_loss_j() > 0.0);
    let refresh = source
        .normal_source()
        .refresh_at_held_indentation(
            before.normal_response().state(),
            data.cards(),
            warm_temperatures(),
        )
        .unwrap();
    let (contact, memory) = before.attached().unwrap();
    let original = contact
        .continue_at_temperature(&memory, refresh.clone())
        .unwrap();
    let step = source.continue_at_temperature(&before, refresh).unwrap();
    assert_eq!(step.previous_state(), &before);
    assert_eq!(step.report(), original.report());
    assert_eq!(
        step.state().attached(),
        Some((original.source(), original.state()))
    );
    let (refreshed, after, report) = step.into_parts();
    assert_eq!(refreshed.normal_source(), original.source().normal_source());
    assert_eq!(refreshed.controls(), source.controls());
    assert_eq!(after.normal_response(), original.source().normal_response());
    let (contact, retained) = after.attached().unwrap();
    assert_eq!(retained.elastic_twist_rad(), memory.elastic_twist_rad());
    assert_eq!(
        retained.accumulated_angle_rad(),
        memory.accumulated_angle_rad()
    );
    assert_eq!(
        after.torsion_storage_j(),
        contact.stored_energy_j(&retained).unwrap()
    );
    assert_eq!(
        after.accumulated_normal_release_j(),
        before.accumulated_normal_release_j()
    );
    assert_eq!(
        after.accumulated_torsion_release_j(),
        before.accumulated_torsion_release_j()
    );
    assert_eq!(
        after.accumulated_torsion_loss_j(),
        before.accumulated_torsion_loss_j()
    );
    assert_eq!(report.before.indentation_m, report.after.indentation_m);
    assert_ne!(report.parameter_work_j, 0.0);
    close(
        report.parameter_work_j,
        report.normal_parameter_work_j + report.torsion_parameter_work_j,
    );
    assert!(matches!(
        source.advance_torsion(&after, 0.0),
        Err(SourceJkrTorsionError::StateMismatch)
    ));
    let separated = refreshed
        .advance_normal(
            &after,
            refreshed.normal_source().minimum_attached_indentation_m(),
        )
        .unwrap();
    assert_eq!(
        separated.state().normal_response().branch,
        JkrBranch::Detached
    );
    assert!(separated.report().normal_released_energy_j > 0.0);
    assert!(separated.report().torsion_released_energy_j > 0.0);
    assert_eq!(
        separated.state().accumulated_torsion_loss_j(),
        before.accumulated_torsion_loss_j()
    );
    let rejoined = refreshed.advance_normal(separated.state(), 0.0).unwrap();
    let (contact, memory) = rejoined.state().attached().unwrap();
    assert_eq!(memory, contact.initial_state());
    assert_eq!(rejoined.state().torsion_storage_j(), 0.0);
    assert_eq!(
        rejoined.state().accumulated_torsion_loss_j(),
        before.accumulated_torsion_loss_j()
    );
    close(
        rejoined.state().accumulated_torsion_release_j(),
        before.accumulated_torsion_release_j() + separated.report().torsion_released_energy_j,
    );
    close(
        rejoined.state().accumulated_normal_release_j(),
        before.accumulated_normal_release_j()
            + separated.report().normal_released_energy_j
            + rejoined.report().normal_released_energy_j,
    );
}

#[test]
fn thermal_torsion_event_refusals_preserve_the_original_complete_history() {
    let data = thermal_fixture();
    let (source, before) = complete_history(&data);
    let saved = before.clone();
    let refresh = source
        .normal_source()
        .refresh_at_held_indentation(
            before.normal_response().state(),
            data.cards(),
            warm_temperatures(),
        )
        .unwrap();
    let mut controls = source.controls();
    controls.critical_twist_angle_rad *= 0.5;
    let foreign = SourceJkrTorsionEvents::new(source.normal_source().clone(), controls).unwrap();
    assert!(matches!(
        foreign.continue_at_temperature(&before, refresh),
        Err(SourceJkrTorsionError::StateMismatch)
    ));
    let held_load = source
        .normal_source()
        .refresh_at_held_load(0.0, data.cards(), warm_temperatures())
        .unwrap();
    assert!(matches!(
        source.continue_at_temperature(&before, held_load),
        Err(SourceJkrTorsionError::StateMismatch)
    ));
    let detached = source
        .advance_normal(
            &before,
            source.normal_source().minimum_attached_indentation_m(),
        )
        .unwrap()
        .into_parts()
        .0;
    let saved_detached = detached.clone();
    let refresh = source
        .normal_source()
        .refresh_at_held_indentation(
            detached.normal_response().state(),
            data.cards(),
            warm_temperatures(),
        )
        .unwrap();
    assert!(matches!(
        source.continue_at_temperature(&detached, refresh),
        Err(SourceJkrTorsionError::Normal(JkrError::InvalidInput {
            field: "attached branch for temperature continuation"
        }))
    ));
    assert_eq!(detached, saved_detached);
    assert_eq!(before, saved);
    assert_eq!(
        source.advance_torsion(&before, 0.0).unwrap().state(),
        &before
    );
}
