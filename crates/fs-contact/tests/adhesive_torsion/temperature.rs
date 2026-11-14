//! Independently refreshed bulk/interface sources retain loaded bonded twist.
#[path = "../adhesive_temperature/common.rs"]
mod thermal;

use super::*;
use fs_contact::normal_patch::adhesive::temperature::{
    AdhesiveTemperatures, SourceJkrThermalRefresh,
};
pub(super) use thermal::Fixture;

pub(super) fn thermal_fixture() -> Fixture {
    Fixture::with_options(10.0, 0.05, false, 320.0)
}
pub(super) fn thermal_source(fixture: &Fixture) -> SourceJkrTorsionContact {
    SourceJkrTorsionContact::new(fixture.jkr(), 0.0, controls()).unwrap()
}
pub(super) fn warm_temperatures() -> AdhesiveTemperatures {
    AdhesiveTemperatures {
        surface_a_k: 310.0,
        surface_b_k: 305.0,
        interface_k: 315.0,
    }
}
fn room_temperatures() -> AdhesiveTemperatures {
    AdhesiveTemperatures {
        surface_a_k: 300.0,
        surface_b_k: 300.0,
        interface_k: 300.0,
    }
}
fn refresh(
    source: &SourceJkrTorsionContact,
    fixture: &Fixture,
    temperatures: AdhesiveTemperatures,
) -> SourceJkrThermalRefresh {
    source
        .normal_source()
        .refresh_at_held_indentation(
            source.normal_response().state(),
            fixture.cards(),
            temperatures,
        )
        .unwrap()
}

#[test]
fn original_bulk_temperatures_and_interface_adhesion_change_live_torsion_and_work() {
    let fixture = thermal_fixture();
    let source = SourceJkrTorsionContact::new(fixture.jkr(), 1e-5, controls()).unwrap();
    let c = controls().critical_twist_angle_rad;
    let yielded = source.advance(&source.initial_state(), 3.0 * c).unwrap();
    let previous = source.advance(&yielded.state(), -0.25 * c).unwrap().state();
    assert!(yielded.torsion().irreversible_loss_j > 0.0);
    let normal = refresh(&source, &fixture, warm_temperatures());
    let before = *normal.before();
    let after = *normal.response();
    let normal_work = normal.parameter_work_j();
    let step = source.continue_at_temperature(&previous, normal).unwrap();
    let updated = step.source();
    let report = step.report();
    // Independent circular-punch compliances in series. Original curves give
    // E_a(310)=9e6, E_b(305)=19e6, nu=1/4 and W(315)=0.0275.
    let ga: f64 = 9e6 / 2.5;
    let gb: f64 = 19e6 / 2.5;
    let g = (ga.recip() + gb.recip()).recip();
    let a = after.contact_radius_m;
    let expected_k = (3.0 / (16.0 * a.powi(3)) * (ga.recip() + gb.recip())).recip();
    let phi = previous.elastic_twist_rad();
    let expected_work = 0.5 * (expected_k - source.law().stiffness_nm_per_rad()) * phi * phi;
    close(updated.law().parameters().shear_moduli_pa[0], ga);
    close(updated.law().parameters().shear_moduli_pa[1], gb);
    close(updated.law().stiffness_nm_per_rad(), expected_k);
    close(
        updated
            .normal_source()
            .normal_state()
            .adhesion_energy_j_per_m2(),
        0.0275,
    );
    close(
        updated.maximum_mode_iii_release_j_m2(),
        4.0 * g * a * c * c / PI,
    );
    close(report.torsion_storage_after_j, 0.5 * expected_k * phi * phi);
    close(report.torsion_parameter_work_j, expected_work);
    close(report.parameter_work_j, normal_work + expected_work);
    assert_eq!(report.normal_parameter_work_j, normal_work);
    assert_eq!(report.before, before);
    assert_eq!(report.after, after);
    assert_eq!(report.normal_storage_before_j, before.reversible_energy_j);
    assert_eq!(report.normal_storage_after_j, after.reversible_energy_j);
    assert_eq!(
        report.torsion_storage_before_j,
        source.stored_energy_j(&previous).unwrap()
    );
    assert_eq!(report.original_model_identity, source.identity());
    assert_eq!(report.refreshed_model_identity, updated.identity());
    assert_ne!(updated.identity(), source.identity());
    assert_ne!(after.contact_radius_m, before.contact_radius_m);
    assert_eq!(updated.controls(), source.controls());
    assert_eq!(step.previous_state(), previous);
    assert_eq!(step.state().elastic_twist_rad().to_bits(), phi.to_bits());
    assert_eq!(
        step.state().accumulated_angle_rad().to_bits(),
        previous.accumulated_angle_rad().to_bits()
    );
    let held = updated.advance(&step.state(), 0.0).unwrap();
    close(held.torsion().endpoint_moment_n_m, -expected_k * phi);
    assert_eq!(held.torsion().irreversible_loss_j, 0.0);
    assert_eq!(held.torsion().plastic_twist_increment_rad, 0.0);
    assert_eq!(held.torsion().applied_work_j, 0.0);
}

#[test]
fn temperature_no_op_and_return_restore_exact_loaded_twist_and_reverse_yield() {
    let fixture = thermal_fixture();
    let source = thermal_source(&fixture);
    let c = controls().critical_twist_angle_rad;
    let yielded = source.advance(&source.initial_state(), 3.0 * c).unwrap();
    let unloaded = source.advance(&yielded.state(), -0.25 * c).unwrap();
    let previous = source
        .advance(&unloaded.state(), 0.125 * f64::EPSILON * c)
        .unwrap()
        .state();
    let same = source
        .continue_at_temperature(&previous, refresh(&source, &fixture, room_temperatures()))
        .unwrap();
    assert_eq!(same.source(), &source);
    assert_eq!(same.state(), previous);
    assert_eq!(same.report().normal_parameter_work_j, 0.0);
    assert_eq!(same.report().torsion_parameter_work_j, 0.0);
    assert_eq!(same.report().parameter_work_j, 0.0);
    let hot = source
        .continue_at_temperature(&previous, refresh(&source, &fixture, warm_temperatures()))
        .unwrap();
    let back = hot
        .source()
        .continue_at_temperature(
            &hot.state(),
            refresh(hot.source(), &fixture, room_temperatures()),
        )
        .unwrap();
    // The way back re-resolves the cards with pinned-claim receipts, so the
    // identities differ from the original's single-claim-only ones; the
    // source's values and the loaded memory must return exactly.
    thermal::same_after_refresh(back.source(), &source);
    // Full private double-double memory.
    thermal::same_after_refresh(&back.state(), &previous);
    assert_eq!(
        back.report().normal_parameter_work_j,
        -hot.report().normal_parameter_work_j
    );
    assert_eq!(
        back.report().torsion_parameter_work_j,
        -hot.report().torsion_parameter_work_j
    );
    let expected = source.advance(&previous, -3.5 * c).unwrap();
    thermal::same_after_refresh(
        &back.source().advance(&back.state(), -3.5 * c).unwrap(),
        &expected,
    );
    assert!(expected.torsion().irreversible_loss_j > 0.0);
    assert_eq!(
        source.stored_energy_j(&hot.state()),
        Err(SourceJkrTorsionError::StateMismatch)
    );
    let (returned, memory, report) = back.into_parts();
    thermal::same_after_refresh(&returned, &source);
    thermal::same_after_refresh(&memory, &previous);
    assert_eq!(report.refreshed_model_identity, returned.identity());
}

#[test]
fn wrong_owner_control_domain_and_new_mode_iii_release_refuse_atomically() {
    let fixture = thermal_fixture();
    let source = thermal_source(&fixture);
    let c = controls().critical_twist_angle_rad;
    let previous = source
        .advance(&source.initial_state(), 2.5 * c)
        .unwrap()
        .state();
    let original = source.clone();
    let foreign_fixture = Fixture::with_options(11.0, 0.05, false, 320.0);
    let foreign = thermal_source(&foreign_fixture);
    assert_eq!(
        source.continue_at_temperature(
            &previous,
            refresh(&foreign, &foreign_fixture, warm_temperatures())
        ),
        Err(SourceJkrTorsionError::StateMismatch)
    );
    assert_eq!(
        source.continue_at_temperature(
            &foreign.initial_state(),
            refresh(&source, &fixture, warm_temperatures())
        ),
        Err(SourceJkrTorsionError::StateMismatch)
    );
    let held_load = source
        .normal_source()
        .refresh_at_held_load(0.0, fixture.cards(), warm_temperatures())
        .unwrap();
    assert_eq!(held_load.mechanical_work_j(), 0.0);
    assert_eq!(
        source.continue_at_temperature(&previous, held_load),
        Err(SourceJkrTorsionError::StateMismatch)
    );
    let wrong_endpoint = source
        .normal_source()
        .refresh_at_held_indentation(
            source
                .normal_source()
                .evaluate_attached(1e-5)
                .unwrap()
                .state(),
            fixture.cards(),
            warm_temperatures(),
        )
        .unwrap();
    assert_eq!(
        source.continue_at_temperature(&previous, wrong_endpoint),
        Err(SourceJkrTorsionError::StateMismatch)
    );
    assert!(
        source
            .normal_source()
            .refresh_at_held_indentation(
                source.normal_response().state(),
                fixture.cards(),
                AdhesiveTemperatures {
                    surface_a_k: 321.0,
                    ..room_temperatures()
                }
            )
            .is_err()
    );
    assert!(
        source
            .normal_source()
            .refresh_at_held_indentation(
                source.normal_response().state(),
                [&fixture.cards[1], &fixture.cards[0]],
                warm_temperatures()
            )
            .is_err()
    );
    // W weakens faster than G and area. Even though the normal JKR endpoint
    // remains valid, the retained frozen-area torsion approximation can fail.
    let hot = source
        .continue_at_temperature(&previous, refresh(&source, &fixture, warm_temperatures()))
        .unwrap();
    let old_ratio = source.maximum_mode_iii_release_j_m2()
        / source
            .normal_source()
            .normal_state()
            .adhesion_energy_j_per_m2();
    let hot_ratio = hot.source().maximum_mode_iii_release_j_m2()
        / hot
            .source()
            .normal_source()
            .normal_state()
            .adhesion_energy_j_per_m2();
    assert!(hot_ratio > old_ratio);
    let strict = SourceJkrTorsionContact::new(
        fixture.jkr(),
        0.0,
        JkrTorsionControls {
            maximum_mode_iii_energy_release_to_adhesion: 0.5 * (old_ratio + hot_ratio),
            ..controls()
        },
    )
    .unwrap();
    let loaded = strict
        .advance(&strict.initial_state(), 2.5 * c)
        .unwrap()
        .state();
    assert!(matches!(
        strict.continue_at_temperature(&loaded, refresh(&strict, &fixture, warm_temperatures())),
        Err(SourceJkrTorsionError::OutsideFrozenArea { .. })
    ));
    assert_eq!(source, original);
    assert_eq!(
        source.advance(&previous, -c).unwrap(),
        original.advance(&previous, -c).unwrap()
    );
    assert!(strict.stored_energy_j(&loaded).unwrap() > 0.0);
}
