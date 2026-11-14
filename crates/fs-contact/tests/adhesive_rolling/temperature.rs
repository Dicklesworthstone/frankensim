//! Original temperature claims change the live spring without resetting memory.
#[path = "../adhesive_temperature/common.rs"]
mod thermal;

use super::*;
use fs_contact::normal_patch::adhesive::temperature::{
    AdhesiveTemperatures, SourceJkrThermalRefresh,
};
use fs_tribo::rolling_loss::adhesive::AdhesiveRollingError;
pub(super) use thermal::Fixture;

pub(super) fn thermal_fixture() -> Fixture {
    Fixture::with_options(10.0, 0.05, false, 320.0)
}
pub(super) fn thermal_source(fixture: &Fixture) -> SourceJkrRollingContact {
    SourceJkrRollingContact::new(fixture.jkr(), 0.0, controls()).unwrap()
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
    source: &SourceJkrRollingContact,
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
fn original_independent_temperatures_change_live_rolling_stiffness_and_parameter_work() {
    let fixture = thermal_fixture();
    // Positive approach avoids the special delta=0 cancellation of E* in k.
    let source = SourceJkrRollingContact::new(fixture.jkr(), 1e-5, controls()).unwrap();
    let c = source.law().critical_rolling_angle_rad();
    let yielded = source.advance(&source.initial_state(), 3.0 * c).unwrap();
    let previous = source.advance(&yielded.state(), -0.25 * c).unwrap().state();
    assert!(yielded.rolling().irreversible_loss_j > 0.0);
    let normal = refresh(&source, &fixture, warm_temperatures());
    let before = *normal.before();
    let after = *normal.response();
    let normal_work = normal.parameter_work_j();
    let step = source.continue_at_temperature(&previous, normal).unwrap();
    let updated = step.source();
    let report = step.report();
    // Independent manufactured claims: E_a(310)=9e6, E_b(305)=19e6,
    // nu_a=nu_b=1/4; W(315)=0.05*(1-0.03*15).
    let reduced_e = ((1.0 - 0.25_f64.powi(2)) * (1.0 / 9e6 + 1.0 / 19e6)).recip();
    let w = 0.0275;
    let r = source.normal_source().effective_radius_m();
    let a0 = (4.5 * PI * w * r * r / reduced_e).cbrt();
    let expected_k = 6.0 * PI * w * (after.contact_radius_m / a0).powf(1.5);
    let xi = previous.elastic_shift_m();
    let expected_work = 0.5 * (expected_k - source.law().stiffness_n_per_m()) * xi * xi;
    close(updated.law().parameters().reduced_modulus_pa, reduced_e);
    close(updated.law().parameters().work_of_adhesion_j_m2, w);
    close(updated.law().stiffness_n_per_m(), expected_k);
    close(report.rolling_storage_after_j, 0.5 * expected_k * xi * xi);
    close(report.rolling_parameter_work_j, expected_work);
    close(report.parameter_work_j, normal_work + expected_work);
    assert_eq!(report.normal_parameter_work_j, normal_work);
    assert_eq!(report.before, before);
    assert_eq!(report.after, after);
    assert_eq!(report.normal_storage_before_j, before.reversible_energy_j);
    assert_eq!(report.normal_storage_after_j, after.reversible_energy_j);
    assert_eq!(
        report.rolling_storage_before_j,
        source.stored_energy_j(&previous).unwrap()
    );
    assert_eq!(report.original_model_identity, source.identity());
    assert_eq!(report.refreshed_model_identity, updated.identity());
    assert_ne!(updated.identity(), source.identity());
    assert_ne!(after.contact_radius_m, before.contact_radius_m);
    assert_eq!(updated.controls(), source.controls());
    assert_eq!(step.previous_state(), previous);
    assert_eq!(step.state().elastic_shift_m().to_bits(), xi.to_bits());
    assert_eq!(
        step.state().accumulated_angle_rad().to_bits(),
        previous.accumulated_angle_rad().to_bits()
    );
    let held = updated.advance(&step.state(), 0.0).unwrap();
    close(held.rolling().endpoint_moment_n_m, -expected_k * r * xi);
    assert_eq!(held.rolling().irreversible_loss_j, 0.0);
    assert_eq!(held.rolling().plastic_rolling_increment_m, 0.0);
    assert_eq!(held.rolling().applied_work_j, 0.0);
}

#[test]
fn same_temperature_and_round_trip_preserve_full_loaded_memory_and_reverse_path() {
    let fixture = thermal_fixture();
    let source = thermal_source(&fixture);
    let c = source.law().critical_rolling_angle_rad();
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
    assert_eq!(same.report().rolling_parameter_work_j, 0.0);
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
    // Includes retained double-double low bits.
    thermal::same_after_refresh(&back.state(), &previous);
    assert_eq!(
        back.report().normal_parameter_work_j,
        -hot.report().normal_parameter_work_j
    );
    assert_eq!(
        back.report().rolling_parameter_work_j,
        -hot.report().rolling_parameter_work_j
    );
    let expected = source.advance(&previous, -3.5 * c).unwrap();
    thermal::same_after_refresh(
        &back.source().advance(&back.state(), -3.5 * c).unwrap(),
        &expected,
    );
    assert!(expected.rolling().irreversible_loss_j > 0.0);
    assert_eq!(
        source.stored_energy_j(&hot.state()),
        Err(SourceJkrRollingError::StateMismatch)
    );
    let (returned, memory, report) = back.into_parts();
    thermal::same_after_refresh(&returned, &source);
    thermal::same_after_refresh(&memory, &previous);
    assert_eq!(report.refreshed_model_identity, returned.identity());
}

#[test]
fn foreign_owner_control_claim_domain_and_hot_small_shift_refuse_without_history_change() {
    let fixture = thermal_fixture();
    let source = thermal_source(&fixture);
    let c = source.law().critical_rolling_angle_rad();
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
        Err(SourceJkrRollingError::StateMismatch)
    );
    assert_eq!(
        source.continue_at_temperature(
            &foreign.initial_state(),
            refresh(&source, &fixture, warm_temperatures())
        ),
        Err(SourceJkrRollingError::StateMismatch)
    );
    // Zero held load does zero actuator work even while approach changes. It
    // must still be refused by this specifically held-indentation owner.
    let held_load = source
        .normal_source()
        .refresh_at_held_load(0.0, fixture.cards(), warm_temperatures())
        .unwrap();
    assert_eq!(held_load.mechanical_work_j(), 0.0);
    assert_eq!(
        source.continue_at_temperature(&previous, held_load),
        Err(SourceJkrRollingError::StateMismatch)
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
        Err(SourceJkrRollingError::StateMismatch)
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
    // A real thermal decrease of area violates the original xi_c/a ceiling.
    let hot = refresh(&source, &fixture, warm_temperatures());
    let old_ratio =
        controls().critical_rolling_displacement_m / source.normal_response().contact_radius_m;
    let hot_ratio = controls().critical_rolling_displacement_m / hot.response().contact_radius_m;
    assert!(hot_ratio > old_ratio);
    let strict = SourceJkrRollingContact::new(
        fixture.jkr(),
        0.0,
        JkrRollingControls {
            maximum_critical_displacement_to_contact_radius: 0.5 * (old_ratio + hot_ratio),
            ..controls()
        },
    )
    .unwrap();
    let loaded = strict
        .advance(&strict.initial_state(), 2.5 * c)
        .unwrap()
        .state();
    assert!(matches!(
        strict.continue_at_temperature(&loaded, hot),
        Err(SourceJkrRollingError::Rolling(
            AdhesiveRollingError::OutsideSmallShift { .. }
        ))
    ));
    assert_eq!(source, original);
    assert_eq!(
        source.advance(&previous, -c).unwrap(),
        original.advance(&previous, -c).unwrap()
    );
    assert!(strict.stored_energy_j(&loaded).unwrap() > 0.0);
}
