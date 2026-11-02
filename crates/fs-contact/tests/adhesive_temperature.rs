//! Independent thermal normal-adhesion equilibria and original-source refusals.
#[path = "adhesive_temperature/common.rs"]
mod common;
#[path = "adhesive_temperature/finite_range.rs"]
mod finite_range;

use core::f64::consts::PI;
use fs_contact::normal_patch::adhesive::fixture::JkrLoadingError;
use fs_contact::normal_patch::adhesive::temperature::{
    AdhesiveTemperatureError, AdhesiveTemperatures, AdhesiveThermalControl,
};
use fs_contact::normal_patch::adhesive::{JkrBranch, JkrError, SourceJkrSphereContact};
use fs_material::state_point::{MaterialPropertySelection, resolve_isotropic_solid_state_point};

fn temperatures(a: f64, b: f64, interface: f64) -> AdhesiveTemperatures {
    AdhesiveTemperatures {
        surface_a_k: a,
        surface_b_k: b,
        interface_k: interface,
    }
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * a.abs().max(b.abs()).max(1e-20),
        "{a:.17e} != {b:.17e}"
    );
}

fn modulus(a: f64, b: f64) -> f64 {
    let ea = 1e6 * (1.0 - 0.01 * (a - 300.0));
    let eb = 2e6 * (1.0 - 0.01 * (b - 300.0));
    1.0 / ((1.0 - 0.25_f64.powi(2)) / ea + (1.0 - 0.25_f64.powi(2)) / eb)
}

// Direct dimensional JKR equations in radius coordinates, independent of the
// production scaled-polynomial inverse and its stored energy coefficients.
fn physical(e: f64, w: f64, r: f64, a: f64) -> (f64, f64, f64) {
    let c = (2.0 * PI * w / e).sqrt();
    let delta = a * a / r - c * a.sqrt();
    let force = 4.0 * e * a.powi(3) / (3.0 * r) - 2.0 * e * c * a.powf(1.5);
    let energy = 8.0 * e * a.powi(5) / (15.0 * r * r) - 4.0 * e * c * a.powf(3.5) / (3.0 * r)
        + PI * w * a * a;
    (delta, force, energy)
}

#[test]
fn unequal_bulk_and_interface_temperatures_change_actual_jkr_area_force_and_potential() {
    let fixture = common::Fixture::new();
    let law = fixture.jkr();
    let original = law.clone();
    let before = law.evaluate_attached(0.0).unwrap();
    let update = law
        .refresh_at_held_indentation(
            before.state(),
            fixture.cards(),
            temperatures(310.0, 290.0, 315.0),
        )
        .unwrap();
    let normal = update.source().normal_state();
    let e = modulus(310.0, 290.0);
    let w = 0.55;
    let radius = (2.0 * PI * w / e).cbrt();
    let (_, force, potential) = physical(e, w, 1.0, radius);
    close(normal.elastic().reduced_modulus_pa(), e, 3e-15);
    close(
        normal.elastic().surface_a_state().young_modulus_pa(),
        9e5,
        3e-15,
    );
    close(
        normal.elastic().surface_b_state().young_modulus_pa(),
        2.2e6,
        3e-15,
    );
    close(normal.adhesion_energy_j_per_m2(), w, 3e-15);
    close(update.response().contact_radius_m, radius, 3e-13);
    close(
        update.response().contact_area_m2,
        PI * radius * radius,
        6e-13,
    );
    close(update.response().normal_force_n, force, 2e-12);
    close(update.response().reversible_energy_j, potential, 3e-12);
    assert!(update.response().normal_force_n < 0.0);
    assert!(update.response().contact_area_m2 < before.contact_area_m2);
    assert_eq!(update.response().indentation_m, before.indentation_m);
    assert_eq!(update.response().branch, before.branch);
    assert_eq!(update.original_model_identity(), law.identity());
    assert_eq!(update.response().normal_state_identity(), normal.identity());
    assert_eq!(update.control(), AdhesiveThermalControl::HeldIndentation);
    assert_eq!(update.mechanical_work_j(), 0.0);
    close(
        update.potential_change_j(),
        potential - before.reversible_energy_j,
        3e-12,
    );
    assert_eq!(update.parameter_work_j(), update.potential_change_j());
    for (old, new) in [
        (
            law.normal_state().elastic().surface_a_state(),
            normal.elastic().surface_a_state(),
        ),
        (
            law.normal_state().elastic().surface_b_state(),
            normal.elastic().surface_b_state(),
        ),
    ] {
        assert_eq!(
            old.resolved().axis_quantities(),
            new.resolved().axis_quantities()
        );
        assert_eq!(new.resolved().state_coordinate("p"), Some(2e5));
        for (old, new) in old
            .resolved()
            .properties()
            .iter()
            .zip(new.resolved().properties())
        {
            assert_eq!(old.answer().receipt.selected, new.answer().receipt.selected);
        }
    }
    assert_eq!(normal.interface_state().state_coordinate("T"), Some(315.0));
    assert_eq!(law, original);
    assert_eq!(update.clone().into_source(), *update.source());
}

#[test]
fn held_load_retains_servo_work_and_noop_retains_original_receipts() {
    let fixture = common::Fixture::new();
    let law = fixture.jkr();
    let before = law.evaluate_attached(0.0).unwrap();
    let noop = law
        .refresh_at_held_indentation(
            before.state(),
            fixture.cards(),
            temperatures(300.0, 300.0, 300.0),
        )
        .unwrap();
    assert_eq!(noop.source(), &law);
    assert_eq!(noop.before(), noop.response());
    assert_eq!(noop.potential_change_j(), 0.0);
    for force in [0.0, -0.5 * PI, 1.0] {
        let update = law
            .refresh_at_held_load(force, fixture.cards(), temperatures(310.0, 290.0, 315.0))
            .unwrap();
        let e = modulus(310.0, 290.0);
        let w = 0.55;
        let radius = ((3.0 / (4.0 * e))
            * (force + 3.0 * PI * w + (6.0 * PI * w * force + (3.0 * PI * w).powi(2)).sqrt()))
        .cbrt();
        let (delta, physical_force, energy) = physical(e, w, 1.0, radius);
        close(update.response().contact_radius_m, radius, 3e-12);
        close(update.response().indentation_m, delta, 3e-12);
        assert!((physical_force - force).abs() < 1e-12);
        assert!((update.response().normal_force_n - force).abs() < 1e-12);
        close(update.response().reversible_energy_j, energy, 3e-12);
        let work = force * (delta - update.before().indentation_m);
        close(update.mechanical_work_j(), work, 3e-12);
        close(
            update.parameter_work_j(),
            energy - update.before().reversible_energy_j - work,
            3e-12,
        );
        assert_eq!(
            update.control(),
            AdhesiveThermalControl::HeldLoad { force_n: force }
        );
        assert!(update.response().contact_area_m2 > 0.0);
    }
}

#[test]
fn original_claims_stay_pinned_and_expiry_never_selects_an_eligible_replacement() {
    let fixture = common::Fixture::with_options(1.0, 1.0, true, 320.0);
    let law = fixture.jkr();
    let updated = law
        .at_temperatures(fixture.cards(), temperatures(315.0, 315.0, 315.0))
        .unwrap();
    close(
        updated
            .normal_state()
            .elastic()
            .surface_a_state()
            .young_modulus_pa(),
        8.5e5,
        3e-15,
    );
    close(
        updated.normal_state().adhesion_energy_j_per_m2(),
        0.55,
        3e-15,
    );
    assert_eq!(
        law.normal_state().interface_state().properties()[0]
            .answer()
            .receipt
            .selected,
        updated.normal_state().interface_state().properties()[0]
            .answer()
            .receipt
            .selected,
    );
    let hot = common::point().with("T", 315.0).unwrap();
    assert!(
        resolve_isotropic_solid_state_point(
            &fixture.cards[0],
            &hot,
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_err()
    );
    let expired = common::Fixture::with_options(1.0, 1.0, true, 310.0);
    let old = expired.jkr();
    let saved = old.clone();
    assert!(matches!(
        old.at_temperatures(expired.cards(), temperatures(315.0, 300.0, 300.0)),
        Err(AdhesiveTemperatureError::Material(_))
    ));
    assert!(matches!(
        old.at_temperatures(expired.cards(), temperatures(300.0, 300.0, 315.0)),
        Err(AdhesiveTemperatureError::Material(_))
    ));
    assert_eq!(old, saved);
}

#[test]
fn thermal_pull_off_refuses_while_detached_and_original_state_tokens_are_preserved() {
    let fixture = common::Fixture::new();
    let law = fixture.jkr();
    let attached = law
        .evaluate_attached(0.9 * law.minimum_attached_indentation_m())
        .unwrap();
    let hot = temperatures(300.0, 300.0, 320.0);
    assert!(matches!(
        law.refresh_at_held_indentation(attached.state(), fixture.cards(), hot),
        Err(AdhesiveTemperatureError::Jkr(
            JkrError::NoAttachedEquilibrium { .. }
        ))
    ));
    assert!(matches!(
        law.refresh_at_held_load(
            0.9 * law.force_controlled_pull_off_n(),
            fixture.cards(),
            hot
        ),
        Err(AdhesiveTemperatureError::JkrLoading(
            JkrLoadingError::NoForceEquilibrium { .. }
        ))
    ));
    let detached = law.detached_state(-1e-4).unwrap();
    let update = law
        .refresh_at_held_indentation(detached, fixture.cards(), hot)
        .unwrap();
    assert_eq!(update.response().branch, JkrBranch::Detached);
    assert_eq!(update.response().contact_area_m2, 0.0);
    assert_eq!(update.potential_change_j(), 0.0);
    let foreign =
        SourceJkrSphereContact::new(fixture.normal.clone(), 0.9, 1e-7, common::envelope()).unwrap();
    assert!(matches!(
        law.refresh_at_held_indentation(
            foreign.evaluate_attached(0.0).unwrap().state(),
            fixture.cards(),
            hot
        ),
        Err(AdhesiveTemperatureError::Jkr(JkrError::StateMismatch))
    ));
    let retry = law
        .refresh_at_held_indentation(
            attached.state(),
            fixture.cards(),
            temperatures(300.0, 300.0, 300.0),
        )
        .unwrap();
    assert_eq!(retry.response(), &attached);
}

#[test]
fn source_domains_original_card_roles_and_the_original_jkr_envelope_still_refuse() {
    let fixture = common::Fixture::new();
    let law = fixture.jkr();
    for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY, 330.0] {
        for t in [
            temperatures(invalid, 300.0, 300.0),
            temperatures(300.0, invalid, 300.0),
            temperatures(300.0, 300.0, invalid),
        ] {
            assert!(law.at_temperatures(fixture.cards(), t).is_err());
        }
    }
    assert!(
        law.at_temperatures(
            [&fixture.cards[1], &fixture.cards[0]],
            temperatures(300.0, 300.0, 300.0)
        )
        .is_err()
    );
    let same_names_different_card = common::Fixture::with_options(1.01, 1.0, false, 320.0);
    assert!(
        law.at_temperatures(
            same_names_different_card.cards(),
            temperatures(300.0, 300.0, 300.0)
        )
        .is_err()
    );
    let e = law.normal_state().elastic().reduced_modulus_pa();
    let range = (1.0 / (e * e)).cbrt() / 6.0;
    let near_limit =
        SourceJkrSphereContact::new(fixture.normal.clone(), 1.0, range, common::envelope())
            .unwrap();
    close(near_limit.tabor_parameter(), 6.0, 3e-14);
    assert!(matches!(
        near_limit.at_temperatures(fixture.cards(), temperatures(300.0, 300.0, 320.0)),
        Err(AdhesiveTemperatureError::Jkr(
            JkrError::OutsideApplicability { .. }
        ))
    ));
}
